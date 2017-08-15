package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.api.services.genomics.model.Read;
import com.google.common.collect.Iterators;
import hdfs.jsr203.HadoopFileSystem;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CigarUtil;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextComparator;
import javafx.collections.transformation.SortedList;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.broadcast.Broadcast;
import org.apache.spark.sql.execution.CollapseCodegenStages$;
import org.bdgenomics.adam.rdd.ADAMContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.ReadsDataSource;
import org.broadinstitute.hellbender.engine.Shard;
import org.broadinstitute.hellbender.engine.ShardBoundary;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.SparkSharder;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SATagBuilder;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.collections.IntervalsSkipList;
import org.broadinstitute.hellbender.utils.gcs.BamBucketIoUtils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.ojalgo.function.BinaryFunction;
import org.seqdoop.hadoop_bam.BAMInputFormat;
import scala.Tuple2;

import javax.annotation.Nullable;
import java.io.Serializable;
import java.util.*;
import java.util.function.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Created by valentin on 7/14/17.
 */
@CommandLineProgramProperties(summary = "test", oneLineSummary = "test",
        programGroup = StructuralVariationSparkProgramGroup.class )
@BetaFeature
public class ComposeStructuralVariantHaplotypesSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    public static final String CONTIGS_FILE_SHORT_NAME = "C";
    public static final String CONTIGS_FILE_FULL_NAME = "contigs";
    public static final String SHARD_SIZE_SHORT_NAME = "sz";
    public static final String SHARD_SIZE_FULL_NAME = "shardSize";
    public static final String PADDING_SIZE_SHORT_NAME = "pd";
    public static final String PADDING_SIZE_FULL_NAME = "paddingSize";

    public static final int DEFAULT_SHARD_SIZE = 10_000;
    public static final int DEFAULT_PADDING_SIZE = 50;


    @Argument(doc = "shard size",
              shortName = SHARD_SIZE_SHORT_NAME,
              fullName = SHARD_SIZE_FULL_NAME,
    optional = true)
    private int shardSize = DEFAULT_SHARD_SIZE;

    @Argument(doc ="padding size",
              shortName = PADDING_SIZE_SHORT_NAME,
              fullName = PADDING_SIZE_FULL_NAME,
              optional = true)
    private int paddingSize = DEFAULT_PADDING_SIZE;

    @Argument(doc = "input variant file",
              shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
              fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME)
    private String variantsFileName;

    @Argument(doc = "aligned contig file",
              fullName = CONTIGS_FILE_FULL_NAME,
              shortName = CONTIGS_FILE_SHORT_NAME
    )
    private String alignedContigsFileName;

    @Argument(doc = "output bam file with contigs per variant",
              fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFileName;

    @Override
    protected void onStartup() {
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        Utils.nonNull(ctx);
        Utils.nonNull(alignedContigsFileName);
        final ReadsSparkSource alignedContigs = new ReadsSparkSource(ctx);
        final VariantsSparkSource variantsSource = new VariantsSparkSource(ctx);

        final JavaRDD<GATKRead> contigs = alignedContigs.getParallelReads(alignedContigsFileName, referenceArguments.getReferenceFileName(), getIntervals());
        final JavaRDD<VariantContext> variants = variantsSource.getParallelVariantContexts(variantsFileName, getIntervals()).filter(ComposeStructuralVariantHaplotypesSpark::supportedVariant);

        final JavaPairRDD<VariantContext, List<GATKRead>> variantOverlappingContigs = composeOverlappingContigs(ctx, contigs, variants);
        processVariants(variantOverlappingContigs, getReferenceSequenceDictionary(), alignedContigs);

    }

    private static boolean supportedVariant(final VariantContext vc) {
        final List<Allele> alternatives = vc.getAlternateAlleles();
        if (alternatives.size() != 1) {
            return false;
        } else {
            final Allele alternative = alternatives.get(0);
            if (!alternative.isSymbolic()) {
                return false;
            } else if (alternative.getDisplayString().equals("<INS>")) {
                return true;
            } else if (alternative.getDisplayString().equals("<DEL>")) {
                return true;
            } else {
                return false;
            }
        }
    }

    private JavaPairRDD<VariantContext,List<GATKRead>> composeOverlappingContigs(final JavaSparkContext ctx, final JavaRDD<GATKRead> contigs, final JavaRDD<VariantContext> variants) {
        final SAMSequenceDictionary sequenceDictionary = getBestAvailableSequenceDictionary();
        final List<SimpleInterval> intervals = hasIntervals() ? getIntervals() : IntervalUtils.getAllIntervalsForReference(sequenceDictionary);
        // use unpadded shards (padding is only needed for reference bases)
        final List<ShardBoundary> shardBoundaries = intervals.stream()
                .flatMap(interval -> Shard.divideIntervalIntoShards(interval, shardSize, 0, sequenceDictionary).stream())
                .collect(Collectors.toList());
        final IntervalsSkipList<SimpleInterval> shardIntervals = new IntervalsSkipList<>(shardBoundaries.stream()
                .map(ShardBoundary::getPaddedInterval)
                .collect(Collectors.toList()));
        final Broadcast<SAMSequenceDictionary> dictionaryBroadcast = ctx.broadcast(sequenceDictionary);

        final Broadcast<IntervalsSkipList<SimpleInterval>> shardIntervalsBroadcast = ctx.broadcast(shardIntervals);

        final JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval,GATKRead>>> contigsInShards =
            groupInShards(contigs, ComposeStructuralVariantHaplotypesSpark::contigIntervals, shardIntervalsBroadcast);
        final int paddingSize = this.paddingSize;

        final JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval, VariantContext>>> variantsInShards =
            groupInShards(variants, (v) -> variantsBreakPointIntervals(v, paddingSize, dictionaryBroadcast.getValue()), shardIntervalsBroadcast);

        final JavaPairRDD<SimpleInterval, Tuple2<List<Tuple2<SimpleInterval, GATKRead>>, List<Tuple2<SimpleInterval, VariantContext>>>> contigAndVariantsInShards =
                contigsInShards.join(variantsInShards);


        final JavaPairRDD<VariantContext, List<GATKRead>> contigsPerVariantInterval =
                contigAndVariantsInShards.flatMapToPair(t -> {
                    final List<Tuple2<SimpleInterval, VariantContext>> vars = t._2()._2();
                    final List<Tuple2<SimpleInterval, GATKRead>> ctgs = t._2()._1();
                    return vars.stream()
                            .map(v -> {
                                final List<GATKRead> cs = ctgs.stream()
                                        .filter(ctg -> v._1().overlaps(ctg._1()))
                                        .map(Tuple2::_2)
                                        .collect(Collectors.toList());

                                return new Tuple2<>(v._2(), cs);})
                            .collect(Collectors.toList()).iterator();
                });

        final Function<VariantContext, String> variantId = (Function<VariantContext, String> & Serializable) ComposeStructuralVariantHaplotypesSpark::variantId;
        final Function2<List<GATKRead>, List<GATKRead>, List<GATKRead>> readListMerger = (a, b) -> Stream.concat(a.stream(), b.stream()).collect(Collectors.toList());

        // Merge contig lists on the same variant-context coming from different intervals
        // into one.
        final JavaPairRDD<VariantContext, List<GATKRead>> contigsPerVariant = contigsPerVariantInterval
                .mapToPair(t -> new Tuple2<>(variantId.apply(t._1()), t))
                .reduceByKey((a, b) -> new Tuple2<>(a._1(), readListMerger.call(a._2(), b._2())))
                .mapToPair(Tuple2::_2);
        return contigsPerVariant;
    }

    private static String variantId(final VariantContext variant) {
        if (variant.getID() != null && !variant.getID().isEmpty()) {
            return variant.getID();
        } else {
            final int length = Math.abs(variantLength(variant));
            return "var_" + variant.getAlternateAllele(0).getDisplayString() + "_" + length;
        }
    }

    private static List<SimpleInterval> contigIntervals(final GATKRead contig) {
        if (contig.isUnmapped()) {
            return Collections.emptyList();
        } else {
            return Collections.singletonList(new SimpleInterval(contig.getContig(), contig.getStart(), contig.getEnd()));
        }
    }

    private static List<SimpleInterval> variantsBreakPointIntervals(final VariantContext variant, final int padding, final SAMSequenceDictionary dictionary) {
        final String contigName = variant.getContig();
        final int contigLength = dictionary.getSequence(contigName).getSequenceLength();
        if (variant.getAlternateAllele(0).getDisplayString().equals("<INS>")) {
            return Collections.singletonList(new SimpleInterval(contigName, Math.max(1, variant.getStart() - padding), Math.min(variant.getStart() + 1 + padding, contigLength)));
        } else { // must be <DEL>
            final int length = - variantLength(variant);
            return Arrays.asList(new SimpleInterval(contigName, Math.max(1, variant.getStart() - padding) , Math.min(contigLength, variant.getStart() + padding)),
                                 new SimpleInterval(contigName, Math.max(1, variant.getStart() + length - padding), Math.min(contigLength, variant.getStart() + length + padding)));
        }
    }

    private static int variantLength(VariantContext variant) {
        final int length = variant.getAttributeAsInt("SVLEN", 0);
        if (length == 0) {
            throw new IllegalStateException("missing SVLEN annotation in " + variant.getContig() + ":" + variant.getStart());
        }
        return length;
    }

    private static SimpleInterval locatableToSimpleInterval(final Locatable loc) {
        return new SimpleInterval(loc.getContig(), loc.getStart(), loc.getEnd());
    }

    private static JavaRDD<GATKRead> readsByInterval(final SimpleInterval interval, final ReadsSparkSource alignedContigs, final String alignedContigFileName, final String referenceName) {

        return alignedContigs.getParallelReads(alignedContigFileName,
                referenceName, Collections.singletonList(interval));
    }

    private <T> JavaPairRDD<SimpleInterval, List<Tuple2<SimpleInterval, T>>> groupInShards(final JavaRDD<T> elements, final org.apache.spark.api.java.function.Function<T, List<SimpleInterval>> intervalsOf,
                                                                  final Broadcast<IntervalsSkipList<SimpleInterval>> shards) {
        final PairFlatMapFunction<T, SimpleInterval, Tuple2<SimpleInterval, T>> flatMapIntervals =
                t -> intervalsOf.call(t).stream().map(i -> new Tuple2<>(i, new Tuple2<>(i,t))).iterator();

        return elements
                .flatMapToPair(flatMapIntervals)
                .flatMapToPair(t -> shards.getValue().getOverlapping(t._1()).stream().map(i -> new Tuple2<>(i, t._2())).iterator())
                .aggregateByKey(new ArrayList<Tuple2<SimpleInterval, T>>(10),
                        (l1, c) -> { l1.add(c); return l1;},
                        (l1, l2) -> {l1.addAll(l2); return l1;});
    }

    protected void processVariants(final JavaPairRDD<VariantContext, List<GATKRead>> variantsAndOverlappingContigs, final SAMSequenceDictionary dictionary, final ReadsSparkSource s) {

        final SAMFileHeader outputHeader = new SAMFileHeader();
        final SAMProgramRecord programRecord = new SAMProgramRecord(getProgramName());
        programRecord.setCommandLine(getCommandLine());
        outputHeader.setSequenceDictionary(dictionary);
        outputHeader.addProgramRecord(programRecord);
        outputHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        outputHeader.addReadGroup(new SAMReadGroupRecord("CTG"));
        final SAMFileWriter outputWriter = BamBucketIoUtils.makeWriter(outputFileName, outputHeader, true);


        final BinaryOperator<GATKRead> readMerger = (BinaryOperator<GATKRead> & Serializable)
                ComposeStructuralVariantHaplotypesSpark::reduceContigReads;

        final VariantContextComparator variantComparator = new VariantContextComparator(dictionary);

        final JavaPairRDD<VariantContext, List<GATKRead>> variantsAndOverlappingUniqueContigs
                = variantsAndOverlappingContigs
                .mapValues(l -> l.stream().collect(Collectors.groupingBy(GATKRead::getName)))
                .mapValues(m -> m.values().stream()
                        .map(l -> l.stream().reduce(readMerger).orElseThrow(IllegalStateException::new))
                        .collect(Collectors.toList()))
                .sortByKey(variantComparator);

        Utils.stream(variantsAndOverlappingUniqueContigs.toLocalIterator())
                .map(t -> resolvePendingContigs(t, s))
                .forEach(t -> {
                    final StructuralVariantContext svc = new StructuralVariantContext(t._1());
                    final List<GATKRead> contigs = t._2();
                    final int maxLength = contigs.stream()
                            .mapToInt(GATKRead::getLength)
                            .max().orElse(paddingSize);
                    final SimpleInterval referenceInterval = new SimpleInterval(
                            t._1().getContig(), (int) Math.floor(t._1().getStart() - maxLength * 2.0), (int) Math.ceil(t._1().getEnd() + maxLength * 2.0));
                    final Haplotype referenceHaplotype = svc.composeHaplotype(0, maxLength * 2, getReference());
                    referenceHaplotype.setGenomeLocation(null);
                    final Haplotype alternativeHaplotype = svc.composeHaplotype(1, maxLength * 2, getReference());
                    alternativeHaplotype.setGenomeLocation(null);
                    final String idPrefix = String.format("var_%s_%d", t._1().getContig(), t._1().getStart());
                    final Consumer<SAMRecord> haplotypeExtraSetup = r -> {
                        r.setReferenceName(t._1().getContig());
                        r.setAlignmentStart(t._1().getStart());
                        r.setAttribute(SAMTag.RG.name(), "HAP");
                    };
                    outputWriter.addAlignment(referenceHaplotype.convertToSAMRecord(outputHeader, idPrefix + ":ref",
                            haplotypeExtraSetup));
                    outputWriter.addAlignment(alternativeHaplotype.convertToSAMRecord(outputHeader, idPrefix + ":alt",
                            haplotypeExtraSetup));
                    contigs.forEach(c -> {
                        c.clearAttributes();
                        c.setName( idPrefix + ":" + c.getName());
                        c.setIsPaired(false);
                        c.setIsDuplicate(false);
                        c.setIsSecondaryAlignment(false);
                        c.setReadGroup("CTG");
                        c.setCigar("*");
                        if (c.isReverseStrand()) {
                            c.setIsReverseStrand(false);
                            final byte[] bases = c.getBases();
                            SequenceUtil.reverseComplement(bases);
                            c.setBases(bases);
                            final byte[] quals = c.getBaseQualities();
                            if (quals != null && quals.length > 0) {
                                SequenceUtil.reverseQualities(quals);
                                c.setBaseQualities(quals);
                            }
                        }
                        c.setPosition(t._1().getContig(), t._1().getStart());
                        c.setMappingQuality(0);
                        c.setMateIsUnmapped();
                        c.setIsUnmapped();
                        outputWriter.addAlignment(c.convertToSAMRecord(outputHeader));
                    });
                });
        outputWriter.close();
    }

    private Tuple2<VariantContext, List<GATKRead>> resolvePendingContigs(final Tuple2<VariantContext, List<GATKRead>> vc, final ReadsSparkSource s) {
        logger.debug("VC " + vc._1().getContig() + ":" + vc._1().getStart() + "-" + vc._1().getEnd());
        final List<GATKRead> updatedList = vc._2().stream()
                .map(r -> {
                    if (!r.getCigar().containsOperator(CigarOperator.H)) {
                        return r;
                    } else {
                        final SATagBuilder saTagBuilder = new SATagBuilder(r);
                        final String targetName = r.getName();
                        final List<SimpleInterval> locations = saTagBuilder.getArtificialReadsBasedOnSATag(s.getHeader(alignedContigsFileName, referenceArguments.getReferenceFileName()))
                                .stream().map(rr -> new SimpleInterval(rr.getContig(), rr.getStart(), rr.getEnd()))
                                .collect(Collectors.toList());
                        final GATKRead candidate = s.getParallelReads(alignedContigsFileName, referenceArguments.getReferenceFileName(), locations)
                                .filter(rr -> rr.getName().equals(targetName))
                                .reduce(((Function2< GATKRead, GATKRead, GATKRead> & Serializable) ComposeStructuralVariantHaplotypesSpark::reduceContigReads));

                        final GATKRead canonicRead = reduceContigReads(candidate, r);
                        if (canonicRead.getCigar().containsOperator(CigarOperator.H)) {
                            logger.warn("Contig " + canonicRead.getName() + " " + readAlignmentString(canonicRead) + " " + canonicRead.getAttributeAsString("SA") + " gave-up!");
                            return null;
                        } else {
                            return canonicRead;
                        }
                    }})
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
        return new Tuple2<>(vc._1, updatedList);
    }

    private static String readAlignmentString(GATKRead r) {
        return r.getContig() + "," + r.getStart() + "," + (r.isReverseStrand() ?  "-" : "+") + "," + r.getCigar() + "," + r.getMappingQuality() + "," + r.getAttributeAsString("NM");
    }

    private static GATKRead reduceContigReads(final GATKRead read1, final GATKRead read2) {
        if (read2.isUnmapped() || read2.isSecondaryAlignment()) {
            return read1;
        } else if (read1.isUnmapped() || read1.isSecondaryAlignment()) {
            return read2;
        } else if (!containsHardclips(read1)) {
            return read1;
        } else if (!containsHardclips(read2)) {
            return read2;
        } else {
            return mergeSequences(read1, read2);
        }
    }

    private static boolean containsHardclips(final GATKRead read1) {
        return read1.getCigar().getCigarElements().stream().anyMatch(e -> e.getOperator() == CigarOperator.HARD_CLIP);
    }

    private static GATKRead mergeSequences(final GATKRead read1, final GATKRead read2) {
        final int contigLength = contigLength(read1);
        final byte[] bases = new byte[contigLength];
        final MutableInt start = new MutableInt(contigLength);
        final MutableInt end = new MutableInt(-1);
        final SATagBuilder saBuilder1 = new SATagBuilder(read1);
        final SATagBuilder saBuilder2 = new SATagBuilder(read2);
        saBuilder1.removeTag(read2);
        saBuilder1.removeTag(read1);
        mergeSequences(bases, start, end, read1);
        mergeSequences(bases, start, end, read2);
        final byte[] mergedBases = Arrays.copyOfRange(bases, start.intValue(), end.intValue());
        final List<CigarElement> elements = new ArrayList<>();
        if (start.intValue() > 0) {
            elements.add(new CigarElement(start.intValue(), CigarOperator.H));
        }
        elements.add(new CigarElement(end.intValue() - start.intValue(), CigarOperator.M));
        if (end.intValue() < contigLength) {
            elements.add(new CigarElement(contigLength - end.intValue(), CigarOperator.H));
        }
        final Cigar mergedCigar = new Cigar(elements);

        final GATKRead result = new ContigMergeGATKRead(read1.getName(), read1.getContig(), read1.getStart(), mergedBases,
                mergedCigar, Math.max(read1.getMappingQuality(), read2.getMappingQuality()), bases.length, read1.getReadGroup(), read1.isSupplementaryAlignment());

        final SATagBuilder resultSABuilder = new SATagBuilder(result);

        resultSABuilder.addAllTags(saBuilder1);
        resultSABuilder.addAllTags(saBuilder2);
        resultSABuilder.setSATag();
        return result;
    }

    private static void mergeSequences(final byte[] bases, final MutableInt start, final MutableInt end, final GATKRead read) {
        final byte[] readBases = read.getBases();
        Cigar cigar = read.getCigar();
        if (read.isReverseStrand()) {
            SequenceUtil.reverseComplement(readBases, 0, readBases.length);
            cigar = CigarUtils.invertCigar(cigar);
        }
        int nextReadBase = 0;
        int nextBase = 0;
        int hardClipStart = 0; // unless any leading H found.
        int hardClipEnd = bases.length; // unless any tailing H found.
        for (final CigarElement element : cigar) {
            final CigarOperator operator = element.getOperator();
            final int length = element.getLength();
            if (operator == CigarOperator.H) {
                    if (nextBase == 0) { // hard-clip at the beginning.
                        hardClipStart = length;
                    } else { // hard-clip at the end.
                        hardClipEnd = bases.length - length;
                    }
                nextBase += length;
            } else if (operator.consumesReadBases()) {
                for (int i = 0; i < length; i++) {
                    bases[nextBase + i] = mergeBase(bases[nextBase + i], readBases[nextReadBase + i], () -> "mismatching bases");
                }
                nextBase += element.getLength();
                nextReadBase += element.getLength();
            }
        }
        if (hardClipStart < start.intValue()) {
            start.setValue(hardClipStart);
        }
        if (hardClipEnd > end.intValue()) {
            end.setValue(hardClipEnd);
        }
    }

    private static byte mergeQual(final byte a, final byte b) {
        return (byte) Math.max(a, b);
    }

    private static byte mergeBase(final byte a, final byte b, final Supplier<String> errorMessage) {
        if (a == 0) {
            return b;
        } else if (b == 0) {
            return a;
        } else if (a == b) {
            return a;
        } else {
            throw new IllegalStateException(errorMessage.get());
        }
    }


    private static int contigLength(final GATKRead contig) {
        return contig.getCigar().getCigarElements().stream()
                .filter(e -> e.getOperator() == CigarOperator.H || e.getOperator().consumesReadBases())
                .mapToInt(CigarElement::getLength)
                .sum();
    }
}

