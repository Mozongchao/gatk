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
import javafx.collections.transformation.SortedList;
import org.apache.commons.lang3.mutable.MutableInt;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
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
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.GoogleGenomicsReadToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.ojalgo.function.BinaryFunction;
import org.seqdoop.hadoop_bam.BAMInputFormat;
import scala.Tuple2;

import javax.annotation.Nullable;
import java.io.Serializable;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
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
    public static final int DEFAULT_PADDING_SIZE = 300;


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
        processVariants(variantOverlappingContigs, alignedContigs);

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

        final Broadcast<IntervalsSkipList<SimpleInterval>> shardIntervalsBroadcast = ctx.broadcast(shardIntervals);

        final JavaPairRDD<SimpleInterval, List<GATKRead>> contigsInShards =
            groupInShards(contigs, shardIntervalsBroadcast);
        final JavaPairRDD<SimpleInterval, List<VariantContext>> variantsInShards =
            groupInShards(variants,
                    (Function<VariantContext, VariantContextVariantAdapter> & Serializable) VariantContextVariantAdapter::new,
                    (Function<VariantContextVariantAdapter, VariantContext> & Serializable) VariantContextVariantAdapter::getVariantContext,
                    shardIntervalsBroadcast);

        final JavaPairRDD<SimpleInterval, Tuple2<List<GATKRead>, List<VariantContext>>> contigAndVariantsInShards =
                contigsInShards.join(variantsInShards);

        final JavaPairRDD<VariantContext, List<GATKRead>> contigsPerVariant =
                contigAndVariantsInShards.flatMapToPair(t -> {
                    final List<VariantContext> vars = t._2()._2();
                    final List<GATKRead> ctgs = t._2()._1();
                    return vars.stream()
                            .map(v -> {
                                final List<GATKRead> cs = ctgs.stream()
                                        .filter(c -> c.getStart() <= v.getEnd() && c.getEnd() >= v.getStart())
                                .collect(Collectors.toList());
                                return new Tuple2<>(v, cs);})
                            .collect(Collectors.toList()).iterator();
                });
        return contigsPerVariant;
    }

    private static SimpleInterval locatableToSimpleInterval(final Locatable loc) {
        return new SimpleInterval(loc.getContig(), loc.getStart(), loc.getEnd());
    }

    private static JavaRDD<GATKRead> readsByInterval(final SimpleInterval interval, final ReadsSparkSource alignedContigs, final String alignedContigFileName, final String referenceName) {

        return alignedContigs.getParallelReads(alignedContigFileName,
                referenceName, Collections.singletonList(interval));
    }

    private <T, U extends Locatable> JavaPairRDD<SimpleInterval, List<T>> groupInShards(final JavaRDD<T> elements,
                                                                                        final Function<T, U> locatableOf,
                                                                                        final Function<U, T> elementOf,
                                                                                        final Broadcast<IntervalsSkipList<SimpleInterval>> shards) {
        return groupInShards(elements.map(e -> locatableOf.apply(e)), shards)
                .mapValues(l -> l.stream().map(elementOf).collect(Collectors.toList()));
    }

    private <T extends Locatable> JavaPairRDD<SimpleInterval, List<T>> groupInShards(final JavaRDD<T> elements, final Broadcast<IntervalsSkipList<SimpleInterval>> shards) {
        return elements
                .filter(c -> c.getContig() != null && c.getStart() > 0 && c.getEnd() >= c.getStart())
                .mapToPair(c -> new Tuple2<>(c, new SimpleInterval(c.getContig(), c.getStart(), c.getEnd())))
                .mapToPair(t -> new Tuple2<>(t._1(), shards.getValue().getOverlapping(t._2())))
                .flatMapValues(l -> l)
                .mapToPair(t -> new Tuple2<>(t._2(), t._1()))
                .aggregateByKey(new ArrayList<T>(10),
                        (l1, c) -> { l1.add(c); return l1;},
                        (l1, l2) -> {l1.addAll(l2); return l1;});

    }


    protected void processVariants(final JavaPairRDD<VariantContext, List<GATKRead>> variantsAndOverlappingContigs, final ReadsSparkSource s) {

        final BinaryOperator<GATKRead> readMerger = (BinaryOperator<GATKRead> & Serializable)
                ComposeStructuralVariantHaplotypesSpark::reduceContigReads;

        final JavaPairRDD<VariantContext, List<GATKRead>> variantsAndOverlappingUniqueContigs
                = variantsAndOverlappingContigs
                .mapValues(l -> l.stream().collect(Collectors.groupingBy(GATKRead::getName)))
                .mapValues(m -> m.values().stream()
                        .map(l -> l.stream().reduce(readMerger).orElseThrow(IllegalStateException::new))
                        .collect(Collectors.toList()));

        variantsAndOverlappingUniqueContigs.toLocalIterator()
                .forEachRemaining(vc -> {
                    logger.info("VC " + vc._1().getContig() + ":" + vc._1().getStart() + "-" + vc._1().getEnd());
                    vc._2().forEach(r -> {
                        if (!r.getCigar().containsOperator(CigarOperator.H)) {
                            logger.info("Contig " + r.getName() + " " + r.getLength() + " bases");
                        } else {
                            final SATagBuilder saTagBuilder = new SATagBuilder(r);
                            final String targetName = r.getName();
                            final List<SimpleInterval> locations = saTagBuilder.getArtificialReadsBasedOnSATag(s.getHeader(alignedContigsFileName, referenceArguments.getReferenceFileName()))
                                    .stream().map(rr -> new SimpleInterval(rr.getContig(), rr.getStart(), rr.getEnd()))
                                    .collect(Collectors.toList());
                            final GATKRead candidate = s.getParallelReads(alignedContigsFileName, referenceArguments.getReferenceFileName(), locations)
                                    .filter(rr -> rr.getName().equals(targetName))
                                    .reduce((a,b) -> readMerger.apply(a, b));

                            final GATKRead canonicRead = readMerger.apply(candidate, r);
                            if (!canonicRead.getCigar().containsOperator(CigarOperator.H)) {
                                logger.info("Contig " + canonicRead.getName() + " " + canonicRead.getLength() + " bases, needed additional search ");
                            } else {
                                logger.info("Contig " + canonicRead.getName() + " gave-up!");
                            }
                        }
                    });
                    logger.info("/VC");

                });
    }

    private GATKRead obtainCannonicAlignment(final GATKRead contig, final Function<SimpleInterval, JavaRDD<GATKRead>> overlappingReads) {
        if (!contig.getCigar().containsOperator(CigarOperator.H)) {
            return contig;
        } else {
            SATagBuilder builder = new SATagBuilder(contig);
            return builder.getArtificialReadsBasedOnSATag(getHeaderForReads()).stream()
                    .map(r -> new SimpleInterval(r.getContig(), r.getStart(), r.getEnd()))
                    .map(i -> overlappingReads.apply(i))
                    .map(rdd -> rdd.filter(r -> r.getName().equals(contig.getName())).first())
                    .filter(Objects::nonNull)
                    .filter(r -> r.getCigar().containsOperator(CigarOperator.H))
                    .findFirst().orElse(null);
        }
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

        resultSABuilder.addTag(saBuilder1);
        resultSABuilder.addTag(saBuilder2);
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

