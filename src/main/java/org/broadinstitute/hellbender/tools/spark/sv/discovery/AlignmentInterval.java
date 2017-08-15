package org.broadinstitute.hellbender.tools.spark.sv.discovery;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.tribble.annotation.Strand;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.MISSING_NM;

/**
 * Each assembled contig should have at least one such accompanying structure, or 0 when it is unmapped.
 */
@DefaultSerializer(AlignmentInterval.Serializer.class)
public final class AlignmentInterval {

    public static final int MISSING_NM = -1;
    public static final int MISSING_AS = -1;

    public final SimpleInterval referenceSpan;
    public final int startInAssembledContig;   // 1-based, inclusive
    public final int endInAssembledContig;     // 1-based, inclusive

    public final Cigar cigarAlong5to3DirectionOfContig;

    public final boolean forwardStrand;
    public final int mapQual;
    public final int mismatches;

    public final int alnScore;
    public final boolean isFromSplitGapAlignment;

    /**
     * Compose an alignment interval instance from a SAM supplementary alignment formatted string.
     * <p>
     *     The input string format is:
     *     <pre>
     *         chr-name,start,strand,cigar,mq,nm,as
     *     </pre>
     *     where mq (mapping-quality), nm (number of mismatches) and as (alignment score) might be absent. The strand
     *     is symbolized as '+' for the forward-strand and '-' for the reverse-strand.
     *     <p>Examples:</p>
     *     <pre>
     *         chr10,1241241,+,10S1313M45I14M100H,30,5,66
     *         chr10,1241241,+,10S1313M45I14M100H,30,5
     *         chr10,1241241,+,10S1313M45I14M100H,30
     *         chr10,1241241,+,10S1313M45I14M100H
     *     </pre>
     * </p>
     * @param str the input string.
     * @throws IllegalArgumentException if {@code str} is {@code null} or it does not look like a valid
     *   SA string.
     */
    public AlignmentInterval(final String str) {
        Utils.nonNull(str, "input str cannot be null");
        final String[] parts = str.replaceAll(";$", "").split(",");
        if (parts.length < 4) {
            throw new IllegalArgumentException("the input SA string at least must contain 4 parts");
        }
        final String referenceContig = parts[0];
        final int start = Integer.parseInt(parts[1]);
        final Strand strand = Strand.toStrand(parts[2]);
        if (strand == Strand.NONE) {
            throw new IllegalArgumentException("the input strand cannot be " + Strand.NONE);
        }
        final boolean forwardStrand = strand == Strand.POSITIVE;
        final Cigar originalCigar;
        try {
            originalCigar = TextCigarCodec.decode(parts[3]);
        } catch (final RuntimeException ex) {
            throw new IllegalArgumentException("bad formatted cigar " + parts[3]);
        }
        final Cigar cigar = forwardStrand ? originalCigar : CigarUtils.invertCigar(originalCigar);
        final int mappingQuality = parts.length >= 5 ? Integer.parseInt(parts[4]) : 0;
        final int mismatches = parts.length >= 6 ? Integer.parseInt(parts[5]) : -1;
        final int alignmentScore = parts.length >= 7 ? Integer.parseInt(parts[6]) : -1;
        this.referenceSpan = new SimpleInterval(referenceContig, start,
                Math.max(start, CigarUtils.referenceBasesConsumed(cigar) + start - 1));
        this.startInAssembledContig = 1 + CigarUtils.leftHardClippedBases(cigar);
        this.endInAssembledContig = CigarUtils.readLength(cigar)
                - CigarUtils.leftHardClippedBases(cigar);
        this.mapQual = mappingQuality;
        this.mismatches = mismatches;
        this.alnScore = alignmentScore;
        this.forwardStrand = forwardStrand;
        this.cigarAlong5to3DirectionOfContig = cigar;
        this.isFromSplitGapAlignment = false;
    }


    @VisibleForTesting
    public AlignmentInterval(final SAMRecord samRecord) {

        final boolean isMappedReverse = samRecord.getReadNegativeStrandFlag();
        this.referenceSpan = new SimpleInterval(samRecord);
        this.startInAssembledContig = getAlignmentStartInOriginalContig(samRecord);
        this.endInAssembledContig = getAlignmentEndInOriginalContig(samRecord);

        this.cigarAlong5to3DirectionOfContig = isMappedReverse ? CigarUtils.invertCigar(samRecord.getCigar())
                                                               : samRecord.getCigar();
        this.forwardStrand = !isMappedReverse;
        this.mapQual = samRecord.getMappingQuality();
        final Integer nMismatches = samRecord.getIntegerAttribute(SAMTag.NM.name());
        this.mismatches = nMismatches==null ? MISSING_NM : nMismatches;

        int alignerScore;
        try {
            alignerScore = samRecord.getIntegerAttribute(SAMTag.AS.name());
        } catch (final RuntimeException ex) {
            alignerScore = MISSING_AS;
        }
        this.alnScore = alignerScore;
        this.isFromSplitGapAlignment = false;
    }

    /**
     * Construct an alignment interval that reflects on the mapping properties of a {@link GATKRead} instance.
     * @param read
     */
    public AlignmentInterval(final GATKRead read) {
        final boolean isMappedReverse = read.isReverseStrand();
        this.referenceSpan = new SimpleInterval(read);
        this.startInAssembledContig = read.getFirstAlignedReadPosition();
        this.endInAssembledContig = read.getLastAlignedReadPosition();
        this.cigarAlong5to3DirectionOfContig = isMappedReverse ? CigarUtils.invertCigar(read.getCigar()) : read.getCigar();
        this.forwardStrand = !isMappedReverse;
        this.mapQual = read.getMappingQuality();
        this.mismatches = read.hasAttribute(SAMTag.NM.name()) ?  read.getAttributeAsInteger(SAMTag.NM.name()) : MISSING_NM;
        this.alnScore = read.hasAttribute(SAMTag.AS.name()) ?  read.getAttributeAsInteger(SAMTag.AS.name()) : MISSING_AS;
        this.isFromSplitGapAlignment = false;
    }

    @VisibleForTesting
    public AlignmentInterval(final BwaMemAlignment alignment, final List<String> refNames, final int unclippedContigLength) {

        // +1 because the BwaMemAlignment class has 0-based coordinate system
        this.referenceSpan = new SimpleInterval(refNames.get(alignment.getRefId()),
                                            alignment.getRefStart() + 1, alignment.getRefEnd());
        this.forwardStrand = 0==(alignment.getSamFlag() & SAMFlag.READ_REVERSE_STRAND.intValue());
        this.cigarAlong5to3DirectionOfContig = forwardStrand ? TextCigarCodec.decode(alignment.getCigar())
                                                             : CigarUtils.invertCigar(TextCigarCodec.decode(alignment.getCigar()));
        Utils.validateArg(
                cigarAlong5to3DirectionOfContig.getReadLength() + SvCigarUtils.getTotalHardClipping(cigarAlong5to3DirectionOfContig)
                        == unclippedContigLength,
                "contig length provided in constructor and inferred length by computation are different: " +
                        unclippedContigLength + "\t" + alignment.toString());

        // BwaMemAlignment has negative mapQ for unmapped sequences, not the same as its SAMRecord conversion
        // (see BwaMemAlignmentUtils.applyAlignment())
        this.mapQual = Math.max(SAMRecord.NO_MAPPING_QUALITY, alignment.getMapQual());
        this.mismatches = alignment.getNMismatches();
        if ( forwardStrand ) {
            this.startInAssembledContig = alignment.getSeqStart() + 1;
            this.endInAssembledContig = alignment.getSeqEnd();
        } else {
            this.startInAssembledContig = unclippedContigLength - alignment.getSeqEnd() + 1;
            this.endInAssembledContig = unclippedContigLength - alignment.getSeqStart();
        }

        this.alnScore = alignment.getAlignerScore();
        this.isFromSplitGapAlignment = false;
    }

    public AlignmentInterval(final SimpleInterval referenceSpan, final int startInAssembledContig, final int endInAssembledContig,
                             final Cigar cigarAlong5to3DirectionOfContig, final boolean forwardStrand,
                             final int mapQual, final int mismatches, final int alignerScore, final boolean isFromSplitGapAlignment) {
        this.referenceSpan = referenceSpan;
        this.startInAssembledContig = startInAssembledContig;
        this.endInAssembledContig = endInAssembledContig;

        this.cigarAlong5to3DirectionOfContig = cigarAlong5to3DirectionOfContig;

        this.forwardStrand = forwardStrand;
        this.mapQual = mapQual;
        this.mismatches = mismatches;
        this.alnScore = alignerScore;
        this.isFromSplitGapAlignment = isFromSplitGapAlignment;
    }

    /**
     * @return the number of bases of overlap between two alignment regions overlap on the locally-assembled contig they originate from.
     *          Mostly useful for computing micro-homologyForwardStrandRep.
     */
    @VisibleForTesting
    static int overlapOnContig(final AlignmentInterval one, final AlignmentInterval two) {
        return Math.max(0,
                Math.min(one.endInAssembledContig + 1, two.endInAssembledContig + 1)
                        - Math.max(one.startInAssembledContig, two.startInAssembledContig)
        );
    }

    static int getAlignmentStartInOriginalContig(final SAMRecord samRecord) {
        return SvCigarUtils.getNumClippedBases(!samRecord.getReadNegativeStrandFlag(), samRecord.getCigar()) + 1;
    }

    static int getAlignmentEndInOriginalContig(final SAMRecord samRecord) {
        final Cigar cigar = samRecord.getCigar();
        return cigar.getReadLength() + SvCigarUtils.getTotalHardClipping(cigar) -
                SvCigarUtils.getNumClippedBases(samRecord.getReadNegativeStrandFlag(), cigar);
    }

    AlignmentInterval(final Kryo kryo, final Input input) {
        final String chr   = input.readString();
        final int refStart = input.readInt(),
                  refEnd   = input.readInt();
        referenceSpan = new SimpleInterval(chr, refStart, refEnd);
        startInAssembledContig = input.readInt();
        endInAssembledContig = input.readInt();
        cigarAlong5to3DirectionOfContig = TextCigarCodec.decode(input.readString());
        forwardStrand = input.readBoolean();
        mapQual = input.readInt();
        mismatches = input.readInt();
        alnScore = input.readInt();
        isFromSplitGapAlignment = input.readBoolean();
    }

    void serialize(final Kryo kryo, final Output output) {
        output.writeString(referenceSpan.getContig());
        output.writeInt(referenceSpan.getStart());
        output.writeInt(referenceSpan.getEnd());
        output.writeInt(startInAssembledContig);
        output.writeInt(endInAssembledContig);
        output.writeString(TextCigarCodec.encode(cigarAlong5to3DirectionOfContig));
        output.writeBoolean(forwardStrand);
        output.writeInt(mapQual);
        output.writeInt(mismatches);
        output.writeInt(alnScore);
        output.writeBoolean(isFromSplitGapAlignment);
    }

    public static final class Serializer extends com.esotericsoftware.kryo.Serializer<AlignmentInterval> {
        @Override
        public void write( final Kryo kryo, final Output output, final AlignmentInterval alignmentInterval){
            alignmentInterval.serialize(kryo, output);
        }

        @Override
        public AlignmentInterval read(final Kryo kryo, final Input input, final Class<AlignmentInterval> clazz ) {
            return new AlignmentInterval(kryo, input);
        }
    }

    static final String PACKED_STRING_REP_SEPARATOR = "_";
    /**
     * @return  A packed String representation of this alignment interval; intended for debugging or annotation usage
     *          (both requires compactified message).
     */
    public String toPackedString() {
        return String.join(PACKED_STRING_REP_SEPARATOR, String.valueOf(startInAssembledContig),
                String.valueOf(endInAssembledContig), referenceSpan.toString(), (forwardStrand ? "+" : "-"),
                TextCigarCodec.encode(cigarAlong5to3DirectionOfContig),
                String.valueOf(mapQual), String.valueOf(mismatches), String.valueOf(alnScore),
                (isFromSplitGapAlignment ? "s" : "o"));
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        AlignmentInterval that = (AlignmentInterval) o;

        if (startInAssembledContig != that.startInAssembledContig) return false;
        if (endInAssembledContig != that.endInAssembledContig) return false;
        if (forwardStrand != that.forwardStrand) return false;
        if (mapQual != that.mapQual) return false;
        if (mismatches != that.mismatches) return false;
        if (alnScore != that.alnScore) return false;
        if (isFromSplitGapAlignment != that.isFromSplitGapAlignment) return false;
        if (!referenceSpan.equals(that.referenceSpan)) return false;
        return cigarAlong5to3DirectionOfContig.equals(that.cigarAlong5to3DirectionOfContig);
    }

    @Override
    public int hashCode() {
        int result = referenceSpan.hashCode();
        result = 31 * result + startInAssembledContig;
        result = 31 * result + endInAssembledContig;
        result = 31 * result + cigarAlong5to3DirectionOfContig.hashCode();
        result = 31 * result + (forwardStrand ? 1 : 0);
        result = 31 * result + mapQual;
        result = 31 * result + mismatches;
        result = 31 * result + alnScore;
        result = 31 * result + (isFromSplitGapAlignment ? 1 : 0);
        return result;
    }
    /**
     * Returns a {@link SAMRecord} instance that reflects this alignment interval given the
     * output {@link SAMFileHeader} and the enclosing {@link AlignedContig}.
     *
     * @param header the returned record header.
     * @param contig the enclosing contig.
     * @param hardClip whether clippings must be hard ones.
     * @return never {@code null}.
     */
    public SAMRecord convertToSAMRecord(final SAMFileHeader header, final AlignedContig contig, final boolean hardClip) {
        Utils.nonNull(header, "the input header cannot be null");
        Utils.nonNull(contig, "the input contig cannot be null");
        final SAMRecord result = new SAMRecord(header);

        result.setReadName(contig.contigName);
        result.setReadPairedFlag(false);
        result.setReadNegativeStrandFlag(!forwardStrand);

        // taking care of the bases;
        final byte[] bases = hardClip ? Arrays.copyOfRange(contig.contigSequence, startInAssembledContig - 1, endInAssembledContig) : contig.contigSequence.clone();
        if (!forwardStrand) {
            SequenceUtil.reverseComplement(bases);
        }
        result.setReadBases(bases);

        // taking care of the cigar.
        final Cigar cigar = forwardStrand ? this.cigarAlong5to3DirectionOfContig :
                CigarUtils.invertCigar(this.cigarAlong5to3DirectionOfContig);

        result.setCigar(hardClip ? CigarUtils.hardClip(cigar) : CigarUtils.softClip(cigar));

        result.setReferenceName(referenceSpan.getContig());
        result.setAlignmentStart(referenceSpan.getStart());
        if (mapQual >= 0) {
            result.setMappingQuality(this.mapQual);
        }
        if (mismatches != -1) {
            result.setAttribute(SAMTag.NM.name(), mismatches);
        }
        if (alnScore != -1) {
            result.setAttribute(SAMTag.AS.name(), alnScore);
        }
        return result;
    }
}
