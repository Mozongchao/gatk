package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.api.services.genomics.model.Read;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by valentin on 7/17/17.
 */
class ContigMergeGATKRead implements GATKRead {

    private final Map<String, String> attributes;
    private final String contig;
    private final int length;
    private final String name;
    private final int start;
    private final byte[] bases;
    private final Cigar cigar;
    private final int mq;
    private final String group;
    private boolean supplementary;

    ContigMergeGATKRead(final String name, final String contig, final int start, final byte[] bases, final Cigar cigar, final int mq, final int length, final String readGroup, final boolean supplementary) {
        this.name = name;
        this.contig = contig;
        this.start = start;
        this.bases = bases;
        this.cigar = cigar;
        this.group = readGroup;
        this.supplementary = supplementary;
        this.attributes = new HashMap<>(1);
        this.mq = mq;
        this.length = length;
    }

    @Override
    public String getContig() {
        return contig;
    }

    @Override
    public int getStart() {
        return start;
    }

    @Override
    public int getEnd() {
        return start + bases.length;
    }

    @Override
    public String getName() {
        return name;
    }

    @Override
    public void setName(String name) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getLength() {
        return length;
    }

    @Override
    public void setPosition(String contig, int start) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setPosition(Locatable locatable) {
        throw new UnsupportedOperationException();
    }

    @Override
    public String getAssignedContig() {
        return contig;
    }

    @Override
    public int getAssignedStart() {
        return start;
    }

    @Override
    public int getUnclippedStart() {
        return start - (cigar.getCigarElements().size() == 1 ? 0 : cigar.getCigarElement(0).getLength());
    }

    @Override
    public int getUnclippedEnd() {
        return start + bases.length + (cigar.getCigarElements().size() == 1? 0 : cigar.getCigarElement(2).getLength());
    }

    @Override
    public String getMateContig() {
        return contig;
    }

    @Override
    public int getMateStart() {
        return start;
    }

    @Override
    public void setMatePosition(String contig, int start) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setMatePosition(Locatable locatable) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getFragmentLength() {
        return length;
    }

    @Override
    public void setFragmentLength(int fragmentLength) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getMappingQuality() {
        return mq;
    }

    @Override
    public void setMappingQuality(int mappingQuality) {
        throw new UnsupportedOperationException();
    }

    @Override
    public byte[] getBases() {
        return bases.clone();
    }

    @Override
    public void setBases(byte[] bases) {
        throw new UnsupportedOperationException();
    }

    @Override
    public byte[] getBaseQualities() {
        return new byte[0];
    }

    @Override
    public void setBaseQualities(byte[] baseQualities) {
        throw new UnsupportedOperationException();
    }

    @Override
    public Cigar getCigar() {
        return cigar;
    }

    @Override
    public void setCigar(Cigar cigar) {
        throw new UnsupportedOperationException();
    }

    @Override
    public void setCigar(String cigarString) {
        throw new UnsupportedOperationException();
    }

    @Override
    public String getReadGroup() {
        return group;
    }

    @Override
    public void setReadGroup(String readGroupID) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isPaired() {
        return false;
    }

    @Override
    public void setIsPaired(boolean isPaired) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isProperlyPaired() {
        return false;
    }

    @Override
    public void setIsProperlyPaired(boolean isProperlyPaired) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isUnmapped() {
        return false;
    }

    @Override
    public void setIsUnmapped() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean mateIsUnmapped() {
        return false;
    }

    @Override
    public void setMateIsUnmapped() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isReverseStrand() {
        return false;
    }

    @Override
    public void setIsReverseStrand(boolean isReverseStrand) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean mateIsReverseStrand() {
        return false;
    }

    @Override
    public void setMateIsReverseStrand(boolean mateIsReverseStrand) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isFirstOfPair() {
        return true;
    }

    @Override
    public void setIsFirstOfPair() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isSecondOfPair() {
        return true;
    }

    @Override
    public void setIsSecondOfPair() {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isSecondaryAlignment() {
        return false;
    }

    @Override
    public void setIsSecondaryAlignment(boolean isSecondaryAlignment) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isSupplementaryAlignment() {
        return supplementary;
    }

    @Override
    public void setIsSupplementaryAlignment(boolean isSupplementaryAlignment) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean failsVendorQualityCheck() {
        return false;
    }

    @Override
    public void setFailsVendorQualityCheck(boolean failsVendorQualityCheck) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean isDuplicate() {
        return false;
    }

    @Override
    public void setIsDuplicate(boolean isDuplicate) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean hasAttribute(final String attributeName) {
        return attributes.containsKey(attributeName);
    }

    @Override
    public Integer getAttributeAsInteger(final String attributeName) {
        final String str = attributes.get(attributeName);
        return Integer.valueOf(str);
    }

    @Override
    public String getAttributeAsString(final String attributeName) {
        return attributes.get(attributeName);
    }

    @Override
    public byte[] getAttributeAsByteArray(final String attributeName) {
        final String str = attributes.get(attributeName);
        return str == null ? null : str.getBytes();
    }

    @Override
    public void setAttribute(final String attributeName, final Integer attributeValue) {
        attributes.put(attributeName, "" + attributeValue);
    }

    @Override
    public void setAttribute(final String attributeName, final String attributeValue) {
        attributes.put(attributeName, attributeValue);
    }

    @Override
    public void setAttribute(String attributeName, byte[] attributeValue) {
        attributes.put(attributeName, new String(attributeValue));
    }

    @Override
    public void clearAttribute(String attributeName) {
        attributes.remove(attributeName);
    }

    @Override
    public void clearAttributes() {
        attributes.clear();
    }

    @Override
    public GATKRead copy() {
        return this;
    }

    @Override
    public GATKRead deepCopy() {
        return this;
    }

    @Override
    public SAMRecord convertToSAMRecord(SAMFileHeader header) {
        final SAMRecord record = new SAMRecord(header);
        record.setAlignmentStart(this.getAssignedStart());
        record.setReferenceName(this.getContig());
        record.setBaseQualities(this.getBaseQualities());
        record.setReadBases(this.getBases());
        record.setReadName(this.getName());
        record.setReadNegativeStrandFlag(this.isReverseStrand());
        record.setMappingQuality(this.getMappingQuality());
        if (this.hasAttribute("SA")) {
            record.setAttribute("SA", this.getAttributeAsString("SA"));
        }
        record.setCigar(this.getCigar());
        return record;
    }

    @Override
    public Read convertToGoogleGenomicsRead() {
        throw new UnsupportedOperationException();
    }

    @Override
    public String getSAMString() { // To much work to bother to provide an implemntation that won't ever be used.
        throw new UnsupportedOperationException();
    }
}
