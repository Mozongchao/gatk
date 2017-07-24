package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.evidence.AlignedAssemblyOrExcuse;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFHeaderLines;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;

/**
 * Created by valentin on 4/26/17.
 */
public class StructuralVariantContext extends VariantContext {

    private static final String ASSEMBLY_ID_PREFIX = "asm";
    private static final String CONTIG_ID_PREFIX = "tig";
    private static final char CONGIT_NAME_PART_SEPARATOR_CHR = ':';

    private static final long serialVersionUID = 1L;

    protected StructuralVariantContext(final VariantContext other) {
        super(other);
    }

    public boolean isStructural() {
        return getNAlleles() > 1 && getAlternateAlleles().stream()
                .anyMatch(StructuralVariantAllele::isStructural);
    }

    /**
     * Returns the assembly ids for this context's structural variant.
     * <p>
     *     The list returned is an immutable list.
     * </p>
     * @throws IllegalStateException if the {@link GATKSVVCFHeaderLines#CONTIG_NAMES} annotation contains
     * undefined contig names (.)
     *
     * @return never {@code null}, an empty list if no structural variant is specified.
     */
    public List<String> contigNames() {
        if (!hasAttribute(GATKSVVCFHeaderLines.CONTIG_NAMES)) {
            return Collections.emptyList();
        } else {
            final List<String> contigNames = getAttributeAsStringList(GATKSVVCFHeaderLines.CONTIG_NAMES, null);
            if (contigNames.contains(null)) {
                throw new IllegalStateException("the contig names annotation contains undefined values");
            }
            return contigNames;
        }
    }

    public List<String> assemblyIDs() {
        final List<String> contigNames = contigNames();
        return contigNames.stream().map(StructuralVariantContext::extractAssemblyIdFromContigName).collect(Collectors.toList());
    }

    private static String extractAssemblyIdFromContigName(final String s) {
        final int partSplitIndex = s.indexOf(CONGIT_NAME_PART_SEPARATOR_CHR);
        if (partSplitIndex == -1) {
            throw new UserException.BadInput(
                    String.format("found an %s annotation with and invalid contig name without parts separator %s", GATKSVVCFHeaderLines.CONTIG_NAMES, CONGIT_NAME_PART_SEPARATOR_CHR));
        } else {
            final String firstPart = s.substring(0, partSplitIndex);
            if (firstPart.startsWith(ASSEMBLY_ID_PREFIX)) {
                return firstPart;
            } else {
                throw new UserException.BadInput(
                        String.format("found an %s annotation with an invalid assembly id part '%s' starting with the wrong prefix (expects '%s' instead)", GATKSVVCFHeaderLines.CONTIG_NAMES,
                                firstPart, ASSEMBLY_ID_PREFIX));
            }
        }
    }

    public Haplotype composeHaplotype(final int index, final int paddingSize, final ReferenceMultiSource reference)  {
        if (index < 0 || index > 1)
            throw new IllegalArgumentException("wrong index must be 0 (ref) or 1 (alt)");
        Utils.nonNull(reference);

        final SimpleInterval referenceInterval = new SimpleInterval(getContig(), getStart() - paddingSize - 1, getStart() + paddingSize + Math.abs(Math.min(0,getStructuralVariantLength())));
        final ReferenceBases bases;
        try {
            bases = reference.getReferenceBases(null, referenceInterval);
        } catch (final IOException ex) {
            throw new UserException.CouldNotReadInputFile("could not read reference file");
        }
        if (index == 0) {
            final Haplotype result = new Haplotype(bases.getBases(), true);
            result.setCigar(new Cigar(Collections.singletonList(new CigarElement(referenceInterval.size(), CigarOperator.M))));
            result.setGenomeLocation(referenceInterval);
            return result;
        } else { //index == 1
            switch (getStructuralVariantAllele()) {
                case INS:
                    return composeInsertionHaplotype(bases);
                case DEL:
                    return composeDeletionHaplotype(bases);
                default: // not jet supported.
                    throw new UnsupportedOperationException("not supported yet");
            }
        }
    }

    private Haplotype composeDeletionHaplotype(final ReferenceBases referenceBases) {
        final int deletionSize = - getStructuralVariantLength();
        final byte[] resultBases = new byte[referenceBases.getInterval().size() - deletionSize];
        final int leftPaddingSize = getStart() - referenceBases.getInterval().getStart() + 1;
        final int rightPaddingSize = referenceBases.getInterval().getEnd() - getStart() - deletionSize;
        final byte[] referenceBaseBytes = referenceBases.getBases();
        System.arraycopy(referenceBaseBytes, 0, resultBases, 0, leftPaddingSize);
        System.arraycopy(referenceBaseBytes, leftPaddingSize + deletionSize, resultBases, leftPaddingSize, rightPaddingSize);
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(leftPaddingSize, CigarOperator.M),
                new CigarElement(deletionSize, CigarOperator.D),
                new CigarElement(rightPaddingSize, CigarOperator.M)));
        final Haplotype result = new Haplotype(resultBases, false);
        result.setCigar(cigar);
        result.setGenomeLocation(referenceBases.getInterval());
        return result;
    }

    private Haplotype composeInsertionHaplotype(final ReferenceBases referenceBases) {
        final byte[] insertedSequence = getInsertedSequence();
        final byte[] referenceBaseBytes = referenceBases.getBases();
        final byte[] resultBases = new byte[referenceBases.getInterval().size() + insertedSequence.length];
        final int leftPaddingSize = getStart() - referenceBases.getInterval().getStart() + 1;
        final int rightPaddingSize = referenceBases.getInterval().getEnd() - getStart();
        System.arraycopy(referenceBaseBytes, 0, resultBases, 0, leftPaddingSize);
        System.arraycopy(insertedSequence, 0, resultBases, leftPaddingSize, insertedSequence.length);
        System.arraycopy(referenceBaseBytes, rightPaddingSize, resultBases, leftPaddingSize + insertedSequence.length, rightPaddingSize);
        final Cigar cigar = new Cigar(Arrays.asList(new CigarElement(leftPaddingSize, CigarOperator.M),
                new CigarElement(insertedSequence.length, CigarOperator.I),
                new CigarElement(rightPaddingSize, CigarOperator.M)));
        final Haplotype result = new Haplotype(resultBases, false);
        result.setCigar(cigar);
        result.setGenomeLocation(referenceBases.getInterval());
        return result;
    }

    public StructuralVariantAllele getStructuralVariantAllele() {
        if (!isStructural()) {
            throw new IllegalStateException("this is not in fact a structural variant context");
        } else {
            return StructuralVariantAllele.valueOf(alleles.get(1));
        }
    }

    public byte[] getInsertedSequence() {
        if (hasAttribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE)) {
            final String asString = getAttributeAsString(GATKSVVCFHeaderLines.INSERTED_SEQUENCE, null);
            if (asString == null) {
                return null;
            } else {
                return asString.getBytes();
            }
        } else {
            return null;
        }
    }

    public int getStructuralVariantLength() {

        final Supplier<List<String>> s = ArrayList<String>::new;
        if (hasAttribute(GATKSVVCFHeaderLines.SVLEN)) {
            return getAttributeAsInt(GATKSVVCFHeaderLines.SVLEN, 0);
        } else {
            throw new IllegalStateException("missing insertion length");
        }
    }
}
