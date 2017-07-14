package org.broadinstitute.hellbender.tools.spark.sv;

import com.sun.tools.internal.xjc.reader.xmlschema.bindinfo.BIConversion;
import htsjdk.samtools.util.Locatable;
import org.apache.hadoop.yarn.webapp.hamlet.Hamlet;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;
import org.broadinstitute.hellbender.utils.tsv.TableReader;
import org.broadinstitute.hellbender.utils.tsv.TableWriter;
import scala.Tuple2;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Aligned Contig Table utility class.
 */
class AlignedContigCoverageTable {

    static final String CONTIG_NAME_COLUMN_NAME = "CONTIG_NAME";

    static final String INTERVAL_COLUMN_NAME = "INTERVAL";

    static final String SUB_INTERVALS_COLUMN_NAME = "SUB_INTERVALS";

    static final String OTHER_INTERVALS_COLUMN_NAME = "OTHER_INTERVALS";

    static final Pattern START_END_PATTERN = Pattern.compile("^\\s*(\\d+)\\s*(\\-\\s*(\\d+))?\\s*$");

    static final TableColumnCollection COLUMNS = new TableColumnCollection(
            INTERVAL_COLUMN_NAME, CONTIG_NAME_COLUMN_NAME, SUB_INTERVALS_COLUMN_NAME,
            OTHER_INTERVALS_COLUMN_NAME);

    /**
     * Represents a record in the aligned-contig table.
     */
    static class Record implements Serializable, Locatable {

        private static long serialVersionID = 1L;

        public final String contigName;
        public final SimpleInterval interval;
        public final List<SimpleInterval> otherIntervals;
        public final List<SimpleInterval> subIntervals;

        private Record(final String contigName, final SimpleInterval interval, final List<SimpleInterval> subIntervals,
                       final List<SimpleInterval> otherIntervals) {
            this.contigName = contigName;
            this.interval = interval;
            this.otherIntervals = otherIntervals;
            this.subIntervals = subIntervals;
        }

        public static List<Record> fromAlignedContig(final AlignedContig ac) {
            final String contigName = ac.contigName;
            final Map<String, List<AlignedAssembly.AlignmentInterval>> byRefContig = ac.alignmentIntervals.stream()
                    .collect(Collectors.groupingBy(ai -> ai.referenceInterval.getContig()));
            final Map<String, Tuple2<SimpleInterval, List<SimpleInterval>>> byRefContigResult = new HashMap<>(byRefContig.size());
            for (final String contig : byRefContig.keySet()) {
                final List<AlignedAssembly.AlignmentInterval> alignedIntervals = byRefContig.get(contig);
                @SuppressWarnings("ConstantConditions")
                final int start = alignedIntervals.stream().mapToInt(ai -> ai.referenceInterval.getStart()).min().getAsInt();
                @SuppressWarnings("ConstantConditions")
                final int end = alignedIntervals.stream().mapToInt(ai -> ai.referenceInterval.getEnd()).max().getAsInt();
                final SimpleInterval interval = new SimpleInterval(contigName, start, end);
                byRefContigResult.put(contig,
                        new Tuple2<>(interval, byRefContig.get(contig).stream()
                                .map(ai -> ai.referenceInterval).collect(Collectors.toList())));
            }
            final List<SimpleInterval> allIntervals = byRefContigResult.values().stream().map(v -> v._1()).collect(Collectors.toList());

            return byRefContigResult.values().stream()
                    .map(v -> new Record(null, v._1(), v._2(),
                            allIntervals.stream().filter(interval -> !interval.equals(v._1())).collect(Collectors.toList())))
                    .collect(Collectors.toList());
        }

        @Override
        public String getContig() {
            return interval.getContig();
        }

        @Override
        public int getStart() {
            return interval.getStart();
        }

        @Override
        public int getEnd() {
            return interval.getEnd();
        }
    }

    /**
     * Table writer.
     */
    static class Writer extends TableWriter<Record> {
        public Writer(final java.io.Writer writer) throws IOException {
            super(writer, COLUMNS);
        }

        @Override
        protected void composeLine(final Record record, final DataLine dataLine) {
            dataLine.append(record.contigName, record.interval.toString());
            dataLine.append(record.subIntervals.stream().map(i -> "" + i.getStart() + "-" + i.getEnd()).collect(Collectors.joining(",")));
            dataLine.append(record.otherIntervals.stream().map(SimpleInterval::toString).collect(Collectors.joining(",")));
        }
    }

    /**
     * Table reader.
     */
    static class Reader extends TableReader<Record> {

        protected Reader(final String sourceName, final java.io.Reader sourceReader) throws IOException {
            super(sourceName, sourceReader);
        }

        @Override
        protected void processColumns(final TableColumnCollection columns) {
            if (!COLUMNS.matchesExactly(columns)) {
                throw formatException("bad header, expected column names : " + String.join(", ", COLUMNS.names()));
            }
        }

        @Override
        protected Record createRecord(final DataLine dataLine) {
            final String contigName = dataLine.nextString();
            final SimpleInterval interval;
            try {
                interval = new SimpleInterval(dataLine.nextString());
            } catch (final IllegalArgumentException ex) {
                throw formatException("bad interval string: " + dataLine.get(2));
            }
            final String subPartsString = dataLine.nextString();
            final String[] subPartsSplitString = subPartsString.split("\\s*,\\s*");
            final List<SimpleInterval> subIntervals = Stream.of(subPartsSplitString).map(p -> {
                final Matcher m = START_END_PATTERN.matcher(p);
                if (!m.find()) {
                    throw formatException("bad range string: " + p);
                } else if (m.group(3) != null || !m.group(3).isEmpty()) { // format is: start-end
                    try {
                        return new SimpleInterval(interval.getContig(), Integer.parseInt(m.group(1)), Integer.parseInt(m.group(3)));
                    } catch (final IllegalArgumentException ex) {
                        throw formatException("invalid range: " + p);
                    }
                } else { // format is: start
                    final int startEnd = Integer.parseInt(m.group(1));
                    return new SimpleInterval(interval.getContig(), startEnd, startEnd);
                }
            }).collect(Collectors.toList());
            final String otherIntervalsString = dataLine.nextString();
            final String[] otherPartsSplitString = otherIntervalsString.split("\\s*,\\s*");
            final List<SimpleInterval> otherIntervals = Stream.of(otherPartsSplitString).map(p -> {
                try {
                    return new SimpleInterval(p);
                } catch (final IllegalArgumentException ex) {
                    throw formatException("bad interval string: " + p);
                }
            }).collect(Collectors.toList());
            return new Record(contigName, interval, subIntervals, otherIntervals);
        }
    }

    static void write(final java.io.Writer where, final String whereName, final Stream<Record> records) {
        try (final Writer writer = new Writer(where)) {
            records.forEach(r -> {
                try {
                    writer.writeRecord(r);
                } catch (final IOException ex) {
                    throw new UncheckedIOException(ex);
                }
            });
        } catch (final Exception ex) {
            throw new UserException.CouldNotCreateOutputFile(whereName, ex.getMessage(), ex);
        }
    }
}
