/**
 * Copyright (c) 2016, Fulcrum Genomics LLC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

package com.fulcrumgenomics.metrics;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import com.fulcrumgenomics.util.SAMRecordDuplicateComparator;
import picard.cmdline.programgroups.Metrics;
import picard.util.IlluminaUtil;
import picard.util.RExecutor;

import java.io.File;
import java.util.*;

/** A tool to compute various unique molecular identifier metrics */
@CommandLineProgramProperties(
        usage = CollectUniqueMolecularIdentifierMetrics.USAGE,
        usageShort = CollectUniqueMolecularIdentifierMetrics.USAGE_SHORT,
        programGroup = Metrics.class
)
public class CollectUniqueMolecularIdentifierMetrics extends CommandLineProgram {

    public final static String USAGE_SHORT = "Calculates metrics for unique molecular identifiers.";

    public final static String USAGE = USAGE_SHORT + "\n" +
            "This tool produces a histogram of the number of reads per UMI. Additionally, this tool\n" +
            "\n" +
            "iterates through each set of duplicates and assigns a single UMI to the reasd in the\n" +
            "duplicate set in three ways:\n" +
            "  1. (FROM_DUPLICATE_SET) randomly choose a UMI based on the set of UMIs observed in a given duplicate set.\n" +
            "  2. (FROM_BACKGROUND) randomly choose a UMI based on the weighted global distribution of observed UMIs. The\n" +
            "     distribution is weighted on how frequently we see each UMI across all reads.\n" +
            "  3. (UNIFORM_RANDOM) randomly choose a UMI uniformly across all possible UMIs.\n" +
            "\n";

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "A SAM or BAM file to process.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The metrics file to write.")
    public File OUTPUT;

    @Option(shortName = "CHART", doc = "A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;

    @Option(doc = "Unique molecular identifier (UMI) SAM tag (ex. RX).")
    public String BARCODE_TAG = "RX";

    @Option(doc = "If multiple sequences are stored in the UMI SAM tag, use the ith one (0-based).", optional = true)
    public int BARCODE_TAG_INDEX = -1;

    @Option(doc = "Use the barcode tag when creating the duplicate sets.")
    public boolean USE_BARCODE_TAG = false;

    @Option(doc = "Assign the UMIs for each read independently, otherwise set the UMI for all reads in a duplicate set to the same UMI.")
    public boolean ASSIGN_UMI_PER_READ = true;

    @Option(doc = "When choosing a random UMI for FROM_DUPLICATE_SET (see usage), weight the distribution by occurrence\n" +
            " in the duplicate set, otherwise choose based on a uniform distribution across UMIs in the duplicate set.")
    public boolean FROM_DUPLICATE_SET_WEIGHTED = false;

    private final Log log = Log.getInstance(CollectUniqueMolecularIdentifierMetrics.class);

    private final Random random = new Random(42);

    private final String histogramBinLabel = "UMI_SEQUENCE";

    private int barcodeTagLength = -1;

    /** Used for generating random barcode tags and includes the delimiter */
    private byte[] barcodeTagTemplate = null;

    public final String R_SCRIPT = "com/fulcrumgenomics/metrics/uniqueMolecularIdentifierMetricsHistogram.R";

    @Override
    public int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(CHART_OUTPUT);

        final UniqueMolecularIdentifierMetrics metrics = new UniqueMolecularIdentifierMetrics();
        final Histogram<String> localRandomDistribution = new Histogram<String>(histogramBinLabel, "FROM_DUPLICATE_SET");
        final Histogram<String> globalRandomUmiDistribution = new Histogram<String>(histogramBinLabel, "FROM_BACKGROUND");
        final Histogram<String> randomDistribution = new Histogram<String>(histogramBinLabel, "UNIFORM_RANDOM");
        final Histogram<String> globalRandomUmiDistributionForReadOne = new Histogram<String>();

        log.info("Getting the global distribution of unique molecular identifiers.");
        final Histogram<String> globalUmiDistribution = getGlobalUniqueMolecularIdentifierDistribution();
        final int globalUmiDistributionSize = (int)globalUmiDistribution.getSumOfValues();

        log.info("Computing the metrics (sorting may take some time).");
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        // Get the appropriate comparator for building duplicate sets
        final htsjdk.samtools.SAMRecordDuplicateComparator comparator;
        if (USE_BARCODE_TAG) comparator = new SAMRecordDuplicateComparator(Collections.singletonList(reader.getFileHeader()), BARCODE_TAG, null, null);
        else comparator = new htsjdk.samtools.SAMRecordDuplicateComparator(Collections.singletonList(reader.getFileHeader()));

        final DuplicateSetIterator iterator = new DuplicateSetIterator(getIterator(reader), reader.getFileHeader(), false, comparator);
        final ProgressLogger progress = new ProgressLogger(log, 10000);
        for (final DuplicateSet duplicateSet : new IterableAdapter<DuplicateSet>(iterator)) {
            final SAMRecord representativeRecord = duplicateSet.getRepresentative();

            // We assume that the barcode tag is per paired-end read, so we care only about one end of a pair.  We choose
            // arbitrarily to look at when the representative is read one.
            if (representativeRecord.getReadPairedFlag() && representativeRecord.getFirstOfPairFlag() && !representativeRecord.isSecondaryOrSupplementary()) {
                final List<SAMRecord> records = new ArrayList<SAMRecord>(duplicateSet.size());
                final Set<String> barcodeTags = new HashSet<String>(duplicateSet.size());

                // do not include secondary or supplementary
                for (final SAMRecord record : duplicateSet.getRecords()) {
                    if (record.getReadPairedFlag() && !record.isSecondaryOrSupplementary()) {
                        records.add(record);
                        if (!FROM_DUPLICATE_SET_WEIGHTED) barcodeTags.add(record.getStringAttribute(BARCODE_TAG));
                    }
                }

                final int numRecords = records.size();
                // NB: could use getTransientAttribute below, but that may be a little brittle only to improve performance.

                if (0 < numRecords) {
                    if (ASSIGN_UMI_PER_READ) { // a new UMI for each read
                        for (final SAMRecord record : records) {
                            // 1. choose a UMI randomly from the set of records in this duplicate set
                            final String localRandomUmi = getLocalRandomUmi(records, barcodeTags);
                            localRandomDistribution.increment(localRandomUmi);

                            // 2. choose a UMI randomly from the global distribution of UMIs for this duplicate set
                            final String globalRandomUmi = getGlobalRandomUmi(globalUmiDistribution, globalUmiDistributionSize);
                            globalRandomUmiDistribution.increment(globalRandomUmi);

                            // 3. for each read in the duplicate set, choose a random UMI
                            randomDistribution.increment(getRandomBarcodeTag());
                        }
                    }
                    else { // same UMI for all reads in the duplicate set
                        // 1. choose a UMI randomly from the set of records in this duplicate set
                        final String localRandomUmi = getLocalRandomUmi(records, barcodeTags);
                        localRandomDistribution.increment(localRandomUmi, numRecords);

                        // 2. choose a UMI randomly from the global distribution of UMIs for this duplicate set
                        final String globalRandomUmi = getGlobalRandomUmi(globalUmiDistribution, globalUmiDistributionSize);
                        globalRandomUmiDistribution.increment(globalRandomUmi, numRecords);

                        // 3. choose a random UMI for the this duplicate set
                        randomDistribution.increment(getRandomBarcodeTag(), numRecords);
                    }
                }
            }

            // update the progress and read pair metrics
            for (final SAMRecord record : duplicateSet.getRecords()) {
                if (record.getReadPairedFlag() && record.getFirstOfPairFlag() && !record.isSecondaryOrSupplementary()) {
                    final String umi = record.getStringAttribute(BARCODE_TAG);
                    globalRandomUmiDistributionForReadOne.increment(umi);
                }
                progress.record(record);
            }
        }
        CloserUtil.close(reader);

        // update the metrics
        metrics.TOTAL_READ_PAIRS              = (long)globalRandomUmiDistributionForReadOne.getSumOfValues();
        metrics.NUM_OBSERVED_UMIS             = globalRandomUmiDistributionForReadOne.size();
        metrics.MEAN_READS_PER_OBSERVED_UMI   = globalRandomUmiDistributionForReadOne.getMeanBinSize();
        metrics.MEDIAN_READS_PER_OBSERVED_UMI = globalRandomUmiDistributionForReadOne.getMedianBinSize();
        metrics.SD_READS_PER_OBSERVED_UMI     = globalRandomUmiDistributionForReadOne.getStandardDeviationBinSize(metrics.MEAN_READS_PER_OBSERVED_UMI);
        metrics.MEAN_READS_PER_POSSIBLE_UMI   = dividByFourRepeatedly(metrics.TOTAL_READ_PAIRS, barcodeTagLength);
        metrics.MEAN_READS_PER_UMI_RATIO      = metrics.MEAN_READS_PER_POSSIBLE_UMI / metrics.MEAN_READS_PER_OBSERVED_UMI;

        // write the metrics file
        final MetricsFile<UniqueMolecularIdentifierMetrics, String> metricsFile = getMetricsFile();
        metricsFile.addMetric(metrics);
        metricsFile.addHistogram(globalUmiDistribution);
        metricsFile.addHistogram(localRandomDistribution);
        metricsFile.addHistogram(globalRandomUmiDistribution);
        metricsFile.addHistogram(randomDistribution);
        metricsFile.write(OUTPUT);

        // plot the histogram
        final int rResult = RExecutor.executeFromClasspath(R_SCRIPT,
                OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                INPUT.getName());
        if (rResult != 0) {
            throw new PicardException(String.format("R script %s failed with return code %d", R_SCRIPT, rResult));
        }

        return 0;
    }

    /** Metrics for unique molecular identifiers. */
    public class UniqueMolecularIdentifierMetrics extends MetricBase {

        /** The total number of read pairs with unique molecular identifiers. */
        public long TOTAL_READ_PAIRS;

        /** The total number of unique molecular identifiers observed in read pairs. */
        public long NUM_OBSERVED_UMIS;

        /** The mean number of read pairs per observed UMI. */
        public double MEAN_READS_PER_OBSERVED_UMI;

        /** The median number of read pairs per observed UMI. */
        public double MEDIAN_READS_PER_OBSERVED_UMI;

        /** The standard deviation of the number of reads per observed UMI. */
        public double SD_READS_PER_OBSERVED_UMI;

        /** The mean number of read pairs per possible UMI.  In this case, we assume
         * a fixed-length UMI with all 4^n possible UMI sequences. */
        public double MEAN_READS_PER_POSSIBLE_UMI;

        /** The ratio between MEAN_READS_PER_POSSIBLE_UMI and MEAN_READS_PER_OBSERVED_UMI, namely
         * MEAN_READS_PER_POSSIBLE_UMI / MEAN_READS_PER_OBSERVED_UMI.  Values closer to 1.0 indicate that the UMIs are
         * evenly distributed across the reads.  Values closer to 0.0 indicate that either the space of UMIs is
         * low and the read count is high (ex. 4bp UMI with 1M reads), or that a small fraction of UMIs are present
         * in the majority of reads.
         */
        public double MEAN_READS_PER_UMI_RATIO;
    }

    /** Gets an iterator of the reader, depending on if we want to split the barcode tag or not */
    private CloseableIterator<SAMRecord> getIterator(final SamReader reader) {
        if (BARCODE_TAG_INDEX < 0) {
            return reader.iterator();
        }
        else {
            return new UpdateBarcodeTagSAMRecordIterator(reader.iterator(), BARCODE_TAG, BARCODE_TAG_INDEX);
        }
    }

    /** Get a random barcode tag. */
    private String getRandomBarcodeTag() {
        final byte[] barcodeTag = new byte[barcodeTagTemplate.length];
        for (int i = 0; i < barcodeTagTemplate.length; i++) {
            if (barcodeTagTemplate[i] == SequenceUtil.N) {
                barcodeTag[i] = SequenceUtil.VALID_BASES_UPPER[random.nextInt(4)];
            }
            else {
                barcodeTag[i] = (byte)IlluminaUtil.BARCODE_DELIMITER.charAt(0);
            }
        }
        return new String(barcodeTag);
    }

    /** Is this the real life, or just paranoia? */
    private double dividByFourRepeatedly(final long numerator, final int numTimes) {
        double value = numerator;
        for (int i = 0; i < numTimes; i++) {
            value /= 4.0;
        }
        return value;
    }

    /** Gets a random UMI tag randomly based on the distribution of tags in barcodeTags */
    private String getLocalRandomUmi(final List<SAMRecord> records, final Set<String> barcodeTags) {
        final String localRandomUmi;
        if (FROM_DUPLICATE_SET_WEIGHTED) {
            localRandomUmi = records.get(random.nextInt(records.size())).getStringAttribute(BARCODE_TAG);
        }
        else {
            // choose a UMI randomly from set of UMIs found in the duplicate set
            int numLeft = random.nextInt(barcodeTags.size());
            final Iterator<String> barcodeTagsIterator = barcodeTags.iterator();
            while (0 < numLeft) { // skip over
                barcodeTagsIterator.next();
                numLeft--;
            }
            localRandomUmi = barcodeTagsIterator.next();
        }
        return localRandomUmi;
    }

    /** Gets a random UMI tag randomly based on the global observed distribution of UMIs */
    private String getGlobalRandomUmi(final Histogram<String> globalUmiDistribution, final int globalUmiDistributionSize) {
        final int r = random.nextInt(globalUmiDistributionSize);
        int sum = 0;
        for (final Histogram<String>.Bin bin : globalUmiDistribution.values()) {
            sum += (int)bin.getValue();
            if (r < sum) return bin.getId();
        }
        throw new PicardException("BUG: Could not get a random UMI");
    }

    /** Gets the distribution of UMIs from read one of primary alignment read pairs. */
    private Histogram<String> getGlobalUniqueMolecularIdentifierDistribution() {
        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        final Histogram<String> histogram = new Histogram<String>(histogramBinLabel, "OBSERVED");

        // Get the appropriate comparator for building duplicate sets
        // TODO: move this to it's own method
        final htsjdk.samtools.SAMRecordDuplicateComparator comparator;
        if (USE_BARCODE_TAG) comparator = new SAMRecordDuplicateComparator(Collections.singletonList(reader.getFileHeader()), BARCODE_TAG, null, null);
        else comparator = new htsjdk.samtools.SAMRecordDuplicateComparator(Collections.singletonList(reader.getFileHeader()));

        final DuplicateSetIterator iterator = new DuplicateSetIterator(getIterator(reader), reader.getFileHeader(), false, comparator);
        final ProgressLogger progress = new ProgressLogger(log, 10000);
        for (final DuplicateSet duplicateSet : new IterableAdapter<DuplicateSet>(iterator)) {
            final SAMRecord representativeRecord = duplicateSet.getRepresentative();
            if (representativeRecord.getReadPairedFlag() && representativeRecord.getFirstOfPairFlag() && !representativeRecord.isSecondaryOrSupplementary()) {
                for (final SAMRecord record : duplicateSet.getRecords()) {
                    if (record.getReadPairedFlag() && !record.isSecondaryOrSupplementary()) {
                        final String umi = record.getStringAttribute(BARCODE_TAG);
                        if (barcodeTagLength == -1) {
                            barcodeTagLength = 0;
                            barcodeTagTemplate = new byte[umi.length()];
                            int i = 0;
                            for (final byte base : umi.getBytes()) {
                                if (SequenceUtil.isValidBase(base) || SequenceUtil.isNoCall(base)) {
                                    barcodeTagLength++;
                                    barcodeTagTemplate[i] = SequenceUtil.N;
                                } else {
                                    barcodeTagTemplate[i] = base;
                                }
                                i++;
                            }
                        }
                        histogram.increment(umi);
                    }
                }
            }
            for (final SAMRecord record : duplicateSet.getRecords()) {
                progress.record(record);
            }
        }

        CloserUtil.close(iterator);
        return histogram;
    }

    /**
     * For the given barcode SAM tag, split the tag by the '-' delimiter, and updated the tag with the 0-based token.
     */
    class UpdateBarcodeTagSAMRecordIterator implements CloseableIterator<SAMRecord> {
        final private String barcodeTag;
        final private int tokenIndex;
        final private CloseableIterator<SAMRecord> iterator;

        public UpdateBarcodeTagSAMRecordIterator(final CloseableIterator<SAMRecord> iterator, final String barcodeTag, final int tokenIndex) {
            this.barcodeTag = barcodeTag;
            this.tokenIndex = tokenIndex;
            this.iterator   = iterator;
        }

        public boolean hasNext() { return this.iterator.hasNext(); }

        public SAMRecord next() {
            final SAMRecord record = this.iterator.next();
            final String tagValue = record.getStringAttribute(barcodeTag);
            final String[] tokens = tagValue.split(IlluminaUtil.BARCODE_DELIMITER);
            if (tokens.length <= tokenIndex) {
                throw new PicardException(String.format("Not enough tokens (%d) for SAM tag (%s) in record: %s", tokens.length, barcodeTag, record.getSAMString()));
            }
            record.setAttribute(barcodeTag, tokens[tokenIndex]);
            return record;
        }

        public void close() { this.iterator.close(); }
    }
}
