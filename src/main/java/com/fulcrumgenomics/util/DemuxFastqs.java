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

package com.fulcrumgenomics.util;

import com.fulcrumgenomics.cmdline.Utilities;
import com.fulcrumgenomics.util.miseq.Sample;
import com.fulcrumgenomics.util.miseq.SampleBarcode;
import com.fulcrumgenomics.util.miseq.SampleSheet;
import htsjdk.samtools.*;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.illumina.ExtractIlluminaBarcodes;
import picard.illumina.ExtractIlluminaBarcodes.BarcodeMetric;
import picard.illumina.parser.ReadStructure;
import picard.util.IlluminaUtil;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * This tool demultiplexes a set of FASTQs based on the given Illumina Experiment Manager Sample Sheet for dual-indexed
 * sequencing runs.
 *
 * See the USAGE for a detailed description.
 *
 * Possible Future Improvements:
 * - adapter trimming (see {@link IlluminaUtil.IlluminaAdapterPair}).
 * - more metadata stored in the SampleSheet.csv
 */
@CommandLineProgramProperties(
        usage = DemuxFastqs.USAGE,
        usageShort = "Demultiplexes FASTQs based on sample barcodes and annotates molecular barcode information.",
        programGroup = Utilities.class
)
public class DemuxFastqs extends CommandLineProgram {

    public final static String USAGE =
            "This tool demultiplexes a set of FASTQs based on the given Illumina Experiment Manager Sample Sheet for dual-indexed\n" +
            "sequencing runs.\n" +
            "\n" +
            "Fastqs and read structures for read one, read two, i7 read, and i5 read should be given.  The read structures\n" +
            "may contain sample barcodes bases ('B'), molecular identifier bases ('M'), and template bases ('T'), with the latter\n" +
            "only for template reads (read one or read two).  Any molecular identifiers will be concatenated using the '-'\n" +
            "delimiter and placed in the given SAM record tag (\"RX\" by default).  The order of concatenation will be read one,\n" +
            "read two, i7 read, and i5 read, only considering reads with molecular identifiers.  Similarly, the sample\n" +
            "barcode bases from the given read will be placed in the \"BC\" tag, using the same rules as molecular identifiers,.\n" +
            "but applied to sample barcodes.\n" +
            "\n" +
            "The output directory will contain one BAM file per sample in the sample sheet, plus a BAM for reads that could\n" +
            "not be assigned to a sample given the criteria.  The output files will be the concatenation of sample id, sample\n" +
            "name, and sample barcode bases (expected not observed), delimited by \".\".  A metrics file will also be output\n" +
            "for sample barcodes.  More information about these metrics can be found here:\n" +
            "https://broadinstitute.github.io/picard/picard-metric-definitions.html#ExtractIlluminaBarcodes.BarcodeMetric\n" +
            "\n" +
            "The read group's sample id, sample name, and library id all correspond to the similarly named values in the sample sheet.\n" +
            "Library id will be the sample id if not found, and the platform unit will be the sample name concatenated with the sample\n" +
            "barcode bases delimited by a \".\".\n" +
            "\n" +
            "The sample section of the sample sheet should contain information related to each sample with the following keys and values:\n" +
            "  - Sample Identifier: " + SampleSheet.SAMPLE_ID + "\n" +
            "  - Sample Name: " + SampleSheet.SAMPLE_NAME + "\n" +
            "  - Library Identifier: " + SampleSheet.LIBRARY_ID + "\n" +
            "  - Sample Project: " + SampleSheet.SAMPLE_PROJCET + "\n" +
            "  - Description: : " + SampleSheet.DESCRIPTION + "\n" +
            "The following are optional values may included in the sample sheet to specify sample barcodes in any read: \n" +
            "  - Read One Inline Sample Barcode Bases: " + SampleSheet.R1_BARCODE_BASES + "\n" +
            "  - Read One Inline Sample Barcode Bases: " + SampleSheet.R2_BARCODE_BASES + "\n" +
            "  - i7 Sample Barcode Bases: " + SampleSheet.I7_BASES + "\n" +
            "  - i5 Sample Barcode Bases: " + SampleSheet.I5_BASES + "\n" +
            "The i7 and i5 Sample Barcode Bases may include extra bases before or after the sample barcode (if present\n" +
            "and specified by the read structure).\n" +
            "\n" +
            "The read structures will be used to extract the observed sample barcode and molecular identifiers from each\n" +
            "read.  The observed sample barcode will be matched to the sample barcodes extracted from the bases in the sample sheet\n" +
            "and associated read structures.\n" +
            "\n" +
            ReadStructure.PARAMETER_DOC;

    @Option(shortName="R1", doc="Input fastq file (optionally gzipped) for the first read of paired end data.")
    public List<File> READ_ONE_FASTQ;

    @Option(shortName="R2", doc="Input fastq file (optionally gzipped) for the second read of paired end data.")
    public List<File> READ_TWO_FASTQ;

    @Option(shortName="I7", doc="Input fastq file (optionally gzipped) for the index read of the Illumina i7 sequencing primer.  This is typically the I1 FASTQ (index read one for dual-indexed data).")
    public List<File> I7_FASTQ;

    @Option(shortName="I5", doc="Input fastq file (optionally gzipped) for the index read of the Illumina i5 sequencing primer. This is typically the I2 FASTQ (index read two for dual-indexed data).")
    public List<File> I5_FASTQ;

    @Option(shortName = "RS1", doc = "Read structure for the first read of paired end data.  Set to \"1000T\" or some large value to have all bases found be template bases. See the DemuxFastqs help message for more details.")
    public String READ_ONE_READ_STRUCTURE = "1000T";

    @Option(shortName = "RS2", doc = "Read structure for the first read of paired end data.  Set to \"1000T\" or some large value to have all bases found be template bases. See the DemuxFastqs help message for more details.")
    public String READ_TWO_READ_STRUCTURE = "1000T";

    @Option(shortName = "RS7", doc = "Read structure for the index read of the Illumina i7 sequencing primer.  The total bases must match the length of the index read. See the DemuxFastqs help message for more details.")
    public String I7_READ_STRUCTURE;

    @Option(shortName = "RS5", doc = "Read structure for the index read of the Illumina i5 sequencing primer.  The total bases must match the length of the index read.See the DemuxFastqs help message for more details.")
    public String I5_READ_STRUCTURE;

    @Option(shortName = "SS", doc = "The Sample Sheet (SampleSheet.csv).  This must include information about the i7/i5 index reads' ids and bases.")
    public File SAMPLE_SHEET;

    @Option(doc="Output directory. ", shortName= StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Option(doc = "Per-barcode and per-lane metrics written to this file.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME)
    public File METRICS_FILE;

    @Option(doc="Output file for the unmatched records.")
    public String UNMATCHED_OUTPUT = "unmatched.bam";

    @Option(doc="The SAM tag for the molecular barcode.")
    public String MOLECULAR_BARCODE_TAG = "RX";

    @Option(shortName="V", doc="A value describing how the quality values are encoded in the fastq.  Either Solexa for pre-pipeline 1.3 " +
            "style scores (solexa scaling + 66), Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard for phred scaled " +
            "scores with a character shift of 33.  If this value is not specified, the quality format will be detected automatically.", optional = true)
    public FastqQualityFormat QUALITY_FORMAT;

    @Option(shortName="PL", doc="The platform type (e.g. illumina, solid) to insert into the read group header", optional=true)
    public String PLATFORM = "ILLUMINA";

    @Option(shortName="CN", doc="The sequencing center from which the data originated", optional=true)
    public String SEQUENCING_CENTER;

    @Option(shortName = "PI", doc = "Predicted median insert size, to insert into the read group header", optional = true)
    public Integer PREDICTED_INSERT_SIZE;

    @Option(shortName = "PG", doc = "Program group to insert into the read group header.", optional=true)
    public String PROGRAM_GROUP;

    @Option(shortName = "PM", doc = "Platform model to insert into the group header (free-form text providing further details of the platform/technology used)", optional=true)
    public String PLATFORM_MODEL;

    @Option(doc="Comment(s) to include in the merged output file's header.", optional=true, shortName="CO")
    public List<String> COMMENT = new ArrayList<>();

    @Option(shortName = "DS", doc = "Inserted into the read group header", optional = true)
    public String DESCRIPTION;

    @Option(shortName = "DT", doc = "Date the run was produced, to insert into the read group header", optional = true)
    public Iso8601Date RUN_DATE;

    @Option(shortName="SO", doc="The sort order for the output sam/bam file.")
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.queryname;

    @Option(doc="Minimum quality allowed in the input fastq.  An exception will be thrown if a quality is less than this value.")
    public int MIN_Q = 0;

    @Option(doc="Maximum quality allowed in the input fastq.  An exception will be thrown if a quality is greater than this value.")
    public int MAX_Q = SAMUtils.MAX_PHRED_SCORE;

    @Option(doc="If true and this is an unpaired fastq any occurance of '/1' will be removed from the end of a read name.")
    public Boolean STRIP_UNPAIRED_MATE_NUMBER = false;

    @Option(doc="Allow (and ignore) empty lines")
    public Boolean ALLOW_AND_IGNORE_EMPTY_LINES = false;

    @Option(doc = "Maximum mismatches for a barcode to be considered a match.")
    public int MAX_MISMATCHES = 1;

    @Option(doc = "Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.")
    public int MIN_MISMATCH_DELTA = 2;

    @Option(doc = "Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.")
    public int MAX_NO_CALLS = 2;

    @Option(shortName = "Q", doc = "Minimum base quality. Any barcode bases falling below this quality will be considered a mismatch even in the bases match.")
    public int MINIMUM_BASE_QUALITY = 0;

    private static final SolexaQualityConverter solexaQualityConverter = SolexaQualityConverter.getSingleton();

    public static final String UNMATCHED_SAMPLE_ID = "unmatched";

    final Log log = Log.getInstance(DemuxFastqs.class);

    @Override
    protected int doWork() {
        // Check all input and output files are readable or writable respectively.
        READ_ONE_FASTQ.stream().forEach(IOUtil::assertFileIsReadable);
        READ_TWO_FASTQ.stream().forEach(IOUtil::assertFileIsReadable);
        I7_FASTQ.stream().forEach(IOUtil::assertFileIsReadable);
        I5_FASTQ.stream().forEach(IOUtil::assertFileIsReadable);
        IOUtil.assertFileIsReadable(SAMPLE_SHEET);
        IOUtil.assertDirectoryIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(METRICS_FILE);

        // Make some data structures to help us on our way
        final SampleSheet sampleSheet = new SampleSheet(SAMPLE_SHEET);
        final ReadStructureInfo readOneReadStructureInfo = new ReadStructureInfo(READ_ONE_READ_STRUCTURE);
        final ReadStructureInfo readTwoReadStructureInfo = new ReadStructureInfo(READ_TWO_READ_STRUCTURE);
        final ReadStructureInfo i7ReadStructureInfo = new ReadStructureInfo(I7_READ_STRUCTURE);
        final ReadStructureInfo i5ReadStructureInfo = new ReadStructureInfo(I5_READ_STRUCTURE);

        if (readOneReadStructureInfo.templateCycles.length == 0) {
            throw new PicardException("No template bases found in READ_ONE_READ_STRUCTURE.");
        }
        if (readTwoReadStructureInfo.templateCycles.length == 0) {
            throw new PicardException("No template bases found in READ_TWO_READ_STRUCTURE.");
        }
        if (i7ReadStructureInfo.templateCycles.length > 0) {
            throw new PicardException("Template bases not allowed in the i7 read.");
        }
        if (i5ReadStructureInfo.templateCycles.length > 0) {
            throw new PicardException("Template bases not allowed in the i5 read.");
        }

        // Set the barcode bases for each sample
        setSampleBarcode(sampleSheet,
                readOneReadStructureInfo,
                readTwoReadStructureInfo,
                i7ReadStructureInfo,
                i5ReadStructureInfo
        );

        // Determine the quality format
        determineQualityFormat();

        // Create the files for writing
        final SAMFileWriter[] writers = createSamFileWriters(sampleSheet);

        // Open the input FASTQs for reading
        final ParallelFastqReader reader = new ParallelFastqReader(READ_ONE_FASTQ, READ_TWO_FASTQ, I7_FASTQ, I5_FASTQ);

        // Read in the quadruples of FASTQ records
        final FastqConverter fastqConverter = new FastqConverter(sampleSheet,
                readOneReadStructureInfo,
                readTwoReadStructureInfo,
                i7ReadStructureInfo,
                i5ReadStructureInfo);
        final ProgressLogger progress = new ProgressLogger(log, 1000000);
        while (reader.hasNext()) {
            final List<FastqRecord> records = reader.next();
            // convert the FASTQ records to SAM records and output
            final List<SAMRecord> samRecords = fastqConverter.convertFastqsAndWriteRecords(records, writers);
            samRecords.forEach(progress::record);
        }

        // close them all
        reader.close();
        Arrays.stream(writers).forEach(SAMFileWriter::close);

        // write the metrics
        final MetricsFile<BarcodeMetric, Integer> metrics = getMetricsFile();
        fastqConverter.getBarcodeMetrics().forEach(metrics::addMetric);
        metrics.write(METRICS_FILE);

        return 0;
    }

    private void determineQualityFormat() {
        final FastqReader[] readers = new FastqReader[READ_ONE_FASTQ.size()];
        for (int i = 0; i < READ_ONE_FASTQ.size(); i++) {
            readers[i] = new FastqReader(READ_ONE_FASTQ.get(i), ALLOW_AND_IGNORE_EMPTY_LINES);
        }
        final QualityEncodingDetector detector = new QualityEncodingDetector();
        detector.add(QualityEncodingDetector.DEFAULT_MAX_RECORDS_TO_ITERATE, readers);
        for (final FastqReader reader : readers) reader.close();

        final FastqQualityFormat qualityFormat =  detector.generateBestGuess(QualityEncodingDetector.FileContext.FASTQ, QUALITY_FORMAT);
        if (detector.isDeterminationAmbiguous()) {
            log.warn("Making ambiguous determination about fastq's quality encoding; more than one format possible based on observed qualities.");
        }
        log.info(String.format("Auto-detected quality format as: %s.", qualityFormat));

        QUALITY_FORMAT = qualityFormat;
    }

    /** Creates SAM file writers for each sample in the Sample Sheet, returning them in the same order as specified
     * by each sample's ordinal.  An extra writer is appended for records that do not match a sample barcode. */
    private SAMFileWriter[] createSamFileWriters(final SampleSheet sampleSheet) {
        final SAMFileWriter[] writers = new SAMFileWriter[sampleSheet.size() + 1];
        for (final Sample sample : sampleSheet) {
            final File outputFile = new File(this.OUTPUT, IOUtil.makeFileNameSafe(String.format("%s.%s.%s.bam", sample.sampleId, sample.sampleName, sample.sampleBarcode)));
            final SAMFileHeader header = createSamFileHeader(sample);
            final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, outputFile);
            writers[sample.sampleOrdinal-1] = writer;
        }
        // for unmatched reads
        {
            final File outputFile = new File(this.OUTPUT, UNMATCHED_OUTPUT);
            final SAMFileHeader header = createSamFileHeader(new Sample(
                    writers.length-1,
                    UNMATCHED_SAMPLE_ID,
                    UNMATCHED_SAMPLE_ID,
                    UNMATCHED_SAMPLE_ID,
                    null,
                    null)
            );
            final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, outputFile);
            writers[writers.length-1] = writer;
        }

        return writers;
    }

    /** Creates a read group for the SAM header with the values provided on the command line.. */
    private SAMReadGroupRecord createSamReadGroupRecord(final String readGroupId,
                                                        final String sampleName,
                                                        final String libraryName,
                                                        final String platformUnit) {
        final SAMReadGroupRecord rgroup = new SAMReadGroupRecord(readGroupId);
        rgroup.setSample(sampleName);
        rgroup.setLibrary(libraryName);
        rgroup.setPlatformUnit(platformUnit);
        if (this.PLATFORM != null) rgroup.setPlatform(this.PLATFORM);
        if (this.SEQUENCING_CENTER != null) rgroup.setSequencingCenter(SEQUENCING_CENTER);
        if (this.PREDICTED_INSERT_SIZE != null) rgroup.setPredictedMedianInsertSize(PREDICTED_INSERT_SIZE);
        if (this.DESCRIPTION != null) rgroup.setDescription(this.DESCRIPTION);
        if (this.RUN_DATE != null) rgroup.setRunDate(this.RUN_DATE);
        if (this.PLATFORM_MODEL != null) rgroup.setPlatformModel(this.PLATFORM_MODEL);
        if (this.PROGRAM_GROUP != null) rgroup.setProgramGroup(this.PROGRAM_GROUP);
        return rgroup;
    }

    /** Creates a simple header with the values provided on the command line. */
    private SAMFileHeader createSamFileHeader(final Sample sample) {
        final SAMFileHeader header = new SAMFileHeader();

        // Update library id and platform unit if they have not been set.
        final String libraryId = (sample.libraryId == null) ? sample.sampleName : sample.libraryId;
        final String platformUnit = (sample.sampleBarcode == null) ? sample.sampleName : (sample.sampleName + "." + sample.sampleBarcode);
        final SAMReadGroupRecord rgroup = createSamReadGroupRecord(
                sample.sampleId,
                sample.sampleName,
                libraryId,
                platformUnit
        );

        header.addReadGroup(rgroup);

        COMMENT.forEach(header::addComment);

        header.setSortOrder(this.SORT_ORDER);
        return header;
    }

    /** Gets all the bases (concatenated) for the given cycle indexes (1-based), ignoring cycles not after the end
     * of `bases`.  The cycles should be in ascending order. */
    private static byte[] getBasesAtCycles(final byte[] bases, final int[] cycles) {
        // first get the number of bases we should collate, since we ignore cycles past the end of `bases`.
        int numBases = 0;
        for (final int cycle : cycles) {
            if (cycle-1 >= bases.length) break;
            numBases++;
        }
        // now fill in the bases
        final byte[] basesAtCycles = new byte[numBases];
        for (int i=0; i<numBases; ++i) {
            basesAtCycles[i] = bases[cycles[i]-1]; // NB: cycles[i] is one-based
        }
        return basesAtCycles;
    }

    /** Sets the sample barcode in the given sample sheet using the given read structures and bases from the sample sheet.
     *
     * The final sample barcode is the concatenation of the sample barcode from the i7 and i5 bases found in the sample
     * sheet in that order.  The concatenation uses `IlluminaUtil.barcodeSeqsToString`.  The specific bases extracted
     * from the i7 and i5 reads are the sample barcode bases specified in the given i7 and i5 read structures.  The
     * final sample barcode is then updated in the given sample sheet.  An exception is thrown if no sample barcode
     * bases are found in the i5 or i7 read structures.
     * */
    private void setSampleBarcode(final SampleSheet sampleSheet,
                                  final ReadStructure r1ReadStructure,
                                  final ReadStructure r2ReadStructure,
                                  final ReadStructure i7ReadStructure,
                                  final ReadStructure i5ReadStructure) {
        final ReadStructure [] readStructures = new ReadStructure[] { r1ReadStructure, r2ReadStructure, i7ReadStructure, i5ReadStructure};
        for (final Sample sample : sampleSheet) {
            final List<String> sampleBarcodes = new ArrayList<>();
            final byte[][] barcodeBases = new byte[][]{ sample.r1SampleBarcodeBases, sample.r2SampleBarcodeBases, sample.i7IndexBases, sample.i5IndexBases };
            for (int i = 0; i < barcodeBases.length; i++) {
                addToSampleBarcodeList(sampleBarcodes, readStructures[i], barcodeBases[i]);
            }
            if (sampleBarcodes.isEmpty()) {
                throw new PicardException("No sample barcodes found in the i7 or i5 read structures.");
            }
            sample.sampleBarcode = new SampleBarcode(sampleBarcodes);
        }
    }

    private void addToSampleBarcodeList(final List<String> sampleBarcodes,
                                        final ReadStructure readStructure,
                                        final byte[] barcodeBases) {
        if (0 < readStructure.sampleBarcodes.getTotalCycles()) {
            sampleBarcodes.add(new String(getBasesAtCycles(barcodeBases, readStructure.sampleBarcodes.getCycles())));
        }
    }

    /** Counts the nucleotide mismatches between two strings of the same length.  Ignores no calls in expectedBases.
     * Observed base qualities less than the minimum base quality are counted as mismatches if not a  no call. Observed
     * qualities may be null. */
    static int countMismatches(final byte[] observedBases,
                               byte[] observedQualities,
                               final byte[] expectedBases,
                               final int minimumBaseQuality) {
        int mismatches = 0;
        for (int baseNumber = 0; baseNumber < observedBases.length; baseNumber++) {
            if (!SequenceUtil.isNoCall(expectedBases[baseNumber])) {
                if (!SequenceUtil.basesEqual(observedBases[baseNumber], expectedBases[baseNumber])) ++mismatches;
                else if (observedQualities[baseNumber] < minimumBaseQuality) ++mismatches;
            }
        }
        return mismatches;
    }

    class ReadStructureInfo extends ReadStructure {
        private final int[] templateCycles;
        private final int[] sampleBarcodeCycles;
        private final int[] molecularBarcodeCycles;

        public ReadStructureInfo(final String readStructureString) {
            super(readStructureString);

            this.templateCycles = this.templates.getCycles();
            this.sampleBarcodeCycles = this.sampleBarcodes.getCycles();
            this.molecularBarcodeCycles = this.molecularBarcode.getCycles();
        }
    }

    /** Helper class to convert fastq records into SAM records and write them out */
    private class FastqConverter {
        private final List<ReadStructureInfo> readStructureInfos = new ArrayList<>();
        private final SampleSheet sampleSheet;

        /** Sample barcode to its associate metric */
        private final Map<String, BarcodeMetric> barcodeToMetrics = new LinkedHashMap<>();
        private final BarcodeMetric noMatchBarcodeMetric;

        /**
         *
         * @param sampleSheet the sample sheet describing all the samples
         * @param readStructureInfos the read structures info, one for each expected fastq record that will be passed
         *                           into `convertFastqsAndWriteRecords`.
         */
        public FastqConverter(final SampleSheet sampleSheet, final ReadStructureInfo... readStructureInfos) {
            this.sampleSheet = sampleSheet;
            this.readStructureInfos.addAll(Arrays.asList(readStructureInfos));

            // Initialize the barcode metrics
            for (final Sample sample : sampleSheet) {
                final String barcode = sample.sampleBarcode.concatenatedBarcode;
                final BarcodeMetric metric = new BarcodeMetric(sample.sampleName, sample.libraryId, barcode, new String[]{barcode});
                barcodeToMetrics.put(barcode, metric);
            }
            final List<String> noMatchBarcodes = new ArrayList<>();
            for (final ReadStructureInfo readStructureInfo : readStructureInfos) {
                final int totalCycles = readStructureInfo.sampleBarcodes.getTotalCycles();
                if (0 < totalCycles) noMatchBarcodes.add(StringUtil.repeatCharNTimes('N', totalCycles));
            }
            final String noMatchBarcode = IlluminaUtil.barcodeSeqsToString(noMatchBarcodes);
            noMatchBarcodeMetric = new BarcodeMetric(UNMATCHED_SAMPLE_ID, UNMATCHED_SAMPLE_ID, noMatchBarcode, new String[]{noMatchBarcode});
            barcodeToMetrics.put(noMatchBarcode, noMatchBarcodeMetric);
        }

        /** Converts the given FASTQ records to SAM records and writes them out to the output file based on the
         * sample barcode. */
        public List<SAMRecord> convertFastqsAndWriteRecords(final List<FastqRecord> fastqRecords,
                                                            final SAMFileWriter[] writers) {
            final byte[] sampleBarcodeReadBases = getSampleBarcode(fastqRecords);
            final byte[] sampleBarcodeReadQualities = getSampleBarcodeQualities(fastqRecords);

            // match the observed sample barcode to those in the sample sheet
            int sampleOrdinal = getSampleOrdinalFromSampleBarcode(sampleBarcodeReadBases, sampleBarcodeReadQualities, sampleSheet);
            if (sampleOrdinal == -1) sampleOrdinal = writers.length;

            // get the writer and header
            final SAMFileWriter writer = writers[sampleOrdinal-1];
            final SAMFileHeader header = writer.getFileHeader();

            // Get the read group identifier to use for the records; we assume there is only one read group in the header.
            final String sampleId = writers[sampleOrdinal-1].getFileHeader().getReadGroups().get(0).getId();

            // Create the molecular barcode sequence tag
            final String molecularBarcodeTag = getMolecularBarcode(fastqRecords);

            // Create the SAM records
            final List<SAMRecord> samRecords = new ArrayList<>();
            for (int i = 0; i < fastqRecords.size(); i++) {
                final FastqRecord record = fastqRecords.get(i);
                final ReadStructureInfo readStructureInfo = readStructureInfos.get(i);

                // ignore records without template sequence
                if (readStructureInfo.templates.isEmpty()) continue;

                // now make sure we do not have more than two ends
                if (samRecords.size() == 2) throw new PicardException("Found more than two ends with template bases");

                final SAMRecord samRecord = makeSAMRecord(header, record, readStructureInfo, sampleId, samRecords.isEmpty());
                samRecord.setAttribute("BC", new String(sampleBarcodeReadBases));
                samRecord.setAttribute(MOLECULAR_BARCODE_TAG, molecularBarcodeTag);
                writer.addAlignment(samRecord);
                samRecords.add(samRecord);
            }
            return samRecords;
        }

        public Collection<BarcodeMetric> getBarcodeMetrics() {
            ExtractIlluminaBarcodes.finalizeMetrics(barcodeToMetrics, noMatchBarcodeMetric);
            return this.barcodeToMetrics.values();
        }

        /** Searches for a matching sample given the observed barcode.  Returns -1 if no match is found given the
         * sample sheet and matching parameters */
        private int getSampleOrdinalFromSampleBarcode(final byte[] sampleBarcodeReadBases,
                                                      final byte[] sampleBarcodeReadQualities,
                                                      final SampleSheet sampleSheet) {


            // Count the # of no calls
            int numNoCalls = 0;
            for (final byte base: sampleBarcodeReadBases) {
                if (SequenceUtil.isNoCall(base)) numNoCalls++;
            }

            int bestSampleOrdinal = -1;
            int bestMismatches = Integer.MAX_VALUE;
            int secondBestMismatches = Integer.MAX_VALUE;

            if (numNoCalls <= MAX_NO_CALLS) { // do not do any matching if we have too many no-calls
                for (final Sample sample : sampleSheet) {
                    final int numMismatches = countMismatches(sampleBarcodeReadBases,
                            sampleBarcodeReadQualities,
                            sample.sampleBarcode.barcodeBytes,
                            MINIMUM_BASE_QUALITY);
                    if (numMismatches < bestMismatches) {
                        bestSampleOrdinal = sample.sampleOrdinal;
                        secondBestMismatches = bestMismatches;
                        bestMismatches = numMismatches;
                    } else if (numMismatches < secondBestMismatches) {
                        secondBestMismatches = numMismatches;
                    }
                }
            }

            // Make sure we are within the parameter limits and update barcode metrics if necessary
            if (MAX_MISMATCHES < bestMismatches ||
                    MAX_NO_CALLS < numNoCalls ||
                    (secondBestMismatches - bestMismatches) < MIN_MISMATCH_DELTA) {
                bestSampleOrdinal = -1;
                ++noMatchBarcodeMetric.READS;
                ++noMatchBarcodeMetric.PF_READS;
            }
            else {
                final BarcodeMetric bestBarcodeMetric = this.barcodeToMetrics.get(sampleSheet.get(bestSampleOrdinal-1).sampleBarcode.concatenatedBarcode);
                ++bestBarcodeMetric.READS;
                ++bestBarcodeMetric.PF_READS;
                if (bestMismatches == 0) {
                    ++bestBarcodeMetric.PERFECT_MATCHES;
                    ++bestBarcodeMetric.PF_PERFECT_MATCHES;
                } else if (bestMismatches == 1) {
                    ++bestBarcodeMetric.ONE_MISMATCH_MATCHES;
                    ++bestBarcodeMetric.PF_ONE_MISMATCH_MATCHES;
                }
            }

            return bestSampleOrdinal;
        }

        /** Creates a SAM record (always paired) */
        private SAMRecord makeSAMRecord(final SAMFileHeader header,
                                        final FastqRecord fastqRecord,
                                        final ReadStructureInfo readStructureInfo,
                                        final String readGroupId,
                                        final boolean firstOfPair) {
            // get read bases and quality values
            final byte[] readBases = getBasesAtCycles(fastqRecord.getReadString().getBytes(), readStructureInfo.templateCycles);
            final byte[] baseQualities = getBasesAtCycles(fastqRecord.getBaseQualityString().getBytes(), readStructureInfo.templateCycles);

            final SAMRecord record = new SAMRecord(header);
            record.setReadName(getReadName(fastqRecord.getReadHeader(), true, STRIP_UNPAIRED_MATE_NUMBER));
            record.setReadBases(readBases);
            record.setReadUnmappedFlag(true);
            record.setAttribute(ReservedTagConstants.READ_GROUP_ID, readGroupId);

            // Set base qualities
            convertQuality(baseQualities, QUALITY_FORMAT);
            for (final byte qual : baseQualities)  {
                final int uQual = qual & 0xff;
                if (uQual < MIN_Q || uQual > MAX_Q) {
                    throw new PicardException("Base quality " + uQual + " is not in the range " + MIN_Q + ".." +
                            MAX_Q + " for read " + fastqRecord.getReadHeader());
                }
            }
            record.setBaseQualities(baseQualities);

            // always paired
            record.setReadPairedFlag(true);
            record.setMateUnmappedFlag(true);

            if (firstOfPair) record.setFirstOfPairFlag(true);
            else record.setSecondOfPairFlag(true);

            return record;
        }

        /** Based on the type of quality scores coming in, converts them to a numeric byte[] in phred scale. */
        void convertQuality(final byte[] quals, final FastqQualityFormat version) {
            switch (version)  {
                case Standard:
                    SAMUtils.fastqToPhred(quals);
                    break ;
                case Solexa:
                    solexaQualityConverter.convertSolexaQualityCharsToPhredBinary(quals);
                    break ;
                case Illumina:
                    solexaQualityConverter.convertSolexa_1_3_QualityCharsToPhredBinary(quals);
                    break ;
            }
        }

        /** Gets an array of byte arrays representating the sample barcode found in each record, according to their
         * individual read structures.  This will include empty byte[]s for any reads that don't contain a sample
         * barcode read.
         */
        private byte[] getSampleBarcode(final List<FastqRecord> records) {
            final int numCycles = readStructureInfos.stream().mapToInt(r -> r.sampleBarcodes.getTotalCycles()).sum();
            final byte[] sampleBarcode = new byte[numCycles];
            int sampleBarcodeOffset = 0;
            for (int i = 0; i < records.size(); i++) {
                final FastqRecord record = records.get(i);
                final ReadStructureInfo readStructureInfo = readStructureInfos.get(i);
                final byte[] readBarcode = getBasesAtCycles(record.getReadString().getBytes(), readStructureInfo.sampleBarcodeCycles);
                System.arraycopy(readBarcode, 0, sampleBarcode, sampleBarcodeOffset, readBarcode.length);
                sampleBarcodeOffset += readBarcode.length;
            }
            return sampleBarcode;
        }

        private byte[] getSampleBarcodeQualities(final List<FastqRecord>  records) {
            final int numCycles = readStructureInfos.stream().mapToInt(r -> r.sampleBarcodes.getTotalCycles()).sum();
            final byte[] sampleBarcodeQualities = new byte[numCycles];
            int sampleBarcodeQualitiesOffset = 0;
            for (int i = 0; i < records.size(); i++) {
                final FastqRecord record = records.get(i);
                final ReadStructureInfo readStructureInfo = readStructureInfos.get(i);
                final byte[] readBarcodeQualities = getBasesAtCycles(record.getBaseQualityString().getBytes(), readStructureInfo.sampleBarcodeCycles);
                System.arraycopy(readBarcodeQualities, 0, sampleBarcodeQualities, sampleBarcodeQualitiesOffset, readBarcodeQualities.length);
                sampleBarcodeQualitiesOffset += readBarcodeQualities.length;
            }
            return sampleBarcodeQualities;
        }

        private String getMolecularBarcode(final List<FastqRecord>  records) {
            final byte[][] molecularBarcode = new byte[records.size()][];
            for (int i = 0; i < records.size(); i++) {
                final FastqRecord record = records.get(i);
                final ReadStructureInfo readStructureInfo = readStructureInfos.get(i);
                molecularBarcode[i] = getBasesAtCycles(record.getReadString().getBytes(), readStructureInfo.molecularBarcodeCycles);
            }
            return IlluminaUtil.barcodeSeqsToString(molecularBarcode);
        }
    }

    /** Reads multiple FASTQs at the same time, ensuring their read names match along the way.  This is
     * useful for when we don't have a single interleaved FASTQ. */
    class ParallelFastqReader implements Iterator<List<FastqRecord>> {
        public final List<SequentialFastqReader> readers = new ArrayList<>();

        @SafeVarargs
        public ParallelFastqReader(final List<File> fastq, final List<File>... fastqs) {
            readers.add(new SequentialFastqReader(fastq));
            for (final List<File> otherFastqs : fastqs) {
                if (otherFastqs.size() != fastq.size()) throw new PicardException("List of list fastqs had differeing lengths.");
                readers.add(new SequentialFastqReader(otherFastqs));
            }
        }

        public List<FastqRecord> next() {
            final List<FastqRecord> records = new ArrayList<>();
            records.addAll(readers.stream().map(SequentialFastqReader::next).collect(Collectors.toList()));
            verifySameReadNames(records);
            return records;
        }

        public boolean hasNext() {
            if (readersHaveNext(true)) return true;
            else if (readersHaveNext(false)) return false;
            else throw new PicardException("Fastqs have differing lengths"); // TODO: better error message
        }

        private boolean readersHaveNext(boolean desiredValue) {
            for (final SequentialFastqReader reader : readers) {
                if (reader.hasNext() != desiredValue) return false;
            }
            return true;
        }

        private void verifySameReadNames(final List<FastqRecord> records) {
            final String readName = getReadName(records.get(0).getReadHeader(), true, STRIP_UNPAIRED_MATE_NUMBER);
            for (int i = 1; i < records.size(); i++) {
                final String otherReadName = getReadName(records.get(i).getReadHeader(), true, STRIP_UNPAIRED_MATE_NUMBER);
                if (readName.compareTo(otherReadName) != 0) {
                    throw new PicardException(String.format("Mismatching read names in FASTQS:\n%s\n%s\n", readName, otherReadName));
                }
            }
        }

        public void close() {
            readers.forEach(SequentialFastqReader::close);
        }
    }

    // TODO: duplicated but private in picard.sam.FastqToSam so make it public and static
    // Read names cannot contain blanks
    private static String getReadName(final String fastqHeader, final boolean paired, final boolean stripUnpairedMateNumber) {
        final int idx = fastqHeader.indexOf(" ");
        String readName = (idx == -1) ? fastqHeader : fastqHeader.substring(0,idx);

        // NOTE: the while loop isn't necessarily the most efficient way to handle this but we don't
        // expect this to ever happen more than once, just trapping pathological cases
        while (stripUnpairedMateNumber && !paired && (readName.endsWith("/1") || readName.endsWith("/2"))) {
            // If this is an unpaired run we want to make sure that "/1" isn't tacked on the end of the read name,
            // as this can cause problems down the road in MergeBamAlignment
            readName = readName.substring(0, readName.length() - 2);
        }

        return readName;
    }
}