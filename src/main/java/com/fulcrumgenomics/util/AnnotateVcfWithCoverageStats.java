/**
 * Copyright (c) 2015, Fulcrum Genomics LLC
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
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFStandardHeaderLines;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Simple tool to annotate a VCF with read counts per read group.
 */
@CommandLineProgramProperties(
        usage = "Reads a VCF/VCF.gz/BCF  and a BAM file (with read groups) and creates a new VCF with per read group (each their own sample) coverage statistics (SNPs only).",
        usageShort = "Annotate a VCF with per read group coverage statistics.",
        programGroup = Utilities.class
)
public class AnnotateVcfWithCoverageStats extends CommandLineProgram {
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input SAM or BAM")
    public File INPUT;

    @Option(shortName= "V", doc="Input VCF with no genyptes")
    public File DBSNP_VCF;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF or BCF to emit with per-sample info.")
    public File OUTPUT;

    @Option(doc="The SAM or BAM output to write overlapping reads (may split across read pairs).", optional = true)
    public File OUTPUT_BAM = null;

    @Option(doc="Output only overlapping reads that match the alternative allele.", optional = true)
    public boolean OUTPUT_RECORDS_ONLY_WITH_ALT_BASES = true;

    @Option(doc="The minimum mapping quality to accept", optional = true)
    public Integer MIN_MAPQ = 20;

    @Option(shortName="MAX_REF_AF", doc="The maximum reference allele frequency as reported in the CAF INFO field.", optional = true)
    public Double MAX_REFERENCE_ALLELE_FREQUENCY = 0.95;

    @Option(doc="The interesting intervals file to which to output", optional = true)
    public File OUTPUT_INTERVALS = null;

    @Option(doc="The maximum distance...", optional = true)
    public int MAX_SAMPLE_DISTANCE_FOR_INTERVALS = 10000;

    @Option(doc="The minimum distance...", optional = true)
    public int MIN_INTERVAL_LENGTH = 1000;

    final Log log = Log.getInstance(AnnotateVcfWithCoverageStats.class);

    // Stock main method
    public static void main(final String[] args) {
        new AnnotateVcfWithCoverageStats().instanceMainWithExit(args);
    }

    public AnnotateVcfWithCoverageStats() {
        CREATE_INDEX = true;
    }

    private class InterestingIntervals {

        private final int maximumDistance;
        private int currentSequenceIndex = -1;
        private final Map<String, Integer> sampleNameToPositionMap = new HashMap<>();
        private final SAMSequenceDictionary dictionary;
        private final List<String> sampleNames;
        private IntervalList intervals;
        private int MAX_INTERVALS_UNTIL_UNIQUING = 10000;
        private final int minimumIntervalLength;

        public InterestingIntervals(final int maximumDistance, final int minimumIntervalLength, final List<String> sampleNames, final SAMFileHeader header) {
            this.maximumDistance = maximumDistance;
            this.minimumIntervalLength = minimumIntervalLength;
            this.dictionary = header.getSequenceDictionary();
            this.sampleNames = sampleNames;
            this.intervals = new IntervalList(header);
            reset(this.currentSequenceIndex);
        }

        void reset(final int sequenceIndex) {
            this.sampleNameToPositionMap.clear();
            for (final String sampleName : sampleNames) {
                sampleNameToPositionMap.put(sampleName, -1);
            }
            this.currentSequenceIndex = sequenceIndex;
        }

        void add(final VariantContext variantContext) {

            final int contextSequenceIndex = this.dictionary.getSequenceIndex(variantContext.getContig());
            final int contextPosition = variantContext.getStart();

            // if we move to a new contig
            if (this.currentSequenceIndex != contextSequenceIndex) {
                this.reset(contextSequenceIndex);
            }

            // Go through each sample's genotype
            for (final Genotype genotype : variantContext.getGenotypes()) {
                if (!genotype.isCalled()) continue;
                final String sampleName = genotype.getSampleName();
                final int previousPosition = this.sampleNameToPositionMap.get(sampleName);
                this.sampleNameToPositionMap.put(sampleName, contextPosition);

                if (-1 == previousPosition) continue; // not interesting to go from the start to here
                if (contextPosition == previousPosition) continue; // ignore

                final int distance = contextPosition - previousPosition;
                if (distance <= maximumDistance) {
                    final Interval interval = new Interval(variantContext.getContig(), previousPosition, contextPosition, false, sampleName + ":" + previousPosition + "-" + contextPosition);
                    intervals.add(interval);
                }
            }

            // Update intervals
            if (intervals.size() > MAX_INTERVALS_UNTIL_UNIQUING) {
                intervals = intervals.uniqued(true);
                MAX_INTERVALS_UNTIL_UNIQUING = 2 * MAX_INTERVALS_UNTIL_UNIQUING;
            }
        }

        // Should only be accessed when all intervals have been crated
        public IntervalList getIntervals() {
            final IntervalList candidates = intervals.uniqued();
            intervals = new IntervalList(intervals.getHeader());
            for (final Interval interval : candidates) {
                if (minimumIntervalLength <= interval.length()) {
                    intervals.add(interval);
                }
            }
            return intervals;
        }
    }

    private List<String> getReadGroupIdListFromHeader(final SAMFileHeader header) {
        final List<String> samples = new ArrayList<>();
        for (final SAMReadGroupRecord record : header.getReadGroups()) {
            samples.add(record.getReadGroupId());
        }
        return samples;
    }

    private static final byte[] BASES = {'A', 'C', 'G', 'T', 'N'};

    private class AlleleCounts {
        private final Map<String, int[]> alleleCounts = new HashMap<>();
        private final List<String> samples;

        public AlleleCounts(final List<String> samples) {
            this.samples = samples;
            for (final String sample : samples) {
                final int[] counts = new int[BASES.length];
                alleleCounts.put(sample, counts);
            }
            clear();
        }

        public void clear() {
            for (final String sample : samples) {
                final int[] counts = alleleCounts.get(sample);
                for (int i = 0; i < BASES.length; i++) {
                    counts[i] = 0;
                }
            }
        }

        private void addBase(final String sample, final byte base) {
            final int[] counts = alleleCounts.get(sample);
            int i;
            for (i = 0; i < BASES.length; i++) {
                if (base == BASES[i]) {
                    counts[i]++;
                    break;
                }
            }
            if (i == BASES.length) {
                throw new PicardException("Could not find count for base: " + (char)base);
            }
        }

        public void addN(final String sample) {
            addBase(sample, (byte)'N');
        }

        public void add(final String sample, final Allele allele) {
            final byte[] bases = allele.getBases();

            if (1 != bases.length) {
                throw new PicardException("Allele did not have exactly one base: " + allele.toString());
            }

            final byte base = bases[0];
            addBase(sample, base);
        }

        public int get(final String sample, final Allele allele) {
            final byte[] bases = allele.getBases();

            if (1 != bases.length) {
                throw new PicardException("Allele did not have exactly one base: " + allele.toString());
            }

            final byte base = bases[0];

            final int[] counts = alleleCounts.get(sample);
            int i;
            for (i = 0; i < BASES.length; i++) {
                if (base == BASES[i]) {
                    return counts[i];
                }
            }
            throw new PicardException("Could not find count for base: " + allele.getBaseString());
        }

        public int getSumOfCounts(final String sample) {
            int sum = 0;
            final int[] counts = alleleCounts.get(sample);
            for (int i = 0; i < BASES.length; i++) {
                sum += counts[i];
            }
            return sum;
        }
    }

    /**
     * @return true if it matches one the allelles (or only non-reference alleles if OUTPUT_RECORDS_ONLY_WITH_ALT_BASES is true), false otherwise
     */
    private boolean updateAlleleCounts(final SAMRecord record, final List<Allele> alleles, final int positionToMatch, final String sampleName, final AlleleCounts alleleCounts) {

        boolean foundAltBase = false;
        boolean foundMatch = false;

        try {
            // Get the read base matching the position
            final byte[] readBases = record.getReadBases();
            byte readBase = -1;
            for (final AlignmentBlock block : record.getAlignmentBlocks()) {
                final int readBlockStart = block.getReadStart() - 1; // 0-based
                final int referenceBlockStart = block.getReferenceStart(); // 1-based
                final int length = block.getLength();

                if (referenceBlockStart <= positionToMatch && positionToMatch < referenceBlockStart + length) {
                    final int offset = positionToMatch - referenceBlockStart;
                    readBase = readBases[readBlockStart + offset];
                    break;
                }
            }
            // Not considered out since we could have deletions
            /*
            if (-1 == readBase) {
                throw new PicardException("Could not find the read base at the given position: " + positionToMatch + " (" + record.getAlignmentStart() + "-" + record.getAlignmentEnd() + ")");
            }
            */

            // To which allele does it match?
            for (final Allele allele : alleles) {
                final byte[] alleleBases = allele.getBases();
                if (1 != alleleBases.length) {
                    throw new PicardException("Allele did not have exactly one base: " + allele.toString());
                }
                if (readBase == alleleBases[0]) {
                    alleleCounts.add(sampleName, allele);
                    foundMatch = true;
                    if (allele.isNonReference()) {
                        foundAltBase = true;
                    }
                    break;
                }
            }
            // N
            if (!foundMatch) {
                alleleCounts.addN(sampleName);
            }
        } catch (final Exception e) {
            throw new SAMException("Exception for read " + record, e);
        }

        if (OUTPUT_RECORDS_ONLY_WITH_ALT_BASES) {
            return foundAltBase;
        }
        else {
            return foundMatch;
        }
    }

    private VariantContext addSamplesToVariantContext(final VariantContext variantContext, final SamReader samReader,
                                                      final SAMFileWriter samFileWriter, final Set<Integer> recordsSeen,
                                                      final AlleleCounts alleleCounts,
                                                      final List<String> sampleNames) {
        final List<Allele> alleles = variantContext.getAlleles();

        // clear any previous counts
        alleleCounts.clear();

        // Get an iterator over overlapping reads
        final SAMRecordIterator iterator = samReader.queryOverlapping(variantContext.getContig(), variantContext.getStart(), variantContext.getEnd());
        final Set<Integer> recordsSeenThisTime = new HashSet<>(); // for the next time this function is called
        // Go through each overlapping record
        for (final SAMRecord record : new IterableAdapter<>(iterator)) {

            // Exclude secondary, supplementary, duplicate, low mapq, etc.
            if (record.isSecondaryOrSupplementary() || record.getDuplicateReadFlag() || record.getMappingQuality() < MIN_MAPQ) continue;

            // Assumes that a read group exists, very naughty
            final SAMReadGroupRecord readGroup = record.getReadGroup();
            if (null == readGroup) {
                throw new PicardException("No read group found for record: " + record.getSAMString());
            }

            // Update the allele counts
            final int code = record.hashCode();
            if (updateAlleleCounts(record, alleles, variantContext.getStart(), readGroup.getReadGroupId(), alleleCounts)) {
                if (null != samFileWriter) {
                    //System.err.print("(" + variantContext.getStart() + "-" + variantContext.getEnd() + ":" + variantContext.getID() + "): " + record.getSAMString());
                    if (!recordsSeen.contains(code)) { // do not add records seen in previous iterations
                        samFileWriter.addAlignment(record);
                    }
                    recordsSeenThisTime.add(code);
                }
            }
            else if (recordsSeen.contains(code)) { // since we may be looking at a different set of alleles
                recordsSeenThisTime.add(code);
            }
        }
        iterator.close();

        // Update the records seen
        recordsSeen.clear();
        recordsSeen.addAll(recordsSeenThisTime);

        // Update the genotype fields
        // DP - read depth for the given sample
        // AC - allele counts for REF, ALTs
        final List<Genotype> genotypes = new ArrayList<>();
        for (final String sampleName : sampleNames) {
            final int dp = alleleCounts.getSumOfCounts(sampleName);
            final List<Integer> counts = new ArrayList<>();
            final List<Allele> genotypeAlleles = new ArrayList<>();
            for (final Allele allele : alleles) {
                final int count = alleleCounts.get(sampleName, allele);
                counts.add(count);
                if (0 < count) {
                    genotypeAlleles.add(allele);
                }
            }
            final Genotype genotype;
            if (genotypeAlleles.isEmpty()) {
                genotype = GenotypeBuilder.createMissing(sampleName, 2);
            }
            else {
                genotype = new GenotypeBuilder(sampleName, genotypeAlleles).DP(dp).attribute("AC", counts).make();
            }
            genotypes.add(genotype);
        }

        // Make the variant context
        final VariantContextBuilder builder = new VariantContextBuilder(variantContext);
        builder.genotypes(genotypes);
        return builder.make();
    }

    @SuppressWarnings("unchecked")
    private boolean referenceAlleleFrequencyIsTooDamnHigh(final VariantContext context) {
        final Object attribute = context.getAttribute("CAF");
        if (attribute == null) {
            throw new PicardException("CAF missing form context: " + context.toString());
        }
        final List<String> alleleFrequencies = (List<String>)attribute;
        final double alleleFrequency = Double.parseDouble(alleleFrequencies.get(0));
        return (alleleFrequency > MAX_REFERENCE_ALLELE_FREQUENCY);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(DBSNP_VCF);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (null != OUTPUT_BAM) IOUtil.assertFileIsWritable(OUTPUT_BAM);
        if (null != OUTPUT_INTERVALS) IOUtil.assertFileIsWritable(OUTPUT_INTERVALS);

        final VCFFileReader vcfReader = new VCFFileReader(DBSNP_VCF, false);
        final VCFHeader inputVcfHeader = new VCFHeader(vcfReader.getFileHeader().getMetaDataInInputOrder());
        SAMSequenceDictionary sequenceDictionary = inputVcfHeader.getSequenceDictionary();

        final SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        final SAMFileWriter samFileWriter;
        if (null == OUTPUT_BAM) {
            samFileWriter = null;
        }
        else {
            // TODO: enforce that the output gets sorted
            samFileWriter = new SAMFileWriterFactory().setCreateIndex(CREATE_INDEX).makeSAMOrBAMWriter(samReader.getFileHeader(), false, OUTPUT_BAM);
        }

        if (null == sequenceDictionary) {
            log.warn("Input VCF sequence dictionary is missing, trying with SAM/BAM sequence dictionary");
            sequenceDictionary = samReader.getFileHeader().getSequenceDictionary();
        }

        if (CREATE_INDEX && sequenceDictionary == null) {
            throw new PicardException("A sequence dictionary must be available (either through the input file or by setting it explicitly) when creating indexed output.");
        }

        // Make sure the sequence dictionaries are the same
        sequenceDictionary.assertSameDictionary(samReader.getFileHeader().getSequenceDictionary());

        final ProgressLogger progress = new ProgressLogger(log, 10000);

        // Setup the VCF file writer
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
                .setOutputFile(OUTPUT)
                .setReferenceDictionary(sequenceDictionary);
        if (CREATE_INDEX) {
            builder.setOption(Options.INDEX_ON_THE_FLY);
        }
        else {
            builder.unsetOption(Options.INDEX_ON_THE_FLY);
        }
        final VariantContextWriter writer = builder.build();

        final List<String> sampleNames = getReadGroupIdListFromHeader(samReader.getFileHeader());

        // Make sure to add the samples to the new VCF, as well as the two extra FORMAT fields
        final Set<VCFHeaderLine> headerLines = new HashSet<>(inputVcfHeader.getMetaDataInInputOrder());
        VCFStandardHeaderLines.addStandardFormatLines(headerLines, false, Genotype.PRIMARY_KEYS); // This is very frustrating to have to add
        final VCFHeader header = new VCFHeader(headerLines, sampleNames);
        header.addMetaDataLine(new VCFFormatHeaderLine("DP", 1, VCFHeaderLineType.Integer, "The read depth for the given genotype"));
        header.addMetaDataLine(new VCFFormatHeaderLine("AC", VCFHeaderLineCount.A, VCFHeaderLineType.Integer, "The read depth per allele for the genotype"));
        writer.writeHeader(header);

        final AlleleCounts alleleCounts = new AlleleCounts(header.getSampleNamesInOrder());

        final InterestingIntervals interestingIntervals;
        if (null == OUTPUT_INTERVALS) {
            interestingIntervals = null;
        }
        else {
            interestingIntervals = new InterestingIntervals(MAX_SAMPLE_DISTANCE_FOR_INTERVALS, MIN_INTERVAL_LENGTH, sampleNames, samReader.getFileHeader());
        }

        // Go through the input
        final CloseableIterator<VariantContext> dbsnpIterator = vcfReader.iterator();
        final Set<Integer> recordsSeen = new HashSet<>(); // only really needed if we write to the SAM or BAM output
        while (dbsnpIterator.hasNext()) {
            final VariantContext full = dbsnpIterator.next();

            if (full.getNSamples() > 0) {
                throw new PicardException("Expecting no samples in the DBSNP VCF, but found " + full.getNSamples() + ".");
            }

            // SNPs only, for now
            if (!full.isSNP()) {
                continue;
            }

            // The reference allele frequency is too damn high
            if (referenceAlleleFrequencyIsTooDamnHigh(full)) {
                continue;
            }

            progress.record(full.getContig(), full.getStart());

            // TODO:
            // output an interval list of regions with potentially interesting linked reads

            // Add the new context
            final VariantContext context = addSamplesToVariantContext(full, samReader, samFileWriter, recordsSeen, alleleCounts, sampleNames);
            writer.add(context);
            if (null != interestingIntervals) interestingIntervals.add(context);
        }

        if (null != OUTPUT_INTERVALS && null != interestingIntervals) interestingIntervals.getIntervals().write(OUTPUT_INTERVALS);

        CloserUtil.close(dbsnpIterator);
        CloserUtil.close(vcfReader);
        CloserUtil.close(samReader);
        if (null != samFileWriter) CloserUtil.close(samFileWriter);

        writer.close();

        return 0;
    }
}
