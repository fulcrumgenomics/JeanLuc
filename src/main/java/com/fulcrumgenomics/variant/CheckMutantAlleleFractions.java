package com.fulcrumgenomics.variant;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.filter.DuplicateReadFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.util.*;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.Metrics;

import java.io.File;
import java.io.PrintWriter;
import java.util.*;

/**
 * Tool to check the mutant allele frequencies at a set of SNP sites described by a VCF. Each site
 * is downsampled (once) to the target coverage if it is above the specified coverage, and then the
 * reads are examined and counts of the ref and alt alleles produced.
 *
 * Output is a simple tab-separated text table with one row per SNP.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage = "Computes the mutant allele fraction at SNP loci after downsampling.",
        usageShort = "Computes the mutant allele fraction at SNP loci after downsampling.",
        programGroup = Metrics.class
)
public class CheckMutantAlleleFractions extends CommandLineProgram {
    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="The input SAM or BAM file.")
    public File INPUT;

    @Option(shortName="S", doc="A VCF containing the set of sites to extract mutant allele fractions at.")
    public File SITES_VCF;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output file to write to.")
    public File OUTPUT;

    @Option(shortName="DS", doc="The target coverage to downsample to.")
    public int DOWNSAMPLING_TARGET = 500;

    @Option(shortName=StandardOptionDefinitions.MINIMUM_MAPPING_QUALITY_SHORT_NAME,
            doc="Exclude reads with mapping quality below this value.")
    public int MINIMUM_MAPPING_QUALITY = 10;

    @Option(shortName="Q", doc="Exclude bases below this quality value.")
    public int MINIMUM_BASE_QUALITY = 10;

    @Option(doc="If true allow bases to be counted from overlapping reds from the same insert..")
    public boolean ALLOW_OVERLAPPING_READS = false;

    @Option(doc="If true include duplicate reads, otherwise exclude duplicate reads.")
    public boolean INCLUDE_DUPLICATES = false;


    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(SITES_VCF);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (DOWNSAMPLING_TARGET < 1) throw new IllegalArgumentException("Downsampling target must be > 0.");

        // Setup the various inputs
        final SamReader in = SamReaderFactory.make().open(INPUT);
        final IntervalList sites = buildIntervalList(SITES_VCF, in.getFileHeader().getSequenceDictionary());
        final SamLocusIterator iterator = buildLocusIterator(in, sites);
        final Iterator<VariantContext> vcfIterator = new VCFFileReader(SITES_VCF, false).iterator();

        // Object needed for creating the output report
        final String sample = in.getFileHeader().getReadGroups().iterator().next().getSample();
        final FormatUtil fmt = new FormatUtil();
        final PrintWriter out = new PrintWriter(IOUtil.openFileForBufferedWriting(OUTPUT));
        out.append("chrom\tposition\tid\tref_allele\talt_allele\tsample\tdepth\tref_count\talt_count\tother_count\tmaf\n");

        for (final LocusInfo info : iterator) {
            final VariantContext ctx = nextSnp(vcfIterator);
            if (!info.getSequenceName().equals(ctx.getContig()) || info.getPosition() != ctx.getStart()) {
                throw new IllegalStateException("BAM and VCF out of sync. Bam @ " + info.toString() +
                                                ", VCF @ " + ctx.getContig() + ":" + ctx.getStart());
            }

            final byte refAllele = ctx.getReference().getBases()[0];
            final byte altAllele = ctx.getAlternateAllele(0).getBases()[0];

            List<RecordAndOffset> records = new ArrayList<>(info.getRecordAndPositions());
            if (!ALLOW_OVERLAPPING_READS) records = filterForOverlaps(records);
            records = downsample(records, DOWNSAMPLING_TARGET);
            int refCount=0, altCount=0, otherCount=0;

            for (final RecordAndOffset rec : records) {
                final byte base = rec.getReadBase();
                if      (base == refAllele) ++refCount;
                else if (base == altAllele) ++altCount;
                else                        ++otherCount;
            }

            final double maf = altCount / (double) (altCount + refCount);

            out.append(info.getSequenceName()).append('\t');
            out.append(fmt.format(info.getPosition())).append('\t');
            out.append(ctx.getID()).append('\t');
            out.append((char) refAllele).append('\t');
            out.append((char) altAllele).append('\t');
            out.append(sample).append('\t');
            out.append(fmt.format(records.size())).append('\t');
            out.append(fmt.format(refCount)).append('\t');
            out.append(fmt.format(altCount)).append('\t');
            out.append(fmt.format(otherCount)).append('\t');
            out.append(fmt.format(maf));
            out.append('\n');
        }

        out.close();
        return 0;
    }

    /**
     * Ensures that the returned list only contains one record per read-name, effectively
     * eliminating any double-counting where reads from the same template overlap the same position.
     */
    private List<RecordAndOffset> filterForOverlaps(List<RecordAndOffset> records) {
        Collections.shuffle(records);
        final List<RecordAndOffset> retval = new ArrayList<>();
        final Set<String> readNames = new HashSet<>();

        for (final RecordAndOffset rec : records) {
            if (readNames.add(rec.getRecord().getReadName())) {
                retval.add(rec);
            }
        }

        return retval;
    }

    /** Builds a SamLocusIterator with appropriate filtering. */
    private SamLocusIterator buildLocusIterator(SamReader in, IntervalList sites) {
        final SamLocusIterator iterator = new SamLocusIterator(in, sites);
        iterator.setEmitUncoveredLoci(true);
        iterator.setMappingQualityScoreCutoff(MINIMUM_MAPPING_QUALITY);
        iterator.setQualityScoreCutoff(MINIMUM_BASE_QUALITY);

        final List<SamRecordFilter> filters = new ArrayList<>();
        if (!INCLUDE_DUPLICATES) filters.add(new DuplicateReadFilter());
        iterator.setSamFilters(filters);
        return iterator;
    }

    /** Builds an interval list out of the entries in a VCF. */
    IntervalList buildIntervalList(final File vcf, final SAMSequenceDictionary dict) {
        final VCFFileReader in = new VCFFileReader(vcf, false);
        final SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(dict);
        final IntervalList ilist = new IntervalList(header);

        for (final VariantContext vc : in) {
            if (!vc.isSNP()) continue;
            final Interval i = new Interval(vc.getContig(), vc.getStart(), vc.getEnd(), false, vc.getID());
            ilist.add(i);
        }

        return ilist.sorted();
    }

    /** Consumes from the iterator until a SNP site is found, which is then returned. Null if no more SNPs. */
    VariantContext nextSnp(final Iterator<VariantContext> iterator) {
        while (iterator.hasNext()) {
            final VariantContext ctx = iterator.next();
            if (ctx.isSNP()) return ctx;
        }

        return null;
    }


    /**
     * Downsamples a list of RecordAndOffset to a given coverage. If the list is already below the requested
     * size, it is returned as is.
     */
    <T> List<T> downsample(final List<T> in, final int target) {
        if (in.size() <= target) return in;

        final List<T> out = new ArrayList<>(in);
        Collections.shuffle(out);
        return out.subList(0, target);
    }
}
