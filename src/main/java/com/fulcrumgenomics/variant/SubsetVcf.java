package com.fulcrumgenomics.variant;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.programgroups.VcfOrBcf;
import picard.vcf.ByIntervalListVariantContextIterator;

import java.io.File;

/**
 * Trivial command line program to produce a new VCF from an input VCF and a set of intervals.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage="Subsets a VCF by genomic region.",
        usageShort="Subsets a VCF by genomic region.",
        programGroup = VcfOrBcf.class
)
public class SubsetVcf extends CommandLineProgram {
    @Option(shortName="I", doc="Input VCF.")
    public File INPUT;

    @Option(shortName="O", doc="Output VCF.")
    public File OUTPUT;

    @Option(shortName="L", doc="Intervals to subset to.")
    public File INTERVALS;

    @Option(shortName="P", doc="Pad intervals by n bases prior to subsetting.")
    public int PADDING = 0;

    @Override protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(INTERVALS);
        IOUtil.assertFileIsWritable(OUTPUT);

        final IntervalList intervals = IntervalList.fromFile(INTERVALS).padded(PADDING, PADDING).uniqued();
        final VCFFileReader in = new VCFFileReader(INPUT, false);
        final ByIntervalListVariantContextIterator iterator = new ByIntervalListVariantContextIterator(in, intervals);

        final VariantContextWriter out = new VariantContextWriterBuilder().setOutputFile(OUTPUT)
                .setReferenceDictionary(in.getFileHeader().getSequenceDictionary())
                .build();
        out.writeHeader(in.getFileHeader());

        while (iterator.hasNext()) {
            out.add(iterator.next());
        }

        out.close();
        return 0;
    }
}
