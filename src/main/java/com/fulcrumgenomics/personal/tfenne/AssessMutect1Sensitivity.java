package com.fulcrumgenomics.personal.tfenne;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.cmdline.programgroups.VcfOrBcf;
import picard.util.TabbedTextFileWithHeaderParser;
import picard.util.TabbedTextFileWithHeaderParser.Row;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.stream.Collectors;

/**
 *
 */
@CommandLineProgramProperties(
        usageShort = "Assesses the sensitivity of MuTect 1 on a synthetic mixture sample",
        usage      = "Assesses the sensitivity of MuTect 1 on a synthetic mixture sample",
        programGroup = VcfOrBcf.class
)
public class AssessMutect1Sensitivity extends CommandLineProgram {
    @Option(doc="The MuTect callstats file that is to be assessed.")
    public File TEST_CALLSTATS;

    @Option(doc="The tumor BAM file used in calling.")
    public File BAM;

    @Option(doc="The Truth VCF, which contains genotype data for at least the tumor, with the AF genotype field.")
    public File TRUTH_VCF;

    @Option(doc="The name of the tumor sample with the TRUTH_VCF. Can be omitted if the VCF has only a single sample.", optional=true)
    public String TUMOR;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Basename of output files to write.")
    public File OUTPUT;

    @Option(doc="The label to output in the generated table for this dataset.")
    public String LABEL;

    @Option(doc="Zero or more MuTect filters to ignore when consuming the call stats file.")
    public Set<String> IGNORE_FILTER = new HashSet<>();

    private final Log log = Log.getInstance(AssessMutect1Sensitivity.class);

    static class MutationResult {
        VariantContext ctx;
        double alleleFraction;
        boolean called = false;
        boolean presentInCallStats = false;
        Set<String> rejectionReasons = null;
        int readDepth;
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(TEST_CALLSTATS);
        IOUtil.assertFileIsReadable(TRUTH_VCF);

        ////////////////////////////////////////////////////////////////////////
        // Go through the mutect calls and grab all the passing ones
        ////////////////////////////////////////////////////////////////////////
        final Map<String,Set<String>> mutationToRejectionReasons = new HashMap<>();
        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(TEST_CALLSTATS);
        for (final Row row : parser) {
            final String chrom = row.getField("contig");
            final int    pos   = row.getIntegerField("position");
            final String refAllele = row.getField("ref_allele");
            final String altAllele = row.getField("alt_allele");

            final Set<String> filters = new HashSet<String>();
            final String filterString = row.getField("failure_reasons");
            if (filterString != null && !filterString.isEmpty()) {
                filters.addAll(CollectionUtil.makeSet(row.getField("failure_reasons").split(",")));
            }

            mutationToRejectionReasons.put(makeMutationString(chrom, pos, refAllele, altAllele), filters);
        }
        parser.close();

        ////////////////////////////////////////////////////////////////////////
        // Then go through the VCF and count up what we called and missed
        ////////////////////////////////////////////////////////////////////////
        final Histogram<Double> mutationsByAf       = new Histogram<>();
        final Histogram<Double> calledMutationsByAf = new Histogram<>();
        final List<MutationResult> results  = new ArrayList<>();

        final VCFFileReader in = new VCFFileReader(TRUTH_VCF);
        final SamReader bam  = SamReaderFactory.make().open(BAM);
        final String sample = determineTumorSampleName(in);

        for (final VariantContext ctx : in) {
            if (ctx.isIndel() || ctx.isFiltered()) continue; // No Indels for MuTect!

            final String key = makeMutationString(ctx.getContig(), ctx.getStart(), ctx.getReference().getBaseString(), ctx.getAlternateAllele(0).getBaseString());
            final double af  = Double.parseDouble((String) ctx.getGenotype(sample).getExtendedAttribute("AF"));
            final Set<String> rejections = mutationToRejectionReasons.get(key);

            final MutationResult result = new MutationResult();
            result.ctx = ctx;
            result.presentInCallStats = rejections != null;
            result.called = rejections != null && IGNORE_FILTER.containsAll(rejections);
            result.alleleFraction = af;
            result.readDepth = calculateDepthAtLocus(ctx, bam);
            result.rejectionReasons = rejections;
            results.add(result);

            mutationsByAf.increment(af);
            calledMutationsByAf.increment(af, result.called ? 1 : 0);
        }
        in.close();

        ////////////////////////////////////////////////////////////////////////
        // Write out the results
        ////////////////////////////////////////////////////////////////////////
        try {
            // First the summary level
            BufferedWriter out = IOUtil.openFileForBufferedWriting(new File(OUTPUT.getParentFile(), OUTPUT.getName() + ".sensitivity_summary.txt"));
            final NumberFormat dfmt = new DecimalFormat("0.00000");
            final NumberFormat ifmt = new DecimalFormat("0");
            out.append("condition\tallele_fraction\ttotal_mutations\tcalled_mutations\tsensitivity\n");
            for (final Double af : mutationsByAf.keySet()) {
                final int totalCount  = (int) mutationsByAf.get(af).getValue();
                final int calledCount = (int) calledMutationsByAf.get(af).getValue();
                final double sensitivity = calledCount / (double) totalCount;

                out.append(LABEL).append('\t');
                out.append(dfmt.format(af)).append('\t');
                out.append(String.valueOf(totalCount)).append('\t');
                out.append(String.valueOf(calledCount)).append('\t');
                out.append(dfmt.format(sensitivity));
                out.newLine();
            }

            // Output a total(s) line
            final double called = calledMutationsByAf.getSumOfValues();
            final double total  = mutationsByAf.getSumOfValues();
            out.append(LABEL).append('\t');
            out.append("all").append('\t');
            out.append(String.valueOf((int) total)).append('\t');
            out.append(String.valueOf((int) called)).append('\t');
            out.append(dfmt.format(called/total));
            out.newLine();
            out.close();

            // Then the details
            out = IOUtil.openFileForBufferedWriting(new File(OUTPUT.getParentFile(), OUTPUT.getName() + ".sensitivity_details.txt"));
            out.append("condition\tchrom\tpos\tmutation\tref_allele\talt_allele\tmaf\tin_callstats\tcalled\tcoverage\trejection_reasons\n");
            for (final MutationResult result : results) {
                out.append(LABEL).append('\t');
                out.append(result.ctx.getContig()).append('\t');
                out.append(ifmt.format(result.ctx.getStart())).append('\t');
                out.append(result.ctx.getID()).append('\t');
                out.append(result.ctx.getReference().getBaseString()).append('\t');
                out.append(result.ctx.getAlternateAllele(0).getBaseString()).append('\t');
                out.append(dfmt.format(result.alleleFraction)).append('\t');
                out.append(result.presentInCallStats ? "Yes" : "No").append('\t');
                out.append(result.called ? "Yes" : "No").append('\t');
                out.append(ifmt.format(result.readDepth)).append('\t');
                out.append(result.rejectionReasons == null ? "" : result.rejectionReasons.stream().collect(Collectors.joining(",")));
                out.newLine();
            }

            out.close();

        }
        catch (IOException ioe) {
            throw new RuntimeIOException(ioe);
        }

        return 0;
    }

    int calculateDepthAtLocus(final VariantContext ctx, final SamReader in) {
        final SAMRecordIterator iterator = in.queryOverlapping(ctx.getContig(), ctx.getStart(), ctx.getStart());
        final Set<String> readNames = new HashSet<>();
        while (iterator.hasNext()) {
            final SAMRecord rec = iterator.next();
            if (!rec.getReadUnmappedFlag() && !rec.getDuplicateReadFlag()) readNames.add(rec.getReadName());
        }
        iterator.close();
        return readNames.size();
    }


    /**
     * Returns the given tumor name if supplied, else the single sample name from the VCF,
     * else throws an exception if the VCF is multi-sample.
     */
    private String determineTumorSampleName(VCFFileReader in) {
        if (TUMOR != null) {
            return TUMOR;
        }
        else if (in.getFileHeader().getNGenotypeSamples() == 1) {
            return in.getFileHeader().getGenotypeSamples().get(0);
        }
        else {
            throw new IllegalStateException("TUMOR must be provided when the TRUTH_VCF contains multiple samples.");
        }
    }

    /** Makes a String that should be unique for any given SNV. */
    String makeMutationString(final String contig, final int pos, final String refAllele, final String altAllele) {
        return contig + ":" + pos + ":" + refAllele + ":" +  ":" + altAllele;
    }
}
