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
package com.fulcrumgenomics.personal.tfenne;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.*;

/**
 * Creates a VCF by mixing two germline samples at a given proportion.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        usage="Creates a new VCF for a 'tumor' and a 'normal' where the tumor is a mixture of\n" +
                "two germline samples, and the normal is one of the same germline samples. Only\n" +
                "outputs loci that are variant in the 'tumor'. Writes the 'AF' attribute in the\n" +
                "format field with the expected allele fraction in the tumor.",
        usageShort = "Creates a VCF by in-silico mixing data from two samples."
)
public class MakeGermlineMixtureVcf extends CommandLineProgram {
    @Option(shortName=StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input VCF file.")
    public File INPUT;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output VCF file to write.")
    public File OUTPUT;

    @Option(doc="The name of the 'tumor'  sample.")
    public String TUMOR;

    @Option(doc="The name of the 'normal' sample.")
    public String NORMAL;

    @Option(doc="What fraction of the mixture is from the tumor sample?")
    public double TUMOR_FRACTION = 0.5d;

    @Option(doc="A set of intervals to include variants from.")
    public File INTERVALS;

    public static final String ALLELE_FRACTION_FIELD = "AF";
    public static final String GERMLINE_FILTER = "alt_allele_in_normal";

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<>();
        if (TUMOR_FRACTION <= 0) errors.add("TUMOR_FRACTION must be greater than zero.");
        if (TUMOR_FRACTION >= 1) errors.add("TUMOR_FRACTION must be less than one.");
        if (TUMOR.equals(NORMAL)) errors.add("TUMOR and NORMAL must be different samples.");

        if (errors.isEmpty()) return null;
        else return errors.toArray(new String[errors.size()]);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(INTERVALS);

        final VCFFileReader in = new VCFFileReader(INPUT, true);
        {
            final Collection<String> samples = in.getFileHeader().getSampleNamesInOrder();
            if (!samples.contains(TUMOR))  throw new IllegalStateException("VCF does not contain sample: " + TUMOR);
            if (!samples.contains(NORMAL)) throw new IllegalStateException("VCF does not contain sample: " + NORMAL);
        }

        final VariantContextWriter out = new VariantContextWriterBuilder().setOutputFile(OUTPUT).setReferenceDictionary(in.getFileHeader().getSequenceDictionary()).build();
        out.writeHeader(buildOutputHeader(in.getFileHeader(), "tumor", "normal"));

        final IntervalList intervals = IntervalList.fromFile(INTERVALS).uniqued();
        for (final Interval interval : intervals) {
            final Iterator<VariantContext> iterator = in.query(interval.getContig(), interval.getStart(), interval.getEnd());
            while (iterator.hasNext()) {
                final VariantContext ctx = iterator.next();
                final Genotype tumor  = ctx.getGenotype(TUMOR);
                final Genotype normal = ctx.getGenotype(NORMAL);

                if (ctx.isFiltered() || ctx.getAlternateAlleles().size() > 1 || tumor.isHomRef()) continue;

                final Allele refAllele = ctx.getReference();
                final Allele altAllele = ctx.getAlternateAllele(0);

                // Build the Tumor genotype
                final GenotypeBuilder tumorBuilder = new GenotypeBuilder("tumor");
                final List<Allele> tumorAlleles = new ArrayList<>();
                if (tumor.countAllele(refAllele) + normal.countAllele(refAllele) > 0) tumorAlleles.add(refAllele);
                if (tumor.countAllele(altAllele) + normal.countAllele(altAllele) > 0) tumorAlleles.add(altAllele);
                if (tumorAlleles.size() == 1) tumorAlleles.add(tumorAlleles.get(0));

                tumorBuilder.alleles(tumorAlleles);
                final double af =  (TUMOR_FRACTION    * tumor.countAllele(altAllele) / 2)
                                + ((1-TUMOR_FRACTION) * normal.countAllele(altAllele) / 2);
                tumorBuilder.attribute(ALLELE_FRACTION_FIELD, af);

                // Build the normal genotype
                final GenotypeBuilder normalBuilder = new GenotypeBuilder("normal");
                normalBuilder.alleles(normal.getAlleles());

                // Build and emit the variant context
                VariantContextBuilder builder = new VariantContextBuilder(ctx.getSource(), ctx.getContig(), ctx.getStart(), ctx.getEnd(), ctx.getAlleles());
                builder.genotypes(tumorBuilder.make(), normalBuilder.make());
                if (!normal.isHomRef()) builder.filter(GERMLINE_FILTER);

                out.add(builder.make());
            }
        }

        out.close();
        return 0;
    }

    /** Builds a header that can be used to write the mixture VCF. */
    VCFHeader buildOutputHeader(final VCFHeader in, String... samples) {
        final VCFHeader out = new VCFHeader(Collections.emptySet(), Arrays.asList(samples));

        in.getFilterLines().stream().forEach(out::addMetaDataLine);
        out.addMetaDataLine(new VCFFilterHeaderLine(GERMLINE_FILTER, "Evidence seen in the normal sample"));

        in.getFormatHeaderLines().stream().forEach(out::addMetaDataLine);
        out.addMetaDataLine(new VCFFormatHeaderLine("AF", 1, VCFHeaderLineType.Float, "Alt allele fraction in the tumor"));

        in.getInfoHeaderLines().stream().forEach(out::addMetaDataLine);
        out.setSequenceDictionary(in.getSequenceDictionary());
        return out;
    }
}
