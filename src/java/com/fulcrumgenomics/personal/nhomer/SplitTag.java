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

package com.fulcrumgenomics.personal.nhomer;

import com.fulcrumgenomics.cmdline.Utilities;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.IlluminaUtil;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(
        usage = "Splits an optional tag in a SAM or BAM into multple optional tags.",
        usageShort = "Splits an optional tag in a SAM or BAM into multple optional tags.",
        programGroup = Utilities.class
)
public class SplitTag extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input SAM or BAM.")
    public File INPUT;

    @Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output SAM or BAM.")
    public File OUTPUT;

    @Option(doc = "The tag to split.")
    public String INPUT_TAG;

    @Option(doc = "The output tags.  There should be one per token.")
    public List<String> OUTPUT_TAGS;

    @Option(doc = "The delimiter used to split the string.")
    public String DELIMITER = IlluminaUtil.BARCODE_DELIMITER;

    @Override
    protected String[] customCommandLineValidation() {
        final List<String> errors = new ArrayList<>();
        if (2 != INPUT_TAG.length()) {
            errors.add("INPUT_TAG must be of length two (was " + INPUT_TAG.length() + ").");
        }
        for (final String outputTag : OUTPUT_TAGS) {
            if (2 != outputTag.length()) {
                errors.add("Tags in OUTPUT_TAGS must be of length two (was " + outputTag.length() + ").");
            }
        }
        if (errors.isEmpty()) return super.customCommandLineValidation();
        else return errors.toArray(new String[errors.size()]);
    }

    @Override
    public int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        final SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(reader.getFileHeader(), true, OUTPUT);

        for (final SAMRecord record : reader) {
            final String value = record.getStringAttribute(INPUT_TAG);
            if (value == null) throw new PicardException(String.format("Record '%s' was missing the tag '%s'", record.getReadName(), INPUT_TAG));
            final String[] tokens = value.split(DELIMITER);
            if (tokens.length != OUTPUT_TAGS.size()) throw new PicardException(String.format("Record '%s' did not have '%d' tokens", record.getReadName(), OUTPUT_TAGS.size()));
            for (int i = 0; i < tokens.length; i++) {
                record.setAttribute(OUTPUT_TAGS.get(i), tokens[i]);
            }
            writer.addAlignment(record);
        }

        CloserUtil.close(reader);
        writer.close();

        return 0;
    }
}
