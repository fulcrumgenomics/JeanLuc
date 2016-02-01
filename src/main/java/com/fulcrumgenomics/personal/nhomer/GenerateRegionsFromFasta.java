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

import com.fulcrumgenomics.cmdline.Personal;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.IOUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;

@CommandLineProgramProperties(
        usage = "Generates a list of freebayes/bamtools region specifiers on stdout.",
        usageShort = "Generates a list of freebayes/bamtools region specifiers on stdout.",
        programGroup = Personal.class
)
public class GenerateRegionsFromFasta extends CommandLineProgram {

    @Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "Input FASTA.")
    public File INPUT;

    @Option(shortName = "W", doc = "The size of regions to output.")
    public int REGION_SIZE = 100000;

    @Override
    public int doWork() {
        IOUtil.assertFileIsReadable(INPUT);

        final ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(INPUT);
        final SAMSequenceDictionary sequenceDictionary = referenceSequenceFile.getSequenceDictionary();

        for (int i = 0; i < sequenceDictionary.size(); i++) {
            final SAMSequenceRecord sequenceRecord = sequenceDictionary.getSequence(i);
            int regionStart = 0; // yes, this must be zero-based
            while (regionStart < sequenceRecord.getSequenceLength()) {
                final int regionEnd = Math.min(regionStart + REGION_SIZE, sequenceRecord.getSequenceLength());
                System.out.println(sequenceRecord.getSequenceName() + ":" + regionStart + "-" + regionEnd);
                regionStart = regionEnd;
            }
        }

        return 0;
    }
}
