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
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

@CommandLineProgramProperties(
        usage = "Re-assigns existing read groups to new read groups.",
        usageShort = "Re-assigns existing read groups to new read groups.",
        programGroup = Utilities.class
)
public class RemapReadGroups extends CommandLineProgram {

    @Option(shortName= StandardOptionDefinitions.INPUT_SHORT_NAME, doc="Input file (bam or sam).")
    public File INPUT = null;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file (bam or sam).")
    public File OUTPUT = null;

    @Option(doc="The read group ID mapping csv file (RGID_before,RGID_after,LB,PL,PU,SM,CN,DS,DT,PI).  Empty nor non-existent columns are ignored.")
    public File READ_GROUP_ID_MAPPING = null;

    @Option(doc="Keep old read group info tags if new read group info tags are not present", optional=true)
    public boolean KEEP_OLD_READ_GROUP_ATTRIBUTES = false;

    private final Log log = Log.getInstance(RemapReadGroups.class);

    /** Required main method implementation. */
    public static void main(final String[] anewReadGroupv) {
        new RemapReadGroups().instanceMainWithExit(anewReadGroupv);
    }


    private Map<String,SAMReadGroupRecord> getReadGroupMap(final File input) {


        final Map<String,SAMReadGroupRecord> readGroupMap = new HashMap<>();
        final LineReader reader;
        try {
            reader = new BufferedLineReader(new FileInputStream(input));
        }
        catch (final IOException e) {
            throw new RuntimeException(e);
        }
        String line;
        while (null != (line = reader.readLine())) {
            final String[] tokens = line.split(",", 2);
            if (tokens.length != 2) {
                throw new PicardException("Incorrect number of columns (" + tokens.length + ") in : " + line);
            }
            final SAMTextHeaderCodec codec = new SAMTextHeaderCodec();
            final LineReader lineReader = new StringLineReader(tokens[1]);
            final SAMFileHeader header = codec.decode(lineReader, input.getAbsolutePath());
            final List<SAMReadGroupRecord> readGroupRecords = header.getReadGroups();
            if (0 == readGroupRecords.size()) throw new PicardException("Could not find a read group in item: " + tokens[1]);
            if (1 < readGroupRecords.size()) throw new PicardException("Multiple read groups found in item: " + tokens[1]);
            lineReader.close();
            final SAMReadGroupRecord rg = readGroupRecords.get(0);
            readGroupMap.put(tokens[0], rg);
        }
        reader.close();
        return readGroupMap;
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(READ_GROUP_ID_MAPPING);

        final Map<String,SAMReadGroupRecord> readGroupMap = getReadGroupMap(READ_GROUP_ID_MAPPING);

        final SamReader reader = SamReaderFactory.makeDefault().referenceSequence(REFERENCE_SEQUENCE).open(INPUT);

        final SAMFileHeader inHeader = reader.getFileHeader();
        final SAMFileHeader outHeader = inHeader.clone();

        // Create new read groups
        final List<SAMReadGroupRecord> newReadGroupsList = new ArrayList<>();
        for (final SAMReadGroupRecord oldReadGroup : inHeader.getReadGroups()) { // original records
            final String readGroupID = oldReadGroup.getReadGroupId(); // old read group
            if (!readGroupMap.containsKey(readGroupID)) {
                throw new PicardException("Read group does not contain key: " + readGroupID);
            }
            // get the new read group
            final SAMReadGroupRecord newReadGroup = readGroupMap.get(readGroupID);
            // keep, if desired, the attributes from the old read group
            if (KEEP_OLD_READ_GROUP_ATTRIBUTES) {
                // NB: read group ID is not an attribute, so we should be ok copying over all attributes
                for (final Map.Entry<String, String> entry : oldReadGroup.getAttributes()) {
                    if (null == newReadGroup.getAttribute(entry.getKey())) {
                        newReadGroup.setAttribute(entry.getKey(), entry.getValue());
                    }
                }
            }
            newReadGroupsList.add(newReadGroup);
            log.info(String.format("Created read group OLDID=%s ID=%s PL=%s LB=%s SM=%s",
                    readGroupID, newReadGroup.getId(), newReadGroup.getPlatform(),
                    newReadGroup.getLibrary(), newReadGroup.getSample()));
        }

        // create the new header and output file
        outHeader.setReadGroups(newReadGroupsList);
        final SAMFileWriter outWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader,
                outHeader.getSortOrder() == inHeader.getSortOrder(),
                OUTPUT);

        // now modify the records
        final ProgressLogger progress = new ProgressLogger(log);
        for (final SAMRecord read : reader) {
            // get the current read group
            final SAMReadGroupRecord readGroup = read.getReadGroup();
            if (null != readGroup) {
                // get the current read group id
                final String readGroupID = readGroup.getReadGroupId();
                if (null != readGroupID) {
                    // update the read group tag in the SAM record
                    if (!readGroupMap.containsKey(readGroupID)) {
                        throw new PicardException(String.format("Did not find read group id '%s' in read group map", readGroupID));
                    }
                    read.setAttribute(SAMTag.RG.name(), readGroupMap.get(readGroupID).getReadGroupId());
                }
            }
            outWriter.addAlignment(read);
            progress.record(read);
        }

        // cleanup
        CloserUtil.close(reader);
        outWriter.close();

        return 0;
    }
}
