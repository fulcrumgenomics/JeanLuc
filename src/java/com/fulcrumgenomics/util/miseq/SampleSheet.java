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

package com.fulcrumgenomics.util.miseq;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import picard.PicardException;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Stores information about samples from an Illumina Experiment Manager Sample Sheet (typically a MiSeq).  The
 * samples may also include derived information.  If library id is not specified, it will be set to the sample name.
 *
 * @author Nils Homer
 */
public class SampleSheet implements Iterable<Sample> {

    /** Required header names for standard fields in the Illumina Experiment Manager Sample Sheet. */
    public static final String SAMPLE_ID      = "Sample_ID";
    public static final String SAMPLE_NAME    = "Sample_Name";
    public static final String LIBRARY_ID     = "Library_ID";
    public static final String SAMPLE_PROJCET = "Sample_Project";
    public static final String DESCRIPTION    = "Description";

    /** Optional header names for standard fields in the Illumina Experiment Manager Sample Sheet. */
    public static final String R1_BARCODE_BASES = "R1_Bases";
    public static final String R2_BARCODE_BASES = "R2_Bases";
    public static final String I7_BASES         = "I7_Bases";
    public static final String I5_BASES         = "I5_Bases";

    private List<Sample> samples = new ArrayList<>();

    /** Create a sample sheet from the given file (typically 'SampleSheet.csv') */
    public SampleSheet(final File sampleSheet) {
        final List<Map<String, String>> sampleData = getSampleData(sampleSheet);

        // Create and add samples
        int sampleOrdinal = 1;
        for (final Map<String, String> sampleDatum : sampleData) {
            this.samples.add(makeSample(sampleOrdinal, sampleDatum));
            sampleOrdinal++;
        }
    }

    /** Reads in the the Sample data only, one per row */
    private List<Map<String, String>> getSampleData(final File sampleSheet) {
        final List<Map<String, String>> sampleData = new ArrayList<>();
        try {
            final BufferedReader in = IOUtil.openFileForBufferedReading(sampleSheet);
            int lineNumber = 0;

            // Skip everything up until the sample section
            while (!in.readLine().startsWith("[Data]")) { lineNumber++; }
            lineNumber++;

            List<String> header = null;
            String line;
            int rowNumber = 1;
            while ((line = in.readLine()) != null) {
                if (header == null) {
                    header = Arrays.asList(line.split(",", -1)); // NB: include trailing empty strings
                }
                else {
                    final List<String> values = Arrays.asList(line.split(",", -1)); // NB: include trailing empty strings
                    // check we have the correct # of columns
                    if (values.size() != header.size()) {
                        throw new PicardException(String.format("# of columns in the header and current row do not match ('%d' != '%d'): row #%d line number #%d",
                                header.size(), values.size(), rowNumber, lineNumber));
                    }
                    // create the map and store it
                    final Map<String, String> map = new HashMap<>();
                    for (int i = 0; i < values.size(); i++) map.put(header.get(i).trim(), values.get(i).trim());
                    sampleData.add(map);
                    rowNumber++;
                }
                lineNumber++;
            }

            in.close();
        }
        catch (final IOException ioe) { throw new RuntimeIOException(ioe); }
        return sampleData;
    }

    /** Attempts to clean special characters from a sample id or name as Illumina's software does. */
    public static String cleanMiseqSampleId(final String id) {
        String out = id;
        // remove '#'
        for (final String ch : new String[] {"#"}) { out = out.replace(ch, ""); }
        // replace '[_+ ]' with '-'
        for (final char   ch : new char[]   {'_', '+', ' '}) { out = out.replace(ch, '-'); }
        return out;
    }

    public Iterator<Sample> iterator() { return this.samples.iterator(); }

    public int size() { return this.samples.size(); }

    public Sample get(final int index) { return this.samples.get(index); }

    /** Creates a sample from the given row data */
    private Sample makeSample(final int sampleOrdinal, final Map<String, String> sampleDatum) {
        final String sampleName = getStringField(sampleDatum, SAMPLE_NAME);
        String libraryId = getStringField(sampleDatum, LIBRARY_ID);
        if (libraryId == null) libraryId = sampleName;
        return new Sample(
                sampleOrdinal,
                getStringField(sampleDatum, SAMPLE_ID),
                sampleName,
                libraryId,
                getStringField(sampleDatum, SAMPLE_PROJCET),
                getStringField(sampleDatum, DESCRIPTION),
                getStringField(sampleDatum, R1_BARCODE_BASES),
                getStringField(sampleDatum, R2_BARCODE_BASES),
                getStringField(sampleDatum, I7_BASES),
                getStringField(sampleDatum, I5_BASES)
        );
    }

    /** Gets the trimmed value of the given key (`name`) from the map.  If empty or not found, returns null */
    private static String getStringField(final Map<String, String> sampleDatum, final String name) {
        String str = sampleDatum.get(name);
        if (null != str) str = str.trim();
        if (null == str || str.isEmpty()) return null;
        return str;
    }
}
