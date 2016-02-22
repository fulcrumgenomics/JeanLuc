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

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Class that wraps ViennaRNAs RNAFold utility which can be used to estimate the minimum free energy
 * secondary structure of DNA and RNA molecules. When constructed a background process is started
 * running RNAFold, calls to predict() then pipe input and output through RNAFold.
 *
 * @author Tim Fennell
 */
public class DnaFoldPredictor {
    private final Process process;
    private final BufferedWriter out;
    private final BufferedReader in;

    /**
     * Constructs an instance that has a running copy of ViennaRNA's RNAFold executable behind it
     * that can then be interactively pushed a new sequence and get back folding information.
     *
     * @param viennaRnaInstallDir the parent dir of Vienna RNA's bin, lib, share directories
     * @param tm the tm at which the calculations should be performed
     */
    public DnaFoldPredictor(final File viennaRnaInstallDir, final double tm) {
        final File binary = new File(viennaRnaInstallDir, "bin/RNAfold");
        final File params = new File(viennaRnaInstallDir, "share/ViennaRNA/dna_mathews2004.par");
        final List<String> args = new ArrayList<>();
        args.add(binary.getAbsolutePath());
        args.add("--noconv");
        args.add("--noPS");
        args.add("--temp=" + tm);
        args.add("--paramFile=" + params.getAbsolutePath());

        try {
            final ProcessBuilder builder = new ProcessBuilder(args);
            builder.redirectError(ProcessBuilder.Redirect.INHERIT);
            this.process = builder.start();
            this.out = new BufferedWriter(new OutputStreamWriter(this.process.getOutputStream()), Defaults.BUFFER_SIZE);
            this.in = new BufferedReader(new InputStreamReader(this.process.getInputStream()), Defaults.BUFFER_SIZE);
        }
        catch (IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    /** Creates the secondary structure prediction for the input sequence. */
    public DnaFoldPrediction predict(final String sequence) {
        try {
            this.out.write(sequence);
            this.out.newLine();
            this.out.flush();

            // Output when using stdin/stdout looks like this:
            // TGAACTCCTCAACCCTCTTCTCATCAGGAGTGATAGTGGCACATTTGACG
            // ((.(((..(..(((((.(....)..))).)).).)))..))......... ( -3.13)

            final String seq2 = readLine();
            final String result = readLine();

            if (!seq2.equals(sequence)) {
                throw new IllegalStateException("Return sequence '" + seq2 + "' does not match entered sequence '" + sequence + "'");
            }

            final String structure = result.substring(0, seq2.length());
            final double dg = Double.parseDouble(result.substring(seq2.length()).replace("(", "").replace(")", ""));

            return new DnaFoldPrediction(seq2, structure, dg);
        }
        catch (IOException ioe) {
            throw new RuntimeIOException(ioe);
        }
    }

    /** Reads a line from the input, ignoring warning lines. */
    private String readLine() throws IOException {
        while (true) {
            final String line = this.in.readLine();
            if (line == null || !line.startsWith("WARNING")) return line;
        }
    }

    /** Kills the underlying process and makes future calls to predict() invalid. */
    public void close() {
        CloserUtil.close(this.out);
        if (this.process.isAlive()) this.process.destroy();
    }
}
