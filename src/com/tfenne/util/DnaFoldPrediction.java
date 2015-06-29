package com.tfenne.util;

/**
 *
 */
public class DnaFoldPrediction {
    private final String sequence;
    private final String structure;
    private final double deltaG;

    public DnaFoldPrediction(String sequence, String structure, double deltaG) {
        this.sequence = sequence;
        this.structure = structure;
        this.deltaG = deltaG;
    }

    public String sequence() { return this.sequence; }
    public String structure() { return this.structure; }
    public double deltaG() { return this.deltaG; }
}
