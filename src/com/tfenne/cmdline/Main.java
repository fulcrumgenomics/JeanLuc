package com.tfenne.cmdline;

import htsjdk.samtools.util.CollectionUtil;
import picard.cmdline.PicardCommandLine;

/**
 * Simple class that causes uses the Picard system for command line parsing and sets
 * us up to have multiple programs accessible.
 *
 * @author Tim Fennell
 */
public class Main extends PicardCommandLine {
    public static void main(final String[] args) {
        System.exit(new Main().instanceMain(args, CollectionUtil.makeList("com.tfenne"), "Main"));
    }
}
