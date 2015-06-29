package com.tfenne.cmdline;

import picard.cmdline.CommandLineProgramGroup;

/**
 * Program group/type for utility programs.
 */
public class Utilities implements CommandLineProgramGroup {
    @Override public String getName() { return "Utilities"; }
    @Override public String getDescription() { return "Various utility programs."; }
}

