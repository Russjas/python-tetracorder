"""
Configuration settings for python-tetracorder
"""

#=============== pre-defined values for the symbolic thresholds in Tetracorder
#comments are preserved from Tetracorder cmd files
##comments are mine
DEFAULT_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.4, 0.6],
"[GLBLFITVEG1]": [0.2, 0.3],
"[GLBLFITVEG2]": [0.5, 0.7],
"[FDFITveg3]": [0.3, 0.5],
#GLBLFITVEGR is for vegetation red-edge detection, case 1
"[GLBLFITVEGR]": [0.3, 0.40],
#GLBLFITVEGT is for vegetation type detection, case 2
"[GLBLFITVEGT]": [0.2, 0.30],
#GLBLFITVEGT is for vegetation water detection, case 3, 4, 5
"[GLBLFITVEGW]": [0.4, 0.55],
#GLBLDPFITREE is set high to reduce false positives
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

#  Continuum thresholds

#NB: ct, rct, lct do not use fuzzy logic
"[CTHRESH1]": 0.01,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.02,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.04,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.05,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.08,          # ct threshold value (continuum threshold)

##fuzzy logic params for feature tests

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.0002, 0.0004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour - specpr also does this
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region

# Output DN scaling, e.g. 8 DN 255 = 0.5000 for 8-bit scaled 0.5 band depth = 255

"[O8DN50]":   0.5000,
"[O8DN100]":  1.0000,
"[O8DN14um]": 0.5000,
"[O8DN19um]": 0.5000,
"[O8DN27um]": 0.5000,
"[O8DNree]":  0.1000,
"[O8DNree2]": 0.0200,   # used in grolup 21


}
DEFAULT_BAK23_VALUES = {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.4, 0.6],
"[GLBLFITVEG1]": [0.2, 0.3],
"[GLBLFITVEG2]": [0.5, 0.7],
"[FDFITveg3]": [0.3, 0.5],
"[GLBLFITVEGR]": [0.3, 0.40],
"[GLBLFITVEGT]": [0.2, 0.30],
"[GLBLFITVEGW]": [0.4, 0.55],
"[GLBLDPFITREE]": [0.65, 0.80],## These are ommited from the bak23 variable cmd file - added here to avoid leaving symbols in the rule
"[GLBLFITREE]": [0.5, 0.6],## These are ommited from the bak23 variable cmd file - added here to avoid leaving symbols in the rule

"[CTHRESH1]": 0.01,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.02,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.04,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.05,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.08,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.001, 0.002],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour - specpr also does this
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region

# Output DN scaling, e.g. 8 DN 255 = 0.5000 for 8-bit scaled 0.5 band depth = 255

"[O8DN50]":   0.5000,
"[O8DN100]":  1.0000,
"[O8DN14um]": 0.5000,
"[O8DN19um]": 0.5000,
"[O8DN27um]": 0.5000,
# There two are not defined in Bak23, added here for completeness
"[O8DNree]":  0.1000,
"[O8DNree2]": 0.0200,   # used in grolup 21
}

EMIT_C_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.4, 0.6],
"[GLBLFITVEG1]": [0.2, 0.3],
"[GLBLFITVEG2]": [0.5, 0.7],
"[FDFITveg3]": [0.3, 0.5],
"[GLBLFITVEGR]": [0.3, 0.40],
"[GLBLFITVEGT]": [0.2, 0.30],
"[GLBLFITVEGW]": [0.4, 0.55],
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

"[CTHRESH1]": 0.01,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.02,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.04,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.05,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.08,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.0002, 0.0004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour - specpr also does this
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
# Output DN scaling, e.g. 8 DN 255 = 0.5000 for 8-bit scaled 0.5 band depth = 255
"[O8DN50]":   0.5000,
"[O8DN100]":  1.0000,
"[O8DN14um]": 0.5000,
"[O8DN19um]": 0.5000,
"[O8DN27um]": 0.5000,
"[O8DNree]":  0.1000,
"[O8DNree2]": 0.0200,   # used in group 21
}
MMM_09C_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.91, 0.93],
"[GLBLFITVEG1]": [0.91, 0.93],
"[GLBLFITVEG2]": [0.91, 0.93],
"[FDFITveg3]": [0.91, 0.93],
"[GLBLFITVEGR]": [0.88, 0.91],
"[GLBLFITVEGT]": [0.88, 0.91],
"[GLBLFITVEGW]": [0.88, 0.91],
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

"[CTHRESH1]": 0.005,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.01,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.02,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.025,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.04,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.0002, 0.0004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour - specpr also does this
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
# Output DN scaling, e.g. 8 DN 255 = 0.5000 for 8-bit scaled 0.5 band depth = 255
"[O8DN50]":   0.1000,
"[O8DN100]":  1.0000,
"[O8DN14um]": 0.1000,
"[O8DN19um]": 0.1000,
"[O8DN27um]": 0.2000,
"[O8DNree]":  0.1000,
"[O8DNree2]": 0.0200,   # used in group 21
}
MMM255T_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.91, 0.93],
"[GLBLFITVEG1]": [0.91, 0.93],
"[GLBLFITVEG2]": [0.91, 0.93],
"[FDFITveg3]": [0.91, 0.93],
"[GLBLFITVEGR]": [0.88, 0.91],
"[GLBLFITVEGT]": [0.88, 0.91],
"[GLBLFITVEGW]": [0.88, 0.91],
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

"[CTHRESH1]": 0.005,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.01,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.02,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.025,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.04,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.002,  0.004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.0002, 0.0004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.002, 0.004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.001, 0.002],     # r*bd for 1-      micron region
"[RBD14]": [0.0008, 0.002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.0007, 0.0011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.0007, 0.0011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.001, 0.002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour - specpr also does this
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
# Output DN scaling, e.g. 8 DN 255 = 0.5000 for 8-bit scaled 0.5 band depth = 255
"[O8DN50]":   0.1000,
"[O8DN100]":  1.0000,
"[O8DN14um]": 0.1000,
"[O8DN19um]": 0.1000,
"[O8DN27um]": 0.2000,
"[O8DNree]":  0.1000,
"[O8DNree2]": 0.0200,   # used in group 21

}
HYB2RYUG_VALUES= {
"[GLBLFITALL]": [0.2, 0.3],
"[GLBLFDFIT]": [0.3, 0.4],
"[GLBLDPFIT]": [0.5, 0.6],
"[GLBLDPFITg2]": [0.65, 0.7],
"[GLBLFITVEG]": [0.91, 0.93],
"[GLBLFITVEG1]": [0.91, 0.93],
"[GLBLFITVEG2]": [0.91, 0.93],
"[FDFITveg3]": [0.91, 0.93],
"[GLBLFITVEGR]": [0.88, 0.91],
"[GLBLFITVEGT]": [0.88, 0.91],
"[GLBLFITVEGW]": [0.88, 0.91],
"[GLBLDPFITREE]": [0.65, 0.80],
"[GLBLFITREE]": [0.5, 0.6],

"[CTHRESH1]": 0.001,          # ct threshold value (continuum threshold)
"[CTHRESH2]": 0.002,         # ct threshold value (continuum threshold)
"[CTHRESH4]": 0.002,          # ct threshold value (continuum threshold)
"[CTHRESH5]": 0.0025,          # ct threshold value (continuum threshold)
"[CTHRESH8]": 0.004,          # ct threshold value (continuum threshold)

"[RBD1b]": [0.0002,  0.0004],     # r*bd for 1-2 micron region broad bands
"[RBDree]": [0.00002, 0.00004],    # r*bd for Vis-NIR region REE bands
"[RBDvrs]": [0.0002, 0.0004],    # r*bd for Veg red-edge shift
"[RBD10]": [0.0001, 0.00002],     # r*bd for 1-      micron region ##sic Tetracorder but right<left could be typo
"[RBD14]": [0.00008, 0.0002],     # r*bd for 1.3-1.5 micron region
"[RBD22]": [0.00007, 0.00011],    # r*bd for 2.2-    micron region
"[RBD23]": [0.00007, 0.00011],    # r*bd for 2.3-    micron region
"[RBD25]": [0.0001, 0.0002],     # r*bd for 2.5     micron region
##sic Tetracorder repeats key. Preserved here in text
##Python behaviour will overwrite first with second
##unsure which is preserved by specpr behaviour - specpr also does this
"[RBD25]": [0.00007, 0.00011],    # r*bd for 2.5-    micron region #sic Tetracorder
"[RBD30]": [0.00005, 0.0001],     # r*bd for 3-      micron region
"[RBD35]": [0.00005, 0.0001],     # r*bd for 3.5-    micron region
# Output DN scaling, e.g. 8 DN 255 = 0.5000 for 8-bit scaled 0.5 band depth = 255
"[O8DN50]":   0.1000,
"[O8DN100]":  1.0000,
"[O8DN14um]": 0.2000,
"[O8DN19um]": 0.2000,
"[O8DN27um]": 0.5000,
"[O8DNree]":  0.1000,
"[O8DNree2]": 0.0200,   # used in group 21
}
#============== Other variables declared by tetracorder cmd files==============

NOGROUP0 = {3, 10, 11, 12,}
RRATIO_REFERENCES = {
            "[RATIOGREENVEG]": ("splib06b", 7260),
            "[RATIOGVEG1]": ("splib06b", 7644),
        }
# Tetracorder output directories, from cmds.start.t6.00a:104-149 (==[DIRgN] / ==[DIRcN])
OUTPUT_DIRS = {
    "group": {
        1: "group.1um",           2: "group.2um",            3: "group.veg",
        4: "group.1.5um-broad",   5: "group.2um-broad",      6: "group.2.5um",
        7: "group.3um",           8: "group.2.8um",          9: "group.zz",
        10: "group.3.5um_curve",  11: "group.3.5um",         12: "group.4um",
        13: "group.1.3-1.4um",    14: "group.1.4um",         15: "group.1.5um",
        16: "group.zz",           17: "group.1.7um",         18: "group.zz",
        19: "group.1.9um",        20: "group.ree",           21: "group.ree_g21",
        22: "group.ree_samar",    23: "group.3.8um",         24: "group.4um",
        25: "group.4.1um",        26: "group.4.25um-co2",    27: "group.4.5um",
        28: "group.5um",          29: "group.6um",           30: "group.7-8um",
        31: "group.8-10um",       32: "group.10-12um",       33: "group.12-14um",
        34: "group.14-16um",      35: "group.16-19um",       36: "group.19-22um",
        37: "group.ch4gas-2.3um", 38: "group.co2gas-2um",
    },
    "case": {
        1: "case.red-edge",  2: "case.veg.type",   3: "group.veg",
        4: "group.veg",      5: "group.veg",       6: "case.ep-cal-chl",
        7: "case.carbonate-2feat",
    },
}
#============== Preparing the rulesets for use by RuleEvaluator================

VARIABLE_PRESETS = {
    "default": DEFAULT_VALUES,
    "default_bak23": DEFAULT_BAK23_VALUES,
    "emit_c": EMIT_C_VALUES,
    "MMM_09c": MMM_09C_VALUES,
    "MMM255t": MMM255T_VALUES,
    "HYB2ryug": HYB2RYUG_VALUES,
}