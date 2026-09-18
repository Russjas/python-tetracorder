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
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
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
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
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
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
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
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
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
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.0007, 0.0011],    # r*bd for 2.5-    micron region
"[RBD30]": [0.0005, 0.001],     # r*bd for 3-      micron region
"[RBD35]": [0.0005, 0.001],     # r*bd for 3.5-    micron region
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
##unsure which is preserved by specpr behaviour
"[RBD25]": [0.00007, 0.00011],    # r*bd for 2.5-    micron region #sic Tetracorder
"[RBD30]": [0.00005, 0.0001],     # r*bd for 3-      micron region
"[RBD35]": [0.00005, 0.0001],     # r*bd for 3.5-    micron region
}
#============== Other variables declared by tetracorder cmd files==============

NOGROUP0 = {3, 10, 11, 12,}
RRATIO_REFERENCES = {
            "[RATIOGREENVEG]": ("splib06b", 7260),
            "[RATIOGVEG1]": ("splib06b", 7644),
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