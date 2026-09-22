# Fidelity test - database references, Python convolver (vendored tetracorder-lite/tetrapy/convolve)

## Run summary

| Item | Value |
|---|---:|
| Package run | 1018 s |
| Materials disabled | 26 |

### Disabled materials

| Material | Type | Number |
|---|---|---:|
| `CO2_ice g0.0300cm.curvc` | group | 10 |
| `CO2_ice g0.0300cm_linc` | group | 11 |
| `CO2_ice g0.1000cm.curvc` | group | 10 |
| `CO2_ice g0.1000cm_linc` | group | 11 |
| `co2_ice-on-moon62231.curvc` | group | 10 |
| `co2_ice-on-moon62231_linc` | group | 11 |
| `h2o-frost+ice-drops1-100k` | group | 0 |
| `h2o-frost+ice-drops3-100k` | group | 0 |
| `snow-ice-10mm+0.3water-273k` | group | 0 |
| `snow-ice-10mm-273k` | group | 0 |
| `snow-ice-1mm-273k` | group | 0 |
| `snow-ice-3mm+0.1water-273k` | group | 0 |
| `snow-ice-3mm+0.2water-273k` | group | 0 |
| `snow-ice-3mm+0.3water-273k` | group | 0 |
| `snow-ice-3mm-273k` | group | 0 |
| `snow.and.ice-170k` | group | 0 |
| `snow.and.ice-248k` | group | 0 |
| `snow.and.ice-77k` | group | 0 |
| `snow.melting.1a` | group | 0 |
| `snow.melting.3` | group | 0 |
| `snow.melting.8` | group | 0 |
| `snow.melting16+0.5veg` | group | 0 |
| `snow.melting1a+0.5veg` | group | 0 |
| `snow.melting9+0.5veg` | group | 0 |
| `snow.slush.16` | group | 0 |
| `snow.slush.9` | group | 0 |

## Group 1


**compared 140 materials (20 without native output); native px 573079, python px 573079, pooled jaccard 1.0000, common-pixel fit exact 0.9999**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `nanohematite.BR34b2b` | 131,357 | 131,359 | 1.000 | 1.000 | 1.000 | 1 |
| `Fe2++goeth+musc` | 85,571 | 85,571 | 1.000 | 1.000 | 1.000 | 1 |
| `generic.Fe2+nrw.cummingtonite` | 60,308 | 60,309 | 1.000 | 1.000 | 1.000 | 1 |
| `hematite-nano+goethite-fg-amix` | 46,030 | 46,030 | 1.000 | 1.000 | 1.000 | 1 |
| `chlor+muscphy` | 45,893 | 45,893 | 1.000 | 1.000 | 1.000 | 1 |
| `fe3+bearing1` | 40,238 | 40,236 | 1.000 | 1.000 | 1.000 | 1 |
| `jarosite-thinfilm-GDS243` | 30,869 | 30,867 | 1.000 | 1.000 | 1.000 | 1 |
| `goethite.thincoat` | 14,095 | 14,095 | 1.000 | 1.000 | 1.000 | 1 |
| `fe3+_hematite.4+goethite.1+qtz` | 13,413 | 13,413 | 1.000 | 1.000 | 1.000 | 0 |
| `jarosite-K-JR2501` | 13,213 | 13,214 | 1.000 | 1.000 | 1.000 | 1 |
| `nanohematite.BR34b2` | 12,724 | 12,724 | 1.000 | 1.000 | 1.000 | 1 |
| `hematite.fine.gr.fe2602` | 11,424 | 11,424 | 1.000 | 1.000 | 1.000 | 1 |
| `jarosite_br34a2` | 11,279 | 11,279 | 1.000 | 1.000 | 1.000 | 0 |
| `fe3+_goethite.nano0.05+quartz` | 10,366 | 10,366 | 1.000 | 1.000 | 1.000 | 0 |
| `chlor+goeth.propylzone` | 9,860 | 9,860 | 1.000 | 1.000 | 1.000 | 1 |
| `goethite+qtz.medgr.gds240` | 7,575 | 7,574 | 1.000 | 1.000 | 1.000 | 1 |
| `nh3jarosite-scr-nhj` | 5,087 | 5,087 | 1.000 | 1.000 | 1.000 | 0 |
| `fe3+_hematite.1+goethite.9-hg1` | 3,074 | 3,074 | 1.000 | 1.000 | 1.000 | 0 |
| `hematite.thincoat` | 2,942 | 2,942 | 1.000 | 1.000 | 1.000 | 0 |
| `fe3+_hematite.1+goethite.4+qtz` | 2,810 | 2,810 | 1.000 | 1.000 | 1.000 | 0 |
| `goethite.coarsegr.ws222` | 2,646 | 2,647 | 1.000 | 1.000 | 1.000 | 0 |
| `kjarosite200` | 2,338 | 2,339 | 1.000 | 1.000 | 1.000 | 0 |
| `Fe2+_hematite_weathering` | 1,269 | 1,269 | 1.000 | 1.000 | 1.000 | 0 |
| `goethite.medcoarsegr.mpc` | 1,238 | 1,237 | 0.999 | 1.000 | 1.000 | 0 |
| `fe3+_hematite.2+goethite.8-hg2` | 1,120 | 1,120 | 1.000 | 1.000 | 1.000 | 0 |
| `AMD.assemb1` | 1,112 | 1,112 | 1.000 | 1.000 | 1.000 | 0 |
| `generic.Fe2+med.jadeite` | 886 | 886 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] vegetation.dry+green` | 823 | 823 | 1.000 | 1.000 | 1.000 | 0 |
| `schwertmannite` | 618 | 618 | 1.000 | 1.000 | 1.000 | 0 |
| `fe2+br33a_bioqtzmonz_epidote` | 511 | 511 | 1.000 | 1.000 | 1.000 | 0 |
| `jarosite-na-gds100` | 408 | 408 | 1.000 | 1.000 | 1.000 | 0 |
| `Fe2+_and_hematite_br5b` | 313 | 313 | 1.000 | 1.000 | 1.000 | 0 |
| `fe3+_hematite.3+goethite.7-hg3` | 293 | 293 | 1.000 | 1.000 | 1.000 | 0 |
| `AMD.assemb2` | 186 | 186 | 1.000 | 1.000 | 1.000 | 0 |
| `generic.Fe2+.broad_br60b` | 174 | 174 | 1.000 | 1.000 | 1.000 | 0 |
| `goethite.fingr.mpcma2` | 168 | 168 | 1.000 | 1.000 | 1.000 | 0 |
| `hematite.lg.gr.br25c` | 133 | 133 | 1.000 | 1.000 | 1.000 | 0 |
| `hematite.fine.gr.gds76` | 116 | 116 | 1.000 | 1.000 | 1.000 | 0 |
| `goethite.medgr.ws222` | 76 | 76 | 1.000 | 1.000 | 1.000 | 0 |
| `fe3+_goethite.nano-70nm` | 68 | 68 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont100gpl+0.5veg` | 53 | 53 | 1.000 | 1.000 | 1.000 | 0 |
| `ferrihydrite` | 51 | 51 | 1.000 | 1.000 | 1.000 | 0 |
| `siderite1` | 47 | 47 | 1.000 | 1.000 | 1.000 | 0 |
| `hematite.med.gr.br25b` | 42 | 42 | 1.000 | 1.000 | 1.000 | 0 |
| `jarosite-na-gds101` | 30 | 30 | 1.000 | 1.000 | 1.000 | 0 |
| `generic.Fe2+nrw.hs-actinolite` | 28 | 28 | 1.000 | 1.000 | 1.000 | 0 |
| `hematite.lg.gr.br25a` | 26 | 26 | 1.000 | 1.000 | 1.000 | 0 |
| `goethite.lepidocrosite` | 25 | 25 | 1.000 | 1.000 | 1.000 | 0 |
| `albite_hs143` | 25 | 25 | 1.000 | 1.000 | 1.000 | 0 |
| `generic.Fe2+basalt_br46b` | 23 | 23 | 1.000 | 1.000 | 1.000 | 0 |
| `jarosite-H3O-SJ-1` | 22 | 22 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont300gpl` | 20 | 20 | 1.000 | 1.000 | 1.000 | 0 |
| `jarosite-na-nmnh95074-1` | 10 | 10 | 1.000 | 1.000 | 1.000 | 0 |
| `desert.varnish1` | 9 | 9 | 1.000 | 1.000 | 1.000 | 0 |
| `generic.Fe2+nrw.actinolite` | 7 | 7 | 1.000 | 1.000 | 1.000 | 0 |
| `generic.Fe2+.vbroad_br20` | 7 | 7 | 1.000 | 1.000 | nan | 0 |
| `hematite.lg.gr.br34c` | 5 | 5 | 1.000 | 1.000 | 1.000 | 0 |
| `pigeonite` | 5 | 5 | 1.000 | 1.000 | nan | 0 |
| `fe3+_hematite.nano.KP58610-400nm` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `nontronite.ng-1.a.g1` | 3 | 3 | 1.000 | 1.000 | nan | 0 |
| `feldspar.bytownite` | 3 | 3 | 1.000 | 1.000 | nan | 0 |
| `[g0] water+mont50gpl+0.5veg` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `chalcopyrite` | 2 | 2 | 1.000 | 1.000 | nan | 0 |
| `hematite.med.gr.gds27` | 1 | 1 | 1.000 | 1.000 | 1.000 | 0 |
| `jarosite-k-90c-gds98` | 1 | 1 | 1.000 | 1.000 | 1.000 | 0 |
| `pyroxene.hypersthene.pyx02.e.g1` | 1 | 1 | 1.000 | 1.000 | nan | 0 |
| `fe2+pyroxene.diopside` | 1 | 1 | 1.000 | 1.000 | nan | 0 |
| `greenslime` | 1 | 1 | 1.000 | 1.000 | nan | 0 |
| `[g0] white_pvc_pipe` | 1 | 1 | 1.000 | 1.000 | 1.000 | 0 |

## Group 2


**compared 227 materials (21 without native output); native px 596794, python px 596794, pooled jaccard 1.0000, common-pixel fit exact 1.0000**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `muscovite-medhigh-Al` | 89,694 | 89,693 | 1.000 | 1.000 | 1.000 | 1 |
| `calcite+0.3muscovite` | 66,854 | 66,854 | 1.000 | 1.000 | 1.000 | 0 |
| `beidellite_gds123` | 65,271 | 65,271 | 1.000 | 1.000 | 1.000 | 1 |
| `calcite.ws272.g2` | 55,267 | 55,267 | 1.000 | 1.000 | 1.000 | 1 |
| `kaolgrp_nacrite` | 34,914 | 34,914 | 1.000 | 1.000 | 1.000 | 0 |
| `montna` | 33,649 | 33,649 | 1.000 | 1.000 | 1.000 | 1 |
| `muscovite-med-Al` | 24,865 | 24,866 | 1.000 | 1.000 | 1.000 | 1 |
| `dolomite` | 21,333 | 21,333 | 1.000 | 1.000 | 1.000 | 0 |
| `kaolin.5+muscov.medhighAl` | 19,948 | 19,948 | 1.000 | 1.000 | 1.000 | 1 |
| `kaolin+musc.intimat` | 19,331 | 19,331 | 1.000 | 1.000 | 1.000 | 1 |
| `aragonite` | 13,168 | 13,168 | 1.000 | 1.000 | 1.000 | 1 |
| `calcite+0.2Na-mont` | 11,025 | 11,025 | 1.000 | 1.000 | 1.000 | 0 |
| `calcite.33+kaol.33+mus-AMXr2` | 10,400 | 10,400 | 1.000 | 1.000 | 1.000 | 0 |
| `chalcedony` | 10,184 | 10,184 | 1.000 | 1.000 | 1.000 | 0 |
| `alunite.5+kaol.5` | 8,697 | 8,698 | 1.000 | 1.000 | 1.000 | 0 |
| `paragonite` | 8,115 | 8,115 | 1.000 | 1.000 | 1.000 | 0 |
| `kaolpxl` | 6,945 | 6,945 | 1.000 | 1.000 | 1.000 | 0 |
| `halloy` | 6,595 | 6,595 | 1.000 | 1.000 | 1.000 | 0 |
| `dolomite.5+Na-mont.5` | 6,483 | 6,483 | 1.000 | 1.000 | 1.000 | 0 |
| `plastic.hdpe.1` | 5,995 | 5,995 | 1.000 | 1.000 | 1.000 | 0 |
| `Kalun+kaol.intmx` | 5,906 | 5,906 | 1.000 | 1.000 | 1.000 | 0 |
| `muscovite+chlorite` | 5,802 | 5,802 | 1.000 | 0.999 | 1.000 | 1 |
| `kaol.75+alun.25` | 5,582 | 5,581 | 1.000 | 1.000 | 1.000 | 0 |
| `beidellite_gds124` | 4,071 | 4,071 | 1.000 | 1.000 | 1.000 | 0 |
| `sulfate-mix_gypsum.trace.dust+debris-WTC01` | 4,071 | 4,071 | 1.000 | 1.000 | 1.000 | 1 |
| `kaolin.5+smect.5` | 3,836 | 3,836 | 1.000 | 1.000 | 1.000 | 0 |
| `hectorite` | 3,638 | 3,638 | 1.000 | 1.000 | 1.000 | 1 |
| `alunNa03` | 2,834 | 2,834 | 1.000 | 1.000 | 1.000 | 0 |
| `calcite.25+dolom.25+Na-mont.5` | 2,520 | 2,520 | 1.000 | 1.000 | 1.000 | 0 |
| `kaolwxl` | 2,442 | 2,442 | 1.000 | 1.000 | 1.000 | 0 |
| `drygrass+.17Na-mont` | 2,241 | 2,241 | 1.000 | 1.000 | 1.000 | 0 |
| `plastic.Vinyl` | 2,086 | 2,086 | 1.000 | 1.000 | 1.000 | 0 |
| `alunite.33+kaol.33+musc.33` | 1,882 | 1,882 | 1.000 | 1.000 | 1.000 | 0 |
| `kaolin.5+muscov.medAl` | 1,788 | 1,788 | 1.000 | 1.000 | 1.000 | 0 |
| `jarosite-lowT` | 1,749 | 1,749 | 1.000 | 1.000 | 1.000 | 0 |
| `Na-alun+kaol.intmx` | 1,645 | 1,645 | 1.000 | 1.000 | 1.000 | 0 |
| `cookeite-CAr-1.c` | 1,642 | 1,642 | 1.000 | 1.000 | 1.000 | 0 |
| `chrysotile.gypsum.wtc01-8` | 1,588 | 1,588 | 1.000 | 1.000 | 1.000 | 0 |
| `muscovite-low-Al` | 1,436 | 1,436 | 1.000 | 1.000 | 1.000 | 0 |
| `hydrated_basaltic_glass` | 1,364 | 1,364 | 1.000 | 1.000 | 1.000 | 0 |
| `phlogopite` | 1,305 | 1,305 | 1.000 | 1.000 | 1.000 | 0 |
| `plastic.tarp1` | 1,235 | 1,235 | 1.000 | 1.000 | 1.000 | 0 |
| `sepiolite` | 1,100 | 1,100 | 1.000 | 0.999 | 1.000 | 1 |
| `kaolin.3+smect.7` | 1,024 | 1,024 | 1.000 | 1.000 | 1.000 | 0 |
| `margarite` | 870 | 870 | 1.000 | 1.000 | 1.000 | 0 |
| `methane-gas-2.37um-experimental` | 852 | 852 | 1.000 | 1.000 | 1.000 | 0 |
| `nontronite.ng-1.a.g2` | 833 | 833 | 1.000 | 1.000 | 1.000 | 0 |
| `calcite+0.5Ca-mont` | 825 | 825 | 1.000 | 1.000 | 1.000 | 0 |
| `montca` | 820 | 820 | 1.000 | 1.000 | 1.000 | 0 |
| `dolo+.5ca-mont` | 799 | 799 | 1.000 | 1.000 | 1.000 | 0 |
| `vermiculite_GDS458` | 799 | 799 | 1.000 | 1.000 | 1.000 | 0 |
| `montmorillonite_fe` | 775 | 775 | 1.000 | 1.000 | 1.000 | 0 |
| `palygorskite` | 708 | 708 | 1.000 | 1.000 | 1.000 | 0 |
| `sulfate-bloedite` | 650 | 650 | 1.000 | 1.000 | 0.992 | 0 |
| `alun66K34Na.low` | 613 | 613 | 1.000 | 1.000 | 1.000 | 0 |
| `alun73K27Na.low` | 533 | 533 | 1.000 | 1.000 | 1.000 | 0 |
| `magnesite` | 533 | 533 | 1.000 | 1.000 | 1.000 | 0 |
| `vermiculite_ALB4SA00` | 508 | 508 | 1.000 | 1.000 | 1.000 | 0 |
| `saponite.or.talc` | 477 | 477 | 1.000 | 1.000 | 1.000 | 0 |
| `chlorite-skarn` | 428 | 428 | 1.000 | 1.000 | 1.000 | 0 |
| `gypsum.trace.dust+debris-WTC01-28` | 344 | 344 | 1.000 | 1.000 | 1.000 | 0 |
| `calcite+dolomite.5` | 297 | 297 | 1.000 | 1.000 | 1.000 | 0 |
| `niter` | 282 | 282 | 1.000 | 1.000 | 1.000 | 0 |
| `dickite` | 263 | 263 | 1.000 | 1.000 | 1.000 | 0 |
| `kalun250c` | 257 | 257 | 1.000 | 1.000 | 1.000 | 0 |
| `clinochlore.fe.sc-cca-1` | 225 | 225 | 1.000 | 1.000 | 1.000 | 0 |
| `calcite+0.2Ca-mont` | 186 | 186 | 1.000 | 1.000 | 1.000 | 0 |
| `portlandite` | 182 | 182 | 1.000 | 1.000 | 1.000 | 0 |
| `calcite+0.2kaolwxl` | 155 | 155 | 1.000 | 1.000 | 1.000 | 0 |
| `cookeite-CAr-1.a` | 148 | 148 | 1.000 | 1.000 | 1.000 | 0 |
| `dryveg.grass.golden` | 129 | 129 | 1.000 | 1.000 | 1.000 | 0 |
| `calcite0.7+kaol0.3` | 115 | 115 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] vegetation.dry+green` | 115 | 115 | 1.000 | 1.000 | 1.000 | 0 |
| `chrysotile.med-fine-fibrous` | 114 | 114 | 1.000 | 1.000 | 1.000 | 0 |
| `vermiculite_WS681` | 111 | 111 | 1.000 | 1.000 | 1.000 | 0 |
| `kalun150c` | 101 | 101 | 1.000 | 1.000 | 1.000 | 0 |
| `hornblende` | 98 | 98 | 1.000 | 1.000 | 1.000 | 0 |
| `prehnite+.75chlorite` | 97 | 97 | 1.000 | 1.000 | 1.000 | 0 |
| `prehnite.gds613` | 97 | 97 | 1.000 | 1.000 | 1.000 | 0 |
| `alunite.5+musc.5` | 86 | 86 | 1.000 | 1.000 | 1.000 | 0 |
| `vermiculite_WS682` | 74 | 74 | 1.000 | 1.000 | 1.000 | 0 |
| `buddington.namont2` | 72 | 72 | 1.000 | 1.000 | 1.000 | 0 |
| `clinochlore.fe.gds157` | 70 | 70 | 1.000 | 1.000 | 1.000 | 0 |
| `clinochlore.nmnh83369` | 57 | 57 | 1.000 | 1.000 | 1.000 | 0 |
| `amphibole` | 56 | 56 | 1.000 | 1.000 | 1.000 | 0 |
| `na40alun400c` | 48 | 48 | 1.000 | 1.000 | 1.000 | 0 |
| `dryveg.grass.long` | 44 | 44 | 1.000 | 1.000 | 1.000 | 0 |
| `strontianite` | 43 | 43 | 1.000 | 1.000 | 1.000 | 0 |
| `natroalun+dickite` | 30 | 30 | 1.000 | 1.000 | 1.000 | 0 |
| `clintonite` | 29 | 29 | 1.000 | 1.000 | 1.000 | 0 |
| `unleaded.gas` | 28 | 28 | 1.000 | 1.000 | 1.000 | 0 |
| `jarosite-K` | 23 | 23 | 1.000 | 1.000 | 1.000 | 0 |
| `tremolite.or.talc` | 23 | 23 | 1.000 | 1.000 | 1.000 | 0 |
| `prehnite+.50chlorite` | 22 | 22 | 1.000 | 1.000 | 1.000 | 0 |
| `elbaite` | 21 | 21 | 1.000 | 1.000 | 1.000 | 0 |
| `prehnite+.67chlorite` | 20 | 20 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont100gpl+0.5veg` | 19 | 19 | 1.000 | 1.000 | 1.000 | 0 |
| `sulfate_kieserite` | 15 | 15 | 1.000 | 1.000 | 1.000 | 0 |
| `chrysotile.fine-fibrous` | 15 | 15 | 1.000 | 1.000 | 1.000 | 0 |
| `biotite-phlogopite-P23-5` | 14 | 14 | 1.000 | 1.000 | nan | 0 |
| `musc+pyroph` | 13 | 13 | 1.000 | 1.000 | 1.000 | 0 |
| `muscoviteFerich` | 12 | 12 | 1.000 | 1.000 | 1.000 | 0 |
| `sulfate_szomolnokite` | 11 | 11 | 1.000 | 1.000 | 1.000 | 0 |
| `talc.fngrnd` | 10 | 10 | 1.000 | 1.000 | 1.000 | 0 |
| `illite` | 9 | 9 | 1.000 | 1.000 | 1.000 | 0 |
| `illite.gds4` | 8 | 8 | 1.000 | 1.000 | 1.000 | 0 |
| `calcite.25+dolom.25+Ca-mont.5` | 7 | 7 | 1.000 | 1.000 | 1.000 | 0 |
| `siderite` | 7 | 7 | 1.000 | 1.000 | 1.000 | 0 |
| `calcite.33+Ca-mont.67` | 6 | 6 | 1.000 | 1.000 | 1.000 | 0 |
| `pyroph.5+alunit.5` | 5 | 5 | 1.000 | 1.000 | 1.000 | 0 |
| `talc+calcite.parkcity` | 5 | 5 | 1.000 | 1.000 | 1.000 | 0 |
| `lepidolite` | 5 | 5 | 1.000 | 1.000 | 1.000 | 0 |
| `musc+jarosite.intimat` | 4 | 4 | 1.000 | 1.000 | 1.000 | 0 |
| `naalun150c` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont50gpl+0.5veg` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `kalun450c` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `na82alun100c` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `pyrophyllite` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `buddington` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `cyanide-CdK` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont300gpl` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `actinolite` | 1 | 1 | 1.000 | 1.000 | 1.000 | 0 |
| `jarosite-Na` | 1 | 1 | 1.000 | 1.000 | nan | 0 |
| `gibbsite` | 1 | 1 | 1.000 | 1.000 | 1.000 | 0 |
| `prehnite+muscovite` | 1 | 1 | 1.000 | 1.000 | nan | 0 |
| `cronstedtite_M3542` | 1 | 1 | 1.000 | 1.000 | 1.000 | 0 |
| `sodium_nitrate` | 1 | 1 | 1.000 | 1.000 | nan | 0 |

## Group 3


**compared 4 materials (0 without native output); native px 550506, python px 550506, pooled jaccard 1.0000, common-pixel fit exact 0.9999**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `vegetation.weak.map` | 548,157 | 548,157 | 1.000 | 1.000 | 1.000 | 1 |
| `vegetation.dry_nonphotosyn` | 1,914 | 1,914 | 1.000 | 1.000 | 1.000 | 0 |
| `vegetation.map` | 427 | 427 | 1.000 | 1.000 | 1.000 | 0 |
| `vegetation.dry+green.g3` | 8 | 8 | 1.000 | 1.000 | 1.000 | 0 |

## Group 4


**compared 48 materials (21 without native output); native px 401157, python px 401157, pooled jaccard 1.0000, common-pixel fit exact 0.9999**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `g4-generic.Fe2+nrw.cummingtonite` | 301,030 | 301,030 | 1.000 | 1.000 | 1.000 | 1 |
| `g4-fe2+generic_chlor+muscphy` | 96,629 | 96,629 | 1.000 | 1.000 | 1.000 | 1 |
| `[g0] vegetation.dry+green` | 3,304 | 3,304 | 1.000 | 1.000 | 1.000 | 0 |
| `microcline_hs103.3b` | 97 | 97 | 1.000 | 0.990 | nan | 1 |
| `[g0] water+mont100gpl+0.5veg` | 70 | 70 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont300gpl` | 23 | 23 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont50gpl+0.5veg` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] white_pvc_pipe` | 1 | 1 | 1.000 | 1.000 | 1.000 | 0 |

## Group 5


**compared 35 materials (20 without native output); native px 7265, python px 7265, pooled jaccard 1.0000, common-pixel fit exact 1.0000**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `[g0] vegetation.dry+green` | 4,087 | 4,087 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont100gpl+0.5veg` | 2,052 | 2,052 | 1.000 | 1.000 | 1.000 | 0 |
| `pyroxene.enstatite` | 761 | 761 | 1.000 | 1.000 | nan | 0 |
| `pyroxene.pigeonite` | 197 | 197 | 1.000 | 1.000 | nan | 0 |
| `pyroxene.hypersthene.pyx02.e.g5` | 99 | 99 | 1.000 | 1.000 | nan | 0 |
| `[g0] water+mont300gpl` | 34 | 34 | 1.000 | 1.000 | 1.000 | 0 |
| `pyroxene.augite` | 29 | 29 | 1.000 | 1.000 | nan | 0 |
| `[g0] water+mont50gpl+0.5veg` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] white_pvc_pipe` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `pyroxene.bronzite.hs9.3b.g5` | 1 | 1 | 1.000 | 1.000 | nan | 0 |

## Group 20


**compared 41 materials (20 without native output); native px 17727, python px 17727, pooled jaccard 1.0000, common-pixel fit exact 0.9987**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `ree_nd_xenotime+monazite` | 9,728 | 9,728 | 1.000 | 0.998 | 1.000 | 1 |
| `[g0] vegetation.dry+green` | 3,987 | 3,987 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont100gpl+0.5veg` | 2,033 | 2,033 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_xenotime` | 1,701 | 1,701 | 1.000 | 0.998 | 0.998 | 1 |
| `ree_nd_bastnaesite` | 168 | 168 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont300gpl` | 34 | 34 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_pyromorphite` | 30 | 30 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_monazite_gds957` | 22 | 22 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_ancylite` | 15 | 15 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont50gpl+0.5veg` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_neodymium_oxide` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_neodymium_oxide-0.6um` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] white_pvc_pipe` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |

## Group 21


**compared 40 materials (20 without native output); native px 17726, python px 17726, pooled jaccard 1.0000, common-pixel fit exact 0.9987**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `ree_nd_xenotime+monazite_g21` | 9,728 | 9,728 | 1.000 | 0.998 | 0.998 | 1 |
| `[g0] vegetation.dry+green` | 3,987 | 3,987 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont100gpl+0.5veg` | 2,034 | 2,034 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_xenotime_g21` | 1,701 | 1,701 | 1.000 | 0.998 | 0.996 | 1 |
| `ree_nd_bastnaesite_g21` | 168 | 168 | 1.000 | 1.000 | 0.992 | 0 |
| `[g0] water+mont300gpl` | 34 | 34 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_pyromorphite_g21` | 30 | 30 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_monazite_gds957_g21` | 22 | 22 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_ancylite_g21` | 15 | 15 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont50gpl+0.5veg` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `ree_nd_neodymium_oxide_g21` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] white_pvc_pipe` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |

## Group 22


**compared 29 materials (20 without native output); native px 6178, python px 6178, pooled jaccard 1.0000, common-pixel fit exact 1.0000**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `[g0] vegetation.dry+green` | 4,087 | 4,087 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont100gpl+0.5veg` | 2,052 | 2,052 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont300gpl` | 34 | 34 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont50gpl+0.5veg` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] white_pvc_pipe` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |

## Group 37


**compared 29 materials (20 without native output); native px 6178, python px 6178, pooled jaccard 1.0000, common-pixel fit exact 1.0000**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `[g0] vegetation.dry+green` | 4,087 | 4,087 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont100gpl+0.5veg` | 2,052 | 2,052 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont300gpl` | 34 | 34 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont50gpl+0.5veg` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] white_pvc_pipe` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |

## Group 38


**compared 29 materials (20 without native output); native px 107050, python px 107050, pooled jaccard 1.0000, common-pixel fit exact 0.9997**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `co2-gas-2um_g38-experimental` | 101,269 | 101,269 | 1.000 | 1.000 | 1.000 | 1 |
| `[g0] vegetation.dry+green` | 3,732 | 3,732 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont100gpl+0.5veg` | 2,019 | 2,019 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont300gpl` | 25 | 25 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] water+mont50gpl+0.5veg` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `[g0] white_pvc_pipe` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |

## Case 1


**compared 2 materials (0 without native output); native px 528969, python px 528969, pooled jaccard 1.0000, common-pixel fit exact 0.9996**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `red.edge.shift.2` | 466,101 | 466,097 | 1.000 | 1.000 | 1.000 | 1 |
| `red.edge.shift.1` | 62,868 | 62,872 | 1.000 | 0.998 | 1.000 | 1 |

## Case 2


**compared 15 materials (0 without native output); native px 466554, python px 466554, pooled jaccard 1.0000, common-pixel fit exact 1.0000**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `vegetation.type13` | 465,846 | 465,845 | 1.000 | 1.000 | 1.000 | 1 |
| `vegetation.type6` | 652 | 653 | 0.998 | 1.000 | 1.000 | 0 |
| `vegetation.type11` | 21 | 21 | 1.000 | 1.000 | 1.000 | 0 |
| `vegetation.type4` | 19 | 19 | 1.000 | 1.000 | 1.000 | 0 |
| `vegetation.type15` | 6 | 6 | 1.000 | 1.000 | 1.000 | 0 |
| `vegetation.type7` | 4 | 4 | 1.000 | 1.000 | 1.000 | 0 |
| `vegetation.type14` | 3 | 3 | 1.000 | 1.000 | 1.000 | 0 |
| `vegetation.type9` | 2 | 2 | 1.000 | 1.000 | 1.000 | 0 |
| `vegetation.type3` | 1 | 1 | 1.000 | 1.000 | 1.000 | 0 |

## Case 3


**compared 1 materials (0 without native output); native px 432034, python px 432034, pooled jaccard 1.0000, common-pixel fit exact 0.9999**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `veg0.9um.band` | 432,034 | 432,034 | 1.000 | 1.000 | 1.000 | 1 |

## Case 4


**compared 1 materials (0 without native output); native px 460509, python px 460509, pooled jaccard 1.0000, common-pixel fit exact 0.9999**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `veg1.2um.band` | 460,509 | 460,509 | 1.000 | 1.000 | 1.000 | 1 |

## Case 5


**compared 1 materials (0 without native output); native px 514984, python px 514984, pooled jaccard 1.0000, common-pixel fit exact 1.0000**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `veg1.4um.band` | 514,984 | 514,984 | 1.000 | 1.000 | 1.000 | 1 |

## Case 6


**compared 13 materials (0 without native output); native px 45759, python px 45759, pooled jaccard 1.0000, common-pixel fit exact 1.0000**

| Material | Native | Python | Jaccard | Fit exact | Depth exact | Fit max |
|---|---:|---:|---:|---:|---:|---:|
| `calcite.ws272.c6` | 45,686 | 45,686 | 1.000 | 1.000 | 1.000 | 1 |
| `epidote.gds26.a.c6` | 49 | 49 | 1.000 | 1.000 | 1.000 | 0 |
| `CEC-D.33Ep66Ca` | 20 | 20 | 1.000 | 1.000 | 1.000 | 0 |
| `CEC-J.33Ch66Ca` | 4 | 4 | 1.000 | 1.000 | 1.000 | 0 |

## Summary

| Set | Materials | No native output | Native pixels | Python pixels | Jaccard | Fit exact |
|---|---:|---:|---:|---:|---:|---:|
| Group 1 | 140 | 20 | 573,079 | 573,079 | 1.0000 | 0.9999 |
| Group 2 | 227 | 21 | 596,794 | 596,794 | 1.0000 | 1.0000 |
| Group 3 | 4 | 0 | 550,506 | 550,506 | 1.0000 | 0.9999 |
| Group 4 | 48 | 21 | 401,157 | 401,157 | 1.0000 | 0.9999 |
| Group 5 | 35 | 20 | 7,265 | 7,265 | 1.0000 | 1.0000 |
| Group 20 | 41 | 20 | 17,727 | 17,727 | 1.0000 | 0.9987 |
| Group 21 | 40 | 20 | 17,726 | 17,726 | 1.0000 | 0.9987 |
| Group 22 | 29 | 20 | 6,178 | 6,178 | 1.0000 | 1.0000 |
| Group 37 | 29 | 20 | 6,178 | 6,178 | 1.0000 | 1.0000 |
| Group 38 | 29 | 20 | 107,050 | 107,050 | 1.0000 | 0.9997 |
| Case 1 | 2 | 0 | 528,969 | 528,969 | 1.0000 | 0.9996 |
| Case 2 | 15 | 0 | 466,554 | 466,554 | 1.0000 | 1.0000 |
| Case 3 | 1 | 0 | 432,034 | 432,034 | 1.0000 | 0.9999 |
| Case 4 | 1 | 0 | 460,509 | 460,509 | 1.0000 | 0.9999 |
| Case 5 | 1 | 0 | 514,984 | 514,984 | 1.0000 | 1.0000 |
| Case 6 | 13 | 0 | 45,759 | 45,759 | 1.0000 | 1.0000 |
