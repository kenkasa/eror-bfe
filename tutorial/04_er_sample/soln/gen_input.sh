#!/bin/bash

exe=../../../packages/eror/tools/GENESIS/gen_input
genesis_log=../../02_genesis_md/host-guest_in_water/genesis_log
dcd=../../02_genesis_md/host-guest_in_water/run.dcd

# 1. Prepare template input files
#    (parameters_er)

python3 $exe                \
  -l $genesis_log           \
  -x $dcd                   \
  -d 1                      \
  --minenergy -60.0

# 2. Modify parameters_er to perform the analysis at the OR state 

# ... Original parameters_er ...
#&ene_param
#        slttype = 1,
#        sltspec = 2,
#        boxshp = 1,
#        estype = 1,
#        inptemp = 298,
#        ljformat = 5,
#        ljswitch = 0,
#        upljcut = 9,
#        lwljcut = 9,
#        cltype = 2,
#        elecut = 9,
#        screen = 0.30768,
#        splodr = 4,
#        ms1max = 64,
#        ms2max = 48,
#        ms3max = 48,
#        engdiv = 1,
#/
#&hist
#      eclbin=5.0e-2, ecfbin=2.0e-3, ec0bin=2.0e-4, finfac=10.0e0,
#      ecdmin=-60.000000, ecfmns=-0.20e0, ecdcen=0.0e0, eccore=20.0e0,
#      ecdmax=1.0e11, pecore=200
#/
# ..............................

cp parameters_er parameters_er.orig

cat << EOF > parameters_er
&ene_param
        slttype = 1,
        sltspec = 2,
        boxshp = 1,
        estype = 1,
        inptemp = 298,
        ljformat = 5,
        ljswitch = 0,
        upljcut = 9,
        lwljcut = 9,
        cltype = 2,
        elecut = 9,
        screen = 0.30768,
        splodr = 4,
        ms1max = 64,
        ms2max = 48,
        ms3max = 48,
        engdiv = 1,
        use_eror = .true. ! ADDED
/
&hist
      eclbin=5.0e-2, ecfbin=2.0e-3, ec0bin=2.0e-4, finfac=10.0e0,
      ecdmin=-60.000000, ecfmns=-0.20e0, ecdcen=0.0e0, eccore=20.0e0,
      ecdmax=1.0e11, pecore=200
/
! <--- ADDED
&eror
      refs_filename = '../refs/engref'
/
! --->
EOF

exit
Usage: gen_input [options]

Options:
  -h, --help            show this help message and exit
  -t TARGETDIRECTORY, --targetdirectory=TARGETDIRECTORY
                        Directory that contains configuration data
  -l LOG, --log=LOG     GENESIS log file
  -x DCD, --dcd=DCD     Solution or reference trajectory file
  -p PDB, --pdb=PDB     (deprecated) same as --rigid
  -r PDB, --rigid=PDB   Solute structure PDB file (for rigid insertion)
  -s SOLUTE, --flexible=SOLUTE
                        Solute trajectory file (for flexible insertion)
  -d DIV, --div=DIV     Number of divided sections
  --minenergy=MINENERGY
                        Minimum energy in kcal/mol
  --maxins=MAXINS       Number of insertions for refs (default: 1000)
