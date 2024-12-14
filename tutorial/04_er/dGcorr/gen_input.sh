#!/bin/bash

exe=../../../packages/eror/tools/GENESIS/gen_input
genesis_log=../../02_genesis_md/host_in_water/genesis_log
dcd=../../02_genesis_md/host_in_water/run.dcd
dcd_gas=../../02_genesis_md/guest_in_gas/run.dcd

# 1. Prepare template input files
#    (parameters_er)

python2 $exe                \
  -l $genesis_log           \
  -x $dcd                   \
  -s $dcd_gas               \
  -d 1                      \
  --minenergy -60.0

# 2. Modify parameters_er to perform the test-particle insertion 
#    to the binding pocket

# ... Original parameters_er ...
#&ene_param
#        slttype = 3,
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
#        maxins = 1000,
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
        slttype = 3,
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
        maxins = 1000,
        engdiv = 1,
        use_anatra  = .true.,   ! ADDED
        calc_energy = .false.   ! ADDED
        
/

&output_param
  fhead = "out" ! header of several output files
/

&output_param
  fhead = "out"
/

&hist
      eclbin=5.0e-2, ecfbin=2.0e-3, ec0bin=2.0e-4, finfac=10.0e0,
      ecdmin=-60.000000, ecfmns=-0.20e0, ecdcen=0.0e0, eccore=20.0e0,
      ecdmax=1.0e11, pecore=200
/

!<--- ADDED
&input_param
  fxyz    = "../../03_reaction_coord/02_xyz/refbcd.xyz" ! used for fitting as reference coordinate
  fdxs    = "../../03_reaction_coord/03_sdf/sdf.dx"     ! DX-formatted sdf 
  fprmtop = "../../01_amber_setup/host-guest_in_water_for_insertion/hg_in_water_for_ins.prmtop" ! used for extracting potential parameters
/

&trajopt_param
  molinfo = "../../03_reaction_coord/01_spec/guest_ins.spec" ! 1 : ins
            "../../03_reaction_coord/01_spec/guest_fit.spec" ! 2 : fit (guest)
            "../../03_reaction_coord/01_spec/host_fit.spec"  ! 3 : fit (host)
            "../../03_reaction_coord/01_spec/guest_uij.spec" ! 4 : uattr (guest)  Dim # 1 
            "../../03_reaction_coord/01_spec/host_uij.spec"  ! 5 : uattr (host)   Dim # 1
/

&option_param
  ndim        = 1
  reac_coords = "ENERGY"
  rljcut      = 9.0
  use_sdf     = .true.
  uid_ins     = 1  ! solute molinfo id for insertion (CoM position) 
  uid_fit     = 2  ! solute molinfo id for fitting 
  vid_fit     = 3  ! guest  molinfo id for fitting
  box_shrink  = 20.0 20.0 20.0
/

&energy_param1
  uid         = 4 
  vid         = 5  
  calc_vdw    = .true.
  calc_elec   = .false.
  vdw         = "ATTRACTIVE"  ! Uattr is selected as reaction coordinate
  umin        = -38.8630197   ! Lower limit of Uattr
  umax        =  -9.3507904   ! Upper limit of Uattr
  ! Inserted configurations that satisfy umin <= uattr <= umax
  ! are used for calculating rho and chi 
  !
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
