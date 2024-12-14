# Tutorial (Binding of aspirin to b-cyclodextrin in water)] 

In this tutorial, "solvation" free energy of aspirin (solute, resname APR) at cavity of b-cyclodextrin (host, resname BCD) (delta mu^B) is computed using ER-OR.
The computation of solvation free energy at the dissociate state (delta mu^D) is typical appication for the original ERmod program, 
and thus we omit the scheme of computing delta mu^D.    

## Protocols

Directory tree:
```
├── 01_amber_setup                           # Amber parameter files (prmtop, inpcrd) are stored
│   ├── guest_in_gas                         #   parameter files for guest in gas
│   │   ├── g.inpcrd
│   │   ├── g.pdb
│   │   └── g.prmtop
│   ├── host-guest_in_water                  #   parameter files for host-guest in water (host, guest, water)
│   │   ├── hg_in_water.inpcrd
│   │   ├── hg_in_water.pdb
│   │   └── hg_in_water.prmtop
│   ├── host-guest_in_water_for_insertion    #   parameter files for host-guest in water used for insertion (host, water, guest)
│   │   ├── hg_in_water_for_ins.inpcrd
│   │   ├── hg_in_water_for_ins.pdb
│   │   └── hg_in_water_for_ins.prmtop
│   └── host_in_water                        #   parameter files for host in water (host, water) 
│       ├── h_in_water.inpcrd
│       ├── h_in_water.pdb
│       └── h_in_water.prmtop
├── 02_genesis_md                            #   Trajectory and GENESIS log files are stored
│   ├── guest_in_gas
│   │   └── run.dcd
│   ├── host-guest_in_water
│   │   ├── genesis_log
│   │   └── run.dcd
│   └── host_in_water
│       ├── genesis_log
│       └── run.dcd
├── 03_reaction_coord                        # Reaction coordinates for defining bound state
│   ├── 01_spec                              #   VMD atomselection data for host and guest are prepared 
│   │   ├── guest_fit.spec
│   │   ├── guest_ins.spec
│   │   ├── guest_uij.spec
│   │   ├── host_fit.spec
│   │   ├── host_uij.spec
│   │   └── run.sh
│   ├── 02_xyz                               #   xyz-formatted host coordinate file is prepared  
│   │   ├── refbcd.pdb
│   │   ├── refbcd.xyz
│   │   └── run.sh
│   └── 03_sdf                               #   Spatial distribution function data is stored 
│       └── sdf.dx
├── 04_er                                    # ER-OR calculation is performed
│   ├── dGcorr                               #   Standard-state correction term is calculated 
│   │   ├── Long_Sample
│   │   │   └── out.popul
│   │   ├── gen_input.sh
│   │   └── run.sh
│   ├── prep_dir.sh
│   ├── refs                                 #   Rho and chi in the reference are calculated 
│   │   ├── Long_Sample
│   │   │   ├── corref.01.zip
│   │   │   ├── engref.01
│   │   │   └── weight_refs
│   │   ├── gen_input.sh
│   │   └── run.sh
│   ├── run.sh
│   └── soln                                 #   Rho in the solution is calculated
│       ├── Long_Sample
│       ├── gen_input.sh
│       └── run.sh
├── 04_er_sample                             # Reference data
```

In ``01_amber_setup``, the Amber parameter files are stored, and there is nothing to do in this directory.

In ``02_genesis_md``, the trajectory and GENESIS log files are sto
red. These data are neccessary for the ER-OR calculation. Due to the file size limitation on GitHub, the trajecotries on the short timescale are provided.

In ``03_reaction_coord``, the parameter files for calculating the reaction coordinate are prepared.

In ``04_er``, the ER-OR calculation is performed.

All the output files in ``04_er_sample`` can be reproduced in ``04_er`` after finishing all the procedures. 


## Reaction coordinate (03_reaction_coord)

In the case of CD-aspirin, the bound state is defined using attractive part of LJ interaction (uattr) between the CD and aspirin, and the spatial distribution function (SDF) of aspirin, g(r).   
Bound state: umin <= uattr <= umax, g(r) > 0

To calculate the uattr and the position of aspirin, the parameter files are prepared in this directory.

```
cd 03_reaction_coord/01_spec
ls 
> run.sh sample

cat run.sh
> #!/bin/bash
> 
> prmtop=../../01_amber_setup/host-guest_in_water_for_insertion/hg_in_water_for_ins.prmtop
> inpcrd=../../01_amber_setup/host-guest_in_water_for_insertion/hg_in_water_for_ins.inpcrd
> 
> anatra specfile                                              \
>   -stype   parm7                                             \
>   -sfile   $prmtop                                           \
>   -tintype rst7                                              \
>   -tin     $inpcrd                                           \
>   -fhead   out                                               \
>   -sel0    resname APR                                       \
>   -sel1    resname BCD                                       \
>   -sel2    resname APR and noh                               \
>   -sel3    resname BCD and name O15 O16 O17 O18 O19 O20 O21
> 
> 
> #Uattr  
> cp out_0.spec guest_uij.spec
> cp out_1.spec host_uij.spec
> 
> #Fitting  
> cp out_2.spec guest_fit.spec
> cp out_3.spec host_fit.spec
> 
> #Ins  
> cp out_2.spec guest_ins.spec
> 
> rm -f out_*.spec

./run.sh
ls
> guest_fit.spec  guest_uij.spec  host_uij.spec  sample
> guest_ins.spec  host_fit.spec   run.sh 
```
By executing run.sh, the sets of atoms required for calculating the reaction coordinates are listed in *.spec files.

```
cd ../02_xyz
ls
> refbcd.pdb   # reference CD structure 
> run.sh  
> sample

cat run.sh
> #!/bin/bash
> 
> anatra trjconv                                 \
>   -stype      pdb                              \
>   -sfile      refbcd.pdb                       \
>   -tintype    pdb                              \
>   -tin        refbcd.pdb                       \
>   -totype     xyz                              \
>   -to         refbcd.xyz                       \
>   -sel0       name O15 O16 O17 O18 O19 O20 O21 \
>   -outselid   0

./run.sh
ls
> refbcd.pdb  
> run.sh  
> sample
> refbcd.xyz # Reference coordinate of CD  
```
To use the SDF as the bound-state criteria, the translation/rotateion of CD to superimpose the reference structure of CD is required. At this step, the atoms in CD used for fitting are selected and stored in ``refbcd.xyz``. 
```
cd ../03_sdf
ls
> sdf.dx
``` 
In this directory, SDF data (``sdf.dx``) is present. 
This can be computed using the solution trajectories at the bound state. Since the sufficiently long trajectory can not be placed in the repository, only the SDF data is placed.
The program for computinng the SDF is implemented in ANATRA (``anatra sd``). Before using this program, the fitting of the trajecotry to the reference should be done using ``anatra trjconv``.

## ER-OR calculation (04_er)

### Setup directories

```
cd ../../04_er
ls 
> prep_dir.sh
> run.sh
> refs
> dGcorr    
> soln

cat prep_dir.sh
> #!/usr/bin/env bash
> 
> exe=../../packages/eror/tools/GENESIS/gen_structure_AMBER
> prmtop=../01_amber_setup/host-guest_in_water/hg_in_water.prmtop
> 
> # 1. Prepare soln and refs directories 
> #   soln : working directory for rho in the solution  system
> #   refs : working directory for rho and chi in the reference system
> python2 $exe  \
>   -t $prmtop  \
>   -s APR
> 
> # 2. Copy refs directories for dGcorr calc. 
> #
> rsync -av --exclude "run.sh" --exclude "gen_input.sh" refs/ dGcorr/

./prep_dir.sh
```
By using ``gen_structure_AMBER`` Python script in ERmod, 
the directories for computing rho in solution (soln) and in reference (refs) are prepared.
The directory for standard-state correction (dGcorr) is also prepared at this stage.

### Refs calculation

```
cd refs
ls
> gen_input.sh  
> run.sh
> Long_Sample  

cat gen_input.sh
> #!/bin/bash
> 
> exe=../../../packages/eror/tools/GENESIS/gen_input
> genesis_log=../../02_genesis_md/host_in_water/genesis_log
> dcd=../../02_genesis_md/host_in_water/run.dcd
> dcd_gas=../../02_genesis_md/guest_in_gas/run.dcd
> 
> # 1. Prepare template input files
> #    (parameters_er)
> 
> python2 $exe                \
>   -l $genesis_log           \
>   -x $dcd                   \
>   -s $dcd_gas               \
>   -d 1                      \
>   --minenergy -60.0
> 
> # 2. Modify parameters_er to perform the test-particle insertion 
> #    to the binding pocket
> 
> # ... Original parameters_er ...
> #&ene_param
> #        slttype = 3,
> #        boxshp = 1,
> #        estype = 1,
> #        inptemp = 298,
> #        ljformat = 5,
> #        ljswitch = 0,
> #        upljcut = 9,
> #        lwljcut = 9,
> #        cltype = 2,
> #        elecut = 9,
> #        screen = 0.30768,
> #        splodr = 4,
> #        ms1max = 64,
> #        ms2max = 48,
> #        ms3max = 48,
> #        maxins = 1000,
> #        engdiv = 1,
> #/
> #&hist
> #      eclbin=5.0e-2, ecfbin=2.0e-3, ec0bin=2.0e-4, finfac=10.0e0,
> #      ecdmin=-60.000000, ecfmns=-0.20e0, ecdcen=0.0e0, eccore=20.0e0,
> #      ecdmax=1.0e11, pecore=200
> #/
> # ..............................
> 
> cp parameters_er parameters_er.orig
> 
> cat << EOF > parameters_er
> &ene_param
>         slttype = 3,
>         boxshp = 1,
>         estype = 1,
>         inptemp = 298,
>         ljformat = 5,
>         ljswitch = 0,
>         upljcut = 9,
>         lwljcut = 9,
>         cltype = 2,
>         elecut = 9,
>         screen = 0.30768,
>         splodr = 4,
>         ms1max = 64,
>         ms2max = 48,
>         ms3max = 48,
>         maxins = 1000,
>         engdiv = 1,
>         use_anatra  = .true.,   ! ADDED
>         insposition = 7         ! ADDED
>         stop_fail_ins = .false. ! ADDED
> /
> 
> &output_param
>   fhead = "out" ! header of several output files
> /
> 
> &output_param
>   fhead = "out"
> /
> 
> &hist
>       eclbin=5.0e-2, ecfbin=2.0e-3, ec0bin=2.0e-4, finfac=10.0e0,
>       ecdmin=-60.000000, ecfmns=-0.20e0, ecdcen=0.0e0, eccore=20.0e0,
>       ecdmax=1.0e11, pecore=200
> /
> 
> !<--- ADDED
> &input_param
>   fxyz    = "../../03_reaction_coord/02_xyz/refbcd.xyz" ! used for fitting as reference coordinate
>   fdxs    = "../../03_reaction_coord/03_sdf/sdf.dx"     ! DX-formatted sdf 
>   fprmtop = "../../01_amber_setup/host-guest_in_water_for_insertion/hg_in_water_for_ins.prmtop" ! used for extracting potential parameters
> /
> 
> &trajopt_param
>   molinfo = "../../03_reaction_coord/01_spec/guest_ins.spec" ! 1 : ins
>             "../../03_reaction_coord/01_spec/guest_fit.spec" ! 2 : fit (guest)
>             "../../03_reaction_coord/01_spec/host_fit.spec"  ! 3 : fit (host)
>             "../../03_reaction_coord/01_spec/guest_uij.spec" ! 4 : uattr (guest)  Dim # 1 
>             "../../03_reaction_coord/01_spec/host_uij.spec"  ! 5 : uattr (host)   Dim # 1
> /
> 
> &option_param
>   ndim        = 1
>   reac_coords = "ENERGY"
>   rljcut      = 9.0
>   use_sdf     = .true.
>   uid_ins     = 1  ! solute molinfo id for insertion (CoM position) 
>   uid_fit     = 2  ! solute molinfo id for fitting 
>   vid_fit     = 3  ! guest  molinfo id for fitting
> /
> 
> &energy_param1
>   uid         = 4 
>   vid         = 5  
>   calc_vdw    = .true.
>   calc_elec   = .false.
>   vdw         = "ATTRACTIVE"  ! Uattr is selected as reaction coordinate
>   umin        = -38.8630197   ! Lower limit of Uattr
>   umax        =  -9.3507904   ! Upper limit of Uattr
>   ! Inserted configurations that satisfy umin <= uattr <= umax
>   ! are used for calculating rho and chi 
>   !
> /
> ! --->
> 
> EOF

./gen_input.sh
```
After executing ``gen_input.sh``, input file for ERmod (``parameters_er``) is generated.
The part in the file sandwitched between ``! <--- ADDED`` and ``!--->`` is related with the test-particle insertion of aspirin to the cavity of CD. The other parameters are described at 
https://sourceforge.net/p/ermod/wiki/parameters-erdst/ 
```
# Execute ERmod for computing rho and chi in reference 
./run.sh
ls
> HISTORY
> LJTable
> MDinfo
> MDinfo.bak
> MolPrm1
> MolPrm2
> SltConf
> SltInfo
> engref.01    # rho
> corref.01    # chi
> ermod.log
> gen_input.sh
> parameters_er
> parameters_er.orig
> progress.tt
> run.sh
> uvrange.tt
> weight_refs  # weight file for the target trajectory
```
``engref.*``, ``corref.*``, and ``weight_refs`` are used for the solvation free-energy calculation.

### dGcorr calculation

dGcorr calculation can be done in a similar way.
```
cd ../dGcorr
./gen_input.sh

> cat parameters_er
> &ene_param
>         slttype = 3,
>         boxshp = 1,
>         estype = 1,
>         inptemp = 298,
>         ljformat = 5,
>         ljswitch = 0,
>         upljcut = 9,
>         lwljcut = 9,
>         cltype = 2,
>         elecut = 9,
>         screen = 0.30768,
>         splodr = 4,
>         ms1max = 64,
>         ms2max = 48,
>         ms3max = 48,
>         maxins = 1000,
>         engdiv = 1,
>         use_anatra  = .true.,   ! ADDED
>         calc_energy = .false.   ! ADDED
>         
> /
> 
> &output_param
>   fhead = "out" ! header of several output files
> /
> 
> &output_param
>   fhead = "out"
> /
> 
> &hist
>       eclbin=5.0e-2, ecfbin=2.0e-3, ec0bin=2.0e-4, finfac=10.0e0,
>       ecdmin=-60.000000, ecfmns=-0.20e0, ecdcen=0.0e0, eccore=20.0e0,
>       ecdmax=1.0e11, pecore=200
> /
> 
> !<--- ADDED
> &input_param
>   fxyz    = "../../03_reaction_coord/02_xyz/refbcd.xyz" ! used for fitting as reference coordinate
>   fdxs    = "../../03_reaction_coord/03_sdf/sdf.dx"     ! DX-formatted sdf 
>   fprmtop = "../../01_amber_setup/host-guest_in_water_for_insertion/hg_in_water_for_ins.prmtop" ! used for extracting potential parameters
> /
> 
> &trajopt_param
>   molinfo = "../../03_reaction_coord/01_spec/guest_ins.spec" ! 1 : ins
>             "../../03_reaction_coord/01_spec/guest_fit.spec" ! 2 : fit (guest)
>             "../../03_reaction_coord/01_spec/host_fit.spec"  ! 3 : fit (host)
>             "../../03_reaction_coord/01_spec/guest_uij.spec" ! 4 : uattr (guest)  Dim # 1 
>             "../../03_reaction_coord/01_spec/host_uij.spec"  ! 5 : uattr (host)   Dim # 1
> /
> 
> &option_param
>   ndim        = 1
>   reac_coords = "ENERGY"
>   rljcut      = 9.0
>   use_sdf     = .true.
>   uid_ins     = 1  ! solute molinfo id for insertion (CoM position) 
>   uid_fit     = 2  ! solute molinfo id for fitting 
>   vid_fit     = 3  ! guest  molinfo id for fitting
>   box_shrink  = 20.0 20.0 20.0
> /
> 
> &energy_param1
>   uid         = 4 
>   vid         = 5  
>   calc_vdw    = .true.
>   calc_elec   = .false.
>   vdw         = "ATTRACTIVE"  ! Uattr is selected as reaction coordinate
>   umin        = -38.8630197   ! Lower limit of Uattr
>   umax        =  -9.3507904   ! Upper limit of Uattr
>   ! Inserted configurations that satisfy umin <= uattr <= umax
>   ! are used for calculating rho and chi 
>   !
> /
> ! --->
```
An important parameter is ``box_shrink``. 
This specify the test-particle insertion region. 
In this case, the volume for this region is Vins=20^3 = 8000 A^3.
This value is used for computing dGcorr later. 
```
# Execute the test-particle insertion for dGcorr calculation
./run.sh
ls
> HISTORY
> LJTable
> Long_Sample
> MDinfo
> MDinfo.bak
> MolPrm1
> MolPrm2
> SltConf
> SltInfo
> engref.tt
> ermod.log
> gen_input.sh
> out.popul     # <--- Output
> parameters_er
> parameters_er.orig
> progress.tt
> run.sh
> uvrange.tt
> weight_refs

cat out.popul
>        250000.00000
>          2698.00000
```
N=250000 is the total number of insertion and 
Nb=2698 is the number of insertion to the bound state.
dGcorr is expressed as dGcorr = -kT*log(Vins/V0 * (N/Nb)) (k: Boltzmann factor, T: temperature (298 K), V0: standard volume (1661 A^3))
Since kT is ~0.592 kcal/mol, dGcorr is ~ 1.8 kcal/mol.
Next step is the rho calculation in solution (soln).
```
cd ../soln
./gen_input.sh
cat parameters_er
> &ene_param
>         slttype = 1,
>         sltspec = 2,
>         boxshp = 1,
>         estype = 1,
>         inptemp = 298,
>         ljformat = 5,
>         ljswitch = 0,
>         upljcut = 9,
>         lwljcut = 9,
>         cltype = 2,
>         elecut = 9,
>         screen = 0.30768,
>         splodr = 4,
>         ms1max = 64,
>         ms2max = 48,
>         ms3max = 48,
>         engdiv = 1,
>         use_eror = .true. ! ADDED
> /
> &hist
>       eclbin=5.0e-2, ecfbin=2.0e-3, ec0bin=2.0e-4, finfac=10.0e0,
>       ecdmin=-60.000000, ecfmns=-0.20e0, ecdcen=0.0e0, eccore=20.0e0,
>       ecdmax=1.0e11, pecore=200
> /
> ! <--- ADDED
> &eror
>       refs_filename = '../refs/engref'
> /
> ! --->
```

### Soln calculation

Mode of using EROR can be turned on with ``use_eror = .true.``. 
Then, the OR state is constructed based on the rho in refs (``refs_filename = '../refs/engref'``). 
```
# Execute the interaction energy calculation in the solution (=> rho calculation)
./run.sh
ls 
> HISTORY
> LJTable
> Long_Sample
> MDinfo
> MDinfo.bak
> MolPrm1
> MolPrm2
> Prob_Result   # <--- Output: Prob. of finding OR state in B state
> SltInfo
> aveuv.tt
> avint.tt      # <--- Output: Average interaction energy of solute with surrounding
> engint.01     # <--- Output: rho  
> engint.tt
> engsln.01
> engsln.tt
> ermod.log
> flcuv.tt
> gen_input.sh
> parameters_er
> parameters_er.orig
> run.sh
> uvrange.tt
> weight_soln  # <--- Output:  # weight file for the target trajectory
``` 
``engint.*`` and ``weight_soln`` are used for the solvation free-energy calculation. The free energy change associated with the transition from the OR to solution is described in ``Prob_Result``
```
cat Prob_Result
> Probability to find the system in the OR state =     0.040000
> Free-energy change from the OR state to the simulated state
                                                =    -1.906187 kcal/mol
```

### Solvation free-energy calculation

Change the directory
```
cd ..  # 04_er directory
```
We have computed all the quantities (rho and chi in reference and rho in solution) required for the free-energy calculation.
However, the trajectory used in this tutorial is too short, the calculation can not be done properly. 
Then, let's use the data obtained using the long timescale trajectories. The location of the data is  
```
ls soln/Long_Sample
> avint.tt  dGrel.09  dGrel.18   engint.02  engint.11  engint.20
> dGrel.01  dGrel.10  dGrel.19   engint.03  engint.12  engint.21
> dGrel.02  dGrel.11  dGrel.20   engint.04  engint.13  engint.22
> dGrel.03  dGrel.12  dGrel.21   engint.05  engint.14  engint.23
> dGrel.04  dGrel.13  dGrel.22   engint.06  engint.15  engint.24
> dGrel.05  dGrel.14  dGrel.23   engint.07  engint.16  engint.25
> dGrel.06  dGrel.15  dGrel.24   engint.08  engint.17  weight_soln
> dGrel.07  dGrel.16  dGrel.25   engint.09  engint.18
> dGrel.08  dGrel.17  engint.01  engint.10  engint.19

ls refs/Long_Sample
> corref.01.zip  engref.01  weight_refs
```
The free-energy calculation using these files can be done using the script ``run.sh`` in ``04_refs`` directory.
```
ls ./run.sh
> #!/usr/bin/env bash
> 
> export exe=../../packages/eror/bin/slvfe
> export OMP_NUM_THREADS=128
> 
> cwd=`pwd`
> cd refs/Long_Sample
> unzip corref.01.zip
> cd $cwd
> 
> cat << EOF > parameters_fe
> &fevars
> inptemp=298.000000
> slndnspf  = 'engint'
> aveuvfile = 'avint.tt'
> solndirec = 'soln/Long_Sample'
> refsdirec = 'refs/Long_Sample'
> numdiv    = 5
> /
> 
> EOF
$exe >& slvfe.out

./run.sh
cat slvfe.out

> Soln digit : 2
> Refs digit : 2
>     1    0.1000000E+01
>     2    0.1000000E+01
>     3    0.1000000E+01
>     4    0.1000000E+01
>     5    0.1000000E+01
>     6    0.1000000E+01
>     7    0.1000000E+01
>     8    0.1000000E+01
>     9    0.1000000E+01
>    10    0.1000000E+01
>    11    0.1000000E+01
>    12    0.1000000E+01
>    13    0.1000000E+01
>    14    0.1000000E+01
>    15    0.1000000E+01
>    16    0.1000000E+01
>    17    0.1000000E+01
>    18    0.1000000E+01
>    19    0.1000000E+01
>    20    0.1000000E+01
>    21    0.1000000E+01
>    22    0.1000000E+01
>    23    0.1000000E+01
>    24    0.1000000E+01
>    25    0.1000000E+01
>     1    0.4000000E-01
>     2    0.4000000E-01
>     3    0.4000000E-01
>     4    0.4000000E-01
>     5    0.4000000E-01
>     6    0.4000000E-01
>     7    0.4000000E-01
>     8    0.4000000E-01
>     9    0.4000000E-01
>    10    0.4000000E-01
>    11    0.4000000E-01
>    12    0.4000000E-01
>    13    0.4000000E-01
>    14    0.4000000E-01
>    15    0.4000000E-01
>    16    0.4000000E-01
>    17    0.4000000E-01
>    18    0.4000000E-01
>    19    0.4000000E-01
>    20    0.4000000E-01
>    21    0.4000000E-01
>    22    0.4000000E-01
>    23    0.4000000E-01
>    24    0.4000000E-01
>    25    0.4000000E-01
>  
>   Number of the   1-th solvent  =            1
>   Number of the   2-th solvent  =         7200
>  
>   Self-energy of the solute   =        -0.0049  kcal/mol
>  
>  
>  cumulative average & 95% error for solvation energy
>               total             1st component         2nd component
>   1   -51.3349              -25.5660              -25.7640
>   2   -51.2051     0.2597   -25.5169     0.0982   -25.6833     0.1616
>   3   -51.1239     0.2209   -25.5351     0.0673   -25.5839     0.2195
>   4   -51.2096     0.2319   -25.5195     0.0569   -25.6852     0.2551
>   5   -51.1210     0.2523   -25.4947     0.0665   -25.6215     0.2352
>  
>  
>  group    solvation free energy     error          difference
>   total solvation free energy
>    1           -14.22629           0.21784          -0.05585
>    2           -14.16333           0.17661           0.00711
>    3           -14.17043           0.17493           0.00000
>    4           -14.10345           0.15803           0.06699
>    5           -14.13441           0.15385           0.03602
>  
>   contribution from  1-th solvent component
>    1            -7.98671           0.07975          -0.05087
>    2            -7.93023           0.04198           0.00561
>    3            -7.93584           0.03879           0.00000
>    4            -7.88481           0.02932           0.05104
>    5            -7.90776           0.02831           0.02808
>  
>   contribution from  2-th solvent component
>    1            -6.23472           0.14139          -0.00499
>    2            -6.22824           0.14104           0.00149
>    3            -6.22973           0.14299           0.00000
>    4            -6.21379           0.13985           0.01595
>    5            -6.22179           0.13729           0.00795
>  
>  
>  group   Estimated free energy: total (kcal/mol)
>    1       -14.3655     -14.1522     -14.0599     -14.5766     -13.9772
>    2       -14.2713     -14.1128     -14.0395     -14.4458     -13.9473
>    3       -14.2746     -14.1218     -14.0489     -14.4513     -13.9557
>    4       -14.1976     -14.0622     -13.9988     -14.3552     -13.9035
>    5       -14.2266     -14.0945     -14.0331     -14.3789     -13.9389
>  
>  group   Estimated free energy: 1-th solvent contribution (kcal/mol)
>    1        -8.0636      -7.9414      -7.9435      -8.0972      -7.8879
>    2        -7.9758      -7.9084      -7.9294      -7.9735      -7.8640
>    3        -7.9767      -7.9160      -7.9392      -7.9743      -7.8729
>    4        -7.9169      -7.8717      -7.9018      -7.8996      -7.8340
>    5        -7.9389      -7.8951      -7.9261      -7.9199      -7.8588
>  
>  group   Estimated free energy: 2-th solvent contribution (kcal/mol)
>    1        -6.2970      -6.2059      -6.1116      -6.4745      -6.0845
>    2        -6.2906      -6.1995      -6.1052      -6.4674      -6.0785
>    3        -6.2930      -6.2009      -6.1048      -6.4721      -6.0779
>    4        -6.2758      -6.1856      -6.0922      -6.4507      -6.0647
>    5        -6.2829      -6.1946      -6.1021      -6.4541      -6.0752
>  
>  
>  cumulative average & 95% error for solvation free energy
>               total             1st component         2nd component
>   1   -14.2746               -7.9767               -6.2930
>   2   -14.1982     0.1528    -7.9464     0.0607    -6.2469     0.0921
>   3   -14.1484     0.1330    -7.9440     0.0354    -6.1995     0.1086
>   4   -14.2241     0.1783    -7.9516     0.0292    -6.2677     0.1564
>   5   -14.1704     0.1749    -7.9358     0.0388    -6.2297     0.1430
```

The last line in ``slvfe.out`` shows the average value of solvation free energy and its 95% error. 
