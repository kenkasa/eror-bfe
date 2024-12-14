# ANAlysis TRAjectory (ANATRA) 

## Introduction

  ANATRA enables us to conduct the various types of the analyses of molecular simulation trajectory.
  This version is the development version and the stable version will be released soon. 
  
## Requirement 

  * VMD (ver. 1.9.3 or higher)
  * Intel OneAPI in which ifort is available.
  * Login shell should be BASH

## Library used in ANATRA

  * xdrfile-1.1.4
	The library is distributed under the BSD License.
	https://ftp.gromacs.org/pub/contrib/

  * xdr.F90
	This routine was originally developed by Wes Barnett under the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
	https://github.com/wesbarnett/libgmxfort
	Modified version of the routine was developed by Kai-Min Tu.
	https://github.com/kmtu/xdrfort 

## Install 

  Execute ``install.sh``.
  ```
  cd /path/to/anatra
  ./install.sh
  ```

  This script registers variable ``ANATRA_PATH`` as a global variable by modifying .bashrc file, followed by the compilation of Fortran programs. After finishing the installation, please source your bashrc file.
  ```
  source ~/.bashrc
  ``` 

## Basic usage (anatra, fanatra) 

  ANATRA provides two interfaces to conduct analysis, ``anatra`` and ``fanatra``.
  In the case of ``anatra``, it first call VMD for reading the analysis setting, followed by the execution of the Fortran programs.
  ``fanatra`` start the analysis by calling the Fortran programs without VMD.

  Usage of anatra command is:  
  ```
  anatra <analysis type> -opt1 parameter1 -opt2 parameter2 ... 
  ```
  Available analyses in anatra can be obtained by 
  ```
  anatra -h
  ===
  Available analysis in ANATRA
  trjconv     : Convert trajectory
  distance    : Distance Analysis
  comcrd      : Center-of-Mass Analysis
  rotation    : Rotational Correlation Analysis
  rdf         : RDF Analysis
  sdf         : SDF Analysis
  sasa        : SASA Analysis
  rmsd        : Root-Mean-Square-Deviation Analysis
  energy      : Energy Analysis
  strucsample : Structural Sampling Analysis
  lipidorder  : Lipid Scd order-parameter Analysis
  area        : Area Per Lipid Analysis
  z_deviation : Z-deviation Analysis
  z_profile   : Z-profile Analysis
  z_orient    : Orientation Analysis along z
  elecond     : Electric conductivity Analysis
  specfile    : Generate Species information files
  === 
  ```
  If you want to perform ``trjconv`` and to know the options, please execute the following.
  You can obtain the help message for the analysis.
  ```
  anatra trjconv -h 
  === 
  Info) VMD for LINUXAMD64, version 1.9.4a48 (October 13, 2020)
  Info) http://www.ks.uiuc.edu/Research/vmd/
  Info) Email questions and bug reports to vmd@ks.uiuc.edu
  Info) Please include this reference in published work using VMD:
  Info)    Humphrey, W., Dalke, A. and Schulten, K., `VMD - Visual
  Info)    Molecular Dynamics', J. Molec. Graphics 1996, 14.1, 33-38.
  Info) -------------------------------------------------------------
  Info) Multithreading available, 20 CPUs detected.
  Info)   CPU features: SSE2 SSE4.1 AVX AVX2 F16 AVX512F AVX512CD HT
  Info) Free system memory: 181GB (96%)
  Info) No CUDA accelerator devices available.
  Info) Dynamically loaded 3 plugins in directory:
  Info) /data2/kasa/packages/vmd/lib/plugins/LINUXAMD64/molfile
  AOOutline
  ============================================================
  
                      Trajectory Convert
  
  ============================================================
  Usage:
  anatra trjconv                                                            \
    -stype      <structure file type>                                       \
    -sfile      <structure file name>                                       \
    -tintype    <input trajectory file type>                                \
    -tin        <input trajectory file name>                                \
    -totype     <output trajectory file type>                               \
    -to         <output trajectory file name>                               \
    -beg        <first frame to be read> (default: 1)                       \
    -end        <last frame to be read>                                     \
                (default: 0, correspoding to last frame))                   \
    -selX       <X-th VMD selection> (X=0,1,2...)                           \
    -fit        <fit is performed or not (true or false)>                   \
    -centering  <centering is performed or not (true or false)>             \
                (default: false)                                            \
    -wrap       <wrap is performed or not (true or false)>                  \
                (default: false)                                            \
    -wrapcenter <wrapping center (origin or com) (default: fragment)>       \
    -wrapcomp   <how to wrap molecules (residue, segid,chain, or fragment)> \
                (residue or segid or chain or fragment or nomp)             \
                (default: fragment)                                         \
    -refpdb     <reference pdb file name>                                   \
                (necessary if fit = true)                                   \
    -outselid   <selection id for output molecules>                         \
    -fitselid   <selection id for fitting or centering>                     \
    -refselid   <selection id for reference>
  
  Remark:
  o If you specify fit = true & wrap = true, wrapcomp is automatically changed to 'com'
  
  Example (fitting):
  anatra trjconv                           \
    -stype      parm7                      \
    -sfile      str.prmtop                 \
    -tintype    dcd                        \
    -tin        inp.dcd                    \
    -totype     dcd                        \
    -to         out.dcd                    \
    -beg        1                          \
    -end        150                        \
    -sel0       not water                  \
    -sel1       resid 1 to 275 and name CA \
    -sel2       resid 1 to 275 and name CA \
    -fit        true                       \
    -centering  false                      \
    -wrap       true                       \
    -wrapcenter origin                     \
    -wrapcomp   fragment                   \
    -refpdb     ref.pdb                    \
    -outselid   0                          \
    -fitselid   1                          \
    -refselid   2
  
  Info) VMD for LINUXAMD64, version 1.9.4a48 (October 13, 2020)
  Info) Exiting normally.
  === 
  ```

## Common options 

  * -stype \<structure file type\>
  
    Structure file type. (ex. pdb, psf, parm7, gro)
  
  * -sfile \<structure file name\>
  
    Structure file 
  
  * -tintype \<input trajectory file type\>
 
    Input trajectory type 
  
  * -tin \<input  trajectory file name\>
 
    Input trajectory filename 
  
  * -totype \<trajectory output file type\>
 
    Output trajectory file type  
  
  * -to \<output trajectory file type\>
 
    Output trajectory filename 
  
  * -selX \<X-th vmd selection\> (X=0, 1, 2, ...)
 
    VMD selection for the target molecuels (ex. ``-sel0 name CA``) 

  * -dt \<time step\>

    Timestep (ns)

