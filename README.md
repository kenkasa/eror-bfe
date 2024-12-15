# Energy-representation theory incorporating a solution state with overlapped distributions with reference (ER-OR) 

## Introduction

  ER-OR method enables us to compute binding free energy using the information on the endpoint states (solution and reference) 
  
## Programs

  * EROR: Modified ver. of ERmod 0.3.7  
    Original version is available at https://sourceforge.net/p/ermod/wiki/Home/    
    Distributed GNU GENERAL PUBLIC LICENSE Version 2, June 1991

  * ANATRA (ANAlyze TRAjectory)  
    Distributed GNU GENERAL PUBLIC LICENSE Version 2, June 1991
    In this program, following external libraries are installed

	* xdrfile-1.1.4  
	The library is distributed under the BSD License.
	https://ftp.gromacs.org/pub/contrib/

	* xdr.F90  
	This routine was originally developed by Wes Barnett and distributed under the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.  
	https://github.com/wesbarnett/libgmxfort  
	Modified version of the routine was developed by Kai-Min Tu.  
	https://github.com/kmtu/xdrfort

	* mt19937.f  
	This routine was developed by Makoto Matsumoto and Takuji Nishimura and was distributed under the GNU Library General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.   
	https://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/VERSIONS/FORTRAN/fortran.html

## Requirement 

  * VMD (ver. 1.9.3 or higher)  (for ANATRA)
  * Intel OneAPI in which ifort is available. (for ERmod and ANATRA) 
  * Login shell should be BASH (for ANATRA)

## Libraries used in ANATRA

  * xdrfile-1.1.4    
        This library is distributed under the BSD License.
	https://ftp.gromacs.org/pub/contrib/

  * xdr.F90  
	This routine was originally developed by Wes Barnett and distributed under the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
	https://github.com/wesbarnett/libgmxfort
	Modified version of the routine was developed by Kai-Min Tu.
	https://github.com/kmtu/xdrfort 

## Install 

  Execute ``install.sh`` on this directry.
  ```
  cd /path/to/eror-bfe
  ./install.sh
  ```
  Compiling EROR and ANATRA will start.
  This script registers variable ``ANATRA_PATH`` as a global variable by modifying .bashrc file, followed by the compilation of Fortran programs. After finishing the installation, please source your bashrc file.
  ```
  source ~/.bashrc
  ``` 

## Basic usage

  * [ANATRA](./packages/anatra/README.md)
  * EROR: see official [ERmod website](https://sourceforge.net/p/ermod/wiki/Home/) for the detailed description of the parameters common to ER and EROR.

  * [Tutorial (Binding of aspirin to b-cyclodextrin in water)](./tutorial/README.md) 
