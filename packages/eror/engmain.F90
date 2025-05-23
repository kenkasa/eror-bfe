! -*- F90 -*-
! ERmod - Energy Representation Module
! Copyright (C) 2000- The ERmod authors
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


! engmain.f90: various parameters about energy calculation
!
! Angstrom and kcal/mol are taken as the units of length and energy.
!
!  names of parameters
!   numtype : number of molcular types
!   nummol : total number of molecules
!   numatm : total number of atomic sites
!   maxcnf : maximum number of configurations in MD
!   engdiv : number of divisions of the total simulation length
!   skpcnf : interval to skip the configuration examined
!   corrcal : calculation of the correlation matrix
!               0 : no calculation  1 : calculation performed
!   selfcal : construction of the self-energy distribution
!               0 (default) : no construction  1 : constructed
!   slttype : type of solute treatment
!               1 : physical
!               2 : test particle (rigid)  3 : test particle (flexible)
!            The file for the solute configuration
!            is SltInfo when slttype = 2 and is SltConf when slttype = 3
!   wgtslf : weighting by the self-energy  --- 0 : no  1 : yes
!   wgtins : weight of the solute intramolecular configuration
!               0 : no  1 : yes (can be = 1 only when slttype = 3)
!   wgtsys : weight of the solution / solvent configuration
!               0 : no  1 : yes
!   boxshp : shape of the unit cell box
!               0 : non-periodic  1 : periodic (parallelepiped or triclinic)
!   estype : type of system
!               1 : constant volume  2 : constant pressure
!
!   max_for_spec : maxinum number of species in the system
!   sltspec : specifying the solute species
!               1 <= sltspec <= numtype (default = 1) when slttype = 1
!               sltspec = numtype when slttype >= 2
!             This parameter is effective as an input only in soln calculation.
!   hostspec : solvent spcies to act as a host and bind the guest solute
!              (micelle, membrane or protein)
!               1 <= hostspec <= numtype      when slttype = 1
!               1 <= hostspec <= numtype - 1  when slttype >= 2
!   refspec : specifying the mixed solvent species for superposition reference
!               1 <= refspec <= numtype       when slttype = 1
!               1 <= refspec <= numtype - 1   when slttype >= 2
!             This parameter is effective as an input only when insorigin = 3
!
!   insorigin : translational origin of the solute position
!               The value of insorigin is set internally from insposition
!               0 : mass weighted center is moved to (0, 0, 0)
!               1 : no COM change from the value in the file to be read
!               2 : mass weighted center is moved to aggregate center
!                   (species forming the aggregate is defined by hostspec)
!               3 : fit to reference structure.
!                   reference structure needs be given as RefInfo in PDB format.
!                   RefInfo contains the structure of the host species
!                   (species forming the reference host is defined by refspec)
!                   and the solute structure, in order.
!   insposition : position for the solute
!               0 (default) : fully random position (within perodic bondary)
!               1 : no position change from the value in the file to be read
!               2 : spherically random position,
!                   with radius specified from lwreg to upreg.
!               3 : slab random position (generic case)
!                   slab geometry specified as z = com(aggregate) + dz with
!                   lwreg < dz < upreg for rectangular box periodic condition.
!                   Positioning is more complicated in parallelpiped cell.
!                   (see insertion.F90)
!               4 : slab random position (symmetric bilayer)
!                   slab geometry specified as z = com(aggregate) + dz with
!                   -upreg < dz < -lwreg or lwreg < dz < upreg
!                   for rectangular box periodic condition.
!                   Positioning is more complicated in parallelpiped cell.
!                   (see insertion.F90)
!               5 : random position relative to a reference structure
!                   solvent species identified with the refspec parameter
!                   is set to the reference structure accompanying
!                   the reference position of solute insertion
!                   and the solute is placed relative to that reference
!                   with condition of lwreg < RMSD < upreg
!               6 : (experimental) Gaussian random position.
!                   Position is given by displacing the reference coordinate,
!                   or coordinate fit to reference (insorigin = 3), with upreg.
!                   Solute weight is automatically adjusted
!               insorigin is set to 0 when insposition = 0
!                                to 1 when insposition = 1
!                                to 2 when insposition = 2, 3 or 4
!                                to 3 when insposition = 5 or 6
!   insorient : orientation for the solute
!               0 (default) : random orientation
!               1 : no orientation change from the value in the file to be read
!   insstructure : intramolecular structure of the solute
!               0 (default) : no restriction, used as is from trajectory or file
!               1 : only the structures with lwstr < RMSD < upstr is counted
!                     RefInfo needs to be prepared to determine RMSD
!
!   inscnd : (deprecated) geometrical condition of the solute configuration
!               0 : random            (insorigin = 0, insposition = 0)  default
!               1 : spherical         (insorigin = 2, insposition = 2)
!               2 : symmetric bilayer (insorigin = 2, insposition = 4)
!               3 : reference         (insorigin = 3, insposition = 5)
!   inscfg : (deprecated) position and orientation for the inserted solute
!               0 : only the intramolecular configuration is from the file.
!                   (insorient = 0)  default
!               1 : orientation is fixed from the file with random position
!                   (insorient = 1)
!               2 : position and orientation are also fixed from the file
!                   (insorient = 1, insposition = 1)
!
!   lwreg : lower bound of the region of solute position
!   upreg : upper bound of the region of solute position
!            effective only when insorigin = 2 or 3 and insposition >= 2
!   lwstr : lower bound of order parameter of solute intramolecular structure
!   upstr : upper bound of order parameter of solute intramolecular structure
!            effective only when insstructure = 1
!
!   ljformat : input-file format for the LJ energy and length parameters
!               0 : epsilon (kcal/mol) and sigma (A)
!               1 (default) : epsilon (kcal/mol) and Rmin/2 (A)
!               2 : epsilon (kJ/mol) and sigma (nm)
!               3 : A (kcal/mol A^12) and C (kcal/mol A^6)
!               4 : C12 (kJ/mol nm^12) and C6 (kJ/mol nm^6)
!               5 : Read from table, LJTable file
!                   (epsilon in kcal/mol and sigma in A)
!   ljswitch : switching function for smooth LJ truncation
!               0 (default) : potential switch in CHARMM form
!               1 : potential switch in GROMACS form
!               2 : force switch in CHARMM form
!               3 : force switch in GROMACS form
!               tapering function is defined by lwljcut and upljcut variables
!   iseed : seed parameter for uniform random number
!   inptemp : temperature of the system in Kelvin
!   temp  : temperature of the system in kcal/mol
!
!   ermax_limit : limiting the size of the distribution functions
!   block_threshold : box size for cell-link list based method in realcal.F90
!   force_calculation: if set to .true.,
!                      the program continues to run even if there is a warning
!   stdout : standard output
!
!
!  names of constants and variables for trajectory generation
!   moltype : type of the molecule numbering between 1 and numtype
!   numsite : number of sites in a molecule
!   sluvid : solvent or solute
!               0 : solvent  1 : solute
!               2 : test particle (rigid)  3 : test particle (flexible)
!   bfcoord : coordinate of a rigid molecule in body-fixed frame
!             used only for test particle and kept NULL for physical particle
!   sitemass : mass of an interaction site in a molecule
!   charge : partial charges on the sites in a molecule
!   ljene : energy parameter for Lennard-Jones potential in a molecule
!   ljlen : length parameter for Lennard-Jones potential in a molecule
!   intprm : whether the intereaction paramters given below
!                    (from elecut to ms1max,ms2max,ms3max)
!                    and boxshp, estype, and inptemp
!                    are read from the parent MD program
!      default = 0 in the case of on-the-fly calculation
!      default = 1 in the case of trajectory reading
!            on-the-fly calculation is not effective in the current version
!            and some modification of the programs is necessary
!   elecut : cutoff of the real-space part of electrostatic interaction
!   lwljcut : lower limit of the LJ cutoff tapering function
!   upljcut : upper limit of the LJ cutoff tapering function
!   cmbrule : combination rule for LJ interaction
!        0 : arithmetic mean is used for LJ sigma as for AMBER and CHARMM
!        1 : geometric mean is used for LJ sigma as for OPLS
!      default = 0
!      geometric mean is always used for LJ epsilon
!   cltype : treatment of Coulomb interaction
!        0 : bare  1 : Ewald  2 : PME  3 : PPPM
!   screen : screening parameter in Ewald, PME or PPPM
!   ewtoler : tolerance in Ewald, PME or PPPM to set the screening parameter
!      when screen is given, screen has the priority
!   splodr : order of spline function used in PME or PPPM
!   scrtype : with 'relative', Ewald tolerance refers to erfc(screen * elecut)
!             with 'distance', to erfc(screen * elecut) / elecut (NAMD type)
!   ew1max,ew2max,ew3max : number of reciprocal vectors along one direction
!   ms1max,ms2max,ms3max : number of meshes in PME or PPPM along one direction
!   specatm : specification of the site, defined as an integer function
!   sitepos : coordiantes of interaction site
!   cell : unit cell vector
!   invcl : inversion of the cell matrix
!   volume : system volume
!
!  names of constants and variables for energy distribution
!   ermax : size of the energy-represented distribution functions
!   numslv : number of solvent species
!   uvmax : number of discretization for each solvent species
!   uvsoft : number of discretization in soft interaction region
!   esmax : number of discretization of the solute self-energy
!   maxins : maximum number of insertions for test solute particle
!             This parameter is effective as an input only in refs calculation.
!   uvspec : assignment to the number representing the solvent species
!   numslt : number of solute molecules
!   sltlist : list of solute molecules
!   engnorm : normalization factor
!   engsmpl : number of samplings
!   voffset : offset value for the self-energy
!
!  constants defining the discretized energy coordinate
!     these appear only in the enginit subroutine of engproc.F90
!     and are used for each of the solute and solvent species
!     ecmns0, ecpls0, and ecfpls are internally set within enginit
!     and cannot be changed in the parameter files
!   ecprread : whether the energy parameters are read from a separate file
!       0 (default) : user-defined parameters are not read from outside
!       1 : parameters are read separately from a file EcdInfo
!   meshread : whether the energy meshes are taken from a separate file
!       0 (default) : user-defined meshes are not taken from outside
!       1 : energy coordinate meshes are read from a file EcdMesh
!   peread : (deprecated) same as ecprread
!   pecore : number of discretization in the core interaction region
!   ecdmin : minimum value of the solute-solvent energy
!   ecfmns : smaller side of the finely dicretized solute-solvent energy
!   ecmns0 : smaller side of the very finely dicretized energy near ecdcen
!   ecdcen : central value of the energy coordinate, typically zero
!   ecpls0 : larger side of the very finely dicretized energy near ecdcen
!   ecfpls : larger side of the finely dicretized solute-solvent energy
!   eccore : the solute-solvent energy at which the dicretization is changed
!   ecdmax : maximum value of the solute-solvent energy
!   eclbin : linear mesh for the solute-solvent energy
!   ecfbin : fine linear mesh for the solute-solvent energy
!   ec0bin : very fine linear mesh for the solute-solvent energy near 0
!   finfac : additional "margin" is set in the low-energy domain
!                       by shifting ecdmin and ecfmns by finfac * ecfbin
!
!
module engmain
!
! ANATRA modules
   use mod_input,        only : s_input               ! str
   use mod_output,       only : s_output              ! str
   use mod_xyz,          only : s_xyz                 ! str
   use mod_traj,         only : s_trajopt, s_traj     ! str
   use mod_anatra_ermod, only : s_option,         &   ! str
                                check_ioparams,   &   ! sbr
                                read_ctrl_option      ! sbr
   use mod_com,          only : s_com                 ! str
   use mod_grid3d,       only : s_func3d, s_spregion  ! str
   use mod_prmtop,       only : s_prmtop              ! str
   use mod_anaparm,      only : s_anaparm             ! str 
   use mod_potential,    only : s_pot                 ! str
   use mod_fitting,      only : s_fit                 ! str

!
   implicit none
   ! Note for optimization: any major compilers shall inline expand "parameter"s
   ! mathematical & physical constants
   real, parameter :: PI = 3.1415926535897932
   real(kind=8), parameter :: PI8 = 3.1415926535897932d0
   real, parameter :: cal_per_joule = 4.1840   ! thermochemical cal / J
   integer, parameter :: MaxCharER  = 512      ! ANATRA
!
   integer :: numtype, nummol, numatm, maxcnf, engdiv, skpcnf, corrcal, selfcal
   integer :: slttype, wgtslf, wgtsys, wgtins, boxshp, estype
   integer, parameter :: max_for_spec = 100    ! maximum number of species
   integer :: sltspec, hostspec(max_for_spec), refspec(max_for_spec)
   integer :: insorigin, insposition, insorient, insstructure
   integer :: sltpick, refpick, inscnd, inscfg           ! deprecated
   real :: lwreg, upreg, lwstr, upstr
   integer :: ljformat, ljswitch
   integer(8) :: iseed
   real(kind=8) :: inptemp, temp
   integer :: ermax_limit
   real :: block_threshold
   logical :: only_lj  = .false.
   logical :: use_eror = .false.
   logical :: force_calculation


   ! IO units and files
   ! system trajectory and setups
   character(len=*), parameter :: trjfile = 'HISTORY'      ! trajectory filename
   character(len=*), parameter :: inffile = 'MDinfo'       ! MD info filename
   integer, parameter :: io_MDinfo = 93                    ! MD info file IO
   character(len=*), parameter :: ene_confname = 'parameters_er'
   integer, parameter :: io_paramfile = 91
   integer, parameter :: io_flcuv = 99      ! IO for flcuv.tt and progress.tt
   integer, parameter :: io_popul = 16      !  
   integer, parameter :: stdout = 6         ! standard output

   ! interaction parameters, see setconf.F90
   character(len=*), parameter :: solute_file = 'SltInfo'  ! solute species
   character(len=*), parameter :: solvent_file = 'MolPrm'  ! solvent species
   character(len=*), parameter :: ljtable_file = 'LJTable' ! table for LJ
   integer, parameter :: mol_io = 79        ! IO for SltInfo and MolPrmX
   integer, parameter :: ljtable_io = 70    ! IO for LJ table

   ! single-solute trajectory
   !      used only when slttype = SLT_REFS_FLEX, see insertion.F90
   character(len=*), parameter :: slttrj = 'SltConf'  ! solute trajectory

   ! configuration-dependent weight for system used when wgtsys = YES
   character(len=*), parameter :: syswgt_file = "SysWght"  ! filename
   integer, parameter :: syswgt_io = 33                    ! file IO
   ! configuration-dependent weight for single solute used when wgtins = YES
   !      effective only when slttype = SLT_REFS_FLEX, see insertion.F90
   character(len=*), parameter :: sltwgt_file = 'SltWght'  ! filename
   integer, parameter :: sltwgt_io = 31                    ! file IO

   ! reference structure, see insertion.F90
   !   insorigin = INSORG_REFSTR: solvent species as superposition reference
   !      refspec : solvent species used as the reference
   !   insstructure = INSSTR_RMSD: solute itself as superposition reference
   character(*), parameter :: refstr_file = 'RefInfo'      ! structure filename
   integer, parameter :: refstr_io = 71                    ! structure file IO

   ! used-defined energy coordinate
   !      effective only when ecprread = YES or meshread = YES, see engproc.F90
   character(*), parameter :: ecdinfo_file = 'EcdInfo'     ! filename
   integer, parameter :: ecdinfo_io = 95                   ! file IO
   character(*), parameter :: ecdmesh_file = 'EcdMesh'     ! filename
   integer, parameter :: ecdmesh_io = 95                   ! file IO

   ! permutation of atom index, see setconf.F90
   character(len=*), parameter :: perm_file = "PermIndex"  ! filename
   integer, parameter :: perm_io = 75                      ! file IO


   ! <--- ANATRA
   integer, parameter :: ftrj_io    = 61
   integer, parameter :: fcv_io     = 62
   integer, parameter :: fcvread_io = 63
   ! --->
 
   ! for combination with ANATRA
   logical                            :: use_anatra    = .false.
   logical                            :: calc_energy   = .true.
   logical                            :: stop_fail_ins = .true.
 
   ! <--- ANATRA
   logical                               :: use_multiple_trj    = .false.
   logical                               :: use_selection       = .false.
   integer                               :: ntrj                = 1
   real                                  :: stat_cv_system      = 0.0
   real                                  :: uprc                = 0.0
   real                                  :: lwrc                = 0.0
   integer                               :: maxtrial            = 10000
   integer                               :: anatra_pot_type     = 0
   integer,                  allocatable :: state_count(:)     
 
   logical                               :: use_energy_criteria    = .false.
   integer                               :: energy_criteria_slvind = 10000
   real(8)                               :: energy_min             =  1.0d10
   real(8)                               :: energy_max             = -1.0d10
                                                                
   character(len=MaxCharER)              :: ftrjlist            = ""
   character(len=MaxCharER)              :: fcvlist             = ""
   character(len=MaxCharER), allocatable :: ftrj(:)
   character(len=MaxCharER), allocatable :: fcv(:)
   ! --->


   integer, dimension(:), allocatable :: moltype, numsite, sluvid
   real, dimension(:,:),  allocatable :: bfcoord
   real, dimension(:),    allocatable :: sitemass, charge, ljene, ljlen

   integer                            :: ljtype_max
   integer, dimension(:), allocatable :: ljtype
   real, dimension(:,:),  allocatable :: ljlensq_mat, ljene_mat

   real, dimension(:,:),  allocatable :: sitepos
   real, dimension(:),    allocatable :: mol_charge
   integer, dimension(:), allocatable :: mol_begin_index, belong_to
   real, dimension(3,3)               :: cell, invcl
   real, dimension(3)                 :: celllen
   real                               :: volume

   real :: elecut, lwljcut, upljcut, screen, ewtoler
   character(len=8) :: scrtype
   integer :: intprm, cmbrule, cltype, splodr
   integer :: ew1max, ew2max, ew3max, ms1max, ms2max, ms3max

   integer :: ermax, numslv, esmax, maxins
   integer, dimension(:), allocatable :: uvmax, uvsoft, uvspec
   real(kind=8), dimension(:),    allocatable :: uvcrd, edens
   real(kind=8), dimension(:,:),  allocatable :: ecorr
   real(kind=8), dimension(:),    allocatable :: escrd, eself
   real(kind=8), dimension(:,:),  allocatable :: aveuv
   real, dimension(:),    allocatable :: slnuv
   real, dimension(:,:),  allocatable :: avediv
   real(kind=8)                       :: avslf
   real, dimension(:),    allocatable :: minuv, maxuv
   integer                            :: numslt
   integer, dimension(:), allocatable :: sltlist
   real :: stat_weight_system
   real(kind=8) :: engnorm, engsmpl, voffset
   logical :: voffset_initialized = .false.

   ! ANATRA variables
   !
   type(s_input)    :: anatra_input
   type(s_output)   :: anatra_output
   type(s_prmtop)   :: prmtop
   type(s_anaparm)  :: anaparm 
   type(s_trajopt)  :: trajopt
   !type(s_traj)     :: refu, refv
   type(s_option)   :: option
   type(s_xyz)      :: vfit_xyz
   type(s_spregion) :: spregion
   type(s_pot)      :: anatra_pot(2)
   type(s_fit)      :: anatra_fitu, anatra_fitv
 
   type(s_func3d), allocatable :: sdf(:)
   type(s_traj),   allocatable :: traj(:)
   type(s_com),    allocatable :: com(:)

   ! numeric constants reference
   integer, parameter :: NO = 0, YES = 1
   integer, parameter :: SYS_NONPERIODIC = 0, SYS_PERIODIC = 1
   integer, parameter :: ES_NVT = 1, ES_NPT = 2
   integer, parameter :: LJFMT_EPS_cal_SGM_nm = 0, LJFMT_EPS_Rminh = 1, &
      LJFMT_EPS_J_SGM_A = 2, LJFMT_A_C = 3, &
      LJFMT_C12_C6 = 4, LJFMT_TABLE = 5
   integer, parameter :: LJSWT_POT_CHM = 0, LJSWT_POT_GMX = 1, &
      LJSWT_FRC_CHM = 2, LJSWT_FRC_GMX = 3
   integer, parameter :: LJCMB_ARITH = 0, LJCMB_GEOM = 1
   integer, parameter :: EL_COULOMB = 0, EL_EWALD = 1, EL_PME = 2, EL_PPPM = 3
   integer, parameter :: SLT_SOLN = 1, SLT_REFS_RIGID = 2, SLT_REFS_FLEX = 3
   integer, parameter :: PT_SOLVENT = 0, &
      PT_SOLUTE = SLT_SOLN, PT_TEST_RIGID = SLT_REFS_RIGID, &
      PT_TEST_FLEX = SLT_REFS_FLEX
   ! PT_SOLUTE to PT_TEST_FLEX should correspond to SLT_SOLN to SLT_REFS_FLEX
!   integer, parameter :: INSORG_ORIGIN = 0, INSORG_NOCHANGE= 1, &
!      INSORG_AGGCEN = 2, INSORG_REFSTR = 3
!   integer, parameter :: INSPOS_RANDOM = 0, INSPOS_NOCHANGE= 1, &
!      INSPOS_SPHERE = 2, &
!      INSPOS_SLAB_GENERIC = 3, INSPOS_SLAB_SYMMETRIC = 4, &
!      INSPOS_RMSD = 5, INSPOS_GAUSS = 6
!   integer, parameter :: INSROT_RANDOM = 0, INSROT_NOCHANGE= 1
!   integer, parameter :: INSSTR_NOREJECT = 0, INSSTR_RMSD = 1
   integer, parameter :: INSORG_ORIGIN = 0, INSORG_NOCHANGE= 1,              &
                         INSORG_AGGCEN = 2, INSORG_REFSTR = 3,               &
                         INSORG_ANATRA = 4, INSORG_ANATRA_RANDOM = 5
 
   integer, parameter :: INSPOS_RANDOM         = 0, &
                         INSPOS_NOCHANGE       = 1, &
                         INSPOS_SPHERE         = 2, &
                         INSPOS_SLAB_GENERIC   = 3, &
                         INSPOS_SLAB_SYMMETRIC = 4, &
                         INSPOS_RMSD           = 5, &
                         INSPOS_GAUSS          = 6, &
                         INSPOS_ANATRA         = 7, &
                         INSPOS_ANATRA_RANDOM  = 8
 
   integer, parameter :: INSROT_RANDOM = 0, INSROT_NOCHANGE= 1
   integer, parameter :: INSSTR_NOREJECT = 0, INSSTR_RMSD = 1
   integer, parameter :: INSSTR_ANATRA   = 2  ! ANATRA
 
   integer, parameter :: ANATRA_POT_AMBER = 0, ANATRA_POT_CHARMM = 1  



   namelist /ene_param/ iseed,                                  &
      skpcnf, corrcal, selfcal,                                 &
      slttype, wgtslf, wgtsys, wgtins, boxshp, estype,          &
      sltspec, hostspec, refspec, lwreg, upreg, lwstr, upstr,   &
      insposition, insorient, insstructure,                     &
      sltpick, refpick, inscnd, inscfg,                         & ! deprecated
      ljformat, ljswitch,                                       &
      inptemp, temp,                                            &
      engdiv, maxins,                                           &
      intprm, elecut, lwljcut, upljcut,                         &
      cmbrule, cltype, screen, ewtoler, splodr, scrtype,        &
      ew1max, ew2max, ew3max, ms1max, ms2max, ms3max,           &
      ermax_limit, block_threshold, force_calculation,          &
      use_anatra, calc_energy, stop_fail_ins,                   & ! for ANATRA
      use_eror, only_lj,                                        & ! for ANATRA
      use_multiple_trj, use_selection, ftrjlist, fcvlist, ntrj, & ! for ANATRA
      uprc, lwrc, maxtrial, anatra_pot_type,                    & ! for ANATRA
      use_energy_criteria, energy_min, energy_max,              & ! for ANATRA
      energy_criteria_slvind

contains
   subroutine init_params(directory)
      implicit none
      integer :: param_err
      character(len=*), optional :: directory
      character(len=1024) :: eneconfbuf

      if(present(directory)) then
        eneconfbuf = trim(directory) // "/" // ene_confname
      else
        eneconfbuf = ene_confname
      end if

      param_err = 0
      open(unit = io_paramfile, file = trim(eneconfbuf), action = "read", iostat = param_err)

      if(param_err == 0) then
         read(io_paramfile, nml = ene_param)
         close(io_paramfile)
      else
         stop "parameters_er file does not exist"
      end if

   end subroutine init_params

   ! ANATRA subroutines
   !
   !-----------------------------------------------------------------------------
   subroutine setup_anatra(calc_energy, myrank)
   !-----------------------------------------------------------------------------
     use mod_input,     only : read_ctrl_input  
     use mod_output,    only : read_ctrl_output
     use mod_traj,      only : read_ctrl_trajopt,    &
                               setup_traj_from_args
     use mod_com,       only : get_com
     use mod_xyz,       only : read_xyz
     use mod_grid3d,    only : read_dx, search_nonzero_region
     use mod_prmtop,    only : read_prmtop
     use mod_anaparm,   only : read_anaparm
     use mod_potential, only : alloc_pot, setup_pot
     use mod_fitting,   only : setup_fit
 
     implicit none
 
     logical, intent(in) :: calc_energy
     integer, intent(in) :: myrank
 
     integer :: i, itraj, uid, vid, ierr
     integer :: param_err
 
 
     if (myrank == 0) then
       write(6,*)
       write(6,'("================================================================================")')
       write(6,*)
       write(6,'("                Conditional Solvation Free Energy Calculation")')
       write(6,*)
       write(6,'("================================================================================")')
       write(6,*)
     end if
 
     ! Read INPUT namelist
     !   - values of fxyz, fdx, and fprmtop are read from this namelist
     !
     param_err = 0
     open(unit   = io_paramfile, &
          file   = ene_confname, &
          action = "read",       &
          iostat = param_err)
 
     if (param_err == 0) then
        call read_ctrl_input(io_paramfile, anatra_input, myrank)
     else
       stop "parameters_er file does not exist"
     end if
     close(io_paramfile)
 
 !
     ! Read OUTPUT namelist
     !
     param_err = 0
     open(unit   = io_paramfile, &
          file   = ene_confname, &
          action = "read",       &
          iostat = param_err)
 
 
     if (param_err == 0) then
       call read_ctrl_output(io_paramfile, anatra_output, myrank)
     else
       stop "parameters_er file does not exist"
     end if
     close(io_paramfile)
 
     ! Check INPUT/OUTPUT parameters specified
     !
     if (myrank == 0) then
       call check_ioparams(anatra_input, anatra_output)
     end if
 
     ! Read TRAJOPT namelist
     !
     param_err = 0
     open(unit   = io_paramfile, &
          file   = ene_confname, &
          action = "read",       &
          iostat = param_err)
 
     if (param_err == 0) then
       call read_ctrl_trajopt(io_paramfile, trajopt, myrank = myrank)
     else
       stop "parameters_er file does not exist"
     end if
     close(io_paramfile)
 
     ! Read OPT_PARAM namelist 
     ! 
     param_err = 0
     open(unit   = io_paramfile, &
          file   = ene_confname, &
          action = "read",       &
          iostat = param_err)
 
     if (param_err == 0) then
       call read_ctrl_option(io_paramfile, calc_energy, option, myrank = myrank)
     else
       stop "parameters_er file does not exist"
     end if
     close(io_paramfile)
 
     ! Allocate memory
     !
     allocate(traj(trajopt%nmolinfo))
     allocate(com(trajopt%nmolinfo))
 
     ! Setup TRAJ structures
     !
     do itraj = 1, trajopt%nmolinfo
       call setup_traj_from_args(trajopt, 1, traj(itraj), trajid = itraj)
       traj(itraj)%box(1:3, 1) = option%box(1:3)
     end do
 
     ! Setup SDF & Fit
     !
     if (option%use_sdf) then
       ! Read Reference solvent xyz file
       !
       param_err = 0
       open(unit   = io_paramfile,               &
            file   = trim(anatra_input%fxyz(1)), &
            action = "read",                     &
            iostat = param_err)
      
       if (param_err == 0) then
         call read_xyz(io_paramfile, vfit_xyz)
       else
         stop "xyz file is not found"
       end if
       close(io_paramfile)
      
       ! Read Spatial distribution function from dx file
       !
       allocate(sdf(1:option%nstate))
 
       do i = 1, option%nstate
         call read_dx(anatra_input%fdxs(i), sdf(i))
       end do
 
       ! Search Nonzero value grids in SDF
       !
       call search_nonzero_region(sdf(1), spregion, threshold = option%sdf_threshold)
 
       ! Load reference local solvent structure
       !
       traj(option%vid_fit)%coord(:, :, 1) = vfit_xyz%coord(:, :)
 
       ! Setup Fit structure
       !
       call setup_fit(traj(option%uid_fit), anatra_fitu)
       call setup_fit(traj(option%vid_fit), anatra_fitv)
     end if
 
     ! Setup state_count
     !
     allocate(state_count(0:option%nstate))
     state_count(0:option%nstate) = 0
 
     ! Read PRMTOP
     !
     if (anatra_pot_type == ANATRA_POT_AMBER) then 
       call read_prmtop(anatra_input%fprmtop, prmtop)
     else if (anatra_pot_type == ANATRA_POT_CHARMM) then
       call read_anaparm(anatra_input%fanaparm, anaparm)
     end if
 
     ! Setup COM structures
     !
     do itraj = 1, trajopt%nmolinfo
       call get_com(2, traj(itraj), com(itraj), myrank = myrank) ! Note: "2" of 1st arg. means WHOLE mode  
     end do
 
     ! Setup POT structure
     !
     do i = 1, option%ndim
       if (option%reac_coords(i) == 2) then ! correspond to "ENERGY"
         uid = option%eneopt(i)%uid
         vid = option%eneopt(i)%vid
         call alloc_pot((/traj(uid)%natm, traj(vid)%natm/), anatra_pot(i))
 
         if (anatra_pot_type == ANATRA_POT_AMBER) then
           call setup_pot((/traj(uid), traj(vid)/), anatra_pot(i), prmtop = prmtop)
         else if (anatra_pot_type == ANATRA_POT_CHARMM) then
           call setup_pot((/traj(uid), traj(vid)/), anatra_pot(i), anaparm = anaparm)
         end if
       end if
     end do
 
   !-----------------------------------------------------------------------------
   end subroutine setup_anatra 
   !-----------------------------------------------------------------------------

   ! returns atom no. for [i]th atom in [mol]th molecule
   integer function specatm(i, mol)
      implicit none
      integer, intent(in) :: i, mol
      specatm = mol_begin_index(mol) + (i - 1)
   end function specatm

   ! helper function that corresponds to mol_begin_index
   integer function mol_end_index(mol)
      implicit none
      integer, intent(in) :: mol
      mol_end_index = mol_begin_index(mol + 1) - 1
   end function mol_end_index
end module engmain
