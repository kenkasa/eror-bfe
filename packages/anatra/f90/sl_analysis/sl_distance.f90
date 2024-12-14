!=======================================================================
module mod_ctrl
!=======================================================================
  use mod_util
  use mod_const
  use mod_input
  use mod_output
  use mod_traj
  use mod_com

  implicit none

  ! constants
  !

  ! structures
  !

  ! subroutines
  !
  !public  :: read_ctrl
  !private :: read_ctrl_option

  contains

!-----------------------------------------------------------------------
    subroutine read_ctrl(input, output, option, trajopt)
!-----------------------------------------------------------------------
      implicit none

      integer, parameter           :: iunit = 10

      type(s_input),   intent(out) :: input
      type(s_output),  intent(out) :: output
      type(s_option),  intent(out) :: option
      type(s_trajopt), intent(out) :: trajopt

      character(len=MaxChar)       :: f_ctrl

      ! get control file name
      !
      call getarg(1, f_ctrl)

      write(iw,*)
      write(iw,'("Read_Ctrl> Reading parameters from ", a)') trim(f_ctrl)
      open(iunit, file=trim(f_ctrl), status='old')
        call read_ctrl_input  (iunit, input)
        call read_ctrl_output (iunit, output)
        call read_ctrl_option (iunit, option)
        call read_ctrl_trajopt(iunit, trajopt)
      close(iunit)

    end subroutine read_ctrl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    subroutine read_ctrl_option(iunit, option)
!-----------------------------------------------------------------------
      implicit none
!
      integer,        intent(in)  :: iunit
      type(s_option), intent(out) :: option 

      real(8)                :: dr              = 0.1d0
      character(len=MaxChar) :: mode(2)         = (/"RESIDUE", "RESIDUE"/)
      logical                :: identical       = .false.
      logical                :: normalize       = .false.
      logical                :: separate_self   = .false.

      logical                :: use_conditional = .false.
      logical                :: cond_symmetric  = .false.
      character(len=MaxChar) :: cond_type       = 'Z'
      real(8)                :: cond_range(2)   = (/0.0d0, 20.0d0/)
      real(8)                :: bond_range      = 3.5d0

      integer                :: itraj, i
      integer                :: iopt, ierr

      namelist /option_param/ dr, mode, identical, normalize, &
                              separate_self, &
                              use_conditional, cond_symmetric, &
                              cond_type, cond_range, bond_range

      rewind iunit
      read(iunit, option_param)

      write(iw,*)
      write(iw,'(">> Option section parameters")')
      write(iw,'("dr            = ", f15.7)')   dr 
      write(iw,'("mode          = ", a,2x,a)')  trim(mode(1)), trim(mode(2))
      write(iw,'("identical     = ", a)')       get_tof(identical)
      write(iw,'("normalize     = ", a)')       get_tof(normalize)
      write(iw,'("separate_self = ", a)')       get_tof(separate_self)

      write(iw,'("use_conditional = ", a)')       get_tof(use_conditional)
      if (use_conditional) then
        write(iw,'("cond_symmetric  = ", a)')       get_tof(cond_symmetric)
        write(iw,'("cond_type       = ", a)')       trim(cond_type)
        write(iw,'("cond_range      = ", 2f20.10)') (cond_range(i), i = 1, 2)
        write(iw,'("bond_range      = ", f20.10)')   bond_range
      end if


      do itraj = 1, 2
        iopt = get_opt(mode(itraj), CoMMode, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("mode = ",a," is not available.")') trim(mode(itraj))
          stop
        end if
        option%mode(itraj) = iopt
      end do
  
      if (use_conditional) then
        iopt = get_opt(cond_type, CondType, ierr)
        if (ierr /= 0) then
          write(iw,'("Read_Ctrl_Option> Error.")')
          write(iw,'("cond_type = ",a," is not available.")') trim(cond_type)
          stop
        end if
        option%cond_type = iopt
      end if

                         
      option%dr              = dr
      option%identical       = identical
      option%normalize       = normalize
      option%separate_self   = separate_self

      option%use_conditional = use_conditional
      option%cond_symmetric  = cond_symmetric
      option%cond_range      = cond_range
      option%bond_range      = bond_range

      ! Combination check
      !




    end subroutine read_ctrl_option
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
end module
!=======================================================================
