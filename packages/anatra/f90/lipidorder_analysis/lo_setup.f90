!=======================================================================
module mod_setup
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_dcdio
  use mod_traj
  implicit none


  ! subroutines
  !
  public :: setup

  contains
!-----------------------------------------------------------------------
    subroutine setup(input, trajopt, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),                  intent(in)  :: input
      type(s_trajopt),                intent(in)  :: trajopt
      type(s_traj),    allocatable,   intent(out) :: traj(:)
  
      integer :: itraj, nmolinfo


      if (input%ftraj(1) == "") then
        write(iw,'("Setup> Error.")')
        write(iw,'("ftraj in input_param is empty...")')
        stop
      end if

      nmolinfo = trajopt%nmolinfo

      allocate(traj(nmolinfo))

      do itraj = 1, nmolinfo 
        call setup_traj_from_args(trajopt,      &
                                  1,            &
                                  traj(itraj),  &
                                  trajid = itraj)
      end do

    end subroutine setup

!-----------------------------------------------------------------------

end module mod_setup
!=======================================================================
