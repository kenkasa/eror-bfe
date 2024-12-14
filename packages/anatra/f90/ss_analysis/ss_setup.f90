!=======================================================================
module mod_setup
!=======================================================================
  use mod_const
  use mod_ctrl
  use mod_dcdio
  implicit none


  ! subroutines
  !
  public :: setup

  contains
!-----------------------------------------------------------------------
    subroutine setup(input, dcd)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)   :: input
      type(s_dcd),     intent(out)  :: dcd


      call read_dcd(input%fdcd(1), 0, dcd)

    end subroutine setup
!-----------------------------------------------------------------------

end module mod_setup
!=======================================================================
