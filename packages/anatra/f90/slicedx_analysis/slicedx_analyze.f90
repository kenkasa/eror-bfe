!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_akima
  use mod_grid3d
  use mod_input
  use mod_output
  use mod_ctrl

  ! constants
  !
  real(8), parameter, private :: EPS = 1.0d-8 

  ! structures
  !

  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, cvinfo, mpl2d, bootopt)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),   intent(in)    :: input
      type(s_output),  intent(in)    :: output
      type(s_option),  intent(in)    :: option
      type(s_cvinfo),  intent(in)    :: cvinfo
      type(s_mpl2d),   intent(in)    :: mpl2d
      type(s_bootopt), intent(inout) :: bootopt

      type(s_func3d)          :: g0,  g1 

      real(8)                 :: del(3), origin(3), zcut
      integer                 :: ng3(3), ngz
      character(len=MaxChar)  :: fhead_out

      integer                 :: igx, igy, igz 

      real(8), allocatable    :: ftmp(:)


      ! Read DX file
      !
      call read_dx(input%fdx, g0) 

      ! Setup Grid info
      !
      ng3 = 1 
      ng3(1:2)    = g0%ng3(1:2)
      ngz         = g0%ng3(3)
      origin(1:3) = g0%origin(1:3)
      del(1:3)    = g0%del(1:3)
      zcut        = option%zval - g0%origin(3)

      ! Allocate memory
      !
      allocate(ftmp(1:ngz))

      write(iw,*)
      write(iw,'("Analyze> Calculate Histogram")')

      ! setup 3d-functions
      !
      call setup_func3d(ng3, del, origin, g1) 

      write(iw,*)
      write(iw,'("Analyze> Spline interpolation")')

      do igy = 1, ng3(2)
        do igx = 1, ng3(1)
          ftmp(1:ngz) = g0%data(igx, igy, 1:ngz)

          call akima(ngz, 0.0d0, del(3), zcut, ftmp, g1%data(igx, igy, 1)) 
        end do
      end do

      write(iw,'(">> Finished")')

      write(fhead_out,'(a)') trim(output%fhead)
      call generate_2dmap_matplotlib(fhead_out, g1) 
      !call generate_2dmap_gnuplot(fhead_out, g1)


    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
