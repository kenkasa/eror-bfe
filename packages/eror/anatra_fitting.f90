!=======================================================================
module mod_fitting
!=======================================================================
  use mod_util
  use mod_const
  use mod_traj

  type :: s_fit
    integer              :: natm
    real(8), allocatable :: mass(:)
    real(8), allocatable :: refcoord(:, :)
    real(8), allocatable :: movcoord(:, :)
    real(8)              :: refcom(3)
    real(8)              :: movcom(3)
    real(8)              :: rot_matrix(3, 3)
    integer, allocatable :: ind(:)
  end type s_fit

  ! subroutines
  !
  public :: setup_fit
  public :: get_trrot
  public :: operate_trrot 

  contains

!-----------------------------------------------------------------------
    subroutine setup_fit(traj, fit, refcoord) 
!-----------------------------------------------------------------------
      implicit none

      type(s_traj),               intent(in)    :: traj
      type(s_fit),                intent(inout) :: fit
      real(8),          optional, intent(in)    :: refcoord(:, :)

      integer :: i
      integer :: natm


      natm = traj%natm

      allocate(fit%mass(natm), fit%refcoord(3, natm), fit%movcoord(3, natm))
      allocate(fit%ind(natm))
      do i = 1, natm
        fit%mass(i) = traj%mass(i)
      end do

      fit%refcoord = 0.0d0
      fit%movcoord = 0.0d0

      if (present(refcoord)) then
        fit%refcoord(:, :) = refcoord(:, :)
      end if

      fit%natm   = traj%natm
      fit%ind(:) = traj%ind(:)
!
    end subroutine setup_fit
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_trrot(fit) 
!-----------------------------------------------------------------------
      implicit none

      integer,           parameter     :: lwork = max(1,4*3-1)

      type(s_fit),       intent(inout) :: fit

      real(8) :: vec(4), eigen_work(4, 4) 

      integer :: i, natm, ierr
      real(8) :: total_mass
      real(8) :: smat(4, 4), eval(1:4), work(1:lwork), eu(1:4)
      real(8) :: dref(3), dmov(3), ds(3), da(3)


      natm      = fit%natm
      fit%refcom = 0.0d0
      fit%movcom = 0.0d0
      total_mass = 0.0d0 

      do i = 1, natm
        fit%refcom(1:3) = fit%refcom(1:3) + fit%mass(i) * fit%refcoord(1:3, i) 
        fit%movcom(1:3) = fit%movcom(1:3) + fit%mass(i) * fit%movcoord(1:3, i) 
        total_mass      = total_mass + fit%mass(i)
      end do

      total_mass = 1.0d0 / total_mass
      fit%refcom(1:3) = fit%refcom(1:3) * total_mass
      fit%movcom(1:3) = fit%movcom(1:3) * total_mass

      smat(1:4, 1:4) = 0.0d0

      do i = 1, natm
        dref(1:3) = fit%refcoord(1:3, i) - fit%refcom(1:3)
        dmov(1:3) = fit%movcoord(1:3, i) - fit%movcom(1:3)
        ds(1:3)   = fit%mass(i) * (dref(1:3) - dmov(1:3))
        da(1:3)   = fit%mass(i) * (dref(1:3) + dmov(1:3))

        smat(1,1) = smat(1,1) + ds(1)*ds(1) + ds(2)*ds(2) + ds(3)*ds(3)  
        smat(1,2) = smat(1,2) + da(2)*ds(3) - ds(2)*da(3)
        smat(1,3) = smat(1,3) + ds(1)*da(3) - da(1)*ds(3)
        smat(1,4) = smat(1,4) + da(1)*ds(2) - ds(1)*da(2)
        smat(2,2) = smat(2,2) + ds(1)*ds(1) + da(2)*da(2) + da(3)*da(3)  
        smat(2,3) = smat(2,3) + ds(1)*ds(2) - da(1)*da(2)
        smat(2,4) = smat(2,4) + ds(1)*ds(3) - da(1)*da(3)
        smat(3,3) = smat(3,3) + da(1)*da(1) + ds(2)*ds(2) + da(3)*da(3)  
        smat(3,4) = smat(3,4) + ds(2)*ds(3) - da(2)*da(3)
        smat(4,4) = smat(4,4) + da(1)*da(1) + da(2)*da(2) + ds(3)*ds(3)  
      end do

      smat(2,1) = smat(1,2)
      smat(3,1) = smat(1,3)
      smat(3,2) = smat(2,3)
      smat(4,1) = smat(1,4)
      smat(4,2) = smat(2,4)
      smat(4,3) = smat(3,4)

      call dsyev('V', 'U', 4, smat, 4, eval, work, lwork, ierr)

      eu(1:4) = smat(1:4, 1)

      fit%rot_matrix = 0.0d0

      fit%rot_matrix(1,1) = eu(1)*eu(1) + eu(2)*eu(2) - eu(3)*eu(3) - eu(4)*eu(4)
      fit%rot_matrix(1,2) = 2.0d0  * (eu(2)*eu(3) + eu(1)*eu(4))
      fit%rot_matrix(1,3) = 2.0d0  * (eu(2)*eu(4) - eu(1)*eu(3))
      fit%rot_matrix(2,1) = 2.0d0  * (eu(2)*eu(3) - eu(1)*eu(4))
      fit%rot_matrix(2,2) = eu(1)*eu(1) - eu(2)*eu(2) + eu(3)*eu(3) - eu(4)*eu(4)
      fit%rot_matrix(2,3) = 2.0d0  * (eu(3)*eu(4) + eu(1)*eu(2))
      fit%rot_matrix(3,1) = 2.0d0  * (eu(2)*eu(4) + eu(1)*eu(3))
      fit%rot_matrix(3,2) = 2.0d0  * (eu(3)*eu(4) - eu(1)*eu(2))
      fit%rot_matrix(3,3) = eu(1)*eu(1) - eu(2)*eu(2) - eu(3)*eu(3) + eu(4)*eu(4)

    end subroutine get_trrot 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine operate_trrot(fit, natm, coord) 
!-----------------------------------------------------------------------
      implicit none

      type(s_fit),       intent(inout)   :: fit
      integer,           intent(in)      :: natm
      real(8),           intent(inout)   :: coord(3, natm)

      integer :: iatm
      real(8) :: dmov(1:3), rot_mov(1:3)


      do iatm = 1, natm
      
        dmov(1:3)       = coord(1:3,iatm) - fit%movcom(1:3)
        rot_mov(1:3)    = fit%rot_matrix(1:3,1) * dmov(1) + &
                          fit%rot_matrix(1:3,2) * dmov(2) + &
                          fit%rot_matrix(1:3,3) * dmov(3)
        coord(1:3,iatm) = rot_mov(1:3) + fit%refcom(1:3)
      
      end do
!
    end subroutine operate_trrot 
!-----------------------------------------------------------------------

end module mod_fitting
!=======================================================================
