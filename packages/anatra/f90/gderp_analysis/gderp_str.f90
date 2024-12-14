!=======================================================================
module mod_gderp_str
!=======================================================================
  implicit none

  type :: s_gde
    real(8) :: dt
    integer :: nstep
    integer :: ncut 
    integer :: nextend

    real(8), allocatable :: p(:) 
    real(8), allocatable :: pdot(:)
    real(8), allocatable :: kernel(:)

    real(8), allocatable :: p_ext(:)
    real(8), allocatable :: pdot_ext(:)

    real(8), allocatable :: tau_p(:)
    real(8), allocatable :: tau_k(:)

  end type s_gde

end module mod_gderp_str
!=======================================================================
