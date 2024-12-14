!=======================================================================
module mod_analyze_str
!=======================================================================
  use mod_util
  use mod_const

  ! constants
  !

  ! structures
  !
  type :: s_state
    integer :: nstep
    integer :: step_final
    integer :: ndim
    integer :: nBevents   = 0
    integer :: nUBevents  = 0
    logical :: is_reacted = .false.
    logical :: is_reacted_final = .false.
    integer, allocatable :: data(:)
    real(8), allocatable :: data_ts(:) 
    real(8), allocatable :: avedata(:) 
    real(8), allocatable :: hist(:, :, :)
    real(8), allocatable :: norm(:, :)
  end type s_state

  type :: s_booteach
    integer              :: ntrial
    integer              :: nsample
    integer              :: nt_range
    integer, allocatable :: rand(:, :)
    real(8), allocatable :: taud(:), kins(:)
    real(8), allocatable :: hist(:, :, :, :)
  end type s_booteach

  type :: s_bootave
    integer              :: nt_range
    real(8)              :: taud, taud_stdev, taud_err
    real(8)              :: kins, kins_stdev, kins_err
    real(8), allocatable :: hist(:, :, :), hist_err(:, :, :), hist_stdev(:, :, :)
    real(8), allocatable :: cumm(:)
  end type s_bootave

  ! subroutines
  !

end module mod_analyze_str
!=======================================================================
