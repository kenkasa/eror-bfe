!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj

  ! constants
  !
  integer :: StatusUNBOUNDED = 0 
  integer :: StatusBOUNDED   = 1 


  ! subroutines
  !
  public :: analyze 

  contains
!-----------------------------------------------------------------------
    subroutine analyze(output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_output), intent(in)   :: output
      type(s_option), intent(in)   :: option
      type(s_traj),   intent(in)   :: traj(2)

      type(s_com) :: com(2)

      integer :: iatm, ires, ires_prev, id, ixyz, iz, ierr
      integer :: istep, jstep, ij, nend, itraj, nmol_large
      integer :: imol, jmol, ig
      integer :: imsta, imend, jmsta, jmend, nself
      integer :: nmol(2), nr, npair, npair_distinct 
      real(8) :: boxave(3), vol, dv, zsta, zmin, zmax,  zmax2, cmp, cmpn, cmpa
      real(8) :: m, d(3), d2, r, fourpi, fact 
      real(8) :: grsum

      logical :: is_outside
      character(len=MaxChar) :: fout

      integer, allocatable :: is_bond(:,:)

!$    write(iw,'("Analyze> OpenMP parallelization is activated")')

      ! get number of moleucles 
      !
      do itraj = 1, 2
        call get_com(option%mode(itraj), traj(itraj), com(itraj))
        nmol(itraj) = com(itraj)%nmol
      end do

      ! get average box size
      !
      do ixyz = 1, 3
        boxave(ixyz) = sum(traj(1)%box(ixyz, :)) / dble(traj(1)%nstep) 
      end do
      zmax = maxval(traj(1)%box(3,:))
      zsta = -0.5d0 * zmax 
      nr   =  zmax / option%dr + 1
      vol  =  boxave(1) * boxave(2) * boxave(3)

      ! get constants for normalization
      !
      fourpi = 4.0d0 * PI
      if (option%identical) then
        npair = nmol(1) * (nmol(1) - 1) / 2
        imsta = 1
        imend = nmol(1) - 1
        jmend = nmol(1)
      else
        npair = nmol(1) * nmol(2)
        imsta = 1
        imend = nmol(1)
        jmend = nmol(2)
      end if

      npair_distinct = nmol(1) * (nmol(1) - 1) / 2

      ! allocate memory
      !
      allocate(is_bond(jmend, traj(1)%nstep))

     ! calculation range (parameter)of coordination molecule
      is_bond                 = StatusUNBOUNDED 

      do istep = 1, traj(1)%nstep        
! $omp parallel private(imol, jmol, jmsta, d, r, ig, cmp, cmpn, is_outside), &
! $omp        & shared(istep, traj, com, is_bond) default(shared)
! $omp do
        do imol = imsta, imend

         is_outside = .false.
          if (option%use_conditional) then
            cmp  =  com(1)%coord(option%cond_type, imol, istep)
            cmpa =  abs(cmp)
            cmpn = -cmp
            is_outside = .true. 
            if (option%cond_symmetric) then
              !if ((cmp  >= option%cond_range(1) .and. cmp  <= option%cond_range(2)) .or. &
              !    (cmpn >= option%cond_range(1) .and. cmpn <= option%cond_range(2))) then
              if (cmpa  >= option%cond_range(1) .and. cmpa  <= option%cond_range(2)) then
                is_outside = .false.
              end if
            else
              if (cmp >= option%cond_range(1) .and. cmp <= option%cond_range(2)) then
                is_outside = .false.
              end if
            end if
          end if

          if (option%identical) then
            jmsta = imol + 1      
          else
            jmsta = 1
          end if

          if (is_outside) cycle

          do jmol = jmsta, jmend
            if (is_bond(jmol, istep) == StatusBOUNDED) &
              cycle

            d(1:3) = com(2)%coord(1:3, jmol, istep) &
                   - com(1)%coord(1:3, imol, istep)
            d(1:3) = d(1:3) - traj(1)%box(1:3, istep) &
                            * nint(d(1:3) / traj(1)%box(1:3, istep))
            r      = sqrt(dot_product(d, d))

            if ( r < option%bond_range ) then
              is_bond(jmol, istep) = StatusBOUNDED
              if (option%identical) &
                is_bond(imol, istep) = StatusBOUNDED
            end if
            
            !if ( r < option%bond_range ) then
            !  is_bond(jmol,istep) = StatusBOUNDED
            !else
            !  is_bond(jmol,istep) = StatusUNBOUNDED 
            !end if
         
          end do
  
        end do
! $omp end do
! $omp end parallel
      end do

      ! generate files
      !
      ! - proximal
      write(fout,'(a,".px")') trim(output%fhead) 
      open(UnitOUT, file=trim(fout))
        do istep = 1 , traj(1)%nstep
          write(UnitOUT,'(i10)',advance='no') istep 
          do jmol = jmsta, jmend
            write(UnitOUT,'(i2)',advance='no')  is_bond(jmol, istep) 
          end do
          write(UnitOUT,*)
        end do
      close(UnitOUT)

      ! - binary format RDF
      !
!      if (trim(output%frdbin) /= "") then
!       open(10, file=trim(output%frdbin),form='unformatted')
!        write(10) option%identical
!        write(10) option%separate_self
!        write(10) option%use_conditional
!        write(10) traj(1)%nstep
!        write(10) nmol
!        write(10) vol
!        write(10) option%dr
!        write(10) nr
!        close(10)
!      end if
      

      deallocate(is_bond)

    end subroutine analyze
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
