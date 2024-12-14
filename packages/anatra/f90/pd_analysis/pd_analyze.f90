!=======================================================================
module mod_analyze
!=======================================================================
!$ use omp_lib
  use mod_util
  use mod_const
  use mod_ctrl
  use mod_traj
  use mod_xtcio
  use mod_dcdio
  use mod_com

  ! subroutines
  !
  public  :: analyze
  private :: get_dist
  private :: get_mindist

  contains
!-----------------------------------------------------------------------
    subroutine analyze(input, output, option, traj)
!-----------------------------------------------------------------------
      implicit none

      type(s_input),  intent(in)    :: input
      type(s_output), intent(in)    :: output
      type(s_option), intent(in)    :: option
      type(s_traj),   intent(inout) :: traj(2)

      type(s_dcd)   :: dcd 
      type(xtcfile) :: xtc 
      type(s_com)   :: com(2)

      integer                :: i
      integer                :: istep, istep_tot, istep_use, nstep_tot
      integer                :: iatm, ires, ires_prev, ixyz
      integer                :: itraj, iunit, trajtype
      integer                :: imol, jmol, natm, nmol(2)
      real(8)                :: m, d2, d(3), ti, tj
      character(len=MaxChar) :: fout
      logical                :: is_end

      real(8), allocatable :: dist(:, :)
      real(8), allocatable :: dist_closest_pair(:)


      do itraj = 1, 2
        call get_com(option%mode(itraj),   &
                     traj(itraj),          &
                     com(itraj),           &
                     setup = .true.,       &
                     calc_coord = .false., &
                     myrank = 0)
        nmol(itraj) = com(itraj)%nmol
      end do

      ! allocate memory
      !
      allocate(dist(nmol(2), nmol(1)))
      if (option%distance_type == DistanceTypeMINIMUM) then
        allocate(dist_closest_pair(nmol(1)))
      end if

      ! calculate Distance 
      !
      if (trim(output%fds) /= "") then 
        open(UnitOut,file=trim(output%fds))
      else
        write(fout,'(a,".dis")') trim(output%fhead)
        open(UnitOut,file=trim(fout))
      end if

      if (option%distance_type == DistanceTypeMINIMUM) then
        write(fout,'(a,".dis_closest")') trim(output%fhead)
        open(UnitOut2,file=trim(fout))
      end if

      ! Get trajectory types
      !
      call get_trajtype(input%ftraj(1), trajtype)

      write(iw,*)
      write(iw,'("Analyze> Start")')

      dist = 0.0d0
      istep_tot = 0
      istep_use = 0
      do itraj = 1, input%ntraj
        call open_trajfile(input%ftraj(itraj), trajtype, iunit, dcd, xtc)
        call init_trajfile(trajtype, iunit, dcd, xtc, natm)


        is_end = .false.
        istep  = 0
        do while (.not. is_end)
          istep     = istep     + 1
          istep_tot = istep_tot + 1

          if (mod(istep_tot, 1) == 100) then
            write(iw,'("Progress: ",i0)') istep_tot
          end if

          call read_trajfile_oneframe(trajtype, iunit, istep, dcd, xtc, is_end)

          if (is_end) exit

          ti = traj(1)%dt * istep_tot

          if (ti > option%t_end) then
            is_end = .true.
            exit
          end if

          if (ti >= option%t_sta .and. ti <= option%t_end) then
            istep_use = istep_use + 1
          else
            cycle  
          end if


          do i = 1, 2
            call send_coord_to_traj(1, trajtype, dcd, xtc, traj(i))

            call get_com(option%mode(i),      &
                         traj(i),             &
                         com(i),              &
                         setup = .false.,     &
                         calc_coord = .true., & 
                         myrank = 1)
          end do

          if (option%distance_type == DistanceTypeSTANDARD) then
         
            ! calculate pair distance
            !
            call get_dist(1, option%pbc, traj(1)%box(1:3, 1), &
                          nmol, com, dist)
         
          else if (option%distance_type == DistanceTypeMINIMUM) then
         
            ! calculate pair minimum distance
            !
            call get_mindist(1, option, traj(1)%box(1:3, 1),   &
                             nmol, com, traj, dist, dist_closest_pair)
             
          else if (option%distance_type == DistanceTypeINTRA) then
         
            ! calculate intramolecular distance
            !
            call get_intradist(1, option%pbc, traj(1)%box(1:3, 1), &
                               nmol, com, dist)
         
          end if
         
          ! write pair distance at istep  
          !
         
          if (     option%distance_type == DistanceTypeSTANDARD      &
              .or. option%distance_type == DistanceTypeMINIMUM) then
            !write(UnitOut,'(f20.10)',advance="no") traj(1)%dt * istep_tot 
            write(UnitOut,'(f20.10)',advance="no") traj(1)%dt * istep_use
            do imol = 1, nmol(1)
              do jmol = 1, nmol(2)
                write(UnitOut,'(f20.10)',advance="no") dist(jmol, imol)
              end do
            end do
            write(UnitOut,*)
         
            if (option%distance_type == DistanceTypeMINIMUM) then
              !write(UnitOut2,'(f20.10)',advance="no") traj(1)%dt * istep_tot 
              write(UnitOut2,'(f20.10)',advance="no") traj(1)%dt * istep_use
              do imol = 1, nmol(1)
                write(UnitOut2,'(f20.10)',advance="no") &
                  dist_closest_pair(imol)
              end do
              write(UnitOut2,*)
            end if
         
          else if (option%distance_type == DistanceTypeINTRA) then
            !write(UnitOut,'(f20.10)',advance="no") traj(1)%dt * istep_tot 
            write(UnitOut,'(f20.10)',advance="no") traj(1)%dt * istep_use
            do imol = 1, nmol(1)
              write(UnitOut,'(f20.10)',advance="no") dist(imol, imol)
            end do
            write(UnitOut,*)
          end if

        end do

        istep_tot = istep_tot - 1

      end do

      close(UnitOut)
      if (option%distance_type == DistanceTypeMINIMUM) &
        close(UnitOut2)

      ! deallocate memory
      !
      deallocate(dist)

      if (option%distance_type == DistanceTypeMINIMUM) &
        deallocate(dist_closest_pair)


    end subroutine analyze
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_dist(istep, pbc, box, nmol, com, dist)
!-----------------------------------------------------------------------
      implicit none

      integer,      intent(in)  :: istep
      logical,      intent(in)  :: pbc
      real(8),      intent(in)  :: box(3)
      integer,      intent(in)  :: nmol(2)
      type(s_com),  intent(in)  :: com(2)
      real(8),      intent(out) :: dist(:,:)

      integer :: imol, jmol
      real(8) :: d(3) 

     
      dist = 0.0d0

      if(pbc) then
        do imol = 1, nmol(1)
          do jmol = 1, nmol(2)
            d(1:3) = com(2)%coord(1:3, jmol, istep) &
                   - com(1)%coord(1:3, imol, istep)
            d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))

            dist(jmol, imol)  = sqrt(dot_product(d, d))
          end do
        end do
      else
        do imol = 1, nmol(1)
          do jmol = 1, nmol(2)
            d(1:3) = com(2)%coord(1:3, jmol, istep) &
                   - com(1)%coord(1:3, imol, istep)

            dist(jmol, imol) = sqrt(dot_product(d, d))
          end do
        end do
      end if

!
    end subroutine get_dist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_mindist(istep, option, box, nmol, com, traj, dist, &
                           dist_closest_pair)
!-----------------------------------------------------------------------
      implicit none

      integer,        intent(in)  :: istep
      type(s_option), intent(in)  :: option
      real(8),        intent(in)  :: box(3)
      integer,        intent(in)  :: nmol(2)
      type(s_com),    intent(in)  :: com(2) 
      type(s_traj),   intent(in)  :: traj(2)
      real(8),        intent(out) :: dist(:,:)
      real(8),        intent(out) :: dist_closest_pair(:)

      integer :: iatm, jatm, imol, jmol
      real(8) :: d(3), r, dist_min 

     
      dist = 0.0d0

      if (option%pbc) then
        if (option%mindist_type(1) == MinDistTypeSITE &
            .and. option%mindist_type(2) == MinDistTypeSITE) then
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do iatm = com(1)%molsta(imol), com(1)%molend(imol)
                do jatm = com(2)%molsta(jmol), com(2)%molend(jmol)
                  d(1:3) = traj(2)%coord(1:3, jatm, istep) &
                         - traj(1)%coord(1:3, iatm, istep)
                  d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))
                  r      = sqrt(dot_product(d, d))
         
                  if (r <= dist_min) &
                    dist_min = r
         
                end do
              end do
              dist(jmol, imol) = dist_min
         
            end do
          end do
        else if (option%mindist_type(1) == MinDistTypeSITE &
            .and. option%mindist_type(2) == MinDistTypeCOM) then
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do iatm = com(1)%molsta(imol), com(1)%molend(imol)
                d(1:3) = com(2)%coord(1:3, jmol, istep) &
                       - traj(1)%coord(1:3, iatm, istep)
                d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))
                r      = sqrt(dot_product(d, d))
         
                if (r <= dist_min) &
                  dist_min = r
         
              end do
              dist(jmol, imol) = dist_min
         
            end do
          end do

        else if (option%mindist_type(1) == MinDistTypeCOM   &
          .and.  option%mindist_type(2) == MinDistTypeSITE) then 
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do jatm = com(2)%molsta(jmol), com(2)%molend(jmol)
                d(1:3) = traj(2)%coord(1:3, jatm, istep) &
                       - com(1)%coord(1:3, imol, istep) 
                d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))
                r      = sqrt(dot_product(d, d))
         
                if (r <= dist_min) &
                  dist_min = r
         
              end do
              dist(jmol, imol) = dist_min
         
            end do
          end do
        end if
      else
        if (option%mindist_type(1) == MinDistTypeSITE &
            .and. option%mindist_type(2) == MinDistTypeSITE) then
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do iatm = com(1)%molsta(imol), com(1)%molend(imol)
                do jatm = com(2)%molsta(jmol), com(2)%molend(jmol)
                  d(1:3) = traj(2)%coord(1:3, jatm, istep) &
                         - traj(1)%coord(1:3, iatm, istep)
                  r      = sqrt(dot_product(d, d))
         
                  if (r <= dist_min) &
                    dist_min = r
         
                end do
              end do
              dist(jmol, imol) = dist_min
            end do
          end do
        else if (option%mindist_type(1) == MinDistTypeSITE &
            .and. option%mindist_type(2) == MinDistTypeCOM) then
          do imol = 1, nmol(1)
            do jmol = 1, nmol(2)
              dist_min = 1.0d10 
              do iatm = com(1)%molsta(imol), com(1)%molend(imol)
                d(1:3) = com(2)%coord(1:3, jmol, istep) &
                       - traj(1)%coord(1:3, iatm, istep)
                r      = sqrt(dot_product(d, d))
         
                if (r <= dist_min) &
                  dist_min = r
         
              end do
              dist(jmol, imol) = dist_min
            end do
          end do
        else if (option%mindist_type(1) == MinDistTypeCOM) then 
          write(iw,'("Get_Mindist> Error.")')
          write(iw,'("mindist_type(1) = COM is not supported")')
          stop
        end if
      end if

      do imol = 1, nmol(1)
        dist_closest_pair(imol) = minval(dist(:, imol))
      end do

!
    end subroutine get_mindist
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
    subroutine get_intradist(istep, pbc, box, nmol, com, dist)
!-----------------------------------------------------------------------
      implicit none

      integer,      intent(in)  :: istep
      logical,      intent(in)  :: pbc
      real(8),      intent(in)  :: box(3)
      integer,      intent(in)  :: nmol(2)
      type(s_com),  intent(in)  :: com(2)
      real(8),      intent(out) :: dist(:, :)

      integer :: imol
      real(8) :: d(3) 

     
      dist = 0.0d0

      if(pbc) then
        do imol = 1, nmol(1)
          d(1:3) = com(2)%coord(1:3, imol, istep) &
                 - com(1)%coord(1:3, imol, istep)
          d(1:3) = d(1:3) - box(1:3) * nint(d(1:3) / box(1:3))

          dist(imol, imol)  = sqrt(dot_product(d, d))
        end do
      else
        do imol = 1, nmol(1)
          d(1:3) = com(2)%coord(1:3, imol, istep) &
                 - com(1)%coord(1:3, imol, istep)

          dist(imol, imol) = sqrt(dot_product(d, d))
        end do
      end if

!
    end subroutine get_intradist
!-----------------------------------------------------------------------

end module mod_analyze
!=======================================================================
