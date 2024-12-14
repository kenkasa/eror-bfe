!=======================================================================
module mod_output
!=======================================================================
  use mod_util
  use mod_const

  implicit none

  ! structures
  !
  type :: s_output
    character(len=MaxChar) :: ftraj
    character(len=MaxChar) :: fdcd
    character(len=MaxChar) :: fds
    character(len=MaxChar) :: fcom, fmsd, fmsdave
    character(len=MaxChar) :: fscd, fscdave
    character(len=MaxChar) :: fpdb
    character(len=MaxChar) :: fzp
    character(len=MaxChar) :: fdist
    character(len=MaxChar) :: frd, frdbin
    character(len=MaxChar) :: fdp, fdpori
    character(len=MaxChar) :: for, fordist
    character(len=MaxChar) :: frp
    character(len=MaxChar) :: fcharge
    character(len=MaxChar) :: fweight
!    character(len=MaxChar) :: fhead_btstrp = "btstrp"
    ! Common
    character(len=MaxChar) :: fts
    character(len=MaxChar) :: fhead  
  end type s_output

  ! subroutines
  !
  public :: read_ctrl_output

  contains
!-----------------------------------------------------------------------
    subroutine read_ctrl_output(iunit, output, myrank)
!-----------------------------------------------------------------------
      implicit none
!
      integer,           intent(in)  :: iunit
      type(s_output),    intent(out) :: output
      integer, optional, intent(in)  :: myrank 

      character(len=MaxChar)      :: ftraj = ""
      character(len=MaxChar)      :: fdcd  = ""
      character(len=MaxChar)      :: fpdb  = ""
      character(len=MaxChar)      :: fds   = ""
      character(len=MaxChar)      :: fcom  = "", fmsd = "", fmsdave = ""
      character(len=MaxChar)      :: fscd  = "", fscdave = "" 
      character(len=MaxChar)      :: fzp   = ""
      character(len=MaxChar)      :: fdist = ""
      character(len=MaxChar)      :: frd   = "", frdbin  = ""
      character(len=MaxChar)      :: fdp   = "", fdpori  = ""
      character(len=MaxChar)      :: for   = "", fordist = ""
 !     character(len=MaxChar)      :: fhead_btstrp = "btstrp"
      character(len=MaxChar)      :: fts   = ""
      character(len=MaxChar)      :: frp   = ""
      character(len=MaxChar)      :: fcharge = ""
      character(len=MaxChar)      :: fweight = ""
      character(len=MaxChar)      :: fhead = "out"

      integer :: irank

      namelist /output_param/ ftraj, fdcd, fpdb, fds, fcom, fmsd, fmsdave, &
                              fscd, fscdave, fzp, fdist, &
                              frd, frdbin, &
                              fdp, fdpori, for, fordist, &
  !                            fhead_btstrp, &
                              fts, frp, fcharge, fweight, fhead


      if (present(myrank)) then
        irank = myrank
      else
        irank = 0
      end if

      ! Initialize
      !
      ftraj   = ""
      fdcd    = ""
      fpdb    = ""
      fds     = ""
      fcom    = ""
      fmsd    = ""
      fmsdave = ""
      fscd    = ""
      fscdave = ""
      fzp     = ""
      fdist   = ""
      frd     = ""
      frdbin  = ""
      fdp     = ""
      fdpori  = ""
      for     = ""
      fordist = ""
      fts     = ""
      frp     = ""
      fcharge = ""
      fweight = ""
      fhead   = "out"

      rewind iunit
      read(iunit, output_param)

      if (irank == 0) then
        write(iw,*)
        write(iw,'(">> Output section parameters")')
        if (trim(ftraj) /= "") &
          write(iw,'("ftraj   = ", a)') trim(ftraj)
        if (trim(fdcd) /= "") &
          write(iw,'("fdcd    = ", a)') trim(fdcd)
        if (trim(fpdb) /= "") &
          write(iw,'("fpdb    = ", a)') trim(fpdb)
        if (trim(fds) /= "") &
          write(iw,'("fds     = ", a)') trim(fds)
        if (trim(fcom) /= "") &
          write(iw,'("fcom    = ", a)') trim(fcom)
        if (trim(fmsd) /= "") &
          write(iw,'("fmsd    = ", a)') trim(fmsd)
        if (trim(fmsdave) /= "") &
          write(iw,'("fmsdave = ", a)') trim(fmsdave) 
        if (trim(fscd) /= "") &
          write(iw,'("fscd    = ", a)') trim(fscd) 
        if (trim(fscdave) /= "") &
          write(iw,'("fscdave = ", a)') trim(fscdave) 
        if (trim(fzp) /= "") &
          write(iw,'("fzp     = ", a)') trim(fzp) 
        if (trim(fdist) /= "") &
          write(iw,'("fdist     = ", a)') trim(fdist) 
        if (trim(frd) /= "") &
          write(iw,'("frd       = ", a)') trim(frd) 
        if (trim(frdbin) /= "") &
          write(iw,'("frdbin    = ", a)') trim(frdbin) 
        if (trim(fdp) /= "") &
          write(iw,'("fdp       = ", a)') trim(fdp) 
        if (trim(fdpori) /= "") &
          write(iw,'("fdpori    = ", a)') trim(fdpori) 
        if (trim(for) /= "") &
          write(iw,'("for       = ", a)') trim(for) 
        if (trim(fordist) /= "") &
          write(iw,'("fordist   = ", a)') trim(fordist) 
        if (trim(fts) /= "") &
          write(iw,'("fts       = ", a)') trim(fts) 
        if (trim(frp) /= "") &
          write(iw,'("frp       = ", a)') trim(frp) 
        if (trim(fcharge) /= "") &
          write(iw,'("fcharge   = ", a)') trim(fcharge) 
        if (trim(fweight) /= "") &
          write(iw,'("fweight   = ", a)') trim(fweight) 
        if (trim(fhead) /= "")   &
          write(iw,'("fhead     = ", a)') trim(fhead)
      end if

      output%ftraj   = ftraj
      output%fdcd    = fdcd
      output%fpdb    = fpdb
      output%fds     = fds
      output%fcom    = fcom
      output%fmsd    = fmsd
      output%fmsdave = fmsdave
      output%fscd    = fscd
      output%fscdave = fscdave
      output%fzp     = fzp 
      output%fdist   = fdist
      output%frd     = frd
      output%frdbin  = frdbin
      output%fdp     = fdp
      output%fdpori  = fdpori
      output%for     = for
      output%fordist = fordist
      output%fts     = fts
      output%frp     = frp
      output%fcharge = fcharge
      output%fweight = fweight
      output%fhead   = fhead

    end subroutine read_ctrl_output
!-----------------------------------------------------------------------
!

end module mod_output
!=======================================================================
