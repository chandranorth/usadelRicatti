! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
module error

  character(len=580) :: errmsg
  integer :: errorflag

  character(len=580) :: warnmsg
  integer :: warningflag
  integer :: nwarnings

contains

  subroutine error_clear()
    implicit none
    errmsg = ""
    errorflag = 0
    warnmsg = ""
    warningflag = 0
    nwarnings = 0
  end subroutine error_clear
  
  function check_nwire(nw) result(success)
    use params
    implicit none
    integer, intent(in) :: nw
    logical :: success
    
    if (nwire .ne. nw) then
       errmsg = "Invalid value given for parameter nwire"
       errorflag = -1
       success = .FALSE.
    else
       errorflag = 0
       success = .TRUE.
    end if
  end function check_nwire


  function check_mstar(mstar, mstar_) result(success)
    use params
    implicit none
    integer, intent(in) :: mstar, mstar_
    logical :: success

    if (mstar .ne. mstar_) then
       errmsg = "Invalid value given for parameter mstar"
       errorflag = -1
       success = .FALSE.
    else
       errorflag = 0
       success = .TRUE.
    end if
  end function check_mstar

  subroutine push_warning(code, msg)
    implicit none
    integer, intent(in) :: code
    character(len=*), intent(in) :: msg

    warnmsg = msg
    warningflag = code
    nwarnings = nwarnings + 1
  end subroutine push_warning

end module error
