! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
! Python interface for the Usadel1 solver core.
!

!------------------------------------------------------------------------------
! Set parameters
!------------------------------------------------------------------------------

subroutine spectral2_set_params( nw, nt, &
     w_length, w_conduct, w_inel, w_sf, w_phase_jump, &
     t_delta, t_phase, t_inel, t_sf, t_t, t_mu)
  use params
  use error
  implicit none
  integer, intent(in) :: nw, nt
  double precision, dimension(nw), intent(in) :: &
       w_length, w_conduct, w_inel, w_sf, w_phase_jump
  double precision, dimension(nt), intent(in) :: &
       t_delta, t_phase, t_inel, t_sf, t_t, t_mu

  call error_clear()

  ! initial check
  if (nw .le. 0 .or. nt .le. 0) then
     errorflag = -1
     errmsg = "no wires or no terminals given"
     return
  end if

  ! fill in parameters
  call resize_params(nw, nt)
  wire_length(1:nw) = w_length(1:nw)
  wire_conduct(1:nw) = w_conduct(1:nw)
  wire_rate_inelastic(1:nw) = w_inel(1:nw)
  wire_rate_spinflip(1:nw) = w_sf(1:nw)
  wire_phase_jump(1:nw) = w_phase_jump(1:nw)

  term_delta(1:nt) = t_delta(1:nt)
  term_phase(1:nt) = t_phase(1:nt)
  term_rate_inelastic(1:nt) = t_inel(1:nt)
  term_rate_spinflip(1:nt) = t_sf(1:nt)
  term_t(1:nt) = t_t(1:nt)
  term_mu(1:nt) = t_mu(1:nt)
 
end subroutine spectral2_set_params

!------------------------------------------------------------------------------
! Set coefficients
!------------------------------------------------------------------------------

subroutine spectral2_set_delta(nw, nx, x, w_delta, w_phase)
  use params
  use error
  implicit none
  integer, intent(in) :: nx, nw
  double precision, dimension(nx), intent(in) :: x
  double precision, dimension(nx, nw), intent(in) :: w_delta, w_phase

  call error_clear()

  if (.not. check_nwire(nw)) return

  call set_delta(nx, x, w_delta, w_phase)
end subroutine spectral2_set_delta

subroutine spectral2_set_kinetic( maxne, nx, nw, e, x, &
     w_dl, w_dt, w_tt, w_js, w_ddl, w_ddt, w_dtt, w_djs, w_ctl, w_ctt)
  use params
  use error
  implicit none
  integer, intent(in) :: nx, maxne, nw
  double complex, dimension(maxne), intent(in) :: e
  double precision, dimension(nx), intent(in) :: x
  double precision, dimension(nx, nw, maxne), intent(in) :: &
       w_dl, w_dt, w_tt, w_js, w_ddl, w_ddt, w_dtt, w_djs, w_ctl, w_ctt
  
  call error_clear()

  if (.not. check_nwire(nw)) return

  call set_energy(maxne, e)
  call set_kinetic(maxne, nx, x, &
       w_dl, w_dt, w_tt, w_js, w_ddl, w_ddt, w_dtt, w_djs, w_ctl, w_ctt)
end subroutine spectral2_set_kinetic

!------------------------------------------------------------------------------
! Initialize
!------------------------------------------------------------------------------

subroutine spectral2_sp_initialize(nw, mstar_, &
     w_type, bctype_new, bcdata_new, bcconnect_new)
  use params
  use sp_equations
  use error
  implicit none
  integer, intent(in) :: nw, mstar_
  integer, dimension(nw), intent(in) :: w_type
  integer, dimension(mstar_), intent(in) :: bctype_new
  double precision, dimension(mstar_*16), intent(in) :: bcdata_new
  integer, dimension(mstar_*(nw+1)), intent(in) :: bcconnect_new

  call error_clear()

  if (.not. check_nwire(nw)) return

  ! initialize equations
  call init_equations(w_type)

  if (.not. check_mstar(mstar, mstar_)) return

  bctype = bctype_new(1:mstar)
  bcdata = reshape(bcdata_new, (/ mstar, 16 /) )
  bcconnect = reshape(bcconnect_new, (/ mstar, nwire+1 /) )

  ! init solver
  if (sp_solver_type == 1) then
     call sp_init_1()
  else if (sp_solver_type == 2) then
     call sp_init_2()
  else
     errorflag = -1
     errmsg = "Invalid spectral solver chosen"
     return
  end if
end subroutine spectral2_sp_initialize

subroutine spectral2_kin_initialize(&
     nw, mstar_, w_type, bctype_new, bcdata_new, bcconnect_new, nrec)
  use params
  use kin_equations
  use error
  implicit none
  integer, intent(in) :: nrec, nw, mstar_
  integer, dimension(nw), intent(in) :: w_type
  integer, dimension(mstar_), intent(in) :: bctype_new
  double precision, dimension(mstar_*16), intent(in) :: bcdata_new
  integer, dimension(mstar_*(nw+1)), intent(in) :: bcconnect_new

  call error_clear()

  if (.not. check_nwire(nw)) return

  ! initialize equations
  call init_equations(w_type)

  if (.not. check_mstar(mstar, mstar_)) return

  bctype = bctype_new(1:mstar)
  bcdata = reshape(bcdata_new, (/ mstar, 16 /) )
  bcconnect = reshape(bcconnect_new, (/ mstar, nwire+1 /) )

  ! init solver
  if (kin_solver_type == 1) then
     call kin_init_1(nrec)
  else if (kin_solver_type == 2) then
     call kin_init_2(nrec)
  else if (kin_solver_type == 3) then
     call kin_init_3(nrec)
  else
     errorflag = -1
     errmsg = "Invalid kinetic solver chosen"
     return
  end if
end subroutine spectral2_kin_initialize

!------------------------------------------------------------------------------
! Solve
!------------------------------------------------------------------------------

subroutine spectral2_sp_solve(ne, nx, nw, ee, x, &
     a, b, da, db)
  use miscmath
  use params
  use error
  use sp_equations
  implicit none
  integer, intent(in) :: ne, nx, nw
  double complex, dimension(ne), intent(in) :: ee
  double precision, dimension(nx), intent(in) :: x
  double complex, dimension(nx, nw, ne), intent(out) :: &
       a, b, da, db
!f2py intent(inout) :: a, b, da, db

  double precision :: nan = 1d99
  
  double complex :: start_energy = 0d0, last_ok_energy = 0d0 
  double precision :: last_ok_relax = 0d0, start_relax = 0d0
  integer, parameter :: max_tries = 15
  integer :: ie, ix, k, j, krelax, jrelax, ier
  logical :: ok, reset

  integer, dimension(6), parameter :: splits = &
       (/ 1, 2, 8, 32, 128, 256 /)

  call error_clear()

  if (.not. check_nwire(nw)) return

  call set_energy(ne, ee)

  if (ne > 0) then
     start_energy = last_E
  else
     start_energy = energy(ne)
  end if

  reset = .false.

  do ie=ne,1,-1
     last_ok_energy = start_energy

     !!
     !! Energy grid loop
     !! 
     do k=1, size(splits)
        do j=1, splits(k)
           current_energy = energy(ie) - (energy(ie) - start_energy) * (splits(k) - j) / splits(k)
           !write(0,*) '% E', current_energy

           !!
           !! Relaxation loop
           !! 
           last_ok_relax = 1d0
           start_relax = 1d0
           do krelax = 1, size(splits)
              do jrelax = 0, splits(krelax)
                 if (krelax == 1 .and. jrelax == 0) cycle

                 relax = 0d0 - (0d0 - start_relax) * (splits(krelax) - jrelax)/splits(krelax)

                 !write(0,*) '  % relax', relax, splits(krelax)

                 call sp_solve_one(nx, nw, current_energy, x, &
                      a(:,:,ie), b(:,:,ie), da(:,:,ie), db(:,:,ie), &
                      sp_solver_type, reset, ier)

                 if (ier .eq. 2) then
                    reset = .true.
                 else
                    reset = .false.
                 end if

                 ok = (ier .eq. 0)
                 if (ok .and. relax .eq. 0d0) then
                    goto 200
                 else if (ok) then
                    ok = .false.
                    last_ok_relax = relax
                    !write(0,*) '  % OK at relax = ', relax, krelax, jrelax
                 else
                    ! Thicken sub-grid
                    start_relax = last_ok_relax
                    !write(0,*) '  % Try again at relax = ', relax, k, j, reset
                    exit
                 end if
              end do
           end do
           !!
           !! END relaxation loop
           !!

200        if (ok .and. current_energy == energy(ie)) then
              goto 100
           else if (ok) then
              !write(0,*) '% OK at E = ', current_energy, k, j
              ok = .false.
              last_ok_energy = current_energy
           else
              ! Thicken sub-grid
              start_energy = last_ok_energy
              !write(0,*) '% Try again at E = ', current_energy, k, j, reset
              exit
           end if
        end do
        !!
        !! END energy grid loop
        !!
     end do

100  start_energy = energy(ie)

     if (.not. ok) then
        a(:, :, ie) = nan
        b(:, :, ie) = nan
        da(:, :, ie) = nan
        db(:, :, ie) = nan
        errorflag = -1
        write(errmsg,'(a,g15.7)') "no convergence at e = ", real(energy(ie))
        return
     end if
  end do
end subroutine spectral2_sp_solve

recursive subroutine sp_solve_one(nx, nw, ee, x, a, b, da, db, &
     solver_type, reset, ier)
  use miscmath
  use params
  use error
  use sp_equations
  implicit none

  integer, intent(in) :: nx, nw
  double complex, intent(in) :: ee
  double precision, dimension(nx), intent(in) :: x
  double complex, dimension(nx, nw), intent(out) :: a, b, da, db
  integer, intent(in) :: solver_type
  integer, intent(out) :: ier
  logical, intent(in) :: reset
!f2py intent(inout) :: a, b, da, db

  logical :: ok
  double precision :: maxv

  ier = 0

  if (solver_type == 1) then
     call sp_solve_1(ok, reset)
     if (ok) then
        call sp_get_vars_1(nx, nw, x, a(:, :), da(:, :), &
             b(:, :), db(:, :))
     else
        ier = 1
     end if
  else if (solver_type == 2) then
     call sp_solve_2(ok, reset)
     if (ok) then
        call sp_get_vars_2(nx, nw,x, a(:, :), da(:, :), &
             b(:, :), db(:, :))
     else
        ier = 1
     end if
  else
     ier = 10
     errorflag = -1
     errmsg = "Invalid spectral solver chosen"
     return
  end if

  maxv = max(maxval(abs(a(:,:))), maxval(abs(b(:,:))))
  if (maxv - 1 > 1e-5) then
     !write(0,*) '% maxv', maxv
     ier = 2
  end if
end subroutine sp_solve_one

subroutine spectral2_kin_solve(ne, nx, nw, e, x, &
     fl, ft, dfl, dft, jL, jT)
  use miscmath
  use kin_equations
  use params
  use error
  implicit none
  integer, intent(in) :: ne, nx, nw
  double precision, dimension(ne), intent(in) :: e
  double precision, dimension(nx), intent(in) :: x
  double precision, dimension(nx, nw, ne), intent(out) :: &
       fl, ft, dfl, dft, jT, jL
!f2py intent(inout) :: fl, ft, dfl, dft, jT, jL

  double precision :: nan = 1d99
  
  integer :: ie, ix
  logical :: ok

  call error_clear()

  if (.not. check_nwire(nw)) return

  do ie=1,ne
     if (e(ie) .lt. real(energy(1)) .or. e(ie) .gt. real(energy(nenergy)) + 1e-6) then
        write(warnmsg,'(a,i6,a,i6,a,g15.7,a,g15.7,a,g15.7,a)') &
             "energy outside spectral range: ", &
             ie, " / ", nenergy, ": ", &
             real(energy(1)), " < ", e(ie), " < ", real(energy(nenergy)), &
             " violated"
        call push_warning(-1, warnmsg)
     end if

     current_energy = e(ie)
     if (kin_solver_type == 1) then
        call kin_solve_1(ok)
     else if (kin_solver_type == 2) then
        call kin_solve_2(ok)
     else if (kin_solver_type == 3) then
        call kin_solve_3(ok)
     else
        errorflag = -1
        errmsg = "Invalid kinetic solver chosen"
        return
     end if

     if (.not. ok) then
        fl(:, :, ie) = nan
        ft(:, :, ie) = nan
        dfl(:, :, ie) = nan
        dft(:, :, ie) = nan
        errorflag = -1
        write(errmsg,'(a,g15.7)') "no convergence at e = ", &
             real(current_energy)
        return
     else
        if (kin_solver_type == 1) then
           call kin_get_vars_1(nx, nw, x, fl(:, :, ie), dfl(:, :, ie), &
                ft(:, :, ie), dft(:, :, ie), &
                jL(:, :, ie), jT(:, :, ie))
        else if (kin_solver_type == 2) then
           call kin_get_vars_2(nx, nw, x, fl(:, :, ie), dfl(:, :, ie), &
                ft(:, :, ie), dft(:, :, ie), &
                jL(:, :, ie), jT(:, :, ie))
        else if (kin_solver_type == 3) then
           call kin_get_vars_3(nx, nw, x, fl(:, :, ie), dfl(:, :, ie), &
                ft(:, :, ie), dft(:, :, ie), &
                jL(:, :, ie), jT(:, :, ie))
        else
           errorflag = -1
           errmsg = "Invalid kinetic solver chosen"
           return
        end if
     end if
  end do
end subroutine spectral2_kin_solve

!------------------------------------------------------------------------------
! Set parameters
!------------------------------------------------------------------------------

subroutine spectral2_set_solvers(sp_type, kin_type, sp_tol, kin_tol)
  use params
  use error
  implicit none
  integer, intent(in) :: sp_type, kin_type
  double precision :: sp_tol, kin_tol

  call error_clear()

  if (sp_type .gt. 0) sp_solver_type = sp_type
  if (sp_tol .ge. 0) sp_solver_tol = sp_tol

  if (kin_type .gt. 0) kin_solver_type = kin_type
  if (kin_tol .ge. 0) kin_solver_tol = kin_tol

end subroutine spectral2_set_solvers

!------------------------------------------------------------------------------
! Error reporting
!------------------------------------------------------------------------------

subroutine spectral2_get_error(errmsg_, errorflag_, warnmsg_, warningflag_, &
     nwarnings_)
  use error
  implicit none
  character(len=380), intent(out) :: errmsg_, warnmsg_
  integer, intent(out) :: errorflag_, warningflag_, nwarnings_
  errmsg_ = errmsg
  errorflag_ = errorflag
  warnmsg_ = warnmsg
  warningflag_ = warningflag
  nwarnings_ = nwarnings
end subroutine spectral2_get_error


!------------------------------------------------------------------------------
! For each solver: solve, init_solver, get_vars
!------------------------------------------------------------------------------

! Spectral solve

subroutine sp_solve_1(ok, reset)
  use miscmath
  use sp_solve
  implicit none
  logical, intent(out) :: ok
  logical, intent(in) ::reset
  call solve(ok, reset)
end subroutine sp_solve_1

subroutine sp_solve_2(ok, reset)
  use miscmath
  use sp_solve2
  implicit none
  logical, intent(out) :: ok
  logical, intent(in) ::reset
  call solve(ok, reset)
end subroutine sp_solve_2

! Spectral init

subroutine sp_init_1()
  use miscmath
  use sp_solve
  implicit none
  call init_solver()
end subroutine sp_init_1

subroutine sp_init_2()
  use miscmath
  use sp_solve2
  implicit none
  call init_solver()
end subroutine sp_init_2

! Spectral get vars

subroutine sp_get_vars_1(nx, nw, x, a, da, b, db)
  use miscmath
  use sp_solve
  implicit none
  integer, intent(in) :: nx, nw
  double precision, dimension(nx), intent(in) :: x
  double complex, dimension(nx, nw), intent(out) :: &
       a, da, b, db
  call get_vars(x, a, da, b, db)
end subroutine sp_get_vars_1

subroutine sp_get_vars_2(nx, nw, x, a, da, b, db)
  use miscmath
  use sp_solve2
  implicit none
  integer, intent(in) :: nx, nw
  double precision, dimension(nx), intent(in) :: x
  double complex, dimension(nx, nw), intent(out) :: &
       a, da, b, db
  call get_vars(x, a, da, b, db)
end subroutine sp_get_vars_2

subroutine sp_set_vars_2(nx, nw, x, a, da, b, db)
  use miscmath
  use sp_solve2
  implicit none
  integer, intent(in) :: nx, nw
  double precision, dimension(nx), intent(in) :: x
  double complex, dimension(nx, nw), intent(out) :: &
       a, da, b, db
  call set_vars(x, a, da, b, db)
end subroutine sp_set_vars_2


! Kinetic solve

subroutine kin_solve_1(ok)
  use miscmath
  use kin_solve
  implicit none
  logical, intent(out) :: ok
  call solve(ok)
end subroutine kin_solve_1

subroutine kin_solve_2(ok)
  use miscmath
  use kin_solve2
  implicit none
  logical, intent(out) :: ok
  call solve(ok)
end subroutine kin_solve_2

subroutine kin_solve_3(ok)
  use miscmath
  use kin_solve3
  implicit none
  logical, intent(out) :: ok
  call solve(ok)
end subroutine kin_solve_3

! Kinetic init

subroutine kin_init_1(nrec)
  use miscmath
  use kin_solve
  implicit none
  integer, intent(in) :: nrec
  call init_solver(nrec)
end subroutine kin_init_1

subroutine kin_init_2(nrec)
  use miscmath
  use kin_solve2
  implicit none
  integer, intent(in) :: nrec
  call init_solver(nrec)
end subroutine kin_init_2

subroutine kin_init_3(nrec)
  use miscmath
  use kin_solve3
  implicit none
  integer, intent(in) :: nrec
  call init_solver(nrec)
end subroutine kin_init_3

! Kinetic get vars

subroutine kin_get_vars_1(nx, nw, x, fl, dfl, ft, dft, jl, jt)
  use miscmath
  use kin_solve
  implicit none
  integer, intent(in) :: nx, nw
  double precision, dimension(nx), intent(in) :: x
  double precision, dimension(nx,nw), intent(out) :: fl, dfl, ft, dft, jl, jt
  call get_vars(x, fl, dfl, ft, dft, jl, jt)
end subroutine kin_get_vars_1

subroutine kin_get_vars_2(nx, nw, x, fl, dfl, ft, dft, jl, jt)
  use miscmath
  use kin_solve2
  implicit none
  integer, intent(in) :: nx, nw
  double precision, dimension(nx), intent(in) :: x
  double precision, dimension(nx,nw), intent(out) :: fl, dfl, ft, dft, jl, jt
  call get_vars(x, fl, dfl, ft, dft, jl, jt)
end subroutine kin_get_vars_2

subroutine kin_get_vars_3(nx, nw, x, fl, dfl, ft, dft, jl, jt)
  use miscmath
  use kin_solve3
  implicit none
  integer, intent(in) :: nx, nw
  double precision, dimension(nx), intent(in) :: x
  double precision, dimension(nx,nw), intent(out) :: fl, dfl, ft, dft, jl, jt
  call get_vars(x, fl, dfl, ft, dft, jl, jt)
end subroutine kin_get_vars_3
