! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
MODULE KIN_SOLVE2
  USE PARAMS
  USE KIN_EQUATIONS
  use lazy_alloc

  PUBLIC SOLVE, INIT_SOLVER, GET_VARS

  PRIVATE

  INTEGER :: ntol
  INTEGER, DIMENSION(:), POINTER :: ltol
  double precision, DIMENSION(:), POINTER :: tol

  INTEGER, DIMENSION(:), POINTER :: iwrk
  double precision, DIMENSION(:), POINTER :: wrk, xx, oxx
  double precision, DIMENSION(:,:), POINTER :: u, ou
  INTEGER :: nmsh, onmsh

  INTEGER :: nmax, nucol, nlbc, lwrkfl, lwrkin

  LOGICAL :: guess_exists = .FALSE.

CONTAINS

  SUBROUTINE INIT_SOLVER(NREC)
    IMPLICIT NONE

    INTEGER :: KD, KDM, NSIZEI, NSIZEF, NMAX, NREC
    INTEGER :: I

    nmax = 500
    nucol = 2*nmax

    ! Number of left boundary conditions
    nlbc = MSTAR - NREC

    ! Boundary condition positions
    call lazy_alloc_r1(ZETA, MSTAR)
    ZETA = 1d0
    ZETA(1:(MSTAR - NREC)) = 0d0

    ! Tolerances (do not specify)
    NTOL = MSTAR
    call lazy_alloc_i1(LTOL, NTOL)
    call lazy_alloc_r1(TOL, NTOL)
    DO I = 1, MSTAR
       LTOL(I) = I
    END DO
    IF (KIN_SOLVER_TOL .LE. 0) THEN
       TOL = 1e-8
    ELSE
       TOL = KIN_SOLVER_TOL
    END IF

    ! Workspace
    lwrkfl = 5*mstar*mstar + 2*ntol + 9*mstar + nmax*( &
         4*mstar*mstar + 12*mstar + 3)
    lwrkin = 100 + mstar + nmax*(mstar + 2)

    call lazy_alloc_r1(xx, nucol)
    call lazy_alloc_r1(oxx, nucol)
    call lazy_alloc_r2(u, mstar,nucol)
    call lazy_alloc_r2(ou, mstar,nucol)
    call lazy_alloc_r1(wrk, 2*lwrkfl)
    call lazy_alloc_i1(iwrk, 2*lwrkin)
    
    guess_exists = .FALSE.
  END SUBROUTINE INIT_SOLVER

  SUBROUTINE SOLVE(CONVERGED)
    USE TWPBVPC_MOD
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: CONVERGED

    double precision :: ALEFT, ARIGHT, EPS, EPSMIN
    double precision :: FIXPNT(1)
    INTEGER :: IFLAG, nfxpnt
    LOGICAL :: giveeps, givemsh, giveu, linear
    
    double precision :: ckappa1,gamma1,ckappa, rpar(1)
    integer :: ipar(1)

    CALL INTERPOLATE_KINETIC_FOR_ENERGY(real(CURRENT_ENERGY))

    ALEFT = 0d0
    ARIGHT = 1d0

    FIXPNT = 0
    nfxpnt = 0

    linear = .TRUE.

    iprint = -1
    pdebug = .FALSE.
    uval0 = 1

    IF (guess_exists) THEN
       givemsh = .TRUE.
       giveu = .TRUE.
    ELSE
       CALL GUESS_TWP()
       givemsh = .TRUE.
       giveu = .TRUE.
    END IF

    !givemsh = .FALSE.
    !giveu = .FALSE.
    
    !CALL TWPBVP(mstar, nlbc, aleft, aright, nfxpnt, fixpnt, &
    !     ntol, ltol, tol, linear, givemsh, giveu, nmsh, &
    !     xx, mstar, u, nmax, lwrkfl, wrk, lwrkin, iwrk, &
    !     fsub_twp, dfsub_twp, gsub_twp, dgsub_twp, iflag)

    CALL TWPBVPC(mstar, nlbc, aleft, aright, nfxpnt, fixpnt, &
         ntol, ltol, tol, linear, givemsh, giveu, nmsh, &
         nucol, xx, mstar, u, nmax, lwrkfl, wrk, lwrkin, iwrk, &
         fsub_twp, dfsub_twp, gsub_twp, dgsub_twp, &
         ckappa1,gamma1,ckappa,rpar,ipar, &
         iflag)
    
    IF (IFLAG .NE. 0) THEN
       CONVERGED = .FALSE.
       
       IF (guess_exists) THEN
          ! Restore old solution
          xx = oxx
          u = ou
          nmsh = onmsh
          GUESS_EXISTS = .TRUE.
       ELSE
          ! Fail
          GUESS_EXISTS = .FALSE.
       END IF
    ELSE
       CONVERGED = .TRUE.
       GUESS_EXISTS = .TRUE.
       
       oxx = xx
       ou = u
       onmsh = nmsh
    END IF
  END SUBROUTINE SOLVE

  SUBROUTINE get_vars(X, fL, dfL, fT, dfT, jL, jT)
    USE INTERPOLATE
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(in) :: X
    double precision, DIMENSION(:,:), INTENT(out) :: fL, fT, dfL, dfT, &
         jL, jT
    INTEGER :: ix, wire
    TYPE(Interpolation), SAVE :: u_interp
         
    double precision, DIMENSION(MSTAR) :: z

    INTEGER :: jjfL, jjfT, jjdfL, jjdfT

    double precision, DIMENSION(NWIRE) :: &
         DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT

    CALL interpolation_new(u_interp, MSTAR, nmsh, xx, u)

    DO ix = 1, SIZE(X)
       z = interpolation_get(u_interp, X(ix))
       CALL GET_KINETIC(X(ix),DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)
       DO wire = 1, nwire
          CALL get_vars_from_Z(X(ix), z, wire, &
               fL(ix,wire), dfL(ix,wire), fT(ix,wire), dfT(ix,wire))

          ! Scale lengths properly
          dfL(ix,wire) = dfL(ix,wire) / WIRE_LENGTH(wire)
          dfT(ix,wire) = dfT(ix,wire) / WIRE_LENGTH(wire)

          ! Evaluate currents
          jS(wire) = jS(wire) / WIRE_LENGTH(wire)
          jT(ix,wire) = DT(wire)*dfT(ix,wire) + TT(wire)*dfL(ix,wire) + jS(wire)*fL(ix,wire)
          jL(ix,wire) = DL(wire)*dfL(ix,wire) - TT(wire)*dfT(ix,wire) + jS(wire)*fT(ix,wire)
       END DO
    END DO

  END SUBROUTINE get_vars


  !! TWPBVP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  SUBROUTINE FSUB_TWP(n, x, u, f,RPAR,IPAR)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    double precision, INTENT(IN) :: x
    double precision, DIMENSION(n), INTENT(IN) :: u
    double precision, DIMENSION(n), INTENT(OUT) :: f
    double precision, dimension(*) :: rpar
    integer, dimension(*) :: ipar

    double precision, DIMENSION(ncomp) :: ff

    INTEGER :: i

    CALL FSUB(x, u, ff)

    ! Rearrange
    DO i = 1, ncomp
       f(2*(i-1)+1) = u(2*(i-1)+2)
       f(2*(i-1)+2) = ff(i)
    END DO
  END SUBROUTINE FSUB_TWP

  SUBROUTINE DFSUB_TWP(n, x, u, df,RPAR,IPAR)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    double precision, INTENT(IN) :: x
    double precision, DIMENSION(n), INTENT(IN) :: u
    double precision, DIMENSION(n,n), INTENT(OUT) :: df
    double precision, dimension(*) :: rpar
    integer, dimension(*) :: ipar

    double precision, DIMENSION(ncomp, mstar) :: dff

    INTEGER :: i

    CALL DFSUB(x, u, dff(1:ncomp, :))

    df = 0

    ! Rearrange
    DO i = 1, ncomp
       df(2*(i-1)+1, 2*(i-1)+2) = 1
       df(2*(i-1)+2, :) = dff(i, :)
    END DO
  END SUBROUTINE DFSUB_TWP


  SUBROUTINE GSUB_TWP(i, n, u, g,RPAR,IPAR)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, n
    double precision, DIMENSION(n), INTENT(in) :: u
    double precision, INTENT(out) :: g
    double precision, dimension(*) :: rpar
    integer, dimension(*) :: ipar

    CALL GSUB(i, u, g)
  END SUBROUTINE GSUB_TWP


  SUBROUTINE DGSUB_TWP(i, n, u, dg,RPAR,IPAR)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, n
    double precision, DIMENSION(n), INTENT(IN) :: u
    double precision, DIMENSION(n), INTENT(OUT) :: dg
    double precision, dimension(*) :: rpar
    integer, dimension(*) :: ipar

    CALL DGSUB(i, u, dg)
  END SUBROUTINE DGSUB_TWP


  !! Guess !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GUESS_TWP()
    IMPLICIT NONE

    double precision, DIMENSION(mstar) :: z
    double precision, DIMENSION(ncomp) :: dm

    double precision, DIMENSION(5), PARAMETER :: x0 = &
         (/ 0d0, .25d0, .5d0, .75d0, 1d0 /)

    INTEGER :: ix

    nmsh = SIZE(x0)

    DO ix = 1, SIZE(x0)
       xx(ix) = x0(ix)
       
       CALL GUESS(x0(ix), z, dm)
       u(:, ix) = z
    END DO
  END SUBROUTINE GUESS_TWP

END MODULE KIN_SOLVE2
