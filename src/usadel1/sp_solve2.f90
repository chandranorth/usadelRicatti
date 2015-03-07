! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
MODULE SP_SOLVE2
  USE PARAMS
  USE SP_EQUATIONS
  use lazy_alloc

  PUBLIC SOLVE, INIT_SOLVER, GET_VARS, SET_VARS

  PRIVATE

  INTEGER :: ntol
  INTEGER, DIMENSION(:), POINTER :: ltol
  double precision, DIMENSION(:), POINTER :: tol

  INTEGER, DIMENSION(:), POINTER :: iwrk
  double precision, DIMENSION(:), POINTER :: wrk, xx, oxx, ooxx
  double precision, DIMENSION(:,:), POINTER :: u, ou, oou
  INTEGER :: nmsh, onmsh, oonmsh

  double complex :: llast_E

  INTEGER :: nmax, nucol, nlbc, lwrkfl, lwrkin

  LOGICAL :: guess_exists = .FALSE.
  LOGICAL :: predictor_exists = .FALSE.

CONTAINS

  SUBROUTINE INIT_SOLVER()
    IMPLICIT NONE

    INTEGER :: KD, KDM, NSIZEI, NSIZEF, NMAX
    INTEGER :: I

    nmax = 200
    nucol = 2*nmax

    ! Number of left boundary conditions
    nlbc = NCOMP
    
    ! Tolerances (do not specify)
    NTOL = MSTAR
    call lazy_alloc_i1(LTOL, NTOL)
    call lazy_alloc_r1(TOL, NTOL)
    DO I = 1, MSTAR
       LTOL(I) = I
    END DO
    IF (SP_SOLVER_TOL .LE. 0) THEN
       TOL = 1e-6
    ELSE
       TOL = SP_SOLVER_TOL
    END IF

    ! Workspace
    lwrkfl = 5*mstar*mstar + 2*ntol + 9*mstar + nmax*( &
         4*mstar*mstar + 12*mstar + 3)
    lwrkin = mstar + nmax*(mstar + 2)

    call lazy_alloc_r1(xx, nucol)
    call lazy_alloc_r1(oxx, nucol)
    call lazy_alloc_r1(ooxx, nucol)
    call lazy_alloc_r2(u, mstar, nucol)
    call lazy_alloc_r2(ou, mstar, nucol)
    call lazy_alloc_r2(oou, mstar, nucol)
    call lazy_alloc_r1(wrk, lwrkfl)
    call lazy_alloc_i1(iwrk, lwrkin)
    
    guess_exists = .FALSE.
    predictor_exists = .FALSE.
  END SUBROUTINE INIT_SOLVER

  SUBROUTINE SOLVE(CONVERGED, RESET)
    USE TWPBVPC_MOD
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: CONVERGED
    LOGICAL, OPTIONAL, INTENT(IN) :: RESET

    double precision :: ALEFT, ARIGHT
    double precision :: FIXPNT(1)
    INTEGER :: IFLAG, nfxpnt, TRY
    LOGICAL :: giveeps, givemsh, giveu, linear

    integer, parameter :: liseries = 6000
    integer :: iseries(liseries)
    integer :: indnms

    double precision :: ckappa1,gamma1,ckappa, rpar(1)
    integer :: ipar(1)

    IF (PRESENT(RESET) .AND. RESET) THEN
       predictor_exists = .false.
       guess_exists = .false.
    END IF

    ALEFT = 0d0
    ARIGHT = 1d0

    FIXPNT = 0
    nfxpnt = 0

    ckappa1 = 0
    gamma1 = 0
    ckappa = 0

    linear = .FALSE.

    iprint = -1
    pdebug = .FALSE.
    uval0 = 0
    
    iwrk = 0
    wrk = 0

    E = CURRENT_ENERGY

    IF (E == 0) E = 1e-8 ! This solver fails for E=0

    IF (predictor_exists) THEN
       givemsh = .TRUE.
       giveu = .TRUE.
       
       CALL SECANT_PREDICTOR(MSTAR, 10d0, llast_E, last_E, E, &
            ooxx, oou, oonmsh, oxx, ou, onmsh, xx, u, nmsh)
       !IF (nmsh > 3*nmax/4) CALL HALVE_MESH(MSTAR, nmsh, xx, u)
    ELSE IF (guess_exists) THEN
       givemsh = .TRUE.
       giveu = .TRUE.
       !IF (nmsh > 3*nmax/4) CALL HALVE_MESH(MSTAR, nmsh, xx, u)
    ELSE
       CALL GUESS_TWP()
       givemsh = .TRUE.
       giveu = .TRUE.
    END IF

    DO TRY = 1, 2
       if (TRY == 2) then
          ! Restart from Ansatz
          CALL GUESS_TWP()
          givemsh = .TRUE.
          giveu = .TRUE.
       end if

       !CALL TWPBVP(mstar, nlbc, aleft, aright, nfxpnt, fixpnt, &
       !     ntol, ltol, tol, linear, givemsh, giveu, nmsh, &
       !     xx, mstar, u, nmax, lwrkfl, wrk, lwrkin, iwrk, &
       !     fsub_twp, dfsub_twp, gsub_twp, dgsub_twp, iflag)
       
       CALL TWPBVPC(mstar, nlbc, aleft, aright, nfxpnt, fixpnt, &
            ntol, ltol, tol, linear, givemsh, giveu, nmsh, &
            nucol, xx, mstar, u, nmax, lwrkfl, wrk, lwrkin, iwrk, &
            fsub_twp, dfsub_twp, gsub_twp, dgsub_twp, &
            ckappa1,gamma1,ckappa,rpar,ipar, &
            iflag)!,&
            !liseries, iseries, indnms)

       IF (IFLAG .EQ. 0) EXIT
    END DO

    IF (IFLAG .NE. 0) THEN
       CONVERGED = .FALSE.
          
       IF (guess_exists) THEN
          ! Restore old solution
          xx = oxx
          u = ou
          nmsh = onmsh
       ELSE
          ! Fail
          GUESS_EXISTS = .FALSE.
          predictor_exists = .false.
       END IF
    ELSE
       CONVERGED = .TRUE.

       if (guess_exists) predictor_exists = .true.
       GUESS_EXISTS = .TRUE.

       ooxx = oxx
       oou = ou
       oonmsh = onmsh
       
       oxx = xx
       ou = u
       onmsh = nmsh
       
       llast_E = last_E
       last_E = E
    END IF
  END SUBROUTINE SOLVE

  SUBROUTINE get_vars(X, a, da, b, db)
    USE INTERPOLATE
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(in) :: X
    double complex, DIMENSION(:,:), INTENT(out) :: &
         a, da, b, db
    INTEGER :: ix, wire

    double precision, DIMENSION(MSTAR) :: z

    TYPE(Interpolation), SAVE :: u_interp

    CALL interpolation_new(u_interp, MSTAR, nmsh, xx, u)

    DO ix = 1, SIZE(X)
       z = interpolation_get(u_interp, X(ix))
       DO wire = 1, nwire
          CALL get_vars_from_Z(X(ix), z, wire, &
                               a(ix,wire), da(ix,wire), &
                               b(ix,wire), db(ix,wire))
       END DO
    END DO

  END SUBROUTINE get_vars


  SUBROUTINE set_vars(X, a, da, b, db)
    USE INTERPOLATE
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(in) :: X
    double complex, DIMENSION(:,:), INTENT(in) :: a, da, b, db
    INTEGER :: ix, wire

    double precision, DIMENSION(MSTAR) :: z

    DO ix = 1, SIZE(X)
       xx(ix) = X(ix)
       DO wire = 1, nwire
          CALL set_vars_to_Z(X(ix), z, wire, &
                               a(ix,wire), da(ix,wire), &
                               b(ix,wire), db(ix,wire))
       END DO
       u(:,ix) = z
    END DO

    guess_exists = .true.
    predictor_exists = .false.
  END SUBROUTINE set_vars
  


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
    !IF (contains_nan_0(g)) WRITE(0,*) 'g>>>', g
  END SUBROUTINE GSUB_TWP

  SUBROUTINE DGSUB_TWP(i, n, u, dg,RPAR,IPAR)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, n
    double precision, DIMENSION(n), INTENT(IN) :: u
    double precision, DIMENSION(n), INTENT(OUT) :: dg
    double precision, dimension(*) :: rpar
    integer, dimension(*) :: ipar

    CALL DGSUB(i, u, dg)
    !IF (contains_nan_1(dg)) WRITE(0,*) 'dg>>>', dg
  END SUBROUTINE DGSUB_TWP


  !! Guess !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GUESS_TWP()
    IMPLICIT NONE

    double precision, DIMENSION(mstar) :: z
    double precision, DIMENSION(ncomp) :: dm

    integer, parameter :: nmsh0 = 16

    INTEGER :: ix

    nmsh = nmsh0

    DO ix = 1, nmsh0
       xx(ix) = (ix-1) * 1.0 / (nmsh0 - 1)
       
       CALL GUESS(xx(ix), z, dm)
       u(:, ix) = z
    END DO
  END SUBROUTINE GUESS_TWP

  SUBROUTINE SECANT_PREDICTOR(m, Delta, p1, p2, p3, &
       x1, u1, nmsh1, x2, u2, nmsh2, x3, u3, nmsh3)
    USE INTERPOLATE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m, nmsh1, nmsh2
    INTEGER, INTENT(INOUT) :: nmsh3
    double complex, INTENT(IN) :: p1, p2, p3
    double precision, INTENT(IN) :: Delta
    double precision, INTENT(IN) :: x1(nmsh1), x2(nmsh2)
    double precision, INTENT(IN) :: u1(m,nmsh1), u2(m,nmsh2)
    double precision, INTENT(OUT) :: x3(nmsh3), u3(m,nmsh3)

    double precision, DIMENSION(m) :: z

    TYPE(Interpolation), SAVE :: u12
    INTEGER :: k

    nmsh3 = nmsh2
    x3(1:nmsh3) = x2(1:nmsh2)

    IF (ABS(p1 - p2) > Delta) THEN
       u3(:, 1:nmsh3) = u2(:, 1:nmsh2)
    ELSE
       CALL interpolation_new(u12, m, nmsh1, x1, u1)
       DO k = 1, nmsh2
          z = interpolation_get(u12, x2(k))
          u3(:,k) = u2(:,k) + (p3 - p2)/(p2 - p1) * (u2(:,k) - z)
       END DO

       !if (contains_nan_2(u3)) write(0,*) 'u3>>>', u3
    END IF

  END SUBROUTINE SECANT_PREDICTOR

  SUBROUTINE HALVE_MESH(m, nmsh, x, u)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m
    INTEGER, INTENT(INOUT) :: nmsh
    double precision, DIMENSION(nmsh), INTENT(INOUT) :: x
    double precision, DIMENSION(m,nmsh), INTENT(INOUT) ::u

    INTEGER :: j

    j = nmsh/2
    WRITE(*,*) 'degrade', nmsh, '->', j
    x(1:j) = x(1:nmsh:2)
    u(:,1:j) = u(:,1:nmsh:2)
    nmsh = j
  END SUBROUTINE HALVE_MESH

END MODULE SP_SOLVE2
