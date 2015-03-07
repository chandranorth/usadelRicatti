! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
MODULE SP_SOLVE
  USE PARAMS
  USE SP_EQUATIONS
  use lazy_alloc

  PUBLIC SOLVE, INIT_SOLVER, GET_VARS

  PRIVATE

  INTEGER :: K, NTOL, NISPACE, NFSPACE

  INTEGER, DIMENSION(:), POINTER :: M, ISPACE, OISPACE, LTOL
  double precision, DIMENSION(:), POINTER :: FSPACE, OFSPACE, ZETA, TOL

  LOGICAL :: guess_exists = .FALSE.

CONTAINS

  SUBROUTINE INIT_SOLVER()
    IMPLICIT NONE

    INTEGER :: KD, KDM, NSIZEI, NSIZEF, NMAX, NREC
    INTEGER :: I

    NMAX = 500
    K = 0     ! Number of collocation points

    ! Number of right boundary conditions
    NREC = NCOMP

    ! Degree of equations
    call lazy_alloc_i1(M, NCOMP)
    M = 2

    ! Boundary condition positions
    call lazy_alloc_r1(ZETA, MSTAR)
    ZETA = 1d0
    ZETA(1:NCOMP) = 0d0

    ! Tolerances (do not specify)
    NTOL = MSTAR
    call lazy_alloc_i1(LTOL, NTOL) 
    call lazy_alloc_r1(TOL, NTOL)
    DO I = 1, MSTAR
       LTOL(I) = I
    END DO
    IF (SP_SOLVER_TOL .LE. 0) THEN
       TOL = 1e-4
    ELSE
       TOL = SP_SOLVER_TOL
    END IF

    !CALL SET_TOL(TOL, 1d-4, 1d-2, 1d-3)
    
    ! Workspace
    KD = K * NCOMP
    KDM = KD + MSTAR
    NSIZEI = 3 + KDM
    NISPACE = NMAX * NSIZEI
    call lazy_alloc_i1(ISPACE, NISPACE)
    call lazy_alloc_i1(OISPACE, NISPACE)

    NSIZEF = 8 + 4*MSTAR + (5+KD)*KDM + (2*MSTAR-NREC)*2*MSTAR + KD
    NFSPACE = NMAX * NSIZEF
    call lazy_alloc_r1(FSPACE, NFSPACE)
    call lazy_alloc_r1(OFSPACE, NFSPACE)
    
    guess_exists = .FALSE.
  END SUBROUTINE INIT_SOLVER

  SUBROUTINE SOLVE(CONVERGED, RESET)
    USE COLNEW_MOD
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: CONVERGED
    LOGICAL, OPTIONAL, INTENT(IN) :: RESET

    double precision :: ALEFT, ARIGHT
    double precision :: FIXPNT(1)
    INTEGER, DIMENSION(13) :: IPAR
    INTEGER :: IFLAG, TRY

    IF (PRESENT(RESET) .AND. RESET) THEN
       guess_exists = .false.
    END IF

    ALEFT = 0d0
    ARIGHT = 1d0

    FIXPNT = 0

    IPAR(1) = 1
    IPAR(2) = K
    IPAR(3) = 0
    IPAR(4) = NTOL
    IPAR(5) = NFSPACE
    IPAR(6) = NISPACE
    IPAR(7) = 1 ! iprint
    IPAR(8) = 0
    IPAR(9) = 0
    IPAR(10) = 0
    IPAR(11) = 0
    IPAR(12) = 0
    IPAR(13) = 0

    E = CURRENT_ENERGY

    IF (guess_exists) THEN
       IPAR(9) = 3
       IPAR(3) = ISPACE(1)
       IPAR(8) = 1
    ELSE
       IPAR(9) = 1
    END IF

    DO TRY = 1, 3
       IF (TRY == 2) THEN
          ! Try changing number of collocation points
          IPAR(2) = 4
          ISPACE = OISPACE
          FSPACE = OFSPACE
       ELSE IF (TRY == 3) THEN
          ! Restart from Ansatz
          IPAR(9) = 1
       END IF

       !! Solve the actual problem
       CALL COLNEW(NCOMP, M, ALEFT, ARIGHT, ZETA, IPAR, LTOL, &
            TOL, FIXPNT, ISPACE, FSPACE, IFLAG, &
            FSUB, DFSUB, GSUB, DGSUB, GUESS)

       IF (IFLAG .EQ. 1) EXIT
    END DO

    IF (IFLAG .NE. 1) THEN
       CONVERGED = .FALSE.
          
       IF (guess_exists) THEN
          ! Restore old solution
          ISPACE = OISPACE
          FSPACE = OFSPACE
          GUESS_EXISTS = .TRUE.
       ELSE
          ! Fail
          GUESS_EXISTS = .FALSE.
       END IF
    ELSE
       CONVERGED = .TRUE.
       GUESS_EXISTS = .TRUE.
       OISPACE = ISPACE
       OFSPACE = FSPACE
       last_E = E
    END IF
  END SUBROUTINE SOLVE

  SUBROUTINE get_vars(X, a, da, b, db)
    USE COLNEW_MOD
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: X
    double complex, DIMENSION(:, :), INTENT(out) :: &
         a, da, b, db
    
    double precision, DIMENSION(MSTAR) :: Z
    INTEGER :: wire, ix

    DO ix = 1, SIZE(X)
       CALL APPSLN(X(ix), Z, FSPACE, ISPACE)
       DO wire = 1, NWIRE
          CALL get_vars_from_Z(X(ix), Z, wire, &
                               a(ix,wire), da(ix,wire), &
                               b(ix,wire), db(ix,wire))
       END DO
    END DO
  END SUBROUTINE get_vars
END MODULE SP_SOLVE
