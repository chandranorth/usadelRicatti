! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
MODULE KIN_SOLVE
  USE PARAMS
  USE KIN_EQUATIONS
  use lazy_alloc

  PUBLIC SOLVE, INIT_SOLVER, GET_VARS

  PRIVATE

  INTEGER :: K, NTOL, NISPACE, NFSPACE

  INTEGER, DIMENSION(:), POINTER :: M, ISPACE, OISPACE, LTOL
  double precision, DIMENSION(:), POINTER :: FSPACE, OFSPACE, TOL

  LOGICAL :: guess_exists = .FALSE.

CONTAINS

  SUBROUTINE INIT_SOLVER(NREC)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NREC

    INTEGER :: KD, KDM, NSIZEI, NSIZEF, NMAX
    INTEGER :: I, J

    NMAX = 500
    K = 2     ! Number of collocation points

    ! Degree of equations
    call lazy_alloc_i1(M, NCOMP)
    M = 2

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
    KD = K * NCOMP
    KDM = KD + MSTAR
    NSIZEI = 3 + KDM
    NISPACE = NMAX * NSIZEI
    call lazy_alloc_i1(ISPACE, NISPACE)
    call lazy_alloc_i1(OISPACE, NISPACE)

    NSIZEF = 4 + 3*MSTAR + (5+KD)*KDM + (2*MSTAR-NREC)*2*MSTAR
    NFSPACE = NMAX * NSIZEF
    call lazy_alloc_r1(FSPACE, NFSPACE)
    call lazy_alloc_r1(OFSPACE, NFSPACE)
    
    guess_exists = .FALSE.
  END SUBROUTINE INIT_SOLVER

  SUBROUTINE SOLVE(CONVERGED)
    USE COLNEW_MOD
    IMPLICIT NONE
    LOGICAL, INTENT(OUT) :: CONVERGED

    double precision :: ALEFT, ARIGHT
    double precision :: FIXPNT(1)
    INTEGER, DIMENSION(11) :: IPAR
    INTEGER :: IFLAG
    INTEGER :: TRY

    CALL INTERPOLATE_KINETIC_FOR_ENERGY(real(CURRENT_ENERGY))
    
    ALEFT = 0d0
    ARIGHT = 1d0

    FIXPNT = 0

    IPAR(1) = 0 ! This is linear!
    IPAR(2) = K
    IPAR(3) = 0
    IPAR(4) = NTOL
    IPAR(5) = NFSPACE
    IPAR(6) = NISPACE
    IPAR(7) = 1
    IPAR(8) = 0
    IPAR(9) = 0
    IPAR(10) = 0
    IPAR(11) = 0

    DO TRY = 1,2

       IF (guess_exists) THEN
          IPAR(9) = 3
          IPAR(3) = ISPACE(1)
       ELSE
          IPAR(9) = 1
       END IF
       !IPAR(9) = 1

       IF (TRY .EQ. 2) THEN
          IPAR(10) = 1
       END IF
       
       CALL COLNEW(NCOMP, M, ALEFT, ARIGHT, ZETA, IPAR, LTOL, &
                   TOL, FIXPNT, ISPACE, FSPACE, IFLAG, &
                   FSUB, DFSUB, GSUB, DGSUB, GUESS)
       
       IF (IFLAG .NE. 1) THEN
          WRITE(0,*) '%% Failed to converge: ', IFLAG, ' try ', TRY
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
          EXIT
       END IF
    END DO
  END SUBROUTINE SOLVE

  SUBROUTINE get_vars(X, fL, dfL, fT, dfT, jL, jT)
    USE COLNEW_MOD
    IMPLICIT NONE
    double precision, DIMENSION(:), INTENT(IN) :: X
    double precision, DIMENSION(:,:), INTENT(out) :: fL, dfL, fT, dfT, jL, jT
    
    double precision, DIMENSION(MSTAR) :: Z
    INTEGER :: wire, ix

    double precision, DIMENSION(NWIRE) :: &
         DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT

    DO ix = 1, SIZE(X)
       CALL APPSLN(X(ix), Z, FSPACE, ISPACE)
       CALL GET_KINETIC(X(ix),DL,DT,TT,jS,dDL,dDT,dTT,djS,cTL,cTT)

       DO wire = 1, NWIRE
          CALL get_vars_from_Z(X(ix), Z, wire, &
                               fL(ix,wire), dfL(ix,wire), &
                               fT(ix,wire), dfT(ix,wire))
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
END MODULE KIN_SOLVE
