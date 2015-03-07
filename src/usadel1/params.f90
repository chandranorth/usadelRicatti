! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
MODULE PARAMS
  USE MISCMATH
  USE INTERPOLATE
  use lazy_alloc


!!!
!!! Evaluation parameters
!!!
  double complex, DIMENSION(:), POINTER :: ENERGY
  INTEGER :: NENERGY
  double complex :: CURRENT_ENERGY

!!!
!!! Wire parameters
!!!
  INTEGER :: NWIRE, NTERM

  INTEGER, PARAMETER :: WIRE_TYPE_N        = 0, &
                        WIRE_TYPE_S        = 1, &
                        WIRE_TYPE_NULL     = 2

  ! Wire lengths
  double precision, DIMENSION(:), POINTER :: WIRE_LENGTH
  ! Relative (cross_section*normal_state_conductivity) of wires
  double precision, DIMENSION(:), POINTER :: WIRE_CONDUCT

  double precision, DIMENSION(:), POINTER :: WIRE_RATE_INELASTIC
  double precision, DIMENSION(:), POINTER :: WIRE_RATE_SPINFLIP
  double precision, DIMENSION(:), POINTER :: WIRE_PHASE_JUMP

  TYPE(Interpolation), SAVE, PRIVATE :: WIRE_DELTAPHASE

!!!
!!! Kinetic wire parameters
!!!
  TYPE(Interpolation), DIMENSION(:), allocatable :: WIRE_coefs

  INTEGER, PARAMETER, PRIVATE :: i_DL = 1
  INTEGER, PARAMETER, PRIVATE :: i_DT = 2
  INTEGER, PARAMETER, PRIVATE :: i_TT = 3
  INTEGER, PARAMETER, PRIVATE :: i_jS = 4
  INTEGER, PARAMETER, PRIVATE :: i_cTL = 5
  INTEGER, PARAMETER, PRIVATE :: i_cTT = 6

  INTEGER, PARAMETER, PRIVATE :: i_dDL = 7
  INTEGER, PARAMETER, PRIVATE :: i_dDT = 8
  INTEGER, PARAMETER, PRIVATE :: i_dTT = 9
  INTEGER, PARAMETER, PRIVATE :: i_djS = 10

!!!
!!! Reservoir parameters.
!!!
  double precision, DIMENSION(:), POINTER :: TERM_DELTA ! Energy gap
  double precision, DIMENSION(:), POINTER :: TERM_PHASE ! Phase

  double precision, DIMENSION(:), POINTER :: TERM_RATE_INELASTIC
  double precision, DIMENSION(:), POINTER :: TERM_RATE_SPINFLIP

  double precision, DIMENSION(:), POINTER :: TERM_T  ! Temperature
  double precision, DIMENSION(:), POINTER :: TERM_MU ! Potential
  

!!!
!!! Boundary condition parameters
!!!

  INTEGER, PARAMETER :: BCTYPE_CLEAN_N_TERMINAL         = 1, &
                        BCTYPE_CLEAN_S_TERMINAL         = 2, &
                        BCTYPE_CLEAN_NODE               = 3, &
                        BCTYPE_CLEAN_S_TERMINAL_INFTY   = 4, &
                        BCTYPE_FREE_INTERFACE           = 5, &
                        BCTYPE_FREE_INTERFACE_FIX_PHASE = 6, &
                        BCTYPE_TUNNEL_N_TERMINAL        = 7, &
                        BCTYPE_TUNNEL_S_TERMINAL        = 8, &
                        BCTYPE_TUNNEL_NODE              = 9

  INTEGER, PARAMETER :: BCSUB_CHARACTERISTIC_ZERO      = 1, &
                        BCSUB_CHARACTERISTIC_L         = 2, &
                        BCSUB_CHARACTERISTIC_T         = 3

!!!
!!! Solver parameters
!!!
  INTEGER :: SP_SOLVER_TYPE = 1, KIN_SOLVER_TYPE = 2
  double precision :: SP_SOLVER_TOL = 0, KIN_SOLVER_TOL = 0

CONTAINS

  SUBROUTINE RESIZE_PARAMS(NW, NT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NW, NT

    NWIRE = NW
    NTERM = NT

    call lazy_alloc_r1(WIRE_LENGTH, NW)
    call lazy_alloc_r1(WIRE_CONDUCT, NW)
    call lazy_alloc_r1(WIRE_RATE_INELASTIC, NW)
    call lazy_alloc_r1(WIRE_RATE_SPINFLIP, NW)
    call lazy_alloc_r1(WIRE_PHASE_JUMP, NW)
    call lazy_alloc_r1(TERM_DELTA, NT)
    call lazy_alloc_r1(TERM_PHASE, NT)
    call lazy_alloc_r1(TERM_T, NT)
    call lazy_alloc_r1(TERM_MU, NT)
    call lazy_alloc_r1(TERM_RATE_INELASTIC, NT)
    call lazy_alloc_r1(TERM_RATE_SPINFLIP, NT)
  END SUBROUTINE RESIZE_PARAMS

  ! Set Delta
  SUBROUTINE SET_DELTA(NX, X, Delta, phase)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NX
    double precision, DIMENSION(NX, NWIRE), INTENT(IN) :: X, Delta, phase

    double precision, DIMENSION(2, NWIRE, NX) :: z

    INTEGER :: k, m

    z(1, :, :) = TRANSPOSE(Delta(:, :))
    z(2, :, :) = TRANSPOSE(Phase(:, :))

    CALL interpolation_new(WIRE_DELTAPHASE, 2*NWIRE, NX, X, &
         RESHAPE(z, (/ 2*NWIRE, NX /)))
  END SUBROUTINE SET_DELTA

  SUBROUTINE SET_ENERGY(NE, E)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NE
    double complex, DIMENSION(NE), INTENT(IN) :: E

    call lazy_alloc_c1(ENERGY, NE)
    ENERGY = E
  END SUBROUTINE SET_ENERGY

  ! Set kinetic parameters
  SUBROUTINE SET_KINETIC(NE, NX, X, DL, DT, TT, jS, &
       dDL, dDT, dTT, djS, cTL, cTT)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NE, NX
    double precision, DIMENSION(NX) :: X
    double precision, DIMENSION(NX, NWIRE, NE), INTENT(IN) :: &
         DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT
    double precision, DIMENSION(10, NWIRE, NX) :: z

    INTEGER :: k, m, n

    IF (NX .LE. 0 .OR. NE .LE. 0) RETURN

    IF (.NOT. ALLOCATED(WIRE_coefs) .OR. NE /= NENERGY) THEN
       IF (ALLOCATED(WIRE_coefs)) THEN
          DO m = 1, SIZE(WIRE_coefs)
             CALL interpolation_destroy(WIRE_coefs(m))
          END DO
          DEALLOCATE(WIRE_coefs)
       END IF
       NENERGY = NE
       ALLOCATE(WIRE_coefs(NENERGY+1))
    END IF
         
    DO m = 1, NE
       DO k = 1, NWIRE
          z(i_DL,k,:) = DL(:,k,m)
          z(i_DT,k,:) = DT(:,k,m)
          z(i_TT,k,:) = TT(:,k,m)
          z(i_jS,k,:) = jS(:,k,m)*WIRE_LENGTH(k)
          z(i_dDL,k,:) = dDL(:,k,m)*WIRE_LENGTH(k)
          z(i_dDT,k,:) = dDT(:,k,m)*WIRE_LENGTH(k)
          z(i_dTT,k,:) = dTT(:,k,m)*WIRE_LENGTH(k)
          z(i_djS,k,:) = djS(:,k,m)*WIRE_LENGTH(k)*WIRE_LENGTH(k)
          z(i_cTT,k,:) = cTT(:,k,m)*WIRE_LENGTH(k)*WIRE_LENGTH(k)
          z(i_cTL,k,:) = cTL(:,k,m)*WIRE_LENGTH(k)*WIRE_LENGTH(k)
       END DO

       CALL interpolation_new(WIRE_coefs(m), 10*NWIRE, NX, X, &
            RESHAPE(z, (/ 10*NWIRE, NX /)))
    END DO
    
    ! (re)allocate some spare space
    CALL interpolation_new(WIRE_coefs(NENERGY+1), 10*NWIRE, NX, X, &
         RESHAPE(z, (/ 10*NWIRE, NX /)))
  END SUBROUTINE SET_KINETIC

  SUBROUTINE INTERPOLATE_KINETIC_FOR_ENERGY(E)
    IMPLICIT NONE
    double precision, INTENT(IN) :: E
    double precision :: x
    INTEGER :: I, k, j

    CALL lower_bound(NENERGY, real(ENERGY), E, I)
    IF (I <= 0) THEN
       I = 1
    ELSE IF (E == ENERGY(NENERGY)) THEN
       I = NENERGY-1
    ELSE IF (I >= NENERGY) THEN
       I = NENERGY-1
    END IF

    x = (E - ENERGY(I)) / (ENERGY(I+1) - ENERGY(I))

    CALL interpolation_linear_sum(WIRE_coefs(I),WIRE_coefs(I+1), &
         WIRE_coefs(NENERGY+1), x)
  END SUBROUTINE INTERPOLATE_KINETIC_FOR_ENERGY

  ! Interpolate Delta
  SUBROUTINE GET_DELTA(X, Delta, phase)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(NWIRE), INTENT(OUT) :: Delta, phase

    double precision, DIMENSION(2, NWIRE) :: z

    IF (X .GE. 0 .AND. X .LE. 1) THEN
       z = RESHAPE(interpolation_get(WIRE_DELTAPHASE, X), &
                   (/ 2, NWIRE /))
       Delta = z(1,:)
       phase = z(2,:)
    ELSE
       Delta = 0
       phase = 0
    END IF

  END SUBROUTINE GET_DELTA

  ! Interpolate kinetic coefficients
  SUBROUTINE GET_KINETIC(X, DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(NWIRE), INTENT(OUT) :: &
         DL, DT, TT, jS, dDL, dDT, dTT, djS, cTL, cTT

    double precision, DIMENSION(10, NWIRE) :: z

    IF (X .GE. 0d0 .AND. X .LE. 1d0) THEN
       z = RESHAPE(interpolation_get(WIRE_coefs(NENERGY+1), X), &
                   (/ 10, NWIRE /))
       DL  = z(i_DL,:)
       DT  = z(i_DT,:)
       TT  = z(i_TT,:)
       jS  = z(i_jS,:)
       dDL = z(i_dDL,:)
       dDT = z(i_dDT,:)
       dTT = z(i_dTT,:)
       djS = z(i_djS,:)
       cTL = z(i_cTL,:)
       cTT = z(i_cTT,:)
    ELSE
       DL = 1
       DT = 1
       TT = 0
       jS = 0
       dDL = 0
       dDT = 0
       dTT = 0
       djS = 0
       cTL = 0
       cTT = 0
    END IF
  END SUBROUTINE GET_KINETIC

  !! Bulk equilibrium fT
  double precision FUNCTION bulk_fT(E, k)
    IMPLICIT NONE
    double precision, INTENT(in) :: E
    INTEGER, INTENT(IN) :: k
    bulk_fT = .5 * (TANH((E + TERM_MU(k))*.5/TERM_T(k)) &
                  - TANH((E - TERM_MU(k))*.5/TERM_T(k)))
  END FUNCTION bulk_fT

  !! Bulk equilibrium fL
  double precision FUNCTION bulk_fL(E, k)
    IMPLICIT NONE
    double precision, INTENT(in) :: E
    INTEGER, INTENT(IN) :: k
    bulk_fL = .5 * (TANH((E + TERM_MU(k))*.5/TERM_T(k)) &
                  + TANH((E - TERM_MU(k))*.5/TERM_T(k)))
  END FUNCTION bulk_fL

END MODULE PARAMS
