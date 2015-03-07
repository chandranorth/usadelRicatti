! -*-f90-*-
!
! :Author: Pauli Virtanen <pauli@ltl.tkk.fi>
! :Organization: Low Temperatory Laboratory, Helsinki University of Technology
! :Date: 2005-2006
!
MODULE SP_EQUATIONS
  USE PARAMS
  USE MISCMATH
  use lazy_alloc
  IMPLICIT NONE

  INTEGER :: NCOMP, MSTAR
  INTEGER, DIMENSION(:), POINTER :: WIREOFFSET  ! [NWIRE]
  double precision, DIMENSION(:), POINTER :: lratio

  INTEGER, DIMENSION(:),   POINTER :: WIRE_TYPE   ! [NWIRE]
  INTEGER, DIMENSION(:),   POINTER :: BCTYPE      ! [MSTAR]
  double precision, DIMENSION(:,:),  POINTER :: BCDATA ! [MSTAR,16]
  INTEGER, DIMENSION(:,:), POINTER :: BCCONNECT   ! [MSTAR, NWIRE+1]

  double complex :: E, last_E
  double precision :: RELAX = 0d0

CONTAINS

  SUBROUTINE INIT_EQUATIONS(WIRE_TYPE_NEW)
    IMPLICIT NONE
    INTEGER, DIMENSION(NWIRE), INTENT(IN) :: WIRE_TYPE_NEW
    
    INTEGER :: I
    
    call lazy_alloc_i1(WIREOFFSET, NWIRE)
    call lazy_alloc_i1(WIRE_TYPE, NWIRE)

    WIRE_TYPE = WIRE_TYPE_NEW

    NCOMP = 0
    DO I = 1, NWIRE
       WIREOFFSET(I) = NCOMP
       NCOMP = NCOMP + 4
    END DO
    MSTAR = 2*NCOMP

    call lazy_alloc_i1(BCTYPE, MSTAR)
    call lazy_alloc_r2(BCDATA, MSTAR, 16)
    call lazy_alloc_i2(BCCONNECT, MSTAR, NWIRE+1)
    call lazy_alloc_r1(lratio, NWIRE)

    lratio = WIRE_LENGTH
  END SUBROUTINE INIT_EQUATIONS
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Equations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine indices(wire, &
     jddra, jddia, jddrb, jddib, jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)
  implicit none

  integer, intent(in) :: wire
  integer, intent(out) :: &
       jddra, jddia, jddrb, jddib, jra, jia, jrb, jib, jdra, jdia, jdrb, jdib

  jddra = (1 + WIREOFFSET(wire))
  jddia = (2 + WIREOFFSET(wire))
  jddrb = (3 + WIREOFFSET(wire))
  jddib = (4 + WIREOFFSET(wire))
  
  jra =  ((0 + WIREOFFSET(wire))*2 + 1)
  jia =  ((1 + WIREOFFSET(wire))*2 + 1)
  jrb =  ((2 + WIREOFFSET(wire))*2 + 1)
  jib =  ((3 + WIREOFFSET(wire))*2 + 1)
  
  jdra = ((0 + WIREOFFSET(wire))*2 + 2)
  jdia = ((1 + WIREOFFSET(wire))*2 + 2)
  jdrb = ((2 + WIREOFFSET(wire))*2 + 2)
  jdib = ((3 + WIREOFFSET(wire))*2 + 2)
end subroutine indices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Usadel equations

subroutine eval_d(a, b, da, db, Delta, phi, EE, Gsf, Gin, &
     dda, ddb)
  implicit none
  double complex, intent(in) :: a, b, da, db
  double precision, intent(in) :: Delta, phi, Gsf, Gin
  double complex, intent(in) :: EE
  double complex, intent(out) :: dda, ddb

  double complex :: ab_1

  ab_1 = a*b - 1

  dda = &
       -(0,2)*EE*a + 2*b*da**2/ab_1 &
       - 2*Gsf*a*(a*b+1)/ab_1 &
       + (0,1)*Delta*(exp(-(0,1)*phi)*a**2 + exp((0,1)*phi))
  
  ddb = &
       -(0,2)*EE*b + 2*a*db**2/ab_1 &
       - 2*Gsf*b*(a*b+1)/ab_1 &
       + (0,1)*Delta*(exp((0,1)*phi)*b**2 + exp(-(0,1)*phi))
end subroutine eval_d

!!!--------------------------------------------------------------------------
!!! Grad of Usadel equations

subroutine eval_dd(a, b, da, db, Delta, phi, EE, Gsf, Gin, &
     dda_a, dda_b, dda_da, dda_db, &
     ddb_a, ddb_b, ddb_da, ddb_db)
  implicit none
  double complex, intent(in) :: a, b, da, db
  double precision, intent(in) :: Delta, phi, Gsf, Gin
  double complex, intent(in) :: EE
  double complex, intent(out) :: &
     dda_a, dda_b, dda_da, dda_db, &
     ddb_a, ddb_b, ddb_da, ddb_db

  double complex :: ab_1

  ab_1 = a*b - 1

  dda_a = &
       - (0,2)*EE - 2*b*da**2/ab_1**2*b &
       - 2*Gsf*( (a*b+1)/ab_1 + a*b/ab_1 - a*(a*b+1)/ab_1**2*b ) &
       + (0,1)*Delta*(exp(-(0,1)*phi)*2*a)
  
  dda_b = &
       + 2*da**2/ab_1 - 2*b*da**2/ab_1**2*a &
       - 2*Gsf*( a*a/ab_1 - a*(a*b+1)/ab_1**2*a )
  
  dda_da = &
       2*b*2*da/ab_1
  
  dda_db = 0
       
  !!---
       
  ddb_a = &
       + 2*db**2/ab_1 - 2*a*db**2/ab_1**2*b &
       - 2*Gsf*( b*b/ab_1 - b*(a*b+1)/ab_1**2*b )
  
  ddb_b = &
       - (0,2)*EE  - 2*a*db**2/ab_1**2*a &
       - 2*Gsf*( (a*b+1)/ab_1 + b*(a)/ab_1 - b*(a*b+1)/ab_1**2*a ) &
       + (0,1)*Delta*(exp(+(0,1)*phi)*2*b)
  
  ddb_da = 0
  ddb_db = 2*a*2*db/ab_1
end subroutine eval_dd

!!!--------------------------------------------------------------------------
!!! Kupriyanov-Lukichev

subroutine eval_kl(r, a1, a2, b2, da1, bc1)
  implicit none
  double precision, intent(in) :: r
  double complex, intent(in) :: a1, a2, b2, da1
  double complex, intent(out) :: bc1

  bc1 = (-a1 + a2)*(-1 + a1*b2) - (-1 + a2*b2)*da1*r
end subroutine eval_kl

subroutine eval_kl_d(r, a1, a2, b2, da1, bc1_a1, bc1_a2, bc1_b2, bc1_da1)
  implicit none
  double precision, intent(in) :: r
  double complex, intent(in) :: a1, a2, b2, da1
  double complex, intent(out) :: bc1_a1, bc1_a2, bc1_b2, bc1_da1

  bc1_a1 = 1 - 2*a1*b2 + a2*b2
  bc1_da1 = r*(1 - a2*b2)

  bc1_a2 = -1 + a1*b2 - b2*da1*r
  bc1_b2 = (-a1 + a2)*a1 - a2*da1*r
end subroutine eval_kl_d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE FSUB(X, Z, F)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    double precision, DIMENSION(NCOMP), INTENT(OUT) :: F

    INTEGER :: wire, m, k
    double precision :: Gsf, Gin
    double precision, DIMENSION(NWIRE) :: Delta, phase
    double precision :: y

    double complex :: a, b, da, db, EE
    double complex :: dda, ddb

    integer :: jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib

    F = 0

    CALL GET_DELTA(X, Delta, phase)

    DO wire=1,NWIRE
       IF (WIRE_TYPE(wire) .EQ. WIRE_TYPE_NULL) CYCLE

       call indices(wire, jddra, jddia, jddrb, jddib, &
            jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)

       Gin = WIRE_RATE_INELASTIC(wire)
       Gsf = WIRE_RATE_SPINFLIP(wire)
       EE = E + imagu*Gin

       EE = EE * lratio(wire)**2
       Gsf = Gsf * lratio(wire)**2
       Delta(wire) = Delta(wire)  * lratio(wire)**2

       a = Z(jra) + imagu*Z(jia)
       da = Z(jdra) + imagu*Z(jdia)

       b = Z(jrb) + imagu*Z(jib)
       db = Z(jdrb) + imagu*Z(jdib)

       call eval_d(a, b, da, db, Delta(wire), phase(wire), EE, Gsf, Gin, &
            dda, ddb)

       !!! Assign
       F(jddra) = REAL (dda)
       F(jddia) = AIMAG(dda)
       
       F(jddrb) = REAL (ddb)
       F(jddib) = AIMAG(ddb)
    END DO
  END SUBROUTINE FSUB

  SUBROUTINE DFSUB(X, Z, DF)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    double precision, DIMENSION(NCOMP, MSTAR), INTENT(OUT) :: DF

    INTEGER :: wire, m, n, k
    double precision :: Gsf, Gin
    double precision, DIMENSION(NWIRE) :: Delta, phase
    double precision :: y

    double complex :: a, b, da, db, EE
    double complex :: dda_a, dda_b, ddb_a, ddb_b, &
         dda_da, dda_db, ddb_da, ddb_db
    integer :: jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib


    DF = 0

    CALL GET_DELTA(X, Delta, phase)

    DO wire=1,NWIRE
       IF (WIRE_TYPE(wire) .EQ. WIRE_TYPE_NULL) CYCLE

       call indices(wire, jddra, jddia, jddrb, jddib, &
            jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)

       Gin = WIRE_RATE_INELASTIC(wire)
       Gsf = WIRE_RATE_SPINFLIP(wire)
       EE = E + imagu*Gin

       EE = EE * lratio(wire)**2
       Gsf = Gsf * lratio(wire)**2
       Delta(wire) = Delta(wire)  * lratio(wire)**2

       a = Z(jra) + imagu*Z(jia)
       da = Z(jdra) + imagu*Z(jdia)

       b = Z(jrb) + imagu*Z(jib)
       db = Z(jdrb) + imagu*Z(jdib)

       call eval_dd(a, b, da, db, Delta(wire), phase(wire), EE, Gsf, Gin, &
            dda_a, dda_b, dda_da, dda_db, &
            ddb_a, ddb_b, ddb_da, ddb_db)

       !!! Assign
       DF(jddra, jra) =  REAL (dda_a)
       DF(jddra, jia) = -AIMAG(dda_a)
       DF(jddia, jra) =  AIMAG(dda_a)
       DF(jddia, jia) =  REAL (dda_a)

       DF(jddra, jrb) =  REAL (dda_b)
       DF(jddra, jib) = -AIMAG(dda_b)
       DF(jddia, jrb) =  AIMAG(dda_b)
       DF(jddia, jib) =  REAL (dda_b)
       
       DF(jddra, jdra) =  REAL (dda_da)
       DF(jddra, jdia) = -AIMAG(dda_da)
       DF(jddia, jdra) =  AIMAG(dda_da)
       DF(jddia, jdia) =  REAL (dda_da)

       DF(jddra, jdrb) =  REAL (dda_db)
       DF(jddra, jdib) = -AIMAG(dda_db)
       DF(jddia, jdrb) =  AIMAG(dda_db)
       DF(jddia, jdib) =  REAL (dda_db)
       
       DF(jddrb, jra) =  REAL (ddb_a)
       DF(jddrb, jia) = -AIMAG(ddb_a)
       DF(jddib, jra) =  AIMAG(ddb_a)
       DF(jddib, jia) =  REAL (ddb_a)
       
       DF(jddrb, jrb) =  REAL (ddb_b)
       DF(jddrb, jib) = -AIMAG(ddb_b)
       DF(jddib, jrb) =  AIMAG(ddb_b)
       DF(jddib, jib) =  REAL (ddb_b)
       
       DF(jddrb, jdra) =  REAL (ddb_da)
       DF(jddrb, jdia) = -AIMAG(ddb_da)
       DF(jddib, jdra) =  AIMAG(ddb_da)
       DF(jddib, jdia) =  REAL (ddb_da)
       
       DF(jddrb, jdrb) =  REAL (ddb_db)
       DF(jddrb, jdib) = -AIMAG(ddb_db)
       DF(jddib, jdrb) =  AIMAG(ddb_db)
       DF(jddib, jdib) =  REAL (ddb_db)
    END DO
  END SUBROUTINE DFSUB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DOUBLE PRECISION FUNCTION g_val(I, Z, wire, id)
    IMPLICIT NONE
    integer, intent(in) :: I, wire
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    character*3, intent(in) :: id
    
    logical :: imag
    integer :: jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib
    integer :: jr, ji
    double precision :: phi

    call indices(wire, jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)

    if (id(1:2) == 'ra' .or. id(1:2) == 'ia') then
       jr = jra
       ji = jia
    else if (id(1:2) == 'rb' .or. id(1:2) == 'ib') then
       jr = jrb
       ji = jib
    else if (id == 'dra' .or. id == 'dia') then
       jr = jdra
       ji = jdia
    else if (id == 'drb' .or. id == 'dib') then
       jr = jdrb
       ji = jdib
    else
       write(0,*) 'internal error in sp_equations.g_val'
       write(0,*) id
       stop
    end if
    imag = (id(1:2) == 'ia' .or. id(1:2) == 'ib' &
         .or. id == 'dia' .or. id == 'dib')
    
    if (I <= NCOMP) then
       ! at the beginning of a wire
       if (.not. imag) then
          g_val = Z(jr)
       else
          g_val = Z(ji)
       end if
    else
       ! at the end of a wire: apply phase jump
       phi = -WIRE_PHASE_JUMP(wire)
       if (id(1:2) == 'rb' .or. id == 'drb' &
            .or. id(1:2) == 'ib' .or. id == 'dib') then
          phi = -phi
       end if
       if (.not. imag) then
          g_val = Z(jr) * cos(phi) - Z(ji) * sin(phi)
       else
          g_val = Z(jr) * sin(phi) + Z(ji) * cos(phi)
       end if
    end if
  END FUNCTION g_val

  SUBROUTINE g_der(I, DG, wire, id, val)
    IMPLICIT NONE
    integer, intent(in) :: I, wire
    double precision, DIMENSION(MSTAR), INTENT(INOUT) :: DG
    character*3, intent(in) :: id
    double precision, intent(IN) :: val
    
    logical :: imag
    integer :: jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib
    integer :: jr, ji
    double precision :: phi

    call indices(wire, jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)

    if (id(1:2) == 'ra' .or. id(1:2) == 'ia') then
       jr = jra
       ji = jia
    else if (id(1:2) == 'rb' .or. id(1:2) == 'ib') then
       jr = jrb
       ji = jib
    else if (id == 'dra' .or. id == 'dia') then
       jr = jdra
       ji = jdia
    else if (id == 'drb' .or. id == 'dib') then
       jr = jdrb
       ji = jdib
    else
       write(0,*) 'internal error in sp_equations.g_der'
       write(0,*) id
       stop
    end if
    imag = (id(1:2) == 'ia' .or. id(1:2) == 'ib' &
         .or. id == 'dia' .or. id == 'dib')
    
    if (I <= NCOMP) then
       ! at the beginning of a wire
       if (.not. imag) then
          DG(jr) = DG(jr) + val
       else
          DG(ji) = DG(ji) + val
       end if
    else
       ! at the end of a wire: apply phase jump
       phi = -WIRE_PHASE_JUMP(wire)
       if (id(1:2) == 'rb' .or. id == 'drb' &
            .or. id(1:2) == 'ib' .or. id == 'dib') then
          phi = -phi
       end if
       if (.not. imag) then
          DG(jr) = DG(jr) + val * cos(phi)
          DG(ji) = DG(ji) - val * sin(phi)
       else
          DG(jr) = DG(jr) + val * sin(phi)
          DG(ji) = DG(ji) + val * cos(phi)
       end if
    end if
  END SUBROUTINE g_der

  SUBROUTINE GSUB(I, Z, G)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    double precision, INTENT(OUT) :: G

    INTEGER :: J, BCI, BCT, BCS
    INTEGER :: wire, wire2, term
    double complex :: bc_a, bc_b, zz, aa, bb, &
         bc1, bc2
    double complex :: Gin, EE
    double precision :: rr

    G = 0

    BCS = IAND(ISHFT(BCTYPE(I), -16), 255)
    BCT = IAND(ISHFT(BCTYPE(I), -8), 255)
    BCI = IAND(BCTYPE(I), 255)

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (BCT .EQ. BCTYPE_CLEAN_S_TERMINAL) THEN
       wire = BCCONNECT(I, 1)
       term = BCCONNECT(I, 2)

       EE = E
       Gin = TERM_RATE_INELASTIC(term)

       zz = TERM_DELTA(term) / (EE + Gin*imagu)

       bc_a = exp(+(0,1)*TERM_PHASE(term)) * zz / (sqrt(1 - zz**2) + 1 + RELAX)
       bc_b = exp(-(0,1)*TERM_PHASE(term)) * zz / (sqrt(1 - zz**2) + 1 + RELAX)

       IF (BCI == 1) THEN
          G = g_val(I,Z,wire,'ra') - real(bc_a)
       ELSE IF (BCI == 2) THEN
          G = g_val(I,Z,wire,'ia') - aimag(bc_a)
       ELSE IF (BCI == 3) THEN
          G = g_val(I,Z,wire,'rb') - real(bc_b)
       ELSE IF (BCI == 4) THEN
          G = g_val(I,Z,wire,'ib') - aimag(bc_b)
       ELSE
          WRITE(*,*) '%% Invalid BC 1: ', BCI
       END IF
    ELSE IF (BCT .EQ. BCTYPE_CLEAN_N_TERMINAL) THEN
       wire = BCCONNECT(I, 1)

       IF (BCI == 1) THEN
          G = g_val(I,Z,wire,'ra') - 0
       ELSE IF (BCI == 2) THEN
          G = g_val(I,Z,wire,'ia') - 0
       ELSE IF (BCI == 3) THEN
          G = g_val(I,Z,wire,'rb') - 0
       ELSE IF (BCI == 4) THEN
          G = g_val(I,Z,wire,'ib') - 0
       ELSE
          WRITE(*,*) '%% Invalid BC 2'
       END IF
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_TUNNEL_S_TERMINAL .OR. &
         BCT .EQ. BCTYPE_TUNNEL_N_TERMINAL) THEN
       wire = BCCONNECT(I, 1)
       term = BCCONNECT(I, 2)

       EE = E
       Gin = TERM_RATE_INELASTIC(term)
       rr = BCDATA(I,1) + RELAX

       IF (I <= NCOMP) THEN
          ! At the beginning of the wire
          rr = -rr
       END IF

       rr = rr * (WIRE_CONDUCT(wire)/lratio(wire))

       zz = TERM_DELTA(term) / (EE + Gin*imagu)

       IF (BCT .EQ. BCTYPE_TUNNEL_S_TERMINAL) THEN
          bc_a = exp(+(0,1)*TERM_PHASE(term)) * zz / (sqrt(1 - zz**2) + 1 + RELAX)
          bc_b = exp(-(0,1)*TERM_PHASE(term)) * zz / (sqrt(1 - zz**2) + 1 + RELAX)
       ELSE
          bc_a = 0.
          bc_b = 0.
       END IF

       IF (BCI == 1 .or. BCI == 2) THEN
          call eval_kl(rr, &
               g_val(I,Z,wire,'ra') + imagu*g_val(I,Z,wire,'ia'), &
               bc_a, bc_b, &
               g_val(I,Z,wire,'dra') + imagu*g_val(I,Z,wire,'dia'), &
               bc1)
          
          IF (BCI == 1) THEN
             G = real(bc1)
          ELSE IF (BCI == 2) THEN
             G = aimag(bc1)
          END IF
       ELSE IF (BCI == 3 .OR. BCI == 4) THEN
          call eval_kl(rr, &
               g_val(I,Z,wire,'rb') + imagu*g_val(I,Z,wire,'ib'), &
               bc_b, bc_a, &
               g_val(I,Z,wire,'drb') + imagu*g_val(I,Z,wire,'dib'), &
               bc1)
          IF (BCI == 3) THEN
             G = real(bc1)
          ELSE IF (BCI == 4) THEN
             G = aimag(bc1)
          END IF
       ELSE
          WRITE(*,*) '%% Invalid BC tunnel: ', BCI
       END IF
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_TUNNEL_NODE) THEN
       wire = BCCONNECT(I, 1)
       wire2 = BCCONNECT(I, 2)

       IF (NWIRE .NE. 2) THEN
          WRITE(*,*) '%% Invalid number of wires (should be 2): ', BCI
          STOP
       END IF

       rr = BCDATA(I,1) + RELAX

       IF (I <= NCOMP) THEN
          ! At the beginning of the wire
          rr = -rr
       END IF

       rr = rr * (WIRE_CONDUCT(wire)/lratio(wire))

       bc_a = g_val(I,Z,wire2,'ra') + imagu*g_val(I,Z,wire2,'ia')
       bc_b = g_val(I,Z,wire2,'rb') + imagu*g_val(I,Z,wire2,'ib')

       IF (BCI == 1 .or. BCI == 2) THEN
          call eval_kl(rr, &
               g_val(I,Z,wire,'ra') + imagu*g_val(I,Z,wire,'ia'), &
               bc_a, bc_b, &
               g_val(I,Z,wire,'dra') + imagu*g_val(I,Z,wire,'dia'), &
               bc1)

          IF (BCI == 1) THEN
             G = real(bc1)
          ELSE IF (BCI == 2) THEN
             G = aimag(bc1)
          END IF
       ELSE IF (BCI == 3 .OR. BCI == 4) THEN
          call eval_kl(rr, &
               g_val(I,Z,wire,'rb') + imagu*g_val(I,Z,wire,'ib'), &
               bc_b, bc_a, &
               g_val(I,Z,wire,'drb') + imagu*g_val(I,Z,wire,'dib'), &
               bc1)
          IF (BCI == 3) THEN
             G = real(bc1)
          ELSE IF (BCI == 4) THEN
             G = aimag(bc1)
          END IF
       ELSE
          WRITE(*,*) '%% Invalid BC tunnel: ', BCI
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_CLEAN_NODE) THEN
       wire = BCCONNECT(I, 1)
       wire2 = BCCONNECT(I, 2)

       IF (BCI == 1) THEN
          G = g_val(I,Z,wire,'ra') - g_val(I,Z,wire2,'ra')
       ELSE IF (BCI == 2) THEN
          G = g_val(I,Z,wire,'ia') - g_val(I,Z,wire2,'ia')
       ELSE IF (BCI == 3) THEN
          G = g_val(I,Z,wire,'rb') - g_val(I,Z,wire2,'rb')
       ELSE IF (BCI == 4) THEN
          G = g_val(I,Z,wire,'ib') - g_val(I,Z,wire2,'ib')
       ELSE IF (BCI == 5) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * g_val(I,Z,wire,'dra')
          END DO
       ELSE IF (BCI == 6) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * g_val(I,Z,wire,'dia')
          END DO
       ELSE IF (BCI == 7) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * g_val(I,Z,wire,'drb')
          END DO
       ELSE IF (BCI == 8) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * g_val(I,Z,wire,'dib')
          END DO
       ELSE
          WRITE(*,*) '%% Invalid BC 3'
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_FREE_INTERFACE) THEN
       wire = BCCONNECT(I, 1)
       
       IF (BCI == 1) THEN
          G = g_val(I,Z,wire,'dra') - 0
       ELSE IF (BCI == 2) THEN
          G = g_val(I,Z,wire,'dia') - 0
       ELSE IF (BCI == 3) THEN
          G = g_val(I,Z,wire,'drb') - 0
       ELSE IF (BCI == 4) THEN
          G = g_val(I,Z,wire,'dib') - 0
       ELSE
          WRITE(*,*) '%% Invalid BC 5'
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_FREE_INTERFACE_FIX_PHASE) THEN
       wire = BCCONNECT(I, 1)
       
       IF (BCI == 1) THEN
          G = g_val(I,Z,wire,'ra') - g_val(I,Z,wire,'rb')
       ELSE IF (BCI == 2) THEN
          G = g_val(I,Z,wire,'ia') - g_val(I,Z,wire,'ib')
       ELSE IF (BCI == 3) THEN
          G = g_val(I,Z,wire,'drb') - 0
       ELSE IF (BCI == 4) THEN
          G = g_val(I,Z,wire,'dib') - 0
       ELSE
          WRITE(*,*) '%% Invalid BC 6'
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_CLEAN_NODE) THEN
       wire = BCCONNECT(I, 1)
       wire2 = BCCONNECT(I, 2)

       IF (BCI == 1) THEN
          G = g_val(I,Z,wire,'ra') - g_val(I,Z,wire2,'ra')
       ELSE IF (BCI == 2) THEN
          G = g_val(I,Z,wire,'ia') - g_val(I,Z,wire2,'ia')
       ELSE IF (BCI == 3) THEN
          G = g_val(I,Z,wire,'rb') - g_val(I,Z,wire2,'rb')
       ELSE IF (BCI == 4) THEN
          G = g_val(I,Z,wire,'ib') - g_val(I,Z,wire2,'ib')
       ELSE IF (BCI == 5) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * g_val(I,Z,wire,'dra')
          END DO
       ELSE IF (BCI == 6) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * g_val(I,Z,wire,'dia')
          END DO
       ELSE IF (BCI == 7) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * g_val(I,Z,wire,'drb')
          END DO
       ELSE IF (BCI == 8) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             G = G + WIRE_CONDUCT(wire)/lratio(wire) * g_val(I,Z,wire,'dib')
          END DO
       ELSE
          WRITE(*,*) '%% Invalid BC 3'
       END IF
    ELSE
       WRITE(*,*) '%% Invalid spectral GSUB BC'
    END IF
  END SUBROUTINE GSUB

  SUBROUTINE DGSUB(I, Z, DG)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: I
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    double precision, DIMENSION(MSTAR), INTENT(OUT) :: DG

    INTEGER :: term, wire, wire2, J, BCI, BCT
    double complex :: Gin, EE, zz
    double complex :: bc1_a1, bc1_b1, bc1_da1, bc1_db1, &
         bc_a, bc_b, bc1_a2, bc1_b2
    double precision :: rr

    BCT = ISHFT(BCTYPE(I), -8)
    BCI = IAND(BCTYPE(I), 255)

    DG = 0

    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (BCT .EQ. BCTYPE_CLEAN_S_TERMINAL .OR. &
         BCT .EQ. BCTYPE_CLEAN_N_TERMINAL) THEN
       wire = BCCONNECT(I, 1)

       IF (BCI == 1) THEN
          call g_der(I,DG,wire,'ra',1d0)
       ELSE IF (BCI == 2) THEN
          call g_der(I,DG,wire,'ia',1d0)
       ELSE IF (BCI == 3) THEN
          call g_der(I,DG,wire,'rb',1d0)
       ELSE IF (BCI == 4) THEN
          call g_der(I,DG,wire,'ib',1d0)
       END IF
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_TUNNEL_S_TERMINAL .OR. &
         BCT .EQ. BCTYPE_TUNNEL_N_TERMINAL) THEN
       wire = BCCONNECT(I, 1)
       term = BCCONNECT(I, 2)

       EE = E
       Gin = TERM_RATE_INELASTIC(term)
       rr = BCDATA(I,1) + RELAX

       IF (I <= NCOMP) THEN
          ! At the beginning of the wire
          rr = -rr
       END IF
       
       rr = rr * (WIRE_CONDUCT(wire)/lratio(wire))

       zz = TERM_DELTA(term) / (EE + Gin*imagu)

       IF (BCT .EQ. BCTYPE_TUNNEL_S_TERMINAL) THEN
          bc_a = exp(+(0,1)*TERM_PHASE(term)) * zz / (sqrt(1 - zz**2) + 1 + RELAX)
          bc_b = exp(-(0,1)*TERM_PHASE(term)) * zz / (sqrt(1 - zz**2) + 1 + RELAX)
       ELSE
          bc_a = 0.
          bc_b = 0.
       END IF

       IF (BCI == 1 .or. BCI == 2) THEN
          call eval_kl_d(rr, &
               g_val(I,Z,wire,'ra') + imagu*g_val(I,Z,wire,'ia'), &
               bc_a, bc_b, &
               g_val(I,Z,wire,'dra') + imagu*g_val(I,Z,wire,'dia'), &
               bc1_a1, bc1_a2, bc1_b2, bc1_da1)

          IF (BCI == 1) THEN
             call g_der(I,DG,wire,'ra',   real(bc1_a1))
             call g_der(I,DG,wire,'ia', -aimag(bc1_a1))
             call g_der(I,DG,wire,'dra',  real(bc1_da1))
             call g_der(I,DG,wire,'dia',-aimag(bc1_da1))
          ELSE IF (BCI == 2) THEN
             call g_der(I,DG,wire,'ra',  aimag(bc1_a1))
             call g_der(I,DG,wire,'ia',   real(bc1_a1))
             call g_der(I,DG,wire,'dra', aimag(bc1_da1))
             call g_der(I,DG,wire,'dia',  real(bc1_da1))
          END IF
       ELSE IF (BCI == 3 .OR. BCI == 4) THEN
          call eval_kl_d(rr, &
               g_val(I,Z,wire,'rb') + imagu*g_val(I,Z,wire,'ib'), &
               bc_b, bc_a, &
               g_val(I,Z,wire,'drb') + imagu*g_val(I,Z,wire,'dib'), &
               bc1_b1, bc1_b2, bc1_a2, bc1_db1)

          IF (BCI == 3) THEN
             call g_der(I,DG,wire,'rb',   real(bc1_b1))
             call g_der(I,DG,wire,'ib', -aimag(bc1_b1))
             call g_der(I,DG,wire,'drb',  real(bc1_db1))
             call g_der(I,DG,wire,'dib',-aimag(bc1_db1))
          ELSE IF (BCI == 4) THEN
             call g_der(I,DG,wire,'rb',  aimag(bc1_b1))
             call g_der(I,DG,wire,'ib',   real(bc1_b1))
             call g_der(I,DG,wire,'drb', aimag(bc1_db1))
             call g_der(I,DG,wire,'dib',  real(bc1_db1))
          END IF
       ELSE
          WRITE(*,*) '%% Invalid BC 1: ', BCI
       END IF
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_TUNNEL_NODE) THEN
       wire = BCCONNECT(I, 1)
       wire2 = BCCONNECT(I, 2)

       rr = BCDATA(I,1) + RELAX

       IF (I <= NCOMP) THEN
          ! At the beginning of the wire
          rr = -rr
       END IF
       
       rr = rr * (WIRE_CONDUCT(wire)/lratio(wire))

       IF (BCI == 1 .or. BCI == 2) THEN
          call eval_kl_d(rr, &
               g_val(I,Z,wire,'ra') + imagu*g_val(I,Z,wire,'ia'), &
               g_val(I,Z,wire2,'ra') + imagu*g_val(I,Z,wire2,'ia'), &
               g_val(I,Z,wire2,'rb') + imagu*g_val(I,Z,wire2,'ib'), &
               g_val(I,Z,wire,'dra') + imagu*g_val(I,Z,wire,'dia'), &
               bc1_a1, bc1_a2, bc1_b2, bc1_da1)

          IF (BCI == 1) THEN
             call g_der(I,DG,wire,'ra',   real(bc1_a1))
             call g_der(I,DG,wire,'ia', -aimag(bc1_a1))
             call g_der(I,DG,wire,'dra',  real(bc1_da1))
             call g_der(I,DG,wire,'dia',-aimag(bc1_da1))

             call g_der(I,DG,wire2,'ra',   real(bc1_a2))
             call g_der(I,DG,wire2,'ia', -aimag(bc1_a2))
             call g_der(I,DG,wire2,'rb',   real(bc1_b2))
             call g_der(I,DG,wire2,'ib', -aimag(bc1_b2))
          ELSE IF (BCI == 2) THEN
             call g_der(I,DG,wire,'ra',  aimag(bc1_a1))
             call g_der(I,DG,wire,'ia',   real(bc1_a1))
             call g_der(I,DG,wire,'dra', aimag(bc1_da1))
             call g_der(I,DG,wire,'dia',  real(bc1_da1))

             call g_der(I,DG,wire2,'ra',  aimag(bc1_a2))
             call g_der(I,DG,wire2,'ia',   real(bc1_a2))
             call g_der(I,DG,wire2,'rb',  aimag(bc1_b2))
             call g_der(I,DG,wire2,'ib',   real(bc1_b2))
          END IF
       ELSE IF (BCI == 3 .OR. BCI == 4) THEN
          call eval_kl_d(rr, &
               g_val(I,Z,wire,'rb') + imagu*g_val(I,Z,wire,'ib'), &
               g_val(I,Z,wire2,'rb') + imagu*g_val(I,Z,wire2,'ib'), &
               g_val(I,Z,wire2,'ra') + imagu*g_val(I,Z,wire2,'ia'), &
               g_val(I,Z,wire,'drb') + imagu*g_val(I,Z,wire,'dib'), &
               bc1_b1, bc1_b2, bc1_a2, bc1_db1)

          IF (BCI == 3) THEN
             call g_der(I,DG,wire,'rb',   real(bc1_b1))
             call g_der(I,DG,wire,'ib', -aimag(bc1_b1))
             call g_der(I,DG,wire,'drb',  real(bc1_db1))
             call g_der(I,DG,wire,'dib',-aimag(bc1_db1))

             call g_der(I,DG,wire2,'ra',   real(bc1_a2))
             call g_der(I,DG,wire2,'ia', -aimag(bc1_a2))
             call g_der(I,DG,wire2,'rb',   real(bc1_b2))
             call g_der(I,DG,wire2,'ib', -aimag(bc1_b2))
          ELSE IF (BCI == 4) THEN
             call g_der(I,DG,wire,'rb',  aimag(bc1_b1))
             call g_der(I,DG,wire,'ib',   real(bc1_b1))
             call g_der(I,DG,wire,'drb', aimag(bc1_db1))
             call g_der(I,DG,wire,'dib',  real(bc1_db1))

             call g_der(I,DG,wire2,'ra',  aimag(bc1_a2))
             call g_der(I,DG,wire2,'ia',   real(bc1_a2))
             call g_der(I,DG,wire2,'rb',  aimag(bc1_b2))
             call g_der(I,DG,wire2,'ib',   real(bc1_b2))
          END IF
       ELSE
          WRITE(*,*) '%% Invalid BC 1: ', BCI
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_CLEAN_NODE) THEN
       wire = BCCONNECT(I, 1)
       wire2 = BCCONNECT(I, 2)

       IF (BCI == 1) THEN
          call g_der(I,DG,wire,'ra',   1d0)
          call g_der(I,DG,wire2,'ra', -1d0)
       ELSE IF (BCI == 2) THEN
          call g_der(I,DG,wire,'ia',   1d0)
          call g_der(I,DG,wire2,'ia', -1d0)
       ELSE IF (BCI == 3) THEN
          call g_der(I,DG,wire,'rb',   1d0)
          call g_der(I,DG,wire2,'rb', -1d0)
       ELSE IF (BCI == 4) THEN
          call g_der(I,DG,wire,'ib',  1d0)
          call g_der(I,DG,wire2,'ib',-1d0)
       ELSE IF (BCI == 5) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             call g_der(I,DG,wire,'dra',WIRE_CONDUCT(wire)/lratio(wire))
          END DO
       ELSE IF (BCI == 6) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             call g_der(I,DG,wire,'dia',WIRE_CONDUCT(wire)/lratio(wire))
          END DO
       ELSE IF (BCI == 7) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             call g_der(I,DG,wire,'drb',WIRE_CONDUCT(wire)/lratio(wire))
          END DO
       ELSE IF (BCI == 8) THEN
          DO J=1, NWIRE+1
             wire = BCCONNECT(I,J)
             IF (wire .LE. 0) EXIT
             call g_der(I,DG,wire,'dib',WIRE_CONDUCT(wire)/lratio(wire))
          END DO
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_FREE_INTERFACE) THEN
       wire = BCCONNECT(I, 1)

       IF (BCI == 1) THEN
          call g_der(I,DG,wire,'dra',1d0)
       ELSE IF (BCI == 2) THEN
          call g_der(I,DG,wire,'dia',1d0)
       ELSE IF (BCI == 3) THEN
          call g_der(I,DG,wire,'drb',1d0)
       ELSE IF (BCI == 4) THEN
          call g_der(I,DG,wire,'dib',1d0)
       END IF
       !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF (BCT .EQ. BCTYPE_FREE_INTERFACE_FIX_PHASE) THEN
       wire = BCCONNECT(I, 1)

       IF (BCI == 1) THEN
          call g_der(I,DG,wire,'ra', 1d0)
          call g_der(I,DG,wire,'rb',-1d0)
       ELSE IF (BCI == 2) THEN
          call g_der(I,DG,wire,'ia', 1d0)
          call g_der(I,DG,wire,'ib',-1d0)
       ELSE IF (BCI == 3) THEN
          call g_der(I,DG,wire,'drb',1d0)
       ELSE IF (BCI == 4) THEN
          call g_der(I,DG,wire,'dib',1d0)
       END IF
    ELSE
       WRITE(*,*) '%% Invalid spectral DGSUB BC', BCT
    END IF
  END SUBROUTINE DGSUB


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initial guess
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE GUESS(X, Z, DM)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(MSTAR), INTENT(OUT) :: Z
    double precision, DIMENSION(NCOMP), INTENT(OUT) :: DM

    INTEGER :: wire
    double precision, DIMENSION(NWIRE) :: delta, phase

    double complex :: c
    integer :: jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib

    DM = 0
    Z = 0

    CALL GET_DELTA(X, delta, phase)

    DO wire = 1, NWIRE
       call indices(wire, jddra, jddia, jddrb, jddib, &
            jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)
       Z(jra) =  0
       Z(jia) =  0
       Z(jrb) =  0
       Z(jib) =  0
    END DO
  END SUBROUTINE GUESS


  SUBROUTINE get_vars_from_Z(X, Z, wire, a, da, b, db)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(MSTAR), INTENT(IN) :: Z
    INTEGER, INTENT(IN) :: wire
    double complex, INTENT(OUT) :: a, da, b, db
    integer :: jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib

    call indices(wire, jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)

    a  = Z(jra)  + imagu*Z(jia)
    da = Z(jdra) + imagu*Z(jdia)

    b  = Z(jrb)  + imagu*Z(jib)
    db = Z(jdrb) + imagu*Z(jdib)

    da = da / lratio(wire)
    db = db / lratio(wire)
  END SUBROUTINE get_vars_from_Z

  SUBROUTINE set_vars_to_Z(X, Z, wire, a, da, b, db)
    IMPLICIT NONE
    double precision, INTENT(IN) :: X
    double precision, DIMENSION(MSTAR), INTENT(OUT) :: Z
    INTEGER, INTENT(IN) :: wire
    double complex, INTENT(IN) :: a, da, b, db

    integer :: jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib

    call indices(wire, jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)

    Z(jra) = real(a)
    Z(jia) = aimag(a)
    Z(jdra) = real(da) * lratio(wire)
    Z(jdia) = aimag(da) * lratio(wire)

    Z(jrb) = real(b)
    Z(jib) = aimag(b)
    Z(jdrb) = real(db) * lratio(wire)
    Z(jdib) = aimag(db) * lratio(wire)
  END SUBROUTINE set_vars_to_Z

  SUBROUTINE get_var_indices(wire, jra, jia, jdra, jdia, &
                                   jrb, jib, jdrb, jdib)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: wire
    INTEGER, INTENT(OUT) :: jra, jia, jrb, jib, jdra, jdia, jdrb, jdib
    INTEGER :: jddra, jddia, jddrb, jddib

    call indices(wire, jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)
  END SUBROUTINE get_var_indices

  SUBROUTINE SET_TOL(TOLS, TOL, DTOL, TOLB, DTOLB)
    IMPLICIT NONE
    double precision, DIMENSION(MSTAR), INTENT(OUT) :: TOLS
    double precision, INTENT(IN) :: TOL, DTOL, TOLB, DTOLB

    INTEGER :: wire
    integer :: jddra, jddia, jddrb, jddib, &
         jra, jia, jrb, jib, jdra, jdia, jdrb, jdib

    DO wire=1,NWIRE
       call indices(wire, jddra, jddia, jddrb, jddib, &
            jra, jia, jrb, jib, jdra, jdia, jdrb, jdib)

       TOLS(jra) = TOL
       TOLS(jia) = TOL
       TOLS(jdra) = DTOL
       TOLS(jdia) = DTOL

       TOLS(jrb) = TOLB
       TOLS(jib) = TOLB
       TOLS(jdrb) = DTOLB
       TOLS(jdib) = DTOLB
    END DO
  END SUBROUTINE SET_TOL
END MODULE SP_EQUATIONS
