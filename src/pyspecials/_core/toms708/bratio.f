C     ALGORITHM 708, COLLECTED ALGORITHMS FROM ACM.
C     THIS WORK WAS PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C     VOL. 18, NO. 3, SEPTEMBER, 1992, PP. 360-373z.
C-----------------------------------------------------------------------
C
C     THIS VERSION WAS MODIFIED BY LEANDRO AUGUSTO LACERDA CAMPOS ON
C     OCTOBER 19, 2023, TO MAKE ALGORITHM 708 WORK WITH DOUBLE-PRECISION
c     VARIABLES.
C
C-----------------------------------------------------------------------
      SUBROUTINE BRATIO(A, B, X, Y, W, W1, IERR)
C-----------------------------------------------------------------------
C
C            EVALUATION OF THE INCOMPLETE BETA FUNCTION IX(A,B)
C
C                     --------------------
C
C     IT IS ASSUMED THAT A AND B ARE NONNEGATIVE, AND THAT X .LE. 1
C     AND Y = 1 - X.  BRATIO ASSIGNS W AND W1 THE VALUES
C
C                      W  = IX(A,B)
C                      W1 = 1 - IX(A,B)
C
C     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
C     IF NO INPUT ERRORS ARE DETECTED THEN IERR IS SET TO 0 AND
C     W AND W1 ARE COMPUTED. OTHERWISE, IF AN ERROR IS DETECTED,
C     THEN W AND W1 ARE ASSIGNED THE VALUE 0 AND IERR IS SET TO
C     ONE OF THE FOLLOWING VALUES ...
C
C        IERR = 1  IF A OR B IS NEGATIVE
C        IERR = 2  IF A = B = 0
C        IERR = 3  IF X .LT. 0 OR X .GT. 1
C        IERR = 4  IF Y .LT. 0 OR Y .GT. 1
C        IERR = 5  IF X + Y .NE. 1
C        IERR = 6  IF X = A = 0
C        IERR = 7  IF Y = B = 0
C
C--------------------
C     WRITTEN BY ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN, VIRGINIA
C     REVISED ... NOV 1991
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B, W, W1, X, Y
         INTEGER IERR
C
         DOUBLE PRECISION A0, B0, EPS, LAMBDA, T, X0, Y0, Z
         INTEGER IERR1, IND, N
C
         DOUBLE PRECISION APSER, BASYM, BFRAC, BPSER, BUP, DPMPAR, FPSER
C-----------------------------------------------------------------------
C
C     ****** EPS IS A MACHINE DEPENDENT CONSTANT. EPS IS THE SMALLEST
C            FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
C
         EPS = DPMPAR(1)
C
C-----------------------------------------------------------------------
         W = 0.0D0
         W1 = 0.0D0
         IF (A .LT. 0.0D0 .OR. B .LT. 0.0D0) GO TO 300
         IF (A .EQ. 0.0D0 .AND. B .EQ. 0.0D0) GO TO 310
         IF (X .LT. 0.0D0 .OR. X .GT. 1.0D0) GO TO 320
         IF (Y .LT. 0.0D0 .OR. Y .GT. 1.0D0) GO TO 330
         Z = ((X + Y) - 0.5D0) - 0.5D0
         IF (ABS(Z) .GT. 3.0D0*EPS) GO TO 340
C
         IERR = 0
         IF (X .EQ. 0.0D0) GO TO 200
         IF (Y .EQ. 0.0D0) GO TO 210
         IF (A .EQ. 0.0D0) GO TO 211
         IF (B .EQ. 0.0D0) GO TO 201
C
         EPS = DMAX1(EPS, 1.0D-15)
         IF (DMAX1(A,B) .LT. 1.0D-3*EPS) GO TO 230
C
         IND = 0
         A0 = A
         B0 = B
         X0 = X
         Y0 = Y
         IF (DMIN1(A0, B0) .GT. 1.0D0) GO TO 30
C
C             PROCEDURE FOR A0 .LE. 1 OR B0 .LE. 1
C
         IF (X .LE. 0.5D0) GO TO 10
         IND = 1
         A0 = B
         B0 = A
         X0 = Y
         Y0 = X
C
   10    IF (B0 .LT. DMIN1(EPS,EPS*A0)) GO TO 80
         IF (A0 .LT. DMIN1(EPS,EPS*B0) .AND. B0*X0 .LE. 1.0D0) GO TO 90
         IF (DMAX1(A0, B0) .GT. 1.0D0) GO TO 20
         IF (A0 .GE. DMIN1(0.2D0, B0)) GO TO 100
         IF (X0**A0 .LE. 0.9D0) GO TO 100
         IF (X0 .GE. 0.3D0) GO TO 110
         N = 20
         GO TO 130
C
   20    IF (B0 .LE. 1.0D0) GO TO 100
         IF (X0 .GE. 0.3D0) GO TO 110
         IF (X0 .GE. 0.1D0) GO TO 21
         IF ((X0*B0)**A0 .LE. 0.7D0) GO TO 100
   21    IF (B0 .GT. 15.0D0) GO TO 131
         N = 20
         GO TO 130
C
C             PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
C
   30    IF (A .GT. B) GO TO 31
         LAMBDA = A - (A + B)*X
         GO TO 32
   31    LAMBDA = (A + B)*Y - B
   32    IF (LAMBDA .GE. 0.0D0) GO TO 40
         IND = 1
         A0 = B
         B0 = A
         X0 = Y
         Y0 = X
         LAMBDA = ABS(LAMBDA)
C
   40    IF (B0 .LT. 40.0D0 .AND. B0*X0 .LE. 0.7D0) GO TO 100
         IF (B0 .LT. 40.0D0) GO TO 140
         IF (A0 .GT. B0) GO TO 50
         IF (A0 .LE. 100.0D0) GO TO 120
         IF (LAMBDA .GT. 0.03D0*A0) GO TO 120
         GO TO 180
   50    IF (B0 .LE. 100.0D0) GO TO 120
         IF (LAMBDA .GT. 0.03D0*B0) GO TO 120
         GO TO 180
C
C            EVALUATION OF THE APPROPRIATE ALGORITHM
C
   80    W = FPSER(A0, B0, X0, EPS)
         W1 = 0.5D0 + (0.5D0 - W)
         GO TO 220
C
   90    W1 = APSER(A0, B0, X0, EPS)
         W = 0.5D0 + (0.5D0 - W1)
         GO TO 220
C
  100    W = BPSER(A0, B0, X0, EPS)
         W1 = 0.5D0 + (0.5D0 - W)
         GO TO 220
C
  110    W1 = BPSER(B0, A0, Y0, EPS)
         W = 0.5D0 + (0.5D0 - W1)
         GO TO 220
C
  120    W = BFRAC(A0, B0, X0, Y0, LAMBDA, 15.0D0*EPS)
         W1 = 0.5D0 + (0.5D0 - W)
         GO TO 220
C
  130    W1 = BUP(B0, A0, Y0, X0, N, EPS)
         B0 = B0 + N
  131    CALL BGRAT(B0, A0, Y0, X0, W1, 15.0D0*EPS, IERR1)
         W = 0.5D0 + (0.5D0 - W1)
         GO TO 220
C
  140    N = INT(B0)
         B0 = B0 - N
         IF (B0 .NE. 0.0D0) GO TO 141
         N = N - 1
         B0 = 1.0D0
  141    W = BUP(B0, A0, Y0, X0, N, EPS)
         IF (X0 .GT. 0.7D0) GO TO 150
         W = W + BPSER(A0, B0, X0, EPS)
         W1 = 0.5D0 + (0.5D0 - W)
         GO TO 220
C
  150    IF (A0 .GT. 15.0) GO TO 151
         N = 20
         W = W + BUP(A0, B0, X0, Y0, N, EPS)
         A0 = A0 + N
  151    CALL BGRAT(A0, B0, X0, Y0, W, 15.0D0*EPS, IERR1)
         W1 = 0.5D0 + (0.5D0 - W)
         GO TO 220
C
  180    W = BASYM(A0, B0, LAMBDA, 100.0D0*EPS)
         W1 = 0.5D0 + (0.5D0 - W)
         GO TO 220
C
C               TERMINATION OF THE PROCEDURE
C
  200    IF (A .EQ. 0.0D0) GO TO 350
  201    W = 0.0D0
         W1 = 1.0D0
         RETURN
C
  210    IF (B .EQ. 0.0D0) GO TO 360
  211    W = 1.0D0
         W1 = 0.0D0
         RETURN
C
  220    IF (IND .EQ. 0) RETURN
         T = W
         W = W1
         W1 = T
         RETURN
C
C           PROCEDURE FOR A AND B .LT. 1.D-3*EPS
C
  230    W = B/(A + B)
         W1 = A/(A + B)
         RETURN
C
C                       ERROR RETURN
C
  300    IERR = 1
         RETURN
  310    IERR = 2
         RETURN
  320    IERR = 3
         RETURN
  330    IERR = 4
         RETURN
  340    IERR = 5
         RETURN
  350    IERR = 6
         RETURN
  360    IERR = 7
         RETURN
      END
      DOUBLE PRECISION FUNCTION FPSER(A, B, X, EPS)
C-----------------------------------------------------------------------
C
C                 EVALUATION OF I (A,B)
C                                X
C
C          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.
C
C-----------------------------------------------------------------------
C
C                  SET  FPSER = X**A
C
         IMPLICIT NONE
         DOUBLE PRECISION A, B, EPS, X
C
         DOUBLE PRECISION AN, C, S, T, TOL
C
         DOUBLE PRECISION EXPARG
C
         FPSER = 1.0D0
         IF (A .LE. 10D-3*EPS) GO TO 10
         FPSER = 0.0D0
         T = A*DLOG(X)
         IF (T .LT. EXPARG(1)) RETURN
         FPSER = EXP(T)
C
C                NOTE THAT 1/B(A,B) = B
C
   10    FPSER = (B/A)*FPSER
         TOL = EPS/A
         AN = A + 1.0D0
         T = X
         S = T/AN
   20    AN = AN + 1.0D0
         T = X*T
         C = T/AN
         S = S + C
         IF (ABS(C) .GT. TOL) GO TO 20
C
         FPSER = FPSER*(1.0D0 + A*S)
         RETURN
      END
      DOUBLE PRECISION FUNCTION APSER(A, B, X, EPS)
C-----------------------------------------------------------------------
C     APSER YIELDS THE INCOMPLETE BETA RATIO I(SUB(1-X))(B,A) FOR
C     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. USED WHEN
C     A IS VERY SMALL. USE ONLY IF ABOVE INEQUALITIES ARE SATISFIED.
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B, EPS, X
C
         DOUBLE PRECISION AJ, BX, C, G, J, S, T, TOL
C
         DOUBLE PRECISION PSI
C--------------------
         DATA G/.577215664901533D0/
C--------------------
         BX = B*X
         T = X - BX
         IF (B*EPS .GT. 2.0D-2) GO TO 10
         C = DLOG(X) + PSI(B) + G + T
         GO TO 20
   10    C = DLOG(BX) + G + T
C
   20    TOL = 5.0D0*EPS*ABS(C)
         J = 1.0D0
         S = 0.0D0
   30    J = J + 1.0D0
         T = T*(X - BX/J)
         AJ = T/J
         S = S + AJ
         IF (ABS(AJ) .GT. TOL) GO TO 30
C
         APSER = -A*(C + S)
         RETURN
      END
      DOUBLE PRECISION FUNCTION BPSER(A, B, X, EPS)
C-----------------------------------------------------------------------
C     POWER SERIES EXPANSION FOR EVALUATING IX(A,B) WHEN B .LE. 1
C     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B, EPS, X
         DOUBLE PRECISION A0, APB, B0, C, N, SUM, T, TOL, U, W, Z
         INTEGER I, M
C
         DOUBLE PRECISION ALGDIV, BETALN, GAM1, GAMLN1
C
         BPSER = 0.0D0
         IF (X .EQ. 0.0D0) RETURN
C-----------------------------------------------------------------------
C            COMPUTE THE FACTOR X**A/(A*BETA(A,B))
C-----------------------------------------------------------------------
         A0 = DMIN1(A,B)
         IF (A0 .LT. 1.0D0) GO TO 10
         Z = A*DLOG(X) - BETALN(A,B)
         BPSER = EXP(Z)/A
         GO TO 70
   10    B0 = DMAX1(A,B)
         IF (B0 .GE. 8.0D0) GO TO 60
         IF (B0 .GT. 1.0D0) GO TO 40
C
C            PROCEDURE FOR A0 .LT. 1 AND B0 .LE. 1
C
         BPSER = X**A
         IF (BPSER .EQ. 0.0D0) RETURN
C
         APB = A + B
         IF (APB .GT. 1.0D0) GO TO 20
         Z = 1.0D0 + GAM1(APB)
         GO TO 30
   20    U = A + B - 1.D0
         Z = (1.0D0 + GAM1(U))/APB
C
   30    C = (1.0D0 + GAM1(A))*(1.0D0 + GAM1(B))/Z
         BPSER = BPSER*C*(B/APB)
         GO TO 70
C
C         PROCEDURE FOR A0 .LT. 1 AND 1 .LT. B0 .LT. 8
C
   40    U = GAMLN1(A0)
         M = INT(B0) - 1
         IF (M .LT. 1) GO TO 50
         C = 1.0D0
         DO I = 1,M
            B0 = B0 - 1.0D0
            C = C*(B0/(A0 + B0))
         END DO
         U = DLOG(C) + U
C
   50    Z = A*DLOG(X) - U
         B0 = B0 - 1.0D0
         APB = A0 + B0
         IF (APB .GT. 1.0D0) GO TO 51
         T = 1.0D0 + GAM1(APB)
         GO TO 52
   51    U = A0 + B0 - 1.0D0
         T = (1.0D0 + GAM1(U))/APB
   52    BPSER = EXP(Z)*(A0/A)*(1.0D0 + GAM1(B0))/T
         GO TO 70
C
C            PROCEDURE FOR A0 .LT. 1 AND B0 .GE. 8
C
   60    U = GAMLN1(A0) + ALGDIV(A0,B0)
         Z = A*DLOG(X) - U
         BPSER = (A0/A)*EXP(Z)
   70    IF (BPSER .EQ. 0.0D0 .OR. A .LE. 0.1D0*EPS) RETURN
C-----------------------------------------------------------------------
C                     COMPUTE THE SERIES
C-----------------------------------------------------------------------
         SUM = 0.0D0
         N = 0.0
         C = 1.0D0
         TOL = EPS/A
  100    N = N + 1.0D0
         C = C*(0.5D0 + (0.5D0 - B/N))*X
         W = C/(A + N)
         SUM = SUM + W
         IF (ABS(W) .GT. TOL) GO TO 100
         BPSER = BPSER*(1.0D0 + A*SUM)
         RETURN
      END
      DOUBLE PRECISION FUNCTION BUP(A, B, X, Y, N, EPS)
C-----------------------------------------------------------------------
C     EVALUATION OF IX(A,B) - IX(A+N,B) WHERE N IS A POSITIVE INTEGER.
C     EPS IS THE TOLERANCE USED.
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B, EPS, X, Y
         INTEGER N
C
         DOUBLE PRECISION AP1, APB, D, L, R, T, W
         INTEGER I, K, KP1, MU, NM1
C
         DOUBLE PRECISION BRCMP1, EXPARG
C
C          OBTAIN THE SCALING FACTOR EXP(-MU) AND
C             EXP(MU)*(X**A*Y**B/BETA(A,B))/A
C
         APB = A + B
         AP1 = A + 1.0D0
         MU = 0
         D = 1.0D0
         IF (N .EQ. 1 .OR. A .LT. 1.0D0) GO TO 10
         IF (APB .LT. 1.1D0*AP1) GO TO 10
         MU = ABS(INT(EXPARG(1)))
         K = INT(EXPARG(0))
         IF (K .LT. MU) MU = K
         T = MU
         D = EXP(-T)
C
   10    BUP = BRCMP1(MU,A,B,X,Y)/A
         IF (N .EQ. 1 .OR. BUP .EQ. 0.0D0) RETURN
         NM1 = N - 1
         W = D
C
C          LET K BE THE INDEX OF THE MAXIMUM TERM
C
         K = 0
         IF (B .LE. 1.0D0) GO TO 40
         IF (Y .GT. 1.D-4) GO TO 20
         K = NM1
         GO TO 30
   20    R = (B - 1.0D0)*X/Y - A
         IF (R .LT. 1.0D0) GO TO 40
         K = NM1
         T = NM1
         IF (R .LT. T) K = INT(R)
C
C          ADD THE INCREASING TERMS OF THE SERIES
C
   30    DO 31 I = 1,K
            L = I - 1
            D = ((APB + L)/(AP1 + L))*X*D
            W = W + D
   31    CONTINUE
         IF (K .EQ. NM1) GO TO 50
C
C          ADD THE REMAINING TERMS OF THE SERIES
C
   40    KP1 = K + 1
         DO 41 I = KP1,NM1
            L = I - 1
            D = ((APB + L)/(AP1 + L))*X*D
            W = W + D
            IF (D .LE. EPS*W) GO TO 50
   41    CONTINUE
C
C               TERMINATE THE PROCEDURE
C
   50    BUP = BUP*W
         RETURN
      END
      DOUBLE PRECISION FUNCTION BFRAC(A, B, X, Y, LAMBDA, EPS)
C-----------------------------------------------------------------------
C     CONTINUED FRACTION EXPANSION FOR IX(A,B) WHEN A,B .GT. 1.
C     IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B, EPS, LAMBDA, X, Y
C
         DOUBLE PRECISION ALPHA, AN, ANP1, BETA, BN, BNP1, C, C0, C1, E,
     +                    N, P, R, R0, S, T, W, YP1
C
         DOUBLE PRECISION BRCOMP
C
         BFRAC = BRCOMP(A,B,X,Y)
         IF (BFRAC .EQ. 0.0D0) RETURN
C
         C = 1.0D0 + LAMBDA
         C0 = B/A
         C1 = 1.0D0 + 1.0D0/A
         YP1 = Y + 1.0D0
C
         N = 0.0D0
         P = 1.0D0
         S = A + 1.0D0
         AN = 0.0D0
         BN = 1.0D0
         ANP1 = 1.0D0
         BNP1 = C/C1
         R = C1/C
C
C        CONTINUED FRACTION CALCULATION
C
   10    N = N + 1.0D0
         T = N/A
         W = N*(B - N)*X
         E = A/S
         ALPHA = (P*(P + C0)*E*E)*(W*X)
         E = (1.0D0 + T)/(C1 + T + T)
         BETA = N + W/S + E*(C + N*YP1)
         P = 1.0D0 + T
         S = S + 2.0D0
C
C        UPDATE AN, BN, ANP1, AND BNP1
C
         T = ALPHA*AN + BETA*ANP1
         AN = ANP1
         ANP1 = T
         T = ALPHA*BN + BETA*BNP1
         BN = BNP1
         BNP1 = T
C
         R0 = R
         R = ANP1/BNP1
         IF (ABS(R - R0) .LE. EPS*R) GO TO 20
C
C        RESCALE AN, BN, ANP1, AND BNP1
C
         AN = AN/BNP1
         BN = BN/BNP1
         ANP1 = R
         BNP1 = 1.0D0
         GO TO 10
C
C                 TERMINATION
C
   20    BFRAC = BFRAC*R
         RETURN
      END
      DOUBLE PRECISION FUNCTION BRCOMP(A, B, X, Y)
C-----------------------------------------------------------------------
C               EVALUATION OF X**A*Y**B/BETA(A,B)
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B, X, Y
C
         DOUBLE PRECISION A0, APB, B0, C, CONST, E, H, LAMBDA, LNX, LNY,
     +                    T, U, V, X0, Y0, Z
         INTEGER I, N
C
         DOUBLE PRECISION ALGDIV, ALNREL, BCORR, BETALN, GAM1, GAMLN1,
     +                    RLOG1
C-----------------
C     CONST = 1/SQRT(2*PI)
C-----------------
         DATA CONST/.398942280401433D0/
C
         BRCOMP = 0.0D0
         IF (X .EQ. 0.0D0 .OR. Y .EQ. 0.0D0) RETURN
         A0 = DMIN1(A,B)
         IF (A0 .GE. 8.0D0) GO TO 100
C
         IF (X .GT. 0.375D0) GO TO 10
         LNX = DLOG(X)
         LNY = ALNREL(-X)
         GO TO 20
   10    IF (Y .GT. 0.375D0) GO TO 11
         LNX = ALNREL(-Y)
         LNY = DLOG(Y)
         GO TO 20
   11    LNX = DLOG(X)
         LNY = DLOG(Y)
C
   20    Z = A*LNX + B*LNY
         IF (A0 .LT. 1.0D0) GO TO 30
         Z = Z - BETALN(A,B)
         BRCOMP = EXP(Z)
         RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .LT. 1 OR B .LT. 1
C-----------------------------------------------------------------------
   30    B0 = DMAX1(A,B)
         IF (B0 .GE. 8.0D0) GO TO 80
         IF (B0 .GT. 1.0D0) GO TO 60
C
C                   ALGORITHM FOR B0 .LE. 1
C
         BRCOMP = EXP(Z)
         IF (BRCOMP .EQ. 0.0D0) RETURN
C
         APB = A + B
         IF (APB .GT. 1.0D0) GO TO 40
         Z = 1.0D0 + GAM1(APB)
         GO TO 50
   40    U = A + B - 1.0D0
         Z = (1.0D0 + GAM1(U))/APB
C
   50    C = (1.0D0 + GAM1(A))*(1.0D0 + GAM1(B))/Z
         BRCOMP = BRCOMP*(A0*C)/(1.0D0 + A0/B0)
         RETURN
C
C                ALGORITHM FOR 1 .LT. B0 .LT. 8
C
   60    U = GAMLN1(A0)
         N = INT(B0) - 1
         IF (N .LT. 1) GO TO 70
         C = 1.0D0
         DO 61 I = 1,N
            B0 = B0 - 1.0D0
            C = C*(B0/(A0 + B0))
   61    CONTINUE
         U = DLOG(C) + U
C
   70    Z = Z - U
         B0 = B0 - 1.0D0
         APB = A0 + B0
         IF (APB .GT. 1.0D0) GO TO 71
         T = 1.0D0 + GAM1(APB)
         GO TO 72
   71    U = A0 + B0 - 1.0D0
         T = (1.0D0 + GAM1(U))/APB
   72    BRCOMP = A0*EXP(Z)*(1.0D0 + GAM1(B0))/T
         RETURN
C
C                   ALGORITHM FOR B0 .GE. 8
C
   80    U = GAMLN1(A0) + ALGDIV(A0,B0)
         BRCOMP = A0*EXP(Z - U)
         RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .GE. 8 AND B .GE. 8
C-----------------------------------------------------------------------
  100    IF (A .GT. B) GO TO 101
         H = A/B
         X0 = H/(1.0D0 + H)
         Y0 = 1.0D0/(1.0D0 + H)
         LAMBDA = A - (A + B)*X
         GO TO 110
  101    H = B/A
         X0 = 1.0D0/(1.0D0 + H)
         Y0 = H/(1.0D0 + H)
         LAMBDA = (A + B)*Y - B
C
  110    E = -LAMBDA/A
         IF (ABS(E) .GT. 0.6D0) GO TO 111
         U = RLOG1(E)
         GO TO 120
  111    U = E - DLOG(X/X0)
C
  120    E = LAMBDA/B
         IF (ABS(E) .GT. 0.6D0) GO TO 121
         V = RLOG1(E)
         GO TO 130
  121    V = E - DLOG(Y/Y0)
C
  130    Z = EXP(-(A*U + B*V))
         BRCOMP = CONST*SQRT(B*X0)*Z*EXP(-BCORR(A,B))
         RETURN
      END
      DOUBLE PRECISION FUNCTION BRCMP1(MU, A, B, X, Y)
C-----------------------------------------------------------------------
C          EVALUATION OF  EXP(MU) * (X**A*Y**B/BETA(A,B))
C-----------------------------------------------------------------------
         IMPLICIT NONE
         INTEGER MU
         DOUBLE PRECISION A, B, X, Y
C
         DOUBLE PRECISION A0, APB, B0, C, CONST, E, H, LAMBDA, LNX, LNY,
     +                    T, U, V, X0, Y0, Z
C
         INTEGER I, N
C
         DOUBLE PRECISION ALGDIV, ALNREL, BCORR, BETALN, ESUM, GAM1,
     +                    GAMLN1, RLOG1
C-----------------
C     CONST = 1/SQRT(2*PI)
C-----------------
         DATA CONST/.398942280401433D0/
C
         A0 = DMIN1(A,B)
         IF (A0 .GE. 8.0D0) GO TO 100
C
         IF (X .GT. 0.375D0) GO TO 10
         LNX = DLOG(X)
         LNY = ALNREL(-X)
         GO TO 20
   10    IF (Y .GT. 0.375D0) GO TO 11
         LNX = ALNREL(-Y)
         LNY = DLOG(Y)
         GO TO 20
   11    LNX = DLOG(X)
         LNY = DLOG(Y)
C
   20    Z = A*LNX + B*LNY
         IF (A0 .LT. 1.0D0) GO TO 30
         Z = Z - BETALN(A,B)
         BRCMP1 = ESUM(MU,Z)
         RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .LT. 1 OR B .LT. 1
C-----------------------------------------------------------------------
   30    B0 = DMAX1(A,B)
         IF (B0 .GE. 8.0D0) GO TO 80
         IF (B0 .GT. 1.0D0) GO TO 60
C
C                   ALGORITHM FOR B0 .LE. 1
C
         BRCMP1 = ESUM(MU,Z)
         IF (BRCMP1 .EQ. 0.0D0) RETURN
C
         APB = A + B
         IF (APB .GT. 1.0D0) GO TO 40
         Z = 1.0D0 + GAM1(APB)
         GO TO 50
   40    U = A + B - 1.0D0
         Z = (1.0D0 + GAM1(U))/APB
C
   50    C = (1.0D0 + GAM1(A))*(1.0D0 + GAM1(B))/Z
         BRCMP1 = BRCMP1*(A0*C)/(1.0D0 + A0/B0)
         RETURN
C
C                ALGORITHM FOR 1 .LT. B0 .LT. 8
C
   60    U = GAMLN1(A0)
         N = INT(B0) - 1
         IF (N .LT. 1) GO TO 70
         C = 1.0D0
         DO 61 I = 1,N
            B0 = B0 - 1.0D0
            C = C*(B0/(A0 + B0))
   61    CONTINUE
         U = DLOG(C) + U
C
   70    Z = Z - U
         B0 = B0 - 1.0D0
         APB = A0 + B0
         IF (APB .GT. 1.0D0) GO TO 71
         T = 1.0D0 + GAM1(APB)
         GO TO 72
   71    U = A0 + B0 - 1.0D0
         T = (1.0D0 + GAM1(U))/APB
   72    BRCMP1 = A0*ESUM(MU,Z)*(1.0D0 + GAM1(B0))/T
         RETURN
C
C                   ALGORITHM FOR B0 .GE. 8
C
   80    U = GAMLN1(A0) + ALGDIV(A0,B0)
         BRCMP1 = A0*ESUM(MU,Z - U)
         RETURN
C-----------------------------------------------------------------------
C              PROCEDURE FOR A .GE. 8 AND B .GE. 8
C-----------------------------------------------------------------------
  100    IF (A .GT. B) GO TO 101
         H = A/B
         X0 = H/(1.0D0 + H)
         Y0 = 1.0D0/(1.0D0 + H)
         LAMBDA = A - (A + B)*X
         GO TO 110
  101    H = B/A
         X0 = 1.0D0/(1.0D0 + H)
         Y0 = H/(1.0D0 + H)
         LAMBDA = (A + B)*Y - B
C
  110    E = -LAMBDA/A
         IF (ABS(E) .GT. 0.6D0) GO TO 111
         U = RLOG1(E)
         GO TO 120
  111    U = E - DLOG(X/X0)
C
  120    E = LAMBDA/B
         IF (ABS(E) .GT. 0.6D0) GO TO 121
         V = RLOG1(E)
         GO TO 130
  121    V = E - DLOG(Y/Y0)
C
  130    Z = ESUM(MU,-(A*U + B*V))
         BRCMP1 = CONST*SQRT(B*X0)*Z*EXP(-BCORR(A,B))
         RETURN
      END
      SUBROUTINE BGRAT(A, B, X, Y, W, EPS, IERR)
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR IX(A,B) WHEN A IS LARGER THAN B.
C     THE RESULT OF THE EXPANSION IS ADDED TO W. IT IS ASSUMED
C     THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED.
C     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B, EPS, W, X, Y
         INTEGER IERR
C
         DOUBLE PRECISION BM1, BP2N, CN, COEF, DJ, J, L, LNX, N2, NU, P,
     +                    Q, R, S, SUM, T, T2, U, V, Z
         INTEGER I, N, NM1
C
         DOUBLE PRECISION C(30), D(30)
C
         DOUBLE PRECISION ALGDIV, ALNREL, GAM1
C
         BM1 = (B - 0.5D0) - 0.5D0
         NU = A + 0.5D0*BM1
         IF (Y .GT. 0.375D0) GO TO 10
         LNX = ALNREL(-Y)
         GO TO 11
   10    LNX = DLOG(X)
   11    Z = -NU*LNX
         IF (B*Z .EQ. 0.0D0) GO TO 100
C
C                 COMPUTATION OF THE EXPANSION
C                 SET R = EXP(-Z)*Z**B/GAMMA(B)
C
         R = B*(1.0D0 + GAM1(B))*EXP(B*DLOG(Z))
         R = R*EXP(A*LNX)*EXP(0.5D0*BM1*LNX)
         U = ALGDIV(B,A) + B*DLOG(NU)
         U = R*EXP(-U)
         IF (U .EQ. 0.0D0) GO TO 100
         CALL GRAT1(B,Z,R,P,Q,EPS)
C
         V = 0.25D0*(1.0D0/NU)**2
         T2 = 0.25D0*LNX*LNX
         L = W/U
         J = Q/R
         SUM = J
         T = 1.0D0
         CN = 1.0D0
         N2 = 0.0D0
         DO 22 N = 1,30
            BP2N = B + N2
            J = (BP2N*(BP2N + 1.0D0)*J + (Z + BP2N + 1.0D0)*T)*V
            N2 = N2 + 2.0D0
            T = T*T2
            CN = CN/(N2*(N2 + 1.0D0))
            C(N) = CN
            S = 0.0D0
            IF (N .EQ. 1) GO TO 21
            NM1 = N - 1
            COEF = B - N
            DO I = 1,NM1
               S = S + COEF*C(I)*D(N-I)
               COEF = COEF + B
            END DO
   21       D(N) = BM1*CN + S/N
            DJ = D(N)*J
            SUM = SUM + DJ
            IF (SUM .LE. 0.0D0) GO TO 100
            IF (ABS(DJ) .LE. EPS*(SUM + L)) GO TO 30
   22    CONTINUE
C
C                    ADD THE RESULTS TO W
C
   30    IERR = 0
         W = W + U*SUM
         RETURN
C
C               THE EXPANSION CANNOT BE COMPUTED
C
  100    IERR = 1
         RETURN
      END
      SUBROUTINE GRAT1(A,X,R,P,Q,EPS)
C-----------------------------------------------------------------------
C        EVALUATION OF THE INCOMPLETE GAMMA RATIO FUNCTIONS
C                      P(A,X) AND Q(A,X)
C
C     IT IS ASSUMED THAT A .LE. 1.  EPS IS THE TOLERANCE TO BE USED.
C     THE INPUT ARGUMENT R HAS THE VALUE E**(-X)*X**A/GAMMA(A).
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, EPS, P, Q, R, X
C
         DOUBLE PRECISION A2N, A2NM1, AM0, AN, AN0, B2N, B2NM1, C, CMA,
     +                    G, H, J, L, SUM, T, TOL, W, Z
C
         DOUBLE PRECISION CUSTOM_ERF, ERFC1, GAM1, REXP
C
         IF (A*X .EQ. 0.0D0) GO TO 130
         IF (A .EQ. 0.5D0) GO TO 120
         IF (X .LT. 1.1D0) GO TO 10
         GO TO 50
C
C             TAYLOR SERIES FOR P(A,X)/X**A
C
   10    AN = 3.0D0
         C = X
         SUM = X/(A + 3.0D0)
         TOL = 0.1*EPS/(A + 1.0D0)
   11    AN = AN + 1.0D0
         C = -C*(X/AN)
         T = C/(A + AN)
         SUM = SUM + T
         IF (ABS(T) .GT. TOL) GO TO 11
         J = A*X*((SUM/6.0D0 - 0.5D0/(A + 2.0D0))*X + 1.0D0/(A + 1.0D0))
C
         Z = A*DLOG(X)
         H = GAM1(A)
         G = 1.0D0 + H
         IF (X .LT. 0.25D0) GO TO 20
         IF (A .LT. X/2.59D0) GO TO 40
         GO TO 30
   20    IF (Z .GT. -.13394D0) GO TO 40
C
   30    W = EXP(Z)
         P = W*G*(0.5D0 + (0.5D0 - J))
         Q = 0.5D0 + (0.5D0 - P)
         RETURN
C
   40    L = REXP(Z)
         W = 0.5D0 + (0.5D0 + L)
         Q = (W*J - L)*G - H
         IF (Q .LT. 0.0D0) GO TO 110
         P = 0.5D0 + (0.5D0 - Q)
         RETURN
C
C              CONTINUED FRACTION EXPANSION
C
   50    A2NM1 = 1.0D0
         A2N = 1.0D0
         B2NM1 = X
         B2N = X + (1.0D0 - A)
         C = 1.0D0
   51    A2NM1 = X*A2N + C*A2NM1
         B2NM1 = X*B2N + C*B2NM1
         AM0 = A2NM1/B2NM1
         C = C + 1.0D0
         CMA = C - A
         A2N = A2NM1 + CMA*A2N
         B2N = B2NM1 + CMA*B2N
         AN0 = A2N/B2N
         IF (ABS(AN0 - AM0) .GE. EPS*AN0) GO TO 51
         Q = R*AN0
         P = 0.5D0 + (0.5D0 - Q)
         RETURN
C
C                SPECIAL CASES
C
  100    P = 0.0D0
         Q = 1.0D0
         RETURN
C
  110    P = 1.0D0
         Q = 0.0D0
         RETURN
C
  120    IF (X .GE. 0.25D0) GO TO 121
         P = CUSTOM_ERF(SQRT(X))
         Q = 0.5D0 + (0.5D0 - P)
         RETURN
  121    Q = ERFC1(0,SQRT(X))
         P = 0.5D0 + (0.5D0 - Q)
         RETURN
C
  130    IF (X .LE. A) GO TO 100
         GO TO 110
      END
      DOUBLE PRECISION FUNCTION BASYM(A, B, LAMBDA, EPS)
C-----------------------------------------------------------------------
C     ASYMPTOTIC EXPANSION FOR IX(A,B) FOR LARGE A AND B.
C     LAMBDA = (A + B)*Y - B  AND EPS IS THE TOLERANCE USED.
C     IT IS ASSUMED THAT LAMBDA IS NONNEGATIVE AND THAT
C     A AND B ARE GREATER THAN OR EQUAL TO 15.
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B, EPS, LAMBDA
         DOUBLE PRECISION BSUM,DSUM, E0, E1, F, H, H2, HN, J0, J1, R,
     +                    R0, R1, S, SUM, T, T0, T1, U, W, W0, Z, Z0,
     +                    Z2, ZN, ZNM1
         INTEGER I, IM1, IMJ, J, M, MM1, MMJ, N, NP1, NUM
C
         DOUBLE PRECISION A0(21), B0(21), C(21), D(21)
C
         DOUBLE PRECISION BCORR, ERFC1, RLOG1
C------------------------
C     ****** NUM IS THE MAXIMUM VALUE THAT N CAN TAKE IN THE DO LOOP
C            ENDING AT STATEMENT 50. IT IS REQUIRED THAT NUM BE EVEN.
C            THE ARRAYS A0, B0, C, D HAVE DIMENSION NUM + 1.
C
         DATA NUM/20/
C------------------------
C     E0 = 2/SQRT(PI)
C     E1 = 2**(-3/2)
C------------------------
         DATA E0/1.12837916709551D0/, E1/.353553390593274D0/
C------------------------
         BASYM = 0.0D0
         IF (A .GE. B) GO TO 10
         H = A/B
         R0 = 1.0D0/(1.0D0 + H)
         R1 = (B - A)/B
         W0 = 1.0D0/SQRT(A*(1.0D0 + H))
         GO TO 20
   10    H = B/A
         R0 = 1.0D0/(1.0D0 + H)
         R1 = (B - A)/A
         W0 = 1.0D0/SQRT(B*(1.0D0 + H))
C
   20    F = A*RLOG1(-LAMBDA/A) + B*RLOG1(LAMBDA/B)
         T = EXP(-F)
         IF (T .EQ. 0.0D0) RETURN
         Z0 = SQRT(F)
         Z = 0.5D0*(Z0/E1)
         Z2 = F + F
C
         A0(1) = (2.0D0/3.0D0)*R1
         C(1) = - 0.5D0*A0(1)
         D(1) = - C(1)
         J0 = (0.5D0/E0)*ERFC1(1,Z0)
         J1 = E1
         SUM = J0 + D(1)*W0*J1
C
         S = 1.0D0
         H2 = H*H
         HN = 1.0D0
         W = W0
         ZNM1 = Z
         ZN = Z2
         DO 50 N = 2, NUM, 2
            HN = H2*HN
            A0(N) = 2.0D0*R0*(1.0D0 + H*HN)/(N + 2.0D0)
            NP1 = N + 1
            S = S + HN
            A0(NP1) = 2.0D0*R1*S/(N + 3.0D0)
C
            DO I = N, NP1
               R = -0.5D0*(I + 1.0D0)
               B0(1) = R*A0(1)
               DO M = 2, I
                  BSUM = 0.0D0
                  MM1 = M - 1
                  DO J = 1, MM1
                     MMJ = M - J
                     BSUM = BSUM + (J*R - MMJ)*A0(J)*B0(MMJ)
                  END DO
                  B0(M) = R*A0(M) + BSUM/M
               END DO
               C(I) = B0(I)/(I + 1.0D0)
C
               DSUM = 0.0D0
               IM1 = I - 1
               DO J = 1, IM1
                  IMJ = I - J
                  DSUM = DSUM + D(IMJ)*C(J)
               END DO
               D(I) = -(DSUM + C(I))
            END DO
C
            J0 = E1*ZNM1 + (N - 1.0D0)*J0
            J1 = E1*ZN + N*J1
            ZNM1 = Z2*ZNM1
            ZN = Z2*ZN
            W = W0*W
            T0 = D(N)*W*J0
            W = W0*W
            T1 = D(NP1)*W*J1
            SUM = SUM + (T0 + T1)
            IF ((ABS(T0) + ABS(T1)) .LE. EPS*SUM) GO TO 60
   50    CONTINUE
C
   60    U = EXP(-BCORR(A,B))
         BASYM = E0*T*U*SUM
         RETURN
      END
      DOUBLE PRECISION FUNCTION DPMPAR(I)
C-----------------------------------------------------------------------
C
C     DPMPAR PROVIDES THE DOUBLE PRECISION MACHINE CONSTANTS FOR
C     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
C     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
C     DOUBLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
C     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
C
C        DPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
C
C        DPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
C
C        DPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
C
C-----------------------------------------------------------------------
C     WRITTEN BY
C        ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN VIRGINIA
C-----------------------------------------------------------------------
         IMPLICIT NONE
         INTEGER I
         DOUBLE PRECISION B, BINV, BM1, ONE, W, Z
         INTEGER EMIN, EMAX, IBETA, M
C
         INTEGER IPMPAR
C
         IF (I .GT. 1) GO TO 10
         B = IPMPAR(4)
         M = IPMPAR(8)
         DPMPAR = B**(1 - M)
         RETURN
C
   10    IF (I .GT. 2) GO TO 20
         B = IPMPAR(4)
         EMIN = IPMPAR(9)
         ONE = DBLE(1)
         BINV = ONE/B
         W = B**(EMIN + 2)
         DPMPAR = ((W * BINV) * BINV) * BINV
         RETURN
C
   20    IBETA = IPMPAR(4)
         M = IPMPAR(8)
         EMAX = IPMPAR(10)
C
         B = IBETA
         BM1 = IBETA - 1
         ONE = DBLE(1)
         Z = B**(M - 1)
         W = ((Z - ONE)*B + BM1)/(B*Z)
C
         Z = B**(EMAX - 2)
         DPMPAR = ((W * Z) * B) * B
         RETURN
      END
      DOUBLE PRECISION FUNCTION EXPARG(L)
C--------------------------------------------------------------------
C     IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
C     EXP(W) CAN BE COMPUTED.
C
C     IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
C     WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.
C
C     NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
C--------------------------------------------------------------------
         IMPLICIT NONE
         INTEGER L
C
         DOUBLE PRECISION LNB
         INTEGER M
C
         INTEGER IPMPAR
C
         LNB = .69314718055995D0
C
         IF (L .EQ. 0) GO TO 10
         M = IPMPAR(9) - 1
         EXPARG = 0.99999D0 * (M * LNB)
         RETURN
   10    M = IPMPAR(10)
         EXPARG = 0.99999D0 * (M * LNB)
         RETURN
      END
      DOUBLE PRECISION FUNCTION ESUM(MU, X)
C-----------------------------------------------------------------------
C                    EVALUATION OF EXP(MU + X)
C-----------------------------------------------------------------------
         IMPLICIT NONE
         INTEGER MU
         DOUBLE PRECISION X
C
         DOUBLE PRECISION W
C
         IF (X .GT. 0.0D0) GO TO 10
C
         IF (MU .LT. 0) GO TO 20
         W = MU + X
         IF (W .GT. 0.0D0) GO TO 20
         ESUM = EXP(W)
         RETURN
C
   10    IF (MU .GT. 0) GO TO 20
         W = MU + X
         IF (W .LT. 0.0D0) GO TO 20
         ESUM = EXP(W)
         RETURN
C
   20    W = MU
         ESUM = EXP(W)*EXP(X)
         RETURN
      END
      DOUBLE PRECISION FUNCTION REXP(X)
C-----------------------------------------------------------------------
C            EVALUATION OF THE FUNCTION EXP(X) - 1
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION X
C
         DOUBLE PRECISION P1, P2, Q1, Q2, Q3, Q4, W
C-----------------------
         DATA P1/ .914041914819518D-09/, P2/ .238082361044469D-01/,
     *        Q1/-.499999999085958D+00/, Q2/ .107141568980644D+00/,
     *        Q3/-.119041179760821D-01/, Q4/ .595130811860248D-03/
C-----------------------
         IF (ABS(X) .GT. 0.15D0) GO TO 10
         REXP = X*(((P2*X + P1)*X + 1.0D0)/((((Q4*X + Q3)*X + Q2)*X
     *                    + Q1)*X + 1.0D0))
         RETURN
C
   10    W = EXP(X)
         IF (X .GT. 0.0D0) GO TO 20
         REXP = (W - 0.5D0) - 0.5D0
         RETURN
   20    REXP = W*(0.5D0 + (0.5D0 - 1.0D0/W))
         RETURN
      END
      DOUBLE PRECISION FUNCTION ALNREL(A)
C-----------------------------------------------------------------------
C            EVALUATION OF THE FUNCTION LN(1 + A)
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A
C
         DOUBLE PRECISION P1, P2, P3, Q1, Q2, Q3, T, T2, W, X
C--------------------------
         DATA P1/-.129418923021993D+01/, P2/.405303492862024D+00/,
     *        P3/-.178874546012214D-01/
         DATA Q1/-.162752256355323D+01/, Q2/.747811014037616D+00/,
     *        Q3/-.845104217945565D-01/
C--------------------------
         IF (ABS(A) .GT. 0.375D0) GO TO 10
         T = A/(A + 2.0D0)
         T2 = T*T
         W = (((P3*T2 + P2)*T2 + P1)*T2 + 1.0D0)/
     *       (((Q3*T2 + Q2)*T2 + Q1)*T2 + 1.0D0)
         ALNREL = 2.0D0*T*W
         RETURN
C
   10    X = 1.0D0 + A
         ALNREL = DLOG(X)
         RETURN
      END
      DOUBLE PRECISION FUNCTION RLOG1(X)
C-----------------------------------------------------------------------
C             EVALUATION OF THE FUNCTION X - LN(1 + X)
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION X
C
         DOUBLE PRECISION A, B, H, P0, P1, P2, Q1, Q2, R, T, W, W1
C------------------------
         DATA A/.566749439387324D-01/
         DATA B/.456512608815524D-01/
C------------------------
         DATA P0/ .333333333333333D+00/, P1/-.224696413112536D+00/,
     *        P2/ .620886815375787D-02/
         DATA Q1/-.127408923933623D+01/, Q2/ .354508718369557D+00/
C------------------------
         IF (X .LT. -0.39D0 .OR. X .GT. 0.57D0) GO TO 100
         IF (X .LT. -0.18D0) GO TO 10
         IF (X .GT.  0.18D0) GO TO 20
C
C              ARGUMENT REDUCTION
C
         H = X
         W1 = 0.0D0
         GO TO 30
C
   10    H = X + 0.3D0
         H = H/0.7D0
         W1 = A - H*0.3D0
         GO TO 30
C
   20    H = 0.75D0*X - 0.25D0
         W1 = B + H/3.0D0
C
C               SERIES EXPANSION
C
   30    R = H/(H + 2.0D0)
         T = R*R
         W = ((P2*T + P1)*T + P0)/((Q2*T + Q1)*T + 1.0D0)
         RLOG1 = 2.0D0*T*(1.0D0/(1.0D0 - R) - R*W) + W1
         RETURN
C
C
  100    W = (X + 0.5D0) + 0.5D0
         RLOG1 = X - DLOG(W)
         RETURN
      END
      DOUBLE PRECISION FUNCTION CUSTOM_ERF(X)
C-----------------------------------------------------------------------
C             EVALUATION OF THE REAL ERROR FUNCTION
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION X
C
         DOUBLE PRECISION AX, BOT, C, T, TOP, X2
C
         DOUBLE PRECISION A(5), B(3), P(8), Q(8), R(5), S(4)
C-------------------------
         DATA C /.564189583547756D0/
C-------------------------
         DATA A(1) /.771058495001320D-04/, A(2)/-.133733772997339D-02/,
     *        A(3) /.323076579225834D-01/, A(4) /.479137145607681D-01/,
     *        A(5) /.128379167095513D+00/
         DATA B(1) /.301048631703895D-02/, B(2) /.538971687740286D-01/,
     *        B(3) /.375795757275549D+00/
C-------------------------
         DATA P(1)/-1.36864857382717D-07/, P(2) /5.64195517478974D-01/,
     *        P(3) /7.21175825088309D+00/, P(4) /4.31622272220567D+01/,
     *        P(5) /1.52989285046940D+02/, P(6) /3.39320816734344D+02/,
     *        P(7) /4.51918953711873D+02/, P(8) /3.00459261020162D+02/
         DATA Q(1) /1.00000000000000D+00/, Q(2) /1.27827273196294D+01/,
     *        Q(3) /7.70001529352295D+01/, Q(4) /2.77585444743988D+02/,
     *        Q(5) /6.38980264465631D+02/, Q(6) /9.31354094850610D+02/,
     *        Q(7) /7.90950925327898D+02/, Q(8) /3.00459260956983D+02/
C-------------------------
         DATA R(1) /2.10144126479064D+00/, R(2) /2.62370141675169D+01/,
     *        R(3) /2.13688200555087D+01/, R(4) /4.65807828718470D+00/,
     *        R(5) /2.82094791773523D-01/
         DATA S(1) /9.41537750555460D+01/, S(2) /1.87114811799590D+02/,
     *        S(3) /9.90191814623914D+01/, S(4) /1.80124575948747D+01/
C-------------------------
         AX = ABS(X)
         IF (AX .GT. 0.5D0) GO TO 10
         T = X*X
         TOP = ((((A(1)*T + A(2))*T + A(3))*T + A(4))*T + A(5)) + 1.0D0
         BOT = ((B(1)*T + B(2))*T + B(3))*T + 1.0D0
         CUSTOM_ERF = X*(TOP/BOT)
         RETURN
C
   10    IF (AX .GT. 4.0D0) GO TO 20
         TOP = ((((((P(1)*AX + P(2))*AX + P(3))*AX + P(4))*AX + P(5))*AX
     *                       + P(6))*AX + P(7))*AX + P(8)
         BOT = ((((((Q(1)*AX + Q(2))*AX + Q(3))*AX + Q(4))*AX + Q(5))*AX
     *                       + Q(6))*AX + Q(7))*AX + Q(8)
         CUSTOM_ERF = 0.5D0 + (0.5D0 - EXP(-X*X)*TOP/BOT)
         IF (X .LT. 0.0D0) CUSTOM_ERF = -CUSTOM_ERF
         RETURN
C
   20    IF (AX .GE. 5.8D0) GO TO 30
         X2 = X*X
         T = 1.0D0/X2
         TOP = (((R(1)*T + R(2))*T + R(3))*T + R(4))*T + R(5)
         BOT = (((S(1)*T + S(2))*T + S(3))*T + S(4))*T + 1.0D0
         CUSTOM_ERF = (C - TOP/(X2*BOT)) / AX
         CUSTOM_ERF = 0.5D0 + (0.5D0 - EXP(-X2)*CUSTOM_ERF)
         IF (X .LT. 0.0D0) CUSTOM_ERF = -CUSTOM_ERF
         RETURN
C
   30    CUSTOM_ERF = SIGN(1.0D0,X)
         RETURN
      END
      DOUBLE PRECISION FUNCTION ERFC1(IND, X)
C-----------------------------------------------------------------------
C         EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION
C
C          ERFC1(IND,X) = ERFC(X)            IF IND = 0
C          ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
C-----------------------------------------------------------------------
         IMPLICIT NONE
         INTEGER IND
         DOUBLE PRECISION X
C
         DOUBLE PRECISION AX, BOT, C, E, T, TOP, W
C
         DOUBLE PRECISION A(5), B(3), P(8), Q(8), R(5), S(4)
C
         DOUBLE PRECISION EXPARG
C-------------------------
         DATA C /.564189583547756D0/
C-------------------------
         DATA A(1) /.771058495001320D-04/, A(2)/-.133733772997339D-02/,
     *        A(3) /.323076579225834D-01/, A(4) /.479137145607681D-01/,
     *        A(5) /.128379167095513D+00/
         DATA B(1) /.301048631703895D-02/, B(2) /.538971687740286D-01/,
     *        B(3) /.375795757275549D+00/
C-------------------------
         DATA P(1)/-1.36864857382717D-07/, P(2) /5.64195517478974D-01/,
     *        P(3) /7.21175825088309D+00/, P(4) /4.31622272220567D+01/,
     *        P(5) /1.52989285046940D+02/, P(6) /3.39320816734344D+02/,
     *        P(7) /4.51918953711873D+02/, P(8) /3.00459261020162D+02/
         DATA Q(1) /1.00000000000000D+00/, Q(2) /1.27827273196294D+01/,
     *        Q(3) /7.70001529352295D+01/, Q(4) /2.77585444743988D+02/,
     *        Q(5) /6.38980264465631D+02/, Q(6) /9.31354094850610D+02/,
     *        Q(7) /7.90950925327898D+02/, Q(8) /3.00459260956983D+02/
C-------------------------
         DATA R(1) /2.10144126479064D+00/, R(2) /2.62370141675169D+01/,
     *        R(3) /2.13688200555087D+01/, R(4) /4.65807828718470D+00/,
     *        R(5) /2.82094791773523D-01/
         DATA S(1) /9.41537750555460D+01/, S(2) /1.87114811799590D+02/,
     *        S(3) /9.90191814623914D+01/, S(4) /1.80124575948747D+01/
C-------------------------
C
C                     ABS(X) .LE. 0.5
C
         AX = ABS(X)
         IF (AX .GT. 0.5D0) GO TO 10
         T = X*X
         TOP = ((((A(1)*T + A(2))*T + A(3))*T + A(4))*T + A(5)) + 1.0D0
         BOT = ((B(1)*T + B(2))*T + B(3))*T + 1.0D0
         ERFC1 = 0.5D0 + (0.5D0 - X*(TOP/BOT))
         IF (IND .NE. 0) ERFC1 = EXP(T) * ERFC1
         RETURN
C
C                  0.5 .LT. ABS(X) .LE. 4
C
   10    IF (AX .GT. 4.0D0) GO TO 20
         TOP = ((((((P(1)*AX + P(2))*AX + P(3))*AX + P(4))*AX + P(5))*AX
     *                       + P(6))*AX + P(7))*AX + P(8)
         BOT = ((((((Q(1)*AX + Q(2))*AX + Q(3))*AX + Q(4))*AX + Q(5))*AX
     *                       + Q(6))*AX + Q(7))*AX + Q(8)
         ERFC1 = TOP/BOT
         GO TO 40
C
C                      ABS(X) .GT. 4
C
   20    IF (X .LE. -5.6D0) GO TO 50
         IF (IND .NE. 0) GO TO 30
         IF (X .GT. 100.0D0) GO TO 60
         IF (X*X .GT. -EXPARG(1)) GO TO 60
C
   30    T = (1.0D0/X)**2
         TOP = (((R(1)*T + R(2))*T + R(3))*T + R(4))*T + R(5)
         BOT = (((S(1)*T + S(2))*T + S(3))*T + S(4))*T + 1.0D0
         ERFC1 = (C - T*TOP/BOT)/AX
C
C                      FINAL ASSEMBLY
C
   40    IF (IND .EQ. 0) GO TO 41
         IF (X .LT. 0.0D0) ERFC1 = 2.0D0*EXP(X*X) - ERFC1
         RETURN
   41    W = X*X
         T = W
         E = W - T
         ERFC1 = ((0.5D0 + (0.5D0 - E)) * EXP(-T)) * ERFC1
         IF (X .LT. 0.0D0) ERFC1 = 2.0D0 - ERFC1
         RETURN
C
C             LIMIT VALUE FOR LARGE NEGATIVE X
C
   50    ERFC1 = 2.0D0
         IF (IND .NE. 0) ERFC1 = 2.0D0*EXP(X*X)
         RETURN
C
C             LIMIT VALUE FOR LARGE POSITIVE X
C                       WHEN IND = 0
C
   60    ERFC1 = 0.0D0
         RETURN
      END
      DOUBLE PRECISION FUNCTION GAM1(A)
C     ------------------------------------------------------------------
C     COMPUTATION OF 1/GAMMA(A+1) - 1  FOR -0.5 .LE. A .LE. 1.5
C     ------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A
C
         DOUBLE PRECISION BOT, D, S1, S2, T, TOP, W
         DOUBLE PRECISION P(7), Q(5), R(9)
C     -------------------
         DATA P(1)/ .577215664901533D+00/, P(2)/-.409078193005776D+00/,
     *        P(3)/-.230975380857675D+00/, P(4)/ .597275330452234D-01/,
     *        P(5)/ .766968181649490D-02/, P(6)/-.514889771323592D-02/,
     *        P(7)/ .589597428611429D-03/
C     -------------------
         DATA Q(1)/ .100000000000000D+01/, Q(2)/ .427569613095214D+00/,
     *        Q(3)/ .158451672430138D+00/, Q(4)/ .261132021441447D-01/,
     *        Q(5)/ .423244297896961D-02/
C     -------------------
         DATA R(1)/-.422784335098468D+00/, R(2)/-.771330383816272D+00/,
     *        R(3)/-.244757765222226D+00/, R(4)/ .118378989872749D+00/,
     *        R(5)/ .930357293360349D-03/, R(6)/-.118290993445146D-01/,
     *        R(7)/ .223047661158249D-02/, R(8)/ .266505979058923D-03/,
     *        R(9)/-.132674909766242D-03/
C     -------------------
         DATA S1  / .273076135303957D+00/, S2  / .559398236957378D-01/
C     -------------------
         T = A
         D = A - 0.5D0
         IF (D .GT. 0.0D0) T = D - 0.5D0
C        IF (T) 30,10,20
         IF (T .LT. 0.0D0) GO TO 30
         IF (T .GT. 0.0D0) GO TO 20
C
         GAM1 = 0.0D0
         RETURN
C
   20    TOP = (((((P(7)*T + P(6))*T + P(5))*T + P(4))*T + P(3))*T
     *                     + P(2))*T + P(1)
         BOT = (((Q(5)*T + Q(4))*T + Q(3))*T + Q(2))*T + 1.0D0
         W = TOP/BOT
         IF (D .GT. 0.0D0) GO TO 21
         GAM1 = A*W
         RETURN
   21    GAM1 = (T/A)*((W - 0.5D0) - 0.5D0)
         RETURN
C
   30    TOP = (((((((R(9)*T + R(8))*T + R(7))*T + R(6))*T + R(5))*T
     *                       + R(4))*T + R(3))*T + R(2))*T + R(1)
         BOT = (S2*T + S1)*T + 1.0D0
         W = TOP/BOT
         IF (D .GT. 0.0D0) GO TO 31
         GAM1 = A*((W + 0.5D0) + 0.5D0)
         RETURN
   31    GAM1 = T*W/A
         RETURN
      END
      DOUBLE PRECISION FUNCTION GAMLN1(A)
C-----------------------------------------------------------------------
C     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A
C
         DOUBLE PRECISION P0, P1, P2, P3, P4, P5, P6, Q1, Q2, Q3, Q4,
     +                    Q5,  Q6, R0, R1, R2, R3, R4, R5, S1, S2, S3,
     +                    S4, S5, W, X
C----------------------
         DATA P0/ .577215664901533D+00/, P1/ .844203922187225D+00/,
     *        P2/-.168860593646662D+00/, P3/-.780427615533591D+00/,
     *        P4/-.402055799310489D+00/, P5/-.673562214325671D-01/,
     *        P6/-.271935708322958D-02/
         DATA Q1/ .288743195473681D+01/, Q2/ .312755088914843D+01/,
     *        Q3/ .156875193295039D+01/, Q4/ .361951990101499D+00/,
     *        Q5/ .325038868253937D-01/, Q6/ .667465618796164D-03/
C----------------------
         DATA R0/.422784335098467D+00/,  R1/.848044614534529D+00/,
     *        R2/.565221050691933D+00/,  R3/.156513060486551D+00/,
     *        R4/.170502484022650D-01/,  R5/.497958207639485D-03/
         DATA S1/.124313399877507D+01/,  S2/.548042109832463D+00/,
     *        S3/.101552187439830D+00/,  S4/.713309612391000D-02/,
     *        S5/.116165475989616D-03/
C----------------------
         IF (A .GE. 0.6D0) GO TO 10
         W = ((((((P6*A + P5)*A + P4)*A + P3)*A + P2)*A + P1)*A + P0)/
     *       ((((((Q6*A + Q5)*A + Q4)*A + Q3)*A + Q2)*A + Q1)*A + 1.0D0)
         GAMLN1 = -A*W
         RETURN
C
   10    X = (A - 0.5D0) - 0.5D0
         W = (((((R5*X + R4)*X + R3)*X + R2)*X + R1)*X + R0)/
     *       (((((S5*X + S4)*X + S3)*X + S2)*X + S1)*X + 1.0D0)
         GAMLN1 = X*W
         RETURN
      END
      DOUBLE PRECISION FUNCTION PSI(XX)
C---------------------------------------------------------------------
C
C                 EVALUATION OF THE DIGAMMA FUNCTION
C
C                           -----------
C
C     PSI(XX) IS ASSIGNED THE VALUE 0 WHEN THE DIGAMMA FUNCTION CANNOT
C     BE COMPUTED.
C
C     THE MAIN COMPUTATION INVOLVES EVALUATION OF RATIONAL CHEBYSHEV
C     APPROXIMATIONS PUBLISHED IN MATH. COMP. 27, 123-127(1973) BY
C     CODY, STRECOK AND THACHER.
C
C---------------------------------------------------------------------
C     PSI WAS WRITTEN AT ARGONNE NATIONAL LABORATORY FOR THE FUNPACK
C     PACKAGE OF SPECIAL FUNCTION SUBROUTINES. PSI WAS MODIFIED BY
C     A.H. MORRIS (NSWC).
C---------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION XX
C
         DOUBLE PRECISION AUG, DEN, DX0, PIOV4, SGN, UPPER, W, X, XMAX1,
     +                    XMX0, XSMALL, Z
         INTEGER I, M, N, NQ
C
         DOUBLE PRECISION P1(7), P2(4), Q1(6), Q2(4)
C
         INTEGER IPMPAR
         DOUBLE PRECISION DPMPAR
C---------------------------------------------------------------------
C
C     PIOV4 = PI/4
C     DX0 = ZERO OF PSI TO EXTENDED PRECISION
C
C---------------------------------------------------------------------
         DATA PIOV4/.785398163397448D0/
         DATA DX0/1.461632144968362341262659542325721325D0/
C---------------------------------------------------------------------
C
C     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
C     PSI(X) / (X - X0),  0.5 .LE. X .LE. 3.0
C
C---------------------------------------------------------------------
         DATA P1(1)/.895385022981970D-02/,  P1(2)/.477762828042627D+01/,
     *        P1(3)/.142441585084029D+03/,  P1(4)/.118645200713425D+04/,
     *        P1(5)/.363351846806499D+04/,  P1(6)/.413810161269013D+04/,
     *        P1(7)/.130560269827897D+04/
         DATA Q1(1)/.448452573429826D+02/,  Q1(2)/.520752771467162D+03/,
     *        Q1(3)/.221000799247830D+04/,  Q1(4)/.364127349079381D+04/,
     *        Q1(5)/.190831076596300D+04/,  Q1(6)/.691091682714533D-05/
C---------------------------------------------------------------------
C
C     COEFFICIENTS FOR RATIONAL APPROXIMATION OF
C     PSI(X) - LN(X) + 1 / (2*X),  X .GT. 3.0
C
C---------------------------------------------------------------------
         DATA P2(1)/-.212940445131011D+01/,P2(2)/-.701677227766759D+01/,
     *        P2(3)/-.448616543918019D+01/,P2(4)/-.648157123766197D+00/
         DATA Q2(1)/ .322703493791143D+02/,Q2(2)/ .892920700481861D+02/,
     *        Q2(3)/ .546117738103215D+02/,Q2(4)/ .777788548522962D+01/
C---------------------------------------------------------------------
C
C     MACHINE DEPENDENT CONSTANTS ...
C
C        XMAX1  = THE SMALLEST POSITIVE FLOATING POINT CONSTANT
C                 WITH ENTIRELY INTEGER REPRESENTATION. ALSO USED
C                 AS NEGATIVE OF LOWER BOUND ON ACCEPTABLE NEGATIVE
C                 ARGUMENTS AND AS THE POSITIVE ARGUMENT BEYOND WHICH
C                 PSI MAY BE REPRESENTED AS DLOG(X).
C
C        XSMALL = ABSOLUTE ARGUMENT BELOW WHICH PI*COTAN(PI*X)
C                 MAY BE REPRESENTED BY 1/X.
C
C---------------------------------------------------------------------
         XMAX1 = IPMPAR(3)
         XMAX1 = DMIN1(XMAX1, 1.0D0/DPMPAR(1))
         XSMALL = 1.0D-9
C---------------------------------------------------------------------
         X = XX
         AUG = 0.0D0
         IF (X .GE. 0.5D0) GO TO 200
C---------------------------------------------------------------------
C     X .LT. 0.5,  USE REFLECTION FORMULA
C     PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
C---------------------------------------------------------------------
         IF (ABS(X) .GT. XSMALL) GO TO 100
         IF (X .EQ. 0.0D0) GO TO 400
C---------------------------------------------------------------------
C     0 .LT. ABS(X) .LE. XSMALL.  USE 1/X AS A SUBSTITUTE
C     FOR  PI*COTAN(PI*X)
C---------------------------------------------------------------------
         AUG = -1.0D0 / X
         GO TO 150
C---------------------------------------------------------------------
C     REDUCTION OF ARGUMENT FOR COTAN
C---------------------------------------------------------------------
  100    W = - X
         SGN = PIOV4
         IF (W .GT. 0.0D0) GO TO 120
         W = - W
         SGN = -SGN
C---------------------------------------------------------------------
C     MAKE AN ERROR EXIT IF X .LE. -XMAX1
C---------------------------------------------------------------------
  120    IF (W .GE. XMAX1) GO TO 400
         NQ = INT(W)
         W = W - DBLE(NQ)
         NQ = INT(W*4.0D0)
         W = 4.0D0 * (W - DBLE(NQ) * .25D0)
C---------------------------------------------------------------------
C     W IS NOW RELATED TO THE FRACTIONAL PART OF  4.0 * X.
C     ADJUST ARGUMENT TO CORRESPOND TO VALUES IN FIRST
C     QUADRANT AND DETERMINE SIGN
C---------------------------------------------------------------------
         N = NQ / 2
         IF ((N+N) .NE. NQ) W = 1.0D0 - W
         Z = PIOV4 * W
         M = N / 2
         IF ((M+M) .NE. N) SGN = - SGN
C---------------------------------------------------------------------
C     DETERMINE FINAL VALUE FOR  -PI*COTAN(PI*X)
C---------------------------------------------------------------------
         N = (NQ + 1) / 2
         M = N / 2
         M = M + M
         IF (M .NE. N) GO TO 140
C---------------------------------------------------------------------
C     CHECK FOR SINGULARITY
C---------------------------------------------------------------------
         IF (Z .EQ. 0.0D0) GO TO 400
C---------------------------------------------------------------------
C     USE COS/SIN AS A SUBSTITUTE FOR COTAN, AND
C     SIN/COS AS A SUBSTITUTE FOR TAN
C---------------------------------------------------------------------
         AUG = SGN * ((COS(Z) / SIN(Z)) * 4.0D0)
         GO TO 150
  140    AUG = SGN * ((SIN(Z) / COS(Z)) * 4.0D0)
  150    X = 1.0D0 - X
  200    IF (X .GT. 3.0D0) GO TO 300
C---------------------------------------------------------------------
C     0.5 .LE. X .LE. 3.0
C---------------------------------------------------------------------
         DEN = X
         UPPER = P1(1) * X
C
         DO 210 I = 1, 5
            DEN = (DEN + Q1(I)) * X
            UPPER = (UPPER + P1(I+1)) * X
  210    CONTINUE
C
         DEN = (UPPER + P1(7)) / (DEN + Q1(6))
         XMX0 = DBLE(X) - DX0
         PSI = DEN * XMX0 + AUG
         RETURN
C---------------------------------------------------------------------
C     IF X .GE. XMAX1, PSI = LN(X)
C---------------------------------------------------------------------
  300    IF (X .GE. XMAX1) GO TO 350
C---------------------------------------------------------------------
C     3.0 .LT. X .LT. XMAX1
C---------------------------------------------------------------------
         W = 1.0D0 / (X * X)
         DEN = W
         UPPER = P2(1) * W
C
         DO 310 I = 1, 3
            DEN = (DEN + Q2(I)) * W
            UPPER = (UPPER + P2(I+1)) * W
  310    CONTINUE
C
         AUG = UPPER / (DEN + Q2(4)) - 0.5D0 / X + AUG
  350    PSI = AUG + DLOG(X)
         RETURN
C---------------------------------------------------------------------
C     ERROR RETURN
C---------------------------------------------------------------------
  400    PSI = 0.0D0
         RETURN
      END
      DOUBLE PRECISION FUNCTION BETALN(A0, B0)
C-----------------------------------------------------------------------
C     EVALUATION OF THE LOGARITHM OF THE BETA FUNCTION
C-----------------------------------------------------------------------
C     E = 0.5*LN(2*PI)
C--------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A0, B0
C
         DOUBLE PRECISION A, B, C, E, H, U, V, W, Z
         INTEGER I, N
C
         DOUBLE PRECISION ALGDIV, ALNREL, BCORR, GAMLN, GSUMLN
C--------------------------
         DATA E /.918938533204673D0/
C--------------------------
         A = DMIN1(A0,B0)
         B = DMAX1(A0,B0)
         IF (A .GE. 8.0D0) GO TO 60
         IF (A .GE. 1.0D0) GO TO 20
C-----------------------------------------------------------------------
C                   PROCEDURE WHEN A .LT. 1
C-----------------------------------------------------------------------
         IF (B .GE. 8.0D0) GO TO 10
         BETALN = GAMLN(A) + (GAMLN(B) - GAMLN(A + B))
         RETURN
   10    BETALN = GAMLN(A) + ALGDIV(A,B)
         RETURN
C-----------------------------------------------------------------------
C                PROCEDURE WHEN 1 .LE. A .LT. 8
C-----------------------------------------------------------------------
   20    IF (A .GT. 2.0D0) GO TO 30
         IF (B .GT. 2.0D0) GO TO 21
         BETALN = GAMLN(A) + GAMLN(B) - GSUMLN(A,B)
         RETURN
   21    W = 0.0D0
         IF (B .LT. 8.0D0) GO TO 40
         BETALN = GAMLN(A) + ALGDIV(A,B)
         RETURN
C
C                REDUCTION OF A WHEN B .LE. 1000
C
   30    IF (B .GT. 1000.0D0) GO TO 50
         N = INT(A) - 1
         W = 1.0D0
         DO 31 I = 1,N
            A = A - 1.0D0
            H = A/B
            W = W * (H/(1.0D0 + H))
   31    CONTINUE
         W = DLOG(W)
         IF (B .LT. 8.0D0) GO TO 40
         BETALN = W + GAMLN(A) + ALGDIV(A,B)
         RETURN
C
C                 REDUCTION OF B WHEN B .LT. 8
C
   40    N = INT(B) - 1
         Z = 1.0D0
         DO 41 I = 1,N
            B = B - 1.0D0
            Z = Z * (B/(A + B))
   41    CONTINUE
         BETALN = W + DLOG(Z) + (GAMLN(A) + (GAMLN(B) - GSUMLN(A,B)))
         RETURN
C
C                REDUCTION OF A WHEN B .GT. 1000
C
   50    N = INT(A) - 1
         W = 1.0D0
         DO 51 I = 1,N
            A = A - 1.0D0
            W = W * (A/(1.0D0 + A/B))
   51    CONTINUE
         BETALN = (DLOG(W) - N*DLOG(B)) + (GAMLN(A) + ALGDIV(A,B))
         RETURN
C-----------------------------------------------------------------------
C                   PROCEDURE WHEN A .GE. 8
C-----------------------------------------------------------------------
   60    W = BCORR(A,B)
         H = A/B
         C = H/(1.0D0 + H)
         U = -(A - 0.5D0)*DLOG(C)
         V = B*ALNREL(H)
         IF (U .LE. V) GO TO 61
         BETALN = (((-0.5D0*DLOG(B) + E) + W) - V) - U
         RETURN
   61    BETALN = (((-0.5D0*DLOG(B) + E) + W) - U) - V
         RETURN
      END
      DOUBLE PRECISION FUNCTION GSUMLN(A, B)
C-----------------------------------------------------------------------
C          EVALUATION OF THE FUNCTION LN(GAMMA(A + B))
C          FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B
C
         DOUBLE PRECISION X
C
         DOUBLE PRECISION ALNREL, GAMLN1
C
         X = A + B - 2.0D0
         IF (X .GT. 0.25D0) GO TO 10
         GSUMLN = GAMLN1(1.0D0 + X)
         RETURN
   10    IF (X .GT. 1.25D0) GO TO 20
         GSUMLN = GAMLN1(X) + ALNREL(X)
         RETURN
   20    GSUMLN = GAMLN1(X - 1.0D0) + DLOG(X*(1.0D0 + X))
         RETURN
      END
      DOUBLE PRECISION FUNCTION BCORR(A0, B0)
C-----------------------------------------------------------------------
C
C     EVALUATION OF  DEL(A0) + DEL(B0) - DEL(A0 + B0)  WHERE
C     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
C     IT IS ASSUMED THAT A0 .GE. 8 AND B0 .GE. 8.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A0, B0
C
         DOUBLE PRECISION A, B, C, C0, C1, C2, C3, C4, C5, H, S11, S3,
     +                    S5, S7, S9, T, W, X, x2
C------------------------
         DATA C0/.833333333333333D-01/, C1/-.277777777760991D-02/,
     *        C2/.793650666825390D-03/, C3/-.595202931351870D-03/,
     *        C4/.837308034031215D-03/, C5/-.165322962780713D-02/
C------------------------
         A = DMIN1(A0, B0)
         B = DMAX1(A0, B0)
C
         H = A/B
         C = H/(1.0D0 + H)
         X = 1.0D0/(1.0D0 + H)
         X2 = X*X
C
C                SET SN = (1 - X**N)/(1 - X)
C
         S3 = 1.0D0 + (X + X2)
         S5 = 1.0D0 + (X + X2*S3)
         S7 = 1.0D0 + (X + X2*S5)
         S9 = 1.0D0 + (X + X2*S7)
         S11 = 1.0D0 + (X + X2*S9)
C
C                SET W = DEL(B) - DEL(A + B)
C
         T = (1.0D0/B)**2
         W = ((((C5*S11*T + C4*S9)*T + C3*S7)*T + C2*S5)*T + C1*S3)*T
         W = W + C0
         W = W*(C/B)
C
C                   COMPUTE  DEL(A) + W
C
         T = (1.0D0/A)**2
         BCORR = (((((C5*T + C4)*T + C3)*T + C2)*T + C1)*T + C0)/A + W
         RETURN
      END
      DOUBLE PRECISION FUNCTION ALGDIV(A, B)
C-----------------------------------------------------------------------
C
C     COMPUTATION OF LN(GAMMA(B)/GAMMA(A+B)) WHEN B .GE. 8
C
C                         --------
C
C     IN THIS ALGORITHM, DEL(X) IS THE FUNCTION DEFINED BY
C     LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).
C
C-----------------------------------------------------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A, B
C
         DOUBLE PRECISION C, C0, C1, C2, C3, C4, C5, D, H, S11, S3, S5,
     +                    S7, S9, T, U, V, W, X, X2
C
         DOUBLE PRECISION ALNREL
C------------------------
         DATA C0/.833333333333333D-01/, C1/-.277777777760991D-02/,
     *        C2/.793650666825390D-03/, C3/-.595202931351870D-03/,
     *        C4/.837308034031215D-03/, C5/-.165322962780713D-02/
C------------------------
         IF (A .LE. B) GO TO 10
         H = B/A
         C = 1.0D0/(1.0D0 + H)
         X = H/(1.0D0 + H)
         D = A + (B - 0.5D0)
         GO TO 20
   10    H = A/B
         C = H/(1.0D0 + H)
         X = 1.0D0/(1.0D0 + H)
         D = B + (A - 0.5D0)
C
C                SET SN = (1 - X**N)/(1 - X)
C
   20    X2 = X*X
         S3 = 1.0D0 + (X + X2)
         S5 = 1.0D0 + (X + X2*S3)
         S7 = 1.0D0 + (X + X2*S5)
         S9 = 1.0D0 + (X + X2*S7)
         S11 = 1.0D0 + (X + X2*S9)
C
C                SET W = DEL(B) - DEL(A + B)
C
         T = (1.0D0/B)**2
         W = ((((C5*S11*T + C4*S9)*T + C3*S7)*T + C2*S5)*T + C1*S3)*T
         W = W + C0
         W = W*(C/B)
C
C                    COMBINE THE RESULTS
C
         U = D*ALNREL(A/B)
         V = A*(DLOG(B) - 1.0D0)
         IF (U .LE. V) GO TO 30
         ALGDIV = (W - V) - U
         RETURN
   30    ALGDIV = (W - U) - V
         RETURN
      END
      DOUBLE PRECISION FUNCTION GAMLN(A)
C-----------------------------------------------------------------------
C            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
C-----------------------------------------------------------------------
C     WRITTEN BY ALFRED H. MORRIS
C          NAVAL SURFACE WARFARE CENTER
C          DAHLGREN, VIRGINIA
C--------------------------
C     D = 0.5*(LN(2*PI) - 1)
C--------------------------
         IMPLICIT NONE
         DOUBLE PRECISION A
C
         DOUBLE PRECISION C0, C1, C2, C3, C4, C5, D, T, W
         INTEGER I, N
C
         DOUBLE PRECISION GAMLN1
C--------------------------
         DATA D/.418938533204673D0/
C--------------------------
         DATA C0/.833333333333333D-01/, C1/-.277777777760991D-02/,
     *        C2/.793650666825390D-03/, C3/-.595202931351870D-03/,
     *        C4/.837308034031215D-03/, C5/-.165322962780713D-02/
C-----------------------------------------------------------------------
         IF (A .GT. 0.8D0) GO TO 10
         GAMLN = GAMLN1(A) - DLOG(A)
         RETURN
   10    IF (A .GT. 2.25D0) GO TO 20
         T = (A - 0.5D0) - 0.5D0
         GAMLN = GAMLN1(T)
         RETURN
C
   20    IF (A .GE. 10.0D0) GO TO 30
         N = INT(A - 1.25D0)
         T = A
         W = 1.0D0
         DO I = 1,N
            T = T - 1.0D0
            W = T*W
         END DO
         GAMLN = GAMLN1(T - 1.0D0) + DLOG(W)
         RETURN
C
   30    T = (1.0D0/A)**2
         W = (((((C5*T + C4)*T + C3)*T + C2)*T + C1)*T + C0)/A
         GAMLN = (D + W) + (A - 0.5D0)*(DLOG(A) - 1.0D0)
      END
      INTEGER FUNCTION IPMPAR(I)
C-----------------------------------------------------------------------
C
C     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
C     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
C     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
C
C  INTEGERS.
C
C     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
C
C               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
C
C     IPMPAR(1) = A, THE BASE.
C
C     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
C
C     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
C     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
C     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
C
C               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
C
C               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
C               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
C
C     IPMPAR(4) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
C
C-----------------------------------------------------------------------
C
C     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
C     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
C     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
C     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
C
C-----------------------------------------------------------------------
         INTEGER I
         INTEGER IMACH(10)
C
         DATA IMACH( 1) /     2 /
         DATA IMACH( 2) /    31 /
         DATA IMACH( 3) / 2147483647 /
         DATA IMACH( 4) /     2 /
         DATA IMACH( 5) /    24 /
         DATA IMACH( 6) /  -125 /
         DATA IMACH( 7) /   128 /
         DATA IMACH( 8) /    53 /
         DATA IMACH( 9) / -1021 /
         DATA IMACH(10) /  1024 /
C
         IPMPAR = IMACH(I)
         RETURN
      END
