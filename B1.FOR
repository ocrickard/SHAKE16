C $NOFLOATCALLS
C $NODEBUG
C************************************
      SUBROUTINE AMP( N1,IN,INT,LL,LT,KPL,IDAMP,NA,DF)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C   THIS ROUTINE COMPUTES THE AMPLIFICATION SPECTRUM BETWEEN ANY TWO
C   LAYERS
C
C        N1      = NUMBER OF SOIL LAYERS EXCLUDING ROCK
C        IN      = NUMBER OF SUBLAYER FROM WHICH AMPLIFICATION IS COMP.
C        INT     = SUBLAYER TYPE
C                       0 - OUTCROPPING LAYER
C                       1 - LAYER WITHIN PROFILE
C        LL      = NUMBER OF SUBLAYER TO WHICH AMPLIFICATION IS COMP.
C        LT      = SUBLAYER TYPE
C                       0 - OUTCROPPING LAYER
C                       1 - LAYER WITHIN PROFILE
C        DF      = FREQUENCY STEPS IN AMP. FUNCTION
C        NA      = CURVE NUMBER IN PLOTTING
C        IDAMP   = IDENTIFICATION
C
C   CODED PER B SCHNABEL  FEB. 1971
C   modified to increase number of sublayers to 50
C   February 1991
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      COMPLEX G, V, PLUS, MINUS
      COMPLEX E, F, EE, FF, A, EX, AIN, IPI2,AA
      CHARACTER*60 ABSIS
      CHARACTER*6 ID,IDNT,IDAMP
C
      DIMENSION IDAMP (27,11),T(200)
      COMMON /JOE4/ ST(27,200)
      COMMON /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      COMMON /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
      COMMON /CCG/ ID(27,11)
C
      ABSIS = ' CYCLES/SEC. '
C
      IPI2 = CMPLX(0., 6.283185307)
      FREQ = 0.
      ST(NA,1) = 1.
      DO 19 I = 2,200
      E = 1.
      FF = 1.
      FREQ = FREQ + DF
      A = FREQ*IPI2
      DO 191 K = 1,N1
      IF (K.NE.IN) GO TO 192
      AIN = E + FF
      IF (INT.EQ.0) AIN = 2.*E
  192 IF (K.NE.LL) GO TO 11
      AA    = E + FF
      IF (LT.EQ.0) AA = 2.*E
   11 EX = CEXP(H(K)*A/V(K))
      EE = E*EX
      F = FF/EX
      E = EE*PLUS(K) + MINUS(K)*F
      FF = PLUS(K)*F + MINUS(K)*EE
  191 CONTINUE
      IF (IN.NE.N1+1) GO TO 193
      AIN = E + FF
      IF (INT.EQ.0) AIN = 2.*E
  193 IF (LL   .NE.N1+1) GO TO 21
      AA    = E + FF
      IF (LT.EQ.0)  AA = 2.*E
   21 ST(NA,I) = CABS(AA/AIN)
   19 CONTINUE
      DO 23 I = 1,200
   23 T(I) = DF*FLOAT(I-1)
      AMAX = 0.
      WRITE(6,2)
      DO 22 I = 1,200
      IF (KPL .GE. 2) WRITE(6,1) T(I), ST(NA,I)
      IF (ST(NA,I) .LT.  AMAX) GO TO 22
      TMAX = T(I)
      AMAX = ST(NA,I)
   22 CONTINUE
      IF (NA.LT.9) NA=NA+1
      PERIOD = 1./TMAX
      IF (TMAX.LT. .0001)  WRITE(6,1001) AMAX, TMAX
      IF (TMAX.GT. .0001)  WRITE(6,1001) AMAX, TMAX,PERIOD
      IF (KPL.EQ.0) RETURN
      WRITE(6,1000)
      N = NA-1
      NA = 1
      RETURN
    1 FORMAT(1X,F10.4, 3X, F10.4)
    2 FORMAT(/2X,'FREQUENCY    AMPLITUDE')
 1000 FORMAT(33H1  PLOT OF AMPLIFICATION SPECTRA  /)
 1001 FORMAT(25H MAXIMUM AMPLIFICATION =   F6.2/
     1 25H FOR FREQUENCY         = F6.2, 7H C/SEC. /
     1 25H     PERIOD            =  F6.2, 5H SEC. )
      END
C******************************************
      SUBROUTINE UTPR(KK,DPTH,LS,K2,LH,LT,X,AX,S,INV)
C***********************************************************************
C
C   THIS ROUTINE TRANSFERS THE VALUES IN AX(LH, ) INTO THE TIME DOMAIN
C   IN X( ), TRANSFERS RESULTS TO OUTPUT FILE
C
C        KK      = 5 TABULATE MAX. ACC.
C                  6 PRINT MAX ACC. SEPARATELY
C        DPTH    = DEPTH OF LAYER
C        X( )    = OBJECT MOTION
C        AX(LS, )= COMPUTED MOTION
C        LS      = COMPUTED MOTION NUMBER
C                  0 IF OBJECT MOTION
C        LH      = SUBLAYER NUMBER
C        LT      = SUBLAYER TYPE
C                  0 - OUTCROPPING
C                  1 - INSIDE
C        S,INV  SCRATCH ARRAYS
C
C   CODED  PER B SCHNABEL  OCT. 1970
C   MODIFIED PBS AUG. 1971
C   modified to increase number of layers to 50
C***********************************************************************
      CHARACTER*6 TITLE,IDNT
      COMPLEX SAVE
      COMPLEX X, AX
C
      DIMENSION XR(8)
      DIMENSION X(300),AX(3,270),S(70),INV(70)
      COMMON /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      COMMON /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
C
      FREQ = 0.
      SFX = 0.
      SXX = 0.
C   TRANSFORM VALUES IN X OR IN AX INTO THE TIME DOMAIN
      DO 24  I = 1,MFOLD
      IF (LS.EQ.0) GO TO 241
      SAVE = X(I)
      X(I) = AX(LS,I)
      AX(LS,I) = SAVE
  241 XA = CABS(X(I))
      SXX= SXX + XA*XA
      SFX = SFX + XA*FREQ*XA
      FREQ = FREQ + DF
   24 CONTINUE
      SFX = SFX/SXX
C
      CALL RFSN(X,MX,INV,S,IFERR,-2)
C
      CALL XMX(X,MA,XMAX,NMAX)
      TMAX = DT*FLOAT(NMAX-1)
      XEND = 0.
      N = MA/20
      NN = 9*N
      N = 8*N
      DO 25 I = N,NN
      XABS = REAL(X(I))
      XABS = ABS(XABS)
      IF (XABS.GT.XEND) XEND = XABS
      XABS = AIMAG(X(I))
      XABS = ABS(XABS)
      IF (XABS.GT.XEND) XEND = XABS
   25 CONTINUE
      XEND = XEND/XMAX
C
C   SAVE OUTPUT
      N = 1
      NN = 4
      NCARDS=MA/8
      NC = NCARDS
      IF (K2.EQ.0) NC = 0
C     IF (KK.EQ.5) GO TO 252
      IF (KK.EQ.6) GO TO 252
      IF (LT.EQ.0)  WRITE(6,2000) LH, (IDNT(I),I=1,6)
      IF (LT.EQ.1)  WRITE(6,2002) LH, (IDNT(I),I=1,6)
      WRITE(6, 2005) SFX
      WRITE(6,2003) XMAX, TMAX
  252 IF (KK.EQ.6.AND.LT.EQ.0)  WRITE(6,2001) DPTH,XMAX,TMAX,SFX,XEND,NC
      IF (KK.EQ.6.AND.LT.EQ.1) WRITE(6,2010) DPTH,XMAX,TMAX,SFX, XEND,NC
      IF (K2.EQ.0) GO TO 262
      WRITE(7,2006) XMAX, (TITLE(I),I=1,5)
      IF (LT.EQ.1)  WRITE(7,2002) LH, (IDNT(I),I=1,6)
      IF (LT.EQ.0)  WRITE(7,2000) LH, (IDNT(I),I=1,6)
      DO 26 I = 1,NCARDS
      K = 0
      DO 261 J = N,NN
      K = K+ 1
      XR(K) = REAL(X(J))
      K = K + 1
      XR(K) = AIMAG(X(J))
  261 CONTINUE
      WRITE(7,2009) (XR(J),J=1,8),I
      IF (K2 .EQ. 2) WRITE(6,2019) (XR(J), J = 1,8), I
      NN = 4 + NN
      N = N + 4
   26 CONTINUE
  262 CALL RFFT(X,MX,INV,S,IFERR,2)
      IF (LS.EQ.0)   RETURN
      DO 27 I = 1,MFOLD
      SAVE = AX(LS,I)
      AX(LS,I) = X(I)
   27 X(I) = SAVE
      RETURN
C
 2000 FORMAT(43H  ACCELERATION VALUES AT OUTCROPPING LAYER I3,3H - 6A6)
 2001 FORMAT(5X,6HOUTCR. F15.1,F15.5,2F15.2,F20.3,I20)
 2010 FORMAT(5X,6HWITHIN F15.1,F15.5,2F15.2,F20.3,I20)
 2002 FORMAT(42H  ACCELERATION VALUES AT THE TOP OF LAYER I3,3H - 6A6)
 2003 FORMAT(/15H  MAX. ACC. =  F9.6,11H AT TIME = F6.3, 5H SEC. /)
 2005 FORMAT(/26H   MEAN SQUARE FREQUENCY = F10.2/)
 2006 FORMAT(21X,6HXMAX= F7.4,5A6)
 2008 FORMAT(2X, I5, 5X, 8F15.6)
 2009 FORMAT(8F9.6,I7)
 2019 FORMAT(8F14.6,I10)
      END
C**************************************
      SUBROUTINE REDUCE(IFR,X,AX,LL)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   THIS ROUTINE INCREASES TIME INTERVAL AND REDUCES NUMBER OF VALUES
C
C        IFR   =   DIVIDING FACTOR ON LENGTH OF RECORD
C                  MULTIPLICATION FACTOR ON TIME STEP
C                  MUST BE A POWER OF 2.
C        DT    =   TIMESTEP IN SEC.
C        DF    =   FREQUENCY STEP IN C/SEC.
C        MA    =   NUMBER OF POINTS USED IN FOURIER TRANSFORM
C        X     =   FOURIER TRANSFORM OF OBJECT MOTION
C        AX    =   FOURIER TRANSFORM OF COMPUTED MOTIONS
C
C
C   CODED BY PER B. SCHNABEL DEC. 1970.
C   MODIFIED SEPT. 1971
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      CHARACTER*6 TITLE
      COMPLEX X, AX
C
      DIMENSION X( 68), AX(3, 64), LL(3)
      COMMON /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
      COMMON/FRCUT/ NCUT,NZERO
C
      F1 = .5/DT
      FR = FLOAT(IFR)
      DT = DT*FR
      MA = MA/IFR
      MMA = MMA/IFR
      MA2 = MA + 2
      MFOLD = MA2/2
      N = MFOLD + 1
      DO 12 I = MFOLD,N
      X(I) = 0.
C
      DO 12 L = 1,3
      IF (LL(L).LE.0) GO TO 12
      AX(L,I) = 0.
   12 CONTINUE
      MFOLD = MFOLD + 1
      F2 = .5/DT
      WRITE(6,1000) F1,F2,DT, MA
      FMA = FLOAT(MA)
      MX = (ALOG10(FMA)/ALOG10(2.))-1.
      IF (MA.GT.2**(MX+1)) MX=MX+1
      IF(NCUT.LE.MFOLD) GO TO 15
      NCUT=MFOLD
   15 CONTINUE
 1000 FORMAT(  20H  FREQUENCIES FROM  F6.2, 3H TO F6.2,14H C/SEC ARE REM
     15HOVED /
     216H NEW TIMESTEP = F5.4/19H NUMBER OF VALUES = I5)
      RETURN
      END
C*****************************
      SUBROUTINE INCR(IFR,X,AX)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   THIS ROUTINE INCREASES NUMBER OF POINTS IN THE RECORD
C   BY DECREASING TIMESTEP
C
C        IFR   =   MULTIPLYING FACTOR ON LENGTH OF RECORD
C                  MUST BE A POWER OF 2.
C        DT    =   TIMESTEP IN SEC.
C        DF    =   FREQUENCY STEP IN C/SEC.
C        MA    =   NUMBER OF POINTS USED IN FOURIER TRANSFORM
C        X     =   FOURIER TRANSFORM OF OBJECT MOTION
C        AX    =   FOURIER TRANSFORM OF COMPUTED MOTIONS
C
C
C   CODED BY PER B. SCHNABEL DEC. 1970.
C   MODIFIED OCT.  1971
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      COMPLEX X, AX
      CHARACTER*6 TITLE
C
      DIMENSION X( 68), AX(3, 64)
      COMMON /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
C
      F1 = .5/DT
      FR = FLOAT(IFR)
      DT = DT/FR
      N = MFOLD- 1
      MA = MA*IFR
      MMA = MMA*IFR
      MA2 = MA + 2
      MFOLD = MA2/2
      MFOLD = MFOLD + 1
      DO 10 I = N, MFOLD
      X(I) = 0.
      DO 10 L = 1,3
   10 AX(L,I) = 0.
      F2 = .5/DT
      WRITE(6,1000) F1,F2,DT, MA
      FMA = FLOAT(MA)
      MX = (ALOG10(FMA)/ALOG10(2.))-1.
      IF (MA.GT.2**(MX+1)) MX=MX+1
 1000 FORMAT(27H    FREQUENCIES ADDED FROM  F6.2,3H TO F6.2/
     216H NEW TIME STEP = F5.4/19H NUMBER OF VALUES = I5/)
      RETURN
      END
C***********************************************************
      SUBROUTINE MOTION(N1,IN,INT,LL,LT, X,AX)
***********************************************************
C   THIS ROUTINE CALCULATES THE MOTION IN ANY TWO SOIL LAYERS OR IN
C   ROCK FROM MOTION GIVEN IN ANY LAYER OR IN ROCK
C
C        N1      = NUMBER OF SOIL LAYERS EXCLUDING ROCK
C        IN      = NUMBER OF LAYER WHERE OBJECT MOTION IS GIVEN
C        INT     = MOTION TYPE
C                  IF EQUEL 0     OUTCROPPING LAYER
C        LL()    = NUMBER OF LAYERS WHERE OUTPUT MOTION IS WANTED
C                  MAX  3 LAYERS
C        LT()    = MOTION TYPE
C                  0 - OUTCROPPING LAYER
C                  1 - LAYER WITHIN PROFILE
C        X()     = OBJECT MOTION
C        AX()    = OUTPUT MOTION
C
C   CODED BY PER B SCHNABEL  OCT 1970
C   modified to increase the number of layers to 50
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      INTEGER LL(3), LT(3)
      CHARACTER*6 TITLE,IDNT
      COMPLEX AA(3)
      COMPLEX X, AX
      COMPLEX G, V, PLUS, MINUS
      COMPLEX E, F, EE, FF, A, EX, AIN, IPI2
C
      DIMENSION X(300),AX(3,270),S(70)
      COMMON /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
      COMMON /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      COMMON /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
      COMMON/FRCUT/ NCUT,NZERO
C
      IPI2 = CMPLX(0., 6.283185307)
      DO 20 L = 1,3
      IF (LL(L) .GT. 0) AX(L,1) = X(1)
      IF(NCUT.EQ.MFOLD) GO TO 20
      DO 30 I=NZERO,MFOLD
      AX(L,I)=CMPLX(0.,0.)
   30 CONTINUE
   20 CONTINUE
      FREQ = 0.
      DO 19 I=2,NCUT
      E = 1.
      FF = 1.
      FREQ = FREQ + DF
      A = FREQ*IPI2
      DO 191 K = 1,N1
      IF (K.NE.IN) GO TO 192
      AIN = E + FF
      IF (INT.EQ.0) AIN = 2.*E
C   FIND SUBLAYER WHERE MOTION IS WANTED
  192 DO 11 L = 1,3
      IF (K.NE.LL(L)) GO TO 11
C   AMPLIFICATION FACTOR FOR SUBLAYER WITHIN PROFILE
      AA(L) = E + FF
C   AMPLIFICATION FACTOR FOR OUTCROPPING SUBLAYER
      IF (LT(L).EQ.0)  AA(L) = 2.*E
   11 CONTINUE
      EX = CEXP(H(K)*A/V(K))
      EE = E*EX
      F = FF/EX
      E  = EE*PLUS(K) + MINUS(K)*F
      FF = PLUS(K)*F + MINUS(K)*EE
  191 CONTINUE
      IF (IN.NE.N1+1) GO TO 193
      AIN = E + FF
      IF (INT.EQ.0) AIN = 2.*E
  193 DO 21 L = 1,3
      IF (LL(L).NE.N1+1) GO TO 21
      AA(L) = E + FF
      IF (LT(L).EQ.0) AA(L) = 2.*E
   21 CONTINUE
      DO 23 L = 1,3
      IF (LL(L) .GT. 0) AX(L,I) = X(I)*AA(L)/AIN
   23 CONTINUE
   19 CONTINUE
      RETURN
      END
C*************************************************
      SUBROUTINE CXSOIL(N1)
C***********************************************************************
C
C   THIS ROUTINE CALCULATES THE COMPLEX SOIL PROPERTIES AND TRANSFER
C   FUNCTIONS FOR THE LAYERS
C
C        N1      = NUMBER OF SOIL LAYERS
C        BL      = RATIO OF CRITICAL DAMPING
C        GL      = SHEAR MODULUS
C        R       = DENSITY
C        G       = COMPLEX SHEAR MODULUS
C        V       = COMPLEX SHEAR WAVE VELOCITY
C        PLUS    = COMPLEX TRANSFER FUNCTION
C        MINUS   = COMPLEX TRANSFER FUNCTION
C
C   CODED BY PER B SCHNABEL OCT 1971
C
C***********************************************************************
C
      COMPLEX G, V, PLUS, MINUS, MU
      CHARACTER*6 IDNT
      COMMON /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      COMMON /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
C
      N = N1 + 1
      DO 1 I = 1,N
      GIMAG=2.*BL(I)*GL(I)*SQRT(1.-BL(I)*BL(I))
      GREAL=GL(I)*(1.-2.*BL(I)*BL(I))
      G(I)=CMPLX(GREAL,GIMAG)
      V(I) = CSQRT(G(I)/R(I))
    1 CONTINUE
      DO 2 I = 1,N1
      J = I + 1
      MU = CSQRT(R(I)/R(J)*G(I)/G(J))
      PLUS(I) = (1. + MU)/2.
      MINUS(I)= (1. - MU)/2.
    2 CONTINUE
      RETURN
      END
C*********************************************
      SUBROUTINE STRAIN( LL, LGS, LPCH, LPL,LNV,X,AX,AA,N1,S,INV)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C   THIS SUBROUTINE COMPUTES STRAIN AND/OR STRESS TIME-HISTORY AT THE
C   TOP OF ANY LAYER FOR ACCELERATION HISTORY KNOWN IN ANY LAYER
C   TWO RESPONSE HISTORIES ARE COMPUTED IN ONE RUN
C
C        LL    =   SUBLAYER NUMBER WHERE RESPONSE IS TO BE COMPUTED
C        LGS   =   SWITCH FOR STRESS OR STRAIN
C        LPCH  =   SWITCH FOR SAVING OUTPUT
C        LPL   =   SWITCH FOR PLOT
C        X     =   FOURIER TRANSFORM OF OBJECT MOTION
C        AX(1, )   FOURIER TRANSFORM OF SURFACE MOTION
C        AX(2, )   FOURIER TRANSFORM OF FIRST COMPUTED RESPONSE
C        AX(3, )   FOURIER TRANSFORM OF SECOND RESPONSE
C        AA(1, )   TIME HISTORY OF FIRST RESPONSE
C        AA(2, )   TIME HISTORY OF SECOND RESPONSE
C
C   CODED BY PER B. SCHNABEL  JULY 1971
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      INTEGER TP
      CHARACTER*6 TITLE,IDNT,ID
      CHARACTER*60 ABSIS
      COMPLEX X, AX
      COMPLEX G, V, PLUS, MINUS
      COMPLEX E,F,EE, A,AH,IPI2, AE,AF,EX,AI
C
      DIMENSION AE(2), AF(2)
      DIMENSION X(1),  AX(3,1), AA(2,1), S(1),  INV(1)
      DIMENSION LL(2), LGS(2),  LPCH(2), LPL(2), LNV(2)
C
      COMMON /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      COMMON /SOILB/ FAC(51), WL(51), TP(51), DEPTH(51), WEIGHT(51)
      COMMON /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
      COMMON /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
      COMMON /CCG/ ID(27,11)
      COMMON /FRCUT/ NCUT,NZERO
      COMMON /TIME/ T(9)
C
      ABSIS = ' TIME IN SEC '
      IPI2 = CMPLX(0.,6.283185307)
      GT = 32.2
      AX(2,1) = 0.
      AX(3,1) = 0.
      FREQ = 0.
      AI = GT/IPI2
C
C  STARTING AT THE SURFACE THE STRAIN IS COMPUTED SUCCESSIVELY DOWNWARDS
C   FOR EACH FREQUENCY
      DO 1 I=2,NCUT
      E = AX(1,I)/2.
      F = E
      FREQ = FREQ + DF
      AH = AI/FREQ
      A = FREQ*IPI2
      DO 11  K = 1,N1
      DO 12 L = 1,2
      IF (K.NE.LL(L)) GO TO 12
      AE(L) = E/V(K)
      AF(L) = F/V(K)
   12 CONTINUE
      EX = CEXP(H(K)*A/V(K))
      E = E*EX
      F = F/EX
      EE = E*PLUS(K) + MINUS(K)*F
      F  = F*PLUS(K) + MINUS(K)*E
      E = EE
   11 CONTINUE
      DO 13 L = 1,2
      IF (LL(L).NE.N1+1) GO TO 13
      AE(L) = E/V(N1+1)
      AF(L) = F/V(N1+1)
   13 CONTINUE
      DO 14 L = 1,2
      IF (LL(L).GT.0) AX(L+1,I) = (AE(L) -AF(L))*AH
   14 CONTINUE
    1 CONTINUE
      DO 2 I = 1,MFOLD
    2 AX(1,I) = X(I)
      DO 3 L = 1,2
      IF (LL(L).EQ.0) GO TO 3
      X(1) = 0.
      DO 31 I=2,NCUT
   31 X(I) = AX(L+1,I)
      IF(NCUT.EQ.MFOLD) GO TO 33
      DO 34 II=NZERO,MFOLD
      X(II)=CMPLX(0.,0.)
   34 CONTINUE
   33 CONTINUE
      CALL RFSN(X,MX,INV,S,IFERR,-2)
      DO 32 I =1,MFOLD
      AA(L,2*I-1) =REAL(X(I))*100.
   32 AA(L,2*I) = AIMAG(X(I))*100.
    3 CONTINUE
C
      DO 4 I = 1,MFOLD
    4 X(I) = AX(1,I)
C   COMPUTE STRESS IF WANTED AND SAVE COMPUTED RESPONSES
      DO 5 L = 1,2
      IF (LL(L) .EQ. 0) GO TO 5
      NVAL = LNV(L)
      IF (NVAL.LE.0) NVAL = MMA
      IF (NVAL.GT.MA)  NVAL = MA
      IF (NVAL.GT.2049) NVAL = 2049
      DO 51 I = 1,5
   51 ID(L,I) = TITLE(I)
      N = LL(L)
      ID(L,6) = 'STRAIN'
      IF (LGS(L) .EQ.0) GO TO 53
      ID(L,6) = 'STRESS'
      DO 52 I = 1,NVAL
   52 AA(L,I) = GL(N)*AA(L,I)/100.
   53 IF (LPCH(L).EQ.0) GO TO 54
      WRITE(7,2000) (ID(L,I), I=1,11),N
      N = 1
      NCARDS = NVAL/8
      DO 55 K = 1,NCARDS
      NN = N + 7
      WRITE(7,2001) (AA(L,I), I = N,NN), K
   55 N = N + 8
   54 IF (LPL(L).EQ.0) GO TO 5
      N = 0
      NSKIP = 1
      DO 56 I = 1,NVAL,NSKIP
      N = N + 1
      IF (NSKIP.GT.1) AA(L,N) = AA(L,I)
   56 T(N) = DT*FLOAT(I-1)
      IF  (LGS(L).EQ.0)  WRITE(6,2002)
      IF (LGS(L).EQ.1)  WRITE(6,2003)
      IF (LPL(L).EQ.0) GO TO 5
      IF (LPL(2).EQ.2) GO TO 5
      IF (L.EQ.1)  GO TO  58
      DO 57 I = 1,N
   57 AA(1,I) = AA(2,I)
      DO 50 I = 1,11
   50 ID(1,I) = ID(2,I)
   58 CONTINUE
      GO TO 5
    5 CONTINUE
 2000 FORMAT(11A6,5HLAYER I5)
 2001 FORMAT(8F9.6,I7)
 2002 FORMAT(41H1  TIME HISTORY OF STRAIN IN PERCENT         )
 2003 FORMAT(41H1  TIME HISTORY OF STRESS IN KSF            )
      RETURN
      END
