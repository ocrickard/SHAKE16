C $NOFLOATCALLS
C $NODEBUG
C************************************
      subroutine AMP( N1,IN,INT,LL,LT,KPL,IDAMP,NA,DF)
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
      complex G, V, PLUS, MINUS
      complex E, F, EE, FF, A, EX, AIN, IPI2,AA
      character*60 ABSIS
      character*6 ID,IDNT,IDAMP
C
      dimension IDAMP (27,11),T(200)
      common /JOE4/ ST(27,200)
      common /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      common /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
      common /CCG/ ID(27,11)
C
      ABSIS = ' CYCLES/SEC. '
C
      IPI2 = CMPLX(0., 6.283185307)
      FREQ = 0.
      ST(NA,1) = 1.
      do 19 I = 2,200
      E = 1.
      FF = 1.
      FREQ = FREQ + DF
      A = FREQ*IPI2
      do 191 K = 1,N1
      if (K.NE.IN) go to 192
      AIN = E + FF
      if (INT.EQ.0) AIN = 2.*E
  192 if (K.NE.LL) go to 11
      AA    = E + FF
      if (LT.EQ.0) AA = 2.*E
   11 EX = CEXP(H(K)*A/V(K))
      EE = E*EX
      F = FF/EX
      E = EE*PLUS(K) + MINUS(K)*F
      FF = PLUS(K)*F + MINUS(K)*EE
  191 continue
      if (IN.NE.N1+1) go to 193
      AIN = E + FF
      if (INT.EQ.0) AIN = 2.*E
  193 if (LL   .NE.N1+1) go to 21
      AA    = E + FF
      if (LT.EQ.0)  AA = 2.*E
   21 ST(NA,I) = CABS(AA/AIN)
   19 continue
      do 23 I = 1,200
   23 T(I) = DF*FLOAT(I-1)
      AMAX = 0.
      write(6,2)
      do 22 I = 1,200
      if (KPL .GE. 2) write(6,1) T(I), ST(NA,I)
      if (ST(NA,I) .LT.  AMAX) go to 22
      TMAX = T(I)
      AMAX = ST(NA,I)
   22 continue
      if (NA.LT.9) NA=NA+1
      PERIOD = 1./TMAX
      if (TMAX.LT. .0001)  write(6,1001) AMAX, TMAX
      if (TMAX.GT. .0001)  write(6,1001) AMAX, TMAX,PERIOD
      if (KPL.EQ.0) return
      write(6,1000)
      N = NA-1
      NA = 1
      return
    1 format(1X,F10.4, 3X, F10.4)
    2 format(/2X,'FREQUENCY    AMPLITUDE')
 1000 format(33H1  PLOT OF AMPLIFICATION SPECTRA  /)
 1001 format(25H MAXIMUM AMPLIFICATION =   F6.2/
     1 25H FOR FREQUENCY         = F6.2, 7H C/SEC. /
     1 25H     PERIOD            =  F6.2, 5H SEC. )
      end
C******************************************
      subroutine UTPR(KK,DPTH,LS,K2,LH,LT,X,AX,S,INV)
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
      character*6 TITLE,IDNT
      complex SAVE
      complex X, AX
C
      dimension XR(8)
      dimension X(300),AX(3,270),S(70),INV(70)
      common /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      common /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
C
      FREQ = 0.
      SFX = 0.
      SXX = 0.
C   TRANSFORM VALUES IN X OR IN AX INTO THE TIME DOMAIN
      do 24  I = 1,MFOLD
      if (LS.EQ.0) go to 241
      SAVE = X(I)
      X(I) = AX(LS,I)
      AX(LS,I) = SAVE
  241 XA = CABS(X(I))
      SXX= SXX + XA*XA
      SFX = SFX + XA*FREQ*XA
      FREQ = FREQ + DF
   24 continue
      SFX = SFX/SXX
C
      call RFSN(X,MX,INV,S,IFERR,-2)
C
      call XMX(X,MA,XMAX,NMAX)
      TMAX = DT*FLOAT(NMAX-1)
      XEND = 0.
      N = MA/20
      NN = 9*N
      N = 8*N
      do 25 I = N,NN
      XABS = REAL(X(I))
      XABS = ABS(XABS)
      if (XABS.GT.XEND) XEND = XABS
      XABS = AIMAG(X(I))
      XABS = ABS(XABS)
      if (XABS.GT.XEND) XEND = XABS
   25 continue
      XEND = XEND/XMAX
C
C   SAVE OUTPUT
      N = 1
      NN = 4
      NCARDS=MA/8
      NC = NCARDS
      if (K2.EQ.0) NC = 0
C     if (KK.EQ.5) go to 252
      if (KK.EQ.6) go to 252
      if (LT.EQ.0)  write(6,2000) LH, (IDNT(I),I=1,6)
      if (LT.EQ.1)  write(6,2002) LH, (IDNT(I),I=1,6)
      write(6, 2005) SFX
      write(6,2003) XMAX, TMAX
  252 if (KK.EQ.6.AND.LT.EQ.0)  write(6,2001) DPTH,XMAX,TMAX,SFX,XEND,NC
      if (KK.EQ.6.AND.LT.EQ.1) write(6,2010) DPTH,XMAX,TMAX,SFX, XEND,NC
      if (K2.EQ.0) go to 262
      write(7,2006) XMAX, (TITLE(I),I=1,5)
      if (LT.EQ.1)  write(7,2002) LH, (IDNT(I),I=1,6)
      if (LT.EQ.0)  write(7,2000) LH, (IDNT(I),I=1,6)
      do 26 I = 1,NCARDS
      K = 0
      do 261 J = N,NN
      K = K+ 1
      XR(K) = REAL(X(J))
      K = K + 1
      XR(K) = AIMAG(X(J))
  261 continue
      write(7,2009) (XR(J),J=1,8),I
      if (K2 .EQ. 2) write(6,2019) (XR(J), J = 1,8), I
      NN = 4 + NN
      N = N + 4
   26 continue
  262 call RFFT(X,MX,INV,S,IFERR,2)
      if (LS.EQ.0)   return
      do 27 I = 1,MFOLD
      SAVE = AX(LS,I)
      AX(LS,I) = X(I)
   27 X(I) = SAVE
      return
C
 2000 format(43H  ACCELERATION VALUES AT OUTCROPPING LAYER I3,3H - 6A6)
 2001 format(5X,6HOUTCR. F15.1,F15.5,2F15.2,F20.3,I20)
 2010 format(5X,6HWITHIN F15.1,F15.5,2F15.2,F20.3,I20)
 2002 format(42H  ACCELERATION VALUES AT THE TOP OF LAYER I3,3H - 6A6)
 2003 format(/15H  MAX. ACC. =  F9.6,11H AT TIME = F6.3, 5H SEC. /)
 2005 format(/26H   MEAN SQUARE FREQUENCY = F10.2/)
 2006 format(21X,6HXMAX= F7.4,5A6)
 2008 format(2X, I5, 5X, 8F15.6)
 2009 format(8F9.6,I7)
 2019 format(8F14.6,I10)
      end
C**************************************
      subroutine REDUCE(IFR,X,AX,LL)
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
      character*6 TITLE
      complex X, AX
C
      dimension X( 68), AX(3, 64), LL(3)
      common /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
      common/FRCUT/ NCUT,NZERO
C
      F1 = .5/DT
      FR = FLOAT(IFR)
      DT = DT*FR
      MA = MA/IFR
      MMA = MMA/IFR
      MA2 = MA + 2
      MFOLD = MA2/2
      N = MFOLD + 1
      do 12 I = MFOLD,N
      X(I) = 0.
C
      do 12 L = 1,3
      if (LL(L).LE.0) go to 12
      AX(L,I) = 0.
   12 continue
      MFOLD = MFOLD + 1
      F2 = .5/DT
      write(6,1000) F1,F2,DT, MA
      FMA = FLOAT(MA)
      MX = (ALOG10(FMA)/ALOG10(2.))-1.
      if (MA.GT.2**(MX+1)) MX=MX+1
      IF(NCUT.LE.MFOLD) go to 15
      NCUT=MFOLD
   15 continue
 1000 format(  20H  FREQUENCIES FROM  F6.2, 3H TO F6.2,14H C/SEC ARE REM
     15HOVED /
     216H NEW TIMESTEP = F5.4/19H NUMBER OF VALUES = I5)
      return
      end
C*****************************
      subroutine INCR(IFR,X,AX)
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
      complex X, AX
      character*6 TITLE
C
      dimension X( 68), AX(3, 64)
      common /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
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
      do 10 I = N, MFOLD
      X(I) = 0.
      do 10 L = 1,3
   10 AX(L,I) = 0.
      F2 = .5/DT
      write(6,1000) F1,F2,DT, MA
      FMA = FLOAT(MA)
      MX = (ALOG10(FMA)/ALOG10(2.))-1.
      if (MA.GT.2**(MX+1)) MX=MX+1
 1000 format(27H    FREQUENCIES ADDED FROM  F6.2,3H TO F6.2/
     216H NEW TIME STEP = F5.4/19H NUMBER OF VALUES = I5/)
      return
      end
C***********************************************************
      subroutine MOTION(N1,IN,INT,LL,LT, X,AX)
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
      integer LL(3), LT(3)
      character*6 TITLE,IDNT
      complex AA(3)
      complex X, AX
      complex G, V, PLUS, MINUS
      complex E, F, EE, FF, A, EX, AIN, IPI2
C
      dimension X(300),AX(3,270),S(70)
      common /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
      common /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      common /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
      common/FRCUT/ NCUT,NZERO
C
      IPI2 = CMPLX(0., 6.283185307)
      do 20 L = 1,3
      if (LL(L) .GT. 0) AX(L,1) = X(1)
      IF(NCUT.EQ.MFOLD) go to 20
      do 30 I=NZERO,MFOLD
      AX(L,I)=CMPLX(0.,0.)
   30 continue
   20 continue
      FREQ = 0.
      do 19 I=2,NCUT
      E = 1.
      FF = 1.
      FREQ = FREQ + DF
      A = FREQ*IPI2
      do 191 K = 1,N1
      if (K.NE.IN) go to 192
      AIN = E + FF
      if (INT.EQ.0) AIN = 2.*E
C   FIND SUBLAYER WHERE MOTION IS WANTED
  192 do 11 L = 1,3
      if (K.NE.LL(L)) go to 11
C   AMPLIFICATION FACTOR FOR SUBLAYER WITHIN PROFILE
      AA(L) = E + FF
C   AMPLIFICATION FACTOR FOR OUTCROPPING SUBLAYER
      if (LT(L).EQ.0)  AA(L) = 2.*E
   11 continue
      EX = CEXP(H(K)*A/V(K))
      EE = E*EX
      F = FF/EX
      E  = EE*PLUS(K) + MINUS(K)*F
      FF = PLUS(K)*F + MINUS(K)*EE
  191 continue
      if (IN.NE.N1+1) go to 193
      AIN = E + FF
      if (INT.EQ.0) AIN = 2.*E
  193 do 21 L = 1,3
      if (LL(L).NE.N1+1) go to 21
      AA(L) = E + FF
      if (LT(L).EQ.0) AA(L) = 2.*E
   21 continue
      do 23 L = 1,3
      if (LL(L) .GT. 0) AX(L,I) = X(I)*AA(L)/AIN
   23 continue
   19 continue
      return
      end
C*************************************************
      subroutine CXSOIL(N1)
C***********************************************************************
C
C   THIS ROUTINE CALCULATES THE complex SOIL PROPERTIES AND TRANSFER
C   FUNCTIONS FOR THE LAYERS
C
C        N1      = NUMBER OF SOIL LAYERS
C        BL      = RATIO OF CRITICAL DAMPING
C        GL      = SHEAR MODULUS
C        R       = DENSITY
C        G       = complex SHEAR MODULUS
C        V       = complex SHEAR WAVE VELOCITY
C        PLUS    = complex TRANSFER FUNCTION
C        MINUS   = complex TRANSFER FUNCTION
C
C   CODED BY PER B SCHNABEL OCT 1971
C
C***********************************************************************
C
      complex G, V, PLUS, MINUS, MU
      character*6 IDNT
      common /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      common /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
C
      N = N1 + 1
      do 1 I = 1,N
      GIMAG=2.*BL(I)*GL(I)*SQRT(1.-BL(I)*BL(I))
      GREAL=GL(I)*(1.-2.*BL(I)*BL(I))
      G(I)=CMPLX(GREAL,GIMAG)
      V(I) = CSQRT(G(I)/R(I))
    1 continue
      do 2 I = 1,N1
      J = I + 1
      MU = CSQRT(R(I)/R(J)*G(I)/G(J))
      PLUS(I) = (1. + MU)/2.
      MINUS(I)= (1. - MU)/2.
    2 continue
      return
      end
C*********************************************
      subroutine STRAIN( LL, LGS, LPCH, LPL,LNV,X,AX,AA,N1,S,INV)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C   THIS subroutine COMPUTES STRAIN AND/OR STRESS TIME-HISTORY AT THE
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
      integer TP
      character*6 TITLE,IDNT,ID
      character*60 ABSIS
      complex X, AX
      complex G, V, PLUS, MINUS
      complex E,F,EE, A,AH,IPI2, AE,AF,EX,AI
C
      dimension AE(2), AF(2)
      dimension X(1),  AX(3,1), AA(2,1), S(1),  INV(1)
      dimension LL(2), LGS(2),  LPCH(2), LPL(2), LNV(2)
C
      common /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      common /SOILB/ FAC(51), WL(51), TP(51), DEPTH(51), WEIGHT(51)
      common /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
      common /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
      common /CCG/ ID(27,11)
      common /FRCUT/ NCUT,NZERO
      common /TIME/ T(9)
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
      do 1 I=2,NCUT
      E = AX(1,I)/2.
      F = E
      FREQ = FREQ + DF
      AH = AI/FREQ
      A = FREQ*IPI2
      do 11  K = 1,N1
      do 12 L = 1,2
      if (K.NE.LL(L)) go to 12
      AE(L) = E/V(K)
      AF(L) = F/V(K)
   12 continue
      EX = CEXP(H(K)*A/V(K))
      E = E*EX
      F = F/EX
      EE = E*PLUS(K) + MINUS(K)*F
      F  = F*PLUS(K) + MINUS(K)*E
      E = EE
   11 continue
      do 13 L = 1,2
      if (LL(L).NE.N1+1) go to 13
      AE(L) = E/V(N1+1)
      AF(L) = F/V(N1+1)
   13 continue
      do 14 L = 1,2
      if (LL(L).GT.0) AX(L+1,I) = (AE(L) -AF(L))*AH
   14 continue
    1 continue
      do 2 I = 1,MFOLD
    2 AX(1,I) = X(I)
      do 3 L = 1,2
      if (LL(L).EQ.0) go to 3
      X(1) = 0.
      do 31 I=2,NCUT
   31 X(I) = AX(L+1,I)
      IF(NCUT.EQ.MFOLD) go to 33
      do 34 II=NZERO,MFOLD
      X(II)=CMPLX(0.,0.)
   34 continue
   33 continue
      call RFSN(X,MX,INV,S,IFERR,-2)
      do 32 I =1,MFOLD
      AA(L,2*I-1) =REAL(X(I))*100.
   32 AA(L,2*I) = AIMAG(X(I))*100.
    3 continue
C
      do 4 I = 1,MFOLD
    4 X(I) = AX(1,I)
C   COMPUTE STRESS IF WANTED AND SAVE COMPUTED RESPONSES
      do 5 L = 1,2
      if (LL(L) .EQ. 0) go to 5
      NVAL = LNV(L)
      if (NVAL.LE.0) NVAL = MMA
      if (NVAL.GT.MA)  NVAL = MA
      if (NVAL.GT.2049) NVAL = 2049
      do 51 I = 1,5
   51 ID(L,I) = TITLE(I)
      N = LL(L)
      ID(L,6) = 'STRAIN'
      if (LGS(L) .EQ.0) go to 53
      ID(L,6) = 'STRESS'
      do 52 I = 1,NVAL
   52 AA(L,I) = GL(N)*AA(L,I)/100.
   53 if (LPCH(L).EQ.0) go to 54
      write(7,2000) (ID(L,I), I=1,11),N
      N = 1
      NCARDS = NVAL/8
      do 55 K = 1,NCARDS
      NN = N + 7
      write(7,2001) (AA(L,I), I = N,NN), K
   55 N = N + 8
   54 if (LPL(L).EQ.0) go to 5
      N = 0
      NSKIP = 1
      do 56 I = 1,NVAL,NSKIP
      N = N + 1
      if (NSKIP.GT.1) AA(L,N) = AA(L,I)
   56 T(N) = DT*FLOAT(I-1)
      IF  (LGS(L).EQ.0)  write(6,2002)
      if (LGS(L).EQ.1)  write(6,2003)
      if (LPL(L).EQ.0) go to 5
      if (LPL(2).EQ.2) go to 5
      if (L.EQ.1)  go to  58
      do 57 I = 1,N
   57 AA(1,I) = AA(2,I)
      do 50 I = 1,11
   50 ID(1,I) = ID(2,I)
   58 continue
      go to 5
    5 continue
 2000 format(11A6,5HLAYER I5)
 2001 format(8F9.6,I7)
 2002 format(41H1  TIME HISTORY OF STRAIN IN PERCENT         )
 2003 format(41H1  TIME HISTORY OF STRESS IN KSF            )
      return
      end
