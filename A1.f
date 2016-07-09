C**********************************************************************
      subroutine  CURVEG(NC, NV, K1, A, B, NN, TSTEP, NT, T,V,X,Y,NSTEP)
C***********************************************************************
C   THE PROGRAM GENERATES NEW POINTS ON A CURVE BY LINEAR INTERPOLATION
C   USING AN ARITHMETIC OR A HALFLOGARITHMIC SCALE
C
C        NV(I)     =    NUMBER OF VALUES ON CURVE I
C        NC        =    NUMBER OF CURVES
C        K1        =    SWITCH    K1 = 1  ARITHMETIC SCALE
C                                 K1 = 2  HALFLOGARITHMIC SCALE
C        A,B       =    PARAMETERS FOR CALCULATING NEW VALUES
C                       Y = A*X + B
C        X,Y       =    KNOWN POINTS ON CURVE
C        T         =    VALUES ON ABSISSA WHERE NEW POINTS ARE GENERATED
C        V         =    NEW ORDINATE VALUES
C
C        ARITHMETIC SCALE  K1 = 1
C        NN        =    NUMBER OF INTERVALS
C        TSTEP     =    LARGEST VALUE IN EACH INTERVAL
C        NT        =    NUMBER OF STEPS IN EACH INTERVAL
C
C        HALFLOGARITHMIC SCALE
C        NN        =    NUMBER OF VALUES IN EACH LOG10
C
C   CODED BY PER B SCHNABEL SEPT 1970
C
C***********************************************************************
C
      dimension X(27,20),Y(27,20),A(27,20),B(27,20),NV(27),TSTEP(27)
      dimension NT(27), T(200), V(27,200)
C
      XMIN = 100000000.
      XMAX = 0.
      do L= 1,NC
        M = NV(L)
        if (XMAX .lt. X(L,M))   XMAX = X(L,M)
        if (XMIN .gt. X(L,1)) XMIN = X(L,1)
        M = M - 1
        do I = 1,M
          X1 = X(L,I)
          X2 = X(L,I+1)
          if (K1 .eq. 2)  X1 = ALOG10(X1)
          if (K1 .eq. 2)  X2 = ALOG10(X2)
          X(L,I) = X(L,I+1)
          A(L,I) = (Y(L,I+1) - Y(L,I))/(X2 - X1)
          B(L,I) = -A(L,I)*X1 + Y(L,I)
        end do
      end do
C
      call STEPG(K1, NN, TSTEP, NT, XMIN, XMAX, T, NSTEP)
C
      do L = 1,NC
        M = NV(L) - 1
        do I = 1,NSTEP
          do J = 1,M
            if (T(I) .lt. X(L,J))  go to  31
          end do
          J = M
   31     TT = T(I)
          if (K1 .eq. 2) TT = ALOG10(TT)
          V(L,I) = A(L,J)*TT  + B(L,J)
        end do
      end do
      return
      end
C*********************************************************************
      subroutine STEPG(KK, NN, TSTEP, NT, T1, TN, T, NSTEP)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  *
C
C   THE ROUTINE GENERATES STEPS IN LINEAR OR LOGARITHMIC INCREMENT
C
C        KK      = SWITCH    KK = 1    STEP INCREASE OF VALUES
C                            KK = 2    LOGARITHMIC INCREASE OF VALUES
C        NN      = NUMBER OF STEPS   OR NUMBER OF VALUES IN EACH 10
C        TSTEP   = LARGEST VALUE IN EACH STEP
C        NT      = NUMBER OF VALUES IN EACH STEP
C        T1      = FIRST VALUE IN LOG-STEP
C        TN      = LAST  VALUE IN LOG-STEP
C        T       = VALUES GENERATED
C        NSTEP   = NUMBER OF VALUES
C
C   CODED PER B SCHNABEL SEPT. 1970
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  *
C
      dimension T(200), TSTEP(27), NT(27)
C
      if (KK .eq. 1) then
        K = 1
        T(K) = 0.
        SAVE = 0.
        do N = 1,NN
          M = NT(N)
          STEP = (TSTEP(N) - SAVE)/FLOAT(M)
          SAVE = TSTEP(N)
          do I = 1,M
            K = K + 1
            T(K) = T(K-1) + STEP
          end do
        end do
        NSTEP = K
        return
      else if (KK .eq. 2) then
        NST = ALOG10(T1)
        if (T1.lt. 1.) NST = NST - 1
        STEP = 1./NN
        K = 1
        TA = 10.**FLOAT(NST)
        T(1) = TA
        do J = 2,NN
          K = K + 1
          T(K) = TA*10.**(STEP*FLOAT(J))
          if (T(K) .gt. T1) go to 221
        end do
  221   TA = T(K-1)
        K = 0
  211   do J = 1,NN
          K = K + 1
          T(K) = TA*10.**(STEP*FLOAT(J))
          if (T(K) .gt. TN)  go to 212
        end do
        TA = TA*10.
        go to 211
  212   NSTEP = K
        return
      end if
      end
C***********************************************************
      subroutine RESP(LN,LS,NN,X,AX,A,S,INV)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   THIS PROGRAM readS DATA FOR RESPONSE SPECTRUM ANALYSIS
C   NECESSARY subroutineS    DRCTSP,  CMPMAX
C
C        NN      = RESPONSE SPECTRUM NUMBER
C        ND      = NUMBER OF DAMPING VALUES
C        X       = FOURIER TRANSFORM OF OBJECT MOTION
C        AX      = FOURIER TRANSFORM OF COMPUTED MOTIONS
C        T       = PERIODS FOR WHICH RESPONSE IS to BE COMPUTED
C
C   CODED PER B SCHNABEL DEC. 1970
C   New Sets of Periods -- included in February 1991
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      character*6 TITLE,ID,IBLANK,IDNT
      character*60 ABSIS
      character*32 FPERIOD
      character*80 headerd
      complex X, AX
C
      dimension X(64), AX(3,64),A(2,64), S(10),INV(10)
      dimension ID(27,11)
C
      common /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      common /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
      common /RVAL/ NND(27), ZLD(6),T(200), SA(5,200),SV(5,200)
C
      IBLANK = '      '
      ABSIS = ' PERIOD IN SEC.'
C
      read(5,4) ND, KPER, GGT
    4 format(2I5,F10.2)
      read(5,5) (ZLD(I), I = 1,ND)
    5 format(6F10.3)
      write(6,9001) LN, (ZLD(I), I = 1,ND)
C -------------------------------------------------------------------
C   if KPER = 0; Periods from 0.03 to 10 sec are included in data block
C                in this subroutine
C    otherwise, periods are specified by user (maximum is 200 periods)
C -------------------------------------------------------------------
      if (KPER.eq. 0) go to 99
        read(5,'(A32)') FPERIOD
        write(6,60) FPERIOD
   60   format(' File from which periods were read: ' A32)
        open(8,FILE=FPERIOD,STATUS='OLD')
        read (8,4) NLINES, NNM
        do I = 1, NLINES
          read(8,*) headerd
          write(6,*) headerd
        end do
        read(8,*) (T(I), I=1, NNM)
        close(8)
        go to 101
C ----------------------------------------------------------------------
C  default periods for calculating response spectra
C  ---------------------------------------------------------------------
   99   NNM=152
      T(1) = .01
      data (t(i), i=2,152)/
     1     0.03,    0.04,    0.05,    0.06,    0.07,    0.08,    0.09,
     2     0.10,    0.11,    0.12,    0.13,    0.14,    0.15,    0.16,
     3     0.17,    0.18,    0.19,    0.20,    0.21,    0.22,    0.23,
     4     0.24,    0.25,    0.26,    0.27,    0.28,    0.29,    0.30,
     5     0.31,    0.32,    0.33,    0.34,    0.35,    0.36,    0.37,
     6     0.38,    0.39,    0.40,    0.41,    0.42,    0.43,    0.44,
     7     0.45,    0.46,    0.47,    0.48,    0.49,    0.50,    0.51,
     8     0.52,    0.53,    0.54,    0.55,    0.56,    0.57,    0.58,
     9     0.60,    0.62,    0.64,    0.66,    0.68,    0.70,    0.72,
     T     0.74,    0.76,    0.78,    0.80,    0.82,    0.84,    0.86,
     1     0.88,    0.90,    0.92,    0.94,    0.96,    0.98,    1.00,
     2     1.05,    1.10,    1.15,    1.20,    1.25,    1.30,    1.35,
     3     1.40,    1.45,    1.50,    1.55,    1.60,    1.65,    1.70,
     4     1.75,    1.80,    1.85,    1.90,    1.95,    2.00,    2.05,
     5     2.10,    2.15,    2.20,    2.25,    2.30,    2.35,    2.40,
     6     2.50,    2.60,    2.70,    2.80,    2.90,    3.00,    3.10,
     7     3.20,    3.30,    3.40,    3.50,    3.60,    3.70,    3.80,
     8     3.90,    4.00,    4.10,    4.20,    4.30,    4.40,    4.50,
     9     4.60,    4.70,    4.80,    4.90,    5.00,    5.10,    5.20,
     T     5.40,    5.60,    5.80,    6.00,    6.20,    6.40,    6.60,
     1     6.80,    7.00,    7.20,    7.40,    7.60,    7.80,    8.00,
     2     8.50,    9.00,    9.50,   10.00/
c ---------------------------------------------------------------------
C   SAVE VALUES OF X IN AA
  101 do I = 1,MFOLD
        A(1,I) = real(X(I))
        A(2,I) = AIMAG(X(I))
        if (LS .ne. 0) then
          X(I) = AX(LS,I)
        end if
      end do
C
C   TRANSFORM VALUES IN X OR AX INto THE TIME DOMAIN
      call RFSN(X,MX,INV,S,IFERR,-2)
      do L = 1,ND
        if (NN .ge. 5)  NN= 0
        NN = NN + 1
        do I = 1,5
          ID(NN,I) = TITLE(I)
        end do
        do I = 6,11
          ID(NN,I) = IDNT(I-5)
          if (LS .eq. 0) ID(NN,I) = IBLANK
        end do
C
C       COMPUTE RESPONSE FOR ACCELERATION VALUES IN AA(1, )FOR THE PERIODS
C       GIVEN IN T( )
        call  DRCTSP(NN,MMA, DT, GGT, ID, ZLD(L),NNM,X)
      end do
C
C   GIVE X BACK ORIGINAL VALUES
      do I = 1,MFOLD
        X(I) = CMPLX(A(1,I),A(2,I))
      end do
C     ==============================================================
  134 NN = 0
      return
 1000 format(10I5)
 9000 format( 8F10.3)
 9001 format( 50H  RESPONSE SPECTRUM ANALYSIS FOR LAYER NUMBER     I4
     1/26H CALCULATED FOR DAMPING    8F10.3)
      end
C**********************************************************************
      subroutine DRCTSP(NN, KG, DT, GGT, ID, D, M, A)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   THIS ROUTINE COMPUTES RESPONSE SPECTRA BY THE STEP BY STEP METHOD
C
C       NN      = RESPONSE SPECTRUM CURVE NUMBER  USED (Canceled)
C       KG      = NUMBER OF ACCELERATION VALUES
C       DT      = TIME STEP BETWEEN EACH ACCELERATION VALUE
C       M       = NUMBER OF PERIODS FOR WHICH RESPONSE IS to BE COMPUTED
C       T       = ARRAY WITH THE PERIODS
C       A       = ACCELERATION VALUES
C       D       = CRITICAL DAMPING RATIO
C       ID      = IDENTIFICATION
C       GGT     = Acceleration of gravity - cm/sec/sec, or in/sec/sec
C                  or ft/sec/sec
C
C   CODED BY I. M. IDRISS 1967
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      character*6 ID
      dimension A(1)
      common /RVAL/ NND(27), ZLD(6),T(200), SA(5,200),SV(5,200)
      dimension PRV(200), PAA(200), RD(200)
      dimension ID(27,11)
C .....................................................................
      zmax =0
      do K = 1, KG
        if (zmax .gt. ABS(A(K))) then
          A(K) = GGT*A(K)
        else
          zmax = ABS(A(K))
        end if
      end do
      PIW = 6.283185307
      SV(NN,1) = zmax*GGT*T(1)/PIW
      SA(NN,1) = zmax
      KUG = KG-1
      RD(1) = zmax*GGT*T(1)*T(1)/(PIW*PIW)
      PRV(1) = zmax*GGT*T(1)/PIW
      PAA(1) = zmax
      write(6,112) D
      N = 1
      YY = SQRT(1.-D*D)
      do LOOP = 2, M
        W = 6.283185307/T(N)
        WD = YY*W
        W2 = W*W
        W3 = W2*W
        call CMPMAX(KUG,T(N),W,W2,W3,WD,D,DT,ZD,ZV,ZA,A)
        SV(NN,N) = ZV
        SA(NN,N) = ZA/GGT
        RD(N) = ZD
        PRV(N) = W*ZD
        PAA(N) = W2*ZD/GGT
        N = N + 1
      end do
      write(6,312) GGT, (ID(NN,I), I = 1,10),D
      SUMSV = 0.
      SUMSA = 0.
      SUMT  = 0.
      SVMAX = 0.
      SAMAX = 0.
      TT1 = .1
      TT2 = 0.
      do N = 1, M
        FREKV = 1./T(N)
        if (T(N) .lt. .0999 .or. TT2.gt.2.4999) then
          continue
        end if
        TT2 = (T(N+1) + T(N))/2.
        if (TT2.gt.2.5)  TT2 = 2.5
        TT = TT2 - TT1
        SUMSA = SA(NN,N)*TT + SUMSA
        SUMSV = SV(NN,N)*TT + SUMSV
        SUMT = SUMT + TT
        TT1 = TT2
        if (SVMAX.lt.SV(NN,N))  SVMAX = SV(NN,N)
        if (SAMAX.lt.SA(NN,N))  SAMAX = SA(NN,N)
        write(6,322) N,T(N),RD(N),SV(NN,N),PRV(N),SA(NN,N),PAA(N),FREKV
      end do
      write(6,2002) SUMSA,SUMSV,SAMAX,SVMAX
      do K = 1,KG
        A(K) = A(K)/GGT
      end do
      return
C
  112 format(/5X,41HTIMES AT WHICH MAX. SPECTRAL VALUES OCCUR /
     1 10X,33HTD = TIME FOR MAX. RELATIVE DISP. /
     2 10X,33HTV = TIME FOR MAX. RELATIVE VEL.  /
     3 10X,33HTA = TIME FOR MAX. ABSOLUTE ACC.   /
     4 5X, 15HDAMPING RATIO =  F5.2)
  312 format(5X,' SPECTRAL VALUES --'/
     15X,' [Acceleration of gravity used =' F8.2 ']'/
     210A6,2X,15HDAMPING RATIO =
     3 F5.2/5X,3HNO.,4X,6HPERIOD,5X,10HREL. DISP.,6X,9HREL. VEL.,3X,
     4 12HPSU.REL.VEL.,6X,9HABS. ACC.,3X,12HPSU.ABS.ACC. 5X,5HFREQ.)
  322 format(I8,F10.2,5F15.5,F10.2)
  402 format(8F9.5)
  412 format(I5,25H ACC. RESPONSE VALUES FOR   , 8A6)
  413 format(I5,25H VEL. RESPONSE VALUES FOR   , 8A6)
 2002 format(10X,40HVALUES IN PERIOD RANGE .1 to 2.5 SEC.        /
     115X,35HAREA OF ACC. RESPONSE SPECTRUM   =   F10.3/
     215X,35HAREA OF VEL. RESPONSE SPECTRUM   =   F10.3/
     315X,35HMAX. ACCELERATION RESPONSE VALUE =   F10.3/
     415X,35HMAX. VELOCITY     RESPONSE VALUE =   F10.3)
      end
C  ********************************************************************
      subroutine CMPMAX (KUG,PR,W,W2,W3,WD,D,DT,ZD,ZV,ZA, UG)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C   THIS ROUTINE COMPUTES RESPONSE VALUES FOR ONE SINGLE DEGREE OF
C   FREEDOM SYSTEM USING STEP BY STEP METHOD
C
C   EXPLANATIOS to PARAMETERS GIVEN IN DRCTSP
C
C   CODED BY I. M. IDRISS  1967
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     TODO: These arrays need to be checked for size
      dimension XD(2), XV(2), T(3)
      dimension UG(1)
C
      ZA = 0.
      ZD = 0.
      ZV = 0.
      XD(1) = 0.
      XV(1) = 0.
      F1 = 2.*D/(W3*DT)
      F2 = 1./W2
      F3 = D*W
      F4 = 1./WD
      F5 = F3*F4
      F6 = 2.*F3
       E = EXP(-F3*DT)
       S = SIN(WD*DT)
       C= COS(WD*DT)
      G1 = E*S
      G2 = E*C
      H1 = WD*G2 - F3*G1
      H2 = WD*G1 + F3*G2
      do K = 1, KUG
        Y = K-1
        DUG = UG(K+1) - UG(K)
        Z1 = F2*DUG
        Z2 = F2*UG(K)
        Z3 = F1*DUG
        Z4 = Z1/DT
        B = XD(1) + Z2 -Z3
        A = F4*XV(1) + F5*B + F4*Z4
        XD(2) = A*G1 + B*G2 + Z3 - Z2 - Z1
        XV(2) = A*H1 - B*H2 - Z4
        XD(1) = XD(2)
        XV(1) = XV(2)
        AA = -F6*XV(1) - W2*XD(1)
        F = ABS(XD(1))
        G = ABS(XV(1))
        H = ABS(AA)
        if (F .gt. ZD) then
          T(1) = Y
          ZD = F
        end if
        if (G .gt. ZV) then
          T(2) = Y
          ZV = G
        end if
        if (H .gt. ZA) then
          T(3) = Y
          ZA = H
        end if
      end do
      do L = 1, 3
        T(L) = DT*T(L)
      end do
      write(6,112) PR, (T(L),L=1,3)
  112 format(5X,5HPER = F5.2,5X,19HTIMES FOR MAXIMA -- ,3X,
     14HTD = F8.4,3X,4HTV = F8.4,3X,4HTA = F8.4)
      return
      end
C    ****************************************
      subroutine FFT (A,M,INV,S,IFSET,IFERR)
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     TODO: Check dimensions on these arrays
      dimension A(1),INV(1),S(1),N(3),M(3),NP(3),W(2),W2(2),W3(2)
      equivalence (N1,N(1)), (N2,N(2)), (N3,N(3))

      M1=M(1)
      M2=M(2)
      M3=M(3)
      MTT=M1-2
      MT=MAX0(2,MTT)
      NT=2**MT
10    if (IABS(IFSET)-1) 610,610,20
610   MT=MAX0(M(1),M(2),M(3))-2
      MT=MAX0(2,MT)
      if (MT-20) 630,630,620
620   IFERR=1
      go to 600
630   IFERR=0
      NT=2**MT
      NTV2=NT/2
      THETA=.7853981634
      JSTEP=NT
      JDIF=NTV2
      S(JDIF)=SIN(THETA)
      do L=2,MT
        THETA=THETA/2.
        JSTEP2=JSTEP
        JSTEP=JDIF
        JDIF=JSTEP/2
        S(JDIF)=SIN(THETA)
        JC1=NT-JDIF
        S(JC1)=COS(THETA)
        JLAST=NT-JSTEP2
        if (JLAST-JSTEP .ge. 0) then
          do J=JSTEP,JLAST,JSTEP
            JC=NT-J
            JD=J+JDIF
            S(JD)=S(J)*S(JC1)+S(JDIF)*S(JC)
          end do
        end if
      end do
C
C     SET UP INV(J) TABLE
      MTLEXP=NTV2
      LM1EXP=1
      INV(1)=0
      do L=1,MT
        INV(LM1EXP+1)=MTLEXP
        do J=2,LM1EXP
          JJ=J+LM1EXP
          INV(JJ)=INV(J)+MTLEXP
        end do
        MTLEXP=MTLEXP/2
        LM1EXP=LM1EXP*2
      end do
      if (IFSET) 20,600,20
20    MTT=MAX0(M(1),M(2),M(3))-2
      ROOT2=SQRT(2.)
      if (MTT-MT .gt. 0) then
        IFERR=1
        write(6,1000)
        stop
1000    format(31H --- ERROR IN FOURIER TRANSFORM  )
      end if
      IFERR=0
C     M1=M(1)
C     M2=M(2)
C     M3=M(3)
      N1=2**M1
      N2=2**M2
      N3=2**M3
      if (IFSET .le. 0) then
        NX=N1*N2*N3
        FN=NX
        do I=1,NX
          A(2*I-1)=A(2*I-1)/FN
          A(2*I)=-A(2*I)/FN
        end do
      end if
      NP(1)=N1*2
      NP(2)=NP(1)*N2
      NP(3)=NP(2)*N3
      do ID=1,3
        IL=NP(3)-NP(ID)
        IL1=IL+1
        MI=M(ID)
        if (MI) 330,330,80
80      IDIF=NP(ID)
        KBIT=NP(ID)
        MEV=2*(MI/2)
        if (MI-MEV .gt. 0) then
          KBIT=KBIT/2
          KL=KBIT-2
          do I=1,IL1,IDIF
            KLAST=KL+I
            do K=I,KLAST,2
              KD=K+KBIT
              T=A(KD)
              A(KD)=A(K)-T
              A(K)=A(K)+T
              T=A(KD+1)
              A(KD+1)=A(K+1)-T
              A(K+1)=A(K+1)+T
            end do
          end do
          if (MI-1 .gt. 0) then
            LFIRST=3
            JLAST=1
          else
            go to 330
          end if
        else
120       LFIRST=2
          JLAST=0
        end if
130     do 320 L=LFIRST,MI,2
        JJDIF=KBIT
        KBIT=KBIT/4
        KL=KBIT-2
        do I=1,IL1,IDIF
          KLAST=I+KL
          do K=I,KLAST,2
            K1=K+KBIT
            K2=K1+KBIT
            K3=K2+KBIT
            T=A(K2)
            A(K2)=A(K)-T
            A(K)=A(K)+T
            T=A(K2+1)
            A(K2+1)=A(K+1)-T
            A(K+1)=A(K+1)+T
            T=A(K3)
            A(K3)=A(K1)-T
            A(K1)=A(K1)+T
            T=A(K3+1)
            A(K3+1)=A(K1+1)-T
            A(K1+1)=A(K1+1)+T
            T=A(K1)
            A(K1)=A(K)-T
            A(K)=A(K)+T
            T=A(K1+1)
            A(K1+1)=A(K+1)-T
            A(K+1)=A(K+1)+T
            R=-A(K3+1)
            T=A(K3)
            A(K3)=A(K2)-R
            A(K2)=A(K2)+R
            A(K3+1)=A(K2+1)-T
            A(K2+1)=A(K2+1)+T
          end do
        end do
        if (JLAST) 310,310,150
150     JJ=JJDIF+1
        ILAST=IL+JJ
        do I=JJ,ILAST,IDIF
          KLAST=KL+I
          do K=I,KLAST,2
            K1=K+KBIT
            K2=K1+KBIT
            K3=K2+KBIT
            R=-A(K2+1)
            T=A(K2)
            A(K2)=A(K)-R
            A(K)=A(K)+R
            A(K2+1)=A(K+1)-T
            A(K+1)=A(K+1)+T
            AWR = A(K1)-A(K1+1)
            AWI=A(K1+1)+A(K1)
            R=-A(K3)-A(K3+1)
            T=A(K3)-A(K3+1)
            A(K3)=(AWR-R)/ROOT2
            A(K3+1)=(AWI-T)/ROOT2
            A(K1)=(AWR+R)/ROOT2
            A(K1+1)=(AWI+T)/ROOT2
            T=A(K1)
            A(K1)=A(K)-T
            A(K)=A(K)+T
            T=A(K1+1)
            A(K1+1)=A(K+1)-T
            A(K+1)=A(K+1)+T
            R=-A(K3+1)
            T=A(K3)
            A(K3)=A(K2)-R
            A(K2)=A(K2)+R
            A(K3+1)=A(K2+1)-T
            A(K2+1)=A(K2+1)+T
          end do
        end do
        if (JLAST-1) 310,310,170
170     JJ=JJ+JJDIF
        do J=2,JLAST
          I=INV(J+1)
          IC=NT-I
          W(1)=S(IC)
          W(2)=S(I)
          I2=2*I
          I2C=NT-I2
          if (I2C) 200,190,180
180       W2(1)=S(I2C)
          W2(2)=S(I2)
          go to 210
190       W2(1)=0.
          W2(2)=1.
          go to 210
200       I2CC=I2C+NT
          I2C=-I2C
          W2(1)=-S(I2C)
          W2(2)=S(I2CC)
210       I3=I+I2
          I3C=NT-I3
          if (I3C) 240,230,220
220       W3(1)=S(I3C)
          W3(2)=S(I3)
          go to 280
230       W3(1)=0.
          W3(2)=1.
          go to 280
240       I3CC=I3C+NT
          if (I3CC) 270,260,250
250       I3C=-I3C
          W3(1)=-S(I3C)
          W3(2)=S(I3CC)
          go to 280
260       W3(1)=-1.
          W3(2)=0.
          go to 280
270       I3CCC=NT+I3CC
          I3CC=-I3CC
          W3(1)=-S(I3CCC)
          W3(2)=-S(I3CC)
280       ILAST=IL+JJ
          do I=JJ,ILAST,IDIF
            KLAST=KL+I
            do K=I,KLAST,2
              K1=K+KBIT
              K2=K1+KBIT
              K3=K2+KBIT
              R=A(K2)*W2(1)-A(K2+1)*W2(2)
              T=A(K2)*W2(2)+A(K2+1)*W2(1)
              A(K2)=A(K)-R
              A(K)=A(K)+R
              A(K2+1)=A(K+1)-T
              A(K+1)=A(K+1)+T
              R=A(K3)*W3(1)-A(K3+1)*W3(2)
              T=A(K3)*W3(2)+A(K3+1)*W3(1)
              AWR=A(K1)*W(1)-A(K1+1)*W(2)
              AWI=A(K1)*W(2)+A(K1+1)*W(1)
              A(K3)=AWR-R
              A(K3+1)=AWI-T
              A(K1)=AWR+R
              A(K1+1)=AWI+T
              T=A(K1)
              A(K1)=A(K)-T
              A(K)=A(K)+T
              T=A(K1+1)
              A(K1+1)=A(K+1)-T
              A(K+1)=A(K+1)+T
              R=-A(K3+1)
              T=A(K3)
              A(K3)=A(K2)-R
              A(K2)=A(K2)+R
              A(K3+1)=A(K2+1)-T
              A(K2+1)=A(K2+1)+T
            end do
          end do
          JJ=JJDIF+JJ
        end do
310     JLAST=4*JLAST+3
320     continue
330     continue
      end do
      NTSQ=NT*NT
      M3MT=M3-MT
      if (M3MT) 350,340,340
340   IGO3=1
      N3VNT=N3/NT
      MINN3=NT
      go to 360
350   IGO3=2
      N3VNT=1
      NTVN3=NT/N3
      MINN3=N3
360   JJD3=NTSQ/N3
      M2MT=M2-MT
      if (M2MT) 380,370,370
370   IGO2=1
      N2VNT=N2/NT
      MINN2=NT
      go to 390
380   IGO2=2
      N2VNT=1
      NTVN2=NT/N2
      MINN2=N2
390   JJD2=NTSQ/N2
      M1MT=M1-MT
      if (M1MT) 410,400,400
400   IGO1=1
      N1VNT=N1/NT
      MINN1=NT
      go to 420
410   IGO1=2
      N1VNT=1
      NTVN1=NT/N1
      MINN1=N1
420   JJD1=NTSQ/N1
      JJ3=1
      J=1
      do JPP3=1,N3VNT
        IPP3=INV(JJ3)
        do JP3=1,MINN3
          go to (430,440), IGO3
430       IP3=INV(JP3)*N3VNT
          go to 450
440       IP3=INV(JP3)/NTVN3
450       I3=(IPP3+IP3)*N2
          JJ2=1
          do JPP2=1,N2VNT
            IPP2=INV(JJ2)+I3
            do JP2=1,MINN2
              go to (460,470), IGO2
460           IP2=INV(JP2)*N2VNT
              go to 480
470           IP2=INV(JP2)/NTVN2
480           I2=(IPP2+IP2)*N1
              JJ1=1
              do JPP1=1,N1VNT
                IPP1=INV(JJ1)+I2
                do JP1=1,MINN1
                  go to (490,500), IGO1
490               IP1=INV(JP1)*N1VNT
                  go to 510
500               IP1=INV(JP1)/NTVN1
510               I=2*(IPP1+IP1)+1
                  if (J-I) 520,530,530
520               T=A(I)
                  A(I)=A(J)
                  A(J)=T
                  T=A(I+1)
                  A(I+1)=A(J+1)
                  A(J+1)=T
530               continue
                  J=J+2
                end do
                JJ1=JJ1+JJD1
              end do
            end do
            JJ2=JJ2+JJD2
          end do
        end do
        JJ3=JJ3+JJD3
      end do
      if (IFSET) 580,600,600
580   do I=1,NX
        A(2*I)=-A(2*I)
      end do
600   return
      end
C**********************************************
      subroutine RFFT (A,M,INV,S,IFERR,IFSET)
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      dimension A(1), L(3), INV(1), S(1)
C
C     IFSET=1
      L(1)=M
      L(2)=0
      L(3)=0
      NTOT=2**M
      NTOT2=2*NTOT
      FN=NTOT
      do I=2,NTOT2,2
        A(I)=-A(I)
      end do
      do I=1,NTOT2
        A(I)=A(I)/FN
      end do
      call FFT (A,L,INV,S,IFSET,IFERR)
C
C     MOVE LAST HALF OF A(J)S DOWN ONE SLOT AND ADD A(N) AT BOTTOM TO
C     GIVE ARRAY FOR A1PRIME AND A2PRIME CALCULATION
C
      do 30 I=1,NTOT,2
      J0=NTOT2+2-I
      A(J0)=A(J0-2)
30    A(J0+1)=A(J0-1)
      A(NTOT2+3)=A(1)
      A(NTOT2+4)=A(2)
C
C     CALCULATE A1PRIMES AND STORE IN FIRST N SLOTS
C     CALCULATE A2PRIMES AND STORE IN SECOND N SLOTS IN REVERSE ORDER
C
      K0=NTOT+1
      do I=1,K0,2
        K1=NTOT2-I+4
        AP1RE=.5*(A(I)+A(K1))
        AP2RE=-.5*(A(I+1)+A(K1+1))
        AP1IM=.5*(-A(I+1)+A(K1+1))
        AP2IM=-.5*(A(I)-A(K1))
        A(I)=AP1RE
        A(I+1)=AP1IM
        A(K1)=AP2RE
        A(K1+1)=AP2IM
      end do
      NTO=NTOT/2
      NT=NTO+1
      DEL=3.1415927/FLOAT(NTOT)
      SS=SIN(DEL)
      SC=COS(DEL)
      SI=0.0
      CO=1.0
C
C     COMPUTE C(J)S FOR J=0 THRU J=N
      do I=1,NT
        K6=NTOT2-2*I+5
        AP2RE=A(K6)*CO+A(K6+1)*SI
        AP2IM=-A(K6)*SI+A(K6+1)*CO
        CIRE=.5*(A(2*I-1)+AP2RE)
        CIIM=.5*(A(2*I)+AP2IM)
        CNIRE=.5*(A(2*I-1)-AP2RE)
        CNIIM=.5*(A(2*I)-AP2IM)
        A(2*I-1)=CIRE
        A(2*I)=CIIM
        A(K6)=CNIRE
        A(K6+1)=-CNIIM
        SIS=SI
        SI=SI*SC+CO*SS
        CO=CO*SC-SIS*SS
      end do
C
C     SHIFT C(J)S FOR J=N/2+1 to J=N UP ONE SLOT
      do I=1,NTOT,2
        K8=NTOT+4+I
        A(K8-2)=A(K8)
        A(K8-1)=A(K8+1)
      end do
      do I=3,NTOT2,2
        A(I)=2.*A(I)
        A(I+1)= 2.*A(I+1)
      end do
      return
      end
C******************************************
      subroutine RFSN (A,M,INV,S,IFERR,IFSET)
C**********************************************************************
      dimension A(1),L(3),INV(1),S(1)
C
      L(1)=M
      L(2)=0
      L(3)=0
      NTOT=2**M
C     IFSET=-1
      NTOT2=NTOT+NTOT
      NN=NTOT2+2
      A(NN+2)=A(NN)
      A(NN+1)=A(NN-1)
      FN=NTOT
      NTOT3=NTOT2+4
      do I=3,NTOT2,2
        A(I)=0.5* A(I)
        A(I+1)= .5*A(I+1)
      end do
      do I=1,NTOT,2
        K8=NTOT2+2-I
        A(K8)= A(K8-2)
        A(K8+1)=A(K8-1)
      end do
      NTO=NTOT/ 2
      NT=NTO+1
      DEL=3.141592654/FN
      SS= SIN(DEL)
      SC= COS(DEL)
      SI=0.
      CO =1.0
      do I=1,NT
        K6=NTOT2-2*I+5
        CIRE= A(2*I-1) + A(K6)
        CIIM=A(2*I)-A(K6+1)
        CNIRE=(-SI*(A(2*I)+A(K6+1))+CO*(A(2*I-1)-A(K6)))
        if (SI)62,61,62
  62    CNIIM=(A(2*I-1)-A(K6)-CO*CNIRE)/SI
        go to 63
  61    CNIIM=0.
  63    A(2*I-1)=CIRE
        A(2*I)=CIIM
        A(K6)=CNIRE
        A(K6+1)=CNIIM
        SIS=SI
        SI=SI*SC+CO*SS
        CO=CO*SC-SIS*SS
      end do
      KO=NTOT+1
      do I=1,KO,2
        K1=NTOT2-I+4
        AP1RE=A(I)-A(K1+1)
        AP2RE=-(A(I+1)+A(K1))
        AP1IM=A(I)+A(K1+1)
        AP2IM=A(I+1)-A(K1)
        A(I)=AP1RE
        A(I+1)=AP2RE
        A(K1)=AP1IM
        A(K1+1)=AP2IM
      end do
      NTOP=NTOT2+2
      NT00=NTOT+1
      A(1)=A(NTOT2+3)
      A(2)=A(NTOT2+4)
   21 do I=NT00,NTOP,2
        A(I)=A(I+2)
        A(I+1)=A(I+3)
      end do
      call  FFT(A,L,INV,S,IFSET,IFERR)
      do I=1,NTOT2
        A(I)=A(I)*FN
      end do
      do I=2,NTOT2,2
        A(I)=-A(I)
      end do
      return
      end
C*************************************
      subroutine XMX(X,MX,XMAX,NXMAX)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C   THIS ROUTINE FIND MAX. VALUE, XMAX, AND NUMBER OF MAX. VALUE, NXMAX.
C   OF ARRAY X WITH MX NUMBER OF VALUES
C
C   CODED PER B SCHNABEL OCT. 1971
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      dimension X(1)
      XMAX = 0.
      do I = 1,MX
        XA = ABS(X(I))
        if (XMAX .le. XA) then
          NXMAX = I
          XMAX = XA
        end if
      end do
      return
      end
