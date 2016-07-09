C $NODEBUG
C ......................................................................
      SUBROUTINE SHAKIT(X,AX,AA,S,INV)
C ......................................................................
      INTEGER TP
      CHARACTER*6 TITLE,ID,IDNT,IDAMP,IBLANK
      CHARACTER*60 ABSIS, ABSPR,ABSCL
      CHARACTER*80 OPHEAD
      CHARACTER*30 FINPEQ
      COMPLEX X, AX
      COMPLEX G, V, PLUS, MINUS
C
      DIMENSION LL(3), LT(3),LNSW(3)
      DIMENSION LLL(2), LLGS(2),LLPCH(2),LLPL(2),LNV(2),SK(2)
      DIMENSION X(300),AX(3,270),AA(2,550),S(70),INV(70)
      DIMENSION LL5(15), LT5(15), LP5(15),LP(3)
      DIMENSION IDAMP(27,11),MMM(3)
C
      COMMON /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
      COMMON /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      COMMON /SOILB/ FAC(51), WL(51), TP(51), DEPTH(51), WEIGHT(51)
      COMMON /SOILC/ MSOIL,MWL
      COMMON /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
      COMMON /CCG/ ID(27,11)
      COMMON /JISCK/ JIS,FINPEQ
      COMMON/FRCUT/ NCUT,NZERO
      COMMON /TIME/ T(9)
C     originallly coded by Per Schnabel in 1970-71
C     modified by Sun, Dirrim & Idriss in 1990-91 to
C     increase number of layers to 50;
C     renumber the Options & other cleanup
C
      IBLANK = '      '
      ABSIS = ' TIME IN SECONDS  '
      ABSCL = '  CYCLES/SEC'
      ABSPR = ' PERIOD IN SEC.   '
c
C * * * * * * * * *
C
      DO 102 I = 1,3
      LL(I) = 0
  102 LT(I) = 0
      DO 103 L = 1,9
      DO 103 I = 1,11
      ID(L,I) = IBLANK
  103 IDAMP(L,I) = IBLANK
      NF = 0
      NR=0
      NP=0
      NA = 1
C
C * * * * * * * * *
      KK = -1
C ...................................
  101 READ (5,1700) OPHEAD
 1700 FORMAT(A80)
C
      IF(KK.GE.1.AND.KK.LE.11) WRITE(*,24) KK
      READ(5,1000,END=999) KK
      IF (KK .EQ. 0) STOP
      WRITE(*,23) KK
C ******************** Options ****************
      IF (KK .EQ. 1) GO TO 8
      IF (KK .EQ. 2) GO TO 2
      IF (KK .EQ. 3) GO TO 1
      IF (KK .EQ. 4) GO TO 3
      IF (KK .EQ. 5) GO TO 4
      IF (KK .EQ. 6) GO TO 5
      IF (KK .EQ. 7) GO TO 16
      IF (KK .EQ. 8) GO TO 6
      IF (KK .EQ. 9) GO TO 9
      IF (KK .EQ. 10) GO TO 15
      IF (KK .EQ. 11) GO TO 13
C * * * * * * * * *
    1 WRITE(6,1002) KK
      CALL EARTHQ(X,AX, S, INV)
      NSN = 0
      GO TO 101
C * * * * * * * * *
    2 WRITE(6,2002) KK
      CALL SOILIN(N1)
      NSN = 1
C
C   FIND FUNDAMENTAL PERIOD OF DEPOSIT FROM AVERAGE SHEAR WAVE VELOCITY
C   AND FROM THE PERIOD WHICH GIVE MAXIMUM AMPLIFICATION
      SH = 0.
      N = N1 + 1
      SHV = 0.
      DO 21 I = 1,N1
      SH = SH + H(I)
   21 SHV = SHV + H(I)*SQRT(GL(I)/R(I))
      VSAV = SHV/SH
      TT = 4.*SH/VSAV
      WRITE(6,4006) TT, VSAV
      DFA = .01/TT
      CALL AMP(N1,N ,1,1,0,0,IDAMP, 9, DFA)
      GO TO 101
C * * * * * * * * *
    3 WRITE(6,3002) KK
      READ(5,1000) IN, INT
      IF (INT .EQ. 0) WRITE(6,3001) IN
      IF (INT .NE. 0) WRITE(6,3000) IN
      GO TO 101
C * * * * * * * * *
    4 WRITE(6,4007) KK
      READ(5,4000) KS,ITMAX,PRMUL
      WRITE(6,4001) ITMAX,PRMUL
      LL(1) = 1
      LT(1) = 0
      JIS = 0
      WRITE(*,2029)
 2029 FORMAT(/)
      DO 41 L = 1,ITMAX
      WRITE(*,2028) L
 2028 FORMAT(1H+,12X,19H  ITERATION NUMBER  , I2)
      IF (IN.EQ.1) GO TO 412
      CALL MOTION(N1,IN, INT, LL, LT, X,AX)
      IF (L .EQ. ITMAX) JIS = 1
  412 CALL STRT( L, N1, DGMAX,PRMUL,X,AX,AA,S,INV)
C     IF (DGMAX.LT.ERR) GO TO 411
   41 CONTINUE
C
C   FIND FUNDAMENTAL PERIOD OF DEPOSIT FROM AVERAGE SHEAR WAVE VELOCITY
C   AND FROM THE PERIOD WHICH GIVE MAXIMUM AMPLIFICATION
  411 SH = 0.
      N = N1 + 1
      SHV = 0.
      DO 43 I = 1,N1
      SH = SH + H(I)
   43 SHV = SHV + H(I)*SQRT(GL(I)/R(I))
      VSAV = SHV/SH
      TT = 4.*SH/VSAV
      WRITE(6,4006) TT, VSAV
      DFA = .01/TT
      CALL AMP(N1,N ,1,1,0,0,IDAMP, 9, DFA)
C
      IF (KS .EQ. 0) GO TO 101
C   SAVE NEW SET OF SOIL DATA BASED ON NEW PROPERTIES
      WRITE(7,4003) MSOIL,N,MWL,(IDNT(I),I=1,6),(TITLE(I),I=1,4)
      DO 42 I = 1,N1
   42 WRITE(7,4004) I,TP(I), H(I), GL(I), BL(I), WL(I), FAC(I) ,BF(I)
      WRITE(7,4005) N,GL(N),BL(N),WL(N)
      GO TO 101
C * * * * * * * * *
    5 WRITE(6,5001) KK
      READ(5,1000) (LL5(L),L=1,15)
      READ(5,1000) (LT5(L),L=1,15)
      READ(5,1000) (LP5(L),L=1,15)
      WRITE(6,5002) FINPEQ,(IDNT(I),I=1,6)
      I = 0
      DO 51 LOOP = 1,5
      DO 511 L = 1,3
      I = I + 1
      LL(L) = LL5(I)
      LT(L) = LT5(I)
      LP(L) = LP5(I)
C
      IF (LL(1).EQ.0) GO TO 101
  511 CONTINUE
C
      CALL MOTION(N1,IN, INT, LL, LT, X,AX)
      DO 51  L = 1,3
      N = LL(L)
      K = L
      IF (N.EQ.0) GO TO 101
      IF (N.LE.N1)  DPTH = DEPTH(N) - H(N)/2.
      IF (N.GT.N1)  DPTH = DEPTH(N-1) + H(N-1)/2.
      CALL UTPR(KK,DPTH,K,LP(L),LL(L),LT(L),X,AX,S,INV)
   51 CONTINUE
      GO TO 101
C * * * * * * * * *
    6 WRITE(6,6002) KK
      READ(5,1000) K2
      LS = 0
      LN = IN
      IF (K2.EQ.0)  WRITE(6,6000) LN
      IF (K2.EQ.1)  WRITE(6,6001) LN
   62 CALL UTPR(KK,DPTH,LS,K2,LN,INT   ,X,AX,S,INV)
      GO TO 101
C * * * * * * * * *
    7 WRITE(6,7002) KK
      READ(5,7001) LL1,LT1,XF,DTNEW
      IF (DTNEW.LT..001) DTNEW=DT
      IF (LL1  .EQ.0) GO TO 71
C   CHECK IF MOTION IN SUBLAYER LL1 IS IN AX()
      DO 72 I = 1,3
      IF (LL1.NE.LL(I) .OR. LT1.NE.LT(I)) GO TO 72
      L = I
      GO TO 720
   72 CONTINUE
      LL(1) = LL1
      LT(1) = LT1
      L = 1
      CALL MOTION(N1,IN, INT, LL, LT, X,AX)
  720 DO 75 I = 1,MFOLD
   75 X(I) = AX(L,I)*XF
      NEW = LL(L)
      INT = LT(L)
      GO TO 73
   71 DO 74 I = 1,MFOLD
   74 X(I) = X(I)*XF
      NEW = IN
   73 IN = NEW
      WRITE(6,7000) NEW , XF,DT, DTNEW
      IF(IN.NE.1) GO TO 76
      DO 77 II=1,MFOLD
      AX(1,II)=X(II)
   77 CONTINUE
   76 CONTINUE
      DT = DTNEW
      DF = 1./(MA*DT)
      GO TO 101
C * * * * * * * * *
    8 WRITE(6,8001) KK
      CALL CG
      GO TO 101
C * * * * * * * * *
    9 WRITE(6,9002) KK
      READ(5,1000) LL1, LT1
      IF (LL1.NE.0) GO TO 171
      WRITE(6,9001)
      LS = 0
      LN = IN
      GO TO 173
  171 DO 170 I = 1,3
      IF (LL1.NE.LL(I) .OR. LT1.NE.LT(I)) GO TO 170
      LS = I
      GO TO 172
  170 CONTINUE
      LS = 1
      LL(1) = LL1
      LT(1) = LT1
      CALL MOTION(N1,IN, INT, LL, LT, X,AX)
  172 LN = LL(LS)
C      WRITE(6,9000) LN
  173 CALL RESP(LN,LS,NR,X,AX,AA,S,INV)
      GO TO 101
C * * * * * * * * *
   10 WRITE(6,1010) KK
      READ(5,1000) IFR
      CALL REDUCE(IFR,X,AX,LL)
      MMM(1)=MX
      MMM(2)=0
      MMM(3)=0
      CALL FFT(X,MMM,INV,S,0,IFERR)
      GO TO 101
C * * * * * * * * *
   11 WRITE(6,1101) KK
      READ(5,1000) IFR
      CALL INCR(IFR,X,AX)
      MMM(1)=MX
      MMM(2)=0
      MMM(3)=0
      CALL FFT(X,MMM,INV,S,0,IFERR)
      GO TO 101
C * * * * * * * * *
   12 WRITE(6,1203) KK
      READ(5,1000) K1, NSW, N
      IF (INT.EQ. 0) WRITE(6,1201) IN
      IF (INT.EQ. 1) WRITE(6,1202) IN
      NF = NF + 1
      IF (N.LE.0)  N= MFOLD - 1
      IF (N.GT.2049)   N=2049
      DO 120 I = 1,N
  120 AA(NF,I) = CABS(X(I))
      DO 121 I = 1,5
  121 ID(NF,I) = TITLE(I)
      DO 126 I = 6,11
  126 ID(NF,I) = IBLANK
      IF (NSW.EQ.0) GO TO 123
      M = N-1
      DO 124 LOOP = 1,NSW
      AA(NF,1) = (3.*AA(NF,1) + AA(NF,2))/4.
      AA(NF,N) = (3.*AA(NF,N) + AA(NF,N-1))/4.
      DO 124 I = 2,M
  124 AA(NF,I) = (AA(NF,I-1) + 2.*AA(NF,I) + AA(NF,I+1))/4.
  123 IF (K1.NE.1) GO TO 101
      DO 122 I = 1,N
      T(I ) = FLOAT(I-1)*DF
  122 CONTINUE
      WRITE(6,1200)
      NF = 0
      GO TO 101
C * * * * * * * * *
   13 WRITE(6,1301) KK
      DO 180 I = 1,2
  180 READ(5,1000) LL(I), LT(I), LP(I), LNSW(I), LLL(I)
      CALL MOTION(N1, IN, INT, LL, LT, X, AX)
      NF = 0
      DO 184 L = 1,2
      IF (LL(L).EQ.0) GO TO 101
      IF (LT(L).EQ.0) WRITE(6,1201) LL(L)
      IF (LT(L).EQ.1) WRITE(6,1202) LL(L)
      N = LLL(L)
      IF (N.LE.0)   N = MFOLD - 1
      IF (N.GT.2049)   N = 2049
      NF = NF + 1
      IF (NF.LE.2)  GO TO 182
      WRITE(6,1800)
 1800 FORMAT(//  24H TOO MANY ARRAYS STORED   /)
      GO TO 101
  182 DO 188 I = 1,5
  188 ID(NF,I) = TITLE(I)
      DO 187 I = 6,11
  187 ID(NF,I) = IDNT(I-5)
      DO 185 I = 1,N
  185 AA(NF,I) = CABS(AX(L,I))
      NSW = LNSW(L)
      IF (NSW.EQ.0) GO TO 181
      M = N-1
      DO 186 LOOP = 1,NSW
      AA(NF,1) = (3.*AA(NF,1) + AA(NF,2))/4.
      AA(NF,N) = (3.*AA(NF,N) + AA(NF,N-1))/4.
      DO 186 I = 2,M
  186 AA(NF,I) = (AA(NF,I-1) + 2.*AA(NF,I) + AA(NF,I+1))/4.
      IF (LP(L).EQ.0) GO TO 184
  181 DO 183 I = 1,N
  183 T(I) = DF*FLOAT(I-1)
C     WRITE(6,1200)
  184 CONTINUE
      WRITE(6,1204)
 1204 FORMAT(1X,'       FREQ       FOURIER AMPLITUDES')
 1205 FORMAT(1X,F10.4,2F15.6)
      DO 133  I=1,N
  133 WRITE(6,1205)  T(I),(AA(NF,I), NF=1, 2)
      GO TO 101
C * * * * * * * * *
   14 WRITE(6,1404) KK
      READ(5,1000) NSKIP, NN, NSW
      NP = NP + 1
      CALL RFSN(X,MX,INV,S,IFERR,-2)
      IF (NN.LE.0)  NN = MMA/NSKIP
      IF (NN.GT.2049)  NN = 2049
      NN =  NN*NSKIP
      N = 0
      DO 136 I=1, NN, NSKIP
      N = N + 1
      T(N) = FLOAT(I-1)*DT
  136 CONTINUE
      N = 0
      M = NN/2
      DO 130 I = 1,M
      N = N + 1
      AA(NP,N) = REAL(X(I))
      N = N + 1
      AA(NP,N) = AIMAG(X(I))
  130 CONTINUE
      IF  (NSKIP.EQ.1) GO TO 135
      N = 0
      DO 134 I = 1,NN ,NSKIP
      N = N + 1
      AA(NP,N) = AA(NP,I)
  134 CONTINUE
  135 CALL RFFT(X,MX,INV,S,IFERR,2)
      DO 131 I = 1,5
  131 ID(NP,I) = TITLE(I)
      DO 132 I = 6,11
      ID(NP,I) = IDNT(I-5)
      IF (NSN.EQ.0)  ID(NP,I)  = IBLANK
  132 CONTINUE
      IF (NSW.EQ.1) GO TO 101
      NP = 0
      GO TO 101
C * * * * * * * * *
   15 WRITE(6,1502) KK
      READ(5,1400) LIN, LINT ,LOUT,LOTP,DFA,(IDAMP(NA,I),I=1,6)
      KP = 2
      WRITE(6,1401) LIN,LOUT
      IF  (LOTP.EQ.0) WRITE(6,1403)
      IF  (LINT.EQ.0) WRITE(6,1402)
      CALL AMP(N1, LIN, LINT, LOUT, LOTP, KP, IDAMP,NA, DFA)
      GO TO 101
C * * * * * * * * *
   16 WRITE(6,1601) KK
      DO 151 L = 1,2
      READ(5,1500) LLL(L),LLGS(L),LLPCH(L),LLPL(L),LNV(L),SK(L),
     1(ID(L,I),I=7,11)
      IF (LLL(L).GT.0) WRITE(6,1501) LLL(L),SK(L),(ID(L,I),I=7,11)
  151 CONTINUE
      DO 152 L = 1,3
  152 LT(L) = 0
      LL(3) = 0
      LL(2) = 0
      LL(1) = 1
      CALL MOTION(N1, IN, INT, LL, LT, X, AX)
      CALL STRAIN(LLL,LLGS,LLPCH,LLPL,LNV,X,AX,AA,N1,S,INV)
      DO 153 I = 1,3
  153 LL(I) = 0
      GO TO 101
C * * * * * * * * *
C
   23 FORMAT(5X,'Option NO.',I5,'   is started.')
   24 FORMAT(5X,'Option NO.',I5,'   has been concluded.')
 1000 FORMAT(15I5)
 1002 FORMAT(/16H1******   OPTION  I3,
     1 58H  ***  READ INPUT MOTION                                   )
 2002 FORMAT(/16H1******   OPTION  I3,
     1 58H  ***  READ SOIL PROFILE                                   )
 3000 FORMAT(/32H  OBJECT MOTION IN LAYER NUMBER   I3/)
 3001 FORMAT(33H   OBJECT MOTION IN LAYER NUMBER I3,12H OUTCROPPING )
 3002 FORMAT(16H1******   OPTION  I3,
     1 58H  ***  READ WHERE OBJECT MOTION IS GIVEN                     )
 4000 FORMAT( 2I5, 7F10.0)
 4001 FORMAT(
     148H  MAXIMUM NUMBER OF ITERATIONS                =    ,I5/
     148H  FACTOR FOR UNIFORM STRAIN IN TIME DOMAIN    =    ,F6.2/)
 4003 FORMAT(3I5,6A6, 4A6)
 4004 FORMAT(2I5,4X,1H1,F10.2,F10.0,2F10.3,10X,F10.3, F5.2)
 4005 FORMAT(I5,9X,1H110X,F10.0,2F10.3)
 4006 FORMAT(/10H  PERIOD = F5.2,' FROM AVERAGE SHEAR VELOCITY =',F8.0/)
 4007 FORMAT(/16H1******   OPTION   I3,
     1 58H  ***  OBTAIN STRAIN COMPATIBLE SOIL PROPERTIES            )
 5001 FORMAT(/16H1******   OPTION  I3,
     1 58H  ***  COMPUTE MOTION IN NEW SUBLAYERS                     )
 5002 FORMAT(/15H   EARTHQUAKE - A30/  17H  SOIL DEPOSIT -  6A6/
     1  5X,6HLAYER ,
     2 10X,5HDEPTH ,8X,9HMAX. ACC. 10X, 4HTIME ,6X, 12HMEAN SQ. FR. ,
     3  9X, 10HACC. RATIO , 6X, 14H TH SAVED		,
     4/22X, 2HFT ,12X, 1HG ,16X, 4HSEC  ,9X, 5HC/SEC ,13X,10HQUIET ZONE 
     5  7X, 11HACC. RECORD   )
 6000 FORMAT(/37H  PRINT        ACCELERATION IN LAYER  ,I3)
 6001 FORMAT(/46H  PRINT AND PUNCH        ACCELERATION IN LAYER ,I3)
 6002 FORMAT(16H1******   OPTION  I3,
     1 58H  ***  PRINT OR PUNCH OBJECT MOTION                        )
 7000 FORMAT(/21H  SET MOTION IN LAYER  ,I3,17H AS OBJECT MOTION /
     141H MULTIPLICATION FACTOR FOR NEW MOTION =  ,F6.3/
     227H  TIMESTEP DT CHANGED FROM , F6.3, 3H TO, F6.3, 5H SEC./)
 7001 FORMAT(2I5, 5F10.0)
 7002 FORMAT(16H1******   OPTION   I3,
     1 58H  ***  CHANGE OBJECT MOTION                                  )
 8001 FORMAT(/16H1******   OPTION   I3,
     1 58H  ***  READ RELATION BETWEEN SOIL PROPERTIES AND STRAIN    )
 9000 FORMAT( 36H  COMPUTE RESPONSE SPECTRUM IN LAYER  I3)
 9001 FORMAT(43H COMPUTE RESPONSE SPECTRUM OF OBJECT MOTION )
 9002 FORMAT(/16H1******   OPTION  I3,
     1 58H  ***  COMPUTE RESPONSE SPECTRUM                             )
 1010 FORMAT(16H1******   OPTION    I3,
     1 58H  ***  INCREASE TIMESTEP                                   //)
 1101 FORMAT(16H1******   OPTION    I3,
     1 58H  ***  DECREASE TIMESTEP                                   //)
 1200 FORMAT(26H1 FOURIER SPECTRA        ')
 1201 FORMAT(15H  LAYER NUMBER I4, 12H OUTCROPPING /)
 1202 FORMAT(14H  LAYER NUMBER I4)
 1203 FORMAT(16H1******   OPTION    I3,
     1 58H  ***  PLOT OF FOURIER SPECTRUM OF OBJECT MOTION           )
 1301 FORMAT(16H1******   OPTION    I3,
     1 58H  ***       FOURIER SPECTRUM OF COMPUTED MOTION            )
 1400 FORMAT(4I5, F10.0, 6A6)
 1401 FORMAT(/41H    AMPLIFICATION SPECTRUM BETWEEN LAYER I4, 4H AND I4)
 1402 FORMAT(26H  INPUT LAYER OUTCROPPING )
 1403 FORMAT(26H  OUTPUT LAYER OUTCROPPING )
 1404 FORMAT(16H1******   OPTION  I3,
     1 58H  ***  PLOT TIME HISTORY OF OBJECT MOTION                  )
 1500 FORMAT(5I5,F10.0, 5A6)
 1501 FORMAT(/  49H  COMPUTE STRESS OR STRAIN HISTORY AT THE TOP OF
     1 6H LAYER  ,I5 /21H  SCALE FOR PLOTTING  ,F10.4/15H IDENTIFICATION
     2 3H -    ,5A6,6X)
 1502 FORMAT(/16H1******   OPTION  I3,
     1 58H  ***  COMPUTE AMPLIFICATION FUNCTION                      )
 1601 FORMAT(/16H1******   OPTION  I3,
     1 58H  ***  COMPUTE STRESS/STRAIN HISTORY                       )
C
  999 STOP
      END
C*********************************************************************
      SUBROUTINE STRT( IT,N1,DGMAX,PRMUL,X,AX,AA,SF,INV)
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  *
C
C   THIS ROUTINE CALCULATES STRAIN IN THE MIDDLE OF EACH LAYER AND FIND
C   NEW SOIL PROPERTIES COMPATIBLE WITH THE STRAINS
C
C        IT      = ITERATION NUMBER
C        N1      = NUMBER OF LAYERS EXCLUDING ROCK
C        DGMAX   = MAX ERROR IN SOIL PARAMETERS B OR G IN PERCENT
C        X       = OBJECT MOTION
C        AX(1, ) = ACCELERATION VALUES AT THE SURFACE
C        AX(2, ) = INCIDENT WAVE-COMPONENT
C        AX(3, ) = REFLECTED WAVE-COMPONENT
C        PRMUL   = RATIO EFF. STRAIN/MAX. STRAIN
C
C   CODED PER B SCHNABEL OCT. 1970
C   MODIFIED PBS  SEPT. 1971
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  *
      INTEGER TP
      CHARACTER*6 TITLE,IDNT
      CHARACTER*30 FINPEQ
      COMPLEX IPI2, EX, E, F, EE, FF
      COMPLEX X, AX
      COMPLEX G, V, PLUS, MINUS
C
      DIMENSION TMAX(51),EMAX(51),STR(51)
      DIMENSION X( 68), AX(3, 64), AA(2,128),SF(10), INV(10), ratio(51)
      COMMON /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF,MX
      COMMON /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      COMMON /SOILB/ FAC(51), WL(51), TP(51), DEPTH(51), WEIGHT(51)
      COMMON /SOILC/ MSOIL,MWL
      COMMON /SOILD/ GLMAX(51)
      COMMON /SOILDG/ S(27,20), AS(27,20), BS(27,20), NV(27)
      COMMON /CSOIL/ G(51), V(51), PLUS(51), MINUS(51)
      COMMON /JISCK/ JIS,FINPEQ
      COMMON/FRCUT/ NCUT,NZERO
C
      DO 43 I = 1,MFOLD
      AA(1,I) = REAL(X(I))
   43 AA(2,I) = AIMAG(X(I))
      DO 1  I = 1,MFOLD
      AX(2,I) = AX(1,I)/2.
    1 AX(3,I) = AX(2,I)
      PI2=6.283185307
      IPI2=CMPLX(0.,PI2)
      GT = 32.2
      DO 2 K = 1,N1
C     WRITE(*,2029) K
      FREQ = 0.
      X(1) = 0.
      FF = GT/(IPI2*V(K))
      EE = H(K)/2.*IPI2/V(K)
      DO 20 I=2,NCUT
      FREQ = FREQ + DF
      EX = CEXP(FREQ*EE)
      X(I) = (AX(2,I)*EX - AX(3,I)/EX)*FF/FREQ
      EX = EX*EX
      E = AX(2,I)*EX
      F = AX(3,I)/EX
      AX(2,I)= PLUS(K)*E + MINUS(K)*F
      AX(3,I)= PLUS(K)*F + MINUS(K)*E
   20 CONTINUE
      EMAX(K) = 0.
      IF(NCUT.EQ.MFOLD) GO TO 22
      DO 122 II=NZERO,MFOLD
      X(II)=CMPLX(0.,0.)
  122 CONTINUE
   22 CONTINUE
C
C   DETERMINE MAX. STRAIN BY INVERTING FOURIER TRANSFORM OF STRAIN
C   INTO THE TIME DOMAIN
C
      CALL RFSN(X,MX,INV,SF,IFERR,-2)
      CALL XMX(X,MA,XMAX,NXMAX)
C
      EMAX(K) = XMAX
      TMAX(K) = FLOAT(NXMAX-1)*DT
    2 CONTINUE
      IF (IT.GT.1) WRITE(6,2002)
      WRITE(6,2017) FINPEQ, (IDNT(I),I=1,6)
      WRITE(6,2027) IT
C     WRITE(6,2037) PRMUL
      WRITE(6,2000)
      DGMAX = 0.
      DO 23 I = 1,N1
      EM = EMAX(I)*PRMUL*100.
      EMAX(I) = EMAX(I)*100.
      IF (TP(I) .NE. 0)  GO TO 231
      STR(I) = EMAX(I)*GL(I)*10.
      WRITE(6,2107) I, TP(I), DEPTH(I), EM , BL(I), GL(I)
      GO TO 23
C
C   USE UNIFORM STRAIN AMPLITUDE (EM) TO GET NEW VALUES FOR DAMPING
C   AND SHEAR MODULUS
  231 IN = TP(I)*2 - 1
      SS= ABS(EM)
      SL = ALOG10(SS)
      LL = NV(IN)
      DO 31 L = 1,LL
      IF (SS.LE. S(IN,L)) GO TO 311
   31 CONTINUE
      L = LL
  311 GN =AS(IN,L)*SL +BS(IN,L)
      GG = GN*FACT(I)/1000.
      IN = IN + 1
      LL = NV(IN)
      DO 32 L = 1,LL
      IF (SS.LE. S(IN,L)) GO TO 321
   32 CONTINUE
      L = LL
  321 B  =AS(IN,L)*SL +BS(IN,L)
      B = B*BF(I)
C ----------------------------------------------------------------------
C  SHEAR STRESSES ARE COMPUTED USING CURRENT MODULI
C
      STR(I) = EMAX(I)*GL(I)*10.
      RATIO(I) = GL(I) / GLMAX(I)
C ----------------------------------------------------------------------
      B = B/100.
      DG = (GG - GL(I))*100./GG
      DB = ( B - BL(I))*100./B
      WRITE(6,2007) I, TP(I), DEPTH(I), EM, B, BL(I), DB, GG, GL(I), DG,
     + RATIO(I)
      IF (ABS(DG) .GT. DGMAX)  DGMAX = ABS(DG)
      IF (ABS(DB) .GT. DGMAX)  DGMAX = ABS(DB)
      IF (JIS .EQ. 1) GO TO 23
      BL(I) = B
      GL(I) = GG
   23 CONTINUE
      IF (JIS .NE. 1) GO TO 53
      WRITE(6,2011)
      WRITE(6,2001) (I,TP(I), H(I),DEPTH(I), EMAX(I), STR(I), TMAX(I),
     1 I = 1,N1)
   53 CALL CXSOIL(N1)
      DO 44 I = 1,MFOLD
   44 X(I) = CMPLX(AA(1,I),AA(2,I))
      RETURN
C
 2000 FORMAT(/23H  VALUES IN TIME DOMAIN //
     1,' NO TYPE DEPTH  UNIFRM. <---- DAMPING ---->   <---- SHEAR',
     2' MODULUS ----->    G/Go'/
     3'         (FT)   STRAIN  NEW    USED   ERROR      NEW       USED',
     4'    ERROR   RATIO'/
     5'--- ---- ----  ------- ----- ------  ------   -------   -------',
     6'   ------   -----')
 2002 FORMAT(1H1)
 2011 FORMAT(/23H  VALUES IN TIME DOMAIN //
     1 2X, 5HLAYER ,2X,4HTYPE ,6X,9HTHICKNESS ,10X, 5HDEPTH ,5X, 
     2 10HMAX STRAIN ,5X, 10HMAX STRESS ,10X, 4HTIME /
     3 23X, 2HFT 14X, 2HFT 9X, 5HPRCNT 12X, 3HPSF 13X, 3HSEC /)
 2001 FORMAT(2I6, 2F15.1, F15.5,2F15.2)
 2007 FORMAT(2I3, F7.1, F8.5, F6.3, 1X, F6.3 ,F9.1, 1X,
     1 F9.1,1X, F9.1, F9.1, F8.3)
 2107 FORMAT(2I3, F8.1, F9.5, 2F7.3, F8.1, 2F10.1, F8.1, F7.3)
 2027 FORMAT(19H  ITERATION NUMBER   I2)
 2029 FORMAT(1H+,20X,'Processing layer no. ',I2)
 2030 FORMAT(1H+'                                                    ')
 2017 FORMAT(5X,15HEARTHQUAKE   -  A30/5X,15HSOIL PROFILE -   6A6/)
 2037 FORMAT(56H  THE CALCULATION HAS BEEN CARRIED OUT IN THE TIME DOMAI
     121HN WITH EFF. STRAIN = ,F3.2, 1H* 12H MAX. STRAIN )
      END
C****************************************************************
      SUBROUTINE SOILIN(N1)
C***********************************************************************
C
C   THIS ROUTINE READS PROPERTIES OF A SOIL PROFILE, ASSIGNS VALUES TO
C   EACH LAYER, CALCULATES TOTAL PRESSURE AND DEPTH IN MIDDLE OF
C   EACH LAYER AND PRINTS THE RESULTS
C
C        IDNT    = IDENTIFIER FOR SOIL PROFILE
C        BL      = RATIO OF CRITICAL DAMPING
C        GL      = SHEAR MODULUS
C        FACT    = FACTOR FOR CALCULATING SHEAR MODULUS FROM STRAIN
C        H       = LAYER THICKNESS
C        R       = DENSITY
C        WL      = UNIT WEIGHT
C        TP      = SOIL TYPE
C        DEPTH   = DEPTH TO MIDDLE OF LAYER
C        WEIGTH  = TOTAL PRESSURE
C        ML      = NUMBER OF LAYERS INCLUDING HALFSPACE
C        N1      = NUMBER OF SUBLAYERS EXCLUDING HALFSPACE
C        NLN     = NUMBER OF SUBLAYERS IN EACH LAYER
C        W       = UNIT WEIGHT
C        VS      = SHEAR WAVE VELOCITY
C        BFAC    = FACTOR ON DAMPING
C        FACTOR  = FACTOR ON SHEAR MODULUS
C        HL      = THICKNESS OF LAYER
C        H       = THICKNESS OF SUBLAYER
C        GMOD    = SHEAR MODULUS
C        B       = CRITICAL DAMPING RATIO
C
C   CODED BY PER B SCHNABEL OCT. 1970
C   MODIFIED APRIL 1972
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      INTEGER TP, TYPE
      CHARACTER*6 IDNT
      DIMENSION SMEAN(51)
C
      COMMON /SOILA/ IDNT(6),BL(51),GL(51),FACT(51),H(51),R(51),BF(51)
      COMMON /SOILB/ FAC(51), WL(51), TP(51), DEPTH(51), WEIGHT(51)
      COMMON /SOILC/ MSOIL,MWL
      COMMON /SOILD/ GLMAX(51)
      COMMON /WGK/  WW, GT, SKO
C
      READ(5,1003) MSOIL, ML, MWL, (IDNT(I),I=1,6)
      WRITE(6,2020) MSOIL, (IDNT(I),I=1,6)
C
C READ SOIL PROPERTIES FOR EACH LAYER AND ASSIGN VALUES TO EACH SUBLAYER
C
      J = 0
      DO 14 N =1, ML
      READ(5,1004) K, TYPE, NLN,   HL, GMOD, B , W, VS
      FACTOR = 1.
       BFAC = 1.
      IF (NLN.EQ. 0) NLN = 1
      IF (K .EQ. N) GO TO 141
      WRITE(6,2004) N
      STOP
C
C   COMPUTE MODULUS FROM SHEAR WAVE VELOCITY
C
  141 IF (GMOD .EQ. 0.) GMOD = VS*VS*W/GT
      DO 14 I = 1,NLN
      J = J+1
      BL(J) = B
      GL(J) = GMOD
      GLMAX(J) = GMOD
      FAC(J) = 1.
      FACT(J) =  1.
      BF(J) = 1.
      WL(J) = W
      H(J) = HL/NLN
      TP(J) = TYPE
   14 R(J) = W/GT
      N1 = J -1
C
C   CALCULATE AVERAGE DEPTH AND TOTAL PRESSURE IN EACH LAYER
C
      W1 = WL(1)
      IF (MWL .EQ. 1) W1 = WL(1) - WW
      DEPTH(1) = H(1)/2.
      WEIGHT(1) = H(1)*W1/2.
      SMEAN(1) = WEIGHT(1)*(1.+2.*SKO)/3.
      IF  (N1 .EQ. 1)  GO TO 151
      DO 15 I = 2,N1
      W2 = WL(I)
      IF (MWL .LT. I+1) W2 = WL(I) - WW
      DEPTH(I) = DEPTH(I-1) + H(I)/2. + H(I-1)/2.
      WEIGHT(I) = WEIGHT(I-1) + H(I)*W2/2. + H(I-1)*W1/2.
      SMEAN(I) = WEIGHT(I)*(1.+2.*SKO)/3.
   15 W1 = W2
  151 TD = DEPTH(N1) + H(N1)/2.
      IF (MWL .LT. N1+1)  WD = DEPTH(MWL)- H(MWL)/2.
      IF (MWL .EQ. N1+1)  WD = DEPTH(MWL-1) + H(MWL-1)/2.
C
C   CALCULATE FACTOR FOR SHEAR MODULUS
      DO 16 I = 1,N1
      IF (TP(I) .EQ. 0) GO TO 16
      IF (BF(I).LT..01) BF(I) = 2.53 - .45*ALOG10(WEIGHT(I)*1000.)
      NTP = TP(I)
C ----------------------------------------------------------------------
C  A total of 13 G/Gmax material types can be used
c-----------------------------------------------------------------------
      FAC(I) = FACT(I)
      FACT(I) = GL(I) * 1000. * FACT(I)
c-----------------------------------------------------------------------
   16 CONTINUE
  131 WRITE(6,2021) ML,TD
      WRITE(6,2015)
      DO 17 I = 1,N1
      VS = SQRT( GL(I)/R(I))
      WRITE(6,2005) I, TP(I), H(I),DEPTH(I)
     1,WEIGHT(I),GL(I),BL(I),WL(I),VS
   17 CONTINUE
      I = N1 + 1
      VS = SQRT(GL(I)/R(I))
      WRITE(6,2105) I, GL(I), BL(I), WL(I), VS
      CALL CXSOIL(N1)
 1003 FORMAT(3I5, 6A6)
 1004 FORMAT(3I5, 6F10.0,F5.0)
 2004 FORMAT(17H   SOIL CARD NO. I4,17H OUT OF SEQUENCE )
 2020 FORMAT(22H NEW SOIL PROFILE NO. ,I3,5X,17H IDENTIFICATION   ,6A6)
 2021 FORMAT(17H NUMBER OF LAYERS ,I20,10X,16HDEPTH TO BEDROCK,F14.2/)
 2015 FORMAT(  '  NO. TYPE  THICKNESS   DEPTH  ',
     1   'Tot. PRESS.  MODULUS  DAMPING   UNIT WT.  SHEAR VEL' /
     3   '              (ft)      (ft)      (ksf)       (ksf)',
     4    '             (kcf)    (fps)')
 2005 FORMAT(I4,I5,F10.2,F10.2,F10.2,F12.0,F8.3,F9.3,F10.1)
 2105 FORMAT(    I4, 3X, 4HBASE  25X, F15.0, F8.3, F9.3, F10.1)
      RETURN
      END
C $NOFLOATCALLS
C $NODEBUG
C*****************************************************
      SUBROUTINE CG
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * * *
C
C   THE SUBROUTINE READ POINTS ON A CURVE AND GENERATES NEW POINTS
C   BETWEEN THE GIVEN POINTS IN ARITHMETIC OR HALFLOGARITMIC SCALE
C   NECESSARY SUBROUTINES    CURVEG(),
C
C        NST     = NUMBER OF SOILTYPES
C        ABSIS   = TITLE ON ORDINATE FOR PLOTTING
C        NN      = NUMBER OF VALUES IN EACH 10 FOR SEMILOGPLOT
C        SC      = SCALE FOR PLOTTING
C        NC      = NUMBER OF CURVES
C        NV      = NUMBER OF VALUES WHERE STRAIN/PROPERTY-RELATION
C                  IS GIVEN
C        FPL     = MULTIPLICATION FACTOR FOR PLOTTING
C        ID      = IDENTIFICATION
C        X       = STRAIN VALUES
C        Y       = PROPERTY VALUES
C
C   CODED BY PER B SCHNABEL SEPT 1970
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  * * *
C
      CHARACTER*6 ID
      CHARACTER*60 ABSIS
C
C     DIMENSION Y(27,20), TSTEP(27),NT(27),FPL(27),V(27,200),T(200)
      INTEGER ECHO(13)
      COMMON /JOE1/ Y(27,20), TSTEP(27),NT(27),FPL(27),V(27,200),T(200)
      COMMON /SOILDG/ S(27,20),AS(27,20),BS(27,20), NV(27)
      COMMON /CCG/ ID(27,11)
C
      ABSIS = ' STRAIN IN PERCENT '
C
      READ(5,*) NST
      NC = 2*NST
      DO 1 L = 1,NC
      READ(5,2001) NV(L), (ID(L,I), I=1,11)
      M = NV(L)
      READ(5,1002) (S(L,I), I = 1,M)
      READ(5,1002) (Y(L,I), I = 1,M)
    1 CONTINUE
C---------------------------------------------
C     ECHO INPUT DYNAMIC PROPERTY CURVES
C---------------------------------------------
      READ (5, 1007) NECHO,(ECHO(I), I=1,NECHO)
C      DO 10, I=1, NC, 2
      DO 10, K=1,NECHO
C      MTYPE=(I+1)/2
      MTYPE=ECHO(K)
      I=ECHO(K)*2-1
      WRITE(6,1003) MTYPE
      WRITE(6,1004) I, (ID(I,J), J=1,10), I+1, (ID(I+1,J),J=1,10)
      WRITE(6,1005) I,I+1
      M=MAX0(NV(I),NV(I+1))
   10 WRITE(6,1006) (S(I,J),Y(I,J),S(I+1,J),Y(I+1,J), J=1,M)
      CALL CURVEG( NC, NV,  2, AS, BS, 10, TSTEP, NT, T, V, S, Y, NSTEP)
C--------------------------
 1000 FORMAT(3I5,F10.0)
 1001 FORMAT(/I5,F10.2, 11A6)
 1002 FORMAT(8F10.3)
 1003 FORMAT (/' **********************'/
     +       '  MATERIAL TYPE NO.',I2,/
     +       ' **********************')
 1004 FORMAT(/2(1X,'CURVE NO. ',I2,': ', 10a6/)/)
 1005 FORMAT(  '      CURVE NO.',I2,'              CURVE NO.',I2,'    '/
     1         '  ===================       =================='/
     2         '    STRAIN    G/Gmax         STRAIN    DAMPING'/
     3         '  --------   -------        --------  --------')
 1006 FORMAT(1X,F9.4,4X,F6.3,5X,1X,F9.4,4X,F6.2,5X)
 1007 FORMAT(16I5)
 2001 FORMAT(I5, 11A6)
 2002 FORMAT( 12F10.4)
 3000 FORMAT(53H  MODULUS AND DAMPING VALUES ARE SCALED FOR PLOTTING   )
 3001 FORMAT(55H    CURVES FOR RELATION STRAIN VERSUS SHEAR MODULUS AND
     1 8H DAMPING /)
      RETURN
      END

