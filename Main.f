C $NOFLOATCALLS
C $NODEBUG
C ...................................................................
      PROGRAM SHAKE91
      CHARACTER*32 FIN,FOUT,PUNCH
      COMMON X(25620)
      COMMON /TIME/ T(9)
      COMMON /WGK/  WW, GT, SKO
C ...................................................................
      WRITE(*,100)
  100 FORMAT(2X,'*****************************************************'/
     + 2X,'* SHAKE  --   A COMPUTER PROGRAM FOR EARTHQUAKE RESPONSE  *'/
     + 2X,'*             ANALYSIS OF HORIZONTALLY LAYERED SITES      *'/
     + 2X,'*             by: Per B. Schnabel & John Lysmer -- 1970   *'/
     + 2X,'* ------------------------------------------------------- *'/
     + 2X,'* shake85     IBM-PC version of SHAKE                     *'/
     + 2X,'*             by: S.S. (Willie) Lai, January 1985         *'/
     + 2X,'* ------------------------------------------------------- *'/
     + 2X,'* shake88   : New modulus reduction curves for clays added*'/
     + 2X,'*             using results from Sun et al (1988)         *'/
     + 2X,'*             by: J. I. Sun & Ramin Golesorkhi            *'/
     + 2X,'*             February 26, 1988                           *'/
     + 2X,'* ------------------------------------------------------- *'/
     + 2X,'* SHAKE90/91: Adjust last iteration; Input now is either  *'/
     + 2X,'*             Gmax or max Vs; up to 13 material types can *'/
     + 2X,'*             be specified by user; up to 50 Layers can   *'/
     + 2X,'*             be specified; object motion can be read in  *'/
     + 2X,'*             from a separate file and can have user      *'/
     + 2X,'*             specified format; Different periods for     *'/
     + 2X,'*             response spectral calculations; options     *'/
     + 2X,'*             are renumbered; and general cleanup         *'/
     + 2X,'*             by: J. I. Sun, I. M. Idriss & P. Dirrim     *'/
     + 2X,'*             June 1990 - February 1991                   *'/
     + 2X,'* ------------------------------------------------------- *'/
     + 2X,'* SHAKE91   : General cleanup and finalization of input/  *'/
     + 2X,'*             output format ... etc                       *'/
     + 2X,'*             by: I. M. Idriss                            *'/
     + 2X,'*             December 1991                               *'/
     + 2X,'***********************************************************')
C
      WRITE(*,200)
  200 FORMAT(4X,'Name of Input File =')
      READ(*,10) FIN
C
      WRITE(*,300)
  300 FORMAT(4X,'Name of Output File #1 (input, peak values .. etc) =')
      READ(*,10) FOUT
C
      WRITE(*,400)
  400 FORMAT(4X,'Name of Output File #2 (time histories .. etc) =')
      READ(*,10) PUNCH
C
   10 FORMAT(A32)
C
      OPEN(5,FILE=FIN,STATUS='OLD')
      OPEN(6,FILE=FOUT,STATUS='NEW')
      OPEN(7,FILE=PUNCH,STATUS='NEW')
C
      WRITE(6,100)
      WW = .0624
      GT = 32.2
      MAMAX=4096
C ......................................................................
C
      NAX = MAMAX + 5
      NAA = NAX + 3*(MAMAX + 4)
      NS = NAA + 2*MAMAX
      NINV = NS + NAX/8 + 1
      NTOT = NINV + NAX/8 + 1
      IF (SKO .LT. .000001)   SKO = .45
      WRITE(6,2000) MAMAX, NTOT
C
      CALL SHAKIT(X(1), X(NAX), X(NAA), X(NS), X(NINV))
C
      STOP
C ****************************************************
 1000 FORMAT(I5, F10.0)
 2000 FORMAT( 45H  MAX. NUMBER OF TERMS IN FOURIER TRANSFORM =   I10/
     1          45H  NECESSARY LENGTH OF BLANK COMMON X        =   I10)
      END
C********************************************************************
      SUBROUTINE EARTHQ(X,AX,S,INV)
C***********************************************************************
C
C   THIS ROUTINE  READS THE MOTION IN THE TIME DOMAIN, ADDS TRAILING
C   ZEROS, SCALES THE VALUES, FIND MAXIMUM VALUE AND VARIOUS PARAMETERS
C   AND TRANSFER THE MOTION INTO THE FREQUENCY DOMAIN.
C
C   CODED BY PBS SEPT. 1970
C
C        X       = INPUT MOTION
C        AX      = TEMPORARY STORAGE OF X
C        TITLE   = IDENTIFICATION FOR MOTION
C        DT      = TIME STEP BETWEEN VALUES IN TIME DOMAIN
C        NV      = NUMBER OF ACC. VALUES TO BE READ
C        MA      = LENGTH OF MOTION INCLUDING TRAILING ZEROS
C        MMA     = LENGTH OF SIGNIFICANT PART OF MOTION
C        XF      = MULTIPLICATION FACTOR FOR ACCELERATION VALUES
C        DF      = FREQUENCY STEPS IN FREQ. DOMAIN
C
C***********************************************************************
C
      CHARACTER*6 TITLE
      CHARACTER*30 FINPEQ
      CHARACTER*80 HEAD
      CHARACTER*12 FMAT
      COMPLEX X, AX
      DIMENSION XR(8),X(300),AX(3,270),S(70),INV(70)
      COMMON /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF, MX
      COMMON/FRCUT/ NCUT,NZERO
      COMMON /JISCK/ JIS,FINPEQ
C
      PI2 = 6.283185307
      READ(5,1001) NV, MA, DT, FINPEQ, FMAT
      READ(5,1004) XF,XMAX,FMAX,NHEAD, NPL
      IF (FMAX.LT. .001) FMAX = 100000.
      IF (FMAT .EQ. '            ') FMAT = '(8F9.6,I7)'
      MA2=2
    2 IF(MA2.GE.MA) GO TO 3
      MA2=MA2*2
      GO TO 2
    3 MA=MA2
C  ....................................................................
      WRITE(6,2012) FINPEQ, NV, MA, NHEAD, NPL, DT, FMAT
      WRITE (*,2026) FINPEQ,FMAT
      OPEN(8,FILE=FINPEQ,STATUS='OLD')
      WRITE (6,2021)
      DO 4 I=1,NHEAD
      READ(8,2022) HEAD
    4 WRITE(6,2022) HEAD
C      WRITE (6,2023)
      MMA = NV + NV/10
      IF (MMA.GT.MA) MMA=MA
      MA2 = MA + 2
      MFOLD = MA2/2
      MFOLD = MFOLD + 1
      DF = 1./(MA*DT)
      FMA = FLOAT(MA)
      MX = (ALOG10(FMA)/ALOG10(2.))- 1
    1 NMX = 2**(MX+1)
      IF (MA  .LE. NMX) GO TO 11
      MX = MX + 1
      GO TO 1
   11 NCARDS = (NV-1)/NPL + 1
      JL=NPL*NCARDS-NV
      NV = NV + 1
      N = 0
      LC = 0
      WRITE(6,2024)
      DO 31 I = 1,NCARDS
      LC = LC + 1
      READ(8,FMAT) (XR(J), J=1,NPL)
      IF(I.NE.NCARDS) GO TO 6
      IF(JL.EQ.0) GO TO 6
      JL=NPL+1-JL
      DO 5  J=JL,NPL
    5 XR(J)=0.
C---------------------------------------------------------------
C  ONLY PRINT OUT A FEW OF THE FIRST AND THE LAST LINES of Input Motion
C---------------------------------------------------------------
    6 ICHECK = NCARDS - I
      IF (I .LE. 5 .OR. ICHECK .LT. 5) WRITE(6,2008) I,(XR(J), J=1,NPL)
        IF (I .EQ. 10) WRITE (6,2009)
 2009 FORMAT(3X,'........ INPUT MOTION READ NOT ECHOED...........')
C     ENDIF
C
C   FIND MAX. INPUT ACC. (XMAX)
C
  311 DO 31 J = 1,NPL,2
      N = N + 1
      X(N) = CMPLX(XR(J),XR(J+1))
   31 CONTINUE
      CLOSE (8)
      N = N + 1
      DO 32 I = N,MFOLD
   32 X(I) = 0.
      CALL XMX(X,MA,XM,NXMAX)
      IF (XMAX.LT. .000001) GO TO 300
      XF = XMAX/XM
  300 DO 30 I = 1,N
   30 X(I) = X(I)*XF
      XMAX = XM*XF
      TMAX = FLOAT(NXMAX-1)*DT
      WRITE(6,2014) XM,TMAX,XF,XMAX
C
      CALL RFFT(X,MX,INV,S,IFERR,1)
C
C   REMOVE FREQUENCIES ABOVE FMAX AND FIND MAX. ACC. OF NEW MOTION
C
      FREQ = 0.
      SXX = 0.
      SFX = 0.
      NCUT=0
      DO 33 I = 1,MFOLD
      IF(FREQ.LE.FMAX) GO TO 34
      NCUT=NCUT+1
      X(I)=0.0
   34 CONTINUE
      XA = CABS(X(I))
      SXX= SXX + XA*XA
      SFX = SFX + FREQ*XA*XA
      AX(1,I) = X(I)
      FREQ = FREQ + DF
   33 CONTINUE
      SFX = SFX/SXX
      NCUT=MFOLD-NCUT
      NZERO=NCUT+1
      WRITE(6,2005) SFX
      IF (FMAX.GT.FREQ) RETURN
      CALL RFSN(X,MX,INV,S,IFERR,-2)
      CALL XMX(X,MA,XM,NXMAX)
      DO 72 I = 1,MFOLD
   72 X(I) = AX(1,I)
      WRITE(6,2001) XM,FMAX
C
 1001 FORMAT(2I5,  F10.3, A30, A12)
 1002 FORMAT(8F9.5,I7)
 1003 FORMAT(8F10.0)
 1004 FORMAT(3F10.0,2I5)
 2001 FORMAT(21H  MAX ACCELERATION =  F10.5, 22H FOR FREQUENCIES REMOV
     19HED ABOVE F10.2, 7H C/SEC.)
 2003 FORMAT(17H   ACC. CARD NO. I4,16H OUT OF SEQUENCE  )
 2005 FORMAT(25H MEAN SQUARE FREQUENCY =  F10.2, 7H C/SEC. )
C2008 FORMAT(2X, I5, 5X, 8F15.6)
 2008 FORMAT(1X,I5,1X, 8F9.6)
 2012 FORMAT(/1X, ' FILE NAME FOR INPUT MOTION = ', A30,/
     +         1X, '   NO. OF INPUT ACC. POINTS = ',I5,/
     +         1X, '  NO. OF POINTS USED IN FFT = ',I5/
     +         1X, '       NO. OF HEADING LINES = ',I5/
     +         1X, '     NO. OF POINTS PER LINE = ',I5/
     +         1X, ' TIME STEP FOR INPUT MOTION = ',F6.4/
     +         1X, ' FORMAT FOR OF TIME HISTORY = ', A12, /)
 2014 FORMAT(/23H MAXIMUM ACCELERATION =  F9.5/
     1        23H AT TIME              =  F6.2, 4H SEC/
     1 44H THE VALUES WILL BE MULTIPLIED BY A FACTOR = F7.3/
     3 44H TO GIVE NEW MAXIMUM ACCELERATION          = F9.5 )
 2021 FORMAT (/1X,'***** H E A D E R   ')
 2022 FORMAT (A80)
 2023 FORMAT (1X,'*************************************************'
     1       ,'*****************')
 2024 FORMAT (' **  FIRST & LAST 5 LINES OF INPUT MOTION *****'/)
 2025 FORMAT (' *****************************************************'/)
 2026 FORMAT(/,1X, ' READING INPUT MOTION FROM ----> ',A30/
     +         1X, ' FORMAT OF INPUT MOTION USED --> ',A12)
      RETURN
      END
