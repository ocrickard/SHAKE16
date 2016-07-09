C $NOFLOATCALLS
C $NODEBUG
C ...................................................................
      program SHAKE91
      character*32 FIN,FOUT,PUNCH
      common X(25620)
      common /TIME/ T(9)
      common /WGK/  WW, GT, SKO
C ...................................................................
      write(*,100)
  100 format(2X,'*****************************************************'/
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
      write(*,200)
  200 format(4X,'Name of Input File =')
      read(*,10) FIN
C
      write(*,300)
  300 format(4X,'Name of Output File #1 (input, peak values .. etc) =')
      read(*,10) FOUT
C
      write(*,400)
  400 format(4X,'Name of Output File #2 (time histories .. etc) =')
      read(*,10) PUNCH
C
   10 format(A32)
C
      open(5,FILE=FIN,STATUS='OLD')
      open(6,FILE=FOUT,STATUS='NEW')
      open(7,FILE=PUNCH,STATUS='NEW')
C
      write(6,100)
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
      if (SKO .LT. .000001)   SKO = .45
      write(6,2000) MAMAX, NTOT
C
      callSHAKIT(X(1), X(NAX), X(NAA), X(NS), X(NINV))
C
      stop
C ****************************************************
 1000 format(I5, F10.0)
 2000 format( 45H  MAX. NUMBER OF TERMS IN FOURIER TRANSFORM =   I10/
     1          45H  NECESSARY LENGTH OF BLANK COMMON X        =   I10)
      end
C********************************************************************
      subroutine EARTHQ(X,AX,S,INV)
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
      character*6 TITLE
      character*30 FINPEQ
      character*80 HEAD
      character*12 FMAT
      complex X, AX
      dimension XR(8),X(300),AX(3,270),S(70),INV(70)
      common /EQ/ MFOLD,MA2,TITLE(5),DT, MA , MMA, DF, MX
      common/FRCUT/ NCUT,NZERO
      common /JISCK/ JIS,FINPEQ
C
      PI2 = 6.283185307
      read(5,1001) NV, MA, DT, FINPEQ, FMAT
      read(5,1004) XF,XMAX,FMAX,NHEAD, NPL
      if (FMAX .LT. .001) FMAX = 100000.
      if (FMAT .EQ. '            ') FMAT = '(8F9.6,I7)'
      MA2=2
    2 if (MA2.GE.MA) go to 3
      MA2=MA2*2
      go to 2
    3 MA=MA2
C  ....................................................................
      write(6,2012) FINPEQ, NV, MA, NHEAD, NPL, DT, FMAT
      write (*,2026) FINPEQ,FMAT
      open(8,FILE=FINPEQ,STATUS='OLD')
      write (6,2021)
      do 4 I=1,NHEAD
      read(8,2022) HEAD
    4 write(6,2022) HEAD
C      write (6,2023)
      MMA = NV + NV / 10
      if (MMA .GT. MA) MMA=MA
      MA2 = MA + 2
      MFOLD = MA2/2
      MFOLD = MFOLD + 1
      DF = 1./(MA*DT)
      FMA = FLOAT(MA)
      MX = (ALOG10(FMA)/ALOG10(2.))- 1
    1 NMX = 2**(MX+1)
      if (MA .LE. NMX) go to 11
      MX = MX + 1
      go to 1
   11 NCARDS = (NV-1)/NPL + 1
      JL=NPL*NCARDS-NV
      NV = NV + 1
      N = 0
      LC = 0
      write(6,2024)
      do 31 I = 1,NCARDS
      LC = LC + 1
      read(8,FMAT) (XR(J), J=1,NPL)
      if (I .NE. NCARDS) go to 6
      if (JL .EQ. 0) go to 6
      JL=NPL+1-JL
      do 5  J=JL,NPL
    5 XR(J)=0.
C---------------------------------------------------------------
C  ONLY PRINT OUT A FEW OF THE FIRST AND THE LAST LINES of Input Motion
C---------------------------------------------------------------
    6 ICHECK = NCARDS - I
      if (I .LE. 5 .OR. ICHECK .LT. 5) write(6,2008) I,(XR(J), J=1,NPL)
        if (I .EQ. 10) write (6,2009)
 2009 format(3X,'........ INPUT MOTION READ NOT ECHOED...........')
C     endIF
C
C   FIND MAX. INPUT ACC. (XMAX)
C
  311 do 31 J = 1,NPL,2
      N = N + 1
      X(N) = CMPLX(XR(J),XR(J+1))
   31 continue
      close (8)
      N = N + 1
      do 32 I = N,MFOLD
   32 X(I) = 0.
      callXMX(X,MA,XM,NXMAX)
      if (XMAX.LT. .000001) go to 300
      XF = XMAX/XM
  300 do 30 I = 1,N
   30 X(I) = X(I)*XF
      XMAX = XM*XF
      TMAX = FLOAT(NXMAX-1)*DT
      write(6,2014) XM,TMAX,XF,XMAX
C
      callRFFT(X,MX,INV,S,IFERR,1)
C
C   REMOVE FREQUENCIES ABOVE FMAX AND FIND MAX. ACC. OF NEW MOTION
C
      FREQ = 0.
      SXX = 0.
      SFX = 0.
      NCUT=0
      do 33 I = 1,MFOLD
      if (FREQ .LE. FMAX) go to 34
      NCUT=NCUT+1
      X(I)=0.0
   34 continue
      XA = CABS(X(I))
      SXX= SXX + XA*XA
      SFX = SFX + FREQ*XA*XA
      AX(1,I) = X(I)
      FREQ = FREQ + DF
   33 continue
      SFX = SFX/SXX
      NCUT=MFOLD-NCUT
      NZERO=NCUT+1
      write(6,2005) SFX
      if (FMAX.GT.FREQ) return
      callRFSN(X,MX,INV,S,IFERR,-2)
      callXMX(X,MA,XM,NXMAX)
      do 72 I = 1,MFOLD
   72 X(I) = AX(1,I)
      write(6,2001) XM,FMAX
C
 1001 format(2I5,  F10.3, A30, A12)
 1002 format(8F9.5,I7)
 1003 format(8F10.0)
 1004 format(3F10.0,2I5)
 2001 format(21H  MAX ACCELERATION =  F10.5, 22H FOR FREQUENCIES REMOV
     19HED ABOVE F10.2, 7H C/SEC.)
 2003 format(17H   ACC. CARD NO. I4,16H OUT OF SEQUENCE  )
 2005 format(25H MEAN SQUARE FREQUENCY =  F10.2, 7H C/SEC. )
C2008 format(2X, I5, 5X, 8F15.6)
 2008 format(1X,I5,1X, 8F9.6)
 2012 format(/1X, ' FILE NAME FOR INPUT MOTION = ', A30,/
     +         1X, '   NO. OF INPUT ACC. POINTS = ',I5,/
     +         1X, '  NO. OF POINTS USED IN FFT = ',I5/
     +         1X, '       NO. OF HEADING LINES = ',I5/
     +         1X, '     NO. OF POINTS PER LINE = ',I5/
     +         1X, ' TIME STEP FOR INPUT MOTION = ',F6.4/
     +         1X, ' format FOR OF TIME HISTORY = ', A12, /)
 2014 format(/23H MAXIMUM ACCELERATION =  F9.5/
     1        23H AT TIME              =  F6.2, 4H SEC/
     1 44H THE VALUES WILL BE MULTIPLIED BY A FACTOR = F7.3/
     3 44H TO GIVE NEW MAXIMUM ACCELERATION          = F9.5 )
 2021 format (/1X,'***** H E A D E R   ')
 2022 format (A80)
 2023 format (1X,'*************************************************'
     1       ,'*****************')
 2024 format (' **  FIRST & LAST 5 LINES OF INPUT MOTION *****'/)
 2025 format (' *****************************************************'/)
 2026 format(/,1X, ' READING INPUT MOTION FROM ----> ',A30/
     +         1X, ' format OF INPUT MOTION USED --> ',A12)
      return
      end
