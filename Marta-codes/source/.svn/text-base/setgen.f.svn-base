      SUBROUTINE SETGEN(F,NDIM,NPOIN,NPRIN,NTREAT)
C
C  AUTHOR      : J. VERMASEREN
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL F
      DIMENSION X(11),N(11)
      common/vgmaxi/MDUM,MBIN,FFMAX,FMAX(80000),NM(80000)
      data mbin/3/
      real ran2
C
c     WRITE(6,*) ' =======> AFTER CALL     NPOIN =',NPOIN
C
      FFMAX=0.
      SUM=0.
      SUM2=0.
      SUM2P=0.
      MAX=MBIN**NDIM
      IF(NPRIN.GE.2)WRITE(NOUTP,200)MBIN,MAX,NPOIN
      DO 5 J=1,MAX
         NM(J)=0
         FMAX(J)=0.
5     CONTINUE
      DO 1 J=1,MAX
         JJ=J-1
         DO 2 K=1,NDIM
            JJJ=JJ/MBIN
            N(K)=JJ-JJJ*MBIN
            JJ=JJJ
2        CONTINUE
         FSUM=0.
         FSUM2=0.
         DO 3 M=1,NPOIN
            DO 4 K=1,NDIM
               X(K)=(ran2(dum)+N(K))/MBIN
4           CONTINUE
            IF(NTREAT.GT.0)Z=TREAT(F,X,NDIM)
            IF(NTREAT.LE.0)Z=F(X)
c            if((z.lt.0).or.(z.gt.0))print *,'z=',z
            IF(Z.GT.FMAX(J))FMAX(J)=Z
            FSUM=FSUM+Z
            FSUM2=FSUM2+Z*Z
3        CONTINUE
c         WRITE(6,*) ' =======> BEFOR DEVISION NPOIN =',NPOIN
         AV=FSUM/NPOIN
         AV2=FSUM2/NPOIN
         SIG2=AV2-AV*AV
         SIG=SQRT(SIG2)
         SUM=SUM+AV
         SUM2=SUM2+AV2
         SUM2P=SUM2P+SIG2
         IF(FMAX(J).GT.FFMAX)FFMAX=FMAX(J)
         EFF=10000.
         IF(FMAX(J).NE.0)EFF=FMAX(J)/AV
         IF(NPRIN.GE.3)WRITE(NOUTP,100)J,AV,SIG,FMAX(J),EFF,
     +                                 (N(KJ),KJ=1,NDIM)
1     CONTINUE
      SUM=SUM/MAX
      SUM2=SUM2/MAX
      SUM2P=SUM2P/MAX
      SIG=SQRT(SUM2-SUM*SUM)
      SIGP=SQRT(SUM2P)
      EFF1=0.
      DO 6 J=1,MAX
         EFF1=EFF1+FMAX(J)
6     CONTINUE
      EFF1=EFF1/(MAX*SUM)
      EFF2=FFMAX/SUM
      IF(NPRIN.GE.1)WRITE(NOUTP,101)SUM,SIG,SIGP,FFMAX,EFF1,EFF2
C
100   FORMAT(I6,3X,G13.6,G12.4,G13.6,F8.2,3X,10I1)
101   FORMAT(' Average function value =',G14.6/
     +       ' Overall standard dev.  =',G14.4/
     +       ' Average standard dev.  =',G14.4/
     +       ' Maximum function value =',G14.6/
     +       ' Average inefficiency   =',G14.3/
     +       ' Overall inefficiency   =',G14.3/)
200   FORMAT(' Subroutine SETGEN uses a',I3,'**NDIM division'/
     + ' this results in ',I7,' cubes'/
     + ' The program put ',I5,' points in each cube to find',
     + ' starting values for the maxima'//)
C
      RETURN
      END
