*-- Author :
      SUBROUTINE GMUGNA
C
C  MODIFICATION FROM VERMASERENS GENERA
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision Y
      EXTERNAL sigma
      integer dum
      real ran2
      common/vggen/iend,ndim,npoin,nprin,ntreat,itmx,ncvg

      data dum/-1/

*KEEP,COMGNA.
      INTEGER NGNA
      COMMON /COMMUP/  NGNA
      SAVE /commup/
*KEND.
C
      COMMON/VGMAXI/MDUM,MBIN,FFMAX,FMAX(80000),NM(80000)
c      SAVE /vgmaxi/
      DIMENSION X(11),N(11)
      SAVE x,n
C
      SAVE WEIGHT,CORREC,J
C
      SAVE CORRE2,FMDIFF,FMOLD,FMAX2
      DATA j/0/

      integer ncorrec
      data ncorrec/0/
      save ncorrec
c
c
      NGNA=NGNA+1
C
      AMI=1.0D0/DBLE(MBIN)
      MAX=MBIN**NDIM
C
C CORRECTION CYCLES ARE STARTED <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF (J .NE. 0) THEN
         ncorrec = ncorrec + 1
c         print *,'num of times the correction was applied :',ncorrec
!          print *,'gmugna: correction',j,correc,corre2
 4       CONTINUE
         IF (CORREC .LT. 1.0)THEN
            IF(ran2(dum) .GE. CORREC) GOTO 7
            CORREC=-1.0
         ELSE
            CORREC=CORREC-1.0
         ENDIF
C
C SEL X VALUES IN VEGAS BIN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         DO 6 K=1,NDIM
            X(K)=(ran2(dum)+N(K))*AMI
!            if(x(k).lt.0.or.x(k).gt.1) then
!               print *,'correction, X',x(k),k
!            endif
 6       CONTINUE
C
C COMPUTE WEIGHT FOR X VALUES <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         IF (NTREAT .GT. 0) WEIGHT=TREAT(sigma,X,NDIM)
         IF (NTREAT .LE. 0) WEIGHT=sigma(X)
C
C PARAMETER FOR CORRECTION OF CORRECTION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         IF(WEIGHT .GT. FMAX(J)) THEN
            IF (WEIGHT .GT. FMAX2) FMAX2=WEIGHT
            CORRE2=CORRE2-1.0
            CORREC=CORREC+1.0
         ENDIF
C
C ACCEPT EVENT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
         IF (WEIGHT .GE. FMDIFF*ran2(dum)+FMOLD) THEN
!            print *,'gmugna: done',weight,' x',(x(ipr),ipr=1,NDIM)
            RETURN
         ENDIF
         GOTO 4
C
C CORRECTION IF TOO BIG WEIGHT IS FOUND WHILE CORRECTION <<<<<<<<<<<<
 7       CONTINUE
         IF (FMAX2 .GT. FMAX(J)) THEN
            FMOLD=FMAX(J)
            FMAX(J)=FMAX2
            FMDIFF=FMAX2-FMOLD
!            print *,'using CORRE2',corre2
            IF(FMAX2 .LE. FFMAX) THEN
               CORREC=(NM(J)-1.0)*FMDIFF/FFMAX-CORRE2
            ELSE
               FFMAX=FMAX2
               CORREC=(NM(J)-1.0)*FMDIFF/FFMAX*FMAX2/FFMAX-CORRE2
            ENDIF
            CORRE2=0.0
            FMAX2=0.0
            GOTO 4
         ENDIF
      ENDIF
C
C  NORMAL GENERATION CYCLE STARTS HERE !!!!!!!       *******************
C
c      print *,'gmugna: generation'

C  SEL A VEGAS BIN AND REJECT IF FMAX IS TOO LITTLE <<<<<<<<<<<<<<<<<<<<

 1    CONTINUE
      J=ran2(dum)*MAX+1.
      Y=ran2(dum)*FFMAX
      NM(J)=NM(J)+1
c      print *,'gmugna: y=',y,', j=',j,', fmax(j)=',fmax(j),
c     +', ffmax=',ffmax,'random:',ran2(dum),'dum:',dum
      IF(Y.GT.FMAX(J)) GOTO 1
C
C  SEL X VALUES IN THIS VEGAS BIN <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      JJ=J-1
      DO 2 K=1,NDIM
c         print *,k,'->',jj,mbin
         JJJ=JJ/MBIN
         N(K)=JJ-JJJ*MBIN
         X(K)=(ran2(dum)+N(K))*AMI
!         if(x(k).lt.0.or.x(k).gt.1) then
!            print *,'correction, X',x(k),k
!            endif
         JJ=JJJ
 2    CONTINUE
C
C  GET WEIGHT FOR SEL X VALUES <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF(NTREAT.GT.0) THEN
         WEIGHT=TREAT(sigma,X,NDIM)
      ELSE
         WEIGHT=sigma(X)
      ENDIF

C  REJECT IF WEIGHT IS TOO LOW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF(Y .GT. WEIGHT) THEN
c         print *,'weight is too low'
         GO TO 1
      ENDIF
      IF(WEIGHT .LE. FMAX(J)) THEN
         J=0
C
C  INIT CORRECTION CYCLES IF WEIGHT IS HIGHER THEN FMAX OR FFMAX <<<<<<
      ELSE IF(WEIGHT .LE. FFMAX) THEN
         FMOLD=FMAX(J)
         FMAX(J)=WEIGHT
         FMDIFF=WEIGHT-FMOLD
         CORREC=(NM(J)-1.0)*FMDIFF/FFMAX-1.0
c         print *,'need correction of FFMAX',weight,ffmax,correc
      ELSE
         FMOLD=FMAX(J)
         FMAX(J)=WEIGHT
         FMDIFF=WEIGHT-FMOLD
         FFMAX=WEIGHT
         CORREC=(NM(J)-1.0)*FMDIFF/FFMAX*WEIGHT/FFMAX-1.0
c         print *,'need correction of FMAX',j,weight,ffmax,correc
      ENDIF
C
C RETURN WITH AN ACCEPTED EVENT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c      print *,'gmugna: done',weight,' x',(x(ipr),ipr=1,NDIM)

      RETURN
      END
