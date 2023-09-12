C==========================================================================
      SUBROUTINE EXPAN3D(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
      DIMENSION ID(6,*)
      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,DCEXP,DISTAL,DIS1AL,RE,SVEXC,VALST
     &,PERIO1,TOL,NUMST,MESH,NUMALV,NBSTAC,NPRVEL,MOVEW,IAKIRA

      TT=T-PERIO1/4.D0
c      TT=T

C==========================================================================
C Only for Akira's animation, Oct. 13, 2006
c      call EXPAN3Dballon(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,
c     &NSTAC,ZADVRE,AMI,KKORAK)
	return
C==========================================================================

      IF (IAKIRA.EQ.0) RETURN

      PI=4.D0*DATAN(1.D0)

      MOVX=0
      MOVY=0
      IF (MOVEW.EQ.1) MOVX=1
      IF (MOVEW.EQ.2) MOVY=1
      IF (MOVEW.EQ.0) THEN
        MOVX=1
        MOVY=1
      ENDIF
 
      P=PERIO1
      AN=2.D0*PI/P
      C=SVEXC
      FI=(1.D0+C)**(1.D0/3.D0)
      AK=(FI-1.D0)/(FI+1.D0)


C EXPANSION CENTER
      CX=0.D0
	CY=0.D0
      CZ=-0.0215D0
c      CX=-0.0400D0
c	CY=0.D0
c      CZ=-0.0430D0
 
C=======================================================================

      DO 10 I=1,NPT
        VX=0.D0
        VY=0.D0
        VZ=0.D0
        X=CORD(1,I)
        Y=CORD(2,I)
	  Z=CORD(3,I)
        X0=CCORD(1,I)
        Y0=CCORD(2,I)
           CCORD(2,I)=CY+(Y-CY)*(1.D0+AK*DSIN(AN*TT))
           VY=(Y-CY)*AK*AN*DCOS(AN*TT)
C Ovde izvrsena izmena za granican uslov skracenja duzine cevi
C odnosno zadavanje brzine na kracoj duzini prave cevi
C           VY=(3.D0*ALEN-2.D0*Y)*AK*AN*DCOS(AN*Time)

           CCORD(1,I)=CX+(X-CX)*(1.D0+AK*DSIN(AN*TT))
           VX=(X-CX)*AK*AN*DCOS(AN*TT)
           CCORD(3,I)=CZ+(Z-CZ)*(1.D0+AK*DSIN(AN*TT))
           VZ=(Z-CZ)*AK*AN*DCOS(AN*TT)


c        VMESH(1,I)=(CCORD(1,I)-X0)/TIME
c        VMESH(2,I)=(CCORD(2,I)-Y0)/TIME
C        VMESH(2,I)=(CCORD(2,I)-Y0)/TIME
       IF(ID(1,I).EQ.0.AND.ID(2,I).EQ.0.AND.ID(3,I).EQ.0) THEN
c	if (cord(2,i).gt.0.d0) goto 10
c        IF (NPRVEL.EQ.1.AND.I.LE.MESH) THEN
c          GNODE(2,2,I)=ZADVRE(I)*(1.D0+AK*DSIN(AN*Time))
C          GOTO 10
C        ENDIF
         GNODE(2,1,I)=VX
         GNODE(2,2,I)=VY
         GNODE(2,3,I)=VZ
       ENDIF
  10   CONTINUE
 

      END
C==========================================================================
