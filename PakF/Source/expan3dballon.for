C==========================================================================
      SUBROUTINE EXPAN3Dballon(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,
     &NSTAC,ZADVRE,AMI,KKORAK)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
      DIMENSION ID(6,*)
      COMMON /EXPTUB/ RD,RA,ALEN,GAMA,DCEXP,DISTAL,DIS1AL,RE,SVEXC,VALST
     &,PERIO1,TOL,NUMST,MESH,NUMALV,NBSTAC,NPRVEL,MOVEW,IAKIRA

      TT=T-PERIO1/4.D0

      IF (IAKIRA.EQ.0) RETURN

      PI=4.D0*DATAN(1.D0)


      P=PERIO1
      AN=2.D0*PI/P
      C=SVEXC
      FI=(1.D0+C)**(1.D0/3.D0)
      AK=(FI-1.D0)/(FI+1.D0)


      ak=ak*4.d0
C==========================================================================
C     Ballon above X axis
C==========================================================================
C EXPANSION CENTER
	nc=6415+5
	cx=cord(1,nc)
	cy=cord(2,nc)
	cz=cord(3,nc)
      do 10 ii=2,10      
      n1=6415+(ii-1)*132
	n2=n1+65
      do 10 i=n1,n2
        VX=0.D0
        VY=0.D0
        VZ=0.D0
        X=CORD(1,I)
        Y=CORD(2,I)
	  Z=CORD(3,I)
           CCORD(2,I)=CY+(Y-CY)*(1.D0+AK*DSIN(AN*TT))
           VY=(Y-CY)*AK*AN*DCOS(AN*TT)
           CCORD(1,I)=CX+(X-CX)*(1.D0+AK*DSIN(AN*TT))
           VX=(X-CX)*AK*AN*DCOS(AN*TT)
           CCORD(3,I)=CZ+(Z-CZ)*(1.D0+AK*DSIN(AN*TT))
           VZ=(Z-CZ)*AK*AN*DCOS(AN*TT)
       IF(ID(1,I).EQ.0.AND.ID(2,I).EQ.0.AND.ID(3,I).EQ.0) THEN
         GNODE(2,1,I)=VX
         GNODE(2,2,I)=VY
         GNODE(2,3,I)=VZ
       ENDIF
  10   CONTINUE
C==========================================================================
C==========================================================================
C     Ballon bellow X axis
C==========================================================================
C EXPANSION CENTER
	nc=6481+5
	cx=cord(1,nc)
	cy=cord(2,nc)
	cz=cord(3,nc)
      do 20 ii=2,10      
      n1=6481+(ii-1)*132
	n2=n1+65
      do 20 i=n1,n2
        VX=0.D0
        VY=0.D0
        VZ=0.D0
        X=CORD(1,I)
        Y=CORD(2,I)
	  Z=CORD(3,I)
           CCORD(2,I)=CY+(Y-CY)*(1.D0+AK*DSIN(AN*TT))
           VY=(Y-CY)*AK*AN*DCOS(AN*TT)
           CCORD(1,I)=CX+(X-CX)*(1.D0+AK*DSIN(AN*TT))
           VX=(X-CX)*AK*AN*DCOS(AN*TT)
           CCORD(3,I)=CZ+(Z-CZ)*(1.D0+AK*DSIN(AN*TT))
           VZ=(Z-CZ)*AK*AN*DCOS(AN*TT)
       IF(ID(1,I).EQ.0.AND.ID(2,I).EQ.0.AND.ID(3,I).EQ.0) THEN
         GNODE(2,1,I)=VX
         GNODE(2,2,I)=VY
         GNODE(2,3,I)=VZ
       ENDIF
  20   CONTINUE
C==========================================================================
 

      END
C==========================================================================
C==========================================================================
      SUBROUTINE TwoBallons(CORD,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)


C==========================================================================
C     Ballon above X axis
      do 100 i=1,10      
      n1=6415+(i-1)*132
	n2=n1+65
	p=-(i-1)*(i-10)/5.d0+1.d0
	nc=n1+5
      do n=n1,n2
	 cord(1,n)=(1.d0-p)*cord(1,nc)+p*cord(1,n)
	 cord(2,n)=(1.d0-p)*cord(2,nc)+p*cord(2,n)
	 cord(3,n)=(1.d0-p)*cord(3,nc)+p*cord(3,n)
	enddo
100   continue
      do 30 i=2,10      
      n1=6415+(i-1)*132
	n2=n1+65
	p=1.5d0
	k=0
      do n=n1,n2
	 nc=6415+k
	 k=k+1
	 cord(1,n)=(1.d0-p)*cord(1,nc)+p*cord(1,n)
	 cord(2,n)=(1.d0-p)*cord(2,nc)+p*cord(2,n)
	 cord(3,n)=(1.d0-p)*cord(3,nc)+p*cord(3,n)
	enddo
 30   continue

C==========================================================================
C==========================================================================
C     Ballon bellow X axis
      do 200 i=1,10      
      n1=6481+(i-1)*132
	n2=n1+65
	p=-(i-1)*(i-10)/10.d0+1.d0
	nc=n1+5
      do n=n1,n2
	 cord(1,n)=(1.d0-p)*cord(1,nc)+p*cord(1,n)
	 cord(2,n)=(1.d0-p)*cord(2,nc)+p*cord(2,n)
	 cord(3,n)=(1.d0-p)*cord(3,nc)+p*cord(3,n)
	enddo
200   continue
C==========================================================================

      END
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE LiftGraft(CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION CORD(3,*)

C Lifting graft
      nstart=7039-66

	d=(cord(1,nstart+66)-cord(1,nstart))/5.0

      do 100 i=1,12      
      n1=nstart+(i-1)*66+1
	n2=n1+65
	p=d*i
	k=0
      do n=n1,n2
	k=k+1
	nc=nstart+k
	 cord(1,n)=cord(1,nc)+p
	 cord(2,n)=cord(2,nc)
	 cord(3,n)=cord(3,nc)
	enddo
100   continue

      do n=7700,7765
	 cord(1,n)=cord(1,7765)
	enddo


      END
C==========================================================================
C==========================================================================
