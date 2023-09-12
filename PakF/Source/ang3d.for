C==========================================================================
      SUBROUTINE ang3d2(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK,NEL,NDIM,ITER,NET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
      DIMENSION ID(6,*),NEL(NDIM+1,*)
      DIMENSION NIZ(100)
	DIMENSION N(3)


C      IF (KKORAK.EQ.1) RETURN
C	IF (ITER.GT.0) RETURN


      PI=4.D0*DATAN(1.D0)


      THRESHOLD=0.02D0

      IFILE=51
      OPEN(IFILE,FILE='ANG3D.TXT')


      IU=NEL(8,1)-NEL(5,1)-1
	NL=NEL(1,1)-NEL(5,1)

      NLAYER=NPT/NL

	 AMIN=1.D10 
	 AMAX=-1.D10
      DO I=1,NPT
	  VX=GNODE(2,1,I)
	  VY=GNODE(2,2,I)
	  VZ=GNODE(2,3,I)
	V=DSQRT(VX**2+VY**2+VZ**2)
	 IF (V.LE.AMIN) AMIN=V
	 IF (V.GT.AMAX) AMAX=V
	ENDDO

      THR=THRESHOLD*(AMAX-AMIN)+AMIN

c     IU=5 

      IN=2*(IU+1)*(IU+1)-(IU+1)
	N(1)=1
	N(2)=2*IN-IU
	N(3)=4*IN-IU-2*(IU+1)


      DO L=1,NLAYER
	WRITE(51,*) 'NLAYER = ',L
      DO K=1,3
	WRITE(51,*) 'PATCH = ',K
	N1=N(K)+(L-1)*NL
	N2=N1+IN-IU*K
      DO I=N1,N2,(IU+1)
	  II=0
	  DO J=I,I+IU
	     II=II+1
	     NIZ (II)=J
	  ENDDO
C	WRITE(51,*) (NIZ(JJ),JJ=1,IU+1)
	  CALL ESTIMATE_RADIUS(NIZ,NEL,GNODE,CORD,CCORD,ID,IU+1,THR,IFILE,
     &NDIM,NET,L,NLAYER)
	ENDDO
	ENDDO
	ENDDO


	CLOSE (IFILE)
c	STOP

      END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE ang3d(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK,NEL,NDIM,ITER,NET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
      DIMENSION ID(6,*),NEL(NDIM+1,*)
      DIMENSION NIZ(100)
	DIMENSION N(3)


      IF (KKORAK.EQ.1) RETURN
	IF (ITER.GT.0) RETURN


      PI=4.D0*DATAN(1.D0)


      THRESHOLD=0.02D0

      IFILE=51
      OPEN(IFILE,FILE='ANG3D.TXT')


      IU=NEL(8,1)-NEL(5,1)-1
	NL=NEL(1,1)-NEL(5,1)

      NLAYER=NPT/NL

	 AMIN=1.D10 
	 AMAX=-1.D10
      DO I=1,NPT
	  VX=GNODE(2,1,I)
	  VY=GNODE(2,2,I)
	  VZ=GNODE(2,3,I)
	V=DSQRT(VX**2+VY**2+VZ**2)
	 IF (V.LE.AMIN) AMIN=V
	 IF (V.GT.AMAX) AMAX=V
	ENDDO

      THR=THRESHOLD*(AMAX-AMIN)+AMIN

c     IU=5 

      IN=2*(IU+1)*(IU+1)-(IU+1)
	N(1)=1
	N(2)=2*IN-IU
	N(3)=4*IN-IU-2*(IU+1)


      DO L=1,NLAYER
	WRITE(51,*) 'NLAYER = ',L
      DO K=1,3
	WRITE(51,*) 'PATCH = ',K
	N1=N(K)+(L-1)*NL
	N2=N1+IN-IU*K
      DO I=N1,N2,(IU+1)
	  II=0
	  DO J=I,I+IU
	     II=II+1
	     NIZ (II)=J
	  ENDDO
C	WRITE(51,*) (NIZ(JJ),JJ=1,IU+1)
	  CALL ESTIMATE_MESH(NIZ,NEL,GNODE,CORD,CCORD,ID,IU+1,THR,IFILE,
     &NDIM,NET)
	ENDDO
	ENDDO
	ENDDO


	CLOSE (IFILE)
c	STOP

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ESTIMATE_MESH(NIZ,NEL,GNODE,CORD,CCORD,ID,N,THR,IFILE,
     &NDIM,NET)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*)
      DIMENSION ID(6,*)
	DIMENSION NIZ(*),NEL(NDIM+1,*)
	DIMENSION KONT(100),CC(3)

      KK=0

      WRITE(IFILE,*)'velocities: '

      DO I=1,N
	  NODE=NIZ(I) 
	  VX=GNODE(2,1,NODE)
	  VY=GNODE(2,2,NODE)
	  VZ=GNODE(2,3,NODE)
	  V=DSQRT(VX**2+VY**2+VZ**2)
	  
      WRITE(IFILE,*)'V= ',V,i,thr 
        IF (V.GE.THR) THEN
	     KK=KK+1
	     KONT(KK)=NODE
	  ENDIF
	ENDDO




      IF (KK.LE.1) THEN
c      du1=0.4d0
c      du2=0.8d0
c      goto 100

	KK=0
C	RETURN
        DO I=1,N
	  NODE=NIZ(I) 
	
	M=0
	VS=0.D0
	DO J=1,NET
	IMA=0
	  DO JJ=1,NDIM
	   IF (NODE.EQ.NEL(JJ,J)) THEN
	      IMA=1
	   ENDIF
	  ENDDO
	IF (IMA.EQ.1) THEN
	VS=0.D0
	 DO JJ=1,NDIM
	  VX=GNODE(2,1,JJ)
	  VY=GNODE(2,2,JJ)
	  VZ=GNODE(2,3,JJ)
	  V=DSQRT(VX**2+VY**2+VZ**2)
	  VS=VS+V
	 ENDDO
	VS=VS/8.D0
	V=V+VS
	M=M+1
	ENDIF
      ENDDO

	IF (M.NE.0) THEN 
	  V=V/(1.D0*M)
        IF (V.GE.THR) THEN
	     KK=KK+1
	     KONT(KK)=NODE
	  ENDIF
	ENDIF

	ENDDO
      
	IF (KK.LE.1) RETURN
	  
      ENDIF
	



      U1=0.D0
	UN=0.D0
      DUU=1.D0/N
      IF (KONT(1).GT.NIZ(1)) THEN
	  NOD=KONT(1)-1
        V0X=GNODE(2,1,NOD)
        V0Y=GNODE(2,2,NOD)
        V0Z=GNODE(2,3,NOD)
	  V0=DSQRT(V0X**2+V0Y**2+V0Z**2)
	  
	  NOD=KONT(1)
        VX1=GNODE(2,1,NOD)
        VY1=GNODE(2,2,NOD)
        VZ1=GNODE(2,3,NOD)
	  V1=DSQRT(VX1**2+VY1**2+VZ1**2)
	  U1=DUU*(1.D0-(THR-V0)/(V1-V0))
	ENDIF

      IF (KONT(KK).LT.NIZ(N)) THEN
	  NOD=KONT(KK)+1
        V0X=GNODE(2,1,NOD)
        V0Y=GNODE(2,2,NOD)
        V0Z=GNODE(2,3,NOD)
	  V0=DSQRT(V0X**2+V0Y**2+V0Z**2)
	  
	  NOD=KONT(KK)
        VX1=GNODE(2,1,NOD)
        VY1=GNODE(2,2,NOD)
        VZ1=GNODE(2,3,NOD)
	  V1=DSQRT(VX1**2+VY1**2+VZ1**2)
	  UN=DUU*(1.D0-(THR-V0)/(V1-V0))
	ENDIF

      WRITE(IFILE,*)'CHOOSEN NODES'
      WRITE(IFILE,*) (KONT(I),I=1,KK)

      DU1=(1.D0*(KONT(1)-NIZ(1)))/(1.D0*(NIZ(N)-NIZ(1)))-U1
      DU2=(1.D0*(KONT(KK)-NIZ(1)))/(1.D0*(NIZ(N)-NIZ(1)))+UN



100   DU=DU1
	DUU=(DU2-DU1)/(N-1)


      WRITE(IFILE,*) 'DU1= DU2= ',DU1,DU2

      DO I=1,N
	CC(1)=0.D0
	CC(2)=0.D0
	CC(3)=0.D0
	 DO K=1,N
	  CC(1)=CC(1)+CORD(1,NIZ(K))*BEZ(K-1,N-1,DU)
	  CC(2)=CC(2)+CORD(2,NIZ(K))*BEZ(K-1,N-1,DU)
	  CC(3)=CC(3)+CORD(3,NIZ(K))*BEZ(K-1,N-1,DU)
	 ENDDO
	DU=DU+DUU

        CCORD(1,NIZ(I))=CC(1)
        CCORD(2,NIZ(I))=CC(2)
        CCORD(3,NIZ(I))=CC(3)

	ENDDO

    
      

      END
C==========================================================================
C==========================================================================
      SUBROUTINE ESTIMATE_radius(NIZ,NEL,GNODE,CORD,CCORD,ID,N,THR,IFILE
     &,NDIM,NET,LAYER,NLAYER)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*)
      DIMENSION ID(6,*)
	DIMENSION NIZ(*),NEL(NDIM+1,*)
	DIMENSION KONT(100),CC(3)

      KK=0
	THR=0.D0

      WRITE(IFILE,*)'velocities: '

      DO I=1,N
	  NODE=NIZ(I) 
	  VX=GNODE(2,1,NODE)
	  VY=GNODE(2,2,NODE)
	  VZ=GNODE(2,3,NODE)
	  V=DSQRT(VX**2+VY**2+VZ**2)
	  
      WRITE(IFILE,*)'V= ',V,i,thr 
        IF (V.GE.THR) THEN
	     KK=KK+1
	     KONT(KK)=NODE
	  ENDIF
	ENDDO




      IF (KK.LE.1) THEN
c      du1=0.4d0
c      du2=0.8d0
c      goto 100

	KK=0
C	RETURN
        DO I=1,N
	  NODE=NIZ(I) 
	
	M=0
	VS=0.D0
	DO J=1,NET
	IMA=0
	  DO JJ=1,NDIM
	   IF (NODE.EQ.NEL(JJ,J)) THEN
	      IMA=1
	   ENDIF
	  ENDDO
	IF (IMA.EQ.1) THEN
	VS=0.D0
	 DO JJ=1,NDIM
	  VX=GNODE(2,1,JJ)
	  VY=GNODE(2,2,JJ)
	  VZ=GNODE(2,3,JJ)
	  V=DSQRT(VX**2+VY**2+VZ**2)
	  VS=VS+V
	 ENDDO
	VS=VS/8.D0
	V=V+VS
	M=M+1
	ENDIF
      ENDDO

	IF (M.NE.0) THEN 
	  V=V/(1.D0*M)
        IF (V.GE.THR) THEN
	     KK=KK+1
	     KONT(KK)=NODE
	  ENDIF
	ENDIF

	ENDDO
      
	IF (KK.LE.1) RETURN
	  
      ENDIF
	



      U1=0.D0
	UN=0.D0
      DUU=1.D0/N
      IF (KONT(1).GT.NIZ(1)) THEN
	  NOD=KONT(1)-1
        V0X=GNODE(2,1,NOD)
        V0Y=GNODE(2,2,NOD)
        V0Z=GNODE(2,3,NOD)
	  V0=DSQRT(V0X**2+V0Y**2+V0Z**2)
	  
	  NOD=KONT(1)
        VX1=GNODE(2,1,NOD)
        VY1=GNODE(2,2,NOD)
        VZ1=GNODE(2,3,NOD)
	  V1=DSQRT(VX1**2+VY1**2+VZ1**2)
	  U1=DUU*(1.D0-(THR-V0)/(V1-V0))
	ENDIF

      IF (KONT(KK).LT.NIZ(N)) THEN
	  NOD=KONT(KK)+1
        V0X=GNODE(2,1,NOD)
        V0Y=GNODE(2,2,NOD)
        V0Z=GNODE(2,3,NOD)
	  V0=DSQRT(V0X**2+V0Y**2+V0Z**2)
	  
	  NOD=KONT(KK)
        VX1=GNODE(2,1,NOD)
        VY1=GNODE(2,2,NOD)
        VZ1=GNODE(2,3,NOD)
	  V1=DSQRT(VX1**2+VY1**2+VZ1**2)
	  UN=DUU*(1.D0-(THR-V0)/(V1-V0))
	ENDIF

      WRITE(IFILE,*)'CHOOSEN NODES'
      WRITE(IFILE,*) (KONT(I),I=1,KK)

      DU1=(1.D0*(KONT(1)-NIZ(1)))/(1.D0*(NIZ(N)-NIZ(1)))-U1
      DU2=(1.D0*(KONT(KK)-NIZ(1)))/(1.D0*(NIZ(N)-NIZ(1)))+UN



100   DU=DU1

      DU1=0.0D0
	DU2=1.0D0

	IF (LAYER.GT.0.5*NLAYER) then
	 DU1=0.004D0*(8-(NLAYER-LAYER))**2
	endif
	
	DU=DU1

	DUU=(DU2-DU1)/(N-1)

      WRITE(IFILE,*) 'DU1= DU2= ',DU1,DU2

      DO I=1,N
	CC(1)=0.D0
	CC(2)=0.D0
	CC(3)=0.D0
	 DO K=1,N
	  CC(1)=CC(1)+CORD(1,NIZ(K))*BEZ(K-1,N-1,DU)
	  CC(2)=CC(2)+CORD(2,NIZ(K))*BEZ(K-1,N-1,DU)
	  CC(3)=CC(3)+CORD(3,NIZ(K))*BEZ(K-1,N-1,DU)
	 ENDDO
	DU=DU+DUU

        CCORD(1,NIZ(I))=CC(1)
        CCORD(2,NIZ(I))=CC(2)
        CCORD(3,NIZ(I))=CC(3)

	ENDDO

    
      

      END
C==========================================================================
C==========================================================================
      DOUBLE PRECISION FUNCTION BEZ(I,N,U)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
  	   BEZ=CCC(N,I)*U**I*(1-U)**(N-I)
      END
C==========================================================================
C==========================================================================
      DOUBLE PRECISION FUNCTION CCC(N,I)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
        CCC=FACT(N)/FACT(I)/FACT(N-I)
       
      END
C==========================================================================
C==========================================================================
      DOUBLE PRECISION FUNCTION FACT(N)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
	
	P=1.D0
      DO I=N,1,-1
	  P=P*I
	ENDDO
	FACT=P
	RETURN

      END
C==========================================================================
C==========================================================================
      SUBROUTINE RESIST_bif(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK,NEL,NDIM,ITER,NET,IDPRIT,PRIT,PRES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
	DIMENSION ID(6,*),NEL(NDIM+1,*),PRIT(IDPRIT,*),PRES(3,*)
      DIMENSION NIZ(100)
	DIMENSION N(3),P(3),NN(3),Q(3),NN2(8)
	DIMENSION NID(8),NWALL(6), ITR(6,4), M(8)


C      IF (IAKIRA.EQ.0) RETURN
C      IF (KKORAK.EQ.1) RETURNC
c	IF (ITER.GT.0) RETURN


      PI=4.D0*DATAN(1.D0)


      THRESHOLD=0.02D0

      IFILE=51
      OPEN(IFILE,FILE='ANG3D.TXT')


      IU=NEL(8,1)-NEL(5,1)-1
	NL=NEL(1,1)-NEL(5,1)

      NLAYER=NPT/(NL-1)-1
	NLE=NET/NLAYER


      
	NN(1)=2*(IU)*(IU)+1
	NN(2)=NN(1)+4*(IU)*(IU)
	NN(3)=NN(2)+4*(IU)*(IU)



      DO KKK=0,2*IU-1
	  
	N(1)=2*(IU)*(IU)+1+IU*KKK
	N(2)=N(1)+4*(IU)*(IU)
	N(3)=N(2)+4*(IU)*(IU)


      P(1)=0.D0    
      P(2)=0.D0    
      P(3)=0.D0    
      Q(1)=0.D0    
      Q(2)=0.D0    
      Q(3)=0.D0    
      DO K=1,3
	TOTA=0.D0

c	WRITE(IFILE,*) 'TUBE = ',K
      DO L=1,NLAYER
C	WRITE(IFILE,*) 'NLAYER = ',L
C      WRITE(IFILE,*) 'PATCH = ',K
	N1=N(K)
	N2=N1+IU-1
	VSRE=0.D0
        DO I=N1,N2
	   II=I+(L-1)*NLE
c	WRITE(IFILE,*) 'ELEMENTS = ',II,NEL(1,II),PRIT(1,II)
c	WRITE(IFILE,*) 'ELEMENTS = ',II
	DA=SURF(NEL,II,CORD,NDIM,NEL(1,II),NEL(4,II),NEL(8,II),NEL(5,II))
	TOTA=TOTA+DA
	
	if (k.eq.2) ALL=DABS(CORD(1,NEL(2,II)))
	
	VSR=0.D0
	DO JJ=1,NDIM
	 VSR=VSR+DSQRT(GNODE(2,1,NEL(JJ,II))**2+GNODE(2,2,NEL(JJ,II))**2+
     & GNODE(2,3,NEL(JJ,II))**2)*(1.D0/NDIM)
	ENDDO
	
	   P(K)=P(K)+PRIT(1,II)
	   Q(K)=Q(K)+VSR*DA
c	   VSRE=VSRE+VSR
	  ENDDO
      ENDDO
	P(K)=P(K)/(1.D0*IU*NLAYER)
c	Q(K)=TOTA*VSRE/(1.D0*IU*NLAYER)
C	WRITE(IFILE,*) P(K)
	ENDDO

      TAU=0.D0
      PIPE=0.D0
	SURFTAU=0.D0
	TAUX=0.D0
	TAUY=0.D0
	TAUZ=0.D0
	IBR=0
C      do 110 nbrel=1,net
	  
	DO 110 III=0,NLAYER-1
      DO 100 NBREL=1+NLE*III,NLE*(III+1)
	  
	   IIDOWN1=NN(1)+NLE*(III)
	   IIDOWN2=NN(2)+NLE*(III)
	   IIDOWN3=NN(3)+NLE*(III)
	  
	   IDOWN1=N(1)+IU+NLE*(III)
	   IUP1=NN(1)+2*IU*IU-1+NLE*(III)
	   IDOWN2=N(2)+IU+NLE*(III)
	   IUP2=NN(2)+2*IU*IU-1+NLE*(III)
	   IDOWN3=N(3)+IU+NLE*(III)
	   IUP3=NN(3)+2*IU*IU-1+NLE*(III)

	  IF (
     & (NBREL.GE.IDOWN1.AND.NBREL.LE.IUP1).OR.
     & (NBREL.GE.IDOWN2.AND.NBREL.LE.IUP2).OR.
     & (NBREL.GE.IDOWN3.AND.NBREL.LE.IUP3)) GOTO 100

        IPIPE=0
	  IF (
     & (NBREL.GE.IIDOWN1.AND.NBREL.LT.IDOWN1).OR.
     & (NBREL.GE.IIDOWN2.AND.NBREL.LT.IDOWN2).OR.
     & (NBREL.GE.IIDOWN3.AND.NBREL.LT.IDOWN3)) IPIPE=1


      xc=0.d0
	DO 11 I=1,NDIM
	  M(I)=NEL(I,NBREL)
	  NODE=M(I)
        NID(I)=0
	xc=xc+cord(1,node)
	   IF (DABS(GNODE(2,1,NODE)).GT.1D-10) GOTO 11
      IF(ID(1,M(I)).EQ.0.AND.ID(2,M(I)).EQ.0.AND.ID(3,M(I)).EQ.0)
     &NID(I)=1
11    CONTINUE 
C      if (xc.gt.0.d0) goto 100

      ITR(1,1)=M(8)
      ITR(1,2)=M(4)
      ITR(1,3)=M(5)
      ITR(1,4)=M(1)

      ITR(2,1)=M(7)
      ITR(2,2)=M(3)
      ITR(2,3)=M(6)
      ITR(2,4)=M(2)

      ITR(3,1)=M(6)
      ITR(3,2)=M(2)
      ITR(3,3)=M(5)
      ITR(3,4)=M(1)

      ITR(4,1)=M(7)
      ITR(4,2)=M(3)
      ITR(4,3)=M(8)
      ITR(4,4)=M(4)

      ITR(5,1)=M(3)
      ITR(5,2)=M(2)
      ITR(5,3)=M(4)
      ITR(5,4)=M(1)

      ITR(6,1)=M(7)
      ITR(6,2)=M(6)
      ITR(6,3)=M(8)
      ITR(6,4)=M(5)



      NWALL(1)=NID(1)*NID(4)*NID(5)*NID(8)
      NWALL(2)=NID(2)*NID(3)*NID(6)*NID(7)
      NWALL(3)=NID(1)*NID(2)*NID(5)*NID(6)
      NWALL(4)=NID(3)*NID(4)*NID(7)*NID(8)
      NWALL(5)=NID(1)*NID(2)*NID(3)*NID(4)
      NWALL(6)=NID(5)*NID(6)*NID(7)*NID(8)

      DO 300 NPOV=1,6
       IF (NWALL(NPOV).EQ.0) GOTO 300

	   DA=SURF(NEL,NBREL,CORD,NDIM,ITR(NPOV,1),
     &ITR(NPOV,2),ITR(NPOV,3),ITR(NPOV,4))
	   
         DO JJ=1,4
	     NN2(JJ)=ITR(NPOV,JJ) 
         ENDDO
         
	   T1=DSQRT(PRES(1,NN2(1))**2+PRES(2,NN2(1))**2+PRES(3,NN2(1))**2)
	   T2=DSQRT(PRES(1,NN2(2))**2+PRES(2,NN2(2))**2+PRES(3,NN2(2))**2)
	   T3=DSQRT(PRES(1,NN2(3))**2+PRES(2,NN2(3))**2+PRES(3,NN2(3))**2)
	   T4=DSQRT(PRES(1,NN2(4))**2+PRES(2,NN2(4))**2+PRES(3,NN2(4))**2)
C	   write(ifile,*)'ss,da',(0.25D0*(T1+T2+T3+T4)),DA
 	   TAU=TAU+0.25D0*(T1+T2+T3+T4)*DA
 	   IF (IPIPE.EQ.1) PIPE=PIPE+0.25D0*(T1+T2+T3+T4)*DA
	   SURFTAU=SURFTAU+DA

300   CONTINUE 


100   CONTINUE
110   CONTINUE
      


      TAUX=TAUX/TOTA
      TAUY=TAUY/TOTA
      TAUZ=TAUZ/TOTA


c	WRITE(IFILE,*) 'KKK = ', KKK
c      RES1=(p(2)-p(1)*DCOS(PI/3.0)-p(3)*DCOS(PI/3.0))*tota
      RES1=(p(2)-p(1)-p(3))*tota
	RES2=TAU
	RES2X=TAUX/Q(2)
	RES2Y=TAUY/Q(2)
	RES2Z=TAUZ/Q(2)
C	WRITE(IFILE,*) '(P(2)-(P(1)+P(3)))/Q(2)= ',(P(2)-(P(1)+P(3)))/Q(2)
C	WRITE(IFILE,*) 'TAU/Q(2)= ', TAU/Q(2)

      F1X=P(1)*DCOS(PI/6.D0)
      F1Y=P(1)*DCOS(PI/3.D0)
      F3X=P(3)*DCOS(PI/6.D0)
      F3Y=P(3)*DCOS(PI/3.D0)

      RES1Y=(-P(1)*DCOS(PI/6.0)/Q(1)+P(3)*DCOS(PI/6.0)/Q(3))
      RES1X=P(2)/Q(2)-P(1)*DCOS(PI/3.0)/Q(1)-P(3)*DCOS(PI/3.0)/Q(3)
      
c      WRITE(IFILE,'(8D13.5)')0.02D0-ALL,2.d0*RES1,2.d0*RES2,p(2),SURFTAU
      WRITE(IFILE,'(8D13.5)')ALL,2.d0*RES1,2.d0*RES2,p(2),2.D0*PIPE
C      WRITE(IFILE,'(7D13.5)')RES1,DSQRT(RES1X**2+RES1Y**2)
C      WRITE(IFILE,'(7D13.5)')RES2,DSQRT(RES2X**2+RES2Y**2)
      ENDDO

	CLOSE (IFILE)
c	STOP

      END
C==========================================================================
C==========================================================================
      SUBROUTINE RESIST_tube(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK,NEL,NDIM,ITER,NET,IDPRIT,PRIT,PRES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
	DIMENSION ID(6,*),NEL(NDIM+1,*),PRIT(IDPRIT,*),PRES(3,*)
      DIMENSION NIZ(100)
	DIMENSION N(3),P(3),NN(8),Q(3),TOTAT(2)
	DIMENSION NID(8),NWALL(6), ITR(6,4), M(8),NN2(8)


C      IF (IAKIRA.EQ.0) RETURN
C      IF (KKORAK.EQ.1) RETURNC
c	IF (ITER.GT.0) RETURN


      PI=4.D0*DATAN(1.D0)


      THRESHOLD=0.02D0

      IFILE=51
      OPEN(IFILE,FILE='ANG3D.TXT')


      IU=NEL(6,1)-NEL(5,1)-1

      IESLOJ=NET/100

      DO KKK=1,99
	  
	N(1)=1
	N(2)=N(1)+IESLOJ*KKK


	ALL=DABS(CORD(3,NEL(8,N(1)))-CORD(3,NEL(1,N(2))))

      P(1)=0.D0    
      P(2)=0.D0    
      P(3)=0.D0    
      Q(1)=0.D0    
      Q(2)=0.D0    
      Q(3)=0.D0    

      DO K=1,2
	TOTAT(K)=0.D0
	TOTA=0.D0

	N1=N(K)
	N2=N1+IESLOJ-1
	VSRE=0.D0
        DO I=N1,N2
	   II=I
c	WRITE(IFILE,*) 'ELEMENTS = ',II,NEL(1,II),PRIT(1,II)
c	WRITE(IFILE,*) 'ELEMENTS = ',II
	DA1=SURF(NEL,II,CORD,NDIM,NEL(1,II),NEL(2,II),NEL(3,II),NEL(4,II))
	DA2=SURF(NEL,II,CORD,NDIM,NEL(5,II),NEL(6,II),NEL(7,II),NEL(8,II))
	DA=0.5D0*(DA1+DA2)
c	IF (K.EQ.1) TOTA=TOTA+DA
      TOTA=TOTA+DA
	TOTAT(K)=TOTAT(K)+DA
	VSR=0.D0
	DO JJ=1,NDIM
	 VSR=VSR+DSQRT(GNODE(2,1,NEL(JJ,II))**2+GNODE(2,2,NEL(JJ,II))**2+
     & GNODE(2,3,NEL(JJ,II))**2)*(1.D0/NDIM)
	ENDDO
	   P(K)=P(K)+PRIT(1,II)
	   Q(K)=Q(K)+VSR*DA
	  ENDDO

	P(K)=P(K)/(1.D0*IESLOJ)
C	Q(K)=4.D0*TOTA*VSRE/(1.D0*IESLOJ)
	Q(K)=4.D0*Q(K)
	ENDDO
C	WRITE(IFILE,*) P(K)



      TAU=0.D0
	TAUX=0.D0
	TAUY=0.D0
	TAUZ=0.D0
	IBR=0
	SURFTAU=0.D0
 
 
C      DO 100 NBREL=1,NET
      DO 100 NBREL=N(1),N(2)+IESLOJ-1

	   IBR=0



	DO 11 I=1,NDIM
	  M(I)=NEL(I,NBREL)
	  NODE=M(I)
        NID(I)=0
	   IF (DABS(GNODE(2,3,NODE)).GT.1D-10) GOTO 11
      IF(ID(1,M(I)).EQ.0.AND.ID(2,M(I)).EQ.0.AND.ID(3,M(I)).EQ.0)
     &NID(I)=1
11    CONTINUE 

      ITR(1,1)=M(8)
      ITR(1,2)=M(4)
      ITR(1,3)=M(5)
      ITR(1,4)=M(1)

      ITR(2,1)=M(7)
      ITR(2,2)=M(3)
      ITR(2,3)=M(6)
      ITR(2,4)=M(2)

      ITR(3,1)=M(6)
      ITR(3,2)=M(2)
      ITR(3,3)=M(5)
      ITR(3,4)=M(1)

      ITR(4,1)=M(7)
      ITR(4,2)=M(3)
      ITR(4,3)=M(8)
      ITR(4,4)=M(4)

      ITR(5,1)=M(3)
      ITR(5,2)=M(2)
      ITR(5,3)=M(4)
      ITR(5,4)=M(1)

      ITR(6,1)=M(7)
      ITR(6,2)=M(6)
      ITR(6,3)=M(8)
      ITR(6,4)=M(5)



      NWALL(1)=NID(1)*NID(4)*NID(5)*NID(8)
      NWALL(2)=NID(2)*NID(3)*NID(6)*NID(7)
      NWALL(3)=NID(1)*NID(2)*NID(5)*NID(6)
      NWALL(4)=NID(3)*NID(4)*NID(7)*NID(8)
      NWALL(5)=NID(1)*NID(2)*NID(3)*NID(4)
      NWALL(6)=NID(5)*NID(6)*NID(7)*NID(8)

      DO 300 NPOV=1,6
       IF (NWALL(NPOV).EQ.0) GOTO 300

	   DA=SURF(NEL,NBREL,CORD,NDIM,ITR(NPOV,1),
     &ITR(NPOV,2),ITR(NPOV,3),ITR(NPOV,4))
         DO JJ=1,4
	     NN2(JJ)=ITR(NPOV,JJ) 
         ENDDO
	   T1=DSQRT(PRES(1,NN2(1))**2+PRES(2,NN2(1))**2+PRES(3,NN2(1))**2)
	   T2=DSQRT(PRES(1,NN2(2))**2+PRES(2,NN2(2))**2+PRES(3,NN2(2))**2)
	   T3=DSQRT(PRES(1,NN2(3))**2+PRES(2,NN2(3))**2+PRES(3,NN2(3))**2)
	   T4=DSQRT(PRES(1,NN2(4))**2+PRES(2,NN2(4))**2+PRES(3,NN2(4))**2)
C	   write(ifile,*)'ss,da',(0.25D0*(T1+T2+T3+T4)),DA
 	   TAU=TAU+0.25D0*(T1+T2+T3+T4)*DA
	   SURFTAU=SURFTAU+DA

300   CONTINUE 
      goto 100




         DO 10 I=1,NDIM
	      NODE=NEL(I,NBREL)
	      IDX=ID(1,NODE)
	      IDY=ID(2,NODE)
	      IDZ=ID(3,NODE)
	   IF (DABS(GNODE(2,3,NODE)).GT.1D-10) GOTO 10
	     IF (IDX.EQ.0.AND.IDY.EQ.0.AND.IDZ.EQ.0) THEN
	        IBR=IBR+1
		    NN(IBR)=NODE
	     ENDIF
10    CONTINUE
        IF (IBR.LT.4) GOTO 100
        
        IF (IBR.EQ.4) THEN  	  
	   DA=SURF(NEL,NBREL,CORD,NDIM,NN(1),NN(2),NN(3),NN(4))
	   T1=DSQRT(PRES(1,NN(1))**2+PRES(2,NN(1))**2+PRES(3,NN(1))**2)
	   T2=DSQRT(PRES(1,NN(2))**2+PRES(2,NN(2))**2+PRES(3,NN(2))**2)
	   T3=DSQRT(PRES(1,NN(3))**2+PRES(2,NN(3))**2+PRES(3,NN(3))**2)
	   T4=DSQRT(PRES(1,NN(4))**2+PRES(2,NN(4))**2+PRES(3,NN(4))**2)
 	   TAU=TAU+0.25D0*(T1+T2+T3+T4)*DA
	   SURFTAU=SURFTAU+DA
	  ENDIF
        IF (IBR.EQ.6) THEN 
	 WRITE(*,*) 'USAO OVDE BRE'
c	   WRITE(IFILE,*)(NN(JJ),JJ=1,6) 	  
	   DA=SURF(NEL,NBREL,CORD,NDIM,NN(1),NN(2),NN(4),NN(5))
	   T1=DSQRT(PRES(1,NN(1))**2+PRES(2,NN(1))**2+PRES(3,NN(1))**2)
	   T2=DSQRT(PRES(1,NN(2))**2+PRES(2,NN(2))**2+PRES(3,NN(2))**2)
	   T3=DSQRT(PRES(1,NN(4))**2+PRES(2,NN(4))**2+PRES(3,NN(4))**2)
	   T4=DSQRT(PRES(1,NN(5))**2+PRES(2,NN(5))**2+PRES(3,NN(5))**2)
 	   TAU=TAU+0.25D0*(T1+T2+T3+T4)*DA
	   DA=SURF(NEL,NBREL,CORD,NDIM,NN(2),NN(3),NN(5),NN(6))
	   SURFTAU=SURFTAU+DA
	   T1=DSQRT(PRES(1,NN(2))**2+PRES(2,NN(2))**2+PRES(3,NN(2))**2)
	   T2=DSQRT(PRES(1,NN(3))**2+PRES(2,NN(3))**2+PRES(3,NN(3))**2)
	   T3=DSQRT(PRES(1,NN(5))**2+PRES(2,NN(5))**2+PRES(3,NN(5))**2)
	   T4=DSQRT(PRES(1,NN(6))**2+PRES(2,NN(6))**2+PRES(3,NN(6))**2)
 	   TAU=TAU+0.25D0*(T1+T2+T3+T4)*DA
	   SURFTAU=SURFTAU+DA
	  ENDIF
C	WRITE(IFILE,*) 'T1,T2,T3,T4 = ', T1,T2,T3,T4,DA
        ZUN=CORD(3,NN(3))-CORD(3,NN(1))
        YUN=CORD(2,NN(3))-CORD(2,NN(1))
        XUN=CORD(1,NN(3))-CORD(1,NN(1))
        ZZUN=ZUN/DSQRT(XUN**2+YUN**2+ZUN**2) 
	  ZZUN=DSIN(ATAN(10.D0/500.D0))
	  TAUZ=TAUZ-PRIT(1,NBREL)*ZZUN*DA
	  
	  T1X=PRES(1,NN(1))
	  T2X=PRES(1,NN(2))
	  T3X=PRES(1,NN(3))
	  T4X=PRES(1,NN(4))
	  TAUX=TAUX+0.25D0*(T1X+T2X+T3X+T4X)*DA

	  T1Y=PRES(2,NN(1))
	  T2Y=PRES(2,NN(2))
	  T3Y=PRES(2,NN(3))
	  T4Y=PRES(2,NN(4))
	  TAUY=TAUY+0.25D0*(T1Y+T2Y+T3Y+T4Y)*DA

	  T1Z=PRES(3,NN(1))
	  T2Z=PRES(3,NN(2))
	  T3Z=PRES(3,NN(3))
	  T4Z=PRES(3,NN(4))
	  TAUZ=TAUZ+0.25D0*(T1Z+T2Z+T3Z+T4Z)*DA
100   CONTINUE
      
C      TAU=TAU/TOTA




c	WRITE(IFILE,*) 'KKK = ', KKK
c      RES1=(P(1)-P(2))/4.084D-7
      RES1=(P(1)-P(2))/Q(2)
	RES2=TAU/Q(2)
c	RES2=TAU/4.084D-7
      
      
c	RR=CORD(1,1)
	RR=10.D-4
	Vmax=0.14d0
	Vsr=Vmax*0.5d0
	Q2=vsr*RR**2*PI
	ANALR=8.D0*AMI*ALL*Q2/(RR**2)
C	WRITE(IFILE,'(8D13.5)')ALL,RES1,RES2,ANALR
      RES11=2.d0*(P(1)*TOTAT(1)-P(2)*TOTAT(2))
c	RES22=TAUZ
	RES22=2.d0*TAU
	TT=DSQRT(TAUX**2+TAUY**2+TAUZ**2)
	WRITE(IFILE,'(8D13.5)')ALL,RES11,RES22,analr,TAUZ,TT,SURFTAU
C	WRITE(IFILE,'(8D13.5)')ALL,p(1),p(2),TOTAT(1),TOTAT(2)
      
      ENDDO

	CLOSE (IFILE)

      END
C==========================================================================
C==========================================================================
      FUNCTION SURF(NEL,NBREL,CORD,NDIM,N1,N2,N3,N4)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION NEL(NDIM+1,*),CORD(3,*)

	  X1=CORD(1,N1) 
	  Y1=CORD(2,N1) 
	  Z1=CORD(3,N1) 
	  
	  X2=CORD(1,N2) 
	  Y2=CORD(2,N2) 
	  Z2=CORD(3,N2) 
	  
	  X3=CORD(1,N3) 
	  Y3=CORD(2,N3) 
	  Z3=CORD(3,N3) 
	  
	  X4=CORD(1,N4) 
	  Y4=CORD(2,N4) 
	  Z4=CORD(3,N4) 

	  A=DSQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
	  B=DSQRT((X2-X3)**2+(Y2-Y3)**2+(Z2-Z3)**2)
	  C=DSQRT((X1-X3)**2+(Y1-Y3)**2+(Z1-Z3)**2)
	  S=(A+B+C)*0.5D0
	  P1=DSQRT(DABS(S*(S-A)*(S-B)*(S-C)))

	  A=DSQRT((X2-X3)**2+(Y2-Y3)**2+(Z2-Z3)**2)
	  B=DSQRT((X3-X4)**2+(Y3-Y4)**2+(Z3-Z4)**2)
	  C=DSQRT((X2-X4)**2+(Y2-Y4)**2+(Z2-Z4)**2)
	  S=(A+B+C)*0.5D0
	  P2=DSQRT(DABS(S*(S-A)*(S-B)*(S-C)))

        SURF=P1+P2

	END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE RESIST_VOL(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK,NEL,NDIM,ITER,NET,IDPRIT,PRIT,PRES,SPSIL,NETIP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
	DIMENSION ID(6,*),NEL(NDIM+1,*),SPSIL(NETIP,*)
	DIMENSION PRIT(IDPRIT,*),PRES(3,*)
	DIMENSION NID(8),NWALL(6), ITR(6,4), M(8),NN2(8)

C
CE Subroutine ZADNOD is used for inclusion prescribed values
C

      PI=4.D0*DATAN(1.D0)

      GANGLE=PI/2.D0
      GANGLE=PI

      IFILE=51
      OPEN(IFILE,FILE='RESISTANCE.TXT')


      ANGLE1=GANGLE/100.D0
      P1=0.D0
      DO NBREL=1,100
	  P1=P1+PRIT(1,NBREL)
	ENDDO
	P1=P1/100.D0




      Q1=0.D0 
	VSRE=0.D0
	TOTA1=0.D0
      DO II=1,100
	DA1=SURF(NEL,II,CORD,NDIM,NEL(1,II),NEL(2,II),NEL(3,II),NEL(4,II))
	DA2=SURF(NEL,II,CORD,NDIM,NEL(5,II),NEL(6,II),NEL(7,II),NEL(8,II))
		DA=0.5*(DA1+DA2)
		TOTA1=TOTA1+DA
		VSR=0.D0
		DO JJ=1,NDIM
		VSR=VSR+DSQRT(GNODE(2,1,NEL(JJ,II))**2+
     &GNODE(2,2,NEL(JJ,II))**2+GNODE(2,3,NEL(JJ,II))**2)*(1.D0/NDIM)
		ENDDO
	   
		Q1=Q1+VSR*DA
	ENDDO


c start of global loop over element layers               
      N1=1
      N2=116+1
	DX=CORD(1,N2)-CORD(1,N1)
	DY=CORD(2,N2)-CORD(2,N1)
	DZ=CORD(3,N2)-CORD(3,N1)
	DL1=SQRT(DX**2+DY**2+DZ**2)
c      DL=0.d0
      DL=DL1

      DO 10 N=1,99
C      DO 10 N=99,1,-1

      ANGLE=ANGLE1*(N+1)
      FX=0.D0 
      FY=0.D0 
      FZ=0.D0 
	FWALL=0.D0
      N1=(N+1)*116
      N2=(N+2)*116
	DX=CORD(1,N2)-CORD(1,N1)
	DY=CORD(2,N2)-CORD(2,N1)
	DZ=CORD(3,N2)-CORD(3,N1)
	DL=DL+SQRT(DX**2+DY**2+DZ**2)


      DO 5 NODE=1+1*116,(N+1)*116
c      DO 5 NODE=1,(N+2)*116
c      DO 5 NODE=1+(50-N/2)*116,(50+N/2)*116
      
	IF (GNODE(2,3,NODE).GT.1.D-10) GOTO 5
      IDX=ID(1,NODE)
	IDY=ID(2,NODE)
	IDZ=ID(3,NODE)
	IF (IDX.EQ.0.AND.IDY.EQ.0.AND.IDZ.EQ.0) THEN
c	 FX=SPSIL(1,NODE)       
c	 FY=SPSIL(2,NODE)       
c	 FZ=SPSIL(3,NODE)       
	 FX=FX+SPSIL(1,NODE)       
	 FY=FY+SPSIL(2,NODE)       
	 FZ=FZ+SPSIL(3,NODE)       
	ENDIF
C	FWALL =FWALL+dsqrt((FX)**2+(FZ)**2)
5     CONTINUE



      P2=0.D0
      DO NBREL=1+N*100,100+N*100
	  P2=P2+PRIT(1,NBREL)
	ENDDO
	P2=P2/100.D0


      Q2=0.D0 
	VSRE=0.D0
	TOTA2=0.D0
      DO II=1+N*100,100+N*100
	DA1=SURF(NEL,II,CORD,NDIM,NEL(1,II),NEL(2,II),NEL(3,II),NEL(4,II))
	DA2=SURF(NEL,II,CORD,NDIM,NEL(5,II),NEL(6,II),NEL(7,II),NEL(8,II))
		DA=0.5*(DA1+DA2)
		TOTA2=TOTA2+DA
		VSR=0.D0
		DO JJ=1,NDIM
		VSR=VSR+DSQRT(GNODE(2,1,NEL(JJ,II))**2+
     &GNODE(2,2,NEL(JJ,II))**2+GNODE(2,3,NEL(JJ,II))**2)*(1.D0/NDIM)
		ENDDO
	   
		Q2=Q2+VSR*DA
	ENDDO


      TAU=0.D0
      DO 100 NBREL=1,(N+1)*100
C      DO 100 NBREL=1+(49-N/2)*100,(51+N/2)*100
	IBR=0



	DO 11 I=1,NDIM
	  M(I)=NEL(I,NBREL)
	  NODE=M(I)
        NID(I)=0
	   IF (DABS(GNODE(2,3,NODE)).GT.1D-10) GOTO 11
      IF(ID(1,M(I)).EQ.0.AND.ID(2,M(I)).EQ.0.AND.ID(3,M(I)).EQ.0)
     &NID(I)=1
11    CONTINUE 

      ITR(1,1)=M(8)
      ITR(1,2)=M(4)
      ITR(1,3)=M(5)
      ITR(1,4)=M(1)

      ITR(2,1)=M(7)
      ITR(2,2)=M(3)
      ITR(2,3)=M(6)
      ITR(2,4)=M(2)

      ITR(3,1)=M(6)
      ITR(3,2)=M(2)
      ITR(3,3)=M(5)
      ITR(3,4)=M(1)

      ITR(4,1)=M(7)
      ITR(4,2)=M(3)
      ITR(4,3)=M(8)
      ITR(4,4)=M(4)

      ITR(5,1)=M(3)
      ITR(5,2)=M(2)
      ITR(5,3)=M(4)
      ITR(5,4)=M(1)

      ITR(6,1)=M(7)
      ITR(6,2)=M(6)
      ITR(6,3)=M(8)
      ITR(6,4)=M(5)



      NWALL(1)=NID(1)*NID(4)*NID(5)*NID(8)
      NWALL(2)=NID(2)*NID(3)*NID(6)*NID(7)
      NWALL(3)=NID(1)*NID(2)*NID(5)*NID(6)
      NWALL(4)=NID(3)*NID(4)*NID(7)*NID(8)
      NWALL(5)=NID(1)*NID(2)*NID(3)*NID(4)
      NWALL(6)=NID(5)*NID(6)*NID(7)*NID(8)

      DO 300 NPOV=1,6
       IF (NWALL(NPOV).EQ.0) GOTO 300

	   DA=SURF(NEL,NBREL,CORD,NDIM,ITR(NPOV,1),
     &ITR(NPOV,2),ITR(NPOV,3),ITR(NPOV,4))
         DO JJ=1,4
	     NN2(JJ)=ITR(NPOV,JJ) 
         ENDDO
	   T1=DSQRT(PRES(1,NN2(1))**2+PRES(2,NN2(1))**2+PRES(3,NN2(1))**2)
	   T2=DSQRT(PRES(1,NN2(2))**2+PRES(2,NN2(2))**2+PRES(3,NN2(2))**2)
	   T3=DSQRT(PRES(1,NN2(3))**2+PRES(2,NN2(3))**2+PRES(3,NN2(3))**2)
	   T4=DSQRT(PRES(1,NN2(4))**2+PRES(2,NN2(4))**2+PRES(3,NN2(4))**2)
C	   write(ifile,*)'ss,da',(0.25D0*(T1+T2+T3+T4)),DA
 	   TAU=TAU+0.25D0*(T1+T2+T3+T4)*DA
C	   SURFTAU=SURFTAU+DA

300   CONTINUE 

100   CONTINUE


c==========================================================================
C     Analytical solution for curved tube, Burger... et al 1987
c==========================================================================
      VMAX=1.8D3
      A=10.D-4
	R=20.D0*A
      AK=(2.D0*A/R)*(VMAX*A/0.035D0)**2
	AN=1.D0-0.0306D0*(AK/576.D0)**2+0.0120*(AK/576.D0)**4
c     data for straight tube  
c      ANGLE=0.D0
c      AN=1.D0
      TOTAA=A**2*PI
C	DLA=(R*PI/100.D0)*(N+1)
	DLA=DL
	ANAL=((VMAX*4.D0*3.675D-2*(DLA)/(A**2))*TOTAA)/AN
c==========================================================================
	Q2=2.D0*Q2
	


c==========================================================================
C     Analytical solution for straight tube
c==========================================================================

C	RR=10.D-4
C	Vsr=0.07d0*0.5d0
C	ANAL=8.D0*AMI*(DL)*PI*Vsr
	
c      RES1Z=P1*TOTA1
c      RES1X=0.D0
c      RES2Z=-P2*TOTA2*DCOS(ANGLE)
c      RES2X=P2*TOTA2*DSIN(ANGLE)
c	RESX=(RES1X+RES2X)
c	RESZ=(RES1Z+RES2Z)
c	ANUMER1=dsqrt((RES1X*1.D1)**2+(RES1Z*1.D1)**2)-
c     &dsqrt((RES2X*1.D1)**2+(RES2Z*1.D1)**2)
c	ANUMER2=dsqrt((FX*1.D1)**2+(FZ*1.D1)**2)
c     
c	WRITE(IFILE,'(9D13.5)')
c     &DL,RESX*1.D1,RESZ*1.D1,FX*1.D1,FZ*1.D1,
c     &ANUMER1,(P1-P2)*TOTA1*1.D1,
c     &ANUMER2,ANAL



      RES1Z=P1*TOTA1*2.D0
      RES1X=0.D0
      RES2Z=-P2*2.D0*TOTA2*DCOS(ANGLE)
      RES2X=P2*2.D0*TOTA2*DSIN(ANGLE)
      RESX=(RES1X+RES2X)
	RESZ=(RES1Z+RES2Z)
	ANUMER1=dsqrt((RES1X)**2+(RES1Z)**2)-
     &dsqrt((RES2X)**2+(RES2Z)**2)
	ANUMER2=dsqrt((2.D0*FX)**2+(2.D0*FZ)**2)
c	ANUMER2=FWALL
c     
c	WRITE(IFILE,'(9D13.5)')
c     &DL,TAU,RESZ,2.D0*FX,2.D0*FZ,
c     &ANUMER1,(P1-P2)*2.D0*TOTA1,
c     &ANUMER2,ANAL

	WRITE(IFILE,'(9D13.5)')DL,DLA,2.D0*TAU,
     &(P1-P2)*2.D0*TOTA1,ANUMER2,ANAL


10    CONTINUE

	CLOSE (IFILE)

      END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE RESIST_bif2(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK,NEL,NDIM,ITER,NET,IDPRIT,PRIT,PRES,SPSIL,NETIP,
     &INDEL)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
	DIMENSION ID(6,*),NEL(NDIM+1,*),INDEL(*)
	DIMENSION PRIT(IDPRIT,*),PRES(3,*)
	DIMENSION SPSIL(NETIP,3)
      DIMENSION NIZ(100),NNIZ(4)
	DIMENSION N(3),P(3),NN(3),Q(3)


C      IF (IAKIRA.EQ.0) RETURN
C      IF (KKORAK.EQ.1) RETURNC
c	IF (ITER.GT.0) RETURN


      PI=4.D0*DATAN(1.D0)




C      DO I=1,NPT
C       INDEL(I)=0
C      ENDDO

C       DO NBREL=1,NET
C        DO I=1,NDIM
C         NODE=NEL(I,NBREL)
C         INDEL(NODE)=INDEL(NODE)+1
C        ENDDO        
C       ENDDO

      THRESHOLD=0.02D0

      IFILE=51
      OPEN(IFILE,FILE='ANG3D.TXT')


      IU=NEL(8,1)-NEL(5,1)-1
	NL=NEL(1,1)-NEL(5,1)

      NLAYER=NPT/(NL-1)-1
	NLE=NET/NLAYER


      
	NN(1)=2*(IU)*(IU)+1
	NN(2)=NN(1)+4*(IU)*(IU)
	NN(3)=NN(2)+4*(IU)*(IU)





      DO KKK=0,2*IU-1
	  
	N(1)=2*(IU)*(IU)+1+IU*KKK
	N(2)=N(1)+4*(IU)*(IU)
	N(3)=N(2)+4*(IU)*(IU)


      P(1)=0.D0    
      P(2)=0.D0    
      P(3)=0.D0    
      Q(1)=0.D0    
      Q(2)=0.D0    
      Q(3)=0.D0    
      DO K=1,3
	TOTA=0.D0

c	WRITE(IFILE,*) 'TUBE = ',K
      DO L=1,NLAYER
C	WRITE(IFILE,*) 'NLAYER = ',L
C      WRITE(IFILE,*) 'PATCH = ',K
	N1=N(K)
	N2=N1+IU-1
	VSRE=0.D0
        DO I=N1,N2
	   II=I+(L-1)*NLE
c	WRITE(IFILE,*) 'ELEMENTS = ',II,NEL(1,II),PRIT(1,II)
c	WRITE(IFILE,*) 'ELEMENTS = ',II
	DA=SURF(NEL,II,CORD,NDIM,NEL(1,II),NEL(4,II),NEL(8,II),NEL(5,II))
	TOTA=TOTA+DA


	
	if (k.eq.2) ALL=2.D0*DABS(CORD(1,NEL(1,II)))
	
	VSR=0.D0
	DO JJ=1,NDIM
	 VSR=VSR+DSQRT(GNODE(2,1,NEL(JJ,II))**2+GNODE(2,2,NEL(JJ,II))**2+
     & GNODE(2,3,NEL(JJ,II))**2)*(1.D0/NDIM)
	ENDDO
	
	   P(K)=P(K)+PRIT(1,II)
	   Q(K)=Q(K)+VSR*DA
c	   VSRE=VSRE+VSR
	  ENDDO
      ENDDO
	P(K)=P(K)/(1.D0*IU*NLAYER)
c	Q(K)=TOTA*VSRE/(1.D0*IU*NLAYER)
C	WRITE(IFILE,*) P(K)
	ENDDO

      TAU=0.D0
	TAUX=0.D0
	TAUY=0.D0
	TAUZ=0.D0
	IBR=0
	 FXWALL =0.D0
	 FYWALL =0.D0
	 FZWALL =0.D0
	 FWALL=0.D0
	do jjj=1,npt
	  id(5,jjj)=0
	enddo
	K=0
	DO 110 III=0,NLAYER-1
      DO 100 NBREL=1+NLE*III,NLE*(III+1)
	   IDOWN=N(1)+IU+NLE*(III)
	   IUP=NN(1)+2*IU*IU-1+NLE*(III)
	  IF (NBREL.GE.IDOWN.AND.NBREL.LE.IUP) GOTO 100
	   IDOWN=N(2)+IU+NLE*(III)
	   IUP=NN(2)+2*IU*IU-1+NLE*(III)
	  IF (NBREL.GE.IDOWN.AND.NBREL.LE.IUP) GOTO 100
	   IDOWN=N(3)+IU+NLE*(III)
	   IUP=NN(3)+2*IU*IU-1+NLE*(III)
	  IF (NBREL.GE.IDOWN.AND.NBREL.LE.IUP) GOTO 100

      K=K+1
	Ewall=0.d0
	kk=0
      DO 5 JJ=1,8
c	 NODE=NNIZ(JJ)
	 NODE=NEL(JJ,NBREL)
c	 IF (INDEL(NODE).GE.4) GOTO 5
c	 WRITE(IFILE,'(I10)') INDEL(NODE)
	 IF (GNODE(2,1,NODE).GT.1.D-10) GOTO 5
       IDX=ID(1,NODE)
	 IDY=ID(2,NODE)
	 IDZ=ID(3,NODE)
	 IF (IDX.EQ.0.AND.IDY.EQ.0.AND.IDZ.EQ.0) THEN
	 id(5,node)=1
	  kk=kk+1
c	 FXWALL=FXWALL+SPSIL(1,NODE)/(1.D0*INDEL(NODE))
c	 FYWALL=FYWALL+SPSIL(2,NODE)/(1.D0*INDEL(NODE))
c	 FZWALL=FZWALL+SPSIL(3,NODE)/(1.D0*INDEL(NODE))
        Ewall=Ewall+
     &dsqrt(SPSIL(1,NODE)**2+SPSIL(2,NODE)**2+SPSIL(3,NODE)**2)
	 ENDIF


	
5      CONTINUE
        
c       if (kk.ne.0) WRITE(IFILE,'(i10)')kk
c	 if (kk.ne.0) FWALL=FWALL+Ewall/(1.d0*kk)


	  
100   CONTINUE
110   CONTINUE
      


      Fwall=0.d0 
      DO 15 node=1,npt
	 IF (GNODE(2,1,NODE).GT.1.D-10) GOTO 15
	 if (id(5,node).eq.0) goto 15
       IDX=ID(1,NODE)
	 IDY=ID(2,NODE)
	 IDZ=ID(3,NODE)
	 IF (IDX.EQ.0.AND.IDY.EQ.0.AND.IDZ.EQ.0) THEN
        Fwall=Fwall+
     &dsqrt(SPSIL(1,NODE)**2+SPSIL(2,NODE)**2+SPSIL(3,NODE)**2)
	 ENDIF
  15  CONTINUE


C      write(ifile,'(8D13.5)')all,tau,
C     &(p(2)-p(1)*DCOS(PI/3.0)-p(3)*DCOS(PI/3.0))*tota,taux,tauy

c      R=12.5D-4 
c	TAU=TAU/(0.5D0*R**2*PI)
C	WRITE(IFILE,*) 'TAU = ', TAU, TOTA,IBR
c      TAU=TAU/TOTA
      TAU=TAU
      TAUX=TAUX/TOTA
      TAUY=TAUY/TOTA
      TAUZ=TAUZ/TOTA


c	WRITE(IFILE,*) 'KKK = ', KKK
c      RES1=(p(2)-p(1)*DCOS(PI/3.0)-p(3)*DCOS(PI/3.0))/Q(2)
      RES1=(p(2)-p(1)*DCOS(PI/3.0)-p(3)*DCOS(PI/3.0))*TOTA
c	RES2=TAU/Q(2)
	RES2=TAU
	RES2X=TAUX/Q(2)
	RES2Y=TAUY/Q(2)
	RES2Z=TAUZ/Q(2)
C	WRITE(IFILE,*) '(P(2)-(P(1)+P(3)))/Q(2)= ',(P(2)-(P(1)+P(3)))/Q(2)
C	WRITE(IFILE,*) 'TAU/Q(2)= ', TAU/Q(2)

      F1X=P(1)*DCOS(PI/6.D0)
      F1Y=P(1)*DCOS(PI/3.D0)
      F3X=P(3)*DCOS(PI/6.D0)
      F3Y=P(3)*DCOS(PI/3.D0)
	F2X=P(2)

      RES1Y=(-P(1)*DCOS(PI/6.0)+P(3)*DCOS(PI/6.0))*TOTA
      RES1X=(P(2)-P(1)*DCOS(PI/3.0)-P(3)*DCOS(PI/3.0))*TOTA
C      RES1X=(P(2))*TOTA
      
c      WRITE(IFILE,'(8D13.5)')ALL,RES1,RES2,p(2),p(1),p(3),q(2),q(1)
C      WRITE(IFILE,'(10D13.5)')ALL,RES1,RES2,p(2),p(1),p(3),
C     &FXWALL,FYWALL,FZWALL
C      WRITE(IFILE,'(10D13.5)')ALL,RES1X,RES1Y,FXWALL,FYWALL
      WRITE(IFILE,'(10D13.5)')ALL,2.d0*RES1X,2.d0*RES1Y,
     &2.d0*FXWALL,2.d0*FYWALL,FWALL
C      WRITE(IFILE,*)'BROJ ELEMENATA= ',K,NLE


C      WRITE(IFILE,'(7D13.5)')RES2,DSQRT(RES2X**2+RES2Y**2)
      ENDDO



	do jjj=1,npt
	  id(5,jjj)=0
	enddo

      GOTO 200

      fxplus=0.d0
      fyplus=0.d0
      fzplus=0.d0
      fxminus=0.d0
      fyminus=0.d0
      fzminus=0.d0
	fx=0.d0
	fx1=0.d0
	fx1p=0.d0
	fx1m=0.d0
	fx2p=0.d0
	fx2m=0.d0
	fx3p=0.d0
	fx3m=0.d0
	k=0
	xmin=1.d10
      do 30 i=1,npt
	  if (cord(1,i).lt.xmin) xmin=cord(1,i)
       IDX=ID(1,i)
	 IDY=ID(2,i)
	 IDZ=ID(3,i)


	  if(spsil(1,i).gt.0.d0) fx1p=fx1p+spsil(1,i)
	  if(spsil(1,i).le.0.d0) fx1m=fx1m+spsil(1,i)
	  if(spsil(2,i).gt.0.d0) fx2p=fx2p+spsil(2,i)
	  if(spsil(2,i).le.0.d0) fx2m=fx2m+spsil(2,i)
	  if(spsil(3,i).gt.0.d0) fx3p=fx3p+spsil(3,i)
	  if(spsil(3,i).le.0.d0) fx3m=fx3m+spsil(3,i)


	 IF (IDX.EQ.0.AND.IDY.EQ.0.AND.IDZ.EQ.0) THEN
	 IF (GNODE(2,1,i).GT.1.D-5) then
	    k=k+1
	    fx=fx+spsil(1,i)
	   GOTO 30
	else
	 fx1=fx1+spsil(1,i)
	endif
	  if(spsil(1,i).gt.0.d0) fxplus=fxplus+spsil(1,i)
	  if(spsil(1,i).le.0.d0) fxminus=fxminus+spsil(1,i)
	  if(spsil(2,i).gt.0.d0) fyplus=fyplus+spsil(2,i)
	  if(spsil(2,i).le.0.d0) fyminus=fyminus+spsil(2,i)
	  if(spsil(3,i).gt.0.d0) fzplus=fzplus+spsil(3,i)
	  if(spsil(3,i).le.0.d0) fzminus=fzminus+spsil(3,i)
	 endif
30     continue
      WRITE(IFILE,'(10D13.5)')ALL,fxplus,fxminus,fyplus,fyminus,
     &fzplus,fzminus

      WRITE(IFILE,'(i10,10d13.5)')k,fx,fx1,fx1p,fx1m,fx2p,fx2m,fx3p,fx3m

      fxplus=0.d0
      fyplus=0.d0
      fzplus=0.d0
      fxminus=0.d0
      fyminus=0.d0
      fzminus=0.d0
      do i=1,npt
       IDX=ID(1,i)
	 IDY=ID(2,i)
	 IDZ=ID(3,i)
	 IF (cord(1,i).eq.xmin) THEN
c	 IF (GNODE(2,1,i).GT.1.D-5) then
c	 IF (IDX.ne.0.or.IDY.ne.0.or.IDZ.ne.0) THEN
	  if(spsil(1,i).gt.0.d0) fxplus=fxplus+spsil(1,i)
	  if(spsil(1,i).le.0.d0) fxminus=fxminus+spsil(1,i)
	  if(spsil(2,i).gt.0.d0) fyplus=fyplus+spsil(2,i)
	  if(spsil(2,i).le.0.d0) fyminus=fyminus+spsil(2,i)
	  if(spsil(3,i).gt.0.d0) fzplus=fzplus+spsil(3,i)
	  if(spsil(3,i).le.0.d0) fzminus=fzminus+spsil(3,i)
	 endif
	enddo
      WRITE(IFILE,'(10D13.5)')ALL,fxplus,fxminus,fyplus,fyminus,
     &fzplus,fzminus

      WRITE(IFILE,'(8D13.5)')ALL,p(2)*tota,p(1)*tota,p(3)*tota


200	CLOSE (IFILE)
c	STOP

      END
C==========================================================================
C==========================================================================
C==========================================================================
      SUBROUTINE RESIST_bif3(GNODE,CORD,CCORD,ID,T,NPT,VMESH,TIME,NSTAC,
     &ZADVRE,AMI,KKORAK,NEL,NDIM,ITER,NET,IDPRIT,PRIT,PRES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),CCORD(3,*),VMESH(3,*),ZADVRE(*)
	DIMENSION ID(6,*),NEL(NDIM+1,*),PRIT(IDPRIT,*),PRES(3,*)
      DIMENSION NIZ(100)
	DIMENSION N(3),P(3),NN(3),Q(3),NN2(8)
	DIMENSION NID(8),NWALL(6), ITR(6,4), M(8)
	DIMENSION T1(2),T2(2),T3(2),T4(2)


C      IF (IAKIRA.EQ.0) RETURN
C      IF (KKORAK.EQ.1) RETURNC
c	IF (ITER.GT.0) RETURN


      PI=4.D0*DATAN(1.D0)


      THRESHOLD=0.02D0

      IFILE=51
      OPEN(IFILE,FILE='ANG3D33.TXT')


      IU=NEL(8,1)-NEL(5,1)-1
	NL=NEL(1,1)-NEL(5,1)

      NLAYER=NPT/(NL-1)-1
	NLE=NET/NLAYER


      
	NN(1)=2*(IU)*(IU)+1
	NN(2)=NN(1)+4*(IU)*(IU)
	NN(3)=NN(2)+4*(IU)*(IU)


      XMAX=-1.D10
      XMIN=1.D10
      DO 50 I=1,NPT
	 IF (DABS(CORD(2,I)).GT.1.D-10) GOTO 50
	 IF (CORD(1,I).GT.XMAX) XMAX=CORD(1,I)
	 IF (CORD(1,I).LT.XMIN) XMIN=CORD(1,I)
50    CONTINUE
      



      DO 150 KKK=0,2*IU-1
C	XX=(KKK+1)*(2.D0*XMIN)/(2.D0*IU)
C	XX=(KKK+1)*(2.01D0*DABS(XMIN))/(2.D0*IU)
	XX=(KKK+1)*(1.01D0*DABS(XMIN))/(2.D0*IU)
	  
	N(1)=2*(IU)*(IU)+1+IU*KKK
	N(2)=N(1)+4*(IU)*(IU)
	N(3)=N(2)+4*(IU)*(IU)


      P(1)=0.D0    
      P(2)=0.D0    
      P(3)=0.D0    
      Q(1)=0.D0    
      Q(2)=0.D0    
      Q(3)=0.D0    
      DO K=1,3
	TOTA=0.D0

c	WRITE(IFILE,*) 'TUBE = ',K
      DO L=1,NLAYER
C	WRITE(IFILE,*) 'NLAYER = ',L
C      WRITE(IFILE,*) 'PATCH = ',K
	N1=N(K)
	N2=N1+IU-1
	VSRE=0.D0
        DO I=N1,N2
	   II=I+(L-1)*NLE
c	WRITE(IFILE,*) 'ELEMENTS = ',II,NEL(1,II),PRIT(1,II)
c	WRITE(IFILE,*) 'ELEMENTS = ',II
	DA=SURF(NEL,II,CORD,NDIM,NEL(1,II),NEL(4,II),NEL(8,II),NEL(5,II))
	TOTA=TOTA+DA
	
	if (k.eq.2) ALL=DABS(CORD(1,NEL(2,II)))
	IF (KKK.EQ.0) XPIPEMIN=2.D0*ALL
	
	VSR=0.D0
	DO JJ=1,NDIM
	 VSR=VSR+DSQRT(GNODE(2,1,NEL(JJ,II))**2+GNODE(2,2,NEL(JJ,II))**2+
     & GNODE(2,3,NEL(JJ,II))**2)*(1.D0/NDIM)
	ENDDO
	
	   P(K)=P(K)+PRIT(1,II)
	   Q(K)=Q(K)+VSR*DA
c	   VSRE=VSRE+VSR
	  ENDDO
      ENDDO
	P(K)=P(K)/(1.D0*IU*NLAYER)
c	Q(K)=TOTA*VSRE/(1.D0*IU*NLAYER)
C	WRITE(IFILE,*) P(K)
	ENDDO

      TAU=0.D0
      TAU1=0.D0
	SURFTAU=0.D0
	TAUX=0.D0
	TAUY=0.D0
	TAUZ=0.D0
	IBR=0
C      do 110 nbrel=1,net
	
c      do rad=1.d-4,200.d-4,10.d-4
C	DO 110 III=0,NLAYER-1
C      DO 100 NBREL=1+NLE*III,NLE*(III+1)

      

      do 100 nbrel=1,net 	  
	  
      xc=0.d0
      yc=0.d0
	DO 11 I=1,NDIM
	  M(I)=NEL(I,NBREL)
	  NODE=M(I)
        NID(I)=0
	xc=xc+cord(1,node)
	yc=yc+cord(2,node)
	   IF (DABS(GNODE(2,1,NODE)).GT.1D-10) GOTO 11
      IF(ID(1,M(I)).EQ.0.AND.ID(2,M(I)).EQ.0.AND.ID(3,M(I)).EQ.0)
     &NID(I)=1
11    CONTINUE 
      xc=xc/(1.d0*ndim)
      yc=yc/(1.d0*ndim)
C	if (sqrt(xc*xc+yc*yc)>all) goto 100
C	if (sqrt(xc*xc+yc*yc)<0.0018d0) goto 100
      T4(1)=XC
	T4(2)=YC
	XX1=XMAX
	T1(1)=XX1
	T1(2)=0.D0

	T2(1)=-0.5D0*XX1
	T2(2)=XX1*DCOS(PI/6.D0)

	T3(1)=-0.5D0*XX1
	T3(2)=-XX1*DCOS(PI/6.D0)
      if((ISSPR(T1,T2,T3,T4)*ISSPR(T2,T3,T1,T4)*ISSPR(T1,T3,T2,T4))
     &.EQ.1) GOTO 100



	T1(1)=XX
	T1(2)=0.D0

	T2(1)=-0.5D0*XX
	T2(2)=XX*DCOS(PI/6.D0)

	T3(1)=-0.5D0*XX
	T3(2)=-XX*DCOS(PI/6.D0)
      if((ISSPR(T1,T2,T3,T4)*ISSPR(T2,T3,T1,T4)*ISSPR(T1,T3,T2,T4))
     &.EQ.0) GOTO 100

	IF (XX.LT.XMAX) GOTO 100


	T1(1)=XpipeMin
	T1(2)=0.D0

	T2(1)=-0.5D0*XpipeMin
	T2(2)=XpipeMin*DCOS(PI/6.D0)

	T3(1)=-0.5D0*XpipeMin
	T3(2)=-XpipeMin*DCOS(PI/6.D0)
	itau1=1
      if((ISSPR(T1,T2,T3,T4)*ISSPR(T2,T3,T1,T4)*ISSPR(T1,T3,T2,T4))
     &.EQ.1) itau1=0

	IF (XX.LT.XMAX) GOTO 100



      ITR(1,1)=M(8)
      ITR(1,2)=M(4)
      ITR(1,3)=M(5)
      ITR(1,4)=M(1)

      ITR(2,1)=M(7)
      ITR(2,2)=M(3)
      ITR(2,3)=M(6)
      ITR(2,4)=M(2)

      ITR(3,1)=M(6)
      ITR(3,2)=M(2)
      ITR(3,3)=M(5)
      ITR(3,4)=M(1)

      ITR(4,1)=M(7)
      ITR(4,2)=M(3)
      ITR(4,3)=M(8)
      ITR(4,4)=M(4)

      ITR(5,1)=M(3)
      ITR(5,2)=M(2)
      ITR(5,3)=M(4)
      ITR(5,4)=M(1)

      ITR(6,1)=M(7)
      ITR(6,2)=M(6)
      ITR(6,3)=M(8)
      ITR(6,4)=M(5)



      NWALL(1)=NID(1)*NID(4)*NID(5)*NID(8)
      NWALL(2)=NID(2)*NID(3)*NID(6)*NID(7)
      NWALL(3)=NID(1)*NID(2)*NID(5)*NID(6)
      NWALL(4)=NID(3)*NID(4)*NID(7)*NID(8)
      NWALL(5)=NID(1)*NID(2)*NID(3)*NID(4)
      NWALL(6)=NID(5)*NID(6)*NID(7)*NID(8)

      DO 300 NPOV=1,6
       IF (NWALL(NPOV).EQ.0) GOTO 300

	   DA=SURF(NEL,NBREL,CORD,NDIM,ITR(NPOV,1),
     &ITR(NPOV,2),ITR(NPOV,3),ITR(NPOV,4))
	   
         DO JJ=1,4
	     NN2(JJ)=ITR(NPOV,JJ) 
         ENDDO
         
	   TT1=DSQRT(PRES(1,NN2(1))**2+PRES(2,NN2(1))**2+PRES(3,NN2(1))**2)
	   TT2=DSQRT(PRES(1,NN2(2))**2+PRES(2,NN2(2))**2+PRES(3,NN2(2))**2)
	   TT3=DSQRT(PRES(1,NN2(3))**2+PRES(2,NN2(3))**2+PRES(3,NN2(3))**2)
	   TT4=DSQRT(PRES(1,NN2(4))**2+PRES(2,NN2(4))**2+PRES(3,NN2(4))**2)
C	   write(ifile,*)'ss,da',(0.25D0*(T1+T2+T3+T4)),DA
 	   TAU=TAU+0.25D0*(TT1+TT2+TT3+TT4)*DA
  	   IF (itau1.eq.1) TAU1=TAU1+0.25D0*(TT1+TT2+TT3+TT4)*DA
	   SURFTAU=SURFTAU+DA

300   CONTINUE 
100   CONTINUE
c 110   CONTINUE
c      enddo
      


      TAUX=TAUX/TOTA
      TAUY=TAUY/TOTA
      TAUZ=TAUZ/TOTA


c	WRITE(IFILE,*) 'KKK = ', KKK
c      RES1=(p(2)-p(1)*DCOS(PI/3.0)-p(3)*DCOS(PI/3.0))*tota
      RES1=(p(2)-p(1)-p(3))*tota
	RES2=TAU
	RES2X=TAUX/Q(2)
	RES2Y=TAUY/Q(2)
	RES2Z=TAUZ/Q(2)
C	WRITE(IFILE,*) '(P(2)-(P(1)+P(3)))/Q(2)= ',(P(2)-(P(1)+P(3)))/Q(2)
C	WRITE(IFILE,*) 'TAU/Q(2)= ', TAU/Q(2)

      F1X=P(1)*DCOS(PI/6.D0)
      F1Y=P(1)*DCOS(PI/3.D0)
      F3X=P(3)*DCOS(PI/6.D0)
      F3Y=P(3)*DCOS(PI/3.D0)

      RES1Y=(-P(1)*DCOS(PI/6.0)/Q(1)+P(3)*DCOS(PI/6.0)/Q(3))
      RES1X=P(2)/Q(2)-P(1)*DCOS(PI/3.0)/Q(1)-P(3)*DCOS(PI/3.0)/Q(3)
      
c      WRITE(IFILE,'(8D13.5)')0.5*XX,2.d0*RES1,2.d0*RES2,p(2),SURFTAU
      WRITE(IFILE,'(8D13.5)')0.5*XX,2.d0*RES1,2.d0*RES2,2.d0*tau1
150   continue

	CLOSE (IFILE)
c	STOP

      END
C==========================================================================
C==========================================================================
C==========================================================================
      FUNCTION ISSPR(T1,T2,T3,T4)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION T1(2),T2(2),T3(2),T4(2)
      


      TI3=(T3(2) - T1(2))*(T2(1) - T1(1)) - 
     &(T2(2) - T1(2)) * (T3(1) - T1(1))
	
	TI4=(T4(2) - T1(2))*(T2(1) - T1(1)) - 
     &(T2(2) - T1(2)) * (T4(1) - T1(1))
	  
      ISSPR=0
      IF ((TI3 * TI4).GT.0) ISSPR=1


      END
C==========================================================================
C==========================================================================
