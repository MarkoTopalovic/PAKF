C=======================================================================
      SUBROUTINE SHELLWRITE(NEL,NBR2,NET,NGE,II,NETIP,NDIMM,ID,CORD,NPT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine SHELLWRITE is used for printing SHELL data for finite elements 
CE to output file *.UNV
C

      DIMENSION NEL(NDIMM+1,*),ID(6,*)
	DIMENSION CORD(3,*)
	DIMENSION NIZ(4,6)
	DATA NIZ/1,2,3,4,5,6,7,8,1,5,6,2,4,3,7,8,1,4,8,5,3,2,6,7/
C
C     E L E M E N T I   3/D
C
C      IF(ISRPS.EQ.0.AND.(NBR2.NE.8.AND.NBR2.NE.20))
C     1WRITE(II,2200) NGE
C      IF(ISRPS.EQ.1.AND.(NBR2.NE.8.AND.NBR2.NE.20))
C     1WRITE(II,6200) NGE

      NBR=NBR2
      IF(NBR2.LT.20) NBR2=8
      IF(NBR2.EQ.21) NBR2=20
C     GRAFICKI OPIS ELEMENTA: SA 8 CVOROVA = 19, SA 20 CVOROVA = 20
      ITYPE=19
      IF(NBR2.EQ.20) ITYPE=20
C     VRSTA 3/D ELEMENTA: 
      IE1=115
      IF(NBR2.EQ.20) IE1=116
C     TABELA FIZICKIH OSOBINA
      IE2=1
C     TABELA MATERIJALA
      IE3=1
C     BOJA  
      ICOL=8
      IND=-1
      ITYP=71
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
      DO 10 I=1,NET
	 DO J=1,6
        N1=NEL(NIZ(1,J),I)
        N2=NEL(NIZ(2,J),I)
        N3=NEL(NIZ(3,J),I)
        N4=NEL(NIZ(4,J),I)
	  IF ( ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(3,N1).EQ.0.AND.
     &ID(1,N2).EQ.0.AND.ID(2,N2).EQ.0.AND.ID(3,N2).EQ.0.AND.
     &ID(1,N3).EQ.0.AND.ID(2,N3).EQ.0.AND.ID(3,N3).EQ.0.AND.
     &ID(1,N4).EQ.0.AND.ID(2,N4).EQ.0.AND.ID(3,N4).EQ.0) THEN
         WRITE(II,1000) I,5,94,1,1,7,4
         WRITE(II,1000) N1,N2,N3,N4
	   ENDIF
	 ENDDO


   10 CONTINUE
      WRITE(II,1100) IND

      NBR2=NBR

C ELEMENTS
C=====================================================
      TOL=1.D-5
      IND=-1
      ITYP=71
	WRITE(II,1100) IND
      WRITE(II,1100) ITYP
      DO 20 I=1,NET
	 DO J=1,6
        N1=NEL(NIZ(1,J),I)
        N2=NEL(NIZ(2,J),I)
        N3=NEL(NIZ(3,J),I)
        N4=NEL(NIZ(4,J),I)
	  IF ( ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(3,N1).EQ.0.AND.
     &ID(1,N2).EQ.0.AND.ID(2,N2).EQ.0.AND.ID(3,N2).EQ.0.AND.
     &ID(1,N3).EQ.0.AND.ID(2,N3).EQ.0.AND.ID(3,N3).EQ.0.AND.
     &ID(1,N4).EQ.0.AND.ID(2,N4).EQ.0.AND.ID(3,N4).EQ.0) THEN
        R1=DSQRT(CORD(1,N1)**2+CORD(2,N1)**2)
        R2=DSQRT(CORD(1,N2)**2+CORD(2,N2)**2)
        R3=DSQRT(CORD(1,N3)**2+CORD(2,N3)**2)
        R4=DSQRT(CORD(1,N4)**2+CORD(2,N4)**2)
        IF (DABS(R1-2.5D-4).LT.TOL.AND.DABS(R2-2.5D-4).LT.TOL.AND.
     &DABS(R3-2.5D-4).LT.TOL.AND.DABS(R4-2.5D-4).LT.TOL) THEN
          WRITE(II,1000) I,5,94,1,1,7,4
          WRITE(II,1000) N1,N2,N3,N4
	   ENDIF
	   ENDIF
	 ENDDO


   20 CONTINUE
      WRITE(II,1100) IND
C=====================================================

 1200 FORMAT(4I10,4E13.5)
 1100 FORMAT(I6)
 1000 FORMAT(8I10)


      RETURN
      END
C=====================================================
C=====================================================
C=======================================================================
      SUBROUTINE SHELL2(NEL,NET,NETIP,NDIM,ID,CORD,NPT,PRES)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine SHELLWRITE is used for printing SHELL data for finite elements 
CE to output file *.UNV
C

      DIMENSION NEL(NDIM+1,*),ID(6,*)
	DIMENSION CORD(3,*),PRES(3,*)
 
      DIMENSION N(3)
	DIMENSION NIZ(4,6)
	DATA NIZ/1,2,3,4,5,6,7,8,1,5,6,2,4,3,7,8,1,4,8,5,3,2,6,7/


      PI=4.D0*DATAN(1.D0)

      TOL=1.D-5
	IND=-1
      
C NODES
C=====================================================
      II=139
	ITYP=15
      OPEN(II,FILE='newcoord.unv')
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP
      
	
	IU=NEL(8,1)-NEL(5,1)-1
	NL=NEL(1,1)-NEL(5,1)



      NLAYER=NPT/NL

      IN=2*(IU+1)*(IU+1)-(IU+1)
	N(1)=1
	N(2)=2*IN-IU
	N(3)=4*IN-IU-2*(IU+1)

      Z=0.D0

	X1=CORD(1,NEL(8,1))
	Y1=CORD(2,NEL(8,1))

	X2=CORD(1,NEL(5,1))
	Y2=CORD(2,NEL(5,1))
   
      DL=DSQRT((X2-X1)**2+(Y2-Y1)**2)
	N1=N(1)+(NLAYER-1)*NL
	DZ=CORD(3,N1)/(NLAYER-1)





      DO L=1,NLAYER
	Y=DZ*L
	X=0.D0
C	WRITE(51,*) 'NLAYER = ',L
       DO K=1,3
C	WRITE(51,*) 'PATCH = ',K
 	 N1=N(K)+(L-1)*NL
	 N2=N1+IN-IU*K
        DO I=N1,N2,(IU+1)
	   X=X+DL
         WRITE(II,1200) I,0,0,8,X,Y,Z
	  ENDDO
	 ENDDO
	ENDDO
      WRITE(II,1100) IND
C=====================================================


C ELEMENTS
C=====================================================
      ITYP=71
	WRITE(II,1100) IND
      WRITE(II,1100) ITYP


	IE=0
      DO L=1,NLAYER-1
	Z=(1.D-4)*L
	Y=0.D0
C	WRITE(51,*) 'NLAYER = ',L
       DO K=1,3
C	WRITE(51,*) 'PATCH = ',K
 	 N1=N(K)+(L-1)*NL
	 N2=N1+IN-IU*K
        DO I=N1,N2-(IU+1),(IU+1)
	   Y=Y+1.D-4
	    IE=IE+1
          WRITE(II,1000) IE,5,94,1,1,7,4
          WRITE(II,1000) I,I+IU+1,I+IU+1+NL,I+NL
	  ENDDO
	 ENDDO


 	 N1=N(1)+(L-1)*NL
	 N2=N1+IN-IU*1-1
	 N3=N(2)+(L-1)*NL
        IE=IE+1
        WRITE(II,1000) IE,5,94,1,1,7,4
        WRITE(II,1000) N2,N3,N3+NL,N2+NL

 	 N1=N(2)+(L-1)*NL
	 N2=N1+IN-IU*2-1-1
	 N3=N(3)+(L-1)*NL
        IE=IE+1
        WRITE(II,1000) IE,5,94,1,1,7,4
        WRITE(II,1000) N2,N3,N3+NL,N2+NL




	ENDDO
      WRITE(II,1100) IND
C=== ==================================================
C
C  FOR PRINTING SHEAR STRESS
C
C NUMBER '12' IS ANALOGY NODAL ACCELERATIONS
      KOR=1

      WRITE(II,5100) IND
      WRITE(II,5100) 55
      WRITE(II,5003) 'Pillar shear stress'
      WRITE(II,7007)
      WRITE(II,2001) KOR
      WRITE(II,3006) 2,1,1,12,2,3
      WRITE(II,3006) 1,1,KOR,KOR
      WRITE(II,5202) 1.D-3
      DO L=1,NLAYER
	Z=1.D-4*L
	Y=0.D0
C	WRITE(51,*) 'NLAYER = ',L
       DO K=1,3
C	WRITE(51,*) 'PATCH = ',K
 	 N1=N(K)+(L-1)*NL
	 N2=N1+IN-IU*K
        DO I=N1,N2,(IU+1)
	   Y=Y+1.D0
	      WRITE(II,5000) I
            WRITE(II,5200) (PRES(J,I),J=1,3)
C            WRITE(II,5200) 1.D0*I,1.D0*I*L,1.D0*I*L*L
	  ENDDO
	 ENDDO
	ENDDO
      WRITE(II,5100) IND
C=====================================================



      CLOSE(II)

      RETURN






      PI=4.D0*DATAN(1.D0)

      TOL=1.D-5
	IND=-1
      
C NODES
C=====================================================
      II=139
	ITYP=15
      OPEN(II,FILE='newcoord.unv')
      WRITE(II,1100) IND
      WRITE(II,1100) ITYP

      DO I=1,NPT
	 X=CORD(1,I)
	 Y=CORD(2,I)
	 Z=CORD(3,I)

	 R=DSQRT(X**2+Y**2)
	IF (DABS(R-2.5D-4).LT.TOL) THEN
	 FI=DACOS(DABS(X)/R)
       IF (X.GE.0.D0.AND.Y.GE.0.D0) THEN
	  FI=FI
	 ELSEIF(Y.LE.0.D0.AND.X.GE.0.D0) THEN 
	  FI=2.D0*PI-FI
	 ELSEIF(X.LT.0.D0.AND.Y.GE.0.D0) THEN
	  FI=PI+FI
	 ELSEIF(X.LT.0.D0.AND.Y.LT.0.D0) THEN
	  FI=PI-FI
	 ENDIF
	 X=100.D-4
	 Y=R*FI
       WRITE(II,1200) I,0,0,8,X,Y,Z
	ENDIF

	ENDDO

      WRITE(II,1100) IND
C=====================================================


C ELEMENTS
C=====================================================
      ITYP=71
	WRITE(II,1100) IND
      WRITE(II,1100) ITYP
      DO 20 I=1,NET
	 DO J=1,6
        N1=NEL(NIZ(1,J),I)
        N2=NEL(NIZ(2,J),I)
        N3=NEL(NIZ(3,J),I)
        N4=NEL(NIZ(4,J),I)
	  IF ( ID(1,N1).EQ.0.AND.ID(2,N1).EQ.0.AND.ID(3,N1).EQ.0.AND.
     &ID(1,N2).EQ.0.AND.ID(2,N2).EQ.0.AND.ID(3,N2).EQ.0.AND.
     &ID(1,N3).EQ.0.AND.ID(2,N3).EQ.0.AND.ID(3,N3).EQ.0.AND.
     &ID(1,N4).EQ.0.AND.ID(2,N4).EQ.0.AND.ID(3,N4).EQ.0) THEN
        R1=DSQRT(CORD(1,N1)**2+CORD(2,N1)**2)
        R2=DSQRT(CORD(1,N2)**2+CORD(2,N2)**2)
        R3=DSQRT(CORD(1,N3)**2+CORD(2,N3)**2)
        R4=DSQRT(CORD(1,N4)**2+CORD(2,N4)**2)
        IF (DABS(R1-2.5D-4).LT.TOL.AND.DABS(R2-2.5D-4).LT.TOL.AND.
     &DABS(R3-2.5D-4).LT.TOL.AND.DABS(R4-2.5D-4).LT.TOL) THEN
          WRITE(II,1000) I,5,94,1,1,7,4
          WRITE(II,1000) N1,N2,N3,N4
	   ENDIF
	   ENDIF
	 ENDDO


   20 CONTINUE
      WRITE(II,1100) IND
C=====================================================
C
C  FOR PRINTING SHEAR STRESS
C
C NUMBER '12' IS ANALOGY NODAL ACCELERATIONS
      KOR=1

      WRITE(II,5100) IND
      WRITE(II,5100) 55
      WRITE(II,5003) NASLOV
      WRITE(II,7007)
      WRITE(II,2001) KOR
      WRITE(II,3006) 2,1,1,12,2,3
      WRITE(II,3006) 1,1,KOR,KOR
      WRITE(II,5202) 1.D-3
      DO I=1,NPT
	 X=CORD(1,I)
	 Y=CORD(2,I)
	 Z=CORD(3,I)

	 R=DSQRT(X**2+Y**2)
	 IF (DABS(R-2.5D-4).LT.TOL) THEN
	      WRITE(II,5000) I
            WRITE(II,5200) (PRES(J,I),J=1,3)
	 ENDIF
	ENDDO
      WRITE(II,5100) IND
C=====================================================



      CLOSE(II)
C
 1200 FORMAT(4I10,4E13.5)
 1100 FORMAT(I6)
 1000 FORMAT(8I10)

 5100 FORMAT(I6)
 5003 FORMAT(A80)
 5000 FORMAT(6I10)
 5200 FORMAT(6(1PE13.5))
 5201 FORMAT(6(1PE20.12))
 5202 FORMAT(1PE13.2)
C-----------------------------------------------------------------------
 2005 FORMAT('CVORNE TRANSLACIJE I ROTACIJE')
 2006 FORMAT('CVORNE BRZINE I UGAONE BRZINE')
 2007 FORMAT('CVORNA UBRZANJA I UGAONA UBRZANJA')
 3005 FORMAT('PRITISCI U FLUIDU')
 3006 FORMAT(6I10)
 2000 FORMAT('DATUM I VREME'/
     1       'PRAZNA'/
     1       'SLUCAJ OPTERECENJA:',I10)
 2001 FORMAT(' DATE'/
     1       ' EMPTY'/
     1       ' LOAD CASE',I10)
C-----------------------------------------------------------------------
 6005 FORMAT('NODAL TRANSLATIONS AND ROTATIONS')
 6006 FORMAT('NODAL VELOCITIES AND ANGLE VELOCITIES')
 6007 FORMAT('NODAL ACCELERATIONS AND ANGLE ACCELERATIONS')
 7006 FORMAT(' NODAL TEMPERATURE')
 7007 FORMAT(' SHEAR STRESSES')
 7008 FORMAT(' PRESSURES AT NODES')
 6000 FORMAT('DATE AND TIME'/
     1       'EMPTY'/
     1       'LOAD CASE         :',I10)
C-----------------------------------------------------------------------
      END
C=======================================================================
C======================================================================
