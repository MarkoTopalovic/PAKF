C===============================================================================
      SUBROUTINE INTERMOVING(GNODE,ID,NEL,NET,NPT,NDIM,NEQ,NWK,MHT,MAXA,
     &CORD,VVREME,IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	DIMENSION GNODE(2,7,*),CORD(3,*)
	DIMENSION ID(7,*),NEL(NDIM+1,*),MHT(*),MAXA(*)

      PI=4.D0*DATAN(1.D0)

C      ICHANGE=0

      T=100.D0
      W=2.D0*PI/T
      CX=VVREME+30.D0
      CY=15.D0*DSIN(W*VVREME)



      DO 100 NODE=102,NPT-100
      X=CORD(1,NODE)
      Y=CORD(2,NODE)
	IF (X.LT.1.D-5) GOTO 100

	IF (ID(1,NODE).EQ.0) THEN
C	 GNODE(1,1,NODE)=0.5*(GNODE(1,1,NODE-101)+GNODE(1,1,NODE+101))
	 WRITE(IIZLAZ,*)'STARI POLOZAJ= ',NODE
	ENDIF
	IF (ID(2,NODE).EQ.0) THEN
C 	 GNODE(1,2,NODE)=0.5*(GNODE(1,2,NODE-101)+GNODE(1,2,NODE+101))
	ENDIF

	DIST=DSQRT((CX-X)**2+(CY-Y)**2)
	   ID(1,NODE)=-1	    
         ID(2,NODE)=-1  
	IF (DIST.LT.10.D0) THEN
	   WRITE(IIZLAZ,*) 'MOVING BODY = ',NODE
	   ID(1,NODE)=0	    
         ID(2,NODE)=0	  
	ENDIF

100   CONTINUE

		 
C ID MATRIX FROM THE GLOBAL CONSTANT BOUNDARY CONDITIONS

C      DO NODE=1,101
C	 ID (1,NODE)=1
C	 ID (2,NODE)=1
C	ENDDO

C      DO NODE=NPT-100,NPT
C	 ID (1,NODE)=1
C	 ID (2,NODE)=1
C	ENDDO

C      DO NODE=1,NPT,101
C	 ID (1,NODE)=1
C	 ID (2,NODE)=1
C	ENDDO


      KK=0
      DO 30 N=1,NPT
      DO JJ=1,2
       IF (ID(JJ,N).NE.0) THEN
C       IF (ID(JJ,N).EQ.0) THEN
        KK=KK+1
        ID(JJ,N)=KK
       ELSE
	  GNODE(2,JJ,N)=0.D0
	  GNODE(1,JJ,N)=0.D0
        ID(JJ,N)=0
       ENDIF
      ENDDO
30    CONTINUE


	NEQ=KK
	DO I=1,NEQ
	  MHT(I)=0
	ENDDO
      CALL MAXATF(MAXA,MHT,ID,NEL,NET,NDIM,NEQ,NWK,7,NDIM+1,iizlaz)
	
1000  FORMAT (3I5)

	END
C===============================================================================
C===============================================================================
