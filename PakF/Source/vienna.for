C===============================================================================
      SUBROUTINE VIENNA(GNODE,ID,NEL,NET,NPT,NDIM,NEQ,NWK,MHT,MAXA,
     &IIZLAZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	DIMENSION GNODE(2,6,*)
	DIMENSION ID(6,*),NEL(NDIM+1,*),MHT(*),MAXA(*)

      

C      ICHANGE=0

	NSLOJ=NEL(4,1)-NEL(1,1)-1


	DO 10 I=1,NET,NSLOJ
	 IF (NEL(NDIM+1,I).EQ.1) THEN
	   GVY=0.D0
	   DO 20 II=I,I+NSLOJ-2
	   IF (NEL(NDIM+1,II).EQ.1) THEN
	    DO J=1,NDIM
	     NODE=NEL(J,II)
	     VY=GNODE(2,2,NODE) 
	     GVY=GVY+0.25D0*VY
	    ENDDO 
	   ELSE IF (NEL(NDIM+1,II).EQ.2) THEN
          IF (GVY.LT.0.D0) THEN
	     DO J=1,NDIM
	      NODE=NEL(J,II)
            ID(1,NODE)=0	    
            ID(2,NODE)=0	    
	     ENDDO
	    ELSE 
	     DO J=1,NDIM
	      NODE=NEL(J,II)
            ID(1,NODE)=1	    
            ID(2,NODE)=1	    
	     ENDDO
	    ENDIF
	   ENDIF
20    CONTINUE       
	 ENDIF
10    CONTINUE


c	GOTO 50

C ADDITIONAL ALGORITHM FOR THE CASE WHEN VALVE IS CLOSED 100%
C THIS DOES NOT WORK

	NNUM2=0

      DO 40 I=1,NET
	 IF (NEL(NDIM+1,I).EQ.1) GOTO 50
	 IF (NEL(NDIM+1,I).EQ.2) NNUM2=1 
40    CONTINUE

      IF (NNUM2.EQ.1) THEN
	      GVY=0.D0
	      DO J=1,NDIM
	        NODE=NEL(J,NET-NSLOJ+1)
	        VY=GNODE(2,2,NODE) 
	        GVY=GVY+0.25D0*VY
	      ENDDO 
	DO I=1,NET
	
      IF (NEL(NDIM+1,I).EQ.2) THEN 
          IF (GVY.LT.0.D0) THEN
	    
		 DO J=1,NDIM,3
	      NODE=NEL(J,I)
c            ID(1,NODE)=0	    
            ID(2,NODE)=0	    
	     ENDDO
	    ELSE 
	     DO J=1,NDIM,3
	      NODE=NEL(J,I)
c            ID(1,NODE)=1	    
            ID(2,NODE)=1
c		  write(iizlaz,*)'node= ',node, ID(2,NODE)
	     ENDDO
	    ENDIF
	ENDIF
	ENDDO
	ENDIF



50	KK=0
      DO 30 N=1,NPT
      DO JJ=1,2
       IF (ID(JJ,N).NE.0) THEN
        KK=KK+1
        ID(JJ,N)=KK
       ELSE
	  GNODE(2,JJ,N)=0.D0
	  GNODE(1,JJ,N)=0.D0
        ID(JJ,N)=0
       ENDIF
      ENDDO
30    CONTINUE

C      IF (GVY.LT.0.D0) THEN
C	DO I=1,NPT
C	  WRITE(IIZLAZ,1000) I,ID(1,I),ID(2,I)
C	ENDDO
c	STOP
C      endif
	 
	NEQ=KK
	DO I=1,NEQ
	  MHT(I)=0
	ENDDO
      CALL MAXATF(MAXA,MHT,ID,NEL,NET,NDIM,NEQ,NWK,6,NDIM+1,iizlaz)
	
1000  FORMAT (3I5)

	END
C===============================================================================
C===============================================================================
