C==========================================================================
      SUBROUTINE GENERcoating(II,NIZ)
c======================================================================
      !USE IFPORT                
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension NIZ(*)

c======================================================================
C WRITING NODES	
      

C      CALL SRAND(123)    
      CALL RANDOM_SEED()
      DO I=1,110011
       NIZ(I)=0
      ENDDO
      
      DO I = 1, 5000000      
       N = DRAND(0)*50000000+1
C       WRITE (II,*) N
       if (mod(n,5000).eq.1) NIZ ((n/5000)*10+mod(n,5000))=1
       if (mod(n,5000).eq.10) NIZ ((n/5000)*10+mod(n,5000))=2
c       if (mod(n,5001).eq.1) NIZ ((n/5001)*11+mod(n,5001))=1
c       if (mod(n,5001).eq.11) NIZ ((n/5001)*11+mod(n,5001))=1
c      if (mod(n,5001).eq.1) WRITE (II,*) N,(n/5001)*11+mod(n,5001)
c      if (mod(n,5001).eq.11) WRITE (II,*) N,(n/5001)*11+mod(n,5001)
      END DO    
c     STOP
      maxsil=0 
      DO I=1,100000
       IF (NIZ(I).NE.0) THEN
        maxsil=maxsil+1
C        WRITE (II,*) I
       ENDIF
      ENDDO
C      STOP

c======================================================================
C WRITING NODES	
      NR=10
      NL=1000*NR
      R=0.1D0
      AL=100.D0


      WRITE(II,*)'C NODES'
      
      DX=R/NR
      DY=AL/NL
      K=0
      DO J=0,NL 
	  Y=J*DY 
        DO I=0,NR
         IDT=0
	   X=I*DX
	   K=K+1
C	   WRITE(II,1000)K,1,1,1,1,NIZ(K),X,Y,0.D0
	   WRITE(II,1000)K,1,1,1,1,0,X,Y,0.D0
	  ENDDO
      ENDDO
	  


c======================================================================
C WRITING ELEMENTS
      K=0
      WRITE(II,*)'C ELEMENTS'
      DO J=0,NL-1 
	  DO I=1,NR
	   K=K+1
	   I1=I+J*(NR+1)
	   WRITE(II,'(5I10)')K,I1,I1+1,I1+NR+2,I1+NR+1
	  ENDDO
      ENDDO
	  
	
1000  FORMAT (I10,5I5,3(X,E13.5))

      WRITE(II,*)'C MAXSIL', MAXSIL
      
      K=0
      DO J=0,NL-1 
	  DO I=1,NR
	   K=K+1
	   I1=I+J*(NR+1)
         IF (NIZ(K).EQ.1) WRITE(II,2000)K,I1,I1+NR+1,1
         IF (NIZ(K).EQ.2) WRITE(II,2000)K,I1+1,I1+NR+2,1
	  ENDDO
      ENDDO
          
2000  FORMAT (4I10)
	close(ii)
	STOP
c======================================================================
      
      end
c======================================================================
c======================================================================
