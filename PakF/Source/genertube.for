C==========================================================================
      SUBROUTINE genertube(II)
c======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

c======================================================================
C WRITING NODES	
      
      NR=30
      NL=1200

      NR=60
      NL=2400

      WRITE(II,*)'C NODES'
      
      R=2*1.0D0
      AL=40.D0
      BIAS=0.9999D0
      
      
      A1=R*(1.D0-BIAS)/(1.D0-BIAS**(NR))
      DX=R/NR
      DY=AL/NL
      K=0
      DO J=0,NL 
	  Y=J*DY 
        DO I=0,NR
          IDX=0
          IDY=0
c         IF (I.EQ.0) IDX=1
c         IF (I.EQ.NR.OR.J.EQ.0) THEN
         IF (I.EQ.NR.OR.J.EQ.0.or.I.EQ.0) THEN
           IDX=1
           IDY=1
         ENDIF
C	   X=I*DX

         R1=R
         X1=0.D0
         
c         IF (Y.GE.20.D0.AND.Y.LE.22.D0) THEN
c           X1=DSQRT((Y-20.D0)/2.D0)
c         ELSEIF (Y.GT.22.D0)THEN
c           X1=1.D0
c         ENDIF
         
c         IF (Y.GE.20.D0.AND.I.EQ.0) THEN
c          IDX=1
c          IDY=1
c         ENDIF


         
         R1=R-X1
         A1=R1*(1.D0-BIAS)/(1.D0-BIAS**(NR))

         X=X1+A1*(1.D0-BIAS**I)/(1.D0-BIAS)
	   K=K+1
	   WRITE(II,1000)K,IDX,IDY,1,1,1,X,Y,0.D0
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
	  
c======================================================================
      WRITE(II,*)'C VELOCITIES'
         A1=R*(1.D0-BIAS)/(1.D0-BIAS**(NR))
      DO I=2,NR 
c	   X=I*DX
         X=A1*(1.D0-BIAS**(I-1))/(1.D0-BIAS)
c        WRITE(II,'(3I5,E13.5)')I,2,1,1.D0-(X/R)**2
        WRITE(II,'(3I5,E13.5)')I,2,1,1.D0-((X-R/2.d0)/(R/2.d0))**2
      ENDDO
c======================================================================
	
1000  FORMAT (I10,5I5,3(X,E13.5))

	
	close(ii)
	STOP
c======================================================================
      END
C==========================================================================
