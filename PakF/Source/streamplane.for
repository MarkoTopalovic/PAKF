C==========================================================================
      SUBROUTINE streamplane(NPT,NET,GNODE,NEL,NDIM,CORD,NEWNUM,KOR)
c======================================================================
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*)
      DIMENSION NELL(8)
	DIMENSION NEWNUM(*)

C STREAMLINES SUBROUTINE FOR PLANE
      CHARACTER *12 STRING
      CHARACTER *4 INTSTR



      II=100+KOR
	STRING='STRE'//INTSTR(KOR)//'.'//'DAT'
      OPEN(II,FILE=STRING)	
	K=0

c======================================================================
C WRITING NODES	
      K=0
      DO I=1,NPT
	  NEWNUM(I)=0
	  IF (CORD(3,I).EQ.0.D0) THEN
	    K=K+1
	    NEWNUM(I)=K
	    WRITE(II,1000)K,1,1,1,1,1,(CORD(J,I),J=1,3)
	  ENDIF
	ENDDO


c======================================================================
C WRITING ELEMENTS
      K=0
      WRITE(II,*)'C ELEMENTS'
	do NBREL=1,NET
	  NN=0
	 DO J=1,NDIM
	  NODE=NEL(J,NBREL)
	  IF (CORD(3,NODE).EQ.0.D0) THEN
	    NN=NN+1
	    NELL(NN)=NEWNUM(NODE)
	  ENDIF
	 ENDDO
	IF (NN.EQ.4) THEN
	   K=K+1
	   WRITE(II,'(5I5)')K,NELL(1),NELL(2),NELL(4),NELL(3)
	ENDIF
	enddo
c======================================================================
c======================================================================
C WRITING VELOCITES
      WRITE(II,*)'C VELOCITES'
      K=0
      DO I=1,NPT
	  IF (NEWNUM(I).NE.0) THEN
	    K=K+2
	    WRITE(II,2000)NEWNUM(I),1,1,GNODE(2,1,I)
	    WRITE(II,2000)NEWNUM(I),2,1,GNODE(2,2,I)
	  ENDIF
	ENDDO
      WRITE(II,*)'TOTAL VELOCITIES',K
c======================================================================
	
1000  FORMAT (6I5,3F10.4)
2000  FORMAT (3I5,F10.7)
	
	close(ii)
c======================================================================
      END
C==========================================================================
