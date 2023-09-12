C===============================================================================
      SUBROUTINE CHECKJACOBIAN3D (CORD,NEL,NETIP,NET,NDIM,JACOB3D,
     &IIZLAZ)
	IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      DIMENSION CORD(3,*)
	DIMENSION NEL(NDIM+1,*)

      DIMENSION NS(4),NE(4),NODE(4)
	DIMENSION VEC(4,2,3),ANGLE(4)



      PI=4.D0*DATAN(1.D0)



	IF (NETIP.NE.3) RETURN
	IF (JACOB3D.EQ.0) RETURN

	DO 100 I=1,NET
	DO 100 II=0,1
       DO 50 J=1+II*4,4+II*4
         NODE(J)=NEL(J,I)
 
 50   CONTINUE	
         NS(1)=NODE(4) 
         NS(2)=NODE(1) 
         NS(3)=NODE(2) 
         NS(4)=NODE(3) 

         NE(1)=NODE(2)
         NE(2)=NODE(3) 
         NE(3)=NODE(4) 
         NE(4)=NODE(1) 
      DO 80 K=1,3
	   VEC(1,1,K)=CORD(K,NODE(2))-CORD(K,NODE(1))
	   VEC(1,2,K)=CORD(K,NODE(4))-CORD(K,NODE(1))
	   VEC(2,1,K)=CORD(K,NODE(3))-CORD(K,NODE(2))
	   VEC(2,2,K)=CORD(K,NODE(1))-CORD(K,NODE(2))
	   VEC(3,1,K)=CORD(K,NODE(4))-CORD(K,NODE(3))
	   VEC(3,2,K)=CORD(K,NODE(2))-CORD(K,NODE(3))
	   VEC(4,1,K)=CORD(K,NODE(1))-CORD(K,NODE(4))
	   VEC(4,2,K)=CORD(K,NODE(3))-CORD(K,NODE(4))
 80   CONTINUE	
      WRITE(IIZLAZ,*) 'ELEMENT NUMBER: ',I
      DO J=1,4
	UP=
     &VEC(J,1,1)*VEC(J,2,1)+VEC(J,1,2)*VEC(J,2,2)+VEC(J,1,3)*VEC(J,2,3)
	A=DSQRT(VEC(J,1,1)**2+VEC(J,1,2)**2+VEC(J,1,3)**2)
	B=DSQRT(VEC(J,2,1)**2+VEC(J,2,2)**2+VEC(J,2,3)**2)
	IF (DABS(A*B).LT.1.D-15) THEN
	  WRITE(*,*) 'ELEMENT NUMBER: ',I
	  WRITE(*,*)'NODE(1) = ',NODE(1)
	  WRITE(*,*)'NODE(2) = ',NODE(2)
	  WRITE(*,*)'NODE(3) = ',NODE(3)
	  WRITE(*,*)'NODE(4) = ',NODE(4)
	  WRITE(*,*) 'A: ',A
	  WRITE(*,*) 'B: ',B
	  STOP ' THE EDGE OF 3D ELEMENT IS EQUAL ZERO ??!!'
	ENDIF
      ANGLE(J)=DACOS(UP/(A*B))  
      WRITE(IIZLAZ,*)'ANGLE(',J,') = ',ANGLE(J)
	
C	IF (ANGLE(J).GT.PI) THEN
	IF (ANGLE(J).lT.0.0D0) write(*,*)' angle less then zero!!!'

	IF (ANGLE(J).GT.3.13D0) THEN
C	 DO K=1,3
C	  CORD(K,NODE(J))=0.5D0*(CORD(K,NS(J))+CORD(K,NE(J)))
C	 ENDDO
	WRITE(IIZLAZ,*) '=========================================='
	WRITE(IIZLAZ,*) 'THE CORRECTION FOR JACOBIAN 3D FINITE ELEMENT 
     &WAS DONE ON ELEMENT : ',I
	WRITE(IIZLAZ,*) '=========================================='
	ENDIF

	ENDDO
100   CONTINUE	

      STOP

	END
C===============================================================================
C===============================================================================
