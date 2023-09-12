c==============================================================================
c==============================================================================
      subroutine jactSurf5node (r,s,t,ck,h,tt21,nbrel,ndim,det1,
     &nel,sf1,sf2,sf3,nizel,cord,pj,hp,ami)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension p(4,10),pj(3,*),xjj(3,3),xj(3,3),ck(21,*),cord(3,*)
      dimension me(16),le(16)
      dimension xjja(20),Xjjj(17)
      dimension h(*),hp(*),TT21(*),NEL(NDIM+1,*)
      dimension V1(10),V2(10),V3(10),pp(10),temp(10),an(3)
      DIMENSION VSHX(3),VSHY(3),VSHZ(3),SHEAR(3),vtang(3,10)
	dimension xxjjj(4,4),nizel(*)

                 

       do i=1,ndim+1
		 h(i)=0.d0 
         do j=1,4
		   p(j,i)=0.d0 
	     enddo
	   enddo


       e1=1.0d0-r-s-t
       e2=r
       e3=s
       e4=t
	 e5=256.d0*e1*e2*e3*e4

       h(1)=e1-0.25d0*e5
       h(2)=e2-0.25d0*e5
       h(3)=e3-0.25d0*e5
       h(4)=e4-0.25d0*e5
	 h(5)=e5

	 hp(1)=e1
	 hp(2)=e2
	 hp(3)=e3
	 hp(4)=e4
      


       rr=0.25*256.d0*(s*t*e1-r*s*t)
       ss=0.25*256.d0*(r*t*e1-r*s*t)
       tt=0.25*256.d0*(r*s*e1-r*s*t)


       p(1,1)=-1.0d0-rr
       p(2,1)=-1.0d0-ss
       p(3,1)=-1.0d0-tt

       p(1,2)=1.0d0-rr
       p(2,2)=0.0d0-ss
       p(3,2)=0.0d0-tt

       p(1,3)=0.0d0-rr
       p(2,3)=1.0d0-ss
       p(3,3)=0.0d0-tt

       p(1,4)=0.0d0-rr
       p(2,4)=0.0d0-ss
       p(3,4)=1.0d0-tt

       p(1,5)=4.d0*rr
       p(2,5)=4.d0*ss
       p(3,5)=4.d0*tt


           do i=1,3
		   do j=1,3
			   xj(i,j)=0.d0
			   xjj(i,j)=0.d0
			   do kk=1,ndim
				   xj(i,j)=xj(i,j)+p(i,kk)*ck(kk,j)
				   xjj(i,j)=xjj(i,j)+p(i,kk)*ck(kk,j)
			   enddo
			enddo
		enddo


      call minv(xjj,3,det1,le,me)

	  
      if (dabs(det1)<1.d-15) then
       write (*,*) 'determinante less then zero for element',nbrel
c       return
      endif
      
      
	   
	xxjjj(1,1)=1.d0  
	xxjjj(1,2)=1.d0  
	xxjjj(1,3)=1.d0  
	xxjjj(1,4)=1.d0  
	  
          do j=1,3
            do i=1,4
	        xxjjj(j+1,i)=ck(i,j)
		   enddo
	   enddo

       call minv(xxjjj,4,det2,le,me);
c
       AJ_1=(1.0/6.0)*det2
c
       det1=AJ_1;




       do i=1,3
		   do jjj=1,ndim
			   pj(i,jjj)=0.d0
              do k=1,3
				  pj(i,jjj)=pj(i,jjj)+xjj(i,k)*p(k,jjj)
			  enddo
			enddo
		enddo




      DO I=1,NDIM
      V1(I)=TT21(I)
      V2(I)=TT21(I+NDIM)
      V3(I)=TT21(I+2*NDIM)
      TEMP(I)=TT21(I+3*NDIM+4)
      IF (I.LE.4) THEN
        PP(I)=TT21(I+3*NDIM)
      ENDIF
      ENDDO
      
C
CS    JAKOBIJAN U TACKI R,S,T
CE    JACOBIAN AT POINT R,S,T
C

      HV1=DOT(H,V1,NDIM)
      HV2=DOT(H,V2,NDIM)
      HV3=DOT(H,V3,NDIM)
  
  


c      IF (DET1.LT.1.D-15) THEN
c       WRITE(*,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
c       WRITE(IIZLAZ,*)'DETERMINANTE LESS THEN ZERO FOR ELEMENT ',NBREL
c       WRITE(IIZLAZ,*)'DETERMINANTE= ',DET1
c       WRITE(IIZLAZ,*)'NODES  COORDINATES'
c       DO I=1,NDIM
c        WRITE(IIZLAZ,1000) NEL(I,NBREL),CK(I,1),CK(I,2),CK(I,3)
c        WRITE(IIZLAZ,1001) NEL(I,NBREL),CK(I,1),CK(I,2),CK(I,3)
c	  ENDDO
c       return 
c      ENDIF
 1000 FORMAT(I5,3(D13.5))
 1001 FORMAT(I5,3(f10.6))

     
      ax=cord(1,nizel(2))-cord(1,nizel(1))
      ay=cord(2,nizel(2))-cord(2,nizel(1))
      az=cord(3,nizel(2))-cord(3,nizel(1))
      
      bx=cord(1,nizel(3))-cord(1,nizel(1))
      by=cord(2,nizel(3))-cord(2,nizel(1))
      bz=cord(3,nizel(3))-cord(3,nizel(1))


      ANX=ay*bz-az*by
      ANY=az*bx-ax*bz
      ANZ=ax*by-ay*bx



      aaa=dsqrt(anx**2+any**2+anz**2)

      ANX=ANX/AAA
      ANY=ANY/AAA
      ANZ=ANZ/AAA

      AN(1)=ANX
      AN(2)=ANY
      AN(3)=ANZ
      
      V1X=0.D0
      V1Y=0.D0
      V1Z=0.D0
      V2X=0.D0
      V2Y=0.D0
      V2Z=0.D0
      V3X=0.D0
      V3Y=0.D0
      V3Z=0.D0

      VSHX(1)=0.D0
      VSHX(2)=0.D0
      VSHX(3)=0.D0
      VSHY(1)=0.D0
      VSHY(2)=0.D0
      VSHY(3)=0.D0
      VSHZ(1)=0.D0
      VSHZ(2)=0.D0
      VSHZ(3)=0.D0
      DO I=1,NDIM
       V1X=V1X+PJ(1,I)*V1(I)
       V1Y=V1Y+PJ(2,I)*V1(I)
       V1Z=V1Z+PJ(3,I)*V1(I)
       V2X=V2X+PJ(1,I)*V2(I)
       V2Y=V2Y+PJ(2,I)*V2(I)
       V2Z=V2Z+PJ(3,I)*V2(I)
       V3X=V3X+PJ(1,I)*V3(I)
       V3Y=V3Y+PJ(2,I)*V3(I)
       V3Z=V3Z+PJ(3,I)*V3(I)
c       VSHX(1)=VSHX(1)+PJ(1,I)*VTANG(1,I)
c       VSHX(2)=VSHX(2)+PJ(2,I)*VTANG(1,I)
c       VSHX(3)=VSHX(3)+PJ(3,I)*VTANG(1,I)
c       VSHY(1)=VSHY(1)+PJ(1,I)*VTANG(2,I)
c       VSHY(2)=VSHY(2)+PJ(2,I)*VTANG(2,I)
c       VSHY(3)=VSHY(3)+PJ(3,I)*VTANG(2,I)
c       VSHZ(1)=VSHZ(1)+PJ(1,I)*VTANG(3,I)
c       VSHZ(2)=VSHZ(2)+PJ(2,I)*VTANG(3,I)
c       VSHZ(3)=VSHZ(3)+PJ(3,I)*VTANG(3,I)
      ENDDO
c      SHEAR(1)=DOT(AN,VSHX,3)
c      SHEAR(2)=DOT(AN,VSHY,3)
c      SHEAR(3)=DOT(AN,VSHZ,3)
      PRIT=DOT(HP,PP,4)
c      SF1=ANX
c      SF2=ANY
c      SF3=ANZ
c      ANZ=-1.D0
      SF1=-PRIT*ANX+AMI*(V1X*ANX+V1Y*ANY+V1Z*ANZ)
      SF2=-PRIT*ANY+AMI*(V2X*ANX+V2Y*ANY+V2Z*ANZ)
      SF3=-PRIT*ANZ+AMI*(V3X*ANX+V3Y*ANY+V3Z*ANZ)
      
	
       end
c==============================================================================
c==============================================================================
