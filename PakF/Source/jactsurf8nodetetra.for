c==============================================================================
      subroutine jactSurf8nodeTetra (r,s,t,ck,h,tt21,nbrel,ndim,det1,
     &nel,sf1,sf2,sf3,nizel,cord,pj,hp,ami,ndimp)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension p(3,10),pj(3,*),xjj(3,3),xj(3,3),ck(21,*),cord(3,*)
      dimension me(16),le(16)
      dimension xjja(20),Xjjj(17)
      dimension h(*),hp(*),TT21(*),NEL(NDIM+1,*)
      dimension V1(10),V2(10),V3(10),pp(10),temp(10),an(3)
      DIMENSION VSHX(3),VSHY(3),VSHZ(3),SHEAR(3),vtang(3,10)
	dimension xxjjj(4,4),nizel(*)

                 

       do i=1,ndim
		 h(i)=0.d0 
         do j=1,3
		   p(j,i)=0.d0 
	     enddo
	   enddo

c       t=1.d0/3.d0
c       s=1.d0/3.d0
c       r=1.d0/3.d0


       e1=r
       e2=s
       e3=t
       e4=1.0d0-r-s-t

	 e5=27.d0*e1*e2*e3
	 e6=27.d0*e1*e2*e4
	 e7=27.d0*e1*e3*e4
	 e8=27.d0*e2*e3*e4

       hp(1)=e1
       hp(2)=e2
       hp(3)=e3
       hp(4)=e4
	if (ndimp.eq.1) hp(1)=1.d0

      
	c=-1.d0/3.d0

      h(1)=e1+c*(e5+e6+e7)
      h(2)=e2+c*(e5+e6+e8)
      h(3)=e3+c*(e5+e7+e8)
      h(4)=e4+c*(e6+e7+e8)
      h(5)=e5
      h(6)=e6
      h(7)=e7
      h(8)=e8



      if (ndimp.eq.8) then
	 do i=1,ndimp
	   hp(i)=h(i)
	 enddo
	endif 

      p(1,5)=27.d0*e2*e3
      p(2,5)=27.d0*e1*e3
      p(3,5)=27.d0*e1*e2

      p(1,6)=27.d0*(e2*e4-e1*e2)
      p(2,6)=27.d0*(e1*e4-e1*e2)
      p(3,6)=27.d0*(-e1*e2)

      p(1,7)=27.d0*(e3*e4-e1*e3)
      p(2,7)=27.d0*(-e1*e3)
      p(3,7)=27.d0*(e1*e4-e1*e3)


      p(1,8)=27.d0*(-e2*e3)
      p(2,8)=27.d0*(e3*e4-e2*e3)
      p(3,8)=27.d0*(e2*e4-e2*e3)

      

      p(1,1)=1.d0+c*(p(1,5)+p(1,6)+p(1,7))
      p(2,1)=0.d0+c*(p(2,5)+p(2,6)+p(2,7))
      p(3,1)=0.d0+c*(p(3,5)+p(3,6)+p(3,7))

      p(1,2)=0.d0+c*(p(1,5)+p(1,6)+p(1,8))
      p(2,2)=1.d0+c*(p(2,5)+p(2,6)+p(2,8))
      p(3,2)=0.d0+c*(p(3,5)+p(3,6)+p(3,8))


      p(1,3)=0.d0+c*(p(1,5)+p(1,7)+p(1,8))
      p(2,3)=0.d0+c*(p(2,5)+p(2,7)+p(2,8))
      p(3,3)=1.d0+c*(p(3,5)+p(3,7)+p(3,8))

      p(1,4)=-1.d0+c*(p(1,6)+p(1,7)+p(1,8))
      p(2,4)=-1.d0+c*(p(2,6)+p(2,7)+p(2,8))
      p(3,4)=-1.d0+c*(p(3,6)+p(3,7)+p(3,8))

c testiranje tezista
c      tez1=0.d0
c	tez2=0.d0
c	do i=1,4
c	 tez1=tez1+hp(i)*ck(i,1)
c	enddo
c	do i=1,8
c	 tez2=tez2+h(i)*ck(i,1)
c	enddo

c      write(iizlaz,*)'tez1 =', tez1
c      write(iizlaz,*)'tez2 =', tez2
c      stop


c      k=1
      do i=1,3
      do j=1,3
       xj(i,j)=0.d0
       xjj(i,j)=0.d0
	   do kk=1,ndim
           xj(i,j)=xj(i,j)+p(i,kk)*ck(kk,j)
           xjj(i,j)=xjj(i,j)+p(i,kk)*ck(kk,j)
	   enddo
c		   xjja(k)=xjj(i,j)
c	    k=k+1
	enddo
	enddo
		   
     
c========================================================== 
c      call minv(xjja,3,det1,le,me);
      call minv(xjj,3,det1,le,me);
c========================================================== 
c      k=1
c       do i=1,3
c        do j=1,3
c	    xjj(i,j)=xjja(k)
c	     k=k+1
c	  enddo
c	enddo
      if (dabs(det1)<1.d-18) then
       write (*,*) 'determinante less then zero for element',nbrel
       return
      endif
	  
  
c      xxjjj(1,1)=1.d0  
c      xxjjj(1,2)=1.d0  
c      xxjjj(1,3)=1.d0  
c      xxjjj(1,4)=1.d0  


c      do j=1,3
c       do i=1,4
c         xxjjj(j+1,i)=ck(i,j)
c	 enddo
c	enddo

c      call minv(xxjjj,4,det2,le,me)
c      call minv(xjj,4,det2,le,me)

      aJ_1=(1.0d0/6.0d0)*det1

      det1=aJ_1



c      if (det1<0.d0) stop

      do i=1,3
       do ijj=1,ndim
         pj(i,ijj)=0.d0
	    do k=1,3
            pj(i,ijj)=pj(i,ijj) + xjj(i,k)*p(k,ijj)
	    enddo
	 enddo
	enddo






      DO I=1,NDIM
      V1(I)=TT21(I)
      V2(I)=TT21(I+NDIM)
      V3(I)=TT21(I+2*NDIM)
      TEMP(I)=TT21(I+3*NDIM+4)
      IF (I.LE.ndimp) THEN
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
      PRIT=DOT(HP,PP,ndimp)
c      SF1=ANX
c      SF2=ANY
c      SF3=ANZ
c      ANZ=-1.D0
      SF1=-PRIT*ANX+AMI*(V1X*ANX+V1Y*ANY+V1Z*ANZ)
      SF2=-PRIT*ANY+AMI*(V2X*ANX+V2Y*ANY+V2Z*ANZ)
      SF3=-PRIT*ANZ+AMI*(V3X*ANX+V3Y*ANY+V3Z*ANZ)
      
      SF1=-sf1
      SF2=-sf2
      SF3=-sf3

       end
c==============================================================================
