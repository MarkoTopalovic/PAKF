c==============================================================================
c==============================================================================
      subroutine jactSurf10node (r,s,t,ck,h,tt21,nbrel,ndim,det1,
     &nel,sf1,sf2,sf3,nizel,cord,pj,hp,ami)
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

c       r=0.0d0
c       s=0.5d0
c       t=0.5d0
  
       e1=r
       e2=s
       e3=t
       e4=1.0d0-r-s-t

       hp(1)=e1
       hp(2)=e2
       hp(3)=e3
       hp(4)=e4


      h(5)=4.0*e1*e2
      h(6)=4.0*e2*e3
      h(7)=4.0*e1*e3
      h(8)=4.0*e1*e4
      h(9)=4.0*e2*e4
      h(10)=4.0*e3*e4

      h(1)=e1*(2.0*e1-1)
      h(2)=e2*(2.0*e2-1)
      h(3)=e3*(2.0*e3-1)
      h(4)=e4*(2.0*e4-1)
	
c      open (13,file='res.txt')
c      do iii=1,10
c       write(13,*)'h( ',iii,')=',h(iii)   
c      enddo
c      stop

      p(1,1)=4.0*e1-1.0
      p(2,2)=4.0*e2-1.0
      p(3,3)=4.0*e3-1.0

      p(1,4)=-(4.0*e4-1.0)
      p(2,4)=-(4.0*e4-1.0)
      p(3,4)=-(4.0*e4-1.0)

      p(1,5)=4.0*e2
      p(2,5)=4.0*e1

      p(2,6)=4.0*e3
      p(3,6)=4.0*e2

      p(1,7)=4.0*e3
      p(3,7)=4.0*e1


      p(1,8)=4.0*(e4-e1)
      p(2,8)=-4.0*e1
      p(3,8)=-4.0*e1

      p(1,9)=-4.0*e2
      p(2,9)=4.0*(e4-e2)
      p(3,9)=-4.0*e2


      p(1,10)=-4.0*e3
      p(2,10)=-4.0*e3
      p(3,10)=4.0*(e4-e3)


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
      if (dabs(det1)<1.d-15) then
       write (*,*) 'determinante less then zero for element',nbrel
       return
      endif
	  
  
c      xxjjj(1,1)=1.d0  
c      xxjjj(1,2)=1.d0  
c      xxjjj(1,3)=1.d0  
c      xxjjj(1,4)=1.d0  

c      Xjj(0)=1.0
c      Xjj(1)=1.0
c      Xjj(2)=1.0
c      Xjj(3)=1.0
      
c	k=4+1
c      do j=1,3
c       do i=1,4
c         xjj(k)=ck(i,j)
c	    k=k+1
c	enddo
c	enddo

c      do j=1,3
c       do i=1,4
c         xxjjj(j+1,i)=ck(i,j)
c	 enddo
c	enddo

c      call minv(xxjjj,4,det2,le,me)

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
Cc      ANZ=-1.D0
      SF1=-PRIT*ANX+AMI*(V1X*ANX+V1Y*ANY+V1Z*ANZ)
      SF2=-PRIT*ANY+AMI*(V2X*ANX+V2Y*ANY+V2Z*ANZ)
      SF3=-PRIT*ANZ+AMI*(V3X*ANX+V3Y*ANY+V3Z*ANZ)

      SF1=-SF1
      SF2=-SF2
      SF3=-SF3
     
	
       end
c==============================================================================
c==============================================================================
      subroutine surfjact10 (r,s,nizel,h,gnode,sf1,sf2,sf3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension h(*),gnode(2,6,*),pp(3)
	dimension nizel(*)

           



       e1=1.0d0-r-s
       e2=r
       e3=s

       h(1)=e1
       h(2)=e2
       h(3)=e3


c       h(1)=e1*(2.d0*e1-1.d0)
c       h(2)=e2*(2.d0*e2-1.d0)
c       h(3)=e3*(2.d0*e3-1.d0)
c       h(4)=4.d0*e2*e3
c       h(5)=4.d0*e3*e1
c       h(6)=4.d0*e1*e2

       do i=1,3
	  pp(i)=gnode(2,4,nizel(i))
	 enddo


      
      
      PRIT=DOT(H,PP,3)
      
 	
	
	SF1=0.D0
      SF2=0.D0
      SF3=-prit*(1.d0)
      return



       end
c==============================================================================
c==============================================================================
      subroutine strana_tetra(nizel,nel,nbrel,r,s,t,wt,kk,ndim)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension nel(ndim+1,*),nizm(3,4),nizel(*)
	data nizm/1,3,4,1,2,4,1,2,3,2,3,4/

       
      do ii=1,4
	 n=0
	  do j=1,3
	  do i=1,3
	    if (nel(nizm(i,ii),nbrel).eq.nizel(j)) n=n+1
	   enddo
	  enddo
	if (n.eq.3) goto 10
	 enddo

  10  continue

      if (ii.eq.1) then
	  call Integ3points(s,t,wt,kk)
	  r=0.d0
      elseif (ii.eq.2) then
	  call Integ3points(r,t,wt,kk)
	  s=0.d0
      elseif (ii.eq.3) then
	  call Integ3points(r,s,wt,kk)
	  t=0.d0
      elseif (ii.eq.4) then
	  call Integ3points(r,s,wt,kk)
	  t=1.d0-r-s
	endif

      end
c==============================================================================
c==============================================================================
      subroutine Integ3points(r,s,wt,kk)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

       if (kk.eq.1) then
	   r=0.1666666666666667d0
	   s=r
	 elseif(kk.eq.2) then
	   s=0.1666666666666667d0
	   r=0.6666666666666667d0
	 elseif(kk.eq.3) then
	   r=0.1666666666666667d0
	   s=0.6666666666666667d0
	 endif
	  
	  wt=0.333333333333333333d0



      end
c==============================================================================
c==============================================================================
C==========================================================================
      FUNCTION SURF1(nizel,CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION nizel(*),CORD(3,*)


        n1=nizel(1)
        n2=nizel(2)
        n3=nizel(3)

	  X1=CORD(1,N1) 
	  Y1=CORD(2,N1) 
	  Z1=CORD(3,N1) 
	  
	  X2=CORD(1,N2) 
	  Y2=CORD(2,N2) 
	  Z2=CORD(3,N2) 
	  
	  X3=CORD(1,N3) 
	  Y3=CORD(2,N3) 
	  Z3=CORD(3,N3) 
	  

	  A=DSQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
	  B=DSQRT((X2-X3)**2+(Y2-Y3)**2+(Z2-Z3)**2)
	  C=DSQRT((X1-X3)**2+(Y1-Y3)**2+(Z1-Z3)**2)
	  S=(A+B+C)*0.5D0
	  P1=DSQRT(DABS(S*(S-A)*(S-B)*(S-C)))

	 

        SURF1=P1

	END
C==========================================================================
C==========================================================================
