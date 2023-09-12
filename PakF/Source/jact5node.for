c==============================================================================
c==============================================================================
      subroutine jact5node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,
     &hp,hv1,hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,nel)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension p(5,12),pj(3,*),xjj(3,3),xj(3,3),ck(21,*)
      dimension me(16),le(16)
      dimension xjja(10),XXjj(17),Jj(5,4)
      dimension h(*),hp(*),TT21(*),NEL(NDIM+1,*)
      dimension V1(4),V2(4),V3(4),pp(4),temp(4),an(3)
      DIMENSION VSHX(3),VSHY(3),VSHZ(3),SHEAR(3),vtang(3,4)
	dimension xxjjj(4,4)

           

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


      iupwin=0
      IF (IUPWIN.EQ.1.AND.INDUP.EQ.0) THEN
c       CALL INTER1(CK,V1,V2,V3,H,PJ,NDIM,AMI)
c      	AK=2.D0/AMI 
      	AK=2.D0/2.d0
       VV=sqrt(HV1**2+HV2**2+HV3**2)
	if (dabs(vv).gt.1.d-15) then
      DO I=1,NDIM
c        PP=AKK*(HV1*PJ(1,I)+HV2*PJ(2,I)+HV3*PJ(3,I))/VV
        P1=AK*(HV1*PJ(1,I)+HV2*PJ(2,I)+HV3*PJ(3,I))/VV
        Hp(I)=Hp(I)+P1
      ENDDO
	endif

       INDUP=1
c      GOTO 100
      ENDIF
      

      HXU=0.D0
      HYU=0.D0
      HZU=0.D0
      HXV=0.D0
      HYV=0.D0
      HZV=0.D0
      HXW=0.D0
      HYW=0.D0
      HZW=0.D0
      ZVXT=0.D0
      ZVYT=0.D0
      ZVZT=0.D0
      DO L =1,NDIM
        HXU=HXU+PJ(1,L)*V1(L)
        HYU=HYU+PJ(2,L)*V1(L)
        HZU=HZU+PJ(3,L)*V1(L)
        HXV=HXV+PJ(1,L)*V2(L)
        HYV=HYV+PJ(2,L)*V2(L)
        HZV=HZV+PJ(3,L)*V2(L)
        HXW=HXW+PJ(1,L)*V3(L)
        HYW=HYW+PJ(2,L)*V3(L)
        HZW=HZW+PJ(3,L)*V3(L)
        ZVXT=ZVXT+PJ(1,L)*TEMP(L)
        ZVYT=ZVYT+PJ(2,L)*TEMP(L)
        ZVZT=ZVZT+PJ(3,L)*TEMP(L)
      ENDDO




       
      
c      IF(KFIX.GT.0) GO TO 70
c      hp(1)=1.d0
      RETURN




       end
c==============================================================================
C==========================================================================
