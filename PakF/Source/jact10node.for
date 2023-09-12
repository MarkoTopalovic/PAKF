c==============================================================================
      subroutine jact10node1 (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,
     &hp,hv1,hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,nel)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension p(3,10),pj(3,*),xjj(3,3),xj(3,3),ck(21,*)
      dimension me(16),le(16)
      dimension xjja(20),Xjjj(17)
      dimension h(*),hp(*),TT21(*),NEL(NDIM+1,*)
      dimension V1(10),V2(10),V3(10),pp(10),temp(10),an(3)
      DIMENSION VSHX(3),VSHY(3),VSHZ(3),SHEAR(3),vtang(3,10)
	dimension xxjjj(4,4)

                 

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
      RETURN



C
CS     DETERMINATA POVRSINSKOG JAKOBIJANA
CE     SURFACE JACOBIAN DETERMINANT
C
  70   GO TO (71,72,73),KFIX
CS     KONSTANTNO KSI
CE     CONSTANT KSI
   71 DET=(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))**2+(XJ(3,1)*XJ(2,3)-
     1XJ(3,3)*XJ(2,1))**2+(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))**2
      ANX=R*(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))
      ANY=R*(XJ(3,1)*XJ(2,3)-XJ(3,3)*XJ(2,1))
      ANZ=R*(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ETA
CE     CONSTANT ETA
   72 DET=(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))**2+(XJ(1,1)*XJ(3,3)-
     1XJ(1,3)*XJ(3,1))**2+(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))**2
      ANX=-S*(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))
      ANY=S*(XJ(1,1)*XJ(3,3)-XJ(1,3)*XJ(3,1))
      ANZ=-S*(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ZETA
CE     CONSTANT ZETA
   73 DET=(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))**2+(XJ(1,1)*
     1XJ(2,3)-XJ(1,3)*XJ(2,1))**2+(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))**2
      ANX=T*(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))
      ANY=-T*(XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1))
      ANZ=T*(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))
   74 DET=DSQRT(DET)
      ANX=ANX/DET
      ANY=ANY/DET
      ANZ=ANZ/DET
      AN(1)=ANX
      AN(2)=ANY
      AN(3)=ANZ
C      WRITE(IIZLAZ,*)'ANX=',ANX
C      WRITE(IIZLAZ,*)'ANY=',ANY
C      WRITE(IIZLAZ,*)'ANZ=',ANZ
c      CALL SHEARS(AN,TT21,VTANG,NDIM)
      
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
       VSHX(1)=VSHX(1)+PJ(1,I)*VTANG(1,I)
       VSHX(2)=VSHX(2)+PJ(2,I)*VTANG(1,I)
       VSHX(3)=VSHX(3)+PJ(3,I)*VTANG(1,I)
       VSHY(1)=VSHY(1)+PJ(1,I)*VTANG(2,I)
       VSHY(2)=VSHY(2)+PJ(2,I)*VTANG(2,I)
       VSHY(3)=VSHY(3)+PJ(3,I)*VTANG(2,I)
       VSHZ(1)=VSHZ(1)+PJ(1,I)*VTANG(3,I)
       VSHZ(2)=VSHZ(2)+PJ(2,I)*VTANG(3,I)
       VSHZ(3)=VSHZ(3)+PJ(3,I)*VTANG(3,I)
      ENDDO
      SHEAR(1)=DOT(AN,VSHX,3)
      SHEAR(2)=DOT(AN,VSHY,3)
      SHEAR(3)=DOT(AN,VSHZ,3)
      PRIT=DOT(HP,PP,8)
      SF1=ANX
      SF2=ANY
      SF3=ANZ
C      SF1=-PRIT*ANX+AMI*(V1X*ANX+V1Y*ANY+V1Z*ANZ)
C      SF2=-PRIT*ANY+AMI*(V2X*ANX+V2Y*ANY+V2Z*ANZ)
C      SF3=-PRIT*ANZ+AMI*(V3X*ANX+V3Y*ANY+V3Z*ANZ)
      IF ( DET.GT.1.D-15) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(3,2000) NBREL,KFIX,R,S,T,DET
      IF(ISRPS.EQ.1)
     *WRITE(3,6000) NBREL,KFIX,R,S,T,DET
      STOP
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
     1       ' FOR ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
C


       end
c==============================================================================
c==============================================================================
      subroutine jact10node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,
     &hp,hv1,hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,nel)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension p(3,10),pj(3,*),xjj(3,3),xj(3,3),ck(21,*)
      dimension me(16),le(16)
      dimension xjja(20),Xjjj(17)
      dimension h(*),hp(*),TT21(*),NEL(NDIM+1,*)
      dimension V1(10),V2(10),V3(10),pp(10),temp(10),an(3)
      DIMENSION VSHX(3),VSHY(3),VSHZ(3),SHEAR(3),vtang(3,10)
	dimension xxjjj(4,4)

                 

       do i=1,ndim
		 h(i)=0.d0 
         do j=1,3
		   p(j,i)=0.d0 
	     enddo
	   enddo



       e1=1.0d0-r-s-t
       e2=r
       e3=s
       e4=t

       hp(1)=e1
       hp(2)=e2
       hp(3)=e3
       hp(4)=e4


      h(5)=4.0*r*(1.0-r-s-t)
      h(6)=4.0*r*s*(1.0-t)
      h(7)=4.0*s*(1.0-r-s-t)
      h(8)=4.0*r*t*(1.0-s)
      h(9)=4.0*s*t*(1.0-r)
      h(10)=4.0*t*(1.0-r-s-t)

      h(1)=1.0-r-s-t-0.5*h(5)-0.5*h(7)-0.5*h(10)
      h(2)=r-0.5*h(5)-0.5*h(6)-0.5*h(8)
      h(3)=s-0.5*h(6)-0.5*h(7)-0.5*h(9)
      h(4)=t-0.5*h(8)-0.5*h(9)-0.5*h(10)
	


      p(1,5)=4.0*(1.0-2.0*r-s-t)
      p(2,5)=-4.0*r 
      p(3,5)=-4.0*r 


      p(1,6)=4.0*s*(1.0-t) 
      p(2,6)=4.0*r*(1.0-t) 
      p(3,6)=-4.0*r*s 
      p(1,7)=-4.0*s 
      p(2,7)=4.0*(1.0-r-2.0*s-t) 
      p(3,7)=-4.0*s 


      p(1,8)=4.0*t*(1.0-s) 
      p(2,8)=-4.0*r*t 
      p(3,8)=4.0*r*(1.0-s) 


      p(1,9)=-4.0*s*t 
      p(2,9)=4.0*t*(1.0-r) 
      p(3,9)=4.0*s*(1.0-r) 



      p(1,10)=-4.0*t 
      p(2,10)=-4.0*t 
      p(3,10)=4.0*(1.0-r-s-2.0*t) 

      p(1,1)=-1.0-0.5*(p(1,5)+p(1,7)+p(1,10)) 
      p(2,1)=-1.0-0.5*(p(2,5)+p(2,7)+p(2,10)) 
      p(3,1)=-1.0-0.5*(p(3,5)+p(3,7)+p(3,10)) 

      p(1,2)=1.0-0.5*(p(1,5)+p(1,6)+p(1,8)) 
      p(2,2)=0.0-0.5*(p(2,5)+p(2,6)+p(2,8)) 
      p(3,2)=0.0-0.5*(p(3,5)+p(3,6)+p(3,8)) 

      p(1,3)=0.0-0.5*(p(1,6)+p(1,7)+p(1,9)) 
      p(2,3)=1.0-0.5*(p(2,6)+p(2,7)+p(2,9)) 
      p(3,3)=0.0-0.5*(p(3,6)+p(3,7)+p(3,9)) 

      p(1,4)=0.0-0.5*(p(1,8)+p(1,9)+p(1,10)) 
      p(2,4)=0.0-0.5*(p(2,8)+p(2,9)+p(2,10)) 
      p(3,4)=1.0-0.5*(p(3,8)+p(3,9)+p(3,10)) 

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
	  
  
      xxjjj(1,1)=1.d0  
      xxjjj(1,2)=1.d0  
      xxjjj(1,3)=1.d0  
      xxjjj(1,4)=1.d0  

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

      do j=1,3
       do i=1,4
         xxjjj(j+1,i)=ck(i,j)
	 enddo
	enddo

      call minv(xxjjj,4,det2,le,me)
c      call minv(xjj,4,det2,le,me)

      aJ_1=(1.0d0/6.0d0)*det2

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
      RETURN



C
CS     DETERMINATA POVRSINSKOG JAKOBIJANA
CE     SURFACE JACOBIAN DETERMINANT
C
  70   GO TO (71,72,73),KFIX
CS     KONSTANTNO KSI
CE     CONSTANT KSI
   71 DET=(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))**2+(XJ(3,1)*XJ(2,3)-
     1XJ(3,3)*XJ(2,1))**2+(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))**2
      ANX=R*(XJ(2,2)*XJ(3,3)-XJ(2,3)*XJ(3,2))
      ANY=R*(XJ(3,1)*XJ(2,3)-XJ(3,3)*XJ(2,1))
      ANZ=R*(XJ(2,1)*XJ(3,2)-XJ(2,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ETA
CE     CONSTANT ETA
   72 DET=(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))**2+(XJ(1,1)*XJ(3,3)-
     1XJ(1,3)*XJ(3,1))**2+(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))**2
      ANX=-S*(XJ(1,2)*XJ(3,3)-XJ(1,3)*XJ(3,2))
      ANY=S*(XJ(1,1)*XJ(3,3)-XJ(1,3)*XJ(3,1))
      ANZ=-S*(XJ(1,1)*XJ(3,2)-XJ(1,2)*XJ(3,1))
      GO TO 74
CS     KONSTANTNO ZETA
CE     CONSTANT ZETA
   73 DET=(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))**2+(XJ(1,1)*
     1XJ(2,3)-XJ(1,3)*XJ(2,1))**2+(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))**2
      ANX=T*(XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2))
      ANY=-T*(XJ(1,1)*XJ(2,3)-XJ(1,3)*XJ(2,1))
      ANZ=T*(XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1))
   74 DET=DSQRT(DET)
      ANX=ANX/DET
      ANY=ANY/DET
      ANZ=ANZ/DET
      AN(1)=ANX
      AN(2)=ANY
      AN(3)=ANZ
C      WRITE(IIZLAZ,*)'ANX=',ANX
C      WRITE(IIZLAZ,*)'ANY=',ANY
C      WRITE(IIZLAZ,*)'ANZ=',ANZ
c      CALL SHEARS(AN,TT21,VTANG,NDIM)
      
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
       VSHX(1)=VSHX(1)+PJ(1,I)*VTANG(1,I)
       VSHX(2)=VSHX(2)+PJ(2,I)*VTANG(1,I)
       VSHX(3)=VSHX(3)+PJ(3,I)*VTANG(1,I)
       VSHY(1)=VSHY(1)+PJ(1,I)*VTANG(2,I)
       VSHY(2)=VSHY(2)+PJ(2,I)*VTANG(2,I)
       VSHY(3)=VSHY(3)+PJ(3,I)*VTANG(2,I)
       VSHZ(1)=VSHZ(1)+PJ(1,I)*VTANG(3,I)
       VSHZ(2)=VSHZ(2)+PJ(2,I)*VTANG(3,I)
       VSHZ(3)=VSHZ(3)+PJ(3,I)*VTANG(3,I)
      ENDDO
      SHEAR(1)=DOT(AN,VSHX,3)
      SHEAR(2)=DOT(AN,VSHY,3)
      SHEAR(3)=DOT(AN,VSHZ,3)
      PRIT=DOT(HP,PP,8)
      SF1=ANX
      SF2=ANY
      SF3=ANZ
C      SF1=-PRIT*ANX+AMI*(V1X*ANX+V1Y*ANY+V1Z*ANZ)
C      SF2=-PRIT*ANY+AMI*(V2X*ANX+V2Y*ANY+V2Z*ANZ)
C      SF3=-PRIT*ANZ+AMI*(V3X*ANX+V3Y*ANY+V3Z*ANZ)
      IF ( DET.GT.1.D-15) RETURN
      IF(ISRPS.EQ.0)
     *WRITE(3,2000) NBREL,KFIX,R,S,T,DET
      IF(ISRPS.EQ.1)
     *WRITE(3,6000) NBREL,KFIX,R,S,T,DET
      STOP
C
 2000 FORMAT(' ** GRESKA **: JAKOBIJAN JEDNAK ILI MANJI OD NULE',
     1       ' ZA ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
 6000 FORMAT(' ** ERROR **: JACOBIAN EQUAL OR LESS THEN ZERO',
     1       ' FOR ELEMENT No.',I5/
     1       9X,'KFIX=',I5/
     2       12X,'R=',F10.5/
     3       12X,'S=',F10.5/
     4       12X,'T=',F10.5/
     5       10X,'DET=',F10.5)
C


       end
c==============================================================================
c==============================================================================
      subroutine Integ10points(r,s,t,wt,kk)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

	if (kk.eq.1) then
	wt=0.013155555555555
	r=0.20d0
	s=0.25d0
	t=0.25d0
	return
	elseif (kk.eq.2) then
	r=0.0714285714285714
	s=0.0714285714285714
	t=0.785714285714286
	wt=0.007622222222222
	return
	elseif (kk.eq.3) then
	r=0.0714285714285714
	t=0.0714285714285714
	s=0.785714285714286
	wt=0.007622222222222
	return
	else if (kk.eq.4) then
	s=0.0714285714285714
	t=0.0714285714285714
	r=0.785714285714286
	wt=0.007622222222222
	return
	elseif (kk.eq.5) then
	r=0.0714285714285714
	s=0.0714285714285714
	t=0.0714285714285714
	wt=0.007622222222222
	return
	elseif (kk.eq.6) then
	r=0.399403576166799
	s=0.399403576166799
	t=0.100596423833201
	wt=0.024888888888888
	return
	elseif (kk.eq.7) then
	r=0.399403576166799
	t=0.399403576166799
	s=0.100596423833201
	wt=0.024888888888888
	return
	elseif (kk.eq.8) then
	s=0.399403576166799
	t=0.399403576166799
	r=0.100596423833201
	wt=0.024888888888888
	return
	elseif (kk.eq.9) then
	t=0.399403576166799
	r=0.100596423833201
	s=0.100596423833201
	wt=0.024888888888888
	return
	elseif (kk.eq.10) then
	s=0.399403576166799
	r=0.100596423833201
	t=0.100596423833201
	wt=0.024888888888888
	return
	elseif (kk.eq.11) then
	r=0.399403576166799
	s=0.100596423833201
	t=0.100596423833201
	wt=0.024888888888888
	return
	endif
      
	end
c==============================================================================
c==============================================================================
c==============================================================================
      subroutine Integ4points(r,s,t,wt,kk)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

c	alfa=0.5854102966249685
c	beta=0.138196601125015

	alfa=0.58541020d0
	beta=0.13819660d0
     
	wt=0.25d0

      if (kk.eq.1) then 
	  r=beta
	  s=beta
	  t=beta
	elseif (kk.eq.2) then
	  s=beta
	  t=beta
	  r=alfa
	elseif (kk.eq.3) then
	  r=beta
	  t=beta
	  s=alfa
	elseif (kk.eq.4) then 
	  r=beta
	  s=beta
	  t=alfa
	endif


	end
c==============================================================================
c==============================================================================
c==============================================================================
      subroutine Integ5points(r,s,t,wt,kk)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      alfa=1.d0/4.d0
	beta=1.d0/6.d0
	gama=1.d0/2.d0

      if (kk.eq.1) then 
	  r=alfa
	  s=alfa
	  t=alfa
	wt=-4.d0/5.d0
	elseif (kk.eq.2) then
	  r=gama
	  s=beta
	  t=beta
	wt=9.d0/20.d0
	elseif (kk.eq.3) then
	  r=beta
	  s=gama
	  t=beta
	wt=9.d0/20.d0
	elseif (kk.eq.4) then 
	  r=beta
	  s=beta
	  t=gama
	wt=9.d0/20.d0
	elseif (kk.eq.5) then 
	  r=beta
	  s=beta
	  t=beta
	wt=9.d0/20.d0
	endif


	end
c==============================================================================
c==============================================================================
C=======================================================================
      SUBROUTINE PENTP3Tetra(GNODE,NEL,CORD,ID,PRIT,NET,NDIM,
     &IDPRIT,PENALT,IIZLAZ,AMI,ISRPS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /PENALL/ PENALT,PRESS
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM

      DIMENSION GNODE(2,6,*),NEL(NDIM+1,*),CORD(3,*),ID(6,*)
      DIMENSION PRIT(IDPRIT,*)
      DIMENSION TT21(92),PJ(3,21),H(21),HP(10),SHEAR(3)
	dimension ck(21,3)

C
CE Subroutine PENTP3 is used for postprocessing pressure calculation
CE for 3D analysis when Penalty method is included
C




      

c      DO 505 NBREL=1,NET
      DO 505 NBREL=1,NET
C inicijalizacija pritiska za primere sa krvnim sudovima
C inicijalni pritisak = 80 mmHg = 10 kPa
C      PRIT(1,NBREL)=0.D0

      DO I=1,92
        TT21(I)=0.D0
      ENDDO


C=========================================================================
      DO  KLM=1,NDIM
      CK(KLM,1)=CORD(1,NEL(KLM,NBREL))
      CK(KLM,2)=CORD(2,NEL(KLM,NBREL))
      CK(KLM,3)=CORD(3,NEL(KLM,NBREL))
C      LM2(KLM)=ID(1,NEL(KLM,NBREL))
C      LM2(KLM+NDIM)=ID(2,NEL(KLM,NBREL))
C      LM2(KLM+2*NDIM)=ID(3,NEL(KLM,NBREL))
C      IF (KLM.LE.8) LM2(KLM+3*NDIM)=ID(4,NEL(KLM,NBREL))
      ENDDO
C=======================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,3
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
C      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
 135  CONTINUE
C      IF (KLM.LE.8.AND.ID(4,NEL(KLM,NBREL)).NE.0) THEN
C      IF (KLM.LE.8) THEN
C        TT21(KLM+3*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
C      ENDIF
        TT21(KLM+3*NDIM+8)=GNODE(2,5,NEL(KLM,NBREL))
C        TT210(KLM+3*NDIM+8)=TT10(ID(5,NEL(KLM,NBREL)))
 140  CONTINUE
C=======================================================================

C
C       IUPWIN=0 
c	  wt=-0.8d0
C      nkk=4
      nkk=4
	ppp=0.d0
	do kk=1,nkk
         pp=0.d0
        IF (nkk.eq.4) then
        call Integ4points(r,s,t,wt,kk)      
	elseif (nkk.eq.11) then
        call Integ10points(r,s,t,wt,kk)      
	else if (nkk.eq.1) then
	  r=0.25d0
	  s=0.25d0
	  t=0.25d0
      endif
      if (ndim.eq.8) then
       call jact8nodeTetra (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,iizlaz,ndimp)
c       CALL JACT(0.D0,0.D0,0.D0,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
c     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
c     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN)
	elseif (ndim.eq.4) then
       call jact4node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL,1,hpv)
	elseif (ndim.eq.10) then
       call jact10node (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,hp,hv1,
     &hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,NEL)
	  
	endif

C       WDT=DET1*WT       
       DO I=1,NDIM
        I1=I
        I2=I+NDIM
        I3=I+2*NDIM
        PPP=-
     &PENALT*(PJ(1,I)*TT21(I1)+PJ(2,I)*TT21(I2)+PJ(3,I)*TT21(I3))
        PRIT(1,NBREL)=PRIT(1,NBREL)+PPP
        pp=pp+ppp 
       ENDDO
	 write(iizlaz,'(2i10,e13.5)')nbrel,kk,pp
	enddo
C        PRIT(1,NBREL)=PPP

505    CONTINUE

C END OF LOOP PER ELEMENT    
      END
C=========================================================================
