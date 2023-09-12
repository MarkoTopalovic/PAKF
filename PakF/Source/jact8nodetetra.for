c==============================================================================
      subroutine jact8nodeTetra (r,s,t,ck,h,tt21,nbrel,ndim,det1,pj,
     &hp,hv1,hv2,hv3,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,nel,
     &iizlaz,ndimp)
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
      PRIT=DOT(HP,PP,ndimp)
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
