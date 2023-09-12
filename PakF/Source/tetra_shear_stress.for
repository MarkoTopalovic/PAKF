C==========================================================================
      SUBROUTINE tetra_SSTRES(NEL,NDIM,ID,CK,PJ,H,HP,TT21,AMI,ISRPS,
     &NUMZAD,IIZLAZ,IUPWIN,XG,WGT,NREF,IBRGT,PRES,INDEL,NBREL,NZAD,
     &ZADVRE,CORD)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION NEL(NDIM+1,*),ID(6,*),INDEL(*),NZAD(3,*)
      DIMENSION CK(21,*),PJ(3,*),H(*),HP(*),TT21(*),PRES(3,*),ZADVRE(*)
      DIMENSION CORD(3,*)
	DIMENSION NID(21),N(21),NWALL(6),ITR(6,4)
      DIMENSION XG(*),WGT(*),NREF(*),SHEAR(3),NIZEL(3)

C
CE Subroutine SSTRES is used for calculation shear stresses
C

C      WRITE(IIZLAZ,*)'ELEMENT= ',NBREL
      
c      do i=1,11716
c        write(iizlaz,*)i,indel(i)
c      enddo
c      stop
      
	DO I=1,NDIM
	  N(I)=NEL(I,NBREL)
        NID(I)=0
      IF(ID(1,N(I)).EQ.1.AND.ID(2,N(I)).EQ.1.AND.ID(3,N(I)).EQ.1)
     &NID(I)=1
      ENDDO




      ITR(1,1)=N(1)
      ITR(1,2)=N(3)
      ITR(1,3)=N(4)

      ITR(2,1)=N(1)
      ITR(2,2)=N(2)
      ITR(2,3)=N(4)

      ITR(3,1)=N(1)
      ITR(3,2)=N(2)
      ITR(3,3)=N(3)

      ITR(4,1)=N(2)
      ITR(4,2)=N(3)
      ITR(4,3)=N(4)


      NWALL(1)=NID(1)*NID(3)*NID(4)
      NWALL(2)=NID(1)*NID(2)*NID(4)
      NWALL(3)=NID(1)*NID(2)*NID(3)
      NWALL(4)=NID(2)*NID(3)*NID(4)

      DO 300 NPOV=1,4
       IF (NWALL(NPOV).EQ.0) GOTO 300

C
CS  PETLJA PO GAUSOVIM TACKAMA
CE  GAUSS POINTS LOOP
C
      DO 210 KK=1,3
C
      NODE=ITR(NPOV,KK)
      NIZEL(1)=ITR(NPOV,1)
      NIZEL(2)=ITR(NPOV,2)
      NIZEL(3)=ITR(NPOV,3)
      
      call strana_tetra2(nizel,nel,nbrel,r,s,t,wt,kk,ndim,npov)
      call ShearStressSurf4node(r,s,t,ck,h,tt21,nbrel,ndim,det1,
     &nel,SHEAR,nizel,pj,hp,ami,ndimp,npov,cord,iizlaz)


	PRES(1,NODE)=PRES(1,NODE)-AMI*SHEAR(1)/INDEL(NODE)
      PRES(2,NODE)=PRES(2,NODE)-AMI*SHEAR(2)/INDEL(NODE)
      PRES(3,NODE)=PRES(3,NODE)-AMI*SHEAR(3)/INDEL(NODE)
c
c	PRES(1,NODE)=PRES(1,NODE)+indel(node)
c      PRES(2,NODE)=PRES(2,NODE)+indel(node)
c      PRES(3,NODE)=PRES(3,NODE)+indel(node)
      

  210 CONTINUE    

  300 CONTINUE    

      END
C==========================================================================
c==============================================================================
      subroutine ShearStressSurf4node (r,s,t,ck,h,tt21,nbrel,ndim,det1,
     &nel,shear,nizel,pj,hp,ami,ndimp,npov,cord,iizlaz)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension p(4,10),pj(3,*),xjj(3,3),xj(3,3),ck(21,*)
      dimension me(16),le(16)
      dimension xjja(20),Xjjj(17)
      dimension h(*),hp(*),TT21(*),NEL(NDIM+1,*)
      dimension cord(3,*)
      dimension V1(10),V2(10),V3(10),pp(10),temp(10),an(3)
      DIMENSION VSHX(3),VSHY(3),VSHZ(3),SHEAR(*),vtang(3,10),TANG(3)
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



       h(1)=e1
       h(2)=e2
       h(3)=e3
       h(4)=e4
	 hp(1)=h(1)
	 hp(2)=h(2)
	 hp(3)=h(3)
	 hp(4)=h(4)

c	 hpv(1)=h(1)
c	 hpv(2)=h(2)
c	 hpv(3)=h(3)
c	 hpv(4)=h(4)

       if (ndimp.eq.1) hp(1)=1.d0

       p(1,1)=-1.0d0
       p(2,1)=-1.0d0
       p(3,1)=-1.0d0

       p(1,2)=1.0d0
       p(2,3)=1.0d0
       p(3,4)=1.0d0




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

c       k=1
c       do j=1,3
c		   do i=1,3
c		   xjja(k)=xjj(i,j) 
c		   k=k+1
c		enddo
c	 enddo


c      call minv(xjja,3,det1,le,me)
      call minv(xjj,3,det1,le,me)

c       k=1
c      do j=1,3
c		   do i=1,3
c			   xjj(i,j)=xjja(k)
c			   k=k+1
c	       enddo
c	   enddo
		  

	  
      if (det1<1.d-15) then
c       write (*,*) 'determinante less then zero for element',nbrel
c       return
      endif
      
      
	   

c	xxjjj(1,1)=1.d0  
c	xxjjj(1,2)=1.d0  
c	xxjjj(1,3)=1.d0  
c	xxjjj(1,4)=1.d0  
	  
c          do j=1,3
c           do i=1,4
c	        xxjjj(j+1,i)=ck(i,j)
c			k=k+1
c		   enddo
c	   enddo

c       call minv(xxjjj,4,det2,le,me);
c
       AJ_1=(1.0/6.0)*det1
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
      
      
c      IF (npov.eq.1) THEN
c       N2=3
c       N1=1
c      ENDIF
c      IF (npov.eq.2.or.npov.eq.3) THEN
c       N2=2
c       N1=1
c      ENDIF
c      IF (npov.eq.4) THEN
c       N2=3
c       N1=2
c      ENDIF
      
      
c      TANG(1)=CK(N2,1)-CK(N1,1)
c      TANG(2)=CK(N2,2)-CK(N1,2)
c      TANG(3)=CK(N2,3)-CK(N1,3)
c      TAN=DSQRT(TANG(1)**2+TANG(2)**2+TANG(3)**2)
c      TANG(1)=TANG(1)/TAN
c      TANG(2)=TANG(2)/TAN
c      TANG(3)=TANG(3)/TAN
      
c      DO I=1,NDIM
c       VVTANG=V1(I)*TANG(1)+V2(I)*TANG(2)+V3(I)*TANG(3)
c       VTANG(1,I)=VVTANG*TANG(1)
c       VTANG(2,I)=VVTANG*TANG(2)
c       VTANG(3,I)=VVTANG*TANG(3)
c      ENDDO

      CALL SHEARS(AN,TT21,VTANG,NDIM)

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


c      return
c==============================================================================
      vx1=0.d0
      vy1=0.d0
      vz1=0.d0
      do i=1,ndim
c      write(iizlaz,'(2i5,3e13.5)') 
c     &i,nbrel,vtang(1,i),vtang(2,i),vtang(3,i)
       vx1=vx1+vtang(1,i)
       vy1=vy1+vtang(2,i)
       vz1=vz1+vtang(3,i)
c       vx1=vx1+v1(i)
c       vy1=vy1+v2(i)
c       vz1=vz1+v3(i)
      enddo
      shear(1)=vx1
      shear(2)=vy1
      shear(3)=vz1
c==============================================================================
c      shear(1)=an(1)
c      shear(2)=an(2)
c      shear(3)=an(3)

       end
c==============================================================================
c==============================================================================
C==========================================================================
      SUBROUTINE Tetra_INDELSSTRES(NEL,INDEL,NET,NPT,NDIM,ID)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
     	
      DIMENSION NEL(NDIM+1,*),INDEL(*),ID(6,*)
	DIMENSION NID(21),N(21),NWALL(4),ITR(4,4)

     
C
CE Subroutine INITIA is used for initialisation all global variable
C


C INITIALISATION:
       DO I=1,NPT
        INDEL(I)=0
       ENDDO



      do nbrel=1,net

	DO I=1,NDIM
	  N(I)=NEL(I,NBREL)
        NID(I)=0
      IF(ID(1,N(I)).EQ.1.AND.ID(2,N(I)).EQ.1.AND.ID(3,N(I)).EQ.1)
     &NID(I)=1
      ENDDO

      ITR(1,1)=N(1)
      ITR(1,2)=N(3)
      ITR(1,3)=N(4)

      ITR(2,1)=N(1)
      ITR(2,2)=N(2)
      ITR(2,3)=N(4)

      ITR(3,1)=N(1)
      ITR(3,2)=N(2)
      ITR(3,3)=N(3)

      ITR(4,1)=N(2)
      ITR(4,2)=N(3)
      ITR(4,3)=N(4)


      NWALL(1)=NID(1)*NID(3)*NID(4)
      NWALL(2)=NID(1)*NID(2)*NID(4)
      NWALL(3)=NID(1)*NID(2)*NID(3)
      NWALL(4)=NID(2)*NID(3)*NID(4)

      DO 300 NPOV=1,4
       IF (NWALL(NPOV).EQ.0) GOTO 300
      do j=1,3
	 node=itr(npov,j)
	 indel(node)=indel(node)+1
      enddo
300   continue
      enddo

      

      END
C==========================================================================
c==============================================================================
      subroutine strana_tetra2(nizel,nel,nbrel,r,s,t,wt,kk,ndim,ii)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      dimension nel(ndim+1,*),nizm(3,4),nizel(*)
	data nizm/1,3,4,1,2,4,1,2,3,2,3,4/

            
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
