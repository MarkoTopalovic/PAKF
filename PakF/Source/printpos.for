c==============================================================================
      subroutine printpos(gnode,npt,nel,net,ndim,cord,penalt,id,pres,
     &spsil)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),id(6,*),spsil(3,*)
      dimension pres(3,*)
      DIMENSION NEL(NDIM+1,*)
      
      if(NDIM.eq.21) then
            call printposvelocity21(gnode,npt,nel,net,ndim,cord,id)
      else
            call printposvelocity(gnode,npt,nel,net,ndim,cord,id)
      endif
      if (ndim.eq.10) then
	  call printpospressure(gnode,npt,nel,net,ndim,cord)
c	elseif(ndim.eq.4.and.penalt.lt.1d0) then
c	  call printpospressureELEM(gnode,npt,nel,net,ndim,cord)
c	elseif(ndim.eq.4.and.penalt.gt.1d0) then
      else
	  call printpospressure(gnode,npt,nel,net,ndim,cord)
	  call printposshearstress(pres,npt,nel,net,ndim,cord,id)
	  call printposwallforce(spsil,npt,nel,net,ndim,cord,id)
	endif

      end
c==============================================================================
c==============================================================================
c==============================================================================
      subroutine printpospressure(gnode,npt,nel,net,ndim,cord)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*)


      ifile=14
      open (ifile,file='Pressure.pos')
      write(ifile,*)'View "PressureTube" {'
	do nbrel=1,net


	  n1=nel(1,nbrel)
	  n2=nel(2,nbrel)
	  n3=nel(3,nbrel)
	  xn1=cord(1,n1)
	  yn1=cord(2,n1)
	  zn1=cord(3,n1)
	  xn2=cord(1,n2)
	  yn2=cord(2,n2)
	  zn2=cord(3,n2)
	  xn3=cord(1,n3)
	  yn3=cord(2,n3)
	  zn3=cord(3,n3)
	  p1=gnode(2,4,n1)
	  p2=gnode(2,4,n2)
	  p3=gnode(2,4,n3)

      if (nel(3,nbrel).eq.nel(4,nbrel)) then
	  n4=nel(5,nbrel)
	  xn4=cord(1,n4)
	  yn4=cord(2,n4)
	  zn4=cord(3,n4)
	  p4=gnode(2,4,n4)
	  
	  n5=nel(6,nbrel)
	  xn5=cord(1,n5)
	  yn5=cord(2,n5)
	  zn5=cord(3,n5)
	  p5=gnode(2,4,n5)
	  
	  n6=nel(7,nbrel)
	  xn6=cord(1,n6)
	  yn6=cord(2,n6)
	  zn6=cord(3,n6)
	  p6=gnode(2,4,n6)
	  write(ifile,
     &'(a3,17(e13.5,a1),e13.5,a2,5(e13.5,a1),e13.5,a2)')'SI(',
     &xn1,',',yn1,',',zn1,',',xn2,',',yn2,',',zn2,',',xn3,',',yn3,',',
     &zn3,',',xn4,',',yn4,',',zn4,',',xn5,',',yn5,',',zn5,
     &',',xn6,',',yn6,',',zn6,
     &'){',p1,',',p2,',',p3,',',p4,',',p5,',',p6,'};'
      else
	  n4=nel(4,nbrel)
	  xn4=cord(1,n4)
	  yn4=cord(2,n4)
	  zn4=cord(3,n4)
	  p4=gnode(2,4,n4)
	  write(ifile,
     &'(a3,11(e13.5,a1),e13.5,a2,3(e13.5,a1),e13.5,a2)')'SS(',
     &xn1,',',yn1,',',zn1,',',xn2,',',yn2,',',zn2,',',xn3,',',yn3,',',
     &zn3,',',xn4,',',yn4,',',zn4,
     &'){',p1,',',p2,',',p3,',',p4,'};'
      endif
	enddo
	  write(ifile,'(a2)') '};'

C
      
      close(ifile)
       end
c==============================================================================
c==============================================================================
      subroutine printpospressureELEM(gnode,npt,nel,net,ndim,cord)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*)


      ifile=14
      open (ifile,file='Pressure.pos')
      write(ifile,*)'View "PressureTube" {'
	do nbrel=1,net
	  n1=nel(1,nbrel)
	  n2=nel(2,nbrel)
	  n3=nel(3,nbrel)
	  n4=nel(4,nbrel)
	  n5=nel(5,nbrel)
	  xn1=cord(1,n1)
	  yn1=cord(2,n1)
	  zn1=cord(3,n1)
	  xn2=cord(1,n2)
	  yn2=cord(2,n2)
	  zn2=cord(3,n2)
	  xn3=cord(1,n3)
	  yn3=cord(2,n3)
	  zn3=cord(3,n3)
	  xn4=cord(1,n4)
	  yn4=cord(2,n4)
	  zn4=cord(3,n4)
	  p1=gnode(2,4,n5)
	  p2=gnode(2,4,n5)
	  p3=gnode(2,4,n5)
	  p4=gnode(2,4,n5)
	  write(ifile,
     &'(a3,11(e13.5,a1),e13.5,a2,3(e13.5,a1),e13.5,a2)')'SS(',
     &xn1,',',yn1,',',zn1,',',xn2,',',yn2,',',zn2,',',xn3,',',yn3,',',
     &zn3,',',xn4,',',yn4,',',zn4,
     &'){',p1,',',p2,',',p3,',',p4,'};'
	enddo
	  write(ifile,'(a2)') '};'

C
      
      close(ifile)
       end
c==============================================================================
c==============================================================================
c==============================================================================
      subroutine printposvelocity(gnode,npt,nel,net,ndim,cord,id)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),id(6,*)
      DIMENSION NEL(NDIM+1,*)


      ifile=14
      open (ifile,file='EffectiveVelocity.pos')
      write(ifile,*)'View "EffectiveVelocityTube" {'
	do nbrel=1,net

	  n1=nel(1,nbrel)
	  n2=nel(2,nbrel)
	  n3=nel(3,nbrel)
	  xn1=cord(1,n1)
	  yn1=cord(2,n1)
	  zn1=cord(3,n1)
	  xn2=cord(1,n2)
	  yn2=cord(2,n2)
	  zn2=cord(3,n2)
	  xn3=cord(1,n3)
	  yn3=cord(2,n3)
	  zn3=cord(3,n3)
	  vx1=gnode(2,1,n1)
	  vy1=gnode(2,2,n1)
	  vz1=gnode(2,3,n1)
	  v1=dsqrt(vx1**2+vy1**2+vz1**2)
	  vx2=gnode(2,1,n2)
	  vy2=gnode(2,2,n2)
	  vz2=gnode(2,3,n2)
	  v2=dsqrt(vx2**2+vy2**2+vz2**2)
	  vx3=gnode(2,1,n3)
	  vy3=gnode(2,2,n3)
	  vz3=gnode(2,3,n3)
	  v3=dsqrt(vx3**2+vy3**2+vz3**2)

      if (nel(3,nbrel).eq.nel(4,nbrel)) then
	  n4=nel(5,nbrel)
	  n5=nel(6,nbrel)
	  n6=nel(7,nbrel)
	  xn4=cord(1,n4)
	  yn4=cord(2,n4)
	  zn4=cord(3,n4)
	  xn5=cord(1,n5)
	  yn5=cord(2,n5)
	  zn5=cord(3,n5)
	  xn6=cord(1,n6)
	  yn6=cord(2,n6)
	  zn6=cord(3,n6)
	  vx4=gnode(2,1,n4)
	  vy4=gnode(2,2,n4)
	  vz4=gnode(2,3,n4)
	  v4=dsqrt(vx4**2+vy4**2+vz4**2)
	  vx5=gnode(2,1,n5)
	  vy5=gnode(2,2,n5)
	  vz5=gnode(2,3,n5)
	  v5=dsqrt(vx5**2+vy5**2+vz5**2)
	  vx6=gnode(2,1,n6)
	  vy6=gnode(2,2,n6)
	  vz6=gnode(2,3,n6)
	  v6=dsqrt(vx6**2+vy6**2+vz6**2)

	  write(ifile,
     &'(a3,17(e13.5,a1),e13.5,a2,5(e13.5,a1),e13.5,a2)')'SI(',
     &xn1,',',yn1,',',zn1,',',xn2,',',yn2,',',zn2,',',xn3,',',yn3,',',
     &zn3,',',xn4,',',yn4,',',zn4,',',xn5,',',yn5,',',zn5,',',xn6,',',
     &yn6,',',zn6,
     &'){',v1,',',v2,',',v3,',',v4,',',v5,',',v6,'};'
      else
	  n4=nel(4,nbrel)
	  xn4=cord(1,n4)
	  yn4=cord(2,n4)
	  zn4=cord(3,n4)
	  vx4=gnode(2,1,n4)
	  vy4=gnode(2,2,n4)
	  vz4=gnode(2,3,n4)
	  v4=dsqrt(vx4**2+vy4**2+vz4**2)
	  write(ifile,
     &'(a3,11(e13.5,a1),e13.5,a2,3(e13.5,a1),e13.5,a2)')'SS(',
     &xn1,',',yn1,',',zn1,',',xn2,',',yn2,',',zn2,',',xn3,',',yn3,',',
     &zn3,',',xn4,',',yn4,',',zn4,
     &'){',v1,',',v2,',',v3,',',v4,'};'
      endif
     	enddo
	  write(ifile,'(a2)') '};'


C
C
      close(ifile)

      end
c==============================================================================
c==============================================================================
c==============================================================================
      subroutine printposvelocity21(gnode,npt,nel,net,ndim,cord,id)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION GNODE(2,6,*),CORD(3,*),id(6,*)
      DIMENSION NEL(NDIM+1,*)


      ifile=14
      open (ifile,file='EffectiveVelocity.pos')
      write(ifile,*)'View "EffectiveVelocity" {'
	do nbrel=1,net
	  n1=nel(1,nbrel)
	  n2=nel(2,nbrel)
	  n3=nel(3,nbrel)
	  n4=nel(8,nbrel)
	  xn1=cord(1,n1)
	  yn1=cord(2,n1)
	  zn1=cord(3,n1)
	  xn2=cord(1,n2)
	  yn2=cord(2,n2)
	  zn2=cord(3,n2)
	  xn3=cord(1,n3)
	  yn3=cord(2,n3)
	  zn3=cord(3,n3)
	  xn4=cord(1,n4)
	  yn4=cord(2,n4)
	  zn4=cord(3,n4)
	  vx1=gnode(2,1,n1)
	  vy1=gnode(2,2,n1)
	  vz1=gnode(2,3,n1)
	  v1=dsqrt(vx1**2+vy1**2+vz1**2)
	  vx2=gnode(2,1,n2)
	  vy2=gnode(2,2,n2)
	  vz2=gnode(2,3,n2)
	  v2=dsqrt(vx2**2+vy2**2+vz2**2)
	  vx3=gnode(2,1,n3)
	  vy3=gnode(2,2,n3)
	  vz3=gnode(2,3,n3)
	  
	  
	  v3=dsqrt(vx3**2+vy3**2+vz3**2)
	  vx4=gnode(2,1,n4)
	  vy4=gnode(2,2,n4)
	  vz4=gnode(2,3,n4)
	  v4=dsqrt(vx4**2+vy4**2+vz4**2)

c        v1=0.d0
c        v2=0.d0
c        v3=0.d0
c       v4=0.d0
c	  if (id(3,n1).ne.0) v1=1.d0
c	  if (id(3,n2).ne.0) v2=1.d0
c	  if (id(3,n3).ne.0) v3=1.d0
c	  if (id(3,n4).ne.0) v4=1.d0

	  write(ifile,
     &'(a3,11(e13.5,a1),e13.5,a2,3(e13.5,a1),e13.5,a2)')'SS(',
     &xn1,',',yn1,',',zn1,',',xn2,',',yn2,',',zn2,',',xn3,',',yn3,',',
     &zn3,',',xn4,',',yn4,',',zn4,
     &'){',v1,',',v2,',',v3,',',v4,'};'
	enddo
	  write(ifile,'(a2)') '};'


C
C
      close(ifile)

      end
c==============================================================================
c==============================================================================
      subroutine printposshearstress(pres,npt,nel,net,ndim,cord,id)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION pres(3,*),CORD(3,*),id(6,*)
      DIMENSION NEL(NDIM+1,*)


      ifile=19
      open (ifile,file='ShearStress.pos')
      write(ifile,*)'View "ShearStress" {'
	do nbrel=1,net
	  n1=nel(1,nbrel)
	  n2=nel(2,nbrel)
	  n3=nel(3,nbrel)
	  xn1=cord(1,n1)
	  yn1=cord(2,n1)
	  zn1=cord(3,n1)
	  xn2=cord(1,n2)
	  yn2=cord(2,n2)
	  zn2=cord(3,n2)
	  xn3=cord(1,n3)
	  yn3=cord(2,n3)
	  zn3=cord(3,n3)
	  sx1=pres(1,n1)
	  sy1=pres(2,n1)
	  sz1=pres(3,n1)
	  s1=dsqrt(sx1**2+sy1**2+sz1**2)
	  sx2=pres(1,n2)
	  sy2=pres(2,n2)
	  sz2=pres(3,n2)
	  s2=dsqrt(sx2**2+sy2**2+sz2**2)
	  sx3=pres(1,n3)
	  sy3=pres(2,n3)
	  sz3=pres(3,n3)
  	  s3=dsqrt(sx3**2+sy3**2+sz3**2)

      if (nel(3,nbrel).eq.nel(4,nbrel)) then
	  n4=nel(5,nbrel)
	  n5=nel(6,nbrel)
	  n6=nel(7,nbrel)
	  xn4=cord(1,n4)
	  yn4=cord(2,n4)
	  zn4=cord(3,n4)
	  xn5=cord(1,n5)
	  yn5=cord(2,n5)
	  zn5=cord(3,n5)
	  xn6=cord(1,n6)
	  yn6=cord(2,n6)
	  zn6=cord(3,n6)
	  sx4=pres(1,n4)
	  sy4=pres(2,n4)
	  sz4=pres(3,n4)
	  s4=dsqrt(sx4**2+sy4**2+sz4**2)
	  sx5=pres(1,n5)
	  sy5=pres(2,n5)
	  sz5=pres(3,n5)
	  s5=dsqrt(sx5**2+sy5**2+sz5**2)
	  sx6=pres(1,n6)
	  sy6=pres(2,n6)
	  sz6=pres(3,n6)
  	  s6=dsqrt(sx6**2+sy6**2+sz6**2)
	  write(ifile,
     &'(a3,17(e13.5,a1),e13.5,a2,5(e13.5,a1),e13.5,a2)')'SI(',
     &xn1,',',yn1,',',zn1,',',xn2,',',yn2,',',zn2,',',xn3,',',yn3,',',
     &zn3,',',xn4,',',yn4,',',zn4,',',xn5,',',yn5,',',zn5,
     &',',xn6,',',yn6,',',zn6,
     &'){',s1,',',s2,',',s3,',',s4,',',s5,',',s6,'};'
      else

	  n4=nel(4,nbrel)
	  xn4=cord(1,n4)
	  yn4=cord(2,n4)
	  zn4=cord(3,n4)
	  sx4=pres(1,n4)
	  sy4=pres(2,n4)
	  sz4=pres(3,n4)
	  s4=dsqrt(sx4**2+sy4**2+sz4**2)


	  write(ifile,
     &'(a3,11(e13.5,a1),e13.5,a2,3(e13.5,a1),e13.5,a2)')'SS(',
     &xn1,',',yn1,',',zn1,',',xn2,',',yn2,',',zn2,',',xn3,',',yn3,',',
     &zn3,',',xn4,',',yn4,',',zn4,
     &'){',s1,',',s2,',',s3,',',s4,'};'


      endif
	enddo
	  write(ifile,'(a2)') '};'


C
C
      close(ifile)

      end
c==============================================================================
c==============================================================================
c==============================================================================
      subroutine printposwallforce(spsil,npt,nel,net,ndim,cord,id)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION spsil(3,*),CORD(3,*),id(6,*)
      DIMENSION NEL(NDIM+1,*)


      ifile=14
      if (ndim.ne.8) return
      open (ifile,file='WallForce.pos')
      write(ifile,*)'View "Wall Force" {'
	do nbrel=1,net
	  n1=nel(1,nbrel)
	  n2=nel(2,nbrel)
	  n3=nel(3,nbrel)
	  n4=nel(4,nbrel)
	  n5=nel(5,nbrel)
	  n6=nel(6,nbrel)
	  n7=nel(7,nbrel)
	  n8=nel(8,nbrel)
	  xn1=cord(1,n1)
	  yn1=cord(2,n1)
	  zn1=cord(3,n1)
	  xn2=cord(1,n2)
	  yn2=cord(2,n2)
	  zn2=cord(3,n2)
	  xn3=cord(1,n3)
	  yn3=cord(2,n3)
	  zn3=cord(3,n3)
	  sx1=spsil(1,n1)
	  sy1=spsil(2,n1)
	  sz1=spsil(3,n1)
	  s1=dsqrt(sx1**2+sy1**2+sz1**2)
	  sx2=spsil(1,n2)
	  sy2=spsil(2,n2)
	  sz2=spsil(3,n2)
	  s2=dsqrt(sx2**2+sy2**2+sz2**2)
	  sx3=spsil(1,n3)
	  sy3=spsil(2,n3)
	  sz3=spsil(3,n3)
  	  s3=dsqrt(sx3**2+sy3**2+sz3**2)

	  xn4=cord(1,n4)
	  yn4=cord(2,n4)
	  zn4=cord(3,n4)
	  xn5=cord(1,n5)
	  yn5=cord(2,n5)
	  zn5=cord(3,n5)
	  xn6=cord(1,n6)
	  yn6=cord(2,n6)
	  zn6=cord(3,n6)
	  sx4=spsil(1,n4)
	  sy4=spsil(2,n4)
	  sz4=spsil(3,n4)
	  s4=dsqrt(sx4**2+sy4**2+sz4**2)
	  sx5=spsil(1,n5)
	  sy5=spsil(2,n5)
	  sz5=spsil(3,n5)
	  s5=dsqrt(sx5**2+sy5**2+sz5**2)
	  sx6=spsil(1,n6)
	  sy6=spsil(2,n6)
	  sz6=spsil(3,n6)
  	  s6=dsqrt(sx6**2+sy6**2+sz6**2)

	  xn7=cord(1,n7)
	  yn7=cord(2,n7)
	  zn7=cord(3,n7)
	  xn8=cord(1,n8)
	  yn8=cord(2,n8)
	  zn8=cord(3,n8)
	  sx7=spsil(1,n7)
	  sy7=spsil(2,n7)
	  sz7=spsil(3,n7)
	  s7=dsqrt(sx7**2+sy7**2+sz7**2)
	  sx8=spsil(1,n8)
	  sy8=spsil(2,n8)
	  sz8=spsil(3,n8)
  	  s8=dsqrt(sx8**2+sy8**2+sz8**2)



	  write(ifile,
     &'(a3,23(e13.5,a1),e13.5,a2,7(e13.5,a1),e13.5,a2)')'SH(',
     &xn1,',',yn1,',',zn1,',',xn2,',',yn2,',',zn2,',',xn3,',',yn3,',',
     &zn3,',',xn4,',',yn4,',',zn4,',',xn5,',',yn5,',',zn5,
     &',',xn6,',',yn6,',',zn6,',',xn7,',',yn7,',',zn7,
     &',',xn8,',',yn8,',',zn8,
     &'){',s1,',',s2,',',s3,',',s4,',',s5,',',s6,',',s7,',',s8,'};'
      

	enddo
	  write(ifile,'(a2)') '};'


C
C
      close(ifile)

      end
c==============================================================================
c==============================================================================
