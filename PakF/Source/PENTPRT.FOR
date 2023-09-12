C=======================================================================
      SUBROUTINE PRIT2DT(GNODE,NEL,CORD,ID,TT10,PRIT,INDAX,NDIM,
     &IDPRIT,PENALT,NUMZAD,NPT,MAXSIL,JEDN,NWK1,NET,AMI,INDAMI,IIZLAZ,
     &IUPWIN,GUSM,AKT)
     
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
C
CE Subroutine PRIT2DT is used for postprocessing pressure calculation
CE for 2D analysis IN TURBULENT MODE when Penalty method is included
C
      DIMENSION GNODE(2,7,*),NEL(NDIM+1,*),ID(7,*),TT10(*)
      DIMENSION CORD(3,*),PRIT(IDPRIT,*)

      DIMENSION H(9),ZVHX(9),ZVHY(9),HP(4)
      DIMENSION R(3,3),S(3,3),W(3,3),X(9),Y(9),TT21(54),TT210(54)
      
      R(3,1)=-0.7745966692415
      R(3,2)=0.0
      R(3,3)=0.77459666924148
      S(3,1)=-0.7745966692415
      S(3,2)=0.0
      S(3,3)=0.77459666924148
      W(3,1)=0.55555555555556
      W(3,2)=0.88888888888889
      W(3,3)=0.55555555555556


      R(2,1)=0.57735026918963
      R(2,2)=-0.5773502691896
      S(2,1)=0.57735026918963
      S(2,2)=-0.5773502691896
      W(2,1)=1.000000000000000
      W(2,2)=1.000000000000000

      R(1,1)=0.0
      S(1,1)=0.0
      W(1,1)=2.0

      DO 400 NBREL=1,NET

C TT21 is vector of unknowns values at element level
C TT210 is vector of unknowns values at element level at start of time step

      DO 125 I=1,4*NDIM
      TT210(I)=0.D0
      TT21(I)=0.D0
 125  CONTINUE
C=========================================================================
      DO 130 KLM=1,NDIM
      X(KLM)=CORD(1,NEL(KLM,NBREL))
      Y(KLM)=CORD(2,NEL(KLM,NBREL))
 130  CONTINUE
C=======================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,2
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
 135  CONTINUE
      IF (KLM.LE.4) THEN
      TT21(KLM+2*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
      TT210(KLM+2*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
      TT21(KLM+2*NDIM+4)=GNODE(2,5,NEL(KLM,NBREL))
      TT210(KLM+2*NDIM+4)=GNODE(1,5,NEL(KLM,NBREL))
 140  CONTINUE
 
C RACUNANJE PRITISAKA 
C calculation of pressure if we have penalty function

      IBRGT=3
      IF (NDIM.EQ.4) IBRGT=2
      K=0
      DO 450 I=1,IBRGT-1
      DO 450 J=1,IBRGT-1
      K=K+1
     
      CALL FNTERPT(R(IBRGT-1,I),S(IBRGT-1,J),0,TT21,H,HP,ZVHX,ZVHY,HV2,
     &HV3,ZVXT,ZVYT,DETJ,DETJS,X,Y,FS2,FS3,ZVXV2,ZVYV3,ZVYV2,ZVXV3,TAU,
     &NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,AMI,INDAMI,PENALT,PRESS,INDAX,
     &AKT,GUSM,IUPWIN,FSK1,FSOMEGA1,HKT,HOM,ZVXK,ZVYK,ZVXOM,ZVYOM,NBREL,
     &IIZLAZ,GNODE,NEL,CORD,ID,TT210)
              
      PRIT(K,NBREL)=PRESS

 450  CONTINUE
								  
 400  CONTINUE      
C======================================================================= 
C End of subroutine
      END
C=======================================================================