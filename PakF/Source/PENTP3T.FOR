C=======================================================================
      SUBROUTINE PENTP3T(GNODE,NEL,CORD,ID,PRIT,NET,NDIM,IDPRIT,PENALT,
     &IIZLAZ,AMI,ISRPS)
     
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      COMMON /NUMNPT/ NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET
C      COMMON /PENALL/ PENALT,PRESS
C      COMMON /NDESUK/ NDES,IDPRIT,IFORM

      DIMENSION GNODE(2,7,*),NEL(NDIM+1,*),CORD(3,*),ID(7,*)
      DIMENSION PRIT(IDPRIT,*)
C      DIMENSION PRITISAK(IDPRIT,NET)
      DIMENSION TT21(92),PJ(3,21),H(21),HP(10),SHEAR(3)
	dimension ck(21,3)

C
CE Subroutine PENTP3 is used for postprocessing pressure calculation
CE for 3D analysis when Penalty method is included
      DO 505 NBREL=1,NET
C inicijalizacija pritiska za primere sa krvnim sudovima
C inicijalni pritisak = 80 mmHg = 10 kPa

      DO I=1,6*NDIM
        TT21(I)=0.D0
        
      ENDDO
C=========================================================================
      DO  KLM=1,NDIM
      CK(KLM,1)=CORD(1,NEL(KLM,NBREL))
      CK(KLM,2)=CORD(2,NEL(KLM,NBREL))
      CK(KLM,3)=CORD(3,NEL(KLM,NBREL))
      ENDDO
C=======================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,3
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
 135  CONTINUE
      TT21(KLM+3*NDIM+8)=GNODE(2,5,NEL(KLM,NBREL))
 140  CONTINUE
C=======================================================================
C	  r=0.25d0
C	  s=0.25d0
C	  t=0.25d0

      IF (NETIP.EQ.3.AND.NDIM.EQ.8) IBRGT=2
      
C      K=0
C      DO 450 
C      K=K+1
      
       CALL JACTT(0.D0,0.D0,0.D0,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21,
     &DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,
     &FSK1,FSOMEGA1,HKT,HOM,ZVXK,ZVYK,ZVZK,ZVXOM,ZVYOM,ZVZOM,
     &ZVXV1,ZVXV2,ZVXV3,ZVYV1,ZVYV2,ZVYV3,ZVZV1,ZVZV2,ZVZV3,CORD,ID,
     &PRITISAK,PENALT,PRIT,IDPRIT)
     
       DO I=1,NDIM
        I1=I
        I2=I+NDIM
        I3=I+2*NDIM
        PPP=-PENALT*(PJ(1,I)*TT21(I1)+PJ(2,I)*TT21(I2)+PJ(3,I)*TT21(I3))
C        PPP=-PENALT*(SF1+SF2+SF3)

C        PRIT(1,NBREL)=PRITISAK
        
        PRIT(1,NBREL)=PRIT(1,NBREL)+PPP
        
c      WRITE(IIZLAZ,*)'PRIT(1,NBREL)= ',PRIT(1,NBREL),  NBREL
      
 450  CONTINUE

       ENDDO

505    CONTINUE



C END OF LOOP PER ELEMENT    
      END
C=========================================================================