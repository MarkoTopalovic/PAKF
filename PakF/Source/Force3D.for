C=========================================================================
      SUBROUTINE FORC3D(GNODE,ALEVO,DESNO,SILE,NEL,
     1ID,NGPSIL,MAXA,CORD,SKEF,PRIT,SPSIL,IBRGT,
     &ITFMAX,TABF,GUSM,CC,AKT,NETIP,IIZLAZ,AMI,INDAMI,BETA,TETAO,FB2,
     &FB3,NUMZAD,NPT,NDIM,MAXSIL,JEDN,NWK,NET,NDES,IDPRIT,IFORM,NPER,
     &NTABFT,PENALT,PRES,NSTAC,VVREME,TIME,INDEL,ZADVRE,NZAD,IUPWIN,
     &VMESH,IALE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
CE Subroutine FORC3D is used for calculation forces for 3D analysis
CE It is used global loop per elements
C

      DIMENSION GNODE(2,6,*),CORD(3,*)
      DIMENSION NEL(NDIM+1,*),ID(6,*),NGPSIL(8,*),MAXA(*)
      DIMENSION TABF(2,NTABFT,*),ITFMAX(*),INDEL(*),NZAD(3,*)

      DIMENSION PRIT(IDPRIT,*),SPSIL(NETIP,*),SKEF(NDES,*)
      DIMENSION ALEVO(*),DESNO(*),SILE(*),PRES(3,*),ZADVRE(*),VMESH(3,*)


      DIMENSION CK(21,3),LM2(92),TT21(92),TT210(92),PJ(3,21),TT21A(3*21)
      DIMENSION H(21),HP(8)
      DIMENSION AKVV1(21,21),AKMIV1(21,21),AMV2(21,21),AJV1V1(21,21)
      DIMENSION AJV1V2(21,21),AJV1V3(21,21),AJV2V1(21,21),AJV2V2(21,21)
      DIMENSION AJV2V3(21,21),AJV3V1(21,21),AJV3V2(21,21),AJV3V3(21,21)
      DIMENSION AKTV1(21,21),AKTV2(21,21),AKTV3(21,21),AKK(21,21)
      DIMENSION AKV1P(21,8),AKV2P(21,8),AKV3P(21,8),RS1(21),RS2(21)
      DIMENSION RS3(21),RB1(21),RB2(21),RB3(21),F92(92),SHEAR(3)
      DIMENSION XG(15),WGT(15),NREF(6)
      DIMENSION PENXX(10,10),PENXY(10,10),PENXZ(10,10)
      DIMENSION PENYX(10,10),PENYY(10,10),PENYZ(10,10)
      DIMENSION PENZX(10,10),PENZY(10,10),PENZZ(10,10)
      DIMENSION NID(8),NID1(8),NN(8)
      DIMENSION C(21,21),F92S(92)
	DIMENSION TOTF(3),TOTSF(3)
      DATA NREF/0,1,3,6,10,15/



 
      DATA WGT/2.,1.00000000000000,1.00000000000000,
     10.55555555555556,0.88888888888889,0.55555555555556,
     20.34785484513745,0.65214515486255,0.65214515486255,
     30.34785484513745,0.23692688505619,0.47862867049937,
     40.56888888888889,0.47862867049937,0.23692688505619/
 
      DATA XG/0.,-0.5773502691896,0.57735026918963,
     1-0.7745966692415,0.,0.77459666924148,-0.8611363115941,
     2-0.3399810435849,0.33998104358486,0.86113631159405,
     3-0.9061798459387,-0.5384693101057,0.,0.53846931010568,
     40.90617984593866/


C      R(3,1)=-0.7745966692415
C      R(3,2)=0.0
C      R(3,3)=0.77459666924148
C      S(3,1)=-0.7745966692415
C      S(3,2)=0.0
C      S(3,3)=0.77459666924148
C      W(3,1)=0.55555555555556
C      W(3,2)=0.88888888888889
C      W(3,3)=0.55555555555556


C      R(2,1)=0.57735026918963
C      R(2,2)=-0.5773502691896
C      S(2,1)=0.57735026918963
C      S(2,2)=-0.5773502691896
C      W(2,1)=1.000000000000000
C      W(2,2)=1.000000000000000

C      R(1,1)=0.0
C      S(1,1)=0.0
C      W(1,1)=2.0
         
C       PENALT=1.D9



      call INDELSSTRES(NEL,INDEL,NET,NPT,NDIM,ID)
      do n=1,numzad
       if (nzad(2,n).le.3) indel(nzad(1,n))=0
      enddo


      DO I=1,NPT
       VMESH(1,I)=0.D0
       VMESH(2,I)=0.D0
       VMESH(3,I)=0.D0
      ENDDO



C GLAVNA PETLJA PO ELEMENTIMA
      DO 400 NBREL=1,NET

      DO 125 I=1,92
      TT210(I)=0.D0
      LM2(I)=0
 125  TT21(I)=0.D0


C=========================================================================
      DO 130 KLM=1,NDIM
      CK(KLM,1)=CORD(1,NEL(KLM,NBREL))
      CK(KLM,2)=CORD(2,NEL(KLM,NBREL))
      CK(KLM,3)=CORD(3,NEL(KLM,NBREL))
      LM2(KLM)=ID(1,NEL(KLM,NBREL))
      LM2(KLM+NDIM)=ID(2,NEL(KLM,NBREL))
      LM2(KLM+2*NDIM)=ID(3,NEL(KLM,NBREL))
      IF (KLM.LE.8) LM2(KLM+3*NDIM)=ID(4,NEL(KLM,NBREL))
      LM2(KLM+3*NDIM+8)=ID(5,NEL(KLM,NBREL))
 130  CONTINUE
C=======================================================================
C=========================================================================
      IF (IALE.EQ.1) THEN
       DO I=1,NDIM
        NODE=NEL(I,NBREL)
        TT21A(I)=VMESH(1,NODE)
        TT21A(I+NDIM)=VMESH(2,NODE)
        TT21A(I+2*NDIM)=VMESH(3,NODE)
       ENDDO
      ENDIF
C=========================================================================
      DO 140 KLM=1,NDIM
      DO 135 NR=1,3
      TT21(KLM+(NR-1)*NDIM)=GNODE(2,NR,NEL(KLM,NBREL))
      TT210(KLM+(NR-1)*NDIM)=GNODE(1,NR,NEL(KLM,NBREL))
 135  CONTINUE
C      IF (KLM.LE.8.AND.ID(4,NEL(KLM,NBREL)).NE.0) THEN
      IF (KLM.LE.8) THEN
        TT21(KLM+3*NDIM)=GNODE(2,4,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM)=GNODE(1,4,NEL(KLM,NBREL))
      ENDIF
        TT21(KLM+3*NDIM+8)=GNODE(2,5,NEL(KLM,NBREL))
        TT210(KLM+3*NDIM+8)=GNODE(1,5,NEL(KLM,NBREL))
 140  CONTINUE
C=======================================================================
C      DO 160 NR=1,NDIM
C      DO 155 NZDT=1,NUMZAD
C      IF(NEL(NR,NBREL).EQ.NZAD(1,NZDT)) THEN
C       CALL TIMFUN (TABF,FK1,VVREME,ITFMAX(NZAD(3,NZDT)),NZAD(3,NZDT),
C     &NTABFT,IIZLAZ)
C      MESTO=NZAD(2,NZDT)
C       IF ((MESTO.EQ.4).AND.(NR.GT.8)) GOTO 160
C       IF (MESTO.EQ.5) THEN
C        TT21(3*NDIM+8+NR)=ZADVRE(NZDT)*FK1
C       ELSE
C        TT21((MESTO-1)*NDIM+NR)=ZADVRE(NZDT)*FK1
C       ENDIF
C      ENDIF  
C 155  CONTINUE
C 160  CONTINUE
C=======================================================================
C======================================================================= 
C RACUNANJE PRITISAKA 
C      IF (INDPR.EQ.1) THEN
C      K=0
C      DO 450 I=1,IBRGT-1
C      DO 450 J=1,IBRGT-1
C      K=K+1
C      CALL FNTERP(R(IBRGT-1,I),S(IBRGT-1,J),1)
C      PRIT(K,NBREL)=PRESS
C450   CONTINUE
C      GOTO 400
C      ENDIF
C
C=======================================================================
C      CALL WRR(TT21,3*NDIM+4,'T211')
      DO 163 K=1,NDIM
      DO 162 N=1,NDIM
  
      IF ((K.LE.8).AND.(N.LE.8)) THEN
       PENXX(K,N)=0.D0
       PENXY(K,N)=0.D0
       PENXZ(K,N)=0.D0
       PENYX(K,N)=0.D0
       PENYY(K,N)=0.D0
       PENYZ(K,N)=0.D0
       PENZX(K,N)=0.D0
       PENZY(K,N)=0.D0
       PENZZ(K,N)=0.D0
      ENDIF
      IF (N.LE.8) THEN
       AKV1P(K,N)=0.D0
       AKV2P(K,N)=0.D0
       AKV3P(K,N)=0.D0
      ENDIF
      AKVV1(K,N)=0.D0
      AKTV1(K,N)=0.D0
      AKTV2(K,N)=0.D0
      AKTV3(K,N)=0.D0
      AKMIV1(K,N)=0.D0
C      AKMIV3(K,N)=0.D0
C      A12(K,N)=0.
C      A21(K,N)=0.
      AMV2(K,N)=0.D0
      C(K,N)=0.D0
      AJV1V1(K,N)=0.D0
      AJV1V2(K,N)=0.D0
      AJV1V3(K,N)=0.D0
      AJV2V1(K,N)=0.D0
      AJV2V2(K,N)=0.D0
      AJV2V3(K,N)=0.D0
      AJV3V1(K,N)=0.D0
      AJV3V2(K,N)=0.D0
      AJV3V3(K,N)=0.D0
      AKK(K,N)=0.D0
C      AJK(K,N)=0.D0
  162 CONTINUE
  163 CONTINUE


C POVRSINSKE SILE I ZAPREMINSKE SILE:
      DO 199 I=1,NDIM
      RS1(I)=0.D0
      RS2(I)=0.D0
      RS3(I)=0.D0
      RB1(I)=0.D0
      RB2(I)=0.D0
      RB3(I)=0.D0
 199  CONTINUE

C===========================================================================

      NGAUSX=IBRGT
      NGAUSY=IBRGT
      NGAUSZ=IBRGT
C INTEGRACIJA U GAUSOVIM TACKAMA
C      DO 180 I=1,IBRGT
C      DO 170 J=1,IBRGT
C      CALL FNTERP(R(IBRGT,I),S(IBRGT,J),0)
C      WDT=W(IBRGT,I)*W(IBRGT,J)*DETJ

      ZAPRE=0.D0

      DO 180 NGX=1,NGAUSX
      JR=NREF(NGAUSX) + NGX
      R = XG(JR)
 
      DO 175 NGY=1,NGAUSY
      JS=NREF(NGAUSY) + NGY
      S = XG(JS)
 
      DO 170 NGZ=1,NGAUSZ
      JT=NREF(NGAUSZ) + NGZ
      T = XG(JT)
 
      WT=WGT(JR)*WGT(JS)*WGT(JT)
 
       CALL JACT(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN)

       IF (IALE.EQ.1) CALL ALEHV3(TT21,TT21A,HV1,HV2,HV3,NDIM,H)
 
      WDT=WT*DET1
      ZAPRE=ZAPRE+WDT
  
      DO 165 K=1,NDIM
      DO 164 N=1,NDIM
      AKVV1(K,N)=AKVV1(K,N)+WDT*((H(K)*HV1*PJ(1,N)+
     1H(K)*HV2*PJ(2,N)+H(K)*HV3*PJ(3,N))*GUSM)
      AKTV1(K,N)=AKTV1(K,N)+WDT*(H(K)*H(N)*ZVXT*GUSM*CC)
      AKTV2(K,N)=AKTV2(K,N)+WDT*(H(K)*H(N)*ZVYT*GUSM*CC)
      AKTV3(K,N)=AKTV3(K,N)+WDT*(H(K)*H(N)*ZVZT*GUSM*CC)
      AKMIV1(K,N)=AKMIV1(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AMI)
      AKK(K,N)=AKK(K,N)+WDT*((PJ(1,K)*PJ(1,N)+
     1PJ(2,K)*PJ(2,N)+PJ(3,K)*PJ(3,N))*AKT)
C      AKMIV3(K,N)=AKMIV1(K,N)
      AMV2(K,N)=AMV2(K,N)+WDT*H(K)*H(N)*GUSM
      AJV1V1(K,N)=AJV1V1(K,N)+H(K)*HXU*H(N)*GUSM*WDT
      AJV1V2(K,N)=AJV1V2(K,N)+H(K)*HYU*H(N)*GUSM*WDT
      AJV1V3(K,N)=AJV1V3(K,N)+H(K)*HZU*H(N)*GUSM*WDT
      AJV2V1(K,N)=AJV2V1(K,N)+H(K)*HXV*H(N)*GUSM*WDT
      AJV2V2(K,N)=AJV2V2(K,N)+H(K)*HYV*H(N)*GUSM*WDT
      AJV2V3(K,N)=AJV2V3(K,N)+H(K)*HZV*H(N)*GUSM*WDT
      AJV3V1(K,N)=AJV3V1(K,N)+H(K)*HXW*H(N)*GUSM*WDT
      AJV3V2(K,N)=AJV3V2(K,N)+H(K)*HYW*H(N)*GUSM*WDT
      AJV3V3(K,N)=AJV3V3(K,N)+H(K)*HZW*H(N)*GUSM*WDT
       IF (N.LE.8.AND.PENALT.LT.1.D0) THEN
       AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
       AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
       AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
      ENDIF
 164  CONTINUE
       RB2(K)=RB2(K)+H(K)*GUSM*FB2*(1.D0+BETA*TETAO)*WDT
       RB3(K)=RB3(K)+H(K)*GUSM*FB3*(1.D0+BETA*TETAO)*WDT
  165 CONTINUE

  170 CONTINUE
  175 CONTINUE
  180 CONTINUE



      IF (PENALT.GT.0.D0) THEN
       NGAUSX=IBRGT-1
       NGAUSY=IBRGT-1
       NGAUSZ=IBRGT-1

       DO  NGX=1,NGAUSX
       JR=NREF(NGAUSX) + NGX
       R = XG(JR)
 
       DO  NGY=1,NGAUSY
       JS=NREF(NGAUSY) + NGY
       S = XG(JS)
 
       DO  NGZ=1,NGAUSZ
       JT=NREF(NGAUSZ) + NGZ
       T = XG(JT)
 
       WT=WGT(JR)*WGT(JS)*WGT(JT)
 
       CALL JACT(R,S,T,DET1,CK,0,PJ,HV1,HV2,HV3,H,HP,TT21
     1,DET,SF1,SF2,SF3,NBREL,HXU,HYU,HZU,HXV,HYV,HZV,HXW,HYW,HZW,
     &AMI,NDIM,ISRPS,IIZLAZ,IUPWIN,NEL,SHEAR,ZVXT,ZVYT,ZVZT,DTEMPDN)

 
        WDT=WT*DET1
  
        DO  K=1,NDIM
        DO  N=1,NDIM
C      IF (N.LE.8.AND.PENALT.LT.1.D0) THEN
C       AKV1P(K,N)=AKV1P(K,N)+WDT*(-PJ(1,K)*HP(N))
C       AKV2P(K,N)=AKV2P(K,N)+WDT*(-PJ(2,K)*HP(N))
C       AKV3P(K,N)=AKV3P(K,N)+WDT*(-PJ(3,K)*HP(N))
C      ENDIF
C         IF (PENALT.GT.0.D0) THEN
          PENXX(K,N)=PENXX(K,N)+WDT*PJ(1,K)*PJ(1,N)*PENALT
          PENXY(K,N)=PENXY(K,N)+WDT*PJ(1,K)*PJ(2,N)*PENALT
          PENXZ(K,N)=PENXZ(K,N)+WDT*PJ(1,K)*PJ(3,N)*PENALT
          PENYX(K,N)=PENYX(K,N)+WDT*PJ(2,K)*PJ(1,N)*PENALT
          PENYY(K,N)=PENYY(K,N)+WDT*PJ(2,K)*PJ(2,N)*PENALT
          PENYZ(K,N)=PENYZ(K,N)+WDT*PJ(2,K)*PJ(3,N)*PENALT
          PENZX(K,N)=PENZX(K,N)+WDT*PJ(3,K)*PJ(1,N)*PENALT
          PENZY(K,N)=PENZY(K,N)+WDT*PJ(3,K)*PJ(2,N)*PENALT
          PENZZ(K,N)=PENZZ(K,N)+WDT*PJ(3,K)*PJ(3,N)*PENALT
C         ENDIF
        ENDDO
        ENDDO
       ENDDO
       ENDDO
       ENDDO
      ENDIF



C       WRITE (IIZLAZ,*)'ZAPREMINA=',ZAPRE
C===========================================================================
C REDUKOVANA INTEGRACIJA CLANOVA SA PENALTY FAKTOROM
C     IF (PENALT.GT.0.D0) THEN
C     DO 200 I=1,IBRGT-1
C     DO 195 J=1,IBRGT-1
C     CALL FNTERP(R(IBRGT-1,I),S(IBRGT-1,J),0)
C     WDT1=W(IBRGT-1,I)*W(IBRGT-1,J)*DETJ
C     DO 190 K=1,NDIM
C     DO 185 N=1,NDIM
C     AKMIV1(K,N)=AKMIV1(K,N)+PENALT*ZVHX(K)*ZVHX(N)*WDT1
C     AKMIV3(K,N)=AKMIV3(K,N)+PENALT*ZVHY(K)*ZVHY(N)*WDT1
C     A12(K,N)=A12(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT1
C     AJV2V3(K,N)=AJV2V3(K,N)+PENALT*ZVHX(K)*ZVHY(N)*WDT1
C     A21(K,N)=A21(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT1
C     AJV3V2(K,N)=AJV3V2(K,N)+PENALT*ZVHY(K)*ZVHX(N)*WDT1
C 185 CONTINUE
C 190 CONTINUE
C 195 CONTINUE
C 200 CONTINUE
C     ENDIF
C===========================================================================


C INCIJALIZACIJA MATRICE SKEF I F92
      DO 260 I=1,NDES
      DO 258 J=1,NDES
      SKEF(I,J)=0.D0
 258  CONTINUE
       F92(I)=0.D0
       F92S(I)=0.D0
 260  CONTINUE
C=========================================================================
C=========================================================================
c       goto 271
C LEVA STRANA:
      DO 270 I=1,NDIM
      DO 265 J=1,NDIM
       SKEF(I,J)=AKMIV1(I,J)+AKVV1(I,J)+AJV1V1(I,J)
       SKEF(I,J+NDIM)=AJV1V2(I,J)
       SKEF(I,J+2*NDIM)=AJV1V3(I,J)
      IF (J.LE.8) THEN
       SKEF(I,J+3*NDIM)=AKV1P(I,J)
       SKEF(I+NDIM,J+3*NDIM)=AKV2P(I,J)
       SKEF(I+2*NDIM,J+3*NDIM)=AKV3P(I,J)
      ENDIF
      IF (I.LE.8) THEN
       SKEF(I+3*NDIM,J)=AKV1P(J,I)
       SKEF(I+3*NDIM,J+NDIM)=AKV2P(J,I)
       SKEF(I+3*NDIM,J+2*NDIM)=AKV3P(J,I)
      ENDIF

       SKEF(I+NDIM,J)=AJV2V1(I,J)
       SKEF(I+NDIM,J+NDIM)=AKMIV1(I,J)+AKVV1(I,J)+AJV2V2(I,J)
       SKEF(I+NDIM,J+2*NDIM)=AJV2V3(I,J)

      SKEF(I+2*NDIM,J)=AJV3V1(I,J)
      SKEF(I+2*NDIM,J+NDIM)=AJV3V2(I,J)
      SKEF(I+2*NDIM,J+2*NDIM)=AKMIV1(I,J)+AKVV1(I,J)+AJV3V3(I,J)
      IF (NSTAC.EQ.0) THEN
       SKEF(I,J)=SKEF(I,J)+AMV2(I,J)/TIME
       SKEF(I+NDIM,J+NDIM)=SKEF(I+NDIM,J+NDIM)+AMV2(I,J)/TIME
       SKEF(I+2*NDIM,J+2*NDIM)=SKEF(I+2*NDIM,J+2*NDIM)+AMV2(I,J)/TIME
      ENDIF
 265  CONTINUE
 270  CONTINUE
 271  continue
       
       
       DO K=0,2
        DO I=1,NDIM
         II=I+K*NDIM
         DO J=1,NDES
          F92S(II)=F92S(II)+SKEF(II,J)*TT21(J)
         ENDDO
       ENDDO
      ENDDO
       
      IF (PENALT.GT.0.D0) THEN
      DO I=1,NDIM
        I1=I
        I2=I+NDIM
        I3=I+2*NDIM
       DO J=1,NDIM
        J1=J
        J2=J+NDIM
        J3=J+2*NDIM
        SKEF(I1,J1)=SKEF(I1,J1)+PENXX(I,J)
        SKEF(I1,J2)=SKEF(I1,J2)+PENXY(I,J)
        SKEF(I1,J3)=SKEF(I1,J3)+PENXZ(I,J)
        SKEF(I2,J1)=SKEF(I2,J1)+PENYX(I,J)
        SKEF(I2,J2)=SKEF(I2,J2)+PENYY(I,J)
        SKEF(I2,J3)=SKEF(I2,J3)+PENYZ(I,J)
        SKEF(I3,J1)=SKEF(I3,J1)+PENZX(I,J)
        SKEF(I3,J2)=SKEF(I3,J2)+PENZY(I,J)
        SKEF(I3,J3)=SKEF(I3,J3)+PENZZ(I,J)
       ENDDO
      ENDDO
      ENDIF

C RACUNANJE SILA INTERAKCIJE KOJIMA FLUID DELUJE NA ZIDOVE
c calculation of interaction forces between fluid and solid

       DO K=0,2
        DO I=1,NDIM
         II=I+K*NDIM
         DO J=1,NDES
          F92(II)=F92(II)+SKEF(II,J)*TT21(J)
         ENDDO
       ENDDO
      ENDDO

      IF (NSTAC.EQ.10) THEN
       DO K=0,2
        DO I=1,NDIM
          II=I+K*NDIM
         DO J=1,NDIM
          JJ=J+K*NDIM
c          F92(II)=F92(II)-AMV2(I,J)*TT210(JJ)/TIME
         ENDDO
        ENDDO
       ENDDO
       ENDIF



      
       DO 299 I=1,NDIM
C       IF (NID1(I).EQ.0) GOTO 299
       N=NEL(I,NBREL)
         P1=F92(I)
         P2=F92(I+NDIM)
         P3=F92(I+2*NDIM)
C        if (dabs(cord(2,n)).lt.1.d-5) P2=0.D0
        if (indel(n).eq.0) goto 299
        
         SPSIL(1,N)=SPSIL(1,N)-4.d0*P1/indel(n)
         SPSIL(2,N)=SPSIL(2,N)-4.d0*P2/indel(n)
         SPSIL(3,N)=SPSIL(3,N)-4.d0*P3/indel(n)
c         SPSIL(1,N)=SPSIL(1,N)-P1
c         SPSIL(2,N)=SPSIL(2,N)-P2
c         SPSIL(3,N)=SPSIL(3,N)-P3

         P1=F92S(I)
         P2=F92S(I+NDIM)
         P3=F92S(I+2*NDIM)

         VMESH(1,N)=VMESH(1,N)-4.d0*P1/indel(n)
         VMESH(2,N)=VMESH(2,N)-4.d0*P2/indel(n)
         VMESH(3,N)=VMESH(3,N)-4.d0*P3/indel(n)

299   CONTINUE      


C=======================================================================
C KRAJ PETLJE PO ELEMENTIMA
C=======================================================================
 400  CONTINUE

      TOTF(1)=0.D0
      TOTF(2)=0.D0
      TOTF(3)=0.D0
 
      DO I=1,NPT
        TOTF(1)=TOTF(1)+SPSIL(1,I)
        TOTF(2)=TOTF(2)+SPSIL(2,I)
        TOTF(3)=TOTF(3)+SPSIL(3,I)
	ENDDO
      write(iizlaz,*) 'total force Fx= ', TOTF(1)
      write(iizlaz,*) 'total force Fy= ', TOTF(2)
      write(iizlaz,*) 'total force Fz= ', TOTF(3)
c      write(49,*) 'total force= ',sqrt(TOTF(1)**2+TOTF(2)**2+TOTF(3)**2)
C End of subroutine
      TOTSF(1)=0.D0
      TOTSF(2)=0.D0
      TOTSF(3)=0.D0
 
      DO I=1,NPT
        TOTSF(1)=TOTSF(1)+VMESH(1,I)
        TOTSF(2)=TOTSF(2)+VMESH(2,I)
        TOTSF(3)=TOTSF(3)+VMESH(3,I)
	ENDDO
      write(iizlaz,*) 'total SHEAR force Fx= ', TOTSF(1)
      write(iizlaz,*) 'total SHEAR force Fy= ', TOTSF(2)
      write(iizlaz,*) 'total SHEAR force Fz= ', TOTSF(3)

      END
C=======================================================================
C=======================================================================    
C=========================================================================  
      SUBROUTINE FORCE3D(GNODE,NEL,ID,CORD,NETIP,IIZLAZ,SPSIL,
     &NUMZAD,NPT,NDIM,NET,PENALT,PRES,NSTAC,INDEL,ZADVRE,NZAD,VMESH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)                                    
C                                                                           
CE Subroutine FORC3D is used for calculation forces for 3D analysis         
CE It is used global loop per elements                                      
C                                                                           
                                                                            
      DIMENSION GNODE(2,6,*),CORD(3,*)                                      
      DIMENSION NEL(NDIM+1,*),ID(6,*)
      DIMENSION INDEL(*),NZAD(3,*)               
                                                                            
      DIMENSION SPSIL(NETIP,*)              
      DIMENSION PRES(3,*),ZADVRE(*),VMESH(3,*)    
                                                                            

      DIMENSION N(8),ITR(6,4),NWALL(6),NN(4),NID(8)
      DIMENSION TOTF(3),TOTSF(3),P(3)
      
      DO I=1,NPT
       VMESH(1,I)=0.D0
       VMESH(2,I)=0.D0
       VMESH(3,I)=0.D0
      ENDDO

      call INDELSSTRES(NEL,INDEL,NET,NPT,NDIM,ID)
      do I=1,numzad
       if (nzad(2,I).le.3) indel(nzad(1,I))=0
      enddo


      do 400 nbrel=1,net


	DO I=1,NDIM
	  N(I)=NEL(I,NBREL)
        NID(I)=0
      IF(ID(1,N(I)).EQ.0.AND.ID(2,N(I)).EQ.0.AND.ID(3,N(I)).EQ.0)
     &NID(I)=1
      ENDDO

      ITR(1,1)=N(8)
      ITR(1,2)=N(5)
      ITR(1,3)=N(1)
      ITR(1,4)=N(4)

      ITR(2,1)=N(7)
      ITR(2,2)=N(3)
      ITR(2,3)=N(2)
      ITR(2,4)=N(6)

      ITR(3,1)=N(6)
      ITR(3,2)=N(2)
      ITR(3,3)=N(1)
      ITR(3,4)=N(5)

      ITR(4,1)=N(7)
      ITR(4,2)=N(8)
      ITR(4,3)=N(4)
      ITR(4,4)=N(3)

      ITR(5,1)=N(3)
      ITR(5,2)=N(4)
      ITR(5,3)=N(1)
      ITR(5,4)=N(2)

      ITR(6,1)=N(7)
      ITR(6,2)=N(6)
      ITR(6,3)=N(5)
      ITR(6,4)=N(8)

      DO III=1,8
       do II=1,numzad                                                         
        if (nzad(2,II).le.3 .AND. NZAD(1,II).EQ.N(III)) NID(III)=0 
       enddo     
      ENDDO                                                            



      NWALL(1)=NID(1)*NID(4)*NID(5)*NID(8)
      NWALL(2)=NID(2)*NID(3)*NID(6)*NID(7)
      NWALL(3)=NID(1)*NID(2)*NID(5)*NID(6)
      NWALL(4)=NID(3)*NID(4)*NID(7)*NID(8)
      NWALL(5)=NID(1)*NID(2)*NID(3)*NID(4)
      NWALL(6)=NID(5)*NID(6)*NID(7)*NID(8)

      DO 300 NPOV=1,6
       IF (NWALL(NPOV).EQ.0) GOTO 300
        NN(1)=ITR(NPOV,1)
        NN(2)=ITR(NPOV,2)
        NN(3)=ITR(NPOV,3)
        NN(4)=ITR(NPOV,4)
	  DA=SURF4(CORD,NN(1),NN(2),NN(3),NN(4))
	   Pout=160000.d0
	   P1=GNODE(2,4,NN(1))+Pout
	   P2=GNODE(2,4,NN(2))+Pout
	   P3=GNODE(2,4,NN(3))+Pout
	   P4=GNODE(2,4,NN(4))+Pout

	   PRESSURE=0.25*(P1+P2+P3+P4)
	   P2X=CORD(1,NN(3))-CORD(1,NN(1))
	   P2Y=CORD(2,NN(3))-CORD(2,NN(1))
	   P2Z=CORD(3,NN(3))-CORD(3,NN(1))

	   P1X=CORD(1,NN(2))-CORD(1,NN(1))
	   P1Y=CORD(2,NN(2))-CORD(2,NN(1))
	   P1Z=CORD(3,NN(2))-CORD(3,NN(1))
	   
          P(1)=P1Y*P2Z-P1Z*P2Y
          P(2)=P1Z*P2X-P1X*P2Z
          P(3)=P1X*P2Y-P1Y*P2X
          
          PP=DSQRT(P(1)**2+P(2)**2+P(3)**2)
          P(1)=P(1)/PP
          P(2)=P(2)/PP
          P(3)=P(3)/PP
        
         	   
         DO K=1,3
  	 TAU=0.25D0*(PRES(K,NN(1))+PRES(K,NN(2))+PRES(K,NN(3))+PRES(K,NN(4)))
 	   PRESS=PRESSURE*P(K)
 	   

 	   SPSIL(K,NN(1))=SPSIL(K,NN(1))+0.25D0*(PRESS+TAU)*DA
 	   SPSIL(K,NN(2))=SPSIL(K,NN(2))+0.25D0*(PRESS+TAU)*DA
 	   SPSIL(K,NN(3))=SPSIL(K,NN(3))+0.25D0*(PRESS+TAU)*DA
 	   SPSIL(K,NN(4))=SPSIL(K,NN(4))+0.25D0*(PRESS+TAU)*DA

c 	   SPSIL(K,NN(1))=SPSIL(K,NN(1))+(PRESS+TAU)*DA/INDEL(NN(1))
c 	   SPSIL(K,NN(2))=SPSIL(K,NN(2))+(PRESS+TAU)*DA/INDEL(NN(2))
c 	   SPSIL(K,NN(3))=SPSIL(K,NN(3))+(PRESS+TAU)*DA/INDEL(NN(3))
c 	   SPSIL(K,NN(4))=SPSIL(K,NN(4))+(PRESS+TAU)*DA/INDEL(NN(4))
 	   
c 	   SPSIL(K,NN(1))=SPSIL(K,NN(1))+0.25D0*(TAU)*DA
c 	   SPSIL(K,NN(2))=SPSIL(K,NN(2))+0.25D0*(TAU)*DA
c 	   SPSIL(K,NN(3))=SPSIL(K,NN(3))+0.25D0*(TAU)*DA
c 	   SPSIL(K,NN(4))=SPSIL(K,NN(4))+0.25D0*(TAU)*DA

 	   VMESH(K,NN(1))=VMESH(K,NN(1))+0.25D0*TAU*DA
 	   VMESH(K,NN(2))=VMESH(K,NN(2))+0.25D0*TAU*DA
 	   VMESH(K,NN(3))=VMESH(K,NN(3))+0.25D0*TAU*DA
 	   VMESH(K,NN(4))=VMESH(K,NN(4))+0.25D0*TAU*DA


c 	   VMESH(K,NN(1))=VMESH(K,NN(1))+TAU*DA/INDEL(NN(1))
c 	   VMESH(K,NN(2))=VMESH(K,NN(2))+TAU*DA/INDEL(NN(2))
c 	   VMESH(K,NN(3))=VMESH(K,NN(3))+TAU*DA/INDEL(NN(3))
c 	   VMESH(K,NN(4))=VMESH(K,NN(4))+TAU*DA/INDEL(NN(4))
 	   ENDDO
 	   
300   continue
400   continue       

                                                                            
                                                                            
      TOTF(1)=0.D0                                                          
      TOTF(2)=0.D0                                                          
      TOTF(3)=0.D0                                                          
                                                                            
      DO I=1,NPT                                                            
        TOTF(1)=TOTF(1)+SPSIL(1,I)                                          
        TOTF(2)=TOTF(2)+SPSIL(2,I)                                          
        TOTF(3)=TOTF(3)+SPSIL(3,I)                                          
	ENDDO                                                                     
      write(iizlaz,*) 'total force Fx= ', TOTF(1)                           
      write(iizlaz,*) 'total force Fy= ', TOTF(2)                           
      write(iizlaz,*) 'total force Fz= ', TOTF(3)                           

      TOTSF(1)=0.D0                                                          
      TOTSF(2)=0.D0                                                          
      TOTSF(3)=0.D0                                                          
                                                                            
      DO I=1,NPT                                                            
        TOTSF(1)=TOTSF(1)+VMESH(1,I)                                          
        TOTSF(2)=TOTSF(2)+VMESH(2,I)                                          
        TOTSF(3)=TOTSF(3)+VMESH(3,I)                                          
	ENDDO                                                                     
      write(iizlaz,*) 'total SHEAR force Fx= ', TOTSF(1)                           
      write(iizlaz,*) 'total SHEAR force Fy= ', TOTSF(2)                           
      write(iizlaz,*) 'total SHEAR force Fz= ', TOTSF(3)                           
C End of subroutine                                                         
                                                                            
      END                                                                   
C=======================================================================    
C=========================================================================  
C==========================================================================
      FUNCTION SURF4(CORD,N1,N2,N3,N4)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
	DIMENSION CORD(3,*)

	  X1=CORD(1,N1) 
	  Y1=CORD(2,N1) 
	  Z1=CORD(3,N1) 
	  
	  X2=CORD(1,N2) 
	  Y2=CORD(2,N2) 
	  Z2=CORD(3,N2) 
	  
	  X3=CORD(1,N3) 
	  Y3=CORD(2,N3) 
	  Z3=CORD(3,N3) 
	  
	  X4=CORD(1,N4) 
	  Y4=CORD(2,N4) 
	  Z4=CORD(3,N4) 

	  A=DSQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
	  B=DSQRT((X2-X3)**2+(Y2-Y3)**2+(Z2-Z3)**2)
	  C=DSQRT((X1-X3)**2+(Y1-Y3)**2+(Z1-Z3)**2)
	  S=(A+B+C)*0.5D0
	  P1=DSQRT(DABS(S*(S-A)*(S-B)*(S-C)))

	  A=DSQRT((X2-X3)**2+(Y2-Y3)**2+(Z2-Z3)**2)
	  B=DSQRT((X3-X4)**2+(Y3-Y4)**2+(Z3-Z4)**2)
	  C=DSQRT((X2-X4)**2+(Y2-Y4)**2+(Z2-Z4)**2)
	  S=(A+B+C)*0.5D0
	  P2=DSQRT(DABS(S*(S-A)*(S-B)*(S-C)))

        SURF4=P1+P2
       
	END
C==========================================================================
C==========================================================================
