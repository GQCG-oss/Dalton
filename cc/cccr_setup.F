!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
#include <single.h>
C
c /* deck cccr_setup */
*=====================================================================*
      SUBROUTINE CCCR_SETUP(MXTRAN2, MXVEC2, MXTRAN3, MXVEC1,
     &                      I0HTRAN, I0HDOTS, N0HTRAN,
     &                      I0GTRAN, I0GDOTS, N0GTRAN,
     &                      IAGTRAN, IAGDOTS, NAGTRAN,
     &                      I0FTRAN, I0FDOTS, N0FTRAN,
     &                      IAFTRAN, IAFDOTS, NAFTRAN,
     &                      I0FATRAN,I0FADOTS,N0FATRAN,
     &                      IAFATRAN,IAFADOTS,NAFATRAN,
     &                      IAEATRAN,IAEADOTS,NAEATRAN,
     &                      IXTRAN,  IXDOTS,  NXTRAN,   
     &                      IOTRAN,  IODOTS,  NOTRAN,
     &                      ILTRAN,  ILDOTS,  NLTRAN    )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CCCR section
*                - list of H^0 matrix transformations 
*                - list of G^0 matrix transformations 
*                - list of G^A matrix transformations 
*                - list of F^0 matrix transformations 
*                - list of F^A matrix transformations 
*                - list of F^0{O} matrix transformations 
*                - list of F^A{O} matrix transformations 
*                - list of ETA^A{O} vector calculations 
*                - list of chi vector dot products
*                - list of xksi vector dot products
*                - list of L2 x O2 vector dot products
*
*     Written by Christof Haettig, februar 1997.
*     turn over rule options (USE_L2BC, USE_LBCD) added in april 1998
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cccrinf.h>
#include <ccroper.h>
#include <cccperm.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCCR_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER MXVEC2, MXTRAN2, MXVEC1, MXTRAN3

      INTEGER I0HTRAN(MXDIM_HTRAN,MXTRAN3)
      INTEGER I0HDOTS(MXVEC1,MXTRAN3)

      INTEGER I0GTRAN(MXDIM_GTRAN,MXTRAN2)
      INTEGER I0GDOTS(MXVEC2,MXTRAN2)

      INTEGER IAGTRAN(MXDIM_GTRAN,MXTRAN3)
      INTEGER IAGDOTS(MXVEC1,MXTRAN3)

      INTEGER I0FTRAN(MXDIM_FTRAN,MXTRAN2)
      INTEGER I0FDOTS(MXVEC2,MXTRAN2)
  
      INTEGER IAFTRAN(MXDIM_FTRAN,MXTRAN2)
      INTEGER IAFDOTS(MXVEC2,MXTRAN2)
  
      INTEGER I0FATRAN(MXDIM_FATRAN,MXTRAN2)
      INTEGER I0FADOTS(MXVEC2,MXTRAN2)

      INTEGER IAFATRAN(MXDIM_FATRAN,MXTRAN3)
      INTEGER IAFADOTS(MXVEC1,MXTRAN3)

C     INTEGER IAEATRAN(3,MXTRAN2)
      INTEGER IAEATRAN(MXDIM_XEVEC,MXTRAN2)
      INTEGER IAEADOTS(MXVEC2,MXTRAN2)

      INTEGER IXTRAN(MXTRAN2)
      INTEGER IXDOTS(MXVEC2,MXTRAN2)

      INTEGER ILTRAN(MXTRAN2)
      INTEGER ILDOTS(MXVEC2,MXTRAN2)

      INTEGER IOTRAN(MXTRAN3)
      INTEGER IODOTS(MXVEC1,MXTRAN3)

      INTEGER N0HTRAN, N0GTRAN, N0FTRAN, N0FATRAN, NXTRAN, NOTRAN
      INTEGER          NAGTRAN, NAFTRAN, NAFATRAN, NAEATRAN, NLTRAN
      INTEGER NCRRESF

      INTEGER IVEC, ITRAN, I
      INTEGER ISYML, ISYM1, ISYM2, ISYM3
      INTEGER IFREQ, IOPER
      INTEGER P, ISIGN
      INTEGER MXV0H, MXV0G, MXVAG, MXV0F, MXVAF, MXV0FA, MXVAFA
      INTEGER MXVAEA, MXX, MXO, MXL

#if defined (SYS_CRAY)
      REAL SIGN
#else
      DOUBLE PRECISION SIGN
#endif

      INTEGER IOP(4), ISY(4), IZT(4), IR1(4), IR2(6), IX2(6)
      INTEGER IO3(4), IO2(6), IL2(6), IE1(4), IKAP(4)

* external functions:
      INTEGER IR3TAMP
      INTEGER IR2TAMP
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER IL2ZETA
      INTEGER IRHSR2
      INTEGER ICHI2
      INTEGER IETA1
      INTEGER IRHSR3


*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      N0HTRAN  = 0
      N0GTRAN  = 0
      NAGTRAN  = 0
      N0FTRAN  = 0
      NAFTRAN  = 0
      N0FATRAN = 0
      NAFATRAN = 0
      NAEATRAN = 0
      NXTRAN   = 0
      NOTRAN   = 0
      NLTRAN   = 0
      NCRRESF  = 0
 
      MXV0H    = 0
      MXV0G    = 0
      MXVAG    = 0
      MXV0F    = 0
      MXVAF    = 0
      MXV0FA   = 0
      MXVAFA   = 0
      MXVAEA   = 0
      MXX      = 0
      MXO      = 0
      MXL      = 0

      DO ITRAN = 1, MXTRAN2
        IAEATRAN(3,ITRAN) = -1
      END DO

*---------------------------------------------------------------------*
* start loop over all requested second hyperpolarizabilities
*---------------------------------------------------------------------*
 
      DO IOPER = 1, NCROPER
        IOP(A) = IACROP(IOPER)
        IOP(B) = IBCROP(IOPER)
        IOP(C) = ICCROP(IOPER)
        IOP(D) = IDCROP(IOPER)

        IKAP(A)= 0
        IKAP(B)= 0
        IKAP(C)= 0
        IKAP(D)= 0

        ISY(A) = ISYOPR(IOP(A))
        ISY(B) = ISYOPR(IOP(B))
        ISY(C) = ISYOPR(IOP(C))
        ISY(D) = ISYOPR(IOP(D))

      IF (MULD2H(ISY(A),ISY(B)).EQ.MULD2H(ISY(C),ISY(D))) THEN

      DO IFREQ = 1, NCRFREQ

        NCRRESF = NCRRESF + 1
      
      DO ISIGN = 1, -1, -2
        SIGN = DBLE(ISIGN)

        IE1(A) = IETA1(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYML)
        IE1(B) = IETA1(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYML)
        IE1(C) = IETA1(LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYML)
        IE1(D) = IETA1(LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYML)

        IZT(A) = IL1ZETA(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYML)
        IZT(B) = IL1ZETA(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYML)
        IZT(C) = IL1ZETA(LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYML)
        IZT(D) = IL1ZETA(LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYML)

        IR1(A) = IR1TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYML)
        IR1(B) = IR1TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYML)
        IR1(C) = IR1TAMP(LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYML)
        IR1(D) = IR1TAMP(LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYML)


        IF (USE_LBCD) THEN ! L2(BC,BD,CD) instead of R2(AD,AC,AB)
         IL2(BC)=IL2ZETA(LBLOPR(IOP(B)),        SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),        SIGN*CCRFR(IFREQ),ISYM2)
         IL2(BD)=IL2ZETA(LBLOPR(IOP(B)),        SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),        SIGN*DCRFR(IFREQ),ISYM2)
         IL2(CD)=IL2ZETA(LBLOPR(IOP(C)),        SIGN*CCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),        SIGN*DCRFR(IFREQ),ISYM2)
C        IX2(AD)=ICHI2(  LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
C    &                   LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
         IO2(AD)=IRHSR2( LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
C        IX2(AC)=ICHI2(  LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
C    &                   LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
         IO2(AC)=IRHSR2( LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
C        IX2(AB)=ICHI2(  LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
C    &                   LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM2)
         IO2(AB)=IRHSR2( LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM2)
        ELSE IF (USE_L2BC) THEN  ! use L2(BC) instead of R2(AD)
         IR2(AB)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM2)
         IR2(AC)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
         IL2(BC)=IL2ZETA(LBLOPR(IOP(B)),        SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),        SIGN*CCRFR(IFREQ),ISYM2)

C        IX2(AD)=ICHI2(  LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
C    &                   LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
         IO2(AD)=IRHSR2( LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
        ELSE  ! 2n+1/2n+2 rule formula symmetric in A,B,C,D
         IR2(AB)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM2)
         IR2(AC)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
         IR2(AD)=IR2TAMP(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
        END IF

        IR2(BC)=IR2TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
        IR2(BD)=IR2TAMP(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
        IR2(CD)=IR2TAMP(LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM1,
     &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)

        IF (L_USE_CHI2) THEN
c         IX2(AB)=ICHI2(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM2)
c         IX2(AC)=ICHI2(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
c         IX2(AD)=ICHI2(LBLOPR(IOP(A)),.FALSE.,SIGN*ACRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
c         IX2(BC)=ICHI2(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM2)
c         IX2(BD)=ICHI2(LBLOPR(IOP(B)),.FALSE.,SIGN*BCRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)
c         IX2(CD)=ICHI2(LBLOPR(IOP(C)),.FALSE.,SIGN*CCRFR(IFREQ),ISYM1,
c    &                  LBLOPR(IOP(D)),.FALSE.,SIGN*DCRFR(IFREQ),ISYM2)

         IX2(AB)=IL2ZETA(LBLOPR(IOP(A)),        SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(B)),        SIGN*BCRFR(IFREQ),ISYM2)
         IX2(AC)=IL2ZETA(LBLOPR(IOP(A)),        SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),        SIGN*CCRFR(IFREQ),ISYM2)
         IX2(AD)=IL2ZETA(LBLOPR(IOP(A)),        SIGN*ACRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),        SIGN*DCRFR(IFREQ),ISYM2)
         IX2(BC)=IL2ZETA(LBLOPR(IOP(B)),        SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(C)),        SIGN*CCRFR(IFREQ),ISYM2)
         IX2(BD)=IL2ZETA(LBLOPR(IOP(B)),        SIGN*BCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),        SIGN*DCRFR(IFREQ),ISYM2)
         IX2(CD)=IL2ZETA(LBLOPR(IOP(C)),        SIGN*CCRFR(IFREQ),ISYM1,
     &                   LBLOPR(IOP(D)),        SIGN*DCRFR(IFREQ),ISYM2)
        END IF

        IF (L_USE_XKS3) THEN
          IO3(ABC) = IR3TAMP(LBLOPR(IOP(A)),SIGN*ACRFR(IFREQ),ISYM1,
     &                      LBLOPR(IOP(B)),SIGN*BCRFR(IFREQ),ISYM2,
     &                      LBLOPR(IOP(C)),SIGN*CCRFR(IFREQ),ISYM3)
          IO3(ABD) = IR3TAMP(LBLOPR(IOP(A)),SIGN*ACRFR(IFREQ),ISYM1,
     &                      LBLOPR(IOP(B)),SIGN*BCRFR(IFREQ),ISYM2,
     &                      LBLOPR(IOP(D)),SIGN*DCRFR(IFREQ),ISYM3)
          IO3(ACD) = IR3TAMP(LBLOPR(IOP(A)),SIGN*ACRFR(IFREQ),ISYM1,
     &                      LBLOPR(IOP(C)),SIGN*CCRFR(IFREQ),ISYM2,
     &                      LBLOPR(IOP(D)),SIGN*DCRFR(IFREQ),ISYM3)
          IO3(BCD) = IR3TAMP(LBLOPR(IOP(B)),SIGN*BCRFR(IFREQ),ISYM1,
     &                      LBLOPR(IOP(C)),SIGN*CCRFR(IFREQ),ISYM2,
     &                      LBLOPR(IOP(D)),SIGN*DCRFR(IFREQ),ISYM3)
        END IF

*---------------------------------------------------------------------*
* set up list of H^0 matrix transformations, 1 permutation
*---------------------------------------------------------------------*
        CALL CC_SETH1111(I0HTRAN,I0HDOTS,MXTRAN3,MXVEC1,
     &                   0,IR1(A),IR1(B),IR1(C),IR1(D),ITRAN,IVEC)
        N0HTRAN = MAX(N0HTRAN,ITRAN)
        MXV0H   = MAX(MXV0H,IVEC) 

*---------------------------------------------------------------------*
* set up list of G^0 matrix transformations, 6 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        IF ((USE_L2BC.AND.P.EQ.4) .OR. (USE_LBCD.AND.P.GT.3) ) THEN
          CONTINUE
        ELSE 
          CALL CC_SETG112(I0GTRAN,I0GDOTS,MXTRAN2,MXVEC2,   
     &                0,IR1(I1(P)),IR1(I2(P)),IR2(IP(P)),ITRAN,IVEC)
          N0GTRAN = MAX(N0GTRAN,ITRAN)
          MXV0G   = MAX(MXV0G,IVEC) 
        END IF
      END DO

*---------------------------------------------------------------------*
* set up list of G^A matrix transformations, 4 permutations
*---------------------------------------------------------------------*
        CALL CCQR_SETG(IAGTRAN,IAGDOTS,MXTRAN3,MXVEC1,  
     &                 IZT(A),IR1(B),IR1(C),IR1(D),ITRAN,IVEC)
        NAGTRAN = MAX(NAGTRAN,ITRAN)
        MXVAG   = MAX(MXVAG,IVEC) 

        CALL CCQR_SETG(IAGTRAN,IAGDOTS,MXTRAN3,MXVEC1,
     &                 IZT(B),IR1(A),IR1(C),IR1(D),ITRAN,IVEC)
        NAGTRAN = MAX(NAGTRAN,ITRAN)
        MXVAG   = MAX(MXVAG,IVEC) 

        CALL CCQR_SETG(IAGTRAN,IAGDOTS,MXTRAN3,MXVEC1,
     &                 IZT(C),IR1(B),IR1(A),IR1(D),ITRAN,IVEC)
        NAGTRAN = MAX(NAGTRAN,ITRAN)
        MXVAG   = MAX(MXVAG,IVEC) 

        CALL CCQR_SETG(IAGTRAN,IAGDOTS,MXTRAN3,MXVEC1,
     &                 IZT(D),IR1(B),IR1(C),IR1(A),ITRAN,IVEC)
        NAGTRAN = MAX(NAGTRAN,ITRAN)
        MXVAG   = MAX(MXVAG,IVEC) 

*---------------------------------------------------------------------*
* set up list of F^0 matrix transformations, 3 permutations
*---------------------------------------------------------------------*
        IF (.NOT. USE_LBCD) THEN
          CALL CCQR_SETF(I0FTRAN,I0FDOTS,MXTRAN2,MXVEC2,
     &                   0,IR2(AB),IR2(CD),ITRAN,IVEC)
          N0FTRAN = MAX(N0FTRAN,ITRAN)
          MXV0F   = MAX(MXV0F,IVEC) 
        END IF

        IF (.NOT. USE_LBCD) THEN
          CALL CCQR_SETF(I0FTRAN,I0FDOTS,MXTRAN2,MXVEC2,
     &                   0,IR2(AC),IR2(BD),ITRAN,IVEC)
          N0FTRAN = MAX(N0FTRAN,ITRAN)
          MXV0F   = MAX(MXV0F,IVEC) 
        END IF

        IF (.NOT. (USE_LBCD .OR. USE_L2BC)) THEN
          CALL CCQR_SETF(I0FTRAN,I0FDOTS,MXTRAN2,MXVEC2,
     &                   0,IR2(AD),IR2(BC),ITRAN,IVEC)
          N0FTRAN = MAX(N0FTRAN,ITRAN)
          MXV0F   = MAX(MXV0F,IVEC) 
        END IF

*---------------------------------------------------------------------*
* set up list of F^A matrix transformations, 12 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        IF ((USE_L2BC.AND.P.EQ.4) .OR. (USE_LBCD.AND.P.GT.3) ) THEN
          CONTINUE
        ELSE 
          CALL CC_SETF12(IAFTRAN,IAFDOTS,MXTRAN2,MXVEC2, ! 1 x 2 x 3,4
     &                   IZT(I1(P)),IR1(I2(P)),IR2(IP(P)),ITRAN,IVEC)
          NAFTRAN = MAX(NAFTRAN,ITRAN)
          MXVAF   = MAX(MXVAF,IVEC) 

          CALL CC_SETF12(IAFTRAN,IAFDOTS,MXTRAN2,MXVEC2, ! 2 x 1 x 3,4
     &                   IZT(I2(P)),IR1(I1(P)),IR2(IP(P)),ITRAN,IVEC)
          NAFTRAN = MAX(NAFTRAN,ITRAN)
          MXVAF   = MAX(MXVAF,IVEC) 
        END IF
      END DO

*---------------------------------------------------------------------*
* set up list of F^0{O} matrix transformations, 12 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        IF ((USE_L2BC.AND.P.EQ.4) .OR. (USE_LBCD.AND.P.GT.3) ) THEN
          CONTINUE
        ELSE 
          CALL CC_SETFA12(I0FATRAN,I0FADOTS,MXTRAN2,MXVEC2, ! 1x2x3,4
     &                 0,IOP(I1(P)),IR1(I2(P)),IR2(IP(P)),ITRAN,IVEC)
          N0FATRAN = MAX(N0FATRAN,ITRAN)
          MXV0FA   = MAX(MXV0FA,IVEC) 

          CALL CC_SETFA12(I0FATRAN,I0FADOTS,MXTRAN2,MXVEC2, ! 2x1x3,4
     &                 0,IOP(I2(P)),IR1(I1(P)),IR2(IP(P)),ITRAN,IVEC)
          N0FATRAN = MAX(N0FATRAN,ITRAN)
          MXV0FA   = MAX(MXV0FA,IVEC) 
        END IF
      END DO

*---------------------------------------------------------------------*
* set up list of F^A{O} matrix transformations, 12 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        CALL CCQR_SETFA(IAFATRAN,IAFADOTS,MXTRAN3,MXVEC1, ! 1 x 2 x 3,4
     &      IZT(I1(P)),IOP(I2(P)),IR1(I3(P)),IR1(I4(P)),ITRAN,IVEC)
        NAFATRAN = MAX(NAFATRAN,ITRAN)
        MXVAFA   = MAX(MXVAFA,IVEC) 

        CALL CCQR_SETFA(IAFATRAN,IAFADOTS,MXTRAN3,MXVEC1, ! 2 x 1 x 3,4
     &      IZT(I2(P)),IOP(I1(P)),IR1(I3(P)),IR1(I4(P)),ITRAN,IVEC)
        NAFATRAN = MAX(NAFATRAN,ITRAN)
        MXVAFA   = MAX(MXVAFA,IVEC) 
      END DO

*---------------------------------------------------------------------*
* set up list of ETA{O} vector calculations, 12 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        IF ((USE_L2BC.AND.P.EQ.4) .OR. (USE_LBCD.AND.P.GT.3) ) THEN
          CONTINUE
        ELSE 
C         CALL CCQR_SETEA(IAEATRAN,IAEADOTS,MXTRAN2,MXVEC2, ! 1x2x3,4
C    &                    IZT(I1(P)),IOP(I2(P)),IR2(IP(P)),ITRAN,IVEC)
          CALL CC_SETXE('Eta',IAEATRAN,IAEADOTS,MXTRAN2,MXVEC2,!1x2x3,4
     &                  IZT(I1(P)),IOP(I2(P)),IKAP(I2(P)),0,0,0,
     &                  IR2(IP(P)),ITRAN,IVEC)
          NAEATRAN = MAX(NAEATRAN,ITRAN)
          MXVAEA   = MAX(MXVAEA,IVEC) 

C         CALL CCQR_SETEA(IAEATRAN,IAEADOTS,MXTRAN2,MXVEC2, ! 2x1x3,4
C    &                    IZT(I2(P)),IOP(I1(P)),IR2(IP(P)),ITRAN,IVEC)
          CALL CC_SETXE('Eta',IAEATRAN,IAEADOTS,MXTRAN2,MXVEC2,!2x1x3,4
     &                  IZT(I2(P)),IOP(I1(P)),IKAP(I1(P)),0,0,0,
     &                  IR2(IP(P)),ITRAN,IVEC)
          NAEATRAN = MAX(NAEATRAN,ITRAN)
          MXVAEA   = MAX(MXVAEA,IVEC) 
        END IF
      END DO

*---------------------------------------------------------------------*
* set up list of L2 x O2 vector dot products, max. 3 permutations
*---------------------------------------------------------------------*
      IF (USE_LBCD .OR. USE_L2BC) THEN
        CALL CC_SETDOT(ILTRAN,ILDOTS,MXTRAN2,MXVEC2,
     &                 IL2(BC),IO2(AD),ITRAN,IVEC)
        NLTRAN = MAX(NLTRAN,ITRAN)
        MXL    = MAX(MXL,IVEC) 
      END IF

      IF (USE_LBCD) THEN
        CALL CC_SETDOT(ILTRAN,ILDOTS,MXTRAN2,MXVEC2,
     &                 IL2(BD),IO2(AC),ITRAN,IVEC)
        NLTRAN = MAX(NLTRAN,ITRAN)
        MXL    = MAX(MXL,IVEC) 

        CALL CC_SETDOT(ILTRAN,ILDOTS,MXTRAN2,MXVEC2,
     &                 IL2(CD),IO2(AB),ITRAN,IVEC)
        NLTRAN = MAX(NLTRAN,ITRAN)
        MXL    = MAX(MXL,IVEC) 
      END IF
*---------------------------------------------------------------------*
* set up list of chi vector dot products, 6 permutations
*---------------------------------------------------------------------*
      IF (L_USE_CHI2) THEN 
        CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                 IX2(AB),IR2(CD),ITRAN,IVEC)
        NXTRAN = MAX(NXTRAN,ITRAN)
        MXX    = MAX(MXX,IVEC) 

        CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                 IX2(AC),IR2(BD),ITRAN,IVEC)
        NXTRAN = MAX(NXTRAN,ITRAN)
        MXX    = MAX(MXX,IVEC) 

        CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                 IX2(AD),IR2(BC),ITRAN,IVEC)
        NXTRAN = MAX(NXTRAN,ITRAN)
        MXX    = MAX(MXX,IVEC) 

        CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                 IX2(CD),IR2(AB),ITRAN,IVEC)
        NXTRAN = MAX(NXTRAN,ITRAN)
        MXX    = MAX(MXX,IVEC) 

        CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                 IX2(BD),IR2(AC),ITRAN,IVEC)
        NXTRAN = MAX(NXTRAN,ITRAN)
        MXX    = MAX(MXX,IVEC) 

        CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                 IX2(BC),IR2(AD),ITRAN,IVEC)
        NXTRAN = MAX(NXTRAN,ITRAN)
        MXX    = MAX(MXX,IVEC) 
      END IF

*---------------------------------------------------------------------*
* set up list of Xksi3 vector dot products, 4 permutations
*---------------------------------------------------------------------*
      IF (L_USE_XKS3) THEN 
        CALL CC_SETDOT(IOTRAN,IODOTS,MXTRAN3,MXVEC1,
     &                 IO3(ABC),IZT(D),ITRAN,IVEC)
        NOTRAN = MAX(NOTRAN,ITRAN)
        MXO    = MAX(MXO,IVEC) 

        CALL CC_SETDOT(IOTRAN,IODOTS,MXTRAN3,MXVEC1,
     &                 IO3(ABD),IZT(C),ITRAN,IVEC)
        NOTRAN = MAX(NOTRAN,ITRAN)
        MXO    = MAX(MXO,IVEC) 

        CALL CC_SETDOT(IOTRAN,IODOTS,MXTRAN3,MXVEC1,
     &                 IO3(ACD),IZT(B),ITRAN,IVEC)
        NOTRAN = MAX(NOTRAN,ITRAN)
        MXO    = MAX(MXO,IVEC) 

        CALL CC_SETDOT(IOTRAN,IODOTS,MXTRAN3,MXVEC1,
     &                 IO3(BCD),IZT(A),ITRAN,IVEC)
        NOTRAN = MAX(NOTRAN,ITRAN)
        MXO    = MAX(MXO,IVEC) 
      END IF


*---------------------------------------------------------------------*
* end loop over all requested hyperpolarizabilities
*---------------------------------------------------------------------*
      END DO
      END DO
      END IF
      END DO

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      WRITE(LUPRI,'(/,/3X,A,I4,A)') 'For the requested',NCRRESF,
     &      ' cubic response functions '
      WRITE(LUPRI,'((8X,A,I3,A))') 
     &   ' - ',N0HTRAN,  ' H^0 matrix transformations ',
     &   ' - ',N0GTRAN,  ' G^0 matrix transformations ',
     &   ' - ',NAGTRAN,  ' G^A matrix transformations ',
     &   ' - ',N0FTRAN,  ' F^0 matrix transformations ',
     &   ' - ',NAFTRAN,  ' F^A matrix transformations ',
     &   ' - ',N0FATRAN, ' F^0{O} matrix transformations ',
     &   ' - ',NAFATRAN, ' F^A{O} matrix transformations ',
     &   ' - ',NAEATRAN, ' ETA^A{O} vector calculations ' 
      WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'


* K^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of H^0 matrix transformations:'
        DO ITRAN = 1, N0HTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I0HTRAN(I,ITRAN),I=1,4),(I0HDOTS(I,ITRAN),I=1,MXV0H)
        END DO
        WRITE (LUPRI,*)
      END IF

* G^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G^0 matrix transformations:'
        DO ITRAN = 1, N0GTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (I0GTRAN(I,ITRAN),I=1,3),(I0GDOTS(I,ITRAN),I=1,MXV0G)
        END DO
        WRITE (LUPRI,*)
      END IF

* G^A matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G^A matrix transformations:'
        DO ITRAN = 1, NAGTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IAGTRAN(I,ITRAN),I=1,3),(IAGDOTS(I,ITRAN),I=1,MXVAG)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^0 matrix transformations:'
        DO ITRAN = 1, N0FTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (I0FTRAN(I,ITRAN),I=1,2),(I0FDOTS(I,ITRAN),I=1,MXV0F)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^A matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^A matrix transformations:'
        DO ITRAN = 1, NAFTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IAFTRAN(I,ITRAN),I=1,2),(IAFDOTS(I,ITRAN),I=1,MXVAF)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^0{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F{O} matrix transformations:'
        DO ITRAN = 1, N0FATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (I0FATRAN(I,ITRAN),I=1,5),(I0FADOTS(I,ITRAN),I=1,MXV0FA)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^A{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F{O} matrix transformations:'
        DO ITRAN = 1, NAFATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (IAFATRAN(I,ITRAN),I=1,5),(IAFADOTS(I,ITRAN),I=1,MXVAFA)
        END DO
        WRITE (LUPRI,*)
      END IF

* ETA{O} vector calculations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of ETA{O} vector calculations:'
        DO ITRAN = 1, NAEATRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IAEATRAN(I,ITRAN),I=1,2),(IAEADOTS(I,ITRAN),I=1,MXVAEA)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

* chi vector dot products
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of chi vector dot products:'
        DO ITRAN = 1, NXTRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &      IXTRAN(ITRAN),(IXDOTS(I,ITRAN),I=1,MXX)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

* L2 x O2 vector dot products
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of L2 x O2 vector dot products:'
        DO ITRAN = 1, NLTRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &      ILTRAN(ITRAN),(ILDOTS(I,ITRAN),I=1,MXL)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

* xksi vector dot products
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of xksi vector dot products:'
        DO ITRAN = 1, NOTRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &      IOTRAN(ITRAN),(IODOTS(I,ITRAN),I=1,MXO)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCCR_SETUP                           *
*---------------------------------------------------------------------*
c /* deck CCCR_DISP_SETUP */
*=====================================================================*
      SUBROUTINE CCCR_DISP_SETUP(
     &                      MXTRAN2, MXVEC2, MXTRAN3, MXVEC1,
     &                      I0KTRAN, I0KDOTS, W0K,  N0KTRAN,
     &                      I0GTRAN, I0GDOTS, W0G,  N0GTRAN,
     &                      IAGTRAN, IAGDOTS, WAG,  NAGTRAN,
     &                      I0FTRAN, I0FDOTS, W0F,  N0FTRAN,
     &                      IAFTRAN, IAFDOTS, WAF,  NAFTRAN,
     &                      I0FATRAN,I0FADOTS,W0FA, N0FATRAN,
     &                      IAFATRAN,IAFADOTS,WAFA, NAFATRAN,
     &                      IAEATRAN,IAEADOTS,WAEA, NAEATRAN,
     &                      IXTRAN,  IXDOTS,  WX,   NXTRAN,   
     &                      IOTRAN,  IODOTS,  WO,   NOTRAN,
     &                      ILTRAN,  ILDOTS,  WL,   NLTRAN,
     &                      EXPCOF,  MXEXPCF, LADD              )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CCCR dispersion coefficients
*                - list of K^0 matrix transformations 
*                - list of G^0 matrix transformations 
*                - list of G^A matrix transformations 
*                - list of F^0 matrix transformations 
*                - list of F^A matrix transformations 
*                - list of F^0{O} matrix transformations 
*                - list of F^A{O} matrix transformations 
*                - list of ETA^A{O} vector calculations 
*                - list of chi vector dot products
*                - list of xksi vector dot products
*
*     Written by Christof Haettig, march 1998.
*     based on CCCR_SETUP routine
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cccrinf.h>
#include <ccroper.h>
#include <cccperm.h>

* local parameters:
      CHARACTER*(25) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCCR_DISP_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LADD

      INTEGER MXVEC2, MXTRAN2, MXVEC1, MXTRAN3, MXEXPCF

      INTEGER I0KTRAN(5,MXTRAN3)
      INTEGER I0KDOTS(MXVEC1,MXTRAN3)

      INTEGER I0GTRAN(4,MXTRAN2)
      INTEGER I0GDOTS(MXVEC2,MXTRAN2)

      INTEGER IAGTRAN(4,MXTRAN3)
      INTEGER IAGDOTS(MXVEC1,MXTRAN3)

      INTEGER I0FTRAN(3,MXTRAN2)
      INTEGER I0FDOTS(MXVEC2,MXTRAN2)
  
      INTEGER IAFTRAN(3,MXTRAN2)
      INTEGER IAFDOTS(MXVEC2,MXTRAN2)
  
      INTEGER I0FATRAN(5,MXTRAN2)
      INTEGER I0FADOTS(MXVEC2,MXTRAN2)

      INTEGER IAFATRAN(5,MXTRAN3)
      INTEGER IAFADOTS(MXVEC1,MXTRAN3)

      INTEGER IAEATRAN(3,MXTRAN2)
      INTEGER IAEADOTS(MXVEC2,MXTRAN2)

      INTEGER IXTRAN(MXTRAN2)
      INTEGER IXDOTS(MXVEC2,MXTRAN2)

      INTEGER IOTRAN(MXTRAN2)
      INTEGER IODOTS(MXVEC2,MXTRAN2)

      INTEGER ILTRAN(MXTRAN2)
      INTEGER ILDOTS(MXVEC2,MXTRAN2)

      INTEGER N0KTRAN, N0GTRAN, N0FTRAN, N0FATRAN, NXTRAN, NOTRAN
      INTEGER          NAGTRAN, NAFTRAN, NAFATRAN, NAEATRAN, NLTRAN
      INTEGER NCREXPCF

      INTEGER IVEC, ITRAN, I
      INTEGER ISYML, ISYM1, ISYM2, ISYM3, ISYM4
      INTEGER IDISP, IOPER
      INTEGER P
      INTEGER MXV0H, MXV0G, MXVAG, MXV0F, MXVAF, MXV0FA, MXVAFA, MXVAEA
      INTEGER MXX, MXO, MXL

#if defined (SYS_CRAY)
      REAL ZERO
      REAL EXPCOF(MXEXPCF)
      REAL W0K(MXVEC1,MXTRAN3)
      REAL W0G(MXVEC2,MXTRAN2)
      REAL WAG(MXVEC1,MXTRAN3)
      REAL W0F(MXVEC2,MXTRAN2)
      REAL WAF(MXVEC2,MXTRAN2)
      REAL W0FA(MXVEC2,MXTRAN2)
      REAL WAFA(MXVEC1,MXTRAN3)
      REAL WAEA(MXVEC2,MXTRAN2)
      REAL WX(MXVEC2,MXTRAN2)
      REAL WO(MXVEC2,MXTRAN2)
      REAL WL(MXVEC2,MXTRAN2)
      REAL K0CON, G0CON(6), GACON(4), F0CON(3), FACON(12)
      REAL F0ACON(12), FAACON(12), EAACON(12), SUM
      REAL XCON(6), OCON(3)
#else
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION EXPCOF(MXEXPCF)
      DOUBLE PRECISION W0K(MXVEC1,MXTRAN3)
      DOUBLE PRECISION W0G(MXVEC2,MXTRAN2)
      DOUBLE PRECISION WAG(MXVEC1,MXTRAN3)
      DOUBLE PRECISION W0F(MXVEC2,MXTRAN2)
      DOUBLE PRECISION WAF(MXVEC2,MXTRAN2)
      DOUBLE PRECISION W0FA(MXVEC2,MXTRAN2)
      DOUBLE PRECISION WAFA(MXVEC1,MXTRAN3)
      DOUBLE PRECISION WAEA(MXVEC2,MXTRAN2)
      DOUBLE PRECISION WX(MXVEC2,MXTRAN2)
      DOUBLE PRECISION WO(MXVEC2,MXTRAN2)
      DOUBLE PRECISION WL(MXVEC2,MXTRAN2)
      DOUBLE PRECISION K0CON, G0CON(6), GACON(4), F0CON(3), FACON(12)
      DOUBLE PRECISION F0ACON(12), FAACON(12), EAACON(12), SUM
      DOUBLE PRECISION XCON(6), OCON(3)
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER*8 LBL1, LBL2, LBL3, LBL4
      INTEGER ICO1, ICO2, ICO3, ICO4, ICP1, ICP2
      INTEGER IOP(4), ICO(4), ISY(4), IZT(4), IR1(4)
      INTEGER IR2(6), IO2(6), IX2(6), IL2(6)
      INTEGER ICM1(6), ICM2(6), ICM3(6), ICM4(6)

* external functions:
      INTEGER ILSTSYM
      INTEGER ILRCAMP
      INTEGER ILC1AMP
      INTEGER ICR2AMP
      INTEGER ICL2AMP
      INTEGER IRHSCR2
      INTEGER IETACL2


*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN
        N0KTRAN  = 0
        N0GTRAN  = 0
        NAGTRAN  = 0
        N0FTRAN  = 0
        NAFTRAN  = 0
        N0FATRAN = 0
        NAFATRAN = 0
        NAEATRAN = 0
        NXTRAN   = 0
        NOTRAN   = 0
        NLTRAN   = 0
      END IF

      MXV0H  = 0
      MXV0G  = 0
      MXVAG  = 0
      MXV0F  = 0
      MXVAF  = 0
      MXV0FA = 0
      MXVAFA = 0
      MXVAEA = 0
      MXX    = 0
      MXO    = 0
      MXL    = 0

      NCREXPCF = 0
 
      CALL DZERO(EXPCOF,MXEXPCF)

*---------------------------------------------------------------------*
* start loop over all requested dispersion coefficients:
*---------------------------------------------------------------------*
 
      DO IOPER = 1, NCROPER
        IOP(A) = IACROP(IOPER)
        IOP(B) = IBCROP(IOPER)
        IOP(C) = ICCROP(IOPER)
        IOP(D) = IDCROP(IOPER)

        ISY(A) = ISYOPR(IOP(A))
        ISY(B) = ISYOPR(IOP(B))
        ISY(C) = ISYOPR(IOP(C))
        ISY(D) = ISYOPR(IOP(D))

      IF (MULD2H(ISY(A),ISY(B)).EQ.MULD2H(ISY(C),ISY(D))) THEN

      DO IDISP = 1, NCRDISP

        NCREXPCF = NCREXPCF + 1
        
        ICO(A) = ICCAUA(IDISP)
        ICO(B) = ICCAUB(IDISP)
        ICO(C) = ICCAUC(IDISP)
        ICO(D) = ICCAUD(IDISP)
      
        IZT(A) = ILC1AMP(LBLOPR(IOP(A)),ICCAUA(IDISP),ISYML)
        IZT(B) = ILC1AMP(LBLOPR(IOP(B)),ICCAUB(IDISP),ISYML)
        IZT(C) = ILC1AMP(LBLOPR(IOP(C)),ICCAUC(IDISP),ISYML)
        IZT(D) = ILC1AMP(LBLOPR(IOP(D)),ICCAUD(IDISP),ISYML)

        IR1(A) = ILRCAMP(LBLOPR(IOP(A)),ICCAUA(IDISP),ISYML)
        IR1(B) = ILRCAMP(LBLOPR(IOP(B)),ICCAUB(IDISP),ISYML)
        IR1(C) = ILRCAMP(LBLOPR(IOP(C)),ICCAUC(IDISP),ISYML)
        IR1(D) = ILRCAMP(LBLOPR(IOP(D)),ICCAUD(IDISP),ISYML)

        IF ( NO_2NP1_RULE ) THEN
          IR2(AB) = ICR2AMP(LBLOPR(IOP(A)),ICCAUA(IDISP),ISYM1,
     &                      LBLOPR(IOP(B)),ICCAUB(IDISP),ISYM2)
          IR2(AC) = ICR2AMP(LBLOPR(IOP(A)),ICCAUA(IDISP),ISYM1,
     &                      LBLOPR(IOP(C)),ICCAUC(IDISP),ISYM2)
          IR2(AD) = ICR2AMP(LBLOPR(IOP(A)),ICCAUA(IDISP),ISYM1,
     &                      LBLOPR(IOP(D)),ICCAUD(IDISP),ISYM2)
          IR2(BC) = ICR2AMP(LBLOPR(IOP(B)),ICCAUB(IDISP),ISYM1,
     &                      LBLOPR(IOP(C)),ICCAUC(IDISP),ISYM2)
          IR2(BD) = ICR2AMP(LBLOPR(IOP(B)),ICCAUB(IDISP),ISYM1,
     &                      LBLOPR(IOP(D)),ICCAUD(IDISP),ISYM2)
          IR2(CD) = ICR2AMP(LBLOPR(IOP(C)),ICCAUC(IDISP),ISYM1,
     &                      LBLOPR(IOP(D)),ICCAUD(IDISP),ISYM2)
        ELSE

          DO P = 1, 3
            LBL1 = LBLOPR(IOP(I1(P)))  ! Labels
            LBL2 = LBLOPR(IOP(I2(P)))
            LBL3 = LBLOPR(IOP(I3(P)))
            LBL4 = LBLOPR(IOP(I4(P)))
            ICO1 = ICO(I1(P))          ! Cauchy orders
            ICO2 = ICO(I2(P))
            ICO3 = ICO(I3(P))
            ICO4 = ICO(I4(P))

            IF      ( (ICO1+ICO2) .GT. (ICO3+ICO4) ) THEN
              IX2(IP1(P)) = IETACL2(LBL1,ICO1,ISYM1, LBL2,ICO2,ISYM2)
              IO2(IP1(P)) = IRHSCR2(LBL1,ICO1,ISYM1, LBL2,ICO2,ISYM2)
              IL2(IP2(P)) = ICL2AMP(LBL3,ICO3,ISYM3, LBL4,ICO4,ISYM4)
              IR2(IP2(P)) = ICR2AMP(LBL3,ICO3,ISYM3, LBL4,ICO4,ISYM4)
            ELSE IF ( (ICO1+ICO2) .EQ. (ICO3+ICO4) ) THEN
              IX2(IP1(P)) = IETACL2(LBL1,ICO1,ISYM1, LBL2,ICO2,ISYM2)
              IR2(IP1(P)) = ICR2AMP(LBL1,ICO1,ISYM1, LBL2,ICO2,ISYM2)
              IX2(IP2(P)) = IETACL2(LBL3,ICO3,ISYM3, LBL4,ICO4,ISYM4)
              IR2(IP2(P)) = ICR2AMP(LBL3,ICO3,ISYM3, LBL4,ICO4,ISYM4)
              IF (ICO1.GT.0) THEN
               ICM1(IP1(P))=ICL2AMP(LBL1,ICO1-1,ISYM1, LBL2,ICO2,ISYM2)
              END IF
              IF (ICO2.GT.0) THEN
               ICM2(IP1(P))=ICL2AMP(LBL1,ICO1,ISYM1, LBL2,ICO2-1,ISYM2)
              END IF
              IF (ICO3.GT.0) THEN
               ICM3(IP2(P))=ICL2AMP(LBL3,ICO3-1,ISYM3, LBL4,ICO4,ISYM4)
              END IF
              IF (ICO4.GT.0) THEN
               ICM4(IP2(P))=ICL2AMP(LBL3,ICO3,ISYM3, LBL4,ICO4-1,ISYM4)
              END IF
            ELSE IF ( (ICO1+ICO2) .LT. (ICO3+ICO4) ) THEN
              IL2(IP1(P)) = ICL2AMP(LBL1,ICO1,ISYM1, LBL2,ICO2,ISYM2)
              IR2(IP1(P)) = ICR2AMP(LBL1,ICO1,ISYM1, LBL2,ICO2,ISYM2)
              IX2(IP2(P)) = IETACL2(LBL3,ICO3,ISYM3, LBL4,ICO4,ISYM4)
              IO2(IP2(P)) = IRHSCR2(LBL3,ICO3,ISYM3, LBL4,ICO4,ISYM4)
            END IF

          END DO ! P

        END IF

*---------------------------------------------------------------------*
* set up list of K^0 matrix transformations, 1 permutation
*---------------------------------------------------------------------*
        CALL CC_SETH1111(I0KTRAN,I0KDOTS,MXTRAN3,MXVEC1,
     &                   0,IR1(A),IR1(B),IR1(C),IR1(D),ITRAN,IVEC)
        N0KTRAN  = MAX(N0KTRAN,ITRAN)
        MXV0H    = MAX(MXV0H,IVEC)
        K0CON    = W0K(IVEC,ITRAN)

*---------------------------------------------------------------------*
* set up list of G^0 matrix transformations, 6 permutations
*---------------------------------------------------------------------*
      IF (NO_2NP1_RULE)  THEN
        DO P = 1, 6
          CALL CC_SETG112(I0GTRAN,I0GDOTS,MXTRAN2,MXVEC2,   
     &                0,IR1(I1(P)),IR1(I2(P)),IR2(IP(P)),ITRAN,IVEC)
          N0GTRAN  = MAX(N0GTRAN,ITRAN)
          MXV0G    = MAX(MXV0G,IVEC)
          G0CON(P) = W0G(IVEC,ITRAN)
        END DO
      ELSE
        DO P = 1, 6
          G0CON(P) = ZERO
        END DO
      END IF

*---------------------------------------------------------------------*
* set up list of G^A matrix transformations, 4 permutations
*---------------------------------------------------------------------*
        CALL CCQR_SETG(IAGTRAN,IAGDOTS,MXTRAN3,MXVEC1,  
     &                 IZT(A),IR1(B),IR1(C),IR1(D),ITRAN,IVEC)
        NAGTRAN  = MAX(NAGTRAN,ITRAN)
        MXVAG    = MAX(MXVAG,IVEC)
        GACON(1) = WAG(IVEC,ITRAN)

        CALL CCQR_SETG(IAGTRAN,IAGDOTS,MXTRAN3,MXVEC1,
     &                 IZT(B),IR1(A),IR1(C),IR1(D),ITRAN,IVEC)
        NAGTRAN  = MAX(NAGTRAN,ITRAN)
        MXVAG    = MAX(MXVAG,IVEC)
        GACON(2) = WAG(IVEC,ITRAN)

        CALL CCQR_SETG(IAGTRAN,IAGDOTS,MXTRAN3,MXVEC1,
     &                 IZT(C),IR1(B),IR1(A),IR1(D),ITRAN,IVEC)
        NAGTRAN  = MAX(NAGTRAN,ITRAN)
        MXVAG    = MAX(MXVAG,IVEC)
        GACON(3) = WAG(IVEC,ITRAN)

        CALL CCQR_SETG(IAGTRAN,IAGDOTS,MXTRAN3,MXVEC1,
     &                 IZT(D),IR1(B),IR1(C),IR1(A),ITRAN,IVEC)
        NAGTRAN  = MAX(NAGTRAN,ITRAN)
        MXVAG    = MAX(MXVAG,IVEC)
        GACON(4) = WAG(IVEC,ITRAN)

*---------------------------------------------------------------------*
* set up list of F^0 matrix transformations, 3 permutations
*---------------------------------------------------------------------*
      IF (NO_2NP1_RULE)  THEN
        CALL CCQR_SETF(I0FTRAN,I0FDOTS,MXTRAN2,MXVEC2,
     &                 0,IR2(AB),IR2(CD),ITRAN,IVEC)
        N0FTRAN  = MAX(N0FTRAN,ITRAN)
        MXV0F    = MAX(MXV0F,IVEC)
        F0CON(1) = W0F(IVEC,ITRAN)

        CALL CCQR_SETF(I0FTRAN,I0FDOTS,MXTRAN2,MXVEC2,
     &                 0,IR2(AC),IR2(BD),ITRAN,IVEC)
        N0FTRAN  = MAX(N0FTRAN,ITRAN)
        MXV0F    = MAX(MXV0F,IVEC)
        F0CON(2) = W0F(IVEC,ITRAN)

        CALL CCQR_SETF(I0FTRAN,I0FDOTS,MXTRAN2,MXVEC2,
     &                 0,IR2(AD),IR2(BC),ITRAN,IVEC)
        N0FTRAN  = MAX(N0FTRAN,ITRAN)
        MXV0F    = MAX(MXV0F,IVEC)
        F0CON(3) = W0F(IVEC,ITRAN)
      END IF

*---------------------------------------------------------------------*
* set up list of F^A matrix transformations, 12 permutations
*---------------------------------------------------------------------*
      IF (NO_2NP1_RULE)  THEN
        DO P = 1, 6
          CALL CC_SETF12(IAFTRAN,IAFDOTS,MXTRAN2,MXVEC2, ! 1 x 2 x 3,4
     &                   IZT(I1(P)),IR1(I2(P)),IR2(IP(P)),ITRAN,IVEC)
          NAFTRAN  = MAX(NAFTRAN,ITRAN)
          MXVAF    = MAX(MXVAF,IVEC)
          FACON(P) = WAF(IVEC,ITRAN)

          CALL CC_SETF12(IAFTRAN,IAFDOTS,MXTRAN2,MXVEC2, ! 2 x 1 x 3,4
     &                   IZT(I2(P)),IR1(I1(P)),IR2(IP(P)),ITRAN,IVEC)
          NAFTRAN    = MAX(NAFTRAN,ITRAN)
          MXVAF      = MAX(MXVAF,IVEC)
          FACON(6+P) = WAF(IVEC,ITRAN)
        END DO
      ELSE
        DO P = 1, 6
          FACON(P)   = ZERO
          FACON(6+P) = ZERO
        END DO
      END IF

*---------------------------------------------------------------------*
* set up list of F^0{O} matrix transformations, 12 permutations
*---------------------------------------------------------------------*
      IF (NO_2NP1_RULE)  THEN
        DO P = 1, 6
          IF ( ICO(I1(P)).EQ.0 ) THEN
            CALL CC_SETFA12(I0FATRAN,I0FADOTS,MXTRAN2,MXVEC2, ! 1x2x3,4
     &                   0,IOP(I1(P)),IR1(I2(P)),IR2(IP(P)),ITRAN,IVEC)
            N0FATRAN  = MAX(N0FATRAN,ITRAN)
            MXV0FA    = MAX(MXV0FA,IVEC)
            F0ACON(P) = W0FA(IVEC,ITRAN)
          ELSE
            F0ACON(P) = ZERO
          END IF

          IF ( ICO(I2(P)).EQ.0 ) THEN
            CALL CC_SETFA12(I0FATRAN,I0FADOTS,MXTRAN2,MXVEC2, ! 2x1x3,4
     &                   0,IOP(I2(P)),IR1(I1(P)),IR2(IP(P)),ITRAN,IVEC)
            N0FATRAN    = MAX(N0FATRAN,ITRAN)
            MXV0FA      = MAX(MXV0FA,IVEC)
            F0ACON(6+P) = W0FA(IVEC,ITRAN)
          ELSE
            F0ACON(6+P) = ZERO
          END IF
        END DO
      ELSE
        DO P = 1, 6
          F0ACON(P)   = ZERO
          F0ACON(6+P) = ZERO
        END DO
      END IF

*---------------------------------------------------------------------*
* set up list of F^A{O} matrix transformations, 12 permutations
*---------------------------------------------------------------------*
      DO P = 1, 6
        IF ( ICO(I2(P)).EQ.0 ) THEN
          CALL CCQR_SETFA(IAFATRAN,IAFADOTS,MXTRAN3,MXVEC1, ! 1x2x3x4
     &        IZT(I1(P)),IOP(I2(P)),IR1(I3(P)),IR1(I4(P)),ITRAN,IVEC)
          NAFATRAN  = MAX(NAFATRAN,ITRAN)
          MXVAFA    = MAX(MXVAFA,IVEC)
          FAACON(P) = WAFA(IVEC,ITRAN)
        ELSE
          FAACON(P) = ZERO
        END IF

        IF ( ICO(I1(P)).EQ.0 ) THEN
          CALL CCQR_SETFA(IAFATRAN,IAFADOTS,MXTRAN3,MXVEC1, ! 2x1x3x4
     &        IZT(I2(P)),IOP(I1(P)),IR1(I3(P)),IR1(I4(P)),ITRAN,IVEC)
          NAFATRAN    = MAX(NAFATRAN,ITRAN)
          MXVAFA      = MAX(MXVAFA,IVEC)
          FAACON(6+P) = WAFA(IVEC,ITRAN)
        ELSE
          FAACON(6+P) = ZERO
        END IF
      END DO

*---------------------------------------------------------------------*
* set up list of ETA{O} vector calculations, 12 permutations
*---------------------------------------------------------------------*
      IF (NO_2NP1_RULE)  THEN
        DO P = 1, 6
          IF ( ICO(I2(P)).EQ.0 ) THEN
            CALL CCQR_SETEA(IAEATRAN,IAEADOTS,MXTRAN2,MXVEC2, ! 1x2x3,4
     &                     IZT(I1(P)),IOP(I2(P)),IR2(IP(P)),ITRAN,IVEC)
            NAEATRAN  = MAX(NAEATRAN,ITRAN)
            MXVAEA    = MAX(MXVAEA,IVEC)
            EAACON(P) = WAEA(IVEC,ITRAN)
          ELSE
            EAACON(P) = ZERO
          END IF

          IF ( ICO(I1(P)).EQ.0 ) THEN
            CALL CCQR_SETEA(IAEATRAN,IAEADOTS,MXTRAN2,MXVEC2, ! 2x1x3,4
     &                     IZT(I2(P)),IOP(I1(P)),IR2(IP(P)),ITRAN,IVEC)
            NAEATRAN    = MAX(NAEATRAN,ITRAN)
            MXVAEA      = MAX(MXVAEA,IVEC)
            EAACON(6+P) = WAEA(IVEC,ITRAN)
          ELSE
            EAACON(6+P) = ZERO
          END IF
        END DO
      ELSE
        DO P = 1, 6
          EAACON(P)   = ZERO
          EAACON(6+P) = ZERO
        END DO
      END IF

*---------------------------------------------------------------------*
* if we use the 2n+1/2n+2 rules for the second-order Cauchy vectors,
* we here set up list of CX2 x CR2 and CL2 x CO2 dot products 
* (max. 3 permutations) and the list for the F transf. (max. 3 perm.)
*---------------------------------------------------------------------*
      IF (.NOT. NO_2NP1_RULE ) THEN

        DO P = 1, 3
          ICP1 = ICO(I1(P))+ICO(I2(P))
          ICP2 = ICO(I3(P))+ICO(I4(P))

          IF ( ICP1.GT.ICP2 ) THEN
            CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                     IX2(IP1(P)),IR2(IP2(P)),ITRAN,IVEC)
            NXTRAN  = MAX(NXTRAN,ITRAN)
            MXX     = MAX(MXX,IVEC)
            XCON(P) = WX(IVEC,ITRAN)

            CALL CC_SETDOT(IOTRAN,IODOTS,MXTRAN2,MXVEC2,
     &                     IO2(IP1(P)),IL2(IP2(P)),ITRAN,IVEC)
            NOTRAN  = MAX(NOTRAN,ITRAN)
            MXO     = MAX(MXO,IVEC)
            OCON(P) = WO(IVEC,ITRAN)

            XCON(P+3) = ZERO
            F0CON(P)  = ZERO

          ELSE IF ( ICP1.LT.ICP2 ) THEN
            CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                     IX2(IP2(P)),IR2(IP1(P)),ITRAN,IVEC)
            NXTRAN  = MAX(NXTRAN,ITRAN)
            MXX     = MAX(MXX,IVEC)
            XCON(P) = WX(IVEC,ITRAN)

            CALL CC_SETDOT(IOTRAN,IODOTS,MXTRAN2,MXVEC2,
     &                     IO2(IP2(P)),IL2(IP1(P)),ITRAN,IVEC)
            NOTRAN  = MAX(NOTRAN,ITRAN)
            MXO     = MAX(MXO,IVEC)
            OCON(P) = WO(IVEC,ITRAN)

            XCON(P+3) = ZERO
            F0CON(P)  = ZERO

          ELSE IF ( ICP1.EQ.ICP2 ) THEN
            CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                     IX2(IP1(P)),IR2(IP2(P)),ITRAN,IVEC)
            NXTRAN  = MAX(NXTRAN,ITRAN)
            MXX     = MAX(MXX,IVEC)
            XCON(P) = WX(IVEC,ITRAN)

            CALL CC_SETDOT(IXTRAN,IXDOTS,MXTRAN2,MXVEC2,
     &                     IX2(IP2(P)),IR2(IP1(P)),ITRAN,IVEC)
            NXTRAN    = MAX(NXTRAN,ITRAN)
            MXX       = MAX(MXX,IVEC)
            XCON(P+3) = WX(IVEC,ITRAN)

            OCON(P) = ZERO

            IF (ICO(I3(P)).GT.0) THEN
              CALL CC_SETDOT(ILTRAN,ILDOTS,MXTRAN2,MXVEC2,
     &                       ICM3(IP2(P)),IR2(IP1(P)),ITRAN,IVEC)
              NLTRAN  = MAX(NLTRAN,ITRAN)
              MXL     = MAX(MXL,IVEC)
              OCON(P) = OCON(P) + WL(IVEC,ITRAN)
            END IF

            IF (ICO(I4(P)).GT.0) THEN
              CALL CC_SETDOT(ILTRAN,ILDOTS,MXTRAN2,MXVEC2,
     &                       ICM4(IP2(P)),IR2(IP1(P)),ITRAN,IVEC)
              NLTRAN  = MAX(NLTRAN,ITRAN)
              MXL     = MAX(MXL,IVEC)
              OCON(P) = OCON(P) + WL(IVEC,ITRAN)
            END IF

            IF (ICO(I1(P)).GT.0) THEN
              CALL CC_SETDOT(ILTRAN,ILDOTS,MXTRAN2,MXVEC2,
     &                       ICM1(IP1(P)),IR2(IP2(P)),ITRAN,IVEC)
              NLTRAN  = MAX(NLTRAN,ITRAN)
              MXL     = MAX(MXL,IVEC)
              OCON(P) = OCON(P) + WL(IVEC,ITRAN)
            END IF

            IF (ICO(I2(P)).GT.0) THEN
              CALL CC_SETDOT(ILTRAN,ILDOTS,MXTRAN2,MXVEC2,
     &                       ICM2(IP1(P)),IR2(IP2(P)),ITRAN,IVEC)
              NLTRAN  = MAX(NLTRAN,ITRAN)
              MXL     = MAX(MXL,IVEC)
              OCON(P) = OCON(P) + WL(IVEC,ITRAN)
            END IF

            CALL CCQR_SETF(I0FTRAN,I0FDOTS,MXTRAN2,MXVEC2,
     &                     0,IR2(IP1(P)),IR2(IP2(P)),ITRAN,IVEC)
            N0FTRAN  = MAX(N0FTRAN,ITRAN)
            MXV0F    = MAX(MXV0F,IVEC)
            F0CON(P) = W0F(IVEC,ITRAN)
          END IF
         
        END DO

      ELSE
        DO P = 1, 3
          XCON(P)   = ZERO
          XCON(P+3) = ZERO
          OCON(P)   = ZERO
        END DO
      END IF

*---------------------------------------------------------------------*
* add contributions to hyperpolarizabilities:
*---------------------------------------------------------------------*
      IF (LADD) THEN

        EXPCOF(NCRDISP*(IOPER-1)+ IDISP) = 
     &         K0CON     +
     &         G0CON(1)  + G0CON(2)  + G0CON(3)   + G0CON(4)   +
     &         G0CON(5)  + G0CON(6)  +
     &         GACON(1)  + GACON(2)  + GACON(3)   + GACON(4)   +
     &         F0CON(1)  + F0CON(2)  + F0CON(3)   +
     &         FACON(1)  + FACON(2)  + FACON(3)   + FACON(4)   +
     &         FACON(5)  + FACON(6)  + FACON(7)   + FACON(8)   +
     &         FACON(9)  + FACON(10) + FACON(11)  + FACON(12)  + 
     &         F0ACON(1) + F0ACON(2) + F0ACON(3)  + F0ACON(4)  +
     &         F0ACON(5) + F0ACON(6) + F0ACON(7)  + F0ACON(8)  +
     &         F0ACON(9) + F0ACON(10)+ F0ACON(11) + F0ACON(12) +
     &         FAACON(1) + FAACON(2) + FAACON(3)  + FAACON(4)  +
     &         FAACON(5) + FAACON(6) + FAACON(7)  + FAACON(8)  +
     &         FAACON(9) + FAACON(10)+ FAACON(11) + FAACON(12) +
     &         EAACON(1) + EAACON(2) + EAACON(3)  + EAACON(4)  +
     &         EAACON(5) + EAACON(6) + EAACON(7)  + EAACON(8)  +
     &         EAACON(9) + EAACON(10)+ EAACON(11) + EAACON(12) +
     &         XCON(1)   + XCON(2)   + XCON(3)    +
     &         XCON(4)   + XCON(5)   + XCON(6)    +
     &         OCON(1)   + OCON(2)   + OCON(3)    

        IF (LOCDBG) THEN
          WRITE(LUPRI,'(A,3I5)') 'IOPER, IDISP:',IOPER,IDISP
          WRITE(LUPRI,'(A,4I5)') 'IOP:',(IOP(I),I=1,4)
          WRITE(LUPRI,'(A,4I5)') 'ICO:',(ICO(I),I=1,4)
          WRITE(LUPRI,'(A,4I5)') 'ISY:',(ISY(I),I=1,4)
          WRITE(LUPRI,'(A,4I5)') 'IZT:',(IZT(I),I=1,4)
          WRITE(LUPRI,'(A,4I5)') 'IR1:',(IR1(I),I=1,4)
          WRITE(LUPRI,'(A,6I5)') 'IR2:',(IR2(I),I=1,6)
          WRITE(LUPRI,*) 'INDEX:', NCRDISP*(IOPER-1)+ IDISP
          WRITE(LUPRI,*) 'EXPCOF:',EXPCOF(NCRDISP*(IOPER-1)+ IDISP)
          WRITE(LUPRI,*) 'K0CON: ',K0CON
          SUM = K0CON
          WRITE(LUPRI,*) 'G0CON: ',(G0CON(I),I=1,6)
          DO I = 1, 6
            SUM = SUM + G0CON(I)
          END DO
          WRITE(LUPRI,*) 'SUM:',SUM
          WRITE(LUPRI,*) 'GACON: ',(GACON(I),I=1,4)
          DO I = 1, 4
            SUM = SUM + GACON(I)
          END DO
          WRITE(LUPRI,*) 'SUM:',SUM
          WRITE(LUPRI,*) 'F0CON: ',(F0CON(I),I=1,3)
          DO I = 1, 3
            SUM = SUM + F0CON(I)
          END DO
          WRITE(LUPRI,*) 'SUM:',SUM
          WRITE(LUPRI,*) 'FACON: ',(FACON(I),I=1,12)
          DO I = 1, 12
            SUM = SUM + FACON(I)
          END DO
          WRITE(LUPRI,*) 'SUM:',SUM
          WRITE(LUPRI,*) 'F0ACON:',(F0ACON(I),I=1,12)
          DO I = 1, 12
            SUM = SUM + F0ACON(I)
          END DO
          WRITE(LUPRI,*) 'SUM:',SUM
          WRITE(LUPRI,*) 'FAACON:',(FAACON(I),I=1,12)
          DO I = 1, 12
            SUM = SUM + FAACON(I)
          END DO
          WRITE(LUPRI,*) 'SUM:',SUM
          WRITE(LUPRI,*) 'EAACON:',(EAACON(I),I=1,12)
          DO I = 1, 12
            SUM = SUM +EAACON(I)
          END DO
          WRITE(LUPRI,*) 'SUM:',SUM
          WRITE(LUPRI,*) 'XCON:',(XCON(I),I=1,6)
          DO I = 1, 6
            SUM = SUM +XCON(I)
          END DO
          WRITE(LUPRI,*) 'SUM:',SUM
          WRITE(LUPRI,*) 'OCON:',(OCON(I),I=1,3)
          DO I = 1, 3
            SUM = SUM +OCON(I)
          END DO
          WRITE(LUPRI,*) 'SUM:',SUM

          SUM = EXPCOF(NCRDISP*(IOPER-1)+ IDISP) - SUM
          WRITE(LUPRI,*) 'DIFFERENCE:', SUM
        END IF

      END IF

*---------------------------------------------------------------------*
* end loop over all requested dispersion coefficients
*---------------------------------------------------------------------*
      END DO
      END IF
      END DO

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN

* general statistics:
      WRITE(LUPRI,'(////,/3X,A,I5,A)') 'For the requested',NCREXPCF,
     &      ' expansion coefficients for cubic response functions '
      WRITE(LUPRI,'((8X,A,I3,A))') 
     &   ' - ',N0KTRAN,  ' H^0 matrix transformations ',
     &   ' - ',N0GTRAN,  ' G^0 matrix transformations ',
     &   ' - ',NAGTRAN,  ' G^A matrix transformations ',
     &   ' - ',N0FTRAN,  ' F^0 matrix transformations ',
     &   ' - ',NAFTRAN,  ' F^A matrix transformations ',
     &   ' - ',N0FATRAN, ' F^0{O} matrix transformations ',
     &   ' - ',NAFATRAN, ' F^A{O} matrix transformations ',
     &   ' - ',NAEATRAN, ' ETA^A{O} vector calculations ',
     &   ' - ',NXTRAN,   ' CX2 x CR2 dot product calculations ',
     &   ' - ',NOTRAN,   ' CL2 x CO2 dot product calculations '
      WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'


* K^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of K^0 matrix transformations:'
        DO ITRAN = 1, N0KTRAN
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') MSGDBG,
     &     (I0KTRAN(I,ITRAN),I=1,4),(I0KDOTS(I,ITRAN),I=1,MXV0H)
        END DO
        WRITE (LUPRI,*)
      END IF

* G^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G^0 matrix transformations:'
        DO ITRAN = 1, N0GTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (I0GTRAN(I,ITRAN),I=1,3),(I0GDOTS(I,ITRAN),I=1,MXV0G)
        END DO
        WRITE (LUPRI,*)
      END IF

* G^A matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G^A matrix transformations:'
        DO ITRAN = 1, NAGTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IAGTRAN(I,ITRAN),I=1,3),(IAGDOTS(I,ITRAN),I=1,MXVAG)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^0 matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^0 matrix transformations:'
        DO ITRAN = 1, N0FTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (I0FTRAN(I,ITRAN),I=1,2),(I0FDOTS(I,ITRAN),I=1,MXV0F)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^A matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F^A matrix transformations:'
        DO ITRAN = 1, NAFTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IAFTRAN(I,ITRAN),I=1,2),(IAFDOTS(I,ITRAN),I=1,MXVAF)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^0{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F{O} matrix transformations:'
        DO ITRAN = 1, N0FATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (I0FATRAN(I,ITRAN),I=1,5),(I0FADOTS(I,ITRAN),I=1,MXV0FA)
        END DO
        WRITE (LUPRI,*)
      END IF

* F^A{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F{O} matrix transformations:'
        DO ITRAN = 1, NAFATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (IAFATRAN(I,ITRAN),I=1,5),(IAFADOTS(I,ITRAN),I=1,MXVAFA)
        END DO
        WRITE (LUPRI,*)
      END IF

* ETA{O} vector calculations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of ETA{O} vector calculations:'
        DO ITRAN = 1, NAEATRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IAEATRAN(I,ITRAN),I=1,2),(IAEADOTS(I,ITRAN),I=1,MXVAEA)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

* eta vector dot products:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of eta vector dot products:'
        DO ITRAN = 1, NXTRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &      IXTRAN(ITRAN),(IXDOTS(I,ITRAN),I=1,MXX)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

* xksi vector dot products:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of xksi vector dot products:'
        DO ITRAN = 1, NOTRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &      IOTRAN(ITRAN),(IODOTS(I,ITRAN),I=1,MXO)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

* CL2 x CR2 vector dot products:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of CL2 x CR2 vector dot products:'
        DO ITRAN = 1, NLTRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &      ILTRAN(ITRAN),(ILDOTS(I,ITRAN),I=1,MXL)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCCR_DISP_SETUP                      *
*---------------------------------------------------------------------*
c /* deck CC_SETH1111 */
*=====================================================================*
      SUBROUTINE CC_SETH1111(IHTRAN,IHDOTS,MXTRAN,MXVEC,IZETAV,
     &                       ITAMPA,ITAMPB,ITAMPC,ITAMPD,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintain a list of dot products of H matrix 
*             transformations with right amplitude vectors:  
*                       (Z*H*T*T*T) * T
*    N.B.: assumes that all four T vectors belong to the same list
*
*             IHTRAN - list of H matrix transformations
*             IHDOTS - list of vectors it should be dottet on
*             
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimesion for IHDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*             ITAMPC - index of amplitude vector C
*             ITAMPD - index of amplitude vector D
*
*             ITRAN - index in IHTRAN list
*             IVEC  - second index in IHDOTS list
*
*    Written by Christof Haettig, february 1997.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
      INTEGER MXVEC, MXTRAN
      INTEGER IHTRAN(5,MXTRAN)
      INTEGER IHDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDA, LFNDB, LFNDC, LFNDD
      INTEGER IZETAV, ITAMPA, ITAMPB, ITAMPC, ITAMPD
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LHTST, LHEND
      INTEGER IL, IA, IB,IC
      LHTST(ITRAN,IL,IA,IB,IC) = IHTRAN(1,ITRAN).EQ.IL .AND. (
     &     (IHTRAN(2,ITRAN).EQ.IA .AND. IHTRAN(3,ITRAN).EQ.IB
     &                            .AND. IHTRAN(4,ITRAN).EQ.IC) .OR. 
     &     (IHTRAN(2,ITRAN).EQ.IB .AND. IHTRAN(3,ITRAN).EQ.IA
     &                            .AND. IHTRAN(4,ITRAN).EQ.IC) .OR.
     &     (IHTRAN(2,ITRAN).EQ.IC .AND. IHTRAN(3,ITRAN).EQ.IA
     &                            .AND. IHTRAN(4,ITRAN).EQ.IB) .OR. 
     &     (IHTRAN(2,ITRAN).EQ.IA .AND. IHTRAN(3,ITRAN).EQ.IC
     &                            .AND. IHTRAN(4,ITRAN).EQ.IB) .OR.
     &     (IHTRAN(2,ITRAN).EQ.IB .AND. IHTRAN(3,ITRAN).EQ.IC
     &                            .AND. IHTRAN(4,ITRAN).EQ.IA) .OR. 
     &     (IHTRAN(2,ITRAN).EQ.IC .AND. IHTRAN(3,ITRAN).EQ.IB
     &                            .AND. IHTRAN(4,ITRAN).EQ.IA)      )
      LHEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &               (IHTRAN(1,ITRAN)+IHTRAN(2,ITRAN)+
     &                IHTRAN(3,ITRAN)+IHTRAN(4,ITRAN) ).LE.0

*---------------------------------------------------------------------*
* set up list of H matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDA = LHTST(ITRAN,IZETAV,ITAMPB,ITAMPC,ITAMPD)
      LFNDB = LHTST(ITRAN,IZETAV,ITAMPC,ITAMPD,ITAMPA)
      LFNDC = LHTST(ITRAN,IZETAV,ITAMPD,ITAMPA,ITAMPB)
      LFNDD = LHTST(ITRAN,IZETAV,ITAMPA,ITAMPB,ITAMPC)

      DO WHILE ( .NOT. (LFNDA.OR.LFNDB.OR.LFNDC.OR.LFNDD
     &                       .OR.LHEND(ITRAN)) )
        ITRAN = ITRAN + 1
        LFNDA = LHTST(ITRAN,IZETAV,ITAMPB,ITAMPC,ITAMPD)
        LFNDB = LHTST(ITRAN,IZETAV,ITAMPC,ITAMPD,ITAMPA)
        LFNDC = LHTST(ITRAN,IZETAV,ITAMPD,ITAMPA,ITAMPB)
        LFNDD = LHTST(ITRAN,IZETAV,ITAMPA,ITAMPB,ITAMPC)
      END DO

      IF (.NOT.(LFNDA.OR.LFNDB.OR.LFNDC.OR.LFNDD)) THEN
        IHTRAN(1,ITRAN) = IZETAV
        IHTRAN(2,ITRAN) = ITAMPA
        IHTRAN(3,ITRAN) = ITAMPB
        IHTRAN(4,ITRAN) = ITAMPC
        ITAMP = ITAMPD
      ELSE 
        IF (LFNDA) ITAMP = ITAMPA
        IF (LFNDB) ITAMP = ITAMPB
        IF (LFNDC) ITAMP = ITAMPC
        IF (LFNDD) ITAMP = ITAMPD
      END IF

      IVEC = 1
      DO WHILE (IHDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &           IHDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IHDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'ITAMPA,ITAMPB,ITAMPC,ITAMPD:',
     &              ITAMPA,ITAMPB,ITAMPC,ITAMPD
        IDX = 1
        DO WHILE (.NOT. LHEND(IDX))
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') 'CC_SETH1111>',
     &       (IHTRAN(I,IDX),I=1,4),(IHDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETH1111.')
      END IF
      
      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETH1111                          *
*---------------------------------------------------------------------*
c /* deck CC_SETG112 */
*=====================================================================*
      SUBROUTINE CC_SETG112(IGTRAN,IGDOTS,MXTRAN,MXVEC,
     &                      IZETAV,ITAMPA,ITAMPB,ITAMPC,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintain a list of dot products of G matrix 
*             transformations with right amplitude vectors:  
*                       (Z*G*T^A*T^B) * T^C
*             assumes that T^A and T^B belong to one list,
*             and T^C belongs to a second list
*
*             IGTRAN - list of G matrix transformations
*             IGDOTS - list of vectors it should be dottet on
*             
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimesion for IGDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*             ITAMPC - index of amplitude vector C
*
*             ITRAN - index in IGTRAN list
*             IVEC  - second index in IGDOTS list
*
*    Written by Christof Haettig, februar 1997.
*
*=====================================================================*
      IMPLICIT NONE  

#include <priunit.h>
      INTEGER MXVEC, MXTRAN
      INTEGER IGTRAN(4,MXTRAN)
      INTEGER IGDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDC
      INTEGER IZETAV, ITAMPA, ITAMPB, ITAMPC
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LGTST, LGEND
      INTEGER IL, IA, IB
      LGTST(ITRAN,IL,IA,IB) = IGTRAN(1,ITRAN).EQ.IL .AND.
     &   ( (IGTRAN(2,ITRAN).EQ.IA .AND. IGTRAN(3,ITRAN).EQ.IB) .OR. 
     &     (IGTRAN(2,ITRAN).EQ.IB .AND. IGTRAN(3,ITRAN).EQ.IA)       )
      LGEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &   (IGTRAN(1,ITRAN)+IGTRAN(2,ITRAN)+IGTRAN(3,ITRAN)).LE.0

*---------------------------------------------------------------------*
* set up list of G matrix transformations
*---------------------------------------------------------------------*
        ITRAN = 1
        LFNDC = LGTST(ITRAN,IZETAV,ITAMPA,ITAMPB)

        DO WHILE ( .NOT. (LFNDC.OR.LGEND(ITRAN)) )
         ITRAN = ITRAN + 1
         LFNDC = LGTST(ITRAN,IZETAV,ITAMPA,ITAMPB)
        END DO

        IF (.NOT.LFNDC) THEN
          IGTRAN(1,ITRAN) = IZETAV
          IGTRAN(2,ITRAN) = ITAMPA
          IGTRAN(3,ITRAN) = ITAMPB
          ITAMP = ITAMPC
        ELSE 
          IF (LFNDC) ITAMP = ITAMPC
        END IF

        IVEC = 1
        DO WHILE (IGDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &             IGDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
          IVEC = IVEC + 1
        END DO

        IGDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
        IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
          WRITE (LUPRI,*) 'IVEC :',IVEC
          WRITE (LUPRI,*) 'ITRAN:',ITRAN
          WRITE (LUPRI,*) 'ITAMPA,ITAMPB,ITAMPC:',ITAMPA,ITAMPB,ITAMPC
          IDX = 1
          DO WHILE (.NOT. LGEND(IDX))
            WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') 'CC_SETG112>',
     &         (IGTRAN(I,IDX),I=1,3),(IGDOTS(I,IDX),I=1,MXVEC)
            IDX = IDX + 1
          END DO
          CALL FLSHFO(LUPRI)
          CALL QUIT('Overflow error in CC_SETG112.')
        END IF
      
        RETURN
        END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETG112                           *
*---------------------------------------------------------------------*
c /* deck CC_SETF12 */
*=====================================================================*
       SUBROUTINE CC_SETF12(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                      IZETAV,ITAMPA,ITAMPB,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintain a list of dot products of F matrix 
*             transformations with right amplitude vectors
*                       (Z*F*T^A) * T^B
*             assumes that T^A and T^B belong to different lists
*
*             IFTRAN - list of F matrix transformations
*             IFDOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IFDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*
*             ITRAN - index in IFTRAN list
*             IVEC  - second index in IFDOTS list
*
*    Written by Christof Haettig, februar 1997.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IFTRAN(3,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDB
      INTEGER IZETAV, ITAMPA, ITAMPB
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LFTST, LFEND
      INTEGER IL, IA
      LFTST(ITRAN,IL,IA) = 
     &      IFTRAN(1,ITRAN).EQ.IL .AND. IFTRAN(2,ITRAN).EQ.IA 
      LFEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &   (IFTRAN(1,ITRAN)+IFTRAN(2,ITRAN)).LE.0

*---------------------------------------------------------------------*
* set up list of F matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDB = LFTST(ITRAN,IZETAV,ITAMPA)

      DO WHILE ( .NOT. (LFNDB.OR.LFEND(ITRAN)) )
       ITRAN = ITRAN + 1
       LFNDB = LFTST(ITRAN,IZETAV,ITAMPA)
      END DO

      IF (.NOT.LFNDB) THEN
        IFTRAN(1,ITRAN) = IZETAV
        IFTRAN(2,ITRAN) = ITAMPA
        ITAMP = ITAMPB
      ELSE 
        IF (LFNDB) ITAMP = ITAMPB
      END IF

      IVEC = 1
      DO WHILE (IFDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &           IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IFDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'ITAMPA,ITAMPB:',ITAMPA,ITAMPB
        IDX = 1
        DO WHILE ( .NOT. LFEND(IDX) )
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') 'CC_SETF12>',
     &       (IFTRAN(I,IDX),I=1,2),(IFDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETF12.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETF12                            *
*---------------------------------------------------------------------*
c /* deck CC_SETAA */
*=====================================================================*
      SUBROUTINE CC_SETAA(IAATRAN,IAADOTS,MXTRAN,MXVEC,
     &                    IZETAV,IOPER,ITAMPB,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintain a list of dot products of A{O} matrix 
*             transformations with right amplitude vectors
*                       Z * A{O} * T^B
*
*             IAATRAN - list of A{O} matrix transformations
*             IAADOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IFDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             IOPER  - index of operator O
*             ITAMPB - index of amplitude vector B
*
*             ITRAN - index in IAATRAN list
*             IVEC  - second index in IAADOTS list
*
*    Written by Christof Haettig, Mai 2003.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <cclists.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IAATRAN(MXDIM_AATRAN,MXTRAN)
      INTEGER IAADOTS(MXVEC,MXTRAN)

      LOGICAL LFNDB
      INTEGER IZETAV, IOPER, ITAMPB
      INTEGER ITRAN, IVEC, I, IDX

* statement  functions:
      LOGICAL LAATST, LFEND
      INTEGER IO, IB
      LAATST(ITRAN,IO,IB) = 
     &      IAATRAN(1,ITRAN).EQ.IO .AND. IAATRAN(2,ITRAN).EQ.IB
      LFEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &   (IAATRAN(1,ITRAN)+IAATRAN(2,ITRAN)).LE.0

*---------------------------------------------------------------------*
* set up list of A{O} matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDB = LAATST(ITRAN,IOPER,ITAMPB)

      DO WHILE ( .NOT. (LFNDB.OR.LFEND(ITRAN)) )
       ITRAN = ITRAN + 1
       LFNDB = LAATST(ITRAN,IOPER,ITAMPB)
      END DO

      IF (.NOT.LFNDB) THEN
        IAATRAN(1,ITRAN) = IOPER
        IAATRAN(2,ITRAN) = ITAMPB
      END IF

      IVEC = 1
      DO WHILE (IAADOTS(IVEC,ITRAN).NE.IZETAV .AND.
     &          IAADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IAADOTS(IVEC,ITRAN) = IZETAV

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'IOPER,ITAMPB:',IOPER,ITAMPB
        IDX = 1
        DO WHILE ( .NOT. LFEND(IDX) )
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') 'CC_SETAA>',
     &       (IAATRAN(I,IDX),I=1,2),(IAADOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETAA.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETAA                             *
*---------------------------------------------------------------------*
c /* deck CC_SETFA12 */
*=====================================================================*
       SUBROUTINE CC_SETFA12(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                       IZETAV,IOPER,ITAMPA,ITAMPB,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintains a list of dot products of F{O} matrix 
*             transformations with right amplitude vectors:
*                        (Z*F{O}*T^A) * T^B
*             assumes that T^A and T^B belong to different lists
*
*             IFTRAN - list of F matrix transformations
*             IFDOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IFDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             IOPER  - index of property operator 
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*
*             ITRAN - index in IFTRAN list
*             IVEC  - second index in IFDOTS list
*
*    Written by Christof Haettig, februar 1997.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IFTRAN(5,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDB
      INTEGER IZETAV, IOPER, ITAMPA, ITAMPB
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LFATST, LFAEND
      INTEGER IL, IA, IO
      LFATST(ITRAN,IL,IO,IA) = IFTRAN(1,ITRAN).EQ.IL 
     &       .AND. IFTRAN(2,ITRAN).EQ.IO .AND. IFTRAN(3,ITRAN).EQ.IA 
      LFAEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &      (IFTRAN(1,ITRAN)+IFTRAN(2,ITRAN)+IFTRAN(3,ITRAN)).LE.0 


*---------------------------------------------------------------------*
* set up list of F{A} matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDB  = LFATST(ITRAN,IZETAV,IOPER,ITAMPA)

      DO WHILE ( .NOT. (LFNDB.OR.LFAEND(ITRAN)))
       ITRAN = ITRAN + 1
       LFNDB  = LFATST(ITRAN,IZETAV,IOPER,ITAMPA)
      END DO

      IF (.NOT.LFNDB) THEN
        IFTRAN(1,ITRAN) = IZETAV
        IFTRAN(2,ITRAN) = IOPER
        IFTRAN(3,ITRAN) = ITAMPA
        IFTRAN(4,ITRAN) = 0
        IFTRAN(5,ITRAN) = 0
        ITAMP = ITAMPB
      ELSE 
        IF (LFNDB) ITAMP = ITAMPB
      END IF

      IVEC = 1
      DO WHILE (IFDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &            IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IFDOTS(IVEC,ITRAN) = ITAMP

C     WRITE (LUPRI,*) 'CC_SETFA12>',IZETAV,IOPER,ITAMPA,ITAMPB,ITRAN,IVEC
*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'ITAMPA,ITAMPB:',ITAMPA,ITAMPB
        IDX = 1
        DO WHILE ( .NOT. LFAEND(IDX) )
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') 'CC_SETFA12>',
     &       (IFTRAN(I,IDX),I=1,5),(IFDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETFA12.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETFA12                           *
*---------------------------------------------------------------------*
c /* deck CC_SETDOT */
*=====================================================================*
       SUBROUTINE CC_SETDOT(IDTRAN,IDDOTS,MXTRAN,MXVEC,
     &                      ICHIA,ITAMPB,ITRAN,IVEC   )
*---------------------------------------------------------------------*
*
*    Purpose: maintain a list of dot products of 'ICHIA' vectors
*             times 'ITAMPB' vectors
*                       X^A * T^B
*             assumes that X^A and T^B belong to different lists
*
*             IDTRAN - list of ICHIA vectors
*             IDDOTS - list of vectors they should be dotted on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IDDOTS
*      
*             ICHIA  - index of ICHIA vector
*             ITAMPB - index of ITAMPB vector
*
*             ITRAN - index in IDTRAN list
*             IVEC  - second index in IDDOTS list
*
*    Written by Christof Haettig, februar 1997.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IDTRAN(MXTRAN)
      INTEGER IDDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDB
      INTEGER ICHIA, ITAMPB
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LFTST, LFEND
      INTEGER IA
      LFTST(ITRAN,IA) = IDTRAN(ITRAN).EQ.IA 
      LFEND(ITRAN) = ITRAN.GT.MXTRAN .OR. IDTRAN(ITRAN).LE.0

*---------------------------------------------------------------------*
* set up list of ICHIA vectors
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDB = LFTST(ITRAN,ICHIA)

      DO WHILE ( .NOT. (LFNDB.OR.LFEND(ITRAN)) )
       ITRAN = ITRAN + 1
       LFNDB = LFTST(ITRAN,ICHIA)
      END DO

      IF (.NOT.LFNDB) THEN
        IDTRAN(ITRAN) = ICHIA
        ITAMP = ITAMPB
      ELSE 
        IF (LFNDB) ITAMP = ITAMPB
      END IF

      IVEC = 1
      DO WHILE (IDDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &           IDDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IDDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC, MXVEC :',IVEC,MXVEC
        WRITE (LUPRI,*) 'ITRAN,MXTRAN:',ITRAN,MXTRAN
        WRITE (LUPRI,*) 'ICHIA,ITAMPB:',ICHIA,ITAMPB
        IDX = 1
        DO WHILE ( .NOT. LFEND(IDX) )
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') 'CC_SETDOT>',
     &       IDTRAN(IDX),(IDDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETDOT.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETDOT                            *
*---------------------------------------------------------------------*
