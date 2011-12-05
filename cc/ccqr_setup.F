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
c /* deck ccqr_setup */
*=====================================================================*
      SUBROUTINE CCQR_SETUP(MXTRAN, MXVEC, 
     &                      IGTRAN, IGDOTS, NGTRAN,
     &                      IFTRAN, IFDOTS, NFTRAN,
     &                      IBTRAN, IBDOTS, NBTRAN,
     &                      IKTRAN, IKDOTS, NKTRAN,
     &                      IFATRAN,IFADOTS,NFATRAN,
     &                      IAATRAN,IAADOTS,NAATRAN,
     &                      IEATRAN,IEADOTS,NEATRAN,
     &                      IR2TRAN,IR2DOTS,NR2TRAN,
     &                      IXETRAN,IXDOTS, IEDOTS, NXETRAN,
     &                      WORK,   LWORK  )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CCQR section
*                - list of G matrix transformations 
*                - list of F/B matrix transformations 
*                - list of F{O} matrix transformations 
*                - list of ETA{O} vector calculations 
*                - list of X1 x R2 dot products
*                - list of Xksi{O2} x L1 and Eta{O2} x R1 dot products
*    CCMM JK+OC  - list of K matrix transformations
*
*     Written by Christof Haettig, november 1996.
*     (option for B matrix transformation introduced in 12/97)
*     (USE_R2 option introduced 3/98)                              
*     (Orbital relaxation for one operator introduced in 6/99)
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccqrinf.h>
#include <ccroper.h>
#include <cclists.h>
#include <ccsections.h>
#include <ccslvinf.h>

* local parameters:
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCQR_SETUP> ')
      LOGICAL LOCDBG, OLD
      PARAMETER (LOCDBG = .FALSE., OLD = .FALSE.)

      INTEGER MXVEC, MXTRAN, LWORK

      INTEGER IGTRAN(MXDIM_GTRAN,MXTRAN)
      INTEGER IGDOTS(MXVEC,MXTRAN)

      INTEGER IBTRAN(MXDIM_BTRAN,MXTRAN)
      INTEGER IBDOTS(MXVEC,MXTRAN)

      INTEGER IKTRAN(MXDIM_KTRAN,MXTRAN)
      INTEGER IKDOTS(MXVEC,MXTRAN)
  
      INTEGER IFTRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)
  
      INTEGER IFATRAN(MXDIM_FATRAN,MXTRAN)
      INTEGER IFADOTS(MXVEC,MXTRAN)

      INTEGER IAATRAN(MXDIM_AATRAN,MXTRAN)
      INTEGER IAADOTS(MXVEC,MXTRAN)

      INTEGER IEATRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IEADOTS(MXVEC,MXTRAN)

      INTEGER IR2TRAN(MXTRAN)
      INTEGER IR2DOTS(MXVEC,MXTRAN)

      INTEGER IXETRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IXDOTS(MXVEC,MXTRAN), IEDOTS(MXVEC,MXTRAN)

      INTEGER NGTRAN,NFTRAN,NR2TRAN,NBTRAN,NFATRAN,NEATRAN,NXETRAN
      INTEGER NKTRAN,NQRRESF,NAATRAN

      INTEGER ISYMA, ISYMB, ISYMC, ISYML, ISYML1, ISYML2
      INTEGER IFREQ, IOPER, ISIGN

      INTEGER IOPERA, IOPERB, IOPERC
      INTEGER IZETVA, IZETVB, IZETVC
      INTEGER ITAMPA, ITAMPB, ITAMPC
      INTEGER IKAPA,  IKAPB,  IKAPC
      INTEGER IR2AB, IR2BC, IR2AC, IETAA, IETAB, IETAC, MXVAA
      INTEGER IVEC, ITRAN, I, ISYSOP, IOPSOP, ISGNSOP, INUM
      INTEGER MXV0G, MXVAF, MXVAEA, MXVXE, MXVAB, MXV0FA, MXVR2
      INTEGER MXVAK

      LOGICAL LORXA, LORXB, LORXC, LPDBSA, LPDBSB, LPDBSC

      CHARACTER*(8) LABSOP

#if defined (SYS_CRAY)
      REAL SIGN, WORK(LWORK)
#else
      DOUBLE PRECISION SIGN, WORK(LWORK)
#endif


* external functions:
      INTEGER IR1TAMP
      INTEGER IR1KAPPA
      INTEGER IL1ZETA
      INTEGER IR2TAMP
      INTEGER IROPER
      INTEGER IETA1


*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      DO ITRAN = 1, MXTRAN
        IR2TRAN(ITRAN) = 0

        DO I = 1, 3
          IFTRAN(I,ITRAN) = 0
          IBTRAN(I,ITRAN) = 0
          IKTRAN(I,ITRAN) = 0
        END DO

        DO I = 1, 4
          IGTRAN(I,ITRAN) = 0
        END DO

        DO I = 1, MXDIM_FATRAN
          IFATRAN(I,ITRAN) = 0
        END DO

        DO I = 1, MXDIM_AATRAN
          IAATRAN(I,ITRAN) = 0
        END DO

        DO I = 1, MXDIM_XEVEC
          IXETRAN(I,ITRAN) = 0
          IEATRAN(I,ITRAN) = 0
        END DO
        IXETRAN(3,ITRAN) = -1
        IXETRAN(4,ITRAN) = -1
        IEATRAN(3,ITRAN) = -1
        IEATRAN(4,ITRAN) = -1

        DO IVEC  = 1, MXVEC
          IGDOTS(IVEC,ITRAN)  = 0
          IFDOTS(IVEC,ITRAN)  = 0
          IBDOTS(IVEC,ITRAN)  = 0
          IKDOTS(IVEC,ITRAN)  = 0
          IFADOTS(IVEC,ITRAN) = 0
          IAADOTS(IVEC,ITRAN) = 0
          IEADOTS(IVEC,ITRAN) = 0
          IR2DOTS(IVEC,ITRAN) = 0
          IXDOTS(IVEC,ITRAN)  = 0
          IEDOTS(IVEC,ITRAN)  = 0
        END DO
      END DO

      NGTRAN  = 0
      NFTRAN  = 0
      NBTRAN  = 0
      NKTRAN  = 0
      NFATRAN = 0
      NAATRAN = 0
      NEATRAN = 0
      NR2TRAN = 0
      NXETRAN = 0
      NQRRESF = 0
 
      MXV0G  = 0
      MXVAF  = 0
      MXVXE  = 0
      MXVAB  = 0
      MXVAK  = 0
      MXV0FA = 0
      MXVR2  = 0
      MXVAEA = 0
      MXVAA  = 0

*---------------------------------------------------------------------*
* start loop over all requested hyperpolarizabilities
*---------------------------------------------------------------------*
 
      DO IOPER = 1, NQROPER
        IOPERA = IAQROP(IOPER)
        IOPERB = IBQROP(IOPER)
        IOPERC = ICQROP(IOPER)

        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        ISYMC  = ISYOPR(IOPERC)

      IF (MULD2H(ISYMB,ISYMC).EQ.ISYMA) THEN

        LORXA  = LAQLRX(IOPER)
        LORXB  = LBQLRX(IOPER)
        LORXC  = LCQLRX(IOPER)

        LPDBSA = LPDBSOP(IOPERA)
        LPDBSB = LPDBSOP(IOPERB)
        LPDBSC = LPDBSOP(IOPERC)

      DO IFREQ = 1, NQRFREQ

        NQRRESF = NQRRESF + 1

      DO ISIGN = +1, -1, -2
        SIGN = DBLE(ISIGN)

        IZETVA = IL1ZETA(LBLOPR(IOPERA),LORXA,SIGN*AQRFR(IFREQ),ISYML)
        IZETVB = IL1ZETA(LBLOPR(IOPERB),LORXB,SIGN*BQRFR(IFREQ),ISYML)
        IZETVC = IL1ZETA(LBLOPR(IOPERC),LORXC,SIGN*CQRFR(IFREQ),ISYML)

        ITAMPA = IR1TAMP(LBLOPR(IOPERA),LORXA,SIGN*AQRFR(IFREQ),ISYML)
        ITAMPB = IR1TAMP(LBLOPR(IOPERB),LORXB,SIGN*BQRFR(IFREQ),ISYML)
        ITAMPC = IR1TAMP(LBLOPR(IOPERC),LORXC,SIGN*CQRFR(IFREQ),ISYML)

        IKAPA = 0
        IKAPB = 0
        IKAPC = 0
        IF(LORXA)IKAPA=IR1KAPPA(LBLOPR(IOPERA),SIGN*AQRFR(IFREQ),ISYML)
        IF(LORXB)IKAPB=IR1KAPPA(LBLOPR(IOPERB),SIGN*BQRFR(IFREQ),ISYML)
        IF(LORXC)IKAPC=IR1KAPPA(LBLOPR(IOPERC),SIGN*CQRFR(IFREQ),ISYML)
 
        IF ( USE_R2 ) THEN
         IR2AB=IR2TAMP(LBLOPR(IOPERA),.FALSE.,SIGN*AQRFR(IFREQ),ISYML1,
     &                 LBLOPR(IOPERB),.FALSE.,SIGN*BQRFR(IFREQ),ISYML2)
         IR2AC=IR2TAMP(LBLOPR(IOPERA),.FALSE.,SIGN*AQRFR(IFREQ),ISYML1,
     &                 LBLOPR(IOPERC),.FALSE.,SIGN*CQRFR(IFREQ),ISYML2)
         IR2BC=IR2TAMP(LBLOPR(IOPERB),.FALSE.,SIGN*BQRFR(IFREQ),ISYML1,
     &                 LBLOPR(IOPERC),.FALSE.,SIGN*CQRFR(IFREQ),ISYML2)

         IETAA = IETA1(LBLOPR(IOPERA),.FALSE.,SIGN*AQRFR(IFREQ),ISYML)
         IETAB = IETA1(LBLOPR(IOPERB),.FALSE.,SIGN*BQRFR(IFREQ),ISYML)
         IETAC = IETA1(LBLOPR(IOPERC),.FALSE.,SIGN*CQRFR(IFREQ),ISYML)
        END IF

*---------------------------------------------------------------------*
* set up list of G matrix transformations, 1 permutation
*---------------------------------------------------------------------*
        CALL CCQR_SETG(IGTRAN,IGDOTS,MXTRAN,MXVEC,
     &                 0,ITAMPA,ITAMPB,ITAMPC,ITRAN,IVEC)
        NGTRAN = MAX(NGTRAN,ITRAN)
        MXV0G  = MAX(MXV0G,IVEC)

*---------------------------------------------------------------------*
* set up list of F/B matrix transformations, 3 permutations
* or if USE_R2 is set a list of XC x CR2 dot products, 3 permutations
* (can be used alternatively from the hyperpolarizability driver)
*---------------------------------------------------------------------*
       IF ( (.NOT.USEBTRAN) .AND. (.NOT.USE_R2) ) THEN
* F matrix transformations:
        CALL CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                 IZETVA,ITAMPB,ITAMPC,ITRAN,IVEC)
        NFTRAN = MAX(NFTRAN,ITRAN)
        MXVAF  = MAX(MXVAF,IVEC)

        CALL CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                 IZETVB,ITAMPC,ITAMPA,ITRAN,IVEC)
        NFTRAN = MAX(NFTRAN,ITRAN)
        MXVAF  = MAX(MXVAF,IVEC)

        CALL CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                 IZETVC,ITAMPA,ITAMPB,ITRAN,IVEC)
        NFTRAN = MAX(NFTRAN,ITRAN)
        MXVAF  = MAX(MXVAF,IVEC)

       ELSE IF ( USEBTRAN .AND. (.NOT.USE_R2) ) THEN
* equivalent B matrix transformations:
        CALL CC_SETB11(IBTRAN,IBDOTS,MXTRAN,MXVEC,
     &                 IZETVA,ITAMPB,ITAMPC,ITRAN,IVEC)
        NBTRAN = MAX(NBTRAN,ITRAN)
        MXVAB  = MAX(MXVAB,IVEC)

        CALL CC_SETB11(IBTRAN,IBDOTS,MXTRAN,MXVEC,
     &                 IZETVB,ITAMPC,ITAMPA,ITRAN,IVEC)
        NBTRAN = MAX(NBTRAN,ITRAN)
        MXVAB  = MAX(MXVAB,IVEC)

        CALL CC_SETB11(IBTRAN,IBDOTS,MXTRAN,MXVEC,
     &                 IZETVC,ITAMPA,ITAMPB,ITRAN,IVEC)
        NBTRAN = MAX(NBTRAN,ITRAN)
        MXVAB  = MAX(MXVAB,IVEC)

       ELSE IF ( USE_R2 ) THEN
* R2 x X1 dot products:
        CALL CC_SETDOT(IR2TRAN,IR2DOTS,MXTRAN,MXVEC,
     &                 IR2AB,IZETVC,ITRAN,IVEC)
        NR2TRAN = MAX(NR2TRAN,ITRAN)
        MXVR2  = MAX(MXVR2,IVEC)

        CALL CC_SETDOT(IR2TRAN,IR2DOTS,MXTRAN,MXVEC,
     &                 IR2AC,IZETVB,ITRAN,IVEC)
        NR2TRAN = MAX(NR2TRAN,ITRAN)
        MXVR2  = MAX(MXVR2,IVEC)
     
        CALL CC_SETDOT(IR2TRAN,IR2DOTS,MXTRAN,MXVEC,
     &                 IR2BC,IZETVA,ITRAN,IVEC)
        NR2TRAN = MAX(NR2TRAN,ITRAN)
        MXVR2  = MAX(MXVR2,IVEC)

       ELSE 
          CALL QUIT('Error in CC_HYPPOL.')
       END IF
* K matrix transformations:
      IF (CCSLV) THEN
        CALL CC_SETB11(IKTRAN,IKDOTS,MXTRAN,MXVEC,
     &                 ITAMPA,IZETVB,IZETVC,ITRAN,IVEC)
        NKTRAN = MAX(NKTRAN,ITRAN)
        MXVAK  = MAX(MXVAK,IVEC)
        CALL CC_SETB11(IKTRAN,IKDOTS,MXTRAN,MXVEC,
     &                 ITAMPB,IZETVC,IZETVA,ITRAN,IVEC)
        NKTRAN = MAX(NKTRAN,ITRAN)
        MXVAK  = MAX(MXVAK,IVEC)
        CALL CC_SETB11(IKTRAN,IKDOTS,MXTRAN,MXVEC,
     &                 ITAMPC,IZETVA,IZETVB,ITRAN,IVEC)

        NKTRAN = MAX(NKTRAN,ITRAN)
        MXVAK  = MAX(MXVAK,IVEC)
      END IF

*---------------------------------------------------------------------*
* set up list of F{O} matrix transformations, 3 permutations
*---------------------------------------------------------------------*
        CALL CC_SETFA(IFATRAN,IFADOTS,MXTRAN,MXVEC,
     &                  0,IOPERA,IKAPA,ITAMPB,ITAMPC,ITRAN,IVEC)
        NFATRAN = MAX(NFATRAN,ITRAN)
        MXV0FA  = MAX(MXV0FA,IVEC)


        CALL CC_SETFA(IFATRAN,IFADOTS,MXTRAN,MXVEC,
     &                  0,IOPERB,IKAPB,ITAMPC,ITAMPA,ITRAN,IVEC)
        NFATRAN = MAX(NFATRAN,ITRAN)
        MXV0FA  = MAX(MXV0FA,IVEC)

        CALL CC_SETFA(IFATRAN,IFADOTS,MXTRAN,MXVEC,
     &                  0,IOPERC,IKAPC,ITAMPA,ITAMPB,ITRAN,IVEC)
        NFATRAN = MAX(NFATRAN,ITRAN)
        MXV0FA  = MAX(MXV0FA,IVEC)

*---------------------------------------------------------------------*
* set up list of ETA{O} vector calculations, 6 permutations
*---------------------------------------------------------------------*
        IF (( .NOT. USE_R2 ) .AND. (.NOT. USE_AAMAT)) THEN
          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IZETVA,IOPERB,IKAPB,0,0,0,ITAMPC,ITRAN,IVEC)
          NEATRAN = MAX(NEATRAN,ITRAN)
          MXVAEA   = MAX(MXVAEA,IVEC)

          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IZETVA,IOPERC,IKAPC,0,0,0,ITAMPB,ITRAN,IVEC)
          NEATRAN = MAX(NEATRAN,ITRAN)
          MXVAEA  = MAX(MXVAEA,IVEC)

          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IZETVB,IOPERC,IKAPC,0,0,0,ITAMPA,ITRAN,IVEC)
          NEATRAN = MAX(NEATRAN,ITRAN)
          MXVAEA  = MAX(MXVAEA,IVEC)

          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IZETVB,IOPERA,IKAPA,0,0,0,ITAMPC,ITRAN,IVEC)
          NEATRAN = MAX(NEATRAN,ITRAN)
          MXVAEA  = MAX(MXVAEA,IVEC)

          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IZETVC,IOPERA,IKAPA,0,0,0,ITAMPB,ITRAN,IVEC)
          NEATRAN = MAX(NEATRAN,ITRAN)
          MXVAEA  = MAX(MXVAEA,IVEC)

          CALL CC_SETXE('Eta',IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                  IZETVC,IOPERB,IKAPB,0,0,0,ITAMPA,ITRAN,IVEC)
          NEATRAN = MAX(NEATRAN,ITRAN)
          MXVAEA  = MAX(MXVAEA,IVEC)
        END IF

*---------------------------------------------------------------------*
* set up list of A{O} matrix transformations, 6 permutations
*---------------------------------------------------------------------*
        IF ((.NOT.USE_R2) .AND. USE_AAMAT) THEN
          IF (IKAPA.NE.0 .OR. IKAPB.NE.0 .OR. IKAPC.NE.0) THEN
            CALL QUIT('A{O} matrix transformation not yet relaxed.')
          END IF
 
          CALL CC_SETAA(IAATRAN,IAADOTS,MXTRAN,MXVEC,
     &                  IZETVA,IOPERB,ITAMPC,ITRAN,IVEC)
          NAATRAN = MAX(NAATRAN,ITRAN)
          MXVAA   = MAX(MXVAA,IVEC)

          CALL CC_SETAA(IAATRAN,IAADOTS,MXTRAN,MXVEC,
     &                  IZETVA,IOPERC,ITAMPB,ITRAN,IVEC)
          NAATRAN = MAX(NAATRAN,ITRAN)
          MXVAA   = MAX(MXVAA,IVEC)

          CALL CC_SETAA(IAATRAN,IAADOTS,MXTRAN,MXVEC,
     &                  IZETVB,IOPERC,ITAMPA,ITRAN,IVEC)
          NAATRAN = MAX(NAATRAN,ITRAN)
          MXVAA   = MAX(MXVAA,IVEC)

          CALL CC_SETAA(IAATRAN,IAADOTS,MXTRAN,MXVEC,
     &                  IZETVB,IOPERA,ITAMPC,ITRAN,IVEC)
          NAATRAN = MAX(NAATRAN,ITRAN)
          MXVAA   = MAX(MXVAA,IVEC)

          CALL CC_SETAA(IAATRAN,IAADOTS,MXTRAN,MXVEC,
     &                  IZETVC,IOPERA,ITAMPB,ITRAN,IVEC)
          NAATRAN = MAX(NAATRAN,ITRAN)
          MXVAA   = MAX(MXVAA,IVEC)

          CALL CC_SETAA(IAATRAN,IAADOTS,MXTRAN,MXVEC,
     &                  IZETVC,IOPERB,ITAMPA,ITRAN,IVEC)
          NAATRAN = MAX(NAATRAN,ITRAN)
          MXVAA   = MAX(MXVAA,IVEC)
        END IF

*---------------------------------------------------------------------*
* set up list of Xksi{O2} and Eta{O2} vector calculations, 3 permut.
*---------------------------------------------------------------------*
        IF (LPDBSA.OR.LORXA .OR. LPDBSB.OR.LORXB) THEN
          CALL CC_FIND_SO_OP(LBLOPR(IOPERA),LBLOPR(IOPERB),
     &                       LABSOP,ISYSOP,ISGNSOP,INUM,WORK,LWORK)
          IOPSOP = IROPER(LABSOP,ISYSOP)

          CALL CC_SETXE('Xi ',IXETRAN,IXDOTS,MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPA,IKAPB,0,0,IZETVC,ITRAN,IVEC)
          NXETRAN = MAX(NXETRAN,ITRAN)
          MXVXE   = MAX(MXVXE,IVEC)

          CALL CC_SETXE('Eta',IXETRAN,IEDOTS,MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPA,IKAPB,0,0,ITAMPC,ITRAN,IVEC)
          NXETRAN = MAX(NXETRAN,ITRAN)
          MXVXE   = MAX(MXVXE,IVEC)
        END IF

        IF (LPDBSA.OR.LORXA .OR. LPDBSC.OR.LORXC) THEN
          CALL CC_FIND_SO_OP(LBLOPR(IOPERA),LBLOPR(IOPERC),
     &                       LABSOP,ISYSOP,ISGNSOP,INUM,WORK,LWORK)
          IOPSOP = IROPER(LABSOP,ISYSOP)

          CALL CC_SETXE('Xi ',IXETRAN,IXDOTS,MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPA,IKAPC,0,0,IZETVB,ITRAN,IVEC)
          NXETRAN = MAX(NXETRAN,ITRAN)
          MXVXE   = MAX(MXVXE,IVEC)

          CALL CC_SETXE('Eta',IXETRAN,IEDOTS,MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPA,IKAPC,0,0,ITAMPB,ITRAN,IVEC)
          NXETRAN = MAX(NXETRAN,ITRAN)
          MXVXE   = MAX(MXVXE,IVEC)
        END IF

        IF (LPDBSC.OR.LORXC .OR. LPDBSB.OR.LORXB) THEN
          CALL CC_FIND_SO_OP(LBLOPR(IOPERC),LBLOPR(IOPERB),
     &                       LABSOP,ISYSOP,ISGNSOP,INUM,WORK,LWORK)
          IOPSOP = IROPER(LABSOP,ISYSOP)

          CALL CC_SETXE('Xi ',IXETRAN,IXDOTS,MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPC,IKAPB,0,0,IZETVA,ITRAN,IVEC)
          NXETRAN = MAX(NXETRAN,ITRAN)
          MXVXE   = MAX(MXVXE,IVEC)

          CALL CC_SETXE('Eta',IXETRAN,IEDOTS,MXTRAN,MXVEC,
     &                  0,IOPSOP,IKAPC,IKAPB,0,0,ITAMPA,ITRAN,IVEC)
          NXETRAN = MAX(NXETRAN,ITRAN)
          MXVXE   = MAX(MXVXE,IVEC)
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
      WRITE(LUPRI,'(/,/3X,A,I3,A)') 'For the requested',NQRRESF,
     &      ' quadratic response functions '
      WRITE(LUPRI,'((8X,A,I3,A))') 
     &   ' - ',NGTRAN,  ' G matrix transformations ',
     &   ' - ',NBTRAN,  ' B matrix transformations ',
     &   ' - ',NKTRAN,  ' K matrix transformations ',
     &   ' - ',NFTRAN,  ' F matrix transformations ',
     &   ' - ',NFATRAN, ' F{O} matrix transformations ',
     &   ' - ',NEATRAN, ' ETA{O} vector calculations ',
     &   ' - ',NXETRAN, ' Xksi{O2} x L1 and Eta{O2} x R1 dot products',
     &   ' - ',NR2TRAN, ' X1 x R2 dot products ' 
      WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'


* G matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G matrix transformations:'
        DO ITRAN = 1, NGTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IGTRAN(I,ITRAN),I=1,MXDIM_GTRAN),
     &     (IGDOTS(I,ITRAN),I=1,MXV0G)
        END DO
        WRITE (LUPRI,*)
      END IF

* B matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of B matrix transformations:'
        DO ITRAN = 1, NBTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IBTRAN(I,ITRAN),I=1,MXDIM_BTRAN),
     &     (IBDOTS(I,ITRAN),I=1,MXVAB)
        END DO
        WRITE (LUPRI,*)
      END IF

* K matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of K matrix transformations:'
        DO ITRAN = 1, NKTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IKTRAN(I,ITRAN),I=1,MXDIM_KTRAN),
     &     (IKDOTS(I,ITRAN),I=1,MXVAK)
        END DO
        WRITE (LUPRI,*)
      END IF

* F matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F matrix transformations:'
        DO ITRAN = 1, NFTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IFTRAN(I,ITRAN),I=1,MXDIM_FTRAN),
     &     (IFDOTS(I,ITRAN),I=1,MXVAF)
        END DO
        WRITE (LUPRI,*)
      END IF

* F{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F{O} matrix transformations:'
        DO ITRAN = 1, NFATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (IFATRAN(I,ITRAN),I=1,5),(IFADOTS(I,ITRAN),I=1,MXV0FA)
        END DO
        WRITE (LUPRI,*)
      END IF

* ETA{O} vector calculations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of ETA{O} vector calculations:'
        DO ITRAN = 1, NEATRAN
          WRITE(LUPRI,'(A,8I5,5X,(12I5,20X))') MSGDBG,
     &     (IEATRAN(I,ITRAN),I=1,MXDIM_XEVEC),
     &     (IEADOTS(I,ITRAN),I=1,MXVAEA)
        END DO
        WRITE (LUPRI,*)
      END IF

* A{O} vector calculations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of A{O} vector calculations:'
        DO ITRAN = 1, NAATRAN
          WRITE(LUPRI,'(A,8I5,5X,(12I5,20X))') MSGDBG,
     &     (IAATRAN(I,ITRAN),I=1,MXDIM_AATRAN),
     &     (IAADOTS(I,ITRAN),I=1,MXVAA)
        END DO
        WRITE (LUPRI,*)
      END IF

* Xksi{O2} x L1 and Eta{O2} x R1 dot products:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) ' Xksi{O2} x L1 and Eta{O2} x R1 dot products:'
        DO ITRAN = 1, NXETRAN
          WRITE(LUPRI,'(A,8I5,5X,(12I5,20X))') MSGDBG,
     &     (IXETRAN(I,ITRAN),I=1,8),(IXDOTS(I,ITRAN),I=1,MXVXE)
          WRITE(LUPRI,'(A,40X,5X,(12I5,20X))') MSGDBG,
     &                              (IEDOTS(I,ITRAN),I=1,MXVXE)
        END DO
        WRITE (LUPRI,*)
      END IF

* X1 x R2 dot products
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of X1 x R2 dot products:'
        DO ITRAN = 1, NR2TRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &      IR2TRAN(ITRAN),(IR2DOTS(I,ITRAN),I=1,MXVR2)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF


      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQR_SETUP                           *
*---------------------------------------------------------------------*
c /* deck CCQR_DISP_SETUP */
*=====================================================================*
      SUBROUTINE CCQR_DISP_SETUP(MXTRAN, MXVEC,
     &                           IGTRAN, IGDOTS, WG, NGTRAN,
     &                           IFTRAN, IFDOTS, WF, NFTRAN,
     &                           IBTRAN, IBDOTS, WB, NBTRAN,
     &                           IFATRAN,IFADOTS,WFA,NFATRAN,
     &                           IEATRAN,IEADOTS,WEA,NEATRAN,
     &                           IR2TRAN,IR2DOTS,WR2,NR2TRAN,
     &                           DISPCF, MXDSPCF,LADD        )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CCQR dispersion coefficients:
*                - list of G matrix transformations 
*                - list of F/B matrix transformations 
*                - list of F{O} matrix transformations 
*                - list of ETA{O} vector calculations 
*                - list of XC x CR2 dot products      
*
*     Written by Christof Haettig, october 1997.
*     based on CCQR_SETUP routine, B matrix transf. added 12/97
*     USE_R2 option introduced 3/98                               
*     odd dispersion coefficients 8/98
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccqrinf.h>
#include <ccroper.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(25) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCQR_DISP_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LADD

      INTEGER MXVEC, MXTRAN, MXDSPCF

      INTEGER IGTRAN(4,MXTRAN)
      INTEGER IGDOTS(MXVEC,MXTRAN)

      INTEGER IFTRAN(3,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)
  
      INTEGER IBTRAN(3,MXTRAN)
      INTEGER IBDOTS(MXVEC,MXTRAN)
  
      INTEGER IFATRAN(MXDIM_FATRAN,MXTRAN)
      INTEGER IFADOTS(MXVEC,MXTRAN)

      INTEGER IEATRAN(3,MXTRAN)
      INTEGER IEADOTS(MXVEC,MXTRAN)

      INTEGER IR2TRAN(MXTRAN)
      INTEGER IR2DOTS(MXVEC,MXTRAN)

      INTEGER NGTRAN,NFTRAN,NBTRAN,NFATRAN,NEATRAN,NQRDSPCF,NR2TRAN

      INTEGER ISYMA, ISYMB, ISYMC, ISYML, ISYML1, ISYML2
      INTEGER IDISP, IOPER, IDX

      INTEGER MXVG, MXVF, MXVB, MXVFA, MXVEA, MXVR2
      INTEGER IOPERA, IOPERB, IOPERC
      INTEGER IZETVA, IZETVB, IZETVC
      INTEGER ITAMPA, ITAMPB, ITAMPC
      INTEGER IR2AB, IR2AC, IR2BC, ICTOT, ISACAU
      INTEGER IVEC, ITRAN, I, ISAMA, ISAMB, ISAMC, ISAPROP

#if defined (SYS_CRAY)
      REAL DISPCF(MXDSPCF)
      REAL WG(MXVEC,MXTRAN)
      REAL WF(MXVEC,MXTRAN)
      REAL WB(MXVEC,MXTRAN)
      REAL WFA(MXVEC,MXTRAN)
      REAL WEA(MXVEC,MXTRAN)
      REAL WR2(MXVEC,MXTRAN)
      REAL GCON, FCON(3), BCON(3), FACON(3), EACON(6), R2CON(3)
      REAL ZERO
#else
      DOUBLE PRECISION DISPCF(MXDSPCF)
      DOUBLE PRECISION WG(MXVEC,MXTRAN)
      DOUBLE PRECISION WF(MXVEC,MXTRAN)
      DOUBLE PRECISION WB(MXVEC,MXTRAN)
      DOUBLE PRECISION WFA(MXVEC,MXTRAN)
      DOUBLE PRECISION WEA(MXVEC,MXTRAN)
      DOUBLE PRECISION WR2(MXVEC,MXTRAN)
      DOUBLE PRECISION GCON,FCON(3),BCON(3),FACON(3),EACON(6),R2CON(3)
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)


* external functions:
      INTEGER ILRCAMP
      INTEGER ILC1AMP
      INTEGER ICR2AMP


*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN
        NGTRAN  = 0
        NFTRAN  = 0
        NBTRAN  = 0
        NFATRAN = 0
        NEATRAN = 0
        NR2TRAN = 0
      END IF

      MXVG  = 0
      MXVF  = 0
      MXVB  = 0
      MXVFA = 0
      MXVEA = 0
      MXVR2 = 0
 
      NQRDSPCF = 0

      CALL DZERO(DISPCF,MXDSPCF)

*---------------------------------------------------------------------*
* start loop over all requested dispersion coefficients:
*---------------------------------------------------------------------*
 
      DO IOPER = 1, NQROPER
        IOPERA = IAQROP(IOPER)
        IOPERB = IBQROP(IOPER)
        IOPERC = ICQROP(IOPER)

        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        ISYMC  = ISYOPR(IOPERC)

        ISAMA  = ISYMAT(IOPERA)
        ISAMB  = ISYMAT(IOPERB)
        ISAMC  = ISYMAT(IOPERC)

        ISAPROP = ISAMA * ISAMB * ISAMC

      IF (MULD2H(ISYMB,ISYMC).EQ.ISYMA) THEN

      DO IDISP = 1, NQRDISP

        ICTOT  = IQCAUA(IDISP) + IQCAUB(IDISP) + IQCAUC(IDISP)
        ISACAU = 2 * ( (ICTOT/2)*2 - ICTOT ) + 1

      IF (ISACAU.EQ.ISAPROP .OR. ISAPROP.EQ.0 .OR. ALLDSPCF) THEN

        NQRDSPCF = NQRDSPCF + 1

        IZETVA = ILC1AMP(LBLOPR(IOPERA),IQCAUA(IDISP),ISYML)
        IZETVB = ILC1AMP(LBLOPR(IOPERB),IQCAUB(IDISP),ISYML)
        IZETVC = ILC1AMP(LBLOPR(IOPERC),IQCAUC(IDISP),ISYML)

        ITAMPA = ILRCAMP(LBLOPR(IOPERA),IQCAUA(IDISP),ISYML)
        ITAMPB = ILRCAMP(LBLOPR(IOPERB),IQCAUB(IDISP),ISYML)
        ITAMPC = ILRCAMP(LBLOPR(IOPERC),IQCAUC(IDISP),ISYML)

        IF ( USE_R2 ) THEN
          IR2AB = ICR2AMP(LBLOPR(IOPERA),IQCAUA(IDISP),ISYML1,
     &                    LBLOPR(IOPERB),IQCAUB(IDISP),ISYML2)
          IR2AC = ICR2AMP(LBLOPR(IOPERA),IQCAUA(IDISP),ISYML1,
     &                    LBLOPR(IOPERC),IQCAUC(IDISP),ISYML2)
          IR2BC = ICR2AMP(LBLOPR(IOPERB),IQCAUB(IDISP),ISYML1,
     &                    LBLOPR(IOPERC),IQCAUC(IDISP),ISYML2)
        END IF



*---------------------------------------------------------------------*
* set up list of G matrix transformations, 1 permutation
*---------------------------------------------------------------------*
        CALL CCQR_SETG(IGTRAN,IGDOTS,MXTRAN,MXVEC,
     &                 0,ITAMPA,ITAMPB,ITAMPC,ITRAN,IVEC)
        NGTRAN = MAX(NGTRAN,ITRAN)
        MXVG   = MAX(MXVG,IVEC)
        GCON   = WG(IVEC,ITRAN)

*---------------------------------------------------------------------*
* set up list of B/F matrix transformations, 3 permutations
* or if USE_R2 is set a list of XC x CR2 dot products, 3 permutations
*---------------------------------------------------------------------*
        IF ( USEBTRAN .AND. (.NOT.USE_R2) ) THEN
          CALL CC_SETB11(IBTRAN,IBDOTS,MXTRAN,MXVEC,
     &                   IZETVA,ITAMPB,ITAMPC,ITRAN,IVEC)
          NBTRAN  = MAX(NBTRAN,ITRAN)
          MXVB    = MAX(MXVB,IVEC)
          BCON(1) = WB(IVEC,ITRAN)

          CALL CC_SETB11(IBTRAN,IBDOTS,MXTRAN,MXVEC,
     &                   IZETVB,ITAMPC,ITAMPA,ITRAN,IVEC)
          NBTRAN  = MAX(NBTRAN,ITRAN)
          MXVB    = MAX(MXVB,IVEC)
          BCON(2) = WB(IVEC,ITRAN)

          CALL CC_SETB11(IBTRAN,IBDOTS,MXTRAN,MXVEC,
     &                   IZETVC,ITAMPA,ITAMPB,ITRAN,IVEC)
          NBTRAN  = MAX(NBTRAN,ITRAN)
          MXVB    = MAX(MXVB,IVEC)
          BCON(3) = WB(IVEC,ITRAN)

          FCON(1)  = 0.0d0
          FCON(2)  = 0.0d0
          FCON(3)  = 0.0d0
          R2CON(1) = 0.0d0
          R2CON(2) = 0.0d0
          R2CON(3) = 0.0d0

        ELSE IF ( (.NOT.USEBTRAN) .AND. (.NOT.USE_R2) ) THEN

          CALL CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                   IZETVA,ITAMPB,ITAMPC,ITRAN,IVEC)
          NFTRAN  = MAX(NFTRAN,ITRAN)
          MXVF    = MAX(MXVF,IVEC)
          FCON(1) = WF(IVEC,ITRAN)

          CALL CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                   IZETVB,ITAMPC,ITAMPA,ITRAN,IVEC)
          NFTRAN  = MAX(NFTRAN,ITRAN)
          MXVF    = MAX(MXVF,IVEC)
          FCON(2) = WF(IVEC,ITRAN)

          CALL CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                   IZETVC,ITAMPA,ITAMPB,ITRAN,IVEC)
          NFTRAN  = MAX(NFTRAN,ITRAN)
          MXVF    = MAX(MXVF,IVEC)
          FCON(3) = WF(IVEC,ITRAN)

          BCON(1) = 0.0d0
          BCON(2) = 0.0d0
          BCON(3) = 0.0d0
          R2CON(1) = 0.0d0
          R2CON(2) = 0.0d0
          R2CON(3) = 0.0d0

        ELSE IF ( USE_R2 ) THEN

          CALL CC_SETDOT(IR2TRAN,IR2DOTS,MXTRAN,MXVEC,
     &                   IR2AB,IZETVC,ITRAN,IVEC)
          NR2TRAN  = MAX(NR2TRAN,ITRAN)
          MXVR2    = MAX(MXVR2,IVEC)
          R2CON(1) = WR2(IVEC,ITRAN)

          CALL CC_SETDOT(IR2TRAN,IR2DOTS,MXTRAN,MXVEC,
     &                   IR2AC,IZETVB,ITRAN,IVEC)
          NR2TRAN  = MAX(NR2TRAN,ITRAN)
          MXVR2    = MAX(MXVR2,IVEC)
          R2CON(2) = WR2(IVEC,ITRAN)

          CALL CC_SETDOT(IR2TRAN,IR2DOTS,MXTRAN,MXVEC,
     &                   IR2BC,IZETVA,ITRAN,IVEC)
          NR2TRAN  = MAX(NR2TRAN,ITRAN)
          MXVR2    = MAX(MXVR2,IVEC)
          R2CON(3) = WR2(IVEC,ITRAN)

          FCON(1)  = 0.0d0
          FCON(2)  = 0.0d0
          FCON(3)  = 0.0d0
          BCON(1) = 0.0d0
          BCON(2) = 0.0d0
          BCON(3) = 0.0d0

        ELSE
          CALL QUIT('Error in CC_HYPPOL.')
        END IF

*---------------------------------------------------------------------*
* set up list of F{O} matrix transformations, 3 permutations
*---------------------------------------------------------------------*
        IF ( IQCAUA(IDISP).EQ.0 ) THEN
          CALL CCQR_SETFA(IFATRAN,IFADOTS,MXTRAN,MXVEC,
     &                    0,IOPERA,ITAMPB,ITAMPC,ITRAN,IVEC)
          NFATRAN  = MAX(NFATRAN,ITRAN)
          MXVFA    = MAX(MXVFA,IVEC)
          FACON(1) = WFA(IVEC,ITRAN)
        ELSE
          FACON(1) = ZERO
        END IF

        IF ( IQCAUB(IDISP).EQ.0 ) THEN
          CALL CCQR_SETFA(IFATRAN,IFADOTS,MXTRAN,MXVEC,
     &                    0,IOPERB,ITAMPC,ITAMPA,ITRAN,IVEC)
          NFATRAN  = MAX(NFATRAN,ITRAN)
          MXVFA    = MAX(MXVFA,IVEC)
          FACON(2) = WFA(IVEC,ITRAN)
        ELSE
          FACON(2) = ZERO
        END IF

        IF ( IQCAUC(IDISP).EQ.0 ) THEN
          CALL CCQR_SETFA(IFATRAN,IFADOTS,MXTRAN,MXVEC,
     &                    0,IOPERC,ITAMPA,ITAMPB,ITRAN,IVEC)
          NFATRAN  = MAX(NFATRAN,ITRAN)
          MXVFA    = MAX(MXVFA,IVEC)
          FACON(3) = WFA(IVEC,ITRAN)
        ELSE
          FACON(3) = ZERO
        END IF

*---------------------------------------------------------------------*
* set up list of ETA{O} vector calculations, 6 permutations
*---------------------------------------------------------------------*
        IF ( IQCAUB(IDISP).EQ.0 .AND. (.NOT.USE_R2) ) THEN
          CALL CCQR_SETEA(IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                    IZETVA,IOPERB,ITAMPC,ITRAN,IVEC)
          NEATRAN  = MAX(NEATRAN,ITRAN)
          MXVEA    = MAX(MXVEA,IVEC)
          EACON(1) = WEA(IVEC,ITRAN)
        ELSE 
          EACON(1) = ZERO
        END IF

        IF ( IQCAUC(IDISP).EQ.0 .AND. (.NOT.USE_R2) ) THEN
          CALL CCQR_SETEA(IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                    IZETVA,IOPERC,ITAMPB,ITRAN,IVEC)
          NEATRAN  = MAX(NEATRAN,ITRAN)
          MXVEA    = MAX(MXVEA,IVEC)
          EACON(2) = WEA(IVEC,ITRAN)
        ELSE 
          EACON(2) = ZERO
        END IF

        IF ( IQCAUC(IDISP).EQ.0 .AND. (.NOT.USE_R2) ) THEN
          CALL CCQR_SETEA(IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                    IZETVB,IOPERC,ITAMPA,ITRAN,IVEC)
          NEATRAN  = MAX(NEATRAN,ITRAN)
          MXVEA    = MAX(MXVEA,IVEC)
          EACON(3) = WEA(IVEC,ITRAN)
        ELSE 
          EACON(3) = ZERO
        END IF

        IF ( IQCAUA(IDISP).EQ.0 .AND. (.NOT.USE_R2) ) THEN
          CALL CCQR_SETEA(IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                    IZETVB,IOPERA,ITAMPC,ITRAN,IVEC)
          NEATRAN  = MAX(NEATRAN,ITRAN)
          MXVEA    = MAX(MXVEA,IVEC)
          EACON(4) = WEA(IVEC,ITRAN)
        ELSE 
          EACON(4) = ZERO
        END IF

        IF ( IQCAUA(IDISP).EQ.0 .AND. (.NOT.USE_R2) ) THEN
          CALL CCQR_SETEA(IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                    IZETVC,IOPERA,ITAMPB,ITRAN,IVEC)
          NEATRAN  = MAX(NEATRAN,ITRAN)
          MXVEA    = MAX(MXVEA,IVEC)
          EACON(5) = WEA(IVEC,ITRAN)
        ELSE 
          EACON(5) = ZERO
        END IF

        IF ( IQCAUB(IDISP).EQ.0 .AND. (.NOT.USE_R2) ) THEN
          CALL CCQR_SETEA(IEATRAN,IEADOTS,MXTRAN,MXVEC,
     &                    IZETVC,IOPERB,ITAMPA,ITRAN,IVEC)
          NEATRAN  = MAX(NEATRAN,ITRAN)
          MXVEA    = MAX(MXVEA,IVEC)
          EACON(6) = WEA(IVEC,ITRAN)
        ELSE 
          EACON(6) = ZERO
        END IF


*---------------------------------------------------------------------*
* add all 13 contributions to the dispersion coefficients:
*---------------------------------------------------------------------*
        IF (LADD) THEN
          IDX = NQRDISP*(IOPER-1) + IDISP
          
          DISPCF(IDX) = GCON + 
     &                  FCON(1)  + FCON(2)  + FCON(3)  +
     &                  BCON(1)  + BCON(2)  + BCON(3)  +
     &                  R2CON(1) + R2CON(2) + R2CON(3) +
     &                  FACON(1) + FACON(2) + FACON(3) +
     &                  EACON(1) + EACON(2) + EACON(3) +
     &                  EACON(4) + EACON(5) + EACON(6)  

          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'ITAMPA,ITAMPB,ITAMPC:',ITAMPA,ITAMPB,ITAMPC
            WRITE(LUPRI,*) 'IZETVA,IZETVB,IZETVC:',IZETVA,IZETVB,IZETVC
            WRITE(LUPRI,*) 'IOPERA,IOPERB,IOPERC:',IOPERA,IOPERB,IOPERC
            WRITE(LUPRI,*) 'GCON :',GCON
            WRITE(LUPRI,*) 'FCON :',FCON(1),FCON(2),FCON(3)
            WRITE(LUPRI,*) 'BCON :',BCON(1),BCON(2),BCON(3)
            WRITE(LUPRI,*) 'R2CON:',R2CON(1),R2CON(2),R2CON(3)
            WRITE(LUPRI,*) 'FACON:',FACON(1),FACON(2),FACON(3)
            WRITE(LUPRI,*) 'EACON:',EACON(1),EACON(2),EACON(3),
     &                              EACON(4),EACON(5),EACON(6)
            WRITE(LUPRI,*) 'DISPCF:',DISPCF(IDX)
          END IF

        END IF

*---------------------------------------------------------------------*
* end loop over all requested dispersion coefficients
*---------------------------------------------------------------------*
      END IF
      END DO
      END IF
      END DO

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
      IF (.NOT. LADD) THEN

* general statistics:
      WRITE(LUPRI,'(/,/3X,A,I4,A)') 'For the requested',NQRDSPCF,
     &      ' expansion coefficients for quadratic response functions'
      WRITE(LUPRI,'((8X,A,I3,A))') 
     &   ' - ',NGTRAN,  ' G matrix transformations ',
     &   ' - ',NBTRAN,  ' B matrix transformations ',
     &   ' - ',NFTRAN,  ' F matrix transformations ',
     &   ' - ',NFATRAN, ' F{O} matrix transformations ',
     &   ' - ',NEATRAN, ' ETA{O} vector calculations ',
     &   ' - ',NR2TRAN, ' XC x CR2 dot products ' 
      WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'


* G matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G matrix transformations:'
        DO ITRAN = 1, NGTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IGTRAN(I,ITRAN),I=1,3),(IGDOTS(I,ITRAN),I=1,MXVG)
        END DO
        WRITE (LUPRI,*)
      END IF

* B matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of B matrix transformations:'
        DO ITRAN = 1, NBTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IBTRAN(I,ITRAN),I=1,2),(IBDOTS(I,ITRAN),I=1,MXVB)
        END DO
        WRITE (LUPRI,*)
      END IF

* F matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F matrix transformations:'
        DO ITRAN = 1, NFTRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IFTRAN(I,ITRAN),I=1,2),(IFDOTS(I,ITRAN),I=1,MXVF)
        END DO
        WRITE (LUPRI,*)
      END IF

* F{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F{O} matrix transformations:'
        DO ITRAN = 1, NFATRAN
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') MSGDBG,
     &     (IFATRAN(I,ITRAN),I=1,5),(IFADOTS(I,ITRAN),I=1,MXVFA)
        END DO
        WRITE (LUPRI,*)
      END IF

* ETA{O} vector calculations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of ETA{O} vector calculations:'
        DO ITRAN = 1, NEATRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IEATRAN(I,ITRAN),I=1,2),(IEADOTS(I,ITRAN),I=1,MXVEA)
        END DO
        WRITE (LUPRI,*)
      END IF

* X1 x R2 dot products
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of XC x CR2 dot products:'
        DO ITRAN = 1, NR2TRAN
          WRITE(LUPRI,'(A,I5,5X,(12I5,20X))') MSGDBG,
     &      IR2TRAN(ITRAN),(IR2DOTS(I,ITRAN),I=1,MXVR2)
        END DO
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

      CALL FLSHFO(LUPRI)

      END IF ! (.NOT. LADD)

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQR_DISP_SETUP                      *
*---------------------------------------------------------------------*
c /* deck ccqr_setg */
*=====================================================================*
      SUBROUTINE CCQR_SETG(IGTRAN,IGDOTS,MXTRAN,MXVEC,
     &                     IZETAV,ITAMPA,ITAMPB,ITAMPC,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintain a list of dot products of G matrix 
*             transformations with right amplitude vectors:  
*                       (Z*G*T*T) * T
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
*    Written by Christof Haettig, november 1996.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IGTRAN(4,MXTRAN)
      INTEGER IGDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDA, LFNDB, LFNDC
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
        LFNDA = LGTST(ITRAN,IZETAV,ITAMPB,ITAMPC)
        LFNDB = LGTST(ITRAN,IZETAV,ITAMPC,ITAMPA)
        LFNDC = LGTST(ITRAN,IZETAV,ITAMPA,ITAMPB)

        DO WHILE ( .NOT. (LFNDA.OR.LFNDB.OR.LFNDC.OR.LGEND(ITRAN)) )
         ITRAN = ITRAN + 1
         LFNDA = LGTST(ITRAN,IZETAV,ITAMPB,ITAMPC)
         LFNDB = LGTST(ITRAN,IZETAV,ITAMPC,ITAMPA)
         LFNDC = LGTST(ITRAN,IZETAV,ITAMPA,ITAMPB)
        END DO

        IF (.NOT.(LFNDA.OR.LFNDB.OR.LFNDC)) THEN
          IGTRAN(1,ITRAN) = IZETAV
          IGTRAN(2,ITRAN) = ITAMPA
          IGTRAN(3,ITRAN) = ITAMPB
          IGTRAN(4,ITRAN) = 0
          ITAMP = ITAMPC
        ELSE 
          IF (LFNDA) ITAMP = ITAMPA
          IF (LFNDB) ITAMP = ITAMPB
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
            WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') 'CCQR_SETG>',
     &         (IGTRAN(I,IDX),I=1,4),(IGDOTS(I,IDX),I=1,MXVEC)
            IDX = IDX + 1
          END DO
          CALL FLSHFO(LUPRI)
          CALL QUIT('Overflow error in CCQR_SETG.')
        END IF
      
        RETURN
        END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQR_SETG                            *
*---------------------------------------------------------------------*
c /* deck ccqr_setf */
*=====================================================================*
       SUBROUTINE CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                      IZETAV,ITAMPA,ITAMPB,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: set up list of F matrix transformations for CCQR section
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
*    Written by Christof Haettig, november 1996.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IFTRAN(3,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDA, LFNDB
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
      LFNDA = LFTST(ITRAN,IZETAV,ITAMPB)
      LFNDB = LFTST(ITRAN,IZETAV,ITAMPA)

      DO WHILE ( .NOT. (LFNDA.OR.LFNDB.OR.LFEND(ITRAN)) )
       ITRAN = ITRAN + 1
       LFNDA = LFTST(ITRAN,IZETAV,ITAMPB)
       LFNDB = LFTST(ITRAN,IZETAV,ITAMPA)
      END DO

      IF (.NOT.(LFNDA.OR.LFNDB)) THEN
        IFTRAN(1,ITRAN) = IZETAV
        IFTRAN(2,ITRAN) = ITAMPA
        IFTRAN(3,ITRAN) = 0
        ITAMP = ITAMPB
      ELSE 
        IF (LFNDA) ITAMP = ITAMPA
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
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') 'CCQR_SETF>',
     &       (IFTRAN(I,IDX),I=1,3),(IFDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CCQR_SETF.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQR_SETF                            *
*---------------------------------------------------------------------*
c /* deck ccqr_setfa */
*=====================================================================*
       SUBROUTINE CCQR_SETFA(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                       IZETAV,IOPER,ITAMPA,ITAMPB,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: set up list of F matrix transformations for CCQR section
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
*    Written by Christof Haettig, november 1996.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IFTRAN(5,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDA, LFNDB
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
      LFNDA  = LFATST(ITRAN,IZETAV,IOPER,ITAMPB)
      LFNDB  = LFATST(ITRAN,IZETAV,IOPER,ITAMPA)

      DO WHILE ( .NOT. (LFNDA.OR.LFNDB.OR.LFAEND(ITRAN)))
       ITRAN = ITRAN + 1
       LFNDA  = LFATST(ITRAN,IZETAV,IOPER,ITAMPB)
       LFNDB  = LFATST(ITRAN,IZETAV,IOPER,ITAMPA)
      END DO

      IF (.NOT.(LFNDA.OR.LFNDB)) THEN
        IFTRAN(1,ITRAN) = IZETAV
        IFTRAN(2,ITRAN) = IOPER
        IFTRAN(3,ITRAN) = ITAMPA
        IFTRAN(4,ITRAN) = 0
        IFTRAN(5,ITRAN) = 0
        ITAMP = ITAMPB
      ELSE 
        IF (LFNDA) ITAMP = ITAMPA
        IF (LFNDB) ITAMP = ITAMPB
      END IF

      IVEC = 1
      DO WHILE (IFDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &            IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IFDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'Overflow error in CCQR_SETFA:'
        WRITE (LUPRI,*) 'IVEC, MXVEC :',IVEC, MXVEC
        WRITE (LUPRI,*) 'ITRAN, MXTRAN:',ITRAN, MXTRAN
        WRITE (LUPRI,*) 'ITAMPA,ITAMPB:',ITAMPA,ITAMPB
        IDX = 1
        DO WHILE ( .NOT. LFAEND(IDX) )
          WRITE(LUPRI,'(A,4I5,5X,(12I5,20X))') 'CCQR_SETFA>',
     &       (IFTRAN(I,IDX),I=1,4),(IFDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CCQR_SETFA')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQR_SETFA                           *
*---------------------------------------------------------------------*
c /* deck ccqr_setea */
*=====================================================================*
       SUBROUTINE CCQR_SETEA(IETRAN,IEDOTS,MXTRAN,MXVEC,
     &                       IZETAV,IOPER,ITAMP,ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: set up list of ETA vectors for CCQR section
*
*             IETRAN - list of F matrix transformations
*             IEDOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum IETRAN list dimension
*             MXVEC  - second maximum dimension maximum for IEDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             IOPER  - index of property operator 
*             ITAMP  - index of amplitude vector
*
*             ITRAN - index in IETRAN list
*             IVEC  - second index in IEDOTS list
*
*    Written by Christof Haettig, november 1996.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IETRAN(3,MXTRAN)
      INTEGER IEDOTS(MXVEC,MXTRAN)

      LOGICAL LFND
      INTEGER IZETAV, IOPER, ITAMP
      INTEGER ITRAN, IVEC
      INTEGER I, IDX

* statement  functions:
      LOGICAL LEATST, LEAEND
      INTEGER IL, IA, IO
      LEATST(ITRAN,IL,IO,IA) = 
     &       IETRAN(1,ITRAN).EQ.IL .AND. IETRAN(2,ITRAN).EQ.IO 
      LEAEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &  (IETRAN(1,ITRAN)+IETRAN(2,ITRAN)).LE.0

*---------------------------------------------------------------------*
* maintain list of ETA{A} vectors
*---------------------------------------------------------------------*
      ITRAN = 1
      LFND  = LEATST(ITRAN,IZETAV,IOPER,ITAMP)

      DO WHILE ( .NOT.(LFND.OR.LEAEND(ITRAN))) 
       ITRAN = ITRAN + 1
       LFND  = LEATST(ITRAN,IZETAV,IOPER,ITAMP)
      END DO

      IF (.NOT.LFND) THEN
        IETRAN(1,ITRAN) = IZETAV
        IETRAN(2,ITRAN) = IOPER
        IETRAN(3,ITRAN) = 0
      END IF

      IVEC = 1
      DO WHILE (IEDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &            IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IEDOTS(IVEC,ITRAN) = ITAMP

*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'ITAMP:',ITAMP
        IDX = 1
        DO WHILE( .NOT. LEAEND(IDX) )
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') 'CCQR_SETEA>',
     &       (IETRAN(I,IDX),I=1,3),(IEDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CCQR_SETEA.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQR_SETEA                           *
*---------------------------------------------------------------------*
c /* deck ccqr_fadrv */
*=====================================================================*
      SUBROUTINE CCQR_FADRV(LISTL, LISTO, LISTB, LISTC,
     &                      NFATRAN, MXVEC,
     &                      IFATRAN, IFADOTS, FACON,
     &                      WORK, LWORK, OPTION)
*---------------------------------------------------------------------*
*
*    Purpose: driver for F{A} matrix transformations
*
*    Written by Christof Haettig, November 1996.
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccroper.h>
#include <ccr1rsp.h>
#include <cclists.h>
#include <ccnoddy.h>
#include <ccfield.h>

* local parameters:
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCQR_FADRV> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER LUFBTA

      INTEGER LWORK, NFATRAN, MXVEC

      CHARACTER*8 FNTOC
      CHARACTER*4 FNDKBC
      CHARACTER*5 FNDKBC3
      CHARACTER*6 FNDELD, FNCKJD, FN3VI, FN3FOPX
      CHARACTER*7 FN3FOP2X
 
      PARAMETER (FNDELD  = 'CKDELD'  , FNCKJD  = 'CKJDEL'  )
      PARAMETER (FNDKBC  = 'DKBC'    , FNTOC   = 'CCSDT_OC')
      PARAMETER (FN3VI   = 'CC3_VI'  )
      PARAMETER (FN3FOPX = 'PTFOPX'  , FN3FOP2X= 'PTFOP2X' )
      PARAMETER (FNDKBC3 = 'DKBC3'   )  

#if defined (SYS_CRAY)
      REAL FACON(MXVEC,NFATRAN)
      REAL WORK(LWORK)
      REAL ZERO
      REAL FACON1, FACON2
      REAL DDOT
#else
      DOUBLE PRECISION FACON(MXVEC,NFATRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION FACON1, FACON2
      DOUBLE PRECISION DDOT
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER LISTL*(*), LISTO*(*), LISTB*(*), LISTC*(*), OPTION*(*)
      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABELA, FILFBTA
      LOGICAL LPDBSA, FADEN_NODDY, NODDY_CCSDT
      INTEGER IFATRAN(MXDIM_FATRAN,NFATRAN)
      INTEGER IFADOTS(MXVEC,NFATRAN)
      INTEGER ISYRES, ISYMTA, ISYMTB, ISYMTC, ISYCTR, IFILE, IRELAX
      INTEGER IVEC, ITRAN, IOPERA, ITAMPB, ITAMPC, IZETAV, IOPT
      INTEGER KGAMMA1, KGAMMA2, KT1AMPC, KT2AMPC, KEND1, LWRK1
      INTEGER KEND1A, LEND1A, IADR,  LEN, IOPTRES
      INTEGER IFBTRAN(5)
      INTEGER MXFADEN, MXINTER, KIDXL_FADEN, KIDXB_FADEN, KIDXC_FADEN,
     &        KIADR_FADEN, KIDXL_INTER, KIDXR_INTER, KLSTR_INTER,
     &        KIADR_INTER, NFADEN, NINTER


* external functions:
      INTEGER ILSTSYM

      CALL QENTER('CCQR_FADRV')

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (MXVEC*NFATRAN.NE.0) CALL DZERO(FACON,MXVEC*NFATRAN)

      IF (CCSDT .AND. NONHF .AND. NFIELD.GT.0 .AND. (.NOT.NODDY_FAMAT))
     &  CALL QUIT('No finite field in CC3 for F{O} module.')

      FADEN_NODDY = .FALSE.
      IF (CCSDT .AND. OPTION(1:4).EQ.'DOTP' .AND. 
     &    NODDY_FAMAT .AND. NODDY_FA_ALTER) THEN
        ! we compute 'dot products' for CC3, which is via densities
        ! and the noddy flag is switched on and we want to use
        ! the noddy code for the densities (not for the transformation)
        ! check for which cases this is possible:
        IF ( (LISTO(1:3).EQ.'o1 ') 
     &       .AND. (LISTL(1:3).EQ.'L0 '.OR.
     &              LISTL(1:3).EQ.'L1 '.OR.LISTL(1:3).EQ.'LE ') 
     &       .AND. (LISTB(1:3).EQ.'R1 '.OR.LISTB(1:3).EQ.'RE '.OR.
     &              LISTB(1:3).EQ.'R2 '.OR.LISTB(1:3).EQ.'ER1'    ) 
     &       .AND. (LISTC(1:3).EQ.'R1 '.OR.LISTC(1:3).EQ.'RE '.OR.
     &              LISTC(1:3).EQ.'R2 '.OR.LISTC(1:3).EQ.'ER1'    ) )
     &    FADEN_NODDY = .TRUE. 
      END IF

*---------------------------------------------------------------------*
* start loop over all requested F{O} matrix transformations
*---------------------------------------------------------------------*
      DO ITRAN = 1, NFATRAN

        IZETAV = IFATRAN(1,ITRAN)
        IOPERA = IFATRAN(2,ITRAN)
        ITAMPB = IFATRAN(3,ITRAN)
        IFILE  = IFATRAN(4,ITRAN)
        IRELAX = IFATRAN(5,ITRAN)

        LABELA = LBLOPR(IOPERA)
        LPDBSA = LPDBSOP(IOPERA)


        ISYCTR = ILSTSYM(LISTL,IZETAV)
        ISYMTA = ILSTSYM(LISTO,IOPERA)
        ISYMTB = ILSTSYM(LISTB,ITAMPB)

        ISYRES = MULD2H(MULD2H(ISYMTA,ISYMTB),ISYCTR)

        KGAMMA1 = 1
        KGAMMA2 = KGAMMA1 + NT1AM(ISYRES)
        KEND1A  = KGAMMA2 + NT2AM(ISYRES)
        LEND1A  = LWORK - KEND1A

        KT1AMPC = KEND1A
        KT2AMPC = KT1AMPC + NT1AM(ISYRES)
        KEND1   = KT2AMPC + NT2AM(ISYRES)
        LWRK1   = LWORK - KEND1


        IF (LWRK1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCQR_FADRV.')
        END IF

*---------------------------------------------------------------------*
* calculate F{A} matrix times a response amplitude vector:
*      F{A} * t^B = <Lambda| [[A,T^B],tau] |CC>
*---------------------------------------------------------------------*
       IF ( (IRELAX.GT.0) .OR. LPDBSA ) THEN

          IF (CC2) THEN

            CALL CC_FDFBMAT2(LISTB,ITAMPB,WORK(KGAMMA1),WORK(KGAMMA2),
     &                       LABELA,IRELAX,WORK(KEND1A),LEND1A)

          ELSE

            IFBTRAN(1) = IOPERA
            IFBTRAN(2) = IZETAV
            IFBTRAN(3) = ITAMPB
            IFBTRAN(4) = IFILE
            IFBTRAN(5) = IRELAX
 
            FILFBTA = 'CC_FBMAT'
            CALL CC_FBTA( IFBTRAN, 1, 0, LISTL, LISTB,
     &                    FILFBTA, IDUMMY,  DUMMY,
     &                       0,   WORK, LWORK)   
 
            IADR = IFBTRAN(4)
            LEN  = NT1AM(ISYRES) + NT2AM(ISYRES)
            LUFBTA = -1
            CALL WOPEN2(LUFBTA,FILFBTA,64,0)
            CALL GETWA2(LUFBTA,FILFBTA,WORK(KGAMMA1),IADR,LEN)
            CALL WCLOSE2(LUFBTA,FILFBTA,'DELETE')

          END IF

          IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'CCQR_FADRV> result of FB transform.:'
            CALL CC_PRP(WORK(KGAMMA1),WORK(KGAMMA2),ISYRES,1,1)
          END IF

          IF (OPTION(1:4).EQ.'FILE') THEN
C            -----------------------
C             add to vector on file:
C            -----------------------
             IF (CCS) THEN
                MODEL = 'CCS       '
                IOPT  = 1
             ELSE IF (CC2) THEN
                MODEL = 'CC2       '
                IOPT  = 3
             ELSE IF (CCSD) THEN
                MODEL = 'CCSD      '
                IOPT  = 3
             ELSE
                CALL QUIT('Unknown coupled cluster model in CCQR_FADRV')
             END IF
            
             CALL CC_WARSP(LISTC, IFILE, ISYRES, IOPT, MODEL,DUMMY,
     &                  WORK(KGAMMA1),WORK(KGAMMA2),WORK(KEND1),LWRK1)

           ELSE IF (OPTION(1:4).EQ.'DOTP') THEN
C            ---------------------------------
C             calculate required dot products:
C            ---------------------------------
             IVEC = 1
             DO WHILE (IFADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
             
               ITAMPC = IFADOTS(IVEC,ITRAN)
               ISYMTC = ILSTSYM(LISTC,ITAMPC)
            
               IF (ISYMTC.NE.ISYRES) THEN
                 CALL QUIT('symmetry mismatch in CCQR_FADRV.')
               END IF
             
               IOPT = 3
               CALL CC_RDRSP(LISTC,ITAMPC,ISYMTC,IOPT,MODEL,
     &                       WORK(KT1AMPC),WORK(KT2AMPC)   )
             
               FACON1  = DDOT(NT1AM(ISYRES),WORK(KT1AMPC),1,
     &                                      WORK(KGAMMA1),1)
               IF (CC2 .OR. CCSD) THEN
                FACON2 = DDOT(NT2AM(ISYRES),WORK(KT2AMPC),1,
     &                                      WORK(KGAMMA2),1)
               ELSE IF (CCS) THEN
                FACON2 = ZERO
               ELSE 
                CALL QUIT('ERROR> unknown coupled cluster model '//
     &                 'in CCQR_FADRV.')
               END IF
             
               FACON(IVEC,ITRAN) = FACON1 + FACON2
             
               IVEC = IVEC + 1
             END DO

          ELSE
             CALL QUIT('Illegal OPTION in CCQR_FADRV')
          END IF

       ELSE

C         CALL CCLR_FA( LABELA, ISYMTA, LISTB, ITAMPB,
C    &                  LISTL,  IZETAV, WORK,  LWORK   )

          IF      (OPTION(1:4).EQ.'FILE') THEN
             IOPTRES = 4
          ELSE IF (OPTION(1:4).EQ.'DOTP') THEN
             IOPTRES = 5
          ELSE
             CALL QUIT('Illegal OPTION in CCQR_FADRV')
          END IF

          ! the noddy code which goes via transformations is only
          ! used if we do transformations or if we do dot products
          ! and FADEN_NODDY is not set.
          NODDY_CCSDT = NODDY_FAMAT .AND. (.NOT.FADEN_NODDY)

          CALL CC_FAMAT( LABELA, ISYMTA,  ! inp: label/symmetry A
     &                   LISTB,  ITAMPB,  ! inp: B resp. amplit.
     &                   LISTL,  IZETAV,  ! inp: C resp. zeta vec.
     &                   IOPTRES,IFATRAN, ! output option,indeces
     &                   IFADOTS,FACON,   ! indeces/dotproducts
     &                   LISTC,  ITRAN,   ! list,index for dotprods
     &                   NFATRAN,MXVEC,   ! dimensions for dotprods
     &                   NODDY_CCSDT,     ! flag for noddy triples code
     &                   WORK,   LWORK   )! work space

       END IF

      END DO

*---------------------------------------------------------------------*
* add triples contributions for contractions (IOPTRES=5):
*---------------------------------------------------------------------*
      IF (CCSDT .AND. IOPTRES.EQ.5 .AND. 
     &    .NOT. ( NODDY_FAMAT .AND. .NOT. FADEN_NODDY ) ) THEN

        MXFADEN = NFATRAN*MXVEC
        MXINTER = 2*MXFADEN

        KIDXL_FADEN = 1
        KIDXB_FADEN = KIDXL_FADEN + MXFADEN
        KIDXC_FADEN = KIDXB_FADEN + MXFADEN
        KIADR_FADEN = KIDXC_FADEN + MXFADEN
        KEND1       = KIADR_FADEN + MXFADEN

        KIDXL_INTER = KEND1
        KIDXR_INTER = KIDXL_INTER + MXINTER  
        KLSTR_INTER = KIDXR_INTER + MXINTER  
        KIADR_INTER = KLSTR_INTER + MXINTER  
        KEND1       = KIADR_INTER + 3 * MXINTER  

        LWRK1 = LWORK - KEND1
 
        IF (LWRK1 .LT. 0) THEN
          CALL QUIT('Insufficient work space CCQR_FADRV')
        END IF 

        CALL CCSDT_FA_SETUP(IFATRAN,IFADOTS,NFATRAN,MXVEC,
     &                      WORK(KIDXL_FADEN),WORK(KIDXB_FADEN),
     &                      WORK(KIDXC_FADEN),NFADEN,MXFADEN,
     &                      WORK(KIDXL_INTER),WORK(KIDXR_INTER),
     &                      WORK(KLSTR_INTER),NINTER,MXINTER,
     &                      LISTL,LISTO,LISTB,LISTC)

        CALL CCSDT_FA_DEN(LISTL,LISTO,LISTB,LISTC,FADEN_NODDY,
     &                    NFATRAN,MXVEC,IFATRAN,IFADOTS,FACON,
     &                    FNDELD,FNCKJD,FNDKBC,FNTOC,
     &                    FN3VI,FNDKBC3,FN3FOPX,FN3FOP2X,
     &                    WORK(KIDXL_INTER),WORK(KIDXR_INTER),
     &                    WORK(KLSTR_INTER),NINTER,
     &                    WORK(KIDXL_FADEN),WORK(KIDXB_FADEN),
     &                    WORK(KIDXC_FADEN),NFADEN,
     &                    WORK(KIADR_INTER),WORK(KIADR_FADEN),
     &                    WORK(KEND1),LWRK1)

      END IF

*---------------------------------------------------------------------*
* print the results: 
*---------------------------------------------------------------------*
* F{A} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Results for F{O} matrix transformations:'
        IF (MXVEC.NE.0) THEN
         DO ITRAN = 1, NFATRAN
          WRITE(LUPRI,'(A,3I5)') MSGDBG, (IFATRAN(I,ITRAN),I=1,3)
          IVEC = 1
          DO WHILE (IFADOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
            WRITE(LUPRI,'(A,20X,I5,2X,E18.12)') MSGDBG,
     &          IFADOTS(IVEC,ITRAN), FACON(IVEC,ITRAN)
            IVEC = IVEC + 1
          END DO
         END DO
        END IF
        WRITE (LUPRI,*)
      END IF

      CALL QEXIT('CCQR_FADRV')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQR_FADRV                           *
*---------------------------------------------------------------------*
c /* deck ccqr_eadrv */
*=====================================================================*
      SUBROUTINE CCQR_EADRV(LISTL, LISTO, LISTB,
     &                      NETRAN, MXVEC,
     &                      IETRAN, IEDOTS, ECON,
     &                      WORK, LWORK, OPTION )
*---------------------------------------------------------------------*
*
*    Purpose: driver for a list of ETA{A} vector calculations
*
*    Written by Christof Haettig, December 1996.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <dummy.h>
#include <ccroper.h>

* local parameters:
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCQR_EADRV> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LWORK, NETRAN, MXVEC

#if defined (SYS_CRAY)
      REAL ECON(MXVEC,NETRAN)
      REAL WORK(LWORK)
      REAL ZERO
      REAL ECON1, ECON2
      REAL DDOT
#else
      DOUBLE PRECISION ECON(MXVEC,NETRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION ECON1, ECON2
      DOUBLE PRECISION DDOT
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER OPTION*(*)
      CHARACTER LISTL*(*), LISTO*(*), LISTB*(*)
      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABELA
      LOGICAL LTWOEL
      INTEGER IETRAN(3,NETRAN)
      INTEGER IEDOTS(MXVEC,NETRAN)
      INTEGER ISYRES, ISYMTB, ISYCTR, ISYMTA
      INTEGER IVEC, ITRAN, ITAMPB, IZETAV, IOPT, IOPERA, IFILE
      INTEGER KETA1, KETA2, KT1AMPB, KT2AMPB, KEND1,LEND1, KEND2,LEND2


* external functions:
      INTEGER ILSTSYM

      CALL QENTER('CCQR_EADRV')

* print input, if locdbg set:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of ETA{O} vector calculations:'
        DO ITRAN = 1, NETRAN
          WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &     (IETRAN(I,ITRAN),I=1,2),(IEDOTS(I,ITRAN),I=1,MXVEC)
        END DO
      END IF


*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      IF (MXVEC*NETRAN.NE.0) CALL DZERO(ECON,MXVEC*NETRAN)

*---------------------------------------------------------------------*
* start loop over all requested ETA vector calculations
*---------------------------------------------------------------------*
      DO ITRAN = 1, NETRAN

        IZETAV = IETRAN(1,ITRAN)
        IOPERA = IETRAN(2,ITRAN)
        IFILE  = IETRAN(3,ITRAN)

        LTWOEL = LPDBSOP(IOPERA)
        LABELA = LBLOPR(IOPERA)
        IF (LTWOEL) THEN  
           WRITE (LUPRI,*) 'List inconsistency in CCQR_EADRV:'
           WRITE (LUPRI,*) 'LISTO, IOPERA :',LISTO,IOPERA
           WRITE (LUPRI,*) 'LABELA, LTWOEL:',LABELA,LTWOEL
           WRITE (LUPRI,*) 'CCQR_EADRV no adapted for '//
     &          'derivative integrals.'
        END IF 


        ISYCTR = ILSTSYM(LISTL,IZETAV)
        ISYMTA = ILSTSYM(LISTO,IOPERA)

        ISYRES = MULD2H(ISYMTA,ISYCTR)

        KETA1   = 1
        KETA2   = KETA1   + NT1AM(ISYRES)
        KEND1   = KETA2   + NT2AM(ISYRES)
        LEND1   = LWORK - KEND1

        KT1AMPB = KEND1
        KT2AMPB = KT1AMPB + NT1AM(ISYRES)
        KEND2   = KT2AMPB + NT2AM(ISYRES)
        LEND2   = LWORK - KEND2


        IF (LEND2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CCQR_EADRV.')
        END IF

*---------------------------------------------------------------------*
* calculate ETA vector:
*      ETA{A} = <Lambda| [A,tau] |CC>
*---------------------------------------------------------------------*
      CALL CC_ETAC(ISYMTA, LABELA, WORK(KETA1),
     &             LISTL,IZETAV,0,DUMMY,WORK(KEND1),LEND1)

      IF (OPTION(1:4).EQ.'FILE') THEN
*---------------------------------------------------------------------*
* add to vector on file:
*---------------------------------------------------------------------*
        IF (CCS) THEN
           MODEL = 'CCS       '
           IOPT  = 1
        ELSE IF (CC2) THEN
           MODEL = 'CC2       '
           IOPT  = 3
        ELSE IF (CCSD) THEN
           MODEL = 'CCSD      '
           IOPT  = 3
        ELSE
           CALL QUIT('Unknown coupled cluster model in CC_EADRV.')
        END IF

        CALL CC_WARSP(LISTB, IFILE, ISYRES, IOPT, MODEL,DUMMY,
     &                WORK(KETA1),WORK(KETA2),WORK(KEND1),LEND1)

      ELSE IF (OPTION(1:4).EQ.'DOTP') THEN
*---------------------------------------------------------------------*
* calculate required dot products:
*---------------------------------------------------------------------*
        IVEC = 1
        DO WHILE (IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)

          ITAMPB = IEDOTS(IVEC,ITRAN)
          ISYMTB = ILSTSYM(LISTB,ITAMPB)
       
          IF (ISYMTB.NE.ISYRES) THEN
            WRITE (LUPRI,*) LISTL, IZETAV, ISYCTR
            WRITE (LUPRI,*) LISTO, IOPERA, ISYMTA
            WRITE (LUPRI,*) LISTB, ITAMPB, ISYMTB
            CALL QUIT('symmetry mismatch in CCQR_EADRV.')
          END IF

          IOPT = 3
          CALL CC_RDRSP(LISTB,ITAMPB,ISYMTB,IOPT,MODEL,
     &                  WORK(KT1AMPB),WORK(KT2AMPB)   )

          ECON1  = DDOT(NT1AM(ISYRES),WORK(KT1AMPB),1,WORK(KETA1),1)
          IF (CC2 .OR. CCSD) THEN
           ECON2 = DDOT(NT2AM(ISYRES),WORK(KT2AMPB),1,WORK(KETA2),1)
          ELSE IF (CCS) THEN
           ECON2 = ZERO
          ELSE
           CALL QUIT('ERROR> unknown CC model in CCQR_EADRV.')
          END IF

          ECON(IVEC,ITRAN) = ECON1 + ECON2

          IVEC = IVEC + 1
        END DO
*---------------------------------------------------------------------*
      ELSE 
        CALL QUIT('Illegal OPTION in CCQR_EADRV.')
      END IF

      END DO

*---------------------------------------------------------------------*
* print the results: 
*---------------------------------------------------------------------*
* ETA vector results:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Results for ETA vector calculations:'
        IF (MXVEC.NE.0) THEN
         DO ITRAN = 1, NETRAN
          WRITE (LUPRI,*)
          IVEC = 1
          DO WHILE (IEDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
           WRITE(LUPRI,'(A,3I5,2X,E18.12)') MSGDBG, 
     &      (IETRAN(I,ITRAN),I=1,2),IEDOTS(IVEC,ITRAN),ECON(IVEC,ITRAN)
           IVEC = IVEC + 1
          END DO
         END DO
        END IF
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

      CALL QEXIT('CCQR_EADRV')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCQR_EADRV                           *
*---------------------------------------------------------------------*
c /* deck CC_DOTDRV */
*=====================================================================*
      SUBROUTINE CC_DOTDRV(LISTL, LISTR, 
     &                     NDTRAN, MXVEC,
     &                     IDTRAN, IDDOTS, DCON,
     &                     WORK, LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: driver for a list of vector dot product calculations
*
*             some special LIST types are used for this routine:
*
*               LISTR = X1B  --  right hand side for L1 equations
*               LISTR = X2B  --  right hand side for L2 equations
*             
*
*    Written by Christof Haettig, April 1997.
*    adapted for CC-R12, Christian Neiss, June 2005
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccroper.h>
#include <ccl1rsp.h>
#include <cclc1rsp.h>
#include <dummy.h>
#include <ccnoddy.h>

* local parameters:
      CHARACTER*(19) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_DOTDRV> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LWORK, NDTRAN, MXVEC

#if defined (SYS_CRAY)
      REAL DCON(MXVEC,NDTRAN)
      REAL WORK(LWORK)
      REAL ZERO
      REAL DCON1, DCON2, DCONR12
      REAL DDOT
#else
      DOUBLE PRECISION DCON(MXVEC,NDTRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO
      DOUBLE PRECISION DCON1, DCON2, DCONR12
      DOUBLE PRECISION DDOT
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER LISTL*(*), LISTR*(*)
      CHARACTER*(8) LABEL
      CHARACTER*(10) MODEL
      INTEGER IDTRAN(NDTRAN)
      INTEGER IDDOTS(MXVEC,NDTRAN)

      INTEGER ITRAN, IVEC, IFILEL, IFILER, ISYML, ISYMR, IOPT, ICAU
      INTEGER KVECL1, KVECL2, KVECR1, KVECR2, KVECLR12, KVECRR12
      INTEGER KEND1, KEND2, LEND1, LEND2, NTAMP, IOPTR12

* external functions:
      INTEGER ILSTSYM
      INTEGER ILRCAMP


      CALL QENTER('CC_DOTDRV')
*---------------------------------------------------------------------*
* initializations & checks:
*---------------------------------------------------------------------*
      IF (MXVEC*NDTRAN.NE.0) THEN
        CALL DZERO(DCON,MXVEC*NDTRAN)
      ELSE
        CALL QEXIT('CC_DOTDRV')
        RETURN
      END IF

* print input, if locdbg set:
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) MSGDBG, 'List of vector dot '//
     &        'product calculations:'
        WRITE(LUPRI,*) MSGDBG, 'LISTL, LISTR:',LISTL,LISTR
        DO ITRAN = 1, NDTRAN
          WRITE(LUPRI,'(A,I5,5X,(25I3,20X))') MSGDBG,
     &     IDTRAN(ITRAN),(IDDOTS(I,ITRAN),I=1,MXVEC)
        END DO
      END IF

*---------------------------------------------------------------------*
* start loop over left vectors:
*---------------------------------------------------------------------*
      DO ITRAN = 1, NDTRAN
        IFILEL = IDTRAN(ITRAN)

        IF (LISTL.EQ.'X2B') THEN
          ISYML  = ILSTSYM('L2 ',IFILEL)
        ELSE
          ISYML  = ILSTSYM(LISTL,IFILEL)
        END IF

        NTAMP = NT1AM(ISYML) + NT2AM(ISYML)
        IF (CCR12) NTAMP = NTAMP + NTR12AM(ISYML)

        KVECL1 = 1
        KVECL2 = KVECL1 + NT1AM(ISYML)
        KEND1  = KVECL2 + NT2AM(ISYML)
        IF (CCR12) THEN
          KVECLR12 = KVECL2 + NT2AM(ISYML)
          KEND1    = KVECLR12 + NTR12AM(ISYML)
        END IF
        LEND1  = LWORK - KEND1

        KVECR1 = KEND1
        KVECR2 = KVECR1 + NT1AM(ISYML)
        KEND2  = KVECR2 + NT2AM(ISYML)
        IF (CCR12) THEN
          KVECRR12 = KVECR2 + NT2AM(ISYML)
          KEND2    = KVECRR12 + NTR12AM(ISYML)
        END IF
        LEND2  = LWORK - KEND2

        IF (CCS) THEN
          KVECR1 = KVECL2
          KEND2  = KVECR1 + NT1AM(ISYML)
          LEND2  = LWORK - KEND2
          NTAMP  = NT1AM(ISYML)
        END IF

        IF (LEND2 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_DOTDRV.')
        END IF

* read left vector:
        IOPT = 3
        IOPTR12 = 32
        IF (LISTL(1:2).EQ.'O3') THEN
          CALL CC_RDRSP(LISTL,IFILEL,ISYML,IOPT,MODEL,
     &                  WORK(KVECL1),WORK(KVECL2)     )
          CALL CCLR_DIASCL(WORK(KVECL2),0.5d0,ISYML)
          IF (CCR12) THEN
            CALL CC_RDRSP(LISTL,IFILEL,ISYML,IOPTR12,MODEL,
     &                    DUMMY,WORK(KVECLR12))
          END IF
        ELSE IF (LISTL(1:3).EQ.'CO2'.OR.LISTL(1:3).EQ.'O2 ') THEN
          CALL CC_RDRSP(LISTL,IFILEL,ISYML,IOPT,MODEL,
     &                  WORK(KVECL1),WORK(KVECL2)     )
          CALL CCLR_DIASCL(WORK(KVECL2),0.5d0,ISYML)
          IF (CCR12) THEN
            CALL CC_RDRSP(LISTL,IFILEL,ISYML,IOPTR12,MODEL,
     &                    DUMMY,WORK(KVECLR12))
          END IF
        ELSE IF (LISTL(1:3).EQ.'X1B') THEN
          CALL CC_LGD(WORK(KVECL1),LRZLBL(IFILEL),ISYLRZ(IFILEL),
     &         FRQLRZ(IFILEL),LORXLRZ(IFILEL),IOPT,WORK(KEND2),LEND2)
        ELSE IF (LISTL(1:3).EQ.'X2B') THEN
          CALL CC_L2GD(WORK(KVECL1),IFILEL,WORK(KEND2),LEND2)
        ELSE
          CALL CC_RDRSP(LISTL,IFILEL,ISYML,IOPT,MODEL,
     &                  WORK(KVECL1),WORK(KVECL2)     )
          IF (CCR12) THEN
            CALL CC_RDRSP(LISTL,IFILEL,ISYML,IOPTR12,MODEL,
     &                    DUMMY,WORK(KVECLR12))
          END IF
        END IF

*---------------------------------------------------------------------*
* calculate singles and doubles contributions of the required dot 
* products (the vectors are mainly read from file)
*---------------------------------------------------------------------*
        IVEC = 1
        DO WHILE (IDDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)

          IFILER = IDDOTS(IVEC,ITRAN)
          IF (LISTR.EQ.'X1B') THEN
            ISYMR  = ILSTSYM('L1 ',IFILER)
          ELSE
            ISYMR  = ILSTSYM(LISTR,IFILER)
          END IF
       
          IF (ISYMR.NE.ISYML) THEN
            WRITE (LUPRI,*) LISTL, IFILEL, ISYML
            WRITE (LUPRI,*) LISTR, IFILER, ISYMR
            CALL QUIT('symmetry mismatch in CC_DOTDRV.')
          END IF

          IOPT = 3
          IOPTR12 = 32

          IF (LISTR(1:2).EQ.'O2') THEN
            CALL CC_RDRSP(LISTR,IFILER,ISYMR,IOPT,MODEL,
     &                    WORK(KVECR1),WORK(KVECR2)   )
            CALL CCLR_DIASCL(WORK(KVECR2),0.5d0,ISYMR)
            IF (CCR12) THEN
              CALL CC_RDRSP(LISTR,IFILER,ISYMR,IOPTR12,MODEL,
     &                      DUMMY,WORK(KVECRR12))
            END IF
          ELSE IF (LISTR(1:3).EQ.'CO2') THEN
            CALL CC_RDRSP(LISTR,IFILER,ISYMR,IOPT,MODEL,
     &                    WORK(KVECR1),WORK(KVECR2)   )
            CALL CCLR_DIASCL(WORK(KVECR2),0.5d0,ISYMR)
            IF (CCR12) THEN
              CALL CC_RDRSP(LISTR,IFILER,ISYMR,IOPTR12,MODEL,
     &                      DUMMY,WORK(KVECRR12))
            END IF
          ELSE IF (LISTR(1:3).EQ.'X1B') THEN
            CALL CC_LGD(WORK(KVECR1),LRZLBL(IFILER),ISYLRZ(IFILER),
     &           FRQLRZ(IFILER),LORXLRZ(IFILER),IOPT,WORK(KEND2),LEND2)
          ELSE IF (LISTR(1:2).EQ.'XC') THEN
            LABEL = LBLLC1(IFILER)
            ICAU  = ILC1CAU(IFILER)
            IF (ICAU.EQ.0) THEN ! return ( F x R1(0) + eta )
              CALL CC_LGD(WORK(KVECR1),LABEL,ISYMR,
     &                    0.0d0,.FALSE.,IOPT,WORK(KEND2),LEND2)
            ELSE ! return F x RC
              IFILER = ILRCAMP(LABEL,ICAU,ISYMR)
              CALL CC_RDRSP('FC',IFILER,ISYMR,IOPT,MODEL,
     &                      WORK(KVECR1),WORK(KVECR2)   )
              IF (CCR12) CALL QUIT('Calculations including Cauchy '//
     &                             'vectors not adapted for CC-R12')
            END IF
          ELSE
            CALL CC_RDRSP(LISTR,IFILER,ISYMR,IOPT,MODEL,
     &                    WORK(KVECR1),WORK(KVECR2)   )
            IF (CCR12) THEN
              CALL CC_RDRSP(LISTR,IFILER,ISYMR,IOPTR12,MODEL,
     &                      DUMMY,WORK(KVECRR12))
            END IF
          END IF

          DCON1  = DDOT(NT1AM(ISYML),WORK(KVECL1),1,WORK(KVECR1),1)
          IF (LOCDBG) THEN
            WRITE(LUPRI,*) 'NORM^2(L1):',
     &        DDOT(NT1AM(ISYML),WORK(KVECL1),1,WORK(KVECL1),1)
            WRITE(LUPRI,*) 'NORM^2(R1):',
     &        DDOT(NT1AM(ISYMR),WORK(KVECR1),1,WORK(KVECR1),1)
            WRITE(LUPRI,*) '(L1|R1):',
     &        DDOT(NT1AM(ISYMR),WORK(KVECL1),1,WORK(KVECR1),1)
          END IF
          IF (CC2 .OR. CCSD .OR. CC3) THEN
           DCON2 = DDOT(NT2AM(ISYML),WORK(KVECL2),1,WORK(KVECR2),1)
          ELSE IF (CCS) THEN
           DCON2 = ZERO
          ELSE
           CALL QUIT('ERROR> unknown CC model in CC_DOTDRV.')
          END IF
          IF (CCR12) THEN
           DCONR12 = DDOT(NTR12AM(ISYML),WORK(KVECLR12),1,
     &                    WORK(KVECRR12),1)
           IF (LOCDBG) THEN
             WRITE(LUPRI,*) 'NORM^2(LR12):',
     &         DDOT(NTR12AM(ISYML),WORK(KVECLR12),1,WORK(KVECLR12),1)
             WRITE(LUPRI,*) 'NORM^2(RR12):',
     &         DDOT(NTR12AM(ISYMR),WORK(KVECRR12),1,WORK(KVECRR12),1)
             WRITE(LUPRI,*) '(LR12|RR12):',
     &         DDOT(NTR12AM(ISYMR),WORK(KVECLR12),1,WORK(KVECRR12),1)
           END IF
          ELSE
           DCONR12 = ZERO
          END IF

          DCON(IVEC,ITRAN) = DCON1 + DCON2 + DCONR12

          IVEC = IVEC + 1
        END DO

      END DO

*---------------------------------------------------------------------*
* calculate the triples contributions by explicitly setting up the 
* triples part of the vectors:
*---------------------------------------------------------------------*
      IF (CCSDT) THEN
        IF (NODDY_OVLP) THEN  
          CALL CCOVLP_NODDY(LISTL,LISTR,2,DUMMY,1,1,1,
     &                      NDTRAN,MXVEC,IDTRAN,IDDOTS,DCON,
     &                      WORK,LWORK) 
        ELSE
          CALL CC3_LR_OVLP(LISTL,LISTR,2,DUMMY,1,1,1,
     &                     NDTRAN,MXVEC,IDTRAN,IDDOTS,DCON,
     &                     WORK,LWORK)  
        END IF
      END IF

*---------------------------------------------------------------------*
* print the results: 
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE(LUPRI,*) MSGDBG,
     &        'Results for vector dot product calculations:'
        IF (MXVEC.NE.0) THEN
         DO ITRAN = 1, NDTRAN
          WRITE (LUPRI,*)
          IVEC = 1
          DO WHILE (IDDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
           WRITE(LUPRI,'(A,2I5,2X,E18.12)') MSGDBG, 
     &      IDTRAN(ITRAN),IDDOTS(IVEC,ITRAN),DCON(IVEC,ITRAN)
           IVEC = IVEC + 1
          END DO
         END DO
        END IF
        WRITE (LUPRI,*)
        CALL FLSHFO(LUPRI)
      END IF

      CALL QEXIT('CC_DOTDRV')
      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_DOTDRV                            *
*---------------------------------------------------------------------*
