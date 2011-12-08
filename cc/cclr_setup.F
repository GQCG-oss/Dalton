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
c /* deck cclr_setup */
*=====================================================================*
      SUBROUTINE CCLR_SETUP(MXTRAN,  MXVEC,
     &                     IFTRAN,  IFDOTS,  FCONS,  NFTRAN,
     &                     IJTRAN,  IJDOTS,  CJCON,  NJTRAN,
     &                     IXITRAN, IXIDOTS, XICONS, NXITRAN,
     &                     IRTRAN,  IRDOTS,  RCONS,  NRTRAN,
     &                     IXETRAN,IXDOTS,IEDOTS,XCONS,ECONS,NXETRAN,
     &                     RESULT,  MXSOP,   LADD,   WORK, LWORK )
*---------------------------------------------------------------------*
*
*    Purpose: set up for CC linear response section:
*         - list of F matrix transformations with Cauchy vectors
*         - list of XKSI and ETA vector calculations 
*         - list of X intermediate contributions
*         - list of second-order reortho./relax. contributions
*
*     Written by Christof Haettig, may 1999 based on CCCM_SETUP 
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cclrinf.h>
#include <ccroper.h>
#include <ccr1rsp.h>
#include <ccsdinp.h>
#include <ccexpfck.h>
#include <cclists.h>
#include <ccsections.h>
#include <ccslvinf.h>

* local parameters:
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCLR_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LADD
      INTEGER MXVEC, MXTRAN, MXSOP

      INTEGER IFTRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)
      INTEGER IJTRAN(MXDIM_JTRAN,MXTRAN)
      INTEGER IJDOTS(MXVEC,MXTRAN)
      INTEGER IXITRAN(1,MXTRAN)
      INTEGER IXIDOTS(MXVEC,MXTRAN)
      INTEGER IRTRAN(1,MXTRAN)
      INTEGER IRDOTS(MXVEC,MXTRAN)
      INTEGER IXETRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IXDOTS(MXVEC,MXTRAN), IEDOTS(MXVEC,MXTRAN)

      INTEGER NFTRAN, NXITRAN, NRTRAN, NXETRAN, LWORK
      INTEGER NJTRAN

#if defined (SYS_CRAY)
      REAL RESULT(MXSOP)
      REAL FCONS(MXVEC,MXTRAN)
      REAL CJCON(MXVEC,MXTRAN)
      REAL XICONS(MXVEC,MXTRAN)
      REAL RCONS(MXVEC,MXTRAN)
      REAL XCONS(MXVEC,MXTRAN), ECONS(MXVEC,MXTRAN)
      REAL WORK(LWORK)
      REAL ZERO, SIGN
      REAL WSTAT, WNUCL, WREO, WXE1, WXE2, WXI1, WXI2, WF, WJ
#else
      DOUBLE PRECISION RESULT(MXSOP)
      DOUBLE PRECISION FCONS(MXVEC,MXTRAN)
      DOUBLE PRECISION CJCON(MXVEC,MXTRAN)
      DOUBLE PRECISION XICONS(MXVEC,MXTRAN)
      DOUBLE PRECISION RCONS(MXVEC,MXTRAN)
      DOUBLE PRECISION XCONS(MXVEC,MXTRAN), ECONS(MXVEC,MXTRAN)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, SIGN
      DOUBLE PRECISION WSTAT, WNUCL, WREO, WXE1, WXE2, WXI1, WXI2, WF,WJ
#endif
      PARAMETER (ZERO = 0.0D0)

      LOGICAL LORXA, LORXB, LPDBSA, LPDBSB
      INTEGER ISYMA,  ISYMB, ITRAN, IVEC, ISYML, IDUM, I, N, IFREQ
      INTEGER IR1VECA,IR1VECB,IOPERA,IOPERB,IEATA1A,IEATA1B, ISGNSOP
      INTEGER IL1VECB, INUM, IOPER, NBSOP, ISIGN, ISYOP, IKAPA, IKAPB
      INTEGER MEAVEC,MFVEC,MXAVEC,MXIVEC,MXRVEC,MXEVEC, IDX, IEXPV
      INTEGER MJVEC, IL1VECA

      CHARACTER LABELA*(8), LABELB*(8), LABSOP*(8)

* external functions:
      INTEGER IR1TAMP
      INTEGER IR1KAPPA
      INTEGER IL1ZETA
      INTEGER IETA1
      INTEGER IEXPECT
#if defined (SYS_CRAY)
      REAL CC_NUCCON
#else
      DOUBLE PRECISION CC_NUCCON
#endif

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      DO ITRAN = 1, MXTRAN
       DO I = 1, MXDIM_XEVEC
        IXETRAN(I,ITRAN)  = 0
       END DO
       DO I = 1, MXDIM_FTRAN
        IFTRAN(I,ITRAN)  = 0
       END DO
       DO I = 1, MXDIM_JTRAN
        IJTRAN(I,ITRAN)  = 0
       END DO
       DO I = 1, 1
        IXITRAN(I,ITRAN) = 0
        IRTRAN(I,ITRAN)  = 0
       END DO

       DO IVEC  = 1, MXVEC
        IFDOTS(IVEC,ITRAN)  = 0
        IJDOTS(IVEC,ITRAN)  = 0
        IXIDOTS(IVEC,ITRAN) = 0
        IRDOTS(IVEC,ITRAN)  = 0
        IXDOTS(IVEC,ITRAN)  = 0
        IEDOTS(IVEC,ITRAN)  = 0
       END DO
      END DO

      NFTRAN  = 0
      NJTRAN  = 0
      NXITRAN = 0
      NRTRAN  = 0
      NXETRAN = 0
 
      MFVEC   = 0 
      MJVEC   = 0
      MXIVEC  = 0
      MXRVEC  = 0
      MXEVEC  = 0

      NBSOP   = 0

*---------------------------------------------------------------------*
* start loop over all requested second-order properties:
*---------------------------------------------------------------------*
      DO IOPER = 1, NLROP
        IOPERA = IALROP(IOPER)
        IOPERB = IBLROP(IOPER)
        LORXA  = LALORX(IOPER)
        LORXB  = LBLORX(IOPER)
        ISYMA  = ISYOPR(IOPERA)
        ISYMB  = ISYOPR(IOPERB)
        LABELA = LBLOPR(IOPERA)
        LABELB = LBLOPR(IOPERB)
        LPDBSA = LPDBSOP(IOPERA)
        LPDBSB = LPDBSOP(IOPERB)


      IF (ISYMA.EQ.ISYMB) THEN

        DO IFREQ = 1, NBLRFR
           NBSOP = NBSOP + 1

           IF (NBSOP.GT.MXSOP) THEN
              CALL QUIT('NBSOP out of range in CCLR_SETUP.')
           END IF

        DO ISIGN = 1, -1, -2
           SIGN = DBLE(ISIGN)

*---------------------------------------------------------------------*
*          in all cases we need Eta{A} x R1^B 
*---------------------------------------------------------------------*
           IR1VECB = IR1TAMP(LABELB,LORXB,SIGN*BLRFR(IFREQ),IDUM)
C          IEATA1A =   IETA1(LABELA,LORXA,SIGN*ALRFR(IFREQ),IDUM)
           IF (LORXA) THEN
             IKAPA = IR1KAPPA(LABELA,SIGN*ALRFR(IFREQ),IDUM)
           ELSE
             IKAPA = 0
           END IF
     
           CALL CC_SETXE('Eta',IXETRAN,IEDOTS,MXTRAN,MXVEC,
     &                   0,IOPERA,IKAPA,0,0,0,IR1VECB,ITRAN,IVEC)
           NXETRAN = MAX(NXETRAN,ITRAN)
           MXEVEC  = MAX(MXEVEC, IVEC)
           WXE1    = ECONS(IVEC,ITRAN)

           IF (.NOT. ASYMSD) THEN
*---------------------------------------------------------------------*
*             symmetric approach: add F * R1^A * R1^B + Eta{B} x R1^A
*---------------------------------------------------------------------*
              IR1VECA = IR1TAMP(LABELA,LORXA,SIGN*ALRFR(IFREQ),IDUM)
C             IEATA1B =   IETA1(LABELB,LORXB,SIGN*BLRFR(IFREQ),IDUM)
              IF (LORXB) THEN
                IKAPB = IR1KAPPA(LABELB,SIGN*BLRFR(IFREQ),IDUM)
              ELSE
                IKAPB = 0
              END IF

              IF (CIS) THEN
                 WF     = ZERO
              ELSE
                 CALL CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                          0,IR1VECA,IR1VECB,ITRAN,IVEC)
                 NFTRAN = MAX(NFTRAN,ITRAN)
                 MFVEC  = MAX(MFVEC, IVEC)
                 WF     = FCONS(IVEC,ITRAN)
C
                 IF (CCSLV) THEN
                 IL1VECA = IL1ZETA(LABELA,LORXA,SIGN*ALRFR(IFREQ),IDUM)
                 IL1VECB = IL1ZETA(LABELB,LORXB,SIGN*BLRFR(IFREQ),IDUM)
C                Here, we substitute amplitudes <-> multipilers 
                   CALL CCQR_SETF(IJTRAN,IJDOTS,MXTRAN,MXVEC,
     &                            0,IL1VECA,IL1VECB,ITRAN,IVEC)
                   NJTRAN = MAX(NJTRAN,ITRAN)
                   MJVEC  = MAX(MJVEC, IVEC)
                   WJ     = CJCON(IVEC,ITRAN)
                 ELSE
                   WJ = ZERO
                 END IF
              END IF

              CALL CC_SETXE('Eta',IXETRAN,IEDOTS,MXTRAN,MXVEC,
     &                      0,IOPERB,IKAPB,0,0,0,IR1VECA,ITRAN,IVEC)
              NXETRAN = MAX(NXETRAN,ITRAN)
              MXEVEC  = MAX(MXEVEC, IVEC)
              WXE2    = ECONS(IVEC,ITRAN)

           ELSE
*---------------------------------------------------------------------*
*             asymmetric approach: add L1^B x Xksi{A}
*---------------------------------------------------------------------*
              WF      = ZERO
              WJ      = ZERO

              IL1VECB = IL1ZETA(LABELB,LORXB,SIGN*BLRFR(IFREQ),IDUM)

              CALL CC_SETXE('Xi ',IXETRAN,IXDOTS,MXTRAN,MXVEC,
     &                      0,IOPERA,IKAPA,0,0,0,IL1VECB,ITRAN,IVEC)
              NXETRAN = MAX(NXETRAN,ITRAN)
              MXEVEC  = MAX(MXEVEC, IVEC)
              WXE2    = XCONS(IVEC,ITRAN)

           END IF

*---------------------------------------------------------------------*
*          for orbital relaxed second-order properties or if we have
*          perturbation-dependent basis sets involved add the contrib.
*          from the first-order effective Fock matrix times the 
*          Q matrix (kappa + R) :
*---------------------------------------------------------------------*
           WXI1 = ZERO
           WXI2 = ZERO

           IF (LORXB.OR.LPDBSB) THEN

              IKAPA = IR1KAPPA(LABELA,SIGN*ALRFR(IFREQ),IDUM)
              IKAPB = IR1KAPPA(LABELB,SIGN*BLRFR(IFREQ),IDUM)

              CALL CC_SETDOT(IXITRAN,IXIDOTS,MXTRAN,MXVEC,
     &                       IKAPA,IKAPB,ITRAN,IVEC)
              NXITRAN = MAX(NXITRAN,ITRAN)
              MXIVEC  = MAX(MXIVEC, IVEC)
              WXI1    = WXI1 + XICONS(IVEC,ITRAN)


           END IF


           IF (LORXA.OR.LPDBSA) THEN

              IKAPA = IR1KAPPA(LABELA,SIGN*ALRFR(IFREQ),IDUM)
              IKAPB = IR1KAPPA(LABELB,SIGN*BLRFR(IFREQ),IDUM)

              CALL CC_SETDOT(IXITRAN,IXIDOTS,MXTRAN,MXVEC,
     &                       IKAPB,IKAPA,ITRAN,IVEC)
              NXITRAN = MAX(NXITRAN,ITRAN)
              MXIVEC  = MAX(MXIVEC, IVEC)
              WXI2    = WXI2 + XICONS(IVEC,ITRAN)

           END IF

*---------------------------------------------------------------------*
*          for derivatives we might need to include coupling between 
*          relaxation and reorthonormalization:
*---------------------------------------------------------------------*
           IF ( (LPDBSA .AND. LORXB) .OR. (LPDBSB .AND. LORXA) ) THEN

              IKAPA = IR1KAPPA(LABELA,SIGN*ALRFR(IFREQ),IDUM)
              IKAPB = IR1KAPPA(LABELB,SIGN*BLRFR(IFREQ),IDUM)

              CALL CC_SETDOT(IRTRAN,IRDOTS,MXTRAN,MXVEC,
     &                       IKAPB,IKAPA,ITRAN,IVEC)
              NRTRAN = MAX(NRTRAN,ITRAN)
              MXRVEC = MAX(MXRVEC,IVEC)
              WREO   = RCONS(IVEC,ITRAN)

           ELSE
              WREO = ZERO
           END IF

*---------------------------------------------------------------------*
*          get "static" and nuclear contribution:
*---------------------------------------------------------------------*
           IF (LPDBSA .OR. LPDBSB) THEN

              CALL CC_FIND_SO_OP(LABELA,LABELB,LABSOP,ISYOP,ISGNSOP,
     &                           INUM,WORK,LWORK)

              IF (INUM.LT.0) CALL QUIT('Operator error in CCLR_SETUP.')

              IEXPV = IEXPECT(LABSOP,ISYOP,1)
              WSTAT = EXPVALUE(1,IEXPV) + EXPVALUE(2,IEXPV)

              WNUCL = CC_NUCCON(LABSOP,ISYOP)

           ELSE
              WSTAT = ZERO
              WNUCL = ZERO
           END IF

*---------------------------------------------------------------------*
*          add contributions together:
*---------------------------------------------------------------------*
           IF (LADD) THEN

              IDX = NBLRFR*(IOPER-1) + IFREQ
              IF (ISIGN.EQ.-1) IDX = IDX + NBLRFR*NLROP

              RESULT(IDX) = WXE1+WXE2+WF-WJ+WXI1+WXI2+WREO-WSTAT-WNUCL

              IF (LOCDBG) THEN
                 WRITE (LUPRI,*) 'LABELA,LABELB:',LABELA,LABELB
                 WRITE (LUPRI,*) 'IFREQ:',IFREQ
                 WRITE (LUPRI,*) 'IDX = ',IDX
                 WRITE (LUPRI,*) 'WSTAT,WNUCL:',WSTAT,WNUCL
                 WRITE (LUPRI,*) 'WXI1,WXI2,WREO:', WXI1, WXI2, WREO
                 WRITE (LUPRI,*) 'WXE1, WXE2, WF, WJ :', 
     &                            WXE1, WXE2, WF, WJ
                 WRITE (LUPRI,*) 'RESULT:',RESULT(IDX)
              END IF

           END IF

*---------------------------------------------------------------------*
*       end loop over second-order properties
*---------------------------------------------------------------------*
        END DO
        END DO

      END IF
      END DO

      IF      (MFVEC.GT.MXVEC) THEN
         CALL QUIT('MFVEC has been out of bounds in CCLR_SETUP.')
      ELSE IF (MXEVEC.GT.MXVEC) THEN
         CALL QUIT('MXEVEC has been out of bounds in CCLR_SETUP.')
      ELSE IF (MXIVEC.GT.MXVEC) THEN
         CALL QUIT('MXIVEC has been out of bounds in CCLR_SETUP.')
      ELSE IF (MXRVEC.GT.MXVEC) THEN
         CALL QUIT('MXRVEC has been out of bounds in CCLR_SETUP.')
      ELSE IF (MJVEC.GT.MXVEC) THEN
          CALL QUIT('MJVEC has been out of bounds in CCLR_SETUP.')
      END IF

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      IF ((.NOT.LADD) .OR. LOCDBG) THEN
       WRITE(LUPRI,'(/,/3X,A,I3,A)') 'For the requested',NBSOP,
     &      ' second-order properties'
       WRITE(LUPRI,'((8X,A,I3,A))') 
     & ' - ',NFTRAN,' F matrix transformations with R1 vectors',
     & ' - ',NJTRAN,' J matrix transformations with L1 vectors',
     & ' - ',NXETRAN,' ETA and XKSI vector calculations ',
     & ' - ',NXITRAN,' X intermediate calculations ',
     & ' - ',NRTRAN, ' 2. order reortho./relax. contributions' 
       WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'
      END IF

      IF (LOCDBG) THEN

         ! F matrix transformations:
         WRITE(LUPRI,*)'List of F matrix transformations:'
         DO ITRAN = 1, NFTRAN
           WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &      (IFTRAN(I,ITRAN),I=1,2),(IFDOTS(I,ITRAN),I=1,MFVEC)
         END DO
         WRITE(LUPRI,*)

         ! J matrix transformations:
         WRITE(LUPRI,*)'List of J matrix transformations:'
         DO ITRAN = 1, NJTRAN
           WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &      (IJTRAN(I,ITRAN),I=1,2),(IJDOTS(I,ITRAN),I=1,MJVEC)
         END DO
         WRITE(LUPRI,*)

         ! Xi{O} and ETA{O} vector calculations:
         WRITE(LUPRI,*) 'List of Xi{O} and ETA{O} vector calculations:'
         DO ITRAN = 1, NXETRAN
           WRITE(LUPRI,'(A,5I5,5X,(25I3,20X))') MSGDBG,
     &      (IXETRAN(I,ITRAN),I=1,5),(IXDOTS(I,ITRAN),I=1,MXEVEC)
           WRITE(LUPRI,'(A,25X,5X,(25I3,20X))') MSGDBG,
     &                               (IEDOTS(I,ITRAN),I=1,MXEVEC)
         END DO
         WRITE(LUPRI,*)

         ! X{O} intermediate calculations:
         WRITE(LUPRI,*) 'List of X{O} intermediate calculations:'
         DO ITRAN = 1, NXITRAN
           WRITE(LUPRI,'(A,2I5,5X,(25I3,20X))') MSGDBG,
     &      (IXITRAN(I,ITRAN),I=1,1),(IXIDOTS(I,ITRAN),I=1,MXIVEC)
         END DO
         WRITE(LUPRI,*)

      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCLR_SETUP                           *
*---------------------------------------------------------------------*
