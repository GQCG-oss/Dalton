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
*=====================================================================*
c /* deck ccmcd_setup */
*=====================================================================*
      SUBROUTINE CCMCD_SETUP(MXTRAN, MXVEC,
     &                       IGTRAN, IGDOTS, NGTRAN,
     &                       IFATRAN,IFADOTS,NFATRAN,
     &                       IEATRAN,IEADOTS,NEATRAN,
     &                       IEATRA1,IEADOT1,NEATRA1,
     &                       IEATRA2,IEADOT2,NEATRA2,
     &                       IFTRAN, IFDOTS, NFTRAN,
     &                       IFTRA1, IFDOT1, NFTRA1,
     &                       IXE2TRA,IX2DOTS,IE2DOTS,NXE2TRA,
     &                       IX2TRAN,IX2DOT1,NX2TRAN,
     &                       WORK, LWORK  )
*---------------------------------------------------------------------*
*
*    Purpose: set up for MCD section
*                - list of G matrix transformations 
*                - list of F{Op} matrix transformations 
*                - list of ETA{Op} vector calculations 
*    if (LUSEPL1) also
*                - list of B transformations
*
*    Written by:       Sonia Coriani 1997-98
*    Restructured by:  Sonia Coriani 26/01-2000
*    Orbital relaxation for B operator introduced March 2000. Sonia
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccmcdinf.h>
#include <ccroper.h>
#include <ccexci.h>
#include <cclists.h>

* local parameters:
      CHARACTER*(23) MSGDBG
      PARAMETER (MSGDBG = '[debug] CCMCD_SETUP> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER MXVEC, MXTRAN, LWORK

#if defined (SYS_CRAY)
      REAL ZERO, D05, TWO
      REAL EIGVAF, WORK(LWORK)
#else
      DOUBLE PRECISION ZERO, D05, TWO
      DOUBLE PRECISION EIGVAF, WORK(LWORK)
#endif
      PARAMETER (ZERO = 0.0D00, D05  = 0.5D00, TWO = 2.0D00)

      INTEGER IGTRAN(MXDIM_GTRAN,MXTRAN)
      INTEGER IGDOTS(MXVEC,MXTRAN)

      INTEGER IFATRAN(MXDIM_FATRAN,MXTRAN)
      INTEGER IFADOTS(MXVEC,MXTRAN)

      INTEGER IEATRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IEADOTS(MXVEC,MXTRAN)

      INTEGER IEATRA1(MXDIM_XEVEC,MXTRAN)
      INTEGER IEADOT1(MXVEC,MXTRAN)

      INTEGER IEATRA2(MXDIM_XEVEC,MXTRAN)
      INTEGER IEADOT2(MXVEC,MXTRAN)

      INTEGER IFTRAN(MXDIM_FTRAN,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)

      INTEGER IFTRA1(MXDIM_FTRAN,MXTRAN)
      INTEGER IFDOT1(MXVEC,MXTRAN)

      INTEGER IXE2TRA(MXDIM_XEVEC,MXTRAN)
      INTEGER IX2DOTS(MXVEC,MXTRAN), IE2DOTS(MXVEC,MXTRAN)
      INTEGER IX2TRAN(MXDIM_XEVEC,MXTRAN)
      INTEGER IX2DOT1(MXVEC,MXTRAN)

      INTEGER NGTRAN, NFATRAN, NEATRAN 
      INTEGER NEATRA1,NEATRA2,NFTRAN,NFTRA1,NXE2TRA,NX2TRAN
      INTEGER NBTERMS

      INTEGER ISYMA, ISYMB, ISYMC, ISYMAB, ISYMSF, ISYSOP, ISGNSOP
      INTEGER IFREQ, IOPER

      CHARACTER*8 LABELA,LABELB,LABELC,LABSOP
      INTEGER IOPERA, IOPERB, IOPERC, IOPSOP
      INTEGER ITAMPA, ITAMPB, ITAMPC
      INTEGER IKAPA,IKAPB
      INTEGER IVEC, ITRAN, I, ISTATE, INUM, K
      INTEGER ISTATF, IEXCIF
      INTEGER IM1F,IZETAA,IPL1A,IZETA0
  
      LOGICAL LORXA,LORXB,LORXC,LPDBSA,LPDBSB,LPDBSC,LRELAX,LPROJ

* external functions:
      INTEGER IR1TAMP
      INTEGER IR1KAPPA
      INTEGER IPL1ZETA
      INTEGER ILRMAMP
      INTEGER IROPER
      INTEGER IRHSR2

*---------------------------------------------------------------------*
* initializations:
*---------------------------------------------------------------------*
      DO ITRAN = 1, MXTRAN
        DO I = 1, MXDIM_GTRAN
          IGTRAN(I,ITRAN) = 0
        END DO
        DO I = 1, MXDIM_FATRAN
          IFATRAN(I,ITRAN) = 0
        END DO
        DO I = 1, MXDIM_XEVEC
          IEATRAN(I,ITRAN) = 0
          IEATRA1(I,ITRAN) = 0
          IEATRA2(I,ITRAN) = 0
          IXE2TRA(I,ITRAN) = 0
          IX2TRAN(I,ITRAN) = 0
        END DO
        IEATRAN(3,ITRAN) = -1
        IEATRA1(3,ITRAN) = -1
        IEATRA2(3,ITRAN) = -1
        IXE2TRA(3,ITRAN) = -1
        IX2TRAN(3,ITRAN) = -1
        IX2TRAN(4,ITRAN) = -1

        DO I = 1, MXDIM_FTRAN
          IFTRAN(I,ITRAN) = 0
          IFTRA1(I,ITRAN) = 0
        END DO

        DO IVEC  = 1, MXVEC
          IGDOTS(IVEC, ITRAN) = 0
          IFADOTS(IVEC,ITRAN) = 0
          IEADOTS(IVEC,ITRAN) = 0
          IEADOT1(IVEC,ITRAN) = 0
          IEADOT2(IVEC,ITRAN) = 0
          IFDOTS(IVEC,ITRAN)  = 0
          IFDOT1(IVEC,ITRAN)  = 0
          IX2DOTS(IVEC,ITRAN)  = 0
          IE2DOTS(IVEC,ITRAN)  = 0
          IX2DOT1(IVEC,ITRAN)  = 0
        END DO
      END DO

      NGTRAN  = 0
      NFATRAN = 0
      NEATRAN = 0
      NEATRA1 = 0
      NEATRA2 = 0
      NFTRAN  = 0
      NFTRA1  = 0
      NXE2TRA = 0
      NX2TRAN = 0

      NBTERMS = 0   !nr operator triplets matching selected exc. state.
      
      LPROJ = .FALSE.
 
*---------------------------------------------------------------------*
* loop # triples A,B,C
*     loop # states 
*---------------------------------------------------------------------*
      
      DO IOPER = 1, NMCDOPER

         IOPERA = IAMCDOP(IOPER)
         IOPERB = IBMCDOP(IOPER)
         IOPERC = ICMCDOP(IOPER)
 
         ISYMA  = ISYOPR(IOPERA)
         ISYMB  = ISYOPR(IOPERB)
         ISYMC  = ISYOPR(IOPERC)

         ISYMAB = MULD2H(ISYMA,ISYMB)

         IF (ISYMAB.EQ.ISYMC) THEN
           
           LABELA = LBLOPR(IOPERA)
           LABELB = LBLOPR(IOPERB)
           LABELC = LBLOPR(IOPERC)

           LORXA  = LAMCDRX(IOPER)
           LORXB  = LBMCDRX(IOPER)
           LORXC  = LCMCDRX(IOPER)

           LPDBSA = LPDBSOP(IOPERA)
           LPDBSB = LPDBSOP(IOPERB)
           LPDBSC = LPDBSOP(IOPERC)

           LRELAX = (LORXA.OR.LPDBSA.OR.
     &               LORXB.OR.LPDBSB.OR.LORXC.OR.LPDBSC)

           DO ISTATE = 1, NMCDST
              ISYMSF = IMCDSTSY(ISTATE)
              ISTATF = IMCDSTNR(ISTATE)
              IEXCIF = ISYOFE(ISYMSF) + ISTATF
              EIGVAF = EIGVAL(IEXCIF)

              IF (ISYMSF.EQ.ISYMC) THEN

                 NBTERMS = NBTERMS + 1
c         
                 ITAMPA = IR1TAMP(LABELA,LORXA,-EIGVAF,ISYMA)
                 ITAMPB = IR1TAMP(LABELB,LORXB,ZERO,ISYMB)
                 IKAPA  = 0
                 IKAPB  = 0
                 IF (LORXA) IKAPA = IR1KAPPA(LABELA,-EIGVAF,ISYMA) 
                 IF (LORXB) IKAPB = IR1KAPPA(LABELB,ZERO,ISYMB) 
    
                 IF (LOCDBG) THEN
                    WRITE (LUPRI,*) 'CCMCD_SETUP> LABELA: ', LABELA, 
     &                      '  ISYMA: ',ISYMA,
     &                      '  LORXA: ', LORXA, ' LPDBSA: ', LPDBSA 
                    WRITE (LUPRI,*) 'CCMCD_SETUP> LABELB: ', LABELB,
     &                      '  ISYMB: ',ISYMB,
     &                      '  LORXB: ', LORXB, ' LPDBSB: ', LPDBSB 
                    WRITE (LUPRI,*) 'LUSEPL1: ', LUSEPL1
                    CALL FLSHFO(LUPRI)
                 ENDIF
                 
                 IZETA0 = 0
                 IM1F   = ILRMAMP(IEXCIF,EIGVAF,ISYMC)

                 IF (LUSEPL1) THEN
                    IF (ISYMB.EQ.1) LPROJ = .TRUE.
                    IPL1A  = IPL1ZETA(LABELA,LORXA,-EIGVAF,ISYMA,
     &                                LPROJ,IEXCIF,EIGVAF,ISYMC)
                 END IF
*---------------------------------------------------------------------*
* set up (asymmetric) list of G matrix transformations:(G*Ta*Tb)*Ef
*---------------------------------------------------------------------*
                  CALL CC_SETG112(IGTRAN,IGDOTS,MXTRAN,MXVEC,IZETA0,
     &                            ITAMPA,ITAMPB,IEXCIF,ITRAN,IVEC)
                  NGTRAN = MAX(NGTRAN,ITRAN)
*---------------------------------------------------------------------*
* set up (asymmetric) list of F{O} matrix transformations 
*---------------------------------------------------------------------*
* (F{A}*TB)*Ef
                  CALL CC_SETFB12(IFATRAN,IFADOTS,MXTRAN,MXVEC,IZETA0,
     &                           IOPERA,IKAPA,ITAMPB,IEXCIF,ITRAN,IVEC)
                  NFATRAN = MAX(NFATRAN,ITRAN)
* (F{B}*TA)*Ef
                  CALL CC_SETFB12(IFATRAN,IFADOTS,MXTRAN,MXVEC,IZETA0,
     &                           IOPERB,IKAPB,ITAMPA,IEXCIF,ITRAN,IVEC)
                  NFATRAN = MAX(NFATRAN,ITRAN)

*---------------------------------------------------------------------*
* set up list of ETA{O} vector calculations
*---------------------------------------------------------------------*
                  CALL CC_SETXE('Eta',IEATRAN,IEADOTS,
     &                           MXTRAN,MXVEC,
     &                           IEXCIF,IOPERA,IKAPA,0,0,0,ITAMPB,
     &                           ITRAN,IVEC)
                  NEATRAN = MAX(NEATRAN,ITRAN)

                  IF (LUSEPL1) THEN
                     ! two permutations
                     CALL CC_SETXE('Eta',IEATRA1,IEADOT1,
     &                              MXTRAN,MXVEC,
     &                              IM1F,IOPERA,IKAPA,0,0,0,ITAMPB,
     &                              ITRAN,IVEC)
                     NEATRA1 = MAX(NEATRA1,ITRAN)
                     CALL CC_SETXE('Eta',IEATRA1,IEADOT1,
     &                              MXTRAN,MXVEC,
     &                              IM1F,IOPERB,IKAPB,0,0,0,ITAMPA,
     &                              ITRAN,IVEC)
                     NEATRA1 = MAX(NEATRA1,ITRAN)
                     ! 
                     CALL CC_SETXE('Eta',IEATRA2,IEADOT2,
     &                             MXTRAN,MXVEC,
     &                             IPL1A,IOPERB,IKAPB,0,0,0,IEXCIF,
     &                             ITRAN,IVEC)
                     NEATRA2 = MAX(NEATRA2,ITRAN)
                  END IF
*---------------------------------------------------------------------*
* set up (asymmetric) list of generalized F matrix transformations 
*---------------------------------------------------------------------*
                  IF (LUSEPL1) THEN

* M^f*(B*T^A*T^B) = (M^f*B*T^A)*T^B = (F[M^f]*T^A)*T^B

                     CALL CCQR_SETF(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                              IM1F,ITAMPA,ITAMPB,ITRAN,IVEC)
                     NFTRAN = MAX(NFTRAN,ITRAN)

* Z^A*(B*T^B*E^f) = (Z^A*B*T^B)*E^f = (F[Z^A]*T^B)*E^f

                     CALL CC_SETF12(IFTRA1,IFDOT1,MXTRAN,MXVEC,
     &                              IPL1A,ITAMPB,IEXCIF,ITRAN,IVEC)
                     NFTRA1 = MAX(NFTRA1,ITRAN)

                  END IF
*---------------------------------------------------------------------*
* set up list of Xksi{O2} and Eta{O2} vector calculations
* Please note: we only use RELAXED B. 
*---------------------------------------------------------------------*
                  IF (LUSEPL1.AND.(LPDBSB.OR.LORXB)) THEN
                     CALL CC_FIND_SO_OP(LBLOPR(IOPERA),LBLOPR(IOPERB),
     &                                  LABSOP,ISYSOP,ISGNSOP,INUM,
     &                                  WORK,LWORK)
                     IOPSOP = IROPER(LABSOP,ISYSOP)
                     !Xksi{O2} dot-multiplied by M^f
                     CALL CC_SETXE('Xi ',IXE2TRA,IX2DOTS,MXTRAN,MXVEC,
     &                            IZETA0,IOPSOP,IKAPA,IKAPB,0,0,IM1F,
     &                             ITRAN,IVEC)
                     NXE2TRA = MAX(NXE2TRA,ITRAN)
                     !Eta{O2} dot-multiplied by E^f 
                     CALL CC_SETXE('Eta',IXE2TRA,IE2DOTS,MXTRAN,MXVEC,
     &                            IZETA0,IOPSOP,IKAPA,IKAPB,0,0,IEXCIF,
     &                             ITRAN,IVEC)
                     NXE2TRA = MAX(NXE2TRA,ITRAN)
                     !Xksi{O2} dot-multiplied by Ebar^f
                     CALL CC_SETXE('Xi ',IX2TRAN,IX2DOT1,MXTRAN,MXVEC,
     &                            IZETA0,IOPSOP,IKAPA,IKAPB,0,0,IEXCIF,
     &                             ITRAN,IVEC)
                     NX2TRAN = MAX(NX2TRAN,ITRAN)
                  ELSE IF (LPDBSB.OR.LORXB) THEN
                    CALL QUIT(
     &                 'PDBS only programmed for B and LUSEPL1 = true')
                  END IF
*---------------------------------------------------------------------*
* end loops
*---------------------------------------------------------------------*
               END IF                ! <F| = C
            END DO                   ! do # selected MCD states 
         END IF                      ! endif matching sym AB = C
      END DO                         !do loop on # triples ABC

*---------------------------------------------------------------------*
* print the lists: 
*---------------------------------------------------------------------*
* general statistics:
      WRITE(LUPRI,'(/,/3X,A,I3,A)') 'For the requested ',NBTERMS,
     &      ' contributions to B terms of magn.circ.dichr. '
      WRITE(LUPRI,'((8X,A,I3,A))') 
     &   ' - ',NGTRAN,  ' G matrix transformations ',
     &   ' - ',NFATRAN, ' F{O} matrix transformations ',
     &   ' - ',NEATRAN, ' ETA{O} vector calculations ', 
     &   ' - ',NEATRA1, ' extra ETA{O} vector calculations (1)', 
     &   ' - ',NEATRA2, ' extra ETA{O} vector calculations (2)', 
     &   ' - ',NFTRAN,  ' extra B matrix transformation (1)', 
     &   ' - ',NFTRA1,  ' extra B matrix transformation (2)(1)' 
      WRITE(LUPRI,'(3X,A,/,/)') 'will be performed.'


* G matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of G matrix transformations:'
        DO ITRAN = 1, NGTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IGTRAN(I,ITRAN),I=1,3),(IGDOTS(I,ITRAN),I=1,MXVEC)
        END DO
        WRITE (LUPRI,*)
      END IF

* F{O} matrix transformations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of F{O} matrix transformations:'
        DO ITRAN = 1, NFATRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IFATRAN(I,ITRAN),I=1,3),(IFADOTS(I,ITRAN),I=1,MXVEC)
        END DO
        WRITE (LUPRI,*)
      END IF

* ETA{O} vector calculations:
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'List of ETA{O} vector calculations:'
        DO ITRAN = 1, NEATRAN
          WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IEATRAN(I,ITRAN),I=1,5),(IEADOTS(I,ITRAN),I=1,MXVEC)
        END DO
        WRITE (LUPRI,*)
        IF (LUSEPL1) THEN
* extra ETA{O} vector calculations:
          WRITE (LUPRI,*) 
     &         'List of additional ETA{O} vector calculations:'
          DO ITRAN = 1, NEATRA1
            WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IEATRA1(I,ITRAN),I=1,5),(IEADOT1(I,ITRAN),I=1,MXVEC)
          END DO
          WRITE (LUPRI,*)
          DO ITRAN = 1, NEATRA2
            WRITE(LUPRI,'(A,2I5,5X,(12I5,20X))') MSGDBG,
     &     (IEATRA2(I,ITRAN),I=1,5),(IEADOT2(I,ITRAN),I=1,MXVEC)
          END DO
          WRITE (LUPRI,*)
        END IF
      END IF

* Additional B (generalized F) matrix transformations:
      IF (LOCDBG.AND.LUSEPL1) THEN
        WRITE (LUPRI,*) 'List of B/F matrix transformations:'
        DO ITRAN = 1, NFTRAN
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IFTRAN(I,ITRAN),I=1,3),(IFDOTS(I,ITRAN),I=1,MXVEC)
        END DO
        WRITE (LUPRI,*)
        DO ITRAN = 1, NFTRA1
          WRITE(LUPRI,'(A,3I5,5X,(12I5,20X))') MSGDBG,
     &     (IFTRA1(I,ITRAN),I=1,3),(IFDOT1(I,ITRAN),I=1,MXVEC)
        END DO
        WRITE (LUPRI,*)
      END IF

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCMCD_SETUP                          *
*---------------------------------------------------------------------*
