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
      SUBROUTINE CCETACOR(ETAAI,ZKDIA,WORK,LWORK)
C
C     Purpose: to calculate the correction to the RHS (eta_ai and eta_aI) 
C     of the z-kappa-bar0_ai equations originating from the
C     "diagonal" z-kappa-bar0_ij and z-kappa-bar0_ab blocks
C     over a loop over delta
C
C     Written by Sonia Coriani, Fall 2001. Based on CC2_DEN
C
C     Current models: CCD, CCSD, CCSD(T)
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION ETAAI(*), ZKDIA(*), WORK(LWORK)
#include <ccorb.h>
#include <ccisao.h>
#include <r12int.h>
#include <inftap.h>
#include <blocks.h>
#include <ccfield.h>
#include <ccsdinp.h>
#include <ccinftap.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <ccfro.h>
C
      CHARACTER MODEL*10
      CHARACTER NAME1*8
      CHARACTER NAME2*8

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.FALSE.)
C
      CALL QENTER('CCETACOR')

C--------------------------------------------------------------
C Initialize work space
C--------------------------------------------------------------

      KCMO = 1
      KEND1 = KCMO + NLAMDS
      LWRK1 = LWORK - KEND1
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
C 
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
      CALL GPCLOSE (LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)

      IF (FROIMP) THEN
         KCMOF = KEND1
         KEND1 = KCMOF + NLAMDS
         LWRK1 = LWORK - KEND1
         CALL CMO_ALL(WORK(KCMOF),WORK(KEND1),LWRK1)
      END IF
C
C---------------------------------------------------------------
C     Start the loop for the 2-electron integrals and densities.
C---------------------------------------------------------------
C
      KENDS2 = KEND1
      LWRKS2 = LWRK1
C
      IF (DIRECT) THEN
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      ICDEL1 = 0
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
c              DTIME  = SECOND()
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRINT)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                        WORK(KODCL1),WORK(KODCL2),
     *                        WORK(KODBC1),WORK(KODBC2),
     *                        WORK(KRDBC1),WORK(KRDBC2),
     *                        WORK(KODPP1),WORK(KODPP2),
     *                        WORK(KRDPP1),WORK(KRDPP2),
     *                        WORK(KCCFB1),WORK(KINDXB),
     *                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
c              DTIME  = SECOND() - DTIME
c              TIMHE2 = TIMHE2   + DTIME
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC2_DEN')
               END IF
C
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
C
               ISYDEN = ISYMD
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
               KXINT  = KEND1
               KEND2  = KXINT  + NDISAO(ISYDIS)
               LWRK2  = LWORK  - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:',KEND2
                  CALL QUIT(
     *         'Insufficient space for allocation 2 in CCETAAICORR')
               ENDIF
C
C--------------------------------------------
C              Read AO integral distribution.
C--------------------------------------------
C
               AUTIME = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C
C--------------------------------------------------------------
C              Calculate correction terms to eta_ai originating
C              from kbar_ij and kbar_ab.
C--------------------------------------------------------------
C
               KDSRHF = KEND2
               K3OINT = KDSRHF + NDSRHF(ISYMD)
               KEND3  = K3OINT + NMAIJK(ISYDIS)
               IF (FROIMP) THEN
                  KDSFRO = KEND3
                  KOFOIN = KDSFRO + NDSFRO(ISYDIS)
                  KEND3  = KOFOIN + NOFROO(ISYDIS)
               ENDIF
               LWRK3  = LWORK  - KEND3
C
               IF (LWRK3 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND3,'Available : ',LWORK
                  CALL QUIT(
     *              'Insufficient core for integrals in CCETAC')
               ENDIF
C
C--------------------------------------------------------------------
C              Transform one index in the integral batch to occupied.
C--------------------------------------------------------------------
C
               CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KCMO),
     *                        ISYMOP,WORK(KEND3),LWRK3,ISYDIS)
csonia
C
C------------------------------------------------------------------
C              Transform one index in the integral batch to frozen.
C------------------------------------------------------------------
C
               IF (FROIMP) THEN
C
                  CALL CC_GTOFRO(WORK(KXINT),WORK(KDSFRO),WORK(KCMOF),
     *                           WORK(KEND3),LWRK3,ISYDIS)
C
C--------------------------------------------------------------
C                 Calculate integral batch (cor fro | cor del).
C--------------------------------------------------------------
C
                  CALL CC_OFROIN(WORK(KDSRHF),WORK(KOFOIN),WORK(KCMOF),
     *                           WORK(KEND3),LWRK3,ISYDIS)
C
C-------------------------------------------------------------------------
C                 Calculate correction to eta_aI from kappabar_ab
C-------------------------------------------------------------------------
C
                  if (.true.) then
                  KAFROI = 1 + NMATIJ(1) + NMATAB(1) + 2*NT1AMX 
                  CALL MP2_EIDV1(ZKDIA(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(NMATIJ(1)+1),WORK(KCMOF),
     *                           WORK(KEND3),LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDV2(ZKDIA(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(NMATIJ(1)+1),WORK(KCMOF),
     *                           WORK(KEND3),LWRK3,IDEL,ISYMD)
C
C-------------------------------------------------------------------------
C                 Calculate correction to eta_aI from kappabar_ij
C-------------------------------------------------------------------------
C
                  CALL MP2_EIDC1(ZKDIA(KAFROI),WORK(KDSFRO),
     *                           ZKDIA(1),WORK(KCMOF),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)
C
                  CALL MP2_EIDC2(ZKDIA(KAFROI),WORK(KOFOIN),
     *                           ZKDIA(1),WORK(KCMOF),WORK(KEND3),
     *                           LWRK3,IDEL,ISYMD)

                  end if
               END IF 
csonia
C
C-------------------------------------------------------------------
C              Calculate integral batch with three occupied indices.
C              Since LUDUMLOCAL < 0 the integrals are not written 
C              to disk.
C-------------------------------------------------------------------
C
               LUDUMLOCAL = -10
               CALL CC_INT3O(WORK(K3OINT),WORK(KDSRHF),WORK(KCMO),
     *                         ISYMOP,WORK(KCMO),WORK(KEND3),LWRK3,
     *                         IDEL,ISYMD,LUDUMLOCAL,'DUMMY')
C
C-------------------------------------------------------
C              Calculate the correction terms to eta_ai
C              from kappabar_ab and kappabar_ij
C-------------------------------------------------------
C
               CALL MP2_YTV(ETAAI,ZKDIA(NMATIJ(1)+1),WORK(KDSRHF),
     *                      WORK(KCMO),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
               CALL MP2_NXY(ETAAI,ZKDIA(1),ZKDIA(NMATIJ(1)+1),
     *                      WORK(K3OINT),WORK(KDSRHF),WORK(KCMO),
     *                      WORK(KEND3),LWRK3,IDEL,ISYMD)
C
               CALL MP2_XTO(ETAAI,ZKDIA(1),WORK(K3OINT),
     *                      WORK(KCMO),WORK(KEND3),LWRK3,IDEL,ISYMD)
C
               AUTIME = SECOND() - AUTIME
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE

C-----------------------------------------------------
C     Some print out
C-----------------------------------------------------
      IF (LOCDBG) THEN
         ETAKAN = DDOT(NALLAI(1),ETAAI,1,ETAAI,1)
         WRITE(LUPRI,*) 'CCETACOR: eta_ai '
         WRITE(LUPRI,*) 'Norm of corrected occ-vir block:', ETAKAN
      ENDIF
C
C-----------------------
C     Regain work space.
C-----------------------
C
      KEND1 = KENDS2
      LWRK1 = LWRKS2
C
C--------------------------------------------------------------- 
      CALL QEXIT('CCETACOR')
      RETURN
      END
