C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
C  /* Deck cc_hpqta */
      SUBROUTINE CC_HPQTA(XHPQ,ISYHPQ,XHCA,XHIK,ISYRES,LISTR,IDLSTR,
     &                    ISYTAM,LE1,LE2,WORK, LWORK)
*
******************************************************************
*
*    PURPOSE:  transform the h^(1)_pq with t1 amplitudes
*              and add to the E1* and E2 intermediates
*      
*    LE1 = true   h_c^a = - sum_k t_ck * h_ka  (for E1* inter.)
*
*    LE2 = true   h_ik^ =   sum_c h_ic * t_ck  (for E2  inter.)
*
*    Input : h_pq stored as Fock-MO of symmetry ISYHPQ
*    Output: h_ik or h_ca stored NMATIJ and NMATAB
*    Local : TAMP read from file according to LISTR,IDLSTR,ISYTAM
*
*    if (LFOCKB) store result as Fock Matrix
*
*    Sonia Coriani 06/09-1999
*     
*    Note: result is CUMULATED (initialize outside!) on XHCA,XHIK
*
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <dummy.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, XMONE= -1.0D0)
      DIMENSION XHPQ(*),XHCA(*),XHIK(*),WORK(LWORK)
      LOGICAL LE1,LE2,LFOCKB
      CHARACTER LISTR*(*), MODEL*(10)
* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
*
      IF (LOCDBG) THEN
         CALL AROUND('IN CC_HPQTA  ')
         CALL FLSHFO(LUPRI)
      ENDIF
*
      IF (LE1 .OR. LE2) THEN
        KTAMP = 1
        KEND0 = KTAMP + NT1AM(ISYTAM)
        LWRK0 = LWORK - KEND0
        IF (LWRK0.LE.NT1AM(ISYHPQ)) 
     &            CALL QUIT(' Insufficient memory in  CC_HPQTA')
        
        IOPTRES = 1
        CALL CC_RDRSP(LISTR,IDLSTR,ISYTAM,IOPTRES,MODEL,
     &                                    WORK(KTAMP),DUMMY)
*
*       Extract h_ia in memory out of h_pq (store ia and not ai!)       
*
        DO ISYMA = 1, NSYM
           ISYMI = MULD2H(ISYMA,ISYHPQ)
           DO A = 1,NVIR(ISYMA)

              KOFFPQ = IFCVIR(ISYMI,ISYMA) + NORB(ISYMI)*(A - 1) + 1
*             KOFFIA = IT1AM(ISYMA,ISYMI)  + NRHF(ISYMI)*(A - 1) + KEND0
              KOFFIA = IT1AMT(ISYMI,ISYMA) + NRHF(ISYMI)*(A - 1) + KEND0

              CALL DCOPY(NRHF(ISYMI),XHPQ(KOFFPQ),1,WORK(KOFFIA),1)
           END DO
        END DO
        
        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_HPQTA: The h_pq matrix'
           CALL CC_PRFCKMO(XHPQ,ISYHPQ)
           WRITE (LUPRI,*) 'CC_HPQTA: The h_ia matrix'
           CALL OUTPUT(WORK(KEND0),1,NRHF(1),1,NVIR(1),NRHF(1),
     &                                             NVIR(1),1,LUPRI)
           CALL FLSHFO(LUPRI)
        END IF
*
      ELSE
        WRITE (LUPRI,*) 'Warning: nothing to do CC_HPQTA. (all false?)'
        RETURN
      END IF
*
*-----------------------------------------
*     Transform h_ka to h_c^a
*-----------------------------------------
*
      IF (LE1) THEN

        DO 100 ISYMC = 1,NSYM
*
           ISYMK  = MULD2H(ISYTAM,ISYMC)
           ISYMA  = MULD2H(ISYHPQ,ISYMK)
*
           NRHFK  = MAX(NRHF(ISYMK),1)
           NVIRC  = MAX(NVIR(ISYMC),1)
*
           KOFF1  = IT1AM(ISYMC,ISYMK)   + KTAMP     !offset t_ck
*           KOFF2  = IT1AM(ISYMA,ISYMK)   + KEND0     !offset h_ka 
           KOFF2  = IT1AMT(ISYMK,ISYMA)   + KEND0     !offset h_ka 
           KOFF3  = IMATAB(ISYMC,ISYMA)  + 1         !offset output h_c^a 

           CALL DGEMM('N','N',NVIR(ISYMC),NVIR(ISYMA),NRHF(ISYMK),
     &              XMONE,WORK(KOFF1),NVIRC,WORK(KOFF2),NRHFK,
     &              ONE,XHCA(KOFF3),NVIRC)
*
  100   CONTINUE
        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_HPQTA: The h_pq matrix'
           CALL CC_PRFCKMO(XHPQ,ISYHPQ)
           WRITE (LUPRI,*) 'CC_HPQTA: The h_c^a matrix'
           CALL OUTPUT(XHCA(1),1,NVIR(1),1,NVIR(1),NVIR(1),
     *                 NVIR(1),1,LUPRI)
           CALL FLSHFO(LUPRI)
        END IF

      END IF
*
*-----------------------------------------
*     Transform h_ic to h_ik^ 
*-----------------------------------------
*
      IF (LE2) THEN
         
         DO 200 ISYMK = 1,NSYM
*
           ISYMC  = MULD2H(ISYTAM,ISYMK)
           ISYMI  = MULD2H(ISYMC,ISYHPQ)
*
           NRHFI = MAX(NRHF(ISYMI),1)
           NVIRC = MAX(NVIR(ISYMC),1)
*
*           KOFF1  = IT1AM(ISYMC,ISYMI)   + KEND0  !offset for h_ic
           KOFF1  = IT1AMT(ISYMI,ISYMC)   + KEND0  !offset for h_ic
           KOFF2  = IT1AM(ISYMC,ISYMK)   + KTAMP   !offset for t_ck
           KOFF3  = IMATIJ(ISYMI,ISYMK)  + 1       !result of transf. h_ik^ 
*
           CALL DGEMM('N','N',NRHF(ISYMI),NRHF(ISYMK),NVIR(ISYMC),
     &              ONE,WORK(KOFF1),NRHFI,WORK(KOFF2),NVIRC,
     &              ONE,XHIK(KOFF3),NRHFI)
*
  200    CONTINUE
        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_HPQTA: The h_pq matrix'
           CALL CC_PRFCKMO(XHPQ,ISYHPQ)
           WRITE (LUPRI,*) 'CC_HPQTA: The h_c^a matrix'
           CALL OUTPUT(XHIK(1),1,NRHF(1),1,NRHF(1),NRHF(1),
     *                 NRHF(1),1,LUPRI)
           CALL FLSHFO(LUPRI)
        END IF
      END IF
*
      RETURN
      END
*------------------------------------------------------------------
