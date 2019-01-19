!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
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
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "dummy.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, XMONE= -1.0D0)
      DIMENSION XHPQ(*),XHCA(*),XHIK(*),WORK(LWORK)
      LOGICAL LE1,LE2,LFOCKB
      CHARACTER LISTR*(*), MODEL*(10)
* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
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
