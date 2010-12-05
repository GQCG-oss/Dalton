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
*=====================================================================*
C  /* Deck cc_e21con */
      SUBROUTINE CC_E21CON(XJILK,ISYJILK,XIJKD0,ISYIJKD0,
     &                     LUIJKD1,FNIJKD1,ISYIJKD1, 
     &                     XLAMDA0,ISYLAM0,XLAMDA1,ISYLAM1,LRELAX,
     &                     ITRAN, WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose:  lead calculation of E2' contribution to FBTA
*               transformed vector
*
*     ISYIJK0 = symmetry of I_{ijk;delta}
*     ISYIJK1 = symmetry of I^(1)_{ijk;delta}
*
*     Sonia Coriani, 14/09-1999
*
* Read derivative/relaxed integrals from file
* call transformation to 4 occupied
* resort to BF ordering
* add to gammaQ(rho^BFQ)
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      LOGICAL LSKIP4O
      PARAMETER (LSKIP4O = .FALSE.)

      INTEGER ISYIJKD0,ISYIJKD1,ISYLAM0,ISYLAM1,ITRAN,LWORK,IOPT
      INTEGER LUIJKD1, ISYJILK, KOFFIJKL
      LOGICAL LRELAX
      CHARACTER*8 FNIJKD1
  
#if defined (SYS_CRAY)
      REAL XJILK(*), XIJKD0(*), XLAMDA0(*), XLAMDA1(*)
      REAL ZERO, ONE, HALF, DNRM2, XNORM, WORK(LWORK)
#else
      DOUBLE PRECISION XJILK(*), XIJKD0(*), XLAMDA0(*), XLAMDA1(*)
      DOUBLE PRECISION ZERO, ONE, HALF, DNRM2, XNORM, WORK(LWORK)
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
*
      INTEGER ISYMD, ISYIJKL, NTOT, LWRK1, IOFF, N3ODMX,ISYM
      INTEGER KXIJKL, KEND1, ISYML0, ISYML1
      INTEGER KXIJK0, ISYIJK0, LENIJK0
      INTEGER KXIJK1, ISYIJK1

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'I am inside CC_E21CON: LRELAX = ', LRELAX
        CALL CC_PRLAM(XLAMDA0,XLAMDA0,ISYLAM0)
        CALL CC_PRLAM(XLAMDA1,XLAMDA1,ISYLAM1)
        CALL FLSHFO(LUPRI)
      END IF
*---------------------------------------------------------------------*
* precalculate some stuff
*---------------------------------------------------------------------*
      N3ODMX = 0
      DO ISYM  = 1, NSYM
        N3ODMX = MAX(N3ODMX,N3ODEL(ISYM))
      END DO
*     --------------------------------------
*     Begin
*     --------------------------------------
*
      ISYIJKL = MULD2H(ISYIJKD0,ISYLAM0)
      IF (LRELAX) THEN
          ISYIJKL = MULD2H(ISYIJKD0,ISYLAM1)
          IF (ISYIJKL .NE. MULD2H(ISYIJKD1,ISYLAM0))
     &       CALL QUIT('Symmetry mismatch in CC_INT4O (i)')
          IF (ISYIJKL.NE.ISYJILK)
     &       CALL QUIT('Symmetry mismatch in CC_INT4O (ii)')
      END IF
*     ------------------------------------------------------
*     Work space allocation for intermediate ijkl output
*     ------------------------------------------------------
  
      KXIJKL = 1
      KEND1  = KXIJKL + N3ORHF(ISYIJKL)
      LWRK1  = LWORK - KEND1
*
      CALL DZERO(WORK(KXIJKL),N3ORHF(ISYIJKL))
*
      IF (LSKIP4O) THEN
         WRITE (LUPRI,*) 
     &        ' CC_E21CON> I am skipping the g_ijkl integrals..'
         GO TO 200
      END IF
*     ------------------------------------------------------
*     Loop over symmetry of delta
*     ------------------------------------------------------
      KXIJK0 = 1
      DO 100 ISYMD = 1, NSYM

         IF (NBAS(ISYMD) .EQ. 0) GOTO 100

         ISYML0  = MULD2H(ISYLAM0,ISYMD)
         IF (LRELAX) THEN
            ISYML1 = MULD2H(ISYLAM1,ISYMD)
         END IF

         ISYIJK1 = MULD2H(ISYIJKD1,ISYMD)
         ISYIJK0 = MULD2H(ISYIJKD0,ISYMD)
         LENIJK0 = NMAIJK(ISYIJK0)*NBAS(ISYMD)
C---------------------------------------------------------------------
C        Work space allocation for derivative/relaxed 
C        3occ integral distribution (if LRELAX) to be read from file.
C---------------------------------------------------------------------
         IF (LRELAX) THEN
            KXIJK1 = KEND1
            KEND1  = KXIJK1 + NMAIJK(ISYIJK1)*NBAS(ISYMD)
         ELSE
            KXIJK1 = KEND1
         END IF
         LWRK1  = LWORK - KEND1
C--------------------------------------------------------------------
C        Read all integrals (ij|kdel)(1) from disc for given ISYMD.
C--------------------------------------------------------------------
         IF (LRELAX) THEN
            NTOT = NMAIJK(ISYIJK1)*NBAS(ISYMD)
            IOFF = N3ODMX*(ITRAN-1) + I3ODEL(ISYIJK1,ISYMD) + 1
            CALL GETWA2(LUIJKD1,FNIJKD1,WORK(KXIJK1),IOFF,NTOT)
         END IF
C---------------------------------------------------------------
C        Transform AO integral index delta to occupied space
C        thru a call to CC_INT4O --> return result in   XIJKL
C--------------------------------------------------------------
         IF (LOCDBG) THEN
           XNORM = DNRM2(LENIJK0,XIJKD0(KXIJK0),1)
           WRITE (LUPRI,*) 
     &        'CC_E21CON>For ISYMD:',ISYMD,' Norm IJK0:',XNORM
         END IF
C
         IF (LRELAX) THEN
            IOPT = 2
         ELSE
            IOPT = 1
         END IF
         CALL CC_INT4O(XIJKD0(KXIJK0),ISYIJK0,WORK(KXIJK1),ISYIJK1,
     &                 XLAMDA0,ISYLAM0,XLAMDA1,ISYLAM1, ISYMD,
     &                 WORK(KXIJKL),LRELAX,WORK(KEND1),LWRK1,IOPT)

         KXIJK0 = KXIJK0 + LENIJK0

 100  CONTINUE 
 200  CONTINUE 

      IF (LOCDBG) THEN
         XNORM = DNRM2(N3ORHF(ISYIJKL),WORK(KXIJKL),1)
         WRITE (LUPRI,*) 'CC_E21CON> Norm I_ijkl before sorting: ',XNORM
      END IF
       
*-----------------------------------------------------
*     Resort result integrals IJKL to M intermediate ordering
*     JILK --> added to GAMMAQ intermediate (IOPT = 2)
*--------------------------------------------------------
      IOPT = 2
      CALL CC_SORT4O(WORK(KXIJKL),ISYIJKL,XJILK,IOPT)
      IF (LOCDBG) THEN
        CALL AROUND('The IJKL integrals (+M) resorted JIL,K')
        XNORM = DNRM2(N3ORHF(ISYIJKL),XJILK,1)
        WRITE (LUPRI,*) 'Norm: ',  XNORM
      END IF

*----------------------------------
*     Close the file.
*----------------------------------

      RETURN
      END
*=====================================================================*
