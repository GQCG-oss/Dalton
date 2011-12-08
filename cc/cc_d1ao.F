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
C  /* Deck cc_d1ao */
      SUBROUTINE CC_D1AO(IPDD,R12PRP,AODEN,ZKAM,WORK,LWORK,MODEL,LLIST,
     &                   ILLNR,NO,FNDPTIA, FNDPTIA2, FNDPTAB, FNDPTIJ)
C
C     Written by Asger Halkier January 1996.
C
C     Version: 1.0
C
C     Purpose: To calculate the one electron Coupled Cluster
C              density matrix and returning it backtransformed
C              to AO-basis in AODEN!
C
C     Current models: MP2, CCD, CCSD, CCSD(T), CC3
C
C     NO = .true.  compute natural occupation numbers
C
C     OC: 15-4-1997: general llist.
C
C     May 1998: MP2 Frozen core density by Asger Halkier.
C
C     Spring 2002 : CC3 section by K. Hald.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      logical r12prp
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION AODEN(*), ZKAM(*), WORK(LWORK)
#include <ccorb.h>
#include <infind.h> 
#include <inftap.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <ccfop.h>
#include <ccfro.h>
C
      CHARACTER MODEL*10,MODDUM*10
      CHARACTER LLIST*(*)
      CHARACTER*(*) FNDPTIA, FNDPTAB, FNDPTIJ, FNDPTIA2
      LOGICAL NO
C
C----------------------------------------
C     Initialization of timing parameter.
C----------------------------------------
C
      CALL QENTER('CC_D1AO')
      TIMTOT = ZERO
      TIMTOT = SECOND()
C
      IF ((CC2) .AND. (RELORB)) THEN
C
         CALL CC2_D1AO(AODEN,ZKAM,WORK,LWORK)
         CALL QEXIT('CC_D1AO')
         RETURN
C
      ELSE IF (MP2) THEN
C
         KONEAI = 1
         KONEAB = KONEAI + NT1AMX
         KONEIJ = KONEAB + NMATAB(1)
         KONEIA = KONEIJ + NMATIJ(1)
         KT1AM  = KONEIA + NT1AMX
         KLAMDH = KT1AM  + NT1AMX
         KLAMDP = KLAMDH + NLAMDT
         KRMAT  = KLAMDP + NLAMDT 
         KEND1  = KRMAT  + NMATIJ(1)
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'MP2 - Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient memory for work '//
     &           'allocation in CC_D1AO for MP2')
         ENDIF
C
C-----------------------------------------------------------------
C        Initialize arrays (note that the t1-amplitudes are zero).
C-----------------------------------------------------------------
C
         CALL DZERO(WORK(KONEAI),NT1AMX)
         CALL DZERO(WORK(KONEAB),NMATAB(1))
         CALL DZERO(WORK(KONEIJ),NMATIJ(1))
         CALL DZERO(WORK(KONEIA),NT1AMX)
         CALL DZERO(WORK(KT1AM),NT1AMX)
         CALL DZERO(WORK(KRMAT),NMATIJ(1))
C
C--------------------------------------------------
C        Set up the relaxation part of the density.
C--------------------------------------------------
C
         CALL DCOPY(NMATIJ(1),ZKAM(1),1,WORK(KONEIJ),1)
         CALL DCOPY(NMATAB(1),ZKAM(NMATIJ(1)+1),1,WORK(KONEAB),1)
         CALL DCOPY(NT1AMX,ZKAM(NMATIJ(1)+NMATAB(1)+1),1,WORK(KONEAI),1)
         CALL DCOPY(NT1AMX,WORK(KONEAI),1,WORK(KONEIA),1)
C
C-------------------------------------
C        Add the Hartree-Fock density.
C        Add R12/A' contribution (Elena)
C-------------------------------------
C
celena
         IF (R12PRP .AND. (IPDD .EQ. 3 .OR. IPDD .EQ. 5)) THEN
            lunit = -1
            IF (IPDD .EQ. 3) THEN
               call gpopen(lunit,'CCR12YIJ','unknown',' ','unformatted',
     &                     idummy,.false.)
             ELSEIF (IPDD .EQ. 5) THEN
               call gpopen(lunit,'CCR12YIJB','unknown',' ',
     &                     'unformatted',idummy,.false.)
             ENDIF
        
           do isymaj = 1,nsym
              isymij = isymaj
              ncvij = 0
              do isyma = 1,nsym
                 isymj = muld2h(isymaj,isyma)
                 isymi = muld2h(isymij,isymj)
                 ncvij = ncvij +  nrhf(isyma)*nrhf(isymi)
              enddo
            enddo
        
            read(lunit)(work(krmat+i-1),i=1,ncvij)
            call gpclose(lunit,'KEEP')

            DO ISYMI = 1,NSYM
               ISYMJ = ISYMI
               DO  J = 1,NRHF(ISYMJ)
                   DO  I = J+1,NRHF(ISYMI)
                       NIJ   = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)
     *                          *(J - 1) + I
                       NJI   = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)
     *                          *(I - 1) + J
                       WORK(KONEIJ + NIJ - 1) = WORK(KONEIJ + NIJ - 1)
     *                          + work(krmat + nij - 1)
                       WORK(KONEIJ + NJI - 1) = WORK(KONEIJ + NJI - 1)
     *                          + work(krmat + nji - 1)
C
                   ENDDO   
               ENDDO   
            ENDDO    
         
         ENDIF

celena

         DO 80 ISYM = 1,NSYM
            DO 90 I = 1,NRHF(ISYM)
               NII = IMATIJ(ISYM,ISYM) + NRHF(ISYM)*(I - 1) + I
               WORK(KONEIJ + NII - 1) = WORK(KONEIJ + NII - 1) + TWO
               if (r12prp .and. (ipdd .eq. 3 .or. ipdd .eq. 5)) then
                  WORK(KONEIJ + NII - 1) = WORK(KONEIJ + NII - 1) + 
     &                   1.0d0*work(krmat + nii - 1)              
               endif
   90       CONTINUE
   80    CONTINUE
C
C----------------------------------------
C        Calculate MO coefficient matrix.
C----------------------------------------
C
         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     *               WORK(KEND1),LWRK1)
C
C--------------------------------------
C        Transform density to AO basis.
C--------------------------------------
C
         CALL DZERO(AODEN,N2BST(1))
C
         ISDEN = 1
         CALL CC_DENAO(AODEN,1,WORK(KONEAI),WORK(KONEAB),
     *                 WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *                 WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
C----------------------------------------------------
C        Add frozen core contributions to AO density.
C----------------------------------------------------
C
         IF (FROIMP) THEN
C
            KOFFAI = NMATIJ(1) + NMATAB(1) + 2*NT1AMX + 1
            KOFFIA = KOFFAI    + NT1FRO(1)
            KOFFIJ = KOFFIA    + NT1FRO(1)
            KOFFJI = KOFFIJ    + NCOFRO(1)
C
            ISDEN = 1
            ICON  = 1
            CALL CC_D1FCB(AODEN,ZKAM(KOFFIJ),ZKAM(KOFFJI),ZKAM(KOFFAI),
     *                    ZKAM(KOFFIA),WORK(KEND1),LWRK1,ISDEN,ICON)
C
         ENDIF
C
      ELSE IF ((CCSD .OR. CCPT .OR. CC3) 
     *          .OR. ((CC2) .AND. (.NOT. RELORB))) THEN
C
C-------------------------------------------------------
C        Both zeta- and t-vectors are totally symmetric.
C-------------------------------------------------------
C
         ISYMTR = 1
         ISYMOP = 1
C
C--------------------------------------
C        Initial work space allocation.
C--------------------------------------
C
         NHELP = NT2AMX
C
         KZ2AM  = 1
         KT2AM  = KZ2AM  + NT2SQ(1)
         KT2AMT = KT2AM  + NT2AMX
         KLAMDP = KT2AMT + NHELP
         KLAMDH = KLAMDP + NLAMDT
         KT1AM  = KLAMDH + NLAMDT
         KZ1AM  = KT1AM  + NT1AMX
         KEND1  = KZ1AM  + NT1AMX
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient core for initial '//
     &           'allocation in CC_D1AO')
         ENDIF
C
C-------------------------------------------
C        Read zero'th order zeta amplitudes.
C-------------------------------------------
C
         ISYML = 1
         IOPT = 3
         CALL CC_RDRSP(LLIST,ILLNR,ISYML,IOPT,MODDUM,
     *                 WORK(KZ1AM),WORK(KT2AM))
         IF ( IPRINT .GT. 10 ) THEN
            RHO1N = DDOT(NT1AM(ISYML),WORK(KZ1AM),1,WORK(KZ1AM),1)
            RHO2N = DDOT(NT2AM(ISYML),WORK(KT2AM),1,WORK(KT2AM),1)
            WRITE(LUPRI,*) 'CC_D1AO: Norm of L. vector :',RHO1N
            WRITE(LUPRI,*) 'CC_D1AO: Norm of L. vector :',RHO2N
         ENDIF

C
C-----------------------------------
C        Square up zeta2 amplitudes.
C-----------------------------------
C
         CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C----------------------------------------------
C        Read zero'th order cluster amplitudes.
C----------------------------------------------
C
         IOPT = 3
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
C
C---------------------------------------------------
C        Zero out single vectors in CCD-calculation.
C---------------------------------------------------
C
         IF (CCD) THEN
            CALL DZERO(WORK(KT1AM),NT1AMX)
            CALL DZERO(WORK(KZ1AM),NT1AMX)
         ENDIF
C
C-------------------------------------
C        Calculate the lamda matrices.
C-------------------------------------
C
         CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *               LWRK1)
C
C------------------------------------------
C        Set up 2C-E of cluster amplitudes.
C------------------------------------------
C
         CALL DCOPY(NT2AMX,WORK(KT2AM),1,WORK(KT2AMT),1)
         ISYOPE = 1
         IOPTTCME = 1
         CALL CCSD_TCMEPK(WORK(KT2AMT),1.0D0,ISYOPE,IOPTTCME)
C
C----------------------------------
C        Work space allocation one.
C        Note that D(ai)=ZETA(ai) &
C        D(ia) is stored transposed
C----------------------------------
C
         KONEAI = KZ1AM
         KONEAB = KONEAI + NT1AMX
         KONEIJ = KONEAB + NMATAB(1)
         KONEIA = KONEIJ + NMATIJ(1)
         KXMAT  = KONEIA + NT1AMX
         KYMAT  = KXMAT  + NMATIJ(1)
         KCMO   = KYMAT  + NMATAB(1)
         KEND1  = KCMO   + NLAMDS
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'CC_D1AO; Available:',LWORK, 'Needed:',KEND1
            CALL QUIT(
     &           'Insufficient memory for first allocation in CC_D1AO')
         ENDIF
C
C---------------------------------------------------------
C        Initialize remaining one electron density arrays.
C---------------------------------------------------------
C
         CALL DZERO(WORK(KONEAB),NMATAB(1))
         CALL DZERO(WORK(KONEIJ),NMATIJ(1))
         CALL DZERO(WORK(KONEIA),NT1AMX)
C
C-----------------------------------------------------------
C        Calculate X-intermediate of zeta- and t-amplitudes.
C-----------------------------------------------------------
C
         CALL CC_XI(WORK(KXMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *              WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------
C        Calculate Y-intermediate of zeta- and t-amplitudes.
C-----------------------------------------------------------
C
         CALL CC_YI(WORK(KYMAT),WORK(KZ2AM),ISYMTR,WORK(KT2AM),ISYMOP,
     *              WORK(KEND1),LWRK1)
C
C-----------------------------------------------------------------
C        Construct three remaining blocks of one electron density.
C-----------------------------------------------------------------
C
         CALL DCOPY(NMATAB(1),WORK(KYMAT),1,WORK(KONEAB),1)
         CALL CC_EITR(WORK(KONEAB),WORK(KONEIJ),WORK(KEND1),LWRK1,1)
         CALL DIJGEN(WORK(KONEIJ),WORK(KXMAT))
         CALL DIAGEN(WORK(KONEIA),WORK(KT2AMT),WORK(KONEAI))
C
C------------------------------------------------------------
C        Get the contribution from CC3.
C------------------------------------------------------------
C
         IF (CC3) THEN
C
            IF (LWRK1 .LT. MAX(NT1AMX,NMATIJ(1),NMATAB(1))) THEN
               CALL QUIT('Out of memory in CC_D1AO (CC3)')
            ENDIF
C
            LUPTIA = -1
            CALL WOPEN2(LUPTIA,FNDPTIA,64,0)
            IOFF = 1 + NT1AMX*(ILLNR-1)
            CALL GETWA2(LUPTIA,FNDPTIA,WORK(KEND1),IOFF,NT1AMX)
            CALL WCLOSE2(LUPTIA,FNDPTIA,'KEEP')
            CALL DAXPY(NT1AMX,1.0D0,WORK(KEND1),1,WORK(KONEIA),1)
C
            LUPTAB = -1
            CALL WOPEN2(LUPTAB,FNDPTAB,64,0)
            IOFF = 1 + NMATAB(1)*(ILLNR-1)
            CALL GETWA2(LUPTAB,FNDPTAB,WORK(KEND1),IOFF,NMATAB(1))
            CALL WCLOSE2(LUPTAB,FNDPTAB,'KEEP')
            CALL DAXPY(NMATAB(1),1.0D0,WORK(KEND1),1,WORK(KONEAB),1)
C
            LUPTIJ = -1
            CALL WOPEN2(LUPTIJ,FNDPTIJ,64,0)
            IOFF = 1 + NMATIJ(1)*(ILLNR-1)
            CALL GETWA2(LUPTIJ,FNDPTIJ,WORK(KEND1),IOFF,NMATIJ(1))
            CALL WCLOSE2(LUPTIJ,FNDPTIJ,'KEEP')
            CALL DAXPY(NMATIJ(1),1.0D0,WORK(KEND1),1,WORK(KONEIJ),1)
C
            LUPTIA2 = -1
            CALL WOPEN2(LUPTIA2,FNDPTIA2,64,0)
            IOFF = 1 + NT1AMX*(ILLNR-1)
            CALL GETWA2(LUPTIA2,FNDPTIA2,WORK(KEND1),IOFF,NT1AMX)
            CALL WCLOSE2(LUPTIA2,FNDPTIA2,'KEEP')
            CALL DAXPY(NT1AMX,1.0D0,WORK(KEND1),1,WORK(KONEIA),1)
         ENDIF
C
C-----------------------------------------------------------
C        Backtransform the one electron density to AO-basis.
C-----------------------------------------------------------
C
         CALL DZERO(AODEN,N2BST(1))
C
         ISDEN = 1
         CALL CC_DENAO(AODEN,1,WORK(KONEAI),WORK(KONEAB),
     *                 WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *                 WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
         IF (RELORB) THEN
C
C----------------------------------------------------------------
C           Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------------
C
            LUSIFC = -1
            CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &                  .FALSE.)
            REWIND LUSIFC
            CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
            READ (LUSIFC)
            READ (LUSIFC)
            READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
            CALL GPCLOSE(LUSIFC,'KEEP')
C
            CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C-------------------------------------
C        Add orbital relaxation terms.
C-------------------------------------
C
c           CALL CC_D1ORRE(AODEN,ZKAM(NMATIJ(1)+NMATAB(1)+1),
c    *                     WORK(KEND1),LWRK1)
C
            KOFDIJ = 1
            KOFDAB = KOFDIJ + NMATIJ(1)
            KOFDAI = KOFDAB + NMATAB(1)
            KOFDIA = KOFDAI + NT1AMX
C
            ISDEN = 1
            CALL CC_DENAO(AODEN,1,ZKAM(KOFDAI),ZKAM(KOFDAB),
     *                    ZKAM(KOFDIJ),ZKAM(KOFDIA),ISDEN,WORK(KCMO),1,
     *                    WORK(KCMO),1,WORK(KEND1),LWRK1)
C
C
C-------------------------------------------------------
C           Add frozen core contributions to AO density.
C-------------------------------------------------------
C
            IF (FROIMP) THEN
C
               KOFFAI = NMATIJ(1) + NMATAB(1) + 2*NT1AMX + 1
               KOFFIA = KOFFAI    + NT1FRO(1)
               KOFFIJ = KOFFIA    + NT1FRO(1)
               KOFFJI = KOFFIJ    + NCOFRO(1)
C
               ISDEN = 1
               ICON  = 2
               CALL CC_D1FCB(AODEN,ZKAM(KOFFIJ),ZKAM(KOFFJI),
     *                       ZKAM(KOFFAI),ZKAM(KOFFIA),WORK(KEND1),
     *                       LWRK1,ISDEN,ICON)
C
            ENDIF
         ENDIF
C
      ENDIF
C
      IF (FROIMP) THEN
         CALL QEXIT('CC_D1AO')
         RETURN
      END IF
C
C-----------------------------------------------------------------
C     Calculate additional terms to the one electron density
C     needed for the evaluation of the natural occupation numbers.
C     Only if NO else skip.
C-----------------------------------------------------------------
C
      IF (.NOT. NO) THEN
         CALL QEXIT('CC_D1AO')
         RETURN
      END IF
C
      IF (.NOT. MP2) THEN
C
         CALL CCD1T1CO(WORK(KONEAI),WORK(KONEAB),WORK(KONEIJ),
     *                 WORK(KONEIA),WORK(KT1AM),WORK(KYMAT),
     *                 WORK(KEND1),LWRK1)
C
      ENDIF
C
C--------------------------------------------------
C     Reorder and diagonalize one electron density.
C--------------------------------------------------
C
      KDENFO = KEND1
      KEND2  = KDENFO + N2BST(1)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insufficient memory for second '//
     &        'allocation in CC_D1AO')
      ENDIF
C
      CALL DZERO(WORK(KDENFO),N2BST(1))
C
      CALL CCD1REOR(WORK(KDENFO),WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA))
C
      IF (IPRINT .GT. 50) THEN
C
         CALL AROUND('Nondiagonalised CC one electron density')
C
         DO 100 ISYM = 1,NSYM
C
            WRITE(LUPRI,333) 'Symmetry block number:', ISYM
  333       FORMAT(3X,A22,2X,I1)
C
            KOFFAI = KONEAI + IT1AM(ISYM,ISYM)
            KOFFAB = KONEAB + IMATAB(ISYM,ISYM)
            KOFFIJ = KONEIJ + IMATIJ(ISYM,ISYM)
            KOFFIA = KONEIA + IT1AM(ISYM,ISYM)
C
            CALL AROUND('DAI')
            CALL OUTPUT(WORK(KOFFAI),1,NVIR(ISYM),1,NRHF(ISYM),
     *                  NVIR(ISYM),NRHF(ISYM),1,LUPRI)
C
            CALL AROUND('DAB')
            CALL OUTPUT(WORK(KOFFAB),1,NVIR(ISYM),1,NVIR(ISYM),
     *                  NVIR(ISYM),NVIR(ISYM),1,LUPRI)
C
            CALL AROUND('DIJ')
            CALL OUTPUT(WORK(KOFFIJ),1,NRHF(ISYM),1,NRHF(ISYM),
     *                  NRHF(ISYM),NRHF(ISYM),1,LUPRI)
C
            CALL AROUND('DIA')
            CALL OUTPUT(WORK(KOFFIA),1,NVIR(ISYM),1,NRHF(ISYM),
     *                  NVIR(ISYM),NRHF(ISYM),1,LUPRI)
C
  100    CONTINUE
C
      ENDIF
C
      CALL CC_DEDIAN(WORK(KDENFO),MODEL,WORK(KEND2),LWRK2)
C
C-----------------------
C     write out timings.
C-----------------------
C
      TIMTOT = SECOND() -TIMTOT
C
      IF (IPRINT .GT. 9) THEN
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'CC one electron density '//
     &        'calculated & diagonalized'
         WRITE(LUPRI,*) 'Total time used in CC_D1AO:', TIMTOT
      ENDIF
C
      CALL QEXIT('CC_D1AO')
      RETURN
      END
C  /* Deck cc_denao */
      SUBROUTINE CC_DENAO(AODEN,ISYDAO,DENAI,DENAB,DENIJ,DENIA,ISYDMO,
     *                    XLAMDP,ISYMLP,XLAMDH,ISYMLH,WORK,LWORK)
C
C     Written by Asger Halkier 26/1 - 1996
C
C     Version: 1.0
C
C     Symmetries:   ISYDAO  --  AODEN
C                   ISYDMO  --  DENAI,DENAB,DENIJ,DENIA
C                   ISYMLP  --  XLAMDP
C                   ISYMLH  --  XLAMDH
C
C     Purpose: To transform all four blocks of the Coupled
C              Cluster one electron density to AO-basis and add
C              the results!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION AODEN(*), DENAI(*), DENAB(*), DENIJ(*), DENIA(*)
      DIMENSION XLAMDP(*), XLAMDH(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC_DENAO')
C
      DO ISYM1 = 1,NSYM
C
         ISYM2 = MULD2H(ISYDMO,ISYM1)
         ISYMA = MULD2H(ISYMLP,ISYM1)
         ISYMB = MULD2H(ISYMLH,ISYM2)
C         
         IF (ISYDAO.NE.MULD2H(ISYMA,ISYMB)) THEN
           CALL QUIT('Symmetry mismatch in CC_DENAO')
         END IF
C
         LVIR  = MAX(NBAS(ISYMA)*NVIR(ISYM2),NBAS(ISYMB)*NVIR(ISYM1))
         LRHF  = MAX(NBAS(ISYMA)*NRHF(ISYM2),NBAS(ISYMB)*NRHF(ISYM2))
         LNEED = MAX(LVIR,LRHF)
C
         IF (LWORK .LT. LNEED) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', LNEED
            CALL QUIT('Insufficient work space in CC_DENAO')
         ENDIF
C
         CALL DZERO(WORK,LNEED)
C
         NBASA  = MAX(NBAS(ISYMA),1)
         NBASB  = MAX(NBAS(ISYMB),1)
         NTOVI1 = MAX(NVIR(ISYM1),1)
         NTOVI2 = MAX(NVIR(ISYM2),1)
         NTORH1 = MAX(NRHF(ISYM1),1)
C
C-----------------------------------------
C        Transform virtual-occupied block.
C-----------------------------------------
C
         KOFF1  = IGLMVI(ISYMA,ISYM1) + 1
         KOFF2  =  IT1AM(ISYM1,ISYM2) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYM2),NVIR(ISYM1),
     *              ONE,XLAMDP(KOFF1),NBASA,DENAI(KOFF2),NTOVI1,
     *              ZERO,WORK,NBASA)
C
         KOFF3  = IGLMRH(ISYMB,ISYM2) + 1
         KOFF4  = IAODIS(ISYMA,ISYMB) + 1
C
         CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYM2),
     *              ONE,WORK,NBASA,XLAMDH(KOFF3),NBASB,ONE,
     *              AODEN(KOFF4),NBASA)
C
C----------------------------------------
C        Transform virtual-virtual block.
C----------------------------------------
C
         KOFF1  = IGLMVI(ISYMA,ISYM1) + 1
         KOFF2  = IMATAB(ISYM1,ISYM2) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMA),NVIR(ISYM2),NVIR(ISYM1),
     *              ONE,XLAMDP(KOFF1),NBASA,DENAB(KOFF2),NTOVI1,
     *              ZERO,WORK,NBASA)
C
         KOFF3  = IGLMVI(ISYMB,ISYM2) + 1
         KOFF4  = IAODIS(ISYMA,ISYMB) + 1
C
         CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NVIR(ISYM2),
     *              ONE,WORK,NBASA,XLAMDH(KOFF3),NBASB,ONE,
     *              AODEN(KOFF4),NBASA)
C
C------------------------------------------
C        Transform occupied-occupied block.
C------------------------------------------
C
         KOFF1  = IGLMRH(ISYMA,ISYM1) + 1
         KOFF2  = IMATIJ(ISYM1,ISYM2) + 1
C
C
         CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYM2),NRHF(ISYM1),
     *              ONE,XLAMDP(KOFF1),NBASA,DENIJ(KOFF2),NTORH1,
     *              ZERO,WORK,NBASA)
C
         KOFF3  = IGLMRH(ISYMB,ISYM2) + 1
         KOFF4  = IAODIS(ISYMA,ISYMB) + 1
C
         CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYM2),
     *              ONE,WORK,NBASA,XLAMDH(KOFF3),NBASB,ONE,
     *              AODEN(KOFF4),NBASA)
C
C-------------------------------------------------------------
C        Transform occupied-virtual block (stored transposed).
C-------------------------------------------------------------
C
         KOFF1  = IGLMVI(ISYMB,ISYM2) + 1
         KOFF2  =  IT1AM(ISYM2,ISYM1) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMB),NRHF(ISYM1),NVIR(ISYM2),
     *              ONE,XLAMDH(KOFF1),NBASB,DENIA(KOFF2),NTOVI2,
     *              ZERO,WORK,NBASB)
C
         KOFF3  = IGLMRH(ISYMA,ISYM1) + 1
         KOFF4  = IAODIS(ISYMA,ISYMB) + 1
C
         CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYM1),
     *              ONE,XLAMDP(KOFF3),NBASA,WORK,NBASB,ONE,
     *              AODEN(KOFF4),NBASA)
C
      END DO
C
      CALL QEXIT('CC_DENAO')
C
      RETURN
      END
C  /* Deck dijgen */
      SUBROUTINE DIJGEN(DIJ,XMAT)
C
C     Written by Asger Halkier 26/1 - 1996
C
C     Version: 1.0
C
C     Purpose: To set up the (occ,occ) part of the Coupled
C              Cluster one electron density.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION DIJ(*), XMAT(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
C--------------------------
C     Set up density block.
C--------------------------
C
      CALL QENTER('DIJGEN')
C
      CALL DAXPY(NMATIJ(1),-ONE,XMAT,1,DIJ,1)
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMJ = ISYMI
C
         DO 110 I = 1,NRHF(ISYMI)
C
            NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(I - 1) + I
C
            DIJ(NIJ) = DIJ(NIJ) + TWO
C
  110    CONTINUE
C
  100 CONTINUE
C
      CALL QEXIT('DIJGEN')
C
      RETURN
      END
C  /* Deck diagen */
      SUBROUTINE DIAGEN(DIA,T2AMT,Z1AM)
C
C     Written by Asger Halkier 26/1 - 1996
C
C     Version: 1.0
C
C     Purpose: To set up the (occ,vir) part of the Coupled
C              Cluster one electron density.
C
C     NOTE: This block is stored transposed, i.e. like a t1-amplitude!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION DIA(*), T2AMT(*), Z1AM(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('DIAGEN')
C
      ISYMDK = 1
      ISYMAI = ISYMDK
C
C--------------------------
C     Set up density block.
C--------------------------
C
      DO 100 NAI = 1,NT1AM(ISYMAI)
C
         DO 110 NDK = 1,NT1AM(ISYMDK)
C
            NDKAI = IT2AM(ISYMDK,ISYMAI) + INDEX(NDK,NAI)
C
            DIA(NAI) = DIA(NAI) + T2AMT(NDKAI)*Z1AM(NDK)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAGEN')
C
      RETURN
      END
C  /* Deck dijk01 */
      SUBROUTINE DIJK01(D2IJK,ISYIJK,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C     Written by Asger Halkier 29/1 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first contribution to D2IJK!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, FOUR = 4.0D0)
      DIMENSION D2IJK(*), XLAMDH(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJK01')
C
      ISYMK  = MULD2H(ISYMD,ISYMLH)
      ISYMIJ = MULD2H(ISYIJK,ISYMK)
C
      IF (ISYMIJ.NE.1) CALL QUIT('Symmetry mismatch in DIJK01')
C
C-----------------------------------------------
C     Calculate special adresses and add result.
C-----------------------------------------------
C
      DO 100 K = 1,NRHF(ISYMK)
C
         DO 110 ISYMJ = 1,NSYM
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               NIJK = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1)
     *              + IMATIJ(ISYMJ,ISYMJ) + NRHF(ISYMJ)*(J - 1) + J
               NDEL = IGLMRH(ISYMD,ISYMK) + NBAS(ISYMD)*(K - 1)
     *              + IDEL - IBAS(ISYMD)
C
               D2IJK(NIJK) = D2IJK(NIJK) + FOUR*XLAMDH(NDEL)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIJK01')
C
      RETURN
      END
C  /* Deck dijk02 */
      SUBROUTINE DIJK02(D2IJK,ISYIJK,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C     Written by Asger Halkier 29/1 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the second contribution to D2IJK!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2IJK(*), XLAMDH(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJK02')
C
      ISYMI  = MULD2H(ISYMD,ISYMLH)
C
      IF (ISYMI.NE.ISYIJK) CALL QUIT('Symmetry mismatch in DIJK02')
C
C-----------------------------------------------
C     Calculate special adresses and add result.
C-----------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMIJ = MULD2H(ISYMI,ISYMK)
C
         DO 110 K = 1,NRHF(ISYMK)
C
            DO 120 I = 1,NRHF(ISYMI)
C
               NIJK = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1)
     *              + IMATIJ(ISYMI,ISYMK) + NRHF(ISYMI)*(K - 1) + I
               NDEL = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I - 1)
     *              + IDEL - IBAS(ISYMD)
C
               D2IJK(NIJK) = D2IJK(NIJK) - TWO*XLAMDH(NDEL)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIJK02')
C
      RETURN
      END
C  /* Deck cc_d2pao */
      SUBROUTINE CC_D2PAO(D2AOIJ,D2AOAI,D2AOAB,D2AOIA,XLAMDP,D2IJK,
     *                    D2AIJ,D2IAJ,D2ABI,G,ISYMG,ISYDEN)
C
C     Written by Asger Halkier 30/1 - 1996
C
C     Version: 1.0
C
C     Purpose: To backtransform the four 2 electron densities
C              d(pq,k;del) to AO basis for a given gamma d(pq;gam;del)!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2AOIJ(*), D2AOAI(*), D2AOAB(*), D2AOIA(*), XLAMDP(*)
      DIMENSION D2IJK(*), D2AIJ(*), D2IAJ(*), D2ABI(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CC_D2PAO')
C
      ISYMK  = ISYMG
      ISYMPQ = MULD2H(ISYMG,ISYDEN)
C
      KOFFGK = IGLMRH(ISYMG,ISYMK) + G
C
C-------------------------------
C     Transform (occ,occ) block.
C-------------------------------
C
      KOFF1  = IMAIJK(ISYMPQ,ISYMK) + 1
C
      NTOTIJ = MAX(NMATIJ(ISYMPQ),1)
C
      CALL DGEMV('N',NMATIJ(ISYMPQ),NRHF(ISYMK),ONE,D2IJK(KOFF1),
     *           NTOTIJ,XLAMDP(KOFFGK),NBAS(ISYMG),ZERO,D2AOIJ,1)
C
C----------------------------------------------
C     Transform (vir,occ) and (occ,vir) blocks.
C----------------------------------------------
C
      KOFF2  = IT2BCD(ISYMPQ,ISYMK) + 1
C
      NTOTAI = MAX(NT1AM(ISYMPQ),1)
C
      CALL DGEMV('N',NT1AM(ISYMPQ),NRHF(ISYMK),ONE,D2AIJ(KOFF2),
     *           NTOTAI,XLAMDP(KOFFGK),NBAS(ISYMG),ZERO,D2AOAI,1)
C
      CALL DGEMV('N',NT1AM(ISYMPQ),NRHF(ISYMK),ONE,D2IAJ(KOFF2),
     *           NTOTAI,XLAMDP(KOFFGK),NBAS(ISYMG),ZERO,D2AOIA,1)
C
C-------------------------------
C     Transform (vir,vir) block.
C-------------------------------
C
      KOFF3  = ICKASR(ISYMPQ,ISYMK) + 1
C
      NTOTAB = MAX(NMATAB(ISYMPQ),1)
C
      CALL DGEMV('N',NMATAB(ISYMPQ),NRHF(ISYMK),ONE,D2ABI(KOFF3),
     *           NTOTAB,XLAMDP(KOFFGK),NBAS(ISYMG),ZERO,D2AOAB,1)
C
      CALL QEXIT('CC_D2PAO')
C
      RETURN
      END
C  /* Deck ccsd_d2en */
      SUBROUTINE CCSD_D2EN(ECCSD,XINT,XLAMDH,XLAMDP,D2AOIJ,D2AOAI,
     *                     D2AOAB,D2AOIA,G,ISYMG,ISYMD,WORK,LWORK)
C
C     Written by Asger Halkier 30/1 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the contribution from the 2 electron
C              Coupled Cluster density to the ccsd energy!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
      DIMENSION ECCSD(*), XINT(*), XLAMDH(*), XLAMDP(*)
      DIMENSION D2AOIJ(*), D2AOAI(*), D2AOAB(*), D2AOIA(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCSD_D2EN')
C
      ISYINT = ISYMD
      ISYDEN = ISYMD
      ISYMPQ = MULD2H(ISYMG,ISYDEN)
      ISALBE = ISYMPQ
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KAOINT = 1
      KEND1  = KAOINT + N2BST(ISALBE)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1. LT. 0) THEN
         WRITE(LUPRI,*) 'Needed:', KEND1, 'Available:', LWORK
         CALL QUIT('Insufficient memory for allocation in CCSD_D2EN')
      ENDIF
C
C-------------------------------------
C     Square up integral distribution.
C-------------------------------------
C
      KOFF1 = IDSAOG(ISYMG,ISYINT) + NNBST(ISALBE)*(G - 1) + 1
C
      CALL CCSD_SYMSQ(XINT(KOFF1),ISALBE,WORK(KAOINT))
C
C-------------------------------
C     Work space allocation two.
C-------------------------------
C
      LINTMO = MAX(NT1AM(ISYMPQ),NMATIJ(ISYMPQ),NMATAB(ISYMPQ))
C
      KINTMO = KEND1
      KEND2  = KINTMO + LINTMO
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2. LT. 0) THEN
         WRITE(LUPRI,*) 'Needed:', KEND2, 'Available:', LWORK
         CALL QUIT('Insufficient memory for allocation in CCSD_D2EN')
      ENDIF
C
C--------------------------------------
C     Calculate (occ,occ) contribution.
C--------------------------------------
C
      CALL DZERO(WORK(KINTMO),NMATIJ(ISYMPQ))
C
      DO 100 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMAL,ISALBE)
         ISYMJ  = ISYMBE
         ISYMI  = ISYMAL
C
         LNEED  = LWORK - KEND2 - NBAS(ISYMAL)*NRHF(ISYMJ)
C
         IF (LNEED .LT. 0) THEN
            WRITE(LUPRI,*) 'Work exceeded by:', -LNEED
            CALL QUIT('Insufficient work space in CCSD_D2EN')
         ENDIF
C
         CALL DZERO(WORK(KEND2),NBAS(ISYMAL)*NRHF(ISYMJ))
C
         KOFF2  = KAOINT + IAODIS(ISYMAL,ISYMBE)
         KOFF3  = IGLMRH(ISYMBE,ISYMJ) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NBAS(ISYMBE),ONE,
     *              WORK(KOFF2),NTOTAL,XLAMDH(KOFF3),NTOTBE,ZERO,
     *              WORK(KEND2),NTOTAL)
C
         KOFF4  = IGLMRH(ISYMAL,ISYMI) + 1
         KOFF5  = KINTMO + IMATIJ(ISYMI,ISYMJ)
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYMAL),HALF,
     *              XLAMDP(KOFF4),NTOTAL,WORK(KEND2),NTOTAL,ZERO,
     *              WORK(KOFF5),NTOTI)
C
  100 CONTINUE
C
      ECCSD(1) = ECCSD(1) + DDOT(NMATIJ(ISYMPQ),WORK(KINTMO),1,
     *                           D2AOIJ,1)
C
C--------------------------------------
C     Calculate (vir,vir) contribution.
C--------------------------------------
C
      CALL DZERO(WORK(KINTMO),NMATAB(ISYMPQ))
C
      DO 110 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMAL,ISALBE)
         ISYMA  = ISYMAL
         ISYMB  = ISYMBE
C
         LNEED  = LWORK - KEND2 - NBAS(ISYMAL)*NVIR(ISYMB)
C
         IF (LNEED .LT. 0) THEN
            WRITE(LUPRI,*) 'Work exceeded by:', -LNEED
            CALL QUIT('Insufficient work space in CCSD_D2EN')
         ENDIF
C
         CALL DZERO(WORK(KEND2),NBAS(ISYMAL)*NVIR(ISYMB))
C
         KOFF6  = KAOINT + IAODIS(ISYMAL,ISYMBE)
         KOFF7  = IGLMVI(ISYMBE,ISYMB) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NVIR(ISYMB),NBAS(ISYMBE),ONE,
     *              WORK(KOFF6),NTOTAL,XLAMDH(KOFF7),NTOTBE,ZERO,
     *              WORK(KEND2),NTOTAL)
C
         KOFF8  = IGLMVI(ISYMAL,ISYMA) + 1
         KOFF9  = KINTMO + IMATAB(ISYMA,ISYMB)
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NVIR(ISYMB),NBAS(ISYMAL),HALF,
     *              XLAMDP(KOFF8),NTOTAL,WORK(KEND2),NTOTAL,ZERO,
     *              WORK(KOFF9),NTOTA)
C
  110 CONTINUE
C
      ECCSD(1) = ECCSD(1) + DDOT(NMATAB(ISYMPQ),WORK(KINTMO),1,
     *                           D2AOAB,1)
C
C------------------------------------------
C     Calculate the (vir,occ) contribution.
C------------------------------------------
C
      CALL DZERO(WORK(KINTMO),NT1AM(ISYMPQ))
C
      DO 120 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMAL,ISALBE)
         ISYMI  = ISYMBE
         ISYMA  = ISYMAL
C
         LNEED  = LWORK - KEND2 - NBAS(ISYMAL)*NRHF(ISYMI)
C
         IF (LNEED .LT. 0) THEN
            WRITE(LUPRI,*) 'Work exceeded by:', -LNEED
            CALL QUIT('Insufficient work space in CCSD_D2EN')
         ENDIF
C
         CALL DZERO(WORK(KEND2),NBAS(ISYMAL)*NRHF(ISYMI))
C
         KOFF10 = KAOINT + IAODIS(ISYMAL,ISYMBE)
         KOFF11 = IGLMRH(ISYMBE,ISYMI) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMBE),ONE,
     *              WORK(KOFF10),NTOTAL,XLAMDH(KOFF11),NTOTBE,ZERO,
     *              WORK(KEND2),NTOTAL)
C
         KOFF12 = IGLMVI(ISYMAL,ISYMA) + 1
         KOFF13 = KINTMO + IT1AM(ISYMA,ISYMI)
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMAL),HALF,
     *              XLAMDP(KOFF12),NTOTAL,WORK(KEND2),NTOTAL,ZERO,
     *              WORK(KOFF13),NTOTA)
C
  120 CONTINUE
C
      ECCSD(1) = ECCSD(1) + DDOT(NT1AM(ISYMPQ),WORK(KINTMO),1,
     *                           D2AOAI,1)
C
C------------------------------------------
C     Calculate the (occ,vir) contribution.
C------------------------------------------
C
      CALL DZERO(WORK(KINTMO),NT1AM(ISYMPQ))
C
      DO 130 ISYMAL = 1,NSYM
C
         ISYMBE = MULD2H(ISYMAL,ISALBE)
         ISYMI  = ISYMAL
         ISYMA  = ISYMBE
C
         LNEED = LWORK - KEND2 - NBAS(ISYMBE)*NRHF(ISYMI)
C
         IF (LNEED .LT. 0) THEN
            WRITE(LUPRI,*) 'Work exceeded by:', -LNEED
            CALL QUIT('Insufficient work space in CCSD_D2EN')
         ENDIF
C
         CALL DZERO(WORK(KEND2),NBAS(ISYMBE)*NRHF(ISYMI))
C
         KOFF14 = KAOINT + IAODIS(ISYMAL,ISYMBE)
         KOFF15 = IGLMRH(ISYMAL,ISYMI) + 1
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         CALL DGEMM('T','N',NBAS(ISYMBE),NRHF(ISYMI),NBAS(ISYMAL),ONE,
     *              WORK(KOFF14),NTOTAL,XLAMDP(KOFF15),NTOTAL,ZERO,
     *              WORK(KEND2),NTOTBE)
C
         KOFF16 = IGLMVI(ISYMBE,ISYMA) + 1
         KOFF17 = KINTMO + IT1AM(ISYMA,ISYMI)
C
         NTOTBE = MAX(NBAS(ISYMBE),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NBAS(ISYMBE),HALF,
     *              XLAMDH(KOFF16),NTOTBE,WORK(KEND2),NTOTBE,ZERO,
     *              WORK(KOFF17),NTOTA)
C
  130 CONTINUE
C
      ECCSD(1) = ECCSD(1) + DDOT(NT1AM(ISYMPQ),WORK(KINTMO),1,
     *                           D2AOIA,1)
C
      CALL QEXIT('CCSD_D2EN')
C
      RETURN
      END
C  /* Deck dijk11 */
      SUBROUTINE DIJK11(D2IJK,ISYIJK,XLAMDH,ISYMLH,DIA,ISYD1,
     *                  WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 2/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first two contribution to D2IJK
C              from projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2IJK(*), XLAMDH(*), DIA(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJK11')
C
      ISYMA  = MULD2H(ISYMD,ISYMLH)
      ISYMK  = MULD2H(ISYMA,ISYD1)
      ISYMIJ = MULD2H(ISYMK,ISYIJK)
C
      IF (ISYMK.NE.ISYIJK) CALL QUIT('Symmetry mismatch in DIJK11 (1)')
C
      IF (LWORK .LT. NRHF(ISYMK)) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'NEEDED:', NRHF(ISYMK)
         CALL QUIT('Insufficient space for allocation in DIJK11')
      ENDIF
C
      CALL DZERO(WORK,NRHF(ISYMK))
C
C-------------------------------------------
C     Contract D(ia) with lamda hole matrix.
C-------------------------------------------
C
      KOFF1 = IT1AM(ISYMA,ISYMK) + 1
      KOFF2 = IGLMVI(ISYMD,ISYMA) + IDEL - IBAS(ISYMD)
C
      NTOTA = MAX(NVIR(ISYMA),1)
C
      CALL DGEMV('T',NVIR(ISYMA),NRHF(ISYMK),ONE,DIA(KOFF1),NTOTA,
     *           XLAMDH(KOFF2),NBAS(ISYMD),ZERO,WORK,1)
C
C---------------------------------------
C     Calculate the DIJK11 contribution.
C---------------------------------------
C
      DO 100 K = 1,NRHF(ISYMK)
C
         DO 110 ISYMI = 1,NSYM
C
            ISYMJ = ISYMI
C
            DO 120 I = 1,NRHF(ISYMI)
C
               J    = I
               NIJK = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1)
     *              + IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
C
               D2IJK(NIJK) = D2IJK(NIJK) + TWO*WORK(K)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C---------------------------------------
C     Calculate the DIJK12 contribution.
C---------------------------------------
C
      ISYMI  = MULD2H(ISYMA,ISYD1)
      ISYMJK = MULD2H(ISYMI,ISYIJK)

      IF (ISYMI.NE.ISYIJK) CALL QUIT('Symmetry mismatch in DIJK11 (2)')
C
      DO 130 ISYMK = 1,NSYM
C
         ISYMJ  = ISYMK
         ISYMIJ = MULD2H(ISYMI,ISYMJ)
C
         DO 140 K = 1,NRHF(ISYMK)
C
            J = K
C
            DO 150 I = 1,NRHF(ISYMI)
C
               NIJK = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1)
     *              + IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
C
               D2IJK(NIJK) = D2IJK(NIJK) - WORK(I)
C
  150       CONTINUE
  140    CONTINUE
  130 CONTINUE
C
      CALL QEXIT('DIJK11')
C
      RETURN
      END
C  /* Deck dijk13 */
      SUBROUTINE DIJK13(D2IJK,ISYIJK,DENAI,ISYD1,T2TIN,ISYTIN)
C
C     Written by Asger Halkier 2/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the third contribution to D2IJK
C              from projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IJK(*), DENAI(*), T2TIN(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJK13')
C
      IF (MULD2H(ISYTIN,ISYD1).NE.ISYIJK)
     &  CALL QUIT('Symmetry mismatch in DIJK13')
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMEI = MULD2H(ISYMK,ISYTIN)
         ISYMIJ = MULD2H(ISYMK,ISYIJK)
C
         DO 110 K = 1,NRHF(ISYMK)
C
            DO 120 ISYMI = 1,NSYM
C
               ISYME = MULD2H(ISYMI,ISYMEI)
               ISYMJ = MULD2H(ISYME,ISYD1)
C
               KOFF1 = IT2BCD(ISYMEI,ISYMK) + NT1AM(ISYMEI)*(K - 1)
     *               + IT1AM(ISYME,ISYMI) + 1
               KOFF2 = IT1AM(ISYME,ISYMJ) + 1
               KOFF3 = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1)
     *               + IMATIJ(ISYMI,ISYMJ) + 1
C
               NTOTE = MAX(NVIR(ISYME),1)
               NTOTI = MAX(NRHF(ISYMI),1)
C
               CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYME),
     *                    -ONE,T2TIN(KOFF1),NTOTE,DENAI(KOFF2),NTOTE,
     *                    ONE,D2IJK(KOFF3),NTOTI)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIJK13')
C
      RETURN
      END
C  /* Deck daij11 */
      SUBROUTINE DAIJ11(D2AIJ,ISYAIJ,DAI,ISYD1,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C     Written by Asger Halkier 2/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first contribution to D2AIJ
C              from projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2AIJ(*), XLAMDH(*), DAI(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DAIJ11')
C
      ISYMJ  = MULD2H(ISYMD,ISYMLH)
      ISYMAI = ISYD1
      IF (MULD2H(ISYMJ,ISYMAI).NE.ISYAIJ) THEN
        CALL QUIT('Symmetry mismatch in DAIJ11')
      END IF
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      DO 100 J = 1,NRHF(ISYMJ)
C
         NDJ = IGLMRH(ISYMD,ISYMJ) + NBAS(ISYMD)*(J - 1) + 
     *         IDEL - IBAS(ISYMD)
C
         FACT = TWO*XLAMDH(NDJ)
C
         KOFF1 = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1) + 1
C
         CALL DAXPY(NT1AM(ISYMAI),FACT,DAI,1,D2AIJ(KOFF1),1)
C
  100 CONTINUE
C
      CALL QEXIT('DAIJ11')
C
      RETURN
      END
C  /* Deck daij12 */
      SUBROUTINE DAIJ12(D2AIJ,ISYAIJ,DAI,ISYD1,XLAMDH,ISYMLH,
     *                  WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 2/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the second contribution to D2AIJ
C              from projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2AIJ(*), XLAMDH(*), DAI(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DAIJ12')
C
      ISYMK  = MULD2H(ISYMD,ISYMLH)
      ISYMA  = MULD2H(ISYMK,ISYD1)
C
      IF (ISYMA.NE.ISYAIJ) THEN
        CALL QUIT('Symmetry mismatch in DAIJ12')
      END IF
C
      IF (LWORK .LT. NVIR(ISYMA)) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'NEEDED:', NVIR(ISYMA)
         CALL QUIT('Insufficient space for allocation in DAIJ12')
      ENDIF
C
      CALL DZERO(WORK,NVIR(ISYMA))
C
C----------------------------------------------------------
C     Calculate contraction of D(ai) and lamda hole matrix.
C----------------------------------------------------------
C
      KOFF1 = IT1AM(ISYMA,ISYMK) + 1
      KOFF2 = IGLMRH(ISYMD,ISYMK) + IDEl - IBAS(ISYMD)
C
      NTOTA = MAX(NVIR(ISYMA),1)
C
      CALL DGEMV('N',NVIR(ISYMA),NRHF(ISYMK),ONE,DAI(KOFF1),NTOTA,
     *           XLAMDH(KOFF2),NBAS(ISYMD),ZERO,WORK,1)
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMJ  = ISYMI
         ISYMAI = MULD2H(ISYMI,ISYMA)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            J    = I
            NAIJ = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     *           + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
            CALL DAXPY(NVIR(ISYMA),-ONE,WORK,1,D2AIJ(NAIJ),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DAIJ12')
C
      RETURN
      END
C  /* Deck diaj11 */
      SUBROUTINE DIAJ11(D2IAJ,ISYIAJ,DIA,ISYD1,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C     Written by Asger Halkier 4/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first and second contribution to D2IAJ
C              from projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2IAJ(*), XLAMDH(*), DIA(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIAJ11')
C
      ISYMJ  = MULD2H(ISYMD,ISYMLH)
      ISYMAI = ISYD1

      IF (MULD2H(ISYMJ,ISYMAI).NE.ISYIAJ) 
     *   CALL QUIT('Symmetry mismatch in DIAJ11. (1)')
C
C--------------------------------------
C     Calculate the first contribution.
C--------------------------------------
C
      DO 100 J = 1,NRHF(ISYMJ)
C
         NDJ = IGLMRH(ISYMD,ISYMJ) + NBAS(ISYMD)*(J - 1) + 
     *         IDEL - IBAS(ISYMD)
C
         FACT = TWO*XLAMDH(NDJ)
C
         KOFF = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1) + 1
C
         CALL DAXPY(NT1AM(ISYMAI),FACT,DIA,1,D2IAJ(KOFF),1)
C
  100 CONTINUE
C
      ISYMI  = MULD2H(ISYMD,ISYMLH)
      ISYMAJ = ISYD1

      IF (MULD2H(ISYMI,ISYMAJ).NE.ISYIAJ) 
     *   CALL QUIT('Symmetry mismatch in DIAJ11. (2)')
C
C---------------------------------------
C     Calculate the second contribution.
C---------------------------------------
C
      DO 110 ISYMJ = 1,NSYM
C
         ISYMA  = MULD2H(ISYMJ,ISYMAJ)
         ISYMAI = MULD2H(ISYMA,ISYMI)
C
         DO 120 J = 1,NRHF(ISYMJ)
C
            DO 130 I = 1,NRHF(ISYMI)
C
               NDI = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I- 1) + IDEL 
     *             - IBAS(ISYMD)
               NAJ  = IT1AM(ISYMA,ISYMJ) + NVIR(ISYMA)*(J - 1) + 1
               NIAJ = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     *              + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I - 1) + 1
C
               CALL DAXPY(NVIR(ISYMA),-XLAMDH(NDI),DIA(NAJ),1,
     *                    D2IAJ(NIAJ),1)
C
  130       CONTINUE
  120    CONTINUE
  110 CONTINUE
C
      CALL QEXIT('DIAJ11')
C
      RETURN
      END
C  /* Deck dabi11 */
      SUBROUTINE DABI11(D2ABI,ISYABI,DAI,ISYD1,T2TDEL,ISYMTI)
C
C     Written by Asger Halkier 4/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the one and only contribution to D2ABI
C              from projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2ABI(*), DAI(*), T2TDEL(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DABI11')
C
      IF (MULD2H(ISYMTI,ISYD1).NE.ISYABI)
     &  CALL QUIT('Symmetry in DABI11')
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMBM = MULD2H(ISYMI,ISYMTI)
         ISYMAB = MULD2H(ISYMI,ISYABI)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            DO 120 ISYMA = 1,NSYM
C
               ISYMM = MULD2H(ISYMA,ISYD1)
               ISYMB = MULD2H(ISYMM,ISYMBM)
C
               KOFF1 = IT1AM(ISYMA,ISYMM)   + 1
               KOFF2 = IT2BCD(ISYMBM,ISYMI) + NT1AM(ISYMBM)*(I - 1)
     *               + IT1AM(ISYMB,ISYMM)   + 1
               KOFF3 = ICKASR(ISYMAB,ISYMI) + NMATAB(ISYMAB)*(I - 1)
     *               + IMATAB(ISYMA,ISYMB)  + 1
C
               NTOTA = MAX(NVIR(ISYMA),1)
               NTOTB = MAX(NVIR(ISYMB),1)
C
               CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMM),
     *                    ONE,DAI(KOFF1),NTOTA,T2TDEL(KOFF2),NTOTB,ONE,
     *                    D2ABI(KOFF3),NTOTA)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DABI11')
C
      RETURN
      END
C  /* Deck dijg11 */
      SUBROUTINE DIJG11(D2AOIJ,XLAMDH,Z1AO,IDEL,ISYMD,G,ISYMG)
C
C     Written by Asger Halkier 4/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first and second contribution to D2AOIJ
C              from projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2AOIJ(*), XLAMDH(*), Z1AO(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJG11')
C
C--------------------------------------
C     Calculate the first contribution.
C--------------------------------------
C
      IF (ISYMG .EQ. ISYMD) THEN
C
         ISYMK = ISYMG
C
         NDK = IGLMRH(ISYMD,ISYMK) + IDEL - IBAS(ISYMD)
         NGK = IGLMRH(ISYMG,ISYMK) + G
C
         FACT = DDOT(NRHF(ISYMK),Z1AO(NGK),NBAS(ISYMG),XLAMDH(NDK),
     *               NBAS(ISYMD))
C
         DO 100 ISYMI = 1,NSYM
C
            ISYMJ = ISYMI
C  
            DO 110 I = 1,NRHF(ISYMI)
C
               J   = I
               NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
C
               D2AOIJ(NIJ) = D2AOIJ(NIJ) + TWO*FACT
C
  110       CONTINUE
  100    CONTINUE
C
      ENDIF
C
C-----------------------------------
C     Calculate second contribution.
C-----------------------------------
C
      ISYMI = ISYMD
      ISYMJ = ISYMG
C
      DO 120 J = 1,NRHF(ISYMJ)
C
         NGJ = IGLMRH(ISYMG,ISYMJ) + NBAS(ISYMG)*(J - 1) + G
C
         DO 130 I = 1,NRHF(ISYMI)
C
            NDI = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I - 1) + IDEL 
     *          - IBAS(ISYMD)
            NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
C
            D2AOIJ(NIJ) = D2AOIJ(NIJ) - XLAMDH(NDI)*Z1AO(NGJ)
C
  130    CONTINUE
  120 CONTINUE
C
      CALL QEXIT('DIJG11')
C
      RETURN
      END
C  /* Deck diaj13 */
      SUBROUTINE DIAJ13(D2IAJ,ISYIAJ,DAI,ISYD1,T2AMT,XLAMDH,ISYMLH,
     *                  WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 5/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the third contribution to D2IAJ from 
C              projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAJ(*), DAI(*), T2AMT(*), XLAMDH(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('DIAJ13')
C
      IF (MULD2H(ISYMLH,ISYMD).NE.MULD2H(ISYIAJ,ISYD1)) THEN
        CALL QUIT('Symmetry mismatch in DIAJ13.')
      END IF
C
      ISYMK  = MULD2H(ISYMD,ISYMLH)
      ISYME  = MULD2H(ISYMK,ISYD1)
C
C-----------------------------------------
C     Spaghetti goto's if no contribution.
C-----------------------------------------
C
      IF (NRHF(ISYMK) .EQ. 0) GOTO 101
      IF (NVIR(ISYME) .EQ. 0) GOTO 101
      IF (NBAS(ISYMD) .EQ. 0) GOTO 101
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KVECE = 1
      KEND1 = KVECE + NVIR(ISYME)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in DIAJ13')
      ENDIF
C
      CALL DZERO(WORK(KVECE),NVIR(ISYME))
C
C-----------------------------------------------
C     Calculate contraction of zeta1 and xlamdh.
C-----------------------------------------------
C
      KOFF1 = IT1AM(ISYME,ISYMK) + 1
      KOFF2 = IGLMRH(ISYMD,ISYMK) + IDEL - IBAS(ISYMD)
C
      NTOTE = MAX(NVIR(ISYME),1)
C
      CALL DGEMV('N',NVIR(ISYME),NRHF(ISYMK),ONE,DAI(KOFF1),NTOTE,
     *           XLAMDH(KOFF2),NBAS(ISYMD),ZERO,WORK(KVECE),1)
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMJ,ISYIAJ)
         ISYMEJ = ISYMAI
C
C----------------------------------
C        Work space allocation two.
C----------------------------------
C
         KT2SUB = KEND1
         KEND2  = KT2SUB + NT1AM(ISYMAI)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
            CALL QUIT('Insufficient work space for '//
     &           'allocation in DIAJ13')
         ENDIF
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 E = 1,NVIR(ISYME)
C
               NEJ = IT1AM(ISYME,ISYMJ) + NVIR(ISYME)*(J - 1) + E
C
C--------------------------------------------------
C              Copy t1-submatrix (ai) out of t2amt.
C--------------------------------------------------
C
               DO 130 NAI = 1,NT1AM(ISYMAI)
C
                  NAIEJ = IT2AM(ISYMAI,ISYMEJ) + INDEX(NAI,NEJ)
C
                  WORK(KT2SUB + NAI - 1) = T2AMT(NAIEJ)
C
  130          CONTINUE
C
C------------------------------------------------------
C              Final contraction and storage in result.
C------------------------------------------------------
C
               FACT  = - WORK(KVECE + E - 1)
C
               KOFF3 = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1) + 1
C
               CALL DAXPY(NT1AM(ISYMAI),FACT,WORK(KT2SUB),1,
     *                    D2IAJ(KOFF3),1)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
  101 CONTINUE
C
      CALL QEXIT('DIAJ13')
C
      RETURN
      END
C  /* Deck diag11 */
      SUBROUTINE DIAG11(D2AOIA,Z1AO,T2TIN,ISYMD,G,ISYMG)
C
C     Written by Asger Halkier 7/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the one and only contribution to D2AOIA
C              from projection against singles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2AOIA(*), Z1AO(*), T2TIN(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIAG11')
C
      ISYAIM = ISYMD
      ISYMM  = ISYMG
      ISYMAI = MULD2H(ISYMM,ISYAIM)
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      KOFF1  = IT2BCD(ISYMAI,ISYMM) + 1
      KOFF2  = IGLMRH(ISYMG,ISYMM) + G
C
      NTOTAI = MAX(NT1AM(ISYMAI),1)
C
      CALL DGEMV('N',NT1AM(ISYMAI),NRHF(ISYMM),ONE,T2TIN(KOFF1),
     *           NTOTAI,Z1AO(KOFF2),NBAS(ISYMG),ONE,D2AOIA,1)
C
      CALL QEXIT('DIAG11')
C
      RETURN
      END
C  /* Deck dijk25 */
      SUBROUTINE DIJK25(D2IJK,ISYIJK,XMINT,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C     Written by Asger Halkier 10/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the fifth contribution to D2IJK
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IJK(*), XMINT(*), XLAMDH(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJK25')
C
      ISYML  = MULD2H(ISYMD,ISYMLH)
      IF (ISYIJK.NE.ISYML) CALL QUIT('Symmetry mismatch in DIJK25')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      KOFF1  = I3ORHF(ISYIJK,ISYML) + 1
      KOFF2  = IGLMRH(ISYMD,ISYML) + IDEL - IBAS(ISYMD)
C
      NTOIJK = MAX(NMAIJK(ISYIJK),1)
C
      CALL DGEMV('N',NMAIJK(ISYIJK),NRHF(ISYML),ONE,XMINT(KOFF1),
     *           NTOIJK,XLAMDH(KOFF2),NBAS(ISYMD),ONE,D2IJK,1)
C
      CALL QEXIT('DIJK25')
C
      RETURN
      END
C  /* Deck dijk24 */
      SUBROUTINE DIJK24(D2IJK,ISYIJK,XMAT,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C     Written by Asger Halkier 10/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the fourth contribution to D2IJK
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IJK(*), XMAT(*), XLAMDH(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJK24')
C
      ISYMI  = MULD2H(ISYMD,ISYMLH)
      ISYMKJ = 1
C
      IF (ISYMI.NE.ISYIJK) CALL QUIT('Symmetry mismatch in DIJK24')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMJ  = MULD2H(ISYMK,ISYMKJ)
         ISYMIJ = MULD2H(ISYMI,ISYMJ)
C
         DO 110 K = 1,NRHF(ISYMK)
C
            DO 120 J = 1,NRHF(ISYMJ)
C
               NKJ  = IMATIJ(ISYMK,ISYMJ)  + NRHF(ISYMK)*(J - 1) + K
               NDEI = IGLMRH(ISYMD,ISYMI) + IDEL - IBAS(ISYMD)
               NIJK = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1)
     *              + IMATIJ(ISYMI,ISYMJ)  + NRHF(ISYMI)*(J - 1) + 1
C
               CALL DAXPY(NRHF(ISYMI),XMAT(NKJ),XLAMDH(NDEI),
     *                    NBAS(ISYMD),D2IJK(NIJK),1)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIJK24')
C
      RETURN
      END
C  /* Deck dijk23 */
      SUBROUTINE DIJK23(D2IJK,ISYIJK,XMAT,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C     Written by Asger Halkier 10/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the third contribution to D2IJK
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2IJK(*), XMAT(*), XLAMDH(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJK23')
C
      ISYMK  = MULD2H(ISYMD,ISYMLH)
      ISYMIJ = 1
C
      IF (ISYMK.NE.ISYIJK) CALL QUIT('Symmetry mismatch in DIJK23')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 K = 1,NRHF(ISYMK)
C
         NDEK = IGLMRH(ISYMD,ISYMK) + NBAS(ISYMD)*(K - 1)
     *        + IDEL - IBAS(ISYMD)
         NIJK = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1) + 1
C
         FACT = -TWO*XLAMDH(NDEK)
C
         CALL DAXPY(NMATIJ(ISYMIJ),FACT,XMAT,1,D2IJK(NIJK),1)
C
  100 CONTINUE
C
      CALL QEXIT('DIJK23')
C
      RETURN
      END
C  /* Deck dijk21 */
      SUBROUTINE DIJK21(D2IJK,ISYIJK,XMAT,XLAMDH,ISYMLH,
     *                  WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 10/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first and second contribution to D2IJK
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2IJK(*), XMAT(*), XLAMDH(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIJK21')
C
      ISYML  = MULD2H(ISYMD,ISYMLH)
      ISYMX1 = ISYML
C
      IF (ISYML.NE.ISYIJK) CALL QUIT('Symmetry mismatch in DIJK21')
C
C-----------------------------------------------------------------
C     Calculate contraction of X-intermediate & lamda hole matrix.
C-----------------------------------------------------------------
C
      IF (LWORK .LT. NRHF(ISYMX1)) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', NRHF(ISYMX1)
         CALL QUIT('Insufficient work space in DIJK21')
      ENDIF
C
      CALL DZERO(WORK,NRHF(ISYMX1))
C
      KOFF1 = IMATIJ(ISYMX1,ISYML) + 1
      KOFF2 = IGLMRH(ISYMD,ISYML) + IDEL - IBAS(ISYMD)
C
      NTOTI = MAX(NRHF(ISYMX1),1)
C
      CALL DGEMV('N',NRHF(ISYMX1),NRHF(ISYML),ONE,XMAT(KOFF1),NTOTI,
     *           XLAMDH(KOFF2),NBAS(ISYMD),ZERO,WORK,1)
C
C----------------------------------
C     Calculate first contribution.
C----------------------------------
C
      ISYMI = ISYMX1
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMJ  = ISYMK
         ISYMIJ = MULD2H(ISYMJ,ISYMI)
C
         DO 110 K = 1,NRHF(ISYMK)
C
            J    = K
            NIJK = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1)
     *           + IMATIJ(ISYMI,ISYMJ)  + NRHF(ISYMI)*(J - 1) + 1
C
            CALL DAXPY(NRHF(ISYMI),ONE,WORK,1,D2IJK(NIJK),1)
C
  110    CONTINUE
  100 CONTINUE
C
C-----------------------------------
C     Calculate second contribution.
C-----------------------------------
C
      ISYMK  = ISYMX1
      ISYMIJ = 1
C
      DO 120 K = 1,NRHF(ISYMK)
C
         DO 130 ISYMJ = 1,NSYM
C
            ISYMI = ISYMJ
C
            DO 140 J = 1,NRHF(ISYMJ)
C
               I    = J
               NIJK = IMAIJK(ISYMIJ,ISYMK) + NMATIJ(ISYMIJ)*(K - 1)
     *              + IMATIJ(ISYMI,ISYMJ)  + NRHF(ISYMI)*(J - 1) + I
C
               D2IJK(NIJK) = D2IJK(NIJK) - TWO*WORK(K)
C
  140       CONTINUE
  130    CONTINUE
  120 CONTINUE
C
      CALL QEXIT('DIJK21')
C
      RETURN
      END
C  /* Deck daij22 */
      SUBROUTINE DAIJ22(D2AIJ,ISYAIJ,DAB,ISYD1,XLAMDH,ISYMLH,
     *                  WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 10/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the second contribution to D2AIJ
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2AIJ(*), DAB(*), XLAMDH(*), WORK(LWORK)
#include <priunit.h>      
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DAIJ22')
C
      ISYMB  = MULD2H(ISYMD,ISYMLH)
      ISYMA  = MULD2H(ISYMB,ISYD1)
      ISYMIJ = 1
C
      IF (ISYMA.NE.ISYAIJ) CALL QUIT('Symmetry mismatch in DAIJ22')
C
C-----------------------------------------------------------------
C     Calculate contraction of Y-intermediate & lamda hole matrix.
C-----------------------------------------------------------------
C
      IF (LWORK .LT. NVIR(ISYMA)) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', NVIR(ISYMA)
         CALL QUIT('Insufficient work space in DAIJ22')
      ENDIF
C
      CALL DZERO(WORK,NVIR(ISYMA))
C
      KOFF1 = IMATAB(ISYMA,ISYMB) + 1
      KOFF2 = IGLMVI(ISYMD,ISYMB) + IDEL - IBAS(ISYMD)
C
      NTOTA = MAX(NVIR(ISYMA),1)
C
      CALL DGEMV('N',NVIR(ISYMA),NVIR(ISYMB),ONE,DAB(KOFF1),NTOTA,
     *           XLAMDH(KOFF2),NBAS(ISYMD),ZERO,WORK,1)
C
C-------------------------------------
C     Calculate contribution to D2AIJ.
C-------------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMI  = ISYMJ
         ISYMAI = MULD2H(ISYMI,ISYMA)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            I    = J
            NAIJ = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     *           + IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I - 1) + 1
C
            CALL DAXPY(NVIR(ISYMA),-ONE,WORK,1,D2AIJ(NAIJ),1)
C
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DAIJ22')
C
      RETURN
      END
C  /* Deck dabi21 */
      SUBROUTINE DABI21(D2ABI,ISYABI,DAB,ISYD1,XLAMDH,ISYMLH,IDEL,ISYMD)
C
C     Written by Asger Halkier 11/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first contribution to D2ABI
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION D2ABI(*), DAB(*), XLAMDH(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DABI21')
C
      ISYMI  = MULD2H(ISYMD,ISYMLH)
      ISYMAB = ISYD1
C
      IF (MULD2H(ISYMI,ISYMAB).NE.ISYABI)
     *  CALL QUIT('Symmetry mismatch in DABI21')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 I = 1,NRHF(ISYMI)
C
         NDEI = IGLMRH(ISYMD,ISYMI) + NBAS(ISYMD)*(I - 1) 
     *        + IDEL - IBAS(ISYMD)
         NABI = ICKASR(ISYMAB,ISYMI) + NMATAB(ISYMAB)*(I - 1) + 1
C
         FACT = TWO*XLAMDH(NDEI)
C
         CALL DAXPY(NMATAB(ISYMAB),FACT,DAB,1,D2ABI(NABI),1)
C
  100 CONTINUE
C
      CALL QEXIT('DABI21')
C
      RETURN
      END
C  /* Deck diaj21 */
      SUBROUTINE DIAJ21(D2IAJ,ISYAIJ,YMAT,T2TINT,ISYMTI)
C
C     Written by Asger Halkier 11/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the first contribution to D2IAJ
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAJ(*), YMAT(*), T2TINT(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIAJ21')
C
      ISYMAE = 1
C
      IF (ISYMTI.NE.ISYAIJ) CALL QUIT('Symmetry mismatch in DIAJ21.')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMEI = MULD2H(ISYMJ,ISYMTI)
         ISYMAI = MULD2H(ISYMJ,ISYAIJ)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMA = 1,NSYM
C
               ISYME = MULD2H(ISYMA,ISYMAE)
               ISYMI = MULD2H(ISYMA,ISYMAI)
C
               KOFF1 = IMATAB(ISYMA,ISYME)  + 1
               KOFF2 = IT2BCD(ISYMEI,ISYMJ) + NT1AM(ISYMEI)*(J - 1)
     *               + IT1AM(ISYME,ISYMI)   + 1
               KOFF3 = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI)   + 1
C
               NTOTA = MAX(NVIR(ISYMA),1)
               NTOTE = MAX(NVIR(ISYME),1)
C
               CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYME),
     *                    -ONE,YMAT(KOFF1),NTOTA,T2TINT(KOFF2),NTOTE,
     *                    ONE,D2IAJ(KOFF3),NTOTA)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAJ21')
C
      RETURN
      END
C  /* Deck diaj22 */
      SUBROUTINE DIAJ22(D2IAJ,ISYAIJ,XMAT,T2TINT,ISYMTI)
C
C     Written by Asger Halkier 11/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the second contribution to D2IAJ
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAJ(*), XMAT(*), T2TINT(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIAJ22')
C
      IF (ISYAIJ.NE.ISYMTI) CALL QUIT('Symmetry mismatch in DIAJ22.')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMAM = MULD2H(ISYMJ,ISYMTI)
         ISYMAI = MULD2H(ISYMJ,ISYAIJ)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 ISYMA = 1,NSYM
C
               ISYMM = MULD2H(ISYMA,ISYMAM)
               ISYMI = MULD2H(ISYMA,ISYMAI)
C
               KOFF1 = IT2BCD(ISYMAM,ISYMJ) + NT1AM(ISYMAM)*(J - 1)
     *               + IT1AM(ISYMA,ISYMM)   + 1
               KOFF2 = IMATIJ(ISYMI,ISYMM)  + 1
               KOFF3 = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI)   + 1
C
               NTOTA = MAX(NVIR(ISYMA),1)
               NTOTI = MAX(NRHF(ISYMI),1)
C
               CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMM),
     *                    -ONE,T2TINT(KOFF1),NTOTA,XMAT(KOFF2),NTOTI,
     *                    ONE,D2IAJ(KOFF3),NTOTA)
C
  120       CONTINUE 
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAJ22')
C
      RETURN
      END
C  /* Deck diaj23 */
      SUBROUTINE DIAJ23(D2IAJ,ISYAIJ,XMAT,T2TINT,ISYMTI)
C
C     Written by Asger Halkier 11/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the third contribution to D2IAJ
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAJ(*), XMAT(*), T2TINT(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('DIAJ23')
C
      ISYMJM = 1
C
      IF (ISYMTI.NE.ISYAIJ) CALL QUIT('Symmetry mismatch in DIAJ23')
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 ISYMAI = 1,NSYM
C
         ISYMM = MULD2H(ISYMAI,ISYMTI)
         ISYMJ = MULD2H(ISYMAI,ISYAIJ)
C
         KOFF1 = IT2BCD(ISYMAI,ISYMM) + 1
         KOFF2 = IMATIJ(ISYMJ,ISYMM)  + 1
         KOFF3 = IT2BCD(ISYMAI,ISYMJ) + 1
C
         NTOAI = MAX(NT1AM(ISYMAI),1)
         NTOTJ = MAX(NRHF(ISYMJ),1)
C
         CALL DGEMM('N','T',NT1AM(ISYMAI),NRHF(ISYMJ),NRHF(ISYMM),-ONE,
     *              T2TINT(KOFF1),NTOAI,XMAT(KOFF2),NTOTJ,ONE,
     *              D2IAJ(KOFF3),NTOAI)
C
  100 CONTINUE
C
      CALL QEXIT('DIAJ23')
C
      RETURN
      END
C  /* Deck diaj24 */
      SUBROUTINE DIAJ24(D2IAJ,ISYIAJ,T2AMT,DAB,ISYD1,XLAMDH,ISYMLH,
     *                  WORK,LWORK,IDEL,ISYMD)
C
C     Written by Asger Halkier 13/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To calculate the fourth contribution to D2IAJ
C              from projection against doubles!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IAJ(*), T2AMT(*), DAB(*), XLAMDH(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      CALL QENTER('DIAJ24')
C
      ISYMB  = MULD2H(ISYMD,ISYMLH)
      ISYMC  = MULD2H(ISYMB,ISYD1)
C
      IF (MULD2H(ISYMLH,ISYMD).NE.MULD2H(ISYD1,ISYIAJ)) THEN
        CALL QUIT('Symmetry mismatch in DIAJ24.')
      END IF
C
      IF (NVIR(ISYMB) .EQ. 0) THEN
        CALL QEXIT('DIAJ24')
        RETURN 
      END IF
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KVECC = 1
      KEND1 = KVECC + NVIR(ISYMC)
      LWRK1 = LWORK - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient work space for allocation in DIAJ24')
      ENDIF
C
      CALL DZERO(WORK(KVECC),NVIR(ISYMC))
C
C-------------------------------------------------------------------
C     Calculate contraction of transposed Y-matrix (DAB) and XLAMDH.
C-------------------------------------------------------------------
C
      KOFF1 = IMATAB(ISYMC,ISYMB) + 1
      KOFF2 = IGLMVI(ISYMD,ISYMB) + IDEL - IBAS(ISYMD)
C
      NTOTC = MAX(NVIR(ISYMC),1)
C
      CALL DGEMV('N',NVIR(ISYMC),NVIR(ISYMB),ONE,DAB(KOFF1),NTOTC,
     *           XLAMDH(KOFF2),NBAS(ISYMD),ZERO,WORK(KVECC),1)
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      DO 100 ISYMJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMJ,ISYIAJ)
         ISYMCJ = ISYMAI
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            DO 120 C = 1,NVIR(ISYMC)
C
               NCJ = IT1AM(ISYMC,ISYMJ) + NVIR(ISYMC)*(J - 1) + C
C
               DO 130 NAI = 1,NT1AM(ISYMAI)
C
                  NAICJ = IT2AM(ISYMAI,ISYMCJ) + INDEX(NAI,NCJ)
                  NAIJ  = IT2BCD(ISYMAI,ISYMJ) + NT1AM(ISYMAI)*(J - 1)
     *                  + NAI
C
                  D2IAJ(NAIJ) = D2IAJ(NAIJ) 
     *                        - T2AMT(NAICJ)*WORK(KVECC + C - 1)
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('DIAJ24')
C
      RETURN
      END
C  /* Deck ccd1reor */
      SUBROUTINE CCD1REOR(DFOCK,DAI,DAB,DIJ,DIA)
C
C     Written by Asger Halkier 28/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To reorder the four separate blocks DAI through DIA
C              of the CC one electron density to one matrix stored
C              as a directly diagonalisable Fock-matrix (DFOCK)!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION DFOCK(*), DAI(*), DAB(*), DIJ(*), DIA(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCD1REOR')
C
      DO 100 ISYM = 1,NSYM
C
C---------------------------------------
C        Do the occupied occupied block.
C---------------------------------------
C
         DO 110 J = 1,NRHF(ISYM)
C
            NIJOR = IMATIJ(ISYM,ISYM) + NRHF(ISYM)*(J - 1) + 1
            NIJNE = IFCKDO(ISYM) + NORB(ISYM)*(J - 1) + 1
C
            CALL DCOPY(NRHF(ISYM),DIJ(NIJOR),1,DFOCK(NIJNE),1)
C
  110    CONTINUE
C
C--------------------------------------
C        Do the virtual occupied block.
C--------------------------------------
C
         DO 120 I = 1,NRHF(ISYM)
C
            NAIOR = IT1AM(ISYM,ISYM)  + NVIR(ISYM)*(I - 1) + 1
            NAINE = IFCKDO(ISYM) + NORB(ISYM)*(I - 1) 
     *            + NRHF(ISYM) + 1
C
            CALL DCOPY(NVIR(ISYM),DAI(NAIOR),1,DFOCK(NAINE),1)
C
  120    CONTINUE
C
C--------------------------------------
C        Do the occupied virtual block.
C--------------------------------------
C
         DO 130 A = 1,NVIR(ISYM)
C
            NIAOR = IT1AM(ISYM,ISYM)  + A
            NIANE = IFCKDV(ISYM) + NORB(ISYM)*(A - 1) + 1
C
            CALL DCOPY(NRHF(ISYM),DIA(NIAOR),NVIR(ISYM),
     *                 DFOCK(NIANE),1)
C
  130    CONTINUE
C
C-------------------------------------
C        Do the virtual virtual block.
C-------------------------------------
C
         DO 140 B = 1,NVIR(ISYM)
C
            NABOR = IMATAB(ISYM,ISYM) + NVIR(ISYM)*(B - 1) + 1
            NABNE = IFCKDV(ISYM) + NORB(ISYM)*(B - 1)
     *            + NRHF(ISYM) + 1
C
            CALL DCOPY(NVIR(ISYM),DAB(NABOR),1,DFOCK(NABNE),1)
C
  140    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CCD1REOR')
C
      RETURN
      END
C  /* Deck ccd1t1co */
      SUBROUTINE CCD1T1CO(DAI,DAB,DIJ,DIA,T1AM,YMAT,WORK,LWORK)
C
C     Written by Asger Halkier 29/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To add contributions from the t1-amplitudes to
C              the CC one electron density in MO-basis.
C                                             --
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION DAI(*), DAB(*), DIJ(*), DIA(*)
      DIMENSION T1AM(*), YMAT(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCD1T1CO')
C
C----------------------------------------------
C     Add contributions to the (occ,vir) block.
C----------------------------------------------
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = ISYMA
         ISYMK = ISYMA
         ISYMC = ISYMA
C
         KOFF1 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF2 = IMATIJ(ISYMI,ISYMK) + 1
         KOFF3 = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTI = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *              ONE,T1AM(KOFF1),NTOTA,DIJ(KOFF2),NTOTI,ONE,
     *              DIA(KOFF3),NTOTA)
C
         KOFF4 = IMATAB(ISYMA,ISYMC) + 1
         KOFF5 = IT1AM(ISYMC,ISYMI)  + 1
         KOFF6 = IT1AM(ISYMA,ISYMI)  + 1
C
         NTOTA = MAX(NVIR(ISYMA),1)
         NTOTC = MAX(NVIR(ISYMC),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NVIR(ISYMC),
     *              -ONE,YMAT(KOFF4),NTOTA,T1AM(KOFF5),NTOTC,ONE,
     *              DIA(KOFF6),NTOTA)
C
C------------------------------
C        Work space allocation.
C------------------------------
C
         KSCR  = 1
         KEND1 = KSCR  + NRHF(ISYMK)*NRHF(ISYMI)
         LWRK1 = LWORK - KEND1
C
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
            CALL QUIT('Insufficient work space for '//
     &           'allocation in CCD1T1CO')
         ENDIF
C
         CALL DZERO(WORK(KSCR),NRHF(ISYMK)*NRHF(ISYMI))
C
         KOFF7 = IT1AM(ISYMC,ISYMK) + 1
         KOFF8 = IT1AM(ISYMC,ISYMI) + 1
C
         NTOTC = MAX(NVIR(ISYMC),1)
         NTOTK = MAX(NRHF(ISYMK),1)
C
         CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMI),NVIR(ISYMC),
     *              ONE,DAI(KOFF7),NTOTC,T1AM(KOFF8),NTOTC,ZERO,
     *              WORK(KSCR),NTOTK)
C
         KOFF9  = IT1AM(ISYMA,ISYMK) + 1
         KOFF10 = IT1AM(ISYMA,ISYMI) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTK  = MAX(NRHF(ISYMK),1)
C
         CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),NRHF(ISYMK),
     *              -ONE,T1AM(KOFF9),NTOTA,WORK(KSCR),NTOTK,ONE,
     *              DIA(KOFF10),NTOTA)
C
  100 CONTINUE
C
C----------------------------------------------
C     Add contributions to the (vir,vir) block.
C----------------------------------------------
C
      DO 110 ISYMA = 1,NSYM
C
         ISYMB  = ISYMA
         ISYMK  = ISYMA
C
         KOFF11 = IT1AM(ISYMA,ISYMK)  + 1
         KOFF12 = IT1AM(ISYMB,ISYMK)  + 1
         KOFF13 = IMATAB(ISYMA,ISYMB) + 1
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTB  = MAX(NVIR(ISYMB),1)
C
         CALL DGEMM('N','T',NVIR(ISYMA),NVIR(ISYMB),NRHF(ISYMK),
     *              ONE,DAI(KOFF11),NTOTA,T1AM(KOFF12),NTOTB,ONE,
     *              DAB(KOFF13),NTOTA)
C
  110 CONTINUE
C
C----------------------------------------------
C     Add contributions to the (occ,occ) block.
C----------------------------------------------
C
      DO 120 ISYMI = 1,NSYM
C
         ISYMJ  = ISYMI
         ISYMC  = ISYMI
C
         KOFF14 = IT1AM(ISYMC,ISYMI)  + 1
         KOFF15 = IT1AM(ISYMC,ISYMJ)  + 1
         KOFF16 = IMATIJ(ISYMI,ISYMJ) + 1
C
         NTOTC  = MAX(NVIR(ISYMC),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
C
         CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NVIR(ISYMC),
     *              -ONE,T1AM(KOFF14),NTOTC,DAI(KOFF15),NTOTC,ONE,
     *              DIJ(KOFF16),NTOTI)
C
  120 CONTINUE
C
      CALL QEXIT('CCD1T1CO')
C
      RETURN
      END
C  /* Deck sortash */
      SUBROUTINE SORTASH(XARR,YARR,N)
C
C     Written by Asger Halkier 4/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To reorder the array xarr after size, and
C              do the same reordering of the concomitant
C              array yarr!
C              (This routine is based on Bubble-sort from NUM.REC.)
C
#include <implicit.h>
      DIMENSION XARR(*), YARR(*)
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('SORTASH')
C
C----------------------
C     Do the resorting.
C----------------------
C
      DO 100 J = 2,N
C
         XA = XARR(J)
         XB = YARR(J)
C
         DO 110 I = J-1,1,-1
C
            IF (XARR(I) .LE. XA) GOTO 120
C
            XARR(I + 1) = XARR(I)
            YARR(I + 1) = YARR(I)
C
  110    CONTINUE
C
         I=0
C
  120    XARR(I + 1) = XA
         YARR(I + 1) = XB
C
  100 CONTINUE
C
      CALL QEXIT('SORTASH')
C
      RETURN
      END 
C  /* Deck ccnaocan */
      SUBROUTINE CCNAOCAN(XARR,YARR)
C
C     Written by Asger Halkier 4/2 - 1996
C
C     Version: 1.0
C
C     Purpose: To analyze and check the natural occupation
C              numbers in xarr and yarr (real and imag)!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION XARR(*), YARR(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      CALL QENTER('CCNAOCAN')
C
      SUMRHF = ZERO
      SUMALL = ZERO
      SUMIMA = ZERO
C
      DO 100 I = 1,NORBT
C
         SUMALL = SUMALL + XARR(I)
         SUMIMA = SUMIMA + YARR(I)
C
  100 CONTINUE
C
      CHECK = 2*NRHFT - SUMALL
      CHEPR = NRHFT*2.0D0
      THR   = 1.0D-06
C
      IF (CHECK .LT. THR) THEN
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,111) 'Total Sum of natural occupation numbers:',
     &        CHEPR
C
      ELSE
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,111) 'Total Sum of natural occupation numbers:',
     &        SUMALL
C
      ENDIF
C
      IF (IPRINT .GT. 20) THEN
C
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,111) 'Total Sum of imaginary parts of nat.occ:',
     &        SUMIMA
C
      ENDIF
C
  111 FORMAT(3X,A40,2X,F9.6)
C
      DO 110 I = 1,NRHFT
C
         SUMRHF = SUMRHF + XARR(NORBT - I + 1)
C
  110 CONTINUE
C
      TOMOVE = 2*NRHFT - SUMRHF
C
      WRITE(LUPRI,*) ' '
      WRITE (LUPRI,222) 'Dynamical correlation move:',
     *                   TOMOVE,'electrons'
C
  222 FORMAT(3X,A27,F9.6,1X,A9)
C
      CALL QEXIT('CCNAOCAN')
C
      RETURN
      END
C  /* Deck cc_fcd1ao */
      SUBROUTINE CC_FCD1AO(AODEN,WORK,LWORK,MODEL)
C
C     Written by Ove Christiansen 2-may 1996.
C
C     Version: 1.0
C
C     Purpose: To calculate the Frozen core contribution to the one electron Coupled Cluster
C              density matrix and add it to AODEN.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION AODEN(*), WORK(LWORK)
#include <ccorb.h>
#include <infind.h>
#include <inftap.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
C
      CHARACTER MODEL*10
C
      CALL QENTER('CC_FCD1AO')
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND( 'Calculate core contribution to density')
      ENDIF
C
      KCMO  = 1
      KEND  = KCMO  + NLAMDS
      LWRK1 = LWORK - KEND
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient spaces in CC_FCD1AO')
      ENDIF
C
C----------------------------------------------
C     Read MO-coefficients from interface file.
C----------------------------------------------
C
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
C
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      ICMOI  = KCMO - 1
      DO ISYMI = 1, NSYM
        ISYMA = ISYMI
        ISYMB = ISYMI
        DO IFR = 1, NRHFFR(ISYMI)
          I = KFRRHF(IFR,ISYMI)            
          DO IAL = 1, NBAS(ISYMA) 
            ICMOIA = ICMOI + NBAS(ISYMA)*(I-1) + IAL
            DO IBE = 1, NBAS(ISYMA) 
              IDAB = IAODIS(ISYMA,ISYMB)+NBAS(ISYMA)*(IBE-1)+IAL
              ICMOIB = ICMOI + NBAS(ISYMB)*(I-1) + IBE
              AODEN(IDAB) = AODEN(IDAB) + TWO*WORK(ICMOIA)*WORK(ICMOIB)
            ENDDO
          ENDDO
        ENDDO
        ICMOI  = ICMOI  + NBAS(ISYMI)*NRHFS(ISYMI)
        ICMOI  = ICMOI  + NBAS(ISYMI)*NVIRS(ISYMI)
      ENDDO
C
      CALL QEXIT('CC_FCD1AO')
C
      END
C  /* Deck cc2_d1ao */
      SUBROUTINE CC2_D1AO(AODEN,ZKAM,WORK,LWORK)
C
C     Written by A. Halkier & S. Coriani 13/01-2000.
C
C     Version: 1.0
C
C     Purpose: To calculate the relaxed one electron CC2
C              density matrix and return it backtransformed
C              to AO-basis in AODEN! ZKAM is kappa-bar-0 
C              (input), which contributes to the relaxed density.
C
C     Important note: NO FROZEN CORE SO FAR!!!
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <ccinftap.h>
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <iratdef.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION INDEXA(MXCORB)
      DIMENSION AODEN(*), ZKAM(*), WORK(LWORK)
      CHARACTER*10 MODEL
#include <ccorb.h>
#include <infind.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <eribuf.h>
#include <ccfop.h>
#include <ccfro.h>
C
      LOGICAL LOCDBG
      PARAMETER ( LOCDBG=.FALSE. )
C
      CALL QENTER('CC2_D1AO')
C
C------------------------------------------------------------------
C     Both t-vectors and tbar-vectors (zeta) are totally symmetric.
C------------------------------------------------------------------
C
      ISYMTR = 1
      ISYMOP = 1
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KT2AM  = 1
      KZ2AM  = KT2AM  + NT2AMX
      KLAMDP = KZ2AM  + NT2SQ(1)
      KLAMDH = KLAMDP + NLAMDT
      KT1AM  = KLAMDH + NLAMDT
      KZ1AM  = KT1AM  + NT1AMX
      KEND1  = KZ1AM  + NT1AMX
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:',
     *                   LWORK, 'Needed:', KEND1
         CALL QUIT(
     *    'Insufficient memory for initial allocation in CC2_D1AO')
      ENDIF
C
C----------------------------------------
C     Read zero'th order zeta amplitudes.
C----------------------------------------
C
      IOPT   = 3
      CALL CC_RDRSP('L0',0,1,IOPT,MODEL,WORK(KZ1AM),WORK(KZ2AM))
C
C--------------------------------
C     Square up zeta2 amplitudes.
C--------------------------------
C
      CALL DCOPY(NT2AMX,WORK(KZ2AM),1,WORK(KT2AM),1)
      CALL CC_T2SQ(WORK(KT2AM),WORK(KZ2AM),1)
C
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
C
C----------------------------------
C     Calculate the lamda matrices.
C----------------------------------
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *            LWRK1)
C
C-----------------------------------------------
C     Set up 2C-E of cluster amplitudes and save
C     in KT2AM, as we only need T(2c-e) below.
C-----------------------------------------------
C
      ISYOPE = 1
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KT2AM),1.0D0,ISYOPE,IOPTTCME)
      KT2AMT = KT2AM                  !for safety
C
c     write(6,*) 'after initial stuff'
c     CALL FLSHFO(LUPRI)
C
C-------------------------------
C     Work space allocation one.
C     Note that D(ai) = ZETA(ai)
C     and both D(ia) and h(ia) 
C     are stored transposed!
C-------------------------------
C
      KONEAI = KZ1AM
      KONEAB = KONEAI + NT1AMX
      KONEIJ = KONEAB + NMATAB(1)
      KONEIA = KONEIJ + NMATIJ(1)
      KCMO   = KONEIA + NT1AMX
      KDUM   = KCMO   + NLAMDS
      KEND1  = KDUM   + NMATIJ(1)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation 1 CC2_D1AO')
      ENDIF
C
C
C------------------------------------------------------
C     Initialize remaining one electron density arrays.
C------------------------------------------------------
C
      CALL DZERO(WORK(KONEAB),NMATAB(1))
      CALL DZERO(WORK(KONEIJ),NMATIJ(1))
      CALL DZERO(WORK(KONEIA),NT1AMX)
C
C--------------------------------------------------------
C     Construct remaining blocks of one electron density.
C--------------------------------------------------------
C
      CALL DZERO(WORK(KDUM),NMATIJ(1))
      CALL DIJGEN(WORK(KONEIJ),WORK(KDUM))
      CALL DIAGEN(WORK(KONEIA),WORK(KT2AMT),WORK(KONEAI))
C
c     write(6,*) 'after D-lambda'
c     CALL FLSHFO(LUPRI)
C
      IF (LOCDBG) THEN
         DIJNO = DDOT(NMATIJ(1),WORK(KONEIJ),1,WORK(KONEIJ),1)
         DAINO = DDOT(NT1AMX,WORK(KONEAI),1,WORK(KONEAI),1)
         DIANO = DDOT(NT1AMX,WORK(KONEIA),1,WORK(KONEIA),1)
         DABNO = DDOT(NMATAB(1),WORK(KONEAB),1,WORK(KONEAB),1)
         WRITE(LUPRI,*) ' '
         WRITE(LUPRI,*) 'Norm of one electron densities in MO-basis:'
         WRITE(LUPRI,*) DIJNO, DAINO, DIANO, DABNO
      ENDIF
C
C--------------------------------------------------------
C     Backtransform the one electron density to AO-basis.
C--------------------------------------------------------
C
      CALL DZERO(AODEN,N2BST(1))
C
      ISDEN = 1
      CALL CC_DENAO(AODEN,1,WORK(KONEAI),WORK(KONEAB),
     *              WORK(KONEIJ),WORK(KONEIA),ISDEN,WORK(KLAMDP),1,
     *              WORK(KLAMDH),1,WORK(KEND1),LWRK1)
C
C----------------------------------------------------------
C     Read MO-coefficients from interface file and reorder.
C----------------------------------------------------------
C
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','UNKNOWN',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
      CALL GPCLOSE (LUSIFC,'KEEP')
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C----------------------------------
C     Add orbital relaxation terms.
C----------------------------------
C
      KOFDIJ = 1
      KOFDAB = KOFDIJ + NMATIJ(1)
      KOFDAI = KOFDAB + NMATAB(1)
      KOFDIA = KOFDAI + NT1AMX
C
      ISDEN = 1
      CALL CC_DENAO(AODEN,1,ZKAM(KOFDAI),ZKAM(KOFDAB),
     *              ZKAM(KOFDIJ),ZKAM(KOFDIA),ISDEN,WORK(KCMO),1,
     *              WORK(KCMO),1,WORK(KEND1),LWRK1)
C
      CALL QEXIT('CC2_D1AO')
      RETURN
      END
