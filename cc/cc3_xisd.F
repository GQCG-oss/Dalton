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
C  /* Deck CC3_XISD */ 
      SUBROUTINE CC3_XISD(XI1,XI2,XI1EFF,XI2EFF,ISYRES,
     *                  FOCKY,ISYFKY,FREQY,ICAU,LCAUCHY,LABELB,
     *                  XLAMDP,XLAMDH,FOCK0,
     *                  LUCKJD,FNCKJD,LUDKBC,FNDKBC,LUDELD,FNDELD,
     *                  LUTOC,FNTOC,LU3VI,FN3VI,LU3VI2,FN3VI2,
     *                  WORK,LWORK)
C
C     Written by F. Pawlowski and P Jorgensen , Spring 2002.
C
C     Partioning the triples part of the right hand side amplitude
C     equation \xi^Y into the singles and doubles space.
C
C     Initially construct 
C
C     t3(ai,Bj,Ck) = S^BC(aikj) + U^BC(aikj) + S^CB(aijk) + U^CB(aijk)
C
C     W^BC(aikj) = sum_d X(ad)*t3(di,Bj,Ck) - sum_l X(li)*t3(al,Bj,Ck)
C                 - sum_ld X(ld)*t2(ai,Cl)*t2(Bj,dk) 
C                 - sum_ld X(ld)*t2(ai,Bl)*t2(Ck,dj) 
C
C     Kasper Hald, Summer 2002
C     Generalized to calculate the entire T3^Y and
C     calculate the contributions for F.
C     add ( -<mu3 | [[H,T1^Y],T2] | HF> -<mu3| [H,T2^Y] |HF> )
C
C
C     Filip Pawlowski, Spring 2004, Aarhus, generalized to treat Cauchy vectors.
C
      IMPLICIT NONE
C
      LOGICAL LCAUCHY
C
      INTEGER ISYRES, ISINT1, ISINT2, ISYMT3, ISYMW
      INTEGER ISYFKY, ISYMT2, IOPT, KFOCKD, KOMG1
      INTEGER KOMG22, KFCKBA, KT2AM, KEND0, LWRK0, LWORK
      INTEGER ISYMC, ISYMK, ISYMD, ISYMB
      INTEGER KOFF1, KOFF2
      INTEGER KTROC0, KTROC02, KEND1, LWRK1
      INTEGER KINTOC, LWRK2, KEND2
      INTEGER IOFF, ISAIJ1, ISAIJ2, KRMAT1, KRMAT2, KFCKY
      INTEGER ISCKB1, ISCKB2 
      INTEGER KTRVI, KTRVI3, KTRVI1, KTRVI2, KTRVI0, KEND3, LWRK3
      INTEGER KINTVI, KEND4, LWRK4
      INTEGER ISYALJ, ISYALJ2, ISYMBD, ISCKIJ, ISCKD2
      INTEGER KSMAT, KSMAT3, KUMAT, KUMAT3, KDIAG, KINDSQ
      INTEGER KINDEX, KINDEX2, KTMAT, KTRVI8, KTRVI9
      INTEGER KTRVI10, KEND5, LWRK5, LENSQ, LUFCK  
      INTEGER LUCKJD, LUDKBC, LUDELD, LUTOC, LU3VI, LU3VI2
      INTEGER KWMAT, ISWMAT, LENSQW, KINDSQW
      INTEGER KXIAJB, KDIAGW, KTROC, KTROC1
      INTEGER LENGTH, ISYOPE, IOPTTCME, KRBJIA
      INTEGER KT1B, KT2B
      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR, LUDKBCR4
      INTEGER KTROC0Y, KTRVI0Y, KTRVI8Y
C
      INTEGER ICAU
C
      integer kx3am
C
      INTEGER IDLSTB,ISYMBB,NCAU,KC1,KC2,MCAU,IDLSTBM,KOCCC1,ISYMBBD
      INTEGER KVIRDC1,ISYMBBB,KVIRBC1,KOFFOCC,KOFFINTD,KOFFINTB
C
      INTEGER ISYMLSTB
C
Cfunctions
      INTEGER ILRCAMP

C
      REAL*8  XI1(*), XI2(*), XI1EFF(*), XI2EFF(*) 
      REAL*8  FOCKY(*), FOCK0(*), XLAMDP(*), XLAMDH(*)
      REAL*8  XINT, XTROC0, XTRVI0, XTRVI2, XDIA, XSMAT, XT2TP
      REAL*8  XTROC02, XUMAT, FREQY, DDOT, XIAJB
      REAL*8  XXI1EFF, XXI2EFF, FREQB, TWO
      REAL*8  WORK(LWORK)
      CHARACTER*10 MODEL
      CHARACTER FNCKJD*8, FNDKBC*8, FNDELD*8
      CHARACTER FNTOC*8, FN3VI2*8, LABELB*8
      CHARACTER FN3VI*6
      CHARACTER CDUMMY*1
      CHARACTER*16 FNCKJDR,FNDKBCR
      CHARACTER*13 FN3SRTR,FNDELDR
C
      PARAMETER(FN3SRTR  = 'CC3_XISD_TMP1', FNDELDR  = 'CC3_XISD_TMP2')
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccr1rsp.h>
#include <ccrc1rsp.h>
C
      PARAMETER(TWO = 2.0D0, CDUMMY = ' ')
C
      CALL QENTER('CC3_XISD')
C
      !Check value of ICAU
      IF (LCAUCHY .AND. ICAU.LT.-1) THEN
         WRITE(LUPRI,*)'ICAU = ', ICAU
         WRITE(LUPRI,*)'No CC3 cauchy implementation for ICAU.LT.-1'
         CALL QUIT('Illegal value of ICAU in CC3_XISD')
      END IF
C
      !Set variables associated with cauchy vectors
      IF (LCAUCHY .AND. ICAU.GE.0) THEN
         ISYMLSTB = 1
         IDLSTB  = ILRCAMP(LABELB,ICAU,ISYMLSTB)
         ISYMBB  = ISYLRC(IDLSTB)
c        FREQB  = ZERO
         LABELB = LRCLBL(IDLSTB)
         NCAU   = ILRCAU(IDLSTB)
C
         !Check symmetry
         IF (ISYMBB .NE. ISYFKY) THEN
          WRITE(LUPRI,*)'Symmetry of the perturbation operator: ',ISYFKY
          WRITE(LUPRI,*)'Symmetry of the Cauchy vector:         ',ISYMBB
          CALL QUIT('Symmetry inconsistency in CC3_XISD (Cauchy vec)')
         END IF
C
      END IF !LCAUCHY
C
C--------------------------
C     Save and set flags.
C--------------------------
C
C     Set symmetry flags.
C
C
C     ISYMT2 is symmetry of T2TP
C     ISYFKY is symmetry of perturbation operator
C     ISINT2 is symmetry of integrals in triples equation (ISINT2=1)
C     ISINT1 is symmetry of integrals in contraction (ISINT1=1)
C     ISYMT3  is symmetry of S and U intermediate
C     ISYMW   is symmetry of W intermediate
C     ISYRES  is symmetry of result vector (here the same as ISYFKY)
C     ISYMW  = ISYFKY*ISYMT3
C     ISYRES = ISYMT3*ISYFKY*ISINT1
C
C-------------------------------------------------------------
C
      IPRCC   = IPRINT
      ISINT1  = 1
      ISINT2  = 1
      ISYMT2  = 1
      ISYMT3  = MULD2H(ISINT2,ISYMT2)
      ISYMW   = MULD2H(ISYMT3,ISYFKY)
      ISYRES  = MULD2H(ISYMW,ISINT1)
C
C---------------------------------------------------------
C     Work allocation 
C---------------------------------------------------------
C
      KT2AM  = 1 
      KRBJIA = KT2AM  + NT2SQ(ISYMT2)
      KFOCKD = KRBJIA + NT2SQ(ISYRES)
      KFCKBA = KFOCKD + NORBTS
      KFCKY  = KFCKBA + N2BST(ISINT1)
      KEND0  = KFCKY  + N2BST(ISYFKY)
      LWRK0  = LWORK  - KEND0
C
      IF (LCAUCHY) THEN
         KC1   = KEND0
         KC2   = KC1   + (NCAU+1)*NT1AM(ISYMBB)
         KEND0 = KC2   + (NCAU+1)*NT2SQ(ISYMBB)
         LWRK0 = LWORK - KEND0
      END IF
C
      IF (LWRK0.LT.0) THEN
         WRITE(LUPRI,*)'Memory available :',LWORK
         WRITE(LUPRI,*)'Memory needed    :',KEND0
         CALL QUIT('Insufficient memory in CC3_XISD(c1)')
      END IF
C
      CALL DZERO(WORK(KRBJIA),NT2SQ(ISYRES))
      CALL DZERO(XI1EFF,NT1AM(ISYRES))
C
C--------------------------------------
C     Read in t2 , square it and reorder 
C--------------------------------------
C
      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYMT2,IOPT,MODEL,DUMMY,WORK(KEND0))
      CALL CC_T2SQ(WORK(KEND0),WORK(KT2AM),ISYMT2)
      IF (LWORK .LT. NT2SQ(ISYMT2)) THEN
         CALL QUIT('Not enough memory to construct T2TP in CC3_XISD') 
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMT2),WORK(KT2AM),1,WORK(KEND0),1)
      CALL CC3_T2TP(WORK(KT2AM),WORK(KEND0),ISYMT2)
C
      IF (IPRINT .GT. 55) THEN
         XT2TP = DDOT(NT2SQ(ISYMT2),WORK(KT2AM),1,WORK(KT2AM),1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XT2TP
      ENDIF
C
C
      IF (LCAUCHY .AND. ICAU.GE.0) THEN
C
         IF (LWRK0.LT.NT2SQ(ISYMBB)) THEN
            WRITE(LUPRI,*)'Memory available :',LWRK0
            WRITE(LUPRI,*)'Memory needed    :',NT2SQ(ISYMBB)
            CALL QUIT('Insufficient memory in CC3_XISD(c2)')
         END IF
C
C--------------------------------------------------------------------------
C        Loop over cauchy vectors in order to
C        1) get C1 and C2 cauchy vectors
C        2) get C1-transformed integrals needed in <mu3|[[H^0,C1],T2^0]|HF>
C--------------------------------------------------------------------------
C
         DO MCAU = 0, NCAU
C
            !Open temporary files
            LU3SRTR  = -1
            LUDELDR  = -1
            CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
            CALL WOPEN2(LUDELDR,FNDELDR,64,0)
C
            !Get the list for current MCAU
            IDLSTBM = ILRCAMP(LABELB,MCAU,ISYMBB)
C
            !Consistency check
            IF (MCAU.EQ.NCAU .AND. IDLSTBM.NE.IDLSTB) THEN
               CALL QUIT('Inconsistency in Cauchy loop in CC3_XISD')
            END IF
C
            !Read in C1 and C2 Cauchy vectors for each MCAU
            IOPT = 3
            KOFF1 = KC1 + MCAU*NT1AM(ISYMBB)
            KOFF2 = KC2 + MCAU*NT2SQ(ISYMBB)
            CALL CC_RDRSP('RC ',IDLSTBM,ISYMBB,IOPT,MODEL,WORK(KOFF1),
     *                    WORK(KOFF2))
            CALL CCLR_DIASCL(WORK(KOFF2),TWO,ISYMBB)
            CALL CC_T2SQ(WORK(KOFF2),WORK(KEND0),ISYMBB)
            CALL CC3_T2TP(WORK(KOFF2),WORK(KEND0),ISYMBB)
C
            !Generate file names for C1-transformed integrals
            CALL GET_FILENAME(FNCKJDR,FNDKBCR,MCAU)
C
            !Open the files for C1-transformed integrals
            LUCKJDR  = -1
            LUDKBCR  = -1
            CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
            CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
C
            !Construct C1-transformed integrals and store on file
            CALL CC3_BARINT(WORK(KOFF1),ISYMBB,XLAMDP,
     *                      XLAMDH,WORK(KEND0),LWRK0,
     *                      LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
            CALL CC3_SORT1(WORK(KEND0),LWRK0,2,ISYMBB,LU3SRTR,FN3SRTR,
     *                     LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *                     IDUMMY,CDUMMY)
C
            CALL CC3_SINT(XLAMDH,WORK(KEND0),LWRK0,ISYMBB,
     *                    LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)

C
            !Close the files for C1-transformed integrals
            CALL WCLOSE2(LUCKJDR,FNCKJDR,'KEEP')
            CALL WCLOSE2(LUDKBCR,FNDKBCR,'KEEP')
C
            !Close and delete temporary files
            CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
            CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
C
         END DO !MCAU
C
      END IF!LCAUCHY

C
C---------------------------------------------------------
C     Read canonical orbital energies.
C---------------------------------------------------------
C
C--------------------------------------
C     Read in orbital energies
C--------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C---------------------------------------------
C     Delete frozen orbitals in Fock diagonal.
C---------------------------------------------
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND0),LWRK0)
C
C If we want to sum the T3 amplitudes
      if (.false.) then
         kx3am  = kend0
         kend0  = kx3am + nt1amx*nt1amx*nt1amx
         call dzero(work(kx3am),nt1amx*nt1amx*nt1amx)
         lwrk0 = lwork - kend0
      endif
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in CC3_XISD')
      END IF
C
C-----------------------------------------------------
C     Construct the T1 transformed Fock matrix
C-----------------------------------------------------
C
      LUFCK = -1
C     This AO Fock matrix is constructed from the T1 transformed density
       CALL GPOPEN(LUFCK,'CC_FCKH','UNKNOWN',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
      REWIND(LUFCK)
      READ(LUFCK)(WORK(KFCKBA + I-1),I = 1,N2BST(ISINT1))
      CALL GPCLOSE(LUFCK,'KEEP' )
C
      IF (IPRINT .GT. 140) THEN
         CALL AROUND( 'Usual Fock AO matrix' )
         CALL CC_PRFCKAO(WORK(KFCKBA),1)
      ENDIF
C
      ! SCF Fock matrix in transformed using Lamda vector
      CALL CC_FCKMO(WORK(KFCKBA),XLAMDP,XLAMDH,
     *              WORK(KEND0),LWRK0,1,1,1)
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'In CC3_XISD: T1 transformed Fock matrix' )
         CALL CC_PRFCKMO(WORK(KFCKBA),ISINT1)
      ENDIF
C
C     Sort the fock matrix
C
      CALL DCOPY(N2BST(ISINT1),WORK(KFCKBA),1,WORK(KEND0),1)
C
      DO ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISINT1)
C
         DO K = 1,NRHF(ISYMK)
C
            DO C = 1,NVIR(ISYMC)
C
               KOFF1 = KEND0 + IFCVIR(ISYMK,ISYMC) + 
     *                 NORB(ISYMK)*(C - 1) + K - 1
               KOFF2 = KFCKBA + IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1) + C - 1
C
               WORK(KOFF2) = WORK(KOFF1)
C
            ENDDO
         ENDDO
      ENDDO
C
      IF (IPRINT .GT. 50) THEN
       CALL AROUND('In CC3_XISD: T1 transf Fock  matrix (sort(ck))')
         CALL CC_PRFCKMO(WORK(KFCKBA),1)
      ENDIF
C
C----------------------------------------------
C     Sort FOCKY to (ck)
C----------------------------------------------
C
      CALL DCOPY(N2BST(ISYFKY),FOCKY,1,WORK(KEND0),1)
C
      DO ISYMC = 1,NSYM
C
         ISYMK = MULD2H(ISYMC,ISYFKY)
C
         DO K = 1,NRHF(ISYMK)
C
            DO C = 1,NVIR(ISYMC)
C
               KOFF1 = KEND0 + IFCVIR(ISYMK,ISYMC) + 
     *                 NORB(ISYMK)*(C - 1) + K - 1
               KOFF2 = KFCKY + IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1) + C - 1
C
               WORK(KOFF2) = WORK(KOFF1)
C
            ENDDO
         ENDDO
      ENDDO
C
      IF (IPRINT .GT. 50) THEN
       CALL AROUND('In CC3_XISD: T1^Y transf Fock matrix (sort(ck))')
         CALL CC_PRFCKMO(WORK(KFCKY),1)
      ENDIF
C
C-----------------------------
C     Read occupied integrals.
C-----------------------------
C
C     Memory allocation.
C
      KTROC  = KEND0
      KTROC1 = KTROC  + NTRAOC(ISINT1)
      KTROC0 = KTROC1 + NTRAOC(ISINT1)
      KTROC02= KTROC0 + NTRAOC(ISINT2)
      KXIAJB = KTROC02+ NTRAOC(ISINT2)
      KEND1  = KXIAJB + NT2AM(ISINT1)
      LWRK1  = LWORK  - KEND1
C
      KINTOC = KEND1
      KEND2  = KINTOC + MAX(NTOTOC(ISINT2),NTOTOC(ISINT1))
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK
         WRITE(LUPRI,*) 'Memory needed    : ',KEND2
         CALL QUIT('Insufficient space in CC3_XISD')
      END IF
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      LENGTH = IRAT*NT2AM(ISINT1)
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,LENGTH,WORK(KXIAJB))
      ISYOPE = ISINT1
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISYOPE,IOPTTCME)
C
C
      IF ( IPRINT .GT. 55) THEN
         XIAJB = DDOT(NT2AM(ISINT1),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of IAJB ',XIAJB
      ENDIF
C
C------------------------
C     Occupied integrals.
C
C     Read in integrals for SMAT etc.
C-----------------------
C
      IOFF = 1
      IF (NTOTOC(ISINT2) .GT. 0) THEN
         CALL GETWA2(LUCKJD,FNCKJD,WORK(KINTOC),IOFF,NTOTOC(ISINT2))
      ENDIF
C
C----------------------------------
C     Write out norms of Integrals.
C----------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XINT
      ENDIF
C
C----------------------------------------------------------------------
C     Transform (ai|j delta) integrals to (ai|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C     here the xlamdp transformation need to be used 
C     (delta is particle index) 
C
      CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),XLAMDP,
     *               WORK(KEND2),LWRK2,ISINT2)

C
C----------------------------------------------------------------------
C     (ai|j k) sorted as (ij,k,a)
C----------------------------------------------------------------------
C
      CALL CCFOP_SORT(WORK(KTROC0),WORK(KTROC02),ISINT2,1)

C
C-------------------------------
C     Write out norms of arrays.
C-------------------------------
C
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',XINT
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTROC0 = DDOT(NTRAOC(ISINT2),WORK(KTROC0),1,
     *                WORK(KTROC0),1)
         WRITE(LUPRI,*) 'Norm of TROC0 ',XTROC0
      ENDIF
C
      IF (IPRINT .GT. 55) THEN
         XTROC02 = DDOT(NTRAOC(ISINT2),WORK(KTROC02),1,
     *                WORK(KTROC02),1)
         WRITE(LUPRI,*) 'Norm of TROC02 ',XTROC02
      ENDIF
C
C
C----------------------------------
C
C------------------------
C     Occupied integrals.
C
C     Read in integrals for WMAT etc.
C-----------------------
C
C
      IOFF = 1
      IF (NTOTOC(ISYMOP) .GT. 0) THEN
         CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISINT1))
      ENDIF
C
C----------------------------------
C     Write out norms of Integrals.
C----------------------------------
C
      IF (IPRINT .GT. 55) THEN
         XINT  = DDOT(NTOTOC(ISYMOP),WORK(KINTOC),1,
     *                WORK(KINTOC),1)
         WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XINT
      ENDIF
C
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (i,j,k,a)
C----------------------------------------------------------------------
C
       CALL CCSDT_TROCC(WORK(KINTOC),WORK(KTROC),XLAMDH,
     *                  WORK(KEND2),LWRK2)
C
C----------------------------------------------------------------------
C     sort (i,j,k,a) as (a,i,j,k)
C----------------------------------------------------------------------
C

      CALL CCSDT_SRTOC2(WORK(KTROC),WORK(KTROC1),ISINT1,
     *                  WORK(KEND2),LWRK2)
C
      IF (LCAUCHY .AND. ICAU.GE.0) THEN
C
         !allocate space for ALL C1-transformed occupied integrals 
         KOCCC1 = KEND1
         KEND1  = KOCCC1 + (NCAU+1)*NTRAOC(ISYMBB)
         LWRK1  = LWORK  - KEND1
C     
         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND1
            CALL QUIT('Insufficient space in CC3_XISD(c3)')
         END IF
C
         DO MCAU = 0,NCAU
C
C--------------------------------------------------------------
C          Read in C1-transformed occupied integrals for each MCAU
C-------------------------------------------------------------
C
            !Generate file name for C1-transformed integrals
            !(FNDKBCR is not needed here)
            CALL GET_FILENAME(FNCKJDR,FNDKBCR,MCAU)
C
            !Open the file for C1-transformed integrals
            LUCKJDR  = -1
            CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
C
            !Get the offset for a given MCAU
            KOFF1 = KOCCC1 + MCAU*NTRAOC(ISYMBB)
C
            !Read in the integrals
            CALL INTOCC_T3X(LUCKJDR,FNCKJDR,XLAMDP,ISYMBB,
     *                      WORK(KOFF1),WORK(KEND1),LWRK1)
C
            !Close and delete the file for C1-transformed occupied integrals
            CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
C
         END DO !MCAU
C
      END IF !LCAUCHY
C
C
C----------------------------
C     General loop structure.
C----------------------------
C
      DO ISYMD = 1,NSYM
C
         ISAIJ1  = MULD2H(ISYMD,ISYRES)
         ISCKB2  = MULD2H(ISINT2,ISYMD)
         ISCKB1  = MULD2H(ISINT1,ISYMD)
C
         IF (IPRINT .GT. 55) THEN
C
            WRITE(LUPRI,*) 'In CC3_XI: ISCKB2 :',ISCKB2
C
         ENDIF

C
C--------------------------
C        Memory allocation.
C--------------------------
C
         IF (LCAUCHY .AND. ICAU.GE.0) THEN
C
            !get symmetry of C1-transformed virtual integrals with fixed D
            ISYMBBD = MULD2H(ISYMBB,ISYMD)
C
            !allocation for ALL C1-transformed virtual integrals with fixed D
            KVIRDC1 = KEND1
            KEND1   = KVIRDC1 + (NCAU+1)*NCKATR(ISYMBBD)
            LWRK1   = LWORK  - KEND1
C     
            IF (LWRK1 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND1
               CALL QUIT('Insufficient space in CC3_XISD(c4)')
            END IF
C
         END IF !LCAUCHY
C
         KTRVI  = KEND1
         KTRVI1 = KTRVI  + NCKATR(ISCKB1)
         KTRVI3 = KTRVI1 + NCKATR(ISCKB1)
         KTRVI2 = KTRVI3 + NCKATR(ISCKB2)
         KRMAT1 = KTRVI2 + NCKATR(ISCKB2)
         KEND2  = KRMAT1 + NCKI(ISAIJ1)
         LWRK2  = LWORK  - KEND2
C
         KTRVI0  = KEND2
         KEND3   = KTRVI0  + NCKATR(ISCKB2)
         LWRK3   = LWORK  - KEND3
C
         KINTVI = KEND3
         KEND4  = KINTVI + MAX(NCKA(ISCKB2),NCKA(ISCKB1))
         LWRK4  = LWORK  - KEND4
C
         IF (LWRK4 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND4
            CALL QUIT('Insufficient space in CC3_XISD')
         END IF
C
C---------------------
C        Sum over D
C---------------------
C
         DO D = 1,NVIR(ISYMD)
C
C---------------------------------------
C           Initialise
C---------------------------------------
C
            CALL DZERO(WORK(KRMAT1),NCKI(ISAIJ1))
C
C------------------------------------------------------------
C           Read and transform integrals used in contraction
C           with W intermediate.
C------------------------------------------------------------
C
            IOFF = ICKAD(ISCKB1,ISYMD) + NCKA(ISCKB1)*(D - 1) + 1
            IF (NCKA(ISCKB1) .GT. 0) THEN
               CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     &                     NCKA(ISCKB1))
            ENDIF
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI),XLAMDP,
     *                       ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
            IOFF = ICKAD(ISCKB1,ISYMD) + NCKA(ISCKB1)*(D - 1) + 1
            IF (NCKA(ISCKB1) .GT. 0) THEN
               CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                     NCKA(ISCKB1))
            ENDIF
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI1),XLAMDP,
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
C-----------------------------------------------
C           Read virtual integrals used in first s3am.
C-----------------------------------------------
C
            IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
            IF (NCKATR(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                        NCKATR(ISCKB2))
            ENDIF
C
C-----------------------------------------------------------
C           Sort the integrals for s3am 
C-----------------------------------------------------------
C
            CALL CCSDT_SRTVIR(WORK(KTRVI0),WORK(KTRVI2),WORK(KEND4),
     *                        LWRK4,ISYMD,ISINT2)
C
C
            IF (IPRINT .GT. 55) THEN
               XTRVI0= DDOT(NCKATR(ISCKB2),WORK(KTRVI0),1,
     *                      WORK(KTRVI0),1)
               WRITE(LUPRI,*) 'Norm of TRVI0 ',XTRVI0
            ENDIF
C
            IF (IPRINT .GT. 55) THEN
               XTRVI2= DDOT(NCKATR(ISCKB2),WORK(KTRVI2),1,
     *                      WORK(KTRVI2),1)
               WRITE(LUPRI,*) 'Norm of TRVI2 ',XTRVI2
            ENDIF


C
C------------------------------------------------------
C           Read virtual integrals used in first U.
C------------------------------------------------------
C
            IOFF = ICKAD(ISCKB2,ISYMD) + NCKA(ISCKB2)*(D - 1) + 1
            IF (NCKA(ISCKB2) .GT. 0) THEN
               CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     &                     NCKA(ISCKB2))

            ENDIF
C
            CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),XLAMDH,
     *                       ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
C

            IF (LCAUCHY .AND. ICAU.GE.0) THEN
C
               DO MCAU = 0, NCAU
C
C---------------------------------------------------------------------------
C                 Read in C1-transformed virt int with fixed D for each MCAU
C---------------------------------------------------------------------------
C
                  !Generate file names for C1-transformed integrals
                  !(FNCKJDR is not needed here)
                  CALL GET_FILENAME(FNCKJDR,FNDKBCR,MCAU)
C
                  !Open the files for C1-transformed integrals
                  LUDKBCR  = -1
                  CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
C
                  !Get the offset for a given MCAU
                  KOFF1 = KVIRDC1 + MCAU*NCKATR(ISYMBBD)
C
                  !Read in the integrals
                  IOFF = ICKBD(ISYMBBD,ISYMD)+NCKATR(ISYMBBD)*(D-1)+1
                  IF (NCKATR(ISYMBBD) .GT. 0) THEN
                     CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KOFF1),IOFF,
     &                           NCKATR(ISYMBBD))
                  ENDIF
C
                  !Close the file for C1-transformed virtual integrals
                  !(and keep it: it will be needed in B-loop)
                  CALL WCLOSE2(LUDKBCR,FNDKBCR,'KEEP')
C
               END DO !MCAU
C
            END IF!LCAUCHY
C
C--------------------------------------------------------
C
            DO ISYMB = 1,NSYM
C
               ISAIJ2  = MULD2H(ISYMB,ISYRES)
               ISYALJ  = MULD2H(ISYMB,ISYMT2)
               ISYALJ2 = MULD2H(ISYMD,ISYMT2)
               ISYMBD  = MULD2H(ISYMB,ISYMD)
               ISCKIJ  = MULD2H(ISYMBD,ISYMT3)
               ISCKD2  = MULD2H(ISINT2,ISYMB)
               ISWMAT  = MULD2H(ISYFKY,ISCKIJ)
C
               IF (IPRINT .GT. 55) THEN
C 
                  WRITE(LUPRI,*) 'In CC3_XISD : ISYMD :',ISYMD
                  WRITE(LUPRI,*) 'In CC3_XISD : ISYMB :',ISYMB
                  WRITE(LUPRI,*) 'In CC3_XISD : ISYALJ:',ISYALJ
                  WRITE(LUPRI,*) 'In CC3_XISD : ISYMBD:',ISYMBD
                  WRITE(LUPRI,*) 'In CC3_XISD : ISCKIJ:',ISCKIJ
                  WRITE(LUPRI,*) 'In CC3_XISD : ISWMAT:',ISWMAT
C
               ENDIF
C
               !We can use KEND3 since we do not need KINTVI array from D-loop
               !anymore
C
               IF (LCAUCHY .AND. ICAU.GE.0) THEN
C
                  !get symmetry of C1-transformed virt int with fixed B
                  ISYMBBB = MULD2H(ISYMBB,ISYMB)
C
                  !allocation for ALL C1-transformed virt int with fixed B
                  KVIRBC1 = KEND3
                  KEND3   = KVIRBC1 + (NCAU+1)*NCKATR(ISYMBBB)
                  LWRK3   = LWORK  - KEND3
C     
                  IF (LWRK3 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Memory available : ',LWORK
                     WRITE(LUPRI,*) 'Memory needed    : ',KEND3
                     CALL QUIT('Insufficient space in CC3_XISD(c5)')
                  END IF
C
               END IF !LCAUCHY
C
               KSMAT   = KEND3
               KSMAT3  = KSMAT   + NCKIJ(ISCKIJ)
               KUMAT   = KSMAT3  + NCKIJ(ISCKIJ)
               KUMAT3  = KUMAT   + NCKIJ(ISCKIJ)
               KDIAG   = KUMAT3  + NCKIJ(ISCKIJ)
               KDIAGW  = KDIAG   + NCKIJ(ISCKIJ)
               KWMAT   = KDIAGW  + NCKIJ(ISWMAT)
               KRMAT2  = KWMAT   + NCKIJ(ISWMAT) 
               KINDSQW = KRMAT2  + NCKI(ISAIJ2) 
               KINDSQ  = KINDSQW + (6*NCKIJ(ISWMAT) - 1)/IRAT + 1
               KINDEX  = KINDSQ  + (6*NCKIJ(ISCKIJ) - 1)/IRAT + 1
               KINDEX2 = KINDEX  + (NCKI(ISYALJ) - 1)/IRAT + 1
               KTMAT   = KINDEX2 + (NCKI(ISYALJ2) - 1)/IRAT + 1
               KTRVI8  = KTMAT   + MAX(NCKIJ(ISCKIJ),NCKIJ(ISWMAT))
               KTRVI9  = KTRVI8  + NCKATR(ISCKD2)
               KTRVI10 = KTRVI9  + NCKATR(ISCKD2)
               KEND4   = KTRVI10 + NCKATR(ISCKD2)
               LWRK4   = LWORK   - KEND4
C
               KINTVI  = KEND4
               KEND5   = KINTVI  + NCKA(ISCKD2)
               LWRK5   = LWORK   - KEND5
C
               IF (LWRK5 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Memory available : ',LWORK
                  WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                  CALL QUIT('Insufficient space in CC3_XISD ')
               END IF
C
C---------------------------------------------
C              Construct part of the diagonal.
C---------------------------------------------
C
               CALL CC3_DIAG(WORK(KDIAG),WORK(KFOCKD),ISCKIJ)
               CALL CC3_DIAG(WORK(KDIAGW),WORK(KFOCKD),ISWMAT)
C
               IF (IPRINT .GT. 55) THEN
                  XDIA  = DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,
     *                    WORK(KDIAG),1)
                  WRITE(LUPRI,*) 'Norm of DIA  ',XDIA
                  XDIA  = DDOT(NCKIJ(ISCKIJ),WORK(KDIAGW),1,
     *                    WORK(KDIAGW),1)
                  WRITE(LUPRI,*) 'Norm of DIA_W',XDIA
               ENDIF
C
C-------------------------------------
C              Construct index arrays.
C-------------------------------------
C
               LENSQ = NCKIJ(ISCKIJ)
               CALL CC3_INDSQ(WORK(KINDSQ),LENSQ,ISCKIJ)
               CALL CC3_INDEX(WORK(KINDEX),ISYALJ)
               CALL CC3_INDEX(WORK(KINDEX2),ISYALJ2)
               LENSQW = NCKIJ(ISWMAT)
               CALL CC3_INDSQ(WORK(KINDSQW),LENSQW,ISWMAT) 
C
               DO B = 1,NVIR(ISYMB)
C
C---------------------------------------
C           Initialise
C---------------------------------------
C
                  CALL DZERO(WORK(KRMAT2),NCKI(ISAIJ2))
                  CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))
C
C-----------------------------------------------
C                 Read virtual integrals used in second s3am.
C-----------------------------------------------
C
                  IOFF = ICKBD(ISCKD2,ISYMB) + 
     &                   NCKATR(ISCKD2)*(B - 1) + 1
                  IF (NCKATR(ISCKD2) .GT. 0) THEN
                     CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI8),IOFF,
     &                           NCKATR(ISCKD2))
                  ENDIF
C
C-----------------------------------------------------------
C           Sort the integrals for s3am 
C-----------------------------------------------------------
C
                  CALL CCSDT_SRTVIR(WORK(KTRVI8),
     *                              WORK(KTRVI9),WORK(KEND4),
     *                              LWRK4,ISYMB,ISINT2)
C
C
C----------------------------------------------------------
C           Read virtual integrals used in second U
C----------------------------------------------------------
C
                  IOFF = ICKAD(ISCKD2,ISYMB) + NCKA(ISCKD2)*(B - 1) + 1
                  IF (NCKA(ISCKD2) .GT. 0) THEN
                     CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     *                           NCKA(ISCKD2))
                  ENDIF
C
                  CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI10),XLAMDH,
     *                             ISYMB,B,ISINT2,WORK(KEND5),LWRK5)
C
C
                  IF (LCAUCHY .AND. ICAU.GE.0) THEN
C
                     DO MCAU = 0, NCAU
C
C---------------------------------------------------------------------------
C                    Read in C1-transformed virt int with fixed B for each MCAU
C---------------------------------------------------------------------------
C
                        !Generate file names for C1-transformed integrals
                        !(FNCKJDR is not needed here)
                        CALL GET_FILENAME(FNCKJDR,FNDKBCR,MCAU)
C
                        !Open the files for C1-transformed integrals
                        LUDKBCR  = -1
                        CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
C
                        !Get the offset for a given MCAU
                        KOFF1 = KVIRBC1 + MCAU*NCKATR(ISYMBBB)
C
                        !Read in the integrals
                        IOFF = ICKBD(ISYMBBB,ISYMB)
     *                       +NCKATR(ISYMBBB)*(B-1)+1
                        IF (NCKATR(ISYMBBB) .GT. 0) THEN
                           CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KOFF1),IOFF,
     &                                 NCKATR(ISYMBBB))
                        ENDIF
C
                        !Close and keep the file for C1-transformed virt int
                        !...or delte it if you don't need it anymore

                        IF (      (ISYMD.EQ.NSYM) .AND. (ISYMB.EQ.NSYM)
     *                      .AND. (D.EQ.NVIR(ISYMD)) 
     *                      .AND. (B.EQ.NVIR(ISYMB))) THEN
C
                               CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
                        ELSE
                               CALL WCLOSE2(LUDKBCR,FNDKBCR,'KEEP')
                        END IF
C
                     END DO !MCAU
C
                  END IF!LCAUCHY

C
C------------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for B,D.
C-------------------------------------------------------------------
C
                  CALL CC3_SMAT(0.0D0,WORK(KT2AM),ISYMT2,
     *                          WORK(KTMAT),WORK(KTRVI0),
     *                          WORK(KTRVI2),WORK(KTROC0),ISINT2,
     *                          WORK(KFOCKD),WORK(KDIAG),
     *                          WORK(KSMAT),WORK(KEND4),LWRK4,
     *                          WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)
C
                  CALL T3_FORBIDDEN(WORK(KSMAT),ISYMT3,ISYMB,B,ISYMD,D)
C
c      call sum_pt3(work(ksmat),isymb,b,isymd,d,
c    *             isckij,work(kx3am),1)
C
                  IF (IPRINT .GT. 55) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                       WORK(KSMAT),1)
                     WRITE(LUPRI,*) 'Norm of SMAT     ',XSMAT
                  ENDIF
C
C-------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for D,B.
C-------------------------------------------------------------------
C
                  CALL CC3_SMAT(0.0D0,WORK(KT2AM),ISYMT2,
     *                          WORK(KTMAT),WORK(KTRVI8),
     *                          WORK(KTRVI9),WORK(KTROC0),ISINT2,
     *                          WORK(KFOCKD),WORK(KDIAG),
     *                          WORK(KSMAT3),WORK(KEND4),LWRK4,
     *                          WORK(KINDEX2),WORK(KINDSQ),LENSQ,
     *                          ISYMD,D,ISYMB,B)
C
                  CALL T3_FORBIDDEN(WORK(KSMAT3),ISYMT3,ISYMD,D,ISYMB,B)
C
c     call sum_pt3(work(ksmat3),isymd,d,isymb,b,
c    *             isckij,work(kx3am),1)
C
                  IF (IPRINT .GT. 55) THEN
                     XSMAT = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT3),1,
     *                       WORK(KSMAT3),1)
                     WRITE(LUPRI,*) 'Norm of SMAT3    ',XSMAT
                  ENDIF
C
C---------------------------------------------------------------------------
C                 Calculate U(ci,jk) for fixed b,d.
C--------------------------------------------------
C
                  CALL CC3_UMAT(0.0D0,WORK(KT2AM),ISYMT2,WORK(KTRVI3),
     *                          WORK(KTROC02),ISINT2,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KUMAT),WORK(KTMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMB,B,ISYMD,D)
C
                  CALL T3_FORBIDDEN(WORK(KUMAT),ISYMT3,ISYMB,B,ISYMD,D)
C
c      call sum_pt3(work(kumat),isymb,b,isymd,d,
c    *             isckij,work(kx3am),3)
C
                  IF (IPRINT .GT. 55) THEN
                     XUMAT = DDOT(NCKIJ(ISCKIJ),WORK(KUMAT),1,
     *                       WORK(KUMAT),1)
                     WRITE(LUPRI,*) 'Norm of UMAT     ',XUMAT
                  ENDIF
C
C--------------------------------------------------
C
C
C--------------------------------------------------
C                 Calculate U(ci,jk) for fixed d,b.
C--------------------------------------------------
C
                  CALL CC3_UMAT(0.0D0,WORK(KT2AM),ISYMT2,WORK(KTRVI10),
     *                          WORK(KTROC02),ISINT2,WORK(KFOCKD),
     *                          WORK(KDIAG),WORK(KUMAT3),WORK(KTMAT),
     *                          WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                          ISYMD,D,ISYMB,B)
C
                  CALL T3_FORBIDDEN(WORK(KUMAT3),ISYMT3,ISYMD,D,ISYMB,B)
C
c     call sum_pt3(work(kumat3),isymd,d,isymb,b,
c    *             isckij,work(kx3am),3)
C
                  IF (IPRINT .GT. 55) THEN
                     XUMAT = DDOT(NCKIJ(ISCKIJ),WORK(KUMAT3),1,
     *                       WORK(KUMAT3),1)
                     WRITE(LUPRI,*) 'Norm of UMAT3    ',XUMAT
                  ENDIF
C
C
C-------------------------------------------------
C     Calculate the term <mu2|[V,T3]|HF>
C-------------------------------------------------
C
                  IF (.NOT.(LCAUCHY .AND. (ICAU.GE.0))) THEN
                    CALL CCFOP_ONEL(XI2EFF,WORK(KRMAT1),WORK(KRMAT2),
     *                              WORK(KFCKY),WORK(KSMAT),WORK(KTMAT),
     *                              ISYMT3,ISYFKY,WORK(KINDSQ),
     *                              LENSQ,WORK(KEND4),LWRK4,
     *                              ISYMB,B,ISYMD,D)
                  END IF
C
C--------------------------------------------------
C Sum up S and U intermediates to get T3 BD amplitudes
C--------------------------------------------------
C
                  CALL CC3_T3BD(ISCKIJ,WORK(KSMAT),WORK(KSMAT3),
     *                          WORK(KUMAT),WORK(KUMAT3),WORK(KTMAT),
     *                          WORK(KINDSQ),LENSQ)
c
c                  call sum_pt3(work(ktmat),isymb,b,isymd,d,
c    *             1,work(kx3am),3)

C
                  IF (IPRINT .GT. 55) THEN
                     XUMAT = DDOT(NCKIJ(ISCKIJ),WORK(KTMAT),1,
     *                       WORK(KTMAT),1)
                     WRITE(LUPRI,*) 'Norm of TMAT^BD  ',XUMAT
                  ENDIF
C
C------------------------------------------------------
C Based on T3 BD amplitudes calculate W BD intermediates
C------------------------------------------------------
C
C------------------------------------------------------
C     Calculate the  term <mu3|[Y,T3]|HF> virtual contribution 
C     added in W^BD(KWMAT)
C------------------------------------------------------
C
                  CALL WBD_V(WORK(KTMAT),ISCKIJ,FOCKY,ISYFKY,
     *                            WORK(KWMAT), 
     *                            ISWMAT,WORK(KEND4),LWRK4)  
C
C------------------------------------------------------
C     Calculate the  term <mu3|[Y,T3]|HF> occupied contribution 
C     added in W^BD(KWMAT)
C------------------------------------------------------
C
                  CALL WBD_O(WORK(KTMAT),ISCKIJ,FOCKY,ISYFKY,
     *                            WORK(KWMAT),
     *                            ISWMAT,WORK(KEND4),LWRK4)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[[Y,T2],T2]|HF> 
C     added in W^BD(KWMAT)
C------------------------------------------------------
C
                  CALL WBD_T2(.FALSE.,B,ISYMB,D,ISYMD,
     *                        WORK(KT2AM),ISYMT2,WORK(KT2AM),ISYMT2,
     *                        FOCKY,ISYFKY,
     *                        WORK(KINDSQW),LENSQW,WORK(KWMAT), 
     *                        ISWMAT,WORK(KEND4),LWRK4) 
C
*                    call sum_pt3(work(kwmat),isymb,b,isymd,d,
*    *                            iswmat,work(kx3am),4)
C
C
                  IF (LCAUCHY .AND. ICAU.GE.0) THEN
C
                     DO MCAU = 0, NCAU
C
                        !Get offset for C2 vec for a given MCAU
                        KOFF2 = KC2 + MCAU*NT2SQ(ISYMBB)
C
                        !Calculate <mu3|[H^0,C2]|HF> cont to C3 cauchy vec
                        CALL WBD_GROUND(WORK(KOFF2),ISYMBB,WORK(KTMAT),
     *                               WORK(KTRVI0),WORK(KTRVI8),
     *                               WORK(KTROC0),ISINT2,WORK(KWMAT),
     *                               WORK(KEND4),LWRK4,
     *                               WORK(KINDSQW),LENSQW,
     *                               ISYMB,B,ISYMD,D)

C
                        !Get offset for C1-trans occ int for a given MCAU
                        KOFFOCC  = KOCCC1 + MCAU*NTRAOC(ISYMBB)
C
                        !Get offset for C1-trans virt int with fixed D for MCAU
                        KOFFINTD = KVIRDC1 + MCAU*NCKATR(ISYMBBD)  
                        !Get offset for C1-trans virt int with fixed B for MCAU
                        KOFFINTB = KVIRBC1 + MCAU*NCKATR(ISYMBBB)
C
                        
                        !Calculate <mu3|[[H^0,C1],T2^0]|HF>
                        CALL WBD_GROUND(WORK(KT2AM),ISYMT2,WORK(KTMAT),
     *                               WORK(KOFFINTD),WORK(KOFFINTB),
     *                               WORK(KOFFOCC),ISYMBB,WORK(KWMAT),
     *                               WORK(KEND4),LWRK4,
     *                               WORK(KINDSQW),LENSQW,
     *                               ISYMB,B,ISYMD,D)
C
*                    call sum_pt3(work(kwmat),isymb,b,isymd,d,
*    *                            iswmat,work(kx3am),4)

                        !Divide by the energy difference 
                        !and remove the forbidden elements
                        CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQY, 
     *                               ISWMAT,WORK(KWMAT), 
     *                               WORK(KDIAGW),WORK(KFOCKD)) 
C
                        CALL T3_FORBIDDEN(WORK(KWMAT),ISYFKY,ISYMB,B,
     *                                    ISYMD,D) 
C
                     END DO !MCAU
C
                     !trun the sign C_T  <-  -C_T
                     CALL DSCAL(NCKIJ(ISWMAT),-1.0D0,WORK(KWMAT),1)
C
                  END IF !LCAUCHY
C

*                    call sum_pt3(work(kwmat),isymb,b,isymd,d,
*    *                            iswmat,work(kx3am),4)

                  !Divide by the energy difference and
                  !remove the forbidden elements
                  CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQY, 
     *                         ISWMAT,WORK(KWMAT), 
     *                         WORK(KDIAGW),WORK(KFOCKD)) 
C
                  CALL T3_FORBIDDEN(WORK(KWMAT),ISYFKY,ISYMB,B,
     *                              ISYMD,D) 
C
*                    call sum_pt3(work(kwmat),isymb,b,isymd,d,
*    *                            iswmat,work(kx3am),4)
C
C------------------------------------------------------
C     Calculate the  term <mu1|[H,W^BD(3)]|HF> 
C     added in XI1EFF 
C------------------------------------------------------
C
                  CALL CC3_CY1(XI1EFF,ISYRES,WORK(KWMAT),ISWMAT,
     *                         WORK(KTMAT),WORK(KXIAJB),
     *                         ISINT1,WORK(KINDSQW),LENSQW,
     *                         WORK(KEND4),LWRK4,
     *                         ISYMB,B,ISYMD,D)

C
C------------------------------------------------------
C     Calculate the  term <mu2|[H,W^BD(3)]|HF> ( Fock matrix cont ) 
C     added in XI2EFF 
C------------------------------------------------------
C
                  CALL CC3_CY2F(XI2EFF,ISYRES,WORK(KWMAT),ISWMAT,
     *                          WORK(KTMAT), FOCK0,ISINT1,WORK(KINDSQW),
     *                          LENSQW,WORK(KEND4),LWRK4,
     *                          ISYMB,B,ISYMD,D)
c                 
C------------------------------------------------------
C     Calculate the  term <mu2|[H,W^BD(3)]|HF> ( Occupied  cont ) 
C     added in XI2EFF 
C------------------------------------------------------
C
                 CALL CC3_CY2O(XI2EFF,ISYRES,WORK(KWMAT),ISWMAT,
     *                          WORK(KTMAT),WORK(KTROC),
     *                          WORK(KTROC1),ISINT1,WORK(KEND4),LWRK4,
     *                          WORK(KINDSQW),LENSQW,
     *                          ISYMB,B,ISYMD,D)
C
C
C------------------------------------------------------
C     Calculate the  term <mu2|[H,W^BD(3)]|HF> ( Vccupied  cont ) 
C     added in XI2EFF 
C------------------------------------------------------
C
                  CALL CC3_CY2V(XI2EFF,ISYRES,WORK(KRBJIA),WORK(KWMAT),
     *                          ISWMAT,WORK(KTMAT),WORK(KTRVI),
     *                          WORK(KTRVI1),ISINT1,WORK(KEND4),LWRK4,
     *                          WORK(KINDSQW),LENSQW,
     *                          ISYMB,B,ISYMD,D)
C 
C--------------------------------------
C     Accumulate RMAT2 in XI2EFF
C--------------------------------------
C
                  IF (.NOT.(LCAUCHY .AND. (ICAU.GE.0))) THEN
                    CALL CC3_RACC(XI2EFF,WORK(KRMAT2),ISYMB,B,ISYRES)
                  END IF
c
               END DO 
            END DO
C
C--------------------------------------
C     Accumulate RMAT1 in XI2EFF
C--------------------------------------
C
            IF (.NOT.(LCAUCHY .AND. (ICAU.GE.0))) THEN
              CALL CC3_RACC(XI2EFF,WORK(KRMAT1),ISYMD,D,ISYRES)
            END IF
C
         END DO
      END DO
*      write(lupri,*) 'TETAXKL in CC3_XISD  : '
*      write(lupri,*) 'C3 in CC3_XISD, ncau = ',ncau
*      call print_pt3(work(kx3am),1,4)
C
C
C------------------------------------------------------
C     Accumulate RBJIA from <mu2|[H,W^BD(3)]|HF> ( Vccupied  cont ) 
C     in XI2EFF 
C------------------------------------------------------
C
      CALL CC3_RBJIA(XI2EFF,ISYRES,WORK(KRBJIA))
c
c     write(lupri,*)'XI1EFF in CC3_XISD'
c     call output(XI1EFF,1,nvir(1),1,nrhf(1),nvir(1),nrhf(1),1,lupri)
c
c     write(lupri,*)'XI2EFF in CC3_XISD'
c     call outpak(XI2EFF,NT1AM(isyres),1,lupri)
C
      IF (IPRINT .GT. 55) THEN
         XXI2EFF = DDOT(NT2AM(ISYRES),XI2EFF,1,XI2EFF,1)
         WRITE(LUPRI,*) 'Norm of XI2EFF final before added  ',XXI2EFF
      ENDIF
C
      DO I = 1,NT2AM(ISYRES)
         XI2EFF(I) = XI2EFF(I) + XI2(I)
      END DO
C
      IF (IPRINT .GT. 55) THEN
         XXI2EFF = DDOT(NT2AM(ISYRES),XI2EFF,1,XI2EFF,1)
         WRITE(LUPRI,*) 'Norm of XI2EFF final, xi2eff + xi2  ',XXI2EFF
      ENDIF
C
      DO I = 1,NT1AM(ISYRES)
         XI1EFF(I) = XI1EFF(I) + XI1(I)
      END DO

*     write(lupri,*)'xi1 at end of CC3_XISD'
*     call PRINT_MATAI(XI1,ISYRES)
C
*     write(lupri,*)'xi1eff at end of CC3_XISD'
*     call PRINT_MATAI(XI1eff,ISYRES)
C
      IF (IPRINT .GT. 55) THEN
         XXI1EFF = DDOT(NT1AM(ISYRES),XI1EFF,1,XI1EFF,1)
         WRITE(LUPRI,*) 'Norm of XI1EFF final, xi1eff + xi1  ',XXI1EFF
      ENDIF
C
C---------------------------------------------------------------------
C     It might have happened that 'CC3_CAUINT_V*' files have not been
C     deleted up there. Do it now!
C---------------------------------------------------------------------
C
      IF (LCAUCHY .AND. ICAU.GE.0) THEN
         CALL DELETE_FILES('CC3_CAUINT_V*')
      END IF
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CC3_XISD ')
C
      RETURN
C
      END
C
C  /* Deck cc3_t3bd */
      SUBROUTINE CC3_T3BD(ISSMAT,SMAT,SMAT2,
     *                           UMAT,UMAT2,TMAT,
     *                           INDSQ,LENSQ)
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C
C     Sum up S and U intermediates for fixed b and d and 
C     write on TMAT
C    
C     tmat^{bd}(ckij) = smat^{bd}(ckij) + umat^{bd}(ckij)
C    *                + smat^{db}(ckji) + umat^{db}(ckji)
C
C     ISSMAT is symmetry of S matrix (=symmetry of TMAT)
C     
C
      IMPLICIT NONE
C
#include <ccsdsym.h>
C
      INTEGER ISSMAT, LENGTH, LENSQ, INDSQ(LENSQ,6)
C
      REAL*8  SMAT(*), SMAT2(*)
      REAL*8  UMAT(*), UMAT2(*), TMAT(*)
C
      CALL QENTER('CC3_T3BD')
C
      LENGTH = NCKIJ(ISSMAT)

C
C---------------------------
C     Sum up intermediates
C---------------------------
C
      DO I = 1, LENGTH
         TMAT(I) = 
     *             SMAT(I)
     *           + UMAT(I)
     *           + SMAT2(INDSQ(I,3))
     *           + UMAT2(INDSQ(I,3))
      ENDDO
C
C---------------------
C     End.
C---------------------
C
      CALL QEXIT('CC3_T3BD')
C
      RETURN
      END


c  /* deck print_t3bd */
c
      subroutine print_t3bd(tmat,isymim,b,c,iopt)
c
c     Print t3^{BC}_{aikj} amplitudes which have been 
c     summed up (outside) to tmat array. (IOPT = 1)
c
c     Print W^{BC}_{aikj} "amplitudes" which have been 
c     summed up (outside) to wmat array. (IOPT = 2)
c
c     isymim is symmetry of (aikj)
c
c     OBS!  B and C are FIXED and coming from OUTSIDE
c
c
c     P. Jorgensen and F. Pawlowski, Spring 2002.
c
      IMPLICIT NONE
c
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      INTEGER ISYMIM, KOFF1, KOFF4, KOFF5, KOFF6, KOFF7
      INTEGER ISYMA, ISYMI, ISYMJ, ISYMK
      INTEGER KH, ISYIJK, ISYMJK, ISYMAI, ISYAIK
      INTEGER IOPT
C
      REAL*8  tmat(*)
C
      CALL QENTER('print_t3bd')
C
C----------------------------------------
C     Print the triples amplitudes.
C----------------------------------------
C
      DO ISYMA = 1, NSYM
C
         KOFF1 = 0
         DO KH = 1, ISYMA-1
           KOFF1 = KOFF1 + NVIR(KH)
         ENDDO
C
         ISYIJK = MULD2H(ISYMIM,ISYMA)
C
         DO ISYMI = 1, NSYM
C
            KOFF4 = 0
            DO KH = 1, ISYMI-1
              KOFF4 = KOFF4 + NRHF(KH)
            ENDDO
C
            ISYMJK = MULD2H(ISYIJK,ISYMI)
            ISYMAI = MULD2H(ISYMA,ISYMI)
C
            DO ISYMJ = 1, NSYM
C
               KOFF5 = 0
               DO KH = 1, ISYMJ-1
                 KOFF5 = KOFF5 + NRHF(KH)
               ENDDO
C
               ISYMK   = MULD2H(ISYMJK,ISYMJ)
               ISYAIK = MULD2H(ISYMAI,ISYMK)
C
               KOFF6 = 0
               DO KH = 1, ISYMK-1
                 KOFF6 = KOFF6 + NRHF(KH)
               ENDDO
C
                DO A = 1, NVIR(ISYMA)
                DO I = 1, NRHF(ISYMI)
                DO J = 1, NRHF(ISYMJ)
                DO K = 1, NRHF(ISYMK)
C
                KOFF7 = ISAIKJ(ISYAIK,ISYMJ)
     *                + NCKI(ISYAIK)*(J - 1)
     *                + ISAIK(ISYMAI,ISYMK) + NT1AM(ISYMAI)*(K-1)
     *                + IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A

C
        IF (ABS(TMAT(KOFF7)) 
     *                                    .GT. 1.0D-12) THEN
           IF (IOPT .EQ. 1) THEN
              write(lupri,1) 'T3BD(',a+koff1,',',b,',',
     *                               c,',',i+koff4,',',
     *                               j+koff5,',',k+koff6,') = ',
     *        TMAT(KOFF7) 
           ELSE IF (IOPT .EQ. 2) THEN
              write(lupri,1) 'WBD(',a+koff1,',',b,',',
     *                               c,',',i+koff4,',',
     *                               j+koff5,',',k+koff6,') = ',
     *        TMAT(KOFF7) 
           ELSE 
              call quit('Illegal value for IOPT in print_t3bd ')
           ENDIF
C
        ENDIF
C
                ENDDO ! K
                ENDDO ! J
                ENDDO ! I
                ENDDO ! A
                ENDDO ! ISYMJ
                ENDDO ! ISYMI
C
      ENDDO           ! ISYMA
C
      CALL QEXIT('print_t3bd')
C
    1 FORMAT(1X,A6,I3,A1,I3,A1,I3,A1,I3,A1,I3,A1,I3,A4,E20.10)
      RETURN
      END
C
C  /* Deck wbd_v */
      SUBROUTINE WBD_V(TMAT,ISTMAT,FOCKY,ISYFKY, 
     *                 WMAT,ISWMAT,WRK,LWRK)
C 
C WBD(a,i,k,j) = WBD(a,i,k,j) + sum (f)  focky(a,f)*tmatBD(f,i,k,j) 
C
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C

C
      IMPLICIT NONE
C
      INTEGER LWRK, KFCAF, KEND0, LWRK0, KOFF1, KOFF2 
      INTEGER NF, KOFFY, KOFFT, KOFFW 
      INTEGER ISTMAT, ISYFKY, ISWMAT, ISFIKJ, ISYFIK 
      INTEGER ISYMA, ISYAI, ISYAIK, NA
      INTEGER ISYMF, ISYMJ, ISYMK, ISYMI, ISYFI
C
      REAL*8  TMAT(*), FOCKY(*), WMAT(*), WRK(*)
      REAL*8  HALF, ONE
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0)
C
      CALL QENTER('WBD_V')
C
C
C RESORT VIR-VIR  FOCKY ELEMENTS (A,F)
C
C
      KFCAF  = 1
      KEND0  = KFCAF + NMATAB(ISYFKY)
      LWRK0  = LWRK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK0
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in WBD_V')
      END IF
C
      DO ISYMF = 1,NSYM
         ISYMA = MULD2H(ISYMF,ISYFKY)
         DO F = 1,NVIR(ISYMF)
            KOFF1 = IFCVIR(ISYMA,ISYMF) + NORB(ISYMA)*(F - 1)
     *                                  + NRHF(ISYMA) + 1
            KOFF2 = KFCAF + IMATAB(ISYMA,ISYMF) + NVIR(ISYMA)*(F - 1) 
            CALL DCOPY(NVIR(ISYMA),FOCKY(KOFF1),1,WRK(KOFF2),1)
         END DO
      END DO
C
C CARRY OUT MATRIX MULTIPLICATION
C WBD(a,i,k,j) = WBD(a,i,k,j) + sum (f)  focky(a,f)*tmatBD(f,i,k,j) 
C 
      ISFIKJ = ISTMAT
      DO ISYMJ = 1,NSYM
         ISYFIK =MULD2H(ISYMJ,ISFIKJ)
         DO J   = 1,NRHF(ISYMJ)
            DO ISYMK = 1,NSYM
               ISYFI = MULD2H(ISYMK,ISYFIK)
               DO K  = 1,NRHF(ISYMK)
                  DO ISYMI = 1,NSYM
                     ISYMF = MULD2H(ISYFI,ISYMI)   
                     ISYMA = MULD2H(ISYFKY,ISYMF)
                     ISYAIK = MULD2H(ISWMAT,ISYMJ)
                     ISYAI = MULD2H(ISYAIK,ISYMK)
                     NA    = MAX(1,NVIR(ISYMA))
                     NF    = MAX(1,NVIR(ISYMF))
                     KOFFY = KFCAF + IMATAB(ISYMA,ISYMF) 
                     KOFFT = ISAIKJ(ISYFIK,ISYMJ)
     *                      + NCKI(ISYFIK)*(J-1)
     *                      + ISAIK(ISYFI,ISYMK) 
     *                      + NT1AM(ISYFI)*(K-1)
     *                      + IT1AM(ISYMF,ISYMI) + 1
                     KOFFW = ISAIKJ(ISYAIK,ISYMJ)
     *                      + NCKI(ISYAIK)*(J-1)
     *                      + ISAIK(ISYAI,ISYMK) 
     *                      + NT1AM(ISYAI)*(K-1)
     *                      + IT1AM(ISYMA,ISYMI) + 1
C
C SYMMETRY BETWEEN BJ AND CK INTRODUCE A FACTOR TWO
C DENOTE t3 IS CALCULATED WITH NEGATIVE SIGN
C
                     CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI), 
     *                          NVIR(ISYMF),-ONE,WRK(KOFFY),NA,
     *                          TMAT(KOFFT),NF,ONE,WMAT(KOFFW),NA) 
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('WBD_V')
C
      RETURN
      END
C
C  /* Deck wbd_o */
      SUBROUTINE WBD_O(TMAT,ISTMAT,FOCKY,ISYFKY,
     *                 WMAT,ISWMAT,WRK,LWRK)
C 
C WBD(a,i,k,j) = WBD(a,i,k,j) - sum (l) tmatBD(a,l,k,j)*focky(l,i)
C
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C

      IMPLICIT NONE
C
      INTEGER LWRK, KFCLI, KEND0, LWRK0, KOFF1, KOFF2 
      INTEGER NL, KOFFY, KOFFT, KOFFW 
      INTEGER ISTMAT, ISYFKY, ISWMAT, ISALKJ
      INTEGER ISYMA, ISYAI, ISYAIK, ISYALK, ISYAL, NA
      INTEGER ISYMJ, ISYMK, ISYMI, ISYML, ISYFI
C
      REAL*8  TMAT(*), FOCKY(*), WMAT(*), WRK(*)
      REAL*8  MHALF, ONE 
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (MHALF = -0.5D0, ONE = 1.0D0)
C
      CALL QENTER('WBD_O')
C

C
C RESORT OCC-OCC  FOCKY ELEMENTS (L,I)
C
C
      KFCLI  = 1
      KEND0  = KFCLI + NMATIJ(ISYFKY)
      LWRK0  = LWRK  - KEND0
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK0
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in WBD_O')
      END IF
C   
      DO ISYMI = 1,NSYM
         ISYML = MULD2H(ISYMI,ISYFKY)
         DO I = 1,NRHF(ISYMI)
             KOFF1 = IFCRHF(ISYML,ISYMI) + NORB(ISYML)*(I - 1) + 1
             KOFF2 = KFCLI + IMATIJ(ISYML,ISYMI) + NRHF(ISYML)*(I - 1)
             CALL DCOPY(NRHF(ISYML),FOCKY(KOFF1),1,WRK(KOFF2),1) 
         END DO
      END DO
C
C CARRY OUT MATRIX MULTIPLICATION
C WBD(a,i,k,j) = WBD(a,i,k,j) - sum (l) tmatBD(a,l,k,j)*focky(l,i)
C 
      ISALKJ = ISTMAT
      DO ISYMJ = 1,NSYM
         ISYALK =MULD2H(ISYMJ,ISALKJ)
         DO J = 1,NRHF(ISYMJ)
            DO ISYMK = 1,NSYM
               ISYAL = MULD2H(ISYMK,ISYALK)
               DO K  = 1,NRHF(ISYMK)
                  DO ISYML = 1,NSYM
                     ISYMA = MULD2H(ISYAL,ISYML)   
                     ISYMI = MULD2H(ISYFKY,ISYML)
                     ISYAIK = MULD2H(ISWMAT,ISYMJ)
                     ISYAI = MULD2H(ISYAIK,ISYMK)
                     NA    = MAX(1,NVIR(ISYMA))
                     NL    = MAX(1,NRHF(ISYML))
                     KOFFY = KFCLI + IMATIJ(ISYML,ISYMI) 
                     KOFFT = ISAIKJ(ISYALK,ISYMJ)
     *                      + NCKI(ISYALK)*(J-1)
     *                      + ISAIK(ISYAL,ISYMK) 
     *                      + NT1AM(ISYAL)*(K-1)
     *                      + IT1AM(ISYMA,ISYML) + 1
                     KOFFW = ISAIKJ(ISYAIK,ISYMJ)
     *                      + NCKI(ISYAIK)*(J-1)
     *                      + ISAIK(ISYAI,ISYMK) 
     *                      + NT1AM(ISYAI)*(K-1)
     *                      + IT1AM(ISYMA,ISYMI) + 1
C
C SYMMETRY BETWEEN BJ AND CK INTRODUCE A FACTOR TWO
C DENOTE t3 IS CALCULATED WITH NEGATIVE SIGN
C
                     CALL DGEMM('N','N',NVIR(ISYMA),NRHF(ISYMI),
     *                          NRHF(ISYML),ONE,TMAT(KOFFT),NA,
     *                          WRK(KOFFY),NL,ONE,WMAT(KOFFW),NA)
                  END DO
               END DO
            END DO
         END DO
      END DO
C
      CALL QEXIT('WBD_O')
C
      RETURN
      END
C  /* Deck wbd_t2 */
      SUBROUTINE WBD_T2(T2XNET2Z,B,ISYMB,D,ISYMD,
     *                 T2TPX,ISYMT2X,T2TPZ,ISYMT2Z,
     *                 FOCKY,ISYFKY,
     *                 INDSQ,LENSQ,WMAT,ISWMAT,WRK,LWRK)
*
**********************************************************************
*
* Calculate the term <mu3|[[Y,T2X],T2Z]|HF> as the contribution to 
* WBC intermmediate.
*
* If T2XNET2Z = .true. then calculate:
*
* what is returned is:
*
* WBC(a,i,k,j) = WBC(a,i,k,j) -
*  P(bj,ck) sum(l,d) focky(l,d)*( t2X(ai,cl)*t2Z(bj,dk)+t2Z(ai,cl)*t2X(bj,dk))
* 
* ELSE
* 
* <mu3|[[Y,T2],T2]|HF> = P(aibjck) 2 sum(l,d) focky(l,d)*(t2X(ai,cl)*t2Z(bj,dk))
*
* what is returned by the routine is WBC contribution corresponding to 
* (ai,bj,ck) + (ai,ck,bj) permutation.
*
* Totally is returned: (note that 2 factor is absorbed usually by 1/2
*                       in the formula; for example in T3X  !!!)
*
* WBC(a,i,k,j) = WBC(a,i,k,j) -
*   P(bj,ck) sum(l,d) focky(l,d)*( t2X(ai,cl)*t2Z(bj,dk) )
*
*
* F. Pawlowski, P. Jorgensen, Aarhus 14-May-2003.
*
* This is just a little driver which calls WBD_T2_1. 
* 
**********************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccsdsym.h>
C
      LOGICAL T2XNET2Z
C
      INTEGER LENSQ, INDSQ(LENSQ,6), LWRK
      INTEGER ISYMB, ISYMD, ISYMT2X, ISYFKY, ISWMAT, ISYMT2Z
C
      REAL*8  T2TPX(*), FOCKY(*), WMAT(*), WRK(*), T2TPZ(*)
      REAL*8  FAC,HALF,ONE
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0 )
C
      CALL QENTER('WBDT2')
C
      IF (T2XNET2Z) THEN
         FAC = ONE
      ELSE
         FAC = ONE
      END IF
C
      CALL WBD_T2_1(FAC,B,ISYMB,D,ISYMD,
     *              T2TPX,ISYMT2X,T2TPZ,ISYMT2Z,
     *              FOCKY,ISYFKY,
     *              INDSQ,LENSQ,WMAT,ISWMAT,WRK,LWRK)
C
      IF (T2XNET2Z) THEN
C
         CALL WBD_T2_1(FAC,B,ISYMB,D,ISYMD,
     *                 T2TPZ,ISYMT2Z,T2TPX,ISYMT2X,
     *                 FOCKY,ISYFKY,
     *                 INDSQ,LENSQ,WMAT,ISWMAT,WRK,LWRK)      
C
      END IF
C
      CALL QEXIT('WBDT2')
C
      RETURN
      END
C  /* Deck wbd_t2_1 */
      SUBROUTINE WBD_T2_1(FAC,B,ISYMB,D,ISYMD,T2TPX,ISYMT2X,T2TPZ,
     *                 ISYMT2Z,
     *                 FOCKY,ISYFKY,
     *                 INDSQ,LENSQ,WMAT,ISWMAT,WRK,LWRK)
C 
C WBD(a,i,k,j) = WBD(a,i,k,j) -
C      sum (f,l) focky(l,f)*( t2X(ai,dl)*t2Z(fk,bj) + t2X(ai,bl)*t2Z(fj,dk) )
C
C
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C
C     Modified to handle two different symmetries of T2, 14-May-2003.

      IMPLICIT NONE
C
#include <priunit.h>
#include <dummy.h>
#include <iratdef.h>
#include <ccsdsym.h>
#include <inftap.h>
#include <ccinftap.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      INTEGER LENSQ
      INTEGER INDSQ(LENSQ,6)
      INTEGER LWRK, KFCLF, KEND0, LWRK0, KOFF1, KOFF2, KTB, KEND1, LWRK1
      INTEGER NL, NF, KOFFY, KOFFT2, KOFFT, KOFFW, KTD, KW
      INTEGER ISYMB, ISYMD, ISYMT2X, ISYFKY, ISWMAT 
      INTEGER ISYAIL, ISYAI, ISYAIK, NA, NAI, LENGTH
      INTEGER ISYMF, ISYML, ISYFKJ, ISYTB, ISYMJ, ISYFK, ISYMK, ISYLK
      INTEGER ISYFJK, ISYTD, ISYLJ, ISYFJ, ISYAIJ 
      INTEGER ISYMT2Z
C
      REAL*8  T2TPX(*), FOCKY(*), WMAT(*), WRK(*)
      REAL*8  HALF, ONE, ZERO, FAC 
      REAL*8  T2TPZ(*)
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0, ZERO = 0.0D0)
C
      CALL QENTER('WT21')
C
C
C RESORT VIR-OCC  FOCKY ELEMENTS (l,f)
C
C
      KW = 1
      KFCLF = KW + NCKIJ(ISWMAT)
      KEND0  = KFCLF + NT1AM(ISYFKY)
      LWRK0  = LWRK  - KEND0
      CALL DZERO(WRK(KW),NCKIJ(ISWMAT))
C
      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK0
         WRITE(LUPRI,*) 'Memory needed    : ',KEND0
         CALL QUIT('Insufficient space in WBD_T2_1 (1)')
      END IF
C
      DO ISYMF = 1,NSYM
         ISYML = MULD2H(ISYMF,ISYFKY)
         DO L = 1,NRHF(ISYML)
            DO F = 1,NVIR(ISYMF)
               KOFF1 = IFCVIR(ISYML,ISYMF) + NORB(ISYML)*(F - 1) + L
               KOFF2 = KFCLF +  IT1AMT(ISYML,ISYMF) 
     *               + NRHF(ISYML)*(F - 1) + L -1 
C
                  WRK(KOFF2) = FOCKY(KOFF1)
C
            END DO
         END DO
      END DO
C
C    calculate first t2 contribution to W matrix
C
C construct tZB(l,k,j) = sum (f) focky(l,f)*t2tpZ(f,k,j,B)
C
      ISYFKJ   = MULD2H(ISYMT2Z,ISYMB) 
      ISYTB    = MULD2H(ISYFKY,ISYFKJ) 
      KTB      = KEND0
      KEND1    = KTB  + NMAIJK(ISYTB)
      LWRK1    = LWRK  - KEND1
C
      CALL DZERO(WRK(KTB),NMAIJK(ISYTB))
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK1
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBD_T2_1 (2)')
      END IF
C
      DO ISYMJ = 1,NSYM
         ISYFK  = MULD2H(ISYFKJ,ISYMJ) 
         DO J  = 1,NRHF(ISYMJ)
            DO ISYMK = 1,NSYM
               ISYMF = MULD2H(ISYFK,ISYMK)
               ISYML = MULD2H(ISYFKY,ISYMF)
               ISYLK  = MULD2H(ISYML,ISYMK)
               NL = MAX(1,NRHF(ISYML))
               NF = MAX(1,NVIR(ISYMF))
               KOFFY  = KFCLF + IT1AMT(ISYML,ISYMF)  
               KOFFT2 = IT2SP(ISYFKJ,ISYMB) + NCKI(ISYFKJ)*(B-1)
     *                 + ISAIK(ISYFK,ISYMJ) + NT1AM(ISYFK)*(J-1)
     *                 + IT1AM(ISYMF,ISYMK) + 1
               KOFFT =  KTB + IMAIJK(ISYLK,ISYMJ) + NMATIJ(ISYLK)*(J-1) 
     *               + IMATIJ(ISYML,ISYMK)
C
               CALL DGEMM('N','N',NRHF(ISYML),NRHF(ISYMK),
     *                    NVIR(ISYMF),ONE,WRK(KOFFY),NL,
     *                    T2TPZ(KOFFT2),NF,ONE,WRK(KOFFT),NL)
C
            END DO
         END DO
      END DO
C
C      WBD(a,i,k,j) = WBD(a,i,k,j) -
C                        sum (f,l) focky(l,f)* t2X(ai,Dl)*t2Z(fk,Bj) 
C                   = WBD(a,i,k,j) - 
C                        sum(l) t2tpX(a,i,l,D) * tZB(l,k,j)
C
      ISYAIL = MULD2H(ISYMT2X,ISYMD)
      DO ISYMJ = 1,NSYM
         ISYLK  = MULD2H(ISYTB,ISYMJ)
         DO J  = 1,NRHF(ISYMJ)
            DO ISYMK = 1,NSYM
               ISYML = MULD2H(ISYLK,ISYMK)
               ISYAI = MULD2H(ISYAIL,ISYML) 
               ISYAIK = MULD2H(ISYAI,ISYMK)
               NAI = MAX(1,NT1AM(ISYAI))
               NL = MAX(1,NRHF(ISYML)) 
               KOFFT2 = IT2SP(ISYAIL,ISYMD) + NCKI(ISYAIL)*(D-1)
     *                 + ISAIK(ISYAI,ISYML) + 1
               KOFFT =  KTB + IMAIJK(ISYLK,ISYMJ) + NMATIJ(ISYLK)*(J-1) 
     *                      + IMATIJ(ISYML,ISYMK)
               KOFFW  = ISAIKJ(ISYAIK,ISYMJ) + NCKI(ISYAIK)*(J-1)
     *                 + ISAIK(ISYAI,ISYMK) + 1
               CALL DGEMM('N','N',NT1AM(ISYAI),NRHF(ISYMK),
     *                    NRHF(ISYML),-FAC,T2TPX(KOFFT2),NAI,
     *                    WRK(KOFFT),NL,ONE,WMAT(KOFFW),NAI)

C
            END DO
         END DO
      END DO
C
C    calculate second t2 contribution to W matrix
C
C construct tZD(l,j,k) = sum (f) focky(l,f)*t2tpZ(f,j,k,D)
C
      ISYFJK   = MULD2H(ISYMT2Z,ISYMD) 
      ISYTD    = MULD2H(ISYFKY,ISYFJK) 
      KTD      = KEND0
      KEND1    = KTD  + NMAIJK(ISYTD)
      LWRK1    = LWRK  - KEND1
C
      CALL DZERO(WRK(KTD),NMAIJK(ISYTD))
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWRK1
         WRITE(LUPRI,*) 'Memory needed    : ',KEND1
         CALL QUIT('Insufficient space in WBD_T2_1 (3)')
      END IF
C

      DO ISYMK = 1,NSYM
         ISYFJ  = MULD2H(ISYFJK,ISYMK) 
         DO K  = 1,NRHF(ISYMK)
            DO ISYMJ = 1,NSYM
               ISYMF = MULD2H(ISYFJ,ISYMJ)
               ISYML = MULD2H(ISYFKY,ISYMF)
               ISYLJ  = MULD2H(ISYML,ISYMJ)
               NL = MAX(1,NRHF(ISYML))
               NF = MAX(1,NVIR(ISYMF))
               KOFFY  = KFCLF + IT1AMT(ISYML,ISYMF)  
               KOFFT2 = IT2SP(ISYFJK,ISYMD) + NCKI(ISYFJK)*(D-1)
     *                 + ISAIK(ISYFJ,ISYMK) + NT1AM(ISYFJ)*(K-1)
     *                 + IT1AM(ISYMF,ISYMJ) + 1
               KOFFT =  KTD + IMAIJK(ISYLJ,ISYMK) + NMATIJ(ISYLJ)*(K-1) 
     *                      + IMATIJ(ISYML,ISYMJ)
               CALL DGEMM('N','N',NRHF(ISYML),NRHF(ISYMJ),
     *                    NVIR(ISYMF),ONE,WRK(KOFFY),NL,
     *                    T2TPZ(KOFFT2),NF,ONE,WRK(KOFFT),NL)
C
            END DO
         END DO
      END DO
C
C      WBD(a,i,k,j) = WBD(a,i,k,j) -
C                        sum (f,l) focky(l,f)*t2X(ai,Bl)*t2Z(fj,Dk) )
C                   = WBD(a,i,k,j) -
C                        sum(l) t2tpX(a,i,l,B) * tZD(l,j,k)
C
      ISYAIL = MULD2H(ISYMT2X,ISYMB)
      DO ISYMK = 1,NSYM
         ISYLJ  = MULD2H(ISYTD,ISYMK)
         DO K  = 1,NRHF(ISYMK)
            DO ISYMJ = 1,NSYM
               ISYML = MULD2H(ISYLJ,ISYMJ)
               ISYAI = MULD2H(ISYAIL,ISYML) 
               ISYAIJ = MULD2H(ISYAI,ISYMJ) 
               NAI = MAX(1,NT1AM(ISYAI))
               NL = MAX(1,NRHF(ISYML)) 
               KOFFT2 = IT2SP(ISYAIL,ISYMB) + NCKI(ISYAIL)*(B-1)
     *                 + ISAIK(ISYAI,ISYML) + 1
               KOFFT =  KTD + IMAIJK(ISYLJ,ISYMK) + NMATIJ(ISYLJ)*(K-1) 
     *                      + IMATIJ(ISYML,ISYMJ)
               KOFFW  = KW  + ISAIKJ(ISYAIJ,ISYMK) + NCKI(ISYAIJ)*(K-1)
     *                 + ISAIK(ISYAI,ISYMJ) 
               CALL DGEMM('N','N',NT1AM(ISYAI),NRHF(ISYMJ),
     *                    NRHF(ISYML),-FAC,T2TPX(KOFFT2),NAI,
     *                    WRK(KOFFT),NL,ONE,WRK(KOFFW),NAI)

C
            END DO
         END DO
      END DO
C
C     change order aijk to aikj
C
      DO I = 1,NCKIJ(ISWMAT)
         WMAT(I) = WMAT(I) + WRK(INDSQ(I,3))
      END DO
C
      CALL QEXIT('WT21')
C
      RETURN
      END


C  /* Deck cc3_cY1 */
      SUBROUTINE CC3_CY1(OMEGA1,ISOMG1,WMAT,ISWMAT,TMAT,
     *                   XIAJB,ISYINT,INDSQ,LENSQ,WORK,LWORK,
     *                   ISYMIB,IB,ISYMID,ID)
C
C     P. Jorgensen and F. Pawlowski.         Spring 2002
C
      IMPLICIT NONE
C
C     Calculate Omega1 
C
C
C     General symmetry: ISWMAT is symmetry of WMAT and TMAT intermdiates.  
C                       (exclusive isymb,isymd) 
C                       ISYINT = XIAJB
C                       ISOMG1 = ISWMAT*ISYINT
C 
C        Omega1(ai)   
C
C        = sum_{bjck} ( W^{abc}_{ijk}  -  W^{abc}_{kji} ) * L_{jbkc}
C  (1):
C        = sum_{jk} ( W^{BC}_{aikj} - W^{BC}_{akij} ) * L_{CkBj} 
C  (2):
C        + sum_{bjk} ( W^{CA}_{bjik} - W^{AC}_{bjki} ) * L_{bjCk} 
C  (3):
C        + sum_{cjk} ( W^{AB}_{ckji} - W^{AB}_{cijk} ) * L_{ckBj} 
C
C
C
      INTEGER ISOMG1, ISWMAT, ISYINT, LENSQ, LWORK
      INTEGER ISYMIB, IB, ISYMID, ID, INDEX, ISYBJK
      INTEGER ISYBC, ISYMB, ISYMC, ISYKJ, ISYAI, ISYMJ, ISYBJ
      INTEGER ISYMK, NKJ, NBJ, NCK, NCKBJ, NTOTAI, KOFF1, NBJK
      INTEGER ISYMA, ISYMJK, ISYMI, ISYCK, NBJCK, NTOTA, NTOBJK
      INTEGER KOFF2, ISYCKJ, NCKJ, NTOCKJ, KSTART, KEND, LEND 
      INTEGER INDSQ(LENSQ,6)
      REAL*8  OMEGA1(*), WMAT(*), TMAT(*), XIAJB(*), WORK(*)
      REAL*8  ONE, ZERO, TWO

C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C 
      CALL QENTER('CC3_CY1')
C
C----------------------------------
C     Calculate (1):
C         sum_{jk} ( W^{BC}_{aikj} - W^{BC}_{akij} ) * L_{CkBj}
C----------------------------------
C
C       
C     maximum work space
C
      KSTART  = 1
      KEND    = KSTART + NRHFT*NRHFT*NVIRT
      LEND    = LWORK - KEND
C
      IF (LEND .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK 
         WRITE(LUPRI,*) 'Memory needed    : ',KEND
         CALL QUIT('Insufficient space in CC3_CY1')
      END IF

      B = IB
      C = ID
      ISYMB = ISYMIB
      ISYMC = ISYMID
C
C      T(aikj) =  W^{BC}_{aikj} - W^{BC}_{akij} )
C
      DO I = 1,LENSQ
         TMAT(I) =   WMAT(I)
     *             - WMAT(INDSQ(I,1))
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENSQ,WORK(KEND),TMAT,INDSQ(1,6))
         CALL DCOPY(LENSQ,WORK(KEND),1,TMAT,1)
      ENDIF
C
      ISYBC  = MULD2H(ISYMB,ISYMC)
      ISYKJ = MULD2H(ISYBC,ISYINT)
      ISYAI = ISOMG1
C
C        Construct LCB(k,j) = L(Ck,Bj)
C        ---------------------------
C
      DO 400 ISYMJ = 1,NSYM
C
         ISYMK = MULD2H(ISYMJ,ISYKJ)
         ISYCK = MULD2H(ISYMC,ISYMK)
         ISYBJ = MULD2H(ISYMB,ISYMJ)
C
         DO 410 J = 1,NRHF(ISYMJ)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            DO 420 K = 1,NRHF(ISYMK)
C
               NKJ = IMATIJ(ISYMK,ISYMJ)+ NRHF(ISYMK)*(J - 1) + K
               NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
               NCKBJ = IT2AM(ISYCK,ISYBJ) + INDEX(NCK,NBJ)
C
               WORK(NKJ) = XIAJB(NCKBJ)
C
  420       CONTINUE
  410    CONTINUE
  400 CONTINUE
C
         NTOTAI = MAX(NT1AM(ISYAI),1)
C
         KOFF1 = ISAIKL(ISYAI,ISYKJ) + 1
C
         CALL DGEMV('N',NT1AM(ISYAI),NMATIJ(ISYKJ),-ONE,TMAT(KOFF1),
     *              NTOTAI,WORK,1,ONE,OMEGA1,1)
C
C
C  Calculate  (2):
C
C  Omega1(ai) = Omega1(ai)  
C        + sum_{jk} ( W^{AC}_{bjki} - W^{AC}_{bjik} ) * L_{bjCk} 
C
C        + sum_{bjk} ( W^{CA}_{bjik} - W^{CA}_{bjki} ) * L_{bjCk} 
C
      C = IB
      A = ID
c     write(lupri,*)' start term 2 a,id',a,id
C
      ISYMC = ISYMIB
      ISYMA = ISYMID
C
C
C      T(bjki) =  W^{CA}_{bjik} - W^{AC}_{bjki} )  
C
      DO I = 1,LENSQ
C
         TMAT(I) = - WMAT(I) 
     *             + WMAT(INDSQ(I,3))
C
      ENDDO
C
C  Construct LC(bj,k) = L(bj,Ck)
C        ---------------------------
C
      ISYBJK = MULD2H(ISYINT,ISYMC) 
      ISYMI  = MULD2H(ISYBJK,ISWMAT)
      DO 100 ISYMK = 1,NSYM
C
         ISYCK = MULD2H(ISYMC,ISYMK)
         ISYBJ = MULD2H(ISYMK,ISYBJK) 
C
         DO 110 K = 1,NRHF(ISYMK)
C
            NCK = IT1AM(ISYMC,ISYMK) + NVIR(ISYMC)*(K - 1) + C
C
            DO 120 NBJ = 1,NT1AM(ISYBJ)
C
               NBJCK = IT2AM(ISYBJ,ISYCK) + INDEX(NBJ,NCK)
               NBJK  = ICKI(ISYBJ,ISYMK)
     *               + NT1AM(ISYBJ)*(K - 1) + NBJ
C
               WORK(NBJK) = XIAJB(NBJCK)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      NTOTA  = MAX(NVIR(ISYMA),1)
      NTOBJK = MAX(NCKI(ISYBJK),1)
C
      KOFF1 = ISAIKJ(ISYBJK,ISYMI) + 1
      KOFF2 = IT1AM(ISYMA,ISYMI) + A
C
      CALL DGEMV('T',NCKI(ISYBJK),NRHF(ISYMI),-ONE,TMAT(KOFF1),
     *           NTOBJK,WORK,1,ONE,OMEGA1(KOFF2),NTOTA)
C
C  Calculate  (3):
C
C  Omega1(ai) = Omega1(ai)  
C        + sum_{jk} ( W^{AB}_{ckji} - W^{AB}_{cijk} ) * L_{ckBj}
C
C
      A = IB
      B = ID
C
      ISYMA = ISYMIB
      ISYMB = ISYMID
C
C
C      T(ckji) =  W^{AB}_{ckji} - W^{AB}_{cijk}
C
      DO I = 1,LENSQ
C
         TMAT(I) =   WMAT(I) 
     *             - WMAT(INDSQ(I,5))
C
      ENDDO
C
C  Construct LB(ck,j) = L(ck,Bj)    
C  ----------------------------
C
      ISYCKJ = MULD2H(ISYINT,ISYMB) 
      ISYMI  = MULD2H(ISYCKJ,ISWMAT)
      DO 200 ISYMJ = 1,NSYM
C
         ISYBJ = MULD2H(ISYMB,ISYMJ)
         ISYCK = MULD2H(ISYMJ,ISYCKJ) 
C
         DO 210 J = 1,NRHF(ISYMJ)
C
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
C
            DO 220 NCK = 1,NT1AM(ISYCK)
C
               NCKBJ = IT2AM(ISYCK,ISYBJ) + INDEX(NCK,NBJ)
               NCKJ  = ICKI(ISYCK,ISYMJ)
     *               + NT1AM(ISYCK)*(J - 1) + NCK
C
               WORK(NCKJ) = XIAJB(NCKBJ)
C
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      NTOTA  = MAX(NVIR(ISYMA),1)
      NTOCKJ = MAX(NCKI(ISYCKJ),1)
C
      KOFF1 = ISAIKJ(ISYCKJ,ISYMI) + 1 
      KOFF2 = IT1AM(ISYMA,ISYMI) + A
C
      CALL DGEMV('T',NCKI(ISYCKJ),NRHF(ISYMI),-ONE,TMAT(KOFF1),
     *           NTOCKJ,WORK,1,ONE,OMEGA1(KOFF2),NTOTA)
C
      CALL QEXIT('CC3_CY1') 
C
      RETURN
      END
C
C  /* Deck wbd_dia */ 
      SUBROUTINE WBD_DIA(B,ISYMB,C,ISYMC,FREQY, 
     *                   ISWMAT,WMAT,DIAG,FOCKD)
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C
C     Divide by orbital energies
C     
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER ISYMB,ISYMC,ISWMAT,NB,NC 
C
      REAL*8  WSMAT,DDOT,FREQY,EPSIBC 
      REAL*8  WMAT(*), DIAG(*), FOCKD(*) 
C
      CALL QENTER('WBD_DIA')
C
C-----------------------------------------
C     Divide by the Fock matrix diagonals.
C-----------------------------------------
C
      NB = IORB(ISYMB) + NRHF(ISYMB) + B
      NC = IORB(ISYMC) + NRHF(ISYMC) + C
C
      EPSIBC = FOCKD(NB) + FOCKD(NC) - FREQY
C
      DO 500 L = 1,NCKIJ(ISWMAT) 
C
         WMAT(L) = WMAT(L)/(DIAG(L) + EPSIBC)
C
  500 CONTINUE
C
c     IF (IPRINT .GT. 55) THEN
c        WSMAT = DDOT(NCKIJ(ISWMAT),WMAT,1,WMAT,1)
c        WRITE(LUPRI,*) 'In wbd_dia: 5. Norm of WMAT ',WSMAT
c     ENDIF
C
C---------------------
C     End.
C---------------------
C
      CALL QEXIT('WBD_DIA')
C
      RETURN
      END

C  /* Deck wbdtest */ 
      SUBROUTINE WBDtest(omega1,IB,IC,
     *                   TMAT,xiajb,iopt,work,lwrk)
C
C     Written by P. Jorgensen and F. Pawlowski, Spring 2002.
C
C    noddy code for the three cont to omega1 
C     
C
      IMPLICIT NONE
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      INTEGER iopt,lwrk,ai,nck,nbj,kckbj,kaikj,kbjki,index,ib,ic
      integer kckji
C
      REAL*8  tMAT(*), xiajb(*), work(*), omega1(*) 
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
c     write(lupri,*)' entering wbdtest iop',iopt
      write(lupri,*)' tmat i test'
C     call output(tmat,1,nvirt*nrhft,nrhft,nrhft,1,1,
C    *            nvirt*nrhft,nrhft,nrhft,1,lupri)
      do i = 1,nvirt*nrhft*nrhft*nrhft
         write(lupri,*) i , tmat(I)
      end do

C
C-----------------------------------------
      call dcopy(nt1amx,omega1,1,work,1) 
      write(lupri,*)' xi1eff wbdtest before cont is added '
      do i = 1,nt1am(1)
         write(lupri,*) i , work(I)
      end do
C     call dzero(work,nt1amx)
C
      if (iopt.eq.1) then
      b = ib
      c = ic
      write(lupri,*)' b,c third contribution',b,c
      DO AI = 1,NT1AMX
        DO K = 1,NRHFT
          DO J= 1,NRHFT
                nck = nvirt*(k-1) + c
                nbj = nvirt*(j-1) + b
                kckbj = index(nck,nbj) 
                
                kaikj = nvirt*nrhft*nrhft*(j-1)
     *                + nvirt*nrhft*(k-1) + ai
c           write(lupri,*)'nck,nbj,kckbj, xckbj'
c           write(lupri,*)nck,nbj,kckbj,XIAJB(KCKBJ) 
c           write(lupri,*)'ai,j,k,kaikj,tmat(kaikj)'
c           write(lupri,*)ai,j,k,kaikj,tmat(kaikj)
c           write(lupri,*)'WORK(AI),XIAJB(KCKBJ)*tmat(kaikj'
c           write(lupri,*)WORK(AI),XIAJB(KCKBJ)*tmat(kaikj)
            WORK(AI) = WORK(AI) + XIAJB(KCKBJ)*tmat(kaikj)
c           write(lupri,*)'WORK(AI)',WORK(AI)
          end do
        end do
      end do

      write(lupri,*)' xi1eff wbdtest first cont b,c',b,c
      do i = 1,nt1am(1)
         write(lupri,*) i , work(I)
      end do
      end if
cc
cc
      if (iopt.eq.2) then
      A = ic
      C = ib
      write(lupri,*)' a,c third contribution',a,c
      do i = 1,nrhft
       DO B = 1,nvirt
       ai = nvirt*(i-1)+a
        DO K = 1,NRHFT
          DO J= 1,NRHFT
                nck = nvirt*(k-1) + c
                nbj = nvirt*(j-1) + b
                kckbj = index(nck,nbj) 
                
                kbjki =nvirt*nrhft*nrhft*(i-1)
     *                + nvirt*nrhft*(k-1) + nbj
c           write(lupri,*)'nck,nbj,kckbj, xckbj'
c           write(lupri,*)nck,nbj,kckbj,XIAJB(kckbj) 
c           WRITE(LUPRI,*)'AI,I,A'
c           WRITE(LUPRI,*)AI,I,A
c           write(lupri,*)'bj,j,k,kbjki,tmat(kbjki)'
c           write(lupri,*)nbj,j,k,kbjki,tmat(kbjki)
c           write(lupri,*)'WORK(AI),XIAJB(kckbj)*tmat(kbjki'
c           write(lupri,*)WORK(AI),XIAJB(kckbj)*tmat(kbjki)
            WORK(AI) = WORK(AI) + XIAJB(kckbj)*tmat(kbjki)
c           write(lupri,*)'WORK(AI)',WORK(AI)
          end do
        end do
       end do
      end do

      write(lupri,*)' xi1eff wbdtest second cont a,c',a,c
      do i = 1,nt1am(1)
         write(lupri,*) i , work(I)
      end do
      end if
C
C
cc
      if (iopt.eq.3) then
      b = ic
      A = ib
      write(lupri,*)' a,b third contribution',a,b
      do i = 1,nrhft
       DO c = 1,nvirt
       ai = nvirt*(i-1)+a
        DO K = 1,NRHFT
          DO J= 1,NRHFT
                nck = nvirt*(k-1) + c
                nbj = nvirt*(j-1) + b
                kckbj = index(nck,nbj) 
                
                kckji =nvirt*nrhft*nrhft*(i-1)
     *                + nvirt*nrhft*(j-1) + nck
c           write(lupri,*)'nck,nbj,kckbj, xckbj'
c           write(lupri,*)nck,nbj,kckbj,XIAJB(kckbj) 
c           WRITE(LUPRI,*)'AI,I,A'
c           WRITE(LUPRI,*)AI,I,A
c           write(lupri,*)'bj,j,k,kbjki,tmat(kbjki)'
c           write(lupri,*)nbj,j,k,kbjki,tmat(kbjki)
c           write(lupri,*)'WORK(AI),XIAJB(kckbj)*tmat(kbjki'
c           write(lupri,*)WORK(AI),XIAJB(kckbj)*tmat(kbjki)
            WORK(AI) = WORK(AI) + XIAJB(kckbj)*tmat(kckji)
c           write(lupri,*)'WORK(AI)',WORK(AI)
          end do
        end do
       end do
      end do

      write(lupri,*)' xi1eff wbdtest third cont a,c',a,c
      do i = 1,nt1am(1)
         write(lupri,*) i , work(I)
      end do
      end if

      RETURN
      END
C  /* Deck CC3_CY2F */
      SUBROUTINE CC3_CY2F(XI2,ISYRES,WMAT,ISWMAT,TMAT,
     *                   FOCK0,ISYINT,INDSQ,LENSQ,WORK,LWORK,
     *                   ISYMIB,IB,ISYMID,ID)
C
C  P. Jorgensen and F. Pawlowski.         Spring 2002
C
      IMPLICIT NONE
C
C  Calculate Fock part to XI2 
C
C
C  General symmetry: ISWMAT is symmetry of WMAT and TMAT intermdiates.  
C                       (exclusive isymb,isymd) 
C                       ISYINT is symmetry of FOCK0 
C                       ISYRES is symmetry of XI2 and XI2EFF 
C 
C        Omega1(aibj)   
C
C        = sum_{ck} ( W^{abc}_{ijk}  -  W^{abc}_{kji} ) * F_{kc} 
C  (1):
C        = sum_{k} ( W^{BC}_{aikj} - W^{BC}_{akij} ) *  F_{kC} 
C  (2):
C        + sum_{k} ( W^{CA}_{bjik} - W^{CA}_{bjki} ) *  F_{kC}
C  (3):
C        + sum_{ck} ( W^{AB}_{ckji} - W^{AB}_{cijk} ) *  F_{kc}
C
C
C
      INTEGER ISYRES, ISWMAT, ISYINT , LENSQ, LWORK, KSTART
      INTEGER ISYMIB, IB, ISYMID, ID, INDEX
      INTEGER INDSQ(LENSQ,6)
      INTEGER KFMAT, KABIJ, KEND, LEND
      INTEGER ISYMB, ISYMC, ISYMK, ISYAIJ, KOFFT, KOFFF, NTOAIJ
      INTEGER ISYMA, ISYBJI, NTOBJI, ISYAB, ISYIJ 
      INTEGER KOFF1, KOFF2, ISYCK, NTOCK 
      REAL*8  XI2(*), WMAT(*), TMAT(*), FOCK0(*), WORK(*)
      REAL*8  ZERO, ONE, TWO

C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C 
      CALL QENTER('CC3_CY2F')

C       
C     work space allocations
C
      KSTART    = 1
      KEND    = KSTART  + NVIRT*NRHFT*NRHFT 
      LEND    = LWORK - KEND
C
      IF (LEND .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK 
         WRITE(LUPRI,*) 'Memory needed    : ',KEND
         CALL QUIT('Insufficient space in CC3_CY2F')
      END IF
C
C 
C----------------------------------
C  (1):
C     Omega2(aibj) = Omega2(aibj)  
C                   + sum_{k} ( W^{BC}_{aikj} - W^{BC}_{akij} ) *  F_{kC} 
C----------------------------------
C
C       
      B = IB
      C = ID
      ISYMB = ISYMIB
      ISYMC = ISYMID
C
C      T(aijk) =  W^{BC}_{aikj} - W^{BC}_{akij} )
C
      DO I = 1,LENSQ
C
         TMAT(I) =   WMAT(INDSQ(I,3))
     *             - WMAT(INDSQ(I,4))
C
      ENDDO
C
      ISYMK  = MULD2H(ISYINT,ISYMC) 
      ISYAIJ = MULD2H(ISYMK,ISWMAT) 
      KOFFT  = ISAIKJ(ISYAIJ,ISYMK) + 1
      KOFFF  = IFCVIR(ISYMK,ISYMC)  + NORB(ISYMK)*(C-1) + 1 
      NTOAIJ = MAX(NCKI(ISYAIJ),1)
C
      IF (NRHF(ISYMK).GT.0) THEN
         CALL DGEMV('N',NCKI(ISYAIJ),NRHF(ISYMK),-ONE,TMAT(KOFFT),
     *               NTOAIJ,FOCK0(KOFFF),1,ZERO,WORK,1)
C
        IF (NCKI(ISYAIJ).GT.0 ) THEN
           CALL CC3_RACC(XI2,WORK,ISYMB,B,ISYRES)
        END IF
      END IF
C
C----------------------------------
C  (2):
C     Omega2(aibj) = Omega2(aibj)  
C        + sum_{k} ( W^{CA}_{bjik} - W^{CA}_{bjki} ) *  F_{kC}
C----------------------------------
C
      C = IB
      A = ID
C
      ISYMC = ISYMIB
      ISYMA = ISYMID
C
C
C      T(bjik) =  W^{CA}_{bjik} - W^{AC}_{bjki} )  
C
      DO I = 1,LENSQ
C
         TMAT(I) =   WMAT(I) 
     *             - WMAT(INDSQ(I,3))
C
      ENDDO
C
C
      ISYMK  = MULD2H(ISYINT,ISYMC) 
      ISYBJI = MULD2H(ISYMK,ISWMAT) 
      KOFFT  = ISAIKJ(ISYBJI,ISYMK) + 1
      KOFFF  = IFCVIR(ISYMK,ISYMC)  + NORB(ISYMK)*(C-1) + 1 
      NTOBJI = MAX(NCKI(ISYBJI),1)
      IF (NRHF(ISYMK).GT.0) THEN
         CALL DGEMV('N',NCKI(ISYBJI),NRHF(ISYMK),-ONE,TMAT(KOFFT),
     *               NTOBJI,FOCK0(KOFFF),1,ZERO,WORK,1)
C
        IF (NCKI(ISYBJI).GT.0) THEN
           CALL CC3_RACC(XI2,WORK,ISYMA,A,ISYRES)
        END IF
      END IF
C
C----------------------------------
C  (3):
C     Omega2(aibj) = Omega2(aibj)  
C        + sum_{ck} ( W^{AB}_{ckji} - W^{AB}_{cijk} ) *  F_{kc}
C----------------------------------
C
      C = IB
C
      A = IB
      B = ID
C
      ISYMA = ISYMIB
      ISYMB = ISYMID
C
C
C      T(ckij) =  W^{AB}_{ckji} - W^{AB}_{cijk}
C
      DO I = 1,LENSQ
C
         TMAT(I) =   WMAT(INDSQ(I,3)) 
     *             - WMAT(INDSQ(I,2))
C
      ENDDO
C
      IF (NSYM .GT. 1) THEN
         CALL CC_GATHER(LENSQ,WORK,TMAT,INDSQ(1,6))
         CALL DCOPY(LENSQ,WORK,1,TMAT,1)
      ENDIF
C

C
C WORK SPACE ALLOCATION
C
      ISYAB = MULD2H(ISYMA,ISYMB)
      ISYIJ = MULD2H(ISYRES,ISYAB)
C
      KFMAT    = 1
      KABIJ   = KFMAT  + NT1AM(ISYINT)
      KEND    = KABIJ  + NMATIJ(ISYIJ) 
      LEND    = LWORK - KEND
C
      IF (LEND .LT. 0) THEN
         WRITE(LUPRI,*) 'Memory available : ',LWORK 
         WRITE(LUPRI,*) 'Memory needed    : ',KEND
         CALL QUIT('Insufficient space in CC3_CY2F')
      END IF
C
C
C   sort fock matrix 
C
C   FMAT(C,K) = FOCK0(K,C)
C
      DO ISYMC = 1,NSYM
         ISYMK = MULD2H(ISYMC,ISYINT)
         DO K = 1,NRHF(ISYMK)
            DO C = 1,NVIR(ISYMC)
               KOFF1 = IFCVIR(ISYMK,ISYMC) + NORB(ISYMK)*(C - 1) + K
               KOFF2 = KFMAT +  IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K - 1) + C -1
C
                  WORK(KOFF2) = FOCK0(KOFF1)
C
            END DO
         END DO
      END DO
C
      ISYCK = ISYINT
      KOFFT = ISAIKL(ISYCK,ISYIJ) + 1
      NTOCK = MAX(NT1AM(ISYCK),1)
      IF (NMATIJ(ISYIJ).GT.0) THEN
C
         CALL DGEMV('T',NT1AM(ISYCK),NMATIJ(ISYIJ),-ONE,TMAT(KOFFT),
     *               NTOCK,WORK(KFMAT),1,ZERO,WORK(KABIJ),1)
         CALL CC3_RABIJ(XI2,WORK(KABIJ),ISYMA,A,ISYMB,B,ISYRES)
C
      END IF
C
      CALL QEXIT('CC3_CY2F') 
C
C
      RETURN
      END
C
C  /* Deck CC3_RABIJ */
      SUBROUTINE CC3_RABIJ(XI2,ABIJ,ISYMA,A,ISYMB,B,ISYRES)
C
C distribute matrix elements in ABIJ(I,J) to result vector XI2(AIBJ)
C
      IMPLICIT NONE 
C
      INTEGER ISYMA, ISYMB, ISYRES, INDEX
      INTEGER ISYAB, ISYIJ, ISYMJ, ISYMI, ISYMBJ, ISYMAI
      INTEGER NBJ, NIJ, NAI, NAIBJ 
      REAL*8  XI2(*), ABIJ(*)
      REAL*8  TWO
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_RABIJ')
C
      ISYAB = MULD2H(ISYMA,ISYMB)
      ISYIJ = MULD2H(ISYAB,ISYRES)
C
      DO 300 ISYMJ = 1,NSYM
         ISYMI  = MULD2H(ISYIJ,ISYMJ)
         ISYMBJ = MULD2H(ISYMB,ISYMJ)
         ISYMAI = MULD2H(ISYMA,ISYMI)
C
         DO 310 J = 1,NRHF(ISYMJ)
            NBJ = IT1AM(ISYMB,ISYMJ) + NVIR(ISYMB)*(J - 1) + B
            IF (ISYMAI .EQ. ISYMBJ) THEN
C
               DO 320 I = 1,NRHF(ISYMI)
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NAI = IT1AM(ISYMA,ISYMI)  + NVIR(ISYMA)*(I - 1) + A
                  IF (NAI .EQ. NBJ) ABIJ(NIJ) = TWO*ABIJ(NIJ)
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(NAI,NBJ)
                  XI2(NAIBJ) = XI2(NAIBJ) + ABIJ(NIJ)
  320          CONTINUE
C
            ELSE IF (ISYMAI .LT. ISYMBJ) THEN
C
               DO 330 I = 1,NRHF(ISYMI)
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NAI = IT1AM(ISYMA,ISYMI)  + NVIR(ISYMA)*(I - 1) + A
                  NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                  + NT1AM(ISYMAI)*(NBJ-1) + NAI
                  XI2(NAIBJ) = XI2(NAIBJ) + ABIJ(NIJ)
  330          CONTINUE
C
            ELSE IF (ISYMBJ .LT. ISYMAI) THEN
C
               DO 340 I = 1,NRHF(ISYMI)
                  NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J - 1) + I
                  NAI = IT1AM(ISYMA,ISYMI)  + NVIR(ISYMA)*(I - 1) + A
                  NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                  + NT1AM(ISYMBJ)*(NAI-1) + NBJ
                  XI2(NAIBJ) = XI2(NAIBJ) + ABIJ(NIJ)
  340          CONTINUE
C
            ENDIF
  310    CONTINUE
C
  300 CONTINUE
C
      CALL QEXIT('CC3_RABIJ') 
C
      RETURN
      END
C  /* Deck cc3_cy2o */
      SUBROUTINE CC3_CY2O(XI2,ISYRES,WMAT,ISWMAT,TMAT,TROCC,
     *                      TROCC1,ISYINT,WORK,LWORK,INDSQ,LENSQ,
     *                      ISYMIB,IB,ISYMID,ID)
C
C
C     General symmetry: ISWMAT is symmetry of WMAT and TMAT intermediates.
C                       (exclusive isymib,isymid)
C                       ISYINT is symmetry of integrals in TROCC and TROCC1.
C                       ISYRES = ISWMAT*ISYINT*ISYMIB*ISYMID
C
C   XI2(aibj) =  XI2(aibj) 
C
C     + sum_{ckl} (2W^{abc}_{ikl}-W^{abc}_{lki}-W^{abc}_{ilk}) * (lc|kj)
C (1):
C     = sum_{kl} (2W^{BC}_{ailk}-W^{BC}_{alik}-W^{BC}_{aikl}) * (lc|kj) 
C (2):
C     + sum_{kl} (2W^{CA}_{bkil}-W^{CA}_{bkli}-W^{CA}_{blik}) * (lc|kj) 
C (3):
C    + sum_{ckl} (2W^{AB}_{clki}-W^{AB}_{cikl}-W^{AB}_{ckli}) * (lc|kj) 
C
      IMPLICIT NONE
C
      INTEGER ISYRES, ISWMAT, ISYINT, LWORK, LENSQ 
      INTEGER ISYMIB, IB, ISYMID, ID 
      INTEGER INDSQ(LENSQ,6)
      INTEGER ISYMC, ISYMB, ISYBC, JSAIKL, LENGTH, ISYMJ, ISYBJ
      INTEGER ISYAI, ISYKL, NTOTAI 
      INTEGER NTOTKL, KOFF1, KOFF2, KOFF3 
      INTEGER ISYMA, ISYAB, JSCKLI, JSBIKL, ISYMI
      INTEGER ISYCKL, NTOCKL, NRHFI 
      INTEGER KRMAT1, KRMAT2, KEND, LEND, INDEX
      INTEGER ISYCK, ISYCJ, ISYAC, ISYBI, ISYKLJ, NTOTBI
      INTEGER ISWB, ISWD, NCKIMX
      INTEGER ISYAIJ, ISYBJI, ISYIJ 
      INTEGER KTMATX
C
      REAL*8  XI2(*), WMAT(*), TMAT(*), TROCC(*), TROCC1(*)
      REAL*8  WORK(*) 
      REAL*8  TWO, ONE, ZERO, XTMAT, XRMAT, DDOT
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_CY2O')
C
      
      ISWB   = MULD2H(ISWMAT,ISYMIB)
      ISWD   = MULD2H(ISWMAT,ISYMID)
      ISYAIJ = MULD2H(ISYMIB,ISYRES)
      ISYBJI = MULD2H(ISYMID,ISYRES)
      NCKIMX = MAX(NCKI(ISWB),NCKI(ISWD),NCKI(ISYAIJ),NCKI(ISYBJI)) 
C
      KRMAT1 = 1
      KRMAT2 = KRMAT1 + NCKIMX
      KTMATX = KRMAT2 + NCKIMX
      KEND   = KTMATX + NCKIJ(ISWMAT)
      LEND   = LWORK  - KEND
      IF (LWORK .LT. KEND) THEN
         CALL QUIT('Insufficient space in CY2O')
      ENDIF
C
C-------------------------
C     First occupied term.
C-------------------------
C (1):
C             XI2(aibj) =  XI2(aibj) 
C
C     + sum_{kl} (2W^{BC}_{ailk}-W^{BC}_{alik}-W^{BC}_{aikl}) * (lc|kj)

C
      B = IB
      C = ID
C
      ISYMB = ISYMIB
      ISYMC = ISYMID
C
      ISYAIJ = MULD2H(ISYMB,ISYRES)
      IF (NCKI(ISYAIJ).GT.0 ) THEN 
        CALL DZERO(WORK(KRMAT1),NCKI(ISYAIJ))
C
        ISYBC = MULD2H(ISYMB,ISYMC)
        JSAIKL = ISWMAT 
C
        LENGTH = NCKIJ(JSAIKL)
C
C---------------------------------------
C
C     TMAT(aikl) = (2W^{BC}_{ailk}-W^{BC}_{alik}-W^{BC}_{aikl}) 
C
C---------------------------------------
C
        DO I = 1,LENGTH
C
           TMAT(I) =   TWO*WMAT(INDSQ(I,3))
     *               -     WMAT(INDSQ(I,4))
     *               -     WMAT(I) 
C
        ENDDO
C
C----------------------------------
C     Symmetry sorting if symmetry.
C----------------------------------
C
      
        IF (NSYM .GT. 1) THEN
           CALL CC_GATHER(LENGTH,WORK(KTMATX),TMAT,INDSQ(1,6))
           CALL DCOPY(LENGTH,WORK(KTMATX),1,TMAT,1)
        ENDIF
C
        IF (IPRINT .GT. 55) THEN
           XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
           WRITE(LUPRI,*) 'In CC3_CY2O: 1. Norm of TMAT = ',XTMAT
        ENDIF
C
C-----------------------
C     First contraction.
C-----------------------
C
C         TROCC(kljC) = (lc|kj)
C
C   XI2(aibj) =  XI2(aibj) 
C             + sum(kl) TMAT(aikl)* TROCC(kljC) 
C
        DO 200 ISYMJ = 1,NSYM
C
           ISYBJ = MULD2H(ISYMB,ISYMJ)
           ISYAI = MULD2H(ISYBJ,ISYRES)
           ISYKL = MULD2H(JSAIKL,ISYAI)
           ISYKLJ = MULD2H(ISYKL,ISYMJ)
C
           NTOTAI = MAX(NT1AM(ISYAI),1)
           NTOTKL = MAX(NMATIJ(ISYKL),1)
C
           KOFF1  = ISAIKL(ISYAI,ISYKL) + 1
           KOFF2  = ISJIKA(ISYKLJ,ISYMC)
     *          + NMAJIK(ISYKLJ)*(C - 1)
     *          + ISJIK(ISYKL,ISYMJ) + 1
           KOFF3  = KRMAT1 + ISAIK(ISYAI,ISYMJ) 
C
           CALL DGEMM('N','N',NT1AM(ISYAI),NRHF(ISYMJ),NMATIJ(ISYKL),
     *                ONE,TMAT(KOFF1),NTOTAI,TROCC(KOFF2),NTOTKL,
     *                ONE,WORK(KOFF3),NTOTAI) 
C
  200   CONTINUE
C
        IF (IPRINT .GT. 55) THEN
           XRMAT = DDOT(NCKI(ISYAIJ),WORK(KRMAT1),1,WORK(KRMAT1),1)
           WRITE(LUPRI,*) '1. In CC3_CY2O: Norm of RMAT1 =  ',XRMAT
        ENDIF
C
        CALL CC3_RACC(XI2,WORK(KRMAT1),ISYMB,B,ISYRES)
C
      END IF
C
C--------------------------
C     Second occupied term.
C--------------------------
C
C (2):
C
C   XI2(aibj) =  XI2(aibj)
C
C     + sum_{kl} (2W^{CA}_{bkil}-W^{CA}_{bkli}-W^{CA}_{blik}) * (lc|kj)
C 
      A = ID
      C = IB
C
      ISYMA = ISYMID
      ISYMC = ISYMIB
C
      ISYBJI = MULD2H(ISYMA,ISYRES)
      IF (NCKI(ISYBJI).GT.0) THEN
        CALL DZERO(WORK(KRMAT1),NCKI(ISYBJI))
C
        ISYAC = MULD2H(ISYMA,ISYMC)
        JSBIKL = ISWMAT 
C
C---------------------------------------
C
C     TMAT(bikl) = (2W^{CA}_{bkil}-W^{CA}_{bkli}-W^{CA}_{blik}) 
C
C---------------------------------------
C
        LENGTH = NCKIJ(JSBIKL)
C
        DO I = 1,LENGTH
C
           TMAT(I) =   TWO*WMAT(INDSQ(I,1)) 
     *               -     WMAT(INDSQ(I,2))
     *               -     WMAT(INDSQ(I,4))
C
        ENDDO
C
C----------------------------------
C     Symmetry sorting if symmetry.
C----------------------------------
C
        IF (NSYM .GT. 1) THEN
           CALL CC_GATHER(LENGTH,WORK(KTMATX),TMAT,INDSQ(1,6))
           CALL DCOPY(LENGTH,WORK(KTMATX),1,TMAT,1)
        ENDIF
C
        IF (IPRINT .GT. 55) THEN
           XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
           WRITE(LUPRI,*) 'In CC3_CY2O: 2. Norm of TMAT = ',XTMAT
        ENDIF
C
C------------------------
C     Second contraction.
C------------------------
C
C         TROCC(kljC) = (lc|kj)
C
C      M^A_(bij) =   sum(kl) TMAT(bikl)* TTROCC(kljC)
C

        DO 400 ISYMJ = 1,NSYM
C
           ISYCJ = MULD2H(ISYMC,ISYMJ)
           ISYKL = MULD2H(ISYINT,ISYCJ)
           ISYBI = MULD2H(ISYKL,ISWMAT)
           ISYKLJ = MULD2H(ISYKL,ISYMJ)
C
           NTOTBI = MAX(NT1AM(ISYBI),1)
           NTOTKL = MAX(NMATIJ(ISYKL),1)
C
           KOFF1  = ISAIKL(ISYBI,ISYKL) + 1
           KOFF2  = ISJIKA(ISYKLJ,ISYMC)
     *            + NMAJIK(ISYKLJ)*(C - 1)
     *            + ISJIK(ISYKL,ISYMJ) + 1
           KOFF3  = KRMAT1 + ISAIK(ISYBI,ISYMJ) 
C

           CALL DGEMM('N','N',NT1AM(ISYBI),NRHF(ISYMJ),NMATIJ(ISYKL),
     *                ONE,TMAT(KOFF1),NTOTBI,TROCC(KOFF2),NTOTKL,
     *                ONE,WORK(KOFF3),NTOTBI)
C
  400   CONTINUE
C

        IF (IPRINT .GT. 55) THEN
           XRMAT = DDOT(NCKI(ISYBJI),WORK(KRMAT1),1,WORK(KRMAT1),1)
           WRITE(LUPRI,*) '2.In CC3_CY2O: Norm of RMAT1 =  ',XRMAT 
        ENDIF
C 
C reorder result vector M^A_(bij) as M^A_(bji)
C
        CALL CC3_MAJI(WORK(KRMAT1),WORK(KRMAT2),ISYMA,A,ISYRES) 
C
        IF (IPRINT .GT. 55) THEN
           XRMAT = DDOT(NCKI(ISYBJI),WORK(KRMAT2),1,WORK(KRMAT2),1)
           WRITE(LUPRI,*) 
     *         'In CC3_CY2O: Reorder :Norm of RMAT2 =  ',XRMAT 
        ENDIF
C
C add vector to XI2
C
        CALL CC3_RACC(XI2,WORK(KRMAT2),ISYMA,A,ISYRES) 
C
      END IF
C
C-------------------------
C     Third occupied term.
C-------------------------
C (3):
C             XI2(aibj) =  XI2(aibj) 
C
C    + sum_{ckl} (2W^{AB}_{clki}-W^{AB}_{cikl}-W^{AB}_{ckli}) * (lc|kj) 
C
      B = ID
      A = IB
C
      ISYMB = ISYMID
      ISYMA = ISYMIB
C
      ISYAB = MULD2H(ISYMA,ISYMB)
      ISYIJ = MULD2H(ISYAB,ISYRES)
      IF (NMATIJ(ISYIJ).GT.0) THEN
        CALL DZERO(WORK,NMATIJ(ISYIJ))
C
        JSCKLI = ISWMAT 
C
        LENGTH = NCKIJ(JSCKLI)
C
C----------------------------------
C
C   TMAT(ckli) = 2W^{AB}_{clki}-W^{AB}_{cikl}-W^{AB}_{ckli}
C
C----------------------------------
C
        DO I = 1,LENGTH
C
           TMAT(I) =  TWO*WMAT(INDSQ(I,1))    
     *               -    WMAT(INDSQ(I,4)) 
     *               -    WMAT(I)
C
        ENDDO
C
C
        IF (IPRINT .GT. 55) THEN
           XTMAT = DDOT(NCKIJ(JSAIKL),TMAT,1,TMAT,1)
           WRITE(LUPRI,*) 'In CC3_CY2O: 3. Norm of TMAT = ',XTMAT
        ENDIF
C
C-----------------------
C     Third contraction.
C-----------------------
C
C
C         TROCC1(cklj) = (lc|kj)
C
C   XI2(aibj) =  XI2(aibj)
C             + sum(kl) TMAT(ckli)* TROCC1(cklj) 
C
C
        DO 600 ISYMJ = 1,NSYM
C
           ISYBJ = MULD2H(ISYMB,ISYMJ)
           ISYAI = MULD2H(ISYBJ,ISYRES)
           ISYMI  = MULD2H(ISYAI,ISYMA)
           ISYCKL = MULD2H(ISYMI,JSCKLI)
C
           IF (LWORK .LT. NRHF(ISYMI)*NRHF(ISYMJ)) THEN
              CALL QUIT('Insufficient memory in CCSDT_CY2O')
           END IF
C
           NTOCKL = MAX(NCKI(ISYCKL),1)
           NRHFI  = MAX(NRHF(ISYMI),1)
C
           KOFF1  = ISAIKJ(ISYCKL,ISYMI) + 1
           KOFF2  = ISAIKJ(ISYCKL,ISYMJ) + 1
           KOFF3  = IMATIJ(ISYMI,ISYMJ) + 1
C

           CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NCKI(ISYCKL),
     *                ONE,TMAT(KOFF1),NTOCKL,TROCC1(KOFF2),NTOCKL,
     *                ONE,WORK(KOFF3),NRHFI)
C
C

  600   CONTINUE
C
        IF (IPRINT .GT. 55) THEN
           XRMAT = DDOT(NMATIJ(ISYIJ),WORK,1,WORK,1)
           WRITE(LUPRI,*) '3.In CC3_CY2O: Norm of RABIJ =  ',XRMAT 
        ENDIF
C
        CALL CC3_RABIJ(XI2,WORK,ISYMA,A,ISYMB,B,ISYRES)
C
      END IF
C
      CALL QEXIT('CC3_CY2O')
C
      RETURN
      END
C  /* Deck CC3_MAJI */
      SUBROUTINE CC3_MAJI(RMAT1,RMAT2,ISYMB,B,ISYRES) 
C
C   order matrix  M^B_(aji) as  M^B_(aij)
C
      IMPLICIT NONE
C
      INTEGER ISYMB,ISYRES
      INTEGER ISYAJI, ISYMI, ISYAJ, ISYMJ, ISYMA, ISYAI
      INTEGER KOFF1, KOFF2, KMAT, KEND0, LWRK0
      REAL*8  RMAT1(*), RMAT2(*)
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>

C
      CALL QENTER('CC3_MAJI')
C
      ISYAJI   = MULD2H(ISYRES,ISYMB)
      DO ISYMI = 1,NSYM
         ISYAJ = MULD2H(ISYMI,ISYAJI)
         DO I = 1,NRHF(ISYMI)
            DO ISYMJ = 1,NSYM
               ISYMA = MULD2H(ISYAJ,ISYMJ)
               ISYAI = MULD2H(ISYMA,ISYMI)
               DO J = 1,NRHF(ISYMJ)
                  DO A = 1,NVIR(ISYMA)
                     KOFF1 = ISAIK(ISYAJ,ISYMI)
     *                     + NT1AM(ISYAJ)*(I-1)
     *                     + IT1AM(ISYMA,ISYMJ)
     *                     + NVIR(ISYMA)*(J-1) + A 
                     KOFF2 = ISAIK(ISYAI,ISYMJ)
     *                     + NT1AM(ISYAI)*(J-1)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1) + A 

                     RMAT2(KOFF2) = RMAT1(KOFF1)
                  END DO
               END DO
            END DO 
         END DO 
      END DO 
C
      CALL QEXIT('CC3_MAJI')
C

      RETURN
      END
C  /* Deck cc3_convir */
      SUBROUTINE CC3_CY2V(XI2,ISYRES,RBJIA,WMAT,ISWMAT,TMAT,TRVIR,
     *                      TRVIR1,ISYINT,WORK,LWORK,INDSQ,LENSQ,
     *                      ISYMIB,IB,ISYMID,ID)
C
C    
C General symmetry: ISWMAT is symmetry of WMAT and TMAT intermediates.
C                   (exclusive isymib,isymid)
C                   ISYINT is symmetry of integrals in TRVIR and TROVIR1.
C                   ISYRES = ISWMAT*ISYINT*ISYMIB*ISYMID
C                   RBJIA intermediate result vector
C
C   XI2(aibj) =  XI2(aibj)
C
C     + sum_{cdk} (2W^{bcd}_{jik}-W^{bcd}_{kij}-W^{bcd}_{jki}) * (ac|kd)
C (1):
C     = sum_{ck} (2W^{DB}_{cijk}-W^{DB}_{cikj}-W^{DB}_{ckji}) * (ac|kD)
C (2):
C     + sum_{dk} (2W^{BC}_{dkij}-W^{BC}_{djik}-W^{BC}_{dikj}) * (aC|kd) 
C (3):
C    + sum_{k} (2W^{CD}_{bjki}-W^{CD}_{bkji}-W^{CD}_{bjik}) * (aC|kD)
C
      IMPLICIT NONE
C
      INTEGER ISYRES, ISWMAT, ISYINT , LENSQ, LWORK
      INTEGER INDSQ(LENSQ,6)
      INTEGER ISYMIB, IB, ISYMID, ID
      INTEGER INDEX ,LENGTH, ISCKIJ, ISYMB, ISYAIJ, KRMAT, KEND, LEND
      INTEGER ISYMJ ,ISYBJ, ISYAI, ISYCKI, ISYCK, ISYMA, NTOTCK, NVIRA
      INTEGER KOFF1 , KOFF2, KOFF3, ISYMD, ISDKIJ, ISYDKI
      INTEGER ISYDK, NTOTDK, ISBJIK, ISYCKA
      INTEGER ISYKA, KINT, ISYBJI, KOFFT, KOFFM, KOFFR
      INTEGER NTOK, NTOBJI, ISYMK, ISYMC, ISYMI
      REAL*8  XI2(*), RBJIA(*), WMAT(*), TMAT(*), TRVIR(*), TRVIR1(*) 
      REAL*8  WORK(*)
      REAL*8  ZERO, ONE, TWO

C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
C
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC3_CY2V')
      LENGTH = NCKIJ(ISWMAT)
C
C------------------------
C (1):
C     = sum_{ck} (2W^{DB}_{cijk}-W^{DB}_{cikj}-W^{DB}_{ckji}) * (ac|kD)
C
C------------------------
C
C   TMAT(ckij) = 2W^{DB}_{cijk}-W^{DB}_{cikj}-W^{DB}_{ckji}
C
C   Use here the symmetry between BJ, DK !!
C
      ISCKIJ = ISWMAT 
C
      ISYMB = ISYMIB
      ISYMD = ISYMID
      B = IB
      D = ID
C
      ISYAIJ = MULD2H(ISYMB,ISYRES)
      KRMAT  = 1
      KEND   = KRMAT + NCKI(ISYAIJ)
      LEND   = LWORK - KEND
      IF (LEND .LT. 0)THEN
         CALL QUIT('1. Insufficient core in CC3_CY2V')
      ENDIF
      CALL DZERO(WORK(KRMAT), NCKI(ISYAIJ))
      DO I = 1,LENGTH
         TMAT(I) =  TWO*WMAT(INDSQ(I,1))
     *           -      WMAT(INDSQ(I,2))
     *           -      WMAT(I)
      ENDDO
C
C---------------------------
C     (ac|kD) = TRVIR^D_{ck,a}
C
C     RMAT^B_(aij) =  TRVIR^D_{ck,a} )^T *  TMAT(ckij)
C---------------------------
C
      DO 200 ISYMJ = 1,NSYM
C
         ISYBJ = MULD2H(ISYMB,ISYMJ)
         ISYAI = MULD2H(ISYBJ,ISYRES)
         ISYCKI = MULD2H(ISCKIJ,ISYMJ)
C
         DO 210 J = 1,NRHF(ISYMJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYCK = MULD2H(ISYCKI,ISYMI)
               ISYMA  = MULD2H(ISYAI,ISYMI)
C 
               NTOTCK = MAX(NT1AM(ISYCK),1)
               NVIRA  = MAX(NVIR(ISYMA),1)
C
               KOFF1 = ICKATR(ISYCK,ISYMA) + 1
               KOFF2 = ISAIKJ(ISYCKI,ISYMJ)
     *               + NCKI(ISYCKI)*(J - 1)
     *               + ISAIK(ISYCK,ISYMI)  + 1
               KOFF3 = KRMAT
     *               + ISAIK(ISYAI,ISYMJ)
     *               + NT1AM(ISYAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI)
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NT1AM(ISYCK),
     *                   -ONE,TRVIR(KOFF1),NTOTCK,TMAT(KOFF2),NTOTCK,
     *                    ONE,WORK(KOFF3),NVIRA) 
C
  220       CONTINUE
C
  210    CONTINUE
C
  200 CONTINUE
C
      ISYAIJ = MULD2H(ISYMB,ISYRES)
      IF ( NCKI(ISYAIJ).GT.0 ) THEN
        CALL CC3_RACC(XI2,WORK(KRMAT),ISYMB,B,ISYRES)
      END IF
C
C
C (2):
C
C   XI2(aibj) =  XI2(aibj)
C
C     + sum_{dk} (2W^{BC}_{dkij}-W^{BC}_{djik}-W^{BC}_{dikj}) * (aC|kd)
C
C  TMAT(dkij) = 2W^{BC}_{dkij}-W^{BC}_{djik}-W^{BC}_{dikj}
C
      ISYMC = ISYMID
      ISYMB = ISYMIB
      B  =  IB
      C  =  ID 
C
      ISDKIJ = ISWMAT
      ISYAIJ = MULD2H(ISYMB,ISYRES)
      KRMAT  = 1
      KEND   = KRMAT + NCKI(ISYAIJ)
      LEND   = LWORK - KEND
      IF (LEND .LT. 0)THEN
         CALL QUIT('2. Insufficient core in CC3_CY2V')
      ENDIF
      CALL DZERO(WORK(KRMAT), NCKI(ISYAIJ))
C 
      DO I = 1,LENGTH
         TMAT(I) =  TWO*WMAT(I) 
     *             -    WMAT(INDSQ(I,5))
     *             -    WMAT(INDSQ(I,1))
      ENDDO
C
C     (ad|kC) = TRVIR1^C(dka)
C
C  RMAT^B_(aij) = ( TRVIR1^C(dka) )^T * TMAT(dkij)
C
      DO 300 ISYMJ = 1,NSYM
C
         ISYBJ = MULD2H(ISYMB,ISYMJ)
         ISYAI = MULD2H(ISYBJ,ISYRES)
         ISYDKI = MULD2H(ISDKIJ,ISYMJ)
C
         DO 310 J = 1,NRHF(ISYMJ)
C
            DO 320 ISYMI = 1,NSYM
C
               ISYDK = MULD2H(ISYDKI,ISYMI)
               ISYMA  = MULD2H(ISYAI,ISYMI)
C
               NTOTDK = MAX(NT1AM(ISYDK),1)
               NVIRA  = MAX(NVIR(ISYMA),1)
C
               KOFF1 = ICKATR(ISYDK,ISYMA) + 1
               KOFF2 = ISAIKJ(ISYDKI,ISYMJ)
     *               + NCKI(ISYDKI)*(J - 1)
     *               + ISAIK(ISYDK,ISYMI)  + 1
               KOFF3 = KRMAT + ISAIK(ISYAI,ISYMJ)
     *               + NT1AM(ISYAI)*(J - 1)
     *               + IT1AM(ISYMA,ISYMI) 
C
               CALL DGEMM('T','N',NVIR(ISYMA),NRHF(ISYMI),NT1AM(ISYDK),
     *                   -ONE,TRVIR1(KOFF1),NTOTDK,TMAT(KOFF2),NTOTDK,
     *                    ONE,WORK(KOFF3),NVIRA)
C
  320       CONTINUE
  310    CONTINUE
  300 CONTINUE
C
      ISYAIJ = MULD2H(ISYMB,ISYRES)
      IF ( NCKI(ISYAIJ).GT.0 ) THEN
        CALL CC3_RACC(XI2,WORK(KRMAT),ISYMB,B,ISYRES)
      END IF
C
C-----------------------------------
C (3):
C    + sum_{k} (2W^{CD}_{bjki}-W^{CD}_{bkji}-W^{CD}_{bjik}) * (ac|kD)
C
C-----------------------------------
C
      ISYMC = ISYMIB
      ISYMD = ISYMID
      C = IB
      D = ID
C
C     TMAT(bjik) = 2W^{CD}_{bjki}-W^{CD}_{bkji}-W^{CD}_{bjik}
C
      ISBJIK = ISWMAT
C
      DO I = 1,LENGTH
         TMAT(I) =  TWO*WMAT(INDSQ(I,3))
     *           -      WMAT(INDSQ(I,4))
     *           -      WMAT(I)
      ENDDO
C
C     (ac|kD) = TRVIR^D_{Cka} = M^DC_{ka}
C
c
      ISYCKA = MULD2H(ISYINT,ISYMD)
      ISYKA  = MULD2H(ISYCKA,ISYMC)
C
      DO  ISYMA = 1,NSYM
         ISYCK  = MULD2H(ISYCKA,ISYMA)
         ISYMK  = MULD2H(ISYMC,ISYCK)
         DO A = 1,NVIR(ISYMA)
            DO K = 1,NRHF(ISYMK) 
               KOFFM = IT1AMT(ISYMK,ISYMA)
     *               + NRHF(ISYMK)*(A-1) + K 
               KINT  = ICKATR(ISYCK,ISYMA)
     *               + NT1AM(ISYCK)*(A-1) 
     *               + IT1AM(ISYMC,ISYMK)
     *               + NVIR(ISYMC)*(K-1) + C 
               WORK(KOFFM) = TRVIR(KINT)
            END DO
         END DO
      END DO
C
C  RBJIA(bjia) = RBJIA(bjia) + TMAT(bjik) *  M^DC_{ka}
C
      DO ISYMA = 1,NSYM
         ISYMK = MULD2H(ISYMA,ISYKA)
         ISYBJI = MULD2H(ISBJIK,ISYMK)
         KOFFT  = ISAIKJ(ISYBJI,ISYMK) + 1 
         KOFFM  = IT1AMT(ISYMK,ISYMA)  + 1
         KOFFR  = IT2SP(ISYBJI,ISYMA)  + 1
         NTOK = MAX(NRHF(ISYMK),1)
         NTOBJI = MAX(NCKI(ISYBJI),1)
C
         CALL DGEMM('N','N',NCKI(ISYBJI),NVIR(ISYMA),
     *               NRHF(ISYMK),-ONE,TMAT(KOFFT),NTOBJI, 
     *               WORK(KOFFM),NTOK,ONE,RBJIA(KOFFR),NTOBJI)
C
      END DO
C
      CALL QEXIT('CC3_CY2V')
      RETURN
      END
C  /* Deck CC3_RBJIA */
      SUBROUTINE CC3_RBJIA(XI2,ISYRES,RBJIA) 
C
C   distribute RBJIA in result vector XI2
C    XI2(ai,bj) = XI2(ai,bj) + RBJIA(bjai) + RBJIA(aibj)
C  
C    XI2 triangular packed
C    RBJIA squared packed
C
      IMPLICIT NONE
C
      INTEGER ISYMA,ISYRES,ISYBJI,KOFFR
C
      REAL*8  XI2(*), RBJIA(*)
C
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>

C
      CALL QENTER('CC3_RBJIA')
C
      DO  ISYMA = 1,NSYM
         ISYBJI = MULD2H(ISYRES,ISYMA)
         DO A = 1,NVIR(ISYMA)
            KOFFR = IT2SP(ISYBJI,ISYMA)
     *            + NCKI(ISYBJI)*(A-1) + 1 
            CALL CC3_RACC(XI2,RBJIA(KOFFR),ISYMA,A,ISYRES)
         END DO
      END DO
C
      CALL QEXIT('CC3_RBJIA')
C
      RETURN
      END
C
C    /* Deck get_filename */
      SUBROUTINE GET_FILENAME(FNCKJDR,FNDKBCR,MCAU)
*
************************************************************
*     Based on MCAU value get the name of the files for C1-transformed
*     integrals in CC3 Cauchy vector calculations
*
*     Currently MCAU may not be larger than 999
*
*     Filip Pawlowski, 11-May-2004, Aarhus
************************************************************
*
      IMPLICIT NONE
C
#include <priunit.h>
C
      LOGICAL LCDBG
      PARAMETER (LCDBG = .FALSE.)
C
      INTEGER MCAU
      INTEGER UNITS,TENTH,HUNDREDS
      CHARACTER*(*) FNCKJDR,FNDKBCR
      CHARACTER*3  NFIL
C
      CALL QENTER('GETFN')
C
      !Convert MCAU to character string NFIL
      IF ((MCAU.GE.0) .AND. (MCAU.LE.9)) THEN
C
         NFIL = '_'//'_'//CHAR(MCAU+48)
C
      ELSE IF ((MCAU.GE.10) .AND. (MCAU.LE.99)) THEN
C
         TENTH = DINT(AINT(DFLOAT(MCAU)/10.0D0))
         UNITS = MCAU - TENTH*10
         NFIL = '_'//CHAR(TENTH+48)//CHAR(UNITS+48)
C
      ELSE IF ((MCAU.GE.100) .AND. (MCAU.LE.999)) THEN
C
         HUNDREDS = DINT(AINT(DFLOAT(MCAU)/100.0D0))
         TENTH    = DINT(AINT(DFLOAT(MCAU-HUNDREDS*100)/10.0D0))
         UNITS    = MCAU - HUNDREDS*100 - TENTH*10
         NFIL = CHAR(HUNDREDS+48)//CHAR(TENTH+48)//CHAR(UNITS+48)
C
      ELSE
         WRITE(LUPRI,*)'MCAU = ', MCAU
         CALL QUIT('Case not implemented in GET_FILENAME')
      END IF
C
      !Generate file name for C1-transformed integrals
      FNCKJDR = 'CC3_CAUINT_O_'//NFIL
      FNDKBCR = 'CC3_CAUINT_V_'//NFIL
C
      IF (LCDBG) THEN
         WRITE(LUPRI,*)'MCAU = ',MCAU
         WRITE(LUPRI,*)'FNCKJDR = ',FNCKJDR
         WRITE(LUPRI,*)'FNDKBCR = ',FNDKBCR
         WRITE(LUPRI,*)' '
      END IF
C
      !End
      CALL QEXIT('GETFN')
C
      RETURN
      END
C  /* Deck delete_file */
      SUBROUTINE DELETE_FILES(FILENAME)
*
*********************************************************************
*     Delete file(s) with name FILENAME
*
*     NOTE: the file name FILENAME can contain wildcard
*     characters, for example CALL DELETE_FILES('CC3_CAUINT_V*')
*
*     Filip Pawlowski, Aarhus, 16-Apr-2004
*********************************************************************
*
      IMPLICIT NONE
C
      INTEGER ILEN,INFO
      CHARACTER*(*)  FILENAME
      CHARACTER*(80) COMMAND
C
      !Determine length of the filename 
      !(cannot be larger than 74 due to the way COMMAND is declared)
      ILEN = MIN(LEN(FILENAME),74)
C
      !Generate the command to delete the file
      WRITE(COMMAND,'(2A)') 'rm -f ', FILENAME(1:ILEN)
C
      !Execute the command
#if defined (SYS_CRAY)
      CALL INFO = ISHELL(COMMAND)
#else
      CALL SYSTEM(COMMAND)
#endif

      RETURN
      END


