C    /* Deck cc3_xi */ 
      SUBROUTINE CCSDT_FBMAT(IBTRAN,IDOTS,DOTPROD,NBTRAN,MXVEC,
     *                       LISTL,LISTB,LISTC,FNCKJD,FNDKBC,FNDELD,
     *                       FNTOC,FN3VI,FN3VI2,WORK,LWORK)
C
C     Written by K. Hald, Summer 2002.
C
C     Calculate <L2|[[H,R_1],R_3]|HF>
C
C     Initially construct 
C
C     t3(ai,Bj,Ck) = S^BC(aikj) + U^BC(aikj) + S^CB(aijk) + U^CB(aijk)
C
C     W^BC(aikj) = sum_d X(ad)*t3(di,Bj,Ck) - sum_l X(li)*t3(al,Bj,Ck)
C                 - sum_ld X(ld)*t2(ai,Cl)*t2(Bj,dk) 
C                 - sum_ld X(ld)*t2(ai,Bl)*t2(Ck,dj) 
C                 + ( -<mu3 | [[H,T1^Y],T2] | HF> -<mu3| [H,T2^Y] |HF> )
C
      IMPLICIT NONE
C
      INTEGER NBTRAN, MXVEC, LWORK
      INTEGER IBTRAN(3,NBTRAN), IDOTS(MXVEC,NBTRAN)
      INTEGER ISINT1, ISINT2, ISYMT2, ISYMT3, ISYMW
      INTEGER LU3SRTR, LUCKJDR, LUDELDR, LUDKBCR, LUDKBCR4
      INTEGER KT1AM, KT2TP, KRBJIA, KFOCKD, KEND0, LWRK0, IOPT
      INTEGER IDLSTB, ISYMB, KT1B, KT2B
      INTEGER ISINTR1, ISINTR2, KLAMDP, KLAMDH, KEND1, LWRK1
      INTEGER KTROC, KTROC1, KTROC0, KTROC02, KTROC0Y
      INTEGER KINTOC, KEND2, LWRK2
      INTEGER IOFF, ISYMD, ISCKB1, ISCKB2, ISCKB2Y
      INTEGER KTRVI, KTRVI1, KTRVI2, KTRVI3, KTRVI0, KTRVI0Y
      INTEGER KEND3, LWRK3, KEND4, LWRK4, KINTVI, LENSQ, LENSQW
      INTEGER ISYALJ, ISYALJ2, ISYMBD, ISCKIJ, ISCKD2, ISCKD2Y
      INTEGER ISWMAT, KSMAT, KSMAT3, KUMAT, KUMAT3, KDIAG, KDIAGW
      INTEGER KWMAT, KINDSQ, KINDSQW, KINDEX, KINDEX2, KTMAT
      INTEGER KTRVI8, KTRVI8Y, KTRVI9, KTRVI10, KEND5, LWRK5
      INTEGER LUCKJD, LUTOC, LU3VI, LU3VI2, LUDKBC, LUDELD
      INTEGER KOMEG2, KSAVE, KEND00, LWRK00, KEND01, LWRK01
      INTEGER ITRAN, IDLSTC, ISYMRB, ISYMRC, ISINT1C, ISYRES
      INTEGER KLAMPB, KLAMHB, KFCKB, LUFCK, IVEC, IDLSTL, ISYML
      INTEGER KLAMPC, KLAMHC, KT1C, ISYFCKB, ISYFCKC, KFCKC
      INTEGER IRREP, IERR, KXIAJB, IOPTTCME, ISYMK, KOFF1, KOFF2
      INTEGER ISYMC, IRREP2, ISCKB1C, ISYCKB
C
C     Functions :
C
      INTEGER ILSTSYM
C
COMMENT COMENT
COMMENT COMENT
      INTEGER kx3am
COMMENT COMENT
COMMENT COMENT
C
#if defined (SYS_CRAY)
#else
      DOUBLE PRECISION DOTPROD(MXVEC,NBTRAN), WORK(LWORK)
      DOUBLE PRECISION FREQB, FREQC, TCON
      DOUBLE PRECISION XNORM, DDOT, HALF, TWO
#endif
      CHARACTER*(*) FNCKJD, FNDKBC, FNDELD, FNTOC, FN3VI, FN3VI2
      CHARACTER*12 FN3SRTR, FNCKJDR, FNDELDR, FNDKBCR, FNDKBCR4
      CHARACTER*10 MODEL
      CHARACTER*8 LABELB, LABELC
      CHARACTER*3 LISTL, LISTB, LISTC
      CHARACTER CDUMMY*1
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
C
      PARAMETER(HALF = 0.5D0, TWO = 2.0D0)
      PARAMETER(FN3SRTR  = 'CCSDT_FBMAT1', FNCKJDR  = 'CCSDT_FBMAT2', 
     *          FNDELDR  = 'CCSDT_FBMAT3', FNDKBCR  = 'CCSDT_FBMAT4',
     *          FNDKBCR4 = 'CCSDT_FBMAT5')
C
      CALL QENTER('CCSDT_FBMAT')
C
C--------------------------
C     Save and set flags.
C--------------------------
C
C     Set symmetry flags.
C
C
C     ISYMT2 is symmetry of T2TP
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
C
C--------------------------------
C     Open files
C--------------------------------
C
      LUCKJD   = 0
      LUTOC    = 0
      LU3VI    = 0
      LU3VI2   = 0
      LUDKBC   = 0
      LUDELD   = 0
C
      CALL WOPEN2(LUCKJD,FNCKJD,64,0)
      CALL WOPEN2(LUTOC,FNTOC,64,0)
      CALL WOPEN2(LU3VI,FN3VI,64,0)
      CALL WOPEN2(LU3VI2,FN3VI2,64,0)
      CALL WOPEN2(LUDKBC,FNDKBC,64,0)
      CALL WOPEN2(LUDELD,FNDELD,64,0)
C
C--------------------------------
C     Open temporary files
C--------------------------------
C
      LU3SRTR  = 0
      LUCKJDR  = 0
      LUDELDR  = 0
      LUDKBCR  = 0
      LUDKBCR4 = 0
C
      CALL WOPEN2(LU3SRTR,FN3SRTR,64,0)
      CALL WOPEN2(LUCKJDR,FNCKJDR,64,0)
      CALL WOPEN2(LUDELDR,FNDELDR,64,0)
      CALL WOPEN2(LUDKBCR,FNDKBCR,64,0)
      CALL WOPEN2(LUDKBCR4,FNDKBCR4,64,0)
C
C----------------------------------------------
C     Calculate the zeroth order stuff once
C----------------------------------------------
C
      KT2TP  = 1 
      KFOCKD = KT2TP  + NT2SQ(ISYMT2)
      KLAMDP = KFOCKD + NORBTS
      KLAMDH = KLAMDP + NLAMDT
      KXIAJB = KLAMDH + NLAMDT
      KEND00 = KXIAJB + NT2AM(ISINT1)
      LWRK00 = LWORK  - KEND00
C
      KT1AM  = KEND00
      KEND01 = KT1AM  + NT1AM(ISYMT2)
      LWRK01 = LWORK  - KEND01
C
      IF (LWRK01 .LT. 0) THEN
         CALL QUIT('Out of memory in CCSDT_FBMAT (zeroth allo.')
      ENDIF
C
C------------------------
C     Construct L(ia,jb).
C------------------------
C
      REWIND(LUIAJB)
      CALL READI(LUIAJB,IRAT*NT2AM(ISINT1),WORK(KXIAJB))
      IOPTTCME = 1
      CALL CCSD_TCMEPK(WORK(KXIAJB),1.0D0,ISINT1,IOPTTCME)
C
      IF ( IPRINT .GT. 55) THEN
         XNORM = DDOT(NT2AM(ISINT1),WORK(KXIAJB),1,
     *                WORK(KXIAJB),1)
         WRITE(LUPRI,*) 'Norm of IAJB ',XNORM
      ENDIF
C
C----------------------------------------------------------------
C     Read t1 and calculate the zero'th order Lambda matrices
C----------------------------------------------------------------
C
      IOPT   = 1
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),DUMMY)
C
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND01),LWRK01)
C
C-------------------------------------------
C     Read in t2 , square it and reorder 
C-------------------------------------------
C
      IOPT = 2
      CALL CC_RDRSP('R0',0,ISYMT2,IOPT,MODEL,DUMMY,WORK(KEND00))
      CALL CC_T2SQ(WORK(KEND00),WORK(KT2TP),ISYMT2)
      IF (LWRK00 .LT. NT2SQ(ISYMT2)) THEN
         CALL QUIT('Not enough memory to construct T2TP in CCSDT_FBMAT') 
      ENDIF
C
      CALL DCOPY(NT2SQ(ISYMT2),WORK(KT2TP),1,WORK(KEND00),1)
      CALL CC3_T2TP(WORK(KT2TP),WORK(KEND00),ISYMT2)
C
      IF (IPRINT .GT. 55) THEN
         XNORM = DDOT(NT2SQ(ISYMT2),WORK(KT2TP),1,WORK(KT2TP),1)
         WRITE(LUPRI,*) 'Norm of T2TP ',XNORM
      ENDIF
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
     *   CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND00),LWRK00)
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT If we want to sum the T3 amplitudes
COMMENT COMMENT COMMENT
      if (.false.) then
         kx3am  = kend00
         kend00 = kx3am + nt1amx*nt1amx*nt1amx
         call dzero(work(kx3am),nt1amx*nt1amx*nt1amx)
         lwrk00 = lwork - kend00
         if (lwrk00 .lt. 0) then
            write(lupri,*) 'Memory available : ',lwork
            write(lupri,*) 'Memory needed    : ',kend00
            call quit('Insufficient space (T3) in CCSDT_FBMAT')
         END IF
      endif
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C----------------------------------------------------
C     Loop over left and right amplitude vectors:
C----------------------------------------------------
C
      DO ITRAN = 1, NBTRAN
C
         IDLSTB = IBTRAN(1,ITRAN)
         IDLSTC = IBTRAN(2,ITRAN)
C
         IF (LISTB(1:3).EQ.'R1 ') THEN
            ISYMRB = ISYLRT(IDLSTB)
            FREQB  = FRQLRT(IDLSTB)
            LABELB = LRTLBL(IDLSTB)
            ISYMRC = ISYLRT(IDLSTC)
            FREQC  = FRQLRT(IDLSTC)
            LABELC = LRTLBL(IDLSTC)
         ELSE
            CALL QUIT('Unkown list in CCSDT_FBC_NODDY.')
         END IF
C
         ISYMW   = MULD2H(ISYMT3,ISYMRB)
         ISINT1C = MULD2H(ISINT1,ISYMRC)
         ISYRES  = MULD2H(ISYMW,ISINT1C)
C
         ISYFCKB = MULD2H(ISYMOP,ISYMRB)
         ISYFCKC = MULD2H(ISYMOP,ISYMRC)
C
C-------------------------------------------------
C        Read T1^B and T2^B
C        Calculate (ck|de)-tilde(B) and (ck|lm)-tilde(B)
C        Used to construct T3^B
C-------------------------------------------------
C
         KT2B   = KEND00
         KFCKB  = KT2B   + NT2SQ(ISYMRB)
         KFCKC  = KFCKB  + N2BST(ISYFCKB)
         KEND0  = KFCKC  + N2BST(ISYFCKC)
         LWRK0  = LWORK  - KEND0
C
         KT1B   = KEND0
         KLAMPB = KT1B   + NT1AM(ISYMRB)
         KLAMHB = KLAMPB + NLAMDT
         KEND1  = KLAMHB + NLAMDT
         LWRK1  = LWORK  - KEND1
C
         IF (LWRK1 .LT. NT2SQ(ISYMRB)) THEN
            CALL QUIT('Out of memory in CC3_XI (TOT_T3Y) ')
         ENDIF
C
C        Readin
C
         IOPT = 3
         CALL CC_RDRSP(LISTB,IDLSTB,ISYMRB,IOPT,MODEL,
     *                 WORK(KT1B),WORK(KT2B))
         CALL CCLR_DIASCL(WORK(KT2B),TWO,ISYMRB)
         CALL CC_T2SQ(WORK(KT2B),WORK(KEND1),ISYMRB)
         CALL CC3_T2TP(WORK(KT2B),WORK(KEND1),ISYMRB)
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NT1AM(ISYMRB),WORK(KT1B),1,WORK(KT1B),1)
            WRITE(LUPRI,*) 'Norm of T1B  ',XNORM
            XNORM = DDOT(NT2SQ(ISYMRB),WORK(KT2B),1,WORK(KT2B),1)
            WRITE(LUPRI,*) 'Norm of T2B  ',XNORM
         ENDIF
C
C        Integrals
C
         ISINTR1 = MULD2H(ISINT1,ISYMRB)
         ISINTR2 = MULD2H(ISINT2,ISYMRB)
C
         CALL CC3_BARINT(WORK(KT1B),ISYMRB,WORK(KLAMDP),
     *                   WORK(KLAMDH),WORK(KEND1),LWRK1,
     *                   LU3SRTR,FN3SRTR,LUCKJDR,FNCKJDR)
C
         CALL CC3_SORT1(WORK(KEND1),LWRK1,2,ISINTR1,LU3SRTR,FN3SRTR,
     *                  LUDELDR,FNDELDR,IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *                  IDUMMY,CDUMMY)
C
         CALL CC3_SINT(WORK(KLAMDH),WORK(KEND1),LWRK1,ISINTR1,
     *                 LUDELDR,FNDELDR,LUDKBCR,FNDKBCR)
C
         CALL CC3_TCME(WORK(KLAMDP),ISINTR1,WORK(KEND1),LWRK1,
     *                 IDUMMY,CDUMMY,LUDKBCR,FNDKBCR,
     *                 IDUMMY,CDUMMY,IDUMMY,CDUMMY,
     *                 IDUMMY,CDUMMY,LUDKBCR4,FNDKBCR4,2)
C
C---------------------------------------
C        Calculate the Lambda^B
C---------------------------------------
C
         CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLAMPB),WORK(KLAMDH),
     *                    WORK(KLAMHB),WORK(KT1B),ISYMRB)
C
C------------------------------------------
C        Calculate the F^B matrix
C------------------------------------------
C
         CALL CCPRPAO(LABELB,WORK(KFCKB),ISYFCKB,IRREP,IERR,
     *                WORK(KEND1),LWRK1)
C
         IF (IERR.GT.0) THEN
           CALL QUIT('CCSDT_FBMAT : error reading operator '//LABELB)
         ELSE IF (IERR.LT.0) THEN
           CALL DZERO(WORK(KFCKB),N2BST(ISYMRB))
         END IF
C
         CALL CC_FCKMO(WORK(KFCKB),WORK(KLAMDP),WORK(KLAMDH),
     *                 WORK(KEND1),LWRK1,ISYMRB,ISYMOP,ISYMOP)
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND( 'In CCSDT_FBMAT: Fock^B MO matrix' )
            CALL CC_PRFCKMO(WORK(KFCKB),ISYFCKB)
         ENDIF
C
C---------------------------------------------
C        Regain work space and create lam^C
C---------------------------------------------
C
         KT1C   = KEND0
         KLAMPC = KT1C   + NT1AM(ISYMRC)
         KLAMHC = KLAMPC + NLAMDT
         KEND0  = KLAMHC + NLAMDT
         LWRK0  = LWORK  - KEND0
C
         IOPT = 1
         CALL CC_RDRSP(LISTC,IDLSTC,ISYMRC,IOPT,MODEL,
     *                 WORK(KT1C),DUMMY)
C
         CALL CCLR_LAMTRA(WORK(KLAMDP),WORK(KLAMPC),WORK(KLAMDH),
     *                    WORK(KLAMHC),WORK(KT1C),ISYMRC)
C
C------------------------------------------
C        Calculate the F^C matrix
C------------------------------------------
C
         CALL CC3LR_MFOCK(WORK(KEND1),WORK(KT1C),WORK(KXIAJB),ISYFCKC)
C
C     Put the F_{kc} part into correct F_{pq}
C
         CALL DZERO(WORK(KFCKC),N2BST(ISYFCKC))
C
         DO ISYMC = 1, NSYM
            ISYMK = MULD2H(ISYFCKC,ISYMC)
            DO K = 1, NRHF(ISYMK)
               DO C = 1, NVIR(ISYMC)
                  KOFF1 = KFCKC - 1
     *                  + IFCVIR(ISYMK,ISYMC)
     *                  + NORB(ISYMK)*(C - 1)
     *                  + K
                  KOFF2 = KEND1 - 1
     *                  + IT1AM(ISYMC,ISYMK)
     *                  + NVIR(ISYMC)*(K - 1)
     *                  + C
C
                  WORK(KOFF1) = WORK(KOFF2)
C
               ENDDO
            ENDDO
         ENDDO
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND( 'In CCSDT_FBMAT: Fock^C MO matrix' )
            CALL CC_PRFCKMO(WORK(KFCKC),ISYFCKC)
         ENDIF
C
C
C-------------------------------------
C        Summation over vectors!
C-------------------------------------
C
         IVEC = 1
         DO WHILE (IDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
C
         IDLSTL = IDOTS(IVEC,ITRAN)
         ISYML  = ILSTSYM(LISTL,IDLSTL)
C
C-----------------------------
C        Read occupied integrals.
C-----------------------------
C
C        Memory allocation.
C
C---------------------------------------------------------
C        Work allocation 
C---------------------------------------------------------
C
         KOMEG2 = KEND0
         KRBJIA = KOMEG2 + NT2AM(ISYRES)
         KTROC  = KRBJIA + NT2SQ(ISYRES)
         KTROC1 = KTROC  + NTRAOC(ISINT1C)
         KTROC0 = KTROC1 + NTRAOC(ISINT1C)
         KTROC02= KTROC0 + NTRAOC(ISINT2)
         KEND1  = KTROC02+ NTRAOC(ISINT2)
         LWRK1  = LWORK  - KEND1
C
         KTROC0Y = KEND1
         KEND1   = KTROC0Y + NTRAOC(ISINTR2)
C
         KINTOC = KEND1
         KEND2  = KINTOC 
     *          + MAX(NTOTOC(ISINT2),NTOTOC(ISINT1),NTOTOC(ISYMOP))
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Memory available : ',LWORK
            WRITE(LUPRI,*) 'Memory needed    : ',KEND2
            CALL QUIT('Insufficient space in CC3_XI')
         END IF
C
C----------------------------------------
C        Initialise the result vectors
C----------------------------------------
C
         CALL DZERO(WORK(KRBJIA),NT2SQ(ISYRES))
         CALL DZERO(WORK(KOMEG2),NT2AM(ISYRES))
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
            XNORM = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                   WORK(KINTOC),1)
            WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XNORM
         ENDIF
C
C----------------------------------------------------------------------
C     Transform (ai|j delta) integrals to (ai|j k) and sort as (ij,k,a)
C----------------------------------------------------------------------
C
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0),WORK(KLAMDP),
     *                  WORK(KEND2),LWRK2,ISINT2)
C
C------------------------------------------
C     B transformed Occupied integrals.
C-----------------------------------------
C
         IOFF = 1
         IF (NTOTOC(ISINTR2) .GT. 0) THEN
            CALL GETWA2(LUCKJDR,FNCKJDR,WORK(KINTOC),IOFF,
     *                  NTOTOC(ISINTR2))
         ENDIF
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NTOTOC(ISINTR2),WORK(KINTOC),1,
     *                   WORK(KINTOC),1)
            WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT (B transformed) ',XNORM
         ENDIF
C
         CALL CC3_TROCC(WORK(KINTOC),WORK(KTROC0Y),WORK(KLAMDP),
     *                  WORK(KEND2),LWRK2,ISINTR2)
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
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NTOTOC(ISINT2),WORK(KINTOC),1,
     *                   WORK(KINTOC),1)
            WRITE(LUPRI,*) 'Norm of CKJDEL-INT  ',XNORM
         ENDIF
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NTRAOC(ISINT2),WORK(KTROC0),1,
     *                   WORK(KTROC0),1)
            WRITE(LUPRI,*) 'Norm of TROC0 ',XNORM
         ENDIF
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NTRAOC(ISINT2),WORK(KTROC02),1,
     *                   WORK(KTROC02),1)
            WRITE(LUPRI,*) 'Norm of TROC02 ',XNORM
         ENDIF
C
C------------------------
C     Occupied integrals.
C
C     Read in integrals for WMAT etc.
C-----------------------
C
         IOFF = 1
         IF (NTOTOC(ISYMOP) .GT. 0) THEN
            CALL GETWA2(LUTOC,FNTOC,WORK(KINTOC),IOFF,NTOTOC(ISYMOP))
         ENDIF
C
C----------------------------------
C     Write out norms of Integrals.
C----------------------------------
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NTOTOC(ISYMOP),WORK(KINTOC),1,
     *                   WORK(KINTOC),1)
            WRITE(LUPRI,*) 'Norm of CCSDT_OC-INT ',XNORM
         ENDIF
C
C----------------------------------------------------------------------
C     Transform (ia|j delta) integrals to (ia|j k) and sort as (i,j,k,a)
C----------------------------------------------------------------------
C
          CALL CCLR_TROCC(WORK(KINTOC),WORK(KTROC),
     *                     WORK(KLAMHC),ISYMRC,
     *                     WORK(KEND2),LWRK2)
C
C----------------------------------------------------------------------
C     sort (i,j,k,a) as (a,i,j,k)
C----------------------------------------------------------------------
C
         CALL CCSDT_SRTOC2(WORK(KTROC),WORK(KTROC1),ISINT1C,
     *                     WORK(KEND2),LWRK2)
C
C----------------------------
C     General loop structure.
C----------------------------
C
         DO ISYMD = 1,NSYM
C
            ISYCKB  = MULD2H(ISYMOP,ISYMD)
            ISCKB2  = MULD2H(ISINT2,ISYMD)
            ISCKB2Y = MULD2H(ISINTR2,ISYMD)
            ISCKB1  = MULD2H(ISINT1,ISYMD)
            ISCKB1C = MULD2H(ISINT1C,ISYMD)
C
            IF (IPRINT .GT. 55) THEN
C
               WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISYCKB :',ISYCKB
               WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISCKB2 :',ISCKB2
               WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISCKB2Y:',ISCKB2Y
               WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISCKB1 :',ISCKB1
               WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISCKB1C:',ISCKB1C
C
            ENDIF

C
C--------------------------
C        Memory allocation.
C--------------------------
C
            KTRVI  = KEND1
            KTRVI1 = KTRVI  + NCKATR(ISCKB1C)
            KTRVI3 = KTRVI1 + NCKATR(ISCKB1C)
            KTRVI2 = KTRVI3 + NCKATR(ISCKB2)
            KEND2  = KTRVI2 + NCKATR(ISCKB2)
            LWRK2  = LWORK  - KEND2
C
            KTRVI0  = KEND2
            KEND3   = KTRVI0  + NCKATR(ISCKB2)
            LWRK3   = LWORK  - KEND3
C
            KTRVI0Y  = KEND3
            KEND3    = KTRVI0Y  + NCKATR(ISCKB2Y)
            LWRK3    = LWORK    - KEND3
C
            KINTVI = KEND3
            KEND4  = KINTVI 
     *             + MAX(NCKA(ISCKB2),NCKA(ISCKB1),NCKA(ISYCKB))
            LWRK4  = LWORK  - KEND4
C
            IF (LWRK4 .LT. 0) THEN
               WRITE(LUPRI,*) 'Memory available : ',LWORK
               WRITE(LUPRI,*) 'Memory needed    : ',KEND4
               CALL QUIT('Insufficient space in CC3_XI')
            END IF
C
C---------------------
C        Sum over D
C---------------------
C
            DO D = 1,NVIR(ISYMD)
C
C------------------------------------------------------------
C           Read and transform integrals used in contraction
C           with W intermediate.
C------------------------------------------------------------
C
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LU3VI2,FN3VI2,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
               ENDIF
               CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI),
     *                          WORK(KLAMPC),ISYMRC,
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
               IOFF = ICKAD(ISYCKB,ISYMD) + NCKA(ISYCKB)*(D - 1) + 1
               IF (NCKA(ISYCKB) .GT. 0) THEN
                  CALL GETWA2(LU3VI,FN3VI,WORK(KINTVI),IOFF,
     &                        NCKA(ISYCKB))
               ENDIF
               CALL CCLR_TRVIR(WORK(KINTVI),WORK(KTRVI1),
     *                          WORK(KLAMPC),ISYMRC,
     *                          ISYMD,D,ISYMOP,WORK(KEND4),LWRK4)
C
C-----------------------------------------------
C           Read virtual integrals used in first s3am.
C-----------------------------------------------
C
               IOFF = ICKBD(ISCKB2,ISYMD) + NCKATR(ISCKB2)*(D - 1) + 1
               IF (NCKATR(ISCKB2) .GT. 0) THEN
                  CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI0),IOFF,
     &                           NCKATR(ISCKB2))
               ENDIF
C
C-----------------------------------------------------------------------
C           Read B transformed virtual integrals used for W for TOT_T3Y
C-----------------------------------------------------------------------
C
               IOFF = ICKBD(ISCKB2Y,ISYMD) + NCKATR(ISCKB2Y)*(D - 1) + 1
               IF (NCKATR(ISCKB2) .GT. 0) THEN
                  CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI0Y),IOFF,
     &                        NCKATR(ISCKB2Y))
               ENDIF
C
C-----------------------------------------------------------
C           Sort the integrals for s3am 
C-----------------------------------------------------------
C
               CALL CCSDT_SRTVIR(WORK(KTRVI0),WORK(KTRVI2),WORK(KEND4),
     *                           LWRK4,ISYMD,ISINT2)
C
               IF (IPRINT .GT. 55) THEN
                  XNORM = DDOT(NCKATR(ISCKB2),WORK(KTRVI0),1,
     *                         WORK(KTRVI0),1)
                  WRITE(LUPRI,*) 'Norm of TRVI0 ',XNORM
               ENDIF
C
               IF (IPRINT .GT. 55) THEN
                  XNORM = DDOT(NCKATR(ISCKB2),WORK(KTRVI2),1,
     *                         WORK(KTRVI2),1)
                  WRITE(LUPRI,*) 'Norm of TRVI2 ',XNORM
               ENDIF
C
C------------------------------------------------------
C           Read virtual integrals used in first U.
C------------------------------------------------------
C
               IOFF = ICKAD(ISCKB2,ISYMD) + NCKA(ISCKB2)*(D - 1) + 1
               IF (NCKA(ISCKB2) .GT. 0) THEN
                  CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     &                        NCKA(ISCKB2))
               ENDIF
C
               CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI3),WORK(KLAMDH),
     *                          ISYMD,D,ISINT2,WORK(KEND4),LWRK4)
C
C--------------------------------------------------------
C              Sum over ISYMB
C--------------------------------------------------------
C
               DO ISYMB = 1,NSYM
C
                  ISYALJ  = MULD2H(ISYMB,ISYMT2)
                  ISYALJ2 = MULD2H(ISYMD,ISYMT2)
                  ISYMBD  = MULD2H(ISYMB,ISYMD)
                  ISCKIJ  = MULD2H(ISYMBD,ISYMT3)
                  ISCKD2  = MULD2H(ISINT2,ISYMB)
                  ISCKD2Y = MULD2H(ISINTR2,ISYMB)
                  ISWMAT  = MULD2H(ISYMRB,ISCKIJ)
C
                  IF (IPRINT .GT. 55) THEN
C 
                     WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISYMD  :',ISYMD
                     WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISYMB  :',ISYMB
                     WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISYALJ :',ISYALJ
                     WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISYALJ2:',ISYALJ2
                     WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISYMBD :',ISYMBD
                     WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISCKIJ :',ISCKIJ
                     WRITE(LUPRI,*) 'In CCSDT_FBMAT: ISWMAT :',ISWMAT
C
                  ENDIF
C
C              Can use kend3 since we do not need the integrals anymore.
                  KSMAT   = KEND3
                  KSMAT3  = KSMAT   + NCKIJ(ISCKIJ)
                  KUMAT   = KSMAT3  + NCKIJ(ISCKIJ)
                  KUMAT3  = KUMAT   + NCKIJ(ISCKIJ)
                  KDIAG   = KUMAT3  + NCKIJ(ISCKIJ)
                  KDIAGW  = KDIAG   + NCKIJ(ISCKIJ)
                  KWMAT   = KDIAGW  + NCKIJ(ISWMAT)
                  KINDSQW = KWMAT   + NCKIJ(ISWMAT)
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
                  KTRVI8Y = KEND4
                  KEND4   = KTRVI8Y + NCKATR(ISCKD2Y)
C
C
                  KINTVI  = KEND4
                  KEND5   = KINTVI  + NCKA(ISCKD2)
                  LWRK5   = LWORK   - KEND5
C
                  IF (LWRK5 .LT. 0) THEN
                     WRITE(LUPRI,*) 'Memory available : ',LWORK
                     WRITE(LUPRI,*) 'Memory needed    : ',KEND5
                     CALL QUIT('Insufficient space in CC3_XI')
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
                     XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KDIAG),1,
     *                       WORK(KDIAG),1)
                     WRITE(LUPRI,*) 'Norm of DIA  ',XNORM
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
                     CALL DZERO(WORK(KWMAT),NCKIJ(ISWMAT))
C
C-----------------------------------------------
C                 Read virtual integrals used in second s3am.
C-----------------------------------------------
C
                     IOFF = ICKBD(ISCKD2,ISYMB) + 
     &                      NCKATR(ISCKD2)*(B - 1) + 1
                     IF (NCKATR(ISCKD2) .GT. 0) THEN
                        CALL GETWA2(LUDKBC,FNDKBC,WORK(KTRVI8),IOFF,
     &                              NCKATR(ISCKD2))
                     ENDIF
C
C--------------------------------------------------------------------
C                 Read B transformed virtual integrals used in W
C--------------------------------------------------------------------
C
                     IOFF = ICKBD(ISCKD2Y,ISYMB) + 
     &                      NCKATR(ISCKD2Y)*(B - 1) + 1
                     IF (NCKATR(ISCKD2) .GT. 0) THEN
                        CALL GETWA2(LUDKBCR,FNDKBCR,WORK(KTRVI8Y),IOFF,
     &                              NCKATR(ISCKD2Y))
                     ENDIF
C
C-----------------------------------------------------------
C           Sort the integrals for s3am 
C-----------------------------------------------------------
C
                     CALL CCSDT_SRTVIR(WORK(KTRVI8),
     *                                 WORK(KTRVI9),WORK(KEND4),
     *                                 LWRK4,ISYMB,ISINT2)
C
C
C----------------------------------------------------------
C           Read virtual integrals used in second U
C----------------------------------------------------------
C
                    IOFF = ICKAD(ISCKD2,ISYMB) 
     *                   + NCKA(ISCKD2)*(B - 1) + 1
                    IF (NCKA(ISCKD2) .GT. 0) THEN
                       CALL GETWA2(LUDELD,FNDELD,WORK(KINTVI),IOFF,
     *                             NCKA(ISCKD2))
                    ENDIF
C
                    CALL CCSDT_TRVIR(WORK(KINTVI),WORK(KTRVI10),
     *                               WORK(KLAMDH),ISYMB,B,ISINT2,
     &                               WORK(KEND5),LWRK5)
C
C------------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for B,D.
C-------------------------------------------------------------------
C
                     CALL CC3_SMAT(WORK(KT2TP),ISYMT2,
     *                             WORK(KTMAT),WORK(KTRVI0),
     *                             WORK(KTRVI2),WORK(KTROC0),ISINT2,
     *                             WORK(KFOCKD),WORK(KDIAG),
     *                             WORK(KSMAT),WORK(KEND4),LWRK4,
     *                             WORK(KINDEX),WORK(KINDSQ),LENSQ,
     *                             ISYMB,B,ISYMD,D)
C
                     CALL T3_FORBIDDEN(WORK(KSMAT),ISYMT3,ISYMB,B,
     *                                 ISYMD,D)
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C       call sum_pt3(work(ksmat),isymb,b,isymd,d,
C     *             isckij,work(kx3am),1)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
                     IF (IPRINT .GT. 55) THEN
                        XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT),1,
     *                          WORK(KSMAT),1)
                        WRITE(LUPRI,*) 'Norm of SMAT     ',XNORM
                     ENDIF
C
C-------------------------------------------------------------------
C                 Calculate the S(ci,bk,dj) matrix for T3 for D,B.
C-------------------------------------------------------------------
C
                     CALL CC3_SMAT(WORK(KT2TP),ISYMT2,
     *                             WORK(KTMAT),WORK(KTRVI8),
     *                             WORK(KTRVI9),WORK(KTROC0),ISINT2,
     *                             WORK(KFOCKD),WORK(KDIAG),
     *                             WORK(KSMAT3),WORK(KEND4),LWRK4,
     *                             WORK(KINDEX2),WORK(KINDSQ),LENSQ,
     *                             ISYMD,D,ISYMB,B)
C
                     CALL T3_FORBIDDEN(WORK(KSMAT3),ISYMT3,
     *                                 ISYMD,D,ISYMB,B)
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C      call sum_pt3(work(ksmat3),isymd,d,isymb,b,
C     *             isckij,work(kx3am),1)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
                     IF (IPRINT .GT. 55) THEN
                        XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KSMAT3),1,
     *                          WORK(KSMAT3),1)
                        WRITE(LUPRI,*) 'Norm of SMAT3    ',XNORM
                     ENDIF
C
C---------------------------------------------------------------------------
C                 Calculate U(ci,jk) for fixed b,d.
C--------------------------------------------------
C
                     CALL CC3_UMAT(WORK(KT2TP),ISYMT2,WORK(KTRVI3),
     *                             WORK(KTROC02),ISINT2,WORK(KFOCKD),
     *                             WORK(KDIAG),WORK(KUMAT),WORK(KTMAT),
     *                             WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                             ISYMB,B,ISYMD,D)
C
                     CALL T3_FORBIDDEN(WORK(KUMAT),ISYMT3,
     *                                 ISYMB,B,ISYMD,D)
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C       call sum_pt3(work(kumat),isymb,b,isymd,d,
C     *             isckij,work(kx3am),3)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
                     IF (IPRINT .GT. 55) THEN
                        XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KUMAT),1,
     *                          WORK(KUMAT),1)
                        WRITE(LUPRI,*) 'Norm of UMAT     ',XNORM
                     ENDIF
C
C--------------------------------------------------
C                 Calculate U(ci,jk) for fixed d,b.
C--------------------------------------------------
C
                     CALL CC3_UMAT(WORK(KT2TP),ISYMT2,WORK(KTRVI10),
     *                             WORK(KTROC02),ISINT2,WORK(KFOCKD),
     *                             WORK(KDIAG),WORK(KUMAT3),WORK(KTMAT),
     *                             WORK(KEND4),LWRK4,WORK(KINDSQ),LENSQ,
     *                             ISYMD,D,ISYMB,B)
C
                     CALL T3_FORBIDDEN(WORK(KUMAT3),ISYMT3,
     *                                 ISYMD,D,ISYMB,B)
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C      call sum_pt3(work(kumat3),isymd,d,isymb,b,
C     *             isckij,work(kx3am),3)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
                     IF (IPRINT .GT. 55) THEN
                        XNORM = DDOT(NCKIJ(ISCKIJ),WORK(KUMAT),1,
     *                          WORK(KUMAT),1)
                        WRITE(LUPRI,*) 'Norm of UMAT3    ',XNORM
                     ENDIF
C
C--------------------------------------------------
C Sum up S and U intermediates to get T3 BD amplitudes
C--------------------------------------------------
C
                     CALL CC3_T3BD(ISCKIJ,WORK(KSMAT),WORK(KSMAT3),
     *                             WORK(KUMAT),WORK(KUMAT3),WORK(KTMAT),
     *                             WORK(KINDSQ),LENSQ)
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
                     CALL WBD_V(WORK(KTMAT),ISCKIJ,WORK(KFCKB),ISYMRB,
     *                               WORK(KWMAT), 
     *                               ISWMAT,WORK(KEND4),LWRK4)  
C
C------------------------------------------------------
C     Calculate the  term <mu3|[Y,T3]|HF> occupied contribution 
C     added in W^BD(KWMAT)
C------------------------------------------------------
C
                     CALL WBD_O(WORK(KTMAT),ISCKIJ,WORK(KFCKB),ISYMRB,
     *                               WORK(KWMAT),
     *                               ISWMAT,WORK(KEND4),LWRK4)
C
C------------------------------------------------------
C     Calculate the  term <mu3|[[Y,T2],T2]|HF> 
C     added in W^BD(KWMAT)
C------------------------------------------------------
C
                     CALL WBD_T2(B,ISYMB,D,ISYMD,WORK(KT2TP),ISYMT2,
     *                           WORK(KFCKB),ISYMRB,
     *                           WORK(KINDSQW),LENSQW,WORK(KWMAT), 
     *                           ISWMAT,WORK(KEND4),LWRK4) 
C
C------------------------------------------------------
C     To get the entire T3^Y add the two terms
C------------------------------------------------------
C
                     CALL WBD_GROUND(WORK(KT2B),ISYMRB,WORK(KTMAT),
     *                               WORK(KTRVI0),WORK(KTRVI8),
     *                               WORK(KTROC0),ISINT2,WORK(KWMAT),
     *                               WORK(KEND4),LWRK4,
     *                               WORK(KINDSQW),LENSQW,
     *                               ISYMB,B,ISYMD,D)
C
                     CALL WBD_GROUND(WORK(KT2TP),ISYMT2,WORK(KTMAT),
     *                               WORK(KTRVI0Y),WORK(KTRVI8Y),
     *                               WORK(KTROC0Y),ISINTR2,WORK(KWMAT),
     *                               WORK(KEND4),LWRK4,
     *                               WORK(KINDSQW),LENSQW,
     *                               ISYMB,B,ISYMD,D)
C
C------------------------------------------------
C     Divide by the energy difference and
C     remove the forbidden elements
C------------------------------------------------
C
                     CALL WBD_DIA(B,ISYMB,D,ISYMD,FREQB, 
     *                            ISWMAT,WORK(KWMAT), 
     *                            WORK(KDIAGW),WORK(KFOCKD)) 
C
                     CALL T3_FORBIDDEN(WORK(KWMAT),ISYMRB,
     *                                 ISYMB,B,ISYMD,D) 
C
C----------------------------------------------------------------------
C     Calculate the  term <mu2|[H,W^BD(3)]|HF> ( Fock matrix cont ) 
C     added in WORK(KOMEG2)
C----------------------------------------------------------------------
C
                     CALL CC3_CY2F(WORK(KOMEG2),ISYRES,WORK(KWMAT),
     *                             ISWMAT,WORK(KTMAT),WORK(KFCKC),
     *                             ISYFCKC,WORK(KINDSQW),LENSQW,
     *                             WORK(KEND4),LWRK4,
     *                             ISYMB,B,ISYMD,D)
C
C------------------------------------------------------
C     Calculate the  term <mu2|[H,W^BD(3)]|HF> ( Occupied  cont ) 
C     added in WORK(KOMEG2) 
C------------------------------------------------------
C
                    CALL CC3_CY2O(WORK(KOMEG2),ISYRES,WORK(KWMAT),
     *                            ISWMAT,WORK(KTMAT),WORK(KTROC),
     *                            WORK(KTROC1),ISINT1C,
     *                            WORK(KEND4),LWRK4,
     *                            WORK(KINDSQW),LENSQW,
     *                            ISYMB,B,ISYMD,D)
C
C------------------------------------------------------
C     Calculate the  term <mu2|[H,W^BD(3)]|HF> ( Virtual cont ) 
C     added in WORK(KOMEG2) 
C------------------------------------------------------
C
                     CALL CC3_CY2V(WORK(KOMEG2),ISYRES,WORK(KRBJIA),
     *                             WORK(KWMAT),ISWMAT,WORK(KTMAT),
     *                             WORK(KTRVI),WORK(KTRVI1),ISINT1C,
     *                             WORK(KEND4),LWRK4,WORK(KINDSQW),
     *                             LENSQW,ISYMB,B,ISYMD,D)
C 
                  END DO  ! B
               END DO     ! ISYMB
C
            END DO        ! D
         END DO           ! ISYMD
C
COMMENT COMMENT
COMMENT COMMENT
C      write(lupri,*) 'The summed W terms : '
C      call print_pt3(work(kx3am),1,4)
C      write(lupri,*) 'The summed S terms : '
C      call print_pt3(work(kx3am),1,1)
COMMENT COMMENT
COMMENT COMMENT
C
C
C------------------------------------------------------
C     Accumulate RBJIA from <mu2|[H,W^BD(3)]|HF> ( Vccupied  cont ) 
C     in WORK(KOMEG2) 
C------------------------------------------------------
C
         CALL CC3_RBJIA(WORK(KOMEG2),ISYRES,WORK(KRBJIA))
C
         IF (IPRINT .GT. 55) THEN
            XNORM = DDOT(NT2AM(ISYRES),WORK(KOMEG2),1,WORK(KOMEG2),1)
            WRITE(LUPRI,*) 'Norm of final WORK(KOMEG2) ',XNORM
         ENDIF
C
C--------------------------------------------
C     Multiply the result vector with L2
C     and put the result into DOTPROD
C--------------------------------------------
C
         IOPT = 2
         CALL CC_RDRSP(LISTL,0,1,IOPT,MODEL,DUMMY,WORK(KEND1))
C
         CALL CCLR_DIASCL(WORK(KEND1),HALF,ISYML)
         TCON = DDOT(NT2AM(ISYRES),WORK(KOMEG2),1,WORK(KEND1),1)
C
         DOTPROD(IVEC,ITRAN) = DOTPROD(IVEC,ITRAN)
     *                       + TCON
C
         IF (IPRINT .GT. 11) THEN
            WRITE(LUPRI,*) 'IVEC, ITRAN, Contribution: ',IVEC,ITRAN,TCON
         ENDIF
C
C-------------------------------------
C     End loop over vectors.
C     End loop over left and right
C-------------------------------------
C
         IVEC = IVEC + 1
C
         ENDDO    ! DO WHILE
C
      ENDDO       ! ITRAN
C
C--------------------------------
C     Close files
C--------------------------------
C
      CALL WCLOSE2(LUCKJD,FNCKJD,'KEEP')
      CALL WCLOSE2(LUTOC,FNTOC,'KEEP')
      CALL WCLOSE2(LU3VI,FN3VI,'KEEP')
      CALL WCLOSE2(LU3VI2,FN3VI2,'KEEP')
      CALL WCLOSE2(LUDKBC,FNDKBC,'KEEP')
      CALL WCLOSE2(LUDELD,FNDELD,'KEEP')
C
C--------------------------------
C     Close files for "response"
C--------------------------------
C
      CALL WCLOSE2(LU3SRTR,FN3SRTR,'DELETE')
      CALL WCLOSE2(LUCKJDR,FNCKJDR,'DELETE')
      CALL WCLOSE2(LUDELDR,FNDELDR,'DELETE')
      CALL WCLOSE2(LUDKBCR,FNDKBCR,'DELETE')
      CALL WCLOSE2(LUDKBCR4,FNDKBCR4,'DELETE')
C
C-------------
C     End
C-------------
C
      CALL QEXIT('CCSDT_FBMAT')
C
      RETURN
      END
