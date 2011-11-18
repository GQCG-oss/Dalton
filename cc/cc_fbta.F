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
       SUBROUTINE CC_FBTA( IFBTRAN, NFBTRAN, IOPTRES, LISTL, LISTR, 
     &                     FILFBTA, IFBDOTS,  FBCONS, 
     &                     MXVEC,   WORK,    LWORK           )
*---------------------------------------------------------------------*
*
*    Purpose: batched loop over F^B T^A transformations
*             (to avoid problems with dynamic allocation)
*
*             for more detailed documentation see: CC_FBTA1
*        
*    Sonia Coriani 07/02/2000 based on CC_XIETA
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <maxorb.h>

* local
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER MAXFBTRAN
      PARAMETER (MAXFBTRAN = 100)
* argument-list
      CHARACTER*(*) LISTL, LISTR, FILFBTA
      INTEGER IOPTRES, MXVEC
      INTEGER NFBTRAN, LWORK
      INTEGER IFBTRAN(5,NFBTRAN)
      INTEGER IFBDOTS(MXVEC,NFBTRAN)
      
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL FBCONS(MXVEC,NFBTRAN)
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION FBCONS(MXVEC,NFBTRAN)
#endif 

      INTEGER NTRAN, IBATCH, NBATCH, IDUM
      INTEGER ISTART, IEND, IADRE0, IADRE1
      INTEGER IADRPQ0, IADRPQ1, IADRPQMO 
      INTEGER IADRPQA0, IADRPQA1, IADPQAMO 
      INTEGER IADRBFA, IADRBFQA, IADRBZA, IADRBZQA
      INTEGER IT2DEL0, IT2DELB, IT2DELA,IT2DELBA
      INTEGER KCMOPQ, KCMOHQ, KLAMDPQ, KLAMDHQ
      INTEGER KLAMDPA, KLAMDHA, KLAMPQA, KLAMHQA
      INTEGER KONEHB, KDENSB, KDPCKB
      INTEGER KDENS0A,KDPCK0A,KDENSBA,KDPCKBA
      INTEGER KFOCKBCC, KFOCKACC, KFCKBACC
      INTEGER KRAIM, KRQAIM, KFAIM, KFQAIM
      INTEGER KTHETA1, KTHETA2, KRHOBFA,KRHOBFQA
      INTEGER KXAINT, KYAINT
      INTEGER KBFZI,  KBZADEN, KBZQADEN
      INTEGER KEND, LEND

      NBATCH = (NFBTRAN+MAXFBTRAN-1)/MAXFBTRAN

      IF (NBATCH.GT.1) THEN
        WRITE (LUPRI,*) 
     *           'Batching for CC_FBTA not debugged for multipass case'
        CALL QUIT(
     *           'Batching for CC_FBTA not debugged for multipass case')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Batching over FbTa vector calculations:'
        WRITE (LUPRI,*) 'nb. of batches needed:', NBATCH
      END IF
  
      DO IBATCH = 1, NBATCH
        ISTART = (IBATCH-1) * MAXFBTRAN + 1
        NTRAN  = MIN(NFBTRAN-(ISTART-1),MAXFBTRAN)

        ISTART      = 1
        IEND        = ISTART    + NTRAN
        IADRE0      = IEND      + NTRAN
        IADRE1      = IADRE0    + MXCORB*NTRAN
        IADRPQ0     = IADRE1    + MXCORB*NTRAN
        IADRPQ1     = IADRPQ0   + MXCORB*NTRAN
        IADRPQMO    = IADRPQ1   + MXCORB*NTRAN
        IADRPQA0    = IADRPQMO  + MXCORB*NTRAN
        IADRPQA1    = IADRPQA0  + MXCORB*NTRAN
        IADPQAMO    = IADRPQA1  + MXCORB*NTRAN
        IADRBFA     = IADPQAMO  + MXCORB*NTRAN
        IADRBFQA    = IADRBFA   + MXCORB*NTRAN
        IADRBZA     = IADRBFQA  + MXCORB*NTRAN
        IADRBZQA    = IADRBZA   + MXCORB*NTRAN
C
        IT2DEL0     = IADRBZQA  + MXCORB*NTRAN
        IT2DELB     = IT2DEL0   + MXCORB*1
        IT2DELA     = IT2DELB   + MXCORB*NTRAN
        IT2DELBA    = IT2DELA   + MXCORB*NTRAN
C
        KCMOPQ      = IT2DELBA  + MXCORB*NTRAN
        KCMOHQ      = KCMOPQ    + NTRAN
        KLAMDPQ     = KCMOHQ    + NTRAN
        KLAMDHQ     = KLAMDPQ   + NTRAN
        KLAMDPA     = KLAMDHQ   + NTRAN
        KLAMDHA     = KLAMDPA   + NTRAN
        KLAMPQA     = KLAMDHA   + NTRAN
        KLAMHQA     = KLAMPQA   + NTRAN
C
        KONEHB      = KLAMHQA   + NTRAN 
        KDENSB      = KONEHB    + NTRAN
        KDPCKB      = KDENSB    + NTRAN
        KDENS0A     = KDPCKB    + NTRAN
        KDPCK0A     = KDENS0A   + NTRAN
        KDENSBA     = KDPCK0A   + NTRAN
        KDPCKBA     = KDENSBA   + NTRAN

        KFOCKBCC    = KDPCKBA  + NTRAN
        KFOCKACC    = KFOCKBCC + NTRAN
        KFCKBACC    = KFOCKACC + NTRAN
C
        KRAIM       = KFCKBACC + NTRAN
        KRQAIM      = KRAIM    + NTRAN
        KFAIM       = KRQAIM   + NTRAN
        KFQAIM      = KFAIM    + NTRAN

        KTHETA1     = KFQAIM   + NTRAN
        KTHETA2     = KTHETA1  + NTRAN
        KRHOBFA     = KTHETA2  + NTRAN
        KRHOBFQA    = KRHOBFA  + NTRAN
        KXAINT      = KRHOBFQA + NTRAN
        KYAINT      = KXAINT   + NTRAN
        KBFZI       = KYAINT   + NTRAN
        KBZADEN     = KBFZI    + NTRAN
        KBZQADEN    = KBZADEN  + NTRAN
        KEND        = KBZQADEN + NTRAN
        LEND        = LWORK    - KEND

        IF (LEND .LT. 0) THEN
           WRITE (LUPRI,*) 'Insufficient work space in CC_FBTA.'
           WRITE (LUPRI,*) 'Available    :',LWORK,' words,'
           WRITE (LUPRI,*) 'Need at least:',KEND, ' words.'
           CALL QUIT('Insufficient work space in CC_FBTA.')
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Batch No.:',IBATCH
          WRITE (LUPRI,*) 'start at :',ISTART
          WRITE (LUPRI,*) '# transf.:',NTRAN
          WRITE (LUPRI,*) 'kend     :',KEND
          IDUM = 0 
          WRITE (LUPRI,*) 'work(0)  :',WORK(IDUM)
        END IF

        CALL CC_FBTA1(IFBTRAN, NFBTRAN, IOPTRES, LISTL, LISTR,
     &                FILFBTA, IFBDOTS, FBCONS, MXVEC,
     &                WORK(ISTART), WORK(IEND), 
     &                WORK(IADRE0), WORK(IADRE1),
     &                WORK(IADRPQ0), WORK(IADRPQ1), WORK(IADRPQMO),
     &                WORK(IADRPQA0),WORK(IADRPQA1),WORK(IADPQAMO),
     &                WORK(IADRBFA),WORK(IADRBFQA),
     &                WORK(IADRBZA),WORK(IADRBZQA),
     &                WORK(IT2DEL0), WORK(IT2DELB),
     &                WORK(IT2DELA), WORK(IT2DELBA),
     &                WORK(KCMOPQ), WORK(KCMOHQ),
     &                WORK(KLAMDPQ), WORK(KLAMDHQ),
     &                WORK(KLAMDPA), WORK(KLAMDHA),
     &                WORK(KLAMPQA), WORK(KLAMHQA),
     &                WORK(KONEHB),
     &                WORK(KDENSB),  WORK(KDPCKB), 
     &                WORK(KDENS0A), WORK(KDPCK0A), 
     &                WORK(KDENSBA), WORK(KDPCKBA), 
     &                WORK(KFOCKBCC),WORK(KFOCKACC),WORK(KFCKBACC), 
     &                WORK(KRAIM),WORK(KRQAIM),
     &                WORK(KFAIM),WORK(KFQAIM),
     &                WORK(KTHETA1),WORK(KTHETA2),
     &                WORK(KRHOBFA),WORK(KRHOBFQA),
     &                WORK(KXAINT),WORK(KYAINT),
     &                WORK(KBFZI),WORK(KBZADEN),WORK(KBZQADEN), 
     &                WORK(KEND), LEND )

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'returned from CC_FBTA1:'
          IDUM = 0
          WRITE (LUPRI,*) 'work(0)  :',WORK(IDUM)
        END IF

      END DO

      RETURN
      END

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_FBTA                              *
*---------------------------------------------------------------------*
*=====================================================================*
c /* deck cc_fbta1 */
*=====================================================================*
       SUBROUTINE CC_FBTA1( IFBTRAN, NFBTRAN, IOPTRES, LISTL, LISTR,
     &                FILFBTA, IFBDOTS, FBCONS, MXVEC,
     &                ISTART, IEND,IADRE0, IADRE1,
     &                IADRPQ0, IADRPQ1, IADRPQMO,
     &                IADRPQA0,IADRPQA1,IADPQAMO,
     &                IADRBFA,IADRBFQA,
     &                IADRBZA,IADRBZQA,
     &                IT2DEL0, IT2DELB, IT2DELA, IT2DELBA,
     &                KCMOPQ, KCMOHQ, KLAMDPQ, KLAMDHQ,
     &                KLAMDPA, KLAMDHA, KLAMPQA, KLAMHQA,
     &                KONEHB, KDENSB,  KDPCKB,
     &                KDENS0A, KDPCK0A, KDENSBA, KDPCKBA,
     &                KFOCKBCC,KFOCKACC,KFCKBACC,
     &                KRAIM,KRQAIM,
     &                KFAIM,KFQAIM,
     &                KTHETA1,KTHETA2,
     &                KRHOBFA,KRHOBFQA,
     &                KXAINT,KYAINT,
     &                KBFZI,KBZADEN,KBZQADEN,
     &                WORK, LWORK)
*---------------------------------------------------------------------*
*
*   Purpose: Calculation of the F^B T^A vector for a generalized
*            Hamiltonian B (i.e. 2-electron operator)
*            describing a (first-order) perturbation including the
*            contributions of derivative integrals, orbital relaxation
*            and reorthonormalization (orbital connection).
*
*    LISTL         --  Type of Left vector (Zeta)
*    LISTR         --  Type of Right vector (T^A)
*    IFBTRAN(1,*)  --  B operator indices in LBLOPR array 
*    IFBTRAN(2,*)  --  indices for left vectors for RES result vector
*    IFBTRAN(3,*)  --  indices for right vector for RES result vector
*    IFBTRAN(4,*)  --  indices for output RES vector on FILFBTA list
*    IFBTRAN(5,*)  --  flag for orbital relaxation
*    FILFBTA       --  list type of FBTA vectors or of the vectors that
*                      are dotted on FBTA vector (IOPTRES = 5)
*    IFBDOTS       --  indices for vectors to be dotted with
*    FBCONS        --  final dot product contributions (scalar) on return
*
*    NFBTRAN  -- number of requested transformations/vectors
*
*    IOPTRES = 0 :  all result vectors are written to a direct
*                   access file, FILFBTA is used as file name
*                   the start addresses of the vectors are
*                   returned in IFBTRAN(4,*)
*
*    IOPTRES = 1 :  the vectors are kept and returned in WORK
*                   if possible, start addresses returned in
*                   IFBTRAN(4,*). N.B.: if WORK is not large
*                   enough iopt is automatically reset to 0!!!
*                   (each result vector is passed back to caller
*                    routine as first element in WORK (debug purposes)
*                   NOT WORKING pga INTERFACE ROUTINE
*
*    IOPTRES = 3 : each result vector is written to its own file by
*                  a CALL to CC_WRRSP using FILFBTA as list type
*                  and IFBTRAN(4,*) as list index
*                  NOTE that IFBTRAN(4,*) is in this case input!
*                  NOT YET TESTED NOR USED
*
*    IOPTRES = 4 : each result vector is added to a vector on file by
*                  a CALL to CC_WARSP using FILFBTA as list type
*                  and IFBTRAN(4,*) as list index
*                  NOTE that IFBTRAN(4,*) is in this case input!
*                  NOT YET TESTED NOR USED
*    
*    IOPTRES = 5 : each result vector is dotted with (an array of) vectors 
*                  in IFBDOTS, the type of the array of vectors to dot on is
*                  given by FILFBTA and indices are taken from the array
*                  IFBDOTS. Scalar results are returned in the FBCONS matrix
*
*    IT2DEL0,IT2DELB  -- integer scratch arrays of dimensions
*                        MXCORB and MXCORB x NFBTRAN  (integrals)
*    IT2DELA,IT2DELBA -- integer scratch arrays of dimensions
*                        MXCORB x NFBTRAN       (integrals)
*    operator labels:
*
*           HAM0     : unperturbed Hamiltonian (1-el & 2-el part)
*                      (for test purposes)
*           1DHAMxxx : geometrical first derivatives (1-el & 2-el part)
*           dh/dB xx : London first derivatives (1-el & 2-el part)
*
*           ELSE the label is interpreted as a one-electron operator
*           and the relaxation/reorthonormalization terms are skipped
*
*    Written by Sonia Coriani, 1999
*    based on cc_xieta, cc_fmat
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <dummy.h>
#include <aovec.h>
#include <maxorb.h>
#include <mxcent.h>
#include <blocks.h>
#include <maxaqn.h>
#include <symmet.h>
#include <nuclei.h>
#include <ccorb.h>
#include <ccfield.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <iratdef.h>
#include <distcl.h>
#include <eribuf.h>
#include <ccisao.h>
#include <ccroper.h>
#include <ccfro.h>
#include <ccr1rsp.h>
#include <chrxyz.h>
#include <r12int.h>
#include <cbieri.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)       
      INTEGER KDUM
      PARAMETER (KDUM = +99 999 999)
      INTEGER IDUM, IDUM1
      PARAMETER (IDUM = 0)         
* files
      INTEGER LUSIRG                                    !CMO(0) coefficients
      INTEGER LUFOCK                                    !FockAO(0) matrix
      INTEGER LU0IAJB,LU0IJCB,LU0CJIB                   
      INTEGER LU1IAJB,LU1IJCB,LU1CJIB                  
      INTEGER LU0IJKL,LU1IJKL
      INTEGER LUBFDA, LUBFDQA           !effec. BF(A) and BF(QA) densities
      INTEGER LUBFDE0, LUBFDE1          !the d and d-bar eff densities (BZ)
      INTEGER LUPQMO,LUPQ0, LUPQ1       !P&Q,P&Q^del,P&Q^del-bar,
      INTEGER LUPQAMO, LUPQA0, LUPQA1   !PA&QA,PA&QA^del,PA&QA^del-bar
      INTEGER LUBFZA, LUBFZQA           !the BZA (total) intermediates
      INTEGER LUBFA, LUBFQA             !the BF(A) and BF(QA) intermediates 
      INTEGER LUFBTA                    !unit for output (iopres =0,1)
C      PARAMETER (LU0IAJB = 71, LU1IAJB = 72)
C      PARAMETER (LU0IJCB = 73, LU1IJCB = 74)
C      PARAMETER (LU0CJIB = 75, LU1CJIB = 76)
C      PARAMETER (LUBFDA  = 77, LUBFDQA = 78, LUBFDE0 = 79, LUBFDE1= 80)
C      PARAMETER (LUPQMO  = 81, LUPQ0   = 82, LUPQ1   = 83)
C      PARAMETER (LUPQAMO = 84, LUPQA0  = 85, LUPQA1  = 86)
C      PARAMETER (LUBFZA  = 87, LUBFZQA = 88)
C      PARAMETER (LUBFA   = 77, LUBFQA  = 78)    !genbrug unit-number
C      PARAMETER (LU0IJKL = 90, LU1IJKL = 91)
C      PARAMETER (LUFBTA  = 89)

      CHARACTER*8 FN0IAJB,  FN1IAJB
      CHARACTER*8 FN0IJCB,  FN0CJIB, FN1IJCB, FN1CJIB
      CHARACTER*8 FN0IJKL,  FN1IJKL
      CHARACTER*8 FNBFDA,   FNBFDQA, FNBFDE0, FNBFDE1
      CHARACTER*8 FILPQMO,  FILPQ0,  FILPQ1
      CHARACTER*8 FILPQAMO, FILPQA0, FILPQA1
      CHARACTER*8 FILBFZA,  FILBFZQA
      CHARACTER*8 FILBFA,   FILBFQA
      PARAMETER (FN0IAJB = 'CCXIAJB0', FN1IAJB = 'CCXIAJB1')
      PARAMETER (FN0IJCB = 'CCXIJCB0', FN1IJCB = 'CCXIJCB1')
      PARAMETER (FN0CJIB = 'CCXCJIB0', FN1CJIB = 'CCXCJIB1')
      PARAMETER (FN0IJKL = 'CCXIJKL0', FN1IJKL = 'CCXIJKL1')
      PARAMETER (FNBFDA  = 'CCBFDEN1', FNBFDQA  = 'CCBFDEN5')     
      PARAMETER (FNBFDE0 = 'CCBFDEN2', FNBFDE1  = 'CCBFDEN6')      
      PARAMETER (FILPQMO = 'CCPQIMMO', FILPQ0   = 'CCPQINT0')    
      PARAMETER (FILPQ1  = 'CCPQINT1')
      PARAMETER (FILPQAMO= 'CCPQAIMO')
      PARAMETER (FILPQA0 = 'CCPQAIN0', FILPQA1  = 'CCPQAIN1')
      PARAMETER (FILBFZA = 'CCBFZAIN', FILBFZQA = 'CCBFZQAI')       
      PARAMETER (FILBFA  = 'CCBFAINT', FILBFQA  = 'CCBFQAIN')

#if defined (SYS_CRAY)
      REAL DDOT, DNRM2
      REAL ONE, TWO, THREE, HALF, ZERO, FAC
#else
      DOUBLE PRECISION DDOT, DNRM2
      DOUBLE PRECISION ONE,TWO,THREE,HALF,ZERO,FAC
#endif
      PARAMETER (ONE = 1.0d0, TWO = 2.0d0, THREE = 3.0d0, HALF = 0.5d0)
      PARAMETER (ZERO = 0.0d0)

* argument-list

      CHARACTER*10 MODEL
      CHARACTER CDUM
      CHARACTER*(*) LISTL, LISTR
      CHARACTER*8 FILFBTA, LABELH, LAB1
      INTEGER LWORK, NFBTRAN, IOPTRES, MAXSIM, MXVEC

* pointers in memory

      INTEGER IFBTRAN(5,NFBTRAN)
      INTEGER IFBDOTS(MXVEC,NFBTRAN)             
      INTEGER ISTART(NFBTRAN), IEND(NFBTRAN)
C
      INTEGER IT2DEL0(MXCORB)                 !the 3-MO integrals g(0)
      INTEGER IT2DELB(MXCORB,NFBTRAN)         !the 3-MO integrals g(1)
      INTEGER IT2DELA(MXCORB,NFBTRAN)         !the 3-MO integrals g(0A) 
      INTEGER IT2DELBA(MXCORB,NFBTRAN)        !the 3-MO integrals g(1A) 

      INTEGER KRAIM(NFBTRAN), KFAIM(NFBTRAN)      !calR,calF
      INTEGER KRQAIM(NFBTRAN), KFQAIM(NFBTRAN)    !calR-bar,calF-bar

      INTEGER KTHETA1(NFBTRAN), KTHETA2(NFBTRAN)      !result
      !the eff. densities BZA and BZQA (MGD)
      INTEGER KBZADEN(NFBTRAN), KBZQADEN(NFBTRAN)     
      INTEGER KFOCKBCC(NFBTRAN)                       !Fock-bar (AO) 
      INTEGER KONEHB(NFBTRAN)                         !h-bar (derivative?)
      INTEGER KFOCKACC(NFBTRAN)                       !calFock(0) matrix
      INTEGER KFCKBACC(NFBTRAN)                       !calFock(1) matrix
      INTEGER KBFZI(NFBTRAN)

      INTEGER KLAMDPQ(NFBTRAN),KLAMDHQ(NFBTRAN)       !LambdaQ-p,LambdaQ-h
      INTEGER KLAMDPA(NFBTRAN),KLAMDHA(NFBTRAN)       !LambdaA-p,LambdaA-h
      INTEGER KLAMPQA(NFBTRAN),KLAMHQA(NFBTRAN)       !LambdaQA-p,LambdaQA-h
      !Derivative MO coefficients (C*Q)
      INTEGER KCMOPQ(NFBTRAN), KCMOHQ(NFBTRAN)        

      INTEGER KDENSB(NFBTRAN), KDPCKB(NFBTRAN)      
      INTEGER KDENS0A(NFBTRAN), KDPCK0A(NFBTRAN)   
      INTEGER KDENSBA(NFBTRAN), KDPCKBA(NFBTRAN)  
      INTEGER KRHOBFA(NFBTRAN), KRHOBFQA(NFBTRAN)

      INTEGER KXAINT(NFBTRAN), KYAINT(NFBTRAN)   

cs pointers for quantities on file (for each del or f, and in some cases ITRAN)
cs integer arrays addresses

      INTEGER IADRE0(MXCORB,NFBTRAN),  IADRE1(MXCORB,NFBTRAN)   
      INTEGER IADRBFA(MXCORB,NFBTRAN), IADRBFQA(MXCORB,NFBTRAN)
      INTEGER IADRBZA(MXCORB,NFBTRAN), IADRBZQA(MXCORB,NFBTRAN)

      INTEGER IADRPQ0(MXCORB,NFBTRAN), IADRPQ1(MXCORB,NFBTRAN) 
      INTEGER IADRPQMO(MXCORB,NFBTRAN)                        

      INTEGER IADRPQA0(MXCORB,NFBTRAN), IADRPQA1(MXCORB,NFBTRAN)
      INTEGER IADPQAMO(MXCORB,NFBTRAN)                          

      INTEGER INDEXA(MXCORB), IGAM(MXCORB)                     

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL FBCONS(MXVEC,NFBTRAN)
      REAL FREQ, XNORM
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION FBCONS(MXVEC,NFBTRAN)
      DOUBLE PRECISION FREQ, XNORM
#endif
      CHARACTER*10 MODELW
      LOGICAL LDERINT(8,3), LNEWTA, LNEWOP, LNEWZE, SYM1ONLY
      LOGICAL LTWOEL, LRELAX, LLTWOEL, LLRELAX, NEWTYPE, DIRINT, LZERO
      LOGICAL LWRDSK, SQRINT, LSQRUP, LX1ISQ, LLAO
      INTEGER MT2BGD, MT2ORT, MT2BCD, MT2AM, MT2SQ, M2BST, MEMAT1,MT1AO
      INTEGER MDISAO, MGAMMA, MSCRATCH0, MSCRATCH1, MWORK, NBATCH,MT1AM
      INTEGER ISYM, IATOM, ITRAN, IOPERB, IBATCH, IOPT, IFIELD, ILLL
      INTEGER ISYHOP, ISYRES, ISYMD1, IVEC, ISYIJK, ISYMD, IOFF
      INTEGER IDR_OLD, IOP_OLD, IDL_OLD, ID_RES
      INTEGER NTOSYM, NTOT, NUMDIS, IDEL2, IDEL
      INTEGER ISY0DIS, ISY1DIS, ISYDEL, ISYMG, NUMG,ICOOR, ICORSY
      INTEGER NGDER, NBDER, NBUFMX, NFILES, KNRECS
      INTEGER KCCFB1,KINDXB, KEND,LWRK, KODCL1, KODCL2, KODBC1, KODBC2
      INTEGER KRDBC1,KRDBC2, KODPP1,KODPP2, KRDPP1,KRDPP2,KFREE,LFREE
      INTEGER KRECNR, ITYPE, NUMD, MXCOMP, JGAM
      INTEGER IFILE, KEND1SV, LWRK1SV, ITST, KT2AMP0, LEN0, LEN1
      INTEGER IDLSTR,IDLSTL,ISYMTA,ISYHTA,ISYZTA, IREAL
cs
      INTEGER IADRI0, IADRIB, IADRI0A, IADRIAB                
      INTEGER IADRBFE0, IADRBFE1, IADBFA, IADBFQA, IADR, IADRCI 
      INTEGER IADRPQ, IADRPQI0, IADRPQI1, IADR0, IADRB
      INTEGER IADRPQA, IADRPQAI0, IADRPQAI1
      INTEGER IADRBZ1, IADRBZ2
      INTEGER IADR1, IADR2, IADR_FBTA
      INTEGER NBSRHF(8), IBSRHF(8,8), ICOUNT, ISYMAK, ISYBET
cs
      INTEGER KONEH0, KFOCK0, IOPTW, ISYMJ
      INTEGER IRELAX, LUDIM, LUCIM, LEN
      INTEGER KGAMMA, KEMAT1, KEMAT2, ISY0IAJ, ISY1IAJ, KXIAJB, KDPK0
      INTEGER K1XINT, K0XINT, KD1SRHF, KD0SRHF, KEND3, LWRK3
      INTEGER KZETA1, KZETA2, KBFZ0, KEND2, LWRK2, ISYCTR
      INTEGER KEND0, LWRK0, KEND1, LWRK1, KF0IM, KG0IM, KR0IM
      INTEGER KEND4, LWRK4, LGIM
      INTEGER ICORE, KLAMH0, KLAMP0, KDENS0, KT1AMP0
      INTEGER KKAPPA, KRMAT, IDXR, IERR, IORDER
      INTEGER KXINT, KYINT, KMINT
      INTEGER KTA1AM, KTA2AM, KMAINT, KSCR
      INTEGER LENBFA, LENBFQA 
      INTEGER KGAMMAT, KRES1, KRES2, LENTEST
      INTEGER KBF0RHF, KBF1RHF, N3ODMX, LX3MO, KX3MO, LENIJK, KXIJK0
      INTEGER JCOOR, JSCOOR, ISCOOR, IOPTYP,KOFF
cnewforcc2
      INTEGER KZDPK0, KZDEN0, KZDPKB, KZDENB, KCHI, KCHIQ
cstest
      INTEGER LENALL
    

* external functions:
      INTEGER ILSTSYM

*---------------------------------------------------------------------*
* begin: check wavefunction model and flush print unit
*---------------------------------------------------------------------*
      IF (CCSDT) THEN
        WRITE (LUPRI,'(/1x,a)') 
     *      'CC_FBTA not yet implemented for triples ...'
        CALL QUIT('Triples not implemented in CC_FBTA.')
      END IF

      IF ( .NOT. (CCS .OR. CC2 .OR. CCSD) ) THEN
        WRITE (LUPRI,'(/1x,a)') 'CC_FBTA called for a Coupled Cluster '
     &          //'method not implemented in CC_FBTA...'
        CALL QUIT('Unknown CC method in CC_FBTA.')
      END IF
      
      IF (CC2) THEN
        WRITE (LUPRI,'(/1x,a)') 'CC_FBTA not YET implemented for CC2 '
        CALL QUIT('CC2 method not yet available in CC_FBTA.')
      END IF

      CALL FLSHFO(LUPRI)

      IF (LOCDBG) THEN
        ITST = 0
        WRITE (LUPRI,'(/1x,a,i15)')  'Work space in CC_FBTA:',LWORK
        WRITE (LUPRI,*) 'FILFBTA  = ',FILFBTA(1:3)
        WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
        CALL FLSHFO(LUPRI)
      END IF
*---------------------------------------------------------------------*
* check return option for the result vectors, initialize output:
*---------------------------------------------------------------------*
      IF (IOPTRES .EQ. 0 .OR. IOPTRES .EQ. 1) THEN
         LUFBTA = -1
         CALL WOPEN2(LUFBTA, FILFBTA, 64, 0)
         IADR_FBTA = 1
      ELSE IF (IOPTRES .EQ. 3 .OR. IOPTRES .EQ. 4) THEN
         CONTINUE
      ELSE IF (IOPTRES .EQ. 5) THEN
         IF (MXVEC*NFBTRAN.GT.0) THEN
            CALL DZERO(FBCONS,MXVEC*NFBTRAN)
         ELSE
            WRITE (LUPRI,*) 
     *         'WARNING: CC_FBTA called, but nothing to do!'
            RETURN
         END IF
      ELSE
         CALL QUIT('Illegal value of IOPTRES in CC_FBTA.')
      END IF

      IF (CCS) THEN
         MODELW = 'CCS       '
         IOPTW  = 1
      ELSE IF (CC2) THEN
         MODELW = 'CC2       '
         IOPTW  = 3
      ELSE IF (CCSD) THEN
         MODELW = 'CCSD      '
         IOPTW  = 3
      ELSE
         CALL QUIT('Unknown coupled cluster model in CC_FBTA.')
      END IF

*---------------------------------------------------------------------*
* precalculate symmetry array for BSRHF:
*---------------------------------------------------------------------*
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBSRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBSRHF(ISYM) = ICOUNT
      END DO
*---------------------------------------------------------------------*
* precalculate max dimensions for IJK,del distributions 
*---------------------------------------------------------------------*
      N3ODMX = 0
      DO ISYM  = 1, NSYM
        N3ODMX = MAX(N3ODMX,N3ODEL(ISYM))
      END DO

      LX3MO = 0
      DO ISYIJK = 1, NSYM
         LX3MO  = MAX(LX3MO,NMAIJK(ISYIJK))
      END DO

*---------------------------------------------------------------------*
* open 2 files for the different integrals:
*---------------------------------------------------------------------*
      LU0IAJB = -1
      LU1IAJB = -1
      CALL WOPEN2(LU0IAJB, FN0IAJB, 64, 0)
      CALL WOPEN2(LU1IAJB, FN1IAJB, 64, 0)
*---------------------------------------------------------------------*
* open additional 4  for the special (combined) integrals (C, D terms):
*---------------------------------------------------------------------*
      LU0IJCB = -1
      LU0CJIB = -1
      LU1IJCB = -1
      LU1CJIB = -1
      CALL WOPEN2(LU0IJCB, FN0IJCB, 64, 0)
      CALL WOPEN2(LU0CJIB, FN0CJIB, 64, 0)
      CALL WOPEN2(LU1IJCB, FN1IJCB, 64, 0)
      CALL WOPEN2(LU1CJIB, FN1CJIB, 64, 0)
*---------------------------------------------------------------------*
* open 2 files for the 3 occupied integrals:
*---------------------------------------------------------------------*
      LU0IJKL = -1
      LU1IJKL = -1
      CALL WOPEN2(LU0IJKL, FN0IJKL, 64, 0)
      CALL WOPEN2(LU1IJKL, FN1IJKL, 64, 0)
*---------------------------------------------------------------------*
* open files for several effective densities: 
*---------------------------------------------------------------------*
      LUBFDA  = -1
      LUBFDQA = -1
      LUBFDE0 = -1
      LUBFDE1 = -1
      CALL WOPEN2(LUBFDA,  FNBFDA,  64, 0)
      CALL WOPEN2(LUBFDQA, FNBFDQA, 64, 0)
      CALL WOPEN2(LUBFDE0, FNBFDE0, 64, 0)
      CALL WOPEN2(LUBFDE1, FNBFDE1, 64, 0)
*---------------------------------------------------------------------*
* open files for BZ(A) intermediates 
*---------------------------------------------------------------------*
      LUBFZA  = -1
      LUBFZQA = -1
      CALL WOPEN2(LUBFZA,  FILBFZA,  64, 0)
      CALL WOPEN2(LUBFZQA, FILBFZQA, 64, 0)
*---------------------------------------------------------------------*
* open files for P and Q intermediates: 
*---------------------------------------------------------------------*
      LUPQMO = -1
      LUPQ0  = -1
      LUPQ1  = -1
      CALL WOPEN2(LUPQMO, FILPQMO, 64, 0)
      CALL WOPEN2(LUPQ0,  FILPQ0,  64, 0)
      CALL WOPEN2(LUPQ1,  FILPQ1,  64, 0)
*---------------------------------------------------------------------*
* open files for calP and calQ intermediates (PQA): 
*---------------------------------------------------------------------*
      LUPQAMO = -1
      LUPQA0  = -1
      LUPQA1  = -1
      CALL WOPEN2(LUPQAMO, FILPQAMO, 64, 0)
      CALL WOPEN2(LUPQA0,  FILPQA0,  64, 0)
      CALL WOPEN2(LUPQA1,  FILPQA1,  64, 0)

*=====================================================================*
* prepare batching:
*=====================================================================*
*     estimate memory requirements: THIS NEEDS TO BE CHANGED
      MT2BGD = 0
      MT2BCD = 0
      MT2ORT = 0
      MT2AM  = 0
      MT2SQ  = 0
      MDISAO = 0
      MEMAT1 = 0
      M2BST  = 0
      MGAMMA = 0
      MT1AO  = 0
      MT1AM  = 0
      DO ISYM = 1, NSYM
        MT2BCD = MAX(MT2BCD,NT2BCD(ISYM))
        MT2BGD = MAX(MT2BGD,NT2BGD(ISYM))
        MT2ORT = MAX(MT2ORT,NT2ORT(ISYM))
        MT2AM  = MAX(MT2AM ,NT2AM(ISYM))
        MT2SQ  = MAX(MT2SQ ,NT2SQ(ISYM))
        MDISAO = MAX(MDISAO,NDISAO(ISYM))
        MEMAT1 = MAX(MEMAT1,NEMAT1(ISYM))
        M2BST  = MAX(M2BST ,N2BST(ISYM) )
        MGAMMA = MAX(MGAMMA,NGAMMA(ISYM))
        MT1AO  = MAX(MT1AO, NT1AO(ISYM))
        MT1AM  = MAX(MT1AM, NT1AM(ISYM))
      END DO

* fixed requirement, independent of batch length:
*     maximum of MT2SQ and MGAMMA (for A, C and D terms)
*     a NT2BGD array (CC_BF) and 6 NT2BCD arrays (transformations)
*     + NT2SQ (T2AMP0, and some other intermediates) + NEMAT1(1) (R0IM)
*     + 2 (**|*del) integral arrays and some extra reserve

      IF (CCS) THEN
        MSCRATCH0 = 2*MDISAO + 10*N2BASX
      ELSE IF (CC2 .OR. CCSD) THEN
        MSCRATCH0 = MAX(MT2SQ,MGAMMA,MT2BGD,6*MT2BCD) + 2*MDISAO + 
     &        MAX(MT2SQ,2*MT2ORT) + NEMAT1(ISYM0) + 10*N2BASX
      ELSE
        CALL QUIT('Unknown CC model in CC_FBTA.')
      END IF

* requirements per simultaneous transformation:
*     BF and R intermediate and Fock matrix and one-electron integrals

      IF (CCS) THEN
        MSCRATCH1 = 2*M2BST
      ELSE IF (CC2 .OR. CCSD) THEN
        MSCRATCH1 = MAX(2*MT2ORT,MT2AM) + MEMAT1 + 2*M2BST
      ELSE
        CALL QUIT('Unknown CC model in CC_FBTA1.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_FBTA> scratch space estimates:'
        WRITE (LUPRI,*) 'CC_FBTA> fixed (MSCRATCH0)   : ',MSCRATCH0
        WRITE (LUPRI,*) 'CC_FBTA> per simult. transf. : ',MSCRATCH1
        WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
      END IF

* prepare batching, make sure that for derivatives all transf. 
* in a batch belong to coordinates of the same atom,
* do not mix derivatives with other perturbations
* -> IOPTYP = 0 : 1el perturbations
*    IOPTYP = 1 : first geom. derivatives
*    IOPTYP = 2 : first magn. derivatives (London integrals)

      NBATCH = 1
      MWORK  = 0
      IATOM  = 0
      IOPTYP = 0
      ISTART(NBATCH) = 1
      DO ITRAN = 1, NFBTRAN

        IOPERB = IFBTRAN(1,ITRAN)
        LABELH = LBLOPR(IOPERB)

        IF (LABELH(1:5).EQ.'1DHAM') THEN

           NEWTYPE = (IATOPR(IOPERB).NE.IATOM .OR. IOPTYP.NE.1)
           IATOM   = IATOPR(IOPERB)
           IOPTYP  = 1
           IF ( (NEWTYPE .AND. ITRAN.GT.1) .OR. 
     &          ((MWORK+MSCRATCH1+MSCRATCH0).GT.LWORK) ) THEN
              IEND(NBATCH)   = ITRAN - 1
              NBATCH         = NBATCH + 1
              ISTART(NBATCH) = ITRAN
              MWORK = 0
           END IF

        ELSE IF (LABELH(1:5).EQ.'dh/dB') THEN

           NEWTYPE = (IOPTYP.NE.2)
           IATOM   = 0
           IOPTYP  = 2
           IF ( (NEWTYPE .AND. ITRAN.GT.1) .OR.
     &          ((MWORK+MSCRATCH1+MSCRATCH0).GT.LWORK) ) THEN
              IEND(NBATCH)   = ITRAN - 1
              NBATCH         = NBATCH + 1
              ISTART(NBATCH) = ITRAN
              MWORK = 0
           END IF

        ELSE

           NEWTYPE = (IOPTYP.NE.0)
           IATOM   = 0
           IOPTYP  = 0
           IF ( (NEWTYPE .AND. ITRAN.GT.1) .OR.
     &          ((MWORK+MSCRATCH1+MSCRATCH0).GT.LWORK) ) THEN
              IEND(NBATCH)   = ITRAN - 1
              NBATCH         = NBATCH + 1
              ISTART(NBATCH) = ITRAN
              MWORK = 0
           END IF

        END IF

        MWORK = MWORK + MSCRATCH1
        IF ((MWORK+MSCRATCH0) .GT. LWORK) THEN
           WRITE (LUPRI,*) 'Insufficient work space in CC_FBTA.'
           WRITE (LUPRI,*) 'Available:',LWORK,' words.'
           WRITE (LUPRI,*) 'Need at least:',MWORK+MSCRATCH0,' words.'
           CALL FLSHFO(LUPRI)
         CALL QUIT('Insufficient work space in CC_FBTA.')
        END IF

      END DO
      IEND(NBATCH) = NFBTRAN

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_FBTA> batching statistics:'
        WRITE (LUPRI,*) 'CC_FBTA> NBATCH : ',NBATCH
        WRITE (LUPRI,*) 
     *       'CC_FBTA> ISTART : ',(ISTART(IBATCH),IBATCH=1,NBATCH)
        WRITE (LUPRI,*) 
     *       'CC_FBTA> IEND   : ',(IEND(IBATCH),  IBATCH=1,NBATCH)
        WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
      END IF

*     -------------------------------------------------------
*     set start address for 0.th-order intermediates on file:
*     (not depending on ITRAN)
*     -------------------------------------------------------
      IADRI0   = 1                             !g(0)_iaj del

*     --------------------------------------------------------------
*     allocate work space for 0.th-order Lambda matrices, MO vector,
*     densities, Fock matrices, etc., which are kept in memory:
*     --------------------------------------------------------------
      
      KDENS0  = 1 
      KDPK0   = KDENS0  + N2BST(ISYM0)   
      KT1AMP0 = KDPK0   + NNBST(ISYM0)  
      KLAMP0  = KT1AMP0 + NT1AMX                     
      KLAMH0  = KLAMP0  + NLAMDT                    
      KEND0   = KLAMH0  + NLAMDT        
      LWRK0   = LWORK   - KEND0

      IF (LWRK0 .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_FBTA. (0)')
      END IF

*     ----------------------------------------------------------------
*     initialize 0.th-order Lambda matrices, densities, etc.:
*     ----------------------------------------------------------------
      IOPT = 1
      CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))        

      CALL LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),                  
     &            WORK(KEND0),LWRK0)

c take care here, concerning frozen core? See fmat

      !CC Fock density 
      ICORE = 1
      CALL CC_AODENS(WORK(KLAMP0),WORK(KLAMH0),WORK(KDENS0),        
     &               ISYM0,ICORE,WORK(KEND0),LWRK0)

      !pack CC Fock density
      CALL CC_DNSPK(WORK(KDENS0),WORK(KDPK0),ISYM0)                

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c I removed all CM0 stuff CC2 by now (fock densities with CMO's etc
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

*     ----------------------------------------------------------------
*     get the zeroth-order one-electron Hamiltonian and Fock matrices:
*     (for CCS we have to calculate the Fock matrix and do only
*      initialize it here with the one-electron part)
*     ----------------------------------------------------------------
      KONEH0 = KEND0                              !h^0_{alpha,beta}
      KFOCK0 = KONEH0 + N2BST(ISYM0)              !AO Fock matrix F^0_{al,be}
      KEND0  = KFOCK0 + N2BST(ISYM0)

      LWRK0 = LWORK - KEND0
      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_FBTA. (FOCK0)')
      END IF

*     read zeroth-order one-electron hamiltonian from file
*     and add finite fields contributions (if any)
*     
      CALL CCRHS_ONEAO(WORK(KONEH0),WORK(KEND0),LWRK0)
      DO IFIELD = 1, NFIELD
        CALL CC_ONEP(WORK(KONEH0),WORK(KEND0),LWRK0,
     &               EFIELD(IFIELD),ISYM0,LFIELD(IFIELD))
      END DO                      
      IF (CCS) THEN
        CALL DCOPY(N2BST(ISYM0),WORK(KONEH0),1,WORK(KFOCK0),1)
      ELSE
*       read zeroth-order AO Fock matrix from file:
        LUFOCK = -1
        CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',
     *              IDUMMY,.FALSE.)
        REWIND(LUFOCK)
        READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
        CALL GPCLOSE(LUFOCK,'KEEP')
      END IF 

*=====================================================================*
* loop over the batches of transformations:
*=====================================================================*
      DO IBATCH = 1, NBATCH

*       ---------------------------------------------
*       set start addresses for intermediate on files (counters?):
*       ---------------------------------------------
        IADRIB   = 1                        !g(1)_iaj del
        IADRI0A  = 1                        !g~(0A)_ijc del, _cji del
        IADRIAB  = 1                        !g~(1A)_ijc del, _cji del
*
        IADRPQ   = 1                        !PQ MO
        IADRPQI0 = 1                        !PQ delta
        IADRPQI1 = 1                        !PQ-bar delta
*
        IADRBFE0 = 1                        !effect. d
        IADRBFE1 = 1                        !d-bar
*
        IADBFA  = 1                         !Delta(BFA) 
        IADBFQA = 1                         !Delta-bar(BFQA)
*
        IADRPQA   = 1                       !calPQ                 
        IADRPQAI0 = 1                       !calPQ_delta
        IADRPQAI1 = 1                       !calPQ-bar delta
*
        IADRBZ1 = 1                         !BZA interm
        IADRBZ2 = 1                         !BZQA interm
*
*---------------------------------------------------------------------*
* loop over transformations allocate work space and 
* initialize calR, calG, F, BF & BFZ intermediates: 
*---------------------------------------------------------------------*
        LLRELAX = .FALSE.
        LLTWOEL = .FALSE.
*
        KEND1 = KEND0
*
        IOP_OLD  = -1
        IDR_OLD  = -1
        IDL_OLD  = -1
*
        DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
*
          IOPERB  = IFBTRAN(1,ITRAN)   ! operator index
          IDLSTL  = IFBTRAN(2,ITRAN)   ! ZETA vector index
          IDLSTR  = IFBTRAN(3,ITRAN)   ! T^A  vector index
          IRELAX  = IFBTRAN(5,ITRAN)   ! flag for orbital relaxation
*
          LABELH = LBLOPR(IOPERB)      ! operator label
          ISYHOP = ISYOPR(IOPERB)      ! symmetry of the B operator
          IREAL  = ISYMAT(IOPERB)      ! flag for real/imag. operators
          LTWOEL = LPDBSOP(IOPERB)     ! two-electron contribution
          LRELAX = LTWOEL .OR. (IRELAX.GE.1)
*
          LLTWOEL = (LLTWOEL .OR. LTWOEL)  ! any two-el. perturbation
          LLRELAX = (LLRELAX .OR. LRELAX)  ! any relaxed perturbation   (updated/Checked at each iter)

          LZERO   = (ITRAN .EQ. 1)         ! test for zero-order interm.

          LNEWTA = (IDLSTR .NE. IDR_OLD)
          LNEWZE = (IDLSTL .NE. IDL_OLD)
          LNEWOP = (IOPERB .NE. IOP_OLD)
*
          ISYCTR  = ILSTSYM(LISTL,IDLSTL) ! symmetry of ZETA vector
          ISYMTA  = ILSTSYM(LISTR,IDLSTR) ! symmetry of T^A vector
          ISYZTA  = MULD2H(ISYCTR,ISYMTA) ! symmetry of Zeta x T^A
          ISYHTA  = MULD2H(ISYMTA,ISYHOP) ! symmetry of Oper x T^A
          ISYRES  = MULD2H(ISYHOP,ISYZTA) ! symmetry result vect.
*
          IF (LOCDBG) THEN
            WRITE (LUPRI,*) 'IOPERB: ', IOPERB, ' IOP_OLD: ', IOP_OLD
            WRITE (LUPRI,*) 'IDLSTR: ', IDLSTR, ' IDR_OLD: ', IDR_OLD
          END IF
 
C         --------------------------------------------------------------
C         Some intermediates (allocations) for RES: 
C         Always a new allocation if IOPERB changes. 
C         If IOPERB = old, check if T^A has changed.
C         There are some redundancies but it's the easiest way right now.
C         --------------------------------------------------------------

          IF (LNEWOP) THEN                      !if new IOPER (and RELAX?)
*
*  only operator dependent
*
             KLAMDPQ(ITRAN)  = KEND1   
             KLAMDHQ(ITRAN)  = KLAMDPQ(ITRAN)  + NGLMDT(ISYHOP) 
             KDENSB(ITRAN)   = KLAMDHQ(ITRAN)  + NGLMDT(ISYHOP)
             KDPCKB(ITRAN)   = KDENSB(ITRAN)   + N2BST(ISYHOP)
             KONEHB(ITRAN)   = KDPCKB(ITRAN)   + NNBST(ISYHOP) 
             KFOCKBCC(ITRAN) = KONEHB(ITRAN)   + N2BST(ISYHOP)
             KCMOPQ(ITRAN)   = KFOCKBCC(ITRAN) + 
     &                         MAX(N2BST(ISYHOP),N2BST(ISYHTA))  !FOCKB reused in CCFBTAF
             KCMOHQ(ITRAN)   = KCMOPQ(ITRAN)   + NGLMDT(ISYHOP)
             KEND1           = KCMOHQ(ITRAN)   + NGLMDT(ISYHOP)
          
          ELSE
            
             KLAMDPQ(ITRAN)  = KLAMDPQ(ITRAN - 1)
             KLAMDHQ(ITRAN)  = KLAMDHQ(ITRAN - 1)
             KDENSB(ITRAN)   = KDENSB(ITRAN - 1)
             KDPCKB(ITRAN)   = KDPCKB(ITRAN - 1)
             KONEHB(ITRAN)   = KONEHB(ITRAN - 1)
             KFOCKBCC(ITRAN) = KFOCKBCC(ITRAN - 1)
             KCMOPQ(ITRAN)   = KCMOPQ(ITRAN - 1)
             KCMOHQ(ITRAN)   = KCMOHQ(ITRAN - 1)

          END IF
          
          IF (LNEWTA) THEN
*
*  only T^A dependent 
*
             KLAMDPA(ITRAN)  = KEND1
             KLAMDHA(ITRAN)  = KLAMDPA(ITRAN)  + NGLMDT(ISYMTA)      
             KDENS0A(ITRAN)  = KLAMDHA(ITRAN)  + NGLMDT(ISYMTA)
             KDPCK0A(ITRAN)  = KDENS0A(ITRAN)  + N2BST(ISYMTA)
             KFOCKACC(ITRAN) = KDPCK0A(ITRAN)  + NNBST(ISYMTA)
             KEND1           = KFOCKACC(ITRAN) + N2BST(ISYMTA)
             IF (.NOT.(CCS)) THEN
                KRAIM(ITRAN) = KEND1
                KEND1        = KRAIM(ITRAN) + NEMAT1(ISYMTA)    
             END IF
          ELSE
             KLAMDPA(ITRAN)  = KLAMDPA(ITRAN - 1)
             KLAMDHA(ITRAN)  = KLAMDHA(ITRAN - 1)
             KDENS0A(ITRAN)  = KDENS0A(ITRAN - 1)
             KDPCK0A(ITRAN)  = KDPCK0A(ITRAN - 1)
             KFOCKACC(ITRAN) = KFOCKACC(ITRAN -1)
             IF (.NOT.(CCS)) THEN
                KRAIM(ITRAN) = KRAIM(ITRAN - 1)
             END IF
          END IF
*          
*  both operator and T^A dependent
*
          IF ((LNEWTA).OR.(LNEWOP)) THEN

             KLAMHQA(ITRAN)  = KEND1
             KLAMPQA(ITRAN)  = KLAMHQA(ITRAN)  + NGLMDT(ISYHTA)
             KDENSBA(ITRAN)  = KLAMPQA(ITRAN)  + NGLMDT(ISYHTA)
             KDPCKBA(ITRAN)  = KDENSBA(ITRAN)  + N2BST(ISYHTA)
             KFCKBACC(ITRAN) = KDPCKBA(ITRAN)  + NNBST(ISYHTA)
             KEND1           = KFCKBACC(ITRAN) + N2BST(ISYHTA)
             IF (.NOT. CCS) THEN
                KRQAIM(ITRAN) = KEND1
                KEND1         = KRQAIM(ITRAN) + NEMAT1(ISYHTA)
             END IF

          ELSE
 
             KLAMHQA(ITRAN)  = KLAMHQA(ITRAN - 1) 
             KLAMPQA(ITRAN)  = KLAMPQA(ITRAN - 1)
             KDENSBA(ITRAN)  = KDENSBA(ITRAN - 1)
             KDPCKBA(ITRAN)  = KDPCKBA(ITRAN - 1)
             KFCKBACC(ITRAN) = KFCKBACC(ITRAN - 1)
             IF (.NOT. CCS) THEN
                KRQAIM(ITRAN) = KRQAIM(ITRAN - 1)
             END IF
 
          END IF
*
* T^A and Zeta dependent intermediates
*
          IF (.NOT.(CCS)) THEN

             IF ((LNEWTA).OR.(LNEWZE)) THEN

                KXAINT(ITRAN)  = KEND1
                KYAINT(ITRAN)  = KXAINT(ITRAN) + NMATIJ(ISYZTA)
                KBZADEN(ITRAN) = KYAINT(ITRAN) + NMATAB(ISYZTA)
                KEND1          = KBZADEN(ITRAN) + NT2AOIJ(ISYZTA)
                IF (CCSD) THEN
                   KFAIM(ITRAN) = KEND1
                   KEND1        = KFAIM(ITRAN) + NT1AO(ISYZTA)
                END IF
             ELSE
                KXAINT(ITRAN)  = KXAINT(ITRAN - 1)
                KYAINT(ITRAN)  = KYAINT(ITRAN - 1)
                KBZADEN(ITRAN) = KBZADEN(ITRAN- 1)
                IF (CCSD) THEN
                   KFAIM(ITRAN) = KFAIM(ITRAN - 1)
                END IF
             END IF
*
* IOPER, T^A and Zeta dependent intermediates
*
             IF ((LNEWTA).OR.(LNEWZE).OR.(LNEWOP)) THEN
                KBZQADEN(ITRAN) = KEND1
                KEND1           = KBZQADEN(ITRAN) + NT2AOIJ(ISYRES)
                IF (CCSD) THEN
                   KFQAIM(ITRAN) = KEND1
                   KEND1         = KFQAIM(ITRAN) + NT1AO(ISYRES)
                END IF
             ELSE
                KBZQADEN(ITRAN) = KBZQADEN(ITRAN - 1)
                IF (CCSD) THEN
                   KFQAIM(ITRAN)   =  KFQAIM(ITRAN-1)
                END IF
             END IF
          END IF

          LWRK1   = LWORK - KEND1
          IF (LWRK1 .LT. 0) THEN
             CALL QUIT('Insufficient work space in CC_FBTA. (A)')
          END IF
*        -----------------------------------------------------------
*         get the orbital relaxation vector and the connection matrix 
*         R, and calculate the derivative Lambda matrices "LambdaQ" 
*         and the derivative density matrices:
*         Note, that we only have a orbital relaxation vector if
*         IRELAX.GE.1, the logical LRELAX is true if we have orbital
*         relaxation or a connection matrix, i.e. if we have nonzero
*         lambdaQ matrices.
*        -----------------------------------------------------------
*
          IF (LRELAX) THEN
*
            KKAPPA = KEND1                             !Kappa^(B) Matrix
            KRMAT  = KKAPPA + 2*NALLAI(ISYHOP)         !R^(B) Matrix
            KEND2  = KRMAT  + N2BST(ISYHOP)
            LWRK2  = LWORK  - KEND2
*       
            IF (LWRK2 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_FBTA. (2)')
            END IF

            IF (LABELH.EQ.'HAM0    ' .AND. IRELAX.GE.1) THEN
               IRELAX = 0
               WRITE (LUPRI,*) 
     &           'Test case "HAM0"... no relaxation vector used.'
            END IF                               

* Note that we keep orbital relaxation vector on the 
* same index  list as R1 vector

            IF (IRELAX.GE.1) THEN
               CALL CC_RDHFRSP('R1 ',IRELAX,ISYHOP,WORK(KKAPPA))
            ELSE
               CALL DZERO(WORK(KKAPPA),2*NALLAI(ISYHOP))
            END IF
*
            IF (LOCDBG) THEN
              WRITE (LUPRI,*) 
     &           'CC_FBTA1> orbital relax. vector:',IRELAX,LABELH
              CALL OUTPUT(WORK(KKAPPA),1,2*NALLAI(ISYHOP),1,1,
     &                                 2*NALLAI(ISYHOP),1,1,LUPRI)
              WRITE (LUPRI,*) 'WORK(0): ', WORK(ITST)
            END IF
*
            IORDER = 1 ! might need some generalization later...

            CALL CC_GET_RMAT(WORK(KRMAT),IOPERB,IORDER,ISYHOP,
     &                       WORK(KEND2),LWRK2)
            IOPT = 1
            CALL CC_LAMBDAQ(WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),
     &                      WORK(KCMOPQ(ITRAN)), WORK(KCMOHQ(ITRAN)),
     &                      ISYHOP,WORK(KT1AMP0),
     &                      WORK(KKAPPA),WORK(KRMAT),IREAL,IOPT,
     &                      WORK(KEND2),LWRK2)
            IF (LOCDBG) THEN
              WRITE (LUPRI,*) 
     &           'CC_FBTA> The perturbed LambdaQs (CC_LAMBDAQ)'
              CALL CC_PRLAM(WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),
     &                                                       ISYHOP)
            END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C         THIS is still unfinished. 
C         frozen core contributions have to be introduced!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

* B-perturbed CC-Fock density
            ICORE = 1                      
            IOPT  = 2
            CALL CC_AODENS2(WORK(KLAMP0),WORK(KLAMH0),ISYM0,
     &              WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),ISYHOP,
     &              WORK(KDENSB(ITRAN)),ICORE,IOPT,WORK(KEND2),LWRK2)

            IF (LOCDBG) THEN 
              XNORM = DNRM2(N2BST(ISYHOP),WORK(KDENSB(ITRAN)),1)
              WRITE (LUPRI,*) 
     *         'CC_FBTA> Norm of B prtbd CC-Fock dens. DENSB: ',XNORM
            END IF
*         ------------------------------------------------------------
*         calculate the T^A dependent Lambdas and the Fock densities
*         ------------------------------------------------------------
            KTA1AM = KEND2
            KEND3  = KTA1AM + NT1AM(ISYMTA)
            LWRK3  = LWORK - KEND3
            IF (LWRK3 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_FBTA. (T^A in)')
            END IF
*
            IOPT=1                          
            CALL CC_RDRSP(LISTR,IDLSTR,ISYMTA,IOPT,MODEL,
     &                                      WORK(KTA1AM),WORK(KDUM))
*
* trasform zero-order lambdas with T^A 
*
            CALL CC_LAMTRA(WORK(KLAMP0),ISYM0,
     &                     WORK(KLAMDPA(ITRAN)),ISYMTA,
     &                     WORK(KLAMH0),ISYM0,
     &                     WORK(KLAMDHA(ITRAN)),ISYMTA,
     &                     WORK(KTA1AM),ISYMTA)
*
* trasform Q-lambdas with T^A
*
            CALL CC_LAMTRA(WORK(KLAMDPQ(ITRAN)), ISYHOP,
     &                     WORK(KLAMPQA(ITRAN)), ISYHTA,
     &                     WORK(KLAMDHQ(ITRAN)), ISYHOP,
     &                     WORK(KLAMHQA(ITRAN)), ISYHTA,
     &                     WORK(KTA1AM),ISYMTA)
*
* The DA(CC-Fock) and DBA(CC-Fock) densities
* DENS0A = sum_k Lambda0P_\gk LamdbaAH_k\d
* DENSBA = sum_k (LambdaQP_\gk LamdbaAH_k\d + Lambda0P_\gk LambdaQAH_k\d)
*
            ICORE = 1                     !with frozen core not working?  
            IOPT  = 1
            CALL CC_AODENS2(WORK(KLAMP0),WORK(KLAMH0),ISYM0,
     &         WORK(KLAMDPA(ITRAN)),WORK(KLAMDHA(ITRAN)),ISYMTA,
     &         WORK(KDENS0A(ITRAN)),ICORE,IOPT,WORK(KEND2),LWRK2)

            IF (LOCDBG) THEN
              XNORM = DNRM2(N2BST(ISYMTA),WORK(KDENS0A(ITRAN)),1)
              WRITE (LUPRI,*) 
     &              'CC_FBTA> Norm of DENS0A(CCFOCKA density): ',XNORM
            END IF
*
            ICORE = 1                      !with frozen core not working...
            IOPT  = 2                      !ORDER IS IMPORTANT
            CALL CC_AODENS3(WORK(KLAMP0), ISYM0,
     &                      WORK(KLAMHQA(ITRAN)), ISYHTA,
     &                      WORK(KLAMDPQ(ITRAN)), ISYHOP,
     &                      WORK(KLAMDHA(ITRAN)), ISYMTA,
     &                      WORK(KDENSBA(ITRAN)), ISYHTA,
     &                      ICORE,IOPT,WORK(KEND2),LWRK2)

            IF (LOCDBG) THEN
              XNORM = DNRM2(N2BST(ISYHTA),WORK(KDENSBA(ITRAN)),1)
              WRITE (LUPRI,*) 
     &              'CC_FBTA> Norm of DENSBA(CCFOCKBA dens): ', XNORM
            END IF
*
          END IF                    !if lrelax 
*         --------------------------------------------------------------
*         initialize the remaining intermediates with zero:
*         --------------------------------------------------------------
          IF (.NOT. CCS) THEN
             CALL DZERO(WORK(KRAIM(ITRAN)),NEMAT1(ISYMTA))
             CALL DZERO(WORK(KRQAIM(ITRAN)),NEMAT1(ISYHTA))
          END IF
          IF (LOCDBG) WRITE (LUPRI,*) 'WORK(0): ', WORK(ITST)
*
          IF (CCSD.AND.(LTWOEL.OR.LRELAX)) THEN
            CALL DZERO(WORK(KFQAIM(ITRAN)),NT1AO(ISYRES))
*            IF (LNEWTA.OR.LNEWZE) THEN
                CALL DZERO(WORK(KFAIM(ITRAN)),NT1AO(ISYZTA))
*            END IF
          END IF
*
          IOP_OLD = IOPERB
          IDR_OLD = IDLSTR
          IDL_OLD = IDLSTL
*
          IF (LOCDBG) THEN
             WRITE (LUPRI,*) 
     *             'CC_FBTA> memory alloc. for AO section finished:'
             WRITE (LUPRI,*) 'CC_FBTA> IOPERB    = ',IOPERB
             WRITE (LUPRI,*) 
     *             'CC_FBTA> IDLSTR, IDLSTL    = ',IDLSTR,IDLSTL
             WRITE (LUPRI,*) 'CC_FBTA> LTWOEL,LRELAX = ',LTWOEL,LRELAX
             WRITE (LUPRI,*) 'CC_FBTA> LABELH,ISYHOP = ',LABELH,ISYHOP
             WRITE (LUPRI,*) 'CC_FBTA> KEND1, LWRK1:',KEND1,LWRK1
             WRITE (LUPRI,*) 'CC_FBTA1> WORK(0): ', WORK(ITST)
          END IF
*
      END DO ! ITRAN
*=====================================================================*
* Calculate preintermediates only depending on ZETA and ZETA/IOPERB
* (and containg T_2(0)):
* -- X, Y, M, P and Q intermediates (MO, AO0, AObar) 
* -- The effective densities d and d-bar for the (E1')' term
* For these intermediates we need the amplituded packed in core
* and the Lagrangian multipliers squared. 
* (The Zeta vector is read in CCETAINT1 into WORK(KZETA2) )
* The X,Y,M intermediates are not needed after this
*=====================================================================*
C      WRITE(LUPRI,*) 'point1:'
      IF (.NOT.CCS) THEN
         KT2AMP0  = KEND1
         KEND2    = KT2AMP0  + NT2AM(ISYM0)
         LWRK2    = LWORK    - KEND2
         IF (LWRK2 .LT. 0) 
     &   CALL QUIT('Insufficient work space in CC_FBTA. (CCETAINT1)')
*
         IOPT = 2
         CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KT2AMP0))

         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
C            WRITE(LUPRI,*) 'point2:'
*
            IOPERB = IFBTRAN(1,ITRAN) ! B operator index
            IDLSTL = IFBTRAN(2,ITRAN) ! ZETA vector index
            IRELAX = IFBTRAN(5,ITRAN) ! flag for orbital relaxation
            ISYHOP = ISYOPR(IOPERB)   ! symmetry of hamiltonian
            LTWOEL = LPDBSOP(IOPERB)  ! two-electron contribution
            LRELAX = LTWOEL .OR. (IRELAX.GE.1)
*
            ISYCTR = ILSTSYM(LISTL,IDLSTL)  ! symmetry of ZETA vector
*
            KZETA1 = KEND2
            KZETA2 = KZETA1 + NT1AM(ISYCTR)
            KXINT  = KZETA2 + NT2SQ(ISYCTR)
            KYINT  = KXINT  + NMATIJ(ISYCTR)
            KMINT  = KYINT  + NMATAB(ISYCTR)
            KEND3  = KMINT  + N3ORHF(ISYCTR)

CCH quick hack to avoid core dump in CCETAINT1 after the call list
CCH was changed for CC2, needs to be fixed...
            KZDPK0 = KEND3
            KZDEN0 = KZDPK0 + NNBST(ISYCTR)
            KZDPKB = KZDEN0 + N2BST(ISYCTR)
            KZDENB = KZDPKB + NNBST(MULD2H(ISYCTR,ISYHOP))
            KCHI   = KZDENB + N2BST(MULD2H(ISYCTR,ISYHOP))
            KCHIQ  = KCHI   + NGLMDT(ISYCTR)
            KEND3  = KCHIQ  + NGLMDT(MULD2H(ISYCTR,ISYHOP))
CCH

            LWRK3  = LWORK  - KEND3
            IF (LWRK3 .LT. 0) 
     &        CALL QUIT(
     *           'Insufficient work space in CC_FBTRAN. (CCETAINT1)')
         
C            WRITE(LUPRI,*) 'enter CCETAINT1...'
            CALL CCETAINT1(ITRAN,ISTART(IBATCH),LISTL,IDLSTL,
     &                  WORK(KZDPK0),WORK(KZDEN0),
     &                  WORK(KZDPKB),WORK(KZDENB),
     &                  WORK(KXINT),WORK(KYINT),WORK(KMINT),
     &                  WORK(KCHI),WORK(KCHIQ),
     &                  WORK(KZETA1),WORK(KZETA2),ISYCTR,
     &                  WORK(KT2AMP0),WORK(KLAMP0), WORK(KLAMH0), ISYM0,
     &                WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)), ISYHOP,
     &                FNBFDE0, LUBFDE0, IADRE0,  IADRBFE0,
     &                FNBFDE1, LUBFDE1, IADRE1,  IADRBFE1,
     &                FILPQMO, LUPQMO,  IADRPQMO,IADRPQ,
     &                FILPQ0,  LUPQ0,   IADRPQ0, IADRPQI0,
     &                FILPQ1,  LUPQ1,   IADRPQ1, IADRPQI1,
     &                LRELAX,  LTWOEL,  WORK(KEND3),LWRK3)
C            WRITE(LUPRI,*) 'returned from CCETAINT1...'
*
         END DO
      END IF
*
*=====================================================================*
* Precalculate intermediates depending on both ZETA and T^A 
* and on ZETA/TA/IOPERB:
* -- XA and YA (in memory), PA and QA intermediates (MO, AO, AObar;disk) 
* -- the effective densities D and D-bar for the rho^BZA,rho^BZQA interm.
*    These effective densities are kept in memory 
*
* We need the amplitudes (TA2) packed in core and the Lagrangian 
* multipliers (Zeta2) squared.  Both are read inside, but allocated here
*=====================================================================*

      IF (.NOT.CCS) THEN
         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)

            IOPERB  = IFBTRAN(1,ITRAN) ! B operator index
            IDLSTL = IFBTRAN(2,ITRAN) ! ZETA vector index
            IDLSTR = IFBTRAN(3,ITRAN) ! T^A  vector index
            IRELAX = IFBTRAN(5,ITRAN) ! flag for orbital relaxation
            ISYHOP = ISYOPR(IOPERB)    ! symmetry of hamiltonian
            LTWOEL = LPDBSOP(IOPERB)   ! two-electron contribution
            LRELAX = LTWOEL .OR. (IRELAX.GE.1)

            ISYCTR = ILSTSYM(LISTL,IDLSTL)  ! symmetry of ZETA vector
            ISYMTA = ILSTSYM(LISTR,IDLSTR)  ! symmetry of T^A vector
            ISYZTA = MULD2H(ISYCTR,ISYMTA)
            ISYHTA = MULD2H(ISYHOP,ISYMTA)

            KZETA1 = KEND2
            KZETA2 = KZETA1 + NT1AM(ISYCTR)
            KTA1AM = KZETA2 + NT2SQ(ISYCTR) !Z2 are squared 
            KTA2AM = KTA1AM + NT1AM(ISYMTA)
            KMAINT = KTA2AM + NT2AM(ISYMTA)
            KEND3  = KMAINT + N3ORHF(ISYZTA)
            LWRK3  = LWORK  - KEND3
*
            CALL CCFBINT1(ITRAN, LISTL, IDLSTL, LISTR, IDLSTR,
     &               WORK(KXAINT(ITRAN)),WORK(KYAINT(ITRAN)),
     &               WORK(KMAINT),
     &               WORK(KZETA1), WORK(KZETA2), ISYCTR,
     &               WORK(KTA1AM), WORK(KTA2AM), ISYMTA,
     &               WORK(KLAMP0), WORK(KLAMH0), ISYM0,               
     &               WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)), ISYHOP,
     &               WORK(KBZADEN(ITRAN)),WORK(KBZQADEN(ITRAN)),
     &               FILPQAMO, LUPQAMO, IADPQAMO, IADRPQA,
     &               FILPQA0,  LUPQA0,  IADRPQA0,  IADRPQAI0,
     &               FILPQA1,  LUPQA1,  IADRPQA1,  IADRPQAI1,
     &               LRELAX,   LTWOEL,  WORK(KEND3), LWRK3)

         END DO
      END IF
*=====================================================================*
* Precalculate intermediates depending only on T^A or/and IOPERB: 
* -- the AO-FOCKB  (only initialized here with h^(1))
* -- the AO-FOCKA  (only initialized here with zeros)
* -- the AO-FOCKBA (only initialized here with zeros)
* -- the effective density (Delta) for the rho^BFA interm.  (if LNEWTA)
* -- the effective density (Delta) for the rho^BFQA interm. (if LNEWTA.or.LNEWOP)
*
* Move maybe calculation of Fock densities inside it
* As it stays now, only Lambdas^h are needed (0,A,QA)
*=====================================================================*
*
      IOP_OLD = -1
      IDR_OLD  = -1
*
      DO ITRAN = ISTART(IBATCH), IEND(IBATCH)

         IOPERB = IFBTRAN(1,ITRAN) ! operator index
         IDLSTR = IFBTRAN(3,ITRAN) ! index for T^A
         IRELAX = IFBTRAN(5,ITRAN) ! flag for orbital relaxation
         ISYHOP = ISYOPR(IOPERB)    ! symmetry of hamiltonian
         LABELH = LBLOPR(IOPERB)    ! label for 1-el integrals
         LTWOEL = LPDBSOP(IOPERB)   ! two-electron contribution
         LRELAX = LTWOEL .OR. (IRELAX.GE.1)
*
         ISYMTA = ILSTSYM(LISTR,IDLSTR)
         ISYHTA = MULD2H(ISYHOP,ISYMTA)
*
         LNEWTA = (IDLSTR .NE. IDR_OLD) 
         LNEWOP = (IOPERB .NE. IOP_OLD)

         KTA2AM = KEND2                         !use squared T^A_2
         KEND3  = KTA2AM + NT2SQ(ISYMTA)
         LWRK3  = LWORK  - KEND3

         IF (LWRK3.LT.NT2AM(ISYMTA)) THEN
           CALL QUIT('Insufficient memory in CC_FBTA (LNEWTA)')
         END IF
*
* read it here for the time being
*
         IF (LNEWTA) THEN
           IOPT = 2
           CALL CC_RDRSP(LISTR,IDLSTR,ISYMTA,IOPT,MODEL,
     &                              WORK(KDUM),WORK(KEND3))
           CALL CCLR_DIASCL(WORK(KEND3),TWO,ISYMTA)
           CALL CC_T2SQ(WORK(KEND3),WORK(KTA2AM),ISYMTA)
         END IF
*
         CALL CCFBINT2(ITRAN, LABELH, WORK(KTA2AM),
     &             WORK(KDPCKB(ITRAN)),WORK(KONEHB(ITRAN)),
     &             WORK(KFOCKBCC(ITRAN)),WORK(KDENSB(ITRAN)),
     &             WORK(KDENS0A(ITRAN)),WORK(KDPCK0A(ITRAN)),
     &             WORK(KFOCKACC(ITRAN)),
     &             WORK(KDENSBA(ITRAN)),WORK(KDPCKBA(ITRAN)),
     &             WORK(KFCKBACC(ITRAN)),
     &             WORK(KLAMH0), ISYM0,
     &             WORK(KLAMDHQ(ITRAN)), ISYHOP,
     &             WORK(KLAMDHA(ITRAN)), ISYMTA,
     &             WORK(KLAMHQA(ITRAN)), ISYHTA,
     &             FNBFDA,  LUBFDA,  IADRBFA, IADBFA,
     &             FNBFDQA, LUBFDQA, IADRBFQA, IADBFQA,
     &             LRELAX,  LTWOEL,  LNEWTA,  LNEWOP,
     &             WORK(KEND3),LWRK3)

        IOP_OLD = IOPERB
        IDR_OLD = IDLSTR

      END DO

* Allocate and initialize rho^BF(A,QA)_{alphal,nm} (move up)

      DO ITRAN = ISTART(IBATCH), IEND(IBATCH)

         IOPERB = IFBTRAN(1,ITRAN) ! B operator index
         IDLSTL = IFBTRAN(2,ITRAN) ! ZETA vector index
         IDLSTR = IFBTRAN(3,ITRAN) ! T^A  vector index
         IRELAX = IFBTRAN(5,ITRAN) ! flag for orbital relaxation
         ISYHOP = ISYOPR(IOPERB)    ! symmetry of hamiltonian
         LTWOEL = LPDBSOP(IOPERB)   ! two-electron contribution
         LRELAX = LTWOEL .OR. (IRELAX.GE.1)

         ISYCTR = ILSTSYM(LISTL,IDLSTL)  ! symmetry of ZETA vector
         ISYMTA = ILSTSYM(LISTR,IDLSTR)  ! symmetry of T^A vector
         ISYZTA = MULD2H(ISYCTR,ISYMTA)
         ISYHTA = MULD2H(ISYHOP,ISYMTA)

         KRHOBFA(ITRAN)  = KEND2
         KRHOBFQA(ITRAN) = KRHOBFA(ITRAN) + NT2AOIJ(ISYMTA)
         KEND2           = KRHOBFQA(ITRAN) + NT2AOIJ(ISYHTA)
         LWRK2           = LWORK  - KEND2
*      
         CALL DZERO(WORK(KRHOBFA(ITRAN)),NT2AOIJ(ISYMTA))
         CALL DZERO(WORK(KRHOBFQA(ITRAN)),NT2AOIJ(ISYHTA))
      END DO

*---------------------------------------------------------------------*
* initialize integral program:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_FBTA: Initialize now integral program:'
         WRITE (LUPRI,*) 'LLTWOEL = ',LLTWOEL
         WRITE (LUPRI,*) 'LABELH  = ',LABELH     !where does it get it from?
      END IF
*
      IF (LLTWOEL .AND.
     &    (LABELH(1:5).EQ.'1DHAM' .OR. LABELH(1:5).EQ.'dh/dB') ) THEN
*
*        derivative integrals
*
         SYM1ONLY = .FALSE.
         CALL CC_SETDORPS(LABELH,SYM1ONLY,IPRINT)

         IOPERB  = IFBTRAN(1,ISTART(IBATCH)) !assume batching according to IOPERB?
         IATOM   = IATOPR(IOPERB)

         IF (.NOT.DIRECT) THEN
           CALL CCDER1(IATOM,LABELH,LDERINT,WORK(KEND2),LWRK2)
         END IF

      END IF

      IF (LLRELAX .OR. LABELH(1:5).EQ.'HAM0 ') THEN
          
         IF (DIRECT) THEN
            NTOSYM = 1
            IF (HERDIR) THEN
              CALL HERDI1(WORK(KEND2),LWRK2,IPRERI)
            ELSE
              KCCFB1 = KEND2
              KINDXB = KCCFB1 + MXPRIM*MXCONT
              KEND   = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
              LWRK   = LWORK  - KEND
              KEND2  = KEND
              CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                    KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                    KFREE,LFREE,KEND,WORK(KCCFB1),WORK(KINDXB),
     *                    WORK(KEND2),LWRK,IPRERI)
            END IF
            KEND2 = KFREE
            LWRK2 = LFREE
         ELSE
            NTOSYM = NSYM
         END IF

      END IF

*---------------------------------------------------------------------*
* loop over integral distributions:
* (nr of delta obtained simultaneously)
*---------------------------------------------------------------------*
      IF (LLTWOEL .OR. LLRELAX) THEN

        DO ISYMD1 = 1, NTOSYM                !=NSYM for nondirect, 1 for direct     

          IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            END IF
          ELSE
            NTOT = NBAS(ISYMD1)
          END IF

          DO ILLL = 1, NTOT
     
            IF (DIRECT) THEN

              IF      (LABELH(1:5).EQ.'HAM0 ' .OR. (.NOT.LLTWOEL)) THEN
                NGDER = 0
                NBDER = 0
                NFILES = 1
              ELSE IF (LABELH(1:5).EQ.'1DHAM') THEN
                NGDER = 1
                NBDER = 0
                NFILES = 1 + 3*NUCDEP
              ELSE IF (LABELH(1:5).EQ.'dh/dB') THEN
                NGDER = 0
                NBDER = 1
                NFILES = 1 + 3*2
              ELSE
                CALL QUIT('Unknown 2e- integral type in CC_FBTA.')
              END IF

              KEND = KEND2
              LWRK = LWRK2

C             IF (HERDIR) THEN
C                CALL HERDI2(WORK(KEND),LWRK,INDEXA,ILLL,NUMDIS,
C    &                       IPRINT)
C             ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,NGDER,NBDER,
     *                       WORK(KODCL1),WORK(KODCL2),
     *                       WORK(KODBC1),WORK(KODBC2),
     *                       WORK(KRDBC1),WORK(KRDBC2),
     *                       WORK(KODPP1),WORK(KODPP2),
     *                       WORK(KRDPP1),WORK(KRDPP2),
     *                       WORK(KCCFB1),WORK(KINDXB),
     *                       WORK(KEND), LWRK,IPRERI)
C             END IF

              NBUFMX = NBUFX(0)
              DO I = 1, NFILES-1
                NBUFMX = MAX(NBUFMX,NBUFX(I))
              END DO

              ! allocate memory for orbital/record labels
              KRECNR = KEND
              KNRECS = KRECNR + (NBUFMX*NFILES - 1)/IRAT + 1
              KEND3  = KNRECS + (NFILES - 1)/IRAT + 1
              LWRK3  = LWORK - KEND3

              IF (LWRK3 .LT. 0) THEN
               CALL QUIT('Insufficient work space in CC_FBTA. (ERIDI2)')
              END IF

              CALL RDERILBS(WORK(KRECNR),WORK(KNRECS),NBUFMX,NFILES)
* 
            ELSE
              NUMDIS = 1
              KEND3  = KEND2
              LWRK3  = LWRK2
            END IF
*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_FBTA> after integral evaluation '//
     &        'for ILLL:',ILLL
        WRITE (LUPRI,*) 'CC_FBTA> KEND3, LWRK3:',KEND3,LWRK3
        WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
      END IF
*
*=====================================================================*
* Calculate intermediates depending on the two-electron AO integrals,
* in the loop over AO integral shells. We put the loop over the result
* vectors here, because the required derivatives integrals can change
* for each ITRAN. We assume that all RES vectors for operators
* coming from the same IATOM are grouped together. In this way we
* avoid recalculation of the derivative integrals for each vector, but
* accept that the AO integrals are reread from file for each ITRAN.
*=====================================================================*
*
      IOP_OLD = -1
      IDR_OLD = -1
      IDL_OLD = -1
*
      DO ITRAN = ISTART(IBATCH), IEND(IBATCH)

        IOPERB = IFBTRAN(1,ITRAN)   ! operator index
        IDLSTL = IFBTRAN(2,ITRAN)   ! ZETA vector index
        IDLSTR = IFBTRAN(3,ITRAN)   ! T^A  vector index
        ID_RES = IFBTRAN(4,ITRAN)   ! index for RES vector
        IRELAX = IFBTRAN(5,ITRAN)   ! flag for orbital relaxation
*
        LABELH = LBLOPR(IOPERB)      ! operator label
        IATOM  = IATOPR(IOPERB)      ! associated atom
        ISYHOP = ISYOPR(IOPERB)      ! operator symmetry
        IREAL  = ISYMAT(IOPERB)      ! flag for real/imag. operators
        LTWOEL = LPDBSOP(IOPERB)     ! pert.-dep. basis set
        LRELAX = LTWOEL .OR. (IRELAX.GE.1)
*
        LNEWTA = (IDLSTR.NE.IDR_OLD)
        LNEWZE = (IDLSTL.NE.IDL_OLD)
*
        ! number of components (directions) in total for this operator,
        ! coordinate index (ICOOR) and coordinate symmetry (ICORSY)
*
        MXCOMP = 0
        IF      ( LABELH(1:5).EQ.'HAM0 ' ) THEN
           MXCOMP = 1
           ICOOR  = 1
           ICORSY = 1
           ISCOOR = 0
        ELSE IF ( LABELH(1:5).EQ.'dh/dB' ) THEN
           MXCOMP = 3
           ICORSY = ISYHOP
           DO JCOOR = 1, MXCOMP
              IF (CHRXYZ(JCOOR).EQ.LABELH(6:6)) THEN
                 ICOOR = JCOOR
              END IF
           END DO
           ISCOOR = ICOOR
        ELSE IF ( LABELH(1:5).EQ.'1DHAM' ) THEN
           MXCOMP = 3
           ICORSY = ISYHOP
           READ(LABELH,'(A5,I3)') LAB1,ISCOOR

           DO JCOOR = 1, MXCOMP
              JSCOOR = IPTCNT(3*(IATOM-1)+JCOOR,ICORSY-1,1)
              IF (JSCOOR.EQ.ISCOOR) THEN
                 ICOOR = JCOOR
              END IF
           END DO
        ELSE IF (LTWOEL) THEN
           WRITE (LUPRI,*) 'Unknown 2el operator label in CC_FBTA.'
           CALL QUIT('Unknown 2el operator label in CC_FBTA.')
        END IF
*
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_FBTA> IOPERB        = ',IOPERB
          WRITE (LUPRI,*) 'CC_FBTA> LTWOEL,LRELAX = ',LTWOEL,LRELAX
          WRITE (LUPRI,*) 'CC_FBTA> ID_RES        = ',ID_RES
          WRITE (LUPRI,*) 'CC_FBTA> LABELH,ISYHOP = ',LABELH,ISYHOP
          WRITE (LUPRI,*) 'CC_FBTA> IATOM,ISCOOR  = ',IATOM,ISCOOR
          WRITE (LUPRI,*) 'CC_FBTA> ICORSY,ICOOR  = ',ICORSY,ICOOR
        END IF
*
*---------------------------------------------------------------------*
*     loop over number of distributions in this shell:
*     (we enter this loop only for those ITRAN where we have
*      a two-electron contribution, i.e. if LRELAX is true. )
*---------------------------------------------------------------------*
*
       IF (LRELAX) THEN

         DO IDEL2 = 1, NUMDIS
*
*         set orbital index and symmetry class for next 
*         AO integral distribution:
*
          IF (DIRECT) THEN              !find  sym and index for given delta
            IDEL   = INDEXA(IDEL2)
            IF (NOAUXB) THEN
               IDUM1 = 1
               CALL IJKAUX(IDEL,IDUM1,IDUM1,IDUM1)
            END IF
            ISYDEL = ISAO(IDEL)
          ELSE
            IDEL   = IBAS(ISYMD1) + ILLL
            ISYDEL = ISYMD1
          END IF
*
          ISY0DIS = MULD2H(ISYDEL,ISYM0)         !sym of I_alp-bet,gam (del)
          ISY1DIS = MULD2H(ISYDEL,ISYHOP)        !sym of I(1)_alp-bet,gam (del)
*
*         read undifferentiated AO integral distribution,
*         and transform gamma index to occupied:
*         note that ISYD0SRHF = ISY0DIS*ISYM0 = ISY0DIS
*         and ISYD1SRHF = ISY1DIS*ISYM0 = ISY0DIS*ISYHOP
*
          K0XINT  = KEND3
          KD0SRHF = K0XINT  + NDISAO(ISY0DIS)     !the I_alp-bet,K;delta
          KEND4   = KD0SRHF + NDSRHF(ISY0DIS)
          LWRK4   = LWORK   - KEND4
          IF (LWRK4 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_FBTA. (CCRDAO)')
          END IF
*
COLD
C         CALL CCRDAO(WORK(K0XINT),IDEL,IDEL2,
C    &                WORK(KEND4),LWRK4,WORK(KRECNR),DIRECT)
COLD
C         read all gamma for given delta (``3-index'' approach)
C
          NUMD = 1
          NUMG = 0
          DO ISYMG = 1, NSYM
            DO G = 1, NBAS(ISYMG)
              NUMG       = NUMG + 1
              IGAM(NUMG) = G + IBAS(ISYMG)
            END DO
          END DO

          DIRINT = DIRECT
          CALL CC_RDAOD(WORK(K0XINT),0,1,1,IGAM,IDEL,NUMG,NUMD,
     &                  WORK(KEND4),LWRK4,WORK(KRECNR),WORK(KNRECS),
     &                  NBUFMX,NFILES,DIRINT,0,1,LDERINT)


          CALL CCTRBT(WORK(K0XINT),WORK(KD0SRHF),WORK(KLAMP0),
     &                        ISYM0,WORK(KEND4),LWRK4,ISY0DIS)

          IF (LTWOEL) THEN
             IF       ( LABELH(1:5).EQ.'HAM0 ' ) THEN
               ITYPE  = 0
               SQRINT = .FALSE.
             ELSE IF ( LABELH(1:5).EQ.'1DHAM' ) THEN
               ITYPE  = 1
               SQRINT = .FALSE.
             ELSE IF ( LABELH(1:5).EQ.'dh/dB' ) THEN
               ITYPE  = 5
               SQRINT = .TRUE.
             ELSE
               CALL QUIT('ITYPE unknown in CC_FBTA.')
             END IF
          ELSE                                       
             SQRINT = .FALSE.     !use SQRINT in BFBSORT
          END IF
          
          KD1SRHF = KEND4         !the I_alp-bet,bar-K;delta (relaxed)
          IF ((LTWOEL).AND.(SQRINT)) THEN
             KEND4   = KD1SRHF + NDSRHFSQ(ISY1DIS)
             LSQRUP  = .TRUE.
          ELSE 
             KEND4   = KD1SRHF + NDSRHF(ISY1DIS)
             LSQRUP  = .FALSE.
          END IF

          IF (LWRK4 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_FBTA. (CCTRBT2)')
          END IF
*

          IOPT = 0
          CALL CCTRBT2(WORK(K0XINT),WORK(KD1SRHF),WORK(KLAMDPQ(ITRAN)),
     &                 ISYHOP,WORK(KEND4),LWRK4,ISY0DIS,IOPT,
     &                 .FALSE.,LSQRUP,ONE)

          IF (LTWOEL) THEN
C
C            read first derivative AO integral distribution:
C            transform gamma index to occupied
C            and add to KD1SRHF
C
             K1XINT  = KEND4                           
             IF (SQRINT) THEN
                KEND4 = K1XINT + NDISAOSQ(ISY1DIS)
             ELSE
                KEND4 = K1XINT + NDISAO(ISY1DIS)
             END IF
             LWRK4   = LWORK   - KEND4
*
             IF (LWRK4 .LT. 0) THEN
               CALL QUIT(
     *           'Insufficient work space in CC_FBTA. (CCRDAOD)')
             END IF                              
*
*            read all gamma for given delta (``3-index'' approach)
*
             NUMD = 1
             NUMG = 0
             DO ISYMG = 1, NSYM
               DO G = 1, NBAS(ISYMG)
                 NUMG       = NUMG + 1
                 IGAM(NUMG) = G + IBAS(ISYMG)
               END DO
             END DO
*
             DIRINT = DIRECT
             CALL CC_RDAOD(WORK(K1XINT),ISCOOR,ICOOR,ICORSY,
     &                     IGAM,IDEL,NUMG,NUMD,
     &                     WORK(KEND4),LWRK4,WORK(KRECNR),WORK(KNRECS),
     &                     NBUFMX,NFILES,DIRINT,ITYPE,MXCOMP,LDERINT)
  
             IF (LOCDBG) THEN
                WRITE (LUPRI,*) 'Norm of undifferentiated integrals:',
     &          DNRM2(NDISAO(ISY0DIS),WORK(K0XINT),1)
              IF (SQRINT) THEN
                WRITE (LUPRI,*) 
     &             'Norm of derivative integrals (squared):',
     &              DNRM2(NDISAOSQ(ISY1DIS),WORK(K1XINT),1)
              ELSE
                WRITE (LUPRI,*) 
     &             'Norm of derivative integrals (packed):',
     &              DNRM2(NDISAO(ISY1DIS),WORK(K1XINT),1)
              END IF
             END IF

*
* transform gamma to K and cumulate result to g_al_bet,bar-K,delta 
*
             IOPT = 1
             CALL CCTRBT2(WORK(K1XINT),WORK(KD1SRHF),WORK(KLAMP0),
     &                ISYM0,WORK(KEND4),LWRK4,ISY1DIS,IOPT,SQRINT,
     &                .FALSE.,ONE)
*
* WORK(K1XINT) should now contain g^(1)_al_bet,K,delta for the given ITRAN
*
          END IF

* Presort integrals for rho^BF(A,QA) intermediates KBF0RHF, KBF1RHF
* Maybe do inside FBTAAO1

          KBF0RHF = KEND4
          KBF1RHF = KBF0RHF + NBSRHF(ISY0DIS)
          KEND4   = KBF1RHF + NBSRHF(ISY1DIS)
          LWRK4   = LWORK - KEND4

          CALL CC_BFBSORT(WORK(KD0SRHF),WORK(KBF0RHF),ISY0DIS)
*
          CALL CC_BFBSORT1(WORK(KD1SRHF),WORK(KBF1RHF),ISY1DIS,SQRINT)

*---------------------------------------------------------------------*
* Calculate the (ij,k;del) integrals from the (al bet|k;del)
* for the E2'_ai contribution. Write result on file
* (can also be kept in memory)
*---------------------------------------------------------------------*
          LZERO = (ITRAN.EQ.1)

          KX3MO = KEND4
          KEND4 = KX3MO + LX3MO
          LWRK4 = LWORK - KEND4

          IF (LZERO) THEN
            LWRDSK = .TRUE.
            FAC    =  ZERO
            CALL CC_INT3O2(WORK(KX3MO),WORK(KD0SRHF),ISY0DIS,
     *                  WORK(KLAMH0),ISYM0,WORK(KLAMP0),ISYM0,
     *                  WORK(KEND4),LWRK4,IDEL,ISYDEL,
     *                  FAC,LWRDSK,LU0IJKL,FN0IJKL,ITRAN,.FALSE.)
          END IF 

          FAC = ZERO
          LWRDSK = .FALSE.
          CALL CC_INT3O2(WORK(KX3MO),WORK(KD0SRHF),ISY0DIS,
     *                  WORK(KLAMDHQ(ITRAN)),ISYHOP,WORK(KLAMP0),ISYM0,
     *                  WORK(KEND4),LWRK4,IDEL,ISYDEL,
     *                  FAC,LWRDSK,IDUM,CDUM,ITRAN,.FALSE.)
          FAC = ONE
          LWRDSK = .FALSE.
          CALL CC_INT3O2(WORK(KX3MO),WORK(KD0SRHF),ISY0DIS,
     *                  WORK(KLAMH0),ISYM0,WORK(KLAMDPQ(ITRAN)),ISYHOP,
     *                  WORK(KEND4),LWRK4,IDEL,ISYDEL,
     *                  FAC,LWRDSK,IDUM,CDUM,ITRAN,.FALSE.)
          LWRDSK = .TRUE.
          CALL CC_INT3O2(WORK(KX3MO),WORK(KD1SRHF),ISY1DIS,
     *                  WORK(KLAMH0),ISYM0,WORK(KLAMP0),ISYM0,
     *                  WORK(KEND4),LWRK4,IDEL,ISYDEL,
     *                  FAC,LWRDSK,LU1IJKL,FN1IJKL,ITRAN,SQRINT)
*
*---------------------------------------------------------------------*
* Calculate intermediates that do NOT contain ZETA: 
*         - some zeroth-order intermediates need only be computed once
*         - the others have to be computed whenever IOPERB(?) changes
*
*cs The BF intermediate, 
*cs the F and G for CC2 (use one-MO tranformed integrals)
* -- the integrals (IA|Jdelta) (LZERO)
* -- the special integrals (IJCdel)^A, (CJIdel)^A (LNEWTA)
* -- the integrals (IA|Jdelta)^(1) (always)
* -- the special integrals (IJCdel)^A(1) (always)
* Il counter viene passato quando il vettore e' calcolato e
* salvato su file la prima volta. Non richiesto in lettura
* (a meno che non sia indice di posizione in scrittura)
*---------------------------------------------------------------------*
*
          LZERO = (ITRAN.EQ.1)

c          IF ( (IOPERB .NE. IOP_OLD) .OR. LZERO ) THEN
*
            CALL CCFBTAAO1(WORK(K0XINT),ISY0DIS,WORK(K1XINT),ISY1DIS,
     &                     WORK(KBF0RHF),WORK(KBF1RHF),
     &                     WORK(KDENS0),WORK(KDPK0),WORK(KFOCK0),
     &                     WORK(KDENSB(ITRAN)),WORK(KDPCKB(ITRAN)),
     &                     WORK(KFOCKBCC(ITRAN)),
     &                     WORK(KDENS0A(ITRAN)),WORK(KDPCK0A(ITRAN)),
     &                     WORK(KFOCKACC(ITRAN)),
     &                     WORK(KDENSBA(ITRAN)),WORK(KDPCKBA(ITRAN)),
     &                     WORK(KFCKBACC(ITRAN)),
     &                     WORK(KLAMP0),WORK(KLAMH0),
     &                     WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),
     &                     WORK(KLAMDPA(ITRAN)),WORK(KLAMDHA(ITRAN)),
     &                     WORK(KLAMPQA(ITRAN)),WORK(KLAMHQA(ITRAN)),
     &                     WORK(KRHOBFA(ITRAN)),WORK(KRHOBFQA(ITRAN)),
     &                     LUBFDA,FNBFDA,IADRBFA(1,ITRAN),
     &                     LUBFDQA,FNBFDQA,IADRBFQA(1,ITRAN),
     &                     LU0IAJB,FN0IAJB,LU1IAJB,FN1IAJB,
     &                     IT2DEL0,IADRI0,IT2DELB(1,ITRAN),IADRIB,
     &                     LU0IJCB,LU0CJIB,FN0IJCB,FN0CJIB,
     &                     LU1IJCB,LU1CJIB,FN1IJCB,FN1CJIB,
     &                     IT2DELA(1,ITRAN),IADRI0A,
     &                     IT2DELBA(1,ITRAN),IADRIAB,
     &                     IDEL,LZERO,LNEWTA,LRELAX,LTWOEL,
     &                     SQRINT,IREAL,
     &                     ISYHOP,ISYMTA,WORK(KEND4), LWRK4)
c          END IF
*---------------------------------------------------------------------*
*  Calculate additional intermediates that do contain ZETA: 
*         - some zeroth-order intermediates need only be computed once?
*         - the others have to be computed whenever IDXETA changes, 
*           which should be the case for every ITRAN
*---------------------------------------------------------------------*
*
          LZERO = (ITRAN.EQ.1)

          ISYCTR = ILSTSYM(LISTL,IDLSTL)  ! symmetry of ZETA vector
          ISYZTA = MULD2H(ISYMTA,ISYCTR)
          IF (CCS) GOTO 222               ! nothing to do for CCS

          CALL CCFBTAAO2(WORK(K0XINT),ISY0DIS,WORK(K1XINT),ISY1DIS,
     &                   WORK(KBZADEN(ITRAN)), WORK(KBZQADEN(ITRAN)), 
     &                   WORK(KBF0RHF),WORK(KBF1RHF),
     &                   WORK(KLAMP0), WORK(KLAMH0), 
     &                   WORK(KLAMDPQ(ITRAN)), WORK(KLAMDHQ(ITRAN)), 
     &                   WORK(KLAMDPA(ITRAN)), WORK(KLAMDHA(ITRAN)), 
     &                   WORK(KLAMPQA(ITRAN)), WORK(KLAMHQA(ITRAN)), 
     &                   WORK(KFAIM(ITRAN)), WORK(KFQAIM(ITRAN)),
     &                   FILBFZA,  LUBFZA,  IADRBZA(1,ITRAN), IADRBZ1, 
     &                   FILBFZQA, LUBFZQA, IADRBZQA(1,ITRAN),IADRBZ2, 
     &                   FNBFDE0, LUBFDE0, IADRE0(1,ITRAN),
     &                   FNBFDE1, LUBFDE1, IADRE1(1,ITRAN),
     &                   FILPQ0,  LUPQ0,   IADRPQ0,  
     &                   FILPQ1,  LUPQ1,   IADRPQ1(1,ITRAN),  
     &                   FILPQA0,  LUPQA0,   IADRPQA0(1,ITRAN),  
     &                   FILPQA1,  LUPQA1,   IADRPQA1(1,ITRAN),  
     &                   LRELAX,  LTWOEL,  LZERO, LNEWTA, LNEWZE,
     &                   SQRINT,  ISYCTR, ISYMTA, ISYHOP, IREAL,
     &                   IDEL, WORK(KEND4), LWRK4)

 222     CONTINUE      
*---------------------------------------------------------------------*
*     close loop over AO two-electron integrals and ITRAN
*---------------------------------------------------------------------*
        END DO ! IDEL2
        END IF ! LRELAX
*
        IOP_OLD = IOPERB
        IDR_OLD = IDLSTR


      END DO ! ITRAN 
      END DO ! ILLL
      END DO ! ISYMD1
      END IF ! (LLTWOEL .OR. LLRELAX) 

*---------------------------------------------------------------------*
* finished with intermediates depending two-electron AO integrals
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_FBTA> after loop over AO'
         CALL FLSHFO(LUPRI)
      END IF
*---------------------------------------------------------------------*
* close & delete the files with effective densities: no longer needed
*---------------------------------------------------------------------*
      CALL WCLOSE2(LUBFDA,  FNBFDA,  'DELETE')
      CALL WCLOSE2(LUBFDQA, FNBFDQA, 'DELETE')
      CALL WCLOSE2(LUBFDE0, FNBFDE0, 'DELETE')
      CALL WCLOSE2(LUBFDE1, FNBFDE1, 'DELETE')
*---------------------------------------------------------------------*
* write the rho^BF(A) and rho^BF(QA) intermediates to newly opened file: 
* but reuse the (1.st row) address-array of the effective densities
* to save some memory
*---------------------------------------------------------------------*
      IF ((LLTWOEL .OR. LLRELAX).AND.(CCSD)) THEN

         LUBFA  = -1
         LUBFQA = -1
         CALL WOPEN2(LUBFA,  FILBFA,  64, 0)
         CALL WOPEN2(LUBFQA, FILBFQA, 64, 0)
*
* the rho^BF(A) (symmetry T^A)
* 
         IADR1 = 1
         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
           IDLSTR  = IFBTRAN(3,ITRAN)
           ISYMTA = ILSTSYM(LISTR,IDLSTR)      ! symmetry of TA
           LENBFA = NT2AOIJ(ISYMTA)            ! length of BFA (at each ITRAN)
C
           CALL PUTWA2(LUBFA,FILBFA,WORK(KRHOBFA(ITRAN)),
     *                                     IADR1,LENBFA)
           IF (LOCDBG) THEN
             XNORM = DNRM2(LENBFA,WORK(KRHOBFA(ITRAN)),1)
             WRITE (LUPRI,*) 'CC_FBTA> Norm of varrho^BFA: ', XNORM
           END IF
           IADRBFA(1,ITRAN) = IADR1                !(reusing array)
           IADR1 = IADR1+ LENBFA                   !use new integer array
         END DO
c
c the rho^BF(QA) (symmetry T^A*Oper)
c 
         IADR2 = 1
         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
           IOPERB = IFBTRAN(1,ITRAN)      ! operator index
           IDLSTR = IFBTRAN(3,ITRAN)      ! TA
           ISYMTA = ILSTSYM(LISTR,IDLSTR) ! symmetry of TA
           ISYHOP = ISYOPR(IOPERB)        ! symmetry of hamiltonian
           ISYHTA = MULD2H(ISYMTA,ISYHOP) ! symmetry TA*Oper (TA*Q)

           LENBFQA = NT2AOIJ(ISYHTA)       ! length of BFQA intermediate

           CALL PUTWA2(LUBFQA,FILBFQA,WORK(KRHOBFQA(ITRAN)),
     *                                      IADR2,LENBFQA)
           IF (LOCDBG) THEN
             XNORM = DNRM2(LENBFQA,WORK(KRHOBFQA(ITRAN)),1)
             WRITE (LUPRI,*) 'CC_FBTA> Norm of varrho^BFQA: ', XNORM
           END IF
           IADRBFQA(1,ITRAN) = IADR2
           IADR2 = IADR2+ LENBFQA
         END DO

      END IF ! (LLTWOEL .OR. LLRELAX) 
*---------------------------------------------------------------------*
*     Calculate RA intermediate from (ia|j del) integrals
*     and RA-bar intermediate from (ia|j del)-bar integrals
*     Not done very smart, ie at each transformation!!!
* The calculation of RA,RQA might be avoided completely for CCSD
*---------------------------------------------------------------------*
C  RA_cdel may change for each IDLSTR (T^A) only (no Zeta)
C  RAbar_cdel may change for each IOPERB (T^A,g^[1]) 
C  Loop over ITRAN first: there are usually many more deltas than NTRAN

      IF (LLTWOEL .OR. LLRELAX) THEN

         DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
 
           IOPERB  = IFBTRAN(1,ITRAN)      ! operator index
           IDLSTR  = IFBTRAN(3,ITRAN)      ! T^A index
           ISYHOP  = ISYOPR(IOPERB)        ! symmetry of hamiltonian
           ISYMTA  = ILSTSYM(LISTR,IDLSTR) ! symmetry of T^A
           ISYHTA  = MULD2H(ISYHOP,ISYMTA)
 
           KTA2AM = KEND1
           KEND2  = KTA2AM + NT2SQ(ISYMTA)
           LWRK2  = LWORK - KEND2
 
           KXIAJB = KEND2                    !reused for both g(0) and g(1)
 
           IF (LWRK2 .LT. MAX(NT2AM(ISYMTA),MT2BCD) ) THEN
             CALL QUIT(
     &         'Insufficient work space in CC_FBTA. (RIM section)')
           END IF

C        ---------------------------------------------------
C        read response cluster amplitudes, square up and
C        calculate 2T^A(ia,jb) - T^A(ja,ib) in place:
C        ---------------------------------------------------
           IOPT = 2
           CALL CC_RDRSP(LISTR,IDLSTR,ISYMTA,IOPT,MODEL,WORK(KDUM),
     &                                                 WORK(KEND2))
c           CALL CC_PRP(WORK(KDUM),WORK(KEND2),ISYMTA,0,1)
           CALL CCLR_DIASCL(WORK(KEND2),TWO,ISYMTA)

           CALL CC_T2SQ(WORK(KEND2),WORK(KTA2AM),ISYMTA)
           CALL CCRHS_T2TR(WORK(KTA2AM),WORK(KEND2),LWRK2,ISYMTA)

           DO IDEL = 1, NBAST                  !loop over all deltas

             ISYDEL  = ISAO(IDEL)
             ISY0IAJ = MULD2H(ISYDEL,ISYM0)
             !length of this portion of integral vect
             LEN0    = NT2BCD(ISY0IAJ)         
             IADR0   = IT2DEL0(IDEL)

C           ------------------------------------------------------
C           calculate RA intermediates from zeroth-order integrals:
C           do only for a new T^A?
C           ------------------------------------------------------
             CALL GETWA2(LU0IAJB,FN0IAJB,WORK(KXIAJB),IADR0,LEN0)
             CALL CC_RIM(WORK(KRAIM(ITRAN)),WORK(KTA2AM),ISYMTA,
     &                  WORK(KXIAJB),ISY0IAJ,IDEL,ISYDEL)

c            XNORM = DDOT(NEMAT1(ISYMTA),WORK(KRAIM(ITRAN)),1,
c     &                                 WORK(KRAIM(ITRAN)),1)
c             WRITE(LUPRI,*) 'FBTA: Check RA intermediate', XNORM
c

             ISY1IAJ = MULD2H(ISYDEL,ISYHOP)
             LEN1    = NT2BCD(ISY1IAJ)
             IADRB   = IT2DELB(IDEL,ITRAN) 

C           ----------------------------------------------------
C           calculate RA-bar intermediates from derivative integrals:
C           do always, integrals can change at each ITRAN (IOPERB)
C           ----------------------------------------------------
             CALL GETWA2(LU1IAJB,FN1IAJB,WORK(KXIAJB),IADRB,LEN1)
             CALL CC_RIM(WORK(KRQAIM(ITRAN)),WORK(KTA2AM),ISYMTA,
     &                        WORK(KXIAJB),ISY1IAJ,IDEL,ISYDEL)

           END DO

c          WRITE (LUPRI,*) 'FBTA> test RA and RQA intemediates'
c          LEN   = NEMAT1(ISYMTA)
c          XNORM = DDOT(LEN,WORK(KRAIM(ITRAN)),1,WORK(KRAIM(ITRAN)),1)
c          WRITE (LUPRI,*) 'Norm of saved RAIM       is ',XNORM
c          LEN   = NEMAT1(ISYHTA)
c          XNORM = DDOT(LEN,WORK(KRQAIM(ITRAN)),1,WORK(KRQAIM(ITRAN)),1)
c          WRITE (LUPRI,*) 'Norm of saved RQAIM       is ',XNORM
*
* Just for testing with YBARA: contract with Lambda0
*
*          CALL CC_EIM(WORK(KEND1),DUMMY,WORK(KRAIM(ITRAN)), DUMMY,
*     &                 DUMMY,DUMMY,DUMMY,DUMMY,WORK(KLAMH0),
*     &                 DUMMY,ISYM0,DUMMY,IDUM,DUMMY,.FALSE.,
*     &                 .TRUE.,.FALSE.,
*     &                 LRELAX,1,ISYMTA)
*          WRITE (LUPRI,*) 'FBTA> test transformed RA'
*          LEN   = NMATAB(ISYMTA)
*          XNORM = DDOT(LEN,WORK(KEND1),1,WORK(KEND1),1)
*          WRITE (LUPRI,*) 'Norm of saved RAIM  (cb)     is ',XNORM
*

         END DO

      END IF ! (LLTWOEL .OR. LLRELAX) 

*---------------------------------------------------------------------*
* start a new loop over the transformations and calculate now the
* complete vectors from the intermediates:
*---------------------------------------------------------------------*
* read in memory those vectors that do not depend on ITRAN:
* the t0_2 amplitudes and the I_ijkdelta (zero order) integrals
* check whether I need them at all (2 e- part)
* The T0_1 are already in memory
*---------------------------------------------------------------------*
* if ccsd maybe
*
      KT2AMP0 = KEND1
      KEND1   = KT2AMP0 + NT2AM(ISYM0)

      LWRK1 = LWORK - KEND1
      IF (LWRK1.LE.0) THEN
        WRITE (LUPRI,*) 
     *         'Insuf. WRK in CC_FBTA: avail',LWRK1,' need ', KEND1
        CALL QUIT('Insufficient wrk in CC_FBTA: T20')
      END IF
      
      IF ((.NOT.CCS).AND.(LLTWOEL.OR.LLRELAX)) THEN
        IOPT = 2
        CALL CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KDUM),WORK(KT2AMP0))           !read in T(0)_2
      END IF
*
* Read in memory the COMPLETE (N3ODEL(ISYM0)) array IJKdelta(0) 
*
      IF (LLTWOEL.OR.LLRELAX) THEN              
        KXIJK0 = KEND1
        KEND1  = KXIJK0 + N3ODEL(ISYM0)
        LWRK1 = LWORK - KEND1
        IF (LWRK1.LE.0) CALL QUIT('Insufficient wrk in CC_FBTA (XIJK0)')
        KOFF   = KXIJK0
        DO ISYMD = 1, NSYM
           ISYIJK = ISYMD
           LENIJK = NMAIJK(ISYIJK)*NBAS(ISYMD)
           IOFF = I3ODEL(ISYIJK,ISYMD) + 1
           IF (LENIJK.GT.0) THEN
             CALL GETWA2(LU0IJKL,FN0IJKL,WORK(KOFF),IOFF,LENIJK)
             IF (LOCDBG) THEN
               XNORM = DNRM2(LENIJK,WORK(KOFF),1)
               WRITE (LUPRI,*)
     &            'CC_FBTA> ISYMD:',ISYMD, ' Norm XIJKdel0:',XNORM
             END IF
           END IF
           KOFF = KOFF + LENIJK
        END DO
      END IF
*
*---------------------------------------------------------------------*
      KEND1SV = KEND1
      LWRK1SV = LWRK1
*---------------------------------------------------------------------*
      DO ITRAN = ISTART(IBATCH), IEND(IBATCH)
          
        IOPERB = IFBTRAN(1,ITRAN)  ! operator index
        IDLSTL = IFBTRAN(2,ITRAN)  ! ZETA vector index
        IDLSTR = IFBTRAN(3,ITRAN)  ! T^A vector index
        IVEC   = IFBTRAN(4,ITRAN)  ! output index
        IRELAX = IFBTRAN(5,ITRAN)  ! flag for relax. contrib.

        ISYHOP = ISYOPR(IOPERB)     ! symmetry of hamiltonian
        LABELH = LBLOPR(IOPERB)     ! operator label
        LTWOEL = LPDBSOP(IOPERB)    ! two-electron contribution
        LRELAX = LTWOEL .OR. (IRELAX.GE.1)

        ISYCTR = ILSTSYM(LISTL,IDLSTL)  ! symmetry of ZETA vector
        ISYMTA = ILSTSYM(LISTR,IDLSTR)  ! symmetry of T^A vector
        ISYZTA = MULD2H(ISYMTA,ISYCTR)  ! symmetry of Z*TA 
        ISYHTA = MULD2H(ISYMTA,ISYHOP)  ! symmetry of Op*TA 
        ISYRES = MULD2H(ISYHOP,ISYZTA)  ! symmetry of result vector

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_FBTA> IOPERB,IVEC   = ',IOPERB,IVEC
          WRITE (LUPRI,*) 'CC_FBTA> LTWOEL,LRELAX = ',LTWOEL,LRELAX
          WRITE (LUPRI,*) 'CC_FBTA> LABELH,ISYHOP = ',LABELH,ISYHOP
          WRITE (LUPRI,*) 'CC_FBTA> LISTL,IDLSTL  = ',LISTL,IDLSTL
          WRITE (LUPRI,*) 'WORK(0) = ',WORK(ITST)
          CALL FLSHFO(LUPRI)
        END IF

        KEND1 = KEND1SV
        LWRK1 = LWRK1SV
*
*=====================================================================*
* calculate now the F^B T^A vector:
* Note that this part is entered both for usual one-el unrelaxed 
* operators and for two-el/relaxed operators
*=====================================================================*
* 
        KZETA1 = KEND1
        KRES1  = KZETA1 + NT1AM(ISYCTR)
        KEND2  = KRES1  + NT1AM(ISYRES)
        IF (CC2 .OR. CCSD) THEN
           KRES2 = KEND2
           KEND2 = KRES2  + NT2AM(ISYRES)
        END IF
        LEN = NT1AM(ISYRES)+NT2AM(ISYRES)
        LWRK2 = LWORK - KEND2
*
        IF (LWRK2 .LT. 0) THEN
          WRITE (LUPRI,*) ' Available:', LWRK2, ' Need', KEND2
          CALL QUIT('Insufficient work space in CC_FBTA. (CCFBTAF)')
          CALL FLSHFO(LUPRI)
        END IF
*
        IOPT = 1
        CALL CC_RDRSP(LISTL,IDLSTL,ISYCTR,IOPT,MODEL,
     &                              WORK(KZETA1),WORK(KDUM))
*
        CALL CCFBTAF(WORK(KRES1),WORK(KRES2),
     &              WORK(KRAIM(ITRAN)),WORK(KRQAIM(ITRAN)),
     &              WORK(KFAIM(ITRAN)),WORK(KFQAIM(ITRAN)),
     &              WORK(KXAINT(ITRAN)),WORK(KYAINT(ITRAN)),
     &              WORK(KZETA1),
     &              WORK(KT1AMP0),WORK(KT2AMP0),
     &              WORK(KXIJK0),
     &              WORK(KONEH0),WORK(KONEHB(ITRAN)),
     &              WORK(KFOCK0),WORK(KFOCKBCC(ITRAN)),
     &              WORK(KFOCKACC(ITRAN)),WORK(KFCKBACC(ITRAN)),
     &              LUBFA, FILBFA, IADRBFA(1,ITRAN),
     &              LUBFQA,FILBFQA,IADRBFQA(1,ITRAN),
     &              LUBFZA, FILBFZA, IADRBZA(1,ITRAN),
     &              LUBFZQA,FILBFZQA,IADRBZQA(1,ITRAN),
     &              LU0IAJB, FN0IAJB, IT2DEL0,
     &              LU1IAJB, FN1IAJB, IT2DELB(1,ITRAN),
     &              LU0IJCB, FN0IJCB, 
     &              LU0CJIB, FN0CJIB, IT2DELA(1,ITRAN),
     &              LU1IJCB, FN1IJCB,
     &              LU1CJIB, FN1CJIB, IT2DELBA(1,ITRAN),
     &              LU1IJKL, FN1IJKL,
     &              LUPQMO, FILPQMO,IADRPQMO(1,ITRAN),
     &              LUPQAMO, FILPQAMO,IADPQAMO(1,ITRAN),
     &              IADRI0,                !to append IAJB(1) to IAJ,delta(0)
     &              WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KLAMDPQ(ITRAN)),WORK(KLAMDHQ(ITRAN)),
     &              WORK(KLAMDPA(ITRAN)),WORK(KLAMDHA(ITRAN)),
     &              WORK(KLAMPQA(ITRAN)),WORK(KLAMHQA(ITRAN)),
     &              LISTL,IDLSTL,LABELH,ISYHOP,
     &              LISTR,IDLSTR,
     &              LTWOEL,LRELAX,ITRAN,WORK(KEND2), LWRK2)
*
*---------------------------------------------------------------------*
* transformation finished: write result vector to file:
*---------------------------------------------------------------------*
*
        IF (CCS) THEN
           MODELW = 'CCS       '
           IOPTW  = 1
        ELSE IF (CC2) THEN
           MODELW = 'CC2       '
           IOPTW  = 3
        ELSE IF (CCSD) THEN
           MODELW = 'CCSD      '
           IOPTW  = 3
        ELSE
           CALL QUIT('Unknown coupled cluster model in CC_FBTA.')
        END IF
*
        IF (IOPTRES.EQ.0 .OR. IOPTRES.EQ.1) THEN
           IFBTRAN(4,ITRAN) = IADR_FBTA
           CALL PUTWA2(LUFBTA,FILFBTA,WORK(KRES1),IADR_FBTA,
     &                                          NT1AM(ISYRES))
           IADR_FBTA = IADR_FBTA + NT1AM(ISYRES)
           IF (.NOT.CCS) THEN
              CALL PUTWA2(LUFBTA,FILFBTA,WORK(KRES2),IADR_FBTA,
     &                    NT2AM(ISYRES))
              IADR_FBTA = IADR_FBTA + NT2AM(ISYRES)
           END IF
        ELSE IF (IOPTRES.EQ.3) THEN
           IFILE  = IFBTRAN(4,ITRAN)
           CALL CC_WRRSP(FILFBTA,IFILE,ISYRES,IOPTW,MODELW,WORK(KDUM),
     &                   WORK(KRES1),WORK(KRES2),WORK(KEND2),LWRK2)
        ELSE IF (IOPTRES.EQ.4) THEN
           IFILE  = IFBTRAN(4,ITRAN)
           CALL CC_WARSP(FILFBTA,IFILE,ISYRES,IOPTW,MODELW,WORK(KDUM),
     &                   WORK(KRES1),WORK(KRES2),WORK(KEND2),LWRK2)
        ELSE IF (IOPTRES.EQ.5) THEN
           IF (.NOT.CCS) CALL CCLR_DIASCL(WORK(KRES2),TWO,ISYRES)
           CALL CCDOTRSP(IFBDOTS,FBCONS,IOPTW,FILFBTA,ITRAN,NFBTRAN,
     &                   MXVEC,WORK(KRES1),WORK(KRES1),ISYRES,
     &                   WORK(KEND2),LWRK2)
        ELSE
           CALL QUIT('Illegal value for IOPTRES in CC_FBTA.')
        END IF
     
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'Final result of CC_FBTA:'
          WRITE (LUPRI,*) 'operator label:      ',LBLOPR(IOPERB)
          WRITE (LUPRI,*) 'output file type:    ',FILFBTA
          WRITE (LUPRI,*) 'index of output file:',IFILE
          WRITE (LUPRI,*) 'two-electron contr.: ',LTWOEL
          WRITE (LUPRI,*) 'relax/reorth. contr.:',LRELAX
          CALL CC_PRP(WORK(KRES1),WORK(KRES2),ISYRES,1,1)
          XNORM = DDOT(LEN, WORK(KRES1),1, WORK(KRES1),1)
          CALL FLSHFO(LUPRI)
        END IF
*
**=====================================================================*
** end of loop over transformations.
**=====================================================================*
*
      END DO ! ITRAN
      END DO ! IBATCH 
*
*---------------------------------------------------------------------*
* close and delete integral & density files:
*---------------------------------------------------------------------*
      CALL WCLOSE2(LU0IAJB, FN0IAJB,  'DELETE')
      CALL WCLOSE2(LU0IJCB, FN0IJCB,  'DELETE')
      CALL WCLOSE2(LU0CJIB, FN0CJIB,  'DELETE')
      CALL WCLOSE2(LU1IAJB, FN1IAJB,  'DELETE')
      CALL WCLOSE2(LU1IJCB, FN1IJCB,  'DELETE')
      CALL WCLOSE2(LU1CJIB, FN1CJIB,  'DELETE')
      CALL WCLOSE2(LU0IJKL, FN0IJKL,  'DELETE')
      CALL WCLOSE2(LU1IJKL, FN1IJKL,  'DELETE')
*
      CALL WCLOSE2(LUPQMO,  FILPQMO,  'DELETE')
      CALL WCLOSE2(LUPQ0,   FILPQ0,   'DELETE')
      CALL WCLOSE2(LUPQ1,   FILPQ1,   'DELETE')
      CALL WCLOSE2(LUPQAMO, FILPQAMO, 'DELETE')
      CALL WCLOSE2(LUPQA0,  FILPQA0,  'DELETE')
      CALL WCLOSE2(LUPQA1,  FILPQA1,  'DELETE')
*
      CALL WCLOSE2(LUBFA,   FILBFA,   'DELETE')
      CALL WCLOSE2(LUBFQA,  FILBFQA,  'DELETE')
*
      CALL WCLOSE2(LUBFZA,  FILBFZA,  'DELETE')
      CALL WCLOSE2(LUBFZQA, FILBFZQA, 'DELETE')
*
*---------------------------------------------------------------------*
* if IOPTRES=1 and enough work space available, read result
* vectors back into memory:
*---------------------------------------------------------------------*

* check size of work space:
      IF (IOPTRES .EQ. 1) THEN
        LENALL = IADR_FBTA - 1
        IF (LENALL .GT. LWORK) IOPTRES = 0
      END IF                                  
* read the result vectors back into memory:
      IF (IOPTRES .EQ. 1) THEN

        CALL GETWA2(LUFBTA,FILFBTA,WORK(1),1,LENALL)

        IF (LOCDBG) THEN
          DO ITRAN = 1, NFBTRAN
            IF (ITRAN.LT.NFBTRAN) THEN
              LEN     = IFBTRAN(4,ITRAN+1)-IFBTRAN(4,ITRAN)
            ELSE
              LEN     = IADR_FBTA-IFBTRAN(4,NFBTRAN)
            END IF
            KRES1 = IFBTRAN(4,ITRAN)
            XNORM = DDOT(LEN, WORK(KRES1),1, WORK(KRES1),1)
            WRITE (LUPRI,*) 'Read FB matrix transformation nb. ',ITRAN
            WRITE (LUPRI,*)
     &         'Adress,length,NORM:',IFBTRAN(4,ITRAN),LEN,XNORM
          END DO
          CALL FLSHFO(LUPRI)
        END IF
      END IF                                        
      IF (IOPTRES.EQ.0 ) THEN
        CALL WCLOSE2(LUFBTA, FILFBTA, 'KEEP')
      ELSE IF (IOPTRES.EQ.1) THEN
        CALL WCLOSE2(LUFBTA, FILFBTA, 'DELETE')
      ELSE IF (IOPTRES.EQ.3 .OR. IOPTRES.EQ.4 .OR. IOPTRES.EQ.5) THEN
        CONTINUE
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_FBTA.')
      END IF                                         
*=======================================================================*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CC_FBTA ended successfully (?) ... return now.'
        WRITE (LUPRI,*) 'IT IS A MIRACLE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        CALL FLSHFO(LUPRI)
      END IF

      RETURN
      END
*=====================================================================*
*                  END OF SUBROUTINE CC_FBTA                          *
*=====================================================================*IM
