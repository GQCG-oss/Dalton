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
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE T_ROW_TO_H(T,H,K,TKK)
*
* Set H integrals
*
*    Column K : H(P,K) = T(P,K)/T(K,K), P.NE.K
*    Other Columns     = 0
* - and return T_{kk} in TKK
*
*
* Jeppe Olsen, Jan 98
* For rotation of CI vectors
*
      IMPLICIT REAL*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "glbbas.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "lucinp.inc"
*. Input ( in blocked form)
      DIMENSION T(*)
*. Output ( also in blocked form)
      DIMENSION H(*)
*
      KSM = ISMFSO(K)
      KOFF = IBSO(KSM)
      KREL = K - KOFF + 1
      NK = NTOOBS(KSM)

*
      ZERO = 0.0D0
      CALL SETVEC(H,ZERO,NTOOB**2)
*
      IOFF = IFRMR(WORK(KPGINT1A(1)),1,KSM)
      CALL COPVEC(T(IOFF+(KREL-1)*NK),H(IOFF+(KREL-1)*NK),NK)
      TKK = H(IOFF-1+(KREL-1)*NK+KREL)
      IF(TKK .NE. 0.0D0) THEN
        FAC = 1.0D0/TKK
        CALL SCALVE(H(IOFF+(KREL-1)*NK),FAC,NK)
C       H(IOFF-1+(K-1)*NK+K) = H(IOFF-1+(K-1)*NK+K) -1.0D0
        H(IOFF-1+(KREL-1)*NK+KREL) = 0.0D0
      ELSE
C       TKK = 1.0D0
        TKK = 0.0D0
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output from T_ROW_TO_H '
        WRITE(6,*) ' Row to be transferred ', KREL
        WRITE(6,*) ' Updated H matrix '
        WRITE(6,*) ' NTOOB, NTOOBS(1)', NTOOB, NTOOBS(1)
        CALL APRBLM2(H,NTOOBS,NTOOBS,NSMOB,0)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE T_TO_NK_VEC(T,KORB,ISM,ISPC,LUCIN,LUCOUT,C)
*
* Evaluate T**(NK_operator) times vector on file LUIN
* to yield vector on file LUOUT
* (NK_operator is number operator for orbital K )
*
* Note LUCIN and LUCOUT are both rewinded before read/write
* Input
* =====
*  T : Input constant
*  KORB : Orbital in symmetry order
*
*  ISM,ISPC : Symmetry and space of state on LUIN
*  C : Scratch block
*
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLASTR, KLBSTR, KLKAOC, KLKBOC,
     &          KLCLBT, KLCLEBT, KLCI1BT, KLCIBT, KLCBLTP
!               for addressing of WORK
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strinp.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "crun.inc"
#include "csm.inc"

*. Scratch block, must hold a batch of blocks
      DIMENSION C(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' T_TO_NK_VEC speaking '
        WRITE(6,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'T_TO_N')
*. Set up block and batch structure of vector
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      CALL Z_BLKFO(ISPC,ISM,IATP,IBTP,KLCLBT,KLCLEBT,
     &            KLCI1BT,KLCIBT,KLCBLTP,NBATCH,NBLOCK)
C           Z_BLKFO(ISPC,ISM,IATP,IBTP,KPCLBT,KPCLEBT,
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
      CALL MEMMAN(KLKAOC,MXNSTR,     'ADDL  ',1,'KLKAOC')
      CALL MEMMAN(KLKBOC,MXNSTR,     'ADDL  ',1,'KLKBOC')
      
*. Orbital K in type ordering
      KKORB = IREOST(KORB)
      CALL T_TO_NK_VECS
     &  (T,KKORB,C,LUCIN,LUCOUT,WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &      NBLOCK,WORK(KLCIBT),
     &      NAEL,NBEL,WORK(KLASTR),WORK(KLBSTR),
     &      WORK(KLCBLTP),NSMST,
     &      ICISTR,NTOOB,WORK(KLKAOC),WORK(KLKBOC))

      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'T_TO_N')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE T_TO_NK_VECS(T,KORB,C,LUCIN,LUCOUT,NSSOA,NSSOB,
     &                 NBLOCK,IBLOCK,
     &                 NAEL,NBEL,IASTR,IBSTR,IBLTP,NSMST,
     &                 ICISTR,NORB,IKAOCC,IKBOCC)
*
* Multiply Vector in LUCIN with t **NK_op to yield vector on LUCOUT
*
* Both files are initially rewinded
*
*
* Jeppe Olsen, Feb. 1998
*

      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
*. Scratch
      DIMENSION C(*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IKAOCC(*),IKBOCC(*)
*. Specific input
      DIMENSION IBLOCK(8,NBLOCK)
      DIMENSION IBLTP(*)
*
      CALL REWINE(LUCIN,-1)
      CALL REWINE(LUCOUT,-1)
*
      T2 = T**2
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
C?      WRITE(6,*) ' IATP IBTP IASM IBSM ', IATP,IBTP,IASM,IBSM
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
*. Occupation of orb KORB
        DO JSTR = 1, NASTR1
          KOCC = 0
          DO JAEL = 1, NAEL
            IF(IASTR(JAEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKAOCC(JSTR) = KOCC
        END DO
C?      WRITE(6,*) ' IKAOCC array '
C?      CALL IWRTMA(IKAOCC,1,NASTR1,1,NASTR1)


*. Obtain Beta  strings of sym IBSM and type IBTP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                           NORB,0,IDUM,IDUM)
C?      WRITE(6,*) ' After GETSTR, NBSTR1=',NBSTR1
*. Occupation of orb KORB
        DO JSTR = 1, NBSTR1
C?        write(6,*) ' JSTR = ', JSTR
          KOCC = 0
          DO JBEL = 1, NBEL
C?          write(6,*) JBEL, IBSTR(JBEL,JSTR)
            IF(IBSTR(JBEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKBOCC(JSTR) = KOCC
        END DO
C?      WRITE(6,*) ' IKBOCC array '
C?      CALL IWRTMA(IKBOCC,1,NBSTR1,1,NBSTR1)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
C?      WRITE(6,*) ' IBLTP ', IBLTP(IASM)
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
C?      WRITE(6,*) ' NIA NIB ', NIA,NIB
*
        IMZERO = 0
        IF( ICISTR.GE.2 ) THEN
*. Read in a Type-Type-symmetry block
          CALL IFRMDS(LDET,1,-1,LUCIN)
          CALL FRMDSC_LUCI(C,LDET,-1,LUCIN,IMZERO,IAMPACK)
        END IF
        IF(IMZERO.NE.1) THEN
*
          IDET = 0
          DO  IB = 1,NIB
            IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
              MINIA = IB
            ELSE
              MINIA = 1
            END IF
            DO  IA = MINIA,NIA
*
              IDET = IDET + 1
C?            WRITE(6,*) ' IA IB IDET',IA,IB,IDET
              KABOCC = IKAOCC(IA)+IKBOCC(IB)
              IF(KABOCC.EQ.1) THEN
                C(IDET) = T*C(IDET)
              ELSE IF(KABOCC.EQ.2) THEN
                C(IDET) = T2 *C(IDET)
              END IF
            END DO
*           ^ End of loop over alpha strings
          END DO
*         ^ End of loop over beta strings
*
        END IF
*       ^ End of if statement for nonvanishing blocks
*. Save result on LUCOUT
        CALL ITODS(LDET,1,-1,LUCOUT)
        CALL TODSC_LUCI(C,LDET,-1,LUCOUT)
      END DO
*     ^ End of loop over blocks
*. This is the end, the end of every file my friend, the end
       CALL ITODS(-1,1,-1,LUCOUT)

      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRA1_SIMPLE(C)
*
* Transform one-electron integrals
*
* Integrals from KINT1O used : Integrals without particle-hole corrections
*
* Input integrals in KINT1
* Output integrals also in KINT1
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLH1, KLSCR
!               for addressing of WORK
*. General Input
#include "mxpdim.inc"
#include "orbinp.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*. Specific Input
      DIMENSION C(*)
*
      NTEST = 100
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA1_S')
*. Largest symmetry block of orbitals
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Largest number of orbitals in sym block ', MXSOB
      END IF
*. Space for 1e- integrals and two sym blocks
      LENH1 = NTOOB ** 2
      CALL MEMMAN(KLH1,LENH1,'ADDL  ',2,'LH1   ')
      
      LENSCR = 2 * MXSOB ** 2
      CALL MEMMAN(KLSCR,LENSCR,'ADDL   ',2,'H1SCR ')
      WRITE(6,*) ' LENH1 LENSCR ', LENH1, LENSCR
*. and do it
C     TRAN_SYM_BLOC_MAT3 (AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
      WRITE(6,*) ' Before TRAN_SYM '
      CALL COPVEC(WORK(KINT1O),WORK(KINT1),NINT1)
      CALL TRAN_SYM_BLOC_MAT3(WORK(KINT1),C,NSMOB,NTOOBS,NTOOBS,
     &                        WORK(KLH1),WORK(KLSCR),1)
      WRITE(6,*) ' After TRAN_SYM '
*. Copy to original array
      CALL COPVEC(WORK(KLH1),WORK(KINT1),NINT1)
*. Flush memory
      WRITE(6,*) ' Returning from TRA1 '
      CALL MEMMAN(IDUM,IDUM,'FLUSM  ',IDUM,'TRA1_S')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRA2_SIMPLE(C)
*
* Trivial 2-electron integral transformation
*
*. Input integrals in KINT2
*. Output integrals in KINT2
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLBL1, KLBL2
!               for addressing of WORK
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "multd2h.inc"
*. Specific input
      DIMENSION C(*)
*
      NTEST  = 0
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRA2_S')
*. Largest symmetry block
      MXSOB = IMNMX(NTOOBS,NSMOB,2)
*. Two symmetry blocks
      LENBL = MXSOB ** 4
      IF(NTEST.GE.100)
     &WRITE(6,*) ' Size of symmetry block ', LENBL
      CALL  MEMMAN(KLBL1,LENBL,'ADDL ',2,'KLBL1 ')
      CALL  MEMMAN(KLBL2,LENBL,'ADDL ',2,'KLBL2 ')
C     WRITE(6,*) ' Two integral blocks allocated '
*. Loop over symmetry blocks
      DO ISM = 1, NSMOB
        DO JSM = 1, ISM
          DO KSM = 1, ISM
            IF(KSM.EQ.ISM) THEN
              LSM_MX = JSM
            ELSE
              LSM_MX = KSM
            END IF
            DO LSM = 1, LSM_MX
*. Ensure that integrals have correct symmetry
              INTSYM = 1
              IJSM = MULTD2H(ISM,JSM)
              IJKSM = MULTD2H(IJSM,KSM)
              IF(INTSYM.EQ.MULTD2H(IJKSM,LSM)) THEN
*. Fetch 2-e integral block (IJ!KL)
C             LGETINT(XINT,ITP,ISM,JTP,JSM,KTP,KSM,LTP,LSM,
C    &               IXCHNG,IKSM,JLSM,ICOUL)
C             WRITE(6,*)  'ISM, JSM, KSM, LSM ',ISM,KSM,KSM,LSM
              CALL LGETINT(WORK(KLBL1),0,ISM,0,JSM,0,KSM,0,LSM,
     &                     0,0,0,1)
*. (Type = 0  => Complete symmetryblock)
*. Offsets and dimensions for symmetryblocks in C
              IOFF = 1
              DO IISM = 1, ISM-1
                IOFF = IOFF + NTOOBS(IISM)**2
              END DO
              NI = NTOOBS(ISM)
*
              JOFF = 1
              DO JJSM = 1, JSM-1
                JOFF = JOFF + NTOOBS(JJSM)**2
              END DO
              NJ = NTOOBS(JSM)
*
              KOFF = 1
              DO KKSM = 1, KSM-1
                KOFF = KOFF + NTOOBS(KKSM)**2
              END DO
              NK = NTOOBS(KSM)
*
              LOFF = 1
              DO LLSM = 1, LSM-1
                LOFF = LOFF + NTOOBS(LLSM)**2
              END DO
              NL = NTOOBS(LSM)
*. Transform 2-electron integral block
              CALL TRA_2EL_BLK_SIMPLE(WORK(KLBL1),
     &        C(IOFF),NI,C(JOFF),NJ,C(KOFF),NK,C(LOFF),NL,WORK(KLBL2))
*. Transfer symmetry block to integral list
              CALL PUTINT(WORK(KLBL1),0,ISM,0,JSM,0,KSM,0,LSM)
            END IF
*           ^ Check if integrals have correct symmetry
            END DO
          END DO
        END DO
      END DO
*     ^ End of loop over symmetries
*. Dump integrals to file LU90
      LU90 = 90
      CALL DMPINT(LU90)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRA2_S')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRA_2EL_BLK_SIMPLE(XINT,CI,NI,CJ,NJ,CK,NK,CL,NL,SCR)
*
* Transform 2-electron integral block
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CI(NI,NI),CJ(NJ,NJ),CK(NK,NK),CL(NL,NL)
*. Input and output
      DIMENSION XINT(*)
* Matrix given in complete form XI(NI,NJ,NK,NL)
*. Scratch
      DIMENSION SCR(*)
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input integral block as X(IJ,KL)'
        CALL WRTMT_LU(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*. Transform first two indeces
      DO K = 1, NK
       DO L = 1, NL
         KLOF = 1 + ((L-1)*NK+K-1)*NI*NJ
         FACTORC = 0.0D0
         FACTORAB = 1.0D0
         CALL MATML7(SCR(KLOF),XINT(KLOF),CJ,NI,NJ,NI,NJ,NJ,NJ,
     &                FACTORC,FACTORAB,0)
         CALL MATML7(XINT(KLOF),CI,SCR(KLOF),NI,NJ,NI,NI,NI,NJ,
     &                FACTORC,FACTORAB,1)
       END DO
      END DO
*. Transpose Block X(IJ,KL)
      NIJ = NI*NJ
      NKL = NK*NL
      CALL TRPMT3(XINT,NIJ,NKL,SCR)
*. Matrix is now SCR(KL,IJ)
*. Transform last two indeces
      DO I = 1, NI
        DO J = 1, NJ
         IJOF = ((J-1)*NI+I-1)*NK*NL+1
         CALL MATML7(XINT(IJOF),SCR(IJOF),CL,NK,NL,NK,NL,NL,NL,
     &               FACTORC,FACTORAB,0)
         CALL MATML7(SCR(IJOF),CK,XINT(IJOF),NK,NL,NK,NK,NK,NL,
     &               FACTORC,FACTORAB,1)
        END DO
      END DO
*. Transpose to obtain XINT(IJ,KL)
      CALL TRPMT3(SCR,NKL,NIJ,XINT)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output integral block as X(IJ,KL)'
        CALL WRTMT_LU(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRACI_LUCI(X,LUCIN,LUCOUT,IXSPC,IXSM,VEC1,VEC2)
*
* A rotation matrix X is defining expansion from
* old to new orbitals
*        PHI(NEW) = PHI(OLD) * X
*
* change CI coefficients(sym IXSM, space IXSPC )
* so they corresponds to PHI(NEW) basis
*
* The input CI vector is on LUCIN and the transformed CI vector
* will be delivered on LUCOUT.
*
* Transformation as conceived by Per-Aake Malmquist
* (I.J.Q.C. vol XXX, p479 ,1986 (OCTOBER ISSUE ))
*
*  Jeppe Olsen 1988
*
* New LUCIA version of Jan 1998
*
* note The transformation matrix X is supposed to be in complete form
* as a matrix over NTOOB orbitals.
*
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*8 KLSCR, KLT, KLTB
!               for addressing of WORK
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "clunit.inc"
*. Common block for communicating with sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
      DIMENSION X(*)
*
      CALL QENTER('TRACI')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACI ')
*
      NTEST = 10
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' Welcome to TRACI '
        WRITE(6,*) ' ================'
        WRITE(6,*)
        WRITE(6,*) ' IXSPC,IXSM = ', IXSPC,IXSM
      END IF
*. Memory allocation
* for a matrix T
      CALL MEMMAN(KLT,NTOOB**2,'ADDL  ',2,'TMAT  ')
      CALL MEMMAN(KLTB,NTOOB**2,'ADDL  ',2,'TMATBL')
*. Scratch in PAMTMT
      LSCR = NTOOB**2 +NTOOB*(NTOOB+1)/2
      CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'KLSCR ')
      
*. Obtain T matrix used for transformation, for each symmetry separately
      DO ISM = 1, NSMOB
        IF(ISM.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISM-1)**2
        END IF
        IF(NTOOBS(ISM).GT.0)
     &  CALL PAMTMT(X(IOFF),WORK(KLT-1+IOFF),WORK(KLSCR),NTOOBS(ISM))
      END DO
*. Transform CI-vector
      ICSPC = IXSPC
      ICSM  = ICSM
      ISSPC = IXSPC
      ISSM  = IXSM
*
      CALL TRACID(WORK(KLT),LUCIN,LUCOUT,LUSC1,LUSC2,LUSC3,
     &            VEC1,VEC2)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACI ')
      CALL QEXIT('TRACI')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRACI_CTL
*
* Master routine for transforming CI vectors to new orbital basis
*
* Input vectors on LUC
* output vectors delivered on LUHC
* (LUSC1,LUSC2,LUSC3) used as scratch files
*
* Jeppe Olsen, January 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KVEC1, KVEC2, KLCMOMO, KLH1SAVE,
     &          KLCLBT, KLCLEBT, KLCI1BT, KLCIBT, KLCBLTP
!               for addressing of WORK
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "clunit.inc"
#include "crun.inc"
#include "cstate.inc"
#include "glbbas.inc"
#include "orbinp.inc"
*
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
      NTEST = 0
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACIC')
*
      WRITE(6,*)
      WRITE(6,*) ' Wellcome to TRACI_CTL'
      WRITE(6,*) ' ====================='
      WRITE(6,*)
*. The three scratch  blocks
C          GET_3BLKS(KVEC1,KVEC2,KC2)
      CALL GET_3BLKS(KVEC1,KVEC2,KVEC3)

*. MO-MO transformation matrix :
      CALL MEMMAN(KLCMOMO,NTOOB**2,'ADDL  ',2,'CMOMO ')
*. Copy of one-electron integrals
      CALL MEMMAN(KLH1SAVE,NTOOB**2,'ADDL  ',2,'H1SAVE')
*. We are going to mess with the one-electron integrals, take a copy
      CALL COPVEC(WORK(KINT1),WORK(KLH1SAVE),NTOOB*NTOOB)
*. Set up block structure of CI space
      IATP = 1
      IBTP = 2
      CALL  Z_BLKFO(ISSPC,ISSM,IATP,IBTP,KLCLBT,KLCLEBT,
     &      KLCI1BT,KLCIBT,KLCBLTP,NBATCH,NBLOCK)
C           Z_BLKFO(ISPC,ISM,IATP,IBTP,KPCLBT,KPCLEBT,
C    &              KPCI1BT,KPCIBT,KPCBLTP,NBATCH,NBLOCK)
*
      LBLK = -1
      CALL REWINE(LUC,-1)
*. Make sure LUDIA corresponds to original def
      LUDIA = 20
      CALL REWINE(LUDIA,-1)
      DO JROOT = 1, NROOT
*. One-electron density for root JROOT
        CALL REWINE(LUSC1,-1)
        CALL COPVCD(LUC,LUSC1,WORK(KVEC1),0,LBLK)
        CALL COPVCD(LUSC1,LUSC2,WORK(KVEC1),1,LBLK)
        CALL DENSI2(1,WORK(KRHO1),WORK(KRHO2),
     &        WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC2,EXPS2)
*. Obtain MO-MO transformation matrix
C             GET_CMOMO(CMOMO)
        CALL GET_CMOMO(WORK(KLCMOMO))
        WRITE(6,*) ' Memcheck after GET_CMOMO'
        CALL LMEMCHK('TRACI_CTL')
*
*. Transform CI vector : Input on LUHC, output on LUDIA (!)
        CALL COPVCD(LUSC1,LUHC,WORK(KVEC1),1,LBLK)
*
        CALL TRACI_LUCI(WORK(KLCMOMO),LUHC,LUDIA,ISSPC,ISSM,
     &             WORK(KVEC1),WORK(KVEC2))
C            TRACI(X,LUCIN,LUCOUT,IXSPC,IXSM,VEC1,VEC2)
       WRITE(6,*)
       WRITE(6,*) ' Analysis of rotated state number ', JROOT
       WRITE(6,*) ' ====================================='
       WRITE(6,*)
       CALL GASANA(WORK(KVEC1),NBLOCK,WORK(KLCIBT),WORK(KLCBLTP),
     &                LUSC1,ICISTR)
      END DO
*     ^ End of loop over roots
      CALL REWINE(LUDIA,-1)
      IF(NTEST.GE.100) THEN
        DO JROOT = 1, NROOT
          CALL WRTVCD(WORK(KVEC1),LUDIA,0,LBLK)
        END DO
      END IF
*
*. clean up time : copy 1-e integrals back in place
      CALL COPVEC(WORK(KLH1SAVE),WORK(KINT1),NTOOB*NTOOB)
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACIC')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRACID(T,LUCIN,LUCOUT,LUSC1,LUSC2,LUSC3,
     &           VEC1,VEC2)
*
* Transform CI vector on LUCIN with T matrix after
* Docent Malmquist's recipe. Place result as next vector on LUOUT
*
* The transformation is done as a sequence of one-electron transformations
*
* with each orbital transformation being
*
* Sum(k=0,2) ( 1/k! sum(n'.ne.n) S(n'n) E_{n'n} ) Tnn^N_n
*
* with Sn'n = T(n'n)/Tnn
*
* each transformation is
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "oper.inc"
#include "intform.inc"
#include "lucinp.inc"
#include "orbinp.inc"
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
      REAL*8 INPRDD
*. Input
      DIMENSION T(*)
*. Scratch blocks ( two of them)
      DIMENSION VEC1(*),VEC2(*)
*
      NTEST = 000
      LBLK = -1
      IDUM = 1
*. Transfer vector on LUCIN to LUSC1
C           COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
      CALL  COPVCD(LUCIN,LUSC1,VEC1,1,LBLK)
*. A bit of info for the sigma routine
      I_RES_AB = 0
*. Do the one-electron update
        I12 = 1
*. With 1-electron integrals in complete block form
        IH1FORM = 2
*. Transform each orbital separately
      DO K = 1, NTOOB
*. Place (T(P,K)/S(K,K)   in one-electron integral list
C                       T_ROW_TO_H(T,H,K)
        CALL T_ROW_TO_H(T,WORK(KINT1),K,TKK)
*. T_{kk}^Nk
C            T_TO_NK_VEC(T,KORB,ISM,ISPC,LUCIN,LUCOUT,C)
        CALL T_TO_NK_VEC(TKK,K,ISSM,ISSPC,LUSC1,LUSC2,VEC1)
        CALL COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' output from T_TO_NK'
          CALL WRTVCD(VEC1,LUSC1,1,LBLK)
        END IF
*. For each orbital calculate (1+T+1/2 T^2)|0>
* + T
        CALL MV7(VEC1,VEC2,LUSC1,LUSC2)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Correction vector'
          CALL WRTVCD(VEC1,LUSC2,1,LBLK)
        END IF
        ONE = 1.0D0
        CALL VECSMDP(VEC1,VEC2,ONE,ONE,LUSC1,LUSC2,LUSC3,1,LBLK)
        CALL COPVCD(LUSC3,LUSC1,VEC1,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Updated vector'
          CALL WRTVCD(VEC1,LUSC1,1,LBLK)
        END IF
*. + 1/2 T^2
        CALL MV7(VEC1,VEC2,LUSC2,LUSC3)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Correction vector'
          CALL WRTVCD(VEC1,LUSC3,1,LBLK)
        END IF
        ONE = 1.0D0
        HALF  = 0.5D0
        CALL VECSMDP(VEC1,VEC2,ONE,HALF,LUSC1,LUSC3,LUSC2,1,LBLK)
*. and transfer back to LUSC1
        CALL COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' Updated vector'
          CALL WRTVCD(VEC1,LUSC1,1,LBLK)
        END IF
      END DO
*. And transfer to LUCOUT
      CNORM = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
      WRITE(6,*) ' Norm of transformed vector', CNORM
C?    WRITE(6,*) ' Transformed vector'
C?    CALL WRTVCD(VEC1,LUSC1,1,LBLK)
      CALL REWINE(LUSC1,-1)
C?    WRITE(6,*) ' LUCOUT LUSC1 = ', LUCOUT,LUSC1
      CALL COPVCD(LUSC1,LUCOUT,VEC1,0,LBLK)
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRAH1(NBAS,NORB,NSM,HAO,C,HMO,IHSM,SCR)
*
*. Transform one-electron integrals from ao's to mo's.
*
*. Symmetry of integrals is IHSM, all integrals blocks assumed complete,
* i.e not packed to lower half
*
* Jeppe Olsen
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION HAO(*),C(*)
      DIMENSION NORB(*),NBAS(*)
#include "multd2h.inc"

*. Output
      DIMENSION  HMO(*)
*. Scratch
      DIMENSION SCR(*)
*. Loop over integral blocks
      IBHAO = 1
      IBHMO = 1
      DO IRSM = 1, NSM
        ICSM = MULTD2H(IRSM,IHSM)
*. Pointers to offsets in transformation matrices
        IBR = 1
        DO ISM = 1, IRSM-1
          IBR = IBR + NORB(ISM)*NBAS(ISM)
        END DO
        IBC = 1
        DO ISM = 1, ICSM-1
          IBC = IBC + NORB(ISM)*NBAS(ISM)
        END DO
*.
        LRMO = NORB(IRSM)
        LRAO = NBAS(IRSM)
*
        LCMO = NORB(ICSM)
        LCAO = NBAS(ICSM)
C       write(6,*) ' TRAH1 : IRSM ICSM ',IRSM,ICSM
C       WRITE(6,*) ' LRAO LRMO LCAO LCMO ',LRAO,LRMO,LCAO,LCMO

*
C            MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
C    &             NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
        ZERO = 0.0D0
        ONE= 1.0D0
*.C(row)T*Hao
        CALL MATML7(SCR,C(IBR),HAO(IBHAO),
     &       LRMO,LCAO,LRAO,LRMO,LRAO,LCAO,ZERO,ONE,1)
*. (C(row)T*Hao)*C(column)
        CALL MATML7(HMO(IBHMO),SCR,C(IBC),
     &       LRMO,LCMO,LRMO,LCAO,LCAO,LCMO,ZERO,ONE,0)
*
        IBHAO = IBHAO + LRAO*LCAO
        IBHMO = IBHMO + LRMO*LCMO
*.
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRAINT
*
* Very simple integral transformation for LUCIA
* - for testing the real one
*
* MO => MO transformation defined by matrix C
*
* Jeppe Olsen, One day in November 98, Magistratsvaegen 37D
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
*
      CALL QENTER('TRAIN')
*
      WRITE(6,*) ' TRACTL(Fusk) entered '
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Transformation matrix C '
        CALL APRBLM2(WORK(KMOMO),NTOOBS,NTOOBS,NSMOB,0)
C            APRBLM2(A,LROW,LCOL,NBLK,ISYM)
      END IF

*
*. Transform 1-electron integrals :
*
      CALL TRA1_SIMPLE(WORK(KMOMO))
*
*. Transform 2-electron integrals
*
      CALL TRA2_SIMPLE(WORK(KMOMO))
*
      CALL QEXIT ('TRAIN')
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRAN_SYM_BLOC_MAT(AIN,X,NBLOCK,LBLOCK,AOUT,SCR)
*
* Transform a blocked symmetric matrix AIN with blocked matrix
*  X to yield blocked matrix AOUT
*
* Aout = X(transposed) A X
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION AIN(*),X(*),LBLOCK(NBLOCK)
*. Output
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block
      DIMENSION SCR(*)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP = 1
         IOFFC = 1
       ELSE
         IOFFP = IOFFP + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
         IOFFC = IOFFC + LBLOCK(IBLOCK-1)** 2
       END IF
       L = LBLOCK(IBLOCK)
       K1 = 1
       K2 = 1 + L **2
*. Unpack block of A
C      TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
       SIGN = 1.0D0
       CALL TRIPAK_LUCI(SCR(K1),AIN(IOFFP),2,L,L,SIGN)
*. X(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
       CALL MATML7(SCR(K2),X(IOFFC),SCR(K1),L,L,L,L,L,L,
     &             ZERO,ONE,1)
*. X(T) (IBLOCK) A(IBLOCK) X (IBLOCK)
       CALL MATML7(SCR(K1),SCR(K2),X(IOFFC),L,L,L,L,L,L,
     &             ZERO,ONE,0)
*. Pack and transfer
       CALL TRIPAK_LUCI(SCR(K1),AOUT(IOFFP),1,L,L,SIGN)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        CALL APRBLM2(AOUT,LBLOCK,LBLOCK,NBLOCK,1)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRAN_SYM_BLOC_MAT2(AIN,X,NBLOCK,LBLOCK,AOUT,SCR,ISYM)
*
* Transform a blocked matrix AIN with blocked matrix
*  X to yield blocked matrix AOUT
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Aout = X(transposed) A X
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION AIN(*),X(*),LBLOCK(NBLOCK)
*. Output
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block
      DIMENSION SCR(*)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP = 1
         IOFFC = 1
       ELSE
         IOFFP = IOFFP + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
         IOFFC = IOFFC + LBLOCK(IBLOCK-1)** 2
       END IF
       L = LBLOCK(IBLOCK)
       K1 = 1
       K2 = 1 + L **2
*. Unpack block of A
C      TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
       SIGN = 1.0D0
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK_LUCI(SCR(K1),AIN(IOFFP),2,L,L,SIGN)
       ELSE
         CALL COPVEC(AIN(IOFFC),SCR(K1),L*L)
       END IF
*. X(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
       CALL MATML7(SCR(K2),X(IOFFC),SCR(K1),L,L,L,L,L,L,
     &             ZERO,ONE,1)
*. X(T) (IBLOCK) A(IBLOCK) X (IBLOCK)
       CALL MATML7(SCR(K1),SCR(K2),X(IOFFC),L,L,L,L,L,L,
     &             ZERO,ONE,0)
*. Pack and transfer
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK_LUCI(SCR(K1),AOUT(IOFFP),1,L,L,SIGN)
       ELSE
         CALL COPVEC(SCR(K1),AOUT(IOFFC),L*L)
       END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        CALL APRBLM2(AOUT,LBLOCK,LBLOCK,NBLOCK,ISYM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TRAN_SYM_BLOC_MAT3
     &   (AIN,X,NBLOCK,LX_ROW,LX_COL,AOUT,SCR,ISYM)
*
* Transform a blocked matrix AIN with blocked matrix
*  X to yield blocked matrix AOUT
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Aout = X(transposed) A X
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION AIN(*),X(*),LX_ROW(NBLOCK),LX_COL(NBLOCK)
*. Output
      DIMENSION AOUT(*)
*. Scratch : At least twice the length of largest block
      DIMENSION SCR(*)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1) THEN
         IOFFP_IN = 1
         IOFFC_IN = 1
         IOFFP_OUT= 1
         IOFFC_OUT= 1
         IOFFX = 1
       ELSE
         IOFFP_IN = IOFFP_IN + LX_ROW(IBLOCK-1)*(LX_ROW(IBLOCK-1)+1)/2
         IOFFC_IN = IOFFC_IN + LX_ROW(IBLOCK-1) ** 2
         IOFFP_OUT= IOFFP_OUT+ LX_COL(IBLOCK-1)*(LX_COL(IBLOCK-1)+1)/2
         IOFFC_OUT= IOFFC_OUT+ LX_COL(IBLOCK-1) ** 2
         IOFFX = IOFFX + LX_ROW(IBLOCK-1)*LX_COL(IBLOCK-1)
       END IF
       LXR = LX_ROW(IBLOCK)
       LXC = LX_COL(IBLOCK)
       K1 = 1
       K2 = 1 + MAX(LXR,LXC) ** 2
*. Unpack block of A
       SIGN = 1.0D0
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK_LUCI(SCR(K1),AIN(IOFFP_IN),2,LXR,LXR,SIGN)
C             TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGN)
       ELSE
         CALL COPVEC(AIN(IOFFC_IN),SCR(K1),LXR*LXR)
       END IF
*. X(T)(IBLOCK)A(IBLOCK)
       ZERO = 0.0D0
       ONE  = 1.0D0
       CALL SETVEC(SCR(K2),ZERO,LXR*LXC)
       CALL MATML7(SCR(K2),X(IOFFX),SCR(K1),LXC,LXR,LXR,LXC,LXR,LXR,
     &             ZERO,ONE,1)
C?     WRITE(6,*) ' Half transformed matrix '
C?     CALL WRTMT_LU(SCR(K2),LXC,LXR,LXC,LXR)

*. X(T) (IBLOCK) A(IBLOCK) X (IBLOCK)
       CALL SETVEC(SCR(K1),ZERO,LXC*LXC)
       CALL MATML7(SCR(K1),SCR(K2),X(IOFFX),LXC,LXC,LXC,LXR,LXR,LXC,
     &             ZERO,ONE,0)
C?     WRITE(6,*) ' Transformed matrix '
C?     CALL WRTMT_LU(SCR(K1),LXC,LXC,LXC,LXC)
*. Pack and transfer
       IF(ISYM.EQ.1) THEN
         CALL TRIPAK_LUCI(SCR(K1),AOUT(IOFFP_OUT),1,LXC,LXC,SIGN)
       ELSE
         CALL COPVEC(SCR(K1),AOUT(IOFFC_OUT),LXC*LXC)
       END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' output matrix TRAN_SYM_BLOC_MAT '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        CALL APRBLM2(AOUT,LX_COL,LX_COL,NBLOCK,ISYM)
      END IF
*
      END

***********************************************************************
      SUBROUTINE SCLH2(XLAMBDA)
*
*. Scale two electron integrals
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      call dscal(nint2,xlambda,work(kint2),1)
*
      END
