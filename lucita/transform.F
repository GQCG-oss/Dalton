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

!**********************************************************************

      SUBROUTINE TRACI_CTL(vec1,vec2,c2,mo2mo_mat,par_dist_block_list)
*
* Master routine for transforming CI vectors to new orbital basis
*
* Input vectors on LUC
* output vectors delivered on LUHC
* (LUSC1,LUSC2) used as scratch files
*
* Jeppe Olsen, January 98
! re-factored and parallelized by Stefan Knecht, November 2011
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLH1SAVE
!               for addressing of WORK
      real(8), intent(inout) :: vec1(*)
      real(8), intent(inout) :: vec2(*)
      real(8), intent(inout) :: c2(*)
      real(8), intent(inout) :: mo2mo_mat(*)
      integer, intent(in)    :: par_dist_block_list(*)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "clunit.inc"
#include "crun.inc"
#include "cstate.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "parluci.h"
      integer, parameter  :: isigden = 2
      real(8)             :: cdummy(2), sdummy(2)
*
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
!     NTEST = 1000 ! debug
      NTEST = 0    ! default
      IDUM  = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACIC')
*
#ifdef LUCI_DEBUG
      WRITE(luwrt,*)
      WRITE(luwrt,*) ' Welcome to TRACI_CTL'
      WRITE(luwrt,*) ' ===================='
      WRITE(luwrt,*)
#endif

*. Copy of one-electron integrals
      CALL MEMMAN(KLH1SAVE,NTOOB**2,'ADDL  ',2,'H1SAVE')
*. We are going to mess with the one-electron integrals, take a copy
      CALL COPVEC(WORK(KINT1),WORK(KLH1SAVE),NTOOB*NTOOB)
*
      CALL REWINE(LUC,-1)
      CALL REWINE(luhc,-1)

      DO JROOT = 1, NROOT
!
!       stefan - nov 2011: the mo2mo transformation matrix is provided from the outer 
!                          CI/MCSCF environment - the following two steps are therefore commented out. 
!
!       step 1: get PV for state JROOT
!       ------------------------------
!       call rewine(lusc1,-1)
!       CALL COPVCD(lusc1,lusc2,VEC1,0,-1)
!       call rewine(lusc1,-1)
!       call rewine(lusc2,-1)
!       CALL SIGDEN_CI(VEC1,VEC2,C2,LUSC1,LUSC2,
!    &                 cdummy,sdummy,ISIGDEN,1)
!       write(luwrt,*) ' density matrix generated '
!       call dcopy(NTOOB**2,mo2mo_mat,1,WORK(KRHO1),1)
!       write(luwrt,*) ' density matrix copied '
!       call wrtmatmn(WORK(KRHO1),1,NTOOB**2,1,NTOOB**2,luwrt)

!       step 2: Obtain MO-MO transformation matrix 
!       ------------------------------------------
!       CALL GET_CMOMO(mo2mo_mat)
*
*. Transform CI vector : Input on lusc1, output on luhc
        call rewine(lusc1,-1)
        CALL COPVCD(LUC,lusc1,VEC1,0,-1)
        call rewine(lusc1,-1)

        CALL TRACID(mo2mo_mat,lusc1,luhc,LUSC2,LUSC3,
     &              VEC1,VEC2,C2,par_dist_block_list)

      END DO ! loop over roots

#ifdef LUCI_DEBUG
      CALL REWINE(luhc,-1)
      IF(NTEST.GE.100) THEN
        DO JROOT = 1, NROOT
          WRITE(luwrt,*)
          WRITE(luwrt,*) ' Analysis of rotated state number ', JROOT
          WRITE(luwrt,*) ' ====================================='
          WRITE(luwrt,*)
          CALL WRTVCD(VEC1,luhc,0,-1)
        END DO
      END IF
#undef LUCI_DEBUG
#endif

!     clean up time : restore 1-e integrals
      CALL COPVEC(WORK(KLH1SAVE),WORK(KINT1),NTOOB*NTOOB)

!     release local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACIC')

      END
!**********************************************************************

      SUBROUTINE TRACID(T,LUCIN,LUCOUT,LUSC1,LUSC2,
     &                  VEC1,VEC2,C2,par_dist_block_list)
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
      INTEGER*8 KLCLBTx, KLCLEBTx, KLCI1BTx, KLCIBTx, KLCBLTPx
!               for addressing of WORK
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "oper.inc"
#include "intform.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "cicisp.inc"
#include "cstate.inc"
#include "parluci.h"
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
      REAL*8 INPRDD
*. Input
      DIMENSION T(*)
      integer, intent(in) :: par_dist_block_list(*)
      integer, parameter  :: isigden = 1
      real(8)             :: cdummy(2), sdummy(2)
*. Scratch blocks ( three of them)
      DIMENSION VEC1(*),VEC2(*), C2(*)

!     NTEST = 1000 ! debug
      NTEST = 0000 ! default

!     a bit of info for the sigma routine
      I_RES_AB = 0

!     1-electron integrals in complete block form (NTOOB,NTOOB)
      IH1FORM  = 2

!     set up block and batch structure of vector
      IDUM     = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRACID')

      CALL MEMMAN(KLCLBTx ,  MXNTTS,'ADDL  ',1,'CLBTx ')
      CALL MEMMAN(KLCLEBTx,  MXNTTS,'ADDL  ',1,'CLEBTx')
      CALL MEMMAN(KLCI1BTx,  MXNTTS,'ADDL  ',1,'CI1BTx')
      CALL MEMMAN(KLCIBTx ,8*MXNTTS,'ADDL  ',1,'CIBTx ')

      IATP     = 1
      IBTP     = 2
      CALL Z_BLKFO_partitioning_parallel(ISSPC,ISSM,iatp,ibtp,
     &                                   work(KLCLBTx),
     &                                   work(KLCLEBTx),
     &                                   work(KLCI1BTx),
     &                                   work(KLCIBTx),
     &                                   NBATCH,NBLOCK,
     &                                   par_dist_block_list)

      CALL MEMMAN(KLCBLTPx,       NSMST,'ADDL  ',2,'CBLTPA')
      CALL ZBLTP_IDC(ISMOST(1,ISSM),NSMST,IDC,WORK(KLCBLTPx))

#ifdef LUCI_DEBUG
      write(luwrt,*) 'ntoob is ==> ',ntoob
      IF(NTEST.GE.1000) THEN
        WRITE(luwrt,*) ' input vector to ntoob loop'
        CALL WRTVCD(VEC1,lucin,1,-1)
      END IF
#endif

      call rewine(lucin,-1)

!     transform each orbital separately
      DO K = 1, NTOOB

!       place (T(P,K)/S(K,K) in one-electron integral list work(kint1)
        CALL T_ROW_TO_H(T,WORK(KINT1),K,TKK)

!       scale elements with T_{kk}^Nk depending on occupation type Nk (Nk=0,1,2)
        CALL T_TO_NK_VEC(TKK,K,ISSM,ISSPC,lucin,LUSC1,VEC1,
     &                   work(KLCIBTx),work(KLCBLTPx),nblock)

        CALL COPVCD(LUSC1,lucin,VEC1,1,-1)

#ifdef LUCI_DEBUG
        IF(NTEST.GE.1000) THEN
          WRITE(luwrt,*) ' output from T_TO_NK'
          CALL WRTVCD(VEC1,lucin,1,-1)
        END IF
#endif
        call rewine(lucin,-1)
        call rewine(lusc1,-1)
!       for each orbital calculate (1+T+1/2 T^2)|0>
!       step 1: + T
        CALL SIGDEN_CI(VEC1,VEC2,C2,lucin,LUSC1,
     &                 cdummy,sdummy,ISIGDEN,-1)
        
#ifdef LUCI_DEBUG
        IF(NTEST.GE.1000) THEN
          WRITE(luwrt,*) ' Correction vector'
          CALL WRTVCD(VEC1,LUSC1,1,-1)
        END IF
#endif

        CALL VECSMD(VEC1,VEC2,1.0d0,1.0d0,lucin,LUSC1,lusc2,1,-1)
        CALL COPVCD(lusc2,lucin,VEC1,1,-1)

#ifdef LUCI_DEBUG
        IF(NTEST.GE.1000) THEN
          WRITE(luwrt,*) ' Updated vector'
          CALL WRTVCD(VEC1,lucin,1,-1)
        END IF
#endif

        CALL rewine(lusc2,-1)
        CALL rewine(lusc1,-1)

!       step 2: + 1/2 T^2
        CALL SIGDEN_CI(VEC1,VEC2,C2,LUSC1,lusc2,
     &                 cdummy,sdummy,ISIGDEN,-1)

        IF(NTEST.GE.1000) THEN
          WRITE(luwrt,*) ' Correction vector'
          CALL WRTVCD(VEC1,lusc2,1,-1)
        END IF

        CALL VECSMD(VEC1,VEC2,1.0d0,0.5d0,lucin,lusc2,LUSC1,1,-1)

!       transfer back to LUSC1
        CALL COPVCD(LUSC1,lucin,VEC1,1,-1)

        IF(NTEST.GE.1000) THEN
          WRITE(luwrt,*) ' Updated vector'
          CALL WRTVCD(VEC1,lucin,1,-1)
        END IF
      END DO

#ifdef LUCI_DEBUG
      CNORM = INPRDD(VEC1,VEC2,lucin,lucin,1,-1)
      WRITE(luwrt,*) ' Norm of transformed vector', CNORM
#endif

!     transfer rotated vector to LUCOUT
      CALL REWINE(lucin,-1)
      CALL COPVCD(lucin,LUCOUT,VEC1,0,-1)

!     release local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRACID')

      END
!**********************************************************************

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
#include "parluci.h"
*. Input ( in blocked form)
      DIMENSION T(*)
*. Output ( also in blocked form)
      DIMENSION H(*)
*
      CALL DZERO(H,NTOOB**2)

      KSM  = ISMFSO(K)
      KOFF = IBSO(KSM)
      KREL = K - KOFF + 1
      NK   = NTOOBS(KSM)
!     symmetry offset in full 1e-matrix block
      IOFF = IFRMR(WORK(KPGINT1A(1)),1,KSM)

#ifdef LUCI_DEBUG
      write(luwrt,*) 'ioff in T and H ==> ',IOFF+(KREL-1)*NK
      write(luwrt,*) 'NK elements     ==> ',NK
      write(luwrt,*) 'absolute ioff   ==> ',ioff
      write(luwrt,*) 'krel, ksm       ==> ',krel,ksm
#endif
!     copy orbital block with symmetry ksm
      CALL DCOPY(NK,T(IOFF+(KREL-1)*NK),1,H(IOFF+(KREL-1)*NK),1)

!     set T_kk
      TKK = H(IOFF-1+(KREL-1)*NK+KREL)
      IF(TKK .NE. 0.0D0) THEN
        FAC = 1.0D0/TKK
        CALL DSCAL(NK,FAC,H(IOFF+(KREL-1)*NK),1)
        H(IOFF-1+(KREL-1)*NK+KREL) = 0.0D0
      END IF

#ifdef LUCI_DEBUG
      WRITE(luwrt,*) ' output from T_ROW_TO_H, TKK ==> ',TKK
      WRITE(luwrt,*) ' Row to be transferred ', KREL
      WRITE(luwrt,*) ' Updated H matrix '
      WRITE(luwrt,*) ' NTOOB, NTOOBS(1)', NTOOB, NTOOBS(1)
      CALL APRBLM2(H,NTOOBS,NTOOBS,NSMOB,0)
#undef LUCI_DEBUG
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE T_TO_NK_VEC(T,KORB,ISM,ISPC,LUCIN,LUCOUT,C,
     &                       block_info,block_type_info,nblock)
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
      INTEGER*8 KLASTR, KLBSTR, KLKAOC, KLKBOC
!               for addressing of WORK
#include "priunit.h"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strinp.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "gasstr.inc"
#include "crun.inc"
#include "csm.inc"
      integer, intent(in) :: nblock
      integer, intent(in) :: block_type_info(*)
      integer, intent(in) :: block_info(8,nblock)

*. Scratch block, must hold a batch of blocks
      DIMENSION C(*)
*
#ifdef LUCI_DEBUG
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(lupri,*) ' T_TO_NK_VEC speaking '
        WRITE(lupri,*) ' ISM, ISPC = ', ISM,ISPC
      END IF
#endif

      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'T_TO_N')

      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
      CALL MEMMAN(KLKAOC,MXNSTR,     'ADDL  ',1,'KLKAOC')
      CALL MEMMAN(KLKBOC,MXNSTR,     'ADDL  ',1,'KLKBOC')
      
!     orbital K in type ordering
      KKORB = IREOST(KORB)
#ifdef LUCI_DEBUG
      write(lupri,*) 'korb, kkorb ==> ',korb, kkorb
#endif

      CALL T_TO_NK_VECS(T,KKORB,C,LUCIN,LUCOUT,
     &                  WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &                  NBLOCK,block_info,
     &                  NAEL,NBEL,WORK(KLASTR),WORK(KLBSTR),
     &                  block_type_info,NSMST,
     &                  ICISTR,NTOOB,WORK(KLKAOC),WORK(KLKBOC))

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
     &                        NBLOCK,IBLOCK,
     &                        NAEL,NBEL,IASTR,IBSTR,IBLTP,NSMST,
     &                        ICISTR,NORB,IKAOCC,IKBOCC)
*
* Multiply Vector in LUCIN with t **NK_op to yield vector on LUCOUT
*
* Both files are initially rewinded
*
*
* Jeppe Olsen, Feb. 1998
*

      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
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

#ifdef LUCI_DEBUG
        WRITE(lupri,*) ' IATP IBTP IASM IBSM ', IATP,IBTP,IASM,IBSM
        WRITE(lupri,*) ' IBLTP(IASM)         ', IBLTP(IASM) 
#endif

!       obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
!       occupation of orb KORB
        DO JSTR = 1, NASTR1
          KOCC = 0
          DO JAEL = 1, NAEL
            IF(IASTR(JAEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKAOCC(JSTR) = KOCC
        END DO

#ifdef LUCI_DEBUG
        WRITE(lupri,*) ' IKAOCC array '
        CALL IWRTMAMN(IKAOCC,1,NASTR1,1,NASTR1,lupri)
#endif

!       Obtain Beta  strings of sym IBSM and type IBTP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                           NORB,0,IDUM,IDUM)
        DO JSTR = 1, NBSTR1
          KOCC = 0
          DO JBEL = 1, NBEL
            IF(IBSTR(JBEL,JSTR).EQ.KORB) KOCC = 1
          END DO
          IKBOCC(JSTR) = KOCC
        END DO

#ifdef LUCI_DEBUG
        WRITE(lupri,*) ' IKBOCC array '
        CALL IWRTMAMN(IKBOCC,1,NBSTR1,1,NBSTR1,lupri)
#endif

        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF

        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)

!       WRITE(lupri,*) ' NIA NIB ', NIA,NIB

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

              IDET = IDET + 1
              KABOCC = IKAOCC(IA)+IKBOCC(IB)

#ifdef LUCI_DEBUG
              WRITE(lupri,*) ' KABOCC, C(IDET), T, T**2',
     &                         KABOCC, C(IDET), T, T**2
              WRITE(lupri,*) ' IA IB IDET',IA,IB,IDET
#endif

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

      END
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
!**********************************************************************

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
!**********************************************************************

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
!**********************************************************************

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
!**********************************************************************

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
!**********************************************************************

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
