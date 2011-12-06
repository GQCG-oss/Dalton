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
c /* deck cc_gethfgd */
*=====================================================================*
      SUBROUTINE CC_GETHFGD(IVEC,TYPE,LABEL,ISYMS,ISTAT,EIGV,ISYMO,
     &                      FREQS,ICAU,NVEC,MAXVEC,IREAL,
     &                      CMO,UDV,XINDX,FRVAL,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: get right hand side vector for CPHF equation
*             vector is returned at the beginning of the work space
*
*    implemented types:  R1  
*
*    Written by Christof Haettig, november 1998
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccfro.h>
#include <ccroper.h>
#include <ccexpfck.h>
#include <dummy.h>

* local parameters:
      CHARACTER*(20) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_GETHFGD> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE. )
      INTEGER LFOCK

      CHARACTER TYPE*(*)

      INTEGER IREAL, IVEC, NVEC, MAXVEC, LWORK
      INTEGER ISYMS(MAXVEC,*), ISYMO(MAXVEC,*)
      INTEGER ISTAT(MAXVEC,*), ICAU(MAXVEC,*)

      CHARACTER*8 LABEL(MAXVEC,*)

#if defined (SYS_CRAY)
      REAL FREQS(MAXVEC,*), EIGV(MAXVEC,*), FRVAL
      REAL CMO(*), UDV(*), XINDX(*), WORK(LWORK)
      REAL ZERO, ONE, TWO, SIGN
      REAL XNORM, DDOT
#else
      DOUBLE PRECISION FREQS(MAXVEC,*), EIGV(MAXVEC,*), FRVAL
      DOUBLE PRECISION CMO(*), UDV(*), XINDX(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, TWO, SIGN
      DOUBLE PRECISION XNORM, DDOT
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0)

      CHARACTER*8 LABL1
      INTEGER ISYM, IOPER, IOPT, IADR1, IADR0
      INTEGER ISYM0, ISYM1, ISYM2, ICOUNT, ICMO(8,8), NCMO(8)
      INTEGER KOFF1, KOFF2, KOFF3, KOFF4, KOFF5, KOFF6, KOFF7
      INTEGER NBASA, NBASB, NBSA0, NBSB0, NVIRA
      INTEGER KFOCK1, KFOCK0, KFCKMO, KSCR1, KSCR0, KRMAT, KKAPPA
      INTEGER IFOCK1, IFOCK0, KEND1, LWRK1
      INTEGER KCMOPQ, KCMOHQ, KLAMDPQ, KLAMDHQ
      INTEGER ISYMA, ISYMI, ISYALP, ISYBET, ISYAL0, ISYBT0, KT1AM

* external functions:
      INTEGER IROPER
      INTEGER IEFFFOCK

      CALL QENTER('GETHFGP')
*---------------------------------------------------------------------*
* initialize some symmetry arrays:
*---------------------------------------------------------------------*
      DO ISYM = 1, NSYM
         ICOUNT = 0
         DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            ICMO(ISYM1,ISYM2) = ICOUNT
            ICOUNT = ICOUNT + NBAS(ISYM1)*NORBS(ISYM2)
         END DO
         NCMO(ISYM) = ICOUNT
      END DO      

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*

      IF (TYPE(1:2).EQ.'R1') THEN

        LABL1 = LABEL(IVEC,1)
        FRVAL = FREQS(IVEC,1)
        IOPER = IROPER(LABL1,ISYM)
        IREAL = ISYMAT(IOPER)

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_GETHFGD> TYPE:',TYPE
           WRITE (LUPRI,*) 'CC_GETHFGD> IVEC:',IVEC
           WRITE (LUPRI,*) 'CC_GETHFGD> LABL1:',LABL1
           WRITE (LUPRI,*) 'CC_GETHFGD> FRVAL:',FRVAL
           WRITE (LUPRI,*) 'CC_GETHFGD> IOPER:',IOPER
           WRITE (LUPRI,*) 'CC_GETHFGD> IREAL:',IREAL
           WRITE (LUPRI,*) 'CC_GETHFGD> LWORK:',LWORK
        END IF

        IF (LWORK.LT.NALLAI(ISYM)) THEN
           WRITE(LUPRI,*) 'have:',LWORK
           WRITE(LUPRI,*) 'need:',NALLAI(ISYM)
           CALL QUIT('Insufficient memory in CC_GETHFGD.')
        END IF

        ISYM0  = 1
        ISYM1  = ISYM

        KFCKMO = 1
        KFOCK0 = KFCKMO + 2*NALLAI(ISYM1)
        KFOCK1 = KFOCK0 + N2BST(ISYM0)
        KCMOPQ = KFOCK1 + N2BST(ISYM1)
        KCMOHQ = KCMOPQ + NCMO(ISYM1)
        KRMAT  = KCMOHQ + NCMO(ISYM1)
        KKAPPA = KRMAT  + N2BST(ISYM1)
        KSCR1  = KKAPPA + 2*NALLAI(ISYM1)
        KSCR0  = KSCR1  + NT1AO(ISYM1)
        KEND1  = KSCR0  + NT1AO(ISYM0)
        LWRK1  = LWORK  - KEND1

        IF (LWRK1.LT.0) THEN
           WRITE(LUPRI,*) 'have:',LWORK
           WRITE(LUPRI,*) 'need:',KEND1
           CALL QUIT('Insufficient memory in CC_GETHFGD.')
        END IF
 
        ! get the derivative part of the AO Fock operator for the rhs
        IFOCK1 = IEFFFOCK(LABL1,     ISYM1,1)
        IFOCK0 = IEFFFOCK('HAM0    ',ISYM0,1)

        IADR1  = IADRFCK(2,IFOCK1)    
        IADR0  = IADRFCK(2,IFOCK0)    

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_GETHFGD> IFOCK1 = ',IFOCK1
           WRITE (LUPRI,*) 'CC_GETHFGD> IFOCK0 = ',IFOCK0
           WRITE (LUPRI,*) 'CC_GETHFGD> IADR1  = ',IADR1 
           WRITE (LUPRI,*) 'CC_GETHFGD> IADR0  = ',IADR0 
           WRITE (LUPRI,*) 'CC_GETHFGD> KFOCK1 = ',KFOCK1 
           WRITE (LUPRI,*) 'CC_GETHFGD> KFOCK0 = ',KFOCK0 
        END IF

        LFOCK = -1
        CALL WOPEN2(LFOCK,FILFCKEFF,64,0)

        CALL GETWA2(LFOCK,FILFCKEFF,WORK(KFOCK0),IADR0,N2BST(ISYM0))
        CALL GETWA2(LFOCK,FILFCKEFF,WORK(KFOCK1),IADR1,N2BST(ISYM1))

        CALL WCLOSE2(LFOCK,FILFCKEFF,'KEEP')

        IF (LOCDBG) THEN
           WRITE (LUPRI,*)
     &          'CC_GETHFGD> The Fock1 matrix as read from file:'
           CALL CC_PRONELAO(WORK(KFOCK1),ISYM1)
           WRITE (LUPRI,*) 
     &          'CC_GETHFGD> The Fock0 matrix as read from file:'
           CALL CC_PRONELAO(WORK(KFOCK0),ISYM0)
        END IF

        ! get the orbital connection matrix:
        CALL CC_GET_RMAT(WORK(KRMAT),IOPER,1,ISYM,WORK(KEND1),LWRK1)
        
        ! calculate Rmat transformed CMO vector

        IOPT  = 0
        CALL DZERO(WORK(KKAPPA),2*NALLAI(ISYM))
        CALL CC_LAMBDAQ(DUMMY,DUMMY,
     &                  WORK(KCMOPQ), WORK(KCMOHQ),
     &                  ISYM,DUMMY,WORK(KKAPPA),WORK(KRMAT),
     &                  IREAL,IOPT,WORK(KEND1),LWRK1)
        
        
        CALL DZERO(WORK(KSCR0), NT1AO(ISYM0))
        CALL DZERO(WORK(KSCR1), NT1AO(ISYM1))
        CALL DZERO(WORK(KFCKMO),2*NALLAI(ISYM1))
CCH
cch     CALL DZERO(WORK(KFOCK1),N2BST(ISYM1))
cch     CALL DZERO(WORK(KFOCK0),N2BST(ISYM0))
CCH

        DO ISYMA = 1, NSYM
           ISYMI  = MULD2H(ISYM1,ISYMA)
           ISYBT0 = MULD2H(ISYM1,ISYMI)
           ISYAL0 = MULD2H(ISYM1,ISYMA)
           ISYALP = ISYMA
           ISYBET = ISYMI

           NBASA = MAX(NBAS(ISYALP),1)
           NBASB = MAX(NBAS(ISYBET),1)
           NBSA0 = MAX(NBAS(ISYAL0),1)
           NBSB0 = MAX(NBAS(ISYBT0),1)
           NVIRA = MAX(NVIRS(ISYMA),1)

C          ------------------------------------------------------
C          transform second index of FOCK1 to occupied using CMO:
C          ------------------------------------------------------
           KOFF1 = IAODIS(ISYALP,ISYBET)
           KOFF2 = ICMO(ISYBET,ISYMI)
           KOFF3 = IT1AO(ISYALP,ISYMI)

           CALL DGEMM('N','N',NBAS(ISYALP),NRHFS(ISYMI),NBAS(ISYBET),
     &                ONE, WORK(KFOCK1+KOFF1), NBASA,
     &                     CMO(1      +KOFF2), NBASB,
     &                ONE, WORK(KSCR1 +KOFF3), NBASA)


C          --------------------------------------------------------
C          transform second index of FOCK0 to occupied using CMOHQ:
C          --------------------------------------------------------
           KOFF4  = IAODIS(ISYALP,ISYBT0)
           KOFF5  = ICMO(ISYBT0,ISYMI)
           KOFF3  = IT1AO(ISYALP,ISYMI)


           CALL DGEMM('N','N',NBAS(ISYALP),NRHFS(ISYMI),NBAS(ISYBT0),
     &                ONE, WORK(KFOCK0+KOFF4), NBASA,
     &                     WORK(KCMOHQ+KOFF5), NBSB0,
     &                ONE, WORK(KSCR1 +KOFF3), NBASA)

C          ------------------------------------------------------------
C          transform first index of halftrf. Fock to virtual using CMO:
C          ------------------------------------------------------------
           KOFF3 = IT1AO(ISYALP,ISYMI)
           KOFF6 = ICMO(ISYALP,ISYMA) + NBAS(ISYALP)*NRHFS(ISYMA)
           KOFF7 = IT1AM(ISYMA,ISYMI)

           CALL DGEMM('T','N',NVIRS(ISYMA),NRHFS(ISYMI),NBAS(ISYALP),
     &                ONE, CMO(1      +KOFF6), NBASA,
     &                     WORK(KSCR1 +KOFF3), NBASA,
     &                ONE, WORK(KFCKMO+KOFF7), NVIRA)


C          ------------------------------------------------------
C          transform second index of FOCK0 to occupied using CMO:
C          ------------------------------------------------------
           KOFF1  = IAODIS(ISYAL0,ISYBET)
           KOFF2  = ICMO(ISYBET,ISYMI)
           KOFF3  = IT1AO(ISYAL0,ISYMI)

           CALL DGEMM('N','N',NBAS(ISYAL0),NRHFS(ISYMI),NBAS(ISYBET),
     &                ONE, WORK(KFOCK0+KOFF1), NBSA0,
     &                     CMO(1      +KOFF2), NBASB,
     &                ONE, WORK(KSCR0 +KOFF3), NBSA0)

C          -------------------------------------------------------------
C          transform first index of halftrf. Fock to virtual with CMOPQ:
C          -------------------------------------------------------------
           KOFF3 = IT1AO(ISYAL0,ISYMI)
           KOFF6 = ICMO(ISYAL0,ISYMA) + NBAS(ISYAL0)*NRHFS(ISYMA)
           KOFF7 = IT1AM(ISYMA,ISYMI)

C          SIGN = DBLE(IREAL)
C          SIGN = 1.0D0
C          CALL DGEMM('T','N',NVIRS(ISYMA),NRHFS(ISYMI),NBAS(ISYAL0),
C    &                SIGN,WORK(KCMOQ +KOFF6), NBSA0,
C    &                     WORK(KSCR0 +KOFF3), NBSA0,
C    &                ONE, WORK(KFCKMO+KOFF7), NVIRA)

           CALL DGEMM('T','N',NVIRS(ISYMA),NRHFS(ISYMI),NBAS(ISYAL0),
     &                ONE,WORK(KCMOPQ+KOFF6), NBSA0,
     &                    WORK(KSCR0 +KOFF3), NBSA0,
     &                ONE,WORK(KFCKMO+KOFF7), NVIRA)

        END DO

        CALL DSCAL(NALLAI(ISYM1),TWO,WORK(KFCKMO),1)
        CALL DCOPY(NALLAI(ISYM1),WORK(KFCKMO),1,
     &                           WORK(KFCKMO+NALLAI(ISYM1)),1)
        IF (IREAL.EQ.1) THEN
           CALL DSCAL(NALLAI(ISYM1),-1.0d0,WORK(KFCKMO+NALLAI(ISYM1)),1)
        END IF

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_GETHFGD> The Fock1 matrix:'
           CALL CC_PRONELAO(WORK(KFOCK1),ISYM1)
           WRITE (LUPRI,*) 'CC_GETHFGD> The Fock0 matrix:'
           CALL CC_PRONELAO(WORK(KFOCK0),ISYM0)
           WRITE (LUPRI,*) 
     &          'CC_GETHFGD> The AI block of the Fock matrix:'
           CALL CC_PRP(WORK(KFCKMO),WORK,ISYM1,1,0)
           WRITE (LUPRI,*) 'CC_GETHFGD> The the kappa rhs vector:'
           WRITE (LUPRI,'(5X,I5,F12.8)') (I,WORK(I),I=1,NALLAI(ISYM))
        END IF

      ELSE
        WRITE (LUPRI,*) 'CPHF rhs vectors ',TYPE(1:2),
     &        ' not implemented.'
        CALL QUIT('required rhs CPHF vectors not implemented.')
      END IF

      CALL QEXIT('GETHFGP')
      RETURN
      END

*=====================================================================*
*                    END OF SUBROUTINE CC_GETHFGD                     *
*=====================================================================*
