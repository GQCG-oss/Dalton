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
C  /* Deck ccsd_triple */
      SUBROUTINE CC_FOPTRIPLES(OMEGA1,DENAB,DENIJ,T1AM,T2AM,
     *                         XLAMDP,XLAMDH,WORK,LWORK)
C
C     Written by K. Hald, Fall 2001 based on 
C     CCSD_TRIPLE() written by Henrik Koch 19-Sep-1994
C
C     Purpose:
C
C     Calculate the iterative triples corrections to the CCSD
C     equations.
C
C     NB! The T2 amplitudes are assumed to be a full square.
C
C
C     NB! It is assumed that the vectors are allocated in the following
C     order:
C           T1AM(*), OMEGA1(*), OMEGA2(*), T2AM(*), SCR(*), WRK(*).
C
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <maxash.h>
#include <aovec.h>
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION OMEGA1(*),DENAB(*),DENIJ(*),T1AM(*),T2AM(*)
      DIMENSION XLAMDP(*),XLAMDH(*),WORK(LWORK)
C
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
CCN #include <infind.h> not compatible with R12!
#include <ccisao.h>
#include <r12int.h>
#include <blocks.h>
#include <inftap.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
COMMENT COMMENT
COMMENT COMMENT
#include <mxcent.h>
#include <eribuf.h>
#include <ccfield.h>
#include <ccfop.h>
COMMENT COMMENT
COMMENT COMMENT
C
      LOGICAL ITERATIVE
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (NSYM .NE. 1) THEN
         CALL QUIT('No symmetry in this part yet')
      ENDIF
C
COMMENT COMMENT
COMMENT COMMENT
      write(lupri,*) 'Entered FOPTRIPLES'
      write(lupri,*) 'NLAMDT = ',nlamdt
      write(lupri,*) 'NLAMDS = ',nlamds
      write(lupri,*) 'NORBT  = ',norbt
      write(lupri,*) 'NORBTS = ',norbts
      call flshfo(lupri)
C------------------------
C     Dynamic Allocation.
C------------------------
C
      KCMO   = 1
      KOME1  = KCMO   + NLAMDS
COMMENT COMMENT
COMMENT COMMENT CHANGED NLAMDT TO NLAMDS
COMMENT COMMENT
      KOME2  = KOME1  + NT1AMX
      KOME11 = KOME2  + NT1AMX*NT1AMX
      KOME22 = KOME11 + NT1AMX
      KSCR1  = KOME22 + NT1AMX*NT1AMX
      KFOCK  = KSCR1  + NT1AMX
      KFOCKD = KFOCK  + N2BST(ISYMOP)
      KINT1S = KFOCKD + NORBTS
COMMENT COMMENT
COMMENT COMMENT CHANGED NORBT TO NORBTS
COMMENT COMMENT
      KINT2S = KINT1S + NT1AMX*NVIRT*NVIRT
      KINT3S = KINT2S + NRHFT*NRHFT*NT1AMX
      KINT4S = KINT3S + NT1AMX*NVIRT*NVIRT
      KINT1T = KINT4S + NRHFT*NRHFT*NT1AMX
      KINT2T = KINT1T + NT1AMX*NVIRT*NVIRT
      KIAJB  = KINT2T + NRHFT*NRHFT*NT1AMX
      KYIAJB = KIAJB  + NT1AMX*NT1AMX
      KT3AM  = KYIAJB + NT1AMX*NT1AMX
      KT3BAR = KT3AM  + NT1AMX*NT1AMX*NT1AMX
      KVIRT  = KT3BAR + NT1AMX*NT1AMX*NT1AMX
      KOCC   = KVIRT  + NVIRT*NVIRT*NVIRT*NRHFT
      KKAPAA = KOCC   + NRHFT*NRHFT*NRHFT*NVIRT
      KKAPII = KKAPAA + NVIRT
      KEND1  = KKAPII + NRHFT
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSD_TRIPLE')
      ENDIF
C
C--------------------------------
C     Initialize integral arrays.
C--------------------------------
C
      CALL DZERO(WORK(KINT1T),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT1S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KINT3S),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT4S),NT1AMX*NRHFT*NRHFT)
      CALL DZERO(WORK(KIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KT3AM),NT1AMX*NT1AMX*NT1AMX)
      CALL DZERO(WORK(KT3BAR),NT1AMX*NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME1),NT1AMX)
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME11),NT1AMX)
      CALL DZERO(WORK(KOME22),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KVIRT),NVIRT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KOCC),NRHFT*NRHFT*NRHFT*NVIRT)
C
C----------------------------------------------
C     Read MO-coefficients from interface file.
C----------------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
C
      READ (LUSIFC)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDS)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT  CHANGED NLAMDT TO NLAMDS
COMMENT COMMENT COMMENT
C
C
      CALL CMO_REORDER(WORK(KCMO),WORK(KEND1),LWRK1)
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBTS)
COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT CHANGED NORBT TO NORBTS
COMMENT COMMENT COMMENT COMMENT
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      IF (FROIMP .OR. FROEXP) THEN
         CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
      ENDIF
C
      CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
C
C--------------------------------------------------
C     Read in AO fock matrix and transform to MO
C--------------------------------------------------
C
      LUFCK = -1
C     The T1 transformed density is used to calculate this AO fock matrix.
C      CALL GPOPEN(LUFCK,'CC_FCKH','OLD',' ','UNFORMATTED',
C     *            IDUMMY,.FALSE.)
C     The CMO transformed density is used to calculate this AO fock matrix.
      CALL GPOPEN(LUFCK,'CC_FCKREF','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
C
      REWIND(LUFCK)
      READ(LUFCK)(WORK(KFOCK + I-1),I = 1,N2BST(ISYMOP))
      CALL GPCLOSE(LUFCK,'KEEP' )
C
      IF (IPRINT .GT. 50) THEN
         CALL AROUND( 'FOP3 : Usual Fock AO matrix' )
         CALL CC_PRFCKAO(WORK(KFOCK),ISYMOP)
      ENDIF
C
      CALL CC_FCKMO(WORK(KFOCK),WORK(KCMO),WORK(KCMO),
     *              WORK(KEND1),LWRK1,1,1,1)
C
COMMENT
C      IF (IPRINT .GT. 50) THEN
COMMENT
         CALL AROUND('OUTPUT OF FOCK')
         CALL OUTPUT(WORK(KFOCK),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
COMMENT
C      ENDIF
COMMENT
C
C----------------------------------------------
C      Calculate the one electron integrals
C      from the field. Add this contri
C      to the fock matrix as well.
C----------------------------------------------
C
      IF ((NONHF) .AND. NFIELD .GT. 0) THEN
C
         KONEP = KEND1
         KEND1 = KONEP + N2BST(ISYMOP)
         LWRK1 = LWORK - KEND1
C
         CALL DZERO(WORK(KONEP),N2BST(ISYMOP))
C
         DO I = 1, NFIELD
C
            KONEP2 = KEND1
            KEND2  = KONEP2 + N2BST(ISYMOP)
            LWRK2  = LWORK - KEND2
C
            IF (LWRK2 .LT. N2BST(ISYMOP)) THEN
               CALL QUIT('Exceeded workmemory in CC_FOP3 (Field)')
            ENDIF
C
            CALL DZERO(WORK(KONEP2),N2BST(ISYMOP))
            FF = EFIELD(I)
            CALL CC_ONEP(WORK(KONEP2),WORK(KEND2),LWRK2,FF,1,LFIELD(I))
C
            CALL DAXPY(N2BST(ISYMOP),1.0D0,WORK(KONEP2),1,WORK(KONEP),1)
         ENDDO
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND( 'FOP3 : Usual AO matrix from field(s)' )
            CALL CC_PRFCKAO(WORK(KONEP),ISYMOP)
         ENDIF
C
         CALL CC_FCKMO(WORK(KONEP),WORK(KCMO),WORK(KCMO),
     *                 WORK(KEND1),LWRK1,1,1,1)
C
COMMENT
C         IF (IPRINT .GT. 50) THEN
COMMENT
            CALL AROUND('OUTPUT OF STRENGTH ONEP')
            CALL OUTPUT(WORK(KONEP),1,NORBT,1,NORBT,NORBT,NORBT,
     *                  1,LUPRI)
COMMENT
C         ENDIF
COMMENT
C
         CALL DAXPY(N2BST(ISYMOP),1.0D0,WORK(KONEP),1,WORK(KFOCK),1)
COMMENT
C         IF (IPRINT .GT. 50) THEN
COMMENT
            CALL AROUND('OUTPUT OF FOCK+STRENGTH ONEP')
            CALL OUTPUT(WORK(KFOCK),1,NORBT,1,NORBT,NORBT,NORBT,
     *                  1,LUPRI)
COMMENT
C         ENDIF
COMMENT
      ENDIF
C
C---------------------------------------------------------
C      Calculate the one electron integrals
C      from the field (no strength).
C---------------------------------------------------------
C
      IF (((NONHF) .AND. (NFIELD .GT. 0)) .OR. (DIPMOM)) THEN
C
         KONEP2 = KEND1
         KEND1  = KONEP2 + N2BST(ISYMOP)
         LWRK1  = LWORK - KEND1
C
         CALL DZERO(WORK(KONEP2),N2BST(ISYMOP))
C
         IF (DIPMOM) THEN
            FF = 1.0D0
            ISYONE = -1
            CALL CC_ONEP(WORK(KONEP2),WORK(KEND1),LWRK1,FF,
     *                   ISYONE,'ZDIPLEN ')
         ELSE
            DO I = 1, NFIELD
C
               KONEP3 = KEND1
               KEND2  = KONEP3 + N2BST(ISYMOP)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. N2BST(ISYMOP)) THEN
                  CALL QUIT('Exceeded workmemory in CC_FOP3 (Field)')
               ENDIF
C
               CALL DZERO(WORK(KONEP3),N2BST(ISYMOP))
               FF = 1.0D0
               ISYONE = 1
               CALL CC_ONEP(WORK(KONEP3),WORK(KEND2),LWRK2,FF,
     *                      ISYONE,LFIELD(I))
C
               CALL DAXPY(N2BST(ISYMOP),1.0D0,WORK(KONEP3),1,
     *                    WORK(KONEP2),1)
            ENDDO
         ENDIF
C
         IF (IPRINT .GT. 50) THEN
            CALL AROUND( 'FOP3 : Usual AO matrix from field(s)' )
            CALL CC_PRFCKAO(WORK(KONEP2),ISYMOP)
         ENDIF
C
         CALL CC_FCKMO(WORK(KONEP2),WORK(KCMO),WORK(KCMO),
     *                 WORK(KEND1),LWRK1,1,1,1)
C
COMMENT
C         IF (IPRINT .GT. 50) THEN
COMMENT
            CALL AROUND('OUTPUT OF ONEP WITHOUT STREGTH')
            CALL OUTPUT(WORK(KONEP2),1,NORBT,1,NORBT,NORBT,NORBT,
     *                  1,LUPRI)
COMMENT
C         ENDIF
COMMENT
C
      ENDIF
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      IF (DIRECT) THEN
         NTOSYM = 1
         DTIME  = SECOND()
c        CALL HERDI1(WORK(KEND1),LWRK1,IPRINT)
         CALL QUIT('Direct not implemented')
         DTIME  = SECOND() - DTIME
         TIMHER1 = TIMHER1 + DTIME
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            NTOT = MAXSHL
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C-----------------------------------------------------------------
C           If direct calculate the integrals and transposed t2am.
C-----------------------------------------------------------------
C
            IF (DIRECT) THEN
               DTIME  = SECOND()
c              CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,IPRINT)
               CALL QUIT('Direct not implemented')
               DTIME  = SECOND() - DTIME
               TIMHER2 = TIMHER2 + DTIME
COMMENT COMMENT
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CC_FOP3')
               END IF
COMMENT COMMENT
            ELSE
               NUMDIS = 1
            ENDIF
C
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
C------------------------------------------
C              Work space allocation no. 2.
C------------------------------------------
C
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient space in CCSD_TRIPLE')
               ENDIF
C
C
C-----------------------------------------
C              Read in batch of integrals.
C-----------------------------------------
C
               DTIME   = SECOND()
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
               DTIME   = SECOND() - DTIME
               TIMRDAO = TIMRDAO  + DTIME
C
C----------------------------------------------------
C              Calculate integrals needed in CCSD[T].
C----------------------------------------------------
C
C             CALL CCSDT_TRAN1(WORK(KINT1),WORK(KINT2),WORK(KCMO),
C     *                          WORK(KCMO),WORK(KXINT),IDEL)
C
COMMENT COMMENT
COMMENT COMMENT
               CALL CCSDT_TRAN1(WORK(KINT1T),WORK(KINT2T),WORK(KCMO),
     *                          WORK(KCMO),WORK(KXINT),IDEL)
C               CALL CCSDT_TRAN1(WORK(KINT1T),WORK(KINT2T),XLAMDP,
C     *                          XLAMDH,WORK(KXINT),IDEL)
C
               CALL CCSDT_TRAN2(WORK(KIAJB),WORK(KYIAJB),WORK(KCMO),
     *                          WORK(KXINT),IDEL)
C
               CALL CCSDT_TRAN3(WORK(KINT1S),WORK(KINT2S),XLAMDP,
     *                          XLAMDH,WORK(KXINT),IDEL)
COMMENT COMMENT
COMMENT COMMENT
               CALL CCSDT_TRAN3PT(WORK(KINT3S),WORK(KINT4S),WORK(KCMO),
     *                            WORK(KXINT),IDEL)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
      call print_integrals(work(kint3s),work(kint4s))
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C------------------------------------------------
C     Calculate the triple cluster amplitudes.
C------------------------------------------------
C
      rho1n = ddot(nt1amx,t1am,1,t1am,1)
      write(lupri,*) 'Norm of T1AM before T3AM   : ',rho1n
      rho1n = ddot(nt1amx*nt1amx,t2am,1,t2am,1)
      write(lupri,*) 'Norm of T2AM before T3AM   : ',rho1n
      rho1n = ddot(nt1amx*nt1amx*nt1amx,work(kt3am),1,work(kt3am),1)
      write(lupri,*) 'Norm of T3AM before T3AM   : ',rho1n
      rho1n = ddot(nt1amx*nvirt*nvirt,work(kint3s),1,work(kint3s),1)
      write(lupri,*) 'Norm of KINT3S before T3AM : ',rho1n
      rho1n = ddot(nt1amx*nrhft*nrhft,work(kint4s),1,work(kint4s),1)
      write(lupri,*) 'Norm of KINT4S before T3AM : ',rho1n
      rho1n = ddot(nt1amx,work(kscr1),1,work(kscr1),1)
      write(lupri,*) 'Norm of KSCR1  before call : ',rho1n
      write(lupri,*) 'Before CCSDT_T3AM'
      call flshfo(lupri)
C
      IF ((NFIELD .GT. 0) .AND. (NONHF)) THEN
         ITERATIVE = .TRUE.
         CALL CCSDT_T3AM(WORK(KT3AM),WORK(KINT3S),WORK(KINT4S),T2AM,
     *                WORK(KSCR1),WORK(KFOCKD),WORK(KT3BAR),WORK(KONEP),
     *                ITERATIVE)
      ELSE
         ITERATIVE = .FALSE.
         CALL CCSDT_T3AM(WORK(KT3AM),WORK(KINT3S),WORK(KINT4S),T2AM,
     *                   WORK(KSCR1),WORK(KFOCKD),DUMMY,DUMMY,ITERATIVE)
      ENDIF
C
      call temppr(work(kt3am),1)
C
      rho1n = ddot(nt1amx*nt1amx*nt1amx,work(kt3am),1,work(kt3am),1)
      write(lupri,*) 'Norm of T3AM after T3AM   : ',rho1n
C
C-----------------------------------------------------
C     For reference calculate the CCSD(T) energy.
C-----------------------------------------------------
C
      CALL AROUND('START OF CCSD(T) ENERGY')
C
COMMENT COMMENT
COMMENT COMMENT
      rho1n = ddot(nt1amx*nt1amx,work(kiajb),1,work(kiajb),1)
      write(lupri,*) 'Norm of IAJB integrals     = ',rho1n
      rho1n = ddot(nt1amx*nvirt*nvirt,work(kint1t),1,work(kint1t),1)
      write(lupri,*) 'Norm of INT1T integrals    = ',rho1n
      rho1n = ddot(nt1amx*nrhft*nrhft,work(kint2t),1,work(kint2t),1)
      write(lupri,*) 'Norm of INT2T integrals    = ',rho1n
      rho1n = ddot(n2bst(isymop),work(kfock),1,work(kfock),1)
      write(lupri,*) 'Norm of FOCK iintermediate = ',rho1n
      rho1n = ddot(nt1amx*nt1amx*nt1amx,work(kt3am),1,work(kt3am),1)
      write(lupri,*) 'Norm of T3AM before calc.  = ',rho1n
      xnorm = ddot(nt1amx,work(kome1),1,work(kome1),1)
      write(lupri,*) 'Omega1 norm before calc.   = ',xnorm
      xnorm = ddot(nt1amx*nt1amx,work(kome2),1,work(kome2),1)
      write(lupri,*) 'Omega2 norm before calc.   = ',xnorm
COMMENT COMMENT
COMMENT COMMENT
      CALL CCSDT_OMEGA1OLD(WORK(KOME1),WORK(KIAJB),WORK(KT3AM))
C
      CALL CCSDT_OMEGA2OLD(WORK(KOME2),WORK(KINT1T),WORK(KINT2T),
     *                     WORK(KFOCK),WORK(KT3AM))
C
COMMENT COMMENT
COMMENT COMMENT
      xnorm = ddot(nt1amx,work(kome1),1,work(kome1),1)
      write(lupri,*) 'Omega1 norm after  calc. = ',xnorm
      xnorm = ddot(nt1amx*nt1amx,work(kome2),1,work(kome2),1)
      write(lupri,*) 'Omega2 norm after  calc. = ',xnorm
      xnorm = ddot(nt1amx,t1am,1,t1am,1)
      write(lupri,*) 'T1 norm                  = ',xnorm
      xnorm = ddot(nt1amx*nt1amx,t2am,1,t2am,1)
      write(lupri,*) 'T2 norm                  = ',xnorm
COMMENT COMMENT
COMMENT COMMENT
      ENERT1 = ddot(nt1amx,t1am,1,work(kome1),1)
      write(LUPRI,*) 'The E5(ST) contribution :',2.0D0*ENERT1
C
      CALL OUTPUT(WORK(KOME2),1,NT1AMX,1,NT1AMX,NT1AMX,NT1AMX,1,LUPRI)
      CALL ECCSDT(T2AM,WORK(KOME2))
C
      CALL AROUND('END OF CCSD(T) ENERGY')
C
C---------------------------------------------
C     Calculate the triple bar amplitudes.
C---------------------------------------------
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
      rho1n = ddot(nt1amx*nt1amx*nt1amx,work(kt3bar),1,work(kt3bar),1)
      write(lupri,*) 'Norm of T3BAR before call  : ',rho1n
      rho1n = ddot(nt1amx*nvirt*nvirt,work(kint3s),1,work(kint3s),1)
      write(lupri,*) 'Norm of KINT3S before call : ',rho1n
      rho1n = ddot(nt1amx*nrhft*nrhft,work(kint4s),1,work(kint4s),1)
      write(lupri,*) 'Norm of KINT4S before call : ',rho1n
      rho1n = ddot(nt1amx*nt1amx,work(kiajb),1,work(kiajb),1)
      write(lupri,*) 'Norm of KIAJB before call  : ',rho1n
      rho1n = ddot(n2bst(isymop),work(kfock),1,work(kfock),1)
      write(lupri,*) 'Norm of KFOCK before call  : ',rho1n
      rho1n = ddot(nt1amx,t1am,1,t1am,1)
      write(lupri,*) 'Norm of T1AM before call   : ',rho1n
      rho1n = ddot(nt1amx*nt1amx,t2am,1,t2am,1)
      write(lupri,*) 'Norm of T2AM before call   : ',rho1n
      rho1n = ddot(norbt,work(kfockd),1,work(kfockd),1)
      write(lupri,*) 'Norm of KFOCKD before call : ',rho1n
      rho1n = ddot(nt1amx,work(kscr1),1,work(kscr1),1)
      write(lupri,*) 'Norm of KSCR1  before call : ',rho1n
      call flshfo(lupri)
C
C           CALL AROUND( 'In CC_FOP3: Triples Fock MO matrix' )
C           CALL CC_PRFCKMO(WORK(KFOCK),ISYMOP)
C      do a = 1, nvirt
C      do i = 1, nrhft
C          nai = nvirt*(I-1) + a
C          if (abs(t1am(nai)) .gt. 1.0D-11) then
C             write(lupri,*) 'T1AM^{',a,'}_{',i,'} = ',t1am(nai)
C          endif
C      enddo
C      enddo
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
C
C
      CALL CCSDT_T3BAR(WORK(KT3BAR),WORK(KINT3S),WORK(KINT4S),
     *                 WORK(KIAJB),WORK(KFOCK),T1AM,T2AM,
     *                 WORK(KSCR1),WORK(KFOCKD))
C
      call temppr(work(kt3bar),2)
C
      rho1n = ddot(nt1amx*nt1amx*nt1amx,work(kt3am),1,work(kt3am),1)
      write(lupri,*) 'Norm of T3AM after T3BAR   = ',rho1n
      rho1n = ddot(nt1amx*nt1amx*nt1amx,work(kt3bar),1,work(kt3bar),1)
      write(lupri,*) 'Norm of T3BAR after T3BAR  = ',rho1n
C
C-----------------------------------------------------------
C     Calculate the CCSD(T) corrections to the density.
C-----------------------------------------------------------
C
      CALL DZERO(WORK(KOME1),NT1AMX)
C ia block
      CALL CCSDT_OMEGA1(WORK(KOME1),T2AM,WORK(KT3AM))
C
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
      call around('1e- density in NODDY')
      call cc_prp(work(kome1),DUMMY,1,1,0)
COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT
      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      ENDDO
C------------------------------------------------------
C     Calculate the term to the unrelaxed gradient.
C------------------------------------------------------
C
COMMENT COMMENT
      write(lupri,*) 'Field before contraction (fop3)'
      call output(work(konep2),1,norbt,1,norbt,norbt,norbt,1,lupri)
COMMENT COMMENT
      CALL CCSDT_UNRE(T2AM,WORK(KT3AM),WORK(KT3BAR),WORK(KOME1),
     *                WORK(KONEP2))
C
C ij and ab block
C      CALL CCSDT_NONREL_TERM(DENAB,DENIJ,WORK(KT3AM),WORK(KT3BAR))
C-----------------------------------------------------------
C     For debugging purposes calculate the right hand side
C     of the equations for the kappa amplitudes.
C-----------------------------------------------------------
C
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
C
      rho1n = ddot(nt1amx*nt1amx*nt1amx,work(kt3am),1,work(kt3am),1)
      write(lupri,*) 'Norm of T3AM before kappa  : ',rho1n
      rho1n = ddot(nt1amx*nt1amx*nt1amx,work(kt3bar),1,work(kt3bar),1)
      write(lupri,*) 'Norm of T3BAR before kappa : ',rho1n
      rho1n = ddot(nt1amx,t1am,1,t1am,1)
      write(lupri,*) 'Norm of T1AM before kappa  : ',rho1n
      rho1n = ddot(nt1amx*nt1amx,t2am,1,t2am,1)
      write(lupri,*) 'Norm of T2AM before kappa  : ',rho1n
      rho1n = ddot(nvirt*nvirt*nvirt*nrhft,work(kvirt),1,work(kvirt),1)
      write(lupri,*) 'Norm of VIRT before call   : ',rho1n
      rho1n = ddot(nvirt*nrhft*nrhft*nrhft,work(kocc),1,work(kocc),1)
      write(lupri,*) 'Norm of OCC  before kappa  : ',rho1n
      rho1n = ddot(nt1amx*nt1amx,work(kome2),1,work(kome2),1)
      write(lupri,*) 'Norm of OMEGA2 before kappa: ',rho1n
C
      CALL CCSDT_KAPPA(WORK(KT3AM),WORK(KT3BAR),T1AM,T2AM,
     *                 WORK(KVIRT),WORK(KOCC),WORK(KOME2),
     *                 WORK(KOME22))
C
C-----------------------------------------------------------
C     For debugging purposes calculate the right hand side
C     of the equations for the t1-bar and t2-bar
C     amplitudes.
C-----------------------------------------------------------
C
C      CALL CCSDT_OMEGA1OLD(WORK(KOME11),WORK(KIAJB),WORK(KT3AM))
C
COMMENT COMMENT
COMMENT COMMENT
C      do i = 1, nt1amx
C         if (abs(work(kome11-1+i)) .gt. 1.0D-9) then
C            write(lupri,*) 'Omega1_eta_debug(',i,') = ',work(kome11-1+i)
C         endif
C      enddo
COMMENT COMMENT
COMMENT COMMENT
C
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KOME22),NT1AMX*NT1AMX)
      CALL CCETA_OMEGA2(WORK(KOME22),WORK(KOME2),WORK(KT3AM),
     *                  WORK(KT3BAR),WORK(KINT1T),WORK(KINT2T),
     *                  WORK(KFOCK))
C
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
C       call around('OMEGA2_ETA in NODDY :')
C       call output(work(kome22),1,nt1amx,1,nt1amx,nt1amx,nt1amx,1,lupri)
COMMENT COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT COMMENT
C-----------------------------------------------------------
C     For debugging purposes calculate the diagonal
C     of the kappa amplitudes
C-----------------------------------------------------------
C
      CALL DZERO(WORK(KKAPAA),NVIRT)
      CALL DZERO(WORK(KKAPII),NRHFT)
      CALL DEBUG_KAPPADIAG(WORK(KKAPAA),WORK(KKAPII),
     *                     WORK(KT3AM),WORK(KT3BAR))
C
      DO I = 1, NVIRT
C         IF (ABS(WORK(KKAPAA+I-1)) .GT. 1.0D-9) THEN
            WRITE(LUPRI,*) 'Kappa_{aa}(',i,') = ',WORK(KKAPAA+I-1)
C         ENDIF
      ENDDO
      DO I = 1, NRHFT
C         IF (ABS(WORK(KKAPII+I-1)) .GT. 1.0D-9) THEN
            WRITE(LUPRI,*) 'Kappa_{ii}(',i,') = ',WORK(KKAPII+I-1)
C         ENDIF
      ENDDO
C
      RETURN
      END
C  /* Deck ccsdt_tran1 */
      SUBROUTINE CCSDT_TRAN1(XINT1,XINT2,XLAMDP,XLAMDH,AOINT,IDEL)
C
C
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
COMMENT COMMENT
COMMENT COMMENT
      write(lupri,*) 'NRHFT, NVIRT, NT1AMX : ',nrhft,nvirt,nt1amx
      write(lupri,*) 'NORBT, NBAST         : ',norbt,nbast
COMMENT COMMENT
COMMENT COMMENT
C
      DO 100 G = 1,NBAST
         DO 110 IB = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,IB)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
               DO 200 D = 1,NVIRT
                  DO 210 B = 1,NVIRT
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT1(NCK,B,D) = XINT1(NCK,B,D)
     *               + AOINT(NAB,G)*XLAMDH(A,NRHFT+C)*XLAMDP(IB,K)
     *                      *XLAMDP(G,NRHFT+B)*XLAMDH(IDEL,NRHFT+D)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
               DO 300 J = 1,NRHFT
                  DO 310 L = 1,NRHFT
                     DO 320 K = 1,NRHFT
                        DO 330 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT2(NCK,L,J) = XINT2(NCK,L,J)
     *                  + AOINT(NAB,G)*XLAMDH(A,NRHFT+C)*XLAMDP(IB,K)
     *                                *XLAMDP(G,L)*XLAMDH(IDEL,J)
C
  330                   CONTINUE
  320                CONTINUE
  310             CONTINUE
  300          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_tran2 */
      SUBROUTINE CCSDT_TRAN2(XIAJB,YIAJB,CMO,AOINT,IDEL)
C
C
C
#include <implicit.h>
      PARAMETER (TWO = 2.0D0)
      DIMENSION XIAJB(NT1AMX,NT1AMX), AOINT(NNBAST,NBAST)
      DIMENSION YIAJB(NT1AMX,NT1AMX)
      DIMENSION CMO(NBAST,NORBT)
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 G = 1,NBAST
         DO 110 B = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,B)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
               DO 200 L = 1,NRHFT
                  DO 210 D = 1,NVIRT
                     NDL = NVIRT*(L-1) + D
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XIAJB(NCK,NDL) = XIAJB(NCK,NDL)
     *                             + AOINT(NAB,G)*
     *        (TWO*CMO(A,NRHFT+C)*CMO(B,K)*CMO(G,NRHFT+D)*CMO(IDEL,L)
     *         - CMO(A,NRHFT+C)*CMO(B,L)*CMO(G,NRHFT+D)*CMO(IDEL,K))
C
                           YIAJB(NCK,NDL) = YIAJB(NCK,NDL)
     *                             + AOINT(NAB,G)*
     *           CMO(A,NRHFT+C)*CMO(B,K)*CMO(G,NRHFT+D)*CMO(IDEL,L)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_t3am */
      SUBROUTINE CCSDT_T3AM(T3AM,XINT1,XINT2,T2AM,SCR1,FOCKD,
     *                      T3AM2,FCKFIELD,ITERATIVE)
C
C
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX),SCR1(NT1AMX),FOCKD(*)
      DIMENSION T3AM2(NT1AMX,NT1AMX,NT1AMX)
      DIMENSION T2AM(NT1AMX,NT1AMX)
      DIMENSION FCKFIELD(NORBT,NORBT)
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccfield.h>
C
      PARAMETER(ONE = 1.0D0, HALF = 0.5D0)
C
      LOGICAL LDBG1, LDBG2, FIRST, NONCONV, ITERATIVE
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 50 I = 1,NRHFT
         DO 60 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
            SCR1(NAI) = FOCKD(NRHFT+A) - FOCKD(I)
   60    CONTINUE
   50 CONTINUE
C
      IF (ITERATIVE) THEN
         CALL DZERO(T3AM2,NT1AMX*NT1AMX*NT1AMX)
      ENDIF
C
      FIRST  = .TRUE.
      NONCONV = .TRUE.
      MAXITE  = 0
C
COMMENT COMMENT
      write(lupri,*) 'ITERATIVE = ',iterative
COMMENT COMMENT
      DO WHILE ((FIRST) .OR. 
     *          ((ITERATIVE) .AND. (NONCONV)))
C
         FIRST  = .FALSE.
C
         DO NCK = 1,NT1AMX
C
            DO J = 1,NRHFT
               DO B = 1,NVIRT
C
                  NBJ = NVIRT*(J-1) + B
C
                  DO NAI = 1,NT1AMX
C
                     AIBJCK = 0.0D0
                     DO D = 1,NVIRT
C
                        NDJ = NVIRT*(J-1) + D
C
                        AIBJCK = AIBJCK + XINT1(NCK,B,D)*T2AM(NDJ,NAI)
C
                     ENDDO
C
                     DO L = 1,NRHFT
C
                        NBL = NVIRT*(L-1) + B
C
                        AIBJCK = AIBJCK - XINT2(NCK,L,J)*T2AM(NBL,NAI)
COMMENT COMMENT
COMMENT COMMENT
C          if (abs(XINT2(NCK,L,J)*T2AM(NBL,NAI)) .gt. 1.0D-12) then
C              write(lupri,'(A,5I3)') 'NAI, B, J, NCK, L : ',
C     *                                NAI, B, J, NCK, L
C              write(lupri,*) 'INT(NCK.L.J) T2AM(NBL,NAI)  : ',
C     *                      XINT2(NCK,L,J),T2AM(NBL,NAI)
C          endif
COMMENT COMMENT
COMMENT COMMENT
C
                     ENDDO
C
C---------------------------------------------------------------------
C       Add field dependent terms
C       N^8 scaling but who cares in a noddy!!!
C       Divide the first multiplication in two to get N^7 scaling
C---------------------------------------------------------------------
C
                     IF (ITERATIVE) THEN
C
                        DO D = 1, NVIRT
                           NDJ = NVIRT*(J-1) + D
                           DO L = 1, NRHFT
                              NBL = NVIRT*(L-1) + B
COMMENT COMMENT
                              AIBJCK = AIBJCK 
     *                               + T2AM(NAI,NBL)
     *                                *T2AM(NCK,NDJ)
     *                                *FCKFIELD(L,NRHFT+D)
COMMENT COMMENT
                           ENDDO
                        ENDDO
C
                        DO D = 1, NVIRT
                           NDJ = NVIRT*(J-1) + D
COMMENT COMMENT
                           AIBJCK = AIBJCK 
     *                            - HALF*T3AM2(NAI,NDJ,NCK)
     *                                  *FCKFIELD(NRHFT+B,NRHFT+D)
COMMENT COMMENT
                        ENDDO
C
                        DO L = 1, NRHFT
                           NBL = NVIRT*(L-1) + B
COMMENT COMMENT
                           AIBJCK = AIBJCK
     *                            + HALF*T3AM2(NAI,NBL,NCK)
     *                               *FCKFIELD(L,J)
COMMENT COMMENT
                        ENDDO
C
                     ENDIF
C
C-----------------------------------
C       Final summation
C-----------------------------------
C
                     AIBJCK = AIBJCK/(SCR1(NAI) + SCR1(NBJ) + SCR1(NCK))
COMMENT COMMENT
COMMENT COMMENT
C       if (abs(aibjck) .gt. 1.0D-12) then
C           write(lupri,'(A,e20.10)') 'Contri  : ',AIBJCK
C       endif
COMMENT COMMENT
COMMENT COMMENT
C
                     T3AM(NAI,NBJ,NCK) = T3AM(NAI,NBJ,NCK) + AIBJCK
                     T3AM(NAI,NCK,NBJ) = T3AM(NAI,NCK,NBJ) + AIBJCK
                     T3AM(NBJ,NAI,NCK) = T3AM(NBJ,NAI,NCK) + AIBJCK
                     T3AM(NCK,NAI,NBJ) = T3AM(NCK,NAI,NBJ) + AIBJCK
                     T3AM(NBJ,NCK,NAI) = T3AM(NBJ,NCK,NAI) + AIBJCK
                     T3AM(NCK,NBJ,NAI) = T3AM(NCK,NBJ,NAI) + AIBJCK
C
                  ENDDO  ! NAI
               ENDDO     ! B
            ENDDO        ! J
         ENDDO           ! NCK
C
C----------------------------------------
C     Remove the forbidden elements.
C----------------------------------------
C
         do a = 1, nvirt
         do i = 1, nrhft
            nai = nvirt*(i-1) + a
            do b = 1, nvirt
               nbi = nvirt*(i-1) + b
               do c = 1, nvirt
                  nci = nvirt*(i-1) + c
C
                  t3am(nai,nbi,nci) = 0.0D0
                  t3am(nai,nci,nbi) = 0.0D0
                  t3am(nbi,nai,nci) = 0.0D0
                  t3am(nbi,nci,nai) = 0.0D0
                  t3am(nci,nai,nbi) = 0.0D0
                  t3am(nci,nbi,nai) = 0.0D0
               enddo
            enddo
C
            do j = 1, nrhft
               naj = nvirt*(j-1) + a
               do k = 1, nrhft
                  nak = nvirt*(k-1) + a
C
                  t3am(nai,naj,nak) = 0.0D0
                  t3am(nai,nak,naj) = 0.0D0
                  t3am(naj,nai,nak) = 0.0D0
                  t3am(naj,nak,nai) = 0.0D0
                  t3am(nak,nai,naj) = 0.0D0
                  t3am(nak,naj,nai) = 0.0D0
               enddo
            enddo
C
         enddo
         enddo
C
C------------------------------------------------------------
C     If iterative field thing compare the new and the old
C     T3 amplitudes and see if they differ.
C------------------------------------------------------------
C
        IF (ITERATIVE) THEN
C
           CALL DAXPY(NT1AMX*NT1AMX*NT1AMX,-ONE,T3AM,1,T3AM2,1)
C
           XNORM1 = DDOT(NT1AMX*NT1AMX*NT1AMX,T3AM2,1,T3AM2,1)
C
           XNORM = SQRT(XNORM1)
COMMENT COMMENT
         write(lupri,*) 'Norm of difference in T3AM = ',XNORM
COMMENT COMMENT
C
           IF (XNORM .LT. 1.0D-10) THEN
              NONCONV = .FALSE.
           ELSE
              CALL DCOPY(NT1AMX*NT1AMX*NT1AMX,T3AM,1,T3AM2,1)
              CALL DZERO(T3AM,NT1AMX*NT1AMX*NT1AMX)
              MAXITE = MAXITE + 1
              IF (MAXITE .GT. 30) THEN
                 CALL QUIT('MAX ITERATIONS EXCEEDED in CCSDT_T3AM')
              ENDIF
           ENDIF
        ENDIF
C
C--------------------------
C     End do while loop
C--------------------------
C
      ENDDO
C
C
      RETURN
      END
C  /* Deck ccsdt_omega1 */
      SUBROUTINE CCSDT_OMEGA1(OMEGA1,T2AM,T3AM)
C
C
C
#include <implicit.h>
      PARAMETER (TWO = 2.0D0)
      DIMENSION OMEGA1(NT1AMX),T2AM(NT1AMX,NT1AMX)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO 115 B = 1,NVIRT
            DO 120 J = 1,NRHFT
               NBJ = NVIRT*(J-1) + B
C
               DO 130 K = 1,NRHFT
                  NAK = NVIRT*(K-1) + A
                  NBK = NVIRT*(K-1) + B
                  DO 140 C = 1,NVIRT
                     NCK = NVIRT*(K-1) + C
                     NCI = NVIRT*(I-1) + C
                     NCJ = NVIRT*(J-1) + C
C
                     XTEMP = TWO*
     *               (TWO*T2AM(NCK,NBJ)-T2AM(NCJ,NBK))*
     *               ( T3AM(NAI,NBJ,NCK) - T3AM(NAK,NBJ,NCI) )
C
                     OMEGA1(NAI) = OMEGA1(NAI) + XTEMP
C
COMMENT COMMENT
COMMENT COMMENT
C        if (abs(xtemp) .gt. 1.0D-9) then
C           write(lupri,*) 'Cont to NAI (',nai,') = ',xtemp
C           if (i .EQ. 1 .AND. a .EQ. 3) then
C               write(lupri,*) 'A, B, C, I, J, K = ',a,b,c,i,j,k
C               write(lupri,*) 'NAI, NBJ, NCK = ',nai,nbj,nck
C               write(lupri,*) 'NAK, NCI = ',nak,nci
C               write(lupri,*) 'TWO*T2AM(CK,BJ) = ',TWO*T2AM(NCK,NBJ)
C               write(lupri,*) 'T2AM(CJ,BK) = ',T2AM(NCJ,NBK)
C               write(lupri,*) 'T3AM(AI,BJ,CK) = ',T3AM(NAI,NBJ,NCK) 
C               write(lupri,*) 'T3AM(AK,BJ,CI) = ', T3AM(NAK,NBJ,NCI)
C           endif
C        endif
COMMENT COMMENT
COMMENT COMMENT
C
C
  140             CONTINUE
  130          CONTINUE
C
  120       CONTINUE
  115       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_tran3 */
      SUBROUTINE CCSDT_TRAN3(XINT1,XINT2,XLAMDP,XLAMDH,AOINT,IDEL)
C
C
C
#include <implicit.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION XLAMDP(NBAST,NORBT), XLAMDH(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      LOGICAL LDEBUG
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      LDEBUG = .FALSE.
C
      DO 100 G = 1,NBAST
         DO 110 IB = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,IB)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
               DO 200 D = 1,NVIRT
                  DO 210 B = 1,NVIRT
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT1(NCK,B,D) = XINT1(NCK,B,D)
     *               + AOINT(NAB,G)*XLAMDP(A,NRHFT+C)*XLAMDH(IB,K)
     *                      *XLAMDP(G,NRHFT+B)*XLAMDH(IDEL,NRHFT+D)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
               DO 300 J = 1,NRHFT
                  DO 310 L = 1,NRHFT
                     DO 320 K = 1,NRHFT
                        DO 330 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT2(NCK,L,J) = XINT2(NCK,L,J)
     *                  + AOINT(NAB,G)*XLAMDP(A,NRHFT+C)*XLAMDH(IB,K)
     *                                *XLAMDP(G,L)*XLAMDH(IDEL,J)
C
  330                   CONTINUE
  320                CONTINUE
  310             CONTINUE
  300          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF (LDEBUG) THEN
         do a = 1, nt1amx
            do b = 1, nrhtf
               do c = 1, nrhtf
                 if (abs(xint1(a,b,c)) .gt. 1.0D-9) then
                    write(lupri,*) 'Lam-int1(',a,',',b,',',c,') = ',
     *                   xint1(a,b,c)
                 endif
               enddo
            enddo
         enddo
         do a = 1, nt1amx
            do b = 1, nvirt
               do c = 1, nvirt
                 if (abs(xint2(a,b,c)) .gt. 1.0D-9) then
                    write(lupri,*) 'Lam-int2(',a,',',b,',',c,') = ',
     *                   xint2(a,b,c)
                 endif
               enddo
            enddo
         enddo
      ENDIF
C
      RETURN
      END
C  /* Deck ccsdt_tran3pt */
      SUBROUTINE CCSDT_TRAN3PT(XINT1,XINT2,CMO,AOINT,IDEL)
C
C
C
#include <implicit.h>
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT), XINT2(NT1AMX,NRHFT,NRHFT)
      DIMENSION CMO(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 G = 1,NBAST
         DO 110 IB = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,IB)
C
               if (aoint(nab,g) .eq. 0.0d0) goto 120
               DO 200 D = 1,NVIRT
                  DO 210 B = 1,NVIRT
                     DO 220 K = 1,NRHFT
                        DO 230 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT1(NCK,B,D) = XINT1(NCK,B,D)
     *               + AOINT(NAB,G)*CMO(A,NRHFT+C)*CMO(IB,K)
     *                      *CMO(G,NRHFT+B)*CMO(IDEL,NRHFT+D)
C
  230                   CONTINUE
  220                CONTINUE
  210             CONTINUE
  200          CONTINUE
C
               DO 300 J = 1,NRHFT
                  DO 310 L = 1,NRHFT
                     DO 320 K = 1,NRHFT
                        DO 330 C = 1,NVIRT
C
                           NCK = NVIRT*(K-1) + C
C
                           XINT2(NCK,L,J) = XINT2(NCK,L,J)
     *                  + AOINT(NAB,G)*CMO(A,NRHFT+C)*CMO(IB,K)
     *                                *CMO(G,L)*CMO(IDEL,J)
C
  330                   CONTINUE
  320                CONTINUE
  310             CONTINUE
  300          CONTINUE
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_t3bar */
      SUBROUTINE CCSDT_T3BAR(T3BAR,XINT1,XINT2,XIAJB,FOCK,T1AM,
     *                       T2AM,SCR1,FOCKD)
C
C
C
#include <implicit.h>
C
      PARAMETER (TWO = 2.0D0)
C
      DIMENSION XINT1(NT1AMX, NVIRT,NVIRT), XINT2(NT1AMX, NRHFT,NRHFT)
      DIMENSION XIAJB(NT1AMX, NT1AMX), FOCK(NORBT,NORBT)
      DIMENSION T3BAR(NT1AMX, NT1AMX,NT1AMX), SCR1(NT1AMX), FOCKD(*)
      DIMENSION T1AM(NT1AMX), T2AM(NT1AMX,NT1AMX)
      LOGICAL LDBG1, LDBG2
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 50 NI = 1,NRHFT
         DO 60 NA = 1,NVIRT
            NAI = NVIRT*(NI-1) + NA
            SCR1(NAI) = FOCKD(NRHFT+NA) - FOCKD(NI)
   60    CONTINUE
   50 CONTINUE
COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT
C       do a = 1, nvirt
C       do i = 1, nrhft
C          nai = nvirt*(i-1) + a
C          do b = 1, nvirt
C          do j = 1, nvirt
C             nbj = nvirt*(j-1) + b
CC
C             if (abs(t2am(nai,nbj)) .gt. 1.0D-10) then
C                write(lupri,*) 'T2AM^{',a,',',b,'}_{',i,',',j,'}) = ',
C     *                       t2am(nai,nbj)
C             endif
CC
C          enddo
C          enddo
C       enddo
C       enddo
C
C       do ni = 1, nt1amx
C           if (abs(scr1(ni)) .gt. 1.0D-9) then
C              write(lupri,*) 'SCR1(',ni,') = ',SCR1(ni)
C           endif
C       enddo
C
C       do na = 1, nvirt
C          do ni = 1, nrhft
C             do nb = 1, nvirt
C             do nc = 1, nvirt
CC
C             nai = NVIRT*(NI-1) + NA
C             nbi = NVIRT*(NI-1) + NB
CC
C             if (abs(TWO*XINT1(nai,nb,nc)-XINT1(nbi,na,nc)) 
C     *              .gt. 1.0D-9) then
C                write(lupri,*)'L_vir(',nai,',',nb,',',nc,') =',
C     *                 TWO*XINT1(nai,nb,nc)-XINT1(nbi,na,nc)
C             endif
C             enddo
C             enddo
C          enddo
C       enddo
C       do na = 1, nvirt
C          do ni = 1, nrhft
C             do nb = 1, nvirt
C             do nc = 1, nvirt
CC
C             nai = NVIRT*(NI-1) + NA
C             nbi = NVIRT*(NI-1) + NB
CC
C             if (abs(XINT1(nai,nb,nc)) .gt. 1.0D-9) then
C                write(lupri,*)'g_vir(',nai,',',nb,',',nc,') =',
C     *                 XINT1(nai,nb,nc)
C             endif
C             enddo
C             enddo
C          enddo
C       enddo
CC
CC
C       xnorm = 0.0D0
C       do na = 1, nvirt
C          do ni = 1, nrhft
C             do nj = 1, nrhft
C             do nk = 1, nrhft
CC
C             nai = NVIRT*(NI-1) + NA
C             nak = NVIRT*(NK-1) + NA
CC
C             xnorm = xnorm + (TWO*XINT2(nai,nj,nk)-XINT2(nak,nj,ni))
C     *                      *(TWO*XINT2(nai,nj,nk)-XINT2(nak,nj,ni))
CC
C             if (abs(TWO*XINT2(nai,nj,nk)-XINT2(nak,nj,ni)) 
C     *              .gt. 1.0D-9) then
C                write(lupri,*)'L_occ(',nai,',',nj,',',nk,') =',
C     *                 TWO*XINT2(nai,nj,nk)-XINT2(nak,nj,ni)
C             endif
C             enddo
C             enddo
C          enddo
C       enddo
CC
C       write(lupri,*) 'Norm squared of the L_occ : ',xnorm
C       do na = 1, nvirt
C          do ni = 1, nrhft
C             do nj = 1, nrhft
C             do nk = 1, nrhft
CC
C             nai = NVIRT*(NI-1) + NA
CC
C             if (abs(XINT2(nai,nj,nk)) .gt. 1.0D-9) then
C                write(lupri,*)'g_occ(',nai,',',nj,',',nk,') =',
C     *                 XINT2(nai,nj,nk)
C             endif
C             enddo
C             enddo
C          enddo
C       enddo
COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT
COMMENT COMMENT COMMENT COMMENT
C
      DO 100 NK = 1,NRHFT
      DO 105 NC = 1,NVIRT
C
         NCK = NVIRT*(NK-1) + NC
C
         DO 110 NJ = 1,NRHFT
            NCJ = NVIRT*(NJ-1) + NC
            DO 120 NB = 1,NVIRT
C
               NBJ = NVIRT*(NJ-1) + NB
               NBK = NVIRT*(NK-1) + NB
C
               DO 130 NI = 1,NRHFT
                  NBI = NVIRT*(NI-1) + NB
                  NCI = NVIRT*(NI-1) + NC
               DO 135 NA = 1,NVIRT
C
                  NAI = NVIRT*(NI-1) + NA
                  NAJ = NVIRT*(NJ-1) + NA
                  NAK = NVIRT*(NK-1) + NA
C
                  AIBJCK = 0.0D0
C
          if (.true.) then
C     T1* = TWO T1 => Factor of two
                  AIBJCK = AIBJCK + TWO*XIAJB(NBJ,NCK)*T1AM(NAI)
COMMENT COMMENT
COMMENT COMMENT
C            if (abs(XIAJB(NBJ,NCK)*T1AM(NAI)) .gt. 1.0D-10) then
C               write(lupri,*) 'Contribution from 1.st term : ',
C     *                             XIAJB(NBJ,NCK)*T1AM(NAI)
C               write(lupri,*) 'With NAI, NBJ, NCK : ',NAI,NBJ,NCK
C            endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
          if (.true.) then
C     T1* = TWO T1 => Factor of two
                  AIBJCK = AIBJCK - TWO*XIAJB(NBJ,NCI)*T1AM(NAK)
COMMENT COMMENT
COMMENT COMMENT
C             if (abs(XIAJB(NBJ,NCK)*T1AM(NAI)) .gt. 1.0D-10) then
C                write(lupri,*) 'Contribution from 2.nd term : ',
C     *                              XIAJB(NBJ,NCK)*T1AM(NAI)
C                write(lupri,*) 'With NAI, NBJ, NCK : ',NAI,NBJ,NCK
C           endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
          if (.true.) then
                  AIBJCK = AIBJCK + TWO*FOCK(NK,NRHFT+NC)*
     *                           (TWO*T2AM(NAI,NBJ)-T2AM(NAJ,NBI))
          endif
C
          if (.true.) then
                  AIBJCK = AIBJCK - TWO*FOCK(NJ,NRHFT+NC)*
     *                           (TWO*T2AM(NAI,NBK)-T2AM(NAK,NBI))
          endif
C
                  DO 140 ND = 1,NVIRT
C
                     NDI = NVIRT*(NI-1) + ND
                     NDJ = NVIRT*(NJ-1) + ND
                     NDK = NVIRT*(NK-1) + ND
C
          if (.true.) then
                     AIBJCK = AIBJCK + TWO*
     *                        (TWO*XINT1(NCK,NB,ND)-XINT1(NBK,NC,ND))*
     *                        (TWO*T2AM(NDJ,NAI)-T2AM(NDI,NAJ))
C
COMMENT COMMENT
COMMENT COMMENT
C              if (abs((TWO*XINT1(NCK,NB,ND)-XINT1(NBK,NC,ND))*
C     *                        (TWO*T2AM(NDJ,NAI)-T2AM(NDI,NAJ)))
C     *                          .gt. 1.0D-10) then
C                   write(lupri,*) 'Contribution from 5.th term :',
C     *                    TWO*(TWO*XINT1(NCK,NB,ND)-XINT1(NBK,NC,ND))*
C     *                        (TWO*T2AM(NDJ,NAI)-T2AM(NDI,NAJ))
C                   write(lupri,*) '2*XINT1, XINT1, 2*T2AM, T2AM',
C     *                    TWO*XINT1(NCK,NB,ND),XINT1(NBK,NC,ND),
C     *                        TWO*T2AM(NDJ,NAI),T2AM(NDI,NAJ)
C              endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
          if (.true.) then
                     AIBJCK = AIBJCK - TWO*
     *                        XINT1(NBI,NC,ND)*
     *                        (TWO*T2AM(NDK,NAJ)-T2AM(NDJ,NAK))
COMMENT COMMENT
COMMENT COMMENT
C              if (abs(TWO*XINT1(NBI,NC,ND)*
C     *                        (TWO*T2AM(NDK,NAJ)-T2AM(NDJ,NAK)))
C     *                          .gt. 1.0D-10) then
C                   write(lupri,*) 'Contribution from 6.th term :',
C     *                        TWO*XINT1(NBI,NC,ND)*
C     *                        (TWO*T2AM(NDK,NAJ)-T2AM(NDJ,NAK))
C              endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
  140             CONTINUE
C
                  DO 150 NL = 1,NRHFT
C
                     NAL = NVIRT*(NL-1) + NA
                     NBL = NVIRT*(NL-1) + NB
C
          if (.true.) then
                     AIBJCK = AIBJCK - TWO*
     *                        (TWO*XINT2(NCK,NL,NJ)-XINT2(NCJ,NL,NK))*
     *                        (TWO*T2AM(NBL,NAI)-T2AM(NBI,NAL))
COMMENT COMMENT
COMMENT COMMENT
C         if (abs(TWO*
C     *           (TWO*XINT2(NCK,NL,NJ)-XINT2(NCJ,NL,NK))*
C     *           (TWO*T2AM(NBL,NAI)-T2AM(NBI,NAL))) .gt. 1.0D-10) then
CC
C            write(lupri,*) 'Contribution from L_occ term :',TWO*
C     *                 (TWO*XINT2(NCK,NL,NJ)-XINT2(NCJ,NL,NK))*
C     *                 (TWO*T2AM(NBL,NAI)-T2AM(NBI,NAL))
C            write(lupri,*) 'T2TCME, L : ',
C     *                 (TWO*T2AM(NBL,NAI)-T2AM(NBI,NAL)),
C     *                 (TWO*XINT2(NCK,NL,NJ)-XINT2(NCJ,NL,NK))
C         endif
COMMENT COMMENT
COMMENT COMMENT
          endif
C
          if (.true.) then
                     AIBJCK = AIBJCK + TWO*
     *                        XINT2(NCJ,NL,NI)*
     *                        (TWO*T2AM(NBK,NAL)-T2AM(NBL,NAK))
          endif
C
  150             CONTINUE
C
C
                  AIBJCK = AIBJCK/(SCR1(NAI) + SCR1(NBJ) + SCR1(NCK))
C
COMMENT COMMENT
COMMENT COMMENT
C       if (abs(aibjck).gt.(1.0D-10/(scr1(nai)+scr1(nbj)+scr1(nck)))) 
C     *                                                              then
C          write(lupri,'(1X,A,I3,I3,I3,I3,I3,I3)') 
C     *                  'Contribution for A,B,C,I,J,K = ',
C     *                                 NA,NB,NC,NI,NJ,NK
C          write(lupri,'(A,3I3,e20.10)') 'NAI,NBJ,NCK,AIBJCK = ',
C     *                                   NAI,NBJ,NCK,AIBJCK
CC          write(lupri,*) 'SCR1(NAI),SCR1(NBJ),SCR1(NCK), SUM = ',
CC     *                    SCR1(NAI),SCR1(NBJ),SCR1(NCK),
CC     *                    SCR1(NAI)+SCR1(NBJ)+SCR1(NCK)
C          write(lupri,*) '       '
C       endif
COMMENT COMMENT
COMMENT COMMENT
                  T3BAR(NAI,NBJ,NCK) = T3BAR(NAI,NBJ,NCK) + AIBJCK
                  T3BAR(NAI,NCK,NBJ) = T3BAR(NAI,NCK,NBJ) + AIBJCK
                  T3BAR(NBJ,NAI,NCK) = T3BAR(NBJ,NAI,NCK) + AIBJCK
                  T3BAR(NCK,NAI,NBJ) = T3BAR(NCK,NAI,NBJ) + AIBJCK
                  T3BAR(NBJ,NCK,NAI) = T3BAR(NBJ,NCK,NAI) + AIBJCK
                  T3BAR(NCK,NBJ,NAI) = T3BAR(NCK,NBJ,NAI) + AIBJCK
C
  135          CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  105 CONTINUE
  100 CONTINUE
C
C------------------------------------------------------
C     Get rid of amplitudes that are not allowed.
C------------------------------------------------------
C
      do a = 1, nvirt
      do i = 1, nrhft
         nai = nvirt*(i-1) + a
         do b = 1, nvirt
            nbi = nvirt*(i-1) + b
            do c = 1, nvirt
               nci = nvirt*(i-1) + c
C
               t3bar(nai,nbi,nci) = 0.0D0
               t3bar(nai,nci,nbi) = 0.0D0
               t3bar(nbi,nai,nci) = 0.0D0
               t3bar(nbi,nci,nai) = 0.0D0
               t3bar(nci,nai,nbi) = 0.0D0
               t3bar(nci,nbi,nai) = 0.0D0
            enddo
         enddo
C
         do j = 1, nrhft
            naj = nvirt*(j-1) + a
            do k = 1, nrhft
               nak = nvirt*(k-1) + a
C
               t3bar(nai,naj,nak) = 0.0D0
               t3bar(nai,nak,naj) = 0.0D0
               t3bar(naj,nai,nak) = 0.0D0
               t3bar(naj,nak,nai) = 0.0D0
               t3bar(nak,nai,naj) = 0.0D0
               t3bar(nak,naj,nai) = 0.0D0
            enddo
         enddo
C
      enddo
      enddo
C
      RETURN
      END
C  /* Deck temppr */
      SUBROUTINE TEMPPR(T3AM,IOPT)
C
#include <implicit.h>
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL AROUND('ENTERED TEMPPR')
      WRITE(LUPRI,*) 'IOPT  = ',IOPT
C
      DO NA = 1,NVIRT
      DO NB = 1,NVIRT
      DO NC = 1,NVIRT
         DO NI = 1,NRHFT
         NAI = NVIRT*(NI-1) + NA
            DO NJ = 1,NRHFT
               NBJ = NVIRT*(NJ-1) + NB
               DO NK = 1,NRHFT
                  NCK = NVIRT*(NK-1) + NC
C
               if (abs(t3am(nai,nbj,nck)) .gt. 1.0D-11) then
                 if (iopt .eq. 1) then
           write(lupri,1) 't3am (^{',na,',',nb,',',nc,'}_{',ni,',',
     *                                nj,',',nk,'}) -> (',
     *                      NAI,',',NBJ,',',NCK,') = ',t3am(nai,nbj,nck)
                 else if (iopt .eq. 2) then
           write(lupri,1) 't3bar(^{',na,',',nb,',',nc,'}_{',ni,',',
     *                                nj,',',nk,'}) -> (',
     *                      NAI,',',NBJ,',',NCK,') = ',t3am(nai,nbj,nck)
                 else
           call quit('Wrong IOPT in temppr')
                 endif
               endif
C
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ENDDO
      ENDDO
C
   1  FORMAT(1x,A8,I3,A1,I3,A1,I3,A3,I3,A1,I3,A1,I3,A7,I3,A1,I3,A1,I3,
     *       A4,E20.10)
C
      RETURN
      END
C  /* Deck ccsdt_nonrelterm */
      SUBROUTINE CCSDT_NONREL_TERM(DENAB,DENIJ,T3AM,T3BAR)
C
C
C
#include <implicit.h>
      PARAMETER (TWO = 2.0D0, THREE = 3.0D0, HALF = 0.5D0)
      DIMENSION DENAB(NVIRT,NVIRT), DENIJ(NRHFT,NRHFT)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
      DIMENSION T3BAR(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO 115 B = 1,NVIRT
            DO 120 J = 1,NRHFT
               NBJ = NVIRT*(J-1) + B
C
               DO 130 K = 1,NRHFT
                  NAK = NVIRT*(K-1) + A
                  NBK = NVIRT*(K-1) + B
                  DO 140 C = 1,NVIRT
                     NCK = NVIRT*(K-1) + C
                     NCI = NVIRT*(I-1) + C
                     NCJ = NVIRT*(J-1) + C
C
                     DO 150 L = 1, NRHFT
C
                        NCL = NVIRT*(L-1) + C
C
                        DENIJ(K,L) = DENIJ(K,L) 
     *                             - HALF*T3BAR(NAI,NBJ,NCL)*
     *                                      T3AM(NAI,NBJ,NCK)
C
  150                CONTINUE
C
                     DO 160 D = 1, NVIRT
C
                        NDK = NVIRT*(K-1) + D
C
                        DENAB(C,D) = DENAB(C,D)
     *                             + HALF*T3BAR(NAI,NBJ,NCK)*
     *                                      T3AM(NAI,NBJ,NDK)
C
  160                CONTINUE
C
  140             CONTINUE
  130          CONTINUE
C
  120       CONTINUE
  115       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE ECCSDT(T2AM,OMEGA2)
C
C
C
#include <implicit.h>
      PARAMETER (TWO=2.0d0)
      DIMENSION T2AM(NT1AMX,NT1AMX),OMEGA2(NT1AMX,NT1AMX)
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ENERGY = 0.0d0
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
            DO 120 J = 1,NRHFT
               NAJ = NVIRT*(J-1) + A
               DO 130 B = 1,NVIRT
                  NBJ = NVIRT*(J-1) + B
                  NBI = NVIRT*(I-1) + B
                  DELTA = OMEGA2(NAI,NBJ)
     *                   *(TWO*T2AM(NAI,NBJ) - T2AM(NAJ,NBI))
                  ENERGY = ENERGY + DELTA
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      WRITE(LUPRI,*) 'Energy contribution E4(TT) :',ENERGY
C
      RETURN
      END
C  /* Deck ccsdt_omega2 */
      SUBROUTINE CCSDT_OMEGA2OLD(OMEGA2,XINT1T,XINT2T,FOCK,T3AM)
C
C
C
#include <implicit.h>
      PARAMETER (TWO = 2.0D0)
      DIMENSION XINT1T(NT1AMX,NVIRT,NVIRT)
      DIMENSION XINT2T(NT1AMX,NRHFT,NRHFT)
      DIMENSION OMEGA2(NT1AMX,NT1AMX)
      DIMENSION FOCK(NORBT,NORBT)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO 120 J = 1,NRHFT
               DO 130 B = 1,NVIRT
                  NBJ = NVIRT*(J-1) + B
C
                  DO 140 K = 1,NRHFT
                     NBK = NVIRT*(K-1) + B
                     NAK = NVIRT*(K-1) + A
                     DO 150 C = 1,NVIRT
C
                        NCK = NVIRT*(K-1) + C
                        NCJ = NVIRT*(J-1) + C
                        NCI = NVIRT*(I-1) + C
C
                        OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) -
     *    (T3AM(NAI,NBJ,NCK)-T3AM(NAK,NBJ,NCI))*FOCK(K,NRHFT+C)
C
                        DO 160 D = 1,NVIRT
C
                           NDJ = NVIRT*(J-1) + D
                           NDK = NVIRT*(K-1) + D
                           NDI = NVIRT*(I-1) + D
C
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) +
     *  (T3AM(NBK,NCI,NDJ) - TWO*T3AM(NBJ,NCI,NDK) + T3AM(NBJ,NCK,NDI))
     *   *XINT1T(NDK,A,C)
C
  160                   CONTINUE
C
                        DO 170 L = 1,NRHFT
C
                           NBL = NVIRT*(L-1) + B
                           NCL = NVIRT*(L-1) + C
                           NAL = NVIRT*(L-1) + A
C
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) -
     *  (T3AM(NBL,NAK,NCJ) - TWO*T3AM(NBJ,NAK,NCL) + T3AM(NBJ,NAL,NCK))
     *   *XINT2T(NCL,K,I)
C
  170                   CONTINUE
C
  150                CONTINUE
  140             CONTINUE
C
  130          CONTINUE
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      DO 200 NAI = 1,NT1AMX
         DO 210 NBJ = 1,NAI
C
            XAIBJ = OMEGA2(NAI,NBJ) + OMEGA2(NBJ,NAI)
            OMEGA2(NAI,NBJ) = XAIBJ
            OMEGA2(NBJ,NAI) = XAIBJ
C
  210    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck ccsdt_omega1old */
      SUBROUTINE CCSDT_OMEGA1OLD(OMEGA1,XIAJB,T3AM)
C
C
C
#include <implicit.h>
      PARAMETER (TWO = 2.0D0)
      DIMENSION OMEGA1(NT1AMX),XIAJB(NT1AMX,NT1AMX)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
#include <priunit.h>
COMMENT COMMENT COMMENT
C#include <inforb.h>
#include <ccorb.h>
COMMENT COMMENT COMMENT
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      DO 100 I = 1,NRHFT
         DO 110 A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO 120 NBJ = 1,NT1AMX
C
               DO 130 K = 1,NRHFT
                  NAK = NVIRT*(K-1) + A
                  DO 140 C = 1,NVIRT
                     NCK = NVIRT*(K-1) + C
                     NCI = NVIRT*(I-1) + C
                     OMEGA1(NAI) = OMEGA1(NAI) - XIAJB(NCK,NBJ)*
     *               ( T3AM(NAI,NBJ,NCK) - T3AM(NAK,NBJ,NCI) )
  140             CONTINUE
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck cceta_omega2 */
      SUBROUTINE CCETA_OMEGA2(OMEGA2,OMEGA22,T3AM,T3BAR,XINT1T,
     *                        XINT2T,FOCK)
C
C     Kasper Hald, Fall 2001.
C
C     Calculate the doubles of the right hand side of the t1-bar and t2-bar
C     equation.
C     NOTE : OMEGA22 will be erased.
C
#include <implicit.h>
C
      INTEGER INDEX, NAI, NBJ, NBK, NAK, NCK, NCJ, NCI, NDJ, NDK, NDI
      INTEGER NBL, NCL, NAL, NAJ, NBI, AI
C
#if defined (SYS_CRAY)
      REAL OMEGA2(NT1AMX,NT1AMX)
      REAL OMEGA22(NT1AMX,NT1AMX)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL T3BAR(NT1AMX,NT1AMX,NT1AMX)
      REAL XINT1T(NT1AMX,NVIRT,NVIRT)
      REAL XINT2T(NT1AMX,NRHFT,NRHFT)
      REAL FOCK(NORBT,NORBT)
      REAL ONE, TWO, XAIBJ
#else
      DOUBLE PRECISION OMEGA2(NT1AMX,NT1AMX)
      DOUBLE PRECISION OMEGA22(NT1AMX,NT1AMX)
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION T3BAR(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION XINT1T(NT1AMX,NVIRT,NVIRT)
      DOUBLE PRECISION XINT2T(NT1AMX,NRHFT,NRHFT)
      DOUBLE PRECISION FOCK(NORBT,NORBT)
      DOUBLE PRECISION ONE, TWO, XAIBJ
#endif
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
C
#include <priunit.h>
C#include <inforb.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C      DO AI = 1, NT1AMX
C         DO A = 1, NVIRT
C         DO B = 1, NVIRT
C            if (abs(xint1t(ai,a,b)) .gt. 1.0D-9) then
C               write(lupri,*) 'int1t(',ai,',',a,',',b,') = ',
C     *              xint1t(ai,a,b)
C            endif
C         ENDDO
C         ENDDO
C      ENDDO
C      DO AI = 1, NT1AMX
C         DO I = 1, NRHFT
C         DO J = 1, NRHFT
C            if (abs(xint2t(ai,i,j)) .gt. 1.0D-9) then
C               write(lupri,*) 'int2t(',ai,',',i,',',j,') = ',
C     *              xint2t(ai,i,j)
C            endif
C         ENDDO
C         ENDDO
C      ENDDO
C
C---------------------------------------------
C     Start by calculating the T3AM terms
C     and sum up.
C---------------------------------------------
C
      DO I = 1,NRHFT
         DO A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO J = 1,NRHFT
               DO B = 1,NVIRT
                  NBJ = NVIRT*(J-1) + B
C
                  DO K = 1,NRHFT
                     NBK = NVIRT*(K-1) + B
                     NAK = NVIRT*(K-1) + A
                     DO C = 1,NVIRT
C
                        NCK = NVIRT*(K-1) + C
                        NCJ = NVIRT*(J-1) + C
                        NCI = NVIRT*(I-1) + C
C
                if (.true.) then
                   OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) +
     *                (
     *                 T3AM(NBJ,NAI,NCK)
     *                -T3AM(NBJ,NAK,NCI)
     *                )*FOCK(K,NRHFT+C)
C
                endif
C
                        DO D = 1,NVIRT
C
                           NDJ = NVIRT*(J-1) + D
                           NDK = NVIRT*(K-1) + D
                           NDI = NVIRT*(I-1) + D
C
                if (.false.) then
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) -
     *  (T3AM(NBK,NCI,NDJ) - TWO*T3AM(NBJ,NCI,NDK) + T3AM(NBJ,NCK,NDI))
     *   *XINT1T(NDK,A,C)
C
                endif
C
                        ENDDO
C
                        DO L = 1,NRHFT
C
                           NBL = NVIRT*(L-1) + B
                           NCL = NVIRT*(L-1) + C
                           NAL = NVIRT*(L-1) + A
C
                 if (.false.) then
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) +
     *  (T3AM(NBL,NAK,NCJ) - TWO*T3AM(NBJ,NAK,NCL) + T3AM(NBJ,NAL,NCK))
     *   *XINT2T(NCL,K,I)
C
                 endif
C
                        ENDDO
C
                     ENDDO
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C------------------------------------------------
C     Make P^{ab}_{ij} and P^{2CME}_{aibj} for
C     the T3AM contributions.
C------------------------------------------------
C
      DO NAI = 1,NT1AMX
         DO NBJ = 1,NAI
C
            XAIBJ = OMEGA2(NAI,NBJ) + OMEGA2(NBJ,NAI)
            OMEGA2(NAI,NBJ) = XAIBJ
            OMEGA2(NBJ,NAI) = XAIBJ
C
         ENDDO
      ENDDO
C
      CALL DZERO(OMEGA22,NT1AMX*NT1AMX)
C
      DO A = 1,NVIRT
      DO I = 1,NRHFT
         DO B = 1,NVIRT
         DO J = 1,NRHFT
C
            NAI = NVIRT*(I-1) + A
            NAJ = NVIRT*(J-1) + A
            NBI = NVIRT*(I-1) + B
            NBJ = NVIRT*(J-1) + B
            OMEGA22(NAI,NBJ) = OMEGA22(NAI,NBJ) + TWO*OMEGA2(NAI,NBJ)
            OMEGA22(NAJ,NBI) = OMEGA22(NAJ,NBI) - OMEGA2(NAI,NBJ)
C
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
C---------------------------------------------
C     Now calculate the T3BAR terms
C---------------------------------------------
C
      CALL DZERO(OMEGA2,NT1AMX*NT1AMX)
C
      DO I = 1,NRHFT
         DO A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO J = 1,NRHFT
               DO B = 1,NVIRT
                  NBJ = NVIRT*(J-1) + B
C
                  DO K = 1,NRHFT
                     NBK = NVIRT*(K-1) + B
                     NAK = NVIRT*(K-1) + A
                     DO C = 1,NVIRT
C
                        NCK = NVIRT*(K-1) + C
                        NCJ = NVIRT*(J-1) + C
                        NCI = NVIRT*(I-1) + C
C
                        DO D = 1,NVIRT
C
                           NDJ = NVIRT*(J-1) + D
                           NDK = NVIRT*(K-1) + D
                           NDI = NVIRT*(I-1) + D
C
                if (.true.) then
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) +
     *              T3BAR(NBJ,NCI,NDK)*XINT1T(NDK,A,C)
C
                endif
C
                        ENDDO
C
                        DO L = 1,NRHFT
C
                           NBL = NVIRT*(L-1) + B
                           NCL = NVIRT*(L-1) + C
                           NAL = NVIRT*(L-1) + A
C
                 if (.true.) then
                           OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI) +
     *                     - T3BAR(NBJ,NAK,NCL)*XINT2T(NCL,K,I)
C
                 endif
C
                        ENDDO
C
                     ENDDO
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------------------------------------
C     Make P^{ab}_{ij} for the T3BAR contributions.
C----------------------------------------------------
C
      DO NAI = 1,NT1AMX
         DO NBJ = 1,NAI
C
            XAIBJ = OMEGA2(NAI,NBJ) + OMEGA2(NBJ,NAI)
            OMEGA2(NAI,NBJ) = XAIBJ
            OMEGA2(NBJ,NAI) = XAIBJ
C
         ENDDO
      ENDDO
C
C-------------------------------
C     Sum to final result.
C-------------------------------
C
      CALL DAXPY(NT1AMX*NT1AMX,ONE,OMEGA22,1,OMEGA2,1)
C
      RETURN
      END
C  /* Deck ccsdt_kappa */
      SUBROUTINE CCSDT_KAPPA(T3AM,T3BAR,T1AM,T2AM,VIRT,OCC,
     *                       OMEGA2,OMEGA22)
C
C     Kasper Hald, Fall 2001.
C
C     Calculate the doubles of the right hand side of the Kappa
C     equation.
C
#include <implicit.h>
C
      INTEGER INDEX, NAI, NBJ, NBK, NAK, NCK, NCJ, NCI, NDJ, NDK, NDI
      INTEGER NBL, NCL, NAL, NAJ, NBI, AI
C
#if defined (SYS_CRAY)
      REAL T3AM(NT1AMX,NT1AMX,NT1AMX)
      REAL T3BAR(NT1AMX,NT1AMX,NT1AMX)
      REAL T1AM(NT1AMX)
      REAL T2AM(NT1AMX,NT1AMX)
      REAL VIRT(NVIRT,NVIRT,NRHFT,NVIRT)
      REAL OCC(NRHFT,NVIRT,NRHFT,NRHFT)
      REAL OMEGA2(NT1AMX,NT1AMX)
      REAL OMEGA22(NT1AMX,NT1AMX)
      REAL ONE, TWO, XAIBJ
#else
      DOUBLE PRECISION T3AM(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION T3BAR(NT1AMX,NT1AMX,NT1AMX)
      DOUBLE PRECISION T1AM(NT1AMX)
      DOUBLE PRECISION T2AM(NT1AMX,NT1AMX)
      DOUBLE PRECISION VIRT(NVIRT,NVIRT,NRHFT,NVIRT)
      DOUBLE PRECISION OCC(NRHFT,NVIRT,NRHFT,NRHFT)
      DOUBLE PRECISION OMEGA2(NT1AMX,NT1AMX)
      DOUBLE PRECISION OMEGA22(NT1AMX,NT1AMX)
      DOUBLE PRECISION ONE, TWO, XAIBJ
#endif
C
      LOGICAL L3ABIC, L3IAJB, L3IAJK1, L3IAJK2, L3IAJK3
      LOGICAL LBARABCI, LBARAIJK, LDEBUG
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
C
#include <priunit.h>
C#include <inforb.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
C-------------------------------------------------------
C     Specify which terms that should be calculated.
C-------------------------------------------------------
C
      LDEBUG    = .false.
      L3ABIC    = .true.
      L3IAJB    = .true.
      L3IAJK1   = .true.
      L3IAJK2   = .true.
      L3IAJK3   = .true.
      LBARABCI  = .true.
      LBARAIJK  = .true.
C
C---------------------------------------------
C     Start by calculating the T3AM terms
C     and sum up.
C---------------------------------------------
C
      DO I = 1,NRHFT
         DO A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO J = 1,NRHFT
               DO B = 1,NVIRT
                  NBI = NVIRT*(I-1) + B
                  NBJ = NVIRT*(J-1) + B
                  NAJ = NVIRT*(J-1) + A
C
                  DO K = 1,NRHFT
                     NBK = NVIRT*(K-1) + B
                     NAK = NVIRT*(K-1) + A
                     DO C = 1,NVIRT
C
                        NCK = NVIRT*(K-1) + C
                        NCJ = NVIRT*(J-1) + C
                        NCI = NVIRT*(I-1) + C
C
                   if (L3IAJB) then
                        OMEGA2(NBJ,NAI) = OMEGA2(NBJ,NAI)  +
     * ( T3AM(NCK,NAI,NBJ) -T3AM(NCI,NAK,NBJ) )*T1AM(NCK)
C
COMMENT COMMENT
COMMENT COMMENT
C       if (abs(T3AM(NCI,NAK,NBJ)*T1AM(NCK)) .gt. 1.0D-11) then
C          write(lupri,'(A,6I3)') 'Cont. from A,B,C,I,J,K : ',a,b,c,i,j,k
C          write(lupri,*) 'Term is equal to ',
C     *               -TWO*T3AM(NCI,NAK,NBJ)*T1AM(NCK)
C          write(lupri,*) 'With T1 and T3 : ',T1AM(NCK),T3AM(NCI,NAK,NBJ)
C       endif
COMMENT COMMENT
COMMENT COMMENT
                   endif
C
                        DO D = 1,NVIRT
C
                           NDJ = NVIRT*(J-1) + D
                           NDK = NVIRT*(K-1) + D
                           NDI = NVIRT*(I-1) + D
C
                if (L3ABIC) then
                           VIRT(A,B,I,C) = VIRT(A,B,I,C)  - TWO*
     *  (
     *       T3AM(NDK,NCJ,NBI) 
     *  -TWO*T3AM(NBJ,NCI,NDK) 
     *  +    T3AM(NBJ,NCK,NDI)
     *  )
     *   *(TWO*T2AM(NDK,NAJ)-T2AM(NDJ,NAK))
C
COMMENT COMMENT
COMMENT COMMENT
C       if (abs(TWO*
C     * (T3AM(NDK,NCJ,NBI)-TWO*T3AM(NBJ,NCI,NDK)+T3AM(NBJ,NCK,NDI))
C     * (-TWO*T3AM(NBJ,NCI,NDK))
C     *  *(TWO*T2AM(NDK,NAJ)-T2AM(NDJ,NAK))
C     *        ) .gt. 1.0D-12) then
C            write(lupri,*) 'ABIC contribution : ',TWO*
C     *  (T3AM(NDK,NCJ,NBI)-TWO*T3AM(NBJ,NCI,NDK)+T3AM(NBJ,NCK,NDI))
C     *  *(TWO*T2AM(NDK,NAJ)-T2AM(NDJ,NAK))
C            write(lupri,*) 'T3_1             = ',T3AM(NDK,NCJ,NBI)
C            write(lupri,*) 'T3_2             = ',-TWO*T3AM(NBJ,NCI,NDK)
C            write(lupri,*) 'T3_3             = ',T3AM(NBJ,NCK,NDI)
C            write(lupri,*) '2CME OF T2       = ',
C     *              (TWO*T2AM(NDK,NAJ)-T2AM(NDJ,NAK))
C            write(lupri,*) '  '
C       endif
COMMENT COMMENT
COMMENT COMMENT
C
                endif
C
                        ENDDO
C
                        DO L = 1,NRHFT
C
                           NAL = NVIRT*(L-1) + A
                           NBL = NVIRT*(L-1) + B
                           NCL = NVIRT*(L-1) + C
C
                 if (L3IAJK1) then
                           OCC(I,A,J,K) = OCC(I,A,J,K) + TWO*
     *  (
     *  +     T3AM(NBJ,NAL,NCI) 
     *  - TWO*T3AM(NBJ,NAI,NCL) 
     *  +     T3AM(NBI,NAJ,NCL)
     *  )
     *   *(TWO*T2AM(NCL,NBK)-T2AM(NCK,NBL))
C
                 endif
C
                 if (L3IAJK2) then
                           OCC(I,A,J,J) = OCC(I,A,J,J) + TWO*TWO*
     *  (T3AM(NCL,NBK,NAI)-T3AM(NCL,NBI,NAK))*
     *  (TWO*T2AM(NCL,NBK)-T2AM(NCK,NBL))
C
                 endif
C
                 if (L3IAJK3) then
                           OCC(I,A,J,I) = OCC(I,A,J,I) - TWO*
     *  (T3AM(NCL,NBK,NAJ)-T3AM(NCL,NBJ,NAK))*
     *  (TWO*T2AM(NCL,NBK)-T2AM(NCK,NBL))
C
                 endif
C
                        ENDDO
C
                     ENDDO
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------------------------------------------------------
C     Symmetrize the omega2 density with both P^{2CME} and normal P
C----------------------------------------------------------------------
C
      IF (L3IAJB) THEN
         DO NAI = 1, NT1AMX
            DO NBJ = 1, NAI
               XAIBJ = OMEGA2(NAI,NBJ) + OMEGA2(NBJ,NAI)
               OMEGA2(NAI,NBJ) = XAIBJ
               OMEGA2(NBJ,NAI) = XAIBJ
            ENDDO
         ENDDO
C
         CALL DZERO(OMEGA22,NT1AMX*NT1AMX)
C
         DO A = 1,NVIRT
         DO I = 1,NRHFT
            DO B = 1,NVIRT
            DO J = 1,NRHFT
C
               NAI = NVIRT*(I-1) + A
               NAJ = NVIRT*(J-1) + A
               NBI = NVIRT*(I-1) + B
               NBJ = NVIRT*(J-1) + B
               OMEGA22(NAI,NBJ) = OMEGA22(NAI,NBJ) + TWO*OMEGA2(NAI,NBJ)
               OMEGA22(NAJ,NBI) = OMEGA22(NAJ,NBI) - OMEGA2(NAI,NBJ)
C
            ENDDO
            ENDDO
         ENDDO
         ENDDO
      ENDIF
C
C--------------------------------------------
C     Print the calculated terms :
C--------------------------------------------
C
      if (L3ABIC .AND. LDEBUG) then
         do a = 1, nvirt
           do b = 1, nvirt
             do i = 1, nrhft
               do c = 1, nvirt
                  if (abs(virt(a,b,i,c)) .gt. 1.0D-9) then
                     write(lupri,'(A,I3,A,I3,A,I3,A,I3,A,e20.10)')
     *        'd2-3-abic(',a,',',b,',',i,',',c,') = ',
     *               virt(a,b,i,c)
                  endif
               enddo
             enddo
           enddo
         enddo
         write(lupri,*) '        '
      endif
C
      if (L3IAJB .AND. LDEBUG) then
         do a = 1, nvirt
           do i = 1, nrhft
             NAI = NVIRT*(I-1) + A
             do b = 1, nvirt
               do j = 1, nrhft
                  NBJ = NVIRT*(J-1) + B
                  if (abs(omega22(nai,nbj)) .gt. 1.0D-11) then
                     write(lupri,'(A,I3,A,I3,A,I3,A,I3,A,e20.10)')
     *        'd2-3-iajb(',i,',',a,',',j,',',b,') = ',
     *               omega22(nai,nbj)
                  endif
               enddo
             enddo
           enddo
         enddo
          write(lupri,*) '       '
          call around('Two electron density from NODDY')
          call output(omega22,1,nt1amx,1,nt1amx,nt1amx,nt1amx,1,lupri)
          write(lupri,*) '       '
      endif
      
C
      if ((L3IAJK1 .OR. L3IAJK2 .OR. L3IAJK3) .AND. LDEBUG) then
         do i = 1, nrhft
           do a = 1, nvirt
             do j = 1, nrhft
               do k = 1, nrhft
                  if (abs(occ(i,a,j,k)) .gt. 1.0D-9) then
                     write(lupri,'(A,I3,A,I3,A,I3,A,I3,A,e20.10)')
     *        'd2-3-iajk(',i,',',a,',',j,',',k,') = ',
     *               occ(i,a,j,k)
                  endif
               enddo
             enddo
           enddo
         enddo
         write(lupri,*) '        '
      endif
C
C---------------------------------------------
C     Now calculate the T3BAR terms
C---------------------------------------------
C
      CALL DZERO(VIRT,NVIRT*NVIRT*NRHFT*NVIRT)
      CALL DZERO(OCC,NRHFT*NVIRT*NRHFT*NRHFT)
C
      DO I = 1,NRHFT
         DO A = 1,NVIRT
            NAI = NVIRT*(I-1) + A
C
            DO J = 1,NRHFT
               DO B = 1,NVIRT
                  NAJ = NVIRT*(J-1) + A
                  NBJ = NVIRT*(J-1) + B
C
                  DO K = 1,NRHFT
                     NBK = NVIRT*(K-1) + B
                     NAK = NVIRT*(K-1) + A
                     DO C = 1,NVIRT
C
                        NCK = NVIRT*(K-1) + C
                        NCJ = NVIRT*(J-1) + C
                        NCI = NVIRT*(I-1) + C
C
                        DO D = 1,NVIRT
C
                           NDJ = NVIRT*(J-1) + D
                           NDK = NVIRT*(K-1) + D
                           NDI = NVIRT*(I-1) + D
C
                if (LBARABCI) then
C
C       Store it as a,b,i,c though it is a,b,c,i
C
                           VIRT(A,B,I,C) = VIRT(A,B,I,C) -
     *              T3BAR(NAJ,NCI,NDK)*T2AM(NDK,NBJ)
COMMENT COMMENT
COMMENT COMMENT
C       if (abs(T3BAR(NAJ,NCI,NDK)*T2AM(NDK,NBJ)) .gt. 1.0D-12) then
C            write(lupri,*) 'ABCI-bar contribution : ',
C     *         T3BAR(NAJ,NCI,NDK)*T2AM(NDK,NBJ)
C            write(lupri,*) 'T3, T2 : ',T3BAR(NAJ,NCI,NDK),T2AM(NDK,NBJ)
C       endif
COMMENT COMMENT
COMMENT COMMENT
C
                endif
C
                        ENDDO
C
                        DO L = 1,NRHFT
C
                           NBL = NVIRT*(L-1) + B
                           NCL = NVIRT*(L-1) + C
                           NAL = NVIRT*(L-1) + A
C
                 if (LBARAIJK) then
C
C      Store as though it is i,a,j,k though it is a,i,j,k
C
                           OCC(I,A,J,K) = OCC(I,A,J,K) -
     *                       T3BAR(NAI,NBK,NCL)*T2AM(NCL,NBJ)
C
                 endif
C
                        ENDDO
C
                     ENDDO
                  ENDDO
C
               ENDDO
            ENDDO
C
         ENDDO
      ENDDO
C
C----------------------------------------------------
C     Print the contributions from t3bar.
C----------------------------------------------------
C
      if (LBARABCI .AND. LDEBUG) then
         do a = 1, nvirt
           do b = 1, nvirt
             do i = 1, nrhft
               do c = 1, nvirt
                  if (abs(virt(a,b,i,c)) .gt. 1.0D-9) then
                     write(lupri,'(A,I3,A,I3,A,I3,A,I3,A,e20.10)')
     *        'd2-bar-abci(',a,',',b,',',c,',',i,') = ',
     *               virt(a,b,i,c)
                  endif
               enddo
             enddo
           enddo
         enddo
         write(lupri,*) '        '
      endif
C
      if (LBARAIJK .AND. LDEBUG) then
         do a = 1, nvirt
           do i = 1, nrhft
             do j = 1, nrhft
               do k = 1, nrhft
                  if (abs(occ(i,a,j,k)) .gt. 1.0D-9) then
                     write(lupri,'(A,I3,A,I3,A,I3,A,I3,A,e20.10)')
     *        'd2-bar-aijk(',a,',',i,',',j,',',k,') = ',
     *               occ(i,a,j,k)
                  endif
               enddo
             enddo
           enddo
         enddo
         write(lupri,*) '        '
      endif
C
C-------------------------------
C     Sum to final result.
C-------------------------------
C
      CALL DAXPY(NT1AMX*NT1AMX,ONE,OMEGA22,1,OMEGA2,1)
C
      RETURN
      END
C  /* Deck debug_kappadiag */
      SUBROUTINE DEBUG_KAPPADIAG(XKAPAA,XKAPII,T3AM,T3BAR)
C
C
#include <implicit.h>
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      DIMENSION XKAPAA(NVIRT), XKAPII(NRHFT)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
      DIMENSION T3BAR(NT1AMX,NT1AMX,NT1AMX)
C
      PARAMETER (HALF = 0.5D0)
C
C----------------------------------------------
C     Kappa_{aa}
C----------------------------------------------
C
      DO A = 1, NVIRT
         DO I = 1, NRHFT
            NAI = NVIRT*(I-1) + A
            DO NBJ = 1, NT1AMX
               DO NCK = 1, NT1AMX
                 XKAPAA(A) = XKAPAA(A)
     *                     + T3AM(NAI,NBJ,NCK)*T3BAR(NAI,NBJ,NCK)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C----------------------------------------------
C     Kappa_{ii}
C----------------------------------------------
C
      DO I = 1, NRHFT
         DO A = 1, NVIRT
            NAI = NVIRT*(I-1) + A
            DO NBJ = 1, NT1AMX
               DO NCK = 1, NT1AMX
                 XKAPII(I) = XKAPII(I)
     *                     + T3AM(NAI,NBJ,NCK)*T3BAR(NAI,NBJ,NCK)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C------------------------------
C     Scale the kappas
C------------------------------
C
      CALL DSCAL(NVIRT,HALF,XKAPAA,1)
      CALL DSCAL(NRHFT,HALF,XKAPII,1)
C
      RETURN
      END
C  /* Deck print_integrals */
      SUBROUTINE PRINT_INTEGRALS(XINT1,XINT2)
C
#include <implicit.h>
C
      DIMENSION XINT1(NT1AMX,NVIRT,NVIRT)
      DIMENSION XINT2(NT1AMX,NRHFT,NRHFT)
      LOGICAL LDEBUG
C
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      LDEBUG = .FALSE.
C
      IF (LDEBUG) THEN
         DO NAI = 1, NT1AMX
            DO NB = 1, NVIRT
               DO NC = 1, NVIRT
                  IF (ABS(XINT1(NAI,NB,NC)) .GT. 1.0D-9) THEN
                     WRITE(LUPRI,'(A,3I3,A,E20.10)') 
     *                     'INT(NAI,B,C) = INT(',NAI,NB,NC,') = ',
     *                     XINT1(NAI,NB,NC)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         WRITE(LUPRI,*) '             '
         DO NAI = 1, NT1AMX
            DO NL = 1, NVIRT
               DO NK = 1, NVIRT
                  IF (ABS(XINT2(NAI,NL,NK)) .GT. 1.0D-9) THEN
                     WRITE(LUPRI,'(A,3I3,A,E20.10)') 
     *                     'INT(NAI,K,L) = INT(',NAI,NL,NK,') = ',
     *                     XINT2(NAI,NL,NK)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
C
C  /* Deck ccsdt_unre */
      SUBROUTINE CCSDT_UNRE(T2AM,T3AM,T3BAR,DENS,FCKFIELD)
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION T2AM(NT1AMX,NT1AMX)
      DIMENSION T3AM(NT1AMX,NT1AMX,NT1AMX)
      DIMENSION T3BAR(NT1AMX,NT1AMX,NT1AMX)
      DIMENSION DENS(NT1AMX)
      DIMENSION FCKFIELD(NORBT,NORBT)
COMMENT COMMENT
      dimension sumtemp(nt1amx,nrhft,nrhft)
      dimension omeg1(nvirt,nrhft)
COMMENT COMMENT
      PARAMETER (HALF = 0.5D0)
C
COMMENT COMMENT
      call dzero(sumtemp,nt1amx*nrhft*nrhft)
      call dzero(omeg1,nvirt*nrhft)
COMMENT COMMENT
      SUM  = 0.0D0
      SUM1 = 0.0D0
      SUM2 = 0.0D0
      SUM3 = 0.0D0
C
      DO NA = 1, NVIRT
         DO NI = 1, NRHFT
            NAI = NVIRT*(NI-1) + NA
            SUM = SUM + DENS(NAI)*FCKFIELD(NI,NRHFT+NA)
         ENDDO
      ENDDO
C
      DO NCK = 1,NT1AMX
C
         DO J = 1,NRHFT
            DO B = 1,NVIRT
C
               NBJ = NVIRT*(J-1) + B
C
               DO NAI = 1,NT1AMX
C
                  DO D = 1, NVIRT
                     NDJ = NVIRT*(J-1) + D
                     DO L = 1, NRHFT
                        NBL = NVIRT*(L-1) + B
                        SUM1 = SUM1
     *                       - (T2AM(NAI,NBL)
     *                        *T2AM(NCK,NDJ)
     *                        *FCKFIELD(L,NRHFT+D))*T3BAR(NAI,NBJ,NCK)
                     ENDDO
                  ENDDO
C
                  DO D = 1, NVIRT
                     NDJ = NVIRT*(J-1) + D
                     SUM2 = SUM2
     *                    + (HALF*T3AM(NAI,NDJ,NCK)
     *                       *FCKFIELD(NRHFT+B,NRHFT+D)
     *                      )*T3BAR(NAI,NBJ,NCK)
                  ENDDO
C
                  DO L = 1, NRHFT
                     NBL = NVIRT*(L-1) + B
                     SUM3 = SUM3
     *                    - (HALF*T3AM(NAI,NBL,NCK)
     *                       *FCKFIELD(L,J))*T3BAR(NAI,NBJ,NCK)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      CALL AROUND('Contribution to the unrelaxed density :')
      write(lupri,*) 't^{*}*H(1)T3         = ',sum
      write(lupri,*) 't3bar*H(1)*T2*T2     = ',sum1
      write(lupri,*) 't3bar*H(1)*T3 (virt) = ',sum2
      write(lupri,*) 't3bar*H(1)*T3 (occ)  = ',sum3
      write(lupri,*) 'Total contribution   = ',sum+sum1+sum2+sum3
C
      RETURN
      END
C -- end of cc_fop3.F --
