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
      SUBROUTINE LUCITA_ALLOC
*
* Dimensions and
* Allocation of static memory
*
* =====
* Input
* =====
*
* KFREE : Pointer to first element of free space
* Information in /LUCINP/,/ORBINP/,/CSYM/
*
* ======
* Output
* ======
* KFREE : First array of free space after allocation of
*         static memory
*
* =======
* Version
* =======
*
* Modified Jan 1997
*           Fall 97 (KPGINT1 added )

*. Input
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cstate.inc"
#include "csm.inc"
#include "crun.inc"
#include "cprnt.inc"
*.CSMPRD
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)

*./CINTFO/
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*.Output
#include "glbbas.inc"
*.1 : One electron integrals( Complete matrix allocated ) + copy of the original UNMODIFIED 1-elecs ints
      CALL MEMMAN(KINT1, NTOOB ** 2,'ADDS  ',2,'INT1  ')
      CALL MEMMAN(KINT1O,NTOOB ** 2,'ADDS  ',2,'INT1O ')
*.1.1 : Inactive fock matrix
      CALL MEMMAN(KFI  ,NTOOB ** 2,'ADDS  ',2,'FI    ')
      CALL MEMMAN(KFIO ,NTOOB ** 2,'ADDS  ',2,'FIO   ')
!
!     dec 2010 - stefan: reactivate or rather allocate when needed and
!     NOT in general.
!     KINT1O = 0
!     KFIO   = 0
*.2 : Two electron integrals
      IF(NOINT.EQ.0.AND.INCORE.EQ.1) THEN
         CALL MEMMAN(KINT2,NINT2,'ADDS  ',2,'INT2  ')
         CALL DZERO(WORK(KINT2) ,NINT2)
      END IF
*.    initialize
      CALL DZERO(WORK(KINT1) ,NTOOB**2)
      CALL DZERO(WORK(KFI)   ,NTOOB**2)
      CALL DZERO(WORK(KFIO)  ,NTOOB**2)
      CALL DZERO(WORK(KINT1O),NTOOB**2)
*. Pointers to symmetry block of integrals
      CALL MEMMAN(KPINT1,NBINT1,'ADDS  ',2,'PINT1 ')
      CALL MEMMAN(KPINT2,NBINT2,'ADDS  ',2,'PINT2 ')
*. Pointers to nonsymmetric one-electron integrals
      DO ISM = 1, NSMOB
*. triangular packed
        CALL MEMMAN(KPGINT1(ISM),NSMOB,'ADDS  ',2,'PGINT1')
*. no packing
        CALL MEMMAN(KPGINT1A(ISM),NSMOB,'ADDS  ',2,'PGIN1A')
      END DO
*. Symmetry of last index as a function of initial index
      CALL MEMMAN(KLSM1,NBINT1,'ADDS  ',2,'LSM1   ')
      CALL MEMMAN(KLSM2,NBINT2,'ADDS  ',2,'LSM2   ')

!     density matrices
!     ----------------
!     one-body density matrix
      CALL MEMMAN(KRHO1,NTOOB ** 2,'ADDS  ',2,'RHO1  ')
!     two-body density matrix
      LRHO2 = NTOOB**2*(NTOOB**2+1)/2
      if(idensi.eq.2)then 
        CALL MEMMAN(KRHO2,LRHO2,'ADDS  ',2,'RHO2  ')
      end if

      KPNIJ = 1
      KIJKK = 1
      IF(MXP1+MXP2+MXQ. NE. 0 .OR. IPROCC.NE.0 ) THEN
*. Allocate memory for explicit hamiltonian and roots
        NSBDET = MXP1 + MXP2 + MXQ
        NSBDETP = MAX(NSBDET,IPROCC)
        MXP = MXP1 + MXP2
        CALL MEMMAN(KSBEVC,NSBDET*NROOT,'ADDS  ',2,'KSBEVC')
        CALL MEMMAN(KSBEVL,       NROOT,'ADDS  ',2,'KSBEVL')
*. KSBIDT must be able to hold list of subspace dets / print dets
        LSCR = NSBDETP+1
        CALL MEMMAN(KSBIDT,LSCR,'ADDS  ',1,'KSBIDT')
*. ( NSBDET + 1 elements allocated since NSBDET+1 values are
* obtained in order to check for degenerencies )

        CALL MEMMAN(KSBCNF,LSCR        ,'ADDS  ',1,'KSBCNF')
        CALL MEMMAN(KSBIA ,LSCR        ,'ADDS  ',1,'KSBIA ')
        CALL MEMMAN(KSBIB ,LSCR        ,'ADDS  ',1,'KSBIB ')
*
        LH0 = MXP*(MXP+1)/2 + MXP1*MXQ + MXQ
        CALL MEMMAN(KH0,LH0,'ADDS  ',2,'KH0   ')
*. Scratch space for manipulating with H0
        LH0SCR =
     &  MAX(6*NSBDETP,4*NSBDETP+4*NOCOB,
     &      MXP1*(MXP1+1)/2+MXP1**2)
        CALL MEMMAN(KH0SCR,LH0SCR      ,'ADDS  ',2,'KH0SCR')
      ELSE
        KSBEVC = 0
        KSBEVL = 0
        KSBIDT = 0
        KSBCNF = 0
        KH0    = 0
        KH0SCR = 0
      END IF
*. indeces for pair of orbitals symmetry ordered
*
      CALL MEMMAN(KINH1,NTOOB*NTOOB,'ADDS  ',2,'KINH1  ')
*
*. Length of MO-MO and MO_AO expansion file
      LMOMO = 0
      LMOAO = 0
      DO ISM = 1, NSMOB
        LMOMO = LMOMO + NMOS_ENV(ISM)*NMOS_ENV(ISM)
        LMOAO = LMOAO + NMOS_ENV(ISM)*NAOS_ENV(ISM)
      END DO
      CALL MEMMAN(KMOMO  ,LMOMO ,'ADDS  ',2,'MOMO  ')
      IF(NOMOFL.EQ.0) THEN
        WRITE(6,*)
     &  ' Alloc : space for MO-AO matrices will be allocated '
        CALL MEMMAN(KMOAOIN,LMOAO ,'ADDS  ',2,'MOAOIN')
        CALL MEMMAN(KMOAOUT,LMOAO ,'ADDS  ',2,'MOAOUT')
      ELSE
        KMOAOIN = 0
        KMOAOUT = 0
      END IF

!     call lmemchk('LUCITA_ALLOC')

      END
***********************************************************************
      SUBROUTINE LMEMCHK(TEXT)

      CHARACTER*(*) TEXT
#include "priunit.h"
*
* Check memory allocated  with the memory manager
      write (lupri,'(/2A)') 'Memory check - ',TEXT
      CALL MEMMAN(IDUM,IDUM,'CHECK ',IDUM,'IDUM  ')
*
      END
***********************************************************************
      SUBROUTINE MEMMAN(KBASE,KADD,TASK,IR,IDENT)
*
* Memory manager routine
*
* KBASE : New base address
*         If TASK = INI, KBASE is offset for memory to be controlled
*         by MEMMAN
* KADD  : Dimension of array to be added
*         If TASK = INI, KADD is total length of array
* TASK  : = INI  : Initialize                 Character*6
*         = ADDS : Add static memory
*         = ADDL : Add Local memory
*         = FLUSH : Flush local memory
*         = CHECK : Check memory paddings
*         = FREE  : Return first Free word in KBASE
*         = MARK  : Set a mark at current free adress
*         = FLUSM : Flush local memory to previous mark
*         = SFREEM: Return the free memory currently available in real*8 words
* IR    : 1 => integer , 2 => real,
*         ratio between integer and real is RtoI
* IDENT : identifier of memory slice,Character*6
*
* Local Memory not flushed before allocation of additional static memory
* is tranferred to static memory
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KBASE, KFREES, KFREEL, NWORD, IMARK, IBASE, JBASE,
     &          JBASEN, IBASE_L
!               for addressing of WORK
      CHARACTER*6 TASK,IDENT,IIDENT,MARKC
      PARAMETER(NPAD = 1 )
      PARAMETER(MAXLVL = 1024 )
      PARAMETER(MAXMRK = 1024)
      COMMON/CMEMO8/NWORD,KFREES,KFREEL,IBASE(MAXLVL),
     &             IMARK(MAXMRK)
      COMMON/CMEMO4/NS,NL,NM,
     &             LENGTH(MAXLVL),IIDENT(MAXLVL),
     &             MARKL(MAXMRK),MARKS(MAXMRK),MARKC(MAXMRK)
*. Two real*8 words, one added NPAD times before each array, another
*. added NPAD times after each array
      DATA PAD1/0.123456789D0/
      DATA PAD2/0.987654321D0/
*
#include "priunit.h"
#include "mxpdim.inc"
#include "SysDef.inc"
#include "wrkspc.inc"
*
      ISTOP = 0
      ITSOK  = 1
      IPRNTMP = 0
*
      IF(TASK(1:3).EQ.'INI') THEN
*
**.Initialize
*
        NS     = 0
        NL     = 0
        NSNLI  = 0
        KFREES = KBASE
        KFREEL = KBASE
*
        NM        = 1
        MARKC(NM) = '-INI--'
        IMARK(NM) = KFREEL
        MARKL(NM) = NL
        MARKS(NM) = NS
*
*       Compute pointer to end element
*
        NWORD = KBASE + MXPWRD-1
        IPRNTMP = 0
        ISTOP = 0

!       write(lupri,*)'RtoI',RtoI,'for IR',IR
*
** First free word
*
      ELSE IF (TASK(1:4).EQ.'FREE') THEN
       KBASE = KFREEL
       NSNLI = NS+NL
*
**. Static memory
*
      ELSE IF(TASK(1:4).EQ.'ADDS') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREES = KFREEL +(KADD+1)/RtoI + 2*NPAD
        ELSE
          KFREES = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREES-1 .GT. NWORD ) THEN
          WRITE(2,*)
          WRITE(2,*) ' You can''t always get what you want'
          WRITE(2,*) ' No, you can''t always get what you want'
          WRITE(2,*) ' But if you try sometime, you may find '
          WRITE(2,*) ' you get what you need '
          WRITE(2,*) '                       Jagger/Richard '
*
          WRITE(2,*) ' MEMMAN : work array too short '
          WRITE(2,*) ' current and required length ',NWORD,KFREES-1
*
          WRITE(2,*) ' Trying to allocate : identifer,length'
          WRITE(2,'(20X,A,I9)')IDENT,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NS = NS + NL + 1
        NL = 0
        NSNLI = NS + NL
        IF(NS.GT.MAXLVL) THEN
          WRITE(2,*) ' Too many levels in MEMMAN '
          WRITE(2,*) ' Increase MAXLVL from ', MAXLVL
          Call Abend1(11)
        END IF
        IIDENT(NS) = IDENT
        KFREEL     = KFREES
        IBASE(NS)  = KBASE
        DO 10 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD)  = PAD1
          WORK(KFREEL-NPAD-1+IPAD) = PAD2
   10   CONTINUE
*
**. Local memory
*
      ELSE IF(TASK(1:4).EQ.'ADDL') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREEL = KFREEL +(KADD+1)/RtoI + 2*NPAD
        ELSE
          KFREEL = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREEL-1 .GT. NWORD ) THEN
          WRITE(2,*)
          WRITE(2,*) ' You can''t always get what you want'
          WRITE(2,*) ' No, you can''t always get what you want'
          WRITE(2,*) ' But if you try sometime, you may find '
          WRITE(2,*) ' you get what you need '
          WRITE(2,*) '                       Jagger/Richard '

          WRITE(2,*) ' MEMMAN : work array too short '
          WRITE(2,*) ' current and required length ',NWORD,KFREEL-1
          WRITE(2,*) ' Trying to allocate : identifer,offset,  length'
          WRITE(2,'(24X,A,2I9)')IDENT,KBASE,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NL =  NL + 1
        NSNLI = NS+NL
        IF(NS+NL.GT.MAXLVL) THEN
          WRITE(2,*) ' Too many levels in MEMMAN '
          WRITE(2,*) ' Increase MAXLVL from ', MAXLVL
          Call Abend1(11)
        END IF
        IIDENT(NS+NL) = IDENT
        IBASE(NS+NL) = KBASE
        DO 20 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD) = PAD1
          WORK(KFREEL-NPAD-1+IPAD) =  PAD2
   20   CONTINUE
*
** Flush local memory
*
      ELSE IF(TASK(1:5).EQ.'FLUSH') THEN
        NSNLI = NS+NL
        KFREEL = KFREES
        NL = 0
      ELSE IF(TASK(1:4).EQ.'MARK') THEN
*. Set a mark at current free address
        NM = NM + 1
        IF(NM.GT.MAXMRK) THEN
          WRITE(2,*) ' Too many marks  in MEMMAN '
          WRITE(2,*) ' Increase MAXMRK from ', MAXMRK
          Call Abend1(11)
        END IF
        MARKC(NM) = IDENT
        IMARK(NM) = KFREEL
        MARKL(NM) = NL
        MARKS(NM) = NS
        NSNLI = NS + NL
      ELSE IF (TASK(1:5).EQ.'FLUSM') THEN
        NSNLI = NS+NL
        KFREEL = IMARK(NM)
        IF(KFREES.GT.IMARK(NM)) KFREES = IMARK(NM)
        NL = MARKL(NM)
        NS = MARKS(NM)
        NM = NM - 1
      ELSE IF( TASK(1:5).EQ.'CHECK') THEN
        NSNLI = NS+ NL
      ELSE IF( TASK(1:6) .EQ. 'SFREEM') THEN
*. compute current free memory
        KBASE = NWORD - KFREEL - 1
        NSNLI = NS+ NL
      ELSE
          WRITE(2,*) ' MEMMAN : Unknown task parameter ',TASK
          WRITE(2,*) ' Too confused to continue  '
          Call Abend1(11)
      END IF
*
**. Check paddings
*
      ICHECK = 1
      ITSOK = 1
      IF(TASK(1:5).EQ.'CHECK'.OR.ICHECK.EQ.1) THEN
        DO 100 IL = 1, NSNLI
          JBASE = IBASE(IL)
          IF(IL.NE.NSNLI) THEN
           JBASEN = IBASE(IL+1)
          ELSE
           JBASEN = KFREEL + 1
          END IF
          L1OK = 1
          L2OK = 1
          DO  IPAD = 1, NPAD
            IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1 .OR.
     &         WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) THEN
               ITSOK = 0
               WRITE(2,*) ' Memory problem for : '
               WRITE(2,*) '   Level (IL) ',IL
               IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1) L1OK = 0
               IF(WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) L2OK = 0
               IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
                WRITE(2,'(4X,A,I20,4X,A)')
     &          IIDENT(IL),IBASE(IL),'    OKAY     OKAY '
               ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
                WRITE(2,'(4X,A,I20,4X,A)')
     &          IIDENT(IL),IBASE(IL),'    OKAY       -  '
               ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
                WRITE(2,'(4X,A,I20,4X,A)')
     &          IIDENT(IL),IBASE(IL),'     -       OKAY '
               ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
                  WRITE(2,'(4X,A,I20,4X,A)')
     &            IIDENT(IL),IBASE(IL),'     -       OKAY '
               ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
                 WRITE(2,'(4X,A,I20,4X,A)')
     &            IIDENT(IL),IBASE(IL),'     -       -    '
               END IF
            END IF
          END DO
  100   CONTINUE
CSK        iprntmp = 1
      END IF
 1001 CONTINUE
*
        IF(ITSOK.EQ.0.OR.IPRNTMP.NE.0) THEN
          WRITE(2,'(A,A,1X,A)') ' Current task : ', TASK, IDENT
          WRITE(2,*) ' NS, NL, NSNLI',NS,NL,NSNLI
          IF(ITSOK.EQ.0)
     &    WRITE(2,*) '  Sorry to say it , but memory is CORRUPTED '
          WRITE(2,*) '  Memory map : '
          WRITE(2,*) '  Identifier    Offset     Pad1 okay Pad2 okay '
          WRITE(2,*) '  ========== ============  ========= ========= '
          DO 200 IL = 1, NSNLI
            JBASE = IBASE(IL)
            IF(IL.NE.NSNLI) THEN
             JBASEN = IBASE(IL+1)
            ELSE
             JBASEN = KFREEL + 1
            END IF
            L1OK = 1
            L2OK = 1
            DO 40 IPAD = 1, NPAD
              IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1) L1OK = 0
   40       CONTINUE
            DO 50 IPAD = 1, NPAD
              IF(WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) L2OK = 0
   50       CONTINUE
C           find local off-set to WORK for print
            IBASE_L = IBASE(IL) - IMARK(1)
            IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
               WRITE(2,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'    OKAY     OKAY '
            ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
               WRITE(2,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'    OKAY       -  '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(2,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(2,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
               WRITE(2,'(4X,A,I12,4X,A)')
     &         IIDENT(IL),IBASE_L,'     -       -    '
            END IF
  200     CONTINUE
*
* Marks
*
          WRITE(2,*)
          WRITE(2,*) '======='
          WRITE(2,*) ' Marks '
          WRITE(2,*) '======='
          WRITE(2,*)
*
          WRITE(2,*) ' Identifier  Start of free memory '
          WRITE(2,*) ' ================================='
          DO JMARK = 1, NM
            WRITE(2,'(3X,A6,10X,I18)') MARKC(JMARK),IMARK(JMARK)
          END DO
*
        IF (ITSOK.EQ.0)
     &     Call Abend2(' Error observed by  memory manager ')
        END IF
*
      IF (ISTOP.NE.0)
     &   Call Abend2(' Error observed by  memory manager ')

      END
