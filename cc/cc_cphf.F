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
c /* deck cc_cphf */
*=====================================================================*
      SUBROUTINE CC_CPHF(TYPE,LABEL,ISYMS,ISTAT,EIGV,
     &                   ISYMO,FREQS,ICAU,NVEC,MAXVEC,
     &                   WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: solve CPHF equations needed in CC program
*
*    implemented types:  R1  
*
*    Written by Christof Haettig, november 1998
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
      INTEGER IPRLIN, LSYMRF, LSYMPT, LSYMST,
     &        NCONRF, NCONST, NWOPPT, NVARPT
#else
#include <implicit.h>
#endif
#include <priunit.h>
#include <dummy.h>
#include <exeinf.h>
#include <maxorb.h>
#include <infvar.h>
#include <inftap.h>
#include <iratdef.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <cclr.h>
#include <ccfro.h>
#include <inflin.h>
Cholesky
#include <ccdeco.h>
Cholesky

* local parameters:
      CHARACTER*(18) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_CPHF>  ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .true. )
      CHARACTER*8 FNGDVE, FNSOVE, FNREVE, FILCPHF
      PARAMETER(FNGDVE='CCPHFRHS',FNSOVE='CCPHFSOL')
      PARAMETER(FNREVE='CCPHFRES')
      INTEGER MXFRVAL
      PARAMETER (MXFRVAL = 100)
C     PARAMETER (MXFRVAL = 1)  ! switch off simultaneous equations

      CHARACTER TYPE*(3)

      INTEGER NVEC, MAXVEC, LWORK
      INTEGER ISYMS(MAXVEC,*), ISYMO(MAXVEC,*)
      INTEGER ISTAT(MAXVEC,*), ICAU(MAXVEC,*)

      CHARACTER*8 LABEL(MAXVEC,*)

      REAL*8  FREQS(MAXVEC,*), EIGV(MAXVEC,*)
      REAL*8  WORK(LWORK)
      REAL*8  ZERO, RDUM
      REAL*8  XNORM, DDOT, DSQRT, GDNORM

      PARAMETER (ZERO = 0.0d0)

      CHARACTER MODEL*(10), MODWR*(10)

      LOGICAL MAYBE_MORE, RELAX, CICLC, HFCLC, TRPCLC,OOTV
      LOGICAL EXCLC, FASAVE, OPTST, EX
      INTEGER IOPT, ISYM, IVEC, NSTAT, ORDER, IDUM, IPERT,LWRK0,LWRK1
      INTEGER NCMOT,NASHT,N2ASHX,LCINDX,KCMO,KUDV,KXINDX,KEND0,KEND1
      INTEGER LUGDVE, LUSOVE, MXRM, MXPHP, NCOSAV, IOPTWR, KSLVEC
      INTEGER IREAL, NFRVAL, KFRVAL, NABAOP, NABATY, IFRVAL, IDX
      INTEGER IPRABA, LUREVE, LUCPHF, MALLAI, INUM

* external functions:
      INTEGER IR1TAMP

*---------------------------------------------------------------------*
* do some checks:
*---------------------------------------------------------------------*
      CALL QENTER('CC_CPHF')

      IF (NVEC.LT.1) THEN
        CALL QEXIT('CC_CPHF')
        RETURN
      END IF

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'Entered CC_CPHF. NVEC =',NVEC
         WRITE(LUPRI,*) 'TYPE ',TYPE
         CALL FLSHFO(LUPRI)
      END IF

* check vector type:
      IF (TYPE(1:3).EQ.'R1 ') THEN
        ORDER = 1
        NSTAT = 0
        MODWR = 'SCF?      '
      ELSE
        WRITE (LUPRI,*) 'CPHF vectors ',TYPE(1:3),' not implemented.'
        CALL QUIT('required CPHF vectors not implemented.')
      END IF

* Refuse to do anything for Cholesky.
      IF (CHOINT) THEN
         NWARN = NWARN + 1
         WRITE(LUPRI,'(/A/A/A/A/)') '*** WARNING from CC_CPHF:',
     &      ' WARNING: Refusing to do CPHF for Cholesky, because',
     &      ' WARNING: ABACUS has not been modified yet!',
     &      ' WARNING: Program continues nevertheless...'
         RETURN
      ENDIF

* get some variables from SIRIUS common blocks
      CALL CC_SIRINF(NCMOT,NASHT,N2ASHX,LCINDX)

* check number of active shells from SIRIUS:
      IF (NASHT.NE.0) THEN
        WRITE (LUPRI,*) 'non-zero number of active shells:',NASHT
        CALL QUIT('Non-zero number of active shells in CC_CPHF.')
      END IF

*---------------------------------------------------------------------*
* print header for rhs vector section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '------------------------------------',
     &                               '-------------------------------',
     & '|                   OUTPUT FROM CPHF',    
     &                               ' SECTION                      |',
     & '------------------------------------',
     &                               '-------------------------------' 

* print some debug/info output
      IF (IPRINT .GT. 10 .OR. LOCDBG) THEN
        WRITE(LUPRI,*) 'CC_CPHF Workspace:',LWORK
        WRITE(LUPRI,*) '1 MODWR ',MODWR
      END IF
  
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* some initilizations:
*---------------------------------------------------------------------*

* maximum of nallai over isym (used as fixed record lengths for LUCPHF:
      MALLAI = NALLAI(1)
      DO ISYM = 2, NSYM
        MALLAI = MAX(MALLAI,NALLAI(ISYM))
      END DO

* arrays for GETGPV and ABARSP:
      KCMO   = 1
      KUDV   = KCMO   + NCMOT
      KXINDX = KUDV   + N2ASHX
      KFRVAL = KXINDX + LCINDX
      KEND0  = KFRVAL + MXFRVAL
      LWRK0  = LWORK  - KEND0

      IF (LWRK0 .LT. 0) THEN
         CALL QUIT('Insufficient memory in CC_CPHF.')
      END IF

* read MO coefficients from file:
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI(LUSIFC,IRAT*NCMOT,WORK(KCMO))
      CALL GPCLOSE(LUSIFC,'KEEP')

* flags for ABARSP:
      CICLC  = .FALSE.
      HFCLC  = .TRUE.
      TRPCLC = .FALSE.
      OOTV   = .FALSE.
      EXCLC  = .FALSE.
      MXRM   = 40
      MXPHP  = 1
 
      FASAVE = FABRHS
      NCOSAV = NCONF

      FABRHS = .TRUE.
      NCONF  = 1

      IF (DIRECT) CALL CCDFFOP

* flags for CC_WRRSP routine:
      IOPTWR = 4

* open file for right hand side and solution vectors:

      LUGDVE = -1
      LUSOVE = -1
      LUREVE = -1
      CALL GPOPEN(LUGDVE,FNGDVE,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL GPOPEN(LUSOVE,FNSOVE,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL GPOPEN(LUREVE,FNREVE,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
* open property integral file:
      IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)

* close twoelectron integral file 'AOTWOINT':
      IF (LUINTA.GT.0) CALL GPCLOSE(LUINTA,'KEEP')

* file for CPHF vectors:
      WRITE(FILCPHF,'(A5,A3)') 'CPHF_',TYPE(1:3)
      DO I = 6,8
        IF (FILCPHF(I:I).EQ.' ') FILCPHF(I:I) = '_'
      END DO

      LUCPHF = -1
      CALL WOPEN2(LUCPHF,FILCPHF,64,0)

*---------------------------------------------------------------------*
* loop over vectors, set up rhs vectors and call ABARSP to solve eqs.:
*---------------------------------------------------------------------*

      IVEC = 0
      DO WHILE(IVEC.LT.NVEC)
        IVEC = IVEC + 1

        ISYM = 1
        DO IDX = 1, NSTAT
          ISYM = MULD2H(ISYM,ISYMS(IVEC,IDX))
        END DO
        DO IDX = 1, ORDER
          ISYM = MULD2H(ISYM,ISYMO(IVEC,IDX))
        END DO

        IF (LWRK0 .LT. NALLAI(ISYM)) THEN
           CALL QUIT('Insufficient memory in CC_CPHF.')
        END IF

        IF (LOCDBG) THEN
          WRITE (LUPRI,*) 'CC_CPHF> GP vector, label, symmetry:',
     &                     LABEL(IVEC,1), ISYM
          WRITE(LUPRI,*) '2 MODWR ',MODWR
          CALL FLSHFO(LUPRI)
        END IF

C       ---------------------------
C       get right hand side vector:
C       ---------------------------
        CALL CC_GETHFGD(IVEC,TYPE,LABEL,ISYMS,ISTAT,EIGV,ISYMO,
     &                  FREQS,ICAU,NVEC,MAXVEC,IREAL,
     &                  WORK(KCMO),WORK(KUDV),WORK(KXINDX),
     &                  WORK(KFRVAL),WORK(KEND0),LWRK0)

        IF (LOCDBG) THEN
          WRITE (LUPRI,'(5X,I5,F12.8)') (I,WORK(KEND0-1+I),
     &                     I=1,NALLAI(ISYM))
          WRITE(LUPRI,*) '3 MODWR ',MODWR
        END IF


C       ------------------------------
C       write gradient vector to file:
C       ------------------------------
        REWIND LUGDVE
        CALL WRITT(LUGDVE,NALLAI(ISYM),WORK(KEND0))


C       ----------------------------------
C       check norm of the gradient vector:
C       ----------------------------------
        GDNORM=DSQRT(DDOT(NALLAI(ISYM),WORK(KEND0),1,WORK(KEND0),1))
        IF (LOCDBG) WRITE (LUPRI,*) 'GDNORM:',GDNORM
        IF (LOCDBG) WRITE (LUPRI,*) '4 MODWR ',MODWR


C       --------------------------------------------------------
C       for 'R1' check if several frequencies for same operator:
C       --------------------------------------------------------
        NFRVAL     = 1
        MAYBE_MORE = .TRUE.
        DO WHILE (MAYBE_MORE) 
           IF (TYPE(1:3).EQ.'R1 ' 
     &         .AND. IVEC.LT.NVEC .AND. NFRVAL.LT.MXFRVAL
     &         .AND. LABEL(IVEC,1).EQ.LABEL(IVEC+1,1) ) THEN
             WORK(KFRVAL+NFRVAL) = FREQS(IVEC+1,1)
             IVEC   = IVEC + 1
             NFRVAL = NFRVAL + 1
           ELSE
             MAYBE_MORE = .FALSE.
           END IF
        END DO


        IF (GDNORM.GT.THRLDPHF) THEN
C
          CALL SETSIR(WORK(KEND0),LWRK0)

          IF (LOCDBG) WRITE (LUPRI,*) 'going to call ABARSP'
          IF (LOCDBG) WRITE (LUPRI,*) '5 MODWR ',MODWR
C         -----------------------
C         solve CPHF equation(s):
C         -----------------------
          IPRABA = IPRINT
          NABAOP = 1
          NABATY = IREAL  ! flag for real/imaginary operators
          CALL ABARSP(CICLC,HFCLC,TRPCLC,OOTV,ISYM,EXCLC,
     &                WORK(KFRVAL),NFRVAL,NABATY,NABAOP,
     &                LABEL(IVEC,1),LUGDVE,LUSOVE,LUREVE,THRLEQ,
     &                MAXITE,IPRABA,MXRM,MXPHP,WORK(KEND0),LWRK0)

          IF (LOCDBG) WRITE (LUPRI,*) 'returned from ABARSP'
          IF (LOCDBG) WRITE (LUPRI,*) '6 MODWR ',MODWR

C         --------------------------------------------
C         clean up `left overs' from ABARSP:
C         --------------------------------------------
          CALL GPCLOSE(LUSIFC,'KEEP')
          IF (LUONEL .GT. 0) CALL GPCLOSE(LUONEL,'KEEP')
          IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
C
C         --------------------------------------------
C         read solution vector(s) and save on CC file:
C         --------------------------------------------
          REWIND LUSOVE

        END IF

        DO IFRVAL = 1, NFRVAL
           KSLVEC = KEND0
           KEND1  = KSLVEC + 2*MALLAI
           LWRK1  = LWORK  - KEND1
           IF (LWRK1.LT.0) THEN
              CALL QUIT('Insufficient memory in CC_CPHF.')
           END IF

           CALL DZERO(WORK(KSLVEC),2*MALLAI)
           IF (GDNORM.GT.THRLDPHF) THEN
             CALL  READT(LUSOVE,2*NALLAI(ISYM),WORK(KSLVEC))
           END IF


           ! save on CPHF vector file 
           IDX = IVEC - NFRVAL + IFRVAL
           CALL PUTWA2(LUCPHF,FILCPHF,WORK(KSLVEC),
     &                 1+2*MALLAI*(IDX-1),2*MALLAI)


           ! check if a corresponding CC vector exists
           INUM = -1  
           IF (TYPE(1:3).EQ.'R1 ') THEN
            INUM=IR1TAMP(LABEL(IDX,1),.TRUE.,WORK(KFRVAL-1+IFRVAL),ISYM)
           END IF
           IF (LOCDBG) WRITE (LUPRI,*) 'IFRVAL, TYPE, MODWR ',
     &         IFRVAL,TYPE,MODWR

           ! if yes put the CPHF also on the CC vector file
           IF (INUM.GT.0) THEN
             CALL CC_WRRSP(TYPE,INUM,ISYM,IOPTWR,MODWR,
     &                   WORK(KSLVEC),DUMMY,DUMMY,WORK(KEND1),LWRK1)
           END IF

           IF (LOCDBG) THEN
              WRITE (LUPRI,*) 
     &           'CC_CPHF> solution vector, label, freq:',
     &           LABEL(IDX,1),WORK(KFRVAL-1+IFRVAL)
              WRITE (LUPRI,'(5X,I5,F12.8)')
     &           (I,WORK(KSLVEC-1+I),I=1,2*NALLAI(ISYM))
              WRITE (LUPRI,*) 
     &           'CC_CPHF> saved CPHF solution for ',TYPE,
     &           ' equation nb. ',IDX,INUM
           END IF
           IF (LOCDBG) WRITE (LUPRI,*) '7 MODWR ',MODWR

        END DO

      END DO
*---------------------------------------------------------------------*
* that's it: close files, restore variables and return
*---------------------------------------------------------------------*
      CALL WCLOSE2(LUCPHF,FILCPHF,'KEEP')
 
      IF (LUINTM .GT. 0) CALL GPCLOSE(LUINTM,'DELETE')
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUSOVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')     

      IF (LUINTA .LE. 0) THEN
        CALL GPOPEN(LUINTA,'AOTWOINT','UNKNOWN',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
      END IF

      FABRHS = FASAVE
      NCONF  = NCOSAV

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'leave CC_CPHF'
         CALL FLSHFO(LUPRI)
      END IF
           IF (LOCDBG) WRITE (LUPRI,*) '8 MODWR ',MODWR

      CALL QEXIT('CC_CPHF')
      RETURN
      END
*=====================================================================*
*              END OF SUBROUTINE CC_CPHF                              *
*=====================================================================*
c /* deck cc_sirinf */
*=====================================================================*
      SUBROUTINE CC_SIRINF(NCMOT1,NASHT1,N2ASHX1,LCINDX1)
*---------------------------------------------------------------------*
*   Purpose: read some variables from SIRIUS common blocks
*---------------------------------------------------------------------*
#include <implicit.h>
#include <inforb.h>
#include <inftap.h>
#include <infdim.h>

      NCMOT1  = NCMOT
      NASHT1  = NASHT
      N2ASHX1 = N2ASHX
      LCINDX1 = LCINDX
      RETURN
      END 
*=====================================================================*
c /* deck cc_rdhfrsp */
*=====================================================================*
      SUBROUTINE CC_RDHFRSP(LIST,IDXLST,ISYM,XKAPPA)
*---------------------------------------------------------------------*
C
C   Purpose:  Read a CPHF response vector from file
C             for explanation of LIST, IDXLIST & MODFIL see CC_RDRSP
C
C  Christof Haettig, summer 2003         
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>    
#include <ccfro.h>    

#if defined (SYS_CRAY)
      REAL XKAPPA(*)
#else
      DOUBLE PRECISION XKAPPA(*)
#endif

      CHARACTER FILCPHF*8, LIST*3
      INTEGER MALLAI, ISYM, JSYM, LUCPHF, IDXLST

      CALL QENTER('CCRDHFRSP')

* set file name: 
      WRITE(FILCPHF,'(A5,A3)') 'CPHF_',LIST(1:3)
      DO I = 6,8
        IF (FILCPHF(I:I).EQ.' ') FILCPHF(I:I) = '_'
      END DO

* calculate record lengths:
      MALLAI = NALLAI(1)
      DO JSYM = 2, NSYM
        MALLAI = MAX(MALLAI,NALLAI(JSYM))
      END DO

* open direct access file: 
      LUCPHF = -1
      CALL WOPEN2(LUCPHF,FILCPHF,64,0)

* read vector number IDXLST from file:
      CALL GETWA2(LUCPHF,FILCPHF,XKAPPA,
     &                 1+2*MALLAI*(IDXLST-1),2*NALLAI(ISYM))

* close file and return:
      CALL WCLOSE2(LUCPHF,FILCPHF,'KEEP')

      CALL QEXIT('CCRDHFRSP')
      RETURN
      END 
*=====================================================================*
*              END OF SUBROUTINE CC_RDHFRSP                           *
*=====================================================================*
