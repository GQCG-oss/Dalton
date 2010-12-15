C*******************************************************************************
C
C This file contains the routines for interfacing the
C
C      Cholesky Decomposition Utility, Version 2.0
C
C to the programs of DALTON (i.e. Hermit, Sirius, CC, etc.)
C
C Written 2004/2005 by
C
C      Thomas Bondo Pedersen
C      Department of Theoretical Chemistry
C      Chemical Center, University of Lund
C      P.O. Box 124, S-221 00 Lund, Sweden.
C
C      Alfredo M. J. Sanchez de Meras
C      Institute of Molecular Science, Department of Physical Chemistry
C      University of Valencia
C      E-46100 Burjassot, Valencia, Spain.
C
C      Henrik Koch
C      Department of Physical Chemistry
C      Norwegian University of Science and Technology
C      N-7491 Trondheim, Norway.
C
C*******************************************************************************
C  /* Deck cho_dal_drv */
      SUBROUTINE CHO_DAL_DRV(WORK,LWORK)
C
C     Purpose: DALTON interface to Cholesky decomposition.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <cholesky.h>

C     Initialize integral program and set #shells.
C     --------------------------------------------

      CALL CHO_DAL_HERINI(WORK,LWORK,NUMSHL)
      NSHELL = NUMSHL

C     Read input.
C     -----------

      CALL CHO_DAL_INP

C     Start the Cholesky decomposition procedure.
C     -------------------------------------------

      CALL CHO_DRV(WORK,LWORK)

      END
C  /* Deck cho_dal_inp */
      SUBROUTINE CHO_DAL_INP
C
C     Purpose: read and process input or set defaults only.
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>

      CHARACTER*7 WORD

C     Input processing.
C     -----------------

      IF (LUCMD .LT. 1) THEN
         CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND(LUCMD)
      ELSE
         REWIND(LUCMD)
      END IF

  100 READ(LUCMD,'(A7)',END=1100,ERR=1200) WORD
      IF (WORD(1:7) .EQ. '*END OF') THEN
         CALL CHO_INP(.TRUE.,LUCMD,LUPRI) ! set defaults
         CALL GPCLOSE(LUCMD,'KEEP')
      ELSE IF (WORD(1:7) .EQ. '**CHOLE') THEN
         CALL CHO_INP(.FALSE.,LUCMD,LUPRI) ! process input
         CALL GPCLOSE(LUCMD,'KEEP')
      ELSE
         GO TO 100
      END IF

C     Exit.
C     -----

      RETURN

C     Error branches.
C     ---------------

 1100 CONTINUE ! end-of-file reached: set defaults and return
      CALL CHO_INP(.TRUE.,LUCMD,LUPRI)
      CALL GPCLOSE(LUCMD,'KEEP')
      RETURN

 1200 CONTINUE ! I/O error
      WRITE(LUPRI,*) SECNAM,': I/O error reading input file, unit ',
     &               LUCMD
      CALL GPCLOSE(LUCMD,'KEEP')
      CALL CHO_QUIT('I/O error while reading input',104)
      RETURN

      END
C  /* Deck cho_dal_getkey */
      SUBROUTINE CHO_DAL_GETKEY(LUNIT,OPTION,LOPTION,NOPTION,IDKEY,
     &                          LUPRI)
C
C     Purpose: get next keyword and convert it to internal IDKEY.
C
#include <implicit.h>
      CHARACTER*(*) OPTION(NOPTION) ! <-- character*(loption)

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_GETKEY')

      INTEGER CHO_TABIND

      PARAMETER (NTABLE = 35, LKEY = 7, NEOINP = 1, LEOINP = 1)
      PARAMETER (NALIAS = 5, NOBSOL = 3)
      CHARACTER*(LKEY)   TABLE(NTABLE)    ! keywords
      CHARACTER*(LKEY)   ALIAS(NALIAS,2)  ! aliases
      CHARACTER*(LKEY)   OBSOL(NOBSOL)    ! obsolete but allowed
      CHARACTER*(LEOINP) EOINP(NEOINP)    ! end-of-input characters
      CHARACTER*(LKEY)   KEY, KEYINP

      DATA TABLE /'.THRCOM','.PRINT ','.DIABUF','.THRDIA','.DAMP 1',
     &            '.DAMP 2','.SPAN  ','.MIN QU','.MAX SH','.SCREEN',
     &            '.NOSCRE','.QUALIF','.THRNEG','.THRWAR','.TOONEG',
     &            '.CHKCNF','.CHK AL','.RSTDIA','.RSTCHO','.RSTMOD',
     &            '.MAX QU','.MAXRED','.MAXVEC','.CHK SP','.CHK MI',
     &            '.USEABS','.NOABSO','.TRCNEG','.IO VEC','.ONESTE',
     &            '.TWOSTE','.NAIVE ','.BUFFER','.DENDEC','.XXXXXX'/
      DATA ALIAS /'.THINDI','.THSUDI','.NOSCDI','.SPANDI','.LENBUF',
     &            '.DAMP 1','.DAMP 2','.NOSCRE','.SPAN  ','.BUFFER'/
      DATA OBSOL /'.NOCOMP','.REDUCE','.THRSTR'/
      DATA EOINP /'*'/

C     Check that we are in sync with caller.
C     --------------------------------------

      IF (NOPTION .NE. NTABLE) THEN
         WRITE(LUPRI,*) SECNAM,': NOPTION = ',NOPTION,
     &                          ' NTABLE = ',NTABLE
         IDKEY = -5
         GO TO 2000
      END IF

C     Read keyword. Disregard comments.
C     ---------------------------------

      NUMOBS = 0  ! count #obsolete keywords encountered

    1 READ(LUNIT,'(A7)',END=1100,ERR=1200) KEY
      DO WHILE (KEY(1:1).EQ.'!' .OR. KEY(1:1).EQ.'#')
         READ(LUNIT,'(A7)',END=1100,ERR=1200) KEY
      END DO
      KEYINP = KEY

C     Check for obsolete keyword.
C     ---------------------------

      IF (NOBSOL .GT. 0) THEN
         IOBSOL = CHO_TABIND(OBSOL,LKEY,NOBSOL,' ',0,0,KEY)
         IF (IOBSOL.GT.0 .AND. IOBSOL.LE.NOBSOL) THEN
            NUMOBS = NUMOBS + 1
            IF (NUMOBS .EQ. 1) THEN
               WRITE(LUPRI,'(/,A,A7,A,A)')
     &         '*** NOTICE: Keyword "',KEY,'" is obsolete in this ',
     &         'version of the Cholesky decomposition. Disregarded.'
            ELSE
               WRITE(LUPRI,'(12X,A,A7,A,A)')
     &         'Keyword "',KEY,'" is obsolete in this ',
     &         'version of the Cholesky decomposition. Disregarded.'
            END IF
            GO TO 1
         END IF
      END IF

C     Check for aliases.
C     ------------------

      IALIAS = 0
      IF (NALIAS .GT. 0) THEN
         IALIAS = CHO_TABIND(ALIAS(1,1),LKEY,NALIAS,' ',0,0,KEY)
         IF (IALIAS.GT.0 .AND. IALIAS.LE.NALIAS) THEN
            KEY = ALIAS(IALIAS,2)
         ELSE
            IALIAS = 0
         END IF
      END IF

C     Table look-up.
C     --------------

      IDKEY = CHO_TABIND(TABLE,LKEY,NTABLE,EOINP,LEOINP,NEOINP,KEY)
      IF (IDKEY .EQ. -1) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,': keyword not recognized:'
         WRITE(LUPRI,*) 'Keyword from input: ',KEYINP
         IF (IALIAS .NE. 0) THEN
            WRITE(LUPRI,*) 'Internal key used : ',KEY,' (aliasing)'
         END IF
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Available keywords (and short explanations):'
         DO IOPTION = 1,NOPTION
            IF (TABLE(IOPTION) .NE. '.XXXXXX') THEN
               WRITE(LUPRI,*) TABLE(IOPTION),': ',OPTION(IOPTION)
            END IF
         END DO
         IF (NALIAS .GT. 0) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Aliases:'
            DO IALIAS = 1,NALIAS
               WRITE(LUPRI,*) ALIAS(IALIAS,1),'<--->',ALIAS(IALIAS,2)
            END DO
         END IF
         IF (NOBSOL .GT. 0) THEN
            WRITE(LUPRI,*)
            WRITE(LUPRI,*) 'Obsolete keywords accepted but disregarded:'
            DO IOBSOL = 1,NOBSOL
               WRITE(LUPRI,*) OBSOL(IOBSOL)
            END DO
         END IF
         WRITE(LUPRI,*)
      END IF

C     Normal exit point.
C     ------------------

 2000 CALL CHO_FLUSH(LUPRI)
      RETURN
c2000 RETURN

C     Error branches.
C     ---------------

 1100 CONTINUE ! end-of-file reached
      WRITE(LUPRI,*) SECNAM,': EOF reached while reading unit ',LUNIT
      IDKEY = -5
      RETURN

 1200 CONTINUE ! I/O error
      WRITE(LUPRI,*) SECNAM,': I/O error while reading unit ',LUNIT
      IDKEY = -5
      RETURN

      END
C  /* Deck cho_dal_herini */
      SUBROUTINE CHO_DAL_HERINI(WORK,LWORK,NUMSHL)
C
C     Purpose: initialize integral program and get #shells.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <blocks.h>

      IPRINT = 0
      CALL HERDI1(WORK,LWORK,IPRINT)
      NUMSHL = MAXSHL

      END
C  /* Deck cho_dal_inimul */
      BLOCK DATA CHO_DAL_INIMUL
C
C     Purpose: initialize group multiplication table in cholesky.h
C
#include <cholesky.h>

      DATA MULD2H/1,2,3,4,5,6,7,8,
     &            2,1,4,3,6,5,8,7,
     &            3,4,1,2,7,8,5,6,
     &            4,3,2,1,8,7,6,5,
     &            5,6,7,8,1,2,3,4,
     &            6,5,8,7,2,1,4,3,
     &            7,8,5,6,3,4,1,2,
     &            8,7,6,5,4,3,2,1/

      END
C  /* Deck cho_dal_init */
      SUBROUTINE CHO_DAL_INIT(WORK,LWORK)
C
C     Purpose: initialization of Cholesky decomposition.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
      EXTERNAL CHO_DAL_INIMUL
#include <cholesky.h>
#include <choorb.h>
#include <chomeman.h>

      INTEGER INDEXA(CHO_MAXBAS), INDEXB(CHO_MAXBAS)

      LOGICAL SETUP
      PARAMETER (SETUP = .TRUE.)

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_INIT')

      ITRI(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

C     Initialize memory manager.
C     (Is this really needed???)
C     --------------------------

      LWORK_TOT = LWORK

C     Check that the number of shells is within limits.
C     -------------------------------------------------

      IF ((NSHELL.LT.1) .OR. (NSHELL.GT.CHO_MAXSHL)) THEN
         WRITE(LUPRI,*) 'NSHELL out of bounds: ',NSHELL
         WRITE(LUPRI,*) 'Max. allowed # shells in Cholesky',
     &                  ' decomposition is CHO_MAXSHL = ',CHO_MAXSHL
         CALL CHO_QUIT('NSHELL out of bounds in '//SECNAM,102)
      END IF

C     Get symmetry and basis set information.
C     NSYM      : #irreps
C     NBAS(ISYM): # basis functions (SOs) in symmetry ISYM
C     IBAS(ISYM): offset to basis functions in symmetry ISYM
C     NBAST     : total number of basis functions
C     Note: the information resides in cholesky.h as well as
C           in ccorb.h after this routine.
C     ------------------------------------------------------

      CALL CHO_DAL_GETBAS(IBAS,NBAS,NBAST,NSYM,LUPRI)

C     Check total number of basis functions.
C     --------------------------------------

      IF ((NBAST.LT.1) .OR. (NBAST.GT.CHO_MAXBAS)) THEN
         WRITE(LUPRI,*) 'NBAST out of bounds: ',NBAST
         WRITE(LUPRI,*) 'Max. allowed #orbitals in Cholesky',
     &                  ' decomposition is CHO_MAXBAS = ',CHO_MAXBAS
         CALL CHO_QUIT('NBAST out of bounds in '//SECNAM,102)
      END IF
      CALL CHO_DAL_CHKBAS(NBAST,CHO_MAXBAS,LUPRI)

C     Make sure chov2io.h is initialized (only used after
C     decomposition).
C     ---------------------------------------------------

      CALL CHO_DAL_CHOV2IO()

C     Initialize commons used in "old" decomposition.
C     I.e., ccsdsym.h etc. Set also flag CHOINT in
C     ccdeco.h (used in integral program) and the 
C     buffer length used by the integral program.
C     -----------------------------------------------

      CALL CHO_DAL_INDEX(LERIBUF)

C     Get shell info from integral program.
C     ISOSHL(I)        : shell to which SO I belongs
C     NBASSH(ISYM,ISHL): #SO of sym. ISYM in shell ISHL
C     -------------------------------------------------

      ISHLB = 1
      LPRNT = 0
      DO ISHLA = 1,NSHELL
         CALL NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHLA,ISHLB,
     &               NUMA,NUMB,LPRNT,SETUP,1)
         DO IA = 1,NUMA
            ISOSHL(INDEXA(IA)) = ISHLA
         END DO
      END DO

      CALL CHO_IZERO(NBASSH,8*NSHELL)
      DO ISYM = 1,NSYM
         DO IA = 1,NBAS(ISYM)
            I = IBAS(ISYM) + IA
            ISHL = ISOSHL(I)
            NBASSH(ISYM,ISHL) = NBASSH(ISYM,ISHL) + 1
         END DO
      END DO

C     NBSTSH(ISHL): total dimension of shell ISHL
C     MXORSH      : max. shell dimension
C     -------------------------------------------

      MXORSH = 0
      DO ISHL = 1,NSHELL
         NBSTSH(ISHL)   = NBASSH(1,ISHL)
         IBASSH(1,ISHL) = 0
         DO ISYM = 2,NSYM
            IBASSH(ISYM,ISHL) = NBSTSH(ISHL)
            NBSTSH(ISHL) = NBSTSH(ISHL) + NBASSH(ISYM,ISHL)
         END DO
         MXORSH = MAX(MXORSH,NBSTSH(ISHL))
      END DO

C     MX2SH: max. dimension of shell pair.
C     ------------------------------------

      MX2SH = 0
      DO ISHL = 1,NSHELL
         DO JSHL = 1,ISHL
            IJSHL = ITRI(ISHL,JSHL)
            IF (ISHL .EQ. JSHL) THEN
               NUMIJ = NBSTSH(ISHL)*(NBSTSH(ISHL) + 1)/2
            ELSE
               NUMIJ = NBSTSH(ISHL)*NBSTSH(JSHL)
            END IF
            MX2SH = MAX(MX2SH,NUMIJ)
         END DO
      END DO

C     ISHLSO(I): address in shell to which SO I belongs.
C     --------------------------------------------------

      DO ISHL = 1,NSHELL
         ICOUNT = 0
         I      = 0
         DO WHILE ((I.LT.NBAST) .AND. (ICOUNT.LT.NBSTSH(ISHL)))
            I = I + 1
            IF (ISOSHL(I) .EQ. ISHL) THEN
               ICOUNT = ICOUNT + 1
               ISHLSO(I) = ICOUNT
            END IF
         END DO
      END DO

C     Set static units for the 64 bit crayio routines.
C     Define file names.
C     ------------------------------------------------

      DO ISYM = 1,NSYM
         LUNCHO(ISYM) = 30 + ISYM
      END DO
      LURSH = -1  ! INDRSH file handled by Fortran I/O
      LURED = -1  ! reduced set index file handled by Fortran I/O
      LURST = -1  ! restart file handled by Fortran I/O
      LUMAP = -1  ! integral map file handled by Fortran I/O

      FNRSH = 'CHRSH'
      FNRED = 'CHRED'
      DO ISYM = 1,NSYM
         WRITE(FNVEC(ISYM),'(A5,I1)') 'CHVEC',ISYM
      END DO
      FRST = 'CHORST'
      FMAP = 'CHOMAP'

#if defined (CHODAL_DEBUG)
C     Set mapping from shell index to global AO index.
C     ------------------------------------------------

      CALL IZERO(ISHL2F,CHO_MAXBAS*CHO_MAXSHL)
      DO I = 1,NBAST
         ISHL = ISOSHL(I)
         IA   = ISHLSO(I)
         IF (ISHL2F(IA,ISHL).NE.0) THEN
            WRITE(LUPRI,*) SECNAM,': ISHL2F(',IA,',',ISHL,') already ',
     &                     'set to: ',ISHL2F(IA,ISHL)
            CALL CHO_QUIT('Setup error for ISHL2F in '//SECNAM,103)
         ELSE
            ISHL2F(IA,ISHL) = I
         END IF
      END DO
#endif

      END
C  /* Deck cho_dal_index */
      SUBROUTINE CHO_DAL_INDEX(LOCBUF)
C
C     Purpose: initializations of Dalton common blocks not used
C              in the decomposition itself but used by integral
C              program (or interface to it). Also set default values
C              for the Cholesky SCF Fock matrix calculation.
C
#include <implicit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <choscf.h>
#include <mxcent.h>
#include <eribuf.h>

      CHOINT = .TRUE.
      LBUF   = LOCBUF
      CALL CHO_INDEX

      END
C  /* Deck cho_dal_chkbas */
      SUBROUTINE CHO_DAL_CHKBAS(LNBAST,LCHOMXB,LUPRI)
C
C     Purpose: check number of basis functions.
C
#include <implicit.h>
#include <maxorb.h>

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_CHKBAS')

      IF (LNBAST.GT.MXCORB .OR. LNBAST.GT.LCHOMXB) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'NBAST out of bounds: ',LNBAST
         WRITE(LUPRI,*) 'Max. allowed #orbitals in Cholesky',
     &                  ' decomposition is CHO_MAXBAS = ',LCHOMXB
         WRITE(LUPRI,*) 'Max. allowed #orbitals in Dalton',
     &                  ' part is MXCORB = ',MXCORB
         IF (LCHOMXB .NE. MXCORB) THEN
            WRITE(LUPRI,*) 'Suggestion: bring these max. values in ',
     &                     'sync and recompile.'
         END IF
         CALL CHO_QUIT('NBAST out of bounds in '//SECNAM,102)
      END IF

      END
C  /* Deck cho_dal_getbas */
      SUBROUTINE CHO_DAL_GETBAS(IBAS,NBAS,NBAST,NSYM,LUPRI)
C
C     Purpose: set up basis set information.
C
#include <implicit.h>
      INTEGER IBAS(8), NBAS(8)
#include <inftap.h>
#include <dummy.h>

      DIMENSION TITMOL(12,2)

C     Read number of basis functions.
C     -------------------------------

      IF (LUONEL .LT. 1) THEN
         CALL GPOPEN(LUONEL,'AOONEINT','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
      END IF
      REWIND(LUONEL)
      READ (LUONEL) ((TITMOL(I,J),I=1,12),J=1,2),
     &               NSYM,(NBAS(ISYM),ISYM=1,NSYM),POTNUC
      CALL GPCLOSE(LUONEL,'KEEP')

C     Set up offset array IBAS.
C     -------------------------

      IBAS(1) = 0
      NBAST   = NBAS(1)
      DO ISYM = 2,NSYM
         IBAS(ISYM) = NBAST
         NBAST = NBAST + NBAS(ISYM)
      END DO

C     Feed info to ccorb.h.
C     ---------------------

      CALL CHO_DAL_CCORB(IBAS,NBAS,NBAST,NSYM)

      END
C  /* Deck cho_dal_ccorb */
      SUBROUTINE CHO_DAL_CCORB(LIBAS,LNBAS,LNBAST,LNSYM)
C
C     Purpose: populate ccorb.h.
C
#include <implicit.h>
      INTEGER LIBAS(8), LNBAS(8)
#include <ccorb.h>

      NSYM  = LNSYM
      NBAST = LNBAST

      DO ISYM = 1,NSYM
         IBAS(ISYM) = LIBAS(ISYM)
         NBAS(ISYM) = LNBAS(ISYM)
      END DO

      END
C  /* Deck cho_dal_mem */
      SUBROUTINE CHO_DAL_MEM(TEXT,JOB,TYPE,KPOINT,LENGTH,KFREE,LFREE)
C
C     Purpose: memory management.
C
C     NOTE: this will only work for sequential allocations!!!
C
#include <implicit.h>
      CHARACTER*(*) TEXT,JOB,TYPE
#include <cholesky.h>
#include <chomeman.h>

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_MEM')

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER CHO_DAL_MEM_INT

      IF (LOCDBG) THEN
         LENTXT = MIN(LEN(TEXT),72)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,': input parameters:'
         WRITE(LUPRI,*) 'TEXT  : ',TEXT(1:LENTXT)
         WRITE(LUPRI,*) 'JOB   : ',JOB(1:4)
         WRITE(LUPRI,*) 'TYPE  : ',TYPE(1:4)
         WRITE(LUPRI,*) 'KPOINT: ',KPOINT
         WRITE(LUPRI,*) 'LENGTH: ',LENGTH
         WRITE(LUPRI,*) 'KFREE : ',KFREE
         WRITE(LUPRI,*) 'LFREE : ',LFREE
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,': global available memory: ',
     &                  LWORK_TOT
         CALL CHO_FLUSH(LUPRI)
      END IF

      IF (JOB(1:4) .EQ. 'MAX ') THEN
         IF (LFREE .GT. 0) THEN
            IF (TYPE(1:4) .EQ. 'INTE') THEN
               WRITE(LUPRI,*) SECNAM,': integer "MAX" ',
     &                        'allocation has not (yet) been ',
     &                        'implemented...sorry!!'
               CALL CHO_DAL_MEM_ERR(' ',TEXT,JOB,TYPE,KPOINT,LENGTH,
     &                              KFREE,LFREE,LUPRI)
            ELSE IF (TYPE .EQ. 'REAL') THEN
               KPOINT = KFREE
               LENGTH = LFREE
               KFREE  = KPOINT + LENGTH
               LFREE  = 0
            ELSE
               CALL CHO_DAL_MEM_ERR(SECNAM//': "TYPE" not recognized',
     &                              TEXT,JOB,TYPE,KPOINT,LENGTH,
     &                              KFREE,LFREE,LUPRI)
            END IF
         ELSE
            CALL CHO_DAL_MEM_ERR(SECNAM//': insufficient memory for '
     &                           //'allocation',TEXT,JOB,TYPE,
     &                           KPOINT,LENGTH,KFREE,LFREE,LUPRI)
         END IF
      ELSE IF (JOB(1:4) .EQ. 'ALLO') THEN
         IF (LENGTH .LT. 0) THEN
            CALL CHO_DAL_MEM_ERR(SECNAM//': negative LENGTH requested',
     &                           TEXT,JOB,TYPE,KPOINT,LENGTH,
     &                           KFREE,LFREE,LUPRI)
         END IF
         IF (LFREE .GE. 0) THEN
            IF (TYPE(1:4) .EQ. 'INTE') THEN
               ILENGT = CHO_DAL_MEM_INT(LENGTH)
               KPOINT = KFREE
               KFREE  = KPOINT + ILENGT
               LFREE  = LFREE  - ILENGT
            ELSE IF (TYPE(1:4) .EQ. 'REAL') THEN
               KPOINT = KFREE
               KFREE  = KPOINT + LENGTH
               LFREE  = LFREE  - LENGTH
            ELSE
               CALL CHO_DAL_MEM_ERR(SECNAM//': "TYPE" not recognized',
     &                              TEXT,JOB,TYPE,KPOINT,LENGTH,
     &                              KFREE,LFREE,LUPRI)
            END IF
         ELSE
            CALL CHO_DAL_MEM_ERR(SECNAM//': insufficient memory for '
     &                           //'allocation',TEXT,JOB,TYPE,
     &                           KPOINT,LENGTH,KFREE,LFREE,LUPRI)
         END IF
         IF (LFREE .LT. 0) THEN
            CALL CHO_DAL_MEM_ERR(SECNAM//': insufficient memory for '
     &                           //'allocation',TEXT,JOB,TYPE,
     &                           KPOINT,LENGTH,KFREE,LFREE,LUPRI)
         END IF
      ELSE IF (JOB(1:4) .EQ. 'FLUS') THEN
         IF (KPOINT.LT.1 .OR. KPOINT.GT.KFREE) THEN
            CALL CHO_DAL_MEM_ERR(SECNAM//': KPOINT error',
     &                           TEXT,JOB,TYPE,
     &                           KPOINT,LENGTH,KFREE,LFREE,LUPRI)
         ELSE
            LFREE = LFREE + (KFREE - KPOINT)
            KFREE = KPOINT
         END IF
      ELSE IF (JOB(1:4) .EQ. 'FREE') THEN
         IF (KPOINT.LT.1 .OR. KPOINT.GT.KFREE) THEN
            CALL CHO_DAL_MEM_ERR(SECNAM//': KPOINT error',
     &                           TEXT,JOB,TYPE,
     &                           KPOINT,LENGTH,KFREE,LFREE,LUPRI)
         ELSE IF (LENGTH.LT.0) THEN
            CALL CHO_DAL_MEM_ERR(SECNAM//': LENGTH error',
     &                           TEXT,JOB,TYPE,
     &                           KPOINT,LENGTH,KFREE,LFREE,LUPRI)
         ELSE
            KFREE = KPOINT
            LFREE = LFREE + LENGTH
         END IF
      ELSE
         WRITE(LUPRI,*) SECNAM,': request not recognized: ',JOB(1:4)
         CALL CHO_QUIT('JOB error in '//SECNAM,104)
      END IF

      IF (LOCDBG) THEN
         LENTXT = MIN(LEN(TEXT),72)
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SECNAM,': output parameters:'
         WRITE(LUPRI,*) 'TEXT  : ',TEXT(1:LENTXT)
         WRITE(LUPRI,*) 'JOB   : ',JOB(1:4)
         WRITE(LUPRI,*) 'TYPE  : ',TYPE(1:4)
         WRITE(LUPRI,*) 'KPOINT: ',KPOINT
         WRITE(LUPRI,*) 'LENGTH: ',LENGTH
         WRITE(LUPRI,*) 'KFREE : ',KFREE
         WRITE(LUPRI,*) 'LFREE : ',LFREE
         CALL CHO_FLUSH(LUPRI)
      END IF

      END
C  /* Deck cho_dal_mem_err */
      SUBROUTINE CHO_DAL_MEM_ERR(STRING,TEXT,JOB,TYPE,KPOINT,LENGTH,
     &                           KFREE,LFREE,LUPRI)
C
C     Purpose: error handling in memory manager.
C
      IMPLICIT NONE
      CHARACTER*(*) STRING, TEXT, JOB, TYPE
      INTEGER KPOINT, LENGTH, KFREE, LFREE, LUPRI
#include <chomeman.h>

      INTEGER LENTXT

      LENTXT = MIN(LEN(STRING),79)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) STRING(1:LENTXT)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) 'Note: the total available memory is: ',LWORK_TOT 

      LENTXT = MIN(LEN(TEXT),73)
      WRITE(LUPRI,*)
      WRITE(LUPRI,*) 'The control arguments of the last call to ',
     &               'CHO_MEM are:'
      WRITE(LUPRI,*) 'TEXT: ',TEXT(1:LENTXT)
      WRITE(LUPRI,*) 'JOB : ',JOB(1:4)
      WRITE(LUPRI,*) 'TYPE: ',TYPE(1:4)
      WRITE(LUPRI,*) 'The pointer arguments of the last call',
     &               ' to CHO_MEM are (or have become):'
      WRITE(LUPRI,*) 'KPOINT: ',KPOINT
      WRITE(LUPRI,*) 'LENGTH: ',LENGTH
      WRITE(LUPRI,*) 'KFREE : ',KFREE
      WRITE(LUPRI,*) 'LFREE : ',LFREE

      CALL CHO_QUIT('Memory error',101)

      END
C  /* Deck cho_dal_mem_int */
      INTEGER FUNCTION CHO_DAL_MEM_INT(LENGTH)
C
C     Purpose: return LENGTH in double precision words
C              (for allocation of integers).
C
#include <implicit.h>
#include <iratdef.h>

      IF (LENGTH .LT. 1) THEN
         CHO_DAL_MEM_INT = 0
      ELSE
         CHO_DAL_MEM_INT = (LENGTH - 1)/IRAT + 1
      END IF

      END
C  /* cho_dal_openvr */
      SUBROUTINE CHO_DAL_OPENVR(IOPT)
C
C     Purpose: open (IOPT=1) or close (IOPT=2) files for vector
C              and index array storage.
C
#include <implicit.h>
#include <cholesky.h>
#include <choorb.h>

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_OPENVR')

      LOGICAL OLDDX

      IF (IOPT .EQ. 1) THEN
         DO ISYM = 1,NSYM
            CALL WOPEN2U(LUNCHO(ISYM),FNVEC(ISYM),64,0)
         END DO
         LUMAP = -1
         IDUM  = -1
         CALL GPOPEN(LUMAP,FMAP,'UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUM,.FALSE.)
         LURED = -1
         LRLEN = NSYM*(NSHELL*(NSHELL+1)/2) + NBAST*(NBAST+1)
         CALL GPOPEN(LURED,FNRED,'UNKNOWN','DIRECT',
     &               'UNFORMATTED',LRLEN,OLDDX)
         LURSH = -1
         IDUM  = -1
         CALL GPOPEN(LURSH,FNRSH,'UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUM,.FALSE.)
         LURST = -1
         IDUM  = -1
         CALL GPOPEN(LURST,FRST,'UNKNOWN','SEQUENTIAL',
     &               'UNFORMATTED',IDUM,.FALSE.)
      ELSE IF (IOPT .EQ. 2) THEN
C        DO ISYM = 1,NSYM
C           CALL WCLOSE2U(LUNCHO(ISYM),FNVEC(ISYM),'KEEP')
C        END DO
C        CALL GPCLOSE(LUMAP,'KEEP')
C        CALL GPCLOSE(LURSH,'KEEP')
C        CALL GPCLOSE(LURED,'KEEP')
         CALL GPCLOSE(LURST,'KEEP')
C        LUMAP = -1
C        LURSH = -1
C        LURED = -1
         LURST = -1
      ELSE
         WRITE(LUPRI,*) SECNAM,': IOPT out of bounds: ',IOPT
         CALL CHO_QUIT('Error in '//SECNAM,104)
      END IF

      END
C  /* Deck cho_dal_gabn */
      SUBROUTINE CHO_DAL_GABN(ISHLA,ISHLB,DIAG,WORK,LWORK)
C
C     Written by Henrik Koch 25-may-2004
C
C     Calculate diagonal elements (AB|AB).
C
#include <implicit.h>
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include <maxash.h>
#include <maxorb.h>
#include <aovec.h>
#include <priunit.h>
#include <iratdef.h>
      LOGICAL   SETUP
      DIMENSION DIAG(*), WORK(LWORK)
      DIMENSION INDEXA(MXCORB), INDEXB(MXCORB)
      DIMENSION JNDEXA(MXCORB), JNDEXB(MXCORB)
#include <ccorb.h>
#include <ccisao.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <cbieri.h>
#include <distcl.h>
#include <mxcent.h>
#include <eribuf.h>
#include <ccfro.h>
#include <ccfop.h>
#include <ccsections.h>
#include <symsq.h>
#include <ccdeco.h>
C
      ITRI(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      DIACAL = .TRUE.
C
      SETUP = .FALSE.
      CALL NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHLA,ISHLB,
     *            NUMDIA,NUMDIB,IPRINT,SETUP,2)
C
C-tbp:
C These lines have been replaced to make storage in agreement with the
C expectations of version 2 of the decomposition.
CORIG DO I = 1,NUMDIA
CORIG    JNDEXA(INDEXA(I)) = I
CORIG ENDDO
CORIG DO I = 1,NUMDIB
CORIG    JNDEXB(INDEXB(I)) = I
CORIG ENDDO
      CALL CHO_DAL_COPY_ISHLSO(JNDEXA,INDEXA,NUMDIA)
      CALL CHO_DAL_COPY_ISHLSO(JNDEXB,INDEXB,NUMDIB)
C
      IF (ISHLA .EQ. ISHLB) THEN
C
         KDIAG = 1
         KIREC = KDIAG + NUMDIA*NUMDIB
         KEND1 = KIREC + (NBUFX(0) - 1)/IRAT + 1
         LWRK1 = LWORK - KEND1 + 1
         IF (LWRK1 .LT. 0) THEN
            CALL CHO_QUIT('Not enough space in CHO_DAL_GABN [1]',101)
         END IF
C
         NEWDIS = .TRUE.
         CALL CC_RDAOX(WORK(KDIAG),NUMDIA,NUMDIB,JNDEXA,JNDEXB,
     *                 ISHLA,ISHLB,WORK(KEND1),LWRK1,WORK(KIREC))
C
         DO IB = 1,NUMDIB
            KOFF = KDIAG + NUMDIA*(IB - 1) - 1
            DO IA = IB,NUMDIB
               IAB = ITRI(IA,IB)
               KAB = KOFF + IA
               DIAG(IAB) = WORK(KAB)
            END DO
         END DO
C
      ELSE
C
         KIREC = 1
         KEND1 = KIREC + (NBUFX(0) - 1)/IRAT + 1
         LWRK1 = LWORK - KEND1 + 1
         IF (LWRK1 .LT. 0) THEN
            CALL CHO_QUIT('Not enough space in CHO_DAL_GABN [2]',101)
         END IF
C
         NEWDIS = .TRUE.
         CALL CC_RDAOX(DIAG,NUMDIA,NUMDIB,JNDEXA,JNDEXB,
     *                 ISHLA,ISHLB,WORK(KEND1),LWRK1,WORK(KIREC))
C
      END IF
C
      DIACAL = .FALSE.
C
      RETURN
      END
C  /* Deck cho_dal_copy_ishlso */
      SUBROUTINE CHO_DAL_COPY_ISHLSO(JNDX,INDX,NUM)
C
C     Purpose: copy section of ISHLSO (in cholesky.h) to JNDX.
C
C     Input:
C            NUM  : dimension (of shell)
C            INDX : mapping from shell to global AO list
C
C
#include <implicit.h>
      DIMENSION JNDX(*), INDX(NUM)
#include <cholesky.h>

      DO I = 1,NUM
         JNDX(INDX(I)) = ISHLSO(INDX(I))
      END DO

      END
C  /* Deck cho_dal_calcint */
      SUBROUTINE CHO_DAL_CALCINT(IIBSTRSH,NNBSTRSH,INDRED,INDRSH,
     &                           MSYM,MMSHL,MMBSTRT,ISHLA,ISHLB,
     &                           WORK,LWORK)
C
C     Purpose: calculate qualified integral columns from
C              shell pair distribution (**|ISHLA ISHLB).
C
#include <implicit.h>
      INTEGER   INDRED(MMBSTRT,3), INDRSH(MMBSTRT)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      DIMENSION WORK(LWORK)
#include <cholesky.h>
#include <iratdef.h>
#if defined (CHODAL_DEBUG)
#include <choorb.h>
      INTEGER CHO_ISAO
      INTEGER LIAODPK(8,8)
#endif

      CHARACTER*15 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_CALCINT')

      INTEGER IQUVEA(CHO_MAXQUAL,8), IQUVEB(CHO_MAXQUAL,8)
      INTEGER LNNBST(8)

#if defined (CHODAL_DEBUG)
      ITRI(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
#endif

C     Initialization.
C     ---------------

      CALL CHO_DAL_NNBST(LNNBST,NSYM)

C     Calculate (**|ISHLA ISHLB) and store on disk.
C     Set up index arrays IQUVEA and IQUVEB needed
C     for reading the integral buffers from Hermit.
C     ---------------------------------------------

      CALL CHO_TIMER(C1,W1)
      CALL CHO_DAL_CALCINT1(ISHLA,ISHLB,WORK,LWORK,IPRINT,
     &                      INDRED,IQUVEA,IQUVEB,
     &                      MMBSTRT,CHO_MAXQUAL)
      CALL CHO_TIMER(C2,W2)
      TINTEG(1,1) = TINTEG(1,1) + C2 - C1
      TINTEG(2,1) = TINTEG(2,1) + W2 - W1

C     Get qualified columns.
C     ----------------------

      DO ISYM = 1,NSYM
         NAB = NQUAL(ISYM) - IOFFQ(ISYM)
         IF (NAB .GT. 0) THEN

C           Read in full storage mode.
C           --------------------------

            CALL CHO_TIMER(C1,W1)

            KINT1 = 1
            KEND1 = KINT1 + LNNBST(ISYM)*NAB
            LWRK1 = LWORK - KEND1 + 1
            IF (LWRK1 .LT. 0) THEN
               CALL CHO_QUIT('Insufficient memory [1] in '//SECNAM,101)
            END IF

            CALL CHO_DAL_RDAO(WORK(KINT1),IQUVEA(1,ISYM),IQUVEB(1,ISYM),
     &                        NAB,ISYM,WORK(KEND1),LWRK1)

            CALL CHO_TIMER(C2,W2)
            TINTEG(1,1) = TINTEG(1,1) + C2 - C1
            TINTEG(2,1) = TINTEG(2,1) + W2 - W1

C           Get columns in current reduced set.
C           -----------------------------------

            KINTQ = KEND1
            KR2FL = KINTQ + NNBSTR(ISYM,2)*NAB
            KEND2 = KR2FL + (NNBSTR(ISYM,2) - 1)/IRAT + 1
            LWRK2 = LWORK - KEND2 + 1
            IF (LWRK2 .LT. 0) THEN
               CALL CHO_QUIT('Insufficient memory [2] in '//SECNAM,101)
            END IF

            CALL CHO_DAL_RS2FULL(WORK(KR2FL),ISYM,
     &                           IIBSTRSH,NNBSTRSH,INDRED,
     &                           MSYM,MMSHL,MMBSTRT)
            CALL CHO_DAL_GETQCOL(WORK(KINT1),WORK(KINTQ),WORK(KR2FL),
     &                           LNNBST(ISYM),NNBSTR(ISYM,2),NAB)

#if defined (CHODAL_DEBUG)
C           Check extraction.
C           -----------------

            CALL CHO_DAL_IAODPK(LIAODPK,NSYM)

            DO I = 1,NAB
               DO J = 1,NNBSTR(ISYM,2)
                  K   = INDRED(IIBSTR(ISYM,2)+J,2) ! addr. in rs1
                  IAG = IRS1TOF(1,K)
                  IBG = IRS1TOF(2,K)
                  ISYMA = CHO_ISAO(IAG)
                  ISYMB = CHO_ISAO(IBG)
                  ISYMAB = MULD2H(ISYMA,ISYMB)
                  IF (ISYMAB .NE. ISYM) THEN
                     CALL CHO_QUIT('Symmetry error [1] in '//SECNAM,103)
                  END IF
                  IF (ISYMA .LT. ISYMB) THEN
                     CALL CHO_QUIT('Symmetry error [2] in '//SECNAM,103)
                  END IF
                  IA  = IAG - IBAS(ISYMA)
                  IB  = IBG - IBAS(ISYMB)
                  IF (ISYMA .EQ. ISYMB) THEN
                     KOFF1 = KINT1 - 1 + LNNBST(ISYM)*(I-1)
     &                     + LIAODPK(ISYMA,ISYMB)
     &                     + ITRI(IA,IB)
                  ELSE
                     KOFF1 = KINT1 - 1 + LNNBST(ISYM)*(I-1)
     &                     + LIAODPK(ISYMA,ISYMB)
     &                     + NBAS(ISYMA)*(IB-1) + IA
                  END IF
                  KOFF2 = KINTQ - 1 + NNBSTR(ISYM,2)*(I-1) + J
                  DIFF = WORK(KOFF1) - WORK(KOFF2)
                  IF (ABS(DIFF) .GT. 1.0D-20) THEN
                     WRITE(LUPRI,*)
     &               SECNAM,': difference full - reduced: ',
     &               DIFF,' is too large!'
                     WRITE(LUPRI,*)
     &               'Value in full storage   : ',WORK(KOFF1)
                     WRITE(LUPRI,*)
     &               'Value in reduced storage: ',WORK(KOFF2)
                     WRITE(LUPRI,*) 'Qualified number ',I,
     &                              ' of symmetry ',ISYM
                     L = IQUAB(I,ISYM)
                     M = INDRED(L,2)
                     WRITE(LUPRI,*) 'ID of qualified in rs1,rs2  : ',M,L
                     WRITE(LUPRI,*) 'Global AO indices of qual.  : ',
     &                              IRS1TOF(1,M),IRS1TOF(2,M)
                     WRITE(LUPRI,*) 'ID of element in rs1,rs2    : ',
     &                              K,IIBSTR(ISYM,2)+J
                     WRITE(LUPRI,*) 'Global AO indices of element: ',
     &                              IAG,IBG
                     CALL CHO_QUIT('Integral extraction error in '
     &                             //SECNAM,103)
                  END IF
               END DO
            END DO
#endif

C           Write columns to disk.
C           ----------------------

            CALL CHO_TIMER(C1,W1)
            DO IAB = 1,NAB
               KOFF = KINTQ + NNBSTR(ISYM,2)*(IAB - 1) - 1
               WRITE(LUSEL(ISYM)) (WORK(KOFF+I),I=1,NNBSTR(ISYM,2))
            END DO
            CALL CHO_TIMER(C2,W2)
            TINTEG(1,2) = TINTEG(1,2) + C2 - C1
            TINTEG(2,2) = TINTEG(2,2) + W2 - W1

         END IF
      END DO

      END
C  /* Deck cho_dal_nnbst */
      SUBROUTINE CHO_DAL_NNBST(LNNBST,NSYM)
C
C     Purpose: read array NNBST from ccsdsym.h.
C
#include <implicit.h>
      DIMENSION LNNBST(8)
#include <ccsdsym.h>

      DO ISYM = 1,NSYM
         LNNBST(ISYM) = NNBST(ISYM)
      END DO

      END
C  /* Deck cho_dal_n2bst */
      SUBROUTINE CHO_DAL_N2BST(LN2BST,NSYM)
C
C     Purpose: read array N2BST from ccsdsym.h.
C
#include <implicit.h>
      DIMENSION LN2BST(8)
#include <ccsdsym.h>

      DO ISYM = 1,NSYM
         LN2BST(ISYM) = N2BST(ISYM)
      END DO

      END
C  /* Deck cho_dal_iaodpk */
      SUBROUTINE CHO_DAL_IAODPK(LIAODPK,NSYM)
C
C     Purpose: read array IAODPK from symsq.h.
C
#include <implicit.h>
      DIMENSION LIAODPK(8,8)
#include <maxorb.h>
#include <symsq.h>

      DO ISYMB = 1,NSYM
         DO ISYMA = 1,NSYM
            LIAODPK(ISYMA,ISYMB) = IAODPK(ISYMA,ISYMB)
         END DO
      END DO

      END
C  /* Deck cho_dal_calcint1 */
      SUBROUTINE CHO_DAL_CALCINT1(ISHLA,ISHLB,WORK,LWORK,IPRINT,
     &                            INDRED,IQUVEA,IQUVEB,
     &                            MMBSTRT,CHO_MAXQUAL)
C
C     Purpose: calculate (**|AB). Result is stored on disk.
C              In addition, set up index arrays needed for reading
C              the integral buffers.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
      INTEGER   CHO_MAXQUAL
      DIMENSION INDRED(MMBSTRT,3)
      DIMENSION IQUVEA(CHO_MAXQUAL,8), IQUVEB(CHO_MAXQUAL,8)
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
      LOGICAL   SETUP
      DIMENSION INDEXA(MXCORB), INDEXB(MXCORB)
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <iratdef.h>

      DIMENSION JNDEXA(MXCORB), JNDEXB(MXCORB)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

C     Calculate integrals (stored on disk).
C     -------------------------------------

      DIACAL = .FALSE.
      SETUP  = .FALSE.
      CALL NERDI2(WORK,LWORK,INDEXA,INDEXB,ISHLA,ISHLB,
     &            NUMDIA,NUMDIB,IPRINT,SETUP,1)

C     Set shell addr. to full mapping according to our standard.
C     ----------------------------------------------------------

      CALL CHO_DAL_NEWINDX(JNDEXA,INDEXA,NUMDIA)
      CALL CHO_DAL_NEWINDX(JNDEXB,INDEXB,NUMDIB)

      IF (LOCDBG) THEN
         NUMA = 0
         DO I = 1,NUMDIA
            if (INDEXA(I) .ne. JNDEXA(I)) NUMA = NUMA + 1
         END DO
         NUMB = 0
         DO I = 1,NUMDIB
            IF (INDEXB(I) .ne. JNDEXB(I)) NUMB = NUMB + 1
         END DO
         IF (NUMA .eq. 0) THEN
          WRITE(LUPRI,*)
     &    'Cholesky and Hermit ordering identical for shell ',
     &    ISHLA
         ELSE
          WRITE(LUPRI,*)
     &     'Cholesky and Hermit ordering differ for shell ',
     &     ISHLA
         END IF
         IF (NUMB .EQ. 0) THEN
          WRITE(LUPRI,*)
     &    'Cholesky and Hermit ordering identical for shell ',
     &    ISHLB
         ELSE
          WRITE(LUPRI,*)
     &     'Cholesky and Hermit ordering differ for shell ',
     &     ISHLB
         END IF
      END IF

C     Set up indices for qualified columns.
C     -------------------------------------

      CALL CHO_DAL_GETQIN(ISHLA,ISHLB,INDRED,
     &                    IQUVEA,IQUVEB,JNDEXA,JNDEXB,
     &                    MMBSTRT,CHO_MAXQUAL,MXCORB)

      END
C  /* Deck cho_dal_newindx */
      SUBROUTINE CHO_DAL_NEWINDX(JNDX,INDX,NUM)
C
C     Purpose: Shift from Hermit to Cholesky indexation of shells.
C
#include <implicit.h>
      DIMENSION JNDX(*), INDX(NUM)
#include <cholesky.h>

      DO I = 1,NUM     ! addr. in shell (Hermit)
         J = INDX(I)   ! global AO
         K = ISHLSO(J) ! addr. in shell (Cholesky)
         JNDX(K) = J
      END DO

      END
C  /* Deck cho_dal_getqin */
      SUBROUTINE CHO_DAL_GETQIN(ISHLA,ISHLB,INDRED,
     &                          IQUVEA,IQUVEB,INDEXA,INDEXB,
     &                          MMBSTRT,LOC_MAXQUAL,LOC_MXCORB)
C
C     Purpose: set up index arrays for qualified columns.
C
#include <implicit.h>
      DIMENSION INDRED(MMBSTRT,3)
      DIMENSION IQUVEA(LOC_MAXQUAL,8), IQUVEB(LOC_MAXQUAL,8)
      DIMENSION INDEXA(LOC_MXCORB), INDEXB(LOC_MXCORB)
#include <cholesky.h>

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_GETQIN')

      INTEGER A,B

      DO ISYM = 1,NSYM

         NAB = NQUAL(ISYM) - IOFFQ(ISYM)
         IF (NAB .GT. LOC_MAXQUAL) THEN
            WRITE(LUPRI,*) SECNAM,': NAB,LOC_MAXQUAL: ',
     &                     NAB,LOC_MAXQUAL
            CALL CHO_QUIT('Dimension problem in '//SECNAM,104)
         END IF

         DO I = 1,NAB

            LAB = IOFFQ(ISYM) + I  ! qualified number
            KAB = IQUAB(LAB,ISYM)  ! addr. in rs2
            JAB = INDRED(KAB,2)    ! addr. in rs1
            IAB = INDRED(JAB,1)    ! addr. in sh. pair ISHLAB

            IF (ISHLA .EQ. ISHLB) THEN
               CALL CHO_INVPCK(IAB,IA,IB,.TRUE.)
            ELSE IF (ISHLA .GT. ISHLB) THEN
               CALL CHO_INVREC(IAB,IA,IB,NBSTSH(ISHLA),NBSTSH(ISHLB))
            ELSE
               CALL CHO_QUIT('ISHLA < ISHLB in '//SECNAM,104)
            END IF

            IQUVEA(I,ISYM) = INDEXA(IA)
            IQUVEB(I,ISYM) = INDEXB(IB)

         END DO

      END DO

      END
C  /* Deck cho_dal_rdao */
      SUBROUTINE CHO_DAL_RDAO(XINT,IQUVEA,IQUVEB,NAB,ISYMAB,
     &                        WORK,LWORK)
C
C     Purpose: read qualified integral columns from disk (full storage).
C
#include <implicit.h>
      DIMENSION XINT(*)
      DIMENSION IQUVEA(*),IQUVEB(*)
      DIMENSION WORK(LWORK)
#include <iratdef.h>
#include <mxcent.h>
#include <eribuf.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_RDAO')

      KIREC = 1
      KEND1 = KIREC + (NBUFX(0) - 1)/IRAT + 1
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL CHO_QUIT('Insufficient memory in '//SECNAM,101)
      END IF

      NEWDIS = .TRUE.
      CALL CC_RDAON(XINT,IQUVEA,IQUVEB,NAB,ISYMAB,WORK(KIREC),
     &              WORK(KEND1),LWRK1)

      END
C  /* Deck cho_dal_rs2full */
      SUBROUTINE CHO_DAL_RS2FULL(IRS2F,ISYM,IIBSTRSH,NNBSTRSH,INDRED,
     &                           MSYM,MMSHL,MMBSTRT)
C
C     Purpose: set up map from current reduced set to full storage
C              for symmetry block ISYM.
C
#include <implicit.h>
      DIMENSION IRS2F(*)
      DIMENSION IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      DIMENSION INDRED(MMBSTRT,3)
#include <choorb.h>
#include <cholesky.h>

      CHARACTER*15 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_RS2FULL')

      INTEGER IAODPK(8,8)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER CHO_DAL_RS2F

      ITRI(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J

      CALL CHO_IZERO(IRS2F,NNBSTR(ISYM,2))
      CALL CHO_DAL_IAODPK(IAODPK,NSYM)

      DO ISYMB = 1,NSYM
         IF (NBAS(ISYMB) .GT. 0) THEN
            ISYMA = MULD2H(ISYMB,ISYM)
            IF (NBAS(ISYMA) .GT. 0) THEN
               IF (ISYMA .GT. ISYMB) THEN
                  DO KB = 1,NBAS(ISYMB)
                     LB = ISHLSO(IBAS(ISYMB)+KB)
                     ISHLB = ISOSHL(IBAS(ISYMB)+KB)
                     DO KA = 1,NBAS(ISYMA)
                        LA = ISHLSO(IBAS(ISYMA)+KA)
                        ISHLA = ISOSHL(IBAS(ISYMA)+KA)
                        ISHLAB = ITRI(ISHLA,ISHLB)
                        IF (ISHLA .LT. ISHLB) THEN
                           LAB = NBSTSH(ISHLB)*(LA - 1) + LB
                        ELSE IF (ISHLA .EQ. ISHLB) THEN
                           LAB = ITRI(LA,LB)
                        ELSE
                           LAB = NBSTSH(ISHLA)*(LB - 1) + LA
                        END IF
                        I = CHO_DAL_RS2F(LAB,ISHLAB,ISYM,
     &                                   IIBSTRSH,NNBSTRSH,INDRED,
     &                                   MSYM,MMSHL,MMBSTRT)
                        IF (I .GT. 0) THEN
#if defined (CHODAL_DEBUG)
                           IF (IRS2F(I).NE.0) THEN
                              WRITE(LUPRI,*) 
     &                                  SECNAM,': ooops: rs element ',i,
     &                                  ' has already been set!'
                              CALL CHO_QUIT('Error detected in '
     &                                      //secnam,103)
                           END IF
#endif
                           IRS2F(I) = IAODPK(ISYMA,ISYMB)
     &                              + NBAS(ISYMA)*(KB - 1) + KA
#if defined (CHODAL_DEBUG)
                           KAG = IBAS(ISYMA) + KA
                           KBG = IBAS(ISYMB) + KB
                           II  = INDRED(I,2) ! addr. in rs1
                           IF (IRS1TOF(1,II).NE.KAG .OR.
     &                         IRS1TOF(2,II).NE.KBG) THEN
                              WRITE(LUPRI,*) SECNAM,': index error!!'
                              WRITE(LUPRI,*)
     &                        'Addr. in rs1 and current: ',II,I
                              WRITE(LUPRI,*)
     &                        'alpha and beta from IRS2F: ',KAG,KBG
                              WRITE(LUPRI,*)
     &                        'alpha and beta should be : ',
     &                        IRS1TOF(1,II),IRS1TOF(2,II)
                              CALL CHO_QUIT('Error detected in '
     &                                      //SECNAM,103)
                           END IF
#endif
                        END IF
                     END DO
                  END DO
               ELSE IF (ISYMA .EQ. ISYMB) THEN
                  DO KA = 1,NBAS(ISYMA)
                     LA = ISHLSO(IBAS(ISYMA)+KA)
                     ISHLA = ISOSHL(IBAS(ISYMA)+KA)
                     DO KB = 1,KA
                        LB = ISHLSO(IBAS(ISYMB)+KB)
                        ISHLB = ISOSHL(IBAS(ISYMB)+KB)
                        ISHLAB = ITRI(ISHLA,ISHLB)
                        IF (ISHLA .LT. ISHLB) THEN
                           LAB = NBSTSH(ISHLB)*(LA - 1) + LB
                        ELSE IF (ISHLA .EQ. ISHLB) THEN
                           LAB = ITRI(LA,LB)
                        ELSE
                           LAB = NBSTSH(ISHLA)*(LB - 1) + LA
                        END IF
                        I = CHO_DAL_RS2F(LAB,ISHLAB,ISYM,
     &                                   IIBSTRSH,NNBSTRSH,INDRED,
     &                                   MSYM,MMSHL,MMBSTRT)
                        IF (I .GT. 0) THEN
#if defined (CHODAL_DEBUG)
                           IF (IRS2F(I).NE.0) THEN
                              WRITE(LUPRI,*) 
     &                                  SECNAM,': ooops: rs element ',i,
     &                                  ' has already been set!'
                              CALL CHO_QUIT('Error detected in '
     &                                      //secnam,103)
                           END IF
#endif
                           IRS2F(I) = IAODPK(ISYMA,ISYMB) + ITRI(KA,KB)
#if defined (CHODAL_DEBUG)
                          KAG = IBAS(ISYMA) + KA
                          KBG = IBAS(ISYMB) + KB
                          II  = INDRED(I,2) ! addr. in rs1
                          IF (IRS1TOF(1,II).NE.KAG .OR.
     &                        IRS1TOF(2,II).NE.KBG) THEN
                             WRITE(LUPRI,*) SECNAM,': index error!!'
                             WRITE(LUPRI,*)
     &                       'Addr. in rs1 and current: ',II,I
                             WRITE(LUPRI,*)
     &                       'alpha and beta from IRS2F: ',KAG,KBG
                             WRITE(LUPRI,*)
     &                       'alpha and beta should be : ',
     &                       IRS1TOF(1,II),IRS1TOF(2,II)
                             CALL CHO_QUIT('Error detected in '
     &                                     //SECNAM,103)
                          END IF
#endif
                        END IF
                     END DO
                  END DO
               END IF
            END IF
         END IF
      END DO

      IF (LOCDBG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) 'Debug info from CHO_DAL_RS2FULL: ',
     &                  'mapping from current red.set to full'
         WRITE(LUPRI,*) 'Symmetry block: ',ISYM,
     &                  ' Reduced dimension: ',NNBSTR(ISYM,2)
         ICOUNT = 0
         DO I = 1,NNBSTR(ISYM,2)
            IF (IRS2F(I) .EQ. 0) ICOUNT = ICOUNT + 1
         END DO
         IF (ICOUNT .EQ. 0) THEN
            WRITE(LUPRI,*) '#elements not set: ',ICOUNT,' (OK)'
            WRITE(LUPRI,*)
            CALL CHO_FLUSH(LUPRI)
         ELSE
            WRITE(LUPRI,*) '#elements not set: ',ICOUNT,' (OH DEAR!!!)'
            CALL CHO_QUIT('Setup error detected in '//SECNAM,103)
         END IF
      END IF

      END
C  /* Deck cho_dal_rs2f */
      INTEGER FUNCTION CHO_DAL_RS2F(LAB,ISHLAB,ISYMAB,
     &                              IIBSTRSH,NNBSTRSH,INDRED,
     &                              MSYM,MMSHL,MMBSTRT)
C
C     Purpose: return index in current reduced set, 0 if not there.
C
C     LAB is the index in shell pair ISHLAB.
C
#include <implicit.h>
      DIMENSION IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      DIMENSION INDRED(MMBSTRT,3)
#include <cholesky.h>

      CHO_DAL_RS2F = 0

      K1 = IIBSTR(ISYMAB,2) + IIBSTRSH(ISYMAB,ISHLAB,2) + 1
      K2 = K1 + NNBSTRSH(ISYMAB,ISHLAB,2) - 1
      DO K = K1,K2
         IF (INDRED(INDRED(K,2),1) .EQ. LAB) THEN
            CHO_DAL_RS2F = K - IIBSTR(ISYMAB,2)
            RETURN
         END IF
      END DO

      END
C  /* Deck cho_dal_getqcol */
      SUBROUTINE CHO_DAL_GETQCOL(XINT,XINTR,IRS2F,NX,NXR,NCOL)
C
C     Purpose: get qualified integral columns sorted in current
C              reduced set.
C
#include <implicit.h>
      DIMENSION XINT(NX,NCOL), XINTR(NXR,NCOL)
      DIMENSION IRS2F(NXR)

      DO ICOL = 1,NCOL
         DO I = 1,NXR
            XINTR(I,ICOL) = XINT(IRS2F(I),ICOL)
         END DO
      END DO

      END
C  /* Deck cho_dal_setrs1tof */
      SUBROUTINE CHO_DAL_SETRS1TOF(IIBSTRSH,NNBSTRSH,INDRED,
     &                             MSYM,MMSHL,MMBSTRT)
C
C     Purpose: set mapping from first reduced set to full storage,
C
C              IRS1TOF(1,I) = global alpha corresponding to element I of
C                             rs1.
C              IRS1TOF(1,I) = global beta corresponding to element I of
C                             rs1.
C
C     NB! This routine requires that the CHODAL_DEBUG precompiler flag
C         has been set.
C
#include <implicit.h>
      DIMENSION IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
      DIMENSION INDRED(MMBSTRT,3)
#include <choorb.h>
#include <cholesky.h>

      CHARACTER*17 SECNAM
      PARAMETER (SECNAM = 'CHO_DAL_SETRS1TOF')

#if defined (CHODAL_DEBUG)
      INTEGER A, B
      INTEGER ICHECK(CHO_MXNNBST)
      ITRI(I,J) = MAX(I,J)*(MAX(I,J)-3)/2+I+J

      CALL IZERO(IRS1TOF,2*CHO_MXNNBST)
      NNSHL = NSHELL*(NSHELL+1)/2
      DO ISHLAB = 1,NNSHL
         CALL CHO_INVPCK(ISHLAB,ISHLA,ISHLB,.TRUE.)
         DO ISYM = 1,NSYM
            I1 = IIBSTR(ISYM,1) + IIBSTRSH(ISYM,ISHLAB,1) + 1
            I2 = I1 + NNBSTRSH(ISYM,ISHLAB,1) - 1
            DO I = I1,I2
               IAB = INDRED(I,1) ! addr in shell pair
               IF (ISHLA .EQ. ISHLB) THEN
                  CALL CHO_INVPCK(IAB,IA,IB,.TRUE.)
               ELSE
                  CALL CHO_INVREC(IAB,IA,IB,NBSTSH(ISHLA),
     &                                      NBSTSH(ISHLB))
               END IF
               A = ISHL2F(IA,ISHLA) ! global AO index, alpha
               B = ISHL2F(IB,ISHLB) ! global AO index, beta
               IRS1TOF(1,I) = A
               IRS1TOF(2,I) = B
            END DO
         END DO
      END DO
      ICOUNT = 0
      DO I = 1,NNBSTRT(1)
         IF (IRS1TOF(1,I).LT.1 .OR. IRS1TOF(2,I).LT.1 .OR.
     &       IRS1TOF(1,I).GT.NBAST .OR. IRS1TOF(2,I).GT.NBAST) THEN
            WRITE(LUPRI,*) SECNAM,': RS1 address, alpha, beta: ',
     &                    I,IRS1TOF(1,I),IRS1TOF(2,I),' (out of bounds)'
            ICOUNT = ICOUNT + 1
         END IF
      END DO
      IF (ICOUNT .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': ',ICOUNT,' setup errors detected'
         WRITE(LUPRI,*) SECNAM,': alpha,beta must be in range 1 to ',
     &                  NBAST
         CALL CHO_QUIT('Setup error [1] in '//SECNAM,103)
      END IF
      CALL IZERO(ICHECK,CHO_MXNNBST)
      ICOUNT = 0
      DO I = 1,NNBSTRT(1)
         A = IRS1TOF(1,I)
         B = IRS1TOF(2,I)
         J = ITRI(A,B)
         IF (ICHECK(J) .EQ. 0) THEN
            ICHECK(J) = 1
         ELSE
            WRITE(LUPRI,*) SECNAM,': alpha,beta ',A,B,' already set!',
     &                     '   RS1 addr: ',I
            ICOUNT = ICOUNT + 1
         END IF
      END DO
      IF (ICOUNT .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': ',ICOUNT,
     &                  ' elements set more than once!'
         CALL CHO_QUIT('Setup error [2] in '//SECNAM,103)
      END IF
#else
      WRITE(LUPRI,*) SECNAM,' was called without the CHODAL_DEBUG ',
     &               ' precompiler flag set. Returning.'
      RETURN
#endif

      END
C  /* Deck cho_getvecelm */
      SUBROUTINE CHO_GETVECELM(ISYM,A,B,IADR,INDRED,MMBSTRT)
C
C     Purpose: return address of specified (global AO indices) element
C              in current reduced set, sym. block ISYM.
C              IADR=0 is returned if not found.
C
C     NB! This routine requires that the CHODAL_DEBUG precompiler flag
C         has been set.
C
#include <implicit.h>
      INTEGER   A, B
      DIMENSION INDRED(MMBSTRT,3)
#include <cholesky.h>

#if defined (CHODAL_DEBUG)
      IADR = 0
      DO I = 1,NNBSTR(ISYM,2)
         J = INDRED(IIBSTR(ISYM,2)+I,2)
         IF ((IRS1TOF(1,J).EQ.A .AND. IRS1TOF(2,J).EQ.B) .OR.
     &       (IRS1TOF(1,J).EQ.B .AND. IRS1TOF(2,J).EQ.A)) THEN
            IADR = I
            RETURN
         END IF
      END DO
#else
      CHARACTER*13 SECNAM
      PARAMETER (SECNAM = 'CHO_GETVECELM')

      IADR = 0
      WRITE(LUPRI,*) SECNAM,' was called without the CHODAL_DEBUG ',
     &               ' precompiler flag set. Returning.'
      RETURN
#endif

      END
C  /* Deck cho_dal_final */
      SUBROUTINE CHO_DAL_FINAL(XLTHRC,DAMP,XLSPAN,IIBSTR,NNBSTR,LNUMCHO,
     &                         LRSTDIA,LRSTCHO,NSYM)
C
C     Purpose: transfer information to ccdeco.h about decomposition.
C
#include <implicit.h>
      DIMENSION DAMP(2)
      DIMENSION IIBSTR(8,3), NNBSTR(8,3), LNUMCHO(8)
      LOGICAL   LRSTDIA, LRSTCHO
#include <maxorb.h>
#include <ccdeco.h>
#include <choscf.h>

      THRCOM = XLTHRC
      THINDI = DAMP(1)
      THSUDI = DAMP(2)
      SPAN   = XLSPAN
      LREDU  = 0
      DO ISYM = 1,NSYM
         IDIAG(ISYM)  = IIBSTR(ISYM,1)
         NREDUC(ISYM) = NNBSTR(ISYM,1)
         LREDU = LREDU + NNBSTR(ISYM,1)
      END DO
      NDIAG = LREDU
      DO ISYM = 1,NSYM
         NUMCHO(ISYM) = LNUMCHO(ISYM)
      END DO
      RSTDIA = LRSTDIA
      RSTCHO = LRSTCHO

      CHOINT = .TRUE.

      END
C  /* Deck cho_dal_chov2io */
      SUBROUTINE CHO_DAL_CHOV2IO()
C
C     Purpose: set info in chov2io.h
C
      IMPLICIT NONE
#include <chov2io.h>

      CALL IZERO(CHOV2I,3)

      END
C  /* Deck set_choscfdef */
      SUBROUTINE SET_CHOSCFDEF()
C
C     Purpose: set defaults for Cholesky SCF.
C
#include <implicit.h>
#include <choscf.h>

      CCMODSK = .TRUE.
      LUCCMO  = 79
      FCCMO   = 'CHOCMO'
      THRDC   = 1.00D-12
      CALL IZERO(NDVCS,8)

      END
