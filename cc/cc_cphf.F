c /* deck cc_cphf */
*=====================================================================*
      SUBROUTINE CC_CPHF(TYPE,LABEL,ISYMS,ISTAT,EIGV,
     &                   ISYMO,FREQS,ICAU,NVEC,MAXVEC,
     &                   LORX,WORK,LWORK)
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
#  include <priunit_dec.h>
#  include <dummy_dec.h>
#  include <exeinf_dec.h>
#  include <maxorb_dec.h>
#  include <infvar_dec.h>
#  include <iratdef_dec.h>
#  include <ccsdsym_dec.h>
#  include <ccorb_dec.h>
#  include <cclr_dec.h>
#  include <inftap_dec.h>
      INTEGER IPRLIN, LSYMRF, LSYMPT, LSYMST,
     &        NCONRF, NCONST, NWOPPT, NVARPT
#else
#  include <implicit.h>
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

* local parameters:
      CHARACTER*(18) MSGDBG
      PARAMETER (MSGDBG = '[debug] CC_CPHF>  ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE. )
      INTEGER KDUMMY
      PARAMETER (KDUMMY = -99 999 999)
      CHARACTER*8 FNGDVE, FNSOVE, FNREVE
      PARAMETER(FNGDVE='CCPHFRHS',FNSOVE='CCPHFSOL')
      PARAMETER(FNREVE='CCPHFRES')
      INTEGER MXFRVAL
      PARAMETER (MXFRVAL = 100)
C     PARAMETER (MXFRVAL = 1)  ! switch off simultaneous equations

      CHARACTER TYPE*(*)

      INTEGER NVEC, MAXVEC, LWORK
      LOGICAL LORX(MAXVEC,*)
      INTEGER ISYMS(MAXVEC,*), ISYMO(MAXVEC,*)
      INTEGER ISTAT(MAXVEC,*), ICAU(MAXVEC,*)

      CHARACTER*8 LABEL(MAXVEC,*)

#if defined (SYS_CRAY)
      REAL FREQS(MAXVEC,*), EIGV(MAXVEC,*)
      REAL WORK(LWORK)
      REAL ZERO, RDUM, DUMMY
      REAL XNORM, DDOT, DSQRT, GDNORM
#else
      DOUBLE PRECISION FREQS(MAXVEC,*), EIGV(MAXVEC,*)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION ZERO, RDUM
      DOUBLE PRECISION XNORM, DDOT, DSQRT, GDNORM
#endif
      PARAMETER (ZERO = 0.0d0)

      CHARACTER MODEL*(10), MODWR*(10)
      LOGICAL MAYBE_MORE, NEED_CPHF, RELAX, CICLC, HFCLC, TRPCLC,OOTV
      LOGICAL EXCLC, FASAVE, OPTST, EX
      INTEGER IOPT, ISYM, IVEC, NSTAT, ORDER, IDUM, IPERT,LWRK0,LWRK1
      INTEGER NCMOT,NASHT,N2ASHX,LCINDX,KCMO,KUDV,KXINDX,KEND0,KEND1
      INTEGER LUGDVE, LUSOVE, MXRM, MXPHP, NCOSAV, IOPTWR, KSLVEC
      INTEGER ISYMAT, NFRVAL, KFRVAL, NABAOP, NABATY, IFRVAL, IDX
      INTEGER IPRABA, LUREVE

* external functions:
      INTEGER ILSTSYM


*---------------------------------------------------------------------*
* do some checks:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE(LUPRI,*) 'Entered CC_CPHF.'
         CALL FLSHFO(LUPRI)
      END IF

* check vector type:
      IF (TYPE(1:2).EQ.'R1') THEN
        ORDER = 1
        NSTAT = 0
        MODWR  = 'SCF?      '
      ELSE
        WRITE (LUPRI,*) 'CPHF vectors ',TYPE(1:2),' not implemented.'
        CALL QUIT('required CPHF vectors not implemented.')
      END IF

* check if there is anything to do:
      NEED_CPHF = .FALSE.
      DO IVEC = 1, NVEC
         RELAX = .TRUE.
         DO IPERT = 1, ORDER
            RELAX = RELAX .AND. LORX(IVEC,IPERT)
         END DO
         NEED_CPHF = NEED_CPHF .OR. RELAX
      END DO   

      IF (.NOT. NEED_CPHF) RETURN

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
        WRITE(LUPRI,*) 'CCRHSVEC Workspace:',LWORK
      END IF
  
      CALL FLSHFO(LUPRI)

*---------------------------------------------------------------------*
* some initilizations:
*---------------------------------------------------------------------*

* arrays for GETGP and ABARSP:
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
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',KDUMMY,
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

* flags for CC_WRCPHF routine:
      IOPTWR = 4

* open file for right hand side and solution vectors:

      LUGDVE = 0
      LUSOVE = 0
      LUREVE = 0
      CALL GPOPEN(LUGDVE,FNGDVE,'UNKNOWN',' ','UNFORMATTED',KDUMMY,
     &            .FALSE.)
      CALL GPOPEN(LUSOVE,FNSOVE,'UNKNOWN',' ','UNFORMATTED',KDUMMY,
     &            .FALSE.)
      CALL GPOPEN(LUREVE,FNREVE,'UNKNOWN',' ','UNFORMATTED',KDUMMY,
     &            .FALSE.)
* open property integral file:
      IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                               'UNFORMATTED',KDUMMY,.FALSE.)

* close twoelectron integral file 'AOTWOINT':
      IF (LUINTA.GT.0) CALL GPCLOSE(LUINTA,'KEEP')

*---------------------------------------------------------------------*
* loop over vectors, set up rhs vectors and call ABARSP to solve eqs.:
*---------------------------------------------------------------------*

      IVEC = 0
      DO WHILE(IVEC.LT.NVEC)
        IVEC  = IVEC + 1

        RELAX = .TRUE.
        DO IPERT = 1, ORDER
           RELAX = RELAX .AND. LORX(IVEC,IPERT)
        END DO

        IF (RELAX) THEN

           ISYM   = ILSTSYM(TYPE,IVEC)
C          CALL LNRVAR(ISYM,IPRINT,WORK(KEND0),LWRK0)


           IF (LWRK0 .LT. NALLAI(ISYM)) THEN
              CALL QUIT('Insufficient memory in CC_CPHF.')
           END IF

C          ---------------------------
C          get right hand side vector:
C          ---------------------------
           CALL CC_GETHFGD(IVEC,TYPE,LABEL,ISYMS,ISTAT,EIGV,ISYMO,
     &                     FREQS,LORX,ICAU,NVEC,MAXVEC,ISYMAT,
     &                     WORK(KCMO),WORK(KUDV),WORK(KXINDX),
     &                     WORK(KFRVAL),WORK(KEND0),LWRK0)

           IF (LOCDBG) THEN
             WRITE (LUPRI,*) 'CC_CPHF> GP vector, label:',
     &                        LABEL(IVEC,1)
             WRITE (LUPRI,'(5X,I5,F12.8)') (I,WORK(KEND0-1+I),
     &                        I=1,NALLAI(ISYM))
           END IF


C          ------------------------------
C          write gradient vector to file:
C          ------------------------------
           REWIND LUGDVE
           CALL WRITT(LUGDVE,NALLAI(ISYM),WORK(KEND0))


C          ----------------------------------
C          check norm of the gradient vector:
C          ----------------------------------
           GDNORM=DSQRT(DDOT(NALLAI(ISYM),WORK(KEND0),1,WORK(KEND0),1))
           IF (LOCDBG) WRITE (LUPRI,*) 'GDNORM:',GDNORM


C          --------------------------------------------------------
C          for 'R1' check if several frequencies for same operator:
C          --------------------------------------------------------
           NFRVAL     = 1
           MAYBE_MORE = .TRUE.
           DO WHILE (MAYBE_MORE) 
              IF (TYPE(1:2).EQ.'R1' 
     &            .AND. IVEC.LT.NVEC .AND. NFRVAL.LT.MXFRVAL
     &            .AND. LABEL(IVEC,1).EQ.LABEL(IVEC+1,1) ) THEN
                WORK(KFRVAL+NFRVAL) = FREQS(IVEC+1,1)
                IVEC   = IVEC + 1
                NFRVAL = NFRVAL + 1
              ELSE
                MAYBE_MORE = .FALSE.
              END IF
           END DO


           IF (GDNORM.GT.THRLDPHF) THEN

C            -----------------------
C            solve CPHF equation(s):
C            -----------------------
             IPRABA = IPRINT
             NABAOP = 1
             NABATY = ISYMAT  ! flag for real/imaginary operators
             CALL ABARSP(CICLC,HFCLC,TRPCLC,OOTV,ISYM,EXCLC,
     &                   WORK(KFRVAL),NFRVAL,NABATY,NABAOP,
     &                   LABEL(IVEC,1),LUGDVE,LUSOVE,LUREVE,THRLEQ,
     &                   MAXITE,IPRABA,MXRM,MXPHP,WORK(KEND0),LWRK0)

C            --------------------------------------------
C            clean up `left overs' from ABARSP:
C            --------------------------------------------
             CALL GPCLOSE(LUSIFC,'KEEP')
             IF (LUONEL .GT. 0) CALL GPCLOSE(LUONEL,'KEEP')
             IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
C
C            --------------------------------------------
C            read solution vector(s) and save on CC file:
C            --------------------------------------------
             REWIND LUSOVE

           END IF

           DO IFRVAL = 1, NFRVAL
              KSLVEC = KEND0
              KEND1  = KSLVEC + 2*NALLAI(ISYM)
              LWRK1  = LWORK  - KEND1
              IF (LWRK1.LT.0) THEN
                 CALL QUIT('Insufficient memory in CC_CPHF.')
              END IF

              IF (GDNORM.GT.THRLDPHF) THEN
                CALL  READT(LUSOVE,2*NALLAI(ISYM),WORK(KSLVEC))
              ELSE
                CALL DZERO(WORK(KSLVEC),2*NALLAI(ISYM))
              END IF

              IDX = IVEC - NFRVAL + IFRVAL

              CALL CC_WRRSP(TYPE,IDX,ISYM,IOPTWR,MODWR,
     &                      WORK(KSLVEC),DUMMY,DUMMY,WORK(KEND1),LWRK1)

              IF (LOCDBG) THEN
                 WRITE (LUPRI,*) 
     &              'CC_CPHF> solution vector, label, freq:',
     &              LABEL(IVEC,1),WORK(KFRVAL)
                 WRITE (LUPRI,'(5X,I5,F12.8)')
     &              (I,WORK(KSLVEC-1+I),I=1,2*NALLAI(ISYM))
                 WRITE (LUPRI,*) 
     &              'CC_CPHF> saved CPHF solution for ',TYPE,
     &              ' equation nb. ',IDX
              END IF

           END DO

        END IF

      END DO
*---------------------------------------------------------------------*
* that's it: close files, restore variables and return
*---------------------------------------------------------------------*
C
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

      RETURN
      END
*=====================================================================*
*              END OF SUBROUTINE CCRHSVEC                             *
*=====================================================================*
c /* deck cc_sirinf */
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
