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
#include "single.h"

C  /* Deck twoint */

      SUBROUTINE TWOINT(WORK,LWORK,HESSEE,FMAT,DMAT,NDMAT,IREPDM,
     &                  IFCTYP,GMAT,
     &                  INDXAB,NUMDIS,MAXDIS,ITYPE,MAXDIF,JATOM,NODV,
     &                  NOPV,NOCONT,TTIME,JPRINT,IPRNTA,IPRNTB,IPRNTC,
     &                  IPRNTD,RETUR,ISHLA,I2TYP,JSTRSH,NPRIMS,NCONTS,
     &                  IORBSH,ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,
     &                  DINTSKP,RELCAL)
C
C     Trygve Helgaker, University of Oslo, Norway, 1984
C
C     General contraction, Feb - Mar 1988, TUH
C     Symmetry processing 880408  PRT & TUH
C     Screening March 1997 T.Saue
C
C     References for calculation of Cartesian integrals:
C     --------------------------------------------------
C
C     L. E. McMurchie & E. R. Davidson, J. Comput. Chem. 26, 218 (1978)
C     V. R. Saunders, in "Methods in Computational Molecular Physics",
C       G. H. F. Diercksen and S. Wilson, eds. (Reidel,Dordrecht,1983)
C     T. U. Helgaker et al., JCP 84, 6266 (1986)
C
C
C     References for symmetry:
C     ------------------------
C
C     E. R. Davidson, JCP 62, 400 (1975)
C     P. R. Taylor,   TCA 69, 447 (1986)
C     J. Almloef, MOLECULE Program Description, USIP Report 74 - 29
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "dummy.h"
      PARAMETER (D0 = 0.00 D00)
      LOGICAL NOPV, NODV, PERTUR, EXPECT, UNDIFF, DIRFCK, RETUR, TTIME,
     &        NOCONT, SPNORB, DIA2SO, ZFS2EL, DISTRI, SQ12EL, LONDON,
     &        SUSCEP, DDFOCK, ADISTR,SOFOCK,RELCAL
      DIMENSION HESSEE(*),DMAT(*), FMAT(*), GMAT(*), IREPDM(*),
     &        IFCTYP(*),DINTSKP(*),GABRAO(*),DMRAO(*),DMRSO(*),
     &        JSTRSH(*), NPRIMS(*), NCONTS(*), IORBSH(*),
     &        INDXAB(MAXDIS), WORK(LWORK)
#include "ccom.h"
#include "symmet.h"
#include "nuclei.h"
#include "twosta.h"
#include "infpar.h"
#include "blocks.h"
C
Cholesky
#include "ccdeco.h"
Cholesky  
C
      CALL QENTER('TWOINT')
      IF (JPRINT .GT. 5) CALL TITLER('Output from TWOINT','*',103)
C
C980827-hjaaj : "old" NOCONT does not work in this version
C
      IF (NOCONT) THEN
         WRITE (LUPRI,*)
     &      'TWOINT ERROR, .NOCONT is not implemented in this version'
         CALL QUIT('TWOINT: .NOCONT is not implemented in this version')
      END IF
C
      TKTIME = TTIME
C
C     CPU time statistics
C
      CALL GETTIM(T2INT,DUMTIM)
      CALL DELSTA(-1,0)
      CALL REDSTA(-1,0,0,0,0)
      IF (TKTIME) THEN
         TR000  = D0
         THERI  = D0
         TODCVE = D0
         TEXCOE = D0
         TC10IN = D0
         TC1EIN = D0
         TC2HIN = D0
         TC2EIN = D0
         TSYM2S = D0
         TSYMOU = D0
         TFCKOU = D0
         TDSOUT = D0
         TSPOSY = D0
         TDRSYM = D0
         TINTEX = D0
         TPATH1 = D0
         TPATH2 = D0
         MAXJ = 4*(NHTYP - 1) + 2
         DO 100 I = 0, MAXJ
            TR000X(I) = D0
            THERIX(I) = D0
  100    CONTINUE
      END IF
C
C     Work space statistics
C
      LWTOT  = 0
      MWTOT  = 0
      MWFCAB = 0
      MWFCCD = 0
      MWPSO  = 0
      MWSOIN = 0
      MWAOIN = 0
      MWHRIN = 0
      MWHRND = 0
      MWHCIN = 0
      MWRJ00 = 0
      MWHRSQ = 0
      MWC1DR = 0
      MWC2DR = 0
      MWC2HI = 0
      MWC2EI = 0
      MWINTE = 0
      MWPPRI = 0
      MWDRSY = 0
C
C     Determine run type
C
      CALL TWORUN(ITYPE,MAXDIF,JATOM,PERTUR,EXPECT,UNDIFF,
     &            DIRFCK,SOFOCK,SPNORB,DIA2SO,ZFS2EL,DISTRI,SQ12EL,
     &            LONDON,SUSCEP,DDFOCK,ADISTR,MAXDER,IATOM,MULE,MULTE,
     &            JPRINT)
C
C     Threshold for integrals
C
      THRESH = MAX(THRS,1.00D-15)
C
Cholesky
      IF (CHOINT) THRESH = 1.0D-40    ! Probably not needed and surely too much
Cholesky
C
C     *******************************
C     ***** Calculate integrals *****
C     *******************************
C
      CALL TWOCAL(WORK,LWORK,HESSEE,FMAT,DMAT,NDMAT,GMAT,INDXAB,MAXDER,
     &            EXPECT,UNDIFF,DDFOCK,DIRFCK,SOFOCK,DISTRI,LONDON,
     &            SPNORB,DIA2SO,ZFS2EL,SUSCEP,PERTUR,IATOM,MULE,MULTE,
     &            NODV,NOPV,NOCONT,
     &            THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,NUMDIS,
     &            MAXDIS,SQ12EL,RETUR,NHTYP,IREPDM,IFCTYP,ISHLA,ADISTR,
     &            ITYPE,JSTRSH,NPRIMS,NCONTS,IORBSH,
     &            I2TYP,ICEDIF,IFTHRS,DINTSKP,GABRAO,DMRAO,DMRSO,
     &            RELCAL)
      CALL GETTIM(TIMEND,DUMTIM)
      T2INT = TIMEND - T2INT
C
C     *******************
C     ***** Timings *****
C     *******************
C
      IF (TKTIME) THEN
         CALL HEADER('Timings for TWOINT',1)
         WRITE (LUPRI, '(18(/,1X,A,F8.2,A))')
     &      ' Time in R000:   ', TR000, ' seconds',
     &      ' Time in HERI:   ', THERI, ' seconds',
     &      ' Time in ODCVEC: ', TODCVE,' seconds',
     &      ' Time in EXCOEF: ', TEXCOE,' seconds',
     &      ' Time in C10INT: ', TC10IN,' seconds',
     &      ' Time in C1EINT: ', TC1EIN,' seconds',
     &      ' Time in C2HINT: ', TC2HIN,' seconds',
     &      ' Time in C2EINT: ', TC2EIN,' seconds',
     &      ' Time in SYM2:   ', TSYM2S,' seconds',
     &      ' Time in DRSYM2: ', TDRSYM,' seconds',
     &      ' Time in SPOSYM: ', TSPOSY,' seconds',
     &      ' Time in UN2OUT: ', TSYMOU,' seconds',
     &      ' Time in FCKOUT: ', TFCKOU,' seconds',
     &      ' Time in DR2OUT: ', TDSOUT,' seconds',
     &      ' Time in INTEXP: ', TINTEX,' seconds',
     &      ' Time for PATH1: ', TPATH1,' seconds',
     &      ' Time for PATH2: ', TPATH2,' seconds'
C        WRITE (LUPRI, '(/,2(1X,A,I2,A,/,1X,10F6.2))')
C    &      ' Detailed timings for R000 integrals (JMAX = 1, ',
C    &        MAXJ,'):', (TR000X(I), I = 1, MAXJ),
C    &      ' Detailed timings for Hermitian integrals (JMAX = 1, ',
C    &        MAXJ,'):', (THERIX(I), I = 1, MAXJ)
      END IF
C
C     *********************************
C     ***** Work space statistics *****
C     *********************************
C
      IF (LWTOT .NE. 0 .OR.
     &    (JPRINT .GT. 2 .AND. .NOT. SLAVE)  .OR.
     &    (JPRINT .GT. 5 .AND. SLAVE) ) THEN
      IF (DISTRI .AND. (JPRINT.LE.4)) THEN
         IF (NUMDIS .NE. -1) THEN
            WRITE (LUPRI,'(/1X,A,I12)')
     &      ' Work space in TWOINT for these distributions:',MWTOT
         END IF
      ELSE
         CALL HEADER('Maximum work space allocations in TWOINT',1)
         WRITE (LUPRI,'(15(1X,A,I12/))')
     &      ' Total allocation:                   ', MWTOT,
     &      ' Expansion coefficients (electron 1):', MWFCAB,
     &      ' Expansion coefficients (electron 2):', MWFCCD,
     &      ' Two-electron densities:             ', MWPSO ,
     &      ' SO integrals:                       ', MWSOIN,
     &      ' AO integrals:                       ', MWAOIN,
     &      ' Hermite integrals:                  ', MWHRIN,
     &      ' Hermite integral pointer:           ', MWHRND,
     &      ' Hermite-Cartesian integrals:        ', MWHCIN,
     &      ' Gamma functions:                    ', MWRJ00,
     &      ' Work space for Hermite integrals:   ', MWHRSQ,
     &      ' Allocations in C1DRIV:              ', MWC1DR,
     &      ' Allocations in C2DRIV:              ', MWC2DR,
     &      ' Allocations in C2HINT:              ', MWC2HI,
     &      ' Allocations in C2EINT:              ', MWC2EI,
     &      ' Allocations in INTEXP:              ', MWINTE,
     &      ' Allocations in PPRIM :              ', MWPPRI,
     &      ' Allocations in DRSYM2:              ', MWDRSY
      END IF
      END IF
      IF (LWTOT .NE. 0) THEN
         WRITE (LUPRI,'(1X,A)')     ' Error in work space statistics.'
         WRITE (LUPRI,'(1X,A,I12)') ' LWTOT at end of TWOINT:',LWTOT
      END IF
      IF ((JPRINT .GT. 2 .AND. .NOT. SLAVE)  .OR.
     &    (JPRINT .GT. 5 .AND. SLAVE) ) THEN
         CALL TWOPER
         CALL REDSTA(1,0,0,0,0)
         CALL DELSTA(1,0)
      END IF
      CALL QEXIT('TWOINT')
      RETURN
      END
C  /* Deck twocal */
      SUBROUTINE TWOCAL(WORK,LWORK,HESSEE,FMAT,DMAT,NDMAT,GMAT,INDXAB,
     &            MAXDER,
     &            EXPECT,UNDIFF,DDFOCK,DIRFCK,SOFOCK,DISTRI,LONDON,
     &            SPNORB,DIA2SO,ZFS2EL,SUSCEP,PERTUR,IATOM,MULE,MULTE,
     &            NODV,NOPV,NOCONT,
     &            THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,NUMDIS,
     &            MAXDIS,SQ12EL,RETUR,NHTYP,IREPDM,IFCTYP,ISHLA,ADISTR,
     &            ITYPE,JSTRSH,NPRIMS,NCONTS,IORBSH,
     &            I2TYP,ICEDIF,IFTHRS,DINTSKP,GABRAO,DMRAO,DMRSO,
     &            RELCAL)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL NOPV, NODV, PERTUR, EXPECT, UNDIFF, DDFOCK, DIRFCK, RETUR,
     &        NOCONT, LONDON, SPNORB, DIA2SO, ZFS2EL, DISTRI, SQ12EL,
     &        SUSCEP, ADISTR, RELCAL
      DIMENSION DMAT(*), FMAT(*), HESSEE(*),GMAT(*), IREPDM(*),
     &          IFCTYP(*), INDXAB(*), JSTRSH(*), NPRIMS(*), NCONTS(*),
     &          IORBSH(*),GABRAO(*),DMRAO(*), DMRSO(*),DINTSKP(*),
     &          WORK(LWORK)
#include "hertop.h"
#include "infpar.h"
C NECgh
#include "nuclei.h"
#include "veclen.h"
#include "drw2el.h"
#include "r12int.h"
C NECgh
C
      IF (JPRINT .GT. 5) CALL TITLER('Output from TWOCAL','*',103)
C
      JTOP  = 4*(NHTYP - 1) + MAXDER
C
C     For [r12,T1+T2] integrals, JTOP must be increased by 1. For
C     r12 integrals, JTOP must be 1 for (ss|ss) integrals, but is
C     not increased otherwise (WK/UniKA/14-11-2002).
C
      IF (U12INT) JTOP = JTOP + 1
C     JTOP must be increased for integrals #5 (WK/UniKA/14-11-2002).
      IF (INTGAC .EQ. 5) JTOP = JTOP + 1
      IF (R12INT) JTOP = MAX(JTOP,1)
C
C     For two-electron integrals for relativistic Direct Perturbation
C     Theory (DPT), JTOP must be increased by 2 (WK/UniKA/14-11-2002).
C
      IF (FINDPT .OR. DPTINT) JTOP = JTOP + 2
C
C     For two-electron Breit terms (orbit-orbit and spin-spin) integrals
C     (tuh/wk/07-04-2003)
C
      IF (BPH2OO) JTOP = JTOP + 4
C
      JTOP3 = (JTOP + 1)**3
      NRTOP = (JTOP + 1)*(JTOP + 2)*(JTOP + 3)/6
      KINDHR = 1
      KINDSQ = KINDHR + (  JTOP3 + 1)/IRAT
      KIODHR = KINDSQ + (  NRTOP + 1)/IRAT
#if defined (VAR_VECTOR)
C  NECgh We allocate temporary FOCK matrices for TWOLOP

      KTFMAT = KIODHR + (8*NRTOP + 1)/IRAT
      JUNK   = MAX (IVECLN/NDMAT,1)
      KLAST  = KTFMAT + NDMAT * (NBASIS + NODD)*NBASIS * JUNK
#else
      KLAST  = KIODHR + (8*NRTOP + 1)/IRAT
      KTFMAT = KLAST
#endif
      IF (KLAST .GT. LWORK) CALL STOPIT('TWOCAL',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL HERPRP(WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR))
C
C     *******************************
C     ***** Calculate integrals *****
C     *******************************
C
      IF (DISTRI) THEN
C
C        Distributions
C        =============
C
         CALL TWODIS(WORK(KLAST),LWRK,HESSEE,FMAT,DMAT,NDMAT,GMAT,
     &               INDXAB,MAXDER,EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,
     &               SOFOCK,DISTRI,LONDON,SPNORB,
     &               DIA2SO,ZFS2EL,PERTUR,IATOM,MULE,MULTE,NODV,
     &               NOPV,NOCONT,THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,
     &               IPRNTD,NUMDIS,MAXDIS,SQ12EL,WORK(KINDHR),
     &               WORK(KINDSQ),WORK(KIODHR),IFCTYP,ADISTR,I2TYP,
     &               JSTRSH,NPRIMS,NCONTS,IORBSH,RELCAL)
      ELSE IF (ADISTR) THEN
C
C        All g_abcd for given index a
C        ============================
C
         CALL TWODSA(WORK(KLAST),LWRK,HESSEE,FMAT,DMAT,NDMAT,GMAT,
     &               MAXDER,EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,
     &               SOFOCK,DISTRI,LONDON,SPNORB,
     &               DIA2SO,ZFS2EL,PERTUR,IATOM,MULE,MULTE,NODV,NOPV,
     &               NOCONT,THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,
     &               RETUR,SQ12EL,WORK(KINDHR),WORK(KINDSQ),
     &               WORK(KIODHR),IREPDM,IFCTYP,ISHLA,ADISTR,I2TYP,
     &               JSTRSH,NPRIMS,NCONTS,IORBSH,RELCAL)
C
      ELSE
C
C        All integrals
C        =============
C
#if defined (VAR_MPI) || defined (VAR_PVM)
         IF (SLAVE) THEN
            KTKCPU = KLAST
            KWHICH = KTKCPU +  MTOTTK
            KIJS   = KWHICH + (MTOTTK + 1)/IRAT
            KLAST  = KIJS   + (NTASK  + 1)/IRAT

            IF (KLAST .GT. LWORK)
     &          CALL STOPIT('TWOCAL','PARLOP',KLAST,LWORK)
            LWRK = LWORK - KLAST + 1
C
            CALL PARLOP(WORK(KLAST),LWRK,HESSEE,FMAT,DMAT,NDMAT,
     &           GMAT,MAXDER,EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,
     &           SOFOCK,DISTRI,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &           PERTUR,IATOM,MULE,MULTE,NODV,NOPV,NOCONT,
     &           THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &           SQ12EL,WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR),
     &           WORK(KTKCPU),WORK(KWHICH),WORK(KIJS),
     &           IREPDM,IFCTYP,ADISTR,ITYPE,JSTRSH,NPRIMS,NCONTS,
     &           IORBSH,I2TYP,ICEDIF,IFTHRS,
     &           GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL)
         ELSE
#endif
            CALL TWOLOP(WORK(KLAST),LWRK,HESSEE,FMAT,DMAT,NDMAT,
     &           GMAT,MAXDER,EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,
     &           SOFOCK,DISTRI,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &           PERTUR,IATOM,MULE,MULTE,NODV,NOPV,NOCONT,
     &           THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &           SQ12EL,WORK(KINDHR),WORK(KINDSQ),WORK(KIODHR),
     &           IREPDM,IFCTYP,ADISTR,I2TYP,JSTRSH,NPRIMS,NCONTS,
     &           IORBSH,ICEDIF,IFTHRS,
     &           GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL,WORK(KTFMAT))
#if defined (VAR_MPI) || defined (VAR_PVM)
         END IF
#endif
      END IF
      RETURN
      END
C  /* Deck twolop */
      SUBROUTINE TWOLOP(WORK,LWORK,HESSEE,FMAT,DMAT,NDMAT,GMAT,MAXDER,
     &                  EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,
     &                  SOFOCK,DISTRI,LONDON,SPNORB,
     &                  DIA2SO,ZFS2EL,PERTUR,IATOM,MULE,MULTE,
     &                  NODV,NOPV,NOCONT,THRESH,JPRINT,IPRNTA,IPRNTB,
     &                  IPRNTC,IPRNTD,RETUR,SQ12EL,INDHER,INDHSQ,IODDHR,
     &                  IREPDM,IFCTYP,ADISTR,I2TYP,JSTRSH,NPRIMS,NCONTS,
     &                  IORBSH,ICEDIF,IFTHRS,
     &                  GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL,TFMAT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "dummy.h"
      LOGICAL PRINTA, PRINTB, PRINTC, PRINTD, NOPV, NODV, PERTUR,
     &        EXPECT, UNDIFF, DDFOCK, DIRFCK, DISTRI, NOCONT, SPNORB,
     &        RETUR, FIRST, SQ12EL, LONDON, SUSCEP, ADISTR,
     &        SOFOCK,RELCAL, DIA2SO, ZFS2EL
      DIMENSION DMAT(*), GMAT(*), HESSEE(*), IREPDM(*),
     &        IFCTYP(*), INDHSQ(*), IODDHR(*), INDHER(*),
     &        IORBSH(MXSHEL,MXAOVC), JSTRSH(*), NPRIMS(*), NCONTS(*),
     &        GABRAO(*), DMRAO(*),DMRSO(*),DINTSKP(*),
     &        WORK(LWORK)
#include "twocom.h"
#include "nuclei.h"
#include "energy.h"
#include "taymol.h"
#include "taysol.h"
#include "suscpt.h"
#include "inforb.h"
#include "blocks.h"
#include "symmet.h"
#include "veclen.h"
#include "inftap.h"
#include "gtensor.h"
#include "r12int.h"
#include "zfs.h"
Cef begin
#include "incore.h"
Cef end
      DIMENSION FMAT(NBASIS,NBASIS,*)
#if defined (VAR_VECTOR)
C NECgh Dimensioning TFMAT adapted to always odd first dimension
      DIMENSION TFMAT(NBASIS + NODD,NBASIS,NDMAT,*)
C
C NECgh We zero out the temporary FOCK matrices
C
      JUNK = MAX (IVECLN/NDMAT,1)
      CALL DZERO(TFMAT,(NBASIS + NODD)*NBASIS*NDMAT*JUNK)
#endif
C
      CALL QENTER('TWOLOP')
C
C
      IF (JPRINT .GT. 5) CALL TITLER('Output from TWOLOP','*',103)
C
      FIRST = .TRUE.
      DIRAC = RELCAL
      IF (EXPECT .AND. .NOT.NOPV) THEN
         CALL REWSPL(LUPAO)
      END IF
      IF (SUSCEP) THEN
         IF (.NOT.NOPV) CALL REWSPL(LUPAO)
         IF (.NOT.NOPV) CALL REWSPL(LUPAS)
COBL/AMT
C for CAMB3LYP, we must not zero the susceptibilities the second run in TWOLOP
         IF (.NOT. DONETWOLOP)THEN
             DONETWOLOP = .TRUE.
             CALL DZERO(SUS2EL,9)
         END IF
      END IF
      IF (DIA2SO) THEN
         IF (.NOT. NOPV) CALL REWSPL(LUPAO)
         CALL DZERO(GC2,9)
      END IF
      IF (ZFS2EL) THEN
         IF (.NOT. NOPV) CALL REWSPL(LUPAO)
         CALL DZERO(ZFS,9)
      END IF
C
C     For direct contributions: some parameters
C       FCKTHR is threshold for screening
C       ICEFLG gives information about separate screening of
C              Coulomb/exchange for each DMAT
C       NCM    is the number of DMAT requiring Coulomb-contributions
C       NEM    is the number og DMAT requiring exchange-contributions
C       Screening proceeds in three steps as documented by DINTSKP:
C         Step 1: Screening on integral batches
C           DINTSKP(1,1) - total number of integrals
C           DINTSKP(2,1) - number of integrals skipped (batchwise)
C         Step 2: Screening on individual integrals
C                 while unpacking indices
C           DINTSKP(1,2) - number of integrals remaining after step 1
C           DINTSKP(2,2) - number of integrals skipped
C         Step 3a: Screening on Coulomb contributions
C           DINTSKP(1,3) - NCM times number of integrals remaining
C                         after step 2
C           DINTSKP(2,3) - NCM times number of integrals skipped
C         Step 3b: Screening on exchange contributions
C           DINTSKP(1,4) - NEM times number of integrals remaining
C                         after step 2
C           DINTSKP(2,4) - NEM times number of integrals skipped
C
      DOSCRN = .FALSE.
      IF(DIRFCK.OR.SOFOCK)  THEN
        CALL DZERO(DINTSKP,8)
        IF(IFTHRS.LT.16) THEN
          DOSCRN = .TRUE.
          FCKTHR = (10.0D0)**DBLE(-IFTHRS)
          ICEFLG = ICEDIF
          NCM = 0
          NEM = 0
          DO I = 1,NDMAT
            IY  = MOD(IFCTYP(I),10)
            IC  = MOD(IY,2)
            NCM = NCM + IC
            IE  = (IY - IC)/2
            NEM = NEM + IE
          ENDDO
        ENDIF
      ENDIF
C
      IF(I2TYP.EQ.0) THEN
        IASTRT = 1
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        IASTRT = 1
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = NLRGBL
        IBSMAX = NLRGBL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.3) THEN
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSE
        WRITE(LUPRI,'(A,I5)') 'TWOLOP: Unknown I2TYP =' ,I2TYP
        CALL QUIT('Unknown I2TYP !!!')
      ENDIF
C
Cef begin
#if defined(VAR_MPI)
      LINTSV = .FALSE.
#else
      LINTSV = LINTMP
#endif
Cef end
      I_SHL = 1
      IPREVA = -1
      IPREVB = 0
      IPREVC = 0
      IPREVD = 0
      
C
C     *****************************
C     ***** First Shell Index *****
C     *****************************
C
      DO 100 ISHELA = IASTRT,IASMAX
C
C        Check for basis-set identifier (WK/UniKA/04-11-2002).
         MBSA = MBSISH(ISHELA)
         IF (ISHELA .NE. IASMAX) THEN
            IF (.NOT. R12TRA .AND. MBSISH(ISHELA) .GT. 1) GOTO 100
            IF (MBSA .GT. MBSMAX) GOTO 100
         END IF
         ICA    = LCLASH(ISHELA)
C
         NHKTA  = NHKTSH(ISHELA)
         KHKTA  = KHKTSH(ISHELA)
         KCKTA  = KCKTSH(ISHELA)
         SPHRA  = SPHRSH(ISHELA)
         NCENTA = NCNTSH(ISHELA)
         MULA   = ISTBSH(ISHELA)
         MULTA  = MULT(MULA)
         NSTRA  = IORBSB(IORBSH(ISHELA,1))
         NUCA   = NUCOSH(ISHELA)
         NORBA  = NORBSH(ISHELA)
         IF (.NOT.BIGVEC) THEN
            CORAX0 = CENTSH(ISHELA,1)
            CORAY0 = CENTSH(ISHELA,2)
            CORAZ0 = CENTSH(ISHELA,3)
         END IF
         PRINTA = .TRUE.
         IF ((ISHELA .NE. IPRNTA).AND.(IPRNTA .NE. 0)) PRINTA = .FALSE.
C
C        ******************************
C        ***** Second Shell Index *****
C        ******************************
C
         DO 200 ISHELB = IBSTRT,ISHELA
C
C           Check for basis-set identifier (WK/UniKA/04-11-2002).
            MBSAB = MBSA + MBSISH(ISHELB)
            IF (ISHELB .NE. IBSMAX) THEN
               IF (MBSAB .GT. MBSMAX) GOTO 200
               IF (.NOT. R12TRA .AND. MBSISH(ISHELB) .GT. 1) GOTO 200
            END IF
            ICB    = LCLASH(ISHELB)
C
            NHKTB  = NHKTSH(ISHELB)
            KHKTB  = KHKTSH(ISHELB)
            KCKTB  = KCKTSH(ISHELB)
            SPHRB  = SPHRSH(ISHELB)
            NCENTB = NCNTSH(ISHELB)
            MULB   = ISTBSH(ISHELB)
            MULTB  = MULT(MULB)
            NSTRB  = IORBSB(IORBSH(ISHELB,1))
            NUCB   = NUCOSH(ISHELB)
            NORBB  = NORBSH(ISHELB)
            IF (.NOT.BIGVEC) THEN
               CORBX0 = CENTSH(ISHELB,1)
               CORBY0 = CENTSH(ISHELB,2)
               CORBZ0 = CENTSH(ISHELB,3)
            END IF
            GENAB  = .NOT.(SEGMSH(ISHELA) .AND. SEGMSH(ISHELB))
            IGENAB = 1
            IF (.NOT.GENAB) IGENAB = 2
            NSETA  = NSETSH(ISHELA,IGENAB)
            NSETB  = NSETSH(ISHELB,IGENAB)
            PRINTB = PRINTA
            IF ((ISHELB.NE.IPRNTB).AND.(IPRNTB.NE.0)) PRINTB = .FALSE.
            IF (PRINTB) THEN
               IPRINT = JPRINT
            ELSE
               IPRINT = 0
            END IF
C
C           *****************************
C           ***** Third Shell Index *****
C           *****************************
C
            ICMAX = ISHELA
            IF (SPNORB .OR. DIA2SO) ICMAX = MAXSHL
            IF(I2TYP.EQ.2) ICMAX = NLRGBL
            DO 300 ISHELC = ICSTRT, ICMAX
C
C              Check for basis-set identifier (WK/UniKA/04-11-2002).
               MBSABC = MBSAB + MBSISH(ISHELC)
               IF (ISHELC .NE. ICSMAX) THEN
                  IF (.NOT. R12TRA .AND. MBSISH(ISHELC) .GT. 1) GOTO 300
                  IF (MBSABC .GT. MBSMAX) GOTO 300
               END IF
               ICC    = LCLASH(ISHELC)
C
               NHKTC  = NHKTSH(ISHELC)
               KHKTC  = KHKTSH(ISHELC)
               KCKTC  = KCKTSH(ISHELC)
               SPHRC  = SPHRSH(ISHELC)
               NCENTC = NCNTSH(ISHELC)
               MULC   = ISTBSH(ISHELC)
               MULTC  = MULT(MULC)
               NSTRC  = IORBSB(IORBSH(ISHELC,1))
               NUCC   = NUCOSH(ISHELC)
               NORBC  = NORBSH(ISHELC)
               IF (.NOT.BIGVEC) THEN
                  CORCX0 = CENTSH(ISHELC,1)
                  CORCY0 = CENTSH(ISHELC,2)
                  CORCZ0 = CENTSH(ISHELC,3)
               END IF
               PRINTC = PRINTB
               IF ((ISHELC.NE.IPRNTC).AND.(IPRNTC.NE.0)) PRINTC=.FALSE.
C
C              ******************************
C              ***** Fourth Shell Index *****
C              ******************************
C
               IDMAX = ISHELC
               IF (.NOT.(SPNORB .OR. DIA2SO) .AND.
     &                  (ISHELA.EQ.ISHELC)) IDMAX = ISHELB
               DO 400 ISHELD = IDSTRT,IDMAX
               
C
C                 Check for basis-set identifier (WK/UniKA/04-11-2002).
                  MBSABCD = MBSABC + MBSISH(ISHELD)
                  IF (ISHELD .NE. IDSMAX) THEN
                     IF (.NOT.R12TRA.AND.MBSISH(ISHELD) .GT. 1) GOTO 400
                     IF (MBSABCD .GT. MBSMAX) GOTO 400
                  END IF
                  ICD    = LCLASH(ISHELD)
C
                  NHKTD  = NHKTSH(ISHELD)
                  KHKTD  = KHKTSH(ISHELD)
                  KCKTD  = KCKTSH(ISHELD)
                  SPHRD  = SPHRSH(ISHELD)
                  NCENTD = NCNTSH(ISHELD)
                  MULD   = ISTBSH(ISHELD)
                  MULTD  = MULT(MULD)
                  NSTRD  = IORBSB(IORBSH(ISHELD,1))
                  NUCD   = NUCOSH(ISHELD)
                  NORBD  = NORBSH(ISHELD)
                  IF (.NOT.BIGVEC) THEN
                     CORDX0 = CENTSH(ISHELD,1)
                     CORDY0 = CENTSH(ISHELD,2)
                     CORDZ0 = CENTSH(ISHELD,3)
                  END IF
                  GENCD = .NOT.(SEGMSH(ISHELC) .AND. SEGMSH(ISHELD))
                  IGENCD = 1
                  IF (.NOT.GENCD) IGENCD = 2
                  NSETC = NSETSH(ISHELC,IGENCD)
                  NSETD = NSETSH(ISHELD,IGENCD)
                  PRINTD = PRINTC
                  IF ((ISHELD .NE. IPRNTD).AND.(IPRNTD .NE. 0))
     &               PRINTD = .FALSE.
                  IF (PRINTD) THEN
                     IPRINT = JPRINT
                  ELSE
                     IPRINT = 0
                  END IF
C
                  SHAEQB = ISHELA .EQ. ISHELB
                  SHCEQD = ISHELC .EQ. ISHELD
                  SHABAB = (ISHELA.EQ.ISHELC) .AND. (ISHELB.EQ.ISHELD)
C
C                 *******************************
C                 ***** Calculate integrals *****
C                 *******************************
#if defined (VAR_VECTOR)
C
C For the moment, only INTFCK is vectorized. Use TFMAT only then.
C
      IF (DIRFCK .AND. .NOT.(EXPECT.OR.DDFOCK.OR.SUSCEP.OR.LONDON.OR.
     &   PERTUR.OR.SPNORB.OR.UNDIFF.OR.DISTRI.OR.SOFOCK.OR.DIA2SO.OR.
     &   ZFS2EL)) THEN

                  CALL TWOODS(TFMAT,DMAT,NDMAT,GMAT,HESSEE,WORK,LWORK,
     &                        UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                        EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,
     &                        IATOM,MULE,MULTE,
     &                        MAXDER,NOCONT,NODV,NOPV,THRESH,IPRINT,
     &                        FIRST,SQ12EL,INDHSQ,IODDHR,INDHER,IFCTYP,
     &                        ADISTR,JSTRSH,NPRIMS,NCONTS,IORBSH,DUMMY,
     &                        ICEDIF,IFTHRS,DINTSKP,GABRAO,DMRAO,
     &                        DMRSO,IREPDM)

       ELSE
#endif
                  CALL TWOODS(FMAT,DMAT,NDMAT,GMAT,HESSEE,WORK,LWORK,
     &                        UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                        EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,
     &                        IATOM,MULE,MULTE,
     &                        MAXDER,NOCONT,NODV,NOPV,THRESH,IPRINT,
     &                        FIRST,SQ12EL,INDHSQ,IODDHR,INDHER,IFCTYP,
     &                        ADISTR,JSTRSH,NPRIMS,NCONTS,IORBSH,DUMMY,
     &                        ICEDIF,IFTHRS,DINTSKP,GABRAO,DMRAO,
     &                        DMRSO,IREPDM)
#if defined (VAR_VECTOR)
       END IF
#endif
C
                  IF (RETUR) THEN
                     IF (ISHELA .EQ. IPRNTA .AND.
     &                   ISHELB .EQ. IPRNTB .AND.
     &                   ISHELC .EQ. IPRNTC .AND.
     &                   ISHELD .EQ. IPRNTD) GOTO 9999
                  END IF
  400          CONTINUE
  300       CONTINUE
  200    CONTINUE
  100 CONTINUE
      IF (DIA2SO) THEN
         IF (JPRINT .GT. 1) THEN
            CALL HEADER('Two-electron diamagnetic spin-orbit '//
     &                  'contribution',-1)
            CALL OUTPUT(GC2,1,3,1,3,3,3,1,LUPRI)
         END IF
         IF (.NOT.NOPV) CALL GPCLOSE(LUPAO,'DELETE')
      END IF
      IF (ZFS2EL) THEN
         IF (JPRINT .GT. 1) THEN
            CALL HEADER('Two-electron spin-spin interaction',-1)
            CALL OUTPUT(ZFS,1,3,1,3,3,3,1,LUPRI)
         END IF
ckr         IF (.NOT.NOPV) CALL GPCLOSE(LUPAO,'DELETE')
      END IF
      IF (LONDON .AND. MAXDER.EQ.2) THEN
         SUS2EL(2,1) = SUS2EL(1,2)
         SUS2EL(3,1) = SUS2EL(1,3)
         SUS2EL(3,2) = SUS2EL(2,3)
         IF (JPRINT .GT. 1) THEN
            CALL HEADER('Two-electron integral susceptibilities',-1)
            CALL OUTPUT(SUS2EL,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
#if defined (VAR_VECTOR)
C
C NECgh
C        Here we have to sum up the temporary Fockmatrices
C
      IF(DIRFCK.AND..NOT.(EXPECT.OR.DDFOCK.OR.SUSCEP.OR.LONDON.OR.
     &   PERTUR.OR.SPNORB.OR.UNDIFF.OR.DISTRI.OR.SOFOCK.OR.DIA2SO.OR.
     &   ZFS2EL)) THEN
C
         JUNK = MAX (IVECLN/NDMAT,1)
         DO L = 1,JUNK
            DO I = 1,NDMAT
               DO K = 1,NBASIS
                  DO J = 1,NBASIS
                     FMAT(J,K,I)=FMAT(J,K,I)+TFMAT(J,K,I,L)
                  END DO
               END DO
            END DO
         END DO
      END IF
#endif
C
C     Print Fock matrices
C
      IF ((DIRFCK.OR.SOFOCK) .AND. IPRINT.GT.3) THEN
         CALL HEADER('Fock matrix in TWOINT',-1)
         DO 700 I = 1, NDMAT
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(1,1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  700    CONTINUE
      END IF
 9999 CALL QEXIT('TWOLOP')
      RETURN
      END
C  /* Deck twodsa */
      SUBROUTINE TWODSA(WORK,LWORK,HESSEE,FMAT,DMAT,NDMAT,GMAT,MAXDER,
     &                  EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,SOFOCK,
     &                  DISTRI,LONDON,
     &                  SPNORB,DIA2SO,ZFS2EL,PERTUR,IATOM,MULE,MULTE,
     &                  NODV,NOPV,NOCONT,THRESH,JPRINT,IPRNTA,IPRNTB,
     &                  IPRNTC,IPRNTD,RETUR,SQ12EL,INDHER,INDHSQ,IODDHR,
     &                  IREPDM,IFCTYP,ISHLA,ADISTR,I2TYP,JSTRSH,NPRIMS,
     &                  NCONTS,IORBSH,RELCAL)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "iratdef.h"
#include "dummy.h"
C      PARAMETER (LBUF = 600)
      LOGICAL PRINTA, PRINTB, PRINTC, PRINTD, NOPV, NODV, PERTUR,
     &        EXPECT, UNDIFF, DDFOCK, DIRFCK, DISTRI, NOCONT, SPNORB,
     &        RETUR, FIRST, SQ12EL, LONDON, SUSCEP, ADISTR, DOINDX,
     &        SOFOCK,RELCAL, DIA2SO, ZFS2EL
      DIMENSION DMAT(*), FMAT(*), GMAT(*), HESSEE(*), IREPDM(*),
     &        IFCTYP(*), INDHSQ(*), IODDHR(*), INDHER(*),
     &        IORBSH(MXSHEL,MXAOVC), JSTRSH(*), NPRIMS(*), NCONTS(*),
     &        WORK(LWORK)
#include "cbisol.h"
#include "twocom.h"
#include "nuclei.h"
#include "energy.h"
#include "suscpt.h"
#include "inforb.h"
#include "blocks.h"
#include "symmet.h"
#include "disbuf.h"
#include "eribuf.h"
#include "r12int.h"
#include "drw2el.h"
C
Cholesky
#include <ccdeco.h>
#include <choles.h>
#include <choskp.h>
#include <ccdeco2.h>
Cholesky
C
      CALL QENTER('TWODSA')
      IF (JPRINT .GT. 5) CALL TITLER('Output from TWODSA','*',103)
C
      FIRST  = .TRUE.
      DIRAC = RELCAL
      IF(I2TYP.EQ.0) THEN
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = ISHLA
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = ISHLA
        IBSMAX = NLRGBL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        IBSTRT = NLRGBL+1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = ISHLA
        IBSMAX = MAXSHL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.3) THEN
        IBSTRT = NLRGBL+1
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        IASMAX = ISHLA
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
C
Cholesky
C
      ELSEIF(I2TYP.EQ.4) THEN
        IBSTRT = ISHLB
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = ISHLA
        IBSMAX = ISHLB
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.5) THEN
        IBSTRT = ISHLB
        ICSTRT = ISHLA
        IDSTRT = ISHLB
        IASMAX = ISHLA
        IBSMAX = ISHLB
        ICSMAX = ISHLA
        IDSMAX = ISHLB
C
Cholesky
C
      ELSE
        WRITE(LUPRI,'(A,I5)') 'TWODSA: Unknown I2TYP =' ,I2TYP
      ENDIF
C
C     Allocation of memory for integral buffers.
C     Addresses are save in COMMON /DISBUF/
C
      DOINDX = .FALSE.
      CALL AINDEX(ISHLA,NDIST,IDUMMY,DOINDX,IORBSH,JPRINT)
      LDSBUF = LBUF
      KDSBF  = 1
      IF (U21INT) THEN
        KDUBF = KDSBF + NDIST*LDSBUF
      ELSE
        KDUBF = KDSBF
      END IF
      KDSIBF = KDUBF + NDIST*LDSBUF
      KDSNCT = KDSIBF + (2*NDIST*LDSBUF - 1)/IRAT + 1
      KDSORB = KDSNCT + (NDIST          - 1)/IRAT + 1
      KORBDS = KDSORB + (NDIST          - 1)/IRAT + 1
      KDSEND = KORBDS + (NBASIS         - 1)/IRAT + 1
C
      KDSBUF = 1
      KLAST = KDSBUF + KDSEND
      LWRK  = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('TWODSA',' ',KLAST,LWORK)
C
C     *****************************
C     ***** First Shell Index *****
C     *****************************
C
      ISHELA = ISHLA
C
      MBSA   = MBSISH(ISHELA)
      NHKTA  = NHKTSH(ISHELA)
      KHKTA  = KHKTSH(ISHELA)
      KCKTA  = KCKTSH(ISHELA)
      SPHRA  = SPHRSH(ISHELA)
      NCENTA = NCNTSH(ISHELA)
      MULA   = ISTBSH(ISHELA)
      MULTA  = MULT(MULA)
      NSTRA  = IORBSB(IORBSH(ISHELA,1))
      NUCA   = NUCOSH(ISHELA)
      NORBA  = NORBSH(ISHELA)
      IF (.NOT.BIGVEC) THEN
         CORAX0 = CENTSH(ISHELA,1)
         CORAY0 = CENTSH(ISHELA,2)
         CORAZ0 = CENTSH(ISHELA,3)
      END IF
      PRINTA = .TRUE.
      IF ((ISHELA .NE. IPRNTA).AND.(IPRNTA .NE. 0)) PRINTA = .FALSE.
C
C     ******************************
C     ***** Second Shell Index *****
C     ******************************
C
Cholesky
      IALDON = 0
      IALSKP = 0
      J2TYP  = I2TYP
Cholesky
C
      DO 200 ISHELB = IBSTRT,IBSMAX
C        Check for basis-set identifier (WK/UniKA/04-11-2002).
         MBSAB = MBSA + MBSISH(ISHELB)
         IF (ISHELB .NE. IBSMAX) THEN
            IF (.NOT. (R12HYB .AND. COMBSS)) THEN
            IF ((.NOT. R12TRA .OR. ONEAUX) 
     &           .AND. MBSISH(ISHELB) .GT. 1) GOTO 200
            END IF
            IF (MBSAB .GT. MBSMAX) GOTO 200
         END IF
         NHKTB  = NHKTSH(ISHELB)
         KHKTB  = KHKTSH(ISHELB)
         KCKTB  = KCKTSH(ISHELB)
         SPHRB  = SPHRSH(ISHELB)
         NCENTB = NCNTSH(ISHELB)
         MULB   = ISTBSH(ISHELB)
         MULTB  = MULT(MULB)
         NSTRB  = IORBSB(IORBSH(ISHELB,1))
         NUCB   = NUCOSH(ISHELB)
         NORBB  = NORBSH(ISHELB)
         IF (.NOT.BIGVEC) THEN
            CORBX0 = CENTSH(ISHELB,1)
            CORBY0 = CENTSH(ISHELB,2)
            CORBZ0 = CENTSH(ISHELB,3)
         END IF
         GENAB  = .NOT.(SEGMSH(ISHELA) .AND. SEGMSH(ISHELB))
         IGENAB = 1
         IF (.NOT.GENAB) IGENAB = 2
         NSETA  = NSETSH(ISHELA,IGENAB)
         NSETB  = NSETSH(ISHELB,IGENAB)
         PRINTB = PRINTA
         IF ((ISHELB.NE.IPRNTB).AND.(IPRNTB.NE.0)) PRINTB = .FALSE.
         IF (PRINTB) THEN
            IPRINT = JPRINT
         ELSE
            IPRINT = 0
         END IF
C
C        *****************************
C        ***** Third Shell Index *****
C        *****************************
C
         DO 300 ISHELC = ICSTRT,ICSMAX
C           Check for basis-set identifier (WK/UniKA/04-11-2002).
            MBSABC = MBSAB + MBSISH(ISHELC)
            IF (ISHELB .NE. IBSMAX .OR. ISHELC .NE. ICSMAX) THEN
               IF (.NOT. (R12HYB .AND. COMBSS)) THEN
               IF ((.NOT. R12TRA .OR. ONEAUX)
     &             .AND. MBSISH(ISHELC) .GT. 1) GOTO 300
               END IF
               IF (MBSABC .GT. MBSMAX) GOTO 300
            END IF
            NHKTC  = NHKTSH(ISHELC)
            KHKTC  = KHKTSH(ISHELC)
            KCKTC  = KCKTSH(ISHELC)
            SPHRC  = SPHRSH(ISHELC)
            NCENTC = NCNTSH(ISHELC)
            MULC   = ISTBSH(ISHELC)
            MULTC  = MULT(MULC)
            NSTRC  = IORBSB(IORBSH(ISHELC,1))
            NUCC   = NUCOSH(ISHELC)
            NORBC  = NORBSH(ISHELC)
            IF (.NOT.BIGVEC) THEN
               CORCX0 = CENTSH(ISHELC,1)
               CORCY0 = CENTSH(ISHELC,2)
               CORCZ0 = CENTSH(ISHELC,3)
            END IF
            PRINTC = PRINTB
            IF ((ISHELC.NE.IPRNTC).AND.(IPRNTC.NE.0)) PRINTC=.FALSE.
C
C           ******************************
C           ***** Fourth Shell Index *****
C           ******************************
C
            ICMAX = ISHELC
            IF (SPNORB .OR. DIA2SO) ICMAX = MAXSHL
Cholesky
            IF (I2TYP.EQ.5) ICMAX=ISHELB
Cholesky
            DO 400 ISHELD = IDSTRT, ICMAX
C              Check for basis-set identifier (WK/UniKA/04-11-2002).
               MBSABCD = MBSABC + MBSISH(ISHELD)
               IF (ISHELB .NE. IBSMAX .OR.
     &             ISHELC .NE. ICSMAX .OR.
     &             ISHELD .NE. ICMAX) THEN
                  IF (.NOT. (R12HYB .AND. COMBSS)) THEN
                  IF ((.NOT. R12TRA .OR. ONEAUX) 
     &                 .AND. MBSISH(ISHELD) .GT. 1) GOTO 400
                  END IF
                  IF (MBSABCD .GT. MBSMAX) GOTO 400
               END IF
               NHKTD  = NHKTSH(ISHELD)
               KHKTD  = KHKTSH(ISHELD)
               KCKTD  = KCKTSH(ISHELD)
               SPHRD  = SPHRSH(ISHELD)
               NCENTD = NCNTSH(ISHELD)
               MULD   = ISTBSH(ISHELD)
               MULTD  = MULT(MULD)
               NSTRD  = IORBSB(IORBSH(ISHELD,1))
               NUCD   = NUCOSH(ISHELD)
               NORBD  = NORBSH(ISHELD)
               IF (.NOT.BIGVEC) THEN
                  CORDX0 = CENTSH(ISHELD,1)
                  CORDY0 = CENTSH(ISHELD,2)
                  CORDZ0 = CENTSH(ISHELD,3)
               END IF
               GENCD = .NOT.(SEGMSH(ISHELC) .AND. SEGMSH(ISHELD))
               IGENCD = 1
               IF (.NOT.GENCD) IGENCD = 2
               NSETC = NSETSH(ISHELC,IGENCD)
               NSETD = NSETSH(ISHELD,IGENCD)
               PRINTD = PRINTC
               IF ((ISHELD .NE. IPRNTD).AND.(IPRNTD .NE. 0))
     &            PRINTD = .FALSE.
               IF (PRINTD) THEN
                  IPRINT = JPRINT
               ELSE
                  IPRINT = 0
               END IF
C
               SHAEQB = .FALSE.
               IF (U21INT .OR. BPH2OO) THEN
                  SHCEQD = .FALSE.
               ELSE
                  SHCEQD = ISHELC .EQ. ISHELD
               END IF
               SHABAB = .FALSE.
C
C_to_do        Cholesky screening.
C              -------------------
C
c               IF (I2TYP .EQ. 4) THEN
C
c                 IF (DIASC1(ISHELA,ISHELB) .EQ. ZERO) THEN
c                    write(LUPRI,*) 'Ishelab screening:',ishela,ishelb
c                    IALSKP = IALSKP + 1
c                    GOTO 400
c                 ENDIF
c
c                 IF (DIASC1(ISHELC,ISHELD) .EQ. ZERO) THEN
c                    write(LUPRI,*) 'Ishelcd screening:',ishelc,isheld
c                    write(LUPRI,*) 'Other indices:',ishela,ishelb
c                    IALSKP = IALSKP + 1
c                    GOTO 400
c                 ENDIF
C
c                 write(LUPRI,*) 'Distribution calculated',ishela,ishelb,
c    &                        ishelc,isheld
c                 IALDON = IALDON + 1
C
c              ENDIF
C
C_to_do_end
C
C              *******************************
C              ***** Calculate integrals *****
C              *******************************
C
               CALL TWOODS(FMAT,DMAT,NDMAT,GMAT,HESSEE,WORK(KLAST),LWRK,
     &                     UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                     EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,
     &                     IATOM,MULE,MULTE,MAXDER,
     &                     NOCONT,NODV,NOPV,THRESH,IPRINT,FIRST,
     &                     SQ12EL,INDHSQ,IODDHR,INDHER,IFCTYP,ADISTR,
     &                     JSTRSH,NPRIMS,NCONTS,IORBSH,WORK(KDSBUF),
     &                     IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                     IDUMMY)
C
               IF (RETUR) THEN
                  IF (ISHELA .EQ. IPRNTA .AND.
     &                ISHELB .EQ. IPRNTB .AND.
     &                ISHELC .EQ. IPRNTC .AND.
     &                ISHELD .EQ. IPRNTD) GOTO 9999
               END IF
  400       CONTINUE
  300    CONTINUE
  200 CONTINUE
 9999 CALL QEXIT('TWODSA')
      RETURN
      END
C  /* Deck twodis */
      SUBROUTINE TWODIS(WORK,LWORK,HESSEE,FMAT,DMAT,NDMAT,GMAT,INDXAB,
     &                  MAXDER,EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,
     &                  SOFOCK,DISTRI,
     &                  LONDON,SPNORB,DIA2SO,ZFS2EL,PERTUR,IATOM,MULE,
     &                  MULTE,NODV,NOPV,NOCONT,THRESH,JPRINT,IPRNTA,
     &                  IPRNTB,IPRNTC,IPRNTD,NUMDIS,MAXDIS,SQ12EL,
     &                  INDHER,INDHSQ,IODDHR,IFCTYP,ADISTR,I2TYP,JSTRSH,
     &                  NPRIMS,NCONTS,IORBSH,RELCAL)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      LOGICAL NOPV, NODV, PERTUR,
     &        EXPECT, UNDIFF, DDFOCK, DIRFCK, NOCONT, SPNORB, DIA2SO,
     &        FIRST, DISTRI, SQ12EL, LONDON, SUSCEP, ADISTR, SOFOCK,
     &        RELCAL, ZFS2EL
      DIMENSION DMAT(*), FMAT(*), GMAT(*), HESSEE(*), IFCTYP(*),
     &        INDXAB(MAXDIS), INDHSQ(*), IODDHR(*),
     &        IORBSH(MXSHEL,MXAOVC), JSTRSH(*), NPRIMS(*), NCONTS(*),
     &        WORK(LWORK)
#include "blocks.h"
#include "twocom.h"
C
Cholesky #include "inforb.h"
#include "ccorb.h"
Cholesky
C
#include "symmet.h"
#include "r12int.h"
      SAVE ISHLAB
C
      CALL QENTER('TWODIS')
      IF (JPRINT .GT. 5) CALL TITLER('Output from TWODIS','*',103)
C
      IF (NUMDIS .EQ. -1) THEN
         FIRST = .TRUE.
         ISHLAB = 1
      ELSE
         FIRST = .FALSE.
      END IF
C
      IF(I2TYP.EQ.0) THEN
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        IASMAX = NLRGBL
        IBSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.3) THEN
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
      ELSE
        WRITE(LUPRI,'(A,I5)') 'TWODIS: Unknown I2TYP =' ,I2TYP
      ENDIF
C
      IF (ISHLAB .GT. MAXSHL*(MAXSHL + 1)/2) THEN
         NUMDIS = -1
         IF (JPRINT .GT. 1) THEN
            WRITE (LUPRI,'(/A/)')
     &         '  All integral distributions have now been calculated.'
         END IF
      ELSE
         CALL UNPKIJ(ISHLAB,ISHELA,ISHELB)
C
C        First Shell Index
C        =================
C
         ICA    = LCLASH(ISHELA)
C
         NHKTA  = NHKTSH(ISHELA)
         KHKTA  = KHKTSH(ISHELA)
         KCKTA  = KCKTSH(ISHELA)
         SPHRA  = SPHRSH(ISHELA)
         NCENTA = NCNTSH(ISHELA)
         MULA   = ISTBSH(ISHELA)
         MULTA  = MULT(MULA)
         NSTRA  = IORBSB(IORBSH(ISHELA,1))
         NUCA   = NUCOSH(ISHELA)
         NORBA  = NORBSH(ISHELA)
         IF (.NOT.BIGVEC) THEN
            CORAX0 = CENTSH(ISHELA,1)
            CORAY0 = CENTSH(ISHELA,2)
            CORAZ0 = CENTSH(ISHELA,3)
         END IF
C
C        Second Shell Index
C        ==================
C
         ICB    = LCLASH(ISHELB)
C
         NHKTB  = NHKTSH(ISHELB)
         KHKTB  = KHKTSH(ISHELB)
         KCKTB  = KCKTSH(ISHELB)
         SPHRB  = SPHRSH(ISHELB)
         NCENTB = NCNTSH(ISHELB)
         MULB   = ISTBSH(ISHELB)
         MULTB  = MULT(MULB)
         NSTRB  = IORBSB(IORBSH(ISHELB,1))
         NUCB   = NUCOSH(ISHELB)
         NORBB  = NORBSH(ISHELB)
         IF (.NOT.BIGVEC) THEN
            CORBX0 = CENTSH(ISHELB,1)
            CORBY0 = CENTSH(ISHELB,2)
            CORBZ0 = CENTSH(ISHELB,3)
         END IF
         GENAB  = .NOT.(SEGMSH(ISHELA) .AND. SEGMSH(ISHELB))
         IGENAB = 1
         IF (.NOT.GENAB) IGENAB = 2
         NSETA  = NSETSH(ISHELA,IGENAB)
         NSETB  = NSETSH(ISHELB,IGENAB)
         SHAEQB = ISHELA .EQ. ISHELB
         DIAGAB = SHAEQB .AND. .NOT.BIGVEC
C
         NUMDIS = NDISTR(INDXAB,MAXDIS,JPRINT)
         CALL DZERO(GMAT,NBAST*NBAST*NUMDIS)
         IF (JPRINT .GT. 1) THEN
            WRITE (LUPRI,'(//,2(1X,A,I5,/))')
     &         ' Calculation of integral distribution set:',ISHLAB,
     &         ' Number of distributions in this set:     ',NUMDIS
            IF (JPRINT .GT. 4) THEN
               CALL HEADER('Pair index   Orbital A    Orbital B',2)
               DO 100 I = 1, NUMDIS
                  CALL UNPKIJ(INDXAB(I),IA,IB)
                  WRITE (LUPRI,'(1X,3(I8,5X))')  INDXAB(I), IA, IB
 100           CONTINUE
            END IF
         END IF
         CALL DISLOP(GMAT,FMAT,DMAT,NDMAT,HESSEE,WORK,LWORK,MAXDER,
     &               EXPECT,
     &               SUSCEP,UNDIFF,DDFOCK,DIRFCK,SOFOCK,DISTRI,LONDON,
     &               SPNORB,DIA2SO,ZFS2EL,
     &               PERTUR,IATOM,MULE,MULTE,NODV,NOPV,NOCONT,THRESH,
     &               JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,SQ12EL,INDHSQ,
     &               IODDHR,INDHER,IFCTYP,ADISTR,JSTRSH,NPRIMS,NCONTS,
     &               IORBSH,I2TYP,MBSAB)
 999     ISHLAB = ISHLAB + 1
      END IF
C
      CALL QEXIT('TWODIS')
      RETURN
      END
C  /* Deck dislop */
      SUBROUTINE DISLOP(GMAT,FMAT,DMAT,NDMAT,HESSEE,WORK,LWORK,MAXDER,
     &                  EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,SOFOCK,
     &                  DISTRI,LONDON,
     &                  SPNORB,DIA2SO,ZFS2EL,PERTUR,IATOM,MULE,MULTE,
     &                  NODV,NOPV,NOCONT,THRESH,JPRINT,IPRNTA,IPRNTB,
     &                  IPRNTC,IPRNTD,SQ12EL,INDHSQ,IODDHR,INDHER,
     &                  IFCTYP,ADISTR,JSTRSH,NPRIMS,NCONTS,IORBSH,I2TYP,
     &                  MBSAB)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "dummy.h"
      LOGICAL PRINTA, PRINTB, PRINTC, PRINTD, NOPV, NODV, PERTUR,
     &        EXPECT, UNDIFF, DDFOCK, DIRFCK, NOCONT, SPNORB, DIA2SO,
     &        FIRST, SQ12EL, DISTRI, LONDON, SUSCEP, ADISTR, SOFOCK,
     &        ZFS2EL
      DIMENSION GMAT(*), DMAT(*), FMAT(*), HESSEE(*), IFCTYP(*),
     &          INDHSQ(*), IODDHR(*), INDHER(*), IORBSH(MXSHEL,MXAOVC),
     &          JSTRSH(*), NPRIMS(*), NCONTS(*), WORK(LWORK)
#include "twocom.h"
#include "nuclei.h"
#include "energy.h"
#include "taymol.h"
#include "inforb.h"
#include "blocks.h"
#include "symmet.h"
#include "r12int.h"
C
      IF (JPRINT .GT. 5) CALL TITLER('Output from DISLOP','*',103)
C
      PRINTA = .TRUE.
      IF ((ISHELA .NE. IPRNTA).AND.(IPRNTA .NE. 0)) PRINTA = .FALSE.
      PRINTB = PRINTA
      IF ((ISHELB.NE.IPRNTB).AND.(IPRNTB.NE.0)) PRINTB = .FALSE.
      IF (PRINTB) THEN
         IPRINT = JPRINT
      ELSE
         IPRINT = 0
      END IF
      IF(I2TYP.EQ.0) THEN
        ICSTRT = 1
        IDSTRT = 1
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        ICSTRT = 1
        IDSTRT = 1
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        ICSTRT = 1
        IDSTRT = 1
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.3) THEN
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSE
        WRITE(LUPRI,'(A,I5)') 'DISLOP: Unknown I2TYP =' ,I2TYP
      ENDIF
C
C     *****************************
C     ***** Third Shell Index *****
C     *****************************
C
      DO 100 ISHELC = ICSTRT, ICSMAX
         NHKTC  = NHKTSH(ISHELC)
         KHKTC  = KHKTSH(ISHELC)
         KCKTC  = KCKTSH(ISHELC)
         SPHRC  = SPHRSH(ISHELC)
         NCENTC = NCNTSH(ISHELC)
         MULC   = ISTBSH(ISHELC)
         MULTC  = MULT(MULC)
         NSTRC  = IORBSB(IORBSH(ISHELC,1))
         NUCC   = NUCOSH(ISHELC)
         NORBC  = NORBSH(ISHELC)
         IF (.NOT.BIGVEC) THEN
            CORCX0 = CENTSH(ISHELC,1)
            CORCY0 = CENTSH(ISHELC,2)
            CORCZ0 = CENTSH(ISHELC,3)
         END IF
         PRINTC = PRINTB
         IF ((ISHELC.NE.IPRNTC).AND.(IPRNTC.NE.0)) PRINTC=.FALSE.
C
C        ******************************
C        ***** Fourth Shell Index *****
C        ******************************
C
         IDMAX = ISHELC
         DO 200 ISHELD = IDSTRT, ISHELC
            NHKTD  = NHKTSH(ISHELD)
            KHKTD  = KHKTSH(ISHELD)
            KCKTD  = KCKTSH(ISHELD)
            SPHRD  = SPHRSH(ISHELD)
            NCENTD = NCNTSH(ISHELD)
            MULD   = ISTBSH(ISHELD)
            MULTD  = MULT(MULD)
            NSTRD  = IORBSB(IORBSH(ISHELD,1))
            NUCD   = NUCOSH(ISHELD)
            NORBD  = NORBSH(ISHELD)
            IF (.NOT.BIGVEC) THEN
               CORDX0 = CENTSH(ISHELD,1)
               CORDY0 = CENTSH(ISHELD,2)
               CORDZ0 = CENTSH(ISHELD,3)
            END IF
            GENCD = .NOT.(SEGMSH(ISHELC) .AND. SEGMSH(ISHELD))
            IGENCD = 1
            IF (.NOT.GENCD) IGENCD = 2
            NSETC = NSETSH(ISHELC,IGENCD)
            NSETD = NSETSH(ISHELD,IGENCD)
            SHCEQD = ISHELC .EQ. ISHELD
            SHABAB = (ISHELA.EQ.ISHELC) .AND. (ISHELB.EQ.ISHELD)
            PRINTD = PRINTC
            IF ((ISHELD .NE. IPRNTD).AND.(IPRNTD .NE. 0))
     &         PRINTD = .FALSE.
            IF (PRINTD) THEN
               IPRINT = JPRINT
            ELSE
               IPRINT = 0
            END IF
C
C           *******************************
C           ***** Calculate integrals *****
C           *******************************
C
            CALL TWOODS(FMAT,DMAT,NDMAT,GMAT,HESSEE,WORK,LWORK,UNDIFF,
     &                  PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,EXPECT,
     &                  SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,IATOM,MULE,
     &                  MULTE,MAXDER,NOCONT,
     &                  NODV,NOPV,THRESH,IPRINT,FIRST,SQ12EL,INDHSQ,
     &                  IODDHR,INDHER,IFCTYP,ADISTR,JSTRSH,NPRIMS,
     &                  NCONTS,IORBSH,DUMMY,
     &                  IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                  IDUMMY)
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck twoods */
      SUBROUTINE TWOODS(FMAT,DMAT,NDMAT,GMAT,HESSEE,WORK,LWORK,UNDIFF,
     &                  PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,EXPECT,
     &                  SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,IATOM,MULE,
     &                  MULTE,MAXDER,NOCONT,
     &                  NODV,NOPV,THRESH,IPRINT,FIRST,SQ12EL,INDHSQ,
     &                  IODDHR,INDHER,IFCTYP,ADISTR,JSTRSH,NPRIMS,
     &                  NCONTS,IORBSH,WRKBUF,ICEDIF,
     &                  IFTHRS,DINTSKP,GABRAO,DMRAO,DMRSO,IREPDM)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "iratdef.h"
#include "maxaqn.h"
#include "mxcent.h"
      LOGICAL NOPV, NODV, PERTUR, EXPECT, UNDIFF, NOCONT, SPNORB, FIRST,
     &        DDFOCK, DIRFCK, DISTRI, SQ12EL, LONDON, SUSCEP, ADISTR,
     &        SOFOCK, DIA2SO, ZFS2EL
      DIMENSION DMAT(*), FMAT(*), GMAT(*), HESSEE(*), IFCTYP(*),
     &          INDHSQ(*), IODDHR(*), INDHER(*), WRKBUF(*), JSTRSH(*),
     &          NPRIMS(*), NCONTS(*), IORBSH(*), WORK(LWORK),
     &          GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*),IREPDM(*)
#include "inftap.h"
#include "symmet.h"
#include "twocom.h"
#include "blocks.h"
#include "twosta.h"
#include "drw2el.h"
#include "expopt.h"
#if defined (VAR_SPLITFILES)
      CHARACTER*80 FNNAME, FNNM2
#include "chrnos.h"
#endif
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from TWOODS','*',103)
      IF (FINDPT.OR.DPTINT) THEN
         JMAXA  = NHKTA + 1
         JMAXB  = NHKTB + 1
         JMAXC  = NHKTC + 1
         JMAXD  = NHKTD + 1
      ELSE IF (BPH2OO) THEN
         JMAXA = NHKTA
         JMAXB = NHKTB
         JMAXC = NHKTC
         JMAXD = NHKTD
      ELSE IF (DIA2SO) THEN
         JMAXA = NHKTA
         JMAXB = NHKTB - 1
         JMAXC = NHKTC
         JMAXD = NHKTD
      ELSE IF (ZFS2EL) THEN
         JMAXA = NHKTA - 1 + MAXDER
         JMAXB = NHKTB - 1 + MAXDER
         JMAXC = NHKTC - 1
         JMAXD = NHKTD - 1
      ELSE
         JMAXA  = NHKTA - 1 + MAXDER
         JMAXB  = NHKTB - 1 + MAXDER
         JMAXC  = NHKTC - 1 + MAXDER
         JMAXD  = NHKTD - 1 + MAXDER
         IF (LONDON) JMAXB = NHKTB - 1
         IF (LONDON) JMAXD = NHKTD - 1
      END IF
      MAXAB  = NHKTA + NHKTB - 2
      MAXCD  = NHKTC + NHKTD - 2
      TCONAB = SHAEQB .AND. MAXAB .EQ. 0
      TCONCD = SHCEQD .AND. MAXCD .EQ. 0
      DIAGAB = SHAEQB .AND. .NOT.BIGVEC
      DIAGCD = SHCEQD .AND. .NOT.BIGVEC
      SPHRAB = SPHRA .OR. SPHRB
      SPHRCD = SPHRC .OR. SPHRD
      DIACAB = DIAGAB .AND. .NOT.SPHRAB
      DIACCD = DIAGCD .AND. .NOT.SPHRCD
C
C     Number of two-electron densities elements
C
C
      IF ((EXPECT .OR. SUSCEP .OR. DIA2SO .OR. ZFS2EL) 
     &     .AND..NOT .NOPV) THEN
#if defined (VAR_SPLITFILES)
         READ (LUPAO,END=414) NPMAT
         GOTO 415
 414     INQUIRE(UNIT=LUPAO,NAME=FNNAME)
         LN = 1
 410     CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 410
         END IF
         LN = LN - 1
         CALL GPCLOSE(LUPAO,'KEEP')
         I = LN - 1
         IF (FNNAME(I:I) .NE. '-') THEN
            FNNM2 = FNNAME(1:LN)//'-0'
            LN = LN + 2
         ELSE
            READ(FNNAME(LN:),'(I1)') INUM
            INUM = INUM + 1
            IF (INUM .GT. 9) THEN
               WRITE (LUPRI,'(/A)') ' DALTON needs to split a '//
     &              ' file more than 11 times.',
     &              ' This is currently not supported'
               CALL QUIT('Too many splittings of a file')
            END IF
            FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
         END IF
         CALL GPOPEN(LUPAO,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &               .FALSE.)
         IF (SUSCEP) THEN
            INQUIRE(UNIT=LUPAS,NAME=FNNAME)
            LN = 1
 411        CONTINUE
            IF (FNNAME(LN:LN) .NE. ' ') THEN
               LN = LN + 1
               GOTO 411
            END IF
            LN = LN - 1
            CALL GPCLOSE(LUPAS,'KEEP')
            I = LN - 1
            IF (FNNAME(I:I) .NE. '-') THEN
               FNNM2 = FNNAME(1:LN)//'-0'
               LN = LN + 2
            ELSE
               READ(FNNAME(LN:),'(I1)') INUM
               INUM = INUM + 1
               IF (INUM .GT. 9) THEN
                  WRITE (LUPRI,'(/A)') ' DALTON needs to split a '//
     &                 ' file more than 11 times.',
     &                 ' This is currently not supported'
                  CALL QUIT('Too many splittings of a file')
               END IF
               FNNM2 = FNNAME(1:I)//CHRNOS(INUM)
            END IF
            CALL GPOPEN(LUPAS,FNNM2(1:LN),'UNKNOWN',' ',' ',IDUMMY,
     &                  .FALSE.)
         END IF
         READ (LUPAO) NPMAT
 415     CONTINUE
#else
         READ (LUPAO) NPMAT
#endif
         IF (SUSCEP) THEN
            READ (LUPAS) NPMATA
            IF (NPMAT .NE. NPMATA) THEN
               WRITE (LUPRI,'(1X,A,I5,A,I5)')
     &         ' Error in TWOODS: NPMAT = ',NPMAT,', NPMATA = ',NPMATA
               CALL QUIT('Program aborted due to error in TWOODS')
            END IF
         ELSE
            NPMATA = 0
         END IF
      ELSE
         NPMAT  = 0
         NPMATA = 0
      END IF
      NODCAB = NODSYM(MAXOPR,MULA,MULB)
      NODCCD = NODSYM(MAXOPR,MULC,MULD)
C
      KPSO   = 1
      KPSA   = KPSO   + NPMAT
      KNPCOA = KPSA   + NPMATA
      KNPCOB = KNPCOA + (2*NSETA*(NODCAB + 1) + 1)/IRAT
      KNPCOC = KNPCOB + (2*NSETB*(NODCAB + 1) + 1)/IRAT
      KNPCOD = KNPCOC + (2*NSETC*(NODCCD + 1) + 1)/IRAT
      KJSTRA = KNPCOD + (2*NSETD*(NODCCD + 1) + 1)/IRAT
      KJSTRB = KJSTRA + (NSETA + 1)/IRAT
      KJSTRC = KJSTRB + (NSETB + 1)/IRAT
      KJSTRD = KJSTRC + (NSETC + 1)/IRAT
      KLAST  = KJSTRD + (NSETD + 1)/IRAT
C
      IF (EXPGRA) THEN
         LIXPAB = (2*NUCA*NUCB*NODCAB + 1)/IRAT
         LIXPCD = (2*NUCC*NUCD*NODCCD + 1)/IRAT
      ELSE
         LIXPAB = 0
         LIXPCD = 0
      END IF
      KCKMAX = MAX(KCKTA,KCKTB,KCKTC,KCKTD)
      KCORAB = KLAST
      KCORCD = KCORAB + 9*NUCA*NUCB*NODCAB
      KEXPAB = KCORCD + 9*NUCC*NUCD*NODCCD
      KEXPCD = KEXPAB + 3*NUCA*NUCB*NODCAB
      KIXPAB = KEXPCD + 3*NUCC*NUCD*NODCCD
      KIXPCD = KIXPAB +   LIXPAB
      KFACAB = KIXPCD +   LIXPCD
      KFACCD = KFACAB +   NUCA*NUCB*NODCAB
      KLMNVL = KFACCD +   NUCC*NUCD*NODCCD
      KLAST  = KLMNVL + (20*KCKMAX  + 1)/IRAT
      IF (GENAB) THEN
         KCONTA = KLAST
         KCONTB = KCONTA + 2*NORBA*NUCA*NODCAB
         KPNTAB = KCONTB + 2*NORBB*NUCB*NODCAB
         KREDAB = KPNTAB + (2*NUCA*NUCB*NODCAB + 1)/IRAT
         KLAST  = KREDAB + (NORBA*NORBB + 1)/IRAT
         KNCSAB = KLAST
      ELSE
         KCONTA = KLAST
         KCONTB = KLAST
         KPNTAB = KLAST
         KREDAB = KLAST
         KNCSAB = KLAST
         KLAST  = KNCSAB + (NORBA*NORBB*NODCAB + 1)/IRAT
      END IF
      IF (GENCD) THEN
         KCONTC = KLAST
         KCONTD = KCONTC + 2*NORBC*NUCC*NODCCD
         KPNTCD = KCONTD + 2*NORBD*NUCD*NODCCD
         KREDCD = KPNTCD + (2*NUCC*NUCD*NODCCD + 1)/IRAT
         KLAST  = KREDCD + (NORBC*NORBD + 1)/IRAT
         KNCSCD = KLAST
      ELSE
         KCONTC = KLAST
         KCONTD = KLAST
         KPNTCD = KLAST
         KREDCD = KLAST
         KNCSCD = KLAST
         KLAST  = KNCSCD + (NORBC*NORBD*NODCCD + 1)/IRAT
      END IF
      KINDAB = KLAST
C.....New code to fix Heike's error on tchplx4.
C     Add 2 in place of 1 (WK/UniKA/18-06-2004).
C     KINDCD = KINDAB + (2*NORBA*NORBB + 1)/IRAT
C     KLAST  = KINDCD + (2*NORBC*NORBD + 1)/IRAT
C
      KINDCD = KINDAB + (2*NORBA*NORBB + 2)/IRAT
      KLAST  = KINDCD + (2*NORBC*NORBD + 2)/IRAT
C....New code
      LWRK = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('TWOODS',' ',KLAST,LWORK)
      LWTOT  = LWTOT + KLAST
      MWTOT  = MAX(MWTOT,LWTOT)
C
      CALL TWOOD1(FMAT,DMAT,NDMAT,GMAT,HESSEE,WORK(KLAST),LWRK,
     &            WORK(KPSO),WORK(KPSA),NPMAT,WORK(KCORAB),WORK(KCORCD),
     &            WORK(KEXPAB),WORK(KEXPCD),WORK(KFACAB),WORK(KFACCD),
     &            WORK(KCONTA),WORK(KCONTB),WORK(KCONTC),WORK(KCONTD),
     &            UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,EXPECT,
     &            SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,IATOM,MULE,MULTE,
     &            MAXDER,NOCONT,NODV,NOPV,THRESH,IPRINT,FIRST,SQ12EL,
     &            WORK(KNPCOA),WORK(KNPCOB),WORK(KNPCOC),WORK(KNPCOD),
     &            WORK(KNCSAB),WORK(KNCSCD),WORK(KJSTRA),WORK(KJSTRB),
     &            WORK(KJSTRC),WORK(KJSTRD),WORK(KINDAB),WORK(KINDCD),
     &            INDHSQ,IODDHR,INDHER,WORK(KLMNVL),
     &            WORK(KPNTAB),WORK(KPNTCD),WORK(KREDAB),WORK(KREDCD),
     &            IFCTYP,ADISTR,JSTRSH,NPRIMS,NCONTS,IORBSH,WRKBUF,
     &            GABRAO,DMRAO,DMRSO,DINTSKP,IREPDM,
     &            WORK(KIXPAB),WORK(KIXPCD))
C
      LWTOT  = LWTOT - KLAST
      RETURN
      END
C  /* Deck twood1 */      
      SUBROUTINE TWOOD1(FMAT,DMAT,NDMAT,GMAT,HESSEE,WORK,LWORK,PSO,PSA,
     &                  NPMAT,COORAB,COORCD,EXPAB,EXPCD,FACAB,FACCD,
     &                  CONTA,CONTB,CONTC,CONTD,UNDIFF,PERTUR,LONDON,
     &                  SPNORB,DIA2SO,ZFS2EL,EXPECT,SUSCEP,DDFOCK,
     &                  DIRFCK,SOFOCK,DISTRI,IATOM,
     &                  MULE,MULTE,MAXDER,NOCONT,NODV,NOPV,THRESH,
     &                  IPRINT,FIRST,SQ12EL,NPCOA,NPCOB,NPCOC,NPCOD,
     &                  NUCSAB,NUCSCD,JSTRA,JSTRB,JSTRC,JSTRD,NINDAB,
     &                  NINDCD,INDHSQ,IODDHR,INDHER,LMNVLS,NPNTAB,
     &                  NPNTCD,NREDAB,NREDCD,IFCTYP,ADISTR,JSTRSH,
     &                  NPRIMS,NCONTS,IORBSH,WRKBUF,
     &                  GABRAO,DMRAO,DMRSO,DINTSKP,IREPDM,IXPAB,IXPCD)
#include "implicit.h"
#include "priunit.h"
#include "r12int.h"
#include "drw2el.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "maxorb.h"
      LOGICAL NOPV, NODV, PERTUR, EXPECT, UNDIFF, DDFOCK, DIRFCK,
     &        DISTRI, NOCONT, LONDON, SPNORB, DIA2SO, ZFS2EL, ONECEN,
     &        FIRST, LAST, SQ12EL, SUSCEP, TPRAB, TPRCD, ADISTR, SOFOCK
      DIMENSION WORK(LWORK), DMAT(*), FMAT(*), GMAT(*), PSO(*), PSA(*),
     &          HESSEE(*), IFCTYP(*), INDHSQ(*), IODDHR(*), INDHER(*),
     &          CONTA(*), CONTB(*), CONTC(*), CONTD(*),
     &          NPCOA(NSETA,2,0:NODCAB), NPCOB(NSETB,2,0:NODCAB),
     &          NPCOC(NSETC,2,0:NODCCD), NPCOD(NSETD,2,0:NODCCD),
     &          JSTRA(NSETA), JSTRB(NSETB), JSTRC(NSETC), JSTRD(NSETD),
     &          COORAB(NUCA*NUCB,3,3,NODCAB),
     &          COORCD(NUCC*NUCD,3,3,NODCCD),
     &          EXPAB(NUCA*NUCB,3,NODCAB), FACAB(NUCA*NUCB,NODCAB),
     &          EXPCD(NUCC*NUCD,3,NODCCD), FACCD(NUCC*NUCD,NODCCD),
     &          NUCSAB(NORBA*NORBB,NODCAB), NUCTAB(8),
     &          NUCSCD(NORBC*NORBD,NODCCD), NUCTCD(8),
     &          NINDAB(NORBA*NORBB,2), NORTAB(8),
     &          NINDCD(NORBC*NORBD,2), NORTCD(8),
     &          LMNVLS(KCKMAX,5,4), NPNTAB(*), NPNTCD(*),
     &          NREDAB(*), NREDCD(*), WRKBUF(*), IORBSH(MXSHEL,MXAOVC),
     &          JSTRSH(MXSHEL,MXAOVC,2), NPRIMS(MXSHEL,MXAOVC,2),
     &          NCONTS(MXSHEL,MXAOVC,2),GABRAO(*),
     &          DMRAO(*),DMRSO(*),DINTSKP(*),IREPDM(*),
     &          IXPAB(*),IXPCD(*)
#include "twocom.h"
#include "inftap.h"
#include "blocks.h"
#include "symmet.h"
#include "dorps.h"
#include "nuclei.h"
#include "twosta.h"
C
Cholesky
#include <ccdeco.h>
Cholesky
C
#include "ibtfun.h"
C
      IF (IPRINT .GE. 5) CALL TITLER('Output from TWOOD1','*',103)
C
      CALL ICOPY(NSETA,NPRIMS(ISHELA,1,IGENAB),MXSHEL,NPCOA(1,1,0),1)
      CALL ICOPY(NSETA,NCONTS(ISHELA,1,IGENAB),MXSHEL,NPCOA(1,2,0),1)
      CALL ICOPY(NSETB,NPRIMS(ISHELB,1,IGENAB),MXSHEL,NPCOB(1,1,0),1)
      CALL ICOPY(NSETB,NCONTS(ISHELB,1,IGENAB),MXSHEL,NPCOB(1,2,0),1)
      CALL ICOPY(NSETC,NPRIMS(ISHELC,1,IGENCD),MXSHEL,NPCOC(1,1,0),1)
      CALL ICOPY(NSETC,NCONTS(ISHELC,1,IGENCD),MXSHEL,NPCOC(1,2,0),1)
      CALL ICOPY(NSETD,NPRIMS(ISHELD,1,IGENCD),MXSHEL,NPCOD(1,1,0),1)
      CALL ICOPY(NSETD,NCONTS(ISHELD,1,IGENCD),MXSHEL,NPCOD(1,2,0),1)
      CALL ICOPY(NSETA,JSTRSH(ISHELA,1,IGENAB),MXSHEL,JSTRA,1)
      CALL ICOPY(NSETB,JSTRSH(ISHELB,1,IGENAB),MXSHEL,JSTRB,1)
      CALL ICOPY(NSETC,JSTRSH(ISHELC,1,IGENCD),MXSHEL,JSTRC,1)
      CALL ICOPY(NSETD,JSTRSH(ISHELD,1,IGENCD),MXSHEL,JSTRD,1)
C
      LAST   = (ISHELA.EQ.IASMAX) .AND. (ISHELB.EQ.IBSMAX) .AND.
     &         (ISHELC.EQ.ICSMAX) .AND. (ISHELD.EQ.IDSMAX)
C
C     Read two-electron densities
C
      IF ((EXPECT .OR. SUSCEP .OR. DIA2SO .OR. ZFS2EL) .AND.
     &     .NOT.NOPV) THEN
         IF (EXPECT) THEN
            READ (LUPAO) (PSO(I),I=1,NPMAT)
         ELSEIF (SUSCEP) THEN
            READ (LUPAO) (PSO(I),I=1,NPMAT)
            READ (LUPAS) (PSA(I),I=1,NPMAT)
         ELSEIF (DIA2SO) THEN
            READ (LUPAO) (PSO(I),I=1,NPMAT)
         ELSEIF (ZFS2EL) THEN   
            READ (LUPAO) (PSO(I),I=1,NPMAT)
         END IF
      END IF
C
C     Overlap distributions for first electron
C     ========================================
C
      TPRAB = .FALSE.
C
      CALL ODCVEC(COORAB,EXPAB,FACAB,CONTA,CONTB,JMAXA,JMAXB,NSETA,
     &            NSETB,NUCA,NUCB,NUCTAB,NORBA,NORBB,NPCOA,NPCOB,NUCSAB,
     &            JSTRA,JSTRB,TCONAB,TPRAB,GENAB,12,THRESH,MAXDER,
     &            MULA,MULB,NODCAB,NORTAB,NINDAB,NPNTAB,NREDAB,
     &            KHKTA,KHKTB,EXPECT,DIRFCK,WORK,LWORK,RPRIAB,RCNTAB,
     &            IXPAB,IPRINT)
C
      IF (ISUM(NODCAB,NUCTAB,1).EQ.0.AND..NOT.LAST) RETURN
C
C     Overlap distributions for second electron
C     =========================================
C
      TPRCD = .FALSE.
      CALL ODCVEC(COORCD,EXPCD,FACCD,CONTC,CONTD,JMAXC,JMAXD,NSETC,
     &            NSETD,NUCC,NUCD,NUCTCD,NORBC,NORBD,NPCOC,NPCOD,NUCSCD,
     &            JSTRC,JSTRD,TCONCD,TPRCD,GENCD,34,THRESH,MAXDER,
     &            MULC,MULD,NODCCD,NORTCD,NINDCD,NPNTCD,NREDCD,
     &            KHKTC,KHKTD,EXPECT,DIRFCK,WORK,LWORK,RPRICD,RCNTCD,
     &            IXPCD,IPRINT)
      IF (ISUM(NODCCD,NUCTCD,1).EQ.0.AND..NOT.LAST) RETURN
C
      CALL GETLMN(LMNVLS,IPRINT)
      KHKTAB = KHKTA*KHKTB
      KHKTCD = KHKTC*KHKTD
      KCKTAB = KCKTA*KCKTB
      KCKTCD = KCKTC*KCKTD
      IF (DIAGAB) KHKTAB = KHKTA*(KHKTA + 1)/2
      IF (DIAGCD) KHKTCD = KHKTC*(KHKTC + 1)/2
      IF (DIACAB) KCKTAB = KCKTA*(KCKTA + 1)/2
      IF (DIACCD) KCKTCD = KCKTC*(KCKTC + 1)/2
      NORBAB = IMXVEC(NORTAB,NODCAB)
      NORBCD = IMXVEC(NORTCD,NODCCD)
      NOABCD = NORBAB*NORBCD
C
C     Allocate work space
C
      KCORBA = 1
      KCORBB = KCORBA + (NORBA + 1)/IRAT
      KCORBC = KCORBB + (NORBB + 1)/IRAT
      KCORBD = KCORBC + (NORBC + 1)/IRAT
      KLAST  = KCORBD + (NORBD + 1)/IRAT
C
      IF (U12INT) THEN
         NCFTYP = 2
      ELSE IF (MAXDER .EQ. 0) THEN
         NCFTYP = 1
      ELSE IF (MAXDER .EQ. 1) THEN
         NCFTYP = 3
      ELSE IF (MAXDER .EQ. 2) THEN
         NCFTYP = 6
      ELSE IF (FINDPT .OR. DPTINT) THEN
         NCFTYP = 6
      ELSE IF (DIA2SO) THEN
         NCFTYP = 9
      END IF
      IF (FINDPT .OR. DPTINT) NCFTYP = 6
      IF (DIA2SO) NCFTYP = 9
      IF (BPH2OO) NCFTYP = 3
      MXUCAB = IMXVEC(NUCTAB,NODCAB)
      MXUCCD = IMXVEC(NUCTCD,NODCCD)
      LCOFAB = MXUCAB*(JMAXA+JMAXB+1)*(JMAXA+1)*(JMAXB+1)*3*NCFTYP
      LCOFCD = MXUCCD*(JMAXC+JMAXD+1)*(JMAXC+1)*(JMAXD+1)*3*NCFTYP
      IF (U12INT) THEN
         KCOFAB = KLAST
         KCOFCD = KCOFAB + LCOFAB + MXUCAB
         KLAST  = KCOFCD + LCOFCD + MXUCCD
         CALL DZERO(WORK(KCOFAB),LCOFAB)
         CALL DZERO(WORK(KCOFCD),LCOFCD)
      ELSE
         KCOFAB = KLAST
         KCOFCD = KCOFAB + LCOFAB
         KLAST  = KCOFCD + LCOFCD
      END IF
C
C     Number of SO integrals
C
      IF (EXPECT .OR. DIRFCK .OR. DDFOCK .OR. DIA2SO .OR. ZFS2EL) THEN
         NINTS  = 0
         NINTMX = 0
      ELSE
         CALL NINTSO(MULE,LONDON,SPNORB,UNDIFF,SOFOCK,DISTRI,SQ12EL,
     &               IPRINT)
      END IF
C
      NSOINT = NOABCD*NINTS
      IF (EXPECT .OR. DIRFCK .OR. SUSCEP .OR. DDFOCK .OR. DIA2SO .OR.
     &           ZFS2EL .OR. LAST .OR. (NSOINT .GT. 0)) THEN
C        NOPP12 is number of different types of operators (WK/UniKA/19-11-2002).
         NOPP12 = 0
         IF (R12EIN) THEN
            NOPP12 = NOPP12 + 1
         ELSE
            IF (V12INT) NOPP12 = NOPP12 + 1
            IF (R12INT) NOPP12 = NOPP12 + 1
            IF (U12INT) NOPP12 = NOPP12 + 1
            IF (U21INT) NOPP12 = NOPP12 + 1
         END IF
         KSOINT = KLAST
         KPNTAO = KSOINT + NSOINT*NOPP12
         KPNTOP = KPNTAO + (  NINTMX*NOPREP + 1)/IRAT
         KPNTNO = KPNTOP + (3*NINTMX*NOPREP + 1)/IRAT
         KPNTRP = KPNTNO + (4*NINTMX*NOPREP + 1)/IRAT
         KPNTLG = KPNTRP + (3*NINTMX*NOPREP + 1)/IRAT
         KLAST  = KPNTLG + (3*NINTMX*NOPREP + 1)/IRAT
         IF(SOFOCK.OR.DIRFCK) THEN
           KDNSBF = KLAST
           KLAST  = KDNSBF + 2*NDMAT
         ELSE
           KDNSBF = KLAST
         ENDIF
C
         LWRK   = LWORK - KLAST + 1
         IF (KLAST .GT. LWORK) CALL STOPIT('TWOOD1',' ',KLAST,LWORK)
         MWFCAB = MAX(MWFCAB,LCOFAB)
         MWFCCD = MAX(MWFCCD,LCOFCD)
         MWPSO  = MAX(MWPSO, NPMAT)
         MWSOIN = MAX(MWSOIN,KPNTAO - KSOINT)
         LWTOT  = LWTOT + KLAST
         MWTOT  = MAX(MWTOT,LWTOT)
C
         CALL SYMLOP(WORK(KSOINT),FMAT,DMAT,NDMAT,GMAT,PSO,PSA,HESSEE,
     &               WORK(KLAST),LWRK,WORK(KCOFAB),WORK(KCOFCD),COORAB,
     &               COORCD,EXPAB,EXPCD,FACAB,FACCD,CONTA,CONTB,CONTC,
     &               CONTD,NSOINT,NPMAT,UNDIFF,PERTUR,LONDON,SPNORB,
     &               DIA2SO,ZFS2EL,EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,
     &               DISTRI,IATOM,MULE,
     &               MULTE,MAXDER,NOCONT,NODV,NOPV,THRESH,IPRINT,FIRST,
     &               LAST,SQ12EL,NPCOA,NPCOB,NPCOC,NPCOD,NUCSAB,NUCSCD,
     &               NINDAB,NINDCD,JSTRA,JSTRB,JSTRC,JSTRD,WORK(KCORBA),
     &               WORK(KCORBB),WORK(KCORBC),WORK(KCORBD),
     &               WORK(KPNTAO),WORK(KPNTOP),WORK(KPNTNO),
     &               WORK(KPNTRP),WORK(KPNTLG),NUCTAB,NUCTCD,INDHSQ,
     &               IODDHR,INDHER,LMNVLS,NPNTAB,NPNTCD,NREDAB,NREDCD,
     &               IFCTYP,ADISTR,IORBSH,WRKBUF,
     &               GABRAO,DMRAO,DMRSO,DINTSKP,WORK(KDNSBF),IREPDM,
     &               IXPAB,IXPCD,CHOINT)
C
         LWTOT = LWTOT - KLAST
      END IF
      RETURN
      END
C  /* Deck symlop */
Cef begin
      SUBROUTINE SYMLOP(SOINT,FMAT,DMAT,NDMAT,GMAT,PSO,PSA,HESSEE,WORK,
     &                  LWORK,COEFAB,COEFCD,COORAB,COORCD,EXPAB,EXPCD,
     &                  FACAB,FACCD,CONTA,CONTB,CONTC,CONTD,NSOINT,
     &                  NPMAT,UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                  EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,IATOM,
     &                  MULE,MULTE,MAXDER,
     &                  NOCONT,NODV,NOPV,THRESH,IPRINT,FIRST,LAST,
     &                  SQ12EL,NPCOA,NPCOB,NPCOC,NPCOD,NUCSAB,NUCSCD,
     &                  NINDAB,NINDCD,JSTRA,JSTRB,JSTRC,JSTRD,ICORBA,
     &                  ICORBB,ICORBC,ICORBD,IPNTAO,IPNTOP,IPNTNO,
     &                  IPNTRP,IPNTLG,NUCTAB,NUCTCD,INDHSQ,IODDHR,
     &                  INDHER,LMNVLS,NPNTAB,NPNTCD,NREDAB,NREDCD,
     &                  IFCTYP,ADISTR,IORBSH,WRKBUF,
     &                  GABRAO,DMRAO,DMRSO,DINTSKP,DNSBUF,IREPDM,
     &                  IXPAB,IXPCD,CHOINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxorb.h"
#include "twocom.h"
#include "symmet.h"
#include "dorps.h"
#include "nuclei.h"
#include "twosta.h"
#include "energy.h"
#include "expopt.h"
#include "drw2el.h"
#include "r12int.h"
#include "doxyz.h"
#include "blocks.h"
#include "dftcom.h"
#include "incore.h"
#include "gnrinf.h"
      PARAMETER (D1 = 1.0D0, D0 = 0.0D0, D2 = 2.0D0, D4 = 4.0D0)
      CHARACTER COMP(2)*1,SPDCAR*1
      LOGICAL AACDX, AACDY, AACDZ, ABCCX, ABCCY, ABCCZ, NOPV, NODV,
     &        PERTUR, EXPECT, INTS, UNDIFF, CAEQCB, CCEQCD, ABADX,
     &        ABADY,  ABADZ, DDFOCK, DIRFCK, DISTRI, NOCONT, SPNORB,
     &        DIA2SO, ZFS2EL, ONECEN, FIRST, LAST, SQ12EL, SOP000,
     &        IPNTLG(3,*), LONDON, SUSCEP, ADISTR,SOFOCK
      DIMENSION SOINT(NSOINT,*), PSO(NPMAT), PSA(NPMAT), HESSEE(*),
     &          WORK(LWORK), DMAT(*), IFCTYP(*),
     &          INDHSQ(*), IODDHR(*), INDHER(*), FMAT(*), GMAT(*),
     &          CONTA(NORBA*NUCA,2,NODCAB), CONTB(NORBB*NUCB,2,NODCAB),
     &          CONTC(NORBC*NUCC,2,NODCCD), CONTD(NORBD*NUCD,2,NODCCD),
     &          NPCOA(NSETA,2,0:NODCAB), NPCOB(NSETB,2,0:NODCAB),
     &          NPCOC(NSETC,2,0:NODCCD), NPCOD(NSETD,2,0:NODCCD),
     &          JSTRA(NSETA), JSTRB(NSETB), JSTRC(NSETC), JSTRD(NSETD),
     &          ICORBA(NORBA),ICORBB(NORBB),ICORBC(NORBC),ICORBD(NORBD),
     &          IPNTAO(*), IPNTOP(3,*), IPNTNO(4,*), IPNTRP(3,*),
     &          COORAB(NUCA*NUCB,3,3,NODCAB),
     &          COORCD(NUCC*NUCD,3,3,NODCCD),
     &          EXPAB(NUCA*NUCB,3,NODCAB), FACAB(NUCA*NUCB,NODCAB),
     &          EXPCD(NUCC*NUCD,3,NODCCD), FACCD(NUCC*NUCD,NODCCD),
     &          NUCSAB(NORBA*NORBB,NODCAB), NUCTAB(8),
     &          NUCSCD(NORBC*NORBD,NODCCD), NUCTCD(8), SIGNT(3),
     &          NINDAB(*), NINDCD(*), LMNVLS(*),
     &          NPNTAB(NUCA*NUCB,2,NODCAB), NPNTCD(NUCC*NUCD,2,NODCCD),
     &          NREDAB(*), NREDCD(*), WRKBUF(*), IORBSH(MXSHEL,MXAOVC),
     &          IEFFB(0:7), IEFFC(0:7), IEFFD(0:7),
     &          GABRAO(NSYMBL,NSYMBL),DMRAO(NSYMBL,NSYMBL,NDMAT),
     &          DMRSO(*),DNSBUF(2,NDMAT),DINTSKP(2,4),IREPDM(*),
     &          IXPAB(NUCA*NUCB,2,NODCAB), IXPCD(NUCC*NUCD,2,NODCCD)
      dimension COEFCD(*)
C
Cholesky
#include "ccdeco2.h"
      LOGICAL CHOINT
Cholesky
C
#include "ibtfun.h"
      IBTAXO(I,J) = IBTAND(I,IBTXOR(I,J))
      XAND(I)     = PT(IBTAND(ISYMAX(1,1),I))
      YAND(I)     = PT(IBTAND(ISYMAX(2,1),I))
      ZAND(I)     = PT(IBTAND(ISYMAX(3,1),I))
      
      CALL QENTER('SYMLOP')
C
      INTS = .FALSE.

C
      LMCORE = MMCORE
      JSCORE = ISCORE
C
      IF (IPRINT .GE. 5) CALL TITLER('Output from SYMLOP','*',103)
C
      IF (.NOT.(EXPECT.OR.DIRFCK.OR.DDFOCK.OR.DIA2SO.OR.ZFS2EL)
     &     .AND. NSOINT.EQ.0 .AND. LAST) THEN
         IF (CHOINT) CALL CPRLOP(IPNTAO,IPNTOP,IPNTNO,IPNTRP,
     &                           IPNTLG,SQ12EL,IPRINT)
         GOTO 400
      END IF
C
Cholesky
C
C
      IF (EXPECT .OR. SUSCEP .OR. DIRFCK .OR. DDFOCK .OR. DIA2SO .OR.
     &    ZFS2EL) THEN
         CALL SETEFF(IEFFB,IEFFC,IEFFD)
      END IF
C
      MULAB  = IBTAND(MULA,MULB)
      MULCD  = IBTAND(MULC,MULD)
C
      IF(DIRAC) THEN
        COMP(1) = 'L'
        COMP(2) = 'S'
      ELSE
        COMP(1) = ' '
        COMP(2) = ' '
      ENDIF
      IF (DIRFCK .OR. SOFOCK) CALL DZERO(DNSBUF,2*NDMAT)
C
C     ***** Both distributions *****
C
      JMAX0  = MAXAB + MAXCD
      HKABCD = FMULT(IBTAND(MULAB,MULCD))
      IF (SPNORB) HKABCD = - HKABCD
      IF (UNDIFF .OR. LONDON .OR. SPNORB .OR. DIRFCK .OR.
     &    DISTRI .OR. SOFOCK .OR. DIA2SO .OR. ZFS2EL .OR. EXPGRA) THEN
         DOX = .TRUE.
         DOY = .TRUE.
         DOZ = .TRUE.
      ELSE
         DOX = DOCOOR(1,NCENTA) .OR. DOCOOR(1,NCENTB)  .OR.
     &         DOCOOR(1,NCENTC) .OR. DOCOOR(1,NCENTD)
         DOY = DOCOOR(2,NCENTA) .OR. DOCOOR(2,NCENTB)  .OR.
     &         DOCOOR(2,NCENTC) .OR. DOCOOR(2,NCENTD)
         DOZ = DOCOOR(3,NCENTA) .OR. DOCOOR(3,NCENTB)  .OR.
     &         DOCOOR(3,NCENTC) .OR. DOCOOR(3,NCENTD)
      END IF
      DOXYZ(1) = DOX
      DOXYZ(2) = DOY
      DOXYZ(3) = DOZ
C
      ICENTA = -1
      ICENTB = -2
      ICENTC = -3
      ICENTD = -4
      SIGNAX = D1
      SIGNAY = D1
      SIGNAZ = D1
      ICENTA = NUCNUM(NCNTSH(ISHELA),1)
      IDEG   = NDEGNM(ICENTA)
      ISYBLA = ISYMBL(ISHELA,IDEG)
      IF (NCENTA .LE. 0) ICENTA = -1
C
C     Symmetrization loop
C     ===================
C
      IF (.NOT.(EXPECT.OR.DIRFCK.OR.DDFOCK.OR.DIA2SO.OR.ZFS2EL)) THEN
         CALL CPRLOP(IPNTAO,IPNTOP,IPNTNO,IPNTRP,IPNTLG,SQ12EL,IPRINT)
      END IF
C
C     **********************************
C     ***** First Symmetry Index R *****
C     **********************************
C
C     Generates distinct overlap distributions A*R(B)
C
      IF (UNDIFF .OR. DIRFCK .OR. DISTRI .OR. SOFOCK) THEN
         SOP000 = .TRUE.
         IF (LAST) CALL DZERO(SOINT,NSOINT*NOPP12)
      ELSE
         SOP000 = .FALSE.
         IF (PERTUR .OR. LONDON .OR. SPNORB .OR. LAST)
     &       CALL DZERO(SOINT,NSOINT)
      END IF
      NSYMR = 0
      DO 100 ISYMR = 0,MAXOPR
C
      IF (IBTAND(ISYMR,IBTOR(MULA,MULB)) .EQ. 0) THEN
         NSYMR = NSYMR + 1
         ICENTB = NUCNUM(NCNTSH(ISHELB),IBTAXO(ISYMR,MULB)+1)
         IDEG   = NDEGNM(ICENTB)
         ISYBLB = ISYMBL(ISHELB,IDEG)
         IF (NCENTB .LE. 0) ICENTB = -2
         NUCAB  = NUCTAB(NSYMR)
         SIGNBX = XAND(ISYMR)
         SIGNBY = YAND(ISYMR)
         SIGNBZ = ZAND(ISYMR)
C
C        ***************************************************
C        ***** Charge Distributions for First Electron *****
C        ***************************************************
C
         TPRIAB = .FALSE.
         CAEQCB = .NOT.BIGVEC .AND. ICENTA .EQ. ICENTB
     &                        .AND. IATOM  .NE. NCENTA
         CALL ODCOEF(COEFAB,COORAB(1,1,1,NSYMR),EXPAB(1,1,NSYMR),WORK,
     &               LWORK,JMAXA,JMAXB,NHKTA,NHKTB,NSETA,NSETB,NUCA,
     &               NUCB,NUCAB,MXUCAB,NORBA,NORBB,NPCOA,NPCOB,
     &               NUCSAB(1,NSYMR),JSTRA,JSTRB,D1,D1,D1,SIGNBX,SIGNBY,
     &               SIGNBZ,CORAX0,CORAY0,CORAZ0,CORBX0,CORBY0,CORBZ0,
     &               AACDX,AACDY,AACDZ,IAB0X,IAB0Y,IAB0Z,CAEQCB,.TRUE.,
     &               .TRUE.,BIGVEC,UNDIFF,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &               12,THRESH,MAXDER,IPRINT)
C
C        **********************************
C        **** Second Symmetry Index S *****
C        **********************************
C
C        Generates distinct overlap distributions C*S(D)
C
         NSYMS = 0
         DO 200 ISYMS = 0, MAXOPR
         IF (IBTAND(ISYMS,IBTOR(MULC,MULD)) .EQ. 0) THEN
            NSYMS = NSYMS + 1
            TPRICD = .FALSE.
C
C           ****************************************************
C           ***** Charge Distributions for Second Electron *****
C           ****************************************************
C
            SIGNDX = XAND(ISYMS)
            SIGNDY = YAND(ISYMS)
            SIGNDZ = ZAND(ISYMS)
            NUCCD = NUCTCD(NSYMS)
C
            CALL ODCOEF(COEFCD,COORCD(1,1,1,NSYMS),EXPCD(1,1,NSYMS),
     &                  WORK,LWORK,JMAXC,JMAXD,NHKTC,NHKTD,NSETC,NSETD,
     &                  NUCC,NUCD,NUCCD,MXUCCD,NORBC,NORBD,NPCOC,NPCOD,
     &                  NUCSCD(1,NSYMS),JSTRC,JSTRD,D1,D1,D1,SIGNDX,
     &                  SIGNDY,SIGNDZ,CORCX0,CORCY0,CORCZ0,CORDX0,
     &                  CORDY0,CORDZ0,ABCCX,ABCCY,ABCCZ,ICD0X,ICD0Y,
     &                  ICD0Z,.FALSE.,.TRUE.,.TRUE.,BIGVEC,UNDIFF,
     &                  LONDON,SPNORB,DIA2SO,ZFS2EL,34,THRESH,
     &                  MAXDER,IPRINT)
            IF (NOCONT) THEN
              NCCINT = NUCAB *NUCCD *KHKTAB*KHKTCD
            ELSE
              NCCINT = NORBAB*NORBCD*KHKTAB*KHKTCD
            END IF
C
C           **********************************
C           ***** Third Symmetry Index T *****
C           **********************************
C
C           Generates distinct quadruplets A*R(B) * T(C*S(D))
C
            DO 300 ISYMT = 0, MAXOPR
            IF (IBTAND(ISYMT,IBTOR(MULAB,MULCD)) .EQ. 0) THEN
              ISYMTS = IBTXOR(ISYMT,ISYMS)
              ICENTC = NUCNUM(NCNTSH(ISHELC),IBTAXO(ISYMT,MULC)+1)
              IDEG   = NDEGNM(ICENTC)
              ISYBLC = ISYMBL(ISHELC,IDEG)
              IF (NCENTC .LE. 0) ICENTC = -3
              ICENTD = NUCNUM(NCNTSH(ISHELD),IBTAXO(ISYMTS,MULD)+1)
              IDEG   = NDEGNM(ICENTD)
              ISYBLD = ISYMBL(ISHELD,IDEG)
              IF (NCENTD .LE. 0) ICENTD = -4
C
C              ****************************************************
C              ***** Charge Distributions for Second Electron *****
C              ****************************************************
C
               SIGNT(1) = XAND(ISYMT)
               SIGNT(2) = YAND(ISYMT)
               SIGNT(3) = ZAND(ISYMT)
               SIGNCX = XAND(ISYMT)
               SIGNCY = YAND(ISYMT)
               SIGNCZ = ZAND(ISYMT)
               SIGNDX = XAND(ISYMTS)
               SIGNDY = YAND(ISYMTS)
               SIGNDZ = ZAND(ISYMTS)
               CCEQCD = .NOT.BIGVEC .AND. (ICENTC .EQ. ICENTD)
     &              .AND. (NCENTC .NE. IATOM)
C     
C              If necessary change sign of expansion coefficients
C              ==================================================
C
               CALL EXCSGN(COEFCD,JMAXC,JMAXD,NHKTC,NHKTD,NUCCD,
     &                     MXUCCD,ICD0X,ICD0Y,ICD0Z,MAXDER,LONDON,
     &                     SPNORB,DIA2SO,EXPGRA,ISYMT,IPRINT)
ckr           Do we need anything fancy changing of sign for ZFS2EL
C
C
C              b) atom to be differentiated does not enter integral
C     
               IF (PERTUR .AND. MAXDER .GT. 0) THEN
                  IF (IATOM .NE. NCENTA .AND. IATOM .NE. NCENTB
     &                                  .AND. IATOM .NE. NCENTC
     &                                  .AND. IATOM .NE. NCENTD) THEN
                     LINTMP = LINTSV
                     GO TO 300
                  END IF
               END IF
C
C              c) one-center integrals
C
               ONECEN = (ICENTA .EQ. ICENTB) .AND.
     &                  (ICENTA .EQ. ICENTC) .AND.
     &                  (ICENTA .EQ. ICENTD) .AND.
     &                  (ICENTA .NE. 0)
               IF (ONECEN) THEN
                  IF (MAXDER .EQ. 0 .OR. EXPGRA) THEN
                     IF (IBTAND(JMAX0,1).EQ.1) THEN
                        LINTMP = LINTSV
                        GO TO 300
                     END IF
                  ELSE
                     IF (PERTUR .OR. EXPECT .OR. LONDON) THEN
                        LINTMP = LINTSV 
                        GO TO 300
                     END IF
                  END IF
               END IF
               IF (LONDON) THEN
                  IF ((ICENTA .EQ. ICENTB) .AND.
     &                (ICENTC .EQ. ICENTD)) THEN
                     LINTMP = LINTSV
                     GO TO 300
                  END IF
               END IF
C
C              d) none of the directions are requested
C
               IF (.NOT. (DOX .OR. DOY .OR. DOZ)) THEN
                  LINTMP = LINTSV
                  GO TO 300
               END IF          
C
C              Check whether this integral gives zero contribution
C              ===================================================
C
C              a) no contribution to direct Fock matrix construction
C
C              Screening: pretty thorough, taking care also of
C                         non-symmetric density matrices
                  IF(DIRFCK.OR.SOFOCK) THEN
                    DINTSKP(1,1) = DINTSKP(1,1) + NCCINT
                    IF(DOSCRN) THEN
                      ULINT = GABRAO(ISYBLA,ISYBLB)*
     &                        GABRAO(ISYBLC,ISYBLD)
                      DNSMAX   = D0
                      DO I = 1,NDMAT
C Largest Coulomb contribution
                        DNSBUF(1,I) = MAX(DMRAO(ISYBLA,ISYBLB,I),
     &                                    DMRAO(ISYBLB,ISYBLA,I),
     &                                    DMRAO(ISYBLC,ISYBLD,I),
     &                                    DMRAO(ISYBLD,ISYBLC,I))
                        DNSMAX = MAX(DNSMAX,D4*DNSBUF(1,I))
C Largest exchange contribution
                        IF (HFXFAC .NE. 0.0D0) THEN
                          DNSBUF(2,I) = MAX(DMRAO(ISYBLC,ISYBLA,I),
     &                                    DMRAO(ISYBLA,ISYBLC,I),
     &                                    DMRAO(ISYBLC,ISYBLB,I),
     &                                    DMRAO(ISYBLB,ISYBLC,I),
     &                                    DMRAO(ISYBLD,ISYBLA,I),
     &                                    DMRAO(ISYBLA,ISYBLD,I),
     &                                    DMRAO(ISYBLB,ISYBLD,I),
     &                                    DMRAO(ISYBLD,ISYBLB,I))
                          DNSMAX = MAX(DNSMAX,DNSBUF(2,I))
                        ELSE
                          DNSBUF(2,I) = 0.0D0
                        END IF
                      ENDDO
                      FCKMAX = DNSMAX*ULINT
                      IF(FCKMAX.LT.FCKTHR) THEN
                        DINTSKP(2,1) = DINTSKP(2,1) + NCCINT
                        IF(IPRINT.GE.5) WRITE(LUPRI,
     &              '(A,1P,D8.2,4(I4,A1,A6,1X,3A1))')
     &        '*Skip: ',FCKMAX,
     &        ISHELA,'(',NAMDEP(ICENTA),COMP(LCLASH(ISHELA)),
     &        SPDCAR(NHKTSH(ISHELA)-1),')',
     &        ISHELB,'(',NAMDEP(ICENTB),COMP(LCLASH(ISHELB)),
     &        SPDCAR(NHKTSH(ISHELB)-1),')',
     &        ISHELC,'(',NAMDEP(ICENTC),COMP(LCLASH(ISHELC)),
     &        SPDCAR(NHKTSH(ISHELC)-1),')',
     &        ISHELD,'(',NAMDEP(ICENTD),COMP(LCLASH(ISHELD)),
     &        SPDCAR(NHKTSH(ISHELD)-1),')'
#if defined (VAR_MPI)
          IF (AOSAVE .AND. LINTSV .AND. MAXDER .EQ. 0 .AND. 
     &        DIRFCK .AND. HFXMU .EQ. D0) THEN
             I_SHL = I_SHL + 1
             IF (I_SHL .GT. MXTSK) THEN
                write(lupri,*) 'Received quit because '//
     &               'I_SHL exceeds MXTSK',I_SHL,MXTSK
                CALL QUIT('AOSAVE: I_SHL exceeds MXTSK.')
             END IF
#else
                  IF (AOSAVE .AND. INDX_SHL1 .GE. ISHELA .AND. 
     &                 INDX_SHL2 .GE. ISHELB .AND.
     &                 INDX_SHL3 .GE. ISHELC .AND. 
     &                 INDX_SHL4 .GE. ISHELD .AND. IPREVA .EQ. -1
     &                 .AND. HFXMU .EQ. D0) THEN
                     NINTYP = 1
                     I_SHL = I_SHL + NINTYP*NCCINT
                     IF (I_SHL .GT. MMCORE) THEN
                        write(lupri,*) 'Received quit because '//
     &                       'I_SHL exceeds MMCORE',I_SHL,MMCORE
                        CALL QUIT('I_SHL exceeds MMCORE.')
                     END IF
#endif
                     LINTMP = .TRUE.
                  END IF
                        GOTO 300
                      ENDIF
                    ENDIF
                  ENDIF
C     
C              Integral contributes
C              ====================
C
               IF (IPRINT .GE. 5) THEN
                  WRITE (LUPRI, '(/,1X,A,3I5)')
     &               ' Symmetry operations:',ISYMR, ISYMS, ISYMT
                  WRITE (LUPRI, '(/,1X,A,I5)') ' ISYMTS ', ISYMTS
                  WRITE (LUPRI, '(1X,A,4I5)')
     &               ' ICENTA-D ',ICENTA, ICENTB, ICENTC, ICENTD
                  WRITE (LUPRI, '(1X,A,4I5)')
     &               ' NCENTA-D ',NCENTA, NCENTB, NCENTC, NCENTD
               END IF
C
C              Local symmetries
C              ================
C
               IF (BIGVEC) THEN
                  CALL CORDIF(NSETA,NSETC,THRESH,ABADX,ABADY,ABADZ,
     &                        IPRINT,NPCOA(1,1,NSYMR),NPCOC(1,1,NSYMS),
     &                        JSTRA,JSTRC)
               ELSE
                  ABADX  = ABS(CORAX0 - XAND(ISYMT)*CORCX0) .LT. THRESH
                  ABADY  = ABS(CORAY0 - YAND(ISYMT)*CORCY0) .LT. THRESH
                  ABADZ  = ABS(CORAZ0 - ZAND(ISYMT)*CORCZ0) .LT. THRESH
               END IF
               ISAMEX = 0
               ISAMEY = 0
               ISAMEZ = 0
               IF (AACDX .AND. ABADX .AND. ABCCX) ISAMEX = 1
               IF (AACDY .AND. ABADY .AND. ABCCY) ISAMEY = 1
               IF (AACDZ .AND. ABADZ .AND. ABCCZ) ISAMEZ = 1
               ISMXYZ = ISAMEX + 2*ISAMEY + 4*ISAMEZ
C
               IF (EXPECT .OR. SUSCEP .OR. DIRFCK .OR. DDFOCK
     &                    .OR. DIA2SO .OR. ZFS2EL) THEN
                  DO 55 I = 1, NORBA
                     ICORBA(I) = IORBSH(ISHELA,I)
   55             CONTINUE
                  DO 65 I = 1, NORBB
                     ICORBB(I) = IORBSH(ISHELB,I) + IEFFB(ISYMR)
   65             CONTINUE
                  DO 75 I = 1, NORBC
                     ICORBC(I) = IORBSH(ISHELC,I) + IEFFC(ISYMT)
   75             CONTINUE
                  DO 85 I = 1, NORBD
                     ICORBD(I) = IORBSH(ISHELD,I) + IEFFD(ISYMTS)
   85             CONTINUE
               END IF
C
#if defined (VAR_MPI)
               IF (AOSAVE .AND. DIRFCK .AND. MAXDER .EQ. 0 .AND.
     &             HFXMU .EQ. D0) THEN
C
                  NINTYP = 1
                  IF (LINTSV) THEN
C
C                      write(lupri,*) 'LINTSV,ISHELA/B/C/D,ISYMR,ISYMS',
C     &                 LINTSV,ISHELA,ISHELB,ISHELC,ISHELD,
C     &                 ISYMR,ISYMS
                     LINTMP = .TRUE.
                     KAOINT = 1
C
                     CALL DCOPY(NCCINT*NINTYP,
     &                    AOINTSCORE(INDX_SHL(I_SHL)),1,WORK(KAOINT),1)
C                     write(lupri,*) 'I_SHL,INDX_SHL(I_SHL)', 
C     &                    I_SHL,INDX_SHL(I_SHL)
C
C                     write(lupri,*) 'WORK(KAOINT): ',(WORK(KAOINT+II),
C     &                    II=0,NCCINT-1)
                     CALL INTFCK(FMAT,WORK(KAOINT),DMAT,NDMAT,
     &                    NCCINT,NINTYP,ICORBA,ICORBB,ICORBC,ICORBD,
     &                    HKABCD,IPRINT,NODV,NINDAB,NINDCD,SUSCEP,
     &                    IFCTYP,DNSBUF,DINTSKP,
     &                    WORK(KAOINT+NCCINT*NINTYP),
     &                    (LWORK-NCCINT*NINTYP))
C
                     I_SHL = I_SHL + 1
                     IF (I_SHL .GT. MXTSK) THEN
                        write(lupri,*) 'Received quit because '//
     &                       'I_SHL exceeds MXTSK',I_SHL,MXTSK
                        CALL QUIT('AOSAVE: I_SHL exceeds MXTSK.')
                     END IF
C
                     GOTO 300
                  END IF
C
               ELSE IF (MAXDER .GT. 0) THEN
C                  write(lupri,*) 'lintmp,lintsv,maxder',lintmp,lintsv,
C     &                 maxder
                  LINTMP = LINTSV
               END IF
#else
               IF (AOSAVE .AND. DIRFCK .AND. MAXDER .EQ. 0 .AND.
     &             HFXMU .EQ. D0) THEN
                  NINTYP = 1

                  IF (INDX_SHL1 .GE. ISHELA .AND. 
     &                INDX_SHL2 .GE. ISHELB .AND.
     &                INDX_SHL3 .GE. ISHELC .AND. 
     &                INDX_SHL4 .GE. ISHELD .AND. IPREVA .EQ. -1) THEN
C
C                  IF (LINTSV) THEN
C                   write(lupri,*)'ISHELA/B/C/D,INDX_SHL1/2/3/4,lintsv',
C     &                    ishela,ishelb,ishelc,isheld,indx_shl1,
C     &                    indx_shl2,indx_shl3,indx_shl4,lintsv
C
                     LINTMP = .TRUE.
                     KAOINT = 1 
C
                     CALL DCOPY(NCCINT*NINTYP,
     &                    AOINTSCORE(I_SHL),1,WORK(KAOINT),1)   
C
                     CALL INTFCK(FMAT,WORK(KAOINT),DMAT,NDMAT,
     &                    NCCINT,NINTYP,ICORBA,ICORBB,ICORBC,ICORBD,
     &                    HKABCD,IPRINT,NODV,NINDAB,NINDCD,SUSCEP,
     &                    IFCTYP,DNSBUF,DINTSKP,
     &                    WORK(KAOINT+NCCINT*NINTYP),
     &                    (LWORK-NCCINT*NINTYP))
C
                     I_SHL = I_SHL + NCCINT*NINTYP
                     IF (I_SHL .GT. MMCORE) THEN
                        write(lupri,*) 'Received quit because '//
     &                       'I_SHL exceeds MMCORE',I_SHL,MMCORE
                        CALL QUIT('I_SHL exceeds MMCORE.')
                     END IF
C
                     GOTO 300
                  END IF
C
               ELSE IF (MAXDER .GT. 0) THEN
                  LINTMP = LINTSV
               END IF
#endif
C
C
C              *******************************
C              ***** Integral Directives *****
C              *******************************
C
               CALL DIRECT(BIGVEC,ICENTA,ICENTB,ICENTC,ICENTD,
     &              NCENTA,NCENTB,NCENTC,NCENTD,
     &              0,     ISYMR, ISYMT, ISYMTS,
     &              SIGNAX,SIGNAY,SIGNAZ,SIGNBX,SIGNBY,SIGNBZ,
     &              SIGNCX,SIGNCY,SIGNCZ,SIGNDX,SIGNDY,SIGNDZ,
     &              NCCINT,MAXDER,EXPECT,LONDON,SPNORB,DIA2SO,
     &              ZFS2EL,EXPGRA,IATOM,MULTE,NINTYP,IPRINT)
C
C              ************************
C              ***** AO integrals *****
C              ************************
C     
               KAOINT = 1
               LAOINT = NCCINT*NINTYP
               KLAST  = KAOINT + LAOINT
               LWRK   = LWORK  - KLAST + 1
               IF (KLAST.GT.LWORK) CALL STOPIT('CAOINT',' ',KLAST,LWORK)
               MWAOIN = MAX(MWAOIN,KLAST - KAOINT)
               LWTOT  = LWTOT + KLAST
               MWTOT  = MAX(MWTOT,LWTOT)
C
               CALL CAOINT(SOINT,DMAT,NDMAT,PSO,PSA,FMAT,WORK(KAOINT),
     &              WORK(KLAST),LWRK,COEFAB,COEFCD,
     &              COORAB(1,1,1,NSYMR),COORCD(1,1,1,NSYMS),
     &              EXPAB(1,1,NSYMR),EXPCD(1,1,NSYMS),
     &              FACAB(1,NSYMR),FACCD(1,NSYMS),
     &              CONTA(1,1,NSYMR),CONTB(1,1,NSYMR),
     &              CONTC(1,1,NSYMS),CONTD(1,1,NSYMS),
     &              NCCINT,NINTYP,UNDIFF,PERTUR,LONDON,SPNORB,
     &              DIA2SO,ZFS2EL,EXPECT,SUSCEP,DDFOCK,DIRFCK,
     &              SOFOCK,DISTRI,IATOM,
     &              MULE,MULTE,MAXDER,BIGVEC,NOCONT,NODV,NOPV,
     &              ISMXYZ,THRESH,ONECEN,IPRINT,ICORBA,ICORBB,
     &              ICORBC,ICORBD,INTS,HKABCD,JMAX0,ISYMR,ISYMS,
     &              ISYMT,ISYMTS,SQ12EL,SOP000,IPNTAO,IPNTOP,
     &              NPCOA(1,1,NSYMR),NPCOB(1,1,NSYMR),
     &              NPCOC(1,1,NSYMS),NPCOD(1,1,NSYMS),
     &              NUCSAB(1,NSYMR), NUCSCD(1,NSYMS),
     &              SIGNT,INDHSQ,IODDHR,INDHER,LMNVLS,NINDAB,
     &              NINDCD,NPNTAB(1,1,NSYMR),NPNTCD(1,1,NSYMS),
     &              NREDAB,NREDCD,IFCTYP,HESSEE,DNSBUF,DINTSKP,
     &              NSOINT,IXPAB(1,1,NSYMR),IXPCD(1,1,NSYMS))
               LWTOT  = LWTOT - KLAST
C               
C               write(lupri,*) 'WORK(KAOINT)2: ',(WORK(KAOINT+II-1),
C     &                    II=1,NCCINT)
               IF (N_SHL .GT. MXTSK) THEN
                  write(lupri,*) 'Received quit because '//
     &                 'N_SHL exceeds MXTSK',N_SHL,MXTSK
                  CALL QUIT('AOSAVE: N_SHL exceeds MXTSK')
               END IF 
               ITOTNT = ITOTNT + NINTYP*NCCINT
               IF (AOSAVE
     &              .AND. MMCORE .GT. (NCCINT*NINTYP)
     &              .AND. N_SHL  .LE.  MXTSK
     &              .AND. MSAVE  .AND. MAXDER .EQ. 0
     &              .AND. HFXMU  .EQ.  D0) THEN
C
C                  write(lupri,*) 'SAV:LINTSV,ISHELA/B/C/D,ISYMR,ISYMS',
C     &                 LINTSV,ISHELA,ISHELB,ISHELC,ISHELD,
C     &                 ISYMR,ISYMS
                  LINTMP = .TRUE.
                  MMCORE = MMCORE - (NCCINT*NINTYP)
                  CALL DCOPY((NCCINT*NINTYP),WORK(KAOINT),1,
     &               AOINTSCORE(ISCORE),1)
C                  INDX_SHL(1,1) = ISHELA
C                  INDX_SHL(2,1) = ISHELB
C                  INDX_SHL(3,1) = ISHELC
C                  INDX_SHL(4,1) = ISHELD
                  IPREVA = INDX_SHL1
                  IPREVB = INDX_SHL2
                  IPREVC = INDX_SHL3
                  IPREVD = INDX_SHL4
                  INDX_SHL1 = ISHELA
                  INDX_SHL2 = ISHELB
                  INDX_SHL3 = ISHELC
                  INDX_SHL4 = ISHELD
#if defined (VAR_MPI)                  
                  INDX_SHL(N_SHL) = ISCORE
C                  write (lupri,*) 'n_shl,indx',n_shl,indx_shl(n_shl)
                  N_SHL = N_SHL + 1
                  ISCORE = ISCORE + (NINTYP*NCCINT)
                  IF (N_SHL .GT. (MXTSK+1)) THEN
                     write(lupri,*) 'Received quit because '//
     &                  'N_SHL exceeds MXTSK+1',N_SHL,MXTSK+1
                     CALL QUIT('AOSAVE: N_SHL exceeds MXTSK+1.')
                  END IF
#else
                  ISCORE = ISCORE + (NINTYP*NCCINT)
#endif                  
C     Lagrer integraler i minne: Oppdaterer array med beregnede ishela/b/c/d
C  =  Store integrals in memory: Update array with calculated ishela/b/c/d
               ELSE
                  IF (MAXDER .GT. 0) THEN
                     LINTMP = LINTSV
                  ELSE
                     LINTMP = .FALSE.
                     INDX_SHL1 =IPREVA
                     INDX_SHL2 =IPREVB
                     INDX_SHL3 =IPREVC
                     INDX_SHL4 =IPREVD
                  END IF
                  MSAVE = .FALSE.
C        Reset MMCORE,ISCORE since we did not save anything useful in this call.
                  MMCORE = LMCORE
                  ISCORE = JSCORE
C
               END IF
            END IF
 300     CONTINUE
      END IF
 200  CONTINUE
      END IF
C
 100  CONTINUE
C
C     ********************************
C     ***** Process SO integrals *****
C     ********************************
C
 400  CONTINUE
      IF (.NOT. (EXPECT.OR.DIRFCK.OR.DDFOCK.OR. DIA2SO .OR. ZFS2EL)
     &    .AND. (INTS.OR.LAST)) THEN
C
C        A) Undifferentiated integrals
C        =============================
C
         IF (UNDIFF) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            IF (.NOT.ADISTR) THEN
               IF (DCCR12) THEN
                  CALL RN2OUT(SOINT,NSOINT,IPNTNO,IPNTRP,IPNTLG,FIRST,
     &                        LAST,THRESH,NINDAB,NINDCD,IPRINT)
               ELSE IF (SRINTS) THEN
                  CALL SR2OUT(SOINT,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,
     &                        THRESH,NINDAB,NINDCD,IPRINT)
               ELSE
                  CALL UN2OUT(SOINT,IPNTNO,IPNTRP,IPNTLG,FIRST,LAST,
     &                        THRESH,NINDAB,NINDCD,IPRINT)
               END IF
            ELSE
C           ... ADISTR true: Distributions (all gabcd for fixed a)
               IF (U21INT) THEN
                  CALL US2OUT(SOINT,NSOINT,WRKBUF,IPNTNO,IPNTRP,IPNTLG,
     &                        FIRST,LAST,THRESH,NINDAB,NINDCD,
     &                        IORBSH,IPRINT)
               ELSE
                  CALL DS2OUT(SOINT,WRKBUF,IPNTNO,IPNTRP,IPNTLG,
     &                        FIRST,LAST,THRESH,NINDAB,NINDCD,
     &                        IORBSH,IPRINT)
               END IF
            END IF
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TSYMOU = TSYMOU + TIMEND - TIMSTR
            END IF
            FIRST = .FALSE.
C
C        B) Differentiated integrals
C        ===========================
C
         ELSE IF (PERTUR .OR. SPNORB) THEN
#ifndef VAR_INT64
           IF (NBASIS .GT. 255) THEN
              CALL QENTER('DR2OUT')
              CALL QUIT('NBASIS .gt. 255 not implemented for 32-bit'//
     &        ' integers - recompile with 64-bit integers.')
           END IF
#endif
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            CALL DR2OUT(SOINT,IATOM,MULE,FIRST,LAST,SPNORB,SQ12EL,
     &                  THRESH,IPNTNO,IPNTRP,IPNTLG,NINDAB,NINDCD,
     &                  IPRINT)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TDSOUT = TDSOUT + TIMEND - TIMSTR
            END IF
            FIRST = .FALSE.
C
C        C) London orbitals
C        ==================
C
         ELSE IF (LONDON) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            CALL MG2OUT(SOINT,FIRST,LAST,THRESH,IPNTNO,IPNTRP,IPNTLG,
     &                  NINDAB,NINDCD,IPRINT)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TDSOUT = TDSOUT + TIMEND - TIMSTR
            END IF
            FIRST = .FALSE.
C
C        D) Direct calculation of Fock matrices in SO basis
C        ==================================================
C
         ELSE IF (SOFOCK) THEN
            CALL QUIT(
     &      'SOFOCK call of FCKOUT not implemented in this version')
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
CSCR            CALL FCKOUT(FMAT,DMAT,NDMAT,SOINT,IPNTNO,IPNTRP,IPNTLG,
CSCR     &                  NINDAB,NINDCD,IFCTYP,DINTSKP,IREPDM,
CSCR     &                  DMRSO,DNSBUF,WORK,LWORK,IPRINT)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TFCKOU = TFCKOU + TIMEND - TIMSTR
            END IF
            FIRST = .FALSE.
C
C        E) Calculation of distributions
C        ===============================
C
         ELSE IF (DISTRI) THEN
            CALL QUIT(
     &      'DISTRI call of FCKOUT not implemented in this version')
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
C     TODO: Make FCKOUT general enough to handle distributions
C           CALL FCKOUT(FMAT,DMAT,NDMAT,GMAT,SOINT,DIRFCK,DISTRI,
C    &                  FIRST,LAST,THRESH,SQ12EL,IPRINT)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TFCKOU = TFCKOU + TIMEND - TIMSTR
            END IF
            FIRST = .FALSE.
         END IF
      ELSE IF (EXPECT) THEN
         IF (IPRINT .GE. 6) 
     & 	    CALL PRIGRD(GRADEE,WORK,WORK(1+MXCOOR*MXCOOR))
      END IF
C
      CALL QEXIT('SYMLOP')
      RETURN
C     
      END
Cef end
C  /* Deck caoint */
      SUBROUTINE CAOINT(SOINT,DMAT,NDMAT,PSO,PSA,FMAT,AOINT,WORK,LWORK,
     &                  COEFAB,COEFCD,COORAB,COORCD,EXPAB,EXPCD,FACAB,
     &                  FACCD,CONTA,CONTB,CONTC,CONTD,NCCINT,NINTYP,
     &                  UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                  EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,IATOM,
     &                  MULE,MULTE,MAXDER,
     &                  BIGVEC,NOCONT,NODV,NOPV,ISMXYZ,THRESH,ONECEN,
     &                  IPRINT,ICORBA,ICORBB,ICORBC,ICORBD,INTS,HKABCD,
     &                  JMAX0,ISYMR,ISYMS,ISYMT,ISYMTS,SQ12EL,SOP000,
     &                  IPNTAO,IPNTOP,NPCOA,NPCOB,NPCOC,NPCOD,NUCSAB,
     &                  NUCSCD,SIGNT,INDHSQ,IODDHR,INDHER,LMNVLS,NINDAB,
     &                  NINDCD,NPNTAB,NPNTCD,NREDAB,NREDCD,IFCTYP,
     &                  HESSEE,DNSBUF,DINTSKP,NSOINT,IXPAB,IXPCD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "r12int.h"
#include "ccom.h"
#include "symmet.h"
#include "twocom.h"
#include "twosta.h"
      LOGICAL NOPV, NODV, PERTUR, EXPECT, NOINT, DDFOCK, DIRFCK, DISTRI,
     &        INTS, UNDIFF, NOCONT, SPNORB, DIA2SO, BIGVEC, ONECEN,
     &        SQ12EL, SOP000, LONDON, SUSCEP, SOFOCK, ZFS2EL
      DIMENSION AOINT(NCCINT,NINTYP), SOINT(NSOINT,*), PSO(*),PSA(*),
     &          COEFAB(*), COEFCD(*), COORAB(*), COORCD(*), EXPAB(*),
     &          EXPCD(*), CONTA(*), CONTB(*), CONTC(*), CONTD(*),
     &          FACAB(*), FACCD(*), NPCOA(*), NPCOB(*), NPCOC(*),
     &          NPCOD(*), NUCSAB(*), NUCSCD(*), ICORBA(*), ICORBB(*),
     &          ICORBC(*), ICORBD(*), IPNTAO(*), IPNTOP(3,*),
     &          SIGNT(3), INDHSQ(*), IODDHR(*), INDHER(*), LMNVLS(*),
     &          NPNTAB(*), NPNTCD(*), NREDAB(*), NREDCD(*),
     &          FMAT(NBASIS,NBASIS,*), DMAT(NBASIS,NBASIS,*), IFCTYP(*),
     &          HESSEE(*), DNSBUF(*),DINTSKP(*), IXPAB(*), IXPCD(*),
     &          WORK(LWORK)
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IF (IPRINT .GE. 10) THEN
         CALL TITLER('Output from CAOINT','*',103)
         WRITE (LUPRI, 1005) NHKTA, NHKTB, NHKTC, NHKTD
         WRITE (LUPRI, 1010) MAXAB, MAXCD, JMAX0
         WRITE (LUPRI, 1020) KHKTA, KHKTB, KHKTC, KHKTD
         WRITE (LUPRI, 1021) KCKTA, KCKTB, KCKTC, KCKTD
         WRITE (LUPRI, 1030) KHKTAB, KHKTCD
         WRITE (LUPRI, 1031) KCKTAB, KCKTCD
         WRITE (LUPRI, 1050) DIAGAB, IAB0X, IAB0Y, IAB0Z
         WRITE (LUPRI, 1060) DIAGCD, ICD0X, ICD0Y, ICD0Z
         WRITE (LUPRI, 1075) TPRIAB, TCONAB
         WRITE (LUPRI, 1076) TPRICD, TCONCD
         WRITE (LUPRI, 1080) ISMXYZ
         WRITE (LUPRI, 1090) NORBA, NORBB, NORBC, NORBD
         WRITE (LUPRI, 1100) NORBAB, NORBCD
         WRITE (LUPRI, 1110) NUCAB, NUCCD
      END IF
C
C     **********************************
C     ***** Calculate AO integrals *****
C     **********************************
C
      CALL INTDER(AOINT,NCCINT,NINTYP,WORK,LWORK,COEFAB,COEFCD,COORAB,
     &            COORCD,EXPAB,EXPCD,FACAB,FACCD,CONTA,CONTB,CONTC,
     &            CONTD,IPRINT,NOINT,NOCONT,SPNORB,DIA2SO,ZFS2EL,MAXDER,
     &            THRESH,JMAX0,ISMXYZ,ONECEN,NPCOA,NPCOB,NPCOC,NPCOD,
     &            NUCSAB,NUCSCD,SIGNT,IODDHR,INDHSQ,INDHER,LMNVLS,
     &            NPNTAB,NPNTCD,NREDAB,NREDCD)
      INTS = INTS .OR. .NOT.NOINT
C
C     ********************************
C     ***** Process AO integrals *****
C     ********************************
C
      IF (.NOT.NOINT) THEN
C
C        A) Contributions to SO integrals
C        ================================
C
C        a) Undifferentiated integrals
C
         IF (UNDIFF .OR. DISTRI. OR. SOFOCK) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            IF (R12EIN) THEN
               CALL SYM2(SOINT,AOINT,IPNTAO,IPNTOP,
     &                   ISYMR,ISYMT,ISYMTS,HKABCD,
     &                   SQ12EL,SOP000,NINTS,IPRINT)
            ELSE
               IX = 1
               IF (V12INT) THEN
                  CALL SYM2(SOINT(1,IX),AOINT(1,IX),IPNTAO,IPNTOP,
     &                      ISYMR,ISYMT,ISYMTS,HKABCD,
     &                      SQ12EL,SOP000,NINTS,IPRINT)
                  IX = IX + 1
               END IF
               IF (R12INT) THEN
                  CALL SYM2(SOINT(1,IX),AOINT(1,IX),IPNTAO,IPNTOP,
     &                      ISYMR,ISYMT,ISYMTS,HKABCD,
     &                      SQ12EL,SOP000,NINTS,IPRINT)
                  IX = IX + 1
               END IF
               IF (U12INT) THEN
                  CALL SYM2(SOINT(1,IX),AOINT(1,IX),IPNTAO,IPNTOP,
     &                      ISYMR,ISYMT,ISYMTS,HKABCD,
     &                      SQ12EL,SOP000,NINTS,IPRINT)
                  IX = IX + 1
               END IF
               IF (U21INT) THEN
                  CALL SYM2(SOINT(1,IX),AOINT(1,IX),IPNTAO,IPNTOP,
     &                      ISYMR,ISYMT,ISYMTS,HKABCD,
     &                      SQ12EL,SOP000,NINTS,IPRINT)
               END IF
            END IF
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TSYM2S = TSYM2S + TIMEND - TIMSTR
            END IF
            SOP000 = .FALSE.
C
C        b) Spin-orbit integrals
C
         ELSE IF (SPNORB) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            IF (.NOT. DDFOCK) THEN
               CALL SPOSYM(SOINT,AOINT,NCCINT,HKABCD,ISYMR,ISYMT,ISYMTS,
     &                  IPNTAO,IPNTOP,SQ12EL,SOP000,IPRINT)
            END IF
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TSPOSY = TSPOSY + TIMEND - TIMSTR
            END IF
            SOP000 = .FALSE.
            IF (DDFOCK) THEN
               IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
               IADRT = 1
               IADRA = 4
               CALL INTEXP(HESSEE,AOINT,DMAT,NDMAT,PSO,PSA,
     &                     FMAT(1,1,IADRT),FMAT(1,1,IADRA),NINDAB,
     &                     NINDCD,NCCINT,NINTYP,WORK,LWORK,ISYMR,ISYMS,
     &                     ISYMT,ICORBA,ICORBB,ICORBC,ICORBD,THRESH,
     &                     HKABCD,IPRINT,NOPV,NODV,EXPECT,LONDON,SUSCEP,
     &                     DIA2SO,ZFS2EL,DDFOCK,SPNORB,IXPAB,IXPCD,
     &                     IFCTYP)
               IF (TKTIME) THEN
                  CALL GETTIM(TIMEND,DUMTIM)
                  TINTEX = TINTEX + TIMEND - TIMSTR
               END IF
            END IF
C
C        c) First derivative integrals
C
         ELSE IF (PERTUR .AND. .NOT.DDFOCK) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            CALL DRSYM2(SOINT,AOINT,WORK,NCCINT,LWORK,HKABCD,ISYMR,
     &                  ISYMT,ISYMTS,IPNTAO,IPNTOP,MULE,SQ12EL,SOP000,
     &                  IPRINT)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TDRSYM = TDRSYM + TIMEND - TIMSTR
            END IF
            SOP000 = .FALSE.
C
C        d) London orbitals
C
         ELSE IF (LONDON) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            IF (.NOT. DDFOCK) THEN
               CALL MGSYM2(SOINT,AOINT,WORK,NCCINT,LWORK,HKABCD,ISYMR,
     &                     ISYMT,ISYMTS,IPNTAO,IPNTOP,SOP000,IPRINT)
            END IF
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TDRSYM = TDRSYM + TIMEND - TIMSTR
            END IF
            SOP000 = .FALSE.
C
            IF (SUSCEP .OR. DDFOCK) THEN
               IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
               IADRT = 1
               IADRA = 4
               CALL INTEXP(HESSEE,AOINT,DMAT,NDMAT,PSO,PSA,
     &                     FMAT(1,1,IADRT),FMAT(1,1,IADRA),NINDAB,
     &                     NINDCD,NCCINT,NINTYP,WORK,LWORK,ISYMR,ISYMS,
     &                     ISYMT,ICORBA,ICORBB,ICORBC,ICORBD,THRESH,
     &                     HKABCD,IPRINT,NOPV,NODV,EXPECT,LONDON,SUSCEP,
     &                     DIA2SO,ZFS2EL,DDFOCK,SPNORB,IXPAB,IXPCD,
     &                     IFCTYP)
               IF (TKTIME) THEN
                  CALL GETTIM(TIMEND,DUMTIM)
                  TINTEX = TINTEX + TIMEND - TIMSTR
               END IF
            END IF
C
C        e) Contributions to expectation values
C        ======================================
C
         ELSE IF (EXPECT .OR. DDFOCK) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            IF (PERTUR) THEN
               IADRT = 1
               IADRA = 1 + 3*NUCDEG(IATOM)
            ELSE
               IADRT = 1
               IADRA = 1 + 3*NUCDEP
            END IF
            CALL INTEXP(HESSEE,AOINT,DMAT,NDMAT,PSO,PSA,FMAT(1,1,IADRT),
     &                  FMAT(1,1,IADRA),NINDAB,NINDCD,NCCINT,
     &                  NINTYP,WORK,LWORK,ISYMR,ISYMS,ISYMT,ICORBA,
     &                  ICORBB,ICORBC,ICORBD,THRESH,HKABCD,IPRINT,
     &                  NOPV,NODV,EXPECT,LONDON,SUSCEP,DIA2SO,ZFS2EL,
     &                  DDFOCK,SPNORB,IXPAB,IXPCD,IFCTYP)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TINTEX = TINTEX + TIMEND - TIMSTR
            END IF
C
C        f) Contributions to Fock matrices
C        =================================
C
         ELSE IF (DIRFCK) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            CALL INTFCK(FMAT,AOINT,DMAT,NDMAT,NCCINT,NINTYP,
     &                  ICORBA,ICORBB,ICORBC,ICORBD,HKABCD,
     &                  IPRINT,NODV,NINDAB,NINDCD,SUSCEP,IFCTYP,
     &                  DNSBUF,DINTSKP,WORK,LWORK)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TFCKOU = TFCKOU + TIMEND - TIMSTR
            END IF
C
C        g) Contributions to expectation values, dia2so
C        ==============================================
C
         ELSE IF (DIA2SO) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            CALL INTEXP(HESSEE,AOINT,DMAT,NDMAT,PSO,PSA,DUM,
     &                  DUM,NINDAB,NINDCD,NCCINT,
     &                  NINTYP,WORK,LWORK,ISYMR,ISYMS,ISYMT,ICORBA,
     &                  ICORBB,ICORBC,ICORBD,THRESH,HKABCD,IPRINT,
     &                  NOPV,NODV,EXPECT,LONDON,SUSCEP,DIA2SO,ZFS2EL,
     &                  DDFOCK,SPNORB,IXPAB,IXPCD,IFCTYP)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TINTEX = TINTEX + TIMEND - TIMSTR
            END IF
C
C        h) Zero-field electronic spin-spin couplings
C        ============================================
C
         ELSE IF (ZFS2EL) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            CALL INTEXP(HESSEE,AOINT,DMAT,NDMAT,PSO,PSA,DUM,
     &                  DUM,NINDAB,NINDCD,NCCINT,
     &                  NINTYP,WORK,LWORK,ISYMR,ISYMS,ISYMT,ICORBA,
     &                  ICORBB,ICORBC,ICORBD,THRESH,HKABCD,IPRINT,
     &                  NOPV,NODV,EXPECT,LONDON,SUSCEP,DIA2SO,ZFS2EL,
     &                  DDFOCK,SPNORB,IXPAB,IXPCD,IFCTYP)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TINTEX = TINTEX + TIMEND - TIMSTR
            END IF
         ELSE
            CALL QUIT('Programming error in CAOINT please report.')
         END IF
      END IF
      RETURN
C
 1005 FORMAT(//,'  NHKTA-D   ',4I5)
 1010 FORMAT(   '  JMAX0     ',I1,' + ',I1,' = ',I2)
 1020 FORMAT(   '  KHKTA-D   ',4I5)
 1021 FORMAT(   '  KCKTA-D   ',4I5)
 1030 FORMAT(   '  KHKTAB-CD ',2I5)
 1031 FORMAT(   '  KCKTAB-CD ',2I5)
 1050 FORMAT(   '  DIAGAB    ',L5,3I5)
 1060 FORMAT(   '  DIAGCD    ',L5,3I5)
 1075 FORMAT(   '  TPRI/CONAB',2L5)
 1076 FORMAT(   '  TPRI/CONCD',2L5)
 1080 FORMAT(   '  ISMXYZ    ', I5)
 1090 FORMAT(   '  NORBA-D   ',4I5)
 1100 FORMAT(   '  NORBAB-CD ',2I5)
 1110 FORMAT(   '  NUCAB-CD  ',2I5)
      END
C  /* Deck paovec */
      SUBROUTINE PAOVEC(JSTRSH,NPRIMS,NCONTS,IORBSH,JORBSH,KORBSH,
     &                  IATOM,NOBLK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "aovec.h"
#include "maxorb.h"
#include "mxcent.h"
      LOGICAL ORBOUT(MXSHEL), ALONE(MXSHEL), ALLSHR, ALLSEG, NOBLK
      DIMENSION IPRVOR(MXSHEL), IPRVPR(MXSHEL),
     *          JORBSH(MXSHEL,MXAOVC), KORBSH(MXSHEL,MXAOVC)
      DIMENSION JSTRSH(MXSHEL,MXAOVC,2), NPRIMS(MXSHEL,MXAOVC,2),
     &          NCONTS(MXSHEL,MXAOVC,2), IORBSH(MXSHEL,MXAOVC)
#include "gnrinf.h"
#include "shells.h"
#include "blocks.h"
#include "symmet.h"
#include "nuclei.h"
      BIGVEC = .FALSE.
      SEGMEN = .FALSE.
      MXAOSH = 0
      NLRGBL = 0
C
C     First assign each shell to a vector
C
      DO 100 I = 1, KMAX
         ORBOUT(I) = .FALSE.
 100  CONTINUE
C
C     Determine ALONE() and ALLSEG
C
      ALLSEG = .TRUE.
      ALLSHR = .TRUE.
      DO 110 I = 1, KMAX
         ALLSEG = ALLSEG .AND. SEGM(I)
         ALLSHR = ALLSHR .AND. .NOT.SHARE(I)
 110  CONTINUE
      DO 130 I = 1, KMAX
         ALONE(I) = .NOT.ALLSHR .AND. .NOT.SHARE(I)
 130  CONTINUE
C
      MAXSHL = 0
      NAOOUT = 0
      NSYMBL = 0
 1000 CONTINUE
         DO 200 IA = 1, KMAX
            IF (.NOT.ORBOUT(IA)) THEN
               NHKTA  = NHKT(IA)
               KHKTA  = KHKT(IA)
               KCKTA  = KCKT(IA)
               NCENTA = NCENT(IA)
               NBCHA  = NBCH(IA)
               NUCA   = NUCO(IA)
               MAXSHL = MAXSHL + 1
               NUCSH  = 0
               LCLASH(MAXSHL) = LCLASS(IA)
               IF (RELCAL .AND. LCLASS(IA).EQ.1) NLRGBL = NLRGBL + 1
               NHKTSH(MAXSHL) = NHKTA
               KHKTSH(MAXSHL) = KHKTA
               KCKTSH(MAXSHL) = KCKTA
               ISTBSH(MAXSHL) = ISTBAO(IA)
               SEGMSH(MAXSHL) = SEGM(IA)
               SPHRSH(MAXSHL) = SPHR(IA)
C              Basis-set identifier (WK/UniKA/04-11-2002).
               MBSISH(MAXSHL) = MBSID(IA)
C
C              Make offsets for symmetry-dependent centers
C
               NDEG = NUCDEG(NCENTA)
               DO IDEG = 1,NDEG
                  NSYMBL = NSYMBL + 1
                  ISYMBL(MAXSHL,IDEG) = NSYMBL
               ENDDO
               IF (BIGVEC) THEN
                  IF (NCENTA .EQ. IATOM) THEN
                     NCNTSH(MAXSHL) = NCENTA
                  ELSE
                     NCNTSH(MAXSHL) = 0
                  END IF
                  IF (MXAOSH .EQ. 0) THEN
                      IF (NHKTA .EQ. 1) MAXVC = 10
                      IF (NHKTA .EQ. 2) MAXVC = 5
                      IF (NHKTA .GE. 3) MAXVC = 3
                  ELSE
                      MAXVC = MXAOSH
                  END IF
               ELSE
                  NCNTSH(MAXSHL) = NCENTA
                  CENTSH(MAXSHL,1) = CENT(IA,1,1)
                  CENTSH(MAXSHL,2) = CENT(IA,2,1)
                  CENTSH(MAXSHL,3) = CENT(IA,3,1)
                  IF (MXAOSH .EQ. 0) THEN
                      MAXVC = NUCA
                  ELSE
                      MAXVC = MXAOSH
                  END IF
               END IF
               IORB = 1
               NSET = 0
               LBCH = 0
               DO 210 IB = 1, KMAX
                  IF (ORBOUT(IB))                              GO TO 210
                  IF (NHKT(IB) .NE. NHKTA)                     GO TO 210
                  IF (KHKT(IB) .NE. KHKTA)                     GO TO 210
                  IF (NCENT(IB) .NE. NCENTA) THEN
                     IF (.NOT.BIGVEC)                          GO TO 210
                     IF (BIGVEC .AND. ALONE(IA))               GO TO 210
                     IF (BIGVEC .AND. ALONE(IB))               GO TO 210
                  END IF
                  IF (NCENTA.EQ.IATOM.AND.NCENT(IB).NE.IATOM)  GO TO 210
                  IF (NCENTA.NE.IATOM.AND.NCENT(IB).EQ.IATOM)  GO TO 210
                  IF (RELCAL.AND.(LCLASS(IB).NE.LCLASS(IA)))   GO TO 210
C
C                 Correct orbital type
C
                  NUCB = NUCO(IB)
                  IF (SEGMEN) THEN
                     IF (NUCSH + NUCB .LE. MAXVC) THEN
                        ORBOUT(IB) = .TRUE.
                        IPTSHL(IB) = MAXSHL
                        JORBSH(MAXSHL,IORB) = IB
                        IORB       = IORB   + 1
                        NAOOUT     = NAOOUT + 1
                        NUCSH      = NUCSH  + NUCB
                        NSET       = 1
                     END IF
                  ELSE
                     NBCHB = NBCH(IB)
                     IF ((NBCHB .EQ. LBCH) .OR.
     *                    (NUCSH + NUCB .LE. MAXVC)) THEN
                        ORBOUT(IB) = .TRUE.
                        IPTSHL(IB) = MAXSHL
                        JORBSH(MAXSHL,IORB) = IB
                        IORB       = IORB   + 1
                        NAOOUT     = NAOOUT + 1
                        IF (NBCHB .NE. LBCH) THEN
                           NSET  = NSET + 1
                           LBCH  = NBCHB
                           NUCSH = NUCSH  + NUCB
                           KORBSH(MAXSHL,NSET) = IB
                        END IF
                     ELSE IF (NUCSH + NUCB .GT. MAXVC .AND. NOBLK
     &                          .AND. NHKT(IB) .EQ. NHKTA .AND.
     &                            NCENT(IB) .EQ. NCENTA .AND.
     &                            IA .NE. IB) THEN
                        WRITE (LUPRI,'(/1X,A)')
     &                      'PAOVEC: Blocking of basis sets not '//
     &                      'permitted in derivative calculations'
                        CALL QUIT('Blocked basis sets not allowed for'//
     &                       ' this type of calculation')
                     END IF
                  END IF
  210          CONTINUE
               NUCOSH(MAXSHL) = NUCSH
               NSETSH(MAXSHL,1) = NSET
               IF (NAOOUT .LT. KMAX) THEN
                   GO TO 1000
               ELSE
                   GO TO 1100
               END IF
            END IF
  200    CONTINUE
      IF (NAOOUT .LT. KMAX) GO TO 1000
 1100 CONTINUE
C
      NSMLBL = MAXSHL - NLRGBL
C
C     All orbitals have now been assigned to a shell
C
C
C     Determine NORBSH()
C
      DO 300 I = 1, MAXSHL
         NORBSH(I) = 0
         DO 310 J = 1, KMAX
            IF (IPTSHL(J) .EQ. I) NORBSH(I) = NORBSH(I) + 1
  310    CONTINUE
  300 CONTINUE
C
C     Determine and NSTRSH()
C
C     NSTRSH(1) = 0
C     DO 500 IA = 1,MAXSHL - 1
C        NSTRSH(IA + 1)   = NSTRSH(IA)   + KHKTSH(IA)*NORBSH(IA)
C 500 CONTINUE
C
C     Determine IORBSH(MAXSHL,NORBSH)
C
      IORB = 0
      IPRVOR(1) = 0
      DO 800 I = 1, KMAX - 1
         IORB = IORB + KHKT(I)*MULT(ISTBAO(I))
         IPRVOR(I + 1) = IORB
  800 CONTINUE
      DO 810 I = 1, MAXSHL
         DO 820 J = 1, NORBSH(I)
            IORBSH(I,J) = IPRVOR(JORBSH(I,J))
  820    CONTINUE
  810 CONTINUE
C
C     IORBSB()
C
      IORBB = 0
      DO 900 I = 1, MAXSHL
         DO 910 J = 1, NORBSH(I)
            IORBS = IORBSH(I,J)
            DO 920 ICMP = 1, KHKTSH(I)
               IORBSB(IORBS) = IORBB
               IORBB = IORBB + 1
               IORBS = IORBS + 1
  920       CONTINUE
  910    CONTINUE
  900 CONTINUE
      NTOTAL = IORBB + 1
      CALL PAOSET(JSTRSH,NPRIMS,NCONTS,KORBSH)
C
C     Print Section
C
      IF (IPRINT .GE. 5) THEN
         CALL HEADER('Output from PAOVEC',-1)
         WRITE (LUPRI,'(1X,A,L5)') ' BIGVEC ', BIGVEC
         WRITE (LUPRI,'(1X,A,I5)') ' MAXSHL ', MAXSHL
         WRITE (LUPRI,'(1X,A,(5I5))') ' NHKTSH ',(NHKTSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' KHKTSH ',(KHKTSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5L5))') ' SPHRSH ',(SPHRSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' KCKTSH ',(KCKTSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' NUCOSH ',(NUCOSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))') ' NORBSH ',(NORBSH(I),I=1,MAXSHL)
C        WRITE (LUPRI,'(1X,A,(5I5))') ' NSTRSH ',(NSTRSH(I),I=1,MAXSHL)
         WRITE (LUPRI,'(1X,A,(5I5))')' IORBSB ',(IORBSB(I-1),I=1,NTOTAL)
         IF (BIGVEC) THEN
            WRITE (LUPRI,'(1X,A,(10I5))')
     *            ' NCNTSH ',(NCNTSH(I),I=1,MAXSHL)
         END IF
         DO 2000 I = 1, MAXSHL
            WRITE (LUPRI,'(1X,A,I5)')'JSTRSH (gen. con.) for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)')(JSTRSH(I,J,1),J=1,NORBSH(I))
            WRITE (LUPRI,'(1X,A,I5)')'JSTRSH (seg. con.) for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)')(JSTRSH(I,J,2),J=1,NORBSH(I))
 2000    CONTINUE
         DO 2001 I = 1, MAXSHL
            WRITE (LUPRI,'(1X,A,I5)')'NPRIMS (gen. con.) for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)')(NPRIMS(I,J,1),J=1,NORBSH(I))
            WRITE (LUPRI,'(1X,A,I5)')'NPRIMS (seg. con.) for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)')(NPRIMS(I,J,2),J=1,NORBSH(I))
 2001    CONTINUE
         DO 2004 I = 1, MAXSHL
            WRITE (LUPRI,'(1X,A,I5)')'NCONTS (gen. con.) for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)')(NCONTS(I,J,1),J=1,NORBSH(I))
            WRITE (LUPRI,'(1X,A,I5)')'NCONTS (seg. con.) for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)')(NCONTS(I,J,2),J=1,NORBSH(I))
 2004    CONTINUE
         DO 2002 I = 1, MAXSHL
            WRITE (LUPRI,'(1X,A,I5)')'CENTSH for ISHELL =',I
            WRITE (LUPRI,'(1X,3F20.10)')
     *             CENTSH(I,1),CENTSH(I,2),CENTSH(I,3)
 2002    CONTINUE
         DO 2003 I = 1, MAXSHL
            WRITE (LUPRI,'(1X,A,I5)')'IORBSH for ISHELL =',I
            WRITE (LUPRI,'(1X,15I5)')(IORBSH(I,J),J=1,NORBSH(I))
 2003    CONTINUE
      END IF
      RETURN
      END
C  /* Deck cordif */
      SUBROUTINE CORDIF(NORBA,NORBB,THRESH,D0X,D0Y,D0Z,
     &                  IPRINT,NUCOA,NUCOB,JSTRA,JSTRB)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
      LOGICAL DIFX, DIFY, DIFZ, D0X, D0Y, D0Z
      DIMENSION NUCOA(*), NUCOB(*), JSTRA(*), JSTRB(*)
#include "primit.h"
C
      DIFX = .FALSE.
      DIFY = .FALSE.
      DIFZ = .FALSE.
C
C     A - A
C
      DO 100 I = 1, NORBA
      DO 100 J = 1, NUCOA(I)
         IJ  = JSTRA(I) + J
         CRX = PRICRX(IJ)
         CRY = PRICRY(IJ)
         CRZ = PRICRZ(IJ)
         DO 200 K = 1, I
         DO 200 L = 1, NUCOA(K)
            KL   = JSTRA(K) + L
            DIFX = DIFX .OR. ABS(PRICRX(KL)-CRX) .GT. THRESH
            DIFY = DIFY .OR. ABS(PRICRY(KL)-CRY) .GT. THRESH
            DIFZ = DIFZ .OR. ABS(PRICRZ(KL)-CRZ) .GT. THRESH
  200    CONTINUE
  100 CONTINUE
C
C     B - B
C
      IF (.NOT.(DIFX .AND. DIFY .AND. DIFZ)) THEN
         DO 300 I = 1, NORBB
         DO 300 J = 1, NUCOB(I)
            IJ  = JSTRB(I) + J
            CRX = PRICRX(IJ)
            CRY = PRICRY(IJ)
            CRZ = PRICRZ(IJ)
            DO 400 K = 1, I
            DO 400 L = 1, NUCOB(K)
               KL   = JSTRB(K) + L
               DIFX = DIFX .OR. ABS(PRICRX(KL)-CRX) .GT. THRESH
               DIFY = DIFY .OR. ABS(PRICRY(KL)-CRY) .GT. THRESH
               DIFZ = DIFZ .OR. ABS(PRICRZ(KL)-CRZ) .GT. THRESH
  400       CONTINUE
  300    CONTINUE
      END IF
C
C     A - B
C
      IF (.NOT.(DIFX .AND. DIFY .AND. DIFZ)) THEN
         DO 500 I = 1, NORBA
         DO 500 J = 1, NUCOA(I)
            IJ  = JSTRA(I) + J
            CRX = PRICRX(IJ)
            CRY = PRICRY(IJ)
            CRZ = PRICRZ(IJ)
            DO 600 K = 1, NORBB
            DO 600 L = 1, NUCOB(K)
               KL   = JSTRB(K) + L
               DIFX = DIFX .OR. ABS(PRICRX(KL)-CRX) .GT. THRESH
               DIFY = DIFY .OR. ABS(PRICRY(KL)-CRY) .GT. THRESH
               DIFZ = DIFZ .OR. ABS(PRICRZ(KL)-CRZ) .GT. THRESH
  600       CONTINUE
  500    CONTINUE
      END IF
      D0X = .NOT.DIFX
      D0Y = .NOT.DIFY
      D0Z = .NOT.DIFZ
C
      IF (IPRINT .LT. 05) RETURN
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      CALL HEADER('SUBROUTINE CORDIF',-1)
      WRITE (LUPRI, 1010) NORBA, NORBB
      WRITE (LUPRI, 1020) (NUCOA(I), I = 1, NORBA)
      WRITE (LUPRI, 1030) (NUCOB(I), I = 1, NORBB)
      WRITE (LUPRI, 1040) (JSTRA(I), I = 1, NORBA)
      WRITE (LUPRI, 1050) (JSTRB(I), I = 1, NORBB)
      WRITE (LUPRI, 1060) D0X, D0Y, D0Z
 1010 FORMAT(  '  NORB     ',2I7)
 1020 FORMAT(  '  NUCOA:   ',15I7)
 1030 FORMAT(  '  NUCOB:   ',15I7)
 1040 FORMAT(  '  JSTRA:   ',15I7)
 1050 FORMAT(  '  JSTRB:   ',15I7)
 1060 FORMAT(  '  D0X/Y/Z: ',3L5)
      RETURN
      END
C  /* Deck nintso */
      SUBROUTINE NINTSO(MULE,LONDON,SPNORB,UNDIFF,SOFOCK,DISTRI,SQ12EL,
     &                  IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      LOGICAL LONDON, SPNORB, UNDIFF, SOFOCK, DISTRI, SQ12EL
#include "twocom.h"
#include "symmet.h"
#include "dorps.h"
#include "ibtfun.h"
C
      IF (UNDIFF .OR. SOFOCK .OR. DISTRI) THEN
         NINTS     = NMBOSO(SQ12EL,0)
         NINTMX    = NINTS
         NINTSR(1) = NINTS
      ELSE
         NINTS  = 0
         NINTMX = 0
         IRPXYZ = IBTXOR(ISYMAX(1,1),IBTXOR(ISYMAX(2,1),ISYMAX(3,1)))
         DO 100 IREPX = 1, NOPREP
            IREPE  = IPTREP(IREPX,1)
            NINTE  = NMBOSO(SQ12EL,IREPE)
            NINTMX = MAX(NINTMX,NINTE)
            NINTSR(IREPX) = NINTE
            NTYPE = 0
            DO 200 ICOOR = 1, 3
               IF (SPNORB) THEN
                 IF(IREPE .EQ. IBTXOR(ISYMAX(ICOOR,1),IRPXYZ)) THEN
                    NTYPE = NTYPE + 1
                 END IF
               ELSE IF (LONDON) THEN
                 IF(IREPE .EQ. ISYMAX(ICOOR,2)) NTYPE = NTYPE + 2
               ELSE
                 IF(IBTAND(MULE,IBTXOR(IREPE,ISYMAX(ICOOR,1))).EQ.0)THEN
                    NTYPE = NTYPE + 1
                 END IF
               END IF
  200       CONTINUE
            NINTS  = NINTS + NTYPE*NINTE
  100    CONTINUE
      END IF
      IF (IPRINT .GT. 5) THEN
         CALL TITLER('Output from NINTSO','*',103)
         WRITE (LUPRI,'(2X,A, I5)') ' NINTS  ', NINTS
         WRITE (LUPRI,'(2X,A,8I5)') ' NINTSR ', (NINTSR(I),I=1,NOPREP)
         WRITE (LUPRI,'(2X,A, I5)') ' NINTMX ', NINTMX
      END IF
      RETURN
      END
C  /* Deck nmboso */
      FUNCTION NMBOSO(SQ12EL,IREPE)
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      LOGICAL SQ12EL, DCMPAC, DCABAB, NOTEST
#include "twocom.h"
#include "symmet.h"
#include "ibtfun.h"
C
      NOTEST = .NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)
      IF (NOTEST) THEN
        NMBOSO = 0
        DO 100 NA = 1,KHKTA
          ISA = ISYMAO(NHKTA,NA)
          DO 110 IREPA = 0, MAXREP
          IF (IBTAND(MULA,IBTXOR(IREPA,ISA)) .EQ. 0) THEN
            IRPAE = IBTXOR(IREPA,IREPE)
            DO 200 NB = 1,KHKTB
              ISB = ISYMAO(NHKTB,NB)
              DO 210 IREPB = 0, MAXREP
              IF (IBTAND(MULB,IBTXOR(IREPB,ISB)) .EQ. 0) THEN
                IRPABE = IBTXOR(IREPB,IRPAE)
                DO 300 NC = 1,KHKTC
                  ISC = ISYMAO(NHKTC,NC)
                  DO 410 IREPC = 0, MAXREP
                  IF (IBTAND(MULC,IBTXOR(IREPC,ISC)) .EQ. 0) THEN
                    IREPD=IBTXOR(IREPC,IRPABE)
                    DO 500 ND = 1,KHKTD
                      ISD = ISYMAO(NHKTD,ND)
                      IF (IBTAND(MULD,IBTXOR(IREPD,ISD)) .EQ. 0) THEN
                         NMBOSO = NMBOSO + 1
                      END IF
  500               CONTINUE
                  END IF
  410             CONTINUE
  300           CONTINUE
              END IF
  210         CONTINUE
  200       CONTINUE
          END IF
  110     CONTINUE
  100   CONTINUE
      ELSE
        NMBOSO = 0
        DO 600 NA = 1,KHKTA
          ISA = ISYMAO(NHKTA,NA)
          KHKTBB = KHKTB
          IF (DIAGAB) KHKTBB = NA
          DO 610 NB = 1,KHKTBB
            ISB = ISYMAO(NHKTB,NB)
            KHKTCC = KHKTC
            IF (.NOT.SQ12EL .AND. SHABAB) KHKTCC = NA
            DO 620 NC = 1,KHKTCC
              ISC = ISYMAO(NHKTC,NC)
              DCMPAC = SHABAB .AND. NA .EQ. NC
              KHKTDD = KHKTD
              IF (DIAGCD) KHKTDD = NC
              IF (.NOT.SQ12EL .AND. DCMPAC) KHKTDD = NB
              DO 630 ND = 1,KHKTDD
                ISD = ISYMAO(NHKTD,ND)
                DCABAB = .NOT.SQ12EL .AND. DCMPAC .AND. NB .EQ. ND
C
                DO 700 IREPA = 0, MAXREP
                IF (IBTAND(MULA,IBTXOR(IREPA,ISA)) .EQ. 0) THEN
                  DO 710 IREPB = 0, MAXREP
                  IF (IBTAND(MULB,IBTXOR(IREPB,ISB)) .EQ. 0) THEN
                   DO 720 IREPC = 0, MAXREP
                   IF (IBTAND(MULC,IBTXOR(IREPC,ISC)) .EQ. 0) THEN
                   IREPD=IBTXOR(IBTXOR(IBTXOR(IREPA,IREPB),IREPC),IREPE)
                    IF (IBTAND(MULD,IBTXOR(IREPD,ISD)) .EQ. 0) THEN
                     IF (.NOT.(DCABAB .AND. (IREPA .LT. IREPC .OR.
     &                  (IREPA.EQ.IREPC .AND. IREPB.LT.IREPD)))) THEN
                         NMBOSO = NMBOSO + 1
                     END IF
                    END IF
                   END IF
  720              CONTINUE
                  END IF
  710             CONTINUE
                END IF
  700           CONTINUE
C
  630         CONTINUE
  620       CONTINUE
  610     CONTINUE
  600   CONTINUE
      END IF
      RETURN
      END
C  /* Deck intder */
      SUBROUTINE INTDER(AOINT,NCCINT,NINTYP,WORK,LWORK,COEFAB,COEFCD,
     &                  COORAB,COORCD,EXPAB,EXPCD,FACAB,FACCD,CONTA,
     &                  CONTB,CONTC,CONTD,IPRINT,NOINT,NOCONT,SPNORB,
     &                  DIA2SO,ZFS2EL,MAXDER,THRESH,JMAX0,ISMXYZ,ONECEN,
     &                  NPCOA,NPCOB,NPCOC,NPCOD,NUCSAB,NUCSCD,SIGNT,
     &                  IODDHR,INDHSQ,INDHER,LMNVLS,NPNTAB,NPNTCD,
     &                  NREDAB,NREDCD)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxaqn.h"
      LOGICAL NOINT, NOCONT, SPNORB, DIA2SO, ONECEN, ZFS2EL
      DIMENSION AOINT(NCCINT,NINTYP), INDHSQ(*), IODDHR(*),
     &          WORK(LWORK), COEFAB(*), COEFCD(*), COORAB(*), COORCD(*),
     &          EXPAB(*), EXPCD(*), FACAB(*), FACCD(*), LMNVLS(*),
     &          CONTA(*), CONTB(*), CONTC(*), CONTD(*),
     &          NPCOA(*), NPCOB(*), NPCOC(*), NPCOD(*),
     &          NUCSAB(*), NUCSCD(*), SIGNT(*), INDHER(*),
     &          NPNTAB(*), NPNTCD(*), NREDAB(*), NREDCD(*)
#include "twocom.h"
#include "twosta.h"
#include "drw2el.h"
#include "r12int.h"
      CALL QENTER('INTDER')
#include "memint.h"
      IF (IPRINT .GT. 5) CALL TITLER('Output from INTDER','*',103)
C
      JMAX   = JMAX0 + MAXDER
C
C     For two-electron integrals for relativistic Direct Perturbation
C     Theory (DPT), JMAX must be increased by 2 (WK/UniKA/14-11-2002).
C
      IF (FINDPT .OR. DPTINT) JMAX = JMAX + 2
      IF (BPH2OO) JMAX = JMAX + 4
C
C     For [r12,T1+T2] integrals, JMAX must be increased by 1. For
C     r12 integrals, JMAX must be 1 for (ss|ss) integrals, but is
C     not increased otherwise (WK/UniKA/14-11-2002).
C
      IF (U12INT) JMAX = JMAX + 1
      IF (INTGAC .EQ. 5) JMAX = JMAX + 1
      IF (R12INT) JMAX = MAX(JMAX,1)
C
      NUABCD = NUCAB*NUCCD
      NRTUV  = (JMAX + 1)**3
      NTUV   = (JMAX + 1)*(JMAX + 2)*(JMAX + 3)/6
C
C     Allocate work space for Hermite integrals
C     More space is needed for R12 method (WK/UniKA/14-11-2002).
C
      IF (R12EIN .OR. R12INT .OR. U12INT .OR. BPH2OO) THEN
         LHRINT = 2*NUABCD*NTUV
      ELSE
         LHRINT = NUABCD*NTUV
      END IF
      CALL MEMGET('REAL',KHRINT,LHRINT,WORK,KFREE,LFREE)
      MWHRIN = MAX(MWHRIN,LHRINT)
      LWTOT  = LWTOT + KFREE
      MWTOT  = MAX(MWTOT,LWTOT)
      CALL INTDR1(AOINT,NCCINT,NINTYP,WORK(KHRINT),INDHER,NTUV,COEFAB,
     &            COEFCD,COORAB,COORCD,EXPAB,EXPCD,FACAB,FACCD,CONTA,
     &            CONTB,CONTC,CONTD,WORK(KFREE),LFREE,NUABCD,IPRINT,
     &            NOINT,NOCONT,SPNORB,DIA2SO,ZFS2EL,MAXDER,THRESH,JMAX,
     &            ISMXYZ,ONECEN,NPCOA,NPCOB,NPCOC,NPCOD,NUCSAB,NUCSCD,
     &            SIGNT,INDHSQ,IODDHR,LMNVLS,NPNTAB,NPNTCD,
     &            NREDAB,NREDCD,V12INT,R12INT,U12INT,
     &            R12EIN,U21INT,INTGAC,IADV12,IADR12,IADU12)
      LWTOT  = LWTOT - KFREE
      CALL MEMREL('INTDER',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('INTDER')
      RETURN
      END
C  /* Deck intdr1 */
      SUBROUTINE INTDR1(AOINT,NCCINT,NINTYP,HERINT,INDHER,NTUV,COEFAB,
     &                  COEFCD,COORAB,COORCD,EXPAB,EXPCD,FACAB,FACCD,
     &                  CONTA,CONTB,CONTC,CONTD,WORK,LWORK,NUABCD,
     &                  IPRINT,NOINT,NOCONT,SPNORB,DIA2SO,ZFS2EL,MAXDER,
     &                  THRESH,JMAX,ISMXYZ,ONECEN,NPCOA,NPCOB,NPCOC,
     &                  NPCOD,NUCSAB,NUCSCD,SIGNT,INDHSQ,IODDHR,LMNVLS,
     &                  NPNTAB,NPNTCD,NREDAB,NREDCD,V12INT,R12INT,
     &                  U12INT,R12EIN,U21INT,INTGAC,
     &                  IADV12,IADR12,IADU12)
C
C     tuh
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      LOGICAL R12INT, U12INT, V12INT, R12EIN
      LOGICAL NOINT,NOCONT,SPNORB,DIA2SO,PQSYM,ONECEN, ZFS2EL, U21INT
C     3-dimensional HERINT for r12 integrals (WK/UniKA/15-11-2002).
      DIMENSION HERINT(NUABCD,NTUV,2), INDHER(*),
     &          AOINT(NCCINT,NINTYP), WORK(LWORK),
     &          COEFAB(*), COEFCD(*), COORAB(*), COORCD(*),
     &          EXPAB(*), EXPCD(*), FACAB(*), FACCD(*),
     &          CONTA(*), CONTB(*), CONTC(*), CONTD(*),
     &          NPCOA(NSETA,2), NPCOB(NSETB,2),
     &          NPCOC(NSETC,2), NPCOD(NSETD,2),
     &          NUCSAB(*), NUCSCD(*), SIGNT(*),
     &          INDHSQ(*), IODDHR(NRTOP,0:7),
     &          LMNVLS(KCKMAX,5,4),
     &          NPNTAB(*), NPNTCD(*), NREDAB(*), NREDCD(*)
#include "hertop.h"
#include "ccom.h"
#include "twocom.h"
#include "twosta.h"
#include "twoao.h"
#include "crossd.h"
#include "subdir.h"
#include "drw2el.h"
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from INTDR1','*',103)
C
C     ***************************************
C     ********** Hermite Integrals **********
C     ***************************************
C
      NOINT = .FALSE.
      IF (.NOT.DPATH1) THEN
         CALL HR2DRV(HERINT,INDHER,COORCD,COORAB,EXPCD,EXPAB,FACCD,
     &               FACAB,NTUV,WORK,LWORK,JMAX,MAXDER,NUABCD,IPQ0X,
     &               IPQ0Y,IPQ0Z,NOINT,ONECEN,NUCC,NUCD,NUCCD,NUCA,NUCB,
     &               NUCAB,THRESH,IPRINT,SIGNT,IODDHR(1,ISMXYZ))
      ELSE
         CALL HR2DRV(HERINT,INDHER,COORAB,COORCD,EXPAB,EXPCD,FACAB,
     &               FACCD,NTUV,WORK,LWORK,JMAX,MAXDER,NUABCD,IPQ0X,
     &               IPQ0Y,IPQ0Z,NOINT,ONECEN,NUCA,NUCB,NUCAB,NUCC,NUCD,
     &               NUCCD,THRESH,IPRINT,SIGNT,IODDHR(1,ISMXYZ))
      END IF
C
C     *****************************************
C     ********** Cartesian Integrals **********
C     *****************************************
C
      IHRSYM = ISMXYZ
      KCKMX  = KCKMAX
      NOCNT  = NOCONT
      IDERIV = MAXDER
      IF (DIA2SO .OR. BPH2OO) IDERIV = 1
      IF (FINDPT.OR.DPTINT) IDERIV = 2
      JM1234 = JMAX
      NU1234 = NUABCD
      NC1234 = NCCINT
      SPIORB = SPNORB
      IF (.NOT.NOINT) THEN
         IF (NOCONT) THEN
            NRBA  = NUCA
            NRBB  = NUCB
            NRBC  = NUCC
            NRBD  = NUCD
         ELSE
            NRBA  = NORBA
            NRBB  = NORBB
            NRBC  = NORBC
            NRBD  = NORBD
         END IF
         CALL DZERO(AOINT,NCCINT*NINTYP)
         IF (DPATH1 .AND. .NOT. (U12INT .AND. .NOT. U21INT)) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            IPATH  = 1
            PATH1  = .TRUE.
            NHKT1  = NHKTA
            NHKT2  = NHKTB
            NHKT3  = NHKTC
            NHKT4  = NHKTD
            KCKT1  = KCKTA
            KCKT2  = KCKTB
            KCKT3  = KCKTC
            KCKT4  = KCKTD
            KCKT12 = KCKTAB
            KCKT34 = KCKTCD
            KHKT1  = KHKTA
            KHKT2  = KHKTB
            KHKT3  = KHKTC
            KHKT4  = KHKTD
            KHKT12 = KHKTAB
            KHKT34 = KHKTCD
            NORB1  = NRBA
            NORB2  = NRBB
            NORB3  = NRBC
            NORB4  = NRBD
            NORR1  = ISUM(NSETA,NPCOA(1,2),1)
            NORR2  = ISUM(NSETB,NPCOB(1,2),1)
            NORR3  = ISUM(NSETC,NPCOC(1,2),1)
            NORR4  = ISUM(NSETD,NPCOD(1,2),1)
            NORB12 = NORBAB
            NORB34 = NORBCD
            NUC1   = NUCA
            NUC2   = NUCB
            NUC3   = NUCC
            NUC4   = NUCD
            NUC12  = NUCAB
            NUC34  = NUCCD
            NUCR1  = ISUM(NSETA,NPCOA(1,1),1)
            NUCR2  = ISUM(NSETB,NPCOB(1,1),1)
            NUCR3  = ISUM(NSETC,NPCOC(1,1),1)
            NUCR4  = ISUM(NSETD,NPCOD(1,1),1)
            MXUC12 = MXUCAB
            MXUC34 = MXUCCD
            NSET1  = NSETA
            NSET2  = NSETB
            NSET3  = NSETC
            NSET4  = NSETD
            MAX12  = MAXAB
            MAX34  = MAXCD
            JMAX1  = JMAXA
            JMAX2  = JMAXB
            JMAX3  = JMAXC
            JMAX4  = JMAXD
            TPRI12 = TPRIAB
            TPRI34 = TPRICD
            TCON12 = TCONAB
            TCON34 = TCONCD
            DIAG12 = DIAGAB
            DIAG34 = DIAGCD
            DIAC12 = DIACAB
            DIAC34 = DIACCD
            I120X  = IAB0X
            I120Y  = IAB0Y
            I120Z  = IAB0Z
            I340X  = ICD0X
            I340Y  = ICD0Y
            I340Z  = ICD0Z
            NO1234 = NOABCD
            GEN12  = GENAB
            GEN34  = GENCD
            RPRI12 = RPRIAB
            RPRI34 = RPRICD
            RCNT12 = RCNTAB
            RCNT34 = RCNTCD
            SPHR1  = SPHRA
            SPHR2  = SPHRB
            SPHR3  = SPHRC
            SPHR4  = SPHRD
            SPHR12 = SPHRAB
            SPHR34 = SPHRCD
            DC10   = DC101
            DC1E   = DC1E1
            DC2H   = DC2H1
            DC2E   = DC2E1
            CROSS  = CROSS1
            MAX34D = MAX34 + IDERIV
C
C           Special case: [T1+T2,r12] integrals (WK/UniKA/14-11-2002).
            IF (U12INT) MAX34D = MAX34 + 1
            IF (INTGAC .EQ. 5) MAX34D = MAX34D + 1
C
            NTUV34 = (MAX34D + 1)*(MAX34D + 2)*(MAX34D + 3)/6
            MAX34D = MAX34 + IDERIV - 1
            KTUV34 = (MAX34D + 1)*(MAX34D + 2)*(MAX34D + 3)/6
            LTUV34 = (MAX34D    )*(MAX34D + 1)*(MAX34D + 2)/6
            IF (DIA2SO) KTUV34 = NTUV34
            NCCPP  = NORB12*NUC34
            CALL CCDRIV(AOINT,NCCINT,HERINT,INDHER,COEFAB,COEFCD,CONTA,
     &                  CONTB,CONTC,CONTD,WORK,LWORK,IPRINT,NPCOA,NPCOB,
     &                  NPCOC,NPCOD,NUCSAB,NUCSCD,INDHSQ,
     &                  IODDHR(1,ISMXYZ),LMNVLS(1,1,1),LMNVLS(1,1,3),
     &                  NPNTAB,NPNTCD,NREDAB,NREDCD,NTUV,NUABCD,
     &                  V12INT,R12INT,U12INT,R12EIN,
     &                  IADV12,IADR12,IADU12)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TPATH1 = TPATH1 + TIMEND - TIMSTR
            END IF
         END IF
         IF (DPATH2) THEN
            IF (TKTIME) CALL GETTIM(TIMSTR,DUMTIM)
            IPATH  = 2
            PATH1  = .FALSE.
            NHKT1  = NHKTC
            NHKT2  = NHKTD
            NHKT3  = NHKTA
            NHKT4  = NHKTB
            KCKT1  = KCKTC
            KCKT2  = KCKTD
            KCKT3  = KCKTA
            KCKT4  = KCKTB
            KCKT12 = KCKTCD
            KCKT34 = KCKTAB
            KHKT1  = KHKTC
            KHKT2  = KHKTD
            KHKT3  = KHKTA
            KHKT4  = KHKTB
            KHKT12 = KHKTCD
            KHKT34 = KHKTAB
            NORB1  = NRBC
            NORB2  = NRBD
            NORB3  = NRBA
            NORB4  = NRBB
            NORR1  = ISUM(NSETC,NPCOC(1,2),1)
            NORR2  = ISUM(NSETD,NPCOD(1,2),1)
            NORR3  = ISUM(NSETA,NPCOA(1,2),1)
            NORR4  = ISUM(NSETB,NPCOB(1,2),1)
            NORB12 = NORBCD
            NORB34 = NORBAB
            NUC1   = NUCC
            NUC2   = NUCD
            NUC3   = NUCA
            NUC4   = NUCB
            NUCR1  = ISUM(NSETC,NPCOC(1,1),1)
            NUCR2  = ISUM(NSETD,NPCOD(1,1),1)
            NUCR3  = ISUM(NSETA,NPCOA(1,1),1)
            NUCR4  = ISUM(NSETB,NPCOB(1,1),1)
            NUC12  = NUCCD
            NUC34  = NUCAB
            MXUC12 = MXUCCD
            MXUC34 = MXUCAB
            NSET1  = NSETC
            NSET2  = NSETD
            NSET3  = NSETA
            NSET4  = NSETB
            MAX12  = MAXCD
            MAX34  = MAXAB
            JMAX1  = JMAXC
            JMAX2  = JMAXD
            JMAX3  = JMAXA
            JMAX4  = JMAXB
            TPRI12 = TPRICD
            TPRI34 = TPRIAB
            TCON12 = TCONCD
            TCON34 = TCONAB
            DIAG12 = DIAGCD
            DIAG34 = DIAGAB
            DIAC12 = DIACCD
            DIAC34 = DIACAB
            I120X  = ICD0X
            I120Y  = ICD0Y
            I120Z  = ICD0Z
            I340X  = IAB0X
            I340Y  = IAB0Y
            I340Z  = IAB0Z
            NO1234 = NOABCD
            GEN12  = GENCD
            GEN34  = GENAB
            RPRI12 = RPRICD
            RPRI34 = RPRIAB
            RCNT12 = RCNTCD
            RCNT34 = RCNTAB
            SPHR1  = SPHRC
            SPHR2  = SPHRD
            SPHR3  = SPHRA
            SPHR4  = SPHRB
            SPHR12 = SPHRCD
            SPHR34 = SPHRAB
            DC10   = DC102
            DC1E   = DC1E2
            DC2H   = DC2H2
            DC2E   = DC2E2
            CROSS  = CROSS2
            IF (DPATH1) THEN
               CALL HERSWP(HERINT(1,1,1),NTUV,NUABCD,
     &                     WORK,LWORK,JMAX,NUCAB,NUCCD,IPRINT)
C
C              Also swap Hermite integrals over r12 (WK/UniKA/14-11-2002).
               IF (R12INT .OR. U12INT)
     &         CALL HERSWP(HERINT(1,1,2),NTUV,NUABCD,
     &                     WORK,LWORK,JMAX,NUCAB,NUCCD,IPRINT)
            END IF
            MAX34D = MAX34 + IDERIV
C
C           Special case: [T1+T2,r12] integrals (WK/UniKA/14-11-2002).
            IF (U12INT) MAX34D = MAX34 + 1
            IF (INTGAC .EQ. 5) MAX34D = MAX34D + 1
C
            NTUV34 = (MAX34D + 1)*(MAX34D + 2)*(MAX34D + 3)/6
            MAX34D = MAX34 + IDERIV - 1
            KTUV34 = (MAX34D + 1)*(MAX34D + 2)*(MAX34D + 3)/6
            LTUV34 = (MAX34D    )*(MAX34D + 1)*(MAX34D + 2)/6
            NCCPP  = NORB12*NUC34
            CALL CCDRIV(AOINT,NCCINT,HERINT,INDHER,COEFCD,COEFAB,CONTC,
     &                  CONTD,CONTA,CONTB,WORK,LWORK,IPRINT,NPCOC,NPCOD,
     &                  NPCOA,NPCOB,NUCSCD,NUCSAB,INDHSQ,
     &                  IODDHR(1,ISMXYZ),LMNVLS(1,1,3),LMNVLS(1,1,1),
     &                  NPNTCD,NPNTAB,NREDCD,NREDAB,NTUV,NUABCD,
     &                  V12INT,R12INT,U12INT,R12EIN,
     &                  IADV12,IADR12,IADU12)
            IF (TKTIME) THEN
               CALL GETTIM(TIMEND,DUMTIM)
               TPATH2 = TPATH2 + TIMEND - TIMSTR
            END IF
          END IF
      END IF
      RETURN
      END
C  /* Deck ccdriv */
      SUBROUTINE CCDRIV(AOINT,NCCINT,HERINT,INDHER,COEF12,COEF34,CONT1,
     &                  CONT2,CONT3,CONT4,WORK,LWORK,IPRINT,NPCO1,NPCO2,
     &                  NPCO3,NPCO4,NUCS12,NUCS34,INDHSQ,IODDHR,
     &                  LMNV12,LMNV34,NPNT12,NPNT34,NRED12,NRED34,
     &                  NTUV,NUABCD,V12INT,R12INT,U12INT,R12EIN,
     &                  IADV12,IADR12,IADU12)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "iratdef.h"
      DIMENSION HERINT(NUABCD,NTUV,2), INDHER(*), AOINT(NCCINT,*),
     &          WORK(LWORK), COEF12(*), COEF34(*),
     &          CONT1(NORB1*NUC1,2), CONT2(NORB2*NUC2,2),
     &          CONT3(NORB3*NUC3,2), CONT4(NORB4*NUC4,2),
     &          NPCO1(*), NPCO2(*), NPCO3(*), NPCO4(*),
     &          INDHSQ(*), IODDHR(*), NPNT12(*), NPNT34(*),
     &          LMNV12(KCKMX,5,2), LMNV34(KCKMX,5,2),
     &          NRED12(*), NRED34(*), NUCS12(*), NUCS34(*)
#include "twoao.h"
#include "drw2el.h"
#include "crsdir.h"
#include "twosta.h"
      LOGICAL V12INT, R12INT, U12INT, R12EIN
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from CCDRIV','*',103)
C
C     Allocations
C
      NHCINT = 0
      IF (DC10) NHCINT = NHCINT + 1
      IF (DC1E) THEN
         IF (DHCEX1) NHCINT = NHCINT + 1
         IF (DHCEX2) NHCINT = NHCINT + 1
         IF (DHCEY1) NHCINT = NHCINT + 1
         IF (DHCEY2) NHCINT = NHCINT + 1
         IF (DHCEZ1) NHCINT = NHCINT + 1
         IF (DHCEZ2) NHCINT = NHCINT + 1
      END IF
      IF (BPH2OO) NHCINT = NHCINT + 3
      LHCINT = NCCPP*NTUV34*KHKT12*NHCINT
      IF (DPTINT) LHCINT = LHCINT * 2
C
      KODD12 = 1
      KODD34 = KODD12 + (2*KCKT12 + 1)/IRAT
      KHCINT = KODD34 + (2*KCKT34 + 1)/IRAT
      KLAST  = KHCINT + LHCINT
C
C     Allocate work space for r12 integrals (WK/UniKA/15-11-2002).
      KHRINT = KHCINT
      IF (.NOT. R12EIN .AND. (R12INT .OR. U12INT)) THEN
         KHRINT = KLAST
         IF (.NOT. V12INT .AND. .NOT. U12INT) KHRINT = KHCINT
         KLAST = KHRINT + LHCINT
      END IF
      IADUUU = IADU12 + 2 - IPATH
C
      IF (KLAST .GT. LWORK) CALL STOPIT('CCDRIV',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      MWHCIN = MAX(MWHCIN,KLAST - KHCINT)
      LWTOT  = LWTOT + KLAST
      MWTOT  = MAX(MWTOT,LWTOT)
C
      CALL GETODA(WORK(KODD12),WORK(KODD34),LMNV12,LMNV34,IPRINT)
C
      IF (.NOT. GEN12) THEN
C     ... 980826-hjaaj: check if no contraction and no omissions
         NOCNT = .TRUE.
         DO IJ = 1,NORB12
            IF (NUCS12(IJ) .NE. 1) THEN
               NOCNT = .FALSE.
               GO TO 100
            END IF
         END DO
  100    IF (IPRINT .GT. 5) WRITE (LUPRI,*) 'NOCNT for C1DRIV:', NOCNT
      ELSE
         NOCNT = .FALSE.
      END IF
      IF ((R12EIN .OR. V12INT .OR. U12INT ) .AND. .NOT.BPH2OO)
     &CALL C1DRIV(HERINT(1,1,1),INDHER,WORK(KHCINT),COEF12,CONT1(1,2),
     &            CONT2(1,2),WORK(KLAST),LWRK,NPCO1,NPCO2,NUCS12,INDHSQ,
     &            IODDHR,IPRINT,LMNV12,WORK(KODD12),NPNT12,NRED12,
     &            NHCINT)
      IF (((R12INT .OR. U12INT) .AND. .NOT. R12EIN) .OR. BPH2OO)
     &CALL C1DRIV(HERINT(1,1,2),INDHER,WORK(KHRINT),COEF12,CONT1(1,2),
     &            CONT2(1,2),WORK(KLAST),LWRK,NPCO1,NPCO2,NUCS12,INDHSQ,
     &            IODDHR,IPRINT,LMNV12,WORK(KODD12),NPNT12,NRED12,
     &            NHCINT)
      IF (.NOT. GEN34) THEN
C     ... 980826-hjaaj: check if no contraction and no omissions
         NOCNT = .TRUE.
         DO IJ = 1,NORB34
            IF (NUCS34(IJ) .NE. 1) THEN
               NOCNT = .FALSE.
               GO TO 200
            END IF
         END DO
  200    IF (IPRINT .GT. 5) WRITE (LUPRI,*) 'NOCNT for C2DRIV:', NOCNT
      ELSE
         NOCNT = .FALSE.
      END IF
      IF (V12INT .OR. (R12EIN .AND. .NOT. U12INT))
     &CALL C2DRIV(AOINT(1,IADV12),
     &            WORK(KHCINT),COEF34,INDHER,CONT3,CONT4,
     &            WORK(KLAST),LWRK,NPCO3,NPCO4,NUCS34,IPRINT,LMNV12,
     &            LMNV34,WORK(KODD12),WORK(KODD34),NPNT34,NRED34)
      IF (R12INT)
     &CALL C2DRIV(AOINT(1,IADR12),
     &            WORK(KHRINT),COEF34,INDHER,CONT3,CONT4,
     &            WORK(KLAST),LWRK,NPCO3,NPCO4,NUCS34,IPRINT,LMNV12,
     &            LMNV34,WORK(KODD12),WORK(KODD34),NPNT34,NRED34)
      IF (U12INT)
     &CALL U2DRIV(AOINT(1,IADUUU),
     &            WORK(KHCINT),WORK(KHRINT),COEF34,INDHER,CONT3,CONT4,
     &            WORK(KLAST),LWRK,NPCO3,NPCO4,NUCS34,IPRINT,LMNV12,
     &            LMNV34,WORK(KODD12),WORK(KODD34),NPNT34,NRED34,R12EIN)
      LWTOT = LWTOT - KLAST
      RETURN
      END
C  /* Deck ndistr */
      FUNCTION NDISTR(INDXAB,MAXDIS,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
      LOGICAL DCMPAB
      DIMENSION INDXAB(MAXDIS)
#include "twocom.h"
#include "symmet.h"
C
#include "ibtfun.h"
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine NDISTR',-1)
C
      IJ = 0
C
C     Loop over components
C
      DO 100 NA = 1,KHKTA
         NSTRNA = NSTRA + NA
         ITYNA  = ISYMAO(NHKTA,NA)
         KHKTBB = KHKTB
         IF (DIAGAB) KHKTBB = NA
         DO 200 NB = 1,KHKTBB
            NSTRNB = NSTRB + NB
            ITYNB  = ISYMAO(NHKTB,NB)
            DCMPAB = SHAEQB .AND. NA .EQ. NB
C
C           Loop over symmetries
C
            DO 300 IREPA = 0, MAXREP
            IF (IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
               DO 400 IREPB = 0, MAXREP
               IF (IBTAND(MULB,IBTXOR(IREPB,ITYNB)) .EQ. 0) THEN
C
C                 Loop over contracted orbitals
C
                  DO 500 IA = 1,NORBA
                     INDA = IPTSYM(NSTRNA + KHKTA*(IA-1),IREPA)
                     NORBBB = NORBB
                     IF (DCMPAB) NORBBB = IA
                     IF (DCMPAB .AND. IREPA.LT.IREPB) NORBBB = IA - 1
                     DO 600 IB = 1,NORBBB
                        IJ = IJ + 1
                        INDB = IPTSYM(NSTRNB+KHKTB*(IB-1),IREPB)
                        MAXAB = MAX(INDA,INDB)
                        INDXAB(IJ) = MAXAB*(MAXAB-1)/2 + MIN(INDA,INDB)
 600                 CONTINUE
 500              CONTINUE
C
               END IF
 400           CONTINUE
            END IF
 300        CONTINUE
C
 200     CONTINUE
 100  CONTINUE
      NDISTR = IJ
      IF (NDISTR .GT. MAXDIS) THEN
         WRITE (LUPRI,'(A)')    ' Too many distributions in NDISTR.'
         WRITE (LUPRI,'(A,I6)') ' Number of distributions:',NDISTR
         WRITE (LUPRI,'(A,I6)') ' Maximum allowed:        ',MAXDIS
         CALL QUIT('Too many distributions required in NDISTR')
      END IF
      RETURN
      END
C  /* Deck paoset */
      SUBROUTINE PAOSET(JSTRSH,NPRIMS,NCONTS,KORBSH)
#include "implicit.h"
#include "maxorb.h"
#include "aovec.h"
C
C     DSM = parameter to determine zero coefficients for segmented contractions
      PARAMETER (DSM = 1.0D-15)
C
      DIMENSION KORBSH(MXSHEL,MXAOVC), IPRVPR(MXSHEL)
      DIMENSION JSTRSH(MXSHEL,MXAOVC,2), NPRIMS(MXSHEL,MXAOVC,2),
     &          NCONTS(MXSHEL,MXAOVC,2)
C
#include "shells.h"
#include "blocks.h"
#include "primit.h"
#include <priunit.h>
C
C     Parameters for sets of primitives
C     =================================
C
C     Determine NPRIMS(MAXSHL,NSETSH,2), NCONTS(MAXSHL,NSETSH,2),
C           and JSTRSH(MAXSHL,NSETSH,2)
C
      INUC = 0
      IPRVPR(1) = 0
      DO 100 I = 1, KMAX - 1
         IF (NUMCF(I + 1) .EQ. 1) INUC = INUC + NUCO(I)
         IPRVPR(I + 1) = INUC
  100 CONTINUE
      DO 200 I = 1, MAXSHL
#include "scalar.h"
         DO 300 J = 1, NSETSH(I,1)
            NPRIMS(I,J,1) = NUCO  (KORBSH(I,J))
            NCONTS(I,J,1) = NRCO  (KORBSH(I,J))
            JSTRSH(I,J,1) = IPRVPR(KORBSH(I,J))
  300    CONTINUE
C
C        Segmented transformation matrix: Each contracted is one set
C
         IF (SEGMSH(I)) THEN
            ISET = 0
            DO 600 J = 1, NSETSH(I,1)
               JSTR = IPRVPR(KORBSH(I,J))
               NPRS = NUCO(KORBSH(I,J))
               NCRS = NRCO(KORBSH(I,J))
               IOFF = JSTR
               DO 700 K = 1, NCRS
                  ISET = ISET + 1
                  NPRI = NDXGTA(NPRS,DSM,PRICCF(JSTR+1,K),1)
                  NPRIMS(I,ISET,2) = NPRI
                  NCONTS(I,ISET,2) = -K
                  JSTRSH(I,ISET,2) = IOFF
                  IOFF = IOFF + NPRI
  700          CONTINUE
  600       CONTINUE
            NSETSH(I,2) = ISET
         END IF
  200 CONTINUE
      RETURN
      END
C  /* Deck cprlop */
      SUBROUTINE CPRLOP(IPNTAO,IPNTOP,IPNTNO,IPNTRP,IPNTLG,SQ12EL,
     &                  IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
      LOGICAL SQ12EL, DCMPAB, DCMPCD, DCABAB, IPNTLG(3,NINTMX,*)
      DIMENSION IPNTAO(NINTMX,*), IPNTOP(3,NINTMX,*),
     &          IPNTNO(4,NINTMX,*), IPNTRP(3,NINTMX,*)
#include "twocom.h"
#include "symmet.h"
#include "ibtfun.h"
C
      IF (IPRINT .GT. 5) THEN
        CALL TITLER('Output from CPRLOP','*',103)
        IF (IPRINT .GT. 10) THEN
          WRITE (LUPRI,'(2X,A,4I5)')'NHKT?     ',NHKTA,NHKTB,NHKTC,NHKTD
          WRITE (LUPRI,'(2X,A,4I5)')'KHKT?     ',KHKTA,KHKTB,KHKTC,KHKTD
          WRITE (LUPRI,'(2X,A,4I5)')'MUL?      ',MULA,MULB,MULC,MULD
          WRITE (LUPRI,'(2X,A,2L5)')'DIAGAB/CD ',DIAGAB,DIAGCD
          WRITE (LUPRI,'(2X,A, L5)')'SQ12EL    ',SQ12EL
          WRITE (LUPRI,'(2X,A, I5)')'NOABCD    ',NOABCD
          WRITE (LUPRI,'(2X,A, I5)')'NOPREP    ',NOPREP
        END IF
      END IF
      DO 100 IREPX = 1, NOPREP
      IREPE = IPTREP(IREPX,1)
      IF (.NOT.(SHAEQB .OR. SHCEQD .OR. SHABAB)) THEN
         IAOFF = 1
         ISOFF = 1
         DO 200 ICOMPA = 1,KHKTA
           ITYNA  = ISYMAO(NHKTA,ICOMPA)
           NSTRNA = NSTRA + ICOMPA
         DO 200 ICOMPB = 1,KHKTB
           ITYNB  = ISYMAO(NHKTB,ICOMPB)
           NSTRNB = NSTRB + ICOMPB
         DO 200 ICOMPC = 1,KHKTC
           ITYNC  = ISYMAO(NHKTC,ICOMPC)
           NSTRNC = NSTRC + ICOMPC
         DO 200 ICOMPD = 1,KHKTD
           ITYND  = ISYMAO(NHKTD,ICOMPD)
           NSTRND = NSTRD + ICOMPD
           DO 300 IREPA = 0, MAXREP
           IF (IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
             IRPAE = IBTXOR(IREPA,IREPE)
             DO 310 IREPB = 0, MAXREP
             IF (IBTAND(MULB,IBTXOR(IREPB,ITYNB)) .EQ. 0) THEN
               IRPABE = IBTXOR(IREPB,IRPAE)
               IRPTYB = IBTXOR(IREPB,ITYNB)
               DO 320 IREPC = 0, MAXREP
               IF (IBTAND(MULC,IBTXOR(IREPC,ITYNC)) .EQ. 0) THEN
                 IREPD = IBTXOR(IREPC,IRPABE)
                 IF (IBTAND(MULD,IBTXOR(IREPD,ITYND)) .EQ. 0) THEN
                   IPNTAO(ISOFF,IREPX)   = IAOFF
                   IPNTOP(1,ISOFF,IREPX) = IRPTYB
                   IPNTOP(2,ISOFF,IREPX) = IBTXOR(IREPC,ITYNC)
                   IPNTOP(3,ISOFF,IREPX) = IBTXOR(IREPD,ITYND)
                   IPNTNO(1,ISOFF,IREPX) = NSTRNA
                   IPNTNO(2,ISOFF,IREPX) = NSTRNB
                   IPNTNO(3,ISOFF,IREPX) = NSTRNC
                   IPNTNO(4,ISOFF,IREPX) = NSTRND
                   IPNTRP(1,ISOFF,IREPX) = IREPA
                   IPNTRP(2,ISOFF,IREPX) = IREPB
                   IPNTRP(3,ISOFF,IREPX) = IREPC
                   ISOFF = ISOFF + 1
                 END IF
               END IF
  320          CONTINUE
             END IF
  310        CONTINUE
           END IF
  300      CONTINUE
           IAOFF = IAOFF + 1
  200    CONTINUE
      ELSE
         IAOFF = 1
         ISOFF = 1
         DO 400 ICOMPA = 1,KHKTA
           ITYNA  = ISYMAO(NHKTA,ICOMPA)
           NSTRNA = NSTRA + ICOMPA
           KHKTBB = KHKTB
           IF (DIAGAB) KHKTBB = ICOMPA
         DO 400 ICOMPB = 1,KHKTBB
           ITYNB  = ISYMAO(NHKTB,ICOMPB)
           NSTRNB = NSTRB + ICOMPB
         DO 400 ICOMPC = 1,KHKTC
           ITYNC  = ISYMAO(NHKTC,ICOMPC)
           NSTRNC = NSTRC + ICOMPC
           KHKTDD = KHKTD
           IF (DIAGCD) KHKTDD = ICOMPC
         DO 400 ICOMPD = 1,KHKTDD
           IF (SQ12EL .OR. .NOT.(SHABAB .AND. (ICOMPA .LT. ICOMPC .OR.
     &         (ICOMPA .EQ. ICOMPC .AND. ICOMPB .LT. ICOMPD)))) THEN
             ITYND  = ISYMAO(NHKTD,ICOMPD)
             NSTRND = NSTRD + ICOMPD
             DCMPAB = SHAEQB.AND.ICOMPA.EQ.ICOMPB
             DCMPCD = SHCEQD.AND.ICOMPC.EQ.ICOMPD
             DCABAB = .NOT.SQ12EL .AND. SHABAB .AND.
     &                ICOMPA.EQ.ICOMPC .AND. ICOMPB.EQ.ICOMPD
             DO 500 IREPA = 0, MAXREP
             IF (IBTAND(MULA,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
               IRPAE = IBTXOR(IREPA,IREPE)
               DO 510 IREPB = 0, MAXREP
               IF (IBTAND(MULB,IBTXOR(IREPB,ITYNB)) .EQ. 0) THEN
                 IRPABE = IBTXOR(IREPB,IRPAE)
                 IRPTYB = IBTXOR(IREPB,ITYNB)
                 DO 520 IREPC = 0, MAXREP
                 IF (IBTAND(MULC,IBTXOR(IREPC,ITYNC)) .EQ. 0) THEN
                   IREPD = IBTXOR(IREPC,IRPABE)
                   IF (IBTAND(MULD,IBTXOR(IREPD,ITYND)) .EQ. 0) THEN
                   IF (.NOT.(DCABAB .AND. (IREPA .LT. IREPC .OR.
     &               (IREPA.EQ.IREPC .AND. IREPB.LT.IREPD)))) THEN
                     IPNTAO(ISOFF,IREPX)   = IAOFF
                     IPNTOP(1,ISOFF,IREPX) = IRPTYB
                     IPNTOP(2,ISOFF,IREPX) = IBTXOR(IREPC,ITYNC)
                     IPNTOP(3,ISOFF,IREPX) = IBTXOR(IREPD,ITYND)
                     IPNTNO(1,ISOFF,IREPX) = NSTRNA
                     IPNTNO(2,ISOFF,IREPX) = NSTRNB
                     IPNTNO(3,ISOFF,IREPX) = NSTRNC
                     IPNTNO(4,ISOFF,IREPX) = NSTRND
                     IPNTRP(1,ISOFF,IREPX) = IREPA
                     IPNTRP(2,ISOFF,IREPX) = IREPB
                     IPNTRP(3,ISOFF,IREPX) = IREPC
                     IPNTLG(1,ISOFF,IREPX) = DCMPAB
                     IPNTLG(2,ISOFF,IREPX) = DCMPCD
                     IPNTLG(3,ISOFF,IREPX) = DCABAB
                     ISOFF = ISOFF + 1
                   END IF
                   END IF
                 END IF
  520            CONTINUE
               END IF
  510          CONTINUE
             END IF
  500        CONTINUE
           END IF
           IAOFF = IAOFF + 1
  400    CONTINUE
      END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck twoper */
      SUBROUTINE TWOPER
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
      PARAMETER (DMEGA = 1.0D6)
      INTEGER*8 N2CART, N2PRIM
#include "twosta.h"
#include "nuclei.h"
C
      CALL HEADER('Performance of TWOINT',1)
      N2CART = NBASIS*(NBASIS+1)/2
      N2CART = N2CART*(N2CART+1)/2
      N2PRIM = NPBAS *(NPBAS +1)/2
      N2PRIM = N2PRIM*(N2PRIM+1)/2
      WRITE (LUPRI,'(2(2X,A,I5,/),/,3(2X,A,I12,/))')
     &   'Number of contracted functions:',NBASIS,
     &   'Number of primitive  functions:',NPBAS,
     &   'Number of integrals written                    :',N2WRIT,
     &   'Number of contracted integrals (counting zeros):',N2CART,
     &   'Number of primitive integrals (counting zeros) :',N2PRIM
      IF (N2WRIT .GT. 0) THEN
         PRFNON = N2WRIT
         PRFNON = DMEGA*T2INT/ PRFNON
         PRFCAR = N2CART
         PRFCAR = DMEGA*T2INT/ PRFCAR
         PRFPRM = N2PRIM
         PRFPRM = DMEGA*T2INT/ PRFPRM
         WRITE (LUPRI,'(3(/,2X,A,F7.1,A))')
     &      'Average time for non-zero (written) integrals:',
     &       PRFNON,' microseconds',
     &      'Average time for zero/non-zero integrals     :',
     &       PRFCAR,' microseconds',
     &      'Average time for primitive integrals         :',
     &       PRFPRM,' microseconds'
      ELSE
         PRFCAR = N2CART
         PRFCAR = DMEGA*T2INT/ PRFCAR
         PRFPRM = N2PRIM
         PRFPRM = DMEGA*T2INT/ PRFPRM
         WRITE (LUPRI,'(2(/,2X,A,F7.1,A))')
     &      'Average time for zero/non-zero integrals     :',
     &       PRFCAR,' microseconds',
     &      'Average time for primitive integrals         :',
     &       PRFPRM,' microseconds'
      END IF
      RETURN
      END
C  /* Deck delsta */
      SUBROUTINE DELSTA(ITYPE,ICOUNT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.D0, D100 = 100.0D0, DMEGA = 1.D6)
      CHARACTER*1 SPDCAR
      DIMENSION TMEINT(MXQN,MXQN,MXQN,MXQN), INTWRT(MXQN,MXQN,MXQN,MXQN)
#include "twocom.h"
#include "twosta.h"
      SAVE TMEINT, INTWRT, TIMOLD
C
      IF (.NOT.TKTIME) RETURN
      IF (ITYPE .EQ. -1) THEN
         CALL GETTIM(TIMOLD,DUMTIM)
         CALL DZERO(TMEINT,MXQN**4)
         CALL IZERO(INTWRT,MXQN**4)
      ELSE IF (ITYPE .EQ. 0) THEN
         CALL GETTIM(TIMNOW,DUMTIM)
         TMEINT(NHKTA,NHKTB,NHKTC,NHKTD)
     & = TMEINT(NHKTA,NHKTB,NHKTC,NHKTD) + TIMNOW - TIMOLD
         INTWRT(NHKTA,NHKTB,NHKTC,NHKTD)
     & = INTWRT(NHKTA,NHKTB,NHKTC,NHKTD) + ICOUNT
         TIMOLD = TIMNOW
      ELSE
         TOTTIM = DSUM(MXQN**4,TMEINT,1)
         NINTOT = ISUM(MXQN**4,INTWRT,1)
         CALL HEADER('               # int'/
     &        /'               time               time/int',3)
         DO 100 I = 0, 4*(MXQN - 1)
            NINTI = 0
            TIMEI = D0
            DO 200 IA = 0, MXQN - 1
            DO 200 IB = 0, MIN(I - IA,MXQN - 1)
            DO 200 IC = 0, MIN(I - IA - IB,MXQN - 1)
               ID = I - IA - IB - IC
               IF (ID .LT. MXQN) THEN
                  NINTI = NINTI + INTWRT(IA+1,IB+1,IC+1,ID+1)
                  TIMEI = TIMEI + TMEINT(IA+1,IB+1,IC+1,ID+1)
               END IF
  200       CONTINUE
            IF (NINTI .GT. 0) THEN
               PERCNT = D100*NINTI
               PERCNT = PERCNT / NINTOT
               PERINT = (DMEGA*TIMEI) / NINTI
               WRITE (LUPRI,2000) ' j =',I,
     &             NINTI, ' (',PERCNT,'%)',
     &             TIMEI, ' (',D100*TIMEI/TOTTIM,'%)',PERINT
               DO 300 IA = 0, MXQN - 1
               DO 300 IB = 0, MIN(I - IA,MXQN - 1)
               DO 300 IC = 0, MIN(I - IA - IB,MXQN - 1)
                  ID = I - IA - IB - IC
                  IF (ID .LT. MXQN) THEN
                     NINTX = INTWRT(IA+1,IB+1,IC+1,ID+1)
                     IF (NINTX .GT. 0) THEN
                        TIME  = TMEINT(IA+1,IB+1,IC+1,ID+1)
                        PERINT = DMEGA*TIME/NINTX
                        WRITE (LUPRI,1000)'(',SPDCAR(IA),SPDCAR(IB),
     &                                    '|',SPDCAR(IC),SPDCAR(ID),')',
     &                                        NINTX, TIME, PERINT
                     END IF
                  END IF
  300          CONTINUE
               WRITE (LUPRI,'()')
            END IF
  100    CONTINUE
      END IF
      RETURN
 1000 FORMAT(5X,7A1,2X,I10,11X,F8.2,14X,F8.2)
 2000 FORMAT(5X,A,I2,3X,I10,A,F5.2,A,2X,F8.2,A,F5.2,A,5X,F8.2)
      END
C  /* Deck redsta */
      SUBROUTINE REDSTA(ITYPE,NREDPR,NTOTPR,NREDCN,NTOTCN)
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.D0, D1 = 1.0D0, D100 = 100.0D0)
      SAVE NODC, REDPR, REDCN
C
      IF (ITYPE .EQ. -1) THEN
         NODC  = 0
         REDPR = D0
         REDCN = D0
      ELSE IF (ITYPE .EQ. 0) THEN
         NODC  = NODC + 1
         XNREDPR = NREDPR
         XNREDCN = NREDCN
         REDPR = REDPR + XNREDPR / NTOTPR
         REDCN = REDCN + XNREDCN / NTOTCN
      ELSE
         WRITE (LUPRI,'(2(/1X,A,F5.2,A))')
     &      ' Average reduction of primitive ODC vectors:   ',
     &      D100*(D1 - REDPR/NODC),' %',
     &      ' Average reduction of contracted ODC vectors:  ',
     &      D100*(D1 - REDCN/NODC),' %'
      END IF
      RETURN
      END
C  /* Deck nodsym */
      FUNCTION NODSYM(MAXOPR,MULA,MULB)
#include "implicit.h"
#include "ibtfun.h"
      NODSYM = 0
      DO 100 I = 0, MAXOPR
         IF (IBTAND(I,IBTOR(MULA,MULB)).EQ.0) NODSYM = NODSYM + 1
  100 CONTINUE
      RETURN
      END
C  /* Deck herprp */
      SUBROUTINE HERPRP(INDHER,INDHSQ,IODDHR)
#include "implicit.h"
#include "priunit.h"
      INTEGER T, U, V, TUV
      DIMENSION INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7)
#include "hertop.h"
#include "ibtfun.h"
C
      TUV = 0
      DO 100 J = 0, JTOP
         DO 110 T = J, 0, -1
            DO 120 U = J - T, 0, -1
               V = J - T - U
               TUV = TUV + 1
               INDHER(T,U,V) = TUV
               INDHSQ(TUV)   = T + U*(JTOP + 1) + V*(JTOP + 1)**2
               IBTUV = MOD(T,2) + 2*MOD(U,2) + 4*MOD(V,2)
               DO 200 ISYM = 0, 7
                  IODDHR(TUV,ISYM) = IBTAND(IBTUV,ISYM)
  200          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck getlmn */
      SUBROUTINE GETLMN(LMNVLS,IPRINT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION LMNVLS(KCKMAX,5,4)
#include "twocom.h"
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from GETLMN','*',103)
C
      CALL GETLM1(NHKTA,KCKTA,KHKTA,LMNVLS(1,1,1),'A',IPRINT)
      CALL GETLM1(NHKTB,KCKTB,KHKTB,LMNVLS(1,1,2),'B',IPRINT)
      CALL GETLM1(NHKTC,KCKTC,KHKTC,LMNVLS(1,1,3),'C',IPRINT)
      CALL GETLM1(NHKTD,KCKTD,KHKTD,LMNVLS(1,1,4),'D',IPRINT)
      RETURN
      END
C  /* Deck getlm1 */
      SUBROUTINE GETLM1(NHKTX,KCKTX,KHKTX,LMN,TYPE,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      CHARACTER*1 TYPE
      DIMENSION LMN(KCKMAX,5)
#include "ccom.h"
#include "sphtrm.h"
#include "twocom.h"
      CALL LMNVAL(NHKTX,KCKTX,LMN(1,1),LMN(1,2),LMN(1,3))
      DO 100 I = 1, KCKTX
         LMN(I,4) = MOD(LMN(I,1),2)+2*MOD(LMN(I,2),2)+4*MOD(LMN(I,3),2)
  100 CONTINUE
      DO 200 I = 1, KHKTX
         IF (DOCART) THEN
            LMN(I,5) = LMN(I,4)
         ELSE
            INDMAX = IDMAX(KCKTX,CSP(ISPADR(NHKTX)+I-1),KHKTX)
            LMN(I,5) = LMN(INDMAX,4)
         END IF
  200 CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Cartesian components for orbital '//TYPE,1)
         WRITE (LUPRI,'(5X,A,/)')'   cmp        x    y    z        odd '
         DO 300 I = 1, KCKTX
            WRITE (LUPRI,'(5X,I5,5X,3I5,5X,I5)') I, (LMN(I,J),J=1,4)
  300    CONTINUE
         WRITE (LUPRI,'(/,A,/)') ' Odd for spherical integrals:'
         WRITE (LUPRI,'(15I5)') (LMN(I,5),I=1,KHKTX)
      END IF
      RETURN
      END
C  /* Deck getoda */
      SUBROUTINE GETODA(IODD12,IODD34,LMNV12,LMNV34,IPRINT)
C
C     TUH
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION LMNV12(KCKMX,5,2), LMNV34(KCKMX,5,2),
     &          IODD12(KCKT12,2),  IODD34(KCKT34,2)
#include "twoao.h"
#include "ibtfun.h"
C
      IF (IHRSYM .GT. 0) THEN
         IJ = 0
         DO 100 I = 1, KCKT1
            IODD = LMNV12(I,4,1)
            MAXJ = KCKT2
            IF (DIAC12) MAXJ = I
            DO 110 J = 1, MAXJ
               IJ = IJ + 1
               IODD12(IJ,1) = IBTAND(IHRSYM,IBTXOR(IODD,LMNV12(J,4,2)))
  110       CONTINUE
  100    CONTINUE
C
         IJ = 0
         DO 200 I = 1, KHKT1
            IODD = LMNV12(I,5,1)
            MAXJ = KHKT2
            IF (DIAG12) MAXJ = I
            DO 210 J = 1, MAXJ
               IJ = IJ + 1
               IODD12(IJ,2) = IBTAND(IHRSYM,IBTXOR(IODD,LMNV12(J,5,2)))
  210       CONTINUE
  200    CONTINUE
C
         IJ = 0
         DO 300 I = 1, KCKT3
            IODD = LMNV34(I,4,1)
            MAXJ = KCKT4
            IF (DIAC34) MAXJ = I
            DO 310 J = 1, MAXJ
               IJ = IJ + 1
               IODD34(IJ,1) = IBTAND(IHRSYM,IBTXOR(IODD,LMNV34(J,4,2)))
  310       CONTINUE
  300    CONTINUE
C
         IJ = 0
         DO 400 I = 1, KHKT3
            IODD = LMNV34(I,5,1)
            MAXJ = KHKT4
            IF (DIAG34) MAXJ = I
            DO 410 J = 1, MAXJ
               IJ = IJ + 1
               IODD34(IJ,2) = IBTAND(IHRSYM,IBTXOR(IODD,LMNV34(J,5,2)))
  410       CONTINUE
  400    CONTINUE
      ELSE
         CALL IZERO(IODD12,2*KCKT12)
         CALL IZERO(IODD34,2*KCKT34)
      END IF
C
      IF (IPRINT .GE. 10) THEN
         CALL TITLER('Output from GETODA','*',103)
         CALL HEADER('IODD12 (Cartesian)',1)
         WRITE (LUPRI,'(1X,(10I7))') (IODD12(I,1),I=1,KCKT12)
         CALL HEADER('IODD12 (spherical)',1)
         WRITE (LUPRI,'(1X,(10I7))') (IODD12(I,2),I=1,KHKT12)
         CALL HEADER('IODD34 (Cartesian)',1)
         WRITE (LUPRI,'(1X,(10I7))') (IODD34(I,1),I=1,KCKT34)
         CALL HEADER('IODD34 (spherical)',1)
         WRITE (LUPRI,'(1X,(10I7))') (IODD34(I,2),I=1,KHKT34)
      END IF
      RETURN
      END
C  /* Deck whtrep */
      SUBROUTINE WHTREP(ITYPE,MULE,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      LOGICAL DOPREP(0:7)
#include "dorps.h"
#include "symmet.h"
#include "twocom.h"
#include "expopt.h"
#include "ibtfun.h"
C
      DO 100 I = 0, MAXREP
         DOPREP(I) = .FALSE.
  100 CONTINUE
C
C     Undifferentiated integrals written on file
C
      IF (ITYPE .EQ. 0 .OR. EXPGRA) THEN
         DOPREP(0) = .TRUE.
C
C     Differentiated integrals written on file (atom JATOM)
C
      ELSE IF (ITYPE .EQ. 1 .OR. ITYPE .EQ. 8) THEN
         DO 200 IREPE = 0, MAXREP
         IF (DOREPS(IREPE)) THEN
            DO 210 ICOOR = 1, 3
            IF (IBTAND(MULE,IBTXOR(IREPE,ISYMAX(ICOOR,1))).EQ.0) THEN
               DOPREP(IREPE) = .TRUE.
            END IF
  210       CONTINUE
         END IF
  200    CONTINUE
C
C     Differentiated integrals for expectation values
C
      ELSE IF (ITYPE .EQ. 2 .OR. ITYPE .EQ. 6) THEN
         DO 220 IREPE = 0, MAXREP
            DOPREP(IREPE) = DOREPS(IREPE)
  220    CONTINUE
C
C     Spin-orbit integrals
C
      ELSE IF (ITYPE .EQ. -2 .OR. ITYPE .EQ. -20) THEN
         DOPREP(ISYMAX(1,2)) = .TRUE.
         DOPREP(ISYMAX(2,2)) = .TRUE.
         DOPREP(ISYMAX(3,2)) = .TRUE.
C
C     Direct calculation of Fock matrix in AO basis
C
      ELSE IF (ITYPE .EQ. 3) THEN
         DOPREP(0) = .TRUE.
C
C     Direct calculations of Fock matrices in SO basis
C
      ELSE IF (ITYPE .EQ. 9) THEN
         DOPREP(0) = .TRUE.
C
C     Calculation of distributions
C
      ELSE IF (ITYPE .EQ. 4) THEN
         DOPREP(0) = .TRUE.
C
C     Magnetic field derivatives
C
      ELSE IF (ABS(ITYPE) .EQ. 5) THEN
         DOPREP(ISYMAX(1,2)) = .TRUE.
         DOPREP(ISYMAX(2,2)) = .TRUE.
         DOPREP(ISYMAX(3,2)) = .TRUE.
      ELSE IF (ITYPE .EQ. 7) THEN
         DOPREP(0) = .TRUE.
      ELSE IF (ITYPE.EQ.10) THEN
         DOPREP(0) = .TRUE.
         DOPREP(ISYMAX(1,2)) = .TRUE.
         DOPREP(ISYMAX(2,2)) = .TRUE.
         DOPREP(ISYMAX(3,2)) = .TRUE.
      ELSE
      END IF
C
      NOPREP = 0
      DO 400 IREP = 0, MAXREP
         IF (DOPREP(IREP)) THEN
            NOPREP = NOPREP + 1
            IPTREP(NOPREP,1) = IREP
            IPTREP(IREP  ,2) = NOPREP
         END IF
  400 CONTINUE
C
      IF (IPRINT .GT. 5) THEN
         CALL TITLER('Output from WHTREP','*',103)
         WRITE (LUPRI,'(/,1X,A,I2,/,1X,A,8I2)')
     &         ' Number of perturbation representations:', NOPREP,
     &         ' Perturbation representations:          ',
     &           (IPTREP(I,1),I=1,NOPREP)
      END IF
      RETURN
      END
C  /* Deck seteff */
      SUBROUTINE SETEFF(IEFFB,IEFFC,IEFFD)
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      DIMENSION IEFFB(0:7), IEFFC(0:7), IEFFD(0:7)
#include "symmet.h"
#include "twocom.h"
#include "ibtfun.h"
C
      IOFF = 0
      DO 100 I = 0, MAXOPR
        IF (IBTAND(I,MULB).EQ.0) THEN
            IEFFB(I) = IOFF
            IOFF = IOFF + KHKTB
         ELSE
            IEFFB(I) = IEFFB(IBTXOR(I,IBTAND(I,MULB)))
         END IF
 100  CONTINUE
      IOFF = 0
      DO 200 I = 0, MAXOPR
         IF (IBTAND(I,MULC).EQ.0) THEN
            IEFFC(I) = IOFF
            IOFF = IOFF + KHKTC
         ELSE
            IEFFC(I) = IEFFC(IBTXOR(I,IBTAND(I,MULC)))
         END IF
 200  CONTINUE
      IOFF = 0
      DO 300 I = 0, MAXOPR
         IF (IBTAND(I,MULD).EQ.0) THEN
            IEFFD(I) = IOFF
            IOFF = IOFF + KHKTD
         ELSE
            IEFFD(I) = IEFFD(IBTXOR(I,IBTAND(I,MULD)))
         END IF
 300  CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Tworun */
      SUBROUTINE TWORUN(ITYPE,MAXDIF,JATOM,PERTUR,EXPECT,UNDIFF,
     &           DIRFCK,SOFOCK,SPNORB,DIA2SO,ZFS2EL,DISTRI,SQ12EL,
     &           LONDON,SUSCEP,DDFOCK,ADISTR,MAXDER,IATOM,MULE,MULTE,
     &           IPRINT)
C*****************************************************************************
C
C     Determine run type for two-electron integrals
C
C     Compiled by T.Saue Oct 7 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      LOGICAL PERTUR, EXPECT, UNDIFF, DIRFCK, SPNORB, DIA2SO, DISTRI,
     &        SQ12EL, LONDON, SUSCEP, DDFOCK, ADISTR, SOFOCK, ZFS2EL
#include "nuclei.h"
#include "symmet.h"
C
C
C     *****************************
C     **** Determine run type *****
C     *****************************
C
      PERTUR = .FALSE.
      EXPECT = .FALSE.
      UNDIFF = .FALSE.
      DIRFCK = .FALSE.
      SOFOCK = .FALSE.
      SPNORB = .FALSE.
      DIA2SO = .FALSE.
      ZFS2EL = .FALSE.
      DISTRI = .FALSE.
      SQ12EL = .FALSE.
      LONDON = .FALSE.
      SUSCEP = .FALSE.
      DDFOCK = .FALSE.
      ADISTR = .FALSE.
      IATOM  = 0
      MULTE  = 1
C
C     a) Undifferentiated integrals (or 2-el.Darwin) written on file
C     ==============================================================
C
      IF (ITYPE .EQ. 0) THEN
         UNDIFF = .TRUE.
         MAXDER = 0
C
C     b) Differentiated integrals written on file (atom JATOM)
C     ========================================================
C
      ELSE IF (ITYPE .EQ. 1) THEN
         PERTUR = .TRUE.
         MAXDER = MAXDIF
         IATOM  = JATOM
         MULE   = ISTBNU(IATOM)
         MULTE  = MULT(MULE)
C
C     c) Expectation values of differentiated integrals (all atoms)
C     =============================================================
C
      ELSE IF (ITYPE .EQ. 2) THEN
         EXPECT = .TRUE.
         MAXDER = MAXDIF
C
C     d) Direct calculation of Fock matrices
C        using skeleton matrix approach (AO basis)
C     ============================================
C
      ELSE IF (ITYPE .EQ. 3) THEN
         DIRFCK = .TRUE.
         MAXDER = MAXDIF
C
C     e) Direct calculation of Fock matrices in SO basis
C     ==================================================
C
      ELSE IF (ITYPE .EQ. 9) THEN
         SOFOCK = .TRUE.
         MAXDER = MAXDIF
C
C     f) Spin-orbit integrals
C     =======================
C
      ELSE IF (ITYPE .EQ. -2 .OR. ITYPE .EQ. -20) THEN
         SPNORB = .TRUE.
         SQ12EL = .TRUE.
         MAXDER = 1
         IF (ITYPE.EQ.-20) DDFOCK = .TRUE.
C
C     g) Integral distributions
C     =========================
C
      ELSE IF (ITYPE .EQ. 4) THEN
         DISTRI = .TRUE.
         SQ12EL = .TRUE.
         MAXDER = 0
C
C     h) Derivatives with respect to magnetic field
C     =============================================
C
      ELSE IF (ITYPE .EQ. 5 .OR. ITYPE .EQ. - 5) THEN
         LONDON = .TRUE.
         MAXDER = MAXDIF
         IF (MAXDER .EQ. 2) SUSCEP = .TRUE.
         IF (ITYPE .EQ. -5) DDFOCK = .TRUE.
C
C     i) Expectation values and Fock matrics of differentiated integrals
C     ==================================================================
C
      ELSE IF (ITYPE .EQ. 6 .OR. ITYPE .EQ. -6) THEN
         EXPECT = .TRUE.
         MAXDER = MAXDIF
         DDFOCK = MAXDIF .GT. 1
         IF (ITYPE .EQ. -6) DDFOCK = .TRUE.
C
C     j) Distributions (all gabcd for fixed a)
C     ======================================
C
      ELSE IF (ITYPE .EQ. 7) THEN
         UNDIFF = .TRUE.
         ADISTR = .TRUE.
         MAXDER = 0
C
C     k) Derivative Fock matrix for specified atom
C     ============================================
C
      ELSE IF (ITYPE .EQ. 8) THEN
         PERTUR = .TRUE.
         MAXDER = MAXDIF
         IATOM  = JATOM
         MULE   = ISTBNU(IATOM)
         MULTE  = MULT(MULE)
         DDFOCK = .TRUE.
C
C     l) Diamagnetic two-electron spin-orbit integrals
C     ================================================
C
      ELSE IF (ITYPE .EQ. 10) THEN
         MAXDER = 2
         DIA2SO = .TRUE.
         SQ12EL = .TRUE.
C
C     m) Zero-field electron spin-spin couplings
C     ==========================================
C
      ELSE IF (ITYPE .EQ. 11) THEN
         MAXDER = 2
         ZFS2EL = .TRUE.
         SQ12EL = .FALSE.
C
C     x) Error
C     ========
C
      ELSE
         WRITE (LUPRI,'(1X,A,I5)')
     &       ' Wrong run type in TWORUN, ITYPE =', ITYPE
         CALL QUIT('ERROR in TWORUN, unknown run type.')
      END IF
C
C     Set symmetry of operator
C
      CALL WHTREP(ITYPE,MULE,IPRINT)
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(A,I5)') ' ITYPE  ', ITYPE
         WRITE (LUPRI,'(A,I5)') ' MAXDER ', MAXDER
         WRITE (LUPRI,'(A,I5)') ' IATOM  ', IATOM
         WRITE (LUPRI,'(A,I5)') ' MULTE  ', MULTE
         WRITE (LUPRI,'(A,L5)') ' UNDIFF ', UNDIFF
         WRITE (LUPRI,'(A,L5)') ' PERTUR ', PERTUR
         WRITE (LUPRI,'(A,L5)') ' EXPECT ', EXPECT
         WRITE (LUPRI,'(A,L5)') ' DDFOCK ', DDFOCK
         WRITE (LUPRI,'(A,L5)') ' DIRFCK ', DIRFCK
         WRITE (LUPRI,'(A,L5)') ' SOFOCK ', SOFOCK
         WRITE (LUPRI,'(A,L5)') ' SPNORB ', SPNORB
         WRITE (LUPRI,'(A,L5)') ' DIA2SO ', DIA2SO
         WRITE (LUPRI,'(A,L5)') ' ZFS2EL ', ZFS2EL
         WRITE (LUPRI,'(A,L5)') ' LONDON ', LONDON
         WRITE (LUPRI,'(A,L5)') ' SUSCEP ', SUSCEP
         WRITE (LUPRI,'(A,L5)') ' DISTRI ', DISTRI
         WRITE (LUPRI,'(A,L5)') ' ADISTR ', ADISTR
         WRITE (LUPRI,'(A,L5)') ' SQ12EL ', SQ12EL
      END IF
C
C
      RETURN
      END
C
C  /* Deck CLEARMEM */
      SUBROUTINE CLEAR_INCOREMEM()
C
#include "maxorb.h"
#include "incore.h"
C
      LINTSV = .FALSE.
      LINTMP = .FALSE.
      INITX = .FALSE.
      MSAVE = .TRUE.
      MMCORE = MMWORK
      LMCORE = MMCORE
      ISCORE = 1
      JSCORE = ISCORE
      N_SHL = 1
      I_SHL = 1
      INDX_SHL1 = 0
      INDX_SHL2 = 0
      INDX_SHL3 = 0
      INDX_SHL4 = 0
C
      RETURN
      END
C -- end of her2drv.F --
