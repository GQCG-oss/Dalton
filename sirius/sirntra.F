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
C
#ifdef OLD_REV_LOG
!===========================================================================
!MAERKE 940727 ITRSIR=10 og USEDRC (som i .HYPCAL) giver fuldt saet
! af Dirac integraler: dvs. alle integraler baade som Coulomb og Dirac !!!!
!MAERKE 940520 Brug DARECL til LUHLFx, hvor recl kendes
!MAERKE 931201 external option for KTRLVL ??
!MAERKE 900620 new DRCCTL which calls TRA2 to add Dirac integrals.
!MAERKE 900322 definer /TRHSO/ ILXYZ, LUAHSO, LUMHSO;
!              definer IDATA (i.e. LADX,NPQ) for LUAHSO
!              anvend IADOUT(1:3,*) for lx, ly, og lz !!
!/* Comdeck notes */
!se efter MAERKE
!NOTE: 891127 LBINTM skal vaere max record paa LUMINT
!             (bruges i DRCCTL til memory estimate for nxth2m)
!NOTE: 900103 TRACHK routine which checks CMO, INTSYM, ITRLVL m.m.
!         MAERKE 900307: naeste linjer er fra motra2.l02.u
!         ?? anvend #include "tralistgs.h"
!NOTE: 900320 use NCOINT,NE1INT,NE2INT to count integrals written to LUMINT
!=======================
!/* Comdeck rev_log */
!--------> lav TRACHK rutine (se ovenover) 940520/hjaaj !!!!!!!!!!
!-- version n06 --
!950102-hjaaj
!TRA2E1: corrected OPHLF3 calculation
!-- version n05 --
!940727-hjaaj
!TRACTL: delete AOORDINT after transformation if ITRSIR .lt. 0
!TRALIM: comments on MISH, MASH reset for JTRLVL .gt. 2
!940520-hjaaj
!TRA2,TRA2E1: OPHLFx parameter, call DARMOV(LUHLFx) if true to delete
!   any large temporary files before exit from TRACTL.
!-- version n04 --
!931201-hjaaj:
!TRALVL: 'MOFCFV' is not defined any more.
!920529-hjaaj:
!TR2CTL: fixed error in LTUPQ calculation (s/NOCR/NOCQ/)
!        use NOCOUL in LTUPQ calculation
!-- version m02 --
!901016-hjaaj: TRA2: a) return immediately if DOTUPQ,DORUPQ,DOUSPQ all
!   false; b) changed DORUPQ and DOUSPQ definitions so now always false
!   if the corresponding integrals not needed for exchange.
!-- TRACTL: USEDRC reset to true if ITRLVL 1,2,3, or 4
!900924-hjaaj: define LUAINT as LUORD = LUDA = 14
!-- changed IBM names from ORDINT,MOLINT to AORDINT,MTWOINT
!900619-hjaaj: new routine TRALVL for translating from Sirius to
!              newtra transformation levels.
!900320-hjaaj: implemented ICDTRA in TRA2E1
!900316-hjaaj: added NXTH2D and NX2H2D routines
!900314-hjaaj: rewritten NX2H2M to outer loop over symmetries.
!===========================================================================
#endif
C /* Deck n_tractl */
      SUBROUTINE N_TRACTL(ITRSIR,CMO,WORK,KWORK)
C
C Mar 1990 Hans Joergen Aa. Jensen
C last revision 100915-hjaaj (no AOORDINT)
C
#include "implicit.h"
      DIMENSION CMO(*),WORK(KWORK)
C
#include "dummy.h"
C
C Used from common blocks:
C   INFTRA : IPRTRA,USEDRC
C CBGETDIS : IAD*
C   INFTAP : LUORDA
C
#include "priunit.h"
#include "inftra.h"
#include "cbgetdis.h"
#include "inftap.h"
C
      LOGICAL NOCOUL,NOEXCH,NOTUVX,MCOPT,DELORD
      LOGICAL AOEXIST,MOEXIST
C
      CALL QENTER('N_TRACTL')
C
C     (Re)calculate AO integrals, if necessary
C     /hjaaj 11-jan-2007
C
      IF (USE_INTSORT) THEN
         JRDAO = 0
         CALL GPINQ('AOORDINT','EXIST',AOEXIST)
      ELSE
         JRDAO = 3
         CALL GPINQ('AOTWOINT','EXIST',AOEXIST)
      END IF
      CALL GPINQ(FNINTM,'EXIST',MOEXIST)
C                default is FNINTM = 'MOTWOINT'
      IF (.NOT.(AOEXIST .AND. MOEXIST)) THEN
C        ... we must maybe (re)generate AOORDINT
C            but first we must test for AOTWOINT
         CALL MAKE_AOTWOINT(WORK,KWORK)
C
         IF (USE_INTSORT) THEN
            WRITE (LUPRI,'(/A)') '---> Generating AOORDINT'
            CALL TIMER('START ',TIMSTR,TIMEND)
            CALL SORTAO(WORK,KWORK)
            CALL TIMER('SORTAO',TIMSTR,TIMEND)
         END IF
         CALL FLSHFO(LUPRI)
C
      END IF
      CALL SIR_INTOPEN_NEWTRA
C
C     ***** Reset /CBGETD/ because any H2AC on disk will be
C           obsolete after transformation.
C
      IADINT = -1
      IADH2  = -1
      IADH2X = -1
      IADH2D = -1
C
C     Delete AOORDINT after transformation if ITRSIR .lt. 0
C
      DELORD = (ITRSIR .LT. 0)
      JTRSIR = ABS(ITRSIR)
      IF (JTRSIR .GE. 1 .AND. JTRSIR .LE. 4) USEDRC = .TRUE.
C     ... USEDRC must be true for second order transformation
C         can be true by input for higher order (e.g. for MP2 and RESPONS)
C
C     Translate from SIRIUS convention for transformation level
C     to the convention used in this program.
C     
      CALL TRALVL(JTRSIR,NTRLVL)
C
      NOTUVX = .TRUE.
      NOCOUL = .FALSE.
      NOEXCH = .NOT.USEDRC
      MCOPT  = .FALSE.
C

      IPRFIO = IPRTRA
C
      IF (IPRFIO .GE. 20) THEN
         CALL FASTIO('TRACE=ON')
C        CALL FASTIO('ERRFIX=FIRM')
C        no, ERRFIX=FIRM changes unit number to error code, if error.
      END IF
C
      IF (IPRFIO .GE. 3) CALL FASTIO('STATUS')
C
      IF (IPRTRA .GT. 0) DST = SECOND()
C
C
C     Set TRINP and TRUNIT
C
      CALL TRASET(JRDAO, KWORK)
C
      CALL TRACT2(NTRLVL,NOTUVX,NOCOUL,NOEXCH,MCOPT,CMO,DUMMY,
     &            WORK,KWORK)
C     CALL TRACT2(NTRLVL,NOTUVX,NOCOUL,NOEXCH,MCOPT,CMO,TUVX,
C    &            WRK,MEMX)
C
      IF (IPRTRA .GT. 0) THEN
         TRATIM = SECOND() - DST
         WRITE (LUPRI,'(/A,F10.2)')
     &      ' Total CPU time used in integral transformation:',TRATIM
      END IF
C
      IF (IPRFIO .GE. 3) CALL FASTIO('STATUS')
C
C     Reset FASTIO options to defaults, if changed:
C
      IF (IPRFIO .GE. 20) THEN
         CALL FASTIO('TRACE=OFF')
C        CALL FASTIO('ERRFIX=HARD')
      END IF
C
      IF (DELORD .AND. USE_INTSORT) THEN
         WRITE (LUPRI,'(/A)') ' AOORDINT will be deleted as requested.'
         CALL DARMOV(LUORDA)
      END IF
C
      IF (LUINTA .GT. 0) CALL GPCLOSE(LUINTA,'KEEP')
      CALL FLSHFO(LUPRI)
      CALL QEXIT('N_TRACTL')
      RETURN
      END
C  /* Deck tract2 */
      SUBROUTINE TRACT2(NTRLVL,NOTUVX,NOCOUL,NOEXCH,MCOPT,CMO,TUVX,
     &                  WRK,MEMX)
C
C     26-Aug-1989 Hans Joergen Aa. Jensen
C
C     Control of second order integral transformation, the work is
C     done in TR2CTL. Interface to Sirius.
C
C    Comments:
C       hj-890106: INTSYM = symmetry of integrals
C       In revised INTSORT KEEP(2) = INTSYM below
C
#include "implicit.h"
      DIMENSION CMO(*),TUVX(*),WRK(*)
      LOGICAL   NOTUVX, NOCOUL,NOEXCH, MCOPT
C
#include "maxorb.h"
C
C Used from common blocks:
C  TRINP  : MMORBX
C  INFTRA : IPRTRA
C
#include "trinp.h"
#include "inftra.h"
C
      CALL QENTER('TRACT2')
C
C     Allocate memory for ICDTRA(MMORBX)
C
      KFREE  = 1
      LFREE  = MEMX
      LCDTRA = MMORBX
      CALL MEMGET('INTE',KCDTRA,LCDTRA,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KTRTYP,MXCORB,WRK,KFREE,LFREE)
C
C     Set MO loop limits for the four integral indices
C
      CALL TRALIM(NTRLVL,WRK(KCDTRA),WRK(KTRTYP))
C
C     Initialize output integral file : LUMINT
C
      CALL TRAINI(NTRLVL,CMO,WRK(KCDTRA),WRK(KTRTYP))
C
C     Now we are ready for the transformation
C
      CALL TR2CTL(NOTUVX,NOCOUL,NOEXCH,MCOPT,CMO,TUVX,
     &            WRK(KCDTRA),WRK(KFREE),LFREE)
C
C
      CALL QEXIT ('TRACT2')
C
      RETURN
      END
C  /* Deck traadr */
      SUBROUTINE TRAADR(IADRM,IADFRE,IADCMO,IADCDT)
C
C  5-Oct-1990 Hans Joergen Aa. Jensen
C
C  Retrieve first free address, address of CMO and
C  address of ICDTRA from the LUMINT address array
C  (here IADRM).
C
#include "implicit.h"
      DIMENSION IADRM(3,1297)
      IADFRE = IADRM(1,1297)
      IADCMO = IADRM(2,1297)
      IADCDT = IADRM(3,1297)
      RETURN
      END
C  /* Deck sir_intopen_newtra */
      SUBROUTINE SIR_INTOPEN_NEWTRA
C
C 900309-hjaaj
C
#include "implicit.h"
C
C Used from common blocks:
C   priunit : LUPRI
C   INFORB  : NBAS()
C   INFTRA  : USE_INTSORT, NIBUF,NBITS,LBUF
C   INFTAP  : LUORDA,LUMINT,LUINTA
C
#include "priunit.h"
#include "inforb.h"
#include "inftra.h"
#include "inftap.h"
C
      LOGICAL FEXIST
      integer nbas_x(8)
C
C     Open old AOORDINT/AOTWOINT if it exists,
C     otherwise it will be attempted generated later.
C
      IF (USE_INTSORT) THEN
         CALL GPINQ('AOORDINT','EXIST',FEXIST)
         IF (FEXIST) THEN
           IF (LUORDA .LE. 0) CALL DAOPEN(LUORDA,'AOORDINT')
         END IF
      ELSE
         CALL GPINQ('AOTWOINT','EXIST',FEXIST)
         IF (FEXIST) THEN
            IF (LUINTA .LE. 0) CALL GPOPEN(LUINTA,AO2INTFILE_LABEL,
     &         'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
            LUORDA = LUINTA

            CALL REWSPL(LUINTA)
            CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
!           consistency checks before reading of AOTWOINT in RDAO2
            READ (LUINTA) NSYM_x, (NBAS_x(I),I=1,8), LBUF, NIBUF, NBITS
            ierr = 0
            IF (NSYM_x .ne. NSYM) THEN
               ierr = ierr + 1
               WRITE(LUPRI,*)
     &         'SIR_INTOPEN_NEWTRA: NSYM_AOTWOINT .ne. NSYM',NSYM_x,NSYM
            END IF
            do i = 1,8
               if (nbas(i) .ne. nbas_x(i)) ierr = ierr + 10
            end do
            IF (ierr.ge.10) THEN
               write(lupri,*)
     &            'SIR_INTOPEN_NEWTRA: NBAS_AOTWOINT .ne. NBAS:'
               write(lupri,'(8I6)') (NBAS_x(i), i=1,NSYM_x)
               write(lupri,'(8I6)') (NBAS  (i), i=1,NSYM  )
            END IF
            IF (LBUF .ne. 600) THEN
               write(lupri,*) 'LBUF .ne. 600 not implemented yet',LBUF
               ! hjaaj sep 2010: must allocate enough for larger buffer
               ! in RDAO2 outside RDAO2!
            END IF
            IF (ierr .gt. 0) THEN
               call quit(
     &            'SIR_INTOPEN_NEWTRA: AOTWOINT info not consistent')
            END IF
         END IF
      END IF
C
      IF (LUMINT .LE. 0) CALL DAOPEN(LUMINT,FNINTM)
      RETURN
      END
C  /* Deck n_nxth2m */
      SUBROUTINE N_NXTH2M(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
C
C  Written by Hans Joergen Aa. Jensen December 1989
C  This version is interface routine for new integral transformation.
C
C NOTE: The space allocated in WRK must not be touched outside
C       until all desired distributions have been read.
C
C Purpose:
C    Read next Mulliken two-electron integral distribution (**|cd)
C    where (cd) distribution is needed according to NEEDTP(ITYPCD)
C
C Usage:
C    Set IDIST = 0 before first call of NXTH2M.
C    Do NOT change IDIST or WRK(KFREE1:KFREE2-1) in calling routine
C    until last distribution has been read (signalled by IDIST .eq. -1)
C    Prototype code:
C     IDIST = 0
C     define NEEDTP(1:6)
C 100 CALL N_NXTH2M(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
C     IF (IDIST .GT. 0) THEN
C        KW1 = KFREE
C        LW1 = LFREE
C        use (**|cd) distribution in H2CD as desired
C        WRK(KW1:KW1-1+LW1) may be used
C        GO TO 100
C     END IF
C
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION H2CD(*),NEEDTP(6),WRK(*)
C
      PARAMETER (LIADRM = 3*36*36 + 3)
C     last 3 elements are used for 1) first free address,
C     2) address of CMO and 3) address of ICDTRA
C
C Used from common blocks:
C   INFORB : NNORBX
C   INFDIM : NORBMA
C
#include "inforb.h"
#include "infdim.h"
#include "infpri.h"
C
      SAVE KFRSAV, KNEXT, KIADRM, KICDTR
      DATA KNEXT /-1/
C
      CALL QENTER('N_NXTH2M')
C
C     If first read (IDIST .eq. 0) then allocate space
C     for buffers.
C
      IF (IDIST .EQ. 0) THEN
         KFRSAV = KFREE
         CALL MEMGET('INTE',KIADRM,LIADRM,WRK,KFREE,LFREE)
         CALL MEMGET('INTE',KICDTR,NNORBX,WRK,KFREE,LFREE)
         KNEXT = KFREE
      ELSE
C        ... check that IADRM has not been destroyed by calling
C            routine.
         IF (KNEXT.EQ. -1  ) THEN
            WRITE (LUERR,*)
     &         'N_NXTH2M error, IDIST must be zero in first call'
            WRITE (LUERR,*) 'IDIST =',IDIST
            CALL QTRACE(LUERR)
           CALL QUIT('N_NXTH2M error, IDIST must be zero in first call')
         END IF
         IF (KFREE.LT.KNEXT) THEN
            WRITE (LUERR,*)
     &         'N_NXTH2M error, KFREE lower than buffer allocation'
            WRITE (LUERR,*) 'KFREE ',KFREE
            WRITE (LUERR,*) 'KIADRM',KIADRM
            WRITE (LUERR,*) 'KICDTR',KICDTR
            WRITE (LUERR,*) 'KNEXT ',KNEXT,
     &         ' ( next avail. address after n_nxth2m alloc.)'
            CALL MEMCHK('KFREE.lt.KNEXT MEMCHK in N_NXTH2M',WRK,KFRSAV)
            CALL QTRACE(LUERR)
            CALL QUIT('N_NXTH2M error: KFREE .lt. buffer allocation')
         END IF
         CALL MEMCHK('IDIST .ne. 0 MEMCHK in N_NXTH2M',WRK,KFRSAV)
      END IF
C
      KREL   = KFREE
      CALL MEMGET('REAL',KBUF,2*NORBMA*NORBMA,WRK,KFREE,LFREE)
      CALL N_NX2H2M(IC,ID,H2CD,NEEDTP,IDIST,WRK(KIADRM),WRK(KICDTR),
     &            WRK(KBUF))
C
C     If finished (IDIST .lt. 0) then release buffer space
C
      IF (IDIST .LT. 0) THEN
         CALL MEMREL('Releasing all buffer space in NXTH2M',WRK,
     &               KFRSAV,KFRSAV,KFREE,LFREE)
         KNEXT = -1
      ELSE
         CALL MEMREL('Releasing BUF(NORBMA,NORBMA) in NXTH2M',WRK,
     &               KFRSAV,KREL,KFREE,LFREE)
      END IF
      CALL QEXIT('N_NXTH2M')
      RETURN
      END
C  /* Deck n_nx2h2m */
      SUBROUTINE N_NX2H2M(IC,ID,H2CD,NEEDTP,IDIST,IADRM,ICDTRA,BUF)
C
C  Written by Hans Joergen Aa. Jensen December 1989
C
C Purpose:
C
C    Read next Mulliken two-electron integral distribution (**|cd)
C    where (cd) distribution is needed according to NEEDTP(ITYPCD)
C
C Input:
C       NEEDTP(6); positive for needed (cd) distribution types
C                  negative if not all distributions needed
C                  zero if no distributions needed for this type
C       IDIST; .eq. 0 first read
C              .gt. 1 intermediate read
C              .lt. 0 end-of-file has been reached previously
C Output:
C       H2CD(NORBT,NORBT); H2CD(a,b) = (ab|cd)
C       IC,ID; value of c and d
C       IDIST; .gt. 0 when next distribution IC,ID available in H2CD
C              = -1 when no more distributions
C Scratch:
C       IADRM() for start addresses
C       ICDTRA(ICD) .ne. 0 if (**|cd) distribution has been transformed.
C       BUF(NORBMA,NORBMA,2)
C
C ****************************************************************
C
#include "implicit.h"
      DIMENSION H2CD(NORBT,NORBT),BUF(*)
      INTEGER   IADRM(3,36,36), NEEDTP(6), ICDTRA(NNORBX)
C
#include "iratdef.h"
c#include "dacodes.h"
#include "dalistgs.h"
C
C Used from common blocks:
C   INFORB : NORBX,NORBT,NSYM
C   INFIND : IROW(*)
C   INFTAP : LUMINT
C   INFPRI : IPRSIR
C
#include "maxash.h"
#include "maxorb.h"
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"
#include "infpri.h"
C
#include "orbtypdef.h"
C
      DIMENSION IADCD(8)
      SAVE      IADCD,INTSY1, JSUM, ICFRST,ICLAST,IDFRST,IDLAST
      SAVE      ICOLD,IDOLD, ISYMCO,ISYMDO,ISYMCD,ISYMAB
C
C
C ****************************************************************
C
C     If first read (IDIST .EQ. 0)
C     then setup for reading MO integrals ...
C
      IF (IDIST .EQ. 0) THEN
         ISYMCO = 1
         ISYMDO = 0
         ICFRST = 1
         ICLAST = 0
         IDFRST = 1
         IDLAST = 0
         ICOLD  = ICFRST
         IDOLD  = IDFRST - 1
C        MAERKE put this in a new NX0H2M ???
         IAD13  = 0
         LIADRM = 3*1297
         CALL DAREAD(LUMINT,IADRM,LIADRM,IAD13)
         NLIST  = 4
#if 1
         CALL DARELIST(LUMINT,IAD13,
     &                 NLIST,NTRLV1,1,INTSY1,1,MSYM,1,LCMO,1)
#else
         CALL GSLIST(LISTGS,NLIST,NTRLV1,1,INTSY1,1,MSYM,1,LCMO,1)
C        NLIST = 10
C        CALL GSLIST(LISTGS,NLIST,NTRLVL,1,INTSYM,1,MSYM,1,LCMO,1,
C    &               MFRO,8,MISH,8,MASH,8,MORB,8,MBAS,8,ITRTYP,MXCORB)
         CALL DAFILE(LUMINT,DASREA,LISTGS,NLIST,IAD13)
#endif
         CALL TRAADR(IADRM,IADFRE,IADCMO,IADCDT)
C        ... retrieve first free address, CMO address,
C            and ICDTRA address
         IAD13  = IADCDT
         CALL DAREAD(LUMINT,ICDTRA,NNORBX,IAD13)
         IF (IPRSIR .GT. 50) THEN
            WRITE (LUPRI,*) 'Test output from NX2H2M for IDIST = 0'
            WRITE (LUPRI,*) 'NTRLVL from LUMINT :',NTRLV1
            WRITE (LUPRI,*) 'NSYM   from LUMINT :',MSYM
            WRITE (LUPRI,*) 'INTSYM from LUMINT :',INTSY1
            WRITE (LUPRI,*) 'ICDTRA matrix:'
            DO 10 I = 1,NORBT
               IR = IROW(I)
               WRITE (LUPRI,'(I5,(T8,10I7))') I,(ICDTRA(IR+J),J=1,I)
   10       CONTINUE
            WRITE (LUPRI,*) 'Positive addresses in IADRM:'
            WRITE (LUPRI,*)
     &      ' index(sympq) index(symrs) IADC  IADE1 IADE2'
            DO 20 J = 1,36
            DO 20 I = 1,36
               IF (IADRM(1,I,J) .GT. 0 .OR.
     &             IADRM(2,I,J) .GT. 0 .OR.
     &             IADRM(3,I,J) .GT. 0 ) THEN
                  WRITE (LUPRI,*) I,J,(IADRM(K,I,J),K=1,3)
               END IF
   20       CONTINUE
            WRITE (LUPRI,*) 'Address for ICDTRA',IADCDT
            WRITE (LUPRI,*) 'Address for CMO   ',IADCMO
            WRITE (LUPRI,*) 'First free address',IADFRE
         END IF
      END IF
C
C *** Initialize H2CD
C
      CALL DZERO(H2CD,N2ORBX)
C
C *** Read next distribution which is needed according to NEEDTP(6)
C     into H2CD
C
C  ITYPCD values:  1=i*i :  2=t*i : 3=t*t : 4=a*i : 5=a*t : 6=a*a
C                  0 for not wanted type.
C
C  The CD distributions are stored by the present transformation
C  program with IC.ge.ID
C
      ICNEW = ICOLD
      IDNEW = IDOLD
  200 CONTINUE
      IDNEW = IDNEW + 1
      IF (IDNEW .GT. IDLAST) THEN
         ICNEW = ICNEW + 1
         IDNEW = IDFRST
         IF (ISYMCO .EQ. ISYMDO) IDLAST = ICNEW
      END IF
      IF (ICNEW .GT. ICLAST) THEN
C        This symmetry block is finished, change symmetry block:
         ISYMDO = ISYMDO + 1
         IF (ISYMDO .GT. ISYMCO) THEN
            ISYMCO = ISYMCO + 1
            ISYMDO = 1
         END IF
         IF (ISYMCO .GT. NSYM) THEN
C           Last distribution has been read
            IDIST = -1
            GO TO 9999
         END IF
C
         ICFRST = IORB(ISYMCO) + 1
         ICLAST = IORB(ISYMCO) + NORB(ISYMCO)
         IDFRST = IORB(ISYMDO) + 1
         IDLAST = IORB(ISYMDO) + NORB(ISYMDO)
         ICNEW  = ICFRST
         IDNEW  = IDFRST - 1
         IF (ISYMCO .EQ. ISYMDO) IDLAST = ICNEW
C
         ISYMCD = MULD2H(ISYMCO,ISYMDO)
         ISYMAB = MULD2H(ISYMCD,INTSY1)
         NSYMCD = IROW(ISYMCO) + ISYMDO
         JSUM   = 0
         DO 300 ISYMA = 1,NSYM
            ISYMB = MULD2H(ISYMA,ISYMAB)
            IF (ISYMB .GT. ISYMA) THEN
               IADCD(ISYMA) = 0
            ELSE
               NSYMAB = IROW(ISYMA) + ISYMB
               IADCD(ISYMA) = IADRM(1,NSYMAB,NSYMCD)
               JSUM = JSUM + MAX(IADCD(ISYMA),0)
            END IF
  300    CONTINUE
         IF (IPRSIR .GT. 50) THEN
            WRITE(LUPRI,*) 'New ISYMC and ISYMD :',ISYMCO,ISYMDO
            WRITE(LUPRI,*) 'IADCD(isyma) :',(IADCD(ISYMA),ISYMA=1,NSYM)
            WRITE(LUPRI,*) 'JSUM =',JSUM
         END IF
         GO TO 200
      END IF
C
      ICDNEW = IROW(ICNEW) + IDNEW
      IDIST  = IDIST + 1
      ITYPC  = IOBTYP(ICNEW)
      ITYPD  = IOBTYP(IDNEW)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF ( NEEDTP(ITYPCD) .EQ. 0 ) THEN
         ITYPCD = 0
      ELSE IF (NEEDTP(ITYPCD) .LT. 0) THEN
         ITYPCD = -ITYPCD
      END IF
C
      IF (IPRSIR .GT. 50) THEN
         WRITE(LUPRI,*) 'ICNEW ,IDNEW   :',ICNEW,IDNEW
         WRITE(LUPRI,*) 'ICDTRA(ICDNEW) :',ICDTRA(ICDNEW)
         WRITE(LUPRI,*) 'IDIST ,ITYPCD  :',IDIST,ITYPCD
         WRITE(LUPRI,*) 'IADCD(isyma) :',(IADCD(ISYMA),ISYMA=1,NSYM)
      END IF
C
      IF (ITYPCD .GT. 0 .AND. ICDTRA(ICDNEW) .EQ. 0) THEN
         WRITE(LUPRI,'(//A/)') 'NX2H2M ERROR: '//
     &     'needed integral distribution has not been calculated'
         WRITE(LUPRI,*) 'IC    ,ID     :',ICNEW,IDNEW
         WRITE(LUPRI,*) 'IYPC  ,ITYPD  : ',COBTYP(ITYPC),
     &                                '  ',COBTYP(ITYPD)
         WRITE (LUPRI,'(/A)') 'ICDTRA matrix:'
         DO I = 1,NORBT
            IR = IROW(I)
            WRITE (LUPRI,'(I5,(T8,10I7))') I,(ICDTRA(IR+J),J=1,I)
         END DO
         CALL QUIT('N_NX2H2M error: needed integrals not calculated')
      END IF
C
      IF (ICDTRA(ICDNEW) .EQ. 0 .OR. JSUM .EQ. 0) GO TO 200
C     If JSUM = 0 then no integrals in /csym dsym) symmetry block
C
      NINT = 0
      DO 400 ISYMA = 1,NSYM
         ISYMB = MULD2H(ISYMA,ISYMAB)
         IF (ISYMB .GT. ISYMA) THEN
            LREC = 0
         ELSE IF (ISYMB .EQ. ISYMA) THEN
            NORBA = NORB(ISYMA)
            NORBB = NORBA
            IA1   = IORB(ISYMA) + 1
            IB1   = IA1
            LREC  = IRAT * IROW( NORBA + 1 )
            IBUF1 = 1 + NORBA*NORBA
            IBUF2 = 1
         ELSE
            NORBA = NORB(ISYMA)
            NORBB = NORB(ISYMB)
            IA1   = IORB(ISYMA) + 1
            IB1   = IORB(ISYMB) + 1
            LREC  = IRAT * NORBA * NORBB
            IBUF1 = 1
            IBUF2 = 1
         END IF
         IF (LREC .GT. 0) THEN
            IF (IADCD(ISYMA) .LT. 0) THEN
               CALL QTRACE(LUERR)
               CALL QUIT('NX2H2M ERROR, IADCD(ISYMA) .lt. 0')
            END IF
            IF (ITYPCD.EQ.0) THEN
C           ... skip this record
               CALL DASKIP(LUMINT,LREC,IADCD(ISYMA))
            ELSE
               NINT = NINT + 1
               CALL DAREAD(LUMINT,BUF(IBUF2),LREC,IADCD(ISYMA))
C              ... BUF(b,a) now contains (ab/cd)
               IF (ISYMB .EQ. ISYMA) THEN
                  CALL DSPTSI(NORBA,BUF(IBUF2),BUF(IBUF1))
                  CALL MCOPY(NORBA,NORBA,BUF(IBUF1),NORBA,
     &                       H2CD(IA1,IA1),NORBT)
               ELSE
                  CALL MCOPY(NORBB,NORBA,BUF(IBUF1),NORBB,
     &                       H2CD(IB1,IA1),NORBT)
                  CALL MTRSP(NORBB,NORBA,BUF(IBUF1),NORBB,
     &                       H2CD(IA1,IB1),NORBT)
               END IF
            END IF
         END IF
  400 CONTINUE
      IF (NINT.EQ.0) GO TO 200
C
C*******************************************************************
C
C End of subroutine NX2H2M
C
 9999 CONTINUE
      ICOLD  = ICNEW
      IDOLD  = IDNEW
      IC     = ICNEW
      ID     = IDNEW
      RETURN
      END
C  /* Deck n_nxth2d */
      SUBROUTINE N_NXTH2D(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
C
C  Written by Hans Joergen Aa. Jensen March 1990
C  This version is interface routine for new integral transformation.
C
C NOTE: The space allocated in WRK must not be touched outside
C       until all desired distributions have been read.
C
C Purpose:
C    Read next Dirac two-electron integral distribution <**|cd>
C    where (cd) distribution is needed according to NEEDTP(ITYPCD)
C
C Usage:
C    Set IDIST = 0 before first call of NXTH2D.
C    DO NOT CHANGE IDIST or WRK(KFREE1:KFREE2-1) in calling routine
C    until last distribution has been read (signalled by IDIST .eq. -1)
C    Prototype code:
C     IDIST = 0
C     define NEEDTP(1:6)
C 100 CALL N_NXTH2D(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
C     IF (IDIST .GT. 0) THEN
C        KW1 = KFREE
C        LW1 = LFREE
C        use <**|cd> distribution in H2CD as desired
C        WRK(KW1:KW1-1+LW1) may be used
C        GO TO 100
C     END IF
C
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION H2CD(*),NEEDTP(6),WRK(*)
C
      PARAMETER (LIADRM = 3*36*36 + 3)
C     last 3 elements are used for 1) first free address,
C     2) address of CMO and 3) address of ICDTRA
C
C Used from common blocks:
C   INFORB : NNORBX
C   INFDIM : NORBMA
C
#include "inforb.h"
#include "infdim.h"
#include "infpri.h"
C
      SAVE KFRSAV, KNEXT, KIADRM, KICDTR
      DATA KNEXT /-1/
C
      CALL QENTER('N_NXTH2D')
C
C     If first read (IDIST .eq. 0) then allocate space
C     for buffers.
C
      IF (IDIST .EQ. 0) THEN
         KFRSAV = KFREE
         CALL MEMGET('INTE',KIADRM,LIADRM,WRK,KFREE,LFREE)
         CALL MEMGET('INTE',KICDTR,NNORBX,WRK,KFREE,LFREE)
         KNEXT = KFREE
      ELSE
C        ... check that IADRM has not been destroyed by calling
C            routine.
         IF (KNEXT.EQ. -1  ) THEN
            WRITE (LUERR,*)
     &         'NXTH2D error, IDIST must be zero in first call'
            WRITE (LUERR,*) 'IDIST =',IDIST
            CALL QTRACE(LUERR)
           CALL QUIT('N_NXTH2D error, IDIST must be zero in first call')
         END IF
         IF (KFREE.LT.KNEXT) THEN
            WRITE (LUERR,*)
     &         'NXTH2D error, KFREE lower than buffer allocation'
            WRITE (LUERR,*) 'KFREE ',KFREE
            WRITE (LUERR,*) 'KIADRM',KIADRM
            WRITE (LUERR,*) 'KICDTR',KICDTR
            WRITE (LUERR,*) 'KNEXT ',KNEXT,
     &         ' ( next avail. address after nxth2d alloc.)'
            CALL MEMCHK('KFREE.lt.KNEXT MEMCHK in NXTH2D',WRK,KFRSAV)
            CALL QTRACE(LUERR)
            CALL QUIT('NXTH2D error: KFREE inside buffer allocation')
         END IF
      END IF
C
      KREL   = KFREE
      CALL MEMGET('REAL',KBUF,NORBMA*NORBMA,WRK,KFREE,LFREE)
      CALL N_NX2H2D(IC,ID,H2CD,NEEDTP,IDIST,WRK(KIADRM),WRK(KICDTR),
     &            WRK(KBUF))
C
C     If finished (IDIST .lt. 0) then release buffer space
C
      IF (IDIST .LT. 0) THEN
         CALL MEMREL('Releasing all buffer space in NXTH2D',WRK,
     &               KFRSAV,KFRSAV,KFREE,LFREE)
         KNEXT = -1
      ELSE
         CALL MEMREL('Releasing BUF(NORBMA,NORBMA) in NXTH2D',WRK,
     &               KFRSAV,KREL,KFREE,LFREE)
      END IF
      CALL QEXIT('N_NXTH2D')
      RETURN
      END
C  /* Deck n_nx2h2d */
      SUBROUTINE N_NX2H2D(IC,ID,H2CD,NEEDTP,IDIST,IADRM,ICDTRA,BUF)
C
C  Written by Hans Joergen Aa. Jensen March 1990
C
C Purpose:
C
C    Read next Dirac two-electron integral distribution <**|cd>
C    where (cd) distribution is needed according to NEEDTP(ITYPCD)
C
C Input:
C       NEEDTP(6); positive for needed (cd) distribution types
C                  negative if not all distributions needed
C                  zero if no distributions needed for this type
C       IDIST; .eq. 0 first read
C              .gt. 1 intermediate read
C              .lt. 0 end-of-file has been reached previously
C Output:
C       H2CD(NORBT,NORBT); H2CD(a,b) = <ab|cd>
C       IC,ID; value of c and d
C       IDIST; .gt. 0 when next distribution IC,ID available in H2CD
C              = -1 when no more distributions
C Scratch:
C       IADRM() for start addresses
C       ICDTRA(ICD) .ne. 0 if <**|cd> distribution has been transformed.
C       BUF(NORBMA,NORBMA,2)
C
C ****************************************************************
C
#include "implicit.h"
      DIMENSION H2CD(NORBT,NORBT),BUF(*)
      INTEGER   IADRM(3,36,36), NEEDTP(6), ICDTRA(NNORBX)
C
#include "iratdef.h"
c#include "dacodes.h"
#include "dalistgs.h"
C
C Used from common blocks:
C   INFORB : NORBX,NORBT,NSYM
C   INFIND : IROW(*)
C   INFTAP : LUMINT
C   INFPRI : IPRSIR
C
#include "maxash.h"
#include "maxorb.h"
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"
#include "infpri.h"
C
#include "orbtypdef.h"
C
      DIMENSION IADCD1(8),IADCD2(8)
      SAVE      IADCD1,IADCD2,INTSY1, JSUM, ICFRST,ICLAST,IDFRST,IDLAST
      SAVE      ICOLD,IDOLD, ISYMCO,ISYMDO,ISYMCD,ISYMAB
C
C
C ****************************************************************
C
C     If first read (IDIST .EQ. 0)
C     then setup for reading Dirac MO integrals ...
C
      IF (IDIST .EQ. 0) THEN
         ISYMCO = 1
         ISYMDO = 0
         ICFRST = 1
         ICLAST = 0
         IDFRST = 1
         IDLAST = 0
         ICOLD  = ICFRST
         IDOLD  = IDFRST - 1
C        MAERKE put this in a new NX0H2D ???
         IAD13  = 0
         LIADRM = 3*1297
         CALL DAREAD(LUMINT,IADRM,LIADRM,IAD13)
         NLIST  = 4
#if 1
         CALL DARELIST(LUMINT,IAD13,
     &                 NLIST,NTRLV1,1,INTSY1,1,MSYM,1,LCMO,1)
#else
         CALL GSLIST(LISTGS,NLIST,NTRLV1,1,INTSY1,1,MSYM,1,LCMO,1)
C        NLIST = 10
C        CALL GSLIST(LISTGS,NLIST,NTRLVL,1,INTSYM,1,MSYM,1,LCMO,1,
C    &               MFRO,8,MISH,8,MASH,8,MORB,8,MBAS,8,ITRTYP,MXCORB)
         CALL DAFILE(LUMINT,DASREA,LISTGS,NLIST,IAD13)
#endif
         CALL TRAADR(IADRM,IADFRE,IADCMO,IADCDT)
C        ... retrieve ICDTRA address in IADCDT
         IAD13  = IADCDT
         CALL DAREAD(LUMINT,ICDTRA,NNORBX,IAD13)
         IF (IPRSIR .GT. 50) THEN
            WRITE (LUPRI,*) 'Test output from NX2H2D for IDIST = 0'
            WRITE (LUPRI,*) 'NTRLVL from LUMINT :',NTRLV1
            WRITE (LUPRI,*) 'NSYM   from LUMINT :',MSYM
            WRITE (LUPRI,*) 'INTSYM from LUMINT :',INTSY1
            WRITE (LUPRI,*) 'ICDTRA matrix:'
            DO 10 I = 1,NORBT
               IR = IROW(I)
               WRITE (LUPRI,'(I5,(T8,10I7))') I,(ICDTRA(IR+J),J=1,I)
   10       CONTINUE
            WRITE (LUPRI,*) 'Positive addresses in IADRM:'
            WRITE (LUPRI,*)
     &      ' index(sympq) index(symrs) IADC  IADE1 IADE2'
            DO 20 J = 1,36
            DO 20 I = 1,36
               IF (IADRM(1,I,J) .GT. 0 .OR.
     &             IADRM(2,I,J) .GT. 0 .OR.
     &             IADRM(3,I,J) .GT. 0 ) THEN
                  WRITE (LUPRI,*) I,J,(IADRM(K,I,J),K=1,3)
               END IF
   20       CONTINUE
            WRITE (LUPRI,*) 'Address for ICDTRA',IADCDT
            WRITE (LUPRI,*) 'Address for CMO   ',IADCMO
            WRITE (LUPRI,*) 'First free address',IADFRE
         END IF
      END IF
C
C *** Initialize H2CD
C
      CALL DZERO(H2CD,N2ORBX)
C
C *** Read next distribution which is needed according to NEEDTP(6)
C     into H2CD
C
C  ITYPCD values:  1=i*i :  2=t*i : 3=t*t : 4=a*i : 5=a*t : 6=a*a
C                  0 for not wanted type.
C
C  The CD distributions are stored by the present transformation
C  program with IC.ge.ID
C
      ICNEW = ICOLD
      IDNEW = IDOLD
  200 CONTINUE
      IDNEW = IDNEW + 1
      IF (IDNEW .GT. IDLAST) THEN
         ICNEW = ICNEW + 1
         IDNEW = IDFRST
         IF (ISYMCO .EQ. ISYMDO) IDLAST = ICNEW
      END IF
      IF (ICNEW .GT. ICLAST) THEN
C        This symmetry block is finished, change symmetry block:
         ISYMDO = ISYMDO + 1
         IF (ISYMDO .GT. ISYMCO) THEN
            ISYMCO = ISYMCO + 1
            ISYMDO = 1
         END IF
         IF (ISYMCO .GT. NSYM) THEN
C           Last distribution has been read
            IDIST = -1
            GO TO 9999
         END IF
C
         ICFRST = IORB(ISYMCO) + 1
         ICLAST = IORB(ISYMCO) + NORB(ISYMCO)
         IDFRST = IORB(ISYMDO) + 1
         IDLAST = IORB(ISYMDO) + NORB(ISYMDO)
         ICNEW  = ICFRST
         IDNEW  = IDFRST - 1
         IF (ISYMCO .EQ. ISYMDO) IDLAST = ICNEW
C
         ISYMCD = MULD2H(ISYMCO,ISYMDO)
         ISYMAB = MULD2H(ISYMCD,INTSY1)
         NSYMCD = IROW(ISYMCO) + ISYMDO
         JSUM   = 0
         DO 300 ISYMA = 1,NSYM
            ISYMB = MULD2H(ISYMA,ISYMAB)
            IF (ISYMB .GT. ISYMA) THEN
               IADCD1(ISYMA) = 0
               IADCD2(ISYMA) = 0
            ELSE
               NSYMAB = IROW(ISYMA) + ISYMB
               IADCD1(ISYMA) = IADRM(2,NSYMAB,NSYMCD)
               IADCD2(ISYMA) = IADRM(3,NSYMAB,NSYMCD)
               JSUM = JSUM + MAX(IADCD1(ISYMA),IADCD2(ISYMA),0)
            END IF
  300    CONTINUE
         IF (IPRSIR .GT. 50) THEN
            WRITE(LUPRI,*) 'Test output from NX2H2D'
            WRITE(LUPRI,*) 'New ISYMC and ISYMD :',ISYMCO,ISYMDO
            WRITE(LUPRI,*) '   ISYMA    IADCD1(isyma) IADCD2(isyma):'
            DO 310 ISYMA =1,NSYM
  310          WRITE (LUPRI,*) ISYMA,IADCD1(ISYMA),IADCD2(ISYMA)
            WRITE(LUPRI,*) 'JSUM =',JSUM
         END IF
         GO TO 200
      END IF
C
      ICDNEW = IROW(ICNEW) + IDNEW
      IDIST  = IDIST + 1
      ITYPC  = IOBTYP(ICNEW)
      ITYPD  = IOBTYP(IDNEW)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF ( NEEDTP(ITYPCD) .EQ. 0 ) THEN
         ITYPCD = 0
      ELSE IF (NEEDTP(ITYPCD) .LT. 0) THEN
         ITYPCD = -ITYPCD
      END IF
C
      IF (IPRSIR .GT. 50) THEN
         WRITE(LUPRI,*) 'Test output from N_NX2H2D'
         WRITE(LUPRI,*) 'ICNEW ,IDNEW   :',ICNEW,IDNEW
         WRITE(LUPRI,*) 'ICDTRA(ICDNEW) :',ICDTRA(ICDNEW)
         WRITE(LUPRI,*) 'IDIST ,ITYPCD  :',IDIST,ITYPCD
         WRITE(LUPRI,*) '   ISYMA    IADCD1(isyma) IADCD2(isyma):'
         DO 320 ISYMA =1,NSYM
  320       WRITE (LUPRI,*) ISYMA,IADCD1(ISYMA),IADCD2(ISYMA)
      END IF
C
      IF (ITYPCD .GT. 0 .AND. ICDTRA(ICDNEW) .EQ. 0) THEN
         WRITE (LUERR,*) ' N_NX2H2D ERROR: needed integral distribution'
         WRITE (LUERR,*) '               has not been calculated'
         WRITE (LUERR,*) 'IC    ,ID     :',ICNEW,IDNEW
         WRITE (LUERR,*) 'IYPC  ,ITYPD  :',COBTYP(ITYPC),COBTYP(ITYPD)
         CALL QTRACE(LUERR)
         CALL QUIT('N_NX2H2D error: needed integrals not calculated')
      END IF
C
      IF (ICDTRA(ICDNEW) .EQ. 0 .OR. JSUM .EQ. 0) GO TO 200
C     If JSUM = 0 then no integrals in /csym dsym) symmetry block
C
      NINT = 0
      DO 400 ISYMA = 1,NSYM
         ISYMB = MULD2H(ISYMA,ISYMAB)
         IF (ISYMB .GT. ISYMA) THEN
            LREC = 0
         ELSE
            NORBA = NORB(ISYMA)
            NORBB = NORB(ISYMB)
            IA1   = IORB(ISYMA) + 1
            IB1   = IORB(ISYMB) + 1
            LREC  = IRAT * NORBA * NORBB
            IBUF1 = 1
         END IF
         IF (LREC .GT. 0) THEN
            IF (IADCD1(ISYMA) .LT. 0 .AND. IADCD2(ISYMA) .LT. 0) THEN
               CALL QTRACE(LUERR)
               CALL QUIT('N_NX2H2D ERROR, IADCD*(ISYMA) .lt. 0')
            END IF
            IF (ITYPCD.EQ.0) THEN
C           ... skip this record
               IF (IADCD1(ISYMA) .GT. 0)
     &            CALL DASKIP(LUMINT,LREC,IADCD1(ISYMA))
               IF (IADCD2(ISYMA) .GT. 0)
     &            CALL DASKIP(LUMINT,LREC,IADCD2(ISYMA))
            ELSE
               NINT = NINT + 1
               IF (IADCD1(ISYMA) .GT. 0) THEN
                  CALL DAREAD(LUMINT,BUF(IBUF1),LREC,IADCD1(ISYMA))
C                 ... BUF(b,a) now contains <ab/cd>
                  CALL MTRSP(NORBB,NORBA,BUF(IBUF1),NORBB,
     &                       H2CD(IA1,IB1),NORBT)
               END IF
               IF (IADCD2(ISYMA) .GT. 0) THEN
                  CALL DAREAD(LUMINT,BUF(IBUF1),LREC,IADCD2(ISYMA))
C                 ... BUF(a,b) now contains <ba/cd>
                  CALL MCOPY(NORBB,NORBA,BUF(IBUF1),NORBB,
     &                       H2CD(IB1,IA1),NORBT)
               END IF
            END IF
         END IF
  400 CONTINUE
      IF (NINT.EQ.0) GO TO 200
C
C*******************************************************************
C
C End of subroutine NX2H2D
C
 9999 CONTINUE
      ICOLD  = ICNEW
      IDOLD  = IDNEW
      IC     = ICNEW
      ID     = IDNEW
      RETURN
      END
      SUBROUTINE N_TRAUTH(IWUNIT)
C
C     dec 89-hjaaj
C     Author of two-electron integral transformation routines.
C
      WRITE(IWUNIT,100)
  100 FORMAT(/T2,'Two-electron integral transformation:'/
     *   T5,'Byorn Roos,              University of Lund,    Sweden'/,
     *   T5,'Hans Joergen Aa. Jensen, University of Odense,  Denmark')
      RETURN
      END
#ifdef UNDEF
C* Comdeck notes and ideas*
!se MAERKE
!ABACUS new transf: soerg for alt i TRASET, TRALIM, TRAINI er i orden!
!900308-hjaaj: fjerne MCOPT og altid blot kraeve mindst to okkup.?
! (forudsat ITRLVL .le. 2); eller teste paa INTSYM? (er det altid ok?).
!900124-hjaaj: indfoert ITUSYM parameter for spin-orbit integraler!
!891227-hjaaj: new input: *TRAINP .PRINT .FIOPRI .PRIUNIT .SKIP
!  .OLDMOLECULE .NOSUPMAT .OPTIONS ...
!  - REMEMBER change abadtrn.u if change in tr2ctl (everything defined
!  in TRACTL, TRASET and TRALIM must be define in abadtrn.u)
!  HUSK IDISKM in IADOUT(*,LIADUT)
!  - LIADUT = 1297 = 1296 + 1 where IADOUT(*,1297) is used for next
!  free address
!891220-Henrik: infoer kun Dirac til CCSD.
!890825-hjaaj: if out-of-core nec. then try to keep two or at least one
!  of TUPQ, RUPQ, USPQ in core instead of taking all three out-of-core.
!  - changed NBAS to MBAS etc.
!890704-hjaaj: option for appending exchange integral to previously
!   calculated Coulomb integrals
!890612-hjaaj TRA2 call of GTUVX does not work if MISH.gt.0 or MFRO.gt.0
!890611-hjaaj CHECK2 does not work if MISH.gt.0 or MFRO.gt.0
!890605-hjaaj PROGRAM MOTRA2 gives conflict with NAMELIST MOTRA2
!/* Comdeck job_log */
!900314-hjaaj: - written code in TRALIM to determine transf. from ITRLVL
!900313-hjaaj: - implemented JMO[PQRS]2 and ICDTRA in TR2CTL and TRA2
!  - major overhaul of TRALIM (modify MFRO,MISH,MASH,
!  define ITRTYP, define ICDTRA)
!900124-hjaaj: - ITUSYM in TRINP cdk; define ITUSYM in TRASET;
!  s/SQUARE/TRASQ/; added IASYM parameter to TRASQ and calls of TRASQ
!900102-hjaaj: - TRAINI: added CMO to param. and LUMINT
!  - TRALIM : defined ICDTRA; - new option NOCOUL;
!  - TRA2: new local var. DORUPQ, DOUSPQ, DOTUPQ
!891229-hjaaj: - revised symmetry checks for exhange integrals so
!  INTSYM.ne.1 cases also can be handled.
!891228-hjaaj: - changed code in TR2CTL to minimize number of temporary
!  files (LUHLF*); - new routine TRAINI
!  - split TRA2 in TRA2 and TRA2E1
!891227-hjaaj: - new routine TRASET
!  - changed test program motra2 to define inforb,inftap(used in TRASET)
!  - GTUVX and allocation of TUVX : NNASHX*NNASHX instead of NNASHY
!891221-hjaaj: - RDDER2 dummy routine
!  - Reordered decks in two groups: newtractl.u and chktra2.u
!  - Worked on deck tractl, was unfinished. Introduced ICDTRA(NNORBX)
!891220-hjaaj: transfer INTSYM and ICDTRA(NNORBX) to LUMINT
!  - #include "dacodes.h"and use symbolic names for DAFILE operations.
!  - #include "dalistgs.h"for gslist.
!890826-hjaaj
!   moved IDATA(352),KEEP(18) to TRUNIT (was local to TR2CTL)
!   changed IADOUT(3888) to IADOUT(3,LIADUT), LIADUT = 1296
!   new IDATA(LIDATA), LIDATA = 1296=36*36 for absolute addressing
!   changed NSYMP from NSYM*(NSYM+1)/2 to 36 = 8*(8+1)/2 for
!     absolute addressing, defined NSYMP in TRINFO
!   changed NCHAIN in TR2CTL to absolute addressing
!890825-hjaaj
!   use local DOEXCH for greater readability in TRA2
!   changed *URPQ,NOUR into *USPQ,NOUS in /TRINFO/,TR2CTL,TRA2
!890824-hjaaj
!   New routine TRALIM; called from TR2CTL
!890823-hjaaj-t13
!   Replaced IEOR with IBTXOR statement function
!890705-hjaaj-t13
!   Skip memory tests for exchange integrals if NOEXCH
!================
!890620-hjaaj-t12.1
!   NOEXCH, MCOPT, and MAXMEM options
!   Test that at least two symmetries occupied in TR2CTL
!================
!890619-hjaaj-t11.2
!   NOTUVX option
!   Check that MEMX .gt. 0
!================
!890610-hjaaj-t10.2
!   Corrected error in CHECK2 (calculation of PQRS address)
!   Count non-zero integrals in test in CHECK2
!   Check that LADX ne 0 in TR2CTL
!890609-hjaaj-t10.1
!   New module GETH2A for reading AO integrals from Molecule LUINTA
!     used when VECTOR=UNIT
!   INTSYM in common TRINP, use INTSYM in CHECK/CHECK2
!================
!890605-hjaaj-t09.1
!   Introduced ITRLVL in parameter lists, value not used yet.
!   Check for INTSYM from LUINTA in KEEP(2)
!   call darecl(LUINTA,8*8*1024) if *DEFINE DARECL
!   Introduced FASTIO test output, starting w. STATUS at IPRTRA=10
!================
#endif
C  /* Deck tr2ctl */
      SUBROUTINE TR2CTL(NOTUVX,NOCOUL,NOEXCH,MCOPT,CMO,TUVX,
     &                  ICDTRA,WORK,MEMX)
C
C     SECOND ORDER TWO-ELECTRON TRANFORMATION PROGRAM. CONTROL SECTION
C
C     THIS SUBROUTINE SETS UP THE MEMORY ALLOCATIONS FOR TRA2 AND LOOPS
C     OVER THE SYMMETRY BLOCKS. TRA2 IS CALLED ONCE FOR EACH SYMMETRY
C     BLOCK OF INTEGRALS. SYMMETRY BLOCKED AO INTEGRALS MUST HAVE BEEN
C     GENERATED BY INTSORT ON UNIT LUORDA.
C
C    written in Garching in september 1987
C    author: Byoern Roos
C            Department of Theoretical Chemistry
C            Chemical Centre
C            P.O.B. 124
C            S-221 00 Lund Sweden     tel: 46-10 82 51
C
C    Last revision 89 06 05 Hans Joergen Aa. Jensen
C    Revisions:
C       INTSYM parameter (890106,890605/hjaaj)
C       ********** IBM-3090 RELELASE 87 09 14 **********
C
C    Comments:
C       hj-890106: INTSYM = symmetry of integrals
C       In revised INTSORT KEEP(2) = INTSYM below
C
#include "implicit.h"
      DIMENSION CMO(*),TUVX(*),WORK(*),ICDTRA(*)
      LOGICAL   NOTUVX, NOCOUL, NOEXCH, MCOPT
C
c#include "dacodes.h"
C
C Used from common blocks:
C  TRINP  : INTSYM
C  TRUNIT : LIDATA,IDATA,KEEP, LIADUT,IADOUT, ?
C  TRINFO : LX2X3, ?
C  INFTRA : IPRTRA
C  INFTAP : LUMINT
C
#include "priunit.h"
#include "trinp.h"
#include "trunit.h"
#include "trinfo.h"
#include "inftra.h"
#include "inftap.h"
C
#include "ibtfun.h"
C
      CALL QENTER('TR2CTL')
C
C
C     SET TIME AT START OF TRANSFORMATION
C
      IF (IPRTRA .GT. 0) CALL SETTIM
C
      IF (IPRTRA .GE. 2) THEN
         IF (NOCOUL) THEN
            WRITE (LUPRI,'(/A)') ' Coulomb  integrals not calculated.'
         ELSE
            WRITE (LUPRI,'(/A)') ' Coulomb  integrals calculated.'
         END IF
         IF (NOEXCH) THEN
            WRITE (LUPRI,'( A)') ' Exchange integrals not calculated.'
C            ' Exchange (Dirac) integrals not calculated.'
         ELSE
            WRITE (LUPRI,'( A)') ' Exchange integrals calculated.'
         END IF
         IF (.NOT.NOTUVX) THEN
            WRITE (LUPRI,'(/A)') ' H2AC extracted.'
         END IF
         IF (MCOPT) THEN
            WRITE (LUPRI,'(/A)')
     &         ' MCSCF optimization: only integrals with all 4'/
     &        /' symmetries occupied needed.'
         END IF
      END IF
      L2M = 0
C
C     READ THE DAFILE ADDRESS LIST ON UNIT LUMINT
C
      IDISKM = 0
      CALL DAREAD(LUMINT,IADOUT,3*LIADUT,IDISKM)
      IDISKM = IADOUT(1,LIADUT)
C
C     Loop over quadruples of symmetries (NSP,NSQ,NSR,NSS)
C     NOTE that the integrals on LUORDA must be sorted in the same
C     order as the loop structure below (use program INTSORT)
C     890826: This is not neccessary now because we have switched
C             to absolute addressing of a given quadruplet.
C
      IF (IPRTRA .GE. 2 .AND. IPRTRA .LT. 5) WRITE(LUPRI,2000)
#if defined (SYS_UNIX) || defined (SYS_CRAY) || defined (SYS_AIX) || defined (SYS_T3D) || defined (SYS_HAL) || defined (SYS_T90)
 2000 FORMAT(/' SYMMETRY       BASIS FUNCTION        ORBITALS',3X,
     &       ' MEMLFT     LTUPQ   USER(SEC) SYS(SEC)')
#else
 2000 FORMAT(/' SYMMETRY       BASIS FUNCTION        ORBITALS',3X,
     &       ' MEMLFT     LTUPQ    CPU(SEC) I/O(SEC)')
#endif
C
      ITP    = 0
      LMOP1  = 1
      JMOP1  = 0
      DO 104 NSP = 1,MSYM
       IF(NSP.NE.1) THEN
         ITP   = ITP   +MASH(NSP-1)
         LMOP1 = LMOP1 +MBAS(NSP-1)*MORB(NSP-1)
         JMOP1 = JMOP1 +MORB(NSP-1)
       END IF
       NBP   = MBAS(NSP)
       LMOP  = LMOP1 +NBP*MFRO(NSP)
       LMOP2 = LMOP  +NBP*MISH(NSP)
       JMOP2 = JMOP1 +MFRO(NSP)+MISH(NSP)
       NOP   = MORB(NSP)
       NOCP  = MASH(NSP)
       KEEPP = KEEP(10+NSP)
       ISP   = NSP
       ITQ   = 0
       LMOQ1 = 1
       JMOQ1 = 0
       DO 103 NSQ = 1,NSP
        IF(NSQ.NE.1) THEN
          ITQ   = ITQ+MASH(NSQ-1)
          LMOQ1 = LMOQ1+MBAS(NSQ-1)*MORB(NSQ-1)
          JMOQ1 = JMOQ1+MORB(NSQ-1)
        END IF
        NBQ   = MBAS(NSQ)
        LMOQ  = LMOQ1+NBQ*MFRO(NSQ)
        LMOQ2 = LMOQ +NBQ*MISH(NSQ)
        JMOQ2 = JMOQ1+MFRO(NSQ)+MISH(NSQ)
        KEEPQ = KEEP(10+NSQ)
        NOQ   = MORB(NSQ)
        NOCQ  = MASH(NSQ)
        ISQ   = NSQ
        NSPQ  = IBTXOR(NSP-1,NSQ-1)+1
        ITR   = 0
        LMOR1 = 1
        JMOR1 = 0
        DO 102 NSR = 1,MSYM
         IF(NSR.NE.1) THEN
           ITR   = ITR+MASH(NSR-1)
           LMOR1 = LMOR1+MBAS(NSR-1)*MORB(NSR-1)
           JMOR1 = JMOR1+MORB(NSR-1)
         END IF
         NBR   = MBAS(NSR)
         LMOR  = LMOR1+NBR*MFRO(NSR)
         LMOR2 = LMOR +NBR*MISH(NSR)
         JMOR2 = JMOR1+MFRO(NSR)+MISH(NSR)
         KEEPR = KEEP(10+NSR)
         NOR   = MORB(NSR)
         NOCR  = MASH(NSR)
         NSPQR = IBTXOR(NSPQ-1,NSR-1)+1
         ISR   = NSR
         ITS   = 0
         LMOS1 = 1
         JMOS1 = 0
         DO 101 NSS = 1,NSR
          IF(NSS.NE.1) THEN
            ITS   = ITS+MASH(NSS-1)
            LMOS1 = LMOS1+MBAS(NSS-1)*MORB(NSS-1)
            JMOS1 = JMOS1+MORB(NSS-1)
          END IF
          NSPQRS= IBTXOR(NSPQR-1,NSS-1) + 1
         IF(NSPQRS .NE. INTSYM) GO TO 101
          NBS   = MBAS(NSS)
          LMOS  = LMOS1+NBS*MFRO(NSS)
          LMOS2 = LMOS+NBS*MISH(NSS)
          JMOS2 = JMOS1+MFRO(NSS)+MISH(NSS)
          NOS   = MORB(NSS)
          NOCS  = MASH(NSS)
          KEEPS = KEEP(10+NSS)
          ISS   = NSS
C
C
          IF (NOP*NOQ*NOR*NOS.EQ.0) GO TO 101
          NOCTOT = MIN(1,NOCP) + MIN(1,NOCQ)
     &           + MIN(1,NOCR) + MIN(1,NOCS)
          IF (NOCTOT .LT. 2) GO TO 101
C         ... at least two occupied symmetries needed
          IF (MCOPT .AND. NOCTOT .LT. 4) GO TO 101
C         ... for MCSCF optimization only occupied symmetries needed
C
          KEEPT = KEEPP+KEEPQ+KEEPR+KEEPS
          IF (KEEPT.NE.0) THEN
C           NO AO INTEGRALS AVAILABLE BUT MO INTEGRALS ARE NEEDED
C           PROBABLY AN ERROR IN INPUT FOR INTSORT
            WRITE (LUPRI,'(/A)')
     &         ' TR2CTL error, integrals needed but KEEPT .ne. 0'
            WRITE(LUPRI,2000)
            WRITE(LUPRI,2100) ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,
     &                         NOP,NOQ,NOR,NOS
            GO TO 901
          END IF
C
          NCHAIN = ( (ISP**2 - ISP)/2 + ISQ - 1 ) * NSYMP
     &           +   (ISR**2 - ISR)/2 + ISS
          LADX = IDATA(NCHAIN,1)
          NPQ  = IDATA(NCHAIN,2)
          IF (LADX .LT. 0 .OR. NPQ .LT. 0) THEN
            WRITE (LUPRI,'(/A,I7,A,I7/A,I7)')
     &         ' TR2CTL error, LADX =',LADX,' for NCHAIN =',NCHAIN,
     &         '               NPQ  =',NPQ
            WRITE(LUPRI,2000)
            WRITE(LUPRI,2100) ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,
     &                         NOP,NOQ,NOR,NOS
            CALL QUIT(' TR2CTL error, LADX .le. 0')
          END IF
C
C         Calling sequence for second order transformation TRA2
C         for symmetry quadruplet (isp,isq;isr,iss)
C         First allocate and check memory
C
          NBPQ = NBP*NBQ
          IF(ISP.EQ.ISQ) NBPQ = (NBP**2+NBP)/2
          NBRS = NBR*NBS
          IF(ISR.EQ.ISS) NBRS = (NBR**2+NBR)/2
          NOTU = NOCR*NOCS
          IF(ISR.EQ.ISS) NOTU = (NOCR**2+NOCR)/2
C
          LW1    = 1
          LW2    = LW1+MAX(NPQ*NBRS,NBP*NOQ,NBQ*NOCP)
          LW3    = LW2
     *           + MAX(NBR*NBS,NBP*NBQ,NOP*NOR,NOQ*NOR,NOP*NOS,NOQ*NOS)
          LW4    = LW3+MAX(NBR*NOCS,NBS*NOCR)
C         LX2X3 is buffer space needed by TRRDAO,
C         LX2X3 is defined in TRASET
          LW4    = MAX(LW4, LW2 + LX2X3)
          MEMLFT = MEMX-LW4+1
C         Core allocation for the sorting areas RUPQ, USPQ, and TUPQ
          IF (NOCOUL .OR. NOTU.EQ.0) THEN
           NBPQMAX= NBPQ
           LTUPQCM= 0
           LTUPQ  = 0
          ELSE
           NBPQMAX= MIN( MEMLFT/NOTU, NBPQ)
           LTUPQCM= NOTU
           LTUPQ  = NOTU*NBPQMAX
          END IF
          IF (NOEXCH) THEN
           LTUPQM = LTUPQCM
           LRUPQ  = 0
           LUSPQ  = 0
          ELSE
           L2M    = MAX(NOCQ*NOCR*NOP,NOCP*NOCR*NOQ)
C          ... L2M is minimum size of LRUPQ+LTUPQ
           LTUPQXM= MAX(NOCQ*NOCS*NOP,NOCP*NOCS*NOQ)
           LTUPQM = MAX(LTUPQCM,LTUPQXM)
           LRUPQM = NBR*NOCS
           LUSPQM = NBS*NOCR
           IF (LRUPQM+LUSPQM .NE. 0)
     &     NBPQMAX  = MIN(NBPQMAX, (MEMLFT-LTUPQM)/(LRUPQM+LUSPQM))
C
C          --------
C          Find optimal NBPQX, number of simultaneous P,Q
  100      IF (NBPQMAX .LE. 0) GO TO 902
           MEMLFTX= MAX(NBPQMAX*LTUPQCM,LTUPQXM,
     &                  MEMLFT - NBPQMAX*(LRUPQM+LUSPQM))
C          Memory left for LTUPQ with min. NBPQMAX elements in each record
C          on any temp. files for RUPQ and USPQ (if possible).
C          We try if we can avoid temp. file for TUPQ.
C          Now find LTUPQ with these constraints:
           LTUPQ  = NBPQMAX*LTUPQCM
           NBRX   = NBR
           IF (LTUPQXM .NE. 0) THEN
              NBRX   = MIN( MEMLFTX/LTUPQXM, NBR)
              LTUPQ  = MAX( LTUPQM, NBRX*LTUPQXM )
           END IF
C          Memory left for RUPQ and USPQ:
           MEMLFTX= MEMLFT - LTUPQ
           NBPQX  = NBPQMAX
           IF (LRUPQM+LUSPQM .NE. 0)
     &     NBPQX  = MIN(NBPQX, MEMLFTX/(LRUPQM+LUSPQM))
           IF (LUSPQM .NE. 0)
     &     NBPQX  = MIN(NBPQX, (MEMLFT-L2M)/LUSPQM)
           LRUPQ  = NBPQX*LRUPQM
           LUSPQ  = NBPQX*LUSPQM
           LTUPQ  = MAX(NBPQX*LTUPQCM, NBRX*LTUPQXM)
           MEMT   = LRUPQ+LUSPQ+LTUPQ
C
           NBPQMAX = NBPQMAX  - 32
           IF (NBPQX .LT. NBPQMAX) GO TO 100
C          -------- end Find optimal NBPQX
          END IF
C
          LTUPQ  = MEMLFT-LRUPQ-LUSPQ
C
          IF (.NOT. NOEXCH) THEN
           IF(LRUPQ.LT.NBR*NOCS)  GO TO 902
           IF(LUSPQ.LT.NBS*NOCR)  GO TO 902
          END IF
          IF(LTUPQ.LT.LTUPQM)    GO TO 902
          IF(LRUPQ+LTUPQ.LT.L2M) GO TO 902
C
          LW5 = LW4+LUSPQ
          LW6 = LW5+LRUPQ
          CALL TRA2(CMO,NOTUVX,TUVX,NOCOUL,NOEXCH,ICDTRA,
     &              WORK(LW1),WORK(LW2),WORK(LW3),
     &              WORK(LW4),WORK(LW5),WORK(LW6))
C         CALL TRA2(CMO,NOTUVX,TUVX,NOCOUL,NOEXCH,ICDTRA,
C    &              X1,X2,X3, USPQ,RUPQ,TUPQ)
          IF (IPRTRA .GE. 2) THEN
             CALL TIMING(CPT,CPE,TIOT,TIOE)
             IF (IPRTRA .GE. 5) WRITE(LUPRI,2000)
             WRITE(LUPRI,2100) ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,
     &                          NOP,NOQ,NOR,NOS,MEMLFT,LTUPQ,CPE,TIOE
 2100        FORMAT(1X,4I2,1X,4I5,1X,4I5,2I12,2F9.3)
          END IF
  101    CONTINUE
  102   CONTINUE
  103  CONTINUE
  104 CONTINUE
C
      IF (IPRTRA .GT. 0) THEN
         CALL TIMING(CPT,CPE,TIOT,TIOE)
         WRITE(LUPRI,2200) CPT,TIOT
      END IF
#if defined (SYS_UNIX) || defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_HAL) || defined (SYS_T90)
 2200 FORMAT(/' TOTAL USER TIME(SEC)',F9.3,
     &        ' TOTAL SYSTEM TIME(SEC)',F9.3)
#else
 2200 FORMAT(/' TOTAL CPU TIME(SEC)',F9.3,' TOTAL I/O TIME(SEC)',F9.3)
#endif
C
C     FINALLY save first free address on LUMINT in IADOUT(1,LIADUT) and
C     WRITE OUT THE DAFILE ADDRESS LIST ON UNIT 13
C
      IADOUT(1,LIADUT) = IDISKM
      IDISKM=0
      CALL DAWRITE(LUMINT,IADOUT,3*LIADUT,IDISKM)
C
      IF (IPRTRA .GT. 2) THEN
         WRITE (LUPRI,'()')
         WRITE (LUPRI,*) 'Positive addresses in IADOUT for LUMINT:'
         WRITE (LUPRI,*) ' index(sympq) index(symrs) IADC  IADE1 IADE2'
         DO 200 J = 1,36
         DO 200 I = 1,36
            NIJ = I + (J-1)*NSYMP
            IF (IADOUT(1,NIJ) .GT. 0 .OR.
     &          IADOUT(2,NIJ) .GT. 0 .OR.
     &          IADOUT(3,NIJ) .GT. 0 ) THEN
               WRITE (LUPRI,*) I,J,(IADOUT(K,NIJ),K=1,3)
            END IF
  200    CONTINUE
         WRITE (LUPRI,*) 'Address for ICDTRA',IADOUT(3,LIADUT)
         WRITE (LUPRI,*) 'Address for CMO   ',IADOUT(2,LIADUT)
         WRITE (LUPRI,*) 'First free address',IADOUT(1,LIADUT)
         WRITE (LUPRI,'()')
      END IF
C
C
      CALL QEXIT ('TR2CTL')
C
      RETURN
C
C     HERE IF interface FROM SORT of AO integrals IN ERROR
C
  901 WRITE(LUPRI,9010) (KEEP(10+I),I=1,MSYM)
      WRITE(LUPRI,9011) (MASH(I),I=1,MSYM)
 9010 FORMAT(/5X,'ERROR IN KEEP PARAMETER FROM INTSORT FILE:  ',8I5)
 9011 FORMAT( 5X,'NOT CONSISTENT WITH OCCUPIED ORBITAL SPACE: ',8I5,
     &       /5X,'PROGRAM STOP IN SUBROUTINE TR2CTL')
      CALL QTRACE(LUPRI)
      CALL QUIT('TR2CTL, error in KEEP parameter from INTSORT file')
C
C     HERE IF NOT ENOUGH CORE SPACE
C
  902 WRITE(LUPRI,9020) MEMLFT,LRUPQ,NBR*NOCS,LUSPQ,NBS*NOCR,
     &                   LTUPQ,LTUPQM,LRUPQ+LTUPQ,L2M
 9020 FORMAT(/' NOT ENOUGH CORE SPACE FOR SORTING IN TRA2'
     &       /' TOTAL SORTING SPACE IS',I12,
     &       /' STEP1: AVAILABLE IS',I12,',  NEEDED IS',I12,
     &       /' STEP2:     "       ',I12,',    "      ',I12,
     &       /' STEP3:     "       ',I12,',    "      ',I12,
     &       /' STEP4:     "       ',I12,',    "      ',I12)
      CALL QTRACE(LUPRI)
      CALL QUIT('TR2CTL, not enough core space for sorting.')
      END
C  /* Deck tra2 */
      SUBROUTINE TRA2(CMO,NOTUVX,TUVX,NOCOUL,NOEXCH,ICDTRA,
     &                X1,X2,X3,USPQ,RUPQ,TUPQ)
C
C     SECOND ORDER TWO-ELECTRON TRANSFORMATION ROUTINE
C
C     THIS ROUTINE IS CALLED FOR EACH SYMMETRY BLOCK OF INTEGRALS
C     (ISP,ISQ,ISR,ISS) WITH ISP.GE.ISQ AND ISR.GE.ISS.
C
C     INTEGRALS (AB/TU) ARE ALWAYS GENERATED if NOCOUL is false
C
#if defined (VAR_BJORN)
C     EXCHANGE INTEGRALS (AT/BU) ARE GENERATED AS FOLLOWS:
C     (AT/BU) IF ISP.GE.ISR
C     (AT/UB) IF ISP.GT.ISS AND ISP.NE.ISQ
C     (TA/BU) IF ISQ.GT.ISR AND ISP.NE.ISQ
C     (TA/UB) IF ISQ.GE.ISS AND ISP.NE.ISQ
--- 891229-hjaaj: modified to handle INTSYM .ne. 1
#else
C     If NOEXCH is false then
C     EXCHANGE INTEGRALS <AB/TU> ARE GENERATED from (PQ/RS) AS FOLLOWS:
C     case 1a: (AT/BU) IF ISP.GE.ISR
C     case 2a: (AT/UB) IF ISP.GE.ISS AND ISR.NE.ISS
C     case 1b: (TA/BU) IF ISQ.GE.ISR AND ISP.NE.ISQ
C     case 2b: (TA/UB) IF ISQ.GE.ISS AND ISP.NE.ISQ AND ISR.NE.ISS
C
C  Notes:
C     sym(A) .ge. sym(B) has been used
C     case 1* are obtained from RUPQ; common condition: ISP.ge.ISR
C     case 2* are obtained from USPQ; common condition: ISP.ge.ISS
C                                                   and ISR.ne.ISS
#endif
C
C     ********** IBM-3090 RELEASE 87 09 14 **********
C
#include "implicit.h"
      DIMENSION CMO(*),X1(*),X2(*),X3(*),RUPQ(*),USPQ(*),TUPQ(*),TUVX(*)
      DIMENSION ICDTRA(*)
      LOGICAL   NOTUVX, NOCOUL, NOEXCH, OPHLF1, OPHLF2, OPHLF3
#include "iratdef.h"
c#include "dacodes.h"
C
C Used from common blocks:
C   TRINP  : ?,ITUSYM,?
C   TRUNIT : ?
C   TRINFO : NSYMP, LX2X3, ?
C   INFTRA : IPRTRA
C   INFTAP : LUORDA,?
C
#include "priunit.h"
#include "trinp.h"
#include "trunit.h"
#include "trinfo.h"
#include "inftra.h"
#include "inftap.h"
C
      LOGICAL   DOTUPQ,DORUPQ,DOUSPQ
C
C
      CALL QENTER('TRA2  ')
C
      LUHLF1 = -1
      LUHLF2 = -1
      LUHLF3 = -1
      DOTUPQ = .NOT. NOCOUL .AND. NOCR*NOCS .GT. 0
      IF (NOEXCH) THEN
         DORUPQ = .FALSE.
         DOUSPQ = .FALSE.
      ELSE
         DORUPQ = (NOCQ*NOCS.GT.0 .AND. ISP.GE.ISR) .OR.
     *            (NOCP*NOCS.GT.0 .AND. ISQ.GE.ISR)
C        DORUPQ = (case 1a integrals) .or. (case 1b integrals)
C        case 1a: (PQ/RS) -> (AT/BU) IF ISP.GE.ISR
C        case 1b: (PQ/RS) -> (TA/BU) IF ISP.GE.ISS AND ISP.NE.ISQ
         IF (ISR .EQ. ISS) THEN
            DOUSPQ = .FALSE.
         ELSE
            DOUSPQ = (NOCQ*NOCR.GT.0 .AND. ISP.GE.ISS) .OR.
     *               (NOCP*NOCR.GT.0 .AND. ISQ.GE.ISS)
C           DOUSPQ = (case 2a integrals) .or. (case 2b integrals)
C           case 2*: ISR .ne. ISS (ISR.eq.ISS done by DORUPQ)
C           case 2a: (PQ/RS) -> (AT/UB) IF ISP.GE.ISS
C           case 2b: (PQ/RS) -> (TA/UB) IF ISQ.GE.ISS AND ISP.NE.ISQ
         END IF
      END IF
      IF (IPRTRA .GE. 3) THEN
         WRITE (LUPRI,'(/A/A,4I3/)')
     &      ' <<<<<< Test output from TRA2 >>>>>>',
     &      ' ISP,ISQ,ISR,ISS =',ISP,ISQ,ISR,ISS
         WRITE (LUPRI,'(3(5X,A,L10))')
     &      'Do TUPQ :',DOTUPQ,
     &      'Do RUPQ :',DORUPQ,
     &      'Do USPQ :',DOUSPQ
      END IF
      IF (.NOT.DOTUPQ .AND. .NOT.DORUPQ .AND. .NOT.DOUSPQ) GO TO 9999
      NORU =NBR*NOCS
      NOUS =NBS*NOCR
      NOTU =NOCR*NOCS
      IF(ISR.EQ.ISS) NOTU=(NOCR**2+NOCR)/2
C
C     *****************************************************************
C     *****************************************************************
C
C     CHECK FOR IN CORE OR OUT OF CORE TRANSFORMATION
C
C     1. SORT OF PARTIALLY TRANSFORMED INTEGRALS (RU/PQ) ON UNIT LUHLF1
C     2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (US/PQ) ON UNIT LUHLF2
C     3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (TU/PQ) ON UNIT LUHLF3
C
      OPHLF1 = .FALSE.
      OPHLF2 = .FALSE.
      OPHLF3 = .FALSE.
      IF (DORUPQ) THEN
       IPQMX1=NBPQ
       IF(NBPQ*NORU.GT.LRUPQ) THEN
        IPQMX1=LRUPQ/NORU
        IF (IPRTRA .GE. 5) WRITE(LUPRI,*)
     *   'OUT OF CORE SORT FOR INTEGRALS (RU/PQ)',IPQMX1
        IAD1S=0
        CALL DAOPEN(LUHLF1,'HALF1.DA')
        CALL DASKIP(LUHLF1,IRAT*IPQMX1,IAD1S)
        OPHLF1 = .TRUE.
       ENDIF
       IAD1=0
      END IF
C
      IF (DOUSPQ) THEN
       IPQMX2=NBPQ
       IF(NBPQ*NOUS.GT.LUSPQ) THEN
        IPQMX2=LUSPQ/NOUS
        IF (IPRTRA .GE. 5) WRITE(LUPRI,*)
     *   'OUT OF CORE SORT FOR INTEGRALS (US/PQ)',IPQMX2
        IAD2S=0
        CALL DAOPEN(LUHLF2,'HALF2.DA')
        CALL DASKIP(LUHLF2,IRAT*IPQMX2,IAD2S)
        OPHLF2 = .TRUE.
       ENDIF
       IAD2=0
      END IF
C
      IF (DOTUPQ) THEN
       IPQMX3=NBPQ
       IF(NBPQ*NOTU.GT.LTUPQ) THEN
        IPQMX3=LTUPQ/NOTU
        IF (IPRTRA .GE. 5) WRITE(LUPRI,*)
     *   'OUT OF CORE SORT FOR INTEGRALS (TU/PQ)',IPQMX3
        IAD3S=0
        CALL DAOPEN(LUHLF3,'HALF3.DA')
        CALL DASKIP(LUHLF3,IRAT*IPQMX3,IAD3S)
        OPHLF3 = .TRUE.
       ENDIF
       IAD3=0
      END IF
C
C     *****************************************************************
C     *****************************************************************
C
C     Begin first partial transformation
C     START LOOP OVER SORTED AO-INTEGRALS: NPQ PQ-PAIRS IN EACH BUFFER
C
      IOUT1=0
      IOUT2=0
      IOUT3=0
      NX1  =NPQ*NBRS
      IDA  =LADX
      IPQ  =0
      LPQ  =NPQ
      IRSST=1-NBRS
      DO 11 NP=1,NBP
       NQM=NBQ
       IF(ISP.EQ.ISQ) NQM=NP
       DO 10 NQ=1,NQM
        IPQ=IPQ+1
        IF(LPQ.EQ.NPQ) THEN
         LPQ=0
         IRSST=1-NBRS
C
C        READ IN A BLOCK OF INTEGRALS FOR NPQ PQ-VALUES
C
Cold     CALL DAFILE(LUORDA,DAREAD,X1,IRAT*NX1,IDA)
         CALL TRRDAO(LUORDA,X1,NX1,X2,LX2X3,IDA)
C        CALL TRRDAO(LUORDA,X1,NX1,WRK ,LWRK ,IDA)
        ENDIF
        LPQ=LPQ+1
        IRSST=IRSST+NBRS
C
C       START TRANFORMATION OF THIS PQ PAIR
C
        IF(ISR.EQ.ISS) THEN
         CALL TRASQ(X1(IRSST),X2,1,NBS,NBS,ITUSYM)
        ELSE
         CALL DCOPY(NBRS,X1(IRSST),1,X2,1)
        ENDIF
        IF (IPRTRA .GT. 50) THEN
         IF(IPRTRA .GT. 100 .OR. NP.EQ.1.AND.NQ.EQ.1)
     &   WRITE(LUPRI,6000) NP,NQ,(X2(I),I=1,NBR*NBS)
 6000    FORMAT(' AO INTEGRALS FOR P,Q =',2I5/(10F10.6))
        END IF
C
C       INTEGRALS (PQ/US)
C
        IF (DOUSPQ) THEN
C        ... exchange integrals case 2*
         CALL MXMA(X2,        1,NBS,
     &             CMO(LMOR2),1,NBR,
     &             X3,        1,NBS,
     &             NBS,NBR,NOCR)
         IF (IPRTRA .GT. 50) THEN
          IF(IPRTRA .GT. 100 .OR. NP.EQ.1.AND.NQ.EQ.1)
     &    WRITE(LUPRI,6001) NP,NQ,(X3(I),I=1,NBS*NOCR)
 6001     FORMAT(' PQUS FOR P,Q =',2I5/(10F10.6))
         END IF
C
C        SORT THESE INTEGRALS AS (US/PQ)
C
         IOUT2 = IOUT2 + 1
         IF(IOUT2.GT.IPQMX2) THEN
          DO 2 I=1,NOUS
           CALL DAWRITE(LUHLF2,USPQ(1+IPQMX2*(I-1)),IRAT*IPQMX2,
     &                 IAD2)
    2     CONTINUE
          IOUT2=1
         ENDIF
         CALL DCOPY(NOUS,X3,1,USPQ(IOUT2),IPQMX2)
        ENDIF
C
C       INTEGRALS (PQ/RU)
C
        IF(NOCS.NE.0) THEN
         CALL MXMA(X2,        NBS,1,
     &             CMO(LMOS2),1,NBS,
     &             X3,        1,NBR,
     &             NBR,NBS,NOCS)
         IF (IPRTRA .GT. 50) THEN
          IF(IPRTRA .GT. 100 .OR. NP.EQ.1.AND.NQ.EQ.1)
     &    WRITE(LUPRI,6002) NP,NQ,(X3(I),I=1,NBR*NOCS)
 6002     FORMAT(' PQRU FOR P,Q =',2I5/(10F10.6))
         END IF
        END IF
C
C       SORT THESE INTEGRALS AS (RU/PQ)
C
        IF (DORUPQ) THEN
C        ... exchange integrals case 1*
          IOUT1=IOUT1+1
          IF(IOUT1.GT.IPQMX1) THEN
           DO 4 I=1,NORU
            CALL DAWRITE(LUHLF1,RUPQ(1+IPQMX1*(I-1)),IRAT*IPQMX1,
     &                   IAD1)
    4      CONTINUE
           IOUT1=1
          ENDIF
          CALL DCOPY(NORU,X3,1,RUPQ(IOUT1),IPQMX1)
        END IF
C
C       INTEGRALS (PQ/TU)
C
        IF(DOTUPQ) THEN
         IF(ISR.EQ.ISS) THEN
          CALL MXMT(X3,        NBR,1,
     &              CMO(LMOR2),1,NBR,
     &              X2,
     &              NOCR,NBR)
         ELSE
          CALL MXMA(X3,        NBR,1,
     &              CMO(LMOR2),1,NBR,
     &              X2,        1,NOCS,
     &              NOCS,NBR,NOCR)
         ENDIF
         IF (IPRTRA .GT. 50) THEN
          IF(IPRTRA .GT. 100 .OR. NP.EQ.1.AND.NQ.EQ.1)
     &    WRITE(LUPRI,6100) NP,NQ,(X2(I),I=1,NOCR*NOCS)
 6100     FORMAT(' INTEGRALS PQTU FOR P,Q =',2I5/(10F10.6))
         END IF
C
C        SORT INTEGRALS (PQ/TU) INTO TUPQ (SORT AFTER PQ INSTEAD OF TU)
C
         IOUT3 = IOUT3 + 1
         IF(IOUT3.GT.IPQMX3) THEN
          DO 6 I=1,NOTU
           CALL DAWRITE(LUHLF3,TUPQ(1+IPQMX3*(I-1)),IRAT*IPQMX3,
     &                 IAD3)
    6     CONTINUE
          IOUT3=1
         ENDIF
         CALL DCOPY(NOTU,X2,1,TUPQ(IOUT3),IPQMX3)
        ENDIF
   10  CONTINUE
C      ... end do NQ
   11 CONTINUE
C     ... end do NP
C
C       WE NOW HAVE THREE SETS OF PARTIALLY TRANSFORMED INTEGRALS
C       IN TUPQ: (TU/PQ)   TRIANGULAR FOR ISR.EQ.ISS
C       IN RUPQ: (RU/PQ)   IF ISP.GE.ISR
C       IN USPQ: (US/PQ)   IF ISR.GT.ISS .AND. ISP.GE.ISS
C
      IF (IPRTRA .GT. 50) THEN
       IF (DORUPQ) THEN
        IF (IOUT1 .EQ. NBPQ) THEN
         WRITE(LUPRI,6003) (RUPQ(I),I=1,NBPQ*NBR*NOCS)
        ELSE
         WRITE(LUPRI,'(A)') ' RUPQ sorted to LUHLF1'
        END IF
       END IF
       IF (DOUSPQ) THEN
        IF (IOUT2 .EQ. NBPQ) THEN
         WRITE(LUPRI,6013) (USPQ(I),I=1,NBPQ*NBS*NOCR)
        ELSE
         WRITE(LUPRI,'(A)') ' USPQ sorted to LUHLF2'
        END IF
       END IF
       IF (DOTUPQ) THEN
        IF (IOUT3 .EQ. NBPQ) THEN
         WRITE(LUPRI,6200) (TUPQ(I),I=1,NOTU*NBPQ)
        ELSE
         WRITE(LUPRI,'(A)') ' TUPQ sorted to LUHLF3'
        END IF
       END IF
 6003  FORMAT(' RUPQ SORTED'/(10F10.6))
 6013  FORMAT(' USPQ SORTED'/(10F10.6))
 6200  FORMAT(' TUPQ SORTED'/(10F10.6))
      END IF
C
C     EMPTY LAST BUFFERS
C
      IF(DORUPQ .AND. IPQMX1.LT.NBPQ) THEN
       DO 12 I=1,NORU
        CALL DAWRITE(LUHLF1,RUPQ(1+IPQMX1*(I-1)),IRAT*IPQMX1,IAD1)
   12  CONTINUE
      ENDIF
      IF(DOUSPQ .AND. IPQMX2.LT.NBPQ) THEN
       DO 13 I=1,NOUS
        CALL DAWRITE(LUHLF2,USPQ(1+IPQMX2*(I-1)),IRAT*IPQMX2,IAD2)
   13  CONTINUE
      ENDIF
      IF(DOTUPQ .AND. IPQMX3.LT.NBPQ) THEN
       DO 14 I=1,NOTU
        CALL DAWRITE(LUHLF3,TUPQ(1+IPQMX3*(I-1)),IRAT*IPQMX3,IAD3)
   14  CONTINUE
      ENDIF
C
C     *****************************************************************
C     *****************************************************************
C
C     FIRST PARTIAL TRANSFORMATION FINISHED
C     SORTED INTEGRALS ARE ON UNITS LUHLF1 (RUPQ), LUHLF2 (USPQ),
C     AND LUHLF3 (TUPQ), CONTROLLED BY THE ADRESSES IAD1,IAD2, AND IAD3,
C     OR IN CORE (RUPQ, USPQ, AND TUPQ)
C
C     *****************************************************************
C     *****************************************************************
C
C     SECOND HALF TRANSFORMATION FOR INTEGRALS (PQ/TU)
C     FIRST SAVE THE START ADDRESS ON LUMINT FOR THIS BLOCK OF INTEGRALS
C     NOTE THAT THE SYMMETRY LABEL ISPQRS ASSUMES THAT SYMMETRY LOOPS
C     IN THE ORDER T,U,A,B FOR ALL INTEGRAL TYPES.
C
      IF(.NOT.DOTUPQ) GO TO 21
      ISPQRS=((ISR**2-ISR)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISQ
      IADOUT(1,ISPQRS)=IDISKM
      IF (IPRTRA .GT. 50) THEN
       WRITE(LUPRI,7001) ISP,ISQ,ISR,ISS,ISPQRS
 7001  FORMAT(' ADDRESS FOR ABTU',4I2,I11)
      END IF
      ITU=0
      DO 20 NT=1,NOCR
       JT = JMOR2 + NT
       JTR=(JT**2-JT)/2
       NUM=NT
       IF(ISS.NE.ISR) NUM=NOCS
      DO 20 NU=1,NUM
       IPQST=1+NBPQ*ITU
       ITU=ITU+1
C
C      Check if this distribution has been requested
C
       JU = JMOS2 + NU
      IF (ICDTRA(JTR+JU) .EQ. 0) GO TO 20
C
C      READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
       IF(IPQMX3.LT.NBPQ) THEN
        IAD14=IAD3S*(ITU-1)
        IPQ=1
   16   CALL DAREAD(LUHLF3,TUPQ(IPQ),IRAT*IPQMX3,IAD14)
        IPQ=IPQ+IPQMX3
        IAD14=IAD14+IAD3S*(NOTU-1)
        IF(IPQ.LE.NBPQ) GO TO 16
        IPQST=1
       ENDIF
C
       IF(ISP.EQ.ISQ) THEN
        CALL TRASQ(TUPQ(IPQST),X2,1,NBQ,NBQ,+1)
        CALL MXMA(X2,       1,NBQ,
     &            CMO(LMOP),1,NBP,
     &            X1,       1,NBP,
     &            NBQ,NBP,NOP)
        CALL MXMT(X1,       NBQ,1,
     &            CMO(LMOQ),1,NBQ,
     &            X2,
     &            NOP,NBQ)
        IX2=(NOP+NOP**2)/2
       ELSE
        CALL MXMA(TUPQ(IPQST),1,NBQ,
     &            CMO(LMOP),  1,NBP,
     &            X1,         1,NBQ,
     &            NBQ,NBP,NOP)
        CALL MXMA(X1,       NBQ,1,
     &            CMO(LMOQ),1,NBQ,
     &            X2,       NOQ,1,
     &            NOP,NBQ,NOQ)
        IX2=NOP*NOQ
       ENDIF
       IF (IPRTRA .GT. 50) THEN
        WRITE(LUPRI,6010) NT,NU,(X2(I),I=1,IX2)
 6010   FORMAT(' ABTU INTEGRALS FOR NT,NU =',2I5/(10F10.6))
       END IF
C
C      WRITE INTEGRALS (AB/TU) ON OUTPUT UNIT LUMINT
C      INTEGRALS FOR SYMMETRY BLOCK (ISP,ISQ,ISR,ISS) ARE STORED
C      ONE BLOCK FOR EACH TU STARTING AT ADDRESS IADOUT(1,ISPQRS).
C      TRIANGULAR IN AB AND TU IF ISP.EQ.ISQ ( AND ISR.EQ.ISS)
C
       CALL DAWRITE(LUMINT,X2,IRAT*IX2,IDISKM)
C
       IF (.NOT. NOTUVX) THEN
C
C         EXTRACT INTEGRALS WITH ALL INDICES ACTIVE INTO TUVX
C
          IF (ISP.GE.ISR .AND. MASH(ISP)*MASH(ISQ).NE.0)
     &       CALL GTUVX(X2,TUVX,NT,NU,ITP,ITQ,ITR,ITS,ISP,ISQ)
      END IF
C
   20 CONTINUE
   21 CONTINUE
C
C     *****************************************************************
C     *****************************************************************
C
C     Finished Coulomb (Mulliken) integrals
C     now the exchange (Dirac) integrals (unless NOEXCH is set):
C
      IF (DORUPQ .OR. DOUSPQ) THEN
C
C
         CALL TRA2E1(CMO,X1,X2,USPQ,RUPQ,TUPQ,ICDTRA,OPHLF3,
     &                  LUHLF1,LUHLF2,LUHLF3)
C
C
      END IF
C
C     *****************************************************************
C     *****************************************************************
C
C     END OF TRANSFORMATION FOR THIS SYMMETRY BLOCK
C
C     IADOUT CONTAINS START ADRESS FOR EACH TYPE OF INTEGRALS:
C     IADOUT(1,ISPQRS)  COULOMB  INTEGRALS (AB/TU)
C     IADOUT(2,ISPQRS)  EXCHANGE INTEGRALS <AB/TU> FOR SYM T ge SYM U
C     IADOUT(3,ISPQRS)  EXCHANGE INTEGRALS <AB/TU> FOR SYM T lt SYM U
C     THE LAST ADRESS IS ZERO IF SYM T = SYM U
C     TO SEE HOW THE INTEGRALS ARE USED LOOK IN RDINT2
C
C     Close and delete opened temporary files
C
      IF (OPHLF1) CALL DARMOV(LUHLF1)
      IF (OPHLF2) CALL DARMOV(LUHLF2)
      IF (OPHLF3) CALL DARMOV(LUHLF3)
C
 9999 CALL QEXIT ('TRA2  ')
C
      RETURN
      END
C  /* Deck tra2e1 */
      SUBROUTINE TRA2E1(CMO,X1,X2,USPQ,RUPQ,TUPQ,ICDTRA,OPHLF3,
     &                  LUHLF1,LUHLF2,LUHLF3)
C
C     SECOND ORDER TWO-ELECTRON TRANSFORMATION module for exchange
C     integrals.  Called from TRA2.
C
C     THIS ROUTINE IS CALLED FOR EACH SYMMETRY BLOCK OF INTEGRALS
C     (ISP,ISQ,ISR,ISS) WITH ISP.GE.ISQ AND ISR.GE.ISS.
C
C     EXCHANGE INTEGRALS <AB/TU> ARE GENERATED from (PQ/RS) AS FOLLOWS:
C     case 1a: (AT/BU) IF ISP.GE.ISR
C     case 2a: (AT/UB) IF ISP.GE.ISS AND ISR.NE.ISS
C     case 1b: (TA/BU) IF ISQ.GE.ISR AND ISP.NE.ISQ
C     case 2b: (TA/UB) IF ISQ.GE.ISS AND ISP.NE.ISQ AND ISR.NE.ISS
C
C  Notes:
C     sym(A) .ge. sym(B) has been used
C     case 1* are obtained from RUPQ; common condition: ISP.ge.ISR
C     case 2* are obtained from USPQ; common condition: ISP.ge.ISS
C                                                   and ISR.ne.ISS
C
C     ********** IBM-3090 RELEASE 87 09 14 **********
C     Last revision Dec 29, 1989 hjaaj
C
#include "implicit.h"
      DIMENSION CMO(*),X1(*),X2(*),RUPQ(*),USPQ(*),TUPQ(*),ICDTRA(*)
      LOGICAL   OPHLF3
#include "iratdef.h"
c#include "dacodes.h"
C
C Used from common blocks:
C   TRUNIT : ?
C   TRINFO : NSYMP, ?
C   INFTRA : IPRTRA
C   INFTAP : LUMINT
C
#include "priunit.h"
#include "trinp.h"
#include "trunit.h"
#include "trinfo.h"
#include "inftra.h"
#include "inftap.h"
C
      CALL QENTER('TRA2E1')
C
C
      IF (IPRTRA .GT. 50) THEN
         WRITE (LUPRI,'(/A)') ' Test output from TRA2E1'
      END IF
C
C
      NORU =NBR*NOCS
      NOUS =NBS*NOCR
C
C     *****************************************************************
C     *****************************************************************
C
C     Now the exchange (Dirac) integrals:
C
C     Exchange integrals case 1a : (pq/rs) -> (at/bu)
C     SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/RU)-> (AT/BU)
C     IF ISP.EQ.ISR THEN T.GE.U BUT ALWAYS ALL A AND B
C
      NOTU =NOCQ*NOCS
      IF(ISQ.EQ.ISS) NOTU=(NOCQ**2+NOCQ)/2
      IF(ISP.GE.ISR.AND.NOTU.NE.0) THEN
       LAR=LTUPQ/NOTU
       LR=LAR/NOP
       IF(LR.GT.NBR) LR=NBR
       LAR=NOP*LR
       IF(LR.LT.NBR) THEN
        IAD3S=0
        IF (.NOT.OPHLF3) CALL DAOPEN(LUHLF3,'HALF3.DA')
        CALL DASKIP(LUHLF3,IRAT*LAR,IAD3S)
        OPHLF3 = .TRUE.
        IAD3=0
       END IF
       IR=0
       DO 30 NR=1,NBR
        IR=IR+1
       DO 30 NU=1,NOCS
        IRU=NBR*(NU-1)+NR
        IPQST=1+NBPQ*(IRU-1)
C
C       READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
        IF(IPQMX1.LT.NBPQ) THEN
         IAD1=IAD1S*(IRU-1)
         IPQ=1
   22    CALL DAREAD(LUHLF1,RUPQ(IPQ),IRAT*IPQMX1,IAD1)
         IPQ=IPQ+IPQMX1
         IAD1=IAD1+IAD1S*(NORU-1)
         IF(IPQ.LE.NBPQ) GO TO 22
         IPQST=1
        ENDIF
        IF(ISP.EQ.ISQ) THEN
         CALL TRASQ(RUPQ(IPQST),X2,1,NBQ,NBQ,+1)
        ELSE
         CALL DCOPY(NBPQ,RUPQ(IPQST),1,X2,1)
        ENDIF
        IF(ISQ.EQ.ISS) THEN
         IF (IPRTRA .GT. 50) THEN
          WRITE(LUPRI,6004) NR,NU,(X2(I),I=1,NBP*NBQ)
 6004     FORMAT(' X2=RUPQ FOR R AND U =',2I5/(10F10.6))
         END IF
         CALL MXMA(X2,                   NBQ,1,
     &             CMO(LMOQ2+NBQ*(NU-1)),1,NBQ,
     &             X1,                   1,NBP,
     &             NBP,NBQ,NOCQ-NU+1)
         CALL MXMA(X1,       NBP,1,
     &             CMO(LMOP),1,NBP,
     &             X2,       1,NOCQ-NU+1,
     &             NOCQ-NU+1,NBP,NOP)
         IF (IPRTRA .GT. 50) THEN
          WRITE(LUPRI,6005) NR,NU,(X2(I),I=1,NOP*(NOCQ-NU+1))
 6005     FORMAT(' RUAT FOR R AND U =',2I5/(10F10.6))
         END IF
        ELSE
         CALL MXMA(X2,        NBQ,1,
     &             CMO(LMOQ2),1,NBQ,
     &             X1,        1,NBP,
     &             NBP,NBQ,NOCQ)
         CALL MXMA(X1,       NBP,1,
     &             CMO(LMOP),1,NBP,
     &             X2,       1,NOCQ,
     &             NOCQ,NBP,NOP)
         ENDIF
C
C       INTEGRALS (AT/RU) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
C       VALUES OF A,T(T.GE.U IF ISQ.EQ.ISS)
C
C       SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
C       OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
C       INTEGRALS (TU/PQ))
C
        IF(IR.GT.LR) THEN
         IR=1
         DO 24 I=1,NOTU
          CALL DAWRITE(LUHLF3,TUPQ(1+LAR*(I-1)),IRAT*LAR,IAD3)
   24    CONTINUE
        ENDIF
        NAT=0
        DO 26 NA=1,NOP
         NTM=1
         IF(ISQ.EQ.ISS) NTM=NU
        DO 26 NT=NTM,NOCQ
         ITU=NOCS*(NT-1)+NU-1
         IF(ISQ.LT.ISS) ITU=NOCQ*(NU-1)+NT-1
         IF(ISQ.EQ.ISS) ITU=(NT**2-NT)/2+NU-1
         NAT=NAT+1
         TUPQ(LAR*ITU+NOP*(IR-1)+NA)=X2(NAT)
   26   CONTINUE
   30  CONTINUE
C
C      EMPTY LAST BUFFER IF LR.LT.NBR
C
       IF(LR.LT.NBR) THEN
        DO 31 I=1,NOTU
         CALL DAWRITE(LUHLF3,TUPQ(1+LAR*(I-1)),IRAT*LAR,IAD3)
   31   CONTINUE
       ENDIF
C
C      NOW TRANSFORM INDEX R TO B
C
C      FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK
C
       IF(ISQ.GE.ISS) THEN
        ISPQRS=((ISQ**2-ISQ)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISR
        IADOUT(2,ISPQRS)=IDISKM
        NTMAX=NOCQ
        NUMAX=NOCS
        JMOT2=JMOQ2
        JMOU2=JMOS2
       ELSE
        ISPQRS=((ISS**2-ISS)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISR
        IADOUT(3,ISPQRS)=IDISKM
        NTMAX=NOCS
        NUMAX=NOCQ
        JMOT2=JMOS2
        JMOU2=JMOQ2
        IF (IPRTRA .GT. 50) THEN
         WRITE(LUPRI,7002) ISS,ISQ,ISP,ISR,ISPQRS
 7002    FORMAT(' ADDRESS FOR ATBU',4I2,I11)
        END IF
       ENDIF
       NTU=0
       IST=1-NOP*NBR
       DO 40 NT=1,NTMAX
        JT = JMOT2 + NT
        JTR=(JT**2-JT)/2
        NUM=NUMAX
        IF(ISQ.EQ.ISS) NUM=NT
       DO 39 NU=1,NUM
        NTU=NTU+1
        IST=IST+NOP*NBR
C
C      Check if this distribution has been requested
C
        JU = JMOU2 + NU
       IF (ICDTRA(JTR+JU) .EQ. 0) GO TO 39
C
        IF(LR.LT.NBR)THEN
         IAD3=(NTU-1)*IAD3S
         IST=1
   32    CONTINUE
         CALL DAREAD(LUHLF3,TUPQ(IST),IRAT*LAR,IAD3)
         IST=IST+LR*NOP
         IAD3=IAD3+(NOTU-1)*IAD3S
         IF(IST.LE.NBR*NOP) GO TO 32
         IST=1
        ENDIF
        CALL MXMA(TUPQ(IST),1,NOP,
     &            CMO(LMOR),1,NBR,
     &            X2,       NOR,1,
     &            NOP,NBR,NOR)
C
C       Write this block of exchange integrals on LUMINT
C
        IF (IPRTRA .GT. 50) THEN
         WRITE(LUPRI,6008) (X2(I),I=1,NOP*NOR)
 6008    FORMAT(' ATBU INTEGRALS'/(10F10.6))
        END IF
        CALL DAWRITE(LUMINT,X2,IRAT*NOP*NOR,IDISKM)
   39  CONTINUE
   40  CONTINUE
      ENDIF
C
C     Exchange integrals case 1b : (pq/rs) -> (ta/bu)
C     SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/RU)-> (TA/BU)
C
      NOTU=NOCP*NOCS
C     ... ISP.ne.ISS because ISP .gt. ISQ .ge. ISR .ge. ISS
      IF((ISP.NE.ISQ.AND.ISQ.GE.ISR).AND.NOTU.NE.0) THEN
       LAR=LTUPQ/NOTU
       LR=LAR/NOQ
       IF(LR.GT.NBR) LR=NBR
       LAR=NOQ*LR
       IF(LR.LT.NBR) THEN
        IAD3S=0
        IF (.NOT.OPHLF3) CALL DAOPEN(LUHLF3,'HALF3.DA')
        CALL DASKIP(LUHLF3,IRAT*LAR,IAD3S)
        OPHLF3 = .TRUE.
        IAD3=0
       ENDIF
       IRU=0
       IR=0
       DO 50 NR=1,NBR
        IR=IR+1
       DO 50 NU=1,NOCS
        IRU=NBR*(NU-1)+NR
        IPQST=1+NBPQ*(IRU-1)
C
C       READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
        IF(IPQMX1.LT.NBPQ) THEN
         IAD1=IAD1S*(IRU-1)
         IPQ=1
   42    CALL DAREAD(LUHLF1,RUPQ(IPQ),IRAT*IPQMX1,IAD1)
         IPQ=IPQ+IPQMX1
         IAD1=IAD1+IAD1S*(NORU-1)
         IF(IPQ.LE.NBPQ) GO TO 42
         IPQST=1
        ENDIF
        CALL MXMA(RUPQ(IPQST),1,NBQ,
     &            CMO(LMOP2), 1,NBP,
     &            X1,         1,NBQ,
     &            NBQ,NBP,NOCP)
        CALL MXMA(X1,       NBQ,1,
     &            CMO(LMOQ),1,NBQ,
     &            X2,       1,NOCP,
     &            NOCP,NBQ,NOQ)
C
C       INTEGRALS (TA/RU) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
C       VALUES OF T,A . NOTE THAT T AND U HERE ALWAYS ARE OF DIFFERENT
C       SYMMETRIES
C
C       SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
C       OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
C       INTEGRALS (TU/PQ))
C
        IF(IR.GT.LR) THEN
         IR=1
         DO 44 I=1,NOTU
          CALL DAWRITE(LUHLF3,TUPQ(1+LAR*(I-1)),IRAT*LAR,IAD3)
   44    CONTINUE
        ENDIF
        NAT=0
        DO 46 NA=1,NOQ
        DO 46 NT=1,NOCP
         ITU=NOCS*(NT-1)+NU-1
         NAT=NAT+1
         TUPQ(LAR*ITU+NOQ*(IR-1)+NA)=X2(NAT)
   46   CONTINUE
   50  CONTINUE
C
C      EMPTY LAST BUFFER IF LR.LT.NBR
C
       IF(LR.LT.NBR) THEN
        DO 51 I=1,NOTU
         CALL DAWRITE(LUHLF3,TUPQ(1+LAR*(I-1)),IRAT*LAR,IAD3)
   51   CONTINUE
       ENDIF
C
C      NOW TRANSFORM INDEX R TO B
C
C      FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK
C
       ISPQRS=((ISP**2-ISP)/2+ISS-1)*NSYMP+(ISQ**2-ISQ)/2+ISR
       IADOUT(2,ISPQRS)=IDISKM
       IF (IPRTRA .GT. 50) THEN
        WRITE(LUPRI,7003) ISP,ISS,ISQ,ISR,ISPQRS
 7003   FORMAT(' ADDRESS FOR TABU',4I2,I11)
       END IF
       NTU=0
       IST=1-NOQ*NBR
       DO 60 NT=1,NOCP
        JT = JMOP2 + NT
        JTR=(JT**2-JT)/2
       DO 60 NU=1,NOCS
        NTU=NTU+1
        IST=IST+NOQ*NBR
C
C      Check if this distribution has been requested
C
        JU = JMOS2 + NU
       IF (ICDTRA(JTR+JU) .EQ. 0) GO TO 60
C
        IF(LR.LT.NBR)THEN
         IAD3=(NTU-1)*IAD3S
         IST=1
   52    CALL DAREAD(LUHLF3,TUPQ(IST),IRAT*LAR,IAD3)
         IST=IST+LR*NOQ
         IAD3=IAD3+(NOTU-1)*IAD3S
         IF(IST.LE.NBR*NOQ) GO TO 52
         IST=1
        ENDIF
        CALL MXMA(TUPQ(IST),1,NOQ,
     &            CMO(LMOR),1,NBR,
     &            X2,       NOR,1,
     &            NOQ,NBR,NOR)
C
C       Write this block of exchange integrals on LUMINT
C
        CALL DAWRITE(LUMINT,X2,IRAT*NOR*NOQ,IDISKM)
   60  CONTINUE
      ENDIF
C
C     Exchange integrals case 2a : (pq/rs) -> (at/ub)
C     SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/US)-> (AT/UB)
C
      NOTU=NOCQ*NOCR
      IF (ISQ.EQ.ISR) NOTU=(NOCQ**2+NOCQ)/2
      IF (INTSYM .NE. 1) THEN
C     MAERKE hvad med ISQ.EQ.ISR; er det ikke muligt for INTSYM.ne.1?JO
         WRITE (LUPRI,*)
     &      ' error in TRA2E1 case 2a: INTSYM.ne.1 not implemented'
         CALL QTRACE(LUPRI)
         CALL QUIT('TRA2E1 error: INTSYM.ne.1 not implemented for '//
     &             'case 2a')
      END IF
      IF((ISR.NE.ISS.AND.ISP.GE.ISS).AND.NOTU.NE.0) THEN
       LAR=(LRUPQ+LTUPQ)/NOTU
       LR=LAR/NOP
       IF(LR.GT.NBS) LR=NBS
       LAR=NOP*LR
       IF(LR.LT.NBS) THEN
        IAD3S=0
        IF (.NOT.OPHLF3) CALL DAOPEN(LUHLF3,'HALF3.DA')
        CALL DASKIP(LUHLF3,IRAT*LAR,IAD3S)
        OPHLF3 = .TRUE.
        IAD3=0
       ENDIF
       IRU=0
       IR=0
       DO 70 NR=1,NBS
        IR=IR+1
       DO 70 NU=1,NOCR
        IRU=NBS*(NU-1)+NR
        IPQST=1+NBPQ*(IRU-1)
C
C       READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
        IF(IPQMX2.LT.NBPQ) THEN
         IAD2=IAD2S*(IRU-1)
         IPQ=1
   62    CALL DAREAD(LUHLF2,USPQ(IPQ),IRAT*IPQMX2,IAD2)
         IPQ=IPQ+IPQMX2
         IAD2=IAD2+IAD2S*(NOUS-1)
         IF(IPQ.LE.NBPQ) GO TO 62
         IPQST=1
        ENDIF
        CALL MXMA(USPQ(IPQST),NBQ,1,
     &            CMO(LMOQ2), 1,NBQ,
     &            X1,         1,NBP,
     &            NBP,NBQ,NOCQ)
        CALL MXMA(X1,       NBP,1,
     &            CMO(LMOP),1,NBP,
     &            X2,       1,NOCQ,
     &            NOCQ,NBP,NOP)
C
C       INTEGRALS (AT/US) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
C       VALUES OF A,T. NOTE THAT T AND U HAVE DIFFERENT SYMMETRIES.
C
C       SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
C       OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
C       INTEGRALS (TU/PQ))
C
        IF(IR.GT.LR) THEN
         IR=1
         DO 64 I=1,NOTU
          CALL DAWRITE(LUHLF3,RUPQ(1+LAR*(I-1)),IRAT*LAR,IAD3)
   64    CONTINUE
        ENDIF
        NAT=0
        DO 66 NA=1,NOP
        DO 66 NT=1,NOCQ
         ITU=NOCR*(NT-1)+NU-1
         IF(ISQ.LT.ISR) ITU=NOCQ*(NU-1)+NT-1
         NAT=NAT+1
         RUPQ(LAR*ITU+NOP*(IR-1)+NA)=X2(NAT)
   66   CONTINUE
   70  CONTINUE
C
C      EMPTY LAST BUFFER IF LR.LT.NBS
C
       IF(LR.LT.NBS) THEN
        DO 71 I=1,NOTU
         CALL DAWRITE(LUHLF3,RUPQ(1+LAR*(I-1)),IRAT*LAR,IAD3)
   71   CONTINUE
       ENDIF
C
C      NOW TRANSFORM INDEX R TO B
C
C      FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK
C
       IF(ISQ.GE.ISR) THEN
        ISPQRS=((ISQ**2-ISQ)/2+ISR-1)*NSYMP+(ISP**2-ISP)/2+ISS
        IADOUT(2,ISPQRS)=IDISKM
        NTMAX=NOCQ
        NUMAX=NOCR
        JMOT2=JMOQ2
        JMOU2=JMOR2
        IF (IPRTRA .GT. 50) WRITE(LUPRI,7004) ISQ,ISR,ISP,ISS,ISPQRS
       ELSE
        ISPQRS=((ISR**2-ISR)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISS
        IADOUT(3,ISPQRS)=IDISKM
        NTMAX=NOCR
        NUMAX=NOCQ
        JMOT2=JMOR2
        JMOU2=JMOQ2
        IF (IPRTRA .GT. 50) WRITE(LUPRI,7004) ISR,ISQ,ISP,ISS,ISPQRS
       ENDIF
 7004  FORMAT(' ADDRESS FOR ATUB',4I2,I11)
       NTU=0
       IST=1-NBS*NOP
       DO 80 NT=1,NTMAX
        JT = JMOT2 + NT
        JTR=(JT**2-JT)/2
       DO 80 NU=1,NUMAX
        NTU=NTU+1
        IST=IST+NBS*NOP
C
C      Check if this distribution has been requested
C
        JU = JMOU2 + NU
       IF (ICDTRA(JTR+JU) .EQ. 0) GO TO 80
C
        IF(LR.LT.NBS)THEN
         IAD3=(NTU-1)*IAD3S
         IST=1
   72    CALL DAWRITE(LUHLF3,RUPQ(IST),IRAT*LAR,IAD3)
         IST=IST+LR*NOP
         IAD3=IAD3+(NOTU-1)*IAD3S
         IF(IST.LE.NBS*NOP) GO TO 72
         IST=1
        ENDIF
        CALL MXMA(RUPQ(IST),1,NOP,
     &            CMO(LMOS),1,NBS,
     &            X2,       NOS,1,
     &            NOP,NBS,NOS)
C
C       Write this block of exchange integrals on LUMINT
C
        CALL DAWRITE(LUMINT,X2,IRAT*NOS*NOP,IDISKM)
   80  CONTINUE
      ENDIF
C
C     Exchange integrals case 2b : (pq/rs) -> (ta/ub)
C     SECOND PARTIAL TRANSFORMATION FOR INTEGRALS (PQ/US)-> (TA/UB)
C
      NOTU=NOCP*NOCR
      IF(ISP.EQ.ISR) NOTU=(NOCP**2+NOCP)/2
      IF((ISP.NE.ISQ.AND.ISR.NE.ISS.AND.ISQ.GE.ISS).AND.NOTU.NE.0) THEN
       LAR=(LRUPQ+LTUPQ)/NOTU
       LR=LAR/NOQ
       IF(LR.GT.NBS) LR=NBS
       LAR=NOQ*LR
       IF(LR.LT.NBS) THEN
        IAD3S=0
        IF (.NOT.OPHLF3) CALL DAOPEN(LUHLF3,'HALF3.DA')
        CALL DASKIP(LUHLF3,IRAT*LAR,IAD3S)
        OPHLF3 = .TRUE.
        IAD3=0
       END IF
       IRU=0
       IR=0
       DO 90 NR=1,NBS
        IR=IR+1
       DO 90 NU=1,NOCR
        IRU=NBS*(NU-1)+NR
        IPQST=1+NBPQ*(IRU-1)
C
C       READ ONE BUFFER OF INTEGRALS BACK INTO CORE
C
        IF(IPQMX2.LT.NBPQ) THEN
         IAD2=IAD2S*(IRU-1)
         IPQ=1
   82    CALL DAREAD(LUHLF2,USPQ(IPQ),IRAT*IPQMX2,IAD2)
         IPQ=IPQ+IPQMX2
         IAD2=IAD2+IAD2S*(NOUS-1)
         IF(IPQ.LE.NBPQ) GO TO 82
         IPQST=1
        ENDIF
        IF(ISP.EQ.ISR) THEN
         CALL MXMA(USPQ(IPQST),          1,NBQ,
     &             CMO(LMOP2+NBP*(NU-1)),1,NBP,
     &             X1,                   1,NBQ,
     &             NBQ,NBP,NOCP-NU+1)
         CALL MXMA(X1,       NBQ,1,
     &             CMO(LMOQ),1,NBQ,
     &             X2,       1,NOCP-NU+1,
     &             NOCP-NU+1,NBQ,NOQ)
        ELSE
         CALL MXMA(USPQ(IPQST),   1,NBQ,
     &             CMO(LMOP2),1,NBP,
     &             X1,            1,NBQ,
     &             NBQ,NBP,NOCP)
         CALL MXMA(X1,       NBQ,1,
     &             CMO(LMOQ),1,NBQ,
     &             X2,       1,NOCP,
     &             NOCP,NBQ,NOQ)
        ENDIF
C
C       INTEGRALS (TA/US) ARE NOW IN X2 FOR ONE VALUE OF R,U AND ALL
C       VALUES OF A,T
C
C       SORT THESE INTEGRALS AFTER PAIR INDEX TU, IN CORE IF POSSIBLE
C       OTHERWISE USE LUHLF3 FOR TEMPORARY STORAGE (SAME POSITION AS FOR
C       INTEGRALS (TU/PQ))
C
        IF(IR.GT.LR) THEN
         IR=1
         DO 84 I=1,NOTU
          CALL DAWRITE(LUHLF3,RUPQ(1+LAR*(I-1)),IRAT*LAR,IAD3)
   84    CONTINUE
        ENDIF
        NAT=0
        DO 86 NA=1,NOQ
         NTM=1
         IF(ISP.EQ.ISR) NTM=NU
        DO 86 NT=NTM,NOCP
         ITU=NOCR*(NT-1)+NU-1
         IF(ISP.LT.ISR) ITU=NOCP*(NU-1)+NT-1
         IF(ISP.EQ.ISR) ITU=(NT**2-NT)/2+NU-1
         NAT=NAT+1
         RUPQ(LAR*ITU+NOQ*(IR-1)+NA)=X2(NAT)
   86   CONTINUE
   90  CONTINUE
C
C      EMPTY LAST BUFFER IF LR.LT.NBS
C
       IF(LR.LT.NBS) THEN
        DO 91 I=1,NOTU
         CALL DAWRITE(LUHLF3,RUPQ(1+LAR*(I-1)),IRAT*LAR,IAD3)
   91   CONTINUE
       ENDIF
C
C      NOW TRANSFORM INDEX R TO B
C
C       FIRST COMPUTE AND SAVE START ADDRESS FOR THIS SYMMETRY BLOCK
C
       IF(ISP.GE.ISR) THEN
        ISPQRS=((ISP**2-ISP)/2+ISR-1)*NSYMP+(ISQ**2-ISQ)/2+ISS
        IADOUT(2,ISPQRS)=IDISKM
        NTMAX=NOCP
        NUMAX=NOCR
        JMOT2=JMOP2
        JMOU2=JMOR2
        IF (IPRTRA .GT. 50) WRITE(LUPRI,7005) ISP,ISR,ISQ,ISS,ISPQRS
       ELSE
        ISPQRS=((ISR**2-ISR)/2+ISP-1)*NSYMP+(ISQ**2-ISQ)/2+ISS
        IADOUT(3,ISPQRS)=IDISKM
        NTMAX=NOCR
        NUMAX=NOCP
        JMOT2=JMOR2
        JMOU2=JMOP2
        IF (IPRTRA .GT. 50) WRITE(LUPRI,7005) ISR,ISP,ISQ,ISS,ISPQRS
       ENDIF
 7005  FORMAT(' ADDRESS FOR TAUB',4I2,I11)
       NTU=0
       IST=1-NOQ*NBS
       DO 100 NT=1,NTMAX
        JT = JMOT2 + NT
        JTR=(JT**2-JT)/2
        NUM=NUMAX
        IF(ISP.EQ.ISR) NUM=NT
       DO 100 NU=1,NUM
        NTU=NTU+1
        IST=IST+NOQ*NBS
C
C      Check if this distribution has been requested
C
        JU = JMOU2 + NU
       IF (ICDTRA(JTR+JU) .EQ. 0) GO TO 100
C
        IF(LR.LT.NBS)THEN
         IAD3=(NTU-1)*IAD3S
         IST=1
   92    CALL DAREAD(LUHLF3,RUPQ(IST),IRAT*LAR,IAD3)
         IST=IST+LR*NOQ
         IAD3=IAD3+(NOTU-1)*IAD3S
         IF(IST.LE.NBS*NOQ) GO TO 92
         IST=1
        ENDIF
        CALL MXMA(RUPQ(IST),1,NOQ,
     &            CMO(LMOS),1,NBS,
     &            X2,       NOS,1,
     &            NOQ,NBS,NOS)
C
C       Write this block of exchange integrals on LUMINT
C
        CALL DAWRITE(LUMINT,X2,IRAT*NOS*NOQ,IDISKM)
  100  CONTINUE
      ENDIF
C
C
C     *****************************************************************
C     *****************************************************************
C
C     END OF TRA2E1 TRANSFORMATION FOR THIS SYMMETRY BLOCK
C
C
      CALL QEXIT ('TRA2E1')
C
      RETURN
      END
C  /* Deck gtuvx */
      SUBROUTINE GTUVX(ABVX,TUVX,NV,NX,ITP,ITQ,ITR,ITS,ISP,ISQ)
C
C     EXTRACT FROM THE INTEGRALS (AB/VX) INTEGRALS (TU/VX) WITH ALL
C     INDICES ACTIVE. THEY ARE STORED CANONICALLY WITHOUT SYMMETRY
C     BLOCKING
C
C     ********** IBM-3090 RELEASE 87 09 15 **********
C
#include "implicit.h"
      DIMENSION ABVX(*),TUVX(*)
#include "trinp.h"
C
      CALL QENTER('GTUVX ')
C
      NVV=NV+ITR
      NXX=NX+ITS
      NVX=(NVV**2-NVV)/2+NXX
      DO 10 NT=1,MASH(ISP)
       NTT=NT+ITP
       NA=NT+MISH(ISP)
       NUM=MASH(ISQ)
       IF(ISP.EQ.ISQ) NUM=NT
      DO 10 NU=1,NUM
       NUU=NU+ITQ
       NTU=(NTT**2-NTT)/2+NUU
#if defined (VAR_PACKH2AC)
       IF(NTU.GT.NVX) GO TO 10
       NTUVX=(NVX**2-NVX)/2+NTU
#else
       NTUVX = (NVX-1)*MMASHX + NTU
#endif
C
       NB=NU+MISH(ISQ)
       NAB=MORB(ISQ)*(NA-1)+NB
       IF(ISP.EQ.ISQ) NAB=(NA**2-NA)/2+NB
       TUVX(NTUVX)=ABVX(NAB)
   10 CONTINUE
C
      CALL QEXIT ('GTUVX ')
C
      RETURN
      END
C  /* Deck trasq */
      SUBROUTINE TRASQ(A, B,ICB,IRB, NROW, IASYM)
C
C 24-Jan-1990 hjaaj, IASYM option
C
C If IASYM .lt. 0 then the matrix is antisymmetric,
C and A is the lower triangle.
C
C DIMENSION B(1:IRB:ICB,*)
C B(1+IR*ICB+IC*IRB) corresponds to B(IR+1,IC+1)
C For normal square "operation", ICB=1, IRB=NROW
C
#include "implicit.h"
      DIMENSION A(*),B(*)
      IF (IASYM .GE. 0) THEN
         IND=0
         DO 10 IROW=0,NROW-1
         DO 10 ICOL=0,IROW
            IND=IND+1
            B(1+IROW*ICB+ICOL*IRB)=A(IND)
            B(1+ICOL*ICB+IROW*IRB)=A(IND)
   10    CONTINUE
      ELSE
         IND=0
         DO 20 IROW=0,NROW-1
         DO 20 ICOL=0,IROW
            IND=IND+1
            B(1+IROW*ICB+ICOL*IRB)= A(IND)
            B(1+ICOL*ICB+IROW*IRB)=-A(IND)
   20    CONTINUE
      END IF
      RETURN
      END
C  /* Deck traini */
      SUBROUTINE TRAINI(NTRLVL,CMO,ICDTRA,ITRTYP)
C
C Dec 1989 Hans Joergen Aa. Jensen
C
#include "implicit.h"
#include "maxorb.h"
      DIMENSION CMO(LCMO),ICDTRA(MMORBX),ITRTYP(MXCORB)
C MAERKE ??      LOGICAL   NOTUVX, NOEXCH, MCOPT
#include "iratdef.h"
c#include "dacodes.h"
#include "dalistgs.h"
C
C Used from common blocks:
C  TRINP  : LCMO, MMORBX, MSYM, MBAS(), INTSYM, ...
C  TRUNIT : LIADUT,IADOUT
C  INFTAP : LUMINT
C
#include "trinp.h"
#include "trunit.h"
#include "inftap.h"
C
      CALL QENTER('TRAINI')
C
C     SET ADDRESS FIELD FOR OUTPUT INTEGRAL FILE
C     Initialize to -2 which is used as error code
C                   -1 will mean no integrals
C
      DO 1 I=1,LIADUT
       IADOUT(1,I)=-2
       IADOUT(2,I)=-2
       IADOUT(3,I)=-2
    1 CONTINUE
      IDISKM = 0
      CALL DASKIP(LUMINT,3*LIADUT,IDISKM)
      NLIST = 10
C     MAERKE skal flere med i GSLIST?
#if 1
      CALL DAWRLIST(LUMINT,IDISKM,
     &              NLIST,NTRLVL,1,INTSYM,1,MSYM,1,LCMO,1,
     &              MFRO,8,MISH,8,MASH,8,MORB,8,MBAS,8,
     &              ITRTYP,MXCORB)
#else
      CALL GSLIST(LISTGS,NLIST,NTRLVL,1,INTSYM,1,MSYM,1,LCMO,1,
     &            MFRO,8,MISH,8,MASH,8,MORB,8,MBAS,8,
     &            ITRTYP,MXCORB)
      CALL DAFILE(LUMINT,DAGWRI,LISTGS,NLIST,IDISKM)
#endif
      IADOUT(3,LIADUT) = IDISKM
      CALL DAWRITE(LUMINT,ICDTRA,MMORBX,IDISKM)
      IADOUT(2,LIADUT) = IDISKM
      CALL DAWRITE(LUMINT,CMO,IRAT*LCMO,IDISKM)
      IADOUT(1,LIADUT) = IDISKM
      IDISKM = 0
      CALL DAWRITE(LUMINT,IADOUT,3*LIADUT,IDISKM)
      CALL QEXIT ('TRAINI')
      RETURN
      END
C  /* Deck traset */
      SUBROUTINE TRASET(JRDAO,LWORK)
C
C Dec 1989 Hans Joergen Aa. Jensen
C Define TRINP and TRUNIT
C
C Input:
C     JRDAO=0 : from AOORDINT, ordered integral file produced by
C               the Lund INTSORT program.
C     JRDAO=1 : derivative integrals from Abacus
C     JRDAO=2 : Spin-orbit integrals
C     JRDAO=3 : from AOTWOINT
C
#include "implicit.h"
#include "iratdef.h"
C
C Used from common blocks:
C   INFORB : NSYM, ..., NBAS(), NNORBX,NCMOT
C   TRINP  : MSYM, ..., MBAS(), MMORBX,LCMO, ITUSYM
C   TRINFO : LX2X3
C   TRBINS : LMSBIN
C   INFTRA : IPRTRA
C
#include "priunit.h"
#include "inforb.h"
#include "trinp.h"
#include "trinfo.h"
#include "inftra.h"
C
C     Local variables
C
      LOGICAL FOPEN
C
      CALL QENTER('TRASET')
C
C     Transfer information from INFORB to TRINP
C
      MSYM  = NSYM
      MORBT = 0
      DO 100 I = 1,8
         MFRO(I) = NFRO(I)
         MISH(I) = NISH(I)
         MASH(I) = NASH(I)
         MORB(I) = NORB(I)
         MBAS(I) = NBAS(I)
         MORBT   = MORBT + MORB(I)
  100 CONTINUE
      MMASHX = NNASHX
      MMORBX = NNORBX
      LCMO   = NCMOT
C
C
C     Set code for reading AO integrals
C     IRDAO=0 : from AOORDINT, ordered integral file produced by
C               the Lund INTSORT program.
C     IRDAO=1 : derivative AO integrals from Abacus
C     IRDAO=2 : Spin-orbit AO integrals.
C     IRDAO=3 : from AOTWOINT
C
C     (AB/TU) = ITUSYM*(AB/UT) = (BA/TU) = ITUSYM*(BA/UT)
C     ITUSYM = -1 for spin-orbit integrals (JRDAO .EQ. 2)
C
C     Transfer unit numbers to TRUNIT
C     CALL TRASE[02] to obtain KEEP and IDATA in TRUNIT
C
      IRDAO  = JRDAO
      IF (IRDAO .EQ. 0) THEN
         ITUSYM =  1
         LX2X3  =  0
         CALL TRASET_0
! hjaaj Sep 2010: derivative integrals have not been
!     implemented here for many years: since abadtra.F was
!     written.
!     ELSE IF (IRDAO.EQ.1) THEN
!        ITUSYM =  1
!        LX2X3  = ( LMSBIN*(IRAT+1)+2 - 1 ) / IRAT + 1
!
      ELSE IF (IRDAO.EQ.2) THEN
         ITUSYM = -1
         CALL TRASET_2(LX2X3, LWORK)
      ELSE IF (IRDAO .EQ. 3) THEN
         ITUSYM = 1
         CALL TRASET_3(LX2X3, LWORK)
      ELSE
         WRITE (LUPRI,*) 'ERROR TRASET: Illegal IRDAO =',IRDAO
         CALL QTRACE(LUPRI)
         CALL QUIT('ERROR TRASET: Illegal IRDAO')
      END IF
C
      CALL QEXIT ('TRASET')
      RETURN
      END
C  /* Deck trase0 */
      SUBROUTINE TRASET_0
C
C 19-Jun-1990 Hans Joergen Aa. Jensen
C
C Obtain KEEP and IDATA information from LUORDA
C (Moved to here from TRACT2)
C
#include "implicit.h"
C
c#include "dacodes.h"
#include "dalistgs.h"
C
C Used from common blocks:
C  TRINP  : MSYM, MBAS(), INTSYM
C  TRUNIT : LIDATA,IDATA,KEEP, ?
C  INFTRA : IPRTRA
C  INFTAP : LUORDA
C
#include "priunit.h"
#include "trinp.h"
#include "trunit.h"
#include "inftra.h"
#include "inftap.h"
C
      CALL QENTER('TRASET_0')
C
C     RETRIEVE BASE DATA FROM UNIT LUORDA
C        KEEP(1)     = MSYM
C        KEEP(2)     = INTSYM (from revised INTSORT /hjaaj)
C        KEEP(3:10)  = MBAS(1:8)
C        KEEP(11:18) = KEEP(1:8)
C
C        IDATA(ichain,1) = last address for this chain
C        IDATA(ichain,2) = number of records this chain
C        (a value of -2 signals error, this chain is not defined)
C
      IDISKA=0
      CALL DAREAD(LUORDA,KEEP,18,IDISKA)
      CALL DAREAD(LUORDA,IDATA,2*LIDATA,IDISKA)
C
      NERROR = 0
      IF (MSYM .NE. KEEP(1)) NERROR = NERROR + 1
      DO 10 I = 1,8
         IF (MBAS(I) .NE. KEEP(2+I)) NERROR = NERROR + 1
   10 CONTINUE
C
      IF (IPRTRA .GE. 10 .OR. NERROR .GT. 0) THEN
         WRITE (LUPRI,'(/A,2(/A,I5),2(/A,8I5))')
     *      ' KEEP information from AOORDINT',
     *      ' -- KEEP(1)     = MSYM     ',KEEP(1),
     *      ' -- KEEP(2)     = INTSYM   ',KEEP(2),
     *      ' -- KEEP(3:10)  = MBAS(1:8)',(KEEP(I),I=3,10),
     *      ' -- KEEP(11:18) = KEEP(1:8)',(KEEP(I),I=11,18)
      END IF
      IF (NERROR .GT. 0) THEN
         WRITE (LUPRI,'(/A)')
     *      ' ERROR, MSYM and/or MBAS(*) not consistent with AOORDINT'
         CALL QUIT(' ERROR, MSYM and/or MBAS(*) not consistent '//
     *             'with AOORDINT')
      END IF
      IF (IPRTRA .GE. 15) THEN
         WRITE (LUPRI,'(3(/A))')
     *      ' Non-zero IDATA information from AOORDINT :',
     *      ' ========================================',
     *      '     NCHAIN        LADX         NPQ'
         DO 20 I = 1,LIDATA
            IF (IDATA(I,1) .GE. -1 .OR. IDATA(I,2) .GE. 0)
     *         WRITE(LUPRI,'(3I12)') I,IDATA(I,1),IDATA(I,2)
   20    CONTINUE
      END IF
C
      IF (KEEP(2) .NE. 0) THEN
         INTSYM = KEEP(2)
         IF (IPRTRA .GT. 1 .OR. INTSYM .NE. 1) WRITE(LUPRI,'(/A,I8)')
     *   ' Symmetry of integrals from AOORDINT :',INTSYM
      ELSE
         INTSYM = 1
         IF (IPRTRA .GT. 1) WRITE(LUPRI,'(/A)')
     *   ' Integrals assumed to be totally symmetric.'
      END IF
      CALL QEXIT('TRASET_0')
      RETURN
      END
C/*deck  traset_3*/
      SUBROUTINE TRASET_3( LX2X3, LWORK )
C
C Sep 2010 Hans Joergen Aa. Jensen (based on TRASET_2 from spin-orbit)
C
C Define KEEP and IDATA for AO->MO transformation from AOTWOINT
C (instead of from AOORDINT)
C
#include <implicit.h>
C
C Used from common blocks:
C  INFDIM : N2BASM
C  INFTRA : IPRTRA
C  TRINP  : MSYM, MBAS(), INTSYM
C  TRHSO  : KSYMSO
C  TRUNIT : KEEP,IDATA
C
#include <priunit.h>
#include <infdim.h>
#include <inftra.h>
#include <trinp.h>
#include <trhso.h>
#include <trunit.h>
C
      CALL QENTER('TRASET_3')
C
      LX2X3  = 2*600 + 4*600
C     for BUF(600), IBUF(600), and INDX(4,600) for both i*4 and i*8
C        /Bugfix Oct07-hjaaj
C MAERKE: remember to change here when new UN2OUT in her2out.F
C TODO: use info on LUINTA for allocation /Sep 2010 hjaaj
C
      KEEP(1) = MSYM
      INTSYM = 1
      IF (IPRTRA .GT. 1) WRITE(LUPRI,'(/A)')
     *   ' AO 2-el integrals assumed to be totally symmetric.'
      KEEP(2) = INTSYM
      DO 100 I = 1,8
         KEEP( 2+I) = MBAS(I)
         KEEP(10+I) = 0
  100 CONTINUE
C
C     980820-hjaaj: reserve NWRK*N2BASM memory for other
C     buffers in integral transformation.
C
      NWRK = (LWORK/3) / N2BASM
C     allocate 1/3 of memory for other buffers of size N2BASM
      NWRK = MAX(6+NORBMA,NWRK)
C     We need at least ca. this number of buffers (as far as I can count
C     today!) for the integral transformation to work.
C
      NPQ = MAX( LWORK / N2BASM - NWRK , 1 )
      NPQ = MIN( NPQ , N2BASM )
      DO 200 I = 1,LIDATA
         IDATA(I,1) = 0
         IDATA(I,2) = NPQ
  200 CONTINUE
      CALL QEXIT('TRASET_3')
      RETURN
      END
C /*deck  rdao2*/
      SUBROUTINE RDAO2(LUINTA,BUF,IBUF,H2_RSPQ,MPQOFF)
C
C 20-Mar-1990 hjaaj
C
C Reads two-electron integrals from file to buffer.
C LUINTA: File to read
C BUF   : Points to start element in integral buffer.
C IBUF  : Points to start element in label buffer.
C H2_RSPQ : H2RSPQ(r,s,pq), pq=(mpqoff+1:mpqoff+npq); dim (nbasr,nbass,npq)
C MPQOFF: Offset of last element returned to calling routine
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION H2_RSPQ(*), BUF(*), IBUF(4,600)
C
C Used from common blocsk:
C   INFORB : IBAS(8)
C   INFIND : IROW(*),ISAO(*)
C   TRINFO : ISP,ISQ,ISR,ISS,NBQ,NPQ,NBRS
C
#include <maxash.h>
#include <maxorb.h>
#include <nuclei.h>
#include <inforb.h>
#include <infind.h>
#include <trinfo.h>

#include <inftra.h>
C
      INTEGER NBAS_x(8)
C
#include <ibtfun.h>
C
      IBASP = IBAS(ISP)
      IBASQ = IBAS(ISQ)
      IBASR = IBAS(ISR)
      IBASS = IBAS(ISS)
C
      CALL DZERO(H2_RSPQ,NBRS*NPQ)
C
C     LBUF, NIBUF, and NBITS have been read from LUINTA in
C     SIR_INTOPEN_NEWTRA
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      KINT = 1
      KIINT = KINT + LBUF
C
      CALL REWSPL(LUINTA)
C
      CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
 150  CONTINUE
#if defined (VAR_SPLITFILES)         
         JBUF = 608
         CALL READSI(LUINTA,LENINT,BUF(KINT),JBUF)
#else
         CALL READI(LUINTA,LENINT,BUF(KINT))
#endif
         CALL AOLAB4(BUF(KIINT),LBUF,NIBUF,NBITS,IBUF,LENGTH)
         IF (LENGTH .GT. 0) THEN
            DO 200 I = 1, LENGTH
                  JP  = IBUF(1,I)
                  JSP = ISAO(JP)
               IF (JSP .NE. ISP) GO TO 100
                  JQ  = IBUF(2,I)
                  JSQ = ISAO(JQ)
               IF (JSQ .NE. ISQ) GO TO 100
                  JR  = IBUF(3,I)
                  JSR = ISAO(JR)
               IF (JSR .NE. ISR) GO TO 100
                  JS  = IBUF(4,I)
                  JSS = ISAO(JS)
               IF (JSS .NE. ISS) GO TO 100
C
                  MP = JP - IBASP
                  MQ = JQ - IBASQ
                  IF (ISP .EQ. ISQ) THEN
                     MPQ = IROW(MP) + MQ - MPQOFF
                  ELSE
                     MPQ = (MP-1)*NBQ + MQ - MPQOFF
                  END IF
                  IF (MPQ .GT. 0 .AND. MPQ .LE. NPQ) THEN
                     MR = JR - IBASR
                     MS = JS - IBASS
                     IF (ISR .EQ. ISS) THEN
                        MRS = IROW(MR) + MS
                     ELSE
                        MRS = (MR-1)*NBS + MS
                     END IF
                     H2_RSPQ((MPQ-1)*NBRS + MRS) = BUF(I)
                  END IF

  100       CONTINUE

                  JP  = IBUF(3,I)
                  JSP = ISAO(JP)
               IF (JSP .NE. ISP) GO TO 200
                  JQ  = IBUF(4,I)
                  JSQ = ISAO(JQ)
               IF (JSQ .NE. ISQ) GO TO 200
                  JR  = IBUF(1,I)
                  JSR = ISAO(JR)
               IF (JSR .NE. ISR) GO TO 200
                  JS  = IBUF(2,I)
                  JSS = ISAO(JS)
               IF (JSS .NE. ISS) GO TO 200
C
                  MP = JP - IBASP
                  MQ = JQ - IBASQ
                  IF (ISP .EQ. ISQ) THEN
                     MPQ = IROW(MP) + MQ - MPQOFF
                  ELSE
                     MPQ = (MP-1)*NBQ + MQ - MPQOFF
                  END IF
                  IF (MPQ .GT. 0 .AND. MPQ .LE. NPQ) THEN
                     MR = JR - IBASR
                     MS = JS - IBASS
                     IF (ISR .EQ. ISS) THEN
                        MRS = IROW(MR) + MS
                     ELSE
                        MRS = (MR-1)*NBS + MS
                     END IF
                     H2_RSPQ((MPQ-1)*NBRS + MRS) = BUF(I)
                  END IF

  200       CONTINUE
         ELSE IF (LENGTH .LT. 0 ) THEN
            GO TO 300
         END IF
         GO TO 150
C
 300  CONTINUE
      MPQOFF = MPQOFF + NPQ
      RETURN
      END
C  /* Deck tralim */
      SUBROUTINE TRALIM(NTRLVL,ICDTRA,ITRTYP)
C
C 24-Aug-1989 Hans Joergen Aa. Jensen
C
C     Set MO loop limits for the four integral indices
C
C     Used both for Coulomb (Mulliken) integrals : (1 2 | 3 4)
C                 and exchange (Dirac) integrals : <1 2 | 3 4 >
C
C     JTRLVL = MOD(NTRLVL,100)
C     KTRLVL = NTRLVL/100
C
C     KTRLVL = 0 : a designates inactive+active orbitals
C                  g designates general non-frozen orbitals
C     KTRLVL = 1 : a designates          active orbitals
C                  g designates general non-frozen orbitals
C     KTRLVL = 2 : a designates frozen+inactive+active orbitals
C                  g designates all orbitals
C
C     KTRLVL = 1 is used when an AO super-matrix file is available,
C     otherwise use KTRLVL = 0.
C
C     JTRLVL  INDEX 1  INDEX 2  INDEX 3  INDEX 4   comment
C        0       a        a        a        a      CI calc. (0. ord.)
C        1       a        a        a        g      first-order
C        2       a        g        a        g      Sirius (2. ord.)
C        3       a        g        g        g      third-order
C       >3       g        g        g        g      full transf.
C
C     ITRTYP(1:NORBT) = number of integral indices in which this orbital
C     enters (i.e. 0,1,2,3, or 4)
C
#include "implicit.h"
      DIMENSION ICDTRA(*), ITRTYP(*)
C
C Used from common blocks:
C  TRINP  : MFRO,MISH,MASH,MORB
C  INFTRA : IPRTRA
C
#include "priunit.h"
#include "trinp.h"
#include "inftra.h"
#include "maxorb.h"
C
      CHARACTER*35 EXTPRM(3)
      DATA EXTPRM/'inactive + active orbitals         ',
     &            'only active orbitals               ',
     &            'frozen + inactive + active orbitals'/
C
      CALL QENTER('TRALIM')
      JTRLVL = MOD(NTRLVL,100)
      KTRLVL = NTRLVL/100
C
      IF (KTRLVL .LT. 0 .OR. KTRLVL .GT. 2) THEN
         WRITE (LUPRI,*) 'Illegal KTRLVL in TRALIM, KTRLVL =',KTRLVL
         CALL QTRACE(LUPRI)
         CALL QUIT('FATAL ERROR TRALIM: Illegal KTRLVL')
      END IF
      IF (IPRTRA .GE. 2) THEN
         WRITE (LUPRI,'(/A,I5/2A/)')
     &      ' Integral transformation order : ',JTRLVL,
     &      ' Extent of primary space       : ',EXTPRM(KTRLVL+1)
      END IF
C
C     Calculate ITRTYP.
C     Redefine MFRO,MISH,MASH according to requested integral
C     transformation.  If KTRLVL .eq. 1 then this is the same
C     as in Sirius.
C
      CALL IZERO(ITRTYP,MORBT)
      IMOL = 0
      DO 300 ISYM = 1, MSYM
         IF (KTRLVL .EQ. 0) THEN
            MASH(ISYM) = MISH(ISYM) + MASH(ISYM)
            MISH(ISYM) = 0
         ELSE IF (KTRLVL .EQ. 2) THEN
            MASH(ISYM) = MFRO(ISYM) + MISH(ISYM) + MASH(ISYM)
            MFRO(ISYM) = 0
            MISH(ISYM) = 0
         END IF
         IMOG = IMOL + MFRO(ISYM)
         IMOA = IMOL + MFRO(ISYM) + MISH(ISYM)
         NMOA = MASH(ISYM)
         IMOL = IMOL + MORB(ISYM)
C        ... if 3. order or full integral transformation
C            then reset MASH(*) to MORB(*) - MFRO(ISYM)
C            for symmetry quadruplet tests in TR2CTL
C            and for limits on (c d) summation indices in TRA2
         IF (JTRLVL .GT. 2) THEN
            MASH(ISYM) = MORB(ISYM) - MFRO(ISYM)
            MISH(ISYM) = 0
         END IF
C
         NG = MIN(4,JTRLVL)
         DO 240 I = IMOG + 1,IMOL
            ITRTYP(I) = NG
  240    CONTINUE
         DO 250 I = IMOA + 1,IMOA + NMOA
            ITRTYP(I) = 4
  250    CONTINUE
  300 CONTINUE
C
C     Define ICDTRA matrix
C
      ICDTRA(1:MMORBX) = 0

      IDIST = 0
      DO 780 I = 1,MORBT
         IF (ITRTYP(I) .GE. 3) THEN
            IR = (I*I-I)/2
            DO 770 J = 1,I
               IF ( (ITRTYP(I)+ITRTYP(J)).GE.7 ) THEN
                  IDIST = IDIST + 1
                  ICDTRA(IR+J) = IDIST
               END IF
  770       CONTINUE
         END IF
  780 CONTINUE
      IF (IPRTRA .GE. 5) THEN
         WRITE (LUPRI,'(/A)')
     &      ' ICDTRA: list of distributions included on LUMINT:'
         DO 810 I = 1,MORBT
            IR = (I*I-I)/2
            IMAX = MAXVAL(ICDTRA(IR+1:IR+I))
            IF (IMAX .GT. 0) THEN
             WRITE(LUPRI,'(I5,A,(T8,10I7))') I,' :',(ICDTRA(IR+J),J=1,I)
            ELSE
             WRITE(LUPRI,'(I5,A)') I,' : all zero'
            END IF
  810    CONTINUE
      END IF
      CALL QEXIT ('TRALIM')
      RETURN
      END
C  /* Deck tralvl */
      SUBROUTINE TRALVL(ITRSIR,NTRLVL)
C
C 19-Jun-1990 Hans Joergen Aa. Jensen
C
C Translate from Sirius transformation level ITRSIR
C to NTRLVL used in this module.
C
C Copy from TRALIM of definition of NTRLVL.
C     JTRLVL = MOD(NTRLVL,100)
C     KTRLVL = NTRLVL/100
C
C     KTRLVL = 0 : a designates inactive+active orbitals
C                  g designates general non-frozen orbitals
C     KTRLVL = 1 : a designates          active orbitals
C                  g designates general non-frozen orbitals
C     KTRLVL = 2 : a designates frozen+inactive+active orbitals
C                  g designates all orbitals
C
C     KTRLVL = 1 is used when an AO super-matrix file is available,
C     otherwise use KTRLVL = 0.
C
C     JTRLVL  INDEX 1  INDEX 2  INDEX 3  INDEX 4   comment
C        0       a        a        a        a      CI calc. (0. ord.)
C        1       a        a        a        g      first-order
C        2       a        g        a        g      Sirius (2. ord.)
C        3       a        g        g        g      third-order
C        4       g        g        g        g      full transf.
C
C
#include "implicit.h"
C
      CALL QENTER('TRALVL')
      JTRLVL = ABS(ITRSIR)
#if defined (VAR_OLDCODE)
c --- 931201-hjaaj: MOFCFV does not exist any more
      IF (MOFCFV) THEN
         KTRLVL = 0
      ELSE
         KTRLVL = 1
      END IF
#else
      KTRLVL = 1
#endif
      IF (JTRLVL .EQ. 3) THEN
         JTRLVL = 2
         KTRLVL = 1
      ELSE IF (JTRLVL .EQ. 4) THEN
         JTRLVL = 2
         KTRLVL = 0
      ELSE IF (JTRLVL .EQ. 5) THEN
         JTRLVL = 3
         KTRLVL = 0
      ELSE IF (JTRLVL .GE. 6) THEN
         IF (JTRLVL .EQ. 10) KTRLVL = 0
         JTRLVL = 4
      END IF
      NTRLVL = 100*KTRLVL + JTRLVL
      CALL QEXIT('TRALVL')
      RETURN
      END
C  /* Deck trrdao */
      SUBROUTINE TRRDAO(LUORDA,X1,NX1,WRK ,LWRK ,IDA)
C
C 25-Aug-1989 Hans Joergen Aa. Jensen
C
C Interface routine to read ao integrals for tra2
C
#include "implicit.h"
      DIMENSION X1(NX1), WRK(LWRK)
C
#include "iratdef.h"
c#include "dacodes.h"
C
C Used from common blocks
C  TRINP  : IRDAO
C
#include "priunit.h"
#include "trinp.h"
#include "inftra.h"
C
      IF (IRDAO .EQ. 0) THEN
         CALL DAREAD(LUORDA,X1,IRAT*NX1,IDA)
      ELSE IF (IRDAO .EQ. 1) THEN
         CALL RDDER2(LUORDA,WRK,WRK,LWRK,X1,NX1,IDA)
      ELSE IF (IRDAO .EQ. 2) THEN
         CALL RDHSO2(LUORDA,WRK,WRK(600*(1+IRAT) + 1),X1,IDA)
      ELSE IF (IRDAO .EQ. 3) THEN
         CALL RDAO2(LUORDA,WRK,WRK(600*(1+IRAT) + 1),X1,IDA)
      ELSE
         WRITE (LUPRI,'(/A,I10)')
     &      ' TRDDAO ERROR, illegal code: IRDAO =',IRDAO
         CALL QTRACE(LUPRI)
         CALL QUIT('ILLEGAL IRDAO CODE IN TRRDAO')
      END IF
      RETURN
      END
C  /* Deck rdder2 */
      SUBROUTINE RDDER2(a,b,c,d,e,f,g)
C
C 21-Dec-1989 hjaaj
C
#include "implicit.h"
#include "priunit.h"
#include "inftra.h"
      WRITE(LUPRI,*)' ERROR, RDDER2 called but program is linked'
      WRITE(LUPRI,*)' with dummy RDDER2, not the correct one in Abacus'
      CALL QTRACE(LUPRI)
      CALL QUIT('ERROR, Abacus not linked with correct RDDER2 routine')
      END
