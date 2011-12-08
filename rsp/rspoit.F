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
#ifdef UNDEF
===========================================================================
/* Comdeck notes */
940728-hjaaj: defined IRIX (recl=*4 words, not bytes)
#include <single.h>
931004-pj+hjaaj:
OITH1X: this routine written by OV is never called currently. It may be
  called in rspoit instead of OITH1 if sec-sec block is not needed.
  940721-hjaaj: added sec-sec code if DEF,SECSEC from another file
OITH1X: this routine written by OV is never called currently.
  It may be called in rspoit instead of OITH1 if sec-sec block is
  not needed.
===========================================================================
#endif
C  /* Deck rspoit */
      SUBROUTINE RSPOIT(IDAXO,IDAXN,JTRLVL,LSYOIT,OITMAT,
     *                  WRK,KFREE,LFREE,IPROIT)
C
C Copyright 25-Sep-1990 Hans Joergen Aa. Jensen
C
C Input:
C  IDAXO  : file index for integrals to be one-index transformed
C           (IDAXO = 0 means untransformed integrals)
C  JTRLVL : range 0-4, number of general indices in transformed ints.
C  LSYOIT : Symmetry of one-index transformation
C  OITMAT : The transformation matrix
C
C Output:
C  IDAXN  : file index for the transformed integrals
C
#include <implicit.h>
      DIMENSION OITMAT(NORBT,NORBT), WRK(*)
C
C Used from common blocks:
C   INFORB : NORBT,N2ORBX
C
#include <priunit.h>
#include <inforb.h>
#include <infpri.h>
C
      EXTERNAL OITBD
C
      CALL QENTER('RSPOIT')
C
      IF (IPROIT .GE. 10) THEN
         WRITE (LUPRI,'(//A/A)')
     *      ' Output from RSPOIT (one-index transformation)',
     *      ' ---------------------------------------------'
         WRITE (LUPRI,'(/A,I5)') ' Print level :',IPROIT
         CALL GETTIM(TSTRT,WSTRT)
      END IF
C
C     Check input
C
      NERR = 0
      IF (IDAXO .LT. 0) THEN
         NERR = NERR + 1
         WRITE (LUPRI,*) 'RSPOIT: Illegal IDAXO  =',IDAXO
      END IF
      IF (JTRLVL .LT. 0 .OR. JTRLVL .GT. 4) THEN
         NERR = NERR + 1
         WRITE (LUPRI,*) 'RSPOIT: Illegal JTRLVL =',JTRLVL
      END IF
      IF (LSYOIT .LT. 1 .OR. LSYOIT .GT. NSYM) THEN
         NERR = NERR + 1
         WRITE (LUPRI,*) 'RSPOIT: Illegal LSYOIT =',LSYOIT
      END IF
      IF (NERR .GT. 0) THEN
         CALL QTRACE(LUERR)
         CALL QUIT('Input errors in RSPOIT')
      END IF
C
C     Open integral files and assign IDAXN
C
      LRDAX = N2ORBT
      IDAXN = 999
      CALL OITOPN(IDAXO,IDAXN,LUDAXO,LUDAXN,JTRLVO,JTRLVL,
     *            LSYMXO,LSYOIT,LRDAX)
C
      IF (IPROIT .GE. 2) CALL FLSHFO(LUPRI)
C
C     Do the one-index transformation
C
      IF (IDAXO .EQ. 0) THEN
         CALL OIT2M(LUDAXN,JTRLVL,LSYOIT,OITMAT,WRK,KFREE,LFREE,IPROIT)
      ELSE
         CALL OIT2X(IDAXO,LUDAXN,JTRLVL,LSYMXO,LSYOIT,OITMAT,
     *              WRK,KFREE,LFREE,IPROIT)
      END IF
C
      IF (IPROIT .GE. 10) THEN
         CALL GETTIM(TEND,WEND)
         WRITE (LUPRI,'(//A,F12.2,A)')
     *      ' CPU time used in RSPOIT :',TEND-TSTRT
         WRITE (LUPRI,'(A,F12.2,A)')
     *      ' Elapsed time in RSPOIT  :',WEND-WSTRT
      END IF
C
      CALL QEXIT('RSPOIT')
      RETURN
      END
C  /* Deck oitbd */
      BLOCK DATA OITBD
#include <cbdax.h>
      DATA LUDAX /51,52,53,54,55/
      DATA IOPEN, IUSED, LSYDAX, LVLDAX
     *     /MXDAX*0, MXDAX*0, MXDAX*0, MXDAX*0/
      END
C  /* Deck oitopn */
      SUBROUTINE OITOPN(IDAXO,IDAXN,LUDAXO,LUDAXN,JTRLVO,JTRLVN,
     *                  LSYMXO,LSYOIT,LRDAX)
C
C 11-Nov-1990 Hans Joergen Aa. Jensen
C
C If IDAXO .gt. 0 then open old unit.
C    IDAXO .eq. 0 refers to untransformed integrals.
C If IDAXN .gt. 0 then open new unit.
C
#include <implicit.h>
#include <iratdef.h>
C
C Used from common blocks:
C   CBDAX  : *
C   INFORB : MULD2H
C   INFTAP : LUINTM
C
#include <cbdax.h>
#include <priunit.h>
#include <inforb.h>
#include <inftap.h>
C
C
      IDAXN1 = IDAXN
      GO TO 5
         ENTRY OITOPO(IDAXO,LUDAXO,JTRLVO,LSYMXO,LRDAX)
         IDAXN1 = -1
    5 CONTINUE
      IF (IDAXO .GT. MXDAX) THEN
         WRITE (LUPRI,*) 'OITOPN ERROR: IDAXO =',IDAXO
         CALL QTRACE(LUPRI)
         CALL QUIT('OITOPN ERROR: ILLEGAL IDAXO')
      END IF
      IF (IDAXN1 .LE. 0) GO TO 11
      DO 10 I = 1,MXDAX
         IF (IUSED(I) .EQ. 0) THEN
            IDAXN = I
            IUSED(IDAXN) = 1
            GO TO 11
         END IF
   10 CONTINUE
      WRITE (LUPRI,*) 'OITOPN ERROR: No more available units'
      CALL QTRACE(LUPRI)
      CALL QUIT('OITOPN : no more available units')
   11 CONTINUE
C
      IF (IDAXO .GT. 0) THEN
         LUDAXO = LUDAX(IDAXO)
         LSYMXO = LSYDAX(IDAXO)
         JTRLVO = LVLDAX(IDAXO)
         IF (IUSED(IDAXO) .EQ. 0) THEN
            WRITE (LUPRI,*) 'OITOPN ERROR: File with IDAXO =',IDAXO,
     *         ' has not been made.'
            CALL QTRACE(LUPRI)
            CALL QUIT('OITOPN : Old unit does not exist')
         END IF
         IF (IOPEN(IDAXO) .EQ. 0) THEN
            CALL GPOPEN(LUDAX(IDAXO),' ','OLD','DIRECT','UNFORMATTED',
     &                  IRAT*LRDAX,OLDDX)
            IOPEN(IDAXO) = 1
         END IF
      ELSE IF (IDAXO .EQ. 0) THEN
         LUDAXO = LUINTM
         LSYMXO = 1
         JTRLVO = 4
      END IF
C
      IF (IDAXN1 .GT. 0) THEN
         LUDAXN = -9999
         CALL GPOPEN(LUDAXN,' ','NEW','DIRECT','UNFORMATTED',IRAT*LRDAX,
     &               OLDDX)
         LUDAX(IDAXN) = LUDAXN
         IOPEN(IDAXN) = 1
         LSYDAX(IDAXN) = MULD2H(LSYMXO,LSYOIT)
         LVLDAX(IDAXN) = JTRLVN
         IF (IDAXO .GE. 0) THEN
            JTRCHK = MIN(4,JTRLVN+1)
            IF (JTRLVO .LT. JTRCHK) THEN
               WRITE (LUPRI,*) 'OITOPN ERROR: File with IDAXO =',
     *            IDAXO,' has too low level.'
               WRITE (LUPRI,*) ' Old level           :',JTRLVO
               WRITE (LUPRI,*) ' New level requested :',JTRLVN
               CALL QTRACE(LUPRI)
               CALL QUIT('OITOPN error,too low level on old LUDAX file')
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck oitclo */
      SUBROUTINE OITCLO(IDAXN,DISPOS)
C
C 16-Nov-1990 Hans Joergen Aa. Jensen
C
C Close unit with index IDAXN, if DISPOS(1:3) = 'DEL' then
C delete the file.
C
#include <implicit.h>
      CHARACTER*(*) DISPOS
C
C Used from common blocks:
C   CBDAX  : LUDAX(*), IOPEN(*), ISUED(*)
C
#include <cbdax.h>
#include <priunit.h>
C
      IF (IDAXN .LT. 1 .OR. IDAXN .GT. MXDAX) THEN
         WRITE (LUPRI,*) 'OITCLO ERROR: IDAXN =',IDAXN
         CALL QTRACE(LUPRI)
         CALL QUIT('OITCLO ERROR: ILLEGAL IDAXN')
      END IF
      IF (IUSED(IDAXN) .EQ. 0) GO TO 9999
      LUDAXN = LUDAX(IDAXN)
      IF (DISPOS(1:3) .EQ. 'DEL') THEN
         IF (IOPEN(IDAXN) .EQ. 0) THEN
            CALL OITOPO(IDAXN,LUDAXN,JTRLVN,LSYMXN,LRDAX)
C           CALL OITOPO(IDAXO,LUDAXO,JTRLVO,LSYMXO,LRDAX)
         END IF
         CALL GPCLOSE(LUDAXN,'DELETE')
         IUSED(IDAXN) = 0
      ELSE IF (IOPEN(IDAXN) .NE. 0) THEN
         CALL GPCLOSE(LUDAXN,'KEEP')
      END IF
 9999 IOPEN(IDAXN) = 0
      CONTINUE
      RETURN
      END
C  /* Deck oith1 */
      SUBROUTINE OITH1(LSYOIT,OITMAT,H1,H1X,IH1SYM)
C
C Copyright 16-Nov-1990 Hans Joergen Aa. Jensen
C
C One-index transform H1(norbt,norbt) of symmetry IH1SYM
C to H1X(norbt,norbt) using OITMAT(norbt,norbt) of symmetry
C LSYOIT.
C
C Input : H1(a,b) of symmetry IH1SYM
C Output: H1X(a,b) = H1X(a,b)
C                  + H1(a,b) one-index transformed with OITMAT
C                  = H1X(a,b)
C                  + sum(c) [ OITMAT(a,c) H1(c,b)
C                           - H1(a,c) OITMAT(c,b) ]
C
#include <implicit.h>
      DIMENSION OITMAT(NORBT,NORBT), H1(NORBT,NORBT), H1X(NORBT,NORBT)
C
C Used from common blocks:
C  INFORB : NSYM, MULD2H, NORBT, NORB(*), IORB(*)
C
#include <inforb.h>
C
C
C One-index transform both indices and add to previous content.
C
C  (a~ b~) =  SUM(c) [ OITMAT(a,c)*(c b)
C                    - (a c) OITMAT(c,b) ]
C
      DO 100 ICSYM = 1,NSYM
         IBSYM = MULD2H(ICSYM,IH1SYM)
         IASYM = MULD2H(ICSYM,LSYOIT)
         NORBC = NORB(ICSYM)
         NORBB = NORB(IBSYM)
         NORBA = NORB(IASYM)
         IF ((NORBA.NE.0) .AND. (NORBC.NE.0) .AND. (NORBB.NE.0)) THEN
            ICST  = IORB(ICSYM) + 1
            IBST  = IORB(IBSYM) + 1
            IAST  = IORB(IASYM) + 1
            CALL DGEMM('N','N',NORBA,NORBB,NORBC,1.D0,
     &                 OITMAT(IAST,ICST),NORBT,
     &                 H1(ICST,IBST),NORBT,1.D0,
     &                 H1X(IAST,IBST),NORBT)
            CALL DGEMM('N','N',NORBB,NORBA,NORBC,-1.D0,
     &                 H1(IBST,ICST),NORBT,
     &                 OITMAT(ICST,IAST),NORBT,1.D0,
     &                 H1X(IBST,IAST),NORBT)
         END IF
  100 CONTINUE
C
C     End of OITH1
C
      RETURN
      END
C  /* Deck oit2m */
      SUBROUTINE OIT2M (LUDAXN,JTRLVL,LSYOIT,OITMAT,
     *                  WRK,KFRSAV,LFRSAV,IPROIT)
C
C Copyright 15. Nov 1990 by Hans Jorgen Aa. Jensen.
C
C Purpose:
C   Construct one-index transformed 2-electron integrals.
C
C Input:
C   LUDAXN: unit number for output direct access file.
C   JTRLVL: Transformation level of transformed integrals.
C   LSYOIT: Symmetry of OITMAT
C   OITMAT: One-index transformation matrix
C
C Scratch:
C   WRK(KFRSAV:KFRSAV-1+LFRSAV)
C
C *********************************************************************
C
#include <implicit.h>
      DIMENSION OITMAT(N2ORBX), WRK(LFRSAV)
C
C Used from common blocks
C   INFORB : NSYM,N2ORBX,N2ORBT,NNORBX,NORBT,...
C   INFDIM : MWORK,NORBMA,...
C
#include <maxorb.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
C
C     Local variables
C
      DIMENSION NEEDMU(6), N2DIS(8)

      CALL QENTER('OIT2M')
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C     Set NEEDMU array
C     ... occupied-occupied distributions always needed
C     ... unocc-occ distrib. needed for level 1 and higher
C     ... all distributions are needed for level 2 and higher
C
      NEEDMU(1) = 1
      NEEDMU(2) = 1
      NEEDMU(3) = 1
      IF (JTRLVL .GE. 1) THEN
         NEEDMU(4) = 1
         NEEDMU(5) = 1
      ELSE
         NEEDMU(4) = 0
         NEEDMU(5) = 0
      END IF
      IF (JTRLVL .GE. 2) THEN
         NEEDMU(6) = 1
      ELSE
         NEEDMU(6) = 0
      END IF
C
C     Determine if in core or out of core:
C
      MWORK1 = MIN(MWORK,LFREE-10*N2ORBX)
C     ... 10 is an arbitrary number, which is hoped to be sufficient
      NH2XCD = MWORK1 / N2ORBT
      IF (JTRLVL .EQ. 0 .OR. JTRLVL .EQ. 1) THEN
         IF (NH2XCD .GE. N2OCCX) THEN
            ICTOIT = 2
            NH2XCD = N2OCCX
         ELSE
            ICTOIT = 3
            IF (JTRLVL .EQ. 0) THEN
               NH2XCD = NSYM*((NISHMA+NASHMA+1)*(NISHMA+NASHMA))/2
               NH2XCD = MIN(NH2XCD,NNOCCX)
C              MAERKE need NNOCCT but that is not defined yet
            ELSE
               NH2XCD = NSYM*(NISHMA + NASHMA)*NORBMA
               NH2XCD = MIN(NH2XCD,NNOCCX + NOCCT*NSSHT)
C              MAERKE need max of all symmetries but that is not defined
            END IF
         END IF
      ELSE IF (JTRLVL .EQ. 2) THEN
         IF (NH2XCD .GE. N2OCCX + 2*NOCCT*(NORBT-NOCCT) ) THEN
            ICTOIT = 2
            NH2XCD = N2OCCX + 2*NOCCT*(NORBT-NOCCT)
         ELSE
            ICTOIT = 3
            NH2XCD = NNORBT
         END IF
      ELSE
         IF (NH2XCD .GE. NNORBX) THEN
            ICTOIT = 1
            NH2XCD = NNORBX
         ELSE
            ICTOIT = 3
            NH2XCD = NNORBT
         END IF
      END IF
C
C     Allocate work memory
C
      IF (ICTOIT .EQ. 3) THEN
         LH2X  = MIN(N2ORBT,MWORK1/NH2XCD)
         LH2X  = MAX(NORBMA,LH2X)
      ELSE
         LH2X  = N2ORBT
      END IF
      LH2XT = LH2X*NH2XCD
      CALL MEMGET('REAL',KH2CD ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2XCD,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDAB,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDCD,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIN2CD,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KICDTR,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2X  ,LH2XT ,WRK,KFREE,LFREE)
C
C     Open auxiliary DA file
C
      IF (ICTOIT .EQ. 3) THEN
         IDAXO = -1
         IDAXT = 999
         JTRLVT= 4
         LSYMXO= 1
         LRDAX = N2ORBT
         CALL OITOPN(IDAXO,IDAXT,LUDAXO,LUDAXT,JTRLVO,JTRLVT,
     *                     LSYMXO,LSYOIT,LRDAX)
C        CALL OITOPN(IDAXO,IDAXN,LUDAXO,LUDAXN,JTRLVO,JTRLVN,
C    *                     LSYMXO,LSYOIT,LRDAX)
      ELSE
         LUDAXT = -999
      END IF
C
C     Set INDAB and ICDTRA arrays
C
      CALL OITIND(WRK(KINDAB),N2DIS)
      CALL OITICD(JTRLVL,WRK(KICDTR),WRK(KINDCD),IPROIT)
C     CALL OITICD(JTRLVL,ICDTRA,ITRTYP,IPROIT)
C
      CALL OIT2M2(LUDAXN,JTRLVL,ICTOIT,LUDAXT,LSYOIT,
     *            OITMAT,WRK(KH2CD),WRK(KH2XCD),WRK(KICDTR),
     *            NEEDMU,N2DIS,WRK(KINDAB),WRK(KINDCD),WRK(KIN2CD),
     *            WRK(KH2X),LH2X,NH2XCD, WRK,KFREE,LFREE,IPROIT)
C
      IF (ICTOIT .EQ. 3) CALL OITCLO(IDAXT,'DELETE')
C
C *** end of subroutine OIT2M
C
      CALL MEMREL('OIT2M',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('OIT2M')
      RETURN
      END
C  /* Deck oit2m2 */
      SUBROUTINE OIT2M2(LUDAXN,JTRLVL,ICTOIT,LUDAXT,LSYOIT,
     *                  OITMAT,H2CD,H2XCD,
     *                  ICDTRA,NEEDMU,N2DIS,INDAB,INDCD,IN2CD,
     *                  H2X,LH2X,NH2XCD,WRK,KFRSAV,LFRSAV,IPROIT)
C
C Copyright 15. Nov 1990 by Hans Jorgen Aa. Jensen.
C
C Purpose:
C   Construct one-index transformed 2-electron integrals.
C
C Input:
C   LUDAXN: unit number for output direct access file.
C   JTRLVL: Transformation level of transformed integrals.
C   ICTOIT: control parameter for transformation
C   LUDAXT: temporary file for out-of-core transformation, if ICTOIT=3
C   LSYOIT: Symmetry of OITMAT
C   OITMAT: One-index transformation matrix
C   N2DIS : number of distributions of each symmetry
C   INDAB : index for (AB) symmetry packed integrals.
C   NH2XCD: Number of allocated buffers in H2X
C
C Scratch:
C   The rest, including
C   WRK(KFRSAV:KFRSAV-1+LFRSAV)
C
C *********************************************************************
C
#include <implicit.h>
      DIMENSION OITMAT(N2ORBX), H2CD(N2ORBX), H2XCD(N2ORBX)
      DIMENSION H2X(LH2X,NH2XCD), WRK(*), N2DIS(8)
      DIMENSION ICDTRA(NORBT,NORBT), NEEDMU(6), INDAB(NORBT,NORBT)
      DIMENSION INDCD(NORBT,NORBT), IN2CD(NORBT,NORBT)
      PARAMETER ( D1 = 1.0D0 )
C
C Used from common blocks
C   INFORB : NSYM,N2ORBX,N2ORBT,NNORBX,NORBT,...
C   INFIND : ISMO(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      LOGICAL CDEQDC

      CALL QENTER('OIT2M2')
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C     Allocate work memory
C
C
      IF (ICTOIT .LT. 1 .OR. ICTOIT .GT. 3) IPROIT = 5
      IF (IPROIT .GE. 5) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Test output from OIT2M2.'
         WRITE (LUPRI,*) 'LUDAXN :',LUDAXN
         WRITE (LUPRI,*) 'JTRLVL :',JTRLVL
         WRITE (LUPRI,*) 'ICTOIT :',ICTOIT
         WRITE (LUPRI,*) 'LSYOIT :',LSYOIT
         WRITE (LUPRI,*) 'NH2XCD :',NH2XCD
         WRITE (LUPRI,*) 'LH2X   :',LH2X
         WRITE (LUPRI,*) 'N2ORBT :',N2ORBT
         WRITE (LUPRI,*) 'NNORBX :',NNORBX
         WRITE (LUPRI,*) 'N2DIS  :',(N2DIS(I),I=1,NSYM)
         CALL FLSHFO(LUPRI)
      END IF
      IF (ICTOIT .LT. 1 .OR. ICTOIT .GT. 3) THEN
         WRITE (LUPRI,*) ' *** ERROR: Illegal ICTOIT'
         CALL QTRACE(LUPRI)
         CALL QUIT('*** ERROR: Illegal ICTOIT in OIT2M2.')
      END IF
C
C ****************************************************************
C     Loop over Mulliken distributions allowed in NEEDMU(6)
C
      CALL IZERO(INDCD,N2ORBX)
      IF (ICTOIT .EQ. 2) THEN
         CALL DZERO(H2X,NH2XCD*N2ORBT)
      END IF
      JDIST = 0
      IDIST = 0
  100 CALL NXTH2M(IC,ID,H2CD,NEEDMU,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
         JDIST = JDIST + 1
         IF (ICTOIT .EQ. 1 .AND. JDIST .GT. NH2XCD) THEN
            WRITE (LUPRI,*)
     *         'OIT2M2 error, insufficient allocation for H2X'
            WRITE (LUPRI,*) ' -- Allocated :',NH2XCD
            CALL QTRACE(LUPRI)
            CALL QUIT('OIT2M2 error, insufficient allocation for H2X')
         END IF
         INDCD(IC,ID) = JDIST
         INDCD(ID,IC) = JDIST
         IF (IC .LT. ID) THEN
            IA = IC
            IC = ID
            ID = IA
         END IF
C
         ICSYM  = ISMO(IC)
         IDSYM  = ISMO(ID)
         ICDSYM = MULD2H(ICSYM,IDSYM)
         IABSYM = ICDSYM
C
         IF (IPROIT .GE. 40) THEN
           WRITE (LUPRI,'(/A,I5,A,2I5)')
     &        ' Mulliken distribution no.',JDIST,', IC and ID:',IC,ID
           WRITE (LUPRI,*) ' IDIST =',IDIST
           WRITE (LUPRI,*) 'ICDSYM =',ICDSYM
           CALL OUTPUT(H2CD,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         END IF
C
         CALL DZERO(H2XCD,N2ORBX)
         CALL OITH1(LSYOIT,OITMAT,H2CD,H2XCD,IABSYM)
C        CALL OITH1(LSYOIT,OITMAT,H1,H1X,IH1SYM)
C
         IF (IPROIT .GE. 50) THEN
           WRITE (LUPRI,'(/A,I5,A,2I5)')
     &        ' OIT2M2 distribution no.',JDIST,
     &        ', IC and ID:',IC,ID
           WRITE (LUPRI,*) 'One-index transformed in IA and IB:'
           WRITE (LUPRI,*) 'IABSYM =',IABSYM
           WRITE (LUPRI,*) 'LSYOIT =',LSYOIT
           CALL OUTPUT(H2XCD,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         END IF
C
C        We now have the one-index transformed integrals in H2XCD(*,*):
C        (X Y / C D) = (X B / C D) + (A Y / C D)
C        where X = A~1 and Y = B~2
C
         IXYSYM = MULD2H(IABSYM,LSYOIT)
         N2DXY  = N2DIS(IXYSYM)
         IF (ICTOIT .EQ. 1) THEN
C
C           Save half-transformed integrals in H2X
C
            CALL OITPAK(H2XCD,H2X(1,JDIST),IXYSYM,1)
C           CALL OITPAK(H2XCD,H2X,IXYSYM,IWAY)
         ELSE IF (ICTOIT .EQ. 2) THEN
C
C           Distribute half-transformed integrals in H2X
C           using H2CD as temporary storage for packed integrals
C
            CALL OITPAK(H2XCD,H2CD,IXYSYM,1)
            IRECCD = ICDTRA(IC,ID)
            IRECDC = ICDTRA(ID,IC)
            ICDX   = INDAB (IC,ID)
            IDCX   = INDAB (ID,IC)
            IF (IRECCD .GT. 0) THEN
               CALL DAXPY(N2DXY,D1,H2CD,1,H2X(1,IRECCD),1)
            END IF
            IF (IRECDC .GT. 0 .AND. IRECDC .NE. IRECCD) THEN
               CALL DAXPY(N2DXY,D1,H2CD,1,H2X(1,IRECDC),1)
            END IF
            DO 1200 IY = 1,NORBT
               ISYMY = ISMO(IY)
               ISYMX = MULD2H(IXYSYM,ISYMY)
               IXST  = IORB(ISYMX) + 1
               IXEND = IORB(ISYMX) + NORB(ISYMX)
#include <ivdep.h>
               DO 1100 IX = IXST,IXEND
                  IRECXY = ICDTRA(IX,IY)
                  IF (IRECXY .GT. 0) THEN
                     IXYX = INDAB(IX,IY)
                     H2X(ICDX,IRECXY) = H2X(ICDX,IRECXY) + H2CD(IXYX)
                     IF (ICDX .NE. IDCX) H2X(IDCX,IRECXY)
     *                                = H2X(IDCX,IRECXY) + H2CD(IXYX)
                  END IF
 1100          CONTINUE
 1200       CONTINUE
         ELSE
C
C           ICTOIT = 3: out of core
C           Save half-transformed integrals on disk
C
            CALL OITPAK(H2XCD,H2CD,IXYSYM,1)
            WRITE (LUDAXT,REC=JDIST) (H2CD(I),I=1,N2DXY)
         END IF
C
C        Go to 100 to get next needed Mulliken distribution
C
      GO TO 100
C
C     arrive at 800 when finished with all needed Mulliken distributions
C
  800 CONTINUE
      IF (IPROIT .GE. 5) THEN
         WRITE (LUPRI,'(//A/,3(/A,I5))')
     &     ' End of test output of Mulliken distributions from OIT2M2.',
     &     ' Total number of distributions treated  :',JDIST,
     &     ' Total number of distributions (NNORBX) :',NNORBX,
     &     ' Total number allocated in H2X          :',NH2XCD
         CALL FLSHFO(LUPRI)
      END IF
      IF (KFREE .NE. KFRSAV)
     &   CALL MEMREL('OIT2M2-1',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C
C ****************************************************************
C
C  Now H2X(ij,kl) = (it jt / k l),
C  finish H2X by symmetrizing
C  (remember (it jt / kt lt) = (it jt / k l) + (kt lt / i j) ).
C
C  If requested, print H2X
C
C
      IF (ICTOIT .EQ. 1 .OR. ICTOIT .EQ. 3) THEN
         CALL ICOPY(N2ORBX,INDCD,1,IN2CD,1)
      END IF
      ISYH2X = LSYOIT
      CDEQDC = .TRUE.
      DO 2800 ISYMCD = 1,NSYM
         ISYMAB = MULD2H(ISYH2X,ISYMCD)
         N2DAB  = N2DIS(ISYMAB)
         IF (N2DAB .EQ. 0) GO TO 2800
         IDEND  = 0
 2000 CONTINUE
         IDST   = IDEND + 1
         IF (ICTOIT .EQ. 3) THEN
            CALL OITCOR(ISYMAB,ISYMCD,IDST,IDEND,ICDXOF,CDEQDC,
     *                  INDCD,IN2CD,INDAB,LUDAXT,H2XCD,H2X,LH2X,NH2XCD)
         ELSE
            ICDXOF = 0
            IDEND  = NORBT
         END IF
      DO 2400 ID = IDST,IDEND
         ISYMD  = ISMO(ID)
         ISYMC  = MULD2H(ISYMD,ISYMCD)
         ICST   = IORB(ISYMC) + 1
         ICEND  = IORB(ISYMC) + NORB(ISYMC)
         DO 2300 IC = ICST,ICEND
            IRECCD = ICDTRA(IC,ID)
         IF (IRECCD .EQ. 0) GO TO 2300
         IF (ICTOIT .EQ. 1 .OR. ICTOIT .EQ. 3) THEN
            ICDX   = INDAB(IC,ID)
            ICDDIS = INDCD(IC,ID)
            IF (ICDDIS .EQ. 0) THEN
               CALL QUIT('ERROR OIT2M2,INDCD .eq. 0 when ICDTRA .ne. 0')
            END IF
            IF (ICTOIT .EQ. 1) THEN
               CALL DCOPY(N2DAB,H2X(1,ICDDIS),1,H2XCD,1)
            ELSE
               IF (ICDX .LE. ICDXOF) THEN
                  CALL QUIT('ERROR OIT2M2, ICDX .le. '//
     &                      'ICDXOF for ICTOIT = 3')
               END IF
               READ (LUDAXT,REC=ICDDIS) (H2XCD(I),I=1,N2DAB)
            END IF
            DO 2200 IB = 1,NORBT
               ISYMB = ISMO(IB)
               ISYMA = MULD2H(ISYMB,ISYMAB)
               IAST  = IORB(ISYMA) + 1
               IAEND = IORB(ISYMA) + NORB(ISYMA)
               DO 2100 IA = IAST,IAEND
                  IABDIS = IN2CD(IA,IB)
                  IF (IABDIS .EQ. 0) GO TO 2100
                  IABX   = INDAB(IA,IB)
                  IF (IABDIS .GT. 0) THEN
                     H2XCD(IABX) = H2XCD(IABX) + H2X(ICDX-ICDXOF,IABDIS)
                  ELSE
                     READ(LUDAXT,REC=-IABDIS) (H2CD(I),I=1,ICDX)
                     H2XCD(IABX) = H2XCD(IABX) + H2CD(ICDX)
                  END IF
 2100          CONTINUE
 2200       CONTINUE
         END IF
            IF (IPROIT .GE. 40) THEN
               WRITE (LUPRI,'(/A,I8,A,2I5)')
     &        ' OIT2M2 record no.',IRECCD,', IC and ID:',IC,ID
               WRITE (LUPRI,*) 'ISYMCD =',ISYMCD
               WRITE (LUPRI,*) 'ISYMAB =',ISYMAB
               WRITE (LUPRI,*) 'LSYOIT =',LSYOIT
               IF (ICTOIT .EQ. 1 .OR. ICTOIT .EQ. 3) THEN
                  CALL OITPAK(H2CD,H2XCD,ISYMAB,-1)
               ELSE
                  CALL OITPAK(H2CD,H2X(1,IRECCD),ISYMAB,-1)
               END IF
               CALL OUTPUT(H2CD,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
            END IF
            IF (ICTOIT .EQ. 1 .OR. ICTOIT .EQ. 3) THEN
               WRITE (LUDAXN,REC=IRECCD) (H2XCD(I),I=1,N2DAB)
            ELSE
               WRITE (LUDAXN,REC=IRECCD) (H2X(I,IRECCD),I=1,N2DAB)
            END IF
 2300    CONTINUE
 2400 CONTINUE
         IF (IDEND .LT. NORBT) GO TO 2000
 2800 CONTINUE
C
      IF (KFREE .NE. KFRSAV)
     &   CALL MEMREL('OIT2M2-2',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C *** end of subroutine OIT2M2
C
 9999 CALL QEXIT('OIT2M2')
      RETURN
      END
C  /* Deck oit2x */
      SUBROUTINE OIT2X (IDAXO,LUDAXN,JTRLVL,LSYMXO,LSYOIT,OITMAT,
     *                  WRK,KFRSAV,LFRSAV,IPROIT)
C
C Copyright 21. Nov 1990 by Hans Jorgen Aa. Jensen.
C
C Purpose:
C   Construct one-index transformed 2-electron integrals.
C
C Input:
C   IDAXO : Identification of input direct access file
C   LUDAXN: unit number for output direct access file.
C   JTRLVL: Transformation level of transformed integrals.
C   LSYMXO: Symmetry of "old" integrals
C   LSYOIT: Symmetry of OITMAT
C   OITMAT: One-index transformation matrix
C
C Scratch:
C   WRK(KFRSAV:KFRSAV-1+LFRSAV)
C
C *********************************************************************
C
#include <implicit.h>
      DIMENSION OITMAT(N2ORBX), WRK(LFRSAV)
C
C Used from common blocks
C   INFORB : NSYM,N2ORBX,N2ORBT,NNORBX,NORBT,...
C   INFDIM : MWORK, NORBMA, ...
C
#include <maxorb.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
C
C     Local variables
C
      DIMENSION NEEDMU(6), N2DIS(8)

      CALL QENTER('OIT2X')
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C     Set NEEDMU array
C     ... occupied-occupied distributions always needed
C     ... unocc-occ distrib. needed for level 1 and higher
C     ... all distributions are needed for level 2 and higher
C
      NEEDMU(1) = 1
      NEEDMU(2) = 1
      NEEDMU(3) = 1
      IF (JTRLVL .GE. 1) THEN
         NEEDMU(4) = 1
         NEEDMU(5) = 1
      ELSE
         NEEDMU(4) = 0
         NEEDMU(5) = 0
      END IF
      IF (JTRLVL .GE. 2) THEN
         NEEDMU(6) = 1
      ELSE
         NEEDMU(6) = 0
      END IF
C
C     Determine if in core or out of core
C
      MWORK1 = MIN(MWORK,LFREE-10*N2ORBX)
C     ... 10 is an arbitrary number, which is hoped to be sufficient
      NH2XCD = MWORK1 / N2ORBT
      IF (JTRLVL .EQ. 0 .OR. JTRLVL .EQ. 1) THEN
         IF (NH2XCD .GE. N2OCCX) THEN
            ICTOIT = 2
            NH2XCD = N2OCCX
         ELSE
            ICTOIT = 3
            IF (JTRLVL .EQ. 0) THEN
               NH2XCD = NSYM*(NISHMA+NASHMA)*(NISHMA+NASHMA)
               NH2XCD = MIN(NH2XCD,N2OCCX)
C              MAERKE need N2OCCT but that is not defined
            ELSE
               NH2XCD = NSYM*(NISHMA + NASHMA)*(NSSHMA + NORBMA)
               NH2XCD = MIN(NH2XCD,N2OCCX + 2*NOCCT*NSSHT)
C              MAERKE need max of all symmetries but that is not defined
            END IF
         END IF
      ELSE IF (JTRLVL .EQ. 2) THEN
         IF (NH2XCD .GE. N2OCCX + 2*NOCCT*(NORBT-NOCCT)) THEN
            ICTOIT = 2
            NH2XCD = N2OCCX + 2*NOCCT*(NORBT-NOCCT)
         ELSE
            ICTOIT = 3
            NH2XCD = N2ORBT
         END IF
      ELSE
         IF (NH2XCD .GE. N2ORBX) THEN
            ICTOIT = 1
            NH2XCD = N2ORBX
         ELSE
            ICTOIT = 3
            NH2XCD = N2ORBT
         END IF
      END IF
C
C     Allocate work memory
C
      IF (ICTOIT .EQ. 3) THEN
         LH2X  = MIN(N2ORBT,MWORK1/NH2XCD)
         LH2X  = MAX(NORBMA,LH2X)
      ELSE
         LH2X  = N2ORBT
      END IF
      LH2XT = LH2X*NH2XCD
      CALL MEMGET('REAL',KH2CD ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2XCD,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDAB,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDCD,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIN2CD,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KICDTR,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2X  ,LH2XT ,WRK,KFREE,LFREE)
C
C     Open auxiliary DA file
C
      IF (ICTOIT .EQ. 3) THEN
         IDAXX = -1
         IDAXT = 999
         JTRLVT= 4
         LRDAX = N2ORBT
         CALL OITOPN(IDAXX,IDAXT,LUDAXX,LUDAXT,JTRLVX,JTRLVT,
     *                     LSYMXO,LSYOIT,LRDAX)
C        CALL OITOPN(IDAXO,IDAXN,LUDAXO,LUDAXN,JTRLVO,JTRLVN,
C    *                     LSYMXO,LSYOIT,LRDAX)
      ELSE
         LUDAXT = -999
      END IF
C
C     Set INDAB and ICDTRA arrays
C
      CALL OITIND(WRK(KINDAB),N2DIS)
      CALL OITICD(JTRLVL,WRK(KICDTR),WRK(KINDCD),IPROIT)
C     CALL OITICD(JTRLVL,ICDTRA,ITRTYP,IPROIT)
      CALL OIT2X2(IDAXO,LUDAXN,JTRLVL,ICTOIT,LUDAXT,LSYMXO,LSYOIT,
     *            OITMAT,WRK(KH2CD),WRK(KH2XCD), WRK(KICDTR),
     *            NEEDMU,N2DIS,WRK(KINDAB),WRK(KINDCD),WRK(KIN2CD),
     *            WRK(KH2X),LH2X,NH2XCD, WRK,KFREE,LFREE,IPROIT)
C
      IF (ICTOIT .EQ. 3) CALL OITCLO(IDAXT,'DELETE')
C
C *** end of subroutine OIT2X
C
      CALL MEMREL('OIT2X',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('OIT2X')
      RETURN
      END
C  /* Deck oit2x2 */
      SUBROUTINE OIT2X2(IDAXO,LUDAXN,JTRLVL,ICTOIT,LUDAXT,LSYMXO,LSYOIT,
     *                  OITMAT,H2CD,H2XCD,
     *                  ICDTRA,NEEDMU,N2DIS,INDAB,INDCD,IN2CD,
     *                  H2X,LH2X,NH2XCD,WRK,KFRSAV,LFRSAV,IPROIT)
C
C Copyright 15. Nov 1990 by Hans Jorgen Aa. Jensen.
C
C Purpose:
C   Construct one-index transformed 2-electron integrals.
C
C Input:
C   IDAXO : Identification of "old" integrals
C   LUDAXN: unit number for output direct access file.
C   JTRLVL: Transformation level of transformed integrals.
C   LSYMXO: Symmetry of "old" integrals
C   LSYOIT: Symmetry of OITMAT
C   ICTOIT: Control parameter for transformation
C   LUDAXT: temporary file for out-of-core transformation, if ICTOIT=3
C   OITMAT: One-index transformation matrix
C   N2DIS : number of distributions of each symmetry
C   INDAB : index for (AB) symmetry packed integrals.
C   NH2XCD: Number of allocated buffers in H2X
C
C Scratch:
C   The rest, including
C   WRK(KFRSAV:KFRSAV-1+LFRSAV)
C
C *********************************************************************
C
#include <implicit.h>
      DIMENSION OITMAT(N2ORBX), H2CD(N2ORBX), H2XCD(N2ORBX)
      DIMENSION H2X(LH2X,NH2XCD), WRK(*), N2DIS(8)
      DIMENSION ICDTRA(NORBT,NORBT), NEEDMU(6), INDAB(NORBT,NORBT)
      DIMENSION INDCD(NORBT,NORBT), IN2CD(NORBT,NORBT)
      PARAMETER ( D1 = 1.0D0 )
C
C Used from common blocks
C   INFORB : NSYM,N2ORBX,N2ORBT,NNORBX,NORBT,...
C   INFIND : ISMO(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>

      LOGICAL CDEQDC

      CALL QENTER('OIT2X2')
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C     Allocate work memory
C
C
      IF (ICTOIT .LT. 1 .OR. ICTOIT .GT. 3) IPROIT = 5
      IF (IPROIT .GE. 5) THEN
         WRITE (LUPRI,'(///A)')
     &      ' Test output from OIT2X2.'
         WRITE (LUPRI,*) ' IDAXO :',IDAXO
         WRITE (LUPRI,*) 'LUDAXN :',LUDAXN
         WRITE (LUPRI,*) 'JTRLVL :',JTRLVL
         WRITE (LUPRI,*) 'LSYMXO :',LSYMXO
         WRITE (LUPRI,*) 'LSYOIT :',LSYOIT
         WRITE (LUPRI,*) 'ICTOIT :',ICTOIT
         WRITE (LUPRI,*) 'NH2XCD :',NH2XCD
         WRITE (LUPRI,*) 'LH2X   :',LH2X
         WRITE (LUPRI,*) 'N2ORBT :',N2ORBT
         WRITE (LUPRI,*) 'N2OCCX :',N2OCCX
         WRITE (LUPRI,*) 'N2ORBX :',N2ORBX
         WRITE (LUPRI,*) 'N2DIS  :',(N2DIS(I),I=1,NSYM)
         CALL FLSHFO(LUPRI)
      END IF
      IF (ICTOIT .LT. 1 .OR. ICTOIT .GT. 3) THEN
         WRITE (LUPRI,*) ' *** ERROR: Illegal ICTOIT'
         CALL QTRACE(LUPRI)
         CALL QUIT('*** ERROR: Illegal ICTOIT in OIT2X2')
      END IF
C
C ****************************************************************
C     Loop over Mulliken distributions allowed in NEEDMU(6)
C
      CALL IZERO(INDCD,N2ORBX)
      IF (ICTOIT .EQ. 2) THEN
         CALL DZERO(H2X,NH2XCD*N2ORBT)
      END IF
      JDIST = 0
      IDIST = 0
  100 CALL NXTH2X(IC,ID,H2CD,IDAXO,NEEDMU,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
         JDIST = JDIST + 1
         IF (ICTOIT .EQ. 1 .AND. JDIST .GT. NH2XCD) THEN
            WRITE (LUPRI,*)
     *         'OIT2X2 error, insufficient allocation for H2X'
            WRITE (LUPRI,*) ' -- Allocated :',NH2XCD
            CALL QTRACE(LUPRI)
            CALL QUIT('OIT2X2 error, insufficient allocation for H2X')
         END IF
         INDCD(IC,ID) = JDIST
C
         ICSYM  = ISMO(IC)
         IDSYM  = ISMO(ID)
         ICDSYM = MULD2H(ICSYM,IDSYM)
         IABSYM = MULD2H(LSYMXO,ICDSYM)
C
         IF (IPROIT .GE. 40) THEN
           WRITE (LUPRI,'(/A,I5,A,2I5)')
     &        ' OIT2X2 Distribution no.',JDIST,', IC and ID:',IC,ID
           WRITE (LUPRI,*) ' IDIST =',IDIST
           WRITE (LUPRI,*) 'ICDSYM =',ICDSYM
           WRITE (LUPRI,*) 'IABSYM =',IABSYM
           CALL OUTPUT(H2CD,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         END IF
C
         CALL DZERO(H2XCD,N2ORBX)
         CALL OITH1(LSYOIT,OITMAT,H2CD,H2XCD,IABSYM)
C        CALL OITH1(LSYOIT,OITMAT,H1,H1X,IH1SYM)
C
         IF (IPROIT .GE. 50) THEN
           WRITE (LUPRI,'(/A,I5,A,2I5)')
     &        ' OIT2X2 distribution no.',JDIST,', IC and ID:',IC,ID
           WRITE (LUPRI,*) 'One-index transformed in IA and IB:'
           WRITE (LUPRI,*) 'IABSYM =',IABSYM
           WRITE (LUPRI,*) 'LSYOIT =',LSYOIT
           CALL OUTPUT(H2XCD,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
         END IF
C
C        We now have the one-index transformed integrals in H2XCD(*,*):
C        (X Y / C D) = (X B / C D) + (A Y / C D)
C        where X = A~1 and Y = B~2
C
         IXYSYM = MULD2H(IABSYM,LSYOIT)
         N2DXY  = N2DIS(IXYSYM)
         IF (ICTOIT .EQ. 1) THEN
C
C           Save half-transformed integrals in H2X
C
            CALL OITPAK(H2XCD,H2X(1,JDIST),IXYSYM,1)
C           CALL OITPAK(H2XCD,H2X,IXYSYM,IWAY)
         ELSE IF (ICTOIT .EQ. 2) THEN
C
C           Distribute half-transformed integrals in H2X
C           using H2CD as temporary storage for packed integrals
C
            CALL OITPAK(H2XCD,H2CD,IXYSYM,1)
            IRECCD = ICDTRA(IC,ID)
            ICDX   = INDAB (IC,ID)
            IF (IRECCD .GT. 0) THEN
               CALL DAXPY(N2DXY,D1,H2CD,1,H2X(1,IRECCD),1)
            END IF
            DO 1200 IY = 1,NORBT
               ISYMY = ISMO(IY)
               ISYMX = MULD2H(IXYSYM,ISYMY)
               IXST  = IORB(ISYMX) + 1
               IXEND = IORB(ISYMX) + NORB(ISYMX)
#include <ivdep.h>
               DO 1100 IX = IXST,IXEND
                  IRECXY = ICDTRA(IX,IY)
                  IF (IRECXY .GT. 0) THEN
                     IXYX = INDAB(IX,IY)
                     H2X(ICDX,IRECXY) = H2X(ICDX,IRECXY) + H2CD(IXYX)
                  END IF
 1100          CONTINUE
 1200       CONTINUE
         ELSE
C
C           ICTOIT = 3: out of core
C           Save half-transformed integrals on disk
C
            CALL OITPAK(H2XCD,H2CD,IXYSYM,1)
            WRITE (LUDAXT,REC=JDIST) (H2CD(I),I=1,N2DXY)
         END IF
C
C        Go to 100 to get next needed Mulliken distribution
C
      GO TO 100
C
C     arrive at 800 when finished with all needed Mulliken distributions
C
  800 CONTINUE
      IF (IPROIT .GE. 5) THEN
         WRITE (LUPRI,'(//A/,3(/A,I5))')
     &     ' End of test output of input distributions from OIT2X2.',
     &     ' Total number of distributions treated  :',JDIST,
     &     ' Total number of distributions (N2ORBX) :',N2ORBX,
     &     ' Total number allocated in H2X          :',NH2XCD
         CALL FLSHFO(LUPRI)
      END IF
      IF (KFREE .NE. KFRSAV)
     &   CALL MEMREL('OIT2X2-1',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C
C ****************************************************************
C
C  Now H2X(ij,kl) = (it jt / k l),
C  finish H2X by symmetrizing
C  (remember (it jt / kt lt) = (it jt / k l) + (kt lt / i j) ).
C
C  If requested, print H2X
C
C
      IF (ICTOIT .EQ. 1 .OR. ICTOIT .EQ. 3) THEN
         CALL ICOPY(N2ORBX,INDCD,1,IN2CD,1)
      END IF
      ISYH2X = MULD2H(LSYOIT,LSYMXO)
      CDEQDC = .FALSE.
      DO 2800 ISYMCD = 1,NSYM
         ISYMAB = MULD2H(ISYH2X,ISYMCD)
         N2DAB  = N2DIS(ISYMAB)
         IF (N2DAB .EQ. 0) GO TO 2800
         IDEND  = 0
 2000 CONTINUE
         IDST   = IDEND + 1
         IF (ICTOIT .EQ. 3) THEN
            CALL OITCOR(ISYMAB,ISYMCD,IDST,IDEND,ICDXOF,CDEQDC,
     *                  INDCD,IN2CD,INDAB,LUDAXT,H2XCD,H2X,LH2X,NH2XCD)
         ELSE
            ICDXOF = 0
            IDEND  = NORBT
         END IF
      DO 2400 ID = IDST,IDEND
         ISYMD  = ISMO(ID)
         ISYMC  = MULD2H(ISYMD,ISYMCD)
         ICST   = IORB(ISYMC) + 1
         ICEND  = IORB(ISYMC) + NORB(ISYMC)
         DO 2300 IC = ICST,ICEND
            IRECCD = ICDTRA(IC,ID)
         IF (IRECCD .EQ. 0) GO TO 2300
         IF (ICTOIT .EQ. 1 .OR. ICTOIT .EQ. 3) THEN
            ICDX   = INDAB(IC,ID)
            ICDDIS = INDCD(IC,ID)
            IF (ICDDIS .EQ. 0) THEN
               CALL QUIT('ERROR OIT2X2, INDCD .eq. 0 when ICDTRA .ne.0')
            END IF
            IF (ICTOIT .EQ. 1) THEN
               CALL DCOPY(N2DAB,H2X(1,ICDDIS),1,H2XCD,1)
            ELSE
               IF (ICDX .LE. ICDXOF) THEN
                  CALL QUIT('ERROR OIT2X2, ICDX .le. '//
     &                      'ICDXOF for ICTOIT = 3')
               END IF
               READ (LUDAXT,REC=ICDDIS) (H2XCD(I),I=1,N2DAB)
            END IF
            DO 2200 IB = 1,NORBT
               ISYMB = ISMO(IB)
               ISYMA = MULD2H(ISYMB,ISYMAB)
               IAST  = IORB(ISYMA) + 1
               IAEND = IORB(ISYMA) + NORB(ISYMA)
               DO 2100 IA = IAST,IAEND
                  IABDIS = IN2CD(IA,IB)
                  IF (IABDIS .EQ. 0) GO TO 2100
                  IABX   = INDAB(IA,IB)
                  IF (IABDIS .GT. 0) THEN
                     H2XCD(IABX) = H2XCD(IABX) + H2X(ICDX-ICDXOF,IABDIS)
                  ELSE
                     READ(LUDAXT,REC=-IABDIS) (H2CD(I),I=1,ICDX)
                     H2XCD(IABX) = H2XCD(IABX) + H2CD(ICDX)
                  END IF
 2100          CONTINUE
 2200       CONTINUE
         END IF
            IF (IPROIT .GE. 40) THEN
               WRITE (LUPRI,'(/A,I8,A,2I5)')
     &        ' OIT2X2 record no.',IRECCD,', IC and ID:',IC,ID
               WRITE (LUPRI,*) 'ISYMCD =',ISYMCD
               WRITE (LUPRI,*) 'ISYMAB =',ISYMAB
               WRITE (LUPRI,*) 'LSYOIT =',LSYOIT
               IF (ICTOIT .EQ. 1 .OR. ICTOIT .EQ. 3) THEN
                  CALL OITPAK(H2CD,H2XCD,ISYMAB,-1)
               ELSE
                  CALL OITPAK(H2CD,H2X(1,IRECCD),ISYMAB,-1)
               END IF
               CALL OUTPUT(H2CD,1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
            END IF
            IF (ICTOIT .EQ. 1 .OR. ICTOIT .EQ. 3) THEN
               WRITE (LUDAXN,REC=IRECCD) (H2XCD(I),I=1,N2DAB)
            ELSE
               WRITE (LUDAXN,REC=IRECCD) (H2X(I,IRECCD),I=1,N2DAB)
            END IF
 2300    CONTINUE
 2400 CONTINUE
         IF (IDEND .LT. NORBT) GO TO 2000
 2800 CONTINUE
C
      IF (KFREE .NE. KFRSAV)
     &   CALL MEMREL('OIT2X2-2',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C *** end of subroutine OIT2X2
C
 9999 CALL QEXIT('OIT2X2')
      RETURN
      END
C  /* Deck oitpak */
      SUBROUTINE OITPAK(H2XCD,H2X,IXYSYM,IWAY)
C
C Copyright 16-Nov-1990 Hans Joergen Aa. Jensen
C
C IWAY .ge. 0 : Pack H2XCD in H2X
C      else   : Unpack H2X in H2XCD
C
#include <implicit.h>
      DIMENSION H2XCD(NORBT,NORBT), H2X(N2ORBT)
C
C Used from common blocks:
C   INFORB : NORBT,N2ORBT,N2ORBX,NSYM,MULD2H
C
#include <inforb.h>
C
      IF (IWAY .LT. 0) CALL DZERO(H2XCD,N2ORBX)
      ISTH2X = 1
      DO 300 IBSYM = 1,NSYM
         IASYM  = MULD2H(IBSYM,IXYSYM)
         NORBA  = NORB(IASYM)
         NORBB  = NORB(IBSYM)
         IAST   = IORB(IASYM) + 1
         IBST   = IORB(IBSYM) + 1
         IF (IWAY .GE. 0) THEN
            CALL MCOPY(NORBA,NORBB,H2XCD(IAST,IBST),NORBT,
     *                 H2X(ISTH2X),NORBA)
         ELSE
            CALL MCOPY(NORBA,NORBB,H2X(ISTH2X),NORBA,
     *                 H2XCD(IAST,IBST),NORBT)
         END IF
         ISTH2X = ISTH2X + NORBA*NORBB
  300 CONTINUE
C
C         MCOPY(nrowa,ncola,A,nrdima,B,nrdimb)
C
      RETURN
      END
C  /* Deck oitind */
      SUBROUTINE OITIND(INDAB,N2DIS)
C
C Copyright 16-Nov-1990 Hans Joergen Aa. Jensen
C
C Find index for integral (AB) in symmetry packed integral list
C
#include <implicit.h>
      DIMENSION INDAB(NORBT,NORBT), N2DIS(8)
C
C Used from common blocks:
C   INFORB : NORBT,N2ORBT,N2ORBX,NSYM,MULD2H
C
#include <inforb.h>
C
      LOGICAL NOIND
C
      NOIND = .FALSE.
      GO TO 10
         ENTRY OITDIS(N2DIS)
         NOIND = .TRUE.
   10 CONTINUE
C
      DO 400 IABSYM = 1,NSYM
         INDXAB = 0
         DO 300 IBSYM = 1,NSYM
            IASYM  = MULD2H(IBSYM,IABSYM)
            IF (NOIND) THEN
               INDXAB = INDXAB + NORB(IASYM)*NORB(IBSYM)
            ELSE
               IAST   = IORB(IASYM) + 1
               IAEND  = IORB(IASYM) + NORB(IASYM)
               IBST   = IORB(IBSYM) + 1
               IBEND  = IORB(IBSYM) + NORB(IBSYM)
               DO 200 IB = IBST,IBEND
                  DO 100 IA = IAST,IAEND
                     INDXAB = INDXAB + 1
                     INDAB(IA,IB) = INDXAB
  100             CONTINUE
  200          CONTINUE
            END IF
  300    CONTINUE
         N2DIS(IABSYM) = INDXAB
  400 CONTINUE
C
      RETURN
      END
C  /* Deck oiticd */
      SUBROUTINE OITICD(JTRLVL,ICDTRA,ITRTYP,IPROIT)
C
C Copyright 19-Nov-1990 Hans Joergen Aa. Jensen
C
C Find file index for (C D) distribution of final integrals.
C
C     JTRLVL = general transformation level
C     KTRLVL = 2 in this version
C
C     KTRLVL = 0 : a designates inactive+active orbitals
C                  g designates general non-frozen orbitals
C     KTRLVL = 1 : a designates          active orbitals
C                  g designates general non-frozen orbitals
C     KTRLVL = 2 : a designates frozen+inactive+active orbitals
C                  g designates all orbitals
C
C     JTRLVL  INDEX 1  INDEX 2  INDEX 3  INDEX 4   comment
C        0       a        a        a        a      CI calc. (0. ord.)
C        1       a        a        a        g      first-order
C        2       a        a        g        g      second-order
C                a        g        a        g
C        3       a        g        g        g      third-order
C                g        g        a        g
C        4       g        g        g        g      full transf.
C     and permutations between indices 1 and 2 and between 3 and 4.
C
C     Indices 1 and 2 define the distributions needed (stored in
C     ICDTRA(1:NORBT,1:NORBT)).
C     The (g g) distributions are needed at level 3 in order to perform
C     a one-index transformation leading to (aa/gg) integrals at level 2
C
C     ITRTYP(1:NORBT) = number of integral indices in which this orbital
C     enters (i.e. 0,1,2,3, or 4)
C
#include <implicit.h>
      DIMENSION ICDTRA(NORBT,NORBT), ITRTYP(NORBT)
      PARAMETER (KTRLVL = 2)
C
C Used from common blocks:
C  INFORB : NSYM,NORBT,...
C  INFIND : ISMO(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
      CHARACTER*35 EXTPRM(3)
      DATA EXTPRM/'inactive + active orbitals         ',
     &            'only active orbitals               ',
     &            'frozen + inactive + active orbitals'/
C
      CALL QENTER('OITICD')
C
      IF (KTRLVL .LT. 0 .OR. KTRLVL .GT. 2) THEN
         WRITE (LUPRI,*) 'Illegal KTRLVL in OITICD, KTRLVL =',KTRLVL
         CALL QTRACE(LUPRI)
         CALL QUIT('FATAL ERROR OITICD: Illegal KTRLVL')
      END IF
      IF (IPROIT .GE. 10) THEN
         WRITE (LUPRI,'(/1X,A,I5/1X,2A/)')
     &      'Integral transformation order : ',JTRLVL,
     &      'Extent of primary space       : ',EXTPRM(KTRLVL+1)
      END IF
C
C     Calculate ITRTYP.
C
      CALL IZERO(ITRTYP,NORBT)
      DO 300 ISYM = 1, NSYM
         IF (KTRLVL .EQ. 0) THEN
            IMOG = IORB(ISYM) + NFRO(ISYM) + 1
            IMOA = IMOG
            NMOA = NISH(ISYM) + NASH(ISYM)
         ELSE IF (KTRLVL .EQ. 1) THEN
            IMOG = IORB(ISYM) + NFRO(ISYM) + 1
            IMOA = IMOG + NISH(ISYM)
            NMOA = NASH(ISYM)
         ELSE IF (KTRLVL .EQ. 2) THEN
            IMOG = IORB(ISYM) + 1
            IMOA = IMOG
            NMOA = NFRO(ISYM) + NISH(ISYM) + NASH(ISYM)
         END IF
         IMOL = IORB(ISYM) + NORB(ISYM)
C
         NG = MIN(4,JTRLVL)
         DO 240 I = IMOG ,IMOL
            ITRTYP(I) = NG
  240    CONTINUE
         DO 250 I = IMOA ,IMOA - 1 + NMOA
            ITRTYP(I) = 4
  250    CONTINUE
  300 CONTINUE
C
C     Determine ICDTRA matrix
C
      CALL IZERO(ICDTRA,N2ORBX)
      IDIST = 0
      DO 790 IJSYM = 1,NSYM
      DO 780 I = 1,NORBT
         IF (ITRTYP(I) .GE. 2) THEN
            ISYM = ISMO(I)
            JSYM = MULD2H(ISYM,IJSYM)
            JST  = IORB(JSYM) + 1
            JEND = IORB(JSYM) + NORB(JSYM)
            DO 770 J = JST,JEND
               IF ( (ITRTYP(I)+ITRTYP(J)) .GE. 6 ) THEN
                  IDIST = IDIST + 1
                  ICDTRA(J,I) = IDIST
               END IF
  770       CONTINUE
         END IF
  780 CONTINUE
  790 CONTINUE
      IF (IPROIT .GE. 25) THEN
         WRITE (LUPRI,'(/A)')
     &      ' ICDTRA: list of distributions included :'
         DO 810 I = 1,NORBT
            WRITE(LUPRI,'(I5,A,(T8,10I7))')
     *        I,' :',(ICDTRA(I,J),J=1,NORBT)
  810    CONTINUE
      END IF
      CALL QEXIT ('OITICD')
      RETURN
      END
C  /* Deck oitcor */
      SUBROUTINE OITCOR(ISYMAB,ISYMCD,IDST,IDEND,ICDXOF,CDEQDC,
     *                  INDCD,IN2CD,INDAB,LUDAXT,H2XCD,H2X,LH2X,NH2XCD)
C
C  Copyright Hans Joergen Aa. Jensen December 1990
C
C  Purpose: read as many half-transformed integrals of symmetry ISYMAB
C           into core as possible.
C
C  Input:
C   CDEQDC: if true then INDCD(i,j) = INDCD(j,i)
C   IDST  : Starting ID
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION INDCD(NORBT,NORBT),IN2CD(NORBT,NORBT),INDAB(NORBT,NORBT)
      DIMENSION H2XCD(N2ORBX),H2X(LH2X,NH2XCD)
      LOGICAL   CDEQDC
C
C  Used from common blocks:
C    INFORB: NORBT, N2ORBT,...
C    INFIND: ISMO(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      CALL QENTER('OITCOR')
C
C     Copy and change sign, negative sign will signify that
C     the distribution is only on disk
C
      DO 2 J = 1,NORBT
      DO 2 I = 1,NORBT
         IN2CD(I,J) = -INDCD(I,J)
   2  CONTINUE
C
C     Find ICDXOF and IDEND
C
      ISYMD  = ISMO(IDST)
  100 ISYMC  = MULD2H(ISYMD,ISYMCD)
      IF (NORB(ISYMC) .EQ. 0 .OR. NORB(ISYMD) .EQ. 0) THEN
         ISYMD = ISYMD + 1
         IF (ISYMD .GT. NSYM) THEN
            ICDXOF = 0
            IDST   = NORBT + 1
            IDEND  = NORBT
            GO TO 9999
         END IF
         IDST  = IORB(ISYMD) + 1
         GO TO 100
      END IF
      ICST   = IORB(ISYMC) + 1
      ICDST  = INDAB(ICST,IDST)
      ICDXOF = ICDST - 1
C
      IDEND  = 0
      LNEED  = 0
      DO 120 ID = IDST, NORBT
         ISYMD  = ISMO(ID)
         ISYMC  = MULD2H(ISYMD,ISYMCD)
         IF (NORB(ISYMC) .EQ. 0) GO TO 120
         LNEED  = LNEED + NORB(ISYMC)
         IF (LNEED .GT. LH2X) GO TO 130
         IDEND = ID
         ICEND = IORB(ISYMC) + NORB(ISYMC)
  120 CONTINUE
  130 CONTINUE
      IF (IDEND .EQ. 0) THEN
         WRITE(LUERR,*) 'OITCOR ERROR: IDEND .eq. 0'
         WRITE(LUERR,*) '  IDST,  IDEND :',IDST,IDEND
         WRITE(LUERR,*) '  LH2X :',LH2X
         CALL QTRACE(LUERR)
         CALL QUIT('OITCOR ERROR: IDEND .eq. 0')
      END IF
      ICDEND = INDAB(ICEND,IDEND)
      NCD    = ICDEND - ICDST + 1
      IF (NCD .GT. LH2X) THEN
         WRITE(LUERR,*) 'OITCOR ERROR: NCD .gt. LH2X'
         WRITE(LUERR,*) '  ICST,  ICEND :',ICST,ICEND
         WRITE(LUERR,*) '  IDST,  IDEND :',IDST,IDEND
         WRITE(LUERR,*) ' ICDST, ICDEND :',ICDST,ICDEND
         WRITE(LUERR,*) '   NCD :',NCD
         WRITE(LUERR,*) '  LH2X :',LH2X
         CALL QTRACE(LUERR)
         CALL QUIT('OITCOR ERROR: NCD .gt. LH2X')
      END IF
C
      IDISXY = 0
      ISYMXY = ISYMAB
C
      DO 220 IY = 1,NORBT
         ISYMY = ISMO(IY)
         ISYMX = MULD2H(ISYMY,ISYMXY)
         IXST  = IORB(ISYMX) + 1
         IXEND = IORB(ISYMX) + NORB(ISYMX)
         DO 210 IX = IXST,IXEND
            IRECXY = -IN2CD(IX,IY)
         IF (IRECXY .LE. 0) GO TO 210
            IDISXY = IDISXY + 1
            IF (ICDST .EQ. 1) THEN
               READ(LUDAXT,REC=IRECXY) (H2X(I,IDISXY),I=1,ICDEND)
            ELSE
               READ(LUDAXT,REC=IRECXY) (H2XCD(I),I=1,ICDEND)
               CALL DCOPY(NCD,H2XCD(ICDST),1,H2X(1,IDISXY),1)
            END IF
            IN2CD(IX,IY) = IDISXY
            IF (CDEQDC) IN2CD(IY,IX) = IDISXY
            IF (IDISXY .EQ. NH2XCD) GO TO 9999
C           .............. EXIT FROM LOOP
  210    CONTINUE
  220 CONTINUE
C
 9999 CALL QEXIT('OITCOR')
      RETURN
      END
C  /* Deck nxth2x */
      SUBROUTINE NXTH2X(IC,ID,H2CD,IDAX,NEEDTP,WRK,KFREE,LFREE,IDIST)
C
C  Copyright Hans Joergen Aa. Jensen November 1990
C
C NOTE: The space allocated in WRK must not be touched outside
C       until all desired distributions have been read.
C
C Purpose:
C    Read next Mulliken two-electron integral distribution (**|cd)
C    where (cd) distribution is needed according to NEEDTP(ITYPCD)
C    (if needtp(itypcd) .gt. 0 all distributions of that type needed;
C     if needtp(itypcd) .lt. 0 at least one distribution needed;
C     if needtp(itypcd) .eq. 0 no distributions of that type needed).
C
C Usage:
C    Set IDIST = 0 before first call of NXTH2X.
C    DO NOT CHANGE IDIST or WRK(KFREE1:KFREE2-1) in calling routine
C    until last distribution has been read (signalled by IDIST .eq. -1)
C    Prototype code:
C     IDIST = 0
C     define NEEDTP(1:6)
C 100 CALL NXTH2X(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
C     IF (IDIST .GT. 0) THEN
C        KW1 = KFREE
C        LW1 = LFREE
C        use (**|cd) distribution in H2CD as desired
C        WRK(KW1:KW1-1+LW1) may be used
C        GO TO 100
C     END IF
C
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION H2CD(*),NEEDTP(6),WRK(LFREE)
C
C Used from common blocks:
C   INFORB : N2ORBX
C   INFPRI : LUERR
C
#include <inforb.h>
#include <infpri.h>
C
      DIMENSION N2DIS(8)
      LOGICAL CDSWAP
      SAVE N2DIS, CDSWAP, IC1, ID1, LUDAX, LSYMX
      SAVE KICDTR, KH2CDP, KNEXT
      DATA KNEXT /-1/
C
      CALL QENTER('NXTH2X')
C
C MAERKE : print level:
C
      IPRINT = 0
C
C     Read untransformed integrals if IDAX = 0
C
      IF (IDAX .EQ. 0) THEN
         IF (IDIST .EQ. 0) CDSWAP = .FALSE.
         IF (CDSWAP) THEN
            IC = ID1
            ID = IC1
            CDSWAP = .FALSE.
         ELSE
            CALL NXTH2M(IC1,ID1,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
            IC = IC1
            ID = ID1
            CDSWAP = (IC1 .NE. ID1)
         END IF
         GO TO 9999
      END IF
C
C     If first read (IDIST .eq. 0) then
C         Open IDAX file and get information about unit and level
C         allocate space for ICDTRA.
C         recalculate ICDTRA
C
      IF (IDIST .EQ. 0) THEN
         CALL OITOPO(IDAX,LUDAX,JTRLVO,LSYMX,LRDAX)
C        CALL OITOPO(IDAXO,LUDAXO,JTRLVO,LSYMXO,LRDAX)
C
         CALL MEMGET('INTE',KICDTR,N2ORBX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KH2CDP,N2ORBT,WRK,KFREE,LFREE)
         KNEXT = KFREE
C        Recalculate index vector for distributions
         CALL OITICD(JTRLVO,WRK(KICDTR),WRK(KFREE),IPRINT)
C        CALL OITICD(JTRLVL,ICDTRA,ITRTYP,IPROIT)
         CALL OITDIS(N2DIS)
      ELSE
C        ... check that work allocation has not been destroyed by
C            calling routine.
         IF (KNEXT.EQ. -1  ) THEN
            WRITE (LUERR,*)
     &         'NXTH2X error, IDIST must be zero in first call'
            WRITE (LUERR,*) 'IDIST =',IDIST
            CALL QTRACE(LUERR)
            CALL QUIT('NXTH2X error, IDIST must be zero in first call')
         END IF
         IF (KFREE.LT.KNEXT) THEN
            WRITE (LUERR,*)
     &         'NXTH2X error, KFREE lower than internal allocation'
            WRITE (LUERR,*) 'KFREE ',KFREE
            WRITE (LUERR,*) 'KICDTR',KICDTR
            WRITE (LUERR,*) 'KNEXT ',KNEXT,
     &         ' ( next avail. address after internal allocation)'
         END IF
         CALL MEMCHK('IDIST .ne. 0 MEMCHK in NXTH2X',WRK,KICDTR)
         IF (KFREE.LT.KNEXT) THEN
            CALL QTRACE(LUERR)
            CALL QUIT('NXTH2X error: KFREE '//
     &                'lower than internal allocation')
         END IF
      END IF
C
      CALL NX2H2X(IC,ID,H2CD,NEEDTP,IDIST, LUDAX,LSYMX,N2DIS,
     *            WRK(KICDTR),WRK(KH2CDP))
C     CALL NX2H2X(IC,ID,H2CD,NEEDTP,IDIST, LUDAX,LSYMX,N2DIS,
C    *            ICDTRA,H2CDPK)
C
C     If no more integrals (IDIST .lt. 0) then release internal space
C
      IF (IDIST .LT. 0) THEN
         CALL MEMREL('Releasing internal space in NXTH2X',WRK,KICDTR,
     &               KICDTR,KFREE,LFREE)
         KNEXT = -1
      END IF
 9999 CALL QEXIT('NXTH2X')
      RETURN
      END
C  /* Deck nx2h2x */
      SUBROUTINE NX2H2X(IC,ID,H2CD,NEEDTP,IDIST,
     *                  LUDAX,LSYMX,N2DIS,ICDTRA,H2CDPK)
C
C  Copyright November 1990 by Hans Joergen Aa. Jensen
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
C Internal:
C       ICDTRA(ICD) .ne. 0 if (**|cd) distribution has been transformed.
C       H2CDPK(*) is used for packed integrals on disk.
C
C ****************************************************************
C
#include <implicit.h>
      DIMENSION H2CD(NORBT,NORBT), H2CDPK(N2ORBT)
      INTEGER   NEEDTP(6), ICDTRA(NORBT,NORBT), N2DIS(8)
C
C
C Used from common blocks:
C   INFORB : N2ORBX,NORBT,NSYM,...
C   INFIND : IOBTYP(*),?
C   INFPRI : IPRSIR
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
C
#include <orbtypdef.h>
C
      SAVE      ICOLD,IDOLD
C
C
C ****************************************************************
C
C     If first read (IDIST .EQ. 0)
C     then initialize ...
C
      IF (IDIST .EQ. 0) THEN
         ICOLD  = 0
         IDOLD  = 1
         IF (IPRSIR .GT. 20) THEN
            WRITE (LUPRI,'(/A//A,8I8)')
     *         ' Test output from NX2H2X for IDIST = 0',
     *         ' N2DIS  array :',(N2DIS(I),I=1,NSYM)
         END IF
         IF (IPRSIR .GT. 50) THEN
            WRITE (LUPRI,'(/A)') ' ICDTRA matrix:'
            DO 10 I = 1,NORBT
               WRITE (LUPRI,'(I5,A,(T8,10I7))') I,' :',
     *            (ICDTRA(I,J),J=1,NORBT)
   10       CONTINUE
         END IF
      END IF
C
C *** Read next distribution which is needed according to NEEDTP(6)
C     into H2CD
C
C  ITYPCD values:  1=i*i :  2=t*i : 3=t*t : 4=a*i : 5=a*t : 6=a*a
C                  0 for not wanted type.
C
C
      ICNEW = ICOLD
      IDNEW = IDOLD
  200 CONTINUE
      ICNEW = ICNEW + 1
      IF (ICNEW .GT. NORBT) THEN
         IDNEW = IDNEW + 1
         ICNEW = 1
      END IF
      IF (IDNEW .GT. NORBT) THEN
C        Last distribution has been read
         IDIST = -1
         GO TO 9999
      END IF
C
      ITYPC  = IOBTYP(ICNEW)
      ITYPD  = IOBTYP(IDNEW)
      ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF ( NEEDTP(ITYPCD) .EQ. 0 ) THEN
         ITYPCD = 0
      ELSE IF (NEEDTP(ITYPCD) .LT. 0) THEN
         ITYPCD = -ITYPCD
      END IF
      IDIST = ICDTRA(ICNEW,IDNEW)
C
      IF (IPRSIR .GT. 50) THEN
         WRITE(LUPRI,*) 'From NX2H2X:'
         WRITE(LUPRI,*) 'ICNEW ,IDNEW        :',ICNEW,IDNEW
         WRITE(LUPRI,*) 'ICDTRA(ICNEW,IDNEW) :',ICDTRA(ICNEW,IDNEW)
         WRITE(LUPRI,*) 'ITYPCD              :',ITYPCD
      END IF
C
      IF (ITYPCD .GT. 0 .AND. IDIST .EQ. 0) THEN
         WRITE (LUERR,*) ' NX2H2X ERROR: needed integral distribution'
         WRITE (LUERR,*) '               has not been calculated'
         WRITE (LUERR,*) 'IC    ,ID     :',ICNEW,IDNEW
         WRITE (LUERR,*) 'IYPC  ,ITYPD  :',COBTYP(ITYPC),COBTYP(ITYPD)
         CALL QTRACE(LUERR)
         CALL QUIT('NX2H2X error: needed integrals not calculated')
      END IF
C
      IF (ITYPCD .EQ. 0) GO TO 200
C
      ISYMC  = ISMO(ICNEW)
      ISYMD  = ISMO(IDNEW)
      ISYMCD = MULD2H(ISYMC,ISYMD)
      ISYMAB = MULD2H(ISYMCD,LSYMX)
      N2DAB  = N2DIS(ISYMAB)
C
      READ(LUDAX,REC=IDIST) (H2CDPK(I),I=1,N2DAB)
C
C     Unpack integrals
C
      CALL OITPAK(H2CD,H2CDPK,ISYMAB,-1)
C
C*******************************************************************
C
C End of subroutine NX2H2X
C
 9999 CONTINUE
      ICOLD  = ICNEW
      IDOLD  = IDNEW
      IC     = ICNEW
      ID     = IDNEW
      RETURN
      END
C  /* Deck oith1x */
#if defined (VAR_OITH1X)
      SUBROUTINE OITH1X(LSYOIT,OITMAT,H1,H1X,IH1SYM)
C
C Copyright 16-Nov-1990 Hans Joergen Aa. Jensen
C Nov 17 1991: secondary-secondary blocks in OITMAT avoided (ov)
C
C One-index transform H1(norbt,norbt) of symmetry IH1SYM
C to H1X(norbt,norbt) using OITMAT(norbt,norbt) of symmetry
C LSYOIT.
C
C Input : H1(a,b) of symmetry IH1SYM
C Output: H1X(a,b) = H1X(a,b)
C                  + H1(a,b) one-index transformed with OITMAT
C                  = H1X(a,b)
C                  + sum(c) [ OITMAT(a,c) H1(c,b)
C                           - H1(a,c) OITMAT(c,b) ]
C
#include <implicit.h>
      DIMENSION OITMAT(NORBT,NORBT), H1(NORBT,NORBT), H1X(NORBT,NORBT)
C
C Used from common blocks:
C  INFORB : NSYM, MULD2H, NORBT, NORB(*), IORB(*)
C
#include <inforb.h>
C
C
C One-index transform both indices and add to previous content.
C
C  (a~ b~) =  SUM(c) [ OITMAT(a,c)*(c b)
C                    - (a c) OITMAT(c,b) ]
C
      DO 100 ICSYM = 1,NSYM
         IBSYM = MULD2H(ICSYM,IH1SYM)
         IASYM = MULD2H(ICSYM,LSYOIT)
         NORBC = NORB(ICSYM)
         NORBB = NORB(IBSYM)
         NORBA = NORB(IASYM)
         IF ((NORBA.NE.0) .AND. (NORBC.NE.0) .AND. (NORBB.NE.0)) THEN
            NOCCC = NOCC(ICSYM)
            NOCCB = NOCC(IBSYM)
            NOCCA = NOCC(IASYM)
            NSSHC = NSSH(ICSYM)
            NSSHB = NSSH(IBSYM)
            NSSHA = NSSH(IASYM)
C
C occ-occ
C
            ICST  = IORB(ICSYM) + 1
            IBST  = IORB(IBSYM) + 1
            IAST  = IORB(IASYM) + 1
            CALL DGEMM('N','N',NOCCA,NOCCB,NORBC,1.D0,
     &                 OITMAT(IAST,ICST),NORBT,
     &                 H1(ICST,IBST),NORBT,1.D0,
     &                 H1X(IAST,IBST),NORBT)
            CALL DGEMM('N','N',NOCCB,NOCCA,NORBC,-1.D0,
     &                 H1(IBST,ICST),NORBT,
     &                 OITMAT(ICST,IAST),NORBT,1.D0,
     &                 H1X(IBST,IAST),NORBT)
C
C occ-sec
C
            ICST  = IORB(ICSYM) + 1
            IBST  = IORB(IBSYM) + 1 + NOCCB
            IAST  = IORB(IASYM) + 1
            CALL DGEMM('N','N',NOCCA,NSSHB,NORBC,1.D0,
     &                 OITMAT(IAST,ICST),NORBT,
     &                 H1(ICST,IBST),NORBT,1.D0,
     &                 H1X(IAST,IBST),NORBT)
            IBST  = IORB(IBSYM) + 1
            IAST  = IORB(IASYM) + 1 + NOCCA
            CALL DGEMM('N','N',NOCCB,NSSHA,NOCCC,-1.D0,
     &                 H1(IBST,ICST),NORBT,
     &                 OITMAT(ICST,IAST),NORBT,1.D0,
     &                 H1X(IBST,IAST),NORBT)
C
C sec-occ
C
            ICST  = IORB(ICSYM) + 1
            IBST  = IORB(IBSYM) + 1
            IAST  = IORB(IASYM) + 1 + NOCCA
            CALL DGEMM('N','N',NSSHA,NOCCB,NOCCC,1.D0,
     &                 OITMAT(IAST,ICST),NORBT,
     &                 H1(ICST,IBST),NORBT,1.D0,
     &                 H1X(IAST,IBST),NORBT)
            IBST  = IORB(IBSYM) + 1 + NOCCB
            IAST  = IORB(IASYM) + 1
            CALL DGEMM('N','N',NSSHB,NOCCA,NORBC,-1.D0,
     &                 H1(IBST,ICST),NORBT,
     &                 OITMAT(ICST,IAST),NORBT,1.D0,
     &                 H1X(IBST,IAST),NORBT)
#if defined (VAR_SECSEC)
C
C sec-sec
C
            ICST  = IORB(ICSYM) + 1
            IBST  = IORB(IBSYM) + 1 + NOCCB
            IAST  = IORB(IASYM) + 1 + NOCCA
            CALL DGEMM('N','N',NSSHA,NSSHB,NOCCC,1.D0,
     &                 OITMAT(IAST,ICST),NORBT,
     &                 H1(ICST,IBST),NORBT,1.D0,
     &                 H1X(IAST,IBST),NORBT)
            CALL DGEMM('N','N',NSSHB,NSSHA,NOCCC,-1.D0,
     &                 H1(IBST,ICST),NORBT,
     &                 OITMAT(ICST,IAST),NORBT,1.D0,
     &                 H1X(IBST,IAST),NORBT)
#endif
         END IF
  100 CONTINUE
C
C     End of OITH1X
C
      RETURN
      END
#endif
