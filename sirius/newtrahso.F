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
#ifdef OLD_REVISION_OLG
c===========================================================================
c950103-hjaaj
cNX2HSO: changed IADRM(x,37,36) to IADRM(x, 1,37)
c        and dimension from IARDM(3,36,36) to IADRM(3,36,37)
c        however, IADRM(1:3,2:36,37) are not defined nor allocated.
c===========================================================================
#endif
C/*deck  trahso*/
      SUBROUTINE TRAHSO(ITRLSO,CMO,WRK,KWORK)
C
C Copyright Mar 1990 Hans Joergen Aa. Jensen
C
#include "implicit.h"
      DIMENSION CMO(*),WRK(KWORK)
C
#include "dummy.h"
#include <maxorb.h>
C
C Used from common blocks:
C   INFTRA : IPRTRA
C   INFTAP : LUAHSO, LUMHSO
C   TRUNIT : LUORDA, LUMINT
C
#include <priunit.h>
#include <inftra.h>
#include <inftap.h>
#include <trhso.h>
#include <trunit.h>
C
C
      LOGICAL NOCOUL,NOEXCH,NOTUVX,MCOPT, FOPEN
C
      CALL QENTER('TRAHSO')
C
      NOCOUL = .FALSE.
      NOEXCH = .TRUE.
      NOTUVX = .TRUE.
      MCOPT  = .FALSE.
      IF (ITRLSO .EQ. 0) THEN
         ITRLVL = 200
      ELSE IF (ITRLSO .EQ. 1) THEN
         ITRLVL = 203
      ELSE
         ITRLVL = 204
      END IF
C
      IF (IPRTRA .GT. 0) DST = SECOND()
C
C
C     Set TRINP and TRUNIT
C     JRDAO = 2 is code for reading spin-orbit integrals
C
      JRDAO = 2
C
C
      IPRFIO = IPRTRA
C
      IF (IPRFIO .GE. 20) THEN
         CALL FASTIO('TRACE=ON')
C        CALL FASTIO('ERRFIX=FIRM')
C        no, ERRFIX=FIRM changes unit number to error code, if error.
      END IF
C
      IF (LUMHSO .LT. 0) CALL DAOPEN(LUMHSO,'MOHSOINT')
C
C     save LUORDA, LUMINT (may point to normal integrals)
      LUORDASV = LUORDA
      LUORDA   = LUAHSO
      LUMINTSV = LUMINT
      LUMINT   = LUMHSO
C
      CALL TRASET(JRDAO,KWORK)
      IF (IPRFIO .GE. 3) CALL FASTIO('STATUS')
C
      IF (IPRTRA .GE. 2) THEN
         WRITE (LUPRI,'(/A)') ' Spin-orbit MO integrals calculated.'
      END IF
C
      CALL TRACT2(ITRLVL,NOTUVX,NOCOUL,NOEXCH,MCOPT,CMO,DUMMY,
     &            WRK,KWORK)
C     CALL TRACT2(ITRLVL,NOTUVX,NOCOUL,NOEXCH,MCOPT,CMO,TUVX,
C    &            WRK,MEMX)
C
C
C     Reset FASTIO options to defaults, if changed:
C
      IF (IPRFIO .GE. 20) THEN
         CALL FASTIO('TRACE=OFF')
C        CALL FASTIO('ERRFIX=HARD')
      END IF
      IF (IPRFIO .GE. 3) CALL FASTIO('STATUS')
C
      IF (IPRTRA .GT. 0) THEN
         TRATIM = SECOND() - DST
         WRITE (LUPRI,'(/A,F10.2)')
     &      ' Total CPU time used in spin-orbit transformation:',TRATIM
      END IF
C     restore LUORDA, LUMINT (may point to normal integrals)
      LUORDA   = LUORDASV
      LUMINT   = LUMINTSV
      CALL QEXIT('TRAHSO')
      RETURN
      END
C/*deck  traset_2*/
      SUBROUTINE TRASET_2( LX2X3,LWORK )
C
C 19-Jun-1990 Hans Joergen Aa. Jensen
C
C Define KEEP and IDATA for spin-orbit transformation
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
      CALL QENTER('TRASET_2')
C
         LX2X3  = 2*600 + 4*600
C        for BUF(600), IBUF(600), and INDX(4,600) for both i*4 and i*8
C        /Bugfix Oct07-hjaaj
C MAERKE: remember to change here when new DR2OUT in her2out.F
C
      KEEP(1) = MSYM
      INTSYM  = KSYMSO
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
      IF (IPRTRA .GT. 1) WRITE(LUPRI,'(/A,I8)')
     *   ' Symmetry of spin-orbit integrals from LUAHSO :',INTSYM
      CALL QEXIT('TRASET_2')
      RETURN
      END
C /*deck  rdhso2*/
      SUBROUTINE RDHSO2(LUHSO2,BUF,IBUF,HSOPQ,MPQOFF)
C
C 20-Mar-1990 hjaaj
C
C Reads two-electron integrals from file to buffer.
C LUHSO2: File to read
C BUF   : Points to start element in integral buffer.
C IBUF  : Points to start element in label buffer.
C HSOPQ : HSOPQ(r,s,pq), pq=(mpqoff+1:mpqoff+npq); dim (nbasr,nbass,npq)
C MPQOFF: Offset of last element returned to calling routine
C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
      DIMENSION HSOPQ(*), BUF(*), IBUF(4,600)
C
C Used from common blocsk:
C   INFORB : IBAS(8)
C   INFIND : IROW(*),ISAO(*)
C   TRHSO  : ILXYZ
C   TRINFO : ISP,ISQ,ISR,ISS,NBQ,NPQ,NBRS
C
#include <maxash.h>
#include <maxorb.h>
#include <nuclei.h>
#include <inforb.h>
#include <infind.h>
#include <trhso.h>
#include <trinfo.h>
C
C     Local variables:
C
      CHARACTER*8 KEY
      LOGICAL DOCOOR
C
#include <ibtfun.h>
C
      IBASP = IBAS(ISP)
      IBASQ = IBAS(ISQ)
      IBASR = IBAS(ISR)
      IBASS = IBAS(ISS)
C
      LBUF = 600
      JBUF = 256
#if defined (VAR_INT64)
      NIBUF = 1
      NBITS = 16
#else
      IF (NBASIS .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         CALL QUIT('NBASIS .gt. 255 not implemented yet here')
         NIBUF = 2
         NBITS = 16
      END IF
#endif
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      KINT = 1
      KIINT = KINT + LBUF
C
      CALL DZERO(HSOPQ,NBRS*NPQ)
C
      CALL REWSPL(LUHSO2)
      KEY = 'AO2SOINT'
      CALL MOLLAB(KEY,LUHSO2,LUPRI)
      DOCOOR = .FALSE.
 150  CONTINUE
#if defined (VAR_SPLITFILES)         
         CALL READSI(LUHSO2,LENINT,BUF(KINT),JBUF)
#else
         CALL READI(LUHSO2,LENINT,BUF(KINT))
#endif
         CALL AOLAB4(BUF(KIINT),LBUF,NIBUF,NBITS,IBUF,LENGTH)
         IF (LENGTH .GT. 0) THEN
            DO 100 I = 1, LENGTH
               JS = IBUF(4,I)
               IF (JS .EQ. 0) THEN
                  ICOOR = IBUF(3,I)
                  DOCOOR = ICOOR.EQ.ILXYZ
               ELSE IF (DOCOOR) THEN
                  JP  = IBUF(1,I)
                  JSP = ISAO(JP)
               IF (JSP .NE. ISP) GO TO 100
                  JQ  = IBUF(2,I)
                  JSQ = ISAO(JQ)
               IF (JSQ .NE. ISQ) GO TO 100
                  JR  = IBUF(3,I)
                  JSR = ISAO(JR)
               IF (JSR .NE. ISR) GO TO 100
                  JSS  = ISAO(JS)
               IF (JSS .NE. ISS) GO TO 100
C
                  MP  = JP - IBASP
                  MQ  = JQ - IBASQ
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
                     HSOPQ((MPQ-1)*NBRS + MRS) = BUF(I)
                  END IF
               ENDIF
 100        CONTINUE
         ELSE IF (LENGTH .LT. 0 ) THEN
            GO TO 300
         END IF
         GO TO 150
C
 300  CONTINUE
      MPQOFF = MPQOFF + NPQ
      RETURN
      END
C /*deck  nxthso*/
      SUBROUTINE NXTHSO(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
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
C    Set IDIST = 0 before first call of NXTHSO.
C    DO NOT CHANGE IDIST or WRK(KFREE1:KFREE2-1) in calling routine
C    until last distribution has been read (signalled by IDIST .eq. -1)
C    Prototype code:
C     IDIST = 0
C     define NEEDTP(1:6)
C 100 CALL NXTHSO(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
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
#include <inforb.h>
#include <infdim.h>
C
      SAVE KFRSAV, KNEXT, KIADRM, KICDTR
      DATA KNEXT /-1/
C
      CALL QENTER('NXTHSO')
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
            WRITE (LUPRI,*)
     &         'NXTHSO error, IDIST must be zero in first call'
            WRITE (LUPRI,*) 'IDIST =',IDIST
            CALL QTRACE(LUPRI)
            CALL QUIT('NXTHSO error, IDIST must be zero in first call')
         END IF
         IF (KFREE.LT.KNEXT) THEN
            WRITE (LUPRI,*)
     &         'NXTHSO error, KFREE lower than buffer allocation'
            WRITE (LUPRI,*) 'KFREE ',KFREE
            WRITE (LUPRI,*) 'KIADRM',KIADRM
            WRITE (LUPRI,*) 'KICDTR',KICDTR
            WRITE (LUPRI,*) 'KNEXT ',KNEXT,
     &         ' ( next avail. address after nxth2m alloc.)'
         END IF
         CALL MEMCHK('IDIST .ne. 0 MEMCHK in NXTHSO',WRK,KFRSAV)
         IF (KFREE.LT.KNEXT) THEN
           CALL QTRACE(LUPRI)
           CALL QUIT('NXTHSO error: KFREE lower than buffer allocation')
         END IF
      END IF
C
      KREL   = KFREE
      CALL MEMGET('REAL',KBUF,2*NORBMA*NORBMA,WRK,KFREE,LFREE)
      CALL NX2HSO(IC,ID,H2CD,NEEDTP,IDIST,WRK(KIADRM),WRK(KICDTR),
     &            WRK(KBUF))
C
C     If finished (IDIST .lt. 0) then release buffer space
C
      IF (IDIST .LT. 0) THEN
         CALL MEMREL('Releasing all buffer space in NXTHSO',WRK,
     &               KFRSAV,KFRSAV,KFREE,LFREE)
         KNEXT = -1
      ELSE
         CALL MEMREL('Releasing BUF(NORBMA,NORBMA) in NXTHSO',WRK,
     &               KFRSAV,KREL,KFREE,LFREE)
      END IF
      CALL QEXIT('NXTHSO')
      RETURN
      END
C /*deck  nx2hso*/
      SUBROUTINE NX2HSO(IC,ID,H2CD,NEEDTP,IDIST,IADRM,ICDTRA,BUF)
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
#include <implicit.h>
      DIMENSION H2CD(NORBT,NORBT),BUF(*)
      INTEGER   IADRM(3,36,37), NEEDTP(6), ICDTRA(NNORBX)
C     ... actual dimension: IARDM(3, (36,36) + 1)
C
#include <iratdef.h>
c#include <dacodes.h>
#include <dalistgs.h>
C
C Used from common blocks:
C   INFORB : NORBX,NORBT,NSYM
C   INFIND : IROW(*)
C   INFPRI : IPRSIR
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <inftap.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
#include <trhso.h>
C
#include <orbtypdef.h>
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
C        MAERKE put this in a new NX0HSO ???
         IAD13  = 0
         LIADRM = 3*1297
         CALL DAREAD(LUMHSO,IADRM,LIADRM,IAD13)
         NLIST  = 4
#if 1
         CALL DARELIST(LUMHSO,IAD13,
     &                 NLIST,ITRLV1,1,INTSY1,1,MSYM,1,LCMO,1)
#else
         CALL GSLIST(LISTGS,NLIST,ITRLV1,1,INTSY1,1,MSYM,1,LCMO,1)
C        NLIST = 10
C        CALL GSLIST(LISTGS,NLIST,ITRLVL,1,INTSYM,1,MSYM,1,LCMO,1,
C    &               MFRO,8,MISH,8,MASH,8,MORB,8,MBAS,8,ITRTYP,MXCORB)
         CALL DAFILE(LUMHSO,DASREA,LISTGS,NLIST,IAD13)
#endif
C MAERKE: hvad med flg. non-standard?
C adressen paa ICDTRA er gemt i IADRM(3,"1297")
C        IADCMO = IADRM(2, 1,37)
         IAD13  = IADRM(3, 1,37)
         CALL DAREAD(LUMHSO,ICDTRA,NNORBX,IAD13)
         IF (IPRSIR .GT. 50) THEN
            WRITE (LUPRI,*) 'Test output from NX2HSO for IDIST = 0'
            WRITE (LUPRI,*) 'ITRLVL from LUMHSO :',ITRLV1
            WRITE (LUPRI,*) 'NSYM   from LUMHSO :',MSYM
            WRITE (LUPRI,*) 'INTSYM from LUMHSO :',INTSY1
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
            WRITE (LUPRI,*) 'Address for ICDTRA',IADRM(3, 1,37)
            WRITE (LUPRI,*) 'Address for CMO   ',IADRM(2, 1,37)
            WRITE (LUPRI,*) 'First free address',IADRM(1, 1,37)
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
         WRITE (LUPRI,*) ' NX2HSO ERROR: needed integral distribution'
         WRITE (LUPRI,*) '               has not been calculated'
         WRITE (LUPRI,*) 'IC    ,ID     :',ICNEW,IDNEW
         WRITE (LUPRI,*) 'ITYPC ,ITYPD  :',COBTYP(ITYPC),COBTYP(ITYPD)
         CALL QTRACE(LUPRI)
         CALL QUIT('NXTHSO error: needed integrals not calculated')
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
               CALL QTRACE(LUPRI)
               CALL QUIT('NX2HSO ERROR, IADCD(ISYMA) .lt. 0')
            END IF
            IF (ITYPCD.EQ.0) THEN
C           ... skip this record
               CALL DASKIP(LUMHSO,LREC,IADCD(ISYMA))
            ELSE
               NINT = NINT + 1
               CALL DAREAD(LUMHSO,BUF(IBUF2),LREC,IADCD(ISYMA))
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
C End of subroutine NX2HSO
C
 9999 CONTINUE
      ICOLD  = ICNEW
      IDOLD  = IDNEW
      IC     = ICNEW
      ID     = IDNEW
      RETURN
      END
