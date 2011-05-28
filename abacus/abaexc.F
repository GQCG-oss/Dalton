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
C  /* Deck excinp */
      SUBROUTINE EXCINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (NTABLE = 16)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <cbiexc.h>
C
#include <gnrinf.h>
#include <abainf.h>
#include <anrinf.h>
#include <dorps.h>
#include <nuclei.h>
#include <symmet.h>
      DATA TABLE /'.SKIP  ', '.INTPRI','.MAX IT','.THRESH',
     *            '.MAXRED', '.MAXPHP','.TRIPLE','.XXXXXX',
     *            '.OPTORB', '.NEXCIT','.DIPSTR','.PRINT ',
     *            '.ROTVEL', '.FNAC  ','.STOP  ','.SUMRUL'/
C
      NEWDEF = (WORD .EQ. '*EXCITA')
      ICHANG = 0
      IF (NEWDEF) THEN
         NSYM = MAXREP + 1 ! NSYM in inforb.h is not always defined when EXCINP is called
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/4A/)') ' Keyword "',WORD,
     *            '" not recognized for ',WORD1
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword for '//WORD1)
    1          CONTINUE
                  SKIP = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPR1IN
                  ICHANG = ICHANG + 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) MAXITE
                  ICHANG = ICHANG + 1
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) THREXC
                  ICHANG = ICHANG + 1
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) MXRM
                  ICHANG = ICHANG + 1
               GO TO 100
    6          CONTINUE
                  READ (LUCMD,*) MXPHP
                  ICHANG = ICHANG + 1
               GO TO 100
    7          CONTINUE
                  EXCTRP = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    8          CONTINUE
               GO TO 100
    9          CONTINUE
                  OOTV   = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   10          CONTINUE
                  READ (LUCMD,*) (NEXCIT(I),I=1,NSYM)
                  ICHANG = ICHANG + 1
               GO TO 100
   11          CONTINUE
                  DIPSTR = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   12          CONTINUE
                  READ (LUCMD,*) IPREXC
                  ICHANG = ICHANG + 1
               GO TO 100
   13          CONTINUE
                  ROTVEL = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   14          CONTINUE
                  FNAC   = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   15          CONTINUE
                  CUT    = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   16          CONTINUE
                  SUMRUL = .TRUE.
CSPAS:23/5-11: second and third moment sum rules
C                 DIPSTR = .TRUE.
CKeinSPASmehr
                  ICHANG = ICHANG + 1
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in EXCINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in EXCINP.')
            END IF
      END IF
  300 CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         THREXC = THREXC*THR_REDFAC
      END IF
C
      NTOT = ISUM(NSYM,NEXCIT,1)
C
      IF (NTOT.EQ.0) THEN
         NEXCIT(:) = 1
      END IF
C
CSPAS : let's try whether we cannot calculate all excitation energies
C
C     maximum number of excitation energies
C
      MXNEXI = 1
      DO ISYM = 1,8
         MXNEXI = MAX(MXNEXI,NEXCIT(ISYM))
      END DO
CKeinSPASmehr
C
      CALL HEADER('Changes of defaults for .EXCITA:',0)
      IF (EXCTRP) WRITE (LUPRI,'(A)') ' Triplet excitation energies'
      WRITE (LUPRI,'(A,8I5)')
     &   ' Number of excitation energies:',(NEXCIT(I),I=1,8)
C
      IF (ICHANG .GT. 0) THEN
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' EXCITA skipped in this run.'
         ELSE
            WRITE (LUPRI,'(A,I5)')
     &         ' Print level          :',IPREXC
            WRITE (LUPRI,'(A,I5)')
     &         ' Integral print level :',IPR1IN
            WRITE (LUPRI,'(A,1P,E9.2)')
     &         ' Threshold            :',THREXC
            WRITE(LUPRI,'(A,I5)')
     &         ' Maximum iterations   :',MAXITE
            IF (ROTVEL)
     &         WRITE (LUPRI,'(A)') ' Conventional rotational strength'
            IF (DIPSTR) WRITE (LUPRI,'(A)') ' Dipole strength'
            IF (FNAC) WRITE (LUPRI,'(A)')
     &         ' First order nonadiabatic coupling elements'
            IF (CUT) WRITE (LUPRI,'(//A)')
     &         ' Program is stopped after EXCITA.'
         END IF
      END IF
      RETURN
      END
C  /* Deck excini */
      SUBROUTINE EXCINI
C
C     Initialize /cbiexc/
C
#include <implicit.h>
#include <mxcent.h>
#include <cbiexc.h>
#include <abainf.h>
C
      IPR1IN = IPRDEF
      IPREXC = IPRDEF
      SKIP   = .FALSE.
      CUT    = .FALSE.
      OOTV   = .FALSE.
      DIPSTR = .FALSE.
      ROTVEL = .FALSE.
      FNAC   = .FALSE.
      EXCTRP = .FALSE.
      ROTSTR = .FALSE.
      CALL IZERO(NEXCIT,8)
      THREXC = 1.D-04
      MAXITE = 60
      MXRM   = 400
      MXPHP  = 0
      NABAPP = 0
      MXNEXI = 1
C
      RETURN
      END
C  /* Deck excita */
      SUBROUTINE EXCITA(TRLEN,TRVEL,TQLEN,TQVEL,
CSPAS:23/5-11: second and third moment sum rules
C    &                  TRLON,TRMAG,BSRLON,EXENG,FONAC,
     &                  TRLON,TRMAG,BSRLON,TTLEN,EXENG,FONAC,
CKeinSPASmehr
     &                  FONA2,SECMAT,WORK,LWORK,PASS)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <inforb.h>
#include <cbiexc.h>
#include <infrsp.h>
#include <nuclei.h>
CSPAS : 6/11-2009 AO-SOPPA included
#include <soppinf.h>
CKeinSPASmehr
      LOGICAL PASS
      DIMENSION WORK(LWORK)
      DIMENSION TRLEN(3,NSYM,MXNEXI,*), TRVEL(3,NSYM,MXNEXI,*)
CSPAS:23/5-11: second and third moment sum rules
C     DIMENSION TQLEN(3,3,NSYM,MXNEXI), TQVEL(3,3,NSYM,MXNEXI)
      DIMENSION TQLEN(3,3,NSYM,MXNEXI,*), TQVEL(3,3,NSYM,MXNEXI,*)
      DIMENSION TTLEN(10,NSYM,MXNEXI,*)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI,*), TRMAG(3,NSYM,MXNEXI,*)
      DIMENSION BSRLON(3,NSYM,MXNEXI,*),EXENG(NSYM,MXNEXI,*)
      DIMENSION FONAC(3*NATOMS,NSYM,MXNEXI),FONA2(3*NATOMS,NSYM,MXNEXI)
      DIMENSION SECMAT(3,MXNEXI,NSYM)
C
      IF (SKIP) RETURN
      CALL QENTER('EXCITA')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPREXC .GT. 0) WRITE (LUPRI,'(A,/)')
     *    '  <<<<<<<<<< Output from EXCITA >>>>>>>>>> '
C
      IPRRSP = IPREXC
C
      IF (AORPA .OR. DCRPA .OR. AOSOP .OR. AOSOC) THEN
CSPAS:23/5-11: second and third moment sum rules
C        CALL SO_EXCIT1(TRLEN,TRVEL,TRLON,TRMAG,BSRLON,SECMAT,EXENG,
C    &                  WORK,LWORK)
         CALL SO_EXCIT1(TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,TRLON,TRMAG,
     &                  BSRLON,SECMAT,EXENG,WORK,LWORK)
CKeinSPASmehr
      ELSE
CSPAS:23/5-11: second and third moment sum rules
C        CALL EXCIT1(TRLEN,TRVEL,TQLEN,TQVEL,
         CALL EXCIT1(TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,
CKeinSPASmehr
     &               TRLON,TRMAG,BSRLON,EXENG,FONAC,FONA2,
     &               WORK,LWORK)
      END IF
C
      CALL TIMER ('EXCITA',TIMEIN,TIMOUT)
      PASS   = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/A)')
     &          ' Program stopped after EXCITA as requested.'
         CALL QUIT(' ***** End of ABACUS (in EXCITA) *****')
      END IF
      CALL QEXIT('EXCITA')
      RETURN
      END
C  /* Deck excit1 */
CSPAS:23/5-11: second and third moment sum rules
C      SUBROUTINE EXCIT1(TRLEN,TRVEL,TQLEN,TQVEL,
      SUBROUTINE EXCIT1(TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,
CKeinSPASmehr
     &                  TRLON,TRMAG,BSRLON,EXENG,FONAC,
     &                  FONA2,WORK,LWORK)
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <cbiexc.h>
#include <inforb.h>
#include <nuclei.h>
#include <inflin.h>
#include <infdim.h>
#include <inftap.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <maxmom.h>
#include <maxaqn.h>
#include <symmet.h>
#include <abainf.h>
#include <gnrinf.h>
#include <infsop.h>
C
      DIMENSION WORK(LWORK)
      DIMENSION TRLEN(3,NSYM,MXNEXI), TRVEL(3,NSYM,MXNEXI)
      DIMENSION TQLEN(3,3,NSYM,MXNEXI), TQVEL(3,3,NSYM,MXNEXI)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TTLEN(10,NSYM,MXNEXI)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI), TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI),EXENG(NSYM,MXNEXI)
      DIMENSION FONAC(3*NATOMS,NSYM,MXNEXI),FONA2(3*NATOMS,NSYM,MXNEXI)
      CHARACTER*8 LABEL, LABINT(3*MXCOOR)
      LOGICAL   LOCDIPLEN, LOCDIPVEL, LOCANGMOM, LOCLONDON, LOCSECMOM,
CSPAS:23/5-11: second and third moment sum rules
C    &          LOCROTSTR, LOCHBDO, LOCHDO
     &          LOCTHIRDM, LOCROTSTR, LOCHBDO, LOCHDO
CKeinSPASmehr
      LOGICAL   FOUND
C
      CALL QENTER('EXCIT1')
C
C     If OECD or ECD true then
C     set flags in cbiexc.h true to calculate needed integrals
C
      IF (OECD) THEN
         DIPSTR = .TRUE.
         ROTVEL = .TRUE.
         ROTSTR = .TRUE.
      ELSE IF (ECD) THEN
         ROTSTR = .TRUE.
      END IF
C
C     Get reference state
C     ===================
C
C     1. Work allocations:
C
      IF (ABASOP) THEN
         LUDV   = NORBT * NORBT
         LPVX   = LPVMAT
         IF (EXCTRP) LPVX = 2*LPVMAT
      ELSE
         LUDV   = N2ASHX
         LPVX   = 0
      ENDIF
      KFREE  = 1
      LFREE  = LWORK
C
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDV  ,LUDV  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVX  ,LPVX  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXINDX,LCINDX,WORK,KFREE,LFREE)
C
      KWORK1 = KFREE
      LWORK1 = LFREE
C
      CALL RD_SIRIFC('CMO',FOUND,WORK(KCMO),WORK(KFREE),LFREE)
      IF (.NOT.FOUND) CALL QUIT('EXCIT1 error: CMO not found on SIRIFC')
      IF (NASHT .GT. 0) THEN
         CALL RD_SIRIFC('DV',FOUND,WORK(KWORK1),WORK(KWORK1),LWORK1)
         IF (.NOT.FOUND)
     &      CALL QUIT('EXCIT1 error: DV not found on SIRIFC')
         CALL DSPTSI(NASHT,WORK(KWORK1),WORK(KUDV))
      END IF
C
      LFREE = LWORK1
C
      CALL LNRVAR(LSYMRF,IPREXC,WORK(KWORK1),LFREE)
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LFREE,0)
C
C     SOPPA :
C
      IF (ABASOP) THEN
C
C        Initialize XINDX
C
         CALL DZERO(WORK(KXINDX),LCINDX)
C
C        Find address array's for SOPPA calculation
C
         CALL SET2SOPPA(WORK(KXINDX+KABSAD-1),WORK(KXINDX+KABTAD-1),
     *                  WORK(KXINDX+KIJSAD-1),WORK(KXINDX+KIJTAD-1),
     *                  WORK(KXINDX+KIJ1AD-1),WORK(KXINDX+KIJ2AD-1),
     *                  WORK(KXINDX+KIJ3AD-1),WORK(KXINDX+KIADR1-1))
C
C
         REWIND (LUSIFC)
         IF (CCPPA) THEN
            CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
         ELSE
            CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
         ENDIF
C
C        reads the MP2 or CCSD correlation coefficients into PV
C
         CALL READT (LUSIFC,LPVMAT,WORK(KPVX))
C
         IF (IPREXC.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' EXCIT1 : CCSD correlation ',
     &                           'coefficients'
            ELSE
               WRITE(LUPRI,'(/A,A)')' EXCIT1 :',
     &                              ' MP2 correlation coefficients'
            ENDIF
            CALL OUTPUT(WORK(KPVX),1,LPVMAT,1,1,LPVMAT,1,1,LUPRI)
         END IF
C
         IF (EXCTRP) THEN
            CALL TRPKAP(WORK(KPVX),WORK(KPVX+LPVMAT),
     &                  WORK(KXINDX+KIADR1-1),WORK(KWORK1))
C
            IF (IPREXC.GT.10) THEN
               IF (CCPPA) THEN
                  WRITE(LUPRI,'(/A)')' RSPMC : CCSD triplet ',
     &                              'correlation coefficients'
               ELSE
                  WRITE(LUPRI,'(/A)')' RSPMC : MP2 triplet correlation',
     &                               ' coefficients'
               ENDIF
               CALL OUTPUT(WORK(KPVX),1,2*LPVMAT,1,1,2*LPVMAT,1,1,LUPRI)
            END IF
         END IF
C
C        reads the MP2 or CCSD second order one particle density matrix 
C
         CALL READT (LUSIFC,NORBT*NORBT,WORK(KUDV))
C
C        UDV contains the MP2 one-density. Remove the diagonal
C        contribution from the zeroth order. (Added in MP2FAC)
C
         IF (IPREXC.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' RSPMC : CCSD density'
            ELSE
               WRITE(LUPRI,'(/A)')' RSPMC : MP2 density'
            END IF
            CALL OUTPUT(WORK(KUDV),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,
     &                  LUPRI)
         END IF
C
         CALL SOPUDV(WORK(KUDV))
C
      END IF
C
C      CALL MCORL(WORK(KWORK1),LWORK1,PASS)
C
C     Construct property-integrals and write to LUPROP
C     ================================================
C
      LOCDIPLEN = .FALSE.
      LOCDIPVEL = .FALSE.
      LOCANGMOM = .FALSE.
      LOCLONDON = .FALSE.
      LOCSECMOM = .FALSE.
      LOCROTSTR = .FALSE.
CSPAS:23/5-11: second and third moment sum rules
      LOCTHIRDM = .FALSE.
CKeinSPASmehr
      LOCHBDO   = .FALSE.
      LOCHDO    = .FALSE.

      IF (DIPSTR) THEN         ! oscillator strength, len+vel
         LOCDIPLEN = .TRUE.
         LOCDIPVEL = .TRUE.
      END IF
      IF (ROTSTR) THEN         ! rot. str., len+London
         LOCDIPLEN = .TRUE.
         LOCLONDON = .TRUE.
         LOCANGMOM = .TRUE.
         IF (NODIFC) LOCHBDO = .TRUE.
      END IF
      IF (ROTVEL) THEN         ! rot. str., len+vel
         LOCDIPLEN = .TRUE.
         LOCDIPVEL = .TRUE.
         LOCANGMOM = .TRUE.
      END IF
      IF (OECD) THEN           ! rot. str. tensors, len+vel
         LOCDIPLEN = .TRUE.
         LOCDIPVEL = .TRUE.
         LOCANGMOM = .TRUE.
         LOCSECMOM = .TRUE.
         LOCROTSTR = .TRUE.
      END IF
CSPAS:23/5-11: second and third moment sum rules
      IF (SUMRUL) THEN           ! sum rules, dip, quad, oct.
         LOCDIPLEN = .TRUE.
         LOCDIPVEL = .TRUE.
         LOCSECMOM = .TRUE.
         LOCTHIRDM = .TRUE.
      END IF
CKeinSPASmehr
      IF (FNAC .AND. NODIFC) THEN
         LOCHDO    = .TRUE.
      END IF
C
C     2. Work allocations:
C
      KIDSYM = KWORK1
      KIDADR = KIDSYM + 9*MXCENT
      KWORK2 = KIDADR + 9*MXCENT
      LWORK2 = LWORK  - KWORK2
C
      NLBTOT = 0
C
      IF (LOCDIPLEN) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'DIPLEN ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,
     &               DUMMY,IPR1IN)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      END IF
      IF (LOCDIPVEL) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'DIPVEL ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &               IPR1IN)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      END IF
      IF (LOCANGMOM) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'ANGMOM ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &               IPR1IN)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,
     &               LABSYM)
      END IF
      IF (LOCLONDON) THEN
         CALL LABCOP(1,NLBTOT,'XLONMAG ',ISYMAX(1,2),LABAPP,LABSYM)
         CALL LABCOP(1,NLBTOT,'YLONMAG ',ISYMAX(2,2),LABAPP,LABSYM)
         CALL LABCOP(1,NLBTOT,'ZLONMAG ',ISYMAX(3,2),LABAPP,LABSYM)
      END IF
      IF (LOCSECMOM) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'SECMOM ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &               IPR1IN)
         NLAB = 6
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,
     &               LABSYM)
      END IF
CSPAS:23/5-11: second and third moment sum rules
      IF (LOCTHIRDM) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'THRMOM ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &               IPR1IN)
         NLAB = 10
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,
     &               LABSYM)
      END IF
CKeinSPASmehr
      IF (LOCROTSTR) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'ROTSTR ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &               IPR1IN)
         NLAB = 6
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,
     &                  LABSYM)
      END IF
      IF (LOCHBDO) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'HBDO   ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,
     &               IPR1IN)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,
     &               LABSYM)
      END IF
      IF (LOCHDO) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'HDO    ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),IDUMMY,
     &               .TRUE.,NPATOM,.TRUE.,DUMMY,.FALSE.,DUMMY,IPR1IN)
         NLAB = 3 * NUCDEP
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      ENDIF
C
      CALL DZERO(TRLEN,3*NSYM*MXNEXI)
      CALL DZERO(TRVEL,3*NSYM*MXNEXI)
      CALL DZERO(TQLEN,3*3*NSYM*MXNEXI)
      CALL DZERO(TQVEL,3*3*NS*MXNEXI)
CSPAS:23/5-11: second and third moment sum rules
      CALL DZERO(TTLEN,10*NSYM*MXNEXI)
CKeinSPASmehr
      CALL DZERO(TRLON,3*NSYM*MXNEXI)
      CALL DZERO(TRMAG,3*NSYM*MXNEXI)
      CALL DZERO(BSRLON,3*NSYM*MXNEXI)
      CALL DZERO(FONAC,3*NUCDEP*NSYM*MXNEXI)
      CALL DZERO(FONA2,3*NUCDEP*NSYM*MXNEXI)
C
C     Loop over all symmetries
C     ========================
C
      LUSOVE = -9002
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
CSPAS : 26/02-2008 Keld Bak's VIB-G 
      NBR = 0
C
CKeinSPASmehr
      DO 100 ISYM = 1, NSYM
C
         CALL ABAVAR(ISYM,EXCTRP,IPREXC,WORK(KWORK2),LWORK2)
C
         NEXVAL = MIN(NEXCIT(ISYM),NVARPT)
C
         IF (NEXVAL.GT.0) THEN
C
C           3. Work allocations:
C
            KEXVAL = KWORK1
            KWRK1  = KEXVAL + NEXVAL
            LWRK1  = LWORK - KWRK1
            KGD1   = KWRK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KGD1
            KSLV   = KGD1 + 2*NVARPT
C           KLAST  = KSLV + 2*NVARPT
            KWRK2  = KSLV + 2*NVARPT
            KLAST  = KWRK2 + 2*NVARPT
            IF (KLAST.GT.LWORK) CALL STOPIT('EXCITA',' ',KLAST,LWORK)
            KWRK = KLAST
            LWRK = LWORK - KLAST + 1
C
C           Find excitation energies and excitation vectors
C           Excitation vectors are written to file
C           ===============================================
C
            CALL EXCIT2(WORK(KEXVAL),NEXVAL,ISYM,LUSOVE,
     &                  WORK(KWRK1),LWRK1)
C
            DO 70 IEXVAL = 1,NEXVAL
               EXENG(ISYM,IEXVAL) = WORK(KEXVAL-1+IEXVAL)
   70       CONTINUE
C
            IF (FNAC .AND. (.NOT. NODIFC)) THEN
C
              DO 72 IATOM = 1,NUCIND
C
                DO 71 ICOOR = 1,3
C
                  ISCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM-1,1)
C
                  IF (ISCOOR.GT.0) THEN
C
                    NBR = NBR + 1
C
                    CALL RDS2X(NBR,WORK(KUDV),WORK(KXINDX),
     &                         WORK(KWRK2),LWRK)
C
                    IF (IPREXC.GT.3) THEN
                       WRITE (LUPRI,'(2A)') 'GP Vector, label: ','d/dR'
                       CALL OUTPUT(WORK(KWRK2),1,2*NVARPT,1,1,2*NVARPT,
     &                             1,1,LUPRI)
                    END IF
C
                    REWIND LUSOVE
C
                    DO 73 IEXVAL = 1,NEXVAL
C
                      CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
                      EXMOM = DDOT(2*NVARPT,WORK(KSLV),1,WORK(KWRK2),1)
C
                      FONAC(NBR,ISYM,IEXVAL) = EXMOM
C
   73               CONTINUE
C
                  END IF
C
   71           CONTINUE
C
   72         CONTINUE
C
            ENDIF
C
C           Loop over property labels
C           =========================
C
            DO 90 IPRLBL = 1, NLBTOT
              LABEL = LABAPP(IPRLBL)
              KSYM  = LABSYM(IPRLBL)
C
C             ISYM = Reference state symmetry times excited state symmetry.
C
              IF (KSYM.EQ.ISYM) THEN
                KSYMOP = ISYM
                CALL GETGPV(LABEL,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &               WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
                IF (IPREXC.GT.3) THEN
                  WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL
                   CALL OUTPUT(WORK(KGD1),1,2*NVARPT,1,1,2*NVARPT,
     &                         1,1,LUPRI)
                ENDIF
C
C               Form transition moments
C               =======================
C
                Rewind LUSOVE
                DO 80 IEXVAL = 1,NEXVAL
                  CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
                  EXMOM = DDOT(2*NVARPT,WORK(KSLV),1,WORK(KGD1),1)
C
                  IF (IPREXC.GT.3) THEN
                     WRITE (LUPRI,'(A,I4)') 'Solution Vector no. ',
     &                      IEXVAL
                     CALL OUTPUT(WORK(KSLV),1,2*NVARPT,1,1,2*NVARPT,
     &                           1,1,LUPRI)
                     WRITE (LUPRI,'(A,F15.8)')
     &               'Excitation energy = ',WORK(KEXVAL-1+IEXVAL)
C-tbp                WRITE (LUPRI,'(3A,F15.8)')
C-tbp&               'Transition moment for ',LABEL,' = ',EXMOM
                     WRITE (LUPRI,'(A,I4,A,I2,3A,F15.8)')
     &               '<',IEXVAL,' (sym.',ISYM,')|',LABEL,'|0> = ',EXMOM
                  ENDIF
C
C                 Write transition momets into relevant arrays.
C                 ==============================================
C                 We have calculated <f|Operator|0> but want to
C                 keep them as <O|Operator|f>. Therefore, for anti-
C                 hermitian operators we multiply with a minus-sign.
C hjaaj 080830 TODO CHECK: sounds strange, p_x = -i d/dx, does that
C                 not give a double sign change ????????
C                 cf. the tbp comments later, changing signs of l_x etc back!!!
C
                  IF (LABEL(2:7).EQ.'DIPLEN') THEN
                    IF (LABEL(1:1).EQ.'X') TRLEN(1,ISYM,IEXVAL) =  EXMOM
                    IF (LABEL(1:1).EQ.'Y') TRLEN(2,ISYM,IEXVAL) =  EXMOM
                    IF (LABEL(1:1).EQ.'Z') TRLEN(3,ISYM,IEXVAL) =  EXMOM
                  ELSE IF (LABEL(2:7).EQ.'DIPVEL') THEN
                    IF (LABEL(1:1).EQ.'X') TRVEL(1,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Y') TRVEL(2,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Z') TRVEL(3,ISYM,IEXVAL) = -EXMOM
                  ELSE IF (LABEL(3:8).EQ.'SECMOM') THEN
                    IF (LABEL(1:2) .EQ. 'XX') THEN
                       TQLEN(1,1,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'XY') THEN
                       TQLEN(1,2,ISYM,IEXVAL) = EXMOM
                       TQLEN(2,1,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'XZ') THEN
                       TQLEN(1,3,ISYM,IEXVAL) = EXMOM
                       TQLEN(3,1,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'YY') THEN
                       TQLEN(2,2,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'YZ') THEN
                       TQLEN(2,3,ISYM,IEXVAL) = EXMOM
                       TQLEN(3,2,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'ZZ') THEN
                       TQLEN(3,3,ISYM,IEXVAL) = EXMOM
                    ELSE
                       CALL QUIT('Logical error (SECMOM) in EXCIT1')
                    END IF
CSPAS:23/5-11: second and third moment sum rules
                  ELSE IF (LABEL(5:8).EQ.'3MOM') THEN
                    IF (LABEL(1:3) .EQ. 'XXX') THEN
                       TTLEN(1,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:3) .EQ. 'YYY') THEN
                       TTLEN(2,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:3) .EQ. 'ZZZ') THEN
                       TTLEN(3,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:3) .EQ. 'XXY') THEN
                       TTLEN(4,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:3) .EQ. 'XXZ') THEN
                       TTLEN(5,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:3) .EQ. 'XYY') THEN
                       TTLEN(6,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:3) .EQ. 'XYZ') THEN
                       TTLEN(7,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:3) .EQ. 'XZZ') THEN
                       TTLEN(8,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:3) .EQ. 'YYZ') THEN
                       TTLEN(9,ISYM,IEXVAL) = EXMOM
                    ELSE IF (LABEL(1:3) .EQ. 'YZZ') THEN
                       TTLEN(10,ISYM,IEXVAL) = EXMOM
                    ELSE
                       CALL QUIT('Logical error (THRMOM) in EXCIT1')
                    END IF
CKeinSPASmehr
                  ELSE IF (LABEL(3:8).EQ.'ROTSTR') THEN
                    IF (LABEL(1:2) .EQ. 'XX') THEN
                       TQVEL(1,1,ISYM,IEXVAL) = -EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'XY') THEN
                       TQVEL(1,2,ISYM,IEXVAL) = -EXMOM
                       TQVEL(2,1,ISYM,IEXVAL) = -EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'XZ') THEN
                       TQVEL(1,3,ISYM,IEXVAL) = -EXMOM
                       TQVEL(3,1,ISYM,IEXVAL) = -EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'YY') THEN
                       TQVEL(2,2,ISYM,IEXVAL) = -EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'YZ') THEN
                       TQVEL(2,3,ISYM,IEXVAL) = -EXMOM
                       TQVEL(3,2,ISYM,IEXVAL) = -EXMOM
                    ELSE IF (LABEL(1:2) .EQ. 'ZZ') THEN
                       TQVEL(3,3,ISYM,IEXVAL) = -EXMOM
                    ELSE
                       CALL QUIT('Logical error (ROTSTR) in EXCIT1')
                    END IF
                  ELSE IF (LABEL(2:7).EQ.'ANGMOM') THEN
                    IF (LABEL(1:1).EQ.'X') TRMAG(1,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Y') TRMAG(2,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Z') TRMAG(3,ISYM,IEXVAL) = -EXMOM
                  ELSE IF (LABEL(2:7).EQ.'LONMAG') THEN
                    IF (LABEL(1:1).EQ.'X') TRLON(1,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Y') TRLON(2,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Z') TRLON(3,ISYM,IEXVAL) = -EXMOM
                  ELSE IF (LABEL(2:6).EQ.'HBDO ') THEN
                    IF (LABEL(7:7).EQ.'X') BSRLON(1,ISYM,IEXVAL)=
     &                                     -EXMOM * EXENG(ISYM,IEXVAL)
                    IF (LABEL(7:7).EQ.'Y') BSRLON(2,ISYM,IEXVAL)=
     &                                     -EXMOM * EXENG(ISYM,IEXVAL)
                    IF (LABEL(7:7).EQ.'Z') BSRLON(3,ISYM,IEXVAL)=
     &                                     -EXMOM * EXENG(ISYM,IEXVAL)
                  ELSE IF (LABEL(1:4).EQ.'HDO ') THEN
                    DO 57 INR = 1, 3*NUCDEP
                       READ(LABEL(4:6),'(I3)') NR
                       IF (INR.EQ.NR) FONA2(NR,ISYM,IEXVAL) = EXMOM
   57               CONTINUE
C
                    CALL RDS2(LABEL,WORK(KUDV),WORK(KXINDX),
     &                        WORK(KWRK2),LWRK)
C
                    EXMOM = DDOT(2*NVARPT,WORK(KSLV),1,WORK(KWRK2),1)
C
                    DO 58 INR = 1, 3*NUCDEP
                       READ(LABEL(4:6),'(I3)') NR
                       IF (INR.EQ.NR) FONAC(NR,ISYM,IEXVAL) = EXMOM
   58               CONTINUE
C
                  ENDIF
   80           CONTINUE
              ENDIF
   90       CONTINUE
         END IF
  100 CONTINUE
      CALL GPCLOSE(LUSOVE,'DELETE')
C
      CALL QEXIT('EXCIT1')
C
      RETURN
      END
C  /* Deck excout */
      SUBROUTINE EXCOUT(TRLEN,TRVEL,TQLEN,TQVEL,TRMAG,TRLON,BSRLON,
CSPAS:23/5-11: second and third moment sum rules
C    &                  EXENG,FONAC,FONA2,RMLEN,RQLEN,RLEN,
     &                  TTLEN,EXENG,FONAC,FONA2,RMLEN,RQLEN,RLEN,
CKeinSPASmehr
     &                  RMVEL,RQVEL,RVEL,SLEN,SVEL,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <codata.h>
#include <gnrinf.h>
#include <cbiexc.h>
#include <inforb.h>
#include <nuclei.h>
#include <abainf.h>
#include <orgcom.h>
      CHARACTER COORDI*21, NMDTXT*9, LAPRPC*8
      DIMENSION TRLEN(3,NSYM,MXNEXI), TRVEL(3,NSYM,MXNEXI)
      DIMENSION TQLEN(3,3,NSYM,MXNEXI), TQVEL(3,3,NSYM,MXNEXI)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TTLEN(10,NSYM,MXNEXI)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI), TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI),EXENG(NSYM,MXNEXI)
      DIMENSION FONAC(3*NATOMS,NSYM,MXNEXI),FONA2(3*NATOMS,NSYM,MXNEXI)
      DIMENSION RMLEN(3,3,NSYM,MXNEXI), RQLEN(3,3,NSYM,MXNEXI),
     &          RLEN(3,3,NSYM,MXNEXI)
      DIMENSION RMVEL(3,3,NSYM,MXNEXI), RQVEL(3,3,NSYM,MXNEXI),
     &          RVEL(3,3,NSYM,MXNEXI)
      DIMENSION SLEN(NSYM,MXNEXI), SVEL(NSYM,MXNEXI)
      DIMENSION WORK(LWORK)
      DIMENSION DSSUML(-6:2,1:4),DLSUML(-6:2,1:4),DISUML(-6:2,1:4)
      DIMENSION DSSUMM(-6:2,1:4),DLSUMM(-6:2,1:4),DISUMM(-6:2,1:4)
      DIMENSION DSSUMV(-6:2,1:4),DLSUMV(-6:2,1:4),DISUMV(-6:2,1:4)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION QSSUM(-6:2,1:4),DTSSUM(-6:2,1:4)
CKeinSPASmehr
      DIMENSION XLEVICI(3,3,3), AVEC(3), XNEW(3)
      LOGICAL   LCOECD, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      PARAMETER (TOLRTT = 1.0D-12)
      PARAMETER (DTWO =  2.0D0, DTHREE = 3.0D0,
     &           ZERO = 0.0D0, ESUDIP = 1D4*DEBYE**2,
     &           ESUECD =  ECHARGE*XTANG*CCM*1D36*ECHARGE*HBAR/EMASS)
C
C-tbp: setup for OECD.
C
      LOCWRN = 0
      LCOECD = .FALSE.

      CALL DZERO(XLEVICI,27)
      XLEVICI(3,2,1) = -1.0D0
      XLEVICI(2,3,1) =  1.0D0
      XLEVICI(3,1,2) =  1.0D0
      XLEVICI(1,3,2) = -1.0D0
      XLEVICI(2,1,3) = -1.0D0
      XLEVICI(1,2,3) =  1.0D0
      IF (OECD) THEN
         CALL DZERO(RMLEN,3*3*NSYM*MXNEXI)
         CALL DZERO(RQLEN,3*3*NSYM*MXNEXI)
         CALL DZERO(RLEN,3*3*NSYM*MXNEXI)
         CALL DZERO(RMVEL,3*3*NSYM*MXNEXI)
         CALL DZERO(RQVEL,3*3*NSYM*MXNEXI)
         CALL DZERO(RVEL,3*3*NSYM*MXNEXI)
         CALL DZERO(SLEN,NSYM*MXNEXI)
         CALL DZERO(RLEN,NSYM*MXNEXI)
         XNRM = DSQRT(DDOT(3,GAGORG,1,GAGORG,1))
         IF (XNRM .GT. 0.0D0) THEN
            WRITE(LUPRI,'(//,10X,A,/,10X,A,/,10X,A,//)')
     &      '*** NOTICE: OECD only at mathematical origin!',
     &      '            Use the .NOCMC option to force this...',
     &      '            OECD calculation is cancelled.'
            LCOECD = .TRUE.
         END IF
      END IF
      IF (LOCDBG) THEN
         AVEC(1) = 0.0D0
         AVEC(2) = 0.0D0
         AVEC(3) = 2.063628D0
      END IF
C
C     Write transition moments, dip. strength, and rot. strength
C     to outputfile.
C     Transitions moments are keept as <O|Operator|f>. Forming
C     products like <O|Operator1|f><f|Operator2|O> we must therefore
C     multiply with a minus sign if Operator2 is antihermitian.
C     Notice: <O|r|f>        =  <f|r|O>
C             <O|grad|f>     = -<f|grad|O>
C             <O|r X grad|f> = -<f|r X grad|O>
C
C
C     Write electric transition dipole moments to output:
C     ====================================================
C
C-tbp IF (.NOT. (DIPSTR .OR. ROTSTR .OR. ROTVEL)) THEN
         IF (EXCTRP) THEN
            CALL HEADER('Triplet electronic excitation energies',15)
         ELSE
            CALL HEADER('Singlet electronic excitation energies',15)
         END IF
         WRITE (LUPRI,'(16X,A/16X,A/16X,A)')
     &      ' Sym.   Mode   Frequency    Frequency',
     &      'ex. st.  No.      (au)          (eV)',
     &      '======================================='
         DO 15 ISYM = 1, NSYM
            DO 14 IEXVAL = 1, NEXCIT(ISYM)
               WRITE (LUPRI,'(I20,I9,2F12.6)')
     &            ISYM,IEXVAL,EXENG(ISYM,IEXVAL),
     &            EXENG(ISYM,IEXVAL)*XTEV
                  LAPRPC = 'Excita  '
                  ISPINEX = 1
                  IF (EXCTRP) ISPINEX = 3
                  CALL WRIPRO(EXENG(ISYM,IEXVAL),"SCF/DFT   ",0,
     *                        LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *                        0.0D0,0.0D0,0.0D0,1,ISYM,ISPINEX,IEXVAL)
 14         CONTINUE
            WRITE(LUPRI,'(16X,A)')
     &      '---------------------------------------'
 15      CONTINUE
C-tbp END IF
      IF (DIPSTR .OR. ROTSTR .OR. ROTVEL) THEN
         CALL HEADER('Electric transition dipole moments (in a.u.)',15)
         WRITE (LUPRI,'(A,2(/A))')
     &   '  Sym.   Mode    Frequency       Velocity/Frequency  '//
     &   '            Length',
     &   ' ex. st.  No.      (au)          x       y       z   '//
     &   '      x       y       z',
     &   ' ===================================================='//
     &   '=========================='
         DO ISYM = 1, NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)
               WRITE (LUPRI,'(I4,I9,F13.6,3X,3F8.4,2X,3F8.4)')
     &            ISYM, IEXVAL, EXENG(ISYM,IEXVAL),
     &            TRVEL(1,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TRVEL(2,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TRVEL(3,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TRLEN(1,ISYM,IEXVAL),  TRLEN(2,ISYM,IEXVAL),
     &            TRLEN(3,ISYM,IEXVAL)
            END DO
            WRITE (LUPRI,'(A)')
     &   ' ----------------------------------------------------'//
     &   '--------------------------'
         END DO
      END IF
C
C     Write electric quadrupole transition moments to output:
C     =======================================================
C
CSPAS:23/5-11: second and third moment sum rules
C     IF (OECD) THEN
      IF (OECD .OR. SUMRUL) THEN
CKeinSPASmehr
         CALL HEADER('Diagonal electric quadrupole transition '
     &               //'moments (in a.u.)',15)
         WRITE (LUPRI,'(A,A,2(/A,A))')
     &   '  Sym.   Mode    Frequency       Velocity/Frequency  ',
     &   '            Length',
     &   ' ex. st.  No.      (au)          xx      yy      zz  ',
     &   '      xx      yy      zz',
     &   ' ----------------------------------------------------',
     &   '--------------------------'
         DO ISYM = 1, NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)
               WRITE (LUPRI,'(I4,I9,F13.6,3X,3F8.4,2X,3F8.4)')
     &            ISYM, IEXVAL, EXENG(ISYM,IEXVAL),
     &            TQVEL(1,1,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TQVEL(2,2,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TQVEL(3,3,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TQLEN(1,1,ISYM,IEXVAL),
     &            TQLEN(2,2,ISYM,IEXVAL),
     &            TQLEN(3,3,ISYM,IEXVAL)
            END DO
         END DO
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(//,A,3F10.5)')
     &      'Velocity after translation along a = ',(AVEC(I),I=1,3)
            DO ISYM = 1,NSYM
               DO IEXVAL = 1,NEXCIT(ISYM)
                  XNEW(1) = TQVEL(1,1,ISYM,IEXVAL)
     &                    - 2.0D0*AVEC(1)*TRVEL(1,ISYM,IEXVAL)
                  XNEW(2) = TQVEL(2,2,ISYM,IEXVAL)
     &                    - 2.0D0*AVEC(2)*TRVEL(2,ISYM,IEXVAL)
                  XNEW(3) = TQVEL(3,3,ISYM,IEXVAL)
     &                    - 2.0D0*AVEC(3)*TRVEL(3,ISYM,IEXVAL)
                  WRITE(LUPRI,'(2X,I2,6X,I3,1X,F12.6,3X,3F12.6)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            XNEW(1)/EXENG(ISYM,IEXVAL),
     &            XNEW(2)/EXENG(ISYM,IEXVAL),
     &            XNEW(3)/EXENG(ISYM,IEXVAL)
               END DO
            END DO
            WRITE(LUPRI,'(/,A,3F10.5)')
     &      'Length after translation along a = ',(AVEC(I),I=1,3)
            DO ISYM = 1,NSYM
               DO IEXVAL = 1,NEXCIT(ISYM)
                  XNEW(1) = TQLEN(1,1,ISYM,IEXVAL)
     &                    - 2.0D0*AVEC(1)*TRLEN(1,ISYM,IEXVAL)
                  XNEW(2) = TQLEN(2,2,ISYM,IEXVAL)
     &                    - 2.0D0*AVEC(2)*TRLEN(2,ISYM,IEXVAL)
                  XNEW(3) = TQLEN(3,3,ISYM,IEXVAL)
     &                    - 2.0D0*AVEC(3)*TRLEN(3,ISYM,IEXVAL)
                  WRITE(LUPRI,'(2X,I2,6X,I3,1X,F12.6,3X,3F12.6)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            XNEW(1),XNEW(2),XNEW(3)
               END DO
            END DO
         END IF
         CALL HEADER('Off-diagonal electric quadrupole transition '
     &               //'moments (in a.u.)',15)
         WRITE (LUPRI,'(1X,A,2(/1X,A))')
     &   ' Sym.   Mode    Frequency       Velocity/Frequency  '//
     &   '            Length       ',
     &   'ex. st.  No.      (au)          xy      xz      yz  '//
     &   '      xy      xz      yz ',
     &   '=================================================='//
     &   '============================'
         DO ISYM = 1, NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)
               WRITE (LUPRI,'(2X,I2,6X,I3,1X,F12.6,3X,3F8.3,2X,3F8.3)')
     &            ISYM, IEXVAL, EXENG(ISYM,IEXVAL),
     &            TQVEL(1,2,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TQVEL(1,3,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TQVEL(2,3,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TQLEN(1,2,ISYM,IEXVAL),
     &            TQLEN(1,3,ISYM,IEXVAL),
     &            TQLEN(2,3,ISYM,IEXVAL)
            END DO
            WRITE(LUPRI,'(1X,A)')
     &   '--------------------------------------------------'//
     &   '----------------------------'
         END DO
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(//,A,3F10.5)')
     &      'Velocity after translation along a = ',(AVEC(I),I=1,3)
            DO ISYM = 1,NSYM
               DO IEXVAL = 1,NEXCIT(ISYM)
                  XNEW(1) = TQVEL(1,2,ISYM,IEXVAL)
     &                    - AVEC(1)*TRVEL(2,ISYM,IEXVAL)
     &                    - AVEC(2)*TRVEL(1,ISYM,IEXVAL)
                  XNEW(2) = TQVEL(1,3,ISYM,IEXVAL)
     &                    - AVEC(1)*TRVEL(3,ISYM,IEXVAL)
     &                    - AVEC(3)*TRVEL(1,ISYM,IEXVAL)
                  XNEW(3) = TQVEL(2,3,ISYM,IEXVAL)
     &                    - AVEC(2)*TRVEL(3,ISYM,IEXVAL)
     &                    - AVEC(3)*TRVEL(2,ISYM,IEXVAL)
                  WRITE(LUPRI,'(2X,I2,6X,I3,1X,F12.6,3X,3F12.6)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            XNEW(1)/EXENG(ISYM,IEXVAL),
     &            XNEW(2)/EXENG(ISYM,IEXVAL),
     &            XNEW(3)/EXENG(ISYM,IEXVAL)
               END DO
            END DO
            WRITE(LUPRI,'(/,A,3F10.5)')
     &      'Length after translation along a = ',(AVEC(I),I=1,3)
            DO ISYM = 1,NSYM
               DO IEXVAL = 1,NEXCIT(ISYM)
                  XNEW(1) = TQLEN(1,2,ISYM,IEXVAL)
     &                    - AVEC(1)*TRLEN(2,ISYM,IEXVAL)
     &                    - AVEC(2)*TRLEN(1,ISYM,IEXVAL)
                  XNEW(2) = TQLEN(1,3,ISYM,IEXVAL)
     &                    - AVEC(1)*TRLEN(3,ISYM,IEXVAL)
     &                    - AVEC(3)*TRLEN(1,ISYM,IEXVAL)
                  XNEW(3) = TQLEN(2,3,ISYM,IEXVAL)
     &                    - AVEC(2)*TRLEN(3,ISYM,IEXVAL)
     &                    - AVEC(3)*TRLEN(2,ISYM,IEXVAL)
                  WRITE(LUPRI,'(2X,I2,6X,I3,1X,F12.6,3X,3F12.6)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            XNEW(1),XNEW(2),XNEW(3)
               END DO
            END DO
            WRITE(LUPRI,*)
            WRITE(LUPRI,*)
         END IF
      END IF
CSPAS:23/5-11: second and third moment sum rules
      IF (SUMRUL) THEN
         CALL HEADER('Electric octapole transition '
     &               //'moments (in a.u.)',15)
         WRITE (LUPRI,'(A,A,2(/A,A))')
     &   '  Sym.   Mode    Frequency                Length     ',
     &   '                  ',
     &   ' ex. st.  No.      (au)         xxx     xxy     xxz  ',
     &   '   xyy     xyz     xzz',
     &   ' ----------------------------------------------------',
     &   '--------------------------'
         DO ISYM = 1, NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)
               WRITE (LUPRI,'(I4,I9,F13.6,3X,6F8.4)')
     &            ISYM, IEXVAL, EXENG(ISYM,IEXVAL),
     &            TTLEN(1,ISYM,IEXVAL),
     &            TTLEN(4,ISYM,IEXVAL),
     &            TTLEN(5,ISYM,IEXVAL),
     &            TTLEN(6,ISYM,IEXVAL),
     &            TTLEN(7,ISYM,IEXVAL),
     &            TTLEN(8,ISYM,IEXVAL)
            END DO
         END DO
         WRITE (LUPRI,'(//,A,A,2(/A,A))')
     &   '  Sym.   Mode    Frequency                Length     ',
     &   '                  ',
     &   ' ex. st.  No.      (au)         yyy     yyz     yzz  ',
     &   '   zzz     ',
     &   ' ----------------------------------------------------',
     &   '----------'
         DO ISYM = 1, NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)
               WRITE (LUPRI,'(I4,I9,F13.6,3X,4F8.4)')
     &            ISYM, IEXVAL, EXENG(ISYM,IEXVAL),
     &            TTLEN(2,ISYM,IEXVAL),
     &            TTLEN(9,ISYM,IEXVAL),
     &            TTLEN(10,ISYM,IEXVAL),
     &            TTLEN(3,ISYM,IEXVAL)
            END DO
         END DO
      END IF
CKeinSPASmehr
C
C     Write magnetic transition dipole moments to output:
C     ===================================================
C
      IF (ROTSTR .OR. ROTVEL) THEN
         CALL HEADER('Magnetic transition dipole moments (in a.u.)',15)
         WRITE (LUPRI,'(A//A/A/A)')
     &   '  ( mu_B*<0|l_i|n>, where mu_B = 0.5 is the Bohr magneton)',
     &   ' Symm.  Mode   Frequency         Conventional     '//
     &   '             London',
     &   ' ex.st.  No.     (eV)        x        y       z   '//
     &   '      x        y        z',
     &   ' =================================================='//
     &   '============================'
C
C      The factor 0.5 is the Bohr-magneton
C      -tbp: correct signs! was FAC_BM = +0.5D0
         FAC_BM = -0.5D0
C
         DO 400 ISYM = 1,NSYM
            DO 300 IEXVAL = 1,NEXCIT(ISYM)
               WRITE (LUPRI,'(I4,I7,F11.4,2X,3F9.5,1X,3F9.5)')
     &    ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &    FAC_BM*TRMAG(1,ISYM,IEXVAL), FAC_BM*TRMAG(2,ISYM,IEXVAL),
     &    FAC_BM*TRMAG(3,ISYM,IEXVAL), FAC_BM*TRLON(1,ISYM,IEXVAL),
     &    FAC_BM*TRLON(2,ISYM,IEXVAL), FAC_BM*TRLON(3,ISYM,IEXVAL)
               IF (NODIFC) THEN
                  WRITE (LUPRI,'(51X,SP,3F9.5/33X,A,S,3F9.5/51X,A)')
     &            FAC_BM*BSRLON(1,ISYM,IEXVAL),
     &            FAC_BM*BSRLON(2,ISYM,IEXVAL),
     &            FAC_BM*BSRLON(3,ISYM,IEXVAL),
     & 'Total:  ',FAC_BM*(TRLON(1,ISYM,IEXVAL)+BSRLON(1,ISYM,IEXVAL)),
     &            FAC_BM*(TRLON(2,ISYM,IEXVAL)+BSRLON(2,ISYM,IEXVAL)),
     &            FAC_BM*(TRLON(3,ISYM,IEXVAL)+BSRLON(3,ISYM,IEXVAL)),
     &       '----------------------------'
               END IF
  300       CONTINUE
            WRITE (LUPRI,'(A)') ' ----------------------------------'//
     &         '--------------------------------------------'
  400    CONTINUE
         IF (LOCDBG) THEN
            WRITE(LUPRI,'(//,A,3F12.5)')
     &      'Conventional after translation along a = ',(AVEC(I),I=1,3)
            DO ISYM = 1,NSYM
               DO IEXVAL = 1,NEXCIT(ISYM)
                  XNEW(1) = TRMAG(1,ISYM,IEXVAL)
     &                    - AVEC(2)*TRVEL(3,ISYM,IEXVAL)
     &                    + AVEC(3)*TRVEL(2,ISYM,IEXVAL)
                  XNEW(2) = TRMAG(2,ISYM,IEXVAL)
     &                    - AVEC(3)*TRVEL(1,ISYM,IEXVAL)
     &                    + AVEC(1)*TRVEL(3,ISYM,IEXVAL)
                  XNEW(3) = TRMAG(3,ISYM,IEXVAL)
     &                    - AVEC(1)*TRVEL(2,ISYM,IEXVAL)
     &                    + AVEC(2)*TRVEL(1,ISYM,IEXVAL)
                  WRITE(LUPRI,'(2X,I2,6X,I3,1X,F12.6,3X,3F12.6)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            FAC_BM*XNEW(1),FAC_BM*XNEW(2),FAC_BM*XNEW(3)
               END DO
            END DO
            WRITE(LUPRI,'(/)')
         END IF
      END IF
C
C     Write oscilator and rotational strengths to output:
C     ===================================================
C
      IF (ROTSTR .OR. ROTVEL) THEN
         CALL HEADER('Oscillator and Scalar Rotational Strengths',
     &   18)
         WRITE (LUPRI,'(2(A/))')
     &      '  Units: 10**(-40) (esu**2)*(cm**2) (rotational strength)',
     &      '         dimensionless              (oscillator strength)'
         WRITE (LUPRI,'(A/A/A)')
     &   ' Symm.  Mode   Frequency    Oscillator-strength '//
     &   '     Rotational-strength   ',
     &   ' ex.st.  No.     (eV)        velocity  length   '//
     &   '  velocity  length   London',
     &   ' ==============================================='//
     &   '==========================='
         DO 600 ISYM = 1,NSYM
            DO 500 IEXVAL = 1,NEXCIT(ISYM)
               DIPSTV = DDOT(3,TRVEL(1,ISYM,IEXVAL),1,
     &                  TRVEL(1,ISYM,IEXVAL),1) / EXENG(ISYM,IEXVAL)
               DIPSTL = DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                  TRLEN(1,ISYM,IEXVAL),1) * EXENG(ISYM,IEXVAL)
C              FAC_BM factor is the Bohr-magneton
               ROTV   = FAC_BM * DDOT(3,TRVEL(1,ISYM,IEXVAL),1,
     &                  TRMAG(1,ISYM,IEXVAL),1) / EXENG(ISYM,IEXVAL)
               ROTL   = FAC_BM * DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                  TRMAG(1,ISYM,IEXVAL),1)
               ROTLON = FAC_BM * DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                  TRLON(1,ISYM,IEXVAL),1)
               IF (NODIFC) THEN
                  ROTLO2 = FAC_BM * DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                     BSRLON(1,ISYM,IEXVAL),1)
                  ROTLON = ROTLON + ROTLO2
               ENDIF
               WRITE (LUPRI,'(I4,I7,F11.4,5X,2F9.4,3X,3F9.3)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            (DTWO/DTHREE)*DIPSTV, (DTWO/DTHREE)*DIPSTL,
     &            ESUECD*ROTV,   ESUECD*ROTL,   ESUECD*ROTLON
               IF (OECD) THEN
                  SVEL(ISYM,IEXVAL) = ESUECD*ROTV
                  SLEN(ISYM,IEXVAL) = ESUECD*ROTL
               END IF
  500       CONTINUE
            WRITE (LUPRI,'(A)')
     &   ' -------------------------------------------------'//
     &   '-------------------------'
  600    CONTINUE
C
       IF (OECD .AND. (.NOT.LCOECD)) THEN
            CALL HEADER('Magnetic Dipole Rotatory Strength Tensor',18)
            WRITE (LUPRI,'(A/)')
     &      '  Units: 10**(-40) (esu**2)*(cm**2)'
            WRITE (LUPRI,'(3(/1X,A))')
     &   ' Sym. Mode Gauge         Magnetic Dipole Rotatory Strength '//
     &   'Tensor          ',
     &   '       No.            xx        xy        xz        yy     '//
     &   '   yz        zz ',
     &   '==========================================================='//
     &   '================'
         DO ISYM = 1,NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)

               PDOTL = DDOT(3,TRVEL(1,ISYM,IEXVAL),1,
     &                        TRMAG(1,ISYM,IEXVAL),1)
               RDOTL = DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                        TRMAG(1,ISYM,IEXVAL),1)

               DO J = 1,3
                  DO I = 1,3
                     RMVEL(I,J,ISYM,IEXVAL) =
     &                        -TRVEL(J,ISYM,IEXVAL)*TRMAG(I,ISYM,IEXVAL)
                     RMLEN(I,J,ISYM,IEXVAL) =
     &                        -TRLEN(J,ISYM,IEXVAL)*TRMAG(I,ISYM,IEXVAL)
                     IF (I .EQ. J) THEN
                        RMVEL(I,I,ISYM,IEXVAL) = RMVEL(I,I,ISYM,IEXVAL)
     &                                         + PDOTL
                        RMLEN(I,I,ISYM,IEXVAL) = RMLEN(I,I,ISYM,IEXVAL)
     &                                         + RDOTL
                     END IF
                  END DO
               END DO

               FACL = 0.375D0*ESUECD
               FACV = FACL/EXENG(ISYM,IEXVAL)
               CALL DSCAL(3*3,FACV,RMVEL(1,1,ISYM,IEXVAL),1)
               CALL DSCAL(3*3,FACL,RMLEN(1,1,ISYM,IEXVAL),1)

               TSTVEL = 0.0D0
               TSTLEN = 0.0D0
               DO J = 1,3
                  DO I = 1,J
                     RMVEL(I,J,ISYM,IEXVAL) = RMVEL(I,J,ISYM,IEXVAL)
     &                                      + RMVEL(J,I,ISYM,IEXVAL)
                     RMVEL(J,I,ISYM,IEXVAL) = RMVEL(I,J,ISYM,IEXVAL)
                     RMLEN(I,J,ISYM,IEXVAL) = RMLEN(I,J,ISYM,IEXVAL)
     &                                      + RMLEN(J,I,ISYM,IEXVAL)
                     RMLEN(J,I,ISYM,IEXVAL) = RMLEN(I,J,ISYM,IEXVAL)
                  END DO
                  TSTVEL = TSTVEL + RMVEL(J,J,ISYM,IEXVAL)
                  TSTLEN = TSTLEN + RMLEN(J,J,ISYM,IEXVAL)
               END DO
               TSTVEL = TSTVEL/3.0D0
               TSTLEN = TSTLEN/3.0D0

               WRITE(LUPRI,'(2X,I2,3X,I4,2X,A3,6(1X,F9.3))')
     &         ISYM,IEXVAL,'Vel',RMVEL(1,1,ISYM,IEXVAL),
     &                           RMVEL(1,2,ISYM,IEXVAL),
     &                           RMVEL(1,3,ISYM,IEXVAL),
     &                           RMVEL(2,2,ISYM,IEXVAL),
     &                           RMVEL(2,3,ISYM,IEXVAL),
     &                           RMVEL(3,3,ISYM,IEXVAL)
               IF (ABS(TSTVEL-SVEL(ISYM,IEXVAL)) .GT. TOLRTT) THEN
                  WRITE(LUPRI,'(1X,A,1P,D16.8,/,1X,A,1P,D16.8)')
     &            'WARNING: wrong average: ',TSTVEL,
     &            '         expected     : ',SVEL(ISYM,IEXVAL)
                  LOCWRN = LOCWRN + 1
               END IF
               WRITE(LUPRI,'(13X,A3,6(1X,F9.3))')
     &                     'Len',RMLEN(1,1,ISYM,IEXVAL),
     &                           RMLEN(1,2,ISYM,IEXVAL),
     &                           RMLEN(1,3,ISYM,IEXVAL),
     &                           RMLEN(2,2,ISYM,IEXVAL),
     &                           RMLEN(2,3,ISYM,IEXVAL),
     &                           RMLEN(3,3,ISYM,IEXVAL)
               IF (ABS(TSTLEN-SLEN(ISYM,IEXVAL)) .GT. TOLRTT) THEN
                  WRITE(LUPRI,'(1X,A,1P,D16.8,/,1X,A,1P,D16.8)')
     &            'WARNING: wrong average: ',TSTVEL,
     &            '         expected     : ',SVEL(ISYM,IEXVAL)
                  LOCWRN = LOCWRN + 1
               END IF

            END DO
            WRITE (LUPRI,'(1X,A)')
     &   '-----------------------------------------------------------'//
     &   '----------------'
         END DO

         CALL HEADER('Electric Quadrupole Rotatory Strength Tensor',18)
         WRITE (LUPRI,'(A/)')
     &      '  Units: 10**(-40) (esu**2)*(cm**2)'
         WRITE (LUPRI,'(3(/1X,A))')
     &   ' Sym. Mode Gauge     Electric Quadrupole Rotatory Strength '//
     &   'Tensor          ',
     &   '       No.            xx        xy        xz        yy     '//
     &   '   yz        zz ',
     &   '==========================================================='//
     &   '================'
         DO ISYM = 1,NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)

               DO K = 1,3
                  DO M = 1,3
                     DO L = 1,3
                        DO J = 1,3
                           RQVEL(J,K,ISYM,IEXVAL) =
     &                     RQVEL(J,K,ISYM,IEXVAL) + XLEVICI(J,L,M)*
     &                       TRVEL(L,ISYM,IEXVAL)*TQVEL(M,K,ISYM,IEXVAL)
                           RQLEN(J,K,ISYM,IEXVAL) =
     &                     RQLEN(J,K,ISYM,IEXVAL) + XLEVICI(J,L,M)*
     &                       TRLEN(L,ISYM,IEXVAL)*TQLEN(M,K,ISYM,IEXVAL)
                        END DO
                     END DO
                  END DO
               END DO

               FACV = -0.375D0*ESUECD/EXENG(ISYM,IEXVAL)
               FACL = -0.375D0*ESUECD*EXENG(ISYM,IEXVAL)
               CALL DSCAL(3*3,FACV,RQVEL(1,1,ISYM,IEXVAL),1)
               CALL DSCAL(3*3,FACL,RQLEN(1,1,ISYM,IEXVAL),1)

               TSTVEL = 0.0D0
               TSTLEN = 0.0D0
               DO J = 1,3
                  DO I = 1,J
                     RQVEL(I,J,ISYM,IEXVAL) = RQVEL(I,J,ISYM,IEXVAL)
     &                                      + RQVEL(J,I,ISYM,IEXVAL)
                     RQVEL(J,I,ISYM,IEXVAL) = RQVEL(I,J,ISYM,IEXVAL)
                     RQLEN(I,J,ISYM,IEXVAL) = RQLEN(I,J,ISYM,IEXVAL)
     &                                      + RQLEN(J,I,ISYM,IEXVAL)
                     RQLEN(J,I,ISYM,IEXVAL) = RQLEN(I,J,ISYM,IEXVAL)
                  END DO
                  TSTVEL = TSTVEL + RQVEL(J,J,ISYM,IEXVAL)
                  TSTLEN = TSTLEN + RQLEN(J,J,ISYM,IEXVAL)
               END DO
               TSTVEL = TSTVEL/3.0D0
               TSTLEN = TSTLEN/3.0D0

               WRITE(LUPRI,'(2X,I2,3X,I4,2X,A3,6(1X,F9.3))')
     &         ISYM,IEXVAL,'Vel',RQVEL(1,1,ISYM,IEXVAL),
     &                           RQVEL(1,2,ISYM,IEXVAL),
     &                           RQVEL(1,3,ISYM,IEXVAL),
     &                           RQVEL(2,2,ISYM,IEXVAL),
     &                           RQVEL(2,3,ISYM,IEXVAL),
     &                           RQVEL(3,3,ISYM,IEXVAL)
               IF (ABS(TSTVEL) .GT. TOLRTT) THEN
                  WRITE(LUPRI,'(1X,A,1P,D16.8)')
     &            'WARNING: wrong average: ',TSTVEL
                  LOCWRN = LOCWRN + 1
               END IF
               WRITE(LUPRI,'(13X,A3,6(1X,F9.3))')
     &                     'Len',RQLEN(1,1,ISYM,IEXVAL),
     &                           RQLEN(1,2,ISYM,IEXVAL),
     &                           RQLEN(1,3,ISYM,IEXVAL),
     &                           RQLEN(2,2,ISYM,IEXVAL),
     &                           RQLEN(2,3,ISYM,IEXVAL),
     &                           RQLEN(3,3,ISYM,IEXVAL)
               IF (ABS(TSTLEN) .GT. TOLRTT) THEN
                  WRITE(LUPRI,'(1X,A,1P,D16.8)')
     &            'WARNING: wrong average: ',TSTLEN
                  LOCWRN = LOCWRN + 1
               END IF

            END DO
            WRITE(LUPRI,'(1X,A)')
     &   '-----------------------------------------------------------'//
     &   '----------------'
         END DO

         CALL HEADER('Total Rotatory Strength Tensor',18)
         WRITE (LUPRI,'(A/)')
     &      '  Units: 10**(-40) (esu**2)*(cm**2)'
         WRITE (LUPRI,'(3(/,1X,A,A))')
     &   ' Sym. Mode Gauge                   Total Rotatory Strength ',
     &   'Tensor          ',
     &   '       No.            xx        xy        xz        yy     ',
     &   '   yz        zz ',
     &   '-----------------------------------------------------------',
     &   '----------------'
         LENTOT = 3*3*NSYM*MXNEXI
         CALL DCOPY(LENTOT,RMVEL,1,RVEL,1)
         CALL DAXPY(LENTOT,1.0D0,RQVEL,1,RVEL,1)
         CALL DCOPY(LENTOT,RMLEN,1,RLEN,1)
         CALL DAXPY(LENTOT,1.0D0,RQLEN,1,RLEN,1)
         DO ISYM = 1,NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)

               TSTVEL = RVEL(1,1,ISYM,IEXVAL)
               TSTLEN = RLEN(1,1,ISYM,IEXVAL)
               DO I = 2,3
                  TSTVEL = TSTVEL + RVEL(I,I,ISYM,IEXVAL)
                  TSTLEN = TSTLEN + RLEN(I,I,ISYM,IEXVAL)
               END DO
               TSTVEL = TSTVEL/3.0D0
               TSTLEN = TSTLEN/3.0D0

               WRITE(LUPRI,'(2X,I2,3X,I4,2X,A3,6(1X,F9.3))')
     &         ISYM,IEXVAL,'Vel',RVEL(1,1,ISYM,IEXVAL),
     &                           RVEL(1,2,ISYM,IEXVAL),
     &                           RVEL(1,3,ISYM,IEXVAL),
     &                           RVEL(2,2,ISYM,IEXVAL),
     &                           RVEL(2,3,ISYM,IEXVAL),
     &                           RVEL(3,3,ISYM,IEXVAL)
               IF (ABS(TSTVEL-SVEL(ISYM,IEXVAL)) .GT. TOLRTT) THEN
                  WRITE(LUPRI,'(1X,A,1P,D16.8,/,1X,A,1P,D16.8)')
     &            'WARNING: wrong average: ',TSTVEL,
     &            '         expected     : ',SVEL(ISYM,IEXVAL)
                  LOCWRN = LOCWRN + 1
               END IF
               WRITE(LUPRI,'(13X,A3,6(1X,F9.3))')
     &                     'Len',RLEN(1,1,ISYM,IEXVAL),
     &                           RLEN(1,2,ISYM,IEXVAL),
     &                           RLEN(1,3,ISYM,IEXVAL),
     &                           RLEN(2,2,ISYM,IEXVAL),
     &                           RLEN(2,3,ISYM,IEXVAL),
     &                           RLEN(3,3,ISYM,IEXVAL)
               IF (ABS(TSTLEN-SLEN(ISYM,IEXVAL)) .GT. TOLRTT) THEN
                  WRITE(LUPRI,'(1X,A,1P,D16.8,/,1X,A,1P,D16.8)')
     &            'WARNING: wrong average: ',TSTLEN,
     &            '         expected     : ',SLEN(ISYM,IEXVAL)
                  LOCWRN = LOCWRN + 1
               END IF

            END DO
         END DO

         IF (LOCWRN .NE. 0) THEN
            WRITE(LUPRI,'(//A,I6,A/A//)') ' WARNING:',LOCWRN,
     &      ' warnings were issued for rotatory strength tensors.',
     &      ' WARNING: Averages are incorrect!!!!'
         END IF

       END IF  ! (OECD .AND. (.NOT.LCOECD))
C
C     next line is end of "IF (ROTSTR .OR. ROTVEL) THEN"
      ELSE IF (DIPSTR .AND. (.NOT. OECD)) THEN
         CALL HEADER('Oscillator strengths',30)
         WRITE (LUPRI,'(A)')
     &      '  Oscillator strengths are dimensionless.'
         WRITE (LUPRI,'(3(/A))')
     &   '  Sym.   Mode        Frequency     Oscillator-strength ',
     &   ' ex. st.  No.           (eV)        velocity   length   ',
     &   ' -------------------------------------------------------'
         DO 800 ISYM = 1,NSYM
            DO 700 IEXVAL = 1,NEXCIT(ISYM)
               DIPSTV = DDOT(3,TRVEL(1,ISYM,IEXVAL),1,
     &                  TRVEL(1,ISYM,IEXVAL),1) / EXENG(ISYM,IEXVAL)
               DIPSTL = DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                  TRLEN(1,ISYM,IEXVAL),1) * EXENG(ISYM,IEXVAL)
               WRITE (LUPRI,'(2X,I2,6X,I3,5X,F12.6,5X,2F9.4)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            (DTWO/DTHREE)*DIPSTV, (DTWO/DTHREE)*DIPSTL
                  LAPRPC = 'OSCI-VEL'
                  CALL WRIPRO((DTWO/DTHREE)*DIPSTV,"SCF/DFT   ",-400,
     *                         LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *                         0.0D0,0.0D0,0.0D0,1,ISYM,1,IEXVAL)
                  LAPRPC = 'OSCI-LEN'
                  CALL WRIPRO((DTWO/DTHREE)*DIPSTL,"SCF/DFT   ",-400,
     *                         LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *                         0.0D0,0.0D0,0.0D0,1,ISYM,1,IEXVAL)

  700       CONTINUE
  800    CONTINUE
      END IF
C
      IF (FNAC) THEN
         NCOOR  = 3*NUCDEP
         KCSTRA = 1
         KSCTRA = KCSTRA + NCOOR*NCOOR
         KLAST  = KSCTRA + NCOOR*NCOOR
         IF (KLAST .GT. LWORK) CALL STOPIT('EXCOUT',' ',KLAST,LWORK)
         CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),1,NCOOR,0)
         CALL HEADER
     &   ('Non-adiabatic coupling matrix elements <O|d/dR|f> (in a.u.)',
     &    9)
         WRITE (LUPRI,'(1X,A,A,2(/,1X,A,A))')
     &   ' Coordinate R          Sym.  Mode Ex. energy    WFR-term ',
     &   ' BSR-term   Total ',
     &   ' Cartesian no.          st.   No.    (Ev)                ',
     &   '                 ',
     &   '----------------------------------------------------------',
     &   '------------------'
         DO 1100 ISYM = 1, NSYM
            DO 1000 ISCOOR = 1, 3*NATOMS
               FX = ZERO
               ICOUNT = 0
               DO 900 IEXVAL = 1,NEXCIT(ISYM)
                  CALL GENCOR(WORK(KCSTRA),NCOOR,ISCOOR,COORDI,ICORSY)
                  IF (ICORSY .EQ. ISYM) THEN
                     WRITE (LUPRI,'(1X,A,2X,I2,4X,I2,2X,F10.6,1X,
     &                      2F10.4,F12.6)')
     &               COORDI, ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &               FONAC(ISCOOR,ISYM,IEXVAL), 
     &               FONA2(ISCOOR,ISYM,IEXVAL),
     &               FONAC(ISCOOR,ISYM,IEXVAL)+FONA2(ISCOOR,ISYM,IEXVAL)
C
                     FX = FX + FONAC(ISCOOR,ISYM,IEXVAL)**2 /
     &                         EXENG(ISYM,IEXVAL)
                     ICOUNT = ICOUNT + 1
                  END IF
  900          CONTINUE
C
               IF (ICOUNT .GT. 0) WRITE(LUPRI,'(/,1X,A,A,22X,1F11.6,/)') 
     &                                  'Vibrational g-factor: ',COORDI,
     &                                  FX*DTWO
C
 1000       CONTINUE
 1100    CONTINUE
      END IF
C
C     Calculate oscillator strength sum rules
C     =======================================
C
      IF (SUMRUL) THEN
         CALL DZERO(DSSUML,4*9)
CSPAS:23/5-11: second and third moment sum rules
         CALL DZERO(QSSUM,4*9)
         CALL DZERO(DTSSUM,4*9)
CKeinSPASmehr
         CALL DZERO(DLSUML,4*9)
         CALL DZERO(DSSUMM,4*9)
         CALL DZERO(DLSUMM,4*9)
         CALL DZERO(DSSUMV,4*9)
         CALL DZERO(DLSUMV,4*9)
         DO ISYM = 1,NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)
               DO ICOM = 1,3
                  DO K = -6,2
                     DSSUML(K,ICOM) = DSSUML(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K+1)
     &                              * TRLEN(ICOM,ISYM,IEXVAL) 
     &                              * TRLEN(ICOM,ISYM,IEXVAL)
CSPAS:23/5-11: second and third moment sum rules
                     QSSUM(K,ICOM) = QSSUM(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K+1)
     &                              * TQLEN(ICOM,ICOM,ISYM,IEXVAL) 
     &                              * TQLEN(ICOM,ICOM,ISYM,IEXVAL)
                     DTSSUM(K,ICOM) = DTSSUM(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K+1)
     &                              * TRLEN(ICOM,ISYM,IEXVAL) 
     &                              * TTLEN(ICOM,ISYM,IEXVAL)
CKeinSPASmehr
                     DLSUML(K,ICOM) = DLSUML(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K+1)
     &                              * DLOG(EXENG(ISYM,IEXVAL))
     &                              * TRLEN(ICOM,ISYM,IEXVAL) 
     &                              * TRLEN(ICOM,ISYM,IEXVAL)
                     DSSUMM(K,ICOM) = DSSUMM(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**K
     &                              * TRLEN(ICOM,ISYM,IEXVAL) 
     &                              * TRVEL(ICOM,ISYM,IEXVAL)
                     DLSUMM(K,ICOM) = DLSUMM(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**K
     &                              * DLOG(EXENG(ISYM,IEXVAL))
     &                              * TRLEN(ICOM,ISYM,IEXVAL) 
     &                              * TRVEL(ICOM,ISYM,IEXVAL)
                     DSSUMV(K,ICOM) = DSSUMV(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K-1)
     &                              * TRVEL(ICOM,ISYM,IEXVAL) 
     &                              * TRVEL(ICOM,ISYM,IEXVAL)
                     DLSUMV(K,ICOM) = DLSUMV(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K-1)
     &                              * DLOG(EXENG(ISYM,IEXVAL))
     &                              * TRVEL(ICOM,ISYM,IEXVAL) 
     &                              * TRVEL(ICOM,ISYM,IEXVAL)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO K = -6,2
            DSSUML(K,4) = (DSSUML(K,1)+DSSUML(K,2)+DSSUML(K,3)) / DTHREE
CSPAS:23/5-11: second and third moment sum rules
            QSSUM(K,4)  = (QSSUM(K,1)+QSSUM(K,2)+QSSUM(K,3)) / DTHREE
            DTSSUM(K,4) = (DTSSUM(K,1)+DTSSUM(K,2)+DTSSUM(K,3)) / DTHREE
CKeinSPASmehr
            DLSUML(K,4) = (DLSUML(K,1)+DLSUML(K,2)+DLSUML(K,3)) / DTHREE
            DSSUMM(K,4) = (DSSUMM(K,1)+DSSUMM(K,2)+DSSUMM(K,3)) / DTHREE
            DLSUMM(K,4) = (DLSUMM(K,1)+DLSUMM(K,2)+DLSUMM(K,3)) / DTHREE
            DSSUMV(K,4) = (DSSUMV(K,1)+DSSUMV(K,2)+DSSUMV(K,3)) / DTHREE
            DLSUMV(K,4) = (DLSUMV(K,1)+DLSUMV(K,2)+DLSUMV(K,3)) / DTHREE
            DO ICOM = 1,4
               IF (DSSUML(K,ICOM).EQ.ZERO) THEN
                  DISUML(K,ICOM) = ZERO
               ELSE
                  DISUML(K,ICOM) = DEXP(DLSUML(K,ICOM)/DSSUML(K,ICOM))
     &                             *XTEV
               ENDIF
               IF (DSSUMM(K,ICOM).EQ.ZERO) THEN
                  DISUMM(K,ICOM) = ZERO
               ELSE
                  DISUMM(K,ICOM) = DEXP(DLSUMM(K,ICOM)/DSSUMM(K,ICOM))
     &                             *XTEV
               ENDIF
               IF (DSSUMV(K,ICOM).EQ.ZERO) THEN
                  DISUMV(K,ICOM) = ZERO
               ELSE
                  DISUMV(K,ICOM) = DEXP(DLSUMV(K,ICOM)/DSSUMV(K,ICOM))
     &                             *XTEV
               ENDIF
            ENDDO
         ENDDO
C
         CALL HEADER('Oscillator strength sum rules',30)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'S(K) Sum Rules : Dipole Length Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DSSUML(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'S(K) Sum Rules : Dipole Mixed Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DSSUMM(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'S(K) Sum Rules : Dipole Velocity Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DSSUMV(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'L(K) Sum Rules : Dipole Length Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DLSUML(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'L(K) Sum Rules : Dipole Mixed Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DLSUMM(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'L(K) Sum Rules : Dipole Velocity Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DLSUMV(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'I(K) Sum Rules : Dipole Length Approximation in eV',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DISUML(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'I(K) Sum Rules : Dipole Mixed Approximation in eV',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DISUMM(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'I(K) Sum Rules : Dipole Velocity Approximation in eV',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DISUMV(K,J),J=1,4),K=-6,2)
CSPAS:23/5-11: second and third moment sum rules
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,2X,A,2X,A,5X,A,/)')
     &   'S(K) Sum Rules : Second Moments in a.u.',
     &   'K','xx,xx component','yy,yy component',
     &   'zz,zz component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(QSSUM(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,4X,A,2X,A,2X,A,6X,A,/)')
     &   'S(K) Sum Rules : Mixed First-third moments in a.u.',
     &   'K','x,xxx component','y,yyy component',
     &   'z,zzz component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DTSSUM(K,J),J=1,4),K=-6,2)
CKeinSPASmehr
      ENDIF
C     
C        *** If numerical derivatives are calculated, then ***
C        *** write properties to file.                     ***
C
      IF (NMWALK) THEN
         DO ISYM = 1, NSYM
            NMDTXT = '         '
            DO IEXVAL=1, NEXCIT(ISYM)
               IF (IEXVAL.LT.10) THEN
                  WRITE (NMDTXT(5:6),'(2I1)') ISYM, IEXVAL
               ELSE
                  WRITE (NMDTXT(5:7),'(I1,I2)') ISYM, IEXVAL
               END IF
               NMDTXT(1:4) = 'EXEN'
               CALL WRAVFL(EXENG(ISYM,IEXVAL),1,1,1,NMDTXT,IPREXC)
               IF (DIPSTR) THEN
                  NMDTXT(1:4) = 'TRDP'
                  CALL WRAVFL(TRLEN(1,ISYM,IEXVAL),3,1,1,NMDTXT,IPREXC)
               END IF
            END DO
         END DO
      END IF

      RETURN
      END
C  /* Deck excit2 */
      SUBROUTINE EXCIT2(EXVAL,NEXENG,ISYM,LUSOVE,WORK,LWORK)
C
C     Triplet excitations now determined in input
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <dummy.h>
#include <priunit.h>
      CHARACTER*8 LABEL
      LOGICAL INTTRS, CICLC, HFCLC, TRIPLE, EXECLC
      DIMENSION WORK(LWORK), EXVAL(NEXENG)
C
#include <cbiexc.h>
#include <infvar.h>
#include <inforb.h>
C
      CICLC  = .FALSE.
      HFCLC  = NASHT .LE. 1
      TRIPLE = EXCTRP
      EXECLC = .TRUE.
      NABATY = 1
      NABAOP = 1
      LABEL  = 'EXCITLAB'
C
C     The following line is needed in order to pass
C     consistency checks in RSPMC
C
      CALL ABAVAR(1,.FALSE.,0,WORK,LWORK)
      IF (HFCLC) NCONF = 1
C
      LUREVE = -9801
      LUGDVE = -9802
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     *            EXVAL,NEXENG,NABATY,NABAOP,LABEL,LUGDVE,LUSOVE,LUREVE,
     *            THREXC,MAXITE,IPREXC,MXRM,MXPHP,WORK,LWORK)
C
C     The following line is needed to reset variables
C
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
      CALL ABAVAR(ISYM,TRIPLE,0,WORK,LWORK)
C
      RETURN
      END
C  /* Deck labcop */
      SUBROUTINE LABCOP(NLAB,NLBTOT,LAB1,ISYM1,LAB2,ISYM2)
C
C     Purpose: To copy a set of NLAB labels and corresponding
C     symmetries from LAB1 and ISYM1 to LAB2 and ISYM2.
C
#include <implicit.h>
#include <priunit.h>
#include <cbiexc.h>
      CHARACTER LAB1(*)*8, LAB2(*)*8
      DIMENSION ISYM1(*), ISYM2(*)
C
      NTEST = NLBTOT + NLAB
      IF (NTEST.GT.MAXPP) THEN
         WRITE(LUPRI,'(A,/A,I6,A,I6)')
     *   ' NUMBER OF SPECIFIED PROPERTIES EXCEEDS THE MAXIMUM ALLOWED',
     *   ' MAXPRP =',MAXPP,' NLBTOT =',NTEST
         CALL QUIT(' LABCOP: TOO MANY PROPERTIES SPECIFIED')
      END IF
C
      DO 100 ILAB = 1, NLAB
         LAB2(NLBTOT + ILAB)  = LAB1(ILAB)
         ISYM2(NLBTOT + ILAB) = ISYM1(ILAB) + 1
  100 CONTINUE
C
C     Update number total number of labels
C
      NLBTOT = NTEST
C
      RETURN
      END
C  /* Deck gencor */
      SUBROUTINE GENCOR(CSTRA,NCOOR,ISCOOR,COORDI,ICORSY)
C
C     Purpose: To generate a readable description of the symmetry
C     coordinate ISCOOR in COORDI.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <nuclei.h>
      CHARACTER COORDI*21, CNBASE*2
      DIMENSION NBASE(8),NSIGN(8), CSTRA(NCOOR,NCOOR)
#include <chrxyz.h>
#include <chrnos.h>
#include <chrsgn.h>
C
      DO 400 ISYM = 0, MAXREP
      IF (NCRREP(ISYM,1) .GT. 0) THEN
         DO 300 IATOM = 1, NUCIND
            DO 200 ICOOR = 1, 3
               JSCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM,1)
               IF (JSCOOR.EQ.ISCOOR) THEN
                  ICORSY = ISYM + 1
                  NB = 0
                  DO 100 I = 1, NCOOR
                     NSGN = NINT(CSTRA(ISCOOR,I))
                     IF (NSGN .NE. 0) THEN
                        NB = NB + 1
                        NBASE(NB) = I
                        NSIGN(NB) = NSGN
                     END IF
  100             CONTINUE
                  IF (NB .EQ. 1) THEN
                     WRITE(CNBASE,'(I2)') NBASE(1)
                     COORDI = NAMEX(3*IATOM)(1:4)//CHRXYZ(-ICOOR)
     &                        //'['//CNBASE//']'
                  ELSE IF (NB .LT. 4) THEN
                     WRITE(CNBASE,'(I2)') NBASE(1)
                     COORDI = NAMEX(3*IATOM)(1:4)//CHRXYZ(-ICOOR)
     &                        //'['//CNBASE
                     J = 10
                     IF ( NB .GE. 2 ) THEN
                        WRITE(CNBASE,'(I2)') NBASE(2)
                        COORDI = COORDI(1:9)//CHRSGN(NSIGN(2))//CNBASE
                        COORDI = COORDI(1:12)//']/'//CHRNOS(2)
                     END IF
                     IF ( NB .GE. 3 ) THEN
                        WRITE(CNBASE,'(I2)') NBASE(3)
                        COORDI = COORDI(1:13)//CHRSGN(NSIGN(3))//CNBASE
                        COORDI = COORDI(1:16)//']/'//CHRNOS(2)
                     END IF
                  ELSE
                     COORDI = '    See above    '
                  END IF
               END IF
  200       CONTINUE
  300    CONTINUE
      END IF
  400 CONTINUE
C
      RETURN
      END
