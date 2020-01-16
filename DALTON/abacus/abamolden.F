!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
C     Original: abamolden  alex a auer 5.2.99
C
C This file contains all output routines that write parts of the molden.inp
C file which can be read by Molden.
C Info about Molden can be obtained at :
C        http://www.caos.kun.nl/~schaft/molden/molden.html
C
C
      SUBROUTINE MOLDEN_HEAD
#include "implicit.h"
#include "dummy.h"
#include "maxorb.h"
#include "inftap.h"
#include "molde.h"
C
      DONEIT=.FALSE.
      DONEIU=.FALSE.
      DONEIV=.FALSE.
      DONEIW=.FALSE.
      DONEIX=.FALSE.
      IXYZPRINTS = 0
      IF (LUMOLDEN .GT. 0) CALL GPCLOSE(LUMOLDEN,'DELETE') ! we may have changed basis set
      CALL GPOPEN(LUMOLDEN,'molden.inp','UNKNOWN',' ',
     &            'FORMATTED',IDUMMY,.FALSE.)
      REWIND LUMOLDEN
      WRITE(LUMOLDEN,'(A)') '[Molden Format]'
      CALL FLSHFO(LUMOLDEN)

!     open temporary file for MO coefficients; this file will always contain the last set saved by DALTON
!     When MOLDEN_TAIL is called, this file will be appended to the molden.inp file
!     (This is a lot easier than checking in SIRIUS and during geometry optimization which set of MO coefficients
!      is the final. We cannot use SIRIUS.RST or SIRIFC as they are now for this, because they do not always contain
!      the orbital energy and occupation information. /hjaaj Aug. 2013)
      IF (LUMOLDEN_MOS .GT. 0) CALL GPCLOSE(LUMOLDEN_MOS,'DELETE')
      CALL GPOPEN(LUMOLDEN_MOS,'molden_MOS.tmp','UNKNOWN',' ',
     &            'FORMATTED',IDUMMY,.FALSE.)
      REWIND LUMOLDEN_MOS

      RETURN
      END
      SUBROUTINE MOLDEN_TAIL
#include "implicit.h"
#include "dummy.h"
#include "maxorb.h"
#include "priunit.h"
#include "inftap.h"
#include "molde.h"
      CHARACTER*80 LINE
C
      IF (LUMOLDEN .LE. 0) THEN
         LUMOLDEN = -1
         CALL GPOPEN(LUMOLDEN,'molden.inp','UNKNOWN',' ',
     &               'FORMATTED',IDUMMY,.FALSE.)
      END IF

      IF (LUMOLDEN_MOS .LE. 0) THEN
         LUMOLDEN_MOS = -1
         CALL GPOPEN(LUMOLDEN_MOS,'molden_MOS.tmp','UNKNOWN',' ',
     &               'FORMATTED',IDUMMY,.FALSE.)
      END IF

!     skip to end of file
      DO WHILE (.TRUE.)
         READ (LUMOLDEN,'()',END=100)
      END DO
  100 CONTINUE
#ifdef VAR_MFDS
      BACKSPACE LUMOLDEN
#endif

      REWIND LUMOLDEN_MOS

!     append lines from LUMOLDEN_MOS to LUMOLDEN
      DO WHILE (.TRUE.)
         READ(LUMOLDEN_MOS, '(A)', END=300, ERR=200) LINE
         LEN_LINE = MAX(1,LNBLNK(LINE))
         WRITE(LUMOLDEN, '(A)') LINE(1:LEN_LINE)
      END DO
  200 CONTINUE
         WRITE(LUPRI,'(//A)')
     &   ' WARNING: error while reading molden_MOS.tmp file,'//
     &   ' molden.inp file may not be complete.'
         GOTO 9000
  300 CONTINUE

 9000 CONTINUE
      WRITE(LUMOLDEN,'(A)') '[End of Molden output from Dalton]'
      CALL GPCLOSE(LUMOLDEN    ,'KEEP')
      CALL GPCLOSE(LUMOLDEN_MOS,'KEEP')
      RETURN
      END
      SUBROUTINE MOLDEN_GTO(NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &                 KATOM,KANG,KBLOCK,KPRIM,CPRIMU,NRMPRI)

#include "implicit.h"
#include "maxorb.h"
#include "molde.h"
#include "inftap.h"
#include "priunit.h"
#include "mxcent.h"
#include "aovec.h"
#include "aosotr.h"
#include "maxaqn.h"
#include "codata.h"
      PARAMETER (D0 = 0.0D0)
C
#include "ccom.h"
#include "cbirea.h"
#include "nuclei.h"
#include "primit.h"


      CHARACTER CHRSEG*10, GTOTMP*1, SPDCAR*1
      LOGICAL SEG, NRMPRI
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK),
     &          SEG(KBLOCK)
      IF (.NOT. DONEIT) THEN

         WRITE(LUMOLDEN,'(A)') '[GTO]'

         ICENT  = 0
c
cPRT   icento is the offset counter for all centres;
cPRT   icent  counts only symmetry-distinct centres
c
         ICENTO = 1
         IBLOCK = 0
         INON   = 0
         IPROLD = 0
         IPRIM = 0
         IBS = 0
         DO 100 I = 1, NONTYP
            DO 110 N = 1, NONT(I)
               ICENT = ICENT + 1
               NDEG  = NUCDEG(ICENT)
               ILL   = 0
               DO, L = 1, NDEG
                  WRITE (LUMOLDEN,'(I5,A2)') ICENTO+L-1,' 0'
                  KBCH  = IBLOCK
                  IF (L.EQ.1) THEN
                     IPROLD  = IPRIM
                  ELSE
                     IPRIM = IPROLD
                  END IF
                  DO 200 J = 1, IQM(I) ! l quantum number s, p, d, ...
                     GTOTMP = SPDCAR(J-1)
                     DO 200 K = 1, JCO(J,I) ! n of blocks for a given AO
                        KBCH = KBCH + 1
                        NNUC  = NUC(KBCH)
                        NNRC  = NRC(KBCH)
                        IF (NNUC .EQ. 0) GO TO 200
                        ITYP = NHKOFF(J)
                        IPSTRT = IPRIM
                        IPRIM = IPRIM + NNUC
                        ITYP = ITYP + 1
                        DO 420 INNRC = 1, NNRC

                           WRITE (LUMOLDEN,'(1X,A1,1X,I3,A5)')
     &                          GTOTMP, NNUC, ' 1.00'
                           DO 410 INNUC = 1, NNUC
                              IF (PRIEXP(IPSTRT+INNUC) .GT. 1.D6) THEN
                                 WRITE (LUMOLDEN,'(1X,F15.2,1X,F15.10)')
     &                             PRIEXP(IPSTRT+INNUC),
     &                             CPRIMU(INNUC,INNRC,KBCH)
                              ELSE
                                 WRITE (LUMOLDEN,'(1X,F15.7,1X,F15.10)')
     &                             PRIEXP(IPSTRT+INNUC),
     &                             CPRIMU(INNUC,INNRC,KBCH)
                              END IF

 410                       CONTINUE
 420                    CONTINUE
 200                 CONTINUE
                     WRITE(LUMOLDEN,'()')
                  END DO ! L = 1, NDEG
                  ICENTO = ICENTO + NDEG
 110           CONTINUE
            IBLOCK = IBLOCK + NBLCK(I)
 100     CONTINUE

      ENDIF

      DONEIT=.TRUE.
      CALL FLSHFO(LUMOLDEN)
      RETURN
      END

      SUBROUTINE MOLDEN_ATOMS(WORD,LU)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "codata.h"

      CHARACTER*6 NAME
      CHARACTER*4 WORD
      LOGICAL     ECP_is_used

C     include character ASYMB(103)*2 with atomic symbols
#include "asymb.h"

#include "molde.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
#include "cbirea.h"
#include "inftap.h"
#include "chrxyz.h"
#include "chrsgn.h"
#include "chrnos.h"


!hjaaj Aug 2013: DONEIU obsolete, instead we use temporary LUMOLDEN_MOS as LU for 'ATOM'
!hjaaj IF (.NOT. DONEIU .OR. WORD .EQ. 'XYZ' .OR. WORD .EQ. 'FREQ') THEN

         NCOOR = 3*NUCDEP
         IF (WORD .EQ. 'XYZ ') WRITE(LU,'(I5/)') NUCDEP
         IF (WORD .EQ. 'ATOM') WRITE(LU,'(A)') '[Atoms] AU'
         IF (WORD .EQ. 'FREQ') WRITE(LU,'(/A)')'[FR-COORD]'

         ECP_is_used = .FALSE.
         IATOM = 0
         DO 100 ICENT = 1, NUCIND
            ICHARGE = IZATOM(ICENT)
            IF (ICHARGE .GT. 0) THEN ! only test for ECP on real nuclei
               IF (ICHARGE .NE. NINT(CHARGE(ICENT)))
     &            ECP_is_used = .TRUE.
            END IF
            MULCNT = ISTBNU(ICENT)
            NAME   = '      '
            J = 0
            DO I = 1,4
               IF (NAMN(ICENT)(I:I) .NE. ' ') THEN
                  J = J + 1
                  NAME(J:J) = NAMN(ICENT)(I:I)
               END IF
            END DO
            IF (MULT(MULCNT) .EQ. 1) THEN
               IF (WORD .EQ. 'ATOM') THEN
C FIXME ?          hjaaj Oct 2003: should we remove point charges ??
                  IATOM = IATOM + 1
                  WRITE (LU,'(A,1X,I5,1X,I5,3(1X,F20.10))')
     &                   NAME,IATOM,ICHARGE,(CORD(K,ICENT),K=1,3)
               ELSE IF (WORD .EQ. 'FREQ') THEN
                  WRITE (LU,'(A,3(1X,F20.10))')
     &                   ASYMB(ICHARGE),(CORD(K,ICENT),K=1,3)
               ELSE IF (WORD .EQ. 'XYZ ') THEN
                  WRITE (LU,'(A,3(1X,F20.10))')
     &                   NAME,(XTANG*CORD(K,ICENT),K=1,3)
               END IF
            ELSE
               J = J + 1
               NAME(J:J) = '_'
               J = J + 1
               JATOM = 0
               DO 200 ISYMOP = 0, MAXOPR
                  IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
                     JATOM = JATOM + 1
                     NAME(J:J) = CHRNOS(JATOM)
                     CRX = PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
                     CRY = PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
                     CRZ = PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
                     IF (WORD .EQ. 'ATOM') THEN
                        IATOM = IATOM + 1
                        WRITE (LU,'(A,1X,I5,1X,I5,3(1X,F20.10))')
     &                        NAME,IATOM,ICHARGE,CRX,CRY,CRZ
                     ELSE IF (WORD .EQ. 'FREQ') THEN
                        WRITE (LU,'(A,3(1X,F20.10))')
     &                        ASYMB(ICHARGE),CRX,CRY,CRZ
                     ELSE IF (WORD .EQ. 'XYZ ') THEN
                        WRITE (LU,'(A,3(1X,F20.10))')
     &                        NAME,XTANG*CRX,XTANG*CRY,XTANG*CRZ
                     END IF
                  END IF
 200           CONTINUE
            END IF
 100     CONTINUE

         IF (ECP_is_used .AND. WORD .EQ. 'ATOM') THEN
            WRITE(LU,'(A)') '[PSEUDO]'
            IATOM = 0
            DO, ICENT = 1, NUCIND
               ICHARGE = NINT(CHARGE(ICENT))
               MULCNT = ISTBNU(ICENT)
               NAME   = '      '
               J = 0
               DO I = 1,4
                  IF (NAMN(ICENT)(I:I) .NE. ' ') THEN
                     J = J + 1
                     NAME(J:J) = NAMN(ICENT)(I:I)
                  END IF
               END DO
               IF (MULT(MULCNT) .EQ. 1) THEN
                  IATOM = IATOM + 1
                  WRITE (LU,'(A,1X,I5,1X,I5)') NAME,IATOM,ICHARGE
               ELSE
                  J = J + 1
                  NAME(J:J) = '_'
                  J = J + 1
                  JATOM = 0
                  DO, ISYMOP = 0, MAXOPR
                  IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
                     JATOM = JATOM + 1
                     NAME(J:J) = CHRNOS(JATOM)
                     CRX = PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
                     CRY = PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
                     CRZ = PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
                     IATOM = IATOM + 1
                     WRITE (LU,'(A,1X,I5,1X,I5)') NAME,IATOM,ICHARGE
                  END IF
                  END DO ! ISYMOP
               END IF
            END DO ! ICENT
         END IF

!hjaajEND IF
!hjaajDONEIU = .TRUE.
      CALL FLSHFO(LU)
      RETURN
      END

      SUBROUTINE MOLDEN_MOS(ITASK,CMO,OCCUP,ORBTRA,UCMO,ORBVEC)
C
C CMO    =  contains MO coefficients, symmetry packed
C UCMO   =  contains MO coefficients, symmetry unpacked
C itask  =  1 : print everything to temporary Molden file
C           2 : save orbital energies in OREN in molde.h

#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxaqn.h"

! gnrinf.h: WFTYPE
#include "gnrinf.h"
! inftap.h: LUMOLDEN_MOS
#include "inftap.h"
! aosotr.h: JAOAO(:,2)
#include "aosotr.h"
! molde.h: DONIV, OREN(:)
#include "molde.h"
! inforb.h: NBAST, NORBT, ...
#include "inforb.h"
! pgroup.h: REP
#include "pgroup.h"
! ccom.h: SPHNRM
#include "ccom.h"

      LOGICAL     WRTELEM,WRTZERO
      CHARACTER*8 MDATE, MTIME
      DIMENSION   CMO(*), OCCUP(*), ORBTRA(*), UCMO(NBAST,NORBT)
      DIMENSION   ORBVEC(*)

#include "chrsgn.h"

      IF (ITASK .EQ. 1) THEN

         REWIND LUMOLDEN_MOS

         CALL GETDAT(MDATE,MTIME)
         WRITE(LUMOLDEN_MOS,'(A/2A,2X,A,2X,A)')
     &      '[TITLE]',
     &      '*** Dalton interface to Molden, wave function type : ',
     &         WFTYPE,MDATE,MTIME

         CALL MOLDEN_ATOMS('ATOM',LUMOLDEN_MOS) ! atomic coordinates are OK, also for Cartesian GTOs ;-)

         IF (.NOT. DONEIV) THEN
            IF (SPHNRM) THEN
               WRITE(LUMOLDEN_MOS,'(A)') '[5D7F]' ! 5d, 7f
               WRITE(LUMOLDEN_MOS,'(A)') '[9G]'   ! 9g
            ELSE
               NINFO = NINFO + 1
               WRITE (LUPRI,'(//A)') ' INFO: Sorry, plot of MOs with '//
     &         'Molden is only implemented for spherical GTOs'
               DONEIV = .TRUE.
            END IF
         END IF

         IF (.NOT. DONEIV) THEN
            WRITE(LUMOLDEN_MOS,'(A)') '[MO]'
C           call traorb to generate SO to AO transformation
            MOPRIN = 0
            CALL TRAORB(ORBTRA,NBAST,-1,MOPRIN)
            CALL DZERO(UCMO,  NBAST * NORBT)
            CALL DZERO(ORBVEC,NBAST)
            CALL UPKCMO(CMO,UCMO)
            jORB = 0
            DO ISYM = 1,NSYM
               NOCCI = NOCC(ISYM)
               DO I = 1,NORB(ISYM)
                  jORB = jORB + 1
                  II = IORB(ISYM) + I
                  WRITE(LUMOLDEN_MOS,'(2A)')     'Sym= ',REP(ISYM-1)
                  WRITE(LUMOLDEN_MOS,'(A,F9.4)') 'Ene= ',OREN(II)
                  WRITE(LUMOLDEN_MOS,'(A)')      'Spin= Alpha'
                  WRITE(LUMOLDEN_MOS,'(A,F7.4)') 'Occup= ',OCCUP(jORB)
                  CALL DGEMV('T',NBAST,NBAST,1.0D0,ORBTRA,NBAST,
     &                 UCMO(1,jORB),1,0.0D0,ORBVEC,1)
                  DO IB = 1,NBAST
                     IF (ABS(ORBVEC(JAOAO(IB,2))) .GE. 5.0D-7)  ! only write out if non-zero in F15.6 format
     &                  WRITE(LUMOLDEN_MOS,'(I5,1X,2F15.6)')
     &                    IB,ORBVEC(JAOAO(IB,2))
                  END DO
               END DO
            END DO
         END IF

         REWIND LUMOLDEN_MOS

      END IF

      IF (ITASK .EQ. 2) THEN
         DO I = 1, NORBT
            OREN(I)=CMO(I)
         ENDDO
      END IF

      RETURN
      END

      SUBROUTINE MOLDEN_FREQ(EVEC,NUMMOD,NCORD,FREQAU)
#include "implicit.h"
#include "maxorb.h"
#include "codata.h"
#include "molde.h"
#include "inftap.h"

      DIMENSION EVEC(NCORD,NCORD), FREQAU(NCORD)
      SXFAMU = SQRT(XFAMU)
      WRITE(LUMOLDEN,'(A)') '[FREQ]'

      DO IMODE = 1, NUMMOD
         WRITE(LUMOLDEN,'(F10.2)') FREQAU(IMODE)*XTKAYS
      END DO

      CALL MOLDEN_ATOMS('FREQ',LUMOLDEN)

      WRITE(LUMOLDEN,'(/A)') '[FR-NORM-COORD]'
      DO IMODE = 1, NUMMOD
         WRITE(LUMOLDEN,'(A,I10)') 'Vibration ',IMODE
         WRITE(LUMOLDEN,'(3(1X,F20.8))')
     &    (SXFAMU*EVEC(I,IMODE),I=1,NCORD)
      END DO
      CALL FLSHFO(LUMOLDEN)
      END

      SUBROUTINE MOLDEN_SCFCON(ITER,EMCSCF,WRITENOW)
#include "implicit.h"
#include "maxorb.h"
#include "codata.h"
#include "molde.h"
#include "inftap.h"
        LOGICAL WRITENOW
        IF(WRITENOW)THEN
           IF (.NOT. DONEIW) THEN

              WRITE(LUMOLDEN,'(A)') '[SCFCONV]'
              WRITE(LUMOLDEN,'(A,I3)') 'scf-first  1  THROUGH ',ITER
              WRITE(LUMOLDEN,'(F20.10)') (OROC(I),I=1,ITER)
              CALL FLSHFO(LUMOLDEN)

           END IF
           DONEIW = .TRUE.

        ELSE
           OROC(ITER)=EMCSCF
        ENDIF

        RETURN
      END

      SUBROUTINE MOLDEN_GECON(WRITENOW,EMCSCF)
#include "implicit.h"
#include "maxorb.h"
#include "codata.h"
#include "molde.h"
#include "inftap.h"
      LOGICAL WRITENOW
      IF(WRITENOW)THEN

         WRITE(LUMOLDEN,'(A)') '[GEOCONV]'
         WRITE(LUMOLDEN,'(A)') 'energy'

         DO 100 I=1,IXYZPRINTS
            WRITE(LUMOLDEN,'(F20.10)') EMCEN(I)
 100     CONTINUE

      ELSE
         IF (.NOT. DONEIX) THEN
            WRITE(LUMOLDEN,'(A)') '[GEOMETRIES] XYZ'
            DONEIX = .TRUE.
         END IF

         CALL MOLDEN_ATOMS('XYZ ',LUMOLDEN)
         IXYZPRINTS = IXYZPRINTS + 1
         EMCEN(IXYZPRINTS) = EMCSCF
      ENDIF
      CALL FLSHFO(LUMOLDEN)
      RETURN
! -- end of abacus/abamolden.F --
      END
