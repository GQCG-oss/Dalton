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
C FILE: sirius/drcctl.F
C
#ifdef old_log
===========================================================================
900108-hjaaj
  LBINTD = optimal buffer length .le. 1024 (instead of N2ORBT)
  new LUINTD label DRCINFO, after which LBINTD,LVLDRC may be read
891227-hjaaj
  DRCCTL : new parameter LVLDRC,
           = 0 then active-active distrib.
           = 1 all occupied-occupied distributions
           else all distributions
891127-hjaaj
  DRCCTL : NCDDIS = MIN(NCDDIS,NNORBX)
  %s/IS(/ISMO(/
  Use IPRTRA instead of IPRI6
  Print number of passes over LUINTM (= NLOAD)
===========================================================================
#endif
C  /* Deck drcctl */
      SUBROUTINE DRCCTL(LVLDRC,CMO,WRK,LWRK)
C
C Last revision Dec 89, Oct 2003 (CMO check) hjaaj
C
C PURPOSE:
C  SET UP MODRCINT CONTAINING DIRAC INTEGTRAL DISTRIBUTIONS
C  <**/CD> C.ge.D
C
C If (LVLDRC .eq. 0)      then C,D both active
C else if (LVLDRC .eq. 1) then C,D both occupied
C else                    C,D general
C
#include "implicit.h"
#include "priunit.h"
C
      REAL*8  CMO(*), WRK(LWRK)
      LOGICAL DRCOLD
C
#include "iratdef.h"
#include "lbmxsq.h"
C
C Used from common blocks:
C   INFORB : N2ORBT, N2ORBX, ...
C   INFTAP : LUINTD, LBINTD
C   INFTRA : IPRTRA
C
#include "inforb.h"
#include "inftap.h"
#include "inftra.h"
#include "infpri.h"
C
      CALL QENTER('DRCCTL')
C
C     Open MODRCINT and check if requested integrals already are there ...
C
      LUINTD = -1
      CALL GPOPEN(LUINTD,'MODRCINT','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      IF ( DRCOLD(CMO,LUINTD,LVLDRC,WRK,LWRK) ) GO TO 9000
C
C Set LBINTD to optimal buffer length .le. LBMXSQ for LUINTD
C ( N2ORBT is the maximum number of non-zero elements in
C any distribution).
C
      NBINTD = (N2ORBT-1)/LBMXSQ + 1
      LBINTD = (N2ORBT-1)/NBINTD + 1
C
C
C DETERMINE NUMBER OF DIRAC DISTRIBUTIONS WHICH CAN BE
C KEPT IN CORE IN ONE LOAD
C
      KOUTB  = 1
      KH2    = KOUTB  + LBINTD + ( LBINTD + 2 - 1)/ 2 + 1 ! integer*4 for IOUTB
      KWRK1  = KH2    + N2ORBX
      LWRK1  = LWRK   - KWRK1
C
C     LNDACD = how much work space needed in DRCACD (for NXTH2M)
      IF (NEWTRA) THEN
         LNDACD = (3*(36*36+1) + NNORBX + 1)/IRAT + NNORBX + 2*N2ORBT
      ELSE
         LNDACD = ((1+IRAT)*LBINTM+2)/IRAT + 1
      END IF
      NCDDIS = (LWRK1-LNDACD)/N2ORBX
      IF ( NCDDIS.LE.0 ) THEN
         WRITE(LUPRI,'(/A,I5)') 'DRCCTL ERROR: '//
     *   ' NOT ENOUGH MEMORY TO HAVE ONE DIRAC DISTRIBUTION. NCDDIS='
     *   ,NCDDIS
         CALL QTRACE(LUERR)
         CALL QUIT('DRCCTL: NOT ENOUGH MEMORY FOR ONE DISTRIBUTION')
      END IF
      NCDDIS = MIN(NCDDIS,NNORBX)
C
      KH2ACD = KWRK1
      KWRK1  = KH2ACD + NCDDIS*N2ORBX
      LWRK1  = LWRK   - KWRK1
C
C INITIALIZE VARIABLES
C
      NLOAD  = 0
      LABCD  = 0
      ICMIN  = 1
      IDMIN  = 1
      IF (LVLDRC .EQ. 0) THEN
         NTOTDI = NNASHX
      ELSE IF (LVLDRC .EQ. 1) THEN
         NTOTDI = NNOCCX
      ELSE
         NTOTDI = NNORBX
      END IF
      ILOWDI = 1
C
C REPEAT UNTILL
C
 100  CONTINUE
         NLOAD  = NLOAD + 1
         IHGHDI = MIN(NTOTDI,NCDDIS+ILOWDI-1)
         ICMAX  = INT(0.49999999D0+SQRT(0.25D0 +2.0D0*IHGHDI))
C
C Retrieve this load of DIRAC INTEGRAL DISTRIBUTIONS
C
         CALL DRCACD(WRK(KH2ACD),WRK(KH2),ICMIN,ICMAX,ILOWDI,IHGHDI,
     *               LVLDRC,WRK(KWRK1),LWRK1)
C
C WRITE these DIRAC DISTRIBUTIONS OUT ON LUINTD
C
         CALL DRCWRT(WRK(KH2ACD),ICMIN,IDMIN,ILOWDI,IHGHDI,NTOTDI,
     *               LVLDRC,LABCD,WRK(KOUTB),WRK(KOUTB),LUINTD,CMO)
C        CALL DRCWRT(H2ACD,ICMIN,IDMIN,ILOWDI,IHGHDI,NTOTDI,
C    *               LVLDRC,LABCD,OUTB,IOUTB,LUINTD,CMO)
C ICMIN AND IDMIN ARE UPDATED IN DRCWRT
         ILOWDI = IHGHDI + 1
      IF (ILOWDI.LE.NTOTDI) GO TO 100
C
C END REPEAT
C
      IF (IPRTRA.GT.0) WRITE(LUPRI,1000) LABCD,NLOAD
 1000 FORMAT(/' DRCCTL: Number of Dirac integrals written on file',
     &    ' "MODRCINT" is',I12,
     &   /'         Number of passes over Mulliken file MOTWOINT',I5)
C
 9000 CONTINUE
      CALL GPCLOSE(LUINTD,'KEEP')
      CALL QEXIT('DRCCTL')
      RETURN
C     ... end of DRCCTL
      END
C  /* Deck drcacd */
      SUBROUTINE DRCACD(H2ACD,H2,ICLOW,ICHGH,ILOWDI,IHGHDI,
     &                  LVLDRC,WRK,LWRK)
C
C Purpose:
C   SORT MULLIKEN INTEGRALS TO GET desired (XY) DIRAC DISTRIBUTIONS
C
C If (LVLDRC .eq. 0)      then X,Y both active
C else if (LVLDRC .eq. 1) then X,Y both occupied
C else                    X,Y general
C
C
C    MULLIKEN INTEGRALS (CD,AB) = (XP,YQ)
C    TO DIRAC INTEGRALS <CA,DB> = <XY,PQ>
C
#include "implicit.h"
      DIMENSION H2ACD(NORBT,NORBT,*),H2(NORBT,NORBT),WRK(LWRK)
C
C Used from common blocks:
C   INFORB : NSYM,NORB(*),IORB(*),?
C   INFIND : ISMO(*),?
C   INFTRA : IPRTRA, THRP
C
#include "maxash.h"
#include "maxorb.h"
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
#include "inftra.h"
C
      DIMENSION NEEDMU(-4:6)
C
      CALL QENTER('DRCACD')
      IF ( IPRTRA.GT.20 ) THEN
         WRITE(LUPRI,'(/4(A,I7))')'  ICLOW :',ICLOW, '  ICHGH :',ICHGH,
     *                            '  ILOWDI:',ILOWDI,'  IHGHDI:',IHGHDI
      ENDIF
      NEEDMU(-4:6) = 0
      NEEDMU(2) = 1 ! ia
      NEEDMU(3) = 1 ! aa
      NEEDMU(5) = 1 ! as
      IF (LVLDRC .GE. 1) THEN
         NEEDMU(1) = 1 ! ii
         NEEDMU(4) = 1 ! is
      END IF
      IF (LVLDRC .EQ. 2) NEEDMU(6) = 1 ! ss
C
C OFFSET FOR (CD) DISTRIBUTIONS
C
      ILOWOF = ILOWDI - 1
      NDIST  = IHGHDI - ILOWOF
C
C     --- zero H2ACD
C
      LENGTH = N2ORBX * (IHGHDI - ILOWDI + 1)
      CALL DZERO(H2ACD,LENGTH)
C
C     --- prepare for reading MO integrals ...
C
      KFREE  = 1
      LFREE  = LWRK
      IDIST  = 0
 200  CALL NXTH2M(IC,ID,H2,NEEDMU,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 9500
         IF ( IPRTRA.GT.150 ) THEN
           WRITE(LUPRI,'(/A,I8)')' DRCACD: Mulliken distribution',IDIST
           WRITE(LUPRI,'(A,2I5)')' IC, ID =',IC,ID
           CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         ENDIF
         ICDSYM = MULD2H(ISMO(IC),ISMO(ID))
         NCW = ISW(IC)
         NDW = ISW(ID)
         IF (LVLDRC .EQ. 0) THEN
            NCW = NCW - NISHT
            NDW = NDW - NISHT
         END IF
C
C        See if we need this distribution:
C         ICDTYP = 0 we do not need this distribution
C                    DISTRIBUTION HAS NO INDICES WITHIN
C                    DESIRED DISTRIBUTION RANGE
C         ICDTYP = 1 DISTRIBUTION WITH ONE INDEX WITHIN
C                    DESIRED DISTRIBUTION RANGE
C         ICDTYP = 2 DISTRIBUTION WITH TWO INDICES WITHIN
C                    DESIRED DISTRIBUTION RANGE
C
         ICDTYP = 0
         IF ( NCW .GE. ICLOW .AND. NCW .LE. ICHGH ) ICDTYP = 1
         IF ( NDW .GE. ICLOW .AND. NDW .LE. ICHGH .AND.
     &        NDW .NE. NCW  ) THEN
            IF (ICDTYP .EQ. 1) THEN
               ICDTYP = 2
               IDOFF  = IROW(NDW) - ILOWOF
            ELSE
               ICDTYP = 1
               ISWAP  = NCW
               NCW    = NDW
               NDW    = ISWAP
               ISWAP  = IC
               IC     = ID
               ID     = ISWAP
            END IF
         END IF

         IF (ICDTYP.EQ.0)     GO TO 200

         ICOFF  = IROW(NCW) - ILOWOF
C
C
C      ... C IS ACTIVE AND WITHIN DESIRED RANGE
C      ... D IS ACTIVE AND WITHIN DESIRED RANGE IF ICDTYP = 2
C
         DO 280 IA = 1,NORBT
            IASYM = ISMO(IA)
            IBSYM = MULD2H(IASYM,ICDSYM)
            IBST  = IORB(IBSYM) + 1
            IBEND = MIN(IA,IORB(IBSYM)+NORB(IBSYM))
            DO 290 IB = IBST,IBEND
               H2AB = H2(IB,IA)
            IF (ABS(H2AB) .LE. THRP) GO TO 290
               NAW  = ISW(IA)
               NBW  = ISW(IB)
               IF (LVLDRC .EQ. 0) THEN
                  NAW = NAW - NISHT
                  NBW = NBW - NISHT
               END IF
               IF ( (NAW.GE.1).AND.(NAW.LE.NCW) ) THEN
                  NCA = ICOFF + NAW
                  IF ( (NCA.GT.0).AND.(NCA.LE.NDIST ) ) THEN
                      H2ACD(ID,IB,NCA) = H2AB
                  END IF
               END IF
               IF ( (NBW.GE.1).AND.(NBW.LE.NCW) ) THEN
                  NCB = ICOFF + NBW
                  IF ( (NCB.GT.0).AND.(NCB.LE.NDIST) ) THEN
                      H2ACD(ID,IA,NCB) = H2AB
                  END IF
               END IF
               IF (ICDTYP.EQ.2) THEN
                  IF ( (NAW.GE.1).AND.(NAW.LE.NDW) ) THEN
                     NDA = IDOFF + NAW
                     IF ( (NDA.GT.0).AND.(NDA.LE.NDIST) ) THEN
                        H2ACD(IC,IB,NDA) = H2AB
                     END IF
                  END IF
                  IF ( (NBW.GE.1).AND.(NBW.LE.NDW) ) THEN
                     NDB = IDOFF + NBW
                     IF ( (NDB.GT.0).AND.(NDB.LE.NDIST) ) THEN
                        H2ACD(IC,IA,NDB) = H2AB
                     END IF
                  END IF
               END IF
  290       CONTINUE
  280    CONTINUE

      GO TO 200
C
C ***
C
 9500 CONTINUE
      IF ( IPRTRA.GT.150 ) THEN
         DO 700 IWR = ILOWDI,IHGHDI,1
            WRITE(LUPRI,'(/A,I8)')' Dirac DISTRIBUTION NUMBER: ',IWR
            CALL OUTPUT(H2ACD(1,1,IWR-ILOWOF),1,NORBT,1,NORBT,NORBT,
     *                  NORBT,1,LUPRI)
 700     CONTINUE
      ENDIF
C
C *** end of subroutine DRCACD
C
      CALL QEXIT('DRCACD')
      RETURN
      END
C  /* Deck drcwrt */
      SUBROUTINE DRCWRT(H2ACD,ICMIN,IDMIN,ILOWDI,IHGHDI,NTOTDI,
     *                  LVLDRC,LABCD,OUTB,IOUTB,LUINTD,CMO)
C
C WRITE DIRAC DISTRIBUTIONS ON FILE LUINTD
C
C If (LVLDRC .eq. 0)      then X,Y both active
C else if (LVLDRC .eq. 1) then X,Y both occupied
C else                    X,Y general
C
C
C NOTE: EQUIVALENCE (OUTB,IOUTB) in call of DRCWRT
C
#include "implicit.h"
C
      REAL*8    H2ACD(NORBT,NORBT,*),OUTB(*),CMO(NCMOT)
      INTEGER*4 IOUTB(*)
#include "iratdef.h"
      PARAMETER ( D0 = 0.0D0 )
C
C Used from common blocks:
C   INFTRA : IPRTRA,THRP
C
#include "maxash.h"
#include "maxorb.h"
#include "priunit.h"
#include "inforb.h"
#include "infind.h"
#include "inftap.h"
#include "inftra.h"
C
      CHARACTER*8 LAB123(3), LABELD(2)
C
      DATA LAB123/'********','********','********'/
      DATA LABELD/'DRCINFO ','DRCTWOEL'/
C
C ************** LENGTH OF BUFFER FOR H2ACD ON LUINTD
C
      LOUT  =LBINTD
      LOUTI =2*LOUT ! "2*" because always INTEGER*4 IOUTB(:)
      LOUT2 =LOUTI+LOUT
      LOUT21=LOUT2+1
      LOUT22=LOUT2+2
C
C WRITE LABEL ON LUINTD
C
      IF ( ICMIN.EQ.1 .AND. IDMIN.EQ.1 ) THEN
         CALL GETDAT(LAB123(2),LAB123(3))
         REWIND LUINTD
         WRITE(LUINTD)LAB123,LABELD(1)
         WRITE(LUINTD)LBINTD,LVLDRC,NCMOT,D0,D0,D0,D0
         WRITE(LUINTD)CMO
         WRITE(LUINTD)LAB123,LABELD(2)
      END IF
C
C INITIALIZE VARIABLES
C
      IOUT = 0
      IDISOF = ILOWDI - 1
      NCW = ICMIN
      NDW = IDMIN
      ICDDIS = IROW(NCW) + NDW
      IF (ICDDIS.NE.ILOWDI) THEN
         WRITE(LUPRI,'(/A,3(A,I6))')
     *   ' *** DRCWRT: INCORRECT INTEGRAL DISTRIBUTION INTERVAL',
     *   ' NCW:',NCW,' NDW:',NDW,' ICDDIS:',ICDDIS
         CALL QTRACE(LUPRI)
         CALL QUIT(' DRCWRT: INCORRECT DISTRIBUTION INTERVAL')
      END IF
      IF ( IPRTRA.GT.20 ) THEN
         WRITE(LUPRI,'(/A,/A,2I8)')' ****DRCWRT*****',
     *   ' ICMIN,IDMIN',ICMIN,IDMIN
      END IF
C
C   ...  REPEAT UNTILL
C
 100  CONTINUE
         IF (LVLDRC .EQ. 0) THEN
            IC = ISX( NISHT + NCW )
            ID = ISX( NISHT + NDW )
         ELSE
            IC = ISX( NCW )
            ID = ISX( NDW )
         END IF
         INDCD = IC*2**16  + ID
         IOUTB(LOUT22) = INDCD
         IF (IPRTRA.GT.20) THEN
            WRITE(LUPRI,'(A,2I10)') 'Writing IC,ID',IC,ID
         END IF
C
C     ***** ALL INTEGRALS <AB/CD> FOR A GIVEN PAIR CD HAVE *****
C     ***** BEEN CREATED. WRITE THEM ON UNIT LUINTD WITH   *****
C     ***** INDICES C,D,A,B                                *****
         DO 170 IA=1,NORBT
            DO 160 IB = 1,NORBT
               P = H2ACD(IA,IB,ICDDIS-IDISOF)
               IF (ABS(P).LT.THRP) GO TO 160
               IALAST=IA
               IBLAST=IB
               LABCD=LABCD+1
               IOUT =IOUT+1
               IF(IOUT.GT.LOUT) THEN
                  IOUT=1
                  IOUTB(LOUT21)=LOUT
                  IF ( IPRTRA.GT.30 ) THEN
                     WRITE(LUPRI,'(/A,/A/,A,4I6,F20.14)')
     *                 ' SEVERAL RECORDS FOR DISTRIBUTION IC,ID',
     *                 ' LAST ELEMENT NOT WRITTEN OUT ',
     *                 ' IC,ID,IA,IB,P',IC,ID,IA,IB,P
                  END IF
C           ***** WRITE THIS BUFFER *****
                  CALL WRITI4(LUINTD,LOUT22,IOUTB)
               END IF
               OUTB(IOUT) = P
               IOUTB(LOUTI+IOUT) = IA*2**16 + IB
  160       CONTINUE
  170    CONTINUE
C   *** Going to next (CD), empty this buffer ***
         IF (IOUT .GT. 0) THEN
            IOUTB(LOUT21) = IOUT
            IF ( IPRTRA.GT.30 ) THEN
               WRITE(LUPRI,'(/A,/A/,A,4I6,F20.14)')
     *            ' LAST RECORD FOR DISTRIBUTION IC,ID',
     *            ' LAST ELEMENT WRITTEN OUT ',
     *            ' IC,ID,IA,IB,P',IC,ID,IALAST,IBLAST,P
            END IF
            CALL WRITI4(LUINTD,LOUT22,IOUTB)
            IOUT = 0
         END IF
         NDW  = NDW + 1
         IF (NDW.GT.NCW) THEN
            NDW = 1
            NCW = NCW+1
         END IF
         ICDDIS = IROW(NCW) + NDW
      IF (ICDDIS.LE.IHGHDI) GO TO 100
C
C ... END REPEAT
C
C     ***** WRITE LAST BUFFER *****
C
      IF ( ICDDIS.GT.NTOTDI) THEN
         IOUTB(LOUT21)=-1
         IOUTB(LOUT22)=-1
         IF ( IPRTRA.GT.20 ) THEN
            WRITE(LUPRI,'(/A,I10)')
     *      ' -1 WRITTEN OUT ON LAST(EMPTY) RECORD, ICDDIS;',ICDDIS
         END IF
         CALL WRITI4(LUINTD,LOUT22,IOUTB)
      END IF
C
C UPDATE C AND D FOR NEXT LOAD
C
      ICMIN = NCW
      IDMIN = NDW
      RETURN
C *** END OF DRCWRT
      END
C  /* Deck drcctl */
      LOGICAL FUNCTION DRCOLD(CMO,LUINTD,LVLDRC,WRK,LWRK)
C
C Check if requested integrals in Dirac format are already available.
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION CMO(*), WRK(LWRK)
      LOGICAL FNDLAB
C
C thrzer.h : THRZER
C inforb.h : NCMOT
C inftra.h : IPRTRA
#include "thrzer.h"
#include "inforb.h"
#include "inftra.h"
C
      DRCOLD = .FALSE.
      REWIND LUINTD
      IF (.NOT.FNDLAB('DRCINFO ',LUINTD)) GOTO 9999
      READ (LUINTD) LBINTD, LVLDRC_OLD, NCMOT_OLD
      IF (LVLDRC_OLD .GE. LVLDRC) THEN
      IF (NCMOT_OLD .EQ. NCMOT) THEN
C        read CMO matrix from LUINTD
C        and subtract from input CMO matrix
         READ  (LUINTD) WRK(1:NCMOT)
         CALL DAXPY(NCMOT,-1.0D0,CMO,1,WRK,1)
         I = IDAMAX(NCMOT,WRK,1)
         DRCOLD = ABS(WRK(I)) .LE. THRZER
      END IF
      END IF
      IF (DRCOLD .AND. IPRTRA.GE.0) WRITE (LUPRI,'(/A/A)')
     &' DRCCTL abandoned: the required MO integrals in Dirac format',
     &' are already available on the MODRCINT file.'
 9999 RETURN
      END
C --- end of drcctl.F ---
