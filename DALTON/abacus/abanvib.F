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
C     /* Deck nvbinp */
      SUBROUTINE NVBINP(WORD)
C     ************************************************
C     **** Input routine for numerical averaging  ****
C     **** of properties. If nothing is already   ****
C     **** in module *NUMDRV, this will take care ****
C     **** finding maximum derivative order. No   ****
C     **** variables are initiated in this routine****
C     **** so one needs to run NMDINI and NMDINP  ****
C     **** first.                                 ****
C     ************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 8)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7,
     &          WORD1*7
#include "numder.h"
#include "prpndr.h"
#include "fcsym.h"
#include "cbinum.h"
#include "abainf.h"
#include "cbiwlk.h"
#include "cbirea.h"
#include "molinp.h"
C
      DATA TABLE /'.EFFECT','.HARM-P','.ANHA-P','.SPIN-S','.MODE A',
     &            '.ONLY-P','.CUBIC ','.P-BASI'/
C
      ICHANG = 0
      WORD1 = WORD
 100  CONTINUE
      READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         ICHANG = ICHANG + 1
         DO 200 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (1,2,3,4,5,6,7,8), I
            END IF
 200     CONTINUE
         IF (WORD .EQ. '.OPTION') THEN
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            GO TO 100
         END IF
         WRITE (LUPRI,'(/4A/)') ' Keyword "',WORD,
     *        '" not recognized in '//WORD1
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal keyword in '//WORD1)
 1       CONTINUE  ! .EFFECT
            NUMVIB = .TRUE.
            ANALZ1 = .TRUE.
            NMORDR = MAX(NMORDR,3-NAORDR)
            GOTO 100
 2       CONTINUE  ! .HARM-P
            NRMCRD = .TRUE.
            NPRPDR = .TRUE.
            FRSTNM = .TRUE.
            PRPVIB = .TRUE.
            ANALZ1 = .TRUE.
            NMRDRP = 2
            NARDRP = 0
            NMORDR = MAX(NMORDR,3-NAORDR) ! Should be 2-NAORDR ??
            FCLASS = 'C1 '
CRF 12/11 - 12 Added GOTO that seemed to have been left out
CRF            Also .HARM-P and .ANHA-P seems to do the same
CRF            only .HARM-P sets symmetry to C1.
CRF            This is somewhat counteri-intuitive and doesn't agree with manual
CRF            Problem seems to be that there's no branch around the anharmonic
CRF            averaging, for harmonic only
            GOTO 100
 3       CONTINUE  ! .ANHA-P
            NRMCRD = .TRUE.
            NPRPDR = .TRUE.
            FRSTNM = .TRUE.
            PRPVIB = .TRUE.
            ANALZ1 = .TRUE.
            NMRDRP = 2
            NARDRP = 0
            NMORDR = MAX(NMORDR,3-NAORDR)
            GOTO 100
 4       CONTINUE  ! .SPIN-S
            NSPNSP = .TRUE.
            GOTO 100
 5       CONTINUE  ! .MODE A
            MDEANA = .TRUE.
            GOTO 100
CRF 12/11 - 12 Test keywords for splitting force field from property derivatives
 6       CONTINUE  ! .ONLY-P
            PREHES = .TRUE.
            NRMCRD = .TRUE.
            NPRPDR = .TRUE.
            ANALZ1 = .TRUE.
            PRPVIB = .TRUE.
            FRSTNM = .TRUE.
            NMRDRP = 2
            NARDRP = 0
            NMORDR = 2 ! Seems like we'll have to keep this at 2 for now
            NAORDR = 0
            PRPONL = .TRUE.
            GOTO 100
 7       CONTINUE  ! .CUBIC
            REUHES = .TRUE. ! Set this here or not?
            NRMCRD = .TRUE.
            ANALZ1 = .TRUE.
            NMORDR = MAX(NMORDR,3-NAORDR)
            NUMVIB = .TRUE. ! Needed in order to print force field
            GOTO 100
 8       CONTINUE  ! .P-BASIS
CRF This only works if basis set library is used
            IF (NMLINE_basis .le. 0) THEN
               WRITE(LUPRI,'(/A/A)')
     &           '*** ERROR *** ' //
     &           'Keyword .P-BASIS only works when ' //
     &           'basis set library is used.',
     &           'The calculation can be split manually ' //
     &           'using the .CUBIC and .ONLY-P keywords'
               CALL QUIT('Input ERROR in '//WORD1)
            END IF
            IF (LMULBS) THEN
               WRITE(LUPRI,'(/A/A)')
     &           '*** ERROR *** ' //
     &           'Keyword .P-BASIS does not work ' //
     &           'for .R12AUX multiple basis sets.'
               CALL QUIT('ERROR: '//
     &         '.P-BASIS does not work for .R12AUX multiple basis sets')
            END IF
            REUHES = .TRUE.
            NRMCRD = .TRUE.
            ANALZ1 = .TRUE.
            NMORDR = MAX(NMORDR,3-NAORDR)
            NUMVIB = .TRUE. ! Needed in order to print force field
            PRPBAS = .TRUE.
            READ (LUCMD,*) PRPBTX
            GOTO 100
CRFend
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 300
      ELSE
         WRITE (LUPRI,'(/4A/)') ' Prompt "',WORD,
     *            '" not recognized in '//WORD1
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in '//WORD1)
      END IF
  300 CONTINUE
CRF We save a few variables in case of PRPBAS
      IF (PRPBAS) THEN
         NMRDBK = NMRDRP
         NARDBK = NARDRP
         NMRDRP = 0 
         NARDRP = 0
         PRPVIB = .FALSE.
         NPRPDR = .FALSE.
         FRSTNM = .FALSE. 
      END IF
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for *PROPAV section:',0)
         IF (.NOT. PRPONL) THEN
             WRITE (LUPRI,'(/5X,A,I4)') 'Numerical derivatives' //
     &             ' of energy calculated to order', NMORDR
             WRITE (LUPRI,'(5X,A,I4,A)') 'from analytical', NAORDR,
     &                              '. derivatives'
             WRITE (LUPRI,'(5X,A,I4)') 'Total order of differentiation:'
     &            , NAORDR+NMORDR
         ELSE 
             WRITE (LUPRI,'(5X,A)') 'Only property derivatives '//
     &             'calculated this run'
             WRITE (LUPRI,'(5X,A)') 'Force field must be availble '//
     &             'from a previous calculation'
         END IF
         IF (NUMVIB) THEN
            WRITE (LUPRI,'(5X,A)') 'Effective geometry calculated.'
         END IF
         IF (PRPVIB) THEN
            WRITE (LUPRI,'(5X,A)') 'Vibrational averages of' //
     &           ' properties calculated.'
            WRITE (LUPRI,'(/5X,A)') 'Vibrational averages of the ' //
     &           'following properties:'
CRF    The NSPNSP flag controls only a few headers like this
CRF    We should probably instead get this info from /ABAINF/: SPNSPN 
            IF (NSPNSP) THEN 
               WRITE (LUPRI,'(5X,A)') 'Spin-spin couplings.'
            END IF
CRF    This doen't seem to match, how the code currently works
C            WRITE (LUPRI,'(/5X,A)') 'Harmonic contribution calculated.'
C            IF ((NMORDR+NAORDR).GT.2) THEN
C               WRITE (LUPRI,'(5X,A)')
C     &              'Anharmonic contribution calculated.'
C            END IF
            WRITE (LUPRI,'(5X,A)') 'Harmonic and anharmonic '//
     &            'contributions calculated'
            WRITE (LUPRI,'(5x,A)') 'Cannot use group theory for this. '
            WRITE (LUPRI,'(5x,A)') 'Group used in calculations: ' //
     &                          FCLASS
         END IF
         IF (PRPBAS) THEN
            WRITE (LUPRI,'(5X,2A)') 'Basis set used in force field '//
     &            'calculation: ' , TRIM( MLINE(NMLINE_basis) ) 
            WRITE (LUPRI,'(5X,2A)') 'Basis set used in calculation '//
     &            'of property derivatives: ', TRIM( PRPBTX )
         END IF

         IF (MDEANA) THEN
            WRITE (LUPRI,'(5X,A)') 'An analysis of the vibrational ' //
     &            'modes will be performed.'
         END IF
      END IF
C
      RETURN
      END
C
C
C  /* Deck nvbini */
      SUBROUTINE NVBINI
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "prpndr.h"
#include "pvibav.h"
#include "abainf.h"
C
C     /PRPNDR/
      NSPNSP = .FALSE.
C
C     /ABAINF/
      SPNSPN = .FALSE.
C
C     /PVIBAV/
      CNMPRP = .FALSE.
      RETURN
      END
C
C
C  /* Deck nvbrst */
      SUBROUTINE NVBPIN(IPRINT)
C     ************************************************
C     *** Subroutine that resets a few variables   ***
C     *** if properties should be calculated using ***
C     *** normal coordinates.                      ***
C     ************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "prpndr.h"
#include "abainf.h"
#include "magone.h"
#include "spnout.h"
C
      IF (IPRINT.GT.20) THEN
         CALL HEADER('No properties calculated in this geometry of ' //
     &               'numerical differentiation',0)
         WRITE (LUPRI,'(A)') '                                    '
         WRITE (LUPRI,'(5X,A)') 'Properties neglected are: '
         IF (SPNSPN) THEN
            WRITE (LUPRI,'(5X,A)') 'No spin-spin couplings.'
         END IF
      END IF
C
      SPNSPN = .FALSE.
      SST    = .FALSE.
      DODSO  = .FALSE.
      DOPSO  = .FALSE.
      DOSD   = .FALSE.
      DOFC   = .FALSE.
C
      RETURN
      END
C
C
C  /* Deck nvbdrv */
      SUBROUTINE NVBDRV(DERIV,SYMCOR,FREQ,RNNORM,CSTART,WORK,LWORK,
     &                  NDERIV,IPRINT)
C     ***********************************************
C     *** Driver routine for vibrational analysis ***
C     *** from numerical derivatives. Property    ***
C     *** derivatives and frequencies.            ***
C     ***********************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "abainf.h"
#include "trkoor.h"
#include "cbinum.h"
#include "numder.h"
      DIMENSION DERIV(NDERIV), SYMCOR(NCOOR,NCOOR), FREQ(NCOOR),
     &          RNNORM(NCOOR), CSTART(NCOOR      ), WORK(LWORK)
C
CRF   Skip if only properties has been  calculated
      IF (NUMVIB .AND. (.NOT. PRPONL) ) THEN
         KTDRS = 1
         KEVEC = KTDRS + NCOOR**3
         KLAST = KEVEC + NCOOR**2
         LWRK = LWORK - KLAST
         CALL AHRMVB(DERIV,SYMCOR,FREQ,RNNORM,CSTART,WORK(KTDRS),
     &               WORK(KEVEC),WORK(KLAST),NDERIV,LWRK,IPRINT)
      END IF
C
      IF (PRPVIB) THEN
         NPPDER = NDCOOR
C
         IF ((NMRDRP.EQ.2).AND.PRPVIB) THEN
            NPPDER = NPPDER + NDCOOR
         ELSE IF (NMRDRP.GE.2) THEN
            NPPDER = NPPDER + NDCOOR*(NDCOOR+1)/2
         END IF
C
         IF (NMRDRP.GE.3) NPPDER = NPPDER
     &                           + NDCOOR*(NDCOOR+1)*(NDCOOR+2)/6
         IF (NMRDRP.GE.4) NPPDER = NPPDER
     &                      + NDCOOR*(NDCOOR+1)*(NDCOOR+2)*(NDCOOR+3)/24
C
         KLAST = 1
         IF (SPNSPN) THEN
            KSPNSP = KLAST
            KDSPSP = KSPNSP + 6*NCOOR**2
            KLAST  = KDSPSP + 6*NCOOR**2*NPPDER
         END IF
C
         LWRK1 = LWORK - KLAST
         IF (LWRK1.LT.1) CALL QUIT('Memory exceeded in PPVBAN')
         CALL PPVBAN(DERIV,FREQ,SYMCOR,RNNORM,WORK(KSPNSP),WORK(KDSPSP),
     &               WORK(KLAST),NPPDER,NDERIV,LWRK1,IPRINT)
      END IF
C
      RETURN
      END
C  /* Deck ahrmvb */
      SUBROUTINE AHRMVB(DERIV,SYMCOR,FREQ,RNNORM,CSTART,TDER,EVEC,WORK,
     &                  NDERIV,LWORK,IPRINT)
C     *********************************************************
C     *** Readying the third derivative matrix for analyses ***
C     *********************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
#include "mxcent.h"
      PARAMETER (D1 = 1.0D0)
#include "trkoor.h"
#include "cbinum.h"
#include "inftap.h"
#include "numder.h"
#include "cbiwlk.h"
      DIMENSION DERIV(NDERIV), SYMCOR(NCOOR,NCOOR), FREQ(NCOOR),
     &          TDER(NCOOR,NCOOR,NCOOR), EVEC(NCOOR,NCOOR),
     &          RNNORM(NCOOR), CSTART(NCOOR), WORK(LWORK)
C
C     *** Initialisation ***
C
      KDIM = NCOOR**3
      CALL DZERO(TDER,KDIM)
C
      DFAC = D1
      IDERIV = 0
      DO 100 ICOOR3 = 1, NCOOR
      DO 100 ICOOR2 = 1, ICOOR3
      DO 100 ICOOR1 = 1, ICOOR2
         IF (NRMCRD) DFAC = RNNORM(ICOOR1)*RNNORM(ICOOR2)*RNNORM(ICOOR3)
         IDERIV = IDERIV + 1
         TDER(ICOOR1,ICOOR2,ICOOR3) = DERIV(IDERIV)*DFAC
         TDER(ICOOR1,ICOOR3,ICOOR2) = DERIV(IDERIV)*DFAC
         TDER(ICOOR2,ICOOR1,ICOOR3) = DERIV(IDERIV)*DFAC
         TDER(ICOOR2,ICOOR3,ICOOR1) = DERIV(IDERIV)*DFAC
         TDER(ICOOR3,ICOOR1,ICOOR2) = DERIV(IDERIV)*DFAC
         TDER(ICOOR3,ICOOR2,ICOOR1) = DERIV(IDERIV)*DFAC
 100  CONTINUE
C
C     *** Transforming the third derivative from symmetry   ***
C     *** to cartesian coordinates (if not in normal coor). ***
C
      IF (.NOT.NRMCRD) THEN
         KSYMCO = 1
         KTMPTD = KSYMCO + NCOOR**2
         KLAST  = KTMPTD + NCOOR**3
         LWRK   = LWORK  - KLAST
         CALL TDRFSC(TDER,CSTART,WORK(KSYMCO),WORK(KTMPTD),WORK(KLAST),
     &               NCOOR,LWRK,IPRINT)
      END IF
C
C     *** Scaling the normal coordinates back. ***
C
      DO 200 ICOOR2 = 1, NDCOOR
      DO 200 ICOOR1 = 1, NCOOR
         EVEC(ICOOR1,ICOOR2) = RNNORM(ICOOR2)*SYMCOR(ICOOR1,ICOOR2)
 200  CONTINUE
C
C     *** Writing the third derivative to file ***
C
      CALL GPOPEN(LUWLK,ABAWLK,'NEW','SEQUENTIAL','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUWLK
      WRITE (LUWLK)
      WRITE (LUWLK)
      WRITE (LUWLK)
      WRITE (LUWLK)
      WRITE (LUWLK) NDCOOR, TDER
      CALL GPCLOSE(LUWLK,'KEEP')
C
      NMODIF = NRMCRD
C
      KFMATF = 1
      KFMATT = KFMATF + NCOOR*NCOOR*NCOOR
      KCORR  = KFMATT + NCOOR*NCOOR*NCOOR
      KPOS   = KCORR  + NCOOR
      KAMASS = KPOS   + NCOOR
      KLAST  = KAMASS + NCOOR
      LWRK  = LWORK   - KLAST
      CALL VIBV3(EVEC,FREQ,WORK(KFMATF),WORK(KFMATT),WORK(KCORR),
     &           WORK(KPOS),WORK(KAMASS),WORK(KLAST),LWRK,NCOOR,NDCOOR)
C
      RETURN
      END
C
C     /* Deck trdfsc */
      SUBROUTINE TDRFSC(TDER,CSTART,SYMCOR,TMPTDR,WORK,NCOOR,LWORK,
     &                  IPRINT)
C     ************************************************************
C     *** Subroutine that transforms the third derivative from ***
C     *** symmetry to cartesian coordinates (if not in normal  ***
C     *** coor).                                               ***
C     ************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "fcsym.h"
      CHARACTER*80 TEXT
      DIMENSION TDER  (NCOOR,NCOOR,NCOOR), TMPTDR(NCOOR,NCOOR,NCOOR),
     &          SYMCOR(NCOOR,NCOOR      ), CSTART(NCOOR            ),
     &          WORK  (LWORK            )
C
      KGRREP = 1
      KCHCTR = KGRREP +   NGORDR*NGVERT
      KICRIP = KCHCTR +   NGORDR*NCVERT
      KLAST  = KICRIP + 2*NCOOR
      LWRK   = LWORK  -   KLAST
      CALL GRPCHR(CSTART,SYMCOR,WORK(KGRREP),WORK(KCHCTR),WORK(KLAST),
     &            WORK(KICRIP),LWRK,IPRINT)
C
      LTXT = 9
      TEXT(1:9) = 'cartesian'
      CALL TRATDR(SYMCOR,TDER,TMPTDR,NCOOR,NCOOR,NCOOR,TEXT,LTXT,IPRINT)
      CALL DCOPY(NCOOR**3,TMPTDR,1,TDER,1)
C
      IF (IPRINT .GT.6) THEN
         CALL HEADER('Third derivative in cartesian coordinates',-1)
C
         IF (MOD(NCOOR,6).EQ.0) THEN
            NLCMAX = NCOOR/6
         ELSE
            NLCMAX = INT(NCOOR/6)+1
         END IF
C
         DO 400 ICOL2 = 1, NCOOR
            WRITE (LUPRI,'(A,I3)') '      Coloumn number', ICOL2
            WRITE (LUPRI,'(A)') '      -----------------'
            INLC = 0
            DO 500 INLCMX = 1, NLCMAX
               INLC2 = 6*(INLCMX-1) + 1
               INLC  = MIN(INLC+6,NCOOR)
               DO 600 ICOL1 = 1, NCOOR
                  WRITE (LUPRI,'(A,6F12.6)') '         ',
     &                             (TDER(I,ICOL1,ICOL2), I=INLC2, INLC)
 600           CONTINUE
               WRITE (LUPRI,'(A)') '                              '
 500        CONTINUE
 400     CONTINUE
      END IF
C
      RETURN
      END
C
C
C     /* Deck wravfl */
      SUBROUTINE WRAVFL(PRMTRX,NDIM1,NDIM2,NDIM3,PRPTXT,IPRINT)
#include "implicit.h"
#include "priunit.h"
C
#include "numder.h"
      CHARACTER*9 PRPTXT
      DIMENSION PRMTRX(NDIM1,NDIM2,NDIM3)
C
C     *********************************************
C     *** Keeping traack on how many properties ***
C     *** we are calculating.                   ***
C     *********************************************
C
      NMDPRP = NMDPRP + 1
C
C     ************************************
C     *** Writing the property to file ***
C     ************************************
C
      NTCOL = NDIM1/3 + 1
      IF (MOD(NDIM1,3).EQ.0) NTCOL = NDIM1/3
C
      WRITE (LUNDPR,'(A)') PRPTXT
      WRITE (LUNDPR,'(4I7)') NDIM1, NDIM2, NDIM3, NMPCAL
      DO IDIM3 = 1, NDIM3
         KDIM = 0
         WRITE (LUNDPR,'(A,I4)') 'The third dimension', IDIM3
         DO ITCOL = 1, NTCOL
            DO IDIM2 = 1, NDIM2
               WRITE (LUNDPR,'(3F36.16)') (PRMTRX(IDIM1,IDIM2,IDIM3),
     &                                 IDIM1 = KDIM+1,MIN(KDIM+3,NDIM1))

            END DO
            WRITE (LUNDPR,'(A)') '                       '
            KDIM = KDIM + 3
         END DO
      END DO
C
C     *******************
C     *** Test print. ***
C     *******************
C
      IF (IPRINT .GT. 20) THEN
         WRITE (LUPRI,'(A)') 'The property text :' // PRPTXT
         WRITE (LUPRI,'(A,I6)')
     &                   'This is property calculation number:',
     &                   NMDPRP
         DO IDIM3 = 1, NDIM3
            KDIM = 0
            WRITE (LUPRI,'(A,I4)') 'The third dimension', IDIM3
            DO ITCOL = 1, NTCOL
               DO IDIM2 = 1, NDIM2
                  WRITE (LUPRI,'(3F36.16)')
     &                  (PRMTRX(IDIM1,IDIM2,IDIM3),
     &                                 IDIM1 = KDIM+1,MIN(KDIM+3,NDIM1))
               END DO
               WRITE (LUPRI,'(A)') '                          '
               KDIM = KDIM + 3
            END DO
         END DO
      END IF
C
      RETURN
      END
C
C
C     /* Deck ndrdpp */
      SUBROUTINE NDRDPP(SPSPFV,TRLNFV,EXENFV,CCPRFV,TRAMAT,SYMCOR,
     &                  EXETMP,TRDTMP,EXCERF,WORK,NREDSM,LWORK,
     &                  IPRINT,CCPRP)
C     **************************************************
C     *** This routine keeps track of which property ***
C     *** that will be read from file next. The      ***
C     *** routine RDAVFL then reads it from file.    ***
C     **************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "cbiexc.h"
#include "inforb.h"
#include "abainf.h"
#include "trkoor.h"
#include "numder.h"
#include "gnrinf.h"
#include "pvibav.h"
#include "prpc.h"
      CHARACTER*9 PRPTXT
      LOGICAL CCPRP
      DIMENSION SPSPFV(NCOOR,NCOOR,6,NMPCAL),
     &          TRLNFV(3,NSYM,MXNEXI,NMPCAL),
     &          EXENFV(  NSYM,MXNEXI,NMPCAL), SYMCOR(NCOOR,NCOOR),
     &          TRAMAT(NCOOR,NCOOR,NMPCAL), CCPRFV(NPRPC,NMPCAL),
     &          EXCERF(NSYM,MXNEXI), EXETMP(NSYM,NTOTEX),
     &          TRDTMP(3,NSYM,NTOTEX),WORK(LWORK)
      DIMENSION NREDSM(NSYM,NSYM)
C
c      IF (DODIPS.OR.EXCIT) THEN
c         CALL IZERO(NREDSM,NSYM**2)
C
C     *** Read reference properties from file. ***
C
c         CALL RDREFG(EXCERF,NSYM,PRPTXT)
c      END IF
C
C     *** Rewind property file. ***
C
      REWIND(LUNDPR)
C
C     *** Number of properties each loop. ***
C
      NMEPRP = NMDPRP/NMPCAL
C
C     *************************************************
C     *** For all molecuar geometries (IDIME) and   ***
C     *** properties in each geometry (IMEPRP) read ***
C     *** from file.                                ***
C     *************************************************
C
      DO IMPCAL  = 2, NMPCAL
      DO IMEPRP = 1, NMEPRP
         READ (LUNDPR,'(A9)') PRPTXT
         READ (LUNDPR,'(4I7)') NDIM1, NDIM2, NDIM3, KGRB
C
C        *** Figures out which property read in, and reads it. ***
C

         CALL CHPRRD(SPSPFV(1,1,1,IMPCAL),TRAMAT(1,1,IMPCAL),
     &               TRLNFV(1,1,1,IMPCAL),EXENFV(1,1,IMPCAL),
     &               CCPRFV(1,IMPCAL),NDIM1,NDIM2,NDIM3,PRPTXT,CCPRP)
C
      END DO
      END DO
C
C     *** Original geometry at the end of the file. ***
C
      DO IMEPRP = 1, NMEPRP
         READ (LUNDPR,'(A9)') PRPTXT
         READ (LUNDPR,'(4I7)') NDIM1, NDIM2, NDIM3, KGRB
C
C        *** Figures out which property read in, and reads it. ***
C
         CALL CHPRRD(SPSPFV(1,1,1,1),TRAMAT(1,1,1),TRLNFV(1,1,1,1),
     &               EXENFV(1,1,1),CCPRFV(1,1),NDIM1,NDIM2,NDIM3,
     &               PRPTXT,CCPRP)
      END DO
C
C     *** Checkin sign on the transition dipole moment ***
C     *** derivatives.                                 ***
C
      IF (DOCCSD) THEN
         CALL CHKCCS(CCPRFV,IPRINT)
      ELSE
         CALL CHKSGN(TRLNFV,IPRINT)
      END IF
C
C     *** Transformation of all the matrices/vectors back ***
C     *** to original cartesian coordinates. Matrix in the***
C     *** original geometry (IMPCAL=1), need not to be    ***
C     *** calculated. Only when force constants are       ***
C     *** calculated from energies.                       ***
C
      IF (SPNSPN) THEN
         KSPTMP = 1
         KLAST  = KSPTMP + NCOOR**2
         IF (KLAST.GT.LWORK) CALL QUIT('Memory exceeded in OTRTEN')
C
         DO IMPCAL = 2, NMPCAL
            IF (IPRINT.GT.20) THEN
               WRITE (LUPRI,'(5X,A)') '                            '
               WRITE (LUPRI,'(5X,A)') 'Transformation of spin-spin.'
               WRITE (LUPRI,'(5X,A,I4)') 'Geometry number: ', IMPCAL
               WRITE (LUPRI,'(5X,A)') '                            '
            END IF
C
            DO ISPIN  = 1, 6
C
               CALL OTRTEN(SPSPFV(1,1,ISPIN,IMPCAL),TRAMAT(1,1,IMPCAL),
     &                     WORK(KSPTMP),NCOOR,NCOOR,NCOOR,IPRINT,'N',
     &                     'T')
            END DO
         END DO
      END IF
C
C     *** Test print. ***
C
      IF (IPRINT.GE.20) THEN
C
C        *** Workaround for common common block variables. ***
C
         CALL STPPVR
C
C        *** If cc-properties. ***
C
         IF (CCPRP) THEN
            CALL HEADER('CC-properties read from file.',0)
            DO IMPCAL = 1, NMPCAL
               CALL AROUND('New geometry:')
               WRITE (LUPRI,'(5X,A,I5)') 'Geometry number:', IMPCAL
               LUPRPCO = -1
               CALL GPOPEN(LUPRPCO,'CC_PRPC_O','UNKNOWN',
     &                     ' ','FORMATTED',IDUMMY,.FALSE.)
               CALL PRPRPC(LUPRPCO,2,CCPRFV(1,IMPCAL),NPRMI)
               CALL GPCLOSE(LUPRPCO,'KEEP')
            END DO
         END IF
C
C        *** If one property is spin-spin. ***
C
         IF (SPNSPN) THEN
            CALL HEADER('Spin-spin read from file.',0)
            DO IMPCAL = 1, NMPCAL
               WRITE (LUPRI,'(A,I6)')'Molecular geometry number:',IMPCAL
               CALL PRSPSP(SPSPFV(1,1,1,IMPCAL),NCOOR,NCOOR,LUPRI)
            END DO
         END IF
C
C        *** Transition dipole moments. ***
C
         IF (DODIPS) THEN
            CALL HEADER('Dipole transition moments read from file.',0)
            DO IMPCAL = 1, NMPCAL
               WRITE (LUPRI,'(1X,A,I6)')'Molecular geometry number:',
     &                IMPCAL
               CALL PRDPTR(TRLNFV(1,1,1,IMPCAL),EXENFV(1,1,IMPCAL),NSYM,
     &                     LUPRI)
            END DO
         END IF
C
C        *** Transformation matrices. ***
C
         CALL HEADER('Transformation matrices read from file.',0)
         DO IMPCAL = 1, NMPCAL
            WRITE (LUPRI,'(5X,A,I6)')'Molecular geometry number:',IMPCAL
            WRITE (LUPRI,'(A)') '                             '
            CALL PRTRMA(TRAMAT(1,1,IMPCAL),NCOOR,NCOOR,NCOOR,NCOOR,
     &                  LUPRI)
         END DO
      END IF
C
      RETURN
      END
C
C
C     /* Deck chprrd*/
      SUBROUTINE CHPRRD(SPSPFV,TRAMAT,TRLNFV,EXENFV,CCPRFV,NDIM1,NDIM2,
     &                  NDIM3,PRPTXT,CCPRP)
#include "implicit.h"
#include "mxcent.h"
C
#include "cbiexc.h"
#include "trkoor.h"
#include "inforb.h"
#include "prpc.h"
      LOGICAL CCPRP
      CHARACTER*9 PRPTXT
      DIMENSION SPSPFV(NCOOR,NCOOR,6), TRLNFV(3,NSYM,MXNEXI),
     &          TRAMAT(NCOOR,NCOOR  ), EXENFV(  NSYM,MXNEXI),
     &          CCPRFV(NPRPC)
C
      IF (PRPTXT(1:9).EQ.'CC PROPER') THEN
         CCPRP = .TRUE.
C
C        *** If excited state properties are calculated. ***
C        *** Find dalton symmetry of this derivative.    ***
C
         CALL RDAVFL(CCPRFV(1),NDIM1,NDIM2,NDIM3)
C
      ELSE IF (PRPTXT(1:9).EQ.'SPIN-SPIN') THEN
C
C        *** If property is spin-spin. ***
C
         CALL RDAVFL(SPSPFV(1,1,1),NDIM1,NDIM2,NDIM3)
      ELSE IF (PRPTXT(1:9).EQ.'CART-TRAN') THEN
C
C        *** If property is the transformation matrix ***
C        *** to old cartesian coordinates.            ***
C
         CALL RDAVFL(TRAMAT(1,1),NDIM1,NDIM2,NDIM3)
      ELSE IF ((PRPTXT(1:4).EQ.'TRDP').OR.(PRPTXT(1:4).EQ.'EXEN')) THEN
         IF (PRPTXT(7:7).EQ.' ') THEN
            READ (PRPTXT(5:6),'(2I1)') ISYM, IEXVAL
         ELSE
            READ (PRPTXT(5:7),'(I1,I2)') ISYM, IEXVAL
         END IF
         IF (PRPTXT(1:4).EQ.'TRDP') THEN
C
C           *** If property is HF transition moments. ***
C
            CALL RDAVFL(TRLNFV(1,ISYM,IEXVAL),3,1,1)
         ELSE
C
C           *** If property is HF excitation energies. ***
C
            CALL RDAVFL(EXENFV(ISYM,IEXVAL),1,1,1)
         END IF
      END IF
C
      RETURN
      END
C
C
C     /* Deck rdavfl*/
      SUBROUTINE RDAVFL(PRMTRX,NDIM1,NDIM2,NDIM3)
C     ****************************************************
C     *** This routine reads a property from vib. ave. ***
C     *** file, and stores it in PRMTRX.               ***
C     ****************************************************
#include "implicit.h"
C
#include "numder.h"
      CHARACTER*23 GRBG
      DIMENSION PRMTRX(NDIM1,NDIM2,NDIM3)
C
      NTCOL = NDIM1/3 + 1
      IF (MOD(NDIM1,3).EQ.0) NTCOL = NDIM1/3
C
      DO IDIM3 = 1, NDIM3
         KDIM = 0
         READ (LUNDPR,'(A)') GRBG
         DO ITCOL = 1, NTCOL
            DO IDIM2 = 1, NDIM2
               READ (LUNDPR,'(3F36.16)') (PRMTRX(IDIM1,IDIM2,IDIM3),
     &                                 IDIM1 = KDIM+1,MIN(KDIM+3,NDIM1))
            END DO
            READ (LUNDPR,'(A)') GRBG
            KDIM = KDIM + 3
         END DO
      END DO
C
      RETURN
      END
C
C
C     /* Deck otrten */
      SUBROUTINE OTRTEN(AMTRIX,TRAMAT,TMPMAT,NDIMA,NDIMT1,NDIMT2,
     &                  IPRINT,T1,T2)
#include "implicit.h"
      PARAMETER (D0 = 0.0D0, D1=1.0D0)
#include "maxorb.h"
C
#include "infpar.h"
#include "priunit.h"
      CHARACTER*1 T1, T2
      DIMENSION AMTRIX(NDIMA,NDIMA), TMPMAT(NDIMT1,NDIMT1),
     &          TRAMAT(NDIMT1,NDIMT2)
C
      CALL DZERO(TMPMAT,NDIMT1*NDIMT2)
      CALL DGEMM(T1,'N',NDIMT2,NDIMT1,NDIMT1,D1,TRAMAT,NDIMT1,AMTRIX,
     &           NDIMA,D0,TMPMAT,NDIMT1)
C
      CALL DZERO(AMTRIX,NDIMA**2)
      CALL DGEMM('N',T2,NDIMT2,NDIMT2,NDIMT1,D1,TMPMAT,NDIMT1,TRAMAT,
     &           NDIMT1,D0,AMTRIX,NDIMA)
C
      IF (IPRINT.GT.20) THEN
         WRITE (LUPRI,'(A)') 'Transformed tensor: '
         CALL PRTRMA(AMTRIX,NDIMA,NDIMA,NDIMT2,NDIMT2,LUPRI)
         WRITE (LUPRI,'(A)') 'Transformation tensor: '
         CALL PRTRMA(TRAMAT,NDIMT1,NDIMT2,NDIMT1,NDIMT2,LUPRI)
      END IF
C
      RETURN
      END
C
C
C     /* Deck ndwtpp */
      SUBROUTINE NDWTPP(SPSPFV,SPSPDR,NPPDER,IPRINT)
C     ***********************************************
C     *** This subroutine organize the writing of ***
C     *** properties and property derivatives for ***
C     *** vibrational averaging.                  ***
C     ***********************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "trkoor.h"
#include "numder.h"
#include "abainf.h"
      DIMENSION SPSPFV(NCOOR,NCOOR,6,NMPCAL),
     &          SPSPDR(NCOOR,NCOOR,6,NPPDER)
C
C     *** Opening file ***
C
      CALL GPOPEN(LUNDPR,'PROPERTY.NDER','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUNDPR)
C
C     *** Writing spin-spin values ***
C
      IF (SPNSPN) THEN
         CALL NVIBFL(SPSPFV(1,1,1,1),SPSPDR,NCOOR,NCOOR,6,NPPDER,IPRINT,
     &               'SPIN-SPIN')
      END IF
C
C     *** Closing file ***
C
      CALL GPCLOSE(LUNDPR,'KEEP')
C
      RETURN
      END
C
C
C     /* Deck nvibfl*/
      SUBROUTINE NVIBFL(FUNVAL,DERVAL,NDIM1,NDIM2,NDIM3,NDERVL,IPRINT,
     &                  PRPTXT)
C     **************************************************************
C     *** This writes the function value of the property in the  ***
C     *** starting geometry (FUNVAL) and the property derivative ***
C     *** (DERVAL) to file (i.e. the values needed for           ***
C     *** vibrational averaging.                                 ***
C     **************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL PRIVAL
      CHARACTER*9 PRPTXT
      DIMENSION FUNVAL(NDIM1,NDIM2,NDIM3       ),
     &          DERVAL(NDIM1,NDIM2,NDIM3,NDERVL)
#include "numder.h"
C
C     *** Function value ***
C
      CALL WRAVFL(FUNVAL,NDIM1,NDIM2,NDIM3,PRPTXT,IPRINT)
C
C     *** Derivatives ***
C
      IDERV = 0
      DO ICOOR = 1, NDCOOR
         IDERV = IDERV + 1
         WRITE (LUNDPR,'(A,I7)') 'Derivative:', ICOOR
         CALL WRAVFL(DERVAL(1,1,1,IDERV),NDIM1,NDIM2,NDIM3,'         ',
     &               IPRINT)
      END DO
      DO ICOOR2 = 1, NDCOOR
      DO ICOOR1 = 1, ICOOR2
         IDERV = IDERV + 1
C
C        *** If only the diagonal is needed. ***
C
         IF (PRPVIB.AND.((NARDRP+NMRDRP).EQ.2)) THEN
            PRIVAL = ICOOR1.EQ.ICOOR2
         ELSE
            PRIVAL = .TRUE.
         END IF
C
C        *** Writing the derivatives ***
C
         IF (PRIVAL) THEN
            WRITE (LUNDPR,'(A,2I6)') 'Derivative:', ICOOR1, ICOOR2
            CALL WRAVFL(DERVAL(1,1,1,IDERV),NDIM1,NDIM2,NDIM3,
     &                  '         ',IPRINT)
         END IF
      END DO
      END DO
C
      RETURN
      END
C
C
C     /* Deck ppvban */
      SUBROUTINE PPVBAN(DERIV,FREQ,SYMCOR,RNNORM,SPSPFV,SPSPDR,WORK,
     &                  NPPDER,NDERIV,LWORK,IPRINT)
C     *****************************************************
C     *** This routine reads nessecary information from ***
C     *** file, before doing the vibrational analysis.  ***
C     *****************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "trkoor.h"
#include "numder.h"
      CHARACTER*9 PRPTXT
      DIMENSION DERIV(NDERIV), SPSPDR(NCOOR,NCOOR,6,NPPDER),
     &          FREQ(NCOOR), SPSPFV(NCOOR,NCOOR,6), RNNORM(NCOOR),
     &          SYMCOR(NCOOR,NCOOR),WORK(LWORK)



C
C     *** Opening property file. ***
C
      LUNDPR = -1
      CALL GPOPEN(LUNDPR,'PROPERTY.NDER','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUNDPR)
C
      KNMPRP = NMEPRP-1
      DO IEPRP = 1, KNMPRP
         READ (LUNDPR,'(A)') PRPTXT
         READ (LUNDPR,'(4I7)') NDIM1, NDIM2, NDIM3, KGRB
C
C        *** If property is spin-spin ***
C
         IF (PRPTXT.EQ.'SPIN-SPIN') THEN
            KTNRDR = 2
            CALL PRPNVBA(DERIV,FREQ,SYMCOR,RNNORM,SPSPFV,SPSPDR,WORK,
     &                   NDIM1,NDIM2,NDIM3,NPPDER,NDERIV,LWORK,KTNRDR,
     &                   IPRINT,PRPTXT)
         END IF
      END DO
C
C     *** Closing property file. ***
C
      CALL GPCLOSE(LUNDPR,'KEEP')
C
      RETURN
      END
C
C
C     /* Deck prpnvba */
      SUBROUTINE PRPNVBA(FDERIV,FREQ,SYMCOR,RNNORM,FUNVAL,DERIVA,WORK,
     &                   NDIM1,NDIM2,NDIM3,NPPDER,NDERIV,LWORK,KTNRDR,
     &                   IPRINT,PRPTXT)
C     **************************************************************
C     *** This sorts out the nescecary information and calculate ***
C     *** vibrational average of property in funval.             ***
C     **************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "trkoor.h"
#include "numder.h"
      CHARACTER*9 PRPTXT
      DIMENSION FDERIV(NDERIV), DERIVA(NDIM1,NDIM2,NDIM3,NPPDER),
     &          FREQ(NCOOR), FUNVAL(NDIM1,NDIM2,NDIM3), RNNORM(NCOOR),
     &          SYMCOR(NCOOR,NCOOR), WORK(LWORK)
C
C     *** Reading nescecary information from file. ***
C
      CALL VBFIRD(FUNVAL,DERIVA,NDIM1,NDIM2,NDIM3,NPPDER,IPRINT)
C
C     *** Calculating the vibrational average. ***
C
CRF   Test print
      IF ( IPRINT .GE. 20 ) THEN
        WRITE(LUPRI,*) 'Used harmonical frequencies'
        DO I=1,NCOOR
          WRITE(LUPRI,'(i5,f18.10)') I, FREQ(I)
        END DO
      END IF

      KHARM  = 1
      KAHRM  = KHARM  + NDIM1*NDIM2*NDIM3
      KVIBAV = KAHRM  + NDIM1*NDIM2*NDIM3
      KTDER  = KVIBAV + NDIM1*NDIM2*NDIM3
      KLAST  = KTDER  + NDCOOR**3
      LWRK1 = LWORK - KLAST + 1
      IF (LWRK1.LT.1) CALL QUIT('Memory exceeded in CLCVBA')
      CALL CLCVBA(FDERIV,FREQ,SYMCOR,RNNORM,FUNVAL,DERIVA,WORK(KHARM),
     &            WORK(KAHRM),WORK(KVIBAV),WORK(KTDER),WORK(KLAST),
     &            NDERIV,NDIM1,NDIM2,NDIM3,NPPDER,KTNRDR,LWRK1,IPRINT,
     &            PRPTXT)
C
      RETURN
      END
C
C
C     /* Deck vbfird */
      SUBROUTINE VBFIRD(FUNVAL,DERIVA,NDIM1,NDIM2,NDIM3,NPPDER,IPRINT)
C     **************************************************************
C     *** This routine reads the nescecary information, from the ***
C     *** file PROPERTY.NDER, for the vibrational averaging.     ***
C     **************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "numder.h"
#include "trkoor.h"
      CHARACTER*9  PRPTXT
      CHARACTER*11 GRBG
      LOGICAL RDDVAL
      DIMENSION FUNVAL(NDIM1,NDIM2,NDIM3       ),
     &          DERIVA(NDIM1,NDIM2,NDIM3,NPPDER)
C
C     *** Function value ***
C
      CALL RDAVFL(FUNVAL,NDIM1,NDIM2,NDIM3)
C
C     *** 1. derivative ***
C
      IDERV = 0
      DO ICOOR1 = 1, NDCOOR
         READ (LUNDPR,'(A,I7)') GRBG, IC1
C
C        *** Test. ***
C
         IF (IC1.NE.ICOOR1) CALL QUIT('Something wrong with ' //
     &        'components in property derivative file.')
C
         READ (LUNDPR,'(A)') PRPTXT
         READ (LUNDPR,'(4I7)') ND1, ND2, ND3, KGRB
C
C        *** Test. ***
C
         IF ((ND1.NE.NDIM1).OR.(ND2.NE.NDIM2).OR.(ND3.NE.NDIM3))
     &        CALL QUIT('Something wrong with the dimensions in' //
     &                  'property derivative file.')
C
         IDERV = IDERV + 1
         CALL RDAVFL(DERIVA(1,1,1,IDERV),NDIM1,NDIM2,NDIM3)
      END DO
C
C     *** 2. derivative ***
C
      DO ICOOR2 = 1, NDCOOR
      DO ICOOR1 = 1, ICOOR2
C
C        *** If only the diagonal is needed. ***
C
         IF (PRPVIB.AND.((NARDRP+NMRDRP).EQ.2)) THEN
            RDDVAL = ICOOR1.EQ.ICOOR2
         ELSE
            RDDVAL = .TRUE.
         END IF
C
         IF (RDDVAL) THEN
            READ (LUNDPR,'(A,2I7)') GRBG, IC1, IC2
C
C           *** Test. ***
C
            IF ((IC1.NE.ICOOR1).OR.(IC2.NE.ICOOR2))
     &           CALL QUIT('Somethings wrong with components in' //
     &                     'property derivative file.')
C
            READ (LUNDPR,'(A)') PRPTXT
            READ (LUNDPR,'(3I7)') ND1, ND2, ND3
C
C           *** Test. ***
C
            IF ((ND1.NE.NDIM1).OR.(ND2.NE.NDIM2).OR.(ND3.NE.NDIM3))
     &           CALL QUIT('Something wrong with the dimensions in ' //
     &                     'property derivative file.')
            IDERV = IDERV + 1
            CALL RDAVFL(DERIVA(1,1,1,IDERV),NDIM1,NDIM2,NDIM3)
         END IF
      END DO
      END DO
C
      IF (IPRINT.GE.2) THEN
         CALL HEADER('Function value at reference geometry',0)
         CALL PRSPSP(FUNVAL,NDIM1,NDIM2,LUPRI)
         CALL HEADER('Numerical derivatives',0)
         CALL HEADER('1. Derivatives',0)
         IDERV = 0
         DO ICOOR1 = 1, NDCOOR
            IDERV = IDERV + 1
            WRITE (LUPRI,'(A,I7)') 'Derivative:', ICOOR1
            CALL PRSPSP(DERIVA(1,1,1,IDERV),NDIM1,NDIM2,LUPRI)
         END DO
         CALL HEADER('2. Derivatives',0)
         DO ICOOR2 = 1, NDCOOR
         DO ICOOR1 = 1, ICOOR2
C
C           *** If only the diagonal is used. ***
C
            IF (PRPVIB.AND.((NARDRP+NMRDRP).EQ.2)) THEN
               RDDVAL = ICOOR1.EQ.ICOOR2
            ELSE
               RDDVAL = .TRUE.
            END IF
C
            IF (RDDVAL) THEN
               IDERV = IDERV + 1
               WRITE (LUPRI,'(A,2I7)') 'Derivative:', ICOOR1, ICOOR2
               CALL PRSPSP(DERIVA(1,1,1,IDERV),NDIM1,NDIM2,LUPRI)
            END IF
         END DO
         END DO
      END IF
C
      RETURN
      END
C
C
C     /* Deck clcvba */
      SUBROUTINE CLCVBA(FDERIV,FREQ,SYMCOR,RNNORM,FUNVAL,DERIVA,HARMCO,
     &                  AHRMCO,VIBAVE,TDER,WORK,NDERIV,NDIM1,NDIM2,
     &                  NDIM3,NPPDER,KTNRDR,LWORK,IPRINT,PRPTXT)
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
#include "mxcent.h"
      PARAMETER (D0=0.0D0, D025 = 0.25D0, D15=1.5D0, D2 = 2.0D0)
#include "trkoor.h"
#include "numder.h"
#include "nuclei.h"
#include "spnout.h"
      CHARACTER*9 PRPTXT
      DIMENSION FDERIV(NDERIV), FREQ(NCOOR), TDER(NDCOOR,NDCOOR,NDCOOR),
     &          DERIVA(NDIM1,NDIM2,NDIM3,NPPDER), RNNORM(NCOOR),
     &          FUNVAL(NDIM1,NDIM2,NDIM3), HARMCO(NDIM1,NDIM2,NDIM3),
     &          AHRMCO(NDIM1,NDIM2,NDIM3), VIBAVE(NDIM1,NDIM2,NDIM3),
     &          SYMCOR(NCOOR,NCOOR), WORK(LWORK)
C
      CALL TITLER('Vibrational average of properties.','*',118)
C
      CALL DZERO(HARMCO,NDIM1*NDIM2*NDIM3)
      CALL DZERO(AHRMCO,NDIM1*NDIM2*NDIM3)
C
C     *** Read or arrange the force constants. ***
C    
      IF (C4FORC .AND. PRPONL) THEN
         KSORT  = 1
         KLAST  = KSORT + NDCOOR
         IF ( KLAST .GT. LWORK ) CALL QUIT ( 'Memory exceeded in '// 
     &                                'RDC4CU ')
         CALL RDC4CU(FREQ,TDER,WORK(KSORT), IPRINT)
      ELSE
         CALL STTDER(FDERIV,TDER,RNNORM,NDERIV,WORK,LWORK,IPRINT)
      END IF
C
C     *** Analysis of the contributing modes. ***
C
      IF (MDEANA) THEN
         KCMPVB = 1
         KLAST  = KCMPVB + NDIM1*NDIM2*NDIM3*NDCOOR
         LWRK1  = LWORK  - KLAST + 1
         CALL MODEAN(DERIVA,FREQ,TDER,VIBAVE,RNNORM,FUNVAL,WORK(KCMPVB),
     &               WORK(KLAST),NDIM1,NDIM2,NDIM3,NPPDER,LWRK1,IPRINT,
     &               PRPTXT)
      END IF
C
C     *** Contribution from the harmonic part ***
C     *** of the vibration.                   ***
C
      IDERV = NDCOOR
      IADD  = 0
      DO ICOOR1 = 1, NDCOOR
         IF (PRPVIB.AND.((NARDRP+NMRDRP).EQ.2)) THEN
            IDERV = IDERV + 1
         ELSE
            IADD  = IADD  + 1
            IDERV = IDERV + IADD
         END IF
C
         PFAC = D025*(RNNORM(ICOOR1)**2)/FREQ(ICOOR1)
C
         DO IDIM3 = 1, NDIM3
         DO IDIM2 = 1, NDIM2
         DO IDIM1 = 1, NDIM1
C
            HARMCO(IDIM1,IDIM2,IDIM3) = HARMCO(IDIM1,IDIM2,IDIM3      )
     &                           + PFAC*DERIVA(IDIM1,IDIM2,IDIM3,IDERV)
         END DO
         END DO
         END DO
      END DO
C
C     *** Contribution from the anharmonic part ***
C     *** of the vibration.                     ***
C
      DO ICOOR2 = 1, NDCOOR
         TMP = D0
         DO ICOOR1 = 1, NDCOOR
            TMP = TMP
     &          + TDER(ICOOR1,ICOOR1,ICOOR2)/FREQ(ICOOR1)
         END DO
         PFAC = -D025*TMP*RNNORM(ICOOR2)/(FREQ(ICOOR2)**D2)
C
         IF (IPRINT.GE.2) THEN
            WRITE (LUPRI,*) 'Prefactor for anharmonic contribution.'
            WRITE (LUPRI,'(A,I4,F8.4)') 'Coordinate:', ICOOR2, PFAC
         END IF
C
         DO IDIM3 = 1, NDIM3
         DO IDIM2 = 1, NDIM2
         DO IDIM1 = 1, NDIM1
            AHRMCO(IDIM1,IDIM2,IDIM3) = AHRMCO(IDIM1,IDIM2,IDIM3       )
     &                           + PFAC*DERIVA(IDIM1,IDIM2,IDIM3,ICOOR2)
         END DO
         END DO
         END DO
      END DO
C
C     *** Calculating the vibrational average of ***
C     *** the property                           ***
C
      DO IDIM3 = 1, NDIM3
      DO IDIM2 = 1, NDIM2
      DO IDIM1 = 1, NDIM1
         VIBAVE(IDIM1,IDIM2,IDIM3) = FUNVAL(IDIM1,IDIM2,IDIM3)
     &                             + HARMCO(IDIM1,IDIM2,IDIM3)
     &                             + AHRMCO(IDIM1,IDIM2,IDIM3)
      END DO
      END DO
      END DO
C
C     *** Results. ***
C
      IF (PRPTXT.EQ.'SPIN-SPIN') THEN
         CALL HEADER('Vibrationally averaged ' // PRPTXT //' couplings',
     &      0)
         KAVEIS = 1
         KANISO = KAVEIS + 5*NUCDEP*(NUCDEP+1)/2
         KASYM  = KANISO +   NUCDEP*(NUCDEP+1)/2
         KSPAR  = KASYM  +   NUCDEP*(NUCDEP+1)/2
         KAPAR =  KSPAR  +   NUCDEP*(NUCDEP+1)/2
         IF (SPNISO) THEN
            CALL ISOSPN(VIBAVE,WORK(KAVEIS),WORK(KANISO),WORK(KASYM),
     &                  WORK(KSPAR),WORK(KAPAR))
         ELSE
            CALL FNLSPN(VIBAVE,WORK(KAVEIS),WORK(KANISO),WORK(KASYM),
     &                  WORK(KSPAR),WORK(KAPAR))
         END IF
      END IF
C
      IF (IPRINT .GT.10) THEN
         DO IDIM3 = 1, NDIM3
            IF (NDIM3.GT.1) THEN
               WRITE (LUPRI,'(//A,I4)') 'Third dimension:   ', IDIM3
            END IF
C
            CALL OUTPUT(VIBAVE(1,1,IDIM3),1,NDIM1,1,NDIM2,NDIM1,NDIM2,
     &                  -1,LUPRI)
C
         END DO
         WRITE (LUPRI,'(/)')
      END IF
C
C     *** Print. ***
C
      IF (PRPTXT.EQ.'SPIN-SPIN') THEN
         CALL HEADER('Harmonic contribution to vibrational average of'//
     &               ' property ' // PRPTXT,0)
         KAVEIS = 1
         KANISO = KAVEIS + 5*NUCDEP*(NUCDEP+1)/2
         KASYM  = KANISO +   NUCDEP*(NUCDEP+1)/2
         KSPAR  = KASYM  +   NUCDEP*(NUCDEP+1)/2
         KAPAR =  KSPAR  +   NUCDEP*(NUCDEP+1)/2
         IF (SPNISO) THEN
            CALL ISOSPN(HARMCO,WORK(KAVEIS),WORK(KANISO),WORK(KASYM),
     &                  WORK(KSPAR),WORK(KAPAR))
         ELSE
            CALL FNLSPN(HARMCO,WORK(KAVEIS),WORK(KANISO),WORK(KASYM),
     &                  WORK(KSPAR),WORK(KAPAR))
         END IF
      END IF
      IF (IPRINT.GT.20) THEN
         DO IDIM3 = 1, NDIM3
            IF (NDIM3.GT.1) THEN
               WRITE (LUPRI,'(//A,I4)') 'Third dimension:   ', IDIM3
            END IF
C
            CALL OUTPUT(HARMCO(1,1,IDIM3),1,NDIM1,1,NDIM2,NDIM1,NDIM2,
     &                  -1,LUPRI)
         END DO
      END IF
C
      IF (PRPTXT.EQ.'SPIN-SPIN') THEN
         CALL HEADER('Anharmonic contribution to vibrational average' //
     &               ' of property ' // PRPTXT,0)
         KAVEIS = 1
         KANISO = KAVEIS + 5*NUCDEP*(NUCDEP+1)/2
         KASYM  = KANISO +   NUCDEP*(NUCDEP+1)/2
         KSPAR  = KASYM  +   NUCDEP*(NUCDEP+1)/2
         KAPAR =  KSPAR  +   NUCDEP*(NUCDEP+1)/2
         IF (SPNISO) THEN
            CALL ISOSPN(AHRMCO,WORK(KAVEIS),WORK(KANISO),WORK(KASYM),
     &                  WORK(KSPAR),WORK(KAPAR))
         ELSE
            CALL FNLSPN(AHRMCO,WORK(KAVEIS),WORK(KANISO),WORK(KASYM),
     &                  WORK(KSPAR),WORK(KAPAR))
         END IF
      END IF
      IF (IPRINT.GT.20) THEN
         DO IDIM3 = 1, NDIM3
            IF (NDIM3.GT.1) THEN
               WRITE (LUPRI,'(//A,I4)') 'Third dimension:   ', IDIM3
            END IF
C
            CALL OUTPUT(AHRMCO(1,1,IDIM3),1,NDIM1,1,NDIM2,NDIM1,NDIM2,1,
     &                  LUPRI)
         END DO
C
         CALL HEADER('Function value of property ' // PRPTXT //
     &               ' at reference geometry',0)
         DO IDIM3 = 1, NDIM3
            IF (NDIM3.GT.1) THEN
               WRITE (LUPRI,'(//A,I4)') 'Third dimension:   ', IDIM3
            END IF
C
            CALL OUTPUT(FUNVAL(1,1,IDIM3),1,NDIM1,1,NDIM2,NDIM1,NDIM2,1,
     &                  LUPRI)
         END DO
C
      END IF
C
      RETURN
      END
C
C
C  /* Deck sttder */
      SUBROUTINE STTDER(DERIV,TDER,RNNORM,NDERIV,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
CRF   added file routine -- appears that we must include inftap.h
#include "inftap.h"
      PARAMETER (D1=1.0D0)
#include "trkoor.h"
#include "numder.h"
#include "cbinum.h"
      INTEGER   NDCOOT
      DIMENSION DERIV(NDERIV), TDER(NDCOOR,NDCOOR,NDCOOR),
     &          RNNORM(NCOOR), WORK(LWORK) 
C     
CRF   The work done here has already been done in ahrmvb 
CRF   If this was a calculation of properties seperate from the calculation
CRF   of the forcefield, we need to read it from there (and perhaps anyway)?

      IF (PRPONL) THEN
C     *** Memory for temporary array ***
        KTEMP  = 1
        KLAST  = KTEMP + NCOOR**3
        IF (KLAST .GT. LWORK) 
     &           CALL QUIT('Memory exceeded in STTDER')
        CALL GPOPEN(LUWLK,ABAWLK,'OLD','SEQUENTIAL','UNFORMATTED',
     &              IDUMMY, .FALSE.)
        REWIND (LUWLK)
        READ (LUWLK)
        READ (LUWLK)
        READ (LUWLK)
        READ (LUWLK)
        READ (LUWLK) NDCOOT, WORK(KTEMP:KLAST-1)
        CALL GPCLOSE(LUWLK,'KEEP')
        IF ( NDCOOT .NE. NDCOOR) WRITE (LUPRI, '(3A)')
     &     'Check that force field in file',
     &     ABAWLK, 'matches current geometry'
C     *** Store the force constants in TDER ***
        DO 100 ICOOR1 = 1, NDCOOR
        DO 100 ICOOR2 = 1, NDCOOR
        DO 100 ICOOR3 = 1, NDCOOR
           ILOC = (ICOOR1-1)*NCOOR**2 + (ICOOR2-1)*NCOOR + ICOOR3
           TDER (ICOOR1,ICOOR2,ICOOR3) = WORK( ILOC )
100     CONTINUE
      ELSE
        DFAC = D1
        IDERIV = 0
        DO 200 ICOOR3 = 1, NDCOOR
        DO 200 ICOOR2 = 1, ICOOR3
        DO 200 ICOOR1 = 1, ICOOR2
           IF (NRMCRD) DFAC =
     &                  RNNORM(ICOOR1)*RNNORM(ICOOR2)*RNNORM(ICOOR3)
           IDERIV = IDERIV + 1
           TDER(ICOOR1,ICOOR2,ICOOR3) = DERIV(IDERIV)*DFAC
           TDER(ICOOR1,ICOOR3,ICOOR2) = DERIV(IDERIV)*DFAC
           TDER(ICOOR2,ICOOR1,ICOOR3) = DERIV(IDERIV)*DFAC
           TDER(ICOOR2,ICOOR3,ICOOR1) = DERIV(IDERIV)*DFAC
           TDER(ICOOR3,ICOOR1,ICOOR2) = DERIV(IDERIV)*DFAC
           TDER(ICOOR3,ICOOR2,ICOOR1) = DERIV(IDERIV)*DFAC
 200    CONTINUE
      END IF
C
      IF (IPRINT.GE.20) THEN
         CALL HEADER('Diagonal of cubic force field, F(I,J,J)',-1)
C
         ISTRT = 1
         KCOL  = 9
         LAST  = MIN(NDCOOR,KCOL)
         KCOOR = NDCOOR
         NCOL  = NDCOOR/KCOL
         IF (MOD(NDCOOR,KCOL).NE.0) NCOL = NCOL + 1
C
         DO ICOL = 1, NCOL
            DO ICOOR = 1, NDCOOR
               WRITE (LUPRI,'(5X,9F12.6)')
     &                        (TDER(ICOOR,I,I),I=ISTRT,LAST)
            END DO
            WRITE (LUPRI,'()')
            ISTRT = ISTRT + KCOL
            LAST  = MIN(NDCOOR,KCOL+LAST)
         END DO
      END IF
      RETURN
      END
C
C
C     /* Deck modean */
      SUBROUTINE MODEAN(DERIVA,FREQ,TDER,VIBAVE,RNNORM,FUNVAL,CMPVIB,
     &                  WORK,NDIM1,NDIM2,NDIM3,NPPDER,LWORK,IPRINT,
     &                  PRPTXT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0=0.0D0, D025 = 0.25D0, D2 = 2.0D0)
#include "trkoor.h"
#include "numder.h"
#include "nuclei.h"
#include "spnout.h"
      CHARACTER*9 PRPTXT
      DIMENSION CMPVIB(NDIM1,NDIM2,NDIM3,NDCOOR),
     &          DERIVA(NDIM1,NDIM2,NDIM3,NPPDER), FREQ(NCOOR),
     &          TDER(NDCOOR,NDCOOR,NDCOOR), VIBAVE(NDIM1,NDIM2,NDIM3),
     &          RNNORM(NCOOR), FUNVAL(NDIM1,NDIM2,NDIM3), WORK(LWORK)
C
      CALL DZERO(CMPVIB,NDIM1*NDIM2*NDIM3*NDCOOR)
C
      IDERV = NDCOOR
      IADD  = 0
      DO ICOOR1 = 1, NDCOOR
         IF (PRPVIB.AND.((NARDRP+NMRDRP).EQ.2)) THEN
            IDERV = IDERV + 1
         ELSE
            IADD  = IADD  + 1
            IDERV = IDERV + IADD
         END IF
C
         PFAC = D025*(RNNORM(ICOOR1)**2)/FREQ(ICOOR1)
C
         DO IDIM3 = 1, NDIM3
         DO IDIM2 = 1, NDIM2
         DO IDIM1 = 1, NDIM1
C
            CMPVIB(IDIM1,IDIM2,IDIM3,ICOOR1) =
     &                                 CMPVIB(IDIM1,IDIM2,IDIM3,ICOOR1)
     &                          + PFAC*DERIVA(IDIM1,IDIM2,IDIM3,IDERV )
         END DO
         END DO
         END DO
      END DO
C
C     *** Contribution from the anharmonic part ***
C     *** of the vibration.                     ***
C
      DO ICOOR2 = 1, NDCOOR
         TMP = D0
         DO ICOOR1 = 1, NDCOOR
            TMP = TMP
     &          + TDER(ICOOR1,ICOOR1,ICOOR2)/FREQ(ICOOR1)
         END DO
         PFAC = -D025*TMP*RNNORM(ICOOR2)/(FREQ(ICOOR2)**D2)
C
         DO IDIM3 = 1, NDIM3
         DO IDIM2 = 1, NDIM2
         DO IDIM1 = 1, NDIM1
            CMPVIB(IDIM1,IDIM2,IDIM3,ICOOR2) =
     &                                  CMPVIB(IDIM1,IDIM2,IDIM3,ICOOR2)
     &                           + PFAC*DERIVA(IDIM1,IDIM2,IDIM3,ICOOR2)
         END DO
         END DO
         END DO
      END DO
cC
cC     *** Calculating the vibrational average of ***
cC     *** the property                           ***
cC
c      DO IDIM3 = 1, NDIM3
c      DO IDIM2 = 1, NDIM2
c      DO IDIM1 = 1, NDIM1
c         VIBAVE(IDIM1,IDIM2,IDIM3) = FUNVAL(IDIM1,IDIM2,IDIM3)
c     &                             + HARMCO(IDIM1,IDIM2,IDIM3)
c     &                             + AHRMCO(IDIM1,IDIM2,IDIM3)
c      END DO
c      END DO
c      END DO
C
C     *** Results. ***
C
      IF (PRPTXT.EQ.'SPIN-SPIN') THEN
         CALL HEADER('Mode analysis for spin-spin coupling const.',0)
         KAVEIS = 1
         KANISO = KAVEIS + 5*NUCDEP*(NUCDEP+1)/2
         KASYM  = KANISO +   NUCDEP*(NUCDEP+1)/2
         KSPAR  = KASYM  +   NUCDEP*(NUCDEP+1)/2
         KAPAR =  KSPAR  +   NUCDEP*(NUCDEP+1)/2
         DO ICOOR = 1, NDCOOR
            WRITE (LUPRI,'(10X,A)') '                   '
            WRITE (LUPRI,'(10X,A)') '                   '
            WRITE (LUPRI,'(10X,A,I4)') 'Mode number:  ', ICOOR
            WRITE (LUPRI,'(10X,A)') '------------------'
            WRITE (LUPRI,'(10X,A)') '                   '
            WRITE (LUPRI,'(10X,A)') '                   '
            IF (SPNISO) THEN
               CALL ISOSPN(CMPVIB(1,1,1,ICOOR),WORK(KAVEIS),
     &              WORK(KANISO),WORK(KASYM),WORK(KSPAR),WORK(KAPAR))
            ELSE
               CALL FNLSPN(CMPVIB(1,1,1,ICOOR),WORK(KAVEIS),
     &              WORK(KANISO),WORK(KASYM),WORK(KSPAR),WORK(KAPAR))
            END IF
         END DO
      END IF
C
C     *** Debug print ***
C
      IF (IPRINT .GT.10) THEN
         DO IDIM4 = 1, NDCOOR
            WRITE (LUPRI,'(A)') '                               '
            WRITE (LUPRI,'(A)') '                               '
            WRITE (LUPRI,'(A,I4)') 'Coordinate:   ', IDIM4
            DO IDIM3 = 1, NDIM3
               IF (NDIM3.GT.1) THEN
                  WRITE (LUPRI,'(A)') '                               '
                  WRITE (LUPRI,'(A)') '                               '
                  WRITE (LUPRI,'(A,I4)') 'Third dimension:   ', IDIM3
               END IF
C
               CALL OUTPUT(CMPVIB(1,1,IDIM3,IDIM4),1,NDIM1,1,NDIM2,
     &                     NDIM1,NDIM2,1,LUPRI)
C
            END DO
            WRITE (LUPRI,'(A)') '                               '
            WRITE (LUPRI,'(A)') '                               '
         END DO
      END IF
C
      RETURN
      END
C
C     /* Deck rdc4cu */
      SUBROUTINE RDC4CU(FREQ,TDER,SORT,IPRINT)
C     **********************************************
C     *** Routine that reads a cubic force field ***
C     *** from the CFOUR style "cubic" file      ***
C     **********************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
#include "mxcent.h"
#include "inftap.h"

#include "trkoor.h"
#include "numder.h"
      PARAMETER (IPRTRH = 20 )
      LOGICAL USED, CUBEXS
      INTEGER SORTAR
      DIMENSION TDER(NDCOOR,NDCOOR,NDCOOR),
     &          SORTAR(NCOOR-NDCOOR+1:NCOOR),
     &          USED(NCOOR), FREQ(NCOOR)
   
      LUC4IF = -399 

      INQUIRE (FILE='cubic',EXIST=CUBEXS)
      IF (.NOT. CUBEXS ) CALL QUIT ('File "cubic" does not exist')      

      CALL GPOPEN(LUC4IF,'cubic','OLD','SEQUENTIAL','FORMATTED',
     &            IDUMMY, .FALSE.)
      USED = .FALSE.
      SORT = 0
      CALL DZERO(TDER, NDCOOR**3) 

C     *** CFOUR orders normal coordinates by increasing harmonic frequency ***
C     *** We need to match the CFOUR indecies to the dalton ones.          ***

      DO IDX = NCOOR,NCOOR-NDCOOR+1,-1
         JTEMP  = 0
         TMPFRQ = 0
c     *** Find the index corresponding to the highest frequency not yet used ***
         DO JDX = 1, NDCOOR
            IF ( USED(JDX) .OR. (TMPFRQ .GE. FREQ(JDX) ) ) CYCLE
            JTEMP  = JDX
            TMPFRQ = FREQ (JDX) 
         END DO
         IF ( JTEMP .LE. 0 ) EXIT ! None found
         SORTAR ( IDX ) = JTEMP 
         USED (JTEMP)   = .TRUE.
      END DO

C     *** Read the cubic forcefield and add it to TDER ***       
      DO 
         READ (LUC4IF,'(3I5,F25.10)',END=100) IDX , JDX , KDX , TMPFF
         IF ( (IDX .GT. NCOOR) .OR. (JDX .GT. NCOOR) .OR. 
     &        (KDX .GT. NCOOR) ) CALL QUIT ('Index out of range in ' //
     &                                'file "cubic"' )
         PRINT *, FREQ ( SORTAR( IDX) )
         TMPFF = - TMPFF *  SQRT( FREQ(SORTAR(IDX)) 
     &            * FREQ(SORTAR(JDX)) * FREQ(SORTAR(KDX)) ) / XTKAYS 
         TDER ( SORTAR(IDX), SORTAR(JDX), SORTAR(KDX ) ) = TMPFF
         TDER ( SORTAR(IDX), SORTAR(KDX), SORTAR(JDX ) ) = TMPFF
         TDER ( SORTAR(JDX), SORTAR(IDX), SORTAR(KDX ) ) = TMPFF
         TDER ( SORTAR(JDX), SORTAR(KDX), SORTAR(IDX ) ) = TMPFF
         TDER ( SORTAR(KDX), SORTAR(IDX), SORTAR(JDX ) ) = TMPFF
         TDER ( SORTAR(KDX), SORTAR(JDX), SORTAR(IDX ) ) = TMPFF
      END DO

 100  CONTINUE
      CALL GPCLOSE (LUC4IF,'KEEP')
 
C     *** Print the cubic forcefield ***
      IF ( IPRINT .GE. IPRTRH ) THEN
         CALL HEADER('Diagonal of cubic force field, F(I,J,J)',-1)
         DO IDX = 1, NDCOOR
            DO JDX = 1, NDCOOR, 6          
               WRITE (LUPRI, '(5X,6F12.6)')  ( TDER(IDX,KDX,KDX),
     &                KDX = JDX , MIN (NDCOOR , JDX + 5) )
            END DO
            WRITE (LUPRI, '(A)') '  ' 
         END DO
      END IF


      RETURN
      END !SUBROUTINE RDC4CU


C
C
#ifdef NOT_USED_IN_THIS_DALTON_VERSION
C     /* Deck rdrefg */
      SUBROUTINE RDREFG(EXCERF,NSYM,PRPTXT)
#include "implicit.h"
#include "priunit.h"
C
#include "cbiexc.h"
      CHARACTER*9 PRPTXT
      DIMENSION EXCERF(NSYM,MXNEXI)
C
C     *** Open ref file. ***
C
      LUREF = -1
      CALL GPOPEN(LUREF,'PROPERTY.REF','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUREF)
C
C     *** Reading properties at reference geometry from file. ***
C
      DO IS = 1, NSYM
      DO IE = 1, NEXCIT(IS)
         READ (LUREF,'(A9)') PRPTXT
         READ (LUREF,'(3I7)') NDIM1, NDIM2, NDIM3
         IF (PRPTXT(7:7).EQ.' ') THEN
            READ (PRPTXT(5:6),'(2I1)') ISYM, IEXVAL
         ELSE
            READ (PRPTXT(5:7),'(I1,I2)') ISYM, IEXVAL
         END IF
         CALL RDAVFL(EXCERF(ISYM,IEXVAL),1,1,1)
      END DO
      END DO
C
C     *** Closing file. ***
C
      CALL GPCLOSE(LUREF,'DELETE')
C
C     *** Debug print. ***
C
      IF (IPRINT.GT.20) THEN
         CALL HEADER('Excitation energies at ref. geom.',0)
         CALL PREXCE(EXCERF,NSYM,LUPRI)
      END IF
C
      RETURN
      END
#endif  /* NOT_USED_IN_THIS_DALTON_VERSION */
