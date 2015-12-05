!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C
c*DECK CC_GETGD
      SUBROUTINE CC_GETGD(GD,ILSTNR,ISIDE,LIST,WORK,LWORK)
C
C----------------------------------------------------------------------
C
C   Purpose: Calculate right hand side vectors for response equations.
C
C     Written by Ove Christiansen 31-10-1996
C     second- and third-order response, Christof Haettig, spring 1997 
C     excited state response, Christof Haettig, summer 1997
C     left Cauchy vectors, Christof Haettig, Oct. 1997
C     second-order Cauchy vectors, Christof Haettig, Feb. 1998
C     CCSLV98,OC mai 98.
C     projected 1st-order lagrang. multipl. (PL1), Sonia Coriani, 2000  
C     changes for CC3 response, Christof Haettig, Spring 2002
C     changes for CCr12 response, Christof Haettig, Summer 2003
C
C----------------------------------------------------------------------
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER( ONE = 1.0D00 ,XMONE = -1.0D0 )
      DIMENSION GD(*),WORK(LWORK)
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
#include "dummy.h"
#include "iratdef.h"
#include "cclr.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccsections.h"
#include "ccexgr.h"
#include "cclrmrsp.h"
#include "ccr1rsp.h"
#include "ccr3rsp.h"
#include "ccr2rsp.h"
#include "ccl1rsp.h"
#include "ccl2rsp.h"
#include "ccl3rsp.h"
#include "ccn2rsp.h"
#include "ccrc1rsp.h"
#include "cclc1rsp.h"
#include "cccr2rsp.h"
#include "cccl2rsp.h"
#include "leinf.h"
#include "ccpl1rsp.h"
Cholesky
#include "ccdeco.h"
Cholesky
C
      CHARACTER MODEL*(10), LIST*(*)
C
      CALL QENTER('CC_GETGD')
C
      IF ( IPRINT .GT. 10 .OR. LOCDBG) THEN
         CALL AROUND( 'CC_GD: Constructing GD vector ')
         WRITE(LUPRI,*) 'LIST  :',LIST
         WRITE(LUPRI,*) 'ILSTNR:',ILSTNR
         CALL FLSHFO(LUPRI)
      ENDIF

C     ------------------------------------
C     set symmetry of the gradient vector:
C     ------------------------------------
      ISYMGD = ILSTSYM(LIST,ILSTNR)

C     -------------------------------------------------
C     set option for readin gradient vector:
C          1 - read singles only (CCS,CIS)
C          3 - read singles and doubles (CC2,CCSD)
C         24 - read effective singles and doubles (CC3)
C     -------------------------------------------------
      IOPTRDG = 3
      IF (CC3) IOPTRDG = 24
      IF (CCS .OR. CIS) IOPTRDG = 1
      IOPTR12 = 32

C     ---------------------------------------------------------------
C     set the start addresses of the single and the double excitation 
C     parts of the gradient vector and the total length of the vector
C     ---------------------------------------------------------------
      IADR1  = 1
      IADR2  = IADR1 + NT1AM(ISYMGD)
      IADR12 = IADR2 + NT2AM(ISYMGD)
      NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
      IF ( CCS ) THEN
        NTAMP = NT1AM(ISYMGD)
      ELSE IF ( CCR12 ) THEN
        NTAMP = NTAMP + NTR12AM(ISYMGD) 
      END IF
 
Cholesky
      IF (CHOINT .AND. CC2) THEN
         NTAMP = NT1AM(ISYMGD)
      END IF
Cholesky

C----------------------------------------------------------------------
C 'R3 ' third-order response t amplitudes:  --> 'O3 '
C----------------------------------------------------------------------
      IF (LIST(1:2) .EQ. 'R3') THEN
         IF (ISIDE .NE. 1) 
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for R3 response:'
            WRITE(LUPRI,'(1x,A,3F12.6,A,3A8,A,4I3)')
     *         'Frequencies',(FRQR3T(ILSTNR,I),I=1,3),
     *         ' and labels: ',(LBLR3T(ILSTNR,I),I=1,3),
     *         ' and symmetries: ',(ISYR3T(ILSTNR,I),I=1,3), ISYMGD
         END IF
         IDLSTO3  = IRHSR3(LBLR3T(ILSTNR,1),FRQR3T(ILSTNR,1),ISYM1,
     &                     LBLR3T(ILSTNR,2),FRQR3T(ILSTNR,2),ISYM2,
     &                     LBLR3T(ILSTNR,3),FRQR3T(ILSTNR,3),ISYM3)
         CALL CC_RDRSP('O3',IDLSTO3,ISYMGD,IOPTRDG,MODEL,
     &                 GD(IADR1),GD(IADR2))
         CALL CCLR_DIASCL(GD(IADR2),0.5d0,ISYMGD)
         IF (CCR12) CALL CC_RDRSP('O3',IDLSTO3,ISYMGD,IOPTR12,MODEL,
     &                            DUMMY,GD(IADR12))
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)
         IF (CCR12) NTAMP = NTAMP + NTR12AM(ISYMGD)

C----------------------------------------------------------------------
C 'R2 ' second-order response t amplitudes:  --> 'O2 '
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'R2') THEN
         IF (ISIDE .NE. 1) 
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = MULD2H(ISYAR2T(ILSTNR),ISYBR2T(ILSTNR))
         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for '
            WRITE(LUPRI,'(1x,2(A,F10.6),2(A,A8),A,3I3)')
     *         'Frequencies',FRQAR2T(ILSTNR),',',FRQBR2T(ILSTNR),
     *         ' and labels: ',LBLAR2T(ILSTNR),',',LBLBR2T(ILSTNR),
     *         ' and symmetries: ',ISYAR2T(ILSTNR),ISYBR2T(ILSTNR),
     *         ISYMGD
         END IF
         CALL CC_R2GD(GD,ILSTNR,WORK,LWORK)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)
         IF (CCR12) NTAMP = NTAMP + NTR12AM(ISYMGD)
 
C----------------------------------------------------------------------
C 'R1 ' first-order response t amplitudes:  --> 'O1 '
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'R1') THEN

         IF (ISIDE.NE.1)
     &      CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         IF (IPRINT.GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for '
            WRITE(LUPRI,'(1x,A,F10.6,A,A8,A,I3)')
     *           'Frequency',FRQLRT(ILSTNR),
     *           ' and label: ',LRTLBL(ILSTNR),
     *           ' and symmetry: ', ISYMGD
         ENDIF

         ILSTGD = IRHSR1(LRTLBL(ILSTNR),LORXLRT(ILSTNR),
     &                   FRQLRT(ILSTNR),ISYM1)

C       For Cholesky CC2 we need to add a frequency-dependent term to
C       get the effective right-hand side. [Which is also why the
C       solver only solves equations for 1 frequency at a time.]
C       Try to read the eff. rhs from disk; if it fails, calculate
C       the eff. rhs and write it to disk for later use.
C       Conventional: simply read from disk.
C       -------------------------------------------------------------

        IF (CHOINT .AND. CC2) THEN
           MODEL = 'CC2       '
           IOPT = 33
           CALL CC_RDRSP('eO1',ILSTNR,ISYMGD,IOPT,MODEL,
     &                   GD,DUMMY)
           IF (IOPT .EQ. 33) THEN
              IF (IPRINT .GT. 2) THEN
                 WRITE(LUPRI,'(1X,A,A)')
     &           'Effective Cholesky gradient not found on file. ',
     &           'Calculating gradient...'
              ENDIF
              IALG = 0
              CALL CC_CHOXI1(GD,WORK,LWORK,ILSTGD,LRTLBL(ILSTNR),
     &                       FRQLRT(ILSTNR),ISYMGD,IALG)
              IOPT = 1
              CALL CC_WRRSP('eO1',ILSTNR,ISYMGD,IOPT,MODEL,
     &                      DUMMY,GD,DUMMY,WORK,LWORK)
           ENDIF

         ELSE         !Conventional

           CALL CC_RDRSP('O1 ',ILSTGD,ISYMGD,IOPTRDG,MODEL,
     &                    GD(IADR1),GD(IADR2))
           IF (CCR12) THEN
             CALL CC_RDRSP('O1 ',ILSTGD,ISYMGD,IOPTR12,MODEL,
     &                     DUMMY,GD(IADR12))
           END IF
         END IF

         IF (CCSLV) THEN
            IF (CCSDT) CALL QUIT('Solvent not implemented for triples')
            IF (CCR12) CALL QUIT('Solvent not implemented for CCR12')
            CALL CC_PTB(GD,LRTLBL(ILSTNR),ISYLRT(ILSTNR),
     *                  FRQLRT(ILSTNR),LORXLRT(ILSTNR),
     *                  WORK,LWORK)
         ENDIF
         IF (LOCDBG) THEN
           WRITE(LUPRI,*)'CC_GETGD> RHS vector "O1":',ISYMGD,NTAMP
           CALL OUTPUT(GD(IADR1),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
         END IF

C----------------------------------------------------------------------
C 'L3 ' third-order response Lagrangian multipliers: --> 'X3 ' + 'F3 '
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'L3') THEN
         IF (ISIDE .NE. -1) 
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for L3 response:'
            WRITE(LUPRI,'(1x,A,3F12.6,A,3A8,A,4I3)')
     *         'Frequencies',(FRQL3(ILSTNR,I),I=1,3),
     *         ' and labels: ',(LBLL3(ILSTNR,I),I=1,3),
     *         ' and symmetries: ',(ISYL3(ILSTNR,I),I=1,3), ISYMGD
         END IF
         CALL CC_L3GD(GD,ILSTNR,WORK,LWORK)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)
         IF ( CCR12 ) CALL QUIT('No R12 for X3+F3')

C----------------------------------------------------------------------
C 'L2 ' second-order response Lagrangian multipliers: --> 'X2 ' + 'F2 '
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'L2') THEN
         IF (ISIDE .NE. -1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for L2 response:'
            WRITE(LUPRI,'(1x,2(A,F10.6),2(A,A8),A,3I3)')
     *         'Frequencies',FRQAL2(ILSTNR),',',FRQBL2(ILSTNR),
     *         ' and labels: ',LBLAL2(ILSTNR),',',LBLBL2(ILSTNR),
     *         ' and symmetries: ',ISYAL2(ILSTNR),ISYBL2(ILSTNR),
     *         ISYMGD
         END IF
         CALL CC_L2GD(GD,ILSTNR,WORK,LWORK)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)
         IF ( CCR12) NTAMP = NTAMP + NTR12AM(ISYMGD)

C----------------------------------------------------------------------
C 'L1 ' first-order response Lagrangian multipliers: --> 'X1 ' + 'F1 '
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'L1') THEN
         IF (ISIDE .NE.-1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(/,1x,A)')
     *           'Finding gradient for '
                 WRITE(LUPRI,'(1x,A,F10.6,A,A8,A,I3)')
     *             'Frequency',FRQLRZ(ILSTNR),
     *           ' and label: ',LRZLBL(ILSTNR),' and symmetry: ',
     *            ISYLRZ(ILSTNR)
         ENDIF
         CALL CC_LGD(GD,LRZLBL(ILSTNR),ISYLRZ(ILSTNR),
     *               FRQLRZ(ILSTNR),LORXLRZ(ILSTNR),IOPTRDG,WORK,LWORK)
         NTAMP = NT1AM(ISYLRZ(ILSTNR)) + NT2AM(ISYLRZ(ILSTNR))
         IF ( CCS ) NTAMP = NT1AM(ISYLRZ(ILSTNR))
         IF ( CCR12 ) NTAMP = NTAMP + NTR12AM(ISYLRZ(ILSTNR))
         ISYMGD = ISYLRZ(ILSTNR)

C----------------------------------------------------------------------
C 'L0 ' zero-order response Lagrangian multipliers:
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'L0') THEN
         IF (ISIDE .NE.-1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')

Cholesky
         IF (CHOINT .AND. CC2) THEN
            CALL CC_CHOETA(GD,WORK,LWORK)
            NTAMP = NT1AM(ISYMOP)
         ELSE
            CALL CC_ETA(GD,WORK,LWORK)
         END IF

C----------------------------------------------------------------------
C 'N2 ' Lagrangian multipliers:
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'N2') THEN
         IF (ISIDE .NE. -1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')

         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)')
     *           'Finding gradient for '
                 WRITE(LUPRI,'(1x,2(A,F10.6),2(A,I3),A,3I3)')
     *             'Frequencies',FRQIN2(ILSTNR),',',FRQFN2(ILSTNR),
     *             ' and states: ',IIN2(ILSTNR),',',IFN2(ILSTNR),
     *             ' and symmetries: ',ISYIN2(ILSTNR),ISYFN2(ILSTNR),
     *             ISYMGD
         END IF

         CALL CC_RDRSP('BR ',ILSTNR,ISYMGD,IOPTRDG,MODEL,
     &                  GD(IADR1),GD(IADR2))
         IF ( CCR12 ) CALL QUIT('No R12 for BR')

         IF (LOCDBG) THEN
           WRITE(LUPRI,*)'CC_GETGD> RHS vector "BR":',ISYMGD,NTAMP
           CALL OUTPUT(GD(IADR1),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
         END IF

         CALL FLSHFO(LUPRI)

C----------------------------------------------------------------------
C 'M1 ' Lagrangian multipliers:
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'M1') THEN
         IF (ISIDE .NE.-1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')

         IF (IPRINT .GT.2) THEN
           WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for M1'
           WRITE(LUPRI,'(1x,A,F10.6,A,I3,A,I3)')
     *       'Frequency ',FRQLRM(ILSTNR),' and state ',ILRM(ILSTNR),
     *       ' and symmetry ',ISYLRM(ILSTNR)
         END IF

         CALL CC_RDRSP('FR ',ILSTNR,ISYMGD,IOPTRDG,MODEL,
     &                 GD(IADR1),GD(IADR2))
         IF ( CCR12 ) CALL QUIT('No R12 for FR')

         IF (LOCDBG) THEN
           WRITE(LUPRI,*)'CC_GETGD> RHS vector "FR":',ISYMGD,NTAMP
           CALL OUTPUT(GD(IADR1),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
         END IF

         CALL FLSHFO(LUPRI)

C----------------------------------------------------------------------
C 'E0 ' Lagrangian multipliers:
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'E0') THEN
         IF (ISIDE .NE.-1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')

         IF (IPRINT .GT.2) THEN
           WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for E0'
           WRITE(LUPRI,'(1x,A,F10.6,A,I3,A,I3)')
     *       'Frequency ',0.0D0,' and state ',IXGRST(ILSTNR),
     *       ' and symmetry ',1
         END IF

         CALL CC_RDRSP('BE ',ILSTNR,ISYMGD,IOPTRDG,MODEL,
     &                 GD(IADR1),GD(IADR2))
         IF ( CCR12 ) CALL QUIT('No R12 for BE')

         IF (LOCDBG) THEN
           WRITE(LUPRI,*)'CC_GETGD> RHS vector "BE":',ISYMGD,NTAMP
           CALL OUTPUT(GD(IADR1),1,NTAMP,1,1,NTAMP,1,1,LUPRI)
         END IF

         CALL FLSHFO(LUPRI)

C----------------------------------------------------------------------
C 'RC ' first-order right Cauchy vectors:
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'RC') THEN
         IF (ISIDE .NE. 1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(/,1x,A)')
     *           'Finding gradient for '
                 WRITE(LUPRI,'(1x,A,I10,A,A8,A,I3)')
     *             'Right Cauchy Eq.',ILRCAU(ILSTNR),
     *           ' Label: ',LRCLBL(ILSTNR),' Symmetry: ', 
     *            ISYLRC(ILSTNR)
         ENDIF
         CALL CC_CAURHS(GD,ILSTNR,LIST,IOPTRDG,WORK,LWORK)
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)
         ! ---> turn sign for right Cauchy vectors <---
         ! this is done because we want to solve A C(n) = +C(n-1),
         ! thus the gradient vector has the opposite sign that in
         ! the standard formula A t = -X
         CALL DSCAL(NTAMP,XMONE,GD,1)

C----------------------------------------------------------------------
C 'CR2' second-order right Cauchy vectors:
C----------------------------------------------------------------------
      ELSE IF (LIST(1:3) .EQ. 'CR2') THEN
         IF (CC3) CALL QUIT('CR2 vectors not yet implemented for CC3!')
         IF (ISIDE .NE. 1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(/,1x,A)')
     *           'Finding gradient for second-order right Cauchy Eq. '
            WRITE(LUPRI,'(1x,2(A,I3),2(A,A8),A,3I3)')
     *         'Cauchy orders',ICR2CAU(ILSTNR,1),',',ICR2CAU(ILSTNR,2),
     *         ' and labels: ',LBLCR2(ILSTNR,1),',',LBLCR2(ILSTNR,2),
     *         ' and symmetries: ',ISYCR2(ILSTNR,1),ISYCR2(ILSTNR,2),
     *         ISYMGD
         ENDIF
         CALL CC_CR2GD(GD,ILSTNR,WORK,LWORK)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)

C----------------------------------------------------------------------
C 'CL2' second-order left Cauchy vectors:
C----------------------------------------------------------------------
      ELSE IF (LIST(1:3) .EQ. 'CL2') THEN
         IF (CC3) CALL QUIT('CR2 vectors not yet implemented for CC3!')
         IF (ISIDE .NE. -1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 
     *           'Finding gradient for second-order left Cauchy Eq. '
            WRITE(LUPRI,'(1x,2(A,I3),2(A,A8),A,3I3)')
     *         'Cauchy orders',ICL2CAU(ILSTNR,1),',',ICL2CAU(ILSTNR,2),
     *         ' and labels: ',LBLCL2(ILSTNR,1),',',LBLCL2(ILSTNR,2),
     *         ' and symmetries: ',ISYCL2(ILSTNR,1),ISYCL2(ILSTNR,2),
     *         ISYMGD
         END IF
         CALL CC_CL2GD(GD,ILSTNR,WORK,LWORK)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)

C----------------------------------------------------------------------
C 'LC ' first-order left Cauchy vectors:
C----------------------------------------------------------------------
      ELSE IF (LIST(1:2) .EQ. 'LC') THEN
         IF (CC3) CALL QUIT('CR2 vectors not yet implemented for CC3!')
         IF (ISIDE .NE. -1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ISYLC1(ILSTNR)
         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(/,1x,A)')
     *           'Finding gradient for '
                 WRITE(LUPRI,'(1x,A,I10,A,A8,A,I3)')
     *             'Left Cauchy Eq.',ILC1CAU(ILSTNR),
     *           ' Label: ',LBLLC1(ILSTNR),' Symmetry: ', ISYMGD
         ENDIF
         CALL CC_LCGD(GD,ILSTNR,WORK,LWORK)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)

C----------------------------------------------------------------------
C 'ER1' first-order right excited state response vector: --> 'EO1'
C----------------------------------------------------------------------
      ELSE IF (LIST(1:3) .EQ. 'ER1') THEN
         IF (ISIDE .NE. 1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for ER1 response:'
         END IF
         IDLSTEO1 = ILSTNR
         CALL CC_RDRSP('EO1',IDLSTEO1,ISYMGD,IOPTRDG,MODEL,
     &                 GD(IADR1),GD(IADR2))
         CALL CCLR_DIASCL(GD(IADR2),0.5d0,ISYMGD)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)
         IF ( CCR12 ) CALL QUIT('No R12 for EO1')

C----------------------------------------------------------------------
C 'ER2' second-order right excited state response vector: --> 'EO2'
C----------------------------------------------------------------------
      ELSE IF (LIST(1:3) .EQ. 'ER2') THEN
         IF (ISIDE .NE. 1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for ER2 response:'
         END IF
         IDLSTEO2 = ILSTNR
         CALL CC_RDRSP('EO2',IDLSTEO2,ISYMGD,IOPTRDG,MODEL,
     &                 GD(IADR1),GD(IADR2))
         CALL CCLR_DIASCL(GD(IADR2),0.5d0,ISYMGD)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)
         IF ( CCR12 ) CALL QUIT('No R12 for EO2')

C----------------------------------------------------------------------
C 'EL1' first-order left excited state response vector: --> 'EX1'
C----------------------------------------------------------------------
      ELSE IF (LIST(1:3) .EQ. 'EL1') THEN
         IF (ISIDE .NE.-1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for EL1 response:'
         END IF
         IDLSTEX1 = ILSTNR
         CALL CC_RDRSP('EX1',IDLSTEX1,ISYMGD,IOPTRDG,MODEL,
     &                 GD(IADR1),GD(IADR2))
C        CALL CCLR_DIASCL(GD(IADR2),0.5d0,ISYMGD)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)
         IF ( CCR12 ) CALL QUIT('No R12 for EX1')

C----------------------------------------------------------------------
C 'EL2' second-order left excited state response vector: --> 'EX2'
C----------------------------------------------------------------------
      ELSE IF (LIST(1:3) .EQ. 'EL2') THEN
         IF (ISIDE .NE.-1)
     &        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         ISYMGD = ILSTSYM(LIST,ILSTNR)
         IF (IPRINT .GT.2) THEN
            WRITE(LUPRI,'(/,1x,A)') 'Finding gradient for EL2 response:'
         END IF
         IDLSTEX2 = ILSTNR
         CALL CC_RDRSP('EX2',IDLSTEX2,ISYMGD,IOPTRDG,MODEL,
     &                 GD(IADR1),GD(IADR2))
C        CALL CCLR_DIASCL(GD(IADR2),0.5d0,ISYMGD)
         NTAMP = NT1AM(ISYMGD) + NT2AM(ISYMGD)
         IF ( CCS ) NTAMP = NT1AM(ISYMGD)
         IF ( CCR12 ) CALL QUIT('No R12 for EX2')

C----------------------------------------------------------------------
C 'PL1 ' projected 1st-order respons Lagrang. multipls: --> 'X1 '+'F1 '
C----------------------------------------------------------------------

      ELSE IF (LIST(1:3) .EQ. 'PL1') THEN
         IF (ISIDE .NE.-1) 
     *        CALL QUIT('Mismatch ISIDE and LIST in CC_GETGD')
         IF (IPRINT .GT. 2) THEN
            WRITE(LUPRI,'(/,1x,A)')
     *           'Finding gradient for '
                 WRITE(LUPRI,'(1x,A,F10.6,A,A8,A,I3)')
     *             'Frequency',FRQPL1(ILSTNR),
     *           ' and label: ',LBLPL1(ILSTNR),' and symmetry: ',
     *            ISYPL1(ILSTNR)
         ENDIF
         CALL CC_LGD(GD,LBLPL1(ILSTNR),ISYPL1(ILSTNR),
     *               FRQPL1(ILSTNR),LORXPL1(ILSTNR),IOPTRDG,WORK,LWORK)
         NTAMP = NT1AM(ISYPL1(ILSTNR)) + NT2AM(ISYPL1(ILSTNR))
         IF ( CCS ) NTAMP = NT1AM(ISYPL1(ILSTNR))
         IF ( CCR12 ) CALL QUIT('No R12 for X1+F1 as used for PL1')
         ISYMGD = ISYPL1(ILSTNR)
C----------------------------------------------------------------------
      ELSE
         WRITE (LUPRI,*) 'Unknown list in CC_GETGD:',LIST(1:3)
         CALL QUIT('Unknown list in CC_GETGD.')
      ENDIF
C
      IF (CCSTST) THEN
         CALL DZERO(GD(1+NT1AM(ISYMGD)),NTAMP-NT1AM(ISYMGD))
      ENDIF
C
      IF (NTAMP.LE.0) 
     &  WRITE (LUPRI,*) 'strange number of amplitudes:',NTAMP

      IF (DEBUG .OR. LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_GETGD> LIST:',LIST
         WRITE (LUPRI,*) 'CC_GETGD> NTAMP:',NTAMP
         XGD   = DDOT(NTAMP,GD,1,GD,1)
         WRITE(LUPRI,'(1X,A,1X,E20.10)') 
     &      'CC_GETGD: Norm of GD vector:       ',XGD
         IF ( XGD .LT. 1.0D-14) THEN
            WRITE(LUPRI,'(1X,A)') 'WARNING: Zero gradient ???? '
         ENDIF
         CALL FLSHFO(LUPRI)
      ENDIF
C
      CALL QEXIT('CC_GETGD')
      RETURN
      END
c*DECK CC_LGD 
      SUBROUTINE CC_LGD(GD,LBLC,ISYMC,FRQ,LORX,IOPTR,WORK,LWORK)
C
C----------------------------------------------------------------------
C
C     Purpose: Calculate GD vector - right hand side for left response
C              equations.
C
C     Written by Ove Christiansen 17-10-1996
C     Changes for orbital relaxed response, Christof Haettig spring '99
C     Changes for CC3 response, Christof Haettig sping '02
C     Adapted for CC-R12 response, Christian Neiss April 2005
C
C----------------------------------------------------------------------
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccorb.h"
#include "iratdef.h"
#include "cclr.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
#include "ccl1rsp.h"
#include "ccr1rsp.h"
#include "ccroper.h"
#include "dummy.h"
C
      PARAMETER( ONE = 1.0D00 )
      DIMENSION GD(*),WORK(LWORK)
      LOGICAL LORX
      CHARACTER LBLC*(*),MODEL*10
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND( 'IN CC_LGD: Constructing left GD vector ')
      ENDIF
C
      TIMEC = SECOND()
C
C-----------------
C     Allocations.
C-----------------
C
      NTAMP = NT1AM(ISYMC) + NT2AM(ISYMC)
      IF ( CCS ) NTAMP = NT1AM(ISYMC)
      IF ( CCR12 ) NTAMP = NTAMP + NTR12AM(ISYMC)
C
      IF (LWORK .LT. 0)
     *      CALL QUIT('Insufficient space for allocation in CC_LGD -1')
C
C-----------------------------------------
C     get eta^C contribution.
C-----------------------------------------
C
      ! find operator index
      ISYM = ISYMC
      IOPER = IROPER(LBLC,ISYM)
      IF (ISYM.NE.ISYMC) CALL QUIT('Symmetry mismatch in CC_LGD.')
C
C     IF ( LORX .OR. LPDBSOP(IOPER) .OR. CCSDT) THEN
C       ! if the orbitals are allowed to relax in the field or if the
C       ! basis set depends on the perturbation, read ETA from file,
C       ! also if we the 'effective' ETA vector for a triples model
C       ! read this vector from file
        ILSTETA = IETA1(LBLC,LORX,FRQ,ISYMC)
        IOPT    = IOPTR
        CALL CC_RDRSP('X1 ',ILSTETA,ISYMC,IOPT,MODEL,GD(1),
     *                GD(1+NT1AM(ISYMC)))
C     ELSE
C       ! if it is a simple unrelaxed one-electron perturbation
C       ! calculate the ETA vector in CC_ETAC
C       CALL CC_ETAC(ISYMC,LBLC,GD,'L0 ',1,0,DUMMY,WORK,LWORK)
C     END IF
C
      IF (CCR12) THEN
        IOPT = 32
        CALL CC_RDRSP('X1 ',ILSTETA,ISYMC,IOPT,MODEL,DUMMY,
     *                GD(1+NTAMP-NTR12AM(ISYMC)))
      END IF
C
      IF ( DEBUG ) THEN
         X    = DDOT(NTAMP,GD,1,GD,1)
         WRITE(LUPRI,1) 'Norm of GD after eta cont.:     ',X
      ENDIF
C
      IF (CIS) RETURN
C
C-----------------------------------------------------------
C     F-contribution section.
C-----------------------------------------------------------
C
      KFRV  = 1
      KEND1 = KFRV  + NTAMP
      LEND1 = LWORK - KEND1
C
      IF (LEND1 .LT. NTAMP)
     *      CALL QUIT('Insufficient space for allocation in CC_LGD -2')
C
      KFRV1 = KFRV  
      KFRV2 = KFRV  + NT1AM(ISYMC)
      KFRVR12 = KEND1 - NTR12AM(ISYMC) 
C
C-----------------------------------------------
C     Calculate list number.
C
C     Find right response vector:
C
C     Symmetry - the same
C     Label    - the same
C     Frequency - the same with opposite sign!!
C-----------------------------------------------
C
      ISAVE = NLRTLBL
      ILTNR = IR1TAMP(LBLC,LORX,FRQ,ISYMC)
      IF (ILTNR .GT. NLRTLBL) THEN
         WRITE(LUPRI,*) 'In CC_LGD: ILTNR, old NLRTLBL',ILTNR,ISAVE
         CALL QUIT('CC_LGD: Did not find t-response vector ')
      ENDIF

      IOPT = IOPTR 
      CALL CC_RDRSP('F1 ',ILTNR,ISYMC,IOPT,MODEL,WORK(KFRV1),
     *              WORK(KFRV2))
C
      IF (CCR12) THEN
        IOPT = 32
        CALL CC_RDRSP('F1 ',ILTNR,ISYMC,IOPT,MODEL,DUMMY,
     *              WORK(KFRVR12))
      END IF
C
      IF ( DEBUG ) THEN
         X = DDOT(NTAMP,WORK(KFRV),1,WORK(KFRV),1)
         WRITE(LUPRI,1) 'CC_LGD: Norm of F*Tx vector:       ',X
         X    = DDOT(NTAMP,GD,1,GD,1)
         WRITE(LUPRI,1) 'CC_LGD: Norm of GD -XXXXXXXX       ',X
      ENDIF
C
C-------------------------
C     3.Add to the vector.
C-------------------------
C
      CALL DAXPY(NTAMP,ONE,WORK(KFRV),1,GD,1)
C
      IF ( DEBUG ) THEN
         X    = DDOT(NTAMP,GD,1,GD,1)
         WRITE(LUPRI,1) 'Norm of left GD after F*Tx cont.:  ',X
      ENDIF
C
   1  FORMAT(1x,A35,1X,E20.10)
      END
c*DECK CC_CAURHS
      SUBROUTINE CC_CAURHS(GD,ILSTNR,LIST,IOPTRDG,WORK,LWORK)
C
C-----------------------------------------------------------------------
C
C     Purpose: Calculate GD vector - right hand side for Cauchy vector
C              equations.
C
C     Written by Ove Christiansen 8-dec.-1996
C
C-----------------------------------------------------------------------
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccorb.h"
#include "iratdef.h"
#include "cclr.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccrc1rsp.h"
C
      PARAMETER( ONE = 1.0D00 )
      DIMENSION GD(*),WORK(LWORK)
      CHARACTER MODEL1*10,LIST*2
C
      IF ( IPRINT .GT. 10 ) THEN
         CALL AROUND('IN CC_CAURHS:Constructing Cauchy right hand side')
      ENDIF
C
      ISYM = ISYLRC(ILSTNR)
      KT1 = 1
      KT2 = KT1 + NT1AM(ISYM)
C
      IF (ILRCAU(ILSTNR).EQ.1) THEN
C
C-----------------------------------------------------------
C       Get first order response vectors for zero frequency.
C-----------------------------------------------------------
C
        FRQ   = 0.0D0
        INUM  = IR1TAMP(LRCLBL(ILSTNR),.FALSE.,FRQ,ISYM)
        CALL CC_RDRSP('R1 ',INUM,ISYM,IOPTRDG,MODEL1,GD(KT1),GD(KT2))
C
      ELSE 
C
C---------------------------
C       Cauchy (k-1) vector.
C---------------------------
C
        KPREV = ILRCAU(ILSTNR)-1
        INUM  = ILRCAMP(LRCLBL(ILSTNR),KPREV,ISYM) 
        CALL CC_RDRSP('RC ',INUM,ISYM,IOPTRDG,MODEL1,GD(KT1),GD(KT2))
C
      ENDIF
      TIMEC = SECOND()
      END
*---------------------------------------------------------------------*
c /* deck cc_r2gd */
*=====================================================================*
      SUBROUTINE CC_R2GD(GD,IDLSTR2,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: build right hand side vector for the second order
*             response equations for the cluster amplitudes (R2)
*
*    input/output:   GD       -- right hand side vector (output)
*                    IDLSTR2  -- list index of the corresponding
*                                right hand side vector
*
*    Written by Christof Haettig, Februar 1997.
*    Adapted for CC-R12, Christian Neiss June 2005
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccroper.h"
#include "ccr1rsp.h"
#include "ccr2rsp.h"
#include "cco2rsp.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IDLSTR2, LWORK

#if defined (SYS_CRAY)
      REAL GD(*), WORK(LWORK)
      REAL DDOT, DUMMY
      REAL ONE
      REAL FREQA, FREQB
#else
      DOUBLE PRECISION GD(*), WORK(LWORK)
      DOUBLE PRECISION DDOT, DUMMY
      DOUBLE PRECISION ONE
      DOUBLE PRECISION FREQA, FREQB
#endif
      PARAMETER (ONE = 1.0d0)

      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABELA, LABELB
      INTEGER ISYMR2, ISYMA, ISYMB, IADR1, IADR2, NTAMP, IADR12
      INTEGER IDLSTR1A, IDLSTR1B, IOPT, ISYMO2, IDLSTO2

* external functions:
      INTEGER IR1TAMP
      INTEGER IRHSR2

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      ISYMA  = ISYAR2T(IDLSTR2)
      ISYMB  = ISYBR2T(IDLSTR2)
      ISYMR2 = MULD2H(ISYMA,ISYMB)

      LABELA = LBLAR2T(IDLSTR2)
      LABELB = LBLBR2T(IDLSTR2)

      FREQA  = FRQAR2T(IDLSTR2)
      FREQB  = FRQBR2T(IDLSTR2)

      IDLSTR1A = IR1TAMP(LABELA,.FALSE.,FREQA,ISYMA)
      IDLSTR1B = IR1TAMP(LABELB,.FALSE.,FREQB,ISYMB)
      IDLSTO2  = IRHSR2(LABELA,.FALSE.,FREQA,ISYMA,
     &                  LABELB,.FALSE.,FREQB,ISYMB)
      ISYMO2   = MULD2H(ISYMA,ISYMB)

      IF (IDLSTO2 .EQ. -1) THEN
        WRITE (LUPRI,*) ' A requested O2 vector is not available.'
        CALL QUIT(' A requested O2 vector was not available.')
      END IF

      IADR1 = 1
      IADR2 = IADR1 + NT1AM(ISYMR2)
      NTAMP = NT1AM(ISYMR2) + NT2AM(ISYMR2)
      IF (CCR12) THEN
        IADR12 = IADR2 + NT2AM(ISYMR2)
        NTAMP = NTAMP + NTR12AM(ISYMR2)
      END IF

      IF (LWORK .LT. NTAMP) THEN
        CALL QUIT('Insufficient memory in CC_GDR2.')
      END IF

*---------------------------------------------------------------------*
* 1. contribution: B * R1(A) * R1(B), should be stored on file        *
*---------------------------------------------------------------------*
      IOPT = 3
      IF (CCSDT) IOPT = 24
      CALL CC_RDRSP('O2',IDLSTO2,ISYMO2,IOPT,MODEL,GD(IADR1),GD(IADR2))
      CALL CCLR_DIASCL(GD(IADR2),0.5D0,ISYMO2)

      IF (CCR12) THEN
        IOPT = 32
        CALL CC_RDRSP('O2',IDLSTO2,ISYMO2,IOPT,MODEL,DUMMY,GD(IADR12))
      END IF

*---------------------------------------------------------------------*
* write debug output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_GDR2> norm^2 of GD-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_GDR2> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYMR2),GD(IADR1),1,GD(IADR1),1)
        WRITE (LUPRI,*) 'DEBUG_CC_GDR2> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYMR2),GD(IADR2),1,GD(IADR2),1)
        IF (CCR12) THEN
          WRITE (LUPRI,*) 'DEBUG_CC_GDR2> norm^2 of r12 doubles part:',
     &       DDOT(NTR12AM(ISYMR2),GD(IADR12),1,GD(IADR12),1)
        END IF
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_R2GD                              *
*---------------------------------------------------------------------*
c /* deck CC_L2GD */
*=====================================================================*
      SUBROUTINE CC_L2GD(GD,IDLSTL2,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: build right hand side vector for the second-order
*             response equations for the lagrange multipliers (L2)
*
*    input/output:   GD       -- right hand side vector (output)
*                    IDLSTL2  -- list index of the corresponding
*                                response vector
*
*    Written by Christof Haettig, April 1997.
*    Adapted for CC-R12, Christian Neiss June 2005
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccr2rsp.h"
#include "ccl2rsp.h"
#include "ccx2rsp.h"
#include "dummy.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IDLSTL2, LWORK

#if defined (SYS_CRAY)
      REAL GD(*), WORK(LWORK)
      REAL DDOT
      REAL ONE
      REAL FREQA, FREQB
#else
      DOUBLE PRECISION GD(*), WORK(LWORK)
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION ONE
      DOUBLE PRECISION FREQA, FREQB
#endif
      PARAMETER (ONE = 1.0d0)

      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABELA, LABELB
      INTEGER ISYML2, ISYMA, ISYMB, IADR1, IADR2, NTAMP, IADR12
      INTEGER IOPT, IDLSTR2, IDLSTX2

* external functions:
      INTEGER ICHI2
      INTEGER IR2TAMP

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      ISYMA  = ISYAL2(IDLSTL2)
      ISYMB  = ISYBL2(IDLSTL2)
      ISYML2 = MULD2H(ISYMA,ISYMB)

      LABELA = LBLAL2(IDLSTL2)
      LABELB = LBLBL2(IDLSTL2)

      FREQA  = FRQAL2(IDLSTL2)
      FREQB  = FRQBL2(IDLSTL2)

      IDLSTR2 = IR2TAMP(LABELA,.FALSE.,FREQA,ISYMA,
     &                  LABELB,.FALSE.,FREQB,ISYMB)
      IDLSTX2 = ICHI2(LABELA,.FALSE.,FREQA,ISYMA,
     &                LABELB,.FALSE.,FREQB,ISYMB)

      IADR1 = 1
      IADR2 = IADR1 + NT1AM(ISYML2)
      NTAMP = NT1AM(ISYML2) + NT2AM(ISYML2)
      IF (CCS) NTAMP = NT1AM(ISYML2)
      IF (CCR12) THEN
        IADR12 = IADR2 + NT2AM(ISYML2)
        NTAMP  = NTAMP + NTR12AM(ISYML2)
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'LABELA, LABELB:',LABELA,LABELB
        WRITE (LUPRI,*) 'FREQA, FREQB:',FREQA,FREQB
        WRITE (LUPRI,*) 'IDLSTL2,IDLSTX2,IDLSTR2:',IDLSTL2,IDLSTX2,
     &                   IDLSTR2
      END IF

      IF (LWORK .LT. NTAMP) THEN
        CALL QUIT('Insufficient memory in CC_GDL2.')
      END IF

*---------------------------------------------------------------------*
* 1. contribution: second-order chi vector
*---------------------------------------------------------------------*
      IOPT = 3
      CALL CC_RDRSP('X2',IDLSTX2,ISYML2,IOPT,MODEL,GD(IADR1),GD(IADR2))
      IF (CCR12) THEN
        IOPT = 32
        CALL CC_RDRSP('X2',IDLSTX2,ISYML2,IOPT,MODEL,DUMMY,GD(IADR12))
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_GDL2> norm^2 of X2-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_GDL2> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYML2),GD(IADR1),1,GD(IADR1),1)
        IF (.NOT.CCS) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_GDL2> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYML2),GD(IADR2),1,GD(IADR2),1)
        IF (CCR12) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_GDL2> norm^2 of r12 doubles part:',
     &     DDOT(NTR12AM(ISYML2),GD(IADR12),1,GD(IADR12),1)
      END IF
*---------------------------------------------------------------------*
* 2. contribution: F matrix transformed R2 vector
*---------------------------------------------------------------------*
      IOPT = 3
      CALL CC_RDRSP('F2',IDLSTR2,ISYML2,IOPT,MODEL,
     &              WORK(IADR1),WORK(IADR2))
      IF (CCR12) THEN
        IOPT = 32
        CALL CC_RDRSP('F2',IDLSTR2,ISYML2,IOPT,MODEL,DUMMY,WORK(IADR12))
      END IF

      CALL DAXPY(NTAMP,ONE,WORK,1,GD,1)
      
*---------------------------------------------------------------------*
* write debug output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_GDL2> norm^2 of GD-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_GDL2> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYML2),GD(IADR1),1,GD(IADR1),1)
        IF (.NOT.CCS) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_GDL2> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYML2),GD(IADR2),1,GD(IADR2),1)
        IF (CCR12)
     &   WRITE (LUPRI,*) 'DEBUG_CC_GDL2> norm^2 of r12 doubles part:',
     &     DDOT(NTR12AM(ISYML2),GD(IADR12),1,GD(IADR12),1)
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_L2GD                              *
c /* deck CC_L3GD */
*=====================================================================*
      SUBROUTINE CC_L3GD(GD,IDLSTL3,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: build right hand side vector for the third-order
*             response equations for the lagrange multipliers (L3)
*
*    input/output:   GD       -- right hand side vector (output)
*                    IDLSTL3  -- list index of the corresponding
*                                response vector
*
*    Written by Christof Haettig, April 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccr3rsp.h"
#include "ccl3rsp.h"
#include "ccx3rsp.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IDLSTL3, LWORK

#if defined (SYS_CRAY)
      REAL GD(*), WORK(LWORK)
      REAL DDOT
      REAL ONE
      REAL FREQA, FREQB, FREQC
#else
      DOUBLE PRECISION GD(*), WORK(LWORK)
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION ONE
      DOUBLE PRECISION FREQA, FREQB, FREQC
#endif
      PARAMETER (ONE = 1.0d0)

      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABELA, LABELB, LABELC
      INTEGER ISYML3, ISYMA, ISYMB, ISYMC, IADR1, IADR2, NTAMP
      INTEGER IOPT, IDLSTR3, IDLSTX3

* external functions:
      INTEGER ICHI3
      INTEGER IR3TAMP
      INTEGER ILSTSYM

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      ISYML3 = ILSTSYM('L3',IDLSTL3)

      LABELA = LBLL3(IDLSTL3,1)
      LABELB = LBLL3(IDLSTL3,2)
      LABELC = LBLL3(IDLSTL3,3)

      FREQA  = FRQL3(IDLSTL3,1)
      FREQB  = FRQL3(IDLSTL3,2)
      FREQC  = FRQL3(IDLSTL3,3)

      IDLSTR3 = IR3TAMP(LABELA,FREQA,ISYMA,LABELB,FREQB,ISYMB,
     &                                     LABELC,FREQC,ISYMC)
      IDLSTX3 = ICHI3(LABELA,FREQA,ISYMA,LABELB,FREQB,ISYMB,
     &                                   LABELC,FREQC,ISYMC)

      IADR1 = 1
      IADR2 = IADR1 + NT1AM(ISYML3)
      NTAMP = NT1AM(ISYML3) + NT2AM(ISYML3)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'LABELA, LABELB:',LABELA,LABELB
        WRITE (LUPRI,*) 'FREQA, FREQB:',FREQA,FREQB
        WRITE (LUPRI,*) 'IDLSTL3,IDLSTX3,IDLSTR3:',IDLSTL3,IDLSTX3,
     &                   IDLSTR3
      END IF

      IF (LWORK .LT. NTAMP) THEN
        CALL QUIT('Insufficient memory in CC_GDL3.')
      END IF

*---------------------------------------------------------------------*
* 1. contribution: third-order chi vector
*---------------------------------------------------------------------*
      IOPT = 3
      CALL CC_RDRSP('X3',IDLSTX3,ISYML3,IOPT,MODEL,GD(IADR1),GD(IADR2))

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_GDL3> norm^2 of X3-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_GDL3> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYML3),GD(IADR1),1,GD(IADR1),1)
        WRITE (LUPRI,*) 'DEBUG_CC_GDL3> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYML3),GD(IADR2),1,GD(IADR2),1)
      END IF
*---------------------------------------------------------------------*
* 2. contribution: F matrix transformed R3 vector
*---------------------------------------------------------------------*
      IOPT = 3
      CALL CC_RDRSP('F3',IDLSTR3,ISYML3,IOPT,MODEL,
     &              WORK(IADR1),WORK(IADR2))

      CALL DAXPY(NTAMP,ONE,WORK,1,GD,1)
      
*---------------------------------------------------------------------*
* write debug output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_GDL3> norm^2 of GD-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_GDL3> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYML3),GD(IADR1),1,GD(IADR1),1)
        WRITE (LUPRI,*) 'DEBUG_CC_GDL3> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYML3),GD(IADR2),1,GD(IADR2),1)
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_L3GD                              *
c /* deck CC_LCGD */
*=====================================================================*
      SUBROUTINE CC_LCGD(GD,IDLSTLC,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: build right hand side vector for the first-order
*             left Cauchy equations (LC)
*
*    input/output:   GD       -- right hand side vector (output)
*                    IDLSTLC  -- list index of the corresponding
*                                Cauchy vector
*
*    Written by Christof Haettig, October 1997.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccrc1rsp.h"
#include "cclc1rsp.h"
#include "dummy.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER IDLSTLC, LWORK

#if defined (SYS_CRAY)
      REAL GD(*), WORK(LWORK)
      REAL DDOT
      REAL ONE
#else
      DOUBLE PRECISION GD(*), WORK(LWORK)
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION ONE
#endif
      PARAMETER (ONE = 1.0d0)

      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABEL
      INTEGER ISYM, IADR1, IADR2, NTAMP, ICAU, ISYM1
      INTEGER IOPT, IDLSTRC, IDLSTMC 

* external functions:
      INTEGER IR1TAMP
      INTEGER IL1ZETA
      INTEGER ILRCAMP
      INTEGER ILC1AMP

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      ISYM   = ISYLC1(IDLSTLC)
      LABEL  = LBLLC1(IDLSTLC)
      ICAU   = ILC1CAU(IDLSTLC)

      IF (ICAU .LT. 0) THEN
        WRITE (LUPRI,*) 'Error in CC_LCGD> rhs vectors for Cauchy ',
     &          'orders < 0 not implemented!'
        CALL QUIT('Error in CC_LCGD: illegal Cauchy order.')
      END IF

      IDLSTRC = ILRCAMP(LABEL,ICAU,ISYM1)
      IF (ISYM1.NE.ISYM) THEN
        CALL QUIT('Symmetry mismatch in CC_LCGD.')
      END IF

      IADR1 = 1
      IADR2 = IADR1 + NT1AM(ISYM)
      NTAMP = NT1AM(ISYM) + NT2AM(ISYM)
      IF (CCS) NTAMP = NT1AM(ISYM)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'LABEL:',LABEL
        WRITE (LUPRI,*) 'ICAU:',ICAU
        WRITE (LUPRI,*) 'IDLSTLC,IDLSTRC:',IDLSTLC,IDLSTRC
        WRITE (LUPRI,*) 'NTAMP:',NTAMP
        WRITE (LUPRI,*) 'LWORK:',LWORK
      END IF

      IF (LWORK .LT. NTAMP) THEN
        CALL QUIT('Insufficient memory in CC_LCGD.')
      END IF

*---------------------------------------------------------------------*
* 1. contribution: Left Cauchy vector of order ICAU-1
*---------------------------------------------------------------------*
      IOPT = 3
      IF (ICAU.GT.1) THEN
        IDLSTMC = ILC1AMP(LABEL,ICAU-1,  ISYM)
        CALL CC_RDRSP('LC',IDLSTMC,ISYM,IOPT,MODEL,GD(IADR1),GD(IADR2))
      ELSE IF (ICAU.EQ.1) THEN
        IDLSTMC = IL1ZETA(LABEL,.FALSE.,0.0d0,ISYM)
        CALL CC_RDRSP('L1',IDLSTMC,ISYM,IOPT,MODEL,GD(IADR1),GD(IADR2))
      ELSE IF (ICAU.EQ.0) THEN
        CALL CC_ETAC(ISYM,LABEL,GD,'L0',1,0,DUMMY,WORK,LWORK)
      ELSE
        CALL QUIT('Error in CC_LCGD.')
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_LCGD> norm^2 of LC(ICAU-1)-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_LCGD> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYM),GD(IADR1),1,GD(IADR1),1)
        IF (.NOT.CCS) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_LCGD> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYM),GD(IADR2),1,GD(IADR2),1)
      END IF

*---------------------------------------------------------------------*
* 2. contribution: F matrix transformed RC vector
*---------------------------------------------------------------------*
      IF (ICAU.GT.0) THEN
        IOPT = 3
        CALL CC_RDRSP('FC',IDLSTRC,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))
      ELSE IF (ICAU.EQ.0) THEN
        IOPT = 3
        IDLSTRC = IR1TAMP(LABEL,.FALSE.,0.0d0,ISYM)
        CALL CC_RDRSP('F1',IDLSTRC,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))
      ELSE
        CALL QUIT('Error in CC_LCGD.')
      END IF

      CALL DAXPY(NTAMP,ONE,WORK(IADR1),1,GD,1)
      
*---------------------------------------------------------------------*
* write debug output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_LCGD> norm^2 of GD-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_LCGD> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYM),GD(IADR1),1,GD(IADR1),1)
        IF (.NOT.CCS) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_LCGD> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYM),GD(IADR2),1,GD(IADR2),1)
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_L2GD                              *
*---------------------------------------------------------------------*
c /* deck cc_fre  */
      SUBROUTINE CC_FRE(WORK,LWORK)
C
C     Ove Christiansen March/april-1997
C
C     Transform eigen-vectors with F-matrix and write to file.
C
C     ISYMV is symmetry of vector and transformed vector.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "dummy.h"
#include "iratdef.h"
#include "cclr.h"
#include "ccexci.h"
#include "ccexgr.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsections.h"
#include "ccsdsym.h"
C
      PARAMETER( TOLFRQ=1.0D-08 )
      PARAMETER( ONE = 1.0D00 , XMONE = -1.0D0)
C
      DIMENSION WORK(LWORK)
      CHARACTER MODEL*10,RSLIST*2
      LOGICAL INCORE
C
#include "leinf.h"
C
      INCORE = .TRUE.
      MODEL           = '----------'
      IF (CC2)  MODEL = 'CC2       '
      IF (CCS)  MODEL = 'CCS       '
      IF (CCSD) MODEL = 'CCSD      '
      IF (CC3)  MODEL = 'CC3       '
      IF (MODEL(1:1).EQ.'-') CALL QUIT('Unknown model in CC_FRE.')
C
      NTAMP = NT1AMX + NT2AMX
      IF (CCS) NTAMP = NT1AMX
C
C------------------
C     Make vectors.
C------------------
C
      DO 1000 IL = 1,NXGRST
         IE = IXGRST(IL) 
CCH
         CALL CC_FTRAN('LE ',IE,'RE ',IE,WORK,LWORK)
CCH      CALL CC_FGL('LE',IE,'RE',IE,'--',0,.TRUE.,WORK,LWORK)
CCH
         IF (.NOT. (CCSLV)) THEN
           KVEC  = 1
           KGD   = KVEC  + NTAMP
           KEND1 = KGD   + NTAMP
           LWRK1 = LWORK - KEND1
           CALL CC_ETA(WORK(KGD),WORK(KEND1),LWRK1)
           CALL DAXPY(NTAMP,ONE,WORK(KGD),1,WORK(KVEC),1)
         END IF
C
         K1    = 1
         K2    = K1 + NT1AMX
         KEND1 = K2 + NT2AMX
         LWRK1 = LWORK - KEND1
         ILSTB = IL 
         IOPT = 3
         CALL CC_WRRSP('BE ',ILSTB,ISYMOP,IOPT,MODEL,DUMMY,
     &                 WORK(K1),WORK(K2),WORK(KEND1),LWRK1)
 1000 CONTINUE
C
      END
c /* deck CC_CR2GD */
*=====================================================================*
      SUBROUTINE CC_CR2GD(GD,IDLSTCR2,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: build right hand side vector for the second-order
*             right Cauchy equations (CR2)
*
*    input/output:   GD        -- right hand side vector (output)
*                    IDLSTCR2  -- list index of the corresponding
*                                 Cauchy vector
*
*    Written by Christof Haettig, March 1998.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccrc1rsp.h"
#include "cccr2rsp.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL GD(*), WORK(LWORK)
      REAL DDOT
      REAL MONE, ZERO, ONE
#else
      DOUBLE PRECISION GD(*), WORK(LWORK)
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION MONE, ZERO, ONE
#endif
      PARAMETER (MONE = -1.0d0, ZERO = 0.0d0, ONE = 1.0d0)

      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABEL1, LABEL2
      INTEGER ISYM, IADR1, IADR2, NTAMP, ICAU1, ICAU2, ISYM1, ISYM2
      INTEGER IOPT, IDLSTCR2, IDLSTO2, IDLSTMC1, IDLSTMC2

* external functions:
      INTEGER ILRCAMP
      INTEGER IR1TAMP
      INTEGER ICR2AMP
      INTEGER IR2TAMP
      INTEGER IRHSCR2
      INTEGER IRHSR2

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      ISYM1   = ISYCR2(IDLSTCR2,1)
      ISYM2   = ISYCR2(IDLSTCR2,2)
      ISYM    = MULD2H(ISYM1,ISYM2)

      LABEL1  = LBLCR2(IDLSTCR2,1)
      LABEL2  = LBLCR2(IDLSTCR2,2)

      ICAU1   = ICR2CAU(IDLSTCR2,1)
      ICAU2   = ICR2CAU(IDLSTCR2,2)

      IF (ICAU1.LT.0 .OR. ICAU2.LT.0) THEN
        WRITE (LUPRI,*) 'Error in CC_CR2GD> rhs vectors for Cauchy ',
     &          'orders < 0 not implemented!'
        CALL QUIT('Error in CC_CR2GD: illegal Cauchy order.')
      END IF

      IADR1 = 1
      IADR2 = IADR1 + NT1AM(ISYM)
      NTAMP = NT1AM(ISYM) + NT2AM(ISYM)
      IF (CCS) NTAMP = NT1AM(ISYM)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'IDLSTCR2:',IDLSTCR2
        WRITE (LUPRI,*) 'LABELs:',LABEL1, LABEL2
        WRITE (LUPRI,*) 'ICAUs:',ICAU1, ICAU2
        WRITE (LUPRI,*) 'NTAMP:',NTAMP
        WRITE (LUPRI,*) 'LWORK:',LWORK
      END IF

      IF (LWORK .LT. NTAMP) THEN
        WRITE (LUPRI,*) 'Insufficient memory in CC_CR2GD:'
        WRITE (LUPRI,*) 'available:',LWORK, '    needed:',NTAMP
        WRITE (LUPRI,*) 'for single part:',NT1AM(ISYM)
        WRITE (LUPRI,*) 'for doubles part:',NT2AM(ISYM)
        CALL QUIT('Insufficient memory in CC_CR2GD.')
      END IF

*---------------------------------------------------------------------*
* 1. contribution: CO2 vector
*---------------------------------------------------------------------*
      IOPT = 3
      IF ((ICAU1+ICAU2).EQ.0) THEN
        IDLSTO2 = IRHSR2(LABEL1,.FALSE.,ZERO,ISYM1,
     &                   LABEL2,.FALSE.,ZERO,ISYM2)
        CALL CC_RDRSP('O2',IDLSTO2,ISYM,IOPT,MODEL,GD(IADR1),GD(IADR2))
      ELSE IF ((ICAU1+ICAU2).GT.0) THEN
        IDLSTO2 = IRHSCR2(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2,ISYM2)
        CALL CC_RDRSP('CO2',IDLSTO2,ISYM,IOPT,
     &                MODEL,GD(IADR1),GD(IADR2))
      END IF
      CALL CCLR_DIASCL(GD(IADR2),0.5d0,ISYM)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_CR2GD> norm^2 of CO2-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_CR2GD> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYM),GD(IADR1),1,GD(IADR1),1)
        IF (.NOT.CCS) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_CR2GD> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYM),GD(IADR2),1,GD(IADR2),1)
      END IF

*---------------------------------------------------------------------*
* 2. contribution: CR2(ICAU1-1,ICAU2) vector
*---------------------------------------------------------------------*
      IF (ICAU1.GT.1 .OR. (ICAU1.EQ.1 .AND. ICAU2.GT.0) ) THEN

        IDLSTMC2 = ICR2AMP(LABEL1,ICAU1-1,ISYM1,LABEL2,ICAU2,ISYM2)

        IOPT = 3
        CALL CC_RDRSP('CR2',IDLSTMC2,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))

        CALL DAXPY(NTAMP,MONE,WORK(IADR1),1,GD,1)

      ELSE IF (ICAU1.EQ.1 .AND. ICAU2.EQ.0) THEN

        IDLSTMC2 = IR2TAMP(LABEL1,.FALSE.,ZERO,ISYM1,
     &                     LABEL2,.FALSE.,ZERO,ISYM2)

        IOPT = 3
        CALL CC_RDRSP('R2',IDLSTMC2,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))

        CALL DAXPY(NTAMP,MONE,WORK(IADR1),1,GD,1)

      ELSE IF (ICAU1.EQ.0) THEN
        CONTINUE
      ELSE
        CALL QUIT('ERROR IN CC_CR2GD.')
      END IF
      
*---------------------------------------------------------------------*
* 3. contribution: CR2(ICAU1,ICAU2-1) vector
*---------------------------------------------------------------------*
      IF (ICAU2.GT.1 .OR. (ICAU2.EQ.1 .AND. ICAU1.GT.0) ) THEN
        IDLSTMC2 = ICR2AMP(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2-1,ISYM2)

        IOPT = 3
        CALL CC_RDRSP('CR2',IDLSTMC2,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))

        CALL DAXPY(NTAMP,MONE,WORK(IADR1),1,GD,1)

      ELSE IF (ICAU2.EQ.1 .AND. ICAU1.EQ.0) THEN

        IDLSTMC2 = IR2TAMP(LABEL1,.FALSE.,ZERO,ISYM1,
     &                     LABEL2,.FALSE.,ZERO,ISYM2)

        IOPT = 3
        CALL CC_RDRSP('R2',IDLSTMC2,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))

        CALL DAXPY(NTAMP,MONE,WORK(IADR1),1,GD,1)

      ELSE IF (ICAU2.EQ.0) THEN
        CONTINUE
      ELSE
        CALL QUIT('ERROR IN CC_CR2GD.')
      END IF
      
*---------------------------------------------------------------------*
* write debug output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_CR2GD> norm^2 of GD-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_CR2GD> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYM),GD(IADR1),1,GD(IADR1),1)
        IF (.NOT.CCS) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_CR2GD> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYM),GD(IADR2),1,GD(IADR2),1)
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_CR2GD                             *
*---------------------------------------------------------------------*
c /* deck CC_CL2GD */
*=====================================================================*
      SUBROUTINE CC_CL2GD(GD,IDLSTCL2,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: build right hand side vector for the second-order
*             left Cauchy equations (CL2)
*
*    input/output:   GD        -- right hand side vector (output)
*                    IDLSTCL2  -- list index of the corresponding
*                                 Cauchy vector
*
*    Written by Christof Haettig, March 1998.
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "cccl2rsp.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LWORK

#if defined (SYS_CRAY)
      REAL GD(*), WORK(LWORK)
      REAL DDOT
      REAL MONE, ZERO, ONE
#else
      DOUBLE PRECISION GD(*), WORK(LWORK)
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION MONE, ZERO, ONE
#endif
      PARAMETER (MONE = -1.0d0, ZERO = 0.0d0, ONE = 1.0d0)

      CHARACTER*(10) MODEL
      CHARACTER*(8)  LABEL1, LABEL2
      INTEGER ISYM, IADR1, IADR2, NTAMP, ICAU1, ICAU2, ISYM1, ISYM2
      INTEGER IOPT, IDLSTCL2, IDLSTX2, IDLSTMC1, IDLSTMC2, IDLSTR2

* external functions:
      INTEGER ICR2AMP
      INTEGER IR2TAMP
      INTEGER ICL2AMP
      INTEGER IL2ZETA
      INTEGER IETACL2
      INTEGER ICHI2

*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      ISYM1   = ISYCL2(IDLSTCL2,1)
      ISYM2   = ISYCL2(IDLSTCL2,2)
      ISYM    = MULD2H(ISYM1,ISYM2)

      LABEL1  = LBLCL2(IDLSTCL2,1)
      LABEL2  = LBLCL2(IDLSTCL2,2)

      ICAU1   = ICL2CAU(IDLSTCL2,1)
      ICAU2   = ICL2CAU(IDLSTCL2,2)

      IF (ICAU1.LT.0 .OR. ICAU2.LT.0) THEN
        WRITE (LUPRI,*) 'Error in CC_CL2GD> rhs vectors for Cauchy ',
     &          'orders < 0 not implemented!'
        CALL QUIT('Error in CC_CL2GD: illegal Cauchy order.')
      END IF

      IADR1 = 1
      IADR2 = IADR1 + NT1AM(ISYM)
      NTAMP = NT1AM(ISYM) + NT2AM(ISYM)
      IF (CCS) NTAMP = NT1AM(ISYM)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'IDLSTCL2:',IDLSTCL2
        WRITE (LUPRI,*) 'LABELs:',LABEL1, LABEL2
        WRITE (LUPRI,*) 'ICAUs:',ICAU1, ICAU2
        WRITE (LUPRI,*) 'NTAMP:',NTAMP
        WRITE (LUPRI,*) 'LWORK:',LWORK
      END IF

      IF (LWORK .LT. NTAMP) THEN
        WRITE (LUPRI,*) 'Insufficient memory in CC_CL2GD:'
        WRITE (LUPRI,*) 'available:',LWORK, '    needed:',NTAMP
        WRITE (LUPRI,*) 'for single part:',NT1AM(ISYM)
        WRITE (LUPRI,*) 'for doubles part:',NT2AM(ISYM)
        CALL QUIT('Insufficient memory in CC_CL2GD.')
      END IF

*---------------------------------------------------------------------*
* 1. contribution: CX2 vector
*---------------------------------------------------------------------*
      IOPT = 3
      IF ((ICAU1+ICAU2).EQ.0) THEN
        IDLSTX2 = ICHI2(LABEL1,.FALSE.,ZERO,ISYM1,
     &                  LABEL2,.FALSE.,ZERO,ISYM2)
        CALL CC_RDRSP('X2',IDLSTX2,ISYM,IOPT,MODEL,GD(IADR1),GD(IADR2))
      ELSE IF ((ICAU1+ICAU2).GT.0) THEN
        IDLSTX2 = IETACL2(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2,ISYM2)
        CALL CC_RDRSP('CX2',IDLSTX2,ISYM,IOPT,
     &                MODEL,GD(IADR1),GD(IADR2))
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_CL2GD> norm^2 of CX2-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_CL2GD> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYM),GD(IADR1),1,GD(IADR1),1)
        IF (.NOT.CCS) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_CL2GD> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYM),GD(IADR2),1,GD(IADR2),1)
      END IF

*---------------------------------------------------------------------*
* 2. contribution: F-matrix transformed CR2 vector
*---------------------------------------------------------------------*
      IOPT = 3
      IF ((ICAU1+ICAU2).EQ.0) THEN
        IDLSTR2 = IR2TAMP(LABEL1,.FALSE.,ZERO,ISYM1,
     &                    LABEL2,.FALSE.,ZERO,ISYM2)
        CALL CC_RDRSP('F2',IDLSTR2,ISYM,IOPT,
     &                MODEL,WORK(IADR1),WORK(IADR2))
      ELSE IF ((ICAU1+ICAU2).GT.0) THEN
        IDLSTR2 = ICR2AMP(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2,ISYM2)
        CALL CC_RDRSP('CF2',IDLSTR2,ISYM,IOPT,
     &                MODEL,WORK(IADR1),WORK(IADR2))
      END IF

      CALL DAXPY(NTAMP,ONE,WORK(IADR1),1,GD,1)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_CL2GD> norm^2 of CR2-vector:',
     &     DDOT(NTAMP,WORK,1,WORK,1)
        WRITE (LUPRI,*) 'DEBUG_CC_CL2GD> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYM),WORK(IADR1),1,WORK(IADR1),1)
        IF (.NOT.CCS) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_CL2GD> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYM),WORK(IADR2),1,WORK(IADR2),1)
      END IF

*---------------------------------------------------------------------*
* 3. contribution: CL2(ICAU1-1,ICAU2) vector
*---------------------------------------------------------------------*
      IF (ICAU1.GT.1 .OR. (ICAU1.EQ.1 .AND. ICAU2.GT.0) ) THEN

        IDLSTMC2 = ICL2AMP(LABEL1,ICAU1-1,ISYM1,LABEL2,ICAU2,ISYM2)

        IOPT = 3
        CALL CC_RDRSP('CL2',IDLSTMC2,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))

        CALL DAXPY(NTAMP,ONE,WORK(IADR1),1,GD,1)

      ELSE IF (ICAU1.EQ.1 .AND. ICAU2.EQ.0) THEN

        IDLSTMC2 = IL2ZETA(LABEL1,ZERO,ISYM1,LABEL2,ZERO,ISYM2)

        IOPT = 3
        CALL CC_RDRSP('L2',IDLSTMC2,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))

        CALL DAXPY(NTAMP,ONE,WORK(IADR1),1,GD,1)

      ELSE IF (ICAU1.EQ.0) THEN
        CONTINUE
      ELSE
        CALL QUIT('ERROR IN CC_CL2GD.')
      END IF
      
*---------------------------------------------------------------------*
* 4. contribution: CL2(ICAU1,ICAU2-1) vector
*---------------------------------------------------------------------*
      IF (ICAU2.GT.1 .OR. (ICAU2.EQ.1 .AND. ICAU1.GT.0) ) THEN
        IDLSTMC2 = ICL2AMP(LABEL1,ICAU1,ISYM1,LABEL2,ICAU2-1,ISYM2)

        IOPT = 3
        CALL CC_RDRSP('CL2',IDLSTMC2,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))

        CALL DAXPY(NTAMP,ONE,WORK(IADR1),1,GD,1)

      ELSE IF (ICAU2.EQ.1 .AND. ICAU1.EQ.0) THEN

        IDLSTMC2 = IL2ZETA(LABEL1,ZERO,ISYM1,LABEL2,ZERO,ISYM2)

        IOPT = 3
        CALL CC_RDRSP('L2',IDLSTMC2,ISYM,IOPT,MODEL,
     &                WORK(IADR1),WORK(IADR2))

        CALL DAXPY(NTAMP,ONE,WORK(IADR1),1,GD,1)

      ELSE IF (ICAU2.EQ.0) THEN
        CONTINUE
      ELSE
        CALL QUIT('ERROR IN CC_CL2GD.')
      END IF
      
*---------------------------------------------------------------------*
* write debug output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_CL2GD> norm^2 of GD-vector:',
     &     DDOT(NTAMP,GD,1,GD,1)
        WRITE (LUPRI,*) 'DEBUG_CC_CL2GD> norm^2 of singles part:',
     &     DDOT(NT1AM(ISYM),GD(IADR1),1,GD(IADR1),1)
        IF (.NOT.CCS) 
     &   WRITE (LUPRI,*) 'DEBUG_CC_CL2GD> norm^2 of doubles part:',
     &     DDOT(NT2AM(ISYM),GD(IADR2),1,GD(IADR2),1)
      END IF

      RETURN
      END
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_CL2GD                             *
*---------------------------------------------------------------------*
c /* deck CC_PTB  */
      SUBROUTINE CC_PTB(GD,LBL,ISY,FRQ,LRLX,WORK,LWORK)
c---------------------------------------------------------------------*
c
c    Purpose: Get solvent contribution to rhs. of t equations.
C             t-bar*P
c    Ove Christiansen, 13-5-1998.
c
c=====================================================================*
#include "implicit.h"
#include "maxorb.h"
#include "ccorb.h"
#include "iratdef.h"
#include "priunit.h"
#include "cclr.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccn2rsp.h"
#include "ccslvinf.h"
#include "qm3.h"
!#include "qmmm.h"

      PARAMETER( ONE = 1.0D00 , XMONE = -1.0D0)
      DIMENSION GD(*),WORK(LWORK)
      CHARACTER*(8)  LBL,LR*1
      CHARACTER MODEL*10,FILEX*10
      LOGICAL LEXIST,LRLX
C
      K1 = 1
      K2 = 1 + NT1AM(ISY)
      KWRK1 = K2 + NT2AM(ISY)
      LWRK1 = LWORK - KWRK1
      IF (LWRK1.LE.0) CALL QUIT('too little work in cc-ptb')
C
C----------------------------------
C     Get t-bar = L1 = zeta vector.
C----------------------------------
C
      IOPT = 33
      IF (IPRINT.GT.10) THEN
        WRITE(LUPRI,*) 'Solvent cont. to gd requires: '
        WRITE(LUPRI,*) 'Left vectors with : ',LBL,LRLX,FRQ,ISY
      ENDIF
      ILSTNR = IL1ZETA(LBL,LRLX,FRQ,ISY)
C
C     Note if CCSAV_LAM does not exist then
C     we have no contribution yet.
C
      IDXFIL = IDXSYM('L1',ISY,ILSTNR)
C
      WRITE(FILEX,'(A5,I1,1X,I3)') 'CCL1_',ISY,IDXFIL
      DO I = 1, 10
        IF ( FILEX(I:I) .EQ. ' ' ) FILEX(I:I) = '_'
      END DO
C
      INQUIRE(FILE=FILEX,EXIST=LEXIST)
C
      IF (IPRINT.GT.15) THEN
         WRITE(LUPRI,*) ' In CC_PTB: calculating t-barP transf.'
         WRITE(LUPRI,*) 'Isy,LBL,FRQ = ',Isy,LBL,FRQ
         WRITE(LUPRI,*) 'Ilstnr = ',ilstnr
         WRITE(LUPRI,*) 'LEXIST,FILEX= ',LEXIST,FILEX
      ENDIF
C
      IF (.NOT.LEXIST) THEN
         IF (IPRINT.GT.1)
     *    WRITE(LUPRI,*) FILEX,' does not exits yet - no t-bar P cont'
         RETURN
      ENDIF
C
      CALL CC_RDRSP('L1',ILSTNR,ISY,IOPT,MODEL,WORK(K1),WORK(K2))
C
C---------------------------
C     Make P-transformation.
C---------------------------
C
       LR = 'P'
       IF (.NOT. CCMM) CALL CCSL_LTRB(GD(1),GD(1+NT1AM(ISY)),WORK(K1),
     *                 WORK(K2),ISY,LR,WORK(KWRK1),LWRK1)
       IF (CCMM) THEN
          IF (.NOT. NYQMMM) THEN 
             CALL CCMM_LTRB(GD(1),GD(1+NT1AM(ISY)),WORK(K1),
     *                   WORK(K2),ISY,LR,WORK(KWRK1),LWRK1)
          ELSE IF (NYQMMM) THEN
             CALL CCMM_TRANSFORMER(GD(1),GD(1+NT1AM(ISY)),WORK(K1),
     *                     WORK(K2),MODEL,ISY,LR,WORK(KWRK1),LWRK1)
C             CALL CCMM_TGTBAR(GD(1),GD(1+NT1AM(ISY)),WORK(K1),WORK(K2),
C     *                       MODEL,ISY,LR,WORK(KWRK1),LWRK1)
          END IF
       END IF
C
C
      END

