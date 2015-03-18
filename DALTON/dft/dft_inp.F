!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
      SUBROUTINE dftINP(WORD,INPERR,ALLOPT)
#include "implicit.h"
#include "gnrinf.h"
#include "dftinf.h"
#include "dftcom.h"
#include "priunit.h"
#include "mxcent.h"
#include "dftd.h"
      PARAMETER ( NTABLE = 26 )
#include "dftacb.h"
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7, LINE*80
      LOGICAL ALLOPT
      ! Checks if the parameters for DFT-D are defined
      EXTERNAL DFT_D2_CHECK, DFT_D3_CHECK, DFT_D3BJ_CHECK
      INTEGER DFT_D2_CHECK, DFT_D3_CHECK, DFT_D3BJ_CHECK
      DATA TABLE /'.PRINT ','.DFTELS','.DFTTHR','.RADINT','.ANGINT',
     &            '.ANGMIN','.NOPRUN','.GRID T','.DFTAC ','.CARPAR',
     &            '.COARSE','.NORMAL','.FINE  ','.ULTRAF','.GRID1 ',
     &            '.GRID2 ','.GRID3 ','.GRID4 ','.GRID5 ','.DFTD2 ',
     &            '.D2PAR ','.DFTD3 ','.DFD3BJ','.3BODY ','.D3PAR ',
     &            '.DFTVXC'/
C
      WORD1 = WORD
C
 100  CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
     &                      17,18,19,20,21,22,23,24,25,26), I
                  END IF
 200           CONTINUE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *              '" not recognized in DFTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in DFTINP.')
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
              WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in DFTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 CALL QUIT('Illegal keyword in DFTINP.')
            END IF
c     keyword handling section
c           >PRINT < keyword
 1    CONTINUE
         READ(LUCMD,*) IPRDFT
      GO TO 100  
c     >.ELECTRONS< check keyword
 2    CONTINUE
         DFTHRS = .TRUE.
         READ(LUCMD,*) DFTELS
      GO TO 100  
c     >.DFTTHR< check keyword
 3    CONTINUE
         DFTHRS = .TRUE.
         READ (LUCMD,*) DFTHR0, DFTHRL, DFTHRI
      GO TO 100  
c     >.RADINT< check keyword
 4    CONTINUE
         READ (LUCMD,*) RADINT
      GO TO 100  
c     >.ANGINT< check keyword
 5    CONTINUE
         READ(LUCMD,*) ANGINT
      GO TO 100  
c     >.ANGMIN< check keyword
 6    CONTINUE
         READ(LUCMD,*) ANGMIN
      GO TO 100  
c     >.NOPRUN< check keyword
 7    CONTINUE
         NOPRUN = .TRUE.
      GO TO 100  
c     >.GRID T< grid type keyword
 8    CONTINUE
         READ(LUCMD,'(A80)') LINE
         CALL DFTGRIDINPUT(LINE)
      GO TO 100  
c     >.DFTAC < asymptotic correction keyword
 9    CONTINUE
      CALL ACINP
      GO TO 100
c     >.CARPAR< cartesian grid parameters.
 10   CONTINUE
         READ(LUCMD,'(A80)') LINE
         CALL DFTCARTESIANINPUT(LINE)
      GO TO 100 
c     >.COARSE< grid parameters.
 11   CONTINUE
         RADINT = 1D-11
         ANGINT = 35
      GO TO 100 
c     >.NORMAL< grid parameters.
 12   CONTINUE
         RADINT = 1D-13
         ANGINT = 35
      GO TO 100 
c     >.FINE  < grid parameters.
 13   CONTINUE
         RADINT = 1D-13
         ANGINT = 42
      GO TO 100 
c     >.ULTRAF< grid parameters.
 14   CONTINUE
         RADINT = 2D-15
         ANGINT = 64
      GO TO 100 
c     >.GRID1 < grid parameters.
 15   CONTINUE
         CALL DFTGRIDINPUT("GC2")
         RADINT = 1D-3
         ANGINT = 17
      GO TO 100 
c     >.GRID2 < grid parameters.
 16   CONTINUE
         CALL DFTGRIDINPUT("GC2")
         RADINT = 1D-5
         ANGINT = 23
      GO TO 100 
c     >.GRID3 < grid parameters.
 17   CONTINUE
         CALL DFTGRIDINPUT("GC2")
         RADINT = 1D-6
         ANGINT = 29
      GO TO 100 
c     >.GRID4 < grid parameters.
 18   CONTINUE
         CALL DFTGRIDINPUT("GC2")
         RADINT = 1D-7
         ANGINT = 35
      GO TO 100 
c     >.GRID5 < grid parameters.
 19   CONTINUE
         CALL DFTGRIDINPUT("GC2")
         RADINT = 1D-9
         ANGINT = 47
      GO TO 100 
c    >.DFTD2 < turns empirical dispersion DFTD2 correction on
 20   DODFTD = .TRUE.
      DO_DFTD2=.TRUE.
      GO TO 100
c    >.D2PAR < Read in custom DFT-D2 parameters
 21   L_INP_D2PAR = .TRUE.
      READ(LUCMD,*) D2_s6_inp, D2_alp_inp, D2_rs6_inp
      GO TO 100
c    >.DFTD3 < turns empirical dispersion DFT-D3 correction on
 22   DODFTD = .TRUE.
      DO_DFTD3 = .TRUE.
      GO TO 100
c    >.DFD3BJ< turns empirical dispersion DFT-D3 correction with
c              Becke-Johnson damping
 23   DODFTD = .TRUE.
      DO_DFTD3 = .TRUE.
      DO_BJDAMP = .TRUE.
      GO TO 100
c   >.3BODY < turns empirical dispersion 3-body terms 
 24   DO_3BODY = .TRUE.
      GO TO 100
c   >.D3PAR < Read in custom DFT-D3 parameters
 25   L_INP_D3PAR = .TRUE.
      READ(LUCMD,*) D3_s6_inp, D3_alp_inp, D3_rs6_inp,
     &              D3_rs18_inp, D3_s18_inp
      GO TO 100
CAMT >.DFTVXC< Consruct V_xc Explicitly and Use in KS Matrix for SCF
 26   CONTINUE
      LDFTVXC = .TRUE.
      GO TO 100

 300  CONTINUE
      IF (THR_REDFAC .GT. 0.0D0) THEN
C arnfinn: ichang is not defined or used...
C         ICHANG = ICHANG + 1
         WRITE (LUPRI,'(3A,1P,D10.2)') '@ INFO ',WORD1,
     &   ' thresholds multiplied with general factor',THR_REDFAC
         RADINT = RADINT*THR_REDFAC
         DFTHR0 = DFTHR0*THR_REDFAC
         DFTHRL = DFTHRL*THR_REDFAC
         DFTHRI = DFTHRI*THR_REDFAC
         DFTELS = DFTELS*THR_REDFAC
         DFTHRS = .TRUE.
      END IF

      IF (DO_DFTD2.AND.(.NOT.L_INP_D2PAR)) THEN
C        empirical disperision DFT-D2 correction only defined for a
C        number of functionals 
C        check this here and quit if not defined  
         ICHK=DFT_D2_CHECK()
         IF (ICHK.NE.1)
     &       CALL QUIT('ERROR: DFT-D2 Func. Dep. Paras. Unknown')
      END IF

      IF (DO_DFTD3.AND.(.NOT.L_INP_D3PAR)) THEN
C        empirical disperision DFT-D3 correction only defined for a
C        number of functionals 
C        check this here and quit if not defined  
         ICHK=DFT_D3_CHECK()
         IF (ICHK.NE.1)
     &       CALL QUIT('ERROR: DFT-D3 Func. Dep. Paras. Unknown')
      END IF

      IF (DO_3BODY.AND.(.NOT.DO_DFTD3)) THEN
         CALL QUIT('ERROR: 3-BODY CORR. APPLY ONLY FOR DFT-D3')
      ENDIF


      CALL HEADER('Settings for DFT calculation:',0)
c     code goes here
         IF (IPRDFT .NE. IPRUSR) THEN
            WRITE (LUPRI,'(4X,A,T160,I5)')
     &        ' Redefined print level for DFT :',IPRDFT
         END IF
         IF (DFTHRS) THEN
            WRITE (LUPRI,'(4X,A,T60,1P,3D12.2)')
     &        ' Redefined thresholds:', DFTHR0, DFTHRL, DFTHRI 
         ELSE
            WRITE (LUPRI,'(4X,A,T60,1P,3D12.2)')
     &        ' Default thresholds:  ', DFTHR0, DFTHRL, DFTHRI 
         END IF
         IF (DFTHRS) THEN
            WRITE (LUPRI,'(4X,A,T60,1P,3D12.2)')
     &        ' Redefined threshold for number of electrons: ', DFTELS 
         ELSE
            WRITE (LUPRI,'(4X,A,T60,1P,3D12.2)')
     &        ' Default threshold for number of electrons:   ', DFTELS 
         END IF
         WRITE (LUPRI,'(4X,A,T74,1P,D10.2,I4)')
     $      ' DFT radial quadrature accuracy and ang. expansion order:',
     &      RADINT,ANGINT

      if (angint .gt. 64) then
!radovan:
!       angint > 64 not implemented
!       higher angint defaults to the poorest angular grid!
!       so better quit here in this case
        WRITE(LUPRI,*) 'FATAL ERROR: angint > 64 not implemented.'
        call quit('FATAL ERROR: angint > 64 not implemented')
      end if

      END


      SUBROUTINE ACINP
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "dftcom.h"
#include "dftacb.h"
#include "gnrinf.h"
      CHARACTER*8 RTYPE, CTYPE
      DFTASC = .TRUE.
      DFTPOT = .TRUE.
      LDFTVXC = .TRUE.
C Initialize
      LGRAC = .FALSE.
      LLIN = .FALSE.
      LTAN = .FALSE.
      DOLB94 = .FALSE.
      DOMPOLE = .FALSE.

      READ (LUCMD,*) RTYPE
      CALL UPCASE(RTYPE)
      READ (LUCMD,*) CTYPE
      CALL UPCASE(CTYPE)

      IF (CTYPE.EQ.'GRAC    ') THEN
        LGRAC = .TRUE.
      ELSE IF (CTYPE.EQ.'LINEAR  ') THEN
        LLIN = .TRUE.
      ELSE IF (CTYPE.EQ.'TANH    ') THEN
        LTAN = .TRUE.
      ELSE
        CALL QUIT('ERROR: Unrecognised AC Interp')
      ENDIF

      IF (RTYPE.EQ.'LB94    ') THEN
        DOLB94 = .TRUE.
      ELSE IF (RTYPE.EQ.'MULTPOLE') THEN
        DOMPOLE = .TRUE.
      ELSE
        CALL QUIT('ERROR: Unrecognised AC potential')
      ENDIF

      READ (LUCMD,*) DFTIPTA, DFTIPTB, DFTBR1, DFTBR2

      CALL AROUND('Asymptotically Corrected XC Potential')

      CALL HEADER('Interpolation Scheme',0)
      IF (LGRAC) THEN
      WRITE(LUPRI,'(4X,A35)')' Gradient Regulated Connection    '
      ENDIF
      IF (LLIN) THEN
      WRITE(LUPRI,'(4X,A35)')' Linear Interpolation             '
      ENDIF
      IF (LTAN) THEN
      WRITE(LUPRI,'(4X,A35)')' tanh Interpolation               '
      ENDIF
      WRITE(LUPRI,'(4X,A35)')'                                  '
      CALL HEADER('Asymptotic Potential',0)
      IF (DOLB94) THEN
      WRITE(LUPRI,'(4X,A35)')' LB94 Exchange-Correlation Form   '
      ENDIF
      IF (DOMPOLE) THEN
      WRITE(LUPRI,'(4X,A35)')' Asymp. Pot. From Multipole Exp.  '
      ENDIF
      CALL HEADER('Input Ionization Potentials',0)
      WRITE(LUPRI,5)' Alpha IP         :',  DFTIPTA
      WRITE(LUPRI,5)' Beta IP          :',  DFTIPTA
5     FORMAT (5x,a18,f6.3)
      IF (.NOT. LGRAC) THEN
      CALL HEADER('Interpolation Region Definition',0)
      WRITE(LUPRI,5)' Bragg Radius 1   :',DFTBR1
      WRITE(LUPRI,5)' Bragg Radius 2   :',DFTBR2
      ENDIF
      IF (LGRAC) THEN
      CALL HEADER('Interpolation Region Definition',0)
      WRITE(LUPRI,5)' Alpha Parameter  :',DFTBR1
      WRITE(LUPRI,5)' Beta Parameter   :',DFTBR2
      ENDIF
      IF (LGRAC .AND.DOMPOLE) THEN
        CALL QUIT('ERROR: Multipole AC not permitted with GRAC')
      ENDIF
      IF (LGRAC .AND. PARCAL) CALL PARQUIT("GRAC")
      RETURN
      END


