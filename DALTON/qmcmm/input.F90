!...
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

module qmcmm_input

   implicit none

   public qmnpinp

   private

contains

   subroutine qmnpinp(word)
!
! Purpose:
!   reads "*QMNPMM" input group and sets up various flags for
!   QM/NP/MM type embedding
!
! Input:
!   WORD - DALTON input group keyword
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

#include "priunit.h"
#include "gnrinf.h"
#include "qmnpmm.h"

      integer, parameter :: NTABLE = 13
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      integer :: ichang, i
!
      DATA TABLE /'.QMNP  ', '.QMNPMM', '.QMMM  ', '.NPPOIN',           &
             '.xxxxxx', '.PRINT ', '.MQITER', '.xxxxxx',                &
             '.NONPCA', '.MMPOLA', '.xxxxxx', '.XXXXXX',                &
             '.DAMPED'/
!
!
      CALL SET_QMNPMM()
!
      NEWDEF = (WORD .EQ. '*QMNPMM')
      ICHANG = 0
!
      IF (NEWDEF) THEN
        WORD1 = WORD
!
  100   CONTINUE
!
        READ (LUCMD, '(A7)') WORD
        CALL UPCASE(WORD)
        PROMPT = WORD(1:1)
!
        IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
          GO TO 100
        ELSE IF (PROMPT .EQ. '.') THEN
          ICHANG = ICHANG + 1
          DO 200 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
              GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13), I
            END IF
  200     CONTINUE
          IF (WORD .EQ. '.OPTION') THEN
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            GO TO 100
          END IF
!
          WRITE (LUPRI,'(/3A/)') ' Keyword "',WORD,                     &
                            '" not recognized for *QMNPMM'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal keyword for *QMNPMM')

!         ".QMNP "
    1     CONTINUE
            QMNPMM  = .TRUE.
            DONPSUB = .TRUE.
            DOMMSUB = .FALSE.
            DONPCAP = .TRUE.
            DONPPOL = .TRUE.
          GO TO 100

!         ".QMNPMM"
    2     CONTINUE
            QMNPMM  = .TRUE.
            DONPSUB = .TRUE.
            DOMMSUB = .TRUE.
            DONPCAP = .TRUE.
            DONPPOL = .TRUE.
          GO TO 100

!         ".QMMM  "
    3     CONTINUE
            QMNPMM  = .TRUE.
            DONPSUB = .FALSE.
            DOMMSUB = .TRUE.
          GO TO 100

!         ".NPPOIN"
    4     CONTINUE
            NPMQGAU = .FALSE.
          GO TO 100
!
    5     CONTINUE
          GO TO 100

!         ".PRINT "
    6     CONTINUE
            READ(LUCMD,*) IPRTLVL
          GO TO 100

!         ".MQITER"
    7     CONTINUE
            MQITER = .TRUE.
          GO TO 100
!
    8     CONTINUE
          GO TO 100

!         ".NONPCA"
    9     CONTINUE
            DONPCAP = .FALSE.
          GO TO 100

!         ".MMPOLA"
   10     CONTINUE
            DOMMPOL = .TRUE.
          GO TO 100
!
   11     CONTINUE
          GO TO 100

!         ".XXXXXX"
   12     CONTINUE
          GO TO 100

!         ".DAMPED"
 13       CONTINUE
            NOVDAMP = .FALSE.
          GO TO 100

        ELSE IF (PROMPT .EQ. '*') THEN
          GO TO 300
        ELSE
          WRITE (LUPRI,'(/3A/)') ' Prompt "',WORD,                      &
                              '" not recognized for *QMNPMM'
          CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
          CALL QUIT('Illegal prompt for *QMNPMM')
        END IF
      END IF
!
  300 CONTINUE
!
   end subroutine

   subroutine set_qmnpmm()
!
! Purpose:
!   sets initial values of all parameters related to
!   QM/NP/MM type embedding
!
! Last updated: 22/03/2013 by Z. Rinkevicius.
!

#include "qmnpmm.h"
!
      DONPSUB = .FALSE.
      DOMMSUB = .FALSE.
      DONPCAP = .FALSE.
      DOMMCAP = .FALSE.
      DONPPOL = .FALSE.
      DOMMPOL = .FALSE.
      NPMQGAU = .TRUE.
      MQITER  = .FALSE.
      NOVDAMP = .TRUE.
      TNPBLK  = 0
      TMMBLK  = 0
      IPRTLVL = 1
      TNPATM  = 0
      TMMATM  = 0
      TPOLATM = 0
      TNPFF   = 0
      TMMFF   = 0
!
      CALL DZERO(NPCORD,3*MXNPATM)
      CALL DZERO(mm_cord,3*MXMMATM)
      CALL DZERO(NPCHRG,MAXBLK)
      CALL DZERO(MMCHRG,MAXBLK)
      CALL IZERO(NPFTYP,MXNPATM)
      CALL IZERO(MMFTYP,MXMMATM)
      CALL IZERO(NPATOM,MAXBLK)
      CALL IZERO(MMATOM,MAXBLK)
      CALL DZERO(NPFPOL,MXNPFF)
      CALL DZERO(NPFCAP,MXNPFF)
      CALL DZERO(NPFOMG1,MXNPFF)
      CALL DZERO(NPFGAM1,MXNPFF)
      CALL DZERO(NPFOMG2,MXNPFF)
      CALL DZERO(NPFGAM2,MXNPFF)
      CALL DZERO(NPFFAC,MXNPFF)
      CALL DZERO(MMFM0,MXMMFF)
      CALL DZERO(MMFPOL,MXMMFF)
      CALL DZERO(MMMOL,MXMMATM)
      CALL IZERO(MMSKIP,MXMMATM)
!
      ENSOLQNP = 0.0D0
      EESOLQNP = 0.0D0
      ENSOLMNP = 0.0D0
      EESOLMNP = 0.0D0
!
      ENSOLQMM = 0.0D0
      EESOLQMM = 0.0D0
      ENSOLMMM = 0.0D0
      EESOLMMM = 0.0D0
!
   end subroutine

end module
