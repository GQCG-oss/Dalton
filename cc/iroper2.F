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
*=====================================================================*
C  /* Deck iroper2 */
*=====================================================================*
      INTEGER FUNCTION IROPER2(NEWLBLA,NEWLBLB,LABSOP,ISIGN,ISYSOP)
*---------------------------------------------------------------------*
*
* maintain the list of second-order operators for perturbations
* involving field-dependent basis sets
*
* if operator is on the list, return list index and set ISYSOP, ISIGN
* and, depending on input, LABSOP or NEWLBLA, NEWLBLB:
*   1) if NEWLBLA/NEWLBLB='?' search for LABSOP and set NEWLBLA/NEWLBLB
*   2) else search for NEWLBLA/NEWLBLB and set LABSOP
* if operator is NOT on the list:
*       LOPR2OPN = .true.  --> extend list, return index
*       LOPR2OPN = .false. --> return -1
*
* Christof Haettig, March 99
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccroper.h>
#include <ccropr2.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LOPNSAVE
      CHARACTER*8 NEWLBLA, NEWLBLB, LABSOP
      INTEGER I, IOP2, ISIGN, ISYSOP, IOPA, IOPB, ISYMA, ISYMB
* external functions:
      INTEGER IROPER

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'IROPER2>',NEWLBLA,NEWLBLB
        CALL FLSHFO(LUPRI)
      END IF
     
      IF (NEWLBLA(1:1).EQ.'?' .OR. NEWLBLB(1:1).EQ.'?') THEN
         DO I = 1, NRSO2LBL
         IF ( LABSOP.EQ.LBLOP2(I,3) ) THEN
           IROPER2 = I
           IF (NEWLBLA(1:1).EQ.'?' .AND. NEWLBLB(1:1).EQ.'?') THEN
             NEWLBLA = LBLOP2(I,1)
             NEWLBLB = LBLOP2(I,2)
           ELSEIF (NEWLBLA(1:1).EQ.'?'.AND.NEWLBLB.EQ.LBLOP2(I,2)) THEN
             NEWLBLA = LBLOP2(I,1)
           ELSE IF (NEWLBLA(1:1).EQ.'?'.AND.NEWLBLB.EQ.LBLOP2(I,1)) THEN
             NEWLBLA = LBLOP2(I,2)
           ELSEIF (NEWLBLB(1:1).EQ.'?'.AND.NEWLBLA.EQ.LBLOP2(I,2)) THEN
             NEWLBLB = LBLOP2(I,1)
           ELSEIF (NEWLBLB(1:1).EQ.'?'.AND.NEWLBLA.EQ.LBLOP2(I,1)) THEN
             NEWLBLB = LBLOP2(I,2)
           ELSE
             CONTINUE
           END IF
           ISYSOP  = ISYOP2(I)
           ISIGN   = ISGNOP2(I)
        
           IF (LOCDBG) THEN
             WRITE (LUPRI,*) 
     *             'IROPER2> ISYSOP,ISIGN,LABSOP:',ISYSOP,ISIGN,LABSOP
             WRITE (LUPRI,*) 
     *             'IROPER2> NEWLBLA,NEWLBLB:',NEWLBLA,NEWLBLB
           END IF

           RETURN
         END IF
         END DO
         IF (.NOT. LQUIET) THEN
            WRITE(LUPRI,'(/5A)') ' WARNING: SECOND-ORDER OPERATOR "',
     *           LABSOP,'" NOT AVAILABLE.'
            WRITE (LUPRI,*) 
     *           'IROPER2> ISYSOP,ISIGN,LABSOP:',ISYSOP,ISIGN,LABSOP
            WRITE (LUPRI,*) 
     *           'IROPER2> NEWLBLA,NEWLBLB:',NEWLBLA,NEWLBLB
         END IF
         IROPER2 = -1
         RETURN
      END IF


      DO I = 1,NRSO2LBL
         IF ((NEWLBLA.EQ.LBLOP2(I,1) .AND. NEWLBLB.EQ.LBLOP2(I,2)) .OR.
     &       (NEWLBLA.EQ.LBLOP2(I,2) .AND. NEWLBLB.EQ.LBLOP2(I,1)))THEN
           IROPER2 = I
           LABSOP  = LBLOP2(I,3)
           ISYSOP  = ISYOP2(I)
           ISIGN   = ISGNOP2(I)

           IF (LOCDBG) THEN
             WRITE (LUPRI,*) 
     *             'IROPER2> ISYSOP,ISIGN,LABSOP:',ISYSOP,ISIGN,LABSOP
             WRITE (LUPRI,*) 
     *             'IROPER2> NEWLBLA,NEWLBLB:',NEWLBLA,NEWLBLB
           END IF

           RETURN
         END IF
      END DO  


      IF (LOPR2OPN) THEN
        NRSO2LBL = NRSO2LBL + 1

        IF (NRSO2LBL.GT.MAXOP2LBL) THEN
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *    ' NUMBER OF SPECIFIED OPERATORS EXCEED THE MAXIMUM ALLOWED',
     *    ' MAXOP2LBL =',MAXOP2LBL,' NRSO2LBL= ',NRSO2LBL
         CALL QUIT(' IROPER2: TOO MANY OPERATORS SPECIFIED')
        END IF

        IOPA = IROPER(NEWLBLA,ISYMA)
        IOPB = IROPER(NEWLBLB,ISYMB)

        LBLOP2(NRSO2LBL,1) = NEWLBLA
        LBLOP2(NRSO2LBL,2) = NEWLBLB
        LBLOP2(NRSO2LBL,3) = LABSOP

        ISYOP2(NRSO2LBL)   = ISYSOP
        ISGNOP2(NRSO2LBL)  = ISIGN
        IATOPR2(NRSO2LBL)  = IATOPR(IOPA) + IATOPR(IOPB)
        LPDBSOP2(NRSO2LBL) = (LPDBSOP(IOPA) .AND. LPDBSOP(IOPB))

        IROPER2 = NRSO2LBL

        ! put the operator also on the usual operator list:
        LOPNSAVE = LOPROPN
        LOPROPN  = .TRUE.

        IOP2 = IROPER(LABSOP,ISYSOP)
        ISYMAT(IOP2)  = 0
        IATOPR(IOP2)  = IATOPR2(IROPER2)
        LPDBSOP(IOP2) = LPDBSOP2(IROPER2)

        LOPROPN = LOPNSAVE

      ELSE
        IF (.NOT. LQUIET) WRITE(LUPRI,'(/5A)')
     &   ' WARNING: SECOND-ORDER OPERATOR FOR PERTURBATION PAIR "',
     &    NEWLBLA,'", "',NEWLBLB,'" NOT AVAILABLE.'
        IROPER2 = -1
      END IF

      IF (LOCDBG)  THEN
        WRITE (LUPRI,*) 'IROPER2> ',IROPER2,LBLOP2(IROPER2,3),
     &                         ISGNOP2(IROPER2),ISYOP2(IROPER2)
      END IF

      RETURN
      END
