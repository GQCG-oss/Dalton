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
c/* Deck cc_rdrsp */
      SUBROUTINE CC_RDRSP(LIST,IDXLST,ISYM,IOPT,MODFIL,VEC1,VEC2)
C---------------------------------------------------------------------*
C
C   Purpose:  Read a response vector from file:
C             - the vector is addressed via a list type LIST and 
C               its index, IDXLST, with in that list. 
C               IDXLST is calculated from IDXLST and symmetry info.
C             - it is checked, that the vector has the symmetry ISYM. 
C             - according to the value of the read option IOPT and
C               its compatibility with the CC model MODFIL for which
C               the vector is stored on file, the single and double
C               excitation parts are returned in VEC1 and VEC2.
C
C      read option IOPT, used as bit wise flags:
C
C                   1  :  read the single excitation part
C                   2  :  read the double excitation part
C                   3  :  read both single and double excitation part
C                   4  :  read the CPHF part (not available here...)
C                   8  :  read the single excitation part 
C                         of the 'effective' vector for triples models
C                  16  :  read the double excitation part 
C                         of the 'effective' vector for triples models
C                  24  :  read both single and double excitation part
C                         of the 'effective' vector for triples models
C                  32  :  read the R12 double excitation part
C
C                  33  :  read both single and double excitation part
C                         SUCCES: IOPT is returned as 1 or 3.
C                         NO SUCCES: opening file fails, 
C                         it will not stop but return IOPT 33.
C                         (out dated, use -IOPT instead)
C
C                -IOPT :  as IOPT, but does not stop on failure,
C                         instead it returns IOPT=33, 
C                         on succes IOPT is returned with positive sign
C
C                  if IOPT=3 is obtained, but the double excitation
C                  part is not available (f.x. because it is a CCS
C                  vector) IOPT is reset to 1
C
C      NOTE: IOPT is input and output!!! do not pass a parameter
C            or integer constant ...
C       
C
C      List types LIST:
C
C          'R0' : right zeroth-order amplitude vector (IDXLST ignored)
C          'R1' : right first-order response vector 
C          'R2' : right second-order response vector  
C          'R3' : right third-order response vector  
C          'R4' : right fourth-order response vector  
C          'RC' : right first-order Cauchy vectors 
C          'CR2': right second-order Cauchy vectors 
C          'RE' : right eigenvectors
C          'ER1': first-order response of right eigenvetors
C          'ER2': second-order response of right eigenvetors
C
C          'L0' : left zeroth-order lamdba vector (IDXLST ignored)
C          'L1' : left first-order response vector
C          'L2' : left second-order response vector
C          'L3' : left third-order response vector
C          'L4' : left fourth-order response vector
C          'LC' : left first-order Cauchy vectors 
C          'CL2': left second-order Cauchy vectors 
C          'LE' : left eigenvectors
C          'EL1': first-order response of left eigenvectors
C          'EL2': second-order response of left eigenvectors
C          'PL1': projected left first-order response vector 
C
C          'F1' : F-transformed first-order right vectors
C          'F2' : F-transformed second-order right vectors
C          'F3' : F-transformed third-order right vectors
C          'F4' : F-transformed fourth-order right vectors
C          'FC' : F-transformed first-order right Cauchy vectors
C          'CF2': F-transformed second-order right Cauchy vectors
C          'FR' : F-transformed right eigenvector
C
C          'XF1': F-transformed right first-order response vectors,
C                 F12*R2 contribution for Cholesky CC2.
C
C
C          'O1 ': rhs for first-order amplitude equations
C          'O2' : rhs for second-order amplitude equations
C          'O3' : rhs for third-order amplitude equations
C          'O4' : rhs for fourth-order amplitude equations
C          'CO2': rhs for second-order cauchy equations
C          'EO1': rhs for right eigenvector first-order response
C          'EO2': rhs for right eigenvector second-order response
C
C          'eO1': effective rhs for first-order amplitude equations for
C                 Cholesky CC2
C
C
C          'X1' : first-order eta vectors
C          'X2': second-order eta vectors
C          'X3': third-order eta vectors
C          'X4': fourth-order eta vectors
C          'CX2': second-order Cauchy eta vectors
C          'EX1': chi vector for left eigenvector first-order response
C          'EX2': chi vector for left eigenvector second-order response
C
C          'M1': zeroth-order Lagrangian multipliers for ground state
C                -- excited state transition moments
C
C          'BE': rhs for zeroth-order excited state Lagrangian 
C                multipliers. (LE*B*RE)
C          'E0': zeroth-order excited state Lagrangian multipliers.
C                (the same as the N2 vectors, but diagonal case only)
C
C          'BR': rhs for N2 vectors
C
C          'N2': zeroth-order Lagrangian multipliers for 
C                excited state -- excited state transition moments
C                and excited state properties
C
C          'D0': dummy vector
C
C
C  Christof Haettig, October 1996, restructered in May 1997
C  PL1 vectors introduced in March 2000 by Sonia Coriani
C  CC3 effective vectors introduced 2002 by Christof Haettig
C  Cholesky extensions by Thomas Bondo Pedersen, March/April 2003.
C  R12 double excitation part introduced Jun 2003 by Christof Haettig
C
C---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <dummy.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccexci.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
C
      CHARACTER SUBRNAME*(8)
      PARAMETER (SUBRNAME = 'CC_RDRSP')
      INTEGER LUSAVE
C
      INTEGER IDXLST, ISYM, IOPT, IDXVEC
      CHARACTER LIST*(*),LISTI*(4), LIST3*(3)
C
      CHARACTER FILEX*(10), MODFIL*(10)
      LOGICAL LEXIST
      INTEGER IOS, IMUL, NVEC, IOPTSV, IOPTIN

#if defined (SYS_CRAY)
      REAL VEC1(*), VEC2(*)
      REAL ZERO
#else
      DOUBLE PRECISION VEC1(*), VEC2(*)
      DOUBLE PRECISION ZERO
#endif
      PARAMETER (ZERO = 0.0d0)

* external functions:
      INTEGER IDXSYM
      INTEGER ILSTSYM
*
      IF (LEN(LIST) .EQ. 2) THEN
         LIST3 = LIST//' '
      ELSE
         LIST3 = LIST(1:3)
      END IF
 
*---------------------------------------------------------------------*
*     get file name and modified LIST label and IDXVEC:
*---------------------------------------------------------------------*
      CALL CC_RWPRE(LIST,IDXLST,ISYM,LISTI,IDXVEC,FILEX)

*---------------------------------------------------------------------*
*     Test if abort under fail 
*     If Cholesky CC2, then IOPT=33 means singles if available, please!
*     [Note: IOPT=3 is actually OK, as IOPT is reset during the check
*      of the header as IOPT = IAND(IOPT,5)]
*---------------------------------------------------------------------*
      IF (IOPT .EQ. 33) THEN
         IOPTSV = 33
Chol     IOPT   = 3
         IF (CHOINT .AND. CC2) THEN
            IOPT = 1
         ELSE
            IOPT   = 3
         ENDIF
      ELSE IF (IOPT.LT.0) THEN
         IOPTSV = 33
         IOPT   = -IOPT
      ELSE
         IOPTSV = 0
      ENDIF
      IOPTIN = IOPT

*---------------------------------------------------------------------*
* return a fake zero vector for CCS R0/L0 vector:
*---------------------------------------------------------------------*
      IF ( (LIST3(1:2).EQ.'R0'.OR.LIST3(1:2).EQ.'L0') .AND. CCS ) THEN
         IOPT = IAND(IOPT,1)
         IF (IOPT.NE.0) CALL DZERO(VEC1,NT1AMX)
         RETURN
      END IF

*---------------------------------------------------------------------*
* open & rewind file:
*---------------------------------------------------------------------*
      LUSAVE = -1

      INQUIRE(FILE=FILEX,EXIST=LEXIST,IOSTAT=IOS,ERR=990)
      IF (.NOT. LEXIST) GOTO 991

      CALL GPOPEN(LUSAVE,FILEX,'OLD','SEQUENTIAL','UNFORMATTED',IDUMMY,
     &            .FALSE.)

      REWIND(LUSAVE,IOSTAT=IOS,ERR=992)

*---------------------------------------------------------------------*
* read and check file header:
*---------------------------------------------------------------------*
      CALL CC_RCRSP(LISTI,IDXVEC,LUSAVE,FILEX,MODFIL,IOPT,IOPTSV)
      IF (IOPT.LT.0) GOTO 993

      IF (IOPT.EQ.0) THEN
        IF (IAND(IOPTIN,3).GT.0) THEN
         WRITE (LUPRI,'(2A)') 
     &    ' no appropriate singles/doubles vector found on file ',FILEX
        ELSE IF (IAND(IOPTIN,24).GT.0) THEN
         WRITE (LUPRI,'(2A)') 
     &    ' no appropriate eff. CC3 rhs vector found on file ',FILEX
        ELSE IF (IAND(IOPTIN,32).GT.0) THEN
         WRITE (LUPRI,'(2A)') 
     &    ' no appropriate R12 doubles vector found on file ',FILEX
        ELSE
         WRITE (LUPRI,'(2A)') 
     &    ' no appropriate vector found on file ',FILEX
        END IF
        GOTO 996
      END IF

*---------------------------------------------------------------------*
* find the multiplicity of the vector:
*---------------------------------------------------------------------*

      IF (LIST3(1:2).EQ.'RE' .OR. LIST3(1:2).EQ.'LE') THEN
      ! for excited states we get it from the common block:
         IMUL = IMULTE(IDXLST)
      ELSE
      ! default is singlet
         IMUL = 1
      END IF
                                                
*---------------------------------------------------------------------*
* read the vectors & and close file:
*---------------------------------------------------------------------*
      CALL CC_RVRSP(ISYM, IMUL, IOPT, DUMMY, LUSAVE, VEC1, VEC2)
      IF (IOPT.LT.0) GOTO 996

      CALL GPCLOSE(LUSAVE,'KEEP')
 
*---------------------------------------------------------------------*
* dirty hack for CCS finite difference:
*---------------------------------------------------------------------*
      IF (CCSTST.AND.(IOPT.GT.1)) THEN
         CALL DZERO(VEC2,NT2AM(ISYM))
      ENDIF


      RETURN
*---------------------------------------------------------------------*
* handle i/o errors:
*---------------------------------------------------------------------*
990   CONTINUE
      IF ((IOPTSV .NE. 33) .OR. (IPRINT .GT. 0)) THEN 
        WRITE (LUPRI,'(2A)')' an error occured while inquiring file ',
     &        FILEX
      END IF
      GOTO 996

991   CONTINUE
      IF ((IOPTSV .NE. 33) .OR. (IPRINT .GT. 0)) THEN 
        WRITE (LUPRI,'(2A)')' could not find required file ',FILEX
      END IF
      GOTO 996

992   CONTINUE
      WRITE (LUPRI,'(2A)') ' i/o error while rewinding file ',FILEX
      GOTO 996

993   CONTINUE
      WRITE (LUPRI,'(A)') ' Warning: error in CC_RCRSP.'
      GOTO 996

995   CONTINUE
      WRITE (LUPRI,'(2A)') ' i/o error while closing file ',FILEX
      GOTO 996

996   CONTINUE
      IF (IOPTSV .EQ. 33) THEN
         IF (IPRINT .GT. 0) THEN
           WRITE (LUPRI,'(A,I5,A)') ' IOPTSV is ',IOPTSV,
     &              ' ... program continues nevertheless.'
         END IF
         IOPT   = IOPTSV
         IF (LUSAVE .GT. 0) CALL GPCLOSE(LUSAVE,'KEEP')
9996     CONTINUE
         RETURN
      ENDIF
      GOTO 999

999   CONTINUE
      WRITE (LUPRI,'(A,A)')       ' Fatal I/O error in ',SUBRNAME
      WRITE (LUPRI,'(A,3A)')      ' LIST / LISTI   :',LIST3,' / ',
     &     LISTI
      WRITE (LUPRI,'(A,I5,A,I5)') ' IOPTSV / IOPT  :',IOPTSV,' / ',
     &     IOPT
      WRITE (LUPRI,'(A,2I5)')     ' IDXLST, IDXVEC :',IDXLST,IDXVEC
      WRITE (LUPRI,'(A,I5)')      ' unit number    :',LUSAVE
      WRITE (LUPRI,'(A,I5)')      ' returned IOSTAT:',IOS
      CALL QUIT ('Fatal i/o error in '//SUBRNAME)

*---------------------------------------------------------------------*
*                    END OF SUBROUTINE CC_RDRSP
*---------------------------------------------------------------------*
      END 
c/* Deck cc_rwpre */
      SUBROUTINE CC_RWPRE(LIST,IDXLST,ISYM,LISTI,IDXVEC,FILEX)
C---------------------------------------------------------------------*
C
C   Purpose:  generate file name and account for special treatment
C             of Cauchy vectors with zero-order Cauchy order
C
C  Christof Haettig, november 1998
C 
C---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cccr2rsp.h>
#include <ccco2rsp.h>
#include <cccl2rsp.h>
#include <cccx2rsp.h>
C
      CHARACTER SUBRNAME*(8)
      PARAMETER (SUBRNAME = 'CC_RWPRE')
C
      INTEGER IDXLST, ISYM, ISYM1, ISYM2, IDXFIL, IDXVEC
      CHARACTER LIST*(*), LISTI*(*), LIST3*(3)
C
      CHARACTER FILEX*(*)
      REAL*8  ZERO
      PARAMETER (ZERO = 0.0d0)


* external functions:
      INTEGER IDXSYM,  ILSTSYM
      INTEGER IR1TAMP, IR2TAMP
      INTEGER IL1ZETA, IL2ZETA
      INTEGER IRHSR2,  ICHI2

 
*---------------------------------------------------------------------*
*     Check symmetry and list index 
*     (ILSTSYM quits if symmetry or index out of range)
*---------------------------------------------------------------------*
      IF (ISYM .NE. ILSTSYM(LIST,IDXLST) ) THEN
        WRITE (LUPRI,'(A,A)') 'symmetry mismatch in ',SUBRNAME
        WRITE (LUPRI,'(A,A5)') 'LIST   = ',LIST
        WRITE (LUPRI,'(A,I5)') 'IDXLST = ',IDXLST
        WRITE (LUPRI,'(A,I5)') 'ISYM   = ',ISYM
        WRITE (LUPRI,'(A,I5)') 'ILSTSYM returned :',
     &       ILSTSYM(LIST,IDXLST)
        CALL QUIT('Symmetry mismatch in CC_RWPRE.')
      END IF

*---------------------------------------------------------------------*
*     Make symmetry reduced list numbers.
*---------------------------------------------------------------------*
      IF (LEN(LIST) .EQ. 2) THEN
         LIST3 = LIST(1:2)//' '
      ELSE
         LIST3 = LIST(1:3)
      END IF
      IF ( LIST(1:2).EQ.'R0' .OR. LIST(1:2).EQ.'L0' .OR.
     &     LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'L1' .OR.
     &     LIST(1:2).EQ.'R2' .OR. LIST(1:2).EQ.'L2' .OR.
     &     LIST(1:2).EQ.'R3' .OR. LIST(1:2).EQ.'L3' .OR.
     &     LIST(1:2).EQ.'R4' .OR. LIST(1:2).EQ.'L4' .OR.
     &     LIST3(1:3).EQ.'O1 '.OR. LIST3(1:3).EQ.'X1 '.OR.
     &     LIST(1:2).EQ.'O2' .OR. LIST(1:2).EQ.'X2' .OR.
     &     LIST(1:2).EQ.'O3' .OR. LIST(1:2).EQ.'X3' .OR.
     &     LIST(1:2).EQ.'O4' .OR. LIST(1:2).EQ.'X4' .OR.
     &     LIST(1:2).EQ.'F1' .OR. LIST(1:2).EQ.'F2' .OR.
     &     LIST(1:2).EQ.'F3' .OR. LIST(1:2).EQ.'F4' .OR.
     &     LIST(1:2).EQ.'FC' .OR. LIST(1:2).EQ.'D0' .OR.
     &     LIST(1:2).EQ.'RC' .OR. LIST(1:2).EQ.'LC' .OR.
     &     LIST(1:2).EQ.'RE' .OR. LIST(1:2).EQ.'LE' .OR.
     &     LIST(1:2).EQ.'E0' .OR. LIST(1:2).EQ.'BE' .OR.
     &     LIST(1:2).EQ.'N2' .OR. LIST(1:2).EQ.'BR' .OR.
     &     LIST(1:2).EQ.'FR' .OR. LIST(1:2).EQ.'M1' 
     &                                                   ) THEN

        LISTI = LIST(1:2)//'_'
      ELSE
        LISTI = LIST(1:3)
      END IF

      IDXVEC = IDXLST
      IDXFIL = IDXSYM(LIST,ISYM,IDXLST)

*---------------------------------------------------------------------*
* special treatment for cauchy vectors with zero Cauchy orders:
*---------------------------------------------------------------------*
      IF      ( LIST(1:2).EQ.'RC') THEN
       IF (ILRCAU(IDXLST).EQ.0 ) THEN
        LISTI(1:3) = 'R1_'
        IDXVEC = IR1TAMP(LRCLBL(IDXLST),.FALSE.,ZERO,ISYM)
        IDXFIL = IDXSYM('R1',ISYM,IDXVEC)
       END IF
      ELSE IF ( LIST(1:2).EQ.'FC') THEN
       IF (ILRCAU(IDXLST).EQ.0 ) THEN
        LISTI(1:3) = 'F1_'
        IDXVEC = IR1TAMP(LRCLBL(IDXLST),.FALSE.,ZERO,ISYM)
        IDXFIL = IDXSYM('F1',ISYM,IDXVEC)
       END IF
      ELSE IF ( LIST(1:2).EQ.'LC') THEN
       IF (ILC1CAU(IDXLST).EQ.0 ) THEN
        LISTI(1:3) = 'L1_'
        IDXVEC = IL1ZETA(LBLLC1(IDXLST),.FALSE.,ZERO,ISYM)
        IDXFIL = IDXSYM('L1',ISYM,IDXVEC)
       END IF
      ELSE IF ( LIST3(1:3).EQ.'CR2') THEN
       IF (ICR2CAU(IDXLST,1).EQ.0 .AND. ICR2CAU(IDXLST,2).EQ.0) THEN
        LISTI(1:3) = 'R2_'
        IDXVEC = IR2TAMP(LBLCR2(IDXLST,1),.FALSE.,ZERO,ISYM1,
     &                   LBLCR2(IDXLST,2),.FALSE.,ZERO,ISYM2)
        IDXFIL = IDXSYM('R2',ISYM,IDXVEC)
       END IF
      ELSE IF ( LIST3(1:3).EQ.'CO2') THEN
       IF (ICO2CAU(IDXLST,1).EQ.0 .AND. ICO2CAU(IDXLST,2).EQ.0) THEN
        LISTI(1:3) = 'O2_'
        IDXVEC = IRHSR2(LBLCO2(IDXLST,1),.FALSE.,ZERO,ISYM1,
     &                  LBLCO2(IDXLST,2),.FALSE.,ZERO,ISYM2)
        IDXFIL = IDXSYM('O2',ISYM,IDXVEC)
       END IF
      ELSE IF ( LIST3(1:3).EQ.'CL2') THEN
       IF (ICL2CAU(IDXLST,1).EQ.0 .AND. ICL2CAU(IDXLST,2).EQ.0) THEN
        LISTI(1:3) = 'L2_'
        IDXVEC = IL2ZETA(LBLCL2(IDXLST,1),ZERO,ISYM1,
     &                   LBLCL2(IDXLST,2),ZERO,ISYM2)
        IDXFIL = IDXSYM('L2',ISYM,IDXVEC)
       END IF
      ELSE IF ( LIST3(1:3).EQ.'CX2') THEN
       IF (ICX2CAU(IDXLST,1).EQ.0 .AND. ICX2CAU(IDXLST,2).EQ.0) THEN
        LISTI(1:3) = 'X2_'
        IDXVEC = ICHI2(LBLCX2(IDXLST,1),.FALSE.,ZERO,ISYM1,
     &                 LBLCX2(IDXLST,2),.FALSE.,ZERO,ISYM2)
        IDXFIL = IDXSYM('X2',ISYM,IDXVEC)
       END IF
      END IF

*---------------------------------------------------------------------*
*     Make symmetry-adapted filenames
*---------------------------------------------------------------------*
      IF (ISYM .EQ. 1) LISTI = LISTI(1:3)//'1'
      IF (ISYM .EQ. 2) LISTI = LISTI(1:3)//'2'
      IF (ISYM .EQ. 3) LISTI = LISTI(1:3)//'3'
      IF (ISYM .EQ. 4) LISTI = LISTI(1:3)//'4'
      IF (ISYM .EQ. 5) LISTI = LISTI(1:3)//'5'
      IF (ISYM .EQ. 6) LISTI = LISTI(1:3)//'6'
      IF (ISYM .EQ. 7) LISTI = LISTI(1:3)//'7'
      IF (ISYM .EQ. 8) LISTI = LISTI(1:3)//'8'
 
      WRITE(FILEX,'(A2,A4,1X,I3)') 'CC',LISTI(1:4), IDXFIL
      DO I = 1, 10
        IF ( FILEX(I:I) .EQ. ' ' ) FILEX(I:I) = '_'
      END DO

*---------------------------------------------------------------------*
* special file name: 
*---------------------------------------------------------------------*
      IF      ( LIST(1:2) .EQ. 'D0' ) THEN
        FILEX = 'CCDUMMY'
      END IF

      RETURN
      END 
*---------------------------------------------------------------------*
*                    END OF SUBROUTINE CC_RWPRE
*---------------------------------------------------------------------*
c/* Deck cc_rcrsp */
      SUBROUTINE CC_RCRSP(LIST,IDXLST,LUSAVE,FILEX,MODFIL,IOPT,IOPTSV)
*=====================================================================*
*
*   Purpose:  check header of a vector file
*
*   Christof Haettig 15-05-1996
*   Introduced PL1, Sonia Coriani March 2000
*=====================================================================*
       IMPLICIT NONE

#include <priunit.h>
#include <ccer1rsp.h>
#include <ccer2rsp.h>
#include <ccel1rsp.h>
#include <ccel2rsp.h>
#include <ccr1rsp.h>
#include <ccr2rsp.h>
#include <ccr3rsp.h>
#include <ccr4rsp.h>
#include <cco1rsp.h>
#include <cco2rsp.h>
#include <cco3rsp.h>
#include <cco4rsp.h>
#include <ccl1rsp.h>
#include <ccl2rsp.h>
#include <ccl3rsp.h>
#include <ccl4rsp.h>
#include <ccx1rsp.h>
#include <ccx2rsp.h>
#include <ccx3rsp.h>
#include <ccx4rsp.h>
#include <ccn2rsp.h>
#include <cclrmrsp.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cccr2rsp.h>
#include <ccco2rsp.h>
#include <cccx2rsp.h>
#include <cccl2rsp.h>
#include <ccexci.h>
#include <ccpl1rsp.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
  
      CHARACTER SUBRNAME*(8)
      PARAMETER (SUBRNAME = 'CC_RCRSP')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      INTEGER LUSAVE, MAXORDER
      PARAMETER (MAXORDER = 4)
 
      CHARACTER LIST*(*), FILEX*(*), LIST3*(3)
      INTEGER IDXLST, IOPT, IOPTSV, NVEC, IOS
 
      CHARACTER MODFIL*(10), TYPE*(3)
      CHARACTER*8 LABFIL(MAXORDER), LABLST(MAXORDER)

      LOGICAL LORXFIL(MAXORDER), LORXLST(MAXORDER)

      INTEGER ISYFIL(MAXORDER), ISYLST(MAXORDER)
      INTEGER INDFIL(MAXORDER), INDLST(MAXORDER)

#if defined (SYS_CRAY)
      REAL FRQFIL(MAXORDER), FRQLST(MAXORDER), ZERO
      REAL THRDIFF
#else
      DOUBLE PRECISION FRQFIL(MAXORDER), FRQLST(MAXORDER), ZERO
      DOUBLE PRECISION THRDIFF
#endif
      PARAMETER (ZERO = 0.0d0)
      PARAMETER (THRDIFF = 1.0d-08)

      LOGICAL LSYM, LLBL, LFRQ, LIND, LORX
      INTEGER ORDER, IOP, IBIT32
  

*=====================================================================*
* set up the information to check the header:
*=====================================================================*

      IF (LEN(LIST) .EQ. 2) THEN
         LIST3 = LIST//' '
      ELSE
         LIST3 = LIST(1:3)
      end if

*.....................................................................
* 'R0' - zeroth-order amplitudes:
*.....................................................................
      IF ( LIST(1:2).EQ.'R0' ) THEN
          ORDER = 0
*.....................................................................
* 'R1' - first-order response t amplitudes:
* 'F1' - F transformation of first-order response t amplitudes:
* 'XF1' - F transformed first-order response t amplitudes,
*         F12*R2 contribution:
* 'eO1' - Effective rhs for Cholesky first-order response t  amplitude
*.....................................................................
Chol  ELSE IF ( LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'F1') THEN
      ELSE IF ( LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'F1' .OR.
     &          LIST(1:3).EQ.'XF1' .OR. LIST(1:3).EQ.'eO1') THEN
          ORDER = 1
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQLRT(IDXLST)
            LABLST(IOP)  = LRTLBL(IDXLST)
            ISYLST(IOP)  = ISYLRT(IDXLST)
            LORXLST(IOP) = LORXLRT(IDXLST)
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'R2' - second-order response t amplitudes:
* 'F2' - F transformation of second-order response t amplitudes:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'R2' .OR. LIST(1:2).EQ.'F2') THEN
          ORDER = 2
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQR2T(IDXLST,IOP)
            LABLST(IOP)  = LBLR2T(IDXLST,IOP)
            ISYLST(IOP)  = ISYR2T(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'R3' - third-order response t amplitudes:
* 'F3' - F transformation of third-order response t amplitudes:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'R3' .OR. LIST(1:2).EQ.'F3') THEN
          ORDER = 3
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQR3T(IDXLST,IOP)
            LABLST(IOP)  = LBLR3T(IDXLST,IOP)
            ISYLST(IOP)  = ISYR3T(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'R4' - fourth-order response t amplitudes:
* 'F4' - F transformation of fourth-order response t amplitudes:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'R4' .OR. LIST(1:2).EQ.'F4') THEN
          ORDER = 4
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQR4T(IDXLST,IOP)
            LABLST(IOP)  = LBLR4T(IDXLST,IOP)
            ISYLST(IOP)  = ISYR4T(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'CR2' -  second-order right Cauchy vectors
* 'CF2' -  F-transformed second-order right Cauchy vectors
*.....................................................................
      ELSE IF (LIST3.EQ.'CR2' .OR. LIST3.EQ.'CF2') THEN
          ORDER = 2
          DO IOP = 1, ORDER
            FRQLST(IOP)  = ZERO
            LABLST(IOP)  = LBLCR2(IDXLST,IOP)
            ISYLST(IOP)  = ISYCR2(IDXLST,IOP)
            INDLST(IOP)  = ICR2CAU(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'CO2' -  RHS vector for second-order right Cauchy vectors
*.....................................................................
      ELSE IF ( LIST3.EQ.'CO2' ) THEN
          ORDER = 2
          DO IOP = 1, ORDER
            FRQLST(IOP)  = ZERO
            LABLST(IOP)  = LBLCO2(IDXLST,IOP)
            ISYLST(IOP)  = ISYCO2(IDXLST,IOP)
            INDLST(IOP)  = ICO2CAU(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'CL2' -  second-order left Cauchy vectors
*.....................................................................
      ELSE IF ( LIST3.EQ.'CL2' ) THEN
          ORDER = 2
          DO IOP = 1, ORDER
            FRQLST(IOP)  = ZERO
            LABLST(IOP)  = LBLCL2(IDXLST,IOP)
            ISYLST(IOP)  = ISYCL2(IDXLST,IOP)
            INDLST(IOP)  = ICL2CAU(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'CX2' -  second-order Cauchy eta vectors
*.....................................................................
      ELSE IF ( LIST3.EQ.'CX2' ) THEN
          ORDER = 2
          DO IOP = 1, ORDER
            FRQLST(IOP)  = ZERO
            LABLST(IOP)  = LBLCX2(IDXLST,IOP)
            ISYLST(IOP)  = ISYCX2(IDXLST,IOP)
            INDLST(IOP)  = ICX2CAU(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'RC' -  Right Cauchy vectors
* 'FC' -  F transformation of right Cauchy vectors
*.....................................................................
      ELSE IF ( LIST(1:2) .EQ. 'RC' .OR. LIST(1:2).EQ.'FC') THEN
          ORDER = 1
          DO IOP = 1, ORDER
            FRQLST(IOP)  = ZERO
            LABLST(IOP)  = LRCLBL(IDXLST)
            ISYLST(IOP)  = ISYLRC(IDXLST)
            INDLST(IOP)  = ILRCAU(IDXLST)
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'RE' -  Right eigenvectors.
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'RE' ) THEN
          ORDER = 0
          DO IOP = 1, ORDER
            FRQLST(IOP)  = ZERO
            LABLST(IOP)  = '        '
            ISYLST(IOP)  = 0
            INDLST(IOP)  = 0
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'ER1' -  Right eigenvectors first-order response
* 'EO1' -  RHS vector for right eigenvectors first-order response
*.....................................................................
      ELSE IF ( LIST3.EQ.'ER1' .OR. LIST3.EQ.'EO1' ) THEN
          ORDER = 2
          FRQLST(1)  = EIGER1(IDXLST)
          LABLST(1)  = '        '
          ISYLST(1)  = ISYSER1(IDXLST)
          LORXLST(1) = .FALSE.
          INDLST(1)  = ISTER1(IDXLST)
          FRQLST(2)  = FRQER1(IDXLST)
          LABLST(2)  = LBLER1(IDXLST)
          ISYLST(2)  = ISYOER1(IDXLST)
          LORXLST(2) = .FALSE.
          INDLST(2)  = 0
*.....................................................................
* 'ER2' -  Right eigenvectors second-order response
* 'EO2' -  RHS vector for right eigenvectors second-order response
*.....................................................................
      ELSE IF ( LIST3.EQ.'ER2' .OR. LIST3.EQ.'EO2' ) THEN
          ORDER = 3
          FRQLST(1)  = EIGER2(IDXLST)
          LABLST(1)  = '        '
          ISYLST(1)  = ISYSER2(IDXLST)
          INDLST(1)  = ISTER2(IDXLST)
          LORXLST(1) = .FALSE.
          DO IOP = 1, ORDER-1
            FRQLST(IOP+1)  = FRQER2(IDXLST,IOP)
            LABLST(IOP+1)  = LBLER2(IDXLST,IOP)
            ISYLST(IOP+1)  = ISYOER2(IDXLST,IOP)
            INDLST(IOP+1)  = 0
            LORXLST(IOP+1) = .FALSE.
          END DO
*.....................................................................
* 'EL1' -  Left eigenvectors first-order response
* 'EX1' -  RHS vector for left eigenvectors first-order response
*.....................................................................
      ELSE IF ( LIST3.EQ.'EL1' .OR. LIST3.EQ.'EX1' ) THEN
          ORDER = 2
          FRQLST(1)  = EIGEL1(IDXLST)
          LABLST(1)  = '        '
          ISYLST(1)  = ISYSEL1(IDXLST)
          INDLST(1)  = ISTEL1(IDXLST)
          LORXLST(1) = .FALSE.
          FRQLST(2)  = FRQEL1(IDXLST)
          LABLST(2)  = LBLEL1(IDXLST)
          ISYLST(2)  = ISYOEL1(IDXLST)
          INDLST(2)  = 0
          LORXLST(2) = .FALSE.
*.....................................................................
* 'EL2' -  Left eigenvectors second-order response
* 'EX2' -  RHS vector for left eigenvectors second-order response
*.....................................................................
      ELSE IF ( LIST3.EQ.'EL2' .OR. LIST3.EQ.'EX2' ) THEN
          ORDER = 3
          FRQLST(1)  = EIGEL2(IDXLST)
          LABLST(1)  = '        '
          ISYLST(1)  = ISYSEL2(IDXLST)
          INDLST(1)  = ISTEL2(IDXLST)
          LORXLST(1) = .FALSE.
          DO IOP = 1, ORDER-1
            FRQLST(IOP+1)  = FRQEL2(IDXLST,IOP)
            LABLST(IOP+1)  = LBLEL2(IDXLST,IOP)
            ISYLST(IOP+1)  = ISYOEL2(IDXLST,IOP)
            LORXLST(IOP+1) = .FALSE.
            INDLST(IOP+1)  = 0
          END DO
*.....................................................................
* 'E0' and 'BE' - zeroth-order lagrangian multipliers and rhs.
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'E0' .OR. LIST(1:2).EQ.'BE' ) THEN
          ORDER = 0
          DO IOP = 1, ORDER
            FRQLST(IOP)  = ZERO
            LABLST(IOP)  = '        '
            ISYLST(IOP)  = 0
            INDLST(IOP)  = 0
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'L0' - zeroth-order lagrangian multipliers:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'L0' ) THEN
          ORDER = 0
*.....................................................................
* 'L1' - first-order response lagrangian multipliers:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'L1' ) THEN
          ORDER = 1
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQLRZ(IDXLST)
            LABLST(IOP)  = LRZLBL(IDXLST)
            ISYLST(IOP)  = ISYLRZ(IDXLST)
            LORXLST(IOP) = LORXLRZ(IDXLST)
Cold        LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'L2' - second-order response lagrangian multiplier:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'L2' ) THEN
          ORDER = 2
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQL2(IDXLST,IOP)
            LABLST(IOP)  = LBLL2(IDXLST,IOP)
            ISYLST(IOP)  = ISYL2(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'L3' - third-order response lagrangian multiplier:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'L3' ) THEN
          ORDER = 3
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQL3(IDXLST,IOP)
            LABLST(IOP)  = LBLL3(IDXLST,IOP)
            ISYLST(IOP)  = ISYL3(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'L4' - fourth-order response lagrangian multiplier:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'L4' ) THEN
          ORDER = 4
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQL4(IDXLST,IOP)
            LABLST(IOP)  = LBLL4(IDXLST,IOP)
            ISYLST(IOP)  = ISYL4(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'LE' -  Left eigenvectors.
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'LE' ) THEN
          ORDER = 0
          DO IOP = 1, ORDER
            FRQLST(IOP)  = ZERO
            LABLST(IOP)  = '        '
            ISYLST(IOP)  = 0
            INDLST(IOP)  = 0
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'LC' - first-order left Cauchy vectors:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'LC' ) THEN
          ORDER = 1
          DO IOP = 1, ORDER
            FRQLST(IOP)  = ZERO
            LABLST(IOP)  = LBLLC1(IDXLST)
            ISYLST(IOP)  = ISYLC1(IDXLST)
            INDLST(IOP)  = ILC1CAU(IDXLST)
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'O1' -  rhs vector for first-order amplitude equations:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'O1' ) THEN
          ORDER = 1
          FRQLST(1)  = FRQO1(IDXLST)
          LABLST(1)  = LBLO1(IDXLST)
          ISYLST(1)  = ISYO1(IDXLST)
          LORXLST(1) = LORXO1(IDXLST)
          INDLST(1)  = 0
          CALL FLSHFO(LUPRI)
*.....................................................................
* 'O2' -  rhs vector for second-order amplitude equations:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'O2' ) THEN
          ORDER = 2
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQO2(IDXLST,IOP)
            LABLST(IOP)  = LBLO2(IDXLST,IOP)
            ISYLST(IOP)  = ISYO2(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'O3' -  rhs vector for third-order amplitude equations:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'O3' ) THEN
          ORDER = 3
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQO3(IDXLST,IOP)
            LABLST(IOP)  = LBLO3(IDXLST,IOP)
            ISYLST(IOP)  = ISYO3(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'O4' -  rhs vector for fourth-order amplitude equations:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'O4' ) THEN
          ORDER = 4
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQO4(IDXLST,IOP)
            LABLST(IOP)  = LBLO4(IDXLST,IOP)
            ISYLST(IOP)  = ISYO4(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'X1_'/'X1e' -  rhs for first-order lagrangian multipliers equations:
*.....................................................................
      ELSE IF ( LIST3.EQ.'X1_'.OR.LIST3.EQ.'X1e' ) THEN
          ORDER = 1
          FRQLST(1)  = FRQX1(IDXLST)
          LABLST(1)  = LBLX1(IDXLST)
          ISYLST(1)  = ISYX1(IDXLST)
          LORXLST(1) = LORXX1(IDXLST)
          INDLST(1)  = 0
          CALL FLSHFO(LUPRI)
*.....................................................................
* 'X2' -  second-order chi vectors:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'X2' ) THEN
          ORDER = 2
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQX2(IDXLST,IOP)
            LABLST(IOP)  = LBLX2(IDXLST,IOP)
            ISYLST(IOP)  = ISYX2(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'X3' -  third-order chi vectors:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'X3' ) THEN
          ORDER = 3
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQX3(IDXLST,IOP)
            LABLST(IOP)  = LBLX3(IDXLST,IOP)
            ISYLST(IOP)  = ISYX3(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'X4' -  fourth-order chi vectors:
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'X4' ) THEN
          ORDER = 4
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQX4(IDXLST,IOP)
            LABLST(IOP)  = LBLX4(IDXLST,IOP)
            ISYLST(IOP)  = ISYX4(IDXLST,IOP)
            LORXLST(IOP) = .FALSE.
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'M1' and 'FR' - first-order response 
*.....................................................................
      ELSE IF ((LIST(1:2).EQ.'M1').OR.(LIST(1:2).EQ.'FR')) THEN
          ORDER = 1
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQLRM(IDXLST)
            LABLST(IOP)  = '        '
            ISYLST(IOP)  = ISYLRM(IDXLST)
            INDLST(IOP)  = ILRM(IDXLST)
            LORXLST(IOP) = .FALSE.
          END DO
*.....................................................................
* 'N2' and 'LBR' - first-order response 
*.....................................................................
      ELSE IF ((LIST(1:2).EQ.'N2').OR.(LIST(1:2).EQ.'BR')) THEN
          ORDER = 2
          FRQLST(1)  = FRQIN2(IDXLST)
          LABLST(1)  = '        '
          ISYLST(1)  = ISYIN2(IDXLST)
          INDLST(1)  = IIN2(IDXLST)
          LORXLST(1) = .FALSE.
          FRQLST(2)  = FRQFN2(IDXLST)
          LABLST(2)  = '        '
          ISYLST(2)  = ISYFN2(IDXLST)
          INDLST(2)  = IFN2(IDXLST)
          LORXLST(2) = .FALSE.
*.....................................................................
* 'PL1' -  Projected left eigenvectors first-order response (sonia)
*.....................................................................
      ELSE IF ( LIST3.EQ.'PL1' ) THEN
          !info expected to be on header as obtained from the list
          ORDER = 1
          DO IOP = 1, ORDER
            FRQLST(IOP)  = FRQPL1(IDXLST)
            LABLST(IOP)  = LBLPL1(IDXLST)
            ISYLST(IOP)  = ISYPL1(IDXLST)
            LORXLST(IOP) = LORXPL1(IDXLST)
            INDLST(IOP)  = 0
          END DO
*.....................................................................
* 'D0' - dummy vector
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'D0' ) THEN
          ORDER = 0
*.....................................................................
* unknown list:
*.....................................................................
      ELSE
         CALL QUIT ('Unknown list '//LIST3//' in '//SUBRNAME)
      ENDIF


*=====================================================================*
* read the header and compare:
*=====================================================================*

*.....................................................................
* R0, L0, D0
*.....................................................................
      IF ( LIST(1:2).EQ.'R0' .OR. LIST(1:2).EQ.'L0' .OR.
     &     LIST(1:2).EQ.'D0'                             ) THEN
        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL
        TYPE   = LIST(1:2)//'_'

*.....................................................................
* M1, FR
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'M1' .OR. LIST(1:2).EQ.'FR'      ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3), 
     &      INDFIL(1), ISYFIL(1), FRQFIL(1)

        LABFIL(1)  = LABLST(1)
        LORXFIL(1) = LORXLST(1)

*.....................................................................
* RE, LE, E0, BE (should be the same case as M1, FR)
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'RE' .OR. LIST(1:2).EQ.'LE' .OR.
     &          LIST(1:2).EQ.'E0' .OR. LIST(1:2).EQ.'BE'      ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC,MODFIL
        TYPE       = LIST(1:2)//'_'

*.....................................................................
* N2, LBR
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'N2' .OR. LIST(1:2).EQ.'BR'      ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC,MODFIL,TYPE(1:3),
     &        INDFIL(1), ISYFIL(1), FRQFIL(1), 
     &        INDFIL(2), ISYFIL(2), FRQFIL(2) 

Cold    TYPE(3:3)  = '_'
        LABFIL(1)  = LABLST(1)
        LABFIL(2)  = LABLST(2)
        LORXFIL(1) = LORXLST(1)
        LORXFIL(2) = LORXLST(2)

*.....................................................................
* 'R1 ', 'L1 ', 'O1 ', 'X1 ', 'F1 '
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'R1' .OR. LIST(1:2).EQ.'L1' .OR.
     &          LIST(1:2).EQ.'O1' .OR. LIST(1:2).EQ.'X1' .OR.
     &          LIST(1:2).EQ.'F1'                             ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3), 
     &      LABFIL(1), ISYFIL(1), FRQFIL(1), LORXFIL(1)

        INDFIL(1) = 0
      
*.....................................................................
* XF1, eO1   (Cholesky)
*.....................................................................
      ELSE IF ( LIST(1:3).EQ.'XF1' .OR. LIST(1:3).EQ.'eO1') THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3),
     &      LABFIL(1), ISYFIL(1), FRQFIL(1), LORXFIL(1)

        INDFIL(1) = 0

*.....................................................................
* ER1, EO1, EL1, EX1
*.....................................................................
      ELSE IF ( LIST3.EQ.'ER1' .OR. LIST3.EQ.'EO1' .OR.
     &          LIST3.EQ.'EL1' .OR. LIST3.EQ.'EX1'      ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3), 
     &      INDFIL(1), ISYFIL(1), FRQFIL(1),
     &      LABFIL(2), ISYFIL(2), FRQFIL(2)


        LABFIL(1)  = LABLST(1)
        INDFIL(2)  = 0
        LORXFIL(1) = LORXLST(1)
        LORXFIL(2) = LORXLST(2)
*.....................................................................
* PL1 (sonia)
*.....................................................................
      ELSE IF ( LIST3.EQ.'PL1' ) THEN
        !read the header of the PL1 vector on file
        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3),
     &      LABFIL(1), ISYFIL(1), FRQFIL(1), LORXFIL(1)

        INDFIL(1) = 0

*.....................................................................
* RC, LC, FC
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'RC' .OR. LIST(1:2).EQ.'LC' .OR.
     &          LIST(1:2).EQ.'FC'                             ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3), 
     &      LABFIL(1), ISYFIL(1), INDFIL(1)

        TYPE(3:3)  = '_'
        FRQFIL(1)  = ZERO
        LORXFIL(1) = LORXLST(1)

*.....................................................................
* CR2, CL2, CO2, CX2, CF2
*.....................................................................
      ELSE IF ( LIST3.EQ.'CR2' .OR. LIST3.EQ.'CL2' .OR.
     &          LIST3.EQ.'CO2' .OR. LIST3.EQ.'CX2' .OR.
     &          LIST3.EQ.'CF2'                              ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3), 
     &      LABFIL(1), ISYFIL(1), INDFIL(1),
     &      LABFIL(2), ISYFIL(2), INDFIL(2)

        FRQFIL(1)  = ZERO
        FRQFIL(2)  = ZERO
        LORXFIL(1) = LORXLST(1)
        LORXFIL(2) = LORXLST(2)

*.....................................................................
* R2, L2, O2, X2, F2
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'R2' .OR. LIST(1:2).EQ.'L2' .OR.
     &          LIST(1:2).EQ.'O2' .OR. LIST(1:2).EQ.'X2' .OR.
     &          LIST(1:2).EQ.'F2'                             ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3), 
     &      LABFIL(1), ISYFIL(1), FRQFIL(1),
     &      LABFIL(2), ISYFIL(2), FRQFIL(2) 

        TYPE(3:3)  = '_'
        INDFIL(1)  = 0
        INDFIL(2)  = 0
        LORXFIL(1) = LORXLST(1)
        LORXFIL(2) = LORXLST(2)

*.....................................................................
* ER2, EO2, EL2, EX2
*.....................................................................
      ELSE IF ( LIST3.EQ.'ER2' .OR. LIST3.EQ.'EO2' .OR.
     &          LIST3.EQ.'EL2' .OR. LIST3.EQ.'EX2'      ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3), 
     &      INDFIL(1), ISYFIL(1), FRQFIL(1),
     &      LABFIL(2), ISYFIL(2), FRQFIL(2),
     &      LABFIL(3), ISYFIL(3), FRQFIL(3)


        LABFIL(1)  = LABLST(1)
        INDFIL(2)  = 0
        INDFIL(3)  = 0
        LORXFIL(1) = LORXLST(1)
        LORXFIL(2) = LORXLST(2)
        LORXFIL(3) = LORXLST(3)

*.....................................................................
* R3, L3, O3, X3, F3
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'R3' .OR. LIST(1:2).EQ.'L3' .OR.
     &          LIST(1:2).EQ.'O3' .OR. LIST(1:2).EQ.'X3' .OR.
     &          LIST(1:2).EQ.'F3'                             ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3), 
     &      LABFIL(1), ISYFIL(1), FRQFIL(1),
     &      LABFIL(2), ISYFIL(2), FRQFIL(2),
     &      LABFIL(3), ISYFIL(3), FRQFIL(3) 

        TYPE(3:3)  = '_'
        INDFIL(1)  = 0
        INDFIL(2)  = 0
        INDFIL(3)  = 0
        LORXFIL(1) = LORXLST(1)
        LORXFIL(2) = LORXLST(2)
        LORXFIL(3) = LORXLST(3)

*.....................................................................
* R4, L4, O4, X4, F4
*.....................................................................
      ELSE IF ( LIST(1:2).EQ.'R4' .OR. LIST(1:2).EQ.'L4' .OR.
     &          LIST(1:2).EQ.'O4' .OR. LIST(1:2).EQ.'X4' .OR.
     &          LIST(1:2).EQ.'F4'                             ) THEN

        READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3), 
     &      LABFIL(1), ISYFIL(1), FRQFIL(1),
     &      LABFIL(2), ISYFIL(2), FRQFIL(2),
     &      LABFIL(3), ISYFIL(3), FRQFIL(3),
     &      LABFIL(4), ISYFIL(4), FRQFIL(4) 


        TYPE(3:3)  = '_'
        INDFIL(1)  = 0
        INDFIL(2)  = 0
        INDFIL(3)  = 0
        INDFIL(4)  = 0
        LORXFIL(1) = LORXLST(1)
        LORXFIL(2) = LORXLST(2)
        LORXFIL(3) = LORXLST(3)
        LORXFIL(4) = LORXLST(4)

*.....................................................................
* unknown list:
      ELSE
         CALL QUIT ('Unknown list '//LIST3//' in '//SUBRNAME)
      ENDIF

*=====================================================================*
* compare information from the header with the information from LIST:
*=====================================================================*
      IF (NVEC .LT.1 ) THEN
        WRITE (LUPRI,*) ' no vector found on file '//FILEX(1:10)
        CALL QUIT (' no vector found on file '//FILEX(1:10))
      END IF

      IF (TYPE(1:3) .NE. LIST3 ) THEN
        WRITE (LUPRI,*) ' wrong type of vector found on file '//
     &        FILEX(1:10)
        WRITE (LUPRI,*) ' on file:',TYPE(1:3),' expected:',LIST3
        CALL QUIT (' wrong type of vector found on file '//FILEX(1:10))
      END IF

      LSYM = .FALSE.
      LLBL = .FALSE.
      LFRQ = .FALSE.
      LIND = .FALSE.
      LORX = .FALSE.
      DO IOP = 1, ORDER
       IF ( ISYLST(IOP) .NE. ISYFIL(IOP) ) LSYM = .TRUE.
       IF ( LABLST(IOP) .NE. LABFIL(IOP) ) LLBL = .TRUE.
       IF ( DABS(FRQLST(IOP)-FRQFIL(IOP)) .GT. THRDIFF ) LFRQ = .TRUE.
       IF ( INDLST(IOP) .NE. INDFIL(IOP) ) LIND = .TRUE.
       IF ( LORXLST(IOP) .NEQV. LORXFIL(IOP) ) LORX = .TRUE.
      END DO

      IF ( LIND .OR. LSYM .OR. LFRQ .OR. LLBL .OR. LORX ) THEN
        IF ( LSYM ) THEN
           WRITE (LUPRI,*) 'vector for wrong symmetries on file ',FILEX
        END IF
        IF ( LFRQ ) THEN
           WRITE (LUPRI,*) 'vector for wrong frequencies on file ',FILEX
        END IF
        IF ( LLBL ) THEN
           WRITE (LUPRI,*) 'vector for wrong operators on file ',FILEX
        END IF
        IF ( LIND ) THEN
           WRITE (LUPRI,*) 'vector for wrong state/order on file ',FILEX
        END IF
        IF ( LORX ) THEN
           WRITE (LUPRI,*) 'vector with wrong relaxation flag on file ',
     &          FILEX
        END IF

        WRITE (LUPRI,'(2A)')       'FILE  :',FILEX
        WRITE (LUPRI,'(A,I5)')     'IDXLST:',IDXLST
        WRITE (LUPRI,'(A,10I5)')   'ISYLST:',(ISYLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10I5)')   'ISYFIL:',(ISYFIL(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10A10)')  'LABLST:',(LABLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10A10)')  'LABFIL:',(LABFIL(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10F10.6)')'FRQLST:',(FRQLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10F10.6)')'FRQFIL:',(FRQFIL(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10I5)')   'INDLST:',(INDLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10I5)')   'INDFIL:',(INDFIL(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10L5)')   'LORXLST:',(LORXLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10L5)')   'LORXFIL:',(LORXFIL(IOP),IOP=1,ORDER)

        IF (IOPTSV.EQ.33) THEN
          WRITE (LUPRI,*) ' IOPTSV=33... read vector nevertheless...'
        ELSE
          CALL QUIT(' wrong vector on file '//FILEX(1:10))
        END IF
      END IF

      IF (LOCDBG) THEN
        WRITE (LUPRI,'(2A)')       'FILE  :',FILEX
        WRITE (LUPRI,'(2A)')       'MODFIL:',MODFIL
        WRITE (LUPRI,'(A,I5)')     'IDXLST:',IDXLST
        WRITE (LUPRI,'(A,I5)')     'IOPT  :',IOPT
        WRITE (LUPRI,'(A,10I5)')   'ISYLST:',(ISYLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10I5)')   'ISYFIL:',(ISYFIL(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10A10)')  'LABLST:',(LABLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10A10)')  'LABFIL:',(LABFIL(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10F10.6)')'FRQLST:',(FRQLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10F10.6)')'FRQFIL:',(FRQFIL(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10I5)')   'INDLST:',(INDLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10I5)')   'INDFIL:',(INDFIL(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10L5)')   'LORXLST:',(LORXLST(IOP),IOP=1,ORDER)
        WRITE (LUPRI,'(A,10L5)')   'LORXFIL:',(LORXFIL(IOP),IOP=1,ORDER)
      END IF

*=====================================================================*
* check compatibility of IOPT with the MODFIL flag, if not reset IOPT:
*=====================================================================*
      IF (MODFIL(4:7)  .EQ.'-R12'.OR.
     *    MODFIL(5:8)  .EQ.'-R12'.OR.
     *    MODFIL(6:9)  .EQ.'-R12'.OR.
     *    MODFIL(7:10) .EQ.'-R12'.OR.
     *    MODFIL(8:10) .EQ.'-R1' .OR.
     *    MODFIL(9:10) .EQ.'-R'  .OR.
     *    MODFIL(10:10).EQ.'-'   .OR.
     *    MODFIL(4:8)  .EQ.'(R12)'.OR.
     *    MODFIL(5:9)  .EQ.'(R12)'.OR.
     *    MODFIL(6:10) .EQ.'(R12)'.OR.
     *    MODFIL(7:10) .EQ.'(R12'.OR.
     *    MODFIL(8:10) .EQ.'(R1' .OR.
     *    MODFIL(9:10) .EQ.'(R'  .OR.
     *    MODFIL(10:10).EQ.'('       ) THEN
        IBIT32=32
      ELSE
        IBIT32 = 0
      END IF

      IF (MODFIL(1:3).EQ.'SCF') THEN
        IOPT = IAND(IOPT,4)
      ELSE IF ((MODFIL(1:3).EQ.'CCS').AND.
     *      (.NOT.(MODFIL(1:4).EQ.'CCSD'))) THEN
        IOPT = IAND(IOPT,1+4)
Cholesky
C     No doubles in Cholesky CC2
      ELSE IF (MODFIL(1:3).EQ.'CC2' .AND. CHOINT) THEN
        IOPT = IAND(IOPT,5)
Cholesky
      ELSE IF (MODFIL(1:3).EQ.'MP2'.OR.MODFIL(1:3).EQ.'CC2'.OR.
     *         MODFIL(1:3).EQ.'CCD'.OR.MODFIL(1:4).EQ.'CCSD'
     *        ) THEN
        IOPT = IAND(IOPT,1+2+4+IBIT32)
      ELSE IF (MODFIL(1:3).EQ.'CC3'      .OR.
     *         MODFIL(1:8).EQ.'CCSDT-1a' .OR.
     *         MODFIL(1:8).EQ.'CCSDT-1B' .OR.
     *         MODFIL(1:7).EQ.'CCSD(T)'  .OR.
     *         MODFIL(1:5).EQ.'CC(3)'    .OR.
     *         MODFIL(1:8).EQ.'CCSDR(T)' .OR.
     *         MODFIL(1:8).EQ.'CCSDR(3)' .OR.
     *         MODFIL(1:9).EQ.'CCSDR(1A)'.OR.
     *         MODFIL(1:9).EQ.'CCSDR(1B)'
     *        ) THEN
        IOPT   = IAND(IOPT,1+2+4+8+16+IBIT32)
      ELSE
        WRITE (LUPRI,*) 'MODFIL:"',MODFIL,'"'
        CALL QUIT('Model not yet implemented in CC_RDRSP.')
      END IF


*---------------------------------------------------------------------*
* that's it:
*---------------------------------------------------------------------*
      RETURN

*---------------------------------------------------------------------*
* handle i/o errors:
*---------------------------------------------------------------------*
993   CONTINUE
      WRITE (LUPRI,'(A)') ' read error on unit LUSAVE in CC_RCRSP:'
      IOPT = -1
      RETURN

      END 
*---------------------------------------------------------------------*
*                    END OF SUBROUTINE CC_RCRSP
*---------------------------------------------------------------------*
c/* Deck cc_rvrsp */
      SUBROUTINE CC_RVRSP(ISYM, IMUL, IOPT, VEC0, LUSAVE, VEC1, VEC2)
C---------------------------------------------------------------------*
C
C   Purpose:  Read a vector with symmetry ISYM from unit LUSAVE
C             for an explanation of IOPT on input see CC_RDRSP
C
C             in case of an error IOPT is set to -IOPT on output
C
C---------------------------------------------------------------------*
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccfro.h>
C
      INTEGER LUSAVE
C
      INTEGER ISYM, IOPT, IMUL

      DIMENSION VEC0(*), VEC1(*), VEC2(*)
C
      INTEGER IOS

* check symmetry:
      IF (ISYM.LT.1 .OR. ISYM.GT.NSYM) THEN
        WRITE (LUPRI,*) 'ERROR IN CC_RVRSP: symmetry out of range:',
     &        ISYM
        CALL QUIT('ERROR IN CC_RVRSP: SYMMETRY OUT OF RANGE.')
      END IF

* read vector(s):
      ! CPHF part for orbital relaxed response 
      IF ( IOPT.EQ.4 ) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) (VEC0(I),I=1,2*NALLAI(ISYM))
      ELSE
         READ(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF

      ! singles cluster amplitudes
      IF ((IOPT.EQ.1).OR.(IOPT.EQ.3)) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) (VEC1(I),I=1,NT1AM(ISYM))
      ELSE IF (IOPT.GT.1) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF

      ! doubles cluster amplitudes
      IF ((IOPT.EQ.2).OR.(IOPT.EQ.3)) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) (VEC2(I),I=1,NT2AM(ISYM))
         IF (IMUL.EQ.3) THEN
           READ(LUSAVE,IOSTAT=IOS,ERR=993)
     &        (VEC2(NT2AM(ISYM)+I),I=1,NT2AMA(ISYM))
         ELSE
           READ(LUSAVE,IOSTAT=IOS,ERR=993)
         END IF                                     
      ELSE IF (IOPT .GT. 3) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) 
         READ(LUSAVE,IOSTAT=IOS,ERR=993)
      END IF

      ! R12 doubles cluster amplitudes
      IF (IOPT.EQ.32) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) (VEC2(I),I=1,NTR12AM(ISYM))
         IF (IMUL.EQ.3) THEN
           CALL QUIT('Triplet not yet implemented for R12 in CC_RVRSP')
C          READ(LUSAVE,IOSTAT=IOS,ERR=993)
C    &        (VEC2(NT2AM(ISYM)+I),I=1,NT2AMA(ISYM))
         ELSE
           READ(LUSAVE,IOSTAT=IOS,ERR=993)
         END IF                                     
      ELSE
         READ(LUSAVE,IOSTAT=IOS,ERR=993) 
         READ(LUSAVE,IOSTAT=IOS,ERR=993)
      END IF

      ! effective singles amplitudes for CC3 RHS vectors
      IF ((IOPT.EQ.8).OR.(IOPT.EQ.24)) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) (VEC1(I),I=1,NT1AM(ISYM))
      ELSE IF (IOPT.GT.8) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) 
      END IF

      ! effective doubles amplitudes for CC3 RHS vectors
      IF ((IOPT.EQ.16).OR.(IOPT.EQ.24)) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) (VEC2(I),I=1,NT2AM(ISYM))
         IF (IMUL.EQ.3) THEN
           READ(LUSAVE,IOSTAT=IOS,ERR=993)
     &        (VEC2(NT2AM(ISYM)+I),I=1,NT2AMA(ISYM))
         ELSE
           CONTINUE
         END IF                                     
      ELSE IF (IOPT .GT. 24) THEN
         CONTINUE
      END IF

      IF (IOPT.LT.1 .OR. IOPT .GT. 32) THEN
         WRITE (LUPRI,*) 'Illegal option in CC_RVRSP: IOPT = ',IOPT
         CALL QUIT('Triples not implemented in CC_RVRSP.')
      END IF

      RETURN

*---------------------------------------------------------------------*
* handle i/o error:
*---------------------------------------------------------------------*
993   CONTINUE
      WRITE(LUPRI,'(A)') ' read error on unit LUSAVE in CC_RVRSP:'
      WRITE(LUPRI,'(A,I5)') ' returned IOSTAT:',IOS
      WRITE(LUPRI,'(A,I5)') ' option   IOPT  :',IOPT
      CALL FLSHFO(LUPRI)
C     CALL QUIT ('fatal i/o error in CC_RVRSP')
      IOPT = -IOPT
      RETURN

      END 
*---------------------------------------------------------------------*
*                    END OF SUBROUTINE CC_RVRSP
*---------------------------------------------------------------------*
c /* deck CCPRPAO */
*=====================================================================*
       SUBROUTINE CCPRPAO(LABEL1,LONLYAO,PRPAO,IRREP,ISYM,IERR,
     &                    WORK,LWORK)
*---------------------------------------------------------------------*
*
*  Purpose: read property one-electron AO integrals from file AOPROPER
*
*       input:   LABEL1 -- search string for AOPROPER file
*                LONLYAO -- property AO integrals only for AO basis
*                           or for AO + auxiliary basis?
*                
*       output:  PRPAO -- property AO integrals in coupled cluster
*                         storage scheme (dimension of ouput vector
*                         is N2BST(IRREP)
*                IRREP -- irrep symmetry of the operator
*                ISYM  -- +1 for a symmetric operator 
*                         -1 for an antisymmetric operator
*                          0 if unknown
*                IERR  -- >0 on I/O error or if integrals not found
*                         <0 if if the absolute value of all integrals
*                            less CKMXPR(1.0d-12) and IRREP symmetry
*                            could not be determined
*                         N.B.: if IERR <> 0, IRREP,ISYM,PRPAO are
*                         undefined on output
*            
*
*  Written by Christof Haettig, April 1997.
*
*  Thomas Bondo Pedersen, April 2003:
*     Correct sign error in *ANGMOM integrals.
*
*  Christian Neiss, Dec. 2004:
*  adapted for R12
*  Note: the meaning of nbas, mbas1, mbas2,... depends on from WHERE
*  this routine is called (these common blocks are overwritten 
*  during R12 with different values):
*  within normal CC: nbas   number of AO basis functions
*                    mbas1  0 or not initialized
*                    mbas2  0 or not initialized
*                    ....
*  within R12 environment: nbas   number of AO + auxiliary basis func.
*                          mbas1  number of AO basis functions
*                          mbas2  number of auxiliary basis functions
*                          ....
*  To ensure the highest possible flexilbility of this subroutine,
*  we use the R12-specific variables throughout the routine. If the 
*  calculation is not R12, the variables are initialized with the 
*  corespondig values; in these cases it does not matter, whether
*  LONLYAO is TRUE or FALSE.
*  
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccorb.h>
#include <inftap.h>
#include <ccisao.h>
#include <dummy.h>
#include <r12int.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

#if defined (SYS_CRAY)
      REAL CKMXPR
#else
      DOUBLE PRECISION CKMXPR
#endif
      PARAMETER (CKMXPR = 1.0d-12)


* input:
      CHARACTER*8 LABEL1
      INTEGER IRREP, ISYM, IERR, LWORK

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL PRPAO(*)
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION PRPAO(*)
#endif

      CHARACTER*8 RTNLBL(2)
      LOGICAL LOPENED, LONLYAO
      INTEGER IDX, IDXI, IDXJ
      INTEGER KEND0, LEND0, KINTEG, ISYMA, ISYMB, IOLD, INEW
      INTEGER nallbas,nsquare,ntriangle,istart,iend
      INTEGER isymi,isymj,isymab,icount1,icount2,ioff
      INTEGER ibasxtd(8),nbasxtd(8),iaodis1(8,8),iaodis2(8,8)

* functions:
      INTEGER IDAMAX

*---------------------------------------------------------------------*
* open file and search integrals on AOPROPER file:
* (if file already open, do not reopen, but get unit number)
*---------------------------------------------------------------------*
      IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUPROP)


      IERR = -1
      CALL MOLLB2(LABEL1,RTNLBL,LUPROP,IERR)

      IF (IERR .NE. 0) THEN
        IERR = ABS(IERR)  ! force IERR to be positive for I/O errors
        WRITE (LUPRI,'(3A)')   'WARNING: property LABEL ',LABEL1,
     &                  ' not found on file AOPROPER.'
        WRITE (LUPRI,'(A,i5)') 'WARNING: IERR =',IERR
        WRITE (LUPRI,'(A)')    'WARNING: no integrals read.'
        CALL GPCLOSE(LUPROP,'KEEP')
        RETURN
      END IF

*---------------------------------------------------------------------*
* set dimensions
*---------------------------------------------------------------------*
      if (.NOT. CCR12) then
        do i = 1, nsym
          mbas1(i) = nbas(i)
          mbas2(i) = 0
        end do
      end if

      nallbas = 0
      do i = 1, nsym
        nbasxtd(i) = mbas1(i) + mbas2(i)
        nallbas = nallbas + nbasxtd(i) 
      end do
      nsquare = nallbas*nallbas
      ntriangle = nallbas*(nallbas+1)/2

      ibasxtd(1) = 0
      do i = 2, nsym
        ibasxtd(i) = ibasxtd(i-1) + nbasxtd(i-1)
      end do
 
      KINTEG = 1
      KEND0  = KINTEG + nsquare
      LEND0  = LWORK - KEND0

      IF (LEND0.LT.0) CALL QUIT('Insufficient memory in CCPRPAO.')

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'CCPRPAO> LABEL:',LABEL1
        WRITE (LUPRI,*) 'CCPRPAO> LONLYAO:',LONLYAO
        WRITE (LUPRI,*) 'CCPRPAO> IERR:',IERR
        WRITE (LUPRI,*) 'CCPRPAO> N2BASX:',N2BASX
        WRITE (LUPRI,*) 'CCPRPAO> NNBASX:',NNBASX
        WRITE (LUPRI,*) 'CCPRPAO> NSQAURE:',nsquare
        WRITE (LUPRI,*) 'CCPRPAO> NTRIANGLE:',ntriangle
        WRITE (LUPRI,*) 'CCPRPAO> RTNLBL(1):',RTNLBL(1)
        WRITE (LUPRI,*) 'CCPRPAO> RTNLBL(2):',RTNLBL(2)
        do i = 1, nsym
C          write(lupri,*) 'ibas(',i,') = ',ibas(i)
          write(lupri,*) 'nbas(',i,') = ',nbas(i)
C          write(lupri,*) 'ibasxtd(',i,') = ',ibasxtd(i)
          write(lupri,*) 'nbasxtd(',i,') = ',nbasxtd(i)
          write(lupri,*) 'mbas1(',i,') = ',mbas1(i)
          write(lupri,*) 'mbas2(',i,') = ',mbas2(i)
        end do
        CALL FLSHFO(LUPRI)
      END IF

*---------------------------------------------------------------------*
* read integrals and resort them to a full square matrix:
*---------------------------------------------------------------------*
      IF (RTNLBL(2).EQ.'SQUARE  ') THEN
        ISYM = 0
        CALL READT(LUPROP,nsquare,WORK(KINTEG))
        CALL GPCLOSE(LUPROP,'KEEP')
      ELSE IF (RTNLBL(2).EQ.'SYMMETRI') THEN
        ISYM = +1
        IF (LEND0.LT.ntriangle)
     &       CALL QUIT('Insufficient memory in CCPRPAO.')
        CALL READT(LUPROP,ntriangle,WORK(KEND0))
        CALL DSPTGE(NALLBAS,WORK(KEND0),WORK(KINTEG))
        CALL GPCLOSE(LUPROP,'KEEP')
      ELSE IF (RTNLBL(2).EQ.'ANTISYMM') THEN
        ISYM = -1
        IF (LEND0.LT.ntriangle)
     &       CALL QUIT('Insufficient memory in CCPRPAO.')
        CALL READT(LUPROP,ntriangle,WORK(KEND0))
        CALL DAPTGE(NALLBAS,WORK(KEND0),WORK(KINTEG))
        CALL GPCLOSE(LUPROP,'KEEP')
      ELSE
        IERR = 1
        WRITE (LUPRI,'(3A)') 'WARNING: no antisymmetry LABEL for ',
     &              LABEL1,' integrals on file AOPROPER.'
        WRITE (LUPRI,'(A)') 'WARNING: no integrals read.'
        CALL GPCLOSE(LUPROP,'KEEP')
        RETURN
      END IF

      IF (IPRINT.GT.99 .OR. LOCDBG) THEN
        CALL AROUND('CCPRPAO: '//LABEL1//'-integrals')
        CALL OUTPUT(WORK(KINTEG),1,NALLBAS,1,NALLBAS,NALLBAS,NALLBAS,
     &              1,LUPRI)
      END IF
*---------------------------------------------------------------------*
* find irrep symmetry of the operator:
*---------------------------------------------------------------------*
      IDX = IDAMAX(nsquare,WORK(KINTEG),1)
      IF ( ABS(WORK(KINTEG-1+IDX)) .GT. CKMXPR ) THEN
        IDXI  = (IDX+NALLBAS-1) / NALLBAS
        IDXJ  = IDX - (IDXI-1)*NALLBAS
C        IRREP = MULD2H(ISAO(IDXI),ISAO(IDXJ))
        do isyma = 1, nsym
c         do i = 1, nbasxtd(isyma)
c           ioff = ibasxtd(i) + i 
c           if (ioff .eq. idxi) isymi = isyma
c           if (ioff .eq. idxj) isymj = isyma
c         end do   
c         a bit simpler / more elegant: 
          istart = ibasxtd(isyma) + 1
          iend   = istart-1 + nbasxtd(isyma)
          if (idxi.ge.istart .and. idxi.le.iend) isymi = isyma
          if (idxj.ge.istart .and. idxj.le.iend) isymj = isyma
        end do
        IRREP = muld2h(isymi,isymj)
        if (locdbg) then
          write(lupri,*) 'CCPRPAO> IRREP = ',IRREP
        end if
      ELSE
        IRREP = 0
        IERR  = -1
        WRITE (LUPRI,'(3A,1P,D15.7)')
     &       'WARNING: integrals for operator with LABEL ',LABEL1,
     &               ' are smaller than ',CKMXPR
        WRITE (LUPRI,'(A)')
     &       'WARNING: irrep symmetry can not be determined.'
        CALL FLSHFO(LUPRI)
        RETURN
      END IF
      
*---------------------------------------------------------------------*
* resort integrals to coupled cluster storage:
*---------------------------------------------------------------------*
      !determine offsets   
      do isymab = 1, nsym
        icount1 = 0
        icount2 = 0
        do isymb = 1, nsym
          isyma = muld2h(isymab,isymb)
          iaodis1(isyma,isymb) = icount1
          iaodis2(isyma,isymb) = icount2
          icount1 = icount1 + mbas1(isyma)*mbas1(isymb)
          icount2 = icount2 + nbasxtd(isyma)*nbasxtd(isymb) 
c         if (locdbg) then
c           write(lupri,*) 'iaodis1(',isyma,',',isymb,') = ',
c    &        iaodis1(isyma,isymb)
c         end if
        end do
      end do

      DO ISYMA = 1, NSYM
        ISYMB = MULD2H(IRREP,ISYMA)
        if (LONLYAO) then
          DO A = 1, mbas1(ISYMA)
          DO B = 1, mbas1(ISYMB)
            IOLD = (ibasxtd(isyma)+A-1)*NALLBAS + ibasxtd(isymb) + B
            INEW = IAODIS1(ISYMB,ISYMA) + MBAS1(ISYMB)*(A-1) + B
            if (locdbg) then
              write(lupri,*) 'isyma, isymb, iold, inew = ',isyma, isymb,
     &          iold, inew
            end if
            PRPAO(INEW) = WORK(KINTEG-1+IOLD)
          END DO
          END DO
        else 
          do A = 1, nbasxtd(isyma)
          do B = 1, nbasxtd(isymb)
            IOLD = (ibasxtd(isyma)+A-1)*NALLBAS +
     &             (ibasxtd(isymb)+B)
            INEW = IAODIS2(ISYMB,ISYMA) + 
     &             (nbasxtd(isymb))*(A-1) + B
            PRPAO(INEW) = WORK(KINTEG-1+IOLD)
          end do
          end do
        end if
      END DO

      IF (LOCDBG) THEN
        call around('CCPRPAO: Resorted '//LABEL1//'-integrals')
        do isyma = 1, nsym
          isymb = muld2h(irrep,isyma)
          write(lupri,*) 'Symmetry block: ',isymb, isyma
          if (nbas(isyma).eq.0 .or. nbas(isymb).eq.0) then
            write(lupri,*) 'This symmetry is empty'
          else
            call output(prpao(1+IAODIS1(ISYMB,ISYMA)),1,nbas(isymb),1,
     &                  nbas(isyma),nbas(isymb),nbas(isyma),1,lupri)
          end if
        end do
      END IF

*---------------------------------------------------------------------*
* Correct sign error in *ANGMOM squaring:
*---------------------------------------------------------------------*
C  Taken out, because we were not sure if it is correct
C  C.H., P.J., F.P. Aarhus, April 04
C  TODO/FIXME: KR has promised to look at the signs of *all* one-el.
C              integrals. For now, the sign of rot. strengths and
C              optical rotations is changed in the CC code.
C
C     IF (LABEL1(2:7) .EQ. 'ANGMOM') THEN
C        CALL DSCAL(N2BST(IRREP),-1.0D0,PRPAO,1)
C     ENDIF

      RETURN
      END
*=====================================================================*
*                     END OF SUBROUTINE CCPRPAO                       *
*=====================================================================*
*=====================================================================*
      SUBROUTINE CC_GET_CMO(CMO)
*---------------------------------------------------------------------*
*     Purpose: read orbital coefficient vector CMO from SIRIFC file
*     Christof Haettig, spring 2001, Oslo
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <inftap.h>
#include <dummy.h>
#include <ccsdinp.h>

      INTEGER I

#if defined (SYS_CRAY) 
      REAL DDOT, CMO(NLAMDS)
#else
      DOUBLE PRECISION DDOT, CMO(NLAMDS)
#endif

*---------------------------------------------------------------------*
*     read MO coefficients from file:
*---------------------------------------------------------------------*
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUSIFC)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC) (CMO(I),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP')

      IF (DEBUG) THEN
        WRITE(LUPRI,*) 'Length of CMO (NLAMDS):',NLAMDS
        WRITE(LUPRI,*) 'Norm^2 of CMO:',DDOT(NLAMDS,CMO,1,CMO,1)
      END IF

      RETURN
      END
*=====================================================================*
*=====================================================================*
      SUBROUTINE CC_GET_ORBENE(ORBENE,DELFRO,LOCDBG,WORK,LWORK)
*---------------------------------------------------------------------*
*     Purpose: read canonical orbital energies from SIRIFC file
*     Christof Haettig, spring 2005, Friedrichstal(Baden)
*=====================================================================*
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <inftap.h>
#include <dummy.h>
#include <ccsdinp.h>

      LOGICAL DELFRO, LOCDBG
      INTEGER IDX, NORBEN, LWORK

#if defined (SYS_CRAY) 
      REAL ORBENE(NORBTS), WORK(LWORK)
#else
      DOUBLE PRECISION ORBENE(NORBTS), WORK(LWORK)
#endif

*---------------------------------------------------------------------*
*     read canonical orbital energies from SIRIFC file:
*---------------------------------------------------------------------*
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(IDX),IDX=1,NORBTS)
      CALL GPCLOSE(LUSIFC,'KEEP')

      NORBEN = NORBTS
      IF ( DELFRO .AND. (FROIMP .OR. FROEXP)) THEN
        CALL CCSD_DELFRO(WORK,WORK(NORBTS+1),LWORK-NORBTS)
        NORBEN = NORBT
      END IF
      CALL FOCK_REORDER(WORK,WORK(NORBTS+1),LWORK-NORBTS)
     
      CALL DCOPY(NORBEN,WORK,1,ORBENE,1)

      IF (LOCDBG) THEN
         CALL AROUND('Canonical orbital energies:')
         CALL OUTPUT(ORBENE,1,NORBEN,1,1,NORBEN,1,1,LUPRI)
      ENDIF

      RETURN
      END
*=====================================================================*
