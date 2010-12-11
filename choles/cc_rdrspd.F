C  /* Deck cc_rdrspd */
      SUBROUTINE CC_RDRSPD(LIST,IDXLST,ISYM,IOPT,MODEL,RLAX,DENAO)
C
C     Thomas Bondo Pedersen, April 2003.
C     - based on CC_RDRSP by C. Hattig.
C
C     Purpose: Read a response AO density matrix from file:
C             - the density is addressed via a list type LIST and
C               its index, IDXLST, within that list. 
C               IDXLST is calculated from IDXLST and symmetry info.
C             - it is checked that the density has the symmetry ISYM.
C             - it is checked that the model in the file header is identical
C               to MODEL.
C             - For the 0th order density the relaxation flag is checked.
C
C     IOPT = 33: try to read the density.
C                SUCCESS: IOPT returned is 1.
C                FAILURE: IOPT returned is 33.
C     For any other value of IOPT, read failure will cause execution abort.
C
C     Implemented LISTs:
C
C     'd00' : 0th order density (IDXLST ignored).
C     'd01' : right first-order density.
C
#include <implicit.h>
      CHARACTER*(*) LIST
      CHARACTER*10  MODEL
      LOGICAL       RLAX
      DIMENSION     DENAO(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccer1rsp.h>
#include <ccer2rsp.h>
#include <ccel1rsp.h>
#include <ccel2rsp.h>
#include <ccr1rsp.h>
#include <ccl1rsp.h>
#include <cco1rsp.h>
#include <ccx1rsp.h>
#include <ccr2rsp.h>
#include <ccl2rsp.h>
#include <cco2rsp.h>
#include <ccx2rsp.h>
#include <ccr3rsp.h>
#include <ccl3rsp.h>
#include <cco3rsp.h>
#include <ccx3rsp.h>
#include <ccr4rsp.h>
#include <ccl4rsp.h>
#include <cco4rsp.h>
#include <ccx4rsp.h>
#include <ccn2rsp.h>
#include <cclrmrsp.h>
#include <ccrc1rsp.h>
#include <cclc1rsp.h>
#include <cccr2rsp.h>
#include <ccco2rsp.h>
#include <cccl2rsp.h>
#include <cccx2rsp.h>
#include <ccpl1rsp.h>
#include <dummy.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_RDRSPD')

      LOGICAL LEXISTS, RFIL

      CHARACTER*3  TYPE
      CHARACTER*4  LISTI
      CHARACTER*10 FILEX, MODFIL

      INTEGER ORDER
      PARAMETER (MAXORDER = 4)
      CHARACTER*8 LABFIL(MAXORDER), LABLST(MAXORDER)
      LOGICAL LORXFIL(MAXORDER), LORXLST(MAXORDER)
      INTEGER ISYFIL(MAXORDER), ISYLST(MAXORDER)
      INTEGER INDFIL(MAXORDER), INDLST(MAXORDER)
      DIMENSION FRQFIL(MAXORDER), FRQLST(MAXORDER)

      PARAMETER (THRDIFF = 1.0D-8)

      LOGICAL LSYM, LLBL, LFRQ, LIND, LORX, LMOD

      INTEGER IDXSYM
      INTEGER ILSTSYM

C     Check symmetry and list index range.
C     ------------------------------------

      ITST = ILSTSYM(LIST,IDXLST)
      IF (ISYM .NE. ITST) THEN
         WRITE(LUPRI,*) 'Symmetry mismatch in ',SECNAM
         WRITE(LUPRI,*) 'LIST: ',LIST
         WRITE(LUPRI,*) 'Symmetry specified on input : ',ISYM
         WRITE(LUPRI,*) 'Symmetry returned by ILSTSYM: ',ITST
         CALL QUIT('Symmetry mismatch in '//SECNAM)
      ENDIF

C     Make symmetry adapted file name.
C     --------------------------------

      LISTI = LIST(1:3)

      IF (ISYM .EQ. 1) LISTI = LISTI(1:3)//'1'
      IF (ISYM .EQ. 2) LISTI = LISTI(1:3)//'2'
      IF (ISYM .EQ. 3) LISTI = LISTI(1:3)//'3'
      IF (ISYM .EQ. 4) LISTI = LISTI(1:3)//'4'
      IF (ISYM .EQ. 5) LISTI = LISTI(1:3)//'5'
      IF (ISYM .EQ. 6) LISTI = LISTI(1:3)//'6'
      IF (ISYM .EQ. 7) LISTI = LISTI(1:3)//'7'
      IF (ISYM .EQ. 8) LISTI = LISTI(1:3)//'8'
      IDXFIL = IDXSYM(LIST,ISYM,IDXLST)

      WRITE(FILEX,'(A2,A4,1X,I3)') 'CC', LISTI(1:4), IDXFIL
      DO I = 1,10
         IF (FILEX(I:I) .EQ. ' ') FILEX(I:I) = '_'
      ENDDO

C     Test if abort under fail.
C     -------------------------

      IF (IOPT .EQ. 33) THEN
         IOPTSV = IOPT
      ELSE
         IOPTSV = 0
      ENDIF
      IOPT = 1

C     Open and rewind file.
C     ---------------------

      INQUIRE(FILE=FILEX,EXIST=LEXISTS,IOSTAT=IOS,ERR=990)
      IF (.NOT. LEXISTS) GO TO 991

      LUSAVE = -1
      CALL GPOPEN(LUSAVE,FILEX,'OLD','SEQUENTIAL','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUSAVE,IOSTAT=IOS,ERR=992)

C     Set up info for checking file header.
C     -------------------------------------

      IF (LIST(1:3) .EQ. 'd00') THEN

         ORDER = 0

      ELSE IF (LIST(1:3) .EQ. 'd01') THEN

          ORDER = 1
          DO IOP = 1,ORDER
            FRQLST(IOP)  = FRQLRT(IDXLST)
            LABLST(IOP)  = LRTLBL(IDXLST)
            ISYLST(IOP)  = ISYLRT(IDXLST)
            LORXLST(IOP) = LORXLRT(IDXLST)
            INDLST(IOP)  = 0
          ENDDO

      ELSE

         WRITE(LUPRI,*) SECNAM,': ERROR: unknown LIST: ',LIST
         CALL QUIT('Unknown list '//LIST(1:3)//' in '//SECNAM)

      ENDIF

C     Read file header.
C     -----------------

      IF (LIST(1:3) .EQ. 'd00') THEN

         READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, RFIL
         TYPE = LIST(1:3)

      ELSE IF (LIST(1:3) .EQ. 'd01') THEN

         READ(LUSAVE,IOSTAT=IOS,ERR=993) NVEC, MODFIL, TYPE(1:3),
     &      LABFIL(1), ISYFIL(1), FRQFIL(1), LORXFIL(1)

        INDFIL(1) = 0

      ELSE

         WRITE(LUPRI,*) SECNAM,': ERROR: unknown LIST: ',LIST(1:3)
         CALL QUIT('Unknown list '//LIST(1:3)//' in '//SECNAM)

      ENDIF

C     Check.
C     ------

      IF (NVEC .LT. 1) THEN
         WRITE(LUPRI,*) SECNAM,':'
         WRITE(LUPRI,*) ' no density found on file '//FILEX(1:10)
         CALL QUIT('No density found on file '//FILEX(1:10)//' in '
     &             //SECNAM)
      ENDIF

      IF (TYPE(1:3) .NE. LIST(1:3) ) THEN
         WRITE(LUPRI,*) SECNAM,':'
         WRITE(LUPRI,*) ' wrong type of density found on file '//
     &                  FILEX(1:10)
         WRITE(LUPRI,*) ' on file:',TYPE(1:3),' expected:',LIST(1:3)
         CALL QUIT('Wrong type of density found on file '
     &             //FILEX(1:10))
      ENDIF

      LMOD = MODFIL(1:10) .NE. MODEL(1:10)
      LSYM = .FALSE.
      LLBL = .FALSE.
      LFRQ = .FALSE.
      LIND = .FALSE.
      LORX = .FALSE.
      IF (ORDER .EQ. 0) THEN
         LORX = RLAX .NEQV. RFIL
      ELSE
         DO IOP = 1, ORDER
            IF (ISYLST(IOP) .NE. ISYFIL(IOP)) LSYM = .TRUE.
            IF (LABLST(IOP) .NE. LABFIL(IOP)) LLBL = .TRUE.
            IF (DABS(FRQLST(IOP)-FRQFIL(IOP)) .GT. THRDIFF)
     &                                        LFRQ = .TRUE.
            IF (INDLST(IOP) .NE. INDFIL(IOP)) LIND = .TRUE.
            IF (LORXLST(IOP) .NEQV. LORXFIL(IOP)) LORX = .TRUE.
         ENDDO
      ENDIF

      IF (LIND .OR. LSYM .OR. LFRQ .OR. LLBL .OR. LORX .OR. LMOD) THEN
         IF (LSYM) THEN
            WRITE(LUPRI,*)
     &      'density for wrong symmetries on file ',FILEX
         ENDIF
         IF (LFRQ) THEN
            WRITE(LUPRI,*)
     &      'density for wrong frequencies on file ',FILEX
         ENDIF
         IF (LLBL) THEN
            WRITE(LUPRI,*)
     &      'density for wrong operators on file ',FILEX
         ENDIF
         IF (LIND) THEN
            WRITE(LUPRI,*)
     &      'density for wrong state/order on file ',FILEX
         ENDIF
         IF (LORX) THEN
            WRITE(LUPRI,*)
     &      'density with wrong relaxation flag on file ',FILEX
         ENDIF
         IF (LMOD) THEN
            WRITE(LUPRI,*)
     &      'density for wrong model on file ',FILEX
         ENDIF

         WRITE(LUPRI,'(2A)')       'FILE   :',FILEX
         WRITE(LUPRI,'(A,I5)')     'IDXLST :',IDXLST
         WRITE(LUPRI,'(A,10I5)')   'ISYLST :',(ISYLST(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(A,10I5)')   'ISYFIL :',(ISYFIL(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(A,10A10)')  'LABLST :',(LABLST(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(A,10A10)')  'LABFIL :',(LABFIL(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(A,10F10.6)')'FRQLST :',(FRQLST(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(A,10F10.6)')'FRQFIL :',(FRQFIL(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(A,10I5)')   'INDLST :',(INDLST(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(A,10I5)')   'INDFIL :',(INDFIL(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(A,10L5)')   'LORXLST:',(LORXLST(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(A,10L5)')   'LORXFIL:',(LORXFIL(IOP),IOP=1,ORDER)
         WRITE(LUPRI,'(2A)')       'MODEL  :',MODEL
         WRITE(LUPRI,'(2A)')       'MODFIL :',MODFIL

         IF (IOPTSV .EQ. 33) THEN
            IF (.NOT. LSYM) THEN
             WRITE(LUPRI,*) ' IOPTSV=33... read density nevertheless...'
            ELSE
               CALL QUIT('Wrong density on file '//FILEX(1:10)//' in '
     &                   //SECNAM)
            ENDIF
         ELSE
            CALL QUIT('Wrong density on file '//FILEX(1:10)//' in '
     &                //SECNAM)
         ENDIF

      ENDIF

C     Read density.
C     -------------

      IF (IOPT .EQ. 1) THEN
         READ(LUSAVE,IOSTAT=IOS,ERR=993) (DENAO(I),I=1,N2BST(ISYM))
      ELSE
         WRITE(LUPRI,*) 'An unexpected error ocurred in ',SECNAM
         WRITE(LUPRI,*) 'IOPT = ',IOPT,' (should be 1)'
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Close file.
C     -----------

      CALL GPCLOSE(LUSAVE,'KEEP')

C     Normal execution: return.
C     -------------------------

      RETURN

C     Error branches.
C     ---------------

990   CONTINUE
      IF ((IOPTSV .NE. 33) .OR. (IPRINT .GT. 0)) THEN
        WRITE (LUPRI,'(2A)')' an error occured while inquireing file ',
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
      WRITE (LUPRI,'(A)') ' i/o error while reading file ',FILEX
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
         IOPT = IOPTSV
         IF (LUSAVE .GT. 0) CALL GPCLOSE(LUSAVE,'KEEP')
9996     CONTINUE
         RETURN
      ENDIF
      GOTO 999

999   CONTINUE
      WRITE (LUPRI,'(A,A)')       ' Fatal I/O error in ',SECNAM
      WRITE (LUPRI,'(A,3A)')      ' LIST / LISTI   :',LIST(1:3),' / ',
     &     LISTI
      WRITE (LUPRI,'(A,I5,A,I5)') ' IOPTSV / IOPT  :',IOPTSV,' / ',
     &     IOPT
      WRITE (LUPRI,'(A,2I5)')     ' IDXLST, IDXLST :',IDXLST,IDXLST
      WRITE (LUPRI,'(A,I5)')      ' Unit number    :',LUSAVE
      WRITE (LUPRI,'(A,I5)')      ' Returned IOSTAT:',IOS
      CALL QUIT ('Fatal i/o error in '//SECNAM)

      END
