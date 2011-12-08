C  /* Deck cc_wrrspd */
      SUBROUTINE CC_WRRSPD(LIST,IDXLST,ISYM,MODFIL,RLAX,DENAO,
     &                     WORK,LWORK)
C
C     Thomas Bondo Pedersen, April 2003.
C     - based on CC_WRRSP by C. Hattig.
C
C     Purpose: Write a response AO density matrix on file:
C              - the LIST and the list index IDXLST are used to 
C                construct a file name and to access the information
C                for the file header
C              - MODFIL is written to file header as a Label to
C                identify the CC model.
C              - RLAX flag is used to identify orbital relaxed 0th order
C                densities. Not used for higher-order densities where the
C                relaxation flag is generated from the list and list-index.
C
C     Implemented LISTs:
C
C     'd00' : 0th order density (IDXLST ignored).
C     'd01' : right first-order density.
C
#include <implicit.h>
      CHARACTER*(*) LIST
      CHARACTER*10  MODFIL
      LOGICAL       RLAX
      DIMENSION     DENAO(*), WORK(LWORK)
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
      PARAMETER (SECNAM = 'CC_WRRSPD')

      CHARACTER*4  LISTI
      CHARACTER*10 FILEX

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


C     Open and rewind file.
C     ---------------------

      LUSAVE = -1
      CALL GPOPEN(LUSAVE,FILEX,'UNKNOWN','SEQUENTIAL','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUSAVE,IOSTAT=IOS,ERR=992)

C     We always write only one density matrix, but for flexibility,
C     set NWRITE to be written to file header as in CC_WRRSP.
C     -------------------------------------------------------------

      NWRITE = 1

C     Write a LIST and index specific header to file.
C     -----------------------------------------------

      IF (LIST(1:3) .EQ. 'd00') THEN  ! 0th order density

         WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL, RLAX

      ELSE IF (LIST(1:3) .EQ. 'd01') THEN  ! right 1st-order density

         WRITE(LUSAVE,IOSTAT=IOS,ERR=993) NWRITE, MODFIL,
     &    LIST(1:3), LRTLBL(IDXLST), ISYLRT(IDXLST), FRQLRT(IDXLST),
     &    LORXLRT(IDXLST)

      ELSE

         WRITE(LUPRI,*) SECNAM,': ERROR: unknown LIST: ',LIST
         CALL QUIT('Unknown list '//LIST(1:3)//' in '//SECNAM)

      ENDIF

C     Write the AO density.
C     ---------------------

      WRITE(LUSAVE,IOSTAT=IOS,ERR=993) (DENAO(I),I=1,N2BST(ISYM))

C     Close file.
C     -----------

      CALL GPCLOSE(LUSAVE,'KEEP')

C     Normal execution: return.
C     -------------------------

      RETURN

C     Error branches.
C     ---------------

991   CONTINUE
      WRITE(LUPRI,'(2A)') ' An error occured while opening file ',FILEX
      GOTO 999

992   CONTINUE
      WRITE(LUPRI,'(2A)') ' I/O error while rewinding file ',FILEX
      GOTO 999

993   CONTINUE
      WRITE(LUPRI,'(2A)') ' Write error, file ',FILEX
      GOTO 999
      
995   CONTINUE
      WRITE(LUPRI,'(2A)') ' I/O error while closing file ',FILEX
      GOTO 999

999   CONTINUE
      WRITE(LUPRI,'(A,I5)') ' Unit number    :',LUSAVE
      WRITE(LUPRI,'(A,I5)') ' Returned IOSTAT:',IOS
      CALL QUIT ('Fatal I/O error in '//SECNAM)

      END
