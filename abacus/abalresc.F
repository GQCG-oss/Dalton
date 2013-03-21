C
C  /* Deck lrscinp */
C This module read nesserary input for LRESC calculations.
C jim-gesc : May-2012
C
C
      SUBROUTINE LRSCINP(WORD)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "lrescinf.h"
#include "nuclei.h"
#include "abainf.h"
      PARAMETER (NTABLE = 11, D0 = 0.0D0)
      LOGICAL NEWDEF 
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION IPOINT(MXCENT)
C#include "dorps.h"
#include "spnout.h"
c#include "lrescinf.h"
C jim  aca van las opctiones para LRESC , aun no las defino todas
      DATA TABLE /'.GAUGEO','.PRINT ','.SELECT','.PRTALL','.xXXXXx',
     &            '.PARA1S','.PARA1T','.PARA3S','.PARA3T','.DIAM0S',
     &            '.DIAM1S'/
C
c#include <dorps.h>
c#include <spnout.h>
cC jim  aca van las opctiones para LRESC , aun no las defino
c      DATA TABLE /'.JIMOPT','.PRINT ','.SELECT','.LINEAR ','.QUADRA',
c     &            '.GESC  ','.NOPSO ','.SD+FC ','.SDxFC','.ISOTOP',
c     &            '.COUPLI'/
cC
c      write(lupri,*) ' LRSCOPTS :  called with : ' , WORD



      NEWDEF = (WORD .EQ. '*LROPTS')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
 100     CONTINUE
            READ (LUCMD,'(A7)') WORD
c            write(lupri,*) ' @LRSCINP : reading : ' , WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GOTO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GOTO (1,2,3,4,5,6,7,8,9,10,11), I
                  END IF
 200           CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GOTO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "', WORD,
     &               '" not recognized in LRESC.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Input keywords in LRESC Input, LRSCINP')
C   READ GAUGEO : DEFINIR ALGUNA VAR PARA GAUGEO
 1             CONTINUE
cv     write(lupri,*)' ..aca leo gaugeo ', WORD
               READ (LUCMD,*) (LRGAUG(IS), IS = 1, 3)
               DO ICENT = 1, NUCIND
c                  NAME =  NAMEX(3*ICENT)(1:4)
                  WRITE (LUPRI,'(2X,A,3X," : ",3(A1,2X,A,F15.10))')
     &            NAMEX(3*ICENT)(1:4), '1' , 'x' , CORD(1,ICENT),
     &                  '2' , 'y' , CORD(2,ICENT),
     &                  '3' , 'z' , CORD(3,ICENT)
               ENDDO
               GAUCHANG =.TRUE.
               ICHANG = ICHANG + 1
               GOTO 100
C   SET PRINT LEVEL 
 2             CONTINUE
               READ (LUCMD,*) JIMPRT
c              write(lupri,*)' ..aca leo print level :', JIMPRT
               ICHANG = ICHANG + 1
               GOTO 100
C   SET WHICH NUC TO DO LRESC 
 3             CONTINUE
cx             write(lupri,*)'@LRINP  antes de leer' 
               READ (LUCMD,*) LRATOM
               ICHANG = ICHANG + 1
cx             write(lupri,*)'LRINP  despues de leer', LRATOM  
cx             write(lupri,*)'Your selection is atom #:',LRATOM,
cx     $         'named : ', NAMN(LRATOM), 'on Molecule.mol file'
               IF (LRATOM.GT.NATOMS)THEN 
                  write(lupri,*) 'Input Error on LRSCINP' 
                  write(lupri,*) 'Your selection exceeds total ', 
     $                          'number of atoms on input file'
                  CALL QUIT('LRATOM grater than NATOMS, LRSCINP')
               END IF
               GOTO 100
C   SET WHICH NUC TO DO LRESC 
 4             CONTINUE
                 PRTALL=.TRUE. 
               GOTO 100
 5             CONTINUE
c                  do algo 
               GOTO 100
 6             CONTINUE
                  SIGMAP1S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 7             CONTINUE
                  SIGMAP1T = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 8             CONTINUE
                  SIGMAP3S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 9             CONTINUE
                  SIGMAP3T = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 10            CONTINUE
                  SIGMAD0S = .TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
 11            CONTINUE
                  SIGMAD1S=.TRUE.
                  LRESCALL = .FALSE.
                  ICHANG = ICHANG + 1
               GOTO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GOTO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     &               '" not recognized in LRESC'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT ('Illegal keyword in LRESC')
            END IF
      END IF
 300  CONTINUE
c jim-dbg ICHANG > 0 cambios en default variables
c jim-dbg ICHANG < 0 no cambios

      IF ((SIGMAD0S).AND.(SIGMAD1S).AND.(SIGMAP1S).AND.(SIGMAP1T).AND.
     &    (SIGMAP3S).AND.(SIGMAP3T)) THEN 
            ICHANG = ICHANG - 6 
            LRESCALL=.TRUE.
            Write(lupri,*) 'You asked for all LRESC corrections.'
      ENDIF
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for LRESC :',0)
         IF (GAUCHANG) THEN
           write(lupri,*)' - Remember to set GAUGEO on INTEGRAL ' //
     &     ' SECTION. GAUGEO MUST be placed on selected NUCLEUS'
           WRITE (LUPRI,*) ' - cambiaste origen de gauge .. ojo  '
         END IF
         IF (JIMPRT.NE.0) THEN
            WRITE (LUPRI,'(A,I5)')
     &       ' - Print level in lresc output routine now is : ', JIMPRT
         END IF
         IF (.NOT. LRESCALL) WRITE (LUPRI,'(A)')
     &      ' -  Not LRESC for all corrections. ' //
     &      'Only the following will be done : '
         IF (SIGMAP1S)  WRITE (LUPRI,'(A)')
     &   '  **Paramagnetic first order singlet corrections to shielding'
         IF (SIGMAP1T) WRITE (LUPRI,'(A)')
     &   '  **Paramagnetic first order triplet corrections to shielding'
         IF (SIGMAP3S) WRITE (LUPRI,'(A)')
     &   '  **Paramagnetic third order singlet corrections to shielding'
         IF (SIGMAP3T) WRITE (LUPRI,'(A)')
     &   '  **Paramagnetic third order triplet corrections to shielding'
         IF (SIGMAD0S) WRITE (LUPRI,*)
     &   '  **Diamagnetic zero order singlet corrections to shielding'
         IF (SIGMAD1S) WRITE (LUPRI,*)
     &   '  **Diamagnetic first order singlet corrections to shielding'
         IF (NUCSPI .GT. 0)
     &      WRITE (LUPRI,'(A)')' XXXXXXXXTn couplings of selected '//
     &           'nuclei with all the others will be calculated'
       WRITE(LUPRI,*)' ' 
      END IF
      IF (LRESCALL) WRITE (LUPRI,'(A)')
     &      ' LRESC Corrections for Shielding on Atom LRATOM. '

      RETURN
      END
C -----------------------------------------------------------------------
C
C  /* Deck lrscini */
      SUBROUTINE LRSCINI 
C
C     Initialize /LRSCINI/
C
#include "implicit.h"
#include "lrescinf.h"
#include "mxcent.h"
#include "cbiqr.h"
#include "abainf.h"
c jim-dbg : sacar esto de abajo
#include "priunit.h"
C
      JIMPRT    = 0
      LRATOM    = 1 
      LRGAUG(1) = 0.000
      LRGAUG(2) = 0.000
      LRGAUG(3) = 0.000
      GAUCHANG = .FALSE.
      SIGMAP1S = .FALSE. 
      SIGMAP1T = .FALSE. 
      SIGMAD1S = .FALSE. 
      SIGMAD0S = .FALSE. 
      SIGMAP3S = .FALSE. 
      SIGMAP3T = .FALSE. 
      LRESCALL = .TRUE.
C
C     Initialize /CUADRA/
C
      IPRINT = IPRDEF
      IPRQR  = IPRINT
      SKIP   = .FALSE.
      CUT    = .FALSE.
      OOTV   = .FALSE.
      THRESH = 1.D-04
      MAXITE = 60
      MXRM   = 400
      MXPHP  = 0
cs    LABAPP = 0 ?????????
      NABAPP = 0
      LBFREQ = 1
      LCFREQ = 1
      CALL DZERO (QBFREQ,NFMAX)
      CALL DZERO (QCFREQ,NFMAX) 
C
C  Init Result Matrices 
C
      Do I=1,3 
      DO J =1, 3
       LRFCAV(I,J)=0.0
       LRDIAK(I,J)=0.0
       LRANGP(I,J)=0.0
C
       LRDIAM(I,J)=0.0
       LRDIAD(I,J)=0.0
C
       LRLKIN(I,J)=0.0
       LRPSOK(I,J)=0.0
       LRPSKI(I,J)=0.0  ! this is just for debugging
       LRFCZK(I,J)=0.0
       LRSDZK(I,J)=0.0
       LRFCBS(I,J)=0.0
       LRSDBS(I,J)=0.0
      ENDDO
      ENDDO
C
      RETURN
      END
c    ---------------------------------------
C  /* Deck lrscdrv */
      SUBROUTINE LRSCDRV(WORK,LWORK)
C
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "priunit.h"
#include "wrkrsp.h"
#include "lrescinf.h"
#include "inforb.h"
C this is for print section
#include "nuclei.h"
C
C

cv      write(lupri,*)'@ LRSCDRV : '

       IF (LRESCALL) THEN
        SIGMAD0S = .TRUE.
        SIGMAD1S = .TRUE.
        SIGMAP1S = .TRUE.
        SIGMAP1T = .TRUE.
      ENDIF
     
C
C Zeroth Order, Singlet of course
C --------------------------------
C
      IF (SIGMAD0S) THEN
           CALL AVELRSC('FCAV',WORK,LWORK) ! ANDA
           CALL AVELRSC('DIAK',WORK,LWORK) ! ANDA 
           CALL AVELRSC('ANGP',WORK,LWORK) ! ANGPSO testing ! should work!
      END IF

C 
C First Order SINLGET 
C ---------------------
C
C     Paramagnetic ones
C

      IF (SIGMAP1S)  THEN
         CALL LINEARLR('PSOK',WORK,LWORK) 
C         WRITE(LUPRI,'(A)')'   Calling to SIGMAP1S for : PSKI  '//
C     &        'using A.1.B routine'  ! con mis integrales
C         WRITE(LUPRI,'(/721A1/)')('*',I=1,72)
C         CALL LINEARLR('PSKI',WORK,LWORK)
C
         CALL LINEARLR('LKIN',WORK,LWORK) ! ANDA
      END IF
C
C    Diamagnetic ones
C
      IF (SIGMAD1S) THEN
           CALL LINEARLR('DIAM',WORK,LWORK) ! ANDA
           CALL LINEARLR('DIAD',WORK,LWORK) ! ANDA
      END IF
C 
C First Order TRIPLET
C----------------------
C
C Paramagnetic, only
C
      IF (SIGMAP1T) THEN
           CALL LINEARLR('FCZK',WORK,LWORK) !ANDA
           CALL LINEARLR('SDZK',WORK,LWORK) !ANDA
           CALL LINEARLR('FCBS',WORK,LWORK) !ANDA 
           CALL LINEARLR('SDBS',WORK,LWORK) !ANDA
      END IF
C 
C Third Order SINLGET 
C ------------------------
C
      IF (SIGMAP3S) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA3S for : Paramagnetic'//
     &         ' Singlet Second Order Shielding corrections'
cv         CALL CUADRAS('LPSO')
      END IF
C 
C Third Order TRIPLET 
C -------------------
C
      IF (SIGMAP3T) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA3T for : '
cv         CALL CUADRAS('LFCO')
      END IF

C *******************************************************************
C *******************************************************************
C
C
C
C -================---------------======================
C                  Print Section
C -================---------------======================
C
C 
cc      IF(JIMPRT.GE.-2) THEN
c         PRTALL = .TRUE.
cc      ENDIF
 
      
      relall = 0.0
      SGP1T = 0.0  ! sigma 1order triplet PARAM : FcKin, SdKin, SdBso, FcBso
      SGP1S = 0.0  ! sigma 1order singlet PARAM : Lkin, Psokin
      SGD0S = 0.0  ! sigma 0rder singlet DIAM : Fc, AngPso, DiaKin
      SGD1S = 0.0  ! sigm`ia 1order singlet DIAM : DiaMv, DiaDw
C    Add constants to print all in ppm's
      Do i=1, 3
       LRFCZK(i,i) = calfa*calfa*calfa*calfa*1.0D+6* CFCZK*LRFCZK(i,i)
c        write(*,*)' fckin' , LRFCZK(i,i)
       LRSDZK(i,i) = calfa*calfa*calfa*calfa*1.0D+6* CSDZK*LRSDZK(i,i)
       LRFCBS(i,i) = calfa*calfa*calfa*calfa*1.0D+6* CFCBS*LRFCBS(i,i)
       LRSDBS(i,i) = calfa*calfa*calfa*calfa*1.0D+6* CSDBS*LRSDBS(i,i)
       LRLKIN(i,i) = calfa*calfa*calfa*calfa*1.0D+6* CLKIN*LRLKIN(i,i)
       LRPSOK(i,i) = calfa*calfa*calfa*calfa*1.0D+6* CPSOK*LRPSOK(i,i)
       LRFCAV(i,i) = calfa*calfa*calfa*calfa*1.0D+6* CFCAV*LRFCAV(i,i)
       LRDIAK(i,i) = calfa*calfa*calfa*calfa*1.0D+6* CDIAK*LRDIAK(i,i)
       LRANGP(i,i) = calfa*calfa*calfa*calfa*1.0D+6* CANGP*LRANGP(i,i)
       LRDIAD(i,i) = calfa*calfa*1.0D+6* CDIAD* LRDIAD(i,i)
       LRDIAM(i,i) = calfa*calfa*1.0D+6* CDIAM* LRDIAM(i,i)
      ENDDO
      
      DO i=1, 3
       SGP1T = SGP1T + LRFCZK(i,i) + LRSDZK(i,i)+LRFCBS(i,i)+LRSDBS(i,i)
       SGP1S = SGP1S + LRLKIN(i,i) + LRPSOK(i,i) 
       SGD0S = SGD0S + LRFCAV(i,i) + LRDIAK(i,i) + LRANGP(i,i)
       SGD1S = SGD1S + LRDIAD(i,i) + LRDIAM(i,i)
       relall= relall + SGD0S + SGD1S + SGP1S + SGP1T
      ENDDO
C
C  shielding constant is the trace of each tensor
C
         relall = relall/3.0
         SGP1T = SGP1T/3.0
         SGP1S = SGP1S/3.0
         SGD1S = SGD1S/3.0
         SGD0S = SGD0S/3.0
C
      CALL TITLER('ABACUS - RELATIVISTIC CORRECTIONS LRESC','*',124)
      write(lupri,*)
      write(lupri,*)' First and Second order relativistic '// 
     & 'corrections to nuclear shielding constant of ATOM ',NAMN(LRATOM)
      WRITE(LUPRI,*)'(For definitions on terms named ABC, refered '//
     &'to relativistic correction of sigma,'
      WRITE(LUPRI,*)' see JCP. 137, 214319(2012))'
      
C
C paramagnetic 1st singlet or triplet
C
      IF ((SIGMAP1S).OR.(SIGMAP1T)) THEN 
         WRITE(LUPRI,*)
         CALL HEADER('Paramagnetic Corrections',-1)
         IF(SIGMAP1S) THEN 
            WRITE(LUPRI,*)
            WRITE (LUPRI,'(5X,A,9X,F17.5,A)')
     &      'Second Order Singlets. Total Value : ', SGP1S,' [ppm]'
            WRITE(LUPRI,*)
            IF(PRTALL) THEN 
              aa =(LRLKIN(1,1)+LRLKIN(2,2)+LRLKIN(3,3))/3.0
              WRITE(LUPRI,*)'       Detailed info :  Diagonal'// 
     &                      ' components, for each operator'
              WRITE(LUPRI,*)'       ---------------            xx  '//
     &                       '          yy             zz          iso'
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'Lkin :  ',LRLKIN(1,1),LRLKIN(2,2),LRLKIN(3,3),aa,' [ppm]'
              aa = (LRPSOK(1,1)+LRPSOK(2,2)+LRPSOK(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'PsoKin :',LRPSOK(1,1),LRPSOK(2,2),LRPSOK(3,3),aa,' [ppm]'
            ENDIF
         ENDIF

         IF (SIGMAP1T) THEN 
            WRITE(LUPRI,*)
            WRITE (LUPRI,'(5X,A,9X,F17.5,A)')
     &      'Second Order Triplets. Total Value : ', SGP1T,'[ppm]'
            WRITE(LUPRI,*)
            IF(PRTALL) THEN
              WRITE(LUPRI,*)'       Detailed info :  Diagonal'// 
     &                      ' components, for each operator'
              WRITE(LUPRI,*)'       ---------------            xx  '//
     &                       '          yy             zz          iso'
              aa = (LRFCZK(1,1) + LRFCZK(2,2) + LRFCZK(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.2,X,F13.2,X,F13.2,F13.5,A)')
     &       'FcKin :',LRFCZK(1,1),LRFCZK(2,2),LRFCZK(3,3),aa,' [ppm]'
              aa = (LRSDZK(1,1) + LRSDZK(2,2) + LRSDZK(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'SdKin :',LRSDZK(1,1),LRSDZK(2,2),LRSDZK(3,3),aa,' [ppm]'
              aa = (LRFCBS(1,1) + LRFCBS(2,2) + LRFCBS(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'FcBso :',LRFCBS(1,1),LRFCBS(2,2),LRFCBS(3,3),aa,' [ppm]'
              aa = (LRSDBS(1,1) + LRSDBS(2,2) + LRSDBS(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'SdBso :',LRSDBS(1,1),LRSDBS(2,2),LRSDBS(3,3),aa,' [ppm]'
            ENDIF
         ENDIF
            WRITE(LUPRI,*)
            WRITE (LUPRI,'(5X,A,9X,F17.5,A)')
     &      'Total Paramagnetic Corrections :     ',
     &       (SGP1S+SGP1T),' [ppm]'

      ENDIF
C
C diamagnetic 0th or 1st order singlet
C
      IF ((SIGMAD0S).OR.(SIGMAD1S)) THEN 
         WRITE(LUPRI,*)
c         WRITE(LUPRI,'(/721A1/)')('=',I=1,35)
         CALL HEADER(' Diamagnetic Corrections :',-1)
         IF(SIGMAD0S) THEN 
            WRITE(LUPRI,*)
            WRITE (LUPRI,'(5X,A,9X,F17.5,A)')
     &      'First Order. Total Value :           ', SGD0S,' [ppm]'
            WRITE(LUPRI,*)
            IF(PRTALL) THEN
              WRITE(LUPRI,*)'       Detailed info :  Diagonal'// 
     &                      ' components, for each operator'
              WRITE(LUPRI,*)'       ---------------            xx  '//
     &                       '          yy             zz          iso'
              aa = (LRFCAV(1,1) + LRFCAV(2,2) + LRFCAV(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'Fc     :',LRFCAV(1,1),LRFCAV(2,2),LRFCAV(3,3),aa,' [ppm]'
              aa = (LRDIAK(1,1) + LRDIAK(2,2) + LRDIAK(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'DiaKin :',LRDIAK(1,1),LRDIAK(2,2),LRDIAK(3,3),aa,' [ppm]'
              aa = (LRANGP(1,1) + LRANGP(2,2) + LRANGP(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'AngPso :',LRANGP(1,1),LRANGP(2,2),LRANGP(3,3),aa,' [ppm]'
            ENDIF
         ENDIF
         IF(SIGMAD1S) THEN 
            WRITE(LUPRI,*)
            WRITE (LUPRI,'(5X,A,9X,F17.5,A)')
     &      'Second Order Singlets. Total Value : ', SGD1S,' [ppm]'
            WRITE(LUPRI,*)
            IF(PRTALL) THEN
              WRITE(LUPRI,*)'       Detailed info :  Diagonal'// 
     &                      ' components, for each operator'
              WRITE(LUPRI,*)'       ---------------            xx  '//
     &                       '          yy             zz          iso'
              aa = (LRDIAM(1,1) + LRDIAM(2,2) + LRDIAM(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'DiaMv  :',LRDIAM(1,1),LRDIAM(2,2),LRDIAM(3,3),aa,' [ppm]'
              aa = (LRDIAD(1,1) + LRDIAD(2,2) + LRDIAD(3,3))/3.0
              WRITE (LUPRI,'(10X,A,9X,F13.5,X,F13.5,X,F13.5,F13.5,A)')
     &       'DiaDw  :',LRDIAD(1,1),LRDIAD(2,2),LRDIAD(3,3),aa,' [ppm]'
            ENDIF
         ENDIF
      ENDIF
            WRITE(LUPRI,*)
            WRITE (LUPRI,'(5X,A,9X,F17.5,A)')
     &      'Total Diamagnetic Corrections :      ',
     &       (SGD0S+SGD1S),' [ppm]'
C
C All of them in detail
C
      IF (LRESCALL) THEN
         WRITE(LUPRI,'(/721A1/)')('=',I=1,35)
         WRITE (LUPRI,'(/5X,A,15X,F17.5,A)') 
     &           'Total relativistic correction: ', relall,
     &           ' [ppm]'
c         WRITE(LUPRI,*)
         WRITE(LUPRI,'(/721A1/)')('=',I=1,35)
      END IF

      WRITE(LUPRI,*) 
      write(LUPRI,*)'Note : Run RESPONSE module for third order '//
     &     'paramagnetic corrections.'
      write(LUPRI,*)'       Can use this as template. Settled for' //
     &     ' the first atom on input file, named X '
      WRITE(LUPRI,*)'          Spin Orbit :' 
      WRITE(LUPRI,*)'                  **RESPONSE'
      WRITE(LUPRI,*)'                  *QUADRATIC'
      WRITE(LUPRI,*)'                  .ISPABC'
      WRITE(LUPRI,*)'                  0    1    1'
      WRITE(LUPRI,*)'                  .APROP'
      WRITE(LUPRI,*)'                  XANGMOM'
      WRITE(LUPRI,*)'                  .APROP'
      WRITE(LUPRI,*)'                  YANGMOM'
      WRITE(LUPRI,*)'                  .APROP'
      WRITE(LUPRI,*)'                  ZANGMOM'
      WRITE(LUPRI,*)'                  .BPROP'
      WRITE(LUPRI,*)'                  FC X  01'
      WRITE(LUPRI,*)'                  .CPROP'
      WRITE(LUPRI,*)'                  X1SPNORB'
      WRITE(LUPRI,*)'                  .CPROP'
      WRITE(LUPRI,*)'                  Y1SPNORB'
      WRITE(LUPRI,*)'                  .CPROP'
      WRITE(LUPRI,*)'                  Z1SPNORB'
      WRITE(LUPRI,*)' ----------------------'
      WRITE(LUPRI,*)'           Massvelo or Darwin :'
      WRITE(LUPRI,*)'                  **RESPONSE'
      WRITE(LUPRI,*)'                  *QUADRATIC'
      WRITE(LUPRI,*)'                  .ISPABC'
      WRITE(LUPRI,*)'                  0    0    0'
      WRITE(LUPRI,*)'                  .APROP'
      WRITE(LUPRI,*)'                  XANGMOM'
      WRITE(LUPRI,*)'                  .APROP'
      WRITE(LUPRI,*)'                  YANGMOM'
      WRITE(LUPRI,*)'                  .APROP'
      WRITE(LUPRI,*)'                  ZANGMOM'
      WRITE(LUPRI,*)'              .    BPROP'
      WRITE(LUPRI,*)'                  PSO 001'
      WRITE(LUPRI,*)'                  .BPROP'
      WRITE(LUPRI,*)'                  PSO 002'
      WRITE(LUPRI,*)'                  .BPROP'
      WRITE(LUPRI,*)'                  PSO 003'
      WRITE(LUPRI,*)'                  .CPROP'
      WRITE(LUPRI,*)'                  MASSVELO   !(or DARWIN)'
       RETURN
      END
     RETURN
C...
c    ---------------------------------------
c    ---------------------------------------

