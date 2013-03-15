C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
C  /* Deck lrscinp */
C This module read nesserary input for LRESC calculations.
C jim-gesc : May-2012
C
C
      SUBROUTINE LRSCINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <lrescinf.h>
#include <nuclei.h>
#include <abainf.h>
      PARAMETER (NTABLE = 11, D0 = 0.0D0)
      LOGICAL NEWDEF 
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION IPOINT(MXCENT)
C#include <dorps.h>
#include <spnout.h>
c#include <lrescinf.h>
C jim  aca van las opctiones para LRESC , aun no las defino todas
      DATA TABLE /'.GAUGEO','.PRINT ','.SELECT','.xXXXXx','.xXXXXx',
     &            '.PARA1S','.PARA1T','.PARA3S','.PARA3T','.DIAM0S',
     &            '.DIAM1S'/
C
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
c                  do algo 
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
#include <implicit.h>
#include <lrescinf.h>
#include <mxcent.h>
#include <cbiqr.h>
#include <abainf.h>
c jim-dbg : sacar esto de abajo
#include <priunit.h>
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
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <priunit.h>
#include <wrkrsp.h>
#include "lrescinf.h"
#include <inforb.h>
C

cv      write(lupri,*)'@ LRSCDRV : '
      IF (LRESCALL) THEN
        SIGMAD0S = .TRUE.
        SIGMAD1S = .TRUE.
        SIGMAP1S = .TRUE.
        SIGMAP1T = .TRUE.
      ENDIF
C 
C First Order SINLGET 
C
      WRITE(LUPRI,*)
      IF (SIGMAP1S)  THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMAP1S for : PsoKin'
         WRITE (LUPRI,'(A)')'   in 2 steps, A.B and using inside integs'
c         write(*,*) 'nbas1' , nbas(1)
c         write(*,*) 'nbas' , nbas
c               nnaos = nbas(1)
c               nmber = nnaos*(nnaos+1)/2
c         CALL angkin(nmber,nnaos)   ! this is with ORIGINAL angkin routine
         WRITE (LUPRI,'(A)')'   Calling to SIGMAP1S for : PSOK  '//
     &        'using dalton integrals' ! con dalton integrales
         write(LUPRI,'(/721A1/)')('*',I=1,72)
         CALL LINEARLR('PSOK',WORK,LWORK) 
C
         WRITE (LUPRI,'(A)')'   Calling to SIGMAP1S for : PSKI  '//
     &        'using A.1.B routine'  ! con mis integrales
         write(LUPRI,'(/721A1/)')('*',I=1,72)
         CALL LINEARLR('PSKI',WORK,LWORK)

         WRITE (LUPRI,'(A)')'   Calling to SIGMAP1S for : Lkin'
         write(LUPRI,'(/721A1/)')('*',I=1,72)
         CALL LINEARLR('LKIN',WORK,LWORK) ! ANDA
      END IF

      IF (SIGMAD1S) THEN
         WRITE (LUPRI,'(A)')'  Calling to SIGMA1S for : Diamagnetic'// 
     &         ' Singlet First Order Shielding corrections'
         CALL LINEARLR('DIAM',WORK,LWORK) ! ANDA
         CALL LINEARLR('DIAD',WORK,LWORK) ! ANDA
      END IF
C 
C First Order TRIPLET
C
      IF (SIGMAP1T) THEN
         WRITE (LUPRI,'(A)')'  Calling to SIGMA1T for : Paramagnetic'//
     &         ' Triplet First Order Shielding corrections'
c         CALL LINEARLR('FCZK',WORK,LWORK) !ANDA
C         CALL LINEARLR('SDZK',WORK,LWORK) !ANDA
c         CALL LINEARLR('FCBS',WORK,LWORK) !ANDA 
         CALL LINEARLR('SDBS',WORK,LWORK) !ANDA
      END IF
C 
C Third Order SINLGET 
C
      IF (SIGMAP3S) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA3S for : Paramagnetic'//
     &         ' Singlet Second Order Shielding corrections'
cv         CALL CUADRAS('LPSO')
      END IF
C 
C Third Order TRIPLET 
C
      IF (SIGMAP3T) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA3T for : '
cv         CALL CUADRAS('LFCO')
      END IF
C
      IF (SIGMAD0S) THEN
         WRITE (LUPRI,'(A)')'   Calling to SIGMA0S for : Diamagnetic '//
     &         ' ZerothORder Shielding corrections ' 
         CALL AVELRSC('FCAV',WORK,LWORK) ! ANDA
         CALL AVELRSC('DIAK',WORK,LWORK) ! ANDA 
         CALL AVELRSC('ANGP',WORK,LWORK) ! ANGPSO testing ! should work!
      END IF

C
C -================---------------======================
C                  Print Section
C -================---------------======================
C
C 
      relall = 0.0
      SGP1T = 0.0  ! sigma 1order triplet PARAM : FcKin, SdKin, SdBso, FcBso
      SGP1S = 0.0  ! sigma 1order singlet PARAM : Lkin, Psokin
      SGD0S = 0.0  ! sigma 0rder singlet DIAM : Fc, AngPso, DiaKin
      SGD1S = 0.0  ! sigma 1order singlet DIAM : DiaMv, DiaDw
      DO i=1, 3
          SGP1T = SG1T + LRFCZK(i,i)+LRSDZK(i,i)+LRFCBS(i,i)+LRSDBS(i,i)
          SGP1S = SG1S  + LRLKIN(i,i)+LRPSOK(i,i) 
          SGD0S = SGD0S + LRFCAV(i,i)+LRDIAK(i,i)+LRANGP(i,i)
          SGD1S = SGD1S + LRDIAD(i,i)+LRDIAM(i,i)
          relall = relall + SGD0S+SGD1S+SGP1S+SGP1T
          write(lupri,*)' lrdiad' , LRDIAD(i,i)
          write(lupri,*)' lrdiaM' , LRDIAM(i,i)
      ENDDO

      CALL HEADER('Relativistic corrections to Nuclear Shielding',-1)
      WRITE(LUPRI,*)' Paramagnetic Corrections:' 
      IF (SIGMAP1S) THEN 
         WRITE (LUPRI,'(5X,A,9X,F17.10,A)')
     &   'First Order Singlets    ', SGP1S,' arb.units'
         WRITE (LUPRI,'(5X,A,9X,F17.10,A)')
     &   'First Order Triplets    ', SGP1T,' arb.units'
      ENDIF
      WRITE(LUPRI,*)' Diamagnetic Corrections:' 
      IF ((SIGMAD0S).OR.(SIGMAD1S)) THEN 
         WRITE (LUPRI,'(5X,A,9X,F17.10,A)')
     &   'Zeroth Order             ', SGD0S,' arb.units'
         WRITE (LUPRI,'(5X,A,9X,F17.10,A)')
     &   'First Order Singlets     ', SGD1S,'arb.units'
      ENDIF
      IF (LRESCALL) THEN
         WRITE (LUPRI,'(/5X,A,F17.10,A,F6.4,A)') 
     &           'Total relativistic correction:         ', relall
c    &           ' au (', D100*RELTOT/EMCSCF,'%)'
c            WRITE (LUPRI,'(5X,A,F17.10,A//)')
c    &           'Non-relativistic + relativistic energy:',
c    &           EMCSCF + RELTOT,' au'
      END IF



c       Paramagnetic Corrections
c            First Order Triplets                    bla bla ppm
c            First Order Singlets                    ble ble ppm
c 
c       Diamagnetic Corrections
c           Zeroth Order                             bli bli ppm
c           First  Order Singlet                     blo blo ppm
c 
c       Total LRESC 0th and 1st Order Corrections    suma ppm
c 
      WRITE(LUPRI,*) 
      write(LUPRI,*)'Note : For 3rd order paramagnetic corrections ' //
     &     'run response module.'
      write(LUPRI,*)'       Can use this as template, it is set for' //
     &     ' the selected atom on INPUT FILE'
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
C...
c    ---------------------------------------
c    ---------------------------------------

