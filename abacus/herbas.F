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
C FILE : herbas.F
C
#include <single.h>
C
C*************************************************************
C* This is the procedures that controls the reading from     *
C* file, and calls all the other subroutines.                *
C* It returns a bunch of variables. It returns:              *
C*     IQM - (Highest angular quantum number) + 1            *
C*     JCO - Redundant variable that gives how many blocks   *
C*           the primitives and contraction coefficients are *
C*           given in (will be set to one) (array)           *
C*     NUC - Number of primitives for a given                *
C*           (quantum number + 1)                            *
C*     NRC - Number of columns with contraction coefficients *
C*           for a given (quantum number + 1).               *
C*     SEG - Gives whether the contraction coefficients are  *
C*           segmented (for a given (quantum number + 1),    *
C*           (array).                                        *
C*     ALPHA - Array with the primitives (for a given        *
C*             (quantum number + 1)).                        *
C*     CPRIMU - Matrix (3D) with the NOT normalized          *
C*              contraction coefficients (for a given        *
C*              (quantum number + 1)).                       *
C*     CPRIM - Like CPRIMU, but this time the contraction    *
C*             coefficients are normalized.                  *
C*     NBLOCK - total number of blocks; here always same as  *
C*              IQM, but in GTOINP the user may divide any   *
C*              IQM value in several blocks.                 *
C*     KAOVEC, KPRIM, Q, DSM are max (angular quantum number *
C*             + 1), max number of primitives, nuclear       *
C*             charge, highest number to be accepted as zero *
C*     BASREF - documentation lines for basis set            *
C*     QEFF is eff. nuclear charge which is .le. Q,          *
C*             Q-QEFF will be described with ECP             *
C*************************************************************
C
      SUBROUTINE BASLIB (IQM, JCO, NUC, NRC, SEG, ALPHA, CPRIM,
     &                   CPRIMU, NBLOCK, ISGEN, KAOVEC, KPRIM, 
     &                   Q, QEFF, DSM, UNCONT, BASNAM, BASREF, IPRBAS)
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      LOGICAL   UNCONT, SEG, READNU_CALLED
      DIMENSION JCO(MXAQN), NUC(KAOVEC), NRC(KAOVEC), SEG(KAOVEC),
     +          ALPHA(KPRIM, KAOVEC), CPRIMU(KPRIM, KPRIM, KAOVEC), 
     +          CPRIM(KPRIM,KPRIM, KAOVEC), ISGEN(KAOVEC)

#ifdef PRG_DIRAC
#include "dcbgen.h"
#include "dcbham.h"
#endif /* PRG_DIRAC */

C*************************************************************
C* Variable declarations:                                    *
C*     BASNAM - A character variable that contains the name  *
C*              of the file with the basis set.              *
C*     NEWEL - Logical variable, that is returned from       *
C*             the subroutine finpos, gives whether the next *
C*             line in basis-file is a new element or not.   *
C*     SEGEJ - This gives whether the contraction            *
C*             coefficients are segmented or not (for        *
C*             (angular quantum-number + 1))                 *
C* External procedures:                                      *
C*     SEGORB - A procedure that checks if the contraction   *
C*              coefficients (for a given (angular quantum   *
C*              number + 1)) are segmented or not.           *
C*     NRMORB - Normalizing a matrix with contraction        *
C*              coefficients                                 *
C*************************************************************
C
      CHARACTER*80 BASNAM,  BASTMP, BASSAV
      CHARACTER*80 BASREF(10)
      LOGICAL NEWEL, SEGEJ, ANO, SADLEJ, POLFUN, FOUND, BLANK, NQVD
#ifdef PRG_DIRAC
C     Define a maximum number of primitives that can go in one block
C     for UNCONTRACTED basis sets. The value given is for s-functions
C     that generate 3 small component p-functions. The maximum for
C     other types of functions is scaled down by KHK(IQM+1) that
C     gives the number of functions in the small component for that type.
      PARAMETER (MAX_IN_BLOCK=3*15)
      LOGICAL SPLIT,TEST
#endif /* PRG_DIRAC */
C
      CALL QENTER('BASLIB')
C

      BASTMP = ' '
      BASSAV = BASNAM
      INTQ   = NINT(Q)
      LUBAS  = -1
      NBLOCK = 0
      IAUG   = 0
      ANO    = .FALSE.
      SADLEJ = .FALSE.
      NQVD   = .FALSE.
      POLFUN = .FALSE.

C     control print in BASLIB, usually IPRBAS = IPREAD = IPRUSR + 1 = 1
C     from calling routine
      IPRINT = IPRBAS
      IF (IPRINT .GE. 5) THEN
         write (lupri,*) 'BASLIB: Q, QEFF, INTQ',Q,QEFF,INTQ
         write (lupri,*) 'BASLIB: BASNAM ',BASNAM
      END IF

      IF (BASNAM(1:6) .EQ. 'HUCKEL') THEN

C        by default do not print which basis set directory
C        is used for huckel (IPRBAS .EQ. 1 by default) /hjaaj
         IPRINT = IPRBAS - 1

C
C hjaaj jan 2000:
C        New Huckel based on ano-4 instead of STO-3G.
C        - Advantages: atomic orbitals are orthonormal,
C                      better description of atomic shells
C        - The number of contracted functions is based
C        on the available Huckel parameters in huckel.F
C        (it would be good to extend with 2p parameters for Li, Be etc.)
C        BASNAM is needed here but contains true basis set name,
C        this is restored from BASSAV just before RETURN.
C panor/johhe 2005:
C        Huckel guess now based on ANO-DK3 instead because it 
C        includes all elements in the periodic table.
C
         ANO = .TRUE.
C hjaaj mar 2006: use ano-4 for Z.le.36 to include
C        2p for Li, Be; 3p for Na, Mg; 4p for K, Ca
C        (these basis fu. does not exist in the minimum basis ANO-DK3)
         IF (INTQ .LE. 36) THEN
            BASTMP = 'ano-4   for Huckel'
         ELSE
            BASTMP = 'ANO-DK3 for Huckel'
         END IF
C        NOTE that ECP Huckel start guess is implemented
C        in REAANO through the Q-QEFF parameter. 
         ISTART = 7
         IF (INTQ .LE. 2) THEN
C           ... H and He
            BASNAM = 'HUCKEL    1 0 0 0'
         ELSE IF (INTQ .LE. 4) THEN
C           ... Li, Be
CHJ         BASNAM = 'HUCKEL    2 0 0 0'
            BASNAM = 'HUCKEL    2 1 0 0'
         ELSE IF (INTQ .LE. 10) THEN
C           ... B, C, N, O, F, Ne
            BASNAM = 'HUCKEL    2 1 0 0'
         ELSE IF (INTQ .LE. 12) THEN
C           ... Na, Mg
CHJ         BASNAM = 'HUCKEL    3 1 0 0'
            BASNAM = 'HUCKEL    3 2 0 0'
         ELSE IF (INTQ .LE. 18) THEN
C           ... Al, Si, P, S, Cl, Ar
CHJ-jun06   BASNAM = 'HUCKEL    3 2 0 0'
CHJ-jun06   added 3d
CKR-march 07 This is dangerous if you run with a minimal basis. Reverted.
            BASNAM = 'HUCKEL    3 2 1 0'
         ELSE IF (INTQ .LE. 20) THEN
C           ...  K, Ca
CHJ         BASNAM = 'HUCKEL    4 2 0 0'
            BASNAM = 'HUCKEL    4 3 0 0'
         ELSE IF (INTQ .LE. 30) THEN
C           ...  Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
            BASNAM = 'HUCKEL    4 2 1 0'
         ELSE IF (INTQ .LE. 36) THEN
C           ...  Ga, Ge, As, Se, Br, Kr
            BASNAM = 'HUCKEL    4 3 1 0'
         ELSE IF (INTQ .LE. 38) THEN
C           ...  Rb, Sr
            BASNAM = 'HUCKEL    5 3 1 0'
         ELSE IF (INTQ .LE. 45) THEN
C           ...  Rh
            BASNAM = 'HUCKEL    5 3 2 0'
         ELSE IF (INTQ .LE. 46) THEN
C           ...  Pd
            BASNAM = 'HUCKEL    4 3 2 0'
         ELSE IF (INTQ .LE. 48) THEN
C           ...  Cd
            BASNAM = 'HUCKEL    5 3 2 0'
         ELSE IF (INTQ .LE. 54) THEN
C           ...  Xe
            BASNAM = 'HUCKEL    5 4 2 0'
         ELSE IF (INTQ .LE. 56) THEN
C           ...  Ba
            BASNAM = 'HUCKEL    6 4 2 0'
         ELSE IF (INTQ .LE. 57) THEN
C           ...  La
            BASNAM = 'HUCKEL    6 4 3 0'
         ELSE IF (INTQ .LE. 58) THEN
C           ...  Ce
            BASNAM = 'HUCKEL    6 4 3 1'
         ELSE IF (INTQ .LE. 63) THEN
C           ...  Eu
            BASNAM = 'HUCKEL    6 4 2 1'
         ELSE IF (INTQ .LE. 64) THEN
C           ...  Gd
            BASNAM = 'HUCKEL    6 4 3 1'
         ELSE IF (INTQ .LE. 70) THEN
C           ...  Yb
            BASNAM = 'HUCKEL    6 4 2 1'
         ELSE IF (INTQ .LE. 80) THEN
C           ...  Hg
            BASNAM = 'HUCKEL    6 4 3 1'
         ELSE IF (INTQ .LE. 86) THEN
C           ...  Rn
            BASNAM = 'HUCKEL    6 5 3 1'
C Manu beg 12-02-2007
         ELSE IF (INTQ .LE. 88) THEN
C           ...  Ra
            BASNAM = 'HUCKEL    7 5 3 1'
         ELSE IF (INTQ .LE. 90) THEN
C           ...  Th
            BASNAM = 'HUCKEL    7 5 4 1'
         ELSE IF (INTQ .LE. 93) THEN
C           ...  Np
            BASNAM = 'HUCKEL    7 5 4 2'
         ELSE IF (INTQ .LE. 95) THEN
C           ...  Am
            BASNAM = 'HUCKEL    7 5 3 2'
C         ELSE IF (INTQ .EQ. 92) THEN
C           ...  Special dirty fix for U
C            BASNAM = 'HUCKEL    7 5 4 2'
         ELSE IF (INTQ .LE. 103) THEN
C           ...  Lr
            BASNAM = 'HUCKEL    7 6 4 2'
C Manu end
         ELSE 
            CALL QUIT('BASLIB error: HUCKEL not defined for Z > 103')
         END IF
C djw Feb 2005
C notation changed from 'd' to 'd-' for consistency with EMSL
      ELSE
         IF (BASNAM(3:10) .EQ. 'aug-cc-p') THEN
            BASTMP(1:78) = BASNAM(3:80)
            IF (BASNAM(1:2) .EQ. 'd-') THEN
               IAUG = 1
            ELSE IF (BASNAM(1:2) .EQ. 't-') THEN
               IAUG = 2
            ELSE IF (BASNAM(1:2) .EQ. 'q-') THEN
               IAUG = 3
            ELSE
               WRITE (LUPRI,'(/A2,A)') BASNAM(1:2),' is an unknown'//
     &              'augmentation level'
               CALL QUIT('Illegal augmentation level in HERBAS')
            END IF
         ELSE IF (BASNAM(1:7) .EQ. 'aug-ecp') THEN
            BASTMP(1:76) = BASNAM(5:80)
            IAUG = 1
         ELSE IF (BASNAM(3:9) .EQ. 'aug-ecp') THEN
            BASTMP(1:74) = BASNAM(7:80)
            IF (BASNAM(1:2) .EQ. 'd-') THEN
               IAUG = 1
            ELSE IF (BASNAM(1:2) .EQ. 't-') THEN
               IAUG = 2
            ELSE IF (BASNAM(1:2) .EQ. 'q-') THEN
               IAUG = 3
            ELSE
               WRITE (LUPRI,'(/A2,A)') BASNAM(1:2),' is an unknown'//
     &              'augmentation level'
               CALL QUIT('Illegal augmentation level in HERBAS')
            END IF
C for backward compatibility the old notation is retained
         ELSE IF (BASNAM(2:9) .EQ. 'aug-cc-p') THEN
            BASTMP(1:79) = BASNAM(2:80)
            IF (BASNAM(1:1) .EQ. 'd') THEN
               IAUG = 1
            ELSE IF (BASNAM(1:1) .EQ. 't') THEN
               IAUG = 2
            ELSE IF (BASNAM(1:1) .EQ. 'q') THEN
               IAUG = 3
            ELSE
               WRITE (LUPRI,'(/A1,A)') BASNAM(1:1),' is an unknown'//
     &              'augmentation level'
               CALL QUIT('Illegal augmentation level in HERBAS')
            END IF 
         ELSE IF (BASNAM(1:7) .EQ. 'aug-ecp') THEN 
            BASTMP(1:76) = BASNAM(5:80)                      
            IAUG = 1 
         ELSE IF (BASNAM(2:8) .EQ. 'aug-ecp') THEN           
            BASTMP(1:75) = BASNAM(6:80)
            IF (BASNAM(1:1) .EQ. 'd') THEN
               IAUG = 1
            ELSE IF (BASNAM(1:1) .EQ. 't') THEN
               IAUG = 2
            ELSE IF (BASNAM(1:1) .EQ. 'q') THEN
               IAUG = 3                                      
            ELSE
               WRITE (LUPRI,'(/A1,A)') BASNAM(1:1),' is an unknown'//
     &              'augmentation level' 
               CALL QUIT('Illegal augmentation level in HERBAS')
            END IF
C
         ELSE IF (BASNAM(1:3) .EQ. 'ano') THEN
            ANO = .TRUE.
            BASTMP(1:5) = BASNAM(1:5)
            ISTART = 6
         ELSE IF (BASNAM(1:3) .EQ. 'ANO') THEN
            ANO = .TRUE.
            BASTMP(1:7) = BASNAM(1:7)
            ISTART = 8
         ELSE IF (BASNAM(1:3) .EQ. 'raf') THEN
            ANO = .TRUE.
            BASTMP(1:5) = 'raf-r'
            ISTART = 6
         ELSE IF (BASNAM(1:11) .EQ. 'Sadlej-pVTZ') THEN
            SADLEJ = .TRUE.
            BASTMP = 'Sadlej-pVTZ'
         ELSE IF (BASNAM(1:4) .EQ. 'NQvD') THEN
            NQVD = .TRUE.
            BASTMP = 'NQvD'
         ELSE IF (BASNAM(6:8) .EQ. 'ano') THEN
            ANO = .TRUE.
            BASTMP(1:8) = BASNAM(1:8)
            ISTART = 9
         ELSE
            ISTART = INDEX(BASNAM,' ')
            BASTMP = BASNAM(1:ISTART - 1)
         END IF
      END IF
C
C     Determine if there are any user added basis functions/polarization 
C     functions
C
      IPOLST = INDEX(BASNAM,'Pol')
      IF (IPOLST .GT. 0) THEN
         POLFUN = .TRUE.
         IPOLST = IPOLST + 3
      ELSE
         POLFUN = .FALSE.
      END IF

      IF (IPRINT.GE.10) THEN
        CALL HEADER('Output from BASLIB',-1)
        WRITE(LUPRI,*) 'ISTART =',ISTART
        WRITE(LUPRI,'(A,A)') ' BASTMP = ',BASTMP
        WRITE(LUPRI,*) 'IPOLST =',IPOLST
        WRITE(LUPRI,*) 'POLFUN =',POLFUN
        WRITE(LUPRI,*) 'ANO    =',ANO
        WRITE(LUPRI,*) 'SADLEJ =',SADLEJ
        WRITE(LUPRI,*) 'NQVD   =',NQVD
      ENDIF
C      
C Finds the right element in the file. 
C
      IF (IPRINT .GT. 0) WRITE (LUPRI,'(A,I4,A)')
     &   '  Basis set file used for this atomic type with Z =',
     &   INTQ,' :'
      IF (NQVD) THEN
         CALL FINNQD(BASNAM,INTQ,NBLOCK,ALPHA,CPRIMU,CPRIM,NUC,NRC,SEG,
     &               KPRIM,KAOVEC,DSM,POLFUN,IPOLST,UNCONT,LUBAS,IPRINT)
      ELSE
         CALL FIND_ELEMENT(BASTMP,INTQ,LUBAS,IPRINT)
C
C Finds the number of primitives and columns of contraction coefficients
C
         CALL FINPOS ( NEWEL, INTEXP, INTORB, INTISG, LUBAS, IPRINT)
C
 10      CONTINUE
C
         IF ( .NOT. NEWEL) THEN
            NBLOCK = NBLOCK + 1
            IF (NBLOCK .GT. KAOVEC) THEN
               CALL QUIT('Too many AO blocks'//
     &            ' for available work memory; increase WRKMEM')
            END IF
C
C Setting ALPHA, CPRIMU AND CPRIM to zero. 
C
            CALL DZERO(ALPHA(1, NBLOCK), KPRIM)
            CALL DZERO(CPRIMU(1,1, NBLOCK), KPRIM*KPRIM)
            CALL DZERO(CPRIM (1,1, NBLOCK), KPRIM*KPRIM)
C
C Reading the primitives and contraction coefficients from file in READNU
C
            IF (ANO .OR. SADLEJ) THEN
               IQCORE = NINT(Q-QEFF)
               CALL REAANO(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU,ALPHA,
     +                     KPRIM,ANO,SADLEJ,BASNAM,ISTART,POLFUN,IPOLST,
     +                     UNCONT,LUBAS,IQCORE,IPRINT)
               IF (INTORB .EQ. 0) THEN
                  NBLOCK = NBLOCK - 1 ! why this? Can we imagine user wants 2 0 1, i.e. 2s 0p 1d ??????
                  CALL FINPOS ( NEWEL, INTEXP, INTORB,
     &                          INTISG, LUBAS, IPRINT)
                  GOTO 10
               END IF
            ELSE
Chj-aug99: if (UNCONT) force uncontracted
               IF (UNCONT) INTORB = -INTORB
               CALL READNU(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU,ALPHA,
     +                KPRIM,IAUG,DSM,BASNAM,POLFUN,IPOLST,LUBAS,IPRINT)

               READNU_CALLED = .TRUE.

            END IF
            IF (INTEXP.GT.KPRIM) GOTO 5010
            IF (INTORB.GT.KPRIM) GOTO 5020          
C   
C Checking if the matrix, with the contraction coefficients, is segmented 
C
            SEGEJ  = .TRUE.
            CALL SEGORB(SEGEJ,INTEXP,INTORB,CPRIMU(1,1,NBLOCK),KPRIM,
     +                  DSM)
            NUC(NBLOCK) = INTEXP
            NRC(NBLOCK) = INTORB
            SEG(NBLOCK) = SEGEJ
            IF (SEGEJ .AND. INTORB .GT. KAOVEC) THEN
C           in her2drv:PAOSET, each segmented is put into its own block
C           thus for UNCONT each primitive is in its own block and we
C           must have at least KAOVEC blocks
               WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *      ' >>> BASLIB error, number of AO-blocks in PAOSET  ',INTORB,
     *      '                   current maximum number MXAOVC =',KAOVEC,
     &      '  Increase MXAOVC; make depend; make'
               CALL QUIT('Too many AO-blocks for segmented basis set')
            END IF
C
C           Reorder primitive orbitals
C           (For Douglas-Kroll we need to keep exponents properly
C           sorted)
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUC(NBLOCK),
     &                  NRC(NBLOCK),SEG(NBLOCK),KPRIM,DSM)
C
C Normalizing CPRIMU for (angular quantum number + 1).
C
             CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +                   ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                   CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
C
C Define the type of the next relevant line in file
C
            CALL FINPOS ( NEWEL, INTEXP, INTORB, INTISG, LUBAS, IPRINT)
            GOTO 10
         END IF
      END IF
C
      IF (POLFUN) THEN
 32      CONTINUE
         INTEXP = 0
         INTORB = 0
         NBLOCK = NBLOCK + 1
         FOUND = .FALSE.
         ISTART = IPOLST
 33      CONTINUE
         CALL CHBLANK(BLANK,BASNAM(ISTART:))
         IF (.NOT. BLANK) THEN
            CALL FREFRM(BASNAM,ISTART,IQUANT,DUMMY,'INT',IERR)
            IF (IERR .EQ. 0) THEN
               CALL FREFRM(BASNAM,ISTART,IDUMMY,EXPON,'REA',IERR2)
               IF (IERR2 .EQ. 0) THEN
                  IF (IQUANT .EQ. NBLOCK) THEN
                     FOUND  = .TRUE.
                     INTEXP = INTEXP + 1
                     INTORB = INTORB + 1
                     ALPHA(INTEXP,NBLOCK) = EXPON
                     CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0
                  END IF
               END IF
            END IF
            IF ((IERR .EQ. 0) .AND. (IERR2 .EQ. 0)) GOTO 33
         END IF
         IF (INTEXP.GT.KPRIM) GOTO 5010
         IF (INTORB.GT.KPRIM) GOTO 5020          
         IF (FOUND) THEN
            SEGEJ  = .TRUE.
            CALL SEGORB( SEGEJ, INTEXP, INTORB, CPRIMU(1,1,NBLOCK),
     +                   KPRIM, DSM)
            NUC(NBLOCK) = INTEXP
            NRC(NBLOCK) = INTORB
            SEG(NBLOCK) = SEGEJ
            IF (SEGEJ .AND. INTORB .GT. KAOVEC) THEN
C           in her2drv:PAOSET, each segmented is put into its own block
C           thus for UNCONT each primitive is in its own block and we
C           must have at least KAOVEC blocks
               WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *      ' >>> BASLIB error, number of AO-blocks for PAOSET ',INTORB,
     *      '                   current maximum number MXAOVC =',KAOVEC,
     &      '  Increase MXAOVC; make depend; make'
               CALL QUIT('Too many AO-blocks for segmented basis set')
         END IF
C
C Normalizing CPRIMU for (angular quantum number + 1).
C
            CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +                  ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                  CPRIMU(1,1, NBLOCK), KPRIM, NBLOCK)
            GOTO 32
         ELSE
            NBLOCK = NBLOCK - 1
         END IF
      END IF
C
      IF (NBLOCK .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I5,A/3X,A)')
     &   'ERROR No basis functions found for element',INTQ,' in basis:',
     &   BASNAM
         IF (ANO) WRITE(LUPRI,'(/A/A)')
     &   'ERROR NB! For ano basis sets, user must specify number of',
     &   'ERROR     AOs in each shell (example: "basis=ano-1 3 2 1")'
         CALL QUIT('No basis functions found for this element')
      END IF
C
      IQM = NBLOCK
C
C Setting JCO to one.
C
      DO I = 1, IQM
         JCO(I) = 1
      END DO
C
C Find reference information for the basis set. 
C
      IF (IPRINT .GT. 0)
     &   CALL FINREF (BASREF,LUBAS,BASSAV)
C
C Close basis set file
C
      BASNAM = BASSAV
#ifdef PRG_DIRAC
      CLOSE (LUBAS, STATUS = 'KEEP')
#else /* PRG_DIRAC */
      CALL GPCLOSE(LUBAS,'KEEP')
#endif /* PRG_DIRAC */
      CALL QEXIT('BASLIB')
      RETURN
C
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/2A/A,I6)')
     &  ' >>> BASLIB error, number of primitives per block  ',INTEXP,
     &  '                   in basis ',BASNAM,
     &  '                   current maximum number          ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/2A/A,I6)')
     &  ' >>> BASLIB error, number of contracted functions  ',INTORB,
     &  '                   in basis ',BASNAM,
     &  '                   current maximum number          ',KPRIM
        CALL QUIT('Too many contracted  functions')
C
      END
C
C*************************************************************
C* This is the subroutine that takes care of the reading of  *
C* the primitives and the contraction coefficients from the  *
C* file (for a given format). The variables that are         *
C* transferred out are:                                      *
C*     CPRIMU - The matrix that the contraction coefficients *
C*              are put in (not normalized).                 *
C*     ALPHA - The matrix where the primitives are put in    *
C*     (Both of the variables are given for a given (quantum *
C*     number + 1))                                          *
C*************************************************************
C
C/* Deck readnu */
      SUBROUTINE READNU(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU, ALPHA, 
     +            KPRIM,IAUG,DSM,BASNAM,POLFUN,IPOLST,LUBAS,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
C Reading the exponents
C
      DIMENSION CPRIMU( KPRIM, KPRIM, KAOVEC), ALPHA(KPRIM, KAOVEC)
C**************************************************************
C* Variable declarations:                                     *
C*     SIGN   - Character that gives what kind of sentence    *
C*              it is.                                        *
C*     STRING - A character string that helps bulletproofing  *
C*              the subroutine                                *
C*     BLANK  - A logical variable that gives whether the line*
C*              is blank or not.                              *
C**************************************************************
      CHARACTER SIGN
      CHARACTER STRING*88, BASNAM*(*)
      LOGICAL BLANK, POLFUN
C
C CONEOR is the maximum number of columns of contraction coefficients on 
C the 1. line in file. If INTORB is .GT. CONEOR then the contraction
C coefficients will continue on the next line. CONORB is the max. number
C of contraction coefficients on the next lines.
C
      INTEGER CONEOR, CONORB
      PARAMETER(CONEOR = 6, CONORB = 7)
      IF (IPRINT.GE.5) THEN
       CALL HEADER('Output from READNU',-1)
       write(lupri,*) 'entering READNU, IPRINT=',IPRINT
       write(lupri,*) 'NBLOCK,INTORB:',NBLOCK,INTORB
       write(lupri,*) 'DSM=',DSM
      ENDIF

C
C Initializing a counting variable that must be 0 before loop
C
      J = 0
C
 20   CONTINUE
      IF ( J .LT. INTEXP) THEN
C
C Reading the primitive and contracted coefficients
C
         READ(LUBAS, '(A88)', IOSTAT = IOERR, ERR = 2000) STRING

         IF (IPRINT.GE.7)
     &      WRITE(LUPRI,'(A,A)') '  READNU: read STRING=',STRING

         READ (STRING, '(A1)') SIGN

         IF (IPRINT.GE.7)
     &      WRITE(LUPRI,'(A,A)') '  READNU: read SIGN=',SIGN

         IF (SIGN .EQ. ' ') THEN
C
C Calling a subroutine that checks whether the line is a blank one. If
C STRING is a blank line BLANK will return as .TRUE.
C
            CALL CHBLANK(BLANK, STRING)
            IF (.NOT. BLANK) THEN
C
C We have found a line with a primitive and coefficients             
C
               J = J + 1
C
C CONLIN returns the number of lines we have the contraction
C coefficients written on.
C
               JNTORB = ABS(INTORB)

               CALL CONLIN(JNTORB, NUMLIN, DSM, IPRINT)

               IF (INTORB .LE. 0) THEN
C                 uncontracted basis set or
C                 forced uncontracted with .UNCONT
                  READ (STRING, '(F16.9)',IOSTAT=IOS) ALPHA(J,NBLOCK)
                  IF (IOS.NE.0) CALL QUIT('Error in reading alpha')

C                 skip continuation lines for contraction coefficients
                  DO I = 2,NUMLIN
                     READ (LUBAS,*)
                  END DO
                   
                  CPRIMU(J,J,NBLOCK) = 1.0D0

                  IF (IPRINT.GE.5) write(lupri,*)
     &              'READNU read 1 - NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)',
     &              NUMLIN,J,NBLOCK,ALPHA(J,NBLOCK)

                  GOTO 20
               ELSE
C
C     Getting the format for the read statement right.
C
                  IF (NUMLIN .EQ. 1) THEN
                     KNTORB = INTORB
                  ELSE
                     KNTORB = CONEOR
                  END IF
C
C     Reading the first line with exponents and contraction coefficients
C
                  READ (STRING, '(F16.9, 6F12.9)',IOSTAT=IOS)
     +                 ALPHA(J, NBLOCK), (CPRIMU(J,I,NBLOCK),
     +                 I = 1, KNTORB)
                  IF (IOS.NE.0) THEN
                    WRITE(LUPRI,*) 'KNTORB=',KNTORB
                    WRITE(LUPRI,'(A,A)') 'STRING=',STRING
                    WRITE(LUPRI,*) 'I,J,NBLOCK:',I,J,NBLOCK
                    WRITE(LUPRI,*) 'ALPHA(J, NBLOCK)=',ALPHA(J, NBLOCK)
                    WRITE(LUPRI,*) 'CPRIMU(J,I,NBLOCK), I = 1, KNTORB:',
     &              ( CPRIMU(J,I,NBLOCK),I = 1, KNTORB )
                    CALL QUIT('Error in reading basis set...')
                  ENDIF

                  IF (IPRINT.GE.6) THEN
                    write(lupri,*)
     &              'READNU read 2 - J,NBLOCK,ALPHA(J,NBLOCK)',
     &              J,NBLOCK,ALPHA(J,NBLOCK)
                    write(lupri,*)
     &              'READNU read 2 - CPRIMU:',
     &              (CPRIMU(J,I,NBLOCK),I=1,KNTORB)
                  ENDIF
C   
C     If there are more lines with contraction coefficients, 
C     they will be read here.
C
                  DO 100 I = 2, NUMLIN
                     NUMNUM = CONEOR + (I-1)*CONORB
C
C Getting the format for the read-stat right.
C
                     KNTORB = MIN0(NUMNUM, INTORB)
 30                  CONTINUE
C
C     Making the usual safety-precautions before we read the 
C     contraction coefficients.
C
                     READ(LUBAS, '(A88)', IOSTAT = IOERR, ERR = 2000)
     +                    STRING
                     READ (STRING, '(A1)') SIGN
                     IF (SIGN .EQ. ' ') THEN
                        CALL CHBLANK(BLANK, STRING)
                        IF (.NOT. BLANK) THEN
C   
C     We now have a line with contraction coefficients.
C   
                           READ (STRING,'(F16.9,6F12.9)',IOSTAT=IOS) 
     &                          (CPRIMU(J,M,NBLOCK),
     +                          M = CONEOR + (I-2)*CONORB +1, KNTORB)
                           IF (IOS.NE.0) THEN
                             CALL QUIT('Error in reading sjjd393jdj')
                           ENDIF
                           IF (IPRINT.GE.6) THEN
                             WRITE(LUPRI,*) 
     &                       J,NBLOCK,' CPRIMU(J,M,NBLOCK):'
                             WRITE(LUPRI,*) ( CPRIMU(J,M,NBLOCK),
     &                         M = CONEOR + (I-2)*CONORB +1, KNTORB )
                           ENDIF
                        END IF
C
C     Blank line read the next one
C
                     ELSE
                        GOTO 30
                     END IF
 100              CONTINUE
                  GOTO 20
               END IF
            ELSE 
C
C blank line, read the next one.
C
               GOTO 20
            END IF
            CALL QUIT('REANU: this line should never be reached')
         ELSE
C
C Found a line with nothing sensible in it
C
            GOTO 20
         END IF
      ELSE   ! IF ( J .LT. INTEXP) THEN
C
C        For uncontracted basis sets, update INTORB.
C
         IF ( INTORB .LE. 0 ) INTORB = INTEXP
C
C   
C Read all of the primitives, and the coefficients
C If this is some kind of augmented basis set, we augment it here, kr-96
C
C Fix jan -01 VB: For certain elements some of the cc basis sets do 
C not contain exponents that are monotonically decreasing.
C Before augmenting we have to determine the two lowest exponents.
C
         IF (IAUG .GT. 0) THEN
            DO 199 KAUG = 1, IAUG
               INTEXP = INTEXP + 1
               INTORB = INTORB + 1
               IF (ALPHA(1,NBLOCK) .GT. ALPHA(2,NBLOCK)) THEN
                  EXMIN1 = ALPHA(2,NBLOCK)
                  EXMIN2 = ALPHA(1,NBLOCK)
               ELSE
                  EXMIN1 = ALPHA(1,NBLOCK)
                  EXMIN2 = ALPHA(2,NBLOCK)
               END IF
               DO 200 IEX = 3, INTEXP - 1
                  IF (ALPHA(IEX,NBLOCK) .LT. EXMIN1) THEN
                     EXMIN2 = EXMIN1
                     EXMIN1 = ALPHA(IEX,NBLOCK)
                  ELSE IF (ALPHA(IEX,NBLOCK) .LT. EXMIN2) THEN
                     EXMIN2 = ALPHA(IEX,NBLOCK)
                  END IF
 200           CONTINUE
               ALPHA(INTEXP,NBLOCK) = EXMIN1*EXMIN1/EXMIN2
               CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0
 199        CONTINUE 
         END IF

         IF (POLFUN) THEN
            ISTART = IPOLST
 33         CONTINUE
            CALL CHBLANK(BLANK,BASNAM(ISTART:))
            IF (.NOT. BLANK) THEN
               CALL FREFRM(BASNAM,ISTART,IQUANT,DUMMY,'INT',IERR)
               IF (IERR .EQ. 0) THEN
                  CALL FREFRM(BASNAM,ISTART,IDUMMY,EXPON,'REA',IERR2)
                  IF (IERR2 .EQ. 0) THEN
                     IF (IQUANT .EQ. NBLOCK) THEN
                        INTEXP = INTEXP + 1
                        INTORB = INTORB + 1
                        ALPHA(INTEXP,NBLOCK) = EXPON
                        CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0
                     END IF
                     IF (IPRINT.GE.5) write(lupri,*)
     &                  'INTEXP,NBLOCK,ALPHA',
     &                   INTEXP,NBLOCK,ALPHA(INTEXP,NBLOCK)
                  END IF
               END IF
               IF ((IERR .EQ. 0) .AND. (IERR2 .EQ. 0)) GOTO 33
            END IF
         END IF

         IF (IPRINT.GE.4) THEN
            write(lupri,*) 'READNU...NBLOCK,J=',NBLOCK,J
            write(lupri,*) (ALPHA(IX,NBLOCK),IX=1,J)
         ENDIF

      END IF
C
      RETURN
C
 2000 CONTINUE
      WRITE (LUPRI,'(/A,I5)') 'Error in reading from basis file.'//
     &     ' IOSTAT =',IOERR
      END
C/* Deck reaano */
      SUBROUTINE REAANO(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU, ALPHA, 
     +                  KPRIM,ANO,SADLEJ,BASNAM,ISTART,POLFUN,IPOLST,
     +                  UNCONT,LUBAS,IQCORE,IPRINT)
C
C     IQCORE is charge of core electrons (.eq. zero if not ECP)
C     NBLOCK = l_quantum_number + 1 (1 for s, 2 for p, etc.)
C
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
C
C Reading the exponents. Format from MOLCAS files
C
      LOGICAL ANO, SADLEJ, POLFUN, BLANK, UNCONT
      DIMENSION CPRIMU( KPRIM, KPRIM, KAOVEC), ALPHA(KPRIM, KAOVEC)
      CHARACTER BASNAM*(*)
C
C Determine number of contracted functions to be read
C
      IF (ANO) THEN
Chjaaj aug00: prepared for adding max # of contracted orbitals
C             in ano-* files (for this consistency check).
         INTORBMAX = INTORB
         IF (INTORBMAX .LE. 0) INTORBMAX = INTEXP
         CALL FREFRM(BASNAM,ISTART,INTORB,DUMMY,'INT',IERR)
         IF (INTORB .GT. INTORBMAX) THEN
            WRITE (LUPRI,'(/2A/A,I5/A,I5)')
     &      'Fatal error in REAANO for basis ',BASNAM,
     &      '- max number of contracted functions',INTORBMAX,
     &      '- you asked for                     ',INTORB
            CALL QUIT('You asked for too many contracted functions')
         END IF
      END IF
C
C If ECP calculate number of core orbitals to exclude for each l-shell
C This time we have l-shell number NBLOCK
C
      CALL ECP_LCORE(IQCORE,NBLOCK,ILOFF)
C
C Reading the primitive and contracted coefficients
C
c     MAXEXP = 8
c     NUMEXP = (INTEXP - 1)/MAXEXP
c     J = 1
      IF (INTORB .GT. 0) THEN
c        DO I = 0, NUMEXP
c           IF (I.LT.NUMEXP .OR. MOD(INTEXP,MAXEXP).EQ.0 ) THEN
c              READ (LUBAS,*) (ALPHA(K,NBLOCK), K = J,J - 1 + MAXEXP)
c           ELSE
c              READ (LUBAS,*) (ALPHA(K,NBLOCK), K = J, 
c    &              J - 1 + MOD(INTEXP,MAXEXP))
c           END IF
c           J = J + MAXEXP
c        END DO
         READ (LUBAS,*) (ALPHA(K,NBLOCK), K = 1, INTEXP)
         IF (UNCONT .AND. BASNAM(1:6) .NE. 'HUCKEL') THEN
C        ... force uncontracted
            DO I = 1, INTEXP
               READ (LUBAS,*)
               CPRIMU(I,I,NBLOCK) = 1.0D0
            END DO
            INTORB = INTEXP
         ELSE IF (ILOFF .GT. 0) THEN
C           skip core orbitals for Huckel guess for ECP
            INTORB = INTORB - ILOFF
            DO I = 1, INTEXP
               READ (LUBAS,*) (XXX, J = 1,ILOFF),
     &              (CPRIMU(I,J,NBLOCK), J = 1, INTORB)
            END DO
         ELSE
            DO I = 1, INTEXP
               READ (LUBAS,*) (CPRIMU(I,J,NBLOCK), J = 1, INTORB)
            END DO
         END IF
  
      ELSE
         READ (LUBAS,*) (XXX, K = 1, INTEXP)
         DO I = 1, INTEXP
            READ (LUBAS,*)
         END DO
      END IF
      IF (POLFUN) THEN
         IREAD = IPOLST
 33      CONTINUE
         CALL CHBLANK(BLANK,BASNAM(IREAD:))
         IF (.NOT. BLANK) THEN
            CALL FREFRM(BASNAM,IREAD,IQUANT,DUMMY,'INT',IERR)
            IF (IERR .EQ. 0) THEN
               CALL FREFRM(BASNAM,IREAD,IDUMMY,EXPON,'REA',IERR2)
               IF (IERR2 .EQ. 0) THEN
                  IF (IQUANT .EQ. NBLOCK) THEN
                     INTEXP = INTEXP + 1
                     INTORB = INTORB + 1
                     ALPHA(INTEXP,NBLOCK) = EXPON
                     CPRIMU(INTEXP,INTORB,NBLOCK) = 1.0D0
                  END IF
               END IF
            END IF
            IF ((IERR .EQ. 0) .AND. (IERR2 .EQ. 0)) GOTO 33
         END IF
      END IF
      RETURN
      END
C/* Deck chblank */
C
C*************************************************************
C* CHBLANK is a subroutine that checks whether a line is a blank*
C* line. It returns a logical variable if it is a blank one. *
C* The in/out-going variables are:                           *
C*     STRING - A character string (the string checked by the*
C*              subroutine).                                 *
C*     BLANK - The logical variable. It will be .TRUE. if the*
C*             line is a blank one.                          *
C*************************************************************
      SUBROUTINE CHBLANK(BLANK,STRING)
C
      LOGICAL BLANK
      CHARACTER*(*) STRING
C
      BLANK = .TRUE.
C
      DO J = 1, LEN(STRING)
         IF (STRING(J:J) .NE. ' ') THEN
            BLANK = .FALSE.
            GO TO 9000
         END IF
      END DO

 9000 RETURN
      END
C/* Deck conlin */
C
C*************************************************************
C* CONLIN is a subroutine that finds out how many lines the  *
C* contraction coefficients are written on.  There are 6     *
C* cc (contraction coefficients) on the 1. line, and then    *
C* there are seven on the rest of the lines. In out going    *
C* variables are :                                           *
C*      INTORB - Total number of columns of contraction      *
C*               coefficients.                               *
C*      NUMLIN - Returns the number of lines the contraction *
C*               coefficients are written on.                *
C*************************************************************
C
      SUBROUTINE CONLIN(INTORB, NUMLIN, DSM, IPRINT)
#include "implicit.h"
#include "priunit.h"
C      
C This is the parameters for the number of contraction coefficients
C on 1. and 2. line.
C
      INTEGER CONEOR, CONORB
      PARAMETER(CONEOR = 6, CONORB = 7)
C
C The intrinsic functions DBLE makes a double precision real number of an integer.
C
      B = DBLE(CONORB)
      C = DBLE(INTORB) - DBLE(CONEOR) 
C
C This finds out how many lines we have, and puts it into numlin.
C
      IF ((INTORB - CONEOR) .LE. 0) THEN
         NUMLIN = 1
      ELSE IF (DMOD(C,B) .LT. DSM) THEN
         NUMLIN = (INTORB - CONEOR)/CONORB + 1
      ELSE
         NUMLIN = (INTORB - CONEOR)/CONORB + 2
      END IF

      IF (IPRINT.GE.10) THEN
       write(lupri,'(/,2X,A)') '*** CONLIN output ***'
       write(lupri,*) 'input INTORB=',INTORB
       write(lupri,*) 'input DSM=',DSM
       write(lupri,*) 'output NUMLIN=',NUMLIN
      ENDIF

      RETURN
      END
C********************************************************
C* This are the subroutines that  reads  the primitives *
C* and the contracted coefficients from file.           *
C*********************************************************
C* This is the subroutine that searches through the file *
C* and finds the element in question.                    *
C*********************************************************
C/* Deck FIND_ELEMENT */
      SUBROUTINE FIND_ELEMENT(BASNAM,INTQ,LUBAS,IPRINT)
#include "implicit.h"
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basis file               *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*********************************************************
C
      CHARACTER*(*) BASNAM
      CHARACTER*88 STRING
      CHARACTER SIGN
C
      CALL FIND_BASFIL(BASNAM,LUBAS,IPRINT)
C
C Searching the file for the element.
C
 20   CONTINUE
         READ(LUBAS,'(A88)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)',IOSTAT=IOS) SIGN
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
            READ (STRING, 110,IOSTAT=IOS) SIGN, NUCEL
            IF (IOS.NE.0) THEN
               CALL QUIT('FIND_ELEMENT: Error in reading NUCEL')
            ENDIF
 110        FORMAT (BN, A1, I4)
            IF (INTQ .EQ. NUCEL) RETURN
         END IF
      GOTO 20
C
C Error messages 
C
 2000 CONTINUE
      INQUIRE (UNIT = LUBAS, NAME = STRING)
      IEND = INDEX(STRING(1:),' ') - 1
      WRITE (LUPRI,'(/2A/A,I5)')
     &   'FIND_ELEMENT: Error when reading from basis file ',BASNAM,
     &   ' IOSTAT =',IOERR
      CALL QUIT('I/O error in FIND_ELEMENT in HERBAS')
 200  CONTINUE
      WRITE (LUPRI,'(/I3,2A)') INTQ,
     &     ' is an unsupported element for basis ',BASNAM
      CALL QUIT('Unsupported element in HERBAS')
      END 
C/* Deck finnqd */
      SUBROUTINE FINNQD(BASNAM, INTQ,NBLOCK,ALPHA,CPRIMU,CPRIM,NUC,
     &                  NRC,SEG,KPRIM,KAOVEC,DSM,POLFUN,IPOLST,
     &                  UNCONT,LUBAS,IPRINT)
#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basis file               *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*********************************************************
C
C
      CHARACTER*(*) BASNAM
      CHARACTER*140 STRING
      CHARACTER SIGN
      LOGICAL   POLFUN, SEG, SEGEJ, BLANK, UNCONT
      DIMENSION ALPHA(KPRIM,KAOVEC), CPRIMU(KPRIM,KPRIM,KAOVEC),
     &          CPRIM(KPRIM,KPRIM,KAOVEC), NUC(KAOVEC), NRC(KAOVEC),
     &          SEG(KAOVEC)
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
C     Determine primitive basis
C
      ISTART = 5
      CALL FREFRM(BASNAM,ISTART,ISPRIM,DUMMY,'INT',IERR)
      IF (INTQ .GT. 4) THEN
         CALL FREFRM(BASNAM,ISTART,IPPRIM,DUMMY,'INT',IERR)
      END IF
      IF (UNCONT) THEN
C     ... force uncontracted irrespective of .mol file
         ISCONT = ISPRIM
         IPCONT = IPPRIM
      ELSE
         CALL FREFRM(BASNAM,ISTART,ISCONT,DUMMY,'INT',IERR)
         IF (INTQ .GT. 4) THEN
            CALL FREFRM(BASNAM,ISTART,IPCONT,DUMMY,'INT',IERR)
         END IF
      END IF
C
      CALL FIND_BASFIL('NQvD',LUBAS,IPRINT)
C
C Searching the file for the element.
C
 20   CONTINUE
         READ(LUBAS,'(A88)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)', IOSTAT = IOS) SIGN
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
            ISTART = 2
            CALL FREFRM(STRING,ISTART,NUCEL,DUMMY,'INT',IERR)
            CALL FREFRM(STRING,ISTART,IS,DUMMY,'INT',IERR)
            IF (NUCEL .EQ. INTQ .AND. IS .EQ. ISPRIM) THEN
               IF (NUCEL .GT. 4) THEN
                  CALL FREFRM(STRING,ISTART,IP,DUMMY,'INT',IERR)
                  IF (IP .NE. IPPRIM) GO TO 20
               END IF
               NBLOCK = NBLOCK + 1
               IF (NBLOCK .GT. KAOVEC) THEN
                  CALL QUIT('Too many AO blocks; increase MXAOVC')
               END IF
C
C Setting ALPHA, CPRIMU AND CPRIM to zero. 
C
               CALL DZERO(ALPHA(1, NBLOCK), KPRIM)
               CALL DZERO(CPRIMU(1,1, NBLOCK), KPRIM*KPRIM)
               CALL DZERO(CPRIM (1,1, NBLOCK), KPRIM*KPRIM)
               IF (NUCEL .GT. 4) THEN
                  CALL DZERO(ALPHA(1, NBLOCK + 1), KPRIM)
                  CALL DZERO(CPRIMU(1,1, NBLOCK + 1), KPRIM*KPRIM)
                  CALL DZERO(CPRIM (1,1, NBLOCK + 1), KPRIM*KPRIM)
               END IF
               IF (NUCEL .LE. 4) THEN
                  ICONT = 1
                  DO ISLOOP = 1, ISPRIM
                     READ (LUBAS,*,IOSTAT=IOS) 
     &                ALPHA(ISLOOP,NBLOCK), CONTS
                     IF (IOS.NE.0) THEN
                        CALL QUIT('FINNQD: Error reading ALPHA, CONTS')
                     ENDIF
                     IF (NUCEL .EQ. 1) ALPHA(ISLOOP,NBLOCK) =
     &                    ALPHA(ISLOOP,NBLOCK)*1.44D0
                     IF (ICONT .LT. ISCONT) THEN
                        CPRIMU(ISLOOP,ICONT,NBLOCK) = 1.0D0
                        ICONT = ICONT + 1
                     ELSE
                        CPRIMU(ISLOOP,ICONT,NBLOCK) = CONTS
                     END IF
                  END DO
               ELSE
                  ICONTS = 1
                  ICONTP = 1
                  DO IPLOOP = 1, IPPRIM
                     READ (LUBAS,*,IOSTAT=IOS) 
     &                     ALPHA(IPLOOP,NBLOCK), CONTS,
     &                     ALPHA(IPLOOP,NBLOCK + 1), CONTP
                     IF (IOS.NE.0) THEN
                       CALL QUIT('FINNQD Error reading ALPHA,CONTS,...')
                     ENDIF
                     IF (ICONTP .LT. IPCONT) THEN
                        CPRIMU(IPLOOP,ICONTP,NBLOCK + 1) = 1.0D0
                        ICONTP = ICONTP + 1
                     ELSE
                        CPRIMU(IPLOOP,ICONTP,NBLOCK + 1) = CONTP
                     END IF
                     IF (ICONTS .LT. ISCONT) THEN
                        CPRIMU(IPLOOP,ICONTS,NBLOCK) = 1.0D0
                        ICONTS = ICONTS + 1
                     ELSE
                        CPRIMU(IPLOOP,ICONTS,NBLOCK) = CONTS
                     END IF
                  END DO
                  DO ISLOOP = IPPRIM+1, ISPRIM
                     READ (LUBAS,*,IOSTAT=IOS)
     &                 ALPHA(ISLOOP,NBLOCK), CONTS
                     IF (IOS.NE.0) THEN
                        CALL QUIT('FINNQD Error reading ALPHA,CONTS')
                     ENDIF
                     IF (ICONTS .LT. ISCONT) THEN
                        CPRIMU(ISLOOP,ICONTS,NBLOCK) = 1.0D0
                        ICONTS = ICONTS + 1
                     ELSE
                        CPRIMU(ISLOOP,ICONTS,NBLOCK) = CONTS
                     END IF
                  END DO
               END IF
               IF (POLFUN) THEN
                  ISTART = IPOLST
 33               CONTINUE
                  CALL CHBLANK(BLANK,BASNAM(ISTART:))
                  IF (.NOT. BLANK) THEN
                     CALL FREFRM(BASNAM,ISTART,IQUANT,DUMMY,'INT',IERR)
                     IF (IERR .EQ. 0) THEN
                     CALL FREFRM(BASNAM,ISTART,IDUMMY,EXPON,'REA',IERR2)
                     IF (IERR2 .EQ. 0) THEN
                     IF (IQUANT .EQ. NBLOCK) THEN
                        ISPRIM = ISPRIM + 1
                        ISCONT = ISCONT + 1
                        ALPHA(ISPRIM,NBLOCK) = EXPON
                        CPRIMU(ISPRIM,ISCONT,NBLOCK) = 1.0D0
                     ELSE IF (IQUANT .EQ. (NBLOCK + 1) 
     &                       .AND. NUCEL .GT. 4) THEN
                        IPPRIM = IPPRIM + 1
                        IPCONT = IPCONT + 1
                        ALPHA(IPPRIM,NBLOCK + 1) = EXPON
                        CPRIMU(IPPRIM,IPCONT,NBLOCK + 1) = 1.0D0
                     END IF
                     END IF
                     END IF
                     IF ((IERR .EQ. 0) .AND. (IERR2 .EQ. 0)) GOTO 33
                  END IF
               END IF
               IF (ISPRIM.GT.KPRIM ) GOTO 5010
               IF (ISCONT.GT.KAOVEC) GOTO 5020          
C
C Checking if the matrix, with the contraction coefficients, is segmented 
C
               SEGEJ  = .TRUE.
               CALL SEGORB(SEGEJ, ISPRIM, ISCONT, CPRIMU(1,1,NBLOCK), 
     +                     KPRIM, DSM)
               NUC(NBLOCK) = ISPRIM
               NRC(NBLOCK) = ISCONT
               SEG(NBLOCK) = SEGEJ
C
C Normalizing CPRIMU for (angular quantum number + 1).
C
               CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +                     ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                     CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
               IF (NUCEL .GT. 4) THEN
C
C Checking if the matrix, with the contraction coefficients, is segmented 
C
                  NBLOCK = NBLOCK + 1
                  SEGEJ = .TRUE.
                  CALL SEGORB(SEGEJ, IPPRIM, IPCONT,
     +                        CPRIMU(1,1,NBLOCK),KPRIM, DSM)
                  NUC(NBLOCK) = IPPRIM
                  NRC(NBLOCK) = IPCONT
                  SEG(NBLOCK) = SEGEJ
C
C Normalizing CPRIMU for (angular quantum number + 1).
C
                  CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +                        ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                        CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
               END IF
               RETURN
            ELSE
               GOTO 20
            END IF
         ELSE 
            GOTO 20
         END IF
C
C Error messages 
C
 2000 CONTINUE
      WRITE (LUPRI,'(/2A/A,I5)')
     &   'FINNQD: Error when reading from basis file ',
     &   BASNAM,
     &   ' IOSTAT =',IOERR
      CALL QUIT('I/O error in FINNQD in HERBAS')
 200  CONTINUE
      WRITE (LUPRI,'(/I3,2A)') INTQ,
     &     ' is an unsupported element for basis ',BASNAM
      CALL QUIT('Unsupported element in specified basis set file')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/2A/A,I6)')
     &  ' >>> FINNQD error, number of primitives per block  ',ISPRIM,
     &  '                   in basis ',BASNAM,
     &  '                   current maximum number          ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/2A/A,I6)')
     &  ' >>> FINNQD error, number of contracted functions  ',ISCONT,
     &  '                   in basis ',BASNAM,
     &  '                   current maximum number          ',KAOVEC
        CALL QUIT('Too many contracted functions')
      END 
C/* Deck finref */
      SUBROUTINE FINREF(BASREF, LUBAS, BASNAM)
C*********************************************************
C* This is the subroutine that searches through the file *
C* and finds reference info if available.                *
C* Based on FINPOS, LV, 14-5-2003                         *
C*********************************************************
#include "implicit.h"
#include "priunit.h"
      CHARACTER*80 LINE,BASREF(10), BASNAM
C*********************************************************
C* Variable declarations:                                *
C*     SIGN - The first sign in a sentence, to find out  *
C*            what kind of sentence it is.               *
C*     STRING - A character string to help bulletproofing*
C*              the subroutine.                          *
C*********************************************************
      CHARACTER SIGN
      CHARACTER*88 STRING
      LOGICAL BLANK
C
      REWIND (LUBAS)
C
C Initialising IREF (number of documentation lines), set to one
C once we find a line containing the string REFERENCE
C
      IREF = 0
C
C     The loops are written f77 style since f90 may still not be
C     available everywhere, better to use while loops here..
C
 20   CONTINUE
      READ (LUBAS, '(A88)', IOSTAT = IOERR, ERR = 1000, END = 200)
     +     STRING
      IF (IOERR .NE. 0) GO TO 200
      READ (STRING, '(A1)',IOSTAT=IOS) SIGN
C
      IF (SIGN .EQ. '$') THEN
         IF (IREF.EQ.0) THEN
C           Check for start of documentation 
            DO I = 2, 80
               IF (STRING(I:I+8).EQ.'REFERENCE') IREF = 1
            ENDDO
         ELSE
C           We have documentation, omit blank lines and store in BASREF
            CALL CHBLANK(BLANK, STRING(2:81))
            IF (.NOT.BLANK) THEN
               BASREF(IREF) = STRING(2:81)
               IREF = IREF + 1
            ENDIF
         END IF
C
C        Maximum number of documentation lines is 10
C
         IF (IREF.GT.10) RETURN
      ELSE 
C        We return if we found a REFERENCE line but no more $ lines.
         IF (IREF.GT.0) RETURN
      END IF
C
      GOTO 20
C
C No documentation when end of file is found
C
 200  CONTINUE
 1000 CONTINUE
!radovan: this error printout below is quite annoying
!         when you use own libraries to organize basis sets
!         then you may get lots (no of processors)
!         of these messages
CMI: well, rather use warning message
CMI   WRITE (LUPRI,'(//A,I5/2A/)')
CMI  &     'Error in reading basis set file, your basis '//
CMI  &     'has no documentation. IOSTAT code:',IOERR,
CMI  &     'Basis set: ',BASNAM
      WRITE (LUPRI,'(/A/A,A)')
     &     '  Info about the basis set file: your basis '//
     &     'has no documentation.',
     &     '  Basis set: ',BASNAM

Chj   CALL QUIT('Incomplete basis set file in HERBAS')
      RETURN
      END
C             
C
C*********************************************************
C* This is the subroutine that finds out whether it is a *
C* new block with primitives and contraction             *
C* coefficients, or if it is a new element. It returns   *
C* three variables:                                      *
C*          NEWEL -Logical variable that gives whether it*
C*                 is a new element or not.              *
C*          INEXP, INTORB - Number of primitives and     *
C*                          columns of contraction       *
C*                          coefficients                 *
C*********************************************************
C
C/* Deck finpos */
      SUBROUTINE FINPOS(NEWEL, INTEXP, INTORB, INTISG, LUBAS, IPRINT)
#include "implicit.h"
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     SIGN - The first sign in a sentence, to find out  *
C*            what kind of sentence it is.               *
C*     STRING - A character string to help bulletproofing*
C*              the subroutine.                          *
C*     BLANK - A logical variable to help bulletproofing *
C*             the subroutine.
C*********************************************************
      LOGICAL NEWEL
      CHARACTER SIGN
      CHARACTER*88 STRING
      LOGICAL BLANK
C
C Initializing NEWEL, because one is innocent until proven guilty
C
      NEWEL = .FALSE.
C
 20   CONTINUE
      READ (LUBAS, '(A88)', IOSTAT = IOERR, ERR = 1000, END = 200)
     +     STRING
      READ (STRING, '(A1)',IOSTAT=IOS) SIGN
C
      IF (SIGN .EQ. ' ') THEN
C
C Calling on a subroutine that checks whether the line is a blank one or not.
C If STRING is a blank line, then BLANK will return as .TRUE.
         CALL CHBLANK(BLANK, STRING)
         IF (BLANK) THEN
            GOTO 20
         ELSE
            READ (STRING, '(BN,3I5)',IOSTAT=IOS) INTEXP, INTORB, INTISG
            IF (IOS.NE.0) GO TO 2000
            IF (IPRINT.GE.10) write(lupri,'(/2X,A,3I10)')
     &        'Output from FINPOS... read INTEXP, INTORB, INTISG',
     &          INTEXP, INTORB, INTISG
            RETURN
         END IF
      ELSE IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
         NEWEL = .TRUE.
         IF (IPRINT.GE.10) write(lupri,'(/2X,A)')
     &      'Output from FINPOS: New element'
         RETURN
      ELSE 
         GOTO 20
      END IF
C
C No more orbitals when end of file, so the orbitals for the last element 
C are read
C
 200  CONTINUE
      NEWEL = .TRUE.
C
C NEWEL must be .TRUE. to break the if-loop.
C
      IF (IPRINT.GE.10) write(lupri,'(/2X,A)')
     &   'Output from FINPOS: Reached end of file'

      RETURN
C
 1000 CONTINUE
      WRITE (LUPRI,'(/A,I5)') 'Error in reading basis set file, '//
     &     'your basis is not complete. IOSTAT code:',IOERR
      CALL QUIT('Incomplete basis set file in FINPOS')
 2000 CONTINUE
      WRITE (LUPRI,'(/A,I5//A/A)')'Error in searching basis set file,'//
     &     ' your basis is not complete. IOSTAT code:',IOERR,
     &     'Last line read:',STRING
      CALL QUIT('Error searching basis set file in FINPOS')
      END
C/* Deck FIND_BASFIL */
      SUBROUTINE FIND_BASFIL(BASNAM,LUBAS,IPRINT)
C
C     Joern Thyssen and Hans Joergen Aa. Jensen.
C     Last revision Feb. 2005
C
#include <implicit.h>
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basis file               *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*     EXST   - Logical variable that helps inquire if   *
C*              there exists a file with that name.      *
C*********************************************************
C
      CHARACTER*(*) BASNAM
      CHARACTER*200 STRING ! long because will contain full directory path
      LOGICAL EXST
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
C     The environment variable is a : separated string, e.g.:
C     .:/programs/Dirac/basis:/programs/Dalton/basis
C
C     We loop over these directories in search for the basis set and
C     element given.
C
      EXST = .FALSE.
      ISTART = 1
C
 10   CONTINUE
C
C     Find next colon in BASDIR
C
      ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:), ':' )
C
C     In there are no colons, we are at the last directory in BASDIR
C
      IF ( ICOLON .EQ. ISTART - 1 ) THEN
C
C        Find end of string
C        (done this way, because user may have specified
C        a directory name with one or more blanks in it,
C        as /home/xxx/My jobs/My basdir/ )
C
         DO I = LBASDIR,ISTART,-1
            IF (BASDIR(I:I) .NE. ' ') THEN
               ICOLON = I+1
C              ... the address of the "virtual" colon
C                  (because IEND = ICOLON - 1 below)
               GO TO 20
            END IF
         END DO
C
C        If end of string = start of string there are no more
C        directories, so we gracefully die...
C
   20    IF ( ICOLON .LE. ISTART ) GOTO 300
      END IF
      IEND = ICOLON - 1
      IF (IEND .LT. ISTART) THEN
C        User has started with ':',
C        or has two consecutive colons as '::',
C        we skip to next item.
         ISTART = ICOLON + 1
         GO TO 10
      END IF
      STRING = BASDIR(ISTART:IEND)
      LSTR   = IEND + 1 - ISTART
      ISTART = ICOLON + 1
C
C     Add trailing slash if not present.
C
      IF ( STRING(LSTR:LSTR) .NE. '/' ) THEN
         LSTR = LSTR + 1
         STRING(LSTR:LSTR) = '/'
      END IF
      LBAS = INDEX(BASNAM,' ') - 1
C     ... cannot use LNBLNK because some basis sets have additional
C         information in BASNAM, e.g. "ano-2 4 3 2'
      IF (LBAS .LE. 0) LBAS = LEN(BASNAM)
C     ... take care of basis set names with no trailing blanks
      STRING = STRING(1:LSTR)//BASNAM
      LSTR = LSTR + LBAS
      IF (IPRINT .GT. 1) WRITE(LUPRI,'(3A)')
     &   '  Trying file: "',STRING(1:LSTR),'"'
C
C     Inquire if "/path/basis-set" exists.
C
      INQUIRE (FILE = STRING(1:LSTR), EXIST = EXST)
      IF (EXST) THEN
C        We print which basis set is used, user may have
C        specified several basis set directories, and
C        without this output there is not documentation
C        for which is used.
C        This text is continuation of text in BASLIB.
C        (By default suppressed for HUCKEL basis by setting
C         IPRINT = IPRBAS - 1 for HUCKEL) /hjaaj
         IF (IPRINT .GE. 1) WRITE (LUPRI,'(3A)')
     &     '     "',STRING(1:LSTR),'"'
         CALL GPOPEN(LUBAS,STRING(1:LSTR),'OLD',' ','FORMATTED',IDUMMY,
     &               .TRUE.)
         RETURN
      ELSE
         GO TO 10
      END IF
C
  300 WRITE (LUPRI,'(/3A//A/A)')
     +      'Basis "',BASNAM,'" does not exist.',
     +      'The following directories have been searched :',BASDIR
      IF (INDEX(BASNAM,'Huckel') .GT. 0) THEN
         WRITE (LUPRI,'(/A)')
     +   'NB! This basis must be available to use Huckel guess!'
      END IF
      CALL QUIT('Non-existing basis set in HERBAS')
      END
C/* Deck ECP_LCORE */
      SUBROUTINE ECP_LCORE(IQCORE,IQM,ILOFF)
C
C     Nov. 2003, Johan, Patrick and Hans Joergen.
C
C     Purpose: calculate ILOFF = number of core orbitals described by ECP
C     for l_quantum_number = IQM - 1
C     when IQCORE electrons are described with the ECP.
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION ILCORE(4)
C
      IF (IQCORE .GT. 0 .AND. IQM .LE. 4) THEN
C        only s, p. d, f (never g, h, ... core orbitals)
         IF (MOD(IQCORE,2) .NE. 0) THEN
            WRITE(LUPRI,'(/A,I5)')
     &      'ECP_LCORE ECP: odd number of core electrons',IQCORE
            CALL QUIT('ECP_LCORE ECP: odd number of core electrons')
         END IF
         DO I = 1, 4
            ILCORE(I) = 0
         END DO
         IF (IQCORE .EQ. 2) THEN
C           1s ECP
            ILCORE(1) = 1
         ELSE IF (IQCORE .EQ. 10) THEN
C           1s,2s,2p ECP
            ILCORE(1) = 2
            ILCORE(2) = 1
         ELSE IF (IQCORE .EQ. 18) THEN
C           1s,2s,2p,3s,3p ECP
            ILCORE(1) = 3
            ILCORE(2) = 2
         ELSE IF (IQCORE .EQ. 28) THEN
C           1s,2s,2p,3s,3p,3d ECP
            ILCORE(1) = 3
            ILCORE(2) = 2
            ILCORE(3) = 1
         ELSE IF (IQCORE .EQ. 36) THEN
C           1s,2s,2p,3s,3p,3d,4s,4p ECP
            ILCORE(1) = 4
            ILCORE(2) = 3
            ILCORE(3) = 1
         ELSE IF (IQCORE .EQ. 46) THEN
C           1s,2s,2p,3s,3p,3d,4s,4p,4d ECP
            ILCORE(1) = 4
            ILCORE(2) = 3
            ILCORE(3) = 2
         ELSE IF (IQCORE .EQ. 60) THEN
C           1s,2s,2p,3s,3p,3d,4s,4p,4d,4f ECP
            ILCORE(1) = 4
            ILCORE(2) = 3
            ILCORE(3) = 2
            ILCORE(4) = 1
         ELSE IF (IQCORE .EQ. 78) THEN
C           1s,2s,2p,3s,3p,3d,4s,4p,4d,4f,5s,5p,5d ECP
            ILCORE(1) = 5
            ILCORE(2) = 4
            ILCORE(3) = 3
            ILCORE(4) = 1
         ELSE
            WRITE(LUPRI,'(/A,I5)') 'REAANO FATAL ERROR: '//
     &         'ECP Huckel not implemented for this core charge:',IQCORE
            CALL QUIT('REAANO ECP not implemented for this core charge')
         END IF
         ILOFF = ILCORE(IQM)
      ELSE
         ILOFF = 0
C        not ECP or g, h, ...
C        i.e. no core orbitals to exclude
      END IF
      RETURN
      END
C  -- end of herbas.F --
