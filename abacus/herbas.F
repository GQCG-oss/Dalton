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
C*************************************************************
C* This is the procedures that controls the reading from     *
C* file, and calls all the other subroutines.                *
C* It returns a bunch of variables. It returns:              *
C*     IQM - (Highest anglar quantumnumber) + 1              *
C*     JCO - Redundant variable that gives how many blocks   *
C*           the primitives and contraction coeffecients are *
C*           given in (will be set to one) (array)           *
C*     NUC - Number of primitives for a given                *
C*           (quantumnuber + 1)                              *
C*     NRC - Number of coloumns with contraction-coeffecients*
C*           for a given (quantumnumber + 1).                *
C*     SEG - Gives whether the contraction-coeffecients are  *
C*           segmented (for a given (quantumnumber + 1),     *
C*           (array).                                        *
C*     ALPHA - Array with the primitives (for a given        *
C*             (quantumnumber + 1)).                         *
C*     CPRIMU - Matrix (3D) with the NOT normalized          *
C*              contraction coeffecients (for a given        *
C*              (quantumnumber + 1)).                        *
C*     CPRIM - Like CPRIMU, but this time the contraction    *
C*             -coeffecients are normalized.                 *
C*     NBLOCK - For some reason is this defined as the same  *
C*              as IQM                                       *
C*     KAOVEC, KPRIM, Q, DSM are max (angular quantumnumber  *
C*             + 1), max number of primitives, nuclear       *
C*             charge, highest number to be accepted as zero *
C*     QEFF is eff. nuclear charge which is .le. Q,          *
C*             Q-QEFF will be described with ECP              *
C*************************************************************
C
      SUBROUTINE BASLIB (IQM, JCO, NUC, NRC, SEG, ALPHA, CPRIM,
     &                   CPRIMU, NBLOCK, KAOVEC, KPRIM, Q, QEFF, DSM)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      LOGICAL SEG
      DIMENSION JCO(MXAQN), NUC(KAOVEC), NRC(KAOVEC), SEG(KAOVEC),
     +          ALPHA(KPRIM, KAOVEC), CPRIMU(KPRIM, KPRIM, KAOVEC), 
     +          CPRIM(KPRIM,KPRIM, KAOVEC)
#include <cbirea.h>
C*************************************************************
C* Variable declarations:                                    *
C*     BASNAM - A character variable that contains the name  *
C*              of the file.                                 *
C*     NEWEL - Logical variable, that is returned from       *
C*             the subroutine finpos, gives whether the next *
C*             line in basis-file is a new element or not.   *
C*     SEGEJ - This gives whether the contraction            *
C*             -coeffecients are segmented or not (for       *
C*             (angular quantum-number + 1))                 *
C* External procedures:                                      *
C*     SEGORB - A procedure that checks if the contraction   *
C*              -coeffecients (for a given (angular quantum  *
C*              number + 1)) are segmented or not.           *
C*     NRMORB - Normalizing a matrix with contraction        *
C*              -coeffecients                                *
C* cbirea.h:
C*     UNCONT - force uncontracted (hjaaj aug99+aug00)
C*************************************************************
C
      CHARACTER BASTMP*80, BASSAV*80
      LOGICAL NEWEL, SEGEJ, ANO, SADLEJ, POLFUN, FOUND, BLANK, NQVD
C
      BASTMP = '                                        '//
     &         '                                        '
      BASSAV = BASNAM
      INTQ   = NINT(Q)
      LUBAS  = -1
      NBLOCK = 0
      IAUG   = 0
      ANO    = .FALSE.
      SADLEJ = .FALSE.
      NQVD   = .FALSE.
      IF (BASNAM(1:6) .EQ. 'HUCKEL') THEN
C jan 2000 hjaaj:
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
C
         ANO = .TRUE.
         BASTMP = 'ano-4'
C        NOTE that ECP Huckel start guess is implemented
C        in REAANO through the Q-QEFF parameter. 
         ISTART = 6
         IF (INTQ .LE. 2) THEN
C           ... H and He
            BASNAM = 'ano-4   1 0 0 0'
         ELSE IF (INTQ .LE. 4) THEN
C           ... Li, Be
            BASNAM = 'ano-4   2 0 0 0'
         ELSE IF (INTQ .LE. 10) THEN
C           ... B, C, N, O, F, Ne
            BASNAM = 'ano-4   2 1 0 0'
         ELSE IF (INTQ .LE. 12) THEN
C           ... Na, Mg
            BASNAM = 'ano-4   3 1 0 0'
         ELSE IF (INTQ .LE. 18) THEN
C           ... Al, Si, P, S, Cl, Ar
            BASNAM = 'ano-4   3 2 0 0'
         ELSE IF (INTQ .LE. 20) THEN
C           ...  K, Ca
            BASNAM = 'ano-4   4 2 0 0'
         ELSE IF (INTQ .LE. 30) THEN
C           ...  Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn
            BASNAM = 'ano-4   4 2 1 0'
         ELSE IF (INTQ .LE. 36) THEN
C           ...  Ga, Ge, As, Se, Br, Kr
            BASNAM = 'ano-4   4 3 1 0'
         ELSE
            CALL QUIT('BASLIB error: Huckel not defined for Z.gt.36')
         END IF
C djw Feb 2005
C notation changed from 'd' to 'd-' for consistency with EMSL
      ELSE
         IF (BASNAM(3:10) .EQ. 'aug-cc-p') THEN
            BASTMP(1:79) = BASNAM(3:80)
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
            BASTMP(1:75) = BASNAM(5:80)
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
            BASTMP(1:75) = BASNAM(5:80)                      
            IAUG = 1 
         ELSE IF (BASNAM(2:8) .EQ. 'aug-ecp') THEN           
            BASTMP(1:74) = BASNAM(6:80)
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
C    
C Finds the right element in the file. 
C
      IF (NQVD) THEN
         CALL FINNQD(BASNAM,INTQ,NBLOCK,ALPHA,CPRIMU,CPRIM,NUC,NRC,SEG,
     &               KPRIM,KAOVEC,DSM,POLFUN,IPOLST,UNCONT,LUBAS)
      ELSE
         CALL FINDEL(BASTMP,INTQ,LUBAS)
C
C Finds the number of primitives and columns of contraction coeffecients
C
         CALL FINPOS ( NEWEL, INTEXP, INTORB, LUBAS)
C   
 10      CONTINUE
         IF ( .NOT. NEWEL) THEN
            NBLOCK = NBLOCK + 1
            IF (NBLOCK .GT. KAOVEC) THEN
               CALL QUIT('Too many AO blocks'//
     &            ' for available memory; increase WRKMEM')
            END IF
C   
C Setting ALPHA, CPRIMU AND CPRIM to zero. 
C
            CALL DZERO(ALPHA(1, NBLOCK), KPRIM)
            CALL DZERO(CPRIMU(1,1, NBLOCK), KPRIM*KPRIM)
            CALL DZERO(CPRIM (1,1, NBLOCK), KPRIM*KPRIM)
C
C Reading the primitives and contraction-coeffecients from file in READNU
C
            IF (ANO .OR. SADLEJ) THEN
               IQCORE = NINT(Q-QEFF)
               CALL REAANO(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU,ALPHA,
     +                     KPRIM,ANO,SADLEJ,BASNAM,ISTART,POLFUN,IPOLST,
     +                     UNCONT,LUBAS,IQCORE)
               IF (INTORB .EQ. 0) THEN
                  NBLOCK = NBLOCK - 1
                  CALL FINPOS(NEWEL,INTEXP,INTORB,LUBAS)
                  GOTO 10
               END IF
            ELSE
Chj-aug99: if (UNCONT) force uncontracted
               IF (UNCONT) INTORB = -INTORB
               CALL READNU(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU,ALPHA,
     +                     KPRIM,IAUG,DSM,BASNAM,POLFUN,IPOLST,LUBAS)
            END IF
            IF (INTEXP.GT.KPRIM) GOTO 5010
            IF (INTORB.GT.KPRIM) GOTO 5020          
C   
C Checking if the matrix, with the contraction-coeffecients, is segmented 
C
            SEGEJ  = .TRUE.
            CALL SEGORB(SEGEJ,INTEXP,INTORB,CPRIMU(1,1,NBLOCK),KPRIM,
     +                  DSM)
            NUC(NBLOCK) = INTEXP
            NRC(NBLOCK) = INTORB
            SEG(NBLOCK) = SEGEJ
C
C     For Douglas-Kroll we need to keep exponents properly sorted
C
            CALL PRIORD(ALPHA(1,NBLOCK),CPRIMU(1,1,NBLOCK),NUC(NBLOCK),
     &                  NRC(NBLOCK),SEG(NBLOCK),KPRIM,DSM)
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
             CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +                   ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                   CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
C
C Define the type of the next relevant line in file
C
            CALL FINPOS ( NEWEL, INTEXP, INTORB, LUBAS )
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
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
            CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +                  ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                  CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
            GOTO 32
         ELSE
            NBLOCK = NBLOCK - 1
         END IF
      END IF
      IQM = NBLOCK
C
C Setting JCO to one.
C
      DO I = 1, IQM
         JCO(I) = 1
      END DO
C
C Closing file
C
      BASNAM = BASSAV
      CALL GPCLOSE(LUBAS,'KEEP')
      RETURN
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error,'//
     *  'number of primitives per block      ',INTEXP,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error, number of contracted functions      ',
     *  INTORB,
     *  '                   current maximum number              ',KPRIM
C
      END
C
C*************************************************************
C* This is the subroutine that takes care of the reading of  *
C* the primitives og the contraction coeffecients from the   *
C* file (for a given format). The variables that are         *
C* transferred out are:                                      *
C*     CPRIMU - The matrix that the contraction coeffecients *
C*              are put in (not normalized).                 *
C*     ALPHA - The matrix where the primitives are put in    *
C*     (Both of the variables are given for a given (quantum-*
C*     number + 1))                                          *
C*************************************************************
C
C/* Deck readnu */
      SUBROUTINE READNU(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU, ALPHA, 
     +                  KPRIM,IAUG,DSM,BASNAM,POLFUN,IPOLST,LUBAS)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
C Reading the exponents
C
      DIMENSION CPRIMU( KPRIM, KPRIM, KAOVEC), ALPHA(KPRIM, KAOVEC)
C*************************************************************
C* Variable declarations:                                    *
C*     SIGN - Character that gives what kind of sentence it  *
C*            is.                                            *
C*     STRING - A character string that helps bulletproofing *
C*             the subroutine                                *
C*     BLANK - A logical variable that gives whether the line*
C*             is blank or not.                              *
C*************************************************************
      CHARACTER SIGN
      CHARACTER STRING*88, BASNAM*80
      LOGICAL BLANK, POLFUN
C
C CONEOR is the maximum number of coloumns of contraction coeffecients on 
C the 1. line in file. If INTORB is .GT. CONEOR then the contraction
C coeffecients will continue on the next line. CONORB is the max. number
C of contraction coeffecients on the next lines.
C
      INTEGER CONEOR, CONORB
      PARAMETER(CONEOR = 6, CONORB = 7)
C
C Initializing a counting variable that must be 0 before loop
C
      J = 0
C
 20   CONTINUE
      IF ( J .LT. INTEXP) THEN
C
C Reading the primitive and contracted coeffecients
C
         READ(LUBAS, '(A88)', IOSTAT = IOERR, ERR = 2000) STRING
         READ (STRING, '(A1)') SIGN
         IF (SIGN .EQ. ' ') THEN
C
C Calling a subroutine that checks whether the line is a blank one. If
C STRING is a blank line BLANK will return as .TRUE.
C
            CALL CHBLANK(BLANK, STRING)
            IF (.NOT. BLANK) THEN
C
C We have found a line with a primitive and coeffecients             
C
               J = J + 1
C
C CONLIN returns the number of lines we have the contraction
C coeffecients written on.
C
               JNTORB = ABS(INTORB)
               CALL CONLIN(JNTORB, NUMLIN, DSM)
               IF (INTORB .LE. 0) THEN
C                 uncontracted basis set or
C                 forced uncontracted with .UNCONT
                  READ (STRING, '(F16.9)') ALPHA(J,NBLOCK)
C                 skip contiuation lines for contraction coeff.
                  DO I = 2,NUMLIN
                     READ (LUBAS,*)
                  END DO
                  CPRIMU(J,J,NBLOCK) = 1.0D0
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
C     Reading the first line with exponents and contraction- coeffecients
C
                  READ (STRING, '(F16.9, 6F12.9)')
     +                 ALPHA(J, NBLOCK), (CPRIMU(J,I,NBLOCK),
     +                 I = 1, KNTORB)
C   
C     If there are more lines with contraction-coeffecients, 
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
C     contraction-coeffecients.
C
                     READ(LUBAS, '(A88)', IOSTAT = IOERR, ERR = 2000)
     +                    STRING
                     READ (STRING, '(A1)') SIGN
                     IF (SIGN .EQ. ' ') THEN
                        CALL CHBLANK(BLANK, STRING)
                        IF (.NOT. BLANK) THEN
C   
C     We now have a line with contraction-coeffecients.
C   
                           READ (STRING,'(F16.9,6F12.9)') 
     &                          (CPRIMU(J,M,NBLOCK),
     +                          M = CONEOR + (I-2)*CONORB +1, KNTORB)
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
         ELSE
C
C Found a line with nothing sensible in it     
C
            GOTO 20
         END IF
      ELSE
C
C        For uncontracted basis sets, update INTORB.
C
         IF ( INTORB .LE. 0 ) INTORB = INTEXP
C
C   
C Read all of the primitives, and the coeffecients
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
                  END IF
               END IF
               IF ((IERR .EQ. 0) .AND. (IERR2 .EQ. 0)) GOTO 33
            END IF
         END IF
         RETURN
      END IF
C
 2000 CONTINUE
      WRITE (LUPRI,'(/A,I5)') 'Error in reading from basis file.'//
     &     ' IOSTAT =',IOERR
      END
C/* Deck reaano */
      SUBROUTINE REAANO(INTEXP, INTORB, NBLOCK, KAOVEC, CPRIMU, ALPHA, 
     +                  KPRIM,ANO,SADLEJ,BASNAM,ISTART,POLFUN,IPOLST,
     +                  UNCONT,LUBAS,IQCORE)
C
C     IQCORE is charge of core electrons (.eq. zero if not ECP)
C     NBLOCK = l_quatum_number + 1 (1 for s, 2 for p, etc.)
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
C
C Reading the exponents. Format from MOLCAS files
C
      LOGICAL ANO, SADLEJ, POLFUN, BLANK, UNCONT
      DIMENSION CPRIMU( KPRIM, KPRIM, KAOVEC), ALPHA(KPRIM, KAOVEC)
      CHARACTER BASNAM*80
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
      CALL ECPLCORE(IQCORE,NBLOCK,ILOFF)
C
C Reading the primitive and contracted coeffecients
C
      MAXEXP = 8
      NUMEXP = (INTEXP - 1)/MAXEXP
      J = 1
      IF (INTORB .GT. 0) THEN
         DO I = 0, NUMEXP
            IF (I .EQ. NUMEXP) THEN
               IF (INTEXP .EQ. MAXEXP) THEN
                  READ (LUBAS,*) (ALPHA(K,NBLOCK), K = J, 
     &                                              J - 1 + MAXEXP)
               ELSE
                  READ (LUBAS,*) (ALPHA(K,NBLOCK), K = J, 
     &                                    J - 1 + MOD(INTEXP,MAXEXP))
               END IF
            ELSE
               READ (LUBAS,*) (ALPHA(K,NBLOCK), K = J, J + MAXEXP - 1)
            END IF
            J = J + MAXEXP
         END DO
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
         DO I = 0, NUMEXP
            READ (LUBAS,*)
         END DO
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
C
C*************************************************************
C* This is a subroutine that checks whether a line is a blank*
C* line. It returns a logical variable if it is a blank one. *
C* The in/out-going variables are:                          *
C*     STRING - A character string (the string checked by the*
C*              subroutine).                                 *
C*     BLANK - The logical variable. It will be .TRUE. if the*
C*             line is a blank one.                          *
C*************************************************************
C/* Deck chblank */
      SUBROUTINE CHBLANK(BLANK,STRING)
C
      CHARACTER*(*) STRING
      LOGICAL BLANK
C
      BLANK = .TRUE.
C
      DO 100 J = 1, LEN(STRING)
         BLANK = BLANK .AND. (STRING(J:J) .EQ. ' ')
 100  CONTINUE
      RETURN
      END
C
C*************************************************************
C* This is a subroutine that finds out how many lines the    *
C* contraction- coeffecients are written on.  There are 6    *
C* cc (contraction-coeffecients) on the 1. line, and then    *
C* there are seven on the rest of the lines. In out going    *
C* variables are :                                           *
C*      INTORB - Total number of coloumns of contraction     *
C*               -coeffecients.                              *
C*      NUMLIN - Returns the number of lines the contaction  *
C*               -coeffecients are written on.               *
C*************************************************************
C/* Deck conlin */
C
      SUBROUTINE CONLIN(INTORB, NUMLIN, DSM)
#include <implicit.h>
C    
C This is the parameters for the number of contraction-coeffecients
C on 1. and 2. line.
C
      INTEGER CONEOR, CONORB
      PARAMETER(CONEOR = 6, CONORB = 7)
C
C The intrisic functions DBLE makes a souble precision reak number of an integer.
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
      END
C********************************************************
C* This are the subroutines that  reads  the primitives *
C* and the contracted coeffecients from file.           *
C*********************************************************    
C* This is the subroutine that searches through the file *
C* and finds the element in question.                    *
C*********************************************************
C/* Deck findel */
      SUBROUTINE FINDEL(BASNAM,INTQ,LUBAS)
#include <implicit.h>
#include <priunit.h>
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basis file               *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*********************************************************
C
      CHARACTER*(*) BASNAM
      CHARACTER*140 STRING
      CHARACTER SIGN
C
      CALL  FIND_BASFIL(BASNAM,LUBAS)
C
C Searching the file for the element.
C
 20   CONTINUE
         READ(LUBAS,'(A88)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)') SIGN
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
            READ (STRING, 110) SIGN, NUCEL
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
     &   'FINDEL: Error when reading from basis file ',STRING(1:IEND),
     &   ' IOSTAT =',IOERR
      CALL QUIT('I/O error in FINDEL in HERBAS')
 200  CONTINUE
      WRITE (LUPRI,'(/I3,2A)') INTQ,
     &     ' is an unsupported element for basis ',BASNAM
      CALL QUIT('Unsupported element in HERBAS')
      END 
C/* Deck finnqd */
      SUBROUTINE FINNQD(BASNAM, INTQ,NBLOCK,ALPHA,CPRIMU,CPRIM,NUC,
     &                  NRC,SEG,KPRIM,KAOVEC,DSM,POLFUN,IPOLST,
     &                  UNCONT,LUBAS)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basisfile                *
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
#include <gnrinf.h>
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
      CALL  FIND_BASFIL('NQvD',LUBAS)
C
C Searching the file for the element.
C
 20   CONTINUE
         READ(LUBAS,'(A88)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)') SIGN
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
                  CALL QUIT('Too many AO blocks'//
     &               ' for available memory; increase WRKMEM')
            END IF
C
C Setting ALPHA, CPRIMU AND CPRIM to zero. 
C
               CALL DZERO(ALPHA(1, NBLOCK), KPRIM)
               CALL DZERO(CPRIMU(1,1, NBLOCK), KPRIM*KPRIM)
               CALL DZERO(CPRIM (1, 1, NBLOCK), KPRIM*KPRIM)
               IF (NUCEL .GT. 4) THEN
                  CALL DZERO(ALPHA(1, NBLOCK + 1), KPRIM)
                  CALL DZERO(CPRIMU(1,1, NBLOCK + 1), KPRIM*KPRIM)
                  CALL DZERO(CPRIM (1, 1, NBLOCK + 1), KPRIM*KPRIM)
               END IF
               IF (NUCEL .LE. 4) THEN
                  ICONT = 1
                  DO ISLOOP = 1, ISPRIM
                     READ (LUBAS,*) ALPHA(ISLOOP,NBLOCK), CONTS
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
                     READ (LUBAS,*) ALPHA(IPLOOP,NBLOCK), CONTS,
     &                           ALPHA(IPLOOP,NBLOCK + 1), CONTP
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
                     READ (LUBAS,*) ALPHA(ISLOOP,NBLOCK), CONTS
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
               IF (ISPRIM.GT.KPRIM) GOTO 5010
               IF (ISCONT.GT.KPRIM) GOTO 5020          
C
C Checking if the matrix, with the contraction-coeffecients, is segmented 
C
               SEGEJ  = .TRUE.
               CALL SEGORB(SEGEJ, ISPRIM, ISCONT, CPRIMU(1,1,NBLOCK), 
     +                     KPRIM, DSM)
               NUC(NBLOCK) = ISPRIM
               NRC(NBLOCK) = ISCONT
               SEG(NBLOCK) = SEGEJ
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
C
               CALL NRMORB(NBLOCK, NRC(NBLOCK), NUC(NBLOCK), 
     +                     ALPHA(1,NBLOCK),CPRIM(1, 1, NBLOCK), 
     +                     CPRIMU(1,1, NBLOCK),  KPRIM, NBLOCK)
               IF (NUCEL .GT. 4) THEN
C
C Checking if the matrix, with the contraction-coeffecients, is segmented 
C
                  NBLOCK = NBLOCK + 1
                  SEGEJ = .TRUE.
                  CALL SEGORB(SEGEJ, IPPRIM, IPCONT,
     +                        CPRIMU(1,1,NBLOCK),KPRIM, DSM)
                  NUC(NBLOCK) = IPPRIM
                  NRC(NBLOCK) = IPCONT
                  SEG(NBLOCK) = SEGEJ
C
C Normalizing CPRIMU for (angular quantumnumber + 1).
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
     &   BASDIR(1:LENBAS)//'NQvD',
     &   ' IOSTAT =',IOERR
      CALL QUIT('I/O error in FINNQD in HERBAS')
 200  CONTINUE
      WRITE (LUPRI,'(/I3,2A)') INTQ,
     &     ' is an unsupported element for basis ',BASNAM
      CALL QUIT('Unsupported element in HERBAS')
 5010 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error,'//
     *  'number of primitives per block      ',ISPRIM,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many primitives')
 5020 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> GTOINP error, number of contracted functions      ',
     *  ISCONT,
     *  '                   current maximum number              ',KPRIM
        CALL QUIT('Too many contracted functions')
      END 
C             
C
C*********************************************************	
C* This is the subroutine that finds out wheter it is a  *
C* new block with primitives and contraction-            *
C* coeffecients, or if it is a new element. It returns   *
C* three variables:                                      *
C*          NEWEL -Logical variable that gives whether it*
C*                 is a new element or not.              *
C*          INEXP, INTORB - Number of primitives and     *
C*                          coloumns of contraction      *
C*                          coeffecients                 *
C*********************************************************
C
C/* Deck finpos */
      SUBROUTINE FINPOS(NEWEL, INTEXP, INTORB,LUBAS)
#include <implicit.h>
#include <priunit.h>
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
C Initialising NEWEL, because one is innocent until proven guilty
C
      NEWEL = .FALSE.
C
 20   CONTINUE
      READ (LUBAS, '(A88)', IOSTAT = IOERR, ERR = 1000, END = 200)
     +     STRING
      READ (STRING, '(A1)') SIGN
C
      IF (SIGN .EQ. ' ') THEN
C
C Calling on a subroutine that checks whether the line is a blank one or not.
C If STRING is a blank line, then BLANK will return as .TRUE.
         CALL CHBLANK(BLANK, STRING)
         IF (BLANK) THEN
            GOTO 20
         ELSE
            READ (STRING, '(BN,2I5)') INTEXP, INTORB
            RETURN
         END IF
      ELSE IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
         NEWEL = .TRUE.
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
      RETURN
C
 1000 CONTINUE
      WRITE (LUPRI,'(/A,I5)') 'Error in reading file, your basis '//
     &     'is not complete. IOSTAT code:',IOERR
      CALL QUIT('Incomplete basis set file in HERBAS')
      END
C/* Deck FIND_BASFIL */
      SUBROUTINE FIND_BASFIL(BASNAM,LUBAS)
C
C     Joern Thyssen and Hans Joergen Aa. Jensen.
C     Last revision Feb. 2005
C
#include <implicit.h>
#include "priunit.h"
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basisfile                *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*     EXST   - Logical variable that helps inquire if   *
C*              there exists a file with that name.      *
C*********************************************************
C
      CHARACTER*(*) BASNAM
      CHARACTER*140 STRING
      LOGICAL EXST
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
C     The environment variable is a : separated string, e.g.:
C     /programs/Dirac/basis:/programs/Dalton/basis:.
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
C     In there are no colons, we are the last directory in BASDIR
C
      IF ( ICOLON .EQ. ISTART - 1 ) THEN
C
C        Find end of string
C
         ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:), ' ' )
C
C        If end of string = start of string there are no more
C        directories, so we gracefully die...
C
         IF ( ICOLON .EQ. ISTART ) GOTO 300
      END IF
      IEND = ICOLON - 1
      STRING = BASDIR(ISTART:IEND)
      ISTART = ICOLON + 1
C
C     Add trailing slash if not present.
C
      IEND = INDEX ( STRING(1:), ' ') - 1
      IF ( STRING(IEND:IEND) .NE. '/' ) THEN
         IEND = IEND + 1
         STRING(IEND:IEND) = '/'
      END IF
      STRING = STRING(1:IEND)//BASNAM
      IEND = INDEX(STRING(1:),' ') - 1
Cluuk Activate when debugging basis set search (should not be on by default)
CLuuk WRITE(LUPRI,'(A/,A1,A,A1)')
CLuuk&     ' Searching basis : trying file: ','"',STRING(1:IEND),'"'
C
C     Inquire if "/path/basis-set" exists.
C
      INQUIRE (FILE = STRING(1:IEND), EXIST = EXST)
      IF (EXST) THEN
         CALL GPOPEN(LUBAS,STRING(1:IEND),'OLD',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
         RETURN
      ELSE
         GO TO 10
      END IF
C
  300 IEND = INDEX ( BASNAM(1:), ' ') - 1
      WRITE (LUPRI,'(/3A)')
     +      'Basis "',BASNAM(1:IEND),'" doesn''t exist'
      IF (BASNAM(1:6) .EQ. 'Huckel') THEN
         WRITE (LUPRI,'(/A)')
     +   'NB! Basis Huckel must be available to use Huckel guess!'
      END IF
      CALL QUIT('Non-existing basis set in HERBAS')

      END
