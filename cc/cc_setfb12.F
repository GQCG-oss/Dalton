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
c /* deck CC_SETFB12 */
*=====================================================================*
       SUBROUTINE CC_SETFB12(IFTRAN,IFDOTS,MXTRAN,MXVEC,
     &                       IZETAV,IOPER,IKAPPA,ITAMPA,ITAMPB,
     &                       ITRAN,IVEC)
*---------------------------------------------------------------------*
*
*    Purpose: maintains a list of dot products of F{O} matrix 
*             transformations with right amplitude vectors:
*                        (Z*F{O}*T^A) * T^B
*             assumes that T^A and T^B belong to different lists
*
*             IFTRAN - list of F matrix transformations
*             IFDOTS - list of vectors it should be dottet on
*        
*             MXTRAN - maximum list dimension
*             MXVEC  - maximum second dimension for IFDOTS
*      
*             IZETAV - index of lagrangian multiplier vector
*             IOPER  - index of property operator 
*             IKAPPA - index of the relaxation vector
*             ITAMPA - index of amplitude vector A
*             ITAMPB - index of amplitude vector B
*
*             ITRAN - index in IFTRAN list
*             IVEC  - second index in IFDOTS list
*
*    Written by Sonia Coriani, Maj 2000. Based on CC_SETFA and CC_SETFA12
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

      INTEGER MXVEC, MXTRAN
      INTEGER IFTRAN(5,MXTRAN)
      INTEGER IFDOTS(MXVEC,MXTRAN)

      LOGICAL LFNDB
      INTEGER IZETAV, IOPER, IKAPPA, ITAMPA, ITAMPB
      INTEGER ITRAN, IVEC
      INTEGER ITAMP, I, IDX

* statement  functions:
      LOGICAL LFATST, LFAEND
      INTEGER IL, IA, IO, IK
      LFATST(ITRAN,IL,IO,IK,IA) = IFTRAN(1,ITRAN).EQ.IL 
     &       .AND. IFTRAN(2,ITRAN).EQ.IO .AND. IFTRAN(3,ITRAN).EQ.IA 
     &       .AND. IFTRAN(5,ITRAN).EQ.IK
      LFAEND(ITRAN) = ITRAN.GT.MXTRAN .OR.
     &      (IFTRAN(1,ITRAN)+IFTRAN(2,ITRAN)+IFTRAN(3,ITRAN)).LE.0 


*---------------------------------------------------------------------*
* set up list of F{B} matrix transformations
*---------------------------------------------------------------------*
      ITRAN = 1
      LFNDB  = LFATST(ITRAN,IZETAV,IOPER,IKAPPA,ITAMPA)

      DO WHILE ( .NOT. (LFNDB.OR.LFAEND(ITRAN)))
       ITRAN = ITRAN + 1
       LFNDB  = LFATST(ITRAN,IZETAV,IOPER,IKAPPA,ITAMPA)
      END DO

      IF (.NOT.LFNDB) THEN
        IFTRAN(1,ITRAN) = IZETAV
        IFTRAN(2,ITRAN) = IOPER
        IFTRAN(3,ITRAN) = ITAMPA
        IFTRAN(4,ITRAN) = 0
        IFTRAN(5,ITRAN) = IKAPPA
        ITAMP = ITAMPB
      ELSE 
        IF (LFNDB) ITAMP = ITAMPB
      END IF

      IVEC = 1
      DO WHILE (IFDOTS(IVEC,ITRAN).NE.ITAMP .AND.
     &            IFDOTS(IVEC,ITRAN).NE.0 .AND. IVEC.LE.MXVEC)
        IVEC = IVEC + 1
      END DO

      IFDOTS(IVEC,ITRAN) = ITAMP

C     WRITE (LUPRI,*) 'CC_SETFB12>',IZETAV,IOPER,ITAMPA,ITAMPB,ITRAN,IVEC
*---------------------------------------------------------------------*
      IF (IVEC.GT.MXVEC .OR. ITRAN.GT.MXTRAN) THEN
        WRITE (LUPRI,*) 'IVEC :',IVEC
        WRITE (LUPRI,*) 'ITRAN:',ITRAN
        WRITE (LUPRI,*) 'IOPER,IKAPPA :',IOPER,IKAPPA
        WRITE (LUPRI,*) 'ITAMPA,ITAMPB:',ITAMPA,ITAMPB
        IDX = 1
        DO WHILE ( .NOT. LFAEND(IDX) )
          WRITE(LUPRI,'(A,5I5,5X,(12I5,20X))') 'CC_SETFA12>',
     &       (IFTRAN(I,IDX),I=1,5),(IFDOTS(I,IDX),I=1,MXVEC)
          IDX = IDX + 1
        END DO
        CALL FLSHFO(LUPRI)
        CALL QUIT('Overflow error in CC_SETFA12.')
      END IF
      
      RETURN
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE CC_SETFB12                           *
*---------------------------------------------------------------------*
