C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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
*=====================================================================*
      SUBROUTINE CC_E1AIM(E1AMAT,RAIM,RBIM,FOCK,
     *                  XLAMDHA,ISYMLA,XLAMDHB,ISYMLB,
     *                  FCKCON,RTRAN,LRELAX,IOPT,ISYE1A)
*---------------------------------------------------------------------*
*
*     Transforms delta index of R intermediates to virtual and
*     calculate the E1A intermediates in the MO basis
*
*     IOPT = 1:     EMAT1 = FOCK   - XLAMDHA * RAIM
*                   
*                   RBIM,XLAMDHB,ISYMB is dummy input
*
*     IOPT = 2:     EMAT1 = FOCK   - XLAMDHA * RBIM - XLAMDHB * RAIM
*
*     RTRAN  = FALSE : skip R intermediate contributions
*     FCKCON = FALSE : skip FOCK matrix contribution
*     LRELAX = FALSE : skip contributions form XLAMDHB/XLAMDPB
*
*     Symmetries:    ISYE1A  --  E1AMAT, FOCK(MO), ONEHAM(MO)
*                    ISYMLA  --  XLAMDHA
*                    ISYMLB  --  XLAMDHB
*
*     Based on Christof's CC_EIM
*     Sonia Coriani 06/09-1999
*
*---------------------------------------------------------------------*
#include <implicit.h>
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DIMENSION E1AMAT(*),FOCK(*)
      DIMENSION RAIM(*),RBIM(*)
      DIMENSION XLAMDHA(*),XLAMDHB(*)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
C
      LOGICAL FCKCON, RTRAN, LRELAX
      INTEGER IOPT
C
      ISYRAIM = MULD2H(ISYMLA,ISYE1A)
      ISYRBIM = 0
      IF (IOPT.EQ.2) THEN
        ISYRAIM = MULD2H(ISYMLB,ISYE1A)
        ISYRBIM = MULD2H(ISYMLA,ISYE1A)
      END IF
C
C---------------------------------------------------------
C     Transform the delta index of R intermediate(s) to c.
C     store result in E1AMAT (R_bc)
C---------------------------------------------------------
C
      CALL DZERO(E1AMAT,NMATAB(ISYE1A))

      IF ( RTRAN ) THEN
C
         IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'CC_E1AIM> norm^2 of RAIM:',
     &       DDOT(NEMAT1(ISYRAIM),RAIM,1,RAIM,1)
           IF (IOPT.EQ.2) 
     &       WRITE (LUPRI,*) 'CC_E1AIM> norm^2 of RBIM:',
     &         DDOT(NEMAT1(ISYRBIM),RBIM,1,RBIM,1)
           CALL FLSHFO(LUPRI)
         END IF
C
         DO ISYMD = 1,NSYM
C
            ISYMC = MULD2H(ISYMD,ISYMLA)
            ISYMB = MULD2H(ISYMC,ISYE1A)
C
            NVIRB = MAX(NVIR(ISYMB),1)
            NBASD = MAX(NBAS(ISYMD),1)
C
            KOFF1 = IEMAT1(ISYMB,ISYMD) + 1
            KOFF2 = IGLMVI(ISYMD,ISYMC) + 1
            KOFF3 = IMATAB(ISYMB,ISYMC) + 1
C
            IF ( IOPT .EQ. 1) THEN

              CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),NBAS(ISYMD),
     *                   -ONE,RAIM(KOFF1),NVIRB,XLAMDHA(KOFF2),NBASD,
     *                   ONE,E1AMAT(KOFF3),NVIRB)

            ELSE

              CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),NBAS(ISYMD),
     *                   -ONE,RBIM(KOFF1),NVIRB,XLAMDHA(KOFF2),NBASD,
     *                   ONE,E1AMAT(KOFF3),NVIRB)

              IF (LRELAX) THEN

                ISYMC = MULD2H(ISYMD,ISYMLB)
                ISYMB = MULD2H(ISYMC,ISYE1A)

                NVIRB = MAX(NVIR(ISYMB),1)
                NBASD = MAX(NBAS(ISYMD),1)
C
                KOFF1 = IEMAT1(ISYMB,ISYMD) + 1
                KOFF2 = IGLMVI(ISYMD,ISYMC) + 1
                KOFF3 = IMATAB(ISYMB,ISYMC) + 1
C
                CALL DGEMM('N','N',NVIR(ISYMB),NVIR(ISYMC),NBAS(ISYMD),
     *                     -ONE,RAIM(KOFF1),NVIRB,XLAMDHB(KOFF2),NBASD,
     *                     ONE,E1AMAT(KOFF3),NVIRB)
              END IF

            END IF
C
         END DO
C
      ELSE
C
         CALL DZERO(E1AMAT,NMATAB(ISYE1A))
C
      END IF 
C
C---------------------------------------------------------
C     Add the Fock matrix contribution:
C---------------------------------------------------------
C
      IF (FCKCON) THEN
 
         DO ISYMC = 1,NSYM
 
            ISYMB = MULD2H(ISYMC,ISYE1A)
 
            DO C = 1,NVIR(ISYMC)
 
               KOFF1 = IFCVIR(ISYMB,ISYMC) + NORB(ISYMB)*(C - 1)
     *                                     + NRHF(ISYMB) + 1
               KOFF2 = IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C - 1) + 1
 
               CALL DAXPY(NVIR(ISYMB),ONE,FOCK(KOFF1),1,E1AMAT(KOFF2),1)
 
            END DO
         END DO
 
      ENDIF
 
      RETURN
      END
*=====================================================================*
