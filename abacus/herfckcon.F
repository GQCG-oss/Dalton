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
C  /* Deck fckcon */
      SUBROUTINE FCKCON(FMAT,DMAT,NDMAT,I,AOINT,IND,
     &        NINT,NIND,NBUF,IX,IY1,FAC)
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
C  This routine tries to replace INTFC1(her2fck.F), FOKDI1(eri2fck.F)
C  and FCKDS1, FCKDS2,FCKDT3 (sirfck.F)
C
C We keep to outer loop over NDMAT(in calling routine), but we vectorize with
C IVECLN / NDMAT temoprary matrices. This will use the same 
C amount of scratch memory we use in fokdi1 in eri. 
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <nuclei.h>
#include <dftcom.h>
#include <veclen.h>
      PARAMETER (D0=0.0D0,DP25=0.25D0,DP5=0.5D0)
      INTEGER A,B,C,D
CNECgh990830
#if defined (VAR_VECTOR)
      DIMENSION FMAT(NBASIS+NODD,NBASIS,NDMAT,*)
#else
      PARAMETER (K = 1)
      DIMENSION FMAT(NBASIS,NBASIS,NDMAT,K)
#endif        
      DIMENSION DMAT(NBASIS,NBASIS,NDMAT), AOINT(NINT),IND(NIND,4)
C
C
      IF (IX.NE.2) THEN
        IY = IY1
      ELSE
C       ... No Coulomb contribution for
C           antisymmetric density matrix
        IY = IY1 - MOD(IY1,2)
      END IF
C
      IF (HFXFAC.EQ.D0 .AND. HFXMU.EQ.D0) THEN
        IY = MOD(IY,2)
C       ... only Coulomb contribution
      ELSE
         DP5X  =  DP5*HFXFAC
         DP25X = DP25*HFXFAC
      END IF
      IF (HFXMU.NE.D0) IY = 2 
      IF (IY .EQ. 0) GO TO 9999
C     ... nothing to do!
C
#if defined (VAR_VECTOR)
        ICHUNK = MAX (IVECLN/NDMAT,1)
#endif
        IF    (IX.EQ.1.OR.IX.EQ.2) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
#if defined (VAR_VECTOR)
            J=0
            DO INT = 1,(NBUF/ICHUNK)
#include <ivdep.h>
              DO K=1,ICHUNK
                J = J + 1
#else
            DO J = 1,NBUF
#endif
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                EINT = -DP25X*DINT
                FMAT(A,B,I,K) = FMAT(A,B,I,K) + DINT*DMAT(C,D,I)
                FMAT(C,D,I,K) = FMAT(C,D,I,K) + DINT*DMAT(A,B,I)
                FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
                FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
                FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
                FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              END DO
#if defined (VAR_VECTOR)
            END DO
#include <ivdep.h>
              DO K=1,MOD(NBUF,ICHUNK)
                J = J + 1
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                EINT = -DP25X*DINT
                FMAT(A,B,I,K) = FMAT(A,B,I,K) + DINT*DMAT(C,D,I)
                FMAT(C,D,I,K) = FMAT(C,D,I,K) + DINT*DMAT(A,B,I)
                FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
                FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
                FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
                FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              END DO
#endif
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (1/4) (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
#if defined (VAR_VECTOR)
            J=0
            DO INT = 1,(NBUF/ICHUNK)
#include <ivdep.h>
              DO K=1,ICHUNK
                J = J + 1
#else
            DO J = 1,NBUF
#endif
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                EINT = -DP25X*DINT
                FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
                FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
                FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
                FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              END DO
#if defined (VAR_VECTOR)
            END DO
#include <ivdep.h>
              DO K=1,MOD(NBUF,ICHUNK)
                J = J + 1
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                EINT = -DP25X*DINT
                FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
                FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
                FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
                FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              END DO
#endif
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
#if defined (VAR_VECTOR)
            J=0
            DO INT = 1,(NBUF/ICHUNK)
#include <ivdep.h>
              DO K=1,ICHUNK
                J = J + 1
#else
            DO J = 1,NBUF
#endif
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                FMAT(A,B,I,K) = FMAT(A,B,I,K) + DINT*DMAT(C,D,I)
                FMAT(C,D,I,K) = FMAT(C,D,I,K) + DINT*DMAT(A,B,I)
              END DO
#if defined (VAR_VECTOR)
            END DO
#include <ivdep.h>
            DO K=1,MOD(NBUF,ICHUNK)
                J = J + 1
                A = IND(J,1)
                B = IND(J,2)
                C = IND(J,3)
                D = IND(J,4)
                DINT = AOINT(J) * FAC
                FMAT(A,B,I,K) = FMAT(A,B,I,K) + DINT*DMAT(C,D,I)
                FMAT(C,D,I,K) = FMAT(C,D,I,K) + DINT*DMAT(A,B,I)
            END DO
#endif
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'FCKCON ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IX*10+IY,
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in FCKCON: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          IF(IY.EQ.3) THEN
#if defined (VAR_VECTOR)
            J=0
            DO INT = 1,(NBUF/ICHUNK)
#include <ivdep.h>
              DO K=1,ICHUNK
                J = J + 1
#else
            DO J = 1,NBUF
#endif
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              GCD  = DINT*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I,K) = FMAT(A,B,I,K) + GCD
              FMAT(B,A,I,K) = FMAT(B,A,I,K) + GCD
              GAB  = DINT*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I,K) = FMAT(C,D,I,K) + GAB
              FMAT(D,C,I,K) = FMAT(D,C,I,K) + GAB
              EINT = -DP5X*DINT
              FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
              FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
              FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
              FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              FMAT(A,C,I,K) = FMAT(A,C,I,K) + EINT*DMAT(B,D,I)
              FMAT(A,D,I,K) = FMAT(A,D,I,K) + EINT*DMAT(B,C,I)
              FMAT(B,C,I,K) = FMAT(B,C,I,K) + EINT*DMAT(A,D,I)
              FMAT(B,D,I,K) = FMAT(B,D,I,K) + EINT*DMAT(A,C,I)
             ENDDO
#if defined (VAR_VECTOR)
            ENDDO
#include <ivdep.h>
             DO K=1,MOD(NBUF,ICHUNK)
              J = J + 1
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              GCD  = DINT*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I,K) = FMAT(A,B,I,K) + GCD
              FMAT(B,A,I,K) = FMAT(B,A,I,K) + GCD
              GAB  = DINT*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I,K) = FMAT(C,D,I,K) + GAB
              FMAT(D,C,I,K) = FMAT(D,C,I,K) + GAB
              EINT = -DP5X*DINT
              FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
              FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
              FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
              FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              FMAT(A,C,I,K) = FMAT(A,C,I,K) + EINT*DMAT(B,D,I)
              FMAT(A,D,I,K) = FMAT(A,D,I,K) + EINT*DMAT(B,C,I)
              FMAT(B,C,I,K) = FMAT(B,C,I,K) + EINT*DMAT(A,D,I)
              FMAT(B,D,I,K) = FMAT(B,D,I,K) + EINT*DMAT(A,C,I)
             ENDDO
#endif
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
#if defined (VAR_VECTOR)
            J=0
            DO INT = 1,(NBUF/ICHUNK)
#include <ivdep.h>
              DO K=1,ICHUNK
                J = J + 1
#else
            DO J = 1,NBUF
#endif
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              EINT = -DP5X*DINT
              FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
              FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
              FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
              FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              FMAT(A,C,I,K) = FMAT(A,C,I,K) + EINT*DMAT(B,D,I)
              FMAT(A,D,I,K) = FMAT(A,D,I,K) + EINT*DMAT(B,C,I)
              FMAT(B,C,I,K) = FMAT(B,C,I,K) + EINT*DMAT(A,D,I)
              FMAT(B,D,I,K) = FMAT(B,D,I,K) + EINT*DMAT(A,C,I)
            END DO
#if defined (VAR_VECTOR)
            END DO
#include <ivdep.h>
            DO K = 1,MOD(NBUF,ICHUNK)
              J = J + 1
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              EINT = -DP5X*DINT
              FMAT(C,A,I,K) = FMAT(C,A,I,K) + EINT*DMAT(D,B,I)
              FMAT(D,A,I,K) = FMAT(D,A,I,K) + EINT*DMAT(C,B,I)
              FMAT(C,B,I,K) = FMAT(C,B,I,K) + EINT*DMAT(D,A,I)
              FMAT(D,B,I,K) = FMAT(D,B,I,K) + EINT*DMAT(C,A,I)
              FMAT(A,C,I,K) = FMAT(A,C,I,K) + EINT*DMAT(B,D,I)
              FMAT(A,D,I,K) = FMAT(A,D,I,K) + EINT*DMAT(B,C,I)
              FMAT(B,C,I,K) = FMAT(B,C,I,K) + EINT*DMAT(A,D,I)
              FMAT(B,D,I,K) = FMAT(B,D,I,K) + EINT*DMAT(A,C,I)
            END DO
#endif
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
#if defined (VAR_VECTOR)
            J=0
            DO INT = 1,(NBUF/ICHUNK)
#include <ivdep.h>
              DO K=1,ICHUNK
                J = J + 1
#else
            DO J = 1,NBUF
#endif
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              GCD  = DINT*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I,K) = FMAT(A,B,I,K) + GCD
              FMAT(B,A,I,K) = FMAT(B,A,I,K) + GCD
              GAB  = DINT*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I,K) = FMAT(C,D,I,K) + GAB
              FMAT(D,C,I,K) = FMAT(D,C,I,K) + GAB
            END DO
#if defined (VAR_VECTOR)
            END DO
#include <ivdep.h>
             DO K=1,MOD(NBUF,ICHUNK)
              J = J + 1
              A = IND(J,1)
              B = IND(J,2)
              C = IND(J,3)
              D = IND(J,4)
              DINT = AOINT(J) * FAC
              GCD  = DINT*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I,K) = FMAT(A,B,I,K) + GCD
              FMAT(B,A,I,K) = FMAT(B,A,I,K) + GCD
              GAB  = DINT*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I,K) = FMAT(C,D,I,K) + GAB
              FMAT(D,C,I,K) = FMAT(D,C,I,K) + GAB
            END DO
#endif
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'FCKCON ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IX*10+IY1,
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in FCKCON: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      ' ERROR, specified IFCTYP not implemented yet',
     &      '        specified IFCTYP was',IX*10+IY1,
     &      '        for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in FCKCON: specified IFCTYP not implemented.')
        ENDIF
C
 9999 CONTINUE
      RETURN
C
      END
