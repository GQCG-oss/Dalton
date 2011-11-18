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
C  /* Deck setsoppa */
      SUBROUTINE SETSOPPA
C
C  Set up SOPPA-variables
C
C  Written 5/4-94 by Erik K. Dalskov
C
C  Modified by K.Ruud in Dec.-96 in order to avoid static allocation of 
C  large particle-hole matrices
C
#include <implicit.h>
C
C Used from common blocks:
C   INFORB : NSYM,MULD2H(8,8),NRHF(),NVIR(),...
C   INFMP2 : ...,NFRMP2(),IFRMP2(NORBT),...
C
#include <maxorb.h>
#include <maxash.h>
#include <infpri.h>
#include <inforb.h>
#include <infmp2.h>
#include <infsop.h>
#include <infind.h>
#include <infdim.h>
#include <iratdef.h>
C
      CALL QENTER('SETSOPPA')
C
C  Find offsets in XINDX for AB and IJ address arrays and for
C  the A(2) matrix. They will be placed in infsop.h.
C
      KABSAD = 1
      KABTAD = KABSAD + (N2ORBX + 1)/IRAT
      KIJSAD = KABTAD + (N2ORBX + 1)/IRAT
      KIJTAD = KIJSAD + (N2ISHX * NSYM + 1)/IRAT
      KIJ1AD = KIJTAD + (N2ISHX * NSYM + 1)/IRAT
      KIJ2AD = KIJ1AD + (N2ISHX * NSYM + 1)/IRAT
      KIJ3AD = KIJ2AD + (N2ISHX * NSYM + 1)/IRAT
      KIADR1 = KIJ3AD + (N2ISHX * NSYM + 1)/IRAT
      KAB2   = KIADR1 + (MAXOCC * MAXVIR + 1)/IRAT
C
C  A(2) matrix does not exist in XINDX yet
C
      A2EXIST = .FALSE.
C
C  Redefine some variables in infdim:
C
      LCINDX = KAB2 + N2ORBX
C  ...no PHP for SOPPA
      MAXPHP = 0
C
C  Find no. of vir/vir-pairs with a>=b (NSVV(SYM)) and a>b (NTVV(SYM))
C  and no. of occ/occ-pairs with i>=j (NSOO(SYM)) and i>j (NTOO(SYM))
C
      DO I=1,NSYM
         NSOO(I) = 0
         NTOO(I) = 0
         NSVV(I) = 0
         NTVV(I) = 0
      ENDDO
C
      DO IJSYM=2,NSYM
         DO ISYM = 1,NSYM
            JSYM = MULD2H(ISYM,IJSYM)
            IF (ISYM .LT. JSYM) GOTO 100
            NHOLEI = NRHF(ISYM) - NFRMP2(ISYM)
            NHOLEJ = NRHF(JSYM) - NFRMP2(JSYM)
            NSOO(IJSYM) = NSOO(IJSYM) + NHOLEI*NHOLEJ
            NTOO(IJSYM) = NTOO(IJSYM) + NHOLEI*NHOLEJ
            NSVV(IJSYM) = NSVV(IJSYM) + NVIR(ISYM)*NVIR(JSYM)
            NTVV(IJSYM) = NTVV(IJSYM) + NVIR(ISYM)*NVIR(JSYM)
 100        CONTINUE
         ENDDO
      ENDDO
C
      DO ISYM=1,NSYM
         NHOLEI = NRHF(ISYM) - NFRMP2(ISYM)
         NSOO(1) = NSOO(1) + NHOLEI*(NHOLEI + 1) / 2
         NTOO(1) = NTOO(1) + NHOLEI*(NHOLEI - 1) / 2
         NSVV(1) = NSVV(1) + NVIR(ISYM)*(NVIR(ISYM) + 1) / 2
         NTVV(1) = NTVV(1) + NVIR(ISYM)*(NVIR(ISYM) - 1) / 2
      ENDDO
C
C  Find total no. of 2p-2h variables in each symmetry
C  (singlet and triplet case)
C
      DO I=1,NSYM
         NS2P2H(I) = 0
         NT2P2H(I) = 0
         N12P2H(I) = 0
         N22P2H(I) = 0
         N32P2H(I) = 0
      ENDDO

C
      DO ISYM=1,NSYM
         DO JSYM=1,NSYM
            IJSYM = MULD2H(ISYM,JSYM)
            NS2P2H(IJSYM) = NS2P2H(IJSYM) + NSVV(ISYM) * NSOO(JSYM)
            NT2P2H(IJSYM) = NT2P2H(IJSYM) + NTVV(ISYM) * NTOO(JSYM)
            N12P2H(IJSYM) = N12P2H(IJSYM) + NTVV(ISYM) * NTOO(JSYM)
            N22P2H(IJSYM) = N22P2H(IJSYM) + NSVV(ISYM) * NTOO(JSYM)
            N32P2H(IJSYM) = N32P2H(IJSYM) + NTVV(ISYM) * NSOO(JSYM)
         ENDDO
      ENDDO
C
      DO I=1,NSYM
         N2P2HS(I) = NS2P2H(I) + NT2P2H(I)
         N2P2HT(I) = N12P2H(I) + N22P2H(I) + N32P2H(I)
      ENDDO
      CALL QEXIT('SETSOPPA')
      RETURN
      END
C
C  /* Deck set2soppa */
      SUBROUTINE SET2SOPPA(ABSADR,ABTADR,IJSADR,IJTADR,IJ1ADR,IJ2ADR,
     &                     IJ3ADR,IADR1)
#include <implicit.h>
      INTEGER A,AA,B,BEND,ASYM,BSYM,ABSYM,ABSNDX,ABTNDX,ABSADR,ABTADR
      PARAMETER (IBIG = -100 000 000)
#include <maxorb.h>
#include <maxash.h>
#include <infpri.h>
#include <inforb.h>
#include <infsop.h>
#include <infmp2.h>
#include <infind.h>
C
      DIMENSION ABSADR(NORBT,NORBT), ABTADR(NORBT,NORBT),
     &          IJSADR(NISHT,NISHT,NSYM), IJTADR(NISHT,NISHT,NSYM),
     &          IJ1ADR(NISHT,NISHT,NSYM), IJ2ADR(NISHT,NISHT,NSYM),
     &          IJ3ADR(NISHT,NISHT,NSYM), IADR1(MAXVIR,MAXOCC)
C
      CALL QENTER('SET2SOPPA')
C
C  Find offset-arrays for the symmetry packed 2p-2h vectors in SOPPA
C  (singlet properties) and new IJ offsets for triplet properties
C
C  Moved out of SETSOPPA by K.Ruud, Dec.-96
C
      DO ABSYM=1,NSYM
         ABSNDX = 0
         ABTNDX = 0
         DO A=1,NP
            AA = IPHORD(NH+A)
            ASYM = ISMO(AA)
            BSYM = MULD2H(ASYM,ABSYM)
            IF (BSYM .LE. ASYM) THEN
               BEND = MIN(ISSH(BSYM)+NVIR(BSYM),A)
               DO B=ISSH(BSYM)+1,BEND
                  IF (A .EQ. B) THEN
                     ABSNDX = ABSNDX + 1
                     ABSADR(A,A) = ABSNDX
                     ABTADR(A,A) = IBIG
                  ELSE
                     ABSNDX = ABSNDX + 1
                     ABTNDX = ABTNDX + 1
                     ABSADR(A,B) = ABSNDX
                     ABTADR(A,B) = ABTNDX
                     ABSADR(B,A) = ABSNDX
                     ABTADR(B,A) = ABTNDX
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
C
C
      DO LSYMOP=1,NSYM
         IJSNDX = 0
         IJTNDX = NS2P2H(LSYMOP)
         IJ1NDX = 0
         IJ2NDX = N12P2H(LSYMOP)
         IJ3NDX = N22P2H(LSYMOP) + N12P2H(LSYMOP)
         DO I=1,NH
            II = IPHORD(I)
            ISYM = ISMO(II)
            DO J=1,I
               JJ = IPHORD(J)
               JSYM = ISMO(JJ)
               IJSYM = MULD2H(ISYM,JSYM)
               ABSYM = MULD2H(IJSYM,LSYMOP)
               IF (I .EQ. J) THEN
                  IF (IFRMP2(II) .EQ. 0) THEN
                     IJSADR(I,I,LSYMOP) = IJSNDX
                     IJTADR(I,I,LSYMOP) = IBIG
                     IJ1ADR(I,I,LSYMOP) = IBIG
                     IJ2ADR(I,I,LSYMOP) = IBIG
                     IJ3ADR(I,I,LSYMOP) = IJ3NDX
                     IJSNDX = IJSNDX + NSVV(ABSYM)
                     IJ3NDX = IJ3NDX + NTVV(ABSYM)
                  ELSE
                     IJSADR(I,I,LSYMOP) = IBIG
                     IJTADR(I,I,LSYMOP) = IBIG
                     IJ1ADR(I,I,LSYMOP) = IBIG
                     IJ2ADR(I,I,LSYMOP) = IBIG
                     IJ3ADR(I,I,LSYMOP) = IBIG
                  END IF
               ELSE
                  IF (IFRMP2(II) .EQ. 0 .AND. IFRMP2(JJ) .EQ. 0) THEN
                     IJSADR(I,J,LSYMOP) = IJSNDX
                     IJSADR(J,I,LSYMOP) = IJSNDX
                     IJTADR(I,J,LSYMOP) = IJTNDX
                     IJTADR(J,I,LSYMOP) = IJTNDX
                     IJ1ADR(I,J,LSYMOP) = IJ1NDX
                     IJ1ADR(J,I,LSYMOP) = IJ1NDX
                     IJ2ADR(I,J,LSYMOP) = IJ2NDX
                     IJ2ADR(J,I,LSYMOP) = IJ2NDX
                     IJ3ADR(I,J,LSYMOP) = IJ3NDX
                     IJ3ADR(J,I,LSYMOP) = IJ3NDX
                     IJSNDX = IJSNDX + NSVV(ABSYM)
                     IJTNDX = IJTNDX + NTVV(ABSYM)
                     IJ1NDX = IJ1NDX + NTVV(ABSYM)
                     IJ2NDX = IJ2NDX + NSVV(ABSYM)
                     IJ3NDX = IJ3NDX + NTVV(ABSYM)
                  ELSE
                     IJSADR(I,J,LSYMOP) = IBIG
                     IJSADR(J,I,LSYMOP) = IBIG
                     IJTADR(I,J,LSYMOP) = IBIG
                     IJTADR(J,I,LSYMOP) = IBIG
                     IJ1ADR(I,J,LSYMOP) = IBIG
                     IJ1ADR(J,I,LSYMOP) = IBIG
                     IJ2ADR(I,J,LSYMOP) = IBIG
                     IJ2ADR(J,I,LSYMOP) = IBIG
                     IJ3ADR(I,J,LSYMOP) = IBIG
                     IJ3ADR(J,I,LSYMOP) = IBIG
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C  Find also IADR1 (the copy from the call to MP2SET has been overwritten)
C
      CALL PHADR1(IADR1,1)
C
C  End of SET2SOPPA
C
      CALL QEXIT('SET2SOPPA')
      RETURN
      END
C  /* Deck a0s2 */
      SUBROUTINE A0S2(FCONE,FC,DONEPT,ZYMAT,TZYMAT,STWO,
     *                ORBE,EPSIL,NSIM,WRK,LWRK)
C
C Copyright by Martin Packer 251193. In formulas below
C i,j,k.. are occupied and a,b,c... are virtual orbitals.
C We must form FCONE = 0.5*(A(0)S(2)ZYMAT + S(2)A(0)ZYMAT)
C IN ORDER TO ACCOUNT FOR THE NON-HERMITICITY.
C
C 940726-hjaaj: corrected code for NSIM .gt. 1
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      DIMENSION FCONE(N2ORBX,*), ORBE(NORBT,NORBT), FC(*)
      DIMENSION WRK(*), ZYMAT(NORBT,NORBT,*), TZYMAT(NORBT,NORBT,*)
      DIMENSION DONEPT(NORBT,NORBT), STWO(NORBT,NORBT), EPSIL(NORBT)
C
      PARAMETER( DMP25 = -0.25D0 )
C
C
C First PUT THE ORBITAL ENERGIES FROM FC INTO THE ARRAY ORBE
C
      IS = 0
      DO IM = 1,NSYM
         IORBM = IORB(IM)
         NORBM = NORB(IM)
         IIORBM = IIORB(IM)
         IF(NORBM.NE.0) THEN
            DO IJ = 1,NORBM
               IS = IS + 1
               IPL = IIORBM + (IJ*(IJ+1))/2
               EPSIL(IS) = FC(IPL)
            ENDDO
         ENDIF
      ENDDO
C
C FORM S(2).ZYMAT MATRIX AND PLACE IN STWO
C
      DO I = 1,NSIM
         CALL DZERO(STWO,N2ORBX)
         DO IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,KSYMOP)
            NOCCA = NOCC(IASYM)
            NSSHB = NSSH(IBSYM)
            IF ((NOCCA.NE.0).AND.(NSSHB.NE.0)) THEN
C S(AI) TERMS
               IAST = IORB(IASYM) + 1
               IBST = IORB(IBSYM) + 1 + NOCC(IBSYM)
               CALL DGEMM('N','N',NSSHB,NOCCA,NOCCA,1.D0,
     &                    ZYMAT(IBST,IAST,I),NORBT,
     &                    DONEPT(IAST,IAST),NORBT,1.D0,
     &                    STWO(IBST,IAST),NORBT)
               CALL DGEMM('N','N',NSSHB,NOCCA,NSSHB,-1.D0,
     &                    DONEPT(IBST,IBST),NORBT,
     &                    ZYMAT(IBST,IAST,I),NORBT,1.D0,
     &                    STWO(IBST,IAST),NORBT)
C S(IA) TERMS
               CALL DGEMM('T','N',NOCCA,NSSHB,NSSHB,1.D0,
     &                    TZYMAT(IBST,IAST,I),NORBT,
     &                    DONEPT(IBST,IBST),NORBT,1.D0,
     &                    STWO(IAST,IBST),NORBT)
               CALL DGEMM('N','T',NOCCA,NSSHB,NOCCA,-1.D0,
     &                    DONEPT(IAST,IAST),NORBT,
     &                    TZYMAT(IBST,IAST,I),NORBT,1.D0,
     &                    STWO(IAST,IBST),NORBT)
            ENDIF
         ENDDO
C
C NOW FORM A(0)S(2).ZYMAT = G, AND STORE IN STWO:
C    G = SUM(JB) (EPSILON(A)-EPSILON(I))*S(AI,JB).ZYMAT(BJ)
C and FORM A(0).ZYMAT = H, AND STORE IN ORBE:
C    H(J,B) =  (EPSILON(B)-EPSILON(J))ZYMAT(BJ)
C
C        CALL DZERO(ORBE,N2ORBX)
         DO ICSYM = 1,NSYM
            IDSYM = MULD2H(ICSYM,KSYMOP)
            IORBC = IORB(ICSYM)
            IORBD = IORB(IDSYM)
            NORBD = NORB(IDSYM)
            NOCCC = NOCC(ICSYM)
            NOCCD = NOCC(IDSYM)
            DO IK = (IORBC + 1),(IORBC + NOCCC)
               DO IL = (IORBD + 1 + NOCCD),(IORBD + NORBD)
                  STWO(IK,IL) =
     *               STWO(IK,IL)*(EPSIL(IL)-EPSIL(IK))
                  STWO(IL,IK) =
     *               STWO(IL,IK)*(EPSIL(IL)-EPSIL(IK))
                  ORBE(IK,IL) =
     *               ZYMAT(IK,IL,I)*(EPSIL(IL)-EPSIL(IK))
                  ORBE(IL,IK) =
     *               ZYMAT(IL,IK,I)*(EPSIL(IL)-EPSIL(IK))
               ENDDO
            ENDDO
         ENDDO
C
C FORM S(2)*ORBE, WHICH IS THE CONJUGATE TERM TO H
         DO IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,KSYMOP)
            NOCCA = NOCC(IASYM)
            NSSHB = NSSH(IBSYM)
            IF ((NOCCA.NE.0).AND.(NSSHB.NE.0)) THEN
C S(AI) TERMS
               IAST = IORB(IASYM) + 1
               IBST = IORB(IBSYM) + 1 + NOCC(IBSYM)
               CALL DGEMM('N','N',NSSHB,NOCCA,NOCCA,1.D0,
     &                    ORBE(IBST,IAST),NORBT,
     &                    DONEPT(IAST,IAST),NORBT,1.D0,
     &                    STWO(IBST,IAST),NORBT)
               CALL DGEMM('N','N',NSSHB,NOCCA,NSSHB,-1.D0,
     &                    DONEPT(IBST,IBST),NORBT,
     &                    ORBE(IBST,IAST),NORBT,1.D0,
     &                    STWO(IBST,IAST),NORBT)
C S(IA) TERMS
               CALL DGEMM('N','N',NOCCA,NSSHB,NSSHB,1.D0,
     &                    ORBE(IAST,IBST),NORBT,
     &                    DONEPT(IBST,IBST),NORBT,1.D0,
     &                    STWO(IAST,IBST),NORBT)
               CALL DGEMM('N','N',NOCCA,NSSHB,NOCCA,-1.D0,
     &                    DONEPT(IAST,IAST),NORBT,
     &                    ORBE(IAST,IBST),NORBT,1.D0,
     &                    STWO(IAST,IBST),NORBT)
            ENDIF
         ENDDO
C PLACE STWO INTO FCONE, WITH THE CORRECT SCALE FACTOR
         CALL DAXPY(N2ORBX,DMP25, STWO,1, FCONE(1,I),1)
      ENDDO
      RETURN
      END
C  /* Deck sopudv */
      SUBROUTINE SOPUDV(DONEPT)
C
C Copyright by Martin Packer 251193. In formulas below
C i,j,k.. are occupied and a,b,c... are virtual orbitals.
C
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
C
      DIMENSION  DONEPT(NORBT,NORBT)
C
      PARAMETER(D2 = 2.0D0)
C
C REMOVE THE DIAGONAL S(0) CONTRIBUTION FROM DONEPT (ADDED IN MP2FAC)
      DO ISY = 1,NSYM
         IORBSY = IORB(ISY)
         NOCCSY = NOCC(ISY)
         NORBSY = NORB(ISY)
         IF (NORBSY.NE.0) THEN
            DO IA = (IORBSY+1),(IORBSY+NOCCSY)
               DONEPT(IA,IA) = DONEPT(IA,IA) - D2
            ENDDO
         ENDIF
      ENDDO
      RETURN
      END
C  /* Deck hrpa */
      SUBROUTINE HRPA (FCONE,DONEPT,COEMP2,H2,ZYMAT,TZYMAT,H2X,H2XP,
     *                 COEUNP,ICI,IDI,ICDSYM,ITYPCD,NSIM,IADR1,WRK,LWRK)
C
C Copyright by Martin Packer 191193. In formulas below
C i,j,k.. are occupied and a,b,c... are virtual orbitals.
C
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
#include <inforb.h>
#include <orbtypdef.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      DIMENSION H2(NORBT,NORBT), FCONE(N2ORBX,*), COEMP2(*)
      DIMENSION DONEPT(NORBT,NORBT)
      DIMENSION WRK(*), ZYMAT(NORBT,NORBT,*), TZYMAT(NORBT,NORBT,*)
      DIMENSION H2X(NORBT,NORBT), H2XP(NORBT,NORBT), COEUNP(*)
C
      CALL QENTER('HRPA  ')
C
C
C     use distribution type ITYPCD =
C     1:inactive-inactive  2:active-inactive  3:active-active
C     4:secondary-inactive 5:secondary-active 6:secondary-secondary
C     We only need types 1, 4 and 6.
C
      IF (ITYPCD .EQ. 1) THEN
         IF (IFRMP2(ICI) .NE. 0 .AND. IFRMP2(IDI) .NE. 0) GO TO 9999
      ELSE IF (ITYPCD .EQ. 4) THEN
         IF (IFRMP2(ICI) .NE. 0) GO TO 9999
      ELSE IF (ITYPCD .NE. 6) THEN
         GO TO 9999
      END IF
C
C One index transform the two-e integrals
      DO I = 1,NSIM
        CALL DZERO(H2X ,N2ORBX)
        IF (ITYPCD .NE. 4) CALL DZERO(H2XP,N2ORBX)
        DO IPSYM = 1,NSYM
         IBSYM = MULD2H(IPSYM,ICDSYM)
         IASYM = MULD2H(IPSYM,KSYMOP)
         NOCCP = NOCC(IPSYM)
         NOCCB = NOCC(IBSYM)
         NOCCA = NOCC(IASYM)
         NSSHP = NSSH(IPSYM)
         NSSHB = NSSH(IBSYM)
         NSSHA = NSSH(IASYM)
cmjp branch to code for each distribution type
         IF (ITYPCD .EQ. 6) GOTO 600
         IF (ITYPCD .EQ. 4) GOTO 400
         IF (ITYPCD .NE. 1) GOTO 9999
C
cmjp occupied-occupied distributions
C
C This is the transform of (CD|**) with C,D occupied
C H2X(K,C) = SUM(M) H2(K,M)ZYMAT(M,C) - SUM(E) ZYMAT(K,E)H2(E,C): F(AI)
C H2X(C,L) = SUM(E)H2(C,E)ZYMAT(E,L) -SUM(M)ZYMAY(C,M)H2(M,L): F(IA)
C
C TRANSFORM OCC-OCC PART OF H2
C
         IF (NOCCB .GT. 0 .AND. NOCCP .GT. 0) THEN
            IAST  = IORB(IASYM) + 1 + NOCCA
            IBST  = IORB(IBSYM) + 1
            IPST  = IORB(IPSYM) + 1
C Now transform for those parts which enter F(ai)
            CALL DGEMM('N','N',NSSHA,NOCCB,NOCCP,1.D0,
     &                 TZYMAT(IAST,IPST,I),NORBT,
     &                 H2(IPST,IBST),NORBT,1.D0,
     &                 H2XP(IAST,IBST),NORBT)
C           CALL AMPAB(H2    (IBST,IPST),NOCCB,NOCCP,NORBT,NORBT,
C    *                 ZYMAT (IPST,IAST,I),NOCCP,NSSHA,NORBT,NORBT,
C    *                 H2X   (IBST,IAST),NORBT,NORBT )
C Now transform for those parts which enter F(ia)
            CALL DGEMM('N','N',NSSHA,NOCCB,NOCCP,-1.D0,
     &                 ZYMAT(IAST,IPST,I),NORBT,
     &                 H2(IPST,IBST),NORBT,1.D0,
     &                 H2X(IAST,IBST),NORBT)
         END IF
C TRANSFORM VIRT-VIRT PART OF H2
         IF (NOCCA .GT. 0) THEN
            IAST  = IORB(IASYM) + 1
            IBST  = IORB(IBSYM) + 1 + NOCCB
            IPST  = IORB(IPSYM) + 1 + NOCCP
C Now transform for those parts which enter F(ia)
            CALL DGEMM('N','N',NSSHB,NOCCA,NSSHP,1.D0,
     &                 H2(IBST,IPST),NORBT,
     &                 ZYMAT(IPST,IAST,I),NORBT,1.D0,
     &                 H2X(IBST,IAST),NORBT)
C Now transform for those parts which enter F(ai)
            CALL DGEMM('N','N',NSSHB,NOCCA,NSSHP,-1.D0,
     &                 H2(IBST,IPST),NORBT,
     &                 TZYMAT(IPST,IAST,I),NORBT,1.D0,
     &                 H2XP(IBST,IAST),NORBT)
C           CALL SMPAB( ZYMAT(IAST,IPST,I),NOCCA,NSSHP,NORBT,NORBT,
C    *                  H2   (IPST,IBST),NSSHP,NSSHB,NORBT,NORBT,
C    *                  H2X  (IAST,IBST),NORBT,NORBT )
         END IF
C
         GOTO 500
C
cmjp occupied-virtual distributions
C Perform two transforms here, which enter different parts
C of h2x, where h2 is an occ-virt distribution.
C H2X(A,D)=SUM(M){ZYMAT(A,M)H2(M,D)-H2(A,M)ZYMAT(M,D)}
C H2X(L,I)=SUM(E){H2(L,E)ZYMAT(E,I)-ZYMAT(L,E)H2(E,I)}
C
C
C Transform occ-virt part of H2
C
 400     CONTINUE
         IF (NOCCP .GT. 0) THEN
            IAST  = IORB(IASYM) + 1 + NOCCA
            IBST  = IORB(IBSYM) + 1 + NOCCB
            IPST  = IORB(IPSYM) + 1
            CALL DGEMM('N','N',NSSHA,NSSHB,NOCCP,1.D0,
     &                 ZYMAT(IAST,IPST,I),NORBT,
     &                 H2(IPST,IBST),NORBT,1.D0,
     &                 H2X(IAST,IBST),NORBT)
C Transform virt-occ part of H2
            CALL DGEMM('N','N',NSSHB,NSSHA,NOCCP,-1.D0,
     &                 H2(IBST,IPST),NORBT,
     &                 ZYMAT(IPST,IAST,I),NORBT,1.D0,
     &                 H2X(IBST,IAST),NORBT)
         END IF
C
C Now do transform which goes into H2X(occ,occ)
C Transform occ-virt part of H2
C
         IF (NOCCB .GT. 0 .AND. NOCCA .GT. 0) THEN
            IAST  = IORB(IASYM) + 1
            IBST  = IORB(IBSYM) + 1
            IPST  = IORB(IPSYM) + 1 + NOCCP
            CALL DGEMM('T','N',NOCCB,NOCCA,NSSHP,1.D0,
     &                 H2(IPST,IBST),NORBT,
     &                 ZYMAT(IPST,IAST,I),NORBT,1.D0,
     &                 H2X(IBST,IAST),NORBT)
C Transform virt-occ part of H2
            CALL DGEMM('T','N',NOCCA,NOCCB,NSSHP,-1.D0,
     &                 TZYMAT(IPST,IAST,I),NORBT,
     &                 H2(IPST,IBST),NORBT,1.D0,
     &                 H2X(IAST,IBST),NORBT)
         END IF
         GOTO 500
C
cmjp virtual-virtual distributions
C transform H2X(L,C)=SUM(E)ZYMAT(L,E)H2(E,C) - SUM(M)H2(L,M)ZYMAT(M,C)
C where h2 is a virt-virt distribution. This contributes to F(ai).
C
C transform H2X(D,L)=SUM(M)ZYMAT(D,M)H2(M,L) - SUM(E)H2(D,E)ZYMAT(E,L)
C where h2 is a virt-virt distribution. This contributes to F(ia).
C
C
C TRANSFORM VIRT-VIRT PART OF H2
C
 600     CONTINUE
         IF (NOCCA .GT. 0) THEN
            IAST  = IORB(IASYM) + 1
            IBST  = IORB(IBSYM) + 1 + NOCCB
            IPST  = IORB(IPSYM) + 1 + NOCCP
C NOW DO F(AI) CONTRIBUTION
            CALL DGEMM('N','N',NSSHB,NOCCA,NSSHP,1.D0,
     &                 H2(IBST,IPST),NORBT,
     &                 TZYMAT(IPST,IAST,I),NORBT,1.D0,
     &                 H2XP(IBST,IAST),NORBT)
C           CALL AMPAB( ZYMAT (IAST,IPST,I),NOCCA,NSSHP,NORBT,NORBT,
C    *                  H2    (IPST,IBST),NSSHP,NSSHB,NORBT,NORBT,
C    *                  H2X   (IAST,IBST),NORBT,NORBT )
C NOW DO F(IA) CONTRIBUTION
            CALL DGEMM('N','N',NSSHB,NOCCA,NSSHP,-1.D0,
     &                 H2(IBST,IPST),NORBT,
     &                 ZYMAT(IPST,IAST,I),NORBT,1.D0,
     &                 H2X(IBST,IAST),NORBT)
         END IF
C TRANSFORM OCC-OCC PART OF H2
         IF (NOCCP .GT. 0 .AND. NOCCB .GT. 0) THEN
            IAST  = IORB(IASYM) + 1 + NOCCA
            IBST  = IORB(IBSYM) + 1
            IPST  = IORB(IPSYM) + 1
C NOW DO F(IA) CONTRIBUTION
            CALL DGEMM('N','N',NSSHA,NOCCB,NOCCP,1.D0,
     &                 ZYMAT(IAST,IPST,I),NORBT,
     &                 H2(IPST,IBST),NORBT,1.D0,
     &                 H2X(IAST,IBST),NORBT)
C NOW DO F(AI) CONTRIBUTION
            CALL DGEMM('N','N',NSSHA,NOCCB,NOCCP,-1.D0,
     &                 TZYMAT(IAST,IPST,I),NORBT,
     &                 H2(IPST,IBST),NORBT,1.D0,
     &                 H2XP(IAST,IBST),NORBT)
C           CALL SMPAB( H2    (IBST,IPST),NOCCB,NOCCP,NORBT,NORBT,
C    *                  ZYMAT (IPST,IAST,I),NOCCP,NSSHA,NORBT,NORBT,
C    *                  H2X   (IBST,IAST),NORBT,NORBT )
C
         END IF
 500     CONTINUE
        ENDDO
C Now call HRPAF to add contributions to the FCONE matrix
C     HRPAF assumes C .ge. D, we know C .le. D therefore
C     we switch ICI and IDI in call of HRPAF
C
         CALL HRPAF(FCONE(1,I),DONEPT,COEMP2,H2X,H2XP,COEUNP,IDI,ICI,
     *                 ICDSYM,ITYPCD,IADR1,WRK,LWRK)
      ENDDO
C
C   End of HRPA
C
 9999 CALL QEXIT('HRPA  ')
      RETURN
      END
C  /* Deck hrpaf */
      SUBROUTINE HRPAF (FCONE,DONEPT,COEMP2,H2X,H2XP,COEUNP,
     *                    ICI,IDI,ICDSYM,ISEL,IADR1,WRK,LWRK)
C
C Copyright by Martin Packer 191193. In the formulas below
C i,j,k.. are occupied and a,b,c... are virtual orbitals.
C Construct the HRPA parts of the transformed Fock matrix.
C ISEL=1 OCC-OCC DISTRIBUTION
C ISEL=4 OCC-VIRT DISTRIBUTION
C ISEL=6 VIRT-VIRT DISTRIBUTION
C
#include <implicit.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      DIMENSION FCONE(NORBT,NORBT), WRK(*), H2XP(NPHMAX)
      DIMENSION COEMP2(*), H2X(NORBT,NORBT)
      DIMENSION COEUNP(NORBT,NORBT), DONEPT(NORBT,NORBT)
C
      IMSYM = MULD2H(ICDSYM,KSYMOP)
      IF (TRPLET) THEN
         MPOFF = LPVMAT
      ELSE
         MPOFF = 0
      END IF
      IF (ISEL .EQ. 1) THEN
C Occ-occ distribution
C Use INDXPH to indicate which occ or virt orbital ICI and IDI
C correspond to.
         ICP = INDXPH(ICI)
         IDP = INDXPH(IDI)
C F(ai) contribution from H2X(occ,virt) stored in H2XP(virt,occ)
C       which we pack into COEUNP
         CALL PHPACK(H2XP,COEUNP,IMSYM,0,-1)
C F(ia) contribution from H2X(virt,occ) packed into H2XP
         CALL PHPACK(H2X ,H2XP,  IMSYM,0,-1)
         DO IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,ISMO(ICI))
            ICSYM = MULD2H(IASYM,ISMO(IDI))
            NOCCA = NOCC(IASYM)
            NSSHA = NSSH(IASYM)
            IAST = IORB(IASYM) + 1 + NOCCA
            IASTM = -INDXPH(IAST)
            IF((KSYMOP.NE.IBSYM).OR.(ICSYM.NE.IMSYM)) GOTO 20
            IF (IFRMP2(IDI).EQ.0) THEN
               DO J = 0,NSSHA-1
                  FCONE(IAST+J,ICI) = FCONE(IAST+J,ICI) +
     *               DDOT(NPHSYM(IMSYM),COEUNP,1,
     *               COEMP2(MPOFF+IADR1(IASTM+J,IDP)+1),1)
                  FCONE(ICI,IAST+J) = FCONE(ICI,IAST+J) +
     *               DDOT(NPHSYM(IMSYM),H2XP  ,1,
     *               COEMP2(MPOFF+IADR1(IASTM+J,IDP)+1),1)
               ENDDO
            ENDIF
 20         IF((KSYMOP.NE.ICSYM).OR.(IBSYM.NE.IMSYM)) GOTO 30
            IF (IFRMP2(ICI).EQ.0) THEN
               IF(ICI .NE. IDI) THEN
                  DO J = 0,NSSHA-1
                    FCONE(IAST+J,IDI) = FCONE(IAST+J,IDI) +
     *                 DDOT(NPHSYM(IMSYM),COEUNP,1,
     *                 COEMP2(MPOFF+IADR1(IASTM+J,ICP)+1),1)
                    FCONE(IDI,IAST+J) = FCONE(IDI,IAST+J) +
     *                 DDOT(NPHSYM(IMSYM),H2XP  ,1,
     *                 COEMP2(MPOFF+IADR1(IASTM+J,ICP)+1),1)
                  ENDDO
               ENDIF
            ENDIF
 30         CONTINUE
          ENDDO
      ENDIF
C NOW FOR VIRT-VIRT DISTRIBUTION
      IF (ISEL .EQ. 6) THEN
         ICP = -INDXPH(ICI)
         IDP = -INDXPH(IDI)
C F(ai) contribution from H2X(occ,virt) stored in H2XP(virt,occ)
         CALL PHPACK(H2XP,COEUNP,IMSYM,0,-1)
C F(ia) contribution from H2X(virt,occ)
         CALL PHPACK(H2X ,H2XP  ,IMSYM,0,-1)
         DO IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,ISMO(ICI))
            ICSYM = MULD2H(IASYM,ISMO(IDI))
            NOCCA = NOCC(IASYM)
            NSSHA = NSSH(IASYM)
            IAST  = IORB(IASYM) + 1
            IASTM = INDXPH(IAST)
            IF((KSYMOP.NE.IBSYM).OR.(ICSYM.NE.IMSYM)) GOTO 60
               DO J = NFRMP2(IASYM),NOCCA-1
                  FCONE(ICI,IAST+J) = FCONE(ICI,IAST+J) +
     *               DDOT(NPHSYM(IMSYM),COEUNP,1,
     *               COEMP2(MPOFF+IADR1(IDP,IASTM+J)+1),1)
                  FCONE(IAST+J,ICI) = FCONE(IAST+J,ICI) +
     *               DDOT(NPHSYM(IMSYM),H2XP  ,1,
     *               COEMP2(MPOFF+IADR1(IDP,IASTM+J)+1),1)
               ENDDO
 60         IF((KSYMOP.NE.ICSYM).OR.(IBSYM.NE.IMSYM)) GOTO 70
               IF (ICI .NE. IDI) THEN
                  DO J = NFRMP2(IASYM),NOCCA-1
                    FCONE(IDI,IAST+J) = FCONE(IDI,IAST+J) +
     *                 DDOT(NPHSYM(IMSYM),COEUNP,1,
     *                 COEMP2(MPOFF+IADR1(ICP,IASTM+J)+1),1)
                    FCONE(IAST+J,IDI) = FCONE(IAST+J,IDI) +
     *                 DDOT(NPHSYM(IMSYM),H2XP  ,1,
     *                 COEMP2(MPOFF+IADR1(ICP,IASTM+J)+1),1)
                  ENDDO
               ENDIF
 70         CONTINUE
         ENDDO
      ENDIF
C NOW DO OCCUPIED-VIRTUAL DISTRIBUTIONS
      IF (ISEL.EQ.4) THEN
C        skip if IDI is frozen in MP2:
         IF (IFRMP2(IDI).NE.0) GO TO 9999
         ICP = -INDXPH(ICI)
         IDP =  INDXPH(IDI)
C        unpack ph block in COEUNP (i.e. COEUNP(i,A) not defined)
         CALL PHPACK(COEUNP,COEMP2(IADR1(ICP,IDP)+1),ICDSYM,-1,-1)
         DO IASYM = 1, NSYM
            IBSYM = MULD2H(IASYM,ICDSYM)
            ICSYM = MULD2H(IASYM,IMSYM)
C CONSTRUCT F(AI) TERMS FROM VIRT-VIRT PART OF H2X
            NOCCB = NOCC(IBSYM)
            IF (NOCCB .GT. 0) THEN
               NSSHA = NSSH(IASYM)
               NSSHC = NSSH(ICSYM)
               IAST = IORB(IASYM) + 1 + NOCC(IASYM)
               IBST = IORB(IBSYM) + 1
               ICST = IORB(ICSYM) + 1 + NOCC(ICSYM)
            CALL DGEMM('N','N',NSSHC,NOCCB,NSSHA,1.D0,
     &                 H2X(ICST,IAST),NORBT,
     &                 COEUNP(IAST,IBST),NORBT,1.D0,
     &                 FCONE(ICST,IBST),NORBT)
C CONSTRUCT F(IA) TERMS FROM VIRT-VIRT PART OF H2X
            CALL DGEMM('T','N',NOCCB,NSSHC,NSSHA,1.D0,
     &                 COEUNP(IAST,IBST),NORBT,
     &                 H2X(IAST,ICST),NORBT,1.D0,
     &                 FCONE(IBST,ICST),NORBT)
            ENDIF
C
C CONSTRUCT F(AI) TERMS FROM OCC-OCC PART OF H2X
C
            NOCCA = NOCC(IASYM)
            NOCCC = NOCC(ICSYM)
            IF (NOCCA .GT. 0 .AND. NOCCC .GT. 0) THEN
               NSSHB = NSSH(IBSYM)
               IAST = IORB(IASYM) + 1
               IBST = IORB(IBSYM) + 1 + NOCCB
               ICST = IORB(ICSYM) + 1
            CALL DGEMM('N','N',NSSHB,NOCCC,NOCCA,1.D0,
     &                 COEUNP(IBST,IAST),NORBT,
     &                 H2X(IAST,ICST),NORBT,1.D0,
     &                 FCONE(IBST,ICST),NORBT)
C CONSTRUCT F(IA) TERMS FROM OCC-OCC PART OF H2X
            CALL DGEMM('N','T',NOCCC,NSSHB,NOCCA,1.D0,
     &                 H2X(ICST,IAST),NORBT,
     &                 COEUNP(IBST,IAST),NORBT,1.D0,
     &                 FCONE(ICST,IBST),NORBT)
            ENDIF
         ENDDO
      ENDIF
 9999 RETURN
C
C   End of HRPAF
C
      END
C  /* Deck rsps2m */
      SUBROUTINE RSPS2M(NSIM,DONEPT,SVEC,ZYMAT,STWO)
C
C Copyright by Martin Packer 251193. In formulas below
C i,j,k.. are occupied and a,b,c... are virtual orbitals.
C
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <infmp2.h>
#include <inforb.h>
#include <infrsp.h>
#include <infvar.h>
#include <wrkrsp.h>
C
      DIMENSION ZYMAT(NORBT,NORBT,*), DONEPT(NORBT,*)
      DIMENSION STWO(NORBT,NORBT,*), SVEC(KZYVAR,*)
C
      PARAMETER(D2 = 2.0D0)
C
      KYCONF = KZCONF + KZVAR
C
      CALL DZERO(STWO,NORBT*NORBT*NSIM)
C FORM S(2).ZYMAT MATRIX AND PLACE IN STWO
      DO I = 1,NSIM
         DO IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,KSYMOP)
            NOCCA = NOCC(IASYM)
            NSSHB = NSSH(IBSYM)
            IF ((NOCCA.NE.0).AND.(NSSHB.NE.0)) THEN
C S(AI) TERMS
               IAST = IORB(IASYM) + 1
               IBST = IORB(IBSYM) + 1 + NOCC(IBSYM)
               CALL DGEMM('N','N',NSSHB,NOCCA,NOCCA,1.D0,
     &                    ZYMAT(IBST,IAST,I),NORBT,
     &                    DONEPT(IAST,IAST),NORBT,1.D0,
     &                    STWO(IBST,IAST,I),NORBT)
               CALL DGEMM('N','N',NSSHB,NOCCA,NSSHB,-1.D0,
     &                    DONEPT(IBST,IBST),NORBT,
     &                    ZYMAT(IBST,IAST,I),NORBT,1.D0,
     &                    STWO(IBST,IAST,I),NORBT)
C S(IA) TERMS
               CALL DGEMM('N','N',NOCCA,NSSHB,NSSHB,1.D0,
     &                    ZYMAT(IAST,IBST,I),NORBT,
     &                    DONEPT(IBST,IBST),NORBT,1.D0,
     &                    STWO(IAST,IBST,I),NORBT)
               CALL DGEMM('N','N',NOCCA,NSSHB,NOCCA,-1.D0,
     &                    DONEPT(IAST,IAST),NORBT,
     &                    ZYMAT(IAST,IBST,I),NORBT,1.D0,
     &                    STWO(IAST,IBST,I),NORBT)
            ENDIF
         ENDDO
      ENDDO
C NOW PACK STWO INTO SVEC USING CODE ADAPTED FROM RSPSOR
      DO IG = 1,KZWOPT
         K     = JWOP(1,IG)
         L     = JWOP(2,IG)
         DO  ISIM = 1 ,NSIM
            SVEC(KZCONF+IG,ISIM) = SVEC(KZCONF+IG,ISIM)
     *              +  STWO(L,K,ISIM)
            SVEC(KYCONF+IG,ISIM) = SVEC(KYCONF+IG,ISIM)
     *              +  STWO(K,L,ISIM)
            SVEC(KZCONF+IG,ISIM) = SVEC(KZCONF+IG,ISIM)
     *              + D2 * ZYMAT(L,K,ISIM)
            SVEC(KYCONF+IG,ISIM) = SVEC(KYCONF+IG,ISIM)
     *              - D2 * ZYMAT(K,L,ISIM)
         ENDDO
      ENDDO
      RETURN
C
C   End of RSPS2M
C
      END
C  /* Deck hrpaa2 */
      SUBROUTINE HRPAA2(COEMP2,H2,A2TEMP,COEUNP,IADR1,ICI,IDI,ICDSYM)
C
C Copyright by Martin Packer 110194. In the formulas below
C i,j,k.. are occupied and a,b,c... are virtual orbitals.
C Add the contribution which symmetrises the A matrix.
C
C Changed to construction of A(2) only. 940916-ekd
C A new routine is then added to symmetrise A.
C
#include <implicit.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      DIMENSION COEMP2(*),H2(NORBT,NORBT),A2TEMP(NORBT,NORBT)
      DIMENSION COEUNP(NORBT,NORBT)
C
C     C occupied, D virtual
C
      IF (IFRMP2(ICI) .NE. 0) RETURN
      ICP = INDXPH(ICI)
      IDP = -INDXPH(IDI)
C     unpack ph block in COEUNP (i.e. COEUNP(i,A) not defined)
      CALL PHPACK(COEUNP,COEMP2(IADR1(IDP,ICP)+1),ICDSYM,-1,-1)
      DO IASYM = 1, NSYM
         IBSYM = MULD2H(IASYM,ICDSYM)
         NSSHA = NSSH(IASYM)
         NOCCB = NOCC(IBSYM)
         IF ((NSSHA.NE.0).AND.(NOCCB.NE.0)) THEN
C CONSTRUCT A2TEMP(IJ) TERMS
            IAST = IORB(IASYM) + 1 + NOCC(IASYM)
            IBST = IORB(IBSYM) + 1
C           we use H2 symmetric to get transpose H2
            CALL DGEMM('T','N',NOCCB,NOCCB,NSSHA,1.D0,
     &                 H2(IAST,IBST),NORBT,
     &                 COEUNP(IAST,IBST),NORBT,1.D0,
     &                 A2TEMP(IBST,IBST),NORBT)
C CONSTRUCT A2TEMP(AB) TERMS
            CALL DGEMM('N','N',NSSHA,NSSHA,NOCCB,1.D0,
     &                 COEUNP(IAST,IBST),NORBT,
     &                 H2(IBST,IAST),NORBT,1.D0,
     &                 A2TEMP(IAST,IAST),NORBT)
         ENDIF
      ENDDO
C
C  END OF HRPAA2
C
      RETURN
      END
C  /* Deck hrpahm */
      SUBROUTINE HRPAHM(FCONE,A2TEMP,ZYMAT,NOSIM)
C
C This routine cancels the nonsymmetric parts of the A matrix
C in SOPPA by using the explicitly constructed A(2) matrix
C which is in A2TEMP (outside called XINDX).
C Adapted by Erik K. Dalskov from HERM_A
C written by Martin J. Packer.
C
#include <implicit.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <infmp2.h>
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      DIMENSION FCONE(NORBT,NORBT,*)
      DIMENSION A2TEMP(NORBT,NORBT), ZYMAT(NORBT,NORBT,*)
C Now cancel the nonsymmetric parts of A(2)
      DO I = 1,NOSIM
         DO IASYM = 1,NSYM
            IBSYM = MULD2H(IASYM,KSYMOP)
            NOCCA = NOCC(IASYM)
            NSSHB = NSSH(IBSYM)
            IF ((NOCCA.NE.0).AND.(NSSHB.NE.0)) THEN
C F(AI) TERMS
               IAST = IORB(IASYM) + 1
               IBST = IORB(IBSYM) + 1 + NOCC(IBSYM)
C
C sum(j)(-b(bj)*A(ai,bj)
               CALL DGEMM('N','N',NSSHB,NOCCA,NOCCA,-1.D0,
     &                    ZYMAT(IBST,IAST,I),NORBT,
     &                    A2TEMP(IAST,IAST),NORBT,1.D0,
     &                    FCONE(IBST,IAST,I),NORBT)
C        +b(bj)*A(bj,ai))
               CALL DGEMM('N','T',NSSHB,NOCCA,NOCCA,1.D0,
     &                    ZYMAT(IBST,IAST,I),NORBT,
     &                    A2TEMP(IAST,IAST),NORBT,1.D0,
     &                    FCONE(IBST,IAST,I),NORBT)
C sum(b)(-b(bj)*A(ai,bj)
               CALL DGEMM('N','N',NSSHB,NOCCA,NSSHB,-1.D0,
     &                    A2TEMP(IBST,IBST),NORBT,
     &                    ZYMAT(IBST,IAST,I),NORBT,1.D0,
     &                    FCONE(IBST,IAST,I),NORBT)
C        +b(bj)*A(bj,ai))
               CALL DGEMM('T','N',NSSHB,NOCCA,NSSHB,1.D0,
     &                    A2TEMP(IBST,IBST),NORBT,
     &                    ZYMAT(IBST,IAST,I),NORBT,1.D0,
     &                    FCONE(IBST,IAST,I),NORBT)
C F(IA) TERMS
C
C sum(b)(-b(jb)*A(jb,ia)
               CALL DGEMM('N','N',NOCCA,NSSHB,NSSHB,-1.D0,
     &                    ZYMAT(IAST,IBST,I),NORBT,
     &                    A2TEMP(IBST,IBST),NORBT,1.D0,
     &                    FCONE(IAST,IBST,I),NORBT)
C       +b(jb)*A(ia,jb))
               CALL DGEMM('N','T',NOCCA,NSSHB,NSSHB,1.D0,
     &                    ZYMAT(IAST,IBST,I),NORBT,
     &                    A2TEMP(IBST,IBST),NORBT,1.D0,
     &                    FCONE(IAST,IBST,I),NORBT)
C sum(j)(-b(jb)*A(jb,ia)
               CALL DGEMM('N','N',NOCCA,NSSHB,NOCCA,-1.D0,
     &                    A2TEMP(IAST,IAST),NORBT,
     &                    ZYMAT(IAST,IBST,I),NORBT,1.D0,
     &                    FCONE(IAST,IBST,I),NORBT)
C       +b(jb)*A(ia,jb))
               CALL DGEMM('T','N',NOCCA,NSSHB,NOCCA,1.D0,
     &                    A2TEMP(IAST,IAST),NORBT,
     &                    ZYMAT(IAST,IBST,I),NORBT,1.D0,
     &                    FCONE(IAST,IBST,I),NORBT)
            ENDIF
         ENDDO
      ENDDO
C
C  END OF HRPAHM
C
      RETURN
      END
C  /* Deck sopdiag */
      SUBROUTINE SOPDIAG(KSYMOP,FC,DIAG,ABSADR,ABTADR,IJSADR,IJTADR,
     &                   IJ1ADR,IJ2ADR,IJ3ADR,WRK,LWRK)
C
C Copyright 11/4-1994 by Erik K. Dalskov
C
C Purpose: Construct D-matrix in SOPPA
C
#include <implicit.h>
      INTEGER ABSADR, ABTADR, A, AA, B, BB, ABSYM, BSYM
C
C
#include <maxorb.h>
#include <maxash.h>
#include <infpri.h>
#include <inforb.h>
#include <infsop.h>
#include <infmp2.h>
#include <infrsp.h>
#include <infind.h>
C
      DIMENSION ABSADR(NORBT,NORBT), ABTADR(NORBT,NORBT), 
     &          IJSADR(NISHT,NISHT,NSYM), IJTADR(NISHT,NISHT,NSYM),
     &          IJ1ADR(NISHT,NISHT,NSYM), IJ2ADR(NISHT,NISHT,NSYM),
     &          IJ3ADR(NISHT,NISHT,NSYM)
      DIMENSION DIAG(*),FC(*),WRK(*)
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KEPSIL,NORBT,WRK,KFREE,LFREE)
      JEPSIL = KEPSIL - 1
      IS = 0
      DO IM=1,NSYM
         IORBM = IORB(IM)
         NORBM = NORB(IM)
         IIORBM = IIORB(IM)
         IF (NORBM .NE. 0) THEN
            DO IJ=1,NORBM
               IL = IIORBM + (IJ*(IJ+1))/2
               WRK((JEPSIL+IS)+IJ) = FC(IL)
            ENDDO
            IS = IS + NORBM
         ENDIF
      ENDDO
C
      DO II=1,NH
         I = IPHORD(II)
         IF (IFRMP2(I) .NE. 0) GO TO 100
         DO JJ=1,II
            J = IPHORD(JJ)
            IF (IFRMP2(J) .NE. 0) GO TO 200
            IF (TRPLET) THEN
               IJ1OFF = IJ1ADR(II,JJ,KSYMOP)
               IJ2OFF = IJ2ADR(II,JJ,KSYMOP)
               IJ3OFF = IJ3ADR(II,JJ,KSYMOP)
            ELSE
               IJSOFF = IJSADR(II,JJ,KSYMOP)
               IJTOFF = IJTADR(II,JJ,KSYMOP)
            END IF
            IJSYM = MULD2H(ISMO(I),ISMO(J))
            ABSYM = MULD2H(IJSYM,KSYMOP)
            ENIJ  = WRK(JEPSIL+I) + WRK(JEPSIL+J)
            DO AA=1,NP
               A = IPHORD(NH+AA)
               ENAIJ = WRK(JEPSIL+A) - ENIJ
               BSYM   = MULD2H(ISMO(A),ABSYM)
               IF (BSYM .GT. ISMO(A)) GO TO 300
               DO BB=1,AA
                  B = IPHORD(NH+BB)
                  IF (ISMO(B) .EQ. BSYM) THEN
                     ENABIJ = ENAIJ + WRK(JEPSIL+B)
                     IF (TRPLET) THEN
                        NABT = ABTADR(AA,BB)
                        IF (I .NE. J) THEN
                           IF (A .NE. B) THEN
                              DIAG(IJ1OFF+NABT)=ENABIJ
                              DIAG(IJ3OFF+NABT)=ENABIJ
                           ENDIF
                           NABS = ABSADR(AA,BB)
                           DIAG(IJ2OFF+NABS)=ENABIJ
                        ELSE
                           IF (A .NE. B) THEN
                              DIAG(IJ3OFF+NABT)=ENABIJ
                           END IF
                        END IF
                     ELSE
                        NABS = ABSADR(AA,BB)
                        DIAG(IJSOFF+NABS) = ENABIJ
                        IF ((I .NE. J) .AND. (A .NE. B)) THEN
                           NABT = ABTADR(AA,BB)
                           DIAG(IJTOFF+NABT) = ENABIJ
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
 300           CONTINUE
            ENDDO
 200        CONTINUE
         ENDDO
 100     CONTINUE
      ENDDO
      CALL MEMREL('SOPDIAG',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C End of SOPDIAG
C
      RETURN
      END
C  /* Deck soppaf */
      SUBROUTINE SOPPAF(FVTD,H2,ZYCVEC,IJ,IB,JBSYM,NCSIM,XINDX,
     *                  WRK,LWRK)
C
C    Purpose: This routine calculates the product of a
C             2p-2h trial vector with a C-matrix for SOPPA.
C
C     Copyright 11/7-1994 by Erik K. Dalskov
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <wrkrsp.h>
#include <infpri.h>
#include <inforb.h>
#include <infmp2.h>
#include <infsop.h>
#include <infind.h>
#include <scbrhf.h>
C
      DIMENSION H2(NORBT,*),FVTD(NORBT,NORBT,*),ZYCVEC(KZCONF,*),WRK(*)
      DIMENSION XINDX(*)
      INTEGER CSYM, ASYM
C
C     skip this distribution if J is frozen in MP2
      IF (IFRMP2(IJ) .NE. 0) GO TO 9999
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KUNP,NP*NH,WRK,KFREE,LFREE)
C
      DO ISIM = 1,NCSIM
         CALL GET2PH(WRK(KUNP),ZYCVEC(1,ISIM),XINDX(KABSAD),
     &               XINDX(KABTAD),XINDX(KIJSAD),XINDX(KIJTAD),
     &               XINDX(KIJ1AD),XINDX(KIJ2AD),XINDX(KIJ3AD),
     &               IJ,IB,JBSYM)
         DO KSYM=1,NSYM
            CSYM  = MULD2H(KSYM,KSYMOP)
            NOCCK = NOCC(KSYM) - NFRRHF(KSYM)
            NSSHC = NSSH(CSYM)
            IF ((NOCCK.NE.0).AND.(NSSHC.NE.0)) THEN
               IKST  = IORB(KSYM) + NFRRHF(KSYM) + 1
               ICST  = IORB(CSYM) + NOCC(CSYM) + 1
C
C   F(c,k) = - SUM(i)[S(c,i)*H2(i,k)]
C             NB! in this summation k may be frozen in MP2
C
               ISYM  = MULD2H(KSYM,JBSYM)
               NOCCI = NOCC(ISYM) - NFRMP2(ISYM)
               IF (NOCCI .GT. 0) THEN
                  IIST  = IORB(ISYM) + NFRMP2(ISYM) + 1
                  ICIST = KUNP + ( IOCC(ISYM) + NFRMP2(ISYM) ) * NP
                  CALL DGEMM('N','N',NSSHC,NOCCK,NOCCI,-1.D0,
     &                       WRK(ICIST),NP,
     &                       H2(IIST,IKST),NORBT,1.D0,
     &                       FVTD(ICST,IKST,ISIM),NORBT)
               END IF
C
C             + SUM(a)[H2(c,a)*S(a,k)]
C             NB! in this summation k must not be frozen in MP2
C
               ASYM  = MULD2H(CSYM,JBSYM)
               NSSHA = NSSH(ASYM)
               NOCCK = NOCC(KSYM) - NFRMP2(KSYM)
               IF (NSSHA .GT. 0 .AND. NOCCK.GT.0) THEN
                  IAST = IORB(ASYM) + NOCC(ASYM) + 1
                  IAKST = KUNP + ( IOCC(KSYM) + NFRMP2(KSYM) ) * NP
                  IKST  = IORB(KSYM) + NFRMP2(KSYM) + 1
                  CALL DGEMM('N','N',NSSHC,NOCCK,NSSHA,1.D0,
     &                       H2(ICST,IAST),NORBT,
     &                       WRK(IAKST),NP,1.D0,
     &                       FVTD(ICST,IKST,ISIM),NORBT)
               END IF
            ENDIF
         ENDDO
      ENDDO
C
C   END OF SOPPAF
C
      CALL MEMREL('SOPPAF',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
 9999 RETURN
      END
C  /* Deck get2ph */
      SUBROUTINE GET2PH(UNPS,ZYCVEC,ABSADR,ABTADR,IJSADR,IJTADR,IJ1ADR,
     &                  IJ2ADR,IJ3ADR,JJ,BB,JBSYM)
C
C     This routine finds the 2p-2h trial vector elements to be used in
C     routine SOPPAF. If a .ne. b and i .ne. j then we have that an
C     unpacked trial matrix are S(1) + ROOT3 * S(2), where the 1 and 2
C     denotes the singlet and triplet coupled elements, respectively.
C
C     Copyright 11/7-1994 by Erik K. Dalskov
C
C     Triplet implemented 22 Feb. 1995 by Thomas Enevoldsen (tec)
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <wrkrsp.h>
#include <infpri.h>
#include <inforb.h>
#include <infsop.h>
#include <infmp2.h>
#include <infind.h>
#include <infrsp.h>
C
      DIMENSION UNPS(NP,*),ZYCVEC(*)
      DIMENSION ABSADR(NORBT,NORBT), ABTADR(NORBT,NORBT)
      DIMENSION IJSADR(NISHT,NISHT,NSYM), IJTADR(NISHT,NISHT,NSYM)
      DIMENSION IJ1ADR(NISHT,NISHT,NSYM), IJ2ADR(NISHT,NISHT,NSYM)
      DIMENSION IJ3ADR(NISHT,NISHT,NSYM)
      INTEGER A,AA,B,BB,ASYM,AOFF,ABSADR,ABTADR
C
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D3 = 3.0D0 , DP5 = 0.5D0 )
C
      LOGICAL FIRST
      SAVE    FIRST,ROOT3H,ROOT2I
      DATA    FIRST /.TRUE./
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         ROOT3H = SQRT(D3) * DP5
         ROOT2I = SQRT(DP5)
      END IF
C
      J = INDXPH(JJ)
      B = -INDXPH(BB)
      IASYM = MULD2H(KSYMOP,JBSYM)
C     Test for triplet or singlet
      IF (.NOT. TRPLET) THEN
         DO ISYM=1,NSYM
            ASYM = MULD2H(ISYM,IASYM)
            IOFF = IORB(ISYM)
            AOFF = IORB(ASYM) + NOCC(ASYM)
            DO II=1,NOCC(ISYM)
               IF (IFRMP2(IOFF + II) .NE. 0) GO TO 100
               I = INDXPH(IOFF + II)
               IJSOFF = IJSADR(I,J,KSYMOP)
               IJTOFF = IJTADR(I,J,KSYMOP)
               IF (I .GE. J) THEN
                  FACIJ = ROOT3H
               ELSE
                  FACIJ = -ROOT3H
               ENDIF
               DO AA=1,NSSH(ASYM)
                  A = -INDXPH(AOFF + AA)
                  NIJABS = IJSOFF + ABSADR(A,B)
                  IF (I .NE. J) THEN
                     IF (A .NE. B) THEN
                        IF (A .GE. B) THEN
                           FAC = FACIJ
                        ELSE
                           FAC = -FACIJ
                        ENDIF
                        NIJABT = IJTOFF + ABTADR(A,B)
                        UNPS(AA,I) = ZYCVEC(NIJABS)*DP5
     *                       + ZYCVEC(NIJABT)*FAC
                     ELSE
                        UNPS(AA,I) = ZYCVEC(NIJABS)*ROOT2I
                     ENDIF
                  ELSE
                     IF (A .NE. B) THEN
                        UNPS(AA,I) = ZYCVEC(NIJABS)*ROOT2I
                     ELSE
                        UNPS(AA,I) = ZYCVEC(NIJABS)
                     ENDIF
                  ENDIF
               ENDDO
 100           CONTINUE
            ENDDO
         ENDDO
      ELSE
C     If TRPLET
         DO ISYM=1,NSYM
            ASYM = MULD2H(ISYM,IASYM)
            IOFF = IORB(ISYM)
            AOFF = IORB(ASYM) + NOCC(ASYM)
            DO II=1,NOCC(ISYM)
               IF (IFRMP2(IOFF + II) .NE. 0) GO TO 110
               I = INDXPH(IOFF + II)
               IJ1OFF = IJ1ADR(I,J,KSYMOP)
               IJ2OFF = IJ2ADR(I,J,KSYMOP)
               IJ3OFF = IJ3ADR(I,J,KSYMOP)
               IF (I .GE. J) THEN
                  TIJSGN = D1
               ELSE
                  TIJSGN = -D1
               ENDIF
               IF (I .EQ. J) THEN
                  C3FAC = ROOT2I
               ELSE
                  C3FAC = DP5
               ENDIF
               DO AA=1,NSSH(ASYM)
                  A = -INDXPH(AOFF + AA)
                  NABS = ABSADR(A,B)
                  NABT = ABTADR(A,B)
                  IF (A .GE. B) THEN
                     TABSGN = D1
                  ELSE
                     TABSGN = -D1
                  ENDIF
                  IF (A .EQ. B) THEN
                     C2FAC =  ROOT2I
                  ELSE
                     C2FAC = DP5
                  ENDIF
C     C(1) CONTRIBUTION
                  IF ((I .NE. J) .AND. (A .NE. B)) THEN
                     UNPS(AA,I) =  TABSGN * TIJSGN
     *                    * ROOT2I * ZYCVEC(IJ1OFF + NABT)
                  ELSE
                     UNPS(AA,I) = D0
                  ENDIF
C     C(2) CONTRIBUTION
                  IF (I .NE. J) THEN
                     UNPS(AA,I) = UNPS(AA,I)
     *                    - C2FAC * TIJSGN * ZYCVEC(IJ2OFF + NABS)
                  ENDIF
C     C(3) CONTRIBUTION
                  IF (A .NE. B) THEN
                     UNPS(AA,I) = UNPS(AA,I)
     *                    - C3FAC * TABSGN * ZYCVEC(IJ3OFF + NABT)
                  ENDIF
               ENDDO
 110        CONTINUE
         ENDDO
         ENDDO
      ENDIF
C
C     END OF GET2PH
C
      RETURN
      END
C  /* Deck soph2x */
      SUBROUTINE SOPH2X(EVECS,ZYMAT,TZYMAT,H2,IJ,IB,JBSYM,NOSIM,
     &                  ABSADR,ABTADR,IJSADR,IJTADR,IJ1ADR,IJ2ADR,
     &                  IJ3ADR,WRK,LWRK)
C
C     This routine calculates ph trial vector times a C-matrix
C     in SOPPA.
C
C     Copyright 11/7-1994 by Erik K. Dalskov
C     Revision 940711/940726 hjaaj
C     Triplet case implementet 21 feb. 1995 by Thomas Enevoldsen (tec)
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <wrkrsp.h>
#include <infpri.h>
#include <inforb.h>
#include <infsop.h>
#include <infmp2.h>
#include <infind.h>
#include <infrsp.h>
C
C
      DIMENSION EVECS(KZYVAR,*),H2(NORBT,*),WRK(*)
      DIMENSION ZYMAT(NORBT,NORBT,*),TZYMAT(NORBT,NORBT,*)
      DIMENSION ABSADR(NORBT,NORBT), ABTADR(NORBT,NORBT)
      DIMENSION IJSADR(NISHT,NISHT,NSYM), IJTADR(NISHT,NISHT,NSYM)
      DIMENSION IJ1ADR(NISHT,NISHT,NSYM), IJ2ADR(NISHT,NISHT,NSYM)
      DIMENSION IJ3ADR(NISHT,NISHT,NSYM)
      INTEGER A,AA,B,ASYM,CSYM, ABSADR, ABTADR
C
      PARAMETER ( D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0 )
C
      LOGICAL FIRST
      SAVE FIRST,ROOT2,ROOT3
      DATA FIRST /.TRUE./
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         ROOT2  = SQRT(D2)
         ROOT3  = SQRT(D3)
      ENDIF
C
C     skip this distribution if J is frozen in MP2
      IF (IFRMP2(IJ) .NE. 0) GO TO 9999
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KB1,NP*NH,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KB2,NP*NH,WRK,KFREE,LFREE)
      LB1 = KB1 - 1
      LB2 = KB2 - 1
C
C
      J = INDXPH(IJ)
      B = -INDXPH(IB)
      IASYM = MULD2H(KSYMOP,JBSYM)
      DO ISYM=1,NSYM
         ASYM = MULD2H(ISYM,IASYM)
         NOCCI = NOCC(ISYM) - NFRMP2(ISYM)
         NSSHA = NSSH(ASYM)
         IF (NOCCI .GT. 0 .AND. NSSHA .GT. 0) THEN
            IIST  = IORB(ISYM) + NFRMP2(ISYM) + 1
            IAST  = IORB(ASYM) + NOCC(ASYM) + 1
C
            KSYM = MULD2H(ISYM,JBSYM)
            NOCCK = NOCC(KSYM)
            CSYM = MULD2H(ASYM,JBSYM)
            NSSHC = NSSH(CSYM)
            IF (NOCCK .GT. 0 .OR. NSSHC .GT. 0) THEN
               DO ISIM = 1,NOSIM
                  IF (NOCCK .EQ. 0) THEN
                     CALL DZERO(WRK(KB1),NSSHA*NOCCI)
                     CALL DZERO(WRK(KB2),NSSHA*NOCCI)
                  ELSE
                     IKST  = IORB(KSYM) + 1
                     CALL DGEMM('N','N',NSSHA,NOCCI,NOCCK,1.D0,
     &                          ZYMAT(IAST,IKST,ISIM),NORBT,
     &                          H2(IKST,IIST),NORBT,0.D0,
     &                          WRK(KB1),NSSHA)
C
                     CALL DGEMM('N','N',NSSHA,NOCCI,NOCCK,1.D0,
     &                          TZYMAT(IAST,IKST,ISIM),NORBT,
     &                          H2(IKST,IIST),NORBT,0.D0,
     &                          WRK(KB2),NSSHA)
                  ENDIF
                  IF (NSSHC .GT. 0) THEN
                     ICST  = IORB(CSYM) + NOCC(CSYM) + 1
                     CALL DGEMM('N','N',NSSHA,NOCCI,NSSHC,-1.D0,
     &                          H2(IAST,ICST),NORBT,
     &                          ZYMAT(ICST,IIST,ISIM),NORBT,1.D0,
     &                          WRK(KB1),NSSHA)
C
                     CALL DGEMM('N','N',NSSHA,NOCCI,NSSHC,-1.D0,
     &                          H2(IAST,ICST),NORBT,
     &                          TZYMAT(ICST,IIST,ISIM),NORBT,1.D0,
     &                          WRK(KB2),NSSHA)
                  ENDIF
                  IF (.NOT. TRPLET) THEN
C     SINGLET CASE
                     DO II=1,NOCCI
                        I = INDXPH(IIST + II - 1)
                        IJSOFF = IJSADR(I,J,KSYMOP)
                        IJTOFF = IJTADR(I,J,KSYMOP)
                        IF (I .GE. J) THEN
                           FACIJ = ROOT3
                        ELSE
                           FACIJ = -ROOT3
                        ENDIF
                        DO AA=1,NSSHA
                           A = -INDXPH(IAST + AA - 1)
                           IF (A .GE. B) THEN
                              FAC = FACIJ
                           ELSE
                              FAC = - FACIJ
                           ENDIF
                           NABS = ABSADR(A,B)
                           NABT = ABTADR(A,B)
                           IF (I .NE. J) THEN
                              IF (A .NE. B) THEN
                                 EVECS(IJSOFF+NABS,ISIM) =
     *                                EVECS(IJSOFF+NABS,ISIM) +
     *                                WRK(LB1+(II-1)*NSSHA+AA)
                                 EVECS(KZVAR+IJSOFF+NABS,ISIM) =
     *                                EVECS(KZVAR+IJSOFF+NABS,ISIM) +
     *                                WRK(LB2+(II-1)*NSSHA+AA)
                                 EVECS(IJTOFF+NABT,ISIM) =
     *                                EVECS(IJTOFF+NABT,ISIM) +
     *                                WRK(LB1+(II-1)*NSSHA+AA) * FAC
                                 EVECS(KZVAR+IJTOFF+NABT,ISIM) =
     *                                EVECS(KZVAR+IJTOFF+NABT,ISIM) +
     *                                WRK(LB2+(II-1)*NSSHA+AA) * FAC
                              ELSE
                                 EVECS(IJSOFF+NABS,ISIM) =
     *                                EVECS(IJSOFF+NABS,ISIM) +
     *                                WRK(LB1+(II-1)*NSSHA+AA) * ROOT2
                                 EVECS(KZVAR+IJSOFF+NABS,ISIM) =
     *                                EVECS(KZVAR+IJSOFF+NABS,ISIM) +
     *                                WRK(LB2+(II-1)*NSSHA+AA) * ROOT2
                              ENDIF
                           ELSE
                              IF (A .NE. B) THEN
                                 EVECS(IJSOFF+NABS,ISIM) =
     *                                EVECS(IJSOFF+NABS,ISIM) +
     *                                WRK(LB1+(II-1)*NSSHA+AA) * ROOT2
                                 EVECS(KZVAR+IJSOFF+NABS,ISIM) =
     *                                EVECS(KZVAR+IJSOFF+NABS,ISIM) +
     *                                WRK(LB2+(II-1)*NSSHA+AA) * ROOT2
                              ELSE
                                 EVECS(IJSOFF+NABS,ISIM) =
     *                                EVECS(IJSOFF+NABS,ISIM) +
     *                                WRK(LB1+(II-1)*NSSHA+AA) * D2
                                 EVECS(KZVAR+IJSOFF+NABS,ISIM) =
     *                                EVECS(KZVAR+IJSOFF+NABS,ISIM) +
     *                                WRK(LB2+(II-1)*NSSHA+AA) * D2
                              ENDIF
                           ENDIF
                        ENDDO
                     ENDDO
                  ELSE
C     Triplet case
                     DO II=1,NOCCI
                        I = INDXPH(IIST + II - 1)
                        IJ1OFF = IJ1ADR(I,J,KSYMOP)
                        IJ2OFF = IJ2ADR(I,J,KSYMOP)
                        IJ3OFF = IJ3ADR(I,J,KSYMOP)
                        IF (I .GE. J) THEN
                           TIJSGN = D1
                        ELSE
                           TIJSGN = -D1
                        ENDIF
                        IF (I .EQ. J) THEN
                           C3FAC = ROOT2
                        ELSE
                           C3FAC = D1
                        ENDIF
                        DO AA=1,NSSHA
                           A = -INDXPH(IAST + AA - 1)
                           IF (A .GE. B) THEN
                              TABSGN = D1
                           ELSE
                              TABSGN = -D1
                           ENDIF
                           IF (A .EQ. B) THEN
                              C2FAC = ROOT2
                           ELSE
                              C2FAC = D1
                           ENDIF
                           NABS = ABSADR(A,B)
                           NABT = ABTADR(A,B)
C     C(1) PART OF VECTOR
                           IF ((I .NE. J) .AND. (A .NE. B)) THEN
                              EVECS(IJ1OFF+NABT,ISIM) =
     *                             EVECS(IJ1OFF+NABT,ISIM) +
     *                             WRK(LB1+(II-1)*NSSHA+AA) * TIJSGN
     *                             * TABSGN * ROOT2
                              EVECS(KZVAR+IJ1OFF+NABT,ISIM) =
     *                             EVECS(KZVAR+IJ1OFF+NABT,ISIM) +
     *                             WRK(LB2+(II-1)*NSSHA+AA) * TIJSGN
     *                             * TABSGN * ROOT2
                           ENDIF
C     C(2) PART OF VECTOR
                           IF (I .NE. J) THEN
                              EVECS(IJ2OFF+NABS,ISIM) =
     *                             EVECS(IJ2OFF+NABS,ISIM) -
     *                             WRK(LB1+(II-1)*NSSHA+AA)*TIJSGN*C2FAC
                              EVECS(KZVAR+IJ2OFF+NABS,ISIM) =
     *                             EVECS(KZVAR+IJ2OFF+NABS,ISIM) -
     *                             WRK(LB2+(II-1)*NSSHA+AA)*TIJSGN*C2FAC
                           ENDIF
C     C(3) PART OF VECTOR
                           IF (A .NE. B) THEN
                              EVECS(IJ3OFF+NABT,ISIM) =
     *                             EVECS(IJ3OFF+NABT,ISIM) -
     *                             WRK(LB1+(II-1)*NSSHA+AA)*TABSGN*C3FAC
                              EVECS(KZVAR+IJ3OFF+NABT,ISIM) =
     *                             EVECS(KZVAR+IJ3OFF+NABT,ISIM) -
     *                             WRK(LB2+(II-1)*NSSHA+AA)*TABSGN*C3FAC
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               END DO
            ENDIF
         ENDIF
      ENDDO
C
      CALL MEMREL('SOPH2X',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C   END OF SOPH2X
C
 9999 RETURN
      END
C  /* Deck sopcli */
      SUBROUTINE SOPCLI(NCSIM,ZYCVEC,EVECS,SOPCL,LWRK)
C
C     This routine calculates D(0)*S, where D(0) is the D-matrix in
C     SOPPA and S is a 2p-2h trial vector
C
C     Copyright 11/7-1994 by Erik K. Dalskov
C
C
#include <implicit.h>
C
C
#include <wrkrsp.h>
#include <inftap.h>
#include <infrsp.h>
C
      DIMENSION EVECS(KZYVAR,*),ZYCVEC(KZCONF,*),SOPCL(*)
C
      IF (KZCONF .GT. LWRK) CALL ERRWRK('SOPCLI',KZCONF,LWRK)
      CALL REWSPL(LURSP4)
#if defined (VAR_SPLITFILES)
      JRSP4 = 0
      CALL READST(LURSP4,LRSP4,SOPCL,JRSP4,.TRUE.)
      IF (LRSP4 .NE. KZCONF) CALL QUIT('Inconsistent record length '//
     &     'on unit LURSP4 in SOPCLI')
#else
      CALL READT(LURSP4,KZCONF,SOPCL)
#endif
C
C
C
      DO ISIM=1,NCSIM
         DO I=1,KZCONF
            EVECS(I,ISIM) = EVECS(I,ISIM) + ZYCVEC(I,ISIM)*SOPCL(I)
         ENDDO
      ENDDO
      RETURN
      END
C  /* Deck sopijab */
      SUBROUTINE SOPIJAB(IG,ISYMV,II,ISYM,JJ,JSYM,AA,ASYM,BB,
     &                   BSYM,STTYP,ABSADR,ABTADR,IJSADR,IJTADR,
     &                   IJ1ADR,IJ2ADR,IJ3ADR)
C
C     This routine finds i,j,a,b-indices corresponding to a specific
C     element IG in a 2p-2h vector
C
C     Copyright 11/7-94 by Erik K. Dalskov
C
#include <implicit.h>
#include <dummy.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
C
      DIMENSION ABSADR(NORBT,NORBT), ABTADR(NORBT,NORBT),
     &          IJSADR(NISHT,NISHT,NSYM), IJTADR(NISHT,NISHT,NSYM),
     &          IJ1ADR(NISHT,NISHT,NSYM), IJ2ADR(NISHT,NISHT,NSYM),
     &          IJ3ADR(NISHT,NISHT,NSYM)
#include <infpri.h>
#include <inforb.h>
#include <infsop.h>
#include <infmp2.h>
#include <infind.h>
#include <infrsp.h>
C
C
      INTEGER A,AA,B,BB,C,CC,D,DD,CSYM,DSYM,ABSYM,DIFF,DISC,ASYM,BSYM,
     &        ABSADR,ABTADR
      CHARACTER*1 STTYP
C
      A = IDUMMY
      B = IDUMMY
      IF (.NOT. TRPLET) THEN
         IF (IG .GT. NS2P2H(ISYMV)) THEN
            STTYP = 'T'
            KOLD = 0
            LOLD = 0
            DO K=2,NH
               IF (IFRMP2(IPHORD(K)) .EQ. 0) THEN
                  DO L=1,K-1
                     IF (IFRMP2(IPHORD(L)) .EQ. 0) THEN
                        DIFF = IG - IJTADR(L,K,ISYMV)
                        IF (DIFF .LE. 0) THEN
                           IF (L .GT. LOLD) THEN
                              I = K
                              J = LOLD
                           ELSE
                              I = KOLD
                              J = LOLD
                           ENDIF
                           GO TO 120
                        ENDIF
                     LOLD = L
                     END IF
                  ENDDO
               KOLD = K
               END IF
            ENDDO
            I = KOLD
            J = LOLD
C
 120        II = IPHORD(I)
            JJ = IPHORD(J)
            IJSYM = MULD2H(ISMO(II),ISMO(JJ))
            ABSYM = MULD2H(IJSYM,ISYMV)
            DISC = IG - IJTADR(I,J,ISYMV)
            DO C =2,NP
               CC = IPHORD(NH+C)
               CSYM = ISMO(CC)
               DO D=1,C-1
                  DD = IPHORD(NH+D)
                  DSYM = ISMO(DD)
                  IF (ABSYM .EQ. MULD2H(DSYM,CSYM)) THEN
                     IF (DISC .EQ. ABTADR(D,C)) THEN
                        A = C
                        B = D
                        GO TO 140
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            WRITE (LUPRI,*) 'SOPIJAB singlet type 2 ERROR CODE 140'
            WRITE (LUPRI,*) 'I,J,II,JJ,IJSYM,ABSYM',I,J,II,JJ,IJSYM,
     &                       ABSYM
            WRITE (LUPRI,*) 'IG,DISC',IG,DISC
            WRITE (LUPRI,*) 'ABTADR array:'
            CALL IWRTMA(ABTADR,NP,NP,NORBT,NORBT)
            WRITE (LUPRI,*) 'IJTADR(*,*,ISYMV) array:'
            CALL IWRTMA(IJTADR(1,1,ISYMV),NH,NH,MAXOCC,MAXOCC)
C           CALL IWRTMA(IMAT,NROW,NCOL,MAXROW,MAXCOL)
            CALL QUIT('SOPIJAB ERROR CODE 140')
 140        CONTINUE
         ELSE
            STTYP = 'S'
            KOLD = 0
            LOLD = 0
            DO K=1,NH
               IF (IFRMP2(IPHORD(K)) .NE. 0) GO TO 333
               DO L=1,K
                  IF (IFRMP2(IPHORD(L)) .NE. 0) GO TO 444
                  DIFF = IG - IJSADR(L,K,ISYMV)
                  IF (DIFF .LE. 0) THEN
                     IF (L .GT. LOLD) THEN
                        I = K
                        J = LOLD
                     ELSE
                        I = KOLD
                        J = LOLD
                     ENDIF
                     GO TO 110
                  ENDIF
                  LOLD = L
 444              CONTINUE
               ENDDO
               KOLD = K
 333           CONTINUE
            ENDDO
            I = KOLD
            J = LOLD
C
 110        II = IPHORD(I)
            JJ = IPHORD(J)
            IJSYM = MULD2H(ISMO(II),ISMO(JJ))
            ABSYM = MULD2H(IJSYM,ISYMV)
            DISC = IG - IJSADR(I,J,ISYMV)
            DO C =1,NP
               CC = IPHORD(NH+C)
               CSYM = ISMO(CC)
               DO D=1,C
                  DD = IPHORD(NH+D)
                  DSYM = ISMO(DD)
                  IF (ABSYM .EQ. MULD2H(DSYM,CSYM)) THEN
                     IF (DISC .EQ. ABSADR(D,C)) THEN
                        A = C
                        B = D
                        GO TO 130
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            WRITE (LUPRI,*) 'SOPIJAB singlet type 1 ERROR CODE 130'
            WRITE (LUPRI,*) 'I,J,II,JJ',I,J,II,JJ
            WRITE (LUPRI,*) 'A,B,IJSYM,ABSYM',A,B,IJSYM,ABSYM
            WRITE (LUPRI,*) 'IG,DISC',IG,DISC
            WRITE (LUPRI,*) 'ABSADR array:'
            CALL IWRTMA(ABSADR,NP,NP,NORBT,NORBT)
            WRITE (LUPRI,*) 'IJSADR(*,*,ISYMV) array:'
            CALL IWRTMA(IJSADR(1,1,ISYMV),NH,NH,MAXOCC,MAXOCC)
C           CALL IWRTMA(IMAT,NROW,NCOL,MAXROW,MAXCOL)
            CALL QUIT('SOPIJAB ERROR CODE 130')
 130        CONTINUE
         END IF
      ELSE
C     TRIPLET CASE
         IF (IG .GT. (N12P2H(ISYMV)+N22P2H(ISYMV))) THEN
            STTYP = '3'
            KOLD = 0
            LOLD = 0
            DO K=1,NH
               IF (IFRMP2(IPHORD(K)) .NE. 0) GO TO 300
               DO L=1,K
                  IF (IFRMP2(IPHORD(L)) .NE. 0) GO TO 400
                  DIFF = IG - IJ3ADR(L,K,ISYMV)
                  IF (DIFF .LE. 0) THEN
                     IF (L.GT.LOLD) THEN
                        I = K
                        J = LOLD
                     ELSE
                        I = KOLD
                        J = LOLD
                     ENDIF
                     GO TO 200
                  ENDIF
                  LOLD = L
 400              CONTINUE
               ENDDO
               KOLD = K
 300           CONTINUE
            ENDDO
            I = KOLD
            J = LOLD
C
 200        II = IPHORD(I)
            JJ = IPHORD(J)
            IJSYM = MULD2H(ISMO(II),ISMO(JJ))
            ABSYM = MULD2H(IJSYM,ISYMV)
            DISC = IG - IJ3ADR(I,J,ISYMV)
            DO C =2,NP
               CC = IPHORD(NH+C)
               CSYM = ISMO(CC)
               DO D=1,C-1
                  DD = IPHORD(NH+D)
                  DSYM = ISMO(DD)
                  IF (ABSYM .EQ. MULD2H(DSYM,CSYM)) THEN
                     IF (DISC .EQ. ABTADR(D,C)) THEN
                        A = C
                        B = D
                        GO TO 100
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            WRITE (LUPRI,*) 'SOPIJAB triplet type 3 ERROR CODE 100'
            WRITE (LUPRI,*) 'I,J,II,JJ,IJSYM,ABSYM',I,J,II,JJ,IJSYM,
     &                       ABSYM
            WRITE (LUPRI,*) 'IG,DISC',IG,DISC
            WRITE (LUPRI,*) 'ABTADR array:'
            CALL IWRTMA(ABTADR,NP,NP,NORBT,NORBT)
            WRITE (LUPRI,*) 'IJ3ADR(*,*,ISYMV) array:'
            CALL IWRTMA(IJ3ADR(1,1,ISYMV),NH,NH,MAXOCC,MAXOCC)
C           CALL IWRTMA(IMAT,NROW,NCOL,MAXROW,MAXCOL)
            CALL QUIT('SOPIJAB ERROR CODE 100')
 100        CONTINUE
         ELSE IF (IG .GT. N12P2H(ISYMV)) THEN
            STTYP = '2'
            KOLD = 0
            LOLD = 0
            DO K=2,NH
               IF (IFRMP2(IPHORD(K)) .EQ. 0) THEN
                  DO L=1,K-1
                     IF (IFRMP2(IPHORD(L)) .EQ. 0) THEN
                        DIFF = IG - IJ2ADR(L,K,ISYMV)
                        IF (DIFF .LE. 0) THEN
                           IF (L .GT. LOLD) THEN
                              I = K
                              J = LOLD
                           ELSE
                              I = KOLD
                              J = LOLD
                           ENDIF
                           GO TO 201
                        ENDIF
                        LOLD = L
                     END IF
                  ENDDO
                  KOLD = K
               END IF
            ENDDO
            I = KOLD
            J = LOLD
C
C
 201        II = IPHORD(I)
            JJ = IPHORD(J)
            IJSYM = MULD2H(ISMO(II),ISMO(JJ))
            ABSYM = MULD2H(IJSYM,ISYMV)
            DISC = IG - IJ2ADR(I,J,ISYMV)
            DO C =1,NP
               CC = IPHORD(NH+C)
               CSYM = ISMO(CC)
               DO D=1,C
                  DD = IPHORD(NH+D)
                  DSYM = ISMO(DD)
                  IF (ABSYM .EQ. MULD2H(DSYM,CSYM)) THEN
                     IF (DISC .EQ. ABSADR(D,C)) THEN
                        A = C
                        B = D
                        GO TO 101
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            WRITE (LUPRI,*) 'SOPIJAB triplet type 2 ERROR CODE 101'
            WRITE (LUPRI,*) 'I,J,II,JJ',I,J,II,JJ
            WRITE (LUPRI,*) 'A,B,IJSYM,ABSYM',A,B,IJSYM,ABSYM
            WRITE (LUPRI,*) 'IG,DISC',IG,DISC
            WRITE (LUPRI,*) 'ABSADR array:'
            CALL IWRTMA(ABSADR,NP,NP,NORBT,NORBT)
            WRITE (LUPRI,*) 'IJ2ADR(*,*,ISYMV) array:'
            CALL IWRTMA(IJ2ADR(1,1,ISYMV),NH,NH,MAXOCC,MAXOCC)
C           CALL IWRTMA(IMAT,NROW,NCOL,MAXROW,MAXCOL)
            CALL QUIT('SOPIJAB ERROR CODE 101')
 101        CONTINUE
         ELSE
            STTYP = '1'
            KOLD = 0
            LOLD = 0
            DO K=2,NH
               IF (IFRMP2(IPHORD(K)) .EQ. 0) THEN
                  DO L=1,K-1
                     IF (IFRMP2(IPHORD(L)) .EQ. 0) THEN
                        DIFF = IG - IJ1ADR(L,K,ISYMV)
                        IF (DIFF .LE. 0) THEN
                           IF (L .GT. LOLD) THEN
                              I = K
                              J = LOLD
                           ELSE
                              I = KOLD
                              J = LOLD
                           ENDIF
                           GO TO 202
                        ENDIF
                        LOLD = L
                     END IF
                  ENDDO
                  KOLD = K
               END IF
            ENDDO
            I = KOLD
            J = LOLD
C
 202        II = IPHORD(I)
            JJ = IPHORD(J)
            IJSYM = MULD2H(ISMO(II),ISMO(JJ))
            ABSYM = MULD2H(IJSYM,ISYMV)
            DISC = IG - IJ1ADR(I,J,ISYMV)
            DO C =2,NP
               CC = IPHORD(NH+C)
               CSYM = ISMO(CC)
               DO D=1,C-1
                  DD = IPHORD(NH+D)
                  DSYM = ISMO(DD)
                  IF (ABSYM .EQ. MULD2H(DSYM,CSYM)) THEN
                     IF (DISC .EQ. ABTADR(D,C)) THEN
                        A = C
                        B = D
                        GO TO 102
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            WRITE (LUPRI,*) 'SOPIJAB triplet type 1 ERROR CODE 102'
            WRITE (LUPRI,*) 'I,J,II,JJ,IJSYM,ABSYM',I,J,II,JJ,IJSYM,
     &                       ABSYM
            WRITE (LUPRI,*) 'IG,DISC',IG,DISC
            WRITE (LUPRI,*) 'ABTADR array:'
            CALL IWRTMA(ABTADR,NP,NP,NORBT,NORBT)
            WRITE (LUPRI,*) 'IJ1ADR(*,*,ISYMV) array:'
            CALL IWRTMA(IJ1ADR(1,1,ISYMV),NH,NH,MAXOCC,MAXOCC)
C           CALL IWRTMA(IMAT,NROW,NCOL,MAXROW,MAXCOL)
            CALL QUIT('SOPIJAB ERROR CODE 102')
 102        CONTINUE
         ENDIF
      ENDIF
      AA = IPHORD(NH+A)
      ASYM = ISMO(AA)
      BB = IPHORD(NH+B)
      BSYM = ISMO(BB)
      ISYM = ISMO(II)
      JSYM = ISMO(JJ)
C
C     END OF SOPIJAB
C
      RETURN
      END
C  /* Deck trpkap */
      SUBROUTINE TRPKAP(COEMP2,TCOEMP,IADR1,IADR2)
C
C     This routine makes the combination:
C     Kt(ai,bj) = 1/3 ( K(ai,bj) + 2 K(aj,bi) )
C     and puts the result in PV(LPVMAT+1)=TCOEMP
C     for use in SOPPA for triplet properties.
C
C     Copyright 941122 Erik K. Dalskov
C
#include <implicit.h>
C
#include <maxorb.h>
#include <maxash.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
#include <inforb.h>
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      INTEGER ASYM,BSYM,AST,AEND,BB1,B1,B1IOFF,AA,A,AIOFF1,BB
C
      DIMENSION COEMP2(*), TCOEMP(*), IADR2(NP,NH)
      PARAMETER ( D3I = 1.0D0/3.0D0 )
C
      CALL DCOPY(LPVMAT,COEMP2(1),1,TCOEMP(1),1)
C
C   call PHADR2 to get the index array IADR2
C
      CALL PHADR2(IADR2,NP,NH)
C
      DO IASYM=1,NSYM
         DO JASYM=1,NSYM
            DO I=1,NH
               II   = IPHORD(I)
               ISYM = ISMO(II)
               ASYM = MULD2H(IASYM,ISYM)
               JSYM = MULD2H(JASYM,ASYM)
               BSYM = MULD2H(IASYM,JSYM)
               AST  = IORB(ASYM) + NOCC(ASYM) + 1
               AEND = IORB(ASYM) + NORB(ASYM)
               JST  = IORB(JSYM) + NFRMP2(JSYM) + 1
               JEND = IORB(JSYM) + NOCC(JSYM)
               BB1  = IORB(BSYM) + NOCC(BSYM) + 1
               B1   = -INDXPH(BB1)
               B1IOFF = IADR2(B1,I) - 1
               NSSHB= NSSH(BSYM)
               DO AA=AST,AEND
                  A = - INDXPH(AA)
                  AIOFF1 = IADR1(A,I) - 1
                  DO JJ=JST,JEND
                     J      = INDXPH(JJ)
                     IAJOFF = IADR1(A,J) + B1IOFF
                     JAIOFF = AIOFF1     + IADR2(B1,J)
                     DO BB=1,NSSHB
                        TCOEMP(JAIOFF + BB) =
     &                  TCOEMP(JAIOFF + BB) + 2 * COEMP2(IAJOFF + BB)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CALL DSCAL(LPVMAT,-D3I,TCOEMP,1)
C
C  END OF TRPKAP
C
      RETURN
      END
C  /* Deck sopdw4 */
      SUBROUTINE SOPDW4(GP,DINVGP,DIAG,EFREQ)
C
C     10-Mar-1995 Hans Joergen Aa. Jensen
C
C     This routine calculates [D(0) +/- EFREQ]**(-1) * GP
C     where D(0) is the D-matrix in SOPPA and
C     GP is the 2p-2h property gradient vector
C
C
#include <implicit.h>
C
C
#include <wrkrsp.h>
#include <inftap.h>
#include <infrsp.h>
C
      DIMENSION GP(KZYVAR), DINVGP(KZYVAR), DIAG(*)
C
      CALL REWSPL(LURSP4)
#if defined (VAR_SPLITFILES)
      JRSP4 = 0
      CALL READST(LURSP4,LRSP4,DIAG,JRSP4,.TRUE.)
      IF (LRSP4 .NE. KZCONF) CALL QUIT('Inconsistent record length '//
     &     'on unit LURSP4 in SOPDW4')
#else
      CALL READT(LURSP4,KZCONF,DIAG)
#endif
C
C     Zero ph part of DINVGP
C
      CALL DZERO(DINVGP(KZCONF+1),KZWOPT)
      CALL DZERO(DINVGP(KZVAR+KZCONF+1),KZWOPT)
C
C     Calculate 2p2h part of DINVGP
C
      DO 400 I = 1,KZCONF
         DINVGP(I) = GP(I) / (DIAG(I) - EFREQ)
         DINVGP(KZVAR+I) = GP(KZVAR + I) / (DIAG(I) + EFREQ)
  400 CONTINUE
      RETURN
      END
C  /* Deck onemp2 */
      SUBROUTINE ONEMP2(PRPMO,DONEPT,ONEACT)
C
C  CALCULATES THE SECOND-ORDER CORRECTION TO AVERAGE VALUE
C
C  Programmed 14/12-1993 by Erik K. Dalskov
C
#include <implicit.h>
C
      DIMENSION PRPMO(NORBT,NORBT), DONEPT(NORBT,NORBT)
C
#include <maxorb.h>
#include <inforb.h>
#include <infdim.h>
#include <infmp2.h>
C
C
         DO 50 IASYM = 1,NSYM
            NORBA = NORB(IASYM)
            IORBA = IORB(IASYM)
            DO 40 I = 1,NORBA
               DO 30 J = 1, NORBA
                 ONEACT = ONEACT + DONEPT(IORBA+J,IORBA+I) *
     *                              PRPMO(IORBA+J,IORBA+I)
30             CONTINUE
40          CONTINUE
50       CONTINUE
C
C  END OF ONEMP2
C
      RETURN
      END
C  /* Deck prpomp */
      SUBROUTINE PRPOMP(PRPMO,PRPSE,DONEPT)
C
C Written by Martin Packer 030394. Calculate the second order
C corrections to the ph transition moments.
C (q|A)_ai = sum_b A_ab D_bi - sum_j D_aj A_ji
C          + sum_j A_aj D_ji - sum_b D_ab A_bi
C (q|A)_ia = sum_j A_ij D_ja - sum_b D_ib P_ba
C          + sum_b A_ib D_ba - sum_j D_ij A_ja
C i,j,k.. are occupied and a,b,c... are virtual orbitals.
C
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <infmp2.h>
#include <infpri.h>
#include <inforb.h>
#include <infrsp.h>
#include <infvar.h>
#include <wrkrsp.h>
C
      DIMENSION DONEPT(NORBT,NORBT), PRPMO(NORBT,NORBT)
      DIMENSION PRPSE(NORBT,NORBT)
C
      CALL DZERO(PRPSE,NORBT*NORBT)
      DO IBSYM = 1,NSYM
         IASYM = MULD2H(IBSYM,KSYMOP)
         NORBA = NORB(IASYM)
         NORBB = NORB(IBSYM)
         NOCCA = NOCC(IASYM)
         NSSHA = NSSH(IASYM)
         NOCCB = NOCC(IBSYM)
         NSSHB = NSSH(IBSYM)
         IF ((NORBA.NE.0).AND.(NORBB.NE.0)) THEN
C PRPSE(AI) a terms of jo.cpr2
            IAST = IORB(IASYM) + 1 + NOCCA
            IBST = IORB(IBSYM) + 1 + NOCCB
            IBST1= IORB(IBSYM) + 1
            CALL DGEMM('N','N',NSSHA,NOCCB,NSSHB,1.D0,
     &                 PRPMO(IAST,IBST),NORBT,
     &                 DONEPT(IBST,IBST1),NORBT,1.D0,
     &                 PRPSE(IAST,IBST1),NORBT)
C PRPSE(IA)  a terms of jo.cpr2
            CALL DGEMM('N','N',NOCCB,NSSHA,NSSHB,-1.D0,
     &                 DONEPT(IBST1,IBST),NORBT,
     &                 PRPMO(IBST,IAST),NORBT,1.D0,
     &                 PRPSE(IBST1,IAST),NORBT)
C PRPSE(AI) a terms of jo.cpr2
            IAST = IORB(IASYM) + 1
            IBST = IORB(IBSYM) + 1 + NOCCB
            IBST1= IORB(IBSYM) + 1
            CALL DGEMM('N','N',NSSHB,NOCCA,NOCCB,-1.D0,
     &                 DONEPT(IBST,IBST1),NORBT,
     &                 PRPMO(IBST1,IAST),NORBT,1.D0,
     &                 PRPSE(IBST,IAST),NORBT)
C PRPSE(IA)  a terms of jo.cpr2
            CALL DGEMM('N','N',NOCCA,NSSHB,NOCCB,1.D0,
     &                 PRPMO(IAST,IBST1),NORBT,
     &                 DONEPT(IBST1,IBST),NORBT,1.D0,
     &                 PRPSE(IAST,IBST),NORBT)
C PRPSE(AI) b terms of jo.cpr2
            IAST = IORB(IASYM) + 1 + NOCCA
            IBST = IORB(IBSYM) + 1
            CALL DGEMM('N','N',NSSHA,NOCCB,NOCCB,1.D0,
     &                 PRPMO(IAST,IBST),NORBT,
     &                 DONEPT(IBST,IBST),NORBT,1.D0,
     &                 PRPSE(IAST,IBST),NORBT)
            CALL DGEMM('N','N',NSSHA,NOCCB,NSSHA,-1.D0,
     &                 DONEPT(IAST,IAST),NORBT,
     &                 PRPMO(IAST,IBST),NORBT,1.D0,
     &                 PRPSE(IAST,IBST),NORBT)
C PRPSE(IA)  b terms of jo.cpr2
            CALL DGEMM('N','N',NOCCB,NSSHA,NSSHA,1.D0,
     &                 PRPMO(IBST,IAST),NORBT,
     &                 DONEPT(IAST,IAST),NORBT,1.D0,
     &                 PRPSE(IBST,IAST),NORBT)
            CALL DGEMM('N','N',NOCCB,NSSHA,NOCCB,-1.D0,
     &                 DONEPT(IBST,IBST),NORBT,
     &                 PRPMO(IBST,IAST),NORBT,1.D0,
     &                 PRPSE(IBST,IAST),NORBT)
         ENDIF
      ENDDO
      RETURN
C
C   End of PRPOMP
C
      END
C  /* Deck prpors */
      SUBROUTINE PRPORS(PRPMO,PRPSE,GP)
C
C WRITTEN 14-FEB 1986
C
C PURPOSE:
C    DISTRIBUTE PROPERTY MO INTEGRALS INTO GP VECTORS IN SOPPA
C
C List of updates
C 21-Jul-1992 Hinne Hettema Average orbital part if necessary.
C 030394 - mjp written for SOPPA case
C          where PRPSE are the second order corrections
C
#include <implicit.h>
C
      DIMENSION PRPMO(NORBT,NORBT),PRPSE(NORBT,NORBT),GP(KZYVAR)
C
C Used from common blocks:
C  INFORB: NORBT
C  INFVAR: JWOP()
C  INFRSP: IPRRSP
C  WRKRSP: KZYVAR
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infvar.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
C -- local constants
C
      PARAMETER ( D2 = 2.0D0 , DP5 = 0.5D0 )
C
C DISTRIBUTE PROPERTY MATRIX IN GP
C
      ROOT2I = SQRT(DP5)
      DO 1400 IG = 1,KZWOPT
         K = JWOP(1,IG)
         L = JWOP(2,IG)
         GP(KZCONF+IG) = GP(KZCONF+IG)
     *      + D2 * PRPMO(L,K) + PRPSE(L,K)
         GP(KZVAR+KZCONF+IG) = GP(KZVAR+KZCONF+IG)
     *      - D2 * PRPMO(K,L) + PRPSE(K,L)
 1400 CONTINUE
C
C      CALL DSCAL(KZWOPT,0.0D0,GP(KZCONF+1),1)
C      CALL DSCAL(KZWOPT,0.0D0,GP(KZVAR+KZCONF+1),1)
C *** Perform averaging
C
      IF (RSPSUP .AND. (KSYMOP .EQ. 1)) THEN
         CALL RSPAVE(GP(KZCONF+1),KZVAR,2)
      END IF
C
      IF (IPRRSP.GT.20) THEN
         WRITE(LUPRI,'(/A)')
     &   ' (PRPORS) Z and Y property GP vector in SOPPA approximation'
         CALL OUTPUT(GP,1,KZVAR,1,2,KZVAR,2,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck prpcmp */
      SUBROUTINE PRPCMP(PRPMO,COEMP2,GP,ABSADR,ABTADR,IJSADR,IJTADR,
     &                  IADR1,WRK,KFREE,LFREE)
C
C  Copyright 8/4-94 by Erik K. Dalskov
C
C  Purpose: Construct 2p-2h corrections to transition moments
C           in SOPPA
C
#include <implicit.h>
C
C  Used from common blocks:
C  INFRSP : IPRRSP
C  WRKRSP : KSYMOP,KZYVAR
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infpri.h>
#include <inforb.h>
#include <infsop.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      INTEGER BSYM,CSYM,A,AA,B,BB,ABSADR,ABTADR
      PARAMETER ( D1 = 1.0D0 , D3 = 3.0D0 , D4 = 4.0D0 )
      PARAMETER ( D2 = 2.0D0 , D8 = 8.0D0 )
C
      DIMENSION PRPMO(NORBT,*), COEMP2(*), GP(KZYVAR)
      DIMENSION ABSADR(NORBT,NORBT), ABTADR(NORBT,NORBT),
     &          IJSADR(NISHT,NISHT,NSYM), IJTADR(NISHT,NISHT,NSYM)
      DIMENSION WRK(*)
C
      ROOT3I = D1 / SQRT(D3)
      ROOT2  = SQRT(D2)
      ROOT8  = SQRT(D8)
C
C Work space allocation
C
      KFRSAV = KFREE
      CALL MEMGET('REAL',KMPUNP,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRESZ,NP*NH,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRESY,NP*NH,WRK,KFREE,LFREE)
      LMPUNP = KMPUNP - 1
      LRESZ  = KRESZ - 1
      LRESY  = KRESY - 1
C
      DO II=1,NH
         I = IPHORD(II)
         IF (IFRMP2(I) .NE. 0) GO TO 9999
         DO AA=1,NP
            A = IPHORD(NH+AA)
            IASYM = MULD2H(ISMO(A),ISMO(I))
            JBSYM = MULD2H(IASYM,KSYMOP)
            CALL PHPACK(WRK(KMPUNP),COEMP2(IADR1(AA,II)+1),IASYM,-1,-1)
            DO JSYM=1,NSYM
               BSYM  = MULD2H(JSYM,JBSYM)
               CSYM  = MULD2H(JSYM,IASYM)
               KSYM  = MULD2H(BSYM,IASYM)
               NOCCJ = NOCC(JSYM) - NFRMP2(JSYM)
               NOCCK = NOCC(KSYM)
               NOCCB = NOCC(BSYM)
               NOCCC = NOCC(CSYM)
               NSSHB = NSSH(BSYM)
               NSSHC = NSSH(CSYM)
               IF ((NOCCJ .NE. 0) .AND. (NSSHB .NE. 0)) THEN
                  IBST = IORB(BSYM) + NOCCB + 1
                  ICST = IORB(CSYM) + NOCCC + 1
                  IJST = IORB(JSYM) + NFRMP2(JSYM) + 1
                  IKST = IORB(KSYM) + 1
                  JCST = ( IORB(JSYM) + NFRMP2(JSYM) ) * NORBT + ICST
                  KBST = IORB(KSYM) * NORBT + IBST
                  IF (NSSHC .EQ. 0) THEN
                     CALL DZERO(WRK(KRESZ),NSSHB*NOCCJ)
                  ELSE
                     CALL DGEMM('N','N',NSSHB,NOCCJ,NSSHC,1.D0,
     &                          PRPMO(IBST,ICST),NORBT,
     &                          WRK(LMPUNP+JCST),NORBT,0.D0,
     &                          WRK(KRESZ),NSSHB)
                  END IF
                  CALL DGEMM('N','N',NSSHB,NOCCJ,NOCCK,-1.D0,
     &                       WRK(LMPUNP+KBST),NORBT,
     &                       PRPMO(IKST,IJST),NORBT,1.D0,
     &                       WRK(KRESZ),NSSHB)
                  IF (NOCCK .EQ. 0) THEN
                     CALL DZERO(WRK(KRESY),NSSHB*NOCCJ)
                  ELSE
                     CALL DGEMM('N','T',NSSHB,NOCCJ,NOCCK,1.D0,
     &                          WRK(LMPUNP+KBST),NORBT,
     &                          PRPMO(IJST,IKST),NORBT,0.D0,
     &                          WRK(KRESY),NSSHB)
                  END IF
                  CALL DGEMM('T','N',NSSHB,NOCCJ,NSSHC,-1.D0,
     &                       PRPMO(ICST,IBST),NORBT,
     &                       WRK(LMPUNP+JCST),NORBT,1.D0,
     &                       WRK(KRESY),NSSHB)
                  DO J=1,NOCCJ
                     JJ = INDXPH(IJST-1+J)
                     IJSOFF = IJSADR(II,JJ,KSYMOP)
                     IJTOFF = IJTADR(II,JJ,KSYMOP)
                     IF (II .GE. JJ) THEN
                        TIJSGN = ROOT3I
                     ELSE
                        TIJSGN = -ROOT3I
                     ENDIF
                     DO B = 1,NSSHB
                        BB = - INDXPH(IBST-1+B)
                        NABS = ABSADR(AA,BB)
                        NABT = ABTADR(AA,BB)
                        IF (II .EQ. JJ) THEN
                           IF (AA .EQ. BB) THEN
                              FAC = D2
                           ELSE
                              FAC = ROOT2
                           ENDIF
                        ELSE
                           IF (AA .EQ. BB) THEN
                              FAC = ROOT2
                           ELSE
                              FAC = D1
                           ENDIF
                        ENDIF
                        GP(IJSOFF + NABS) =
     *                     GP(IJSOFF + NABS)
     *                     - WRK(LRESZ+(J-1)*NSSHB+B) * FAC
                        GP(KZVAR + IJSOFF + NABS) =
     *                     GP(KZVAR + IJSOFF + NABS)
     *                     - WRK(LRESY+(J-1)*NSSHB+B) * FAC
                        IF ((II .NE. JJ) .AND. (AA .NE. BB)) THEN
                           IF (AA .GT .BB) THEN
                              TSIGN = -TIJSGN
                           ELSE
                              TSIGN = TIJSGN
                           ENDIF
                           GP(IJTOFF + NABT) = GP(IJTOFF + NABT)
     *                        + TSIGN * WRK(LRESZ+(J-1)*NSSHB+B)
                           GP(KZVAR + IJTOFF + NABT) =
     *                        GP(KZVAR + IJTOFF + NABT)
     *                        + TSIGN * WRK(LRESY+(J-1)*NSSHB+B)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
 9999    CONTINUE
      ENDDO
      IF ( IPRRSP.GT.110) THEN
         WRITE(LUPRI,'(/A)')
     *   ' 2p-2h transition moments for SOPPA (Z-part and Y-part)'
         CALL OUTPUT(GP,1,KZCONF,1,2,KZVAR,2,1,LUPRI)
      END IF
C
C  End of PRPCMP
C
      CALL MEMREL('PRPCMP',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C  /* Deck prpcmpt */
      SUBROUTINE PRPTCMP(PRPMO,COEMP2,GP,ABSADR,ABTADR,IJ1ADR,IJ2ADR,
     &                   IJ3ADR,IADR1,WRK,KFREE,LFREE)
C
C Written 14-Feb-1995 by Thomas Enevoldsen (tec)
C
C Purpose: Construct 2p-2h contributions to transition moments
C in triplet SOPPA
C
C
C     Notice COEMP2 = Kappa / 2
C
C
#include <implicit.h>
C
C Used from common blocks:
C INFRSP : IPRRSP
C WRKRSP : KSYMOP,KZYVAR
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infpri.h>
#include <inforb.h>
#include <infsop.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
#include <infind.h>
#include <infrsp.h>
#include <wrkrsp.h>
C
      INTEGER BSYM,CSYM,A,AA,B,BB,ABSADR,ABTADR
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 , D3 = 3.0D0 , D6 = 6.0D0 )
C
      DIMENSION PRPMO(NORBT,*), COEMP2(*), GP(KZYVAR)
      DIMENSION ABSADR(NORBT,NORBT), ABTADR(NORBT,NORBT),
     &          IJ1ADR(NISHT,NISHT,NSYM), IJ2ADR(NISHT,NISHT,NSYM),
     &          IJ3ADR(NISHT,NISHT,NSYM)
      DIMENSION WRK(*)
C
      ROOT2  = SQRT(D2)
      D3I = D1 / D3
      T1FAC =  ROOT2 * D3I
C
C Work space allocation
C
      KFRSAV = KFREE
      CALL MEMGET('REAL',KMPUNP,NORBT*NORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRESZ1,NP*NH,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRESZ2,NP*NH,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRESY1,NP*NH,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KRESY2,NP*NH,WRK,KFREE,LFREE)
      LMPUNP = KMPUNP - 1
      LRESZ1  = KRESZ1 - 1
      LRESY1  = KRESY1 - 1
      LRESZ2  = KRESZ2 - 1
      LRESY2  = KRESY2 - 1
C
      DO II=1,NH
         I = IPHORD(II)
         IF (IFRMP2(I) .NE. 0) GO TO 9999
         DO AA=1,NP
            A = IPHORD(NH+AA)
            IASYM = MULD2H(ISMO(A),ISMO(I))
            JBSYM = MULD2H(IASYM,KSYMOP)
            CALL PHPACK(WRK(KMPUNP),COEMP2(IADR1(AA,II)+1),IASYM,-1,-1)
            DO JSYM=1,NSYM
               BSYM  = MULD2H(JSYM,JBSYM)
               CSYM  = MULD2H(JSYM,IASYM)
               KSYM  = MULD2H(BSYM,IASYM)
               NOCCJ = NOCC(JSYM) - NFRMP2(JSYM)
               NOCCK = NOCC(KSYM)
               NOCCB = NOCC(BSYM)
               NOCCC = NOCC(CSYM)
               NSSHB = NSSH(BSYM)
               NSSHC = NSSH(CSYM)
               IF ((NOCCJ .NE. 0) .AND. (NSSHB .NE. 0)) THEN
                  IBST = IORB(BSYM) + NOCCB + 1
                  ICST = IORB(CSYM) + NOCCC + 1
                  IJST = IORB(JSYM) + NFRMP2(JSYM) + 1
                  IKST = IORB(KSYM) + 1
                  JCST = ( IORB(JSYM) + NFRMP2(JSYM) ) * NORBT + ICST
                  KBST = IORB(KSYM) * NORBT + IBST
                  IF (NSSHC .EQ. 0) THEN
                     CALL DZERO(WRK(KRESZ1),NSSHB*NOCCJ)
                     CALL DZERO(WRK(KRESY2),NSSHB*NOCCJ)
                  ELSE
                     CALL DGEMM('N','N',NSSHB,NOCCJ,NSSHC,1.D0,
     &                          PRPMO(IBST,ICST),NORBT,
     &                          WRK(LMPUNP+JCST),NORBT,0.D0,
     &                          WRK(KRESZ1),NSSHB)
                     CALL DGEMM('T','N',NSSHB,NOCCJ,NSSHC,1.D0,
     &                          PRPMO(ICST,IBST),NORBT,
     &                          WRK(LMPUNP+JCST),NORBT,0.D0,
     &                          WRK(KRESY2),NSSHB)
                  END IF
                  IF (NOCCK .EQ. 0) THEN
                     CALL DZERO(WRK(KRESZ2),NSSHB*NOCCJ)
                     CALL DZERO(WRK(KRESY1),NSSHB*NOCCJ)
                  ELSE
                     CALL DGEMM('N','N',NSSHB,NOCCJ,NOCCK,1.D0,
     &                          WRK(LMPUNP+KBST),NORBT,
     &                          PRPMO(IKST,IJST),NORBT,0.D0,
     &                          WRK(KRESZ2),NSSHB)
                     CALL DGEMM('N','T',NSSHB,NOCCJ,NOCCK,1.D0,
     &                          WRK(LMPUNP+KBST),NORBT,
     &                          PRPMO(IJST,IKST),NORBT,0.D0,
     &                          WRK(KRESY1),NSSHB)
                  ENDIF
                  DO J=1,NOCCJ
                     JJ = INDXPH(IJST-1+J)
                     IJ1OFF = IJ1ADR(II,JJ,KSYMOP)
                     IJ2OFF = IJ2ADR(II,JJ,KSYMOP)
                     IJ3OFF = IJ3ADR(II,JJ,KSYMOP)
                     IF (II .GT. JJ) THEN
                        TIJSGN = D1
                     ELSE
                        TIJSGN = -D1
                     ENDIF
                     IF (II .EQ. JJ) THEN
                        T3FAC = ROOT2
                     ELSE
                        T3FAC = D1
                     ENDIF
                     DO B = 1,NSSHB
                        BB = - INDXPH(IBST-1+B)
                        NABS = ABSADR(AA,BB)
                        NABT = ABTADR(AA,BB)
                        IF (AA .GT. BB) THEN
                           TABSGN = D1
                        ELSE
                           TABSGN = -D1
                        ENDIF
                        IF (AA .EQ. BB) THEN
                           T2FAC = ROOT2
                        ELSE
                           T2FAC = D1
                        ENDIF
C T(1) part of GP-vector
                        IF ((II .NE. JJ) .AND. (AA .NE. BB)) THEN
                           GP(IJ1OFF + NABT) = GP(IJ1OFF + NABT)
     *                          - TABSGN * TIJSGN * T1FAC *
     *                     (WRK(LRESZ1+(J-1)*NSSHB+B)-WRK(LRESZ2+(J-1)
     *                          *NSSHB+B))
                           GP(KZVAR + IJ1OFF + NABT) =
     *                          GP(KZVAR + IJ1OFF + NABT)
     *                          - TABSGN * TIJSGN * T1FAC *
     *                     (WRK(LRESY1+(J-1)*NSSHB+B)-WRK(LRESY2+(J-1)
     *                          *NSSHB+B))
                        ENDIF
C T(2) part of GP-vector
                        IF (II .NE. JJ)  THEN
                           GP(IJ2OFF + NABS) = GP(IJ2OFF + NABS)
     *                          + TIJSGN * T2FAC *
     *                     (WRK(LRESZ2+(J-1)*NSSHB+B) - WRK(LRESZ1
     *                          +(J-1)*NSSHB+B) * D3I)
                           GP(KZVAR + IJ2OFF + NABS) = GP(KZVAR + IJ2OFF
     *                          + NABS) + TIJSGN * T2FAC *
     *                     (WRK(LRESY2+(J-1)*NSSHB+B) * D3I - WRK(LRESY1
     *                          +(J-1)*NSSHB+B) )
                        ENDIF
C T(3) part of GP-vector
                        IF (AA .NE. BB)  THEN
                           GP(IJ3OFF + NABT) = GP(IJ3OFF + NABT)
     *                          + TABSGN * T3FAC *
     *                     (WRK(LRESZ2+(J-1)*NSSHB+B) * D3I - WRK(LRESZ1
     *                          +(J-1)*NSSHB+B) )
                           GP(KZVAR + IJ3OFF + NABT) = GP(KZVAR + IJ3OFF
     *                          + NABT) + TABSGN * T3FAC *
     *                     (WRK(LRESY2+(J-1)*NSSHB+B) - WRK(LRESY1
     *                          +(J-1)*NSSHB+B) * D3I )
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
 9999    CONTINUE
      ENDDO
      IF ( IPRRSP.GT.110) THEN
         WRITE(LUPRI,'(/A)')
     *        ' 2p-2h transition moments for SOPPA (Z-part and Y-part)'
         CALL OUTPUT(GP(1),1,KZCONF,1,2,KZVAR,2,1,LUPRI)
      END IF
C
C End of PRPTCMP
C
      CALL MEMREL('PRPTCMP',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C  /* Deck rspact */
      SUBROUTINE RSPACT
C
C 16-Nov-1993 mjp+hjaaj
C
#include <implicit.h>
C
C Used from common blocks:
C  INFRSP: SOPPA,LPVMAT,NACTT,NACT(8),IACT(8),JACT(8)
C  INFORB: NSYM,*ASH*,*ORB*,NVIRT,NRHFT
C  INFMP2: NPHTOT(1)
C
#include <maxorb.h>
#include <infrsp.h>
#include <inforb.h>
#include <infmp2.h>
C
      IF (SOPPA) THEN
         LPVMAT = NPHTOT(1)
         NACTT = NORBT
         DO ISYM = 1,NSYM
            NACT(ISYM) = NORB(ISYM)
            IACT(ISYM) = IORB(ISYM)
            JACT(ISYM) = IORB(ISYM)
         END DO
      ELSE
         LPVMAT = N2ASHX*N2ASHX
         NACTT = NASHT
         DO ISYM = 1,NSYM
            NACT(ISYM) = NASH(ISYM)
            IACT(ISYM) = IASH(ISYM)
            JACT(ISYM) = IORB(ISYM) + NISH(ISYM)
         END DO
      END IF
      RETURN
      END
C  /* Deck sopij */
      SUBROUTINE SOPIJ(ISYMV,IOFFY,IUNIT,RSPVEC,ABSADR,ABTADR,IJSADR,
     &                 IJTADR,IJ1ADR,IJ2ADR,IJ3ADR,THPSOP)
C
C     This routine sums up all excitation contributions to a specific
C     I,J pair in the 2p2h solution vector.
C
C     Copyright 23/2-96 by Erik K. Dalskov and Thomas Enevoldsen
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <infpri.h>
#include <inforb.h>
#include <infsop.h>
#include <infmp2.h>
#include <infind.h>
#include <infrsp.h>
C
      INTEGER ABSYM,BASYM,ASYM,BSYM,A,B,ABSOFF,ABTOFF,ABSADR,ABTADR
      DIMENSION RSPVEC(*), ABSADR(NORBT,NORBT), ABTADR(NORBT,NORBT),
     &          IJSADR(NISHT,NISHT,NSYM), IJTADR(NISHT,NISHT,NSYM),
     &          IJ1ADR(NISHT,NISHT,NSYM), IJ2ADR(NISHT,NISHT,NSYM),
     &          IJ3ADR(NISHT,NISHT,NSYM)
      PARAMETER ( ZERO = 0.0D0 )
C
      IKKEPR = 0
      EXCITOT = ZERO
      THPR = THPSOP * THPSOP
      IF (.NOT.TRPLET) THEN
         WRITE(IUNIT,1000)
         DO I=1,NH
            II = IPHORD(I)
            IF (IFRMP2(II) .EQ. 0) THEN
               ISYM=ISMO(II)
               DO J=1,I
                  JJ = IPHORD(J)
                  IF (IFRMP2(JJ).EQ.0) THEN
                     JSYM=ISMO(JJ)
                     IJSYM=MULD2H(ISYM,JSYM)
                     ABSYM=MULD2H(ISYMV,IJSYM)
                     IJSOFF=IJSADR(I,J,ISYMV)
                     IJTOFF=IJTADR(I,J,ISYMV)
                     EXCISZ = ZERO
                     EXCITZ = ZERO
                     EXCISY = ZERO
                     EXCITY = ZERO
                     DO A=1,NP
                        ASYM=ISMO(IPHORD(NH+A))
                        DO B=1,A
                           BSYM=ISMO(IPHORD(NH+B))
                           BASYM=MULD2H(BSYM,ASYM)
                           IF (ABSYM .EQ. BASYM) THEN
                              ABSOFF = ABSADR(A,B)
                              ABTOFF = ABTADR(A,B)
                              EXCISZ = EXCISZ + RSPVEC(IJSOFF+ABSOFF)
     &                                        * RSPVEC(IJSOFF+ABSOFF)
                              EXCISY = EXCISY
     &                                  + RSPVEC(IOFFY+IJSOFF+ABSOFF) 
     &                                  * RSPVEC(IOFFY+IJSOFF+ABSOFF)
                              IF (I.NE.J .AND. A.NE.B) THEN
                                 EXCITZ = EXCITZ + RSPVEC(IJTOFF+ABTOFF)
     &                                           * RSPVEC(IJTOFF+ABTOFF)
                                 EXCITY = EXCITY
     &                                     + RSPVEC(IOFFY+IJTOFF+ABTOFF)
     &                                     * RSPVEC(IOFFY+IJTOFF+ABTOFF)
                              ENDIF
                           END IF
                        ENDDO
                     ENDDO
                     EXCIMAX = MAX(EXCISZ,EXCISY,EXCITZ,EXCITY)
                     IF (EXCIMAX .GT. THPR) THEN
                        IKKEPR = IKKEPR + 1
                        EXCISZ = SQRT(EXCISZ)
                        EXCISY = SQRT(EXCISY)
                        EXCITZ = SQRT(EXCITZ)
                        EXCITY = SQRT(EXCITY)
                        WRITE(IUNIT,1001) II,ISYM,JJ,JSYM,
     &                                    EXCISZ,EXCISY,EXCITZ,EXCITY
                     ELSE
                        EXCITOT = EXCITOT +
     &                            EXCISZ + EXCISY + EXCITZ + EXCITY
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ELSE
         WRITE(IUNIT,1100)
         DO I=1,NH
            II = IPHORD(I)
            IF (IFRMP2(II) .EQ. 0) THEN
               ISYM=ISMO(II)
               DO J=1,I
                  JJ = IPHORD(J)
                  IF (IFRMP2(JJ) .EQ. 0) THEN
                     JSYM=ISMO(JJ)
                     IJSYM=MULD2H(ISYM,JSYM)
                     ABSYM=MULD2H(ISYMV,IJSYM)
                     IJ1OFF=IJ1ADR(I,J,ISYMV)
                     IJ2OFF=IJ2ADR(I,J,ISYMV)
                     IJ3OFF=IJ3ADR(I,J,ISYMV)
                     EXCI1Z = ZERO
                     EXCI2Z = ZERO
                     EXCI3Z = ZERO
                     EXCI1Y = ZERO
                     EXCI2Y = ZERO
                     EXCI3Y = ZERO
                     DO A=1,NP
                        ASYM=ISMO(IPHORD(NH+A))
                        DO B=1,A
                           BSYM=ISMO(IPHORD(NH+B))
                           BASYM=MULD2H(BSYM,ASYM)
                           IF (ABSYM .EQ. BASYM) THEN
                              ABSOFF = ABSADR(A,B)
                              ABTOFF = ABTADR(A,B)
                              IF (I.NE.J) THEN
                                 IF(A.NE.B) THEN
                                    EXCI1Z = EXCI1Z +
     &                                       RSPVEC(IJ1OFF+ABTOFF) *
     &                                       RSPVEC(IJ1OFF+ABTOFF)
                                    EXCI2Z = EXCI2Z +
     &                                       RSPVEC(IJ2OFF+ABSOFF) *
     &                                       RSPVEC(IJ2OFF+ABSOFF)
                                    EXCI3Z = EXCI3Z +
     &                                       RSPVEC(IJ3OFF+ABTOFF) *
     &                                       RSPVEC(IJ3OFF+ABTOFF)
                                    EXCI1Y = EXCI1Y +
     &                                 RSPVEC(IOFFY+IJ1OFF+ABTOFF) *
     &                                 RSPVEC(IOFFY+IJ1OFF+ABTOFF)
                                    EXCI2Y = EXCI2Y +
     &                                 RSPVEC(IOFFY+IJ2OFF+ABSOFF) *
     &                                 RSPVEC(IOFFY+IJ2OFF+ABSOFF)
                                    EXCI3Y = EXCI3Y +
     &                                 RSPVEC(IOFFY+IJ3OFF+ABTOFF) *
     &                                 RSPVEC(IOFFY+IJ3OFF+ABTOFF)
                                 ELSE
                                    EXCI2Z = EXCI2Z +
     &                                       RSPVEC(IJ2OFF+ABSOFF) *
     &                                       RSPVEC(IJ2OFF+ABSOFF)
                                    EXCI2Y = EXCI2Y +
     &                                 RSPVEC(IOFFY+IJ2OFF+ABSOFF) *
     &                                 RSPVEC(IOFFY+IJ2OFF+ABSOFF)
                                 ENDIF
                              ELSE
                                 IF (A.NE.B) THEN
                                    EXCI3Z = EXCI3Z +
     &                                       RSPVEC(IJ3OFF+ABTOFF) *
     &                                       RSPVEC(IJ3OFF+ABTOFF)
                                    EXCI3Y = EXCI3Y +
     &                                 RSPVEC(IOFFY+IJ3OFF+ABTOFF) *
     &                                 RSPVEC(IOFFY+IJ3OFF+ABTOFF)
                                 ENDIF
                              ENDIF
                           END IF
                        ENDDO
                     ENDDO
                     EXCIMAX = MAX(EXCI1Z,EXCI1Y,EXCI2Z,
     &                             EXCI2Y,EXCI3Z,EXCI3Y)
                     IF (EXCIMAX .GT. THPR) THEN
                        IKKEPR = IKKEPR + 1
                        EXCI1Z = SQRT(EXCI1Z)
                        EXCI2Z = SQRT(EXCI2Z)
                        EXCI3Z = SQRT(EXCI3Z)
                        EXCI1Y = SQRT(EXCI1Y)
                        EXCI2Y = SQRT(EXCI2Y)
                        EXCI3Y = SQRT(EXCI3Y)
                        WRITE(IUNIT,1101) II,ISYM,JJ,JSYM,
     &                         EXCI1Z,EXCI1Y,EXCI2Z,EXCI2Y,EXCI3Z,EXCI3Y
                     ELSE
                        EXCITOT = EXCITOT + EXCI1Z + EXCI2Z + EXCI3Z
     &                                    + EXCI1Y + EXCI2Y + EXCI3Y
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      EXCITOT = SQRT(EXCITOT)
      IKKEPR = NH*(NH+1)/2 - IKKEPR
      WRITE(IUNIT,1500) IKKEPR,THPSOP,EXCITOT
C
C  End of SOPIJ
C
      RETURN
 1000 FORMAT(/'   i      j     S  Z oper.   Y oper.',
     &            '   T  Z oper.   Y oper.',
     &       /'  ---    ---   --- --------  --------',
     &        ' --- --------  --------')
 1001 FORMAT(I4,'(',I1,')',I4,'(',I1,')',3X,2F10.6,3X,2F10.6)
 1100 FORMAT(/'   i      j     1  Z oper.  Y oper.',
     &            '  2  Z oper.  Y oper.',
     &            '  3  Z oper.  Y oper.',
     &       /'  ---    ---   --- -------  -------',
     &        ' --- -------  -------',
     &        ' --- -------  -------')
 1101 FORMAT(I4,'(',I1,')',I4,'(',I1,')',3X,2F9.5,3X,2F9.5,3X,2F9.5)
 1500 FORMAT(/5X,I5,' elements with norm less than',1P,D9.2,
     &              ' not printed',
     &       //' The total norm of the excitation out of these',
     &         ' elements is',1P,D11.4)
      END
