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
C  /* Deck sirfck */
      SUBROUTINE SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRECT,
     &                  WRK,LWORK)
C
C     March 1997 - tsaue Screening
C     921201-hjaaj; partly based on GETFCK
C     written by Henrik Koch and Trygve Helgaker 22-November-1991.
C     DFT modifications T. Helgaker
C
C PURPOSE : Driver routine for the calculation of the two-electron
C           part of the fock matrices from AO integrals.
C           We assume the densities and fock matrices are full
C           squares and without symmtry reduction.
C           The subroutine offers the possibility of calculating
C           the fock matrices directly or by reading the integrals
C           from disk.
C
C output:
C     FMAT   :  NDMAT Fock matrices with 2-el. cont. added
C  input:
C     FMAT   :  NDMAT Fock matrices
C     DMAT   :  NDMAT density matrices
C     NDMAT  :  number of density matrices and Fock matrices
C     ISYMDM :  point group symmetry of DMAT and FMAT
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
C     DIRECT = .TRUE.   Direct calculation using hermit.
C              .FALSE.  Read integrals from AOTWOINT.
#include <implicit.h>
      LOGICAL   DIRECT
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,NDMAT)
      DIMENSION ISYMDM(NDMAT), IFCTYP(NDMAT), MBAS(8)
      DIMENSION WRK(*)
C
#include <dummy.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D4 = 4.0D0, DMP5 = -0.5D0)
C
C Used from common blocks
C   INFORB: NSYM, NBAS(:), N2BASX
C   INFTAP: LUINTA
C   INFPRI: IPRFCK
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <dftcom.h>
C
      LOGICAL SFIRST
      SAVE SFIRST
      DATA SFIRST /.TRUE./
C
C
#include <memint.h>
      CALL QENTER('SIRFCK')
      IF (DIRECT) THEN
         IF (SFIRST) THEN
            IPRHER = 2*IPRFCK
            CALL SETHER(IPRHER,.FALSE.,WRK,LWORK)
C           CALL SETHER(IPRHER,NEWGEO,WORK,LWORK)
            SFIRST = .FALSE.
         END IF
         DO I = 1,NDMAT
            ISYMDM(I) = ISYMDM(I) - 1
         END DO
         IF (NSYM.GT.1 .AND. DFTADD) THEN
C           HERFCK destroys DMAT when called with symmetry
            KDMAT = 1
            KLAST = KDMAT + N2BASX
            LWRK  = LWORK - KLAST + 1
            IF (KLAST.GT.LWORK) CALL STOPIT('SIRFCK','DMAT',KLAST,LWORK)
            CALL DCOPY(N2BASX,DMAT,1,WRK(KDMAT),1)
            CALL HERFCK(FMAT,WRK(KDMAT),NDMAT,ISYMDM,IFCTYP,IPRFCK,
     &                  WRK(KLAST),LWRK)
         ELSE
            CALL HERFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,IPRFCK,WRK,LWORK)
         END IF
         DO I = 1,NDMAT
            ISYMDM(I) = ISYMDM(I) + 1
         END DO
      ELSE
         CALL SIRFCK2(LUINTA,'AOTWOINT',
     &                FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,WRK,LWORK)
C
C        Add DFT contribution
C
         IF (DFTADD) THEN
            IF (NDMAT.GT.1) CALL QUIT('DFTADD and NDMAT.gt.1')
            CALL DFTEXC(DMAT,1,FMAT,1,.FALSE.,.TRUE.,.FALSE.,.FALSE.,
     &           .FALSE.,.FALSE.,DUMMY,.FALSE.,WRK,LWORK,IPRFCK)
         END IF
C
C
      END IF
C
#if defined (VAR_OLDCODE)
--- keep print option here ??? 921201-hjaaj
C
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
      IF (IPRFCK.GT.2) THEN
         CALL HEADER('Density and Fock matrices in SIRFCK',-1)
         DO 300 I = 1, NDMAT
            ISTR = N2BASX*(I - 1) + 1
            WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
            CALL OUTPUT(DMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
#endif
C
      CALL QEXIT('SIRFCK')
      RETURN
C
      END
C  /* Deck sirfck2 */
      SUBROUTINE SIRFCK2(LU2INT,FN2INT,
     &                   FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,WRK,LWORK)
C
C     Extracted non-Direct code from SIRFCK into SIRFCK2,
C     with file unit and name in parameter list.
C     This means that this can also be used e.g. for
C     the short range Hartree terms in the MC-DFT model.
C                \jkp+hjaaj
C
C output:
C     FMAT   :  NDMAT Fock matrices with 2-el. cont. added
C  input:
C     FMAT   :  NDMAT Fock matrices
C     DMAT   :  NDMAT density matrices
C     NDMAT  :  number of density matrices and Fock matrices
C     ISYMDM :  symmetry of DMAT and FMAT
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
#include <implicit.h>
      CHARACTER*(*) FN2INT
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,NDMAT)
      DIMENSION ISYMDM(NDMAT), IFCTYP(NDMAT), MBAS(8)
      DIMENSION WRK(LWORK)
C
#include <iratdef.h>
#include <dummy.h>
      PARAMETER (D1 = 1.0D0, D4 = 4.0D0, DMP5 = -0.5D0)
C
C Used from common blocks
C   INFORB: NSYM, NBAS(:), N2BASX
C   INFPRI: IPRFCK
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
#include <dftcom.h>
C
C
#include <memint.h>
      CALL QENTER('SIRFCK2')
C
C     Transfer constant factors to DMAT
C     IFCTYP 11,12,13,22; symmetrization or antisymmetrization:
C          factor 2 because factor 0.5 in FCKEN1 and FCKEN2 in order
C          to preserve previous content of FMAT
C     IFCTYP 11,13: factor  2 (see note in FCKD11/FCKD13)
C     IFCTYP 12,22: factor  0.5 from exchange integral
C     IFCTYP 02   : factor -0.5 from exchange integral
      DO I = 1,NDMAT
         IF (IFCTYP(I) .EQ. 11 .OR. IFCTYP(I) .EQ. 13) THEN
            CALL DSCAL(N2BASX,D4,DMAT(1,I),1)
         ELSE IF (IFCTYP(I) .EQ. 2) THEN
            CALL DSCAL(N2BASX,DMP5,DMAT(1,I),1)
         END IF
      END DO
C
C
      IF (LU2INT .GT. 0) CALL GPCLOSE(LU2INT,'KEEP')
      CALL GPOPEN(LU2INT,FN2INT,'OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
C
      CALL REWSPL(LU2INT)
      CALL MOLLAB('BASINFO ',LU2INT,LUPRI)
      READ (LU2INT) MSYM, MBAS, LBUF, NIBUF, NBITS
      IF (MSYM .NE. NSYM) THEN
         CALL QUIT('SIRFCK2: NSYM from AOxxxINT .ne. NSYM in common')
      END IF
#if !defined (SYS_CRAY) && !defined (SYS_T3D) && !defined (SYS_T90)
      IF (NBITS .NE. NIBUF*8) THEN
         CALL QUIT('SIRFCK2: NBITS from AOxxxINT .ne. NIBUF*8')
      END IF
#endif
      DO I = 1,MSYM
         IF (NBAS(I) .NE. MBAS(I)) GO TO 110
      END DO
      GO TO 120
  110    WRITE (LUPRI,'(/A/3A,8I5)')
     &      ' SIRFCK2 FATAL ERROR:',
     &      ' NBAS from ',FN2INT,' :',(MBAS(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)')
     &      ' NBAS from /INFORB/ :',(NBAS(I),I=1,NSYM)
         CALL QUIT('SIRFCK2: NBAS(:) from AOxxxINT not consistent')
  120 CONTINUE
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('INTE',KINT ,LENINT,WRK,KFREE,LFREE)
      KIINT = KINT + LBUF
      CALL MEMGET('INTE',KIIN4,4*LBUF,WRK,KFREE,LFREE)
C
      CALL MOLLAB('BASTWOEL',LU2INT,LUPRI)
      JBUF = 608
C
    1 CONTINUE
C
#if defined (VAR_SPLITFILES)         
         CALL READSI(LU2INT,LENINT,WRK(KINT),JBUF)
#else
         CALL READI(LU2INT,LENINT,WRK(KINT))
#endif
         CALL AOLAB4(WRK(KIINT),LBUF,NIBUF,NBITS,WRK(KIIN4),NINT)
         IF (NINT .EQ. 0)  GOTO 1
         IF (NINT .EQ. -1) GOTO 2
C
         DO 200 I = 1,NDMAT
            IF (IFCTYP(I) .EQ. 13) THEN
C           ... e.g. singlet, symmetric density matrix
              CALL FCKD13(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 22
     &          .OR. IFCTYP(I) .EQ. 23      
     &          .OR. IFCTYP(I) .EQ. 12) THEN
C           ... Symmetric or antisymmetric, Only Exchange
C               (type 23: Coulomb term is zero because antisymmetric)
C           ... e.g. singlet, antisymmetric density matrix OR
C           ... or   triplet,     symmetric density matrix OR
C           ... or   triplet, antisymmetric density matrix
              CALL FCKD12(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 11) THEN
C           ... Symmetric, Only Coulomb
              CALL FCKD11(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 3) THEN
C           ... Non-symmetric, Coulomb+Exchange
C           ... e.g. singlet, non-symmetric (general) density matrix
              CALL FCKD03(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 2) THEN
C           ... Non-symmetric, Only Exchange
C           ... e.g. triplet, non-symmetric (general) density matrix
              CALL FCKD02(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 1) THEN
C           ... Non-symmetric, Only Coulomb
              CALL FCKD01(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 0 .OR. IFCTYP(I) .EQ. 21) THEN
C              00:  nothing
C              21:  antisymmetric, only Coulomb gives a zero matrix
               GOTO 200
            ELSE
              GO TO 8000
            END IF
  200    CONTINUE
C
      GOTO 1
C
    2 CONTINUE
      CALL GPCLOSE(LU2INT,'KEEP')
C
C     Finish symmetric and antisymmetric Fock matrices
C     Restore original density matrices
C
      DO 300 I = 1,NDMAT
C
         IF (IFCTYP(1) .EQ. 11 .OR. IFCTYP(I) .EQ. 12 .OR.
     &       IFCTYP(I) .EQ. 13) THEN
            CALL FCKEN1(FMAT(1,I),NBAST)
         ELSE IF (IFCTYP(I) .EQ. 22 .OR. IFCTYP(I) .EQ. 23) THEN
            CALL FCKEN2(FMAT(1,I),NBAST)
         END IF
C
         IF (IFCTYP(I) .EQ. 11 .OR. IFCTYP(I) .EQ. 13) THEN
            DFAC = D1 / D4
            CALL DSCAL(N2BASX,DFAC,DMAT(1,I),1)
         ELSE IF (IFCTYP(I) .EQ. 2) THEN
            DFAC = D1 / DMP5
            CALL DSCAL(N2BASX,DFAC,DMAT(1,I),1)
         END IF
C
  300 CONTINUE
C
      CALL QEXIT('SIRFCK2')
      RETURN
C
 8000 CONTINUE
      WRITE (LUPRI,'(/A/,(8I10))')
     &   ' FATAL ERROR SIRFCK2: Unknown IFCTYP; dump of IFCTYP:',
     &   (IFCTYP(I),I=1,NDMAT)
      CALL QTRACE(LUPRI)
      CALL QUIT('SIRFCK2 ERROR: Unknown IFCTYP')
      END
C  /* Deck fckd13 */
      SUBROUTINE FCKD13(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C Hans Joergen Aa. Jensen October 1994: general routine for
C symmetric singlet Fock matrix from symmetric singlet density matrix
C (based on FCKD03)
C
C This subroutine adds (derivative) two-electron integrals to
C Fock matrices. The Fock matrices are assumed
C to be square matrices in full dimension without symmetry reduction
C in size. Remember to zero out the Fock matrices before starting
C to accumulate.
C
C tsaue Feb 1997: indices permuted for exchange
C HJAaJ October 1994: NOTE:
C     Symmetry is implicitly taken into account by the fact that
C     only non-zero integrals are considered.  The symmetry of DMAT
C     is not taken into account.
C DFT modifications T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP25 = 0.25D0)
      INTEGER P, Q, R, S
#include <inforb.h>
#include <dftcom.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      EFAC = DP25
      IF (DFTRUN) EFAC = DP25*HFXFAC
C     ... because DMAT is multiplied by extra 2 in SIRFCK
C         in order to get full Coulomb contribution:
C         DINT*(DMAT(R,S) + DMAT(S,R)) = 2*DINT*DMAT(R,S)
C         Compare general code in FCKD03.
      DO INT = 1, LENGTH
C        factor "4" is multiplied on DMAT in SIRFCK
C        DINT = 4*BUF(INT)
         DINT =   BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DINT / 2
         IF (R.EQ.S)              DINT = DINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DINT / 2
         EINT = EFAC*DINT
         FMAT(P,Q) = FMAT(P,Q) + DINT*DMAT(R,S)
         FMAT(P,R) = FMAT(P,R) - EINT*DMAT(Q,S)
         FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(P,S)
         FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(P,R)
         FMAT(P,S) = FMAT(P,S) - EINT*DMAT(Q,R)
         FMAT(R,S) = FMAT(R,S) + DINT*DMAT(P,Q)
C        ... DMAT(i,j) = DMAT(j,i) has been used above
      END DO
      RETURN
      END
C  /* Deck fckd12 */
      SUBROUTINE FCKD12(FMAT,XMAT,BUF,IBUF,LENGTH)
C
C     Hans Joergen Aa. Jensen and Rika Kobayashi 24-Apr-1992
C     (based on FCKDIR by Henrik Koch and Trygve Helgaker)
C     970304-tsaue: index permuted
C     941011-hjaaj: renamed from LFKDIR to FCKD12
C                   more than halved the work by using
C                   the antisymmetry also in the exchange part
C     DFT modifications T. Helgaker
C
C     This subroutine adds two-electron integrals multiplied
C     by antisymmetric effective density matrices to Fock matrices
C     (e.g. from London reorthonormalization one-index transformations).
C     The Fock matrices are assumed to be square matrices
C     in full dimension without symmetry reduction in size.
C     Remember to zero out the fock matrices before starting
C     to accumulate.
C
C     The matrix is finished by antisymmetrization in FCKES2.
C
C     FMAT(P,R)  = - 0.5 SUM(q,s) (p q | r s) X(s,q)
C                = + 0.5 SUM(q,s) (p q | r s) X(q,s)
C     Note that X(q,s) = -X(s,q)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D00, DP5 = 0.50D00)
      INTEGER P, Q, R, S
#include <inforb.h>
#include <dftcom.h>
      DIMENSION FMAT(NBAST,NBAST), XMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      EFAC = D1
      IF (DFTRUN) EFAC = HFXFAC
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = EFAC*BUF(INT)
         IF (P.EQ.Q)              GINT = DP5*GINT
         IF (R.EQ.S)              GINT = DP5*GINT
         IF (P.EQ.R .AND. S.EQ.Q) GINT = DP5*GINT
         FMAT(P,R) = FMAT(P,R) - GINT*XMAT(Q,S)
         FMAT(P,S) = FMAT(P,S) - GINT*XMAT(Q,R)
         FMAT(Q,R) = FMAT(Q,R) - GINT*XMAT(P,S)
         FMAT(Q,S) = FMAT(Q,S) - GINT*XMAT(P,R)
CFCKD03 MAERK 941013: bemaerk 0.5 faktor mangler sammenlignet med FCKD03
CFCKD03               and XMAT transposed using XMAT(i,j) = -XMAT(j,i)
CFCKD03  FMAT(P,R) = FMAT(P,R) - EINT*DMAT(S,Q)
CFCKD03  FMAT(P,S) = FMAT(P,S) - EINT*DMAT(R,Q)
CFCKD03  FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(S,P)
CFCKD03  FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckd11 */
      SUBROUTINE FCKD11(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C Hans Joergen Aa. Jensen March 2003: add to FMAT for
C IX = 1 (symmetric DMAT), IY = 1 (only Coulomb)
C (based on FCKD13)
C
C The Fock matrices are assumed to be square matrices
C in full dimension without symmetry reduction in size.
C Remember to zero out the fock matrices before starting
C to accumulate.
C
C HJAaJ October 1994: NOTE:
C     Symmetry is implicitly taken into account by the fact that
C     only non-zero integrals are considered.  The symmetry of DMAT
C     is not taken into account.
C
#include <implicit.h>
#include <priunit.h>
      INTEGER P, Q, R, S
#include <inforb.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      DO INT = 1, LENGTH
C        factor "4" is multiplied on DMAT in SIRFCK
C        DINT = 4*BUF(INT)
         DINT =   BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DINT / 2
         IF (R.EQ.S)              DINT = DINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DINT / 2
         FMAT(P,Q) = FMAT(P,Q) + DINT*DMAT(R,S)
         FMAT(R,S) = FMAT(R,S) + DINT*DMAT(P,Q)
      END DO
      RETURN
      END
C  /* Deck fckd03 */
      SUBROUTINE FCKD03(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     Henrik Koch and Trygve Helgaker 18-NOV-1991.
C     970303-tsaue : index permutation
C     941011-hjaaj: renamed from FCKDI1 to FCKD03
C     DFT modifications T. Helgaker
C
C     This subroutine adds derivative two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.50D00)
      INTEGER P, Q, R, S
#include <inforb.h>
#include <dftcom.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
#include <ibtfun.h>
C
      EFAC = DP5
      IF (DFTRUN) EFAC = DP5*HFXFAC
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = BUF(INT)
         IF (P.EQ.Q)              GINT = DP5*GINT
         IF (R.EQ.S)              GINT = DP5*GINT
         IF (P.EQ.R .AND. S.EQ.Q) GINT = DP5*GINT
C           coulomb:
         FADD = GINT*(DMAT(R,S)+DMAT(S,R))
         FMAT(P,Q) = FMAT(P,Q) + FADD
         FMAT(Q,P) = FMAT(Q,P) + FADD
         FADD = GINT*(DMAT(P,Q)+DMAT(Q,P))
         FMAT(R,S) = FMAT(R,S) + FADD
         FMAT(S,R) = FMAT(S,R) + FADD
C           exchange:
         GINT = EFAC*GINT
         FMAT(P,R) = FMAT(P,R) - GINT*DMAT(Q,S)
         FMAT(R,P) = FMAT(R,P) - GINT*DMAT(S,Q)
         FMAT(P,S) = FMAT(P,S) - GINT*DMAT(Q,R)
         FMAT(S,P) = FMAT(S,P) - GINT*DMAT(R,Q)
         FMAT(Q,R) = FMAT(Q,R) - GINT*DMAT(P,S)
         FMAT(R,Q) = FMAT(R,Q) - GINT*DMAT(S,P)
         FMAT(Q,S) = FMAT(Q,S) - GINT*DMAT(P,R)
         FMAT(S,Q) = FMAT(S,Q) - GINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckd02 */
      SUBROUTINE FCKD02(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     Hans Joergen Aa. Jensen 11-Oct-1994
C     General routine for triplet Fock matrices, based on FCKD03.
C     970304 - tsaue: Index permutation
C
C     This subroutine adds two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP5 = 0.50D00)
      INTEGER P, Q, R, S
#include <inforb.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
#include <ibtfun.h>
C
      DO 100 INT = 1, LENGTH
C        factor "-0.5" is multiplied on DMAT in SIRFCK
C        GINT = -DP5*BUF(INT)
         GINT = BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              GINT = DP5*GINT
         IF (R.EQ.S)              GINT = DP5*GINT
         IF (P.EQ.R .AND. S.EQ.Q) GINT = DP5*GINT
         FMAT(P,R) = FMAT(P,R) + GINT*DMAT(Q,S)
         FMAT(Q,R) = FMAT(Q,R) + GINT*DMAT(P,S)
         FMAT(R,P) = FMAT(R,P) + GINT*DMAT(S,Q)
         FMAT(S,P) = FMAT(S,P) + GINT*DMAT(R,Q)
         FMAT(P,S) = FMAT(P,S) + GINT*DMAT(Q,R)
         FMAT(Q,S) = FMAT(Q,S) + GINT*DMAT(P,R)
         FMAT(R,Q) = FMAT(R,Q) + GINT*DMAT(S,P)
         FMAT(S,Q) = FMAT(S,Q) + GINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckd01 */
      SUBROUTINE FCKD01(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     HJAaJ Mar 2003: based on FCKD03
C
C     The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include <implicit.h>
#include <priunit.h>
      INTEGER P, Q, R, S
      PARAMETER ( DP5 = 0.5D0 )
#include <inforb.h>
#include <dftcom.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = BUF(INT)
         IF (P.EQ.Q)              GINT = DP5*GINT
         IF (R.EQ.S)              GINT = DP5*GINT
         IF (P.EQ.R .AND. S.EQ.Q) GINT = DP5*GINT
C           coulomb:
         FADD = GINT*(DMAT(R,S)+DMAT(S,R))
         FMAT(P,Q) = FMAT(P,Q) + FADD
         FMAT(Q,P) = FMAT(Q,P) + FADD
         FADD = GINT*(DMAT(P,Q)+DMAT(Q,P))
         FMAT(R,S) = FMAT(R,S) + FADD
         FMAT(S,R) = FMAT(S,R) + FADD
  100 CONTINUE
      RETURN
      END
C  /* Deck fcken1 */
      SUBROUTINE FCKEN1(FMAT,NBAST)
C
C     Hans Joergen Aa. Jensen October 1994.
C
C     This subroutine completes symmetric Fock matrices
C     accumulated by FCKD13 or FCKD02 by symmetrization.
C
C     Note: 0.5 factor is required in order to preserve
C     symmetric content of FMAT from before SIRFCK call.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER ( DP5 = 0.5D0 )
      DIMENSION FMAT(NBAST,NBAST)
      INTEGER   P, Q
C
      DO P = 1,NBAST
         DO Q = 1,P
            FPQ = DP5*( FMAT(P,Q) + FMAT(Q,P) )
            FMAT(P,Q) = FPQ
            FMAT(Q,P) = FPQ
         END DO
      END DO
      RETURN
      END
C  /* Deck fcken2 */
      SUBROUTINE FCKEN2(FMAT,NBAST)
C
C     Hans Joergen Aa. Jensen October 1994.
C
C     This subroutine completes antisymmetric Fock matrices
C     accumulated by FCKD12 by antisymmetrization.
C
C     Note: 0.5 factor is required in order to preserve
C     symmetric content of FMAT from before SIRFCK call.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION FMAT(NBAST,NBAST)
      PARAMETER (DP5 = 0.5D0)
      INTEGER   P, Q
C
      DO P = 1,NBAST
         DO Q = 1,P
            FPQ = DP5 * ( FMAT(P,Q) - FMAT(Q,P) )
            FMAT(P,Q) =  FPQ
            FMAT(Q,P) = -FPQ
         END DO
      END DO
      RETURN
      END
