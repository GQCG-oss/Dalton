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
C
C$Id: sirfck.F,v 1.5 2001-10-01 13:22:20 vebjornb Exp $
C
#include <single.h>
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
C     DIRECT = .TRUE.   Direct calculation using hermit.
C              .FALSE.  Read integrals from AOTWOINT.
#include <implicit.h>
      LOGICAL   DIRECT
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,NDMAT)
      DIMENSION ISYMDM(NDMAT), IFCTYP(NDMAT), MBAS(8)
      DIMENSION WRK(*)
C
#include <iratdef.h>
#include <dummy.h>
      PARAMETER (D1 = 1.0D0, D4 = 4.0D0, DMP5 = -0.5D0)
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
C
C        Transfer constant factors to DMAT
C        IFCTYP 12,13,22; symmetrization or antisymmetrization:
C          factor 2 because factor 0.5 in FCKEN1 and FCKEN2 in order
C          to preserve previous content of FMAT
C        IFCTYP 13: factor 2 (see note in FCKDS1)
C        IFCTYP 12,22: factor 0.5 from exchange integral
C        IFCTYP 02: factor -0.5 from exchange integral
         DO 100 I = 1,NDMAT
            IF (IFCTYP(I) .EQ. 13) THEN
               CALL DSCAL(N2BASX,D4,DMAT(1,I),1)
            ELSE IF (IFCTYP(I) .EQ. 2) THEN
               CALL DSCAL(N2BASX,DMP5,DMAT(1,I),1)
            END IF
  100    CONTINUE
C
C
         IF (LUINTA .GT. 0) CALL GPCLOSE(LUINTA,'KEEP')
         CALL GPOPEN(LUINTA,'AOTWOINT','OLD',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
C
         CALL REWSPL(LUINTA)
         CALL MOLLAB('BASINFO ',LUINTA,LUPRI)
         READ (LUINTA) MSYM, MBAS, LBUF, NIBUF, NBITS
         IF (MSYM .NE. NSYM) THEN
            CALL QUIT('SIRFCK: NSYM from AOTWOEL .ne. NSYM in common')
         END IF
#if !defined (SYS_CRAY) && !defined (SYS_T3D) && !defined (SYS_T90)
         IF (NBITS .NE. NIBUF*8) THEN
            CALL QUIT('SIRFCK: NBITS from AOTWOEL .ne. NIBUF*8')
         END IF
#endif
         DO I = 1,MSYM
            IF (NBAS(I) .NE. MBAS(I)) GO TO 110
         END DO
         GO TO 120
  110    WRITE (LUPRI,'(/A/A,8I5)')
     &      ' SIRFCK FATAL ERROR:',
     &      ' NBAS from AOTWOEL :',(MBAS(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)')
     &      ' NBAS from /INFORB/:',(NBAS(I),I=1,NSYM)
         CALL QUIT('SIRFCK: NBAS(:) from AOTWOEL not consistent')
  120    CONTINUE
C
         LENINT = IRAT*LBUF + NIBUF*LBUF + 1
         KFREE = 1
         LFREE = LWORK
         CALL MEMGET('INTE',KINT ,LENINT,WRK,KFREE,LFREE)
	 KIINT = KINT + LBUF
         CALL MEMGET('INTE',KIIN4,4*LBUF,WRK,KFREE,LFREE)
C
         CALL MOLLAB('BASTWOEL',LUINTA,LUPRI)
         JBUF = 608
C
    1    CONTINUE
C
#if defined (VAR_SPLITFILES)         
         CALL READSI(LUINTA,LENINT,WRK(KINT),JBUF)
#else
         CALL READI(LUINTA,LENINT,WRK(KINT))
#endif
         CALL AOLAB4(WRK(KIINT),LBUF,NIBUF,NBITS,WRK(KIIN4),NINT)
         IF (NINT .EQ. 0)  GOTO 1
         IF (NINT .EQ. -1) GOTO 2
C
         DO 200 I = 1,NDMAT
            IF (IFCTYP(I) .EQ. 13) THEN
C           ... singlet, symmetric density matrix
              CALL FCKDS1(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 22
     &          .OR. IFCTYP(I) .EQ. 12) THEN
C           ... singlet, antisymmetric density matrix OR
C           ... triplet,     symmetric density matrix OR
C           ... triplet, antisymmetric density matrix
              CALL FCKDS2(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 3) THEN
C           ... singlet, non-symmetric (general) density matrix
              CALL FCKDS3(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 2) THEN
C           ... triplet, non-symmetric (general) density matrix
              CALL FCKDT3(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 0) THEN
               GOTO 200
            ELSE
              GO TO 8000
            END IF
  200    CONTINUE
C
         GOTO 1
C
    2    CONTINUE
         CALL GPCLOSE(LUINTA,'KEEP')
C
C        Finish symmetric and antisymmetric Fock matrices
C        Restore original density matrices
C
         DO 300 I = 1,NDMAT
            IF (IFCTYP(I) .EQ. 13) THEN
               DFAC = D1/D4
               CALL DSCAL(N2BASX,DFAC,DMAT(1,I),1)
               CALL FCKEN1(FMAT(1,I),NBAST)
            ELSE IF (IFCTYP(I) .EQ. 12) THEN
               CALL FCKEN1(FMAT(1,I),NBAST)
            ELSE IF (IFCTYP(I) .EQ. 22) THEN
               CALL FCKEN2(FMAT(1,I),NBAST)
            ELSE IF (IFCTYP(I) .EQ. 2) THEN
               DFAC = D1/DMP5
               CALL DSCAL(N2BASX,DFAC,DMAT(1,I),1)
            END IF
  300    CONTINUE
C
      END IF
C
C     Add DFT contribution
C
      IF (DFTADD) THEN
         CALL DFTEXC(DMAT,1,FMAT,1,.FALSE.,.TRUE.,.FALSE.,.FALSE.,
     &        .FALSE.,.FALSE.,DUMMY,.FALSE.,WRK,LWORK,IPRFCK)
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
 8000 CONTINUE
      WRITE (LUPRI,'(/A/,(8I10))')
     &   ' FATAL ERROR SIRFCK: Unknown IFCTYP; dump of IFCTYP:',
     &   (IFCTYP(I),I=1,NDMAT)
      CALL QTRACE(LUPRI)
      CALL QUIT('SIRFCK ERROR: Unknown IFCTYP')
      END
C  /* Deck fckds1 */
      SUBROUTINE FCKDS1(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C Hans Joergen Aa. Jensen October 1994: general routine for
C symmetric singlet Fock matrix from symmetric singlet density matrix
C (based on FCKDS3)
C
C This subroutine adds derivative two-electron integrals to
C Fock matrices. The Fock matrices are assumed
C to be square matrices in full dimension without symmetry reduction
C in size. Remember to zero out the fock matrices before starting
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
      PARAMETER (DP5 = 0.50D00, DP25 = 0.25D0)
      INTEGER P, Q, R, S
#include <inforb.h>
#include <dftcom.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
#include <ibtfun.h>
C
      EFAC = DP25
      IF (DFTRUN) EFAC = DP25*HFXFAC
      DO 100 INT = 1, LENGTH
         DINT = BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DP5*DINT
         IF (R.EQ.S)              DINT = DP5*DINT
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DP5*DINT
         EINT = EFAC*DINT
C        ... because DMAT is multiplied by extra 2 in SIRFCK
C            in order to get full Coulomb contribution:
C            DINT*(DMAT(R,S) + DMAT(S,R)) = 2*DINT*DMAT(R,S)
C            Compare general code in FCKDS3.
         FMAT(P,Q) = FMAT(P,Q) + DINT*DMAT(S,R)
         FMAT(P,R) = FMAT(P,R) - EINT*DMAT(Q,S)
         FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(P,S)
         FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(P,R)
         FMAT(P,S) = FMAT(P,S) - EINT*DMAT(Q,R)
         FMAT(R,S) = FMAT(R,S) + DINT*DMAT(Q,P)
C        ... DMAT(i,j) = DMAT(j,i) has been used above
  100 CONTINUE
      RETURN
      END
C  /* Deck fckds2 */
      SUBROUTINE FCKDS2(FMAT,XMAT,BUF,IBUF,LENGTH)
C
C     Hans Joergen Aa. Jensen and Rika Kobayashi 24-Apr-1992
C     (based on FCKDIR by Henrik Koch and Trygve Helgaker)
C     970304-tsaue: index permuted
C     941011-hjaaj: renamed from LFKDIR to FCKDS2
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
#include <ibtfun.h>
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
CFCKDS3 MAERK 941013: bemaerk 0.5 faktor mangler sammenlignet med FCKDS3
CFCKDS3               and XMAT transposed using XMAT(i,j) = -XMAT(j,i)
CFCKDS3  FMAT(P,R) = FMAT(P,R) - EINT*DMAT(S,Q)
CFCKDS3  FMAT(P,S) = FMAT(P,S) - EINT*DMAT(R,Q)
CFCKDS3  FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(S,P)
CFCKDS3  FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckds3 */
      SUBROUTINE FCKDS3(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     Henrik Koch and Trygve Helgaker 18-NOV-1991.
C     970303-tsaue : index permutation
C     941011-hjaaj: renamed from FCKDI1 to FCKDS3
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
C  /* Deck fckdt3 */
      SUBROUTINE FCKDT3(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     Hans Joergen Aa. Jensen 11-Oct-1994
C     General routine for triplet Fock matrices, based on FCKDS3.
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
         FMAT(R,P) = FMAT(R,P) + GINT*DMAT(S,Q)
         FMAT(P,S) = FMAT(P,S) + GINT*DMAT(Q,R)
         FMAT(S,P) = FMAT(S,P) + GINT*DMAT(R,Q)
         FMAT(Q,R) = FMAT(Q,R) + GINT*DMAT(P,S)
         FMAT(R,Q) = FMAT(R,Q) + GINT*DMAT(S,P)
         FMAT(Q,S) = FMAT(Q,S) + GINT*DMAT(P,R)
         FMAT(S,Q) = FMAT(S,Q) + GINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fcken1 */
      SUBROUTINE FCKEN1(FMAT,NBAST)
C
C     Hans Joergen Aa. Jensen October 1994.
C
C     This subroutine completes symmetric Fock matrices
C     accumulated by FCKDS1 or FCKDT1 by symmetrization.
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
C     accumulated by FCKDS2 by antisymmetrization.
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
