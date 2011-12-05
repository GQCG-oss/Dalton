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
C File : herdkh.F
C
C Subroutines related to implementation of DKH2 - 2. order
C Douglas-Kroll-Hess - in Dalton.
C
C Original implementation by Kenneth Ruud, September 2000.
C
C  /* Deck dkintF */
      SUBROUTINE DKINTF(WORK,LWORK)
C
C     Driver for calculating Douglas-Kroll-Hess integrals in the primitive
C     basis and transforming them to the appropriate form before 
C     writing them to file
C
C     K.Ruud, September-00
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
C
      DIMENSION WORK(LWORK)
      LOGICAL HUCKELBKP
C
#include <inftap.h>
#include <cbiher.h>
#include <cbihr1.h>
#include <cbihr2.h>
#include <cbieri.h>
#include <cbihrs.h>
#include <gnrinf.h>
#include <huckel.h>
C
      CALL QENTER('DKINTF')
      CALL GETTIM(TIMHER,WALHER)
C
      CALL HERINI
      HAMILT = .FALSE.
      PVPINT = .TRUE.
      POTENE = .TRUE.
      KINENE = .TRUE.
      PROPRI = .FALSE.
      HUCKELBKP = DOHUCKEL
      DOHUCKEL = .FALSE.
      ONEPRP = .TRUE.
      RUNONE = .TRUE.
      RUNSUP = .FALSE.
      RUNTWO = .FALSE.
C
C     Process basis set, although now as a decontracted set without
C     normalization of the primitive orbitals (UNCONT .true.)
C     (READIN will be called again in HERINP and generate the
C     proper contracted basis and the info necessary to
C     transform Douglas-Kroll-Hess integrals to contracted basis)
C      
      CALL DKREIN(WORK,LWORK)
C
      TIMSTR = SECOND()
      TIMINP = SECOND() - TIMSTR
      CALL FLSHFO(LUPRI)
      CALL GPOPEN(LUPROP,'DKINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL HERINT(WORK,LWORK)
      CALL GPOPEN(LUONEL,'AOONEINT',' ',' ',' ',IDUMMY,.FALSE.)
      CALL GPCLOSE(LUONEL,'DELETE')
      IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
C
C     We now calculate all the Douglas-Kroll-Hess quantities in the primitive
C     basis
C
      CALL DKDRIV(WORK,LWORK,IPRDEF)
C
      DOHUCKEL = HUCKELBKP
      RDINPC = .FALSE.
      RUNSUP = SUPMAT
      PVPINT = .FALSE.
      RUNTWO = .NOT.NOTWO
      HAMILT = .TRUE.
      CALL GETTIM(TEND,WEND)
      TIMHER = TEND - TIMHER
      WALHER = WEND - WALHER
      CALL TIMTXT('>>>> Total CPU  time used in DKINT:',TIMHER,LUPRI)
      CALL TIMTXT('>>>> Total wall time used in DKINT:',WALHER,LUPRI)
      CALL QEXIT('DKINTF')
      RETURN
      END
C/* Deck dkrein */
      SUBROUTINE DKREIN(WORK,LWORK)
C
C     Read specification of primitive basis (i.e. decontracted)
C     for calculation of Douglas-Kroll-Hess integrals.
C
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>

#include <priunit.h>
      DIMENSION WORK(LWORK)
      LOGICAL UNCONTBKP
      CHARACTER*7 WORD
#include <ccom.h>
#include <mxsymm.h>
#include <cbirea.h>
C
      CALL QENTER('DKREIN')
C
C     Define cbirea.h for call of READIN below
C
      THRS = 1.0D-15
C
      LUMLCL = -1
      MAXPRI = MAXPRD
      BIGVC  = .FALSE.
      DIRAC  = .FALSE.
      GENCON = .TRUE.
      BASIS  = .FALSE.
      ATOMBA = .FALSE.
      TOLLRN = 0.0D0
      ZCMVAL = 1.0D0
      LCMMAX = -1
C
C     We also need two variables from cbiher.h
C
      CALL READ_CBIHER(IPREAD, INPTST)
C
C     An important input for the READIN procedure is the number of primitives
C     per block in the input file. First check if this has been reset.
C
      MAXPRI = MAXPRD
      LUTMP = LUCMD
      IF (LUCMD .LE. 0) CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ',
     &     'FORMATTED',IDUMMY,.FALSE.)
      REWIND (LUCMD,IOSTAT=IOS)
 5000 READ (LUCMD,'(A7)',END=5010) WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. '.MAXPRI') THEN
         READ (LUCMD,*) MAXPRI
         GOTO 5010
      ELSE
         GOTO 5000
      END IF
 5010 CONTINUE
      IF (LUTMP .LE.0) CALL GPCLOSE(LUCMD,'KEEP')
C
C     Now call READIN to get the specification of the primitive basis.
C     Note that UNCONT = .true., to force decontracted specification
C     for Douglas-Kroll-Hess integrals.
C
      UNCONTBKP = UNCONT
      UNCONT = .TRUE.
      CALL READIN(WORK,LWORK,.TRUE.)
      UNCONT = UNCONTBKP
      CALL SETDCH
      CALL QEXIT('DKREIN')
      RETURN
      END
C/* Deck dkdriv */
      SUBROUTINE DKDRIV(WORK,LWORK,IPRINT)
C
C     Driver routine for calculation of the various Douglas-Kroll-Hess matrices
C     in the primitive basis
C
C     K.Ruud, September 2000
C
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      DIMENSION WORK(LWORK)
#include <shells.h>
#include <symmet.h>
#include <inftap.h>
#include <ibtfun.h>
C
      CALL QENTER('DKDRIV')
C
      CALL GPOPEN(LUPROP,'DKINT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
      NBAST = 0
      DO 400 ISHELL = 1, KMAX
         DO 410 KB = 0,MAXREP
            IF (IBTAND(KB,ISTBAO(ISHELL)) .EQ. 0)
     *            NBAST = NBAST + KHKT(ISHELL)
  410    CONTINUE
  400 CONTINUE
      NNBASX = NBAST*(NBAST + 1)/2
      N2BASX = NBAST*NBAST
C
C     Memory allocation
C
      KOVERL = 1
      KXMAT  = KOVERL + NNBASX
      KEIGVC = KXMAT  + N2BASX
      KSQMT1 = KEIGVC + N2BASX
      KSQMT2 = KSQMT1 + N2BASX
      KTRMT1 = KSQMT2 + N2BASX
      KTRMT2 = KTRMT1 + NNBASX
      KE1    = KTRMT2 + NNBASX
      KA1    = KE1    + NBAST
      KEK    = KA1    + NBAST
      KD1    = KEK    + NBAST
      KWEWMT = KD1    + NBAST
      KVMAT  = KWEWMT + N2BASX
      KPVPMT = KVMAT  + NNBASX
      KLAST  = KPVPMT + NNBASX
      IF (KLAST .GT. LWORK) CALL STOPIT('DKDRIV','DKDRI1',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
C
C     Call the driver routine
C     
      CALL DKDRI1(WORK(KOVERL),WORK(KXMAT),WORK(KEIGVC),WORK(KSQMT1),
     &            WORK(KSQMT2),WORK(KTRMT1),WORK(KTRMT2),WORK(KWEWMT),
     *            WORK(KE1),WORK(KA1),WORK(KEK),WORK(KD1),WORK(KVMAT),
     &            WORK(KPVPMT),WORK(KLAST),LWRK,NBAST,NNBASX,N2BASX,
     &            IPRINT)
C
      CALL QEXIT('DKDRIV')
      RETURN
      END
C/* Deck dkdri1 */
      SUBROUTINE DKDRI1(OVERLP,XMAT,EIGVEC,SQMAT1,SQMAT2,TRMAT1,TRMAT2,
     &                  WEWMAT,E1DIAG,A1DIAG,EKDIAG,D1DIAG,VMAT,PVPMAT,
     &                  WORK,LWORK,NBAST,NNBASX,N2BASX,IPRINT)
C
C     K.Ruud, Sept.-2000
C
#include <implicit.h>
#include <codata.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, CVEL2 = CVEL*CVEL, DM1 = -D1)
      CHARACTER*8 RTNLBL(2)
      LOGICAL FNDLAB
      DIMENSION OVERLP(NNBASX), XMAT(NBAST,NBAST), EIGVEC(NBAST,NBAST),
     &          SQMAT1(NBAST,NBAST), SQMAT2(NBAST,NBAST),
     &          TRMAT1(NNBASX), TRMAT2(NNBASX), E1DIAG(NBAST),
     &          A1DIAG(NBAST), EKDIAG(NBAST), D1DIAG(NBAST),
     &          WEWMAT(N2BASX), VMAT(NNBASX), PVPMAT(NNBASX)
      DIMENSION WORK(LWORK)
#include <dkinf.h>      
#include <priunit.h>
#include <inftap.h>
#include <aosotr.h>
#include <symmet.h>
C
C     Allocate temporary storage
C
      KSCR1 = 1
      KSCR2 = KSCR1 + N2BASX
      KLAST = KSCR2 + N2BASX
      IF (KLAST .GT. LWORK) CALL STOPIT('DKDRI1','START ',KLAST,LWORK)
C
C     Read in overlap matrix
C
      IF (FNDLAB('OVERLAP ',LUPROP)) THEN
         CALL READT(LUPROP,NNBASX,OVERLP)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Overlap matrix read from LUPROP in DKDRI1',-1)
            CALL OUTPAK(OVERLP,NBAST,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKDRI1 error: Overlap matrix'//
     *      ' not found on LUPROP.'
         CALL QUIT('HUCDR1 error: property not found on LUPROP')
      END IF
C
C     Diagonalize the overlap matrix
C
      CALL DCOPY(NNBASX,OVERLP,1,TRMAT1,1)
      CALL DUNIT(EIGVEC,NBAST)
      CALL JACO(TRMAT1,EIGVEC,NBAST,NBAST,NBAST,WORK(KSCR1),WORK(KSCR2))
C
      CALL DZERO(XMAT,N2BASX)
      DO I = 1, NBAST
         II = I*(I + 1)/2
         IF (TRMAT1(II) .LT. 1.0D-6) THEN
            WRITE (LUPRI,*) ' >>> ERROR DKDRIV: Linear dependency in '//
     &      'primitive basis used in Douglas-Kroll-Hess ', I, TRMAT1(II)
            write (lupri,*) 'The corresponding eigenvector:'
            call output(eigvec(1,i),1,1,1,nbast,1,nbast,-1,lupri)
            DO J = 1,NBAST
               II = J*(J + 1)/2
               E1DIAG(J) = TRMAT1(II)
            END DO
            write (lupri,*) 'All eigenvalues of the overlap matrix:'
            call output(E1DIAG,1,1,1,nbast,1,nbast,-1,lupri)
            CALL QUIT('Linear dependency in primitive Douglas-Kroll'//
     &           '-Hess basis')
         END IF
         E1DIAG(I) = D1/SQRT(TRMAT1(II))
      END DO
      DO I = 1, NBAST
         CONST =E1DIAG(I)
         DO J = 1, NBAST
            XMAT(J,I) = EIGVEC(J,I)*CONST
         END DO
      END DO
C
C     Read in the kinetic energy and transform to orthonormal basis
C     Note that we reuse OVERLP to save memory
C     
      REWIND (LUPROP)
      IF (FNDLAB('KINENERG',LUPROP)) THEN
         CALL READT(LUPROP,NNBASX,OVERLP)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Kinetic energy matrix read from LUPROP '//
     &                  'in DKDRI1',-1)
            CALL OUTPAK(OVERLP,NBAST,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKDRI1 error: Kinetic energy'//
     *      ' matrix not found on LUPROP.'
         CALL QUIT('HUCDR1 error: property not found on LUPROP')
      END IF
C     
      CALL UTHU(OVERLP,TRMAT1,XMAT,WORK(KSCR1),NBAST,NBAST)
C
C     Diagonalize the kinetic energy matrix in orthonormal basis
C
      CALL DUNIT(EIGVEC,NBAST)
      CALL JACO(TRMAT1,EIGVEC,NBAST,NBAST,NBAST,WORK(KSCR1),WORK(KSCR2))
C
C     Extract kinetic energy and construct Epi, Api and Dpi
C     
      DO I = 1, NBAST
         II = I*(I + 1)/2
         DIAG = D2*TRMAT1(II)
         EKDIAG(I) = DIAG
         E1DIAG(I) = CVEL*SQRT(DIAG + CVEL2)
         A1DIAG(I) = SQRT((CVEL2 + E1DIAG(I))/(D2*E1DIAG(I)))
         D1DIAG(I) = CVEL*A1DIAG(I)/(CVEL2 + E1DIAG(I))
         E1DIAG(I) = E1DIAG(I) - CVEL2
      END DO
      IF (IPRINT .GT. 3) THEN
         CALL HEADER('Diagonal E1 matrix',-1)
         CALL OUTPUT(E1DIAG,1,NBAST,1,1,NBAST,1,-1,LUPRI)
         CALL HEADER('Diagonal A1 matrix',-1)
         CALL OUTPUT(A1DIAG,1,NBAST,1,1,NBAST,1,-1,LUPRI)
         CALL HEADER('Diagonal Ek matrix',-1)
         CALL OUTPUT(EKDIAG,1,NBAST,1,1,NBAST,1,-1,LUPRI)
         CALL HEADER('Diagonal D1 matrix',-1)
         CALL OUTPUT(D1DIAG,1,NBAST,1,1,NBAST,1,-1,LUPRI)
      END IF
C
C
C     Read in nuclear attraction and pVp integrals from file LUPROP (DKINT)
C     and transform to orthonormal basis. Then we transform them to the 
C     diagonal basis of the kinetic operator.
C
      REWIND (LUPROP)
      IF (FNDLAB('POTENERG',LUPROP)) THEN
         CALL READT(LUPROP,NNBASX,VMAT)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Potential energy matrix read from LUPROP '//
     &                  'in DKDRI1',-1)
            CALL OUTPAK(VMAT,NBAST,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKDRI1 error: potential energy'//
     *      ' matrix not found on LUPROP.'
         CALL QUIT('DKDRI1 error: property not found on LUPROP')
      END IF
      CALL UTHU(VMAT,OVERLP,XMAT,WORK(KSCR1),NBAST,NBAST)
      CALL UTHU(OVERLP,VMAT,EIGVEC,WORK(KSCR1),NBAST,NBAST)
C     
      REWIND (LUPROP)
      IF (FNDLAB('pVpINTEG',LUPROP)) THEN
         CALL READT(LUPROP,NNBASX,PVPMAT)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('pVp integrals read from LUPROP '//
     &                  'in DKDRI1',-1)
            CALL OUTPAK(PVPMAT,NBAST,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKDRI1 error: pVp integral'//
     *      ' matrix not found on LUPROP.'
         CALL QUIT('DKDRI1 error: property not found on LUPROP')
      END IF
C
      CALL UTHU(PVPMAT,WORK(KSCR2),XMAT,WORK(KSCR1),NBAST,NBAST)
      CALL UTHU(WORK(KSCR2),PVPMAT,EIGVEC,WORK(KSCR1),NBAST,NBAST)
C
C     Create second-order Douglas-Kroll-Hess contributions
C
      CALL DZERO(WEWMAT,N2BASX)
      CALL WWEWEW(VMAT,PVPMAT,E1DIAG,A1DIAG,EKDIAG,TRMAT1,TRMAT2,
     &            SQMAT1,SQMAT2,WEWMAT,WORK(KSCR1),NBAST,NNBASX,N2BASX,
     &            IPRINT)
C
C     Backtransform the second-order Douglas-Kroll-Hess matrix to 
C     orthonormal basis and make it square
C
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           EIGVEC,NBAST,
     &           WEWMAT,NBAST,0.D0,
     &           SQMAT1,NBAST)
      CALL DGEMM('N','T',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           EIGVEC,NBAST,0.D0,
     &           WEWMAT,NBAST)
      CALL DSCAL(N2BASX,-1.0D0,WEWMAT,1)
C
C     -1/2 WWE - 1/2 EWW - WEW
C
C
C     We now backtransform the various matrices to primitive basis
C     and construct the full effective Douglas-Kroll-Hess Hamiltonian
C
      CALL DZERO(SQMAT1,N2BASX)
      DO I = 1, NBAST
         SQMAT1(I,I) = E1DIAG(I) - CVEL2
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           EIGVEC,NBAST,
     &           SQMAT1,NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','T',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT2,NBAST,
     &           EIGVEC,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
C     AHA
C
      CALL DZERO(SQMAT1,N2BASX)
      DO I = 1, NBAST
         SQMAT1(I,I) = A1DIAG(I)
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           EIGVEC,NBAST,
     &           SQMAT1,NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','T',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT2,NBAST,
     &           EIGVEC,NBAST,0.D0,
     &           SQMAT1,NBAST)
      CALL DSPTSI(NBAST,OVERLP,WORK(KSCR1))
      CALL DGEMM('T','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           WORK(KSCR1),NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,-1.D0,
     &           SQMAT2,NBAST,
     &           SQMAT1,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
C     DWD
C
      CALL DZERO(SQMAT1,N2BASX)
      DO I = 1, NBAST
         SQMAT1(I,I) = D1DIAG(I)
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           EIGVEC,NBAST,
     &           SQMAT1,NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','T',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT2,NBAST,
     &           EIGVEC,NBAST,0.D0,
     &           SQMAT1,NBAST)
      CALL DSPTSI(NBAST,WORK(KSCR2),WORK(KSCR1))
      CALL DGEMM('T','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           WORK(KSCR1),NBAST,0.D0,
     &           SQMAT2,NBAST)
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,-1.D0,
     &           SQMAT2,NBAST,
     &           SQMAT1,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
C     Pack and backtransform the final effect Hamiltonian to primitive basis
C
      CALL DSITSP(NBAST,WEWMAT,TRMAT1)
      CALL DGEINV(NBAST,XMAT,SQMAT1,WORK(KSCR1),WORK(KSCR2),INFO)
      CALL UTHU(TRMAT1,TRMAT2,SQMAT1,WORK(KSCR1),NBAST,NBAST)
C
C     Write the final effective Hamiltonian to DKINT and keep until we
C     know what kind of contraction coefficients we have
C
      IF (IPRINT .GT. 1) THEN
         CALL HEADER('Final Douglas-Kroll-Hess effective one-electron'//
     &               ' potential (SO basis)',-1)
         CALL OUTPAK(TRMAT2,NBAST,1,LUPRI)
      END IF
C
C     We transform to AO basis in order to more easily contract it
C
      IF (MAXREP .GT. 0) THEN
         IATOM = 0
         CALL TRSOAO(TRMAT2,WORK,LWORK,NBAST,IATOM,IPRINT)
         CALL DCOPY(NNBASX,WORK,1,TRMAT2,1)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Final Douglas-Kroll-Hess effective '//
     &                 'one-electron potential (AO basis)',-1)
            CALL OUTPAK(TRMAT2,NBAST,1,LUPRI)
         END IF
      END IF
      NPRIM  = NBAST
      NNPRMX = NNBASX
      REWIND(LUPROP)
      CALL WRTPRO(TRMAT2,NNBASX,'DKINTPRI',RTNLBL)
      CALL GPCLOSE(LUPROP,'KEEP')
C
C     End of DKDRIV
C     
      RETURN
      END
C/* Deck wwewew */
      SUBROUTINE WWEWEW(VMAT,PVPMAT,E1DIAG,A1DIAG,EKDIAG,VTMP,GTMP,
     &                  SQMAT1,SQMAT2,WEWMAT,RDIAG,NBAST,NNBASX,N2BASX,
     &                  IPRINT)
C
C     Subroutine for calculating second-order Douglas-Kroll-Hess contributions
C
#include <implicit.h>
#include <codata.h>
      PARAMETER (CVEL2 = CVEL*CVEL, DP5 = 0.50D0)
#include <priunit.h>
      DIMENSION VMAT(NNBASX), PVPMAT(NNBASX), E1DIAG(NBAST),
     &          A1DIAG(NBAST), EKDIAG(NBAST), WEWMAT(NBAST,NBAST),
     &          SQMAT1(NBAST,NBAST), SQMAT2(NBAST,NBAST),
     &          VTMP(NNBASX), GTMP(NNBASX), RDIAG(NBAST)
C
      DO I = 1, NBAST
         E1DIAG(I) = E1DIAG(I) + CVEL2
         RDIAG(I) = CVEL/(CVEL2 + E1DIAG(I))
      END DO
C
      IJ = 1
      DO I = 1, NBAST
         DO J = 1, I
            VTMP(IJ) = VMAT(IJ)/(E1DIAG(I) + E1DIAG(J))
            GTMP(IJ) = PVPMAT(IJ)/(E1DIAG(I) + E1DIAG(J))
            SQMAT1(I,J) = A1DIAG(I)*RDIAG(I)*GTMP(IJ)*A1DIAG(J)
     &                   *A1DIAG(J)
            SQMAT1(J,I) = A1DIAG(J)*RDIAG(J)*GTMP(IJ)*A1DIAG(I)
     &                   *A1DIAG(I)
            SQMAT2(I,J) = RDIAG(I)*VTMP(IJ)*A1DIAG(J)
     &                  - GTMP(IJ)*A1DIAG(J)*RDIAG(J)/EKDIAG(I)
            SQMAT2(J,I) = RDIAG(J)*VTMP(IJ)*A1DIAG(I)
     &                  - GTMP(IJ)*A1DIAG(I)*RDIAG(I)/EKDIAG(J)
            IJ = IJ + 1
         END DO
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           SQMAT2,NBAST,0.D0,
     &           WEWMAT,NBAST)
      IF (IPRINT .GT. 6) THEN
         CALL AROUND('1. transformation of W-matrix')
         CALL OUTPUT(WEWMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
C
      IJ = 1
      DO I = 1, NBAST
         DO J = 1, I
           SQMAT1(I,J) = A1DIAG(I)*VTMP(IJ)*A1DIAG(J)*A1DIAG(J)*RDIAG(J)
           SQMAT1(J,I) = A1DIAG(J)*VTMP(IJ)*A1DIAG(I)*A1DIAG(I)*RDIAG(I)
           SQMAT2(I,J) = -EKDIAG(I)*RDIAG(I)*VTMP(IJ)*A1DIAG(J)
     &                 + GTMP(IJ)*A1DIAG(J)*RDIAG(J)
           SQMAT2(J,I) = -EKDIAG(J)*RDIAG(J)*VTMP(IJ)*A1DIAG(I)
     &                 + GTMP(IJ)*A1DIAG(I)*RDIAG(I)
           IJ = IJ + 1
        END DO
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           SQMAT2,NBAST,1.D0,
     &           WEWMAT,NBAST)
      IF (IPRINT .GT. 6) THEN
         CALL AROUND('2. transformation of W-matrix')
         CALL OUTPUT(WEWMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
C
C     Make 1/2 EW*W + 1/2 W*WE
C
      DO I = 1, NBAST
         DO J = 1, NBAST
            WEWMAT(I,J) = DP5*WEWMAT(I,J)*(E1DIAG(I)+E1DIAG(J))
         END DO
      END DO
      IF (IPRINT .GT. 6) THEN
         CALL AROUND('3. transformation of W-matrix')
         CALL OUTPUT(WEWMAT,1,NBAST,1,NBAST,NBAST,NBAST,-1,LUPRI)
      END IF
C
C     Start constructing WEW
C
      IJ = 1
      DO I = 1, NBAST
         DO J = 1, I
           SQMAT1(I,J) = A1DIAG(I)*RDIAG(I)*GTMP(IJ)*A1DIAG(J)
     &                  *A1DIAG(J)*E1DIAG(J)
           SQMAT1(J,I) = A1DIAG(J)*RDIAG(J)*GTMP(IJ)*A1DIAG(I)
     &                  *A1DIAG(I)*E1DIAG(I)
           SQMAT2(I,J) = RDIAG(I)*VTMP(IJ)*A1DIAG(J)
     &                 - GTMP(IJ)*A1DIAG(J)*RDIAG(J)/EKDIAG(I)
           SQMAT2(J,I) = RDIAG(J)*VTMP(IJ)*A1DIAG(I)
     &                 - GTMP(IJ)*A1DIAG(I)*RDIAG(I)/EKDIAG(J)
           IJ = IJ + 1
        END DO
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           SQMAT2,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
      IJ = 1
      DO I = 1, NBAST
         DO J = 1, I
            SQMAT1(I,J) = A1DIAG(I)*VTMP(IJ)*A1DIAG(J)*A1DIAG(J)
     &                   *RDIAG(J)*E1DIAG(J)
            SQMAT1(J,I) = A1DIAG(J)*VTMP(IJ)*A1DIAG(I)*A1DIAG(I)
     &                   *RDIAG(I)*E1DIAG(I)
            SQMAT2(I,J) = -EKDIAG(I)*RDIAG(I)*VTMP(IJ)*A1DIAG(J)
     &                  + GTMP(IJ)*A1DIAG(J)*RDIAG(J)
            SQMAT2(J,I) = -EKDIAG(J)*RDIAG(J)*VTMP(IJ)*A1DIAG(I)
     &                  + GTMP(IJ)*A1DIAG(I)*RDIAG(I)
            IJ = IJ + 1
         END DO
      END DO
      CALL DGEMM('N','N',NBAST,NBAST,NBAST,1.D0,
     &           SQMAT1,NBAST,
     &           SQMAT2,NBAST,1.D0,
     &           WEWMAT,NBAST)
C
C     The entire matrix is now constructed in the diagonal T basis
C
      RETURN
      END
C/* Deck dkcont */
      SUBROUTINE DKCONT(PROPRI,IPRINT,WORK,LWORK)
C
C     This subroutine reads in the Douglas-Kroll-Hess potential in the primitive
C     basis and contracts it with the contraction coefficients of the 
C     basis set used in the calculation
C
C     K.Ruud, September 2000
C     
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <dkinf.h>      
      LOGICAL FNDLAB, ANTI, PROPRI
      CHARACTER*8 LABINT
      DIMENSION WORK(LWORK)
#include <shells.h>
#include <symmet.h>
C
C     
      NBAST = ISUM(MAXREP+1,NAOS,1)
      NNBASX = NBAST*(NBAST + 1)/2
C
      KWMAT  = 1
      KCONTC = KWMAT + NPRIM*NPRIM 
      KLAST  = KCONTC + NPRIM*NBAST
      KTMP   = KLAST + NPRIM*NPRIM
      LWRK   = LWORK - KLAST + 1
      IF (KTMP .GT. LWORK) CALL STOPIT('DKCONT','PRI2CN',KLAST,LWORK)
      LUDKIN = -9000
      CALL GPOPEN(LUDKIN,'DKINT','OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      IF (FNDLAB('DKINTPRI',LUDKIN)) THEN
         IF (NNPRMX .GT. LWORK) CALL STOPIT('DKCONT','WMAT  ',
     &                                      NNPRMX,LWORK)
         CALL READT(LUDKIN,NNPRMX,WORK(KWMAT))
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Primitive Douglas-Kroll-Hess integrals read '//
     &                  'from LUDKIN in DKCONT',-1)
            CALL OUTPAK(WORK(KWMAT),NPRIM,1,LUPRI)
         END IF
         CALL DCOPY(NNPRMX,WORK(KWMAT),1,WORK(KLAST),1)
         CALL DSPTSI(NPRIM,WORK(KLAST),WORK(KWMAT))
         NOPTYP = 1
         LABINT = 'DKPOTINT'
         ANTI   = .FALSE.
         INTREP = 0
         CALL PRI2CN(WORK(KWMAT),WORK(KCONTC),WORK(KLAST),
     &               LWRK,ANTI,NOPTYP,NBAST,NPRIM,LABINT,
     &               INTREP,PROPRI,IPRINT)
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKCONT error: Primitive Douglas'//
     &      '-Kroll-Hess integrals not found on DKINT file.'
         CALL QUIT('DKCONT error: Integrals not found on DKINT file')
      END IF
      CALL GPCLOSE(LUDKIN,'DELETE')
      RETURN
      END
C/* Deck pri2cn */
      SUBROUTINE PRI2CN(AOINT,CONTC,WORK,LWORK,ANTI,
     &                  NOPTYP,NBAST,NPRIM,LABINT,INTREP,PROPRI,IPRINT)
C
C     Subroutine to transform primitive AO's to contracted SO basis
C
#include <implicit.h>
#include <dummy.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <aovec.h>
#include <pi.h>
      PARAMETER (D2 = 2.0D0, DP5 = 0.50D0, DP75 = 0.75D0, D4 = 4.0D0)
#include <priunit.h>
      LOGICAL ANTI, PROPRI
      CHARACTER*8 RTNLBL(2), LABINT(NOPTYP)
      DIMENSION AOINT(NPRIM,NPRIM), CONTC(NPRIM,NBAST), WORK(LWORK)
      DIMENSION INTREP(NOPTYP)
#include <pgroup.h>
#include <nuclei.h>
#include <primit.h>
#include <shells.h>
#include <symmet.h>
#include <aosotr.h>
#include <ibtfun.h>
C
      ICNT   = 0
      IADR   = 0
      IORBA  = 0
      IF (NOPTYP.GT. 1) CALL QUIT('More than one integral type in'//
     &                            'PRI2CN')
      NELMNT = NBAST*(NBAST + 1)/2
C
      LUCNMT = -9000
      CALL GPOPEN(LUCNMT,'CNTMAT','OLD',' ','FORMATTED',IDUMMY,.FALSE.)
      READ (LUCNMT,*) ((CONTC(I,J),I=1,NPRIM),J=1,NBAST)
      CALL GPCLOSE(LUCNMT,'DELETE')
      IF (IPRINT .GT. 6) THEN
         CALL HEADER('Contraction matrix read from '//
     &               'LUCNMT in PRI2CN',-1)
         CALL OUTPUT(CONTC,1,NPRIM,1,NBAST,NPRIM,NBAST,-1,LUPRI)
      END IF
      CALL DGEMM('T','N',NBAST,NPRIM,NPRIM,1.0D0,CONTC,NPRIM,
     &           AOINT,NPRIM,0.0D0,WORK,NBAST)
      CALL DGEMM('N','N',NBAST,NPRIM,NPRIM,1.0D0,WORK,NBAST,
     &           CONTC,NPRIM,0.0D0,AOINT,NBAST)
      CALL DSITSP(NBAST,AOINT,WORK)
      CALL DCOPY(NELMNT,WORK,1,AOINT,1)
C
C     First we transform to SO basis, then for consistency with WRTUND,
C     we need to symmetry pack the integrals
C
      IF (IPRINT .GE. 6) THEN
         CALL AROUND('Contracted Douglas-Kroll-Hess '//
     &               'integrals in AO basis')
         CALL OUTPAK(AOINT,NBAST,1,LUPRI)
      END IF
      IF (MAXREP .GT. 0) THEN
         IATOM = -1
         CALL TRAOSO(AOINT,WORK,LWORK,NBAST,IATOM,IPRINT)
         CALL DCOPY(NELMNT,WORK,1,AOINT,1)
         CALL PKSYM1(AOINT,AOINT,NAOS,MAXREP+1,2)
      END IF
C
C     Write integrals to file AOPROPER and print them if people are
C     interested
C     
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      IF (ANTI) THEN
         RTNLBL(2)='ANTISYMM'
      ELSE
         RTNLBL(2)='SYMMETRI'
      END IF
C     
      DO 600 I = 1, NOPTYP
         IF (PROPRI .OR. IPRINT .GE. 4) THEN
            CALL AROUND('Integrals of operator: '//LABINT(I))
            WRITE (LUPRI,'(A,2X,A3,A1,I2,A1)')
     &           ' Symmetry of operator:',
     &           REP(INTREP(I)),'(',(INTREP(I) + 1),')'
            CALL OUTPAK(AOINT,NBAST,1,LUPRI)
         END IF
         CALL WRTPRO(AOINT,NELMNT,LABINT(I),RTNLBL)
 600  CONTINUE
C     
      RETURN
      END
C/* Deck READ_CBIHER */
      SUBROUTINE READ_CBIHER(IPRDEF_HER, TSTINP_HER)
C
C Fetch two definitions from cbiher.h
C
#include <implicit.h>
#include <mxcent.h>
#include <cbiher.h>
      LOGICAL TSTINP_HER
      IPRDEF_HER = IPRDEF
      TSTINP_HER = TSTINP
      RETURN
      END
C /* Deck dktrf */
      SUBROUTINE DKTRF(STDER0,NBAST,NNBAST,IPRINT)
C
C     Read in Douglas-Kroll-Hess and overlap integrals and
C     transfer to standard arrays
C
#include <implicit.h>
#include <priunit.h>
      LOGICAL FNDLAB
#include <inftap.h>
      DIMENSION STDER0(NNBAST)
C
      REWIND (LUPROP)
      IF (FNDLAB('DKPOTINT',LUPROP)) THEN
         CALL READT(LUPROP,NNBAST,STDER0)
         IF (IPRINT .GT. 6) THEN
            CALL HEADER('Douglas-Kroll-Hess matrix read from LUPROP '//
     &                  'in DKTRF',-1)
            CALL OUTPAK(STDER0,NBAST,1,LUPRI)
         END IF
      ELSE
         WRITE (LUPRI,'(//3A)') ' >>> DKTRF error: Douglas-Kroll-Hess'//
     &      ' integrals not found on LUPROP.'
         CALL QUIT('DKTRF error: Integrals not found on LUPROP')
      END IF
      RETURN
      END
C  /* Deck dkpro */
      SUBROUTINE DKPRO(NONTYP,NONT,IQM,NBLCK,JCO,NUC,NRC,SEG,
     &           ALPHA,CPRIM,LPRIM,CONTC,KATOM,KANG,KBLOCK,KPRIM,
     &           IPRIMD,IORBD)
C******************************************************************************
C
C     Set up a contraction matrix for use in contraction of 
C     Douglas-Kroll-Hess integrals. Called from BASPRO(herrdn.F).
C
C******************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <pi.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4 = 4.0D0, DP5 = 0.50D0,
     &           DP75 = 0.75D0)
C
#include "cbirea.h"
      LOGICAL SEG,SPHER
      DIMENSION NONT(KATOM),IQM(KATOM),NBLCK(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          SEG(KBLOCK),ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK)
      DIMENSION CONTC(IPRIMD,IORBD)
      DIMENSION LPRIM(MAXPRD)
#include "ccom.h"
#include "nuclei.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "aosotr.h"
#ifdef PRG_DIRAC
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
C
C
      ICENT  = 0
      IBLOCK = 0
      IPR    = 1
      IRB    = 1
      CALL DZERO(CONTC,IPRIMD*IORBD)
      DO 10 I = 1,NONTYP
         DO 20 N = 1,NONT(I)
            ICENT = ICENT + 1
            KBCH = IBLOCK
            NDEG = NUCDEG(ICENT)
            DO 30 J = 1,IQM(I)
               KKK = 0
               NCOMP = KHK(J)
               NCCMP = KCK(J)
               SPHER = SPH(J)
               DO 40 K = 1, JCO(J,I)
                  KBCH = KBCH + 1
                  NCONT = NRC(KBCH)
C
C     Quick sort of primitive exponents if necessary
C
                  DO L1 = 1, NUC(KBCH)
                     LPRIM(L1) = L1
                  END DO
                  DO L1 = 1, NUC(KBCH) - 1
                     DO L2 = L1 + 1, NUC(KBCH)
                        IF(ALPHA(L2,KBCH) .GT. ALPHA(L1,KBCH)) THEN
                           IDUM = LPRIM(L2)
                           LPRIM(L2) = LPRIM(L1)
                           LPRIM(L1) = IDUM
                        END IF
                     END DO
                  END DO
C
C     Core for setting up contraction matrix
C
                  DO ICOMP = 1, NCOMP
                     IF (ICOMP .GT. 1) THEN
                        IPR = IPR - NDEG*(NCOMP*NUC(KBCH) - 1)
                        IRB = IRB - NDEG*NCONT*NCOMP + 1
                     END IF
                     DO LJ = 1, NUC(KBCH)
                        L = LPRIM(LJ)
                        IF (L .GT. 1) IRB = IRB - NCONT*NCOMP*NDEG
     &                       - (NDEG - 1)
                        PRFA  = ((D2*PI)**(DP75))
     &                       *((D4*ALPHA(L,KBCH))**(-DP5*(J + DP5)))
                        DO IDEG = 1, NDEG
                           IF (IDEG .GT. 1) IRB =IRB
     &                          -NCONT*NCOMP*NDEG + 1
                           DO M = 1, NCONT
                              CONTC(IPR,IRB) = CPRIM(L,M,KBCH)*PRFA
                              IRB = IRB + NCOMP*NDEG
                           END DO
                           IPR = IPR + 1
                        END DO
                        IPR = IPR + NDEG*NCOMP - NDEG
                     END DO
                  END DO
                  IPR = IPR - NCOMP*NDEG + NDEG
                  IRB = IRB - NCOMP*NDEG + 1
C
   40          CONTINUE
   30       CONTINUE
   20    CONTINUE
         IBLOCK = IBLOCK + NBLCK(I)
 10   CONTINUE
C
C     Write contraction matrix to file
C
      LUCNMT = -9000
      CALL GPOPEN(LUCNMT,'CNTMAT','NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      WRITE (LUCNMT,*) ((CONTC(I,J),I=1,IPRIMD),J=1,IORBD)
      CALL GPCLOSE(LUCNMT,'KEEP')
C
      RETURN
      END
C -- end of herdkh.F --
