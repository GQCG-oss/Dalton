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
C /* Deck dftdns */
      SUBROUTINE DFTDNS(DMAT,WORK,LWORK,IPRINT)
C
C     T. Helgaker Feb 01 
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <iratdef.h>
#include <inforb.h>
C 
      DIMENSION DMAT(NBAST,NBAST), WORK(LWORK)
C
#include <inftap.h>
#include <infvar.h>
#include <nuclei.h>
#include <dftcom.h>
C
      IF (NCMOT .GT. LWORK) CALL STOPIT('DFTDNS',' ',NCMOT,LWORK)
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) 
      CALL READI(LUSIFC,IRAT*NCMOT,WORK)
      JKEEP = JWOPSY
      JWOPSY = 1
      CALL FCKDEN(.TRUE.,.FALSE.,DMAT,DUMMY,WORK,DUMMY,WORK,LWORK)
      JWOPSY = JKEEP
      IF (IPRINT.GT.100) THEN
         CALL HEADER('AO density matrix in DFTDNS ',-1)
         CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C  /* Deck getden */
      SUBROUTINE GETDEN(DMAT,WORK,LWORK,NODC,NODV,IPRINT)
C
C     T. Helgaker Sep 1999
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (DP5 = 0.50D0, D1 = 1.0D0)
C
#include <inforb.h>
C
      LOGICAL NODC, NODV
      DIMENSION DMAT(NBAST,NBAST), WORK(LWORK)
C
      IF (IPRINT .GE. 5) CALL TITLER('Output from GETDEN','*',103)
C
      KDSO  = 1
      KDNS  = KDSO + NNBAST
      KLAST = KDNS + NNBASX
      LWRK  = LWORK - KLAST + 1
      IF (KLAST.GT.LWORK) CALL STOPIT('GETDEN','QDRDSO',KLAST,LWORK)
      CALL QDRDSO(WORK(KDSO),WORK(KLAST),LWRK,IPRINT,NODC,NODV)
      CALL QDRDAO(WORK(KDNS),WORK(KDSO),NBAST,IPRINT)
C
      IJ = 0
      DO 100 I = 1, NBAST
      DO 100 J = 1, I
        FAC = D1
        IF (I.NE.J) FAC = DP5
        DMAT(I,J) = FAC*WORK(KDNS+IJ)
        DMAT(J,I) = FAC*WORK(KDNS+IJ)
        IJ = IJ + 1
  100 CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Total density matrix from GETDEN',-1)
         CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C
C  /* Deck qdrdso */
      SUBROUTINE QDRDSO(DSO,WORK,LWORK,IPRINT,NODC,NODV)
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
      LOGICAL NODC, NODV
      DIMENSION DSO(NNBAST), WORK(LWORK)
C
      KCMO  = 1
      KDV   = KCMO + NCMOT
      KLAST = KDV  + NNASHX
      IF (KLAST .GT. LWORK) CALL STOPIT('QDRDSO',' ',KLAST,LWORK)
      CALL QDRDS1(WORK(KCMO),WORK(KDV),DSO,IPRINT,NODC,NODV)
      RETURN
      END
C  /* Deck qdrds1 */
      SUBROUTINE QDRDS1(CMO,DV,DSO,IPRINT,NODC,NODV)
C
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, TWO = 2.0D0, HALF = 0.5D0)
#include <maxorb.h>
#include <inforb.h>
      LOGICAL NODC, NODV
      DIMENSION CMO(*), DV(*), DSO(NNBAST)
#include <iratdef.h>
#include <mxcent.h>
#include <priunit.h>
C
#include <abainf.h>
#include <inftap.h>
      INTEGER R, S, RS, U, V, UV
C
C
C     ***** Read input from LUSIFC *****
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,
     *              NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT
      NSSHT  = NORBT - NOCCT
      CALL READI(LUSIFC,IRAT*NCMOT,CMO)
      READ (LUSIFC)
      IF (NASHT .GT. 0) THEN
         CALL READI(LUSIFC,IRAT*NNASHX,DV)
      ELSE
         READ (LUSIFC)
      END IF
C
C     ***** Print Section *****
C
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, '(//A/)') ' <<<<< SUBROUTINE QDRDS1 >>>>>'
         WRITE (LUPRI, '(A,8I5)') ' NISH ', (NISH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NASH ', (NASH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NOCC ', (NOCC(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NORB ', (NORB(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NBAS ', (NBAS(I),I = 1,NSYM)
         IF (IPRINT .GE. 10) THEN
            CALL HEADER('Occupied molecular orbitals',0)
            IEND = 0
            DO 1000 ISYM = 1,NSYM
               IF (NBAS(ISYM) .EQ. 0) GOTO 1000
               IF (NOCC(ISYM) .EQ. 0) GOTO 1100
               WRITE (LUPRI, '(//,A,I5,/)') ' Symmetry ', ISYM
               IENDI = IEND
               DO 1200 I = 1, NOCC(ISYM)
                  WRITE (LUPRI,'(/,A,I5,/)') ' Molecular orbital ', I
                  WRITE (LUPRI,'(6F12.6)') (CMO(IENDI+J),J=1,NBAS(ISYM))
                  IENDI = IENDI + NBAS(ISYM)
1200           CONTINUE
1100           CONTINUE
               IEND = IEND + NORB(ISYM)*NBAS(ISYM)
1000        CONTINUE
            CALL HEADER('Active density matrix (MO basis)',-1)
            CALL OUTPAK(DV,NASHT,1,LUPRI)
         END IF
      END IF
C
C     ***** Construct contravariant SO matrices *****
C
      ISEND = 0
      ICEND = 0
      DO 110 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
         NISHI = NISH(ISYM)
         NASHI = NASH(ISYM)
         IASHI = IASH(ISYM)
         NBASI = NBAS(ISYM)
         IF (NBASI .EQ. 0) GOTO 120
         IF (NOCC(ISYM) .EQ. 0) THEN
            CALL DZERO(DSO(ISEND+1),NNBAS(ISYM))
            GO TO 120
         END IF
         RS = 0
         DO 100 R = 1, NBASI
            DO 200 S = 1,R
               RS = RS + 1
C
               DTRS = D0
C
C              (I) Inactive contribution
C
               IF (NISHI .GT. 0) THEN
                  ICENDI = ICEND
                  DO 300 I = 1, NISHI
                     DTRS = DTRS + CMO(ICENDI+R)*CMO(ICENDI+S)
                     ICENDI = ICENDI + NBASI
  300             CONTINUE
                  DTRS = DTRS + DTRS
               END IF
               IF (NODC) DTRS = D0
C
C              (II) Active contribution
C
               IF (.NOT. NODV) THEN
                  IF (NASHI .GT. 0) THEN
                     UV = ((IASHI + 1)*(IASHI + 2))/2
                     IDVEND = ICEND + NISHI*NBASI
                     ICENDU = IDVEND
                     DO 400 U = 1,NASHI
                        ICENDV = IDVEND
                        DO 410 V = 1, U
                           DUV = DV(UV)
                           IF (ABS(DUV) .GT. D0) THEN
                              TEMP = CMO(ICENDU+R)*CMO(ICENDV+S)
                              IF (U .NE. V) TEMP = TEMP
     *                             + CMO(ICENDU+S)*CMO(ICENDV+R)
                              DTRS = DTRS + DUV*TEMP
                           END IF
                           UV = UV + 1
                           ICENDV = ICENDV + NBASI
  410                   CONTINUE
                        UV = UV + IASHI
                        ICENDU = ICENDU + NBASI
  400                CONTINUE
                  END IF
               END IF
               IF (R .NE. S) DTRS = DTRS + DTRS
               DSO(ISEND+RS) = DTRS
C
  200       CONTINUE
  100    CONTINUE
C
C        ***** Print Section *****
C
         IF (IPRINT .GE. 10) THEN
            WRITE (LUPRI,'(1X,A,I5)') ' Symmetry', ISYM
            CALL HEADER('Total density matrix (SO basis)',-1)
            CALL OUTPAK(DSO(ISEND+1),NBASI,1,LUPRI)
         END IF
120      CONTINUE
         ISEND = ISEND + (NBASI*(NBASI + 1))/2
         ICEND = ICEND + NORBI*NBASI
110   CONTINUE
      RETURN
      END
C  /* Deck qdrdao */
      SUBROUTINE QDRDAO(DENMAT,DSO,NBAST,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      DIMENSION DSO(*), DENMAT(*)
#include <shells.h>
#include <pincom.h>
#include <symmet.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 10) CALL HEADER('Subroutine QDRDAO',-1)
C
C     Loop over all irreps in molecule
C
      ISOFF = 0
      ISTR = 1
      NNBASX = NBAST*(NBAST + 1)/2
      CALL DZERO(DENMAT,NNBASX)
      DO 100 IREP = 0, MAXREP
         NORBI = NAOS(IREP+1)
         IF (NORBI .EQ. 0) GOTO 110
         DO 200 I = ISTR,ISTR + NORBI - 1
            IA   = IBTAND(IBTSHR(IPIND(I),16),65535)
            NA   = IBTAND(IBTSHR(IPIND(I), 8),  255)
            IOFF = KSTRT(IA)
            MULA = ISTBAO(IA)
            INDA = IOFF + NA
            DO 300 J = ISTR,I
               IB   = IBTAND(IBTSHR(IPIND(J),16),65535)
               NB   = IBTAND(IBTSHR(IPIND(J), 8),  255)
               JOFF   = KSTRT(IB)
               NHKTB  = NHKT(IB)
               KHKTB  = KHKT(IB)
               MULB   = ISTBAO(IB)
               MAB    = IBTOR(MULA,MULB)
               KAB    = IBTAND(MULA,MULB)
               HKAB   = FMULT(KAB)
               ISOFF  = ISOFF + 1
               DSYMIJ = DSO(ISOFF)
               INDB   = JOFF + NB - KHKTB
               DO 400 ISYMOP = 0, MAXOPR
                  IF (IBTAND(ISYMOP,MAB) .NE. 0) GOTO 400
                  INDB = INDB + KHKTB
C
C                 Weight and parity factor
C
                  FAC = HKAB*
     *                  PT(IBTAND(ISYMOP,IBTXOR(IREP,ISYMAO(NHKTB,NB))))
                  INDM = MAX(INDA,INDB)
                  IND  = (INDM*(INDM - 3))/2 + INDA + INDB
                  DENMAT(IND) = DENMAT(IND) + FAC*DSYMIJ
400            CONTINUE
300         CONTINUE
200      CONTINUE
110      CONTINUE
         ISTR = ISTR + NORBI
100   CONTINUE
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Total density matrix (sym. distinct AO basis)',-1)
         CALL OUTPAK(DENMAT,NBAST,1,LUPRI)
      END IF
      RETURN
      END
C /* Deck dftdnsab */  
      SUBROUTINE DFTDNSAB(DMATA,DMATB,WORK,LWORK,IPRINT)
C
C     DFTDNS adaptation for open shell case ZR
C
#include <implicit.h>
      DIMENSION DMATA(*), DMATB(*)
      DIMENSION WORK(LWORK)
      INTEGER IPRINT
#include <inforb.h>
#include <inftap.h>
#include <dummy.h>
#include <priunit.h>
      PARAMETER ( D1 = 1.0D0, DP5 =0.5D0)
      LOGICAL CLOSED
C
C used from infvar.h: JWOPSY
#include <infvar.h>
C
      CALL QENTER('DFTDNSAB')
C
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KCMO,NCMOT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDV,NNASHX,WORK,KFREE,LFREE)
C    
      CLOSED=LUSIFC.LT.0
      IF (CLOSED) CALL GPOPEN(
     &   LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND(LUSIFC)
      CALL MOLLAB(LBSIFC,LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      CALL READT(LUSIFC,NCMOT,WORK(KCMO))
      READ(LUSIFC) 
      CALL READT(LUSIFC,NNASHX,WORK(KUDV))
C
      JKEEP = JWOPSY
      JWOPSY = 1
      IF (NASHT.EQ.0) THEN
          CALL FCKDEN(.TRUE.,.FALSE.,DMATB,DUMMY,WORK(KCMO),DUMMY,
     &                WORK(KFREE),LFREE)
          CALL DZERO(DMATA,N2BASX)
      ELSE
          CALL FCKDEN(.TRUE.,.TRUE.,DMATB,DMATA,WORK(KCMO),WORK(KUDV),
     &                WORK(KFREE),LFREE)
      END IF
      JWOPSY = JKEEP
      CALL DSCAL(N2BASX,DP5,DMATB,1)
      CALL DAXPY(N2BASX,D1,DMATB,1,DMATA,1)
C
      IF (CLOSED) CALL GPCLOSE(LUSIFC,'KEEP')
C
      IF (IPRINT.GT.100) THEN
         CALL HEADER('AO alpha density matrix in DFTDNSAB ',-1)
         CALL OUTPUT(DMATA,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL HEADER('AO beta density matrix in DFTDNSAB ',-1)
         CALL OUTPUT(DMATB,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      END IF
C        
      CALL MEMREL('DFTDNSAB',WORK,1,1,KFREE,LFREE)
C
      CALL QEXIT('DFTDNSAB')
C      
      END

C -- end of dft_den.F --
