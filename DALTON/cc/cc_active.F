!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C
C  /* Deck cc_selact */
      SUBROUTINE CC_SELACT(WORK,LWORK)
C
C     Alfredo Sanchez de Meras 2008
C
C     Purpose: Prepare selected orbitals to cheat CC code
C
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxash.h"
#include "dummy.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "inftap.h"
#include "infinp.h"
#include "infind.h"
#include "ccsdsym.h"
#include "ccfro.h"
#include "center.h"
#include "nuclei.h"
#include "molde.h"
#include "iratdef.h"
C
      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_SELACT')
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, HALF = 0.5D0)
C
      LOGICAL ATISAC, PROCC
      DATA PROCC /.FALSE./
C
      CHARACTER*8 STARS(3)
      CHARACTER*8 LABCC, LABIP, LAB21, LABEN, LINE
C
      DATA STARS /'********','********','********'/
      DATA LABCC /'TRCCINT '/
      DATA LABIP /'SIR IPH '/
      DATA LAB21 /'NEWORB  '/
      DATA LABEN /'EODATA  '/
C
      DIMENSION NOCC(8)
      DIMENSION IATORB(MXACBS,8)
      DIMENSION IOFCMO(8), IOFDEN(8), JFOTYP(2), IOFDE2(8)
C
      DIMENSION WORK(LWORK)
C
      LOGICAL DONE
      SAVE DONE
      DATA DONE /.FALSE./
C
      LOGICAL DBGPRT
      DATA DBGPRT  /.FALSE./
c
      logical aoexist
C
C
      IF (DONE) RETURN                ! Run selection only once
      DONE = .TRUE.
C
      CALL QENTER('SELACT')
      CALL GETTIM(TSTART,WSTART)
C
      WRITE(LUPRI,'(///,9X,A,/,9X,A,/,9X,A,/)')
     &            '======================================',
     &            'Entering selection of orbitals section',
     &            '======================================'
      CALL FLSHFO(LUPRI)
C
C     Information from Sirius to CC code
C     ----------------------------------
C
      IF (LUSIFC .GT. 0) CALL GPCLOSE(LUSIFC,'KEEP')

      CALL GPOPEN(LUSIFC,FNSIFC,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUSIFC)
C
      CALL MOLLAB(LABCC,LUSIFC,LUPRI)
      READ(LUSIFC) NSYM, NORBT, NBAST, NLAM, (NOCC(I),I=1,NSYM),
     &           (NORB(I),I=1,NSYM), (NBAS(I),I=1,NSYM), POTNUC, ESCF
      READ(LUSIFC) (WORK(I), I =1,NORBT)
      READ(LUSIFC) (WORK(I), I=NORBT+1,NLAM+NORBT)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C     Print RHF orbitals
C
C
C     Some index arrays
C
      ICOUN1 = 0
      ICOUN2 = 0
      ICOUN3 = 0
      ICOUN4 = 0
      ICOUN5 = 0
      ICOUN6 = 0
      ICOUN7 = 0
      DO I = 1,NSYM
C
         NRHF(I) = NOCC(I) 
         NVIR(I) = NORB(I) - NOCC(I)
C
         IOFCMO(I) = ICOUN1
         IOFDEN(I) = ICOUN2
         IORB(I)   = ICOUN3
         IBAS(I)   = ICOUN6
         IOFDE2(I) = ICOUN7
C
         ICOUN1 = ICOUN1 + NORB(I)*NBAS(I)
         ICOUN2 = ICOUN2 + NBAS(I)*NBAS(I)
         ICOUN3 = ICOUN3 + NORB(I)
         ICOUN4 = ICOUN4 + NBAS(I) * (NBAS(I)+1) / 2
         ICOUN5 = ICOUN5 + NRHF(I)
         ICOUN6 = ICOUN6 + NBAS(I)
         ICOUN7 = ICOUN7 + NBAS(I)*NBAST
C
         NACTOC(I) = 0
         NACTVI(I) = 0
C
         DO IAT1 = 1,NUCIND
            NOCVEC(IAT1,I) = 0
            NVIVEC(IAT1,I) = 0
         END DO
C
      END DO
C
      NLAMDA = ICOUN1
      N2BAST = ICOUN2
      NNBAST = ICOUN4
      NRHFT  = ICOUN5
      N2BASX = NBAST*NBAST
C
      IF (NLAMDA .NE. NLAM) 
     &   CALL QUIT('Something is rotten in the state of Norway.1')
C
      IF (DBGPRT) THEN
C
         write(lupri,*) 'sirius info read'
         write(lupri,*) 'nsym :',nsym
         write(lupri,*) 'norbt:',norbt
         write(lupri,*) 'nbast:',nbast
         write(lupri,*) 'nlam :',nlam
         write(lupri,*) 'nocc :',(nocc(i),i=1,nsym)
         write(lupri,*) 'nbas :',(nbas(i),i=1,nsym)
         write(lupri,*) 'norb :',(norb(i),i=1,nsym)
         write(lupri,*) 'nrhf :',(nrhf(i),i=1,nsym)
         write(lupri,*) 'nvir :',(nvir(i),i=1,nsym)
         call flshfo(lupri)
c
         IEND = 0
         WRITE(LUPRI,'(A,/,A,/)') 'Orbital energies','----------------'
         DO ISYM = 1,NSYM
            ISTART = IEND + 1
            IEND   = IEND + NORB(ISYM)
            WRITE(LUPRI,'(/,A,I3)') 'Symmetry block :',ISYM
            WRITE(LUPRI,'(4(I4,F13.8,3X))')
     &            (I,WORK(I),I=ISTART,IEND)
         END DO
         WRITE(LUPRI,'(/,A,/)')
     &               'RHF orbitals read from SIRIUS in CC_ACTIVE'
         CALL PRORB(WORK(NORBT+1),PROCC,LUPRI)
C
      END IF
C
C     Main memory allocation
C     ----------------------
C
      LOCCHO = 0
      LVICHO = 0
      LOPCHO = 0
      DO ISYM = 1,NSYM
         LOCCHO = LOCCHO + NRHF(ISYM)*NBAS(ISYM)
         LVICHO = LVICHO + NVIR(ISYM)*NBAS(ISYM)
      END DO
C
      KEPSYL = 1
      KCMO   = KEPSYL + NORBT
      KOCDEN = KCMO   + NLAMDA
      KVIDEN = KOCDEN + N2BAST
      KFCKAO = KVIDEN + N2BAST
      IF (DOMC) THEN
         KH1AO = KFCKAO
      ELSE
         KH1AO = KFCKAO + N2BAST
      END IF
      KH1PK  = KH1AO  + N2BAST
      KOCCHO = KH1PK  + NNBAST
      KVICHO = KOCCHO + LOCCHO
      KEND1  = KVICHO + LVICHO
      LWRK1  = LWORK  - KEND1
      IF (LWRK1. LT. 0) CALL QUIT('Not enough space in CC_SELACT.1')
C
      KVICHL = KEND1
c
      IF (DBGPRT) THEN
         write(lupri,*) 'First memory allocation'
         write(lupri,*) 'KCMO   :', KCMO 
         write(lupri,*) 'KOCDEN :', KOCDEN
         write(lupri,*) 'KVIDEN :', KVIDEN
         write(lupri,*) 'KFCKAO :', KFCKAO
         write(lupri,*) 'KH1AO  :', KH1AO 
         write(lupri,*) 'KH1PK  :', KH1PK 
         write(lupri,*) 'KOCCHO :', KOCCHO
         write(lupri,*) 'KVICHO :', KVICHO
         write(lupri,*) 'KVICHL :', KVICHL
         call flshfo(lupri)
      END IF
C
      KOVLP = KH1AO
      KOVPK = KH1PK
C
C     Spread of basis set and cmo
C
      IF (NSYM .NE. 1) DOSPREAD = .FALSE.
      IF (DOSPREAD) THEN
         KTMP  = KEND1
         KEND2 = KEND1 + NBAST*NBAST
         LWRK2 = LWORK - KEND2
         IF (LWRK2 .LT. 0) CALL QUIT('Not enough space for spread')
         CALL DZERO(WORK(KTMP),NBAST*NBAST)
         DO I = 1,NBAST
            II = NBAST*(I-1) + I
            WORK(KTMP+II-1) = 1.0D0
         END DO
         WRITE(LUPRI,*) 'Basis set'
         CALL CHO_SPREAD(WORK(KTMP),NBAST,WORK(KEND2),LWRK2,0,.TRUE.)
         WRITE(LUPRI,*) 'CMO set'
         CALL CHO_SPREAD(WORK(KCMO),NORBT,WORK(KEND2),LWRK2,NRHFT,
     &                   .TRUE.)
      END IF
  
C     Densities
C     ---------
C
      IOPT = 0
      CALL CHO_DENSI(NSYM,NRHF,NBAS,WORK(KCMO),WORK(KOCDEN),
     &               NRHF,IOPT,IOFCMO,IOFDEN)
C
      IF (.NOT. DOMC) CALL DSCAL(N2BAST,TWO,WORK(KOCDEN),1)
C
      IOPT = 1
      CALL CHO_DENSI(NSYM,NVIR,NBAS,WORK(KCMO),WORK(KVIDEN),
     &               NOCC,IOPT,IOFCMO,IOFDEN)
C
      IF (.NOT. DOMC) THEN
C
C        AO Fock matrix
C        --------------
C
         KFCKSQ = KEND1
         KSCR   = KFCKSQ + NBAST*NBAST
         KEND2  = KSCR   + NBAST*NBAST
         LWRK2  = LWORK  - KEND2
         IF (LWRK2 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.2a')
C
         CALL DZERO(WORK(KSCR),NBAST*NBAST)
         CALL CHO_DENUPK(WORK(KOCDEN),WORK(KSCR),IOFDEN)
C
         IF (DBGPRT) THEN
            write(lupri,*) 'Occupied density matrix squared'
            call output(work(kscr),1,nbast,1,nbast,nbast,nbast,
     &                  1,lupri)
         END IF
C
         NDMAT  = 1
         ISYDEN = 1
         JFOTYP(1) = 13
         CALL DZERO(WORK(KFCKSQ),NBAST*NBAST)
         CALL SIRFCK(WORK(KFCKSQ),WORK(KSCR),NDMAT,ISYDEN,
     &               JFOTYP,DIRFCK,WORK(KEND2),LWRK2)
C
         IF (DBGPRT) THEN
            write(lupri,*) 'AO Fock matrix squared'
            call output(work(kfcksq),1,nbast,1,nbast,nbast,nbast,
     &                  1,lupri)
         END IF
C
         CALL CHO_FCKPK(WORK(KFCKSQ),WORK(KFCKAO),IOFDE2)
C
         IF (DBGPRT) THEN
            write(lupri,*) '2-electron AO Fock matrix packed'
            koff = kfckao
            do isym = 1,nsym
               write(lupri,*) 'Symmetry block :',isym
               call output(work(koff),1,nbas(isym),1,nbas(isym),
     &                     nbas(isym),nbas(isym),1,lupri)
               koff = koff + nbas(isym)*nbas(isym)
            end do
         END IF
C
         CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KH1PK),NNBAST)
         IF (DBGPRT) THEN
            WRITE(LUPRI,*) 'One-electron integrals paked'
            CALL OUTPKB(WORK(KH1PK),NBAS,NSYM,1,LUPRI)
         END IF
C
         KOFF1 = KH1PK
         KOFF2 = KH1AO
         DO ISYM = 1,NSYM
            CALL DSPTSI(NBAS(ISYM),WORK(KOFF1),WORK(KOFF2))
            IF (DBGPRT) THEN
               WRITE(LUPRI,*) 'One-electron integrals squared'
               WRITE(LUPRI,*) 'Symmetry block :',ISYM
               CALL OUTPUT(WORK(KOFF2),1,NBAS(ISYM),1,NBAS(ISYM),
     &                     NBAS(ISYM),NBAS(ISYM),1,LUPRI)
            END IF
            KOFF1 = KOFF1 + NBAS(ISYM)*(NBAS(ISYM)+1)/2
            KOFF2 = KOFF2 + NBAS(ISYM)*NBAS(ISYM)
         END DO
C
         CALL DAXPY(N2BAST,ONE,WORK(KH1AO),1,WORK(KFCKAO),1)
C
         IF (DBGPRT) THEN
            write(lupri,*)
            write(lupri,*) 'AO Fock matrix'
            koff = kfckao
            do isym = 1,nsym
               nbasi = nbas(isym)
               call output(work(koff),1,nbasi,1,nbasi,nbasi,nbasi,
     &                     1,lupri)
               koff = koff + nbasi*nbasi
            end do
         END IF
C
         CALL DSCAL(N2BAST,HALF,WORK(KOCDEN),1)
C
      END IF
C
C     Check orthonormality of the orbitals
C
      CALL DZERO(WORK(KOVPK),NNBAST)
      CALL DZERO(WORK(KOVLP),N2BAST)
      CALL RDONEL('OVERLAP',.TRUE.,WORK(KOVPK),NNBAST)
C
      KOFF1 = KOVPK
      KOFF2 = KOVLP
      DO ISYM = 1,NSYM
         CALL DSPTSI(NBAS(ISYM),WORK(KOFF1),WORK(KOFF2))
         KOFF1 = KOFF1 + NBAS(ISYM)*(NBAS(ISYM)+1)/2
         KOFF2 = KOFF2 + NBAS(ISYM)*NBAS(ISYM)
      END DO
C
      DO ISYM = 1,NSYM
C
         KSCR1 = KEND1 
         KSCR2 = KSCR1 + NBAS(ISYM)*NORB(ISYM)
         KEND2 = KSCR2 + NORB(ISYM)*NORB(ISYM)
         LWRK2 = LWORK - KEND2
         IF (LWRK2 .LT. 0) CALL QUIT('Not enough memory in CC_SELACT')
C
         NBASI = MAX(NBAS(ISYM),1)
         NORBI = MAX(NORB(ISYM),1)
C
         KOFF1 = KCMO  + IOFCMO(ISYM)
         KOFF2 = KOVLP + IOFDEN(ISYM)
         KOFF3 = KSCR1
         CALL DGEMM('T','N',NORB(ISYM),NBAS(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NBASI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NORBI)
C
         KOFF1 = KSCR1
         KOFF2 = KCMO + IOFCMO(ISYM)
         KOFF3 = KSCR2
         CALL DGEMM('N','N',NORB(ISYM),NORB(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NORBI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NORBI)
C
         IOVER = 0
         NORBI = NORB(ISYM)
         DO I = 1,NORBI
            KII = KSCR2 + NORBI*(I-1) + I - 1
            DIF = ABS(WORK(KII) - ONE)
            IF (DIF .GT. 1.0D-8) THEN
               WRITE(LUPRI,'(A,I4,A,I2,A,D18.10)') 
     &               'Norm of Hartree-Fock orbital',I,
     &               ' of symmetry',ISYM, ' is :',WORK(KII)
               IOVER = IOVER + 1
            END IF
C
            DO J = 1,I-1
               KIJ = KSCR2 + NORBI*(J-1) + I - 1
               KJI = KSCR2 + NORBI*(I-1) + J - 1
               DI1 = ABS(WORK(KIJ))
               DI2 = ABS(WORK(KIJ))
               IF ((DI1 .GT. 1.0D-8) .OR. (DI2 .GT. 1.0D-8)) THEN
                  WRITE(LUPRI,'(A,2I4,A,I2,2(A,D18.10))') 
     &                  'Overlap of Hartree-Fock orbitals',I,J,
     &                  ' of symmetry',ISYM, ' are :',WORK(KIJ),
     &                  ' and',WORK(KJI)
                  IOVER = IOVER + 1
               END IF
            END DO
         END DO
C
         IF (IOVER .NE. 0) THEN
            WRITE(LUPRI,*) 'Program stops because of previous errors'
            CALL QUIT('Hartree-Fock MOs not orthonormal')
         END IF
C
      END DO
C
C
C
C     Orbitals
C     ========
C
C     One of two possible methods:
C     1) Decompose first in the selected space and later the remaining space.
C        Code follows just below.
C     2) Decompose all in atom-by-atom basis and select afterwards.
C        Code after 9998 label
C     
C     Orthonormality is checked after 9999 label.
C
C
      IF (.NOT. SELDIR) GOTO 9998      
C
C     Get active functions
C
      IF (FULDEC) THEN
         NACATM = NUCIND
         DO ICENT = 1,NACATM
            LACTAT(ICENT) = ICENT
         END DO
      ELSE
         NACATM = NACINP
         DO ICENT = 1,NACATM
            LACTAT(ICENT) = LACINP(ICENT)
         END DO
      END IF
C
      CALL CHO_ACTIVE(.TRUE.)
C
C     First algorithm
C     ---------------
C

      IOFF1 = 1
      IOFF2 = 1
      IATDUM = 0
      DO ISYM = 1,NSYM
C
         NRED = NACTBS(ISYM)
         NDIM = NBAS(ISYM)
         IF (NDIM .EQ. 0) GOTO 9100
C
C        Occupied part
C     
         IF (NRHF(ISYM) .EQ. 0) GOTO 9110
C     
         NVECS = 0
         KSCR1 = KEND1
         KEND2 = KSCR1 + NBAS(ISYM)*NRHF(ISYM)
         LWRK2 = LWORK - KEND2
         IF (LWRK2 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.2b')
C        
         IF (NRED .EQ. 0) GOTO 9120
C
         IOPT  = 0
         THRS = THACOC
         MXVECL = NRHF(ISYM)
         IF (LIMLOC) THEN
            THRS = 1.0D-5
            MXVECL = MXOCC(ISYM)
            IF (MXVECL .LE. 0) MXVECL = NRHF(ISYM)
         END IF
         KOFF1 = KOCDEN + IOFDEN(ISYM)
         WRITE(LUPRI,'(//,A,I3,/)') 'ACTIVE OCCUPIED. SYMMETRY :',ISYM
         CALL CHO_DENDECO(SELDIR,WORK(KOFF1),NDIM,LACBAS(IOFF1),
     &                    NRED,WORK(KSCR1),THRS,NVECS,IOPT,
     &                    MXVECL,IATDUM,MINSPR,WORK(KEND2),LWRK2)

         NACTOC(ISYM) = NVECS
         IF (NVECS .EQ. NRHF(ISYM)) GOTO 9130
C
 9120    CONTINUE
C
         IOPT  = 1
         THRS = 1.0D-6
         KOFF1 = KOCDEN + IOFDEN(ISYM)
         KSCR2 = KSCR1 + NBAS(ISYM)*NVECS
         WRITE(LUPRI,'(//,A,I3,/)') 'INACTIVE OCCUPIED. SYMMETRY :',ISYM
         CALL CHO_DENDECO(SELDIR,WORK(KOFF1),NDIM,LACBAS(IOFF1),
     &                    NRED,WORK(KSCR2),THRS,NVECS,IOPT,
     &                    MXVECL,IATDUM,MINSPR,WORK(KEND2),LWRK2)
         NINAOC(ISYM) = NVECS
C
 9130    CONTINUE
         NTEST = NRHF(ISYM) - NACTOC(ISYM) - NINAOC(ISYM)
         IF (NTEST .NE. 0) THEN
            CALL QUIT('Something is rotten in the state of Norway.3a')
c        ELSE
c           WRITE(LUPRI,'(//,2A,3I5,//)') 'Occupied symmetry.',
c    &            ' Active, inactive :',
c    &            ISYM, NACTOC(ISYM), NINAOC(ISYM)
         END IF
C
         KOFF2  = KCMO + IOFCMO(ISYM)
         LENGTH = NINAOC(ISYM) * NBAS(ISYM)
         IF (LENGTH .NE. 0)
     &       CALL DCOPY(LENGTH,WORK(KSCR2),1,WORK(KOFF2),1)
         KOFF3 = KOFF2 + LENGTH
         LENGTH = NACTOC(ISYM) * NBAS(ISYM)
         IF (LENGTH .NE. 0)
     &       CALL DCOPY(LENGTH,WORK(KSCR1),1,WORK(KOFF3),1)
C
 9110    CONTINUE
C
C        Virtual part
C
         IF (NVIR(ISYM) .EQ. 0) GOTO 9100
C
C
         NVECS = 0
         KSCR1 = KEND1
         KEND2 = KSCR1 + NBAS(ISYM)*NVIR(ISYM)
         LWRK2 = LWORK - KEND2
         IF (LWRK2 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.2c')
C
         NRED = NACBSV(ISYM)
         IF (NRED .EQ. 0) GOTO 9140
C
         IOPT  = 0
         THRS = THACVI
         MXVECL = NVIR(ISYM)
         IF (LIMLOC) THEN
            THRS = 1.0D-5
            MXVECL = MXVIR(ISYM)
            IF (MXVECL .LE. 0) MXVECL = NVIR(ISYM)
         END IF
         KOFF1 = KVIDEN + IOFDEN(ISYM)
         WRITE(LUPRI,'(//,A,I3,/)') 'ACTIVE VIRTUAL. SYMMETRY :',ISYM
         CALL CHO_DENDECO(SELDIR,WORK(KOFF1),NDIM,LACBSV(IOFF2),
     &                    NRED,WORK(KSCR1),THRS,NVECS,IOPT,
     &                    MXVECL,IATDUM,MINSPR,WORK(KEND2),LWRK2)

         NACTVI(ISYM) = NVECS
         IF (NVECS .EQ. NVIR(ISYM)) GOTO 9150
C
 9140    CONTINUE
C
         IOPT  = 1
         THRS = 1.0D-6
         KOFF1 = KVIDEN + IOFDEN(ISYM)
         KSCR2 = KSCR1 + NBAS(ISYM)*NVECS
         WRITE(LUPRI,'(//,A,I3,/)') 'INACTIVE VIRTUAL. SYMMETRY :',ISYM
         CALL CHO_DENDECO(SELDIR,WORK(KOFF1),NDIM,LACBSV(IOFF2),
     &                    NRED,WORK(KSCR2),THRS,NVECS,IOPT,
     &                    MXVECL,IATDUM,MINSPR,WORK(KEND2),LWRK2)
         NINAVI(ISYM) = NVECS
C
 9150    CONTINUE
         NTEST = NVIR(ISYM) - NACTVI(ISYM) - NINAVI(ISYM)
         IF (NTEST .NE. 0) THEN
            CALL QUIT('Something is rotten in the state of Norway.4a')
c        ELSE
c           WRITE(LUPRI,'(//,2A,3I5,//)') 'Virtual symmetry.',
c    &            ' Active, inactive :',
c    &            ISYM, NACTVI(ISYM), NINAVI(ISYM)
         END IF
C
         KOFF2  = KCMO + NRHF(ISYM)*NBAS(ISYM) + IOFCMO(ISYM)
         LENGTH = NACTVI(ISYM) * NBAS(ISYM)
         IF (LENGTH .NE. 0)
     &       CALL DCOPY(LENGTH,WORK(KSCR1),1,WORK(KOFF2),1)
         KOFF3 = KOFF2 + LENGTH
         LENGTH = NINAVI(ISYM) * NBAS(ISYM)
         IF (LENGTH .NE. 0)
     &       CALL DCOPY(LENGTH,WORK(KSCR2),1,WORK(KOFF3),1)
C

 9100    CONTINUE
C
         IOFF1 = IOFF1 + NACTBS(ISYM)
         IOFF2 = IOFF2 + NACBSV(ISYM)
C
      END DO
C
      GOTO 9999            ! After second algorithm
C
 9998 CONTINUE             ! First algorithm has been skipped
C
C     The second algorithm
C     --------------------
C
      IF (NSYM .NE. 1) THEN
         WRITE(LUPRI,*) 'Decomposition in an atom-by-atom'
         WRITE(LUPRI,*) 'basis works only without symmetry'
         CALL QUIT('Atom-by-atom only if NSYM = 1')
      END IF
C
C     Occupied part
C
      WRITE(LUPRI,'(/,A,/,A,/)') 
     &            '>> Decomposition of occupied density <<',
     &            '======================================='
      IOFF = 1
      KORB = KOCCHO
      DO ISYM = 1,NSYM
C
         NDIM = NBAS(ISYM)
         MXVECL = NDIM      ! No limit of vectors 
C
         IF (NDIM .EQ. 0) GOTO 100
         IF (NRHF(ISYM) .EQ. 0) GOTO 100
C
         THRS = 4.0D0
  998    CONTINUE
C
         THRS = THRS / 5.0D0
         IF (THRS .LT. THACOC) THRS = THACOC
         WRITE(LUPRI,'(//,2A,D9.1,/)') 'Threshold for decomposition',
     &                ' of active occupied block of density',THRS
         DO IAT1 = 1,NUCIND
C
            NACATM = 1
            LACTAT(1) = IAT1     
            CALL CHO_ACTIVE(.FALSE.)
C
            NVECS = 0
            NRED  = NACTBS(ISYM)
            IF (NRED .EQ. 0) GOTO 110
C
            IOPT  = 0
            KOFF1 = KOCDEN + IOFDEN(ISYM)
            CALL CHO_DENDECO(SELDIR,WORK(KOFF1),NDIM,LACBAS(IOFF),
     &                       NRED,WORK(KORB),THRS,NVECS,IOPT,
     &                       MXVECL,IAT1,MINSPR,WORK(KEND1),LWRK1)
C
            NOCVEC(IAT1,ISYM) = NOCVEC(IAT1,ISYM) + NVECS
            DO IVEC = 1,NVECS
               IIVEC = NACTOC(ISYM) + IVEC
               IATORB(IIVEC,ISYM) = IAT1
            END DO
            NACTOC(ISYM) = NACTOC(ISYM) + NVECS
C
            KORB = KORB + NVECS*NBAS(ISYM)
C
  110       CONTINUE
         END DO                         !! Centers loop
C
         IF (THRS .GT. THACOC) THEN
            GOTO 998
         else
            write(lupri,*) 'Threshold for occupied ',
     &                     'decomposition satisfied'
         END IF
C
         IOFF = IOFF + NACTBS(ISYM)
  100    CONTINUE
C
         IF (NACTOC(ISYM) .NE. NRHF(ISYM)) 
     &       CALL QUIT('Something is rotten in the state of Norway.3b')
C
      END DO                            !! Symmetry loop
C
C     Reorder occupied orbitals according to atoms
C
      KORB = KOCCHO
      KNEW = KEND1
      LEN  = 0
      DO ISYM = 1, NSYM
         DO IAT1 = 1,NUCIND
            KOLD = KORB
            DO IVEC = 1,NACTOC(ISYM)
               IAT2 = IATORB(IVEC,ISYM)
               IF (IAT2 .EQ. IAT1) THEN
                  CALL DCOPY(NBAS(ISYM),WORK(KOLD),1,WORK(KNEW),1)
                  KNEW = KNEW + NBAS(ISYM)
                  LEN  = LEN  + NBAS(ISYM)
               END IF
               KOLD = KOLD + NBAS(ISYM)
            END DO
         END DO
         KORB = KORB + NACTOC(ISYM)*NBAS(ISYM)
      END DO
C
      KEND2 = KEND1 + LEN
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) CALL QUIT('Not enough space to resort occupied')
C
      IF (LEN .NE. LOCCHO) THEN
         WRITE(LUPRI,*) 'Error checking lenght of occupied orbitals'
         WRITE(LUPRI,*) 'Length allocated :',LOCCHO
         WRITE(LUPRI,*) 'Length computed  :',LEN
         CALL QUIT('Error in length')
       ELSE
         CALL DCOPY(LOCCHO,WORK(KEND1),1,WORK(KOCCHO),1)
       END IF
C
C
C     Virtual part
C
      WRITE(LUPRI,'(/,A,/,A,/)') 
     &            '>> Decomposition of virtual density <<',
     &            '======================================'
      IOFF = 1
      KORB = KVICHO
      DO ISYM = 1,NSYM
C
         NDIM = NBAS(ISYM)
         IF (NDIM .EQ. 0) GOTO 200
         IF (NVIR(ISYM) .EQ. 0) GOTO 200
C
         THRS = 50.0D0
  999    CONTINUE
C
         THRS = THRS / 5.0D0
         IF (THRS .LT. THACVI) THRS = THACVI
         WRITE(LUPRI,'(//,2A,D9.1,/)') 'Threshold for decomposition',
     &               ' of active virtual block of density',THRS
         DO IAT1 = NUCIND,1,-1
            NACATM = 1
            LACTAT(1) = IAT1
            CALL CHO_ACTIVE(.FALSE.)
C
            NVECS = 0
            NRED  = NACTBS(ISYM)
            IF (NRED .EQ. 0) GOTO 210
C
            IOPT  = 0
            KOFF1 = KVIDEN + IOFDEN(ISYM)
            CALL CHO_DENDECO(SELDIR,WORK(KOFF1),NDIM,LACBAS(IOFF),
     &                       NRED,WORK(KORB),THRS,NVECS,IOPT,
     &                       MXVECL,IAT1,MINSPR,WORK(KEND1),LWRK1)
C
            NVIVEC(IAT1,ISYM) = NVIVEC(IAT1,ISYM) + NVECS
            DO IVEC = 1,NVECS
               IIVEC = NACTVI(ISYM) + IVEC
               IATORB(IIVEC,ISYM) = IAT1
            END DO
            NACTVI(ISYM) = NACTVI(ISYM) + NVECS
C
            KORB = KORB + NVECS*NBAS(ISYM)
C
  210       CONTINUE
         END DO             !! Centers loop
c
         IF (THRS .GT. THACVI) THEN
            GOTO 999
         else
            write(lupri,*) 'Threshold for virtual ',
     &                     'decomposition satisfied'
         END IF
C
         IOFF = IOFF + NACTBS(ISYM)
  200    CONTINUE
C
         IF (NACTVI(ISYM) .NE. NVIR(ISYM)) 
     &       CALL QUIT('Something is rotten in the state of Norway.4b')
      END DO                !! Symmetry loop
C
      WRITE(LUPRI,'(//,A,/)')'Number of occupied orbitals in each atom'
      DO I = 1,NUCIND
         WRITE(LUPRI,'(I3,3X,8I6)') (I,NOCVEC(I,ISYM), ISYM=1,NSYM)
      END DO
      WRITE(LUPRI,'(//,A,/)') 'Number of virtual orbitals in each atom'
      DO I = 1,NUCIND
         WRITE(LUPRI,'(I3,3X,8I6)') (I,NVIVEC(I,ISYM), ISYM=1,NSYM)
      END DO
C
C     Reorder virtual orbitals according to atoms
C
      KORB = KVICHO
      KNEW = KEND1
      LEN  = 0
      DO ISYM = 1, NSYM
         DO IAT1 = NUCIND,1,-1
            KOLD = KORB
            DO IVEC = 1,NACTVI(ISYM)
               IAT2 = IATORB(IVEC,ISYM)
               IF (IAT2 .EQ. IAT1) THEN
                  CALL DCOPY(NBAS(ISYM),WORK(KOLD),1,WORK(KNEW),1)
                  KNEW = KNEW + NBAS(ISYM)
                  LEN  = LEN  + NBAS(ISYM)
               END IF
               KOLD = KOLD + NBAS(ISYM)
            END DO
         END DO
         KORB = KORB + NACTVI(ISYM)*NBAS(ISYM)
      END DO
C
      KEND2 = KEND1 + LEN
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) CALL QUIT('Not enough space to resort virtuals')
C
      IF (LEN .NE. LVICHO) THEN
         WRITE(LUPRI,*) 'Error checking lenght of virtual orbitals'
         WRITE(LUPRI,*) 'Length allocated :',LVICHO
         WRITE(LUPRI,*) 'Length computed  :',LEN
         CALL QUIT('Error in length')
       ELSE
         CALL DCOPY(LVICHO,WORK(KEND1),1,WORK(KVICHO),1)
       END IF
C
      WRITE(LUPRI,'(//,A,/)')'Results from decomposition of AO density'
      WRITE(LUPRI,'(A)') 'Sym     Occ     Virt    Tot'
      DO ISYM = 1,NSYM
         WRITE(LUPRI,'(I2,3I8)') ISYM, NACTOC(ISYM), NACTVI(ISYM),
     &                           NACTOC(ISYM) + NACTVI(ISYM)
      END DO
      WRITE(LUPRI,*)
      WRITE(LUPRI,*)
      CALL FLSHFO(LUPRI)
C
C     Sort orbitals according to selected atoms in input
C
      DO ISYM = 1,NSYM
C
         NACTOC(ISYM) = 0
         NACTVI(ISYM) = 0
C
         LACTOC = 0
         LACTVI = 0
         DO IAT1 = 1,NACINP
            IAT11  = LACINP(IAT1)
C
            NACTOC(ISYM) = NACTOC(ISYM) + NOCVEC(IAT11,ISYM)
            NACTVI(ISYM) = NACTVI(ISYM) + NVIVEC(IAT11,ISYM)
C
            LACTOC = LACTOC + NBAS(ISYM)*NOCVEC(IAT11,ISYM)
            LACTVI = LACTVI + NBAS(ISYM)*NVIVEC(IAT11,ISYM)
C
         END DO
         LINAOC = NBAS(ISYM)*NRHF(ISYM) - LACTOC
C
         NINAOC(ISYM) = NRHF(ISYM) - NACTOC(ISYM)
         NINAVI(ISYM) = NVIR(ISYM) - NACTVI(ISYM)
C
         KOINOC = KCMO   + IOFCMO(ISYM)
         KOACOC = KOINOC + LINAOC
         KOACVI = KOACOC + LACTOC
         KOINVI = KOACVI + LACTVI
C
         KFROMO = KOCCHO
         KFROMV = KVICHL
C
         NACINC = 0
         DO IAT1 = 1,NUCIND
C
            LOCAT1 = NBAS(ISYM)*NOCVEC(IAT1,ISYM)
            LVIAT1 = NBAS(ISYM)*NVIVEC(IAT1,ISYM)
C
            KFROMV = KFROMV - LVIAT1
C
            ATISAC = .FALSE.
C
            IF (NACINC .LT. NACINP) THEN
               DO IAT2 = 1,NACINP
                  IF (IAT1 .EQ. LACINP(IAT2)) THEN
                     ATISAC = .TRUE.
                     NACINC = NACINC + 1
                     GOTO 300
                  END IF
               END DO
            END IF
C
  300       CONTINUE
            IF (ATISAC) THEN
               CALL DCOPY(LOCAT1,WORK(KFROMO),1,WORK(KOACOC),1)
               CALL DCOPY(LVIAT1,WORK(KFROMV),1,WORK(KOACVI),1)
C
               KOACOC = KOACOC + LOCAT1
               KOACVI = KOACVI + LVIAT1
            ELSE
               CALL DCOPY(LOCAT1,WORK(KFROMO),1,WORK(KOINOC),1)
               CALL DCOPY(LVIAT1,WORK(KFROMV),1,WORK(KOINVI),1)
C
               KOINOC = KOINOC + LOCAT1
               KOINVI = KOINVI + LVIAT1
            END IF
C
            KFROMO = KFROMO + LOCAT1
C
         END DO             !! IAT1 loop 
      END DO                !! Symmetry loop
C
C     One way or another, we have now a selected set of Cholesky orbitals
 9999 CONTINUE
C
      WRITE(LUPRI,'(//,A)') ' '
C
      IF (DBGPRT) THEN
         WRITE(LUPRI,'(//,A,/)') 'Cholesky Molecular Orbitals'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,*) 'Symmetry block :',ISYM
            WRITE(LUPRI,*) 'Occupied active/inactive :',
     &                     NACTOC(ISYM), NINAOC(ISYM)
            WRITE(LUPRI,*) 'Virtual active/inactive :',
     &                     NACTVI(ISYM), NINAVI(ISYM)
            CALL OUTPUT(WORK(KCMO+IOFCMO(ISYM)),1,NBAS(ISYM),
     &                  1,NORB(ISYM),NBAS(ISYM),NORB(ISYM),1,LUPRI)
         END DO
      END IF
C
C     Orthonormality of Cholesky orbitals
C
      DO ISYM = 1,NSYM
C
         KSCR1 = KEND1
         KSCR2 = KSCR1 + NBAS(ISYM)*NORB(ISYM)
         KEND2 = KSCR2 + NORB(ISYM)*NORB(ISYM)
         LWRK2 = LWORK - KEND2
         IF (LWRK2 .LT. 0) CALL QUIT('Not enough memory in CC_SELACT.2')
C
         NBASI = MAX(NBAS(ISYM),1)
         NORBI = MAX(NORB(ISYM),1)
C
         KOFF1 = KCMO  + IOFCMO(ISYM)
         KOFF2 = KOVLP + IOFDEN(ISYM)
         KOFF3 = KSCR1
         CALL DGEMM('T','N',NORB(ISYM),NBAS(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NBASI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NORBI)
C
         KOFF1 = KSCR1
         KOFF2 = KCMO + IOFCMO(ISYM)
         KOFF3 = KSCR2
         CALL DGEMM('N','N',NORB(ISYM),NORB(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NORBI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NORBI)
C
         DNOMX = ZERO
         DOFMX = ZERO
         IOVER = 0
         NORBI = NORB(ISYM)
         DO I = 1,NORBI
            KII = KSCR2 + NORBI*(I-1) + I - 1
            DIF = ABS(WORK(KII) - ONE)
            IF (DIF .GT. 1.0D-6) THEN
               WRITE(LUPRI,'(A,I4,A,I2,A,D18.10)') 
     &               'Norm of Cholesky orbital',I,' of symmetry',
     &               ISYM, ' is :',WORK(KII)
               IOVER = IOVER + 1
            END IF
C
            IF (DIF .GT. DNOMX) DNOMX = DIF
C
            DO J = 1,I-1
               KIJ = KSCR2 + NORBI*(J-1) + I - 1
               KJI = KSCR2 + NORBI*(I-1) + J - 1
               DI1 = ABS(WORK(KIJ))
               DI2 = ABS(WORK(KIJ))
               IF ((DI1 .GT. 1.0D-6) .OR. (DI2 .GT. 1.0D-6)) THEN
                  WRITE(LUPRI,'(A,2I4,A,I2,2(A,D18.10))') 
     &                  'Overlap of Cholesky orbitals',I,J,
     &                  ' of symmetry',ISYM, ' are :',WORK(KIJ),
     &                  ' and',WORK(KJI)
                  IOVER = IOVER + 1
               END IF
C
               IF (DI1 .GT. DOFMX) DOFMX = DI1
               IF (DI2 .GT. DOFMX) DOFMX = DI2
C
            END DO
         END DO
C
         IF (IOVER .NE. 0) THEN
            WRITE(LUPRI,*) 'Program stops because of previous errors'
            CALL QUIT('Cholesky MOs not orthonormal')
         ELSE
           WRITE(LUPRI,'(A,I2)') 'Symmetry block',ISYM
           WRITE(LUPRI,'(2A,D12.6)') 'Maximum deviation from one ',
     &                            'in diagonal elements      :', DNOMX
           WRITE(LUPRI,'(2A,D12.6,/)') 'Maximum deviation from zero ',
     &                            'in off-diagonal elements :', DOFMX
         END IF
C
         IF (DBGPRT) THEN
            WRITE(LUPRI,*) 'Checking orthonormality of Cholesky CMO'
            WRITE(LUPRI,*) 'Symmetry block',ISYM
            CALL OUTPUT(WORK(KSCR2),1,NORB(ISYM),1,NORB(ISYM),
     &                  NORB(ISYM),NORB(ISYM),1,LUPRI)
         END IF
C
      END DO
C
C     Spread out of Cholesky orbitals 
C
      IF (DOSPREAD) THEN
         WRITE(LUPRI,*) 'Spread of Cholesky orbitals'
         CALL CHO_SPREAD(WORK(KCMO),NORBT,WORK(KEND1),LWRK1,NRHFT,
     &                  .TRUE. )
      END IF
C
C     In MCSCF we do not need orbital energies.
C
      IF (DOMC) THEN
C
C        Declare inactive orbitals as frozen
C
         LNOROT = .TRUE.
C
         DO ISYM = 1,NSYM
C
            NFRACT(ISYM) = NINAOC(ISYM)
C
            DO I = 1,NINAOC(ISYM)
               NOROT(I) = 1
            END DO
C
            I1 = IORB(ISYM) + NRHF(ISYM) + NACTVI(ISYM) + 1
            DO I = I1,NORB(ISYM)
               NOROT(I) = 1
            END DO
C
            IF (DBGPRT) THEN
               write(lupri,'(//,a)') 'Frozen orbitals'
               write(lupri,*) 'nfract :',nfract(isym)
               write(lupri,'(10(i6,i2))') (i,norot(i),i=1,norb(isym))
            END IF
         END DO
C
C        Write fake SIRIUS.RST
C
         LUF21 = -1
         CALL GPOPEN(LUF21,'SIRIUS.RST','UNKNOWN',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         REWIND(LUF21)
C
         NCMOT4 = MAX(4,NLAMDA)
         CALL NEWIT1
         WRITE(LUF21) STARS, LAB21
         CALL WRITT(LUF21,NCMOT4,WORK(KCMO))
C
         CALL GETDAT(STARS(2),STARS(3))
         WRITE (LUF21) STARS, LABEN
C
         CALL GPCLOSE(LUF21,'KEEP')
C
C        Print final results before exiting
C
         GOTO 500
C
      END IF
C
C     Orbital energies
C     ----------------
C
      KFKDIA = KEND1
      KEND2  = KFKDIA + NORBT
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.2d')
C
      DO ISYM = 1,NSYM
C
C        Inactive occupied part
C
         KFCKIJ = KEND2
         KSCRIQ = KFCKIJ + NINAOC(ISYM)*NINAOC(ISYM)
         KEND3  = KSCRIQ + NINAOC(ISYM)*NBAS(ISYM)
         LWRK3  = LWORK  - KEND3
         IF (LWRK3 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.3a')
C
         NBASI = MAX(NBAS(ISYM),1)
         NACTI = MAX(NINAOC(ISYM),1)
C
         KOFF1 = KCMO   + IOFCMO(ISYM)
         KOFF2 = KFCKAO + IOFDEN(ISYM)
         KOFF3 = KSCRIQ
         CALL DGEMM('T','N',NINAOC(ISYM),NBAS(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NBASI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NACTI)
C
         KOFF1 = KSCRIQ
         KOFF2 = KCMO   + IOFCMO(ISYM)
         KOFF3 = KFCKIJ
         CALL DGEMM('N','N',NINAOC(ISYM),NINAOC(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NACTI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NACTI)
C
         IF (DBGPRT) THEN
            write(lupri,*)
     &            'Fock matrix projected on inactive occupied orbitals'
            call output(work(kfckij),1,NINAOC(isym),1,NINAOC(isym),
     &                  NINAOC(isym),NINAOC(isym),1,lupri)
         END IF
C
C        Diagonalize
C
         KEVAL = KEND3
         KEVEC = KEVAL + NINAOC(ISYM)
         KTMP1 = KEVEC + NINAOC(ISYM)*NINAOC(ISYM)
         KTMP2 = KTMP1 + NINAOC(ISYM)
         KEND4 = KTMP2 + NINAOC(ISYM)
         LWRK4 = LWORK - KEND4
         IF (LWRK4 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.4a')
C
         NACTI = NINAOC(ISYM)
         MATZ  = 1
         CALL RS(NACTI,NACTI,WORK(KFCKIJ),WORK(KEVAL),MATZ,
     &           WORK(KEVEC),WORK(KTMP1),WORK(KTMP2),IERR)
         IF (IERR .NE. 0) CALL QUIT('Error diagonalizing occupied part')
C
C        Transform back and save
C
         KTRAN = KTMP1
         KEND5 = KTRAN + NINAOC(ISYM)*NBAS(ISYM)
         LWRK5 = LWORK - KEND5
         IF (LWRK5 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.5a')
C
         NBASI = MAX(NBAS(ISYM),1)
         NACTI = MAX(NINAOC(ISYM),1)
C
         KOFF = KCMO + IOFCMO(ISYM)
         CALL DGEMM('N','N',NBAS(ISYM),NINAOC(ISYM),NINAOC(ISYM),
     &              ONE,WORK(KOFF),NBASI,WORK(KEVEC),NACTI,
     &              ZERO,WORK(KTRAN),NBASI)
C
         LENGTH = NBAS(ISYM) * NINAOC(ISYM)
         IF (LENGTH .NE. 0)
     &       CALL DCOPY(LENGTH,WORK(KTRAN),1,WORK(KOFF),1)
C
         DO I = 1,NINAOC(ISYM)
            KOFF1 = KFKDIA + IORB(ISYM) + I - 1
            KOFF2 = KEVAL  + I - 1
            WORK(KOFF1) = WORK(KOFF2)
         END DO
C
C        Active occupied part
C
         KFCKIJ = KEND2
         KSCRIQ = KFCKIJ + NACTOC(ISYM)*NACTOC(ISYM)
         KEND3  = KSCRIQ + NACTOC(ISYM)*NBAS(ISYM)
         LWRK3  = LWORK  - KEND3
         IF (LWRK3 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.3b')
C
         NBASI = MAX(NBAS(ISYM),1)
         NACTI = MAX(NACTOC(ISYM),1)
C
         KOFF1 = KCMO   + IOFCMO(ISYM) + NBAS(ISYM)*NINAOC(ISYM)
         KOFF2 = KFCKAO + IOFDEN(ISYM)
         KOFF3 = KSCRIQ
         CALL DGEMM('T','N',NACTOC(ISYM),NBAS(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NBASI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NACTI)
C
         KOFF1 = KSCRIQ
         KOFF2 = KCMO   + IOFCMO(ISYM) + NBAS(ISYM)*NINAOC(ISYM)
         KOFF3 = KFCKIJ
         CALL DGEMM('N','N',NACTOC(ISYM),NACTOC(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NACTI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NACTI)
C
         IF (DBGPRT) THEN
            write(lupri,*) 
     &            'Fock matrix projected on active occupied orbitals'
            call output(work(kfckij),1,nactoc(isym),1,nactoc(isym),
     &                  nactoc(isym),nactoc(isym),1,lupri)
         END IF
C
C        Diagonalize
C
         KEVAL = KEND3
         KEVEC = KEVAL + NACTOC(ISYM)
         KTMP1 = KEVEC + NACTOC(ISYM)*NACTOC(ISYM)
         KTMP2 = KTMP1 + NACTOC(ISYM)
         KEND4 = KTMP2 + NACTOC(ISYM)
         LWRK4 = LWORK - KEND4
         IF (LWRK4 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.4b')
C
         NACTI = NACTOC(ISYM)
         MATZ  = 1
         CALL RS(NACTI,NACTI,WORK(KFCKIJ),WORK(KEVAL),MATZ,
     &           WORK(KEVEC),WORK(KTMP1),WORK(KTMP2),IERR)
         IF (IERR .NE. 0) CALL QUIT('Error diagonalizing occupied part')
C
C        Transform back and save
C
         KTRAN = KTMP1
         KEND5 = KTRAN + NACTOC(ISYM)*NBAS(ISYM)
         LWRK5 = LWORK - KEND5
         IF (LWRK5 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.5b')
C
         NBASI = MAX(NBAS(ISYM),1)
         NACTI = MAX(NACTOC(ISYM),1)
C
         KOFF = KCMO + IOFCMO(ISYM) + NINAOC(ISYM)*NBAS(ISYM)
         CALL DGEMM('N','N',NBAS(ISYM),NACTOC(ISYM),NACTOC(ISYM),
     &              ONE,WORK(KOFF),NBASI,WORK(KEVEC),NACTI,
     &              ZERO,WORK(KTRAN),NBASI)
C
         LENGTH = NBAS(ISYM) * NACTOC(ISYM)
         IF (LENGTH .NE. 0) 
     &       CALL DCOPY(LENGTH,WORK(KTRAN),1,WORK(KOFF),1)
C
         DO I = 1,NACTOC(ISYM)
            KOFF1 = KFKDIA + IORB(ISYM) + NINAOC(ISYM) + I - 1
            KOFF2 = KEVAL  + I - 1
            WORK(KOFF1) = WORK(KOFF2)
         END DO
C
         IF (DBGPRT) THEN
            WRITE(LUPRI,'(//,A,/)') 'Transformed Cholesky Occupied MOs'
            CALL OUTPUT(WORK(KCMO+IOFCMO(ISYM)),1,NBAS(ISYM),
     &                  1,NRHF(ISYM),NBAS(ISYM),NRHF(ISYM),1,LUPRI)
         END IF
C
C        Active virtual part
C
         KFCKAB = KEND2
         KSCRAQ = KFCKAB + NACTVI(ISYM)*NACTVI(ISYM)
         KEND3  = KSCRAQ + NACTVI(ISYM)*NBAS(ISYM)
         LWRK3  = LWORK  - KEND3
         IF (LWRK3 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.3c')
C
         NBASI = MAX(NBAS(ISYM),1)
         NACTI = MAX(NACTVI(ISYM),1)
C
         KOFF1 = KCMO   + IOFCMO(ISYM) + NRHF(ISYM)*NBAS(ISYM)
         KOFF2 = KFCKAO + IOFDEN(ISYM)
         KOFF3 = KSCRAQ
         CALL DGEMM('T','N',NACTVI(ISYM),NBAS(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NBASI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NACTI)
C
         KOFF1 = KSCRAQ
         KOFF2 = KCMO   + IOFCMO(ISYM) + NRHF(ISYM)*NBAS(ISYM)
         KOFF3 = KFCKAB
         CALL DGEMM('N','N',NACTVI(ISYM),NACTVI(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NACTI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NACTI)
C
C        Diagonalize
C
         KEVAL = KEND3
         KEVEC = KEVAL + NACTVI(ISYM)
         KTMP1 = KEVEC + NACTVI(ISYM)*NACTVI(ISYM)
         KTMP2 = KTMP1 + NACTVI(ISYM)
         KEND4 = KTMP2 + NACTVI(ISYM)
         LWRK4 = LWORK - KEND4
         IF (LWRK4 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.4c')
C
         NACTI = NACTVI(ISYM)
         MATZ  = 1
         CALL RS(NACTI,NACTI,WORK(KFCKAB),WORK(KEVAL),MATZ,
     &           WORK(KEVEC),WORK(KTMP1),WORK(KTMP2),IERR)
         IF (IERR .NE. 0) CALL QUIT('Error diagonalizing virtual part')
C
C        Transform back and save
C
         KTRAN = KTMP1
         KEND5 = KTRAN + NACTVI(ISYM)*NBAS(ISYM)
         LWRK5 = LWORK - KEND5
         IF (LWRK5 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.5c')
C
         NBASI = MAX(NBAS(ISYM),1)
         NACTI = MAX(NACTVI(ISYM),1)
C
         KOFF = KCMO + IOFCMO(ISYM) + NRHF(ISYM)*NBAS(ISYM)
         CALL DGEMM('N','N',NBAS(ISYM),NACTVI(ISYM),NACTVI(ISYM),
     &              ONE,WORK(KOFF),NBASI,WORK(KEVEC),NACTI,
     &              ZERO,WORK(KTRAN),NBASI)
C
         LENGTH = NBAS(ISYM) * NACTVI(ISYM)
         IF (LENGTH .NE. 0)
     &       CALL DCOPY(LENGTH,WORK(KTRAN),1,WORK(KOFF),1)
C
         DO I = 1,NACTVI(ISYM)
            KOFF1 = KFKDIA + IORB(ISYM) + NRHF(ISYM) + I - 1
            KOFF2 = KEVAL  + I - 1
            WORK(KOFF1) = WORK(KOFF2)
         END DO
C
C        Inactive virtual part
C
         KFCKAB = KEND2
         KSCRAQ = KFCKAB + NINAVI(ISYM)*NINAVI(ISYM)
         KEND3  = KSCRAQ + NINAVI(ISYM)*NBAS(ISYM)
         LWRK3  = LWORK  - KEND3
         IF (LWRK3 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.3d')
c
         NBASI = MAX(NBAS(ISYM),1)
         NACTI = MAX(NINAVI(ISYM),1)

         KOFF1 = KCMO   + IOFCMO(ISYM) + NRHF(ISYM)*NBAS(ISYM)
     &         + NACTVI(ISYM)*NBAS(ISYM)
         KOFF2 = KFCKAO + IOFDEN(ISYM)
         KOFF3 = KSCRAQ
         CALL DGEMM('T','N',NINAVI(ISYM),NBAS(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NBASI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NACTI)

         KOFF1 = KSCRAQ
         KOFF2 = KCMO   + IOFCMO(ISYM) + NRHF(ISYM)*NBAS(ISYM)
     &         + NACTVI(ISYM)*NBAS(ISYM)
         KOFF3 = KFCKAB
         CALL DGEMM('N','N',NINAVI(ISYM),NINAVI(ISYM),NBAS(ISYM),
     &              ONE,WORK(KOFF1),NACTI,WORK(KOFF2),NBASI,
     &              ZERO,WORK(KOFF3),NACTI)
C
c        Diagonalize
C
         KEVAL = KEND3
         KEVEC = KEVAL + NINAVI(ISYM)
         KTMP1 = KEVEC + NINAVI(ISYM)*NINAVI(ISYM)
         KTMP2 = KTMP1 + NINAVI(ISYM)
         KEND4 = KTMP2 + NINAVI(ISYM)
         LWRK4 = LWORK - KEND4
         IF (LWRK4 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.4d')

         NACTI = NINAVI(ISYM)
         MATZ  = 1
         CALL RS(NACTI,NACTI,WORK(KFCKAB),WORK(KEVAL),MATZ,
     &           WORK(KEVEC),WORK(KTMP1),WORK(KTMP2),IERR)
         IF (IERR .NE. 0) CALL QUIT('Error diagonalizing virtual part')
c
c        Transform back and save
c
         KTRAN = KTMP1
         KEND5 = KTRAN + NINAVI(ISYM)*NBAS(ISYM)
         LWRK5 = LWORK - KEND5
         IF (LWRK5 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.5d')

         NBASI = MAX(NBAS(ISYM),1)
         NACTI = MAX(NINAVI(ISYM),1)

         KOFF = KCMO + IOFCMO(ISYM) + NRHF(ISYM)*NBAS(ISYM)
     &         + NACTVI(ISYM)*NBAS(ISYM)
         CALL DGEMM('N','N',NBAS(ISYM),NINAVI(ISYM),NINAVI(ISYM),
     &              ONE,WORK(KOFF),NBASI,WORK(KEVEC),NACTI,
     &              ZERO,WORK(KTRAN),NBASI)

         LENGTH = NBAS(ISYM) * NINAVI(ISYM)
         IF (LENGTH .NE. 0)
     &       CALL DCOPY(LENGTH,WORK(KTRAN),1,WORK(KOFF),1)

         DO I = 1,NINAVI(ISYM)
            KOFF1 = KFKDIA + IORB(ISYM) + NRHF(ISYM)
     &            + NACTVI(ISYM) + I - 1
            KOFF2 = KEVAL  + I - 1
            WORK(KOFF1) = WORK(KOFF2)
         END DO

C
      END DO
C
  500 CONTINUE
C
      WRITE(LUPRI,'(//,A,/,A,//)') 
     &            'Final results of selection of orbitals',
     &            '--------------------------------------'
C
      WRITE(LUPRI,'(A,I4)') 'Number of active atoms :',NACINP
      WRITE(LUPRI,'(A)') 'Active atoms :'
      WRITE(LUPRI,'(16I5)') (LACINP(I), I=1,NACINP)
      WRITE(LUPRI,'(2A,2F6.2)')'Decomposition thresholds of ',
     &            'active occupied and virtual spaces :',
     &            THACOC, THACVI
      WRITE(LUPRI,*)
      WRITE(LUPRI,'(A,8I5)') 'Occupied inactive :',(NINAOC(I),I=1,NSYM)
      WRITE(LUPRI,'(A,8I5)') 'Occupied active   :',(NACTOC(I),I=1,NSYM)
      WRITE(LUPRI,'(A,8I5)') 'Virtual active    :',(NACTVI(I),I=1,NSYM)
      WRITE(LUPRI,'(A,8I5)') 'Virtual inactive  :',(NINAVI(I),I=1,NSYM)
C
      IF (.NOT. DOMC) THEN
C
         WRITE(LUPRI,'(/,A)') 'Orbital energies of Cholesky orbitals'
         KOFF0 = KFKDIA
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,2A,I3)') 'Inactive occupied block ',
     &                               'of symmetry :',ISYM
            KOFF1 = KOFF0
            WRITE(LUPRI,'(4(I3,A,F11.6,5X))')
     &                  (I,')',WORK(KOFF1+I-1),I=1,NINAOC(ISYM))
            WRITE(LUPRI,'(/,2A,I3)') 'Active occupied block ',
     &                               'of symmetry :',ISYM
            KOFF1 = KOFF0 + NINAOC(ISYM)
            WRITE(LUPRI,'(4(I3,A,F11.6,5X))')
     &                  (I,')',WORK(KOFF1+I-1),I=1,NACTOC(ISYM))
            WRITE(LUPRI,'(/,2A,I3)') 'Active virtual block ',
     &                               'of symmetry :',ISYM
            KOFF2 = KOFF0 + NRHF(ISYM)
            WRITE(LUPRI,'(4(I3,A,F11.6,5X))')
     &                  (I,')',WORK(KOFF2+I-1),I=1,NACTVI(ISYM))
            WRITE(LUPRI,'(/,2A,I3)') 'Inactive virtual block ',
     &                               'of symmetry :',ISYM
            KOFF2 = KOFF0 + NRHF(ISYM) + NACTVI(ISYM)
            WRITE(LUPRI,'(4(I3,A,F11.6,5X))')
     &                  (I,')',WORK(KOFF2+I-1),I=1,NINAVI(ISYM))
            KOFF0 = KOFF0 + NORB(ISYM)
         END DO
C
      END IF
C
C     Orbital spread of active orbitals
C
      NUMORB = NACTOC(1) + NACTVI(1)
      KOFF = KCMO + IOFCMO(1) + NINAOC(1)*NBAS(1)
      IF (DOMC) THEN 
         KEND2 = KEND1
         LWRK2 = LWRK1
      END IF
      IF (DOSPREAD) THEN
         WRITE(LUPRI,*) 'Spread of active Cholesky orbitals'
         CALL CHO_SPREAD(WORK(KOFF),NUMORB,WORK(KEND2),LWRK2,NACTOC(1),
     &                   .TRUE.)
      END IF
C
c     IF (NSYM .NE. 1) GOTO 666
C
C     Prepare MOLDEN2.INP
C
      DONEIV = .FALSE.
      DONEIU = .FALSE.
C
      LUMOSV = LUMOLDEN
      LUMOLDEN = -1
      CALL GPOPEN(LUMOLDEN,'MOLDEN2.INP','OLD',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
 
  555 CONTINUE
      READ(LUMOLDEN,'(A8)',END=556) LINE
      GOTO 555
  556 CONTINUE
#if defined (SYS_DARWIN)
      BACKSPACE(LUMOLDEN)
#endif
C
C_To do: Fix symmetry for molden
C
      DO I = 1,NORBT
         KOFF = KFKDIA + I - 1 
         OREN(I) = WORK(KOFF) 
      END DO
C
      KOCCU = KEND2
      KAOSO = KOCCU + NORBT
      KUNPK = KAOSO + N2BASX
      KOVEC = KUNPK + NORBT*NBAST
      KEND6 = KOVEC + NBAST
      LWRK6 = LWORK - KEND6
      IF (LWRK6 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.6')
C
      CALL DZERO(WORK(KOCCU),NORBT)
      DO I = 1,NRHFT
         WORK(KOCCU+I-1) = TWO
      END DO
C
      IF (DBGPRT) THEN
         WRITE(LUPRI,'(//,A,/)') 'Information passed to MOLDEN_MOS'
         DO I = 1,NORBT
            WRITE(LUPRI,'(A,F12.6)') 'Energy',OREN(I)
c    &                 ,'   Occupation',WORK(KOCCU+I-1)
         END DO
      END IF
C
      CALL MOLDEN_MOS(1,WORK(KCMO),WORK(KOCCU),WORK(KAOSO),
     &           WORK(KUNPK),WORK(KOVEC))
 
      CALL GPCLOSE(LUMOLDEN,'KEEP')
      LUMOLDEN = LUMOSV
C
  666 CONTINUE
C
      WRITE(LUPRI,'(//,A,/)') 'Active Cholesky Molecular Orbitals'
      CALL PRORB(WORK(KCMO),PROCC,LUPRI)
C
      IF (DOMC) GOTO 1000
C
C     Write fake SIRIFC
C
      IONE = 1
      IZER = 0
C
      CALL GPOPEN(LUSIFC,FNSIFC,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUSIFC)
C
      WRITE(LUSIFC) STARS,LABIP
      WRITE(LUSIFC) POTNUC, ESCF, ZERO, ESCF,
     &              IONE, IONE, IZER, IONE, IZER
!     WRITE(LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,
!    &              ISTATE,ISPIN,NACTEL,LSYM,MS2
C
      WRITE(LUSIFC) STARS, LABCC
      WRITE(LUSIFC) NSYM, NORBT, NBAST, NLAMDA, (NRHF(I),I=1,NSYM),
     &            (NORB(I),I=1,NSYM), (NBAS(I),I=1,NSYM), POTNUC, ESCF
      WRITE(LUSIFC) (WORK(KFKDIA+I-1), I=1,NORBT), (ISMO(I), I=1,NORBT),
     &            DUMMY, DUMMY, DUMMY
      WRITE(LUSIFC) (WORK(KCMO+I-1), I=1,NLAMDA)
C
      CALL GETDAT(STARS(2),STARS(3))
      WRITE(LUSIFC) STARS, LABEN
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C     Freeze inactive
C
      IF (FREEZE) FREEZE = .FALSE.
      IF (FROEXP) FROEXP = .FALSE.
      FROIMP = .TRUE.
C
      ICOUNT = 0
      NTOTAC = 0
      DO I = 1,NSYM
         NRHFFR(I) = NRHF(I) - NACTOC(I)
         NVIRFR(I) = NVIR(I) - NACTVI(I)     
C
         NTOTAC = NTOTAC + NACTOC(I)
C
         ICOUNT    = ICOUNT + NINAOC(I)
         IORACT(I) = ICOUNT
         ICOUNT    = ICOUNT + NACTOC(I) + NVIR(I)
      END DO
C
C     Freeze inside active space (only occupied orbitals)
C
      KFODI2 = KEND2
      KFOSYM = KFODI2 + NTOTAC
      KEND3  = KFOSYM + (NTOTAC - 1) / IRAT + 1
      LWRK3  = LWORK  - KEND3
      IF (LWRK3 .LT. 0) CALL QUIT('Not enough space in CC_SELACT.7')
C
      IF (ACTFRE) THEN
         IF (NSYM .EQ. 1) THEN
            NRHFFR(1) = NRHFFR(1) + NACTFR
         ELSE
            KOFF2 = KFODI2 - 1
            DO ISYM = 1,NSYM
               DO I = 1,NACTOC(ISYM)
                  KOFF1 = KFKDIA + IORACT(ISYM) + I - 1
                  KOFF2 = KOFF2  + 1
                  WORK(KOFF2) = WORK(KOFF1)
                END DO
             END DO
C
             CALL ACT_FREEZE(NSYM,WORK(KFODI2),WORK(KFOSYM),
     &                       NACTOC,NRHFFR,NACTFR)
C
         END IF
      END IF
C
       WRITE(LUPRI,*)
C
C     Over and done
C     -------------
C
 1000 CONTINUE
C
      CALL GETTIM(TEND,WEND)
      TCPU = TEND -TSTART
      TWAL = WEND -WSTART
      WRITE (LUPRI,'(//A,2F10.2,A//)')
     *      ' CPU and wall time for SELACT :',TCPU,TWAL,' seconds'

      WRITE(LUPRI,'(/,9X,A,/,9X,A,/,9X,A,///)')
     &            '=====================================',
     &            'Leaving selection of orbitals section',
     &            '====================================='
C
      CALL QEXIT('SELACT')
C
      RETURN
      END
C
C
      SUBROUTINE CHO_DENSI(NSYM,NVEC,NBAS,CMO,DMAT,NJUMP,
     &                     IOPT,IOFCMO,IOFDEN)
C
C     asm 2008
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include "priunit.h"
C
      DIMENSION CMO(*), DMAT(*)
      DIMENSION NVEC(8), NBAS(8), NJUMP(8)
      DIMENSION IOFCMO(8), IOFDEN(8)
C
      LOGICAL DBGPRT
      DATA DBGPRT /.FALSE./
C
C
      DO ISYM = 1,NSYM
C
         KOFF2 = IOFDEN(ISYM) + 1
         N2BSI = NBAS(ISYM) * NBAS(ISYM)
C
         IF (NBAS(ISYM) .EQ. 0) GOTO 100
         IF (NVEC(ISYM) .EQ. 0) THEN
            CALL DZERO(DMAT(KOFF2),N2BSI)
            GOTO 100
         END IF
C
         NBASI = MAX(NBAS(ISYM),1)
         IF (IOPT . EQ. 0) THEN
             KOFF1 = IOFCMO(ISYM) + 1
         ELSE IF (IOPT .EQ. 1) THEN
             KOFF1 = IOFCMO(ISYM) + NJUMP(ISYM)*NBAS(ISYM) + 1
         ELSE
             CALL QUIT('Unknown IOPT in CHO_DENSI')
         END IF
C
         CALL DGEMM('N','T',NBAS(ISYM),NBAS(ISYM),NVEC(ISYM),
     &               ONE,CMO(KOFF1),NBASI,CMO(KOFF1),NBASI,
     &               ZERO,DMAT(KOFF2),NBASI)
C
         IF (DBGPRT) THEN
            write(lupri,*)
            write(lupri,*)
            if (iopt. eq. 0) then
               write(lupri,*) 'AO Occupied density matrix. Sym:',isym
            else
               write(lupri,*) 'AO virtual density matrix. Sym:',isym
            end if 
            call output(dmat(koff2),1,nbas(isym),1,nbas(isym),
     &                  nbas(isym),nbas(isym),1,lupri)
         END IF
C
  100    CONTINUE
      END DO
C
      RETURN
      END
C
C
C  /* Deck cho_active */
      SUBROUTINE CHO_ACTIVE(INFO)
C
C     asm 2008
C
C     Purpose: Construct list of selected basis
C
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "center.h"
#include "symmet.h"
#include "ccorb.h"
#include "maxash.h"
#include "infind.h"
C
      LOGICAL INFO
C
C
      MSYM = MAXREP + 1
      IF (MSYM .NE.NSYM) 
     &    CALL QUIT('Something is rotten in the state of Norway.2')
C
      IF (DIALST) THEN
         DO ISYM = 1,NSYM
            NACTBS(ISYM) = 0
            NACBSV(ISYM) = 0
         END DO
         DO I = 1, NABSOC
            ISYM = ISAO(LACBAS(I))
            NACTBS(ISYM) = NACTBS(ISYM) + 1
         END DO
         DO I = 1, NABSVI
            ISYM = ISAO(LACBSV(I))
            NACBSV(ISYM) = NACBSV(ISYM) + 1
         END DO
         GOTO 100
      END IF
C
C     Include basis functions centered in active atoms
C
      NABSTO = 0
      DO ISYM = 1,NSYM
         NACTBS(ISYM) = 0
         DO ICENT = 1,NACATM
            IAT = LACTAT(ICENT)
            NACTBS(ISYM) = NACTBS(ISYM) + NBCENT(IAT,ISYM)
            IF (DIFADD) NACTBS(ISYM) = NACTBS(ISYM) + NEXTBS(ISYM)
         END DO
         NABSTO = NABSTO + NACTBS(ISYM) 
      END DO
C
      IF (NABSTO .GT. MXACBS) THEN
         WRITE(LUPRI,*) 'ERROR: From center.h, maximum ',
     &                  'number of active basis is',MXACBS
         CALL QUIT('Too many active basis')
      END IF
C
      ICOUN1 = 0
      DO ISYM = 1,NSYM
C
C        Include basis functions centered in active atoms
C
         DO ICENT = 1,NACATM
            IAT = LACTAT(ICENT)
            ICOUN2 = IBCENT(IAT,ISYM)
            DO I = 1,NBCENT(IAT,ISYM)
               ICOUN1 = ICOUN1 + 1
               ICOUN2 = ICOUN2 + 1
               LACBAS(ICOUN1) = ICOUN2
            END DO
         END DO
C
C        Add extra basis 
C
         IF (DIFADD) THEN
            DO I = 1,NEXTBS(ISYM)
               ICOUN1 = ICOUN1 + 1
               LACBAS(ICOUN1) = IEXTBS(I,ISYM)
            END DO
         END IF
C
      END DO
C
C     Error check
C
      IF (ICOUN1 .NE. NABSTO) THEN
         WRITE(LUPRI,'(//,A,//)') ' --- ERROR in CHO_ACTIVE'
         WRITE(LUPRI,'(A,I5)') 'Computed number of active basis',NABSTO
         WRITE(LUPRI,'(A,I5)') 'Dimension of LACBAS            ',ICOUN1
         CALL QUIT('Inconsistent number of active basis in CHO_ACTIVE')
      END IF
C
  100 CONTINUE   ! Come directly here if diagonal list given in input
C
      IF (INFO) THEN
         WRITE(LUPRI,'(A,/,A,/)') 'Information from CHO_ACTIVE',
     &                               '---------------------------'
         IF (.NOT. DIALST) THEN
            WRITE(LUPRI,'(I4,A)') NACATM, ' selected atoms :'
            WRITE(LUPRI,'(16I5)') (LACTAT(I), I=1,NACATM)
            WRITE(LUPRI,*)
            WRITE(LUPRI,'(/,I5,A)') NABSTO, ' selected basis :'
            WRITE(LUPRI,'(16I5)') (LACBAS(I), I=1,NABSTO)
        ELSE
            WRITE(LUPRI,'(I5,A)') NABSOC, 
     &                          ' selected basis for occupied space:'
            WRITE(LUPRI,'(16I5)') (LACBAS(I), I=1,NABSOC)
            WRITE(LUPRI,'(I5,A)') NABSVI, 
     &                          ' selected basis for virtual space:'
            WRITE(LUPRI,'(16I5)') (LACBSV(I), I=1,NABSVI)
        END IF
      END IF
C
C     Selected basis in virtual space are the same than in occupied 
C     space, unless if list of diagonals given in input
C
      IF (.NOT. DIALST) THEN
         KOFF = 0
         DO ISYM = 1,NSYM
            NACBSV(ISYM) = NACTBS(ISYM)
            DO I = 1, NACTBS(ISYM)
               LACBSV(KOFF+I) = LACBAS(KOFF+I)
            END DO
            KOFF = KOFF + NACTBS(ISYM)
         END DO
      END IF
C
      IF (NSYM .GT. 1) THEN
         KOFF1 = 0
         KOFF2 = 0
         DO ISYM = 1,NSYM
            DO I = 1,NACTBS(ISYM)
               LACBAS(KOFF1+I) = LACBAS(KOFF1+I) - IBAS(ISYM)
            END DO
            KOFF1 = KOFF1 + NACTBS(ISYM)
C
            DO I = 1,NACBSV(ISYM)
               LACBSV(KOFF2+I) = LACBSV(KOFF2+I) - IBAS(ISYM)
            END DO
            KOFF2 = KOFF2 + NACBSV(ISYM)
         END DO
      END IF
C
      RETURN
      END
C
C
C  /* Deck cho_dendeco */
      SUBROUTINE CHO_DENDECO(SELDIR,XMAT,NDIM,IVEC,NRED,XORB,THRS,NVECS,
     &                       IOPT,MXVECL,IAT1,MINSPR,WORK,LWORK)
C
C     Henrik Koch   3-may-2008
C
#include "implicit.h"
#include "infpri.h"
#include "priunit.h"
C
      LOGICAL MINSPR
      LOGICAL SELDIR
C
      DIMENSION XMAT(NDIM,NDIM),XORB(NDIM,*)
      DIMENSION WORK(LWORK)
      DIMENSION IVEC(NRED)
C
      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_DENDECO')
      PARAMETER (TOL = 1.00D-9)
C
C
      INEG = 0
      XSUM = 0.0D0
      DO I = 1,NDIM
C
         IF ((I .GT. MXVECL) . AND.  (IOPT .EQ. 0)) RETURN
C
         NVECS = I
C
         IF (.NOT. MINSPR) THEN
C
C           Select maximum diagonal 
C
            IF (IOPT .EQ. 0) THEN
               XMAX = -1.0D10
               DO J = 1,NRED
                  INDX = IVEC(J)
                  IF (ABS(XMAT(INDX,INDX)) .GT. XMAX) THEN
                     XMAX = XMAT(INDX,INDX)
                     IMAX = INDX
                  ENDIF
               ENDDO
            ELSE
               XMAX = -1.0D10
               DO J = 1,NDIM
                  IF (ABS(XMAT(J,J)) .GT. XMAX) THEN
                     XMAX = XMAT(J,J)
                     IMAX = J 
                  ENDIF
               ENDDO
            ENDIF
C
         ELSE
C
C           Select orbitals to give minimum spread
C
            KCMO = 1
            KEND = KCMO + NDIM*NDIM
            LWRK = LWORK - KEND
            IF (LWRK .LT. 0)
     &         CALL QUIT('Not enough memory in CHO_DENDECO')
            DO J = 1,NDIM
               XDIAG = SQRT(XMAT(J,J))
               IF (XDIAG .GT. THRS) THEN
                  DO K = 1,NDIM
                     KOFF = NDIM*(J-1) + K
                     WORK(KOFF) = XMAT(K,J) / XDIAG
                  END DO
               END IF
            END DO
C
            CALL CHO_SPREAD(WORK(KCMO),NDIM,WORK(KEND),LWRK,0,.true.)
C
            XMIN = 999.9D0
            IF (IOPT .EQ. 0) THEN
               DO J = 1,NRED
                  INDX = IVEC(J)
                  KOFF = KEND + INDX - 1
                  IF ((XMAT(INDX,INDX) .GT. THRS) .AND. 
     &            (WORK(KOFF) .LT. XMIN)) THEN
                     XMIN = WORK(KOFF)
                     IMIN = INDX
                  END IF
               END DO
            ELSE
               DO J = 1,NDIM
                  KOFF = KEND + J - 1
                  IF ((XMAT(J,J) .GT. THRS) .AND. 
     &            (WORK(KOFF) .LT. XMIN)) THEN
                     XMIN = WORK(KOFF)
                     IMIN = J 
                  END IF
               END DO
            END IF
C
            IF (XMIN .EQ. 999.9D0) THEN    ! All diagonal .lt. threshold
               NVECS = I - 1
               RETURN
            ELSE
               IMAX = IMIN
               XMAX = XMAT(IMAX,IMAX)
            END IF
C
         END IF
C 
         IF (XMAX .LT. 0.0D0) THEN
            IF (ABS(XMAX) .GT. TOL) THEN
               WRITE(LUPRI,'(A,A,A)')
     &         '*** ',SECNAM,': NOTICE ***'
               WRITE(LUPRI,'(A,1P,D30.15)')
     &         '    Negative diagonal encountered: ',XMAX
               WRITE(LUPRI,'(A,/)')
     &         '    - unable to continue!'
               CALL QUIT('ERROR IN '//SECNAM)
            ENDIF
         ENDIF
C
         IF (ABS(XMAX) .LT. THRS) THEN
            NVECS= I-1
            RETURN
         ENDIF
C
         DO J = 1,NDIM
            XORB(J,I) = XMAT(J,IMAX)/SQRT(XMAX)
         ENDDO
C
C        Get the norm of the Cholesky vector
C
         XNORM = DDOT(NDIM,XORB(1,I),1,XORB(1,I),1)
         XSUM  = XSUM + XNORM
         IF (SELDIR) THEN
            WRITE(LUPRI,'(A,I4,A,D15.7,A,I6,A,D15.7)') 
     &                  'Cholesky MO :', I, 
     &                  ' Diagonal',XMAX,' (',IMAX,') Norm',XNORM
         ELSE
            WRITE(LUPRI,'(A,I4,A,I5,A,D12.4,A,I6,A,D12.4)') 
     &                  'Cholesky MO :', I, '  Atom', IAT1,
     &                  ' Diagonal',XMAX,' (',IMAX,') Norm',XNORM
         END IF
C
         DO J = 1,NDIM
            DO K = 1,NDIM
               XMAT(K,J) = XMAT(K,J) - XORB(K,I)*XORB(J,I)
            ENDDO
         ENDDO
C
         DO J = 1,NDIM
            XMAT(J,IMAX) = 0.0D0
            XMAT(IMAX,J) = 0.0D0
         ENDDO
C
      ENDDO
C
      RETURN
      END
C
C
C /* Deck cho_denupk */
      SUBROUTINE CHO_DENUPK(DENPK,DENFUL,IOFDEN)
C
C     asm 2008
C
#include "implicit.h"
#include "inforb.h"
C
      DIMENSION DENPK(*), DENFUL(*)
      DIMENSION IOFDEN(8)
C
      INTEGER A, B
C
C
      DO ISYMB = 1,NSYM
         ISYMA = ISYMB
         DO B = 1,NBAS(ISYMB)
            IB = IBAS(ISYMB) + B
            DO A = 1,NBAS(ISYMA)
               IA = IBAS(ISYMA) + A
               KOFPK = IOFDEN(ISYMA) + NBAS(ISYMA)*(B-1) + A
               KOFFU = NBAST*(IB-1) + IA
               DENFUL(KOFFU) = DENPK(KOFPK)
            END DO
         END DO
      END DO
C
      RETURN
      END
C
C /* Deck cho_spread */
      SUBROUTINE CHO_SPREAD(CMO,NUMORB,WORK,LWORK,IOCVI,PRTSUM)
C
C     asm 2010
C
#include "implicit.h"
#include "priunit.h"
#include "inftap.h"
#include "inforb.h"
#include "mxcent.h"
#include "center.h"
C
      LOGICAL VIROCC,PRTSUM
C
      CHARACTER*8 LABEL
      CHARACTER*8 RTNLBL(2)
C
      DIMENSION CMO(*), WORK(LWORK)
C
C
      IF (NSYM .GT. 1) CALL QUIT('No symmetry implemented yet')
C
      VIROCC = .FALSE.
      IF (IOCVI .GT. 0) VIROCC = .TRUE.
C
      KSPRE = 1
      KINTX = KSPRE + NUMORB
      KINTY = KINTX + NBAST*NBAST
      KINTZ = KINTY + NBAST*NBAST
      KINXX = KINTZ + NBAST*NBAST
      KINYY = KINXX + NBAST*NBAST
      KINZZ = KINYY + NBAST*NBAST
      KREAD = KINZZ + NBAST*NBAST
      KEND1 = KREAD + NBAST*(NBAST+1)/2
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) CALL QUIT('Not enough space in CHO_SPREAD.1')
C
C     Read integrals
C
      LENINT = NBAST*(NBAST+1)/2
      LENSQ  = NBAST*NBAST
      IF (LUPROP .LE. 0) CALL GPOPEN(LUPROP,'AOPROPER','UNKNOWN',' ',
     &                               'UNFORMATTED',IDUMMY,.FALSE.)
C
C
      REWIND(LUPROP)
      LABEL = 'XDIPLEN'
      IERR = - 1
      CALL MOLLB2(LABEL,RTNLBL,LUPROP,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,*) 'Error when reading ', LABEL, ' integrals'
         WRITE(LUPRI,*) 'Calculation of orbital spread aborted'
         DOSPREAD = .FALSE.
         IF (MINSPR) THEN
            WRITE(LUPRI,*) 'Some required integrals not available',
     &                     ' in LUPROP.'
            WRITE(LUPRI,*) 'Diagonals will be selected',
     &                     ' according to their sizes'
            WRITE(LUPRI,*) ' and not to minimize orbital spreads'
            MINSPR = .FALSE.
         END IF
         RETURN
      END IF
C
      CALL DZERO(WORK(KREAD),LENINT)
      CALL READT(LUPROP,LENINT,WORK(KREAD))
      CALL DSPTGE(NBAST,WORK(KREAD),WORK(KINTX))
C
C
      REWIND(LUPROP)
      LABEL = 'YDIPLEN'
      IERR = - 1
      CALL MOLLB2(LABEL,RTNLBL,LUPROP,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,*) 'Error when reading ', LABEL, ' integrals'
         WRITE(LUPRI,*) 'Calculation of orbital spread aborted'
         DOSPREAD = .FALSE.
         IF (MINSPR) THEN
            WRITE(LUPRI,*) 'Some required integrals not available',
     &                     ' in LUPROP.'
            WRITE(LUPRI,*) 'Diagonals will be selected',
     &                     ' according to their sizes'
            WRITE(LUPRI,*) ' and not to minimize orbital spreads'
            MINSPR = .FALSE.
         END IF
         RETURN
      END IF
C
      CALL DZERO(WORK(KREAD),LENINT)
      CALL READT(LUPROP,LENINT,WORK(KREAD))
      CALL DSPTGE(NBAST,WORK(KREAD),WORK(KINTY))
C
C
      REWIND(LUPROP)
      LABEL = 'ZDIPLEN'
      IERR = - 1
      CALL MOLLB2(LABEL,RTNLBL,LUPROP,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,*) 'Error when reading ', LABEL, ' integrals'
         WRITE(LUPRI,*) 'Calculation of orbital spread aborted'
         DOSPREAD = .FALSE.
         IF (MINSPR) THEN
            WRITE(LUPRI,*) 'Some required integrals not available',
     &                     ' in LUPROP.'
            WRITE(LUPRI,*) 'Diagonals will be selected',
     &                     ' according to their sizes'
            WRITE(LUPRI,*) ' and not to minimize orbital spreads'
            MINSPR = .FALSE.
         END IF
         RETURN
      END IF
C
      CALL DZERO(WORK(KREAD),LENINT)
      CALL READT(LUPROP,LENINT,WORK(KREAD))
      CALL DSPTGE(NBAST,WORK(KREAD),WORK(KINTZ))
C
C
      REWIND(LUPROP)
      LABEL = 'XXSECMOM'
      IERR = - 1
      CALL MOLLB2(LABEL,RTNLBL,LUPROP,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,*) 'Error when reading ', LABEL, ' integrals'
         WRITE(LUPRI,*) 'Calculation of orbital spread aborted'
         DOSPREAD = .FALSE.
         IF (MINSPR) THEN
            WRITE(LUPRI,*) 'Some required integrals not available',
     &                     ' in LUPROP.'
            WRITE(LUPRI,*) 'Diagonals will be selected',
     &                     ' according to their sizes'
            WRITE(LUPRI,*) ' and not to minimize orbital spreads'
            MINSPR = .FALSE.
         END IF
         RETURN
      END IF
C
      CALL DZERO(WORK(KREAD),LENINT)
      CALL READT(LUPROP,LENINT,WORK(KREAD))
      CALL DSPTGE(NBAST,WORK(KREAD),WORK(KINXX))
C
C
      REWIND(LUPROP)
      LABEL = 'YYSECMOM'
      IERR = - 1
      CALL MOLLB2(LABEL,RTNLBL,LUPROP,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,*) 'Error when reading ', LABEL, ' integrals'
         WRITE(LUPRI,*) 'Calculation of orbital spread aborted'
         DOSPREAD = .FALSE.
         IF (MINSPR) THEN
            WRITE(LUPRI,*) 'Some required integrals not available',
     &                     ' in LUPROP.'
            WRITE(LUPRI,*) 'Diagonals will be selected',
     &                     ' according to their sizes'
            WRITE(LUPRI,*) ' and not to minimize orbital spreads'
            MINSPR = .FALSE.
         END IF
         RETURN
      END IF
C
      CALL DZERO(WORK(KREAD),LENINT)
      CALL READT(LUPROP,LENINT,WORK(KREAD))
      CALL DSPTGE(NBAST,WORK(KREAD),WORK(KINYY))
C
C
      REWIND(LUPROP)
      LABEL = 'ZZSECMOM'
      IERR = - 1
      CALL MOLLB2(LABEL,RTNLBL,LUPROP,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(LUPRI,*) 'Error when reading ', LABEL, ' integrals'
         WRITE(LUPRI,*) 'Calculation of orbital spread aborted'
         DOSPREAD = .FALSE.
         IF (MINSPR) THEN
            WRITE(LUPRI,*) 'Some required integrals not available',
     &                     ' in LUPROP.'
            WRITE(LUPRI,*) 'Diagonals will be selected',
     &                     ' according to their sizes'
            WRITE(LUPRI,*) ' and not to minimize orbital spreads'
            MINSPR = .FALSE.
         END IF
         RETURN
      END IF
C
      CALL DZERO(WORK(KREAD),LENINT)
      CALL READT(LUPROP,LENINT,WORK(KREAD))
      CALL DSPTGE(NBAST,WORK(KREAD),WORK(KINZZ))
C
C     Calculate orbitals spreads
C     SPREAD(P) = SQRT( <RR> - <R>^2 )
C
      KSCR  = KREAD
      KEND2 = KSCR  + NBAST
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) CALL QUIT('Not enough space in CHO_SPREAD.2')
C
      XMAX = 0.0D0
      XSUM = 0.0D0
      OMAX = 0.0D0
      OSUM = 0.0D0
      VMAX = 0.0D0
      VSUM = 0.0D0
C
      DO I = 1,NUMORB
C
         KOFF = 1 + NBAST * (I-1)
C
         XVAL  = AVEVAL(CMO(KOFF),WORK(KINTX),WORK(KSCR),NBAST)
         YVAL  = AVEVAL(CMO(KOFF),WORK(KINTY),WORK(KSCR),NBAST)
         ZVAL  = AVEVAL(CMO(KOFF),WORK(KINTZ),WORK(KSCR),NBAST)
         XXVAL = AVEVAL(CMO(KOFF),WORK(KINXX),WORK(KSCR),NBAST)
         YYVAL = AVEVAL(CMO(KOFF),WORK(KINYY),WORK(KSCR),NBAST)
         ZZVAL = AVEVAL(CMO(KOFF),WORK(KINZZ),WORK(KSCR),NBAST)
C
         WORK(KSPRE+I-1) = SQRT(ABS(XXVAL + YYVAL + ZZVAL 
     &                        - XVAL**2 - YVAL**2 - ZVAL**2))
C
         IF (WORK(KSPRE+I-1) .GT. XMAX) XMAX = WORK(KSPRE+I-1)
         XSUM = XSUM + WORK(KSPRE+I-1)
C
         IF (VIROCC) THEN
            IF (I .LE. IOCVI) THEN
               IF (WORK(KSPRE+I-1) .GT. OMAX) OMAX = WORK(KSPRE+I-1)
               OSUM = OSUM + WORK(KSPRE+I-1)
            ELSE
               IF (WORK(KSPRE+I-1) .GT. VMAX) VMAX = WORK(KSPRE+I-1)
               VSUM = VSUM + WORK(KSPRE+I-1)
            END IF
         END IF
      END DO
C
      XAVE = XSUM/DFLOAT(NUMORB)
      OAVE = OSUM/DFLOAT(IOCVI)
      VAVE = VSUM/DFLOAT(NUMORB-IOCVI)
C
      IF (.NOT. PRTSUM) RETURN
C
      WRITE(LUPRI,'(//,A,/)') 'Orbital spread'
      WRITE(LUPRI,'(4(I4,D12.4,4X))') (I,WORK(KSPRE+I-1),I=1,NUMORB)
      WRITE(LUPRI,'(/,A,D12.4)') 'Maximum spread :',XMAX
      IF (VIROCC) WRITE(LUPRI,'(A,D12.4)')   '      Occupied :',OMAX
      IF (VIROCC) WRITE(LUPRI,'(A,D12.4,/)')   '       Virtual :',VMAX
      WRITE(LUPRI,'(A,D12.4)') 'Average spread :',XAVE
      IF (VIROCC) WRITE(LUPRI,'(A,D12.4)')   '      Occupied :',OAVE
      IF (VIROCC) WRITE(LUPRI,'(A,D12.4)')   '       Virtual :',VAVE
      WRITE(LUPRI,*)
C
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION AVEVAL(VECTOR,XMAT,SCR,NDIM)
C
#include "implicit.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION VECTOR(NDIM),SCR(NDIM),XMAT(NDIM,NDIM)
C
C
      CALL DGEMV('N',NDIM,NDIM,ONE,XMAT,NDIM,VECTOR,1,ZERO,SCR,1)
      AVEVAL = DDOT(NDIM,VECTOR,1,SCR,1)
C
      RETURN
      END
c
      SUBROUTINE CHO_FCKPK(FCKSQ,FCKAO,IOFDE2)
C
C     asm. July 2011
C
C      Symmetry packed totally symmetric Fock matrix
C  
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
C
      DIMENSION FCKSQ(*), FCKAO(*)
      DIMENSION IOFDE2(*)
C
C
      KOFF1 = 1
      KOFF2 = 1
      DO ISYM = 1,NSYM
C
         KOFF1 = IOFDE2(ISYM) + IBAS(ISYM) + 1
C
         DO I = 1,NBAS(ISYM)
C
            LENGTH = NBAS(ISYM)
            CALL DCOPY(LENGTH,FCKSQ(KOFF1),1,FCKAO(KOFF2),1)
            KOFF1  = KOFF1 + NBAST 
            KOFF2  = KOFF2 + NBAS(ISYM)
C
         END DO
      END DO
C
      RETURN
      END
C
      SUBROUTINE ACT_FREEZE(NSYM,FCKDIA,ISYDIA,NACTOC,NRHFFR,NACTFR)
C
C  
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION FCKDIA(*)
      DIMENSION ISYDIA(*), NACTOC(*),NRHFFR(*)
C
      LOGICAL DBGPRT
      DATA DBGPRT  /.FALSE./
C
      NDIM = 0
      ICOU = 0
      DO ISYM = 1,NSYM
         NDIM = NDIM + NACTOC(ISYM)
         DO I = 1,NACTOC(ISYM)
            ICOU = ICOU + 1
            ISYDIA(ICOU) = ISYM
         END DO
      END DO
C
      IF (DBGPRT) THEN
         write(lupri,*) 'Output from ACT_FREEZE'
         WRITE(LUPRI,'(4(I3,A,F11.6,5X))')
     &                  (ISYDIA(I),')',FCKDIA(I),I=1,NDIM)
      END IF
C
C     Sort FCKDIA
C
      DO I = 1,NDIM
         FOCKI = FCKDIA(I)
         DO J = 1,I-1
            FOCKJ = FCKDIA(J)
            IF (FOCKI .LT. FOCKJ) THEN
               TMP = FCKDIA(I)
               FCKDIA(I) = FCKDIA(J)
               FCKDIA(J) = TMP
C
               ITMP = ISYDIA(I)
               ISYDIA(I) = ISYDIA(J)
               ISYDIA(J) = ITMP
            END IF
         END DO
      END DO
C
      IF (DBGPRT) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(4(I3,A,F11.6,5X))')
     &                 (ISYDIA(I),')',FCKDIA(I),I=1,NDIM)
      END IF
C
      WRITE(LUPRI,'(//,A,/)') ' >> Active orbitals frozen :'
      DO I = 1,NACTFR
         ISYM = ISYDIA(I)
         NRHFFR(ISYM) = NRHFFR(ISYM) + 1
         WRITE(LUPRI,'(A,I3,A,I2,A,F12.6)') 'Freezing orbital',I,
     &         ' with symmetry',ISYM,' and orbital energy',FCKDIA(I)
      END DO
C
      RETURN
       END
