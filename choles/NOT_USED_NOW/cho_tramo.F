      SUBROUTINE CHO_TRAMO(XONEIN,XTWOIN,WORK,LWORK)
C
C     asm  June 2005
C
C     Transform Cholesky AO vectors to MO basis
C
#include <implicit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <alftmp.h>
#include <choio.h>
#include <inftap.h>
C
C
      LOGICAL LDUMMY
C
      DIMENSION MMBST(8)
      DIMENSION XONEIN(*),XTWOIN(*)
      DIMENSION WORK(LWORK)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C
C---------------
C     Initialize
C---------------
C
      CALL CCSD_INIT1(WORK,LWORK)
C
      ICOUNT = 0
      NNBAST = 0
      DO ISYM = 1,NSYM
         ICMO(ISYM) = ICOUNT
         ICOUNT     = ICOUNT + NBAS(ISYM)*NORB(ISYM)
         NNBAST     = NNBAST + NBAS(ISYM)*(NBAS(ISYM)+1)/2
      END DO
C
      DO ISYM = 1,NSYM
C
         NNMOPQ(ISYM) = 0
         N2MOPQ(ISYM) = 0
C
         DO ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYM)
C
            N2MOPQ(ISYM) = N2MOPQ(ISYM) + NORB(ISYMA)*NORB(ISYMB)
C
            IF (ISYMB .GT. ISYMA) THEN
               NNMOPQ(ISYM) = NNMOPQ(ISYM) + NORB(ISYMA)*NORB(ISYMB)
            ELSE IF (ISYMB .EQ. ISYMA) THEN
               NNMOPQ(ISYM) = NNMOPQ(ISYM)
     &                      + NORB(ISYMA)*(NORB(ISYMA)+1)/2
            END IF
C
         END DO
      END DO
C
      ICOUNT = 0
      DO ISYM = 1,NSYM
         IMODIS(ISYM) = ICOUNT
         ICOUNT = ICOUNT + N2MOPQ(ISYM)*N2MOPQ(ISYM)
      END DO
      LE2MOB = ICOUNT
C
      DO ISYMAB = 1,NSYM
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMAB,ISYMB)
C
            IOFAP(ISYMA,ISYMB) = ICOUN1
            IOFPQ(ISYMA,ISYMB) = ICOUN2
            IPQPK(ISYMA,ISYMB) = ICOUN3
            IPQPK(ISYMB,ISYMA) = ICOUN3
C
            ICOUN1 = ICOUN1 + NBAS(ISYMA)*NORB(ISYMB)
            ICOUN2 = ICOUN2 + NORB(ISYMA)*NORB(ISYMB)
C
            IF (ISYMB .GT. ISYMA) THEN
               ICOUN3 = ICOUN3 + NORB(ISYMA)*NORB(ISYMB)
            ELSE IF (ISYMAB .EQ. 1) THEN
               ICOUN3 = ICOUN3 + NORB(ISYMB)*(NORB(ISYMB)+1)/2
            END IF
C
        END DO
      END DO
C
      DO ISYMAB = 1,NSYM
         DO ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)
C
            DO B = 1,NORB(ISYMB)
               DO A = 1,NORB(ISYMA)
C
                  NABSQ = IOFPQ(ISYMA,ISYMB) + NORB(ISYMA)*(B - 1) + A
C
                  IF (ISYMA .EQ. ISYMB) THEN
                     NAB = INDEX(A,B)
                  ELSE IF (ISYMB .GT.ISYMA) THEN
                     NAB = NORB(ISYMA)*(B - 1) + A
                  ELSE
                     NAB = NORB(ISYMB)*(A - 1) + B
                  END IF
C
                  NABPK = IPQPK(ISYMA,ISYMB) + NAB
C
                  IADMOP(NABPK) = NABSQ
C
               END DO
            END DO
         END DO
      END DO
C
C------------------
C     Cholesky info
C------------------
C
      LREDU = 1
      IF (REDUCE) THEN
         DO I = 1,NSYM
            IF (LREDU .LT. NNBST(I)) LREDU = NNBST(I)
         END DO
      END IF
C
      LURST = 0
      CALL GPOPEN(LURST,'CHOLESKY.RST','UNKNOWN','SEQUENTIAL',
     &            'UNFORMATTED',IDUM,.FALSE.)
C
      REWIND(LURST)
      READ(LURST) MSYM,THRCOM,LDUMMY,LDUMMY,LRDTOT,
     &            (MMBST(ISYM), ISYM=1,MSYM),
     &            (NUMCHO(ISYM), ISYM=1,MSYM),
     &            (IADRTO(ISYM), ISYM=1,MSYM),
     &            (NUMFIL(ISYM), ISYM=1,MSYM)
C
      READ(LURST) IDUM
C
      DO ISYM = 1,NSYM
         READ(LURST) (LENCHO(I,ISYM),I=1,NUMCHO(ISYM))
      END DO
      DO ISYM = 1,NSYM
         READ(LURST) (IDNTCH(I,ISYM),I=1,NUMCHO(ISYM))
      END DO
      DO ISYM = 1,NSYM
         READ(LURST) (IADRLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
      END DO
      DO ISYM = 1,NSYM
         READ(LURST) (IDNTLU(I,ISYM),I=0,NUMFIL(ISYM)-1)
      END DO
C
C
C-----------------------
C     Memory allocation
C-----------------------
C
      KIND = 1
      KCMO = KIND  + NSYM * ((LREDU-1)/IRAT + 1)
      KEND = KCMO  + NLAMDS
      LWRK = LWORK - KEND
      IF (LWRK .LT. 0) STOP 'Memory allocation failure.1'
C
C---------------------------
C     Get CMO and SCF energy
C---------------------------
C
      LUSIFC = 0
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',
     *            IDUMMY,.FALSE.)
C
      REWIND(LUSIFC)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
C
      REWIND(LUSIFC)
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,
     *              LSYM
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C---------------------------
C     One electron integrals
C---------------------------
C
      CALL CHO_TRAM0(XONEIN,WORK(KCMO),WORK(KEND),LWORK)
c     write(lupri,*) 'one electron integrals (spin orbital basis)'
c     ndim = 2*norbt
c     call output(xonein,1,ndim,1,ndim,ndim,ndim,1,lupri)
c     call flush(lupri)
C
C--------------
C     Transform
C--------------
C
      CALL CHO_TRAM1(XTWOIN,WORK(KIND),WORK(KCMO),
     &               WORK(KEND),LWRK,LURST)
C
C---------------
C     Write file
C---------------
C
      lumoin = 55
c     call gpopen(lumoin,'CHOMOVEC','UNKNOWN','UNFORMATTED',
c    &            IDUMMY,.FALSE.)
      
      open(lumoin,status='unknown',access='sequential',
     &     form='unformatted',file='CHOMOVEC')
C
      rewind (lumoin)
      write(lumoin) norbt,potnuc,emcscf
      ndim = norbt*norbt*4
      write(lumoin) (xonein(i), i=1,ndim)
      ndim = norbt*norbt*norbt*norbt*16
      write(lumoin) (xtwoin(i), i=1,ndim)
c
      write(lupri,*)
      write(lupri,*) 'CHOMOVEC written by cho_tramo'
      write(lupri,*)
c
c     call gpclose(lumoin,'keep')
C
      RETURN
      END
C
C=================================
C
      SUBROUTINE CHO_TRAM1(XMOINT,INDRED,CMO,WORK,LWORK,LURST)
C
C     asm  June 2005
C
C     Transform Cholesky AO vectors to MO basis
C
#include <implicit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
#include <alftmp.h>
C
C
      DIMENSION XMOINT(*),CMO(*)
      DIMENSION INDRED(*)
      DIMENSION WORK(LWORK)
C
C
C---------------
C     Get INDRED
C---------------
C
      REWIND(LURST)
      READ(LURST) IDUM
      READ(LURST) (INDRED(I), I=1,LRDTOT)
C
      CALL GPCLOSE(LURST,'KEEP')
C
C
      DO ISYMAB = 1,NSYM
C
         IF (NUMCHO(ISYMAB) .EQ. 0) GOTO 999
C
C-------------------------
C        Memory allocation
C        (Orbital basis)
C-------------------------
C
         LENFUL = N2BST(ISYMAB)
         LENPAC = NNBST(ISYMAB)
         LENRED = NREDUC(ISYMAB)
C
         KVECMO = 1
         KVECAO = KVECMO + LENFUL*4
         KVECHA = KVECAO + LENFUL
         KEND1  = KVECHA + LENFUL
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) STOP 'Memory allocation failure.2'
C
         LENVEC = LENPAC + LENRED
         NREAD  = MIN(NUMCHO(ISYMAB),LWRK1/LENVEC)
         IF (NREAD .LT. 1) STOP 'Memory allocation failure.3'
C
         KCHO  = KEND1
         KEND2 = KCHO  + NREAD*LENPAC
         LWRK2 = LWORK - KEND2
C
         NBATCH = (NUMCHO(ISYMAB)-1)/NREAD + 1
C
         ILAST = 0
         DO IBATCH = 1,NBATCH
C
C----------------------------------------------------------
C           Read and transform batch of AO Cholesky vectors
C----------------------------------------------------------
C
            IFIRST = ILAST + 1
            CALL CC_GETCHO(WORK(KCHO),NREAD,ISYMAB,IFIRST,INDRED,
     &                     WORK(KEND2),LWRK2)
            ILAST = ILAST + NREAD
            IF (ILAST .GT. NUMCHO(ISYMAB)) STOP 'Strange ILAST'
C
            IOPACK = 0
            CALL CHO_TRAM2(XMOINT,CMO,WORK(KCHO),ISYMAB,NREAD,
     &                     WORK(KVECMO),WORK(KVECAO),WORK(KVECHA),
     &                     WORK(KEND2),LWRK2,IOPACK)
C
         END DO
C
  999    CONTINUE
      END DO
C
C
      RETURN
      END
C
C=================================
C
      SUBROUTINE CHO_TRAM2(XMOINT,CMO,CHOVEC,ISYMAB,NREAD,VECMO,
     &                     VECAO,VECTMP,WORK,LWORK,IOPACK)
C
C     asm  June 2005
C
C     Transform Cholesky AO vectors to MO basis
C
C     If IOPACK .EQ. 0, vectors are kept in core in XMOINT
C     If IOPACK .EQ. 1, vectors are writen to disk squared up
C     If IOPACK .EQ. 2, vectors are writen to disk packed p .GE. q
C     N.B. Right now Cholesky vectors are in spatial orbital basis,
C     but if (IOPACK .EQ. 0), the full integral matrix is in spin-orbital
C
#include <implicit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
#include <alftmp.h>
C
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION XMOINT(*)
      DIMENSION CMO(*), CHOVEC(*), VECMO(*), VECAO(*), VECTMP(*)
      DIMENSION WORK(LWORK)
C
C
      IF (IOPACK .EQ. 0) THEN
         LENGTH = N2MOPQ(ISYMAB)
         CALL DZERO(VECMO,LENGTH)
         LENGTH = 0
      ELSE IF (IOPACK .EQ. 1) THEN
         LENGTH = N2MOPQ(ISYMAB)
      ELSE IF (IOPACK .EQ. 2) THEN
         LENGTH = NNMOPQ(ISYMAB)
      ELSE
         STOP 'Unknown IOPACK in CHO_TRAM2'
      END IF
C
      IF (LWORK .LT. LENGTH) STOP 'Memory allocation failure.4'
C
      LDUMP = 0
      KOFDU = 1
C
      DO I = 1,NREAD
C
C----------------------------
C        Square the vector up
C----------------------------
C
         KOFF = NNBST(ISYMAB)*(I-1) + 1
c        nx = norbt*(norbt+1)/2
c        write(lupri,*) 'before loop'
c        call flush(lupri)
c        do ni = 1,nx
c           write(lupri,*) 'chomo :',ni,chovec(koff+ni-1)
c           call flush(lupri)
c        enddo
         CALL CCSD_SYMSQ(CHOVEC(KOFF),ISYMAB,VECAO)
c        write(lupri,*) 'vector before transformation'
c        call output(vecao,1,norbt,1,norbt,norbt,norbt,
c    *               1,lupri)
C
C-------------------------------------
C        Transform each symmetry block
C-------------------------------------
C
         DO ISYMA = 1,NSYM
C
            ISYMB = MULD2H(ISYMA,ISYMAB)
            ISYMP = ISYMA
            ISYMQ = ISYMB
C
            NBASA = MAX(NBAS(ISYMA),1)
            NBASB = MAX(NBAS(ISYMB),1)
            NORBP = MAX(NORB(ISYMP),1)
C
            KOFF1 = IAODIS(ISYMA,ISYMB) + 1
            KOFF2 = ICMO(ISYMB) + 1
            KOFF3 = IOFAP(ISYMA,ISYMB) + 1
C
            CALL DGEMM('N','N',NBAS(ISYMA),NORB(ISYMB),NBAS(ISYMB),
     &                 ONE,VECAO(KOFF1),NBASA,CMO(KOFF2),NBASB,
     &                 ZERO,VECTMP(KOFF3),NBASA)
C
            KOFF1 = ICMO(ISYMA) + 1
            KOFF2 = IOFAP(ISYMA,ISYMB) + 1
            KOFF3 = IOFPQ(ISYMA,ISYMB) + 1
C
            CALL DGEMM('T','N',NORB(ISYMA),NORB(ISYMB),NBAS(ISYMA),
     &                 ONE,CMO(KOFF1),NBASA,VECTMP(KOFF2),NBASA,
     &                 ZERO,VECMO(KOFF3),NORBP)
C
         END DO
c        call output(vecmo,1,norbt,1,norbt,norbt,norbt,
c    *               1,lupri)
C
         IF (IOPACK. EQ. 0) THEN
            CALL CHO_MOINT(ISYMAB,VECMO,XMOINT,WORK,LWORK)
c           write(lupri,*) 'vector',isymab
c           call output(vecmo,1,2*norbt,1,2*norbt,2*norbt,2*norbt,
c    *                  1,lupri)
c           write(lupri,*) 
            GOTO 100
         END IF
C
C-----------------------------------
C        Pack the transformed vector
C-----------------------------------
C
         IF (IOPACK .EQ. 2) THEN
            DO II = 1,LENGTH
               KOFF = IADMOP(I) 
               VECAO(II) = VECMO(KOFF)
            END DO
         END IF
C
C----------------------------------------------
C        Move to work or dump to disk if needed
C----------------------------------------------
C
         LTEMP = LDUMP + LENGTH
C
         IF (LTEMP .GT. LWORK) THEN
            CALL CHO_DUMPQ(ISYMAB,WORK,LDUMP)
            KOFDU = 1
            LDUMP = LENGTH
         ELSE
            LDUMP = LTEMP
         END IF
C
         IF (IOPACK .EQ. 1) THEN
            CALL DCOPY(LENGTH,VECMO,1,WORK(KOFDU),1)
         ELSE
            CALL DCOPY(LENGTH,VECAO,1,WORK(KOFDU),1)
         END IF
         KOFDU = KOFDU + LENGTH
C
  100    CONTINUE
C
      END DO
C
C-----------------------------------
C     Dump to disk remaining vectors
C-----------------------------------
C
      IF (IOPACK .NE. 0) CALL CHO_DUMPQ(ISYMAB,WORK,LDUMP)
C
C
      RETURN
      END
C
C=================================
C
      SUBROUTINE CHO_DUMPQ(ISYMAB,CHOPQ,NDIM)
C
C     asm  June 2005
C
C     Transform Cholesky AO vectors to MO basis
C
#include <implicit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
#include <alftmp.h>
C
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION CHOPQ(NDIM)
C
      CHARACTER*11 NAME
      CHARACTER*10 ROOT
      CHARACTER*1  CSYM
      DATA ROOT    /'CHOLES_MO_'/
C
      INTEGER*8 JJADR(8), ITMP
      SAVE JJADR
      DATA JJADR /1,1,1,1,1,1,1,1/
C
C
      WRITE(CSYM,'(I1)') ISYMAB
      NAME = ROOT(1:10)//CSYM
      LUCHO = 70 + ISYMAB
      CALL WOPEN2U(LUCHO,NAME,64,0)
C
      ITMP = JJADR(ISYMAB)
      CALL PUTWA2U(LUCHO,NAME,CHOPQ,ITMP,NDIM)
      JJADR(ISYMAB) = ITMP + NDIM
C
      CALL WCLOSE2(LUCHO,NAME,'KEEP')
C
C
      RETURN
      END
C
C=================================
C
      SUBROUTINE CHO_MOINT(ISYMAB,VECMO,XMOINT,WORK,LWORK)
C
C     asm  June 2005
C
C     Transform Cholesky AO vectors to MO basis
C
#include <implicit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <priunit.h>
#include <alftmp.h>
C
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION XMOINT(*), VECMO(*)
      DIMENSION WORK(LWORK)
C
      LOGICAL FIRST
      SAVE FIRST
      DATA FIRST /.TRUE./
C
C
      IF (ISYMAB .NE. 1)STOP 'Construction of XMOINT only w/o symmetry'
C
c     LENGTH = 16 * LE2MOB
      LENGTH = 16 * norbt*norbt*norbt*norbt 
      IF (FIRST) THEN
         CALL DZERO(XMOINT,LENGTH)
         FIRST = .FALSE.
      END IF
C
C     Construct VECMO in spin orbit basis
C
      IF (ISYMAB .NE. 1) STOP 'No symmetry if spin-orbit'
C
      KORB = 1
      LENGTH = NORBT*NORBT
      CALL DCOPY(LENGTH,VECMO,1,WORK(KORB),1)
      CALL DZERO(VECMO,4*LENGTH)
C
C     Alpha-alpha part
C
      KOFF1 = KORB
      KOFF2 = 1
      DO I = 1,NORBT
         CALL DCOPY(NORBT,WORK(KOFF1),1,VECMO(KOFF2),1)
         KOFF1 = KOFF1 + NORBT
         KOFF2 = KOFF2 + 2*NORBT
      END DO
C
C     Beta-Beta part
C
      KOFF1 = KORB
      KOFF2 = 2*NORBT*NORBT + NORBT + 1
      DO I = 1,NORBT
         CALL DCOPY(NORBT,WORK(KOFF1),1,VECMO(KOFF2),1)
         KOFF1 = KOFF1 + NORBT
         KOFF2 = KOFF2 + 2*NORBT
      END DO
C
C     Build XMOINT in spin orbital basis
C     Does not work with symmetry !!!!
C
      NIJ = IMODIS(ISYMAB)
      DO I = 1,N2MOPQ(ISYMAB)*4
         DO  J = 1,N2MOPQ(ISYMAB)*4
C
             NIJ = NIJ + 1
             XMOINT(NIJ) = XMOINT(NIJ) + VECMO(I) * VECMO(J)
C
         END DO
      END DO
C
C
      RETURN
      END
C
C=================================
C
      SUBROUTINE CHO_TRAM0(XONEIN,CMO,WORK,LWORK)
C
C     asm  June 2005
C
C     Transform Cholesky AO vectors to MO basis
C
#include <implicit.h>
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <alftmp.h>
C
C
      LOGICAL LDUMMY
C
      DIMENSION XONEIN(*),CMO(*)
      DIMENSION WORK(LWORK)
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
C
      ISYMAB = 1
C
C----------------------------------------------
C     Read and square up one-electron integrals
C----------------------------------------------
C
      KOFSQ = 1
      KOFPK = KOFSQ + N2BAST
      KEND1 = KOFPK + NNBAST
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) STOP 'Memory allocation failure.0'
C
      CALL RDONEL('ONEHAMIL',.TRUE.,WORK(KOFPK),NNBAST)
      CALL CCSD_SYMSQ(WORK(KOFPK),ISYMAB,WORK(KOFSQ))
C
C
C----------------------------------
C     Transform each symmetry block
C----------------------------------
C
      KOFHA = KOFPK
      KEND2 = KOFHA + N2BAST
      LWRK2 = LWORK - KEND2
      IF (LWRK2 .LT. 0) STOP 'Memory allocation failure.0'
C
      DO ISYMA = 1,NSYM
C
         ISYMB = ISYMA
         ISYMP = ISYMA
         ISYMQ = ISYMB
C
         NBASA = MAX(NBAS(ISYMA),1)
         NBASB = MAX(NBAS(ISYMB),1)
         NORBP = MAX(NORB(ISYMP),1)
C
         KOFF1 = KOFSQ + IAODIS(ISYMA,ISYMB)
         KOFF2 = ICMO(ISYMB) + 1
         KOFF3 = KOFHA + IOFAP(ISYMA,ISYMB)
C
         CALL DGEMM('N','N',NBAS(ISYMA),NORB(ISYMB),NBAS(ISYMB),
     &              ONE,WORK(KOFF1),NBASA,CMO(KOFF2),NBASB,
     &              ZERO,WORK(KOFF3),NBASA)
C
         KOFF1 = ICMO(ISYMA) + 1
         KOFF2 = KOFHA + IOFAP(ISYMA,ISYMB)
         KOFF3 = IOFPQ(ISYMA,ISYMB) + 1
C
         CALL DGEMM('T','N',NORB(ISYMA),NORB(ISYMB),NBAS(ISYMA),
     &              ONE,CMO(KOFF1),NBASA,WORK(KOFF2),NBASA,
     &              ZERO,XONEIN(KOFF3),NORBP)
C
      END DO
c     write(lupri,*) 'One electron integral (orbital basis)'
c     call flshfo(lupri)
c     call output(xonein,1,norbt,1,norbt,norbt,norbt,1,lupri)
C
C     Construct XONEIN in spin orbit basis
C
      IF (ISYMAB .NE. 1) STOP 'No symmetry when constructing XONEINT'
C
      KORB = 1
      LENGTH = NORBT*NORBT
      CALL DCOPY(LENGTH,XONEIN,1,WORK(KORB),1)
      CALL DZERO(XONEIN,4*LENGTH)
C
C     Alpha-alpha part
C
      KOFF1 = KORB
      KOFF2 = 1
      DO I = 1,NORBT
         CALL DCOPY(NORBT,WORK(KOFF1),1,XONEIN(KOFF2),1)
         KOFF1 = KOFF1 + NORBT
         KOFF2 = KOFF2 + 2*NORBT
      END DO
C
C     Beta-Beta part
C
      KOFF1 = KORB
      KOFF2 = 2*NORBT*NORBT + NORBT + 1
      DO I = 1,NORBT
         CALL DCOPY(NORBT,WORK(KOFF1),1,XONEIN(KOFF2),1)
         KOFF1 = KOFF1 + NORBT
         KOFF2 = KOFF2 + 2*NORBT
      END DO
C
C
      RETURN
      END
C
C=================================
C
c     program test
c
C
C     asm  June 2005
C
C     Transform Cholesky AO vectors to MO basis
C
c#include <implicit.h>
c#include <maxorb.h>
c#include <ccdeco.h>
c#include <ccsdsym.h>
c#include <ccorb.h>
c#include <iratdef.h>
c#include <priunit.h>
c#include <alftmp.h>
c#include <choio.h>
C
c     parameter (morb = 11)
c     parameter (lwork = 1000000)
c
c     dimension xonein(morb,morb), xtwoin(4*morb*morb,morb*morb)
c     dimension work(lwork)
c
c     call cho_tramo(xonein,xtwoin,work,lwork)
c
c     ndim = morb*morb
c     call output(xtwoin,1,ndim,1,ndim,ndim,ndim,1,lupri)
c
c     stop
c     end
