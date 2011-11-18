C  /* Deck cc_chofop */
      SUBROUTINE CC_CHOFOP(WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Drive the calculation of CC2 first-order one-electron
C              properties with Cholesky decomposed two-electron integrals.
C
C     Overall structure is shamelessly copied from CC_FOP by A. Halkier.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <ccdeco.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <ccfop.h>
#include <chocc2.h>
#include <dummy.h>
#include <cclr.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHOFOP')

      INTEGER ORDER
      INTEGER ISYOF(8)

      LOGICAL RLORBS, LDUM

      CHARACTER*1  CDUM
      CHARACTER*3  LIST
      CHARACTER*3  APROXR12
      CHARACTER*10 MODEL,MODFIL

      LOGICAL DEBGER, LOCDBG
      PARAMETER (DEBGER = .FALSE., LOCDBG = .FALSE.)

C     Initializations.
C     ================

C     Check that this is a Cholesky run.
C     ----------------------------------

      IF (CHOINT) THEN
         WRITE(LUPRI,'(10X,A,/)')
     &   '(Calculation using Cholesky decomposed integrals)'
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,/)')
     &   'Error in ',SECNAM,':',
     &   'The two-electron integrals must be Cholesky decomposed!'
         CALL QUIT('Error in '//SECNAM)
      ENDIF

C     Print info, CC_FOP style.
C     -------------------------

      IF (MP2) THEN
         CALL AROUND('Coupled Cluster model is: MP2')
         MODEL  = 'MP2       '
         MODFIL = '          '
         WRITE(LUPRI,'(5X,A)')
     &   '- MP2 has not yet been implemented...returning.'
         RETURN
      ELSE IF (CC2) THEN
         CALL AROUND('Coupled Cluster model is: CC2')
         MODEL  = 'CC2       '
         MODFIL = 'CC2       '
      ELSE
         WRITE(LUPRI,'(//,5X,A,A,A,/,5X,A,/)')
     &   'Error in ',SECNAM,':',
     &   'Only model implemented is CC2'
         CALL QUIT('Error in '//SECNAM)
      ENDIF

      CALL FLSHFO(LUPRI)

C     Orbital relaxation is not implemented.
C     And it will not be until Cholesky RPA has been implemented...
C     -------------------------------------------------------------

      IF (RELORB) THEN
         WRITE(LUPRI,'(5X,A,A)')
     &   '*** WARNING: Cholesky CC2 relaxed first order properties not',
     &   ' implemented.'
         WRITE(LUPRI,'(5X,A)')
     &  '             Orbital relaxation switched off for Cholesky CC2.'
         RLORBS = RELORB
         RELORB = .FALSE.
      ENDIF

C     Hopefulle obsolete call to debugger.
C     ====================================

      IF (DEBGER) THEN
         CALL CC_CHOFOPDBG(WORK,LWORK,DEBGER)
      ENDIF

C     Solve equations for 0th order multipliers.
C     Solution vector stored on disk.
C     ==========================================

      NSTAT = 0
      ORDER = 0
      ISIDE = -1
      LIST  = 'L0 '

      ISYOF(1) = 0
      DO I = 2,NSYM
         ISYOF(I) = 1
      ENDDO

      IF (.NOT. L0SKIP) THEN
C
         APROXR12 = '   '
         CALL CC_SOLDRV(LIST,NSTAT,ORDER,ISIDE,APROXR12,
     *                  IDUMMY,IDUMMY,DUMMY,LDUM,
     *                  IDUMMY,CDUM,DUMMY,IDUMMY,
     *                  ISYOF,1,1,WORK,LWORK)
      ELSE
         WRITE(LUPRI,*) 'Skipped 0th-order multiplier (L0) equation'
         L0SKIP = .FALSE.
      ENDIF

C     Read in solution vectors.
C     -------------------------

      KDEN  = 1
      KTBAR = KDEN  + N2BST(1)
      KEND1 = KTBAR + NT1AM(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: sol. vec.'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      IOPT = 1
      CALL CC_RDRSP(LIST,0,1,IOPT,MODFIL,WORK(KTBAR),DUMMY)

      IF (LOCDBG) THEN
         TB1NRM = DDOT(NT1AM(1),WORK(KTBAR),1,WORK(KTBAR),1)
         WRITE(LUPRI,'(/,A,A,1P,D22.15,/,A,A,1P,D22.15,/)')
     &   SECNAM,': Sq. norm of TBAR1: ',TB1NRM,
     &   SECNAM,': 2---norm of TBAR1: ',DSQRT(TB1NRM)
      ENDIF

C     Calculate CC2 density matrix.
C     =============================

      CALL CC_CHODEN(WORK(KTBAR),WORK(KDEN),WORK(KEND1),LWRK1)

C     Calculate FOPs.
C     ===============

      KEND1 = KTBAR
      LWRK1 = LWORK - KEND1 + 1
      CALL CC_CHOFOP1(WORK(KDEN),WORK(KEND1),LWRK1,MODEL)

C     Restore RELORB.
C     ---------------

      IF (CC2) RELORB = RLORBS

      RETURN
      END
C  /* Deck cc_chofop1 */
      SUBROUTINE CC_CHOFOP1(DENAO,WORK,LWORK,MODEL)
C
C     Thomas Bondo Pedersen, February 2003.
C     Almost exact copy of the analogous CC_FOP section by A. Halkier.
C     This seperate routine has been written merely to avoid messing
C     around too much with allocations in CC_FOP [no tiene gracia, pero
C     funciona!]
C
C     Purpose: Calculate first order properties from the AO density
C              matrix.
C
#include <implicit.h>
      DIMENSION DENAO(*), WORK(LWORK)
      CHARACTER*10 MODEL
#include <mxcent.h>
#include <maxorb.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfop.h>
#include <cclr.h>
#include <ccrspprp.h>
#include <ccsections.h>
#include <ccroper.h>
#include <ccfield.h>
#include <dipole.h>
#include <quadru.h>
#include <nqcc.h>
#include <ccfdgeo.h>

      DIMENSION ELSEMO(3,3), SKODE(3,3), SKODN(3,3)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHOFOP1')

      PARAMETER (XMONE = -1.0D0, HALF = 0.50D0, ZERO = 0.0D0)
      PARAMETER (ONE = 1.0D0)

      CHARACTER*8 LABEL, LBLDP(3), LBLQD(6), LBLSM(6)
      CHARACTER*3 LIST

      DATA LBLDP /'XDIPLEN ', 'YDIPLEN ', 'ZDIPLEN '/
      DATA LBLQD /'XXTHETA ', 'XYTHETA ', 'XZTHETA ',
     &            'YYTHETA ', 'YZTHETA ', 'ZZTHETA '/
      DATA LBLSM /'XXSECMOM', 'XYSECMOM', 'XZSECMOM',
     &            'YYSECMOM', 'YZSECMOM', 'ZZSECMOM'/

C     Initialize.
C     -----------

      TIMTOT = SECOND()
      ICOUNT = 0

C     Allocation for property integrals.
C     ----------------------------------

      KONEP = 1
      KEND1 = KONEP + N2BASX
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',KEND1-1,
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Electric dipole moment.
C     =======================

      IF (DIPMOM) THEN

         ICOUNT = ICOUNT + 1

         CALL AROUND('Electric Dipole Moment')

C        Calculate nuclear contribution.
C        -------------------------------

         IPRLOC = IPRINT - 9
         CALL DIPNUC(WORK(KEND1),WORK(KEND1),IPRLOC,.FALSE.)

C        Calculate electronic contribution.
C        ----------------------------------

         DO IDIP = 1,3

            CALL DZERO(WORK(KONEP),N2BASX)
            FF   = ONE
            ISYM = -1
            CALL CC_ONEP(WORK(KONEP),WORK(KEND1),LWRK1,FF,ISYM,
     &                   LBLDP(IDIP))

            IF (ISYM .EQ. 1) THEN
               DIPME(IDIP) = -DDOT(N2BST(ISYM),WORK(KONEP),1,DENAO,1)
            ELSE
               DIPME(IDIP) = ZERO
            ENDIF

            DIPMN(IDIP) = DIPMN(IDIP) + DIPME(IDIP)

         ENDDO

C        Print.
C        ------

         IF (IPRLOC .GT. 0) THEN
            CALL HEADER('Electronic Contribution to Dipole Moment',-1)
            CALL DP0PRI(DIPME)
         ENDIF
         CALL HEADER('Total Molecular Dipole Moment',-1)
         CALL DP0PRI(DIPMN)

         CALL FLSHFO(LUPRI)

      ENDIF

C     Electric quadrupole moment.
C     ===========================

      IF (QUADRU) THEN

         ICOUNT = ICOUNT + 1

         CALL AROUND('Electric Quadrupole Moment (Buckingham)')

C        Calculate nuclear contribution.
C        -------------------------------

         IOPT   = 1
         IPRLOC = -1
         CALL CCNUCQUA(WORK(KEND1),LWRK1,IOPT,IPRLOC)
         CALL DZERO(QDREL,9)

C        Calculate electronic contribution.
C        ----------------------------------

         IJ = 0
         DO I = 1,3
            DO J = I,3

               IJ = IJ + 1

               CALL DZERO(WORK(KONEP),N2BASX)
               FF   = ONE
               ISYM = -1
               CALL CC_ONEP(WORK(KONEP),WORK(KEND1),LWRK1,FF,ISYM,
     &                      LBLQD(IJ))

               IF (ISYM .EQ. 1) THEN
                  LENAB = N2BST(ISYM)
                  CALL CCELQUA(WORK(KONEP),DENAO,LENAB,I,J,QDREL)
               ENDIF

            ENDDO
         ENDDO

C        Reorder.
C        --------

         CALL CC_QUAREO(QDREL,SKODE)
         CALL CC_QUAREO(QDRNUC,SKODN)

C        Print.
C        ------

         CALL DSCAL(9,XMONE,SKODE,1)
         IF (IPRINT .GT. 9) THEN
            CALL HEADER('Nuclear Contribution to Quadrupole Moment',-1)
            WRITE(LUPRI,1001) 'X','Y','Z'
            CALL OUTPUT(SKODN,1,3,1,3,3,3,1,LUPRI)
            CALL HEADER('Electronic Contribution to Quadrupole Moment',
     &                  -1)
            WRITE(LUPRI,1001) 'X','Y','Z'
            CALL OUTPUT(SKODE,1,3,1,3,3,3,1,LUPRI)
         ENDIF

         CALL DAXPY(9,ONE,SKODE,1,SKODN,1)

         CALL HEADER('Total Molecular Quadrupole Moment',-1)
         WRITE(LUPRI,1001) 'X','Y','Z'
         CALL OUTPUT(SKODN,1,3,1,3,3,3,1,LUPRI)

         CALL FLSHFO(LUPRI)

      ENDIF

C     Electronic 2nd moment of charge.
C     --------------------------------

      IF (SECMOM) THEN

         ICOUNT = ICOUNT + 1

         CALL AROUND('Electronic Second Moment of Charge')

C        Calculate electronic contribution.
C        ----------------------------------

         CALL DZERO(ELSEMO,9)

         IJ = 0
         DO I = 1,3
            DO J = I,3

               IJ = IJ + 1

               CALL DZERO(WORK(KONEP),N2BASX)
               FF   = ONE
               ISYM = -1
               CALL CC_ONEP(WORK(KONEP),WORK(KEND1),LWRK1,FF,ISYM,
     &                      LBLSM(IJ))

               IF (ISYM .EQ. 1) THEN
                  LENAB = N2BST(1)
                  CALL CCELQUA(WORK(KONEP),DENAO,LENAB,I,J,ELSEMO)
               ENDIF

            ENDDO
         ENDDO

C        Reorder.
C        --------

         CALL CC_QUAREO(ELSEMO,SKODE)

C        Print.
C        ------

         WRITE(LUPRI,1001) 'X','Y','Z'
         CALL OUTPUT(SKODE,1,3,1,3,3,3,1,LUPRI)

         CALL FLSHFO(LUPRI)

      ENDIF

C     Electric field gradient.
C     ========================

      IF (NQCC) THEN

         ICOUNT = ICOUNT + 1

         CALL AROUND('Electric Field Gradients')

C        Calculate nuclear contribution.
C        -------------------------------

         IOPT   = 2
         IPRLOC = IPRINT - 5
         CALL CCNUCQUA(WORK,LWORK,IOPT,IPRLOC)

C        Calculate electronic contribution.
C        ----------------------------------

         LENAB = N2BST(1)
         CALL CCELEFG(DENAO,LENAB,WORK,LWORK,IPRLOC)

C        Print result.
C        -------------

         KDIAG = KEND1
         KAXIS = KDIAG + 3*MXCENT
         KEND2 = KAXIS + 9*MXCENT
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND2-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         IPRLOC = 2
         ICCPRI = 2
         CALL NQCRES(IPRLOC,WORK(KDIAG),WORK(KAXIS),ICCPRI)

         CALL FLSHFO(LUPRI)

      ENDIF

C     Relativistic corrections to ground state energy.
C     ================================================

      IF (RELCOR) THEN

         ICOUNT = ICOUNT + 1

         CALL AROUND('Pauli Relativistic Corrections to the'
     &               //' Ground State Energy')

         DO IRC = 1,2

            IF (IRC .EQ. 1) LABEL = 'DARWIN  '
            IF (IRC .EQ. 2) LABEL = 'MASSVELO'

            CALL DZERO(WORK(KONEP),N2BASX)
            FF   = ONE
            ISYM = 1
            CALL CC_ONEP(WORK(KONEP),WORK(KEND1),LWRK1,FF,ISYM,LABEL)

            IF (IRC .EQ. 1) THEN
               DARW = DDOT(N2BST(1),WORK(KONEP),1,DENAO,1)
            ELSE IF (IRC .EQ. 2) THEN
               VELO = DDOT(N2BST(1),WORK(KONEP),1,DENAO,1)
            ENDIF

         ENDDO

         VEDAT = VELO + DARW
         WRITE(LUPRI,131) '1-elec. Darwin term:',DARW
         WRITE(LUPRI,132) '------------------- '
         WRITE(LUPRI,131) 'Mass-Velocity term: ',VELO
         WRITE(LUPRI,132) '------------------  '
         WRITE(LUPRI,133) 'Mass-Velocity + 1-elec. Darwin terms:',VEDAT
         WRITE(LUPRI,134) '------------------------------------ '

  131    FORMAT(/,9X,A20,F17.9)
  132    FORMAT(9X,A20)
  133    FORMAT(/,9X,A37,1X,F17.9)
  134    FORMAT(9X,A37)

      ENDIF

C     Relativistic two-electron Darwin term correction.
C     ***** NOT IMPLEMENTED YET FOR CHOLESKY CC2 ******
C     =================================================

      IF (DAR2EL) THEN
         ICOUNT = ICOUNT + 1
         WRITE(LUPRI,'(//,5X,A,A,/,5X,A,/,5X,A,A,A,//)')
     &   '*** NOTICE: Relativistic two-electron Darwin correction has ',
     &               'not been',
     &   '            implemented with Cholesky decomposed integrals.',
     & '            ',SECNAM,' cowardly refuses to process the request.'
      ENDIF

C     Electronic expectation value of arbitrary operator(s).
C     ======================================================

      IF (NAFOP .GT. 0) THEN
         ICOUNT = ICOUNT + 1
         CALL CC_FOPVAL(DENAO,WORK,LWORK)
      ENDIF

C     Electronic expectation value of commutators.
C     ============================================

      IF (NABFOP .GT. 0) THEN
         ICOUNT = ICOUNT + 1
         CALL CC_FOPCOM(DENAO,WORK,LWORK,MODFIL)
      ENDIF

C     Print end message.
C     ------------------

      IF (ICOUNT .GT. 0) THEN
         TIMTOT = SECOND() - TIMTOT
         WRITE(LUPRI,'(/,1X,A,F10.2,A,/)')
     &   'Requested first order one-electron properties calculated in ',
     &   TIMTOT,' seconds.'
      ENDIF

      RETURN
 1001 FORMAT(20X,A1,14X,A1,14X,A1)
      END
C  /* Deck cc_choden */
      SUBROUTINE CC_CHODEN(TBAR,DENAO,WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate and analyze (natural occupations) the CC2
C              ground state density matrix.
C              The density matrix is returned backtransformed to AO
C              basis in DENAO, including frozen core contributions.
C              The AO density is stored on disk (list 'd00')
C
#include <implicit.h>
      DIMENSION TBAR(*), DENAO(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <chocc2.h>
#include <ccfop.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_CHODEN')

      CHARACTER*3 LIST
      PARAMETER (LIST = 'd00')

      INTEGER ADR

C     Calculate MO density matrix in T1-transformed basis.
C     ----------------------------------------------------

      KDNIJ = 1
      KDNIA = KDNIJ + NMATIJ(1)
      KDNAB = KDNIA + NT1AM(1)
      KEND1 = KDNAB + NMATAB(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: den.'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL CC_CHODEN1(TBAR,WORK(KDNIJ),WORK(KDNIA),WORK(KDNAB),
     &                WORK(KEND1),LWRK1)

C     Write MO blocks to disk.
C     ------------------------

      ADR   = 1
      LDEN  = NMATIJ(1) + NT1AM(1) + NMATAB(1)
      LUDEN = 97
      CALL WOPEN2(LUDEN,FDEN00,64,0)
      CALL PUTWA2(LUDEN,FILDEN,WORK(KDNIJ),ADR,LDEN)
      CALL WCLOSE2(LUDEN,FILDEN,'KEEP')

C     Analysis: CC2 natural occupations.
C     ----------------------------------

      CALL CC_CHODENO(WORK(KDNIJ),TBAR,WORK(KDNIA),WORK(KDNAB),
     &                WORK(KEND1),LWRK1)

C     Back-transform from T1-transformed to AO basis.
C     -----------------------------------------------

      KLAMDP = KEND1
      KLAMDH = KLAMDP + NGLMDT(1)
      KT1AM  = KLAMDH + NGLMDT(1)
      KEND2  = KT1AM  + NT1AM(1)
      LWRK2  = LWORK  - KEND2 + 1

      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: Lambda'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND2-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      IFAIL = -1
      CALL CHO_RDRST(DUM1,WORK(KT1AM),DUM2,.FALSE.,.TRUE.,.FALSE.,IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND2),
     &            LWRK2)

      KEND2 = KT1AM
      LWRK2 = LWORK - KEND2 + 1
      CALL DZERO(DENAO,N2BST(1))
      CALL CC_DENAO(DENAO,1,TBAR,WORK(KDNAB),WORK(KDNIJ),WORK(KDNIA),1,
     &              WORK(KLAMDP),1,WORK(KLAMDH),1,WORK(KEND2),LWRK2)

C     Add frozen core contributions (if any).
C     ---------------------------------------

      IF (FROIMP .OR. FROEXP) THEN
         CALL CC_FCD1AO(DENAO,WORK(KEND2),LWRK2,'CC2       ')
      ENDIF

C     Write AO density to disk.
C     -------------------------

      CALL CC_WRRSPD(LIST,1,1,'CC2       ',RELORB,DENAO,
     &               WORK(KEND2),LWRK2)

      RETURN
      END
C  /* Deck cc_chodeno */
      SUBROUTINE CC_CHODENO(DENIJ,DENAI,DENIA,DENAB,WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate natural occupations. Skipped for calculations
C              with frozen orbitals.
C
#include <implicit.h>
      DIMENSION DENIJ(*), DENAI(*), DENIA(*), DENAB(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHODENO')

      PARAMETER (ZERO = 0.0D0)

C     Print the non-symmetric MO (p,h) density matrix.
C     ------------------------------------------------

      IF (IPRINT .GT. 50) THEN
         TRCO = ZERO
         TRCV = ZERO
         DO ISYM = 1,NSYM
            DO I = 1,NRHF(ISYM)
               KII  = IMATIJ(ISYM,ISYM) + NRHF(ISYM)*(I - 1) + I
               TRCO = TRCO + DENIJ(KII)
            ENDDO
            DO A = 1,NVIR(ISYM)
               KAA  = IMATAB(ISYM,ISYM) + NVIR(ISYM)*(A - 1) + A
               TRCV = TRCV + DENAB(KAA)
            ENDDO
         ENDDO
         CALL AROUND('CC2 Non-Symmetric One-Electron Density Matrix')
         WRITE(LUPRI,'(10X,A)')
     &   '- Density matrix in T1-transformed (Lamda-p/Lambda-h) basis'
         CALL HEADER('Occupied Part',-1)
         CALL NOCC_PRT(DENIJ,1,'IJ  ')
         CALL HEADER('Single Excitation Part',-1)
         CALL NOCC_PRT(DENAI,1,'AI  ')
         CALL HEADER('Single De-Excitation Part',-1)
         CALL NOCC_PRT(DENIA,1,'AI  ')
         CALL HEADER('Virtual Part',-1)
         CALL NOCC_PRT(DENAB,1,'AB  ')
         WRITE(LUPRI,'(/,5X,A,1P,D22.15)')
     &   'Trace, occupied part: ',TRCO
         WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &   'Trace, virtual  part: ',TRCV
         WRITE(LUPRI,'(5X,A,1P,D22.15,/)')
     &   'Trace, total        : ',TRCO+TRCV
      ENDIF

      IF (.NOT. FROIMP) THEN

C        Add extra T1-dependent terms for natural occupations.
C        -----------------------------------------------------

         KDNIJ = 1
         KDNAI = KDNIJ + NMATIJ(1)
         KDNIA = KDNAI + NT1AM(1)
         KDNAB = KDNIA + NT1AM(1)
         KYMAT = KDNAB + NMATAB(1)
         KT1AM = KYMAT + NMATAB(1)
         KDEN  = KT1AM + NT1AM(1)
         KEND1 = KDEN  + N2BST(1)
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      'Insufficient memory in ',SECNAM
            WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &      'Need     : ',KEND1-1,
     &      'Available: ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         CALL DCOPY(NMATIJ(1),DENIJ,1,WORK(KDNIJ),1)
         CALL DCOPY(NT1AM(1),DENAI,1,WORK(KDNAI),1)
         CALL DCOPY(NT1AM(1),DENIA,1,WORK(KDNIA),1)
         CALL DCOPY(NMATAB(1),DENAB,1,WORK(KDNAB),1)
         CALL DCOPY(NMATAB(1),DENAB,1,WORK(KYMAT),1)
         CALL CC_EITRV(WORK(KYMAT),WORK(KEND1),LWRK1,1)
         IFAIL = -1
         CALL CHO_RDRST(DUM1,WORK(KT1AM),DUM2,.FALSE.,.TRUE.,.FALSE.,
     &                  IFAIL)

         CALL CCD1T1CO(WORK(KDNAI),WORK(KDNAB),WORK(KDNIJ),WORK(KDNIA),
     &                 WORK(KT1AM),WORK(KYMAT),WORK(KEND1),LWRK1)

C        Print the non-symmetric MO (CMO) density matrix.
C        ------------------------------------------------

         IF (IPRINT .GT. 50) THEN
            TRCO = ZERO
            TRCV = ZERO
            DO ISYM = 1,NSYM
               DO I = 1,NRHF(ISYM)
                  KII  = KDNIJ + IMATIJ(ISYM,ISYM) + NRHF(ISYM)*(I - 1)
     &                 + I - 1
                  TRCO = TRCO + WORK(KII)
               ENDDO
               DO A = 1,NVIR(ISYM)
                  KAA  = KDNAB + IMATAB(ISYM,ISYM) + NVIR(ISYM)*(A - 1)
     &                 + A - 1
                  TRCV = TRCV + WORK(KAA)
               ENDDO
            ENDDO
            CALL AROUND('CC2 Non-Symmetric One-Electron Density Matrix')
            WRITE(LUPRI,'(10X,A)')
     &      '- Density matrix in molecular orbital (CMO) basis'
            CALL HEADER('Occupied Part',-1)
            CALL NOCC_PRT(WORK(KDNIJ),1,'IJ  ')
            CALL HEADER('Single Excitation Part',-1)
            CALL NOCC_PRT(WORK(KDNAI),1,'AI  ')
            CALL HEADER('Single De-Excitation Part',-1)
            CALL NOCC_PRT(WORK(KDNIA),1,'AI  ')
            CALL HEADER('Virtual Part',-1)
            CALL NOCC_PRT(WORK(KDNAB),1,'AB  ')
            WRITE(LUPRI,'(/,5X,A,1P,D22.15)')
     &      'Trace, occupied part: ',TRCO
            WRITE(LUPRI,'(5X,A,1P,D22.15)')
     &      'Trace, virtual  part: ',TRCV
            WRITE(LUPRI,'(5X,A,1P,D22.15,/)')
     &      'Trace, total        : ',TRCO+TRCV
         ENDIF

C        Reorder density.
C        ----------------

         CALL DZERO(WORK(KDEN),N2BST(1))
         CALL CCD1REOR(WORK(KDEN),WORK(KDNAI),WORK(KDNAB),WORK(KDNIJ),
     &                 WORK(KDNIA))

C        Natural occupations.
C        --------------------

         CALL CC_DEDIAN(WORK(KDEN),'CC2 ',WORK(KEND1),LWRK1)

      ENDIF

      RETURN
      END
C  /* Deck cc_choden1 */
      SUBROUTINE CC_CHODEN1(TBAR,DENIJ,DENIA,DENAB,WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the CC2 one-electron density matrix (block-wise).
C
C     Formulas:
C
C        DENIJ(ij) = 2 delta(ij) - sum(abk) t(ai,bk) * tb(aj,bk)
C
C        DENIA(ai) = sum(bj) [2*t(ai,bj)-t(aj,bi)] * TBAR(bj)
C
C        DENAB(ab) = sum(cij) tb(ai,cj) * t(bi,cj)
C
C     The excitation block is simply equal to TBAR(ai).
C
C     The calculation is based on a decomposition of the t2-amplitudes.
C
#include <implicit.h>
      DIMENSION TBAR(*), DENIJ(*), DENIA(*), DENAB(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <chocc2.h>

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CC_CHODEN1')

      PARAMETER (ZERO = 0.0D0, TWO = 2.0D0)

C     Read canonical orbital energies.
C     --------------------------------

      KFOCKD = 1
      KEND1  = KFOCKD + NORBTS
      LWRK1  = LWORK  - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: orb. en.'
         WRITE(LUPRI,'(//,5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL CHO_RDSIR(DUM1,DUM2,WORK(KFOCKD),DUM3,WORK(KEND1),LWRK1,
     &               .FALSE.,.TRUE.,.FALSE.)

C     Calculate DENIA. The amplitudes will be decomposed in CC_CIM1
C     (if it hasn't been done yet).
C     -------------------------------------------------------------

      CALL DZERO(DENIA,NT1AM(1))
      CALL CC_CHOCIM1(WORK(KFOCKD),TBAR,DENIA,WORK(KEND1),LWRK1,1,1)

C     Allocation.
C     -----------

      KT1AM  = KEND1
      KLAMDP = KT1AM  + NT1AM(1)
      KLAMDH = KLAMDP + NGLMDT(1)
      KEND2  = KLAMDH + NGLMDT(1)
      LWRK2  = LWORK  - KEND2 + 1

      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(//,5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND2-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Read singles amplitudes and calculate Lambda matrices.
C     ------------------------------------------------------

      IFAIL = -1
      CALL CHO_RDRST(DUM1,WORK(KT1AM),DUM2,.FALSE.,.TRUE.,.FALSE.,IFAIL)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND2),
     &            LWRK2)

C     Calculate perturbed Cholesky vectors.
C     -------------------------------------

      ISIDE = -1
      CALL CC_CHOTG(WORK(KLAMDP),WORK(KLAMDH),TBAR,
     &              WORK(KEND2),LWRK2,1,1,ISIDE)

C     Read the F(ia) matrix in place of T1AM. Regain work space
C     occupied by Lambda matrices (no longer needed).
C     ----------------------------------------------------------

      KFIA  = KT1AM
      KEND2 = KFIA  + NT1AM(1)
      LWRK2 = LWORK - KEND2 + 1

      CALL ONEL_OP(-1,3,LUFIA)
      CALL CHO_MOREAD(WORK(KFIA),NT1AM(1),1,1,LUFIA)
      CALL ONEL_OP(1,3,LUFIA)

C     Calculate Y intermediates,
C        Y(ai,J) = sum(bj) tb(ai,bj) * M(bj,J)
C     where M(bj,J) are the Cholesky vectors representing the
C     (minus) t2 amplitudes. Delete the perturbed Cholesky
C     vector files after use.
C     -------------------------------------------------------

c      lwsav = lwrk2
c      nwrk  = nt2sq(1)/5 + 2*nt1am(1)
c      lwrk2 = min(lwrk2,nwrk)
c      write(lupri,*) '   ...calling cc_cyi with lwrk2 = ',lwrk2

c     ISIDE = -3
c     FREQ  = ZERO
c     CALL CC_CYI(WORK(KFOCKD),DUM1,TBAR,WORK(KFIA),WORK(KEND2),LWRK2,
c    &            ISIDE,1,1,FREQ,XINRM,TBNRM,.TRUE.)

      FREQ = ZERO
      CALL CC_CYID0(WORK(KFOCKD),TBAR,WORK(KFIA),WORK(KEND2),LWRK2,
     &              FREQ,XINRM,TBNRM,.FALSE.,.TRUE.)

c      lwrk2 = lwsav

C     Calculate occupied and virtual blocks of the density.
C     Delete the Y intermediate files after use.
C     -----------------------------------------------------

c      lwsav = lwork
c      nwrk  = 5*nt1am(1)
c      lwork = min(lwork,nwrk)
c      write(lupri,*) '   ...calling cc_choey with lwork = ',lwork

      CALL DZERO(DENIJ,NMATIJ(1))
      DO ISYMI = 1,NSYM
         DO I = 1,NRHF(ISYMI)
            II = IMATIJ(ISYMI,ISYMI) + NRHF(ISYMI)*(I - 1) + I
            DENIJ(II) = TWO
         ENDDO
      ENDDO
      CALL DZERO(DENAB,NMATAB(1))

      ISIDE = -3
      ISYMY = 1
      ITYP  = 6
      ISYCH = 1
      CALL CC_CHOEY(DENIJ,DENAB,WORK,LWORK,ISYMY,ISIDE,ISYCH,ITYP,
     &              NCHMOC,.TRUE.)

c      lwork = lwsav

      RETURN
      END
C  /* Deck cc_cyid0 */
      SUBROUTINE CC_CYID0(FOCKD,TBAR,FIA,WORK,LWORK,FREQ,XINRM,TBNRM,
     &                    DELP1,DELP2)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate Y intermediates for ground state density matrix.
C
#include <implicit.h>
      DIMENSION FOCKD(*), TBAR(*), FIA(*), WORK(LWORK)
      LOGICAL DELP1, DELP2
#include <maxorb.h>
#include <ccdeco.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <dummy.h>

C     Make sure that ground state amplitudes have been decomposed.
C     ------------------------------------------------------------

      IF (.NOT. CHOT2C) THEN
         CALL CC_CHOCIM1(FOCKD,DUMMY,DUMMY,WORK,LWORK,0,0)
      ENDIF

C     Set up amplitude assembly.
C     --------------------------

      FAC1   = 1.0D0
      FAC2   = 1.0D0
      SCD    = 1.0D0
      SCDG   = 1.0D0

      NUMP12 = 1
      JTYP1  = 1
      JTYP2  = 8
      ISYMP1 = 1
      ISYMP2 = 1

      IOPTDN = 1
      IOPTCE = 1

C     Set info for Y intermediates.
C     -----------------------------

      NUMZ  = 1
      JTYPZ = 6
      ISYMZ = 1
      JTYPY = -1

C     Add L(ia,J) to LT(ia,J).
C     ------------------------

      ALPHA = 1.0D0
      BETA  = 1.0D0
      CALL CC_CHOADD(JTYP1,JTYP2,1,NUMCHO,WORK,LWORK,ALPHA,BETA)

C     Calculate Y intermediates.
C     --------------------------

      CALL CC_CYI(JTYP1,ISYMP1,JTYP2,ISYMP2,NUMCHO,FAC1,NUMP12,
     &            TBAR,1,FIA,1,FAC2,1,
     &            SCD,SCDG,IOPTDN,FOCKD,FREQ,IOPTCE,
     &            JTYPZ,ISYMZ,NCHMOC,NUMZ,JTYPY,
     &            DUMMY,IDUMMY,IDUMMY,DUMMY,
     &            WORK,LWORK,XINRM,TBNRM,DELP1,DELP2)

      RETURN
      END
C  /* Deck cc_eitrv */
      SUBROUTINE CC_EITRV(E1INT,WORK,LWORK,ISYINT)
C
C     Thomas Bondo Pedersen, February 2003.
C     - copy of the virtual part of CC_EITR by A. Halkier.
C
C     Purpose: Transpose indices of virtual array E1INT,
C
C              E1INT(b,a) <--- E1INT(a,b)
C
#include <implicit.h>
      DIMENSION E1INT(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_EITRV')

      IF (LWORK .LT. NMATAB(ISYINT)) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need     : ',NMATAB(ISYINT),
     &   'Available: ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

      CALL DCOPY(NMATAB(ISYINT),E1INT,1,WORK,1)

      DO ISYMA = 1,NSYM
         ISYMB = MULD2H(ISYMA,ISYINT)
         DO B = 1,NVIR(ISYMB)
            NORIG = IMATAB(ISYMA,ISYMB) + NVIR(ISYMA)*(B - 1) + 1
            NTRAN = IMATAB(ISYMB,ISYMA) + B
            CALL DCOPY(NVIR(ISYMA),WORK(NORIG),1,E1INT(NTRAN),
     &                 NVIR(ISYMB))
         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck cc_chorsdtst */
      SUBROUTINE CC_CHORSDTST(TBAR,WORK,LWORK)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Calculate the residual vector for ground state multiplier
C              effective equation (for debugging purposes).
C
#include <implicit.h>
      DIMENSION TBAR(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHORSDTST')

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0)

C     Allocation.
C     -----------

      KETA  = 1
      KTRF  = KETA  + NT1AM(1)
      KEND1 = KTRF  + NT1AM(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: eta'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Calculate the effective right-hand side.
C     ----------------------------------------

      CALL CC_CHOETA(WORK(KETA),WORK(KEND1),LWRK1)
      CALL DSCAL(NT1AM(1),XMONE,WORK(KETA),1)

C     Calculate transformed vector.
C     -----------------------------

      FREQ  = ZERO
      ISYMX = 1
      NUMX  = 1
      ISIDE = -1

      CALL CC_CHOATR(WORK(KTRF),TBAR,FREQ,WORK(KEND1),LWRK1,ISYMX,NUMX,
     &               ISIDE)

C     Analysis.
C     ---------

c      CALL HEADER('Minus ETA',-1)
c      CALL NOCC_PRT(WORK(KETA),1,'AI  ')

c      CALL HEADER('Transformed Vector',-1)
c      CALL NOCC_PRT(WORK(KTRF),1,'AI  ')

      CALL DBGDIFAI(WORK(KETA),WORK(KTRF),ETANRM,TRFNRM,ERRMIN,ERRMAX,1)
      CALL DAXPY(NT1AM(1),XMONE,WORK(KTRF),1,WORK(KETA),1)
      RSDNRM = DSQRT(DDOT(NT1AM(1),WORK(KETA),1,WORK(KETA),1))

      TBRNRM = DSQRT(DDOT(NT1AM(1),TBAR,1,TBAR,1))
      CALL HEADER('Testing LH Transformation',-1)
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Norm of singles multipliers: ',TBRNRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Norm of transformed vector : ',TRFNRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Norm of effective RH side  : ',ETANRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Residual norm              : ',RSDNRM
      WRITE(LUPRI,'(5X,A,1P,D22.15)')
     & 'Minimum absolute error     : ',ERRMIN
      WRITE(LUPRI,'(5X,A,1P,D22.15,//)')
     & 'Maximum absolute error     : ',ERRMAX
      CALL HEADER('Residual Vector',-1)
      CALL NOCC_PRT(WORK(KETA),1,'AI  ')

      RETURN
      END
C  /* Deck cc_chofopdbg */
      SUBROUTINE CC_CHOFOPDBG(WORK,LWORK,SAFTST)
C
C     Thomas Bondo Pedersen, February 2003.
C
C     Purpose: Read in singles multipliers from conventional calculation
C              and test the density matrix calculation by calculation
C              first-order properties. Test also the left- and right-hand
C              side transformations by transforming with unit vectors.
C
C     IF (SAFTST): kill job after tests....
C
C     Warning: requires that a file 'AM_FOR_CHOTST' has been written by the
C              conventional code.....
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
      LOGICAL SAFTST
#include <maxorb.h>
#include <ccdeco.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <chocc2.h>
#include <priunit.h>

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOFOPDBG')

C     Make sure that amplitude deco. is part of test.
C     -----------------------------------------------

      CHOT2  = .FALSE.
      CHOT2C = .FALSE.

C     Allocation.
C     -----------

      KSAVE = 1
      KSAVT = KSAVE + 1
      KSAVO = KSAVT + NT1AM(1)
      KDEN  = KSAVO + NT1AM(1)
      KTBAR = KDEN  + MAX(N2BST(1),NT1AM(1))
      KLAMP = KTBAR + NT1AM(1)
      KLAMH = KLAMP + NGLMDT(1)
      KTRAC = KLAMH + NGLMDT(1)
      KEND1 = KTRAC + NUMCHO(1)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   'Insufficient memory in ',SECNAM
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND1-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Get data.
C     ---------

      IFAIL = -1
      CALL CHO_RDRST(WORK(KSAVE),WORK(KSAVT),WORK(KSAVO),.TRUE.,.TRUE.,
     &               .TRUE.,IFAIL)

      LUTBAR = -1
      CALL GPOPEN(LUTBAR,'AM_FOR_CHOTST','OLD',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND(LUTBAR)
      READ(LUTBAR) (WORK(KTBAR+I-1), I = 1,NT1AM(1))
      CALL GPCLOSE(LUTBAR,'KEEP')

c      DELTAE = 0.0D0
c      CALL CHO_WRRST(WORK(KSAVE),WORK(KDEN),WORK(KSAVO))
c      CALL LAMMAT(WORK(KLAMP),WORK(KLAMH),WORK(KDEN),WORK(KEND1),LWRK1)
c      CALL CC2_TRFL(WORK(KDEN),WORK(KLAMP),WORK(KLAMH),WORK(KTRAC),
c     &              WORK(KEND1),LWRK1,DELTAE)

c      T1NRM = DDOT(NT1AM(1),WORK(KDEN),1,WORK(KDEN),1)
      T1NRM = DDOT(NT1AM(1),WORK(KSAVT),1,WORK(KSAVT),1)
      TBNRM = DDOT(NT1AM(1),WORK(KTBAR),1,WORK(KTBAR),1)
      WRITE(LUPRI,'(//,5X,A,A)')
     & SECNAM,': Using singles multipliers from conv. calc. for FOPs'
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15)')
     & 'Square norm of singles amplitudes : ',T1NRM,
     & '2-norm      of singles amplitudes : ',DSQRT(T1NRM)
      WRITE(LUPRI,'(5X,A,1P,D22.15,/,5X,A,1P,D22.15,//)')
     & 'Square norm of singles multipliers: ',TBNRM,
     & '2-norm      of singles multipliers: ',DSQRT(TBNRM)

c      CALL DBGDIFAI(WORK(KDEN),WORK(KSAVT),TRGNRM,TSTNRM,ERRMIN,ERRMAX,
c     &              1)
c      WRITE(LUPRI,'(5X,A,1P,D22.15)')
c     & 'Norm of conventional T1AM: ',TRGNRM
c      WRITE(LUPRI,'(5X,A,1P,D22.15)')
c     & 'Norm of Cholesky     T1AM: ',TSTNRM
c      WRITE(LUPRI,'(5X,A,1P,D22.15)')
c     & 'Difference               : ',TRGNRM-TSTNRM
c      WRITE(LUPRI,'(5X,A,1P,D22.15)')
c     & 'Min. absolute error      : ',ERRMIN
c      WRITE(LUPRI,'(5X,A,1P,D22.15,//)')
c     & 'Max. absolute error      : ',ERRMAX

      CALL DZERO(WORK(KDEN),NT1AM(1))
      CALL CC_CHOATRLHDBG(WORK(KDEN),WORK(KTBAR),WORK(KEND1),LWRK1,1)

C     Check LH and RH transformations (trf. by unit vectors).
C     -------------------------------------------------------

      CALL CC_CHOAEFFD(WORK(KEND1),LWRK1,.TRUE.,WORK(KTBAR))

C     Check residual with conv. vector.
C     ---------------------------------

      CALL CC_CHORSDTST(WORK(KTBAR),WORK(KEND1),LWRK1)

C     Kill if requested.
C     ------------------

      IF (SAFTST) THEN
         WRITE(LUPRI,'(//,5X,A,A)')
     &   SECNAM,': aborts execution as requested...'
         CALL QUIT('Killed by request in '//SECNAM)
      ENDIF

      RETURN
      END
C  /* Deck cc_comrl */
      SUBROUTINE CC_COMRL(COM,WORK,LWORK,MODEL)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate the expectation values of the commutators
C              [r(j),L(j)] where L(j) = [r x grad](j) is the imaginary
C              part of the angular momentum operator.
C
C     The ground state AO density matrix must be available on disk.
C
#include <implicit.h>
      DIMENSION COM(3), WORK(LWORK)
      CHARACTER*10 MODEL
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>
#include <dummy.h>

      CHARACTER*8 SECNAM, LABEL1(3), LABEL2(3)
      PARAMETER (SECNAM = 'CC_COMRL')

      DATA LABEL1 /'XDIPLEN ','YDIPLEN ','ZDIPLEN '/
      DATA LABEL2 /'XANGMOM ','YANGMOM ','ZANGMOM '/

C     Allocation.
C     -----------

      KDNPH  = 1
      KLAMDP = KDNPH  + MAX(N2BST(1),NT1AM(1))
      KLAMDH = KLAMDP + NGLMDT(1)
      KEND1  = KLAMDH + NGLMDT(1)
      LWRK1  = LWORK  - KEND1 + 1

      KT1AM  = KDNPH

      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//' [1]')
      ENDIF

C     Read 0th order singles amplitudes.
C     ----------------------------------

      IOPT = 1
      CALL CC_RDRSP('R0 ',1,1,IOPT,MODEL,WORK(KT1AM),DUMMY)

C     Calculate Lambda matrices.
C     --------------------------

      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)

C     Calculate particle-hole density.
C     IFR = 1: Include frozen orbital contributions.
C     ----------------------------------------------

      IFR = 1
      CALL CC_PHDEN(WORK(KLAMDP),WORK(KLAMDH),WORK(KDNPH),
     &              WORK(KEND1),LWRK1,1,1,IFR)

C     Re-allocate.
C     ------------

      KDNAO = KDNPH + N2BST(1)
      KINT1 = KDNAO + N2BST(1)
      KINT2 = KINT1 + N2BASX
      KEND2 = KINT2 + N2BASX
      LWRK2 = LWORK - KEND2 + 1

      IF (LWRK2 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//' [2]')
      ENDIF

C     Read AO density.
C     NOTE: Will fail if orbital relaxation included in density.
C     ----------------------------------------------------------

      IOPT = 1
      CALL CC_RDRSPD('d00',1,1,IOPT,MODEL,.FALSE.,WORK(KDNAO))

C     Loop over commutators.
C     ----------------------

      DO I = 1,3

C        Read AO integrals.
C        ------------------

         IERR1 = -1
         CALL CCPRPAO(LABEL1(I),WORK(KINT1),ISYM1,IHER1,IERR1,
     &                WORK(KEND2),LWRK2)
         IF (IERR1 .GT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      SECNAM,': I/O error ocurred in CCPRPAO reading ',LABEL1(I)
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR1 = ',IERR1
            CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
         ELSE IF (IERR1 .LT. 0) THEN
            WRITE(LUPRI,'(/,5X,A,A,A)')
     &      SECNAM,': irrep not determined for operator ',LABEL1(I)
            WRITE(LUPRI,'(5X,A,I10)')
     &      'CCPRPAO returned IERR1 = ',IERR1
            CALL QUIT('Irrep not determined in CCPRPAO ('//SECNAM//')')
         ENDIF

         IERR2 = -1
         CALL CCPRPAO(LABEL2(I),WORK(KINT2),ISYM2,IHER2,IERR2,
     &                WORK(KEND2),LWRK2)
         IF (IERR2 .GT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      SECNAM,': I/O error ocurred in CCPRPAO reading ',LABEL2(I)
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR2 = ',IERR2
            CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
         ELSE IF (IERR2 .LT. 0) THEN
            WRITE(LUPRI,'(/,5X,A,A,A)')
     &      SECNAM,': irrep not determined for operator ',LABEL2(I)
            WRITE(LUPRI,'(5X,A,I10)')
     &      'CCPRPAO returned IERR2 = ',IERR2
            CALL QUIT('Irrep not determined in CCPRPAO ('//SECNAM//')')
         ENDIF

C        Calculate commutator expectation value.
C        IFAIL = -1: abort if Hermiticity error.
C        ---------------------------------------

         IFAIL = -1
         CALL CC_COMEXP(WORK(KINT1),LABEL1(I),ISYM1,IHER1,
     &                  WORK(KINT2),LABEL2(I),ISYM2,IHER2,
     &                  WORK(KDNPH),1,WORK(KDNAO),1,
     &                  COM(I),WORK(KEND2),LWRK2,IFAIL)

      ENDDO

      RETURN
      END
C  /* Deck cc_comexp */
      SUBROUTINE CC_COMEXP(XINTA,LABELA,ISYMA,IHERA,
     &                     XINTB,LABELB,ISYMB,IHERB,
     &                     DENPH,ISYMPH,DENAO,ISYDEN,
     &                     COM,WORK,LWORK,IFAIL)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate the expectation value of the commutator
C              [A,B].
C
C     If on input IFAIL = -1, then a failure will be fatal.
C     Else the program will continue with COM=0.0D0.
C
C     On exit, IFAIL = 0 means "all OK" and IFAIL = 1 means "failure".
C
#include <implicit.h>
      DIMENSION XINTA(*), XINTB(*), DENPH(*), DENAO(*), WORK(LWORK)
      CHARACTER*8 LABELA, LABELB
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_COMEXP')

      LOGICAL ZBYSYM, ZBYHER, LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0)

C     Debug print.
C     ------------

      IF (LOCDBG) THEN
         XANRM = DSQRT(DDOT(N2BST(ISYMA),XINTA,1,XINTA,1))
         XBNRM = DSQRT(DDOT(N2BST(ISYMB),XINTB,1,XINTB,1))
         XPHNM = DSQRT(DDOT(N2BST(ISYMPH),DENPH,1,DENPH,1))
         XAONM = DSQRT(DDOT(N2BST(ISYDEN),DENAO,1,DENAO,1))
         WRITE(LUPRI,*) SECNAM,':'
         WRITE(LUPRI,*) 'Calculating <[',LABELA,',',LABELB,']>'
         WRITE(LUPRI,*) 'ISYMA , IHERA  : ',ISYMA,IHERA
         WRITE(LUPRI,*) 'ISYMB , IHERB  : ',ISYMB,IHERB
         WRITE(LUPRI,*) 'ISYMPH, ISYDEN : ',ISYMPH,ISYDEN
         WRITE(LUPRI,*) '2-norm of XINTA: ',XANRM
         WRITE(LUPRI,*) '2-norm of XINTB: ',XBNRM
         WRITE(LUPRI,*) '2-norm of DENPH: ',XPHNM
         WRITE(LUPRI,*) '2-norm of DENAO: ',XAONM
         CALL FLSHFO(LUPRI)
      ENDIF

C     Save and set IFAIL for output.
C     ------------------------------

      IFAILS = IFAIL
      IFAIL  = 0

C     Set commutator symmetry and flags for "zero by symmetry"
C     and "zero by Hermoticity".
C     --------------------------------------------------------

      ISYM   = MULD2H(ISYMA,ISYMB)
      ISYCOM = MULD2H(ISYM,ISYMPH)
      ZBYSYM = ISYCOM .NE. ISYDEN

      IHCOM  = -IHERA*IHERB
      ZBYHER = IHCOM .EQ. -1

C     Do calculation according to symmetry and Hermiticity.
C     -----------------------------------------------------

      IF (ZBYSYM) THEN

         COM = ZERO

      ELSE

         IF (IHCOM .EQ. -1) THEN

C           Real anti-Hermitian.
C           --------------------

            COM = ZERO

         ELSE IF (IHCOM .EQ. 1) THEN

C           Real Hermitian.
C           ---------------

            KCMAO = 1
            KEND1 = KCMAO + N2BST(ISYCOM)
            LWRK1 = LWORK - KEND1 + 1

            IF (LWRK1 .LT. 0) THEN
               CALL QUIT('Insufficient memory in '//SECNAM//' [1]')
            ENDIF

            DO ISYMBE = 1,NSYM

               ISYMDE = MULD2H(ISYMBE,ISYMB)
               ISYMGA = MULD2H(ISYMDE,ISYMPH)
               ISYMAL = MULD2H(ISYMGA,ISYMA)

               NEED = NBAS(ISYMGA)*NBAS(ISYMBE)

               IF (LWRK1 .LT. NEED) THEN
                  CALL QUIT('Insufficient memory in '//SECNAM//' [2]')
               ENDIF

               NTOTA = MAX(NBAS(ISYMAL),1)
               NTOTG = MAX(NBAS(ISYMGA),1)
               NTOTD = MAX(NBAS(ISYMDE),1)

               KOFF1 = IAODIS(ISYMDE,ISYMGA) + 1
               KOFF2 = IAODIS(ISYMDE,ISYMBE) + 1
               CALL DGEMM('T','N',
     &                    NBAS(ISYMGA),NBAS(ISYMBE),NBAS(ISYMDE),
     &                    ONE,DENPH(KOFF1),NTOTD,XINTB(KOFF2),NTOTD,
     &                    ZERO,WORK(KEND1),NTOTG)

               KOFF1 = IAODIS(ISYMAL,ISYMGA) + 1
               KOFF2 = KCMAO + IAODIS(ISYMAL,ISYMBE)
               CALL DGEMM('N','N',
     &                    NBAS(ISYMAL),NBAS(ISYMBE),NBAS(ISYMGA),
     &                    ONE,XINTA(KOFF1),NTOTA,WORK(KEND1),NTOTG,
     &                    ZERO,WORK(KOFF2),NTOTA)

            ENDDO

            IF (LOCDBG) THEN
               XCMNM = DSQRT(DDOT(N2BST(ISYCOM),WORK(KCMAO),1,
     &                                          WORK(KCMAO),1))
               WRITE(LUPRI,*) '2-norm of CMAO : ',XCMNM,' [1]'
            ENDIF

            DO ISYMBE = 1,NSYM

               ISYMDE = MULD2H(ISYMBE,ISYMA)
               ISYMGA = MULD2H(ISYMDE,ISYMPH)
               ISYMAL = MULD2H(ISYMGA,ISYMB)

               NEED = NBAS(ISYMGA)*NBAS(ISYMBE)

               IF (LWRK1 .LT. NEED) THEN
                  CALL QUIT('Insufficient memory in '//SECNAM//' [3]')
               ENDIF

               NTOTA = MAX(NBAS(ISYMAL),1)
               NTOTG = MAX(NBAS(ISYMGA),1)
               NTOTD = MAX(NBAS(ISYMDE),1)
         
               KOFF1 = IAODIS(ISYMDE,ISYMGA) + 1
               KOFF2 = IAODIS(ISYMDE,ISYMBE) + 1
               CALL DGEMM('T','N',
     &                    NBAS(ISYMGA),NBAS(ISYMBE),NBAS(ISYMDE),
     &                    ONE,DENPH(KOFF1),NTOTD,XINTA(KOFF2),NTOTD,
     &                    ZERO,WORK(KEND1),NTOTG)

               KOFF1 = IAODIS(ISYMAL,ISYMGA) + 1
               KOFF2 = KCMAO + IAODIS(ISYMAL,ISYMBE)
               CALL DGEMM('N','N',
     &                    NBAS(ISYMAL),NBAS(ISYMBE),NBAS(ISYMGA),
     &                    XMONE,XINTB(KOFF1),NTOTA,WORK(KEND1),NTOTG,
     &                    ONE,WORK(KOFF2),NTOTA)

            ENDDO

            IF (LOCDBG) THEN
               XCMNM = DSQRT(DDOT(N2BST(ISYCOM),WORK(KCMAO),1,
     &                                          WORK(KCMAO),1))
               WRITE(LUPRI,*) '2-norm of CMAO : ',XCMNM,' [2]'
            ENDIF

            COM = DDOT(N2BST(ISYDEN),WORK(KCMAO),1,DENAO,1)

         ELSE

C           Hermiticity out of bounds.
C           --------------------------

            IF ((IPRINT.GT.15) .OR. (IFAILS.EQ.-1)) THEN
               WRITE(LUPRI,'(//,5X,A,A,A,A,A)')
     &         'Expectation value requested for commutator [',
     &         LABELA,',',LABELB,']:'
               WRITE(LUPRI,'(5X,A,A,A,I2,1X,I2)')
     &         'Symmetry and Hermiticity of ',LABELA,': ',ISYMA,IHERA
               WRITE(LUPRI,'(5X,A,A,A,I2,1X,I2)')
     &         'Symmetry and Hermiticity of ',LABELB,': ',ISYMB,IHERB
            ENDIF
            IF (IFAILS .EQ. -1) THEN
               WRITE(LUPRI,'(5X,A,I4)')
     &         'Unable to handle commutator Hermiticity: ',
     &          IHCOM
               CALL QUIT('Error in '//SECNAM)
            ELSE
               IF (IPRINT .GT. 15) THEN
                  WRITE(LUPRI,'(5X,A,I4)')
     &            'WARNING: Unable to handle commutator Hermiticity: ',
     &             IHCOM
                  WRITE(LUPRI,'(5X,A)')
     &   'Program continues nevetheless using zero expectation value...'
               ENDIF
               COM   = ZERO
               IFAIL = 1
            ENDIF

         ENDIF

      ENDIF

C     Print.
C     ------

      IF (IPRINT .GT. 50) THEN
         WRITE(LUPRI,'(/,5X,A,A,A,A,A,1P,D22.15,A)')
     &   '<[',LABELA,',',LABELB,']> = ',COM,' au'
         IF (IFAIL .EQ. 1) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '(Hermiticity failure - see above)'
         ELSE IF (ZBYHER) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '(Zero by Hermiticity)'
         ELSE IF (ZBYSYM) THEN
            WRITE(LUPRI,'(5X,A)')
     &      '(Zero by symmetry)'
         ENDIF
      ENDIF

      RETURN
      END
C  /* Deck cc_phden */
      SUBROUTINE CC_PHDEN(XLAMDP,XLAMDH,DENPH,WORK,LWORK,ISYMLP,ISYMLH,
     &                    IFR)
C
C     Thomas Bondo Pedersen, April 2003.
C     - based on CC_AODENS by H. Koch and A. Sanchez.
C
C     Purpose: Calculate the full particle-hole density matrix, summing
C              over ALL orbitals (not just occupied ones). For IFR=1,
C              contributions from frozen orbitals (occupied as well as
C              virtual) are included.
C
C     NB! IFR = 1 will be ignored for non-tot. symmetric Lambda matrices
C         and IFR=0 is returned. (=> don't pass an integer constant for IFR).
C
#include <implicit.h>
      DIMENSION XLAMDP(*), XLAMDH(*), DENPH(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
#include <dummy.h>

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CC_PHDEN')

      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)

C     Lambda(p,h) occupied contributions.
C     -----------------------------------

      DO ISYMK = 1,NSYM

         ISYMA = MULD2H(ISYMK,ISYMLP)
         ISYMB = MULD2H(ISYMK,ISYMLH)

         NTOTA = MAX(NBAS(ISYMA),1)
         NTOTB = MAX(NBAS(ISYMB),1)

         KOFFP = IGLMRH(ISYMA,ISYMK) + 1
         KOFFH = IGLMRH(ISYMB,ISYMK) + 1
         KOFFD = IAODIS(ISYMA,ISYMB) + 1

         CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NRHF(ISYMK),
     &              ONE,XLAMDP(KOFFP),NTOTA,XLAMDH(KOFFH),NTOTB,
     &              ZERO,DENPH(KOFFD),NTOTA)

      ENDDO

C     Lambda(p,h) virtual contributions.
C     ----------------------------------

      DO ISYMC = 1,NSYM

         ISYMA = MULD2H(ISYMC,ISYMLP)
         ISYMB = MULD2H(ISYMC,ISYMLH)

         NTOTA = MAX(NBAS(ISYMA),1)
         NTOTB = MAX(NBAS(ISYMB),1)

         KOFFP = IGLMVI(ISYMA,ISYMC) + 1
         KOFFH = IGLMVI(ISYMB,ISYMC) + 1
         KOFFD = IAODIS(ISYMA,ISYMB) + 1
      
         CALL DGEMM('N','T',NBAS(ISYMA),NBAS(ISYMB),NVIR(ISYMC),
     &              ONE,XLAMDP(KOFFP),NTOTA,XLAMDH(KOFFH),NTOTB,
     &              ONE,DENPH(KOFFD),NTOTA)

      ENDDO

C     Include frozen orbitals, if requested.
C     --------------------------------------

      IF ((FROIMP.OR.FROEXP) .AND. (IFR.EQ.1)) THEN

         IF ((ISYMLP.NE.1) .OR. (ISYMLH.NE.1)) THEN

            WRITE(LUPRI,'(/,5X,A,/,5X,A,/)')
     &   SECNAM,': WARNING: Frozen orbital contributions requested for',
     &      'non-tot. symmetric Lambda matrices. Request is ignored....'
            IFR = 0

         ELSE

C           Read MO coefficients.
C           ---------------------

            IF (LWORK .LT. NLAMDS) THEN
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            LUSIR = -1
            CALL GPOPEN(LUSIR,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &                  .FALSE.)
            REWIND LUSIR

            CALL MOLLAB('TRCCINT ',LUSIR,LUPRI)
            READ(LUSIR)
            READ(LUSIR)
            READ(LUSIR) (WORK(I), I = 1,NLAMDS)

            CALL GPCLOSE(LUSIR,'KEEP')

C           Frozen occupied orbitals.
C           -------------------------

            KOFF1 = 0
            DO ISYMK = 1,NSYM

               ISYMA = ISYMK
               ISYMB = ISYMK

               DO KK = 1,NRHFFR(ISYMK)

                  K = KFRRHF(KK,ISYMK)

                  DO B = 1,NBAS(ISYMB)
                     DO A = 1,NBAS(ISYMA)

                        NAK = KOFF1 + NBAS(ISYMA)*(K - 1) + A
                        NBK = KOFF1 + NBAS(ISYMB)*(K - 1) + B
                        NAB = IAODIS(ISYMA,ISYMB) + NBAS(ISYMA)*(B - 1)
     &                      + A

                        DENPH(NAB) = DENPH(NAB) + WORK(NAK)*WORK(NBK)

                     ENDDO
                  ENDDO

               ENDDO

               KOFF1 = KOFF1 + NBAS(ISYMK)*NORBS(ISYMK)

            ENDDO

C           Frozen virtual orbitals.
C           ------------------------

            KOFF1 = 0
            DO ISYMC = 1,NSYM

               ISYMK = ISYMC
               ISYMA = ISYMC
               ISYMB = ISYMC

               DO KC = 1,NVIRFR(ISYMC)

                  C = NRHFS(ISYMK) + KFRVIR(KC,ISYMC)

                  DO B = 1,NBAS(ISYMB)
                     DO A = 1,NBAS(ISYMA)

                        NAC = KOFF1 + NBAS(ISYMA)*(C - 1) + A
                        NBC = KOFF1 + NBAS(ISYMB)*(C - 1) + B
                        NAB = IAODIS(ISYMA,ISYMB) + NBAS(ISYMA)*(B - 1)
     &                      + A

                        DENPH(NAB) = DENPH(NAB) + WORK(NAC)*WORK(NBC)

                     ENDDO
                  ENDDO

               ENDDO

               KOFF1 = KOFF1 + NBAS(ISYMC)*NORBS(ISYMC)

            ENDDO

         ENDIF

      ENDIF

      RETURN
      END
C  /* Deck cc_fopval */
      SUBROUTINE CC_FOPVAL(DENAO,WORK,LWORK)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate and print electronic expectation values of
C              one-electron operators.
C
#include <implicit.h>
      DIMENSION DENAO(*), WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfop.h>
#include <ccrspprp.h>
#include <priunit.h>
#include <dummy.h>

      CHARACTER*9 SECNAM 
      PARAMETER (SECNAM = 'CC_FOPVAL')

      INTEGER IDFAIL(MAFOP)

      CHARACTER*8  LABEL, LBFAIL(MAFOP)
      CHARACTER*27 MESSG(2)

      PARAMETER (ZERO = 0.0D0)

      DATA MESSG /'Hermiticity not determined.',
     &            'Symmetry not determined.   '/

C     Return if nothing to do.
C     ------------------------

      IF (NAFOP .LE. 0) RETURN

C     Initializations.
C     ----------------

      NTFAIL = 0
      CALL IZERO(IDFAIL,NAFOP)

C     Print title and header of table.
C     --------------------------------

      CALL AROUND('Electronic Expectation Value of Requested Operators')

      WRITE(LUPRI,'(16X,A,16X,A)')
     & 'Operator','Electronic Expectation Value'
      WRITE(LUPRI,'(5X,A,4X,A,1X,A,14X,A)')
     & 'Label','Hermiticity','Symmetry','Atomic Units'
      WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'

C     Allocation.
C     -----------

      KONEP = 1
      KEND  = KONEP + N2BASX
      LWRK  = LWORK - KEND + 1

      IF (LWRK .LT. 0) CALL QUIT('Insufficient memory in '//SECNAM)

      DO IOP = 1,NAFOP

         LABEL = PRPLBL(IAFOP(IOP))

C        Read integrals. Get the symmetry and
C        Hermiticity of the operator.
C        ------------------------------------

         IERR = -1
         CALL CCPRPAO(LABEL,WORK(KONEP),ISYM,IHERM,IERR,
     &                WORK(KEND),LWRK)

C        Check for I/O and symmetry errors.
C        ----------------------------------

         IF (IERR .GT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A)')
     &      SECNAM,': I/O error ocurred in CCPRPAO'
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR = ',IERR
            CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
         ELSE IF (IERR .LT. 0) THEN
            NTFAIL = NTFAIL + 1
            LBFAIL(NTFAIL) = LABEL
            IDFAIL(NTFAIL) = 2
            IF (IPRINT .GT. 15) THEN
               WRITE(LUPRI,'(/,5X,A,A,A)')
     &         SECNAM,': irrep not determined for operator ',LABEL
               WRITE(LUPRI,'(5X,A,I10)')
     &         'CCPRPAO returned IERR = ',IERR
               WRITE(LUPRI,'(5X,A,A,/)')
     &         SECNAM,' cowardly refuses to compute expectation value!'
            ENDIF
            GO TO 999
         ENDIF

C        Calculate expectation value.
C        ----------------------------

         IF (IHERM .EQ. -1) THEN

C           Real anti-Hermitian operator.
C           -----------------------------

            VAL = ZERO

         ELSE IF (IHERM .EQ. 0) THEN

C           Real general operator.
C           Is there any such thing implemented in the integral program?
C           - in any case, here's the section for dealing appropriately
C             with this situation (i.e. calculate expectation value of
C             Hermitian part only).
C           ------------------------------------------------------------

            IF (ISYM .EQ. 1) THEN
               DO ISYMB = 1,NSYM
                  ISYMA = ISYMB
                  DO B = 1,NBAS(ISYMB)
                     DO A = 1,B
                        KAB = KONEP + IAODIS(ISYMA,ISYMB)
     &                      + NBAS(ISYMA)*(B - 1) + A - 1
                        KBA = KONEP + IAODIS(ISYMB,ISYMA)
     &                      + NBAS(ISYMB)*(A - 1) + B - 1
                        WORK(KAB) = WORK(KAB) + WORK(KBA)
                        WORK(KBA) = WORK(KAB)
                     ENDDO
                  ENDDO
               ENDDO
               VAL = 0.5D0 * DDOT(N2BST(1),WORK(KONEP),1,DENAO,1)
            ELSE
               VAL = ZERO
            ENDIF

         ELSE IF (IHERM .EQ. 1) THEN

C           Real Hermitian operator.
C           ------------------------

            IF (ISYM .EQ. 1) THEN
               VAL = DDOT(N2BST(1),WORK(KONEP),1,DENAO,1)
            ELSE
               VAL = ZERO
            ENDIF

         ELSE

            NTFAIL = NTFAIL + 1
            LBFAIL(NTFAIL) = LABEL
            IDFAIL(NTFAIL) = 1
            IF (IPRINT .GT. 15) THEN
               WRITE(LUPRI,'(/,5X,A,A,A)')
     &         SECNAM,': Hermiticity not determined for operator ',LABEL
               WRITE(LUPRI,'(5X,A,I10)')
     &         'CCPRPAO returned IHERM = ',IHERM
               WRITE(LUPRI,'(5X,A,A,/)')
     &         SECNAM,' cowardly refuses to compute expectation value!'
            ENDIF
            GO TO 999

         ENDIF

C        Print result.
C        -------------

         WRITE(LUPRI,'(5X,A8,4X,I3,9X,I2,12X,1P,D22.15)')
     &   LABEL,IHERM,ISYM,VAL

C        Escape point for undetermined symmetry/Hermiticity.
C        ---------------------------------------------------

  999    CONTINUE

      ENDDO

C     Finish table.
C     -------------

      WRITE(LUPRI,'(5X,A)')
     & '---------------------------------------------------------------'

C     Failure reports.
C     ----------------

      IF (NTFAIL .GT. 0) THEN
         WRITE(LUPRI,'(/,5X,A,A,/,6X,A,I3,A)')
     &  SECNAM,' cowardly refused to compute expectation values of the',
     &   'following',NTFAIL,' operators:'
         DO IFAIL = 1,NTFAIL
            WRITE(LUPRI,'(5X,A8,A,A)')
     &      LBFAIL(IFAIL),': ',MESSG(IDFAIL(IFAIL))
         ENDDO
      ENDIF

      CALL FLSHFO(LUPRI)

      RETURN
      END
C  /* Deck cc_fopcom */
      SUBROUTINE CC_FOPCOM(DENAO,WORK,LWORK,MODEL)
C
C     Thomas Bondo Pedersen, April 2003.
C
C     Purpose: Calculate and print electronic expectation values of
C              commutators of one-electron operators.
C
#include <implicit.h>
      DIMENSION DENAO(*), WORK(LWORK)
      CHARACTER*10 MODEL
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfop.h>
#include <ccrspprp.h>
#include <priunit.h>
#include <dummy.h>

      CHARACTER*9 SECNAM
      PARAMETER (SECNAM = 'CC_FOPCOM')

      INTEGER IDFAIL(MABFOP)

      CHARACTER*8  LABELA, LABELB, LBFAIL(MABFOP,2)
      CHARACTER*37 MESSG(3)

      PARAMETER (ZERO = 0.0D0)

      DATA MESSG /'Hermiticity not determined.          ',
     &            'Symmetry not determined (operator A).',
     &            'Symmetry not determined (operator B).'/

C     Return if nothing to do.
C     ------------------------

      IF (NABFOP .LE. 0) RETURN

C     Initializations.
C     ----------------

      NTFAIL = 0
      CALL IZERO(IDFAIL,NABFOP)

C     Print title and header of table.
C     --------------------------------

      CALL AROUND('Electronic Expectation Value of Requested'
     &            //' Commutators')
      WRITE(LUPRI,'(9X,A,10X,A,6X,A)')
     & 'Operator A','Operator B','Electronic Expectation Value'
      WRITE(LUPRI,'(5X,A,1X,A,8X,A)')
     & 'Label    Herm. Sym.','Label    Herm. Sym.','Atomic Units'
      WRITE(LUPRI,'(5X,A,A)')
     & '--------------------------------------------------------',
     & '------------'

C     Allocation.
C     -----------

      KDNPH  = 1
      KLAMDP = KDNPH  + MAX(N2BST(1),NT1AM(1))
      KLAMDH = KLAMDP + NGLMDT(1)
      KEND1  = KLAMDH + NGLMDT(1)
      LWRK1  = LWORK  - KEND1 + 1

      KINTA = KDNPH + N2BST(1)
      KINTB = KINTA + N2BASX
      KEND2 = KINTB + N2BASX
      LWRK2 = LWORK  - KEND2 + 1

      KT1AM = KDNPH

      IF ((LWRK1.LT.0) .OR. (LWRK2.LT.0)) THEN
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Read 0th order singles amplitudes.
C     ----------------------------------

      IOPT = 1
      CALL CC_RDRSP('R0 ',1,1,IOPT,MODEL,WORK(KT1AM),DUMMY)

C     Calculate Lambda matrices.
C     --------------------------

      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &            WORK(KEND1),LWRK1)

C     Calculate particle-hole density.
C     IFR = 1: Include frozen orbital contributions.
C     ----------------------------------------------

      IFR = 1
      CALL CC_PHDEN(WORK(KLAMDP),WORK(KLAMDH),WORK(KDNPH),
     &              WORK(KEND1),LWRK1,1,1,IFR)

C     Compute and print expectation value.
C     ------------------------------------

      DO IABFOP = 1,NABFOP

         LABELA = PRPLBL(IACFOP(IABFOP))
         LABELB = PRPLBL(IBCFOP(IABFOP))

         IERRA = -1
         CALL CCPRPAO(LABELA,WORK(KINTA),ISYMA,IHERA,IERRA,
     &                WORK(KEND2),LWRK2)
         IF (IERRA .GT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      SECNAM,': I/O error ocurred in CCPRPAO reading ',LABELA
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR = ',IERRA
            CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
         ELSE IF (IERRA .LT. 0) THEN
            NTFAIL = NTFAIL + 1
            LBFAIL(NTFAIL,1) = LABELA
            LBFAIL(NTFAIL,2) = LABELB
            IDFAIL(NTFAIL) = 2
            IF (IPRINT .GT. 15) THEN
               WRITE(LUPRI,'(/,5X,A,A,A)')
     &         SECNAM,': irrep not determined for operator ',LABELA
               WRITE(LUPRI,'(5X,A,I10)')
     &         'CCPRPAO returned IERR = ',IERRA
               WRITE(LUPRI,'(5X,A,A,/)')
     &         SECNAM,' cowardly refuses to compute expectation value!'
            ENDIF
            GO TO 1000
         ENDIF

         IERRB = -1
         CALL CCPRPAO(LABELB,WORK(KINTB),ISYMB,IHERB,IERRB,
     &                WORK(KEND2),LWRK2)
         IF (IERRB .GT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      SECNAM,': I/O error ocurred in CCPRPAO reading ',LABELB
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'CCPRPAO returned IERR = ',IERRB
            CALL QUIT('I/O error in CCPRPAO ('//SECNAM//')')
         ELSE IF (IERRB .LT. 0) THEN
            NTFAIL = NTFAIL + 1
            LBFAIL(NTFAIL,1) = LABELA
            LBFAIL(NTFAIL,2) = LABELB
            IDFAIL(NTFAIL) = 3
            IF (IPRINT .GT. 15) THEN
               WRITE(LUPRI,'(/,5X,A,A,A)')
     &         SECNAM,': irrep not determined for operator ',LABELB
               WRITE(LUPRI,'(5X,A,I10)')
     &         'CCPRPAO returned IERR = ',IERRB
               WRITE(LUPRI,'(5X,A,A,/)')
     &         SECNAM,' cowardly refuses to compute expectation value!'
            ENDIF
            GO TO 1000
         ENDIF

         IERR = 0
         CALL CC_COMEXP(WORK(KINTA),LABELA,ISYMA,IHERA,
     &                  WORK(KINTB),LABELB,ISYMB,IHERB,
     &                  WORK(KDNPH),1,DENAO,1,
     &                  COM,WORK(KEND2),LWRK2,IERR)
         IF (IERR .NE. 0) THEN
            NTFAIL = NTFAIL + 1
            LBFAIL(NTFAIL,1) = LABELA
            LBFAIL(NTFAIL,2) = LABELB
            IDFAIL(NTFAIL) = 1
            GO TO 1000
         ENDIF

       WRITE(LUPRI,'(5X,A8,1X,I3,3X,I2,3X,A8,1X,I3,3X,I2,6X,1P,D22.15)')
     &   LABELA,IHERA,ISYMA,LABELB,IHERB,ISYMB,COM

C        Escape point for failure.
C        -------------------------

 1000    CONTINUE

      ENDDO

C     Print end of table.
C     -------------------

      WRITE(LUPRI,'(5X,A,A)')
     & '--------------------------------------------------------',
     & '------------'

C     Failure reports.
C     ----------------

      IF (NTFAIL .GT. 0) THEN
         WRITE(LUPRI,'(/,5X,A,A,/,6X,A,I3,A)')
     &  SECNAM,' cowardly refused to compute expectation values of the',
     &   'following',NTFAIL,' commutators:'
         DO IFAIL = 1,NTFAIL
            WRITE(LUPRI,'(5X,A,A8,A,A8,A,A)')
     &      '[',LBFAIL(IFAIL,1),',',LBFAIL(IFAIL,2),']: ',
     &      MESSG(IDFAIL(IFAIL))
         ENDDO
      ENDIF

      CALL FLSHFO(LUPRI)

      RETURN
      END
