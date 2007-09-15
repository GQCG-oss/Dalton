C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
#ifdef OLD_REV_LOG
!===========================================================================
!NOTE (961114/hjaaj): if solvent is implemented in MP2 then
!   IF (DOMP2) FLAG(16) = .FALSE.
!in sirctl.F should be removed.
!/* Comdeck mp2log */
!951130-hjaaj
!removed IADR2 from comdeck INFMP2; moved calculation of IADR2 to new
!  subroutine PHADR2, to be called from RESPONSE.
!corrected dimension for initialization if IADR1.
!renamed PHADDR to PHADR1
!950517-hjaaj
!new call to READMO
!941115-hjaaj
!PHPACK: rewritten code completely (old version used 10% of CPU time on
!  Cray for SOPPA!)
!941028-hjaaj
!MP2CTL: JTRLVL = -10 for SOPPA (was JTRLVL = 10)
!941005-hjaaj
!FCKIFR: new routine for "not frozen" index arrays
!MP2FCK: use NRHF instead of IFRMP2 so MP2SET can be called after FCKFRO
!--- Major overhaul of SOPPA code and other code, in particular FCKFRO.
!December 1993-mjp (Martin J. Packer)
!--- Revised module for SOPPA
!===========================================================================
#endif
C  /* Deck mp2ctl */
      SUBROUTINE MP2CTL(WRK,LFRSAV)
C
C MP2 routines written by Poul Joergensen and Hans Joergen Aa. Jensen
C June 1987.
C MODIFIED OCT 1989 TO INCLUDE FROZEN ORBITALS.
C 011293 mjp Correlation coefficients saved to LUSIFC if SOPPA follows.
C
C PURPOSE:
C  1)  CALCULATE ONE ELECTRON DENSITY MATRIX CORRECT THROUGH
C      SECOND ORDER
C  1A) CALCULATE SECOND ORDER ENERGY
C  2)  DETERMINE THE NATURAL ORBITALS
C
#include <implicit.h>
C
      DIMENSION WRK(LFRSAV)
C
      LOGICAL FCKFRO, MP2FRZ
      PARAMETER ( D0 = 0.0D0 )
C
C Used from common blocks:
C     INFORB: NORBT,NRHF,?
C     SCBRHF: NFRRHF()
C     INFMP2: NPHSYM, MP2FRO
C     INFINP: FLAG(25)
C     INFIND: ISX, ISMO
C     INFTRA: USEDRC
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infopt.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <infmp2.h>
#include <scbrhf.h>
#include <infinp.h>
#include <infind.h>
#include <inftra.h>
#include <dftcom.h>
C
      CALL QENTER('MP2CTL')
      CALL GETTIM(TIMMP2,DUMMY)
      IF (IPRMP2 .GT. 0) WRITE (LUW4,'(//A//A/A)')
     * ' >>>>> Output from SIRIUS MP2 module <<<<<'
     *,' Reference: H.J.Aa.Jensen, P.Jorgensen, H.Agren, and J.Olsen,'
     *,'            J. Chem. Phys. 88, 3834 (1988); 89, 5354 (1988)'
C
C
      KFRSAV= 1
      KFREE = KFRSAV
      LFREE = LFRSAV
      CALL MEMGET('REAL',KORBEN,NORBT ,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDONE ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WRK,KFREE,LFREE)
      CALL MEMGET('INTE',KIADR1,MAXVIR*MAXOCC,WRK,KFREE,LFREE)
C
C     READ IN MO'S AND CHECK FOR CANONICAL HARTREE FOCK ORBITALS
C     Use DONE as scratch for FC in MP2FCK.
C
      KFC = KDONE
      CALL MP2FCK(FCKFRO,EMY,WRK(KCMO),WRK(KFC),WRK(KORBEN),
     &            WRK(KFREE),LFREE)
C     CALL MP2FCK(FCKFRO,EMY,CMO,FC,ORBEN,SCRA,LSCRA)
C
C
C     DETERMINE NUMBER OF PARTICLES AND HOLES AND SET UP AN
C     INDEX ARRAY FOR PARTICLES AND HOLES
C
C     Variables MP2FRO, IPHORD, IFRMP2, INDXPH, NP, NH placed in INFMP2
C
C     MP2SET called after MP2FCK because if (FCKFRO) then
C     the number of frozen MP2 orbitals may have been changed.
C
      CALL MP2SET(WRK(KIADR1))
      MP2FRZ = MP2FRO
C     Disable "MP2FRO" if result not to be printed
      IF (IPRMP2 .LE. 3) MP2FRO = .FALSE.
      CALL FLSHFO(LUW4)
      IF (LUPRI.NE.LUW4) CALL FLSHFO(LUPRI)
C
cmjp If response calculation then make space for the correlation coefficients.
      IF (FLAG(25)) THEN
         CALL MEMGET('REAL',KKAPPA,NPHTOT(1),WRK,KFREE,LFREE)
         CALL DZERO(WRK(KKAPPA),NPHTOT(1))
         JTRLVL = -10
         USEDRC = .TRUE.
      ELSE
         CALL MEMGET('REAL',KKAPPA,0,WRK,KFREE,LFREE)
         JTRLVL = 5
      ENDIF
C
C
C     PERFORM INTEGRAL TRANSFORMATION
C     JTRLVL =  5 : code for MP2 transformation.
C     JTRLVL = 10 : code for full transformation if SOPPA response follows.
C     JTRLVL negative: delete AO integral file after transformation.
C
      CALL TRACTL(JTRLVL,WRK(KCMO),WRK(KFREE),LFREE)
C     CALL TRACTL(ITRLVL,CMO,WRK,LFREE)
C
C
      IF (MP2FRO) THEN
         CALL MEMGET('REAL',KDONEF,N2ORBX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KCMO1 ,NCMOT ,WRK,KFREE,LFREE)
         CALL MEMGET('INTE',KIFRHF,NORBT ,WRK,KFREE,LFREE)
         CALL FCKIFR(WRK(KIFRHF))
      ELSE
         CALL MEMGET('REAL',KDONEF,0,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KCMO1 ,0,WRK,KFREE,LFREE)
         CALL MEMGET('INTE',KIFRHF,0,WRK,KFREE,LFREE)
      END IF
      CALL MEMGET('REAL',KH2M  ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2D  ,N2ORBX,WRK,KFREE,LFREE)
C
C     DETERMINE NUMBER OF PARTICLE-HOLE DIRAC DISTRIBUTIONS WHICH CAN BE
C     KEPT IN CORE IN ONE LOAD
C
C     We assume max 4*N2ORBX additional work space needed :
      NPSTP  = (LFREE-4*N2ORBX)/(NP*NH*NH)
C     Make load per batch even :
      IF (NPSTP .GT. 0) THEN
         NBATCH = (NP-1)/NPSTP  + 1
         NPSTP  = (NP-1)/NBATCH + 1
         NBATCH = (NP-1)/NPSTP  + 1
      ELSE
         NBATCH = 0
      END IF
      IF ( IPRMP2 .GE. 6 .OR. NPSTP .LE. 0) THEN
         WRITE (LUW4,'(3(/A,I8))')
     &      ' Number of PARTICLE-HOLE distributions ',NP*NH,
     &      ' Number of BATCHES                     ',NBATCH,
     &      ' Number of P-H distributions each BATCH',NPSTP*NH
      IF ( NPSTP .LE. 0 ) THEN
         WRITE (LUW4,'(/A)')
     &      ' ERROR: Insufficient work memory for one p-h distribution.'
         CALL ERRWRK('MP2CTL Dirac p-h distributions',(-NP*NH*NH),LFREE)
      END IF
      END IF
      CALL MEMGET('REAL',KH2PHD,NPSTP*NH*NH*NP,WRK,KFREE,LFREE)
C
      EMP2   = D0
      EMP2FL = D0
      CALL DZERO(WRK(KDONE),N2ORBX)
      IF (MP2FRO) THEN
         CALL DCOPY(NCMOT,WRK(KCMO),1,WRK(KCMO1),1)
         CALL DZERO(WRK(KDONEF),N2ORBX)
      END IF
C
      IPLOW  = 1
C
C     ... DO WHILE (IPLOW .LE. NP)
C
 100  CONTINUE
         IPHIGH  =  MIN(NP,IPLOW+NPSTP-1)
C
C        Retrieve DIRAC particle-hole integral distributions
C        for this load
C
         CALL MP2PHD(WRK(KH2PHD),WRK(KH2M),IPLOW,IPHIGH,WRK,KFREE,LFREE)
C        CALL MP2PHD(H2PHD,H2CD,IPLOW,IPHIGH,WRK,KFREE,LFREE)
C
C        Determine second order contributions to one electron
C        density matrix
C
         CALL MP2MUL(IPLOW,IPHIGH, WRK(KDONE),WRK(KKAPPA),WRK(KDONEF),
     *               WRK(KORBEN),WRK(KH2PHD),WRK(KH2M),WRK(KH2D),
     *               WRK(KIFRHF),EMP2,EMP2FL,WRK(KIADR1),WRK,
     *               KFREE,LFREE)
C        CALL MP2MUL(IPLOW,IPHIGH,DONE,COEMP2,DONEFL,ORBEN,H2PHD,
C    *               H2M,H2D,IFRRHF,EMP2,EMP2FL,IADR1,WRK,KFRSAV,LFRSAV)
         IPLOW   =  IPHIGH + 1
      IF (IPLOW.LE.NP) GO TO 100
C     ... END DO WHILE
C
C     Modify DONE to obtain correct one electron density matrix
C     ( factors are defined in MP2DEN ) and add zeroth order
C     contributions
C
      IF ( IPRMP2 .GE. 20 ) THEN
         WRITE(LUPRI,'(/A)')
     *     ' MP2 one-electron density matrix, not modified with factors'
         CALL OUTPUT(WRK(KDONE),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      CALL MP2FAC(WRK(KDONE),WRK(KORBEN))
      IF ( IPRMP2 .GE. 10 ) THEN
         WRITE(LUPRI,'(/A)')
     *     ' Final MP2 one-electron density matrix'
         CALL OUTPUT(WRK(KDONE),1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      IF (IPRMP2 .GE. 1) THEN
         IF (MP2FRZ) THEN
            WRITE (LUW4,'(/A,8I5)')
     *      ' MP2 frozen orbitals per symmetry :',(NFRMP2(I),I=1,NSYM)
            WRITE (LUW4,'(/A)')
     *      ' Total energy and NO:s calculated with frozen orbitals.'
         END IF
         IF (.NOT.DODFT) THEN
            WRITE (LUW4,'(/,3(/A,F25.10))')
     *       ' Hartree-Fock total energy :',EMY,
     *       ' + MP2 contribution        :',EMP2,
     *       ' = MP2 second order energy :',EMY+EMP2
         ELSE
            WRITE (LUW4,'(/,(/A,F25.10),(/A,F25.10,A,F8.5,A,F15.10),
     *        (/A,F25.10))')
     *       ' Kohn--Sham contribution   :',EMY,
     *       ' + scaled PT2 contribution :',EMP2,' *',WDFTMP,
     *       ' = ',EMP2*WDFTMP,
     *       ' = KS-PT2 energy           :',EMY+WDFTMP*EMP2
            EMCSCF = EMY + WDFTMP*EMP2
         END IF
      END IF
      ECORR = EMY + EMP2
C
C     write the correlation coefficients and second order one-density to
C     the interface file LUSIFC
C
      IF (FLAG(25)) THEN
         CALL MP2IFC(WRK(KDONE),WRK(KKAPPA))
      END IF
      CALL NATORB(IPRMP2,WRK(KDONE),WRK(KCMO),WRK(KFREE),LFREE)
C
C     Save MP2 orbitals for MCSCF/CI if not SOPPA
C
      IF (.NOT.FLAG(25)) THEN
         CALL NEWORB('MP2SAVE ',WRK(KCMO),.TRUE.)
      END IF
C
      IF (MP2FRO) THEN
C     If (mp2fro) then also calculate full MP2 energy and NO occupations
C     If iprmp2 .le. 0 then mp2fro will be false
         CALL MP2FAC(WRK(KDONEF),WRK(KORBEN))
         IF ( IPRMP2 .GE. 10 ) THEN
            IF ( FCKFRO ) THEN
               WRITE(LUPRI,'(/A)')
     &       ' MP2 one-electron density matrix with RHF frozen orbitals'
            ELSE
               WRITE(LUPRI,'(/A)')
     &       ' MP2 one-electron density matrix with no frozen orbitals'
            END IF
            CALL OUTPUT(WRK(KDONEF),1,NORBT,1,NORBT,NORBT,NORBT,
     &         1,LUPRI)
         END IF
C
         IF (FCKFRO) THEN
            WRITE (LUW4,'(/A,8I5//A)')
     &      ' RHF frozen orbitals per symmetry :',(NFRRHF(I),I=1,NSYM),
     &      ' Total MP2 energy and NO occupations calculated with'//
     &      ' frozen orbitals from RHF.'
         ELSE
            WRITE (LUW4,'(/A)')
     *      ' Full MP2 energy and NO occupations obtained without'//
     *      ' frozen orbitals.'
         END IF
         WRITE (LUW4,'(/A/3(/A,F25.10))')
     *     ' NOTE that the orbitals saved on disk are those obtained'//
     *     ' with frozen orbitals (see above).',
     *    ' Hartree-Fock total energy :',EMY,
     *    ' + MP2 contribution        :',EMP2FL,
     *    ' = MP2 second order energy :',EMY+EMP2FL
         CALL NATORB(IPRMP2,WRK(KDONEF),WRK(KCMO1),WRK(KFREE),LFREE)
      END IF
C     END IF (MP2FRO .AND. .NOT.FCKFRO)
C
      CALL GETTIM(TIMEND,DUMMY)
      TIMMP2 = TIMEND - TIMMP2
      WRITE (LUERR,1860) TIMMP2
      IF (IPRMP2 .GE. 2) WRITE (LUW4,1860) TIMMP2
 1860 FORMAT (/' Time used for MP2 natural orbitals :',F12.3,
     *        ' CPU seconds.')
      CALL MEMREL('MP2CTL',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL FLSHFO(LUERR)
      CALL FLSHFO(LUW4)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('MP2CTL')
      RETURN
      END
C  /* Deck natorb */
      SUBROUTINE NATORB(IPRNO,DONE,CMO,WRK,LFRSAV)
C
C PURPOSE: DETERMINE NATURAL ORBITALS
C
C INPUT  : DONE IS ONE ELECTRON DENSITY MATRIX
C
C
#include <implicit.h>
      DIMENSION DONE(NORBT,NORBT),CMO(*),WRK(*)
C
      PARAMETER ( D0 = 0.0D0 )
C
C Used from common blocks:
C
C INFIND : ISSMO
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
C
      CALL QENTER('NATORB')
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
      CALL MEMGET('REAL',KORBNO,N2ORBT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KORBEN,NORBT ,WRK,KFREE,LFREE)
      CALL MEMGET('WORK',KWRK1 ,LWRK1 ,WRK,KFREE,LFREE)
      IJ     = 0
      DO 150 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
         IORBI = IORB(ISYM)
         DO 160 I = 1,NORBI
            DO 170 J = 1,I
               IJ = IJ + 1
               WRK(KORBNO-1+IJ) = DONE(IORBI+J,IORBI+I)
 170        CONTINUE
 160     CONTINUE
 150  CONTINUE
C
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
         IORBI = IORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         NBASI = NBAS(ISYM)
         ISSYM = IIORB(ISYM)
         ICSYM = ICMO(ISYM)
C
         CALL JACO(WRK(KORBNO+ISSYM),CMO(1+ICSYM),NORBI,NORBI,
     *             NBASI,WRK(KWRK1),WRK(KWRK1+NORBI))
C        CALL JACO(F,VEC,NB,NMAX,NROWV,BIG,JBIG)
C
         OCCSUM = D0
         DO 175 I=1,NORBI
            II = ISSYM + (I*I+I)/2
            WRK(KORBEN-1+IORBI+I)=WRK(KORBNO-1+II)
            OCCSUM = OCCSUM + WRK(KORBEN-1+IORBI + I)
  175    CONTINUE
#if !defined (VAR_MP2SUPSYM)
         CALL ORDER2(CMO(1+ICSYM),WRK(KORBEN+IORBI),NORBI,NBASI)
#else
C     HJ : implement SUPSYM in MP2 ????? 940823 hjaaj
C          Done, Aug. 2004 /hjaaj
         CALL ORD2SS(CMO(1+ICSYM),WRK(KORBEN+IORBI),
     &               ISSMO(1+IORBI),NORBI,NBASI)
#endif
         IF (IPRNO .GE. 1) THEN
            OCCRHF = 2 * NRHF(ISYM)
#if !defined (VAR_MP2SUPSYM)
            WRITE (LUW4,1840) ISYM
            WRITE (LUW4,1850) (WRK(KORBEN-1+IORBI+I),I=1,NORBI)
#else
            IF (.NOT.SUPSYM) THEN
               WRITE (LUW4,1840) ISYM
               WRITE (LUW4,1850) (WRK(KORBEN-1+IORBI+I),I=1,NORBI)
            ELSE
               WRITE (LUW4,1841) ISYM
               WRITE (LUW4,1851)
     &            (WRK(KORBEN-1+IORBI+I),ISSMO(IORBI+I),I=1,NORBI)
            END IF
#endif
            WRITE (LUW4,1860) OCCSUM,OCCRHF,OCCSUM-OCCRHF
         END IF
  200 CONTINUE
#if defined (VAR_MP2SUPSYM)
      IF (SUPSYM) CALL AVEORD()
C     ... remake ISSORD() as ISSMO() may have changed in ORD2SS
#endif
C
 1840 FORMAT(/' Natural orbital occupation numbers, symmetry',I2/)
 1850 FORMAT(4X,5F15.8)
#if defined (VAR_MP2SUPSYM)
 1841 FORMAT(/' Natural orbital occupation numbers, symmetry',I2,
     &        ' (supersymmetry in parenthesis)')
 1851 FORMAT(4(F14.9,' (',I2,')'))
#endif
 1860 FORMAT(/' Sum',F15.8,/' RHF',F15.8,/' Diff.',F13.8)
C
      CALL MEMREL('NATORB',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('NATORB')
      RETURN
C     ... end of NATORB
      END
C  /* Deck mp2den */
      SUBROUTINE MP2DEN(NFROZ,IFROZ,IAH,IMP,IAHMP,ISYMAM,H2PHD,DONE,
     *                  COEMP2,IADR1,ORBEN,H2M,H2D,EMP2,KFRO)
C
C 18-Oct-1989 PJ
C THE PURPOSE OF THIS ROUTINE IS THE SAME AS OLD MP2DEN EXCEPT THAT
C SOME OF THE ORBITALS IN THE MP2 CALCULATIONS MAY BE FROZEN
C AND THEREFORE WILL NOT CONTRIBUTE TO THE MP2 ENERGY AND DENSITY
C MATRIX.  NFROZ(NSYM) specifies the number of frozen orbitals
C per symmetry and IFROZ(NORBT) is an index array telling if
C a specific orbital is frozen or not.
C 011293 mjp KFRO indicates if this is the first or second call
C            to MP2DEN during a frozen calculation.
C
C CALCULATE SECOND ORDER PERTURBATION CONTRIBUTION TO ONE ELECTRON
C DENSITY MATRIX WHICH ORIGINATE FROM <**,MA> AND (**,MA)
C
C  A,B,C,D OCCUPIED AND M,N,P,Q UNOCCUPIED HF ORBITAL INDICES
C
C    D(M,A) =  2/(E(A)-E(M) *
C                ( SUM(P,Q,C) (CP,MQ)*[QA,CP] -
C                  SUM(P,C,D) (CP,DA)*[MD,CP] )
C
C    D(A,B) = -2 * SUM(P,Q,C) (1/(E(A)+E(C)-E(P)-E(Q))*(AQ,CP)*[QB,CP]
C
C    D(M,N) =  2 * SUM(P,C,D) (1/(E(C)+E(D)-E(P)-E(N))*(DN,CP)*[MD,CP]
C
C      WHERE
C
C   [MA,BN] = (1/(E(A)+E(B)-E(M)-E(N)) * (2*(MA,BN)-(MB,AN))
C           = (1/(E(A)+E(B)-E(M)-E(N)) * (2*(MA,BN)-<MA,BN>)
C
C ( FACTORS BETWEEN = AND * ARE MULTIPLIED ONTO D(*,*) IN MP2FAC )
C
C INPUT:
C  MULLIKEN INTEGRAL DISTRIBUTION (**,MA) = H2M(*,*)
C  DIRAC PARTICLE HOLE INTEGRAL DISTRIBUTION H2PHD(*,*,IAHMP)
C
C OUTPUT:
C  SECOND ORDER CONTRIBUTION ADDED IN DONE(*,*)
C
#include <implicit.h>
C
      DIMENSION NFROZ(8),IFROZ(NORBT),H2PHD(NP,NH,*)
      DIMENSION DONE(NORBT,*),ORBEN(*),H2M(NORBT,*),H2D(NORBT,*)
      DIMENSION COEMP2(*)
C
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0 )
C
C Used from common blocks :
C    INFINP : FLAG(25)
C
#include <maxorb.h>
#include <priunit.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <infinp.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
C
      SAVE    TVALUE
      DATA    TVALUE/D0/
C
      CALL QENTER('MP2DEN')
C
C CONSTRUCT   (2*(MA,BN)-(MB,AN))/(E(A)+E(B)-E(M)-E(N)) =
C             (2*(MA,BN)-<MA,BN>)/(E(A)+E(B)-E(M)-E(N))
C
      JMP  = IPHORD(NH+IMP)
      JAH  = IPHORD(IAH)
      ENMA = - ORBEN(JMP) + ORBEN(JAH)
      IF (IPRMP2 .GE. 30) THEN
         WRITE (LUPRI,'(//A/,2(/A,2I5))')
     *      ' >>> TEST OUTPUT FROM MP2DEN',
     *      ' particle indices JMP, IMP     : ',JMP,IMP,
     *      ' hole     indices JAH, IAH     : ',JAH,IAH
         IF (IPRMP2 .GE. 40) WRITE (LUPRI,'(A,2I5,3(/A,1P,D11.4))')
     *      ' IAHMP index, symmetry(A M)    : ',IAHMP,ISYMAM,
     *      ' EMP2 before this distribution :', EMP2,
     *      ' Orbital energy of MP orbital  :', ORBEN(JMP),
     *      ' Orbital energy of AH orbital  :', ORBEN(JAH)
      END IF
      EMP2BF = EMP2
      DO 100 IH = 1,NH
         IBH = IPHORD(IH)
      IF ( IFROZ(IBH) .EQ. 0 ) THEN
         DO 200 IP = 1,NP
            INP = IPHORD(NH+IP)
            ENMANB = ENMA - ORBEN(INP) + ORBEN(IBH)
            H2D(INP,IBH) = ( D2*H2M(INP,IBH) - H2PHD(IP,IH,IAHMP))
     *                      / ENMANB
            IF (IPRMP2 .GE. 30) THEN
               TVALUE = TVALUE + H2D(INP,IBH)*H2D(INP,IBH)
               IF (IPRMP2 .GE. 50 .AND.
     *             H2D(INP,IBH).NE.D0 .AND. H2M(INP,IBH).NE.D0) THEN
                  WRITE (LUPRI,'(/A,4I5,/A,1P,3D15.6)')
     *            ' IH,IBH,IP,INP         :',IH,IBH,IP,INP,
     *            ' ENMANB, a(MA,NB), H2M :',
     *            ENMANB,H2D(INP,IBH),H2M(INP,IBH)
               END IF
            END IF
C
C Accumulate MP2 energy :
C
            EMP2 = EMP2 + H2M(INP,IBH) * H2D(INP,IBH)
C
C PARTICLE HOLE PART OF MULLIKEN DISTRIBUTIONS ARE DIVIDED WITH ORBITAL
C ENERGY DIFFERENCES TO BE ABLE TO CARRY OUT MATRIX MULTIPLICATIONS.
C (THE MODIFIED MATRIX ELEMENTS ENTER ONLY IN DOUBLES CONTRIBUTIONS,
C  NOT IN SINGLES)
C
            H2M(INP,IBH) = H2M(INP,IBH)/ENMANB
            H2M(IBH,INP) = H2M(INP,IBH)
 200     CONTINUE
      END IF
 100  CONTINUE
cmjp Pack the correlation coefficients into COEMP2 if SOPPA follow
      IF (FLAG(25) .AND. KFRO .GE. 0) THEN
         CALL PHPACK(H2D,COEMP2(IADR1(IMP,IAH)+1),ISYMAM,1,-1)
      ENDIF
C
C
C  [*,*] = H2D(*,*)
C  (*,*) = H2M(*,*)
C
C ADD CONTRIBUTION TO DONE (FACTORS ARE TAKEN CARE OF LATER)
C
C   DONE(A,M) = (A,D)*[M,D]
C   DONE(M,N) = (N,D)*[M,D]
C
      DO 1000 ISSYM = 1,NSYM
         IOFFI = IORB(ISSYM) + 1
         NORBI = NORB(ISSYM)
         NRHFI = NRHF(ISSYM)
         NFRMPI = NFROZ(ISSYM)
         NVIRI = NORBI - NRHFI
         JSYM  = MULD2H(ISYMAM,ISSYM)
         IOFFJ = IORB(JSYM) + 1
         NRHFJ = NRHF(JSYM)
         NFRMPJ = NFROZ(JSYM)
         NVIRJ = NORB(JSYM) - NRHFJ
C
C   DONE(A,B) = (A,Q)*[Q,B]
C   DONE(M,B) = (M,Q)*[Q,B]
C        I,I     I,J   J,I
C
         CALL DGEMM('N','N',NORBI-NFRMPI,NRHFI-NFRMPI,NVIRJ,1.D0,
     &              H2M(IOFFI+NFRMPI,IOFFJ+NRHFJ),NORBT,
     &              H2D(IOFFJ+NRHFJ,IOFFI+NFRMPI),NORBT,1.D0,
     &              DONE(IOFFI+NFRMPI,IOFFI+NFRMPI),NORBT)
C
C   DONE(A,M) = (A,D)*[M,D]
C   DONE(M,N) = (N,D)*[M,D]
C        I,I     I,J   I,J
C
         CALL DGEMM('N','T',NORBI-NFRMPI,NVIRI,NRHFJ-NFRMPJ,1.D0,
     &              H2M(IOFFI+NFRMPI,IOFFJ+NFRMPJ),NORBT,
     &              H2D(IOFFI+NRHFI,IOFFJ+NFRMPJ),NORBT,1.D0,
     &              DONE(IOFFI+NFRMPI,IOFFI+NRHFI),NORBT)
 1000 CONTINUE
      IF (IPRMP2 .GE. 30) THEN
         WRITE (LUPRI,'(1P,2(/A,2D11.4))')
     *   ' del(EMP2), EMP2 after this distribution :',EMP2-EMP2BF,EMP2,
     *   ' T = sum (a(MA,NB)**2) after this distr. :',TVALUE
      END IF
      CALL QEXIT('MP2DEN')
      RETURN
C     ... return from MP2DEN
      END
C  /* Deck mp2fac */
      SUBROUTINE MP2FAC(DONE,ORBEN)
C
C MODIFY DONE TO OBTAIN CORRECT ONE ELECTRON DENSITY MATRIX
C ( FACTORS ARE DEFINED IN MP2DEN )
C
#include <implicit.h>
C
      DIMENSION DONE(NORBT,*),ORBEN(*)
C
      PARAMETER (D0 = 0.0D0 , D2 = 2.0D0 ,TOLELE = 1.0D-4 )
C
#include <priunit.h>
#include <inforb.h>
#include <infpri.h>
C
      CALL QENTER('MP2FAC')
C
C     Scale with factor 2 given in formulas in comments in MP2DEN
C
      CALL DSCAL(N2ORBX,D2,DONE,1)
      ELERMV = D0
      ELERCV = D0
      DO 1500 ISYM = 1,NSYM
         IORBI = IORB(ISYM)
         NRHFI = NRHF(ISYM)
         NORBI = NORB(ISYM)
         DO 1600 I = (IORBI+1),(IORBI+NRHFI)
            DO 1700 J = (IORBI+1+NRHFI),(IORBI+NORBI)
               DONE(J,I) = ( DONE(J,I) - DONE(I,J))/(ORBEN(I)-ORBEN(J) )
               DONE(I,J) = DONE(J,I)
 1700       CONTINUE
            DO 1800 J = (IORBI+1),(IORBI+NRHFI)
               DONE(J,I) = -DONE(J,I)
 1800       CONTINUE
            ELERMV = ELERMV - DONE(I,I)
C
C           ADD ZERO'TH ORDER CONTRIBUTION
C
            DONE(I,I) = D2 + DONE(I,I)
 1600    CONTINUE
         DO 1900 IUNOCC = (IORBI+1+NRHFI),(IORBI+NORBI)
            ELERCV = ELERCV + DONE(IUNOCC,IUNOCC)
 1900    CONTINUE
 1500 CONTINUE
C
      IF ( ABS(-ELERMV+ELERCV).GT.TOLELE ) THEN
         NWARN = NWARN + 1
         WRITE(LUW4,'(/A,/A,/A,F12.6,A)') ' *******WARNING*********',
     *   ' MP2 CALCULATION CARRIED OUT INCORRECTLY ',
     *   ' TRACE OF ONE ELECTRON DENSITY MATRIX GIVE ',
     *   (-ELERMV+ELERCV),' ELECTRONS  IN ADDITION TO SPECIFIED NUMBER'
         I = 1
      ELSE
         I = 0
      END IF
      IF (ELERMV .LT. 0 .OR. ELERMV .GT. 2*NRHFT) THEN
         NWARN = NWARN + 1
         WRITE(LUW4,'(/A/A,I5/A,F12.6/A/A)')
     &   ' *******WARNING from MP2*********',
     &   ' Number of electrons in Hartree-Fock orbitals :',2*NRHFT,
     &   ' Number of electrons removed from HF orbitals :',ELERMV,
     &   ' Probable error: this molecule is not suitable for MP2',
     &   ' because E(LUMO) - E(HOMO) small or negative.'
         I = 1
      ELSE
         I = 0
      END IF
      IF (IPRMP2 .GE. 2 .OR. I .EQ. 1) THEN
         WRITE(LUW4,'(/A,F11.6,A)')
     *      ' MP2 move',ELERCV,' electrons to unoccupied HF orbitals'
         IF (IPRMP2 .GE. 6 .OR. I .EQ. 1)
     &      WRITE(LUW4,'(/A,F11.6,A/A,1P,D10.2,A)')
     *      ' MP2 MOVE',ELERMV,' ELECTRONS FROM OCCUPIED HF ORBITALS',
     *      ' MP2 CHANGE IN TRACE OF ONE ELECTRON DENSITY MATRIX :',
     *      -ELERMV+ELERCV,' ELECTRONS.'
      END IF
C
      CALL QEXIT('MP2FAC')
      RETURN
C
C     ... END OF MP2FAC
C
      END
C  /* Deck mp2fck */
      SUBROUTINE MP2FCK(FCKFRO,EMY,CMO,FC,ORBEN,SCRA,LSCRA)
C
C Purpose:
C  CONSTRUCT INACTIVE FOCK MATRIX AND CHECK IF ORBITALS ARE
C  CANONICAL HARTREE-FOCK ORBITALS
C
C Input:
C  CMO; initial molecular orbitals used to build Fock matrix,
C       assumed to be orthonormal.
C
C Output:
C  EMY; Fock energy
C  CMO; molecular orbitals diagonalizing Fock matrix
C
C Scratch:
C  FC; the inactive Fock matrix and scratch area for overlap matrix
C  SCRA; general scratch area
C
#include <implicit.h>
      DIMENSION CMO(*),FC(*),ORBEN(*),SCRA(*)
      LOGICAL   FCKFRO
C
C
      PARAMETER (D0 = 0.0D0,    EMYCNV = 1.D-4,
     &           DBIG = 1.D+12, GAPMIN = 0.1D0)
#include <dummy.h>
C
C  Used from common blocks:
C     INFINP : POTNUC
C     INFORB : NRHF(*), NRHFT, ...
C     SCBRHF : IOPRHF,NFRRHF(*)
C     INFMP2 : NFRMP2(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <scbrhf.h>
#include <infmp2.h>
#include <infind.h>
#include <inftap.h>
#include <infpri.h>
C
      CALL QENTER('MP2FCK')
      IF (IOPRHF .GT. 0) THEN
         WRITE (LUW4,'(//A)')
     *   ' FATAL ERROR in MP2 module: MP2 is not implemented for'//
     *   ' open-shell Hartree-Fock'
         CALL QUIT('ERROR: MP2 is not implemented for open-shell RHF')
      END IF
      CALL READMO(CMO,9,SCRA,LSCRA)
C
      IF (IPRMP2 .GE. 3) WRITE (LUW4,'(///A//A,I5/A,8I5)')
     *   ' Check that orbitals are canonical HARTREE-FOCK orbitals',
     *   ' Number of electrons :',2*NRHFT,
     *   ' Orbital occupations :',(NRHF(I),I=1,NSYM)
C
C
      EMY    = D0
C
C        Step 1: Construct inactive Fock matrix
C
      CALL FCKMAT(.TRUE.,DUMMY,CMO,EMY,FC,DUMMY,SCRA,LSCRA)
C     CALL FCKMAT(ONLYFC,DV,CMO,EMY,FC,FV,WRK,LFREE)
      IF (IPRMP2 .GE. 3) WRITE (LUW4,1730) EMY
      EMY = EMY + POTNUC
      IF (IPRMP2 .GE. 3) WRITE (LUW4,1732) EMY
C
C         Step 2: Check if orbitals are canonical orbitals
C
      CALL DZERO(ORBEN,NORBT)
      FCKFRO = .FALSE.
      VALCON =  D0
      EHOMO  = -DBIG
      ELUMO  =  DBIG
      DO 200 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 200
         NFRHFI= NFRRHF(ISYM)
         IF (NFRHFI .GT. 0) FCKFRO = .TRUE.
         IF (NFRHFI .GT. NFRMP2(ISYM)) THEN
            NWARN = NWARN + 1
            WRITE (LUW4,'(3(/A,I5))')
     *      ' ***** WARNING from MP2FCK for symmetry',ISYM,
     *      ' Number of frozen MP2 orbitals this symmetry',NFRMP2(ISYM),
     *      ' is reset to number of frozen RHF orbitals  ',NFRHFI
            NFRMP2(ISYM) = NFRHFI
         END IF
         NRHFI = NRHF(ISYM)
         IORBI = IORB(ISYM)
         ISSYM = IIORB(ISYM)
         DO 300 I = NFRHFI+1,NORBI
            IJOFF = ISSYM + IROW(I)
            DO 400 J = NFRHFI+1,I-1
               VALCON = VALCON + ABS(FC(IJOFF+J))
400         CONTINUE
            ORBEN(IORBI+I) = FC(IJOFF+I)
            IF (I .LE. NRHFI) THEN
               EHOMO = MAX(EHOMO,ORBEN(IORBI+I))
            ELSE
               ELUMO = MIN(ELUMO,ORBEN(IORBI+I))
            END IF
300      CONTINUE
         IF (IPRMP2 .GE. 3) THEN
            WRITE (LUW4,1740) ISYM
            IF (NFRHFI .GT. 0) WRITE (LUW4,'(A,I3,A)')
     *         ' Note: the first',NFRHFI,' orbitals are frozen'
            WRITE (LUW4,1750) (ORBEN(IORBI+I),I=1,NORBI)
         END IF
200   CONTINUE
      IF (ELUMO-EHOMO .LT. GAPMIN) THEN
         NWARN = NWARN + 1
         WRITE(LUW4,'(/A/A/A)')
     &   ' *******WARNING from MP2*********',
     &   ' This molecule is probably not suitable for MP2',
     &   ' because E(LUMO) - E(HOMO) small or negative.'
         I = 1
      ELSE
         I = 0
      END IF
      IF (IPRMP2 .GE. 1 .OR. I .EQ. 1) THEN
         WRITE (LUW4,'(2(/A,F15.8)/A/A,F15.8)')
     &      '    E(LUMO) :',ELUMO,
     &      '  - E(HOMO) :',EHOMO,
     &      '  --------------------------',
     &      '    gap     :',ELUMO-EHOMO
      END IF
      IF ( VALCON .GT. EMYCNV ) THEN
         WRITE(LUW4,'(/A,//A,F12.8//A)')
     *   ' MOLECULAR ORBITALS ARE NOT CANONICAL HARTREE-FOCK ORBITALS',
     *   ' ABS SUM OF OFF-DIAGONAL FOCK MATRIX ELEMENTS ARE :',VALCON,
     *   ' The Fock matrix :'
         CALL OUTPKB(FC(1),NORB,NSYM,1,LUW4)
         CALL QTRACE(LUW4)
         CALL QUIT('MP2FCK: MO:S ARE NOT CANONICAL HF ORBITALS')
      ELSE IF (IPRMP2 .GE. 2) THEN
         WRITE(LUW4,'(/A,1P,D10.2)')
     *      ' ABS SUM OF OFF-DIAGONAL FOCK CORE ELEMENTS ARE :',VALCON
      END IF
      CALL QEXIT('MP2FCK')
      RETURN
C
 1730 FORMAT(//' Hartree-Fock electronic energy:',F25.12)
 1732 FORMAT( /' Hartree-Fock total      energy:',F25.12)
 1740 FORMAT( /' Hartree-Fock orbital energies, symmetry',I2/)
 1750 FORMAT(4X,5F15.8)
C
C
C *** end of subroutine MP2FCK
C
      END
C  /* Deck mp2ifc */
        SUBROUTINE MP2IFC(DONE,COEMP2)
C
C This routine prepares the SOPPA interface for RESPONSE storing
C the correlation coefficients and the second order one-density
C matrix on the end of the LUSIFC file.
C MJP 181093
C
#include <implicit.h>
#include <priunit.h>
        DIMENSION DONE(*), COEMP2(*)
C Used from common blocks:
C    INFORB: N2ORBT, NRHFT, NVIRT
C    INFMP2: NPHMAX
C    INFTAP: LUSIFC
C
#include <maxorb.h>
#include <dummy.h>
#include <inforb.h>
#include <infmp2.h>
#include <inftap.h>
      CHARACTER*8 LAB123(3), LABEOD, LABELN, LMP2INFO
      DATA LAB123/'********','********','********'/
      DATA LABEOD /'EODATA  '/
      DATA LMP2INFO /'MP2INFO '/
      CALL GETDAT(LAB123(2),LAB123(3))
C     ... place date in LAB123(2) and time in LAB123(3)
C
C Open the file LUSIFC (we must write to the end so check it is open!
      CALL GPOPEN(LUSIFC,FNSIFC,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      CALL MOLLAB(LABEOD,LUSIFC,LUERR)
      BACKSPACE (LUSIFC)
      WRITE (LUSIFC) LAB123,LMP2INFO
C Write the correlation coefficients
      CALL WRITT (LUSIFC,NPHTOT(1),COEMP2)
C
C Write the second order one-density to LUSIFC
C
      CALL WRITT (LUSIFC,N2ORBX,DONE)
      WRITE (LUSIFC) LAB123,LABEOD
      CALL GPCLOSE(LUSIFC,'KEEP')
      RETURN
      END
C  /* Deck mp2mul */
      SUBROUTINE MP2MUL(IPLOW,IPHIGH,DONE,COEMP2,DONEFL,ORBEN,H2PHD,
     *                  H2M,H2D,IFRRHF,EMP2,EMP2FL,IADR1,
     *                  WRK,KFRSAV,LFRSAV)
C
C WRITTEN APRIL-1-87
C   Revision 891017-hjaaj: use NXTH2M
C 011293 mjp call to MP2DEN altered since IPHPORD etc. now in /INFMP2/
C
C PURPOSE:
C     READ IN INTEGRAL DISTRIBUTIONS TO BE USED TO CALCULATE
C     SECOND ORDER CONTRIBUTIONS TO THE ONE ELECTRON DENSITY MATRIX
C
C Scratch:
C     WRK(kfrsav:kfrsav-1+lfrsav)
C
#include <implicit.h>
C
      DIMENSION DONEFL(NORBT,NORBT),DONE(NORBT,NORBT),ORBEN(NORBT)
      DIMENSION H2PHD(NP,NH,*),H2M(NORBT,NORBT),H2D(NORBT,NORBT)
      DIMENSION WRK(*), COEMP2(*), IFRRHF(*)
C
C -- local constants
C
      PARAMETER ( D2 = 2.0D0 )
C
C Used from common blocks:
C   INFORB : N2ORBX, ...
C   SCBRHF : NFRRHF(*)
C   INFMP2 : IPHORD, NP, NH, IFRMP2, MP2FRO
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infpri.h>
#include <scbrhf.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
C     
      LOGICAL CHOLE,DHOLE
#include <line.h>
      DIMENSION NEEDTP(6)
      DATA NEEDTP /0,1,1,1,1,0/
C     ... only particle-hole distributions needed
C         assumed: NRHF(isym) .ge. NISH(isym)
C              and NVIR(isym) .le. NSSH(isym)
#include <line.h>
C
      CALL QENTER('MP2MUL')
      KFREE = KFRSAV
      LFREE = LFRSAV
      IF (MP2FRO) THEN
         CALL MEMGET('REAL',KH2M1 ,N2ORBX,WRK,KFREE,LFREE)
      ELSE
         CALL MEMGET('REAL',KH2M1 ,0,WRK,KFREE,LFREE)
      END IF
      CALL DZERO(H2D,N2ORBX)
C
C ****************************************************************
C     Loop over Mulliken distributions allowed in NEEDTP(6)
C
      IDIST = 0
  100 CALL NXTH2M(IC,ID,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
C
         CHOLE = ( INDXPH(IC) .GT. 0 )
         DHOLE = ( INDXPH(ID) .GT. 0 )
C
C        Check if we need this distribution:
C        HOLE-HOLE AND PARTICLE-PARTICLE  DISTRIBUTIONS NOT NEEDED
C
      IF ( CHOLE .EQV. DHOLE ) GO TO 100
C
         IF (DHOLE) THEN
            ISWAP = IC
            IC    = ID
            ID    = ISWAP
         END IF
C
C        IC IS A HOLE INDEX
C
         IDPA = - INDXPH(ID)
      IF ( IDPA .LT. IPLOW .OR. IDPA .GT. IPHIGH ) GO TO 100
C     ... we do not need this distribution this time
C
         ICHO  = INDXPH(IC)
         ICHDP = ( IDPA - IPLOW ) * NH + ICHO
C
         ISYMC = ISMO(IC)
         ISYMD = ISMO(ID)
         ISYMCD= MULD2H(ISYMC,ISYMD)
C
C        MULLIKEN DISTRIBUTION (**,CD) IN H2M(*,*)
C        DIRAC PARTICLE-HOLE DISTRIBUTION IN H2PHD(*,*,ICHPH) C<D
C
C        ADD SECOND ORDER CONTRIBUTIONS FROM INTEGRAL DISTRIBUTIONS CD
C        TO ONE ELECTRON DENSITY MATRIX; copy H2M to WRK(KH2M1) because
C        H2M is modified in MP2DEN and it will be needed below if
C        MP2FRO true.
C
         IF (MP2FRO) CALL DCOPY(N2ORBX,H2M,1,WRK(KH2M1),1)
         IF ( IFRMP2(IC).EQ.0 ) THEN
C        If ( IC not frozen   ) then
            CALL MP2DEN(NFRMP2,IFRMP2,ICHO,IDPA,ICHDP,ISYMCD,H2PHD,DONE,
     *                  COEMP2,IADR1,ORBEN,H2M,H2D,EMP2,0)
C           CALL MP2DEN(NFROZ,IFROZ,ICHO,IDPA,ICHDP,ISYMCD,H2PHD,DONE,
C     *                 COEMP2,IADR1,ORBEN,H2M),H2D,EMP2,KFRO)
         END IF
C
C        IF MP2FRO, that is, if some orbitals are frozen above,
C        then also collect full MP2 density matrix and energy.
C        IFRRHF is defined to correspond to only RHF frozen orbitals.
C
Cmjp  If the last argument of MP2DEN < 0, then we will form the energy
Cmjp  in MP2DEN assuming no frozen orbitals. If the argument is
Cmjp  >= 0 then we take account of any frozen orbitals.
C
         IF ( MP2FRO ) THEN
         IF ( IFRRHF(IC) .EQ. 0 ) THEN
            CALL MP2DEN(NFRRHF,IFRRHF,ICHO,IDPA,ICHDP,ISYMCD,H2PHD,
     *                  DONEFL,COEMP2,IADR1,ORBEN,WRK(KH2M1),H2D,
     *                  EMP2FL,-1)
         END IF
         END IF
C
C        Go get next Mulliken distribution
C
      GO TO 100
C
C     arrive at 800 when finished with all Mulliken distributions
C
  800 CONTINUE
C
      CALL MEMREL('MP2MUL',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('MP2MUL')
      RETURN
      END
C  /* Deck fckifr */
      SUBROUTINE FCKIFR(IFRRHF)
C
C  941005-hjaaj
C  Set IFRRHF(NORBT) to 1 if orbital is frozen in RHF
C  otherwise to 0
C
#include <implicit.h>
      DIMENSION IFRRHF(NORBT)
C
C Used from common blocks:
C   INFORB : NSYM, NORBT, IORB(8)
C   SCBRHF : NFRRHF(8)
C
#include <inforb.h>
#include <scbrhf.h>
C
      CALL IZERO(IFRRHF,NORBT)
      DO 100 ISYM = 1,NSYM
         IF (NFRRHF(ISYM) .GT. 0) THEN
            IORBI = IORB(ISYM)
            DO 200 I = IORBI+1,IORBI+NFRRHF(ISYM)
               IFRRHF(I) = 1
  200       CONTINUE
         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck mp2phd */
      SUBROUTINE MP2PHD(H2PHD,H2CD,IPLOW,IPHIGH,WRK,KFRSAV,LFRSAV)
C
C (retrieve Particle-Hole integrals in Dirac format)
C Written  18. May 1987 by Hans Jorgen Aa. Jensen.
C Revisions:
C   891017-hjaaj: use NXTH2M
C
C Purpose:
C   For MP-2 density matrix we need particle-hole integrals in Dirac
C   format.
C   Sort Mulliken integrals (i a / j b) = ( A B / C D )
C   to   Dirac    format    <i b / a j> = < A D / B C >
C   (in H2PHD(i,b;a,j) ; PH for particle-hole, D for Dirac)
C
C Scratch:
C   WRK(kfrsav:kfrsav-1+lfrsav)
C   H2CD(NORBT,NORBT)
C
#include <implicit.h>
      DIMENSION H2PHD(NP,NH,*), H2CD(NORBT,NORBT)
      DIMENSION WRK(*)
C
C Used from common blocks:
C   INFORB : NSYM,NRHF(*),NORB(*),IORB(*),?
C   INFIND : ISMO(*),?
C   INFMP2 : NP,NH,?
C   INFPRI : ?
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infmp2.h>
#include <infpri.h>
C
      LOGICAL AHOLE, BHOLE, CHOLE, DHOLE
#include <line.h>
      DIMENSION NEEDTP(6)
      DATA NEEDTP /0,1,1,1,1,0/
C     ... only particle-hole distributions needed
C         assumed: NRHF(isym) .ge. NISH(isym)
C              and NVIR(isym) .le. NSSH(isym)
#include <line.h>
C
      CALL QENTER('MP2PHD')
C
C     >>> zero H2PHD
C
      LENGTH = NP * NH * NH * (IPHIGH - IPLOW + 1)
      CALL DZERO(H2PHD,LENGTH)
C
C
C ****************************************************************
C     Loop over Mulliken distributions allowed in NEEDTP(6)
C
      KFREE = KFRSAV
      LFREE = LFRSAV
      IDIST = 0
  100 CALL NXTH2M(IC,ID,H2CD,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
C
         ICSYM = ISMO(IC)
         IDSYM = ISMO(ID)
         CHOLE = ( INDXPH(IC) .GT. 0 )
         DHOLE = ( INDXPH(ID) .GT. 0 )
C
C        See if we need this distribution
C        (We only need particle-hole distributions)
C
      IF (CHOLE .EQV. DHOLE) GO TO 100
C     ... we do not need this distribution
         IF (DHOLE) THEN
            ISWAP = IC
            IC    = ID
            ID    = ISWAP
         END IF
C        ... now c is hole, d is particle
C            Do we want this distribution this time ?
         IDPA   = - INDXPH(ID)
      IF ( IDPA .LT. IPLOW .OR. IDPA .GT. IPHIGH ) GO TO 100
C     ... No, we don't
         ICHO   = INDXPH(IC)
         IDPOFF = (IDPA - IPLOW) * NH
         ICDSYM = MULD2H(ICSYM,IDSYM)
C
C        loop over ( A B / C D ) of ( p h / h p ) type
C          sort as < A C / B D > in / B D > distribution.
C                ( < A C / B D > is thus < p h / h p > type )
C
         DO 300 IA = 1,NORBT
            AHOLE = ( INDXPH(IA) .GT. 0 )
            IF (AHOLE) THEN
               IAH   =   INDXPH(IA)
               IASYM = ISMO(IA)
               IBSYM = MULD2H(IASYM, ICDSYM)
               IBST  = IORB(IBSYM) + NRHF(IBSYM) + 1
               IBEND = IORB(IBSYM) + NORB(IBSYM)
               DO 200 IB = IBST,IBEND
                  IBP = - INDXPH(IB)
                  H2PHD(IBP,ICHO,IDPOFF+IAH) = H2CD(IB,IA)
  200          CONTINUE
            END IF
  300    CONTINUE
C
C        Go get next Mulliken distribution
C
      GO TO 100
C
C     arrive at 800 when finished with all Mulliken distributions
C
  800 CONTINUE
C
C *** end of subroutine MP2PHD
C
      CALL MEMREL('MP2PHD',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('MP2PHD')
      RETURN
      END
C  /* Deck mp2set */
      SUBROUTINE MP2SET(IADR1)
C
C 19-May-1987 hjaaj
C 20-MAY-1987 PJ     INSERTED IPHORD
Cold  SUBROUTINE MP2SET(MP2FRO,INDXPH,IPHORD,IFRMP2,NP,NH)
C 011293 mjp MP2FRO,IPHORD,INDXPH,IFRMP2,NP,NH now put into /INFMP2/
C
C Calculate NP, number of particles,
C           NH, number of holes, and
C           INDXPH(NORBT), particle or hole index for each orbital.
C           IPHORD(NORBT), GIVE SYMMETRY ORDERED INDEX FOR PARTICLE OR
C                          HOLE INDEX
C           IFRMP2(NORBT), .EQ.1 FOR FROZEN MP2 ORBITAL ELSE .EQ.0
C
#include <implicit.h>
C
C Used from common blocks:
C  INFORB : NSYM, NRHF(NSYM), NVIR(NSYM), NFRMP2(NSYM), ...
C  INFMP2 : IPHORD, MP2FRO, INDXPH, IFRMP2,NP,NH
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infmp2.h>
      DIMENSION IADR1(MAXVIR,MAXOCC)
#include <infind.h>
C
      CALL QENTER('MP2SET')
C
      NH = NRHFT
      NP = NORBT - NRHFT
C
      IPH = 0
      IH  = 0
      IP  = 0
      IORBT = 0
      CALL IZERO(IFRMP2,NORBT)
      IMP2 = 0
      DO 300 ISYM = 1,NSYM
         NRHFI = NRHF(ISYM)
         NVIRI = NVIR(ISYM)
         NFRMPI= NFRMP2(ISYM)
         DO 100 I = 1,NRHFI
            IPH = IPH + 1
            IH  = IH  + 1
            IPHORD(IH)  = IPH
            INDXPH(IPH) = IH
            IF (I .LE. NFRMPI) THEN
               IMP2        = IMP2 + 1
               IFRMP2(IPH) = 1
            END IF
  100    CONTINUE
         DO 200 I = 1,NVIRI
            IPH = IPH + 1
            IP  = IP  + 1
            IPHORD(NH+IP) = IPH
            INDXPH(IPH)   = -IP
  200    CONTINUE
  300 CONTINUE
C
C     Any frozen orbital in MP2?
C
      MP2FRO = (IMP2 .GT. 0)
C
C Calculate the number of ph excitations of each symmetry, ignoring
C orbitals frozen in MP2.
C
         CALL IZERO(NPHSYM,NSYM)
         DO I = 1,NH
            II = IPHORD(I)
            IF (IFRMP2(II).EQ.1) GOTO 500
            DO J = 1,NP
               JA = IPHORD(NH+J)
               KSYM = MULD2H(ISMO(II),ISMO(JA))
               NPHSYM(KSYM) = NPHSYM(KSYM) + 1
            ENDDO
 500     CONTINUE
         ENDDO
         NPHMAX = NPHSYM(1)
         DO I = 2,NSYM
            NPHMAX = MAX(NPHMAX,NPHSYM(I))
         ENDDO
C
C Determine ph offset array IADR1(VIRT,OCC) for symmetry 1
C
         CALL IZERO(NPHTOT,NSYM)
         CALL PHADR1(IADR1,1)
C
      CALL QEXIT('MP2SET')
      RETURN
      END
C  /* Deck phpack */
      SUBROUTINE PHPACK (XK, XKP, JPHSYM, IWAY, IPHLAB)
C
C Pack or unpack ph-array for MP2/SOPPA.
C
C If IWAY.LT.0 then unpack the linear array XKP into XK. Else pack the
C   the array XK into XKP. 
C If IPHLAB>0, pack XK(OCC,VIRT); IPHLAB<0 C pack XK(VIRT,OCC),
C   and IPHLAB=0, pack both (XK(OCC,VIRT) and XK(VIRT,OCC))
C
#include <implicit.h>
      DIMENSION XK(NORBT,NORBT), XKP(*)
C
C Used from common blocks:
C   INFORB: NSYM,NORBT,MULD2H, N2ORBX , NOCC(), ...
C   INFMP2: NFRMP2(*)
C
#include <maxorb.h>
#include <inforb.h>
#include <infmp2.h>
C
      IF (IWAY .LT. 0) CALL DZERO(XK,N2ORBX)
C
      IA = 0
      DO ISYM = 1,NSYM
         JSYM = MULD2H(ISYM,JPHSYM)
         IF (NOCC(ISYM) .EQ. 0) GO TO 600
         IIST = IORB(ISYM) + NFRMP2(ISYM) + 1
         IIEND= IORB(ISYM) + NOCC(ISYM)
         JAST = IORB(JSYM) + NOCC(JSYM) + 1
         JAEND= IORB(JSYM) + NOCC(JSYM) + NSSH(JSYM)
         IF (IWAY .GE. 0) THEN
            IF (IPHLAB .GE. 0) THEN
               DO II = IIST,IIEND
                  DO JA = JAST, JAEND
                     IA = IA + 1
                     XKP(IA) = XK(II,JA)
                  END DO
               END DO
            ELSE
               DO II = IIST,IIEND
                  DO JA = JAST, JAEND
                     IA = IA + 1
                     XKP(IA) = XK(JA,II)
                  END DO
               END DO
            END IF
         ELSE
            IF (IPHLAB .GE. 0) THEN
               DO II = IIST,IIEND
                  DO JA = JAST,JAEND
                     IA = IA + 1
                     XK(II,JA) = XKP(IA)
                  END DO
               END DO
               IF (IPHLAB .EQ. 0) THEN
                  DO II=IIST,IIEND
                     DO JA=JAST,JAEND
                        XK(JA,II) = XK(II,JA)
                     ENDDO
                  ENDDO
               ENDIF
            ELSE
               DO II = IIST,IIEND
                  DO JA = JAST,JAEND
                     IA = IA + 1
                     XK(JA,II) = XKP(IA)
                  END DO
               END DO
            END IF
         END IF
 600  CONTINUE
      END DO
      RETURN
      END
C  /* Deck phadr1 */
      SUBROUTINE PHADR1(IADR1,KSYMOP)
C
C Written by Martin Packer 090394.
C
C IADR1(A,I) is an address C array giving the offsets to
C   the set of ph excitations associated with A,I.
C
#include <implicit.h>
C
C
#include <maxorb.h>
#include <maxash.h>
#include <infind.h>
#include <inforb.h>
#include <infmp2.h>
C
      INTEGER A, AA
      DIMENSION IADR1(MAXVIR,MAXOCC)
C
C
      CALL ISETVC(IADR1,-100 000 000,MXCORB*MAXOCC)
C     call isetvc(ivector,ivalue,ndim)
      IAB = 0
      NPHTOT(KSYMOP) = 0
      DO II = 1,NH
        I = IPHORD(II)
        IF (IFRMP2(I) .EQ. 0) THEN
           DO AA = 1,NP
              A = IPHORD(NH+AA)
              IASYM = MULD2H(ISMO(I),ISMO(A))
              JBSYM = MULD2H(IASYM,KSYMOP)
              NPHTOT(KSYMOP) = NPHTOT(KSYMOP) + NPHSYM(JBSYM)
              IADR1(AA,II) = IAB
              IAB = IAB + NPHSYM(JBSYM)
           ENDDO
         END IF
      ENDDO
      RETURN
      END
C  /* Deck phadr2 */
      SUBROUTINE PHADR2(IADR2,NCOL,NROW)
C
C 30-Nov-1995 Hans Joergen Aa. Jensen.
C
C The address of the MP2 first order 2p2h coefficients
C   kappa(ia,jb) er IADR1(B,J) + IADR2(A,I)
C
C Based on code written by Erik Kragh Dalskov fall 1994.
C
#include <implicit.h>
      DIMENSION IADR2(NCOL,NROW)
C
#include <maxorb.h>
#include <maxash.h>
#include <infind.h>
#include <inforb.h>
#include <infmp2.h>
C
      INTEGER A,AA,ASYM,A1OFF,AST
C
      IF (NOCCT .GT. NROW .OR. NSSHT .GT. NCOL) THEN
         CALL QUIT('prog.error: dimension of IADR2 in PHADR2 too small')
      END IF
C
      CALL ISETVC(IADR2,-100 000 000,NCOL*NROW)
C     call isetvc(ivector,ivalue,ndim)
C
      DO IASYM = 1,NSYM
         IA = 0
         DO I = 1,NH
            II = IPHORD(I)
         IF (IFRMP2(II).EQ.1) GOTO 500
            ISYM = ISMO(II)
            ASYM = MULD2H(ISYM,IASYM)
            AST  = IORB(ASYM) + NOCC(ASYM) + 1
            A1OFF= -INDXPH(AST) - 1
            DO A = 1,NSSH(ASYM)
               IADR2(A1OFF+A,I) = IA+A
            ENDDO
            IA = IA + NSSH(ASYM)
 500     CONTINUE
         ENDDO
      ENDDO
      RETURN
      END
C --- end of sirius/sirmp2.F ---
