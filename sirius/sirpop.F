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
C File: sirius/sirpop.F
C
#include <single.h>
C
C  /* Deck sirpop */
      SUBROUTINE SIRPOP(INP_POPKEY,DV,WRK,LFREE)
C
C Nov-10-85 Hans Agren
C l.r. 910708-hjaaj/0710-hjaaj
C
C Purpose :
C        To control calculations of certain properties from SIRIUS
C        wave function
C
C MOTECC-90: This module, SIRPOP, is described in the input/ouput
C            Documentation of MOTECC-90.
C
C (Implemented 10-Nov-85: Mulliken population analysis, Virial analysis
C  Dipole moments, Quadrupole moments. These four subroutines have been
C  taken from the Stockholm CASSCF program (release 82 09 01) and
C  modified for SIRIUS.)
C
C
#include <implicit.h>
      CHARACTER*(*)   INP_POPKEY
      DIMENSION DV(*), WRK(*)
C
#include <iratdef.h>
C
C     MXCGTO = Max number of cgto's contributing to a symmetry orbital
C
      PARAMETER ( MXCGTO = 8 )
C
C   INFINP : FLAG(*)
C   INFORB ? NBAST,
C   INFPOP ?
C   INFTAP : LUIT1
C
#include <maxorb.h>
#include <mxcent.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include "infpop.h"
#include <inftap.h>
#include <infpri.h>
C
C RESTRICTIONS: -THE MAXIMUM NUMBER OF ATOMS IS SET TO MXCENT
C               -The highest L quantum number allowed for is F
C
C               To alter these limits, change the data statements
C               (in MULLIK) and the dimensions in common block INFPOP
C
C
      LOGICAL DIPOLM, QUADRM, MULLKM
      LOGICAL FNDLAB
      CHARACTER*(6) POPKEY
C
C
      CALL QENTER('SIRPOP')
C 
      POPKEY = INP_POPKEY
      IPRMUL_SAVE = IPRMUL
      IF (POPKEY(1:4) .EQ. 'DIIS') THEN
         iprmul = -1
         DIPOLM = .false.
         QUADRM = .false.
         MULLKM = .true.
      ELSE IF (POPKEY(1:6) .EQ. 'MCITER') THEN
         iprmul = -1
         DIPOLM = .false.
         QUADRM = .false.
         MULLKM = .true.
      ELSE IF (POPKEY(1:5) .EQ. 'FINAL') THEN
         WRITE (LUW4,'(//A/)')
     &   ' <<< SIRIUS Population analysis. >>>'
         IF (LUPRI .NE. LUW4) WRITE (LUPRI,'(//A/)')
     &   ' <<< SIRIUS Population analysis. >>>'
         DIPOLM = .NOT.FLAG(71)
         QUADRM = .NOT.FLAG(72)
         MULLKM = .NOT.FLAG(73)
      ELSE
         CALL QUIT('Unrecognized POPKEY : '//POPKEY)
      END IF
C
      IF (NASHT .GT. 0) THEN
         write (LUPRI,*) ' SIRPOP skipped because NASHT > 0'
         go to 9999
      END IF
C
C     *****   core allocation  *****
C
      MATCMO = 0
      DO 10 ISYM = 1,NSYM
         MATCMO = MATCMO + NBAS(ISYM)*(NISH(ISYM)+NASH(ISYM))
   10 CONTINUE
C
C  ***---VIRIAL---POPTRA---***
C
C     The structure of the working area will be
C
C       1) CMO      MO coefficients (occupied only)
C       2) OCCUP    occupation numbers (inactive and active only)
C       3) IAMN     label of the centre to which this cgto belongs
C       4) IPR      index for the type of the cgto (s,x,y,...)
C       5) T OR S   kinetic energy matrix, overlap matrix
C       6) JTRAN    number of cgto's that contribute to the i'th
C                   symmetry orbital
C       7) ITRAN    numbers of the JTRAN(I) cgto's
C       8) CTRAN    coefficients of the cgto's
C       9) CMOCGT   MO coefficients in cgto basis
C      10) SS       overlap matrix in cgto basis
C
      KCMO   = 1
      KOCC   = KCMO   + MATCMO
      KCMO1  = KOCC   + NOCCT
      KAOCC  = KCMO1  + NCMOT
      KWRK   = KAOCC  + NASHT
C
      KIAMN  = KCMO1
      KIPR   = KIAMN  + NBAST
      KOVLP  = KIPR   + NBAST
      KJTRAN = KOVLP  + NNBAST
      KITRAN = KJTRAN + NBAST
      KCTRAN = KITRAN + NBAST*MXCGTO/IRAT
      KWRK   = MAX(KCTRAN + NBAST*MXCGTO,KWRK)
      KCMOCG = KWRK
      KSS    = KCMOCG + NBAST*NOCCT
      KWRK1  = KSS    + NBAST*(NBAST + 1)/2
C
      IF (KWRK1 .GT. KWRK) KWRK = KWRK1
      IF (KWRK .GE. LFREE) CALL ERRWRK('SIRPOP-1',KWRK,LFREE)
C
C     ***---MULLIK---***
C
C     5)  NUMBER   For each cgto basis function number of the
C                  centre to which the function belongs
C     6)  POPCGT   population of each cgto - sum over MO's
C     7)  POPC     population of each cgto - contribution from
C                  each MO separately
C                  this may partly overlap with
C     9')  CMOCGT  as previously
C     10') SS          - " -
C     11') POPNET  net and overlap populations over cgto basis
C                  functions - sum over all MO's
C
      KNUMB  = KOVLP
      KPOPNAB= KNUMB  + NBAST
      KPOPABI= KPOPNAB+ MXCENT*MXCENT
      KPOPCG = KPOPABI+ MXCENT*MXCENT
      KPOPC  = KWRK1
      KPPNET = KPOPC  + NBAST*NOCCT
      KWRK2  = KPPNET + NBAST*(NBAST + 1)/2
      IF (KWRK2 .GE. LFREE) CALL ERRWRK('SIRPOP-2',KWRK2,LFREE)
C
C  ***---DIPMOM---QUADRP---***
C
C     3)  DIPMAT dipole moment matrix
C     3') QUMAT  quadrupole moment matrix
C
      KDIPM  = KIAMN
      KWRK3D = KDIPM  + NNBAST*3
      KQUMAT = KIAMN
      KWRK3Q = KQUMAT + NNBAST*6
      KWRK3  = MAX(KWRK3D,KWRK3Q)
      IF (KWRK3 .GE. LFREE) CALL ERRWRK('SIRPOP-3',KWRK3,LFREE)
C
C     ***** Get the MO coefficients and NO occupation numbers
C
      REWIND LUIT1
      CALL READMO (WRK(KCMO1),9,DUMMY,0)
      IF (NASHT .LE. 1) THEN
         IF (NASHT .EQ. 1) WRK(KAOCC) = NACTEL
      ELSE IF ( FNDLAB('NATOCCAC',LUIT1) ) THEN
         CALL READT(LUIT1,NASHT,WRK(KAOCC))
      ELSE
         NWARN = NWARN + 1
         WRITE (LUW4,'(//A//A/A/A/)')
     &      ' *** Label "NATOCCAC" not found on unit LUIT1 in SIRPOP',
     &      ' WARNING, occupation of active orbitals not available',
     &      '          population analysis aborted.',
     &      '          Use ".TRACI" option.'
         GO TO 9999
      END IF
      IF (MCTYPE .EQ. 2) THEN
         NWARN = NWARN + 1
         WRITE (LUW4,'(/A/A/A)')
     &       ' WARNING, population analysis for RAS wave functions'//
     &       ' uses pseudo-natural orbitals,',
     &       '          that is, the off-diagonal elements of the'//
     &       ' one-electron density matrix which',
     &       '          couple the different RAS spaces are ignored.'
      END IF
C
      IEOCC = 0
      JCMO  = KCMO
      DO 40 ISYM = 1,NSYM
         NOCCI = NOCC(ISYM)
      IF (NOCCI.EQ.0)   GO TO 40
         NBASI = NBAS(ISYM)
         CALL DCOPY(NOCCI*NBASI,WRK(KCMO1+ICMO(ISYM)),1,WRK(JCMO),1)
         JCMO = JCMO + NOCCI*NBASI
C
C
C     Get the occupation numbers for occupied orbitals
C     (2 for inactive and from WRK(KAOCC) for active orbitals)
C
         ISTOCC = IEOCC + 1
         IEOCC  = IEOCC + NOCCI
         ISTASH = ISTOCC+ NISH(ISYM)
         DO 210 I = ISTOCC,(ISTASH-1)
            WRK(KOCC-1+I) = 2
  210    CONTINUE
         NASHI = NASH(ISYM)
         IF (NASHI.GT.0)  THEN
            CALL DCOPY(NASHI,WRK(KAOCC+IASH(ISYM)),1,
     &                       WRK(KOCC -1+ISTASH),1)
         END IF
   40 CONTINUE
C
C     *****     Perform energy analysis   ****
C
      CALL VIRIAL (NNBAST,MXCGTO,MATCMO,
     *             WRK(KCMO),WRK(KOCC),WRK(KIAMN),WRK(KIPR),
     *             WRK(KOVLP),WRK(KJTRAN),WRK(KITRAN),WRK(KCTRAN))
C
C     *****     Transform from symmetry orbital basis to cgto basis
C
      CALL POPTRA (NNBAST,MXCGTO,MATCMO,
     *             WRK(KCMO),WRK(KOVLP),WRK(KJTRAN),WRK(KITRAN),
     *             WRK(KCTRAN),WRK(KCMOCG),WRK(KSS))
C
C     *****     Perform Mulliken population analysis  *****
C
      IF (MULLKM) THEN
         CALL MULLIK (WRK(KOCC),WRK(KIAMN),WRK(KIPR),WRK(KNUMB),
     &                WRK(KPOPNAB),WRK(KPOPABI),
     *                WRK(KPOPCG),WRK(KPOPC),WRK(KCMOCG),WRK(KSS),
     *                WRK(KPPNET))
      END IF
C
C     *****     Dipole moments  *****
C
      IF (DIPOLM) THEN
         CALL DIPMOM (MATDIM,MATCMO,WRK(KCMO),WRK(KOCC),WRK(KDIPM))
      END IF
C
C     *****     Quadrupole moments  *****
C
      IF (QUADRM) THEN
         CALL QUADRP (MATDIM,MATCMO,WRK(KCMO),WRK(KOCC),WRK(KQUMAT))
      END IF
 9999 CALL QEXIT('SIRPOP')
      IPPRMUL = IPRMUL_SAVE
      RETURN
      END
C  /* Deck virial */
      SUBROUTINE VIRIAL (MATDIM,MXCGTO,MATCMO,CMO,OCCUP,IAMN,IPR,
     *                   T,JTRAN,ITRAN,CTRAN)
C
C Purpose:
C      Form the kinetic energy and print an energy analysis
C      (based on RELEASE 82 08 30)
C
#include <implicit.h>
#include <dummy.h>
C
#if defined (VAR_ITRAN2)
      INTEGER*2 ITRAN
#endif
      DIMENSION CMO(MATCMO),OCCUP(NOCCT),T(MATDIM),JTRAN(NBAST),
     *          ITRAN(NBAST,MXCGTO),CTRAN(NBAST,MXCGTO),
     *          IAMN(NBAST),IPR(NBAST)
C
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
C
C
C   INFINP : POTNUC,
C   INFORB : NBAST,
C   INFTAP : LUONEL
C   INFOPT : EMCSCF,
C
#include <maxorb.h>
#include <mxcent.h>
#include <priunit.h>
#include "infpop.h"
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infopt.h>
#include <infpri.h>
C
      DIMENSION NB(8)
      LOGICAL FNDLAB,FOUND
C
C     ***** Read the transformation info from unit LUONEL
C
      CALL GPOPEN(LUONEL,FNONEL,'UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
#if !defined (VAR_FASTIO)
      REWIND LUONEL
      IF ( .NOT. FNDLAB('SCFINP  ',LUONEL)) THEN
#endif
         WRITE (LUW4,'(//A/2A/)')
     *     ' *** Label "SCFINP" not found on unit LUONEL,',
     *     '     basic information for VIRIAL and MULLIKEN routines',
     *     ' is missing.'
         RETURN
#if !defined (VAR_FASTIO)
      END IF
      READ (LUONEL) TITLEM, NS,(NB(I),I=1,NS), DUMM,
     *              KMA,(MNUC,I=1,KMA), NBRO,(JTRAN(I),I=1,NBRO),
     *              ((ITRAN(I,J),J=1,MXCGTO),I=1,NBRO),
     *              ((CTRAN(I,J),J=1,MXCGTO),I=1,NBRO),
     *              IABAS,(IAMN(I),I=1,IABAS),(IPR(I),I=1,IABAS),
     *              DIPOL,NNUC,((COOO,J=1,3),I=1,NNUC),IFXYZ,TPROP,
     *              QPOL,QQFAC,JFXYZ
      CALL GPCLOSE(LUONEL,'KEEP')
C
C     **********     Kinetic energy
C     Read kinetic energy matrix from unit LUONEL
C
      FOUND = .FALSE.
      CALL RDONEL('KINETINT',FOUND,T,NNBAST)
      IF ( .NOT. FOUND) THEN
         WRITE (LUW4,'(//A/A/)')
     *     ' *** Label "KINETINT" not found on unit LUONEL,',
     *     '     basic information for VIRIAL routine is missing.'
         RETURN
      END IF
C
C
C     **********     Compute kinetic energy     **********
C
C              TK = Sum N  Sum C  * C  * <j/T/k>
C                    p   p j,k  pj   pk
C
C     where p goes over occupied mo's and j and k over symmetry
C     orbitals. note that only occupied orbitals are included in
C     the matrix cmo, so no skip over secondary mo's is needed.
C
C     MATBLK indicates where the current symmetry block in T matrix
C            starts.
C     IST    tells where the appropriate line in the T matrix starts.
C
      TK = D0
      IP = 0
      IPSTRT = 0
      MATBLK = 0
C
      DO 50 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         NOCCI = NOCC(ISYM)
C
         DO 40 ICC = 1,NOCCI
            IP  = IP + 1
            IST = MATBLK
C
            DO 30 J = 1,NBASI
               JORB = IPSTRT + J
C
               DO 20 K = 1,J
                  TK = TK + D2*CMO(JORB)*CMO(IPSTRT+K)
     *                    * OCCUP(IP)*T(IST+K)
   20          CONTINUE
C
C              Take care of the case j = k
C
               TK = TK - CMO(JORB)*CMO(JORB)*OCCUP(IP)*T(IST+J)
               IST = IST + J
   30       CONTINUE
C
            IPSTRT = IPSTRT + NBASI
   40    CONTINUE
C
         MATBLK = MATBLK + NBASI*(NBASI + 1)/2
   50 CONTINUE
C
C     **********     Energy analysis and virial theorem     **********
C
      V   = EMCSCF - TK
      VIR = -V/TK
c     WRITE(LUPRI,1001) EMCSCF,TK,V,POTNUC,VIR
c     IF (LUW4 .NE. LUPRI) WRITE(LUW4,1001) EMCSCF,TK,V,POTNUC,VIR
 1001 FORMAT(//' Energy analysis'
     *        /' ---------------'
     *        /' Total energy',T30,F18.10
     *        /' Kinetic energy',T30,F18.10
     *        /' Potential energy',T30,F18.10
     *        /' Nuclear repulsion energy',T30,F18.10
     *        /' -V/T',T30,F14.6)
      RETURN
#endif
      END
C  /* Deck poptra */
      SUBROUTINE POPTRA (MATDIM,MXCGTO,MATCMO,CMO,S,JTRAN,ITRAN,CTRAN,
     *                   CMOCGT,SS)
C
C Purpose: transform the mo coefficients and overlap matrix
C          from symmetry orbital basis to cgto basis.
C          Called from MULLIK
C          Based on release 82 08 31
C
#include <implicit.h>
C
#if defined (VAR_ITRAN2)
      INTEGER*2 ITRAN
#endif
      DIMENSION CMO(MATCMO),S(MATDIM),JTRAN(NBAST),
     *          ITRAN(NBAST,MXCGTO),CTRAN(NBAST,MXCGTO),
     *          CMOCGT(*),SS(*)
C
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C   INFORB ? NBAST,NNBAST,
C
#include <mxcent.h>
#include <priunit.h>
#include <inforb.h>
#include "infpop.h"
#include <infpri.h>
C
      LOGICAL FOUND
C
C     *****   Read in overlap integrals from unit LUONEL
C
      FOUND = .FALSE.
      CALL RDONEL('OVERLAP ',FOUND,S,NNBAST)
      IF ( .NOT. FOUND) THEN
         WRITE (LUW4,'(//A/A/)')
     *     ' *** Label "OVERLAP" not found on unit LUONEL,',
     *     '     basic information for POPTRA is missing.'
         RETURN
      END IF
C
C
C     *****     TRANSFORM THE MO COEFFICIENTS
C
C     for each mo you have nbas(isym) symmetry coefficients (cmo).
C     for each cmo you have a number of cgto coefficients (cmocgt)
C     so the transformation reads
C
C               C'  = C  * D
C                pa    pk   ka
C
C     where p indicates the mo and k the symmetry orbital. a indexes
C     the cgto's. For each k, all values of a as found in ITRAN are
C     examined.
C
      KSTART = 0
      IPSTC  = 0
      IPSTS  = 0
C
      CALL DZERO(CMOCGT,NORBT*NBAST)
      DO 60 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         NOCCI = NOCC(ISYM)
C
         DO 50 ICC = 1,NOCCI
C
            DO 40 K = 1,NBASI
               KK = KSTART + K
               J  = JTRAN(KK)
C
               DO 30 JJ = 1,J
                  IA = ITRAN(KK,JJ)
                  CMOCGT(IPSTC+IA) =
     1               CMO(IPSTS+K)*SIGN(D1,CTRAN(KK,JJ))
   30          CONTINUE
   40       CONTINUE
            IPSTS = IPSTS + NBASI
            IPSTC = IPSTC + NBAST
   50    CONTINUE
         KSTART = KSTART + NBASI
   60 CONTINUE
C
C     *****     Transform the integrals
C
C     Here we have
C
C             LOOP <a/b> = SUM D  * D  * <k/l>
C              a,b         k,l  ka   lb
C
C     The order of the loops is k,a,l,b.
C
      IDIM = NBAST*(NBAST+1)/2
      CALL DZERO(SS,IDIM)
C
      INTBAS = 0
      KLSTRT = 0
C
      DO 120 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
C
         DO 110 K=1,NBASI
            KROW = INTBAS + K*(K - 1)/2
            KK   = KLSTRT + K
            J    = JTRAN(KK)
            FACKA= J
            FACKA= D1/FACKA
C
            DO 100 JJ = 1,J
               IA  = ITRAN(KK,JJ)
               IIA = IA*(IA - 1)/2
C
               DO 90 L = 1,K
                  IKL = KROW + L
                  LL = KLSTRT + L
                  M  = JTRAN(LL)
                  FACLB = M
                  FACLB = D1/FACLB
C
                  DO 80 MM = 1,M
                     IB  = ITRAN(LL,MM)
                  IF (IB.GT.IA) GO TO 80
                     IND = IIA + IB
                     FAC = SIGN(FACKA,CTRAN(KK,JJ))
     1                   * SIGN(FACLB,CTRAN(LL,MM))
                     SS(IND) = SS(IND) + S(IKL)*FAC
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
         KLSTRT = KLSTRT + NBASI
         INTBAS = INTBAS + NBASI*(NBASI + 1)/2
  120 CONTINUE
C
      RETURN
      END
C  /* Deck mullik */
      SUBROUTINE MULLIK (OCCUP,IAMN,IPR,NUMBER,
     *                   POPNAB,POPABI,POPCGT,POPC,CMOCGT,SS,POPNET)
C
C  Purpose:
C         Given the mo coefficients and overlap integrals in cgto basis,
C         to compute and print out a population analysis.
C         Based on release 82 08 31
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
C
      DIMENSION OCCUP(NOCCT),IAMN(NBAST),IPR(NBAST),NUMBER(NBAST),
     &          POPNAB(MXCENT,MXCENT),POPABI(MXCENT,MXCENT),
     &          POPCGT(NBAST),POPC(NBAST,NOCCT),CMOCGT(NBAST,NOCCT),
     &          SS(*),POPNET(*)
C
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D2 = 2.0D0)
C
C Used from common blocks:
C   INFPOP ? MAXTYP, ...
C   INFORB ? NBAST,
C   INFTAP : LUTEMP
C
#include <maxaqn.h>
#include <ccom.h>
#include "infpop.h"
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
C
      DIMENSION   POPTSM(MXCENT)
      CHARACTER*1 SPDCAR
CMAERKE-hj-910708:
C     This information should be read from integral program
C     (HERMIT, MOLECULE, or whatever), because this only works
C     for old MOLECULE order, not e.g. spherical Gaussians.
C     Various matrices in INFPOP should then be allocated acc.
C     to actual dimensions.
C
C
C     *****     Number of the centres
C
C     given the label of the centre to which the cgto belongs, assign
C     the centre a number. also collect the different labels to a
C     separate list
C
      NNUC = 1
      NAME(1) = IAMN(1)
      NUMBER(1) = 1
      NTYP = IPR(1)
      DO 20 IA = 2,NBAST
         ITYP = IPR(IA)
         NTYP = MAX(NTYP, ITYP)
         DO 10 NUC = 1,NNUC
            IF (IAMN(IA).EQ.NAME(NUC)) THEN
C
C              This centre is already known to us
C
               NUMBER(IA) = NUC
               GO TO 20
            END IF
   10    CONTINUE
C
C        A new label is encountered. Enter it to the list
C
         NNUC = NNUC + 1
         NAME(NNUC) = IAMN(IA)
         NUMBER(IA) = NNUC
   20 CONTINUE
      IF (NNUC.GT.MXCENT)   GO TO 910
      IF (NTYP.GT.MAXTYP)   GO TO 920
C
C   Initialize the matrices where populations will be cumulated
C
      DO 30 NUCB = 1,NNUC
         DO 30 NUCA = 1,NNUC
            POPNAB(NUCA,NUCB) = D0
   30 CONTINUE
C
      CALL DZERO(POPCGT,NBAST)
C
      IND = NNUC * (NNUC + 1) / 2
      CALL DZERO(POPNET,IND)
C
      CALL DZERO(POPC, NOCCT*NBAST )
C
C     *****     Prepare the raw material
C
C     at this stage
C
C         POPNAB will contain the net and overlap populations between
C                nuclei a and b, summed over all mo's
C         POPABI will be the same as popnabi, but resolved for each
C                mo separately (these are written out on tape)
C         POPC   will be the same as previous, but resolved for each
C                mo separately. one line of length nbast is reserved
C                for each mo in core
C         POPNET will contain net and overlap populations between
C                cgto basis functions, summed over all MO's
C
      IF (IPRMUL .GE. 5) THEN
         LUTEMP = -99
         CALL GPOPEN(LUTEMP,' ','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &               .FALSE.)
         REWIND LUTEMP
      END IF
      DO 80 IP = 1,NOCCT
         DO 60 NUCB = 1,NNUC
            DO 60 NUCA = 1,NNUC
               POPABI(NUCA,NUCB) = D0
   60    CONTINUE
C
         IND = 0
         DO 70 IA = 1,NBAST
            DO 70 IB = 1,IA
               IND = IND + 1
C
               NUCA = NUMBER(IA)
               NUCB = NUMBER(IB)
               IF ( NUCA.LT.NUCB ) THEN
                  NUC  = NUCA
                  NUCA = NUCB
                  NUCB = NUC
               END IF
C
               POPAB = CMOCGT(IA,IP)*CMOCGT(IB,IP)
     *               * OCCUP(IP)*SS(IND)
               IF (IA.NE.IB)   POPAB = D2*POPAB
C
C              Distribute
C
               POPNAB(NUCA,NUCB) = POPNAB(NUCA,NUCB) + POPAB
               POPABI(NUCA,NUCB) = POPABI(NUCA,NUCB) + POPAB
               POPNET (IND)      = POPNET(IND)       + POPAB
C
C              (gross populations follow)
C
               POPAB = DP5 * POPAB
               POPC(IA,IP) = POPC(IA,IP) + POPAB
               POPC(IB,IP) = POPC(IB,IP) + POPAB
   70    CONTINUE
C
         IF (IPRMUL .GE. 5) WRITE (LUTEMP) POPABI
C
   80 CONTINUE
C
C     *****     compute some refined results
C
C     Now   POPTYP will contain gross populations on the cgto's
C                  of type t on centre a, summed over all mo's
C           POPCGT will contain the gross populations on each cgto,
C                  summed over all mo's
C           POPNUC will contain the gross populations on each atom,
C                  summed over all mo's
C
      DO 90 NUC = 1,NNUC
         DO 90 ITYP = 1,MAXTYP
            POPTYP(ITYP,NUC) = D0
   90 CONTINUE
C
      DO 100 IA = 1,NBAST
         NUC = NUMBER(IA)
         ITYP = IPR(IA)
         DO 100 IP = 1,NOCCT
            POPTYP(ITYP,NUC) = POPTYP(ITYP,NUC) + POPC(IA,IP)
            POPCGT(IA)       = POPCGT(IA)       + POPC(IA,IP)
  100 CONTINUE
C
      DO 110 NUCA = 1,NNUC
         POPAB = D0
         DO 105 NUCB = 1,NNUC
            POPAB = POPAB + POPNAB(NUCA,NUCB)
     1                    + POPNAB(NUCB,NUCA)
  105    CONTINUE
         POPNUC(NUCA) = DP5 * POPAB
  110 CONTINUE
C
C     *****     Now just write everything out in neat order
C
C
C     POPULATIONS BETWEEN CENTRES
C
      IF ( IPRMUL.GE.5 ) THEN
         WRITE (LUPRI,1001)
 1001    FORMAT (//6X,'Mulliken net and overlap populations',
     *           //6X,'Net and overlap between centres A and B',
     *           //6X,'Contributions from individual MOs')
         REWIND LUTEMP
         DO 130 IP = 1,NOCCT
            WRITE (LUPRI,1002) IP
 1002       FORMAT (//6X,'Molecular orbital number',I4//)
            READ (LUTEMP) POPABI
C
            IEND = 0
    3          IST  = IEND + 1
C              IEND = MIN(NNUC, IEND + 10)
C901028: use 7 instead of 10 (to fit in 80 columns)
               IEND = MIN(NNUC, IEND + 7)
               WRITE (LUPRI,1009) (NAME(NUC),NUC=IST,IEND)
               WRITE (LUPRI,1009)
               DO 120 NUCA = IST,NNUC
                  IENDB = MIN(NUCA,IEND)
                  WRITE (LUPRI,1012)
     *               NAME(NUCA),(POPABI(NUCA,NUCB),NUCB=IST,IENDB)
  120          CONTINUE
            IF (IEND.LT.NNUC)   GO TO 3
  130    CONTINUE
         CALL GPCLOSE(LUTEMP,'DELETE')
      END IF
C
      if (iprmul .ge. 0) then
      WRITE (LUPRI,1004)
      IF (LUW4 .NE. LUPRI) WRITE (LUW4,1004)
 1004 FORMAT(//6X,'Mulliken total overlap populations',
     *        /6X,'Net and overlap populations between centres A and B',
     *        /6X,'Sum over all MOs P',//)
C
      IEND = 0
    5    IST  = IEND + 1
         IEND = MIN(NNUC, IEND + 10)
C901028: use 7 instead of 10 (to fit in 80 columns)
C        IEND = MIN(NNUC, IEND + 7)
         WRITE (LUPRI,1009) (NAME(NUC),NUC=IST,IEND)
         WRITE (LUPRI,1009)
         IF (LUW4 .NE. LUPRI) THEN
            WRITE (LUW4,1009) (NAME(NUC),NUC=IST,IEND)
            WRITE (LUW4,1009)
         END IF
         DO 140 NUCA = IST,NNUC
            IENDB = MIN(NUCA,IEND)
            WRITE (LUPRI,1012)
     *         NAME(NUCA),(POPNAB(NUCA,NUCB),NUCB=IST,IENDB)
            IF (LUW4 .NE. LUPRI) THEN
               WRITE (LUW4,1012)
     *            NAME(NUCA),(POPNAB(NUCA,NUCB),NUCB=IST,IENDB)
            END IF
  140    CONTINUE
      IF (IEND.LT.NNUC)   GO TO 5
      end if
C
C     Gross populations of cgto's
C
      IF (IPRMUL.GE.2) THEN
         WRITE (LUPRI,1005)
         IF (LUW4 .NE. LUPRI) WRITE (LUW4,1005)
 1005    FORMAT (//6X,'Partial orbital gross populations',
     *            /6X,'Gross populations of CGTO basis functions A',
     *            /6X,'Contribution from MO no. P'//)
         DO 160 IP = 1,NOCCT
            WRITE (LUPRI,1002) IP
            WRITE (LUPRI,1006) (IA,POPC(IA,IP),IA=1,NBAST)
            IF (LUW4 .NE. LUPRI) THEN
               WRITE (LUW4,1002) IP
               WRITE (LUW4,1006) (IA,POPC(IA,IP),IA=1,NBAST)
            END IF
  160    CONTINUE
 1006    FORMAT (6(1X,I3,F8.4))
      END IF
C
      if (iprmul .gt. 0) then
      WRITE (LUPRI,1007)
      WRITE (LUPRI,1006) (IA,POPCGT(IA),IA = 1,NBAST)
      IF (LUW4 .NE. LUPRI) THEN
         WRITE (LUW4,1007)
         WRITE (LUW4,1006) (IA,POPCGT(IA),IA = 1,NBAST)
      END IF
      end if
 1007 FORMAT(//6X,
     *   'Mulliken total gross populations of CGTO basis functions A',
     *        /6X,'Sum over all MOs'//)
C
C     GROSS POPULATIONS OF ORBITAL TYPES
C
      IF (IPRMUL.GE.2) THEN
         WRITE (LUPRI,1008)
         IF (LUW4 .NE. LUPRI) WRITE (LUW4,1008)
 1008    FORMAT(//6X,'Partial gross population of CGTO basis',
     *               ' functions of type T',/6X,'for each centre',
     *           /6X,'Contribution from MO no. P'//)
         DO 200 IP = 1,NOCCT
            WRITE (LUPRI,1002) IP
            IF (LUW4 .NE. LUPRI) WRITE (LUW4,1002) IP
C
            DO 170 NUC = 1,NNUC
               DO 170 ITYP = 1,NTYP
                  POPTY(ITYP,NUC) = D0
  170       CONTINUE
C
            DO 180 IA = 1,NBAST
               ITYP = IPR(IA)
               NUC = NUMBER(IA)
               POPTY(ITYP,NUC) = POPTY(ITYP,NUC) + POPC(IA,IP)
  180       CONTINUE
C
            IEND = 0
    7       IST = IEND + 1
            IEND = IEND + 10
C901028: use 7 instead of 10 (to fit in 80 columns)
C           IEND = IEND + 7
            IF (IEND.GT.NNUC) IEND = NNUC
            WRITE (LUPRI,1009) (NAME(NUC),NUC = IST,IEND)
            WRITE (LUPRI,1011)
            IF (LUW4 .NE. LUPRI) THEN
               WRITE (LUW4,1009) (NAME(NUC),NUC = IST,IEND)
               WRITE (LUW4,1011)
            END IF
 1009       FORMAT (' Centre',4X,A4,(9(6X,A4)))
 1011       FORMAT ('  Type'/)
            DO 190 ITYP = 1,NTYP
               WRITE (LUPRI,1012)
     &         GTOTYP(ITYP),(POPTY(ITYP,NUC), NUC=IST,IEND)
               IF (LUW4 .NE. LUPRI) WRITE (LUW4,1012)
     &         GTOTYP(ITYP),(POPTY(ITYP,NUC), NUC=IST,IEND)
  190       CONTINUE
 1012       FORMAT (1X,A4,2X,10F9.4)
            IF (IEND.LT.NNUC)   GO TO 7
  200    CONTINUE
      END IF
C
      if (iprmul .gt. 0) then
      WRITE (LUPRI,1013)
      IF (LUW4 .NE. LUPRI) WRITE (LUW4,1013)
 1013 FORMAT (//6X,'Mulliken total gross population of CGTO basis ',
     *             'functions of type T',/6X,'for each centre'
     *         /6X,'Sum over all MOs'//)
      IEND = 0
    9 IST  = IEND + 1
      IEND = MIN(NNUC, IEND + 10)
C901028: use 7 instead of 10 (to fit in 80 columns)
C     IEND = MIN(NNUC, IEND + 7)
      WRITE (LUPRI,1009) (NAME(NUC),NUC = IST,IEND)
      WRITE (LUPRI,1011)
      IF (LUW4 .NE. LUPRI) THEN
         WRITE (LUW4,1009) (NAME(NUC),NUC = IST,IEND)
         WRITE (LUW4,1011)
      END IF
      LQ     = 0
      DO NUC = IST,IEND
         POPTSM(NUC) = D0
      END DO
      DO 210 ITYP = 1,NTYP
         WRITE (LUPRI,1012) GTOTYP(ITYP),(POPTYP(ITYP,NUC),NUC=IST,IEND)
         IF (LUW4 .NE. LUPRI)
     &     WRITE(LUW4,1012) GTOTYP(ITYP),(POPTYP(ITYP,NUC),NUC=IST,IEND)
C
         IF (LQ.GT.0) THEN
            DO 211 NUC = IST,IEND
               POPTSM(NUC) = POPTSM(NUC) + POPTYP(ITYP,NUC)
  211       CONTINUE
         END IF
      IF (GTOTYP(ITYP+1)(1:1) .EQ. GTOTYP(ITYP)(1:1)) GO TO 210
C     for cartesian GTOTYP(*)(1:1) will be 's', 'p', etc.
C     for spherical GTOTYP(*)(1:1) will be '1', '2', etc.
C     thus: when char 1 changes so does l quantum number and
C     using GTOTYP to decide works both for cartesian and spherical
         IF (LQ.GT.0) THEN
            WRITE (LUPRI,1012) SPDCAR(LQ),(POPTSM(NUC), NUC=IST,IEND)
            IF (LUW4 .NE. LUPRI) THEN
               WRITE (LUW4,1012) SPDCAR(LQ),(POPTSM(NUC), NUC=IST,IEND)
            END IF
         END IF
         WRITE (LUPRI,1099)
         IF (LUW4 .NE. LUPRI) WRITE (LUW4,1099)
         DO NUC = IST,IEND
            POPTSM(NUC) = D0
         END DO
         LQ = LQ + 1
  210 CONTINUE
 1099 FORMAT (1X)
      IF (IEND.LT.NNUC)   GO TO 9
      end if
C
C     POPULATION AT EACH CENTRE
C
      IF (IPRMUL.GE.1) THEN
        WRITE (LUPRI,1014)
        IF (LUW4 .NE. LUPRI) WRITE (LUW4,1014)
 1014 FORMAT (//6X,'Partial gross population of centre A',
     *         /6X,'Contribution from MO no. P'//)
C
        DO 240 IP = 1,NOCCT
          DO 220 NUC = 1,NNUC
            POPABI(NUC,1) = D0
  220     CONTINUE
          DO 230 IA = 1,NBAST
            NUC = NUMBER(IA)
            POPABI(NUC,1) = POPABI(NUC,1) + POPC(IA,IP)
  230     CONTINUE
C
          WRITE (LUPRI,1002) IP
          WRITE (LUPRI,1015) (NAME(NUC),POPABI(NUC,1),NUC=1,NNUC)
c1015     FORMAT (4(3X,A4,F10.4))
 1015     FORMAT (5(1X,A4,F8.2,' ;'))
          IF (LUW4 .NE. LUPRI) THEN
            WRITE (LUW4,1002) IP
            WRITE (LUW4,1015) (NAME(NUC),POPABI(NUC,1),NUC=1,NNUC)
          END IF
  240 CONTINUE
      END IF
C
      WRITE (LUPRI,1016)
      WRITE (LUPRI,1015) (NAME(NUC),POPNUC(NUC),NUC=1,NNUC)
      WRITE (LUPRI,'()')
      IF (LUW4 .NE. LUPRI) THEN
         WRITE (LUW4,1016)
         WRITE (LUW4,1015) (NAME(NUC),POPNUC(NUC),NUC=1,NNUC)
      END IF
 1016 FORMAT (/' Mulliken total gross populations')
c1016 FORMAT (//6X,'Mulliken total gross population of centre A',
c    *         /6X,'Sum over all MOs'/)
C
C
C
C     NET and overlap populations between cgto basis functions - sum
C
      IF (IPRMUL.GE.4) THEN
         WRITE (LUW4,1017)
 1017    FORMAT(//6X,'Total net and overlap populations between CGTO',
     *          ' basis functions'/10X,'Sum over all MOs'//)
C
         IEND = 0
         INDEND = 0
   12    IST = IEND + 1
         IEND = MIN(NBAST, IEND + 10)
C901028: use 7 instead of 10 (to fit in 80 columns)
C        IEND = MIN(NBAST, IEND + 7)
         WRITE (LUW4,1018) (I,I=IST,IEND)
 1018    FORMAT (5X,'CGTO',10I9)
         WRITE (LUW4,1018)
         DO 250 IA=IST,NBAST
            INDSTA = INDEND + IST
            INDE = INDEND + MIN (IA,IEND)
            INDEND = INDEND + IA
            WRITE (LUW4,1019) IA,(POPNET(IND),IND=INDSTA,INDE)
 1019       FORMAT (I7,5X,10F9.4)
  250    CONTINUE
         IF (IEND.LT.NBAST)   GO TO 12
      END IF
C
      GO TO 9999
C
C     Print error messages
C
  910 WRITE (LUW4,9001) NNUC,MXCENT
 9001 FORMAT (/10X,'YOU REALLY MADE IT. THE NUMBER OF NUCLEI,',I5,','
     *        /10X,'EXCEEDS THE PRESENT DIMENSION,',I5/)
      IF (NTYP.LE.MAXTYP)   GO TO 9999
C
  920 WRITE (LUW4,9002) NTYP,MAXTYP
 9002 FORMAT(/10X,'YOU REALLY MADE IT. THE HIGHEST ORBITAL TYPE,',I4,','
     *       /10X,'EXCEEDS THE PRESENT DIMENSION,',I3/)
C     GO TO 9999
C
 9999 CONTINUE
      RETURN
      END
C  /* Deck dipmom */
      SUBROUTINE DIPMOM (MATDIM,MATCMO,CMO,OCCUP,DIPMAT)
C
C PURPOSE:
C       Compute the electronic and total dipole moments
C       Based on release 82 09 01
C
#include <implicit.h>
C
      DIMENSION CMO(MATCMO),DIPMAT(MATDIM,3),OCCUP(NOCCT)
C
C
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0)
C
C   INFORB ? NBAST,
C
#include <mxcent.h>
#include <priunit.h>
#include "infpop.h"
#include <inforb.h>
#include <infpri.h>
C
      LOGICAL FOUND
C
      WRITE (LUPRI,1001)
      IF (LUW4 .NE. LUPRI) WRITE (LUW4,1001)
 1001 FORMAT (///,22X,'Dipole moment analysis')
C
C     *****     READ IN THE DIPOLE MOMENT INTEGRALS
C
      FOUND = .FALSE.
      CALL RDONEL('DIPOLMOM',FOUND,DIPMAT,MATDIM)
      IF ( .NOT. FOUND) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(//A/A/)')
     *     ' *** INFO: Label "DIPOLMOM" not found on unit LUONEL,',
     *     '     basic information for DIPMOM is missing.'
         RETURN
      END IF
C
C
C     *****     Compute the electronic contribution to dipole moment
C
      WRITE (LUPRI,1002)
      IF (LUW4 .NE. LUPRI) WRITE (LUW4,1002)
 1002 FORMAT(//25X,'Orbital dipole moments'//28X,'Vector component'
     *       //10X,'Orbital',12X,'X',11X,'Y',11X,'Z',5X,'Occ.number'/)
      IP = 0
      IPSTRT = 0
      MATBLK = 0
C
      DO 50 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         NOCCI = NISH(ISYM) + NASH(ISYM)
C
         DO 40 ICC = 1,NOCCI
            IP = IP + 1
            IST = MATBLK
            DX = D0
            DY = D0
            DZ = D0
C
            DO 30 J=1,NBASI
               JORB = IPSTRT + J
C
               DO 20 K = 1,(J-1)
                  DENS = D2*CMO(JORB)*CMO(IPSTRT+K)*OCCUP(IP)
                  DX = DX + DENS*DIPMAT(IST+K,1)
                  DY = DY + DENS*DIPMAT(IST+K,2)
                  DZ = DZ + DENS*DIPMAT(IST+K,3)
   20          CONTINUE
C
C              Case j = k, divide by two
C              because no (j,k) permutation exist
C
               DENS = CMO(JORB)*CMO(JORB)*OCCUP(IP)
               DX = DX + DENS*DIPMAT(IST+J,1)
               DY = DY + DENS*DIPMAT(IST+J,2)
               DZ = DZ + DENS*DIPMAT(IST+J,3)
               IST= IST + J
   30       CONTINUE
C
            IPSTRT = IPSTRT + NBASI
C
            IF (IFXYZ(1).EQ.1)   DX = D0
            IF (IFXYZ(2).EQ.1)   DY = D0
            IF (IFXYZ(3).EQ.1)   DZ = D0
            DIPOL(1) = DIPOL(1) - DX
            DIPOL(2) = DIPOL(2) - DY
            DIPOL(3) = DIPOL(3) - DZ
            IF (IPRMUL .GE. 3) THEN
               WRITE (LUPRI,1003) IP,DX,DY,DZ,OCCUP(IP)
               IF (LUW4 .NE. LUPRI)
     1         WRITE (LUW4,1003) IP,DX,DY,DZ,OCCUP(IP)
 1003          FORMAT(10X,I4,7X,4F12.6)
            END IF
C
   40    CONTINUE
C
         MATBLK = MATBLK + NBASI*(NBASI + 1)/2
   50 CONTINUE
C
C     To obtain exactly zero results in cases where symmetry
C     dictates that the dipole moment must vanish, explicitly
C     set the components to zero
C
      IF (IFXYZ(1).EQ.1)   DIPOL(1) = D0
      IF (IFXYZ(2).EQ.1)   DIPOL(2) = D0
      IF (IFXYZ(3).EQ.1)   DIPOL(3) = D0
C
C     *****     Print out the total dipole moments
C
C
      WRITE(LUPRI,1005) (DIPOL(I),I=1,3)
      IF (LUW4 .NE. LUPRI) WRITE(LUW4,1005) (DIPOL(I),I=1,3)
 1005 FORMAT(/11X,'TOTAL',5X,3F12.6)
C
C     MAGNITUDE
C
      DIPITO = DIPOL(1)**2 + DIPOL(2)**2 + DIPOL(3)**2
      DIPITO = SQRT(DIPITO)
      DIPITD = DIPITO * 2.5415D0
      WRITE (LUPRI,1006) DIPITO,DIPITD
      IF (LUW4 .NE. LUPRI) WRITE (LUW4,1006) DIPITO,DIPITD
 1006 FORMAT (/6X,' TOTAL DIPOLE MOMENT (A.U.) = ',F10.6,
     *       //6X,' TOTAL DIPOLE MOMENT (D)    = ',F10.6)
C
C
      RETURN
      END
C  /* Deck quadrp */
      SUBROUTINE QUADRP (MATDIM,MATCMO,CMO,OCCUP,QUMAT)
C
C  Purpose: Compute the quadrupole moments
C
C based on release 82 09 01
C 900709-hjaaj: corrected error in diagonalization of quadrupole tensor
C    relative to mass center (off-diagonal elements were from tensor
C    relative to origo).  Print <r2>.
C
#include <implicit.h>
C
      DIMENSION CMO(MATCMO),QUMAT(MATDIM,6),OCCUP(NOCCT)
C
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0)
C
C   INFORB ? NBAST,
C
#include <mxcent.h>
#include <priunit.h>
#include <inforb.h>
#include "infpop.h"
#include <infpri.h>
C
      DIMENSION QUPOL(6),A(6),EIVR(9),WO(3),IWO(3),ELQPOL(6),
     1          QOPOL(6)
C
      LOGICAL FOUND
C
      WRITE (LUPRI,1001)
      IF (LUW4 .NE. LUPRI) WRITE (LUW4,1001)
 1001 FORMAT (///,21X,' QUADRUPOLE MOMENT ANALYSIS')
C
C     *****     Read in the quadrupole moment integrals
C
      FOUND = .FALSE.
      CALL RDONEL('QUADRP  ',FOUND,QUMAT,MATDIM)
      IF ( .NOT. FOUND) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(//A/A/)')
     *     ' *** INFO: Label "QUADRP" not found on unit LUONEL,',
     *     '     basic information for QUADRP  is missing.'
         RETURN
      END IF
C
C     *****     Compute the electronic contribution to quadrupole moment
C
      IF (IPRMUL .GE. 3) THEN
         WRITE(LUPRI,1002)
         IF (LUW4 .NE. LUPRI) WRITE(LUW4,1002)
      ELSE
         WRITE(LUPRI,2002)
         IF (LUW4 .NE. LUPRI) WRITE(LUW4,2002)
      END IF
 1002 FORMAT(//23X,'Occupation weighted orbital quadrupole moments',
     *        /23X,'(except <r2> which shows spatial extent of orbital)'
     *       //' Orbital',5X,'XX',7X,'YY',7X,'ZZ',7X,'XY',7X,'XZ',
     *         7X,'YZ',7X,'RR     <r2>',/)
 2002 FORMAT(//23X,'Occupation weighted orbital quadrupole moments',
     *       //' Orbital',5X,'XX',7X,'YY',7X,'ZZ',7X,'XY',7X,'XZ',
     *         7X,'YZ',7X,'RR',/)
C
      IP = 0
      IPSTRT = 0
      MATBLK = 0
      DO 20 ICOMP = 1,6
         ELQPOL(ICOMP) = D0
   20 CONTINUE
C
      DO 90 ISYM = 1,NSYM
         NBASI = NBAS(ISYM)
         NOCCI = NISH(ISYM) + NASH(ISYM)
C
         DO 80 ICC = 1,NOCCI
            IP = IP + 1
            DO 30 ICOMP = 1,6
               QUPOL(ICOMP) = D0
   30       CONTINUE
            IST = MATBLK
C
            DO 60 J=1,NBASI
               JORB = IPSTRT + J
C
               DO 40 K = 1,(J-1)
                  DENS = D2*CMO(JORB)*CMO(IPSTRT+K)
                  DO 40 ICOMP = 1,6
                     QUPOL(ICOMP) = QUPOL(ICOMP)
     1                            + DENS*QUMAT(IST+K,ICOMP)
   40          CONTINUE
C
C              CASE J = K, DIVIDE BY TWO BECAUSE
C              NO (J,K) PERMUTATION EXIST
C
               DENS = CMO(JORB)*CMO(JORB)
               DO 50 ICOMP = 1,6
                  QUPOL(ICOMP) = QUPOL(ICOMP)
     1                         + DENS*QUMAT(IST+J,ICOMP)
   50          CONTINUE
               IST= IST + J
   60       CONTINUE
C
            IPSTRT = IPSTRT + NBASI
C           R2 = < r**2 > for this orbital
            R2 = QUPOL(1) + QUPOL(2) + QUPOL(3)
C
            DO 70 ICOMP = 1,6
               QUPOL(ICOMP)  = QUPOL(ICOMP)*OCCUP(IP)
               ELQPOL(ICOMP) = ELQPOL(ICOMP) - QUPOL(ICOMP)
   70       CONTINUE
            RR = QUPOL(1) + QUPOL(2) + QUPOL(3)
            IF (JFXYZ(1).EQ.1)   QUPOL(4) = D0
            IF (JFXYZ(2).EQ.1)   QUPOL(5) = D0
            IF (JFXYZ(3).EQ.1)   QUPOL(6) = D0
            IF (IPRMUL .GE. 3) THEN
               WRITE (LUPRI,1003) IP,(QUPOL(ICOMP),ICOMP=1,6),RR,R2
               IF (LUW4 .NE. LUPRI)
     1         WRITE (LUW4,1003) IP,(QUPOL(ICOMP),ICOMP=1,6),RR,R2
 1003          FORMAT(1X,I4,3X,7F9.4,F7.2)
            END IF
C
   80    CONTINUE
C
         MATBLK = MATBLK + NBASI*(NBASI + 1)/2
   90 CONTINUE
C
C     *****     Print out the total electronic quadrupole moments
C
      SUM = ELQPOL(1) + ELQPOL(2) + ELQPOL(3)
      IF (JFXYZ(1).EQ.1) ELQPOL(4) = D0
      IF (JFXYZ(2).EQ.1) ELQPOL(5) = D0
      IF (JFXYZ(3).EQ.1) ELQPOL(6) = D0
      WRITE(LUPRI,1004) ELQPOL,SUM
      IF (LUW4 .NE. LUPRI) WRITE(LUW4,1004) ELQPOL,SUM
 1004 FORMAT(/' SUMMED ',7F9.4)
C
C     *****     Print out the nuclear quadrupole moment
C
      WRITE (LUPRI,1005) QPOL
      IF (LUW4 .NE. LUPRI) WRITE (LUW4,1005) QPOL
 1005 FORMAT(/' NUCLEAR',6F9.4)
C
C     *****     Print out the total quadrupole moment
C
      DO 100 I=1,6
100      QPOL(I) = QPOL(I) + ELQPOL(I)
      SUM = QPOL(1) + QPOL(2) + QPOL(3)
      WRITE (LUPRI,1006) QPOL,SUM
      IF (LUW4 .NE. LUPRI) WRITE (LUW4,1006) QPOL,SUM
 1006 FORMAT(/'  TOTAL ',7F9.4)
C
C     *****     Quadrupole moments relative to mass centre
C
      QOPOL(1)=QPOL(1)-2*QQFAC(1)*DIPOL(1)
      QOPOL(2)=QPOL(2)-2*QQFAC(2)*DIPOL(2)
      QOPOL(3)=QPOL(3)-2*QQFAC(3)*DIPOL(3)
      QOPOL(4)=QPOL(4)-  QQFAC(1)*DIPOL(2)-QQFAC(2)*DIPOL(1)
      QOPOL(5)=QPOL(5)-  QQFAC(1)*DIPOL(3)-QQFAC(3)*DIPOL(1)
      QOPOL(6)=QPOL(6)-  QQFAC(2)*DIPOL(3)-QQFAC(3)*DIPOL(2)
      SUMR    =QOPOL(1) + QOPOL(2) + QOPOL(3)
      WRITE(LUPRI,1007) QOPOL,SUMR
      IF (LUW4 .NE. LUPRI) WRITE(LUW4,1007) QOPOL,SUMR
 1007 FORMAT(//16X,'QUADRUPOLE MOMENTS RELATIVE MASS CENTRE',
     *       //8X,7F9.4)
C
C     *****     Diagonalize quadrupole moment tensor
C
      ITURN=0
      WRITE(LUPRI,1008)
      IF (LUW4 .NE. LUPRI) WRITE(LUW4,1008)
 1008 FORMAT(//,6X,
     *   'DIAGONALIZATION OF QUADRUPOLE TENSOR RELATIVE ORIGO')
      A(1)=0.5D0*(3.0D0*QPOL(1)-SUM)
      A(3)=0.5D0*(3.0D0*QPOL(2)-SUM)
      A(6)=0.5D0*(3.0D0*QPOL(3)-SUM)
      A(2)=QPOL(4)*1.5D0
      A(4)=QPOL(5)*1.5D0
      A(5)=QPOL(6)*1.5D0
C
  145 CALL DUNIT(EIVR,3)
      CALL JACO(A,EIVR,3,3,3,WO,IWO)
      A(2)=A(3)
      A(3)=A(6)
      WRITE(LUPRI,140)
      IF (LUW4 .NE. LUPRI) WRITE(LUW4,140)
  140 FORMAT(//,18X,'EIGENVALUES',16X,'EIGENVECTORS')
      IEND=0
      DO 141 I=1,3
         IST =IEND+1
         IEND=IEND+3
         WRITE(LUPRI,142)A(I),(EIVR(J),J=IST,IEND)
         IF (LUW4 .NE. LUPRI) WRITE(LUW4,142)A(I),(EIVR(J),J=IST,IEND)
  142    FORMAT(/,20X,F10.6,10X,3F10.6)
  141 CONTINUE
      IF (ITURN.NE.1) THEN
         ITURN = 1
C
         WRITE(LUPRI,1009)
         IF (LUW4 .NE. LUPRI) WRITE(LUW4,1009)
 1009    FORMAT(//'   Diagonalization of quadrupole tensor relative',
     *          ' to mass centre')
         A(1)=0.5D0*(3.0D0*QOPOL(1)-SUM)
         A(3)=0.5D0*(3.0D0*QOPOL(2)-SUM)
         A(6)=0.5D0*(3.0D0*QOPOL(3)-SUM)
         A(2)=QOPOL(4)*1.5D0
         A(4)=QOPOL(5)*1.5D0
         A(5)=QOPOL(6)*1.5D0
         GO TO 145
C     ^-----------
      END IF
C
      RETURN
C -- end of sirpop.F --
      END
