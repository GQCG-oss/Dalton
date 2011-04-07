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
#ifdef NOTES
!===========================================================================
!>>> lontr1.u >>> London 1-index transformation of regular integrals
!941207-hjaaj:
!c: removed SYMTD (=NODIFC) everywhere, not needed!
!   cut CPU time in half in TR1LH1 for SYMTD false
!LFKDEN: corrected code to correspond to formulas, annotated with
! explicit formulas, streamlined, and removed SYMTD extra code
! (old versions had extra factors and signs in LFKDIR, now the result
! should be the same if IFCTYP = 2 or the general IFCTYP = 3.
!9410-hjaaj:
!LFKDIR: removed from this module (now FCKDS2 in SIRFCK module)
!c: s/NODIFC/SYMTD/; SYMTD in parameter list (used to be NODIFC from
!  /ABAINF/, i.e. /ABAINF/ not used any more)
!ABATR1: UDV added to parameter list
!ABARHS: allocate and calculate UDV for ABATR1
!LONTR1: integrated into ABATR1
!930412-kr: New changes due to the new way of calculating right-hand sides.
!931116-hjaaj:
!    LFKDEN: corrected DX[CV]AO(NORBT,NORBT) to DX[CV]AO(NBAST,NBAST)
!    LFK1TR: corrected NBAST to NORBT in 2 AMPAB calls (dim FXC & FXV)
!930929-hjaaj: "cat lontr1.u tr1lfk.u >> siraba.u; rm lontr1.u tr1lfk.u"
!930929-hjaaj: extracted as siraba.u from sirn04.su; replaced abarhs,
!    abaogr, abatr1 with E92 revision for London orbitals (london_sir.u)
!900409-hjaaj: abarhs rev: skip cisigo if rhf (and similar checks).
!891214-hjaaj: this module contains special files for abacus
!===========================================================================
#endif
C  /* Deck abarhs */
      SUBROUTINE ABARHS(NGD,DV,PV,FC,FV,CREF,CMO,
     &                  TD,FDC,FDV,FDQ,H2DAC,LH2DAC,GRDICT,GRDACT,
     &                  GD,XNDXCI,LONDON,DFTFD,WRK,KFRSAV,LFRSAV)
C
C Dec 1989 Hans Joergen Aa. Jensen
C Apr 1992 HJAaJ + Rika Kobayashi (London orbitals)
C
C Purpose: Calculate Abacus GD vectors for right hand side (RHS)
C of response equations.
C
C Input:
C   LH2DAC   dimension of H2DAC
C   LONDON   true if London orbitals
C   etc.
C
C Output:
C   FDC    = FDC(input)   + FSC
C   FDV    = FDV(input)   + FSV
C   FDQ    = FDQ(input)   + FSQ
C   H2DAC  = H2DAC(input) + H2SAC
C   GRDICT = inactiv contribution to gradient = inactive energy
C   GRDACT = active contribution to gradient from CI part of GD
C   GD     = right hand side for ABACUS response calculation
C
C
C Symmetry and number of variational parameters are defined through
C INFLIN.
C
#include <implicit.h>
      DIMENSION DV(*), PV(*), FC(*), FV(*), CREF(*), CMO(*)
      DIMENSION FDC(N2ORBX,*),FDV(N2ORBX,*),FDQ(*)
      DIMENSION H2DAC(LH2DAC,*), DFTFD(N2ORBX,*)
      DIMENSION TD(*), GRDICT(NGD),  GRDACT(NGD),  GD(NVARPT,NGD)
      DIMENSION XNDXCI(*),    WRK(*)
      LOGICAL   LONDON
C
      PARAMETER ( DM1 = -1.0D0, D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0,
     *            D2 = 2.0D0)
C
C  Local:
C
      LOGICAL NOH2, IH8SM
      PARAMETER (NOH2 = .FALSE.)
C
C Used from common blocks:
C   INFORB : NNASHX,N2ASHX,N2ORBX,NASHT
C   INFLIN : NVARPT,LSYMRF,LSYMPT,NCONRF,NCONST,NWOPPT
C CBGETDIS : DISTYP,IADINT
C
#include <maxorb.h>
#include <inforb.h>
#include <inflin.h>
#include <cbgetdis.h>
#include <dftcom.h>
#include <priunit.h>
C
C (0) Initialize
      CALL QENTER('ABARHS')
      CALL DZERO(GD,NGD*NVARPT)
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C     Unpack DV into UDV ( SP to SI format)
      CALL MEMGET('REAL',KUDV,N2ASHX,WRK,KFREE,LFREE)
      IF (NASHT .GT. 0) CALL DSPTSI(NASHT,DV,WRK(KUDV))
C
C (1) Add TD one-index transformation to derivative matrices
C     to obtain the total derivatives.
C
      KFRTR1 = KFREE
      LFRTR1 = LFREE
      CALL ABATR1(NGD,TD,CMO,WRK(KUDV),PV,FC,FV,FDC,FDV,FDQ,
     &            H2DAC,LH2DAC,.TRUE.,LONDON,WRK,KFRTR1,LFRTR1)
C     CALL ABATR1(NDERIV,TR1MAT,CMO,UDV,PV,FC,FV,FXC,FXV,FXQ,
C    &            H2XAC,LH2XAC,ADDTR1,LONDON,WRK,KFRSAV,LFRSAV)
C
C
C (2) Calculate CI part of GD vectors.
C
      IF (NASHT .GT. 1 .AND. NCONST .GT. 0) THEN
         KFRCSI = KFREE
         LFRCSI = LFREE
         CALL MEMGET('REAL',KFDCAC,N2ASHX,WRK,KFRCSI,LFRCSI)
C
C        Set CI control parameters for CISIGD
C
         IADINT = -1
C        ... h2dac in core
         ISPIN1 = 0
         ISPIN2 = 0
C        ... singlet-singlet coupling of 2-electron integrals
         IF (LONDON) THEN
            IH8SM  = .FALSE.
            DISTYP = 9
C           ... DISTYP = 9 for (ij|kl) = H2AC(i,j,k,l) -- no packing !
C           ... HJMAERK 5.oct.92: skal checkes !
C           ... (ij|kl) = (kl|ij) = -(ji|lk) = -(lk|ji)
            GDCFAC = D2
C           950426-hjaaj: old code had DISTYP=10 and GDCFAC=-D2,
C           but that is equivalent to DISTYP=9,GDCFAC=D2
C           (when FDCAC scaled by -1 below).
C           DISTYP=9 is more efficient than DISTYP=10.
         ELSE
            IH8SM  = .TRUE.
            DISTYP = 9
C           ... (ij|kl) = (kl|ij) = (ji|kl) = ... 8-fold symmetry
            GDCFAC = D2
         END IF
C
         DO 200 IGD = 1,NGD
            CALL GETAC1(FDC(1,IGD),WRK(KFDCAC))
            IF (LONDON) CALL DSCAL(N2ASHX,DM1,WRK(KFDCAC),1)
            CALL DZERO(GD(1,IGD),NCONST)
            CALL CISIGD(LSYMRF,LSYMST,NCONRF,NCONST,
     &               CREF,GD(1,IGD),WRK(KFDCAC),H2DAC(1,IGD),
     &               NOH2,IH8SM,XNDXCI,ISPIN1,ISPIN2,
     &               WRK(KFRCSI),LFRCSI)
C           CALL CISIGD(ICSYM,IHCSYM,NCDET,NHCDET, C,HC, FCAC,H2AC,
C    &               NOH2,IH8SM,XNDXCI,ISPIN1,ISPIN2,WRK,LFREE)
            CALL DSCAL(NCONST,GDCFAC,GD(1,IGD),1)
C           ...IF (LSYMPT .eq. 1)
C                 GD(i,IGD) = GD(i,IGD) + D2*EMY(IGD)*CREF(i)
C              but that is removed below when CREF is projected out.
  200    CONTINUE
C
         CALL MEMREL('ABARHS.CISIGD',WRK,KFRSAV,KFREE,KFRCSI,LFRCSI)
      END IF
C
C (3) Calculate orbital part of GD vectors.
C
C     Adding last term of dft-contribution to fock derivative matrix
C
      IF (DFTADD) THEN
         CALL DAXPY(N2ORBX,D1,DFTFD,1,FDC,1)
      END IF
C
      IF (NWOPPT .GT. 0) THEN
         DO 300 IGD = 1,NGD
            JFDQ = 1 + NASHT*NORBT*(IGD-1)
            CALL ABAOGR(WRK(KUDV),FDC(1,IGD),FDV(1,IGD),
     &                  FDQ(JFDQ),LONDON,GD(NCONST+1,IGD))
  300    CONTINUE
      END IF
C     Removing last term of dft-contribution to fock derivative matrix,
C     now the matrix is the generalized fock derivative matrix
      IF (DFTADD) THEN
         CALL DAXPY(N2ORBX,DM1,DFTFD,1,FDC,1)
      END IF
C
C (4) Project out CREF component of GD vectors and
C     calculate GRDICT and GRDACT (the inactive and active
C     contributions to the molecular gradient, calculated from
C     the GD vector).
C
C     GRDICT = inactive energy = sum(k) FDC(kk)
C     If London then GRDICT = 0 because London FDC antisymmetric
C
      IF (.NOT.LONDON .AND. LSYMPT .EQ. 1 .AND. NCONST .GT. 0) THEN
         CALL TR1EMY(NGD,FDC,GRDICT)
      ELSE
         DO IGD = 1,NGD
            GRDICT(IGD) = D0
         END DO
      END IF
      DO 400 IGD = 1,NGD
         IF (LSYMPT .EQ. 1 .AND. NCONST .GT. 0) THEN
            CREFGD = DDOT(NCONST,CREF,1,GD(1,IGD),1)
            CALL DAXPY(NCONST,(-CREFGD),CREF,1,GD(1,IGD),1)
ckr            GRDACT(IGD) = DP5*CREFGD - GRDICT(IGD)
            GRDACT(IGD) = DP5*CREFGD
         ELSE
            GRDACT(IGD) = D0
         END IF
  400 CONTINUE
C
      CALL MEMREL('ABARHS',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('ABARHS')
      RETURN
      END
C  /* Deck abaogr */
      SUBROUTINE ABAOGR(UDV,FDC,FDV,FDQ,IMAGOP,GDORB)
C
C     Jan 90 hjaaj
C     Apr 92 (revised to also handle London orbitals)
C
#include <implicit.h>
      DIMENSION UDV(NASHDI,*)
      DIMENSION FDC(NORBT,*),FDV(NORBT,*),FDQ(NORBT,*), GDORB(*)
      LOGICAL   IMAGOP
C
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0, D4 = 4.0D0)
C
C Used from common blocks:
C   INFVAR : NWOPT,JWOP(2,:)
C   INFORB : NORBT,NNASHX,NNORBT,?
C   INFIND : IROW(:),ISMO(:),ISW(:),ISX(:),IOBTYP(:)
C   INFDIM : NASHDI
C
#include <maxash.h>
#include <maxorb.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>

#include <inflin.h>
C
      IF (IMAGOP) THEN
C     ... Ekl + Elk operator (e.g. B-field)
         FAC =  D2
      ELSE
C     ... Ekl - Elk operator (e.g. geom.deriv.)
         FAC = -D2
      END IF
C
      DO 1300 IG = 1,NWOPPT
         K = JWOP(1,IG)
         L = JWOP(2,IG)
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
C
         IF (ITYPK .EQ. JTINAC) THEN
C **        first index inactive:
            IF (NASHT .GT. 0) THEN
               GDORB(IG) = GDORB(IG) + D4*(FDC(L,K) + FDV(L,K))
            ELSE
               GDORB(IG) = GDORB(IG) + D4*FDC(L,K)
            END IF
         ELSE
C **        first index active:
            ISYMK = ISMO(K)
            NKW   = ISW(K) - NISHT
            TEMP = D0
#include <scalar.h>
            DO 1100 NXW = IASH(ISYMK)+1,IASH(ISYMK)+NASH(ISYMK)
               IX   = ISX(NISHT + NXW)
               TEMP = TEMP + UDV(NKW,NXW)*FDC(L,IX)
 1100       CONTINUE
            GDORB(IG) = GDORB(IG) + D2*(TEMP + FDQ(L,NKW))
         END IF
C
         IF (ITYPL .EQ. JTACT) THEN
C **        second index active:
            ISYML = ISMO(L)
            NLW   = ISW(L) - NISHT
            TEMP = D0
#include <scalar.h>
            DO 1200 NXW = IASH(ISYML)+1,IASH(ISYML)+NASH(ISYML)
               IX   = ISX(NISHT + NXW)
               TEMP = TEMP + UDV(NLW,NXW)*FDC(K,IX)
 1200       CONTINUE
            GDORB(IG) = GDORB(IG) + FAC*(TEMP + FDQ(K,NLW))
         END IF
 1300 CONTINUE
C
      RETURN
      END
C  /* Deck abatr1 */
      SUBROUTINE ABATR1(NDERIV,TR1MAT,CMO,UDV,PV,FC,FV,FXC,FXV,FXQ,
     &                  H2XAC,LH2XAC,ADDTR1,LONDON,
     &                  WRK,KFRSAV,LFRSAV)
C
C Dec 1989 Hans Joergen Aa. Jensen (CALL SIRTR1)
C Revised Apr 1992 for London orbitals (CALL SIRTR1 or LONTR1)
C 23-Apr-1992 Hans Joergen Aa. Jensen and Rika Kobayashi (LONTR1)
C Apr 1993 Kenneth Ruud (always CALL new generalized LONTR1)
C Oct 1994 hjaaj (integrated LONTR1 into ABATR1)
C
C
C Purpose: Calculate one-index transformations needed in ABACUS
C   in FXC,FXV,FXQ,H2XAC.  If ADDTR1 true the result is added
C   to previous contents of these matrices (e.g. for obtaining total
C   derivatives FDC+FSC etc.).
C
C PARAMETERS:
C
C   input : NDERIV,CMO (always needed)
C           LH2XAC = dimension of H2XAC
C           ADDTR1 true : Add FXC,FXV,FXQ,H2XAC to previous
C                         contents
C                  false: Initialize matrices to zero
C   if (NDERIV .gt. 0) then
C     input : TR1MAT,PV,FC,FV
C     output: FXC,FXV,FXQ, H2XAC
C   endif
C   scratch: WRK(KFRSAV:KFRSAV-1+LFRSAV)
C
C
C Symmetry and number of variational parameters are defined through
C INFLIN.
C
C If (LONDON) then perform one-index transformation for
C             London gauge-invariant orbitals
C
#include <implicit.h>
      DIMENSION UDV(*), PV(*), FC(*), FV(*), CMO(*)
      DIMENSION FXC(N2ORBX,*), FXV(N2ORBX,*), FXQ(NASHT*NORBT,*),
     &          TR1MAT(N2ORBX,*), H2XAC(LH2XAC,*)
      DIMENSION WRK(*)
      LOGICAL   ADDTR1, LONDON
C
C *** local constants
C
      PARAMETER ( D0 = 0.0D0, DP5 = 0.5D0 )
C
C Used from common blocks:
C   INFINP : FLAG(*)
C   INFORB : N2ORBX,NASHT,NORBT,...
C   INFVAR : NWOPT,JWOPSY
C   INFDIM : NWOPDI
C   INFLIN : NCONRF,IPRLIN
C   INFTRA : USEDRC
C   INFPRI : P6FLAG(*), ?
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <infdim.h>
#include <inflin.h>
#include <inftra.h>
#include <infpri.h>
C
C
      CALL QENTER('ABATR1')
C
C *** Check input
C
      IF (NDERIV .LE. 0) GO TO 9999
      IF (NWOPT .NE. NWOPDI .AND. NDERIV .GT. 1) THEN
C        SIRIUS/ABACUS program error, correct code.
         WRITE (LUPRI,'(//2A,2I8)') ' >>> FATAL ERROR (ABATR1)',
     *      ' dimension error; NWOPT and NWOPDI =',NWOPT,NWOPDI
         CALL QTRACE(LUPRI)
         CALL QUIT('FATAL ERROR (ABATR1) NWOPDI DIMENSION ERROR')
      END IF
C
C *** Allocate work area for TR1LFK,TR1L2*
C
      KWRK1  = KFRSAV
      LWRK1  = LFRSAV
C
C *** Print one-index transformation matrices, if desired
C
      IF (P6FLAG(19) .OR. IPRLIN .GE. 19) THEN
         DO 210 IDERIV = 1,NDERIV
            WRITE (LUPRI,2110) IDERIV,NDERIV
            CALL OUTPUT(TR1MAT(1,IDERIV),1,NORBT,1,NORBT,NORBT,NORBT,
     &                  1,LUPRI)
  210    CONTINUE
      END IF
 2110 FORMAT (/,' ABATR1: one-index transformation matrix',
     *  ' (no.',I3,' of',I3,')')
C
C     Initialization of output matrices
C
      IF (ADDTR1) THEN
         IF (NASHT .GT. 1) THEN
            CALL DSCAL((NDERIV*LH2XAC),DP5,H2XAC,1)
C           HJMAERK: afh. af om H2DAC er symmetriseret
C           921002-hj: tror jeg er OK (H2DAC antages symmetriseret)
C           ... divide by two because H2XAC and its transposed
C               are added in TR1L2M
         END IF
      ELSE
         CALL DZERO(FXC,NDERIV*N2ORBX)
         CALL DZERO(FXV,NDERIV*N2ORBX)
         CALL DZERO(FXQ,(NDERIV*NASHT*NORBT))
         CALL DZERO(H2XAC,(NDERIV*LH2XAC))
      END IF
C
C
C ***  Calculate contributions from MO two-electron integrals
C
C
      IF (NASHT .GT. 1) THEN
         IF (LONDON) THEN
            LPV = N2ASHX
         ELSE
            LPV = NNASHX
         END IF
         CALL TR1L2M(NDERIV, TR1MAT, PV,LPV, H2XAC,FXQ,
     *               .NOT.USEDRC,LONDON,WRK(KWRK1),LWRK1)
         IF (USEDRC) THEN
            CALL TR1L2D(NDERIV, TR1MAT,PV,LPV,FXQ, LONDON,
     &                  WRK(KWRK1),LWRK1)
         END IF
C
C        CALL TR1L2M(NDERIV,TR1MAT,PV,H2XAC,FXQ, NODRC,LONDON,
C    &               WRK,LFRSAV)
      END IF
C
C
C *** Calculate various transformed Fock matrices from the atomic
C     two-electron integrals and untransformed
C     Fock matrices FC and FV.
C
      CALL TR1LFK(NDERIV,TR1MAT,CMO,UDV,FXC,FXV,
     &            LONDON,WRK(KWRK1),LWRK1)
C
      DO 1200 IDERIV = 1,NDERIV
         IF (P6FLAG(33) .OR. IPRLIN.GE.21 .AND. NASHT.GT.0) THEN
            WRITE (LUPRI,3011) IDERIV,NDERIV
            CALL OUTPUT(FXQ(1,IDERIV),1,NORBT,1,NASHT,NORBT,NASHT,
     &                  1,LUPRI)
         END IF
         IF (P6FLAG(21) .OR. IPRLIN.GE.21) THEN
            WRITE (LUPRI,5020) IDERIV,NDERIV
            CALL OUTPUT(FXC(1,IDERIV),1,NORBT,1,NORBT,
     &                  NORBT,NORBT,1,LUPRI)
            IF (NASHT.GT.0) THEN
               WRITE (LUPRI,5030) IDERIV,NDERIV
               CALL OUTPUT(FXV(1,IDERIV),1,NORBT,1,NORBT,
     &                     NORBT,NORBT,1,LUPRI)
            END IF
         END IF
 1200 CONTINUE
 3011 FORMAT(/' FXQ contribution to one-index transformed Fock matrix ',
     *  '(no.',I3,' of',I3,')')
 5020 FORMAT (/' FXC matrix without FC transf. (no.',I3,' of',I3,')')
 5030 FORMAT (/' FXV matrix without FV transf. (no.',I3,' of',I3,')')
C
C
      CALL TR1LFD(NDERIV,FXC,FXV,TR1MAT,FC,FV,
     &            LONDON,WRK(KWRK1),LWRK1)
C
C ***
C
      IF (KWRK1 .NE. KFRSAV)
     &   CALL MEMREL('ABATR1',WRK,KFRSAV,KFRSAV,KWRK1,LWRK1)
 9999 CALL QEXIT('ABATR1')
      RETURN
      END
C  /* Deck abalin */
      SUBROUTINE ABALIN(NCSIM,NOSIM,BCVECS,BOVECS,
     *                  CMO,CREF,GORB,DV,PV,
     *                  FC,FV,FCAC,H2AC,INDXCI,WRK,KFRSAV,LFRSAV)
C
C Dec 1989 Hans Joergen Aa. Jensen
C
C Purpose: Calculate linear transformations for solving response
C equations in ABACUS.
C
C Symmetry and number of variational parameters are defined through
C INFLIN.
C
#include <implicit.h>
      DIMENSION BCVECS(*), BOVECS(*)
      DIMENSION CREF(*), GORB(*), CMO(*), DV(*), PV(*)
      DIMENSION FC(*), FV(*), FCAC(*), H2AC(*)
      DIMENSION INDXCI(*),WRK(*)
C
#include <thrzer.h>
C
C Used from common blocks:
C  INFOPT : EACTIV
C  INFLIN : LSYMPT,NCONST,NVARPT
C
#include <infopt.h>
#include <inflin.h>
C
#include <sirbkd.h>
C
C (0) Initialize
C
      CALL QENTER('ABALIN')
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C (1) call lintrn, sigma vectors are returned in WRK
C
      KSVECS = KFREE
      CALL LINTRN(NCSIM,NOSIM,BCVECS,BOVECS,
     *            CMO,CREF,EACTIV,GORB,DV,PV,
     *            FC,FV,FCAC,H2AC,INDXCI,WRK,KFREE,LFREE)
C     CALL LINTRN(NCSIM,NOSIM,BCVECS,BOVECS,
C    *            CMO,CREF,EACTIV,GORB,DV,PV,
C    *            FC,FV,FCAC,H2AC,INDXCI,WRK,KFRSAV,LFRSAV)
C
C
C (2) remove CREF component in sigma vectors
C
      IF (LSYMPT .EQ. 1) THEN
         NTSIM = MAX(0,NCSIM) + MAX(0,NOSIM)
         DO 200 ITSIM = 1,NTSIM
            JSVECS = KSVECS + (ITSIM-1)*NVARPT
            CREFSV = DDOT(NCONST,CREF,1,WRK(JSVECS),1)
            IF (ABS(CREFSV) .GT. THRZER) THEN
               CALL DAXPY(NCONST,(-CREFSV),CREF,1,WRK(JSVECS),1)
            END IF
  200    CONTINUE
      END IF
C
      IF (KFREE .NE. KFRSAV)
     &   CALL MEMREL('ABALIN',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('ABALIN')
      RETURN
      END
C  /* Deck tr1lfd */
      SUBROUTINE TR1LFD(NOSIM,FXC,FXV, TMAT,FC,FV,
     &                  LONDON,WRK,LFRSAV)
C
C 23. Apr. 1992 HJAaJ + RK (based on TR1FD)
C
C Purpose:
C
C  To complete the construction of NOSIM one-index transformed,
C  inactive and active Fock matrices (FXC and FXV) by
C  adding the one-index transformed FC and FV matrices to
C  FXC and FXV, resp.
C
#include <implicit.h>
      LOGICAL LONDON
      DIMENSION FXC(N2ORBX,*),FXV(N2ORBX,*)
      DIMENSION TMAT(N2ORBX,*), FC(*),FV(*)
      DIMENSION WRK(LFRSAV)
C
      PARAMETER ( DP5 = 0.5D0 )
C
C Used from common blocks:
C   INFORB : N2ORBX,N2ORBT,NNBASX,NNASHX
C   INFLIN : IPRLIN
C   INFPRI : P6FLAG(*)
C
#include <priunit.h>
#include <inforb.h>
#include <inflin.h>
#include <infpri.h>
C
      CALL QENTER('TR1LFD ')
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C
C ***** Calculate one-index transformation of FC and FV, the inactive
C ***** and active Fock matrices, resp., and add to FXC and FXV.
C ***** This completes FXC and FXV.
C
C
      CALL MEMGET('REAL',KFTMP,NNORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFC  ,N2ORBX,WRK,KFREE,LFREE)
      CALL PKSYM1(WRK(KFTMP),FC,NORB,NSYM,-1)
      CALL DSPTSI(NORBT,WRK(KFTMP),WRK(KFC))
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KFV  ,N2ORBX,WRK,KFREE,LFREE)
         CALL PKSYM1(WRK(KFTMP),FV,NORB,NSYM,-1)
         CALL DSPTSI(NORBT,WRK(KFTMP),WRK(KFV))
      ELSE
         KFV = KFC
C        ... to avoid compiler messages
      END IF
      DO 5900 IOSIM = 1,NOSIM
         CALL DSCAL(N2ORBX,DP5,FXC(1,IOSIM),1)
         CALL TR1LH1(TMAT(1,IOSIM),WRK(KFC),FXC(1,IOSIM),1,
     &               LONDON)
         IF (NASHT .GT. 0) THEN
            CALL DSCAL(N2ORBX,DP5,FXV(1,IOSIM),1)
            CALL TR1LH1(TMAT(1,IOSIM),WRK(KFV),FXV(1,IOSIM),1,
     &                  LONDON)
         END IF
C
         IF (P6FLAG(21) .OR. IPRLIN.GE.21) THEN
            WRITE (LUPRI,5020) IOSIM,NOSIM
            CALL OUTPUT(FXC(1,IOSIM),1,NORBT,1,NORBT,
     &                  NORBT,NORBT,1,LUPRI)
            IF (NASHT.GT.0) THEN
               WRITE (LUPRI,5030) IOSIM,NOSIM
               CALL OUTPUT(FXV(1,IOSIM),1,NORBT,1,NORBT,
     &                     NORBT,NORBT,1,LUPRI)
            END IF
         END IF
C
 5020 FORMAT (/' FXC matrix (no.',I3,' of',I3,')')
 5030 FORMAT (/' FXV matrix (no.',I3,' of',I3,')')
C
 5900 CONTINUE
      CALL MEMREL('TR1LFD.FXC/V',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C *** End of subroutine TR1LFD
C
      CALL QEXIT('TR1LFD ')
      RETURN
      END
C  /* Deck tr1l2m */
      SUBROUTINE TR1L2M(NOSIM,TMAT,PV,LPV,H2XAC,FXQ, NODRC,LONDON,
     &                  WRK,LFRSAV)
C
C Written 23. Apr. 1992 by Hans Jorgen Aa. Jensen and Rika Kobayashi
C Revised 14-Sep-1992 hjaaj
C (based on TR1H2M)
C
C Purpose:
C   Construct one-index transformed 2-electron integrals (*u/vx)
C   needed for ABARHS with London orbitals.
C   To construct FXQ :
C     FXQ(p,t) = SUM(uvx): PV(tu:vx)*(pXuX:vXxX)
C   NOTE: FXQ and H2XAC are added to input
C
C Input:
C   If (LONDON) Full two-electron density matrices PV(N2ASHX,N2ASHX)
C               ( PV(pq;rs) = <e(pq;rs)> )
C   Else        symmtrized PV(NNASHX,NNASHX)
C               ( PV(pq;rs) = .5 * <e(pq;rs)+e(qp;rs)> )
C   Abacus T matrix as a full matrix in TMAT
C
C Output:
C   "FXQ" matrix in FXQ(NORBT,NASHT,NOSIM)
C   H2XAC One-index transformed integrals over active orbitals
C
C Scratch:
C   WRK(KFRSAV:LFRSAV)
C
C
C Externals:
C   DDOT, DZERO, MOLLAB, QTRACE
#include <implicit.h>
      DIMENSION PV(LPV,*),TMAT(N2ORBX,*)
      DIMENSION H2XAC(NASHT,NASHT,NASHT,NASHT,*)
C     HJMAERK: N2ASHX,NNASHX eller NNASHX,NNASHX half-transf. int. ?
C     921002: use N2ASHX,N2ASHX so we may use old GETINT option
      DIMENSION FXQ(NASHT*NORBT,*)
      DIMENSION WRK(LFRSAV)
      LOGICAL   NODRC, LONDON
      PARAMETER ( D2 = 2.0D0 )
C
C Used from common blocks
C   INFINP : FLAG(*)
C   INFORB : NSYM,NNASHX,NORBT,...
C   INFIND : IROW(*),ICH(*),IOBTYP(*),?
C   INFVAR : JWOPSY
C   INFLIN : IPRLIN
C   INFPRI : P6FLAG(*),?
C CBGETDIS : IADH2X
C
#include <maxash.h>
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
#include <inflin.h>
#include <infpri.h>
#include <cbgetdis.h>
C
#include <orbtypdef.h>
C
C     Local variables
C
      DIMENSION NEEDMU(6)

      CALL QENTER('TR1L2M')
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LFRSAV
C
C     Check input
C

      IF (NASHT .LE. 1) GO TO 9999
C     (FXQ zero matrix and H2XAC not used if NASHT .le. 1)
C
      IF (NOSIM .LE. 0) GO TO 9999
C     ... only one-index transformation if orbital vectors
C
C     Set NEEDMU array
C
      DO 10 I = 1,6
         NEEDMU(I) = 0
   10 CONTINUE
      NEEDMU(3) = 1
C     ... active-active distributions always needed
      IF (NODRC) THEN
C        ... all distributions involving active orbitals needed
C            for FXQ
         NEEDMU(2) = 1
         NEEDMU(5) = 1
      END IF
C
C     Allocate work memory
C
      IF (LONDON) THEN
         LPVFCD = 0
      ELSE
         LPVFCD = N2ASHX
      END IF
      CALL MEMGET('REAL',KPVFCD,LPVFCD,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2CD ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KH2XCD,N2ORBX,WRK,KFREE,LFREE)
C
      JDIST = 0
      IF (P6FLAG(36) .OR. IPRLIN .GE. 40) THEN
         WRITE (LUPRI,'(///A)')
     &      ' Test output of Mulliken distributions from TR1L2M.'
      END IF
C
C ****************************************************************
C     Loop over Mulliken distributions allowed in NEEDMU(6)
C
      IDIST = 0
  100 CALL NXTH2M(IC,ID,WRK(KH2CD),NEEDMU,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
      IF (P6FLAG(36) .OR. IPRLIN .GE. 40) THEN
         JDIST = JDIST + 1
         WRITE (LUPRI,'(/A,I5,A,2I5)')
     &      ' Mulliken distribution no.',JDIST,', IC and ID:',IC,ID
         WRITE (LUPRI,'(A,I5)') ' IDIST form NXTH2M :',IDIST
         CALL OUTPUT(WRK(KH2CD),1,NORBT,1,NORBT,NORBT,NORBT,-1,LUPRI)
      END IF
C
C     Find distribution type (only distributions with at least
C     one active index are needed):
C
         ITYPC  = IOBTYP(IC)
         ITYPD  = IOBTYP(ID)
         ITYPCD = IDBTYP(ITYPC,ITYPD)
      IF (ITYPC .NE. JTACT .AND. ITYPD .NE. JTACT) GO TO 100
C     ... only distributions with at least one active index needed
C         (previous statement should never become
C          active because of NEEDMU)
         ISWAP = 0
         IF (ITYPCD .EQ. JTACAC) THEN
            IF (ISW(IC) .LT. ISW(ID)) ISWAP = 1
C           ... we want NCW .ge. NDW for active-active
         ELSE
            IF (IOBTYP(IC) .EQ. JTACT) ISWAP = 1
C           ... we want IC to be non-active index
         END IF
         IF (ISWAP .EQ. 1) THEN
            ISWAP = IC
            IC = ID
            ID = ISWAP
         END IF
         IF (ITYPCD .NE. JTACAC) THEN
            ICXSYM = MULD2H(JWOPSY,ISMO(IC))
            IF (NASH(ICXSYM).EQ.0) GO TO 100
C     ^-------------------------------------
C           if no active orbitals in ICXSYM (new symmetry after
C           one-index transformation of IC), we don't need this
C           distribution
         END IF
C
         ICSYM  = ISMO(IC)
         IDSYM  = ISMO(ID)
         ICDSYM = MULD2H(ICSYM,IDSYM)
         IF (ITYPCD .EQ. JTACAC) THEN
            NCW    = ICH(IC)
            NDW    = ICH(ID)
            NCDW   = (NDW-1)*NASHT + NCW
            NDCW   = (NCW-1)*NASHT + NDW
            IF (.NOT. LONDON) THEN
               CALL SIRPVD(NCW,NDW,WRK(KPVFCD),PV,1)
               IF (NCW.NE.NDW) CALL DSCAL(N2ASHX,D2,WRK(KPVFCD),1)
C              ... C,D and D,C contributions are identical
            END IF
C
            DO 399 IOSIM = 1,NOSIM
               CALL DZERO(WRK(KH2XCD),N2ORBX)
               CALL TR1LH1(TMAT(1,IOSIM),WRK(KH2CD),WRK(KH2XCD),ICDSYM,
     &                     LONDON)
C              CALL TR1LH1(UBOVEC,H1,H1X,IH1SYM,LONDON)
C
C We now have the one-index transformed integrals in H2XCD(*,*):
C (AX BX / C D) = (AX B / C D) - (BX A / C D) = -(BX AX / C D)
C (AX B  / C D) = SUM(r)  TMAT(a,r) * H2CD(r,b)
C We have used TMAT(i,j) = -TMAT(j,i)
C
C              Add contributions to FXQ
C
               IF (LONDON) THEN
                  CALL ADDFQ(NCW,NDW,FXQ(1,IOSIM),
     &                       WRK(KH2XCD),PV(1,NCDW),JWOPSY)
                  IF (NCW .NE. NDW) THEN
                     CALL ADDFQ(NDW,NCW,FXQ(1,IOSIM),
     &                          WRK(KH2XCD),PV(1,NDCW),JWOPSY)
                  END IF
               ELSE
                  CALL ADDFQ(NCW,NDW,FXQ(1,IOSIM),
     &                       WRK(KH2XCD),WRK(KPVFCD),JWOPSY)
               END IF
C              CALL ADDFQ(NXW,NYW,FQ,H2XY,PVXY,JH2SYM)
C
C              add active-active elements to H2XAC
C
               IUVSYM = MULD2H(ICDSYM,JWOPSY)
C              CALL ADH2AC(H2XAC(1,1,NCW,NDW,IOSIM),WRK(KH2XCD),IUVSYM)
               CALL ADDAC1(WRK(KH2XCD),H2XAC(1,1,NCW,NDW,IOSIM),IUVSYM)
               IF (NCW .NE. NDW) THEN
                  CALL ADDAC1(WRK(KH2XCD),H2XAC(1,1,NDW,NCW,IOSIM),
     &                        IUVSYM)
               END IF
C
C              HJMAERK : ADH2AC du'r ikke !
C                        921003: CALL ADDAC1 baseret paa GETAC1
C
  399       CONTINUE
C           end do 399 iosim = 1,nosim
C
         END IF
C        end of if ( (cd) active-active ) then ... end if
C
C        ... If NODRC then
C            now the (p~ v~ | x y) contribution to FXQ
C
C        ... (CD) inactive-active, active-active
C            or secondary-active distribution
C
         IF (NODRC .AND.
     &      (ITYPC.EQ.JTACT .OR. ITYPD.EQ.JTACT)) THEN
            DO 499 IOSIM = 1,NOSIM
               CALL DZERO(WRK(KH2XCD),N2ORBX)
               CALL TR1LH1(TMAT(1,IOSIM),WRK(KH2CD),WRK(KH2XCD),ICDSYM,
     &                     LONDON)
C              CALL TR1LH1(UBOVEC,H1,H1X,IH1SYM,LONDON)
               IF (LONDON) THEN
                  CALL ADFLQM(IC,ID,FXQ(1,IOSIM),WRK(KH2XCD),PV,
     &                        WRK(KFREE),LFREE)
               ELSE
                  CALL ADFXQM(IC,ID,FXQ(1,IOSIM),WRK(KH2XCD),PV,
     &                        WRK(KFREE),LFREE)
               END IF
C
  499       CONTINUE
C           end do 499 iosim = 1,nosim
C
         END IF
C
C        Go to 100 to get next needed Mulliken distribution
C
      GO TO 100
C
C     arrive at 800 when finished with all needed Mulliken distributions
C
  800 CONTINUE
      IF (P6FLAG(36) .OR. IPRLIN .GE. 40) THEN
         WRITE (LUPRI,'(//A/,2(/A,I5))')
     &     ' End of test output of Mulliken distributions from TR1L2M.',
     &     ' Total number of distributions treated  :',JDIST,
     &     ' Total number of distributions (NNORBX) :',NNORBX
      END IF
C
C
C ****************************************************************
C
C  Now H2XAC(ij,kl,iosim) = (it jt / k l),
C  finish H2XAC by symmetrizing
C  (remember (it jt / kt lt) = (it jt / k l) + (kt lt / i j) ).
C
C  If requested, print H2XAC
C
      DO 1800 IOSIM = 1,NOSIM
C     HJMAERK : skal H2XAC symmetriseres for London ????
C     eller skal vi ordne det i CI rutinerne ????
C     921002-hj: vi ordner det her i denne version.
C
         CALL TR1SMT(H2XAC(1,1,1,1,IOSIM),N2ASHX)
         IF (P6FLAG(28) .OR. IPRLIN.GE.28) THEN
            WRITE (LUPRI,1020) IOSIM,NOSIM
            CALL OUTPUT(H2XAC(1,1,1,1,IOSIM),1,N2ASHX,1,N2ASHX,
     *                  N2ASHX,N2ASHX,-1,LUPRI)
         END IF
 1020 FORMAT (/' H2XAC matrix (no.',I3,' of',I3,')')
C
 1800 CONTINUE
C
C *** IADH2X .lt. 0 means that H2XAC is in core.
C
      IADH2X = -1
C
C *** end of subroutine TR1L2M
C
 9999 IF (KFRSAV .NE. KFREE)
     &   CALL MEMREL('TR1L2M',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('TR1L2M')
      RETURN
      END
C  /* Deck tr1smt */
      SUBROUTINE TR1SMT(H2XAC,N2ASHX)
C
C Symmetrize one-index transformation
C
#include <implicit.h>
      DIMENSION H2XAC(N2ASHX,N2ASHX)
      DO 200 KL = 1,N2ASHX
         DO 100 IJ = 1,KL
            X = H2XAC(IJ,KL) + H2XAC(KL,IJ)
            H2XAC(IJ,KL) = X
            H2XAC(KL,IJ) = X
  100    CONTINUE
  200 CONTINUE
C
      RETURN
      END
C  /* Deck adflqm */
      SUBROUTINE ADFLQM(IC,ID,FXQ,H2XCD,PVF,WRK,LWRK)
C
C Sep 92 hjaaj : version for London orbitals
C Based on ADFXQM from Jan 90 hjaaj
C
C Add (p~ v~ | x y) contributions to FXQ from Mulliken integrals.
C
#include <implicit.h>
      DIMENSION FXQ(NORBT,NASHT), H2XCD(NORBT,NORBT),
     &          PVF(N2ASHX,N2ASHX),WRK(LWRK)
C
C
C Used from common blocks:
C  INFORB : NORBT, NASHT, N2ASHX, NNASHX, MULD2H(8,8)
C  INFIND : IROW(:),ISMO(:),NSM(:),ICH(:),IOBTYP(:)
C  INFVAR : JWOPSY
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
#include <infvar.h>
C
      KH2AC1 = 1
      KNEED  = KH2AC1 + N2ASHX
      IF (KNEED .GT. LWRK) CALL ERRWRK('ADFLQM',KNEED,LWRK)
C
      CALL GETAC1(H2XCD,WRK(KH2AC1))
C     CALL GETAC1(HH, HHAC)
C
C We now have the triangular packed one-index transformed integrals
C in H2ACF:
C   H2ACF(uv) = (UX VX / C D)
C where
C   (AX BX / C D) = (AX B / C D) - (A BX / C D)
C   (A  BX / C D) = SUM(r)  TMAT(b,r) * H2CD(a,r)
C
C
C Add PVF(t,d;u,v)*(c d;ux vx) contributions to FXQ(c,t)
C
      IC1 = IC
      ID1 = ID
      IROUND = 0
  420    IROUND = IROUND + 1
         ITYPD  = IOBTYP(ID1)
      IF (ITYPD .EQ. JTACT) THEN
         ICSYM  = ISMO(IC1)
         NDW    = ICH(ID1)
         DO 450 NTW = 1,NASHT
            ITSYM = NSM(NTW)
         IF (MULD2H(ITSYM,ICSYM) .NE. JWOPSY) GO TO 450
            NTDW = (NDW-1)*NASHT + NTW
            FXQ(IC1,NTW) = FXQ(IC1,NTW) +
     *         DDOT(N2ASHX,WRK(KH2AC1),1,PVF(1,NTDW),1)
C       HJMAERK: skal nok vaere N2ASHX og PVF; changed 921003
  450    CONTINUE
      END IF
C
C        If (cd) active-active and c .ne. d we need to
C        consider (dc) distribution also.
C
      IF ( IROUND .EQ. 1 .AND. IC .NE. ID ) THEN
         IC1 = ID
         ID1 = IC
         GO TO 420
      END IF
      RETURN
      END
C  /* Deck tr1l2d */
      SUBROUTINE TR1L2D(A,B,C,D,E,F,G,H)
C
C Apr 92 Hans Joergen Aa. Jensen
C
C Purpose:
C   Calculate contributions from Dirac format integrals
C   of London orbitals
C
#include <implicit.h>
#include <priunit.h>
      CALL QENTER('TR1L2D')
C
      CALL QTRACE (LUPRI)
      CALL QUIT('TR1L2D ERROR: MO integrals Dirac format not '//
     &          'implemented yet')
C
      CALL QEXIT('TR1L2D')
      RETURN
      END
C  /* Deck tr1lh1 */
      SUBROUTINE TR1LH1(TMAT,H1,H1X,IH1SYM,LONDON)
C
C 920423-hjaaj (based on TR1UH1)
C
C Antisymmetric one-index transformation of unpacked H1(norbt,norbt)
C of symmetry IH1SYM to H1X(norbt,norbt) using TMAT(norbt,norbt)
C of symmetry JWOPSY.
C
C Input : H1(a,b) of symmetry IH1SYM
C Output: H1X(a,b) = H1X(a,b)
C                  + H1(a,b) one-index transformed with TMAT
C
C General formula:
C H1X(a,b) = H1X(a,b) + sum(c) [ T(a,c)  H1(c,b) + T*(b,c) H1(a,c) ]
C          = H1X(a,b) + sum(c) [ T*(c,a) H1(c,b) + T(c,b)  H1(a,c) ]
C (First line corresponds to the original code, second line
C is the new modified formulas by K.Ruud, May 1994.  Both forms
C are equivalent for symmetric T matrix, i.e. NODIFC true, but
C only the second form is correct if DIFC reorthonormalization /hjaaj)
C
C Specific formulas for real and imaginary perturbations:
C (real and imag. correspond to not london and london, resp.)
C real       H1X(a,b) =   sum(c) ( T(c,a) H1(c,b) + T(c,b) H1(a,c))
C imag.    i H1X(a,b) = i sum(c) (-T(c,a) H1(c,b) + T(c,b) H1(a,c))
C
C Thus for derivative/London orbitals, i.e. real/imag. :
C  H1X(a,b) = sum(c) H1(a,c) T(c,b)
C  H1X(a,b) = H1X(a,b) +/- H1X(b,a)
C
#include <implicit.h>
      LOGICAL LONDON
      DIMENSION TMAT(NORBT,NORBT), H1(NORBT,NORBT), H1X(NORBT,NORBT)
C
C Used from common blocks:
C  INFORB : NSYM, MULD2H, NORBT, ...
C  INFVAR : JWOPSY
C
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
C
C
C ONE INDEX TRANSFORM FIRST INTEGRAL INDEX
C
C  H1X(a,b) = sum(c) H1(a,c) T(c,b)
C
      DO 1100 IASYM = 1,NSYM
         ICSYM = MULD2H(IASYM,IH1SYM)
         IBSYM = MULD2H(ICSYM,JWOPSY)
         NORBA = NORB(IASYM)
         NORBB = NORB(IBSYM)
         NORBC = NORB(ICSYM)
         IF (NORBA*NORBB*NORBC.NE.0) THEN
            IAST  = IORB(IASYM) + 1
            IBST  = IORB(IBSYM) + 1
            ICST  = IORB(ICSYM) + 1
            CALL DGEMM('N','N',NORBA,NORBB,NORBC,1.D0,
     &                 H1(IAST,ICST),NORBT,
     &                 TMAT(ICST,IBST),NORBT,1.D0,
     &                 H1X(IAST,IBST),NORBT)
         END IF
 1100 CONTINUE
C
C  ADD ONE-INDEX TRANSFORM FROM SECOND INTEGRAL INDEX
C
C           H1X(A,B) = H1X(A,B) +/- H1X(B,A)
C
      IH1XSY = MULD2H(IH1SYM,JWOPSY)
      DO 1200 IASYM = 1,NSYM
         IBSYM = MULD2H(IASYM,IH1XSY)
      IF (IASYM .LT. IBSYM) GO TO 1200
         NORBB = NORB(IBSYM)
         NORBA = NORB(IASYM)
      IF (NORBB.EQ.0 .OR. NORBA.EQ.0) GO TO 1200
         IORBB = IORB(IBSYM)
         IORBA = IORB(IASYM)
         DO 1250 IA = IORBA+1,IORBA+NORBA
            IBEND = MIN(IA,IORBB+NORBB)
            DO 1260 IB = IORBB+1,IBEND
               IF (LONDON) THEN
                  H1X(IA,IB) = H1X(IA,IB) - H1X(IB,IA)
                  H1X(IB,IA) = - H1X(IA,IB)
               ELSE
                  H1X(IA,IB) = H1X(IA,IB) + H1X(IB,IA)
                  H1X(IB,IA) = H1X(IA,IB)
               END IF
 1260       CONTINUE
 1250    CONTINUE
 1200 CONTINUE
C
C     End of TR1LH1
C
      RETURN
      END
C  /* Deck addac1 */
      SUBROUTINE ADDAC1(HH,HHAC,NIJSYM)
C
C  3-Oct-1992 Hans Joergen Aa. Jensen
C
C Purpose: Extract block with active-active orbital indices of symmetry
C          NIJSYM out of full matrix and add to HHAC.
C
#include <implicit.h>
      DIMENSION HH(NORBT,NORBT), HHAC(NASHT,NASHT)
C
C INFORB : NORBT, NASHT, N2ASHX
C INFIND : IOBTYP(*), ICH(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infind.h>
C
      CALL QENTER('ADDAC1')
      IF (NASHT .EQ. 0) GO TO 9999
C
      DO 300 NJSYM = 1,NSYM
         NISYM = MULD2H(NJSYM,NIJSYM)
         J1 = IORB(NJSYM) + 1
         JL = IORB(NJSYM) + NORB(NJSYM)
         I1 = IORB(NISYM) + 1
         IL = IORB(NISYM) + NORB(NISYM)
         DO 200 J = J1,JL
         IF (IOBTYP(J) .EQ. JTACT) THEN
            NJ = ICH(J)
            DO 100 I = I1,IL
            IF (IOBTYP(I) .EQ. JTACT) THEN
               NI = ICH(I)
               HHAC(NI,NJ) = HHAC(NI,NJ) + HH(I,J)
            END IF
  100       CONTINUE
         END IF
  200    CONTINUE
  300 CONTINUE
C
 9999 CALL QEXIT('ADDAC1')
      RETURN
C
C     End of ADDAC1.
C
      END
C  /* Deck tr1lfk */
      SUBROUTINE TR1LFK (NTMAT,TMAT,CMO,UDV,FXC,FXV,
     &                   LONDON,WRK,LWRK)
C
C     Hans Joergen Aa. Jensen and Rika Kobayashi 24-Apr-1992
C     (based on GETFCK by Henrik Koch and Trygve Helgaker Nov 91)
C
C     DFT modifications T. Helgaker
C
C     PURPOSE : Driver routine for the calculation of the two-electron
C               part of the fock matrices.
C               We assume the densities and fock matrices are full
C               squares and without symmtry reduction .
C               The subroutine offers the possibility of calculating
C               the fock matrices directly or by reading the integrals
C               from disk.
C
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
      PARAMETER (D2 = 2.0D0, DM1 = -1.0D0)
C
      LOGICAL LONDON, DFTADX
      DIMENSION TMAT(N2ORBX,NTMAT),CMO(*),UDV(*)
      DIMENSION FXC(NORBT,NORBT,NTMAT), FXV(NORBT,NORBT,NTMAT)
      DIMENSION WRK(LWRK)
C
C Used from common blocks:
C   INFINP : DIRFCK
C   INFORB : N2ORBX,...
C   INFVAR : JWOPSY
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <dftcom.h>
C
      PARAMETER (MXTMAT = 50)
      DIMENSION ISYMDM(MXTMAT), IFCTYP(MXTMAT)
C
      IF (NTMAT .GT. MXTMAT) THEN
         CALL QUIT('FATAL ERROR: increase MXTMAT parameter in TR1LFK')
      END IF
C
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
C
      NXMAT = 0
      IF (NISHT .GT. 0) NXMAT = NXMAT + NTMAT
      IF (NASHT .GT. 0) NXMAT = NXMAT + NTMAT
      LDFAO = NXMAT*N2BASX
      CALL MEMGET('REAL',KDXCAO,LDFAO,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KFXCAO,LDFAO,WRK,KFREE,LFREE)
C     Note: DXCAO and DXVAO must be contiguous.
C     Note: FXCAO and FXVAO must be contiguous.
      IF (NISHT .GT. 0) THEN
         KDXVAO = KDXCAO + NTMAT*N2BASX
         KFXVAO = KFXCAO + NTMAT*N2BASX
      ELSE
         KDXVAO = KDXCAO
         KFXVAO = KFXCAO
      END IF
C
      CALL LFKDEN(NTMAT,CMO,TMAT,UDV, WRK(KDXCAO),WRK(KDXVAO),
     &            LONDON,WRK(KFREE),LFREE)
C
      CALL DZERO(WRK(KFXCAO),LDFAO)
C
      IF (LONDON) THEN
CTROND         IFCSET = 2
C     Scale anti-symmetric density to "get things right..."
         IFCSET = 22
         CALL DSCAL(LDFAO,DM1,WRK(KDXCAO),1)
      ELSE
CTROND         IFCSET = 1
         IFCSET = 13
      END IF
      DO 100 I = 1,NXMAT
         ISYMDM(I) = JWOPSY
         IFCTYP(I) = IFCSET
  100 CONTINUE
      DFTADX = DFTADD
      IF (LONDON) DFTADD = .FALSE.
      CALL SIRFCK(WRK(KFXCAO),WRK(KDXCAO),NXMAT,ISYMDM,IFCTYP,
     &            DIRFCK,WRK(KFREE),LFREE)
      DFTADD = DFTADX
#if defined (VAR_OLDCODE)
C
CHJAAJ: this is from FCKDIR; change it !
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
      IF (IPRINT.GT.2) THEN
         CALL HEADER('AO density and Fock matrices in TR1LFK',-1)
         DO 300 I = 1, NXMAT
            JDXAO = ?
            JFXAO = ?
            WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
            CALL OUTPUT(WRK(JDXAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(WRK(JFXAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
#endif
C
C     Transform Fock matrices to MO basis
C     FC(J,I) = FC(J,I) + sum(a,b) CMO(a,J) FCAO(a,b) CMO(b,I)
C
      DO 4000 ISYM = 1,NSYM
         JSYM  = MULD2H(ISYM,JWOPSY)
         ICMOI = ICMO(ISYM)
         IORBI = IORB(ISYM)
         NORBI = NORB(ISYM)
         IBASI = IBAS(ISYM)
         NBASI = NBAS(ISYM)
         ICMOJ = ICMO(JSYM)
         IORBJ = IORB(JSYM)
         NORBJ = NORB(JSYM)
         IBASJ = IBAS(JSYM)
         NBASJ = NBAS(JSYM)
         JFXCAO = KFXCAO + NBAST*IBASI + IBASJ
         JFXVAO = KFXVAO + NBAST*IBASI + IBASJ
         IF (NORBI*NORBJ*NBASI*NBASJ .NE. 0) THEN
         DO 3800 ITMAT = 1,NTMAT
         IF (NISHT .GT. 0) THEN
            CALL DGEMM('T','N',NORBJ,NBASI,NBASJ,1.D0,
     &                 CMO(ICMOJ+1),NBASJ,
     &                 WRK(JFXCAO),NBAST,0.D0,
     &                 WRK(KDXCAO),NORBJ)
            CALL DGEMM('N','N',NORBJ,NORBI,NBASI,1.D0,
     &                 WRK(KDXCAO),NORBJ,
     &                 CMO(ICMOI+1),NBASI,1.D0,
     &                 FXC(IORBJ+1,IORBI+1,ITMAT),NORBT)
            JFXCAO = JFXCAO + N2BASX
         END IF
         IF (NASHT .GT. 0) THEN
            CALL DGEMM('T','N',NORBJ,NBASI,NBASJ,1.D0,
     &                 CMO(ICMOJ+1),NBASJ,
     &                 WRK(JFXVAO),NBAST,0.D0,
     &                 WRK(KDXVAO),NORBJ)
            CALL DGEMM('N','N',NORBJ,NORBI,NBASI,1.D0,
     &                 WRK(KDXVAO),NORBJ,
     &                 CMO(ICMOI+1),NBASI,1.D0,
     &                 FXV(IORBJ+1,IORBI+1,ITMAT),NORBT)
            JFXVAO = JFXVAO + N2BASX
         END IF
 3800    CONTINUE
         END IF
C **  this symmetry block finished
 4000 CONTINUE
C
      CALL MEMREL('TR1LFK',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C  /* Deck lfkden */
      SUBROUTINE LFKDEN(NTMAT,CMO,TMAT,UDV,DXCAO,DXVAO,
     &                  LONDON,WRK,LWRK)
C
C Hans Joergen Aa. Jensen and Rika Kobayashi 24-Apr-1992
C (based on DEQ27 JAN-20 1992, H.A.)
C P.Joergensen and K.Ruud, Nov.-93 (non-symmetric TD-matrix)
C Revised 14-Oct-1994 HJAaJ
C   (finished the generalization initiated by p.j. & k.r.,
C    in particular removed factor 0.5 and redundant calc.s)
C
C Purpose:
C To form the reorthonormalisation part of the derivative
C Fock matrices
C
C It follows TR1DEN in SIRIUS, which in turn
C follows appendix C in Chem. Phys. 104 (1986) 229.
C
C DXCAO and DXVAO must be allocated as square (unfolded) matrices
C
C General formula:   DX(p,q) = sum(t) (T(t,p)  D(t,q) + T*(t,q) D(p,t))
C                            = sum(t) (T*(p,t) D(t,q) + T(q,t)  D(p,t))
C (First line corresponds to the original code, second line
C is the new modified formulas by K.Ruud, May 1994.  Both forms
C are equivalent for symmetric T matrix, i.e. NODIFC true, but
C only the second form is correct if DIFC reorthonormalization /hjaaj)
C
C Specific formulas for real and imaginary perturbations:
C real (not london)  DX(p,q) =   sum(t) ( T(p,t) D(t,q) + T(q,t) D(p,t))
C imag. (london)   i DX(p,q) = i sum(t) (-T(p,t) D(t,q) + T(q,t) D(p,t))
C Thus for derivative/London orbitals, i.e. real/imag. :
C       DXC(i,q) = 2*T(q,i)
C       DXC(p,q) = DXC(p,q) +/- DXC(q,p)
C       DXV(u,q) = sum(v) T(q,v) UDV(u,v)
C       DXV(p,q) = DXV(p,q) +/- DXV(q,p)
C
C Input:
C   CMO(*)  molecular orbital coefficients
C   TMAT    Reorthonormalisation matrix
C   UDV(*)  active part of one-electron density matrix (over MO's)
C
C Output:
C   DXCAO(*) contravariant one-index transformed DC matrix
C   DXVAO(*) contravariant one-index transformed DV matrix
C
C Scratch:
C   WRK(LWRK)
C
#include <implicit.h>
      LOGICAL LONDON
      DIMENSION CMO(*),TMAT(NORBT,NORBT,NTMAT),UDV(NASHT,*),
     &          DXCAO(NBAST,NBAST,NTMAT),DXVAO(NBAST,NBAST,NTMAT),WRK(*)
C
      PARAMETER ( D2 = 2.0D0, D4 = 4.0D0 )
C
C Used from common blocks:
C  INFORB : NSYM,NASHT,...
C  INFVAR : JWOPSY,??????
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
C
#if defined (VAR_DEBUG)
#include <idbg.h>
#endif
C
      CALL QENTER('LFKDEN')
C
C
      KWRK = 1
      KFREE  = KWRK
      LFREE  = LWRK
      LDXAO1 = MAX(NASHT,NISHT)*NBAST
      CALL MEMGET('REAL',KDXAO1,LDXAO1,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KDXAO2,NNBASX,WRK,KFREE,LFREE)
C
C     *************************************************
C        Active matrices
C     *************************************************
C
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDXVT,NASHT*NORBT,WRK,KFREE,LFREE)
         CALL DZERO(DXVAO,NTMAT*N2BASX)
      DO 2000 ISYM = 1,NSYM
         JSYM  = MULD2H(ISYM,JWOPSY)
         NASHI = NASH(ISYM)
         NORBJ = NORB(JSYM)
         IF (NASHI .EQ. 0 .OR. NORBJ .EQ. 0) GO TO 2000
         NISHI = NISH(ISYM)
         IASHI = IASH(ISYM)
         ICMOI = ICMO(ISYM)
         IORBI = IORB(ISYM)
         NBASI = NBAS(ISYM)
         ICMOJ = ICMO(JSYM)
         IORBJ = IORB(JSYM)
         NBASJ = NBAS(JSYM)
C **     step 1 of active dm:
C *****  first calculate one-index transformation of first index
C *****  of UDV(uv), the active density matrix.
C        DXV(u,q) = sum(v) T(q,v) UDV(u,v) = sum(v) T(q,v) UDV(v,u)
C        DXVT(q,u) = DXV(u,q) = sum(v) T(q,v) UDV(v,u)
      IF (NASHI*NORBJ*NBASJ*NBASI .NE. 0) THEN
      DO 1800 ITMAT = 1,NTMAT
         CALL DGEMM('N','N',NORBJ,NASHI,NASHI,1.D0,
     &              TMAT(IORBJ+1,IORBI+NISHI+1,ITMAT),NORBT,
     &              UDV(IASHI+1,IASHI+1),NASHT,0.D0,
     &              WRK(KDXVT),NORBJ)
C **     step 2 of active dm:
C        DXVAO(a,b) = sum(u) CMO(a,u) (sum(q) CMO(b,q) DXVT(q,u))
         CALL DGEMM('N','N',NBASJ,NASHI,NORBJ,1.D0,
     &              CMO(ICMOJ+1),NBASJ,
     &              WRK(KDXVT),NORBJ,0.D0,
     &              WRK(KDXAO1),NBASJ)
         IOFMOV = ICMOI + 1 + NISHI*NBASI
         CALL DGEMM('N','T',NBASI,NBASJ,NASHI,1.D0,
     &              CMO(IOFMOV),NBASI,
     &              WRK(KDXAO1),NBASJ,0.D0,
     &              DXVAO(IBAS(ISYM)+1,IBAS(JSYM)+1,ITMAT),NBAST)
 1800 CONTINUE
      END IF
C **  this symmetry block finished
 2000 CONTINUE
C
C **     step 3 of active dm:
C        DXV(p,q) = DXV(p,q) +/- DXV(q,p)
C        Scale by two to get DXVAO(i,j) = (DXVAO(i,j) - DXVAO(j,i))
C        below
C
         CALL DSCAL(NTMAT*N2BASX,D2,DXVAO,1)
C
C        Obtain DXVAO(i,j) = 0.5*(DXVAO(i,j) -/+ DXVAO(j,i))
C
         DO 2200 ITMAT = 1,NTMAT
            IF (LONDON) THEN
               CALL DGETAP(NBAST,DXVAO(1,1,ITMAT),WRK(KDXAO2))
               CALL DAPTGE(NBAST,WRK(KDXAO2),DXVAO(1,1,ITMAT))
            ELSE
               CALL DGETSP(NBAST,DXVAO(1,1,ITMAT),WRK(KDXAO2))
               CALL DSPTSI(NBAST,WRK(KDXAO2),DXVAO(1,1,ITMAT))
            END IF
 2200    CONTINUE
C
      END IF
C
C     *************************************************
C        Inactive matrices
C     *************************************************
C
      IF (NISHT .GT. 0) THEN
         CALL DZERO(DXCAO,NTMAT*N2BASX)
      DO 4000 ISYM = 1,NSYM
         JSYM  = MULD2H(ISYM,JWOPSY)
         NISHI = NISH(ISYM)
         NORBJ = NORB(JSYM)
      IF (NISHI .EQ. 0 .OR. NORBJ .EQ. 0) GO TO 4000
         ICMOI = ICMO(ISYM)
         IORBI = IORB(ISYM)
         NBASI = NBAS(ISYM)
         ICMOJ = ICMO(JSYM)
         IORBJ = IORB(JSYM)
         NBASJ = NBAS(JSYM)
C
C **     the inactive one-index transformed dm:
C        DXC(i,q) = 2*T(q,i)
C        DXCAO(a,b) = sum(i) CMO(a,i) (sum(q) CMO(b,q) DXC(i,q))
C                 = 2 sum(i) CMO(a,i) (sum(q) CMO(b,q) TMAT(q,i))
C        The factor 2 is postponed to the DSCAL below
         IF (NBASJ*NISHI*NORBJ*NBASI .NE. 0) THEN
      DO 3800 ITMAT = 1,NTMAT
         CALL DGEMM('N','N',NBASJ,NISHI,NORBJ,1.D0,
     &              CMO(ICMOJ+1),NBASJ,
     &              TMAT(IORBJ+1,IORBI+1,ITMAT),NORBT,0.D0,
     &              WRK(KDXAO1),NBASJ)
C        ... DXAO1(alpha,i) = CMO(alpha,p) TMAT(p,i)
         CALL DGEMM('N','T',NBASI,NBASJ,NISHI,1.D0,
     &              CMO(ICMOI+1),NBASI,
     &              WRK(KDXAO1),NBASJ,0.D0,
     &              DXCAO(IBAS(ISYM)+1,IBAS(JSYM)+1,ITMAT),NBAST)
C        ... DXCAO(beta,alpha) = CMO(beta,i) DXAO1(alpha,i)
 3800 CONTINUE
      END IF
C **  this symmetry block finished
 4000 CONTINUE
C
C        Scale by four to get final DXCAO(i,j) below
C        Factor 2 from the second one-index transformation (as for
C        DXVAO), and factor 2 from DC(i,i) = 2
C
         CALL DSCAL(NTMAT*N2BASX,D4,DXCAO,1)
C
C        DXC(p,q) = DXC(p,q) +/- DXC(q,p)
C        Obtain second one-index transformation from
C           DXCAO(i,j) = 0.5*(DXCAO(i,j) -/+ DXCAO(j,i))
C
         DO 4200 ITMAT = 1,NTMAT
            IF (LONDON) THEN
               CALL DGETAP(NBAST,DXCAO(1,1,ITMAT),WRK(KDXAO2))
               CALL DAPTGE(NBAST,WRK(KDXAO2),DXCAO(1,1,ITMAT))
            ELSE
               CALL DGETSP(NBAST,DXCAO(1,1,ITMAT),WRK(KDXAO2))
               CALL DSPTSI(NBAST,WRK(KDXAO2),DXCAO(1,1,ITMAT))
            END IF
 4200    CONTINUE
C
      END IF
C
C
      CALL MEMREL('LFKDEN',WRK,KWRK,KWRK,KFREE,LFREE)
      CALL QEXIT('LFKDEN')
      RETURN
C
C *** end of subroutine LFKDEN
C
      END
