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
#ifdef UNDEF
===========================================================================
901023-hjaaj: new parameter list to CIDIAG in ABADIA
900216-hjaaj: abalr: exit if ingd true (not reimplemented yet)
 -- abalr: set nload = 1 after resst if no trial vectors found
 (will make anrctl check if already converged without calling abalin)
      END
===========================================================================
#endif
C  /* Deck abalr */
      SUBROUTINE ABALR(LMAXVE,INGD,REDH,REDGD,SOLEQ,IBNDX,CMO,UDV,PV,
     *                 FC,FV,FCAC,H2AC,CREF,GORB,DIAOR,DIACI,INDXCI,
     *                 GD,WRK,LWRK)
C
C Purpose: solve linear response equations in Abacus.
C If INGD true, then solve for the gradient vector in GD,
C else solve for the gradient vectors on LUGDR.
C Note: GD is only used if INGD true
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
#include <iratdef.h>
C
      DIMENSION REDH(*),REDGD(*),SOLEQ(*),IBNDX(*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION CREF(*),GORB(*),DIAOR(*),DIACI(*),INDXCI(*),WRK(*)
      DIMENSION GD(1)
C
#include <abainf.h>
#include <cbires.h>
#include <inftap.h>
#include <nuclei.h>
#include <gdvec.h>
#include <infpri.h>
#include <inflin.h>
#include <infdim.h>
#include <inforb.h>
C
C     Declare logical external functions
C
      LOGICAL ANRRST
C
C     Local variables
C
      LOGICAL INGD, LRSTNR, OPTORB
C
      CALL QENTER('ABALR')
      IF (INGD) THEN
         CALL QUIT('ERROR in ABALR: INGD not implemented.')
      END IF
      LRSTNR = RSTNR
      IF (KAPTST .GT. 0) THEN
         OPTORB = .TRUE.
      ELSE IF (KAPTST .LT. 0) THEN
         OPTORB = .FALSE.
      ELSE
         OPTORB = NCONMA .GT. 1000
      END IF
      KDV    = 1
      KKNR   = KDV  + NNASHX
      KGD    = KKNR + LMAXVE
      LGD    = LWRK - KGD
C
      IF (NASHT .GT. 0) CALL DSITSP(NASHT,UDV,WRK(KDV))
C
      NREDH  = 0
      IF (IPRLIN.GT.100) THEN
         WRITE(LUPRI,'(/A)')' **** ABALR **** '
         IF (NASHT.GT.0) THEN
            WRITE(LUPRI,'(/A)')
     *      ' ONE ELECTRON ACTIVE DENSITY MATRIX  PACKED'
            CALL OUTPAK(WRK(KDV),NASHT,1,LUPRI)
            WRITE(LUPRI,'(/A)')
     *      ' ONE ELECTRON ACTIVE DENSITY MATRIX  UNPACKED'
            CALL OUTPUT(UDV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(/A)') ' FOCK VALENCE MATRIX'
            CALL OUTPKB(FV,NORB,NSYM,1,LUPRI)
         END IF
         WRITE(LUPRI,'(/A)') ' FOCK CORE MATRIX'
         CALL OUTPKB(FC,  NORB,NSYM,1,LUPRI)
      END IF
C
      IF ( IPRLIN.GT.5 ) WRITE(LUPRI,'(/A,L8)')' LRSTNR:',LRSTNR
      IF (LRSTNR) THEN
         LRSTNR = ANRRST(NREDH, IBNDX, REDH)
         IF (LRSTNR) THEN
           IF (IPRINT .GE. 1) WRITE(LUPRI,'(/A,I4,A/)')
     *        ' Restart with',NREDH,' "old" trial vectors.'
           NLOAD = 1
         ELSE
            WRITE (LUPRI,'(/A)')
     *      ' Restart not possible for response vector calculation'
              WRITE(LUPRI,'(A)') ' ".NRRESTART" option ignored.'
         END IF
      END IF
      CALL GPOPEN(LUNR3,ABANR3,'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUNR5,ABANR5,'UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      IF (.NOT. LRSTNR) THEN
C        IF (INGD) CALL DCOPY(NVARPT,WRK(KGD),1,WRK(KWRK1),1)
         NRSTV = NGDVEC(LSYMPT,1)
         CALL RESST(INGD,NRSTV,CREF,DIACI,DIAOR,IBNDX,WRK(KGD),LGD)
C        CALL RESST(INGD,NRSTV,CREF,DIACI,DIAOR,IBNDX,WRK,LWRK)
         NREDH = NRSTV
         NLOAD = 0
         IF (NREDH .EQ. 0) NLOAD = 1
C        ... no start trial vectors found, skip abalin in anrctl
C            and check if converged.
      END IF
      LNRDT = NGDVEC(LSYMPT,1)
      DO 100 ISIM = 1, LNRDT, LMAXVE
         NLOAD = NLOAD + 1
         IRDC = ISIM - 1
         NRDC = MIN(LMAXVE,(LNRDT - IRDC))
C
         KWRK1 = KGD  + NRDC*NVARPT
         LWRK1 = LWRK - KWRK1
C
C           Read gradient vectors:
C
         IOFFGD = KGD
         INOT = 0
         DO 200 I = 1, NRDC
            IREC = IGDREC(IRDC + I,LSYMPT,1)
            IF (MOLHES .OR. (POLAR .AND. IREC .GT. 3*NUCDEP) .OR.
     &         (QPGRAD .AND. IREC .GT. (3*NUCDEP + 3))) THEN
               CALL READDX (LUGDR,IREC,IRAT*NVARPT,WRK(IOFFGD))
               IOFFGD = IOFFGD + NVARPT
            ELSE
               INOT = INOT + 1
            END IF
200      CONTINUE
         NRDC = NRDC - INOT
         IF ( IPRLIN.GT.40 ) THEN
            CALL HEADER('GD vectors in ANRCTL',-1)
            CALL OUTPUT(WRK(KGD),1,NVARPT,1,NRDC,NVARPT,NRDC,1,LUPRI)
         END IF
C 
C In cases where many responses need not be solved (e.g. APT's)
C a batch, as this causes the program to crash. K.Ruud Oct.-00 
C
         IF (NRDC .EQ. 0) THEN
            NLOAD = NLOAD - 1
            GOTO 100
         END IF
C
C        Solve linear equations:
C
         CALL ANRCTL(LMAXVE,REDH,REDGD,NREDH,NRDC,NLOAD,
     &               OPTORB,SOLEQ,CMO,CREF,
     &               GORB,WRK(KDV),PV,FC,FV,FCAC,H2AC,
     &               DIACI,DIAOR,INDXCI,IBNDX,WRK(KKNR),WRK(KGD),
     &               WRK(KWRK1),LWRK1,MAXRED)
C        CALL ANRCTL(MAXSIM,REDH,REDGD,NREDH,NGD,NLOAD,DOXKAP,
C    &               SOLEQ,CMO,CREF,GORB,DV,PV,FC,FV,FCAC,H2AC,
C    &               DIACI,DIAOR,INDXCI,IBNDX,KNR,GD,WRK,LWRK,MAXRED)
C
C        Construct residual vectors in GD and write on file
C
         CALL ABARES(SOLEQ,WRK(KGD),NRDC,NREDH,WRK(KWRK1),LWRK1)
C        CALL ABARES(SOLEQ,GD,NBX,NREDH,WRK,LWRK)
         JRD = KGD
         DO 400 I = 1,NRDC + INOT
            IREC = 2*IGDREC(IRDC + I,LSYMPT,1)
            IF (MOLHES .OR. (POLAR .AND. IREC .GT. 6*NUCDEP).OR.
     &         (QPGRAD .AND. IREC .GT. (6*NUCDEP + 6)))
     &         CALL WRITDX (LURDR,IREC,IRAT*NVARPT,WRK(JRD))
            JRD = JRD + NVARPT
 400     CONTINUE
C
C        Construct solution vectors in GD
C
         CALL ANRCON(NREDH,NRDC,SOLEQ,WRK(KGD),IBNDX,WRK(KWRK1))
C        CALL ANRCON(NREDH,NRD,SOLEQ,RD,IBNDX,WRK)
C
         IF (IPRINT .GE. 10) THEN
            CALL AROUND('Solution vectors written on LURDR in ABALR')
         END IF
         JRD = KGD
         DO 300 I = 1, NRDC + INOT
            IREC = 2*IGDREC(IRDC + I,LSYMPT,1) - 1
            IF (MOLHES .OR. (POLAR .AND. IREC .GT. 6*NUCDEP).OR.
     &         (QPGRAD .AND. IREC .GT. (6*NUCDEP + 6))) THEN
               CALL WRITDX (LURDR,IREC,IRAT*NVARPT,WRK(JRD))
               IF (IPRINT .GE. 10) THEN
                  WRITE (LUPRI,'(//A,I5)') ' Symmetry coordinate: ',
     *                                      IGDCOR(IRDC + I,LSYMPT,1)
                  CALL HEADER('Orbital Solution',-1)
                  IF (IPRINT .GE. 20) THEN
                     PRFAC = 0.0D0
                  ELSE
                     PRFAC = 0.1D0
                  END IF
                  CALL PRKAP(NWOPPT,WRK(JRD + NCONST),PRFAC,LUPRI)
                  IF ((IPRINT .GE. 40) .AND. (NCONST .GT. 0)) THEN
                   CALL HEADER('CI Solution',-1)
                   WRITE (LUPRI,'(1P,6D12.4)') (WRK(JRD-1+K),K=1,NCONST)
                  END IF
               END IF
               JRD = JRD + NVARPT
               IREC = IGDREC(IRDC + I,LSYMPT,1)
               IDORCI(IREC,1) = 2
            END IF
 300     CONTINUE
C     
 100  CONTINUE
      CALL GPCLOSE(LUNR3,'DELETE')
      CALL GPCLOSE(LUNR5,'DELETE')
C
C *** end of ABALR --
C
      CALL QEXIT('ABALR')
      RETURN
      END
C  /* Deck abares */
      SUBROUTINE ABARES(SOLEQ,GD,NBX,NREDH,WRK,LWRK)
C
C     Purpose: Construct residual RD(I) = H(2)*X(I) + GD
C              Residual is returned in GD
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION SOLEQ(NREDH,*)
      DIMENSION GD(*),WRK(*)
      PARAMETER (  DM1=-1.0D0  )

C Used from common blocks:
C
C  INFLIN: NVARPT,NVARPT
C  INFTAP: LUNR5
C  LINABA: LOFFRD
C
#include <linaba.h>
#include <infpri.h>
#include <inflin.h>
#include <inftap.h>
C
#include <ibndxdef.h>
C
C ADD  L[2]*X(I) to GD
C
      IF (LWRK .LT. NVARPT) CALL ERRWRK('ABARES',NVARPT,LWRK)
      REWIND LUNR5
      IF (LOFFRD.EQ.1) THEN
         IF (NCONST.GT.0) READ(LUNR5)
         IF (NWOPPT.GT.0) READ(LUNR5)
      END IF
      DO 600 K = 1,NREDH
         IF (NCONST.GT.0 ) CALL READT(LUNR5,NCONST,WRK(1))
         IF (NWOPPT.GT.0 ) CALL READT(LUNR5,NWOPPT,WRK(1+NCONST))
         DO 700 JR = 1,NBX
            JRDOFF  = (JR-1)*NVARPT + 1
            FAC1   = SOLEQ(K,JR)
            CALL DAXPY(NVARPT,FAC1,WRK(1),1,GD(JRDOFF),1)
 700     CONTINUE
 600  CONTINUE
      IF (IPRLIN.GT.100) THEN
         WRITE (LUPRI,*) ' NBX RESIDUAL VECTORS',NBX
         CALL OUTPUT(GD,1,NVARPT,1,NBX,NVARPT,NBX,1,LUPRI)
      END IF
      RETURN
C
C     End of ABARES
C
      END
C  /* Deck abavar */
      SUBROUTINE ABAVAR(JSYMPT,TRPCLC,JPRINT,WRK,LWRK)
C
C  l.r. 900106
C  Set variables and index arrays to perturbation symmetry JSYMPT
C  970815sps+hjaaj: added TRPCLC, JPRINT parameters and soppa code
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
C
      LOGICAL TRPCLC
      DIMENSION WRK(LWRK)
C
C Used from common blocks:
C  INFORB : MULD2H
C  INFDIM : NCONDI, NWOPDI
C  INFTAP : LUINDX
C  INFVAR : JWOPSY, NCONF, NWOPT, NWOPH, NVAR, NVARH
C  INFLIN : LSYMRF, LSYMST, LSYMPT, NCONST, NWOPPT, NVARPT
C  LINABA : ABAHF, ABACI, LOFFRD
C  ABAINF : ABASOP
C  INFSOP : N2P2HS, N2P2HT, NS2P2H, NT2P2H, N12P2H, N22P2H, N32P2H
C
#include <inforb.h>
#include <infdim.h>
#include <inftap.h>
#include <infvar.h>
#include <inflin.h>
#include <linaba.h>
#include <abainf.h>
#include <infsop.h>
C
      LOGICAL GETWOP, GOTWOP
C
      CHARACTER*8 OPLBL(8)
      DATA OPLBL/'EXOPSYM1','EXOPSYM2','EXOPSYM3','EXOPSYM4',
     *           'EXOPSYM5','EXOPSYM6','EXOPSYM7','EXOPSYM8'/
C
C     *******************************************
C     ***** Number of variatonal parameters *****
C     *******************************************
C
      LSYMPT = JSYMPT
      LSYMST = MULD2H(LSYMRF,LSYMPT)
C
C     Configuration variables
C     -----------------------
C
      IF (ABAHF) THEN
         NCONF  = 0
         IF (ABASOP) THEN
            IF (TRPCLC) THEN
               NCONST = N2P2HT(LSYMPT)
            ELSE
               NCONST = N2P2HS(LSYMPT)
            ENDIF
         ELSE
            NCONST = 0
         ENDIF
      ELSE
         CALL SETCI2(NCONF,LSYMST,TRPCLC,0)
         NCONST = NCONF
      END IF
      NCONDI = MAX(1,NCONF)
C
C     Orbital rotation variables
C     --------------------------
C
      IF ( LSYMPT .NE. JWOPSY ) THEN
         IF (ABACI) THEN
            NWOPT  = 0
         ELSE
            GOTWOP = GETWOP(OPLBL(LSYMPT))
            IF (.NOT.GOTWOP) CALL QUIT
     &         ('ABAVAR error: GETWOP cannot find OPLBL(LSYMPT)')
         END IF
         JWOPSY = LSYMPT
         NWOPDI = MAX(1,NWOPT)
      END IF
      NWOPPT = NWOPT
C
C     Total number of variables
C     -------------------------
C
      NVAR    = NCONF  + NWOPT
      NVARH   = NCONF  + NWOPH
      NVARPT  = NCONST + NWOPPT
      IF (JPRINT .GE. 3) THEN
         IF (ABASOP) THEN
            WRITE(LUPRI,'(/4(/A,I8))')
     *       ' Perturbation symmetry ............ :',LSYMPT,
     *       ' p-h variables.            - NWOPPT :',NWOPPT,
     *       ' 2p-2h variables.          - NCONST :',NCONST,
     *       ' Total number of variables - NVARPT :',NVARPT
            IF (TRPCLC) THEN
               IF (JPRINT .GE. 6) WRITE(LUPRI,'(3(/A,I8))')
     *         ' T(1) 2p-2h variables.         :',N12P2H(LSYMPT),
     *         ' T(2) 2p-2h variables.         :',N22P2H(LSYMPT),
     *         ' T(3) 2p-2h variables.         :',N32P2H(LSYMPT)
            ELSE
               IF (JPRINT .GE. 6) WRITE(LUPRI,'(2(/A,I8))')
     *         ' R(1) 2p-2h variables.         :',NS2P2H(LSYMPT),
     *         ' R(2) 2p-2h variables.         :',NT2P2H(LSYMPT)
            ENDIF
         ELSE
            WRITE(LUPRI,'(/4(/A,I8))')
     *       ' Perturbation symmetry ............ :',LSYMPT,
     *       ' Orbital variables         - NWOPPT :',NWOPPT,
     *       ' Configuration variables   - NCONST :',NCONST,
     *       ' Total number of variables - NVARPT :',NVARPT
         ENDIF
      ENDIF
C
C     Check if the reference state function has to be read in as
C     the first trial vector on LUNR3.
C
      IF (( NCONF .GT. 0 ) .AND. ( LSYMPT .EQ. 1 )) THEN
C           ^- was NCONST, now NCONF because of ABASOP
         LOFFRD = 1
      ELSE
         LOFFRD = 0
      ENDIF
C
C *** End of ABAVAR
C
      RETURN
      END
C  /* Deck abadia */
      SUBROUTINE ABADIA(UDV,FOCK,FC,FV,FCAC,H2AC,XNDXCI,WRK,LWRK)
C
C l.r. 901023 hjaaj
C Calculate diagonal of MCSCF Hessian
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
      LOGICAL NOH2
      PARAMETER (NOH2 = .FALSE.)
C
      DIMENSION UDV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XNDXCI(*), WRK(*)
C
C Used from common blocks:
C  INFORB : NSYM, NASHT, NORB(), ?
C  INFOPT : EACTIV
C  INFLIN : NCONST, NWOPPT
C  LINABA : AVEODG
C
#include <inforb.h>
#include <infopt.h>
#include <inftap.h>
#include <inflin.h>
#include <linaba.h>
C
C     ****************************************
C     ***** Diagonal elements of Hessian *****
C     ****************************************
C
C
C     Diagonal elements of configuration Hessian
C     ------------------------------------------
C
      KDIACI = 1
      IF (NCONST.GE.1) THEN
         KWRK1  = KDIACI + NCONST
         LWRK1  = LWRK   - KWRK1
         LUIT2S = LUIT2
         LUIT2  = -1
C        ... set LUIT2 = 0 so that diagonal isn't written to LUIT2
         CALL CIDIAG(LSYMST,NOH2,FCAC,H2AC,XNDXCI,WRK(KDIACI),
     &               WRK(KWRK1),LWRK1)
         LUIT2  = LUIT2S
         DO 50 I = 1, NCONST
            WRK(KDIACI - 1 + I) = WRK(KDIACI - 1 + I) - EACTIV
 50      CONTINUE
         IF (IPRLIN.GT.100) THEN
            CALL AROUND('Diagonal elements of configuration Hessian')
            CALL OUTPUT(WRK(KDIACI),1,NCONST,1,1,NCONST,1,1,LUPRI)
         END IF
      END IF
C
C     Fock matrix contribution to diagonal elements of orbital Hessian
C     ----------------------------------------------------------------
C
      IF (NWOPPT.GT.0) THEN
         IF (IPRLIN.GT.100) THEN
            IF (NASHT.GT.0) THEN
               CALL AROUND('One-electron active density '/
     *                    /'matrix unpacked in ABADIA')
               CALL OUTPUT(UDV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
               CALL AROUND('Fock valence matrix in ABADIA')
               CALL OUTPKB(FV,NORB,NSYM,1,LUPRI)
            END IF
            CALL AROUND('Fock core matrix in ABADIA')
            CALL OUTPKB(FC,  NORB,NSYM,1,LUPRI)
            CALL AROUND('Total Fock matrix in ABADIA')
            CALL OUTPTB(FOCK,NORB,NSYM,1,LUPRI)
         END IF
         KDIAOR = KDIACI + NCONST
         CALL ABAODG(FOCK,FC,FV,UDV,WRK(KDIAOR))
C        CALL ABAODG(FOCK,FC,FV,UDV,EODIA)
C
         IF (IPRLIN .GT. 100) THEN
            CALL AROUND('Fock matrix contribution to diagonal '/
     *                 /'orbital Hessian')
            IF (AVEODG) WRITE (LUPRI,'(A)')
     *      ' (Fock type decoupling of two el. density matrix used)'
            IF (IPRLIN .GT. 120) THEN
               PRFAC = 0.0D0
            ELSE
               PRFAC = 0.1D0
            END IF
            CALL PRKAP(NWOPPT,WRK(KDIAOR),PRFAC,LUPRI)
         END IF
      END IF
C
C
C *** End of ABADIA
C
      RETURN
      END
C  /* Deck abaodg */
      SUBROUTINE ABAODG(FOCK,FC,FV,UDV,EODIA)
C
C     WRITTEN 13-MAY 1987
C
C     Purpose:
C       Calculate Fock contributions to diagonal elements of Hessian
C
C       EODIA(L,K) = < [E(K,L),H,E(L,K)] >     L > K
C
C    SECONDARY - INACTIVE
C     (L,K) = (A,I)
C    EODIA(A,I) = -FOCK(I,I) + 2*FC(A,A) +2*FV(A,A)
C
C    SECONDARY - ACTIVE
C             (A,M)
C    EODIA(A,M) = UDV(M,M)*FC(A,A) - FOCK(M,M)
C
C    ACTIVE - INACTIVE
C             (M,I)
C    EODIA(M,I) = 2*FC(M,M) + UDV(M,M)*FC(I,I) - FOCK(I,I) - FOCK(M,M)
C               + 2*FV(M,M)
C
C    ACTIVE - ACTIVE
C             (M,N)
C    EODIA(M,N) = UDV(N,N)*FC(M,M) + UDV(M,M)*FC(N,N)
C               - FOCK(M,M) - FOCK(N,N)
C
C    AVEODG = .TRUE.:
C
C       Add FV contributions to EODIA, originating from Fock-type
C       decoupling of the two-electron density matrix. All UDV*FC
C       contributions then become UDV*(FC + FV)
C
C                ******************************
C
#include <implicit.h>
#include <priunit.h>
C
C
      DIMENSION FC(*),FV(*),FOCK(*),UDV(NASHDI,NASHDI)
      DIMENSION EODIA(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (D0 = 0.0D0, D2 =2.0D0)
      PARAMETER (ORBSFT = 1.0D-6)
C  ABAINF : NWNABA
C  INFDIM : NASHDI
C
#include <abainf.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infpri.h>
#include <inflin.h>
#include <linaba.h>
C
      IF (IPRLIN .GT. 140) THEN
         IF (NASHT.GT.0) THEN
            CALL AROUND('One-el. act. density matrix unpacked - ABAODG')
            CALL OUTPUT(UDV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
            CALL AROUND('Fock valence matrix - ABAODG')
            CALL OUTPKB(FV,NORB,NSYM,1,LUPRI)
         END IF
         CALL AROUND('Fock core matrix - ABAODG')
         CALL OUTPKB(FC,  NORB,NSYM,1,LUPRI)
         CALL AROUND('Total Fock matrix - ABAODG')
         CALL OUTPTB(FOCK,NORB,NSYM,1,LUPRI)
      END IF
      CALL DZERO(EODIA,NWOPPT)
      KSYM1 = 0
      DO 100 IG = 1,NWOPPT
        K     = JWOP(1,IG)
        L     = JWOP(2,IG)
        ISMOK = ISMO(K)
        ISMOL = ISMO(L)
        IF (ISMOK .NE. KSYM1) THEN
           KSYM1  = ISMOK
           NORBK  = NORB (ISMOK)
           IORBK  = IORB (ISMOK)
           IIORBK = IIORB(ISMOK)
           I2ORBK = I2ORB(ISMOK)
           IORBL  = IORB (ISMOL)
           NORBL  = NORB (ISMOL)
           IIORBL = IIORB(ISMOL)
           I2ORBL = I2ORB(ISMOL)
        END IF
        NK    = K - IORBK
        NL    = L - IORBL
        ITYPK = IOBTYP(K)
        ITYPL = IOBTYP(L)
        IIKK  = IIORBK + NK*(NK+1)/2
        IILL  = IIORBL + NL*(NL+1)/2
        I2KK  = I2ORBK + NORBK*(NK-1)+NK
        I2LL  = I2ORBL + NORBL*(NL-1)+NL
        NWL   = ISW(L) - NISHT
        NWK   = ISW(K) - NISHT
        IF (ITYPK .EQ. JTINAC) THEN
           EODIA(IG) = EODIA(IG) - FOCK(I2KK) + D2*(FC(IILL) + FV(IILL))
           IF (ITYPL .EQ. JTACT) THEN
              EODIA(IG) = EODIA(IG) - FOCK(I2LL) + UDV(NWL,NWL)*FC(IIKK)
              IF (AVEODG) EODIA(IG) = EODIA(IG)  + UDV(NWL,NWL)*FV(IIKK)
           END IF
        ELSE
           IF (ITYPL .EQ. JTACT) THEN
              EODIA(IG) = EODIA(IG) + UDV(NWK,NWK)*FC(IILL)
     *                              + UDV(NWL,NWL)*FC(IIKK)
     *                              - FOCK(I2LL) - FOCK(I2KK)
              IF (AVEODG) EODIA(IG) = EODIA(IG) + UDV(NWK,NWK)*FV(IILL)
     *                                          + UDV(NWL,NWL)*FV(IIKK)
           ELSE
              EODIA(IG) = EODIA(IG) + UDV(NWK,NWK)*FC(IILL) - FOCK(I2KK)
              IF (AVEODG) EODIA(IG) = EODIA(IG) + UDV(NWK,NWK)*FV(IILL)
           ENDIF
        ENDIF
 100  CONTINUE
C
C     Check if orbital Hessian has negative diagonal elements.
C     Add ORBSFT to remove spurious negative elements.
C
      NNEG = 0
      DO 300 I = 1,NWOPPT
         IF ( EODIA(I) .LT. D0 ) NNEG = NNEG + 1
         IF ( EODIA(I) .LT. ORBSFT ) EODIA(I) = ORBSFT
 300  CONTINUE
      IF (NNEG.GT.0) THEN
         WRITE(LUPRI,'(/A,/A,I5,A)') '  ********* WARNING *********** ',
     *   '  Diagonal orbital Hessian has ',NNEG,' negative elements.'
         NWNABA = NWNABA + 1
      END IF
      IF (IPRLIN.GT.100) THEN
         CALL AROUND('Diagonal elements of orbital Hessian')
         CALL OUTPUT(EODIA,1,NWOPPT,1,1,NWOPPT,1,1,LUPRI)
      END IF
C
C END OF ABAODG
C
      RETURN
      END
