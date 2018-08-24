!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C
C  /* Deck srestr */
      SUBROUTINE SRESTR
C
C     Several SIRIUS-variables has to be restored, because they are overwritten
C     by ABACUS in SETSIR module. Called from SIRRDI.
C
#include "implicit.h"
#include "maxorb.h"
#include "gnrinf.h"
#include "infinp.h"
#include "inpbak.h"
#include "infpri.h"
#include "pribak.h"
#include "pcmlog.h"
C
C     Restoring should only be done if SIRIUS has been through its
C     input processing.
C
      IF (.NOT. SIR_INPPRC) RETURN
C
C     /INFINP/ values are restored
C
      ISPIN  = ISPINB
      ISTATE = ISTATB
      LSYM   = LSYMB
      NACTEL = NACTLB
      MCTYPE = MCTYPB
      LSOLMX = LSLMXB
      NLMSOL = NLMSLB
      NELMN1 = NLMN1B
      NELMX1 = NLMX1B
      NELMN3 = NLMN3B
      NELMX3 = NLMX3B
      LROOTS = LROTSB
      NROOTS = NROTSB
      IORTO  = IORTOB
      ICI0   = ICI0B
      KDEL   = KDELB
      ICHECK = ICHCKB
      NTIT   = NTITB
      MAXMAC = MXMACB
      MAXMIC = MXMICB
      MAXJT  = MAXJTB
      MAXCIT = MXCITB
      MAXUIT = MXUITB
      MAXAPM = MXAPMB
      MAXABS = MXABSB
      ITRLVL = ITRLVB
      ITRFIN = ITRFNB
      JCHSYM = JCHSMB
      JCHORB = JCHRBB
      NROOCI = NROCIB
      ISTACI = ISTCIB
      NFIELD = NFILDB
      MXCIMA = MXCMAB
      ICICNO = ICCNOB
      IMCCNO = IMCNOB
      IF (NEWSYM) THEN
         DOSCF  = DOSCFB
         DOMP2  = DOMP2B
         DOCINO = DOCNOB
         DOCI   = DOCIB
         FCVORB = FCVRBB
      ELSE
         DOSCF  = (DOSCFB .AND. .NOT. DOMCB)
         DOMP2  = (DOMP2B .AND. .NOT. DOMCB)
         DOCINO = (DOCNOB .AND. .NOT. DOMCB)
         DOCI   = (DOCIB  .AND. .NOT. DOMCB)
         FCVORB = (FCVRBB .AND. .NOT. DOMCB)
      END IF
      DOMC   = DOMCB
      DORSP  = DORSPB
      LNOROT = LNROTB
      LMOORD = LMORDB
      DIRFCK = DRFCKB
      CORHOL = CORHLB
      CORRLX = CRRLXB
      RESPHP = RSPHPB
      JOLSEN = JOLSNB
      ABAIPH = ABIPHB
      INERSI = INRSIB
      INERSF = INRSFB
      SUPSYM = SUPSMB
      SPIN   = SPINB
      POTNUC = POTNCB
      EPSOL  = EPSOLB
      EPSTAT = EPSTTB
      PCM    = PCMB
      RSOL(1)= RSOLB(1)
      RSOL(2)= RSOLB(2)
      RSOL(3)= RSOLB(3)
      THRGRD = TRGRDB
      THRPWF = TRPWFB
      THRCI  = TRCIB
      THRMC  = TRMCB
      THRCGR = TRCGRB
      CMAXMO = CMXMOB
      THROVL = TROVLB
      THRSSY = TRSSYB
      DO I = 1, MAXRTS
         IROOT(I) = IROOTB(I)
      END DO
      DO I = 1, MXCORB
         NOROT(I)  = NOROTB(I)
         IMOORD(I) = IMORDB(I)
         CENT(I)   = CENTB(I)
         TYPE(I)   = TYPEB(I)
      END DO
      DO I = 1, NFLAG
         FLAG(I) = FLAGB(I)
      END DO
      DO I = 1, 6
         TITLE(I)      = TITLEB(I)
      END DO
      TITMOL(1)   = TITMLB(1)
      TITMOL(2)   = TITMLB(2)
      DO I = 1, MXFELT
         EFIELD(I) = EFILDB(I)
         LFIELD(I) = LFILDB(I)
      END DO
C
C     /INFPRI/ values are restored
C
      DO I = 1, NPFLAG
         P4FLAG(I) = P4FLGB(I)
         P6FLAG(I) = P6FLGB(I)
      END DO
      IPRI4  = IPRI4B - 2
      IPRI6  = IPRI6B
      IPRSIR = IPRSRB
      IPRCNO = IPRCNB
      IPRDIA = IPRDIB
      IPRSIG = IPRSGB
      IPRDNS = IPDNSB
      IPRSOL = IPRSLB
      IPRKAP = IPRKPB
      IPRMUL = IPRMLB
      IPRCIX = IPRCXB
      IPRRHF = IPRHFB
      IPRFCK = IPRFCB
      IPRAVE = IPRAVB
      MPRI4  = MPRI4B
      MPRI6  = MPRI6B
      MPRSIR = MPRSRB
      RETURN
      END
C  /* Deck ssave */
      SUBROUTINE SSAVE
C
C     Several SIRIUS-variables has to be saved, because they are overwritten
C     by ABACUS in SETSIR module. Called from SIRRDI.
C
#include "implicit.h"
#include "maxorb.h"
#include "gnrinf.h"
#include "infinp.h"
#include "inpbak.h"
#include "infpri.h"
#include "pribak.h"
#include "pcmlog.h"
C
C        /INFINP/ values are saved to /INFBAK/
C
      ISPINB = ISPIN
      ISTATB = ISTATE
      LSYMB  = LSYM
      NACTLB = NACTEL
      MCTYPB = MCTYPE
      LSLMXB = LSOLMX
      NLMSLB = NLMSOL
      NLMN1B = NELMN1
      NLMX1B = NELMX1
      NLMN3B = NELMN3
      NLMX3B = NELMX3
      LROTSB = LROOTS
      NROTSB = NROOTS
      IORTOB = IORTO
      ICI0B  = ICI0
      KDELB  = KDEL
      ICHCKB = ICHECK
      NTITB  = NTIT
      MXMACB = MAXMAC
      MXMICB = MAXMIC
      MAXJTB = MAXJT
      MXCITB = MAXCIT
      MXUITB = MAXUIT
      MXAPMB = MAXAPM
      MXABSB = MAXABS
      ITRLVB = ITRLVL
      ITRFNB = ITRFIN
      JCHSMB = JCHSYM
      JCHRBB = JCHORB
      NROCIB = NROOCI
      ISTCIB = ISTACI
      NFILDB = NFIELD
      MXCMAB = MXCIMA
      ICCNOB = ICICNO
      IMCNOB = IMCCNO
      DOSCFB = DOSCF
      DOMP2B = DOMP2
      DOCNOB = DOCINO
      DOCIB  = DOCI
      DOMCB  = DOMC
      DORSPB = DORSP
      FCVRBB = FCVORB
      LNROTB = LNOROT
      LMORDB = LMOORD
      DRFCKB = DIRFCK
      CORHLB = CORHOL
      CRRLXB = CORRLX
      RSPHPB = RESPHP
      JOLSNB = JOLSEN
      ABIPHB = ABAIPH
      INRSIB = INERSI
      INRSFB = INERSF
      SUPSMB = SUPSYM
      SPINB  = SPIN
      POTNCB = POTNUC
      EPSOLB = EPSOL
      EPSTTB = EPSTAT
      PCMB   = PCM
      RSOLB(1)= RSOL(1)
      RSOLB(2)= RSOL(2)
      RSOLB(3)= RSOL(3)
      TRGRDB = THRGRD
      TRPWFB = THRPWF
      TRCIB  = THRCI
      TRMCB  = THRMC
      TRCGRB = THRCGR
      CMXMOB = CMAXMO
      TROVLB = THROVL
      TRSSYB = THRSSY
      TITLEB(1) = TITLE(1)
      TITLEB(2) = TITLE(2)
      TITLEB(3) = TITLE(3)
      TITLEB(4) = TITLE(4)
      TITLEB(5) = TITLE(5)
      TITLEB(6) = TITLE(6)
      DO 100 I = 1, MAXRTS
         IROOTB(I) = IROOT(I)
 100  CONTINUE
      DO 110 I = 1, MXCORB
         NOROTB(I) = NOROT(I)
         IMORDB(I) = IMOORD(I)
         CENTB(I)  = CENT(I)
         TYPEB(I)  = TYPE(I)
 110  CONTINUE
      DO 130 I = 1, NFLAG
         FLAGB(I) = FLAG(I)
 130  CONTINUE
      TITMLB(1) = TITMOL(1)
      TITMLB(2) = TITMOL(2)
      DO 150 I = 1, MXFELT
         EFILDB(I) = EFIELD(I)
         LFILDB(I) = LFIELD(I)
 150  CONTINUE
C
C        /INFPRI/ values are saved
C
      DO I = 1, NPFLAG
         P4FLGB(I) = P4FLAG(I)
         P6FLGB(I) = P6FLAG(I)
      END DO
      IPRI4B = IPRI4
      IPRI6B = IPRI6
      IPRSRB = IPRSIR
      IPRCNB = IPRCNO
      IPRDIB = IPRDIA
      IPRSGB = IPRSIG
      IPDNSB = IPRDNS
      IPRSLB = IPRSOL
      IPRKPB = IPRKAP
      IPRMLB = IPRMUL
      IPRCXB = IPRCIX
      IPRHFB = IPRRHF
      IPRFCB = IPRFCK
      IPRAVB = IPRAVE
      MPRI4B = MPRI4
      MPRI6B = MPRI6
      MPRSRB = MPRSIR
      RETURN
      END
C  /* Deck prestr */
      SUBROUTINE PRESTR
C
C     ABACUS-variables concerning execution of modules are restored.
C     Called from ABADRV. Necessary because modules should be run
C     in each geometry iteration.
C
#include "implicit.h"
#include "past.h"
#include "pastbk.h"
         PASTWO = PSTWOB
         PASORT = PSORTB
         PASDIP = PSDIPB
         PASONE = PSONEB
         PASRES = PSRESB
         PASREL = PSRELB
         PASDPR = PSDPRB
         PASCRS = PSCRSB
         PASCZR = PSCZRB
         PASCTR = PSCTRB
         PASCRL = PSCRLB
         PASMAG = PSMAGB
         PASAAT = PSAATB
         PASRTR = PSRTRB
         PASLRS = PSLRSB
         PASTRP = PSTRPB
         PASLNR = PSLNRB
         PASEXC = PSEXCB
         RETURN
      END
C  /* Deck psave */
      SUBROUTINE PSAVE
C
C     ABACUS-variables concerning execution of modules are saved.
C     Called from ABADRV.
C
#include "implicit.h"
#include "past.h"
#include "pastbk.h"
         PSTWOB = PASTWO
         PSORTB = PASORT
         PSDIPB = PASDIP
         PSONEB = PASONE
         PSRESB = PASRES
         PSRELB = PASREL
         PSDPRB = PASDPR
         PSCRSB = PASCRS
         PSCZRB = PASCZR
         PSCTRB = PASCTR
         PSCRLB = PASCRL
         PSMAGB = PASMAG
         PSAATB = PASAAT
         PSRTRB = PASRTR
         PSLRSB = PASLRS
         PSTRPB = PASTRP
         PSLNRB = PASLNR
         PSEXCB = PASEXC
         RETURN
      END
C  /* Deck bndchk */
      SUBROUTINE BNDCHK(WORK, LMWORK2, WRKDLM, PROG)
#include "implicit.h"
      DIMENSION WORK(LMWORK2)
      CHARACTER*(*) PROG
#include "priunit.h"
C
C     Check memory traps. Gives error message if any of the programs
C     have been outside the declared memory area.
C
      IF (WORK(1) .NE. WRKDLM .OR.
     &    WORK(LMWORK2) .NE. WRKDLM) THEN
         WRITE (LUPRI,'(//3A)')
     *      ' WARNING, ',PROG,' module has been out of bounds.'
         IF (WORK(1) .NE. WRKDLM) WRITE (LUPRI,'(/A)')
     *      ' WORK(0) has been destroyed'
         IF (WORK(LMWORK2) .NE. WRKDLM) WRITE (LUPRI,'(/A)')
     *      ' WORK(LMWORK+1) has been destroyed'
         CALL QUIT('WARNING, ' // PROG // ' has been out of bounds.')
      END IF
      RETURN
      END
C  /* Deck exeher */
      SUBROUTINE EXEHER(WORK, LMWORK, WRKDLM)

      use qmcmm, only: read_pot_qmnpmm

#include "implicit.h"
#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION WORK(0:LMWORK)
#include "priunit.h"
#include "gnrinf.h"
! cbiher.h : TWOBYTEPACKING
! nuclei.h : NBASIS
#include "cbiher.h"
#include "nuclei.h"
#include "dftcom.h"
#include "pcmlog.h"
#include "qm3.h"
C
C     Run Integral section
C
      CALL QENTER('HERMIT')
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'ini HERMIT')
      NEWPRP = .TRUE.
      NEWMAT = .TRUE.
      GRDONE = .FALSE.
C     ... we need new grid for numerical integrations when new geometry
C
C     ***************************************************************
C     Reading in the interaction parameters for the MM system types:
C     ***************************************************************
C
      IF (QM3 .AND. .NOT.ONLYOV) CALL QM3_RDPOT()
C
C     ***************************************************************
C     Reading in the interaction parameters for the NP and MM
C     subsystems needed by QM/NP/MM embedding scheme
C     ***************************************************************
C
      IF (QMNPMM .AND. .NOT.ONLYOV) THEN
         CALL READ_POT_QMNPMM()
      END IF
C
C     ***************************************************************
C     Read initial geometry, orbital spec., etc. (MOLECULE format).
C     ***************************************************************
C
C     If we use the Douglas-Kroll transformation, we first generate
C     the necessary integrals in the primitive basis.
C     See further notes in DKH_INTF routine.
C
      IF (DKTRAN .AND. RNHERM) THEN
         DKHINT = .TRUE. ! now we are calculating DKH 1-el. integrals
         CALL DKH_INTF(WORK(1),LMWORK)
         CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'DKH_INTF ')
         DKHINT = .FALSE.
         NEWPRP = .TRUE.
      END IF
C
C     Third parameter .TRUE. in call indicates that LUONEL will be written.
C
      CALL READIN(WORK(1),LMWORK,RNHERM)
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'READIN')
C
C     Initialize Hermit variables, also when RNHERM false (e.g. for ABACUS)
C
      CALL HERINI
      TWOBYTEPACKING = TWOBYTEPACKING .OR. NBASIS .GT. 255
      ! we need two bytes to pack indices if NBASIS .gt. 255; we also need
      ! the .OR. test, because TWOBYTEPACKING might have been requested in input
C
C     ***************************************************************
C     Setup of MM information.
C     ***************************************************************
C
      IF (QM3 .AND. .NOT.ONLYOV) CALL QM3_SETUP()
      CALL FLSHFO(LUPRI)
C
C     ***************************************************************
C     Now we generate all requested integrals in the contracted basis
C     (if DKTRAN: including transformation of the Douglas-Kroll
C     integrals to contracted basis).
C
      IF (RNHERM) THEN
         CALL TITLER('Starting in Integral Section (HERMIT)',' ',200)
         CALL HERCTL(WORK(1),LMWORK)
C        ... Note: if RUNTWO or RUNERI then HERCTL sets FTRCTL = .TRUE.
C        for transformation control, old AO/MO files cannot be used any more.
C        ***************************************************************
         CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'HERMIT')
         CALL TITLER('End of Integral Section (HERMIT)',' ',200)
         CALL FLSHFO(LUPRI)
      END IF
      CALL QEXIT('HERMIT')
      RETURN
      END
C  /* Deck exesir */
      SUBROUTINE EXESIR(WORK, LMWORK, WRKDLM)
#include "implicit.h"
      DIMENSION WORK(0:LMWORK)
#include "priunit.h"
#include "maxorb.h"
#include "gnrinf.h"
#include "inftra.h"
#include "cbihr2.h"
#include "inftap.h"
#include "huckel.h"
C
C
C     Run SIRIUS
C
      CALL QENTER('SIRIUS')
C
      CALL TITLER('Starting in Wave Function Section (SIRIUS)',' ',200)
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'ini SIRIUS') ! boundary checks on WORK
C
C     SIRINI defines buffer lengths for I/O and the IROW() array
C
      CALL SIRINI ! SIRius INItialization
      WRINDX = .TRUE. ! new sirius call, make sure orbital rotation info is written to file in SETWOP
C
      CALL SIRCTL(ICONV,WORK(1),LMWORK)
      ! ICONV > 0: converged
      !       = 0: not converged
      !       < 0: convergence not requested (e.g. restart from old MCSCF)

      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'SIRIUS') ! boundary checks on WORK
C
C     If not numerical differentiation:
C     We should now have MO's on file, no need for DOHUCKEL anymore
C
CHJ   IF (.NOT.NMWALK) THEN
CHJ: May09: we try to set DOHUCKEL true later if NEWSYM true instead
CHJ         HUCKEL guess is only used in SIRIUS if first call or NEWSYM
         DOHUCKEL = .FALSE.
CHJ   END IF
C
C     Several files left open, have to be closed manually
C
C     With the present version of the code, AOSUPINT cannot be used either
C     in RESPONSE or ABACUS, and for this reason we delete the file here
C
      IF (LUSUPM.LE.0) CALL GPOPEN(LUSUPM,
     &   FNSUPM,'UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      CALL GPCLOSE(LUSUPM,'DELETE')
      IF (LUIT2  .GT. 0) CALL GPCLOSE(LUIT2,'DELETE')
      IF (LUIT3  .GT. 0) CALL GPCLOSE(LUIT3,'DELETE')
      IF (LUIT5  .GT. 0) CALL GPCLOSE(LUIT5,'DELETE')
      IF (LUINTM .GT. 0) CALL GPCLOSE(LUINTM,'KEEP')
      CALL TITLER('End of Wave Function Section (SIRIUS)',' ',200)

      IF (ICONV .EQ. 0 .AND. .NOT.WLKREJ) THEN
         WRITE(LUPRI, '(/A)')
     &      '@ - DALTON aborted because wave function not converged!'
         CALL QUIT('DALTON aborted because wave function not converged')
      END IF
      CALL FLSHFO(LUPRI)
      CALL QEXIT('SIRIUS')
      RETURN
      END
C  /* Deck execc */
      SUBROUTINE EXECC(WORK, LMWORK, WRKDLM)
#include "implicit.h"
      DIMENSION WORK(0:LMWORK)
#include "priunit.h"
#include "gnrinf.h"
#include "maxorb.h"
#include "mxcent.h"
#include "inftra.h"
#include "optinf.h"
#include "nuclei.h"
C
C     Run Coupled Cluster program
C
      CALL QENTER('CC')
C
      CALL TITLER('Starting in Coupled Cluster Section (CC)',' ',200)
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'ini CC')
C
      IF (RNSIRI .AND. .NOT. RNHERM) THEN
         CALL SETHER(IPRUSR,.FALSE.,WORK(1),LMWORK)
         CALL ER2INI
      END IF

!
!     The CC module uses a lot of static allocation for pointer arrays.
!     To minimize the static allocation in Dalton we use a smaller
!     max number of orbitals in CC than the global MXCORB. /hjaaj Sep 2013
!
      IF (NBASIS .GT. MXCORB_CC) THEN
         WRITE (LUPRI,'(//A,I6/A,I6/A)')
     &  ' ERROR Too many orbitals for CC module:',NBASIS,
     &  ' ERROR Limit is MXCORB_CC in maxorb.h :',MXCORB_CC,
     %  ' ERROR To proceed increase MXCORB_CC in maxorb.h and recompile'
         CALL QUIT('Too many orbitals for CC module.'//
     &      'Increase MXCORB_CC in maxorb.h and recompile.')
      END IF

      ! ITRNMR > 0 required for proper restart of CPHF solver
      ISAVITR = ITRNMR
      ITRNMR = 1

      CALL CC_DRV(WORK(1),LMWORK)

      ITRNMR = ISAVITR

      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'CC')

      CALL TITLER('End of Coupled Cluster Section (CC)',' ',200)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('CC')
      RETURN
      END
C  /* Deck exeaba */
      SUBROUTINE EXEABA(WORK, LMWORK, WRKDLM)

      use pelib_interface, only: use_pelib, pelib_ifc_activate,
     &                           pelib_ifc_deactivate

#include "implicit.h"
      DIMENSION WORK(0:LMWORK)
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "inftra.h"
#include "inftap.h"
#include "huckel.h"
#include "gnrinf.h"
#include "abainf.h"
#include "pcmlog.h"

      LOGICAL SAVE_PELIB
C
C     Run ABACUS
C
      IF (SKIPAB) RETURN
      CALL QENTER('ABACUS')
      CALL TITLER('Starting in Static Property Section (ABACUS) -',
     &   ' ',200)
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'ini ABACUS')

C     Test if a requested molecular derivative calculation is implemented
C     NAORDR = max order of analytical geometry derivatives

      SAVE_PELIB = .FALSE.
      IF (USE_PELIB()) THEN
         SAVE_PELIB = .TRUE.
         CALL PELIB_IFC_DEACTIVATE()
      END IF
      CALL FNDANA(NAORDR)
      IF (SAVE_PELIB) THEN
         IF (MOLHES .AND. NAORDR .EQ. 2) THEN
            WRITE(LUPRI,'(//A)')
     &      '@ WARNING. Analytical molecular Hessian is approximate'//
     &      ' because PEQM contributions are neglected.'
         END IF
         CALL PELIB_IFC_ACTIVATE()
      END IF
      IF (MOLHES .AND. NAORDR .LT. 2) THEN
         WRITE(LUPRI,'(//A/A,I5/A)')
     &   'ERROR, analytical calculation of molecular Hessian is '//
     &      'not available for this type of wave function.',
     &   '       Analytical derivative order available:',NAORDR,
     &   '       (You might want to choose numerical Hessian if your '//
     &      'molecule is not too big.)'
         CALL QUIT(
     &   'ERROR, analytical calculation of molecular Hessian is '//
     &      'not available for this type of wave function')
      END IF
      IF (MOLGRD .AND. NAORDR .LT. 1) THEN
         WRITE(LUPRI,'(//A/A,I5)')
     &   'ERROR, analytical calculation of molecular gradient is '//
     &      'not available for this type of wave function.',
     &   '       Analytical derivative order available:',NAORDR
         CALL QUIT(
     &   'ERROR, analytical calculation of molecular gradient is '//
     &      'not available for this type of wave function')
      END IF

C     Now execute abacus

      CALL ABACTL(WORK(1),LMWORK)
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'ABACUS')
C
C     Several files left open, have to be closed manually.
C
      NEWPRP = .TRUE.
      NEWMAT = .TRUE.
      IF (LUINTM .GT. 0) CALL GPCLOSE(LUINTM,'KEEP')
C
      CALL TITLER('End of Static Property Section (ABACUS) -',
     &   ' ',200)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('ABACUS')
      RETURN
      END
C  /* Deck dalton_exedrv */
      SUBROUTINE DALTON_EXEDRV(LMWORK, WRKDLM)

#ifdef BUILD_OPENRSP
          use openrsp_input, only: read_openrsp_input
#endif
C
C     Top driver for execution of Dalton.
C     Called by PROGRAM DALTON after reading environment
C     variables.
C
#ifdef HAS_PCMSOLVER
      use pcm_config, only: pcm_configuration, pcm_cfg, pcm_finalize
#endif
      use pelib_interface, only: use_pelib, pelib_ifc_finalize,
     &                           pelib_ifc_deactivate,
     &                           pelib_ifc_dolf, pelib_ifc_lf

#include "implicit.h"
#include "dummy.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0)
C
#include "gnrinf.h"
#include "nuclei.h"
#include "infinp.h"
#include "inftap.h"
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include "infpar.h"
#include "inftra.h"
#include "abainf.h"
#include "molinp.h"
#include "cbiwlk.h"
#include "cbirea.h"
#include "cbisol.h"
#include "pcmlog.h"
Cef begin
#include "incore.h"
Cef end
C
Cholesky
#include "ccdeco.h"
#include "center.h"
Cholesky
C
!
! Programmers please note:
!   WORK(0) and WORK(LMOWRK+1) are allocated and contain the value
!   WRKDLM (set in main program). This is a simple (not fool-proof) test
!   that Dalton has not been out of bounds. Thus all subroutines called
!   here needs to be called with WORK(1) (and not just WORK), __except__
!   the subroutines EXE* (here: work(0)). In these driver routines for some of the
!   major modules in Dalton, the WRKDLM is checked with BNDCHK (which is
!   in this file). /HJAaJ, Jan. 2011
!
      DIMENSION CORBKP(3*MXCENT)
      LOGICAL NEWWLK, EX
      CHARACTER FILENM*12, WORD*7

#include "trkoor.h"
#include "taymol.h"
      REAL*8, allocatable ::  GRDMOL(:), HESMOL(:,:)

      real(8), allocatable :: work(:)
      logical :: openrsp_input_found

      CALL QENTER('DALTON')
      CALL GETTIM(CSTR,WSTR)

#ifdef BUILD_OPENRSP
      ! radovan: we need to read the **OPENRSP input before we run anything
      call gpopen(lucmd, 'DALTON.INP','OLD',' ','FORMATTED',
     &            idummy, .false.)
!     search **OPENRSP in input
      rewind (lucmd)
      openrsp_input_found = .false.
      use_openrsp = .false.
      do
         read (lucmd, '(a)', end=2) word
         if (word == '**OPENR') then
            call read_openrsp_input('DALTON.INP', lucmd)
            openrsp_input_found = .true.
            use_openrsp = .true.
            exit
         end if
      end do
2     rewind(lucmd)
      call gpclose(lucmd, 'KEEP')
#else
      use_openrsp = .false.
#endif

C
!     allocate memory using f90 utilities
      allocate(work(0:lmwork+1),stat=i)
!     Set memory traps
      work(0)        = wrkdlm
      work(1+lmwork) = wrkdlm

      if(i /= 0)then
         write (luerr,*) mynum,': ALLOCATE command failed to allocate'//
     &        ' the requested memory. Error code:',I
         write (luerr,*) 'Reduce the memory demands and be welcome back'
         call quit('Failed to allocate memory')
      end if

      IPRINT = IPRUSR
      IPRUSR_orig = IPRUSR
C
      WLKREJ = .FALSE.
      IF (OPTNEW) THEN
         IF (QM3)
     &      CALL QUIT('Geometry optimization not implemented for .QM3')
         RNHERM = .TRUE.
         RNABAC = .TRUE.
         DOWALK = .FALSE.
         CALL OPTMIN(WORK(1),LMWORK,WRKDLM)
      ELSE IF (OPTWLK) THEN
         IF (QM3)
     &      CALL QUIT('ERROR: .WALK is not implemented for .QM3')
         RNHERM = .TRUE.
         RNABAC = .TRUE.
         IF (NMWALK) THEN
            CALL TITLER('Commencing numerical differentiation '//
     &         'based on .NMDDRV',' ',200)
            MOLGRD = .FALSE.
         ELSE IF (IWKTYP .EQ. 6) THEN
            CALL TITLER('Commencing numerical differentiation '//
     &         'based on .WALK',' ',200)
            DOWALK = .TRUE.
            MOLGRD = .TRUE.
            MOLHES = .TRUE.
         ELSE
            CALL TITLER(
     &      'Commencing geometry optimization based on .WALK',' ',200)
            DOWALK = .TRUE.
            MOLGRD = .TRUE.
            MOLHES = .TRUE.
         END IF
         START  = .TRUE.
         GEOCNV = .FALSE.
         NEWWLK = .TRUE.
         IPREAD_orig = IPREAD
         IPREAD_reduced = -2
C        ---------------------------------------
C        loop over geometries
C        ---------------------------------------
 10      CONTINUE
            IF (LUCME.GT.0) REWIND (LUCME)
            SOLVNT = .FALSE.
            WRITE (LUPRI,'(//5X,A,I5,A//)')
     &      '@@---------- Geometry walk iteration number:',ITERNR,
     &       ' ----------'
            CALL EXEHER(WORK(0),LMWORK,WRKDLM)
            IPREAD = IPREAD_reduced
C           ... reduce output from READIN in following geometries
            WLKREJ = .FALSE.
            NEWGEO = .TRUE.
Cef begin
            IF (AOSAVE) THEN
               LINTSV = .FALSE.
               LINTMP = .FALSE.
               INITX = .FALSE.
               MSAVE = .TRUE.
               MMCORE = MMWORK
               LMCORE = MMCORE
               ISCORE = 1
               JSCORE = ISCORE
               N_SHL = 1
               I_SHL = 1
               INDX_SHL1 = 0
               INDX_SHL2 = 0
               INDX_SHL3 = 0
               INDX_SHL4 = 0
C     CALL CLEAR_INCOREMEM()
            END IF
Cef end
            CALL EXESIR(WORK(0),LMWORK,WRKDLM)
            IF (DOCCSD) CALL EXECC(WORK(0),LMWORK,WRKDLM)
C
C     If this is a solvent run, update information for abacus
C     Add cavity center if solvent model (921014-kvm/hjaaj)
C     =====================================================
C
            SOLVNT = FLAG(16)
            IF (SOLVNT) THEN
               NEWGEO = .TRUE.
               NUCIND = NUCIND + 1
               NCTOT  = NCTOT + 1
               NUCDEP = NUCDEP + 1
               NATOMS = NATOMS + 1
               NCNTCV = NCTOT
               NCLINE(NUCIND) = 0
               NAMN(NUCIND)       = 'cav '
               NAMEX(3*NUCIND-2)  = 'cav  x'
               NAMEX(3*NUCIND-1)  = 'cav  y'
               NAMEX(3*NUCIND)    = 'cav  z'
               NAMDEP(NUCDEP)     = 'cavity'
               NAMDPX(3*NUCDEP-2) = 'cavity x'
               NAMDPX(3*NUCDEP-1) = 'cavity y'
               NAMDPX(3*NUCDEP  ) = 'cavity z'
               IF (NUCDEP .GT. MXCENT) THEN
                  WRITE (LUPRI,'(//2A,/A,I5)')
     &         ' Too many atomic centers: MXCENT exceeded for',
     &         ' solvent cavity,',' Current limit:',MXCENT
                  CALL QUIT('*** ERROR *** MXCENT exceeded')
               END IF
               CORD(1,NUCIND) = D0
               CORD(2,NUCIND) = D0
               CORD(3,NUCIND) = D0
               ISTBNU(NUCIND) = 7
               CHARGE(NUCIND) = D0
               CALL NUCPRO(WORK(1),LMWORK)
CLF
            ELSE IF(PCM) THEN
               NPCMIN =.TRUE.
            END IF
C
            IF (WLKREJ) THEN
               WRITE (LUPRI,'(/A)')
     &         '@ Geometry step was rejected, must backstep ...'
C
C     We have to backstep, get old geometry and call for an update
C
               REJECT = .TRUE.
               CALL DCOPY(3*MXCENT,CORBKP,1,CORD,1)
               CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',' ',IDUMMY,
     &                     .FALSE.)
               CALL WLKDRV(DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                     WORK(1),LMWORK)
               CALL GPCLOSE(LUSIFC,'KEEP')
               REJECT = .FALSE.
               RDINPC = .FALSE.
               RDMLIN = .TRUE.
C
C     Punch MOLECULE input with updated coordinates to XXXX_mol.inp
C
               CALL PNCMOL(ITERNR,IPRINT)
               IF (LUSIFC .GT. 0) CALL GPCLOSE(LUSIFC,'KEEP')
               GOTO 10
            ELSE
C
C     We take a backup of the present geometry in case of backsteps
C
               CALL DCOPY(3*MXCENT,CORD,1,CORBKP,1)
            END IF
C
C     We have to check if properties are to be calculated in first run
C
C
C           *** Walk to make numerical derivatives. ***
C
            IF (NMWALK) THEN
               CALL TITLER(
     &            'Geometry steps to make numerical derivatives',
     &            ' ',200)
               CALL NMDINI(IPRINT)
               CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',
     &                     IDUMMY,.FALSE.)
               REWIND (LUCMD,IOSTAT=IOS)
 1100          READ (LUCMD,'(A7)',IOSTAT=IOS) WORD
               IF (IOS .NE. 0) THEN
                  WRITE(LUPRI,'(//A)') 'ERROR: **NMDDRV input'//
     &            ' required when .NMDDRV requested'
                  CALL QUIT('ERROR: **NMDDRV input'//
     &            ' required when .NMDDRV requested')
               END IF
               CALL UPCASE(WORD)
               IF (WORD .NE. '**NMDDR') GOTO 1100
               CALL NMDINP(WORD,IPRINT)
               CALL GPCLOSE(LUCMD,'KEEP')
               CALL NUMDRV(WORK(1),LMWORK,IPRINT,WRKDLM)
               GEOCNV = .TRUE.
c#if defined (VAR_MPI)
c               IF (NODTOT .GE. 1 .AND. MYNUM .EQ. 0) THEN
c                  IJOB = 0
c                  CALL MPI_BCAST(IJOB,1,my_MPI_INTEGER,NODEID(MASTER),
c     &                           MPI_COMM_WORLD,IERR)
c                  GOTO 87
c               END IF
c#endif
            ELSE
               WRINDX = .TRUE.
               IF (NEWWLK) THEN
                  CALL WLKINI
                  IF (WALKIN) THEN
C
C     We made some changes to walk default. Need to process these
C
                     CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ',
     &                           'FORMATTED',IDUMMY,.FALSE.)
                     REWIND (LUCMD,IOSTAT=IOS)
 901                 READ (LUCMD,'(A7)') WORD
                     CALL UPCASE(WORD)
                     IF (WORD .NE. '*WALK  ') GOTO 901
                     CALL WLKINP(WORD)
                     CALL GPCLOSE(LUCMD,'KEEP')
                  END IF
               END IF
               IF (ITERNR .EQ. 0) THEN
                  CALL ABAINP('**START',WORK(1),LMWORK)
               ELSE
                  CALL ABAINP('**EACH ',WORK(1),LMWORK)
               END IF
               IF (NEWWLK) THEN
                  NEWWLK = .FALSE.
                  CALL WLKINI
                  IF (WALKIN) THEN
C
C     We made some changes to walk default. Need to process these
C
                     CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ',
     &                           'FORMATTED',IDUMMY,.FALSE.)
                     REWIND (LUCMD,IOSTAT=IOS)
 900                 READ (LUCMD,'(A7)') WORD
                     CALL UPCASE(WORD)
                     IF (WORD .NE. '*WALK  ') GOTO 900
                     CALL WLKINP(WORD)
                     CALL GPCLOSE(LUCMD,'KEEP')
                  END IF
               END IF
C     We may have more DOSYM under **START than under **PROPE
               DO I = 0, 7
                  DOREPW(I) = DOSYM(I + 1) .AND. DOREPW(I)
               END DO
               IF (REUSED .AND. ITERNR .GT. 0) THEN
                  allocate ( GRDMOL(NCOOR), HESMOL(NCOOR,NCOOR) )
                  CALL ABAREAD_TAYMOL(ERGMOL,GRDMOL,HESMOL,NCOOR)
                  REWIND (LUSIFC)
                  CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
                  READ (LUSIFC) POTNUC,EMY,EACTIV,ERGMOL
                  CALL ABAWRIT_TAYMOL(ERGMOL,GRDMOL,HESMOL,NCOOR)
                  deallocate ( GRDMOL, HESMOL )
                  CALL WLKDRV(DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                        WORK(1),LMWORK)
               ELSE
                  CALL EXEABA(WORK(0),LMWORK,WRKDLM)
               END IF
               IF (.NOT. USRIPR .AND. ((IWKTYP .EQ. 3)
     &              .OR. (IWKTYP .EQ. 6))) THEN
                  IF (MOD(ITERNR+2,10) .EQ. 0) THEN
                     IPRUSR = 0
                  ELSE
                     IPRUSR = -2
                  END IF
               END IF
               WRITE (LUPRI,'(/A,I4,A,F14.6,A,F10.6)')
     &            '@  Geometry walk step no.',ITERNR,
     &            ':  Molecular energy =', ERGMOL,
     &            '   Molecular gradient =', GRADML
            END IF
C
C     Punch MOLECULE input with updated coordinates to XXXX_mol.inp
C
            CALL PNCMOL(ITERNR,IPRINT)
C
            START = .FALSE.
            ITERNR = ITERNR + 1
            RDMLIN = .TRUE.
         IF ((.NOT. GEOCNV).AND.(ITERNR .LE. ITERMX)) THEN
            RDINPC = .FALSE.
            CALL GPINQ('RSPVEC','EXIST',EX)
            IF (EX) THEN
               LU = -9982
               CALL GPOPEN(LU,'RSPVEC','OLD',' ','UNFORMATTED',IDUMMY,
     &                     .FALSE.)
               CALL GPCLOSE(LU,'DELETE')
            END IF
            GOTO 10
         END IF
C
C        ---------------------------------------------------------------
C
         IF (GEOCNV) THEN
            IPRUSR = IPRUSR_orig
            IF (IWKTYP .EQ. 6) THEN
              CALL TITLER('Nummerical differentiation finished',' ',200)
            ELSE
              CALL TITLER('Geometry walk has converged',' ',200)
            END IF
C
C     We need to read in starting geometry before doing the vibrational
C     analysis in a ROA calculation, K.Ruud and G.Hangartner, Oct.-96
C
            WRINDX = .TRUE.
            IF (IWKTYP .EQ. 6 .AND. NUMHES) THEN
               NEWGEO = .TRUE.
               CALL ABAINP('**PROPE',WORK(1),LMWORK)
               CALL VIBCTL(WORK(1),LMWORK)
            ELSE IF (IWKTYP .EQ. 6) THEN
               RDINPC = .FALSE.
               SOLVNT = .FALSE.
               IPREAD = IPREAD_orig
               CALL EXEHER(WORK(0),LMWORK,WRKDLM)
               NEWGEO = .TRUE.
               CALL EXESIR(WORK(0),LMWORK,WRKDLM)
               SOLVNT = FLAG(16)
              IF (SOLVNT) THEN
               NEWGEO = .TRUE.
               NUCIND = NUCIND + 1
               NUCDEP = NUCDEP + 1
               NATOMS = NATOMS + 1
               NCNTCV = NUCIND
               NCLINE(NUCIND) = 0
               NAMN(NUCIND)       = 'cav '
               NAMEX(3*NUCIND-2)  = 'cav  x'
               NAMEX(3*NUCIND-1)  = 'cav  y'
               NAMEX(3*NUCIND)    = 'cav  z'
               NAMDEP(NUCDEP)     = 'cavity'
               NAMDPX(3*NUCDEP-2) = 'cavity x'
               NAMDPX(3*NUCDEP-1) = 'cavity y'
               NAMDPX(3*NUCDEP  ) = 'cavity z'
               IF (NUCDEP .GT. MXCENT) THEN
                  WRITE (LUPRI,'(//2A,/A,I5)')
     &         ' Too many atomic centers: MXCENT exceeded for',
     &         ' solvent cavity,',' Current limit:',MXCENT
                  CALL QUIT('*** ERROR *** MXCENT exceeded')
               END IF
               CORD(1,NUCIND) = D0
               CORD(2,NUCIND) = D0
               CORD(3,NUCIND) = D0
               ISTBNU(NUCIND) = 7
               CHARGE(NUCIND) = D0
               CALL NUCPRO(WORK(1),LMWORK)
              END IF
            END IF
            NEWGEO = .TRUE. ! EXESIR resets NEWGEO to .false.
                            ! but we need to copy abacus info to slaves
            CALL ABAINP('**PROPE',WORK(1),LMWORK)
            CALL EXEABA(WORK(0),LMWORK,WRKDLM)
C
C     We also check to see if we have requested a RESPONSE calculations
C
            CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND (LUCMD)
C
C     Let us first see if we can find any **RESPONSE input
C
 133        CONTINUE
            READ (LUCMD,'(A7)',END=134,ERR=134) WORD
            CALL UPCASE(WORD)
            IF (WORD .EQ. '**RESPO') THEN
               CALL TITLER(
     &         'Starting in Dynamic Property Section (RESPONS)',' ',200)
               REWIND (LUCMD)
               RNRESP = .TRUE.
               CALL GPCLOSE(LUCMD,'KEEP')
               CALL RSPDRV(WORK(1),LMWORK)
               CALL TITLER(
     &         'End of Dynamic Property Section (RESPONS)',' ',200)
            ELSE
               GOTO 133
            END IF
 134        CONTINUE

#ifdef BUILD_OPENRSP
C added by Bin Gao, April 2, 2009
C for response theory in atomic orbital basis
            REWIND(LUCMD)
C search for $RESPON in input (legacy driver)
            DO
               READ (LUCMD,'(A7)',END=135) WORD
               IF (WORD .EQ. '$RESPON') THEN
                  USE_OPENRSP = .TRUE.
                  CALL OPENRSP_DALDRV_OLD(WORK(10),LMWORK-9)
                  EXIT
               END IF
            END DO
 135        CONTINUE
 137        CONTINUE
#endif

            IF (LUCMD .GT. 0) CALL GPCLOSE(LUCMD,'KEEP')
         ELSE IF (ITERNR .GT. ITERMX) THEN
            CALL TITLER('WARNING - Geometry not converged',' ',200)
            WRITE (LUPRI,'(/A//)') '@ *** WARNING:  '//
     &         'Maximum number of geometry iterations exceeded.'
         END IF
      ELSE
C     ... not OPTWLK and not OPTNEW, i.e. not geometry optimization
         NEWGEO = .TRUE.
Chj      IF (RNHERM .OR. TESTIN) CALL EXEHER(WORK,LMWORK,WRKDLM)
         CALL EXEHER(WORK(0),LMWORK,WRKDLM)
Chj always call EXEHER to read MOLECULE.INP (not read in abactl any more)
Chj instead "IF (RNHERM) CALL HERCTL" inside EXEHER /April 2009 hjaaj
C
C Request effective dipole integrals using PE library
         IF (USE_PELIB()) THEN
             IF (PELIB_IFC_DOLF()) THEN
!                 WRITE(LUPRI,'(/A)') 'PEQM calculation includes '//
!     &                               'effective external field effects.'
!                 WRITE(LUPRI,'(A)') 'Reference: '//
!     &                        'Phys. Chem. Chem. Phys. 2016, 18, 10070.'
                 CALL PELIB_IFC_LF()
             END IF
         END IF

Cholesky
         IF (CHOINT) CALL EXECHO(WORK(0),LMWORK,WRKDLM)
Cholesky
         IF (RNSIRI) CALL EXESIR(WORK(0),LMWORK,WRKDLM)
Cact
         IF (ACTSEL .AND. DOCCSD) CALL EXEACT(WORK(0),LMWORK,WRKDLM)
Cact
         IF (DOCCSD) CALL EXECC (WORK(0),LMWORK,WRKDLM)
C
C     If this is a solvent run, update information for abacus
C     Add cavity center if solvent model (921014-kvm/hjaaj)
C     =====================================================
C
         SOLVNT = FLAG(16)
         IF (SOLVNT) THEN
            NEWGEO = .TRUE.
            NUCIND = NUCIND + 1
            NUCDEP = NUCDEP + 1
            NATOMS = NATOMS + 1
            NCNTCV = NUCIND
            NCLINE(NUCIND) = 0
            NAMN(NUCIND)       = 'cav '
            NAMEX(3*NUCIND-2)  = 'cav  x'
            NAMEX(3*NUCIND-1)  = 'cav  y'
            NAMEX(3*NUCIND)    = 'cav  z'
            NAMDEP(NUCDEP)     = 'cavity'
            NAMDPX(3*NUCDEP-2) = 'cavity x'
            NAMDPX(3*NUCDEP-1) = 'cavity y'
            NAMDPX(3*NUCDEP  ) = 'cavity z'
            IF (NUCDEP .GT. MXCENT) THEN
               WRITE (LUPRI,'(//2A,/A,I5)')
     &          ' Too many atomic centers: MXCENT exceeded for',
     &          ' solvent cavity,',' Current limit:',MXCENT
               CALL QUIT('*** ERROR *** MXCENT exceeded')
            END IF
            CORD(1,NUCIND) = D0
            CORD(2,NUCIND) = D0
            CORD(3,NUCIND) = D0
            ISTBNU(NUCIND) = 7
            CHARGE(NUCIND) = D0
            CALL NUCPRO(WORK(1),LMWORK)
         END IF
         WRINDX = .TRUE.
         IF (RNRESP) THEN
            CALL TITLER(
     &         'Starting in Dynamic Property Section (RESPONS)',' ',200)
            CALL RSPDRV(WORK(1),LMWORK)
            CALL TITLER(
     &         'End of Dynamic Property Section (RESPONS)',' ',200)
         END IF

#ifdef BUILD_OPENRSP
C added by Bin Gao, April 2, 2009
C for response theory in atomic orbital basis
         IF (LUCMD .LT. 0) THEN
            CALL GPOPEN(LUCMD, 'DALTON.INP', 'OLD', ' ', 'FORMATTED',
     &                  IDUMMY, .FALSE.)
         END IF
         REWIND(LUCMD)
C search for $RESPON in input
         DO
            READ (LUCMD,'(A7)',END=145) WORD
            IF (WORD .EQ. '$RESPON') then
               USE_OPENRSP = .TRUE.
               CALL OPENRSP_DALDRV_OLD(WORK(10),LMWORK-9)
               EXIT
            END IF
         END DO

 145     CONTINUE
         ! radovan: it's really bad to call the driver while reading the input
         !          we should first read all input, then sanity check, then run the code
         if (openrsp_input_found) then
            call openrsp_driver(work(10), lmwork-9)
         end if
 147     if (lucmd > 0) call gpclose(lucmd, 'KEEP')
#endif

         IF (RNABAC) THEN
            CALL ABAINP('**PROPE',WORK(1),LMWORK)
            CALL WLKINI
            IF (WALKIN) THEN
C
C     User has made some changes to walk default. Need to process these
C
               CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',
     &                     IDUMMY,.FALSE.)
               REWIND (LUCMD,IOSTAT=IOS)
 940           READ (LUCMD,'(A7)') WORD
               CALL UPCASE(WORD)
               IF (WORD .NE. '*WALK  ') GOTO 940
               CALL WLKINP(WORD)
               CALL GPCLOSE(LUCMD,'KEEP')
            END IF
            CALL EXEABA(WORK(0),LMWORK,WRKDLM)
         END IF
      END IF
 87   CONTINUE

C
C     We finally punch information for Gamess-US graphic packages
C
      IF (SOLVNT) NUCIND = NUCIND - 1
      KKIND = 1
      KKOL  = KKIND + MXCENT
      KLAST = KKOL  + MXCENT
      LLEFT = LMWORK - KLAST
C
      CALL MOLPLT(WORK(KKIND),WORK(KKOL),TITMOL,WORK(KLAST),LLEFT)
      CALL PLTORB(WORK(KKIND),LMWORK)

C
!     release memory used in PElib and deactivate it
      IF (USE_PELIB()) THEN
        CALL PELIB_IFC_FINALIZE()
        CALL PELIB_IFC_DEACTIVATE()
      END IF

      CALL GETTIM(CEND,WEND)
      CTOT = CEND - CSTR
      WTOT = WEND - WSTR
C
! Finalize pcm_scf if an EXTPCM calculation (PCM using PCMSolver
! library) was requested.
#ifdef HAS_PCMSOLVER
      if (pcm_cfg%do_pcm) then
        call pcm_finalize
        write(lupri, '(//A/)') 'PCMSolver interface correctly finalized'
      end if
#endif

      CALL TIMTXT(' Total CPU  time used in DALTON:',CTOT,LUPRI)
      CALL TIMTXT(' Total wall time used in DALTON:',WTOT,LUPRI)
C
C     Stamp date and time and hostname to output
C
      CALL TSTAMP(' ',LUPRI)

!     release dynamically allocated memory
      deallocate(work)
C
      CALL QEXIT('DALTON')
      RETURN
      END
C  /* Deck gnrlin */
      SUBROUTINE DALTON_GNRLINP
C
C     GENERAL input for Dalton
C
#ifdef HAS_PCMSOLVER
      use pcm_config
      use pcm_write
#endif
      use qmcmm_input, only: qmnpinp
      use pelib_interface, only: use_pelib, pelib_ifc_activate,
     &                           pelib_ifc_input_reader
      use qfitlib_interface, only : qfitlib_ifc_input_reader
      use fde_mod, only : fde_input_init, fde_dalton_input

#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "priunit.h"
#include "gnrinf.h"
C
#include "molde.h"
#include "abainf.h"
#include "exeinf.h"
#include "inftra.h"
#include "huckel.h"
#include "orgcom.h"
#include "clsfmm.h"
#include "cbieri.h"
#include "veclen.h"
#include "numder.h"
#include "pcmdef.h"
#include "pcm.h"
#include "pcmlog.h"
CSONIA
#include "grdccpt.h"
#include "infpar.h"
C
Cholesky
#include "ccdeco.h"
#include "center.h"
C
      PARAMETER (NDIR = 15, NTABLE = 36, D0 = 0.D0)
#if defined (VAR_VECTOR)
      PARAMETER (IVECDF = 128)
#else
      PARAMETER (IVECDF = 1)
#endif
C
      LOGICAL ALLOPT

#ifdef HAS_PCMSOLVER
      logical extpcm_input_provided
#endif

C
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     *          REWORD*12, RWORD*6, WORD1*7
      DATA TABDIR/'*END OF', '*PARALL', '*WALK  ', '*OPTIMI', 'XXXXXXX',
     &            '*MOLBAS', '*PCM   ', '*QM3   ', '*QMMM  ', '*PEQM  ',
     &            '*READIN', '*QMNPMM', '*QFIT  ', '*EXTPCM', '*FDE   '/
      DATA TABLE /'.PRINT ', '.ITERAT', '.PRIERR',
     &            '.INTEGR', '.WAVE F', '.PROPER', !  6
     &            '.INPTES', '.PARALL', '.DIRECT',
     &            '.WALK  ', '.MAX IT', '.RESPON', ! 12
     &            '.RUN PR', '.RUN RE', '.RUN WA',
     &            '.PRESOR', '.WESTA ', '.CHOLES', ! 18
     &            '.EXTPCM', '.FCKTRA', '.2TRATY',
     &            'xXXXXXX', '.OPTIMI', '.TOTSYM', ! 24
     &            '.VECLEN', '.RUNERI', '.NEWTRA',
     &            '.RUN AL', '.NMDDRV', '.PARNMD', ! 30
     &            '.THRRED', '.DOUGLA', '.PEQM  ',
     &            '.QFIT  ', '.LSLIB ', '.FDE   '/ ! 36
C
      CALL QENTER('DALTON_GNRLINP')
C
C     Initialize /VECLEN/
C
      IVECLN = IVECDF
C
C     Initialize /ABAINF/
C
      MOLGRD = .FALSE.
      MOLHES = .FALSE.
      HELFEY = .FALSE.
      DOWALK = .FALSE.
C
C     Initialize /MOLDEN/
C
      MOLDEN = .TRUE.
C
C     Initialize /GNRINF/
C
      PANAS  = 0.0D0
      THR_REDFAC = -1.0D0 ! negative number signals "not set" /hjaaj
      ITERNR = 0
      ITERMX = 20
      IPRUSR = 0
      DOFDE  = .FALSE.
      SEGBAS = .TRUE.
      SEGAUX = .TRUE.
      WALKIN = .FALSE.
      HRINPC = .FALSE.
      SRINPC = .FALSE.
      RDINPC = .FALSE.
      SIR_INPPRC = .FALSE.
      RDMLIN = .FALSE.
      TESTIN = .FALSE.
      OPTWLK = .FALSE.
      USRIPR = .FALSE.
      RNHERM = .FALSE.
      RUNERI = .FALSE.
      RNSIRI = .FALSE.
      RNABAC = .FALSE.
      RNRESP = .FALSE.
      DOCCSD = .FALSE.
      WESTA  = .FALSE.
      QMNPMM = .FALSE.

C     FCKTRA_TYPE and NEWTRA are in /INFTRA/ for practical reasons
C     NEWTRA and NEWTRA_USEDRC are in /INFTRA/ for practical reasons
      FCKTRA_TYPE = -1
      NEWTRA = .FALSE.         ! if true, use "new" transformation for Mulliken integrals
      NEWTRA_USEDRC = .FALSE.  ! if NEWTRA and NEWTRA_USEDRC both true, also use "new" transformation for Dirac integrals

      OPTNEW = .FALSE.
      NEWSYM = .FALSE.
      NEWBAS = .TRUE.
      NEWPRP = .TRUE.
      TOTSYM = .FALSE.
      RELCAL = .FALSE.
      NMWALK = .FALSE.
      DKTRAN = .FALSE.
      DKHINT = .FALSE.
      GEOALL = .FALSE.
      NPCMIN = .TRUE.
      QFIT   = .FALSE.
      USE_LSLIB = .FALSE.
      IF (NODTOT .GE. 1) THEN
         PARCAL = .TRUE.
         DIRCAL = .TRUE.
         PARHER = .TRUE. ! in infpar.h - calculate 2-el. integrals in parallel
      ELSE
         PARCAL = .FALSE.
         DIRCAL = .FALSE.
         PARHER = .FALSE.
      END IF
C
C ...............
C     Flags for calculation of long and short range integrals
C     (mostly programmed by Jesper Kielberg Pedersen) /hjaaj
C     For ERFEXP false:
C       g(R12) = g_erf(R12) = erf(CHIVAL*R12) / R12
C       where CHIVAL .ge. 0, and CHIVAL.eq.0 is normal 1/R12
C     For ERFEXP true:
C       g(R12) = g_erfexp(R12)
C              = g_erf(R12) - 2*CHIVAL/SQRT(PI) * EXP(-CHIVAL**2*R12**2/3.0D0)
      ERFEXP = .FALSE.
C
C     IF (CHIVAL .EQ. -1.0D0) THEN calculate 1/R12
C                             ELSE calculate g(R12) listed above
C
      CHIVAL = -1.0D0
      CHI1ST = .FALSE.
C
C     SRINTS: calculate 1/R12 - g(R12)
C
      DOSRIN = .FALSE.
      SRINTS = .FALSE.
C ...............
C
C     Initialize /EXEINF/
C
      FTWOXP = .TRUE.
      FTRCTL = .FALSE.
C     ... FTRCTL true forces AO sort and integral transformation,
C         but we may restart with old integrals. FTRCTL is set
C         true by EXEHER signalling that new AO integrals have
C         been calculated.
      NEWCMO = .TRUE.
      ! force integral transformation in abaset.F / abarspn.F;
      ! NEWCMO is then set to .false. after integral transformation.
      ! Will be reset to .true. each time ABACTL is called (typically: new geometry)
      ITRLVL_LAST = -999 ! Mulliken (Coulomb) MO integrals not available
      LVLDRC_LAST = -999 ! Dirac (exchange) MO integrals not available
C
C     Initialize /PCMLOG/ and /EXTPCM/
C
#ifdef HAS_PCMSOLVER
      extpcm_input_provided = .false.
#endif
      PCM    = .FALSE.
      OUTFLD = .FALSE.
      NEWMAT = .TRUE.
      LOCFLD = .FALSE.
      NONEQ  = .FALSE.
      NEQRSP = .FALSE.
Clf   NPCMIN = DEFINED ABOVE.....
      OLDCEN = .FALSE.
C
C     Initialize /HUCKEL/
C
      DOHUCKEL = .TRUE.
      EWMO     = .TRUE.
!     HUCCNT   = 1.75D0
      HUCCNT   = 0.75D0 ! hjaaj: because S is ignored ...
      CALL IZERO(NHUCAO,8)
      CALL IZERO(IHUCPT,MXSHEL)
C
Cholesky
      CHOINT = .FALSE.
Cholesky
C
C     Initialize /GRDCCPT/
C
      IGRDCCPT = 1
      LGRDCCPT = .FALSE.

C
C
      CALL TITLER('Output from DALTON general input processing','*',111)
C
C     **** Find General input *****
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
      J = 0
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
         J = J + 1
         CALL UPCASE(WORD)
         IF ((WORD .EQ. '**GENER') .OR. (WORD .EQ. '*GENERA') .OR.
     &       (WORD .EQ. '**DALTO') .OR. (WORD .EQ. '*DALTON')) THEN
            GO TO 930
         ELSE
            GO TO 900
         END IF
  910 CONTINUE
         WRITE (LUPRI,'(A,I10,A)')
     &   'End of file on DALTON.INP, no **DALTON input found in',
     &   J,' lines.'
         CALL QUIT('End of file on DALTON.INP, no **DALTON input found')
  920 CONTINUE
         WRITE (LUPRI,'(A,I10,A)')
     &   'Error reading line',J+1,' of DALTON.INP'//
     &   ' before finding any **DALTON input.'
         CALL QUIT('Error reading DALTON.INP, no **DALTON input found')
  930 CONTINUE
      WORD1 = WORD
C
C     ***** Process input for COMMON  /GENINF/  *****
C
C Table(01:35) :
C 01-05: '.PRINT ', '.ITERAT', '.PRIERR', '.INTEGR', '.WAVE F',
C 06-10: '.PROPER', '.INPTES', '.PARALL', '.DIRECT', '.WALK  ',
C 11-15: '.MAX IT', '.RESPON', '.RUN PR', '.RUN RE', '.RUN WA',
C 16-20: '.PRESOR', '.WESTA ', '.CHOLES', '.EXTPCM', '.FCKTRA',
C 21-25: '.2TRATY', 'xXXXXXX', '.OPTIMI', '.TOTSYM', '.VECLEN',
C 26-30: '.RUNERI', '.NEWTRA', '.RUN AL', '.NMDDRV', '.PARNMD',
C 31-35: '.THRRED', '.DOUGLA', '.PEQM  ', '.QFIT  ', '.LSLIB ',
C 36-40: '.FDE   '.
C
      IPRSTAT = -10
  100 READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114,115,116,117,118,119,120,
     &                121,122,123,124,125,126,127,128,129,129,
     &                131,132,133,134,135,136), I
            END IF
   99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
              CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
              GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     &         ' Keyword ',WORD,' not recognized in DALTON_GNRLINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in DALTON_GNRLINP.')
  101    CONTINUE  ! .PRINT
            READ (LUCMD,*) IPRUSR
            USRIPR = .TRUE.
            GO TO 100
  102    CONTINUE  ! .ITERATions
            READ (LUCMD,*) ITERNR
            GO TO 100
  103    CONTINUE  ! .PRIERR
            READ (LUCMD,*) IPRSTAT
            GO TO 100
  104    CONTINUE  ! .INTEGRals
            RNHERM = .TRUE.
            GO TO 100
  105    CONTINUE  ! .WAVE Function
            RNSIRI = .TRUE.
            GO TO 100
  106    CONTINUE  ! .PROPERties
            RNABAC = .TRUE.
            GO TO 100
  107    CONTINUE  ! .INPTESt
            TESTIN = .TRUE.
            GO TO 100
  108    CONTINUE  ! .PARALLel calculation
            PARCAL = .TRUE.
            DIRCAL = .TRUE.
            PARHER = .TRUE.
            GO TO 100
  109    CONTINUE  ! .DIRECT calculation (AO direct, AOTWOINT file is not generated unless needed by a module)
            DIRCAL = .TRUE.
         GO TO 100
  110    CONTINUE  ! .WALK
            OPTWLK = .TRUE.
            DOWALK = .TRUE.
         GO TO 100
  111    CONTINUE  ! .MAX ITerations
            READ (LUCMD,*) ITERMX
         GO TO 100
  112    CONTINUE  ! .RESPONs
            RNRESP = .TRUE.
         GO TO 100
  113    CONTINUE  ! .RUN PRoperties
            RNHERM = .TRUE.
            RNSIRI = .TRUE.
            RNABAC = .TRUE.
         GO TO 100
 114     CONTINUE  ! .RUN REspons
            RNHERM = .TRUE.
            RNSIRI = .TRUE.
            RNRESP = .TRUE.
         GO TO 100
 115     CONTINUE  ! .RUN WAve function
            RNHERM = .TRUE.
            RNSIRI = .TRUE.
         GO TO 100
 116     CONTINUE   ! .PRESORt
            WRITE(LUPRI,'(//A/)')
     &      'INFO: .PRESORT deprecated, use .NEWTRA'
            NEWTRA = .TRUE.
         GO TO 100
 117     CONTINUE   ! .WESTA
            WESTA = .TRUE.
         GO TO 100
 118     CONTINUE   ! .CHOLESky
            RNHERM = .TRUE.
            CHOINT = .TRUE.
            DIRCAL = .TRUE.
         GO TO 100
 119     CONTINUE  ! .EXTPCM
#ifdef HAS_PCMSOLVER
            extpcm_input_provided = .true.
#else
            WRITE(LUPRI,'(//A/)')
     $         'ERROR: EXTPCM requested but PCM Module not compiled'
            CALL QUIT('Recompile including the PCM Module')
#endif
         goto 100
 120     CONTINUE  ! .FCKTRA
            IF (FCKTRA_TYPE .LT. 0) FCKTRA_TYPE = 1
            NEWTRA = .TRUE.
         GO TO 100
 121     CONTINUE  !.2TRATYpe
         ! "Secret option" for debugging/enhancing 2-electron integral transformation
         !  negative: old integral transformation
         !  zero: new integral transformation
         !  1: standard .FCKTRA
         !  2: .FCKTRA not parallel, even when PARHER
         !  3: .FCKTRA not parallel with ONLY_J on AO2_JINT file
         !  4: .FCKTRA not parallel, do not use ONLY_J on AO2_JINT file
            READ(LUCMD, *) FCKTRA_TYPE
            NEWTRA = FCKTRA_TYPE .GE. 0
         GO TO 100
 122     CONTINUE  ! not used
         GO TO 100
 123     CONTINUE  ! .OPTIMIzation of geometry
            OPTNEW = .TRUE.
         GO TO 100
 124     CONTINUE  ! .TOTSYM
            TOTSYM = .TRUE.
         GO TO 100
 125     CONTINUE  ! .VECLEN
            READ (LUCMD,*) IVECLN
         GO TO 100
 126     CONTINUE  ! .RUNERI
            RUNERI = .TRUE.
         GO TO 100
 127     CONTINUE  ! .NEWTRA
            NEWTRA = .TRUE.
            FCKTRA_TYPE = 0
         GO TO 100
 128     CONTINUE ! .RUN ALl four modules
            RNHERM = .TRUE.
            RNSIRI = .TRUE.
            RNABAC = .TRUE.
            RNRESP = .TRUE.
         GO TO 100
 129     CONTINUE ! '.NMDDRV' and '.PARNMD'
            NMWALK = .TRUE.
            NOMOVE = .TRUE.
            OPTWLK = .TRUE.
         GOTO 100
 131     CONTINUE ! .THRRED
            READ (LUCMD,*) THR_REDFAC
         GOTO 100
 132     CONTINUE ! .DOUGLAs-kroll-hess scalar second order
            DKTRAN = .TRUE.
         GOTO 100
 133     CONTINUE ! .PEQM
            CALL PELIB_IFC_ACTIVATE()
         GOTO 100
 134     CONTINUE ! .QFIT
#if defined (BUILD_QFITLIB)
            QFIT = .TRUE.
#else
            WRITE(LUPRI,*) 'ERROR for .QFIT: QFITLIB not enabled.'
            CALL QUIT('ERROR for .QFIT: QFITLIB not enabled.')
#endif
         GOTO 100
 135     CONTINUE ! .LSLIB 
#if defined (BUILD_LSLIB)
            USE_LSLIB  = .TRUE.
#else
            WRITE(LUPRI,*) 'ERROR for .LSLIB : LSLIB not enabled.'
            CALL QUIT('ERROR for .LSLIB : LSLIB not enabled.')
#endif
         GOTO 100
 136     CONTINUE ! .FDE
            DOFDE = .TRUE.
         GOTO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 999
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in DALTON_GNRLINP.')
      END IF
C
C
  999 CONTINUE
      WRITE (LUPRI,'(1X,A)') SEPARATOR
      WRITE (LUPRI,'(3X,A,I5)')  'Overall default print level:',IPRUSR
      IF (IPRSTAT .EQ. -10) IPRSTAT = IPRUSR + 1
      WRITE (LUPRI,'(3X,A,I5/)') 'Print level for DALTON.STAT:',IPRSTAT
!     WRITE (LUERR,'(A,I5)')     ' Print level for DALTON.STAT',IPRSTAT
      WRITE (LUSTAT,'(A,I5)')    ' Print level for DALTON.STAT',IPRSTAT
C     ... default print level for statistics /hjaaj apr 2000
      IF (TESTIN) WRITE (LUPRI,'(/A/)') ' *** Input test run only ***'
      IF (THR_REDFAC .GT. 0.0D0) THEN
         WRITE (LUPRI,'(A,1P,D10.2)')
     &   '@   INFO: All thresholds are multiplied by',THR_REDFAC
         IF (THR_REDFAC .LT. 1.D-10) THEN
            WRITE (LUPRI,'(A)') ' FATAL INPUT ERROR: Illegal value'
            CALL QUIT('Illegal input value for .THRRED')
         END IF
      END IF
      IF (PARCAL .AND. NODTOT .EQ. 0) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(4X,A)') SEPARATOR,'INFO: '//
     &   'Request for parallel calculation (.PARALL) is ignored '//
     &   'because only one CPU process available.',SEPARATOR
         PARCAL = .FALSE.
         PARHER = .FALSE.
      END IF
#if defined (VAR_MPI)
      IF (PARCAL) WRITE (LUPRI,'(4X,A)')'Parallel calculation using MPI'
#endif
#if defined (VAR_PVM)
      IF (PARCAL) WRITE (LUPRI,'(4X,A)')'Parallel calculation using PVM'
#endif
#if ! defined (VAR_MPI) && ! defined (VAR_PVM)
      IF (PARCAL) THEN
         WRITE (LUPRI,'(//4X,A/4X,A)')
     &      'Parallel calculation requested, but this is not a',
     &      'parallel version. Please recompile for MPI'
         CALL QUIT('Parallel calculation requested w/o MPI activated')
!    &      'parallel version. Please recompile for MPI or PVM'
!        CALL QUIT('Parallel calc. requested w/o MPI or PVM activated')
      END IF
#endif
      IF (OPTWLK) THEN
         IF (NMWALK) THEN
            WRITE (LUPRI,'(A/A,I5)')
     &      '@   Numerical derivatives will be calculated.',
     &      '      Number of different geometries:', ITERMX
         ELSE
            WRITE (LUPRI,'(A/6X,A,I5/6X,A,I5)')
     &        '@   Geometry optimization',
     &        'Starting at iteration:', ITERNR,
     &        'Maximum number of steps:', ITERMX
         END IF
      ELSE IF (ITERNR .GT. 0) THEN
         NINFO = NINFO + 1
         WRITE (LUPRI,'(4X,A/4X,A/6X,A/4X,A)') SEPARATOR,
     &   'INFO: Input specification of geometry iteration number '//
     &     'is only valid for ".WALK  ", ".NMDDRV", or ".PARNMD"',
     &   'Geometry iteration number has been reset to 0',SEPARATOR
         ITERNR = 0
      END IF
      IF (TOTSYM .AND. OPTNEW) THEN
         WRITE(LUPRI,'(//A/4X,A/6X,A)') SEPARATOR,
     &     'ERROR: .TOTSYM is incompatible with .OPTIMIZE',
     &     'Use instead the .WALK module together with .TOTSYM'
         CALL QUIT('.TOTSYM and .OPTIMIZE is incompatible options')
      END IF
      IF (TOTSYM) WRITE (LUPRI,'(4X,A)')
     &   'Only totally symmetric part of molecular Hessian calculated'
      IF (DIRCAL) WRITE(LUPRI,'(4X,A)')
     &   'AO-direct calculation (in sections where implemented)'
      IF (DKTRAN) WRITE (LUPRI,'(4X,A)')
     &   'The second order Douglas-Kroll-Hess transformation is applied'
      IF (RNHERM) WRITE (LUPRI,'(4X,A)')
     &   'HERMIT 1- and 2-electron integral sections will be executed'
      IF (RUNERI) THEN
         WRITE (LUPRI,'(4X,A)')
     &     '2-elctron integrals are calculated using ERI '//
     &     'instead of HERMIT (where possible)'
!        RUNERI = RUNERI .AND. .NOT. PARCAL
!        IF (.NOT. RUNERI) THEN
!           NINFO = NINFO + 1
!           WRITE (LUPRI,'(3X,A)')
!    & ' INFO: NO! ERI is only implemented for direct and not parallel'
!        END IF
      END IF
      IF (IVECLN .NE. 1) WRITE (LUPRI,'(4X,A,I5)')
     &     'Vector length used in direct Fock matrix calculations '//
     &     '(assuming vector machine)',IVECLN
      IF (FCKTRA_TYPE .gt. 0) THEN
         WRITE(LUPRI,'(4X,A,I5)')
     &     'Fock-matrix based 2-el integral transformation (2016) used,'
     &     //' type',FCKTRA_TYPE
      ELSE IF (NEWTRA) THEN
         WRITE(LUPRI,'(4X,A)')
     &     '"New" (from 1988!) integral transformation used'
      ELSE
         WRITE(LUPRI,'(4X,A)') '"Old" integral transformation used'//
     &     ' (limited to max 255 basis functions)'
      END IF
      IF (RNSIRI) WRITE (LUPRI,'(4X,A)')
     &     'Wave function sections will be executed (SIRIUS module)'
      IF (RNRESP) WRITE (LUPRI,'(4X,A)')
     &     'Dynamic molecular response properties section will be'//
     &          ' executed (RESPONSE module)'
      IF (RNABAC) WRITE (LUPRI,'(4X,A)')
     &     'Static molecular property section will be executed'//
     &          ' (ABACUS module)'
      IF (NMWALK) WRITE (LUPRI,'(4X,A)')
     &   'Numerical derivatives will be calculated'
      IF (WESTA) WRITE (LUPRI,'(4X,A)')
     &  'Information for WESTA will be calculated and written to files.'
      IF (USE_PELIB()) WRITE(LUPRI,'(4X,A)')
     &   'Environment is modeled using polarizable embedding scheme'//
     &   ' (PE library)'
      IF(QFIT) WRITE(LUPRI,'(4X,A)')
     &   'Potential derived multipole moments will be calculated'//
     &   ' (QFITLIB)'
C
Cholesky
      IF (CHOINT) THEN
        IF (PARCAL) CALL PARQUIT('CHOLESKY')
        WRITE (LUPRI,'(4X,A)')
     &  'Cholesky decomposition-based calculation will be done'
      ENDIF
Cholesky
C
      WRITE (LUPRI,'(1X,A)') SEPARATOR
C
C     Initialize for REAINP and READIN (reading of .mol input file)
      IPREAD_ini = IPRUSR + 1
      CALL REAINI(IPREAD_ini,RELCAL,TSTINP)
C
  200 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         READ (LUCMD,'(A7)') WORD
         CALL UPCASE(WORD)
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GOTO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in DALTON_GNRLINP.')
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in DALTON_GNRLINP, illegal prompt.')
      END IF
C     "*PARALL" input
 2    CONTINUE
#if defined (VAR_MPI) || defined (VAR_PVM)
        CALL PRLINP(WORD)
#else
        WRITE (LUPRI,'(//4X,A/4X,A)')
     &    '"*PARALL" input sprecified, but this is not a parallel',
     &    'version. Please recompile DALTON for MPI or PVM.'
         CALL QUIT('*PARALL input specified w/o MPI or PVM activated')
#endif
      GO TO 200
C     "*WALK  " input
    3 WALKIN = .TRUE.
C
C     We have to wait with the input processing of "*WALK  "
C     until we know more about the molecule
C
 203    READ (LUCMD,'(A7)') WORD
        IF (WORD(1:1) .NE. '*') GOTO 203
      GO TO 200
C     "*OPTIMI" input
    4 CONTINUE
        CALL OPINPU(WORD)
      GO TO 200
C     "XXXXXXX" input
    5 CONTINUE
      GO TO 200
C
C     "*MOLBAS" input
    6 CONTINUE
        CALL REAINP(WORD,RELCAL,TSTINP)
      GO TO 200
C
Clf:  "*PCM   " input processing
C
    7 CONTINUE
C       INPERR is defined as a precaution (not used during PCM input processing)
        INPERR = 0
        ALLOPT = .FALSE.
        CALL PCMINP(WORD,INPERR,ALLOPT)
      GO TO 200
C
C     "*QM3   " input
    8 CONTINUE
        CALL QM3INP(WORD)
      GO TO 200
C
C     "*QMMM  " input
    9 CONTINUE
        CALL QMMMIN(WORD)
      GO TO 200

!     "*PEQM  " input
   10 CONTINUE
        CALL PELIB_IFC_INPUT_READER(WORD)
      GO TO 200
C
C     "*READIN" input ( DEPRECATED !!! Use "*MOLBAS")
   11 CONTINUE
        WRITE(LUPRI,'(//A/)')
     &    ' INFO : "*READIN" is deprecated, use instead "*MOLBAS"'
      GO TO 6
C
C     "*QMNPMM" input
   12 CONTINUE
        CALL QMNPINP(WORD)
      GO TO 200
   13 CONTINUE
        CALL QFITLIB_IFC_INPUT_READER(WORD)
      GOTO 200
C     "*END OF" or "**something"

!     *EXTPCM input
   14 CONTINUE
#ifdef HAS_PCMSOLVER
            call move_to_next_star(word,lucmd)
#else
            WRITE(LUPRI,'(//A/)')
     $           'ERROR: EXTPCM requested but PCM Module not compiled'
            CALL QUIT('Recompile including the PCM Module')
#endif
      GOTO 200
!     *FDE  input
   15 CONTINUE
        call fde_input_init(LUPRI,'DALTON.OUT')
        call fde_dalton_input(word,.true.)
      GOTO 200
C     "*END OF" or "**something"
    1 CONTINUE
C
#ifdef HAS_PCMSOLVER
!     read EXTPCM input if applicable
      if(extpcm_input_provided)then
        call read_menu_input('DALTON.INP', lucmd,
     &                       '*EXTPCM', extpcm_input_provided)
        pcm_cfg%do_pcm = .true.
        if (.not. extpcm_input_provided) then
        call quit('Oops, read_menu_input could not find *EXTPCM input')
        end if
      end if
#endif

#if defined (VAR_MPI) || defined (VAR_PVM)
      CALL PRLINP(WORD)
#endif
      CALL OPINPU(WORD)
C
C     In order to see if there is any change in starting orbitals (i.e.
C     not Huckel guess, we have to search for *ORBITA
C
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
 5000 READ (LUCMD,'(A7)',END=8000) WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. '*ORBITA') THEN
 400     READ (LUCMD,'(A7)') WORD
         CALL UPCASE(WORD)
         IF (WORD .EQ. '.MOSTAR') THEN
            READ (LUCMD,'(A7)') WORD
            CALL UPCASE(WORD)
            IF (INDEX(WORD,'EWMO') .GT. 0) THEN
               EWMO = .TRUE.
               DOHUCKEL = .TRUE.
            ELSE IF (INDEX(WORD,'HUCKEL') .GT. 0
     &          .OR. INDEX(WORD,'EHT')    .GT. 0) THEN
               EWMO = .FALSE.
               DOHUCKEL = .TRUE.
            ELSE
               DOHUCKEL = .FALSE.
            END IF
         ELSE IF (WORD(1:1) .EQ. '*') THEN
            GOTO 8000
         ELSE
            GOTO 400
         END IF
      ELSE
         GOTO 5000
      END IF
 8000 CALL FLSHFO(LUPRI)
C
C     *** Doing just a survey of the wavefunctions we use. ***
C
c      INPERR = 0
c      NUMRUN = 0
c      IRDMO  = 8
c      IREST  = 0
c      NSYM   = 1
c      REWIND (LUCMD,IOSTAT=IOS)
c 8100 READ (LUCMD,'(A7)') WORD
c      IF (.NOT.((WORD.EQ.'**WAVE ').OR.(WORD.EQ.'**WAVEF')
c     &                             .OR.(WORD.EQ.'**SIRIUS'))) GOTO 8100
c      CALL NEWINP(INPERR,NUMRUN,IRDMO,IREST,NSYM,NORB)
C
      CALL GPCLOSE(LUCMD,'KEEP')
      CALL QEXIT('DALTON_GNRLINP')
      RETURN
      END
#if defined (VAR_PARIO)
C  /* Deck pariot */
      SUBROUTINE PARIOT
C
C     Master routine for checking whether we will do parallell I/O by
C     parallelizing over nuclear geometries
C
#include "implicit.h"
      INCLUDE 'mpif.h'
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
      CHARACTER WORD*7
#include "priunit.h"
#include "inftap.h"
#include "molinp.h"
#include "infpar.h"
#include "pario.h"
C
      PARIO = .FALSE.
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
 5000 READ (LUCMD,'(A7)',END=4070) WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. '.PARNMD') THEN
         PARIO = .TRUE.
         GOTO 4070
      ELSE
         GOTO 5000
      END IF
 4070 CONTINUE
C
C     We are going to do parallel I/O. Inform the slaves about this
C
      CALL MPI_BCAST(PARIO,1,MPI_LOGICAL,NODEID(MASTER),MPI_COMM_WORLD,
     &               IERR)
C
C     Read in Dalton input. Save it temporarily MLINE
C
      IF (PARIO) THEN
#if defined (VAR_MPI2)
C
C     The following piece of code should be replaced by an initialization
C     of the counter for the numerical differentiation which will be addressed
C     using RMA operations. It should probably be in a common block
C
         I0 = 1
         LUNMCL = -9056
#if defined (SYS_T3D) || defined (SYS_T90)
         CALL PXFGETENV('WRKDIR',6,WRKDIR,LWRKMM,IERR)
#else
         CALL GETENV('WRKDIR',WRKDIR)
#endif
         LENWRK = 0
         DO 43 I = 1, 60
            IF (WRKDIR(I:I) .EQ. ' ') GO TO 44
            LENWRK = LENWRK + 1
 43      CONTINUE
 44      CONTINUE
         IF (WRKDIR(LENWRK:LENWRK) .NE. '/') THEN
            LENWRK = LENWRK + 1
            WRKDIR(LENWRK:LENWRK) = '/'
         END IF
         CALL GPOPEN(LUNMCL,WRKDIR(1:LENWRK)//'NUMCAL','NEW',' ',
     &               'FORMATTED',IDUMMY,.FALSE.)
         WRITE (LUNMCL,'(I5)') I0
         CALL GPCLOSE(LUNMCL,'KEEP')
#endif
         DO ILINE = 1, KMLINE
            MLINE(ILINE) = '                                        '//
     &           '                                        '
         END DO
         REWIND (LUCMD)
         JLINE = 1
 5030    READ (LUCMD,'(A8)',END=5010) MLINE(JLINE)
         JLINE = JLINE + 1
         GOTO 5030
 5010    CALL MPI_BCAST(JLINE,1,my_MPI_INTEGER,NODEID(MASTER),
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MLINE,80*JLINE,MPI_CHARACTER,NODEID(MASTER),
     &                  MPI_COMM_WORLD,IERR)
         CALL GPCLOSE(LUCMD,'KEEP')
C
C     Now we read and send MOLECULE.INP
C
         DO ILINE = 1, JLINE
            MLINE(ILINE)(1:8) = '        '
         END DO
         CALL GPOPEN(LUMOL,'MOLECULE.INP',' ',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
         JLINE = 1
 6030    READ (LUMOL,'(A80)',END=6010) MLINE(JLINE)
         JLINE = JLINE + 1
         GOTO 6030
 6010    CALL MPI_BCAST(JLINE,1,my_MPI_INTEGER,NODEID(MASTER),
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MLINE,80*JLINE,MPI_CHARACTER,NODEID(MASTER),
     &                  MPI_COMM_WORLD,IERR)
         CALL GPCLOSE(LUMOL,'KEEP')
      ELSE
         REWIND (LUCMD)
      END IF
C
C     We should have sent everything we need know. Return to the
C     main program.
C
      RETURN
      END
C  /* Deck parion */
      SUBROUTINE PARION
C
C     Slave routine for learning whether we are going to do parallell I/O
C     by parallelizing over nuclear geometries
C
#include "implicit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
      INCLUDE 'mpif.h'
#include "priunit.h"
#include "inftap.h"
#include "molinp.h"
#include "infpar.h"
#include "pario.h"
C
C     We start by waiting for the master to tell us whether we will do
C     parallel I/O or not
C
      CALL MPI_BCAST(PARIO,1,MPI_LOGICAL,NODEID(MASTER),MPI_COMM_WORLD,
     &               IERR)
      IF (PARIO) THEN
#if defined (SYS_T3D) || defined (SYS_T90)
         CALL PXFGETENV('WRKDIR',6,WRKDIR,LWRKMM,IERR)
#else
         CALL GETENV('WRKDIR',WRKDIR)
#endif
         LENWRK = 0
         DO 43 I = 1, 60
            IF (WRKDIR(I:I) .EQ. ' ') GO TO 44
            LENWRK = LENWRK + 1
 43      CONTINUE
 44      CONTINUE
         IF (WRKDIR(LENWRK:LENWRK) .NE. '/') THEN
            LENWRK = LENWRK + 1
            WRKDIR(LENWRK:LENWRK) = '/'
         END IF
         CALL GPOPEN(LUCMD,'DALTON.INP','NEW',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
         CALL GPOPEN(LUMOL,'MOLECULE.INP','NEW',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
         CALL MPI_BCAST(JLINE,1,my_MPI_INTEGER,NODEID(MASTER),
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MLINE,80*JLINE,MPI_CHARACTER,NODEID(MASTER),
     &                  MPI_COMM_WORLD,IERR)
         REWIND (LUCMD)
         DO ILINE = 1, JLINE
            WRITE (LUCMD,'(A80)') MLINE(ILINE)
         END DO
         CALL GPCLOSE(LUCMD,'KEEP')
C
C     We have written DALTON.INP. Now we collect and write MOLECULE.INP
C
         CALL MPI_BCAST(JLINE,1,my_MPI_INTEGER,NODEID(MASTER),
     &                  MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(MLINE,80*JLINE,MPI_CHARACTER,NODEID(MASTER),
     &                  MPI_COMM_WORLD,IERR)
         REWIND (LUMOL)
         DO ILINE = 1, JLINE
            WRITE (LUMOL,'(A80)') MLINE(ILINE)
         END DO
         CALL GPCLOSE(LUMOL,'KEEP')
      END IF
      RETURN
      END
#endif /* VAR_PARIO WARNING: DO NOT activate VAR_PARIO without fixing the PARIO code */
C
Cholesky
C  /* Deck execho */
      SUBROUTINE EXECHO(WORK, LMWORK, WRKDLM)
#include "implicit.h"
#include "dummy.h"
      DIMENSION WORK(0:LMWORK)
#include "priunit.h"
#include "gnrinf.h"
#include "dftcom.h"
C
C     Run Cholesky stuff
C
      CALL QENTER('CHOLESKY')
C
      CALL TITLER('Starting in Cholesky Integral Section',
     &  ' ',200)
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'ini CHOLESKY')
C
      CALL CC_CHODRV(WORK(1),LMWORK)
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'CHOLESKY')
C
ctst  CLOSE(24)
      CALL TITLER('End of Cholesky Integral Section',' ',200)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('CHOLESKY')
      RETURN
      END
C
C /* Deck exeact */
C
      SUBROUTINE EXEACT(WORK, LMWORK, WRKDLM)
#include "implicit.h"
#include "dummy.h"
      DIMENSION WORK(0:LMWORK)
#include "priunit.h"
#include "gnrinf.h"
#include "mxcent.h"
#include "maxorb.h"
#include "center.h"
casm
#include "ccorb.h"
C
C
C     Select orbitals and basis fucntions for hybrid CC models.
C
      CALL QENTER('EXEACT')
C
      CALL TITLER('Starting in Choact Section',' ',200)
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'ini CHOACT')
C
      CALL CC_SELACT(WORK(1),LMWORK)
      CALL BNDCHK(WORK(0), LMWORK+2, WRKDLM, 'CHOACT')
C
ctst  CLOSE(24)
      CALL TITLER('End of Choact Section',' ',200)
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('EXEACT')
      RETURN
      END
C --- end of dalgnr.F ---
