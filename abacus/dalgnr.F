C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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
C  /* Deck srestr */
      SUBROUTINE SRESTR
C
C     Several SIRIUS-variables has to be restored, because they are overwritten
C     by ABACUS in SIRSET module. Called from SIRRDI.
C
#include <implicit.h>
#include <maxorb.h>
#include <gnrinf.h>
#include <siripc.h>
#include <infinp.h>
#include <inpbak.h>
#include <infpri.h>
#include <pribak.h>
C
C     Restoring should only be done if SIRIUS has been through its
C     input processing.
C
      IF (.NOT. INPPRC) RETURN
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
         DORHF  = DORHFB
         DOMP2  = DOMP2B
         DOCINO = DOCNOB
         DOCI   = DOCIB
         FCVORB = FCVRBB
      ELSE
         DORHF  = (DORHFB .AND. .NOT. DOMCB)
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
      TITLE(1) = TITLEB(1)
      TITLE(2) = TITLEB(2)
      TITLE(3) = TITLEB(3)
      TITLE(4) = TITLEB(4)
      TITLE(5) = TITLEB(5)
      TITLE(6) = TITLEB(6)
      DO 100 I = 1, MAXRTS
         IROOT(I) = IROOTB(I)
 100  CONTINUE
      DO 110 I = 1, MXCORB
         NOROT(I)  = NOROTB(I)
         IMOORD(I) = IMORDB(I)
         CENT(I)   = CENTB(I)
         TYPE(I)   = TYPEB(I)
 110  CONTINUE
      DO 130 I = 1, NFLAG
         FLAG(I) = FLAGB(I)
 130  CONTINUE
      DO 140 I = 1, 12
         TITMOL(I,1) = TITMLB(I,1)
         TITMOL(I,2) = TITMLB(I,2)
 140  CONTINUE
      DO 150 I = 1, MXFELT
         EFIELD(I) = EFILDB(I)
         LFIELD(I) = LFILDB(I)
 150  CONTINUE
C     
C     /INFPRI/ values are restored
C     
      DO 200 I = 1, NPFLAG
         P4FLAG(I) = P4FLGB(I)
         P6FLAG(I) = P6FLGB(I)
 200  CONTINUE
      IPRI4  = IPRI4B - 2
      IPRI6  = IPRI6B
      IPRSIR = IPRSRB
      IPRCNO = IPRCNB
      IPRDIA = IPRDIB
      IPRSIG = IPRSGB
      IPRDNS = IPDNSB
      IPRMP2 = IPMP2B
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
C     by ABACUS in SIRSET module. Called from SIRRDI.
C
#include <implicit.h>
#include <maxorb.h>
#include <siripc.h>
#include <infinp.h>
#include <inpbak.h>
#include <infpri.h>
#include <pribak.h>
C
C        /INFINP/ values are saved
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
      DORHFB = DORHF
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
      DO 140 I = 1, 12
         TITMLB(I,1) = TITMOL(I,1)
         TITMLB(I,2) = TITMOL(I,2)
 140  CONTINUE
      DO 150 I = 1, MXFELT
         EFILDB(I) = EFIELD(I)
         LFILDB(I) = LFIELD(I)
 150  CONTINUE
C     
C        /INFPRI/ values are saved
C
      DO 200 I = 1, NPFLAG
         P4FLGB(I) = P4FLAG(I)
         P6FLGB(I) = P6FLAG(I)
 200  CONTINUE
      IPRI4B = IPRI4
      IPRI6B = IPRI6
      IPRSRB = IPRSIR
      IPRCNB = IPRCNO
      IPRDIB = IPRDIA
      IPRSGB = IPRSIG
      IPDNSB = IPRDNS
      IPMP2B = IPRMP2
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
C     Called from ABADRV. Necessary because modules should be runned
C     in each iteration.
C
#include <implicit.h>
#include <past.h>
#include <pastbk.h>
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
#include <implicit.h>
#include <past.h>
#include <pastbk.h>
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
      SUBROUTINE BNDCHK(WORK, LMWORK, WRKDLM, PROG)
#include <implicit.h>
      DIMENSION WORK(LMWORK)
      CHARACTER*6 PROG
#include <priunit.h>
C
C     Check memory traps. Gives error message if any of the programs
C     have been outside the declared memory area.
C
      IF (WORK(1) .NE. WRKDLM .OR.
     &    WORK(LMWORK + 2) .NE. WRKDLM) THEN
         WRITE (LUPRI,'(//A,A6,A)')
     *      ' >>> WARNING, ',PROG,' has been out of bounds.'
         IF (WORK(1) .NE. WRKDLM) WRITE (LUPRI,'(/A)')
     *      ' >>> WORK(0) has been destroyed'
         IF (WORK(LMWORK + 2) .NE. WRKDLM) WRITE (LUPRI,'(/A)')
     *      ' >>> WORK(LMWORK+1) has been destroyed'
         CALL QUIT('WARNING, ' // PROG // ' has been out of bounds.')
      END IF
      RETURN
      END
C  /* Deck exeher */
      SUBROUTINE EXEHER(WORK, LMWORK, WRKDLM)
#include <implicit.h>
#include <dummy.h>
      DIMENSION WORK(LMWORK)
#include <priunit.h>
#include <gnrinf.h>
#include <dftcom.h>
C
C     Run Integral section
C
      CALL QENTER('HERMIT')
      NEWPRP = .TRUE.
      GRDONE = .FALSE.
C
      WRITE(LUPRI, '(/A/)') 'Starting in Integral Section -'
C
C     If we use the Douglas-Kroll transformation, we first generate
C     the necessary integrals in the primitive basis.
C     See further notes in DKINTF routine.
C
      IF (DKTRAN) THEN
         CALL DKINTF(WORK(2),LMWORK)
         NEWPRP = .TRUE.
      END IF
C
C     Now we generate all requested integrals in the contracted
C     basis (if DKTRAN: including transformation of the Douglas-
C     Kroll integrals to contracted basis).
C
      CALL HERCTL(WORK(2),LMWORK)
C     ... Note: if RUNTWO or RUNERI then HERCTL sets FTRCTL = .TRUE.
C     for transformation control, old AO/MO files cannot be used
C     any more.
      CALL BNDCHK(WORK, LMWORK, WRKDLM, 'HERMIT')
      WRITE(LUPRI, '(/A/)') '- End of Integral Section'
      CALL FLSHFO(LUPRI)
      CALL QEXIT('HERMIT')
      RETURN
      END
C  /* Deck exesir */
      SUBROUTINE EXESIR(WORK, LMWORK, WRKDLM)
#include <implicit.h>
      DIMENSION WORK(LMWORK)
#include <priunit.h>
#include <maxorb.h>
#include <siripc.h>
#include <inftra.h>
#include <cbihr2.h>
#include <inftap.h>
#include <huckel.h>
#include <gnrinf.h>
C
C     Run SIRIUS
C
      CALL QENTER('SIRIUS')
      WRITE(LUPRI, '(/A/)') 'Starting in Wave Function Section -'
C     The dollar sign in the previous line may not be standard.
C     It causes the End of Record to be skipped.
C
C     SIRINI defines buffer lengths for I/O and IROW()
C
      CALL SIRINI
C
      WRINDX = .TRUE.
      CALL SIRCTL(ICONV,WORK(2),LMWORK)
      CALL BNDCHK(WORK, LMWORK, WRKDLM, 'SIRIUS')
C
C     If not numerical differentiation:
C     We should now have MO's on file, no need for DOHUCKEL anymore
C
      IF (.NOT.NMWALK) THEN
         DOHUCKEL = .FALSE.
      END IF
C
C     Several files left open, have to be closed manually
C
C     With the present version of the code, AOSUPINT cannot be used either
C     in RESPONSE or ABACUS, and for this reason we delete the file here
C     
      IF (LUSUPM.LE.0) CALL GPOPEN(LUSUPM,
     &   FNSUPM,'UNKNOWN',' ','UNFORMATTED',IDUMMY,.FALSE.)
      CALL GPCLOSE(LUSUPM,'DELETE')
      CALL GPCLOSE(LUIT2,'DELETE')
      CALL GPCLOSE(LUIT3,'DELETE')
      CALL GPCLOSE(LUIT5,'DELETE')
      IF (LUINTM .GT. 0) CALL GPCLOSE(LUINTM,'KEEP')
      IF (.NOT. USRSCR) IFTHRS = 10
      WRITE(LUPRI, '(/A/)') '- End of Wave Function Section'
      IF (ICONV .LT. 1 .AND. .NOT.WLKREJ) THEN
         WRITE(LUPRI, '(/A)')
     &      '- DALTON aborted because wave function not converged!'
         CALL QUIT('DALTON aborted because wave function not converged')
      END IF
      CALL FLSHFO(LUPRI)
      CALL QEXIT('SIRIUS')
      RETURN
      END
C  /* Deck execc */
      SUBROUTINE EXECC(WORK, LMWORK, WRKDLM)
#include <implicit.h>
      DIMENSION WORK(LMWORK)
#include <priunit.h>
#include <siripc.h>
#include <inftra.h>
#include <gnrinf.h>
#include <mxcent.h>
#include <optinf.h>
C
C     Run Coupled Cluster program
C
      CALL QENTER('CC')
      WRITE(LUPRI, '(/A/)') 'Starting Coupled Cluster -'
C     The dollar sign in the previous line may not be standard.
C     It causes the End of Record to be skipped.
      IF (RNSIRI .AND. .NOT. RNHERM) THEN
         CALL SETHER(IPRUSR,.FALSE.,WORK(2),LMWORK)
#if defined (VAR_MPI) || defined (VAR_PVM)
         CALL QUIT('No CC in parallel code')
#else
         CALL ER2INI
#endif
      END IF

      ! ITRNMR > 0 required for proper restart of CPHF solver
      ISAVITR = ITRNMR
      ITRNMR = 1

      CALL CC_DRV(WORK(2),LMWORK)

      ITRNMR = ISAVITR 

      CALL BNDCHK(WORK, LMWORK, WRKDLM, 'CC    ')
      WRITE(LUPRI, '(/A/)') '- End of Coupled Cluster'
      CALL FLSHFO(LUPRI)
      CALL QEXIT('CC')
      RETURN
      END
C  /* Deck exeaba */
      SUBROUTINE EXEABA(WORK, LMWORK, WRKDLM)
#include <implicit.h>
      DIMENSION WORK(LMWORK)
#include <priunit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <inftra.h>
#include <inftap.h>
#include <huckel.h>
#include <gnrinf.h>
#include <abainf.h>
C
C     Run ABACUS
C
      IF (SKIPAB) RETURN
      CALL QENTER('ABACUS')
      WRITE(LUPRI, '(/A/)') 'Starting in Static Property Section -'
C
      CALL ABACTL(WORK(2),LMWORK)
      CALL BNDCHK(WORK, LMWORK, WRKDLM, 'ABACUS')
C
C     Several files left open, have to be closed manually.
C
      NEWPRP = .TRUE.
      IF (LUINTM .GT. 0) CALL GPCLOSE(LUINTM,'KEEP')
      WRITE(LUPRI, '(/A/)') '- End of Static Property Section'
      CALL FLSHFO(LUPRI)
      CALL QEXIT('ABACUS')
      RETURN
      END
C  /* Deck exedrv */
      SUBROUTINE EXEDRV(WORK, LMWORK, WRKDLM,MAST1,MYNUM1)
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0)
#include <siripc.h>
#include <nuclei.h>
#include <infinp.h>
#include <inftap.h>
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include <infpar.h>
#include <inftra.h>
#include <gnrinf.h>
#include <abainf.h>
#include <molinp.h>
#include <cbiwlk.h>
#include <cbisol.h>
#include <taymol.h>
      DIMENSION WORK(LMWORK), CORBKP(3*MXCENT)
      LOGICAL NEWWLK, EX
      CHARACTER FILENM*12, WORD*7
      CALL QENTER('DALTON')
      CALL GETTIM(CSTR,WSTR)
      WLKREJ = .FALSE.
      IF (OPTNEW) THEN
         RNABAC = .TRUE.
         DOWALK = .FALSE.
         CALL OPTMIN(WORK,LMWORK,WRKDLM)
      ELSE IF (OPTWLK) THEN
         RNABAC = .TRUE.
         IF (.NOT. NMWALK) THEN
            DOWALK = .TRUE.
            MOLGRD = .TRUE.
            MOLHES = .TRUE.
         ELSE
            MOLGRD = .FALSE.
         END IF
         START  = .TRUE.
         GEOCNV = .FALSE.
         NEWWLK = .TRUE.
         CALL TITLER('Commencing geometry optimization using .WALK',
     &   '*',111)
 10      CONTINUE
            IF (LUCME.GT.0) REWIND (LUCME)
            SOLVNT = .FALSE.
            WRITE (LUPRI,'(//5X,A,I5,A//)')
     &      '@@>>>>>>>>>> Geometry walk iteration number:',ITERNR,
     &       ' <<<<<<<<<<'
            CALL EXEHER(WORK,LMWORK,WRKDLM)
            WLKREJ = .FALSE.
            NEWGEO = .TRUE.
            CALL EXESIR(WORK,LMWORK,WRKDLM)
            IF (DOCCSD) CALL EXECC(WORK,LMWORK,WRKDLM)
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
     &         ' Too many atomic centers: MXCENT exceed in READIN for',
     &         ' solvent cavity,',' Current limit:',MXCENT
                  CALL QUIT('*** ERROR *** MXCENT exceeded in READIN')
               END IF
               CORD(1,NUCIND) = D0
               CORD(2,NUCIND) = D0
               CORD(3,NUCIND) = D0
               ISTBNU(NUCIND) = 7
               CHARGE(NUCIND) = D0
               CALL NUCPRO(WORK(2),LMWORK)
            END IF
C
            IF (WLKREJ) THEN
               WRITE (LUPRI,'(/A)')
     &         '@ >>> Geometry step was rejected, we backstep ... >>>'
C
C     We have to backstep, get old geometry and call for an update
C
               REJECT = .TRUE.
               CALL DCOPY(3*MXCENT,CORBKP,1,CORD,1)
               CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',' ',IDUMMY,
     &                     .FALSE.)
               CALL WLKDRV(DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                     WORK(2),LMWORK)
               CALL GPCLOSE(LUSIFC,'KEEP')
               REJECT = .FALSE.
               RDINPC = .FALSE.
               RDMLIN = .TRUE.
C
C     Punch MOLECULE input with updated coordinates to MOLECULE.INP
C
               CALL GPOPEN(LUMOL,'MOLECULE.INP','OLD',' ','FORMATTED',
     &                     IDUMMY,.FALSE.)
               DO 230 IMLINE = 1,NMLINE
                  WRITE (LUMOL,'(A)') MLINE(IMLINE)
 230           CONTINUE
               CALL GPCLOSE(LUMOL,'KEEP')
C
C     Punch MOLECULE input with updated coordinates to XXX_mol.inp
C
               FILENM = 'XXXX_mol.inp'
               WRITE(FILENM(1:4),'(I4)') ITERNR + 1
               IF (ITERNR .LT. 999) WRITE(FILENM(1:1),'(A1)') '0'
               IF (ITERNR .LT. 99) WRITE(FILENM(2:2),'(A1)') '0'
               IF (ITERNR .LT. 9) WRITE(FILENM(3:3),'(A1)') '0'
               CALL GPOPEN(LUMOL,FILENM,' ',' ','FORMATTED',IDUMMY,
     &                     .FALSE.)
               DO 330 IMLINE = 1,NMLINE
                  WRITE (LUMOL,'(A)') MLINE(IMLINE)
 330           CONTINUE
               CALL GPCLOSE(LUMOL,'KEEP')
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
               CALL NMDINI(IPRINT)
               CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',
     &                     IDUMMY,.FALSE.)
               REWIND (LUCMD,IOSTAT=IOS)
 1100          READ (LUCMD,'(A7)') WORD
               IF (WORD .NE. '**NMDDR') GOTO 1100
               CALL NMDINP(WORD,IPRINT)
               CALL GPCLOSE(LUCMD,'KEEP')
               CALL NUMDRV(WORK,LMWORK,0,IPRINT,WRKDLM)
               GEOCNV = .TRUE.
c#if defined (VAR_MPI)
c               IF (NODTOT .GT. 1 .AND. MYNUM .EQ. 0) THEN
c                  IJOB = 0
c                  CALL MPI_BCAST(IJOB,1,MPI_INTEGER,NODEID(MASTER),
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
                     IF (WORD .NE. '*WALK  ') GOTO 901
                     CALL WLKINP(WORD)
                     CALL GPCLOSE(LUCMD,'KEEP')
                  END IF
               END IF
               IF (ITERNR .EQ. 0) THEN
                  CALL ABAINP('**START',WORK(2),LMWORK)
               ELSE
                  CALL ABAINP('**PROPE',WORK(2),LMWORK)
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
                  REWIND (LUSIFC)
                  CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
                  READ (LUSIFC) POTNUC,EMY,EACTIV,ERGMOL
                  CALL WLKDRV(DUMMY,DUMMY,DUMMY,DUMMY,DUMMY,
     &                        WORK(2),LMWORK)
               ELSE
                  CALL EXEABA(WORK,LMWORK,WRKDLM)
               END IF
               IF (.NOT. USRIPR .AND. ((IWKTYP .EQ. 3) 
     &              .OR. (IWKTYP .EQ. 6))) THEN
                  IF (MOD(ITERNR+2,10) .EQ. 0) THEN
                     IPRUSR = 0
                  ELSE
                     IPRUSR = -2
                  END IF
               END IF
               WRITE (LUPRI,'(/4X,A,F12.6,A,F13.10)')'Energy: ', ERGMOL,
     &                                           '   Gradient:', GRADML
            END IF
C
C     Punch MOLECULE input with updated coordinates to MOLECULE.INP
C
            CALL GPOPEN(LUMOL,'MOLECULE.INP','OLD',' ','FORMATTED',
     &                  IDUMMY,.FALSE.)
            REWIND (LUMOL)
            DO 200 IMLINE = 1,NMLINE
               WRITE (LUMOL,'(A)') MLINE(IMLINE)
 200        CONTINUE
            CALL GPCLOSE(LUMOL,'KEEP')
C
C     Punch MOLECULE input with updated coordinates to XXXX_mol.inp
C
            FILENM = 'XXXX_mol.inp'
            WRITE(FILENM(1:4),'(I4)') ITERNR + 1
            IF (ITERNR .LT. 999) WRITE(FILENM(1:1),'(A1)') '0'
            IF (ITERNR .LT. 99) WRITE(FILENM(2:2),'(A1)') '0'
            IF (ITERNR .LT. 9) WRITE(FILENM(3:3),'(A1)') '0'
            CALL GPOPEN(LUMOL,FILENM,' ',' ','FORMATTED',IDUMMY,
     &                  .FALSE.)
            DO 300 IMLINE = 1,NMLINE
               WRITE (LUMOL,'(A)') MLINE(IMLINE)
 300        CONTINUE
            CALL GPCLOSE(LUMOL,'KEEP')
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
C
C
         IF (GEOCNV) THEN
C            WRITE (LUPRI,'(/A//)')
C     &                ' >>>>>>>>>> Geometry has converged <<<<<<<<<<'
            WRITE (LUPRI,'(//5X,A//)')
     &                ' >>>>>>>>>> Geometry has converged <<<<<<<<<<'
C
C     We need to read in starting geometry before doing the vibrational
C     analysis in a ROA calculation, K.Ruud and G.Hangartner, Oct.-96
C
            WRINDX = .TRUE.
            IF (IWKTYP .EQ. 6 .AND. NUMHES) THEN
               CALL ABAINP('**FINAL',WORK(2),LMWORK)
               CALL VIBCTL(WORK,LMWORK)
            ELSE IF (IWKTYP .EQ. 6) THEN
ckr            IF (IWKTYP .EQ. 6) THEN
               RDINPC = .FALSE.
               SOLVNT = .FALSE.
               CALL EXEHER(WORK,LMWORK,WRKDLM)
               CALL EXESIR(WORK,LMWORK,WRKDLM)
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
     &         ' Too many atomic centers: MXCENT exceed in READIN for',
     &         ' solvent cavity,',' Current limit:',MXCENT
                  CALL QUIT('*** ERROR *** MXCENT exceeded in READIN')
               END IF
               CORD(1,NUCIND) = D0
               CORD(2,NUCIND) = D0
               CORD(3,NUCIND) = D0
               ISTBNU(NUCIND) = 7
               CHARGE(NUCIND) = D0
               CALL NUCPRO(WORK(2),LMWORK)
            END IF
            END IF
            CALL ABAINP('**FINAL',WORK(2),LMWORK)
            CALL EXEABA(WORK,LMWORK,WRKDLM)
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
            IF (WORD .EQ. '**RESPO') THEN
               REWIND (LUCMD)
               RNRESP = .TRUE.
               CALL SETSIR(WORK(2),LMWORK)
               IF (LUSIFC.GT.0) CALL GPCLOSE(LUSIFC,'KEEP')
               CALL GPCLOSE(LUCMD,'KEEP')
               CALL RSPDRV(WORK(2),LMWORK)
            ELSE
               GOTO 133
            END IF
 134        CONTINUE 
            IF (LUCMD .GT. 0) CALL GPCLOSE(LUCMD,'KEEP')
         ELSE IF (ITERNR .GT. ITERMX) THEN
C            WRITE (LUPRI,'(//3X,A//)')
C     &       '*** WARNING:  Maximum number of iterations exceeded ***'
            WRITE (LUPRI,'(//3X,A//)')
     &       '*** WARNING:  Maximum number of iterations exceeded ***'
         END IF
      ELSE
C     ... not OPTWLK and not OPTNEW, i.e. not geometry optimization
         NEWGEO = .TRUE.
         IF (RNHERM) CALL EXEHER(WORK,LMWORK,WRKDLM)
         IF (RNSIRI) CALL EXESIR(WORK,LMWORK,WRKDLM)
         IF (DOCCSD) CALL EXECC (WORK,LMWORK,WRKDLM)
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
     &          ' Too many atomic centers: MXCENT exceed in READIN for',
     &          ' solvent cavity,',' Current limit:',MXCENT
               CALL QUIT('*** ERROR *** MXCENT exceeded in READIN')
            END IF
            CORD(1,NUCIND) = D0
            CORD(2,NUCIND) = D0
            CORD(3,NUCIND) = D0
            ISTBNU(NUCIND) = 7
            CHARGE(NUCIND) = D0
            CALL NUCPRO(WORK(2),LMWORK)
         END IF
         WRINDX = .TRUE.
         IF (RNRESP) THEN
            IF (.NOT. RNSIRI) THEN
               CALL SETSIR(WORK(2),LMWORK)
               IF (LUSIFC.GT.0) CALL GPCLOSE(LUSIFC,'KEEP')
            END IF
            CALL RSPDRV(WORK(2),LMWORK)
         END IF
C
         IF (RNABAC) THEN
            CALL ABAINP('**PROPE',WORK(2),LMWORK)
            CALL WLKINI
            IF (WALKIN) THEN
C
C     We made some changes to walk default. Need to process these
C
               CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',
     &                     IDUMMY,.FALSE.)
               REWIND (LUCMD,IOSTAT=IOS)
 940           READ (LUCMD,'(A7)') WORD
               IF (WORD .NE. '*WALK  ') GOTO 940
               CALL WLKINP(WORD)
               CALL GPCLOSE(LUCMD,'KEEP')
            END IF
            CALL EXEABA(WORK,LMWORK,WRKDLM)
         END IF
      END IF
 87   CONTINUE
C
C     We finally punch information for Gamess-US graphic packages
C
      IKIND = 2
      IKOL  = IKIND + MXCENT
      KLAST = IKOL + MXCENT
      KLEFT = 2 + LMWORK - KLAST
C
      CALL MOLPLT(WORK(IKIND),WORK(IKOL),WORK(KLAST),KLEFT)
      CALL PLTORB(WORK(IKIND),LMWORK)
C
      CALL GETTIM(CEND,WEND)
      CTOT = CEND - CSTR
      WTOT = WEND - WSTR
C
      CALL TIMTXT('>>>> Total CPU  time used in DALTON:',CTOT,LUPRI)
      CALL TIMTXT('>>>> Total wall time used in DALTON:',WTOT,LUPRI)
C
C     Stamp date and time and hostname to output
C
      CALL TSTAMP(' ',LUPRI)
C
      CALL QEXIT('DALTON')
      RETURN
      END
C  /* Deck gnrlin */
      SUBROUTINE GNRLIN
C
C     GENERAL input
C
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <priunit.h>
#include <siripc.h>
#include <gnrinf.h>
#include <molde.h>
#include <abainf.h>
#include <exeinf.h>
#include <inftra.h>
#include <huckel.h>
#include <orgcom.h>
#include <clsfmm.h>
#include <cbieri.h>
#include <veclen.h>
#include <numder.h>
CSONIA
#include <grdccpt.h>
      PARAMETER (NDIR = 7, NTABLE = 32, D0 = 0.D0)
#if defined (VAR_VECTOR)
      PARAMETER (IVECDF = 128)
#else
      PARAMETER (IVECDF = 1)
#endif
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     *          REWORD*12, RWORD*6, WORD1*7
      DATA TABDIR /'*END OF', '*PARALL', '*WALK  ',
     &             '*MINIMI', '*NMDDRV', '*OPTIMI',
     &             '*XXXXXX'/
      DATA TABLE  /'.PRINT ', '.ITERAT', '.PRIERR',
     &             '.INTEGR', '.WAVE F', '.PROPER',
     &             '.INPTES', '.PARALL', '.DIRECT',
     &             '.WALK  ', '.MAX IT', '.RESPON',
     &             '.RUN PR', '.RUN RE', '.RUN WA',
     &             '.PRESOR', '.XXXXXX', '.XXXXXX',
     &             '.XXXXXX', '.XXXXXX', '.XXXXXX',
     &             '.XXXXXX', '.XXXXXX', '.TOTSYM',
     &             '.VECLEN', '.RUNERI', '.XXXXXX',
     &             '.RUN AL', '.NMDDRV', '.PARNMD',
     &             '.OPTIMI', '.DOUGLA'/
C
      CALL QENTER('GNRLIN')
C
C     Initialize /SIRIPC/
C
      INPPRC = .FALSE.
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
      SEGBAS = .TRUE.
      WALKIN = .FALSE.
      HRINPC = .FALSE.
      SRINPC = .FALSE.
      RDINPC = .FALSE.
      RDMLIN = .FALSE.
      TESTIN = .FALSE.
      OPTWLK = .FALSE.
      ITERNR = 0
      ITERMX = 20
      IPRUSR = 0
      USRIPR = .FALSE.
      RNHERM = .FALSE.
      RUNERI = .FALSE.
      RNSIRI = .FALSE.
      RNABAC = .FALSE.
      RNRESP = .FALSE.
      PARCAL = .FALSE.
      DIRCAL = .FALSE.
      DOCCSD = .FALSE.
C     NEWTRA is in /INFTRA/ for practical reasons
      NEWTRA = .FALSE.
      OPTNEW = .FALSE.
      NEWSYM = .FALSE.
      NEWBAS = .TRUE.
      NEWPRP = .TRUE.
      TOTSYM = .FALSE.
      RELCAL = .FALSE.
      NMWALK = .FALSE.
      DKTRAN = .FALSE.
      GEOALL = .FALSE.
      PANAS  = 0.0D0
C
C     Initialize /EXEINF/
C
      FTWOXP = .TRUE.
      FTRCTL = .FALSE.
C     ... FTRCTL true forces AO sort and integral transformation,
C         but we may restart with old integrals. FTRCTL is set
C         true by EXEHER signalling that new AO integrals have
C         been calculated.
C
C     Initialize /HUCKEL/
C
      DOHUCKEL = .TRUE.
      HUCCNT = 2.0D0
      CALL IZERO(NHUCAO,8)
      CALL IZERO(IHUCPT,MXSHEL)
C
C     Initialize /GRDCCPT/
C
      IGRDCCPT = 1
      LGRDCCPT = .FALSE.

C
C
      WRITE (LUPRI,'(/A/)')
     *   ' <<<<<<<<<< OUTPUT FROM GENERAL INPUT PROCESSING >>>>>>>>>>'
C
C     **** Find General input *****
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      IF ((WORD .EQ. '**GENER') .OR. (WORD .EQ. '*GENERA') .OR.
     &    (WORD .EQ. '**DALTO') .OR. (WORD .EQ. '*DALTON')) THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
         CALL QUIT('End of file on LUCMD, no General input found')
  920 CONTINUE
         CALL QUIT('Error reading LUCMD, no General input found')
  930 CONTINUE
      WORD1 = WORD
C
C     ***** Process input for COMMON  /GENINF/  *****
C
      IPRERR = -10
  100 READ (LUCMD, '(A7)') WORD
      PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114,115,116,104,105,106,113,
     &                115,114,123,124,125,126,127,128,129,129,
     &                123,130), I
            END IF
   99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
              CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
              GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     &         ' Keyword ',WORD,' not recognized in GNRLIN.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in GNRLIN.')
  101    CONTINUE
            READ (LUCMD,*) IPRUSR
            USRIPR = .TRUE.
            GO TO 100
  102    CONTINUE
            READ (LUCMD,*) ITERNR
            GO TO 100
  103    CONTINUE
            READ (LUCMD,*) IPRERR
            GO TO 100
  104    CONTINUE
            RNHERM = .TRUE.
            GO TO 100
  105    CONTINUE
            RNSIRI = .TRUE.
            GO TO 100
  106    CONTINUE
            RNABAC = .TRUE.
            GO TO 100
  107    CONTINUE
            TESTIN = .TRUE.
            GO TO 100
  108    CONTINUE
            PARCAL = .TRUE.
            DIRCAL = .TRUE.
            GO TO 100
  109    CONTINUE
            DIRCAL = .TRUE.
         GO TO 100
  110    CONTINUE
            OPTWLK = .TRUE.
            DOWALK = .TRUE.
         GO TO 100
  111    CONTINUE
            READ (LUCMD,*) ITERMX
         GO TO 100
  112    CONTINUE
            RNRESP = .TRUE.
         GO TO 100
  113    CONTINUE
            RNHERM = .TRUE.
            RNSIRI = .TRUE.
            RNABAC = .TRUE.
         GO TO 100
  114    CONTINUE
            RNHERM = .TRUE.
            RNSIRI = .TRUE.
            RNRESP = .TRUE.
         GO TO 100
 115     CONTINUE
            RNHERM = .TRUE.
            RNSIRI = .TRUE.
         GO TO 100
 116     CONTINUE
            NEWTRA = .TRUE.
         GO TO 100
 117     CONTINUE
         GO TO 100
 118     CONTINUE
         GO TO 100
 123     CONTINUE
            OPTNEW = .TRUE.
         GO TO 100
 124     CONTINUE
            TOTSYM = .TRUE.
         GO TO 100
 125     CONTINUE
            READ (LUCMD,*) IVECLN
         GO TO 100
 126     CONTINUE
            RUNERI = .TRUE.
         GO TO 100
 127     CONTINUE
         GO TO 100
 128     CONTINUE
            RNHERM = .TRUE.
            RNSIRI = .TRUE.
            RNABAC = .TRUE.
            RNRESP = .TRUE.
         GO TO 100
 129     CONTINUE
            NMWALK = .TRUE.
            NOMOVE = .TRUE.
            OPTWLK = .TRUE.
         GOTO 100
 130     CONTINUE
            DKTRAN = .TRUE.
         GOTO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 999
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in GNRLIN.')
      END IF
  999 CONTINUE
      WRITE (LUPRI,'(/A,I5/)') ' Default print level:    ',IPRUSR
      IF (IPRERR .EQ. -10)IPRERR = IPRUSR + 1
      WRITE (LUERR,'(/A,I5)') ' Error file print level: ',IPRERR
C     ... default print level for statistics /hjaaj apr 2000
      IF (TESTIN) WRITE (LUPRI,'(A/)') ' *** Input test run only ***'
      IF (OPTWLK) THEN
         WRITE (LUPRI,'(3X,A/5X,A,I5/5X,A,I5)')
     &      ' Geometry optimization',
     &        ' Starting at iteration:', ITERNR,
     &        ' Maximum number of iterations:', ITERMX
      ELSE IF (ITERNR .GT. 0) THEN
         WRITE (LUPRI,'(3X,A/3X,A)')
     &' Specification of iteration number only valid for optimization',
     &' Iteration number has been reset to 0'
         ITERNR = 0
      END IF
      IF (TOTSYM .AND. OPTNEW) THEN
         WRITE(LUPRI,'(4X,A,/A)') 
     &        '.TOTSYM is incompatible with .OPTIMIZE',
     &        'Use instead the .WALK module together with .TOTSYM'
         CALL QUIT('.TOTSYM and .OPTIMIZE is incompatible options')
      END IF
      IF (TOTSYM) WRITE (LUPRI,'(4X,A)') 
     &     'Only totally symmetric part of molecular Hessian calculated'
      IF (DIRCAL .AND. .NOT. PARCAL) WRITE(LUPRI,'(4X,A)')
     &         'Direct calculation'
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
     &      'parallel version. Please recompile for MPI or PVM'
         CALL QUIT('Parallel calc. requested w/o MPI or PVM activated')
      END IF
#endif
      CALL MOHEAD
      IF (DKTRAN) WRITE (LUPRI,'(3X,A)') 
     &   ' The second order Douglas-Kroll transformation is applied'
      IF (RNHERM) WRITE (LUPRI,'(3X,A)')
     &     ' Integral sections will be executed'
      IF (RUNERI) THEN
         WRITE (LUPRI,'(3X,A)')
     &     ' Integrals are calculated using ERI '//
     &     'instead of HERMIT (where possible)'
         RUNERI = RUNERI .AND. .NOT. PARCAL
         IF (.NOT. RUNERI) THEN
            NINFO = NINFO + 1
            WRITE (LUPRI,'(3X,A)')
     & ' INFO: NO! ERI is only implemented for direct and not parallel'
         END IF
      END IF
      IF (IVECLN .NE. 1) WRITE (LUPRI,'(4X,A,I5)')
     &     'Vector length used in direct Fock matrix calculations '//
     &     '(assuming vector machine)',IVECLN
      IF (NEWTRA) THEN
         WRITE(LUPRI,'(4X,A)')
     &     '"New" integral transformation used'
      ELSE IF (.NOT. DIRCAL) THEN
         WRITE(LUPRI,'(4X,A)') '"Old" integral transformation used'//
     &     ' (limited to max 255 basis functions)'
      END IF
      IF (RNSIRI) WRITE (LUPRI,'(3X,A)')
     &     ' Wave function sections will be executed'
      IF (RNRESP) WRITE (LUPRI,'(3X,A)')
     &     ' Dynamic molecular property section will be executed'
      IF (RNABAC) WRITE (LUPRI,'(3X,A)')
     &     ' Static molecular property section will be executed'
      IF (NMWALK) WRITE (LUPRI,'(3X,A)') 
     &   'Numerical derivatives will be calculated'
  200 PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         READ (LUCMD,'(A7)') WORD
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GOTO (1,2,3,4,5,4,7), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in GNRLIN.')
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in GNRLIN, error in prompt.')
      END IF
 2    CONTINUE 
#if defined (VAR_MPI) || defined (VAR_PVM)
        CALL PRLINP(WORD)
#else
         WRITE (LUPRI,'(//4X,A/4X,A)')
     &      'Parallel input module found, but this is not a',
     &      'parallel version. Please recompile for MPI or PVM'
         CALL QUIT('Parallel input module w/o MPI or PVM activated')
#endif
      GO TO 200
    3   WALKIN = .TRUE.
C
C     We have to wait with the input processing until we know 
C     more about the molecule
C
 203    READ (LUCMD,'(A7)') WORD
        IF (WORD(1:1) .NE. '*') GOTO 203
      GO TO 200
    4   CONTINUE
        CALL OPINPU(WORD)
      GO TO 200
    5   CONTINUE
C
C     We have to wait with the input processing until we know 
C     more about the molecule
C       
 204    READ (LUCMD,'(A7)') WORD
        IF (WORD(1:1) .NE. '*') GOTO 204        
      GO TO 200
    6   CONTINUE
      GO TO 200
    7   CONTINUE
      GO TO 200
    1 CONTINUE
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
      IF (WORD .EQ. '*ORBITA') THEN
 400     READ (LUCMD,'(A7)') WORD
         IF (WORD .EQ. '.MOSTAR') THEN
            READ (LUCMD,'(A7)') WORD
            IF (WORD(1:6) .NE. 'HUCKEL' .AND. WORD(2:7) .NE. 'HUCKEL')
     &         DOHUCKEL = .FALSE.
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
c      CALL NEWINP(INPERR,NUMRUN,IRDMO,IREST,NSYM)
C
      CALL GPCLOSE(LUCMD,'KEEP')
      CALL QEXIT('GNRLIN')
      RETURN
      END
#if defined (VAR_MPI)
C  /* Deck pariot */
      SUBROUTINE PARIOT
C
C     Master routine for checking whether we will do parallell I/O by 
C     parallelizing over nuclear geometries
C
#include <implicit.h>
      INCLUDE 'mpif.h'
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
      CHARACTER WORD*7
#include <priunit.h>
#include <inftap.h>
#include <molinp.h>
#include <infpar.h>
#include <pario.h>
C
      PARIO = .FALSE.
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
 5000 READ (LUCMD,'(A7)',END=4070) WORD
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
         DO ILINE = 1, 500
            MLINE(ILINE) = '                                        '//
     &           '                                        '
         END DO
         REWIND (LUCMD)
         JLINE = 1
 5030    READ (LUCMD,'(A8)',END=5010) MLINE(JLINE)
         JLINE = JLINE + 1
         GOTO 5030
 5010    CALL MPI_BCAST(JLINE,1,MPI_INTEGER,NODEID(MASTER),
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
         CALL GPOPEN(LUMOL,'MOLECULE.INP','OLD',' ','FORMATTED',IDUMMY,
     &               .FALSE.)
         JLINE = 1
 6030    READ (LUMOL,'(A80)',END=6010) MLINE(JLINE)
         JLINE = JLINE + 1
         GOTO 6030
 6010    CALL MPI_BCAST(JLINE,1,MPI_INTEGER,NODEID(MASTER),
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
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
      INCLUDE 'mpif.h'
#include <priunit.h>
#include <inftap.h>
#include <molinp.h>
#include <infpar.h>
#include <pario.h>
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
         CALL MPI_BCAST(JLINE,1,MPI_INTEGER,NODEID(MASTER),
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
         CALL MPI_BCAST(JLINE,1,MPI_INTEGER,NODEID(MASTER),
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
#endif
