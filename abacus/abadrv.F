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
C  /* Deck abactl */
      SUBROUTINE ABACTL(WORK,LWORK)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
#include <dummy.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION WORK(LWORK)
C
#include <abainf.h>
#include <cbisol.h>
#include <magone.h>
#include <energy.h>
#include <moldip.h>
#include <dipole.h>
#include <difsec.h>
#include <sigma.h>
#include <cbilnr.h>
#include <cbiexc.h>
#include <nuclei.h>
#include <inftap.h>
#include <taysol.h>
#include <inforb.h>
#include <spnout.h>
#include <infinp.h>
#include <numder.h>
#include <dftcom.h>
#include <expopt.h>
#include <pcmlog.h>
#include <ccinftap.h>
CSPAS : 6/11-2009 AO-SOPPA included
#include <soppinf.h>
CKeinSPASmehr
C
C
      LOGICAL PROPTY, OLDDX, DFTADX
#include <past.h>
C
      CHARACTER LABEL1*8, STHELP*10
C---------------------------------------------------------------------
C
      IF (SKIPAB) RETURN
      CALL QENTER('ABACTL')
      CALL TSTAMP(' ',LUPRI)
      CALL GETTIM(TIMABA,WALABA)
      TSTART = TIMABA
C
C     Tag this as an ABACUS run.
C
      CALL ABASET
C     
C     *************************
C     ***** Input Section *****
C     *************************
C
      TIMSTR = TSTART
      CALL FLSHFO(LUPRI)
      CALL GETTIM(TIMINP,TIMDUM)
      TIMINP = TIMINP - TIMSTR
C
C     *******************************************
C     ***** Tra/Rot Symmetry Initialization *****
C     *******************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      CALL TROINV(WORK,LWORK)
      CALL FLSHFO(LUPRI)
      CALL GETTIM(TIMTRO,TIMDUM)
      TIMTRO = TIMTRO - TIMSTR
C
C     *****************************
C     ***** Nuclear Repulsion *****
C     *****************************
C
      CALL NUCREP(WORK,WORK(MXCOOR*MXCOOR+1),WORK(2*MXCOOR*MXCOOR+1))
      CALL FLSHFO(LUPRI)
C
C     *********************************************************
C     ***** One-Electron Expectation Values and Integrals *****
C     *********************************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (.NOT.PASONE) THEN
         PROPTY = .TRUE.
         EXPGRA = EXPGRD
         CALL ONEINT(WORK,LWORK,PASONE,PROPTY,PCM)
         IF (PASONE) CALL ABAWRITRS
         CALL FLSHFO(LUPRI)
      END IF
      CALL GETTIM(TIMONE,TIMDUM)
      TIMONE = TIMONE - TIMSTR
C
C     ************************************************************************
C     ** Magnetic properties and relativistic correction expectation values **
C     ************************************************************************
C
      IF (SPNSPN) THEN
         KDSO  = 1
         KPSO  = KDSO  + MXCOOR*MXCOOR
         KSD   = KPSO  + MXCOOR*MXCOOR
         KFC   = KSD   + MXCOOR*MXCOOR
         KSDFC = KFC   + MXCOOR*MXCOOR
         KLAST = KSDFC + MXCOOR*MXCOOR
         CALL DZERO(WORK,5*MXCOOR*MXCOOR)
         LWRK  = LWORK - KLAST + 1
      ELSE
         KDSO  = 1
         KPSO  = 1
         KSD   = 1 
         KFC   = 1
         KSDFC = 1
         KLAST = 1
         LWRK  = LWORK
      END IF
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (.NOT. PASMAG) THEN
         CALL MAGINT(WORK(KDSO),WORK(KLAST),LWRK,PASMAG)
         IF (PASMAG) CALL ABAWRITRS
         CALL FLSHFO(LUPRI)
      END IF
      CALL GETTIM(TIMMAG,TIMDUM)
      TIMMAG = TIMMAG - TIMSTR
C
C     Open file for derivative fock matrices
C     ======================================
C
      IF (FCKDDR) CALL GPOPEN(LUDFCK,ABADFK,'UNKNOWN','DIRECT',' ',
     &                        IRAT*NNBASX,OLDDX)
C
C     *******************************************
C     ***** Two-electron expectation values *****
C     *******************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (.NOT. PASTWO) THEN
         CALL TWOEXP(WORK(KLAST),LWRK,PASTWO)
         IF (PASTWO) CALL ABAWRITRS
         CALL FLSHFO(LUPRI)
      END IF
      CALL GETTIM(TIMTEX,TIMDUM)
      TIMTEX = TIMTEX - TIMSTR
C
C     *********************************
C     ***** Fast Multipole Method *****
C     *********************************
C
      IF (.FALSE.) THEN
         CALL FMMDRV(WORK(KLAST),LWRK)
      END IF
C
C     *****************************************************************
C     ***** Sort undifferentiated integrals as Dirac distribution *****
C     *****************************************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
#if defined (VAR_ABADRC)
      IF (.NOT. (PASORT .AND. PASRES .AND. PASDIP)) THEN
         IF (H2MO) CALL ABADRC(WORK(KLAST),LWRK)
         CALL FLSHFO(LUPRI)
      END IF
#endif
      CALL GETTIM(TIMDRC,TIMDUM)
      TIMDRC = TIMDRC - TIMSTR
C
C     ***********************************************************
C     ***** Right-hand sides for (MC)SCF response equations *****
C     ***** and (MC)SCF lowest-order reorthonormalization   *****
C     ***********************************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (.NOT. PASORT) THEN
         CALL MCORL(WORK(KLAST),LWRK,PASORT)
         IF (PASORT) CALL ABAWRITRS
         CALL FLSHFO(LUPRI)
      END IF
      CALL GETTIM(TIMORT,TIMDUM)
      TIMORT = TIMORT - TIMSTR
C
C     Close file for derivative fock matrices
C     ======================================
C
      IF (FCKDDR) CALL GPCLOSE(LUDFCK,'DELETE')
C
C     **********************************
C     ***** Dipole moment gradient *****
C     **********************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (.NOT. PASDIP) THEN
         CALL DIPCTL(WORK(KLAST),LWRK,PASDIP)
         IF (PASDIP) CALL ABAWRITRS
         CALL FLSHFO(LUPRI)
      END IF
      CALL GETTIM(TIMDIP,TIMDUM)
      TIMDIP = TIMDIP - TIMSTR
      IF (LUDA1 .GT. 0) THEN
         IF (DFTRUN .AND. MOLHES) THEN
            CALL GPCLOSE(LUDA1,'KEEP')
         ELSE
            CALL GPCLOSE(LUDA1,'DELETE')
         END IF
      END IF
C
C     **************************************
C     ***** Quadrupole moment gradient *****
C     **************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (.NOT. PASQPG) THEN
         CALL QPGCTL(WORK(KLAST),LWRK,PASQPG)
         IF (PASQPG) CALL ABAWRITRS
         CALL FLSHFO(LUPRI)
      END IF
      CALL GETTIM(TIMQPG,TIMDUM)
      TIMQPG = TIMQPG - TIMSTR
C
C     ******************************************
C     ***** Solve MCSCF response equations *****
C     ******************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (.NOT. PASRES) THEN
         CALL RESPON(WORK(KLAST),LWRK,PASRES)
         CALL FLSHFO(LUPRI)
      END IF
      CALL GETTIM(TIMRES,TIMDUM)
      TIMRES = TIMRES - TIMSTR
C
C     *****************************************************************
C     ***** Construct triplet operator right-hand sides and solve *****
C     *****                response equations                     *****
C     *****************************************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (SPNSPN .AND. .NOT. PASTRP) THEN
         CALL TRPDRV(WORK(KLAST),LWRK,PASTRP)
         CALL FLSHFO(LUPRI)
      END IF
      CALL GETTIM(TIMTRP,TIMDUM)
      TIMTRP = TIMTRP - TIMSTR
C
C     ***************************
C     ***** Linear response *****
C     ***************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF ((DOLRES .OR. VCD .OR. SPNSPN .OR. MAGSUS .OR. SHIELD .OR.
     &    SPINRO .OR. MOLGFA) .AND. .NOT. PASLRS) THEN
         CALL LRSDRV(WORK(KLAST),LWRK,PASLRS)
         CALL FLSHFO(LUPRI)
      END IF
      CALL GETTIM(TIMLRS,TIMDUM)
      TIMLRS = TIMLRS - TIMSTR
cLig >> added the call to CTOABA
C
C     *****************************
C     ***** CTOCD-DZ response *****
C     *****************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (CTOCD) THEN
        CALL CTOABA(WORK(KLAST),LWRK) 
        CALL FLSHFO(LUPRI)  
      END IF
      CALL GETTIM(TIMCTO,TIMDUM)
      TIMCTO = TIMCTO - TIMSTR   
cLig <<
C
C     **************************************
C     ***** Relaxation terms (singlet) *****
C     **************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (.NOT. PASREL) THEN
      IF (MOLHES .OR. SHIELD .OR. MAGSUS .OR. VCD .OR.
     &    DIPDER .OR. POLAR .OR. QPGRAD .OR. (SPNSPN .AND. DOPSO) .OR.
     &    ECD .OR. VROA .OR. MOLGFA .OR. SPINRO .OR. OPTROT) THEN
         CALL RELAX(WORK(KPSO),WORK(KSD),WORK(KFC),WORK(KSDFC),
     &              WORK(KLAST),LWRK,PASREL,.FALSE.)
         IF (PASREL) CALL ABAWRITRS
         CALL FLSHFO(LUPRI)
      END IF
      END IF
C
C     **************************************
C     ***** Relaxation terms (triplet) *****
C     **************************************
C
      IF (.NOT. PASRTR) THEN
      IF (SPNSPN .AND. (DOSD .OR. DOFC .OR. DOSDFC) ) THEN
         CALL RELAX(WORK(KPSO),WORK(KSD),WORK(KFC),WORK(KSDFC),
     &              WORK(KLAST),LWRK,PASRTR,.TRUE.)
         IF (PASRTR) CALL ABAWRITRS
         CALL FLSHFO(LUPRI)
      END IF
      END IF
      CALL GETTIM(TIMREL,TIMDUM)
      TIMREL = TIMREL - TIMSTR
C
C     **************************************
C     ***** Atomic Axial Tensors (AAT) *****
C     **************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (.NOT. PASAAT) THEN
         CALL AATDRV(WORK(KLAST),LWRK,PASAAT)
         CALL FLSHFO(LUPRI)
         IF (PASAAT) CALL ABAWRITRS
      END IF
      CALL GETTIM(TIMAAT,TIMDUM)
      TIMAAT = TIMAAT - TIMSTR
C
C     ****************************************
C     **** Do linear response calculation ****
C     ****************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (ABALNR) THEN
C
C     Added an extra dimension to keep complex numbers in case of
C     absorption, panor+kr-02
C
         KPOLDD = KLAST
         KPOLDQ = KPOLDD + 2*9*MXFR
         KPOLDL = KPOLDQ + 2*27*MXFR
         KPOLDA = KPOLDL + 2*9*MXFR
         KPOLVL = KPOLDA + 2*9*MXFR
         KPOLVV = KPOLVL + 2*9*MXFR
         KFOVIB = KPOLVV + 2*9*MXFR
CSPAS:6/10-08: trying to add mass independent vibrational g-factor
C        KLST1  = KFOVIB + 3*NATOMS*3*NATOMS*NSYM*MXFR
         KLST1  = KFOVIB + (3*NATOMS+3)*(3*NATOMS+3)*NSYM*MXFR
CKeinSPASmehr
         LWRK1  = LWORK - KLST1 + 1
         IF (.NOT. PASLNR) THEN
            CALL LNRABA(WORK(KPOLDD),WORK(KPOLDQ),WORK(KPOLDL),
     &                  WORK(KPOLDA),WORK(KPOLVL),WORK(KPOLVV),
     &                  WORK(KFOVIB),WORK(KLST1),LWRK1,PASLNR)
            IF (PASLNR) CALL ABAWRITRS
            CALL FLSHFO(LUPRI)
         END IF
      ELSE
         KPOLDD = KLAST
         KPOLDQ = KLAST
         KPOLDL = KLAST
         KPOLDA = KLAST
         KPOLVL = KLAST
         KPOLVV = KLAST
         KFOVIB = KLAST
         KLST1  = KLAST
         LWRK1  = LWRK
      END IF
      CALL GETTIM(TIMLNR,TIMDUM)
      TIMLNR = TIMLNR - TIMSTR
C
C     ******************************************
C     ***** Electronic excitation energies *****
C     ******************************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
CSPAS : 6/11-2009 AO-SOPPA included
      NMODEL = 1
      IF (AORPA .OR. DCRPA .OR. AOSOP .OR. AOSOC) NMODEL = 4
CKeinSPASmehr
      IF (DOEXCI .AND. .NOT. PASEXC) THEN
         KTRLEN = KLST1
CSPAS : 6/11-2009 AO-SOPPA included
C         KTRVEL = KTRLEN + 3*NSYM*MXNEXI
C         KTQLEN = KTRVEL + 3*NSYM*MXNEXI
C         KTQVEL = KTQLEN + 3*3*NSYM*MXNEXI
C         KTRLON = KTQVEL + 3*3*NSYM*MXNEXI
C         KTRMAG = KTRLON + 3*NSYM*MXNEXI
C         KBSRLO = KTRMAG + 3*NSYM*MXNEXI
C         KEXENG = KBSRLO + 3*NSYM*MXNEXI
C         KFONAC = KEXENG + NSYM*MXNEXI
         KTRVEL = KTRLEN + 3*NSYM*MXNEXI*NMODEL
         KTQLEN = KTRVEL + 3*NSYM*MXNEXI*NMODEL
         KTQVEL = KTQLEN + 3*3*NSYM*MXNEXI
CSPAS:23/5-11: second and third moment sum rules
C        KTQVEL = KTQLEN + 3*3*NSYM*MXNEXI
         KTQVEL = KTQLEN + 3*3*NSYM*MXNEXI*NMODEL
         KTTLEN = KTQVEL + 3*3*NSYM*MXNEXI*NMODEL
C        KTRLON = KTQVEL + 3*3*NSYM*MXNEXI
         KTRLON = KTTLEN + 10*NSYM*MXNEXI*NMODEL
CKeinSPASmehr
         KTRMAG = KTRLON + 3*NSYM*MXNEXI*NMODEL
         KBSRLO = KTRMAG + 3*NSYM*MXNEXI*NMODEL
         KEXENG = KBSRLO + 3*NSYM*MXNEXI*NMODEL
         KSECMT = KEXENG + NSYM*MXNEXI*NMODEL
         KFONAC = KSECMT + 3*NSYM*MXNEXI
CKeinSPASmehr
         KFONA2 = KFONAC + 3*NUCDEP*NSYM*MXNEXI
         KRMLEN = KFONA2 + 3*NUCDEP*NSYM*MXNEXI
         KRQLEN = KRMLEN + 3*3*NSYM*MXNEXI
         KRLEN  = KRQLEN + 3*3*NSYM*MXNEXI
         KRMVEL = KRLEN  + 3*3*NSYM*MXNEXI
         KRQVEL = KRMVEL + 3*3*NSYM*MXNEXI
         KRVEL  = KRQVEL + 3*3*NSYM*MXNEXI
         KSLEN  = KRVEL  + 3*3*NSYM*MXNEXI
         KSVEL  = KSLEN  + NSYM*MXNEXI
         KLST2  = KSVEL  + NSYM*MXNEXI
         LWRK2  = LWORK  - KLST2 + 1
         CALL EXCITA(WORK(KTRLEN),WORK(KTRVEL),WORK(KTQLEN),
     &               WORK(KTQVEL),WORK(KTRLON),WORK(KTRMAG),
CSPAS:23/5-11: second and third moment sum rules
C    &               WORK(KBSRLO),WORK(KEXENG),WORK(KFONAC),
C    &               WORK(KFONA2),WORK(KSECMT),
     &               WORK(KBSRLO),WORK(KTTLEN),WORK(KEXENG),
     &               WORK(KFONAC),WORK(KFONA2),WORK(KSECMT),
CKeinSPASmehr
     &               WORK(KLST2),LWRK2,PASEXC)
         IF (PASEXC) CALL ABAWRITRS
         CALL FLSHFO(LUPRI)
      ELSE
         KTRLEN = KLST1
         KTRVEL = KTRLEN
         KTQLEN = KTRVEL
         KTQVEL = KTQLEN
CSPAS:23/5-11: second and third moment sum rules
C        KTRLON = KTQVEL
         KTTLEN = KTQVEL
         KTRLON = KTTLEN
CKeinSPASmehr
         KTRLON = KTQVEL
         KTRMAG = KTRLON
         KBSRLO = KTRMAG
         KEXENG = KBSRLO
CSPAS : 6/11-2009 AO-SOPPA included
         KSECMT = KEXENG
         KFONAC = KSECMT
C         KFONAC = KEXENG
CKeinSPASmehr
         KFONA2 = KFONAC
         KRMLEN = KFONA2
         KRQLEN = KRMLEN
         KRLEN  = KRQLEN
         KRMVEL = KRLEN
         KRQVEL = KRMVEL
         KRVEL  = KRQVEL
         KSLEN  = KRVEL
         KSVEL  = KSLEN
         KLST2  = KSVEL
         LWRK2  = LWRK1
      END IF
      CALL GETTIM(TIMEXC,TIMDUM)
      TIMEXC = TIMEXC - TIMSTR
C
C     ******************************************
C     ***** Quadratic Response Properties  *****
C     ******************************************
C
cs february 1996
cs
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (HYPER .OR. VERDET .OR. MCD) THEN
         CALL AQRDRV(WORK(KLST2),LWRK2)
      END IF
C
C     ****************************************
C     ***** Add contributions calculated *****
C     ****************************************
C
C     Calculated elements of molecular gradient
C
      CALL ZERGRD
      IF (MOLGRD) THEN
C
C     If this is a finite-field molecular gradient calculation, we now add
C     the electric field term, K.Ruud-120397
C
         IF (MOLGRD .AND. NFIELD .GT. 0) THEN
            DO IFIELD = 1, NFIELD
               IF (EFIELD(IFIELD) .NE. D0) THEN
                  IF (LFIELD(IFIELD) .EQ. 'XDIPLEN ') THEN
                     ICOOR = 1
                  ELSE IF (LFIELD(IFIELD) .EQ. 'YDIPLEN ') THEN
                     ICOOR = 2
                  ELSE IF (LFIELD(IFIELD) .EQ. 'ZDIPLEN ') THEN
                     ICOOR = 3
                  ELSE
                     WRITE (LUPRI,'(/,3A,/)') 'Field type ',
     &                    LFIELD(IFIELD),
     &                    ' not implemented for molecular gradient'
                     CALL QUIT('Illegal field type for mag. properties')
                  END IF
                  DO JCOOR = 1, MXCOOR
C
C     Electronic part
C
                     GRADKE(JCOOR) = GRADKE(JCOOR) 
     &                    - EFIELD(IFIELD)*(DDIPE(ICOOR,JCOOR) + 
     &                                      DDIPN(ICOOR,JCOOR))
                  END DO
C
C     Nuclear part. We assume molecule to be in center-of charge.
C     As it will only contribute in the field direction and we cannot
C     have symmetry in that direction(s), it is very simple
C
C                  DO IATOM = 1, NUCIND
C                     JCOOR = 3*(IATOM - 1) + ICOOR
C                        GRADKE(JCOOR) = GRADKE(JCOOR) 
C     &                       + EFIELD(IFIELD)*CHARGE(IATOM)
C                     END DO
C
               END IF
            END DO
         END IF
         CALL ADDGRD(GRADKE)
         CALL ADDGRD(GRADNA)
         CALL ADDGRD(GRADEE)
         CALL ADDGRD(GRADNN)
         CALL ADDGRD(GRADFS)
         CALL ADDGRD(GRADFT)
         IF (SOLVNT.OR.PCM) CALL ADDGRD(GSOLTT)
         IF (SOLVNT.OR.PCM) CALL ADDGRD(GSOLNN)
         IF (PCM) CALL ADDGRD(GSOLCV)
      END IF
C
C     Calculated dipole moment
C
      CALL DP0SUM
C
C     Calculated elements of dipole gradient
C
      CALL DZERO(DIP1,9*NUCDEP)
Clf to be adjusted for PCM IR intensities calculation
      IF (DIPDER) THEN
         CALL DIPADD(DDIPN)
         CALL DIPADD(DDIPE)
         CALL DIPADD(DDIPS)
         CALL DIPADD(DDIPR)
      END IF
C
C     Calculated elements of second moment gradient
C
      CALL DZERO(SEC1,27*NUCDEP)
      IF (QPGRAD) THEN
         CALL SECADD(DSECN)
         CALL SECADD(DSECE)
         CALL SECADD(DSECS)
         CALL SECADD(DSECR)
      END IF
C
C     Calculated elements of shielding tensors
C
      CALL DZERO(SIGMAT,9*NUCDEP)
      IF (SHIELD) THEN
         CALL SIGADD(SIGMAD)
         CALL SIGADD(SIGMAS)
         CALL SIGADD(SIGMAR)
      END IF
C
C     ****************************
C     ***** Tra/rot symmetry *****
C     ****************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      CALL TROINV(WORK(KLST2),LWRK1)
      CALL GETTIM(TIMTR2,TIMDUM)
      TIMTRO = TIMTRO + TIMTR2 - TIMSTR
      CALL FLSHFO(LUPRI)
C
C     *******************************************************
C     ***** Fold in contrbutions from floating orbitals *****
C     *******************************************************
C
      CALL FLTORB(WORK(KLST2),LWRK1)
      CALL FLSHFO(LUPRI)
C
C     **************************
C     ***** Output section *****
C     **************************
C
      KCSTRA = KLST2
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
      LWRK3  = LWORK - KLAST + 1
      CALL ABARESULT(WORK(KPOLDD),WORK(KPOLDL),WORK(KPOLDA),
     &            WORK(KPOLVL),WORK(KPOLVV),WORK(KCSTRA),
     &            WORK(KSCTRA),WORK(KDSO),WORK(KPSO),WORK(KSD),
     &            WORK(KFC),WORK(KSDFC),WORK(KTRLEN),
     &            WORK(KTRVEL),WORK(KTQLEN),WORK(KTQVEL),
CSPAS:23/5-11: second and third moment sum rules
C    &            WORK(KTRMAG),WORK(KTRLON),WORK(KBSRLO),
     &            WORK(KTRMAG),WORK(KTRLON),WORK(KBSRLO),WORK(KTTLEN),
CKeinSPASmehr
     &            WORK(KEXENG),WORK(KFONAC),WORK(KFONA2),WORK(KFOVIB),
     &            WORK(KSECMT),
     &            WORK(KRMLEN),WORK(KRQLEN),WORK(KRLEN),
     &            WORK(KRMVEL),WORK(KRQVEL),WORK(KRVEL),
     &            WORK(KSLEN),WORK(KSVEL),WORK(KLAST),LWRK3)
      IF (IPRDEF .GE. 0)
     &   CALL GEOANA(CORD,.TRUE.,.FALSE.,NBONDS,-1,WORK(KLAST),LWRK3)
      CALL FLSHFO(LUPRI)
C
C     *************************
C     ***** Geometry walk *****
C     *************************
C
      IF (DOWALK.AND..NOT.NWPROP) THEN
         TIMSTR = SECOND()
         IF (SPNSPN) THEN
C           ... save "SPNTOT" in "SPNDSO" for WLKDRV
            DO I = 0, MXCOOR*MXCOOR - 1
               WORK(KDSO + I) = WORK(KDSO + I) + WORK(KPSO + I) 
     &                        + WORK(KSD + I)  + WORK(KFC + I)
     &                        + WORK(KSDFC + I)
            END DO
         END IF
C
         CALL WLKDRV(WORK(KPOLDD),WORK(KPOLDQ),WORK(KPOLDL),
     &               WORK(KPOLDA),WORK(KDSO),WORK(KLAST),LWRK3)
         TIMWLK = SECOND() - TIMSTR
         CALL FLSHFO(LUPRI)
      ELSE
         TIMWLK = D0
      END IF
C     
C     ********************************
C     ***** Vibrational analysis *****
C     ********************************
C
      CALL GETTIM(TIMSTR,TIMDUM)
      IF (VIB) CALL VIBCTL(WORK(KLAST),LWRK3)
      CALL GETTIM(TIMVIB,TIMDUM)
      TIMVIB = TIMVIB - TIMSTR
      CALL FLSHFO(LUPRI)
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C     ******************************
C     ***** End of calculation *****
C     ******************************
C
      CALL GETTIM(TIMALL,TIMDUM)
      TIMALL = TIMALL - TSTART

#if !defined (SYS_T3D)
cLig <> added the TIMCTO
      CALL ABATIMSTA(TIMALL,TIMINP,TIMTRO,TIMONE,TIMMAG,TIMTEX,TIMDRC,
     &     TIMORT,TIMDIP,TIMQPG,TIMRES,TIMREL,TIMTST,TIMVIB,TIMWLK,
     &     TIMLRS,TIMTRP,TIMEXC,TIMAAT,TIMLNR,TIMCTO)
#endif
      CALL GETTIM(TEND,WEND)
      TIMABA = TEND - TIMABA
      WALABA = WEND - WALABA
      CALL TIMTXT('>>>> Total CPU  time used in ABACUS:',TIMABA,LUPRI)
      CALL TIMTXT('>>>> Total wall time used in ABACUS:',WALABA,LUPRI)
      CALL QEXIT('ABACTL')
      RETURN
      END
cs --------------------february 1996-----------------------------------
C  /* Deck abainp */
      SUBROUTINE ABAINP(WRDSRC,WORK,LWORK)
C
C     General input for ABACUS
C     input:  WRDSRC - which ABACUS label to look for in input file
C             (Feb 2003: one of **START, **PROPE, or **FINAL)
C
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <priunit.h>
      PARAMETER (NDIR = 28, NTABLE = 58, D0 = 0.D0)
      LOGICAL   POPULA, NOHESS, FCKEXP, GAGT
      LOGICAL   GRDINP
      CHARACTER WORD*7, PROMPT*1, TABDIR(NDIR)*7, TABLE(NTABLE)*7,
     *          REWORD*12, RWORD*6, RSTLBL(10)*6, WORD1*7, WRDSRC*7,
     *          WORDL*12
      DIMENSION IPOINT(MXCOOR), IDOSYM(8), WORK(LWORK)
C
#include <abainf.h>
#include <gnrinf.h>
#include <energy.h>
#include <orgcom.h>
#include <dorps.h>
#include <nuclei.h>
#include <symmet.h>
#include <inforb.h>
#include <cbisol.h>
#include <past.h>
#include <exeinf.h>
#include <inftap.h>
#include <rspprp.h>
#include <esg.h>
#include <cbiwlk.h>
#include <dftcom.h>
#include <cbilnr.h>
#include <qmmm.h>
C
C Common blocks used for SOPPA :
C use from COMMON /INFDIM/ :  NCONMA,NWOPMA,NVARMA
C use from COMMON /INFSOP/ :  N2P2HS,N2P2HT
C use from COMMON /INFRSP/ :  SOPPA
C     Aug. 97/spas+hjaaj
#include <infdim.h>
#include <infmp2.h>
#include <infsop.h>
#include <infrsp.h>
#include <infinp.h>
CSPAS : 6/11-2009 AO-SOPPA included
#include <soppinf.h>
CKeinSPASmehr
C
      DATA TABDIR /'*END OF', 'xxxxxxx', '*NUCREP',
     &             '*ONEINT', '*TWOEXP', '*GETSGY',
     &             '*EXCITA', '*REORT ', '*RESPON',
     &             '*TROINV', '*FMM   ', '*DIPCTL',
     &             '*RELAX ', '*GEOANA', '*VIBANA',
     &             '*FLOAT ', 'xxxxxxx', '*QPGCTL',
     &             'xXXXXXX', 'xXXXXXX', '*SOPPA ',
     &             '*LINRES', '*EXPECT', '*SPIN-S',
     &             '*QRPROP', '*AAT   ', '*TRPRSP',
     &             '*ABALNR'/
Caug04-hjaaj: disabled currently inactive .NACME
Csep08-spas: .NACME reactivated for Vib_g
C                      1          2          3
      DATA TABLE  /'.PRINT ', '.MOLGRA', '.SELECT',
C                      4          5          6
     &             '.DIPGRA', '.POLARI', '.MOLHES',
C                      7          8          9
     &             '.INPTES', '.VIBANA', '.RESTAR',
C                     10         11         12
     &             '.OPTROT', '.LINEAR', '.REPS  ',
C                     13         14         15
     &             '.NUMHES', '.POPANA', '.OECD  ',
C                     16         17         18
     &             '.EXCITA', '.SHIELD', '.VCD   ',
C                     19         20         21
     &             '.MAGNET', '.SPIN-S', '.NMR   ',
C                     22         23         24
     &             '.NACME ', '.NOLOND', '.NOHESS',
C                     25         26         27
     &             '.WRTINT', '.ECD   ', '.NODIFC',
C                     28         29         30
     &             '.MOLGFA', '.SPIN-R', '.OR    ',
C                     31         32         33
     &             '.PHASEO', '.NOMASV', '.NODARW',
C                     34         35         36
     &             '.ALPHA ', '.VROA  ', '.OR MVE',
C                     37         38         39
     &             '.SOPPA(', '.HELLMA', '.NOCMC ',
C                     40         41         42
     &             '.EXPFCK', '.GAUGEO', '.DIPORG',
C                     43         44         45
     &             '.RAMAN ', '.QUADRU', '.NQCC  ',
C                     46         47         48
     &             '.CAVORG', '.HYPER ', '.VERDET',
C                     49         50         51
     &             '.MCD   ', '.SOPPA ', '.EXPGRA',
C                     52         53         54
     &             '.CTOCD ', '.D2DQ2 ', '.QPGRAD',
C                     55         56         57
     &             '.DISDIP', '.SECMOM', '.THIRDM',
C                     58         59         60
     &             '.VIB_G '/ 
      DATA RSTLBL /'ORTREL', 'DIPCTL', 'RESPON', 'RELAX ',
     &             'DIPREL', 'TROINV', 'CCRHS ', 'CCTRES',
     &             'CCRELA', 'QPGCTL'/
      CALL QENTER('ABAINP')
C
C     Initialize /PAST/
C
      PASTWO = .FALSE.
      PASORT = .FALSE.
      PASDIP = .FALSE.
      PASQPG = .FALSE.
      PASONE = .FALSE.
      PASRES = .FALSE.
      PASTRP = .FALSE.
      PASLRS = .FALSE.
      PASREL = .FALSE.
      PASRTR = .FALSE.
      PASCRS = .FALSE.
      PASCZR = .FALSE.
      PASCTR = .FALSE.
      PASCRL = .FALSE.
      PASMAG = .FALSE.
      PASAAT = .FALSE.
      PASLNR = .FALSE.
      PASTRP = .FALSE.
C
C     Initialize /ABAINF/
C
      IPRDEF = IPRUSR
      NWNABA = 0
      POLAR  = .FALSE.
      TSTINP = .FALSE.
      RESTAR = .FALSE.
C
C     Initialization of MOLGRD, MOLHES, HELFEY, and DOWALK in /ABAINF/
C     is now done in GNRLIN, because control of the walk/geometry
C     optimization has been moved to GENERAL input processing.
C     For .WALK (OPTWLK) and some other cases the earlier settings have
C     to be changed to new defaults here:
C
      IF (.NOT.OPTNEW .AND. .NOT.NMWALK) THEN
C     ... if (not .OPTIMIZE and not (.NMDDRV or .PARNMD) ) then ...
         DOWALK = OPTWLK .AND. WRDSRC .NE. '**PROPE'
C        ... check if converged geometry ('**PROPE') in a .WALK
         MOLGRD = DOWALK
         MOLHES = DOWALK .OR. (TOTSYM .AND. WRDSRC .EQ. '**PROPE')
         IF (((V3CAL .AND. NMODIF) .OR. VIBAVE) .AND. 
     &       (WRDSRC .NE. '**START')) MOLHES = .FALSE.
         IF (VIBAVE .AND. REUSED) MOLHES = .FALSE.
      END IF
C
C     Initializing some properties depends on if this is a numerical 
C     differentiation. If so initialization is done elsewhere.
C      
      IF (.NOT.NMWALK .OR. (WRDSRC .EQ. '**PROPE')) THEN
         SPNSPN = .FALSE.
      END IF
      DODRCT = DIRCAL
      GRDINP = .FALSE.
      FCKEXP = .FALSE.
      DIPDER = .FALSE.
      QPGRAD = .FALSE.
      VIB    = .FALSE.
      SHIELD = .FALSE.
      MAGSUS = .FALSE.
      GDALL  = .FALSE.
      POPULA = .FALSE.
      DOLRES = .FALSE.
      VCD    = .FALSE.
      NACME  = .FALSE.
      NOLOND = .FALSE.
      NOHESS = .FALSE.
      ECD    = .FALSE.
      OECD   = .FALSE.
      SUPMAT = .FALSE.
      VROA   = .FALSE.
      NODIFC = .FALSE.
      MOLGFA = .FALSE.
      SPINRO = .FALSE.
      NUMHES = .FALSE.
      MASSVE = .TRUE.
      DARWIN = .TRUE.
      NOCMC  = .FALSE.
      RAMAN  = .FALSE.
      NQCC   = .FALSE.
      QUADRU = .FALSE.
      SECNDM = .FALSE.
      GAGT   = .FALSE.
      HYPER  = .FALSE.
      VERDET = .FALSE.
      MCD    = .FALSE.
      ABA_ALPHA = .FALSE.
      OPTROT = .FALSE.
      MVEOR  = .FALSE.
      LINCPL = .FALSE.
      SKIPAB = .FALSE.
      CTOCD  = .FALSE.
      EXPGRD = .FALSE.
      THIRDM = .FALSE.
      ESG    = .FALSE.
      VIB_G  = .FALSE.
      CALL DZERO(ORIGIN,3)
      CALL DZERO(CMXYZ,3)
      CALL DZERO(DIPORG,3)
      CALL DZERO(GAGORG,3)
      CALL DZERO(CAVORG,3)
      ABALNR = .FALSE.
      IF (TOTSYM .AND. WRDSRC .NE. '**PROPE') THEN
         DOSYM(1) = .TRUE.
         DO 40 I = 2, 8
            DOSYM(I) = .FALSE.
 40      CONTINUE
      ELSE
         DO I = 1, 8
            DOSYM(I) = .TRUE.
         END DO
      END IF
      FCKDDR = .TRUE.
      EXPFCK = .FALSE.
      ABASOP = .FALSE.
      CCPPA  = .FALSE.
      DOD2DQ2 = .FALSE.
C
C     Initialize /NUCLEI/
C
      NTRACO = 0
      DO 50 I = 1, 3
         ITRACO(I) = 0
   50 CONTINUE
      DO 60 I = 1, 3*MXCENT
         DOPERT(I,1) = .TRUE.
         DOPERT(I,2) = .TRUE.
   60 CONTINUE
C
C     Set FTRONV and FABRHS in /EXEINF/:
C
C     TROINV is run twice. FTRONV keeps track of it, and setting it TRUE
C     indicates first run.
C
      FTRONV = .TRUE.
C     
C     ABARHS has a variable of the same kind:
C
      FABRHS = .TRUE.
C
C     **** Find Abacus input *****
C
      CALL GPOPEN(LUCMD,'DALTON.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)   
      REWIND (LUCMD,IOSTAT=IOS)
C     ... IOSTAT to avoid program abort on some systems
C         if reading input from a terminal
      WORD1 = WRDSRC
  900 READ (LUCMD,'(A7)',END=910,ERR=920) WORD
      CALL UPCASE(WORD)
      IF (WORD .EQ. WORD1) THEN
         GO TO 930
      ELSE
         GO TO 900
      END IF
  910 CONTINUE
  920 CONTINUE
C     ... WRDSRC not found in input file;
C         if WRDSRC .eq. '**START' we also try for '**PROPER'
C         such that any '**PROPER' input will now also apply for initial
C         geometry if not overridden with '**START' /Feb 2003 hjaaj
         IF (WORD1 .EQ. '**START') THEN
            REWIND (LUCMD,IOSTAT=IOS)
            WORD1 = '**EACH '
            GO TO 900
         END IF
         CALL ABAINIALL
         GOTO 1
  930 CONTINUE
      WORD1 = WORD
C
      CALL TITLER('Output from '//WRDSRC//' input processing for ABACUS'
     &            ,'*',108)
C
C
C     ***** Process input for COMMON  /ABAINF/  *****
C
      INPSYM = 0
  100 READ (LUCMD, '(A7)') WORD
      CALL UPCASE(WORD)
 211  PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 100
      ELSE IF (PROMPT .EQ. '.') THEN
         DO 99 I = 1, NTABLE
            IF (TABLE(I) .EQ. WORD) THEN
               GO TO (101,102,103,104,105,106,107,108,109,110,
     &                111,112,113,114,115,116,117,118,119,120,
     &                121,122,123,124,125,126,127,128,129,130,
     &                131,132,133,134,135,136,137,138,139,140,
     &                141,142,143,144,145,146,147,148,149,150,
     &                151,152,153,154,155,156,157,158),I
            END IF
   99    CONTINUE
            IF (WORD .EQ. '.OPTION') THEN
              CALL PRTAB(NDIR,TABDIR, WORD1//' input keywords',LUPRI)
              CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
              GO TO 100
            END IF
            WRITE (LUPRI,'(/,3A,/)')
     *         ' Keyword ',WORD,' not recognized in ABAINP.'
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal keyword in ABAINP.')

  101    CONTINUE !        .PRINT
            READ (LUCMD,*) IPRDEF
            GO TO 100
  102    CONTINUE !        .MOLGRA
            MOLGRD = .TRUE.
            GRDINP = .TRUE.
            GO TO 100
  103    CONTINUE !        .SELECT
            READ (LUCMD,*) NPERT
            READ (LUCMD, *) (IPOINT(I),I=1,NPERT)
            CALL HEADER('In this run only the following '//
     *                  'perturbations are considered:',1)
            WRITE (LUPRI,'(5I5,5X,5I5)') (IPOINT(I), I = 1, NPERT)
            WRITE (LUPRI,'()')
            DOPERT(:,1) = .FALSE.
            DO I = 1, NPERT
               DOPERT(IPOINT(I),1) = .TRUE.
            END DO
            GO TO 100
  104    CONTINUE !        .DIPTRA
            DIPDER = .TRUE.
            GO TO 100
  105    CONTINUE !        .POLARI
            POLAR = .TRUE.
            GO TO 100
  106    CONTINUE !        .MOLHES
            MOLGRD = .TRUE.
            MOLHES = .TRUE.
            GO TO 100
  107    CONTINUE !        .INPTES
            TSTINP = .TRUE.
            GO TO 100
  108    CONTINUE !        .VIBANA
            VIB    = .TRUE.
            IF (IWKTYP .EQ. 6) THEN
               MOLGRD = .TRUE.
               MOLHES = .NOT.NUMHES
            ELSE
               MOLGRD = .NOT.((WRDSRC .EQ. '**PROPE') .AND. OPTWLK)
               MOLHES = .NOT.((WRDSRC .EQ. '**PROPE') .AND. OPTWLK 
     &                         .AND. .NOT. TOTSYM)
            END IF
            IF (MOLHES) DIPDER = .TRUE.
            GO TO 100
  109    CONTINUE !        .RESTAR
            RESTAR = .TRUE.
            READ (LUCMD,'(A12)') REWORD
            CALL UPCASE(REWORD)
            IF ((REWORD(1:6) .NE. '.FROM ') .AND.
     *          (REWORD(1:4) .NE. '.AT ')) THEN
               WRITE(LUPRI,'(4A)')
     *           '"',REWORD(1:6),'" is incorrect ',
     *           'Input, ".FROM " or ".AT " expected after ".RESTAR"'
               CALL QUIT('Input error in ABAINP.')
            ELSE
               IF (REWORD(1:6) .EQ. '.FROM ') THEN
                  RWORD = REWORD(7:12)
               ELSE IF (REWORD(1:4) .EQ. '.AT ') THEN
                  RWORD = REWORD(5:10)
               END IF
            END IF
         GO TO 100
  110    CONTINUE !        .OPTROT
            OPTROT = .TRUE.
            ABALNR = .TRUE.
         GO TO 100
  111    CONTINUE !        .LINEAR
            VIB    = .TRUE.
            LINCPL = .TRUE.
            MOLGRD = .TRUE.
            MOLHES = .FALSE.
         GO TO 100
  112    CONTINUE !        .REPS  
            READ (LUCMD,*) NREPS
            READ (LUCMD,*) (IDOSYM(I),I=1,NREPS)
            DOSYM(:) = .FALSE.
            DO ISYM = 1, NREPS
               DOSYM(IDOSYM(ISYM)) = .TRUE.
            END DO
            INPSYM = INPSYM + 1
         GO TO 100
  113    CONTINUE !        .NUMHES
            NUMHES = .TRUE.
         GO TO 100
  114    CONTINUE !        .POPANA
            POPULA = .TRUE.
            DIPDER = .TRUE.
            POLAR  = .TRUE.
         GO TO 100
  115    CONTINUE !        .OECD
            OECD   = .TRUE.
            NOCMC  = .TRUE.
            ECD    = .TRUE.
            DOEXCI = .TRUE.
         GO TO 100
  116    CONTINUE !        .EXCITA
            DOEXCI = .TRUE.
         GO TO 100
  117    CONTINUE !        .SHIELD
            SHIELD = .TRUE.
         GO TO 100
  118    CONTINUE !        .VCD
            VCD    = .TRUE.
            MOLGRD = .TRUE.
            MOLHES = .TRUE.
            DIPDER = .TRUE.
            VIB    = .TRUE.
            NOCMC  = .TRUE.
         GO TO 100
  119    CONTINUE !        .MAGNET
            MAGSUS = .TRUE.
            MOLGFA = .TRUE.
         GO TO 100
  120    CONTINUE !        .SPIN-S
            SPNSPN = .TRUE.
         GO TO 100
  121    CONTINUE !        .NMR
            SPNSPN = .TRUE.
            SHIELD = .TRUE.
         GO TO 100
  122    CONTINUE !        .NACME
            NACME  = .TRUE.
            MOLGRD = .TRUE.
            MOLHES = .TRUE.
         GO TO 100
  123    CONTINUE !        .NOLONDON
            NOLOND = .TRUE.
         GO TO 100
  124    CONTINUE !        .NOHESS
            NOHESS = .TRUE.
         GO TO 100
  125    CONTINUE
            FCKDDR = .FALSE.
         GO TO 100
  126    CONTINUE
            ECD    = .TRUE.
            DOEXCI = .TRUE.
         GO TO 100
  127    CONTINUE
            NODIFC = .TRUE.
         GOTO 100
  128    CONTINUE
            MOLGFA = .TRUE.
            MAGSUS = .TRUE.
         GOTO 100
  129    CONTINUE
            SPINRO = .TRUE.
         GOTO 100
  130    CONTINUE ! .OR = .OR MVE = .OPTROT + mod. vel. OR
            MVEOR  = .TRUE.
            OPTROT = .TRUE.
            ABALNR = .TRUE.
         GOTO 100
  131    CONTINUE
            READ (LUCMD,*) (ORIGIN(ICOOR), ICOOR = 1,3)
         GOTO 100
  132    CONTINUE
            MASSVE = .FALSE.
         GOTO 100
  133    CONTINUE
            DARWIN = .FALSE.
         GOTO 100
  134    CONTINUE
            ABA_ALPHA = .TRUE.
            ABALNR = .TRUE.
         GOTO 100
  135    CONTINUE
            VROA   = .TRUE.
            ABALNR = .TRUE.
         GOTO 100
  136    CONTINUE
            MVEOR  = .TRUE.
            OPTROT = .TRUE.
            ABALNR = .TRUE.
         GOTO 100
  137    CONTINUE
            ABASOP = .TRUE.
            CCPPA  = .TRUE.
         GOTO 100
  138    CONTINUE
            HELFEY = .TRUE.
         GOTO 100
  139    CONTINUE
            NOCMC = .TRUE.
         GO TO 100
  140    CONTINUE
            FCKEXP = .TRUE.
         GO TO 100
  141    CONTINUE!          .GAUGEO
            READ (LUCMD,*) (GAGORG(ICOOR), ICOOR = 1,3)
            GAGT = .TRUE.
         GOTO 100
 142     CONTINUE
            READ (LUCMD,*) (DIPORG(ICOOR), ICOOR = 1, 3)
         GOTO 100
 143     CONTINUE
            RAMAN  = .TRUE.
            ABALNR = .TRUE.
         GOTO 100
 144     CONTINUE
            QUADRU = .TRUE.
         GOTO 100
 145     CONTINUE
            NQCC = .TRUE.
         GOTO 100
 146     CONTINUE
            READ (LUCMD,*) (CAVORG(ICOOR), ICOOR = 1, 3)
         GOTO 100
cs
cs
 147     CONTINUE
            HYPER = .TRUE.
         GOTO 100
 148     CONTINUE
            VERDET = .TRUE.
         GOTO 100
 149     CONTINUE
            MCD = .TRUE.
         GOTO 100
 150     CONTINUE 
            ABASOP = .TRUE.
         GOTO 100
 151     CONTINUE 
            EXPGRD = .TRUE.
         GOTO 100
 152     CONTINUE 
            NOLOND = .TRUE.
            CTOCD  = .TRUE.
         GOTO 100
C .D2DQ2
 153     CONTINUE 
            DOD2DQ2 = .TRUE.
            MOLGRD  = .TRUE.
            MOLHES  = .TRUE.
         GOTO 100
 154     CONTINUE
            QPGRAD = .TRUE.
         GOTO 100
 155     CONTINUE
            POLAR  = .TRUE.
            DIPDER = .TRUE.
            QPGRAD = .TRUE.
         GOTO 100

 156     CONTINUE
            SECNDM = .TRUE.
         GOTO 100

 157     CONTINUE
            THIRDM = .TRUE.
         GOTO 100

 158     CONTINUE
            VIB_G  = .TRUE.
            NACME  = .TRUE.
            VIB    = .TRUE.
            MOLGRD = .TRUE.
            MOLHES = .TRUE.
            ABALNR = .TRUE.
         GOTO 100
      ELSE IF (PROMPT .EQ. '*') THEN
         GO TO 999
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal'
         CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal prompt in ABAINP.')
      END IF
  999 CONTINUE
      IF (MOLHES .AND. .NOT. DODRCT) EXPFCK = .TRUE.
      IF (NUMHES) THEN
         MOLHES = .FALSE.
         IF (.NOT. POLAR) DIPDER = .FALSE.
      END IF
      IF (FCKEXP) EXPFCK = .TRUE.
      IF (.NOT.(SHIELD.OR.SPNSPN)) THEN
         DO 65 I = 1, 3*MXCENT
            DOPERT(I,2) = .FALSE.
   65    CONTINUE
      END IF
C
C-----RESTART-------------------------------------------------
C
      IF (RESTAR) THEN
         CALL HEADER('This is a restart run!',0)
         CALL ABAREADRS
         IF (RWORD .EQ.'      ') THEN
            IF (DIPDER .AND. .NOT.PASDIP) PASONE = .FALSE.
            IF (QPGRAD .AND. .NOT.PASQPG) PASONE = .FALSE.
            IF (MOLHES .AND. .NOT.PASORT) PASONE = .FALSE.
            WRITE (LUPRI,'(/A/)')' Restart table used:'
            WRITE (LUPRI,'(2(/5X,A/,3X,8L7))')
     *        ' PASTWO PASORT PASDIP PASONE PASRES PASREL'//
     *        ' PASCRS PASCTR',
     *          PASTWO, PASORT, PASDIP, PASONE, PASRES,
     *          PASREL, PASCRS, PASCTR,
     *        ' PASMAG PASRTR PASAAT PASLRS PASTRP PASLNR'//
     *        ' PASEXC PASQPG',
     *          PASMAG, PASRTR, PASAAT, PASLRS, PASTRP,
     *          PASLNR, PASEXC, PASQPG
         ELSE
            WRITE (LUPRI,'(/2A/)') ' Non-default starting point: ',RWORD
            PASTWO = .FALSE.
            PASORT = .FALSE.
            PASDIP = .FALSE.
            PASQPG = .FALSE.
            PASONE = .FALSE.
            PASRES = .FALSE.
            PASTRP = .FALSE.
            PASREL = .FALSE.
            PASRTR = .FALSE.
            PASCRS = .FALSE.
            PASCTR = .FALSE.
            PASCZR = .FALSE.
            PASCRL = .FALSE.
            IF (RWORD .EQ. RSTLBL(1)) THEN
C              'ORTREL'
               PASTWO = .TRUE.
               PASCZR = .TRUE.
            ELSE IF (RWORD .EQ. RSTLBL(2)) THEN
C              'DIPCTL'
               PASTWO = .TRUE.
               PASORT = .TRUE.
               PASCZR = .TRUE.
            ELSE IF (RWORD .EQ. RSTLBL(3)) THEN
C              'RESPON'
               PASTWO = .TRUE.
               PASORT = .TRUE.
               PASDIP = .TRUE.
               PASONE = .TRUE.
               PASCZR = .TRUE.
            ELSE IF (RWORD .EQ. RSTLBL(4)) THEN
C              'RELAX '
               PASTWO = .TRUE.
               PASORT = .TRUE.
               PASDIP = .TRUE.
               PASONE = .TRUE.
               PASRES = .TRUE.
               PASCZR = .TRUE.
            ELSE IF (RWORD .EQ. RSTLBL(5)) THEN
C              'DIPREL'
               PASTWO = .TRUE.
               PASORT = .TRUE.
               PASDIP = .TRUE.
               PASONE = .TRUE.
               PASRES = .TRUE.
               PASREL = .TRUE.
               PASCZR = .TRUE.
            ELSE IF (RWORD .EQ. RSTLBL(6)) THEN
C              'TROINV'
               PASTWO = .TRUE.
               PASORT = .TRUE.
               PASCRS = .TRUE.
               PASCTR = .TRUE.
               PASDIP = .TRUE.
               PASONE = .TRUE.
               PASRES = .TRUE.
               PASREL = .TRUE.
               PASCZR = .TRUE.
               PASCRL = .TRUE.
            ELSE IF (RWORD .EQ. RSTLBL(7)) THEN
C              'CCRHS '
               PASTWO = .TRUE.
               PASORT = .TRUE.
               PASDIP = .TRUE.
               PASONE = .TRUE.
               PASRES = .TRUE.
               PASREL = .TRUE.
               PASCZR = .TRUE.
            ELSE IF (RWORD .EQ. RSTLBL(8)) THEN
C              'CCTRES'
               PASTWO = .TRUE.
               PASORT = .TRUE.
               PASDIP = .TRUE.
               PASONE = .TRUE.
               PASRES = .TRUE.
               PASREL = .TRUE.
               PASCZR = .TRUE.
               PASCRS = .TRUE.
            ELSE IF (RWORD .EQ. RSTLBL(9)) THEN
C              'CCTREL'
               PASTWO = .TRUE.
               PASORT = .TRUE.
               PASDIP = .TRUE.
               PASONE = .TRUE.
               PASRES = .TRUE.
               PASREL = .TRUE.
               PASCZR = .TRUE.
               PASCRS = .TRUE.
               PASCTR = .TRUE.
            ELSE IF (RWORD .EQ. RSTLBL(10)) THEN
C              'QPGCTL'
               PASTWO = .TRUE.
               PASORT = .TRUE.
               PASCZR = .TRUE.
            ELSE
               WRITE (LUPRI,'(/3A//A//,(T6,A))')
     *            ' Starting point for restart "',RWORD,
     *            '" unknown in ABAINP.',' Known starting points:',
     *            RSTLBL
               WRITE (LUPRI,'(A)') ' or "      " (blank)'
               CALL QUIT('ERROR in restart input in ABAINP.')
            END IF
         END IF
      END IF
cs -------------------------------endofrestart---------------------
C
cs    write session
cs
      IF (.NOT. ABASOP .AND. DOMP2 .AND. (NASHT .EQ. 0)) THEN
         WRITE(LUPRI,'(A)') ' Properties not ' //
     &        'available for MP2 wave functions!'
         WRITE(LUPRI,'(A)') ' You may consider using SOPPA instead.'
         SKIPAB = .TRUE.
         GOTO 789
      END IF
      IF (NOHESS .AND. MOLHES) THEN
         WRITE (LUPRI,'(5X,A)')
     &      ' No calculation of Hessian due to .NOHESS.'
         MOLHES = .FALSE.
      END IF
      IF ((VROA .OR. RAMAN) .AND. DOWALK) THEN
         MOLHES = .FALSE.
         MOLGRD = GRDINP .OR. NUMHES
         IF (DIPDER) POLAR = .TRUE.
      END IF
      IF (ABASOP) THEN
         IF (CCPPA) THEN
            WRITE(LUPRI,'(A,A,/,18X,A,/,A,A)')
     &          ' SOPPA(CCSD) :    ',
     &          'Second Order Polarization Propagator Approximation',
     &          'with Coupled Cluster Singles and Doubles Amplitudes',
     &          '    (Ref.: S.P.A. Sauer, ',
     &          'J. Phys. B: At. Mol. Opt. Phys. 30, 3773-3780 (1997))'
CSPAS : 6/11-2009 AO-SOPPA included
            IF (DODRCT) THEN
               AOSOC = .TRUE.
               WRITE(LUPRI,'(/A,/,A,A,/,5x,A)')
     &              ' AO integral driven SOPPA(CCSD) :',
     &              '    (Ref.: H. H. Falden, K. R. Falster-Hansen, ',
     &              'K. L. Bak, S. Rettrup and S. P. A. Sauer,',
     &              'J. Phys. Chem. A 113, 11995-12012 (2009))'
            ENDIF
CKeinSPASmehr
         ELSE 
CSPAS : 6/11-2009 AO-SOPPA included
            IF (DODRCT) THEN
               AOSOP = .TRUE.
               WRITE(LUPRI,'(//5x,A,/,5x,A,A,/,5x,A)')
     &              ' AO integral driven SOPPA :',
     &              '    (Ref.: K. L. Bak, H. Koch, J. Oddershede,',
     &              ' O. Christiansen and S. P. A. Sauer',
     &              '     J. Chem. Phys. 112, 4173-4185 (2000))'
            ELSE
               WRITE(LUPRI,'(//5x,A,A,/)')
     *             ' SOPPA :',
     *             ' Second Order Polarization Propagator Approximation'
            END IF
CKeinSPASmehr
         END IF
         IF (SHIELD .OR. MAGSUS .OR. SPINRO .OR. MOLGFA) THEN
            NOLOND = .TRUE.
         END IF
C
C     No mass-velocity or Darwin for SOPPA calc.
C
         MASSVE = .FALSE.
         DARWIN = .FALSE.
      END IF

      IF (QMMM) WRITE(LUPRI,*) 'QMMM calculation'

      CALL HEADER('The following molecular properties will be'//
     *            ' calculated in this run:',0)
      WRITE (LUPRI,'(A,I5/)') ' Default print level:    ',IPRDEF
      IF (TSTINP) WRITE (LUPRI,'(/A/)')  ' *** Input test run only ***'
      IF (MOLGRD) WRITE (LUPRI,'(5X,A)') ' The molecular gradient'
      IF (MOLHES) WRITE (LUPRI,'(5X,A)') ' The molecular Hessian'
      IF (DOWALK) WRITE (LUPRI,'(5X,A)') ' Geometry optimization'
      IF (VIB)    WRITE (LUPRI,'(5X,A)') ' Vibrational analysis'
      IF (SECNDM) WRITE (LUPRI,'(5X,A)') ' Second order moments'
      IF (QUADRU) WRITE (LUPRI,'(5X,A)') ' Quadrupole moments'
      IF (THIRDM) WRITE (LUPRI,'(5X,A)') ' Third order moments'
      IF (DIPDER .AND. .NOT.POPULA)
     &            WRITE (LUPRI,'(5X,A)') ' Dipole moment gradient'
      IF (QPGRAD) WRITE (LUPRI,'(5X,A)') ' Quadrupole moment'//
     &     ' gradient'
      IF (HELFEY) THEN
         WRITE (LUPRI,'(5X,A)') ' The Hellman-Feynman '//
     &     ' theorem will be used'
         WRITE (LUPRI,'(5X,A)') 'This option is currently not working'//
     &        ' correctly, program will stop'
         CALL QUIT('Hellmann-Feynman approximation not working')
      END IF
      IF (POLAR .AND. .NOT.POPULA)
     &            WRITE (LUPRI,'(5X,A)') ' Static polarizabilities'
      IF (HYPER)  WRITE (LUPRI,'(5X,A)') 
     &          ' First hyperpolarizabilities (beta)'
      IF (NQCC)   WRITE (LUPRI,'(5X,A)') ' Nuclear quadrupole moments'
      IF (POPULA) WRITE (LUPRI,'(5X,A)') ' Cioslowski population'
     &               //' analysis and static polarizabilities only.'
      IF (DOLRES) WRITE (LUPRI,'(5X,A)') ' Linear response calculation'
      IF (DOEXCI) WRITE (LUPRI,'(5X,A)')
     &               ' Electronic excitation energies '
      IF (SHIELD) WRITE (LUPRI,'(5X,A)') ' Nuclear magnetic shieldings'
      IF (SPNSPN) WRITE (LUPRI,'(5X,A)')
     &               ' Nuclear spin-spin coupling constants'
      IF (MAGSUS) WRITE (LUPRI,'(5X,A)') ' Magnetic susceptibilities'
cLig >> added what to print when the CTOCD option is set true
      IF (CTOCD) THEN
         WRITE (LUPRI,'(5X,A)') ' CTOCD-DZ: '
         WRITE (LUPRI,'(5X,A)')
     &'  Continuous Transformation of the Origin of the Current Density'
         IF(MAGSUS) WRITE(LUPRI,'(5X,A)')
     &     ' CTOCD-DZ Magnetic susceptibilities'
         IF(SHIELD) WRITE(LUPRI,'(5X,A)')
     &     ' CTOCD-DZ Nuclear magnetic shieldings' 
         WRITE (LUPRI,'(5x,A,A,/,5x,A,A)')
     &          '    (Ref.: A. Ligabue, S. P. A. Sauer',
     &          ' and P. Lazzeretti,',
     &          '     J. Chem. Phys. 118, 6830 (2003)',
     &          ' & J. Chem. Phys. 126, 154111 (2007))' 
      END IF
cLig <<
      IF (VCD)    WRITE (LUPRI,'(5X,A)')
     &               ' Vibrational circular dichroism (VCD)'
      IF (NACME)  WRITE (LUPRI,'(5X,A)')
     &   ' First-order non-adiabatic coupling matrix elements (NACMEs)'
      IF (VIB_G)  WRITE (LUPRI,'(5X,A,/,5x,A,A,/,5x,A)')
     &            ' Vibrational g-factors :',
     &              '    (Ref.: K. L. Bak, S. P. A. Sauer,',
     &              ' J. Oddershede and J. F. Ogilvie,',
     &              '     Phys. Chem. Chem. Phys. 7, 1747 (2005))' 
      IF (DOD2DQ2)  WRITE (LUPRI,'(5X,A)')
     &   ' < d Psi/dQ | d Psi/dQ > (Q = nuclear coordinates)'
      IF (NOLOND) WRITE (LUPRI,'(5X,A)')
     &   ' INFO - London orbitals are NOT used'
      IF (DODRCT) WRITE (LUPRI,'(5X,A)')
     &   ' Direct calculation - no two-electron integrals '//
     &   'written or read.'
      IF (EXPFCK .AND. DODRCT) WRITE (LUPRI,'(5X,A)')
     &   ' Two-electron expectation values and derivative Fock ',
     &   ' matrices are calculated simultaneously (for direct SCF).'
      IF (SUPMAT) WRITE (LUPRI,'(5X,A)')
     &   ' Two-electron supermatrix integrals are used.'
      IF (PARCAL) WRITE (LUPRI,'(5X,A)')
     &   ' Two-electron integrals will be calculated in parallel'
      IF (LINCPL) THEN
         WRITE (LUPRI,'(A)') ' Analyze the gradient in normal '//
     &        'coordinates read from Hessian on file'
      END IF
      IF (ABA_ALPHA) WRITE (LUPRI,'(5X,A)')
     &            ' Dynamic electric dipole polarizability'
      IF (OECD)   WRITE (LUPRI,'(5X,A)')
     &            ' Oriented electronic circular dichroism'
      IF (ECD)    WRITE (LUPRI,'(5X,A)')
     &            ' Electronic circular dichroism'
      IF (MCD)    WRITE (LUPRI,'(5X,A)')
     &            ' Magnetic circular dichroism'
      IF (VERDET) WRITE (LUPRI,'(5X,A)')
     &            ' Verdet constants'
      IF (VROA)   WRITE (LUPRI,'(5X,A)')
     &            ' Vibrational Raman intensities and optical activity'
      IF (RAMAN)  WRITE (LUPRI,'(5X,A)')
     &            ' Vibrational Raman intensities'
      IF (NUMHES .AND. (VROA .OR. RAMAN))  WRITE (LUPRI,'(5X,A)')
     &            ' A numerically calculated Hessian will be used '//
     &            'for the frequency analysis'
      IF (OPTROT) THEN
         IF (MVEOR) THEN
            WRITE (LUPRI,'(5X,A)')
     &      ' Optical rotation strengths (London, length, velocity)'
         ELSE
            WRITE (LUPRI,'(5X,A)')
     &      ' Optical rotation strengths (London, length)'
         END IF
      END IF
      IF (SPINRO) WRITE(LUPRI,'(5X,A)') ' Spin-rotation constants'
      IF (MOLGFA) WRITE(LUPRI,'(5X,A)') ' Molecular rotational g-factor'
      IF (EXPGRD) WRITE(LUPRI,'(5X,A)') ' Orbital-exponent gradient' 
      IF (.NOT. MASSVE) WRITE(LUPRI,'(5X,A)')
     &     ' No mass-velocity correction calculated'
      IF (.NOT. DARWIN) WRITE(LUPRI,'(5X,A)')
     &     ' No Darwin correction calculated'
      IF (NODIFC) THEN
         WRITE(LUPRI,'(6X,A)')
     &   'Differentiatied creation operators not used (.NODIFC),',
     &   'i.e. symmetric orbital connection is used',
     &   'for perturbation dependent basis sets.'
      ELSE
         WRITE(LUPRI,'(6X,A)')
     &   'Natural orbital connection is used',
     &   'for perturbation dependent basis sets.'
      END IF
      IF (GAGT) THEN
         WRITE (LUPRI,'(A,3F15.10)')
     &        ' User defined gauge origin       :',
     &        (GAGORG(I), I = 1, 3)
         IF (MOLGFA .OR. SPINRO .OR. ((MAGSUS .OR. SHIELD)
     &              .AND. NOLOND)) THEN
            WRITE(LUPRI,'(2A)')
     &           ' WARNING! Gauge origin changed from default. ',
     &           'This may give wrong results for some properties'
         END IF
         IF ((MAGSUS .OR. SHIELD) .AND. .NOT. NODIFC)
     &      WRITE(LUPRI,'(1X,2A)') ' WARNING! User defined gauge '//
     &           'may give wrong results when using symmetry'
      END IF
      IF (DNRM2(3,CAVORG,1) .GT. D0) THEN
         WRITE (LUPRI,'(A,3F15.10)')
     &        ' User defined cavity origin      :',
     &        (CAVORG(I), I = 1, 3)
      END IF
      IF (DNRM2(3,DIPORG,1) .GT. D0) THEN
         WRITE (LUPRI,'(A,3F15.10/A)')
     &        ' User defined dipole origin      :',
     &        (DIPORG(I), I = 1, 3),
     &        ' (Also used for .QUADRU, .SECMOM, .THIRDM multipoles)'
      END IF
      IF (DNRM2(3,ORIGIN,1) .GT. D0) THEN
        WRITE(LUPRI,'(A,3F15.10)')
     &         ' Origin for London phase factor :',
     &        (ORIGIN(I), I = 1, 3)
      END IF
      IF (QPGRAD .AND. .NOT.NOCMC) THEN
         WRITE(LUPRI,'(A)') ' Incompatible input. Quadrupole '//
     &        'gradient requires NOCMC'
         CALL QUIT('Inconsistent input in ABACUS')
      END IF
      IF (ABALNR) WRITE (LUPRI,'(5X,A)') ' Linear response properties'
      IF (INPSYM .GT. 0) THEN
         IF (NREPS .EQ. 0) THEN
            WRITE (LUPRI,'(/2(/5X,A))')
     &           ' Input error in ABAINP.',
     &           ' No symmetries specified under .REPS'
            CALL QUIT('Input error in ABAINP.')
         ELSE IF (NREPS .EQ. 1) THEN
            WRITE (LUPRI,'(//A,I2,A)')
     &           ' Only distortions belonging to symmetry',IDOSYM(1),
     &           ' are considered.'
         ELSE
            WRITE (LUPRI,'(//A,I1,A,I1,7(A,I1))')
     &           ' Only distortions of the following ',
     &           NREPS,' symmetries are considered: ',IDOSYM(1),
     &           (',',IDOSYM(I),I=2,NREPS)
         END IF
      ELSE IF (TOTSYM) THEN
         WRITE (LUPRI,'(//A)')
     &        ' Only totally symmetric distortions are considered.'
      END IF
      IF ((POPULA .OR. POLAR .OR. DIPDER .OR. QPGRAD .OR. VCD) .AND. 
     &    (TOTSYM .AND. WRDSRC .NE. '**PROPE')) THEN
         NWARN = NWARN + 1
         WRITE(LUPRI,'(/A)') 'WARNING, incompatible input options:'
         IF (POPULA) WRITE (LUPRI,'(2(5X,A))') 
     &      'Cioslowski population analysis incompatible with .TOTSYM.',
     &      '  WARNING: Population analysis is not carried out.'
         IF (POLAR) WRITE (LUPRI,'(2(5X,A))') 
     &      'Static polarizability incompatible with .TOTSYM.',
     &      '  WARNING: No calculation of static polarizabilities.'
         IF (DIPDER) WRITE (LUPRI,'(2(5X,A))') 
     &      'Dipole gradient incompatible with .TOTSYM.',
     &      '  WARNING: Dipole gradient is not calculated.'
         IF (QPGRAD) WRITE (LUPRI,'(2(5X,A))') 
     &      'Quadrupole gradient incompatible with .TOTSYM.',
     &      '  WARNING: Quadupole gradient is not calculated.'
         IF (VCD) WRITE (LUPRI,'(2(5X,A))') 
     &      'VCD incompatible with .TOTSYM.',
     &      '  WARNING: VCD tensors are not calculated.'
         VCD    = .FALSE.
         POPULA = .FALSE.
         DIPDER = .FALSE.
         QPGRAD = .FALSE.
         POLAR  = .FALSE.
      END IF
      IF (DIPDER .AND. .NOT. (MOLHES .OR. POLAR)) THEN
         WRITE (LUPRI,'(//2A/A)')
     *      ' FATAL ERROR .DIPOLE: Dipole moment gradient requires',
     *      ' either molecular Hessian or polarizabilities. ',
     *      ' Set ".POLARI" or delete ".NOHESS" in input!'
         CALL QUIT
     *      ('Inconsistency in input for ".DIPOLE", cannot proceed.')
      END IF
      IF (SPNSPN .AND. NASHT .LE. 1 .AND. .NOT.(ABASOP.OR.DFTRUN)) THEN
         NWARN = NWARN + 1
         WRITE (LUPRI,'(/A/A)')
     &      ' >>>>> WARNING - RHF calculations of spin-spin coupling'//
     &      ' constants are likely to give', 
     &      ' >>>>> WARNING - qualitatively incorrect results !!!'
      END IF
cLig >> added the control to right CTOCD input
      IF ((CTOCD).AND.(.NOT.(MAGSUS .OR. SHIELD))) THEN
         WRITE(LUPRI,'(/A/A)')
     &      ' Inconsistent input - CTOCD-DZ method',
     &      ' the .CTOCD option only works if the .MAGNET and/or'//
     &      ' .SHIELD options are also used'
         CALL QUIT('Inconsistent input in ABAINP')
      ENDIF
      IF ((VIB_G).AND.
     &    ((OPTROT .OR. VROA .OR. RAMAN .OR. ABA_ALPHA))) THEN
         WRITE(LUPRI,'(/A,A/A,A)')
     &      ' Unsupported input - calculation of vibrational g-factor',
     &      ' is not tested to work',' if .RAMAN, .OPTROT, .VROA,'//
     &      ' .OR, .OR MVE or .ALPHA options are also used'
         CALL QUIT('Not supported input in ABAINP')
      ENDIF
cLig <<
      IF (ABASOP) THEN
         IERR = 0
         IF (MOLGRD) THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. MOLGRD '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (MOLHES) THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. MOLHES '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (DIPDER .AND. .NOT.POPULA) THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. DIPDER '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (QPGRAD .AND. .NOT.POPULA) THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. QPGRAD '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (HYPER)  THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. HYPER '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (QUADRU) THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. QUADRU '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (SECNDM) THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. SECMOM '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (NQCC)   THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. NQCC '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (THIRDM)   THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. THIRDM '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (DOWALK) THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. DOWALK '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (VIB)    THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. VIB '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (POPULA) THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. POPULA '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (VCD)    THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. VCD '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (NACME)  THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. NACME '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (OECD)    THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. OECD '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (ECD)    THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. ECD '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (MCD)    THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. MCD '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (VERDET) THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. VERDET '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (VROA)   THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. VROA '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
         IF (RAMAN)  THEN
            WRITE(LUPRI,'(A)') ' Incompatible input. RAMAN '//
     &                         'cannot be calculated with SOPPA'
            IERR = IERR + 1
         END IF
CHJ feb 2007 let us try if it now works ... /hjaaj
CHJ      IF (DODRCT) THEN
CHJ         WRITE (LUPRI,'(2X,A,/,2X,A)')
CHJ  &      ' Inconsistent input - SOPPA calculations cannot yet',
CHJ  &      ' be carried out "directly" (using .DIRECT and .SOPPA).'
CHJ         IERR = IERR + 1
CHJ      END IF
         IF (IERR .GT. 0) THEN
            CALL QUIT('Inconsistent input in ABACUS for SOPPA')
         END IF
      END IF
CAMT--------------
CAMT DFT Contributions for Magnetizabilities do not work with symmetry
CAMT g-tensors and spin-rotation constants also do not work, put a 
CAMT quit here to catch this...
      IF ((NSYM.GT.1).AND.
     &    (MAGSUS.OR.MOLGFA.OR.SPINRO).AND.(DFTRUN)) THEN
            IF (MAGSUS) THEN
              WRITE(LUPRI,'(5X,A50)')
     &           'WARNING: Magnetic Suceptibilities requested with '
              WRITE(LUPRI,'(5X,A50)')
     &           '         DFT and symmetry but not yet implemented'
            ENDIF
            IF (MOLGFA) THEN
              WRITE(LUPRI,'(5X,A50)')
     &           'WARNING: Rotational g-tensors requested with     '
              WRITE(LUPRI,'(5X,A50)')
     &           '         DFT and symmetry but not yet implemented'
            ENDIF
            IF (SPINRO) THEN
              WRITE(LUPRI,'(5X,A50)')
     &           'WARNING: Spin-rotation constants requested with  '
              WRITE(LUPRI,'(5X,A50)')
     &           '         DFT and symmetry but not yet implemented'
            ENDIF
              WRITE(LUPRI,'(5X,A50)')
     &           '         Please re-run with NoSymmetry...        '
            CALL QUIT('Specified Mag. Props. NYI with Symm for DFT')
      ENDIF
CAMT--------------

C
C
C
C     **** PROCESS INPUT FOR VARIOUS PROGRAM SECTIONS  *****
C
cs    introdotta *QRPROP (posizione 25 di TABDIR)
cs    legge input e vede cosa e' attivo
cs
      CALL ABAINIALL
 200  PROMPT = WORD(1:1)
      IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
         GO TO 200
      ELSE IF (PROMPT .EQ. '*') THEN
         DO 210 I = 1, NDIR
            IF (WORD .EQ. TABDIR(I)) THEN
               GO TO
     *         (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     *              22,23,24,25,26,27,28), I
            END IF
  210    CONTINUE
         IF (WORD(1:2) .EQ. '**') GO TO 1
         WRITE (LUPRI,'(/,3A,/)') ' Directory ',WORD,' nonexistent.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Illegal directory in ABAINP.')
      ELSE
         WRITE (LUPRI,'(/,3A,/)') ' Prompter "',PROMPT,'" illegal or',
     *                        ' out of order.'
         CALL PRTAB(NDIR,TABDIR,WORD1//' input keywords',LUPRI)
         CALL QUIT('Program stopped in ABAINP, error in prompt.')
      END IF
    2   CONTINUE
      GO TO 200
    3   CALL NUCINP(WORD)
      GO TO 200
    4   CALL ONEINP(WORD)
      GO TO 200
    5   CALL TWOINP(WORD)
      GO TO 200
    6   CALL RHSINP(WORD)
      GO TO 200
    7   CALL EXCINP(WORD)
      GO TO 200
    8   CALL ORTINP(WORD)
      GO TO 200
    9   CALL RESINP(WORD)
      GO TO 200
   10   CALL TROINP(WORD)
      GO TO 200
   11   CALL FMMINP(WORD)
      GO TO 200
   12   CALL DIPINP(WORD)
      GO TO 200
   13   CALL RLXINP(WORD)
      GO TO 200
   14   CALL ANAINP(WORD)
      GO TO 200
   15   CALL VIBINP(WORD)
      GO TO 200
   16   CALL FLTINP(WORD)
      GO TO 200
   17   CONTINUE
      GO TO 200
   18   CONTINUE
         CALL QPGINP(WORD)
      GO TO 200
   19   CONTINUE
      GO TO 200
   20   CONTINUE
      GO TO 200
   21   CALL SOPINP(WORD)
      GO TO 200
   22   CALL LRSINP(WORD)
      GO TO 200
   23   CALL MAGINP(WORD)
      GO TO 200
   24   CALL SPNINP(WORD)
      GO TO 200
   25   CALL AQRINP(WORD)
      GO TO 200
   26   CALL AATINP(WORD)
      GO TO 200
   27   CALL TRPINP(WORD)
      GO TO 200
   28   CALL LNRINP(WORD)
      GO TO 200
C
    1 CONTINUE
C
C     Properties not available for MP2, ABACUS run will be skipped
C
      IF (.NOT. ABASOP .AND. DOMP2 .AND. (NASHT .EQ. 0)) THEN
         SKIPAB = .TRUE.
         GOTO 789
      END IF
C
C
C     Read geometry, orbital spec., etc. (MOLECULE format)
C
C     Logical argument indicates whether information on LUONEL should be
C     written or not.
C
C     Initialize and define wave function variables in SIRIUS common
C     blocks. Print level in sirius routines (IPRIN4, IPRIN6) and for
C     integral transformation in TRACTL (IPRTRA also specifies print
C     level in DERTRA) can be specified in RHSINP and are here trans-
C     ferred from /CBIRHS/ to Sirius common blocks.
C
Chj   CALL READIN(WORK,LWORK,.FALSE.)
Chj   this call had survived from when ABACUS was a separate program.
Chj   I have now made sure that READIN always is called in the beginning,
Chj   also if RNHERM is false. /April 2009 hjaaj

      CALL SETSIR(WORK,LWORK)
C
C     If SOPPA calculation, then
C     1) set SOPPA true in /INFRSP/
C     2) define other variables for SOPPA in SETSOPPA and RSPACT
C     3) reset NCONMA and NVARMA for SOPPA calculation
C     This must follow after the 'CALL SETSIR' above, which (with
C     correct input) sets variables to HF values when SOPPA.
C     Aug. 97/spas+hjaaj
C
      IF (ABASOP) THEN
         SOPPA = ABASOP
         KIADR1 = 1
         KLAST  = KIADR1 + (MAXOCC * MAXVIR + 1)/IRAT
         IF (KLAST .GT. LWORK) CALL STOPIT('ABADRV','MP2SET',KLAST,
     &                                     LWORK)
         CALL MP2SET(WORK(KIADR1))
         CALL SETSOPPA
         NCONMA = 0
         DO ISYM = 1, NSYM
            NCONMA = MAX(NCONMA,MAX(N2P2HS(ISYM),N2P2HT(ISYM)))
         ENDDO
         NVARMA = NCONMA + NWOPMA
      ENDIF
      CALL RSPACT
C
C     The correct ABACUS execution variables are saved for later iterations.
C     ekd+hjaaj: PRESTR never called, so PSAVE commented out. (960306)
C
C     CALL PSAVE
      IF (MOLHES) CALL ZERHES
C
      IF ( (NASHT .LT. 1) .AND.
     &     (MOLHES .OR. DIPDER .OR. QPGRAD .OR. VCD .OR. MAGSUS .OR.
     &     ECD .OR. VROA .OR. OPTROT .OR. SHIELD .OR. SPINRO)) THEN
C
C	We add this line when direct, non-symmetric SKLFCK is fixed,KR
C
C     &           .OR. MCD .OR. VERDET )
         IF (.NOT. FCKDDR) WRITE (LUPRI,'(/A)')
     &     ' Derivative integrals will be written to disk as requested.'
      ELSE
C        ... FCKDDR cannot be used in this calculation
         FCKDDR = .FALSE.
      END IF
C
      IF ((VROA .OR. RAMAN) .AND. MAXREP .GT. 0) THEN
         WRITE (LUPRI,'(2X,A,I5/A/A)')
     &   ' Number of symmetries:',MAXREP + 1,
     &   ' Vibrational Raman calculations can only be run without '//
     &   'symmetry.',' Calculation aborted.'
         CALL QUIT('ERROR: Inconsistent input (Raman with symmetry)')
      END IF
C
      IF (DKTRAN .AND. MOLGRD) THEN
         WRITE(LUPRI,'(/A)') 'ERROR: molecular gradient not implemented'
     &   //' for Douglas-Kroll-Hess Hamiltonian'
         CALL QUIT('ERROR: mol. grad. not implemented for DKH(2)')
      END IF
      IF (DKTRAN .AND. .NOT.NOLOND) THEN
         WRITE(LUPRI,'(/A)') 'ERROR: London orbitals not implemented'
     &   //' for Douglas-Kroll-Hess Hamiltonian'
         CALL QUIT('ERROR: London orbitals not implemented for DKH(2)')
      END IF
C
C     LUSUPM = -1 means AOSUPMAT not used.
C     AOSUPMAT is only used for non-direct calculations when
C     requested by .SUPMAT under *ABACUS input.
C     ekd+hjaaj: LUSUPM may be reset by SIRIUS, and the assignment of
C     LUSUPM must therefore be repeated for each ABACUS call in geometry
C     iterations, also if no '**PROPER' input has been given.
C     This should be OK with this location.
C
      LUSUPM = - 1
      IF (SUPMAT .AND. .NOT.DODRCT) LUSUPM = -1
C
C     We need to change the origin to the center of mass when calculating
C     spin-rotation constants and molecular g-factors.
C     This is now default for all molecular properties, unless turned off.
C
      KGEOM = 1
      KMASS = KGEOM + 3*(NATOMS + NFLOAT)
      KNAT  = KMASS + NATOMS + NFLOAT
      KNUMIS= KNAT  + (NATOMS + NFLOAT + 1)/IRAT
      KLAST = KNUMIS+ (NATOMS + NFLOAT + 1)/IRAT
      IF (KLAST .GT. LWORK) CALL STOPIT('ABADRV','CMMASS',KLAST,
     &                                   LWORK)
      CALL CMMASS(WORK(KGEOM),WORK(KMASS),WORK(KNAT),WORK(KNUMIS),
     &            IPRDEF-2)
      IF (SHIELD .OR. MAGSUS .OR. VCD .OR. ECD .OR. OECD .OR. VROA .OR.
     &    SPINRO .OR. MOLGFA .OR. MCD .OR. VERDET .OR. OPTROT) THEN
         IPRGAG = 1
      ELSE
         IPRGAG = 0
      END IF
      IF (.NOT.NOCMC) THEN
         IF (DNRM2(3,DIPORG,1) .EQ. D0) THEN
C           ... i.e. user hasn't changed DIPORG with .DIPORG
            CALL DCOPY(3,CMXYZ,1,DIPORG,1)
            WRITE(LUPRI,'(/A,3F12.6)')
     &         ' Center of mass dipole origin  :',(DIPORG(I), I = 1, 3)
         END IF
         IF (GAGT) THEN
            IF (IPRGAG .GT. 0) WRITE(LUPRI,'(/A,3F12.6)') 
     &         ' User defined gauge origin used:',(GAGORG(I), I = 1, 3)
         ELSE
            CALL DCOPY(3,CMXYZ,1,GAGORG,1)
            IF (IPRGAG .GT. 0)
     &         WRITE(LUPRI,'(/A,3F12.6)')
     &         ' Center of mass gauge origin   :',(GAGORG(I), I = 1, 3)
         END IF
         IF (SOLVNT) THEN
            CALL DCOPY(3,CMXYZ,1,CAVORG,1)
            WRITE(LUPRI,'(/,A,3F12.6)')
     &         ' Center of mass cavity origin  :',(CAVORG(I), I = 1, 3)
         END IF
      ELSE
         CALL DZERO(GAGORG,3)
         IF (IPRGAG .GT. 0)
     &   WRITE(LUPRI,'(/,A,3F12.6)') ' Gauge  origin              :',
     &        (GAGORG(I), I = 1, 3)
         IF (SOLVNT) THEN
            CALL DZERO(CAVORG,3)
            WRITE(LUPRI,'(/,A,3F12.6)') ' Cavity origin              :',
     &           (CAVORG(I), I = 1, 3)
         END IF
      END IF
      IF (SOLVNT) THEN
         CORD(1,NUCIND) = CAVORG(1)
         CORD(2,NUCIND) = CAVORG(2)
         CORD(3,NUCIND) = CAVORG(3)
      END IF
      IF (SPINRO) CALL NUCSPR(WORK(KGEOM),CMXYZ,IPRDEF)
      IF (MOLGFA) CALL NUCMGF(WORK(KGEOM),GAGORG,LUPRI,IPRDEF)
      IF (SECNDM) CALL NUCSCM(WORK(KGEOM),DIPORG,LUPRI,IPRDEF)
      IF (QUADRU) CALL NUCQDR(WORK(KGEOM),DIPORG,LUPRI,IPRDEF)
      IF (NQCC)   CALL NUCNQC(WORK(KGEOM),LUPRI,IPRDEF)
      IF (THIRDM) CALL NUC3RD(WORK(KGEOM),DIPORG,LUPRI,IPRDEF)
C     
 789  CONTINUE
      CALL GPCLOSE(LUCMD,'KEEP')
      IF (TSTINP) CALL QUIT('*** End of input test for ABACUS ***')
      CALL FLSHFO(LUPRI)
      CALL QEXIT('ABAINP')
      RETURN
C
      END
C  /* Deck result */
      SUBROUTINE ABARESULT(POLDD,POLDL,POLDA,POLVL,POLVV,
     &                  CSTRA,SCTRA,SPNDSO,SPNPSO,
     &                  SPNSD,SPNFC,SPSDFC,TRLEN,TRVEL,TQLEN,TQVEL,
CSPAS:23/5-11: second and third moment sum rules
C    &                  TRMAG,TRLON,BSRLON,EXENG,FONAC,FONA2,FOVIBG,
C    &                  SECMAT,
     &                  TRMAG,TRLON,BSRLON,TTLEN,EXENG,FONAC,FONA2,
     &                  FOVIBG,SECMAT,
CKeinSPASmehr
     &                  RMLEN,RQLEN,RLEN,RMVEL,RQVEL,RVEL,SLEN,SVEL,
     &                  WORK,LWORK)
C
C <<< Final results from ABACUS >>>
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <codata.h>
#include <aovec.h>
#include <ccinftap.h>
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0,
     *            D3 = 3.0D0, D100 = 100.0D0)
      CHARACTER NAME*6, SPDCAR*1, LAPRPC*8, mlab1*8, mlab2*8, mlab3*8
      DIMENSION PRAPOL(3,3), PRVPOL(6), CSTRA(*), SCTRA(*), POLDD(*),
     *          POLDL(*), POLDA(*), POLVL(*), POLVV(*)
      DIMENSION WORK(LWORK)
      DIMENSION SPNDSO(MXCOOR,MXCOOR), SPNPSO(MXCOOR,MXCOOR),
     &          SPNSD (MXCOOR,MXCOOR), SPNFC (MXCOOR,MXCOOR),
     &          SPSDFC(MXCOOR,MXCOOR)
CSPAS:23/5-11: second and third moment sum rules
C     DIMENSION TRLEN(*), TRVEL(*), TQLEN(*), TQVEL(*),
      DIMENSION TRLEN(*), TRVEL(*), TQLEN(*), TQVEL(*), TTLEN(*),
CKeinSPASmehr
     &          TRLON(*), TRMAG(*), BSRLON(*),
     &          EXENG(*), FONAC(*), FONA2(*), FOVIBG(*),
CSPAS : 6/11-2009 AO-SOPPA included
     &          SECMAT(*),
CKeinSPASmehr
     &          RMLEN(*), RQLEN(*), RLEN(*),
     &          RMVEL(*), RQVEL(*), RVEL(*), SLEN(*), SVEL(*)
      LOGICAL OD, FIRST
#include <abainf.h>
#include <taymol.h>
#include <nuclei.h>
#include <inftap.h>
#include <shells.h>
#include <symmet.h>
#include <helfey.h>
#include <moldip.h>
#include <sigma.h>
#include <relcor.h>
#include <aatens.h>
#include <pcmlog.h>
#include <difsec.h>
#include <expopt.h>
#include <primit.h>
CSPAS : 6/11-2009 AO-SOPPA included
#include <soppinf.h>
CKeinSPASmehr
C
C     Heading
C
      CALL TITLER('FINAL RESULTS from ABACUS','*',124)
C
C     Time and date
C
      CALL TSTAMP(' ',LUPRI)
      IF (NWNABA .NE. 0) THEN
         WRITE (LUPRI,'(/A)') ' WARNING: warnings have been issued,'//
     *      ' please check the output above !'
         NWNABA = 0
      END IF
C
      IF (LUSIFC .LE. 0) CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ',' ',
     &                               IDUMMY,.FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC) POTNUC,EMY,EACTIV,EMCSCF,ISTATE,ISPIN,NACTEL,
     *            LSYM
      READ (LUSIFC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT
C
      ERGMOL = EMCSCF
C
C     Geometry, Energy, Gradient, and Hessian output is
C     controlled by IPRDEF, such that it can been suppressed
C     in geometry optimizations.
C     For other properties the user has explicit asked for them
C     so they are always printed (changed Jan 2003 hjaaj).
C
      IF (IPRDEF .GE. 0) THEN
C
C     Geometry
C
         CALL HEADER('Molecular geometry (au)',-1)
         CALL PRIGEO(CORD)
C
C     Energy
C
         CALL HEADER('Molecular wave function and energy',-1)
         CALL MOLCHR(ICHRGE)
         WRITE (LUPRI,'(3(5X,A,I3))')
     *        'Spin multiplicity',ISPIN,
     *        'State number     ',ISTATE,
     *        'Total charge     ',ICHRGE
         WRITE (LUPRI,'(/T6,A,F17.10,A/T6,A,F17.8,A/T6,A,F17.4,A)')
     *        'Total energy     ',       EMCSCF,' au (Hartrees)',
     *        '                 ',  XTEV*EMCSCF,' eV',
     *        '                 ',XKJMOL*EMCSCF,' kJ/mol'
csonia 04/10/95
         IF (LUCME.GT.0) THEN
            WRITE (LUCME,'(A)') 'Molecular wave function and energy'
            WRITE (LUCME,'(3(5X,A,I3))')
     *        'Spin multiplicity',ISPIN,
     *        'State number     ',ISTATE,
     *        'Total charge     ',ICHRGE
            WRITE (LUCME,'(/T6,A,D24.15,A/T6,A,F17.8,A/T6,A,F17.4,A/)')
     *        'Total energy     ',       EMCSCF,' au (Hartrees)',
     *        '                 ',  XTEV*EMCSCF,' eV',
     *        '                 ',XKJMOL*EMCSCF,' kJ/mol'
         END IF
csonia 04/10/95
C
C     Relativistic corrections, K.Ruud, Dec.-93
C
         IF (DARWIN .OR. MASSVE) THEN
            CALL HEADER('Relativistic corrections',-1)
            IF (DARWIN) WRITE (LUPRI,'(5X,A,9X,F17.10,A)')
     &           'Darwin correction:            ',       DARWN,' au'
            IF (MASSVE) WRITE (LUPRI,'(5X,A,9X,F17.10,A)')
     &           'Mass-velocity correction:     ',       RMASSV,' au'
            IF (DARWIN .AND. MASSVE) THEN
               RELTOT = DARWN + RMASSV
               WRITE (LUPRI,'(/5X,A,F17.10,A,F6.4,A)')
     &              'Total relativistic correction:         ',  RELTOT,
     &              ' au (', D100*RELTOT/EMCSCF,'%)'
               WRITE (LUPRI,'(5X,A,F17.10,A//)')
     &              'Non-relativistic + relativistic energy:',
     &              EMCSCF + RELTOT,' au'
            END IF
         END IF
C
C     Orbital exponent gradient or molecular gradient
C
         IF (EXPGRD) THEN
            CALL HEADER(
     &         'Orbital-exponent gradient for all atoms (au)',-1)
            WRITE (LUPRI,'(4X,2A/)') ' #   atom Z   deg l       ',
     &          '        exponent           grad'//
     &          '           grad ln          grad ln'
            JCNT = 1
            ALPNRM = D0
            ALPNRL = D0
            DO I = 1, KMAX
               ISHELL = IPRSHL(I)
               LVAL = NHKT(ISHELL) - 1
               ICNT = NCENT(ISHELL)
               MLT = MULT(ISTBAO(ISHELL))
               ICHRG = IZATOM(ICNT)
               IF (JCNT .NE. ICNT) WRITE (LUPRI,'()')
               JCNT = ICNT
               GRAD = ALPGRD(I)/FLOAT(MLT)
               GRDL = GRAD*PRIEXP(I)
C              WRITE (LUPRI,'(5X,A,I3,3X,A,I2,I5,4X,A,3F17.10)')
               WRITE (LUPRI,
     &              '(A,I3,3X,A,I2,I5,2X,A,F23.10,2F17.10,1P,D17.8)')
     &              ' ##',I,NAMDEP(ICNT)(1:4),ICHRG,
     &               MLT,SPDCAR(LVAL),PRIEXP(I), GRAD,GRDL,GRDL
               ALPNRM = ALPNRM + FLOAT(MLT)*GRAD*GRAD 
               ALPNRL = ALPNRL + FLOAT(MLT)*GRDL*GRDL 
            END DO
            ALPNRM = SQRT(ALPNRM)
            ALPNRL = SQRT(ALPNRL)
            WRITE (LUPRI,'(/5X,A,19X,1P,D12.6,5X,D12.6)')
     &         '##     Norm of gradient:   ',ALPNRM,ALPNRL
C
            MAXZ = 0
            DO I = 1, KMAX
               MAXZ = MAX(MAXZ,IZATOM(NCENT(IPRSHL(I))))
            END DO
            DO K = 1, MAXZ
               FIRST = .TRUE.
               DO I = 1, KMAX
               IF (IZATOM(NCENT(IPRSHL(I))).EQ.K) THEN
C
                  IPRV = -1 
                  IF (I.GT.1) IPRV = NCENT(IPRSHL(I-1))
                  ICNT = NCENT(IPRSHL(I))
                  INXT = -1
                  IF (I.LT.KMAX) INXT = NCENT(IPRSHL(I+1))
C
                  IF (FIRST) THEN
                     IF (IPRV.NE.ICNT) IREF = I 
                     IF (INXT.NE.ICNT) FIRST = .FALSE.
                  ELSE
                     IF (IPRV.NE.ICNT) IOFF = I
                     J = IREF - IOFF + I
                     ALPGRD(J) = ALPGRD(J) + ALPGRD(I)
                     ALPGRD(I) = D0 
                  END IF
               END IF 
               END DO
            END DO
            CALL HEADER('Orbital-exponent gradient (au)',-1)
            WRITE (LUPRI,'(7X,A/)')  'Z    l            '//
     &       '    exponent                grad             grad ln'//
     &       '             grad ln'
            ALPNRM = D0
            ALPNRL = D0
            DO K = 1, MAXZ
               FIRST = .TRUE.
               DO I = 1, KMAX
               IF (IZATOM(NCENT(IPRSHL(I))).EQ.K) THEN
C
                  ICNT = NCENT(IPRSHL(I))
                  INXT = -1
                  IF (I.LT.KMAX) INXT = NCENT(IPRSHL(I+1))
C
                  IF (FIRST) THEN
                     ISHELL = IPRSHL(I)
                     LVAL = NHKT(ISHELL) - 1
                     GRAD = ALPGRD(I)
                     GRDL = GRAD*PRIEXP(I)
C                    WRITE (LUPRI,'(5X,A,I5,4X,A,3F20.10)')
                     WRITE (LUPRI,
     &                    '(1X,A,I5,4X,A,F24.10,2F20.10,1P,D20.8)')
     &                    '@@',K,SPDCAR(LVAL),PRIEXP(I), GRAD,GRDL,GRDL
                     ALPNRM = ALPNRM + GRAD*GRAD 
                     ALPNRL = ALPNRL + GRDL*GRDL 
                     IF (INXT.NE.ICNT) THEN
                        FIRST = .FALSE.
                        WRITE (LUPRI,'()')
                     END IF
                  END IF
               END IF 
               END DO
            END DO
            ALPNRM = SQRT(ALPNRM)
            ALPNRL = SQRT(ALPNRL)
            WRITE (LUPRI,'(5X,A,17X,1P,E12.6,8X,E12.6)')
     &         '@@    Norm of gradient:',ALPNRM,ALPNRL
C
         ELSE IF (MOLGRD) THEN
            CALL HEADER('Molecular gradient (au)',-1)
            CALL PRIGRD(GRDMOL,CSTRA,SCTRA)
         END IF
C
C     Molecular Hessian
C
C
C     *****************************************************************
C     ***** Translation invariance for Hellmann-Feynmann Hessians *****
C     *****************************************************************
C
         IF (HELFEY .AND. MOLHES) THEN
            DO I = 1, NUCDEP
               DO K = 1, 3
                  IPOS = 3*(I - 1) + K
                  HESMOL(IPOS,IPOS) = D0
                  DO J = 1, NUCDEP
                     IF (I .NE. J) THEN
                        JPOS = 3*(J - 1) + K
                        IF (I .LE. J) THEN
                           HESMOL(IPOS,IPOS) = HESMOL(IPOS,IPOS)
     &                                       - HESMOL(JPOS,IPOS)
                        ELSE
                           HESMOL(IPOS,IPOS) = HESMOL(IPOS,IPOS)
     &                                       - HESMOL(IPOS,JPOS)
                        END IF
                     END IF
                  END DO
               END DO
            END DO
            DO I = 1, 3*NUCDEP
               DO J = 1, I - 1
                  HESMOL(J,I) = HESMOL(I,J)
               END DO
            END DO
         END IF
C
         IF (MOLHES) THEN
C            IF (NFLOAT .EQ. 0) THEN
               CALL HEADER('Molecular Hessian (au)',-1)
               CALL PRIHES(HESMOL,'CENTERS',CSTRA,SCTRA)
C            ELSE
C               CALL HEADER('Hessian with respect to nuclear and '//
C     *              'orbital displacements (au)',-1)
C               CALL PRIHES(HESMOL,'CENTERS',CSTRA,SCTRA)
C               CALL HEADER('Molecular Hessian (au)',-1)
C               CALL PRIHES(HESFLT,'ATOMS',CSTRA,SCTRA)
C            END IF
         END IF
C
C     Dipole moments (disable printing for SOPPA)
Cdjw added SI units
C
         IF (.NOT. ABASOP) THEN
            CALL HEADER('Dipole moment',-1)
            DIPMOM = SQRT(DIP0(1)*DIP0(1) + DIP0(2)*DIP0(2) + 
     &           DIP0(3)*DIP0(3))
C
            WRITE (LUPRI,'(17X,A,15X,A,10X,A/3X,3F19.6)') 
     *             'au','Debye','C m (/(10**-30)',
     *             DIPMOM, DEBYE*DIPMOM, DIPSI*DIPMOM
            CALL HEADER('Dipole moment components',-1)
            CALL DP0PRI(DIP0)
C
         WRITE (LUPRI,'(2X,A,2X,A,F9.5,A)') ' Units:',
     *          ' 1 a.u. = ',DEBYE,' Debye '
         WRITE (LUPRI,'(11X,A,F9.5,A//)') ' 1 a.u. = ',DIPSI,
     *           ' (10**-30) C m (SI)'
C
            DO IPRPC = 1,3 
               IF (IPRPC .EQ. 1) LAPRPC = 'XDIPLEN '
               IF (IPRPC .EQ. 2) LAPRPC = 'YDIPLEN '
               IF (IPRPC .EQ. 3) LAPRPC = 'ZDIPLEN '
               CALL WRIPRO(DIP0(IPRPC),"SCF/DFT   ",1,
     *                     LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *                     0.0D0,0.0D0,0.0D0,1,0,0,0)
            END DO
Clf Printout of the LF dipole moment
            IF (PCM .AND. LOCFLD) THEN
               CALL HEADER('Local-field corrected dipole moment',-1)
               DIPLF = SQRT(DIPLF0(1)*DIPLF0(1) + DIPLF0(2)*DIPLF0(2) + 
     &              DIPLF0(3)*DIPLF0(3))
            WRITE (LUPRI,'(17X,A,15X,A,10X,A/3X,3F19.6)') 
     *             'au','Debye','C m (/(10**-30)',
     *             DIPLF, DEBYE*DIPLF, DIPSI*DIPLF
               CALL HEADER(
     $              'Local-field corrected dipole moment components',-1)
               CALL DP0PRI(DIPLF0)
            END IF
         END IF
      END IF
C     ... end if (iprdef .ge. 0)
C
C     APT population analysis
C
      IF (DIPDER) THEN
         IF (NFLOAT .EQ. 0) THEN
            CALL HEADER('Dipole moment gradient (au)',-1)
            CALL FCPRI(DIPFLT,'APT',CSTRA,SCTRA)
         ELSE
            CALL HEADER
     *           ('Unfolded dipole moment gradient (au) ',-1)
            CALL FCPRI(DIP1,'APT',CSTRA,SCTRA)
            CALL HEADER('Dipole moment gradient (APTs) (au)',-1)
            CALL FCPRI(DIPFLT,'APT',CSTRA,SCTRA)
         END IF
         ITEMP = MXCOOR*MXCOOR + 1
         CALL APTPOP(DIPFLT,HESMOL,QAPT,CSTRA,SCTRA,WORK,WORK(ITEMP),
     &               LWORK-ITEMP)
      END IF
C
C     Polarizabilities
C
      IF (POLAR) THEN
         IF (NFLOAT .EQ. 0) THEN
            CALL HEADER('Static polarizabilities (au)',-1)
            CALL POLPRI(POLARS,'   ',1)
            CALL HEADER('Static polarizabilities (angstroms**3)',-1)
            CALL POLPRI(POLARS,'EXP',1)
            DO 200 I = 1, 3
               DO 210 J = 1, 3
                  POLFLT(I,J) = POLARS(I,J)
 210           CONTINUE
 200        CONTINUE
         ELSE
            CALL HEADER('Unfolded static polarizabilities (au)',-1)
            CALL POLPRI(POLARS,'   ',1)
            CALL HEADER
     *         ('Unfolded static polarizabilities (angstroms**3)',-1)
            CALL POLPRI(POLARS,'EXP',1)
            CALL HEADER('Static polarizabilities (au)',-1)
            CALL POLPRI(POLFLT,'   ',1)
            CALL HEADER('Static polarizabilities (angstroms**3)',-1)
            CALL POLPRI(POLFLT,'EXP',1)
         END IF
C
C     Principal values and axes of polarizability
C
         CALL DIAPOL
      END IF
C
C     Frequency-dependent linear response properties
C     Polarizabilities and G tensors
C
      IF (ABALNR) CALL LNROUT(POLDD,POLDL,POLDA,POLVL,POLVV,FOVIBG,
     &                        IPRDEF,WORK,LWORK)
C
C     Excitation energies and (oriented) electronic circular dichroism
C
      IF (DOEXCI) THEN
         IF (AORPA .OR. DCRPA .OR. AOSOP .OR. AOSOC) THEN
CSPAS:23/5-11: second and third moment sum rules
C           CALL SO_EXCOUT(TRLEN,TRVEL,TRMAG,TRLON,BSRLON,EXENG,
C    &                     FONAC,FONA2,SECMAT,WORK,LWORK)
            CALL SO_EXCOUT(TRLEN,TRVEL,TQLEN,TQVEL,TRMAG,TRLON,BSRLON,
     &                     TTLEN,EXENG,FONAC,FONA2,SECMAT,WORK,LWORK)
CKeinSPASmehr
         ELSE
            CALL EXCOUT(TRLEN,TRVEL,TQLEN,TQVEL,
CSPAS:23/5-11: second and third moment sum rules
C    &                  TRMAG,TRLON,BSRLON,EXENG,
     &                  TRMAG,TRLON,BSRLON,TTLEN,EXENG,
CKeinSPASmehr
     &                  FONAC,FONA2,RMLEN,RQLEN,
     &                  RLEN,RMVEL,RQVEL,RVEL,
     &                  SLEN,SVEL,WORK,LWORK)
         END IF
      END IF
C
C     Molecular second moments
C
      IF (SECNDM) CALL SECRES(IPRDEF)
C
C     Molecular quadrupole moments
C
      IF (QUADRU) THEN
         CALL QDRRES(IPRDEF)
      END IF
C
C     Second moment gradient
C
      IF (QPGRAD) THEN
         CALL HEADER('Second moment gradient (au)',-1)
         CALL PRISEC(SEC1,'SECDER',CSTRA,SCTRA)
chs
         IF (DIPDER) CALL ATMDIP(CSTRA,SCTRA)
chs
      END IF
C
C     Molecular third moments (cartesian, not spherical octupole)
C
      IF (THIRDM) CALL THRDRES(IPRDEF)
C
C     Nuclear quadrupole moments
C
      IF (NQCC) THEN
         KAXIS = 1
         KLAST = KAXIS + MXCENT*9
         LWRK  = LWORK - KLAST + 1
         IF (KLAST .GT. LWORK) CALL
     &        STOPIT('ABARESULT','NQCRES',KLAST,LWORK)
         CALL NQCRES(IPRDEF,WORK(KAXIS))
      END IF
cLig >> added the call to SUSRES and SHIRES for CTOCD
C
C     CTOCD-DZ Susceptibilities and Nuclear Shieldings
C
      IF (IPRDEF .GE. 0) THEN
        IF (CTOCD) THEN
           IF(MAGSUS) THEN
              CALL SUSRES(IPRDEF)
           ENDIF
           IF(SHIELD) THEN
              DONS = .TRUE.
              CALL SHIRES(WORK,LWORK,IPRDEF)
              DONS = .FALSE.
              CALL SHIRES(WORK,LWORK,IPRDEF)
           ENDIF
           CTOCD = .FALSE.
        ENDIF
      ENDIF
cLig <<
C
C     Magnetic Susceptibilities
C
      IF (MAGSUS) THEN
         CALL SUSRES(IPRDEF)
      END IF
C
C     Molecular g-factor
C
      IF (MOLGFA) THEN
         JATOM  = NATOMS + NFLOAT
         KGEOM  = 1
         KAMASS = KGEOM + 3*JATOM
         KLAST  = KAMASS + JATOM
         LWRK   = LWORK - KLAST + 1
         IF (KLAST .GT. LWORK) CALL
     &        STOPIT('ABARESULT','MGFRES',KLAST,LWORK)
         CALL MGFRES(WORK(KGEOM),WORK(KAMASS),WORK(KLAST),LWRK,JATOM,
     &        IPRDEF)
      END IF
C
C     Atomic axial tensor
C
      IF (VCD) THEN
         CALL HEADER('Atomic axial tensors (AATs)',-1)
         CALL FCPRI(AATTOT,'AAT',CSTRA,SCTRA)
      END IF
C
C     Nuclear shieldings
C
      IF (SHIELD) CALL SHIRES(WORK,LWORK,IPRDEF)
C
C     Spin-rotation constants
C
      IF (SPINRO) THEN
         JATOM  = NATOMS + NFLOAT
         KGEOM  = 1
         KAMASS = KGEOM + 3*JATOM
         KGVAL  = KAMASS + JATOM
         KDIAMA = KGVAL + JATOM
         KPARMA = KDIAMA + 9*JATOM
         KLAST  = KPARMA + 9*JATOM
         LWRK   = LWORK - KLAST + 1
         IF (KLAST .GT. LWORK) CALL
     &        STOPIT('ABARESULT','SPRRES',KLAST,LWORK)
         CALL SPRRES(WORK(KGEOM),WORK(KAMASS),WORK(KGVAL),
     &               WORK(KDIAMA),WORK(KPARMA),WORK(KLAST),LWRK,
     &               JATOM,IPRDEF)
      END IF
C
C     Spin-spin-couplings constants
C
      IF (SPNSPN) THEN
         CALL SPIRES(SPNDSO,SPNPSO,SPNSD,SPNFC,SPSDFC,WORK,LWORK)
      END IF
      RETURN
 1000 FORMAT (1X,A6,F17.10,2F24.10)
      END
C  /* Deck readrs */
      SUBROUTINE ABAREADRS
C
C     Read ABACUS restart information
C     TUH 1985
C
C     Revised 12-Nov-1989 tuh
C
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <mxcent.h>
C
      CHARACTER*8 LABEL(4)
      LOGICAL AROUND
#include <abainf.h>
#include <past.h>
#include <energy.h>
#include <taymol.h>
#include <dipole.h>
#include <moldip.h>
#include <sigma.h>
#include <suscpt.h>
#include <spinro.h>
#include <molgfa.h>
#include <quadru.h>
#include <nqcc.h>
#include <aatens.h>
#include <gdvec.h>
C
      LUSTAR = -1
      CALL GPINQ('ABACUS.RESTART','EXIST',AROUND)
      IF (AROUND) THEN
         CALL GPOPEN(LUSTAR,'ABACUS.RESTART','OLD',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
      ELSE
         CALL GPOPEN(LUSTAR,'ABACUS.RESTART','NEW',' ','UNFORMATTED',
     &               IDUMMY,.FALSE.)
         GO TO 900
      END IF
      REWIND LUSTAR
      READ (LUSTAR) LABEL
      READ (LUSTAR) PASTWO, PASORT, PASDIP, PASONE, PASRES, PASREL,
     *              PASCRS, PASCZR, PASCTR, PASCRL, PASMAG, PASRTR,
     &              PASAAT, PASLRS, PASTRP, PASLNR, PASEXC, PASQPG
      READ (LUSTAR) NCOR
      READ (LUSTAR) ENERKE, ENERNA, ENEREE, ENERNN,
     *              (GRADKE(I), I=1,NCOR),
     *              (GRADNA(I), I=1,NCOR),
     *              (GRADEE(I), I=1,NCOR),
     *              (GRADNN(I), I=1,NCOR),
     *              (GRADFS(I), I=1,NCOR),
     *              ((HESMOL(I,J), I=1,NCOR), J=1,NCOR)
      READ (LUSTAR) DIPMN,  DIPME,
     *              ((DDIPN(I,J), I=1,3), J=1,NCOR),
     *              ((DDIPE(I,J), I=1,3), J=1,NCOR),
     *              ((DDIPS(I,J), I=1,3), J=1,NCOR),
     *              ((DDIPR(I,J), I=1,3), J=1,NCOR),
     *              POLARS
      READ (LUSTAR) ((SIGMAD(I,J), I=1,3), J=1,NCOR),
     &              ((SIGMAS(I,J), I=1,3), J=1,NCOR),
     &              ((SIGMAR(I,J), I=1,3), J=1,NCOR),
     &              ((SIGMAT(I,J), I=1,3), J=1,NCOR)
      READ (LUSTAR) ((SUSDIA(I,J), I=1,3), J=1,3), 
     &              ((SUS2EL(I,J), I=1,3), J=1,3),
     &              ((SUSDFT(I,J), I=1,3), J=1,3),
     &              ((SUSFS (I,J), I=1,3), J=1,3), 
     &              ((SUSREL(I,J), I=1,3), J=1,3),
     &              ((SUSTOT(I,J), I=1,3), J=1,3),
     &              ((SUSFSY(I,J), I=1,3), J=1,3)
C      READ (LUSTAR) ((SPNDSO(I,J), I=1,NCOR), J=1,NCOR),
C     &              ((SPNPSO(I,J), I=1,NCOR), J=1,NCOR),
C     &              ((SPNSD (I,J), I=1,NCOR), J=1,NCOR),
C     &              ((SPNFC (I,J), I=1,NCOR), J=1,NCOR),
C     &              ((SPSDFC(I,J), I=1,NCOR), J=1,NCOR)
      READ (LUSTAR)  (((TOTSPR(I,J,K), I=1,3), J=1,3), K=1,NCOR/3),
     &               (((SPRNUC(I,J,K), I=1,3), J=1,3), K=1,NCOR/3),
     &               ((ELSPRD(I,J), I=1,3), J=1, NCOR),
     &               ((ELSPRP(I,J), I=1,3), J=1, NCOR),
     &               (((SPRDNL(I,J,K), I=1,3), J=1,3), K=1,NCOR/3)
      READ (LUSTAR) ((TOTMGF(I,J), I=1,3), J=1,3),
     &              ((GFANUC(I,J), I=1,3), J=1,3),
     &              ((ELMGF(I,J), I=1,3), J=1,3),
     &              ((GFACDI(I,J), I=1,3), J=1,3)
      READ (LUSTAR) ((QDRNUC(I,J), I=1,3), J=1,3),
     &              ((QDREL(I,J), I=1,3), J=1,3),
     &              ((QDRTOT(I,J), I=1,3), J=1,3)
      READ (LUSTAR) (((UCNNQC(I,J,K), I=1,3), J=1,3), K=1,NCOR/3),
     &              (((ELNQC(I,J,K), I=1,3), J=1,3), K=1,NCOR/3)
      READ (LUSTAR) ((AATNUC(I,J), I=1,3), J=1,NCOR),
     &              ((AATORB(I,J), I=1,3), J=1,NCOR),
     &              ((AATCI (I,J), I=1,3), J=1,NCOR),
     &              ((AAT2ND(I,J), I=1,3), J=1,NCOR),
     &              ((AATTOT(I,J), I=1,3), J=1,NCOR)
      READ (LUSTAR) (IDORCT(I), I=1,80*MXCENT)
      READ (LUSTAR) ((IDORCI(I,J), I=1,24*(MXCENT + 1)), J=1,2)
      WRITE (LUPRI,'(4(A,1X),A8)')  ' Restart file used was created',
     *      LABEL(2)(1:2), LABEL(2)(3:5), LABEL(2)(6:7), LABEL(3)
      CALL GPCLOSE(LUSTAR,'KEEP')
      RETURN
C
  900 CONTINUE
      NWNABA = NWNABA + 1
      WRITE (LUPRI,'(//A/A/)')
     *   ' WARNING: No restart file found.',
     *   ' WARNING: Calculation proceeds without restart.'
      RETURN
      END
C  /* Deck abawritrs */
      SUBROUTINE ABAWRITRS
C
C     Write ABACUS restart information
C     tuh 1985
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
C
#include <past.h>
#include <energy.h>
#include <taymol.h>
#include <dipole.h>
#include <moldip.h>
#include <nuclei.h>
#include <sigma.h>
#include <suscpt.h>
#include <spinro.h>
#include <molgfa.h>
#include <quadru.h>
#include <nqcc.h>
#include <aatens.h>
#include <gdvec.h>
C
      LUSTAR = -1
      CALL GPOPEN(LUSTAR,'ABACUS.RESTART','UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
      NCOR = 3*NUCDEP
      REWIND LUSTAR
      CALL NEWLAB('ABARSTR',LUSTAR,LUPRI)
      WRITE (LUSTAR) PASTWO, PASORT, PASDIP, PASONE, PASRES, PASREL,
     *               PASCRS, PASCZR, PASCTR, PASCRL, PASMAG, PASRTR,
     &               PASAAT, PASLRS, PASTRP, PASLNR, PASEXC, PASQPG
      WRITE (LUSTAR) NCOR
      WRITE (LUSTAR) ENERKE, ENERNA, ENEREE, ENERNN,
     *               (GRADKE(I), I=1,NCOR),
     *               (GRADNA(I), I=1,NCOR),
     *               (GRADEE(I), I=1,NCOR),
     *               (GRADNN(I), I=1,NCOR),
     *               (GRADFS(I), I=1,NCOR),
     *               ((HESMOL(I,J), I=1,NCOR), J=1,NCOR)
      WRITE (LUSTAR) DIPMN,  DIPME,
     *               ((DDIPN(I,J), I=1,3), J=1,NCOR),
     *               ((DDIPE(I,J), I=1,3), J=1,NCOR),
     *               ((DDIPS(I,J), I=1,3), J=1,NCOR),
     *               ((DDIPR(I,J), I=1,3), J=1,NCOR),
     *               POLARS
      WRITE (LUSTAR) ((SIGMAD(I,J), I=1,3), J=1,NCOR),
     &               ((SIGMAS(I,J), I=1,3), J=1,NCOR),
     &               ((SIGMAR(I,J), I=1,3), J=1,NCOR),
     &               ((SIGMAT(I,J), I=1,3), J=1,NCOR)
      WRITE (LUSTAR) ((SUSDIA(I,J), I=1,3), J=1,3), 
     &               ((SUS2EL(I,J), I=1,3), J=1,3),
     &               ((SUSDFT(I,J), I=1,3), J=1,3),
     &               ((SUSFS (I,J), I=1,3), J=1,3), 
     &               ((SUSREL(I,J), I=1,3), J=1,3),
     &               ((SUSTOT(I,J), I=1,3), J=1,3),
     &               ((SUSFSY(I,J), I=1,3), J=1,3)
C      WRITE (LUSTAR) ((SPNDSO(I,J), I=1,NCOR), J=1,NCOR),
C     &               ((SPNPSO(I,J), I=1,NCOR), J=1,NCOR),
C     &               ((SPNSD (I,J), I=1,NCOR), J=1,NCOR),
C     &               ((SPNFC (I,J), I=1,NCOR), J=1,NCOR),
C     &               ((SPSDFC(I,J), I=1,NCOR), J=1,NCOR)
      WRITE(LUSTAR)  (((TOTSPR(I,J,K), I=1,3), J=1,3), K=1,NUCDEP),
     &               (((SPRNUC(I,J,K), I=1,3), J=1,3), K=1,NUCDEP),
     &               ((ELSPRD(I,J), I=1,3), J=1, NCOR),
     &               ((ELSPRP(I,J), I=1,3), J=1, NCOR),
     &               (((SPRDNL(I,J,K), I=1,3), J=1,3), K=1,NUCDEP)
      WRITE (LUSTAR) ((TOTMGF(I,J), I=1,3), J=1,3),
     &               ((GFANUC(I,J), I=1,3), J=1,3),
     &               ((ELMGF(I,J), I=1,3), J=1,3),
     &               ((GFACDI(I,J), I=1,3), J=1,3)
      WRITE (LUSTAR) ((QDRNUC(I,J), I=1,3), J=1,3),
     &               ((QDREL(I,J), I=1,3), J=1,3),
     &               ((QDRTOT(I,J), I=1,3), J=1,3)
      WRITE (LUSTAR) (((UCNNQC(I,J,K), I=1,3), J=1,3), K=1,NUCDEP),
     &               (((ELNQC(I,J,K), I=1,3), J=1,3), K=1,NUCDEP)
      WRITE (LUSTAR) ((AATNUC(I,J), I=1,3), J=1,NCOR),
     &               ((AATORB(I,J), I=1,3), J=1,NCOR),
     &               ((AATCI (I,J), I=1,3), J=1,NCOR),
     &               ((AAT2ND(I,J), I=1,3), J=1,NCOR),
     &               ((AATTOT(I,J), I=1,3), J=1,NCOR)
      WRITE (LUSTAR) (IDORCT(I), I=1,80*MXCENT)
      WRITE (LUSTAR) ((IDORCI(I,J), I=1,24*(MXCENT + 1)), J=1,2)
C
      CALL GPCLOSE(LUSTAR,'KEEP')
      RETURN
      END
C  /* Deck timsta */
      SUBROUTINE ABATIMSTA(TIMALL,TIMINP,TIMTRO,TIMONE,TIMMAG,TIMTEX,
     *                     TIMDRC,TIMORT,TIMDIP,TIMQPG,TIMRES,TIMREL,
     *                     TIMTST,TIMVIB,TIMWLK,
     *                     TIMLRS,TIMTRP,TIMEXC,TIMAAT,TIMLNR,TIMCTO)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <infpri.h>
#include <abainf.h>
      CALL HEADER('CPU time statistics for ABACUS',0)
      CALL TIMPRI('INPUT ',TIMINP,TIMALL)
      CALL TIMPRI('TRAROT',TIMTRO,TIMALL)
      CALL TIMPRI('ONEINT',TIMONE,TIMALL)
      CALL TIMPRI('MAGINT',TIMMAG,TIMALL)
      CALL TIMPRI('TWOTEX',TIMTEX,TIMALL)
      CALL TIMPRI('DRCCTL',TIMDRC,TIMALL)
      CALL TIMPRI('RHSIDE',TIMORT,TIMALL)
      CALL TIMPRI('DIPOLE',TIMDIP,TIMALL)
      CALL TIMPRI('QPGRAD',TIMQPG,TIMALL)
      CALL TIMPRI('RESPON',TIMRES,TIMALL)
      CALL TIMPRI('EXCITA',TIMEXC,TIMALL)
      CALL TIMPRI('AATDRV',TIMAAT,TIMALL)
      CALL TIMPRI('LNRABA',TIMLNR,TIMALL)
      CALL TIMPRI('LINRES',TIMLRS,TIMALL)
      CALL TIMPRI('TRP LR',TIMTRP,TIMALL)
      CALL TIMPRI('RELAX ',TIMREL,TIMALL)
      CALL TIMPRI('WALK  ',TIMWLK,TIMALL)
      CALL TIMPRI('VIBCTL',TIMVIB,TIMALL)
      CALL TIMPRI('CTOCD ',TIMCTO,TIMALL)
      CALL TIMPRI('REST  ',TIMWLK,TIMALL)
      WRITE (LUPRI,'()')
      CALL TIMPRI('TOTAL ',TIMALL,TIMALL)
      WRITE (LUPRI,'(/)')
      RETURN
      END
C  /* Deck abadrc */
      SUBROUTINE ABADRC
C Feb 90 hjaaj -- temp. replacement for DRCCTL
C 900216: call of abadrc is disabled by *IF DEF,ABADRC
C
      CALL QUIT('ABACUS error: ABADRC called')
      RETURN
      END
C  /* Deck qdrres */
      SUBROUTINE QDRRES(IPRINT)
C Feb 05 djw - principal axis displays included.
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <codata.h>
      PARAMETER (DM1 = -1.0D0, D1 = 1.0D0,
     &           CGS = (ECHARGE*XTANG**2*CCM*1.D11),
     &           SI  = ECHARGE*XTANGM10**2*1.0D+40)
      DIMENSION QUADRT(3,3) 
      DIMENSION WRK(3), IWRK(3), AXES(3,3), PVAL(6)
#include <symmet.h>
#include <quadru.h>
#include <abainf.h>
      CHARACTER*6 ABC(3)
      DATA ABC /'Q_aa','Q_bb','Q_cc'/ 
C
      CALL TITLER('ABACUS - Molecular quadrupole moments','*',113)
      CALL DCOPY(9,QDREL,1,QUADRT,1)
      CALL DSCAL(9,DM1,QUADRT,1)
C
C     Principal values
C     ================
C
      CALL DUNIT(AXES,3)
      IJ = 1
      DO 100 I = 1, 3
         DO 110 J = 1, I
            PVAL (IJ) = QUADRT(I,J) + QDRNUC(I,J)
            IJ = IJ + 1
 110     CONTINUE
 100  CONTINUE
      CALL JACO(PVAL,AXES,3,3,3,WRK,IWRK)
      PVAL(1) = PVAL(1)
      PVAL(2) = PVAL(3)
      PVAL(3) = PVAL(6)
      CALL ORDER2(AXES,PVAL,3,3)
C
      WRITE (LUPRI,'(2(/A))')
     &      '  Units:  C m**2/(10**(-40))   (SI) ',
     &      '          Debye*AAngstroem     (cgs)'
C
      CALL HEADER('Principal values (au, SI and cgs) and axes:',1)
      DO 200 I = 1, 3
         WRITE (LUPRI,'(2X,A,1X,3F12.4,2X,3F10.4)')
     &          ABC(I),PVAL(I),SI*PVAL(I),CGS*PVAL(I),
     &          (AXES(IPTAX(J,2),I),J=1,3)
 200  CONTINUE
C
      IF (IPRINT .GE. 2) THEN
         CALL HEADER('Nuclear contribution (au)',-1)
         CALL POLPRI(QDRNUC,'   ',1)
         CALL HEADER('Electronic contribution (au)',-1)
         CALL POLPRI(QUADRT,'   ',1)
      END IF
      CALL DAXPY(9,D1,QDRNUC,1,QUADRT,1)
      CALL HEADER('Total quadrupole moments (au)',-1)
      CALL POLPRI(QUADRT,'   ',1)
      CALL HEADER('Total quadrupole moments (cgs, Debye*AAngstroem)',-1)
      CALL POLPRI(QUADRT,'CGS',1)
      CALL HEADER('Total quadrupole moments (SI)',-1)
      CALL POLPRI(QUADRT,'SIU',1)
      RETURN
      END
C  /* Deck secres */
      SUBROUTINE SECRES(IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (DM1 = -1.0D0, D1 = 1.0D0)
      DIMENSION SECNDT(3,3)
#include <symmet.h>
#include <secmom.h>
#include <abainf.h>
C
      CALL TITLER('ABACUS - Molecular second order moments','*',113)
      CALL DCOPY(9,SCMEL,1,SECNDT,1)
      CALL DSCAL(9,DM1,SECNDT,1)
      IF (IPRINT .GE. 2) THEN
         CALL HEADER('Nuclear contribution (au)',-1)
         CALL POLPRI(SCMNUC,'   ',1)
         CALL HEADER('Electronic contribution (au)',-1)
         CALL POLPRI(SECNDT,'   ',1)
      END IF
      CALL DAXPY(9,D1,SCMNUC,1,SECNDT,1)
      CALL HEADER('Total second order moments (au)',-1)
      CALL POLPRI(SECNDT,'   ',1)
      CALL HEADER('Total second order moments (cgs, Debye*AAngstroem)',
     &            -1)
      CALL POLPRI(SECNDT,'CGS',1)
      RETURN
      END
C  /* Deck octres */
      SUBROUTINE THRDRES(IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
      PARAMETER (DM1 = -1.0D0, D1 = 1.0D0)
      DIMENSION SECNDT(3,3)
#include <symmet.h>
#include <thimom.h>
#include <abainf.h>
C
      CALL TITLER('ABACUS - Molecular third order moments','*',113)
      CALL DCOPY(27,THDMEL,1,THDTOT,1)
      CALL DSCAL(27,DM1,THDTOT,1)
      IF (IPRINT .GE. 2) THEN
         CALL HEADER('Nuclear contribution (au)',-1)
         CALL PRIOCT(THDNUC)
         CALL HEADER('Electronic contribution (au)',-1)
         CALL PRIOCT(THDTOT)
      END IF
      CALL DAXPY(27,D1,THDNUC,1,THDTOT,1)
      CALL HEADER('Total third order moments (au)',-1)
      CALL PRIOCT(THDTOT)
c     CALL HEADER('Total third order moments (cgs, Debye*AAngstroem)',
c    &             -1)
c     CALL PRIOCT(THDTOT,'CGS',1)
      RETURN
      END
C
C  /* Deck susres */
      SUBROUTINE SUSRES(IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
#include <dftcom.h>
      PARAMETER (D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0, THRESH = 1.0D-5,
     &           SI = (XTANG*XTANG*ECHARGE*ECHARGE/EMASS)*1.0D-21,
     &           CGS = SI*XFMOL*1.0D6)
      LOGICAL D12, D23
      DIMENSION WRK(3), IWRK(3), AXES(3,3), PVAL(6), DIASUS(3,3),
     &          PARSUS(3,3), PVAL1(6), PVAL2(6)
#include <symmet.h>
#include <suscpt.h>
#include <abainf.h>
#include <inftap.h>
C
CTOCD
#include <ctocdcc.h> 
CTOCD
C
      CHARACTER*6 ABC(3)
      CHARACTER*8 LAPRPC
      DATA ABC /'khi_aa','khi_bb','khi_cc'/

C
      IF (CTOMAG) THEN
         IF (ALCCS) THEN
            CALL TITLER('CCS CTOCD-DZ MAGNETIZABILITIES','*',122)
         ELSE IF (ALCC2) THEN
            CALL TITLER('CC2 CTOCD-DZ MAGNETIZABILITIES','*',122)
         ELSE IF (ALCC3) THEN
            CALL TITLER('CC3 CTOCD-DZ MAGNETIZABILITIES','*',122)
         ELSE IF (ALCCSD) THEN
            CALL TITLER('CCSD CTOCD-DZ MAGNETIZABILITIES','*',122)
         END IF
      ELSE
         CALL TITLER('ABACUS - MAGNETIZABILITIES','*',124)
      END IF
C

cLig <> added the .NOT. CTOCD
      IF ((.NOT.NOLOND).AND.(.NOT.CTOCD) .AND. (.NOT.CTOMAG)) THEN
C
         CALL DZERO(SUSTOT,9)
         CALL DZERO(DIASUS,9)
         CALL DZERO(PARSUS,9)
         DO 100 I = 1, 3
         DO 100 J = 1, 3
            DIASUS(I,J) = DIASUS(I,J) + SUSCOM(I,J)
            SUSTOT(I,J) = SUSTOT(I,J) + SUSREL(I,J)
     &                                + SUSDIA(I,J)
     &                                + SUS2EL(I,J)
     &                                + SUSFS (I,J)
     &                                + SUSFSY(I,J)
            IF (DFTRUN) THEN
               SUSTOT(I,J) = SUSTOT(I,J) + SUSDFT(I,J)
            END IF
            PARSUS(I,J) = SUSTOT(I,J) - DIASUS(I,J)
  100    CONTINUE
C
         IF (IPRINT .GT. 2) THEN
            CALL HEADER('One-electron expectation values',-1)
            CALL POLPRI(SUSDIA,'   ',-2)
            AVERAG = (SUSDIA(1,1) + SUSDIA(2,2) + SUSDIA(3,3))/D3
            WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
C
            CALL HEADER('Two-electron expectation values',-1)
            CALL POLPRI(SUS2EL,'   ',-2)
            AVERAG = (SUS2EL(1,1) + SUS2EL(2,2) + SUS2EL(3,3))/D3
            WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
C
            IF (DFTRUN) THEN
               CALL HEADER('DFT static contribution',-1)
               CALL POLPRI(SUSDFT,'   ',-2)
               AVERAG = (SUSDFT(1,1) + SUSDFT(2,2) + SUSDFT(3,3))/D3
               WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
            END IF
C
            CALL HEADER('Highest-order reorthonormalization',-1)
            CALL POLPRI(SUSFS,'   ',-2)
            AVERAG = (SUSFS (1,1) + SUSFS (2,2) + SUSFS (3,3))/D3
            WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
C
            CALL HEADER('Lowest-order reorthonormalization',-1)
            CALL POLPRI(SUSFSY,'   ',-2)
            AVERAG = (SUSFSY(1,1) + SUSFSY(2,2) + SUSFSY(3,3))/D3
            WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
C
            CALL HEADER('Relaxation',-1)
            CALL POLPRI(SUSREL,'   ',-2)
            AVERAG = (SUSREL(1,1) + SUSREL(2,2) + SUSREL(3,3))/D3
            WRITE (LUPRI,'(/,6X,A,F12.6)') ' Average value:',AVERAG
C
         END IF
cLig >> what to do for CTOCD
      ELSE IF(CTOCD .OR. CTOMAG) THEN
cDEBUG
c 	write(LUPRI,*) 'Dentro SUSRES'
c	write(LUPRI,*) 'SUSDZD= ',SUSDZD
        CALL DZERO(SUSTOT,9)
        CALL DZERO(DIASUS,9)
        CALL DZERO(PARSUS,9)
        DO 110 I = 1, 3
        DO 110 J = 1, 3
          DIASUS(I,J) = SUSDZD(I,J)
          PARSUS(I,J) = SUSREL(I,J)
          SUSTOT(I,J) = DIASUS(I,J) + PARSUS(I,J)
 110    CONTINUE
cLig <<
      ELSE
C
         DO 200 J = 1, 3
         DO 200 I = 1, 3
            DIASUS(I,J) = SUSDIA(I,J)
            PARSUS(I,J) = SUSREL(I,J)
            SUSTOT(I,J) = SUSDIA(I,J) + SUSREL(I,J)
  200    CONTINUE
C
      END IF
C
C     Principal values
C     ================
C
      CALL DUNIT(AXES,3)
      IJ = 1
      DO 300 I = 1, 3
         DO 310 J = 1, I
            PVAL1(IJ) = DIASUS(I,J)
            PVAL2(IJ) = PARSUS(I,J)
            PVAL (IJ) = SUSTOT(I,J)
            IJ = IJ + 1
 310     CONTINUE
 300  CONTINUE
      CALL JACO(PVAL,AXES,3,3,3,WRK,IWRK)
      PVAL(1) = -PVAL(1)
      PVAL(2) = -PVAL(3)
      PVAL(3) = -PVAL(6)
      CALL ORDER2(AXES,PVAL,3,3)
      D12 = (ABS(PVAL(1)-PVAL(2)) .LT. THRESH)
      D23 = (ABS(PVAL(2)-PVAL(3)) .LT. THRESH)
C
      CALL JACO(PVAL1,DUMMY,3,3,0,WRK,IWRK)
      CALL JACO(PVAL2,DUMMY,3,3,0,WRK,IWRK)
      DIAMAG = -(PVAL1(1) + PVAL1(3) + PVAL1(6))/D3
      PARAMA = -(PVAL2(1) + PVAL2(3) + PVAL2(6))/D3
C
      WRITE (LUPRI,'(2(/1X,A))')
     &      ' Units:  JT**(-2)/(10**(-30))   (SI) ',
     &      '         ppm cm**(3) mol**(-1)  (cgs)'
      IF (NOLOND) THEN
         WRITE (LUPRI,'(/A)')
     &   ' Calculation without London orbitals (GIAOs).'
      ELSE
         IF (.NOT. (CTOCD .OR. CTOMAG))
     &      WRITE (LUPRI,'(/A)')' London orbitals (GIAOs) used.'
      END IF
cLig >> added a print statement
      IF (CTOCD .OR. CTOMAG) THEN
        WRITE (LUPRI,'(/,1X,A)')'CTOCD-DZ results:'   
      ENDIF 
cLig <<
      IF (D12 .AND. D23) THEN
         AVE = (PVAL(1) + PVAL(2) + PVAL(3))/D3
         WRITE (LUPRI,'(3(/,1X,A,3(F12.4,A)))')
     &      ' Magnetizability:           ',    AVE,' (au), ',
     &        SI*AVE*1.0D31,' (SI), ', CGS*AVE,' (cgs)',
     &      ' Diamagnetic contribution:  ',DIAMAG,' (au), ',
     &        SI*DIAMAG*1.0D31,' (SI), ',CGS*DIAMAG,' (cgs)',
     &      ' Paramagnetic contribution: ',PARAMA,' (au), ',
     &        SI*PARAMA*1.0D31,' (SI), ',CGS*PARAMA,' (cgs)'
         WRITE (LUPRI,'(/,1X,A)') ' Magnetizability is spherical.'
      ELSE IF (D12 .OR. D23) THEN
         AVE = (PVAL(1) + PVAL(2) + PVAL(3))/D3
         IF (D12) THEN
            PAR = PVAL(3)
            PER = (PVAL(1) + PVAL(2))/D2
         ELSE
            PAR = PVAL(1)
            PER = (PVAL(2) + PVAL(3))/D2
         END IF
         ANI = PAR - PER
         CALL HEADER('Magnetizabilities           '
     &             //'         au          SI         cgs',1)
         WRITE (LUPRI,'(6(1X,A,3F12.4,/))')
     &      ' Isotropic magnetizability: ',AVE,SI*AVE*1.0D31,CGS*AVE,
     &      ' Diamagnetic contribution:  ',DIAMAG,SI*DIAMAG*1.0D31,
     &                                     CGS*DIAMAG,
     &      ' Paramagnetic contribution: ',PARAMA,SI*PARAMA*1.0D31,
     &                                     CGS*PARAMA,
     &      ' Parallel component:        ',PAR,SI*PAR*1.0D31,CGS*PAR,
     &      ' Perpendicular component:   ',PER,SI*PER*1.0D31,CGS*PER,
     &      ' Anisotropy:                ',ANI,SI*ANI*1.0D31,CGS*ANI
         WRITE (LUPRI,'(1X,A)') ' Magnetizability is cylindrical.'
         CALL HEADER('Principal values (au, SI, and cgs) and axes:',1)
         DO 700 I = 1, 3
            WRITE (LUPRI,'(2X,A,1X,3F12.4,2X,3F10.4)')
     &         ABC(I),PVAL(I),SI*PVAL(I)*1.0D31,CGS*PVAL(I),
     &         (AXES(IPTAX(J,2),I),J=1,3)
 700     CONTINUE
         LAPRPC = 'ISO_MAGNE'
         CALL WRIPRO(SI*AVE*1.0D31,"SCF/DFT   ",403,
     *               LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *               0.0D0,0.0D0,0.0D0,1,0,0,0)
         LAPRPC = 'ANI_MAGNE'
         CALL WRIPRO(SI*ANI*1.0D31,"SCF/DFT   ",403,
     *               LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *               0.0D0,0.0D0,0.0D0,1,0,0,0)
         LAPRPC = 'DIA_MAGNE'
         CALL WRIPRO(SI*DIAMAG*1.0D31,"SCF/DFT   ",403,
     *               LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *               0.0D0,0.0D0,0.0D0,1,0,0,0)
         LAPRPC = 'PAR_MAGNE'
         CALL WRIPRO(SI*PARAMA*1.0D31,"SCF/DFT   ",403,
     *               LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *               0.0D0,0.0D0,0.0D0,1,0,0,0)
      ELSE
         CALL DSWAP(1,PVAL(2),1,PVAL(1),1)
         CALL DSWAP(3,AXES(1,2),1,AXES(1,1),1)
         AVE  = (PVAL(1) + PVAL(2) + PVAL(3))/D3
         ANI1 = PVAL(3) - (PVAL(1) + PVAL(2))/D2
         ANI2 = PVAL(2) - (PVAL(1) + PVAL(3))/D2
         CALL HEADER('Magnetizabilities           '
     &             //'         au          SI         cgs',1)
         WRITE (LUPRI,'(5(1X,A,3F12.4,/))')
     &      ' Isotropic magnetizability: ',AVE, SI*AVE*1.0D31, CGS*AVE,
     &      ' Diamagnetic contribution:  ',DIAMAG,SI*DIAMAG*1.0D31,
     &                                     CGS*DIAMAG,
     &      ' Paramagnetic contribution: ',PARAMA,SI*PARAMA*1.0D31,
     &                                     CGS*PARAMA,
     &      ' 1st anisotropy:            ',ANI1,SI*ANI1*1.0D31,CGS*ANI1,
     &      ' 2nd anisotropy:            ',ANI2,SI*ANI2*1.0D31,CGS*ANI2
         CALL HEADER('Principal values (au, SI, and cgs) and axes:',1)
         DO 800 I = 1, 3
            WRITE (LUPRI,'(2X,A,1X,3F12.4,2X,3F10.4)')
     &         ABC(I),PVAL(I),SI*PVAL(I)*1.0D31,CGS*PVAL(I),
     &         (AXES(IPTAX(J,2),I),J=1,3)
 800     CONTINUE
         LAPRPC = 'ISO_MAGNE'
         CALL WRIPRO(SI*AVE*1.0D31,"SCF/DFT   ",403,
     *               LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *               0.0D0,0.0D0,0.0D0,1,0,0,0)
         LAPRPC = 'ANI1MAGNE'
         CALL WRIPRO(SI*ANI1*1.0D31,"SCF/DFT   ",403,
     *               LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *               0.0D0,0.0D0,0.0D0,1,0,0,0)
         LAPRPC = 'ANI2MAGNE'
         CALL WRIPRO(SI*ANI2*1.0D31,"SCF/DFT   ",403,
     *               LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *               0.0D0,0.0D0,0.0D0,1,0,0,0)
         LAPRPC = 'DIA_MAGNE'
         CALL WRIPRO(SI*DIAMAG*1.0D31,"SCF/DFT   ",403,
     *               LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *               0.0D0,0.0D0,0.0D0,1,0,0,0)
         LAPRPC = 'PAR_MAGNE'
         CALL WRIPRO(SI*PARAMA*1.0D31,"SCF/DFT   ",403,
     *               LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *               0.0D0,0.0D0,0.0D0,1,0,0,0)
      END IF
C
      CALL HEADER('Total magnetizability tensor (au)',-1)
csonia 04/10/1995
      IF (LUCME.GT.0)
     &  WRITE (LUCME,'(A)') 'Total magnetizability tensor (au)'
csonia 04/10/1995                                                              
      CALL DSCAL(9,-D1,SUSTOT,1)
      CALL POLPRI(SUSTOT,'   ',-2)
      CALL DSCAL(9,-D1,SUSTOT,1)
C
Cspas: also print diamagnetic and paramagnetic tensors
C
      CALL HEADER('Diamagnetic magnetizability tensor (au)',-1)
      CALL DSCAL(9,-D1,DIASUS,1)
      CALL POLPRI(DIASUS,'   ',-2)
      CALL DSCAL(9,-D1,DIASUS,1)
C
      CALL HEADER('Paramagnetic magnetizability tensor (au)',-1)
      CALL DSCAL(9,-D1,PARSUS,1)
      CALL POLPRI(PARSUS,'   ',-2)
      CALL DSCAL(9,-D1,PARSUS,1)
C
      RETURN
      END
C  /* Deck nqcres */
      SUBROUTINE NQCRES(IPRINT,AXIS)
C
C     Output routine for nuclear quadrupole coupling constants
C     Based on O.Christiansens output routines, K.Ruud, Nov.-94
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <codata.h>
C
      PARAMETER (XMUB =  ECHARGE*HBAR/(2.0D0*EMASS))
      PARAMETER (AUFGR = XTJ*1D10*1D10/(ECHARGE*XTANG**2))
      PARAMETER (CONV = ECHARGE*AUFGR*1D-28*1D-6/(HBAR*2.0*PI),
     &           CONV2 = 1.0D-10*2.0D0*XMUB/(2.0D0*PI*HBAR))
      DIMENSION PVAL(6), AXIS(3,3,MXCENT), WRK(3), IWRK(3)
      CHARACTER*33 TEXT
cmbh string variable
      CHARACTER*8 laname 
cmbh end
#include <nuclei.h>
#include <symmet.h>
#include <qm3.h>
#include <nqcc.h>
#include <ibtfun.h>
C
      IF (IPRINT .GE. 0) THEN
         CALL TITLER('ABACUS - Nuclear Quadrupole moments','*',116)
         IF (IPRINT .GT. 3) THEN
            TEXT = 'Nuclear EFG Tensor of nucleus    '
            CALL EFGPRI(UCNNQC,TEXT)
            TEXT = 'Electronic EFG Tensor of nucleus '
            CALL EFGPRI(ELNQC,TEXT)
         END IF
      END IF
C
C     Add nuclear and electronic contribution. Keep in electronic part
C
      DO 70 I = 1, NUCDEP
chj-s-090516 : defined ISUBSI, NSISY(0) for .NOT.QM3 so test is passed
c       IF (( (QM3) .AND. (ISUBSY(I) .EQ. 0) .AND.
c    &     (ISUBSI(I) .LE. NSISY(0)) ) .OR. (.NOT. QM3)) THEN
        IF ( (ISUBSY(I) .EQ. 0) .AND. (ISUBSI(I) .LE. NSISY(0)) ) THEN
          DO 60 J = 1, 3
            DO 50 K = 1, 3
               ELNQC(J,K,I) = ELNQC(J,K,I) + UCNNQC(J,K,I)
 50         CONTINUE
 60       CONTINUE
        END IF
 70   CONTINUE
      IF (IPRINT .GE. 2) THEN
         TEXT = 'Total non-diagonalized EFG Tensor'
         CALL EFGPRI(ELNQC,TEXT)
      END IF
cmbh print EFG tensor for midas here
      do iatom=1,nucdep
          laname=NAMDEP(iatom)//'  '
          call stripblanks(laname)
          call wripro(ELNQC(1,1,IATOM),"SCF/DFT   ",502,
     *            ' XEFGCAR',' XEFGCAR',laname,laname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
          call wripro(ELNQC(1,2,IATOM),"SCF/DFT   ",502,
     *            ' XEFGCAR',' YEFGCAR',laname,laname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
          call wripro(ELNQC(1,3,IATOM),"SCF/DFT   ",502,
     *            ' XEFGCAR',' ZEFGCAR',laname,laname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
          call wripro(ELNQC(2,1,IATOM),"SCF/DFT   ",502,
     *            ' YEFGCAR',' XEFGCAR',laname,laname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
          call wripro(ELNQC(2,2,IATOM),"SCF/DFT   ",502,
     *            ' YEFGCAR',' YEFGCAR',laname,laname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
          call wripro(ELNQC(2,3,IATOM),"SCF/DFT   ",502,
     *            ' YEFGCAR',' ZEFGCAR',laname,laname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
          call wripro(ELNQC(3,1,IATOM),"SCF/DFT   ",502,
     *            ' ZEFGCAR',' XEFGCAR',laname,laname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
          call wripro(ELNQC(3,2,IATOM),"SCF/DFT   ",502,
     *            ' ZEFGCAR',' YEFGCAR',laname,laname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
          call wripro(ELNQC(3,3,IATOM),"SCF/DFT   ",502,
     *            ' ZEFGCAR',' ZEFGCAR',laname,laname,
     *            0.0D0,0.0D0,0.0D0,1,0,0,0)
      enddo

cmbh end
C
C     Diagonalize each atomic EFG tensor
C
      DO 10 IATOM = 1, NUCDEP
        IF( (ISUBSY(I) .EQ. 0) .AND. (ISUBSI(I) .LE. NSISY(0)) ) THEN
          CALL DUNIT(AXIS(1,1,IATOM),3)
          IJ = 1
          DO 20 I = 1, 3
            DO 30 J = 1, I
              PVAL(IJ) = ELNQC(I,J,IATOM)
              IJ = IJ + 1
 30         CONTINUE
 20       CONTINUE
          CALL JACO(PVAL,AXIS(1,1,IATOM),3,3,3,WRK,IWRK)
          DIANQC(1,IATOM) = PVAL(1)
          DIANQC(2,IATOM) = PVAL(3)
          DIANQC(3,IATOM) = PVAL(6)
        END IF
 10   CONTINUE
      IF (IPRINT .GE. 2) THEN
         TEXT = 'EFG Principal axis               '
         CALL EFGPRI(AXIS,TEXT)
      END IF
C
C     Calculate the nuclear quadrupole moments
C
      NATOM = 0
      DO 100 IATOM = 1, NUCIND
        IF( (ISUBSY(I) .EQ. 0) .AND. (ISUBSI(I) .LE. NSISY(0)) ) THEN
          DO 110 ISYM = 0, MAXOPR
            IF (IBTAND(ISTBNU(IATOM),ISYM) .EQ. 0) THEN
              NATOM = NATOM + 1
              NUCCHA = IZATOM(IATOM)
              IF (NUCCHA  .EQ. 0) GOTO 100
              IF (IPRINT .GE. 0) CALL HEADER('Nuclear quadrupole '//
     &                         'moment for '//NAMDEP(NATOM),-1)
              DO 120 ISO = 1, 5
                QMOM = DISOTP(NUCCHA,ISO,'QMOM')
                IF (QMOM .NE. 0.D0) THEN
                  SPIN = DISOTP(NUCCHA,ISO,'SPIN')
                  MASS = NINT(DISOTP(NUCCHA,ISO,'A'))
                  FACT = CONV*QMOM
                  FACT2 = FACT/CONV2
                  IF (IPRINT .GE. 0) THEN
                    WRITE (LUPRI,'(10X,A10,F5.1,10X,A18,F10.8)')
     &                   'Spin     :',SPIN,'Quadrupole moment:',QMOM
                    WRITE (LUPRI,'(28X,A2,11X,A3,8X,A5)')
     &                   'au','MHz','Gauss'
                  END IF
                  IZ = IDAMAX(3,DIANQC(1,NATOM),1)
                  CALL DSWAP(1,DIANQC(3,NATOM),1,DIANQC(IZ,NATOM),1)
                  IX = IDAMIN(3,DIANQC(1,NATOM),1)
                  CALL DSWAP(1,DIANQC(1,NATOM),1,DIANQC(IX,NATOM),1)
                  IF (IPRINT .GE. 0) THEN
                    WRITE(LUPRI,'(1X,A6,I3,2X,A11,F10.6,3X,F10.6,3X,
     &                            F10.6)')
     &              NAMDEP(NATOM),MASS,'Vx =     ',
     &              DIANQC(1,NATOM),DIANQC(1,NATOM)*FACT,
     &              DIANQC(1,NATOM)*FACT2
                    WRITE(LUPRI,'(1X,A6,I3,2X,A11,F10.6,3X,F10.6,3X,
     &                            F10.6)')
     &              NAMDEP(NATOM),MASS,'Vy =     ',
     &              DIANQC(2,NATOM),DIANQC(2,NATOM)*FACT,
     &              DIANQC(2,NATOM)*FACT2
                    WRITE(LUPRI,'(1X,A6,I3,2X,A11,F10.6,3X,F10.6,3X,
     &                            F10.6)')
     &              NAMDEP(NATOM),MASS,'Vz =     ',
     &              DIANQC(3,NATOM),DIANQC(3,NATOM)*FACT,
     &              DIANQC(3,NATOM)*FACT2
                   END IF
                   IZ = IDAMAX(3,DIANQC(1,NATOM),1)
                   CALL DSWAP(1,DIANQC(3,NATOM),1,DIANQC(IZ,NATOM),1)
                   IX = IDAMIN(3,DIANQC(1,NATOM),1)
                   CALL DSWAP(1,DIANQC(1,NATOM),1,DIANQC(IX,NATOM),1)
                   ETA = (DIANQC(1,NATOM) - DIANQC(2,NATOM))/
     &                    DIANQC(3,NATOM)
                   IF (IPRINT .GE. 0) THEN
                     WRITE (LUPRI,'(1X,A39,F14.6)')
     &                     'Asymmetry constant eta:', ETA
                     IF (SPIN .NE. 0.D0.AND.SPIN.NE.0.5) THEN
                       AUXILM = 3.D0 * (2.D0*SPIN + 3.D0)*QMOM*QMOM*
     &                        ( 1.D0 + ETA*ETA/3.D0 ) *
     &                 DIANQC(3,NATOM)*DIANQC(3,NATOM)/
     &                        ( 40.D0 * SPIN*SPIN *(2.D0*SPIN - 1.D0))
                       WRITE (LUPRI,'(1X,A39,D14.6)')
     &                        'Auxiliary value au  1/T', AUXILM
                     END IF
                   END IF
                 END IF
 120           CONTINUE
             END IF
 110       CONTINUE
         END IF
 100  CONTINUE
      RETURN
      END
C  /* Deck efgpri */
      SUBROUTINE EFGPRI(VMAT,TEXT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <qm3.h>
#include <maxaqn.h>
#include <maxorb.h>
      DIMENSION VMAT(3,3,MXCENT)
      CHARACTER*33 TEXT
#include <nuclei.h>
#include <symmet.h>
#include <ibtfun.h>
C
      NATOM = 1
      DO 10 I = 1, NUCIND
chj-s-090516 : defined ISUBSI, NSISY(0) for .NOT.QM3 so test is passed
        IF( (ISUBSY(I) .EQ. 0) .AND. (ISUBSI(I) .LE. NSISY(0)) ) THEN
          DO 20 J = 0, MAXOPR
            IF (IBTAND(ISTBNU(I),J) .EQ. 0) THEN
              WRITE (LUPRI,'(/A,3X,A)') TEXT,NAMDEP(NATOM)
              CALL OUTPUT(VMAT(1,1,NATOM),1,3,1,3,3,3,1,LUPRI)
              NATOM = NATOM + 1
            END IF
 20       CONTINUE
        END IF
 10   CONTINUE
      RETURN
      END
C  /* Deck iniall */
      SUBROUTINE ABAINIALL
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
C
      CALL NUCINI
      CALL ONEINI
      CALL TWOINI
      CALL RHSINI
      CALL ORTINI
      CALL RESINI
      CALL TROINI
      CALL FMMINI
      CALL DIPINI
      CALL QPGINI
      CALL RLXINI
      CALL ANAINI
      CALL VIBINI
      CALL FLTINI
      CALL EXCINI
      CALL LRSINI
      CALL MAGINI
      CALL SPNINI
      CALL AATINI
      CALL TRPINI
      CALL LNRINI
      CALL AQRINI
      RETURN
      END
