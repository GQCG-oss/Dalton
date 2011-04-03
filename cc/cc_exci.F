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
c*DECK CC_EXCI
       SUBROUTINE CC_EXCI(WORK,LWORK,LIST,APROXR12)
C
C----------------------------------------------------------------------
C
C     Purpose: Direct calculation of Coupled Cluster
C              excitation energies.
C
C        Singlet:
C
C              CCS, CC2, CCSD, CC3
C
C              CCSDT-1a, CCSDT-1b
C
C              CC(2), CCSDR(T), CCSDR(3), CCSDR(1A), CCSDR(1B)
C
C              CCS   = TDA    = CIS
C              CC(2) = CIS(D)
C
C        Triplet:
C
C              CCS, CC2, CCSD
C
C     The first section solves iteratively for CC excitation energies.
C     Both for left and right excitation energies.
C     The next section calculates perturbational corrections CC(2) and
C     CCSDR().
C
C     Written by Ove Christiansen August/November 1994
C
C     Version 3, December 1996.
C
C
C-----------------------------------------------------------------------
C
#include <implicit.h>
#include <pgroup.h>
#include <dummy.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <codata.h>
#include <ccsections.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccexci.h> 
#include <ccexgr.h>
#include <ccfdgeo.h>
#include <ccgr.h>
#include <ccfro.h>
#include <ccinftap.h>
C
      LOGICAL LOCDBG
      PARAMETER(LOCDBG=.FALSE.)
C
      LOGICAL KAJ,LINQCC,RSPIM2,LTRIP(4),LRST1,LRST2,LPROJECT,TRIPLET
      LOGICAL LREDS
      DIMENSION WORK(LWORK),NCCEXSAV(8,3),CONT(28)
      CHARACTER MODEL*24,MODELP*24,MODEL1*24,CHSYM*4,FILEX*10,FILOLD*10
      CHARACTER MOPRPC*10
      CHARACTER CDIP*1,CDUM*3, SPACE*32, MODELSCR*24
      CHARACTER*8 FC1AM,FC2AM,FRHO1,FRHO2,FRHO12,FC12AM,FS12AM,FS2AM
      CHARACTER*(3) LIST
      CHARACTER*(3) APROXR12
      PARAMETER (FC1AM ='CCR_C1AM',FC2AM ='CCR_C2AM')
      PARAMETER (FRHO1 ='CCR_RHO1',FRHO2 ='CCR_RHO2')
      PARAMETER (FRHO12='CCR_RH12',FC12AM='CCR_C12M',FS12AM='CCR_S12M')
      PARAMETER (FS2AM ='CCR_S2DM')
      PARAMETER (TWO = 2.0D00,TF=1.0D-04, TLOVLP = 0.5D0 )
      PARAMETER (XMONE = -1.0D00 , THRLDP = 1.0D-20)
      PARAMETER (THRDEGEN = 1.0D-8)
      CHARACTER*8 LABEL,LABEL1
C
#include <leinf.h>
C
      CALL QENTER('CC_EXCI')
C
      LUFC1  = -1
      LUFC2  = -1
      LUFC12 = -1
      LUFR1  = -1
      LUFR2  = -1
      LUFR12 = -1
      LUFS12 = -1
      LUFS2  = -1
C
      IF (CCP2)   CCS = .TRUE.
      TRIPLET = .FALSE.
      LABEL  = "Excita  "
      LABEL1 = "Exctot  "
C
C----------------------------------------------------------
C     Reinit NEXCI: Necessary if symmetry has been reduced.
C----------------------------------------------------------
C
      NEXCI  = 0
      NTRIP  = 0
      DO ISYM = 1,NSYM
         ISYOFE(ISYM) = NEXCI
         ITROFE(ISYM) = ISYOFE(ISYM) + NCCEXCI(ISYM,1)
         NEXCI        = ITROFE(ISYM) + NCCEXCI(ISYM,3)
         NTRIP        = NTRIP        + NCCEXCI(ISYM,3)
         DO IEX = ISYOFE(ISYM)+1, NEXCI
            ISYEXC(IEX) = ISYM
         END DO
         DO IEX = ISYOFE(ISYM)+1, ITROFE(ISYM)
            IMULTE(IEX) = 1
         END DO
         DO IEX = ITROFE(ISYM)+1, NEXCI
            IMULTE(IEX) = 3
         END DO
      END DO
      IF (IPRINT.GT.15) THEN
         WRITE(LUPRI,*) 'IN CC_EXCI after Reinit'
         WRITE(LUPRI,*) 'NEXCI: ',NEXCI
         WRITE(LUPRI,*) 'Singlet: ',(NCCEXCI(J,1),J=1,NSYM)
         WRITE(LUPRI,*) 'Triplet: ',(NCCEXCI(J,3),J=1,NSYM)
         WRITE(LUPRI,*) 'ISYOFE:',(ISYOFE(J), J=1,NSYM)
         WRITE(LUPRI,*) 'ITROFE:',(ISYOFE(J), J=1,NSYM)
         WRITE(LUPRI,*) 'ISYEXC:',(ISYEXC(J), J=1,NEXCI)
         WRITE(LUPRI,*) 'IMULTE:',(IMULTE(J), J=1,NEXCI)
         WRITE(LUPRI,*) 'EIGVAL:',(EIGVAL(J), J=1,NEXCI)
      ENDIF
C
C-----------------------
C     Initialize Leinfi.
C-----------------------
C
      CALL CCLR_LEINFI(TRIPLET)
      NTAMP  = NT1AMX + NT2AMX
      IF (CCS.OR.CCP2) NTAMP = NT1AMX
      LINQCC = .FALSE.
      NLOAD  = 1
C
C-------------------------------------------------------------------
C     Test of Finite difference jacobian and linear transformations.
C     or test calculation by finite difference jacobian.
C-------------------------------------------------------------------
C
      IF (LIST(1:2) .EQ. 'RE') ISIDE  = +1
      IF (LIST(1:2) .EQ. 'LE') ISIDE  = -1

      IF (JACTST.OR.FDJAC) THEN
         CALL CCLR_JAC(ISIDE,WORK,LWORK,APROXR12)
         IF (JACTST) THEN
            CALL AROUND( ' END OF JACTST ' )
            CALL QEXIT('CC_EXCI')
            RETURN
         ENDIF
      ENDIF
      IF ((.NOT.FDJAC).AND.(.NOT.JACTST).AND.FDEXCI) THEN
         RSPIM2 = RSPIM
         RSPIM  = .FALSE.
         FDJAC  = .TRUE.
         CALL CCLR_FDJAC(NT1AMX,NT2AMX,WORK,LWORK,APROXR12)
         FDJAC  = .FALSE.
         RSPIM  = RSPIM2
      ENDIF
C
C-----------------------------------------------------------------------
C     Calculation of excit. energies from finite difference Jacobian:
C-----------------------------------------------------------------------
C
      IF (FDEXCI) THEN
         CALL AROUND('CC_EXCI: START OF CALCULATION OF FD EXCI.')
         CALL CCLR_FDEXCI(WORK,LWORK)
         CALL AROUND('CC_EXCI: END OF CALCULATION OF FD EXCI. ')
         IF (.NOT.CCEXCI) THEN
            CALL AROUND( ' ONLY CALCULATION FROM FD JACOBIAN' )
            CALL QEXIT('CC_EXCI')
            RETURN
         ENDIF
      ENDIF
C
C---------------------------------------------
C     Header of Excitation Energy calculation.
C---------------------------------------------
C
      CALL TIMER('START ',TIMEIN,TIMOUT)
      WRITE (LUPRI,'(1X,A,/)') '  '
      WRITE (LUPRI,'(1X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1X,A)')
     *'*<<<<<<<<<< OUTPUT FROM COUPLED CLUSTER LINEAR RESPONSE >'//
     *'>>>>>>>>>*'
      IF ( CCEXCI ) THEN
         WRITE (LUPRI,'(1X,A)')
     *   '*                                                        '//
     *   '         *'
         WRITE (LUPRI,'(1X,A)')
     *   '*<<<<<<<<<<     CALCULATION OF EXCITATION ENERGIES      >'//
     *   '>>>>>>>>>*'
      ENDIF
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1X,A,/)')
     *'*********************************************************'//
     *'**********'
C
      MODEL = 'CCSD'
      MOPRPC = 'CCSD      '
      IF (CC2) THEN
         CALL AROUND( ' CC2 Excitation Energies ')
         MODEL = 'CC2'
         MOPRPC = 'CC2       '
      ENDIF
      IF (CCP2) THEN
         CALL AROUND( ' CCS and CC(2) Excitation Energies ')
         MODEL = 'CC(2)'
         MOPRPC ='CC(2)     '
      ENDIF
      IF (CCS.AND.(.NOT.CIS)) THEN
         CALL AROUND( ' CCS Excitation Energies ')
         MODEL = 'CCS'
         MOPRPC ='CCS       '
      ENDIF
      IF (CCS.AND.CIS) THEN
         CALL AROUND( ' CIS Excitation Energies ')
         MODEL = 'CCS'
         MOPRPC ='CCS       '
      ENDIF
      IF (CC3  ) THEN
         CALL AROUND( ' CC3 Excitation Energies ')
         MODEL = 'CC3'
         MOPRPC ='CC3       '
      ENDIF
      IF (CC1B .AND. ( .NOT. CC1a) ) THEN
         CALL AROUND( ' CCSDT-1b Excitation Energies ')
         MODEL = 'CCSDT-1b'
         MOPRPC ='CCSDT-1b  '
      ENDIF
      IF (CC1A ) THEN
         CALL AROUND( ' CCSDT-1a Excitation Energies ')
         MODEL = 'CCSDT-1a'
         MOPRPC ='CCSDT-1a  '
      ENDIF
      IF (CCR3) THEN
         CALL AROUND( ' CCSD and CCSDR(3) Excitation Energies')
         MODEL = 'CCSD'
         MOPRPC ='CCSD      '
      ENDIF
      IF (CCRT) THEN
         CALL AROUND( ' CCSD and CCSDR(T) Excitation Energies')
         MODEL = 'CCSD'
         MOPRPC ='CCSD      '
      ENDIF
      IF (CCR1A) THEN
         CALL AROUND( ' CCSD and CCSDR(1a) Excitation Energies')
         MODEL = 'CCSD'
         MOPRPC ='CCSD      '
      ENDIF
      IF (CCR1B) THEN
         CALL AROUND( ' CCSD and CCSDR(1b) Excitation Energies')
         MODEL = 'CCSD'
         MOPRPC ='CCSD      '
      ENDIF
      IF (CCSD) THEN
         CALL AROUND( ' CCSD Excitation Energies ')
         MODEL = 'CCSD'
         MOPRPC ='CCSD      '
      ENDIF
      IF (CIS) THEN
        MODELP = 'CIS'
        MOPRPC ='CIS       '
      ELSE
        MODELP = MODEL
        IF (CCR12 .AND. (.NOT.CCS)) THEN
          CALL CCSD_MODEL(MODELP,LENMOD,24,MODEL,24,APROXR12)
        END IF
      ENDIF
C
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CC_EXCI-1: Workspace:',LWORK
C
C     ********************************************
C     ***** Determine CC excitation energies *****
C     ********************************************
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      KEXCI  = 1
      KEXCP  = KEXCI  + NEXCI
      KEXCP2 = KEXCP  + NEXCI
      KEXCP3 = KEXCP2 + NEXCI
      KEXCP4 = KEXCP3 + NEXCI
      KT1P   = KEXCP4 + NEXCI
      KWRK0  = KT1P   + NEXCI
      LWRK0  = LWORK  - KWRK0
      IF (LWRK0.LT. 0 ) CALL QUIT(' TOO LITTLE WORKSPACE IN CC_EXCI')
      IF (IPRINT.GT.10) WRITE(LUPRI,*) 'CC_EXCI-1b: Workspace:',LWRK0
C
      KEXCPS  = KEXCP
      KEXCPS2 = KEXCP2
      KEXCPS3 = KEXCP3
      KEXCPS4 = KEXCP4
      KEXCIS  = KEXCI
      KT1PS   = KT1P
C
C---------------------------------
C     test with explicit jacobian.
C---------------------------------
C
      IF (JACEXP) THEN
         ISYMTR = 1    
         CALL CCLR_LEINFI(TRIPLET)
         CALL CC_JACEXP(WORK(KWRK0),LWRK0)
      ENDIF
C
C-------------------------------------
C     Start loop over symmetriclasses.
C-------------------------------------
C
      DO 9000 ISYM = 1, NSYM
C
C
C--------------------------------------
C       Start loop over multiplicities: 
C--------------------------------------
C
        DO 8000 IMULT = 1, 3, 2
C
C-----------------------------
C        Initialize variables.
C-----------------------------
C
         IF      (IMULT.EQ.1) THEN
            TRIPLET = .FALSE.
         ELSE IF (IMULT.EQ.3) THEN
            TRIPLET = .TRUE.
         ELSE
            CALL QUIT('Illegal multiplicity in CC_EXCI.')
         END IF
C
         IF (TRIPLET .AND. NCCEXCI(ISYM,IMULT).GT.0) THEN
           IF (.NOT. (CCS.OR.CC2.OR.CCSD.OR.CC3)) THEN
             WRITE (LUPRI,*) 'WARNING: ',MODELP,
     &             ' TRIPLET EXCITATION ENERGIES',
     &               ' ARE NOT (YET) AVAILABLE.'
           END IF
         END IF
C
         LREDS  = CCR12 .OR. (CCSDT .AND. CCSDT_DIIS) 
         LINQCC = .FALSE.
C
         IF (LIST(1:2) .EQ. 'RE') ISIDE  = +1
         IF (LIST(1:2) .EQ. 'LE') ISIDE  = -1
         IF ((LIST(1:2) .EQ. 'LE').AND.CCSDT) THEN
            IF (.NOT. CC3) THEN
               CALL QUIT('Left excitation energies and '//
     *                   '(non cc3) iterative triples not allowed')
            ENDIF
         ENDIF
C
         ISYMTR = ISYM
         IF (CCS) THEN
           NCCVAR = NT1AM(ISYMTR)
         ELSE
           NCCVAR = NT1AM(ISYMTR) + NT2AM(ISYMTR)
           IF (CCR12)   NCCVAR = NCCVAR + NTR12AM(ISYMTR)
           IF (TRIPLET) NCCVAR = NCCVAR + NT2AMA(ISYMTR)
         END IF
         CALL CCLR_LEINFI(TRIPLET)
C
C---------------
C        OUTPUT.
C---------------
C
         IF (ISYM .GT. 1) WRITE(LUPRI,'(//A)')
     *     ' ***********************************************'
     *    //'********************************'
         WRITE (LUPRI,'(/A)')   ' --------------------------'
         WRITE (LUPRI,'(A,I5)') ' Symmetry class Nr.: ',ISYM
         WRITE (LUPRI,'(A,I5)') ' Multiplicity      : ',IMULT
         WRITE (LUPRI,'(A,/)')  ' --------------------------'
         WRITE (LUPRI,'(A,I10)')
     *   ' Length of Excitation vectors in this class is:',NCCVAR
C
         IF (NCCVAR .LT.NCCEXCI(ISYM,IMULT)) THEN
            WRITE(LUPRI,*) 'Demand for more excitation energies than '
     *           //'amplitudes with this symmetry/multiplicity!!!!'
            WRITE(LUPRI,*) 'Nr. of amplitudes with symmetry ',ISYM,
     *                 ' and multiplicity ',IMULT,' is ',NCCVAR
            WRITE (LUPRI,*) 'Nr. of exci.e. requested with sym ',ISYM,
     *                          ' IS ',NCCEXCI(ISYM,IMULT)
            NCCEXCI(ISYM,IMULT) = NCCVAR
         ENDIF
         NCCEXSAV(ISYM,IMULT) = NCCEXCI(ISYM,IMULT)
C        NCC(ISYM,IMULT) = NCCEXCI(ISYM,IMULT)
C
         IF (NCCEXCI(ISYMTR,IMULT).EQ.0) THEN
            WRITE (LUPRI,'(A)')
     *   ' No Excitation Energies calculated in this symmetry class'
            GOTO 8000
         ENDIF
C
         IF (IMULT.EQ.1) THEN
           IBOFF = ISYOFE(ISYM)
         ELSE
           IBOFF = ITROFE(ISYM)
         END IF
         NSIM    = NCCEXCI(ISYM,IMULT)
         IRST    = IBOFF + 1
         NSIMR   = NSIM
         IREND   = IRST  + NSIMR - 1
C
C---------------------------------------------------------------
C        Allocation of work space.
C        NB: The gradient vectors is not in memory at this time.
C---------------------------------------------------------------
C
C
         KIPLAC = KWRK0 
         KREDH  = KIPLAC + MAXRED
         KREDGD = KREDH  + MAXRED*MAXRED
         KEIVAL = KREDGD + MAXRED*NSIMR
         KSOLEQ = KEIVAL + MAXRED
         KWRK1  = KSOLEQ + MAXRED*MAXRED

         KREDS  = KWRK1 
         IF (LREDS) KWRK1  = KREDS  + MAXRED*MAXRED

         IF (CCSDT .AND. CCSDT_DIIS) THEN
            KAMAT  = KWRK1
            KITRAN = KAMAT  + (MXDIIS+1)*(MXDIIS+1)*NSIM
            KCONV  = KITRAN + NSIM
            KRNORM = KCONV  + NSIM
            KWRK1  = KRNORM + NSIM
         END IF

         LWRK1  = LWORK  - KWRK1
C
         IF (LWRK1.LT. 0 )
     *      CALL QUIT(' TOO LITTLE WORKSPACE IN CC_EXCI')
         CALL DZERO(WORK(KEIVAL),MAXRED)
C
C--------------------------------------------------------------
C        If CCR3 do spaghetti goto to pert. corrections section.
C--------------------------------------------------------------
C
         IF ( CCR3 ) GOTO 1200
C
         WRITE (LUPRI,'(A,I3,A)')
     *   ' Converging for',NCCEXCI(ISYM,IMULT),' roots.'
C
C-----------------------------------------------------------------------
C        Choose start omega for CC3: EOMINP contains input choice, 
C        CCSD or CCSDR excitations energies.
C        Choose the first, for omesc and mxtomn this means the highest,
C        and thus probably the hightest shift.
C-----------------------------------------------------------------------
C
         IF (CCSDT) THEN
            ECURR = EOMINP(1,ISYMTR,IMULT)
         ENDIF
C
         IF ( STSD.AND.CCSDT ) THEN
            CCSDT = .FALSE.
            KAJ   = .TRUE.
         ENDIF
C
C-------------------
C        Open files.
C-------------------
C
         CALL CC_FILOP(FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                 FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 FS12AM,LUFS12,FS2AM,LUFS2)
C
C-----------------------------
C        Create start vectors.
C-----------------------------
C
         LPROJECT = .FALSE.
         CALL CCEQ_STR(FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 FS12AM,LUFS12,FS2AM,LUFS2,LPROJECT,ISTATPRJ,
     *                 TRIPLET,ISIDE,IRST,NSIMR,NUPVEC,
     *                 NREDH,WORK(KEIVAL),WORK(KIPLAC),
     *                 WORK(KWRK1),LWRK1,LIST)
         CALL FLSHFO(LUPRI)
C
C------------------------------------------
C        Solve equations by call to solver.
C------------------------------------------
C
         IF (OMESC .AND. CCSDT .AND. CCSDT_DIIS) THEN

           DO ICOUNT = 1, NSIM
             IF (TRIPLET) THEN
               WORK(KEIVAL-1+ICOUNT) = EIGVAL(ITROFE(ISYM)+ICOUNT)
             ELSE
               WORK(KEIVAL-1+ICOUNT) = EIGVAL(ISYOFE(ISYM)+ICOUNT)
             END IF
           END DO

           CALL CCDIIS_SOL(FRHO1,LUFR1,FRHO2,LUFR2,
     *                     FC1AM,LUFC1,FC2AM,LUFC2,LIST,IRST,
     *                     TRIPLET,ISIDE,NSIMR,NUPVEC,
     *                     NREDH,WORK(KREDH),WORK(KREDS),WORK(KEIVAL),
     *                     WORK(KSOLEQ),WORK(KAMAT),
     *                     WORK(KITRAN),WORK(KCONV),WORK(KRNORM),
     *                     WORK(KWRK1),LWRK1)
         ELSE
           CALL CCEQ_SOL(LIST,LPROJECT,ISTATPRJ,ECURR,
     *                   FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                   FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                   LREDS,WORK(KREDS),FS12AM,LUFS12,FS2AM,LUFS2, 
     *                   LINQCC,TRIPLET,ISIDE,IRST,NSIMR,NUPVEC,
     *                   NREDH,WORK(KREDH),WORK(KEIVAL),
     *                   WORK(KSOLEQ),WORK(KWRK1),LWRK1,APROXR12)
         END IF

         ! test eigenvector: in case of degeneracies orthogonalize
         ! eigenvectors to the same eigenvalue
         CALL CCTSTSOL(WORK(KSOLEQ),WORK(KEIVAL),THRDEGEN,NREDH,NSIMR)
 
         IF ( STSD.AND.KAJ ) THEN
            CCSDT = .TRUE.
         ENDIF
C
C-------------------------------------
C       Write out Excitation energies.
C-------------------------------------
C
         WRITE (LUPRI,'(//A,I3)') ' SYMMETRY CLASS NR.  ',ISYM
         WRITE (LUPRI,'(A,I3)')   ' MULTIPLICITY        ',IMULT
         IF (.NOT. CCSDT) THEN
            IF (LIST(1:2) .EQ. 'RE') THEN
               WRITE (LUPRI,'(//1X,A10,1X,A26)')
     *                MODELP,'right excitation energies:'
            ELSE IF (LIST(1:2) .EQ. 'LE') THEN
               WRITE (LUPRI,'(//1X,A10,1X,A25)')
     *                MODELP,'left excitation energies:'
            ENDIF
            WRITE (LUPRI,'(A)')
     *      ' ===================================='
         ELSE
            IF ( STSD ) THEN
               WRITE (LUPRI,'(//A11,A10,A12)')
     *         ' CCSD with ',MODELP,' amplitudes:'
            WRITE (LUPRI,'(A)')
     *      ' ==============================='
            ELSE
               WRITE (LUPRI,'(//,1X,A10,A33,F10.6)')
     *         MODELP,' excitation energies with Omega= ',ECURR
            WRITE (LUPRI,'(A)')
     *      ' ====================================================='
            ENDIF
C
         ENDIF
C
         WRITE (LUPRI,'(A//A/A)')
     *   ' (conversion factor used: 1 au = 27.2113957 eV)',
     *   ' Excitation no.       Hartree               eV ',
     *   ' --------------       -------               --'
         DO 400 I = 1,NCCEXCI(ISYM,IMULT)
            WRITE (LUPRI,'(I10,4X,2F20.10)')
     *         I,WORK(KEIVAL-1+I),WORK(KEIVAL-1+I)*XTEV
  400    CONTINUE
C
         WRITE (LUPRI,'(//A,2I5,/A/A)')
     *   ' Total excited state energies for states of symmetry/spin',
     *    ISYM,IMULT,
     *   ' Excitation no.       Energy (Hartree) ',
     *   ' ------------------------------------- '
         DO 401 I = 1,NCCEXCI(ISYM,IMULT)
            WRITE (LUPRI,'(A2,2I5,F30.15)')
     *         '@@',ISYM,I,WORK(KEIVAL-1+I)+ECCGRS
  401    CONTINUE
C
C------------------------------------
C       Analysis of solution vectors.
C------------------------------------
C
         NVARPT = NCCVAR+ 2*NALLAI(ISYM)
         KWRK2  = KWRK1 + NVARPT
         LWRK2  = LWORK - KWRK2
         NSIMUL = MIN(NSIMR,LWRK2/NCCVAR)
         NBATCH = (NSIMR -1)/NSIMUL + 1
         IOFF1  = 1
         ICOUNT = 0
         DO 500 I = 1,NBATCH
            IOFF2 = MIN(NSIMUL,NCCEXCI(ISYM,IMULT) - (I-1)*NSIMUL)
            CALL CCCONV(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,
     *                  TRIPLET,CCR12,NREDH,IOFF1,IOFF2,WORK(KSOLEQ),
     *                  WORK(KWRK2),WORK(KWRK1))

            IF (LOCDBG) THEN
               WRITE(LUPRI,*) 'Overlap matrix for solution vectors:'
               CALL PROVLP(WORK(KWRK2),WORK(KWRK2),NCCVAR,IOFF2,
     *                     WORK(KWRK2+IOFF2*NCCVAR),LUPRI)
            END IF
            IF ( IPRINT .GT. 30 ) THEN
               CALL AROUND('CC_EXCI: VECTORS IN MO BASIS' )
               CALL OUTPUT(WORK(KWRK2),1,NCCVAR,1,NCCEXCI(ISYM,IMULT),
     *                     NCCVAR,NCCEXCI(ISYM,IMULT),1,LUPRI)
            ENDIF
C
            DO 510 J = 1,IOFF2
               ICOUNT = ICOUNT + 1
               IF (TRIPLET) THEN
                 ILSTNR = ITROFE(ISYM) + ICOUNT
               ELSE
                 ILSTNR = ISYOFE(ISYM) + ICOUNT
               END IF

               WRITE(LUPRI,'(//1X,2A,I3)') 'Analysis of the Coupled ',
     *          'Cluster Excitation Vector Number :', ICOUNT
               WRITE(LUPRI,'(1X,2A)') '-----------------------------',
     *                             '--------------------------------'
               WRITE(LUPRI,'(/10X,A,23X,F10.4,A)')
     *          'Excitation Energy :  ',WORK(KEIVAL+ICOUNT-1)*XTEV,' eV'
               WORK(KEXCIS+(ICOUNT-1)) = WORK(KEIVAL+ICOUNT-1)
               IF (TRIPLET) THEN
                 KCAM1 = KWRK2 + (J-1)*NCCVAR
                 KCAMP = KCAM1 + NT1AM(ISYMTR)
                 KCAMM = KCAMP + NT2AM(ISYMTR)
                 CALL CC_PRAM3(WORK(KCAM1), WORK(KCAMP), WORK(KCAMM),
     *                         WORK(KT1PS+(ICOUNT-1)),PTP,PTM,ISYMTR,
     *                         .FALSE.)
               ELSE
                 CALL CC_PRAM(WORK(KWRK2+(J-1)*NCCVAR),
     *                        WORK(KT1PS+(ICOUNT-1)),ISYMTR,
     *                        .FALSE.)
               END IF
               CALL FLSHFO(LUPRI)
 
C              -----------------------------------------------------
C              Save solution vectors on file and excitation energies
C              in the result array EIGVAL:
C              (if we compute both left and right eigenvector, 
C               don't overwrite with left eigenvalues the right ones)
C              -----------------------------------------------------
               IOPT = 3
               CALL CC_WRRSP(LIST,ILSTNR,ISYM,IOPT,MODELP,DUMMY,
     *                       WORK(KWRK2+(J-1)*NCCVAR),
     *                       WORK(KWRK2+(J-1)*NCCVAR+NT1AM(ISYM)),
     *                       WORK(KWRK1),NVARPT)
               IF (CCR12) THEN
                 IOPT = 32
                 CALL CC_WRRSP(LIST,ILSTNR,ISYM,IOPT,MODELP,DUMMY,DUMMY,
     *                WORK(KWRK2+(J-1)*NCCVAR+NT1AM(ISYM)+NT2AM(ISYM)),
     *                WORK(KWRK1),NVARPT)
               END IF
 
               IF (LIST(1:2).EQ.'RE'.OR.LHTR) THEN
                 EIGVAL(ILSTNR) = WORK(KEIVAL+ICOUNT-1)
               ENDIF

  510       CONTINUE
            IOFF1 = IOFF1 + NSIMUL
  500    CONTINUE
C
         IF (EXCI_CONT .AND. LIST(1:2).EQ.'RE') THEN
           DO I = 1, NSIMR
              KC1AM = KWRK0
              KWRK1 = KC1AM + NT1AM(ISYM)
              IF (.NOT.CCS) THEN
                KC2AM = KWRK1
                KWRK1 = KC2AM + NT2AM(ISYM)
              END IF
              IF (CCR12) THEN
                KC12AM = KWRK1
                KWRK1  = KC12AM + NTR12AM(ISYM)
              END IF
              LWRK1 = LWORK - KWRK1
              IF(LWRK1.LT.0) CALL QUIT('OUT OF WORKSPACE IN CC_EXCI')

              IF (TRIPLET) THEN
                ILSTNR = ITROFE(ISYM) + I
              ELSE
                ILSTNR = ISYOFE(ISYM) + I
              END IF

              IOPT = 3
              CALL CC_RDRSP(LIST,ILSTNR,ISYM,IOPT,MODEL,
     *                      WORK(KC1AM),WORK(KC2AM))
              IF (CCR12) THEN
                IOPT = 32
                CALL CC_RDRSP(LIST,ILSTNR,ISYM,IOPT,MODEL,
     *                        DUMMY,WORK(KC12AM))
              END IF

              CALL CC_ATRR(ECURR,ISYM,1,WORK(KWRK0),LWRK0,.TRUE.,CONT,
     &                     APROXR12)
              CALL CC_EXCIT_CONT(I,ISYM,TRIPLET,MODELP,
     *                           CONT(4),CONT(1))
           END DO
         END IF
C
         CALL FLSHFO(LUPRI)
C
C----------------------------------------------------------
C==========================================================
C       IF CCSDT we may solve until self consistency:
C==========================================================
C----------------------------------------------------------
C
         IF (OMESC.AND.CCSDT.AND.(.NOT.CCSDT_DIIS)) THEN
C
         CALL DZERO(WORK(KEXCIS),NCCEXCI(ISYMTR,IMULT))
         CALL DZERO(WORK(KT1PS),NCCEXCI(ISYMTR,IMULT))
C
C-----------------------------------------------------------
C        Solve for all NOMINP eigenvalues self-consistently.
C-----------------------------------------------------------
C
         DO 1000 IOM = 1, NOMINP(ISYMTR,IMULT)
C
C-----------------------
C        Solve equtions.
C-----------------------
C
            ISTATE  = IOMINP(IOM,ISYMTR,IMULT)
C
            NCCEXCI(ISYM,IMULT) = ISTATE
C           NCC(ISYM,IMULT) = ISTATE
            NSIMR = ISTATE
            NSIDE = 1
C
C--------------------------------
C           Create start vectors.
C--------------------------------
C
            LPROJECT = .FALSE.
            CALL CCEQ_STR(FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                    FS12AM,LUFS12,FS2AM,LUFS2,LPROJECT,ISTATPRJ,
     *                    TRIPLET,ISIDE,IRST,NSIMR,NUPVEC,
     *                    NREDH,WORK(KEIVAL),WORK(KIPLAC),
     *                    WORK(KWRK1),LWRK1,LIST)
            CALL FLSHFO(LUPRI)
C
C-------------------------------------------------------------------
C           If EOMINP read in is different from zero, use this
C           omega.
C-------------------------------------------------------------------
C
            ITSC = 1
            IF ( STSD ) ITSC = 0
C
            IF (TRIPLET) THEN
               ILSTNR = ITROFE(ISYM) + ISTATE
            ELSE
               ILSTNR = ISYOFE(ISYM) + ISTATE
            ENDIF
C            ECURR1 = EIGVAL(ISTATE) 
            ECURR1 = EIGVAL(ILSTNR) 
            IF (ITSC .EQ. 0 ) THEN
               EITOL = 1.0D-3
               IF (ABS(EOMINP(IOM,ISYMTR,IMULT)).GT.EITOL) THEN
                  ECURR1 = EOMINP(IOM,ISYMTR,IMULT)
                  WRITE(LUPRI,'(/,1X,A24,F10.6)')
     *                 'Omega put equal to input',ECURR1
               ENDIF
            ENDIF
C
C-------------------------------------------------
C
C           LOOP OMESC
C
C           Loop until solution is selfconsistent.
C
C-------------------------------------------------
C
 2000       CONTINUE
C
            ITSC   = ITSC + 1
            ECURR  = ECURR1
C
C
C---------------------------------------------
C           Solve equations by call to solver.
C---------------------------------------------
C
            CALL CCEQ_SOL(LIST,LPROJECT,ISTATPRJ,ECURR,
     *                    FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                    FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                    LREDS,WORK(KREDS),FS12AM,LUFS12,FS2AM,LUFS2,
     *                    LINQCC,TRIPLET,ISIDE,IRST,NSIMR,NUPVEC,
     *                    NREDH,WORK(KREDH),WORK(KEIVAL),
     *                    WORK(KSOLEQ),WORK(KWRK1),LWRK1,APROXR12)
C
C-------------------------------------
C       Write out Excitation energies.
C-------------------------------------
C
            WRITE (LUPRI,'(//,1X,A10,A33,F10.6)')
     *         MODELP,' excitation energies with Omega= ',ECURR
            WRITE (LUPRI,'(A/A//A/A/)')
     *      ' =====================================================',
     *      ' (conversion factor used: 1 au = 27.2113957 eV)',
     *      ' Excitation no.       Hartree               eV',
     *      ' --------------       -------               --'
C
            DO 900 I = 1,NCCEXCI(ISYM,IMULT)
C
               WRITE (LUPRI,'(I10,2F20.6)')
     *            I,WORK(KEIVAL-1+I),WORK(KEIVAL-1+I)*XTEV
C
  900       CONTINUE
C
            ECURR1 = WORK(KEIVAL-1+ISTATE)
C
            IF (IPRINT.GT. 1) THEN
               WRITE(LUPRI,'(/1x,A13,I2,A15,2X,F10.6)')
     *              'Omega in the ',ITSC,'-th iteration: ',ECURR1
               WRITE(LUPRI,'(1x,A32,F10.6)')
     *              'Omega in previous iteration:    ',ECURR
               WRITE(LUPRI,'(1x,A32,F10.6,/)')
     *              'Tolerance for self consistency: ',TOLSC
            ENDIF
C
            IF (ABS(ECURR1-ECURR).LT. TOLSC) THEN
C
               WRITE(LUPRI,'(/,1x,A23,F10.6,A11,I3,A4,I3,A11/)')
     *                'Converged root to diff.'
     *                ,ECURR-ECURR1,' for root ',ISTATE
     *                ,' in ', ITSC ,' iterations'
               EOMINP(IOM,ISYMTR,IMULT) = ECURR1
C
               WRITE (LUPRI,'(//A,2I3,/A/A)')
     *   ' Total excited state energies for states of symmetry/spin ',
     *     ISYM,IMULT,
     *   ' Excitation no.       Energy (Hartree) ',
     *   ' ------------------------------------- '
                WRITE (LUPRI,'(A3,2I5,F30.15)')
C
CKH   *         '@@@',ISYM,I,WORK(KEIVAL-1+I)+ECCGRS
C
     *         '@@@',ISYM,I-1,EOMINP(IOM,ISYMTR,IMULT)+ECCGRS
C
C
            ENDIF
C
C----------------------------------------------------------------------
C           If equations not are selfconsistent with respect to omega,
C           then take solution vectors and energy and start again.
C           Else write out the analysis of the vectors.
C----------------------------------------------------------------------
C
            NVARPT = NCCVAR+ 2*NALLAI(ISYM)
            KWRK2  = KWRK1 + NVARPT
            LWRK2  = LWORK - KWRK2

            THRESH = 0.05
            MAXLIN = 100
            NSIMUL = MIN(NCCEXCI(ISYM,IMULT),LWRK2/NCCVAR -1 )
            NBATCH = (NCCEXCI(ISYM,IMULT)-1)/NSIMUL + 1
            IOFF1  = 1
            ICOUNT = 0
C
            NSIMA = 1
            NSIMB = 0
            NBLE3 = 0
C
            DO 600 I = 1,NBATCH
               IOFF2 = MIN(NSIMUL,NCCEXCI(ISYM,IMULT) - (I-1)*NSIMUL)
               CALL CCCONV(LUFC1,FC1AM,LUFC2,FC2AM,LUFC12,FC12AM,
     *                     TRIPLET,CCR12,NREDH,IOFF1,IOFF2,WORK(KSOLEQ),
     *                     WORK(KWRK2),WORK(KWRK1))
C
               IF ( IPRINT .GT. 30 ) THEN
                  CALL AROUND('CC_EXCI: VECTORS IN MO BASIS' )
                  CALL OUTPUT(WORK(KWRK2),1,NCCVAR,1,
     *                        NCCEXCI(ISYM,IMULT),
     *                        NCCVAR,NCCEXCI(ISYM,IMULT),1,LUPRI)
               ENDIF
C
C-----------------------------------------------------------------
C              If converged write out the analysis of the vectors.
C-----------------------------------------------------------------
C
               DO 610 J = 1,IOFF2
                  ICOUNT = ICOUNT + 1
C                 -------------------------------
C                 save vectors on standard files:
C                 -------------------------------
                  IF (TRIPLET) THEN
                    ILSTNR = ITROFE(ISYM) + ICOUNT
                  ELSE
                    ILSTNR = ISYOFE(ISYM) + ICOUNT
                  END IF
                  IOPT = 3
                  CALL CC_WRRSP(LIST,ILSTNR,ISYM,IOPT,MODEL,DUMMY,
     *                       WORK(KWRK2+(ICOUNT-1)*NCCVAR),
     *                       WORK(KWRK2+(ICOUNT-1)*NCCVAR+NT1AM(ISYM)),
     *                       WORK(KWRK1),NVARPT)
C
                  IF (ABS(ECURR1-ECURR).LE. TOLSC) THEN
C
                     WRITE(LUPRI,'(//1X,A,I2)')
     *'Analysis of the Coupled Cluster Excitation Vector Number : ',
     * ICOUNT
                  WRITE(LUPRI,'(1X,A)')
     *'-------------------------------------------------------------'
                  WRITE(LUPRI,'(/10X,A,23X,F10.4,A)')
     *   'Excitation Energy :  ',WORK(KEIVAL+ICOUNT-1)*XTEV,
     *   ' eV'
C
                   IF (TRIPLET) THEN
                      KCAM1 = KWRK2 + (ICOUNT-1)*NCCVAR
                      KCAMP = KCAM1 + NT1AM(ISYMTR)
                      KCAMM = KCAMP + NT2AM(ISYMTR)
                      CALL CC_PRAM3(WORK(KCAM1), WORK(KCAMP),
     *                              WORK(KCAMM),PT1,
     *                              PTP,PTM,ISYMTR,.FALSE.)
                    ELSE
                      CALL CC_PRAM(WORK(KWRK2 + (ICOUNT-1)*NCCVAR),
     *                             PT1,ISYMTR,.FALSE.)
                    ENDIF
                    IF (ICOUNT .EQ. ISTATE) THEN
                       WORK(KEXCIS+(ICOUNT-1)) = WORK(KEIVAL+ICOUNT-1)
                       WORK(KT1PS+(ICOUNT-1)) = PT1
                    ENDIF

C                   -------------------------------
C                   save / check excitation energy:
C                   -------------------------------
                    IF ((LIST(1:2) .EQ. 'LE').AND.(.NOT.LHTR)) THEN
                     IF(ABS(EIGVAL(ILSTNR)-WORK(KEIVAL+ICOUNT-1)).GT.TF)
     &                 WRITE(LUPRI,'(////,1X,A,//)') 
     &                   'Warning - Large difference between '//
     &                   'left and right excitation energies'
                    ELSE
                     EIGVAL(ILSTNR) = WORK(KEIVAL+ICOUNT-1)
                    END IF
C
C-----------------------------------------------------
C                   Restart if not selfconsistent yet.
C-----------------------------------------------------
C
                  ELSE
C
                    LPROJECT = .FALSE.
                    CALL CCEQ_STR(FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                            FS12AM,LUFS12,FS2AM,LUFS2,
     *                            LPROJECT,ISTATPRJ,
     *                            TRIPLET,ISIDE,IRST,NSIMR,NUPVEC,
     *                            NREDH,WORK(KEIVAL),WORK(KIPLAC),
     *                            WORK(KWRK1),LWRK1,LIST)
                    CALL FLSHFO(LUPRI)
C
                  ENDIF
C
  610          CONTINUE
C
               IOFF1 = IOFF1 + NSIMUL
  600       CONTINUE
C
C---------------------------------------------------------
C           If not self-consistent then return to 2000 and
C           solved eigenvalue equations with new omega.
C---------------------------------------------------------
C
            IF (ABS(ECURR1-ECURR).GT. TOLSC) THEN
C
               GOTO 2000
C
            ENDIF
C
 1000       CONTINUE
C
         ENDIF
C
C-------------------------------
C        Close and delete files.
C-------------------------------
C
         CALL CC_FILCL(FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *                 FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                 FS12AM,LUFS12,FS2AM,LUFS2)

C-----------------------------------------------------------
C        restore number of excitation within symmetry class:
C-----------------------------------------------------------

         NCCEXCI(ISYM,IMULT) = NCCEXSAV(ISYM,IMULT) 

C---------------------------------------------------------------------
C        First, check if left and right eigenvalues agree. Then
C        normalize Eigenvectors such that (RE|RE) = 1  and (LE|RE) = 1
C---------------------------------------------------------------------
        IF (LIST(1:2).EQ.'LE' .AND. (.NOT.LHTR)) THEN
         ! set threshold within which left and right eigenvalues
         ! should agree, else complain 
         THREIG = 10 * THREXC
         IF (OMESC.AND.CCSDT.AND.(.NOT.CCSDT_DIIS)) THREIG = TOLSC

         DO ICOUNT = 1, NCCEXCI(ISYM,IMULT)
           IF (TRIPLET) THEN
             ILSTNR = ITROFE(ISYM) + ICOUNT
           ELSE
             ILSTNR = ISYOFE(ISYM) + ICOUNT
           END IF
           IF(ABS(EIGVAL(ILSTNR)-WORK(KEXCIS+ICOUNT-1)).GT.THREIG)THEN
             WRITE(LUPRI,'(//,1X,A)') 'Warning - Large differ'
     *           //'ence between left and right excitation energies'
             WRITE(LUPRI,'(1X,A,I5)') 'Symmetry class:',ISYM
             WRITE(LUPRI,'(1X,A,I5)') 'Multiplicity  :',IMULT
             WRITE(LUPRI,'(1X,A,I5)') 'State number  :',ICOUNT
             WRITE(LUPRI,'(1X,A,2F20.10)') 'Eigenvalues   :',
     *         EIGVAL(ILSTNR),WORK(KEXCIS+ICOUNT-1)
             WRITE(LUPRI,'(1X,A,2F20.10)') 'Diff/Threshold:',
     *         EIGVAL(ILSTNR)-WORK(KEXCIS+ICOUNT-1),THREIG
             WRITE(LUPRI,'(1X,A,//)') 
     *         'Warning - Check your input and output carefully!!'
           END IF
         END DO

         CALL CCEXNORM(ISYM,IMULT,THREIG,WORK,LWORK)

        END IF
C
C====================================================================
C--------------------------------------------------------------------
C
C        PERTURBATIVE CORRECTIONS TO CCS OR CCSD EXCITATION ENERGIES.
C
C--------------------------------------------------------------------
C====================================================================
C
 1200 CONTINUE
C
      IF ((LIST(1:2).EQ.'LE').AND. (.NOT.TRIPLET) .AND. 
     *    (.NOT.CCR12) .AND.
     *    (CCP2.OR.((CCR3.OR.CCRT).OR.(CCR1A.OR.CCR1B)))) THEN
C
C---------------------------------------------------------------
C        If CCSDR(3) and first time (only calculating left CCSD)
C        then we skip the rest.
C---------------------------------------------------------------
C
         IF ( CCR3 .AND. ( .NOT. CCT)) GOTO 8000
C
         KE1    = KWRK1
         KE2    = KE1  + NCCEXCI(ISYM,IMULT)
         KE3    = KE2  + NCCEXCI(ISYM,IMULT)
         KE4    = KE3  + NCCEXCI(ISYM,IMULT)
         KWRK2  = KE4  + NCCEXCI(ISYM,IMULT)
         LWRK2  = LWORK - KWRK2
C
         CALL CC_PCEXCI(WORK(KEIVAL),WORK(KE1),WORK(KE2),WORK(KE3),
     *                  WORK(KE4),WORK(KWRK2),LWRK2)
C
         IF (CCP2 ) THEN
            WRITE (LUPRI,'(//,1X,A)')
     *      'Excitation energies with doubles corrections'
         ELSE IF (CCRT.OR.CCR3.OR.CCR1A.OR.CCR1B) THEN
            WRITE (LUPRI,'(//,1X,A)')
     *      'Excitation energies with triples corrections'
         ENDIF
C
         IF (CCP2) THEN
               WRITE (LUPRI,'(//A)')
     *   ' CC-model      Excitation no.         Hartree             eV'
         ENDIF
C
         IF (CCR3.OR.CCRT.OR.CCR1A.OR.CCR1B) THEN
              WRITE (LUPRI,'(//A)')
     *   ' CC-model               Excitation no.          Hartree'
     *         //'               eV'
         ENDIF
C
C------------------------------------------------
C        Loop over excitations in symmetry class.
C------------------------------------------------
C
         DO 1300 IST = 1,NCCEXCI(ISYM,IMULT)
C
            OME1 = WORK(KE1-1+IST)
            OME2 = WORK(KE2-1+IST)
            OME3 = WORK(KE3-1+IST)
            OME4 = WORK(KE4-1+IST)
C
C-----------------------------------------
C           Write out Excitation energies.
C-----------------------------------------
C
            IF (CCP2) THEN
               WRITE (LUPRI,'(A)')
     *   ' ----------------------------         --------            --'
     *         //'------'
               WRITE (LUPRI,'(A16,I10,2F20.6)') ' CCS            ',
     *         IST,WORK(KEIVAL-1+IST),WORK(KEIVAL-1+IST)*XTEV
               WRITE (LUPRI,'(A)')
     *   ' ----------------------------         --------            --'
     *         //'------'
C
               WRITE (LUPRI,'(A16,I10,2F20.6)') ' C(0)A(0+2)C(0) ',
     *         IST,OME1,OME1*XTEV
               WRITE (LUPRI,'(A16,I10,2F20.6)') ' C(0)A(1)C(1)   ',
     *         IST,OME2,OME2*XTEV
C
               OMECOR = OME1 + OME2
C
               WRITE (LUPRI,'(A)')
     *   ' ----------------------------         --------            --'
     *         //'------'
               WRITE (LUPRI,'(A16,I10,2F20.6)') ' CC(2)          ',
     *         IST,OMECOR,OMECOR*XTEV
               WRITE (LUPRI,'(A/)')
     *   ' ============================         ========            =='
     *         //'======'
C
            ENDIF
C
            IF (CCR3.OR.CCRT.OR.CCR1A.OR.CCR1B) THEN
              ITST = 0
              DO 1310 IOM = 1, NOMINP(ISYMTR,IMULT)
                 IF (IEX .EQ. IOMINP(IOM,ISYMTR,IMULT)) ITST = ITST + 1
 1310         CONTINUE
              IF ((ITST .EQ. 0 ).AND.CCSDT) GOTO 1300
C
              WRITE (LUPRI,'(A)')
     *   ' --------               --------------          ----'
     *         //'----           ---------'
              IF (.NOT. CCR3 ) THEN
               WRITE (LUPRI,'(A16,5X,I10,5X,2F20.6)')' CCSD           ',
     *         IST,WORK(KEIVAL-1+IST),WORK(KEIVAL-1+IST)*XTEV
               WRITE (LUPRI,'(A)')
     *   ' --------               --------------          ----'
     *         //'----           ---------'
              ELSE IF (IPRINT .GT. 10) THEN
               WRITE (LUPRI,'(A16,5X,I10,5X,2F20.6)')' CCSD*(T*ampl.) ',
     *         IST,WORK(KEIVAL-1+IST),WORK(KEIVAL-1+IST)*XTEV
               WRITE (LUPRI,'(A)')
     *   ' --------               --------------          ----'
     *         //'----           ---------'
              ENDIF
C
              IF (CCR3) WRITE (LUPRI,'(A16,5X,I10,5X,2F20.6)')
     *          ' CCSDR(3)       ',IST,OME1,OME1*XTEV
              IF (CCRT) WRITE (LUPRI,'(A16,5X,I10,5X,2F20.6)')
     *          ' CCSDR(T)       ',IST,OME2,OME2*XTEV
              IF (CCR1A) WRITE (LUPRI,'(A16,5X,I10,5X,2F20.6)')
     *          ' CCSDR(1a)       ',IST,OME3,OME3*XTEV
              IF (CCR1B) WRITE (LUPRI,'(A16,5X,I10,5X,2F20.6)')
     *          ' CCSDR(1b)       ',IST,OME4,OME4*XTEV
C
              OMECOR = OME1
C
              WRITE (LUPRI,'(A/)')
     *   ' --------               --------------          ----'
     *         //'----           ---------'
C
            ENDIF
C
            IF (CCP2)  WORK(KEXCPS+IST-1)  = OMECOR
            IF (CCR3)  WORK(KEXCPS+IST-1)  = OME1
            IF (CCRT)  WORK(KEXCPS2+IST-1) = OME2
            IF (CCR1A) WORK(KEXCPS3+IST-1) = OME3
            IF (CCR1B) WORK(KEXCPS4+IST-1) = OME4
C
C-----------------------------------
C           End of loop over states.
C-----------------------------------
C
 1300    CONTINUE
C
         ENDIF
C
C----------------------------------------------
C        End of symmetry and multiplicity loop.
C----------------------------------------------
C
         KEXCPS  = KEXCPS  + NCCEXSAV(ISYM,IMULT)
         KEXCPS2 = KEXCPS2 + NCCEXSAV(ISYM,IMULT)
         KEXCPS3 = KEXCPS3 + NCCEXSAV(ISYM,IMULT)
         KEXCPS4 = KEXCPS4 + NCCEXSAV(ISYM,IMULT)
         KEXCIS  = KEXCIS  + NCCEXSAV(ISYM,IMULT)
         KT1PS   = KT1PS   + NCCEXSAV(ISYM,IMULT)
C
 8000   CONTINUE
C
 9000 CONTINUE
C
      KEXCIS = KEXCI
      KT1PS  = KT1P
C
      IF ((.NOT.CCR3).AND.(LHTR .OR. (.NOT.LIST(1:2).EQ.'LE'))) THEN
C
      WRITE(LURES,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
      WRITE(LURES,'(1X,A26,A10,A)')
     *     '|  sym. | Exci.  |        ',MODELP,' Excitation energies'
     *     //'            | ||T1||  |'
      WRITE(LURES,'(A)')
     *     ' |(spin, |        +-----------------------------'
     *    //'-------------------------------+'
      WRITE(LURES,'(1X,A)')
     *     '| spat) |        |     Hartree    |       eV.      |'
     *     //'     cm-1       |    %    |'
      WRITE(LURES,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      DO 9100 ISYM = 1, NSYM
C
        DO 9050 IMULT = 1, 3, 2
C
         DO 9200 IEX = 1, NCCEXSAV(ISYM,IMULT)
C
            IF ( OMESC ) THEN
               ITST = 0
               DO 9210 IOM = 1, NOMINP(ISYM,IMULT)
                  IF ((IOMINP(IOM,ISYM,IMULT).EQ.IEX).AND.
     *                (.NOT.(OMEINP)))
     *               EOMINP(IOM,ISYM,IMULT)=WORK(KEXCIS+IEX-1)
                  IF (IEX .EQ. IOMINP(IOM,ISYM,IMULT)) ITST = ITST + 1
 9210          CONTINUE
               IF ((ITST .EQ. 0 ).AND.CCSDT) GOTO 9200
            ENDIF
C
c           IF ( IEX .EQ. 1) THEN
               WRITE(LURES,9990) IMULT,REP(ISYM-1),IEX,
     *            WORK(KEXCIS+IEX-1),
     *            WORK(KEXCIS+IEX-1)*XTEV,WORK(KEXCIS+IEX-1)*XTKAYS,
     *            WORK(KT1PS+IEX-1)
c           ELSE
c              WRITE(LURES,9991) IEX,WORK(KEXCIS+IEX-1),
c    *            WORK(KEXCIS+IEX-1)*XTEV,WORK(KEXCIS+IEX-1)*XTKAYS,
c    *            WORK(KT1PS+IEX-1)
c           ENDIF
            OME   = WORK(KEXCIS+IEX-1)
            ETOT  = ECCGRS+OME
            CALL WRIPRO(OME,MOPRPC,0,
     *                  LABEL,LABEL,LABEL,LABEL,
     *                  ECCGRS,OME,ETOT,1,ISYM,IMULT,IEX)
            CALL WRIPRO(ETOT,MOPRPC,0,
     *                  LABEL1,LABEL1,LABEL1,LABEL1,
     *                  ECCGRS,OME,ETOT,1,ISYM,IMULT,IEX)
 9200    CONTINUE
C
         DO 9220 IEX =  1, NCCEXSAV(ISYM,IMULT)
            IF ((ISYM.EQ.IXSTSY).AND.(IEX.EQ.IXSTAT)) THEN
               OMECCX = WORK(KEXCIS+IEX-1)
               ECCXST = ECCGRS + OMECCX
            ENDIF
 9220    CONTINUE
C
         IF (.NOT.(NCCEXSAV(ISYM,IMULT).EQ.0)) THEN
            NREST = 0
            IF (IMULT.EQ.1) NREST = NCCEXCI(ISYM,3)
            DO ISYM2 = ISYM+1,NSYM
              DO IMULT2 = 1, 3, 2
                NREST = NREST + NCCEXCI(ISYM2,IMULT2)
              END DO
            END DO
            IF (NREST.EQ.0) GOTO 9100
            WRITE(LURES,'(A)')
     *        ' +----------------------------------------------'
     *       //'-------------------------------+'
         ENDIF
         KEXCIS = KEXCIS + NCCEXSAV(ISYM,IMULT)
         KT1PS  = KT1PS + NCCEXSAV(ISYM,IMULT)
C
 9050   CONTINUE
 9100 CONTINUE
C
         WRITE(LURES,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      WRITE(LURES,'(//5X,A)') 'Total energies in Hartree:'
       KEXCIS = KEXCI
       KT1PS  = KT1P
       DO ISYM = 1, NSYM
        DO IMULT = 1, 3, 2
         DO IEX = 1, NCCEXSAV(ISYM,IMULT)
            ITST = 1
            IF ( OMESC ) THEN
               ITST = 0
               DO IOM = 1, NOMINP(ISYM,IMULT)
                  IF (IEX .EQ. IOMINP(IOM,ISYM,IMULT)) ITST = ITST + 1
               END DO
            ENDIF
            IF (.NOT.(ITST.EQ.0.AND.CCSDT)) THEN
               WRITE(LURES,9989) IEX,IMULT,REP(ISYM-1),
     *            WORK(KEXCIS+IEX-1)+ECCGRS
            END IF
         END DO
         KEXCIS = KEXCIS + NCCEXSAV(ISYM,IMULT)
         KT1PS  = KT1PS + NCCEXSAV(ISYM,IMULT)
        END DO
       END DO

      ENDIF
C
 9989 FORMAT(10X,I4,' ^',I1,A3,F20.10)
 9990 FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',F13.7,
     *       '  | ',F13.5,'  | ',F13.3,'  | ',F6.2,'  |')
 9991 FORMAT(1X,'|       | ',I4,'   | ',F13.7,
     *       '  | ',F13.5,'  | ',F13.3,'  | ',F6.2,'  |')
C
      IF ((CCP2.OR.((CCR3.AND.CCT).OR.CCRT.OR.CCR1A.OR.CCR1B))
     *   .AND.(LIST(1:2).EQ.'LE')) THEN
C
         KEXCPS  = KEXCP
         KEXCPS2 = KEXCP2
         KEXCPS3 = KEXCP3
         KEXCPS4 = KEXCP4
C
         WRITE(LURES,'(//A/)')
     *     ' Excitation energies with perturbational corrections '
C
C--------------------------------------
C        Loop over triples corrections.
C--------------------------------------
C
         LTRIP(1)   = CCR3
         LTRIP(2)   = CCRT
         LTRIP(3)   = CCR1A
         LTRIP(4)   = CCR1B
C
         CCR3  = .FALSE.
         CCRT  = .FALSE.
         CCR1A = .FALSE.
         CCR1B = .FALSE.
C
         IF ( CCP2 ) THEN
            NTRIP =  1
         ELSE
            NTRIP = 4
         ENDIF
C
         DO 9299 ITRIP = 1, NTRIP
C
            IF ((.NOT. LTRIP(ITRIP)).AND.(.NOT.CCP2)) GOTO 9299
C
            IF ((ITRIP.EQ.1).AND.(.NOT.CCP2))  CCR3   = .TRUE.
            IF (ITRIP.EQ.2)  CCRT   = .TRUE.
            IF (ITRIP.EQ.3)  CCR1A  = .TRUE.
            IF (ITRIP.EQ.4)  CCR1B  = .TRUE.
C
            IF (CCP2 ) MODEL1 = 'CIS(D)'
            IF (CCR3 ) MODEL1 = 'CCSDR(3)'
            IF (CCRT ) MODEL1 = 'CCSDR(T)'
            IF (CCR1A) MODEL1 = 'CCSDR(1a)'
            IF (CCR1B) MODEL1 = 'CCSDR(1b)'

            IF (CCR12) THEN
              MODELSCR = MODEL1 
              CALL CCSD_MODEL(MODEL1,LENMOD,24,MODELSCR,24,APROXR12)
            END IF
C
            WRITE(LURES,'(/A)')
     *        ' +=============================================='
     *       //'=====================+'
            WRITE(LURES,'(1X,A26,A10,A)')
     *       '|  sym. | Exci.  |        ',MODEL1,' Excitation energies'
     *       //'            |'
            WRITE(LURES,'(A)')
     *        ' |(spin, |        +-----------------------------'
     *       //'---------------------+'
            WRITE(LURES,'(1X,A)')
     *        '| spat) |        |     Hartree    |       eV.      |'
     *        //'     cm-1       |'
            WRITE(LURES,'(A)')
     *        ' +=============================================='
     *       //'=====================+'
C
            DO 9300 ISYM = 1, NSYM
C
C            DO 9310 IMULT = 1, 3, 2
             DO 9310 IMULT = 1, 1, 2
C
               DO 9400 IEX = 1, NCCEXCI(ISYM,IMULT)
C
                  IF (CCP2)  OME = WORK(KEXCPS +IEX-1)
                  IF (CCR3)  OME = WORK(KEXCPS +IEX-1)
                  IF (CCRT)  OME = WORK(KEXCPS2+IEX-1)
                  IF (CCR1A) OME = WORK(KEXCPS3+IEX-1)
                  IF (CCR1B) OME = WORK(KEXCPS4+IEX-1)
C
                  ITST = 0
                  DO 9410 IOM = 1, NOMINP(ISYM,IMULT)
                     IF ((IOMINP(IOM,ISYM,IMULT).EQ.IEX).AND.
     *                   (.NOT.(OMEINP)).AND.OMESC)
     *                   EOMINP(IOM,ISYM,IMULT)=OME
                     IF (IEX .EQ. IOMINP(IOM,ISYM,IMULT)) ITST=ITST+1
 9410             CONTINUE
C
                  IF ((ITST .EQ. 0 ).AND.(.NOT.CCP2)) GOTO 9400
C
c                 IF ( IEX .EQ. 1) THEN
                     WRITE(LURES,9992) IMULT,REP(ISYM-1),IEX,OME,
     *                  OME*XTEV,OME*XTKAYS
c                 ELSE
c                    WRITE(LURES,9993) IEX,OME,
c    *                  OME*XTEV,OME*XTKAYS
c                 ENDIF
                  ETOT  = ECCGRS+OME
                  CALL WRIPRO(OME,MODEL1,0,
     *                        LABEL,LABEL,LABEL,LABEL,
     *                        0.0D0,0.0D0,0.0D0,1,ISYM,IMULT,IEX)
                  CALL WRIPRO(ETOT,MODEL1,0,
     *                        LABEL1,LABEL1,LABEL1,LABEL1,
     *                        0.0D0,0.0D0,0.0D0,1,ISYM,IMULT,IEX)

                  IF ((ISYM .EQ. IXSTSY).AND.(IEX.EQ.IXSTAT)) THEN
                     OMECCX = OME
                     ECCXST = ECCGRS + OMECCX
                  ENDIF
C
 9400          CONTINUE
C
               IF (.NOT.((ISYM .EQ. NSYM).OR.
     *                (NCCEXCI(ISYM,IMULT).EQ.0))) THEN
                 NREST = 0
                 DO 9350 ISYM2 = ISYM+1,NSYM
                    NREST = NREST + NCCEXCI(ISYM2,IMULT)
 9350            CONTINUE
                 IF (NREST.EQ.0) GOTO 9300
                 WRITE(LURES,'(A)')
     *              ' +----------------------------------------------'
     *             //'---------------------+'
               ENDIF
C
               IF (CCP2)  KEXCPS  = KEXCPS  + NCCEXCI(ISYM,IMULT)
               IF (CCR3)  KEXCPS  = KEXCPS  + NCCEXCI(ISYM,IMULT)
               IF (CCRT)  KEXCPS2 = KEXCPS2 + NCCEXCI(ISYM,IMULT)
               IF (CCR1A) KEXCPS3 = KEXCPS3 + NCCEXCI(ISYM,IMULT)
               IF (CCR1B) KEXCPS4 = KEXCPS4 + NCCEXCI(ISYM,IMULT)
C
 9310        CONTINUE
 9300       CONTINUE
C
            WRITE(LURES,'(A)')
     *        ' +=============================================='
     *       //'=====================+'
C
             KEXCPS  = KEXCP
             KEXCPS2 = KEXCP2
             KEXCPS3 = KEXCP3
             KEXCPS4 = KEXCP4
             WRITE(LURES,'(//5X,A)') 'Total energies in Hartree:'
             DO ISYM = 1, NSYM
C             DO IMULT = 1, 3, 2
              DO IMULT = 1, 1, 2
               DO IEX = 1, NCCEXCI(ISYM,IMULT)
                  IF (CCP2)  OME = WORK(KEXCPS +IEX-1)
                  IF (CCR3)  OME = WORK(KEXCPS +IEX-1)
                  IF (CCRT)  OME = WORK(KEXCPS2+IEX-1)
                  IF (CCR1A) OME = WORK(KEXCPS3+IEX-1)
                  IF (CCR1B) OME = WORK(KEXCPS4+IEX-1)
                  ITST = 0
                  DO IOM = 1, NOMINP(ISYM,IMULT)
                     IF (IEX .EQ. IOMINP(IOM,ISYM,IMULT)) ITST=ITST+1
                  END DO
                  IF (.NOT.(ITST.EQ.0 .AND. .NOT.CCP2)) THEN
                     WRITE(LURES,9989) IEX,IMULT,REP(ISYM-1),
     *                  OME+ECCGRS
                  END IF
               END DO
               IF (CCP2)  KEXCPS  = KEXCPS  + NCCEXCI(ISYM,IMULT)
               IF (CCR3)  KEXCPS  = KEXCPS  + NCCEXCI(ISYM,IMULT)
               IF (CCRT)  KEXCPS2 = KEXCPS2 + NCCEXCI(ISYM,IMULT)
               IF (CCR1A) KEXCPS3 = KEXCPS3 + NCCEXCI(ISYM,IMULT)
               IF (CCR1B) KEXCPS4 = KEXCPS4 + NCCEXCI(ISYM,IMULT)
              END DO
             END DO
C
            IF (CCR3)  CCR3  = .FALSE.
            IF (CCRT)  CCRT  = .FALSE.
            IF (CCR1A) CCR1A = .FALSE.
            IF (CCR1B) CCR1B = .FALSE.
C
 9299    CONTINUE
C
         CCR3   = LTRIP(1)
         CCRT   = LTRIP(2)
         CCR1A  = LTRIP(3)
         CCR1B  = LTRIP(4)
C
 9992    FORMAT(1X,'| ^',I1,A3,' | ',I4,'   | ',F13.7,
     *       '  | ',F13.5,'  | ',F13.3,'  |')
 9993    FORMAT(1X,'|       | ',I4,'   | ',F13.7,
     *       '  | ',F13.5,'  | ',F13.3,'  |')
C
      ENDIF
C
      IF (IPRINT .GT.10) THEN
         CALL AROUND( ' END OF CC_EXCI   ' )
      ENDIF
C
      CALL QEXIT('CC_EXCI')
C
      RETURN
      END
c*DECK CCLR_LEINFI
      SUBROUTINE CCLR_LEINFI(TRIPLET)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Purpose: set common Leinf for calculation of
C              Coupled Cluster excitation energies.
C
C     Written by Ove Christiansen November 1994
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <leinf.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
C
      LOGICAL TRIPLET
C
      CALL QENTER('CCLR_LEINFI')
C
C--------------------------
C     Set COMMON /LEINF/
C--------------------------
C
      IF (IPRINT .GT.10) THEN
         CALL AROUND( ' START OF CCLR_LEINFI ' )
      ENDIF
C
      THRLE  = THREXC
      MAXLE  = MAXITE
      IPRLE  = IPRINT
C
Cholesky
C
      IF (CHOINT) THEN
         IF (TRIPLET)
     &   CALL QUIT('CCLR_LEINFI: Triplet Cholesky not implemented')
         IF (CC2) THEN
            LETYPA = NT1AM(ISYMTR)
         ELSE
            CALL QUIT('CCLR_LEINFI: Only Cholesky CC2 implemented')
         ENDIF
      ELSE
         LETYPA = NT1AM(ISYMTR) + NT2AM(ISYMTR)
         IF (  TRIPLET  ) LETYPA = LETYPA + NT2AMA(ISYMTR)
         IF (   CCR12   ) LETYPA = LETYPA + NTR12AM(ISYMTR)
         IF (CCS.OR.CCP2) LETYPA = NT1AM(ISYMTR)
      END IF
C
      LETYPB = 0
      LETOT  = LETYPA + LETYPB
      NBASIS = 1
      NCREF  = 0
      LULEA  = 45
      LULEB  = 46
      LULEC  = 47
C
C--------------------------------------------------------------------
C     B space has zero dimension.
C     A space : integral distribution, and vectors and intermediates.
C--------------------------------------------------------------------
C
      LLEWB = 0
      LLEWA = 0
      NCCEX = 0
C
      DO 10 ISYM = 1, NSYM
         IF (NDISAO(ISYM).GT. LLEWA) LLEWA = NDISAO(ISYM)
  10  CONTINUE
C
      NRHO2 = 2*NT2ORT(1) + NT2SQ(1)
      IF (CCS) NRHO2 = 0
C
      LLEWA = LLEWA + NRHO2
C
      LLEWA = LLEWA + (10+NRHFT)*N2BST(1)
      LLEWA = LLEWA + 2*NT2BGD(1)
      LLEWA = LLEWA + NDSRHF(1)
      LLEWA = LLEWA + 2*NT2MAO(1,1)
C
      IF (CCSDT) LLEWA = LLEWA  +  NT2SQ(1)
C
Cholesky
      IF (CHOINT .AND. CC2) THEN
         NRHO2 = 0
         LLEWA = 0
      END IF
Cholesky
C
C---------------------------------------------------------------
C     t2tcor is taken care of in routine, that is; it is forced
C     to be false if there is problems.
C---------------------------------------------------------------
C
C     IF (T2TCOR) LLEWA = LLEWA + NT2SQ(1)
C
      IF (IPRINT .GT.60) THEN
         CALL AROUND('COMMON CCLR in CCEXCI')
         WRITE (LUPRI,'(A32,F24.12)') 'IN CCLR_LEINFI:  THREXC', THREXC
         IMULT = 1
         IF (TRIPLET) IMULT = 3
         WRITE(LUPRI,1) 'NCCEXCI:',(NCCEXCI(I,IMULT),  I=1,NSYM)
      ENDIF
      IF (IPRINT .GT.60) THEN
         CALL AROUND('COMMON /LEINF/ in CCEXCI')
         WRITE (LUPRI,'(A32,F24.12)')  'IN CCLR_LEINFI:  THRLE ', THRLE
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  MAXLE ', MAXLE
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  IPRLE ', IPRLE
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  LETYPA', LETYPA
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  LETYPB', LETYPB
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  LLEWA ', LLEWA
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  LLEWB ', LLEWB
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  LETOT ', LETOT
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  NSIDE ', NSIDE
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  NCREF ', NCREF
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  LULEA ', LULEA
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  LULEB ', LULEB
         WRITE (LUPRI,'(A32,I12)')     'IN CCLR_LEINFI:  LULEC ', LULEC
      END IF
C
      IF (IPRINT .GT.10) THEN
         CALL AROUND( ' END OF CCLR_LEINFI ' )
      ENDIF
C
    1 FORMAT(3X,A8,8I8)
C
      CALL QEXIT('CCLR_LEINFI')
C
      END
c*DECK CC_PCEXCI 
      SUBROUTINE CC_PCEXCI(OMES,OME1,OME2,OME3,OME4,WORK,LWORK)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Purpose: 
C              Direct calculation of perturbative corrections
C              to Coupled Cluster excitation energies.
C
C              CC(2) correction to CCS  - identical to CIS(D)
C              (CCS =TDA=CIS ),(CC(2)=CIS(D))
C
C              CCSDR(3),CCSDR(1a),CCSDR(1b),CCSDR(T)
C
C     Written by Ove Christiansen 10-3-1995 
C                                 8-10-1995 
C                                 26-2-1996
C                                 10-12-1996
C                                 3-2-1997
C     Version 3, February 1996
C                                 
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <ccexci.h>
C
      PARAMETER(NTRIP = 4,MXISPT = 3)
C
      LOGICAL TRIPLET
      LOGICAL CCSSAV, CC2SAV, CC3SAV, LINQCC, LTRIP(NTRIP),LHSOLD
      DIMENSION OME1(*),OME2(*),OME3(*),OME4(*),WORK(LWORK),OMES(*)
      DIMENSION IADR(MXISPT)
C
      CHARACTER CHSYM*4,MODEL1*10,MODEL*10,APROXR12*3
      CHARACTER LBLPT(MXISPT)*8
      CHARACTER*8 FC1AM,FC2AM,FRHO1,FRHO2
      PARAMETER (TWO = 2.0D00,TF=1.0D-04, TLOVLP = 0.5D0 )
      PARAMETER (FC1AM='CCR_C1AM',FC2AM='CCR_C2AM')
      PARAMETER (FRHO1='CCR_RHO1',FRHO2='CCR_RHO2')
C
#include <leinf.h>
C
      CALL QENTER('CC_PCEXCI')
C
      LUFC1 = -1
      LUFC2 = -1
      LUFR1 = -1
      LUFR2 = -1
C
      WRITE (LUPRI,'(1X,A,/)') '  '
      WRITE (LUPRI,'(1X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      IF ( CCR1A.OR.CCR1B.OR.CCRT.OR.CCR3.OR.CCP2) THEN
         WRITE (LUPRI,'(1X,A)')
     *   '*                                                        '//
     *   '         *'
         WRITE (LUPRI,'(1X,A)')
     *   '*<<<<<<<<<<  CALCULATION OF PERTURBATIONAL CORRECTIONS  >'//
     *   '>>>>>>>>>*'
      ENDIF
      WRITE (LUPRI,'(1X,A)')
     *'*                                                        '//
     *'         *'
      WRITE (LUPRI,'(1X,A,/)')
     *'*********************************************************'//
     *'**********'
C
      TRIPLET = .FALSE.
      IMULT   = 1
C
      MODEL = 'CCSD      '
      IF (CCS)   MODEL = 'CC2       '
      IF (CCP2)  MODEL = 'CC(2)     '
      IF (CCR1A) MODEL = 'CCSDR(1A) '
      IF (CCR1B) MODEL = 'CCSDR(1B) '
      IF (CCRT)  MODEL = 'CCSDR(T) '
      IF (CCR3)  MODEL = 'CCSDR(3) '
C
      IF ( CCR3.OR.CCRT.OR.CCR1A.OR.CCR1B ) THEN
         CC3SAV = CCSDT
         CCSDT  = .FALSE.
      ENDIF
C
      IF ( CCP2 ) THEN
         CCSSAV = CCS
         CC2SAV = CC2
         CCS    = .FALSE.
         CC2    = .TRUE.
      ENDIF
C
      IF (CCP2.OR.((CCR3.OR.CCRT).OR.(CCR1A.OR.CCR1B))) THEN
C
         WRITE (LUPRI,'(//A)')
     *   ' ========================================================'
     *   //'=============='
         WRITE (LUPRI,'(A)')
     *   ' ###                                                   '
     *   //'             ###'
         WRITE (LUPRI,'(A)')
     *   ' ###                                                   '
     *   //'             ###'
         WRITE (LUPRI,'(A)')
     *   ' ###      Perturbational Corrections to Excitation '
     *   //' energies.       ###'
         WRITE (LUPRI,'(A)')
     *   ' ###                                                   '
     *   //'             ###'
         WRITE (LUPRI,'(A)')
     *   ' ###                                                   '
     *   //'             ###'
         IF (CCR3.OR.CCRT.OR.CCR1A.OR.CCR1B) THEN
            WRITE(LUPRI,'(A)')
     *   ' ###    Calculating triples corrected CCSD excitation e'
     *   //'nergies.     ###'
         ENDIF
         IF (CCP2) THEN
            WRITE(LUPRI,'(A)')
     *   ' ###     Calculating doubles corrected CCS excitation e'
     *   //'nergies.     ###'
         ENDIF
         WRITE (LUPRI,'(A)')
     *   ' ###                                                   '
     *   //'             ###'
         WRITE (LUPRI,'(A)')
     *   ' ###                                                   '
     *   //'             ###'
         WRITE (LUPRI,'(A)')
     *   ' ========================================================'
     *   //'=============='
C
      ENDIF
C
      IF (CCR3.OR.CCRT.OR.CCR1A.OR.CCR1B ) THEN
C
C--------------------------
C        Work allocation 2.
C--------------------------
C
         NTAMP = NT1AM(ISYMTR) + NT2AM(ISYMTR)
C
         KRV1   = 1    
         KEND1  = KRV1  + NTAMP
         LEND1  = LWORK - KEND1
C
         IF (LEND1.LE.0) CALL QUIT('Too little work space in cc_pc')
C
         DO 50 IV = 1, NCCEXCI(ISYMTR,IMULT)
C
            ISTATE = ISYOFE(ISYMTR) + IV    
            EIGV   = EIGVAL(ISTATE)
            IF (IPRINT .GT. 5) THEN
               WRITE(LUPRI,'(/,1x,A,I3,/1X,A,I3,A,F16.8)')
     *         'Calculating triples corrections for state ',ISTATE,
     *         'of symmetry ',ISYMTR,' and eigenvalue: ',EIGV
            ENDIF
C
C------------------------------------------------------------
C           Calculate (A(0)+A(1)*(A(0)-E(ccsd)1)-1*A(1))*C(0)
C           Loop over various triple correction approaches.
C------------------------------------------------------------
C
            ECURR = EIGV
C
            ITST = 0
            DO 55 IOM = 1, NOMINP(ISYMTR,IMULT)
               IF (IOMINP(IOM,ISYMTR,IMULT).EQ.IV) ITST = ITST + 1
   55       CONTINUE
C
            IF (ITST .EQ. 0 ) THEN
               IF (CCR3)  OME1(IV) = 0.0
               IF (CCRT)  OME2(IV) = 0.0
               IF (CCR1A) OME3(IV) = 0.0
               IF (CCR1B) OME4(IV) = 0.0
               GOTO 50
            ENDIF
C
            IF ( IPRINT. GT.  5) THEN
               WRITE(LUPRI,'(A,F10.6)') ' Doing partioned triples '
     *             //'linear transformation with ECURR= ',ECURR
            ENDIF
C
C--------------------------------------------------
C           loop over different triple corrections.
C--------------------------------------------------
C
            LTRIP(1) = CCR3
            LTRIP(2) = CCRT
            LTRIP(3) = CCR1A
            LTRIP(4) = CCR1B
C
            CCR3     = .FALSE.
            CCRT     = .FALSE.
            CCR1A    = .FALSE.
            CCR1B    = .FALSE.
C
            DO 60 ITRIP = 1, NTRIP
C
             IF (.NOT. LTRIP(ITRIP)) GOTO 60
C
             CCSDT = .TRUE.
             IF (ITRIP.EQ.1)  CC3   = .TRUE.
             IF (ITRIP.EQ.2)  CCRT  = .TRUE.
             IF (ITRIP.EQ.3)  CC1A  = .TRUE.
             IF (ITRIP.EQ.4)  CC1B  = .TRUE.
C
             IF ( IPRINT .GT. 0 ) THEN
                WRITE(LUPRI,'(/,A,F10.6)') ' ECURR: ',ECURR
                WRITE(LUPRI,'(A,F10.6)') ' L*R Norm is ',
     &               XNORM(ISTATE)**2
             ENDIF
             IF ( IPRINT .GT. 5 ) THEN
                WRITE(LUPRI,*) 'CC3,CC1A,CC1B,CCRT,CC3LR: ',
     *                   CC3,CC1A,CC1B,CCRT,CC3LR
             ENDIF
C
C------------------------------------------------------------
C            Readin right solution and find excitation energy.
C------------------------------------------------------------
C
             IOPT   = 3
             CALL CC_RDRSP('RE',ISTATE,ISYMTR,IOPT,MODEL,WORK(KRV1),
     *                     WORK(KRV1+NT1AM(ISYMTR)))
             IF (IPRINT .GT. 40 ) THEN
                CALL AROUND( 'In CC_PCEXCI:  right eigen vector ' )
                CALL CC_PRP(WORK(KRV1),WORK(KRV1+NT1AM(ISYMTR)),
     *                      ISYMTR,1,1)
             ENDIF
C
C--------------------------------------------------------------------
C            Calculate linear transformation.
C            Input vector is first elements in work - 
C            Output vector replaces vector as first elements in work.
C--------------------------------------------------------------------
C
             CALL CC_ATRR(ECURR,ISYMTR,1,WORK,LWORK,.FALSE.,DUMMY,
     &                     APROXR12)
C
             IF (IPRINT .GT.25) THEN
                CALL AROUND('CC_PCEXCI: partioned A*R(0) vector. ')
                CALL CC_PRP(WORK(KRV1),WORK(KRV1+NT1AM(ISYMTR)),
     *                      ISYMTR,1,1)
             ENDIF
C
C--------------------------------
C            Readin left solution.
C--------------------------------
C
             KLV1   = KEND1
             KEND2  = KLV1  + NTAMP
             LEND2  = LWORK - KEND2
             IF (LEND2.LE.0) CALL QUIT('Too little work '//
     &            'space in cc_pc-2')
C           
             IOPT   = 3
             CALL CC_RDRSP('LE',ISTATE,ISYMTR,IOPT,MODEL,WORK(KLV1),
     *                     WORK(KLV1+NT1AM(ISYMTR)))
             IF (IPRINT .GT. 40 ) THEN
                CALL AROUND( 'In CC_PCEXCI:  Left Eigen vector ' )
                CALL CC_PRP(WORK(KLV1),WORK(KLV1+NT1AM(ISYMTR)),
     *                      ISYMTR,1,1)
             ENDIF
C
             IF (CC3) OME1(IV) = 
     *          DDOT(NTAMP,WORK(KLV1),1,WORK(KRV1),1)/XNORM(ISTATE)**2
             IF (CCRT) OME2(IV) = 
     *          DDOT(NTAMP,WORK(KLV1),1,WORK(KRV1),1)/XNORM(ISTATE)**2
             IF (CC1A) OME3(IV) = 
     *          DDOT(NTAMP,WORK(KLV1),1,WORK(KRV1),1)/XNORM(ISTATE)**2
             IF (CC1B) OME4(IV) = 
     *          DDOT(NTAMP,WORK(KLV1),1,WORK(KRV1),1)/XNORM(ISTATE)**2
C
             IF ( IPRINT .GT. 1) THEN
                IF (CC3) WRITE(LUPRI,'(1X,A9,I3,1X,A15,1X,F10.6)')
     *                   'Exci. nr.',IV,'in CCSDR(3) is',OME1(IV)
                IF (CCRT) WRITE(LUPRI,'(1X,A9,I3,1X,A15,1X,F10.6)')
     *                   'Exci. nr.',IV,'in CCSDR(T) is',OME2(IV)
                IF (CC1A) WRITE(LUPRI,'(1X,A9,I3,1X,A15,1X,F10.6)')
     *                   'Exci. nr.',IV,'in CCSDR(1A) is',OME3(IV)
                IF (CC1B) WRITE(LUPRI,'(1X,A9,I3,1X,A15,1X,F10.6)')
     *                   'Exci. nr.',IV,'in CCSDR(1B) is',OME4(IV)
             ENDIF
             CCSDT = .FALSE.
             IF (CC3)  CC3   = .FALSE.
             IF (CCRT) CCRT  = .FALSE.
             IF (CC1A) CC1A  = .FALSE.
             IF (CC1B) CC1B  = .FALSE.
C
   60       CONTINUE
C
            CCR3     = LTRIP(1)
            CCRT     = LTRIP(2)
            CCR1A    = LTRIP(3)
            CCR1B    = LTRIP(4)
C
   50    CONTINUE
C
      ENDIF
C
C-------------------
C     CC(2) section.
C   NB!! kan laves smatere ved at lave logisk flag
C   saa kun C1 transformeres altsaa spare allokering
C   til C2SQ.
C-------------------
C
      IF (CCP2) THEN
C
C--------------------------
C        Work allocation 3.
C--------------------------
C
         NTAMP = NT1AM(ISYMTR) + NT2AM(ISYMTR)
         NT    = NT1AM(ISYMTR)
C
         DO 100 IV = 1, NCCEXCI(ISYMTR,IMULT)
C
           ISTATE = ISYOFE(ISYMTR) + IV    
           EIGV   = EIGVAL(ISTATE)
           IF (IPRINT .GT. 5) THEN
              WRITE(LUPRI,'(/,1x,A,I3,/1X,A,I3,A,F16.8)')
     *        'Calculating doubles corrections for state ',ISTATE,
     *        'of symmetry ',ISYMTR,' and eigenvalue: ',EIGV
           ENDIF
C
           KRV1   = 1    
           KEND1  = KRV1  + NTAMP
           LEND1  = LWORK - KEND1
           IF (LEND1.LE.0) CALL QUIT('Too little work space in cc_pc')
           CALL DZERO(WORK(KRV1),NTAMP)
C
           IOPT   = 1
           CALL CC_RDRSP('RE',ISTATE,ISYMTR,IOPT,MODEL,WORK(KRV1),
     *                   WORK(KRV1+NT1AM(ISYMTR)))
           IF (IPRINT .GT.15) THEN
              CALL AROUND('CC_PCEXCI: C(0) vector. ')
              CALL CC_PRP(WORK(KRV1),WORK(KRV1+NT1AM(ISYMTR)),
     *                    ISYMTR,1,1)
           ENDIF
C
           CALL CC_ATRR(ECURR,ISYMTR,1,WORK,LWORK,.FALSE.,DUMMY,
     &                     APROXR12)
           IF (IPRINT .GT.15) THEN
              CALL AROUND('CC_PCEXCI: A*C(0) vector. ')
              CALL CC_PRP(WORK(KRV1),WORK(KRV1+NT1AM(ISYMTR)),
     *                    ISYMTR,1,1)
           ENDIF
C
C-----------------------------------------------------------------------
C          Scale Vector with orbital energy diff. and exci. energy.
C-----------------------------------------------------------------------
C
           KLV1   = KEND1
           KEND2  = KLV1  + NT
           LEND2  = LWORK - KEND2
           IF (LEND2.LE.0) CALL QUIT('Too little work space in cc_pc-2')
           IOPT   = 1
           CALL CC_RDRSP('LE',ISTATE,ISYMTR,IOPT,MODEL,WORK(KLV1),
     *                   WORK(KRV1+NT1AM(ISYMTR)))
           OME1(IV) = DDOT(NT1AM(ISYMTR),WORK(KLV1),1,WORK(KRV1),1)
           CALL CC_OMEC(OME2(IV),WORK(KRV1+NT1AM(ISYMTR)),EIGV,
     *                   WORK(KEND1),LEND1,ISYMTR)
  100   CONTINUE
C
      ENDIF
C
      IF (CCP2) THEN
C
         CCS    = CCSSAV
         CC2    = CC2SAV
C
      ENDIF
C
      IF ( CCR3.OR.CCRT.OR.CCR1A.OR.CCR1B )  THEN
         CC3SAV = .FALSE.
         CCSDT  = CC3SAV
      ENDIF
C
      CALL QEXIT('CC_PCEXCI')
C
      END
C  /* Deck cc_vscal */
      SUBROUTINE CC_VSCAL(OMEGA1,OMEGA2,OME,WORK,LWORK,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Ove Christiansen 10-5-1995
C
C     Purpose: Scale OMEGA with diagonals.
C              eps(a) - eps(i)
C              eps(a) + eps(b) - eps(i) - eps(j)
C              Used for calculating pert. corr. amplitudes in
C              CCSDR(3)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <dummy.h>
C
      DIMENSION OMEGA1(*),OMEGA2(*),WORK(LWORK)
C
#include <priunit.h>
#include <inftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
Cholesky
#include <maxorb.h>
#include <ccdeco.h>
Cholesky
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_VSCAL')
C
C-----------------------
C     Memory allocation.
C-----------------------
C
      KSCR1 = 1
      KEND  = KSCR1 + NORBTS
      LWRK  = LWORK - KEND
C
      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient space in CC2_FCK')
      ENDIF
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(I), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KSCR1),WORK(KEND),LWRK)
C
      IF (IPRINT .GT. 80) THEN
         CALL AROUND('CC_VSCAL - Orbital energies. ')
         WRITE (LUPRI,*) (WORK(I), I=1,NORBTS)
         CALL AROUND('CC_VSCAL - start - : RHO1,RHO2 ')
         CALL CC_PRP(OMEGA1,OMEGA2,ISYMTR,1,1)
      ENDIF
C
C----------------------
C     Transform vector.
C----------------------
C
      DO 10 ISYMI = 1, NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMTR)
C
         DO 20 I = 1, NRHF(ISYMI)
C
            MI = IORB(ISYMI) + I
C
            DO 30 A = 1, NVIR(ISYMA)
C
               MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
               NAI = IT1AM(ISYMA,ISYMI)
     *             + NVIR(ISYMA)*(I - 1) + A
C
               DEN =  (WORK(MA) - WORK(MI) - OME  )
C
               XIDEN = 1/DEN
C
               OMEGA1(NAI) = OMEGA1(NAI)*XIDEN
C
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
C
C     Skip doubles part if CCS or Cholesky CC2
C
      IF (.NOT. (CCS .OR. (CHOINT .AND. CC2))) THEN
        DO 100 ISYMBJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMBJ,ISYMTR)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  MJ = IORB(ISYMJ) + J
C
                  DO 140 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
C
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        MI = IORB(ISYMI) + I
C
                        DO 160 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
C
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                           IF (((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .GT. NBJ)).OR.(ISYMAI.GT.ISYMBJ))
     *                          GOTO 160
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ENDIF
C
                           DEN =  (WORK(MA) + WORK(MB)
     *                         - WORK(MI) - WORK(MJ) - OME  )
C
                           XIDEN = 1/DEN
C
                           OMEGA2(NAIBJ) = OMEGA2(NAIBJ)*XIDEN
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100   CONTINUE
      ENDIF
C
      IF (IPRINT .GT. 80) THEN
         IF (CCS .OR. (CHOINT.AND.CC2)) THEN
            I = 0
            CALL AROUND('CC_VSCAL - end - : RHO1')
         ELSE
            I = 1
            CALL AROUND('CC_VSCAL - end - : RHO1,RHO2')
         ENDIF
         CALL CC_PRP(OMEGA1,OMEGA2,ISYMTR,1,I)
      ENDIF
C
      CALL QEXIT('CC_VSCAL')
C
      RETURN
      END
C  /* Deck cc_omec */
      SUBROUTINE CC_OMEC(OMEC,OMEGA2,OME,WORK,LWORK,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Ove Christiansen 10-5-1995
C
C     Purpose: Used in CC(2) for scaling
C              in this case OME is different
C              from zero.
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
      DIMENSION OMEGA2(*),WORK(LWORK)
C
#include <inftap.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
      CALL QENTER('CC_OMEC')
C
C-----------------------
C     Memory allocation.
C-----------------------
C
      KSCR1 = 1
      KEND  = KSCR1 + NORBTS
      LWRK  = LWORK - KEND
C
      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient space in CC_OMEC')
      ENDIF
C
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
C
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
C
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(I), I=1,NORBTS)
C
      CALL GPCLOSE(LUSIFC,'KEEP')
C
      IF (FROIMP .OR. FROEXP)
     *   CALL CCSD_DELFRO(WORK(KSCR1),WORK(KEND),LWRK)
C
C----------------------------
C     Calculate contribution.
C----------------------------
C
      OMEC = 0.0D00
C
      DO 100 ISYMBJ = 1,NSYM
C
         ISYMAI = MULD2H(ISYMBJ,ISYMTR)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 120 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
               ISYMBI = MULD2H(ISYMB,ISYMI)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  MJ = IORB(ISYMJ) + J
C
                  DO 140 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
C
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
C
                     DO 150 I = 1,NRHF(ISYMI)
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
C
                        MI = IORB(ISYMI) + I
C
                        DO 160 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
C
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                           IF (((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .GT. NBJ)).OR.(ISYMAI.GT.ISYMBJ))
     *                          GOTO 160
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ENDIF
C
                           IF (ISYMAJ .EQ. ISYMBI) THEN
C
                              NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                              + INDEX(NAJ,NBI)
C
                           ELSE IF (ISYMAJ .LT. ISYMBI) THEN
C
                              NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                              + NT1AM(ISYMAJ)*(NBI - 1) + NAJ
C
                           ELSE IF (ISYMBI .LT. ISYMAJ) THEN
C
                              NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                              + NT1AM(ISYMBI)*(NAJ - 1) + NBI
C
                           ENDIF
C
                           DEN = - (  WORK(MA) + WORK(MB)
     *                              - WORK(MI) - WORK(MJ) - OME  )
C
                           XIDEN = 1/DEN
                           XIDEN2 = XIDEN
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              IF (NAI.EQ.NBJ) XIDEN2 = XIDEN2*2
                           ENDIF
C
                           OMEC = OMEC +
     *                            (2*OMEGA2(NAIBJ)-OMEGA2(NAJBI))*
     *                             OMEGA2(NAIBJ)*XIDEN2
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      CALL QEXIT('CC_OMEC')
C
      RETURN
      END
c*DECK CC_ATRR
      SUBROUTINE CC_ATRR(ECURR,ISYMV,ISIDE,WORK,LWORK,LCONT,CONT,
     &                   APROXR12)
C
C----------------------------------------------------------------------
C     Ove Christiansen December 1996.
C
C     Jacobian transformation with ONE vector.
C     Vector is first element in WORK and on output is replaced
C     by its linear transformed.
C     For CCR12 on input R_12 is passed after the conventional trial vector
C     and while on output A * R_12 and S * R_12 passed at this place
C
C     removed problem with CCS left transformation, C.H., October 1997
C----------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <iratdef.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <ccsdinp.h>
#include <r12int.h>
C
      LOGICAL LCONT,LOCDBG
      CHARACTER*8 FC1AM,FC2AM,FRHO1,FRHO2,FR2SD,FRHO12,FC12AM,FS12AM,
     &            FS2AM
      PARAMETER (FC1AM ='CCR_C1AM',FC2AM ='CCR_C2AM',FC12AM='CCR_C12M')
      PARAMETER (FRHO1 ='CCR_RHO1',FRHO2 ='CCR_RHO2',FRHO12='CCR_RH12')
      PARAMETER (                  FS2AM ='CCR_S2DM',FS12AM='CCR_S12M')
      PARAMETER (TWO = 2.0D00 )
      PARAMETER (LOCDBG = .FALSE.)
      CHARACTER*3 APROXR12
C
      DIMENSION WORK(LWORK), CONT(28)
C
      CALL QENTER('CC_ATRR')
      IF (CCS) THEN
         NCCVAR = NT1AM(ISYMV)
      ELSE
         NCCVAR = NT1AM(ISYMV) + NT2AM(ISYMV)
         IF (CCR12) NCCVAR = NCCVAR + NTR12AM(ISYMV)
      END IF

C     IF (CCS .OR. CC2 .OR. CCSD) THEN
C       ECURR = 0.0D0
C     ELSE
C       CALL QUIT('Entered for unknown model CC_ATRR')
C     END IF
C
      NSIMTR = 1
      ISYMTR = ISYMV
C
      LUFR1  = -1
      LUFR2  = -1
      LUFR12 = -1
      LUFC1  = -1
      LUFC2  = -1
      LUFC12 = -1
      LUFS12 = -1
      LUFS2  = -1
C
C----------------
C     Open files.
C----------------
C
      CALL CC_FILOP(FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *              FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *              FS12AM,LUFS12,FS2AM,LUFS2)
C
C----------------------------------------------------------------------
C     Make rho2 file name.
C     For CCSD rho2 has to be stored on different file 
C     due to different length.
C----------------------------------------------------------------------
C
      IF (.NOT. (CCS.OR.CC2)) THEN
         LUFSD = -1
         FR2SD = 'CC_TRA2_'
      ELSE
         LUFSD = LUFR2
         FR2SD = FRHO2
      ENDIF
C
C-----------------------------------------------------------------------
C     Cheat and do CCS left transformation by right hand transformation.
C-----------------------------------------------------------------------
C
      MYSIDE = ISIDE
C     IF ((ISIDE .EQ. -1 ) .AND. CCS ) ISIDE = 1
      IF ((ISIDE .EQ. -1 ) .AND. CCS ) MYSIDE = 1

C
C-----------------------------------------------------------------------
C
      K1   = 1
      IVEC = K1
      IF (.NOT. (CCS .OR. CC2)) THEN
         ITR  = 1
      ELSE
         ITR  = K1
      ENDIF
C
C--------------------
C     Allocations.
C--------------------
C
      NRHO2 = MAX(NT2AM(ISYMTR),NT2AM(1),2*NT2ORT(ISYMTR))
CCH   IF (ISIDE  .EQ. -1) NRHO2 = MAX(NRHO2,2*NT2ORT(1))
      IF (MYSIDE .EQ. -1) NRHO2 = MAX(NRHO2,2*NT2ORT(1))

      IF ( CC2 ) NRHO2 = MAX(NT2AM(ISYMTR),NT2AM(1))
      IF ( CCS ) NRHO2 = 2
C
      NC2AM = MAX(NT2SQ(ISYMTR),NT2SQ(1),
     *           NT2AM(ISYMTR)+2*NT2ORT(1),NT2R12(1))
      IF ( CC2 ) NC2AM = MAX(NT2SQ(ISYMTR),NT2SQ(1))
      IF ( CCS ) NC2AM = 2
C
      NRHO1 = NT1AM(ISYMTR)*NSIMTR
C
      KRHO1 = 1
      KRHO2 = KRHO1 + NRHO1
      KC1AM = KRHO2 + NRHO2
      KC2AM = KC1AM + NRHO1
      KEND1 = KC2AM + NC2AM
      LWRK1 = LWORK - KEND1
      IF ( LWRK1 .LE. 0 ) CALL QUIT('Insufficient workspace in CC_ATRR')
C
C
C---------------------------------------------------------------------
C     Prepare the C-amplitudes.
C---------------------------------------------------------------------
C
      IF ( DEBUG .OR. LOCDBG ) THEN
         RHO1N = DDOT(NT1AM(ISYMTR),WORK(KRHO1),1,WORK(KRHO1),1)
         WRITE(LUPRI,1) 'Norm of C1AM -first in CC_ATRR:  ',RHO1N
         IF (.NOT. CCS) THEN
            RHO2N = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
            WRITE(LUPRI,1) 'Norm of C2AM -first in CC_ATRR:  ',RHO2N
         ENDIF
      ENDIF
C
      IF (.NOT. CCS) NRHO2 = NT2AM(ISYMTR)
C
      DO 70 IV = 1, NSIMTR
         NR1 = IV + K1 - 1
         IF (.NOT. (CCS.OR.CC2)) THEN
            NR2 = IV
         ELSE
            NR2 = NR1
         ENDIF
         CALL CC_WVEC(LUFC1,FC1AM,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                NR1,WORK(KRHO1))
         IF (.NOT.CCS) THEN
            CALL CC_WVEC(LUFC2,FC2AM,NRHO2,NRHO2,NR2,WORK(KRHO2))
         ENDIF
         IF (LCONT) THEN
           CONT(1) = DDOT(NT1AM(ISYMTR),WORK(KRHO1),1,WORK(KRHO1),1)
           IF (CCS) THEN
             CONT(2) = 0.0D0
           ELSE IF (CCR12.AND.IANR12.EQ.2) THEN
             CONT(2) = 0.0D0
           ELSE
             CONT(2) = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KRHO2),1)
           END IF
           CONT(3) = 0.0D0
         END IF
         IF (CCR12) THEN
            KRHO12 = NT1AM(ISYMTR)+NRHO2+1
            CALL CC_WVEC(LUFC12,FC12AM,NTR12AM(ISYMTR),NTR12AM(ISYMTR),
     *                   NR1,WORK(KRHO12))
         END IF
  70  CONTINUE
C
      CALL DCOPY(NT1AM(ISYMTR),WORK(KRHO1),1,WORK(KC1AM),1)
      IF ( .NOT. CCS ) THEN
CCH      IF ( ISIDE .GE. 1) THEN
         IF ( MYSIDE .GE. 1) THEN
            CALL CCLR_DIASCL(WORK(KRHO2),TWO,ISYMTR)
         ENDIF
         CALL CC_T2SQ(WORK(KRHO2),WORK(KC2AM),ISYMTR)
      ENDIF
C
C---------------
C     File open.
C---------------
C
      IF (.NOT. (CCS.OR.CC2)) THEN
         CALL WOPEN2(LUFSD,FR2SD,64,0)
      ENDIF
C
C---------------------
C     Zero rho vector.
C---------------------
C
      NRHO2 = MAX(NT2AM(ISYMTR),2*NT2ORT(ISYMTR))          
      IF ( CC2 ) NRHO2 = NT2AM(ISYMTR)
      IF ( CCS ) NRHO2 = 2
      CALL DZERO(WORK(KRHO1),NRHO1)
      CALL DZERO(WORK(KRHO2),NRHO2)
      DO 80 IV = 1, NSIMTR
         NR1 = IV + K1 - 1
         IF (.NOT. (CCS.OR.CC2)) THEN
            NR2 = IV
         ELSE
            NR2 = NR1
         ENDIF
         CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                NR1,WORK(KRHO1))
         IF (.NOT.CCS) THEN
            CALL CC_WVEC(LUFSD,FR2SD,NRHO2,NRHO2,NR2,WORK(KRHO2))
         ENDIF
  80  CONTINUE
C
C----------------------------------
C     Calculate transformed vectors.
C-----------------------------------
C
      LRHO1 = NT1AM(ISYMTR)
C
      IF (MYSIDE .EQ. 1) THEN
c        get A*R
         CALL CC_RHTR(ECURR,FRHO1,LUFR1,FR2SD,LUFSD,FRHO12,LUFR12,
     *                FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                WORK(KRHO1),WORK(KRHO2),
     *                WORK(KC1AM),WORK(KC2AM),
     *                WORK(KEND1),LWRK1,NSIMTR,
     *                IVEC,ITR,LRHO1,LCONT,CONT(7),APROXR12)
      ELSE IF (MYSIDE .EQ. -1) THEN
c        get R*A
         CALL CC_LHTR(ECURR,
     *                FRHO1,LUFR1,FR2SD,LUFSD,FRHO12,LUFR12,
     *                FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *                WORK(KRHO1),WORK(KRHO2),
     *                WORK(KC1AM),WORK(KC2AM),
     *                WORK(KEND1),LWRK1,NSIMTR,
     *                IVEC,ITR,LRHO1,APROXR12)
      ELSE 
         CALL QUIT('CC_ATRR; ISIDE should be -1 or +1 ')
      ENDIF

      IF ( CCR12 ) THEN
         IF (MYSIDE .EQ. 1) THEN
            CALL CC_R12METRIC(ISYMTR,BRASCL,KETSCL,WORK(KEND1),LWRK1,
     *                        FC2AM,LUFC2,FC12AM,LUFC12,FS12AM,LUFS12,
     *                        FS2AM,LUFS2,IVEC,.FALSE.,DUMMY)
         ELSE IF (MYSIDE .EQ. -1) THEN
            CALL CC_R12METRIC(ISYMTR,0.5D0*KETSCL,2.0D0*BRASCL,
     *                        WORK(KEND1),LWRK1,FC2AM,LUFC2,
     *                        FC12AM,LUFC12,FS12AM,LUFS12,FS2AM,LUFS2,
     *                        IVEC,.FALSE.,DUMMY)
         END IF

         KS12AM = NT1AM(ISYMTR)+NT2AM(ISYMTR)+NTR12AM(ISYMTR)+1
         CALL CC_RVEC(LUFS12,FS12AM,NTR12AM(ISYMTR),NTR12AM(ISYMTR),
     *                IVEC,WORK(KS12AM))
         IF (IANR12.EQ.2) THEN
           KS2AM = KS12AM + NTR12AM(ISYMTR)
           KEND1 = KS2AM + NT2AM(ISYMTR)
           LWRK1 = LWORK - KEND1
           IF (LWRK1.LT.0) THEN
             CALL QUIT('Insufficient work space in cc_atrr')
           END IF
           CALL CC_RVEC(LUFS2,FS2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                  IVEC,WORK(KS2AM))
         END IF
      END IF

      NRHO2 = MAX(NT2AM(ISYMTR),2*NT2ORT(ISYMTR))
      IF (CC2 ) NRHO2 = NT2AM(ISYMTR)
      DO 90 IV = 1, NSIMTR
         NR1 = IV + K1 - 1
         IF (.NOT. (CCS.OR.CC2)) THEN
            NR2 = IV
         ELSE
            NR2 = NR1
         ENDIF
         CALL CC_RVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                NR1,WORK(KRHO1))
         IF (.NOT.CCS) THEN
            CALL CC_RVEC(LUFSD,FR2SD,NRHO2,NT2AM(ISYMTR),
     *                   NR2,WORK(KRHO2))
         END IF
         IF (CCR12) THEN
            KRHO12 = NT1AM(ISYMTR)+NT2AM(ISYMTR)+1
            CALL CC_RVEC(LUFR12,FRHO12,NTR12AM(ISYMTR),NTR12AM(ISYMTR),
     *                   NR1,WORK(KRHO12))
         END IF

         IF (LCONT) THEN
           IF (CCS) THEN
             KSCR12 = NT1AM(ISYMTR)+NT2AM(ISYMTR)+1 
             KEND2  = KSCR12 + NT1AM(ISYMTR)
           ELSE IF (CCR12) THEN
             IF (IANR12.EQ.2) THEN
               KSCR12 = KS2AM + NT2AM(ISYMTR)
             ELSE
               KSCR12 = KS12AM + NTR12AM(ISYMTR)
             END IF
             KEND2  = KSCR12 + 
     &                MAX(NT1AM(ISYMTR),NT2AM(ISYMTR),NTR12AM(ISYMTR))
           ELSE
             KSCR12 = NT1AM(ISYMTR)+NT2AM(ISYMTR)+1 
             KEND2  = KSCR12 + MAX(NT1AM(ISYMTR),NT2AM(ISYMTR))
           END IF
           LWRK2  = LWORK - KEND2
           IF(LWRK2.LE.0)CALL QUIT('Insufficient workspace in CC_ATRR')

           CALL CC_RVEC(LUFC1,FC1AM,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                  1,WORK(KSCR12))
           CONT(4) = DDOT(NT1AM(ISYMTR),WORK(KRHO1),1,WORK(KSCR12),1)

           IF (CCS) THEN
             CONT(5) = 0.0d0
           ELSE
             CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     &                    1,WORK(KSCR12))
             CONT(5) = DDOT(NT2AM(ISYMTR),WORK(KRHO2),1,WORK(KSCR12),1)
chf
             IF (CCR12.AND.IANR12.EQ.2) THEN
               CONT(2)= CONT(2)+ 
     &                  DDOT(NT2AM(ISYMTR),WORK(KS2AM),1,WORK(KSCR12),1)
             END IF
chf
           END IF

           IF (CCR12) THEN
             CALL CC_RVEC(LUFC12,FC12AM,NTR12AM(ISYMTR),NTR12AM(ISYMTR),
     *                     NR1,WORK(KSCR12))
             CONT(6)=DDOT(NTR12AM(ISYMTR),WORK(KSCR12),1,WORK(KRHO12),1)
             CONT(3)=DDOT(NTR12AM(ISYMTR),WORK(KSCR12),1,WORK(KS12AM),1)
           ELSE
             CONT(6) = 0.0D0
             CONT(3) = 0.0D0
           END IF

         END IF

         IF ((IPRINT.GT.45) .OR. LOCDBG) THEN
            CALL AROUND('CC_ATRR: RHO = trans. Vector ')
            CALL CC_PRP(WORK(KRHO1),WORK(KRHO2),ISYMTR,1,1)
            IF (CCR12) THEN
              CALL CC_PRPR12(WORK(KRHO12),ISYMTR,1,.TRUE.)
            END IF
         ENDIF
         IF (.NOT.(CCS.OR.CC2)) THEN
            CALL CC_WVEC(LUFR2,FRHO2,NT2AM(ISYMTR),
     *                   NT2AM(ISYMTR),NR1,WORK(KRHO2))
         ENDIF
  90  CONTINUE
C
C----------------------------
C     Close and delete files.
C----------------------------
C
      CALL CC_FILCL(FRHO1,LUFR1,FRHO2,LUFR2,FRHO12,LUFR12,
     *              FC1AM,LUFC1,FC2AM,LUFC2,FC12AM,LUFC12,
     *              FS12AM,LUFS12,FS2AM,LUFS2)
      IF (.NOT. (CCS.OR.CC2)) THEN
         CALL WCLOSE2(LUFSD,FR2SD,'DELETE')
      ENDIF      
C
   1  FORMAT(1x,A35,1X,E20.10)
C
      CALL QEXIT('CC_ATRR')
      RETURN
      END
c*DECK CC_REDEIG
       SUBROUTINE CC_REDEIG(WORK,LWORK,OMEHIT)
C
C-----------------------------------------------------------------------------
C
C     Purpose: Find exci with OMEHIT and skip further calculation of 
C              as many of the other as possible.
C
C     Written by Ove Christiansen 230899 
C
C-----------------------------------------------------------------------------

#include <implicit.h>
#include <priunit.h>
#include <cclr.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccsections.h> 
#include <ccexci.h> 
#include <ccexgr.h>
#include <ccfdgeo.h>
#include <ccgr.h>
#include <cclres.h>
C
      DIMENSION WORK(LWORK)
C
      CALL QENTER('CC_REDEIG')
C
      WRITE (LUPRI,'(/,1X,A)')
     *'*********************************************************'//
     *'**********'
      WRITE(LUPRI,'(/,1X,A,F20.10)') 
     *   'Search for reference excitation energy: ',OMEHIT
      WRITE(LUPRI,'(A)') ' List of excitation energies:'
      WRITE(LUPRI,'(A)') ' Nr.   Sym.  Multi Nr.in sym/mul Exci(au.)'
      DO IEXCI=1,NEXCI
         NR  = IEXCI
         IS  = ILSTSYM('RE ',IEXCI)
         IMUL = IMULTE(IEXCI)
         IF (IMUL.EQ.1) IST = IEXCI-ISYOFE(ILSTSYM('RE ',IEXCI))
         IF (IMUL.EQ.3) IST = IEXCI-ITROFE(ILSTSYM('RE ',IEXCI))
         OM  = EIGVAL(IEXCI)
         WRITE(LUPRI,'(4I5,F20.10)') NR,IS,IMUL,IST,OM
      ENDDO
C
      OMEREF = 1.0D10
      IOMREF = 1
      DO IEXCI=1,NEXCI
         IF (ABS(OMEHIT-EIGVAL(IEXCI)).LT.OMEREF) THEN
            IOMREF = IEXCI
            OMEREF = ABS(OMEHIT-EIGVAL(IEXCI))
         ENDIF 
      ENDDO
C
      ISYHIT = ILSTSYM('RE ',IOMREF)
      IMUHIT = IMULTE(IOMREF)
      IF (IMUHIT.EQ.1) ISTHIT = IOMREF - ISYOFE(ISYHIT)
      IF (IMUHIT.EQ.3) ISTHIT = IOMREF - ITROFE(ISYHIT)
C
C     We will, we will, MUH It!!!
C     Let me hear it again
C     We will, we will, MUH It!!!
C
      WRITE(LUPRI,'(/,F15.10,4(A,I3))')
     *        OMEHIT,' is closest to exci. nr.',IOMREF,
     *        ' which is nr. ',ISTHIT,
     *        ' of symmetry ',ISYHIT,
     *        ' and multiplicity ',IMUHIT
C
C     Test if closest correspondence is ok.
C
      IF (MARGIN) THEN
         WRITE(LUPRI,'(2(A,F15.6))') ' Margin: ',XMARGIN,
     *                            ' correspondence:',OMEREF
         IF (OMEREF.LT.XMARGIN) THEN
           WRITE(LUPRI,'(/,A)') ' Closest correspondence is acceptable'
         ELSE
           WRITE(LUPRI,'(/,A)') 
     *     ' Closest correspondence is NOT acceptable'
           CALL QUIT(
     *    ' Search for specific excitation energy was not satisfactory')
         ENDIF
      ENDIF
C
      DO IMULT = 1, 3, 2
        DO ISYM = 1, NSYM
          IF ((ISYM.EQ.ISYHIT).AND.(IMULT.EQ.IMUHIT)) THEN
              NCCEXCI(ISYM,IMULT) = ISTHIT
          ELSE 
              NCCEXCI(ISYM,IMULT) = 0
          ENDIF
        ENDDO
      ENDDO
C
      DO I=1,ISTHIT
        IF (IMUHIT.EQ.1) EIGVAL(I) = EIGVAL(ISYOFE(ISYHIT)+I)
        IF (IMUHIT.EQ.3) EIGVAL(I) = EIGVAL(ITROFE(ISYHIT)+I)
      ENDDO
C     ----------------------------
C     set up NEXCI + sym into et.
C     ----------------------------
      NEXCI  = 0
      NTRIP  = 0
      DO ISYM = 1,NSYM
         ISYOFE(ISYM) = NEXCI
         ITROFE(ISYM) = ISYOFE(ISYM) + NCCEXCI(ISYM,1)
         NEXCI        = ITROFE(ISYM) + NCCEXCI(ISYM,3)
         NTRIP        = NTRIP        + NCCEXCI(ISYM,3)
         DO IEX = ISYOFE(ISYM)+1, NEXCI
            ISYEXC(IEX) = ISYM
         END DO
         DO IEX = ISYOFE(ISYM)+1, ITROFE(ISYM)
            IMULTE(IEX) = 1
         END DO
         DO IEX = ITROFE(ISYM)+1, NEXCI
            IMULTE(IEX) = 3
         END DO
      END DO
C
      CALL FLSHFO(LUPRI)
C
      IF (IPRINT.GT.15) THEN
         WRITE(LUPRI,*) 'IN CC_REDEIG after Reinit'
         WRITE(LUPRI,*) 'NEXCI: ',NEXCI
         WRITE(LUPRI,*) 'Singlet: ',(NCCEXCI(J,1),J=1,NSYM)
         WRITE(LUPRI,*) 'Triplet: ',(NCCEXCI(J,3),J=1,NSYM)
         WRITE(LUPRI,*) 'ISYOFE:',(ISYOFE(J), J=1,NSYM)
         WRITE(LUPRI,*) 'ITROFE:',(ISYOFE(J), J=1,NSYM)
         WRITE(LUPRI,*) 'ISYEXC:',(ISYEXC(J), J=1,NEXCI)
         WRITE(LUPRI,*) 'IMULTE:',(IMULTE(J), J=1,NEXCI)
         WRITE(LUPRI,*) 'EIGVAL:',(EIGVAL(J), J=1,NEXCI)
      ENDIF
C
C-------------------------------
C     Overwriting IXSTAT IXSTSY.
C-------------------------------
C
      IXSTSY = ISYHIT
      IXSTAT = ISTHIT
      IXSTMU = IMUHIT
      OMECCX = EIGVAL(ISTHIT)
      ECCXST = ECCGRS + OMECCX
C
C-----------------------------------------------------
C     Calculate only residues for this specific state.
C     Input overwritten.
C-----------------------------------------------------
C
      IF (CCLRSD) THEN
         SELLRS =.TRUE.
         NSELRS = 1
         ISELRS(NSELRS,1) = ISYHIT
         ISELRS(NSELRS,2) = ISTHIT
      ENDIF
C
      WRITE(LUPRI,*) ' Note: Optimization and residue calculation '//
     *   ' only carried out for this state'
C
      WRITE (LUPRI,'(/,1X,A)')
     *'*********************************************************'//
     *'**********'
C
      CALL QEXIT('CC_REDEIG')
C
      END
*=====================================================================*
      subroutine cctstsol(solvec,eigval,thrld,nred,nvec)
*---------------------------------------------------------------------*
*     Purpose: test for degenerate eigenvalues/vectors.
*              for degenerate solutions orthogonalize eigenvectors
*              all eigenvectors are normalized to one
*
*     Written by Christof Haettig, Februar 2003, Aarhus
*---------------------------------------------------------------------*
      implicit none
#include <priunit.h>
#include <cclr.h>

      logical locdbg
      parameter (locdbg = .false.)

      integer ivec, jvec, ndeg, nvec, nred
#if defined (SYS_CRAY) 
      real fac, ddot, thrld, solvec(maxred,*), eigval(*)
#else
      double precision fac, ddot, thrld, solvec(maxred,*), eigval(*)
#endif

      if (locdbg) then
        write(lupri,*) 'Matrix with solution vectors in reduced space:'
        call output(solvec,1,maxred,1,maxred,maxred,maxred,1,lupri)
      end if

      do ivec = 1, nvec
        ndeg = 1
        do jvec = 1, nvec
          if (ivec.ne.jvec. and. 
     &          abs(eigval(ivec)-eigval(jvec)).lt.thrld) then
            ndeg = ndeg+1
            fac = -ddot(nred,solvec(1,jvec),1,solvec(1,ivec),1)
            call daxpy(nred,fac,solvec(1,jvec),1,solvec(1,ivec),1)
          end if
        end do
        if (ndeg.gt.1) then
          write(lupri,'(/1x,a,i3,f20.10)') 
     &     'degenerate eigenvector found:',ivec,eigval(ivec)
          write(lupri,'(1x,a,i2,1x,a,/)') 
     &     'eigenvector was orthogonalized against ',ndeg-1,
     &     'degenerate vectors'
        end if
        fac = dsqrt(ddot(nred,solvec(1,ivec),1,solvec(1,ivec),1))
        if (dabs(fac-1.0d0).gt.1.0d-14) then
          call dscal(nred,1.0d0/fac,solvec(1,ivec),1)
          write(lupri,'(/1x,a,i3,f16.8,a,g16.8)')
     &      'eigenvector for state',ivec,eigval(ivec),
     &      'has been renormalized... scaling factor:',fac
        end if
      end do
 
      return
      end 
*=====================================================================*
      SUBROUTINE CC_EXCIT_CONT(ISTATE,ISYM,TRIPLET,MODEL,CONT,SCONT)
      IMPLICIT NONE
#include <priunit.h>
#include <ccsdinp.h>
#include <r12int.h>

      CHARACTER*(*) MODEL
      CHARACTER*(8) MULTIPLICITY
      LOGICAL TRIPLET
      INTEGER ISTATE,ISYM
     
#if defined (SYS_CRAY)
      REAL CONT(*), SCONT(*), EAE, ESE
#else
      DOUBLE PRECISION CONT(*), SCONT(*), EAE, ESE
#endif

      MULTIPLICITY = 'singlet '
      IF (TRIPLET) MULTIPLICITY = 'triplet '

      WRITE(LUPRI,'(//A)') 'Analysis of excitation energy:'
      WRITE(LUPRI,'(30("="),/)')
      WRITE(LUPRI,'(A,I3,5X,2A,5X,A,I3)')'Symmetry:',ISYM,
     *   'Multiplicity: ',MULTIPLICITY,'State nr.:',ISTATE
      WRITE(LUPRI,'(/72("-")/)') 
      WRITE(LUPRI,'(3(A,F10.6,5X))') 'C1 * R1  :',CONT(1),
     *   'C1  * C1                :',SCONT(1),'Ratio:',CONT(1)/SCONT(1)
      IF (.NOT.CCS) THEN
        IF (CCR12.AND.(IANR12.EQ.1.OR.IANR12.EQ.3)) THEN
         WRITE(LUPRI,'(3(A,F10.6,5X))') 'C2 * R2  :',CONT(2),
     *    'C2  * C2                :',SCONT(2),'Ratio:',
     *     CONT(2)/SCONT(2)
        ELSE IF (CCR12.AND.IANR12.EQ.2) THEN
         WRITE(LUPRI,'(3(A,F10.6,5X))') 'C2 * R2  :',CONT(2),
     *    'C2*C2 + C2*S23*C12      :',SCONT(2),'Ratio:',
     *    CONT(2)/SCONT(2)
        ELSE
         WRITE(LUPRI,'(3(A,F10.6,5X))') 'C2 * R2  :',CONT(2),
     *    'C2  * C2                :',SCONT(2),'Ratio:',
     *     CONT(2)/SCONT(2)
        END IF
      END IF
      IF (CCR12) THEN
        IF (IANR12.EQ.1.OR.IANR12.EQ.3) THEN
          WRITE(LUPRI,'(3(A,F10.6,5X))') 'C12 * R12:',CONT(3),
     *    'C12 * S * C12           :',SCONT(3),'Ratio:',CONT(3)/SCONT(3)
        ELSE IF (IANR12.EQ.2) THEN
          WRITE(LUPRI,'(3(A,F10.6,5X))') 'C12 * R12:',CONT(3),
     *    'C12*S32*C2 + C12*S12*C12:',SCONT(3),'Ratio:',
     *    CONT(3)/SCONT(3)
        END IF
      END IF
      EAE = CONT(1)+CONT(2)+CONT(3) 
      ESE = SCONT(1)+SCONT(2)+SCONT(3) 
      WRITE(LUPRI,'(72("."))')
      WRITE(LUPRI,'(3(A,F10.6,5X),/)') 'C   * R  :', EAE,
     *   '  C   * C      :',ESE, 'Ratio:',EAE/ESE

      WRITE(LUPRI,'(A)') 'singles only contributions:'
      WRITE(LUPRI,'(13X,A,F10.6)') "J(T1)       :",CONT(4)/ESE
      WRITE(LUPRI,'(13X,A,F10.6)') "J(E1)       :",CONT(5)/ESE
      WRITE(LUPRI,'(A)') "...................................."
      WRITE(LUPRI,'(A,2F10.6,/)') "<R1|[Hhat,R1]|HF>        :",
     *             (CONT(4)+CONT(5))/ESE, (CONT(4)+CONT(5))/ESE

      IF (.NOT.CCS) THEN
        WRITE(LUPRI,'(A)') 'doubles contributions:'
        WRITE(LUPRI,'(13X,A,F10.6)') "G(T2)+H(T2) :",CONT(7)/ESE
        WRITE(LUPRI,'(13X,A,F10.6)') "I(T2,E1)    :",CONT(9)/ESE
        WRITE(LUPRI,'(A)') "...................................."
        WRITE(LUPRI,'(A,2F10.6)') "<R1|[[H,T2],R1]|HF>      :",
     *             (CONT(7)+CONT(9))/ESE, (CONT(7)+CONT(9))/ESE
        WRITE(LUPRI,'(13X,A,F10.6)') "G(E2)+H(E2) :",CONT(6)/ESE
        WRITE(LUPRI,'(13X,A,F10.6)') "I(E2,T1)    :",CONT(8)/ESE
        WRITE(LUPRI,'(A)') "...................................."
        WRITE(LUPRI,'(A,F10.6)') "<R1|[Hhat,R2]|HF>        :",
     *             (CONT(6)+CONT(8))/ESE

        IF (CC2) THEN
          WRITE(LUPRI,'(3X,A,F10.6)') "<R2|[Hhat,R1]|HF>     :",
     *            CONT(10)/ESE
          WRITE(LUPRI,'(3X,A,F10.6)') "<R2|[F,R2]|HF>        :",
     *            CONT(11)/ESE
          WRITE(LUPRI,'(A)') "...................................."
          WRITE(LUPRI,'(A,2F10.6)') "sum of doubles blocks    :",
     *     (CONT(6)+CONT(8)+CONT(10)+CONT(11))/ESE, 
     *     (CONT(6)+CONT(8)+CONT(10)+CONT(11))/ESE
          WRITE(LUPRI,'(A,F10.6)') "doubles total            :",
     *     (CONT(7)+CONT(9)+CONT(6)+CONT(8)+CONT(10)+CONT(11))/ESE
        END IF
      END IF

      IF (CCSD) THEN
          WRITE(LUPRI,'(A)') 'Analysis incomplete for CCSD!'
      END IF

      IF (CCR12) THEN
        WRITE(LUPRI,'(/A)') 'R12 contributions:'
        WRITE(LUPRI,'(4X,A,2F10.6)') "<R1|[[H,T2'],R1]|HF> :",
     *       CONT(18)/ESE, CONT(18)/ESE
        IF (IANR12.EQ.2.OR.IANR12.EQ.3) THEN
          WRITE(LUPRI,'(4X,A,2F10.6)') "1HP and 1IP contrib  :",
     *         CONT(22)/ESE
        END IF
        WRITE(LUPRI,'(4X,A,F10.6)') "<R1|[Hhat,R2']|HF>   :",
     *       CONT(19)/ESE
        IF (IANR12.EQ.2.OR.IANR12.EQ.3) THEN
          WRITE(LUPRI,'(4X,A,F10.6)') "2HP and 2IP contrib  :",
     *         CONT(23)/ESE
        END IF
        WRITE(LUPRI,'(4X,A,F10.6)') "<R2'|[Hhat,R1]|HF>   :",
     *       CONT(20)/ESE
        WRITE(LUPRI,'(4X,A,F10.6)') "<R2'|[F,R2']|HF>     :",
     *       CONT(21)/ESE
        IF (IANR12.EQ.2.OR.IANR12.EQ.3) THEN
          WRITE(LUPRI,'(4X,A,F10.6)') "<R2'|[F,R2]|HF>      :",
     *       CONT(12)/ESE
          WRITE(LUPRI,'(4X,A,F10.6)') "<R2 |[F,R2']|HF>     :",
     *       CONT(13)/ESE
          WRITE(LUPRI,'(A)') "...................................."
          WRITE(LUPRI,'(A,2F10.6)') "sum of R12 doubles blocks:",
     *        (CONT(19)+CONT(20)+CONT(21)+CONT(12)+CONT(13))/ESE, 
     *        (CONT(19)+CONT(20)+CONT(21)+CONT(12)+CONT(13))/ESE
          WRITE(LUPRI,'(A,F10.6)') "R12 doubles total        :",
     *     (CONT(18)+CONT(19)+CONT(20)+CONT(21)+CONT(12)+CONT(13))/ESE
          WRITE(LUPRI,'(A)') "...................................."
          WRITE(LUPRI,'(A,10X,F10.6)') "total excitation energy  :",
     *     (CONT(4)+CONT(5)+CONT(6)+CONT(7)+CONT(8)+CONT(9)+
     *      CONT(10)+CONT(11)+CONT(18)+CONT(19)+CONT(20)+
     *      CONT(21)+CONT(12)+CONT(13))/ESE
        ELSE IF (IANR12.EQ.1) THEN
          WRITE(LUPRI,'(A)') "...................................."
          WRITE(LUPRI,'(A,2F10.6)') "sum of R12 doubles blocks:",
     *        (CONT(19)+CONT(20)+CONT(21))/ESE,
     *        (CONT(19)+CONT(20)+CONT(21))/ESE
          
          WRITE(LUPRI,'(A,F10.6)') "R12 doubles total        :",
     *          (CONT(18)+CONT(19)+CONT(20)+CONT(21))/ESE
          WRITE(LUPRI,'(A)') "...................................."
          WRITE(LUPRI,'(A,10X,F10.6)') "total excitation energy  :",
     *     (CONT(4)+CONT(5)+CONT(6)+CONT(7)+CONT(8)+CONT(9)+
     *      CONT(10)+CONT(11)+CONT(18)+CONT(19)+CONT(20)+
     *      CONT(21))/ESE
        END IF
      END IF

      WRITE(LUPRI,'(72("-"))') 
      WRITE(LUPRI,'(//)')

      RETURN
      END 
*=====================================================================*
