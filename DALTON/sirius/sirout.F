!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
C
C  /* Deck sirout */
      SUBROUTINE SIROUT(ICONV,WRK,LWRK)
C
C  4-June-1991 hjaaj: driver routine for core allocation
C
#include "implicit.h"
      DIMENSION WRK(LWRK)
#include "litinfo.h"
#include "iratdef.h"
C
C Used from common blocks:
C  INFINP : MAXMAC
#include "maxorb.h"
#include "infinp.h"
C
      KIINFO = 1
      KDINFO = KIINFO + (MAXMAC+10)*((LIINFO-1)/IRAT + 1)  ! +10: saving space for up to 10 backsteps
      KWRK1  = KDINFO + (MAXMAC+10)*  LIINFO
      LWRK1  = LWRK   - KWRK1
      CALL SIROUT_1(ICONV,WRK(KIINFO),WRK(KDINFO),WRK(KWRK1),LWRK1)
      RETURN
      END
C/* Deck sirout_1 */
      SUBROUTINE SIROUT_1(ICONV,IINFO,DINFO,WRK,LWRK)
C
C Hans Agren 8-apr-1984
C Revised 4-Jun-1991 hjaaj (check DOCINO)
C
C Purpose: to control output of summary of SIRIUS calculation
C
C MOTECC-90: This module, SIROUT, is described in the input/ouput
C            Documentation of MOTECC-90.
C
#ifdef HAS_PCMSOLVER
      use pcm_config, only: pcm_configuration, pcm_cfg
#endif
      use qmcmm_fock, only: scf_qmnpmm_out
      use pelib_interface, only: use_pelib
      use qfitlib_interface, only: qfitlib_ifc_sirius_fit,
     * qfitlib_ifc_initialize, qfitlib_ifc_results
#include "implicit.h"
C
#include "litinfo.h"
      DIMENSION IINFO(LIINFO,*),DINFO(LDINFO,*),WRK(LWRK)
C
C Used from common blocks:
C   gnrinf : WFTYPE, ?
C   pgroup : REP,GROUP
C   INFINP : ISPIN,LSYM,ISTATE, NACTEL,POTNUC, THRPWF, DOCINO, MCTYPE
C   INFVAR : NCONF,NWOPT
C   INFORB : NCMOT,NORBT,NNASHX,NISHT,?
C   INFDIM : LCINDX
C
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "priunit.h"
#include "gnrinf.h"
#include "pgroup.h"
#include "infinp.h"
#include "infvar.h"
#include "inforb.h"
#include "infind.h"
#include "infdim.h"
#include "infopt.h"
#include "inftap.h"
#include "infpri.h"
! MM, PCM
#include "pcmdef.h"
#include "pcm.h"
#include "pcmlog.h"
#include "qm3.h"
#include "qmmm.h"
! MOLDEN
#include "molde.h"
C
C
C EXTERNAL FUNCTIONS
C
      LOGICAL FNDLAB
C
      LOGICAL CICALC, MOISNO, LDUMP
      PARAMETER (D0 = 0.0D0)
C
      CHARACTER   CONV_text*17
      CHARACTER*8 TABLE(4), CMO_RTNLBL(2), RTNLBL(2), LAPRPC
      DATA TABLE/'********','OLDORB  ','NEWORB  ','STARTVEC'/
C
      CALL QENTER('SIROUT')
      CICALC = (FLAG(4) .AND. .NOT.FLAG(2))
C
      LDUMP = .TRUE.
      IF (IPRI4 .GE. 3) P4FLAG(3) = .TRUE.
      IF (P4FLAG(3) .AND. FLAG(2)) THEN
C
C        statistics for MCSCF optimization :
C
C
C   write out ITINFO
C
C    itm    IINFO(itm)    DINFO(itm)
C   -----   ----------    ----------
C     1       ITMAC         EMY
C     2       ITMIC         EACTIV
C     3       NREDL         EMCSCF
C     4       NLIN          DEPRED
C     5       NCLIN         DEACT
C     6       NOLIN         RATIO
C     7                     BETA
C     8                     STPLNG
C     9                     RTRUST
C    10                     GAMMA
C    11       INDGCM        GCINRM
C    12       INDGOM        GOBNRM
C    13       IEVA1         GRDNRM
C    14       NEVAL         GCIMAX
C    15                     GOBMAX
C    16                     TIMMAC
C    17                     TIMMIC
C    18                     TIMITR
C    19                     TIMLIN
C    20
C
C  DINFO(21-30):   EVAL(1-10)    eigenvalues of red. L
C  DINFO(31-40):   EVEC(1,1-10)  CREF-coef. of eigenvectors
C  DINFO(41-50):   EVEC(2,1-10)  CI gradient-coef. of eigenvectors
C  DINFO(51-60):   EVEC(3,1-10)  orb gradient-coef. of eigenvectors
C                  (all four have NEVAL entries, starting w. /IEVAL1)
C
        CALL TITLER('--- SIRIUS OPTIMIZATION STATISTICS ---',' ',200)
C       stamp date, time, and hostname to LUPRI
        CALL TSTAMP(' ',LUPRI)
C
        REWIND LUINF
        ITM = 0
 1910   CONTINUE
           ITM = ITM + 1
           READ (LUINF,ERR=1915,END=1920)
     &          (DINFO(I,ITM),I=1,LDINFO), (IINFO(I,ITM),I=1,LIINFO)
        GO TO 1910
 1915   CONTINUE
        WRITE(LUPRI,'(//A,I5)')' I/O error reading info file record',ITM
 1920   REWIND LUINF
        ITM = ITM - 1
        IF (ITM .EQ. 0) THEN
           WRITE (LUPRI,'(/A)')
     &        ' INFO: No statistics found on SIRIUS.ITINFO file.'
           GO TO 2000
        END IF
C
        EMCSCF = DINFO(3,ITM)
        GRDNRM = DINFO(13,ITM)
C
        WRITE(LUPRI,1100) (IINFO(1,I),IINFO(2,I),DINFO(3,I),
     *       DINFO(13,I),DINFO(6,I),DINFO(8,I),I=1,ITM)
C
        WRITE(LUPRI,1200) (IINFO(1,I),
     *                    IINFO(11,I),DINFO(14,I),DINFO(11,I),
     *                    IINFO(12,I),DINFO(15,I),DINFO(12,I),
     *                    DINFO(13,I),I=1,ITM)
C
        WRITE (LUPRI,1300)
        DO 200 I = 1,ITM
           IF (IINFO(2,I).NE.0) THEN
              WRITE(LUPRI,1310) IINFO(1,I),IINFO(2,I),
     *           IINFO(5,I),IINFO(6,I),
     *           DINFO(16,I),DINFO(18,I),DINFO(17,I),DINFO(19,I),
     *           (DINFO(17,I)/IINFO(2,I))
           ELSE
              WRITE(LUPRI,1310) IINFO(1,I),IINFO(2,I),
     *           IINFO(5,I),IINFO(6,I),
     *           DINFO(16,I),DINFO(18,I),DINFO(17,I),DINFO(19,I)
           END IF
  200   CONTINUE
        WRITE(LUPRI,2010) (IINFO(1,I),DINFO(1,I),DINFO(2,I),
     *    DINFO(3,I),I=1,ITM)
        WRITE(LUPRI,2040) (IINFO(1,I),DINFO(4,I),DINFO(5,I),
     *    DINFO(6,I),I=1,ITM)
        WRITE(LUPRI,2070) (IINFO(1,I),DINFO(7,I),DINFO(10,I),
     *                    DINFO(8,I),DINFO(9,I),I=1,ITM)
        IST = IINFO(13,1)
        IEND = IST - 1 + IINFO(14,1)
        IF (IEND .LT. IST) THEN
           ! no optimization in this run; initialization: IST=0, IEND=-1
        ELSE IF (IST .LT. 1 .OR. IEND .GT. 10) THEN
           WRITE (LUPRI,'(/A,2I20/A)')
     *        ' SIROUT stat info, IST and IEND =',IST,IEND,
     *        ' IST or IEND out of bounds (1:10)'
        ELSE
           DO 300 J = IST,IEND
             WRITE(LUPRI,4100) J,(IINFO(1,I),DINFO(20+J,I),
     &           DINFO(30+J,I),DINFO(40+J,I),DINFO(50+J,I), I = 1,ITM)
  300      CONTINUE
        END IF
      END IF
 2000 CONTINUE
C
 1100 FORMAT(//'  ITER ITMIC',T18,'EMCSCF',T35,'GRDNRM',
     *       T49,'RATIO',T60,'STPLNG',/,1X,69('-'),
     *       /,(2I5,F20.12,F15.10,F10.6,F15.10))
 1200 FORMAT(//'  ITER',T9,'INDGCM',T17,'GCIMAX',T29,'GCINRM',
     *       T40,'INDGOM',T48,'GOBMAX',T60,'GOBNRM',T72,'GRDNRM',
     *       /,1X,78('-'),/,(I5,I7,2F12.6,I7,3F12.6))
C     maybe
C        (I5,I10,2F10.6,I10,3F10.6)
C     or (I5,I10,1P,2D10.3,I10,3D10.3)
C     in 1200 FORMAT.
C
 1300 FORMAT(//,'  ITER ITMIC NCLIN NOLIN',T28,'TIMMAC',
     *       T38,'TIMITR',T48,'TIMMIC',T58,'TIMLIN',
     *       T68,'TIMMIC/ITMIC',/,1X,78('-'),/)
 1310 FORMAT(I5,3I6,5F10.2)
 2010 FORMAT(//,' ITER',T15,'EMY',T35,'EACTIV',T55,'EMCSCF',
     *       //,(I5,3F20.12))
 2040 FORMAT(//,' ITER',T15,'DEPRED',T35,'DEACT',T55,'RATIO',
     *       //,(I5,3F20.12))
 2070 FORMAT(//,' ITER',T10,'BETA',T25,'GAMMA',T43,'STPLNG',
     *       T63,'RTRUST',//,(I5,F16.8,F12.8,2F20.12))
 4100 FORMAT(//,' Reduced L root no.',I3,/,' ITER',T15,'EVAL',
     *       T33,'EVEC(1)',T51,'EVEC(2)',T69,'EVEC(3)',/1X,76('-'),
     *       /,(I5,4F18.12))
C
C
C *** If FLAG(3) then call SIRCAN
C     (unless MC optimization, i.e. FLAG(2), which isn't converged,
C      SIRCAN would destroy restart information).
C
      KFRSAV= 1
      KFREE = KFRSAV
      LFREE = LWRK
      CALL MEMGET2('REAL','CMO',  KCMO  ,NCMOT, WRK,KFREE,LFREE)
      CALL MEMGET2('REAL','OCCNO',KOCCNO,NORBT, WRK,KFREE,LFREE)
      CALL DZERO(WRK(KOCCNO),NORBT)
      CALL MEMGET2('REAL','CREF', KCREF ,NCONF, WRK,KFREE,LFREE)
      IF (FLAG(3)) THEN
         IF ( FLAG(2) .AND. ICONV.EQ.0 ) THEN
            NWARN = NWARN + 1
            FLAG(3) = .FALSE.
            WRITE (LUPRI,'(//A/A,I3/)')
     *         ' *** WARNING, SIRCAN not called because MCSCF'//
     *         ' optimization not converged','     ICONV =',ICONV
         ELSE IF ( FLAG(2) .AND. FLAG(25) ) THEN
            NWARN = NWARN + 1
            FLAG(3) = .FALSE.
            WRITE (LUPRI,'(//A/A/)')
     *         ' *** WARNING, SIRCAN not called because orbitals would',
     *         '     then become different from those just written on'//
     *         ' ABACUS/RESPONSE interface.'
         END IF
         IF (FLAG(3)) THEN
            KREL1 = KFREE
            CALL SIRCAN(IPRI4,WRK(KCREF),WRK(KCMO),WRK(KOCCNO),
     *                  WRK,KFREE,LFREE)
C           CALL SIRCAN(IPRCAN,CREF,CMO,OCC,WRK,KFREE,LFREE)
         END IF
      END IF
C
C     Now read the final orbitals and occupation numbers from SIRIUS.RST
C     (if label 'NATOCC  ' not found, then we assume any occupation numbers
C      already in WRK(KOCCNO) are OK)
C
      REWIND LUIT1
      CALL MOLLB2('NEWORB  ',CMO_RTNLBL,LUIT1,LUPRI)
      CALL READT(LUIT1,NCMOT,WRK(KCMO))
      REWIND LUIT1
      IF ( FNDLAB('NATOCC  ',LUIT1) ) THEN
         CALL READT(LUIT1,NORBT,WRK(KOCCNO))
      ELSE IF (MCTYPE .LE. 0) THEN
         DO, ISYM = 1,NSYM
            JOCC1  = KOCCNO -1 + IORB(ISYM)
            NISHI  = NISH(ISYM)
            WRK(JOCC1+1:JOCC1+NISHI) = 2.0D0
            NASHI  = NASH(ISYM)
         IF (NASHI.EQ.0) CYCLE
            JOCC1  = JOCC1 + NISHI
            WRK(JOCC1+1:JOCC1+NASHI) = 1.0D0
         END DO
      END IF
C
C *** FINAL RESULTS
C
      CALL MOLCHR(ICHRG)
C     ... get total molecular charge
      CALL TITLER('--- Final results from SIRIUS ---',' ',200)
      WRITE (LUPRI,5000) ISPIN,LSYM,REP(LSYM-1),GROUP,ICHRG
      IF (NFIELD .GT. 0) THEN
         WRITE (LUPRI,'(/A/A/A)')
     *      ' The molecule is placed in a static field.',
     *      ' Field strength (a.u.)          Field operator',
     *      ' ---------------------          --------------'
csonia 04/10/95
         IF (LUCME.GT.0) WRITE (LUCME,'(/A/A/A)')
     &         ' The molecule is placed in a static field.',
     &         ' Field strength (a.u.)          Field operator',
     &         ' ---------------------          --------------'
csonia 04/10/95
         DO 100 IFIELD = 1,NFIELD
            WRITE (LUPRI,'(1P,G18.8,18X,A8)')
     *         EFIELD(IFIELD),LFIELD(IFIELD)
csonia 04/10/95
            IF (LUCME.GT.0) WRITE (LUCME,'(1P,G18.8,18X,A8)')
     &         EFIELD(IFIELD), LFIELD(IFIELD)
csonia 04/10/95
  100    CONTINUE
         WRITE (LUPRI,'()')
      END IF
      IF (FLAG(16)) THEN
         WRITE (LUPRI,'(/5X,A,2(/10X,A,F12.6)/10X,A,3F12.6/10X,A,I12)')
     *'SOLVATION MODEL: molecule is in a cavity in a dielectric medium,'
     *  ,'dielectric constant =',EPSOL,
     *   'static diel. const. =',EPSTAT,
     *   'cavity dimensions   =',RSOL,
     *   'max l value         =',LSOLMX
         IF (INERSF) WRITE (LUPRI,'(/5X,A)')
     *      'This is the final state in a calculation with'//
     *      ' inertial polarization.'
      ELSE IF (PCM .AND. (.NOT.MMPCM)) THEN
         IF (NONEQ) THEN
            WRITE (LUPRI,'(/5X,A,2(/10X,A,F12.6))')
     *           'SOLVATION MODEL: polarizable continuum model (PCM),'
     *           ,'dielectric constant =',EPSINF,
     *           'static diel. const. =',EPS
         ELSE
            WRITE (LUPRI,'(/5X,A,(/10X,A,F12.6))')
     *           'SOLVATION MODEL: polarizable continuum model (PCM),'
     *           ,'dielectric constant =',EPS
         END IF
#ifdef HAS_PCMSOLVER
      ELSE IF (pcm_cfg%do_pcm) THEN
         WRITE (LUPRI,'(/5X,A,(/10X,A,F12.6))')
     *        'SOLVATION MODEL: polarizable continuum model (PCM),'
         WRITE(LUPRI,'(5X,A24,F20.12)') 'Solvation energy:     ',ESOLT
#endif
      ELSE IF (QM3 .AND. (.NOT.MMPCM)) THEN
        WRITE(LUPRI,'(/5X,A/)') 'QM/MM "QM3" calculation converged :'
        WRITE(LUPRI,'(8X,A22,F20.12)') 'Electrostatic energy: ', ENSQM3
        WRITE(LUPRI,'(8X,A22,F20.12)') 'Polarization energy:  ', EPOQM3
        WRITE(LUPRI,'(8X,A22,F20.12)') 'van der Waals energy: ', ECLVDW
        WRITE(LUPRI,'(8X,A22,F20.12)') 'Total QM/MM energy:   ', ESOLT
      ELSE IF (QMMM) THEN
        WRITE(LUPRI,'(/5X,A/)') 'QM/MM "QMMM" calculation converged :'
        WRITE(LUPRI,'(5X,A24,F20.12)') 'Charge contribution:    ',ECHART
        WRITE(LUPRI,'(5X,A24,F20.12)') 'Dipole contribution:    ',EDIPT
        WRITE(LUPRI,'(5X,A24,F20.12)') 'Quadrupole contribution:',EQUADT
        WRITE(LUPRI,'(5X,A24,F20.12)') 'Electronic Pol. energy: ',EDELD
        WRITE(LUPRI,'(5X,A24,F20.12)') 'Nuclear pol. energy:    ',EDNUC
        WRITE(LUPRI,'(5X,A24,F20.12)') 'Multipole Pol. energy:  ',EDMULT
        WRITE(LUPRI,'(5X,A24,F20.12)') 'Total QM/MM energy:     ',ESOLT
      ELSE IF (MMPCM) THEN
         WRITE(LUPRI,'(/5X,A/)') 'QM/MM/PCM calculation converged :'
         WRITE(LUPRI,'(8X,A22,F20.12)') 'Electrostatic energy: ', ENSQM3
         WRITE(LUPRI,'(8X,A22,F20.12)') 'Polarization energy:  ', EPOQM3
         WRITE(LUPRI,'(8X,A22,F20.12)') 'van der Waals energy: ', ECLVDW
         WRITE(LUPRI,'(8X,A22,F20.12)') 'QM/MM/PCM energy:     ', ESOLT
      ELSE IF (QMNPMM) THEN
         CALL SCF_QMNPMM_OUT()
      ELSE IF (USE_PELIB()) THEN
         CONTINUE
      ELSE
        ESOLT = D0
      END IF
C
      IF (ICONV .EQ. 0) THEN
         CONV_text = ' (NOT CONVERGED!)'
      ELSE
         CONV_text = '                 '
      END IF
C
      IF (CICALC) THEN
C        CI calculation
         IF (DOCINO) THEN
            WRITE (LUPRI,5010)
         ELSE
            WRITE (LUPRI,5012)
         END IF
      ELSE IF (FLAG(2) .AND. MCTYPE .GT. 0) THEN
C        MCSCF calculation
         IF (DOMCSRDFT) THEN
            WRITE (LUPRI,5024) ISTATE,EMCSCF,CONV_text,
     &         POTNUC,(EMCSCF-POTNUC-ESOLT)
         ELSE
            WRITE (LUPRI,5020) ISTATE,EMCSCF,CONV_text,
     &         POTNUC,(EMCSCF-POTNUC-ESOLT)
         END IF
         IF (FLAG(16) .OR. PCM) WRITE (LUPRI,5025) ESOLT
         IF (USE_PELIB()) WRITE(LUPRI,5026) ESOLT
         WRITE (LUPRI,5030) GRDNRM
      ELSE IF (DOSCF) THEN
         IF (DOHFSRDFT) THEN
           WRITE (LUPRI,5023) EMCSCF,POTNUC,(EMCSCF-POTNUC-ESOLT)
         ELSE IF (.NOT.DODFT) THEN
C        ... DIIS-HF or QC-HF optimization
             WRITE (LUPRI,5021) EMCSCF,CONV_text,
     &          POTNUC,(EMCSCF-POTNUC-ESOLT)
         ELSE
C        ... DIIS-DFT or QC-DFT optimization
             WRITE (LUPRI,5022) EMCSCF,CONV_text,
     &          POTNUC,(EMCSCF-POTNUC-ESOLT)
         END IF
         IF (DOCCSD .AND. .NOT.(QM3)) LDUMP = .FALSE.
         IF (DOCCSD .AND. .NOT.(QMMM)) LDUMP = .FALSE.
         IF (LDUMP) THEN
            LAPRPC = 'ENERGY  '
            CALL WRIPRO(EMCSCF,"SCF/DFT   ",0,
     *                  LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *                  0.0D0,0.0D0,0.0D0,LSYM,0,ISPIN,0)
         END IF
         IF (FLAG(16)) WRITE (LUPRI,5025) ESOLT
         IF (USE_PELIB()) WRITE(LUPRI,5026) ESOLT
         WRITE (LUPRI,5030) GRDNRM
      END IF

      IF (QMMM .AND. MMPROP) THEN
        CALL MM_PROPS(WRK(KFREE),LFREE,IPQMMM)
      ENDIF

 5000 FORMAT(/'@    Spin multiplicity:',T31,I5,
     *       /'@    Spatial symmetry:',T31,I5,
     &        ' ( irrep  ',A,' in ',A,' )',
     *       /'@    Total charge of molecule:',T31,I5)
 5010 FORMAT(/'@    This was a CI-NO calculation.')
 5012 FORMAT(/'@    This was a CI calculation.')
 5020 FORMAT( '@    State number:',T31,I5,
     *      //'@    Final MCSCF energy:',T31,F20.12,A,
     *       /'@    Nuclear repulsion:',T31,F20.12,
     *       /'@    Electronic energy:',T31,F20.12)
 5021 FORMAT(/'@    Final HF energy:',T31,F20.12,A,
     *       /'@    Nuclear repulsion:',T31,F20.12,
     *       /'@    Electronic energy:',T31,F20.12)
 5022 FORMAT(/'@    Final DFT energy:',T31,F20.12,A,
     *       /'@    Nuclear repulsion:',T31,F20.12,
     *       /'@    Electronic energy:',T31,F20.12)
 5023 FORMAT(/'@    Final HF-SRDFT energy:',T31,F20.12,
     *       /'@    Nuclear repulsion:',T31,F20.12,
     *       /'@    Electronic energy:',T31,F20.12)
 5024 FORMAT( '@    State number:',T31,I5,
     *      //'@    Final MC-SRDFT energy:',T31,F20.12,A,
     *       /'@    Nuclear repulsion:',T31,F20.12,
     *       /'@    Electronic energy:',T31,F20.12)
 5025 FORMAT( '@    Solvation  energy:',T31,F20.12)
 5026 FORMAT( '@    Embedding energy:',T31,F20.12)
 5030 FORMAT(/'@    Final gradient norm:',T31,F20.12)
C
C       stamp date, time, and hostname to LUPRI
        CALL TSTAMP(' ',LUPRI)
C
C     Natural orbital analysis of MCSCF, also when orbitals are not natural orbitals
C
      IF (MCTYPE .LE. 0) THEN
         MOISNO = .TRUE.
      ELSE IF (.NOT.(FLAG(3) .AND. MCTYPE.EQ.1) .AND.
     &         .NOT.(CICALC  .AND. DOCINO) ) THEN
C     if MC/CI and NO's not already printed
         KREL1 = KFREE
         CALL MEMGET2('REAL','OCC  ',KOCC,  NORBT,WRK,KFREE,LFREE) ! we do not want to overwrite WRK(KOCCNO)
         CALL MEMGET2('REAL','UNO',  KUNO,  N2ASHX,WRK,KFREE,LFREE)
         CALL MEMGET2('REAL','CINDX',KCINDX,LCINDX,WRK,KFREE,LFREE)
         IPRNO = MAX(1,IPRI4-5)
C        980903-hjaaj: IPRI4 = IPRUSR + 5 as default in sirinp.F
         IF (NASHT .GT. 1) THEN
            CALL GETCIX(WRK(KCINDX),LSYM,LSYM,WRK(KFREE),LFREE,0)
            CALL RDCREF(WRK(KCREF))
         END IF
         CALL GETNO(WRK(KCREF),LSYM,WRK(KOCC),WRK(KUNO),WRK(KCMO),
     *              .FALSE.,.FALSE.,.TRUE.,MOISNO,WRK(KCINDX),
     *              WRK,KFREE,LFREE,LUPRI,IPRNO)
C        CALL GETNO(CVEC,ICSYM,OCCNO,UNO,CMO,ORDER,NATORB,NOAVER,
C    *              MOISNO,INDXCI,WRK,KFRSAV,LFRSAV,LUPRI,IPRINT)
         IF (.NOT. MOISNO) WRITE (LUPRI,'(3(/A))')
     *     '  --- NOTE THAT THE OCCUPATION NUMBERS ARE FOR ANALYSIS'//
     &     ' PURPOSE ONLY',
     *     '  --- The orbital coefficients have NOT been transformed',
     *     '  --- to ordered natural orbitals.'
         CALL MEMREL('SIROUT.GETNO ',WRK,KFRSAV,KREL1,KFREE,LFREE)
      END IF
C
      IF (DOMEP) THEN
         CALL CALC_MEP(WRK(KCMO),WRK,KFREE,LFREE)
      ENDIF

      IF (QFIT) THEN
         CALL QFITLIB_IFC_INITIALIZE
         CALL QFITLIB_IFC_SIRIUS_FIT(WRK(KCMO), WRK(KFREE), LFREE)
         CALL QFITLIB_IFC_RESULTS
      ENDIF

      IF (MOLDEN) THEN
         KREL2 = KFREE
         CALL MEMGET('REAL', KAOSO, N2BASX,      WRK, KFREE, LFREE)
         CALL MEMGET('REAL', KUCMO, NORBT*NBAST, WRK, KFREE, LFREE)
         CALL MEMGET('REAL', KOVEC, NBAST,       WRK, KFREE, LFREE)
         CALL MOLDEN_MOS(1,WRK(KCMO), WRK(KOCCNO), WRK(KAOSO),
     &                   WRK(KUCMO), WRK(KOVEC))
         CALL MEMREL('SIROUT.MOLDEN_MOS',WRK,KREL2,KREL2,KFREE,LFREE)
      END IF ! (MOLDEN)
C
C ***  Print final MCSCF orbitals ***
C
      IF (IPRI4 .GE. 5) P4FLAG(5) = .TRUE.
      IF (P4FLAG(5)) THEN
         WRITE(LUPRI,'(/A,2X,A,2X,A)') 'File label for MO orbitals:',
     &      CMO_RTNLBL
         IF ((IPRI4 .GT. 10) .OR. CMOPRI) THEN
            CALL PRORB(WRK(KCMO), .FALSE., LUPRI)
         ELSE
            CALL PRORB(WRK(KCMO), .TRUE.,  LUPRI)
         END IF
C        CALL PRORB(CMO, PROCC, IOUT [,THRPRI_in])
      END IF
C
C *** Print out of wavefunction in MCSCF orbital basis ***
C     (first work space allocation)
C
C
      IF (IPRI4 .GE. 4) P4FLAG(4) = .TRUE.
c      IF (P4FLAG(4) .AND. NCONF.GT.1 .AND. .NOT.DOCINO) THEN
c--renzo modif
      IF (P4FLAG(4) .AND.((NCONF.eq.1.and.DONEVPT).or.(NCONF.gt.1)).AND.
     $     .NOT.(CICALC .AND. DOCINO)) THEN
         IF (LCINDX .GT. LFREE) THEN
            WRITE (LUPRI,1077) LCINDX-LFREE
            GO TO 110
         END IF
         CALL MEMGET('REAL',KCINDX,LCINDX,WRK,KFREE,LFREE)
         CALL MEMGET('WORK',KWRK2,LWRK2,WRK,KFREE,LFREE)
C
         REWIND LUIT1
         CALL MOLLB2(TABLE(4),RTNLBL,LUIT1,LUPRI)
         IF (RTNLBL(1) .EQ. TABLE(1)) THEN
C           this is an old LUIT1 file (before 5-Aug-1986)
            I = ISTATE
         ELSE
            READ (RTNLBL(1),'(2I4)') NRS,I
c--renzo modif
            if(nconf.eq.1) i=1   !just a patch (hopefully harmless)
            IF (RTNLBL(2) .EQ. 'CISAVE ') THEN
               IF (.NOT. CICALC) THEN
                  NWARN = NWARN + 1
                  WRITE (LUPRI,'(//A/A/)')
     &               ' WARNING, internal error in SIROUT subroutine',
     &               ' CICALC false, but CISAVE label on SIRIUS.RST'
               END IF
               CICALC = .TRUE.
               IF (I .GT. 0 .AND. I .NE. ISTACI) THEN
                  NWARN = NWARN + 1
                  WRITE (LUPRI,'(//,(A,I5/))')
     *         ' WARNING from SIROUT: ISTACI specied in input:',ISTACI,
     *         '                      ISTACI read from LUIT1 :',I
               END IF
            ELSE
               IF (CICALC) THEN
                  NWARN = NWARN + 1
                  WRITE (LUPRI,'(//A/A/)')
     &               ' WARNING, internal error in SIROUT subroutine',
     &               ' CICALC true, but no CISAVE label on SIRIUS.RST'
                  CICALC = .FALSE.
               END IF
               IF (ABS(I) .NE. ISTATE) THEN
                  NWARN = NWARN + 1
                  WRITE (LUPRI,'(//,(A,I5/))')
     *         ' WARNING from SIROUT: ISTATE specied in input:',ISTATE,
     *         '                      ISTATE read from LUIT1 :',ABS(I)
               END IF
            END IF
         END IF
         IF (I .GT. 0) THEN
C            ... CI with ISTACI .gt. 0 or MC calculation
            DO 102 I = 1,(I-1)
               READ (LUIT1)
 102        CONTINUE
            ILOW = I
            IHGH = I
         ELSE IF (CICALC) THEN
            ILOW = 1
            IHGH = MIN(NRS,NROOCI)
         ELSE
C        else MCSCF with only reference vector saved (NRSAVE)
            ILOW = ABS(I)
            IHGH = ABS(I)
         END IF
C
         CALL GETCIX(WRK(KCINDX),LSYM,LSYM,WRK(KWRK2),LWRK2,0)
         DO 105 I = ILOW,IHGH
            WRITE (LUPRI,1050) THRPWF,I
            IF (I .EQ. ISTATE .AND. .NOT. CICALC) WRITE (LUPRI,1051)
            IF (I .EQ. ISTACI .AND.       CICALC) WRITE (LUPRI,1052)
            CALL READT(LUIT1,NCONF,WRK(KCREF))
            CALL PRWF(WRK(KCREF),WRK(KCINDX),LUPRI,THRPWF)
C           CALL PRWF(C,INDXCI,IOUT,THRPWF)
  105    CONTINUE
         CALL MEMREL('SIROUT.PRWF',WRK,KFRSAV,KCINDX,KFREE,LFREE)
      END IF
  110 CONTINUE
 1050 FORMAT(/' Printout of CI-coefficients abs greater than',F8.5,
     *         ' for root',I3)
 1051 FORMAT(' *** NOTE: this root is the reference state ***')
 1052 FORMAT(' *** NOTE: this root is the CI reference state ***')
 1077 FORMAT(/' *** SIROUT, Wave function not printed as requested',
     *  /,15X,'Reason: insufficient work memory space',
     *  /,15X,'need more than',I12,' additional real*8 work space')
C
C If FLAG(18):
C punch final orbitals (and, if FLAG(3): call sircan, also
C occupation numbers) formatted on unit 7
C
      IF (FLAG(18)) THEN
         REWIND LUIT1
         IF ( FNDLAB(TABLE(3),LUIT1) ) THEN
            WRITE (LUPRI,'(/A/)') ' "NEWORB " orbitals punched.'
         ELSE
            NWARN = NWARN + 1
            WRITE (LUPRI,'(//A/)')
     *         ' WARNING, new orbitals not found on SIRIUS.RST,'//
     *         ' will try "OLDORB "'
            REWIND LUIT1
            CALL MOLLAB(TABLE(2),LUIT1,lupri)
            WRITE (LUPRI,'(/A/)') ' "OLDORB " orbitals punched.'
         END IF
         CALL READT(LUIT1,NCMOT,WRK(KCMO))
         IF (FLAG(3)) THEN
            IPCTL = 1
         ELSE
            IPCTL = 0
            KOCCNO  = KCMO
         END IF
         CALL PUNMO(IPCTL,WRK(KCMO),WRK(KOCCNO))
C        CALL PUNMO(IPCTL,CMO,OCC)
      END IF
C
C end of SIROUT
C
      CALL MEMREL('SIROUT',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('SIROUT')
      RETURN
      END
C --- end of sirius/sirout.F ---
