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
C  /* Deck sirout */
      SUBROUTINE SIROUT(ICONV,WRK,LWRK)
C
C  4-June-1991 hjaaj: driver routine for core allocation
C
#include <implicit.h>
      DIMENSION WRK(LWRK)
#include <litinfo.h>
#include <iratdef.h>
C
C Used from common blocks:
C  INFINP : MAXMAC
#include <maxorb.h>
#include <infinp.h>
C
      KIINFO = 1
      KDINFO = 1 + MAXMAC*((LIINFO-1)/IRAT + 1)
      CALL SIROU1(ICONV,WRK(KIINFO),WRK(KDINFO),WRK,LWRK)
      RETURN
      END
C/* Deck sirou1 */
      SUBROUTINE SIROU1(ICONV,IINFO,DINFO,WRK,LWRK)
C
C Hans Agren 8-apr-1984
C Last revision 4-Jun-1991 hjaaj (check DOCINO)
C
C Purpose: to control print output of summary of calculation
C
C IINFO,DINFO and WRK are assumed to overlap in core so LWRK is
C valid for both.
C
C To change dimension of IINFO and/or DINFO change in /ITINFO/,
C change here, in SETUP, and in SIRCTL's call of SIROUT. NEW PROCEDURE:
C This is done automatically by changing comdeck LITINFO.
C
C MOTECC-90: This module, SIROUT, is described in the input/ouput
C            Documentation of MOTECC-90.
C
#include <implicit.h>
C
#include <litinfo.h>
      DIMENSION IINFO(LIINFO,*),DINFO(LDINFO,*),WRK(LWRK)
#include <iratdef.h>
      PARAMETER (D0 = 0.0D0)
C
C Used from common blocks:
C   INFINP : ISPIN,LSYM,ISTATE, NACTEL,POTNUC, THRPWF, DOCINO, MCTYPE
C   INFVAR : NCONF,NWOPT
C   INFORB : NCMOT,NORBT,NNASHX,NISHT,?
C   INFDIM : LCINDX
C
#include <maxash.h>
#include <maxorb.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
#include <priunit.h>
#include <molde.h>
#include <infinp.h>
#include <infvar.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infopt.h>
#include <inftap.h>
#include <infpri.h>
#include <gnrinf.h>
C----------------
C CBN+JK 03.01.06
C----------------
#include <qm3.h>
C----------------
C CBN+JK 03.01.06
C----------------
C
C
C EXTERNAL FUNCTIONS
C
      LOGICAL FNDLAB
C
      LOGICAL CICALC, MOISNO, LDUMP
C
      CHARACTER*8 TABLE(4), RTNLBL(2), LAPRPC
      DATA TABLE/'********','OLDORB  ','NEWORB  ','STARTVEC'/
C
      CALL QENTER('SIROUT')
      CICALC = (FLAG(4) .AND. .NOT.FLAG(2))
C
C Workspace check for some sections:
C  ITMAC*LBINFO for iteration information
C  NCMOT for PRORB
C
      LBIMAC = LBINFO*ITMAC
      LWRK1  = MAX(LBIMAC,NCMOT)
      IF (LWRK.LT.LWRK1) CALL ERRWRK('SIROUT',LWRK1,LWRK)
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
C                  (all four have NEVAL entries, starting w. IEVAL1)
C
        WRITE(LUW4,1000)
C       stamp date, time, and hostname to LUW4
        CALL TSTAMP(' ',LUW4)
C
        REWIND LUINF
        ITM = 0
 1910   CONTINUE
           ITM = ITM + 1
           READ (LUINF,ERR=1915,END=1920)
     *          (DINFO(I,ITM),I=1,LDINFO), (IINFO(I,ITM),I=1,LIINFO)
        GO TO 1910
 1915   CONTINUE
        WRITE(LUW4,'(//A,I5)') ' I/O error reading info file record',ITM
 1920   REWIND LUINF
        ITM = ITM - 1
        IF (ITM .EQ. 0) THEN
           WRITE (LUW4,'(//A/)') ' No statistics found on info file.'
           GO TO 2000
        END IF
C
        EMCSCF = DINFO(3,ITM)
        GRDNRM = DINFO(13,ITM)
C
        WRITE(LUW4,1100) (IINFO(1,I),IINFO(2,I),DINFO(3,I),DINFO(13,I),
     *                    DINFO(6,I),DINFO(8,I),I=1,ITM)
C
        WRITE(LUW4,1200) (IINFO(1,I),
     *                    IINFO(11,I),DINFO(14,I),DINFO(11,I),
     *                    IINFO(12,I),DINFO(15,I),DINFO(12,I),
     *                    DINFO(13,I),I=1,ITM)
C
        WRITE (LUW4,1300)
        DO 200 I = 1,ITM
           IF (IINFO(2,I).NE.0) THEN
              WRITE(LUW4,1310) IINFO(1,I),IINFO(2,I),
     *           IINFO(5,I),IINFO(6,I),
     *           DINFO(16,I),DINFO(18,I),DINFO(17,I),DINFO(19,I),
     *           (DINFO(17,I)/IINFO(2,I))
           ELSE
              WRITE(LUW4,1310) IINFO(1,I),IINFO(2,I),
     *           IINFO(5,I),IINFO(6,I),
     *           DINFO(16,I),DINFO(18,I),DINFO(17,I),DINFO(19,I)
           END IF
  200   CONTINUE
        WRITE(LUW4,2010) (IINFO(1,I),DINFO(1,I),DINFO(2,I),
     *    DINFO(3,I),I=1,ITM)
        WRITE(LUW4,2040) (IINFO(1,I),DINFO(4,I),DINFO(5,I),
     *    DINFO(6,I),I=1,ITM)
        WRITE(LUW4,2070) (IINFO(1,I),DINFO(7,I),DINFO(10,I),
     *                    DINFO(8,I),DINFO(9,I),I=1,ITM)
        IST = IINFO(13,1)
        IEND = IST - 1 + IINFO(14,1)
        IF (IST .LT. 1 .OR. IEND .GT. 10) THEN
           WRITE (LUERR,'(/A,2I20/2A)')
     *        ' SIROUT stat info, IST and IEND =',IST,IEND,
     *        ' IST or IEND out of bounds',
     *        ' - probably no optimization in this run.'
        ELSE
           DO 300 J = IST,IEND
             WRITE(LUW4,4100) J,(IINFO(1,I),DINFO(20+J,I),DINFO(30+J,I),
     *                        DINFO(40+J,I),DINFO(50+J,I), I = 1,ITM)
  300      CONTINUE
        END IF
      END IF
 2000 CONTINUE
C
 1000 FORMAT(//T16,'>>> SIRIUS OPTIMIZATION STATISTICS <<<')
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
      CALL MEMGET('REAL', KCMO,NCMOT, WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KCREF,NCONF, WRK,KFREE,LFREE)
      CALL MEMGET('REAL', KOCC,NORBT, WRK,KFREE,LFREE)
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
            CALL SIRCAN(IPRI4,WRK(KCREF),WRK(KCMO),WRK(KOCC),
     *                  WRK,KFREE,LFREE)
C           CALL SIRCAN(IPRCAN,CREF,CMO,OCC,WRK,KFREE,LFREE)
         END IF
      END IF
C
C *** FINAL RESULTS
C
      CALL MOLCHR(ICHRG)
C     ... get total molecular charge
      WRITE (LUW4,5000) ISPIN,LSYM,ICHRG
      IF (NFIELD .GT. 0) THEN
         WRITE (LUW4,'(/A/A/A)')
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
            WRITE (LUW4,'(1P,G18.8,18X,A8)')
     *         EFIELD(IFIELD),LFIELD(IFIELD)
csonia 04/10/95
            IF (LUCME.GT.0) WRITE (LUCME,'(1P,G18.8,18X,A8)')
     &         EFIELD(IFIELD), LFIELD(IFIELD)
csonia 04/10/95                                                            
  100    CONTINUE
         WRITE (LUW4,'()')
      END IF
      IF (FLAG(16)) THEN
         WRITE (LUW4,'(/5X,A,2(/10X,A,F12.6)/10X,A,3F12.6/10X,A,I12)')
     *'SOLVATION MODEL: molecule is in a cavity in a dielectric medium,'
     *  ,'dielectric constant =',EPSOL,
     *   'static diel. const. =',EPSTAT,
     *   'cavity dimensions   =',RSOL,
     *   'max l value         =',LSOLMX
         IF (INERSF) WRITE (LUW4,'(/5X,A)')
     *      'This is the final state in a calculation with'//
     *      ' inertial polarization.'
      ELSE IF (PCM) THEN
         IF (NONEQ) THEN
            WRITE (LUW4,'(/5X,A,2(/10X,A,F12.6))')
     *           'SOLVATION MODEL: polarizable continuum model (PCM),'
     *           ,'dielectric constant =',EPSINF,
     *           'static diel. const. =',EPS
         ELSE
            WRITE (LUW4,'(/5X,A,(/10X,A,F12.6))')
     *           'SOLVATION MODEL: polarizable continuum model (PCM),'
     *           ,'dielectric constant =',EPS
         END IF
C------------------
C CBN+JK 03.01.06
C------------------
      ELSE IF (QM3) THEN
         WRITE(LUW4,'(/5X,A/)') 'QM/MM calculation converged :'
         WRITE(LUW4,'(8X,A22,F20.12)') 'Electrostatic energy: ', ENSQM3
         WRITE(LUW4,'(8X,A22,F20.12)') 'Polarization energy:  ', EPOQM3
         WRITE(LUW4,'(8X,A22,F20.12)') 'van der Waals energy: ', ECLVDW
         WRITE(LUW4,'(8X,A22,F20.12)') 'Total QM/MM energy:   ', ESOLT
      ELSE
C------------------
C CBN+JK 03.01.06
C------------------
         ESOLT = D0
      END IF
      IF (CICALC) THEN
C        CI calculation
         IF (DOCINO) THEN
            WRITE (LUW4,5010)
         ELSE
            WRITE (LUW4,5012)
         END IF
      ELSE IF (FLAG(2) .AND. MCTYPE .GT. 0) THEN
C        MCSCF calculation
         WRITE (LUW4,5020) ISTATE,EMCSCF,POTNUC,(EMCSCF-POTNUC-ESOLT)
         IF (FLAG(16) .OR. PCM) WRITE (LUW4,5025) ESOLT
         WRITE (LUW4,5030) GRDNRM
      ELSE IF (DORHF) THEN
         IF (.NOT.DODFT) THEN
C        ... DIIS-HF or QC-HF optimization
             WRITE (LUW4,5021) EMCSCF,POTNUC,(EMCSCF-POTNUC-ESOLT)
         ELSE
C        ... DIIS-DFT or QC-DFT optimization
             WRITE (LUW4,5022) EMCSCF,POTNUC,(EMCSCF-POTNUC-ESOLT)
         END IF
         IF (DOCCSD .AND. .NOT.(QM3)) LDUMP = .FALSE.
         IF (LDUMP) THEN
            LAPRPC = 'ENERGY  '
            CALL WRIPRO(EMCSCF,"SCF/DFT   ",0,
     *                  LAPRPC,LAPRPC,LAPRPC,LAPRPC,
     *                  0.0D0,0.0D0,0.0D0,LSYM,0,ISPIN,0)
         END IF
         IF (FLAG(16)) WRITE (LUW4,5025) ESOLT
         WRITE (LUW4,5030) GRDNRM
      END IF

 5000 FORMAT(//T21,'>>> FINAL RESULTS FROM SIRIUS <<<',
     *       //T6,'Spin multiplicity:',T31,I5,
     *        /T6,'Spatial symmetry:',T31,I5,
     *        /T6,'Total charge of molecule:',T31,I5)
 5010 FORMAT(/T6,'This was a CI-NO calculation.')
 5012 FORMAT(/T6,'This was a CI calculation.')
 5020 FORMAT( T6,'State number:',T31,I5,
     *      //T6,'Final MCSCF energy:',T31,F20.12,
     *       /T6,'Nuclear repulsion:',T31,F20.12,
     *       /T6,'Electronic energy:',T31,F20.12)
 5021 FORMAT(/T6,'Final HF energy:',T31,F20.12,
     *       /T6,'Nuclear repulsion:',T31,F20.12,
     *       /T6,'Electronic energy:',T31,F20.12)
 5022 FORMAT(/T6,'Final DFT energy:',T31,F20.12,
     *       /T6,'Nuclear repulsion:',T31,F20.12,
     *       /T6,'Electronic energy:',T31,F20.12)
 5025 FORMAT( T6,'Solvation  energy:',T31,F20.12)
 5030 FORMAT(/T6,'Final gradient norm:',T31,F20.12)
C
C       stamp date, time, and hostname to LUW4
        CALL TSTAMP(' ',LUW4)
C
      IF (MCTYPE .GT. 0 .AND. .NOT.(FLAG(3) .AND. MCTYPE.EQ.1) .AND.
     *    .NOT.(CICALC.AND.DOCINO) ) THEN
C     if MC/CI and NO's not already printed
         KREL1 = KFREE
         CALL MEMGET('REAL',KUNO,N2ASHX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KCINDX,LCINDX,WRK,KFREE,LFREE)
         IPRNO = MAX(1,IPRI4-5)
C        980903-hjaaj: IPRI4 = IPRUSR + 5 as default in sirinp.F
         CALL GETCIX(WRK(KCINDX),LSYM,LSYM,WRK(KFREE),LFREE,0)
         CALL RDCREF(WRK(KCREF))
         CALL GETNO(WRK(KCREF),LSYM,WRK(KOCC),WRK(KUNO),WRK(KCMO),
     *              .FALSE.,.FALSE.,.TRUE.,MOISNO,WRK(KCINDX),
     *              WRK,KFREE,LFREE,LUW4,IPRNO)
C        CALL GETNO(CVEC,ICSYM,OCCNO,UNO,CMO,ORDER,NATORB,NOAVER,
C    *              MOISNO,INDXCI,WRK,KFRSAV,LFRSAV,LUPRI,IPRINT)
         IF (.NOT. MOISNO) WRITE (LUW4,'(3(/A))')
     *     '  >>> NOTE THAT THE OCCUPATION NUMBERS ARE FOR ANALYSIS'//
     &     ' PURPOSE ONLY',
     *     '  >>> The orbital coefficients have NOT been transformed',
     *     '  >>> to ordered natural orbitals.'
         CALL MEMREL('SIROUT.GETNO ',WRK,KFRSAV,KREL1,KFREE,LFREE)
      END IF
C
      IF (MOLDEN) THEN
         REWIND LUIT1
         IF ( .NOT. FNDLAB(TABLE(3),LUIT1) ) THEN
            NWARN = NWARN + 1
            WRITE (LUPRI,'(//A/)')
     *         ' WARNING, new orbitals not found on SIRIUS.RST,'//
     *         ' will try "OLDORB "'
            REWIND LUIT1
            CALL MOLLAB(TABLE(2),LUIT1,LUW4)
         END IF
         CALL READT(LUIT1,NCMOT,WRK(KCMO))
       IF (CICALC .AND. DOCINO) THEN
          CALL DZERO(WRK(KOCC),NORBT)
C         ... we do not know the orbital occupation for CINO
C             any more and cannot calculate them again 
C             (CREF was not saved on file)
       ELSE
         KREL1 = KFREE
         CALL MEMGET('REAL',KUNO,N2ASHX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KCINDX,LCINDX,WRK,KFREE,LFREE)
         IPRNO = -1
         CALL GETCIX(WRK(KCINDX),LSYM,LSYM,WRK(KFREE),LFREE,0)
         CALL RDCREF(WRK(KCREF))
         CALL GETNO(WRK(KCREF),LSYM,WRK(KOCC),WRK(KUNO),WRK(KCMO),
     *              .TRUE.,.TRUE.,.TRUE.,MOISNO,WRK(KCINDX),
     *              WRK,KFREE,LFREE,LUW4,IPRNO)
C        CALL GETNO(CVEC,ICSYM,OCCNO,UNO,CMO,ORDER,NATORB,NOAVER,
C    *              MOISNO,INDXCI,WRK,KFRSAV,LFRSAV,LUPRI,IPRINT)
         CALL MEMREL('SIROUT.MOLDEN.GETNO',WRK,KFRSAV,KREL1,KFREE,LFREE)
       END IF
       KREL2 = KFREE
       CALL MEMGET('REAL', KAOSO,  N2BASX,      WRK, KFREE, LFREE)
       CALL MEMGET('REAL', KUNPAK, NORBT*NBAST, WRK, KFREE, LFREE)
       CALL MEMGET('REAL', KOVEC,  NBAST,       WRK, KFREE, LFREE)
       CALL MOMOS(1,WRK(KCMO), WRK(KOCC), WRK(KAOSO),
     $            WRK(KUNPAK), WRK(KOVEC))
       CALL MEMREL('SIROUT.MOMOS',WRK,KREL2,KREL2,KFREE,LFREE)
      END IF
C     
C ***  Print final MCSCF orbitals ***
C
      IF (IPRI4 .GE. 5) P4FLAG(5) = .TRUE.
      IF (P4FLAG(5)) THEN
         IF (.NOT.FLAG(3)) THEN
            REWIND LUIT1
            IF ( .NOT. FNDLAB(TABLE(3),LUIT1) ) THEN
               NWARN = NWARN + 1
               WRITE (LUPRI,'(//A/)')
     *            ' WARNING, new orbitals not found on SIRIUS.RST,'//
     *            ' will try "OLDORB "'
               REWIND LUIT1
               CALL MOLLAB(TABLE(2),LUIT1,LUERR)
            END IF
            CALL READT(LUIT1,NCMOT,WRK(KCMO))
         ELSE
            IF (NASHT .GT. 1) THEN
               WRITE(LUW4,'(//A/A/)') ' Natural orbitals :',
     *                                ' =================='
            ELSE
               WRITE(LUW4,'(//A/A/)') ' Canonical orbitals :',
     *                                ' ===================='
            END IF
         END IF
         IF (IPRI4 .GT. 10) THEN
            CALL PRORB(WRK(KCMO),.FALSE.,LUW4)
         ELSE
            CALL PRORB(WRK(KCMO),.TRUE.,LUW4)
         END IF
C        CALL PRORB(CMO,PROCC,IOUT)
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
            WRITE (LUW4,1077) LCINDX-LFREE
            GO TO 110
         END IF
         CALL MEMGET('REAL',KCINDX,LCINDX,WRK,KFREE,LFREE)
         CALL MEMGET('WORK',KWRK2,LWRK2,WRK,KFREE,LFREE)
C
         REWIND LUIT1
         CALL MOLLB2 (TABLE(4),RTNLBL,LUIT1,LUERR)
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
                  WRITE (LUW4,'(//A/A/)')
     &               ' WARNING, internal error in SIROUT subroutine',
     &               ' CICALC false, but CISAVE label on SIRIUS.RST'
               END IF
               CICALC = .TRUE.
               IF (I .GT. 0 .AND. I .NE. ISTACI) THEN
                  NWARN = NWARN + 1
                  WRITE (LUW4,'(//,(A,I5/))')
     *         ' WARNING from SIROUT: ISTACI specied in input:',ISTACI,
     *         '                      ISTACI read from LUIT1 :',I
               END IF
            ELSE
               IF (CICALC) THEN
                  NWARN = NWARN + 1
                  WRITE (LUW4,'(//A/A/)')
     &               ' WARNING, internal error in SIROUT subroutine',
     &               ' CICALC true, but no CISAVE label on SIRIUS.RST'
                  CICALC = .FALSE.
               END IF
               IF (ABS(I) .NE. ISTATE) THEN
                  NWARN = NWARN + 1
                  WRITE (LUW4,'(//,(A,I5/))')
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
            WRITE (LUW4,1050) THRPWF,I
            IF (I .EQ. ISTATE .AND. .NOT. CICALC) WRITE (LUW4,1051)
            IF (I .EQ. ISTACI .AND.       CICALC) WRITE (LUW4,1052)
            CALL READT(LUIT1,NCONF,WRK(KCREF))
            CALL PRWF(WRK(KCREF),WRK(KCINDX),LUW4,THRPWF)
C           CALL PRWF(C,INDXCI,IOUT,THRPWF)
  105    CONTINUE
         CALL MEMREL('SIROUT.PRWF',WRK,KFRSAV,KCINDX,KFREE,LFREE)
      END IF
  110 CONTINUE
 1050 FORMAT(//' Printout of CI-coefficients larger than',F8.5,
     *         ' for root',I3,/)
 1051 FORMAT(' (this is the reference state)'/)
 1052 FORMAT(' (this is the CI reference state)'/)
 1077 FORMAT(/' *** SIROUT, Wave function not printed as requested',
     *  /,15X,'Reason: insufficient work memory space',
     *  /,15X,'need more than',I10,' additional real*8 work space')
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
            CALL MOLLAB(TABLE(2),LUIT1,LUERR)
            WRITE (LUPRI,'(/A/)') ' "OLDORB " orbitals punched.'
         END IF
         CALL READT(LUIT1,NCMOT,WRK(KCMO))
         IF (FLAG(3)) THEN
            IPCTL = 1
         ELSE
            IPCTL = 0
            KOCC  = KCMO
         END IF
         CALL PUNMO(IPCTL,WRK(KCMO),WRK(KOCC))
C        CALL PUNMO(IPCTL,CMO,OCC)
      END IF
C
C end of SIROUT
C
      CALL MEMREL('SIROUT',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('SIROUT')
      RETURN
      END
C  /* Deck prorb */
      SUBROUTINE PRORB(CMO,PROCC,IOUT)
C
C Written 11-Dec-1983 by ha/hjaaj
C Revised  8-Jul-1992 hjaaj
C
C Purpose:
C  Print molecular orbital coefficients on unit IOUT.
C
C Input:
C  CMO, MO orbital coefficients (symmetry blocked)
C  PROCC, print only occupied orbitals
C  IOUT, output file unit
C
#include <implicit.h>
      DIMENSION CMO(*)
      LOGICAL   PROCC
      PARAMETER (THRPRI = 1.0D-4)
C
C Used from common blocks:
C   INFINP : CENT,TYPE,SUPSYM,?
C   INFORB : NSYM,...
C   INFIND : ISSMO(),
C
#include <maxorb.h>
#include <maxash.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <r12int.h>
C
C *** FORMAT statements
C
 2400 FORMAT(/' Molecular orbitals for symmetry species',I3,
     &       /' ------------------------------------------')
 2600 FORMAT(/' Orbital   ',7I9)
 2620 FORMAT( ' Sup.sym.  ',7I9)
 2800 FORMAT(I4,1X,A4,':',A4,7F9.4)
C
C *** Print molecular orbitals
C
      ISTBAS = 0
      DO 400 ISYM = 1,NSYM
         IF (PROCC) THEN
C           NENDI = NOCC(ISYM)
            NENDI = NOCC(ISYM)
     *            + MIN(5,NORB(ISYM)-NOCC(ISYM) - NRXR12(ISYM))
         ELSE
            NENDI = NORB(ISYM)
         END IF
         NBASI = NBAS(ISYM)
         MBASI = NBAS(ISYM)
         IF (LMULBS) MBASI = MBAS1(ISYM)
      IF (NENDI.EQ.0) GO TO 300
         WRITE(IOUT,2400) ISYM
C
         ICMOI = ICMO(ISYM)
         ISTORB = IORB(ISYM)
         IEND=0
  100       IST =IEND+1
            ISTMO=IEND*NBASI+ICMOI
            IEND=IEND+7
            IF(IEND.GT.NENDI) IEND=NENDI
            IEMO=NBASI*(IEND-1)+ICMOI
            WRITE(IOUT,2600) (I,I=IST,IEND)
            IF (SUPSYM) WRITE(IOUT,2620) (ISSMO(ISTORB+I),I=IST,IEND)
            DO 200 I=1,MBASI
               JSMO=ISTMO+I
               JEMO=IEMO+I
               JJ = 0
               DO J=JSMO,JEMO,NBASI
                  IF ( ABS(CMO(J)) .GE. THRPRI ) JJ = 1
               END DO
               IF (JJ .EQ. 1) WRITE(IOUT,2800) I,CENT(I+ISTBAS),
     *            TYPE(I+ISTBAS),(CMO(J),J=JSMO,JEMO,NBASI)
  200       CONTINUE
         IF (IEND.NE.NENDI) GO TO 100
C
  300 CONTINUE
        ISTBAS = ISTBAS + NBASI
  400 CONTINUE
C
C *** End of subroutine PRORB
C
      RETURN
      END
C  /* Deck prkap */
      SUBROUTINE PRKAP (NKAPT,XKAP,PRFAC,LUPRIN)
C
C Jan 90 hjaaj (descendant of PRWOP)
C
C Purpose:
C  Write orbital operator array XKAP to unit LUPRIN.
C
C  if PRFAC  .gt. 0, write only elements with absolute value
C                    .ge. PRFAC * (max. absolute value)
C  if PRFAC  .eq. 0, write all elements;
C  if PRFAC  .lt. 0, write only elements with
C                    value .le. PRFACX * min. pos value
C                    where PRFACX = MAX(-PRFAC, -1.0/PRFAC)
C
#include <implicit.h>
      DIMENSION XKAP(*)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C Used from common blocks:
C   INFVAR : NWOPT,NWOPH,JWOPSY,JWOP(2,*)
C
#include <maxorb.h>
#include <infvar.h>
C
C Check if NKAPT legal
C
      IF (NKAPT .NE. NWOPT .AND. NKAPT .NE. NWOPH) THEN
         WRITE (LUPRIN,'(//A/A,3I10)')
     &      ' ERROR in PRKAP, NKAPT .ne. NWOPT,NWOPH',
     &      ' NKAPT, NWOPT, NWOPH :',NKAPT,NWOPT,NWOPH
         CALL QUIT(' ERROR in PRKAP, NKAPT .ne. NWOPT,NWOPH')
      END IF
C
C Any elements?
C
      IF (NKAPT.LE.0) THEN
         WRITE (LUPRIN,8000) JWOPSY
         GO TO 9999
      END IF
 8000 FORMAT(/5X,'Orbital operator symmetry =',I2,
     *       /5X,'>> NO ELEMENTS <<')
C
C Print XKAP as requested:
C
      IF (PRFAC .GT. D0) THEN
         KMAX   = IDAMAX(NKAPT,XKAP,1)
         THPKAP = PRFAC*ABS(XKAP(KMAX))
         WRITE (LUPRIN,1010) JWOPSY,PRFAC
         INPR = 0
         DO 100 IG = 1,NKAPT
            IF ( ABS(XKAP(IG)) .GT. THPKAP ) THEN
               K = JWOP(1,IG)
               L = JWOP(2,IG)
               WRITE (LUPRIN,1100) IG,K,L,XKAP(IG)
            ELSE
               INPR = INPR + 1
            END IF
  100    CONTINUE
         IF (INPR.GT.0) WRITE (LUPRIN,1210) INPR,THPKAP
      ELSE IF (PRFAC .EQ. D0) THEN
         WRITE (LUPRIN,1020) JWOPSY
         INPR = 0
         DO 200 IG = 1,NKAPT
            IF ( XKAP(IG) .NE. D0 ) THEN
               K = JWOP(1,IG)
               L = JWOP(2,IG)
               WRITE (LUPRIN,1100) IG,K,L,XKAP(IG)
            ELSE
               INPR = INPR + 1
            END IF
  200    CONTINUE
         IF (INPR.GT.0) WRITE (LUPRIN,1220) INPR
      ELSE
C     (PRFAC .LT. D0)
         KMIN   = IDAMIN(NKAPT,XKAP,1)
         PRFACX = MAX( -PRFAC, -D1/PRFAC)
         THPKAP = PRFACX*ABS(XKAP(KMIN))
         WRITE (LUPRIN,1030) JWOPSY,THPKAP
         INPR = 0
         DO 300 IG = 1,NKAPT
            IF ( XKAP(IG) .LE. THPKAP ) THEN
               K = JWOP(1,IG)
               L = JWOP(2,IG)
               WRITE (LUPRIN,1100) IG,K,L,XKAP(IG)
            ELSE
               INPR = INPR + 1
            END IF
  300    CONTINUE
         IF (INPR.GT.0) WRITE (LUPRIN,1230) INPR,THPKAP
      END IF
C
C
 9999 CONTINUE
      RETURN
C
 1010 FORMAT(/5X,'Orbital operator symmetry =',I2,/,
     *  5X,'(only elements abs greater than',1P,D8.1,
     *     ' times max abs value)'//,
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1020 FORMAT(/5X,'Orbital operator symmetry =',I2,//,
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1030 FORMAT(/5X,'Orbital operator symmetry =',I2,/,
     *  5X,'(only elements less than',1P,D10.2,')'//
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1100 FORMAT(I15,I7,I3,F20.10)
 1210 FORMAT(/1X,I8,' elements with absolute value less than',
     *       1P,D9.2,' not printed.')
 1220 FORMAT(/1X,I8,' zero elements not printed.')
 1230 FORMAT(/1X,I8,' elements with value greater than',
     *       1P,D9.2,' not printed.')
C
C -- end of subroutine PRKAP
C
      END
C  /* Deck prwop */
      SUBROUTINE PRWOP (WOP,IUNITX)
C
C Written 11-Apr-1984 -- hjaaj
C Last revision 17-Oct-1984 hjaaj (IUNITX.lt.0 option)
C                9-May-1984 hjaaj
C
C Purpose:
C  Write orbital operator array WOP to unit IUNIT.
C
C  IUNIT = IABS(IUNITX)
C
C  If IUNITX .lt. 0, write all elements;
C  if IUNITX .gt. 0, write only elements with absolute value
C                    .ge. PRFAC * (max. absolute value)
C
#include <implicit.h>
      DIMENSION WOP(*)
      PARAMETER (PRFAC = 0.1D0, D0 = 0.0D0)
C
C Used from common blocks:
C   INFVAR : NWOPT,JWOPSY,JWOP(2,*)
C
#include <maxorb.h>
#include <infvar.h>
C
      IUNIT = IABS(IUNITX)
C
C Any elements?
C
      IF (NWOPT.LE.0) GO TO 300
C
C 1. find maximum absolute value in WOP
C
      IF (IUNITX.GT.0) THEN
         MAX = IDAMAX(NWOPT,WOP,1)
         THPWOP = PRFAC*ABS(WOP(MAX))
         WRITE (IUNIT,1000) JWOPSY,PRFAC
      ELSE
         WRITE (IUNIT,1010) JWOPSY
         THPWOP = D0
      END IF
      INPR = 0
      DO 200 IG = 1,NWOPT
         IF (ABS(WOP(IG)).GE.THPWOP) THEN
            K = JWOP(1,IG)
            L = JWOP(2,IG)
            WRITE (IUNIT,1100) IG,K,L,WOP(IG)
         ELSE
            INPR = INPR + 1
         END IF
  200 CONTINUE
      IF (INPR.GT.0) WRITE (IUNIT,1200) INPR,THPWOP
C
      RETURN
C
  300 CONTINUE
      WRITE (IUNIT,3000) JWOPSY
      RETURN
C
 1000 FORMAT(/5X,'Orbital operator symmetry =',I2,/,
     *  5X,'(only elements abs greater than',1P,D8.1,
     *     ' times max abs value)'//,
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1010 FORMAT(/5X,'Orbital operator symmetry =',I2,//,
     *  5X,' Index(r,s)     r  s      operator value',/,
     *  5X,' ----------    -- --      --------------')
 1100 FORMAT(I15,I7,I3,F20.10)
 1200 FORMAT(/1X,I8,' elements with absolute value less than',
     *       1P,D9.2,' not printed.')
 3000 FORMAT(/5X,'Orbital operator symmetry =',I2,
     *       /5X,'>> NO ELEMENTS <<')
C
C -- end of subroutine PRWOP
C
      END
C  /* Deck punmo */
        SUBROUTINE PUNMO(IPCTL,CMO,OCC)
C
C 25-May-1985 hjaaj
C l.r. 910715-hjaaj (added formats 7010-7050)
C
C Purpose:
C   punch mo coefficients and, conditionally, occupation numbers
C   to LUPNCH
C
C Control:
C   IPCTL .eq. 1: also punch occupation numbers
C   else        : do not punch occupation numbers
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION CMO(*),OCC(*)
      PARAMETER (D0 = 0.0D0)
C
C Used from common blocks:
C   INFORB : NSYM,NORB(*),NBAS(*)
C
#include <inforb.h>
C
      CHARACTER*8 LBLDAT(2)
C
      CALL QENTER('PUNMO')
      CALL GETDAT(LBLDAT(1),LBLDAT(2))
      LUPNCH = -1
      CALL GPOPEN(LUPNCH,'DALTON.MOPUN','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      REWIND LUPNCH
      IF (IPCTL.EQ.1) THEN
         WRITE (LUPNCH,'(A,A8,2X,A8,A)')
     *      '**NATORB   ( punched by SIRIUS ',LBLDAT,' )'
      ELSE
         WRITE (LUPNCH,'(A,A8,2X,A8,A)')
     *      '**MOLORB   ( punched by SIRIUS ',LBLDAT,' )'
      END IF
C
      IEND =0
      DO 700 ISYM=1,NSYM
         NORBI=NORB(ISYM)
      IF (NORBI.EQ.0) GO TO 700
         NBASI=NBAS(ISYM)
         DO 705 NI=1,NORBI
            IST=IEND+1
            IEND=IEND+NBASI
            IMOMXV = IDAMAX(NBASI,CMO(IST),1)
            CMOMXV = ABS(CMO(IST-1+IMOMXV))
            IF (CMOMXV .LT. 1.0D2) THEN
               WRITE (LUPNCH,7000) (CMO(NP),NP=IST,IEND)
            ELSE IF (CMOMXV .LT. 1.0D3) THEN
               WRITE (LUPNCH,7010) (CMO(NP),NP=IST,IEND)
            ELSE IF (CMOMXV .LT. 1.0D4) THEN
               WRITE (LUPNCH,7020) (CMO(NP),NP=IST,IEND)
            ELSE IF (CMOMXV .LT. 1.0D5) THEN
               WRITE (LUPNCH,7030) (CMO(NP),NP=IST,IEND)
            ELSE IF (CMOMXV .LT. 1.0D6) THEN
               WRITE (LUPNCH,7040) (CMO(NP),NP=IST,IEND)
            ELSE
               WRITE (LUPNCH,7050) (CMO(NP),NP=IST,IEND)
            END IF
  705    CONTINUE
         NDELI = NBASI-NORBI
      IF (NDELI.EQ.0) GO TO 700
         DO 706 NI=1,NDELI
            WRITE (LUPNCH,7000) (D0,NP=1,NBASI)
  706    CONTINUE
  700 CONTINUE
      CALL GPCLOSE(LUPNCH,'KEEP')
C
      IF (IPCTL.EQ.1) THEN
         WRITE (LUPNCH,'(A,A8,2X,A8,A)')
     *   '**NATOCC   ( punched by SIRIUS ',LBLDAT,' )'
         WRITE(LUPNCH,7000) (OCC(NP),NP=1,NORBT)
      END IF
C
      CALL QEXIT('PUNMO')
      RETURN
C
 7000 FORMAT(4F18.14)
 7010 FORMAT(4F18.13)
 7020 FORMAT(4F18.12)
 7030 FORMAT(4F18.11)
 7040 FORMAT(4F18.10)
 7050 FORMAT(4G18.11)
C
C end of PUNMO
C
      END
