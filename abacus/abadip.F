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
C  /* Deck dipinp */
      SUBROUTINE DIPINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 7)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
#include <cbidip.h>
      DATA TABLE /'.SKIP  ', '.PRINT ', '.NODC  ', '.NODV  ',
     *            'xXXXXXX', 'xXXXXXX', '.STOP  '/
C
      NEWDEF = (WORD .EQ. '*DIPCTL')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in DIPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in DIPINP')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  NODC = .TRUE.
               GO TO 100
    4          CONTINUE
                  NODV = .TRUE.
               GO TO 100
    5             CONTINUE
               GO TO 100
    6             CONTINUE
               GO TO 100
    7             CUT   = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in DIPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in DIPINP')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for DIPCTL:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' DIPCTL skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in DIPCTL:',IPRINT
            END IF
            IF (NODC) WRITE (LUPRI,'(/,2A)') ' Inactive one-electron',
     *      ' density matrix neglected in DIPCTL.'
            IF (NODV) WRITE (LUPRI,'(/,2A)') ' Active one-electron',
     *      ' density matrix neglected in DIPCTL.'
            IF (TEST) WRITE (LUPRI,'(/,2A)') ' Test for dipole moments',
     *      ' and dipole reorthonormalization.'
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after DIPCTL.'
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck dipini */
      SUBROUTINE DIPINI
C
C     Initialize /CBIDIP/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbidip.h>
      IPRINT = IPRDEF
      NODC   = .FALSE.
      NODV   = .FALSE.
      TEST   = .FALSE.
      SKIP = .NOT.(DIPDER .OR. POLAR)
      CUT  = .FALSE.
C
      RETURN
      END
C  /* Deck dipctl */
      SUBROUTINE DIPCTL(WORK,LWORK,PASS)
C
C     tuh 1985
C     Rewritten for symmetry January 1990, tuh
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      LOGICAL PASS
      DIMENSION WORK(LWORK)
#include <cbidip.h>
#include <nuclei.h>
#include <inforb.h>
#include <infdim.h>
C
      IF (SKIP) RETURN
      CALL TIMER('START ',TIMSTR,TIMEND)
      IF (IPRINT .GT. 0) CALL TITLER('Output from DIPCTL','*',103)
      KCMO    = 1
      KDVSP   = KCMO   + NCMOT
      KDV     = KDVSP  + NNASHX
      KDPSOP  = KDV    + N2ASHX
      KDPSO   = KDPSOP + 3*NNBASX
      KDPRHS  = KDPSO  + N2BASX
      KDPREO  = KDPRHS + NVARMA
      KWRK    = KDPREO + 3*NUCDEP
      LWRK    = LWORK  - KWRK + 1
      IF (KWRK .GE. LWORK) CALL STOPIT('DIPCTL',' ',KWRK,LWORK)
      CALL DIPCT1(WORK(KCMO),WORK(KDVSP),WORK(KDV),
     *            WORK(KDPSOP),WORK(KDPSO),
     *            WORK(KDPRHS),WORK(KDPREO),
     *            WORK(KWRK),LWRK)
      IF (IPRINT .GT. 1) CALL TIMER('DIPCTL',TIMSTR,TIMEND)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/,A)')
     &          ' Program stopped after DIPCTL as required.'
         WRITE (LUPRI,'(A)') ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in DIPCTL) *****')
      END IF
      RETURN
      END
C  /* Deck dipct1 */
      SUBROUTINE DIPCT1(CMO,DVSP,DV,SOPACK,SOINT,DIPRHS,DIPREO,
     *                  WORK,LWORK)
C
C     December 1989, tuh
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <iratdef.h>
      PARAMETER (DM1 = -1.0D0)
#include <dipole.h>
#include <orgcom.h>
#include <moldip.h>
#include <cbidip.h>
C
#include <symmet.h>
#include <dorps.h>
#include <abainf.h>
#include <linaba.h>
#include <inftap.h>
#include <infvar.h>
#include <inforb.h>
#include <infdim.h>
#include <infinp.h>
#include <inflin.h>
C
      CHARACTER*8 LABINT(3*MXCOOR)
      DIMENSION CMO(NCMOT), DVSP(NNASHX), DV(NASHT,NASHT),
     *          SOPACK(NNBASX,3), DIPRHS(NVARMA), DIPREO(3*NUCDEP),
     *          WORK(LWORK), SOINT(NBAST,NBAST)
C
      LOGICAL OLDDX
#include <chrxyz.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from DIPCT1','*',103)
C
C     ***** Read orbitals *****
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC)
      CALL READI(LUSIFC,IRAT*NCMOT,CMO)
C
C     ***** Read one-electron density *****
C
      IF (NASHT .GT. 0) THEN
         READ (LUSIFC)
         CALL READI(LUSIFC,IRAT*NNASHX,DVSP)
         CALL DSPTSI(NASHT,DVSP,DV)
         IF (IPRINT .GT. 10) THEN
            CALL HEADER('Active density matrix in DIPCT1',-1)
            CALL OUTPUT(DV,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         END IF
      END IF
C
C     ***** Read SO integrals *****
C
      NCOMP  = 3
      NPATOM = 0
      KINT   = 1
      KREP   = KINT + (9*MXCENT + 1)/IRAT
      KLAST  = KREP + (9*MXCENT + 1)/IRAT
      LWRK   = LWORK - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('DIPCT1','GET1IN',KLAST,LWORK)
      CALL GET1IN(SOPACK,'DIPLEN ',NCOMP,WORK(KLAST),LWRK,LABINT,
     &            WORK(KINT),WORK(KREP),IDUMMY,.FALSE.,NPATOM,.TRUE.,
     &            DUMMY,.FALSE.,DUMMY,IPRINT)
      IF (NNBASX .GT. LWORK) CALL STOPIT('DIPCT1','SYMUPK',NNBASX,LWORK)
      IF (IPRINT .GT. 10) THEN
         CALL HEADER('Dipole (x) SO  matrix in DIPCT1',-1)
         CALL OUTPAK(SOPACK(1,1),NBAST,1,LUPRI)
         CALL HEADER('Dipole (y) SO  matrix in DIPCT1',-1)
         CALL OUTPAK(SOPACK(1,2),NBAST,1,LUPRI)
         CALL HEADER('Dipole (z) SO  matrix in DIPCT1',-1)
         CALL OUTPAK(SOPACK(1,3),NBAST,1,LUPRI)
      END IF
C
C     ***** Open file for right-hand side of response equations *****
C
      CALL GPOPEN(LUGDR,ABAGDR,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,OLDDX)
C
C     ********************************
C     ***** Loop over components *****
C     ********************************
C
      IOFFAX = 0
      DO 100 IREP = 0, MAXREP
         NAXIS = NAXREP(IREP,1)
         IF (DOREPS(IREP) .AND. (NCRREP(IREP,1).GT.0)) THEN
            IF (IPRINT.GT.5) WRITE (LUPRI,'(1X,A,I5)')' Symmetry ',IREP
C
C           Number of orbital and configuration variables
C
            CALL ABAVAR(IREP+1,.FALSE.,IPRINT,WORK(KLAST),LWRK)
            IF (NVARPT .GT. 0) THEN
               DO 200 IAX = 1, NAXIS
                  ICOOR = IPTXYZ(IAX,IREP,1)
                  CALL DSCAL(NNBASX,DM1,SOPACK(1,ICOOR),1)
                  CALL DSPTSI(NBAST,SOPACK(1,ICOOR),SOINT)
                  IF (IPRINT .GT. 5) THEN
                     CALL AROUND('Component of dipole moment:'
     &                           //CHRXYZ(ICOOR))
                     IF (IPRINT .GT. 10) THEN
                        CALL HEADER('SO matrix in DIPCT1',-1)
                        CALL OUTPUT(SOINT,1,NBAST,1,NBAST,NBAST,NBAST,
     &                              1,LUPRI)
                     END IF
                  END IF
C
C                 ***** Reorthonormalization and right-hand side *****
                  CALL OPGCTL(DIPREO,DIPRHS,CMO,DV,SOINT,DIPACT,DIPICT,
     &                           WORK,LWORK,IREP,NODC,NODV,IPRINT)
                  CALL CHKPRP(' Testing property: '//CHRXYZ(-ICOOR)//
     &                        ' component of dipole moment',
     &                         DIPME(ICOOR),DIPICT,DIPACT,IPRINT)
C
C                 ***** Add reorthonormalization to DDIPS *****
C
                  IF (.NOT. HELFEY) CALL DCOPY(3*NUCDEP,
     &                                   DIPREO,1,DDIPS(IOFFAX+IAX,1),3)
C
C                 ***** Write right-hand side on file *****
C
                  IDISK = 3*NUCDEP + IOFFAX + IAX
                  CALL WRITDX(LUGDR,IDISK,IRAT*NVARPT,DIPRHS)
  200          CONTINUE
            END IF
         END IF
         IOFFAX = IOFFAX + NAXIS
  100 CONTINUE
C
C     ***** Print static constribution to dipole gradient *****
C
      CALL GPCLOSE(LUGDR,'KEEP')
      IF (IPRINT .GT. 1 .AND. DIPDER) THEN
         KCSTRA = 1
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         IF (KLAST .GT. LWORK)
     &        CALL STOPIT('DIPCT1','TRANUC',KLAST,LWORK)
         CALL HEADER('Reorthonormalization part of dipole gradient',-1)
         CALL FCPRI(DDIPS,'APT',WORK(KCSTRA),WORK(KSCTRA))
         CALL HEADER('Static contribution to dipole gradient',-1)
         CALL DZERO(DIP1,9*NUCDEP)
         CALL DIPADD(DDIPN)
         CALL DIPADD(DDIPE)
         CALL DIPADD(DDIPS)
         CALL FCPRI(DIP1,'APT',WORK(KCSTRA),WORK(KSCTRA))
         CALL DZERO(DIP1,9*NUCDEP)
      END IF
      RETURN
      END
C  /* Deck dipadd */
      SUBROUTINE DIPADD(AMAT)
#include <implicit.h>
#include <mxcent.h>
      DIMENSION AMAT(3,MXCOOR)
#include <nuclei.h>
#include <moldip.h>
      NCOORD = 3*NUCDEP
      DO 100 I = 1, 3
         DO 200 J = 1, NCOORD
            DIP1(I,J) = DIP1(I,J) + AMAT(I,J)
  200    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck polpri */
      SUBROUTINE POLPRI(AMAT,SPC,ITYPE1)
C
C   ITYPE = 1 : electric polarizability
C   ITYPE = 2 : magnetic susceptibility
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <codata.h>
      PARAMETER (D1 = 1.0 D00)
      LOGICAL ALLREP, INERT, LARGE
      DIMENSION AMAT(3,3), IND(3)
      CHARACTER*(*) SPC
      CHARACTER LAB(3)*2, CHRABC(3)*1
#include <inirep.h>
#include <abainf.h>
#include <inftap.h>
#include <dorps.h>
#include <symmet.h>
#include <chrxyz.h>
      DATA CHRABC /'A','B','C'/
C
      INERT = INDEX(SPC,'PRI') .NE. 0
      LARGE = ITYPE1 .LT. 0
      ITYPE = ABS(ITYPE1)
C
C     *****************************************
C     ***** Units (atomic or Angstrom**3) *****
C     *****************************************
C
      IF (INDEX(SPC,'EXP') .NE. 0) THEN
         FAC = XTANG*XTANG*XTANG
      ELSE IF (INDEX(SPC,'CGS') .NE. 0) THEN
         FAC = XTANG*DEBYE
      ELSE IF (INDEX(SPC,'SIU') .NE. 0) THEN
         FAC = XTANGM10*XTANGM10*ECHARGE*1.0D+40
      ELSE
         FAC = D1
      END IF
C
C     ***************************************
C     ***** Pick up components to print *****
C     ***************************************
C
      IF (.NOT.INERT) THEN
         NCOMP = 0
         DO 100 ICOOR = 1, 3
            IF (DOSYM(ISYMAX(ICOOR,ITYPE) + 1)
     &          .OR. ITYPE .EQ. 2) THEN
               NCOMP = NCOMP + 1
               IND(NCOMP) = IPTAX(ICOOR,ITYPE)
               IF (ITYPE .EQ. 1) THEN
                  LAB(NCOMP) = 'E'//CHRXYZ(-ICOOR)
               ELSE
                  LAB(NCOMP) = 'B'//CHRXYZ(-ICOOR)
               END IF
            END IF
  100    CONTINUE
      ELSE
         NCOMP = 0
         DO 200 ICOOR = 1, 3
            ALLREP = .TRUE.
            DO 210 IREPS = 1, NREPPI(ICOOR)
               IF (.NOT.DOREPS(IREPPI(ICOOR,IREPS))) ALLREP = .FALSE.
  210       CONTINUE
            IF (ALLREP) THEN
               NCOMP = NCOMP + 1
               IND(NCOMP) = ICOOR
               LAB(NCOMP) = 'E'//CHRABC(ICOOR)
            END IF
  200    CONTINUE
      END IF
C
C     *****************
C     ***** Print *****
C     *****************
C
      IF (NCOMP .GT. 0) THEN
         IF (INERT) WRITE (LUPRI,'(15X,A/)')
     &      '(Along principal axes of moments of inertia)'
         IF (LARGE) THEN
            WRITE (LUPRI,'(5X,3(18X,A2)/)') (LAB(I),I=1,NCOMP)
         ELSE
            WRITE (LUPRI,'(15X,3(10X,A2)/)') (LAB(I),I=1,NCOMP)
         END IF
         DO 300 I = 1, NCOMP
            IF (LARGE) THEN
               WRITE (LUPRI, '(2X,A2,3X,3F20.12)') LAB(I),
     &                  (FAC*AMAT(IND(I),IND(J)),J=1,NCOMP)
            ELSE
               WRITE (LUPRI, '(12X,A2,3X,3F12.6)') LAB(I),
     &                  (FAC*AMAT(IND(I),IND(J)),J=1,NCOMP)
            END IF
  300    CONTINUE
csonia 04/10/95
         IF (LUCME.GT.0) THEN
            DO I = 1, NCOMP
               WRITE (LUCME, '(1X,A2,1X,3D23.15)') LAB(I),
     &                  (FAC*AMAT(IND(I),IND(J)),J=1,NCOMP)
            END DO
         END IF
csonia 04/10/95
      ELSE
         WRITE (LUPRI,'(2X,A)') ' Polarizabilities not calculated - '//
     &                       ' appropriate symmetries not requested.'
      END IF
      WRITE (LUPRI, '()')
      RETURN
      END
C  /* Deck aptpop */
      SUBROUTINE APTPOP(AMAT,SHESS,QAPT,CSTRA,SCTRA,CHESS,WRK,LWRK)
C
C     hjaaj+tuh 151289
C
C Population Analysis based on Atomic Polar Tensors
C
C according to J. Cioslowski, J.Am.Chem.Soc. 111 (1989) 8333-8336
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (THIRD = 1.0D0/3.0D0)
      DIMENSION AMAT(3,MXCOOR), CMAT(3,MXCOOR)
      DIMENSION HELP1(3,3), HELP2(3,3), HELP3(3,3), CORR(3,3)
      DIMENSION SHESS(MXCOOR,MXCOOR), CHESS(MXCOOR,MXCOOR), WRK(LWRK)
      DIMENSION QAPT(*), CSTRA(*), SCTRA(*)
      LOGICAL ALL
#include <abainf.h>
#include <nuclei.h>
#include <symmet.h>
#include <cbisol.h>
#include <chrxyz.h>
      ALL = DOSYM(ISYMAX(1,1)+1) .AND. DOSYM(ISYMAX(2,1)+1)
     &                           .AND. DOSYM(ISYMAX(3,1)+1)
      IF (ALL) THEN
         NCOOR = 3*NUCDEP
         CALL AROUND('APT Population Analysis')
         WRITE (LUPRI,'(/6X,A/)')
     *  '[ Reference : J. Cioslowski, J.Am.Chem.Soc. 111 (1989) 8333 ]'
         CALL TRADIP(AMAT,CMAT,CSTRA,SCTRA,NCOOR,1,1)
         MXCNT = NUCDEP
         IF (SOLVNT) THEN
            CALL TRAHES(SHESS,MXCOOR,CHESS,CSTRA,SCTRA,
     &                  MXCOOR,NCOOR,1)
            MXCNT = MXCNT
            ICAV  = 3*MXCNT - 3
            DO 20 I = 1, 3
            DO 20 J = 1, 3
               HELP1(I,J) = CHESS(ICAV + I,ICAV + J)
               CALL DGEINV(3,HELP1,HELP2,WRK,WRK(20),INFO)
 20         CONTINUE 
         END IF
         DO 100 IATOM = 1, MXCNT
            ICOOR = (IATOM-1)*3
            TOTCR = 0.0D0
            IF (SOLVNT) THEN
               DO 10 I = 1, 3
               DO 10 J = 1, 3
                  HELP1(I,J) = CHESS(ICOOR + I,ICAV + J)
                  HELP3(I,J) = CMAT(I,ICAV + J)
 10            CONTINUE 
               CALL DGEMM('N','N',3,3,3,1.D0,
     &                    HELP1,3,
     &                    HELP2,3,0.D0,
     &                    WRK,3)
               CALL DGEMM('N','N',3,3,3,1.D0,
     &                    WRK,3,
     &                    HELP3,3,0.D0,
     &                    CORR,3)
               TOTCR = -(CORR(1,1) + CORR(2,2) + CORR(3,3))*THIRD
            END IF
            QAPT(IATOM)=THIRD*(CMAT(1,ICOOR+1)+CMAT(2,ICOOR+2)
     &                        +CMAT(3,ICOOR+3))
            QAPT(IATOM) = QAPT(IATOM) + TOTCR
            WRITE (LUPRI,'(28X,A6,F12.6)') NAMDEP(IATOM), QAPT(IATOM)
 100     CONTINUE
         WRITE (LUPRI,'()')
      END IF
      IF (SOLVNT .AND. .NOT. MOLHES) THEN
         WRITE (LUPRI,'(/,2A)') ' WARNING: APT population incorrect '//
     &        'because the full molecular Hessian',
     &        '          is needed to project the charge on the '//
     &        'cavity center'
         END IF
      RETURN
      END
C  /* Deck dp0sum */
      SUBROUTINE DP0SUM
#include <implicit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0)
#include <dipole.h>
#include <moldip.h>
#include <symmet.h>
#include <pcmlog.h>
      DO 100 I = 1, 3
         IF (ISYMAX(I,1) .EQ. 0) THEN
            DIP0(I) = DIPMN(I) + DIPME(I)
         ELSE
            DIP0(I) = D0
         END IF
  100 CONTINUE
Clf Local field corrected dipole moment
      IF (PCM.AND.LOCFLD) THEN
         DO 200 I = 1, 3
            IF (ISYMAX(I,1) .EQ. 0) THEN
               DIPLF0(I) = DLFN(I) + DLFE(I)
            ELSE
               DIPLF0(I) = D0
            END IF
 200     CONTINUE
      END IF
      RETURN
      END
C  /* Deck chkprp */
      SUBROUTINE CHKPRP(TEXT,PRPDEN,PRPICT,PRPACT,IPRINT)
C
C     Test CI part of gradient:
C
C     Compare one-electron property calculated from contraction of
C     densities with integrals (PRPDEN) and from contraction of CI
C     reference with CI property gradient
C
C     tuh 080190
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (THRS = 1.0D-10)
      CHARACTER TEXT*(*)
#include <abainf.h>
C
      IF (IPRINT .GE. 5) THEN
         CALL TITLER('Output from CHKPRP','*',103)
      END IF
      PRPCI  = PRPACT + PRPICT
      DIFFER = PRPCI  - PRPDEN
      ABSDIF = ABS(DIFFER)
      IF (ABSDIF .GE. THRS) THEN
         CALL HEADER('WARNING: Dif. between prop. calculated'//
     *               ' from density and CI gradient!',0)
         NWNABA = NWNABA + 1
      END IF
      IF ((ABSDIF .GE. THRS) .OR. (IPRINT .GE. 4)) THEN
         CALL HEADER(TEXT,-1)
         WRITE (LUPRI,'(/,A,//)')
     *    '         Active (CI)          Inactive             Total'
         WRITE (LUPRI,'(3F20.10)') PRPACT,PRPICT,PRPCI
         WRITE (LUPRI,'(/,A,//)')
     *    '         CI property      Density property       Difference'
         WRITE (LUPRI,'(3F20.10)') PRPCI, PRPDEN, DIFFER
      END IF
      RETURN
      END
C  /* Deck diapol */
      SUBROUTINE DIAPOL
C
C     Principal values and axes of polarizability
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER ( D0 = 0.0D0 , D2 = 2.0D0, D3 = 3.0D0)
      DIMENSION AXES(3,3), VALUES(6), WORK(3), IWORK(3)
#include <symmet.h>
#include <dorps.h>
#include <moldip.h>
#include <chrxyz.h>
      NAXES = 0
      DO 100 IREP = 0, MAXREP
         IF (DOREPS(IREP) .AND. (NAXREP(IREP,1).GT.0)) THEN
            NAXES = NAXES + NAXREP(IREP,1)
         END IF
  100 CONTINUE
      IF (NAXES .GT. 0) THEN
         CALL AROUND('Principal values and axes of polarizability (au)')
         IF (MAXREP .GT. 0) THEN
            CALL HEADER('sym         value             '/
     &                 /'x          y          z  ',8)
         ELSE
            CALL HEADER
     &        ('   value              x          y          z  ',12)
         END IF
         CALL DUNIT(AXES,3)
         ISTR = 1
         DO 200 IREP = 0, MAXREP
            NAXIS = NAXREP(IREP,1)
            IF (DOREPS(IREP) .AND. (NAXIS.GT.0)) THEN
               IJ = ISTR
               DO 300 I = ISTR, ISTR + NAXIS - 1
                  DO 310 J = ISTR, I
                     VALUES(IJ) = POLFLT(I,J)
                     IJ = IJ + 1
  310             CONTINUE
  300          CONTINUE
               CALL JACO(VALUES(ISTR),AXES(ISTR,ISTR),NAXIS,3,3,
     &                   WORK,IWORK)
               DO 400 I = 1, NAXIS
                  VALUES(ISTR + I - 1) = VALUES(ISTR + (I*(I+1)/2) - 1)
  400          CONTINUE
               CALL ORDER2(AXES(1,ISTR),VALUES(ISTR),NAXIS,3)
               DO 500 I = ISTR, ISTR + NAXIS - 1
                  IF (MAXREP .GT. 0) THEN
                     WRITE (LUPRI,'(9X,I2,3X,F12.4,5X,3F11.4)')
     &                  IREP+1,VALUES(I),(AXES(IPTAX(J,1),I),J=1,3)
                  ELSE
                     WRITE (LUPRI,'(10X,F12.4,5X,3F11.4)')
     &                  VALUES(I), (AXES(IPTAX(J,1),I),J=1,3)
                  END IF
  500          CONTINUE
            END IF
            ISTR = ISTR + NAXIS
  200    CONTINUE
      END IF
      IF (NAXES .EQ. 3) THEN
         PMEAN = (VALUES(1) + VALUES(2) + VALUES(3))/D3
         ANIS1 = (VALUES(1) - VALUES(2))**2
         ANIS2 = (VALUES(2) - VALUES(3))**2
         ANIS3 = (VALUES(3) - VALUES(1))**2
         ANIS  = SQRT((ANIS1 + ANIS2 + ANIS3)/D2)
         WRITE (LUPRI,'(//11X,A,F10.4,7X,A,F10.4)')
     *         ' Mean value:', PMEAN, ' Anisotropy:', ANIS
      END IF
      WRITE (LUPRI,'(//)')
      RETURN
      END
C  /* Deck trapol */
      SUBROUTINE TRAPOL(POLAR,KEY)
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      CHARACTER*(*) KEY
      DIMENSION POLAR(3,3), AMAT(3,3)
#include <symmet.h>
      DO 100 I = 1, 3
        DO 200 J = 1, 3
           IF (KEY .EQ. 'TOSYM') THEN
              AMAT(IPTAX(I,1),IPTAX(J,1)) = POLAR(I,J)
           ELSE IF (KEY .EQ. 'FROMSYM') THEN
              AMAT(I,J) = POLAR(IPTAX(I,1),IPTAX(J,1))
           ELSE
               WRITE (LUPRI,'(//,3A,/,A)')
     *            ' Keyword ',KEY,' unknown in TRAPOL.',
     *            ' Program cannot continue.'
               CALL QUIT('Illegal keyword in TRAPOL')
           END IF
  200   CONTINUE
  100 CONTINUE
      CALL DCOPY(9,AMAT,1,POLAR,1)
      RETURN
      END
