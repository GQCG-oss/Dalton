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
C  /* Deck ortrel */
      SUBROUTINE ORTREL(WORK,LWORK,PASS)
C
C      7-Nov-1988 Hans Joergen Aa. Jensen
C
C     Purpose:
C      Calculate geometry response gradients and
C      reorthonormalization terms (apart from tr(Sab*F))
C
C      This new driver routine calls MCORL or CCORL,
C      depending on wave function type.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
      LOGICAL PASS
      DIMENSION WORK(LWORK)
C
C      Used from common blocks:
C       /ABAINF/: MOLHES,DIPDER
C
#include <abainf.h>
C
C     Return immediately if not (molhes or dipder or qpgrad),
C        we do not want to open LUSFDA and LUGDR without using them.
C
      IF (.NOT. (MOLHES .OR. DIPDER .OR. QPGRAD)) RETURN
C
      CALL MCORL(WORK,LWORK,PASS)
C
      RETURN
      END
C  /* Deck mcorl */
      SUBROUTINE MCORL(WORK,LWORK,PASS)
C
C     Written 23-jan-1985 Hans Joergen Aa. Jensen
C     Modified 14-jun-1985 TUH
C     Modified for symmetry 02-nov-1988 tuh
C
C     Purpose:
C      Calculate geometry response gradients and
c      reorthonormalization terms (apart from tr(Sab*F))
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
C
      LOGICAL RUNORT, OLDRD, OLDDX, PASS, DOREAL, DOIMAG, DOSFD,
     &        DOBRHS, DOPRHS
      DIMENSION WORK(LWORK)
C
C      Used from common blocks:
C       /ABAINF/: POLAR
C       /INFORB/: NNORBT
C       /INFDIM/: NVARMA
C       /NUCLEI/: NUCDEP
C       /INFTAP/: LUGDR,LURDR, LRGDR,LRRDR, NBGDR,NBRDR
C
#include <abainf.h>
#include <dummy.h>
#include <cbisol.h>
#include <inflin.h>
#include <infvar.h>
#include <infdim.h>
#include <inforb.h>
#include <nuclei.h>
#include <inftap.h>
C
      DOSFD  = ((MOLHES .OR. DIPDER .OR. QPGRAD) .AND. .NOT. HELFEY)
     &         .OR. MAGSUS
      DOREAL = MOLHES .OR. DIPDER .OR. QPGRAD .OR. VCD 
      DOIMAG = SHIELD .OR. MAGSUS .OR. VCD .OR. SPNSPN .OR. ECD .OR.
     &         VROA .OR. MOLGFA .OR. SPINRO .OR. OPTROT
      DOBRHS = SHIELD .OR. MAGSUS .OR. VCD .OR. ECD .OR. MOLGFA .OR.
     &         VROA .OR. SPINRO .OR. OPTROT
      DOPRHS = SHIELD .OR. SPNSPN .OR. SPINRO
      IF (.NOT. (DOREAL .OR. DOIMAG)) RETURN
      CALL QENTER('MCORL')
      KTEST = 1
      WORK(KTEST) = -999.9D0
C
C     Unit for reorthonormalization
C     =============================
C
      KLAST = KTEST + 1
      IF (DOSFD .AND. (MOLHES .OR. DIPDER .OR. QPGRAD)) THEN
         CALL GPOPEN(LUSFDA,ABASF,'UNKNOWN','DIRECT',' ',IRAT*N2ORBX,
     &               OLDDX)
         KHESFS = KLAST
         KLAST  = KHESFS + MXCOOR*MXCOOR
      ELSE
         KHESFS = KTEST
      END IF
      LWRK = LWORK - KLAST + 1
C
C     Unit LUGDR for real singlet RHS
C     ===============================
C
      IF (DOREAL) THEN
         NREC = 3*NUCDEP
         IF (POLAR) NREC = NREC + 3
         IF (QPGRAD .AND. .NOT. MOLHES) NREC = NREC + 6
         CALL GPOPEN(LUGDR,ABAGDR,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,
     &               OLDDX)
         IF (OLDDX) THEN
            WRITE (LUPRI,'(/A)') ' Old LUGDR file opened in ORTREL.'
         ELSE
            CALL DZERO(WORK(KLAST),NVARMA)
            DO 100 I = 1, NREC
               CALL WRITDX (LUGDR,I,IRAT*NVARMA,WORK(KLAST))
  100       CONTINUE
         END IF
      END IF
C
C     Unit LUGDI for imaginary singlet RHS
C     ====================================
C
      IF (DOIMAG) THEN
         NREC = 0
         IF (.TRUE. .OR. SHIELD .OR. SPNSPN .OR. SPINRO) 
     &       NREC = NREC + 3*NUCDEP
         IF (SHIELD .OR. MAGSUS .OR. VCD .OR. ECD .OR. MOLGFA .OR.
     &       VROA .OR. SPINRO .OR. OPTROT) NREC = NREC + 3
         CALL GPOPEN(LUGDI,ABAGDI,'UNKNOWN','DIRECT',' ',IRAT*NVARMA,
     &               OLDDX)
         IF (OLDDX) THEN
            WRITE (LUPRI,'(/A)') ' Old LUGDI file opened in ORTREL.'
         ELSE
            CALL DZERO(WORK(KLAST),NVARMA)
            DO 200 I = 1, NREC
               CALL WRITDX (LUGDI,I,IRAT*NVARMA,WORK(KLAST))
  200       CONTINUE
         END IF
      END IF
C
C     Unit for solvent TLM
C     ====================
C
      IF (SOLVNT .AND. DOSFD) CALL GPOPEN(LUTLM,ABATLM,'UNKNOWN',
     &                                    'DIRECT',' ',IRAT*LMTOT,OLDDX)
C
C     ************************************************************
C     ***** Calculate differentiated gradients and Y vectors *****
C     ***** Get differentiated overlap matrices in MO basis  *****
C     ************************************************************
C
      IF (MOLHES .OR. DIPDER .OR. QPGRAD) CALL HESZER(WORK(KHESFS))
      CALL RHSIDE(DOSFD,DOREAL,DOIMAG,DOBRHS,DOPRHS,WORK(KHESFS),
     &            WORK(KLAST),LWRK)
C
      IF (DOREAL) CALL GPCLOSE(LUGDR,'KEEP')
      IF (DOIMAG) CALL GPCLOSE(LUGDI,'KEEP')
C
C     ******************************************************************
C     ***** Calculate reorthonormalization contribution to Hessian *****
C     ******************************************************************
C
C     Evaluate tr(Sa*Yb) contributions to HESMOL(a,b)
C
      IF (DOSFD) THEN
         IF (MOLHES .OR. DIPDER) THEN
            KY     = KLAST
            KSD    = KY     + N2ORBX
            KCSTRA = KSD    + N2ORBX
            KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
            KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
            KGLM   = KLAST
            IF (SOLVNT) THEN
               KTLMA  = KGLM  + LMTOT
               KTLMB  = KTLMA + LMTOT
               KHSOLT = KTLMB + LMTOT
               KLAST  = KHSOLT+ MXCOOR*MXCOOR
            ELSE
               KTLMA  = KTEST
               KTLMB  = KTEST
               KHSOLT = KTEST
            END IF
            IF (KLAST .GT. LWRK) CALL STOPIT('MCORL',' ',KLAST,LWRK)
            CALL REORT(PASS,WORK(KSD),WORK(KY),WORK(KGLM),WORK(KHESFS),
     &                 WORK(KHSOLT),WORK(KTLMA),WORK(KTLMB),
     &                 WORK(KCSTRA),WORK(KSCTRA),LUSFDA,LUTLM,N2ORBX)
            CALL GPCLOSE(LUSFDA,'DELETE')
            IF (SOLVNT) CALL GPCLOSE(LUTLM,'DELETE')
         END IF
      END IF
      IF (WORK(KTEST) .NE. -999.9D0) THEN
         CALL QUIT('programming error, WORK(KTEST) has been modified')
      END IF
      CALL QEXIT('MCORL')
      RETURN
      END
C  /* Deck ortinp */
      SUBROUTINE ORTINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 3)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <abainf.h>
#include <cbiort.h>
      DATA TABLE /'.SKIP  ', '.PRINT ', '.STOP  '/
C
      NEWDEF = (WORD .EQ. '*REORT ')
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
                     GO TO (1,2,3), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ORTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ORTINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3          CONTINUE
                  CUT  = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ORTINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ORTINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for REORT:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' REORT skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in REORT :',IPRINT
            END IF
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after REORT.'
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck ortini */
      SUBROUTINE ORTINI
C
C     Initialize /CBIORT/
C
#include <implicit.h>
#include <mxcent.h>
#include <abainf.h>
#include <cbiort.h>
C
      IPRINT = IPRDEF
      SKIP   = .NOT.MOLHES
      CUT    = .FALSE.
      RETURN
      END
C  /* Deck reort */
      SUBROUTINE REORT(PASS,SKMAT,YXMAT,GLM,HESFS1,HSOLTT,TLMA,TLMB,
     &                 CSTRA,SCTRA,LUDA,LUTLM,LEN)
C
C     tuh 021188
C
C     Purpose:
C      Evaluate reorthonormalization term (apart from tr(Sab F) ) to
C      HESMOL(a,b) = sum(ij) SKMAT(ij,a) * YXMAT(ij,b)
C      Contribution to solvent Hessian Gl*Tlma*Tlmb
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
      PARAMETER (D2 = 2.0D0)
      DIMENSION SKMAT(LEN), YXMAT(LEN), GLM(LMTOT),
     &          TLMA(LMTOT), TLMB(LMTOT), HESFS1(MXCOOR,MXCOOR),
     &          HSOLTT(MXCOOR,MXCOOR), CSTRA(*), SCTRA(*)
      LOGICAL PASS
#include <taysol.h>
#include <cbisol.h>
#include <cbiort.h>
#include <energy.h>
#include <taymol.h>
#include <inforb.h>
#include <nuclei.h>
#include <symmet.h>
#include <dorps.h>
C
      IF (SKIP) RETURN
      CALL QENTER('REORT')
      IF (IPRINT .GT. 0) THEN
         CALL TIMER('START ',TIMSTR,TIMEND)
         WRITE (LUPRI,'(A)')
     *      '  <<<<<<<<<< OUTPUT FROM REORT >>>>>>>>>> '
      END IF
C
      IF (SOLVNT) THEN
         CALL HESZER(HSOLTT)
         CALL SOLFL(GLM,EPDIEL,RCAV,LCAVMX)
         IF (IPRINT .GT. 5) THEN
            CALL HEADER('GLM constructed in REORT',-1)
            CALL OUTPUT(GLM,1,1,1,LMTOT,1,LMTOT,1,LUPRI)
         END IF
      END IF
C
C     ***** Loop over symmetries *****
C
      NOFF = 0
      DO 100 ISYM = 1, NSYM
         NDIR = NCRREP(ISYM - 1,1)
         IF (DOREPS(ISYM - 1) .AND. (NDIR.GT.0)) THEN
            IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,2I5)')
     *            ' Symmetry and number of coordinates:', ISYM, NDIR
C
C           ***** Loop over YXMAT matrices *****
C
            DO 300 ICOOR = 1, NDIR
               ICRD = NOFF + ICOOR
               IF (DOPERT(ICRD,1)) THEN
                  IREC = 2*ICRD
                  CALL READDX(LUDA,IREC,IRAT*LEN,YXMAT)
                  IF (IPRINT .GE. 10) THEN
                     WRITE (LUPRI,'(/A,2I5)') ' ICOOR, ICRD', ICOOR,ICRD
                     CALL AROUND('YMAT in REORT')
                     CALL OUTPUT(YXMAT,1,NORBT,1,NORBT,NORBT,
     &                           NORBT,1,LUPRI)
                  END IF
C
C                 ***** Loop over S/K matrices *****
C
                  DO 400 JCOOR = 1, NDIR
                     JCRD = NOFF + JCOOR
                     IF (DOPERT(JCRD,1)) THEN
                        IMAT = 1
                        IREC = 2*JCRD - 1
                        CALL READDX(LUDA,IREC,IRAT*LEN,SKMAT)
                        HESFS1(ICRD,JCRD) = HESFS1(ICRD,JCRD)
     &                                    - DDOT(LEN,SKMAT,1,
     &                                           YXMAT,1)
                     END IF
  400             CONTINUE
               END IF
  300       CONTINUE
C
C           Tlm*Tlm - contribution to solvent Hessian
C           =========================================
C
C
C     Check dimension (compare CC), Jan.-94, Kr.
C
            IF (SOLVNT) THEN
               DO 311 ICOOR = 1, NDIR
                  ICRD = NOFF + ICOOR
                  IF (DOPERT(ICRD,1)) THEN
                     CALL READDX(LUTLM,ICRD,IRAT*LMTOT,TLMA)
                     IF (IPRINT .GE. 5) THEN
                        CALL AROUND('TLMA in REORT')
                        WRITE (LUPRI,'(A,2I5)')' ICOOR,ICRD',ICOOR,ICRD
                        CALL OUTPUT(TLMA,1,1,1,LMTOT,1,LMTOT,1,LUPRI)
                     END IF
                     DO 411 JCOOR = 1, NDIR
                        JCRD = NOFF + JCOOR
                        IF (DOPERT(JCRD,1)) THEN
                          CALL READDX(LUTLM,JCRD,IRAT*LMTOT,TLMB)
                          HSOLTT(ICRD,JCRD) =
     &                                 D2*DV3DOT(LMTOT,GLM,TLMA,TLMB)
                        END IF
 411                 CONTINUE
                  END IF
 311           CONTINUE
            END IF
         END IF
         NOFF = NOFF + NDIR
  100 CONTINUE
C
C     ***** Print of matrix before symmetrization *****
C
      IF (IPRINT .GE. 15) THEN
         CALL AROUND('HESFS1 before symmetrization')
         CALL OUTPUT(HESFS1,1,3*NUCDEP,1,3*NUCDEP,MXCOOR,MXCOOR,1,LUPRI)
         IF (SOLVNT) THEN
            CALL AROUND('Solvent Hessian HSOLTT')
            CALL OUTPUT(HSOLTT,1,3*NUCDEP,1,3*NUCDEP,MXCOOR,MXCOOR,1,
     &                  LUPRI)
         END IF
      END IF
C
C     ***** Symmetrize and add to HESMOL *****
C
      DO 500 I = 1, 3*NUCDEP
         DO 510 J = 1, I
            HESFS1(I,J) = HESFS1(I,J) + HESFS1(J,I)
            HESFS1(J,I) = HESFS1(I,J)
  510    CONTINUE
  500 CONTINUE
C
C     ***** PRINT SECTION *****
C
      IF (IPRINT .GT. 1) THEN
         CALL HEADER('Lowest-order reorthonormalization Hessian',-1)
         CALL PRIHES(HESFS1,'CENTERS',CSTRA,SCTRA)
         IF (SOLVNT) THEN
            CALL HEADER('Product (HSOLTT) contribution to Hessian',-1)
            CALL PRIHES(HSOLTT,'CENTERS',CSTRA,SCTRA)
         END IF
      END IF
      CALL ADDHES(HESFS1)
      IF (SOLVNT) CALL ADDHES(HSOLTT)
      IF (IPRINT .GT. 0) CALL TIMER('REORT ',TIMSTR,TIMEND)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/,A)')
     &          ' Program stopped after REORT as required.'
         CALL QUIT(' ***** End of ABACUS (in REORT) *****')
      END IF
      CALL QEXIT('REORT')
      RETURN
      END
