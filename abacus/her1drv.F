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
C  /* Deck oneinp */
      SUBROUTINE ONEINP(WORD)
C
C     Trygve Helgaker, University of Oslo, Norway
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 6)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "abainf.h"
#include "cbione.h"
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ','.NCLONE',
     *            '.NODC  ', '.NODV  ','.STOP  '/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      SET = .TRUE.
      CALL ONEINI
C
      NEWDEF = WORD .EQ. '*ONEINT'
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
                     GO TO (1,2,3,4,5,6), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ONEINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ONEINP')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3             NCLONE = .TRUE.
               GO TO 100
    4             NODC = .TRUE.
               GO TO 100
    5             NODV = .TRUE.
               GO TO 100
    6             CUT  = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ONEINP.'
            END IF
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal prompt in ONEINP')
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for ONEINT:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' ONEINT skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in ONEINT:',IPRINT
            END IF
            IF (NODC) WRITE (LUPRI,'(/,2A)') ' Inactive one-electron',
     *      ' density matrix neglected in ONEINT.'
            IF (NODV) WRITE (LUPRI,'(/,2A)') ' Active one-electron',
     *      ' density matrix neglected in ONEINT.'
            IF (NCLONE) WRITE (LUPRI,'(/,2A)') ' Only classical',
     *      ' contributions to nuclear-attraction integrals calculated.'
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after ONEINT.'
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck oneini */
      SUBROUTINE ONEINI
C
C     Initialize /CBIONE/
C
#include "implicit.h"
#include "mxcent.h"
#include <maxorb.h>

#include "abainf.h"
#include "cbione.h"
#include <infinp.h>
C
      IPRINT = IPRDEF
      SKIP   = .FALSE.
      CUT    = .FALSE.
      IF (MOLHES) THEN
         MAXDIF = 2
      ELSE IF (DIPDER .OR. MOLGRD .OR. QPGRAD) THEN
         MAXDIF = 1
      ELSE IF (EXPGRD) THEN
         MAXDIF = 0
      ELSE
         SKIP = .TRUE.
      END IF
      HFONLY = HELFEY
      NODC   = .FALSE.
      NODV   = .FALSE.
      DIFDIP = DIPDER .OR. (MOLGRD .AND. NFIELD .GT. 0)
      DIFQDP = QPGRAD
      NCLONE = .FALSE.
      RETURN
      END
C  /* Deck oneint */
      SUBROUTINE ONEINT(WORK,LWORK,PASS,PROPTY,PCM)
C
C     TUH
C
#include "implicit.h"
#include "priunit.h"
      LOGICAL PASS, PROPTY, PCM
      DIMENSION WORK(LWORK)
#include "cbione.h"
      IF (SKIP) RETURN
      IF (IPRINT .GE. 2) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (IPRINT .GT. 2) CALL TITLER('Output from ONEINT','*',103)
      PROPTY = .TRUE.
      CALL ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,NODV,
     &            DIFDIP,DIFQDP,HFONLY,NCLONE,PCM)
      IF (IPRINT .GE. 2) CALL TIMER ('ONEINT',TIMSTR,TIMEND)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/,A)')
     &          ' Program stopped after ONEINT as required.'
         WRITE (LUPRI,'(A)') ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in ONEINT) *****')
      END IF
      RETURN
      END
C  /* Deck onedrv */
      SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,
     &                  NODV,DIFDIP,DIFQDP,HFONLY,NCLONE,PCM)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "iratdef.h"
C
      LOGICAL DIFINT, NODC, NODV, DIFDIP, PROPTY, DOINT(2,2),
     &        HFONLY, NCLONE, DIFQDP, PCM
      CHARACTER*4 OMITVNUC(2)
      DIMENSION WORK(LWORK)
#include "nuclei.h"
#include "shells.h"
#include "symmet.h"
#include "cbisol.h"

#ifdef PRG_DIRAC
      LOGICAL    FINDPT
      PARAMETER (FINDPT = .FALSE.)
#else
#include <drw2el.h>
#include <inforb.h>
#endif

C
C
#include "ibtfun.h"
      CALL QENTER('ONEDRV')
C
      DOINT(1,1) = .TRUE.
      DOINT(2,1) = .TRUE.
      DOINT(1,2) = .TRUE.
      DOINT(2,2) = .TRUE.
      OMITVNUC(1)   = 'FFFF'
      OMITVNUC(2)   = 'FFFF'
C
C     ***** Number of basis functions *****
C
      NBAST  = 0
      NNBAST = 0
      DO 100 KB = 0, MAXREP
         NBASI  = NAOS(KB+1)
         NBAST  = NBAST  + NBASI
         NNBAST = NNBAST + NBASI*(NBASI + 1)/2
  100 CONTINUE
      NNBASX = NBAST*(NBAST + 1)/2
C
      KSTHMA = 1
      IF (FINDPT) THEN
         N_STHMA = 4
      ELSE
         N_STHMA = 3
      END IF
      KDENMA = KSTHMA + N_STHMA*NNBASX
      KFOCMA = KDENMA +   NNBASX
      KFACIN = KFOCMA +   NNBASX
      KCOORC = KFACIN + 2*NUCDEP
C     allocate for nuclear charge AND an electronic core charge
C     for modified nuclear attraction (e.g. for the small charge
C     in relativistic calculations). /Mar2001 hjaaj
      KSIGNC = KCOORC + 3*NUCDEP
      KNCENT = KSIGNC + 3*NUCDEP
      KJSYMC = KNCENT +  (NUCDEP + 1)/IRAT
      KJCENT = KJSYMC +  (NUCDEP + 1)/IRAT
      KGEXP  = KJCENT +  (NUCDEP + 1)/IRAT
      KCSTRA = KGEXP  +  NUCDEP
      KSCTRA = KCSTRA + MXCOOR*MXCOOR
      KTLMD  = KSCTRA + MXCOOR*MXCOOR
      IF (SOLVNT) THEN
         KHESKE = KTLMD + LMNTOT*6*NUCDEP
         IF (MAXDIF .GE. 2) THEN
            KHESNA = KHESKE + MXCOOR*MXCOOR
            KHESFS = KHESNA + MXCOOR*MXCOOR
            KHSOLT = KHESFS + MXCOOR*MXCOOR
            KHSOLN = KHSOLT + MXCOOR*MXCOOR
            KHEDFT = KHSOLN + MXCOOR*MXCOOR
            KLAST  = KHEDFT + MXCOOR*MXCOOR
         ELSE
            KLAST = KHESKE
         END IF
      ELSE
         IF (MAXDIF .GE. 2) THEN
            KHESKE = KTLMD
            KHESNA = KHESKE + MXCOOR*MXCOOR
            KHESFS = KHESNA + MXCOOR*MXCOOR
            KHEDFT = KHESFS + MXCOOR*MXCOOR
            KLAST  = KHEDFT + MXCOOR*MXCOOR
         ELSE
            KLAST  = KTLMD
            KHESKE = KLAST
            KHESNA = KLAST
            KHESFS = KLAST
            KHEDFT = KLAST
         END IF
         KHSOLT = KLAST
         KHSOLN = KLAST
      END IF
      LWRK   = LWORK  - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ONEDRV',' ',KLAST,LWORK)
      CALL ONEDR1(WORK(KSTHMA),WORK(KDENMA),WORK(KFOCMA),WORK(KFACIN),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KNCENT),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KTLMD),
     &            WORK(KCSTRA),WORK(KSCTRA),
     &            WORK(KLAST),LWRK,IPRINT,PROPTY,MAXDIF,NODC,NODV,
     &            DIFDIP,DIFQDP,NBAST,NNBASX,NNBAST,DOINT,.TRUE.,
     &            WORK(KHESKE),WORK(KHESNA),WORK(KHESFS),
     &            WORK(KHSOLT),WORK(KHSOLN),WORK(KHEDFT),
     &            OMITVNUC,HFONLY,NCLONE,PCM)
      CALL QEXIT('ONEDRV')
      RETURN
      END
C  /* Deck onedr1 */
      SUBROUTINE ONEDR1(STHMAT,DENMAT,FOCMAT,FACINT,COORC,SIGNC,NCENTC,
     &                  JSYMC,JCENTC,GEXP,TLMD,CSTRA,SCTRA,WORK,LWORK,
     &                  IPRINT,PROPTY,MAXDIF,NODC,NODV,DIFDIP,DIFQDP,
     &                  NBAST,NNBASX,NNBAST,DOINT,TOFILE,
     &                  HESSKE,HESSNA,HESFS2,HSOLT2,HSOLNN,HESDFT,
     &                  OMITVNUC,HFONLY,NCLONE,PCM)
C
C     TUH
C
C     This program calculates overlap and one-electron Hamiltonian
C     matrix elements and the first and second derivatives of these
C     elements using the McMurchie/Davidson scheme.  See L. E. McMurchie
C     & E. R. Davidson, J. Comp. Phys. 26 (1978) 218, and also V. R.
C     Saunders in "Methods in Computational Molecular Physics", Reidel
C     1983.
C
C     Symmetry included  880406  TUH & PRT
C
#include "implicit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "iratdef.h"
#include <dummy.h>
#include <drw2el.h>
#include <codata.h>
#include "priunit.h"
      PARAMETER (D0 = 0.00D00)
C
#include "nuclei.h"
      LOGICAL NODC, NODV, FRSDER, SECDER, DIFDIP, PROPTY, DOINT(2,2),
     &        TOFILE, HFONLY, NCLONE, FOUND, DIFQDP, FNDLB2, PCM,
     &        LOMITVNUC(2,5)
      CHARACTER*4 OMITVNUC(2)
      DIMENSION DENMAT(NNBASX), FOCMAT(NNBASX), STHMAT(NNBASX,*),
     &          FACINT(NUCDEP,*), WORK(LWORK), COORC(3,NUCDEP),
     &          SIGNC(3,NUCDEP), NCENTC(NUCDEP), JSYMC(NUCDEP),
     &          JCENTC(NUCDEP), GEXP(NUCDEP),TLMD(3*NUCDEP*LMNTOT,2),
     &          HESSKE(MXCOOR,MXCOOR), HESSNA(MXCOOR,MXCOOR),
     &          HESFS2(MXCOOR,MXCOOR), HSOLT2(MXCOOR,MXCOOR),
     &          HSOLNN(MXCOOR,MXCOOR), HESDFT(MXCOOR,MXCOOR),
     &          CSTRA(*), SCTRA(*)
      CHARACTER*8 RTNLBL(2)
#ifdef PRG_DIRAC
#include "dcbgen.h"
#include "dcbgrd.h"
#else
#include "gnrinf.h"
#include "energy.h"
#include <inftap.h>
#endif
#include "onecom.h"
#include "lmns.h"
#include "dipole.h"
#include "taysol.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "dorps.h"
#include "symind.h"
#include "csym1.h"
#include "cbisol.h"
#include <difsec.h>
#include <dftcom.h>
#include "expopt.h"
#include <rspprp.h>
#include <esg.h>
C
#include "ibtfun.h"
C... added by Bin Gao, Feb. 16, 2011
C... for ECP
#include "ecpinf.h"
      IF (IPRINT .GE. 5) CALL TITLER('Output from ONEDR1','*',103)
C
      TIMHER = D0
      TIMINT = D0
      FRSDER = MAXDIF .GE. 1
      SECDER = MAXDIF .EQ. 2
      IF(THRS.EQ.0D0) THRS = 1.D-15 ! to avoid log(0.0d0) below
      TOLS   = THRS*THRS
      TOLOG  = - LOG(TOLS)
      IF (PROPTY .AND. (SECDER .OR. DIFDIP .OR. DIFQDP)) THEN
         CALL GPOPEN(LUITMP,' ','NEW',' ','UNFORMATTED',IDUMMY,.FALSE.)
         REWIND LUITMP
         INDMAX = 0
         LENGTH = 0
      END IF
      DO I = 1,2
         LOMITVNUC(I,1) = (OMITVNUC(I)(1:1) .EQ. 'T')
         LOMITVNUC(I,2) = (OMITVNUC(I)(2:2) .EQ. 'T')
         LOMITVNUC(I,3) = (OMITVNUC(I)(3:3) .EQ. 'T')
         LOMITVNUC(I,4) = (OMITVNUC(I)(4:4) .EQ. 'T')
         LOMITVNUC(I,5) = LOMITVNUC(I,1) .OR. LOMITVNUC(I,2)
     &             .OR. LOMITVNUC(I,3) .OR. LOMITVNUC(I,4)
      ENDDO
C
Cjth
#ifndef PRG_DIRAC
C
C     ***** Nuclear contributions *****
C
      IF (PROPTY.AND..NOT.EXPGRA) CALL DIPNUC(CSTRA,SCTRA,IPRINT,DIFDIP)
      IF (PROPTY.AND.DIFQDP) CALL SECNUC(CSTRA,SCTRA,IPRINT)
C
#endif
C     Gradient and Hessian elements equal to zero
C
      IF (EXPGRA) THEN
         CALL DZERO(ALPGRD,KMAX)
      ELSE IF (PROPTY) THEN
         ENERKE = D0
         ENERNA = D0
#ifdef PRG_DIRAC
C
         CALL DZERO(GRADKN(1,0),3*NUCDEP)
         CALL DZERO(GRADKN(1,1),3*NUCDEP)
         CALL DZERO(GRADKN(1,2),3*NUCDEP)
C
         CALL DZERO(GRADNU(1,0),3*NUCDEP)
         CALL DZERO(GRADNU(1,1),3*NUCDEP)
         CALL DZERO(GRADNU(1,2),3*NUCDEP)
C
C
         CALL DZERO(GRADRO(1,0),3*NUCDEP)
         CALL DZERO(GRADRO(1,1),3*NUCDEP)
         CALL DZERO(GRADRO(1,2),3*NUCDEP)
C
#else 
         CALL DZERO(GRADKE,3*NUCDEP)
         CALL DZERO(GRADNA,3*NUCDEP)
         CALL DZERO(GRADFS,3*NUCDEP)
#endif
         IF (SECDER) THEN
            CALL HESZER(HESSKE)
            CALL HESZER(HESSNA)
            CALL HESZER(HESFS2)
         END IF
         IF (DIFDIP) THEN
            CALL DZERO(DDIPE,9*NUCDEP)
            CALL DZERO(DDIPS,9*NUCDEP)
         END IF
         IF (DIFQDP) THEN
            CALL DZERO(DSECE,27*NUCDEP)
            CALL DZERO(DSECS,27*NUCDEP)
         END IF
         IF (SOLVNT) THEN
            ESOLTT = D0
            ESOLNN = D0
            CALL DZERO(GSOLTT,3*NUCDEP)
            CALL DZERO(GSOLNN,3*NUCDEP)
            IF (SECDER) THEN
               CALL HESZER(HSOLT2)
               CALL HESZER(HSOLNN)
            END IF
         ELSEIF (PCM) THEN
            ESOLTT = D0
            ESOLNN = D0
            CALL DZERO(GSOLTT,3*NUCDEP)
            CALL DZERO(GSOLNN,3*NUCDEP)
            CALL DZERO(GSOLCV,3*NUCDEP)
Clf            IF (SECDER) THEN
            IF (.FALSE.) THEN
               CALL QUIT('ONEINT: '// 
     $              'SECOND DERIVATIVES NOT IMPLEMENTED WITH PCM')
            END IF
         END IF
      END IF
C
Cjth
C     The contravariant Fock matrix and the density matrix are
C     allready given from DIRAC, so skip the next.
C
#if !defined(PRG_DIRAC)
C
C     **************************************************************
C     ***** Set up total density and Fock matrices in AO basis *****
C     **************************************************************
C
      IF (PROPTY) THEN
         KDSO  = 1
         KFSO  = KDSO  + NNBAST
         KLAST = KFSO  + NNBAST
         LWRK  = LWORK - KLAST + 1
         IF (KLAST.GT.LWORK) CALL STOPIT('ONEDR1','DSOFSO',KLAST,LWORK)
         IF (ESG) THEN
            CALL GPOPEN(LUESG,'ESG_AOMAT',' ',' ',' ',IDUMMY,.FALSE.)
            REWIND(LUESG)
            READ (LUESG) MMBAST
            IF (MMBAST .NE. NNBAST)
     &      CALL QUIT('ESG error: NNBAST on LUESG not coorect')
            CALL GETESG_DENMAT_FOCMAT(WORK(KDSO),WORK(KFSO),NBAST,
     &                                NNBAST)
         ELSE
           CALL DSOFSO(WORK(KDSO),WORK(KFSO),WORK(KLAST),LWRK,IPRINT,
     &                 NODC,NODV)
         END IF
         CALL DSYM1(DENMAT,FOCMAT,WORK(KDSO),WORK(KFSO),.TRUE.,
     &              NBAST,IPRINT)
      END IF
C
C     ********************************************************
C     ***** Preparation for solvent contributions        *****
C     ***** Calculation of solvent nuclear contributions *****
C     ********************************************************
C
      KFCM  = 1
      KLAST = 1
      IF (PROPTY .AND. SOLVNT) THEN
         IF (DIFQDP) THEN
            WRITE (LUPRI,'(1X,A)') 'Second-moment gradient not yet '//
     &           'implemented for spherical cavity model'
            CALL QUIT('Second moment gradients not implemented '//
     &           'for spherical cavity model')
         END IF
         KFRSAV = 1
         KFREE  = KFRSAV
         LFREE  = LWORK
         CALL MEMGET('REAL',KFCM,LMNTOT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KGLM, LMTOT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KTLM, LMTOT,WORK,KFREE,LFREE)
         LWRK  = LWORK - KFREE + 1
         CALL SOLFL(WORK(KGLM),EPDIEL,RCAV,LCAVMX)
         CALL RD_SIRIFC('SOLTLM',FOUND,WORK(KTLM),WORK(KFREE),LWRK)
         IF (.NOT. FOUND)
     &    CALL QUIT('ONEDR1 ERROR: "SOLTLM" not found on SIRIFC.')
         CALL FCMFAC(LCAVMX,LMNTOT,LMTOT,WORK(KGLM),WORK(KTLM),
     &               WORK(KFCM),WORK(KFREE),LWRK,IPRINT)
         CALL MEMREL('ONEDR1.FCMFAC',WORK,KFRSAV,KTLM,KFREE,LFREE)
         IF (MAXDIF .GE. 1) CALL DZERO(TLMD,6*NUCDEP*LMNTOT)
         CALL SOLNUC(.TRUE.,MAXDIF,DIFDIP,TLMD(1,2),HSOLNN,
     &               WORK(KFCM),WORK,KFREE,LFREE,IPRINT)
         KLAST = KFREE
      ELSE IF(PROPTY .AND. PCM) THEN 
         CALL PCMNUC
      END IF
#endif
      LWRK  = LWORK - KLAST + 1
C
C     ************************************************************
C     ***** Triangular loop over symmetry independent shells *****
C     ************************************************************
C
      IF (FINDPT) THEN
         N_STHMA = 4
      ELSE
         N_STHMA = 3
      END IF
      CALL DZERO(STHMAT,N_STHMA*NNBASX)
      CALL IZERO(ISOFRA, 8)

#ifdef BUILD_GEN1INT
C...  added by Bin Gao, Feb. 15, 2011
      if (ECP) then
C...    computes expectation values
        if (PROPTY) then
C...      expectation values from ECP integrals
          call gen1int_driver('.ECP', 0, 0, .false., .false., .false.,
     &                        (/DUMMY/), .true., 1, DENMAT,
     &                        .true., .false., WORK(KLAST:KLAST),
     &                        LWRK-1, WORK(KLAST+1), IPRINT)
          ENERNA = ENERNA+WORK(KLAST)
C...      expectation values and writing the integrals of first-order geometric derivatives
          if (SECDER .OR. DIFDIP .OR. DIFQDP) then
            call gen1int_driver('.ECP', MAXDIF, 0, .false., .false.,
     &                          .true., (/DUMMY/), .true., 1, DENMAT,
     &                          .true., .false., GRADNA,
     &                          LWRK, WORK(KLAST), IPRINT)
C...      expectation values of first-order geometric derivatives
          else
            call gen1int_driver('.ECP', MAXDIF, 0, .false., .false.,
     &                          .false., (/DUMMY/), .true., 1, DENMAT,
     &                          .true., .false., GRADNA,
     &                          LWRK, WORK(KLAST), IPRINT)
          end if
C...      expectation values of second-order geometric derivatives
          if (SECDER) then
            call gen1int_driver('.ECP', MAXDIF, 0, .false., .false.,
     &                          .false., (/DUMMY/), .true., 1, DENMAT,
     &                          .true., .false., HESSNA,
     &                          LWRK, WORK(KLAST), IPRINT)
          end if
C...    ECP integrals, do not write the integrals (FIXME)?
        else
          call gen1int_driver('.ECP', MAXDIF, 0, .false., .true.,
     &                        .false., STHMAT(:,2), .false., 1, DENMAT,
     &                        .false., .false., (/DUMMY/),
     &                        LWRK, WORK(KLAST), IPRINT)
        end if
      end if
#endif

      IDENA = 0
      DO 100 ISHELA = 1,KMAX
         NHKTA = NHKT(ISHELA)
         KHKTA = KHKT(ISHELA)
         KCKTA = KCKT(ISHELA)
         ICA   = LCLASS(ISHELA)
         SPHRA = SPHR(ISHELA)
         CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
         NCENTA = NCENT(ISHELA)
         ICENTA = NUCNUM(NCENTA,1)
         MULA   = ISTBAO(ISHELA)
         MULTA  = MULT(MULA)
         NUCA   = NUCO(ISHELA)
         NUMCFA = NUMCF(ISHELA)
         JSTA   = JSTRT(ISHELA)
         CORAX  = CENT(ISHELA,1,1)
         CORAY  = CENT(ISHELA,2,1)
         CORAZ  = CENT(ISHELA,3,1)
         IDENB0 = 0
C
C        Compute symmetry integral pointers for contributions
C        from this block.  Note that at present this assumes all
C        components from a shell are included.
C
         DO 600 I = 1, 8
            ISOFRB(I) = 0
            DO 610 J = 1, MXAQN
              INDFA(I,J) = -10 000 000
610         CONTINUE
600      CONTINUE
         DO 620 NA = 1, KHKTA
            DO 630 IREP = 0, MAXREP
            IF (IBTAND(MULA,IBTXOR(IREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
               ISOFRA(IREP+1)    = ISOFRA(IREP+1) + 1
               INDFA (IREP+1,NA) = ISOFRA(IREP+1)
            END IF
630         CONTINUE
620      CONTINUE
         IF (IPRINT .GT. 20) THEN
            WRITE(LUPRI,'(A,I4)')' IA address offsets for shell ',ISHELA
            DO 640 NA = 1,KHKTA
               WRITE(LUPRI,'(8(1X,I5))') (INDFA(I,NA), I = 1,MAXREP+1)
640         CONTINUE
         END IF
      DO 110 ISHELB = 1,ISHELA
         LDIAG = ISHELA .EQ. ISHELB
         NHKTB = NHKT(ISHELB)
         KHKTB = KHKT(ISHELB)
         KCKTB = KCKT(ISHELB)
         ICB   = LCLASS(ISHELB)
         SPHRB = SPHR(ISHELB)
         CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
         NCENTB = NCENT(ISHELB)
         NHKTAB = NHKTA + NHKTB
         MULB   = ISTBAO(ISHELB)
         MULTB  = MULT(MULB)
         NUCB   = NUCO(ISHELB)
         NUMCFB = NUMCF(ISHELB)
         JSTB   = JSTRT(ISHELB)
         CORBX0 = CENT(ISHELB,1,1)
         CORBY0 = CENT(ISHELB,2,1)
         CORBZ0 = CENT(ISHELB,3,1)
         KHKTAB = KHKTA*KHKTB
         KCKTAB = KCKTA*KCKTB
         MAB    = IBTOR(MULA,MULB)
         KAB    = IBTAND(MULA,MULB)
         HKAB   = FMULT(KAB)
C
         SPHRAB = SPHRA .OR. SPHRB
C
C        Compute symmetry integral pointers for contributions
C        from this block.  Note that at present this assumes all
C        components from a shell are included
C
         DO 700 I = 1, 8
            DO 710 J = 1, MXAQN
              INDFB(I,J) = -10 000 000
710         CONTINUE
700      CONTINUE
         DO 720 NB = 1, KHKTB
            DO 730 IREP = 0, MAXREP
            IF (IBTAND(MULB,IBTXOR(IREP,ISYMAO(NHKTB,NB))).EQ.0) THEN
               ISOFRB(IREP+1)    = ISOFRB(IREP+1) + 1
               INDFB (IREP+1,NB) = ISOFRB(IREP+1)
            END IF
730         CONTINUE
720      CONTINUE
         IF (IPRINT .GT. 20) THEN
            WRITE(LUPRI,'(A,I4)')' IB address offsets for shell ',ISHELB
            DO 740 NB = 1, KHKTB
               WRITE(LUPRI,'(8(1X,I5))') (INDFB(I,NB), I = 1,MAXREP+1)
740         CONTINUE
         ENDIF
         IF(.NOT.DOINT(ICA,ICB)) GOTO 110
         IF (IPRINT .GE. 05) WRITE (LUPRI, 1000) ISHELA, ISHELB
         IF (IPRINT .GE. 10) THEN
             WRITE (LUPRI,'(A,2I10)') ' NHKT   ', NHKTA, NHKTB
             WRITE (LUPRI,'(A,2I10)') ' KHKT   ', KHKTA, KHKTB
             WRITE (LUPRI,'(A,2I10)') ' KCKT   ', KCKTA, KCKTB
             WRITE (LUPRI,'(A,2I10)') ' NCENT  ', NCENTA, NCENTB
             WRITE (LUPRI,'(A,2I10)') ' ISTBAO ', MULA, MULB
             WRITE (LUPRI,'(A,2I10)') ' MULT   ', MULTA, MULTB
             WRITE (LUPRI,'(A,2I10)') ' NUC    ', NUCA, NUCB
             WRITE (LUPRI,'(A,2I10)') ' NUMCF  ', NUMCFA, NUMCFB
             WRITE (LUPRI,'(A,2I10)') ' JST    ', JSTA, JSTB
             WRITE (LUPRI,'(A,2F12.6)') ' CORAX    ', CORAX, CORBX0
             WRITE (LUPRI,'(A,2F12.6)') ' CORAY    ', CORAY, CORBY0
             WRITE (LUPRI,'(A,2F12.6)') ' CORAZ    ', CORAZ, CORBZ0
         END IF
C
C        Initialization for nuclear attraction integrals
C
         JMAX = NHKTAB - 2
         IF (EXPGRA) THEN
            JMAX = JMAX + 2
         ELSE IF (PROPTY) THEN
            JMAX = JMAX + MAXDIF
         END IF
         ISTEPU = JMAX + 1
         ISTEPV = ISTEPU*ISTEPU
         NAHGTF = ISTEPU*ISTEPV
         NATOMC = 0
         IMXATM = NUCIND
         IF (QM3) IMXATM = NCTOT
         DO 120 IATOMC = 1,IMXATM
            MULC   = ISTBNU(IATOMC)
            MABC   = IBTOR(MULC,KAB)
            CORCX0 = CORD(1,IATOMC)
            CORCY0 = CORD(2,IATOMC)
            CORCZ0 = CORD(3,IATOMC)
            CHARG1 = CHARGE(IATOMC)
            FACTOR = - FMULT(IBTAND(MULC,KAB))/HKAB
            DO 130 ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,MABC) .EQ. 0) THEN
                  NATOMC = NATOMC + 1
                  JSYMC(NATOMC)   = ISYMOP
                  JCENTC(NATOMC)  = IATOMC
                  SIGNC(1,NATOMC) = PT(IBTAND(ISYMAX(1,1),ISYMOP))
                  SIGNC(2,NATOMC) = PT(IBTAND(ISYMAX(2,1),ISYMOP))
                  SIGNC(3,NATOMC) = PT(IBTAND(ISYMAX(3,1),ISYMOP))
                  COORC(1,NATOMC) = SIGNC(1,NATOMC)*CORCX0
                  COORC(2,NATOMC) = SIGNC(2,NATOMC)*CORCY0
                  COORC(3,NATOMC) = SIGNC(3,NATOMC)*CORCZ0
                  GEXP(NATOMC)    = GNUEXP(IATOMC)
                  FACINT(NATOMC,1)= FACTOR*CHARG1
                  FACINT(NATOMC,2)= D0
                  NCENTC(NATOMC)  = NUCNUM(IATOMC,ISYMOP+1)
               END IF
  130       CONTINUE
  120    CONTINUE
C         IF (NATOMC.NE.NUCDEP)
C     &      print *,'ONEDR1 NATOMC, NUCDEP:',NATOMC,NUCDEP
C
         CALL ONESOP(STHMAT,DENMAT,FOCMAT,FACINT,COORC,
     &               WORK(KLAST),LWRK,
     &               IPRINT,PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,CORBZ0,
     &               DIFDIP,DIFQDP,SECDER,NATOMC,TOLOG,TOLS,JSYMC,
     &               JCENTC,NCENTC,SIGNC,GEXP,NNBASX,WORK(KFCM),TLMD,
     &               HESSKE,HESSNA,HESFS2,HSOLT2,HSOLNN,NCLONE,PCM)
  110    IDENB0 = IDENB0 + KHKTB*MULTB
         IDENA  = IDENA  + KHKTA*MULTA
  100 CONTINUE
C
C     ***** End loop over symmetry independent orbitals *****
C
C     ****************************************************
C     ***** Write final buffers and dipole integrals *****
C     ***** ELSE write undifferentiated integrals    *****
C     ****************************************************
C
      IF (PROPTY) THEN
         IF (SECDER .OR. DIFDIP .OR. DIFQDP) THEN
C
C           Write final buffer of first derivative integrals
C
            IF (LENGTH .GT. 0) WRITE (LUITMP) BUF, IBUF, LENGTH
            WRITE (LUITMP) BUF, IBUF, -LENGTH
            IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     *       ' Last buffer of length',LENGTH,
     *       ' has been written on LUITMP in ONEDR1.'
         END IF
      ELSE
         IF (FINDPT) THEN
C
C           Add DPT correction to nuclear attraction integrals (WK/UniKA/26-03-2004).
C
            IF (LUPROP .GT. 0) CALL GPCLOSE(LUPROP,'KEEP')
            CALL GPOPEN(LUPROP,'AOPROPER',
     &                  'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
            IF (FNDLB2('DERXXPVP',RTNLBL,LUPROP)) THEN
                CALL READT(LUPROP,NNBASX,STHMAT(1,4))
            ELSE
                CALL QUIT('Integrals with label ''DERXXPVP'' '//
     &                    'not found on file ''AOPROPER''.')
            END IF 
            FF = DPTFAC * ALPHA2
            CALL PKSYM1(STHMAT(1,4),STHMAT(1,4),NAOS,MAXREP+1,2)
            CALL DAXPY(NNBASX,FF,STHMAT(1,4),1,STHMAT(1,2),1)
         END IF 
C
C        Write undifferentiated one-electron integrals
C
        IF(TOFILE) THEN
C
C     We use the modified Douglas-Kroll integrals instead
C
           IF (DKTRAN) THEN
              CALL DKTRF(STHMAT(1,2),NBAST,NNBAST,IPRINT)
           ELSE
              DO 800 I = 1, NNBAST
                 STHMAT(I,2) = STHMAT(I,2) + STHMAT(I,3)
 800          CONTINUE
           END IF
           CALL WRTUND(STHMAT,NNBAST,NNBASX,IPRINT)
        ENDIF
      END IF
C
C     *******************************
C     ***** Symmetrize Hessians *****
C     *******************************
C
      IF (PROPTY .AND. SECDER) THEN
         DO 500 I = 2,3*NUCDEP
            DO 510 J = 1,I-1
               HESSKE(I,J) = HESSKE(I,J) + HESSKE(J,I)
               HESSKE(J,I) = HESSKE(I,J)
               HESSNA(I,J) = HESSNA(I,J) + HESSNA(J,I)
               HESSNA(J,I) = HESSNA(I,J)
               HESFS2(I,J) = HESFS2(I,J) + HESFS2(J,I)
               HESFS2(J,I) = HESFS2(I,J)
               IF (SOLVNT) THEN
                  HSOLT2(I,J) = HSOLT2(I,J) + HSOLT2(J,I)
                  HSOLT2(J,I) = HSOLT2(I,J)
               END IF
  510       CONTINUE
  500    CONTINUE
      END IF
C
C     ********************************************************
C     ***** Sort gradient elements on direct access unit *****
C     ********************************************************
C
      IF (PROPTY) THEN
         IF (SECDER .OR. DIFDIP .OR. DIFQDP) THEN
            IF (SOLVNT) THEN
               CALL TLMWRT(TLMD,WORK(KFCM),WORK(KLAST),LWRK,IPRINT)
            END IF
            NMATS = 3*NUCDEP*(MAXREP+1)
            CALL SORONE(WORK,LWORK,NMATS,INDMAX,IPRINT)
            IF (IPRINT .GE. 5) CALL SHDPRI(WORK,LWORK)
            IF (DOREPS(0)) CALL SHDCHK(WORK,LWORK,NODC,NODV,IPRINT)
            CALL GPCLOSE(LUITMP,' ')
         END IF
      END IF
C
C     ****************************************
C     ***** DFT contribution to gradient *****
C     ****************************************
C
      IF (PROPTY .AND. DFTRUN) THEN
         IF (EXPGRA) THEN
            CALL DFTEXPGRAD(WORK,LWORK,1)
         ELSE
            CALL DFTMOLGRAD(WORK,LWORK,1)
         END IF
C        CALL HEADER('GRADFT in DFT_MOLGRAD ',-1)
C        CALL OUTPUT(GRADFT,1,3,1,NATOMS,3,NATOMS,1,LUPRI)
      END IF
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IF (IPRINT .GT. 0) THEN
         IF (EXPGRA) THEN
            CALL HEADER('One-electron part of alpha gradient',-1)
            WRITE (LUPRI,'(2X,5F12.8)') (ALPGRD(I),I=1,KMAX)
         ELSE IF (PROPTY) THEN
            KCSTRA = 1
            KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
            KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
            IF (KLAST .GT. LWORK) 
     &           CALL STOPIT('ONEDR1','TRACOR',KLAST,LWORK)
            CALL HEADER('Kinetic energy integral gradient',-1)
            CALL PRIGRD(GRADKE,WORK(KCSTRA),WORK(KSCTRA))
            CALL HEADER('Nuclear attraction integral gradient',-1)
            CALL PRIGRD(GRADNA,WORK(KCSTRA),WORK(KSCTRA))
            CALL HEADER('Reorthonormalization gradient',-1)
            CALL PRIGRD(GRADFS,WORK(KCSTRA),WORK(KSCTRA))
            IF (SOLVNT) THEN
               CALL HEADER('Solvent energy contributions',-1)
               WRITE(LUPRI,'(3(/A,F15.10))')
     &            ' Nuclear solvent energy:',ESOLNN,
     &            ' Electr. solvent energy:',ESOLTT,
     &            ' Total   solvent energy:',ESOLNN+ESOLTT
               CALL HEADER('Electronic solvent part of gradient',-1)
               CALL PRIGRD(GSOLTT,WORK(KCSTRA),WORK(KSCTRA))
               CALL HEADER('Nuclear solvent part of gradient',-1)
               CALL PRIGRD(GSOLNN,WORK(KCSTRA),WORK(KSCTRA))
            ELSE IF (PCM) THEN
               CALL HEADER('PCM Solvent energy contributions',-1)
               WRITE(LUPRI,'(3(/A,F15.10))')
     &              ' Nuclear solvent energy:',ESOLNN,
     &              ' Electr. solvent energy:',ESOLTT,
     &              ' Total   solvent energy:',ESOLNN+ESOLTT
               CALL HEADER('Electronic solvent part of gradient',-1)
               CALL PRIGRD(GSOLTT,WORK(KCSTRA),WORK(KSCTRA))
               CALL HEADER('Nuclear solvent part of gradient',-1)
               CALL PRIGRD(GSOLNN,WORK(KCSTRA),WORK(KSCTRA))
               CALL HEADER('Cavity solvent part of the gradient',-1)
               CALL PRIGRD(GSOLCV,WORK(KCSTRA),WORK(KSCTRA))
            END IF
            IF (DFTRUN) THEN
               CALL HEADER('DFT contribution to gradient',-1)
               CALL PRIGRD(GRADFT,WORK(KCSTRA),WORK(KSCTRA))
            END IF
            IF (SECDER) THEN
               CALL HEADER('Kinetic energy integral Hessian',-1)
               CALL PRIHES(HESSKE,'CENTERS',WORK(KCSTRA),WORK(KSCTRA))
               CALL HEADER('Nuclear attraction integral Hessian',-1)
               CALL PRIHES(HESSNA,'CENTERS',WORK(KCSTRA),WORK(KSCTRA))
               CALL HEADER('Highest order reorthonormalization Hessian',
     &                      -1)
               CALL PRIHES(HESFS2,'CENTERS',WORK(KCSTRA),WORK(KSCTRA))
               IF (SOLVNT) THEN
                 CALL HEADER('HSOLT2 part of solvent Hessian',-1)
                 CALL PRIHES(HSOLT2,'CENTERS',WORK(KCSTRA),WORK(KSCTRA))
                 CALL HEADER('HSOLNN part of solvent Hessian',-1)
                 CALL PRIHES(HSOLNN,'CENTERS',WORK(KCSTRA),WORK(KSCTRA))
               END IF
            END IF
            IF (DIFDIP) THEN
               CALL HEADER('Electronic contributions to static part'//
     &                     ' of dipole gradient',-1)
               CALL FCPRI(DDIPE,'APT',WORK(KCSTRA),WORK(KSCTRA))
            END IF
            IF (DIFQDP) THEN
               CALL HEADER('Electronic contributions to static part'//
     &                     ' of second moment gradient',-1)
               CALL PRISEC(DSECE,'SECDER',WORK(KCSTRA),WORK(KSCTRA))
            END IF
         END IF
      END IF
C
C     Add Hessian contributions to total Hessian
C
      IF (PROPTY .AND. SECDER) THEN
         CALL ADDHES(HESSKE)
         CALL ADDHES(HESSNA)
         CALL ADDHES(HESFS2)
         IF (SOLVNT) THEN
            CALL ADDHES(HSOLT2)
            CALL ADDHES(HSOLNN)
         END IF
      END IF
      RETURN
 1000 FORMAT (//,2X,'***************************************',
     *         /,2X,'********** ISHELA/B =',I3,',',I3,' **********',
     *         /,2X,'***************************************',/)
      END
C  /* Deck dsofso */
      SUBROUTINE DSOFSO(DSO,FSO,WORK,LWORK,IPRINT,NODC,NODV)
C
C     This subroutine calculates the folded total one-electron density
C     and Fock matrices in SO basis (contravariant).  Input is
C     one-electron active density matrix and total Fock matrix in
C     MO basis.
C
C     Symmetry included        880418  PRT & TUH
C
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
! abainf.h for ABASOP
#include "mxcent.h"
#include "abainf.h"
      LOGICAL NODC, NODV
      DIMENSION DSO(NNBAST), FSO(NNBAST), WORK(LWORK)
C
      IF (ABASOP) THEN
         LKDV  = NORBT*NORBT
         LKUDV = NORBT*NORBT
      ELSE
         LKDV  = NNASHX
         LKUDV = 0
      ENDIF
      KCMO  = 1
      KFT   = KCMO + NCMOT
      KDV   = KFT  + N2ORBT
      KUDV  = KDV  + LKDV
      KLAST = KUDV + LKUDV
      IF (KLAST .GT. LWORK) CALL STOPIT('DSOFSO',' ',KLAST,LWORK)
      CALL DSOFSO_1(WORK(KCMO),WORK(KFT),WORK(KDV),WORK(KUDV),
     &            DSO,FSO,IPRINT,NODC,NODV)
      RETURN
      END
C  /* Deck dsofso_1 */
      SUBROUTINE DSOFSO_1(CMO,FT,DV,UDV,DSO,FSO,IPRINT,NODC,NODV)
C
C     This subroutine calculates the folded total one-electron density
C     and Fock matrices in AO basis (contravariant).  Input is
C     one-electron active density matrix and total Fock matrix in
C     MO basis.
C
C     Symmetry included        880418  PRT & TUH
C
#include "implicit.h"
      PARAMETER (D0 = 0.0D0, TWO = 2.0D0, HALF = 0.5D0)
      LOGICAL NODC, NODV
      DIMENSION CMO(*), FT(*), DV(*), DSO(NNBAST), FSO(NNBAST)
      DIMENSION UDV(NORBT,NORBT)
#include "iratdef.h"
#include "mxcent.h"
#include "priunit.h"
#include "maxorb.h"
C
C Used from common blocks:
C  INFRSP : CCPPA
C
#include <abainf.h>
#include <inforb.h>
#include <inftap.h>
#include <infrsp.h>
      INTEGER R, S, RS, U, V, UV
C
C
C     ***** Read input from LUSIFC *****
C
      REWIND LUSIFC
      CALL MOLLAB('SIR IPH ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,
     *              NCDETS,NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT
      NSSHT  = NORBT - NOCCT
      CALL READI(LUSIFC,IRAT*NCMOT,CMO)
      READ (LUSIFC)
      IF (NASHT .GT. 0) THEN
         CALL READI(LUSIFC,IRAT*NNASHX,DV)
      ELSE
         READ (LUSIFC)
      END IF
      CALL READI(LUSIFC,IRAT*N2ORBT,FT)
C     
      IF (ABASOP) THEN
         IF (CCPPA) THEN
            CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
         ELSE
            CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
         ENDIF
         READ (LUSIFC)
         CALL READT (LUSIFC,NORBT*NORBT,UDV)
         CALL DSITSP(NORBT,UDV,DV)
         CALL PKSYM1(DV,DV,NORB,NSYM,+1)
      END IF
C
C     ***** Print Section *****
C
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, '(//A/)') ' <<<<< SUBROUTINE DSOFSO_1 >>>>>'
         WRITE (LUPRI, '(A,8I5)') ' NISH ', (NISH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NASH ', (NASH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NOCC ', (NOCC(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NORB ', (NORB(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NBAS ', (NBAS(I),I = 1,NSYM)
         IF (IPRINT .GE. 10) THEN
            CALL HEADER('Occupied molecular orbitals',0)
            IEND = 0
            DO 1000 ISYM = 1,NSYM
               IF (NBAS(ISYM) .EQ. 0) GOTO 1000
               IF (NOCC(ISYM) .EQ. 0) GOTO 1100
               WRITE (LUPRI, '(//,A,I5,/)') ' Symmetry ', ISYM
               IENDI = IEND
               DO 1200 I = 1, NOCC(ISYM)
                  WRITE (LUPRI,'(/,A,I5,/)') ' Molecular orbital ', I
                  WRITE (LUPRI,'(6F12.6)') (CMO(IENDI+J),J=1,NBAS(ISYM))
                  IENDI = IENDI + NBAS(ISYM)
1200           CONTINUE
               CALL HEADER('Total Fock matrix (MO basis)',-1)
               CALL OUTPUT(FT(I2ORB(ISYM)+1),1,NORB(ISYM),1,NORB(ISYM),
     &                     NORB(ISYM),NORB(ISYM),1,LUPRI)
1100           CONTINUE
               IEND = IEND + NORB(ISYM)*NBAS(ISYM)
1000        CONTINUE
            CALL HEADER('Active density matrix (MO basis)',-1)
            CALL OUTPAK(DV,NASHT,1,LUPRI)
         END IF
      END IF
C
C     ***** Construct contravariant SO matrices *****
C
      ISEND = 0
      ICEND = 0
      IFEND = 0
      DO 110 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
         NISHI = NISH(ISYM)
         NASHI = NASH(ISYM)
         IASHI = IASH(ISYM)
         NBASI = NBAS(ISYM)
         IF (NBASI .EQ. 0) GOTO 120
         IF (.NOT. ABASOP) THEN
            IF (NOCC(ISYM) .EQ. 0) THEN
               CALL DZERO(DSO(ISEND+1),NNBAS(ISYM))
               CALL DZERO(FSO(ISEND+1),NNBAS(ISYM))
               GO TO 120
            END IF
         END IF
         RS = 0
C
C        ***** for SOPPA *****
C
         IF (ABASOP) THEN
C
C           ***************************************
C           ***** One-electron density matrix *****
C           ***************************************
C
            CALL UHUT(DV(IIORB(ISYM)+1), DSO(IIBAS(ISYM)+1),
     *            CMO(ICMO(ISYM)+1), NBAS(ISYM), NORB(ISYM))
            CALL DSCAL(NNBAS(ISYM),TWO,DSO(IIBAS(ISYM)+1),1)
            DO I = 1,NBAS(ISYM)
               II = I*(I+1)/2
               DSO(IIBAS(ISYM)+II) = DSO(IIBAS(ISYM)+II)*HALF
            END DO
C
C           ***********************
C           ***** Fock matrix *****
C           ***********************
C
            DO 10 R = 1, NBASI
               DO 20 S = 1,R
                  RS = RS + 1

                  FTRS = D0
                  IJ = 0
                  ICENDI = 0
                  DO 50 I = 1, NORBI
                     ICENDJ = 0
                     DO 51 J = 1, NOCC(ISYM)
                        IJ = IJ + 1
                        FTIJ = FT(IFEND+IJ)
                        IF (ABS(FTIJ) .GT. D0) THEN
                           TEMP = CMO(ICEND+ICENDI+R)
     &                           *CMO(ICEND+ICENDJ+S)
                           IF (R .NE. S) TEMP = TEMP +
     *                         CMO(ICEND+ICENDI+S)*CMO(ICEND+ICENDJ+R)
                           FTRS = FTRS + FTIJ*TEMP
                        END IF
                        ICENDJ = ICENDJ + NBASI
   51                CONTINUE
                     IJ = IJ + NSSH(ISYM)
                     ICENDI = ICENDI + NBASI
   50             CONTINUE
                  FSO(ISEND+RS) = FTRS
   20          CONTINUE
   10       CONTINUE
C
         ELSE
C
C        ***** for all other methods *****
C
            DO 100 R = 1, NBASI
               DO 200 S = 1,R
                  RS = RS + 1
C
C              ***************************************
C              ***** One-electron density matrix *****
C              ***************************************
C
                  DTRS = D0
C
C              (I) Inactive contribution
C
                  IF (NISHI .GT. 0) THEN
                     ICENDI = ICEND
                     DO 300 I = 1, NISHI
                        DTRS = DTRS + CMO(ICENDI+R)*CMO(ICENDI+S)
                        ICENDI = ICENDI + NBASI
  300                CONTINUE
                     DTRS = DTRS + DTRS
                  END IF
                  IF (NODC) DTRS = D0
C
C              (II) Active contribution
C
                  IF (.NOT. NODV) THEN
                     IF (NASHI .GT. 0) THEN
                        UV = ((IASHI + 1)*(IASHI + 2))/2
                        IDVEND = ICEND + NISHI*NBASI
                        ICENDU = IDVEND
                        DO 400 U = 1,NASHI
                           ICENDV = IDVEND
                           DO 410 V = 1, U
                              DUV = DV(UV)
                              IF (ABS(DUV) .GT. D0) THEN
                                 TEMP = CMO(ICENDU+R)*CMO(ICENDV+S)
                                 IF (U .NE. V) TEMP = TEMP
     *                               + CMO(ICENDU+S)*CMO(ICENDV+R)
                                 DTRS = DTRS + DUV*TEMP
                              END IF
                              UV = UV + 1
                              ICENDV = ICENDV + NBASI
  410                      CONTINUE
                           UV = UV + IASHI
                           ICENDU = ICENDU + NBASI
  400                   CONTINUE
                     END IF
                  END IF
                  IF (R .NE. S) DTRS = DTRS + DTRS
                  DSO(ISEND+RS) = DTRS
C
C              ***********************
C              ***** Fock matrix *****
C              ***********************
C
                  FTRS = D0
                  IJ = 0
                  ICENDI = 0
                  DO 500 I = 1, NORBI
                     ICENDJ = 0
                     DO 510 J = 1, NOCC(ISYM)
                        IJ = IJ + 1
                        FTIJ = FT(IFEND+IJ)
                        IF (ABS(FTIJ) .GT. D0) THEN
                           TEMP = CMO(ICEND+ICENDI+R)
     &                           *CMO(ICEND+ICENDJ+S)
                           IF (R .NE. S) TEMP = TEMP +
     *                         CMO(ICEND+ICENDI+S)*CMO(ICEND+ICENDJ+R)
                           FTRS = FTRS + FTIJ*TEMP
                        END IF
                        ICENDJ = ICENDJ + NBASI
  510                CONTINUE
                     IJ = IJ + NSSH(ISYM)
                     ICENDI = ICENDI + NBASI
  500             CONTINUE
                  FSO(ISEND+RS) = FTRS
  200          CONTINUE
  100       CONTINUE
         ENDIF
C
C        ***** Print Section *****
C
         IF (IPRINT .GE. 10) THEN
            WRITE (LUPRI,'(1X,A,I5)') ' Symmetry', ISYM
            CALL HEADER('Total density matrix (SO basis)',-1)
            CALL OUTPAK(DSO(ISEND+1),NBASI,1,LUPRI)
            CALL HEADER('Total Fock matrix (SO basis)',-1)
            CALL OUTPAK(FSO(ISEND+1),NBASI,1,LUPRI)
         END IF
120      CONTINUE
         IFEND = IFEND + NORBI*NORBI
         ISEND = ISEND + (NBASI*(NBASI + 1))/2
         ICEND = ICEND + NORBI*NBASI
110   CONTINUE
      RETURN
      END
C  /* Deck dipnuc */
      SUBROUTINE DIPNUC(CSTRA,SCTRA,IPRINT,DIFDIP)
C
C     Calculates nuclear contributions to electric dipole moments
C     and dipole gradients
C
C     1985 tuh
C     symmetry Dec 1988 tuh
C
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (D0 = 0.0D0, D100 = 100.0D0)
      LOGICAL DIFDIP
      DIMENSION CSTRA(*), SCTRA(*)
#include <orgcom.h>
#include <abainf.h>
#include <nuclei.h>
#include <dipole.h>
#include <symmet.h>
#include <dorps.h>
#include <ibtfun.h>
C
C     Dipole moment
C
      CALL DZERO(DIPMN,3)
C
      DO 100 IATOM = 1, NUCIND
         CHA = CHARGE(IATOM)
         IF (ABS(CHA) .LT. D100) THEN
            FAC = FMULT(ISTBNU(IATOM))*CHA
            DO 200 ICOOR = 1, 3
               IF (ISYMAX(ICOOR,1) .EQ. 0) THEN
                  DIPMN(ICOOR) = DIPMN(ICOOR) 
     &                         + FAC*(CORD(ICOOR,IATOM) - DIPORG(ICOOR))
               END IF
 200        CONTINUE
         END IF
 100  CONTINUE
C
C     Dipole moment gradient
C
      IF (DIFDIP) THEN
         CALL DZERO(DDIPN,9*NUCDEP)
         DO 300 IREP = 0, MAXREP
            IF (DOREPS(IREP) .AND. (NAXREP(IREP,1).GT.0)) THEN
               DO 400 IATOM = 1, NUCIND
                  CHA = CHARGE(IATOM)
                  IF (ABS(CHA) .LT. D100) THEN
                     FAC = FMULT(ISTBNU(IATOM))*CHA
                     DO 500 ICOOR = 1, 3
                        IF (ISYMAX(ICOOR,1) .EQ. IREP) THEN
                           ISC = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                           IF (ISC.GT.0) DDIPN(IPTAX(ICOOR,1),ISC) = FAC
                        END IF
  500                CONTINUE
                  END IF
  400          CONTINUE
            END IF
  300    CONTINUE
      END IF
      IF (IPRINT .GT. 0) THEN
         CALL HEADER('Nuclear contribution to dipole moments',-1)
         CALL DP0PRI(DIPMN)
         IF (DIFDIP) THEN
            CALL HEADER('Nuclear contribution to dipole gradient',-1)
            CALL FCPRI(DDIPN,'APT',CSTRA,SCTRA)
         END IF
      END IF
      RETURN
      END
C  /* Deck secnuc */
      SUBROUTINE SECNUC(CSTRA,SCTRA,IPRINT)
C
C     Calculates nuclear contributions to electric second moment gradient
C
C     2003 hs
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.0D0, D100 = 100.0D0)
      DIMENSION CSTRA(*), SCTRA(*)
#include <orgcom.h>
#include <abainf.h>
#include <nuclei.h>
#include <difsec.h>
#include <symmet.h>
#include <dorps.h>
#include <ibtfun.h>
C
C     Second moment gradient
C
      CALL DZERO(DSECN,27*NUCDEP)
      DO IATOM = 1, NUCIND
         CHA = CHARGE(IATOM)
         IF (ABS(CHA) .LT. D100) THEN
            FAC = FMULT(ISTBNU(IATOM))*CHA
            DO IXCOOR = 1, 3
               DO IYCOOR = 1, IXCOOR
                  IREP = IBTXOR(ISYMAX(IXCOOR,1),ISYMAX(IYCOOR,1))
                  IF (IXCOOR .EQ. IYCOOR) THEN
                     ISC = IPTCNT(3*(IATOM - 1) + IXCOOR,IREP,1)
                     IF (ISC .GT. 0) THEN
                        DSECN(IPTAX(IXCOOR,1),IPTAX(IXCOOR,1),ISC) = 
     &                       2*FAC*(CORD(IXCOOR,IATOM))
                     END IF
                  ELSE
                     ISC = IPTCNT(3*(IATOM - 1) + IXCOOR,IREP,1)
                     IF (ISC .GT. 0) THEN
                        DSECN(IPTAX(IXCOOR,1),IPTAX(IYCOOR,1),ISC) = 
     &                       FAC*(CORD(IYCOOR,IATOM))
                        DSECN(IPTAX(IYCOOR,1),IPTAX(IXCOOR,1),ISC) = 
     &                       DSECN(IPTAX(IXCOOR,1),IPTAX(IYCOOR,1),ISC) 
                     END IF
                     ISC = IPTCNT(3*(IATOM - 1) + IYCOOR,IREP,1)
                     IF (ISC .GT. 0) THEN
                        DSECN(IPTAX(IXCOOR,1),IPTAX(IYCOOR,1),ISC) = 
     &                       FAC*(CORD(IXCOOR,IATOM))
                        DSECN(IPTAX(IYCOOR,1),IPTAX(IXCOOR,1),ISC) = 
     &                       DSECN(IPTAX(IXCOOR,1),IPTAX(IYCOOR,1),ISC) 
                     END IF
                  END IF
               END DO
            END DO
         END IF
      END DO
      IF (IPRINT .GT. 0) THEN
         CALL HEADER
     &        ('Nuclear contribution to second moment gradient',-1)
         CALL PRISEC(DSECN,'SECDER',CSTRA,SCTRA)
      END IF
      RETURN
      END
C  /* Deck tlmwrt */
      SUBROUTINE TLMWRT(TLMD,FCM,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      DIMENSION TLMD(*), FCM(*), WORK(LWORK)
#include "nuclei.h"
      KGRADE = 1
      KGRADN = KGRADE + 3*NUCDEP
      KGDIFE = KGRADN + 3*NUCDEP
      KGDIFN = KGDIFE + 3*NUCDEP
      KLAST  = KGDIFN + 3*NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('TLMWRT',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL TLMWR1(TLMD,FCM,WORK(KGRADE),WORK(KGRADN),WORK(KGDIFE),
     &            WORK(KGDIFN),WORK(KLAST),LWRK,IPRINT)
      RETURN
      END
C  /* Deck tlmwr1 */
      SUBROUTINE TLMWR1(TLMD,FCM,GRADE,GRADN,GDIFE,GDIFN,WORK,LWORK,
     &                  IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "iratdef.h"
      PARAMETER (DM1 = -1.0D00, THRSH = 1.0D-07)
      LOGICAL ERRORE, ERRORN, OLDDX
#include "nuclei.h"
      DIMENSION TLMD(LMNTOT,3*NUCDEP,2), FCM(LMNTOT), GRADE(*),
     &          GRADN(*), GDIFE(*), GDIFN(*), WORK(LWORK)
#include "inftap.h"
#include "cbisol.h"
#include "taysol.h"
C
C     Test TLMD
C
      IF (IPRINT .GT. 6) THEN
         CALL HEADER('Electronic TLMED in Cart. basis',-1)
         WRITE(LUPRI,'(A)') ' Rows and columns: lmn and ICOOR'
         CALL OUTPUT(TLMD(1,1,1),1,LMNTOT,1,3*NUCDEP,
     &               LMNTOT,3*NUCDEP,1,LUPRI)
         CALL HEADER('Nuclear TLMND in Cart. basis',-1)
         WRITE(LUPRI,'(A)') ' Rows and columns: lmn and ICOOR'
         CALL OUTPUT(TLMD(1,1,2),1,LMNTOT,1,3*NUCDEP,
     &               LMNTOT,3*NUCDEP,1,LUPRI)
      END IF
      DO 100 I = 1, 3*NUCDEP
         GRADE(I) =   DDOT(LMNTOT,FCM,1,TLMD(1,I,1),1)
         GRADN(I) = - DDOT(LMNTOT,FCM,1,TLMD(1,I,2),1)
         GDIFE(I) = GSOLTT(I) - GRADE(I)
         GDIFN(I) = GSOLNN(I) - GRADN(I)
  100 CONTINUE
      ERRORE = DNORM2(3*NUCDEP,GDIFE,1) .GT. THRSH
      ERRORN = DNORM2(3*NUCDEP,GDIFN,1) .GT. THRSH
      IF (IPRINT .GT. 2 .OR. ERRORE) THEN
         IF (ERRORE) WRITE (LUPRI,'(A)')
     &     ' Error in electronic solvent part of gradient '
         CALL HEADER('Difference GSOLTT and GRADE in TLMSTS',-1)
         WRITE (LUPRI,'(11X,A)')
     &      'GSOLTT             GRADE               difference'
         WRITE (LUPRI,'((2X,3F20.10))')
     &         (GSOLTT(I), GRADE(I), GDIFE(I),I=1,3*NUCDEP)
      END IF
      IF (IPRINT .GT. 2 .OR. ERRORN) THEN
         IF (ERRORN) WRITE (LUPRI,'(A)')
     &     ' Error in nuclear solvent part of gradient '
         CALL HEADER('Difference GSOLNN and GRADN in TLMSTS',-1)
         WRITE (LUPRI,'(11X,A)')
     &      'GSOLNN             GRADN               difference'
         WRITE (LUPRI,'((2X,3F20.10))')
     &         (GSOLNN(I), GRADN(I), GDIFN(I),I=1,3*NUCDEP)
      END IF
C
C     Sum electronic and nuclear contributions
C     ========================================
C
      DO 200 I = 1, 3*NUCDEP
         CALL DAXPY(LMNTOT,DM1,TLMD(1,I,1),1,TLMD(1,I,2),1)
  200 CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('TLMD in Cart. basis',-1)
         WRITE(LUPRI,'(A)') ' Rows and columns: lmn and ICOOR'
         CALL OUTPUT(TLMD(1,1,2),1,LMNTOT,1,3*NUCDEP,
     &               LMNTOT,3*NUCDEP,1,LUPRI)
      END IF
C
C     Transform TLMD to spherical basis
C     =================================
C
      CALL TLMTRA(3*NUCDEP,TLMD(1,1,2),TLMD(1,1,1),WORK,LWORK,IPRINT)
C
C     Write TLMD to disk
C     ==================
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('TLMD in spher. basis as written on LUTLM',-1)
         WRITE(LUPRI,'(A)') ' Rows and columns: l,m and ICOOR'
         CALL OUTPUT(TLMD(1,1,1),1,LMTOT,1,3*NUCDEP,
     &               LMNTOT,3*NUCDEP,1,LUPRI)
      END IF
      CALL GPOPEN(LUTLM,ABATLM,'UNKNOWN','DIRECT',' ',IRAT*LMTOT,OLDDX)
      DO 300 I = 1, 3*NUCDEP
         CALL WRITDX(LUTLM,I,IRAT*LMTOT,TLMD(1,I,1))
  300 CONTINUE
      CALL GPCLOSE(LUTLM,'KEEP')
      RETURN
      END
C  /* Deck tlmtra */
      SUBROUTINE TLMTRA(NMAT,CARMAT,SPHMAT,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION CARMAT(*), SPHMAT(*), WORK(LWORK)
#include "cbisol.h"
C
      MXLM  = 2*LCAVMX + 1
      MXXYZ = (LCAVMX+1)*(LCAVMX+2)/2
      KTRAMA = 1
      KWRK1  = KTRAMA + MXLM*MXXYZ
      IF (KWRK1 .GT. LWORK) CALL STOPIT('TLMTRA',' ',LWORK,KWRK1)
      LWRK1  = LWORK + 1 - KWRK1
      CALL TLMTR1(NMAT,CARMAT,SPHMAT,WORK(KTRAMA),MXXYZ,MXLM,
     &            WORK(KWRK1),LWRK1,IPRINT)
      RETURN
      END
C  /* Deck tlmtr1 */
      SUBROUTINE TLMTR1(NMAT,CARMAT,SPHMAT,TRAMAT,MXXYZ,MXLM,
     &                  WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (MORDER = 1, MINTEG = 0)
      DIMENSION CARMAT(LMNTOT,*), SPHMAT(LMNTOT,*), TRAMAT(MXXYZ,MXLM),
     &          WORK(LWORK)
#include "cbisol.h"
C
      JLM  = 0
      JXYZ = 0
      DO 100 IO = 0, LCAVMX
         NLM  = 2*IO + 1
         NXYZ = (IO+1)*(IO+2)/2
         CALL SPHCOM(IO,TRAMAT,MXLM,MXXYZ,MORDER,MINTEG,WORK,LWORK,
     &               IPRINT)
         DO 200 ICOOR = 1, NMAT
            DO 300 I = 1,NLM
               SPHIJ = D0
               DO 400 J = 1,NXYZ
                  SPHIJ = SPHIJ + CARMAT(JXYZ+J,ICOOR) * TRAMAT(J,I)
  400          CONTINUE
               SPHMAT(JLM+I,ICOOR) = SPHIJ
  300       CONTINUE
  200    CONTINUE
         JLM  = JLM  + NLM
         JXYZ = JXYZ + NXYZ
  100 CONTINUE
      RETURN
      END
C /* Deck PCMNUC */
      SUBROUTINE PCMNUC
C
C     Purpose: calculate nuclear contributions to PCM
C
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0)
#include <pi.h>
#include <mxcent.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
#include <nuclei.h>
#include <pcmnuclei.h>
#include <taysol.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <symmet.h>
#include <ibtfun.h>

      FACT = 4.0D0 * PI * EPS/(EPS - D1)
      ISPP = 1
      DO INUC = 1, NUCIND
         DO ICOOR = 1,3
            SESE=D0
            SNSN=D0
            SESN=D0
            IF(ICESPH.NE.1) CALL PCMNUC_OVER(ISPP,ICOOR,SESE,SNSN,SESN)
            DERU = D0
            DERY = D0
            XN = CORD(1,INUC)
            YN = CORD(2,INUC) 
            ZN = CORD(3,INUC)
            DO JSYM = 0, MAXREP
               SIGNX = PT(IBTAND(ISYMAX(1,1),JSYM))
               SIGNY = PT(IBTAND(ISYMAX(2,1),JSYM))
               SIGNZ = PT(IBTAND(ISYMAX(3,1),JSYM))
               DO ITS = 1, NTSIRR
                  XI = SIGNX * XTSCOR(ITS)
                  YI = SIGNY * YTSCOR(ITS)
                  ZI = SIGNZ * ZTSCOR(ITS)
                  DIST = DSQRT ((XN-XI)**2 + (YN-YI)**2 + 
     $                 (ZN-ZI)**2)
                  IF(ICOOR .EQ. 1) PROD = XI - XN
                  IF(ICOOR .EQ. 2) PROD = YI - YN
                  IF(ICOOR .EQ. 3) PROD = ZI - ZN
                  DVNUC = CHARGE(INUC) * PROD / DIST**3
                  DERY = DERY + DVNUC * QSE(ITS)
                  DERU = DERU + DVNUC * QSN(ITS)
               ENDDO
            ENDDO
            IACOOR  = IPTCNT(3*(INUC-1)+ICOOR ,0,1)
            IF(IACOOR .GT. 0) THEN
               DOVERS = FACT * ((SESE + SNSN) * DP5 + SESN)
               GSOLCV(IACOOR) = FMULT(ISTBNU(INUC))*DOVERS
               GSOLNN(IACOOR) = FMULT(ISTBNU(INUC))*(DERY + DERU)
            ENDIF
         ENDDO
         ISPP = ISPP + MULT(ISTBNU(INUC))
      ENDDO
      RETURN
      END
C/*DECK PCMNUC_OVER*/
      SUBROUTINE PCMNUC_OVER(NSJ,IC,SESE,SNSN,SESN)
C
#include <implicit.h>
      PARAMETER (D0 = 0.0D0)
#include <mxcent.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <symmet.h>
#include <ibtfun.h>
C
      DO ISYM = 0, MAXREP
         MTS = NTSIRR * ISYM
         SIGNX = PT(IBTAND(ISYMAX(1,1),ISYM))
         SIGNY = PT(IBTAND(ISYMAX(2,1),ISYM))
         SIGNZ = PT(IBTAND(ISYMAX(3,1),ISYM))
         DO ITS=1,NTSIRR
            L=ISPHE(ITS + MTS)
            XNI = - (XE(L) - SIGNX * XTSCOR(ITS)) / RE(L)
            YNI = - (YE(L) - SIGNY * YTSCOR(ITS)) / RE(L)
            ZNI = - (ZE(L) - SIGNZ * ZTSCOR(ITS)) / RE(L)
            IF(L.LE.NESFP)THEN
               IF(INA(L).EQ.NSJ)THEN
                  IF(IC.EQ.1)DN=XNI
                  IF(IC.EQ.2)DN=YNI
                  IF(IC.EQ.3)DN=ZNI
                  SESE=SESE+DN*(QSE(ITS)**2)/AS(ITS)
                  SNSN=SNSN+DN*(QSN(ITS)**2)/AS(ITS)
                  SESN=SESN+DN*QSE(ITS)*QSN(ITS)/AS(ITS)
               ENDIF
            ELSE
               DCENTN=XNI*Dercen(L,NSJ,IC,1)+
     *                YNI*Dercen(L,NSJ,IC,2)+
     *                ZNI*Dercen(L,NSJ,IC,3)
               DN=DERRAD(L,NSJ,IC)+DCENTN
               SESE=SESE+DN*(QSE(ITS)**2)/AS(ITS)
               SNSN=SNSN+DN*(QSN(ITS)**2)/AS(ITS)
               SESN=SESN+DN*QSE(ITS)*QSN(ITS)/AS(ITS)
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END
