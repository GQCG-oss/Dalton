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
C
C$Id: abaexc.F,v 1.6 2002-02-12 08:17:21 vebjornb Exp $
C
#include <single.h>
C  /* Deck excinp */
      SUBROUTINE EXCINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 16)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <inforb.h>
#include <cbiexc.h>
#include <cbilrs.h>
C
#include <abainf.h>
#include <anrinf.h>
#include <dorps.h>
#include <nuclei.h>
      DATA TABLE /'.SKIP  ', '.INTPRI','.MAX IT','.THRESH',
     *            '.MAXRED', '.MAXPHP','.TRIPLE','.XXXXXX',
     *            '.OPTORB', '.NEXCIT','.DIPSTR','.PRINT ',
     *            '.ROTVEL', '.FNAC  ','.STOP  ','.SUMRUL'/
C
      NEWDEF = (WORD .EQ. '*EXCITA')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in EXCINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in EXCINP.')
    1          CONTINUE
                  SKIP = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    2          CONTINUE
                  READ (LUCMD,*) IPRINT
                  ICHANG = ICHANG + 1
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) MAXITE
                  ICHANG = ICHANG + 1
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) THREXC
                  ICHANG = ICHANG + 1
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) MXRM
                  ICHANG = ICHANG + 1
               GO TO 100
    6          CONTINUE
                  READ (LUCMD,*) MXPHP
                  ICHANG = ICHANG + 1
               GO TO 100
    7          CONTINUE
                  EXCTRP = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    8          CONTINUE
               GO TO 100
    9          CONTINUE
                  OOTV   = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   10          CONTINUE
C                 NSYM is set to 8 since NSYM has not yet been set !!
                  READ (LUCMD, '(8I5)') (NEXCIT(I),I=1,8)
                  ICHANG = ICHANG + 1
               GO TO 100
   11          CONTINUE
                  DIPSTR = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   12          CONTINUE
                  READ (LUCMD,*) IPREXC
                  ICHANG = ICHANG + 1
               GO TO 100
   13          CONTINUE
                  ROTVEL = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   14          CONTINUE
                  FNAC   = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   15          CONTINUE
                  CUT    = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   16          CONTINUE
                  SUMRUL = .TRUE.
                  DIPSTR = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in EXCINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in EXCINP.')
            END IF
      END IF
  300 CONTINUE
C
      IF (ECD) THEN
         ROTSTR = .TRUE.
      ELSE
         ROTSTR = .FALSE.
      END IF
C
C      NSYM is set to 8 since NSYM has not yet been set !!
C      NTOT = ISUM(NSYM,NEXCIT,1)
      NTOT = ISUM(8,NEXCIT,1)
C
CSPAS : let's try whether we cannot calculate all excitation energies
C
C     maximum number of excitation energies
C
      MXNEXI = 1
C     NSYM is set to 8 since NSYM has not yet been set !!
      DO ISYM = 1,8
         MXNEXI = MAX(MXNEXI,NEXCIT(ISYM))
      END DO
CKeinSPASmehr
C
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for EXCITA:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' EXCITA skipped in this run.'
         ELSE
            IF (NTOT .GT. 0) THEN
C               NSYM is set to 8 since NSYM has not yet been set !!
               WRITE (LUPRI,'(A,8I5)')
     &            ' Number of excitation energies:',(NEXCIT(I),I=1,8)
            END IF
            WRITE (LUPRI,'(A,I5)')
     &         ' Print level in EXCITA        :',IPREXC
            WRITE (LUPRI,'(A,I5)')
     &         ' Integral print level in EXCITA        :',IPRINT
            WRITE (LUPRI,'(A,1P,E9.2)')
     &           ' Threshold in EXCITA          :',THREXC
            WRITE(LUPRI,'(A,I5)')' Maximum iterations in EXCITA :',
     &                           MAXITE
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after EXCITA.'
            END IF
            IF (ROTVEL) THEN
               WRITE (LUPRI,'(A)') ' Conventional rotational strength'
            END IF
            IF (DIPSTR) THEN
               WRITE (LUPRI,'(A)') ' Dipole strength'
            END IF
            IF (EXCTRP) THEN
               WRITE (LUPRI,'(A)') ' Triplet excitation energies'
            END IF
            IF (FNAC) THEN
               WRITE (LUPRI,'(A)') ' First order nonadiabatic '//
     &                             'coupling elements'
            END IF
         END IF
      END IF
      IF (NTOT.EQ.0) THEN
C         NSYM is set to 8 since NSYM has not yet been set !!
C         DO 400 I = 1, NSYM
         DO 400 I = 1, 8
            NEXCIT(I) = 1
  400    CONTINUE
      END IF
      RETURN
      END
C  /* Deck excini */
      SUBROUTINE EXCINI
C
C     Initialize /LRSRES/
C
#include <implicit.h>
#include <mxcent.h>
#include <cbiexc.h>
#include <abainf.h>
C
      IPRINT = IPRDEF
      IPREXC = IPRINT
      SKIP   = .FALSE.
      CUT    = .FALSE.
      OOTV   = .FALSE.
      DIPSTR = .FALSE.
      ROTVEL = .FALSE.
      FNAC   = .FALSE.
      EXCTRP = .FALSE.
      CALL IZERO(NEXCIT,8)
      THREXC = 1.D-04
      MAXITE = 60
      MXRM   = 400
      MXPHP  = 0
      NABAPP = 0
C
      RETURN
      END
C  /* Deck excita */
      SUBROUTINE EXCITA(TRLEN,TRVEL,TRLON,TRMAG,BSRLON,EXENG,FONAC,
     &                  FONA2,WORK,LWORK,PASS)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <inforb.h>
#include <cbiexc.h>
#include <infrsp.h>
#include <cbilrs.h>
#include <nuclei.h>
      LOGICAL PASS
      DIMENSION WORK(LWORK)
      DIMENSION TRLEN(3,NSYM,MXNEXI), TRVEL(3,NSYM,MXNEXI)
      DIMENSION TRLON(3,NSYM,MXNEXI), TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI),EXENG(NSYM,MXNEXI)
      DIMENSION FONAC(3*NATOMS,NSYM,MXNEXI),FONA2(3*NATOMS,NSYM,MXNEXI)
C
      IF (SKIP) RETURN
      CALL QENTER('EXCITA')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPREXC .GT. 0) WRITE (LUPRI,'(A,/)')
     *    '  <<<<<<<<<< Output from EXCITA >>>>>>>>>> '
C
CSPAS : I am trying to separate the print level in the integral routines
C       from the overall print level in EXCITA
C
C      IPRINT = IPREXC
CKeinSPASmehr
      IPRRSP = IPREXC
C
C     Work allocations:
C
C
      CALL EXCIT1(TRLEN,TRVEL,TRLON,TRMAG,BSRLON,EXENG,FONAC,FONA2,
     &            WORK,LWORK)
C
      CALL TIMER ('EXCITA',TIMEIN,TIMOUT)
      PASS   = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/A)')
     &          ' Program stopped after EXCITA as requested.'
         CALL QUIT(' ***** End of ABACUS (in EXCITA) *****')
      END IF
      CALL QEXIT('EXCITA')
      RETURN
      END
C  /* Deck excit1 */
      SUBROUTINE EXCIT1(TRLEN,TRVEL,TRLON,TRMAG,BSRLON,EXENG,FONAC,
     &                  FONA2,WORK,LWORK)
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <cbilrs.h>
      DIMENSION WORK(LWORK)
      DIMENSION TRLEN(3,NSYM,MXNEXI), TRVEL(3,NSYM,MXNEXI)
      DIMENSION TRLON(3,NSYM,MXNEXI), TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI),EXENG(NSYM,MXNEXI)
      DIMENSION FONAC(3*NATOMS,NSYM,MXNEXI),FONA2(3*NATOMS,NSYM,MXNEXI)
      CHARACTER*8 LABEL, LABINT(3*MXCOOR)
      LOGICAL   FOUND
C
#include <cbiexc.h>
#include <inforb.h>
#include <nuclei.h>
#include <inflin.h>
#include <infdim.h>
#include <inftap.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <maxmom.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <abainf.h>
#include <gnrinf.h>
#include <infsop.h>
C
C     Get reference state
C     ===================
C
C     1. Work allocations:
C
      IF (ABASOP) THEN
         LUDV   = NORBT * NORBT
         LPVX   = LPVMAT
         IF (EXCTRP) LPVX = 2*LPVMAT
      ELSE
         LUDV   = N2ASHX
         LPVX   = 0
      ENDIF
      KFREE  = 1
      LFREE  = LWORK
C      
      CALL MEMGET('REAL',KCMO  ,NCMOT ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDV  ,LUDV  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPVX  ,LPVX  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXINDX,LCINDX,WORK,KFREE,LFREE)
C
      KWORK1 = KFREE
      LWORK1 = LFREE
C
      CALL RDSIFC('CMO',FOUND,WORK(KCMO),WORK(KFREE),LFREE)
      IF (.NOT.FOUND) CALL QUIT('EXCIT1 error: CMO not found on SIRIFC')
      IF (NASHT .GT. 0) THEN
         CALL RDSIFC('DV',FOUND,WORK(KWORK1),WORK(KWORK1),LWORK1)
         IF (.NOT.FOUND)
     &      CALL QUIT('EXCIT1 error: DV not found on SIRIFC')
         CALL DSPTSI(NASHT,WORK(KWORK1),WORK(KUDV))
      END IF
C
      LFREE = LWORK1
C
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LFREE,0)
C
C     SOPPA :
C
      IF (ABASOP) THEN
C
C        Initialize XINDX
C
         CALL DZERO(WORK(KXINDX),LCINDX)
C
C        Find address array's for SOPPA calculation
C
         CALL SET2SOPPA(WORK(KXINDX+KABSAD-1),WORK(KXINDX+KABTAD-1),
     *                  WORK(KXINDX+KIJSAD-1),WORK(KXINDX+KIJTAD-1),
     *                  WORK(KXINDX+KIJ1AD-1),WORK(KXINDX+KIJ2AD-1),
     *                  WORK(KXINDX+KIJ3AD-1),WORK(KXINDX+KIADR1-1))
C
C
         REWIND (LUSIFC)
         IF (CCPPA) THEN
            CALL MOLLAB('CCSDINFO',LUSIFC,LUPRI)
         ELSE
            CALL MOLLAB('MP2INFO ',LUSIFC,LUPRI)
         ENDIF
C
C        reads the MP2 or CCSD correlation coefficients into PV
C
         CALL READT (LUSIFC,LPVMAT,WORK(KPVX))
C
         IF (IPREXC.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' EXCIT1 : CCSD correlation ',
     &                           'coefficients'
            ELSE
               WRITE(LUPRI,'(/A,A)')' EXCIT1 :',
     &                              ' MP2 correlation coefficients'
            ENDIF
            CALL OUTPUT(WORK(KPVX),1,LPVMAT,1,1,LPVMAT,1,1,LUPRI)
         END IF
C
         IF (EXCTRP) THEN
            CALL TRPKAP(WORK(KPVX),WORK(KPVX+LPVMAT),
     &                  WORK(KXINDX+KIADR1-1),WORK(KWORK1))
C
            IF (IPREXC.GT.10) THEN
               IF (CCPPA) THEN
                  WRITE(LUPRI,'(/A)')' RSPMC : CCSD triplet ',
     &                              'correlation coefficients'
               ELSE
                  WRITE(LUPRI,'(/A)')' RSPMC : MP2 triplet correlation',
     &                               ' coefficients'
               ENDIF
               CALL OUTPUT(WORK(KPVX),1,2*LPVMAT,1,1,2*LPVMAT,1,1,LUPRI)
            END IF
         END IF
C
C        reads the MP2 or CCSD second order one particle density matrix 
C
         CALL READT (LUSIFC,NORBT*NORBT,WORK(KUDV))
C
C        UDV contains the MP2 one-density. Remove the diagonal
C        contribution from the zeroth order. (Added in MP2FAC)
C
         IF (IPREXC.GT.10) THEN
            IF (CCPPA) THEN
               WRITE(LUPRI,'(/A)')' RSPMC : CCSD density'
            ELSE
               WRITE(LUPRI,'(/A)')' RSPMC : MP2 density'
            END IF
            CALL OUTPUT(WORK(KUDV),1,NORBT*NORBT,1,1,NORBT*NORBT,1,1,
     &                  LUPRI)
         END IF
C
         CALL SOPUDV(WORK(KUDV))
C
      END IF
C
C      CALL MCORL(WORK(KWORK1),LWORK1,PASS)
C
C     Construct property-integrals and write to LUPROP
C     ================================================
C
C     2. Work allocations:
C
      KIDSYM = KWORK1
      KIDADR = KIDSYM + 9*MXCENT
      KWORK2 = KIDADR + 9*MXCENT
      LWORK2 = LWORK  - KWORK2
C
      NLBTOT = 0
C
      IF (DIPSTR .OR. ROTSTR) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'DIPLEN ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,IPRINT)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
C
         IF (ROTSTR) THEN
            CALL LABCOP(1,NLBTOT,'XLONMAG ',ISYMAX(1,2),LABAPP,LABSYM)
            CALL LABCOP(1,NLBTOT,'YLONMAG ',ISYMAX(2,2),LABAPP,LABSYM)
            CALL LABCOP(1,NLBTOT,'ZLONMAG ',ISYMAX(3,2),LABAPP,LABSYM)
C
            IF (NODIFC) THEN
               NCOMP  = 0
               NPATOM = 0
               CALL GET1IN(DUMMY,'HBDO   ',NCOMP,WORK(KWORK2),LWORK2,
     &                     LABINT,WORK(KIDSYM),WORK(KIDADR),
     &                     IDUMMY,.TRUE.,NPATOM,.TRUE.,IPRINT)
               NLAB = 3
               CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,
     &                     LABSYM)
            ENDIF
         END IF
      ENDIF
C
      IF (ROTVEL) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'DIPVEL ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,IPRINT)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'ANGMOM ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,IPRINT)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      ENDIF
C
      IF (FNAC  .AND. NODIFC) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'HDO    ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),IDUMMY,
     &               .TRUE.,NPATOM,.TRUE.,IPRINT)
         NLAB = 3 * NUCDEP
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      ENDIF
C
      CALL DZERO(TRLEN,3*NSYM*MXNEXI)
      CALL DZERO(TRVEL,3*NSYM*MXNEXI)
      CALL DZERO(TRLON,3*NSYM*MXNEXI)
      CALL DZERO(TRMAG,3*NSYM*MXNEXI)
      CALL DZERO(BSRLON,3*NSYM*MXNEXI)
      CALL DZERO(FONAC,3*NUCDEP*NSYM*MXNEXI)
      CALL DZERO(FONA2,3*NUCDEP*NSYM*MXNEXI)
C
C     Loop over all symmetries
C     ========================
C
      LUSOVE = -9002
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      DO 100 ISYM = 1, NSYM
C
         CALL ABAVAR(ISYM,EXCTRP,IPREXC,WORK(KWORK2),LWORK2)
C
         NEXVAL = MIN(NEXCIT(ISYM),NVARPT)
C
         IF (NEXVAL.GT.0) THEN
C
C           3. Work allocations:
C
            KEXVAL = KWORK1
            KWRK1  = KEXVAL + NEXVAL
            LWRK1  = LWORK - KWRK1
            KGD1   = KWRK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KGD1
            KSLV   = KGD1 + 2*NVARPT
            KLAST  = KSLV + 2*NVARPT
            IF (KLAST.GT.LWORK) CALL STOPIT('EXCITA',' ',KLAST,LWORK)
            KWRK = KLAST
            LWRK = LWORK - KLAST + 1
C
C           Find excitation energies and excitation vectors
C           Excitation vectors are written to file
C           ===============================================
C
            CALL EXCIT2(WORK(KEXVAL),NEXVAL,ISYM,LUSOVE,
     &                  WORK(KWRK1),LWRK1)
C
            DO 70 IEXVAL = 1,NEXVAL
               EXENG(ISYM,IEXVAL) = WORK(KEXVAL-1+IEXVAL)
   70       CONTINUE
C
C           Loop over property labels
C           =========================
C
            DO 90 IPRLBL = 1, NLBTOT
              LABEL = LABAPP(IPRLBL)
              KSYM  = LABSYM(IPRLBL)
C
C             ISYM = Reference state symmetry times excited state symmetry.
C
              IF (KSYM.EQ.ISYM) THEN
                KSYMOP = ISYM
                CALL GETGPV(LABEL,DUMMY,WORK(KCMO),WORK(KUDV),
     &               WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
                IF (IPREXC.GT.3) THEN
                  WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL
                   CALL OUTPUT(WORK(KGD1),1,2*NVARPT,1,1,2*NVARPT,
     &                         1,1,LUPRI)
                ENDIF
C
C               Form transition moments
C               =======================
C
                Rewind LUSOVE
                DO 80 IEXVAL = 1,NEXVAL
                  CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
                  EXMOM = DDOT(2*NVARPT,WORK(KSLV),1,WORK(KGD1),1)
C
                  IF (IPREXC.GT.3) THEN
                     WRITE (LUPRI,'(A,I4)') 'Solution Vector no. ',
     &                      IEXVAL
                     CALL OUTPUT(WORK(KSLV),1,2*NVARPT,1,1,2*NVARPT,
     &                           1,1,LUPRI)
                     WRITE (LUPRI,'(A,F15.8)')
     &               'Excitation energy = ',WORK(KEXVAL-1+IEXVAL)
                     WRITE (LUPRI,'(3A,F15.8)')
     &               'Transition moment for ',LABEL,' = ',EXMOM
                  ENDIF
C
C                 Write transition momets into relevant arrays.
C                 ==============================================
C                 We have calculated <f|Operator|0> but want to
C                 keep them as <O|Operator|f>. Therefore, for anti-
C                 hermitian operators we multiply with a minus-sign.
C
                  IF (LABEL(2:7).EQ.'DIPLEN') THEN
                    IF (LABEL(1:1).EQ.'X') TRLEN(1,ISYM,IEXVAL) =  EXMOM
                    IF (LABEL(1:1).EQ.'Y') TRLEN(2,ISYM,IEXVAL) =  EXMOM
                    IF (LABEL(1:1).EQ.'Z') TRLEN(3,ISYM,IEXVAL) =  EXMOM
                  ELSE IF (LABEL(2:7).EQ.'DIPVEL') THEN
                    IF (LABEL(1:1).EQ.'X') TRVEL(1,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Y') TRVEL(2,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Z') TRVEL(3,ISYM,IEXVAL) = -EXMOM
                  ELSE IF (LABEL(2:7).EQ.'ANGMOM') THEN
                    IF (LABEL(1:1).EQ.'X') TRMAG(1,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Y') TRMAG(2,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Z') TRMAG(3,ISYM,IEXVAL) = -EXMOM
                  ELSE IF (LABEL(2:7).EQ.'LONMAG') THEN
                    IF (LABEL(1:1).EQ.'X') TRLON(1,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Y') TRLON(2,ISYM,IEXVAL) = -EXMOM
                    IF (LABEL(1:1).EQ.'Z') TRLON(3,ISYM,IEXVAL) = -EXMOM
                  ELSE IF (LABEL(2:6).EQ.'HBDO ') THEN
                    IF (LABEL(7:7).EQ.'X') BSRLON(1,ISYM,IEXVAL)=
     &                                     -EXMOM * EXENG(ISYM,IEXVAL)
                    IF (LABEL(7:7).EQ.'Y') BSRLON(2,ISYM,IEXVAL)=
     &                                     -EXMOM * EXENG(ISYM,IEXVAL)
                    IF (LABEL(7:7).EQ.'Z') BSRLON(3,ISYM,IEXVAL)=
     &                                     -EXMOM * EXENG(ISYM,IEXVAL)
                  ELSE IF (LABEL(1:4).EQ.'HDO ') THEN
                    DO 57 INR = 1, 3*NUCDEP
                       READ(LABEL(5:6),'(I1)') NR
                       IF (INR.EQ.NR) FONAC(NR,ISYM,IEXVAL) = -EXMOM
   57               CONTINUE
                  ENDIF
   80           CONTINUE
              ENDIF
   90       CONTINUE
         END IF
  100 CONTINUE
      CALL GPCLOSE(LUSOVE,'DELETE')
C
      RETURN
      END
C  /* Deck excout */
      SUBROUTINE EXCOUT(TRLEN,TRVEL,TRMAG,TRLON,BSRLON,EXENG,FONAC,
     &                  FONA2,WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <gnrinf.h>
#include <cbiexc.h>
#include <cbilrs.h>
#include <inforb.h>
#include <nuclei.h>
#include <abainf.h>
      CHARACTER COORDI*21, NMDTXT*9
      DIMENSION TRLEN(3,NSYM,MXNEXI), TRVEL(3,NSYM,MXNEXI)
      DIMENSION TRLON(3,NSYM,MXNEXI), TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI),EXENG(NSYM,MXNEXI)
      DIMENSION FONAC(3*NATOMS,NSYM,MXNEXI),FONA2(3*NATOMS,NSYM,MXNEXI)
      DIMENSION WORK(LWORK)
      DIMENSION DSSUML(-6:2,1:4),DLSUML(-6:2,1:4),DISUML(-6:2,1:4)
      DIMENSION DSSUMM(-6:2,1:4),DLSUMM(-6:2,1:4),DISUMM(-6:2,1:4)
      DIMENSION DSSUMV(-6:2,1:4),DLSUMV(-6:2,1:4),DISUMV(-6:2,1:4)
      PARAMETER (D05 = 0.5D0, DTWO =  2.0D0, DTHREE = 3.0D0,
     &           ZERO = 0.0D0,
     &          ESUDIP = 64604.885D0, ESUECD = 471.44360D0)
#include <codata.h>
C
C     Write transition moments, dip. strength, and rot. strength
C     to outputfile.
C     Transitions moments are keept as <O|Operator|f>. Forming
C     products like <O|Operator1|f><f|Operator2|O> we must therefore
C     multiply with a minus sign if Operator2 is antihermitian.
C     Notice: <O|r|f>        =  <f|r|O>
C             <O|grad|f>     = -<f|grad|O>
C             <O|r X grad|f> = -<f|r X grad|O>
C
C
C     Write electric transition dipole moments to output:
C     ====================================================
C
      IF (.NOT. (DIPSTR .OR. ROTSTR .OR. ROTVEL)) THEN
         WRITE (LUPRI,'()')
         IF (EXCTRP) THEN
            CALL HEADER('Triplet electronic excitation energies',15)
         ELSE
            CALL HEADER ('Singlet electronic excitation energies',15)
         END IF
         WRITE (LUPRI,'(14X,A,/,14X,A,/,14X,A)')
     &      ' Sym.   Mode   Frequency    Frequency',
     &      'ex. st.  No.      (au)          (eV)',
     &      '---------------------------------------'
         DO 15 ISYM = 1, NSYM
            DO 14 IEXVAL = 1, NEXCIT(ISYM)
               WRITE (LUPRI,'(16X,I2,6X,I3,2F12.6)')
     &            ISYM,IEXVAL,EXENG(ISYM,IEXVAL),
     &            EXENG(ISYM,IEXVAL)*XTEV
 14         CONTINUE
 15      CONTINUE
      END IF
      IF (DIPSTR .OR. ROTSTR .OR. ROTVEL) THEN
         WRITE (LUPRI,'()')
         CALL HEADER('Electric transition dipole moments (in a.u.)',15)
         WRITE (LUPRI,'(1X,A,A,2(/,1X,A,A))')
     &   ' Sym.   Mode    Frequency       Velocity/Frequency  ',
     &   '            Length       ',
     &   'ex. st.  No.      (au)          x       y       z   ',
     &   '      x       y       z  ',
     &   '--------------------------------------------------',
     &   '----------------------------'
         DO 200 ISYM = 1, NSYM
            DO 100 IEXVAL = 1,NEXCIT(ISYM)
               WRITE (LUPRI,'(2X,I2,6X,I3,1X,F12.6,3X,3F8.3,2X,3F8.3)')
     &            ISYM, IEXVAL, EXENG(ISYM,IEXVAL),
     &            TRVEL(1,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TRVEL(2,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TRVEL(3,ISYM,IEXVAL) / EXENG(ISYM,IEXVAL),
     &            TRLEN(1,ISYM,IEXVAL),  TRLEN(2,ISYM,IEXVAL),
     &            TRLEN(3,ISYM,IEXVAL)
  100       CONTINUE
  200    CONTINUE
      END IF
C
C     Write magnetic transition dipole meoments to output:
C     ====================================================
C
      IF (ROTSTR .OR. ROTVEL) THEN
         WRITE (LUPRI,'()')
         CALL HEADER('Magnetic transition dipole moments (in a.u.)',15)
         WRITE (LUPRI,'(1X,A,A,2(/,1X,A,A))')
     &   ' Sym.   Mode    Frequency        Conventional     ',
     &   '            London       ',
     &   'ex. st.  No.      (eV)          x       y       z   ',
     &   '      x       y       z  ',
     &   '--------------------------------------------------',
     &   '----------------------------'
         DO 400 ISYM = 1,NSYM
            DO 300 IEXVAL = 1,NEXCIT(ISYM)
C              The factor 0.5 is the Bohr-magneton
               WRITE (LUPRI,'(2X,I2,6X,I3,1X,F12.6,3X,3F8.3,2X,3F8.3)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            D05*TRMAG(1,ISYM,IEXVAL), D05*TRMAG(2,ISYM,IEXVAL),
     &            D05*TRMAG(3,ISYM,IEXVAL), D05*TRLON(1,ISYM,IEXVAL),
     &            D05*TRLON(2,ISYM,IEXVAL), D05*TRLON(3,ISYM,IEXVAL)
               IF (NODIFC) THEN
                  WRITE (LUPRI,'(51X,3(A,F7.3))')
     &            '+', D05*BSRLON(1,ISYM,IEXVAL),
     &            '+', D05*BSRLON(2,ISYM,IEXVAL),
     &            '+', D05*BSRLON(3,ISYM,IEXVAL)
                  WRITE (LUPRI,'(51X,3(A,F7.3))')
     &            '=', D05*(TRLON(1,ISYM,IEXVAL)+BSRLON(1,ISYM,IEXVAL)),
     &            '=', D05*(TRLON(2,ISYM,IEXVAL)+BSRLON(2,ISYM,IEXVAL)),
     &            '=', D05*(TRLON(3,ISYM,IEXVAL)+BSRLON(3,ISYM,IEXVAL))
               END IF
  300       CONTINUE
  400    CONTINUE
      END IF
C
C     Write oscilator and rotational strengths to output:
C     ===================================================
C
      IF (ROTSTR .OR. ROTVEL) THEN
         WRITE (LUPRI,'()')
         CALL HEADER('Oscillator and Rotational Strengths',
     &   18)
         WRITE (LUPRI,'(2X,A,/2X,A)')
     &      'Units: 10**(-40) (esu**2)*(cm**2) (rotational strength)',
     &      '       dimensionless              (oscillator strength)'
         WRITE (LUPRI,'(3(/,1X,A,A))')
     &   ' Sym.   Mode     Frequency       Oscillator-strength ',
     &   '     Rotational-strength   ',
     &   'ex. st.  No.        (eV)          velocity  length   ',
     &   '  velocity  length   London',
     &   '------------------------------------------------',
     &   '-------------------------------'
         DO 600 ISYM = 1,NSYM
            DO 500 IEXVAL = 1,NEXCIT(ISYM)
               DIPSTV = DDOT(3,TRVEL(1,ISYM,IEXVAL),1,
     &                  TRVEL(1,ISYM,IEXVAL),1) / EXENG(ISYM,IEXVAL)
               DIPSTL = DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                  TRLEN(1,ISYM,IEXVAL),1) * EXENG(ISYM,IEXVAL)
C              The factor 0.5 is the Bohr-magneton
               ROTV   = - D05 * DDOT(3,TRVEL(1,ISYM,IEXVAL),1,
     &                  TRMAG(1,ISYM,IEXVAL),1) / EXENG(ISYM,IEXVAL)
               ROTL   = - D05 * DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                  TRMAG(1,ISYM,IEXVAL),1)
               ROTLON = - D05 * DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                  TRLON(1,ISYM,IEXVAL),1)
               IF (NODIFC) THEN
                  ROTLO2 = - D05 * DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                     BSRLON(1,ISYM,IEXVAL),1)
                  ROTLON = ROTLON + ROTLO2
               ENDIF
               WRITE (LUPRI,'(2X,I2,6X,I3,2X,F12.6,5X,2F9.4,3X,3F9.3)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            (DTWO/DTHREE)*DIPSTV, (DTWO/DTHREE)*DIPSTL,
     &            ESUECD*ROTV,   ESUECD*ROTL,   ESUECD*ROTLON
  500       CONTINUE
  600    CONTINUE
C
      ELSE IF (DIPSTR .AND. (.NOT. ROTSTR) .AND. (.NOT. ROTVEL)) THEN
         WRITE (LUPRI,'()')
         CALL HEADER('Oscillator strengths',30)
         WRITE (LUPRI,'(2X,A)')
     &      'Oscillator strengths are dimensionless.'
         WRITE (LUPRI,'(3(/,1X,A))')
     &   ' Sym.   Mode        Frequency     Oscillator-strength ',
     &   'ex. st.  No.           (eV)        velocity   length   ',
     &   '-------------------------------------------------------'
         DO 800 ISYM = 1,NSYM
            DO 700 IEXVAL = 1,NEXCIT(ISYM)
               DIPSTV = DDOT(3,TRVEL(1,ISYM,IEXVAL),1,
     &                  TRVEL(1,ISYM,IEXVAL),1) / EXENG(ISYM,IEXVAL)
               DIPSTL = DDOT(3,TRLEN(1,ISYM,IEXVAL),1,
     &                  TRLEN(1,ISYM,IEXVAL),1) * EXENG(ISYM,IEXVAL)
               WRITE (LUPRI,'(2X,I2,6X,I3,5X,F12.6,5X,2F9.4)')
     &            ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            (DTWO/DTHREE)*DIPSTV, (DTWO/DTHREE)*DIPSTL
  700       CONTINUE
  800    CONTINUE
      END IF
C
      IF (FNAC) THEN
         KCSTRA = 1
         KSCTRA = KCSTRA + 9*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP
         IF (KLAST .GT. LWORK) CALL STOPIT('EXCOUT',' ',KLAST,LWORK)
         CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),1,3*NUCDEP,0)
         WRITE (LUPRI,'()')
         CALL HEADER
     &   ('Non adiabatic coupling matrix elements <O|d/dR|f> (in a.u.)',
     &    9)
         WRITE (LUPRI,'(1X,A,A,2(/,1X,A,A))')
     &   '    Coordinate R       Sym.  Mode  Ex. energy   WFR-term ',
     &   ' BSR-term   Total ',
     &   'No.   Cartesian no.     st.   No.     (Ev)               ',
     &   '                 ',
     &   '---------------------------------------------------------',
     &   '------------------'
         DO 1100 ISYM = 1, NSYM
            DO 1000 IEXVAL = 1,NEXCIT(ISYM)
               DO 900 ISCOOR = 1, 3*NATOMS
                  CALL GENCOR(WORK(KCSTRA),3*NUCDEP,ISCOOR,COORDI)
                  WRITE (LUPRI,'(1X,I2,2X,A,2X,I2,4X,I2,1X,F10.6,1X,
     &                   3F10.4)')
     &            ISCOOR, COORDI, ISYM, IEXVAL, XTEV*EXENG(ISYM,IEXVAL),
     &            FONAC(ISCOOR,ISYM,IEXVAL), FONA2(ISCOOR,ISYM,IEXVAL),
     &            FONAC(ISCOOR,ISYM,IEXVAL) + FONA2(ISCOOR,ISYM,IEXVAL)
  900          CONTINUE
               WRITE (LUPRI,'(A)') ' '
 1000       CONTINUE
 1100    CONTINUE
      END IF
C
C     Calculate oscillator strength sum rules
C     =======================================
C
      IF (SUMRUL) THEN
         CALL DZERO(DSSUML,4*9)
         CALL DZERO(DLSUML,4*9)
         CALL DZERO(DSSUMM,4*9)
         CALL DZERO(DLSUMM,4*9)
         CALL DZERO(DSSUMV,4*9)
         CALL DZERO(DLSUMV,4*9)
         DO ISYM = 1,NSYM
            DO IEXVAL = 1,NEXCIT(ISYM)
               DO ICOM = 1,3
                  DO K = -6,2
                     DSSUML(K,ICOM) = DSSUML(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K+1)
     &                              * TRLEN(ICOM,ISYM,IEXVAL) 
     &                              * TRLEN(ICOM,ISYM,IEXVAL)
                     DLSUML(K,ICOM) = DLSUML(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K+1)
     &                              * DLOG(EXENG(ISYM,IEXVAL))
     &                              * TRLEN(ICOM,ISYM,IEXVAL) 
     &                              * TRLEN(ICOM,ISYM,IEXVAL)
                     DSSUMM(K,ICOM) = DSSUMM(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**K
     &                              * TRLEN(ICOM,ISYM,IEXVAL) 
     &                              * TRVEL(ICOM,ISYM,IEXVAL)
                     DLSUMM(K,ICOM) = DLSUMM(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**K
     &                              * DLOG(EXENG(ISYM,IEXVAL))
     &                              * TRLEN(ICOM,ISYM,IEXVAL) 
     &                              * TRVEL(ICOM,ISYM,IEXVAL)
                     DSSUMV(K,ICOM) = DSSUMV(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K-1)
     &                              * TRVEL(ICOM,ISYM,IEXVAL) 
     &                              * TRVEL(ICOM,ISYM,IEXVAL)
                     DLSUMV(K,ICOM) = DLSUMV(K,ICOM) 
     &                              + DTWO * EXENG(ISYM,IEXVAL)**(K-1)
     &                              * DLOG(EXENG(ISYM,IEXVAL))
     &                              * TRVEL(ICOM,ISYM,IEXVAL) 
     &                              * TRVEL(ICOM,ISYM,IEXVAL)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         DO K = -6,2
            DSSUML(K,4) = (DSSUML(K,1)+DSSUML(K,2)+DSSUML(K,3)) / DTHREE
            DLSUML(K,4) = (DLSUML(K,1)+DLSUML(K,2)+DLSUML(K,3)) / DTHREE
            DSSUMM(K,4) = (DSSUMM(K,1)+DSSUMM(K,2)+DSSUMM(K,3)) / DTHREE
            DLSUMM(K,4) = (DLSUMM(K,1)+DLSUMM(K,2)+DLSUMM(K,3)) / DTHREE
            DSSUMV(K,4) = (DSSUMV(K,1)+DSSUMV(K,2)+DSSUMV(K,3)) / DTHREE
            DLSUMV(K,4) = (DLSUMV(K,1)+DLSUMV(K,2)+DLSUMV(K,3)) / DTHREE
            DO ICOM = 1,4
               IF (DSSUML(K,ICOM).EQ.ZERO) THEN
                  DISUML(K,ICOM) = ZERO
               ELSE
                  DISUML(K,ICOM) = DEXP(DLSUML(K,ICOM)/DSSUML(K,ICOM))
     &                             *XTEV
               ENDIF
               IF (DSSUMM(K,ICOM).EQ.ZERO) THEN
                  DISUMM(K,ICOM) = ZERO
               ELSE
                  DISUMM(K,ICOM) = DEXP(DLSUMM(K,ICOM)/DSSUMM(K,ICOM))
     &                             *XTEV
               ENDIF
               IF (DSSUMV(K,ICOM).EQ.ZERO) THEN
                  DISUMV(K,ICOM) = ZERO
               ELSE
                  DISUMV(K,ICOM) = DEXP(DLSUMV(K,ICOM)/DSSUMV(K,ICOM))
     &                             *XTEV
               ENDIF
            ENDDO
         ENDDO
C
         WRITE (LUPRI,'()')
         CALL HEADER('Oscillator strength sum rules',30)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'S(K) Sum Rules : Dipole Length Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DSSUML(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'S(K) Sum Rules : Dipole Mixed Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DSSUMM(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'S(K) Sum Rules : Dipole Velocity Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DSSUMV(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'L(K) Sum Rules : Dipole Length Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DLSUML(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'L(K) Sum Rules : Dipole Mixed Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DLSUMM(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'L(K) Sum Rules : Dipole Velocity Approximation in a.u.',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DLSUMV(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'I(K) Sum Rules : Dipole Length Approximation in eV',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DISUML(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'I(K) Sum Rules : Dipole Mixed Approximation in eV',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DISUMM(K,J),J=1,4),K=-6,2)
         WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,3X,A,3X,A,6X,A,/)')
     &   'I(K) Sum Rules : Dipole Velocity Approximation in eV',
     &   'K','xx - component','yy - component','zz - component','total'
         WRITE (LUPRI,'(9(5X,I3,4(4X,G13.6)/))') 
     &         (K,(DISUMV(K,J),J=1,4),K=-6,2)
      ENDIF
C     
C        *** If numerical derivatives are calculated, then ***
C        *** write properties to file.                     ***
C
      IF (NMWALK) THEN
         DO ISYM = 1, NSYM
            NMDTXT = '         '
            DO IEXVAL=1, NEXCIT(ISYM)
               IF (IEXVAL.LT.10) THEN
                  WRITE (NMDTXT(5:6),'(2I1)') ISYM, IEXVAL
               ELSE
                  WRITE (NMDTXT(5:7),'(I1,I2)') ISYM, IEXVAL
               END IF
               NMDTXT(1:4) = 'EXEN'
               CALL WRAVFL(EXENG(ISYM,IEXVAL),1,1,1,NMDTXT,IPRINT)
               IF (DIPSTR) THEN
                  NMDTXT(1:4) = 'TRDP' // CHAR(ISYM) // CHAR(IEXVAL)
                  CALL WRAVFL(TRLEN(1,ISYM,IEXVAL),3,1,1,NMDTXT,IPRINT)
               END IF
            END DO
         END DO
      END IF

      RETURN
      END
C  /* Deck excit2 */
      SUBROUTINE EXCIT2(EXVAL,NEXENG,ISYM,LUSOVE,WORK,LWORK)
C
C     Triplet excitations now determined in input
C
#include <implicit.h>
#include <maxorb.h>
#include <mxcent.h>
#include <dummy.h>
#include <priunit.h>
      CHARACTER*8 LABEL
      LOGICAL INTTRS, CICLC, HFCLC, TRIPLE, EXECLC
      DIMENSION WORK(LWORK), EXVAL(NEXENG)
C
#include <cbiexc.h>
#include <infvar.h>
#include <inforb.h>
C
      CICLC  = .FALSE.
      HFCLC  = NASHT .LE. 1
      TRIPLE = EXCTRP
      EXECLC = .TRUE.
      NABATY = 1
      NABAOP = 1
      LABEL  = 'EXCITLAB'
C
C     The following line is needed in order to pass
C     consistency checks in RSPMC
C
      CALL ABAVAR(1,.FALSE.,0,WORK,LWORK)
      IF (HFCLC) NCONF = 1
C
      LUREVE = -9801
      LUGDVE = -9802
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     *            EXVAL,NEXENG,NABATY,NABAOP,LABEL,LUGDVE,LUSOVE,LUREVE,
     *            THREXC,MAXITE,IPREXC,MXRM,MXPHP,WORK,LWORK)
C
C     The following line is needed to reset variables
C
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
      CALL ABAVAR(ISYM,TRIPLE,0,WORK,LWORK)
C
      RETURN
      END
C  /* Deck labcop */
      SUBROUTINE LABCOP(NLAB,NLBTOT,LAB1,ISYM1,LAB2,ISYM2)
C
C     Purpose: To copy a set of NLAB labels and corresponding
C     symmetries from LAB1 and ISYM1 to LAB2 and ISYM2.
C
#include <implicit.h>
#include <priunit.h>
#include <cbiexc.h>
      CHARACTER LAB1(*)*8, LAB2(*)*8
      DIMENSION ISYM1(*), ISYM2(*)
C
      DO 100 ILAB = 1, NLAB
         LAB2(NLBTOT + ILAB)  = LAB1(ILAB)
         ISYM2(NLBTOT + ILAB) = ISYM1(ILAB) + 1
  100 CONTINUE
C
C     Update number total number of labels
      NLBTOT = NLBTOT + NLAB
C
      IF (NLBTOT.GT.MAXPP) THEN
         WRITE(LUPRI,'(A,/A,I5,A,I5)')
     *   ' NUMBER OF SPECIFIED PROPERTIES EXCEED THE MAXIMUM ALLOWED',
     *   ' MAXPRP = ',MAXPP,' NLBTOT = ',NLBTOT
         CALL QUIT(' LABCOP: TOO MANY PROPERTIES SPECIFIED')
      END IF
C
      RETURN
      END
C  /* Deck gencor */
      SUBROUTINE GENCOR(CSTRA,NCOOR,ISCOOR,COORDI)
C
C     Purpose: To generate a readable description of the symmetry
C     coordinate ISCOOR in COORDI.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
#include <symmet.h>
#include <nuclei.h>
      CHARACTER COORDI*21, CNBASE*2
      DIMENSION NBASE(8),NSIGN(8), CSTRA(NCOOR,NCOOR)
#include <chrxyz.h>
#include <chrnos.h>
#include <chrsgn.h>
C
      DO 400 ISYM = 0, MAXREP
      IF (NCRREP(ISYM,1) .GT. 0) THEN
         DO 300 IATOM = 1, NUCIND
            DO 200 ICOOR = 1, 3
               JSCOOR = IPTCNT(3*(IATOM - 1) + ICOOR,ISYM,1)
               IF (JSCOOR.EQ.ISCOOR) THEN
                  NB = 0
                  DO 100 I = 1, 3*NUCDEP
                     NSGN = NINT(CSTRA(ISCOOR,I))
                     IF (NSGN .NE. 0) THEN
                        NB = NB + 1
                        NBASE(NB) = I
                        NSIGN(NB) = NSGN
                     END IF
  100             CONTINUE
                  IF (NB .EQ. 1) THEN
                     WRITE(CNBASE,'(I2)') NBASE(1)
                     COORDI = NAMEX(3*IATOM)(1:4)//' '//CHRXYZ(-ICOOR)
     &                        //CNBASE
C                     WRITE (LUPRI,'(1X,I3,3X,A,2X,A,3X,I2)') ISCOOR,
C     &                  NAMEX(3*IATOM)(1:4), CHRXYZ(-ICOOR), NBASE(1)
                  ELSE IF (NB .LT. 4) THEN
                     WRITE(CNBASE,'(I2)') NBASE(1)
                     COORDI = NAMEX(3*IATOM)(1:4)//' '//CHRXYZ(-ICOOR)
     &                        //'['//CNBASE
                     J = 10
                     DO 150 IB = 2, NB
                        WRITE(CNBASE,'(I2)') NBASE(IB)
                        COORDI(J:J+2) = CHRSGN(NSIGN(IB))//CNBASE
                        J = J + 3
  150                CONTINUE
                     COORDI(J:J+2) = ']/'//CHRNOS(NB)
C                     WRITE (LUPRI,'(1X,I3,3X,A,2X,A,3X,A,I2,'
C     &                //CHRNOS(NB-1)//'(2X,A,1X,I2),A,I1)')
C     &                ISCOOR,NAMEX(3*IATOM)(1:4),CHRXYZ(-ICOOR),' [',
C     &                NBASE(1),(CHRSGN(NSIGN(I)),NBASE(I),I=2,NB),
C     &                ' ]/',NB
                  ELSE
                     COORDI = '    See above    '
                  END IF
               END IF
  200       CONTINUE
  300    CONTINUE
      END IF
  400 CONTINUE
C
      RETURN
      END
