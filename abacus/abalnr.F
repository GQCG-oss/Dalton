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
C  /* Deck lnrinp */
      SUBROUTINE LNRINP(WORD)
C
C     Nov. and Dec. 93
C     Written by K.L.Bak and P.Joergensen using EXCITA as template.
C     Purpose: To enable calculations of frequency dependent second
C     order response properties in ABACUS. In particular polariza-
C     bilities and vibrational raman optical activity (VROA).
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 14)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <cbiexc.h>
#include <cbilnr.h>
#include <abainf.h>
#include <anrinf.h>
#include <dorps.h>
#include <nuclei.h>
      DATA TABLE /'.SKIP  ', '.XXXXXX','.MAX IT','.THRESH',
     *            '.MAXRED', '.MAXPHP','.XXXXXX','.XXXXXX',
     *            '.OPTORB', '.FREQUE','.ALFA  ','.IPRINT',
     *            '.XXXXXX','.STOP  '/
C
      NEWDEF = (WORD .EQ. '*ABALNR')
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
                     GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in LNRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in LNRINP.')
    1          CONTINUE
                  SKIP = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
    2          CONTINUE
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) INPVAL
                  ICHANG = ICHANG + 1
                  IF (INPVAL .EQ. MAXITE) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     MAXITE = INPVAL
                  END IF
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) DTHCLN
                  ICHANG = ICHANG + 1
                  IF (DTHCLN .EQ. THCLNR) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     THCLNR = DTHCLN
                  END IF
               GO TO 100
    5          CONTINUE
                  READ (LUCMD,*) IMXRM
                  ICHANG = ICHANG + 1
                  IF (IMXRM .EQ. MXRM) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     MXRM = IMXRM
                  END IF
               GO TO 100
    6          CONTINUE
                  READ (LUCMD,*) IMXPHP
                  ICHANG = ICHANG + 1
                  IF (IMXPHP .EQ. MXPHP) THEN
                     ICHANG = ICHANG - 1
                  ELSE
                     MXPHP = IMXPHP
                  END IF
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
    9             OOTV   = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   10             CONTINUE
                  READ (LUCMD,*) NFRVAL
                  READ (LUCMD,*) (FRVAL(I),I=1,NFRVAL)
                  ICHANG = ICHANG + 1
               GO TO 100
   11             CONTINUE
                  ALFA = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
   12             CONTINUE
                  READ (LUCMD,*) IPRLNR
                  ICHANG = ICHANG + 1
                  IF (IPRLNR .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   13             CONTINUE
               GO TO 100
   14             CONTINUE
                  CUT    = .TRUE.
                  ICHANG = ICHANG + 1
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in LNRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in LNRINP.')
            END IF
      END IF
  300 CONTINUE
C
      IF ((VROA .OR. RAMAN .OR. OPTROT) .AND. NODIFC) THEN
         WRITE (LUPRI,'(/,3A,/)')
     &   ' Raman properties calculation is NOT ',
     &   ' implemented with the NODIFC keyword active. ',
     &   ' Try calculation again without specifying NODIFC.'
         ROAA = .FALSE.
         ROAG = .FALSE.
      ELSEIF (VROA) THEN
         ROAG = .TRUE.
         ROAA = .TRUE.
      ELSE IF (OPTROT) THEN
         ROAG = .TRUE.
      ELSE IF (RAMAN) THEN
         ALFA = .TRUE.
      ELSE
         ROAA = .FALSE.
         ROAG = .FALSE.
      END IF
C
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for LNRABA:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' LNRABA skipped in this run.'
         ELSE
            IF (NFRVAL .GT. 0) THEN
               WRITE (LUPRI,'(A,I4/A,5F10.6:/,(28X,5F10.6))')
     &            ' Number of frequencies     :',NFRVAL,
     &            ' Frequencies               :',
     &            (FRVAL(I), I=1,NFRVAL)
            END IF
            WRITE (LUPRI,'(A,I4)')
     &         ' Print level in LNRABA     :',IPRLNR
            WRITE (LUPRI,'(A,1P,D10.2)')
     &         ' Threshold in LNRABA       :',THCLNR
            WRITE(LUPRI,'(A,I4)')' Max. iterations in LNRABA :',
     &                           MAXITE
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after LNRINP.'
            END IF
         END IF
      END IF
C
      RETURN
      END
C  /* Deck lnrini */
      SUBROUTINE LNRINI
C
C     Initialize /LNREXC/
C
#include <implicit.h>
#include <mxcent.h>
#include <cbiexc.h>
#include <cbilnr.h>
#include <abainf.h>
C
      IPRINT   = IPRDEF
      IPRLNR   = IPRINT
      SKIP     = .FALSE.
      CUT      = .FALSE.
      OOTV     = .FALSE.
      ALFA     = .FALSE.
      NFRVAL   = 1
      FRVAL(1) = 0.0D0
      THCLNR   = 5.D-05
      MAXITE   = 60
      MXRM     = 400
      MXPHP    = 0
      RETURN
      END
C  /* Deck lnraba */
      SUBROUTINE LNRABA(POLDD,POLDQ,POLDL,POLDA,WORK,LWORK,PASS)
#include <implicit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <iratdef.h>
#include <priunit.h>
#include <cbilnr.h>
      LOGICAL PASS
      LOGICAL CICLC, HFCLC, TRIPLE, EXECLC, FOUND
      DIMENSION WORK(LWORK)
      DIMENSION POLDD(3,3,MXFR), POLDQ(3,3,3,MXFR)
      DIMENSION POLDL(3,3,MXFR), POLDA(3,3,MXFR)
      CHARACTER*8 LABEL1, LABEL2, LABINT(3*MXCOOR)
      CHARACTER*1 XYZLAB
      PARAMETER (DP5=0.5D0)
C
#include <cbiexc.h>
#include <inflin.h>
#include <infvar.h>
#include <infdim.h>
#include <inforb.h>
#include <nuclei.h>
#include <inftap.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include <maxmom.h>
#include <maxaqn.h>
#include <symmet.h>
#include <abainf.h>
#include <gnrinf.h>
#include <infsop.h>
C
      IF (SKIP) RETURN
      CALL QENTER('LNRABA')
      CALL TIMER('START ',TIMEIN,TIMOUT)
      IF (IPRLNR .GT. 0) WRITE (LUPRI,'(/A/)')
     *    ' <<<<<<<<<< Output from LNRABA >>>>>>>>>> '
C
      IPRRSP = IPRLNR
C
C     Get reference state
C     ===================
C
C     1. Work Allocations:
C
      IF (ABASOP) THEN
         LUDV   = NORBT * NORBT
         LPVX   = LPVMAT
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
      IF (.NOT.FOUND) CALL QUIT('LNRABA error: CMO not found on SIRIFC')
      IF (NASHT .GT. 0) THEN
         CALL RDSIFC('DV',FOUND,WORK(KWORK1),WORK(KWORK1),LWORK1)
         IF (.NOT.FOUND)
     &      CALL QUIT('LNRABA error: DV not found on SIRIFC')
         CALL DSPTSI(NASHT,WORK(KWORK1),WORK(KUDV))
      END IF
C
      ISYM = 1
      CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK1),LWORK1)
C
      CALL GETCIX(WORK(KXINDX),IREFSY,IREFSY,WORK(KWORK1),LWORK1,0)
C
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
         IF (IPRLNR.GT.10) THEN
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
C        reads the MP2 or CCSD second order one particle density matrix 
C
         CALL READT (LUSIFC,NORBT*NORBT,WORK(KUDV))
C
C        UDV contains the MP2 one-density. Remove the diagonal
C        contribution from the zeroth order. (Added in MP2FAC)
C
         IF (IPRLNR.GT.10) THEN
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
      END IF
C
C     Construct property-integrals and write to LUPROP
C     ================================================
C
C     2. Work Allocations:
C
      KIDSYM = KWORK1
      KIDADR = KIDSYM + 9*MXCENT
      KWORK2 = KIDADR + 9*MXCENT
      LWORK2 = LWORK  - KWORK2
C
      NLBTOT = 0
C
      IF (ALFA .OR. ROAA .OR. ROAG) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'DIPLEN ',NCOMP,WORK(KWORK2),LWORK2,
     &               LABINT,WORK(KIDSYM),WORK(KIDADR),
     &               IDUMMY,.TRUE.,NPATOM,.TRUE.,IPRINT)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      ENDIF
      IF (ROAA) THEN
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'THETA  ',NCOMP,WORK(KWORK2),LWORK2,
     &            LABINT,WORK(KIDSYM),WORK(KIDADR),
     &            IDUMMY,.TRUE.,NPATOM,.TRUE.,IPRINT)
         NLAB = 6
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
      END IF
C
      IF (ROAG) THEN
         CALL LABCOP(1,NLBTOT,'XLONMAG ',ISYMAX(1,2),LABAPP,LABSYM)
         CALL LABCOP(1,NLBTOT,'YLONMAG ',ISYMAX(2,2),LABAPP,LABSYM)
         CALL LABCOP(1,NLBTOT,'ZLONMAG ',ISYMAX(3,2),LABAPP,LABSYM)
C
         NCOMP  = 0
         NPATOM = 0
         CALL GET1IN(DUMMY,'ANGMOM ',NCOMP,WORK(KWORK2),LWORK2,
     &            LABINT,WORK(KIDSYM),WORK(KIDADR),
     &            IDUMMY,.TRUE.,NPATOM,.TRUE.,IPRINT)
         NLAB = 3
         CALL LABCOP(NLAB,NLBTOT,LABINT,WORK(KIDSYM),LABAPP,LABSYM)
C
      ENDIF
C
C     Set variables and logicals
C
      CICLC  = .FALSE.
      HFCLC  = NASHT .LE. 1
      TRIPLE = .FALSE.
      EXECLC = .FALSE.
      NABATY = 1
      NABAOP = 1
C
C     Assign symmetry of x-, y-, and z-dipole operator to I?SYM
C
      IXSYM  = ISYMAX(1,1) + 1
      IYSYM  = ISYMAX(2,1) + 1
      IZSYM  = ISYMAX(3,1) + 1
C
C     Zero the property tensors
C
      CALL DZERO(POLDD,9*MXFR)
      CALL DZERO(POLDQ,27*MXFR)
      CALL DZERO(POLDL,9*MXFR)
      CALL DZERO(POLDA,9*MXFR)
C
C     Loop over the first operators which are the
C     x, y, and, z electric-dipole operators
C     ===========================================
C
      LUSOVE = 0
      LUGDVE = 0
      LUREVE = 0
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      DO 300 IDIP = 1,3
         IF (IDIP.EQ.1) LABEL1 = 'XDIPLEN '
         IF (IDIP.EQ.2) LABEL1 = 'YDIPLEN '
         IF (IDIP.EQ.3) LABEL1 = 'ZDIPLEN '
         ISYM = ISYMAX(IDIP,1)+1
C
         CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK2),LWORK2)
C
         IF (NFRVAL.GT.0) THEN
C
C           3. Work Allocations:
C
            KGD1   = KWORK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KWRKG1
            KSLV   = KGD1 + 2*NVARPT
            KLAST  = KSLV + 2*NVARPT
            IF (KLAST.GT.LWORK) CALL STOPIT('LNRABA',' ',KLAST,LWORK)
            KWRK = KLAST
            LWRK = LWORK - KLAST + 1
C
C           Find right hand side for first operator and write to file
C           =========================================================
C
            KSYMOP = ISYM
            TRPLET = .FALSE.
C
            CALL GETGPV(LABEL1,DUMMY,WORK(KCMO),WORK(KUDV),WORK(KPVX),
     &                  WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
            REWIND LUGDVE
            CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
            IF (IPRLNR.GT.3) THEN
               WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL1
               CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,
     &                           2,1,LUPRI)
            ENDIF
C
C           Calculate eigenvector and write to file
C           =======================================
C
            CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     &            FRVAL,NFRVAL,NABATY,NABAOP,LABEL1,LUGDVE,LUSOVE,
     &            LUREVE,THCLNR,MAXITE,IPRRSP,MXRM,MXPHP,
     &            WORK(KWRK),LWRK)
C
C           Loop over the second property operators
C           =======================================
C
            DO 200 IPRLBL = 1, NLBTOT
C
C              Find label and symmetry of second operator
C
               LABEL2 = LABAPP(IPRLBL)
               KSYM   = LABSYM(IPRLBL)
C
C              If symmetry of first operator equals symmetry of
C              second operator, that is if ISYM = KSYM, then
C              ================================================
C
               IF (KSYM.EQ.ISYM) THEN
                  KSYMOP = ISYM
                  TRPLET = .FALSE.
C
C                 Find right hand side for second operator
C                 ========================================
C
                  CALL GETGPV(LABEL2,DUMMY,WORK(KCMO),WORK(KUDV),
     &                        WORK(KPVX),WORK(KXINDX),ANTSYM,
     &                        WORK(KWRKG1),LWRKG1)
C
                  IF (IPRLNR.GT.3) THEN
                     WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL2
                     CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,
     &                           2,1,LUPRI)
                  ENDIF
C
C                 Form second order properties SNDPRP
C                 ===================================
C
                  REWIND LUSOVE
                  DO 100 IFRVAL = 1,NFRVAL
                     CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
C
                     SNDPRP = DDOT(2*NVARPT,WORK(KSLV),1,WORK(KGD1),1)
C
                     IF (IPRLNR.GT.3) THEN
                        WRITE (LUPRI,'(A,I4)') 'Solution Vector no. ',
     &                         IFRVAL
                        CALL OUTPUT(WORK(KSLV),1,NVARPT,1,2,NVARPT,
     &                              2,1,LUPRI)
                     ENDIF
                     IF (IPRLNR.GT.2) THEN
                        WRITE (LUPRI,'(/,A,F15.8)')
     &                  ' Frequency = ',FRVAL(IFRVAL)
                        WRITE (LUPRI,'(4A,F15.8)')
     &                  ' Second order property for ',LABEL1,LABEL2,
     &                  ' = ',SNDPRP
                     ENDIF
C
C                    Write properties into the various property matrices
C                    ===================================================
C
                     IF (LABEL2(2:7).EQ.'DIPLEN') THEN
                      IF (LABEL2(1:1).EQ.'X')POLDD(IDIP,1,IFRVAL)=SNDPRP
                      IF (LABEL2(1:1).EQ.'Y')POLDD(IDIP,2,IFRVAL)=SNDPRP
                      IF (LABEL2(1:1).EQ.'Z')POLDD(IDIP,3,IFRVAL)=SNDPRP
                     ELSE IF (LABEL2(3:8).EQ.'THETA ') THEN
                      IF(LABEL2(1:2).EQ.'XX')POLDQ(IDIP,1,1,IFRVAL)=
     &                                       SNDPRP
                      IF(LABEL2(1:2).EQ.'XY')POLDQ(IDIP,1,2,IFRVAL)=
     &                                       SNDPRP
                      IF(LABEL2(1:2).EQ.'XY')POLDQ(IDIP,2,1,IFRVAL)=
     &                                       SNDPRP
                      IF(LABEL2(1:2).EQ.'XZ')POLDQ(IDIP,1,3,IFRVAL)=
     &                                       SNDPRP
                      IF(LABEL2(1:2).EQ.'XZ')POLDQ(IDIP,3,1,IFRVAL)=
     &                                       SNDPRP
                      IF(LABEL2(1:2).EQ.'YY')POLDQ(IDIP,2,2,IFRVAL)=
     &                                       SNDPRP
                      IF(LABEL2(1:2).EQ.'YZ')POLDQ(IDIP,2,3,IFRVAL)=
     &                                       SNDPRP
                      IF(LABEL2(1:2).EQ.'YZ')POLDQ(IDIP,3,2,IFRVAL)=
     &                                       SNDPRP
                      IF(LABEL2(1:2).EQ.'ZZ')POLDQ(IDIP,3,3,IFRVAL)=
     &                                       SNDPRP
                     ELSE IF (LABEL2(2:7).EQ.'LONMAG') THEN
C
C                     Multiply with minus the Bohr-magneton (-0.5) to 
C                     create the magnetic dipole operator from the 
C                     angular momentum operator
C
                      IF (LABEL2(1:1).EQ.'X')POLDL(IDIP,1,IFRVAL)=
     &                                       -DP5*SNDPRP
                      IF (LABEL2(1:1).EQ.'Y')POLDL(IDIP,2,IFRVAL)=
     &                                       -DP5*SNDPRP
                      IF (LABEL2(1:1).EQ.'Z')POLDL(IDIP,3,IFRVAL)=
     &                                       -DP5*SNDPRP
                     ELSE IF (LABEL2(2:7).EQ.'ANGMOM') THEN
                      IF (LABEL2(1:1).EQ.'X')POLDA(IDIP,1,IFRVAL)=
     &                                       -DP5*SNDPRP
                      IF (LABEL2(1:1).EQ.'Y')POLDA(IDIP,2,IFRVAL)=
     &                                       -DP5*SNDPRP
                      IF (LABEL2(1:1).EQ.'Z')POLDA(IDIP,3,IFRVAL)=
     &                                       -DP5*SNDPRP
                     ENDIF
  100             CONTINUE
               ENDIF
  200       CONTINUE
         END IF
  300 CONTINUE
      CALL GPCLOSE(LUSOVE,'DELETE')
      CALL GPCLOSE(LUGDVE,'DELETE')
      CALL GPCLOSE(LUREVE,'DELETE')
C
      CALL TIMER ('LNRABA',TIMEIN,TIMOUT)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/,A)')
     &          ' Program stopped after LNRABA as required.'
         WRITE (LUPRI,'(A)') ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in LNRABA) *****')
      END IF
C
      CALL QEXIT('LNRABA')
      RETURN
      END
C
C  /* Deck lnrvar */
      SUBROUTINE LNRVAR(ISYM,IPRINT,WORK,LWORK)
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
      DIMENSION WORK(LWORK)
#include <inflin.h>
#include <infvar.h>
#include <infrsp.h>
#include <inforb.h>
#include <wrkrsp.h>
#include <abainf.h>
C
      LOGICAL TRIPLE
C
C     Set variables for response modules
C     ==================================
C
      TRIPLE = .FALSE.
      CALL ABAVAR(ISYM,TRIPLE,IPRINT,WORK,LWORK)
      IREFSY = LSYMRF
      NCREF  = NCONRF
      KSYMST = LSYMST
      KSYMOP = LSYMPT
      KZWOPT = NWOPT
      IF ((NASHT .EQ. 0).AND.(.NOT.ABASOP)) THEN
         KZCONF = 0
      ELSE IF (NASHT .EQ. 1) THEN
         KZCONF = 1
      ELSE
         KZCONF = NCONST
      END IF
      KZVAR  = KZWOPT + KZCONF
      KZYWOP = 2*KZWOPT
      KZYCON = 2*KZCONF
      KZYVAR = 2*KZVAR
C
      RETURN
      END
C  /* Deck betagm */
      FUNCTION BETAGM(ALFA,GM)
C
C Calculate Beta(G')**2 =
C BETAGM = 0.5*(3*ALFA(I,J)*GM(I,J) - ALFA(I,I)*GM(J,J))
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0  , DP5 = 0.5D0 , D3 = 3.0D0 )
      DIMENSION ALFA(3,3),GM(3,3)
C
      BETAGM = D0
      DO 100 I = 1,3
         DO 200 J = 1,3
            BETAGM = BETAGM +
     &               DP5*(D3*ALFA(I,J)*GM(I,J)-ALFA(I,I)*GM(J,J))
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C  /* Deck betaal */
      FUNCTION BETAAL(ALFA)
C
C Calculate Beta(alfa)**2 =
C BETAAL =  0.5 * (3*ALFA(I,J)*ALFA(I,J) - ALFA(I,I)*ALFA(J,J))
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0  , DP5 = 0.5D0 , D3 = 3.0D0 )
      DIMENSION ALFA(3,3)
C
      BETAAL = D0
      DO 100 I = 1,3
         DO 200 J = 1,3
            BETAAL = BETAAL + DP5*(D3*ALFA(I,J)*ALFA(I,J)
     *                             -ALFA(I,I)*ALFA(J,J))
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C  /* Deck betaa */
      FUNCTION BETAA(ALFA,A,OMEGA)
C
C Calculate Beta(A)**2 =
C BETAA =  0.5*OMEGA*ALFA(I,J)*EPSILON(I,K,L)*A(K,L,J)
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0 , DP5 = 0.5D0 , D1 = 1.0D0 , DM1 = -1.0D0)
      DIMENSION ALFA(3,3), A(3,3,3)
      DIMENSION EPSILO(3,3,3)
      DATA EPSILO / D0, D0,  D0, D0, D0, DM1, D0,  D1, D0,
     &              D0, D0,  D1, D0, D0, D0,  DM1, D0, D0,
     &              D0, DM1, D0, D1, D0, D0,  D0,  D0, D0/
C
      BETAA = D0
      DO 100 I = 1,3
         DO 200 J = 1,3
            DO 300 K = 1,3
               DO 400 L = 1,3
                  BETAA = BETAA + DP5*OMEGA*ALFA(I,J)*EPSILO(I,K,L)*
     *                                A(K,L,J)
 400           CONTINUE
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C  /* Deck alfmn */
      FUNCTION ALFMN(ALFA)
C
C Calculate AlfaMean =  ALFMN = ( (1/3) * ALFA(I,I) )
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D3 = 3.0D0 )
      DIMENSION ALFA(3,3)
C
      ALFMN = D0
      DO 100 I = 1,3
         ALFMN = ALFMN + ALFA(I,I)
 100  CONTINUE
C
      ALFMN = ALFMN  / D3
C
      RETURN
      END
C  /* Deck gmmn */
      FUNCTION GMMN(GM)
C
C Calculate G' =  GMMN = ( (1/3) * GM(I,I) )
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D3 = 3.0D0 )
      DIMENSION GM(3,3)
C
      GMMN = D0
      DO 100 I = 1,3
         GMMN = GMMN + GM(I,I)
 100  CONTINUE
C
      GMMN = GMMN  / D3
C
      RETURN
      END
C  /* Deck raminl */
      FUNCTION RAMINL(ALFMN,BETAAL)
C
C Calculate RamanIntensity =  RAMIN = (45 * ALFMN**2) + (4 * BETA**2)
C for linear pol. incident radiation perpendicular to scattering plane
C 
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D4 = 4.0D0, D45 = 45.0D0 )
C
      RAMINL = (D45 * ALFMN**2) + (D4 * BETAAL)
C
      RETURN
      END
C  /* Deck depoll */
      FUNCTION DEPOLL(ALFMN,BETAAL)
C
C Calculate the Depolarization Ratio for rightangle scattering,
C linear polarized incident light, perpendicular (=parall/perpend)
C plane = DEPOLR = (3 * BETA**2) / (45 * ALFMN**2 + 4 * BETA**2)
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D3 = 3.0D0, D4 = 4.0D0, D45 = 45.0D0 )
C
      A      = D3 * BETAAL
      B      = (D45 * ALFMN**2) + (D4 * BETAAL)
C
      IF ( B .GT. 0.D-6 ) THEN
         DEPOLL = A / B
      ELSE
         DEPOLL = D0
      ENDIF
C
      RETURN
      END
C  /* Deck raminn */
      FUNCTION RAMINN(ALFMN,BETAAL)
C
C Calculate RamanIntensity =  RAMIN = (45 * ALFMN**2) + (7 * BETA**2)
C for natural & circular pol. incident radiation 
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D7 = 7.0D0, D45 = 45.0D0 )
C
      RAMINN = (D45 * ALFMN**2) + (D7 * BETAAL)
C
      RETURN
      END
C  /* Deck depoln */
      FUNCTION DEPOLN(ALFMN,BETAAL)
C
C Calculate the Depolarization Ratio for rightangle scattering, 
C naturel or circ. pol. incident light (=parall/perpend).
C plane = DEPOLR = (6 * BETA**2) / (45 * ALFMN**2 + 7 * BETA**2)
C
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D6 = 6.0D0, D7 = 7.0D0, D45 = 45.0D0 )
C
      A      = D6 * BETAAL
      B      = (D45 * ALFMN**2) + (D7 * BETAAL)
C
      IF ( B .GT. 0.D-6 ) THEN
         DEPOLN = A / B
      ELSE
         DEPOLN = D0
      ENDIF
C
      RETURN
      END
C  /* Deck deltaz */
      FUNCTION DELTAZ(BETAGM,BETAA)
C
C Calculate the difference differential scattering cross section  
C in depolarized rightangle scattering =
C DELTAZ = 4/c ( 6*BETAGM - 2*BETAA)
C
#include <implicit.h>
#include <alphac.h>
      PARAMETER ( D2 = 2.0D0, D4 = 4.0D0, D6 = 6.0D0 )
C
      DELTAZ = D4 / CVEL * (D6 * BETAGM - D2 * BETAA)
C
      RETURN
      END
C  /* Deck deltax */
      FUNCTION DELTAX(BETAGM,BETAA,ALFAMN,GMMN)
C
C Calculate the difference differential scattering cross section
C in polarized rightangle scattering =
C DELTAX = 4/c ( 45*ALFAMN*GMMN+7*BETAGM+BETAA)
C
#include <implicit.h>
#include <alphac.h>
      PARAMETER ( D4 = 4.0D0, D7 = 7.0D0, D45 = 45.0D0)
C
      DELTAX = D4 / CVEL * (D45*ALFAMN*GMMN+D7*BETAGM+BETAA)
C
      RETURN
      END
C  /* Deck delta0 */
      FUNCTION DELTA0(BETAGM,BETAA,ALFAMN,GMMN)
C
C Calculate the difference differential scattering cross section
C in forward scattering with no analyzer =
C DELTA0 = 4/c (180*ALFAMN*GMMN+4*BETAGM-4*BETAA)
C
#include <implicit.h>
#include <alphac.h>
      PARAMETER ( D4 = 4.0D0, D180 = 18.0D1 )
C
      DELTA0 = D4 / CVEL * (D180*ALFAMN*GMMN+D4*BETAGM-D4*BETAA)
C
      RETURN
      END
C  /* Deck deltab */
      FUNCTION DELTAB(BETAGM,BETAA)
C
C Calculate the difference differential scattering cross section
C in backward scattering with no analyzer =
C DELTAB = 4/c ( 24*BETAGM+8*BETAA)
C
#include <implicit.h>
#include <alphac.h>
      PARAMETER ( D4 = 4.0D0, D8 = 8.0D0, D24 = 24.0D0 )
C
      DELTAB = D4 / CVEL * (D24*BETAGM + D8*BETAA)
C
      RETURN
      END
C  /* Deck cid */
      FUNCTION CID(DELTA,RMINN)
C
C Calculate the Circular Intensity Difference (CID) for all
C arragements. CID = - DELTA  / (2*RMINN)
C RMIN(0) and (180) with no analyzer = 2*RMIN.
C Rmin (depolarized) = DEPLN * RMINN
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, DP5 = 0.5D0 )
C
      IF ( RMINN .GT. 0.D-6 ) THEN
         CID = DELTA / RMINN * DP5
      ELSE
         CID = D0
      ENDIF
C
      RETURN
      END
C  /* Deck boltzk */
      FUNCTION BOLTZK(FREQ)
C
C Calculate the Boltzmannfactor for Spectra simulation 
C BOLTZK = 1/(1-EXP(-h*freq/k/T)), FREQ is in Hartree = 2Pi*freq
C
#include <implicit.h>
#include <codata.h>
      PARAMETER ( TSTAND = 298.15D0 )
C
      BOLTZK = 1-EXP(-(FREQ*AUTK)/TSTAND)
      RETURN
      END
C  /* Deck crossk */
      FUNCTION CROSSK(FRVAL,FREQ)
C
C Calculate the Constant for the differential scattering cross section
C K = 16/90*Pi**4*((freq0-freq)/C)**3*freq0/C
C Combined with the freq. factor out of PLACZEK approximation
C B = SQRT(1/(4Pi*freq))
C CROSSK = K*B**2
C ****************** FREQ is in Hartree = 2Pi*freq
C
#include <implicit.h>
#include <alphac.h>
      PARAMETER ( D180 = 18.0D1)
C
      CROSSK = (FRVAL-FREQ)**3*FRVAL/FREQ/D180/CVEL**3
      RETURN
      END
C  /* Deck lnrout */
      SUBROUTINE LNROUT(POLDD,POLDL,POLDA,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pi.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <codata.h>
      PARAMETER (FACTOR = (288.0D-30)*(PI**2)*XFMOL*(XTANG**4),
     &           D3 = 3.0D0)
#include <nuclei.h>
#include <symmet.h>
#include <abainf.h>
#include <cbilnr.h>
      DIMENSION POLDD(3,3,MXFR), POLDL(3,3,MXFR), POLDA(3,3,MXFR)
#include <ibtfun.h>
C
      IF (ALFA .AND. IPRINT .GE. 0) THEN
         WRITE(LUPRI,9000)
         DO 100 IFRVAL = 1,NFRVAL
            WRITE(LUPRI,9001) FRVAL(IFRVAL)
            WRITE(LUPRI,9002)
            WRITE(LUPRI,9003) 'Ex','Ey','Ez'
            WRITE(LUPRI,9004) 'Ex',(POLDD(1,JDIP,IFRVAL), JDIP=1,3)
            WRITE(LUPRI,9004) 'Ey',(POLDD(2,JDIP,IFRVAL), JDIP=1,3)
            WRITE(LUPRI,9004) 'Ez',(POLDD(3,JDIP,IFRVAL), JDIP=1,3)
  100    CONTINUE
      ENDIF
      IF (ROAG) THEN
         CALL TITLER('Optical rotation','+',118)
         WRITE(LUPRI,9000)
         IF (IPRINT .GE. 0) THEN
            DO 110 IFRVAL = 1,NFRVAL
               WRITE(LUPRI,9005) FRVAL(IFRVAL)
               WRITE(LUPRI,9007)
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(POLDL(1,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(POLDL(2,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(POLDL(3,JDIP,IFRVAL), JDIP=1,3)
 110        CONTINUE
            DO 120 IFRVAL = 1,NFRVAL
               WRITE(LUPRI,9006) FRVAL(IFRVAL)
               WRITE(LUPRI,9007)
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(POLDA(1,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(POLDA(2,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(POLDA(3,JDIP,IFRVAL), JDIP=1,3)
 120        CONTINUE
         END IF
         TMASS = 0.0D0
         JATOM = 0
         DO IATOM = 1, NUCIND
            DO ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  JATOM = JATOM + 1
                  NATTYP = NINT(CHARGE(IATOM))
                  IF (NATTYP .NE. 0 .AND. .NOT. NOORBT(IATOM)) THEN
                     AMASS = DISOTP(NATTYP,ISOTOP(JATOM),'MASS')
                     TMASS = TMASS + AMASS
                  END IF
               END IF
            END DO
         END DO
         FACTOT = FACTOR*XTKAYS*XTKAYS
         DO IFRVAL = 1, NFRVAL
            BETAL  = -(POLDL(1,1,IFRVAL) + POLDL(2,2,IFRVAL) +
     &                 POLDL(3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            BETANL = -(POLDA(1,1,IFRVAL) + POLDA(2,2,IFRVAL) +
     &                 POLDA(3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            ALPHAL = FACTOT*BETAL*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            ALPHNL = FACTOT*BETANL*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            WRITE (LUPRI,'(//1X,A,F12.6,A,5X,F12.6,A)') 
     &           'Frequency: ',FRVAL(IFRVAL),' au',
     &           XTNM/FRVAL(IFRVAL), ' nm'
            WRITE (LUPRI,'(/25X,A,F12.6)') 'Beta(London)     :', BETAL
            WRITE (LUPRI,'(25X,A,F12.6)') 'Beta(No-London)  :', BETANL
            WRITE (LUPRI,'(/25X,A,F12.6)') 'Optical rotation '//
     &           '(London)    :', ALPHAL
            WRITE (LUPRI,'(25X,A,F12.6)') 'Optical rotation '//
     &           '(No-London) :', ALPHNL
         END DO
      ENDIF
C
C
 9000 FORMAT (/)
 9001 FORMAT (/,10X,'Polarizability tensor for frequency ',F12.6,' au')
 9002 FORMAT (10X,'----------------------------------------------------'
     &       )
 9003 FORMAT (/,24X,3(A,5X)/)
 9004 FORMAT (16X,A,2X,3F12.7)
 9005 FORMAT (/,10X,'London    G tensor for frequency ',F12.6,' au')
 9006 FORMAT (/,10X,'No-London G tensor for frequency ',F12.6,' au')
9007  FORMAT (10X,'-------------------------------------------------')
9008  FORMAT (16X,A,2X,3F12.7)
9009  FORMAT (/,27X,3(A,10X)/)
C
      RETURN
      END
