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
#include <maxorb.h>
#include <cbiexc.h>
#include <cbilnr.h>
#include <abainf.h>
#include <anrinf.h>
#include <dorps.h>
#include <nuclei.h>
#include <absorp.h>
#include <maxaqn.h>
#include <symmet.h>
C
      PARAMETER (NTABLE = 14, CM_TO_AU = 4.556333D-6)
      LOGICAL NEWDEF
      CHARACTER*8 DIPLEN(3)
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      DATA DIPLEN/'XDIPLEN','YDIPLEN','ZDIPLEN'/
C
      DATA TABLE /'.SKIP  ', '.XXXXXX','.MAX IT','.THRESH',
     *            '.MAXRED', '.MAXPHP','.XXXXXX','.XXXXXX',
     *            '.OPTORB', '.FREQUE','.XXXXXX','.PRINT ',
     *            '.LIFETI', '.STOP  '/
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
               CALL QUIT('Illegal keyword under *ABALNR.')
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
                  J = MIN(NFRVAL,MXFR)
                  READ (LUCMD,*) (FRVAL(I),I=1,J)
                  ICHANG = ICHANG + 1
               GO TO 100
   11             CONTINUE
               GO TO 100
   12             CONTINUE
                  READ (LUCMD,*) IPRLNR
                  ICHANG = ICHANG + 1
                  IF (IPRLNR .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
   13             CONTINUE
                  ABSORP = .TRUE.
                  ICHANG = ICHANG + 1
                  READ (LUCMD,*) LIFETIME
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
               CALL QUIT('Illegal prompt under *ABALNR.')
            END IF
      END IF
  300 CONTINUE
C
      CALL IZERO(NOPER,8)
      DO ILABEL=1,3
         ISYM = ISYMAX(ILABEL,1) + 1
         NOPER(ISYM) = NOPER(ISYM) + 1
         LABOP(NOPER(ISYM),ISYM) = DIPLEN(ILABEL)
      END DO
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
               J = MIN(MXFR,NFRVAL)
               WRITE (LUPRI,'(A,I4/A,5F10.6:/,(28X,5F10.6))')
     &            ' Number of frequencies     :',NFRVAL,
     &            ' Frequencies               :',
     &            (FRVAL(I), I=1,J)
             IF (NFRVAL .GT. MXFR) THEN
               WRITE(LUPRI,'(/A/A,I5/A)')
     &         'You have asked for too many frequencies under *ABALNR',
     &         'Current maximum is MXFR =',MXFR,
     &         'Reduce number of frequencies or increase'//
     &         ' MXFR in cbilnr.h and recompile.'
               CALL QUIT(
     &         'You have asked for too many frequencies under *ABALNR')
             END IF
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
            IF (ABSORP) THEN
               WRITE(LUPRI,'(A,F8.2,A)') 
     &              ' Finite lifetime of the excited states :',
     &              LIFETIME, ' cm-1'
               LIFETIME = LIFETIME*CM_TO_AU                  
               NFREQ_ALPHA = NFRVAL
               DO I=1,NFREQ_ALPHA
                  FREQ_ALPHA(I) = FRVAL(I)
               END DO
               THCLR_ABSORP = THCLNR
               THCPP_ABSORP = 1.0D-3
               IPRABS = IPRLNR
               NEXCITED_STATES = 2
               MAX_MACRO  = 5
               MAX_MICRO  = MAXITE
               MAX_ITORB = 5
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
      ALFA     = ABA_ALPHA
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
      LOGICAL CICLC, HFCLC, TRIPLE, EXECLC, FOUND, CONV
      DIMENSION WORK(LWORK),SNDPRP(2)
      DIMENSION POLDD(2,3,3,MXFR), POLDQ(2,3,3,3,MXFR)
      DIMENSION POLDL(2,3,3,MXFR), POLDA(2,3,3,MXFR)
      CHARACTER*8 LABEL1, LABEL2, LABINT(3*MXCOOR), BLANK
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
#include <absorp.h>
C
      IF (SKIP) RETURN
      CALL QENTER('LNRABA')
      CALL TIMER('START ',TIMEIN,TIMOUT)
C
      IF (IPRLNR .GE. 0) 
     &     CALL TITLER('Solving Linear Response Equations','¤',118)
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
      IPRINT=IPRLNR
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
      TRPLET = .FALSE.
      LUSOVE = 0
      LUGDVE = 0
      LUREVE = 0
      NABAOP = 1
      BLANK  = '        '
C
      CALL GPOPEN(LUSOVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUGDVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
      CALL GPOPEN(LUREVE,' ','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
C
      IF (NFRVAL.LE.0) GOTO 999
C
C     Zero the property complex tensors and start the calculations
C
      CALL DZERO(SNDPRP,2)
      CALL DZERO(POLDD,2*9*MXFR)
      CALL DZERO(POLDQ,2*27*MXFR)
      CALL DZERO(POLDL,2*9*MXFR)
      CALL DZERO(POLDA,2*9*MXFR)
C     
      DO ISYM=1,NSYM
         KSYMOP = ISYM
         DO IOPER=1,NOPER(KSYMOP)
            KOPER=IOPER
            LABEL1 = LABOP(IOPER,KSYMOP)
            IF (LABEL1 .EQ.'XDIPLEN') IDIP = 1
            IF (LABEL1 .EQ.'YDIPLEN') IDIP = 2
            IF (LABEL1 .EQ.'ZDIPLEN') IDIP = 3
C     
            CALL LNRVAR(ISYM,IPRLNR,WORK(KWORK2),LWORK2)
C     
            KGD1   = KWORK1
            KWRKG1 = KGD1
            LWRKG1 = LWORK - KWRKG1
            KSLV   = KGD1 + 2*NVARPT
            KLAST  = KSLV + 4*NVARPT
            IF (KLAST.GT.LWORK) CALL STOPIT('LNRABA',' ',KLAST,LWORK)
            KWRK = KLAST
            LWRK = LWORK - KLAST + 1
C     
C     Find right hand side (gradient) of first operator and write to file
C     ===================================================================
C     
            CALL GETGPV(LABEL1,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &           WORK(KPVX),WORK(KXINDX),ANTSYM,WORK(KWRKG1),LWRKG1)
            REWIND LUGDVE
            CALL WRITT(LUGDVE,2*NVARPT,WORK(KWRKG1))
            IF (IPRLNR.GT.10) THEN
               WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL1
               CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,
     &              2,1,LUPRI)
            ENDIF
C     
C     Calculate linear response vector and write to file
C     ==================================================
C     
            CALL ABARSP(CICLC,HFCLC,TRIPLE,OOTV,ISYM,EXECLC,
     &           FRVAL,NFRVAL,NABATY,NABAOP,LABEL1,LUGDVE,LUSOVE,
     &           LUREVE,THCLNR,MAXITE,IPRRSP,MXRM,MXPHP,
     &           WORK(KWRK),LWRK)
C     
C     Loop over the second property operators
C     =======================================
C     
            DO 200 IPRLBL = 1, NLBTOT
C     
C     Find label and symmetry of second operator
C     
               LABEL2 = LABAPP(IPRLBL)
               KSYM   = LABSYM(IPRLBL)
C     
C     If symmetry of first operator equals symmetry of
C     second operator, that is if ISYM = KSYM, then
C     ================================================
C     
               IF (KSYM.EQ.ISYM) THEN
C     
C     Find right hand side (gradient) for second operator
C     ===================================================
C     
                  CALL GETGPV(LABEL2,DUMMY,DUMMY,WORK(KCMO),WORK(KUDV),
     &                 WORK(KPVX),WORK(KXINDX),ANTSYM,
     &                 WORK(KWRKG1),LWRKG1)
C     
                  IF (IPRLNR.GT.10) THEN
                     WRITE (LUPRI,'(2A)') 'GP Vector, label: ',LABEL2
                     CALL OUTPUT(WORK(KGD1),1,NVARPT,1,2,NVARPT,
     &                    2,1,LUPRI)
                  ENDIF
C     
C     Form second order properties SNDPRP
C     ===================================
C     
                  IF (.NOT.ABSORP) REWIND LUSOVE
                  DO 100 IFRVAL = 1,NFRVAL
                     IF (ABSORP) THEN
                        CALL GPOPEN(LURSP,'RSPVEC','OLD',' ',
     &                       ' ',IDUMMY,.FALSE.)
                        CALL REARSP(LURSP,4*NVARPT,WORK(KSLV),LABEL1,
     &                       BLANK,FREQ_ALPHA(IFRVAL),0.0D0,ISYM,0,
     &                       THCLR_ABSORP,FOUND,CONV,ANTSYM)
                        CALL GPCLOSE(LURSP,'KEEP')
                        SNDPRP(1)=DDOT(2*NVARPT,WORK(KSLV),1,
     &                       WORK(KGD1),1)
                        SNDPRP(2)=DDOT(2*NVARPT,WORK(KSLV+2*NVARPT),1,
     &                       WORK(KGD1),1)
                     ELSE
                        CALL READT(LUSOVE,2*NVARPT,WORK(KSLV))
                        SNDPRP(1)=DDOT(2*NVARPT,WORK(KSLV),1,
     &                       WORK(KGD1),1)
                     END IF
C     
                     IF (IPRLNR.GT.2) THEN
                        WRITE (LUPRI,'(/A,F15.8,A/4A,F15.8)')
     &                       ' Frequency = ',FRVAL(IFRVAL),' au',
     &                       ' Second order property for ',
     &                       LABEL1,LABEL2,' = ',SNDPRP(1)
                     ENDIF
                     IF (IPRLNR.GT.10) THEN
                        WRITE(LUPRI,'(A,I4)') 'Solution Vector no. ',
     &                       IFRVAL
                        CALL OUTPUT(WORK(KSLV),1,NVARPT,1,2,NVARPT,
     &                       2,1,LUPRI)
                     ENDIF
C     
C     Write properties into the various property matrices
C     ===================================================
C     
                     IF (LABEL2(2:7).EQ.'DIPLEN') THEN
                        IF (LABEL2(1:1).EQ.'X') THEN
                           POLDD(1,IDIP,1,IFRVAL) = SNDPRP(1)
                           POLDD(2,IDIP,1,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:1).EQ.'Y') THEN
                           POLDD(1,IDIP,2,IFRVAL) = SNDPRP(1)
                           POLDD(2,IDIP,2,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:1).EQ.'Z') THEN
                           POLDD(1,IDIP,3,IFRVAL) = SNDPRP(1)
                           POLDD(2,IDIP,3,IFRVAL) = SNDPRP(2)
                        END IF
                     ELSE IF (LABEL2(3:8).EQ.'THETA ') THEN
                        IF (LABEL2(1:2).EQ.'XX') THEN
                           POLDQ(1,IDIP,1,1,IFRVAL) = SNDPRP(1)
                           POLDQ(2,IDIP,1,1,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:2).EQ.'XY') THEN
                           POLDQ(1,IDIP,1,2,IFRVAL) = SNDPRP(1)
                           POLDQ(2,IDIP,1,2,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:2).EQ.'XY') THEN
                           POLDQ(1,IDIP,2,1,IFRVAL) = SNDPRP(1)
                           POLDQ(2,IDIP,2,1,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:2).EQ.'XZ') THEN
                           POLDQ(1,IDIP,1,3,IFRVAL) = SNDPRP(1)
                           POLDQ(2,IDIP,1,3,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:2).EQ.'XZ') THEN
                           POLDQ(1,IDIP,3,1,IFRVAL) = SNDPRP(1)
                           POLDQ(2,IDIP,3,1,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:2).EQ.'YY') THEN
                           POLDQ(1,IDIP,2,2,IFRVAL) = SNDPRP(1)
                           POLDQ(2,IDIP,2,2,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:2).EQ.'YZ') THEN
                           POLDQ(1,IDIP,2,3,IFRVAL) = SNDPRP(1)
                           POLDQ(2,IDIP,2,3,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:2).EQ.'YZ') THEN
                           POLDQ(1,IDIP,3,2,IFRVAL) = SNDPRP(1)
                           POLDQ(2,IDIP,3,2,IFRVAL) = SNDPRP(2)
                        END IF
                        IF (LABEL2(1:2).EQ.'ZZ') THEN
                           POLDQ(1,IDIP,3,3,IFRVAL) = SNDPRP(1)
                           POLDQ(2,IDIP,3,3,IFRVAL) = SNDPRP(2)
                        END IF
                     ELSE IF (LABEL2(2:7).EQ.'LONMAG') THEN
C     
C     Multiply with minus the Bohr-magneton (-0.5) to create the magnetic 
C     dipole operator from the angular momentum operator
C     
                        IF (LABEL2(1:1).EQ.'X') THEN
                           POLDL(1,IDIP,1,IFRVAL) = -DP5*SNDPRP(1)
                           POLDL(2,IDIP,1,IFRVAL) = -DP5*SNDPRP(2)
                        END IF
                        IF (LABEL2(1:1).EQ.'Y') THEN
                           POLDL(1,IDIP,2,IFRVAL) = -DP5*SNDPRP(1)
                           POLDL(2,IDIP,2,IFRVAL) = -DP5*SNDPRP(2)
                        END IF
                        IF (LABEL2(1:1).EQ.'Z') THEN
                           POLDL(1,IDIP,3,IFRVAL) = -DP5*SNDPRP(1)
                           POLDL(2,IDIP,3,IFRVAL) = -DP5*SNDPRP(2)
                        END IF
                     ELSE IF (LABEL2(2:7).EQ.'ANGMOM') THEN
                        IF (LABEL2(1:1).EQ.'X') THEN
                           POLDA(1,IDIP,1,IFRVAL) = -DP5*SNDPRP(1)
                           POLDA(2,IDIP,1,IFRVAL) = -DP5*SNDPRP(2)
                        END IF
                        IF (LABEL2(1:1).EQ.'Y') THEN
                           POLDA(1,IDIP,2,IFRVAL) = -DP5*SNDPRP(1)
                           POLDA(2,IDIP,2,IFRVAL) = -DP5*SNDPRP(2)
                        END IF
                        IF (LABEL2(1:1).EQ.'Z') THEN
                           POLDA(1,IDIP,3,IFRVAL) = -DP5*SNDPRP(1)
                           POLDA(2,IDIP,3,IFRVAL) = -DP5*SNDPRP(2)
                        END IF
                     ELSE IF (IPRLNR.LE.2) THEN
C                    ... hjaaj mar 2004: this property is not
C                        predefined, should never happen, but
C                        we make sure the value is printed.
C                        (the value was printed above if IPRLNR.gt.2)
                        WRITE (LUPRI,'(/A,F15.8,A/4A,F15.8)')
     &                       ' Frequency = ',FRVAL(IFRVAL),' au',
     &                       ' Second order property for ',
     &                       LABEL1,LABEL2,' = ',SNDPRP(1)
                     ENDIF
 100              CONTINUE
               ENDIF
 200        CONTINUE
C     
C     End loop over operator in this symmetry
C     
         END DO
C     
C     End loop over symmetries
C     
      END DO
C
 999  CONTINUE
C
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
#include <codata.h>
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
#include <codata.h>
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
#include <codata.h>
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
#include <codata.h>
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
#include <codata.h>
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
#include <absorp.h>
#include <inforb.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <codata.h>
      PARAMETER (FACTOR = (288.0D-30)*(PI**2)*XFMOL*(XTANG**4),
     &           D3 = 3.0D0)
      PARAMETER ( CM_TO_AU = 4.556333D-6, AU_TO_EV=27.2114 )
#include <nuclei.h>
#include <symmet.h>
#include <abainf.h>
#include <cbilnr.h>
      DIMENSION POLDD(2,3,3,MXFR), POLDL(2,3,3,MXFR), POLDA(2,3,3,MXFR),
     &          BETANL(2), BETAL(2), ALPHNL(2), ALPHAL(2)
#include <ibtfun.h>
C
      IF (ALFA .AND. IPRINT .GE. 0) THEN
         WRITE(LUPRI,9000)
         DO 100 IFRVAL = 1,NFRVAL
            WRITE(LUPRI,9001) FRVAL(IFRVAL)
            WRITE(LUPRI,9002)
            IF (ABSORP) WRITE(LUPRI,'(/,10X,A)') 'Real part:'
            WRITE(LUPRI,9003) 'Ex','Ey','Ez'
            WRITE(LUPRI,9004) 'Ex',(POLDD(1,1,JDIP,IFRVAL), JDIP=1,3)
            WRITE(LUPRI,9004) 'Ey',(POLDD(1,2,JDIP,IFRVAL), JDIP=1,3)
            WRITE(LUPRI,9004) 'Ez',(POLDD(1,3,JDIP,IFRVAL), JDIP=1,3)
            IF (ABSORP) THEN
               WRITE(LUPRI,'(/,10X,A)') 'Imaginary part:'
               WRITE(LUPRI,9003) 'Ex','Ey','Ez'
               WRITE(LUPRI,9004) 'Ex',(POLDD(2,1,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9004) 'Ey',(POLDD(2,2,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9004) 'Ez',(POLDD(2,3,JDIP,IFRVAL), JDIP=1,3)
            END IF
  100    CONTINUE
      ENDIF
      IF (ROAG) THEN
         CALL TITLER('Optical rotation','+',118)
         IF (IPRINT .GE. 0) THEN
            DO IFRVAL = 1,NFRVAL
               WRITE(LUPRI,9005) FRVAL(IFRVAL)
               WRITE(LUPRI,9007)
               IF (ABSORP) WRITE(LUPRI,'(/,10X,A)') 'Real part:'
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(POLDL(1,1,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(POLDL(1,2,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(POLDL(1,3,JDIP,IFRVAL), JDIP=1,3)
               IF (ABSORP) THEN
                  WRITE(LUPRI,'(/,10X,A)') 'Imaginary part:'
                  WRITE(LUPRI,9009) 'Bx','By','Bz'
                  WRITE(LUPRI,9008) 'Ex',(POLDL(2,1,JDIP,IFRVAL), 
     &                 JDIP=1,3)
                  WRITE(LUPRI,9008) 'Ey',(POLDL(2,2,JDIP,IFRVAL), 
     &                 JDIP=1,3)
                  WRITE(LUPRI,9008) 'Ez',(POLDL(2,3,JDIP,IFRVAL), 
     &                 JDIP=1,3)
               END IF
            END DO
            DO IFRVAL = 1,NFRVAL
               WRITE(LUPRI,9006) FRVAL(IFRVAL)
               WRITE(LUPRI,9007)
               IF (ABSORP) WRITE(LUPRI,'(/,10X,A)') 'Real part:'
               WRITE(LUPRI,9009) 'Bx','By','Bz'
               WRITE(LUPRI,9008) 'Ex',(POLDA(1,1,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ey',(POLDA(1,2,JDIP,IFRVAL), JDIP=1,3)
               WRITE(LUPRI,9008) 'Ez',(POLDA(1,3,JDIP,IFRVAL), JDIP=1,3)
               IF (ABSORP) THEN
                  WRITE(LUPRI,'(/,10X,A)') 'Imaginary part:'
                  WRITE(LUPRI,9009) 'Bx','By','Bz'
                  WRITE(LUPRI,9008) 'Ex',(POLDA(2,1,JDIP,IFRVAL), 
     &                 JDIP=1,3)
                  WRITE(LUPRI,9008) 'Ey',(POLDA(2,2,JDIP,IFRVAL), 
     &                 JDIP=1,3)
                  WRITE(LUPRI,9008) 'Ez',(POLDA(2,3,JDIP,IFRVAL), 
     &                 JDIP=1,3)
               END IF
            END DO
         END IF
         TMASS = 0.0D0
         JATOM = 0
         DO IATOM = 1, NUCIND
            DO ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  JATOM = JATOM + 1
                  NATTYP = IZATOM(IATOM)
                  IF (NATTYP .NE. 0 .AND. .NOT. NOORBT(IATOM)) THEN
                     AMASS = DISOTP(NATTYP,ISOTOP(JATOM),'MASS')
                     TMASS = TMASS + AMASS
                  END IF
               END IF
            END DO
         END DO
         FACTOT = FACTOR*XTKAYS*XTKAYS
         WRITE(LUPRI,'(2(/,10X,A))')
     &        'Optical rotation summary',
     &        '------------------------'
         DO IFRVAL = 1, NFRVAL
            BETAL(1) =-(POLDL(1,1,1,IFRVAL) + POLDL(1,2,2,IFRVAL) +
     &                 POLDL(1,3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            BETANL(1)=-(POLDA(1,1,1,IFRVAL) + POLDA(1,2,2,IFRVAL) +
     &                 POLDA(1,3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            ALPHAL(1)=FACTOT*BETAL(1)*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            ALPHNL(1)=FACTOT*BETANL(1)*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            BETAL(2) =-(POLDL(2,1,1,IFRVAL) + POLDL(2,2,2,IFRVAL) +
     &                 POLDL(2,3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            BETANL(2)=-(POLDA(2,1,1,IFRVAL) + POLDA(2,2,2,IFRVAL) +
     &                 POLDA(2,3,3,IFRVAL))/(D3*FRVAL(IFRVAL))
            ALPHAL(2)=FACTOT*BETAL(2)*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            ALPHNL(2)=FACTOT*BETANL(2)*FRVAL(IFRVAL)*FRVAL(IFRVAL)/TMASS
            WRITE (LUPRI,'(/10X,A,F12.6,A,F12.6,A)') 
     &           'Frequency: ',FRVAL(IFRVAL),' au  = ',
     &           XTNM/FRVAL(IFRVAL), ' nm'
            IF (ABSORP) THEN
               WRITE(LUPRI,'(10X,A,F12.2,A)')
     &              'Lifetime : ',LIFETIME/CM_TO_AU,' cm-1'
               WRITE(LUPRI,'(4(/10X,2(A,F12.6)))') 
     &              'Beta (London)                :',
     &              BETAL(1), '  + i ',BETAL(2),
     &              'Beta (No-London)             :',
     &              BETANL(1),'  + i ',BETANL(2),
     &              'Optical rotation (London)    :', 
     &              ALPHAL(1),'  + i ',ALPHAL(2),
     &              'Optical rotation (No-London) :', 
     &              ALPHNL(1),'  + i ',ALPHNL(2)
            ELSE
               WRITE(LUPRI,'(4(/10X,A,F12.6))') 
     &              'Beta (London)                :',BETAL(1),
     &              'Beta (No-London)             :',BETANL(1),
     &              'Optical rotation (London)    :',ALPHAL(1),
     &              'Optical rotation (No-London) :',ALPHNL(1)
            END IF
         END DO
      ENDIF
C
      IF (ABSORP .AND. IPRINT.GE.0) THEN
         WRITE(LUPRI,'(/10X,A4,A6,A21,/10X,A)') 
     &        'Sym.','State','Excitation energy',
     &        '---------------------------------------------------'
         DO ISYM=1,NSYM
            IF (NOPER(ISYM).GT.0) THEN
               DO ISTATE=1,NEXCITED_STATES
                  WRITE(LUPRI,'(6X,2I6,F12.4,A,F9.3,A,F10.2,A)') 
     &                 ISYM,ISTATE,EXC_ENERGY(ISTATE,ISYM),' a.u.',
     &                 EXC_ENERGY(ISTATE,ISYM)*AU_TO_EV,' eV',
     &                 XTNM/EXC_ENERGY(ISTATE,ISYM),' nm'
               END DO
            END IF
         END DO
      END IF
C
 9000 FORMAT (/)
 9001 FORMAT (/,10X,'Polarizability tensor for frequency ',F12.6,' au')
 9002 FORMAT (10X,'----------------------------------------------------'
     &       )
 9003 FORMAT (16X,3(A12),/)
 9004 FORMAT (16X,A,2X,3F12.7)
 9005 FORMAT (/,10X,'London    G tensor for frequency ',F12.6,' au')
 9006 FORMAT (/,10X,'No-London G tensor for frequency ',F12.6,' au')
9007  FORMAT (10X,'-------------------------------------------------')
9008  FORMAT (16X,A,2X,3F12.7)
9009  FORMAT (27X,3(A,10X)/)
C
      RETURN
      END
