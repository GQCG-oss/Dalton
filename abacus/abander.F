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
C  /* Deck nmdinp */
      SUBROUTINE NMDINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 9)
      CHARACTER PROMPT*1, WORD*7, GRPTMP*15, TABLE(NTABLE)*7, WORD1*7
#include <numder.h>
#include <fcsym.h>
#include <cbinum.h>
      LOGICAL NEWDEF
C
      DATA TABLE /'.DORDR ', '.SYMMET', '.SDRTST', '.RESTRT', '.DRYRUN',
     *            '.PARTST', '.NORMAL', '.XXXXXX', '.XXXXXX'/
C
      NEWDEF = (WORD .EQ. '*NMDDRV')
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
                     GO TO (1,2,3,4,5,6,7,8,9), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in NMDINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in NMDINP.')
 1             CONTINUE
                  READ (LUCMD,'(2I4)') NMORDR, NAORDR
                  GOTO 100
 2             CONTINUE
                  READ (LUCMD,'(A)') GRPTMP
                  IJ = 0
                  DO II = 1, 15
                     IF (GRPTMP(II:II).NE.' ') THEN
                        IJ = IJ + 1
                        FCLASS(IJ:IJ) = GRPTMP(II:II)
                     END IF
                  END DO
                  GOTO 100
 3             CONTINUE
                  SDRTST = .TRUE.
                  GOTO 100
 4             CONTINUE
                  RESTRT = .TRUE.
                  GOTO 100
 5             CONTINUE
                  DRYRUN = .TRUE.
                  READ (LUCMD,'( I4)')  NMREDU
                  READ (LUCMD,'(6I4)') (KDRYRN(II),II=1,NMREDU)
                  GOTO 100
 6             CONTINUE
                  DRYRUN = .TRUE.
                  READ (LUCMD,'( I4)')  NMREDU
                  READ (LUCMD,'(6I4)') (KDRYRN(II),II=1,NMREDU)
 7             CONTINUE
                  NRMCRD = .TRUE.
 8             CONTINUE
 9             CONTINUE
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in NMDINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in NMDINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for NMDINP:',0)
         WRITE (LUPRI,'(A)') '                                    '
         WRITE (LUPRI,'(5X,A,I4,A,I4,a)') 'Numerical derivatives' //
     &             'calculated to order', NMORDR, ' from analytical', 
     &                                    NAORDR, '. derivatives'
         WRITE (LUPRI,'(5x,A)') 'Group used for force constants ' // 
     &                           GRPTMP
         IF (SDRTST) THEN
            WRITE (LUPRI,'(A)') '                                    '
            WRITE (LUPRI,'(5X,A)') 'Comparison of numerical Hessian ' //
     &                           'with analytical Hessian is performed'
         END IF
         IF (DRYRUN) THEN
            WRITE (LUPRI,'(A)') '                                    '
            WRITE (LUPRI,'(5X,A)') 'Numerical derivatives will be' //
     &                 'conducted as a dry run.'
            WRITE (LUPRI,'(5X,A)') 'No actual derivatives will be' //
     &                             'calculated.'
            WRITE (LUPRI,'(A,I4)') 'Number of rendundant coordinates:', 
     &                             NMREDU
         END IF
      END IF
      RETURN
      END
C
C
C
      SUBROUTINE NMDINI
C
C     Initialize /NUMDER/, /NMDPRP/ and /FCSYM/
C
#include <implicit.h>
#include <numder.h>
#include <fcsym.h>
#include <cbinum.h>
C
C     /CBINUM/
      NRMCRD = .FALSE.
C
C     /NUMDER/
      NMORDR = 0
      NAORDR = 0
      NWPROP = .FALSE.
      FSTPRP = .FALSE.
      NOMOVE = .TRUE.
      NUMELC = .FALSE.
      CMPARE = .FALSE.
      SDRTST = .FALSE.
      DRYRUN = .FALSE.
C      
C     /FCSYM/
      FCLASS = 'C1             '
      MROTAX = .FALSE.
      VPLANE = .FALSE.
      HPLANE = .FALSE.
      ROTAX2 = .FALSE.
      DPLANE = .FALSE.
      ICNTR  = .FALSE.
      ROTARE = .FALSE.
C
      RETURN
      END
C
C
C     /* Deck numdrv */
      SUBROUTINE NUMDRV(WORK,LWORK,IWOFF,WRKDLM)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <infpar.h>
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
      LOGICAL FINISH
#include <gnrinf.h>
#endif
#include <cbiwlk.h>
#include <cbinum.h>
#include <trkoor.h>
#include <nuclei.h>
#include <inftap.h>
#include <cbiexc.h>
#include <symmet.h>
#include <numder.h>
#include <molinp.h>
#include <fcsym.h>
      DIMENSION WORK(LWORK)
      LOGICAL MOLECU
      CHARACTER WORD*7
C
 8    CONTINUE
C
      MOLECU = .TRUE.
      PGMTST = .TRUE.
      MDGTST = .TRUE.
      HTEST  = .TRUE.
      NCOOR  = 3*NUCDEP
C
      NMEXCI = 0
      DO 100 ISYM = 0, MAXREP
         NMEXCI = NMEXCI + NEXCIT(ISYM+1)
 100  CONTINUE
C
      NDERIV = 0
      NDIME  = 1
      NINNER = NCOOR**NAORDR
C     
      MAXADR = 1
      DO 200 J = NMORDR, 1, -2
      DO 200 I = 1, J
         MAXADR = MAXADR + 1
 200  CONTINUE
C
      DO 300 IORDR = 3, NMORDR+NAORDR
         IKDRV = 1
         DO 400 IIORDR =1, IORDR
            IKDRV = IKDRV*(NCOOR+IIORDR-1)/IIORDR
 400     CONTINUE
         NDERIV = NDERIV + IKDRV
 300  CONTINUE
C
      DO 500 IORDR = 1, NMORDR
         IKDIME = 1
         DO 600 IIORDR = 1, IORDR
            IKDIME = IKDIME*(NCOOR+1-IIORDR)/IIORDR
 600     CONTINUE
         NDIME = NDIME + IKDIME*2**IORDR
 500  CONTINUE
      IF (NMORDR .GE. 3) THEN
         NDIME = NDIME + 2*NCOOR
      END IF
      IF (NMORDR .GE. 4) THEN
         NDIME = NDIME + 4*NCOOR*(NCOOR-1)
      END IF
      IF (NMORDR .GE. 5) THEN
         NDIME = NDIME + 2*NCOOR
      END IF
      NDIME = NDIME*NINNER
C
      KDERIV = IWOFF + 2
      KENVAL = KDERIV +   NDERIV
      KCOOR  = KENVAL +   NDIME
      KCSTAR = KCOOR  + 3*NCOOR
      KSYMCO = KCSTAR + 3*NCOOR
      KLAST = KSYMCO +   NCOOR**2
C         
C     *** Memory needed for test on Hessian. ***
C
      IF (SDRTST) THEN
         KTSTGD = KLAST
         KTSTSD = KTSTGD + NCOOR
         KLAST  = KTSTSD + NCOOR**2
      END IF
      IF (PGMTST) THEN
         LTHTST = 2**NMORDR
         KENTST = KLAST
         KLAST  = KENTST + LTHTST
      END IF
C
C     *** Memory for screening of redundant force constants. ***
C
      LDPMTX = 0
      IF (NMORDR .GE. 4) THEN
         LDPMTX = LDPMTX + (NCOOR*(NCOOR+1)*(NCOOR+2)*(NCOOR+3))/24
      END IF
      IF (NMORDR .GE. 3) THEN
         LDPMTX = LDPMTX + (NCOOR*(NCOOR+1)*(NCOOR+2))/6
      END IF
      IF (NMORDR .GE. 2) THEN
         IFRSTD = 2**NMORDR 
         LDPMTX = LDPMTX + (NCOOR*(NCOOR+1))/2
         KDPMTX = KLAST
         KDCOEF = KDPMTX +   IFRSTD*NMORDR*LDPMTX
         KNIDPC = KDCOEF +   IFRSTD       *LDPMTX
         KLAST  = KNIDPC +                 LDPMTX
      END IF
C
C     *** Space for backup av MOLECULE.INP ***
C
      KMBKLN = KLAST
      KLAST  = KMBKLN + 80*NMLINE
C
C     *** Symmetry initialization and symmetry ***
C     ***       related memory allocation      ***
C
      CALL FCSINI
      KGRIRP = KLAST
      KCHRCT = KGRIRP + NGORDR*NGVERT
      KICRIR = KCHRCT + NGORDR*NCVERT
      KLAST  = KICRIR + 2*NCOOR
C
      LWRK1 = LWORK - (KLAST - IWOFF) + 1
C
      DISPLC = 1.0D-2
C
      IWSAV = IWOFF
      IWOFF = 0
      CALL NUMDR1(WORK(KDERIV),WORK(KENVAL),WORK(KCOOR),WORK(KCSTAR),
     &            WORK(KSYMCO),WORK(KDCOEF),WORK(KTSTGD),WORK(KTSTSD),
     &            WORK(KENTST),WORK(KGRIRP),WORK(KCHRCT),WORK(KLAST),
     &            WORK(KDPMTX),WORK(KNIDPC),WORK(KICRIR),WORK(KMBKLN),
     &            LWRK1,IWOFF,NDERIV,NDIME,MAXADR,LTHTST,LDPMTX,IFRSTD,
     &            WRKDLM)
      IWOFF = IWSAV
C
#if defined (VAR_MPI)
C
C     We let the slaves wait for the Master to tell them whether to pick up
C     a new geometry or to end this calculation
C
      IF (MYNUM .GT. 0) THEN
         CALL MPI_BCAST(NTASK,1,MPI_INTEGER,NODEID(MASTER),
     &                 MPI_COMM_WORLD,IERR)
         IF (NTASK .EQ. 1) THEN
            CALL PARION
            RDINPC = .FALSE.
            CALL READIN(WORK(IWOFF),LWORK,.FALSE.)
            GOTO 8
         ELSE IF (NTASK .EQ. 0) THEN
            CALL MPI_BCAST(FINISH,1,MPI_LOGICAL,NODEID(MASTER),
     &                  MPI_COMM_WORLD,IERR)
            CALL MPI_FINALIZE(IERR)
#if defined (SYS_CRAY)
            CALL INFO = ISHELL('rm -f $SCRATCHDIR/*')
#else
            CALL SYSTEM('rm -f $SCRATCHDIR/*')
#endif     
            STOP '*** End of DALTON calculation ***'
         ELSE
            WRITE (LUPRI,'(/A)') 'Unknown message received by slave'
            CALL QUIT('Slave received unknown message from master')
         END IF
      END IF
#endif
      RETURN
      END
C
C  /* Deck numdr1 */
      SUBROUTINE NUMDR1(DERIV,ENVAL,COOR,CSTART,SYMCOR,DCOEFF,TSTGDR,
     &                  TSTSDR,ENTST,GRIREP,CHRCTR,WORK,KDPMTX,NMIDPC,
     &                  ICRIRP,MBKLIN,LWORK,IWOFF,NDERIV,NDIME,MAXADR,
     &                  LTHTST,LDPMTX,IFRSTD,WRKDLM)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <dummy.h>
#include <mxorb.h>
      PARAMETER (D0 = 0.0D0)
#include <cbiwlk.h>
#include <cbiexc.h>
#include <cbinum.h>
#include <nuclei.h>
#include <symmet.h>
#include <exeinf.h>
#include <abainf.h>
#include <infpar.h>
#include <taymol.h>
#include <trkoor.h>
#include <numder.h>
#include <past.h>
#include <gnrinf.h>
#include <inftap.h>
#include <molinp.h>
#include <fcsym.h>
      LOGICAL EXHER, EXSIR, EXABA, GTCOMP, SYMDET, RSTDON
      CHARACTER*80 MBKLIN, PRTEXT
      CHARACTER*8 ANDER
      DIMENSION MBKLIN(NMLINE)
      DIMENSION DERIV(NDERIV), ENVAL(NDIME), CSTART(3*NCOOR), 
     &          COOR(3*NCOOR), SYMCOR(NCOOR,NCOOR), 
     &          DCOEFF(LDPMTX,IFRSTD), TSTGDR(NCOOR), ENTST(LTHTST), 
     &          TSTSDR(NCOOR,NCOOR), GRIREP(NGORDR,NGVERT), 
     &          CHRCTR(NGORDR,NCVERT),WORK(LWORK)
      DIMENSION ICRIRP(NCOOR,2), KDPMTX(LDPMTX,NMORDR,IFRSTD), 
     &          NMIDPC(LDPMTX)
#include <ibtfun.h>
C
C     *** Restart initialization ***
C
      RSTDON = .FALSE.
C
C     *************************************
C     ***Backup of original MOLECULE.INP***
C     ***   To finish off correctly.    ***
C     *************************************
C
      NMBKLN = NMLINE
      DO 100 IMLINE = 1, NMLINE
         MBKLIN(IMLINE) = MLINE(IMLINE)
 100  CONTINUE
C     
C     ****************************************
C     ***Unrolling the symmetry coordinates***
C     ***  In order to take proper steps   ***
C     ****************************************
C
      ICOOR = 0
      IATOM = 0
      DO 200 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         DO 300 IOP = 0, MAXOPR
            IF (IBTAND(IOP,MULCNT) .EQ. 0) THEN
               IATOM = IATOM + 1
               DO 400 I = 1, 3
                  ICOOR = ICOOR + 1
                  CSTART(ICOOR) =
     &                 PT(IBTAND(ISYMAX(I,1),IOP))*CORD(I,ICENT)
 400           CONTINUE
            END IF
 300     CONTINUE
 200  CONTINUE
C
C     *** Restart ***
C     
      LURSTR = 0
      CALL GPOPEN(LURSTR,'RSTRT.FC','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LURSTR)
C     
      IF (RESTRT) THEN
C
C        *** Restart, find which round it ended  ***
C        *** in this run of nmder. Reread enval  ***
C        *** values from file.                   ***
C     
         KEND = 0
         IDIMAX = 0
         IDIMIN = 2
         CALL RERSTR(ENVAL,SYMCOR,ICRIRP,NDIME,KEND,IDIMAX,IDIMIN,
     &               LURSTR,IPRINT,RSTDON)
      ELSE
C
C        *** Open restart file. A zero first in the file ***
C        *** means that the calculation ended here.      ***
C
         WRITE (LURSTR,'(I2)') 0
      END IF
C
      NUMCAL = 0
      SYMDET = .TRUE.
      NDCOOR = NCOOR
      IF (NRMCRD) THEN
         NTMPDR = NMORDR
         NMORDR = 2 - NAORDR
         IF ((NAORDR+NMORDR).LT.2) NDORDR = NMORDR
      END IF
      MAXINR = 2**NMORDR
C
      KIADRS = 1
      KINDST = KIADRS + MAXADR
      KNPRTN = KINDST + NMORDR
      KIDCMP = KNPRTN + MAXINR
      KIEQVG = KIDCMP + NCOOR
      KICIN  = KIEQVG + 2*NMORDR
      KIRPID = KICIN  +   NMORDR
#if defined (VAR_MPI)
      KFTVAL = KIRPID +   NMORDR
      KLAST  = KFTVAL +   NDIME
#else
      KLAST  = KIRPID +   NMORDR
#endif
C
      LWRK  = LWORK - KLAST + 1
C
      CALL NMDER(DERIV,ENVAL,COOR,CSTART,SYMCOR,ENTST,DCOEFF,GRIREP,
     &           CHRCTR,WORK(KLAST),ICRIRP,KDPMTX,NMIDPC,WORK(KIADRS),
     &           WORK(KINDST),WORK(KNPRTN),WORK(KIDCMP),WORK(KIEQVG),
     &           WORK(KICIN),WORK(KIRPID),MBKLIN,NMBKLN,NDERIV,LWRK,
     &           NDIME,MAXADR,LDPMTX,IFRSTD,MAXINR,IWOFF,LTHTST,IDIMAX,
     &           IDIMIN,LURSTR,WRKDLM,
#if defined (VAR_MPI)
     &              WORK(KFTVAL),SYMDET,RSTDON)
#else
     &              SYMDET,RSTDON)
#endif
      IF (NRMCRD) THEN
         NMORDR = NTMPDR
      END IF
C
C        *** Close the restart file ***
C
      CALL GPCLOSE(LURSTR,'KEEP')
C
C        *** Reevaluate restart parameter ***
C
      IF (RESTRT.AND.RSTDON) THEN
         RESTRT = .FALSE.
         RSTDON = .FALSE.
      END IF
C
C        *** Printing the derivatives ***
C
      IF (MYNUM .EQ. 0) THEN
         PRTEXT(1:9) = 'symmetric'
         LTEXT  = 9
         IF (NRMCRD) THEN
            NPRRDR = 2
         ELSE
            NPRRDR = NMORDR
         END IF
         KTDRS  = 1
         NDIMT  = NCOOR*(NCOOR+1)*(NCOOR+2)/6
         KFDRS  = NDIMT  + 1
         NDIMF  = NCOOR*(NCOOR+1)*(NCOOR+2)*(NCOOR+3)/24
         CALL PRDERV(DERIV(KTDRS),DERIV(KFDRS),TSTGDR,TSTSDR,SYMCOR,
     &               WORK,ICRIRP,LWORK,NPRRDR,NDIMT,NDIMF,LTEXT,IPRINT,
     &               PRTEXT)
C
         IF (SDRTST) THEN
            KTMPGD = 1
            KTMPHS = KTMPGD + MXCOOR
            KBKPGD = KTMPHS + NCOOR**2
            KBKPHS = KBKPGD + NCOOR
            KLAST  = KBKPHS + NCOOR**2
            LWRK = LWORK - KLAST
            IF (KLAST.GT.LWORK) CALL QUIT('Memory exceeded in SDERTT')
            CALL SDERTT(TSTSDR,TSTGDR,SYMCOR,WORK(KTMPGD),WORK(KTMPHS),
     &                  WORK(KBKPGD),WORK(KBKPHS),WORK(KLAST),LWRK,
     &                  IPRINT)
         END IF
C     
         IF ((NAORDR+NMORDR).GT.1) THEN
            IF (.NOT.DRYRUN) THEN
               KEIGNV = 1
               KEGNVC = KEIGNV + NCOOR
               KHSMWT = KEGNVC + NCOOR**2
               KMT1TP = KHSMWT + NCOOR**2
               KMT2TP = KMT1TP + NCOOR**2
               KAMASS = KMT2TP + NCOOR**2
               KDKIN  = KAMASS + NATOMS
               KHTSTM = KDKIN  + NCOOR
               KFREQ  = KHTSTM + NCOOR**2
               KNATYP = KFREQ  + NCOOR
               KNMSSP = KNATYP + NATOMS
               KLAST  = KNMSSP + NCOOR
               LWRK   = LWORK  - KLAST 
               IF (KLAST.GT.LWORK)CALL QUIT('Memory exceeded in MKNRMC')
               CALL MKNRMC(SYMCOR,CSTART,WORK(KEIGNV),WORK(KEGNVC),
     &                     WORK(KHSMWT),WORK(KMT1TP),WORK(KMT2TP),
     &                     WORK(KAMASS),WORK(KDKIN),WORK(KHTSTM),
     &                     WORK(KFREQ),WORK(KLAST),ICRIRP,WORK(KNATYP),
     &                     WORK(KNMSSP),LWRK)
            ELSE
               CALL DRNRMC(SYMCOR,ICRIRP,IPRINT)
            END IF
         END IF
      END IF
C
      IF (NRMCRD .AND. ((NAORDR+NMORDR).GT.2)) THEN
C
C        *** Restart ***
C
         LURSTR = 0
         CALL GPOPEN(LURSTR,'RSTRT.FC','UNKNOWN',' ','FORMATTED',
     &               IDUMMY,.FALSE.)
C
         IF (RESTRT) THEN
C
C           *** Restart, find which round it ended  ***
C           *** in this run of nmder. Reread enval  ***
C           *** values and normal coordinates from  ***
C           *** file.                               ***
C     
            KEND = 1
            IDIMAX = 0
            IDIMIN = 2
            CALL RERSTR(ENVAL,SYMCOR,ICRIRP,NDIME,KEND,IDIMAX,IDIMIN,
     &                  LURSTR,IPRINT,RSTDON)
         ELSE
C
C           *** Open restart file. A 1 first in the file ***
C           *** means that the calculation ended here.   ***
C           *** Normal coordinates are also written to   ***
C           *** file for restart purposes.               ***
C
            WRITE (LURSTR,'(I2)') 1
            CALL WRICOR(SYMCOR,ICRIRP,LURSTR,IPRINT)
         END IF
C
         MAXINR = 2**NMORDR
C
         KIADRS = 1
         KINDST = KIADRS + MAXADR
         KNPRTN = KINDST + NMORDR
         KIDCMP = KNPRTN + MAXINR
         KIEQVG = KIDCMP + NCOOR
         KICIN  = KIEQVG + 2*NMORDR
         KIRPID = KICIN   +   NMORDR
#if defined (VAR_MPI)
         KFTVAL = KIRPID +   NMORDR
         KLAST  = KFTVAL +   NDIME
#else
         KLAST  = KIRPID +   NMORDR
#endif
C
         LWRK   = LWORK - KLAST + 1
C     
         CALL NMDER(DERIV,ENVAL,COOR,CSTART,SYMCOR,ENTST,DCOEFF,GRIREP,
     &              CHRCTR,WORK(KLAST),ICRIRP,KDPMTX,NMIDPC,
     &              WORK(KIADRS),WORK(KINDST),WORK(KNPRTN),WORK(KIDCMP),
     &              WORK(KIEQVG),WORK(KICIN),WORK(KIRPID),MBKLIN,NMBKLN,
     &              NDERIV,LWRK,NDIME,MAXADR,LDPMTX,IFRSTD,MAXINR,IWOFF,
     &              LTHTST,WIDIMAX,IDIMIN,LURSTR,WRKDLM,
#if defined (VAR_MPI)
     &              WORK(KFTVAL),SYMDET,RSTDON)
#else
     &              SYMDET,RSTDON)
#endif
C
C        *** Close the restart file ***
C
         CALL GPCLOSE(LURSTR,'KEEP')
C
C        *** Reevaluate restart parameter ***
C
         IF (RESTRT.AND.RSTDON) THEN
            RESTRT = .FALSE.
         END IF
C     
         PRTEXT(1:6) = 'normal'
C
         LTEXT  = 6
         NPRRDR = NAORDR + NMORDR
         KTDRS  = 1
         KFDRS  = (NDCOOR*(NDCOOR+1)*(NDCOOR+2))/6  + 1
C
         NDIMT  = NDCOOR*(NDCOOR+1)*(NDCOOR+2)/6
         NDIMF  = NDCOOR*(NDCOOR+1)*(NDCOOR+2)*(NCOOR+3)/24
         CALL PRDERV(DERIV(KTDRS),DERIV(KFDRS),TSTGDR,TSTSDR,SYMCOR,
     &               WORK,ICRIRP,LWORK,NPRRDR,NDIMT,NDIMF,LTEXT,IPRINT,
     &               PRTEXT)
      END IF
C
      IF (NAORDR.EQ.0) THEN
         ANDER = 'energy  '
      ELSE IF (NAORDR .EQ. 1) THEN
         ANDER = 'gradient'
      ELSE
         ANDER = 'hessian '
      END IF
C         
      IF (DRYRUN) THEN
         WRITE (LUPRI,'(5X,A,I15)') 'Number of ' // ANDER // 
     &              ' calculations needed:', NUMCAL
      ELSE
         WRITE (LUPRI,'(5X,A,I15)') 'Number of ' // ANDER // 
     &        ' calculations done:', NUMCAL
      END IF
C
      RETURN 
      END
C
C  /* Deck nmder */
      SUBROUTINE NMDER(DERIV,FUVAL,COOR,CSTART,SYMCOR,ENTST,DCOEFF,
     &                 GRIREP,CHRCTR,WORK,ICRIRP,KDPMTX,NMIDPC,IADRSS,
     &                 INDSTP,NPRTNR,IDCOMP,IEQVGM,ICIN,IRPIND,MBKLIN,
     &                 NMBKLN,NDERIV,LWORK,NDIME,MAXADR,LDPMTX,IFRSTD,
     &                 MAXINR,IWOFF,LTHTST,IDIMAX,IDIMIN,LURSTR,WRKDLM,
#if defined (VAR_MPI)
     &                 FTVAL,SYMDET,RSTDON)
#else
     &                 SYMDET,RSTDON)
#endif
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
      PARAMETER (D1 = 1.0D0, DM1 = -1.0D0, D0 = 0.0D0, DMAX=1.0D-15) 
#if defined (VAR_MPI) 
      INCLUDE 'mpif.h'      
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      DIMENSION FTVAL(NDIME)
#endif
#if defined (VAR_MPI2)
C
C     pario.h will no longer be needed as an include file when locking RMA
C     operations become available.
C
#include <dummy.h>
#include <pario.h>
#endif
#include <infpar.h>
C      
#include <trkoor.h>
#include <nuclei.h>
#include <numder.h>
#include <molinp.h>
#include <cbiwlk.h>
#include <cbinum.h>
#include <optinf.h>
#include <taymol.h>
#include <fcsym.h>
C
      LOGICAL CALCMP, CLNRGY, LASTE, NOSYM, EXESIR, EXEHER, EXEABA, 
     &        TOTSYM, SYMDET, RSTDON, PRTNR, ALRCAL
      CHARACTER*80 MBKLIN
      DIMENSION MBKLIN(NMBKLN)
      DIMENSION FUVAL(NDIME), COOR(NCOOR), CSTART(NCOOR),WORK(LWORK), 
     &          DERIV(NDERIV), DCOEFF(LDPMTX,IFRSTD),
     &          SYMCOR(NCOOR,NCOOR), GRIREP(NGORDR,NGVERT), 
     &          CHRCTR(NGORDR,NCVERT), EGRAD(MXCOOR), ENTST(LTHTST) 
      DIMENSION INDSTP(NMORDR), ICIN(NMORDR), IDCOMP(NCOOR), 
     &          IADRSS(MAXADR), IRPIND(NMORDR), ICRIRP(NCOOR,2), 
     &          KDPMTX(LDPMTX,NMORDR,IFRSTD), NMIDPC(LDPMTX), 
     &          IEQVGM(NMORDR,2), NPRTNR(MAXINR)
#include <ibtfun.h>
C
C     *** Finite field header. ***
C
      CALL TITLER('Numerical derivatives.','*',118)
      CALL HEADER('Derivatives calculated:',0)
      WRITE (LUPRI,'(A,I3)') 
     &     '           Derivatives calculated to order', NMORDR + NAORDR
      WRITE (LUPRI,'(A,I3,A)') 
     &     '           Analytical derivatives from energies to ',
     &     NAORDR, ' order.'                  
      WRITE (LUPRI,'(A,I3,A,I3,A)') 
     &     '         ', NMORDR, '. numerical derivatives from',NAORDR,
     &     '. order analytical derivatives'
C
C     *** Symmetry adapted coordinates ***
C
      IF (SYMDET) THEN
         CALL GRPCHR(CSTART,SYMCOR,GRIREP,CHRCTR,WORK,ICRIRP,LWORK,
     &               IPRINT)
         SYMDET = .FALSE.
      END IF
C
C     *** Finding force constants that are dependent on each-other ***
C
      NLDPMX = 0
      KDIM = IFRSTD*NMORDR*LDPMTX
      CALL IZERO(KDPMTX,KDIM)
      CALL FSDCST(SYMCOR,GRIREP,DCOEFF,WORK,KDPMTX,NMIDPC,ICRIRP,LDPMTX,
     &            IFRSTD,NLDPMX,LWORK,IPRINT)
c      if (nmordr.gt.2) stop '2'
C     
C     *** Memory allocations for future use ***
C
      NTYPE  = 3
      LASTE  = .FALSE.
      NOSYM  = .FALSE.
      IF (MLINE(NMLAU)(10:10).EQ.'0') NOSYM = .TRUE.
      KNATOM = 1
      KNLINE = KNATOM + NUCIND
      KNATIN = KNLINE + NUCIND
C
      KIDTMP = 1
      KIRPDG = KIDTMP + NMORDR
      KIRPST = KIRPDG + NMORDR
C
      ITYPE = 1
      NINTIN = 1
      IF (NAORDR .GT. 0) NINTIN = NINTIN*NCOOR
      IF (NAORDR .GT. 1) NINTIN = NINTIN*(NCOOR+1)/2
C      
C     Commented out to make optimizations work
C
C      NOMOVE = .FALSE.
      IDIME = NINTIN + 1
      IF (PGMTST) EMAX = D0
      NUMCAL = 0
      CALL DZERO(FUVAL,NDIME)
#if defined (VAR_MPI)
      CALL DZERO(FTVAL,NDIME)
#endif
#if defined (VAR_MPI2)
      LUNMCL = -9056
#endif
C
C     *** Order for derivatives ***
C
      DO 100 IORDR  = 1, NMORDR
         IHORDR = INT((IORDR+1)/2)
C
C        *** IHORDR -> The maximum order in one direction for this ***
C        ***           numerical derivative.                       ***
C
         DO 200 IMXRDR = 1, IHORDR, 1
C
            ITYPE         = ITYPE + 1
            IADRSS(ITYPE) = IDIME - 1
C
C           ***   The first component IX1 has always the largest order,   ***
C           ***          and are then independent of the others           ***
C           *** The order of the other components are not larger than one ***
C
            IMINCR = 1
            IF (IMXRDR .EQ. 1) IMINCR = IORDR
            IRSRDR    = IORDR - (2*IMXRDR-1)
            DO 300 IX1    = IMINCR, NDCOOR       
C
C              *** Starting values for the component-vector. ***
C
               INDSTP(1) = IX1
               DO 400 IC = IRSRDR+1, 2, -1
                  INDSTP(IC) = IRSRDR+2-IC
 400           CONTINUE
               IF (IRSRDR .GT. 0) INDSTP(IRSRDR+1) = INDSTP(IRSRDR+1)-1
C
               NSTP = 1
               IF (IMXRDR .EQ. 1) THEN
                  DO 500 I = 1, IRSRDR
                     NSTP = NSTP*(IX1-I)/I
 500              CONTINUE
               ELSE
                  DO 600 I = 1, IRSRDR
                     NSTP = NSTP*(NDCOOR-I+1)/I
 600              CONTINUE
               END IF
C     
C              *** NSTP -> Number of components for this IORDR, IMXRDR and IX1***
C
               DO 700 ISTP = 1, NSTP
C
C                 *** Finding the other components. ***
C     
                  CALCMP = .TRUE.
                  IF (IMXRDR .EQ. 1) THEN
                     DO 800 IC = IRSRDR+1, 2, -1
                        IF (INDSTP(IC) .LT. (INDSTP(IC-1)-1)) THEN
                           INDSTP(IC) = INDSTP(IC) + 1
                           ITEMP = ITEMP - 1
                           DO 900 I = IC+1, IRSRDR+1
                              INDSTP(I) = (IRSRDR+2) - I
 900                       CONTINUE
                           GOTO 1300
                        END IF
 800                 CONTINUE
                  ELSE
                     DO 1000 IC = IRSRDR+1, 2, -1
                        IF (IC .EQ. 2) THEN
                           INDSTP(2) = INDSTP(2) + 1
                           DO 1100 I = 3, IRSRDR+1
                              INDSTP(I) = 1
 1100                      CONTINUE
                           DO 1150 ICN = 2, IRSRDR+1
                              IF (INDSTP(ICN) .EQ. IX1) CALCMP = .FALSE.
 1150                      CONTINUE
                           GOTO 1300
                        ELSE IF (INDSTP(IC) .LT. (INDSTP(IC-1)-1)) THEN
                           INDSTP(IC) = INDSTP(IC) + 1
                           DO 1200 I = IC+1, IRSRDR
                              INDSTP(I) = 1
 1200                      CONTINUE
                           DO 1250 ICN = 2, IRSRDR+1
                              IF (INDSTP(ICN) .EQ. IX1) CALCMP = .FALSE.
 1250                      CONTINUE
                           GOTO 1300
                        END IF
 1000                CONTINUE                        
                  END IF
C     
 1300             CONTINUE
C
C                 *** Have we calculated this function-value before? ***
C
                  IF (CALCMP) THEN
C
C                    *** IDCOMP(INDSTP(IC)) -> gives the length of the steps we need ***
C                    ***                       to do in INDSTP(IC) direction         ***
C
                     CALL IZERO(IDCOMP,NDCOOR)
                     IDCOMP(INDSTP(1)) = IMXRDR
                     DO 1400 IC = 2, IRSRDR+1
                        IDCOMP(INDSTP(IC)) = IDCOMP(INDSTP(IC)) + 1
 1400                CONTINUE
C     
C                    *** NINNER -> Number of different steps needed ***
C                    *** If numerical derivatives from energy is    ***
C                    *** calculated, we need to check whether the   ***
C                    *** steps are all totally symmetric.           ***
C
                     IF (NAORDR .EQ. 0) THEN
                        IJ = 1
                        TOTSYM = .FALSE.
                        DO IRDR = 1, IRSRDR+1
                           IJ = IJ*ICRIRP(INDSTP(IRDR),1)
                        END DO
                        IF (IJ.EQ.1) TOTSYM = .TRUE.
                     END IF
C
                     NINNER = 2**(IRSRDR+1)
C
                     NMPRTN = 0
                     CALL IZERO(NPRTNR,MAXINR)
C
                     DO 1500 IINNER = 1, NINNER
C     
C                       *** Finding the appropriate step-possibility. ***
C
                        IC   = 0
                        IMOD = 1
                        IDIV = 1
                        DO 1600 I = 1, IRSRDR+1
                           ICIN(I) = MOD(INT((IINNER-1)/IDIV),2) + 1 
                           IDIV = IDIV*2
 1600                   CONTINUE
C
                        KDIM = 3*NCOOR
                        CALL DCOPY(KDIM,CSTART,1,COOR,1)
C
C                       *** Making the appropriate step. ***
C
                        SCLFCK = D1
                        DO 1700 IC = 1, IRSRDR+1
                        DO 1700 IMXN = 1, IDCOMP(INDSTP(IC))
                           CALL STPCOR(COOR,COOR,SYMCOR,DISPLC,NCOOR,
     &                                 ICIN(IC),INDSTP(IC),IPRINT) 
 1700                   CONTINUE
C
C                       *** Symmetry of derivatives calculated ***
C
                        CLNRGY = .FALSE.
                        KIDTMP = 1
                        KIDDBT = KIDTMP + NMORDR
                        KIRPDG = KIDDBT + NMORDR
                        KIRPST = KIRPDG + NMORDR
                        KLAST  = KIRPST + NMORDR
                        LWRK   = LWORK  - KLAST + 1
                        CALL FCSCRN(GRIREP,WORK(KLAST),KDPMTX,INDSTP,
     &                              ICRIRP,IRPIND,WORK(KIDTMP),
     &                              WORK(KIDDBT),WORK(KIRPDG),
     &                              WORK(KIRPST),NPRTNR,LWRK,NLDPMX,
     &                              LDPMTX,IFRSTD,IORDR,IRSRDR,MAXINR,
     &                              IINNER,NMPRTN,IPRINT,CLNRGY,PRTNR,
     &                              ALRCAL,.FALSE.)
C
C                       *** Calculate the energy, gradient or ***
C                       *** hessian. If this is a paralell    ***
C                       *** job, we need to find the proper   ***
C                       *** processor.                        ***
C                       
                        IF (CLNRGY) THEN
C     
C                          *** Update MOLECULE.INP file. ***
C
#if defined (VAR_MPI)                           
#if defined (VAR_MPI2)
C
C     The following code should be replaced with RMA operations following
C     the MPI-2 standard and as exemplified in for example Fig.6-9 in
C     "Using MPI-2" by Gropp, Lusk, and Thakur, and WIN_LOCK becomes
C     available in an MPI-2 implementation
C
 199                       CONTINUE
                           OPEN(UNIT=99,FILE=WRKDIR(1:LENWRK)//'LOCK',
     &                          STATUS='NEW',FORM='FORMATTED',
     &                          ERR=199)
                           CALL GPOPEN(LUNMCL,
     &                                 WRKDIR(1:LENWRK)//'NUMCAL',
     &                                 'OLD',' ','FORMATTED',IDUMMY,
     &                                 .FALSE.)
                           READ (LUNMCL,'(I5)') NUMCL
                           NUMCAL = NUMCAL + 1
                           IF (NUMCL .EQ. NUMCAL) THEN
                              REWIND (LUNMCL)
                              WRITE (LUNMCL,'(I5)') NUMCL + 1
                              CALL GPCLOSE(LUNMCL,'KEEP')
                              CLOSE(UNIT=99,STATUS='DELETE',ERR=107)
 107                          CONTINUE
#else
                           NUMCAL = NUMCAL + 1
                           IF (MYNUM.EQ.MOD(NUMCAL,(NODTOT+1))) THEN
#endif
#endif
                              IF ((.NOT.DRYRUN).AND.((.NOT.RESTRT).OR.
     &                            (     RESTRT .AND.        RSTDON.AND.
     &                            (IDIME.GT.IDIMAX)))) THEN
                                 CALL UPCOMB(COOR,WORK(KNATOM),
     &                                       WORK(KNLINE),WORK(KNATIN),
     &                                       MBKLIN,NMBKLN,LASTE)
                              END IF
C
C                             *** Reset nescesary variables ***
C     
                              CALL NUMRES(NMORDR+NAORDR,EXHER,EXSIR,
     &                                    EXABA)
C     
C                             *** Find the energy, gradient or hessian ***
C     
                              IF (NAORDR .EQ. 0) THEN
                                IF ((.NOT.DRYRUN).AND.((.NOT.RESTRT).OR.
     &                              (     RESTRT .AND.       RSTDON.AND.
     &                                         (IDIME.GT.IDIMAX)))) THEN
                                   CALL GTNRGY(EXHER,EXSIR,EXABA,WORK,
     &                                         LWORK,IWOFF,WRKDLM)
#if defined (VAR_MPI)
                                   FTVAL(IDIME) = ENERGY
#else
                                   FUVAL(IDIME) = ENERGY
#endif
                                   WRITE(LURSTR,'(I8,F24.16)') IDIME,
     &                                                         ENERGY
                                   PRTNR = .FALSE.
                                   IF (PRTNR) THEN
                                      IMDIME = IDIME - IINNER 
     &                                               + NPRTNR(NMPRTN)
#if defined (VAR_MPI)                                      
                                      FTVAL(IMDIME) = ENERGY
#else
                                      FUVAL(IMDIME) = ENERGY
#endif
                                      WRITE(LURSTR,'(I8,F24.16)')IMDIME,
     &                                                           ENERGY
                                   END IF
                                   CALL FLSHFO(LURSTR)
                                END IF
                                IDIME = IDIME + 1
                             ELSE IF (NAORDR .EQ. 1) THEN
                                CALL GTGRAD(EGRAD,EXHER,EXSIR,EXABA,
     &                                      WORK,LWORK,IWOFF,WRKDLM)
                                DO ISCOOR = 1, NCOOR
#if defined (VAR_MPI)
                                   FTVAL(IDIME) = GRDMOL(ISCOOR)
#else
                                   FUVAL(IDIME) = GRDMOL(ISCOOR)
#endif
                                   IDIME = IDIME + 1
                                END DO
                             ELSE IF (NAORDR .EQ. 2) THEN
                                CALL GTHESS(DUMMY,FUVAL(IDIME),DUMMY,
     &                                      EXHER,EXSIR,EXABA,WORK,
     &                                      LWORK,IWOFF,WRKDLM)
                                DO ISCOO2 = 1, NCOOR
                                DO ISCOO1 = 1, ISCOO2
                                   IDIME = IDIME + 1
#if defined (VAR_MPI)
                                   FTVAL(IDIME) = HESMOL(ISCOO1,ISCOO2)
#else
                                   FUVAL(IDIME) = HESMOL(ISCOO1,ISCOO2)
#endif
                                END DO
                                END DO
                             END IF
#if defined (VAR_MPI)
                          ELSE
#if defined (VAR_MPI2)
                             CALL GPCLOSE(LUNMCL,'KEEP')
                             CLOSE(UNIT=99,STATUS='DELETE',ERR=108)
 108                         CONTINUE
#endif
                             DO INTIN = IDIME, IDIME+NINTIN
                                FUVAL(INTIN) = D0
                             END DO
                             IDIME = IDIME + NINTIN
                          END IF
#endif
                       ELSE 
                          IF (NAORDR .EQ. 0) THEN
                             ALRCAL = .FALSE.
                             IF (((.NOT.RESTRT).OR.(RESTRT.AND.
     &                             RSTDON.AND.(IDIME.GT.IDIMAX))).AND.
     &                            (.NOT.ALRCAL)) THEN
                                FUVAL(IDIME) = D0
                                WRITE(LURSTR,'(I8,F24.16)') IDIME, D0
                                CALL FLSHFO(LURSTR)
                             END IF
                             IDIME = IDIME + 1
                          ELSE
                             DO 2100 INTIN = IDIME, IDIME+NINTIN
                                FUVAL(INTIN) = D0
 2100                        CONTINUE
                             IDIME = IDIME + NINTIN
                          END IF
                       END IF
C     
 1500                CONTINUE
                  END IF
 700           CONTINUE
C     
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
C     *** Returning to the original geometry, only for master. ***
C
      IF (MYNUM .EQ. 0) THEN
C      
C        ***      Update MOLECULE.INP file.     ***
C     
         LASTE = .TRUE.
         CALL UPCOMB(CSTART,WORK(KNATOM),WORK(KNLINE),WORK(KNATIN),
     &               MBKLIN,NMBKLN,LASTE)
C
C        *** Reset nescesary variables ***
C
         CALL NUMRES(NMORDR+NAORDR,EXHER,EXSIR,EXABA)
C
C        *** Find the energy, gradient or hessian ***
C
         NUMCAL = NUMCAL + 1
         IF ((.NOT.DRYRUN).AND.((.NOT.RESTRT).OR.
     &                 (RSTDON.AND.RESTRT.AND.(IDIMIN.GT.1)))) THEN
            IF (NAORDR .EQ. 0) THEN
               CALL GTNRGY(EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,WRKDLM)
#if defined (VAR_MPI)
               FTVAL (1) = ENERGY
#else
               FUVAL (1) = ENERGY
#endif
               WRITE(LURSTR,'(I8,F24.16)') 1, ENERGY
            ELSE IF (NAORDR .EQ. 1) THEN
               CALL GTGRAD(EGRAD,EXHER,EXSIR,EXABA,WORK,LWORK,IWOFF,
     &                     WRKDLM)
C     
               DO 2200 I = 1, NCOOR
#if defined (VAR_MPI)
                  FTVAL(I) = EGRAD(I)
#else
                  FUVAL(I) = EGRAD(I)
#endif
 2200          CONTINUE
            ELSE IF (NAORDR .EQ. 2) THEN
               KHESS = 1
               KLAST = KHESS + MXCOOR**2
               CALL GTHESS(DUMMY,WORK(KHESS),DUMMY,EXHER,EXSIR,EXABA,
     &                     WORK,LWORK,IWOFF,WRKDLM)
C         
               KHINT = 0
               KFINT = 0
               DO 2300 J = 1, MXCOOR
               DO 2300 I = 1, J
                  KHINT = KHINT + 1
                  IF (WORK(KHESS+KHINT-1) .GT. DMAX) THEN
                     KFINT = KFINT + 1
#if defined (VAR_MPI)
                     FTVAL(KFINT) = WORK(KHESS+KHINT-1)
#else
                     FUVAL(KFINT) = WORK(KHESS+KHINT-1)
#endif
                  END IF
 2300          CONTINUE
            END IF
         END IF
         IADRSS(1) = 0
#if defined (VAR_MPI)
      ELSE
         DO INTIN = 1, NINTIN
            FUVAL(INTIN) = D0
         END DO
#endif         
      END IF
C
#if defined (VAR_MPI)
C
C     *** If paralell calculation all energies are ***
C     *** collected into one array.                ***
C
#if defined (VAR_MPI2)
C
C     However, in the case of "simulated" MPI2 behaviour (RMA operations),
C     there can occur a "glitch" in the NFS lock file, and points may happen
C     to be calculated on several processors.
C     
C     We collect results from one processor at a time, checking for double
C     counting
C
      IF (MYNUM .EQ. 0) THEN
         CALL DCOPY(NDIME,FTVAL,1,FUVAL,1)
         DO IWHO = 1, NODTOT
            CALL MPI_RECV(NWHO,1,MPI_INTEGER,MPI_ANY_SOURCE,65,
     &                    MPI_COMM_WORLD,ISTAT,IERR)
            CALL MPI_RECV(FTVAL,NDIME,MPI_DOUBLE_PRECISION,NWHO,65,
     &                    MPI_COMM_WORLD,ISTAT,IERR)
            DO IPOS = 1, NDIME
C
C     Molecular energies ought to be negative
C     
               IF (.NOT. (FUVAL(IPOS) .LT. D0))
     &              FUVAL(IPOS) = FTVAL(IPOS)
            END DO
         END DO
      ELSE
         CALL MPI_SEND(MYNUM,1,MPI_INTEGER,NODEID(MASTER),
     &                 65,MPI_COMM_WORLD,IERR)
         CALL MPI_SEND(FTVAL,NDIME,MPI_DOUBLE_PRECISION,NODEID(MASTER),
     &                 65,MPI_COMM_WORLD,IERR)
      END IF
#else
      CALL MPI_REDUCE(FTVAL,FUVAL,NDIME,MPI_DOUBLE_PRECISION,MPI_SUM,0,
     &                MPI_COMM_WORLD,IERR)
#endif
#endif      
C
C     *** Preliminary constants to derivatives ***
C
      IF (MYNUM .EQ. 0) THEN
         MXCOEF = INT(NMORDR/2) + 1
C
C     *** Calculating the derivatives ***
C
         KCOEF  = 1
         KIMAX  = KCOEF  + (2*MXCOEF+1)*(NMORDR+1)
         KIMIN  = KIMAX  +               NMORDR
         KICNT  = KIMIN  +               NMORDR
         KNCVAL = KICNT  +               NTYPE
         KIDDCP = KNCVAL +               NCOOR
         KLAST  = KIDDCP +               NCOOR
         CALL NMNDER(DERIV,WORK(KCOEF),FUVAL,GRIREP,WORK(KLAST),IADRSS,
     &               KDPMTX,ICRIRP,INDSTP,IDCOMP,WORK(KIMAX),
     &               WORK(KIMIN),WORK(KICNT),WORK(KNCVAL),WORK(KIDDCP),
     &               MXCOEF,NDIME,NTYPE,NDERIV,LDPMTX,IFRSTD,NLDPMX,
     &               LWORK)
C
C     *** Assigning values to the dependent force constants. ***
C     
         CALL ADDPFC(DERIV,DCOEFF,KDPMTX,NMIDPC,LDPMTX,IFIRSTD,NDERIV,
     &               NLDPMX,IPRINT)
C
         IF (PGMTST) THEN
            WRITE (LUPRI,'(/A)') 
     &           'Test "equal energy for partner '//
     &           'geometries" is complete.'
            WRITE (LUPRI,'(A,F20.12)') 'Maximum error in energy is:',
     &                                 EMAX
            WRITE (LUPRI,'(A,F20.12)') 'Relative error: '
            DO IRDR = 1, NMORDR
               WRITE (LUPRI,'(I2,A,F14.8)') IRDR, '. derivative: ', 
     &              EMAX/(DISPLC**(DBLE(IRDR)))
            END DO
         END IF
      END IF
C
      RETURN
      END
C     
C  /* Deck upcomb */
      SUBROUTINE UPCOMB(COOR,NATOM,NLINE,NATIND,MBKLIN,NMBKLN,LASTE)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
C
      PARAMETER (D100 = 100.0D0, THRSH = 1.0D-12)
#include <molinp.h>
#include <nuclei.h>
#include <trkoor.h> 
#include <symmet.h>
#include <inftap.h>
      LOGICAL BIG, USED, LASTE
      CHARACTER CHR
      CHARACTER*80 MBKLIN
      DIMENSION MBKLIN(NMBKLN)
      DIMENSION COOR(3,NCOOR/3), NATOM(NUCIND), NLINE(NUCIND), 
     &          NATIND(NUCIND)
#include <ibtfun.h>
C
C     Updates geometry in common block
C     
      IF (LASTE) THEN
         NMLINE = NMBKLN
         DO 100 IBKLIN = 1, NMBKLN
            MLINE(IBKLIN) = MBKLIN(IBKLIN)
 100     CONTINUE
      ELSE
         NADD   = 0
         IATOM  = 0
         NCLAST = 0
         IF (MLINE(NMLAU)(10:10).EQ.'0') THEN
            MLINE(NMLAU)(20:20) = ' '
         ELSE
            MLINE(NMLAU)(10:20) = ' '
         END IF
         DO 200 ICENT = 1, NUCIND
            ISYM   = 0
            NRLINE = NCLINE(ICENT)
            NC     = NCLINE(ICENT)
            MULCNT = ISTBNU(ICENT)
            IF (NC .NE. 0) THEN
               READ (MLINE(NC),9100) NAME
               DO 300 IOP = 0, MAXOPR
                  IF (IBTAND(IOP,MULCNT) .EQ. 0) THEN
                     IATOM = IATOM + 1
                     CRX = COOR(1,IATOM)
                     CRY = COOR(2,IATOM)
                     CRZ = COOR(3,IATOM)
                     BIG = (ABS(CRX) .GT. D100 .OR.
     *                      ABS(CRY) .GT. D100 .OR.
     *                      ABS(CRZ) .GT. D100)
                     IF (ISYM .GT. 0) THEN
                        DO 400 I = NMLINE, NC+1, -1
                           MLINE(I+1) = MLINE(I)                        
 400                    CONTINUE
                        DO 450 IC2 = ICENT+1, NUCIND
                           NCLINE(IC2) = NCLINE(IC2) + 1
 450                    CONTINUE
                        NRLINE = NRLINE + 1
                        NMLINE = NMLINE + 1
                        NC     = NC     + 1
                     END IF
                     IF (BIG) THEN
                        WRITE (MLINE(NC),9200) NAME,CRX,CRY,CRZ,'*'
                     ELSE
                        WRITE (MLINE(NC),9300) NAME,CRX,CRY,CRZ,'*'
                     END IF
                     ISYM = ISYM + 1
                  END IF
 300           CONTINUE
            END IF
 200     CONTINUE
C     
C        Do a count of each type of atom....
C
         KCENT2 = 0
         ICENT1 = 1
 500     CONTINUE
         ICENT1 = ICENT1 + KCENT2
         DO 700 ICENT2 = ICENT1+1, NUCIND
            READ (MLINE(NCLINE(ICENT2)-1),'(A)') CHR
            IF (CHR .EQ. ' ') THEN
               NAT = 0
               DO 800 IC = ICENT1, ICENT2-1
                  NAT = NAT + NUCDEG(IC)
 800           CONTINUE
               WRITE (MLINE(NCLINE(ICENT1)-1)(13:15),'(I3)') NAT
               KCENT2 = ICENT2-ICENT1
               GOTO 500
            ELSE
               IF (ICENT2 .EQ. NUCIND) THEN
                  KCENT1 = ICENT1
                  GOTO 900
               END IF
            END IF
 700     CONTINUE
         KCENT1 = NUCIND
 900     CONTINUE
         NAT = 0
         DO 1100 IC = KCENT1, NUCIND
            NAT = NAT + NUCDEG(IC)
 1100    CONTINUE
         WRITE (MLINE(NCLINE(KCENT1)-1)(13:15),'(I3)') NAT
      END IF
C     
C     Punch MOLECULE input with updated coordinates to LUMOL
C
      CALL GPOPEN(LUMOL,'MOLECULE.INP','OLD',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUMOL)
      DO 1300 IMLINE = 1,NMLINE
         WRITE (LUMOL,'(A)') MLINE(IMLINE)
 1300 CONTINUE
      CALL GPCLOSE(LUMOL,'KEEP')
C     
 9100 FORMAT (A4)
 9200 FORMAT (A4,3F20.10,7X,A1)
 9300 FORMAT (A4,3F20.16,7X,A1)
C     
      RETURN
      END
C
C  /*Deck stpcor*/
      SUBROUTINE STPCOR(COOR,CSTART,SYMCOR,DISPLC,NCOOR,KPM,KSCOOR,
     &                  IPRINT)
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION COOR(NCOOR), CSTART(NCOOR), SYMCOR(NCOOR,NCOOR)
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('CSTART in STPCOR',1)
         CALL OUTPUT(CSTART,1,1,1,NCOOR,1,NCOOR,1,LUPRI)
      END IF 
C 
      IF (KPM.EQ.1) THEN
        FAC =  DISPLC
      ELSE
        FAC = -DISPLC
      END IF

      DO 100 ICOOR = 1, NCOOR
         COOR(ICOOR) = CSTART(ICOOR) + FAC*SYMCOR(ICOOR,KSCOOR)
 100  CONTINUE
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('COOR in STPCOR',1)
         CALL OUTPUT(COOR,1,1,1,NCOOR,1,NCOOR,1,LUPRI)
      END IF
C
      RETURN
      END
C
C  /* Deck numres */
      SUBROUTINE NUMRES(NORDR,EXHER,EXSIR,EXABA)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
C      
#include <optinf.h>
#include <symmet.h>
#include <nuclei.h>
#include <gnrinf.h>
#include <huckel.h>
      LOGICAL EXHER,EXSIR,EXABA
C
C     This routine resets a few variables, to be able to calculate 
C     energy again for a new geometry (and symmetry). Only for 
C     derivatives higher then 1.
C
      EXHER  = .FALSE.
      EXSIR  = .FALSE.
      EXABA  = .FALSE.
      RDMLIN = .FALSE.
      RDINPC = .FALSE.
      IF (NORDR .GT. 1) THEN
         NEWSYM = .TRUE.
         ADDSTO = .TRUE.
         KEEPHE = .FALSE.
         RSTARR = .TRUE.
         BRKSYM = .FALSE.
         INDOLD = INDTOT
         GECONV = .FALSE.
         HRINPC = .FALSE.
         CALL IZERO(NUCNUM, MXCENT*8)
         CALL IZERO(NCRREP, 16)
         CALL IZERO(IPTCNT, MXCENT*48)
         CALL IZERO(NAXREP, 16)
         CALL IZERO(INDHES, 8)
      END IF
      RETURN
      END
C
C     
      SUBROUTINE PRNHST(KSCOO1,KSCOO2,DISPLC,I,J,NC)
C
C     Prints out headers to the different steps in finite field.
C
#include <implicit.h>
#include <priunit.h>
      CHARACTER SIGNS*4
C     
      SIGNS = '    '
C
      CALL TITLER('Finite field geometry','*',118)
C
      IF (NC .EQ. 1) THEN
         WRITE (LUPRI,'(A,I4)') '        Distortions in coordinate:',
     &                           KSCOO2
      ELSE
         WRITE (LUPRI,'(A,2I4)') '        Distortions in coordinates:',
     &                           KSCOO1, KSCOO2
      END IF
C
      WRITE (LUPRI,'(A)') '           '
      WRITE (LUPRI,'(A,F6.4)') '        Length of displacement:', DISPLC
      WRITE (LUPRI,'(A)') '           '
C
      IF (J .EQ. 1) THEN
         SIGNS(2:2) = '+'
      ELSE
         SIGNS(2:2) = '-'
      END IF
      IF (NC .NE. 1) THEN
         IF (I .EQ. 1) THEN
            SIGNS(4:4) = '+'
         ELSE
            SIGNS(4:4) = '-'
         END IF
      END IF
      WRITE (LUPRI,'(A24,A4)') '        Step directions:', SIGNS
C
      WRITE (LUPRI,'(A)') '           '
      WRITE (LUPRI,'(A)') '           '
      WRITE (LUPRI,'(A)') '           '
C     
      RETURN
      END
C /* Deck addnuc */
c      SUBROUTINE ADDNUC(WORK,LWORK)
c#include <implicit.h>
c#include <priunit.h>
c#include <mxcent.h>
cC
c      DIMENSION WORK(LWORK)
cC
c#include <energy.h>
c#include <taymol.h>
c#include <trkoor.h>
c#include <cbiwlk.h>
cC
cC     Calculate nuclear contributions
cC 
c      NDIM = 3*MXCOOR*MXCOOR
c      IF (LWORK .LT. NDIM) CALL STOPIT('ADDNUC',' ',LWORK,NDIM) 
c      CALL NUCREP(WORK,WORK(MXCOOR*MXCOOR+1),WORK(2*MXCOOR*MXCOOR+1))
C
C     Add nuclear contributions to gradient and Hessian
C
c      IF (NUMGRD) THEN
c         DO 100 I = 1, NCOOR
c            GRDMOL(I) = GRDMOL(I) + GRADNN(I)
c  100    CONTINUE
c      END IF
C     IF (NUMHES) THEN
C        DO 200 I = 1, NCOOR
C        DO 200 J = 1, NCOOR
C           HESMOL(J,I) = HESMOL(J,I) + WORK((I-1)*MXCOOR + J)
C 200    CONTINUE
C     END IF
c      RETURN
c      END
C /* Deck GTCOMP */
      LOGICAL FUNCTION GTCOMP()
#include <implicit.h>
#include <cbitro.h>
      GTCOMP = COMPAR
      RETURN
      END
C
C     
C /* Deck NMCOEF */
      SUBROUTINE NMCOEF(COEFF,TCOEFF,WORK,MXCOEF,NMNMDR,LWORK)
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.0D0, D05 = 0.5D0)
C
      DIMENSION COEFF (-MXCOEF:MXCOEF,0:NMNMDR), 
     &          TCOEFF(-NMNMDR:NMNMDR,0:NMNMDR)
C
      KDIM1 = (2*MXCOEF+1)*(NMNMDR+1)
      KDIM2 = (2*NMNMDR+1)*(NMNMDR+1)
      CALL DZERO(COEFF ,KDIM1)
      CALL DZERO(TCOEFF,KDIM2)
      COEFF(0,0) = D1
C
C     *** Temporary coeffesients used to generate coeffisients ***
C     ***            for even-numbered derivatives.            ***
C
      IF (NMNMDR .GT. 1) THEN
         NCOR = 1
         TCOEFF(0,0) = D1
         DO 100 IDR  = 1, NMNMDR
            DO 200 ICOR = -IDR+1, IDR-1
               TCOEFF(ICOR,IDR) = TCOEFF(ICOR-1,IDR-1) 
     &                          - TCOEFF(ICOR+1,IDR-1)
 200        CONTINUE
C
            TCOEFF(-IDR,IDR) = (-D1)**IDR
            TCOEFF( IDR,IDR) =   D1
C
 100     CONTINUE
C
C        *** Coefficients for even-numbered derivatives ***
C
         COEFF(0,0) = D1
         NEVEN = INT(NMNMDR/2)
         DO 300 IEVEN = 1, NEVEN
            IDR = 2*IEVEN
            DO 400 IECOR = -IEVEN,IEVEN
               ICOR = 2*IECOR
               COEFF(IECOR,IDR) = TCOEFF(ICOR,IDR)
 400        CONTINUE
 300     CONTINUE
      END IF
C
C     *** Coefficients for odd-numbered derivatives ***
C
      NODD = INT((NMNMDR+1)/2)
      DO 500 IODD = 1, NODD
         IDR = 2*IODD - 1
        DO 600 IOCOR = -IODD+1, IODD-1
            IF (IOCOR .NE. 0) THEN
               COEFF(IOCOR,IDR) = D05*(COEFF(IOCOR-1,IDR-1) 
     &                               - COEFF(IOCOR+1,IDR-1))
            END IF
 600     CONTINUE
C 
         COEFF(-IODD,IDR) = -D05
         COEFF( IODD,IDR) =  D05
C
 500  CONTINUE
C
      RETURN 
      END
C
C
C /* Deck NMNDER */
      SUBROUTINE NMNDER(DERIV,COEFF,FUNVAL,GRIREP,WORK,IADRSS,KDPMTX,
     &                  ICRIRP,INDSTP,IDCOMP,IMAX,IMIN,ICNT,NCVAL,
     &                  IDDCMP,MXCOEF,NDIME,NTYPE,NDERIV,LDPMTX,IFRSTD,
     &                  NLDPMX,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
      PARAMETER (DMIN = 1.0D-12, D1=1.0D0, D0=0.0D0)
#include <taymol.h>
#include <trkoor.h>
#include <cbiwlk.h>
#include <nuclei.h>
#include <symmet.h>
#include <numder.h>
#include <fcsym.h> 
#include <dummy.h>      
      LOGICAL CLFVAL
      DIMENSION COEFF(-MXCOEF:MXCOEF,0:NMORDR), DERIV(NDERIV), 
     &          FUNVAL(NDIME), GRIREP(NGORDR,NGVERT), WORK(LWORK)
      DIMENSION ICNT(NTYPE), IADRSS(NTYPE), IMAX(NMORDR), IMIN(NMORDR), 
     &          INDSTP(NMORDR), IDCOMP(NCOOR), IDDCMP(NCOOR), 
     &          NCVAL(NCOOR), KDPMTX(LDPMTX,NMORDR,IFRSTD), 
     &          ICRIRP(NCOOR,2)
#include <ibtfun.h>
C
      KDIM = MXCOOR**2
      IF (NAORDR .LT. 1) THEN
         CALL DZERO(GRDMOL,MXCOOR)
      END IF
      IF (NAORDR .LT.2 ) THEN
         CALL DZERO(HESMOL,  KDIM)
      END IF
      CALL DZERO(DERIV ,NDERIV)
C
C     ***************************************************
C     *** This subroutine calculates the coeffisients ***
C     ***      for the numerical differentiation.     ***
C     ***************************************************
C
      KTCOEF = 1
      KLAST  = KTCOEF + (2*NMORDR+1)*(NMORDR+1)
      CALL NMCOEF(COEFF,WORK(KTCOEF),WORK,MXCOEF,NMORDR,LWORK)
C
      IDERIV = 0
      DO 100 IORDR = 1, NMORDR
C
         MAXIND = 3*NUCIND
         IF ((IORDR+NAORDR) .EQ. 1) THEN
C
C           ****************************************
C           *** Numerical gradient from energies ***
C           ****************************************
C
            IMAX(1)  =  1
            IMIN(1)  = -1
            IDIME    =  1
            HINV     =  D1/DISPLC
            DO 200 ICOOR = 1, MAXIND
               IF (ICOOR .NE. 0) THEN
                  DO 250 I   = IMAX(1), IMIN(1), -2
                     IF (COEFF(I,1)**2 .GT. DMIN) THEN
                        IDIME = IDIME + 1
                        GRDMOL(ICOOR) = GRDMOL(ICOOR) 
     &                                + COEFF(I,1)*FUNVAL(IDIME)*HINV
                     END IF
 250              CONTINUE
                  sclfck = 1.0d0
                  GRDMOL(ICOOR) = SCLFCK*GRDMOL(ICOOR)
               END IF
 200        CONTINUE
         ELSE IF ((IORDR+NAORDR) .EQ. 2) THEN
C
C           ******************************
C           *** Numerical Hessian from *** 
C           ***  analytical gradients  ***
C           ******************************
C
            IF (NAORDR .EQ. 1) THEN
               IDIME    =  7
               HINV     =  D1/DISPLC
               DO 260 ICOOR2 = 1, MAXIND
                  IF (ICOOR2 .GT. 0) THEN
                     DO 270 ICOOR1 = 1, NCOOR
                        HESMOL(ICOOR1,ICOOR2) =
     &                            (COEFF( 1,1)*FUNVAL(IDIME      )
     &                           + COEFF(-1,1)*FUNVAL(IDIME+NCOOR))*HINV
                        IDIME = IDIME + 1
 270                 CONTINUE
                     IDIME = IDIME + NCOOR
                  END IF
                  sclfck = 1.0d0
                  GRDMOL(ICOOR2) = SCLFCK*GRDMOL(ICOOR2)
 260           CONTINUE
            ELSE
C
C           ******************************
C           *** Numerical Hessian from *** 
C           ***        energies.       ***
C           ******************************
C
               CALL IZERO(ICNT,NTYPE)
C
               IMAX(1)  =  1
               IMIN(1)  = -1
               HINV    =  D1/(DISPLC**2)
               DO 300 IX2 = 1, NDCOOR
               DO 300 IX1 = 1, IX2
                  IF (IX1 .LE. IX2) THEN
                     MX = IMIN(1)
                     IF (IX1 .EQ. IX2) MX = IMAX(1)
                     DO 400 I2 = IMAX(1), IMIN(1),-1
                     DO 400 I1 = IMAX(1), MX, -1
                        IF (IX1 .EQ. IX2) THEN 
                           IF (I2.EQ.0) THEN
                              ITYPE       = 1
                              ICNT(ITYPE) = 0 
                              BCOEFF      = COEFF(0,2)
                           ELSE
                              ITYPE       = 2
                              BCOEFF      = COEFF(I2,2)
                           END IF
                        ELSE 
                           ITYPE       = 3
                           BCOEFF      = COEFF(I1,1)*COEFF(I2,1)
                        END IF
C
                        IF (BCOEFF**2 .GT. DMIN) THEN
                           ICNT(ITYPE) = ICNT(ITYPE) + 1
                           HESMOL(IX2,IX1) = HESMOL(IX2,IX1) 
     &                               + BCOEFF*HINV
     &                               * FUNVAL(IADRSS(ITYPE)+ICNT(ITYPE))
                        END IF
 400                 CONTINUE
                     sclfck = 1.0d0
                     HESMOL(IX2,IX1) = SCLFCK*HESMOL(IX2,IX1)
                  END IF
 300           CONTINUE
            END IF
         ELSE
C
C          *********************************************
C          *** Numerical N'th derivative from NAORDR ***
C          ***         analytival derivative         ***
C          *********************************************
C
            NINNR2 = 0
C
            NSTP = 1
            DO 500 I = 1, IORDR
               NSTP = NSTP*(NDCOOR+I-1)/I
 500        CONTINUE
            POWER  = DBLE(IORDR)
            DIVDIS = D1/(DISPLC**POWER)
C
            CALL IZERO(INDSTP,NMORDR)
C
C           *** NSTP -> Number of components in the ***
C           ***      numerical differentiation      ***
C
            DO 600 ISTP = 1, NSTP
C
               
 700           CONTINUE
C
C              *** Finding which component this is ***
C
               DO 800 IC = IORDR, 1, -1
                  IF (IC .EQ. 1) THEN
                     INDSTP(1) = INDSTP(1) + 1
                     DO 900 I = 2, IORDR
                        INDSTP(I) = 1
 900                 CONTINUE
                     GOTO 1100
                  ELSE IF (INDSTP(IC) .LE. INDSTP(IC-1)-1) THEN
                     DO 1000 I = IC+1, IORDR
                        INDSTP(I) = 1
 1000                CONTINUE
                     INDSTP(IC) = INDSTP(IC) + 1
                     GOTO 1100
                  END IF
 800           CONTINUE
C
 1100          CONTINUE
C
C              *** Checking whether this component should be something. ***
C
               KIDTMP = 1
               KIDDBT = KIDTMP + NMORDR
               KIRPDG = KIDDBT + NMORDR
               KIRPST = KIRPDG + NMORDR
               KLAST  = KIRPST + NMORDR
               LWRK   = LWORK - KLAST + 1
               CALL FCSCRN(GRIREP,WORK(KLAST),KDPMTX,INDSTP,ICRIRP,
     &                     IDUMMY,WORK(KIDTMP),WORK(KIDDBT),
     &                     WORK(KIRPDG),WORK(KIRPST),IDUMMY,LWORK,
     &                     NLDPMX,LDPMTX,IFRSTD,IORDR,IORDR-1,IDUMMY,
     &                     IDUMMY,IDUMMY,IPRINT,CLFVAL,.FALSE.,.FALSE.,
     &                     .TRUE.)
C
               IF (CLFVAL) THEN
                  CALL IZERO(NCVAL,NDCOOR)
                  DO 1200 IC = 1, IORDR
                     NCVAL(INDSTP(IC)) = NCVAL(INDSTP(IC)) + 1
 1200             CONTINUE
C
C                 *** IDCOMP -> Maks steporder to get the derivative             ***
C                 *** IDDCMP -> Counting array, maks to min steporder            ***
C                 *** NTTYPE -> Number of function values needed for 1 component ***
C
                  NTTYPE = 1
                  CALL IZERO(IDCOMP,NDCOOR)
                  DO 1300 IC = 1, NDCOOR
                     IF (NCVAL(IC) .NE. 0) THEN
                        IDCOMP(IC) =   INT((NCVAL(IC)+1)/2)
                        NTTYPE     =   NTTYPE*(2*IDCOMP(IC) + 1)
                     END IF
 1300             CONTINUE
                  CALL ICOPY(NDCOOR,IDCOMP,1,IDDCMP,1)
C     
                  NINNR1 = 1
                  NINNR2 = 1
                  IF (NAORDR .EQ. 1) THEN
                     NINNR1 = NDCOOR
                     NINNR2 = INDSTP(IORDR)
                  END IF
                  IF (NAORDR .EQ. 2) THEN
                     NINNR1 = NINNR1*(NDCOOR+1)/2
                     NINNR2 = NINNR2*(INDSTP(IORDR)+1)/2
                  END IF
C
                  DO 1400 ITTYPE = 1, NTTYPE
C     
C                    *** Finding the right indices to identify ***
C                    ***        the right function value       ***
C     
                     DO 1500 IC = 1, NDCOOR
                        IF ((IDDCMP(IC) .GT. -IDCOMP(IC)) 
     &                              .AND. (ITTYPE .NE. 1)) THEN
                           IDDCMP(IC) = IDDCMP(IC) - 1
                           DO 1600 ICT = 1, IC-1
                              IDDCMP(ICT) = IDCOMP(ICT)
 1600                      CONTINUE
                           GOTO 1700
                        END IF
 1500                CONTINUE
 1700                CONTINUE
C     
C                    *** Calculate the coeffisient for this function value ***
C
                     BCOEFF = D1
                     NUMCOF = 0
                     DO 1800 IC = 1, NDCOOR
                        IF (NCVAL(IC) .NE. 0) THEN
                           BCOEFF = BCOEFF*COEFF(IDDCMP(IC),NCVAL(IC))
                           NUMCOF = NUMCOF + 1
                        END IF
 1800                CONTINUE
                     IF (NUMCOF .EQ. 0) BCOEFF = D0
C
C                    *** Does the function value contribute? ***
C
                     IF (BCOEFF**2 .GT. DMIN) THEN
C
C                       ************************************************
C                       *** This subroutine finds the adress for the ***
C                       ***      function value, from the indices    ***
C                       ***   NEIND - The adress in the FUVAL-array  **
C                       ************************************************
C
                        KITCMP = 1
                        CALL GTEIND(IADRSS,IDDCMP,IDCOMP,INDSTP,NCVAL,
     &                              WORK(KITCMP),NEIND,ITTYPE,IORDR)
C
C                       *** The derivative is calculated. ***
C
                        DO 1810 I = 1, NINNR2
                           ID    = IDERIV + I
                           DERIV(ID) = DERIV(ID) 
     &                               + BCOEFF*FUNVAL(NEIND)
 1810                   CONTINUE
C     
                     END IF
 1400             CONTINUE
C
C              *** The derivative is correctly scaled. *** 
C
                  SCLFCK = D1
                  DO 1900 I = 1, NINNR2
                     DERIV(IDERIV+I) = DERIV(IDERIV+I)*DIVDIS*SCLFCK
 1900             CONTINUE
               END IF
               IDERIV = IDERIV + NINNR2
 600        CONTINUE
         END IF
C
 100  CONTINUE
C
      RETURN 
      END
C
C
C   /*Deck gteind*/
      SUBROUTINE GTEIND(IADRSS,INDCMP,IDCOMP,INDSTP,NCVAL,ITCMP,NEIND,
     &                  NMTYPE,IORDR)
C     
C     *************************************************************************
C     *** This routine finds the adress of the function-value (in numerical ***
C     ***                differentiation) and returns it.                   ***
C     *** NEIND  -> The adress of the function value.                       ***
C     *** IADRSS -> Adress of the start of the steporder                    ***
C     *** NSTP   -> Number of function-component within the order           ***
C     *** NDIST  -> Number value within the component.                      ***
C     *************************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
      DIMENSION IADRSS(NMTYPE), INDCMP(NCOOR), INDSTP(NMORDR), 
     &          ITCMP(NMORDR), NCVAL(NCOOR)
C
      CALL IZERO(ITCMP,NMORDR)
C
C     *** Preliminary components ***
C
      IRSRDR = 0
      NMX = 0
      ITOT = 0
      ITWOTT = 0
      DO 100 I = NDCOOR, 1, -1
         IF (INDCMP(I) .NE. 0) THEN
            IRSRDR = IRSRDR + 1
            ITOT = ITOT + ABS(INDCMP(I))
            ITWOTT = ITWOTT + (2*ABS(INDCMP(I))-1)
            ITCMP(IRSRDR) = I
            IF (ABS(INDCMP(I)) .GT. NMX) THEN
               NMX   = ABS(INDCMP(I))
               NICMX = I
            END IF
         END IF
 100  CONTINUE
      NRSRDR = IRSRDR
C
      ITYPE = 1
      DO 300 IRDR  = 1, ITWOTT-1
         IHORDR = INT((IRDR+1)/2)
         DO 400 IMXRDR = 1, IHORDR
            ITYPE = ITYPE + 1
 400     CONTINUE
 300  CONTINUE
      NTYPE  = ITYPE + NMX
C
      IF (NMX .EQ. 0) THEN
         ISTP = 1
C
      ELSE IF (NMX .EQ. 1) THEN
         ISTP = 0
         DO 500 IRS1 = 1, NRSRDR - 1
            ITISTP = 1
            DO 600 I  = 1, NRSRDR-IRS1+1
               ITISTP = ITISTP*(ITCMP(IRS1)-I)/I
 600        CONTINUE
            ISTP = ISTP + ITISTP
 500     CONTINUE
         ISTP = ISTP + ITCMP(NRSRDR)
C
C
C
      ELSE
         IF (ITOT .EQ. NMX) THEN
            ISTP = ITCMP(1)
         ELSE
            IF (ITCMP(1) .NE. NICMX) THEN
               ITMP1 = ITCMP(1)
               ITCMP(1) = NICMX
               DO 650 I = 2, NRSRDR
                  ITMP2    = ITCMP(I)
                  ITCMP(I) = ITMP1
                  ITMP1    = ITMP2
 650           CONTINUE
            END IF
C
            ISTP = 0
            ISTP = (ITCMP(1)-1)*(NDCOOR-1) + ITCMP(2)
            IF (NICMX .LT. ITCMP(2)) ISTP = ISTP - 1
C     
            DO 700 IRS1 = 3, NRSRDR-1
               ITISTP = 1
               DO 800 I  = 1, NRSRDR-IRS1+1
                  ITISTP = ITISTP*(ITCMP(IRS1)-I)/I
 800           CONTINUE
               ISTP = ISTP + ITISTP
 700        CONTINUE
            IF (NRSRDR .GT. 2) ISTP = ISTP + ITCMP(NRSRDR)
         END IF
      END IF
      NINNER = 2**NRSRDR
      NSTP = (ISTP-1)*NINNER
C
      NDISP = 1
      IORD = NRSRDR-1
      DO 900 IRS = NRSRDR, 2, -1
         IDISP = 0
         IF (INDCMP(ITCMP(IRS)) .LT. 0) THEN
            IDISP = 2**IORD
         END IF
         IORD = IORD - 1
         NDISP = NDISP + IDISP
 900  CONTINUE
      IF (ITCMP(1) .NE. 0) THEN
         IF (INDCMP(ITCMP(1)) .LE. 0) NDISP = NDISP + 1
      END IF
C
      NEIND = IADRSS(NTYPE) + NSTP + NDISP
C
      RETURN
      END
C
C    /*Deck pritdr*/  
      SUBROUTINE PRITDR(TMPTDR,SYMCOR,TDER,SYMTDR,SCTRA,CSTRA,NTMPDM,
     &                  LTXT,TEXT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
C
#include <trkoor.h>
#include <symmet.h>
#include <numder.h>
#include <fcsym.h>
      character*80 TEXT
      DIMENSION TDER(NCOOR,NCOOR,NCOOR), TMPTDR(NTMPDM), 
     &          SYMTDR(NCOOR,NCOOR,NCOOR), CSTRA(NCOOR,NCOOR),
     &          SCTRA(NCOOR,NCOOR), SYMCOR(NCOOR,NCOOR)
C
C     *****************************************
C     *** Assigning values to proper places ***
C     ***  according to permutational sym.  ***
C     *****************************************
C
      ITMP = 0
      DO 100 K = 1, NDCOOR
      DO 100 J = 1, K
      DO 100 I = 1, J
         ITMP = ITMP + 1
         SYMTDR(I,J,K) = TMPTDR(ITMP)
         SYMTDR(I,K,J) = TMPTDR(ITMP)
         SYMTDR(J,I,K) = TMPTDR(ITMP)
         SYMTDR(J,K,I) = TMPTDR(ITMP)
         SYMTDR(K,I,J) = TMPTDR(ITMP)
         SYMTDR(K,J,I) = TMPTDR(ITMP)
 100  CONTINUE
C
C     **************************************
C     *** Printing symmetric coordinates ***
C     ***    Cartesian if no symmetry    ***
C     **************************************
C
      CALL HEADER('Third derivative of energy in ' // TEXT(1:LTXT) 
     &            // ' coordinates',-1)
C      
      IF (MOD(NDCOOR,6).EQ.0) THEN
         NLCMAX = NCOOR/6
      ELSE
         NLCMAX = INT(NCOOR/6)+1
      END IF
C
      DO 200 ICOL2 = 1, NDCOOR
         WRITE (LUPRI,'(A,I3)') '      Coloumn number', ICOL2
         WRITE (LUPRI,'(A)') '      -----------------'
         INLC = 0
         DO 400 INLCMX = 1, NLCMAX
            INLC2 = 6*(INLCMX-1) + 1
            INLC  = MIN(INLC+6,NDCOOR)
            DO 300 ICOL1 = 1, NDCOOR
               WRITE (LUPRI,'(A,6F15.6)') '         ', 
     &                           (SYMTDR(I,ICOL1,ICOL2), I=INLC2, INLC)
 300        CONTINUE
            WRITE (LUPRI,'(A)') '                              '
 400     CONTINUE
 200  CONTINUE
C
C     ************************************
C     *** Transformation to cartesian  ***
C     ***  coordinstes, and printing.  ***
C     ************************************
C
      IF ((FCLASS(1:3) .NE. 'C1 ').AND.(TEXT(1:LTXT).NE.'normal')) THEN
C
         CALL HEADER('Cartesian coordinates',-1)
C
         CALL TRATDR(CSTRA,SCTRA,SYMTDR,TDER,NDCOOR)
C
         IF (MOD(NCOOR,6).EQ.0) THEN
            NLCMAX = NCOOR/6
         ELSE
            NLCMAX = INT(NCOOR/6)+1
         END IF
C
         DO 500 ICOL2 = 1, NCOOR
            WRITE (LUPRI,'(A,I3)') '      Coloumn number', ICOL2
            WRITE (LUPRI,'(A)') '      -----------------'
            INLC = 0
            DO 600 INLCMX = 1, NLCMAX
               INLC2 = 6*(INLCMX-1) + 1
               INLC  = MIN(INLC+6,NCOOR)
               DO 700 ICOL1 = 1, NCOOR
                  WRITE (LUPRI,'(A,6F15.6)') '         ', 
     &                             (TDER(I,ICOL1,ICOL2), I=INLC2, INLC)
 700           CONTINUE
               WRITE (LUPRI,'(A)') '                              '
 600        CONTINUE
 500     CONTINUE
      END IF
C
      RETURN
      END
C
C    /*Deck tratdr*/  
      SUBROUTINE TRATDR(CSTRA,SCTRA,SYMTDR,TDER,NCOOR)
C     *********************************************
C     *** Transforming to cartesian coordinates ***
C     *********************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION CSTRA(NCOOR,NCOOR), SCTRA(NCOOR,NCOOR), 
     &          SYMTDR(NCOOR,NCOOR,NCOOR), TDER(NCOOR,NCOOR,NCOOR)
C
      CALL DZERO(TDER,NCOOR**3)
C
      ITYPE = 1
      CALL TRACOR(CSTRA,SCTRA,ITYPE,NCOOR,0)
C
      DO 100 ICART3 = 1, NCOOR
      DO 100 ICART2 = 1, NCOOR
      DO 100 ICART1 = 1, NCOOR
         DO 200 ISYCO3 = 1, NCOOR
         DO 200 ISYCO2 = 1, NCOOR
         DO 200 ISYCO1 = 1, NCOOR
            TDER(ICART1,ICART2,ICART3) = TDER(ICART1,ICART2,ICART3)
     &            + SCTRA(ICART1,ISYCO1)
     &             *SCTRA(ICART2,ISYCO2)*SCTRA(ICART3,ISYCO3)
     &             *SYMTDR(ISYCO1,ISYCO2,ISYCO3)
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C
C    /*Deck prifdr*/  
      SUBROUTINE PRIFDR(TMPFDR,FDER,SYMFDR,SCTRA,CSTRA,NTMPDM,LTXT,TEXT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
C
#include <trkoor.h>
#include <symmet.h>
#include <numder.h>
#include <fcsym.h>
      CHARACTER*80 TEXT
      DIMENSION FDER(NCOOR,NCOOR,NCOOR,NCOOR), TMPFDR(NTMPDM), 
     &          SYMFDR(NCOOR,NCOOR,NCOOR,NCOOR), SCTRA(NCOOR,NCOOR),
     &          CSTRA(NCOOR,NCOOR)
C
C     *****************************************
C     *** Assigning values to proper places ***
C     ***  according to permutational sym.  ***
C     *****************************************
C
      ITMP = 0
      DO 100 L = 1, NDCOOR
      DO 100 K = 1, L
      DO 100 J = 1, K
      DO 100 I = 1, J
         ITMP = ITMP + 1
C
         SYMFDR(I,J,K,L) = TMPFDR(ITMP)
         SYMFDR(I,J,L,K) = TMPFDR(ITMP)
         SYMFDR(I,K,J,L) = TMPFDR(ITMP)
         SYMFDR(I,K,L,J) = TMPFDR(ITMP)
         SYMFDR(I,L,J,K) = TMPFDR(ITMP)
         SYMFDR(I,L,K,J) = TMPFDR(ITMP)
C
         SYMFDR(J,I,K,L) = TMPFDR(ITMP)
         SYMFDR(J,I,L,K) = TMPFDR(ITMP)
         SYMFDR(J,K,I,L) = TMPFDR(ITMP)
         SYMFDR(J,K,L,I) = TMPFDR(ITMP)
         SYMFDR(J,L,I,K) = TMPFDR(ITMP)
         SYMFDR(J,L,K,I) = TMPFDR(ITMP)
C
         SYMFDR(K,I,J,L) = TMPFDR(ITMP)
         SYMFDR(K,I,L,J) = TMPFDR(ITMP)
         SYMFDR(K,J,I,L) = TMPFDR(ITMP)
         SYMFDR(K,J,L,I) = TMPFDR(ITMP)
         SYMFDR(K,L,I,J) = TMPFDR(ITMP)
         SYMFDR(K,L,J,I) = TMPFDR(ITMP)
C
         SYMFDR(L,I,J,L) = TMPFDR(ITMP)
         SYMFDR(L,I,L,J) = TMPFDR(ITMP)
         SYMFDR(L,J,I,L) = TMPFDR(ITMP)
         SYMFDR(L,J,L,I) = TMPFDR(ITMP)
         SYMFDR(L,K,I,J) = TMPFDR(ITMP)
         SYMFDR(L,K,J,I) = TMPFDR(ITMP)
 100  CONTINUE
C
C     **************************************
C     *** Printing symmetric coordinates ***
C     ***    Cartesian if no symmetry    ***
C     **************************************
C
      CALL HEADER('Fourth derivative of energy with respect to ' // 
     &             TEXT(1:LTXT) // ' nuclear coordinates',-1)
C         
      IF (MOD(NDCOOR,6).EQ.0) THEN
         NLCMAX = NDCOOR/6
      ELSE
         NLCMAX = INT(NDCOOR/6)+1
      END IF
C
      DO 200 ICOL3 = 1, NDCOOR
         WRITE (LUPRI,'(A,I4)') '      The fourth dimension', ICOL3
         WRITE (LUPRI,'(A)') '      ------------------------'
         WRITE (LUPRI,'(A)') '                                  '
         DO 300 ICOL2 = 1, NDCOOR
            WRITE (LUPRI,'(A,I4)') '        The third Dimension', ICOL2
            INLC = 0
            DO 400 INLCMX = 1, NLCMAX
               INLC2 = 6*(INLCMX-1) + 1
               INLC  = MIN(INLC+6,NDCOOR)
               DO 500 ICOL1 = 1, NDCOOR
                  WRITE (LUPRI,'(A,6F15.6)') '           ', 
     &                 (SYMFDR(I,ICOL1,ICOL2,ICOL3), I=INLC2, INLC)
 500           CONTINUE
               WRITE (LUPRI,'(A)') '                              '
 400        CONTINUE
 300     CONTINUE
 200  CONTINUE
C
C     ************************************
C     *** Transformation to cartesian  ***
C     ***  coordinstes, and printing.  ***
C     ************************************
C
      IF ((FCLASS(1:3).NE.'C1 ').AND.(TEXT(1:LTXT).NE.'normal')) THEN
C
         CALL HEADER('Cartesian coordinates',-1)
C
         CALL TRAFDR(CSTRA,SCTRA,SYMFDR,FDER,NCOOR)
C
         IF (MOD(NCOOR,6).EQ.0) THEN
            NLCMAX = NCOOR/6
         ELSE
            NLCMAX = INT(NCOOR/6)+1
         END IF
C     
         DO 600 ICOL3 = 1, NCOOR
            WRITE (LUPRI,'(A,I4)') '      The fourth dimension', ICOL3
            WRITE (LUPRI,'(A)') '      ------------------------'
            WRITE (LUPRI,'(A)') '                                  '
            DO 700 ICOL2 = 1, NCOOR
               WRITE (LUPRI,'(A,I4)')'        The third Dimension',ICOL2
               INLC = 0
               DO 800 INLCMX = 1, NLCMAX
                  INLC2 = 6*(INLCMX-1) + 1
                  INLC  = MIN(INLC+6,NCOOR)
                  DO 900 ICOL1 = 1, NCOOR
                     WRITE (LUPRI,'(A,6F15.6)') '           ', 
     &                    (SYMFDR(I,ICOL1,ICOL2,ICOL3), I=INLC2, INLC)
 900              CONTINUE
                  WRITE (LUPRI,'(A)') '                              '
 800           CONTINUE
 700        CONTINUE
 600     CONTINUE
      END IF
      RETURN
      END
C
C    /*Deck traFdr*/  
      SUBROUTINE TRAFDR(CSTRA,SCTRA,SYMFDR,FDER,NCOOR)
C     *********************************************
C     *** Transforming to cartesian coordinates ***
C     *********************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION CSTRA(NCOOR,NCOOR), SCTRA(NCOOR,NCOOR), 
     &          SYMFDR(NCOOR,NCOOR,NCOOR,NCOOR), 
     &          FDER(NCOOR,NCOOR,NCOOR,NCOOR)
C
      ITYPE = 1
      CALL TRACOR(CSTRA,SCTRA,ITYPE,NCOOR,0)
C
      DO 100 ICART4 = 1, NCOOR
      DO 100 ICART3 = 1, NCOOR
      DO 100 ICART2 = 1, NCOOR
      DO 100 ICART1 = 1, NCOOR
         DO 200 ISYCO4 = 1, NCOOR
         DO 200 ISYCO3 = 1, NCOOR
         DO 200 ISYCO2 = 1, NCOOR
         DO 200 ISYCO1 = 1, NCOOR
            FDER(ICART1,ICART2,ICART3,ICART4) 
     &           =  FDER(ICART1,ICART2,ICART3,ICART4) 
     *           +  SCTRA(ICART1,ISYCO1)*SCTRA(ICART2,ISYCO2)
     &             *SCTRA(ICART3,ISYCO3)*SCTRA(ICART4,ISYCO4)
     &                   *SYMFDR(ISYCO1,ISYCO2,ISYCO3,ISYCO4)
 200     CONTINUE
 100  CONTINUE
C
      RETURN
      END
C
C
C     /* Deck mknrmc */
      SUBROUTINE MKNRMC(SYMCOR,CSTART,EIGNVL,EGNVCT,HESMWT,TM1TMP,
     &                  TM2TMP,AMASS,DKIN,HTESTM,FREQ,WORK,ICRIRP,
     &                  NATTYP,NMSYSP,LWORK)
************************************************************
*** Makes normal coordinates from the molecular hessian, ***
*** and writes out the  harmonic frequensies according   ***
***             general symmetry species.                ***
************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <codata.h>
      PARAMETER (DMTHR = 1.0D-6, D0 = 0.0D0, D1 = 1.0D0)
#include <taymol.h>
#include <trkoor.h>      
#include <nuclei.h>
#include <cbinum.h>
#include <numder.h>
      DIMENSION SYMCOR(NCOOR,NCOOR), EIGNVL(NCOOR), EGNVCT(NCOOR,NCOOR), 
     &          HESMWT(NCOOR*(NCOOR+1)/2), AMASS(NATOMS), CSTART(NCOOR), 
     &          DKIN(NCOOR), TM1TMP(NCOOR,NCOOR), 
     &          TM2TMP(NCOOR,NCOOR), 
     &          FREQ(NCOOR), HTESTM(NCOOR,NCOOR), WORK(LWORK)
      DIMENSION ICRIRP(NCOOR,2), NATTYP(NATOMS), ISOTOP(MXCENT), 
     &          NMSYSP(NCOOR)
C
      DO 100 I = 1, NATOMS
         ISOTOP(I) = 1
 100  CONTINUE
C
C     *** Calculating center of mass, and mass of each center ***
C
      CALL CMMASS(CSTART,AMASS,NATTYP,ISOTOP,IPRINT)
C
C     *** Diagonal sqrt(mass)^(-1/2) matrix ***
C
      DO 200 IC = 1, NCOOR
         DKIN(IC) = D1/SQRT(XFAMU*AMASS((IC+2)/3))
 200  CONTINUE
C
C     *** The (mass)^(-1/2) matrix for symmetry coordinates. ***
C
      DO 300 IC2 = 1, NCOOR
      DO 300 IC1 = 1, NCOOR
         TM1TMP(IC1,IC2) = DKIN(IC1)*SYMCOR(IC1,IC2)
 300  CONTINUE
C
C     *** TM2TMP is the (mass)^(-1/2) matrix. ***
C
      KDIM = NCOOR**2
      CALL DZERO(TM2TMP,KDIM)
      DO 400 IC3 = 1, NCOOR
      DO 400 IC2 = 1, NCOOR
      DO 400 IC1 = 1, NCOOR
         TM2TMP(IC1,IC3) = TM2TMP(IC1,IC3)  
     &                   + SYMCOR(IC2,IC1)*TM1TMP(IC2,IC3)
 400  CONTINUE
C     
C     *** Test if TM2TMP is a diagonal matrix (if there are ***
C     ***              different isotopes)                  ***
C
      IF (MDGTST) THEN
         DO 500 IC2 = 1, NCOOR
         DO 500 IC1 = 1, NCOOR
            IF ((IC1 .NE. IC2).AND.(ABS(TM2TMP(IC1,IC2)).GT.DMTHR)) 
     &            CALL QUIT('Diagonal mass test failed. Off-diagonal' //
     &                      'elements present.')
 500     CONTINUE
      END IF
C     
C     *** Calculating the mass-weighted Hessian ***
C
      CALL DZERO(TM1TMP,KDIM)
      DO 600 IC3 = 1, NCOOR
      DO 600 IC2 = 1, NCOOR
      DO 600 IC1 = 1, NCOOR
         TM1TMP(IC1,IC3) = TM1TMP(IC1,IC3) 
     &                   + TM2TMP(IC1,IC2)*HESMOL(IC2,IC3)
 600  CONTINUE
C
C     *** HESMWT is the mass-weighted hessian. ***
C
      IC12 = 0
      CALL DZERO(HESMWT,KDIM)
      DO 700 IC2 = 1, NCOOR
      DO 700 IC1 = 1, IC2
         IC12 = IC12 + 1
         HESMWT(IC12) = HESMWT(IC12) + TM1TMP(IC1,IC2)*TM2TMP(IC2,IC2)
 700  CONTINUE
C
C     *** Test to check if mass-weighted hessian is symmetric. ***
C
      IF (HTEST) THEN
         CALL DZERO(HTESTM,KDIM)
         DO 800 IC3 = 1, NCOOR
         DO 800 IC2 = 1, NCOOR
         DO 800 IC1 = 1, NCOOR
            HTESTM(IC1,IC3) = HTESTM(IC1,IC3) 
     &                      + TM1TMP(IC1,IC2)*TM2TMP(IC2,IC3)
 800     CONTINUE
C
         DO 900 IC2 = 1, NCOOR
         DO 900 IC1 = 1, IC2
            IF (ABS(HTESTM(IC1,IC2)-HTESTM(IC2,IC1)).GT.DMTHR) 
     &           CALL QUIT('Mass-weighted hessian is not symmetric.')
 900     CONTINUE
C
         CALL HEADER('Mass weighted Hessian in symmetry coordinates',-1)
         NUMTIM = (NCOOR-1)/9 + 1
         DO 1100 ITIM = 1, NUMTIM
            ISTART =     9*(ITIM-1) + 1
            IEND   = MIN(9* ITIM  ,NCOOR)
            DO 1200 IC1 = 1, NCOOR
               WRITE (LUPRI,'(9F12.7)')(HTESTM(IC1,IC2),IC2=ISTART,IEND)
 1200       CONTINUE
            WRITE (LUPRI,'(A)') '                                  '
 1100    CONTINUE
      END IF
C
C     *** Diagonalizing the mass weighted Hessian. ***
C
      KWRK  = 1
      KIWRK = KWRK + NCOOR
      CALL DUNIT(EGNVCT,NCOOR)
      CALL JACO(HESMWT,EGNVCT,NCOOR,NCOOR,NCOOR,WORK(KWRK),WORK(KIWRK))
C
C     *** Mass-weighting the normal coordinates. ***
C
      KDIM = NCOOR**2
      CALL DCOPY(KDIM,EGNVCT,1,TM1TMP,1)
      CALL DZERO(EGNVCT,KDIM)
      DO 1300 IC3 = 1, NCOOR
      DO 1300 IC2 = 1, NCOOR
      DO 1300 IC1 = 1, NCOOR
         EGNVCT(IC1,IC3) = EGNVCT(IC1,IC3) 
     &                   + TM2TMP(IC1,IC2)*TM1TMP(IC2,IC3)
 1300 CONTINUE
C
C     *** Normalizing the normal coordinates. ***
C
      DO 1400 IC2 = 1, NCOOR
         RLENGTH = D0
         DO 1500 IC1 = 1, NCOOR
            RLENGTH = RLENGTH + EGNVCT(IC1,IC2)**2
 1500    CONTINUE

         DRINV = D1/SQRT(RLENGTH)
         DO 1600 IC1 = 1, NCOOR
            EGNVCT(IC1,IC2) = EGNVCT(IC1,IC2)*DRINV
 1600    CONTINUE
 1400 CONTINUE
C
C     *** Removing the redundant normal coordinates, and setting some ***
C     *** common variables according to this                          ***
C
      IC12 = 0
      IFREQ = 0
      NUMZRO = 0
      DO 1700 IC = 1, NCOOR
         IC12 = IC12 + IC
         IF (ABS(HESMWT(IC12)).GT.DMTHR) THEN
            IFREQ = IFREQ + 1
            FREQ(IFREQ) = SQRT(ABS(HESMWT(IC12)))
         ELSE
            NUMZRO = NUMZRO + 1
            DO 1800 IC2 = IC-NUMZRO+1, NCOOR-1
               ICRIRP(IC2,1) = ICRIRP(IC2+1,1)
               ICRIRP(IC2,2) = ICRIRP(IC2+1,2)
               DO 1900 IC1 =  1, NCOOR
                  EGNVCT(IC1,IC2) = EGNVCT(IC1,IC2+1)
 1900          CONTINUE
C
 1800       CONTINUE
C
         END IF
 1700 CONTINUE
      NDCOOR = NCOOR - NUMZRO
C
C     *** Calculating the normal coordinates in cartesian coordinates ***
C     
      CALL DZERO(TM1TMP,KDIM)
      DO 2100 IC3 = 1, NCOOR-NUMZRO
      DO 2100 IC2 = 1, NCOOR
      DO 2100 IC1 = 1, NCOOR
         TM1TMP(IC1,IC3) = TM1TMP(IC1,IC3) 
     &                   + SYMCOR(IC1,IC2)*EGNVCT(IC2,IC3)
 2100 CONTINUE
C
C     *** Using these coordinates in future differentiations. ***
C
      DO 2200 IC2 = 1, NCOOR-NUMZRO
      DO 2200 IC1 = 1, NCOOR
         SYMCOR(IC1,IC2) = TM1TMP(IC1,IC2)
 2200 CONTINUE
C
C     *** Printing the frequencies. ***
C
      WRITE (LUPRI, '(10X,A,I3)') 
     &                  'Number of modes with zero frequency: ', NUMZRO
      CALL HEADER('Vibrational frequencies harmonic approximation: ',-1)
      WRITE (LUPRI,'(20X,A)') ' Mode     cm-1       hartrees  ' 
      DO 2300 IFREQ = 1, NCOOR-NUMZRO
         WRITE (LUPRI,'(20X,I4,F12.2,F12.6)') ICRIRP(IFREQ,1), 
     &                                XTKAYS*FREQ(IFREQ), FREQ(IFREQ)
 2300 CONTINUE
C
C     *** Printing the cartesian components of the normal coordinates ***
C
      NONZRO = NCOOR-NUMZRO
      CALL HEADER('Normal coordinates:  ',0)
      NUMTIM = (NONZRO-1)/9 + 1
      DO ITIM = 1, NUMTIM
         ISTART =     9*(ITIM-1)+1
         IEND   = MIN(9* ITIM     ,NONZRO)
         WRITE (LUPRI,'(I11,8I13)') (ICRIRP(I,1),I=ISTART,IEND)
         DO ICOOR = 1, NCOOR
            WRITE (LUPRI,'(9F13.7)') (SYMCOR(ICOOR,I),I=ISTART,IEND)
         END DO
         WRITE (LUPRI,'(A)') '                                        '
      END DO
C
C     *** Test printing. ***
C
      IF (IPRINT .GE. 20) THEN
         CALL HEADER('Eigen-vectors of the symmetry adapted hessian',-1)
         DO IC1 = 1, NCOOR
            WRITE (LUPRI,'(10X,12F9.4)') (EGNVCT(IC1,IC2),IC2=1,NCOOR)
         END DO
C
         CALL HEADER ('Atomic masses used',-1)
         WRITE (LUPRI,'(3X,A,10I9)') 'Atom number:', (I,I=1,NATOMS)
         WRITE (LUPRI,'(17X, 10F9.4 )') (AMASS(I), I = 1, NATOMS)
C
         CALL HEADER('Diagonal elements of (sqrt(mass))^-1 matrix',-1)
         WRITE (LUPRI,'(24F9.4)') (DKIN(IC), IC=1,NCOOR)
C
         CALL HEADER('Transformed sqrt(mass)^-1 matrix',-1)
         DO IC1 = 1, NCOOR
            WRITE (LUPRI,'(24F9.7)') (TM2TMP(IC1,IC2),IC2= 1, NCOOR)
         END DO
C
         CALL HEADER('Diagonalized mass weighted hessian.',-1)
         DO IC1 = 1, NCOOR
            WRITE (LUPRI,'(24F9.7)') 
     &            (HESMWT(IC12),IC12= IC1*(IC1-1)/2+1, IC1*(IC1+1)/2)
         END DO
      END IF
      RETURN
      END
C
C     
C     /*Deck prderv*/
      SUBROUTINE PRDERV(TDER,FDER,TSTGDR,TSTSDR,SYMCOR,WORK,ICRIRP,
     &                  LWORK,NPRRDR,NDIMT,NDIMF,LTXT,IPRINT,TEXT)
C     ***************************************************************
C     **** This routine prints out the derivatives of the energy ****
C     **** to NPRRDR order. These are done in 'TEXT' coordinates ****
C     ***************************************************************
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
C
#include <nuclei.h>
#include <trkoor.h>
#include <taymol.h>
#include <cbiwlk.h>
#include <numder.h>
      CHARACTER*80 TEXT
      DIMENSION TDER(NDIMT), FDER(NDIMF), TSTGDR(NCOOR), 
     &          TSTSDR(NCOOR,NCOOR), SYMCOR(NCOOR,NCOOR), WORK(LWORK)
      DIMENSION ICRIRP(NCOOR,2)
C
C     *** Print gradient ***
C
      IF (NPRRDR.GT.0) THEN
C     
         CALL HEADER('Numerical gradient in ' // TEXT(1:LTXT) // 
     &               ' coordinates',-1)
C
         KCGRAD = 1
         KWRK   = KCGRAD + MXCOOR
         LWRK   = LWORK - KWRK
         IF ((LWRK - 2*MXCOOR*MXCOOR) .LT. 0)
     &        CALL STOPIT('Too little memory in PRDERV')
         CALL GSPGRD(SYMCOR,WORK(KCGRAD),WORK(KWRK),LWRK,ICRIRP,LTXT,
     &               IPRINT,TEXT)
C
C        *** If comparing with the analytical gradient. ***
C            
         IF (SDRTST) THEN
            DO 200 IC1 = 1, NCOOR
               TSTGDR(IC1) = WORK(KCGRAD-1+IC1)
 200        CONTINUE
         END IF
      END IF
C
C     *** Print hessian ***
C
      IF (NPRRDR.GT.1) THEN
C
         CALL HEADER('Numerical hessian in ' // TEXT(1:LTXT) //
     &               ' coordinates',-1)  
C
         KSMCIN = 1
         KCHES1 = KSMCIN + NCOOR**2
         KCHES2 = KCHES1 + NCOOR**2
         KLAST  = KCHES2 + NCOOR**2
         LWRK   = LWORK  - KLAST
         CALL GSPHES(SYMCOR,WORK(KSMCIN),WORK(KCHES1),WORK(KCHES2),
     &               WORK(KLAST),ICRIRP,NDCOOR,LWRK,LTXT,IPRINT,TEXT)
C
C        *** If comparing with the analytical hessian. ***
C
         IF (SDRTST) THEN
            IC12 = 0
            DO 300 IC2 = 1, NCOOR
            DO 300 IC1 = 1, NCOOR
               IC12 = IC12 + 1
               TSTSDR(IC1,IC2) = WORK(KCHES2-1+IC12)
 300        CONTINUE
         END IF
      END IF
C
C     *** Print third derivative of energy ***
C      
      IF (NPRRDR.GT.2) THEN
         KTDER  = 1
         KSYMTD = KTDER  + NCOOR**3
         KSCTRA = KSYMTD + NCOOR**3
         KCSTRA = KSCTRA + NCOOR**2
         KLAST  = KCSTRA + NCOOR**2
         CALL HEADER('Numerical third derivative of energy in ' // 
     &                TEXT(1:LTXT) // 'coordinates',-1)
         CALL PRITDR(TDER,SYMCOR,WORK(KTDER),WORK(KSYMTD),WORK(KSCTRA),
     &               WORK(KCSTRA),NDIMT,LTXT,TEXT)
      END IF
C
C     *** Print fourth derivative of energy ***
C
      IF (NPRRDR.GT.3) THEN
         KFDER  = 1
         KSYMFD = KFDER  + NCOOR**4
         KSCTRA = KSYMFD + NCOOR**4
         KCSTRA = KSCTRA + NCOOR**2
         KLAST  = KCSTRA + NCOOR**2
         CALL HEADER('Numerical fourth derivative of energy in ' //
     &                TEXT(1:LTXT) // 'coordinates',-1)
         CALL PRIFDR(FDER,WORK(KFDER),WORK(KSYMFD),WORK(KSCTRA),
     &               WORK(KCSTRA),NDIMF,LTXT,TEXT)
      END IF
C
      IF ((IPRINT .GT. 20).AND.(NPRRDR.EQ.2)) THEN
         CALL HEADER ('Copy of Hessian for test', -1)
         DO IC1 = 1, NCOOR
            WRITE (LUPRI,'(12F12.8)') (TSTSDR(IC1,IC2),IC2=1,NCOOR)
         END DO
      END IF
      RETURN
      END
C     
C
C
C     /*Deck gspgrd*/
      SUBROUTINE GSPGRD(SYMCOR,CGRAD,WORK,LWORK,ICRIRP,LTXT,IPRINT,TEXT)
#include <implicit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <priunit.h>
C
#include <symmet.h>
#include <nuclei.h>
#include <trkoor.h>
#include <taymol.h>
#include <cbiwlk.h>
#include <numder.h>
      CHARACTER TEXT*80
      DIMENSION SYMCOR(NCOOR,NCOOR), CGRAD(NCOOR), ICRIRP(NCOOR,2),
     &          WORK(LWORK)
C
      CALL HEADER('Gradient in ' // TEXT(1:LTXT) // ' coordinates',-1)
C
      DO 100 ICOOR = 1, NDCOOR
         IF (ICRIRP(ICOOR,1).EQ.1) THEN
            WRITE (LUPRI,'(A,F16.8)') '                               ',
     &                                GRDMOL(ICOOR)
         END IF
 100  CONTINUE
C     
      CALL DZERO(CGRAD,NCOOR)
      DO 200 ICOOR2 = 1, NDCOOR
      DO 200 ICOOR1 = 1, NCOOR
         CGRAD(ICOOR1) = CGRAD(ICOOR1) 
     &                 + SYMCOR(ICOOR1,ICOOR2)*GRDMOL(ICOOR2)
 200  CONTINUE
C
      CALL HEADER('Gradient in cartesian coordinates',-1)
C
      IOFF = 0
      DO 300 ICENT = 1, NUCDEP
         WRITE (LUPRI,'(1X,A6,F17.10,2F24.10)') NAMDEP(ICENT), 
     &                                         (CGRAD(IOFF+J), J=1,3)
         IOFF = IOFF + 3
 300  CONTINUE
C
C     Transform to symmetry basis used in Dalton in case of 
C     geometry optimization
C
      IF (MAXREP .GT. 0) THEN
         KCSTRA = 1
         KSCTRA = KCSTRA + MXCOOR*MXCOOR
         CALL TRACOR(WORK(KCSTRA),WORK(KSCTRA),1,MXCOOR,IPRINT)
         CALL TRACTS(CGRAD,3*NUCDEP,WORK(KCSTRA))
         CALL DCOPY(3*NUCDEP,CGRAD,1,GRDMOL,1)
C
         CALL HEADER('Gradient in Dalton symmetry coordinates',-1)
C
         DO 202 I = 1, NCRREP(0,1)
            WRITE (LUPRI,'(25X,A6,F17.10)') NAMEX(IPTCOR(I,1)),GRDMOL(I)
 202     CONTINUE
      END IF
C
C     *** If testing ***
C
      IF (SDRTST) THEN
         DO IC1 = 1, NCOOR
            GRDMOL(IC1) = CGRAD(IC1)
         END DO
      END IF
C
C     *** Print ***
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Symcor matrix in GSPGRD',-1)
         DO 400 I = 1, NCOOR
            WRITE (LUPRI,'(24F12.7)') (SYMCOR(I,J),J=1,NCOOR)
 400     CONTINUE
      END IF
C
      RETURN
      END 
C
C
C     /*Deck gsphes*/
      SUBROUTINE GSPHES(SYMCOR,SMCINV,CHESS1,CHESS2,WORK,ICRIRP,NDCOOR,
     &                  LWORK,LTXT,IPRINT,TEXT)
#include <implicit.h>
#include <mxcent.h>
#include <priunit.h>
      PARAMETER (KCOL=6)
#include <nuclei.h>
#include <trkoor.h>
#include <taymol.h>
#include <cbiwlk.h>
      INTEGER BEGIN, LAST
      CHARACTER TEXT*80
      DIMENSION SYMCOR(NCOOR,NCOOR), SMCINV(NCOOR,NCOOR), 
     &          CHESS1(NCOOR,NCOOR), CHESS2(NCOOR,NCOOR),
     &          WORK(LWORK), ICRIRP(NCOOR,    2)
C
C     *** Print hessian in symmetric coordinates ***
C
      CALL HEADER('Hessian in ' // TEXT(1:LTXT) // ' coordinates',-1)
C
      BEGIN = 1
      LAST  = MIN(NDCOOR,KCOL)
      KCOOR = NDCOOR
      NCOL  = INT(DBLE(NDCOOR)/DBLE(KCOL)) + 1
C
      DO 100 ICOL = 1, NCOL
         WRITE (LUPRI,1000) (ICRIRP(I,1),I = BEGIN,LAST)
C
         DO 200 ICOOR = BEGIN, NDCOOR
            WRITE (LUPRI,2000) ICRIRP(ICOOR,1), 
     &                        (HESMOL(ICOOR,I),I=BEGIN,MIN(LAST,ICOOR))
 200     CONTINUE
         WRITE (LUPRI,'(A)') '                                       '
         BEGIN = BEGIN + LAST
         LAST  = MIN(NDCOOR,KCOL+LAST)
 100  CONTINUE
C
      DO 300 J = 1, NDCOOR
      DO 300 I = 1, J
         HESMOL(I,J) = HESMOL(J,I)
 300  CONTINUE
C
C     *** Transpose of coordinate transformation matrix ***
C
      IF (TEXT(1:6) .NE. 'normal') THEN
         DO 350 J = 1, NDCOOR
         DO 350 I = 1, NCOOR
            SMCINV(J,I) = SYMCOR(I,J)
 350     CONTINUE
C
C     *** Transform to cartesian hessian ***
C
         KDIM = NCOOR**2
         CALL DZERO(CHESS1,KDIM)
         DO 400 K = 1, NCOOR
         DO 400 J = 1, NDCOOR
         DO 400 I = 1, NDCOOR
            CHESS1(I,K) = CHESS1(I,K) + HESMOL(I,J)*SMCINV(J,K)
 400     CONTINUE
C
         KDIM = NCOOR**2
         CALL DZERO(CHESS2,KDIM)
         DO 500 K = 1, NCOOR
         DO 500 J = 1, NDCOOR
         DO 500 I = 1, NCOOR
            CHESS2(I,K) = CHESS2(I,K) + SYMCOR(I,J)*CHESS1(J,K)
 500     CONTINUE
C
C        *** Print cartesian hessian ***
C
         CALL HEADER('Cartesian Hessian in GSPHES',-1)
C      
         BEGIN = 1
         LAST  = MIN(NCOOR,KCOL)
         KCOOR = NCOOR
         NCOL  = INT(DBLE(NCOOR)/DBLE(KCOL)) + 1
C
         DO 600 ICOL = 1, NCOL
            WRITE (LUPRI,3000) (NAMEX(I),I = BEGIN,LAST)
            DO 700 ICOOR = BEGIN, NCOOR
               WRITE (LUPRI,4000) NAMEX(ICOOR), 
     &                        (CHESS2(ICOOR,I),I=BEGIN,MIN(LAST,ICOOR)) 
 700        CONTINUE
            WRITE (LUPRI,'(A)') '                                      '
            BEGIN = BEGIN + LAST
            LAST  = MIN(NCOOR,KCOL+LAST)
 600     CONTINUE
C
C        *** Print ***
C
         IF (IPRINT .GT. 20) THEN
C     
            KDIM = NCOOR**2
            CALL DZERO(CHESS1,KDIM)
C     
            CALL HEADER('Symcor matrix',-1)
            DO 1100 I = 1, NCOOR
               WRITE (LUPRI,'(24F9.6)') (SYMCOR(I,J),J=1,NDCOOR)
 1100       CONTINUE
            WRITE (LUPRI,'(A)') '                                  '
C     
            CALL HEADER('Inverse of symcor matrix',-1)
            DO 1200 I = 1, NDCOOR
               WRITE (LUPRI,'(24F9.6)') (SMCINV(I,J),J=1,NCOOR)
 1200       CONTINUE
            WRITE (LUPRI,'(A)') '                                  '
C
            DO 1300 K = 1, NCOOR
            DO 1300 J = 1, NDCOOR
            DO 1300 I = 1, NCOOR
               CHESS1(I,K) = CHESS1(I,K) + SYMCOR(I,J)*SMCINV(J,K)
 1300       CONTINUE
C
            CALL HEADER('Should be unit matrix',-1)
            DO 1400 J = 1, NCOOR
               WRITE (LUPRI,'(24F9.6)') (CHESS1(I,J),I=1,NCOOR)
 1400       CONTINUE
         END IF
      END IF
C         
 1000 FORMAT (8X,6(3X,I4,5X),(3X,I4,5X))
 2000 FORMAT (1X,I4,2X,6F12.6)
 3000 FORMAT (8X,6(3X,A6,3X),(3X,A6,5X))
 4000 FORMAT (1X,A6,2X,6F12.6)
      RETURN
      END
C
C
C     /* Deck tsths1*/
      SUBROUTINE TSTHS1(SYMCOR,HESMWT,EGNVCT,TM1TMP,TM2TMP,DKIN,WORK,
     &                  ICRIRP,LWORK)
************************************************************
*** Tests the molecular hessian by makeing the cartesian ***
*** hessian, diagonalizing the mass-weighted hessian and ***
*** prints the eigenvalues, and harmonic frequancies     ***
************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <codata.h>
      PARAMETER (DMTHR = 1.0D-9, D1 = 1.0D0)
#include <taymol.h>
#include <trkoor.h>
#include <nuclei.h>
#include <cbinum.h>
      DIMENSION SYMCOR(NCOOR,NCOOR), EGNVCT(NCOOR,NCOOR),
     &          HESMWT(NCOOR*(NCOOR+1)/2), DKIN(NCOOR,NCOOR), 
     &          TM1TMP(NCOOR,NCOOR),TM2TMP(NCOOR,NCOOR),
     &          WORK(LWORK)
C
      KDIM = NCOOR**2
      CALL DZERO(TM1TMP,KDIM)
      DO 100 IC3 = 1, NCOOR
      DO 100 IC2 = 1, NCOOR
      DO 100 IC1 = 1, NCOOR
         TM1TMP(IC1,IC3) = TM1TMP(IC1,IC3) 
     &                   + HESMOL(IC1,IC2)*SYMCOR(IC3,IC2)
 100  CONTINUE
C      
      CALL DZERO(TM2TMP,KDIM)
      DO 200 IC3 = 1, NCOOR
      DO 200 IC2 = 1, NCOOR
      DO 200 IC1 = 1, NCOOR
         TM2TMP(IC1,IC3) = TM2TMP(IC1,IC3)
     &                   + SYMCOR(IC1,IC2)*TM1TMP(IC2,IC3)
 200  CONTINUE
C
      CALL HEADER('Cartesian hessian',-1)
      DO 300 IC1 = 1, NCOOR
         WRITE (LUPRI,'(12F8.5)') (HESMOL(IC1,IC2),IC2=1,NCOOR)
 300  CONTINUE
C
      CALL DZERO(TM1TMP,KDIM)
      DO 400 IC3 = 1, NCOOR
      DO 400 IC2 = 1, NCOOR
      DO 400 IC1 = 1, NCOOR
         TM1TMP(IC1,IC3) = TM1TMP(IC1,IC3) 
     &                   + TM2TMP(IC1,IC2)*DKIN(IC2,IC3)
 400  CONTINUE
C
      CALL DZERO(TM2TMP,KDIM)
      DO 500 IC3 = 1, NCOOR
      DO 500 IC2 = 1, NCOOR
      DO 500 IC1 = 1, NCOOR
         TM2TMP(IC1,IC3) = TM2TMP(IC1,IC3)
     &                   + DKIN(IC1,IC2)*TM1TMP(IC2,IC3)
 500  CONTINUE
C
      IC12 = 0
      DO 600 IC2 = 1, NCOOR
      DO 600 IC1 = 1, IC2
         IC12 = IC12 + 1
         HESMWT(IC12) = TM2TMP(IC1,IC2)
 600  CONTINUE
C
      CALL HEADER('Mass-weighted hessian',-1)
      DO 700 IC1 = 1, NCOOR
         WRITE (LUPRI,'(12F8.5)') (TM2TMP(IC1,IC2),IC2=1,IC1)
 700  CONTINUE
C
      KWRK  = 1
      KIWRK = KWRK + NCOOR
      CALL DUNIT(EGNVCT,NCOOR)
      CALL JACO(HESMWT,EGNVCT,NCOOR,NCOOR,NCOOR,WORK(KWRK),WORK(KIWRK))
C
      CALL HEADER('Diagonalized hessian',-1)
      DO 800 IC1 = 1, NCOOR
         ISTART = (IC1*(IC1-1))/2 + 1
         IEND   = (IC1*(IC1+1))/2
         WRITE (LUPRI,'(12F8.5)') (HESMWT(IC12),IC22=ISTART,IEND)
 800  CONTINUE
C
      RETURN
      END
C
C
C     /*Deck sdertt*/
      SUBROUTINE SDERTT(TSTSDR,TSTGDR,SYMCOR,TMPGRD,TMPHES,BKPGRD,
     &                  BKPHES,WORK,LWORK,IPRINT)
C     *******************************************************************
C     *** This routine tests the numerical derivatives with available ***
C     ***                   analytical derivatives.                   ***
C     *** NOTE: TMPGRD has dimension MXCOOR due to old code.          ***
C     *******************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
C
#include <numder.h>
#include <trkoor.h>
#include <taymol.h>
#include <symmet.h>
#include <abainf.h>
#include <exeinf.h>
#include <gnrinf.h>
#include <past.h>
#include <inftap.h>
      DIMENSION TSTSDR(NCOOR,NCOOR), TSTGDR(NCOOR ),
     &          TMPHES(NCOOR,NCOOR), TMPGRD(MXCOOR),
     &          BKPHES(NCOOR,NCOOR), BKPGRD(NCOOR ),
     &          SYMCOR(NCOOR,NCOOR), WORK  (LWORK )
C
      DO 100 IC2 = 1, NCOOR
         BKPGRD(IC2) = GRDMOL(IC2)
         DO 200 IC1 = 1, NCOOR
            BKPHES(IC1,IC2) = HESMOL(IC1,IC2)
 200     CONTINUE
 100  CONTINUE
C
      DO 300 IC2 = 1, NCOOR
      DO 300 IC1 = 1, NCOOR
         TSTGDR(IC1) = TSTGDR(IC1) + SYMCOR(IC1,IC2)*GRDMOL(IC2)
 300  CONTINUE
C
      CALL DGEMM('N','N',NCOOR,NCOOR,NCOOR,1.D0,
     &           SYMCOR,NCOOR,
     &           HESMOL,MXCOOR,0.D0,
     &           TMPHES,NCOOR)
C
      CALL DGEMM('N','T',NCOOR,NCOOR,NCOOR,1.D0,
     &           TMPHES,NCOOR,
     &           SYMCOR,NCOOR,0.D0,
     &           TSTSDR,NCOOR)
C
      KDIM = MXCOOR**2
      CALL DZERO(GRDMOL,MXCOOR)
      CALL DZERO(HESMOL,KDIM  )
C
      MOLGRD = .TRUE.
      MOLHES = .TRUE.
      PASEXC = .FALSE.
      RNABAC = .TRUE.
      WRINDX = .TRUE. 
      FTRONV = .TRUE.
      DOWALK = .FALSE.
      LUSUPM = -1
      WORK(IWOFF+1) = WRKDLM
      CALL ABAINP('**PROPE',WORK(IWOFF+2),LWORK)
      CALL EXEABA(WORK,LWORK,IWOFF,WRKDLM)
C     
      
C     *** Transforming the analytical derivatives into ***
C     ***               cartesian basis.               ***
C
      IF (MAXREP.GT.0) THEN
         KSCTR = 1
         KCSTR = KSCTR + NCOOR**2
         KLAST  = KCSTR + NCOOR**2
         IF (KLAST.GT.LWORK) CALL QUIT('Memory exceeded inside SDERTT')
         CALL TRAGRD(GRDMOL,TMPGRD,WORK(KCSTR),WORK(KSCTR),NCRREP(0,1),
     &               NCOOR)
         CALL TRAHES(HESMOL,MXCOOR,TMPHES,WORK(KCSTR),WORK(KSCTR),NCOOR,
     &               NCOOR,1)
      ELSE
         DO 400 IC2 = 1, NCOOR
            TMPGRD(IC2) = GRDMOL(IC2)
            DO 500 IC1 = 1, NCOOR
               TMPHES(IC1,IC2) = HESMOL(IC1,IC2)
 500        CONTINUE
 400     CONTINUE
      END IF
C
      RMAXGD = 0.D0
      DO 600 J = 1, NCOOR
         RGRDJ = ABS(TMPGRD(J)-TSTGDR(J))
         IF (RGRDJ .GT. RMAXGD) THEN
            RMAXGD = RGRDJ
            NMG    = J
         END IF
 600  CONTINUE
      CALL HEADER('Comparison of numerical and analytical gradients',-1)
      WRITE (LUPRI,'(//A,E13.7,A,I4)')  'Largest difference ', RMAXGD,
     &                                  ' for element:', NMG
      WRITE (LUPRI,'(A,2E13.7//)') 'The values of these elements are: ',
     &                             TMPGRD(NMG),TSTGDR(NMG)

      RLRGST = 0.D0
      DO 700 J = 1, NCOOR
      DO 700 I = 1, J
         RINTMD = (TMPHES(I,J)-TSTSDR(I,J))**2
         IF ( RINTMD .GT. RLRGST) THEN
            RLRGST = RINTMD
            NMI    = I
            NMJ    = J
            HVALC  = TMPHES(I,J)
            HVALN  = TSTSDR(I,J)
         END IF
 700  CONTINUE
C     
      CALL HEADER('Comparison of numerical and analytical Hessians',-1)
      WRITE (LUPRI,'(//A,E13.7,A,2I4)') 'Largest difference ', 
     &                 SQRT(RLRGST), ' for elements:', NMI, NMJ
      WRITE (LUPRI,'(A,2E13.7//)') 'The values of these elements are: ',
     &                              HVALC, HVALN
C
      DO 800 IC2 = 1, NCOOR
         GRDMOL(IC2) = BKPGRD(IC2)
         DO 900 IC1 = 1, NCOOR
            HESMOL(IC1,IC2) = BKPHES(IC1,IC2)
 900     CONTINUE
 800  CONTINUE
C
      RETURN
      END
C
C     
C     /* Deck drnrmc*/
      SUBROUTINE DRNRMC(SYMCOR,ICRIRP,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
#include <fcsym.h>
      LOGICAL FOUND
      DIMENSION SYMCOR(NCOOR,NCOOR), ICRIRP(NCOOR,2)
C
      DO 100 IC1 = 1, NMREDU
         FOUND = .FALSE.
         IIREP = KDRYRN(IC1)
C         
         DO 200 IC2 = 1, NDCOOR
            IF ((ICRIRP(IC2,1).EQ.IIREP).AND.(.NOT.FOUND)) THEN
               FOUND = .TRUE.
               DO 300 II  = 1,     2
               DO 300 IC3 = IC2+1, NDCOOR
                  ICRIRP(IC3-1,II) = ICRIRP(IC3,II)
 300           CONTINUE
               NDCOOR = NDCOOR - 1
C
               IF (IIREP.GT.N1DIME) THEN
                  FOUND = .FALSE.
                  DO 400 IC3 = IC2-1, NDCOOR
                     IF ((ICRIRP(IC3,1).EQ.IIREP).AND.
     &                   (ICRIRP(IC3,2).EQ.    1).AND.(.NOT.FOUND)) THEN
                        FOUND = .TRUE.
                        DO 500 II  =     1,      2
                        DO 500 IC4 = IC3+1, NDCOOR
                           ICRIRP(IC4-1,II) = ICRIRP(IC4,II)
 500                    CONTINUE
                     END IF
 400              CONTINUE
                  NDCOOR = NDCOOR - 1
               END IF
            END IF
 200     CONTINUE
 100  CONTINUE
C
      IF (IPRINT.GT.20) THEN
         WRITE (LUPRI,'(5X,A)') 'Removed translational and rotational'//
     &         'redundencies.'
         WRITE (LUPRI,'(A)') '                                    '
         WRITE (LUPRI,'(5X,A)') 'Symmetry of coordinates left:' 
         WRITE (LUPRI,'(5X,24I4)') (ICRIRP(II,1),II=1,NDCOOR)
         WRITE (LUPRI,'(5X,24I4)') (ICRIRP(II,2),II=1,NDCOOR)
      END IF
      RETURN
      END
C
C
C     /* Deck wricor*/
      SUBROUTINE WRICOR(SYMCOR,ICRIRP,LURSTR,IPRINT)
C     *************************************************
C     *** This is a routine that writes out normal  ***
C     *** coordinates to file, in case of a restart.***
C     *************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
      DIMENSION SYMCOR(NCOOR,NCOOR)
      DIMENSION ICRIRP(NCOOR,2)
C
      WRITE (LURSTR,'(I8)') NDCOOR
C
      DO 100 IC2 = 1, NDCOOR
         WRITE(LURSTR,'(2I5)') (ICRIRP(IC2,I),I=1,2)
         DO 200 IC1 = 1, NCOOR   
            WRITE(LURSTR,'(F24.16)') SYMCOR(IC1,IC2)
 200     CONTINUE
 100  CONTINUE
      CALL FLSHFO(LURSTR)
C
      RETURN
      END
C
C
C     /*Deck rerstr*/
      SUBROUTINE RERSTR(ENVAL,SYMCOR,ICRIRP,NDIME,KEND,IDIMAX,IDIMIN,
     &                  LURSTR,IPRINT,RSTDON)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
C
#include <trkoor.h>
#include <numder.h>
      LOGICAL RSTDON
      DIMENSION SYMCOR(NCOOR,NCOOR), ENVAL(NDIME)
      DIMENSION ICRIRP(NCOOR,2)
C
      READ(LURSTR,'(I2)') II
C
C     *** The program ended the first time around. ***
C
      IF ((II.EQ.0) .AND. (II.EQ.KEND)) THEN
         RSTDON = .TRUE.
C
C        *** Reading the function values ***
C
 100     CONTINUE
         READ(LURSTR,FMT='(I8,F24.16)',IOSTAT=IOS) IDIME, ENERGY
         IF (IOS.GE.0) THEN
            IDIMAX = MAX(IDIMAX,IDIME)
            IDIMIN = MIN(IDIMIN,IDIME)
            ENVAL(IDIME) = ENERGY
            GOTO 100
         END IF
      ELSE IF ((II.EQ.1) .AND. (II.EQ.KEND)) THEN
         RSTDON = .TRUE.
C
C        *** Number of normal coordinates ***
C
         READ(LURSTR,FMT='(I8)') NDCOOR
C
C        *** The normal coordinates. ***
C
         DO 200 IC2 = 1, NDCOOR
            READ(LURSTR,'(2I5)') (ICRIRP(IC2,I),I=1,2)
            DO 300 IC1 = 1, NCOOR
               READ(LURSTR,'(F24.16)') SYMCOR(IC1,IC2)
 300        CONTINUE
 200     CONTINUE
C
C        *** Read the function values. ***
C
 400     CONTINUE
         READ(LURSTR,FMT='(I8,F24.16)',IOSTAT=IOS) IDIME, ENERGY
         IF (IOS.GE.0) THEN
            IDIMAX = MAX(IDIMAX,IDIME)
            IDIMIN = MIN(IDIMIN,IDIME)
            ENVAL(IDIME) = ENERGY
            GOTO 400
         END IF
      END IF
C
      RETURN
      END
