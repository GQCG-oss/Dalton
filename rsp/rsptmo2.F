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
      SUBROUTINE TM1IND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT ONE-INDEX
C LINEAR RESPONSE EQUATIONS THAT ARE USED IN A
C A CALCULATION OF THIRD ORDER TRANSITION MOMENTS
C
#include <implicit.h>
C
#include <priunit.h>
#include <infsmo.h>
#include <indcr.h>
#include <inforb.h>
#include <infrsp.h>
#include <infspi.h>
#include <inftmo.h>
#include <infcr.h>
C
      DO 200 IDSYM = 1,NSYM
      DO 300 ICSYM = 1,NSYM
      DO 400 IBSYM = 1,NSYM
         IASYM = MULD2H(IDSYM,MULD2H(ICSYM,IBSYM))
         IF ( (NTMCNV(IDSYM).GT.0) .AND. (NCTMOP(ICSYM).GT.0) .AND.
     *       (NBTMOP(IBSYM).GT.0) .AND. (NATMOP(IASYM).GT.0) ) THEN
            DO 450 IDFR = 1,NTMCNV(IDSYM)
               INUM = INCRLR('EXCITLAB',EXCIT2(IDSYM,IDFR),IDSYM)
 450        CONTINUE
            DO 500 ICOP = 1,NCTMOP(ICSYM)
            DO 600 ICFR = 1,NCTMFR
               INUM = INCRLR(CTMLB(ICSYM,ICOP),-CTMFR(ICFR),ICSYM)
 600        CONTINUE
 500        CONTINUE
            DO 700 IBOP = 1,NBTMOP(IBSYM)
            DO 800 IBFR = 1,NBTMFR
               INUM = INCRLR(BTMLB(IBSYM,IBOP),-BTMFR(IBFR),IBSYM)
 800        CONTINUE
 700        CONTINUE
            DO 900  IDFR = 1,NTMCNV(IDSYM)
            DO 1000 ICFR = 1,NCTMFR
            DO 1100 IBFR = 1,NBTMFR
C
C  If harmonic generation is specified we only accept an equal
C  number of b,c frequencies
C
               IF (CTMOHG) THEN
                  IF (NBTMFR.NE.NCTMFR) THEN
                      WRITE (LUPRI,'(2(/A))')
     * '*** ERROR -- If harmonic generation is specified we only',
     * '             accept an equal number of b,c frequencies.'
                      CALL QUIT('wrong frequencies')
                   END IF
                  IF (IBFR.NE.ICFR) GO TO 1100
               END IF
               ATMFR = EXCIT2(IDSYM,IDFR) - CTMFR(ICFR) - BTMFR(IBFR)
               DO 1200 IAOP = 1,NATMOP(IASYM)
                  INUM = INCRLR(ATMLB(IASYM,IAOP),ATMFR,IASYM)
 1200          CONTINUE
 1100       CONTINUE
 1000       CONTINUE
 900        CONTINUE
         END IF
 400  CONTINUE
 300  CONTINUE
 200  CONTINUE
      RETURN
      END
      SUBROUTINE TM2IND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT
C TWO-INDEX LINEAR RESPONSE EQUATIONS THAT NEED TO BE SOLVED
C IN A SECOND HYPERPOLARIZABILITY CALCULATION
C
#include <implicit.h>
C
#include <infcr.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
#include <infspi.h>
#include <inftmo.h>
#include <indcr.h>
C
      CHARACTER*8 DTMLB
C
C Put label EXCITLABin list for two-index vectors
C for vectors of the type N^BX and N^CX
C
      DATA DTMLB/'EXCITLAB'/
C
      DO 300 IDSYM = 1,NSYM
      DO 200 ICSYM = 1,NSYM
      DO 100 IBSYM = 1,NSYM
         IASYM = MULD2H(IDSYM,MULD2H(ICSYM,IBSYM))
         IF ( (NTMCNV(IDSYM).GT.0) .AND. (NCTMOP(ICSYM).GT.0)
     *     .AND. (NBTMOP(IBSYM).GT.0) .AND. (NATMOP(IASYM).GT.0)) THEN
            DO 210 ICOP = 1,NCTMOP(ICSYM)
            DO 110 IBOP = 1,NBTMOP(IBSYM)
               DO 320 IDFR = 1,NTMCNV(IDSYM)
               DO 220 ICFR = 1,NCTMFR
               DO 120 IBFR = 1,NBTMFR
                  DTMFR = EXCIT2(IDSYM,IDFR)
C
C  If harmonic generation is specified we only need two index
C  response vectors with equal B and C frequencies.
C
                  IF (CTMOHG) THEN
                     IF (IBFR.NE.ICFR) GO TO 120
                  END IF
C
                  INUM = INCR2(BTMLB(IBSYM,IBOP),CTMLB(ICSYM,ICOP),
     *                          -BTMFR(IBFR),-CTMFR(ICFR),IBSYM,ICSYM)
                  INUM = INCR2(BTMLB(IBSYM,IBOP),DTMLB,
     *                          -BTMFR(IBFR),DTMFR,IBSYM,IDSYM)
                  INUM = INCR2(CTMLB(ICSYM,ICOP),DTMLB,
     *                          -CTMFR(ICFR),DTMFR,ICSYM,IDSYM)
120            CONTINUE
220            CONTINUE
320            CONTINUE
110         CONTINUE
210         CONTINUE
         END IF
100   CONTINUE
200   CONTINUE
300   CONTINUE
      RETURN
      END

