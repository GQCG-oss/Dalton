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
      SUBROUTINE TP1IND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT ONE-INDEX
C LINEAR RESPONSE EQUATIONS THAT ARE USED IN A
C A CALCULATION OF TWO-PHOTON EXCITATION
C
#include <implicit.h>
C
#include <rspprp.h>
#include <infsmo.h>
#include <indcr.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
#include <infspi.h>
#include <inftpa.h>
#include <infcr.h>
C
C If excited state alpha (TPALP) is specified we only compute certain
C components.
C
      IF (TPALP) THEN
C
      DO 600 ICSYM = 1,NSYM
      DO 620 IBSYM = 1,NSYM
         IDSYM = ICSYM
         IASYM = IBSYM
         IF ( (NTPCN2(IDSYM).GT.0) .AND. (NTPCN1(ICSYM).GT.0) .AND.
     *    (NBTPOP(IBSYM).GT.0) .AND. (NATPOP(IASYM).GT.0) ) THEN
            DO ICFR = 1,NTPCN1(ICSYM)
               INUM = INCRLR('EXCITLAB',EXCIT2(ICSYM,ICFR),ICSYM)
            ENDDO
            DO IBOP = 1,NBTPOP(IBSYM)
            DO IBFR = 1,NBTPFR
               INUM = INCRLR(BTPLB(IBSYM,IBOP),-BTPFR(IBFR),IBSYM)
            ENDDO
            ENDDO
            DO ICFR = 1,NTPCN1(ICSYM)
            DO IBFR = 1,NBTPFR
               IDFR = ICFR
               ATPFR = EXCIT2(IDSYM,IDFR)-EXCIT2(ICSYM,ICFR)-BTPFR(IBFR)
               DO IAOP = 1,NATPOP(IASYM)
                  INUM = INCRLR(ATPLB(IASYM,IAOP),ATPFR,IASYM)
               ENDDO
            ENDDO
            ENDDO
         END IF
 620  CONTINUE
 600  CONTINUE
C
      ELSE
C
      DO 200 IDSYM = 1,NSYM
      DO 300 ICSYM = 1,NSYM
      DO 400 IBSYM = 1,NSYM
         IASYM = MULD2H(IDSYM,MULD2H(ICSYM,IBSYM))
         IF ( (NTPCN2(IDSYM).GT.0) .AND. (NTPCN1(ICSYM).GT.0) .AND.
     *       (NBTPOP(IBSYM).GT.0) .AND. (NATPOP(IASYM).GT.0) ) THEN
            DO 450 IDFR = 1,NTPCN2(IDSYM)
               INUM = INCRLR('EXCITLAB',EXCIT2(IDSYM,IDFR),IDSYM)
 450        CONTINUE
            DO 460 ICFR = 1,NTPCN1(ICSYM)
               INUM = INCRLR('EXCITLAB',EXCIT2(ICSYM,ICFR),ICSYM)
 460        CONTINUE
            DO 700 IBOP = 1,NBTPOP(IBSYM)
            DO 800 IBFR = 1,NBTPFR
               INUM = INCRLR(BTPLB(IBSYM,IBOP),-BTPFR(IBFR),IBSYM)
 800        CONTINUE
 700        CONTINUE
            DO 900  IDFR = 1,NTPCN2(IDSYM)
            DO 1000 ICFR = 1,NTPCN1(ICSYM)
            DO 1100 IBFR = 1,NBTPFR
               ATPFR = EXCIT2(IDSYM,IDFR)-EXCIT2(ICSYM,ICFR)-BTPFR(IBFR)
               DO 1200 IAOP = 1,NATPOP(IASYM)
                  INUM = INCRLR(ATPLB(IASYM,IAOP),ATPFR,IASYM)
 1200          CONTINUE
 1100       CONTINUE
 1000       CONTINUE
 900        CONTINUE
         END IF
 400  CONTINUE
 300  CONTINUE
 200  CONTINUE
C
      END IF
C
      RETURN
      END
      SUBROUTINE TP2IND
C
C CALCULATE A POINTER TO THE NUMBER OF DIFFERENT
C TWO-INDEX LINEAR RESPONSE EQUATIONS THAT NEED TO BE SOLVED
C IN A CUBIC RESPONSE TWO-PHOTON CALCULATION
C
#include <implicit.h>
C
#include <rspprp.h>
#include <infcr.h>
#include <inforb.h>
#include <infrsp.h>
#include <infpri.h>
#include <infspi.h>
#include <inftpa.h>
#include <indcr.h>
C
      CHARACTER*8 CTPLB, DTPLB
C
C Put label EXCITLAB in list for two-index vectors
C for vectors of the type N^XX, N^BX and N^CX
C
      DATA CTPLB/'EXCITLAB'/
      DATA DTPLB/'EXCITLAB'/
C
C If excited state alpha (TPALP) is specified we only compute certain
C components.
C
      IF (TPALP) THEN
C
      DO 600 ICSYM = 1,NSYM
      DO 620 IBSYM = 1,NSYM
         IDSYM = ICSYM
         IASYM = IBSYM
         IF ( (NTPCN2(IDSYM).GT.0) .AND. (NTPCN1(ICSYM).GT.0) .AND.
     *    (NBTPOP(IBSYM).GT.0) .AND. (NATPOP(IASYM).GT.0) ) THEN
            DO IBOP = 1,NBTPOP(IBSYM)
               DO ICFR = 1,NTPCN1(ICSYM)
               DO IBFR = 1,NBTPFR
                  IDFR = ICFR
                  CTPFR = -EXCIT2(ICSYM,ICFR)
                  DTPFR = EXCIT2(IDSYM,IDFR)
                  INUM = INCR2(BTPLB(IBSYM,IBOP),CTPLB,
     *                          -BTPFR(IBFR),CTPFR,IBSYM,ICSYM)
                  INUM = INCR2(BTPLB(IBSYM,IBOP),DTPLB,
     *                          -BTPFR(IBFR),DTPFR,IBSYM,IDSYM)
                  INUM = INCR2(CTPLB,DTPLB,
     *                          CTPFR,DTPFR,ICSYM,IDSYM)
               ENDDO
               ENDDO
            ENDDO
         END IF
620   CONTINUE
600   CONTINUE
C
      ELSE
C
      DO 300 IDSYM = 1,NSYM
      DO 200 ICSYM = 1,NSYM
      DO 100 IBSYM = 1,NSYM
         IASYM = MULD2H(IDSYM,MULD2H(ICSYM,IBSYM))
         IF ( (NTPCN2(IDSYM).GT.0) .AND. (NTPCN1(ICSYM).GT.0)
     *     .AND. (NBTPOP(IBSYM).GT.0) .AND. (NATPOP(IASYM).GT.0)) THEN
            DO 110 IBOP = 1,NBTPOP(IBSYM)
               DO 320 IDFR = 1,NTPCN2(IDSYM)
               DO 220 ICFR = 1,NTPCN1(ICSYM)
               DO 120 IBFR = 1,NBTPFR
                  CTPFR = -EXCIT2(ICSYM,ICFR)
                  DTPFR = EXCIT2(IDSYM,IDFR)
                  INUM = INCR2(BTPLB(IBSYM,IBOP),CTPLB,
     *                          -BTPFR(IBFR),CTPFR,IBSYM,ICSYM)
                  INUM = INCR2(BTPLB(IBSYM,IBOP),DTPLB,
     *                          -BTPFR(IBFR),DTPFR,IBSYM,IDSYM)
                  INUM = INCR2(CTPLB,DTPLB,
     *                          CTPFR,DTPFR,ICSYM,IDSYM)
120            CONTINUE
220            CONTINUE
320            CONTINUE
110         CONTINUE
         END IF
100   CONTINUE
200   CONTINUE
300   CONTINUE
C
      END IF
C
      RETURN
      END

