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
      PROGRAM NODPRG_DALTON
C
C     Node program (slave program) for PVM
C
#include <implicit.h>
#include <priunit.h>
#include <mtags.h>
      PARAMETER (MWORK = 10 000 000)
      DIMENSION WORK(MWORK)
C
C     Join PVM and return my name to master
C
      CALL PVMXMYTID(MYTID)
      CALL PVMXPARENT(MASTER)
      CALL SENDNAME(MASTER,MYTID)
C
C     Get calculation type (iprtyp)
C
      CALL PVMXRECV(MASTER,MTAG1)
      CALL PVMXUNPACK(IPRTYP,1,'INTEGER')
C
      DO WHILE(IPRTYP.GT.0)         
         CALL PVMXUNPACK(IPRINT,1,'INTEGER')
C
         IF (IPRTYP.EQ.1) THEN
            CALL HER_NODSTR(WORK,MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.2) THEN
            CALL ERI_NODSTR(WORK,MWORK,IPRINT)
C         ELSE IF (IPRTYP.EQ.3) THEN
C            CALL GEMNOD(WORK,MWORK,MASTER,MYTID,MYNUM,IPRINT)
         END IF
C
         CALL PVMXRECV(MASTER,MTAG1)
         CALL PVMXUNPACK(IPRTYP,1,'INTEGER')
      END DO
      CALL PVMXEXIT
      STOP
      END
