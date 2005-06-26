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
      SUBROUTINE RSPESG(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
C
C  Purpose:
C     CONTROL CALCULATION OF EXCITED STATE GRADIENT
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)

      LOGICAL PROPTY
C
C
#include <codata.h>
C
#include <priunit.h>
#include <infopt.h>
#include <infrsp.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include <infpp.h>
#include <inflr.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <inftap.h>
#include <mxcent.h>
#include <nuclei.h>
#include <energy.h>
#include <past.h>
#include <gnrinf.h>
#include <esg.h>
#include <taymol.h>
#include <abainf.h>

      CALL HEADER('Excited state gradient calculation',-1)
      CALL GETTIM(ESGTIM_0,DUMTIM)        
      WRITE(LUPRI,'(A,I3)') 'ESG calculation for state number:', IESG
      WRITE(LUPRI,'(A,I3,/)') 'Symmetry of ESG state           :', ISYME

C     ================================================================
C     ***  Calculate and save all the generalised matrices needed  ***
C     ***  for the gradient calculation                            ***
C     ================================================================
C
C     The response vectors for the excited states are calculated 
C     already. Here they are only read from file and the apropriate
C     matrices are constructed.
C
      CALL HEADER('Generalised density matrices for ES lagrangean',-1)

      CALL ESGLAG_AOMAT(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)

      CALL GETTIM(ESGTIM_1,DUMTIM)   
   
      ESGTIM_LAG = ESGTIM_1 - ESGTIM_0

C     =================
C     Initialize Abacus 
C     =================
C

      IPRUSR = 0
      IPRINT = IPRDEF
      IPRESG = IPRINT 
      CALL ABAINP('**PROPE',WRK,LWRK)
      CALL ONEINI
      CALL SETDCR('ABACUS')

      MAXDIF = 1
      MOLGRD = .TRUE. 
      PROPTY = .TRUE.
      ESG    = .TRUE.

      KCSTRA = 1
      KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
      KWRK2  = KSCTRA + 9*NUCDEP*NUCDEP
      LWRK2  = LWRK - KWRK2 + 1 

      CALL HEADER('Calculation of excited state gradient',-1)

      CALL NUCREP(WRK,WRK(MXCOOR*MXCOOR+1),WRK(2*MXCOOR*MXCOOR+1))

      CALL DZERO(GRDMOL,3*NUCDEP)

      CALL GETTIM(ESGTIM_2,DUMTIM)        


C     =========================================================
C     get the one-electron terms ( + reorhonormalization term )
C     =========================================================
C

      CALL ONEDRV(
     &        WRK(KWRK2),LWRK2,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,
     &        NODV,DIFDIP,.FALSE.,HFONLY,NCLONE)

      CALL GETTIM(ESGTIM_3,DUMTIM)        

C     =========================
C     get the two-electron term
C     =========================

      CALL TWOEXP(WRK(KWRK2),LWRK2,PASTWO)

      CALL GETTIM(ESGTIM_4,DUMTIM)        
      ESGTIM_ONE = ESGTIM_3 - ESGTIM_2 
      ESGTIM_TWO = ESGTIM_4 - ESGTIM_3

      CALL ADDGRD(GRADNN)
      CALL ADDGRD(GRADNA)
      CALL ADDGRD(GRADKE)
      CALL ADDGRD(GRADEE)
      CALL ADDGRD(GRADFS)

C     =========
C     printouts
C     =========
C

       IF ( IPRRSP .GE. 4 ) THEN
         CALL HEADER('Molecular gradient (au)',-1)
         WRITE (LUPRI,'(A,I3,/)') 'For excited state number : ', IESG 

         CALL HEADER('Molecular gradient (au) - KE',-1)
         CALL PRIGRD(GRADKE,WRK(KCSTRA),WRK(KSCTRA))
         CALL HEADER('Molecular gradient (au) - NA',-1)
         CALL PRIGRD(GRADNA,WRK(KCSTRA),WRK(KSCTRA))
         CALL HEADER('Molecular gradient (au) - NN',-1)
         CALL PRIGRD(GRADNN,WRK(KCSTRA),WRK(KSCTRA))
         CALL HEADER('Molecular gradient (au) - FS',-1)
         CALL PRIGRD(GRADFS,WRK(KCSTRA),WRK(KSCTRA))
         CALL HEADER('Molecular gradient (au) - EE',-1)
         CALL PRIGRD(GRADEE,WRK(KCSTRA),WRK(KSCTRA))
       END IF

      CALL HEADER('Excited state gradient (au) ',-1)
      WRITE (LUPRI,'(A,I3,/)') 'For excited state number : ', IESG
      CALL PRIGRD(GRDMOL,WRK(KCSTRA),WRK(KSCTRA))

      CALL FLSHFO(LUPRI)
 
      CALL HEADER('Timings for excited state calculation ',-1)
      WRITE (LUPRI,1000) 
     &        ' LAGRANGEAN : ',ESGTIM_LAG,
     &        '    - XVECS : ',ESGTIM_XVECS, 
     &        '    - KVECS : ',ESGTIM_KVECS, 
     &        '    - MOMAT : ',ESGTIM_MOMAT, 
     &        ' ONEINT     : ',ESGTIM_ONE, 
     &        ' TWOINT     : ',ESGTIM_TWO 

 1000 FORMAT(6(/,2X,A,F8.2,' seconds '))

      RETURN
      END

C
C  END OF RSPESG
C














