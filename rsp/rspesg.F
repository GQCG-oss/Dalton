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
#include <pcmlog.h>

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
     &        NODV,DIFDIP,.FALSE.,HFONLY,NCLONE,PCM)

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














