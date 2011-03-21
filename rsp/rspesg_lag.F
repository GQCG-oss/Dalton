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
      SUBROUTINE ESGLAG_AOMAT(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     $                        XINDX,WRK,LWRK)
#include <implicit.h>
#include <dummy.h>
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*),WRK(*)
#include <priunit.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <inftap.h>
#include "rspprp.h"
#include <esg.h>

      KD1MO   = 1
      KD2MO   = KD1MO   + N2ORBX
      KDXMO   = KD2MO   + N2ORBX
      KDXSMO  = KDXMO   + N2ORBX
      KFRMO   = KDXSMO  + N2ORBX
      KTMPAO  = KFRMO   + N2ORBX
      KTMPAOF = KTMPAO  + N2BASX

      KWRK1     = KTMPAOF + NNBASX
      LWRK1     = LWRK    - KWRK1

      IF ( LWRK1 .LT. 0 ) THEN 
        CALL STOPIT('ESG','RSPESG',KWRK1,LWRK)
      END IF 

      CALL ESGLAG_MOMAT(WRK(KD1MO),WRK(KD2MO),
     &                  WRK(KDXMO),WRK(KDXSMO),WRK(KFRMO),
     &        CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX, WRK(KWRK1),LWRK1)

      IF ( IPRRSP .GE. 10 ) THEN
        CALL PRINTMAT('D1MO:          ',1,NORBT,WRK(KD1MO))
        CALL PRINTMAT('D2MO:          ',1,NORBT,WRK(KD2MO))
        CALL PRINTMAT('DXMO:          ',1,NORBT,WRK(KDXMO))
        CALL PRINTMAT('DXSMO:         ',1,NORBT,WRK(KDXSMO))
        CALL PRINTMAT('FRMO:          ',1,NORBT,WRK(KFRMO))
      END IF
C
C  For one-electron part we save folded (triangular) matrices D1AO, FRAO 
C
      CALL GPOPEN(LUESG,'ESG_AOMAT','NEW',' ','UNFORMATTED',
     &                   IDUMMY,.FALSE.)
      REWIND (LUESG)

      WRITE(LUESG) NNBAST
      CALL MO2AO(WRK(KD1MO),WRK(KTMPAO),CMO,WRK(KWRK1),LWRK1)
      CALL FOLD(WRK(KTMPAO),WRK(KTMPAOF),NBAST)
      CALL PKSYM1(WRK(KTMPAOF),WRK(KD1MO),NBAS,NSYM,1)
      CALL WRITT(LUESG,NNBAST,WRK(KD1MO))

      CALL MO2AO(WRK(KFRMO),WRK(KTMPAO),CMO,WRK(KWRK1),LWRK1)
      CALL FOLD(WRK(KTMPAO),WRK(KTMPAOF),NBAST)
      CALL PKSYM1(WRK(KTMPAOF),WRK(KFRMO),NBAS,NSYM,1)
      CALL WRITT(LUESG,NNBAST,WRK(KFRMO))

      CALL GPCLOSE(LUESG,'KEEP')
C
C For two-electron part we save full (square) AO matrices D2AO, DXAO, KDXSAO
C To be read later in TWOEXP
C
      CALL GPOPEN(LUESG2,'ESG_DMAT','NEW',' ','UNFORMATTED',
     &                   IDUMMY,.FALSE.)
      REWIND (LUESG2)

      CALL MO2AO(WRK(KD2MO),WRK(KTMPAO),CMO,WRK(KWRK1),LWRK1)
      CALL WRITT(LUESG2,N2BASX,WRK(KTMPAO))
      CALL MO2AO(WRK(KDXMO),WRK(KTMPAO),CMO,WRK(KWRK1),LWRK1)
      CALL WRITT(LUESG2,N2BASX,WRK(KTMPAO))
      CALL MO2AO(WRK(KDXSMO),WRK(KTMPAO),CMO,WRK(KWRK1),LWRK1)
      CALL WRITT(LUESG2,N2BASX,WRK(KTMPAO))

      CALL GPCLOSE(LUESG2,'KEEP')

      RETURN
      END

C
C  END OF ESGLAG_AOMAT
C
  
      SUBROUTINE ESGLAG_MOMAT(D1MO,D2MO,DXMO,DXSMO,FRMO,
     &        CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK,LWRK)
#include <implicit.h>
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION D1MO(*),D2MO(*),DXMO(*),DXSMO(*),FRMO(*)
      DIMENSION XINDX(*),WRK(*)
      LOGICAL FOUND, CONV
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
#include <esg.h>

C   KXVECS - eigenvectors corresponding to eigenstates
C            of the excited states
C   KKVECS - lagrange multipliers = inv(E2)*TVEC 
C
      KXMATS = 1 
      KKMATS = KXMATS + N2ORBX
      KWRK1  = KKMATS + N2ORBX
      LWRK1  = LWRK - KWRK1

      IF ( LWRK1 .LT. 0 ) THEN 
        CALL STOPIT('ESG','RSPESG',KWRK1,LWRK)
      END IF 

      CALL GETTIM(ESGTIM_0,DUMTIM)        
      CALL GETTIM(ESGTIM_1,DUMTIM)        

c      CALL ESG_XKMATS(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,1,WRK(KXMATS),
c     &             WRK(KKMATS),XINDX,WRK(KWRK1),LWRK1)

      CALL ESGXKMATS_NEW(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,1,WRK(KXMATS),
     &             WRK(KKMATS),XINDX,WRK(KWRK1),LWRK1)

      CALL GETTIM(ESGTIM_2,DUMTIM)        

      CALL ESG_DMOMAT(D1MO,D2MO,DXMO,DXSMO,FRMO,
     &             CMO,UDV,PV,FC,FV,FCAC,H2AC,WRK(KXMATS),
     &             WRK(KKMATS),XINDX,WRK(KWRK1),LWRK1)

      CALL GETTIM(ESGTIM_3,DUMTIM)        

      ESGTIM_KVECS = ESGTIM_2 - ESGTIM_1
      ESGTIM_MOMAT = ESGTIM_3 - ESGTIM_2

      RETURN
      END

      SUBROUTINE ESG_XKMATS(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,NSIM,XMATS,
     &             KMATS,XINDX,WRK,LWRK)
C
#include <implicit.h>
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XMATS(*),KMATS(*),XINDX(*),WRK(*)
      PARAMETER (MXDMAT = 50)
      DIMENSION IFCTYP(MXDMAT), ISYMDM(MXDMAT)
      CHARACTER*8 BLANK
      LOGICAL FOUND, CONV
      PARAMETER ( BLANK = '        ', D0=0.0D0, D1 = 1.0D0 )
#include <priunit.h>
#include <infopt.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <maxorb.h>
#include <infinp.h>
#include <inftap.h>
#include "rspprp.h"
#include <inflr.h>
#include <wrkrsp.h>
#include <esg.h>

#include <qrinf.h>
#include <infvar.h>
C
C     SUBROUTINE FOR CALCUATING THE T-VECTOR FOR 
C     THE LAGRANGE MULTIPLIER IN EXCITED STATE 
C     GRADIENT CALCULATION
C
C     INPUT - excited state eigenvectors X 
C
C     DENSITIES USED
C    
C     D0  = ordinary density for HF ( in MO - 2*delta_{ij} )
C     D1  = X * D0 - D0 * X
C     D2  = X^T * D1 - D1 * X^T 
C
C
C    intermediate two-electron contributions to 
C    Fock matrices used
C
C    FV0_{pq} = D0_{pq} * (g_{pqrs}-0.5g_{psrq})
C    FV1_{pq} = D1_{pq} * (g_{pqrs}-0.5g_{psrq})
C    FV2_{pq} = D2_{pq} * (g_{pqrs}-0.5g_{psrq})
C
C    The generalized Fock matrices are then constructed 
C    as 
C         F1 = D2*FC - D1*FV1 - D1'*FV1' + D0*FV2'
C
C     Resulting T vector (in matrix form) 
C     is than T = F1 - F1^T
C 
      WRITE(LUPRI,*) 'IN ESG_XKMATS NSIM, KZYVAR',NSIM,KZYVAR
C
C     Allocate the workspace for generalized density 
C     matrices, unpacked X vector and generalize Fock
C     matrices
C

      KXVECS   = 1
      KD1MO    = KXVECS + KZYVAR
      KD2MO    = KD1MO  + N2ORBX
      KD1AO    = KD2MO  + N2ORBX
      KD2AO    = KD1AO  + N2ORBX

      KFV1AO   = KD2AO  + N2ORBX 
      KFV2AO   = KFV1AO + N2ORBX 
      KFV1MO   = KFV2AO + N2ORBX 
      KFV2MO   = KFV1MO + N2ORBX 
      KFXXMO   = KFV2MO + N2ORBX 

      KFCDIAG  = KFXXMO + N2ORBX 
      KD0DIAG  = KFCDIAG+ NORBT

      KWRK1    = KD0DIAG+ NORBT
      LWRK1    = LWRK - KWRK1

      IF ( LWRK1 .LT. 0 ) THEN 
        CALL STOPIT('ESG','ESG_XKMATS',KWRK1,LWRK)
      END IF 

C
C     Clear up the place for all the matrices
C
      CALL DZERO(WRK,KWRK1)

C
C READ EIGENVECTORS TO THE EXCITATION ENERGIES FROM FILE
C

      CALL REARSP(LURSP,KZYVAR,WRK(KXVECS),
     &         'EXCITLAB',BLANK, EXCITA(ISYME,IESG,1),D0,
     &          ISYME,0,THCLR, FOUND,CONV,ANTSYM)

      CALL PRINTVEC('XVECS:         ',1,KZYVAR,WRK(KXVECS))
      CALL PRINTVEC2('XVECS:         ',1,KZVAR,WRK(KXVECS))
C
C  Diagonal elements of ordinary one-electron density 
C  matrix and inactive Fock matrix (for HF ground state)
C

      CALL D0FCDIAG(WRK(KD0DIAG),WRK(KFCDIAG),FC)

C
C     Unpack the X-vector to the matrix
C

      CALL RSPZYM(1,WRK(KXVECS),XMATS)
C
C     Construct the one index transformed (by X) one electron 
C     density matrix D1 (in MO basis)
C
C     D1  = X * D0 - D0 * X
C
C
      CALL TRANSFORMED_DENSITY(1,XMATS,WRK(KD1MO))
C
C     Construct the double transformed (by X) one electron 
C     density matrix D2 (in MO basis)
C
C     D2  = X^T * D1 - D1 * X^T C
C
      CALL DGEMM('t','n',NORBT,NORBT,NORBT,1.0D0,XMATS,NORBT,
     &                    WRK(KD1MO),NORBT,0.0D0,WRK(KD2MO),NORBT)
      CALL DGEMM('n','t',NORBT,NORBT,NORBT,-1.0D0, WRK(KD1MO),NORBT,
     &                    XMATS,NORBT,1.0D0, WRK(KD2MO),NORBT)

      IF (IPRRSP.GE.10) THEN    
        CALL PRINTMAT('(ESG_XKMATS)XMAT',1,NORBT,XMATS)
        CALL PRINTMAT('(ESG_XKMATS)D1MO',1,NORBT,WRK(KD1MO))
        CALL PRINTMAT('(ESG_XMATS)D2MO',1,NORBT,WRK(KD2MO))
      END IF

C
C     Transform the generalized density matrices to AO basis
C
      CALL MO2AO(WRK(KD1MO),WRK(KD1AO),CMO,WRK(KWRK1),LWRK1)
      CALL MO2AO(WRK(KD2MO),WRK(KD2AO),CMO,WRK(KWRK1),LWRK1)

      IF (IPRRSP.GE.10) THEN    
        CALL PRINTMAT('D1AO           ',1,NORBT,WRK(KD1AO))
        CALL PRINTMAT('D2AO           ',1,NORBT,WRK(KD2AO))
      END IF
C
C     Calculate the active fock matrices for the generalized 
C     density matrices
C
C     Both in one call (note: corresponding matrices must be consecutive
C     in memory)
C
      NDMAT=2
      ISYMDM(1)=ISYME
      ISYMDM(2)=1
      IFCTYP(1)=03
      IFCTYP(2)=03

      CALL SIRFCK(WRK(KFV1AO), WRK(KD1AO), NDMAT, ISYMDM,
     &                         IFCTYP,DIRFCK,
     &                         WRK(KWRK1),LWRK1)
c      CALL SIRFCK(WRK(KFV2AO), WRK(KD2AO), 1,ISYMDM,
c     &                         IFCTYP,DIRFCK,
c     &                         WRK(KWRK1),LWRK1)

      IF (IPRRSP.GE.10) THEN    
        CALL PRINTMAT('FV1AO          ',1,NORBT,WRK(KFV1AO))
        CALL PRINTMAT('FV2AO          ',1,NORBT,WRK(KFV2AO))
      END IF
C
C    Transform generalized Fock matrices to MO basis
C

      CALL AO2MO(WRK(KFV1AO),WRK(KFV1MO),CMO,WRK(KWRK1),LWRK1)
      CALL AO2MO(WRK(KFV2AO),WRK(KFV2MO),CMO,WRK(KWRK1),LWRK1)

      IF (IPRRSP.GE.10) THEN
        CALL PRINTMAT('FV1MO          ',1,NORBT,WRK(KFV1MO))
        CALL PRINTMAT('FV2MO          ',1,NORBT,WRK(KFV2MO))
        CALL PRINTMAT('FXXMO: (before ',1,NORBT,WRK(KFXXMO))
      END IF
C
C     Construct T vector
C
C
C    first part :  D2*FC we get simply by multiplying the diagonal 
C      elements of FC with corresponding vectors in D2
C      (since FC is diagonal for HF)
C
      IOFF = 0
      DO I=1,NORBT
         CALL DAXPY(NORBT,WRK(KFCDIAG+I-1),
     &        WRK(KD2MO+IOFF),1,WRK(KFXXMO+IOFF),1)
         IOFF=IOFF+NORBT
      END DO
      
      IF (IPRRSP.GE.12) THEN    
         CALL PRINTMAT('FXXMO: step 1.:',1,NORBT,WRK(KFXXMO))
      END IF
C
C    second part : -D1*FV1'
C
      CALL DGEMM('N','T',NORBT,NORBT,NORBT,-1.0D0,WRK(KD1MO),NORBT,
     &           WRK(KFV1MO),NORBT,1.0D0,WRK(KFXXMO),NORBT)

      IF (IPRRSP.GE.12) THEN    
         CALL PRINTMAT( 'FXXMO: step 2.:',1,NORBT,WRK(KFXXMO))
      END IF
C
C  Third part : -D1'*FV1
C
      CALL DGEMM('T','N',NORBT,NORBT,NORBT,-1.0D0,WRK(KD1MO),NORBT,
     &           WRK(KFV1MO),NORBT,1.0D0,WRK(KFXXMO),NORBT)

      IF (IPRRSP.GE.12) THEN    
         CALL PRINTMAT( 'FXXMO: step 3.:',1,NORBT,WRK(KFXXMO))
      END IF
C
C  fourth part: D0*FV2 we get by multiplying diagonal 
C        elements of D0 with the corresponding vectors 
C        in FV2.
C     
      IOFF = 0
      DO ISYM = 1, NSYM
         DO I=1,NISH(ISYM)
            IOFF = IORB(ISYM)+I-1
            CALL DAXPY(NORBT,WRK(KD0DIAG+IOFF),
     &           WRK(KFV2MO+IOFF),NORBT,WRK(KFXXMO+IOFF),NORBT)
         END DO
      END DO

      IF (IPRRSP.GE.12) THEN    
         CALL PRINTMAT( 'FXXMO: step 4.:',1,NORBT,WRK(KFXXMO))              
      END IF 
C
C     Need to reset variables to reference state symmetry for FMAT2VEC
C     to work
C 
      KSYMOP = IREFSY
      CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KWRK1),LWRK1)

      KTVECS = KWRK1 
      KKVECS = KTVECS + KZYVAR 
      KWRK2  = KKVECS + KZYVAR 
      LWRK2  = LWRK1 - KWRK2

      CALL FMAT2VEC(1,WRK(KTVECS),WRK(KFXXMO))
      CALL PRINTVEC2( 'TVECS          ',1,KZVAR,WRK(KTVECS))
c      CALL PRPORB(WRK(KFXXMO),WRK(KFXXMO),WRK(KTVECS))

      KEXSIM = 1
      KEXCNV = KEXSIM
      CALL INVE2VEC(CMO,UDV,PV,FC,FV,FCAC,H2AC,1,WRK(KTVECS),
     &              WRK(KKVECS),XINDX,WRK(KWRK2),LWRK2)

      CALL PRINTVEC2( 'KVECS           ',1,KZVAR,WRK(KKVECS))

      CALL RSPZYM(1,WRK(KKVECS),KMATS)

      RETURN
      END
C 
C     END OF ESG_XKMATS
C
      SUBROUTINE ESGXKMATS_NEW(CMO,UDV,PV,FOCK,FC,FV,FCAC,H2AC,
     $                          NSIM,XMATS,KMATS,XINDX,WRK,LWRK)
C
#include <implicit.h>
#include <dummy.h>
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XMATS(*),KMATS(*),XINDX(*),WRK(*)
      PARAMETER (MXDMAT = 50)
      DIMENSION IFCTYP(MXDMAT), ISYMDM(MXDMAT)
      CHARACTER*8 BLANK
      LOGICAL FOUND, CONV
      PARAMETER ( BLANK = '        ', D0=0.0D0, D1 = 1.0D0 )
#include <priunit.h>
#include <infopt.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <infpri.h>
#include <maxorb.h>
#include <infinp.h>
#include <infspi.h>
#include <inftap.h>
#include "rspprp.h"
#include <inflr.h>
#include <wrkrsp.h>
#include <esg.h>

#include <qrinf.h>
#include <infvar.h>

C     Konstruct T = E[3] X X^T and solve K E[2] = T

C     INPUT - excited state eigenvectors X 
C     OUTPUT - Lagrange multipliers K

      WRITE(LUPRI,*) 'IN ESGXKMATS_NEW NSIM, KZYVAR',NSIM,KZYVAR
      IF (NSIM.GT.1) WRITE(LUPRI,*) 'Warning: NSIM larger than 1'

C     Make sure the symmetry variables in MJWOP have been set (better to
C     set above and pass MJWOP through?)

      KFREE = 1
      LFREE = LWRK
      CALL MEMGET('INTE',KMJWOP,16*MAXWOP,WRK,KFREE,LFREE)
      KSYMOP = IREFSY
      CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KFREE),LFREE)
      call SETZY(WRK(KMJWOP))
      KSYMOP = ISYME
      CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KFREE),LFREE)
      call SETZY(WRK(KMJWOP))

      ISYMT = IREFSY
      ISYMX = ISYME

      KZYVT  = MZYVAR(ISYMT)
      KZYVX  = MZYVAR(ISYMX)
      KZVT  = MZVAR(ISYMT)
      KZVX  = MZVAR(ISYMX)

      CALL MEMGET('REAL',KVECT,KZYVT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECK,KZYVT,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECXN,KZYVX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KVECXT,KZYVX,WRK,KFREE,LFREE)

C      CALL MEMGET('REAL',KFT   ,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KZYMXN,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KZYMXT,N2ORBX,WRK,KFREE,LFREE)

C     Read in eigenvector X and make X^dag which is also needed
cdj  minus sign in triplet case? check ANTSYM

      CALL REARSP(LURSP,KZYVX,WRK(KVECXN),
     &         'EXCITLAB',BLANK, EXCITA(ISYMX,IESG,1),D0,
     &          ISYMX,0,THCLR, FOUND,CONV,ANTSYM)
      CALL DCOPY(KZVX,WRK(KVECXN),1,WRK(KVECXT+KZVX),1)
      CALL DCOPY(KZVX,WRK(KVECXN+KZVX),1,WRK(KVECXT),1)
      CALL RSPZYM(1,WRK(KVECXN),XMATS)
      IF (IPRRSP.GE.12) THEN
         CALL PRINTVEC2('XVECS:         ',1,KZVX,WRK(KVECXN))
         CALL PRINTVEC2('XTVEC:         ',1,KZVX,WRK(KVECXT))
      END IF

C     Do the E[3] X X^T contraction to get T
      
      IBEQC = 0
c      write(lupri,*)  kzyvt,kzyvx,isymt,isymx
c      write(lupri,*)  kvecxn,kvecxt,kvect,kft,kfree,kmjwop
      CALL DZERO(WRK(KVECT),KZYVT)
      CALL E3INIT(WRK(KVECXN),WRK(KVECXT),DUMMY,.FALSE.,IBEQC,
     *            WRK(KVECT),XINDX,UDV,PV,WRK(KFREE),LFREE,
     *           KZYVT,KZYVX,KZYVX,ISYMT,ISYMX,ISYMX,ISPINA,ISPINB,
     *            ISPINC,CMO,FC,FV,WRK(KMJWOP))
C     Scale by -1
      CALL DSCAL(KZYVT,-D1,WRK(KVECT),1)
      CALL PRINTVEC2( 'TVECS E3INIT   ',1,KZVT,WRK(KVECT))

C     Need to reset symmetry variables before solving the multiplier
C     equation

      KSYMOP = ISYMT
      CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KFREE),LFREE)

      KEXSIM = 1
      KEXCNV = KEXSIM
      CALL INVE2VEC(CMO,UDV,PV,FC,FV,FCAC,H2AC,1,WRK(KVECT),
     &              WRK(KVECK),XINDX,WRK(KFREE),LFREE)

C      IF (IPRRSP.GE.10) THEN
         CALL PRINTVEC2( 'KVEC E3INIT     ',1,KZVT,WRK(KVECK))
C      END IF

      CALL RSPZYM(1,WRK(KVECK),KMATS)

      RETURN
      END
C 
C     END OF ESGKVECS_NEW
C

      SUBROUTINE ESG_DMOMAT(D1MO,D2MO,DXMO,DXSMO,FRMO,
     &             CMO,UDV,PV,FC,FV,FCAC,H2AC,XMATS,
     &             KMATS,XINDX,WRK,LWRK)
C
#include <implicit.h>
      DIMENSION D1MO(*),D2MO(*),DXMO(*),DXSMO(*),FRMO(*)
      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XMATS(*),KMATS(*),XINDX(*),WRK(*)
      INTEGER   IFCTYP(1), ISYMDM(1)
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
#include <infsop.h>

C
C   Calculate the excited state energy gradient
C
C   Input : eigenvectors for excited states XVECS
C           Lagrange multipliers KVECS
C
C   Generalized densities used
C
C    D0 = normal HF density 
C    DX  = X   * D0 - D0 * X
C    DXX = X^T * DX - DX * X^T    (^T = transpose)
C    DK  = K   * D0 - D0 * K
C
C    D1  = D0 + Dxx + DK
C    D2  = 0.25*D0 + 0.5*DXX + 0.5*DK 
C
C    Excited state lagrangean than is
C
C    L = D1_pq h_pq 
C      + ( D2_pq*D_rs + D_pq*D2_rs -0.5*(  D2_ps*D_rq + D_ps*D2_rq )
C           - ( DX_pq*DX_rs - 0.5*DX_ps*DX_rq ) )g_pqrs
C 
C    The gradient is calculated with the standard procedure 
C    for the evaluation for derivatives of OMO integrals
C

      KD1    = 1 
      KD2    = KD1   + N2ORBX 
      KDX    = KD2   + N2ORBX 
      KDXX   = KDX   + N2ORBX 
      KDK    = KDXX  + N2ORBX 
      KD0    = KDK   + N2ORBX 
      KD0DIAG= KD0   + N2ORBX 
     
      KFGMO  = KD0DIAG + NORBT

      KWRK1 = KFGMO + N2ORBX
      LWRK1 = LWRK - KWRK1

      IF ( LWRK1 .LT. 0 ) THEN 
        CALL STOPIT('ESG','ESG_DMOMAT',KWRK1,LWRK)
      END IF 

      CALL DZERO(WRK,KWRK1)
C
C     Set the diagonal elements of D0 to 2.0D0
C
      DO ISYM = 1, NSYM
         DO I=1,NISH(ISYM)
            WRK(KD0DIAG+IORB(ISYM) + I - 1) =2.0D0
         END DO
      END DO
C
C     Unpack the X-vector and K-vector 
C     to the matrix form
C
      IF (IPRRSP .GE. 10 ) THEN
        CALL PRINTMAT( 'ESG_MOMAT: XMAT:  ',1,NORBT,XMATS)              
        CALL PRINTMAT( 'ESG_MOMAT: KMAT:  ',1,NORBT,KMATS)              
      END IF
C
C     Construct the transformed (by X and K vectors) one electron 
C     density matrix DX, DK (in MO basis)
C
C     DX  = X * D0 - D0 * X 
C     DK  = K * D0 - D0 * K 
C
      CALL TRANSFORMED_DENSITY(1,XMATS,WRK(KDX))      
      CALL TRANSFORMED_DENSITY(1,KMATS,WRK(KDK))      
C
C     Construct the double transformed (by X) one electron 
C     density matrix D2 (in MO basis)
C
C     DXX  = X^T * DX - DX * X^T C
C
      CALL DGEMM('t','n',NORBT, NORBT, NORBT,1.0D0, XMATS,NORBT,
     &                   WRK(KDX),NORBT,1.0D0, WRK(KDXX),NORBT)
      CALL DGEMM('n','t',NORBT, NORBT, NORBT,-1.0D0, WRK(KDX),NORBT,
     &                   XMATS,NORBT,1.0D0, WRK(KDXX),NORBT)

      IF (IPRRSP.GE.10) THEN    
        CALL PRINTMAT( '(ESG_MOMAT)DX :',1,NORBT,WRK(KDX)) 
        CALL PRINTMAT( '(ESG_MOMAT)DXX:',1,NORBT,WRK(KDXX)) 
        CALL PRINTMAT( '(ESG_MOMAT)DK :',1,NORBT,WRK(KDK)) 
      END IF 
C
C     Construct the generalized density matrices
C
C     D1  = D0 - Dxx - DK
C     D2  = 0.25*D0 - 0.5*DXX - 0.5*DK 
C
      DO ISYM = 1, NSYM
         DO I = 1, NISH(ISYM)
            IOFF = (IORB(ISYM) + I - 1)*NORBT + IORB(ISYM) + I - 1
            WRK(KD1 + IOFF) = WRK(KD0DIAG + IORB(ISYM) + I - 1)
            WRK(KD2 + IOFF) = 0.25D0*WRK(KD0DIAG + IORB(ISYM) + I - 1)
         END DO
      END DO

      CALL DAXPY(NORBT*NORBT,-1.0D0,WRK(KDXX),1,WRK(KD1),1)
      CALL DAXPY(NORBT*NORBT,-1.0D0,WRK(KDK), 1,WRK(KD1),1)

      CALL DAXPY(NORBT*NORBT,-0.5D0,WRK(KDXX),1,WRK(KD2),1)
      CALL DAXPY(NORBT*NORBT,-0.5D0,WRK(KDK), 1,WRK(KD2),1)

      IF (IPRRSP.GE.10) THEN    
        CALL PRINTMAT( '(ESG_MOMAT) D1:',1,NORBT,WRK(KD1))              
        CALL PRINTMAT( '(ESG_MOMAT) D2:',1,NORBT,WRK(KD2))              
      END IF

C 
C   We need to construct the following generalized 
C   Fock Matrix now: 
C
C     FG = 0.5*D0*hMO + 2*D2*F + 2*D0*(FV2^T) - 2*DX*(FX^T)
C

      CALL REORTH_FOCK(CMO,UDV,PV,FC,FV,FCAC,H2AC,1,
     &      WRK(KFGMO),WRK(KDX),WRK(KD2),XINDX,WRK(KWRK1),LWRK1)

      CALL SYMMETRIZE(WRK(KDX),DXSMO,1,NORBT)

      CALL DCOPY(N2ORBX,WRK(KD1),  1,D1MO,1)
      CALL DCOPY(N2ORBX,WRK(KD2),  1,D2MO,1)
      CALL DCOPY(N2ORBX,WRK(KDX),  1,DXMO,1)
      CALL DCOPY(N2ORBX,WRK(KFGMO),1,FRMO,1)

      RETURN
      END
C 
C     END OF ESG_DMOMAT
C

      SUBROUTINE TRANSFORMED_DENSITY(NSIM,X,DX)
#include <implicit.h>
#include <priunit.h> 
#include <inforb.h> 
C
C     Construct the one index transformed (by X) one electron 
C     density matrix DX (in MO basis)
C
C     DX  = X * D0 - D0 * X
C
C     With D0 (in MO) for optimized HF wavefunction this gives
C
C     DX(i,a) = -2 X(i,a) 
C     DX(a,i) = 2 X(a,i)
C

      DIMENSION X(*),DX(*)

      CALL DZERO(DX,NSIM*N2ORBX)

      DO 300 ISIM=1,NSIM
       DO 200 ISYM=1,NSYM
        DO 100 I=1,NISH(ISYM)
         KOFF=I+IORB(ISYM)+(ISIM-1)*N2ORBX
         CALL DAXPY(NORBT,-2.0D0,X(KOFF),NORBT,DX(KOFF),NORBT)

         KOFF=1+(IORB(ISYM)+I-1)*NORBT+(ISIM-1)*N2ORBX
         CALL DAXPY(NORBT,2.0D0,X(KOFF),1,DX(KOFF),1)
 100    CONTINUE
 200   CONTINUE
 300  CONTINUE

      RETURN
      END
C
C     END OF TRANSFORMED DENSITY
C

      SUBROUTINE REORTH_FOCK(CMO,UDV,PV,FC,FV,FCAC,H2AC,NSIM,
     &            FGMO,DXMO,D2MO,XINDX,WRK,LWRK)
C
#include <implicit.h>
#include <priunit.h>
#include <infopt.h>
#include <inforb.h>
#include <maxorb.h>
#include <infinp.h>
#include <infrsp.h>
#include "rspprp.h"
#include <esg.h>

      DIMENSION CMO(*),UDV(*),PV(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION FGMO(*),DXMO(*),D2MO(*),XINDX(*),WRK(*)
      INTEGER   IFCTYP(2), ISYMDM(2)

      KD0DIAG = 1 
      KFCDIAG = KD0DIAG +  NORBT 
      KH1MO   = KFCDIAG +  NORBT 
      KFRMO   = KH1MO   + N2ORBX 
      KFXAO   = KFRMO   + N2ORBX 
      KF2AO   = KFXAO   + N2ORBX 
      KF2MO   = KF2AO   + N2ORBX 
      KFXMO   = KF2MO   + N2ORBX 
      KDXAO   = KFXMO   + N2ORBX 
      KD2AO   = KDXAO   + N2ORBX 

      KWRK1   = KD2AO   + N2ORBX 
      LWRK1   = LWRK - KWRK1 
C
C     Clear up the place for all the matrices
C
      CALL DZERO(WRK,KWRK1)

      IF (IPRRSP.GE.10) THEN    
        CALL PRINTMAT('DXMO:          ',1,NORBT,DXMO)
        CALL PRINTMAT('D2MO           ',1,NORBT,D2MO)
      END IF 

      CALL D0FCDIAG(WRK(KD0DIAG),WRK(KFCDIAG),FC)

      CALL MO2AO(DXMO,WRK(KDXAO),CMO,WRK(KWRK1),LWRK1)
      CALL MO2AO(D2MO,WRK(KD2AO),CMO,WRK(KWRK1),LWRK1)

      IF (IPRRSP .GE. 10 ) THEN
        CALL PRINTMAT('DXAO           ',1,NORBT,WRK(KDXAO))
        CALL PRINTMAT('D2AO           ',1,NORBT,WRK(KD2AO))
      END IF
     
      CALL GETH1MO(WRK(KH1MO),CMO,WRK(KWRK1),LWRK1)

      IF (IPRRSP .GE. 10 ) THEN
        CALL PRINTMAT('H1MO           ',1   ,NORBT,WRK(KH1MO))
      END IF 

      NDMAT = 2
      ISYMDM(1)=ISYME
      ISYMDM(2)=1
      IFCTYP(1)=03
      IFCTYP(2)=03

      CALL SIRFCK(WRK(KFXAO), WRK(KDXAO), NDMAT,ISYMDM,
     &                         IFCTYP,DIRFCK,
     &                         WRK(KWRK1),LWRK1)

      IF (IPRRSP .GE. 10 ) THEN
        CALL PRINTMAT('FXAO           ',1,NORBT,WRK(KFXAO))
        CALL PRINTMAT('F2AO           ',1,NORBT,WRK(KF2AO))
      END IF 

      CALL AO2MO(WRK(KFXAO),WRK(KFXMO),CMO,WRK(KWRK1),LWRK1)
      CALL AO2MO(WRK(KF2AO),WRK(KF2MO),CMO,WRK(KWRK1),LWRK1)

      IF (IPRRSP .GE. 10 ) THEN
        CALL PRINTMAT('FXMO           ',1,NORBT,WRK(KFXMO))
        CALL PRINTMAT('F2MO           ',1,NORBT,WRK(KF2MO))
      END IF 

C
C  Construct the final Fock matrix for the reorthonormalization 
C  contribution
C
C    FGMO = 0.5*D0*hMO + 2*D2*FC + 2*D0*F2^T - DX^T*FX - DX*FX^T 
C

      IOFF = 0
      DO ISYM = 1, NSYM
         DO I=1,NISH(ISYM)
            IOFF = (IORB(ISYM)+I-1)*NORBT
            CALL DAXPY(NORBT,0.5D0*WRK(KD0DIAG+IORB(ISYM)+I-1),
     &             WRK(KH1MO +IOFF),1,FGMO(IORB(ISYM)+I),NORBT)
         END DO
      END DO

      IF (IPRRSP .GE.12 ) THEN
        CALL PRINTMAT('FGMO:  step 1.:',1,NORBT,FGMO)
      END IF

      IOFF = 0
      DO I=1,NORBT
        CALL DAXPY(NORBT,2.0D0*WRK(KFCDIAG+I-1),
     &             D2MO(I),NORBT,FGMO(1+IOFF),1)
        IOFF=IOFF+NORBT
      END DO

      IF (IPRRSP .GE.12 ) THEN
        CALL PRINTMAT('FGMO:  step 2.:',1,NORBT,FGMO)              
      END IF 

      IOFF = 0
      DO ISYM = 1, NSYM
         DO I=1,NISH(ISYM)
            IOFF = (IORB(ISYM)+I-1)*NORBT
            CALL DAXPY(NORBT,2.0D0*WRK(KD0DIAG+IORB(ISYM)+I-1),
     &           WRK(KF2MO+IOFF),1,FGMO(IORB(ISYM)+I),NORBT)
         END DO
      END DO

      IF (IPRRSP .GE.12 ) THEN
        CALL PRINTMAT('FGMO:  step 3.:',1,NORBT,FGMO)              
      END IF

      CALL DGEMM('N','T',NORBT,NORBT,NORBT,1.0D0,DXMO,NORBT,
     &              WRK(KFXMO),NORBT,1.0D0,FGMO,NORBT)
      CALL DGEMM('T','N',NORBT,NORBT,NORBT,1.0D0,DXMO,NORBT,
     &              WRK(KFXMO),NORBT,1.0D0,FGMO,NORBT)

      IF (IPRRSP .GE.10 ) THEN 
        CALL PRINTMAT('FGMO:  step 4.:',1,NORBT,FGMO)              
      END IF

      RETURN
      END
C
C     END OF REORT_FOCK
C
