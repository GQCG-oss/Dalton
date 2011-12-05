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
#include <single.h>
C
C     Parallel-routines, by Alf Christian Hennum 
C 
      SUBROUTINE DFTPAR(DMAT,DTRMAT,EXCMAT,NXCMAT,DOERG,DOGRD,DOLND,
     &     DOATR,CORX,CORY,CORZ,WEIGHT,NBUF,GSO,NCNT,DST,
     &     C0,C1,C2,DMAGAO,HES,DOHES,TRPLET,FACDRC,FACVWN,
     &     FACBCK,FACLYP,FACKT,DODRC,DOVWN,DOBCK,DOLYP,DOGGA,DOKT,
     &     WORK,LWORK,IPRINT)
C     
C     T. Helgaker sep 99
C
#include <implicit.h>
#if defined (VAR_MPI)
      INCLUDE "mpif.h"
#endif
#include <priunit.h>
#include <maxorb.h>
#include <infpar.h>
#include <mxcent.h>
#include <pi.h>
#include <dummy.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
      PARAMETER (THRSH = 1.0D-8, THRSH2 = 1.0D-9) 
C
#include <inforb.h>
#include <nuclei.h>
#include <dftinf.h>
#include <dfterg.h>
#include <energy.h>
#include <dftcom.h>
#include <orgcom.h>
#include <maxaqn.h>
#include <symmet.h>
C
      INTEGER A, B
      LOGICAL DOERG, DOGRD, DODRC, DOVWN, DOBCK, DOLYP, DOKT, 
     &        DOLND, DOATR, FROMVX, DOGGA, TRPLET, DOHES, FINNISH 
      DIMENSION CORX(NBUF), CORY(NBUF), CORZ(NBUF), WEIGHT(NBUF), 
     &          DMAT(NBAST,NBAST), DTRMAT(NBAST,NBAST),
     &          GSO(NBAST*NTYPSO), EXCMAT(NBAST,NBAST,NXCMAT), 
     &          NCNT(NBAST), DST(NATOMS), C0(NORBT), C1(NORBT,3), 
     &          C2(NORBT),DMAGAO(NBAST), HES(NOCCT,NVIRT,NOCCT,NVIRT), 
     &          RHG(3), WORK(LWORK)     
C
C     start up the nodes:  
C
      IPRTYP = 4
      FINNISH = .FALSE.      
#if defined (VAR_MPI)
      CALL MPI_BCAST(IPRTYP,1,my_MPI_INTEGER,MASTER,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IPRINT,1,my_MPI_INTEGER,MASTER,
     &     MPI_COMM_WORLD,IERR)
C
C     transfer data for memory allocation
C           
      CALL MPI_BCAST(NBAST,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(NBUF,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(NATOMS,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(NXCMAT,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)      
C
C     Some hessian stuff that's not used. 
C
      CALL MPI_BCAST(NOCCT,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NVIRT,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOHES,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NORBT,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
C     DFT constants : 
C
      CALL MPI_BCAST(WDFTX,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(WDFTC,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(WDFTB,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(WDFTL,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(WDFTK,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
C      
C     Some other nedeed variables :  
C
      CALL MPI_BCAST(NSYM,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOLND,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOATR,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOGRD,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
      NCO = 0
      IF (DOHES) NC0  = NORBT
C     
C     Move everything else: 
C     
      CALL DFTMPIDAT(DMAT,DTRMAT,EXCMAT,NXCMAT,
     &     DOERG,DOVWN,
     &     NCNT,DST,C0,C1,C2,HES,DOHES,TRPLET,
     &     IPRINT)
C     
C     Call the dft-driver :  
C     
      CALL DFTDRV(DMAT,DTRMAT,EXCMAT,NXCMAT,DOERG,DOGRD,DOLND,
     &     DOATR,CORX,CORY,CORZ,WEIGHT,NBUF,GSO,NCNT,DST,
     &     C0,C1,C2,DMAGAO,HES,DOHES,TRPLET,FACDRC,FACVWN,
     &     FACBCK,FACLYP,FACKT,DODRC,DOVWN,DOBCK,DOLYP,DOGGA,
     &     DOKT,WORK,LWORK,IPRINT)
#endif
      RETURN
      END
C
      SUBROUTINE DFT_NODSTR(WORK,LWORK,IPRINT)
C
C     ACH 
C     Takes info from master and allocate memory. 
C
#include <implicit.h>
#if defined (VAR_MPI)
      INCLUDE "mpif.h"
#endif
#include <maxorb.h>
#include <infpar.h>
#include <priunit.h>
#include <inforb.h>      
#include <mxcent.h>
C
#include <nuclei.h>
#include <dfterg.h>
#include <dftinf.h>
#include <energy.h>
#include <maxaqn.h>
#include <symmet.h>
#include <dftcom.h>
#include <orgcom.h>

      PARAMETER (FACTHR = 1.0D-8)
C
      DIMENSION WORK(LWORK) 
C      
      LOGICAL DOHES, DOERG, DOGRD, DOLND, DOATR, TRPLET,
     &        DOGGA, DOVWN, DODRC, DOBCK, DOLYP, DOKT
C     
#if defined (VAR_MPI)
      CALL MPI_BCAST(NBAST,1,my_MPI_INTEGER,MASTER,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBUF,1,my_MPI_INTEGER,MASTER,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NATOMS,1,my_MPI_INTEGER,MASTER,
     &               MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(NXCMAT,1,my_MPI_INTEGER,MASTER,
     &     MPI_COMM_WORLD,IERR)      
C     
C     hessian stuff
C     
      CALL MPI_BCAST(NOCCT,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NVIRT,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOHES,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NORBT,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
C     some needed constants and other things 
C     
      CALL MPI_BCAST(WDFTX,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(WDFTC,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(WDFTB,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(WDFTL,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(WDFTK,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NSYM,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOLND,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOATR,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOGRD,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
C     decide what to do and allocate: as in dftexe: 
C     
      FACDRC = WDFTX
      FACVWN = WDFTC
      FACBCK = WDFTB
      FACLYP = WDFTL
      FACKT  = WDFTK
C
      DODRC = DABS(FACDRC) .GT. FACTHR
      DOVWN = DABS(FACVWN) .GT. FACTHR 
      DOBCK = DABS(FACBCK) .GT. FACTHR 
      DOLYP = DABS(FACLYP) .GT. FACTHR 
      DOKT  = DABS(FACKT)  .GT. FACTHR
      DOGGA = DOBCK .OR. DOLYP .OR. DOKT
C     
      NC0  = 0
      IF (DOHES) NC0  = NORBT
C     
      NOSYM = NSYM.EQ.1
C     
      NDER = 0
      IF (DOGRD) NDER = 1
      IF (DOGGA) THEN
         NDER = NDER + 1 
         IF (DFTPOT) NDER = NDER + 1 
      END IF
      IF (NDER.EQ.0) NTYPSO =  1
      IF (NDER.EQ.1) NTYPSO =  4
      IF (NDER.EQ.2) NTYPSO = 10
      NSO0 = 1
      NSO1 = 2
      NSO2 = 5
      IF (DOLND) THEN
         NTYPSO = NTYPSO + 3
         NSOB   = NTYPSO - 2 
         IF (DOGGA) THEN
            NTYPSO = NTYPSO + 9
            NSOB1  = NTYPSO - 8 
         END IF
      END IF
      KSO0 = (NSO0-1)*NBAST + 1
      KSO1 = (NSO1-1)*NBAST + 1
      KSO2 = (NSO2-1)*NBAST + 1
      KSOB = (NSOB-1)*NBAST + 1
      KSOB1 = (NSOB1-1)*NBAST + 1
C
      IDMAT = 1 
      IDTRMAT = IDMAT    + NBAST*NBAST
      IEXCMAT = IDTRMAT  + NBAST*NBAST      
      KX    =   IEXCMAT  + NBAST*NBAST*NXCMAT
      KY    = KX    + NBUF
      KZ    = KY    + NBUF
      KW    = KZ    + NBUF
      KGSO  = KW    + NBUF
      KCNT  = KGSO  + NTYPSO*NBAST
      KDST  = KCNT  + NBAST
      KC0   = KDST  + NATOMS
      KC1   = KC0   + NC0
      KC2   = KC1   + 3*NC0
      KDGA  = KC2   + NC0
      KLST  = KDGA  + NBAST
      KHES = KLST
CKR      KHES  = KDGA  + NBAST
CKR      KLST  = KHES  + NOCCT*NVIRT*NOCCT*NVIRT
      LWORK  = LWORK - KLST + 1
C
C     Set Zero: only the masters excmat is going to be updated 
C     
      CALL DZERO(WORK(IEXCMAT),NBAST*NBAST*NXCMAT)
C
C     The nodes do not want to write anything !! 
C     
      iprint = -1000
C
C     Get the rest of the data from master: 
C
      CALL DFTMPIDAT(WORK(IDMAT),WORK(IDTRMAT),WORK(IEXCMAT),NXCMAT,
     &     DOERG,DOVWN,
     &     WORK(KCNT),WORK(KDST),WORK(KC0),WORK(KC1),WORK(KC2),
     &     WORK(KHES),DOHES,TRPLET,IPRINT) 
C
C     Loop
C
      CALL DFTDRV(WORK(IDMAT),WORK(IDTRMAT),WORK(IEXCMAT),NXCMAT,
     &            DOERG,DOGRD,DOLND,
     &            DOATR,WORK(KX),WORK(KY),WORK(KZ),WORK(KW),NBUF,
     &            WORK(KGSO),WORK(KCNT),WORK(KDST),
     &            WORK(KC0),WORK(KC1),
     &            WORK(KC2),WORK(KDGA),HES,DOHES,TRPLET,FACDRC,FACVWN,
     &            FACBCK,FACLYP,FACKT,DODRC,DOVWN,DOBCK,DOLYP,DOGGA,
     &            DOKT,WORK(KLST),LWORK,IPRINT)
      LWORK=LWORK+KLST+1
#endif      
      RETURN
      END
C     
      SUBROUTINE DFTMPIDAT(DMAT,DTRMAT,EXCMAT,NXCMAT,
     &     DOERG,DOVWN,
     &     NCNT,DST,C0,C1,C2,HES,DOHES,TRPLET,IPRINT)
C     
C     Moves everything that is needed. Something that is not needed is also 
C     moved, but what ? 
C
#include <implicit.h>
#if defined (VAR_MPI)
      INCLUDE "mpif.h"
#include <priunit.h>
#include <aovec.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <infpar.h>
#include <mxcent.h>
#include <pi.h>
#include <dummy.h>
#include <inforb.h>
#include <nuclei.h>
#include <dfterg.h>
#include <energy.h>
#include <dftcom.h>
#include <orgcom.h>
#include <onecom.h>
#include <shells.h>
#include <symmet.h>
#include <primit.h>
#include <maxmom.h>
#include <xyzpow.h>
#include <pincom.h>
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
      PARAMETER (THRSH = 1.0D-8, THRSH2 = 1.0D-9) 
C     
      LOGICAL TRPLET,DOERG, DOHES,DOVWN
      DIMENSION EXCMAT(NBAST,NBAST,NXCMAT),
     &          DMAT(NBAST,NBAST), DTRMAT(NBAST,NBAST),
     &          NCNT(NBAST), DST(NATOMS),
     &          C0(NORBT), C1(NORBT,3), C2(NORBT),
     &          HES(NOCCT,NVIRT,NOCCT,NVIRT)
C
#include <wrkrsp.h>
#include <lmns.h>
#include <sphtrm.h>
C
C     Things needed by dftgrd: 
C
C     NBASIS, CHARGE(NATOMS),CORD(MXCENT), NUCIND, MAXOPR
C     ISYMAX,PT      
C
      CALL MPI_BCAST(NBASIS,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NUCIND,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ISTBNU,MXCENT,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(MAXOPR,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ISYMAX,3*2,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(PT,8,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(FMULT,8,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(MULT,8,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(CHARGE,NATOMS,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(CORD,MXCENT*3,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
C     
C     Things needed by getaos: 
C
C     KMAX, NHKT(MXSHEL), KHKT(MXSHEL), 
C     KCKT(),SPHR(),
C     NUCO(), NUMCF(), JSTRT(),
C     CENT(MXSHEL,2,3)     
C     ISTEP(MXAQNM),  MVAL(MXAQNM),   NVAL(MXAQNM)      
C     CSP(NCSP), ISPADR(MXQN), maxrep 
C     NAOS(8) ,MAXREP, MAXOPR, IPIND(MXCORB),NHKT(MXSHEL),   KHKT(MXSHEL), ISTBAO(IA)
C     KSTRT(MXSHEL) , ISYMAO, PT(0:7),ISYMAO(MXQN,MXAQN)
C     ksymop
C     
C     INFORMATIONS FOR DIFFERENT DFT-ROUTINES
C     I have lost the overview.
C     
      CALL MPI_BCAST(KSYMOP,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ISTEP,MXAQNM,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ISYMAO,MXQN*MXAQN,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(MAXREP,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ISTBAO,MXSHEL,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NBAS,8,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IBAS,8,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IPIND,MXCORB,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(NAOS,8,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(KSTRT,MXSHEL,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(CSP,NCSP,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ISPADR,MXQN,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)     
      CALL MPI_BCAST(KMAX,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NHKT,MXSHEL,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(KHKT,MXSHEL,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(KCKT,MXSHEL,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(SPHR,MXSHEL,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(NUCO,MXSHEL,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(NUMCF,MXSHEL,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(JSTRT,MXSHEL,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(CENT,2*3*MXSHEL,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)            
C
C     for GETGAO
C     PRICCF, PRIEXP
C
      CALL MPI_BCAST(PRICCF,MXCONT*MXPRIM,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(PRIEXP,MXPRIM,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
C     Some symmetry related things
C     
C     GETGA1
C     JSTA,NUCA,LMNVALUA
c     IPTCNT(3*MXCENT,0:7,2)
C     IPTSYM(MXCORB,0:7)
C     MULD2H(8,8)
C
      CALL MPI_BCAST(MVAL,MXAQNM,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NVAL,MXAQNM,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)

      CALL MPI_BCAST(LVALUA,MXAQN,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(MVALUA,MXAQN,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NVALUA,MXAQN,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
      CALL MPI_BCAST(MULD2H,8*8,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IPTSYM,MXCORB*8,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(IPTCNT,3*MXCENT*8*2,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(JSTA,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NUCA,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
C     The arrays: 
C 
      CALL MPI_BCAST(DMAT,NBAST*NBAST,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DTRMAT,NBAST*NBAST,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DST,NATOMS,MPI_DOUBLE_PRECISION,MASTER,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(NCENT,MXSHEL,my_MPI_INTEGER,MASTER,
     &     MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(C0,NORBT,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
C     Hessian ? 
C
      IF (DOHES) THEN 
         CALL MPI_BCAST(C1,NORBT*3,MPI_DOUBLE_PRECISION,
     &        MASTER,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(C2,NORBT,MPI_DOUBLE_PRECISION,
     &        MASTER,MPI_COMM_WORLD,IERR)
         CALL MPI_BCAST(HES,NOCCT*NVIRT*NOCCT*NVIRT,
     &        MPI_DOUBLE_PRECISION,MASTER,MPI_COMM_WORLD,IERR)
      ENDIF
C
C     THE common blocks
C
C     dfterg 3 stk double prec      
C
      CALL MPI_BCAST(EDFTX,1,MPI_DOUBLE_PRECISION,
     &     MASTER, MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(EDFTC,1,MPI_DOUBLE_PRECISION,
     &     MASTER, MPI_COMM_WORLD,IERR)      
      CALL MPI_BCAST(EDFTY,1,MPI_DOUBLE_PRECISION,
     &     MASTER, MPI_COMM_WORLD,IERR)      
C     dftcom: 
C     All of it is prob not nescecary
      isize =  1
      CALL MPI_BCAST(DFTADD,ISIZE,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(GRDONE,ISIZE,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTRUN,ISIZE,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTPOT,ISIZE,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTORD,ISIZE,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTASC,ISIZE,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTHES,ISIZE,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(RADINT,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ANGINT,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(ANGMIN,1,my_MPI_INTEGER,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(HFXFAC,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTHRI,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTHR0,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DFTHRL,1,MPI_DOUBLE_PRECISION,
     &     MASTER,MPI_COMM_WORLD,IERR)
C     Some booleans
      CALL MPI_BCAST(TRPLET,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOERG,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(DOVWN,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
      CALL MPI_BCAST(HFXFAC,1,MPI_LOGICAL,
     &     MASTER,MPI_COMM_WORLD,IERR)
C
#endif
      RETURN
      END
C
      SUBROUTINE UPDTENE(ELCTRN,ERGDRC,ERGVWN,ERGBCK,
     &     ERGLYP,ERGKT,EXCMAT,NBAST,NXCMAT,WORK,LWORK)
C
C     Ach 
C     Updates the different energies and the exchange mat: 
C     
#include <implicit.h>       
#if defined (VAR_MPI)
      INCLUDE "mpif.h"
#endif
#include <maxorb.h>
#include <infpar.h>
Ca
      DIMENSION EXCMAT(NBAST*NBAST*NXCMAT),WORK(LWORK)
C
#if defined (VAR_MPI)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
C
      CALL MPI_REDUCE(ELCTRN,EELCTRN,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR) 
      CALL MPI_REDUCE(ERGDRC,EERGDRC,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR) 
      CALL MPI_REDUCE(ERGVWN,EERGVWN,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR) 
      CALL MPI_REDUCE(ERGBCK,EERGBCK,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR) 
      CALL MPI_REDUCE(ERGLYP,EERGLYP, 1,MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR) 
      CALL MPI_REDUCE(ERGKT, EERGKT, 1,MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR) 
      CALL MPI_REDUCE(EXCMAT,WORK,NBAST*NBAST*NXCMAT,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR)
      IF (MYNUM .EQ. 0) then 
         ELCTRN=EELCTRN
         ERGDRC=EERGDRC
         ERGVWN=EERGVWN
         ERGBCK=EERGBCK
         ERGLYP=EERGLYP
         ERGKT =EERGKT
         DO I=1,NBAST*NBAST*NXCMAT
            EXCMAT(I)=WORK(I) 
         ENDDO
      ENDIF
#endif
      RETURN
      END
C
      SUBROUTINE UPDTGRD(WORK,LWORK)
C
C     Update gradient 
C     Ach 
C
#include <implicit.h>       
#if defined (VAR_MPI)
      INCLUDE "mpif.h"
#endif
#include <maxorb.h> 
#include <mxcent.h>
#include <nuclei.h>
#include <energy.h>
#include <infpar.h>
      DIMENSION WORK(LWORK)
C     
#if defined (VAR_MPI)      
      CALL MPI_REDUCE(GRADFT,WORK,MXCOOR,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,
     &     0, MPI_COMM_WORLD, IERR)
      DO I=1,MXCOOR
         GRADFT(I) = WORK(I) 
      ENDDO
#endif
      RETURN 
      END









