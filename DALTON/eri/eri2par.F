!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2016 (2015), see http://daltonprogram.org"
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
C
C  /* Deck pardrv */
      SUBROUTINE PARDRV(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,WORK,NSTAT,
     &                  HESSEE,LWORK,ITYPE,MAXDIF,IATOM,NODV,NOPV,
     &                  NOCONT,TKTIME,RETUR,IPRINT,IPRTYP,
     &                  ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,DINTSKP,GENCTR)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mtags.h"
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,GENCTR
      DIMENSION FMAT(*),DMAT(*),ISYMDM(*),IFCTYP(*),WORK(LWORK),
     &          NSTAT(NODTOT), HESSEE(*)
CTROND
      DIMENSION GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*)
CTROND
C
C Used from common blocks
C
C MXORB  : MXSHEL (for INFPAR)
C INFPAR : MASTER, NODTOT, IPRPAR,  NCODE, DBG
C GNRINF : BASDIR, LENBAS
C
#include "infpar.h"
#include "gnrinf.h"
      DATA TOTWAL /0./
      SAVE TOTWAL
C
      CALL QENTER('PARDRV')
C
!     IPRPAR = 11
      IPRINT = MAX(IPRINT,IPRPAR)
C
      IF (IPRINT .GT. 3) THEN
         CALL TITLER('Output from PARDRV','*',103)
         WRITE(LUPRI,'(/A,I5)')
     &      'Parallel type (ITYPE) is',ITYPE
      END IF
      IF (IPRPAR .GT. 0) CALL GETTIM(CPU1,WALL1)
C
      IF (INFPAR_DEBUG) THEN
         IPRNOD = IPRINT
      ELSE
         IPRNOD = 0
      END IF
C
C     Tell nodes whether this is an ERI or HERMIT run. Send IPRNOD (PVM: NCODE)
C
#if defined (VAR_MPI)
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
      CALL MPIXBCAST(IPRNOD,1,'INTEGER',MASTER)
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(IPRTYP,1,'INTEGER')
      CALL PVMXPACK(IPRNOD,1,'INTEGER')
      CALL PVMXMCAST(NODTOT,NODEID(1),MTAG1)
#endif
      IF (IPRTYP.EQ.1) THEN
         if (iprint .gt. 3)
     &      write (lupri,*) 'Calling parallel HERMIT routines.'
         CALL HER_PARDRV(FMAT,DMAT,HESSEE,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                   WORK,LWORK,ITYPE,MAXDIF,IATOM,NODV,
     &                   NOPV,NOCONT,TKTIME,RETUR,IPRINT,
     &                   ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,
     &                   DINTSKP,RELCAL,GENCTR)
      ELSE
         if (iprint .gt. 3)
     &      write (lupri,*) 'Calling parallel ERI routines.'
         CALL ERI_PARDRV(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                   WORK,LWORK,IPRINT)
      END IF
C
      KCPUS  = 1
      KWALLS = KCPUS  + NODTOT
      KLAST  = KWALLS + NODTOT
      IF (KLAST .GT. LWORK)
     &   CALL STOPIT('PARDRV','RVTIM',KLAST,LWORK)
C
      CALL RVTIM(NSTAT,WORK(KCPUS),WORK(KWALLS),NODNAM(1),NODTOT,IPRINT)
C
      IF (IPRPAR .GT. 0) THEN
         CALL GETTIM(CPU2,WALL2)
         CPU    = CPU2 - CPU1
         WALL   = WALL2 - WALL1
         TOTWAL = TOTWAL + WALL
C
         IMINS  = INT(NINT(TOTWAL)/60)
         IHOURS = INT(IMINS/60)
         IMINS  = IMINS - 60*IHOURS
         ISECS  = NINT(TOTWAL) - 3600*IHOURS - 60*IMINS
C
         WRITE(LUPRI,'(A,I5.2,A,I2.2,A,I2.2)')
     &     ' Total wall time used in PARDRV so far   :',
     &      IHOURS,':',IMINS,':',ISECS
         IF (IPRPAR .GT. 2) THEN
            WRITE(LUPRI,'(A,F11.2,/A,F11.2)')
     &        ' CPU  time used in PARDRV last iteration :',
     &          CPU,
     &        ' Wall time used in PARDRV last iteration :',
     &          WALL
         END IF
      END IF
C
      CALL QEXIT('PARDRV')
      RETURN
      END
C  /* Deck rvtim */
      SUBROUTINE RVTIM(NSTAT,CPUS,WALLS,NODNAM,NODTOT,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "mtags.h"
      CHARACTER*20 NAME, NODNAM(NODTOT)
      DIMENSION NSTAT(NODTOT),CPUS(NODTOT),WALLS(NODTOT)
      DIMENSION TMGARR(2)
C
      SAVE TSOFAR
      DATA TSOFAR/0.0D0/
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from RVTIM','*',103)
C

      DO 100 I = 1,NODTOT
#if defined (VAR_MPI)
         NWHO = -1
         CALL MPIXRECV(TMGARR,2,'DOUBLE',NWHO,MTAG8)
#elif defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG8)
         CALL PVMXUNPACK(NWHO  ,1,'INTEGER')
         CALL PVMXUNPACK(TMGARR,2,'DOUBLE')
#else
         NWHO = I
#endif
         CPUS(NWHO)   = TMGARR(1)
         WALLS(NWHO)  = TMGARR(2)
         TSOFAR = TSOFAR + CPUS(NWHO)
C
  100 CONTINUE
C
      IF (IPRINT .LE. 0) GOTO 9000
      IF (IPRINT .GT. 1) THEN
         CALL AROUND('Overall statistics for distribution of batches')
C
         WRITE(LUPRI,'(9X,5(A,4X))')   'Nodename',
     &                                 'Batches',
     &                                 'CPU time',
     &                                 'Wall time',
     &                                 'Efficiency'
         WRITE(LUPRI,'(9X,5(A,4X)/)')  '--------',
     &                                 '-------',
     &                                 '--------',
     &                                 '---------',
     &                                 '----------'
C
         DO 200 NWHO = 1, NODTOT
            IF (WALLS(NWHO) .NE. 0.0D0) THEN
               EFFI = (CPUS(NWHO)/WALLS(NWHO))*100
            ELSE
               EFFI = 0.0D0
            END IF
            WRITE(LUPRI,'(9X,A10,I7,6X,F8.2,5X,F8.2,6X,F6.2)')
     &           NODNAM(NWHO),NSTAT(NWHO),CPUS(NWHO),WALLS(NWHO),EFFI
 200     CONTINUE
      END IF
C
      IMINS  = INT(NINT(TSOFAR)/60)
      IHOURS = INT(IMINS/60)
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(TSOFAR) - 3600*IHOURS - 60*IMINS
C
      WRITE(LUPRI,'(/A,I5.2,A,I2.2,A,I2.2)')
     &     ' Total CPU  time used in NODES  so far   :',
     &     IHOURS,':',IMINS,':',ISECS
      CALL FLSHFO(LUPRI)
C
 9000 RETURN
      END
C  /* Deck eri_pardrv */
      SUBROUTINE ERI_PARDRV(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                      WORK,LWORK,IPRINT)
C
C     **********************************************************
C     *    Driver routine for the parallel ERI calculation.    *
C     **********************************************************
C
C
C      A short survey of the different messagetags (MTAGX) :
C
C      10 - Tell nodes this is an ERI run, send NCODE and IPRINT
C      30 - Send initialization, new input information and DMAT
C      40 - Receive NODCLS and NTASK from node 1
C      50 - Receive (I|J) request from node
C      60 - Send new (I|J) labels to node
C      70 - Receive fock matrix from nodes if MEMOK
C      80 - Receive overall timings
C
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "iratdef.h"
C
      LOGICAL FIRST
      DIMENSION FMAT(*), DMAT(*), ISYMDM(*), IFCTYP(*),
     &          NSTAT(NODTOT), WORK(LWORK)
C
C Used from common blocks
C
C MXORB  : MXSHEL (for INFPAR)
C GNRINF : NEWBAS
C SIRIPC : NEWGEO
C INFORB : N2BASX
C INFPAR : NODEID(), NODTOT, NCODE, NTASK
C
#include "gnrinf.h"
#include "siripc.h"
#include "inforb.h"
#include "infpar.h"
C
      SAVE FIRST
      DATA FIRST/.TRUE./
C
      CALL QENTER('ERI_PARDRV')
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_PARDRV','*',103)
C
C     Consistancy check
C
      CALL ERI_CONCHK(NDMAT,ISYMDM,IFCTYP)
C
C     Initialization in ER2INI (Set /CBIERI/)
C
      CALL ER2INI
C
C     Send initialization to nodes
C
      CALL ERI_SDINIT(DMAT,ISYMDM,NDMAT,IFCTYP,FIRST,IPRINT)
C
C     Receive number of OD classes from node 1. (Letting Master
C     calculate it is waste of CPU time). Also get NTASK
C
      CALL ERI_RVCLS(NODCLS,IPRINT)
C
C     Workspace allocation for parallel ERI calculation
C     =================================================
C
C     Number of (class|class) integrals
C
      MXCLS2 = NODCLS*(NODCLS + 1)/2
C
      IF (MXCLS2 .GT. MAXTSK) THEN
         WRITE(LUPRI,'(/1X,A,/1X,A,/1X,2(A,I5))')
     &      'Number of actual (class1|class2) integrals too large.',
     &      'Recompilation necessary!',
     &      'MXCLS2 =',MXCLS2,' MAXTSK =',MAXTSK
         CALL RECOMPI
         CALL QUIT('ERROR: MXCLS2 greater than MAXTSK in ERI_PARDRV')
      END IF
C
C     Dimension of array INDEX
C
      IBLOCK = NODTOT*NTASK
      NINDEX = (IBLOCK + 1)*INT(MXCLS2/IBLOCK + 1)
C
      KTMPMT = 1
      KINDEX = KTMPMT +  N2BASX
      KTIMES = KINDEX + (NINDEX + 1)/IRAT
      KWHICH = KTIMES +  MXCLS2
      KCLCPU = KWHICH + (MXCLS2 + 1)/IRAT
      KLAST  = KCLCPU +  MXCLS2
C
      IF (KLAST .GT. LWORK) CALL STOPIT('ERI_PARDRV',' ',KLAST,LWORK)
      LWRK = LWORK - KLAST
C
      CALL ERI_PARDR1(FMAT,DMAT,NSTAT,WORK(KTMPMT),WORK(KINDEX),
     &                WORK(KTIMES),WORK(KWHICH),WORK(KCLCPU),
     &                WORK(KLAST),LWRK,NINDEX,MXCLS2,
     &                FIRST,IPRINT)
C
      IF (FIRST)  FIRST  = .FALSE.
      IF (NEWGEO) NEWGEO = .FALSE.
      IF (NEWBAS) NEWBAS = .FALSE.
C
      CALL QEXIT('ERI_PARDRV')
      RETURN
      END
C  /* Deck eri_conchk */
      SUBROUTINE ERI_CONCHK(NDMAT,ISYMDM,IFCTYP)
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION ISYMDM(NDMAT), IFCTYP(NDMAT)
C
CTROND      IF (IFCTYP(1).NE.1) THE
      IF (IFCTYP(1).NE.13) THEN
         WRITE (LUPRI,'(/,A,I5,/A)')
CTROND     &        ' IFCTYP different from 1 or 3 in ERI_CONCHK',
     &        ' IFCTYP different from 13 or 3 in ERI_CONCHK',
     &          IFCTYP, ' Process aborted.'
CTROND         CALL QUIT('IFCTYP .ne. 1 in ERI_CONCHK.')
         CALL QUIT('IFCTYP .ne. 13 in ERI_CONCHK.')
      END IF
      IF (NDMAT .GT. 1) THEN
         WRITE (LUPRI,'(/A,I5,/A)')
     &        ' NDMAT .gt. 1 in ERI_CONCHK.',
     &          NDMAT, ' Process aborted'
         CALL QUIT('NDMAT .gt. 1 in ERI_CONCHK')
      END IF
      IF (ISYMDM(1) .NE. 0) THEN
         WRITE (LUPRI, '(/A,I5,/A)')
     &        ' ISYMDM .ne. 1 in ERI_CONCHK',
     &          ISYMDM(1), ' Process aborted'
         CALL QUIT('ISYMDM .ne. 1 in ERI_CONCHK')
      END IF
C
      RETURN
      END
C  /* Deck eri_sdinit */
      SUBROUTINE ERI_SDINIT(DMAT,ISYMDM,NDMT,IFCTYP,FIRST,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "mtags.h"
C
      LOGICAL FIRST
      DIMENSION DMAT(*), ISYMDM(NDMT), IFCTYP(NDMT)
C
C Used from common blocks
C
C MXCENT : MXCENT (for MOLINP)
C MXORB  : MXSHEL (for INFPAR)
C GNRINF : NEWBAS (PVM: BASDIR, LENBAS)
C SIRIPC : NEWGEO
C MOLINP : NMLINE
C MOLINC : MLINE()
C INFORB : *
C INFPAR : *
C CBIERI : *
C
#include "gnrinf.h"
#include "siripc.h"
#include "molinp.h"
#include "inforb.h"
#include "infpar.h"
#include "r12int.h"
#if defined (SYS_T3D)
C
C     This kludge fixes the current non-implementation of MPI_CHARACTER on
C     the Cray-T3E, K.Ruud-Nov.96
C
      DIMENSION ITMP(4000)
      EQUIVALENCE (ITMP,MLINE)
#endif
#include "cbieri.h"
#include "cbirea.h"
C
      NDMAT = NDMT
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_SDINIT','*',103)
C
#if defined (VAR_MPI)
C
C     Tell nodes whether they will receive new input.
C
      CALL MPIXBCAST(FIRST ,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWGEO,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWBAS,1,'LOGICAL',MASTER)
C
C     Send  /MOLINC/, /MOLINP/, /INFORB/, /INFPAR/ and /CBIERI/
C
      IF (FIRST .OR. NEWBAS .OR. NEWGEO) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &        ' Sending new orbital inf.',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C   /molinp/
         CALL MPIXBCAST(NMLINE,1,'INTEGER',MASTER)
C   /molinc/
#if defined (SYS_T3D)
         CALL MPIXBCAST(ITMP,10*NMLINE,'INTEGER',MASTER)
#else
         CALL MPIXBCAST(MLINE,len_MLINE*NMLINE,'STRING',MASTER)
#endif
C   /inforb/
         CALL MPIXBCAST(MULD2H,NINFI,'INTEGER',MASTER)
C   /infpar/
         CALL MPIXBCAST(NODTOT,NPARI,'INTEGER',MASTER)
C   /cbieri/
         CALL MPIXBCAST(NDMAT ,NCBII,'INTEGER',MASTER)
         CALL MPIXBCAST(RUNERI,NCBIL,'LOGICAL',MASTER)
C   /cbirea/
         CALL MPIXBCAST(LCMMAX,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NCMSTR,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NCMEND,1,'INTEGER',MASTER)
         CALL MPIXBCAST(LCNTNUUM,1,'LOGICAL',MASTER)
C  comr12
         CALL MPIXBCAST(R12INT,NR12L,'LOGICAL',MASTER)
         CALL MPIXBCAST(LMULBS,1,'LOGICAL',MASTER)
         CALL MPIXBCAST(INTGAC,NR12I,'INTEGER',MASTER)
      END IF
      CALL MPIXBCAST(ISYMDM,NDMAT,'INTEGER',MASTER)
      CALL MPIXBCAST(IFCTYP,NDMAT,'INTEGER',MASTER)
      CALL MPIXBCAST(DMAT,N2BASX ,'DOUBLE' ,MASTER)
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(FIRST ,1,'LOGICAL')
      CALL PVMXPACK(NEWBAS,1,'LOGICAL')
      CALL PVMXPACK(NEWGEO,1,'LOGICAL')
      CALL PVMXPACK(BASDIR,60,'STRING')
      CALL PVMXPACK(LENBAS,1,'INTEGER')
C
C     Send  /MOLINC/, /MOLINP/, /INFORB/, /INFPAR/ and /CBIERI/
C
      IF (FIRST .OR. NEWBAS .OR. NEWGEO) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &        ' Sending new orbital inf.',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C   /molinp/
         CALL PVMXPACK(NMLINE,1,'INTEGER')
C   /molinc/
         CALL PVMXPACK(MLINE ,len_MLINE*NMLINE,'STRING')
C   /inforb/
         CALL PVMXPACK(MULD2H,NINFI,'INTEGER')
C   /infpar/
         CALL PVMXPACK(NODTOT,NPARI,'INTEGER')
C   /cbieri/
         CALL PVMXPACK(NDMAT ,NCBII,'INTEGER')
         CALL PVMXPACK(RUNERI,NCBIL,'LOGICAL')
C  comr12
         CALL PVMXBCAST(R12INT,NR12L,'LOGICAL')
         CALL PVMXBCAST(INTGAC,NR12I,'INTEGER')
      END IF
      CALL PVMXPACK(ISYMDM,NDMAT,'INTEGER')
      CALL PVMXPACK(IFCTYP,NDMAT,'INTEGER')
      CALL PVMXPACK(DMAT,N2BASX,'DOUBLE')
      CALL PVMXMCAST(NODTOT,NODEID(1),MTAG3)
#endif
C
      RETURN
      END
C  /* Deck eri_rvcls */
      SUBROUTINE ERI_RVCLS(NODCLS,IPRINT)
C
#include "implicit.h"
#include "maxorb.h"
#include "mtags.h"
C
C Used from common blocks
C
C  MXORB : MXSHEL (for INFPAR)
C  INFPAR: NODEID(1), NTASK
C
#include "infpar.h"
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_RVCLS','*',103)
C
#if defined (VAR_MPI)
      CALL MPIXRECV(NODCLS,1,'INTEGER',NODEID(1),MTAG4)
      CALL MPIXRECV(NTASK ,1,'INTEGER',NODEID(1),MTAG4)
#endif
#if defined (VAR_PVM)
      CALL PVMXRECV(NODEID(1),MTAG4)
      CALL PVMXUNPACK(NODCLS,1,'INTEGER')
      CALL PVMXUNPACK(NTASK ,1,'INTEGER')
#endif
C
      RETURN
      END
C  /* Deck eri_pardr1 */
      SUBROUTINE ERI_PARDR1(FMAT,DMAT,NSTAT,TMPMAT,INDEX,TIMES,
     &                      IWHICH,CLSCPU,WORK,LWORK,NINDEX,
     &                      MXCLS2,FIRST,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mtags.h"
C
Cef Dimension of array INDEX changed from INDEX(NINDEX) to INDEX(3,NINDEX)
      LOGICAL FIRST
      DIMENSION FMAT(*), DMAT(*), NSTAT(NODTOT), TMPMAT(*),
     &          TIMES(MXCLS2), INDEX(3,NINDEX), IWHICH(MXCLS2),
     &          CLSCPU(MXCLS2), WORK(LWORK)
C
C Used from common blocks
C
C MXORB  : MXSHEL (for INFPAR)
C GNRINF : NEWBAS
C SIRIPC : NEWGEO
C INFPAR : NODTOT, NTASK, TIMING, MAXTSK
C
#include "gnrinf.h"
#include "siripc.h"
#include "infpar.h"
Cef begin
#include "incore.h"
Cef end
C
C     NPOS must be saved and therefore statically allocated.
C
      DIMENSION NPOS(0:MAXTSK+1)
      SAVE NPOS
C
      CALL QENTER('ERI_PARDR1')
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_PARDR1','*',103)
C
C     Take calculation time of each (I|J) pair?
C
      TIMING = FIRST .OR. NEWBAS
C
C--------------------------
C     Make the index array.
C--------------------------
C
      CALL HER_INDEKS(TIMING,INDEX,NINDEX,MXCLS2,NODTOT,
     &            NTASK,NBATCH,NPOS,IPRINT)
C
Cef begin
C  10.mai    INITX = .TRUE.
Cef end
C-----------------------
C     Start calculation.
C-----------------------
C
      CALL DOPAR(INDEX,MXCLS2,NBATCH,NSTAT,IPRINT)
C
C--------------------------------------
C     Receive final results from nodes.
C--------------------------------------
C
      CALL ERI_RVRES(FMAT,TMPMAT,MXCLS2,TIMES,IWHICH,CLSCPU,IPRINT)
C
C-----------------------------------------
C     Sort the integral-calculation times.
C-----------------------------------------
C
      IF (TIMING) CALL PARSRT(MXCLS2,NPOS,TMPMAT,TIMES,IPRINT)
C
C     NB! Array TMPMAT is only used as
C         workspace for array SORTED.
C
      CALL QEXIT('ERI_PARDR1')
      RETURN
      END
C     /* Deck dopar */
Cef begin
      SUBROUTINE DOPAR(INDEX,MXBTCH,NBATCH,NSTAT,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mtags.h"
#include "incore.h"
C
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
C
C
      LOGICAL DONE, TDONE
      DIMENSION INDEX(3,MXBTCH), NSTAT(NODTOT), IQUEUE(2,NODTOT,1000),
     &     NQUEUE(NODTOT)
C
C Used from common blocks
C
C  MXORB  : MXSHEL (for INFPAR)
C  INFPAR : NODTOT, NODEID(), NTASK, NCODE
C
#include "infpar.h"
C
C For MPE event logging
C      evIDb1 = MPE_Log_get_event_number()
C      evIDe1 = MPE_Log_get_event_number()
C      call MPE_Describe_state(evIDb1,evIDe1,
C     &     "DOPAR","magenta")
C
C Begin logging the event
C      call MPE_Log_event(evIDb1,0,'')
C
      CALL QENTER('DOPAR')
      IF (IPRINT .GT. 3) CALL TITLER('Output from DOPAR','*',103)
C
      CALL IZERO(NSTAT,NODTOT)
C
      CALL IZERO(IQUEUE,(2*1000*NODTOT))
C
      CALL IZERO(NQUEUE,NODTOT)

C     Start loop over integrals (class1|class2).
C
      IQSIZE = 1000
      IPLACE = 1
      IPLACQ = 0
      ITOTNT = 0
      DONE   = .FALSE.
      TDONE  = .FALSE.
      IDONE = 0
C
      IF (AOSAVE) THEN
C
C         do ii=1,mxbtch
C            write(lupri,*) 'i,INDEX(1..3,i):',ii,INDEX(1,ii),
C     &           INDEX(2,ii),  INDEX(3,ii)
C         end do
         IF (NTASK .GT. 1) THEN
            CALL QUIT('NTASK > 1 is not ported to AOSAVE = TRUE.')
         END IF
C
 110     CONTINUE
C
#if defined (VAR_MPI)
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MTAG5)
         CALL MPIXRECV(IPLAC2,1,'INTEGER',NWHO,MTAG51)
C         CALL MPI_BARRIER(MPI_COMM_WORLD)
         IF (IPLAC2 .NE. 0) THEN
            CALL MPIXRECV(ISAVED,1,'INTEGER',NWHO,MTAG52)
            INDEX(3,IPLAC2) = ISAVED
         END IF
C
C     Send (I|J) batch to NODEID
C
         IF (IPLACE .GT. NBATCH .AND. NQUEUE(NWHO) .EQ. 0) THEN
            TDONE = .TRUE.
            CALL MPIXSEND(TDONE,1,'LOGICAL',NWHO,MTAG6)
            TDONE = .FALSE.
            IDONE = IDONE + 1
         ELSE
            CALL MPIXSEND(DONE,1,'LOGICAL',NWHO,MTAG6)
            IF (NQUEUE(NWHO) .GT. 0) THEN
               CALL MPIXSEND(IQUEUE(1,NWHO,1),NTASK,
     &              'INTEGER',NWHO,MTAG61)
               CALL MPIXSEND(INDEX(3,IQUEUE(2,NWHO,1)),NTASK,'INTEGER',
     &              NWHO,MTAG64)
               CALL MPIXSEND(IQUEUE(2,NWHO,1),1,'INTEGER',
     &              NWHO,MTAG62)
               NQUEUE(NWHO) = NQUEUE(NWHO) - 1
               CALL ICOPY((NQUEUE(NWHO)+1),IQUEUE(1,NWHO,2),(2*NODTOT),
     &              IQUEUE(1,NWHO,1),(2*NODTOT))
               CALL ICOPY((NQUEUE(NWHO)+1),IQUEUE(2,NWHO,2),(2*NODTOT),
     &              IQUEUE(2,NWHO,1),(2*NODTOT))
C
               IPLACQ = IPLACQ - 1
            ELSE IF ((INDEX(3,IPLACE) .EQ. 1
     &          .AND. INDEX(2,IPLACE) .EQ. NWHO)
     &          .OR. (INDEX(3,IPLACE) .EQ. 0)) THEN
C
               CALL MPIXSEND(INDEX(1,IPLACE),NTASK,'INTEGER',NWHO,
     &              MTAG61)
               CALL MPIXSEND(INDEX(3,IPLACE),NTASK,'INTEGER',NWHO,
     &              MTAG64)
C
               INDEX(2,IPLACE) = NWHO
C
               CALL MPIXSEND(IPLACE,1,'INTEGER',NWHO,MTAG62)
               IPLACE = IPLACE + NTASK
            ELSE
C
               NQUEUE(INDEX(2,IPLACE)) = NQUEUE(INDEX(2,IPLACE)) + 1
               IF (NQUEUE(INDEX(2,IPLACE)) .GT. IQSIZE) THEN
                  CALL QUIT('The queue array IQUEUE must be enlarged.')
               END IF
               IQUEUE(1,INDEX(2,IPLACE),NQUEUE(INDEX(2,IPLACE))) =
     &              INDEX(1,IPLACE)
               IQUEUE(2,INDEX(2,IPLACE),NQUEUE(INDEX(2,IPLACE))) =
     &              IPLACE
               IPLACQ = IPLACQ + 1
               IPLACE = IPLACE + NTASK
               ITMP = 0
C     Tell the unemployed node to try again
               CALL MPIXSEND(ITMP,1,'INTEGER',NWHO,MTAG61)
            END IF
         END IF
C
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG5)
         CALL PVMXUNPACK(NWHO,1,'INTEGER')
C
C        Send (I|J) batch to NODEID
C
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(DONE,1,'LOGICAL')
         CALL PVMXPACK(INDEX(1,IPLACE),NTASK,'INTEGER')
         CALL PVMXSEND(NODEID(NWHO),MTAG6)
#endif
C
C        Who-got-the-last-batch statistics
C
         NSTAT(NWHO) = NSTAT(NWHO) + 1
C
         IF (IPRINT .GT. 4) THEN
            WRITE(LUPRI,'(5X,A10,I5,A8,I4)')
     &              'Send batch', I, 'to node', NWHO
         END IF
C
         IF (IPLACE .LE. NBATCH .OR. IPLACQ .GT. 0) THEN
            GOTO 110
         END IF
C
C     Instruct nodes to return the final results.
C
      DONE = .TRUE.
C
      DO 210 I = 1,(NODTOT-IDONE)
C
#if defined (VAR_MPI)
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MTAG5)
C
C        Send (I|J) batch to NODEID
C
         CALL MPIXRECV(IPLAC2,1,'INTEGER',NWHO,MTAG51)
C
         IF (IPLAC2 .NE. 0) THEN
            CALL MPIXRECV(ISAVED,1,'INTEGER',NWHO,MTAG52)
            INDEX(3,IPLAC2) = ISAVED
         END IF
         CALL MPIXSEND(DONE,1,'LOGICAL',NWHO,MTAG6)
C         CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG5)
         CALL PVMXUNPACK(NWHO,1,'INTEGER')
C
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(DONE,1,'LOGICAL')
         CALL PVMXSEND(NODEID(NWHO),MTAG6)
#endif
 210  CONTINUE
C
      IF (MXBTCH .GT. MXTSK) THEN
          CALL QUIT('INDX_C in incore.h must be enlarged.')
       END IF
       CALL ICOPY((3*MXBTCH),INDEX(1,1),1,INDX_C(1,1),1)
C      write(lupri,*) 'INDX_C'
C      DO J=1,MXBTCH
C         write(lupri,*) j,(INDX_C(I,J),i=1,3)
C      END DO
C
C  ELSE below is for the compability to old code: AOSAVE = FALSE
      ELSE
C
         DO 100, I = 1,NBATCH
C
#if defined (VAR_MPI)
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MTAG5)
C
C        Send (I|J) batch to NODEID
C
         CALL MPIXSEND(DONE,1,'LOGICAL',NWHO,MTAG6)
         CALL MPIXSEND(INDEX(MOD(IPLACE-1,3)+1,INT((IPLACE-1)/3)+1),
     &        NTASK,'INTEGER',NWHO,MTAG6)
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG5)
         CALL PVMXUNPACK(NWHO,1,'INTEGER')
C
C     Send (I|J) batch to NODEID
C
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(DONE,1,'LOGICAL')
         CALL PVMXPACK(INDEX(1,IPLACE),NTASK,'INTEGER')
         CALL PVMXSEND(NODEID(NWHO),MTAG6)
#endif
C
C        Who-got-the-last-batch statistics
C
         NSTAT(NWHO) = NSTAT(NWHO) + 1
C
         IF (IPRINT .GT. 4) THEN
            WRITE(LUPRI,'(5X,A10,I5,A8,I4)')
     &              'Send batch', I, 'to node', NWHO
         END IF
C
         IPLACE = IPLACE + NTASK
C
  100 CONTINUE
C
C     Instruct nodes to return the final results.
C
      DONE = .TRUE.
C
      DO 200 I = 1,NODTOT
C
#if defined (VAR_MPI)
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MTAG5)
         CALL MPIXSEND(DONE,1,'LOGICAL',NWHO,MTAG6)
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG5)
         CALL PVMXUNPACK(NWHO,1,'INTEGER')
C
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(DONE,1,'LOGICAL')
         CALL PVMXSEND(NODEID(NWHO),MTAG6)
#endif
  200 CONTINUE
C  ENDIF below ends if-test on AOSAVE
      END IF
C
C End logging the event
C     call MPE_Log_event(evIDe1,0,'')
      CALL QEXIT('DOPAR')
      RETURN
      END
Cef end
C  /* Deck eri_rvres */
      SUBROUTINE ERI_RVRES(FMAT,TMPMAT,MXCLS2,TIMES,IWHICH,CLSCPU,
     &                     IPRINT)
C
#include "implicit.h"
#include "mtags.h"
#include "priunit.h"
#include "maxorb.h"
      PARAMETER (D0=0.0D0,D1 = 1.0D0)
      LOGICAL NOTMEM
      DIMENSION FMAT(N2BASX), TMPMAT(N2BASX), TIMES(MXCLS2),
     &          IWHICH(MXCLS2), CLSCPU(MXCLS2)
C
C Used from common blocks
C
C MXORB  : MXSHEL (for INFPAR)
C INFPAR : NODTOT
C INFORB : N2BASX
C ERIMEM : *
C
#include "infpar.h"
#include "inforb.h"
#include "erimem.h"
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_RVRES','*',103)
C
C     Clean matrices.
C
      CALL DZERO(FMAT,N2BASX)
      CALL DZERO(TIMES,MXCLS2)
C
      IJSTOT = 0
      MAXADD = 0
      NOTMEM = .FALSE.
C
      DO 100 I = 1,NODTOT
#if defined (VAR_MPI)
         NWHO = -1
         CALL MPIXRECV(MEMOK,1,'LOGICAL',NWHO,MTAG7)
         IF (MEMOK) THEN
            CALL DZERO(TMPMAT,N2BASX)
            CALL MPIXRECV(TMPMAT,N2BASX,'DOUBLE',NWHO,MTAG7)
            CALL MPIXRECV(IJS,1,'INTEGER',NWHO,MTAG7)
            IF (TIMING) THEN
               CALL MPIXRECV(IWHICH,IJS,'INTEGER',NWHO,MTAG7)
               CALL MPIXRECV(CLSCPU,IJS,'DOUBLE' ,NWHO,MTAG7)
               DO 200, J=1,IJS
                  TIMES(IWHICH(J)) = CLSCPU(J)
 200           CONTINUE
            END IF
            CALL DAXPY(N2BASX,D1,TMPMAT,1,FMAT,1)
            IJSTOT = IJSTOT + IJS
         ELSE
            CALL MPIXRECV(MEMADD,3,'INTEGER',NWHO,MTAG7)
            CALL MPIXRECV(LWORK ,1,'INTEGER',NWHO,MTAG7)
            IF (MEMADD .GT. MAXADD) THEN
               IABMAX = MODAB
               ICDMAX = MODCD
               IWORK  = LWORK
               MAXADD = MEMADD
            END IF
            NOTMEM = .TRUE.
         END IF
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG7)
         CALL PVMXUNPACK(MEMOK,1,'LOGICAL')
         IF (MEMOK) THEN
            CALL DZERO(TMPMAT,N2BASX)
            CALL PVMXUNPACK(TMPMAT,N2BASX,'DOUBLE')
            CALL PVMXUNPACK(IJS,1,'INTEGER')
            IF (TIMING) THEN
               CALL PVMXUNPACK(IWHICH,IJS,'INTEGER')
               CALL PVMXUNPACK(CLSCPU,IJS,'DOUBLE')
               DO 200, J=1,IJS
                  TIMES(IWHICH(J)) = CLSCPU(J)
  200           CONTINUE
            END IF
            CALL DAXPY(N2BASX,D1,TMPMAT,1,FMAT,1)
            IJSTOT = IJSTOT + IJS
         ELSE
            CALL PVMXUNPACK(MEMADD,3,'INTEGER')
            CALL PVMXUNPACK(LWORK ,1,'INTEGER')
            IF (MEMADD .GT. MAXADD) THEN
               IABMAX = MODAB
               ICDMAX = MODCD
               IWORK  = LWORK
               MAXADD = MEMADD
            END IF
            NOTMEM = .TRUE.
         END IF
#endif
 100  CONTINUE
C
      IF (NOTMEM) THEN
         WRITE (LUPRI,'(//,1X,A,3(/,1X,A,I10))')
     &      ' Not enough memory for this parallel run of ERIFCK.',
     &      ' Available memory in ERIFCK:',IWORK,
     &      ' Required memory for ERIFCK:',IWORK + MAXADD,
     &      ' Increase memory (LWORK) by:',MAXADD
         WRITE (LUPRI,'(/,1X,A,2I5)')
     &      ' Memory requirements largest for OD classes :',
     &        IABMAX,ICDMAX
         CALL QUIT('Insufficient memory in ERIFCK.')
      ELSE
         IF (IPRINT.GT.5) THEN
            CALL HEADER('Fock matrices in ERI_RVRES',-1)
            CALL OUTPUT(FMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL FLSHFO(LUPRI)
         END IF
      END IF
C
      IF (IJSTOT .NE. MXCLS2) THEN
        WRITE(LUPRI,'(/5X,A)') 'Error in parallel calculation!'
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks to calculate  : ',
     &       MXCLS2
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks    calculated : ',
     &       IJSTOT
        CALL QUIT('Inconsistence in ERI_RVRES -> MXCLS2 .NE. IJSTOT')
      ENDIF
C
      RETURN
      END
C  /* Deck dalton_nodedriver */
#if defined (VAR_MPI)
      SUBROUTINE DALTON_NODEDRIVER(MWORK,WRKDLM)
C
C Slaves are sent to this nodedriver from abacus/dalton.F where they 
C idle until the keyword IPRTYP is set with an MPI_BCAST (this is 
C mediated via the subroutine MPIXBCAST from /gp/mpimacro.F)
C Once the slaves have received the IPRTYP.ne.0 broadcast, they need a 
C broadcast with the print level.
C IPRTYP is set in iprtyp.h. This file associates every parallel 
C computation in Dalton with an specific integer and a character string.
C This iprtyp.h file has to be imported in every subroutine that wants to 
C communicate with this nodedriver.
C 
C If you add a parallel routine to Dalton it needs to be added to 
C iprtyp.h
C 
C  This idling works because MPIXBCAST in a collective and blocking 
C communication. Every slave will be blocked by the broadcast at label 
C 100 until the master participates in the broadcast. As soon as an 
C IPRTYP broadcast is made, every slave executes the computation 
C associated with the broadcast once, exits the if-statement and is 
C blocked/idles once again by the MPI_COMM_WORLD broadcast at label 100.
C 
C If IPRTYP is set to 0 the nodes exit the loop.
C 
C As of Aug 2013, you have to broadcast the IPRTYP from the master node 
C using MPIXBCAST or MPI_BCAST, to communicate with the nodedriver. 
C This mean you can change IPRTYP locally if you need it for whatever 
C reason. To start a parallel calculation the master HAS to broadcast 
C the IPRTYP and the needed parameters for the calculation.
C 
C IMPORTANT: The communication is always via the communicator 
C MPI_COMM_WORLD, so if another communicator is made, only those nodes 
C included in the new communicator will receive the broadcast IPRTYP. 
C This is because the nodedriver uses the MPIXBCAST subroutine to 
C communicate and as of Nov 2014 this is hardcoded to MPI_COMM_WORLD.
C 
#ifdef VAR_IFORT
!     subroutine sleep is in IFPORT module
      use IFPORT, ONLY: SLEEPQQ
#endif
#ifdef MOD_CGTO_DIFF_ERI
      use eri_basis_loops
#endif
#ifdef ENABLE_XCINT
      use dalton_xcint_interface, only: xcint_integrate_worker
#endif
      use pelib_interface, only: pelib_ifc_slave

#include "implicit.h"
#include "priunit.h"
#include "mtags.h"
C infpar.h: MASTER,MYNUM
#include "maxorb.h"
#include "infpar.h"
C iprtyp.h: define parallel calculation types
#include "iprtyp.h"
#include "mpif.h"

      real(8), allocatable  :: work(:)
      integer(4), parameter :: TEN_4 = 10

! ======================================================

      CALL QENTER('DALTON_NODEDRIVER')
      CALL GETTIM(CSTR,WSTR)
      N_TASKS = 0
C
C     Get calculation type (IPRTYP)
  100 CONTINUE
#ifdef VAR_IFORT
!     Nov 2013 hjaaj: sleep 10 millisecond before next check for a task
!     - this saves a lot of cpu time on other cores in a sequential section
!     - may make it possible to go into turbo mode
!     - may make it possible to exploit e.g. parallel DGEMM in MKL in serial sections
C      CALL SLEEP(1)
c
!     - only intel ifort has a sleepqq counting in milliseconds (that I know of),
!       the more common sleep(n) sleeps in at least one second, which is too long.
!       (One could use an interface to the C routine usleep.)
      CALL SLEEPQQ(TEN_4)
#endif
      CALL MPIXBCAST(IPRTYP,1,'INTEGER',MASTER)
C
      IF (IPRTYP.GT.0) THEN
         N_TASKS = N_TASKS + 1

         CALL MPIXBCAST(IPRINT,1,'INTEGER',MASTER)
         IF (IPRINT .GT. 1)
     &      write (lupri,'(/A,I5,A,I8,A,I5,A,I5/A)')
     &      'Node',MYNUM,' received task',N_TASKS,' IPRTYP =',IPRTYP,
     &      ' Print level is',IPRINT,
     &      '(See DALTON/include/iprtyp.h for IPRTYP codes)'

!        allocate memory here if we do not call hermit (more to follow...).
!        idea: move allocation down to main driver routines of modules as we might want to use
!              dynamic memory allocation down there (for MPI performance/MPI standard compliance reasons)
!              which is severely hampered by a "static" WORK for everything.
!              sknecht - feb 2013.
         if(iprtyp /= HER_WORK .and. .not.(allocated(work)))then
!          allocate memory using f90 utilities
           allocate(work(0:mwork+1),stat=i)
!          Set memory traps
           work(0)       = wrkdlm
           work(1+mwork) = wrkdlm

           if(i /= 0)then
             write (luerr,*) mynum,
     &           ': ALLOCATE command failed to allocate'//
     &           ' the requested memory on slave. Error code:',i
             write (luerr,*)
     &           'Reduce the memory demands and be welcome back'
             call quit('Failed to allocate memory on slave')
           end if
         end if

C  new calculation type needs to be added in include/iprtyp.h
         IF      (IPRTYP.EQ.HER_WORK) THEN
            CALL HER_NODSTR(MWORK,wrkdlm,IPRINT)
         ELSE IF (IPRTYP.EQ.ERI_WORK) THEN
            CALL ERI_NODSTR(WORK(1),MWORK,IPRINT)
C        ELSE IF (IPRTYP.EQ.GEM_WORK) THEN
C           CALL GEMNOD(WORK(1),MWORK,MASTER,MYTID,MYNUM,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_KSM_WORK) THEN
            CALL DFT_KSMSLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_C_WORK) THEN
            CALL DFT_CSLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_GRAD_WORK) THEN
            CALL DFT_GRAD_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_EXPGRAD_WORK) THEN
            CALL DFT_EXPG_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_BRHS_WORK) THEN
            CALL DFT_BRHS_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_HESSTAT_WORK) THEN
            CALL DFT_HESSTAT_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFTHED_WORK) THEN
            CALL DFTHED_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_SUSCEP_WORK) THEN
             CALL DFT_SUSCEP_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.NPETES2_WORK .OR.
     &            IPRTYP.EQ.PCMBSOL_WORK) THEN
            CALL J1INTS(WORK(1),MWORK,IPRTYP,IPRINT)
         ELSE IF (IPRTYP.EQ.J1XP_WORK) THEN
            CALL J1XS(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.J2XP_WORK) THEN
            CALL J2XS(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.DFT_MOMO_WORK) THEN
            CALL DFT_MOMO_SLAVE(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3_NSP_WORK) THEN
            CALL QM3_NSPS(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3LNO_WORK) THEN
            CALL QM3LNO_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3QRO_WORK) THEN
            CALL QM3QRO_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3FIRST_1_WORK) THEN
            CALL QM3FIRST_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QM3FIRST_2_WORK) THEN
            CALL QM3FIRST_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.PARQMMM__WORK) THEN
            CALL PARQMMM_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.MM_FIELD_1_WORK) THEN
            CALL MM_FIELD_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.MM_FIELD_2_WORK) THEN
            CALL MM_FIELD_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.MM_POLAR_CONTR_WORK) THEN
            CALL MM_POLAR_CONTR_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.MMITER_INDDIP_WORK) THEN
            CALL MMITER_INDDIP_S(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMLNO_1_WORK) THEN
            CALL QMMMLNO_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMLNO_2_WORK) THEN
            CALL QMMMLNO_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMQRO_1_WORK) THEN
            CALL QMMMQRO_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMQRO_2_WORK) THEN
            CALL QMMMQRO_S2(WORK(1),MWORK,IPRINT)
#ifdef MOD_UNRELEASED
         ELSE IF (IPRTYP.EQ.QMMMCRO_1_WORK) THEN
            CALL QMMMCRO_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMCRO_2_WORK) THEN
            CALL QMMMCRO_S2(WORK(1),MWORK,IPRINT)
#endif
         ELSE IF (IPRTYP.EQ.QMMMFIRST_1_WORK) THEN
            CALL QMMMFIRST_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMFIRST_2_WORK) THEN
            CALL QMMMFIRST_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMFIRST_3_WORK) THEN
            CALL QMMMFIRST_S3(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMFIRST_4_WORK) THEN
            CALL QMMMFIRST_S4(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMB2_1_WORK) THEN
            CALL QMMMB2_S1(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMB2_2_WORK) THEN
            CALL QMMMB2_S2(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMB2_3_WORK) THEN
            CALL QMMMB2_S3(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.QMMMB2_4_WORK) THEN
            CALL QMMMB2_S4(WORK(1),MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.LUCITA_WORK) THEN
            call lucita_coworker_main(work(1),mwork,iprint)
            lucita_models_enabled = .true.
         ELSE IF (IPRTYP.EQ.MEP_WORK) THEN
            CALL MEP_SLAVE(WORK(1),MWORK)
         ELSE IF (IPRTYP == POLARIZABLE_EMBEDDING) THEN
            CALL PELIB_IFC_SLAVE(IPRINT)
#if defined(BUILD_GEN1INT)
C added by Bin Gao, May 12, 2012
         ELSE IF (IPRTYP.EQ.GEN1INT_INIT) THEN
            CALL gen1int_worker_init()
         ELSE IF (IPRTYP.EQ.GEN1INT_GET_INT) THEN
            CALL gen1int_worker_get_int(MWORK,WORK(1),LUPRI,IPRINT)
         ELSE IF (IPRTYP.EQ.GEN1INT_GET_EXPT) THEN
            CALL gen1int_worker_get_expt(MWORK,WORK(1),LUPRI,IPRINT)
         ELSE IF (IPRTYP.EQ.GEN1INT_GET_CUBE) THEN
            CALL gen1int_worker_get_cube(MWORK,WORK(1),LUPRI,IPRINT)
#endif
         ELSE IF (IPRTYP.EQ.PARA_RP_ERES) THEN
            CALL PAR_RP_ERES(work(1), mwork, mytid )
         ELSE IF (IPRTYP.EQ.PARA_SO_ERES) THEN
            CALL PAR_SO_ERES(work(1), mwork, mytid )
#ifdef MOD_CGTO_DIFF_ERI
         ! added by Bin Gao, Oct. 27, 2012
         ! cgto-diff-eri code to calculate arbitray order geometric derivatives
         ! of two-electron integrals
         ELSE IF (IPRTYP.EQ.CGTO_DIFF_ERI_INIT) THEN
            CALL unopt_geodiff_loop_worker(LUPRI,IPRINT)
#endif

#ifdef ENABLE_XCINT
         else if (iprtyp == XCINT_MPI_WAKEUP_SIGNAL) then
            call xcint_integrate_worker(MPI_COMM_WORLD)
#endif
         ELSE IF (IPRTYP.EQ.QMCMM_WORK) THEN
#ifdef ENABLE_VPOTDAMP
            CALL vpotdamped_slave()
#else
            call quit('VPOTDAMP not compiled in this version')
#endif

#if defined(ENABLE_EFS)
         else if (IPRTYP.eq.efs_add_atomtype_work) then
            call efs_add_atomtype_slave()
         else if (IPRTYP.eq.efs_init_work) then
            call efs_init_slave()
         else if (IPRTYP.eq.efs_generate_basis_work) then
            call efs_generate_basis_product_slave()
         else if (IPRTYP.eq.efs_init_2efock_work) then
            call efs_init_2efock_slave()
         else if (IPRTYP.eq.efs_fock_update_work) then
            call efs_fock_update_slave()
#endif

         ELSE IF (IPRTYP.EQ.QMCMM_WORK) THEN
#ifdef ENABLE_VPOTDAMP
            CALL vpotdamped_slave()
#else
            call quit('VPOTDAMP not compiled in this version')
#endif

#if defined(ENABLE_EFS)
         else if (IPRTYP.eq.efs_add_atomtype_work) then
            call efs_add_atomtype_slave()
         else if (IPRTYP.eq.efs_init_work) then
            call efs_init_slave()
         else if (IPRTYP.eq.efs_generate_basis_work) then
            call efs_generate_basis_product_slave()
         else if (IPRTYP.eq.efs_init_2efock_work) then
            call efs_init_2efock_slave()
         else if (IPRTYP.eq.efs_fock_update_work) then
            call efs_fock_update_slave()
#endif

         ELSE
            WRITE(LUPRI,'(/A,I5,A,I10)')
     &         '(DALTON_NODEDRIVER) ERROR on slave',MYTID,
     &         ': Unknown calculation type request',IPRTYP
            WRITE(*,'(/A,I5,A,I10)')
     &         '(DALTON_NODEDRIVER) ERROR on slave',MYTID,
     &         ': Unknown calculation type request',IPRTYP
            CALL QUIT('DALTON_NODEDRIVER: Unknown calculation type')
         END IF
C
         GO TO 100
      END IF ! IPRTYP .gt. 0

!     release memory
      if(allocated(work)) deallocate(work)
C
      IF (MYNUM .LE. 20) THEN
         CALL GETTIM(CEND,WEND)
         CTOT = CEND - CSTR
         WTOT = WEND - WSTR
C
         WRITE (LUPRI,'(/A,I10/)') 'Number of tasks performed:',N_TASKS
         CALL TIMTXT(' Total CPU  time used:',CTOT,LUPRI)
         CALL TIMTXT(' Total wall time used:',WTOT,LUPRI)
C
C        Stamp date and time and hostname to output
C
         CALL TSTAMP(' ',LUPRI)
      END IF
C
      CALL QEXIT('DALTON_NODEDRIVER')
      RETURN
      END
#endif
C  /* Deck sendname */
      SUBROUTINE SENDNAME(NASTER,NYTID)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mtags.h"
C
#include "infpar.h"
C
      MASTER = NASTER
      MYTID  = NYTID
C
#if defined (VAR_PVM)
C
C     NCODE is the first variable PVM-slaves must receive
C
      CALL PVMXRECV(MASTER,MTAG0)
      CALL PVMXUNPACK(NCODE,1,'INTEGER')
#endif
C     Get my name
C
      CALL GETNAME(MYNAME)
C
C     Return name to master
C
#if defined (VAR_MPI)
      CALL MPIXSEND(MYNAME,20,'STRING',MASTER,MTAG8)
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(MYTID,1,'INTEGER')
      CALL PVMXPACK(MYNAME,20,'STRING')
      CALL PVMXSEND(NASTER,MTAG8)
#endif
C
      RETURN
      END
C  /* Deck recvnames */
      SUBROUTINE RECVNAMES
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mtags.h"
      CHARACTER NAME*20
#include "infpar.h"
C
#if defined (VAR_PVM)
C
C     NCODE must be the first variable PVM-slaves must receive
C
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(NCODE,1,'INTEGER')
      CALL PVMXMCAST(NODTOT,NODEID(1),MTAG0)
#endif
C     Get my name
C
      CALL GETNAME(MYNAME)
      NODNAM(0) = MYNAME
C
C     Return name to master
C
      DO I = 1,NODTOT
#if defined (VAR_MPI)
         NWHO = -1
         CALL MPIXRECV(NAME,20 ,'STRING',NWHO,MTAG8)
         NODNAM(NWHO) = NAME
#endif
#if defined (VAR_PVM)
         CALL PVMXRECV(-1,MTAG8)
         CALL PVMXUNPACK(NTID,1,'INTEGER')
         CALL PVMXUNPACK(NAME,20,'STRING')
         DO J=1, NODTOT
            IF (NODEID(J) .EQ. NTID) NWHO = J
         END DO
         NODNAM(NWHO) = NAME
#endif
      END DO
C
      IF (IPRPAR.GT.0) THEN
         WRITE (LUPRI,'(3X,A,A)') ' Slaves are : ',NODNAM(1)
         WRITE (LUPRI,'(17X,A)') (NODNAM(I),I=2,NODTOT)
      END IF
C
      RETURN
      END
C  /* Deck eri_nodstr */
      SUBROUTINE ERI_NODSTR(WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "maxorb.h"
#include "aovec.h"
C
      PARAMETER (MXDMAT = 100)
      DIMENSION WORK(LWORK), ISYMDM(MXDMAT), IFCTYP(MXDMAT)
C
C Used from common blocks
C
C MXORB  : MXSHEL (for INFPAR)
C INFPAR : SLAVE
C INFORB : N2BASX
C
#include "infpar.h"
#include "inforb.h"
#include "veclen.h"
C
#include "memint.h"
C
      CALL QENTER('ERI_NODSTR')
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_NODSTR','*',103)
C
C     Start timing
C
      CALL GETTIM(CPU1,WALL1)
C
C     Tell integral section this is a slave in a parallel run.
C
      SLAVE = .TRUE.
C
C     Receive initialization from Master
C
      NPAO = MXSHEL*MXAOVC
      CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)
      CALL ERI_RVINIT(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &                WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),
     &                ISYMDM,IFCTYP,WORK,LWORK,IPRINT)
C
C     Allocate workspace
C
      NDMAT  = 1
#if defined (VAR_VECTOR)
      ICHUNK = MAX(IVECLN/NDMAT,1)
#endif
C
      KDMAT  = 1
      KFMAT  = KDMAT  + N2BASX
      KCCFBT = KFMAT  + N2BASX
      KINDXB = KCCFBT + MXPRIM*MXCONT
      KFCKPT = KINDXB + (MXSHEL*MXCONT*8 - 1)/IRAT + 1
#if defined (VAR_VECTOR)
      KLAST  = KFCKPT + ICHUNK*NDMAT*(NBASE + NODD)*NBASE
#else
      KLAST  = KFCKPT
#endif
      IF (KLAST .GT. LWORK) CALL STOPIT('ERI_NODSTR',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST
C
C     Zero Fock-matrix
C
      CALL DZERO(WORK(KFMAT),N2BASX)
C
      CALL ERIFCK(WORK(KFMAT),WORK(KDMAT),NDMAT,ISYMDM,IFCTYP,IPRINT,
     &            WORK(KCCFBT),WORK(KINDXB),WORK(KFCKPT),
     &            WORK(KLAST),LWRK)
C
C     Send overall timings if requested.
C
      IF (IPRPAR .GT. 0) CALL SDTIM(CPU1,WALL1,IPRINT)
C
      CALL QEXIT('ERI_NODSTR')
      RETURN
      END
C  /* Deck eri_rvinit */
      SUBROUTINE ERI_RVINIT(JSTRSH,NPRIMS,NCONTS,IORBSH,JORBSH,KORBSH,
     &                      ISYMDM,IFCTYP,WORK,LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "mtags.h"
C
      LOGICAL FIRST, TSTINP
      DIMENSION JSTRSH(*), NPRIMS(*), NCONTS(*), IORBSH(*),
     &          JORBSH(*), KORBSH(*), ISYMDM(*), IFCTYP(*),
     &          WORK(LWORK)
C
C Used from common blocks
C
C MXCENT : MXCENT (for MOLINP)
C MXORB  : MXSHEL (for INFPAR)
C GNRINF : NEWBAS, (PVM: BASDIR, BASLEN), RELCAL
C SIRIPC : NEWGEO
C MOLINP : NMLINE
C MOLINC : MLINE()
C INFORB : *
C INFPAR : *
C CBIERI : *
C
#include "gnrinf.h"
#include "siripc.h"
#include "molinp.h"
#include "infpar.h"
#include "inforb.h"
#include "r12int.h"
#include "cbieri.h"
#include "cbirea.h"
#if defined (SYS_T3D)
C
C     This kludge fixes the current non-implementation of MPI_CHARACTER
C     on the Cray-T3E, K.Ruud-Nov.96
C     This Kludge fix is soon old enough to start at University...
      DIMENSION ITMP(4000)
      EQUIVALENCE (ITMP,MLINE)
#endif
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_RVINIT','*',103)
C
#if defined (VAR_MPI)
C
C     Master tells nodes whether to read new input information.
C
      CALL MPIXBCAST(FIRST ,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWGEO,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWBAS,1,'LOGICAL',MASTER)
C
      TIMING = FIRST .OR. NEWBAS
C
C     Receive /MOLINC/,/MOLINP/, /INFORB/, /INFPAR/ and /CBIERI/
C
      IF (TIMING .OR. NEWGEO) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &        ' Receiving new orbital inf.',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C   /molinp/
         CALL MPIXBCAST(NMLINE,1,'INTEGER',MASTER)
C   /molinc/
#if defined (SYS_T3D)
         CALL MPIXBCAST(ITMP,10*NMLINE,'INTEGER',MASTER)
#else
         CALL MPIXBCAST(MLINE,len_MLINE*NMLINE,'STRING',MASTER)
#endif
C   /inforb/
         CALL MPIXBCAST(MULD2H,NINFI,'INTEGER',MASTER)
C   /infpar/
         CALL MPIXBCAST(NODTOT,NPARI,'INTEGER',MASTER)
C   /cbieri/
         CALL MPIXBCAST(NDMAT ,NCBII,'INTEGER',MASTER)
         CALL MPIXBCAST(RUNERI,NCBIL,'LOGICAL',MASTER)
C   /cbirea/
         CALL MPIXBCAST(LCMMAX,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NCMSTR,1,'INTEGER',MASTER)
         CALL MPIXBCAST(NCMEND,1,'INTEGER',MASTER)
         CALL MPIXBCAST(LCNTNUUM,1,'LOGICAL',MASTER)
C  comr12
         CALL MPIXBCAST(R12INT,NR12L,'LOGICAL',MASTER)
         CALL MPIXBCAST(LMULBS,1,'LOGICAL',MASTER)
         CALL MPIXBCAST(INTGAC,NR12I,'INTEGER',MASTER)
      END IF
      CALL MPIXBCAST(ISYMDM,NDMAT,'INTEGER',MASTER)
      CALL MPIXBCAST(IFCTYP,NDMAT,'INTEGER',MASTER)
#elif defined (VAR_PVM)
C
C     Receive initialization.
C
      CALL PVMXRECV(MASTER,MTAG3)
      CALL PVMXUNPACK(FIRST ,1,'LOGICAL')
      CALL PVMXUNPACK(NEWBAS,1,'LOGICAL')
      CALL PVMXUNPACK(NEWGEO,1,'LOGICAL')
      CALL PVMXUNPACK(BASDIR,60,'STRING')
      CALL PVMXUNPACK(LENBAS,1,'INTEGER')
C
      TIMING = FIRST .OR. NEWBAS
C
C     Receive /MOLINC/,/MOLINP/, /INFORB/, /INFPAR/ and /CBIERI/
C
      IF (TIMING .OR. NEWGEO) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &        ' Receiving new orbital inf.',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C   /molinp/
         CALL PVMXUNPACK(NMLINE,1,'INTEGER')
C   /molinc/
         CALL PVMXUNPACK(MLINE,len_MLINE*NMLINE,'STRING')
C   /inforb/
         CALL PVMXUNPACK(MULD2H,NINFI,'INTEGER')
C   /infpar/
         CALL PVMXUNPACK(NODTOT,NPARI,'INTEGER')
         DO 100, I=1, NODTOT
            IF (NODEID(I) .EQ. MYTID) MYNUM = I
 100     CONTINUE
C   /cbieri/
         CALL PVMXUNPACK(NDMAT ,NCBII,'INTEGER')
         CALL PVMXUNPACK(RUNERI,NCBIL,'LOGICAL')
C  comr12
         CALL PVMXBCAST(R12INT,NR12L,'LOGICAL')
         CALL PVMXBCAST(INTGAC,NR12I,'INTEGER')
      END IF
      CALL PVMXUNPACK(ISYMDM,NDMAT,'INTEGER')
      CALL PVMXUNPACK(IFCTYP,NDMAT,'INTEGER')
#else
!     to remove compiler warnings for sequential builds
      FIRST = .TRUE.
#endif
C
C     Set hermit
C
      IF (NEWGEO) RDINPC = .FALSE.
C
C     Since REAINI is not initialized outside of SETHER, make sure slaves
C     initialize it properly
C
      IF (FIRST) THEN
          RELCAL = .FALSE.
          TSTINP = .FALSE.
         CALL REAINI(0,RELCAL,TSTINP)
      END IF
C
      CALL SETHER(0,NEWGEO,WORK,LWORK)
C
C     Set common-block BLOCKS
C
      CALL PAOVEC(JSTRSH,NPRIMS,NCONTS,IORBSH,JORBSH,KORBSH,0,.FALSE.,0)
C
C     Get DMAT
C
      CALL DZERO(WORK,N2BASX)
C
#if defined (VAR_MPI)
      CALL MPIXBCAST(WORK,N2BASX,'DOUBLE',MASTER)
#endif
#if defined (VAR_PVM)
      CALL PVMXUNPACK(WORK,N2BASX,'DOUBLE')
#endif
C
      RETURN
      END
C  /* Deck sdtim */
      SUBROUTINE SDTIM(CPU1,WALL1,IPRINT)
C
C Used from common blocks:
C
C MXORB  : MXSHEL (for INFPAR)
C INFPAR : MASTER, MYNUM, NCODE
C
#include "implicit.h"
#include "maxorb.h"
#include "mtags.h"
      DIMENSION TMGARR(2)
#include "infpar.h"
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from SDTIM','*',103)
C
      CALL GETTIM(CPU2,WALL2)
C
      TMGARR(1) = CPU2  - CPU1
      TMGARR(2) = WALL2 - WALL1
C
#if defined (VAR_MPI)
      CALL MPIXSEND(TMGARR,2 ,'DOUBLE',MASTER,MTAG8)
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(MYNUM ,1,'INTEGER')
      CALL PVMXPACK(TMGARR,2,'DOUBLE')
      CALL PVMXSEND(MASTER,MTAG8)
#endif
C
      RETURN
      END
C  /* Deck eri_noddrv */
      SUBROUTINE ERI_NODDRV(IODCL1,IODCL2,
     &                      IODBC1,IODBC2,RODBC1,RODBC2,
     &                      IODPP1,IODPP2,RODPP1,RODPP2,
     &                      FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,
     &                      NACTAB,NACTCD,SCRNAB,SCRNCD,
     &                      INDHER,INDHSQ,IODDHR,
     &                      LMNPWR,LMNSYM,IPNTUV,
     &                      BUF,IBUF,
     &                      CLSCPU,IWHICH,IJS,MXCLS2,
     &                      CCFBT,INDXBT,NCOUNT,CSQ,WORK,LWORK,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
#include "mxcent.h"
#include "aovec.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "eribuf.h"
#include "dummy.h"
#include "mtags.h"
#include "cbieri.h"
#include "erithr.h"
#include "odclss.h"
#include "ccom.h"
#include "aobtch.h"
#include "doxyz.h"
#include "ericom.h"
#include "erimem.h"
#include "erista.h"
#include "hertop.h"
#include "symmet.h"
C
C Used from common blocks
C
C INFPAR: MYNUM, NTASK, TIMING
C
#include "infpar.h"
      LOGICAL DONE
      DIMENSION IODCL1(NODCL1,NITCL), IODCL2(NODCL2,NITCL),
     &          IODBC1(NODBC1,NITBC), IODBC2(NODBC2,NITBC),
     &          RODBC1(NODBC1,NRTBC), RODBC2(NODBC2,NRTBC),
     &          IODPP1(NODPP1,NITPP), IODPP2(NODPP2,NITPP),
     &          RODPP1(NODPP1,NRTPP), RODPP2(NODPP2,NRTPP),
     &          NACTAB(NODCL1), NACTCD(NODCL2),
     &          SCRNAB(NODBC1), SCRNCD(NODBC2),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP),
     &          INDHSQ(NRTOP), IODDHR(NRTOP,0:7),
     &          LMNPWR(KCKMAX,NHKMAX,3), LMNSYM(KCKMAX,NHKMAX,2,0:7),
     &          IPNTUV(KC2MAX,2,2),
     &          FMAT(*), DMAT(*), IFCTYP(NDMT), D2MAT(*), ID2MAT(*),
     &          BUF(LBUF,NBUFS), IBUF(LBUF,NIBUF,NBUFS),
     &          CCFBT(*), INDXBT(*),
     &          CLSCPU(MXCLS2), IWHICH(MXCLS2), IJS(NTASK),
     &          NCOUNT(*), CSQ(NCSQ1,NCSQ2), WORK(LWORK)

C
      CALL QENTER('ERI_NODDRV')
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_NODDRV','*',103)
C
C     Send NODCL1 and NTASK to master.
C
      IF (MYNUM .EQ. 1) THEN
         CALL ERI_SDCLS(NODCL1,IPRINT)
      END IF
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('Output from subroutine ODCDR1',-1)
         WRITE (LUPRI,'(2X,A,2I5)') ' NODCL1, NODCL2 ', NODCL1, NODCL2
         WRITE (LUPRI,'(2X,A,2I5)') ' NODPP1, NODPP2 ', NODPP1, NODPP2
         WRITE (LUPRI,'(2X,A,2L5)') ' ODTRI1, ODTRI2 ', ODTRI1, ODTRI2
         IF (IPRINT .GT. 8) THEN
            CALL HEADER('IODCL1 in ODCDR1',1)
            DO 100 I = 1, NODCL1
               WRITE (LUPRI,'(5X,I5,5X,14I5)') I,(IODCL1(I,J),J=1,NITCL)
  100       CONTINUE
            CALL HEADER('IODCL2 in ODCDR1',1)
            DO 110 I = 1, NODCL2
               WRITE (LUPRI,'(5X,I5,5X,14I5)') I,(IODCL2(I,J),J=1,NITCL)
  110       CONTINUE
            CALL HEADER('IODBC1 in ODCDR1',1)
            DO 120 I = 1, NODBC1
               WRITE (LUPRI,'(5X,I5,5X,7I5)') I,(IODBC1(I,J),J=1,NITBC)
  120       CONTINUE
            CALL HEADER('IODBC2 in ODCDR1',1)
            DO 130 I = 1, NODBC2
               WRITE (LUPRI,'(5X,I5,5X,7I5)') I,(IODBC2(I,J),J=1,NITBC)
  130       CONTINUE
            CALL HEADER('RODBC1 in ODCDR1',-1)
            CALL OUTPUT(RODBC1,1,NODBC1,1,NRTBC,NODBC1,NRTBC,1,LUPRI)
            CALL HEADER('RODBC2 in ODCDR1',-1)
            CALL OUTPUT(RODBC2,1,NODBC2,1,NRTBC,NODBC2,NRTBC,1,LUPRI)
            CALL HEADER('RODPP1 in ODCDR1',-1)
            CALL OUTPUT(RODPP1,1,NODPP1,1,NRTPP,NODPP1,NRTPP,1,LUPRI)
            CALL HEADER('RODPP2 in ODCDR1',-1)
            CALL OUTPUT(RODPP2,1,NODPP2,1,NRTPP,NODPP2,NRTPP,1,LUPRI)
         END IF
      END IF
C
C     COMMON /DOXYZ/
C
      DOX = .TRUE.
      DOY = .TRUE.
      DOZ = .TRUE.
C
      FIRST = .TRUE.
      LAST  = .FALSE.
C
      NPPA  = 0
      NCPA  = 0
      NPCA  = 0
      NCCA  = 0
      NODPQ = 0
C
      NTPAS = 0
      NPPXA = 0
      NPCXA = 0
      NCPXA = 0
      NCCXA = 0
      NPPSA = 0
      NPCSA = 0
      NCPSA = 0
      NCCSA = 0
C
      NPPX1 = 1000000
      NPCX1 = 1000000
      NCPX1 = 1000000
      NCCX1 = 1000000
      NPPS1 = 1000000
      NPCS1 = 1000000
      NCPS1 = 1000000
      NCCS1 = 1000000
C
      NPPX2 = 0
      NPCX2 = 0
      NCPX2 = 0
      NCCX2 = 0
      NPPS2 = 0
      NPCS2 = 0
      NCPS2 = 0
      NCCS2 = 0
C
      IPASSA = 0
      IPASS1 = 1000
      IPASS2 = 0
      IPASSN = 0
C
C     Arrays for Hermite integrals
C     ============================
C
      CALL HERPRP(INDHER,INDHSQ,IODDHR,IPRINT)
      CALL LMNPRP(LMNPWR,LMNSYM,IPRINT)
C
C     Active/inactive batches
C     =======================
C
      IELCTR = 1
      CALL ERIACT(NACTAB,SCRNAB,RODBC1,IODBC1,NODBC1,IODCL1,
     &            NODCL1,ODTRI1,IELCTR,IPRINT)
C
      IELCTR = 2
      CALL ERIACT(NACTCD,SCRNCD,RODBC2,IODBC2,NODBC2,IODCL2,
     &            NODCL2,ODTRI2,IELCTR,IPRINT)
C
      IPRKEP = IPRINT
C
      IJSTOT = 0
      CALL IZERO(IWHICH,MXCLS2)
      CALL DZERO(CLSCPU,MXCLS2)
C
  200 CONTINUE
C
C     Send ID number and ask for something to do.
C
#if defined (VAR_MPI)
      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MTAG5)
C
C     Receive IJ labels for NTASK new (class|class) integrals
C
      CALL IZERO(IJS,NTASK)
C
      CALL MPIXRECV(DONE,1,'LOGICAL',MASTER,MTAG6)
      IF (DONE) GOTO 500
      CALL MPIXRECV(IJS,NTASK,'INTEGER',MASTER,MTAG6)
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(MYNUM,1,'INTEGER')
      CALL PVMXSEND(MASTER,MTAG5)
C
C     Receive IJ labels for NTASK new (class|class) integrals
C
      CALL IZERO(IJS,NTASK)
C
      CALL PVMXRECV(MASTER,MTAG6)
      CALL PVMXUNPACK(DONE,1,'LOGICAL')
      IF (DONE) GOTO 500
      CALL PVMXUNPACK(IJS,NTASK,'INTEGER')
#endif
      IF (IPRINT .GT. 4) THEN
         WRITE(LUPRI,'(5X,A,I5,A)')
     &      'Receiving batch',IJSTOT,'from master'
      END IF
C
      DO 300, K = 1,NTASK
C
      IJCLS = IJS(K)
C
      IF (IJCLS .EQ. 0) THEN
         IJSTOT = IJSTOT + K - 1
         GOTO 200
      ELSE
         CALL UNPKIJ(IJCLS,IODAB,IODCD)
      END IF
C
      IF (TIMING) CPU1IJ = SECOND()
C
C     Here starts the integral calculation
C
      IF (NACTAB(IODAB) .GT. 0) THEN
         NODSAB = NACTAB(IODAB)
C
         NODTAB = IODCL1(IODAB, 1)
         KODSAB = IODCL1(IODAB, 2)
         SCRMAB = SCRNAB(KODSAB+IDAMAX(NODTAB,SCRNAB(KODSAB),1)-1)
C
         NHKTA  = IODCL1(IODAB, 3)
         NHKTB  = IODCL1(IODAB, 4)
         KHKTA  = IODCL1(IODAB, 5)
         KHKTB  = IODCL1(IODAB, 6)
         KCKTA  = IODCL1(IODAB, 7)
         KCKTB  = IODCL1(IODAB, 8)
         ISTBLA = IODCL1(IODAB, 9)
         ISTBLB = IODCL1(IODAB,10)
         NPRFAB = IODCL1(IODAB,11)
         NCTFAB = IODCL1(IODAB,12)
         NCNTAB = IODCL1(IODAB,13)
         NGTOAB = IODCL1(IODAB,14)
C
         MAXCD = NODCL2
         IF (ODTR12) MAXCD = IODAB
C
         IF (NACTCD(IODCD) .GT. 0) THEN
            NODSCD = NACTCD(IODCD)
C
            NODTCD = IODCL2(IODCD, 1)
            KODSCD = IODCL2(IODCD, 2)
            SCRMCD = SCRNCD(KODSCD+IDAMAX(NODTCD,SCRNCD(KODSCD),1)-1)
C
         IF (SCRMAB*SCRMCD .GT. THRSH) THEN
C
            IF (EXTPRI .AND. IODAB.EQ.IPROD1 .AND. IODCD.EQ.IPROD2) THEN
               IPRINT = 1000
            ELSE
               IPRINT = IPRKEP
            END IF
C
            NHKTC  = IODCL2(IODCD, 3)
            NHKTD  = IODCL2(IODCD, 4)
            KHKTC  = IODCL2(IODCD, 5)
            KHKTD  = IODCL2(IODCD, 6)
            KCKTC  = IODCL2(IODCD, 7)
            KCKTD  = IODCL2(IODCD, 8)
            ISTBLC = IODCL2(IODCD, 9)
            ISTBLD = IODCL2(IODCD,10)
            NPRFCD = IODCL2(IODCD,11)
            NCTFCD = IODCL2(IODCD,12)
            NCNTCD = IODCL2(IODCD,13)
            NGTOCD = IODCL2(IODCD,14)
C
            NODSPQ = NODSAB*NODSCD
            NPRFPQ = NPRFAB*NPRFCD
            NCTFPQ = NCTFAB*NCTFCD
C
C           Statistics
C
            NODPQ  = NODPQ + 1
            NPPI   = NODSPQ*NPRFPQ
            NCPI   = NODSPQ*NCTFAB*NPRFCD
            NPCI   = NODSPQ*NPRFAB*NCTFCD
            NCCI   = NODSPQ*NCTFAB*NCTFCD
            NPPA   = NPPA  + NPPI
            NCPA   = NCPA  + NCPI
            NPCA   = NPCA  + NPCI
            NCCA   = NCCA  + NCCI
C
            JMAXA  = NHKTA - 1
            JMAXB  = NHKTB - 1
            JMAXC  = NHKTC - 1
            JMAXD  = NHKTD - 1
            JMAXAB = JMAXA + JMAXB + MAXDER
            JMAXCD = JMAXC + JMAXD + MAXDER
            JMAX   = JMAXA + JMAXB + JMAXC + JMAXD + MAXDER
C
            NTUVAB = (JMAXAB + 1)*(JMAXAB + 2)*(JMAXAB + 3)/6
            NTUVCD = (JMAXCD + 1)*(JMAXCD + 2)*(JMAXCD + 3)/6
            NTUV   = (JMAX   + 1)*(JMAX   + 2)*(JMAX   + 3)/6
C
            SPHRA  = KHKTA .NE. KCKTA
            SPHRB  = KHKTB .NE. KCKTB
            SPHRC  = KHKTC .NE. KCKTC
            SPHRD  = KHKTD .NE. KCKTD
C
            SPHRAB = SPHRA .OR. SPHRB
            SPHRCD = SPHRC .OR. SPHRD
C
            ISTBLR = IOR(ISTBLA,ISTBLB)
            ISTBLS = IOR(ISTBLC,ISTBLD)
            ISTBLT = IOR(IAND(ISTBLA,ISTBLB),IAND(ISTBLC,ISTBLD))
C
            MLTPA = MULT(ISTBLA)
            MLTPB = MULT(ISTBLB)
            MLTPC = MULT(ISTBLC)
            MLTPD = MULT(ISTBLD)
C
C
            MLTPR = MULT(ISTBLR)
            MLTPS = MULT(ISTBLS)
            MLTPT = MULT(ISTBLT)
            MLTPX = MLTPR*MLTPS*MLTPT
C
            TKMPAB = NCNTAB.EQ.1 .AND. KHKTA.EQ.KHKTB .AND..NOT.DODIST
            TKMPCD = NCNTCD.EQ.1 .AND. KHKTC.EQ.KHKTD .AND..NOT.DODIST
            TCMPAB = TKMPAB .AND. .NOT.SPHRAB
            TCMPCD = TKMPCD .AND. .NOT.SPHRCD
C
            DIAGAB = NGTOAB .EQ. 1 .AND. ODTRI1
            DIAGCD = NGTOCD .EQ. 1 .AND. ODTRI2
C
            KHKTAB = KHKTA*KHKTB
            KHKTCD = KHKTC*KHKTD
            KCKTAB = KCKTA*KCKTB
            KCKTCD = KCKTC*KCKTD
            IF (TKMPAB) KHKTAB = KHKTA*(KHKTA + 1)/2
            IF (TKMPCD) KHKTCD = KHKTC*(KHKTC + 1)/2
            IF (TCMPAB) KCKTAB = KCKTA*(KCKTA + 1)/2
            IF (TCMPCD) KCKTCD = KCKTC*(KCKTC + 1)/2
C
C           Memory requirements
C
            CALL MEMPQB(MEMBCH,NALLOI,IPRINT)
            IF (MEMBCH .GT. LWORK) THEN
               IF (IPRINT .GT. 5) THEN
                  WRITE (LUPRI,'(/,1X,A,2I5,2(/1X,A,I10))')
     &              ' Not enough memory for one PQ batch for'//
     &              ' OD classes: ',IODAB, IODCD,
     &              ' Available memory in ODCPAR:',LWORK,
     &              ' Required memory for ODCPAR:',MEMBCH
               END IF
               MEMOK = .FALSE.
               NEEDMR = MEMBCH - LWORK
               IF (NEEDMR .GT. MEMADD) THEN
                  MEMADD = NEEDMR
                  MODAB  = IODAB
                  MODCD  = IODCD
               END IF
            ELSE IF (MEMOK) THEN
               IJ = 0
               DO 400 I = KODSAB, KODSAB + NODTAB - 1
                  MAXJ = KODSCD + NODTCD - 1
                  IF (ODTR12 .AND. IODAB.EQ.IODCD) MAXJ = I
                  DO 410 J = KODSCD, MAXJ
                     IF (SCRNAB(I)*SCRNCD(J).GT.THRSH) IJ = IJ + 1
  410             CONTINUE
  400          CONTINUE
               MAXBCH = MIN(IJ,(IRAT*LWORK-(IRAT-1)*NALLOI)/MEMBCH)
            END IF
C
C           Print level zero for parallel calculations
C
            IF (IPRINT .GT. 10) THEN
               CALL TITLER('New pair of OD classes','*',103)
            END IF
            IF (IPRINT .GT. 5) THEN
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' New OD classes:        ',IODAB,IODCD
            END IF
            IF (IPRINT .GT. 10) THEN
               WRITE (LUPRI,'(1X,A,4I5)')
     &            ' Angular momenta (NHKT):',NHKTA,NHKTB,NHKTC,NHKTD
               WRITE (LUPRI,'(1X,A,3I5)')
     &            ' Number of active OD batches:  ',NODSAB,NODSCD,NODSPQ
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' NODTAB, NODTCD:        ',NODTAB,NODTCD
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' KODSAB, KODSCD:        ',KODSAB,KODSCD
               WRITE (LUPRI,'(1X,A,3I5)')
     &            ' NPRFAB, NPRFCD:        ',NPRFAB,NPRFCD,NPRFPQ
               WRITE (LUPRI,'(1X,A,3I5)')
     &            ' NCTFAB, NCTFCD:        ',NCTFAB,NCTFCD,NCTFPQ
               WRITE (LUPRI,'(1X,A,4I5)')
     &            ' NPP,NCPI,NPCI,NCCI:    ',NPPI,NCPI,NPCI,NCCI
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' NCNTAB, NCNTCD:        ',NCNTAB,NCNTCD
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' NGTOAB, NGTOCD:        ',NGTOAB,NGTOCD
               WRITE (LUPRI,'(1X,A,2L5)')
     &            ' TKMPAB, TKMPCD:        ',TKMPAB,TKMPCD
               WRITE (LUPRI,'(1X,A,2I5)')
     &            ' MEMBCH, MAXBCH:        ',MEMBCH,MAXBCH
            END IF
C
            IF (MEMOK) THEN
C
C              Allocate ODCPAR
C
               NITPQ  = 6
               KPNTCR = 1
               KPNTPP = KPNTCR + (4*MAXBCH - 1)/IRAT + 1
               KLAST  = KPNTPP + (2*MAXBCH - 1)/IRAT + 1
               IF (KLAST.GT.LWORK) CALL STOPIT('ODCDR1',' ',KLAST,LWORK)
               LWRK  = LWORK - KLAST + 1
               CALL ODCPAR(IODBC1,IODBC2,
     &                     IODPP1,IODPP2,RODPP1,RODPP2,
     &                     SCRNAB,SCRNCD,
     &                     INDHER,INDHSQ,IODDHR,
     &                     LMNPWR,LMNSYM,IPNTUV,
     &                     WORK(KPNTCR),WORK(KPNTPP),
     &                     FMAT,DMAT,NDMT,IFCTYP,D2MAT,ID2MAT,BUF,IBUF,
     &                     CCFBT,INDXBT,NCOUNT,CSQ,WORK(KLAST),
     &                     LWRK,IPRINT)
            END IF
         END IF
         END IF
      END IF
C
      IF (TIMING) THEN
         CPU2IJ = SECOND()
         CPUIJ  = CPU2IJ - CPU1IJ
         IPLACE = IJSTOT + K
C
         IWHICH(IPLACE) = IJCLS
         CLSCPU(IPLACE) = CPUIJ
      END IF
C
 300  CONTINUE
C
      IJSTOT = IJSTOT + NTASK
      GOTO 200
C
 500  CONTINUE
C
      CALL ERI_SDRES(DMAT,FMAT,NBASE,CLSCPU,IWHICH,IJSTOT,LWORK,
     &               IPRINT)
C
      IPRINT = IPRKEP
C
      CALL QEXIT('ERI_NODDRV')
      RETURN
      END
C  /* Deck eri_sdcls */
      SUBROUTINE ERI_SDCLS(NODCLS,IPRINT)
C
#include "implicit.h"
#include "maxorb.h"
#include "mtags.h"
C
C Used from common blocks
C
C  MXORB : MXSHEL (for INFPAR)
C  INFPAR: MASTER
C
#include "infpar.h"
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_SDCLS','*',103)
C
#if defined (VAR_MPI)
      CALL MPIXSEND(NODCLS,1,'INTEGER',MASTER,MTAG4)
      CALL MPIXSEND(NTASK ,1,'INTEGER',MASTER,MTAG4)
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(NODCLS,1,'INTEGER')
      CALL PVMXPACK(NTASK ,1,'INTEGER')
      CALL PVMXSEND(MASTER,MTAG4)
#endif
C
      RETURN
      END
C  /* Deck sdres */
      SUBROUTINE ERI_SDRES(DMAT,FMAT,NBASE,CLSCPU,IWHICH,IJSTOT,
     &                     LWORK,IPRINT)
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mtags.h"
C
      DIMENSION DMAT(*), FMAT(*), CLSCPU(IJSTOT), IWHICH(IJSTOT)
C
C Used from common blocks
C
C  MXORB  : MXSHEL (for INFPAR)
C  INFPAR : NCODE, MASTER, MYNUM, TIMING
C  ERIMEM : *
C  CBIERI : NDMAT
C
#include "infpar.h"
#include "erimem.h"
#include "cbieri.h"
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from ERI_SDRES','*',103)
C
      IF (MEMOK) THEN
         IF (IPRINT.GT.5) THEN
            CALL HEADER('Density and Fock matrices in ERI_SDRES',-1)
            KSTR = 1
            DO 300 I = 1, NDMAT
               WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
               CALL OUTPUT(DMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                     NBASE,1,LUPRI)
               WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
               CALL OUTPUT(FMAT(KSTR),1,NBASE,1,NBASE,NBASE,
     &                     NBASE,1,LUPRI)
               KSTR = NBASE*NBASE
  300       CONTINUE
            CALL FLSHFO(LUPRI)
         END IF
C
C        Send MEMOK and Fock matrix
C
#if defined (VAR_MPI)
         CALL MPIXSEND(MEMOK,1,'LOGICAL',MASTER,MTAG7)
         CALL MPIXSEND(FMAT,NBASE*NBASE,'DOUBLE',MASTER,MTAG7)
         CALL MPIXSEND(IJSTOT,1,'INTEGER',MASTER,MTAG7)
         IF (TIMING) THEN
            CALL MPIXSEND(IWHICH,IJSTOT,'INTEGER',MASTER,MTAG7)
            CALL MPIXSEND(CLSCPU,IJSTOT,'DOUBLE' ,MASTER,MTAG7)
         END IF
      ELSE
         CALL MPIXSEND(MEMOK ,1,'LOGICAL',MASTER,MTAG7)
         CALL MPIXSEND(MEMADD,3,'INTEGER',MASTER,MTAG7)
         CALL MPIXSEND(LWORK ,1,'INTEGER',MASTER,MTAG7)
#endif
#if defined (VAR_PVM)
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(MEMOK,1,'LOGICAL')
         CALL PVMXPACK(FMAT,NBASE*NBASE,'DOUBLE')
         CALL PVMXPACK(IJSTOT,1,'INTEGER')
         IF (TIMING) THEN
            CALL PVMXPACK(IWHICH,IJSTOT,'INTEGER')
            CALL PVMXPACK(CLSCPU,IJSTOT,'DOUBLE')
         END IF
         CALL PVMXSEND(MASTER,MTAG7)
      ELSE
         CALL PVMXINITSEND(NCODE)
         CALL PVMXPACK(MEMOK ,1,'LOGICAL')
         CALL PVMXPACK(MEMADD,3,'INTEGER')
         CALL PVMXPACK(LWORK ,1,'INTEGER')
         CALL PVMXSEND(MASTER,MTAG7)
#endif
      END IF
C
      RETURN
      END
