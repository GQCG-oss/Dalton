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
#ifdef UNDEF
===========================================================================
/* Comdeck notes */
950620-kr: Ported to MPI, based on previous efforts by Bjarne Herland at
           Para//ab
950531-kr: Reorganized the code so that the node programs are
started only once for a given calculation
(1)  Note: in old HERFCK 'IPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD' was
 transferred to TWOINT, do we want to do that in CALINT for ITYPE=3 ???
(2)  removed JSYMDM,JFCTYP,ITRIPL from /PARINT/ in aba.cdk
===========================================================================
#endif
C  /* Deck her_pardrv */
      SUBROUTINE HER_PARDRV(FMAT,DMAT,HESSEE,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                      WORK,LWORK,ITYPE,MAXDER,IATOM,NODV,
     &                      NOPV,NOCONT,TKTIME,RETUR,IPRINT,
     &                      ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,
     &                      DINTSKP,RELCAL)

C
#include <implicit.h>
#include <iratdef.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h>
C
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,RELCAL
      DIMENSION FMAT(*), DMAT(N2BASX,NDMAT), ISYMDM(*), IFCTYP(*),
     &          NSTAT(NODTOT), WORK(LWORK)
CTROND
      DIMENSION GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*)
CTROND
C
C Used from common blocks:
C
C AOVEC : MXAOVC (for BLOCKS)
C MXCENT: MXCOOR
C MXORB : MXSHEL (for INFPAR)
C INFPAR: NDEGDI,NTASK,MAXTSK
C BLOCKS: MAXSHL
C NUCLEI: NUCDEG
C INFORB: NNBASX,N2BASX
C
#include <infpar.h>
#include <blocks.h>
#include <nuclei.h>
#include <inforb.h>
#include <memint.h>
C
      CALL QENTER('HER_PARDRV')
      IF (IPRINT .GT. 3) CALL TITLER('Output from HER_PARDRV','*',103)
C
C     Determine total number of tasks - MTOTTK.
C
      MTOTTK = MAXSHL*(MAXSHL+1)/2
C
      IF (MTOTTK .GT. MAXTSK) THEN
         WRITE(LUPRI,'(//5X,A,/5X,A,/5X,2(A,I5))')
     &    'ERROR : Number of tasks greater than maximum.',
     &    'Recompilation will be necessary',
     &    'MTOTTK =',MTOTTK, 'MAXTSK =',MAXTSK
         CALL RECOMPI
         CALL QUIT('ERROR: MTOTTK greater than MAXTSK in HER_PARDRV')
      ENDIF
C
C     Determine dimension of array TMPMAT
C
C     Expectation values of differentiated integrals (all atoms)
      IF (ITYPE .EQ. 2) THEN
         NUMGRA = MXCOOR
         NUMHES = NUMGRA*NUMGRA
         NTMPMT = NUMGRA + NUMHES
C     Direct calculation of Fock matrices 
C        using skeleton matrix approach (AO basis)
      ELSE IF (ITYPE .EQ. 3) THEN
         NUMFCK = NDMAT*N2BASX
         NUMSKP = 8
#if defined(_CRAYT3E) && !defined(NO_BINSUM)
         NTMPMT = NUMSKP
#else
         NTMPMT = NUMFCK + NUMSKP
#endif
C     Derivatives with respect to magnetic field(DDFOCK=T)
      ELSE IF (ITYPE .EQ. -5) THEN
         NUMSUS = 9
         NTMPMT = 3*N2BASX + NUMSUS
      ELSE IF (ITYPE .EQ. 8) THEN
C     Derivative Fock matrix for specified atom
         NTMPMT = 3*NUCDEG(IATOM)*N2BASX
      END IF
C
C     Determine number of tasks pr. batch - NTASK.
C
      NTASK = MAX(1,MIN(MTOTTK,INT(NDEGDI*MTOTTK/(100*NODTOT))))
C
C     Determine dimension of array INDEX
C
      IBLOCK = NODTOT*NTASK
      NINDEX = (IBLOCK + 1)*INT(MTOTTK/IBLOCK + 1)
C
C     Be sure to allocate enough workspace for array SORTED
C
      IF (NTMPMT .LT. (MTOTTK + 2)) NTMPMT = MTOTTK + 2
C
      KTMPMT = 1
      KINDEX = KTMPMT +  NTMPMT
      KTIMES = KINDEX + (NINDEX + 1)/IRAT
      KTKCPU = KTIMES +  MTOTTK
      KWHICH = KTKCPU +  MTOTTK
      KLAST  = KWHICH + (MTOTTK + 1)/IRAT 
C
      IF (KLAST .GT. LWORK)
     &   CALL STOPIT('HER_PARDRV',' ',KLAST,LWORK)
      LWRK = LWORK - KLAST
C
      CALL HER_PARDR1(FMAT,DMAT,HESSEE,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                WORK(KTMPMT),WORK(KINDEX),WORK(KTIMES),
     &                WORK(KTKCPU),WORK(KWHICH),WORK(KLAST),LWRK,
     &                NINDEX,NTMPMT,NUMGRA,NUMHES,NUMSUS,NUMFCK,
     &                ITYPE,MAXDER,IATOM,NODV,NOPV,NOCONT,
     &                TKTIME,RETUR,IPRINT,
     &                ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,
     &                DINTSKP)
C
      CALL QEXIT('HER_PARDRV')
      RETURN
      END
C  /* Deck her_pardr1 */
      SUBROUTINE HER_PARDR1(FMAT,DMAT,HESSEE,NDMAT,ISYMDM,IFCTYP,NSTAT,
     &                      TMPMAT,INDEX,TIMES,TSKCPU,IWHICH,WORK,LWORK,
     &                      NINDEX,NTMPMT,NUMGRA,NUMHES,NUMSUS,NUMFCK,
     &                      ITYPE,MAXDIF,IATOM,NODV,NOPV,
     &                      NOCONT,TKTIME,RETUR,IPRINT,
     &                      ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO,
     &                      DINTSKP)
C
C
C     **************************************************************
C     *     Driver routine for the parallel HERMIT calculation.    *
C     **************************************************************
C
C
C      A short survey of the different messagetags (MTAGX) :
C
C
C      NEW LABELS:
C
C      10 - Tell nodes this is a HERMIT run, send NCODE and IPRINT
C      30 - Send initialization to nodes.
C      40   >Label used in ERI only<
C      50 - Receive (IJ|**) request from node
C      60 - Send new (IJ|**) labels to node
C      70 - Receive results from nodes
C      80 - Receive overall timings
C
C
C      OLD LABELS:
C
C      10 - Send initialization to nodes.
C      20 - Receive request from node.
C      30 - Send batch to node.
C      40 - Receive final results from nodes.
C      50 - Get timing results from nodes.
C
C           NCODE  - PVM-message encoding option.
C           ITYPE  - Calc. type: 2 = expectation values.
C                                3 = fock matrices.
C                               -5 = Integrals derivated with respect
C                                    to magnetic field.
C                                8 = Fock matrices derivated with respect
C                                    to atomic coordinates.
C           NTASK  - Number of tasks
C           IPRPAR - Print level during parallelization.
C           MAXDIF -
C           MAXREP -
C           ISYMDM - Symmetries of density matrices
C           IFCTYP - Fock matrix types (see twoint).
C           IATOM  -
C           NODV   - Density matrix neglected in TWOEXP
C           NOPV   -           -   "   -
C           NOCONT -           -   "   -
C           RETUR  - Program will exit after spec. shells
C           TKTIME - Take time in TWOINT
C           TIMING - Take time for each integralbatch IJ
C           DOREPS -
C           DOSYM  -
C           DOCOOR -
C           MULD2H - Start address for common block INFORB.
C           NMLINE - Number of lines in MOLECULE.INP
C           MLINE  - MOLECULE.INP as internal file
C           NDMAT  - Number of density-matrices
C           DMAT   - Density-matrices.
C           HFXFAC - HF exchange factor             
C
C
C
C NOTE: If the program get short of memory in calculation of
C       right hand sides (itype = 8) one can reduce the temporary
C       matrix TMPMAT used in HER_RVRES from 3*N2BASX to NNBASX by
C       looping over Fock matrices in HER_SDRES and HER_RVRES  
C       and then performing a send each time rather than transferring
C       all matrices in one batch.
C
C
C
#include <implicit.h>
#include <maxorb.h>
C
C Used from common blocks
C 
C MXORB : MXSHEL (used for INFPAR)
C GNRINF: NEWBAS
C SIRIPC: NEWGEO
C INFPAR: MAXTSK, NODTOT, NTASK, IPRPAR
C INFORB: N2BASX
C
#include <gnrinf.h>
#include <siripc.h>
#include <inforb.h>
#include <infpar.h>
C
      LOGICAL FIRST,NODV,NOPV,NOCONT,RETUR,TKTIME
      DIMENSION FMAT(*), DMAT(N2BASX,NDMAT), ISYMDM(*), IFCTYP(*),
     &          NSTAT(NODTOT), TMPMAT(NTMPMT), INDEX(NINDEX),
     &          TIMES(MTOTTK), TSKCPU(MTOTTK), IWHICH(MTOTTK),
     &          WORK(LWORK), NPOS(0:MAXTSK+1), HESSEE(*)
CTROND
      DIMENSION GABRAO(*),DMRAO(*),DMRSO(*),DINTSKP(*)
CTROND
C
C     NB! NPOS cannot be allocated from WORK due to the SAVEing below.
C
      SAVE FIRST,NPOS
      DATA FIRST/.TRUE./
C
      CALL QENTER('HER_PARDR1')
      IF (IPRINT .GT. 3) CALL TITLER('Output from HER_PARDR1','*',103)
C
C     Take calculation time of each IJ pair?
C
      TIMING = FIRST .OR. NEWBAS .OR. NEWGEO
C
C--------------------------
C     Make the index array.
C--------------------------
C
      CALL INDEKS(TIMING,INDEX,NINDEX,MTOTTK,NODTOT,NTASK,NBATCH,
     &            NPOS,IPRINT)
C
C---------------------------------------
C     Send initialization data to nodes.
C---------------------------------------
C
      CALL HER_SDINIT(DMAT,NDMAT,ISYMDM,IFCTYP,ITYPE,MAXDIF,IATOM,
     &                NODV,NOPV,NOCONT,TKTIME,RETUR,FIRST,IPRINT,
     &                ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO)
C
C-----------------------
C     Start calculation.
C-----------------------
C
      CALL DOPAR(INDEX,MTOTTK,NBATCH,NSTAT,IPRINT)
C
C--------------------------------------
C     Receive final results from nodes.
C--------------------------------------
C
      CALL HER_RVRES(FMAT,TMPMAT,HESSEE,TIMES,TSKCPU,IWHICH,NTMPMT,
     &               NUMGRA,NUMHES,NUMSUS,NUMFCK,ITYPE,IPRINT,DINTSKP)
C
C-----------------------------------------
C     Sort the integral-calculation times.
C-----------------------------------------
C
      IF (TIMING) CALL PARSRT(MTOTTK,NPOS,TMPMAT,TIMES,IPRINT)
C
C     NB! Array TMPMAT is only used as
C         workspace for array SORTED.
C
      IF (FIRST)  FIRST  = .FALSE.
      IF (NEWGEO) NEWGEO = .FALSE.
      IF (NEWBAS) NEWBAS = .FALSE.
C
      CALL QEXIT('HER_PARDR1')
      RETURN
      END
C  /* Deck indeks */
      SUBROUTINE INDEKS(TIMING,INDEX,NINDEX,MTOTTK,NODTOT,NTASK,NBATCH,
     &                  NPOS,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION INDEX(NINDEX), NPOS(0:MTOTTK+1)
      LOGICAL TIMING
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from INDEKS','*',103)
C
C     This subroutine determines in which order integrals will be calculated.
C
      CALL IZERO(INDEX,NINDEX)
C
      ICOUNT = 0
      IBLOCK = NODTOT*NTASK
C
      DO 100 I = MTOTTK, 1, -1
         MORE1  = (ICOUNT/IBLOCK)*IBLOCK
         MORE2  = MOD(ICOUNT,NODTOT)*NTASK
         MORE3  = ((ICOUNT-MORE1)/NODTOT) + 1
         IPLACE = MORE1+MORE2+MORE3
C
         IF (TIMING) THEN
            INDEX(IPLACE) = I
         ELSE
            INDEX(IPLACE) = NPOS(I)
        END IF
C
         ICOUNT=ICOUNT+1
 100  CONTINUE
C
C     Determine number of batches (NBATCH).
C
      NOTALL = (MTOTTK/IBLOCK)*IBLOCK
      NREST  = MTOTTK - NOTALL
      IF (NREST .GT. NODTOT) NREST = NODTOT
      NBATCH = NOTALL/NTASK + NREST
C
C     Print current settings.
C
      NEACH = NBATCH/NODTOT
      NREST = NBATCH - NEACH*NODTOT
C
      IF (IPRINT .GT. 4) THEN
         WRITE(LUPRI,'(5(/5X,A,I5),2X,A7,I3,A1/)')
     &        'Number of nodes   :',NODTOT,
     &        'Number of tasks   :',MTOTTK,
     &        'Tasks pr. batch   :',NTASK,
     &        'Number of batches :',NBATCH,
     &        'Batches pr. node  :',NEACH,
     &        '(rest =',NREST,')'
C
C        Print the index array.
C
         WRITE(LUPRI,'(/5X,A/)') 'The order of task-distribution :'
         DO 400 I=1, NBATCH*NTASK
            WRITE(LUPRI,'(5X,A3,I5,A15,I5)')
     &           'I =',I,'index(i) =',INDEX(I)
 400      CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck her_sdinit */
      SUBROUTINE HER_SDINIT(DMAT,NDMAT,ISYMDM,IFCTYP,ITYPE,MAXDIF,
     &                      IATOM,NODV,NOPV,NOCONT,TKTIME,RETUR,
     &                      FIRST,IPRINT,
     &                      ICEDIF,IFTHRS,GABRAO,DMRAO,DMRSO)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mtags.h>
      LOGICAL NODV,NOPV,NOCONT,TKTIME,RETUR,FIRST
      DIMENSION DMAT(N2BASX,NDMAT), ISYMDM(*), IFCTYP(*)
CTROND
CKR      DIMENSION GABRAO(*), DMRAO(*), DMRSO(*), DINTSKP(*)
      DIMENSION GABRAO(*), DMRAO(*), DMRSO(*)
CTROND
#include <infpar.h>
C
C Used from common blocks:
C
C MXCENT : MXCENT (for DORPS, etc.)
C MXORB  : MXSHEL (for INFPAR)
C MAXAQN : *      (for SYMMET)
C SIRIPC : NEWGEO (PVM: BASDIR, LENBAS)
C INFORB : N2BASX
C MOLINP : NMLINE, MLINE()
C DORPS  : DOREPS(), DOCOOR()
C ABAINF : DOSYM()
C SYMMET : MAXREP
C NUCLEI : NUCDEG()
C INFPAR : *
C HRUNIT : *
C PARINT : *
C
#include <gnrinf.h>
#include <siripc.h>
#include <inforb.h>
#include <molinp.h>
#include <dorps.h>
#include <abainf.h>
#include <symmet.h>
#include <nuclei.h>
#include <inftap.h>
#include <dftcom.h>
#include <parint.h>
#include <blocks.h>
#include <r12int.h>
#if defined (SYS_T3e)
C
C     This kludge fixes the current non-implementation of MPI_CHARACTER on
C     the Cray-T3E, K.Ruud-Nov.96
C
      DIMENSION ITMP(4000)
      EQUIVALENCE (ITMP,MLINE)
#endif
C
#include <cbisol.h>
C
      CALL QENTER('HER_SDINIT')
      IF (IPRINT .GT. 3) CALL TITLER('Output from HER_SDINIT','*',103)
C
C     Calculate the number of Fock matrices that will be
C     returned from slaves to master.
C
      IF (ITYPE .EQ. 2) THEN
         NFMAT = 0
      ELSE IF (ITYPE .EQ. 3) THEN
         NFMAT = NDMAT
      ELSE IF (ITYPE .EQ. -5) THEN
         NFMAT = 3*NDMAT
      ELSE IF (ITYPE .EQ. 8) THEN
         NFMAT = 3*NDMAT*NUCDEG(IATOM)
      ELSE
         WRITE (LUPRI,'(//1X,A,/1X,A,I2)')
     &    'ERROR: specified ITYPE for TWOINT not defined in HER_SDINIT',
     &    'ITYPE =',ITYPE
         CALL QUIT('Specified ITYPE not defined in HER_SDINIT.')
      END IF
C
C     Set common-block PARINT
C
      JATOM  = IATOM
      JLUDAS = LUDASP
      JLUINT = LUINTA
      JLUONE = LUONEL
      JLUSOL = LUSOL
      JLUSUP = LUSUPM
      JMXDIF = MAXDIF
      JMXREP = MAXREP
      JNFMAT = NFMAT
      JTASK  = NTASK
      JTYPE  = ITYPE
CTROND: Modify when merging with DIRAC...
CTROND      J2TYP  = I2TYP
      J2TYP  = 0
      JCEDIF = ICEDIF
      JFTHRS = IFTHRS
      JNODV  = NODV 
      JNOPV  = NOPV 
      JNOCNT = NOCONT
      JRETUR = RETUR
      JTKTIM = TKTIME
      JSOLVN = SOLVNT
      JRELCL = RELCAL
C
#if defined (VAR_MPI)
      CALL MPIXBCAST(FIRST ,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWBAS,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWGEO,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(HFXFAC,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(HFXMU, 1,'DOUBLE',MASTER)
C
      IF (FIRST .OR. NEWBAS .OR. NEWGEO) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &        ' Sending new orbital inf.',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C  /molinp/
         CALL MPIXBCAST(NMLINE,1,'INTEGER',MASTER)
C  /molinc/
#if defined (SYS_T3e)
         CALL MPIXBCAST(ITMP,10*NMLINE,'INTEGER',MASTER)
#else
         CALL MPIXBCAST(MLINE,80*NMLINE,'STRING',MASTER)
#endif
C  /inforb/
         CALL MPIXBCAST(MULD2H,NINFI,'INTEGER',MASTER)
C  /infpar/
         CALL MPIXBCAST(NODTOT,NPARI,'INTEGER',MASTER)
C  /dorps/
         CALL MPIXBCAST(DOREPS(0),NDORL,'LOGICAL',MASTER)
C  dosym()
         CALL MPIXBCAST(DOSYM,NSYML,'LOGICAL',MASTER)
C  comr12
         CALL MPIXBCAST(R12INT,NR12L,'LOGICAL',MASTER)
         CALL MPIXBCAST(INTGAC,NR12I,'INTEGER',MASTER)
      END IF
C
C  /parint/
      CALL MPIXBCAST(JATOM,NINTI,'INTEGER',MASTER)
      CALL MPIXBCAST(JNODV,NINTL,'LOGICAL',MASTER)
C
      CALL MPIXBCAST(NDMAT,1,'INTEGER',MASTER)
      CALL MPIXBCAST(DMAT,N2BASX*NDMAT,'DOUBLE',MASTER)
      CALL MPIXBCAST(ISYMDM,NFMAT,'INTEGER',MASTER)
      CALL MPIXBCAST(IFCTYP,NFMAT,'INTEGER',MASTER)
      CALL MPIXBCAST(HFXATT,1,'DOUBLE',MASTER)
C
      IF(ITYPE.EQ.3) THEN
        N2RED = NSYMBL*NSYMBL
        CALL MPIXBCAST(N2RED,1,'INTEGER',MASTER)
        CALL MPIXBCAST(DMRAO,N2RED*NDMAT,'DOUBLE',MASTER)
        CALL MPIXBCAST(GABRAO,N2RED,'DOUBLE',MASTER)
      ENDIF
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(FIRST ,1,'LOGICAL')
      CALL PVMXPACK(NEWBAS,1,'LOGICAL')
      CALL PVMXPACK(NEWGEO,1,'LOGICAL')
      CALL PVMXPACK(BASDIR,60,'STRING')
      CALL PVMXPACK(LENBAS,1,'INTEGER')
C
      IF (FIRST .OR. NEWBAS .OR. NEWGEO) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &        ' Sending new orbital inf.',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C  /molinc/
         CALL PVMXPACK(NMLINE,1,'INTEGER')
C  /molinp/
         CALL PVMXPACK(MLINE,80*NMLINE,'STRING')
C  /inforb/
         CALL PVMXPACK(MULD2H,NINFI,'INTEGER')
C  /infpar/
         CALL PVMXPACK(NODTOT,NPARI,'INTEGER')      
C  /dorps/
         CALL PVMXPACK(DOREPS(0),NDORL,'LOGICAL')      
C  dosym()
         CALL PVMXPACK(DOSYM,NSYML,'LOGICAL')      
C  comr12
         CALL PVMXBCAST(R12INT,NR12L,'LOGICAL')
         CALL PVMXBCAST(INTGAC,NR12I,'INTEGER')
      END IF
C  /parint/
      CALL PVMXPACK(JATOM,NINTI,'INTEGER')
      CALL PVMXPACK(JNODV,NINTL,'LOGICAL')
C
C     Pack the density matrices, their symmetries and the Fock matrix types.
C
      CALL PVMXPACK(NDMAT,1,'INTEGER')
      CALL PVMXPACK(DMAT,N2BASX*NDMAT,'DOUBLE')
      CALL PVMXPACK(ISYMDM,NFMAT,'INTEGER')
      CALL PVMXPACK(IFCTYP,NFMAT,'INTEGER')
C
      IF(ITYPE.EQ.3) THEN
        N2RED = NSYMBL*NSYMBL
        CALL PVMXPACK(N2RED,1,'INTEGER')
        CALL PVMXPACK(DMRAO,N2RED*NDMAT,'DOUBLE')
        CALL PVMXPACK(GABRAO,N2RED,'DOUBLE')
      ENDIF
C
C     -----------------------------------
C     Multicast initialization to nodes.
C     -----------------------------------
C
      CALL PVMXMCAST(NODTOT,NODEID(1),MTAG3)
C
#endif
      CALL QEXIT('HER_SDINIT')
      RETURN
      END
C  /* Deck her_rvres */
      SUBROUTINE HER_RVRES(FMAT,TMPMAT,HESSEE,TIMES,TSKCPU,IWHICH,
     &                     NTMPMT,NUMGRA,NUMHES,NUMSUS,NUMFCK,ITYPE,
     &                     IPRINT,DINTSKP)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <mtags.h>
      PARAMETER (D1 = 1.0D0)
      DIMENSION FMAT(*), TMPMAT(NTMPMT), TIMES(MTOTTK),
     &          TSKCPU(MTOTTK), IWHICH(MTOTTK),DSKPBF(8),
     &          DINTSKP(8), HESSEE(*)
C
C Used from common blocks:
C
C MXCENT : MXCENT (for ENERGY)
C MXORB  : MXSHEL (for INFPAR)
C ENERGY : GRADEE,HESSEE
C SUSCPT : SUS2EL
C INFORB : NNBASX,N2BASX
C INFPAR : TIMING, ???
C
#include <energy.h>
#include <suscpt.h>
#include <inforb.h>
#include <infpar.h>
C
      CALL QENTER('HER_RVRES')
      IF (IPRINT .GT. 3) CALL TITLER('Output from HER_RVRES','*',103)
C
C     Clean matrices.
C
      IF (TIMING) CALL DZERO(TIMES,MTOTTK)
C
      IF (ITYPE .EQ. 2) THEN
         CALL DZERO(GRADEE,NUMGRA)
         CALL DZERO(HESSEE,NUMHES)
      ELSE IF (ITYPE .EQ. 3) THEN
         CALL DZERO(FMAT,NUMFCK)
         CALL DZERO(DINTSKP,8)
#if defined(_CRAYT3E) && !defined(NO_BINSUM)
         CALL BIN_SUM(NUMFCK, FMAT)
         CALL BIN_SUM(8, DINTSKP)
#endif
      ELSE IF (ITYPE .EQ. -5) THEN
         CALL DZERO(FMAT,3*N2BASX)
         CALL DZERO(SUS2EL,NUMSUS)
      ELSE IF (ITYPE .EQ. 8) THEN
         CALL DZERO(FMAT,NTMPMT)
      END IF
C
      ICOUNT =  0
C
#if defined (VAR_MPI)
      DO 100 I = 1,NODTOT
         IWHO = -1
         CALL MPIXRECV(NWHO,1,'INTEGER',IWHO,MTAG7)
         IF (ITYPE .EQ. 2) THEN
            CALL MPIXRECV(TMPMAT(1),NUMGRA,'DOUBLE',NWHO,MTAG7)
            CALL MPIXRECV(TMPMAT(NUMGRA+1),NUMHES,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(NUMGRA,D1,TMPMAT,1,GRADEE,1)
            CALL DAXPY(NUMHES,D1,TMPMAT(NUMGRA+1),1,HESSEE,1)
#if !defined(_CRAYT3E) && !defined(NO_BINSUM)
         ELSE IF (ITYPE .EQ. 3) THEN
            CALL MPIXRECV(TMPMAT,NUMFCK,'DOUBLE',NWHO,MTAG7)
            CALL MPIXRECV(TMPMAT(1+NUMFCK),8,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(NUMFCK,D1,TMPMAT,1,FMAT,1)
            CALL DAXPY(8,D1,TMPMAT(NUMFCK+1),1,DINTSKP,1)
#endif
         ELSE IF (ITYPE .EQ. 8) THEN
            CALL MPIXRECV(TMPMAT(1),NTMPMT,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(NTMPMT,D1,TMPMAT,1,FMAT,1)
         ELSE IF (ITYPE .EQ. -5) THEN
            CALL MPIXRECV(TMPMAT(1),3*N2BASX,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(3*N2BASX,D1,TMPMAT,1,FMAT,1)
            CALL MPIXRECV(TMPMAT(1),NUMSUS,'DOUBLE',NWHO,MTAG7)
            CALL DAXPY(NUMSUS,D1,TMPMAT,1,SUS2EL,1)
         END IF
C
         CALL MPIXRECV(ITASKS,1,'INTEGER',NWHO,MTAG7)
C
         IF (TIMING) THEN
            CALL MPIXRECV(IWHICH(ICOUNT+1),ITASKS,'INTEGER',NWHO,MTAG7)
            CALL MPIXRECV(TSKCPU(ICOUNT+1),ITASKS,'DOUBLE' ,NWHO,MTAG7)
            DO 200 J=1,ITASKS
               TIMES(IWHICH(J)) = TSKCPU(J)
 200        CONTINUE
         END IF
C
         ICOUNT = ICOUNT + ITASKS
  100 CONTINUE
#endif
#if defined (VAR_PVM)
      DO 100 I = 1,NODTOT
         CALL PVMXRECV(-1,MTAG7)
         CALL PVMXUNPACK(TMPMAT,NTMPMT,'DOUBLE')
C
         IF (ITYPE .EQ. 2) THEN
            CALL DAXPY(NUMGRA,D1,TMPMAT(1),1,GRADEE,1)
            CALL DAXPY(NUMHES,D1,TMPMAT(NUMGRA+1),1,HESSEE,1)
C
         ELSE IF (ITYPE .EQ. 3) THEN
            CALL DAXPY(NUMFCK,D1,TMPMAT,1,FMAT,1)
            CALL DAXPY(8,D1,TMPMAT(NUMFCK+1),1,DINTSKP,1)
C
         ELSE IF (ITYPE .EQ. 8) THEN
            CALL DAXPY(NTMPMT,D1,TMPMAT(1),1,FMAT,1)
         ELSE IF (ITYPE .EQ. -5) THEN
            CALL DAXPY(3*N2BASX,D1,TMPMAT(1),1,FMAT,1)
            CALL DAXPY(NUMSUS,D1,TMPMAT(3*N2BASX+1),1,SUS2EL,1)
         END IF
C
C        Unpack the number of tasks received by the node.
C
         CALL PVMXUNPACK(ITASKS,1,'INTEGER')
C
         IF (TIMING) THEN
            CALL PVMXUNPACK(IWHICH(ICOUNT+1),ITASKS,'INTEGER')
            CALL PVMXUNPACK(TSKCPU(ICOUNT+1),ITASKS,'DOUBLE')
C
            DO 200 J=1,ITASKS
               TIMES(IWHICH(J)) = TSKCPU(J)
 200        CONTINUE
         END IF
C
         ICOUNT = ICOUNT + ITASKS
 100  CONTINUE
#endif
C
      IF (ICOUNT .NE. MTOTTK) THEN
        WRITE(LUPRI,'(/5X,A)') 'Error in parallel calculation!'
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks to calculate  : ',
     &       MTOTTK
        WRITE(LUPRI,'(5X,A,I5)') 'Number of tasks    calculated : ',
     &       ICOUNT
c        CALL QUIT('Inconsistence in RVRES -> MTOTTK .NE. ICOUNT')
      ENDIF
C
      CALL QEXIT('HER_RVRES')
      RETURN
      END
C  /* Deck parsrt */
      SUBROUTINE PARSRT(MTOTTK,NPOS,SORTED,TIMES,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION  TIMES(MTOTTK), NPOS(0:MTOTTK+1), SORTED(0:MTOTTK+1)        
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from PARSRT','*',103)
C
C------------------------------------------------------
C     The array TIMES containing CPU-times is sorted by
C     increasing CPUs and stored in the array SORTED.
C------------------------------------------------------
C
C     SORTED(0) must be less than any possible CPU-time,
C     and TIMMAX must be greater than any possible CPU-time
C
      ITEMS     = MTOTTK
      TIMMAX    = 99999D00
      SORTED(0) =  -1.0D00
      SORTED(1) = TIMES(1)
      NPOS  (0) = 0
      NPOS  (1) = 1
C
      DO 100,ITEM = 2,ITEMS
C
         LASTLO       = 0
         LASTHI       = ITEM
         NUMBER       = INT(ITEM/2)
         SORTED(ITEM) = TIMMAX
C
 200     CONTINUE
C
         IF (TIMES(ITEM) .EQ. SORTED(NUMBER) .OR.
     &      (TIMES(ITEM) .GT. SORTED(NUMBER) .AND.
     &       TIMES(ITEM) .LT. SORTED(NUMBER+1))) THEN
C
            DO 300, NCOUNT = ITEM,NUMBER,-1
               NPOS  (NCOUNT+1) = NPOS  (NCOUNT)
               SORTED(NCOUNT+1) = SORTED(NCOUNT)
 300        CONTINUE
C
            NPOS  (NUMBER+1) = ITEM
            SORTED(NUMBER+1) = TIMES(ITEM)
            GO TO 100
C
         ELSE IF (TIMES(ITEM) .LT. SORTED(NUMBER)) THEN
C
            LASTHI = NUMBER
            NUMBER = INT((LASTLO+LASTHI)/2)
            GO TO 200
C
         ELSE IF (TIMES(ITEM) .GT. SORTED(NUMBER)) THEN
C
            LASTLO = NUMBER
            NUMBER = INT((LASTLO+LASTHI)/2)
            GO TO 200
C
         END IF
 100  CONTINUE
C
      IF (IPRINT .GT. 4) THEN
         DO 400 I = 1, ITEMS
            WRITE(LUPRI,'(15X,A4,I5,A9,F8.3)')
     &           'IJ =', NPOS(I), '   time =', SORTED(I)
 400      CONTINUE
      END IF
C
      RETURN
      END
C  /* Deck her_nodstr */
      SUBROUTINE HER_NODSTR(WORK,LWORK,IPRINT)
C
C    *****************************************************************
C    *    This is the node program for the construction of fock-     *
C    *   matrices, derivated fock-matrices and expectation values.   *
C    *****************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <dummy.h>
C
      PARAMETER (MXDMAT = 2*MXCOOR)
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME
      DIMENSION WORK(LWORK), ISYMDM(MXDMAT), IFCTYP(MXDMAT),
     &          DINTSKP(2,4)
C
C Used from common blocks
C
C MXORB  : MXSHEL (for INFPAR)
C MXCENT : MXCOOR (for ENERGY)
C ENERGY : GRADEE, HESSEE
C SUSCPT : SUS2EL
C INFORB : N2BASX
C INFPAR : SLAVE
C
#include <energy.h>
#include <suscpt.h>
#include <inforb.h>
#include <infpar.h>
#include <gnrinf.h>
C
#include <memint.h>
C
      CALL QENTER('HER_NODSTR')
C
C     Start timing
C
      CALL GETTIM(CPU1,WALL1)
C
C     This is a slave in a parallel run
C
      SLAVE = .TRUE.
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from HER_NODSTR','*',103)
C
C     Receive initialization from master
C     ==================================
C
      NPAO = MXSHEL*MXAOVC
      CALL MEMGET('INTE',KJSTRS,NPAO*2,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KNPRIM,NPAO*2,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KNCONT,NPAO*2,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KIORBS,NPAO  ,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KJORBS,NPAO  ,WORK,KFREE,LFREE)     
      CALL MEMGET('INTE',KKORBS,NPAO  ,WORK,KFREE,LFREE)     
      CALL HER_RVINIT(WORK,KFREE,LFREE,
     &                NDMAT,ISYMDM,IFCTYP,ITYPE,IATOM,
     &                MAXDIF,NODV,NOPV,NOCONT,RETUR,TKTIME,
     &                WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),
     &                WORK(KIORBS),WORK(KJORBS),WORK(KKORBS),IPRINT,
     &                I2TYP,ICEDIF,IFTHRS,
     &                N2RED,KDMAT,KFMAT,KHESEE,KDMRAO,KGABAO,KDMRSO)
C
C     Initialize matrices
C     ===================
C
      CALL DZERO(WORK(KFMAT),NFMAT*N2BASX)
C
      IF (ITYPE .EQ. 2) THEN
         CALL DZERO(GRADEE,MXCOOR)
         CALL DZERO(WORK(KHESEE),MXCOOR*MXCOOR)
      ELSE IF (ITYPE .EQ. -5) THEN
         CALL DZERO(SUS2EL,9)
      END IF
C
C     Calculate the integrals
C     =======================
C
      GMAT   = 0
      INDXAB = 0
      NUMDIS = 0
      MAXDIS = 1
      IPRNTA = 0
      IPRNTB = 0
      IPRNTC = 0
      IPRNTD = 0
      ISHLA  = 0
C
      CALL TWOINT(WORK(KFREE),LFREE,WORK(KHESEE),WORK(KFMAT),
     &     WORK(KDMAT),NDMAT,
     &     ISYMDM,IFCTYP,GMAT,INDXAB,NUMDIS,MAXDIS,ITYPE,MAXDIF,
     &     IATOM,NODV,NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,
     &     IPRNTB,IPRNTC,IPRNTD,RETUR,ISHLA,I2TYP,WORK(KJSTRS),
     &     WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &     ICEDIF,IFTHRS,WORK(KGABAO),WORK(KDMRAO),
     &     WORK(KDMRSO),DINTSKP,RELCAL)
      CALL MEMREL('HER_NODSTR',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Send overall timings if requested.
C
      IF (IPRPAR .GT. 0) CALL SDTIM(CPU1,WALL1,IPRINT)
C
      CALL QEXIT('HER_NODSTR')
C
      RETURN
      END
C  /* Deck her_rvinit */
      SUBROUTINE HER_RVINIT(WORK,KFREE,LFREE,
     &                      NDMAT,ISYMDM,IFCTYP,ITYPE,
     &                      IATOM,MAXDIF,NODV,NOPV,NOCONT,
     &                      RETUR,TKTIME,JSTRSH,NPRIMS,NCONTS,IORBSH,
     &                      JORBSH,KORBSH,IPRINT,I2TYP,ICEDIF,IFTHRS,
     &                      N2RED,KDMAT,KFMAT,KHESEE,KDMRAO,KGABAO,
     &                      KDMRSO)
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
#include <mtags.h>
      PARAMETER (D0 = 0.0D0)
      LOGICAL NODV,NOPV,NOCONT,RETUR,TKTIME,FIRST
      DIMENSION WORK(*), ISYMDM(*), IFCTYP(*), JSTRSH(*), NPRIMS(*),
     &          NCONTS(*), IORBSH(*), JORBSH(*), KORBSH(*)
#include <infpar.h>
C
C Used from common blocks:
C
C MXCENT : MXCENT (for DORPS, etc.)
C MXORB  : MXSHEL (for INFPAR)
C MAXAQN : *      (for SYMMET)
C GNRINF : NEWBAS, (PVM: BASDIR, LENBAS)
C SIRIPC : NEWGEO
C INFORB : N2BASX
C MOLINP : NMLINE, MLINE()
C DORPS  : DOREPS(), DOCOOR()
C ABAINF : DOSYM()
C SYMMET : MAXREP
C NUCLEI : NUCDEG()
C INFPAR : *
C HRUNIT : *
C PARINT : *
C
#include <gnrinf.h>
#include <siripc.h>
#include <inforb.h>
#include <molinp.h>
#include <dorps.h>
#include <abainf.h>
#include <symmet.h>
#include <nuclei.h>
#include <dftcom.h>
#include <inftap.h>
#include <parint.h>
#include <r12int.h>
C
#if defined (SYS_T3E)
C
C     This kludge fixes the current non-implementation of MPI_CHARACTER on
C     the Cray-T3E, K.Ruud-Nov.96
C
      DIMENSION ITMP(4000)
      EQUIVALENCE (ITMP,MLINE)
#endif
#include <cbisol.h>
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from HER_RVINIT','*',103)
C
C-----------------------------------------
C     Receive initializations from master.
C-----------------------------------------
C
#if defined (VAR_MPI)
      CALL MPIXBCAST(FIRST ,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWBAS,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(NEWGEO,1,'LOGICAL',MASTER)
      CALL MPIXBCAST(HFXFAC,1,'DOUBLE',MASTER)
      CALL MPIXBCAST(HFXMU, 1,'DOUBLE',MASTER)
C
      TIMING = FIRST .OR. NEWBAS .OR. NEWGEO
C
      IF (TIMING) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &        ' Receiving new orbital inf.',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C  /molinp/
            CALL MPIXBCAST(NMLINE,1,'INTEGER',MASTER)
C  /molinc/
#if defined (SYS_T3e)
            CALL MPIXBCAST(ITMP,10*NMLINE,'INTEGER',MASTER)
#else
            CALL MPIXBCAST(MLINE,80*NMLINE,'STRING',MASTER)
#endif
C  /inforb/
            CALL MPIXBCAST(MULD2H,NINFI,'INTEGER',MASTER)
C  /infpar/
            CALL MPIXBCAST(NODTOT,NPARI,'INTEGER',MASTER)
C  /dorps/
            CALL MPIXBCAST(DOREPS(0),NDORL,'LOGICAL',MASTER)
C  dosym()
            CALL MPIXBCAST(DOSYM,NSYML,'LOGICAL',MASTER)
C  /comr12/
            CALL MPIXBCAST(R12INT,NR12L,'LOGICAL',MASTER)
            CALL MPIXBCAST(INTGAC,NR12I,'INTEGER',MASTER)
      END IF
C  /parint/
      CALL MPIXBCAST(JATOM,NINTI,'INTEGER',MASTER)
      CALL MPIXBCAST(JNODV,NINTL,'LOGICAL',MASTER)
#endif
#if defined (VAR_PVM)
      CALL PVMXRECV(MASTER,MTAG3)
      CALL PVMXUNPACK(FIRST ,1,'LOGICAL')
      CALL PVMXUNPACK(NEWBAS,1,'LOGICAL')
      CALL PVMXUNPACK(NEWGEO,1,'LOGICAL')
      CALL PVMXUNPACK(BASDIR,60,'STRING')
      CALL PVMXUNPACK(LENBAS,1,'INTEGER')
C
C     --------------------------------------------------
C     Unpack common block MOLINP and INFORB (only once).
C     --------------------------------------------------
C
      TIMING = FIRST .OR. NEWBAS .OR. NEWGEO
C
      IF (TIMING) THEN
         IF (IPRINT .GT. 3) WRITE(LUPRI,'(/1X,A,/,1X,A,3L3)')
     &        ' Receiving new orbital inf.',
     &        ' FIRST,NEWBAS,NEWGEO :',FIRST,NEWBAS,NEWGEO
C  /molinp/
         CALL PVMXUNPACK(NMLINE,1,'INTEGER')
C  /molinc/
         CALL PVMXUNPACK(MLINE,80*NMLINE,'STRING')
C  /inforb/
         CALL PVMXUNPACK(MULD2H,NINFI,'INTEGER')
C  /infpar/
         CALL PVMXUNPACK(NODTOT,NPARI,'INTEGER')
C
         DO 100 I = 1, NODTOT
            IF (NODEID(I) .EQ. MYTID) MYNUM = I
 100     CONTINUE
C  /dorps/
         CALL PVMXUNPACK(DOREPS(0),NDORL,'LOGICAL')
C  dosym()
         CALL PVMXUNPACK(DOSYM,NSYML,'LOGICAL')
C  comr12
         CALL PVMXBCAST(R12INT,NR12L,'LOGICAL')
         CALL PVMXBCAST(INTGAC,NR12I,'INTEGER')
      END IF
C  /parint/
      CALL PVMXUNPACK(JATOM,NINTI,'INTEGER')
      CALL PVMXUNPACK(JNODV,NINTL,'LOGICAL')
#endif
C
C     Retain information from PARINT.
C
      IATOM  = JATOM
      ITYPE  = JTYPE
      LUDASP = JLUDAS
      LUINTA = JLUINT
      LUONEL = JLUONE
      LUSOL  = JLUSOL
      LUSUPM = JLUSUP
      MAXDIF = JMXDIF
      MAXREP = JMXREP
      NFMAT  = JNFMAT
      NTASK  = JTASK
      I2TYP  = J2TYP
      ICEDIF = JCEDIF
      IFTHRS = JFTHRS
      NOCONT = JNOCNT
      NODV   = JNODV 
      NOPV   = JNOPV 
      RETUR  = JRETUR
      TKTIME = JTKTIM
      SOLVNT = JSOLVN
      RELCAL = JRELCL
C
C     Set hermit.
C
      IF (NEWGEO .AND. SOLVNT) THEN
         NEWGEO = .TRUE.
         NUCIND = NUCIND + 1
         NUCDEP = NUCDEP + 1
         NATOMS = NATOMS + 1
         NCNTCV = NUCIND
         NCLINE(NUCIND) = 0
         NAMN(NUCIND)       = 'cav '
         NAMEX(3*NUCIND-2)  = 'cav  x'
         NAMEX(3*NUCIND-1)  = 'cav  y'
         NAMEX(3*NUCIND)    = 'cav  z'
         NAMDEP(NUCDEP)     = 'cavity'
         NAMDPX(3*NUCDEP-2) = 'cavity x'
         NAMDPX(3*NUCDEP-1) = 'cavity y'
         NAMDPX(3*NUCDEP  ) = 'cavity z'
         IF (NUCDEP .GT. MXCENT) THEN
            WRITE (LUPRI,'(//2A,/A,I5)')
     &         ' Too many atomic centers: MXCENT exceed in READIN for',
     &         ' solvent cavity,',' Current limit:',MXCENT
            CALL QUIT('*** ERROR *** MXCENT exceeded in READIN')
         END IF
         CORD(1,NUCIND) = D0
         CORD(2,NUCIND) = D0
         CORD(3,NUCIND) = D0
         ISTBNU(NUCIND) = 7
         CHARGE(NUCIND) = D0
         CALL NUCPRO(WORK(KFREE),LFREE)
         NEWGEO = .FALSE.   
      END IF
C
      IF (NEWGEO) RDINPC = .FALSE.
C
      CALL SETHER(0,NEWGEO,WORK(KFREE),LFREE)
C
C     Set common-block BLOCKS
C
      CALL PAOVEC(JSTRSH,NPRIMS,NCONTS,IORBSH,JORBSH,KORBSH,0,.FALSE.,0)
C
#if defined (VAR_MPI)
      CALL MPIXBCAST(NDMAT,1,'INTEGER',MASTER)
      CALL MEMGET('REAL',KDMAT,NDMAT*N2BASX,WORK,KFREE,LFREE)     
      CALL MEMGET('REAL',KFMAT,NFMAT*N2BASX,WORK,KFREE,LFREE)     
      CALL MEMGET('REAL',KHESEE,MXCOOR*MXCOOR,WORK,KFREE,LFREE)     
      CALL MPIXBCAST(WORK(KDMAT),NDMAT*N2BASX,'DOUBLE',MASTER)
      CALL MPIXBCAST(ISYMDM,NFMAT,'INTEGER',MASTER)
      CALL MPIXBCAST(IFCTYP,NFMAT,'INTEGER',MASTER)
      CALL MPIXBCAST(HFXATT,1,'DOUBLE',MASTER)
      IF(ITYPE.EQ.3) THEN
        CALL MPIXBCAST(N2RED,1,'INTEGER',MASTER)
        CALL MEMGET('REAL',KDMRAO,NDMAT*N2RED,WORK,KFREE,LFREE)     
        CALL MEMGET('REAL',KGABAO,N2RED,WORK,KFREE,LFREE)     
        CALL MPIXBCAST(WORK(KDMRAO),NDMAT*N2RED,'DOUBLE',MASTER)
        CALL MPIXBCAST(WORK(KGABAO),N2RED,'DOUBLE',MASTER)
      ELSE
        KDMRAO = KFREE
        KGABAO = KFREE
      ENDIF
      KDMRSO = KFREE
#endif
#if defined (VAR_PVM)
C     Unpack density matrices.
C
      CALL PVMXUNPACK(NDMAT,1,'INTEGER')
      CALL MEMGET('REAL',KDMAT,NDMAT*N2BASX,WORK,KFREE,LFREE)     
      CALL MEMGET('REAL',KFMAT,NFMAT*N2BASX,WORK,KFREE,LFREE)     
      CALL MEMGET('REAL',KHESEE,MXCOOR*MXCOOR,WORK,KFREE,LFREE)     
      CALL PVMXUNPACK(WORK(KDMAT),N2BASX*NDMAT,'DOUBLE')
      CALL PVMXUNPACK(ISYMDM,NFMAT,'INTEGER')
      CALL PVMXUNPACK(IFCTYP,NFMAT,'INTEGER')
C
      IF(ITYPE.EQ.3) THEN
        CALL PVMXUNPACK(N2RED,1,'INTEGER')
        CALL MEMGET('REAL',KDMRAO,NDMAT*N2RED,WORK,KFREE,LFREE)     
        CALL MEMGET('REAL',KGABAO,N2RED,WORK,KFREE,LFREE)     
        CALL PVMXUNPACK(WORK(KDMRAO),N2RED*NDMAT,'DOUBLE')
        CALL PVMXUNPACK(WORK(KGABAO),N2RED,'DOUBLE')
        KLAST  = KGABAO + N2RED
      ELSE
        KDMRAO = KFREE
        KGABAO = KFREE
      ENDIF
      KDMRSO = KFREE
#endif
C
      RETURN
      END
C  /* Deck parlop */
      SUBROUTINE PARLOP(WORK,LWORK,HESSEE,FMAT,DMAT,NDMAT,GMAT,
     &                  MAXDER,EXPECT,SUSCEP,UNDIFF,DDFOCK,DIRFCK,
     &                  SOFOCK,DISTRI,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                  PERTUR,IATOM,MULE,MULTE,NODV,NOPV,NOCONT,
     &                  THRESH,JPRINT,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RETUR,
     &                  SQ12EL,INDHER,INDHSQ,IODDHR,TSKCPU,IWHICH,IJS,
     &                  ISYMDM,IFCTYP,ADISTR,ITYPE,JSTRSH,NPRIMS,NCONTS,
     &                  IORBSH,I2TYP,ICEDIF,IFTHRS,
     &                  GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL)
C
C     Copied from TWOLOP and rewritten by Paal Dahle Nov.1994
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxorb.h> 
#include <dummy.h>
#include <mtags.h>
      LOGICAL PRINTA, PRINTB, PRINTC, PRINTD, NOPV, NODV, PERTUR,
     &        EXPECT, UNDIFF, DDFOCK, DIRFCK, DIA2SO, ZFS2EL,
     &        DISTRI, NOCONT, SPNORB,
     &        RETUR, FIRST, SQ12EL, LONDON, SUSCEP, ADISTR, DONE,
     &        SOFOCK,RELCAL
      DIMENSION DMAT(*), FMAT(*), GMAT(*), INDHSQ(*), IODDHR(*),
     &          INDHER(*), WORK(LWORK), ISYMDM(*), IFCTYP(*),
     &          JSTRSH(*), NPRIMS(*), NCONTS(*), IORBSH(MXSHEL,MXAOVC),
     &          TSKCPU(MTOTTK), IWHICH(MTOTTK), IJS(NTASK),
     &          GABRAO(*), DMRAO(*),DMRSO(*),DINTSKP(2,4), HESSEE(*)
#include <cbisol.h>
#include <twocom.h>
#include <nuclei.h>
#include <energy.h>
#include <taymol.h>
#include <taysol.h>
#include <suscpt.h>
#include <blocks.h>
#include <symmet.h>
#include <infpar.h>
#include <inftap.h>
#include <inforb.h>
C
      CALL QENTER('PARLOP')
      IF (JPRINT .GT. 3) CALL TITLER('Output from PARLOP','*',103)
C
      FIRST = .TRUE.
      DIRAC = RELCAL
      IF (EXPECT .AND. .NOT.NOPV) THEN
         CALL REWSPL(LUPAO)
      END IF
      IF (SUSCEP) THEN
         IF (.NOT.NOPV) CALL REWSPL(LUPAO)
         IF (.NOT.NOPV) CALL REWSPL(LUPAS)
         CALL DZERO(SUS2EL,9)
      END IF
C
C
C     For direct contributions: some parameters
C       FCKTHR is threshold for screening
C       ICEFLG gives information about separate screening of
C              Coulomb/exchange for each DMAT
C       NCM    is the number of DMAT requiring Coulomb-contributions
C       NEM    is the number og DMAT requiring exchange-contributions
C       Screening proceeds in three steps as documented by DINTSKP:
C         Step 1: Screening on integral batches
C           DINTSKP(1,1) - total number of integrals
C           DINTSKP(2,1) - number of integrals skipped (batchwise)
C         Step 2: Screening on individual integrals 
C                 while unpacking indices
C           DINTSKP(1,2) - number of integrals remaining after step 1
C           DINTSKP(2,2) - number of integrals skipped
C         Step 3a: Screening on Coulomb contributions
C           DINTSKP(1,3) - NCM times number of integrals remaining 
C                         after step 2
C           DINTSKP(2,3) - NCM times number of integrals skipped
C         Step 3b: Screening on exchange contributions
C           DINTSKP(1,4) - NEM times number of integrals remaining 
C                         after step 2
C           DINTSKP(2,4) - NEM times number of integrals skipped
C
      DOSCRN = .FALSE.
      CALL DZERO(DINTSKP,8)
      IF(DIRFCK.OR.SOFOCK)  THEN
        IF(IFTHRS.LT.16) THEN
          DOSCRN = .TRUE.
          FCKTHR = -IFTHRS
          FCKTHR = 10.0D0**FCKTHR
          ICEFLG = ICEDIF
          NCM = 0
          NEM = 0
          DO I = 1,NDMAT
            IY  = MOD(IFCTYP(I),10)
            IC  = MOD(IY,2)
            NCM = NCM + IC
            IE  = (IY - IC)/2
            NEM = NEM + IE
          ENDDO
        ENDIF
      ENDIF
C
      IF(I2TYP.EQ.0) THEN
        IASTRT = 1
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSEIF(I2TYP.EQ.1) THEN
        IASTRT = 1
        IBSTRT = 1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = NLRGBL
        IBSMAX = NLRGBL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.2) THEN
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        ICSTRT = 1
        IDSTRT = 1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = NLRGBL
        IDSMAX = NLRGBL
      ELSEIF(I2TYP.EQ.3) THEN
        IASTRT = NLRGBL+1
        IBSTRT = NLRGBL+1
        ICSTRT = NLRGBL+1
        IDSTRT = NLRGBL+1
        IASMAX = MAXSHL
        IBSMAX = MAXSHL
        ICSMAX = MAXSHL
        IDSMAX = MAXSHL
      ELSE
        WRITE(LUPRI,'(A,I5)') 'TWOLOP: Unknown I2TYP =' ,I2TYP
        CALL QUIT('Unknown I2TYP !!!')
      ENDIF
C
C     ***********************************
C     ***** Send request for a task *****
C     ***********************************
C
      ICOUNT = 0
      CALL IZERO(IWHICH,MTOTTK)
      CALL DZERO(TSKCPU,MTOTTK)
C
  100 CONTINUE
C
#if defined (VAR_MPI)
      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MTAG5)
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      CALL PVMXPACK(MYNUM,1,'INTEGER')
      CALL PVMXSEND(MASTER,MTAG5)
#endif
C
C     Receive IJ labels for NTASK new (IJ|**) integrals 
C
      CALL IZERO(IJS,NTASK)
C
#if defined (VAR_MPI)
      CALL MPIXRECV(DONE,1,'LOGICAL',MASTER,MTAG6)
      IF (DONE) GOTO 300
      CALL MPIXRECV(IJS,NTASK,'INTEGER',MASTER,MTAG6)
#endif
#if defined (VAR_PVM)
      CALL PVMXRECV(MASTER,MTAG6)
      CALL PVMXUNPACK(DONE,1,'LOGICAL')
      IF (DONE) GOTO 300
      CALL PVMXUNPACK(IJS,NTASK,'INTEGER')
#endif
C
      IF (JPRINT .GT. 4) THEN
         WRITE(LUPRI,'(5X,A,I4)')
     &      'Receiving IJ-batches: ',(IJS(KK),KK=1,NTASK)
      END IF
C
      DO 200, I = 1,NTASK
C
         IJSHEL = IJS(I)
C
         IF (IJSHEL .EQ. 0) THEN
            ICOUNT = ICOUNT + I - 1
            GO TO 100
         END IF
C
C        The compressed index IJSHEL is split to indices ISHELA and ISHELB.
C
         CALL UNPKIJ(IJSHEL,ISHELA,ISHELB)
C
C        Check if this distribution is requested, otherwise
C        go for next task [ FOR EFFICIENCY THIS OUGHT NOT TO
C        HAPPEN! ]  950428-K.Ruud+H.J.Aa.Jensen
C
         ICA = LCLASH(ISHELA)
         ICB = LCLASH(ISHELB)
C
         IF (TIMING) CPUTK1 = SECOND()
C
C        *****************************
C        ***** First Shell Index *****
C        *****************************
C
         NHKTA  = NHKTSH(ISHELA)
         KHKTA  = KHKTSH(ISHELA)
         KCKTA  = KCKTSH(ISHELA)
         SPHRA  = SPHRSH(ISHELA)
         NCENTA = NCNTSH(ISHELA)
         MULA   = ISTBSH(ISHELA)
         MULTA  = MULT(MULA)
         NSTRA  = IORBSB(IORBSH(ISHELA,1))
         NUCA   = NUCOSH(ISHELA)
         NORBA  = NORBSH(ISHELA)
         IF (.NOT.BIGVEC) THEN
            CORAX0 = CENTSH(ISHELA,1)
            CORAY0 = CENTSH(ISHELA,2)
            CORAZ0 = CENTSH(ISHELA,3)
         END IF
         PRINTA = .TRUE.
         IF ((ISHELA .NE. IPRNTA).AND.(IPRNTA .NE. 0)) PRINTA = .FALSE.
C
C        ******************************
C        ***** Second Shell Index *****
C        ******************************
C
         NHKTB  = NHKTSH(ISHELB)
         KHKTB  = KHKTSH(ISHELB)
         KCKTB  = KCKTSH(ISHELB)
         SPHRB  = SPHRSH(ISHELB)
         NCENTB = NCNTSH(ISHELB)
         MULB   = ISTBSH(ISHELB)
         MULTB  = MULT(MULB)
         NSTRB  = IORBSB(IORBSH(ISHELB,1))
         NUCB   = NUCOSH(ISHELB)
         NORBB  = NORBSH(ISHELB)
         IF (.NOT.BIGVEC) THEN
            CORBX0 = CENTSH(ISHELB,1)
            CORBY0 = CENTSH(ISHELB,2)
            CORBZ0 = CENTSH(ISHELB,3)
         END IF
         GENAB  = .NOT.(SEGMSH(ISHELA) .AND. SEGMSH(ISHELB))
         IGENAB = 1
         IF (.NOT.GENAB) IGENAB = 2
         NSETA  = NSETSH(ISHELA,IGENAB)
         NSETB  = NSETSH(ISHELB,IGENAB)
         PRINTB = PRINTA
         IF ((ISHELB.NE.IPRNTB).AND.(IPRNTB.NE.0)) PRINTB = .FALSE.
C
C        *****************************
C        ***** Third Shell Index *****
C        *****************************
C
         ICMAX = ISHELA
         IF (SPNORB) ICMAX = MAXSHL
         IF(I2TYP.EQ.2) ICMAX = NLRGBL
         DO 400 ISHELC = ICSTRT, ICMAX
C
            ICC    = LCLASH(ISHELC)
C
            NHKTC  = NHKTSH(ISHELC)
            KHKTC  = KHKTSH(ISHELC)
            KCKTC  = KCKTSH(ISHELC)
            SPHRC  = SPHRSH(ISHELC)
            NCENTC = NCNTSH(ISHELC)
            MULC   = ISTBSH(ISHELC)
            MULTC  = MULT(MULC)
            NSTRC  = IORBSB(IORBSH(ISHELC,1))
            NUCC   = NUCOSH(ISHELC)
            NORBC  = NORBSH(ISHELC)
            IF (.NOT.BIGVEC) THEN
               CORCX0 = CENTSH(ISHELC,1)
               CORCY0 = CENTSH(ISHELC,2)
               CORCZ0 = CENTSH(ISHELC,3)
            END IF
            PRINTC = PRINTB
            IF ((ISHELC.NE.IPRNTC).AND.(IPRNTC.NE.0)) PRINTC=.FALSE.
C
C           ******************************
C           ***** Fourth Shell Index *****
C           ******************************
C
            IDMAX = ISHELC
            IF (.NOT.SPNORB.AND.(ISHELA.EQ.ISHELC)) IDMAX = ISHELB
            DO 500 ISHELD = IDSTRT,IDMAX
C
               ICD    = LCLASH(ISHELD)
C
               NHKTD  = NHKTSH(ISHELD)
               KHKTD  = KHKTSH(ISHELD)
               KCKTD  = KCKTSH(ISHELD)
               SPHRD  = SPHRSH(ISHELD)
               NCENTD = NCNTSH(ISHELD)
               MULD   = ISTBSH(ISHELD)
               MULTD  = MULT(MULD)
               NSTRD  = IORBSB(IORBSH(ISHELD,1))
               NUCD   = NUCOSH(ISHELD)
               NORBD  = NORBSH(ISHELD)
               IF (.NOT.BIGVEC) THEN
                  CORDX0 = CENTSH(ISHELD,1)
                  CORDY0 = CENTSH(ISHELD,2)
                  CORDZ0 = CENTSH(ISHELD,3)
               END IF
               GENCD = .NOT.(SEGMSH(ISHELC) .AND. SEGMSH(ISHELD))
               IGENCD = 1
               IF (.NOT.GENCD) IGENCD = 2
               NSETC = NSETSH(ISHELC,IGENCD)
               NSETD = NSETSH(ISHELD,IGENCD)
               PRINTD = PRINTC
               IF ((ISHELD .NE. IPRNTD).AND.(IPRNTD .NE. 0))
     &              PRINTD = .FALSE.
C
               SHAEQB = ISHELA .EQ. ISHELB
               SHCEQD = ISHELC .EQ. ISHELD
               SHABAB = (ISHELA.EQ.ISHELC) .AND. (ISHELB.EQ.ISHELD)
C
C              *******************************
C              ***** Calculate integrals *****
C              *******************************
C
c               IPRINT = 0
               IPRINT = JPRINT
C
C
C     Note DUMMY should be HESSEE, K.Ruud-96
C
               CALL TWOODS(FMAT,DMAT,NDMAT,GMAT,HESSEE,WORK,LWORK,
     &                 UNDIFF,PERTUR,LONDON,SPNORB,DIA2SO,ZFS2EL,
     &                 EXPECT,SUSCEP,DDFOCK,DIRFCK,SOFOCK,DISTRI,
     &                 IATOM,MULE,MULTE,
     &                 MAXDER,NOCONT,NODV,NOPV,THRESH,IPRINT,
     &                 FIRST,SQ12EL,INDHSQ,IODDHR,INDHER,IFCTYP,
     &                 ADISTR,JSTRSH,NPRIMS,NCONTS,IORBSH,DUMMY,
     &                 ICEDIF,IFTHRS,DINTSKP,GABRAO,DMRAO,
     &                 DMRSO,IREPDM)
C
               IF (RETUR) THEN
                  IF (ISHELA .EQ. IPRNTA .AND.
     &                 ISHELB .EQ. IPRNTB .AND.
     &                 ISHELC .EQ. IPRNTC .AND.
     &                 ISHELD .EQ. IPRNTD) GOTO 9999
               END IF
 500        CONTINUE
 400     CONTINUE
C
         IF (TIMING) THEN
            CPUTK2 = SECOND()
            CPUTSK = CPUTK2 - CPUTK1
            IPLACE = ICOUNT + I
C
            IWHICH(IPLACE) = IJSHEL
            TSKCPU(IPLACE) = CPUTSK
         END IF
 200  CONTINUE
C
      ICOUNT = ICOUNT + NTASK
      GOTO 100
C
 300  CONTINUE
C
C     Note that DUMMY here should be HESSEE,K.Ruud, Jan.-97
C
CTROND      IF (DIRFCK .OR. DDFOCK .OR. EXPECT) THEN
CTROND         CALL SKLFCK(FMAT,DUMMY,WORK,LWORK,JPRINT,DIRFCK,DDFOCK,EXPECT,
CTROND     &               PERTUR,NODV,MAXDER,LONDON,NDMAT,ISYMDM,IFCTYP,
CTROND     &               IATOM,.FALSE.)
CTROND      END IF
C
C     <<<<< Print Section - Gradient and Hessian Elements >>>>>
C
      IF (EXPECT) THEN
         IF (JPRINT .GT. 5) THEN
            KCSTRA = 1
            KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
            KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
            IF (KLAST .GT. LWORK) CALL STOPIT('PARLOP',' ',KLAST,LWORK)
            CALL HEADER('Two-electron integral gradient',-1)
            CALL PRIGRD(GRADEE,WORK(KCSTRA),WORK(KSCTRA))
            CALL HEADER('Potential energy (NN + NE + EE) gradient',-1)
            CALL ZERGRD
            CALL ADDGRD(GRADNN)
            CALL ADDGRD(GRADNA)
            CALL ADDGRD(GRADEE)
            CALL PRIGRD(GRDMOL,WORK(KCSTRA),WORK(KSCTRA))
            CALL HEADER('Molecular gradient',-1)
            CALL ADDGRD(GRADFS)
            CALL ADDGRD(GRADKE)
            IF (SOLVNT) THEN
               CALL ADDGRD(GSOLTT)
               CALL ADDGRD(GSOLNN)
            END IF
            CALL PRIGRD(GRDMOL,WORK(KCSTRA),WORK(KSCTRA))
            NCDEP3 = 3*NUCDEP
            GRDNRM = DDOT(NCDEP3,GRDMOL,1,GRDMOL,1)
            GRDNRM = SQRT(GRDNRM)
            WRITE (LUPRI,'(/19X,A,1P,E10.2)')
     *         'Molecular gradient norm:', GRDNRM
            CALL ZERGRD
         END IF
      END IF
      IF (LONDON .AND. MAXDER.EQ.2) THEN
         SUS2EL(2,1) = SUS2EL(1,2)
         SUS2EL(3,1) = SUS2EL(1,3)
         SUS2EL(3,2) = SUS2EL(2,3)
         IF (JPRINT .GT. 5) THEN
            CALL HEADER('Two-electron integral susceptibilities',-1)
            CALL OUTPUT(SUS2EL,1,3,1,3,3,3,1,LUPRI)
         END IF
      END IF
C
C     Print Fock matrices
C
      IF (DIRFCK .AND. JPRINT.GT.5) THEN
         CALL HEADER('Fock matrix in PARLOP',-1)
         DO 600 I = 1, NDMAT
            ISTR = NBAST*NBAST*(I - 1) + 1
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(ISTR),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  600    CONTINUE
      END IF
C
C     Send final results to master
C
      CALL HER_SDRES(FMAT,HESSEE,NDMAT,DINTSKP,TSKCPU,IWHICH,ICOUNT,
     &               ITYPE,JPRINT)
C
 9999 CALL QEXIT('PARLOP')
      RETURN
      END
C  /* Deck her_sdres */
      SUBROUTINE HER_SDRES(FMAT,HESSEE,NDMAT,DINTSKP,TSKCPU,IWHICH,
     &                     ICOUNT,ITYPE,IPRINT)
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <mtags.h>
      DIMENSION FMAT(NFMAT*N2BASX), TSKCPU(ICOUNT), IWHICH(ICOUNT),
     &          DINTSKP(8), HESSEE(*)
C
C Used from common blocks:
C
C MXCENT: MXCOOR (for ENERGY)
C MXORB : MXSHEL (for INFPAR)
C ENERGY: GRADEE,HESSEE
C SUSCPT: SUS2EL
C INFORB: N2BASX
C INFPAR: MTOTTK, MASTER, TIMING
C
#include <energy.h>
#include <suscpt.h>
#include <inforb.h>
#include <infpar.h>
C
      IF (IPRINT .GT. 3) CALL TITLER('Output from HER_SDRES','*',103)
C
      NUMGRA = MXCOOR
      NUMHES = MXCOOR*MXCOOR
      NUMSUS = 9
C
#if defined (VAR_MPI)
      CALL MPIXSEND(MYNUM,1,'INTEGER',MASTER,MTAG7)
C
      IF (ITYPE .EQ. 2) THEN
         CALL MPIXSEND(GRADEE,NUMGRA,'DOUBLE',MASTER,MTAG7)
         CALL MPIXSEND(HESSEE,NUMHES,'DOUBLE',MASTER,MTAG7)
      ELSE IF (ITYPE .EQ. 3) THEN
#if defined(_CRAYT3E) && !defined(NO_BINSUM)
         CALL BIN_SUM(NDMAT*N2BASX, FMAT)
         CALL BIN_SUM(8, DINTSKP)
#else
         CALL MPIXSEND(FMAT,NDMAT*N2BASX,'DOUBLE',MASTER,MTAG7)
         CALL MPIXSEND(DINTSKP,8,'DOUBLE',MASTER,MTAG7)
#endif
      ELSE IF (ITYPE .EQ. -5 .OR. ITYPE .EQ. 8) THEN
         CALL MPIXSEND(FMAT,NFMAT*N2BASX,'DOUBLE',MASTER,MTAG7)
         IF (ITYPE .EQ. -5) THEN
            CALL MPIXSEND(SUS2EL,NUMSUS ,'DOUBLE',MASTER,MTAG7)
         END IF
      END IF
C
      CALL MPIXSEND(ICOUNT,1,'INTEGER',MASTER,MTAG7)
C
      IF (TIMING) THEN
         CALL MPIXSEND(IWHICH,ICOUNT,'INTEGER',MASTER,MTAG7)
         CALL MPIXSEND(TSKCPU,ICOUNT,'DOUBLE' ,MASTER,MTAG7)
      END IF
#endif
#if defined (VAR_PVM)
      CALL PVMXINITSEND(NCODE)
      IF (ITYPE .EQ. 2) THEN
         CALL PVMXPACK(GRADEE,NUMGRA,'DOUBLE')
         CALL PVMXPACK(HESSEE,NUMHES,'DOUBLE')
      ELSE IF (ITYPE .EQ. 3) THEN
         CALL PVMXPACK(FMAT,NDMAT*N2BASX,'DOUBLE')
         CALL PVMXPACK(DINTSKP,8,'DOUBLE')
      ELSE IF (ITYPE .EQ. -5 .OR. ITYPE .EQ. 8) THEN
         CALL PVMXPACK(FMAT,NFMAT*N2BASX,'DOUBLE')
         IF (ITYPE .EQ. -5) THEN
            CALL PVMXPACK(SUS2EL,NUMSUS,'DOUBLE')
         END IF
      END IF
C
C     Pack the number of tasks received.
C
      CALL PVMXPACK(ICOUNT,1,'INTEGER')
C
C     Pack the integral-calculation times.
C
      IF (TIMING) THEN
         CALL PVMXPACK(IWHICH,ICOUNT,'INTEGER')
         CALL PVMXPACK(TSKCPU,ICOUNT,'DOUBLE')
      END IF
C
      CALL PVMXSEND(MASTER,MTAG7)
#endif
C
      RETURN
      END
C  /* Deck recompi */
      SUBROUTINE RECOMPI
C
#include <implicit.h>
#include <priunit.h>
C
      WRITE(LUPRI,'(//4(/5X,A))')
     &  'FOR RECOMPILATION:',
     &  'Change parameter in infpar.h and execute:',
     &  '   make depend',
     &  '   make'
C
      RETURN
      END
C  /* Deck prlinp */
      SUBROUTINE PRLINP(WORD)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <mtags.h>
#if defined (VAR_PVM)
      INCLUDE 'fpvm3.h'
#endif
      PARAMETER (NTABLE = 5)
      LOGICAL FIRST, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include <infpar.h>
      SAVE FIRST
      DATA TABLE /'.NODES ', '.DEGREE', '.PRINT ', '.ENCODE', '.DEBUG '/
      DATA FIRST /.TRUE./
C
C-------------------------
C     Initialize /PRLINP/.
C-------------------------
C
      IF (.NOT. FIRST) THEN
         IF (WORD .NE. '*END OF' .OR. WORD(1:2) .NE. '**') THEN
 969        READ (LUCMD, '(A7)') WORD
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      FIRST = .FALSE.
C
C     For infpar.h
C
      SLAVE  = .FALSE.
C     ... this is the master in a parallel run
      DEBUG  = .FALSE.
      NTASK  = 1
      IPRPAR = 0
      NDEGDI = 5
#if defined (VAR_PVM)
      NODTOT = 0
      NCODE  = PVMDATADEFAULT
#endif
C
      NEWDEF = (WORD .EQ. '*PARALL')
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
                     GO TO (1,2,3,4,5), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in PRLINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in PRLINP')
    1          CONTINUE
                  READ(LUCMD,*) NODS
#if defined (VAR_MPI)
                  IF (NODS .NE. NODTOT) THEN
                     WRITE (LUPRI,'(/A,I3,A,I3)')
     &                 'Number of nodes in input ', NODS,
     &                 ' .ne. nodes in MPI calculation ',NODTOT
                     CALL QUIT('Wrong number of nodes in input')
                  END IF
#endif
#if defined (VAR_PVM)
                  NODTOT = NODS
#endif
               GO TO 100
    2          CONTINUE
                  READ(LUCMD,*) NDEGDI
               GO TO 100
    3          CONTINUE
                  READ(LUCMD,*) IPRPAR
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,'(A7)') WORD
                  IF (WORD.EQ.' DATADE') THEN
                     NCODE = PVMDATADEFAULT
                  ELSE IF (WORD.EQ.' DATARA') THEN
                     NCODE = PVMDATARAW
                  ELSE IF (WORD.EQ.' DATAIN') THEN
                     NCODE = PVMDATAINPLACE
                  ELSE
                     WRITE(LUPRI,'(/1X,A,A,A//1X,A/15X,A/15X,A/)')
     &              '.ENCODE type ',word,' unknown!!!',
     &              'Options are : DATADEFAULT','DATARAW','DATAINPLACE'
                     CALL QUIT('Nonexistent encoding type in input.')
                  END IF
               GO TO 100
    5          CONTINUE
                  DEBUG = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in PRLINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in PRLINP')
            END IF
      END IF
  300 CONTINUE
C
      IF (NODTOT .GT. MAXNOD) THEN
         WRITE(LUPRI,'(/5X,A,/12X,A,2(/5X,A20,I8))')
     &    'ERROR: Number of nodes allocated has exceeded hardcoded',
     &    'limit. Reduce number of nodes or recompile program. ',
     &    'Allocated (NODTOT) :',NODTOT,'Limit (MAXNOD)     :',MAXNOD
         CALL RECOMPI
         CALL QUIT('ERROR: Number of nodes requested is too large.')        
      ENDIF
C
      RETURN
      END
#if defined(_CRAYT3E) && !defined(NO_BINSUM)
      SUBROUTINE BIN_SUM(N, SUM)
C     ******************************************************************
C     *  
C     *  Performs a binary tree sum of SUM(0:N-1) across all PE's and
C     *  leaves the result in SUM at node 0. This is equivalent to
C     *
C     *     MPI_Reduce(SUM, SUM, N, MPI_DOUBLE_PRECISION, MPI_SUM, 0,
C     *                MPI_COMM_WORLD, IERR)
C     *
C     * Input:
C     *  N      - dimension of input array SUM.
C     *  SUM    - REAL*8 array to be summed across nodes.
C     *
C     * Output:
C     *  SUM   - sum of SUM() across all nodes is returned on PE 0,
C     *          the arrays on other PEs are not modified.
C     *
C     * NOTES:
C     *        - stream safe SHMEM implementation.
C     *        - use only 1072 words of buffer space.
C     *        - calls BARRIER() which implies this is a globally
C     *          synchronizing call.
C     *        - tested with CRAY T3E only. At least some udcflushes or
C     *          some other magics will be necessary on the T3D.
C     *
C     * -BEGIN-TERMS-
C     *  Redistribution is disallowed without notifying the author. In
C     *  no event shall the lines between and including -BEGIN-TERMS-
C     *  and --END-TERMS-- be modified.
C     *
C     *  (C) 1997, Jorn Amundsen, University Administration, NTNU.
C     * --END-TERMS--
C     *
C     ******************************************************************

      INTEGER N
      REAL*8 SUM(0:N)

C---  Cray SHMEM routines
      INTEGER SHMEM_MY_PE, SHMEM_N_PES, SHMEM_INT8_CSWAP

      INTEGER NS, NR, NPROC, ME
      LOGICAL DOSUM

C---  Communication control data (NPAD is for streams padding)
      INTEGER NSEG, NPAD, M_RD, M_WRT
      PARAMETER (NSEG = 512, NPAD = 24, M_RD = 0, M_WRT = -1)
      INTEGER ISWP, ISEG, I, J
CDIR $SUPPRESS STAT
      INTEGER STAT(2)
      REAL TMP(0:2*(NSEG+NPAD))
CDIR$ CACHE_ALIGN /BINSUM__INTERNAL/
      COMMON /BINSUM__INTERNAL/ TMP,STAT



C---  Initialization
      ME = SHMEM_MY_PE()
      NPROC = SHMEM_N_PES()
      IF (NPROC .LE. 1) RETURN
      DOSUM = .TRUE.
      NS = 1


      DO WHILE (NS .LT. NPROC)

         NR = NS * 2
         CALL BARRIER() ! if DOSUM is false do barriers only

         IF (DOSUM) THEN
            IF (MOD(ME,NR) .EQ. 0 .AND. ME+NS .LT. NPROC) THEN
C---           Receiver: add up SUM and TMP in size NSEG chunks
               ISWP = 0
               DO ISEG = 0, N-1, NSEG
C                 Busy wait on STAT(ISWP) == M_WRT
 100              CONTINUE
                  IF (SHMEM_INT8_CSWAP(STAT(ISWP), M_WRT, M_WRT, ME)
     $                 .NE. M_WRT) GOTO 100

                  J = ISWP*(NSEG+NPAD)
CDIR$             UNROLL 8
                  DO I = ISEG, MIN(ISEG+NSEG-1,N-1)
                     SUM(I) = SUM(I) + TMP(J)
                     J = J + 1
                  ENDDO

C                 Update the status word
                  CALL SHMEM_INT8_SWAP(STAT(ISWP), M_RD, ME)
                  ISWP = XOR(ISWP,1)
               ENDDO
            ELSE IF (MOD(ME,NR) .NE. 0) THEN
C---           Sender: write SUM in segments of size NSEG to node ME-NS
               ISWP = 0
               DO ISEG = 0, N-1, NSEG
C                 Busy wait on STAT(ISWP) == M_RD
 200              CONTINUE
                  IF (SHMEM_INT8_CSWAP(STAT(ISWP), M_RD, M_RD, ME-NS)
     $                 .NE. M_RD) GOTO 200

                  J = ISWP*(NSEG+NPAD)
                  CALL SHMEM_PUT(TMP(J), SUM(ISEG), MIN(NSEG,N-ISEG),
     $                 ME-NS)

C                 Update the status word
                  CALL SHMEM_INT8_SWAP(STAT(ISWP), M_WRT, ME-NS)
                  ISWP = XOR(ISWP,1)
               ENDDO
               DOSUM = .FALSE.
            ENDIF
         ENDIF

         NS = NS * 2
      ENDDO

      CALL BARRIER()
      END
#endif
