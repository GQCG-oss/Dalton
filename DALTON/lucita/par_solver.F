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

#ifdef VAR_MPI
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_IMAT_SC(IBLINF,SCALFAC,ICCTOS,ICWEIGHT,IBTOTW,
     &                        ICWEIGHTF,NBLOCK,IABSOLUTE_WEIGHT)
*
* Compute the connection matrix ICCTOS between the sigma and c vector 
* Store the contribution of one c-block to sigma-blocks in ICWEIGHT.
*
* written by S. Knecht  - Feb. 2007
*
#include "implicit.h"
      INTEGER*8 KSIOIO, KSBLTP, KSVST, KLSLBT, KLSLEBT, KLSI1BT, KLSIBT,
     &          KCONSPA, KCONSPB, KCBLTP, KCIOIO, KLLBT, KLLEBT, KLI1BT,
     &          KLIBT
!               for addressing of WORK
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION IBLINF(NBLOCK),SCALFAC(NBLOCK),ICCTOS(NBLOCK,NBLOCK)
      DIMENSION ICWEIGHT(NBLOCK),ICWEIGHTF(NBLOCK)
*.Definition of c and sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
      INTEGER*8 IST8NULL
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
*. Used : NSMOB
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
*
      IDUM = 0
      IST8NULL = 0
      CALL ISETVC(IBTOTW,IST8NULL,NBLOCK)
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'IMAT_S')
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Arrays giving allowed type combinations
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KSIOIO),IPRDIA)
*. Arrays for additional symmetry operation
      KSVST = 1
*. Arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
*. Arrays for partitioning of sigma
      NTTS = MXNTTS
      CALL MEMMAN(KLSLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLSLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLSI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLSIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
*. Batches  of S vector
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &     WORK(KSIOIO),ISMOST(1,ISSM),
     &     NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &     WORK(KLSI1BT),WORK(KLSIBT),0,ITTSS_ORD)

      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
      IF(IDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE
        MXEXC = 1
      END IF
* Info for this internal space
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*. connection matrices for supergroups
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB),IPRCIX)
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. Arrays giving block type
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      KSVST = 1
*. Arrays giving allowed type combinations
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL IAIBCM(ICSPC,WORK(KCIOIO))
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*.Some TTS arrays
      NTTS = MXNTTS
*. for partitioning of vector
      CALL MEMMAN(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      CALL MEMMAN(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      CALL MEMMAN(KLI1BT,NTTS  ,'ADDL  ',1,'I1BTC ')
      CALL MEMMAN(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*
*
*. Find batches of C - strings
      ITTSS_ORD = 2
      NCBATCH = 0
      CALL PART_CIV2(IDC,WORK(KCBLTP),WORK(KNSTSO(IATP)),
     &               WORK(KNSTSO(IBTP)),
     &               NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KCIOIO),
     &               ISMOST(1,ICSM),NCBATCH,
     &               WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &               0,ITTSS_ORD)
                      
*
*. Initialize SCALFAC as it should look like after a few iterations ...
      DO II = 1, NBLOCK
        LBL = IBLINF(II)
        IF(LBL.NE.0) THEN 
          SCALFAC(II) = 1.0D0
        ELSE
          SCALFAC(II) = 0.0D0
        END IF
      ENDDO
CSK      IF(LUCI_MYPROC.EQ.LUCI_MASTER) CALL WRTMT_LU(SCALFAC,1,NBLOCK,1,NBLOCK)
*. start the search for the INTERACT-MATRIX...
      CALL FIND_INTERACT_MAT(NBLOCK,SCALFAC,ICCTOS,ICWEIGHT,ICWEIGHTF,
     &                      IDC,WORK(KSVST),WORK(KCONSPA),WORK(KCONSPB),
     &                      WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &                      NOCTPA,NOCTPB,WORK(KLIBT),WORK(KLSIBT),
     &                      PSSIGN,MXEXC,NSMST)
      NTEST = 100
CSK      IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
        IF(NTEST.GE.100)THEN
           IF(NTEST.GE.500.AND.LUCI_MYPROC.EQ.LUCI_MASTER) THEN
             WRITE(6,*)'CONNECTION MATRIX BETWEEN C- AND SIGMA-BLOCKS'
             CALL IWRTMAS(ICCTOS,NBLOCK,NBLOCK,NBLOCK,NBLOCK)
           END IF
          CALL PRINT_BATCH_INFO(NCBATCH,WORK(KLLBT),WORK(KLLEBT),
     &         WORK(KLI1BT),WORK(KLIBT),ICCTOS,ICWEIGHT,NBLOCK,IBTOTW,
     &         ICWEIGHTF,IABSOLUTE_WEIGHT)
        END IF
CSK        IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
CSK          WRITE(6,*)'CONNECTION FACTOR MATRIX'
CSK          CALL IWRTMA(ICWEIGHTF,1,NBLOCK,1,NBLOCK)
CSK        END IF
CSK      ENDIF
*. Set SCALFAC back to default 
      CALL DZERO(SCALFAC,NBLOCK)
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'IMAT_S')
            
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_INTERACT_MAT(NBLOCK,SCALFAC,ICCTOS,ICWEIGHT,
     &           ICWEIGHTF,IDC,ISTRFL,ICONSPA,ICONSPB,NSSOA,NSSOB,
     &           NOCTPA,NOCTPB,ICBLOCK,ISBLOCK,PS,MXEXC,NSMST)
*
* written by S. Knecht  - Feb. 2007
*
      IMPLICIT REAL*8           (A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCALFAC(NBLOCK),ICCTOS(NBLOCK,NBLOCK)
      DIMENSION ICWEIGHT(NBLOCK),ICWEIGHTF(NBLOCK)
      DIMENSION ICONSPA(NOCTPA,NOCTPA),ICONSPB(NOCTPB,NOCTPB)
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION ICBLOCK(8,*),ISBLOCK(8,*)
      INTEGER NSSOA(NSMST ,*),NSSOB(NSMST ,*)
      INTEGER NTEST
      NTEST = 0
      
      DO JBLOCK = 1, NBLOCK
C
        INTERACT = 0
        JJJ = 0
C
        IF(SCALFAC(JBLOCK).EQ. 1.0D0) THEN
          JATP = ICBLOCK(1,JBLOCK)
          JBTP = ICBLOCK(2,JBLOCK)
          JASM = ICBLOCK(3,JBLOCK)
          JBSM = ICBLOCK(4,JBLOCK)
          JOFF = ICBLOCK(5,JBLOCK)
          CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
          DO IPERM = 1, NPERM
            LLASM = LASM(IPERM)
            LLBSM = LBSM(IPERM)
            LLATP = LATP(IPERM)
            LLBTP = LBTP(IPERM)
C           loop over Sigma blocks in batch
            DO JSBLOCK = 1, NBLOCK
              IF(ISBLOCK(1,JSBLOCK).GT.0) THEN
                INTERACT = 0
                IATP = ISBLOCK(1,JSBLOCK)
                IBTP = ISBLOCK(2,JSBLOCK)
                IASM = ISBLOCK(3,JSBLOCK)
                IBSM = ISBLOCK(4,JSBLOCK)
                NIA = NSSOA(IASM,IATP)
                NIB = NSSOB(IBSM,IBTP)
                IF(NIA*NIB.NE.0) THEN        
C                 are the two blocks connected by allowed excitation?
                  IF(MXEXC.EQ.2) THEN
                    IF((ICONSPA(IATP,LLATP).LE.1.AND.
     &                  ICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &                 (ICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &                  IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &                 (ICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &                  IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                           INTERACT = 1
                    END IF
                  END IF
                  IF(INTERACT.NE.0) THEN
                   JJJ = JJJ + 1
                   ICCTOS(JSBLOCK,JBLOCK) = 1
                   ICWEIGHT(JJJ) = JSBLOCK
                  END IF
                END IF
              END IF
            END DO
C           ^ loop over s-blocks
          END DO
C         ^loop over nperm
        END IF
        IF(LUCI_MYPROC.EQ.LUCI_MASTER.AND.JJJ.GT.0.AND.NTEST.GE.100)THEN
          WRITE(6,'(A,1X,I5,1X,A,1X,I5,1X,A)') 
     & 'C BLOCK',JBLOCK,'IS CONNECTING TO',JJJ,'S-BLOCKS:'
          CALL IWRTMA(ICWEIGHT,1,JJJ,1,JJJ)
        END IF
      END DO 
C     ^loop over c-blocks
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE GATHER_LOW_PAR(NSUB,NSUBMX,SUBVAL,ISCAT,
     &                          RECVARRAY,TESTARR,
     &                          IRECVARRAY,ITESTARRAY,IRANKARR,
     &                          NTESTG)
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SUBVAL(*),ISCAT(*),RECVARRAY(*),TESTARR(*)
      DIMENSION IRECVARRAY(*),ITESTARRAY(*),IRANKARR(NSUBMX)
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 0
      IF(NTEST.GE.100) THEN
        WRITE(6,*)'SUBROUTINE GATHER_LOW_PAR entered'
      END IF        
      VMAX = 1.0D99
      VMIN = -1.0D99
      IF(NTEST.GE.100) 
     & WRITE(6,*)'SET VMAX TO 1.0D99 and VMIN TO -1.0D99',VMAX,VMIN
      CALL SETVEC(TESTARR,VMAX,NSUBMX*LUCI_NMPROC)
      CALL SETVEC(RECVARRAY,VMAX,NSUBMX*LUCI_NMPROC)
      CALL IZERO(IRECVARRAY,NSUBMX*LUCI_NMPROC)
      CALL IZERO(ITESTARRAY,NSUBMX*LUCI_NMPROC)
      ITESTPURPOS = 0
      IF(LUCI_MYPROC.NE.LUCI_MASTER.AND.ITESTPURPOS.EQ.1)THEN
        CALL SETVEC(SUBVAL,VMAX,NSUBMX)
        CALL IZERO(ISCAT,NSUBMX)
      END IF
      ICOUNT = 1
      DO IL = 1, NSUBMX
        IRANKARR(IL) = ICOUNT
        ICOUNT = ICOUNT + 1
      ENDDO
*
*. gather lowest elements from nodes on all nodes - we need no master
      CALL MPI_ALLGATHER(SUBVAL,NSUBMX,MPI_REAL8,RECVARRAY,NSUBMX,
     &                  MPI_REAL8,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(ISCAT,NSUBMX,my_MPI_INTEGER,IRECVARRAY,NSUBMX,
     &                  my_MPI_INTEGER,MPI_COMM_WORLD,IERR)
      CALL MPI_ALLGATHER(IRANKARR,NSUBMX,my_MPI_INTEGER,ITESTARRAY,
     &                   NSUBMX,my_MPI_INTEGER,MPI_COMM_WORLD,IERR)
*
      IF(LUCI_MYPROC.EQ.LUCI_MASTER.AND.NTEST.GE.10)THEN
        WRITE(6,*)'before search: RECVARRAY and IRECVARRAY'
        CALL WRTMT_LU(RECVARRAY,
     &     1,NSUBMX*LUCI_NMPROC,1,NSUBMX*LUCI_NMPROC)
        CALL IWRTMA(IRECVARRAY,
     &     1,NSUBMX*LUCI_NMPROC,1,NSUBMX*LUCI_NMPROC)
        WRITE(6,*)'before search: ITESTARRAY'
        CALL IWRTMA(ITESTARRAY,
     &     1,NSUBMX*LUCI_NMPROC,1,NSUBMX*LUCI_NMPROC)
      END IF
                    
*. we need the NSUBMX lowest elements from all over the world - sorted
      ISEARCH = 0
      VMINLAST = VMIN
1000  CONTINUE
*
      VMINTMP = VMAX
      IMINRANK = 0
      IMINPLACE = NSUBMX + 1
      
      DO 300 II = 1,NSUBMX*LUCI_NMPROC
*
         TEMPMIN   = RECVARRAY(II)
         ITEMPMIN  = IRECVARRAY(II)
         ITEMPRANK = ITESTARRAY(II) 
         ITEMPN    = II
*
         IF(ITEMPMIN.GT.0)THEN
           IF(TEMPMIN.LE.VMINTMP.AND.TEMPMIN.NE.VMAX)THEN
             IF(ITEMPRANK.EQ.1.AND.IMINRANK.EQ.1)THEN
               IF(TEMPMIN.LT.VMINTMP)THEN
*               
                 IMINPLACE = ITEMPMIN
                 VMINTMP   = TEMPMIN
                 IMINNUMB  = ITEMPN
                 IMINRANK  = ITEMPRANK
*                 
               ELSE IF(TEMPMIN.EQ.VMINTMP)THEN
                 IF(ITEMPMIN.LT.IMINPLACE)THEN
*                 
                 IMINPLACE = ITEMPMIN
                 VMINTMP   = TEMPMIN
                 IMINNUMB  = ITEMPN
                 IMINRANK  = ITEMPRANK
*                 
                 END IF
               ELSE 
                 GOTO 300
               END IF
             ELSE IF(ITEMPRANK.EQ.1.AND.IMINRANK.NE.1)THEN
*                 
               IMINPLACE = ITEMPMIN
               VMINTMP   = TEMPMIN
               IMINNUMB  = ITEMPN
               IMINRANK  = ITEMPRANK
*                 
             ELSE IF(ITEMPRANK.NE.1.AND.IMINRANK.EQ.1)THEN
               IF(TEMPMIN.LT.VMINTMP) THEN
*                 
                 IMINPLACE = ITEMPMIN
                 VMINTMP   = TEMPMIN
                 IMINNUMB  = ITEMPN
                 IMINRANK  = ITEMPRANK
*                 
               END IF
             ELSE IF(ITEMPRANK.NE.1.AND.IMINRANK.NE.1)THEN
               IF(ITEMPRANK.EQ.IMINRANK) THEN
                 IF(TEMPMIN.LT.VMINTMP) THEN
*                 
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*                 
                 ELSE IF(TEMPMIN.EQ.VMINTMP)THEN
                   IF(ITEMPMIN.GT.IMINPLACE)THEN
*                 
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*                 
                   END IF
                 ENDIF
               ELSE IF(ITEMPRANK.LT.IMINRANK) THEN
                 IF(TEMPMIN.LE.VMINTMP) THEN
*                 
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*                 
                 END IF
               ELSE IF(ITEMPRANK.GT.IMINRANK) THEN
                 IF(TEMPMIN.LT.VMINTMP) THEN
*                 
                   IMINPLACE = ITEMPMIN
                   VMINTMP   = TEMPMIN
                   IMINNUMB  = ITEMPN
                   IMINRANK  = ITEMPRANK
*                 
                 END IF                 
               END IF
             END IF
           END IF
         END IF
*
*
 300  CONTINUE
*
      ISEARCH = ISEARCH + 1
*.    VMINTMP should be the lowest value w.r.t. ISEARCH 
      VMINLAST = VMINTMP
*. test writing
      IF(LUCI_MYPROC.EQ.LUCI_MASTER.AND.NTEST.GE.10) 
     & WRITE(6,*)'VMINLAST',VMINLAST
      IF(LUCI_MYPROC.EQ.LUCI_MASTER.AND.NTEST.GE.10) 
     & WRITE(6,*)'PLACE FOR VMINLAST',IMINPLACE
      IF(LUCI_MYPROC.EQ.LUCI_MASTER.AND.NTEST.GE.10) 
     & WRITE(6,*)'NUMBER',IMINNUMB
      IF(LUCI_MYPROC.EQ.LUCI_MASTER.AND.NTEST.GE.10) 
     & WRITE(6,*)'IMINRANK',IMINRANK
*. end of test writing
*. put the result to permanent storage in memory
      SUBVAL(ISEARCH) = VMINTMP
      ISCAT(ISEARCH) = IMINPLACE 
      RECVARRAY(IMINNUMB) = VMAX
      IRECVARRAY(IMINNUMB) = 0
      ITESTARRAY(IMINNUMB) = -1

      IF(ISEARCH.LT.NSUBMX) GOTO 1000
1001  CONTINUE
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MICDV6_PAR(VEC1,VEC2,RNRM,EIG,FINEIG,MAXIT,
     &                      NVAR,NROOT,MAXVEC,NINVEC,APROJ,AVEC,WORK,
     &                      IPRTXX,NPRDIM,H0,IPNTR,
     &                      NP1,NP2,NQ,H0SCR,LBLK,EIGSHF,THRES_E,
     &                      IROOTHOMING,LUWRTOUT,IBLOCKL,NBLOCKDN,
     &                      RCCTOS,LU1LIST,LU2LIST,
     &                      LU3LIST,LU4LIST,LU5LIST,LU6LIST,LU7LIST,
     &                      LUCLIST,NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                      IPROCLIST,IGROUPLIST)
*
* Iterative eigen solver, requires two blocks in core
*
* Multiroot version
*
* From MICDV6 
*
* parallel adaption, Stefan Knecht and Hans Joergen Aa. Jensen, 
* winter 2007 and improved (!) in March 2008. 
* main revision for parallel MPI file I/O 
*                          - SK (working in Odense) March 2008
*
* Input :
* =======
*        LU1 : Initial set of vectors
*        VEC1,VEC2 : Two vectors,each must be dimensioned to hold
*                    largest blocks
*        LU3,LU4   : Scratch files
*        LUDIA     : File containing diagonal of matrix
*        NROOT     : Number of eigenvectors to be obtained
*        MAXVEC    : Largest allowed number of vectors
*                    must atleast be 2 * NROOT
*        NINVEC    : Number of initial vectors ( atleast NROOT )
*        NPRDIM    : Dimension of subspace with
*                    nondiagonal preconditioning
*                    (NPRDIM = 0 indicates no such subspace )
*   For NPRDIM .gt. 0:
*          PEIGVC  : EIGENVECTORS OF MATRIX IN PRIMAR SPACE
*                    Holds preconditioner matrices
*                    PHP,PHQ,QHQ in this order !!
*          PEIGVL  : EIGENVALUES  OF MATRIX IN PRIMAR SPACE
*          IPNTR   : IPNTR(I) IS ORIGINAL ADRESS OF SUBSPACE ELEMENT I
*          NP1,NP2,NQ : Dimension of the three subspaces
*
* H0SCR : Scratch space for handling H0, at least 2*(NP1+NP2) ** 2 +
*         4 (NP1+NP2+NQ)
*           LBLK : Defines block structure of matrices
*
      use dalton_mpi
      use par_mcci_io
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
      INTEGER IERR
#include "parluci.h"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
      DIMENSION APROJ(MAXVEC*(MAXVEC+1)/2),AVEC(*),WORK(*)
      DIMENSION H0(*),IPNTR(1)
      DIMENSION H0SCR(*)
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*), LU7LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LU6LIST(*), LUCLIST(*)
      DIMENSION LBATV(*), LEBATV(*), I1BATV(*), IBATV(8,*)
      real(8), allocatable :: tmp_space0(:)
      real(8), allocatable :: tmp_space1(:)
      real(8), allocatable :: tmp_space2(:)
*
* Dimensioning required of local vectors
*     APROJ  : MAXVEC*(MAXVEC+1)/2
*     AVEC   : MAXVEC ** 2
*     WORK   : MAXVEC*(MAXVEC+1)/2
*     H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*
*     IROOTHOMING : Do roothoming, i.e. select the
*     eigenvectors in iteration n+1 as the approximations
*     with largest overlap with the previous space
*
      DIMENSION FINEIG(NROOT)
      LOGICAL CONVER
      REAL*8 INPRDD, INPROD 
      DOUBLE PRECISION tottime,endtime,starttime
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
*     characters used for timing
      CHARACTER WALLTSTEP*12,  WPART22*12
*     XJEP is used for ROOTHOMING
      real(8), allocatable ::  xjep(:)
      integer, allocatable :: ixjep(:)
      INTEGER RCCTOS(*), IGROUPLIST(LUCI_NMPROC)
      DIMENSION IBLOCKL(*),NBLOCKDN(*), IPROCLIST(LUCI_NMPROC)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      character (len=13)   :: convergence
      logical, allocatable :: rtcnv(:)
 
!     allocate scratch space
      allocate(tmp_space0(maxvec*(maxvec+1)/2))
      allocate(tmp_space1(maxvec**2))
      allocate(tmp_space2(maxvec**2))
      allocate(rtcnv(nroot))
      tmp_space0 = 0
      tmp_space1 = 0
      tmp_space2 = 0
      rtcnv      = .false.
*
C
C     set check point parameters - after every 6th iteration
C     of MAXIT
      CHECKPOINT_LUCIX = .TRUE.
      ITER_CHECKP = 0
      ICHPARAM    = 0
      IF( CHECKPOINT_LUCIX )THEN
        ICHPARAM    = 6
      END IF

      WRITE(luwrt,'(/7X,A)') 
     &'  ***************************************************'
      WRITE(luwrt,'(7X,A)') 
     &'     entering MICDV6_par (parallel solver routine)   '
      WRITE(luwrt,'(7X,A)') 
     &'  ***************************************************'
*
C
C     transfer root information to common block in parluci.h
      NROOT_INFO = NROOT
C
*
      IPICO = 0
      IPRT  = 0000
!     IPRT  = 9999 ! debug
       
      IOLSTM = 1
      IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     &WRITE(luwrt,*) ' Inverse iteration modified Davidson '
      IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     &WRITE(luwrt,*) ' Normal Davidson method '
      IF( MAXVEC .LT. 2 * NROOT ) THEN
        WRITE(luwrt,*) ' Sorry MICDV6 wounded , MAXVEC .LT. 2*NROOT '
        WRITE(luwrt,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
        WRITE(luwrt,*) ' Raise MXCIV to be at least 2 * Nroot '
        WRITE(luwrt,*) ' Enforced stop on MICDV6 '
        Call Abend1( 20 )
      END IF
*
      IF(IROOTHOMING.EQ.1) THEN
        WRITE(luwrt,*) ' Root homing performed '
        allocate(xjep(maxvec*3))
        allocate(ixjep(maxvec*3))
      END IF
      KAPROJ = 1
      KFREE  = KAPROJ+ MAXVEC*(MAXVEC+1)/2
      TEST   = 1.0D-4
      CONVER = .FALSE.
C
C     ===================
C      Initial iteration
C     ===================
C
C     start timing of initial iteration
      WALLITR1 = MPI_WTIME()

      call dzero(aproj,maxvec*(maxvec+1)/2)
C
      CALL IZERO(LU2LIST,IALL_LU2)
C
      ITER = 1
      ITER_CHECKP = 1
C
      DO 10 IVEC = 1, NINVEC
*
*       copy c-vector to working-file ILUC
*
        CALL MPI_BARRIER(MYNEW_COMM,IERR)
*
*       reset LUCLIST; MY_ACT_BLK_ALL = NUMBLOCKS!
*
        CALL IZERO(LUCLIST,NUM_BLOCKS2)
*
csk     WRITE(luwrt,*) 'LUCLIST for the 1st time:'
csk     CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,luwrt)
csk     WRITE(luwrt,*) 'LU1LIST for the 1st time:'
csk     CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,luwrt)
*
        call mcci_cp_vcd_batch(ilu1,iluc,vec1,
     &                         nbatv,lbatv,lebatv,i1batv,
     &                         ibatv,iblockl,my_lu1_off,my_luc_off,
     &                         lu1list,luclist,my_vec2_ioff,my_act_blk2,
     &                         ivec-1)
!       CALL COPVCD_PP_CC_B(ILU1,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
!    &                      IBATV,MY_LU1_OFF,MY_LUC_OFF,LU1LIST,
!    &                      LUCLIST,IBLOCKL,IVEC-1)
*
!       WRITE(luwrt,*) 'LUCLIST for the 2nd time:'
!       CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,luwrt)
!       WRITE(luwrt,*) 'LU1LIST for the 2nd time:'
!       CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,luwrt)
*
*       set offset for sigma-file
*
        JVEC_SF = IVEC - 1
*
*       start calculation: sigma = H x C
*
*       timing this sigma-vector computation
*
        sigmatime = MPI_WTIME()
*
*=======================================================================
        CALL MV7(VEC1,VEC2,ILUC,ILU2 
#ifdef VAR_MPI
     &           ,LUCLIST,LU2LIST,IBLOCKL,NBLOCKDN,RCCTOS,IGROUPLIST,
     &           IPROCLIST
#endif
     &           )
*=======================================================================
*
*       end of timing
        sigmatime2 = MPI_WTIME()
        WALLTID = SECTID(sigmatime2-sigmatime)
        if(IVEC .eq. 1)then
          WRITE(luwrt,9777) WALLTID
        else
          WRITE(luwrt,9778) WALLTID
        end if
*
*       projected matrix using batch structure of CI vector(s)
*
        CALL INPROD_B_PAR_RL(ILU1,ILU2,VEC1,VEC2,APROJ,
     &                       NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                       MY_LU1_OFF,MY_LU2_OFF,LU1LIST,
     &                       LU2LIST,IVEC)
 10   CONTINUE
*     ^ End of loop over IVEC
*
*     timing of initial iteration
*
      WALLITR2 = MPI_WTIME()
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(luwrt,9888) WALLTID
 
!     synchronize MPI_COMM_WORLD
      iredl = ninvec*(ninvec-1)/2 + ninvec
      call mpi_allreduce(mpi_in_place,aproj,iredl,mpi_real8,
     &                   mpi_sum,mpi_comm_world,ierr)
*
*
      IF( IPRT .GE.5 ) THEN
        WRITE(luwrt,*) ' INITIAL PROJECTED MATRIX  '
        CALL PRSYM(APROJ,NINVEC)
      END IF

C     Diagonalize initial projected matrix
      CALL DCOPY(NINVEC*(NINVEC+1)/2,APROJ,1,WORK(KAPROJ),1)
      CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NINVEC,0,1)
      DO IROOT = 1, NROOT
        EIG(1,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2 )
      END DO
C
      IF(IPRT .GE. 3 ) THEN
        WRITE(luwrt,*) ' Eigenvalues of initial iteration (with
     &  shift)'
        WRITE(luwrt,'(5F18.13)')
     &  ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
      END IF
      IF( IPRT  .GE. 5 ) THEN
        WRITE(luwrt,*) ' Eigenvalues of initial iteration (no shift)'
        WRITE(luwrt,'(5F18.13)')
     &  ( EIG(1,IROOT),IROOT=1,NROOT)
      END IF
C      
C     transform vectors: C and sigma
C
      CALL IZERO(LU3LIST,IALL_LU3)
C
      CALL TRAVC_B_RL_DRV(VEC1,VEC2,AVEC,LU1LIST,LU3LIST,
     &                    NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                    MY_LU1_OFF,MY_LU3_OFF,
     &                    NINVEC,NROOT,ILU1,ILU3,IALL_LU1)
C
      CALL IZERO(LU3LIST,IALL_LU3)
C
      CALL TRAVC_B_RL_DRV(VEC1,VEC2,AVEC,LU2LIST,LU3LIST,
     &                    NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                    MY_LU2_OFF,MY_LU3_OFF,
     &                    NINVEC,NROOT,ILU2,ILU3,IALL_LU2)
C
      CALL IZERO(LU3LIST,IALL_LU3)
      CALL IZERO(LU4LIST,IALL_LU4)
C
C     and the corresponding Hamiltonian matrix, no problems
C     with numerical stabilities, so just use eigenvalues
C
      CALL DZERO(APROJ,NROOT*(NROOT+1)/2)
C
      DO IROOT = 1, NROOT
        APROJ(IROOT*(IROOT+1)/2) = EIG(1,IROOT)
      END DO
C
C
      NVEC = NROOT
C
      IF (MAXIT .EQ. 1 ) GOTO  1001
*
* ======================
*. Loop over iterations
* ======================
*
 1000 CONTINUE
*     start timing of iteration
      WALLITR1 = MPI_WTIME()
      starttime = MPI_WTIME()
*
      write(luwrt,'(/A21,3X,I3)') ' Info from iteration ',ITER
      write(luwrt,*) '_______________________'
      ITER = ITER + 1
      ITER_CHECKP = ITER_CHECKP + 1
*
*
*===========================================================
*                       PART 1                             =
*                                                          =
*              New directions to be included               =
*                                                          =
*===========================================================
*
*     1.1 : R = H*X - EIGAPR*X
*
      IADD = 0
      CONVER = .TRUE.
*
      DO 100 IROOT = 1, NROOT
*
*       reset scratch file lists ...
*
        CALL IZERO(LU5LIST,IALL_LU5)
        CALL IZERO(LU6LIST,IALL_LU6)
        CALL IZERO(LU7LIST,IALL_LU7)
*
        EIGAPR = EIG(ITER-1,IROOT)
*
*       calculate residues ...
*
        CALL P1_B_PAR_RL_LUCI1(VEC1,VEC2,EIGAPR,RNRM,EIGSHF,
     &                         EIG,TEST,THRES_E,RTCNV,CONVER,ITER,MAXIT,
     &                         IROOT,LU2LIST,LU1LIST,LU5LIST,NBATV,
     &                         LBATV,LEBATV,I1BATV,IBATV,
     &                         MY_LU2_OFF,MY_LU1_OFF,MY_LU5_OFF,
     &                         tmp_space0,ILU2,ILU1,ILU5)
*
*
        IF( ITER .GT. MAXIT) GOTO 100
*
*       new direction needed?
*
*       1.2 : multiply with inverse Hessian approximation
*             to get new direction
*
        IF( .NOT. RTCNV(IROOT) ) THEN
*
*         (D-E)-1 *( HX - EX )
*
          IADD = IADD + 1
*
          CSCREEN = .TRUE.
csk       CSCREEN = .FALSE.
*
          CALL H0M1TD_REL_PAR(ILU6,IDIA,ILU5,-EIGAPR,VEC1,VEC2,
     &                        LU6LIST,LU5LIST,NBATV,LBATV,LEBATV,
     &                        I1BATV,IBATV,MY_LU6_OFF,MY_DIA_OFF,
     &                        MY_LU5_OFF,1,THRES_E)
*              H0M1TD_REL_PAR(LUOUT,LUDIA,LUIN,SHIFT,VEC1,VEC2,LISTOUT,
*    &                        LISTIN,NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
*    &                        OFFSET_OUT,OFFSET_DIAG,OFFSET_IN,INV)
*
*
*         hardwired to 1!
*
          IF(IOLSTM .NE. 0 ) THEN
*
*           add Olsen correction if neccessary
*
            CSCREEN = .FALSE.
*
            CALL IZERO(LU5LIST,IALL_LU5)
            CALL IZERO(LU7LIST,IALL_LU7)
*
            CALL P1_B_PAR_RL_LUCI2(VEC1,VEC2,-EIGAPR,IROOT,
     &                             LU1LIST,LU5LIST,LU7LIST,LU6LIST,
     &                             NBATV,LBATV,LEBATV,
     &                             I1BATV,IBATV,MY_LU1_OFF,MY_LU5_OFF,
     &                             MY_LU7_OFF,MY_LU6_OFF,MY_DIA_OFF,
     &                             ILU1,ILU5,ILU7,ILU6,IDIA,1)
*
          END IF
*
*         1.3 orthogonalize to all previous vectors
*         1.4 normalize vector
*
          CALL IZERO(LU5LIST,IALL_LU5)
*
          CALL P1_B_PAR_RL_LUCI3(VEC1,VEC2,WORK,LU1LIST,LU6LIST,
     &                           LU3LIST,LU5LIST,
     &                           NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                           MY_LU1_OFF,MY_LU6_OFF,MY_LU3_OFF,
     &                           MY_LU5_OFF,tmp_space0,NVEC,IADD,
     &                           ILU1,ILU6,ILU3,ILU5)
*
        END IF
*       ^ converged?
  100 CONTINUE
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      IF( TIMING_par )
     &WRITE(LUWRT,9250)WALLTID
      IF( CONVER ) THEN
         GOTO  1001
      END IF
      IF( ITER.GT. MAXIT) THEN
         ITER = MAXIT
         GOTO 1001
      END IF
*
*
*===========================================================
*                       PART 2                             =
*                                                          =
*         Optimal combination of new and old directions    =
*                                                          =
*===========================================================
*
*     2.1: multiply new directions with matrix
*
      starttime = MPI_WTIME()
      xixidletime = 0.0D0
*
      call dzero(tmp_space0,maxvec*(maxvec+1)/2)

      IMUSTRED = 0
      ISTRED   = 0
*      
      DO 150 IVEC = 1, IADD
*
*       copy c-vector to working-file ILUC
*
*
        xidletime = MPI_WTIME()
        CALL MPI_BARRIER(MYNEW_COMM,IERR)
        xixidletime = xixidletime - xidletime + MPI_WTIME()
*
csk     WRITE(LUWRT,*) 'LU3LIST again:'
csk     CALL IWRTMAMN(LU3LIST,1,IALL_LU3,1,IALL_LU3,LUWRT)
*
*       reset LUCLIST
*
        CALL IZERO(LUCLIST,MY_ACT_BLK_ALL)
*
        call mcci_cp_vcd_batch(ilu3,iluc,vec1,
     &                         nbatv,lbatv,lebatv,i1batv,
     &                         ibatv,iblockl,my_lu3_off,my_luc_off,
     &                         lu3list,luclist,my_vec2_ioff,my_act_blk2,
     &                         nvec+ivec-1-nroot)
!       CALL COPVCD_PP_CC_B(ILU3,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
!    &                      IBATV,MY_LU3_OFF,MY_LUC_OFF,LU3LIST,
!    &                      LUCLIST,IBLOCKL,NVEC+IVEC-1-NROOT)
*
!       WRITE(LUWRT,*) 'LUCLIST for the 2nd time:'
!       CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,LUWRT)
!       WRITE(LUWRT,*) 'LU3LIST for the 2nd time:'
!       CALL IWRTMAMN(LU3LIST,1,IALL_LU3,1,IALL_LU3,LUWRT)
*
*       set offset for sigma-file
*
        JVEC_SF = NVEC + IVEC - 1 - NROOT
*
*       start calculation: sigma = H x C
*
*       timing this sigma-vector computation
        sigmatime = MPI_WTIME()
*
*==================================================================
        CALL MV7(VEC1,VEC2,ILUC,ILU4
#ifdef VAR_MPI
     &           ,LUCLIST,LU4LIST,IBLOCKL,NBLOCKDN,RCCTOS,IGROUPLIST,
     &           IPROCLIST
#endif
     &            )
*==================================================================
*
*       end of timing
        sigmatime2 = MPI_WTIME()
        WALLTID = SECTID(sigmatime2-sigmatime)
        if(ivec.eq.1) then 
          WRITE(LUWRT,9400) WALLTID
        else
          WRITE(LUWRT,9401) WALLTID
        end if
*
*       augment projected matrix using batch structure of CI vector(s)
*
        CALL INPROD_B_PAR_RL_LUCI2(ILU4,ILU1,ILU3,VEC1,VEC2,tmp_space0,
     &                             NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                             MY_LU4_OFF,MY_LU1_OFF,MY_LU3_OFF,
     &                             LU4LIST,LU1LIST,LU3LIST,IVEC,NVEC,
     &                             IMUSTRED,ISTRED)
*
  150 CONTINUE
*
*     synchronize MPI_COMM_WORLD
*
*
      xidletime = MPI_WTIME()
      call dzero(aproj(istred),imustred)
      call mpi_allreduce(tmp_space0(istred),aproj(istred),imustred,
     &                   mpi_real8,mpi_sum,mpi_comm_world,ierr)
      xixidletime = xixidletime - xidletime + MPI_WTIME()
 
#ifdef LUCI_DEBUG
      write(luwrt,*)'APROJ is (MAXVEC*(MAXVEC+1)/2)',
     &                         MAXVEC*(MAXVEC+1)/2
      CALL WRTMATMN(APROJ,1,MAXVEC*(MAXVEC+1)/2,1,
     &              MAXVEC*(MAXVEC+1)/2,LUWRT)
#endif

*
*     2.2: diagonalize projected matrix
*
      NVEC = NVEC + IADD
      CALL DCOPY(NVEC*(NVEC+1)/2,APROJ,1,WORK(KAPROJ),1)
      CALL EIGEN_LUCI(WORK(KAPROJ),AVEC,NVEC,0,1)
*
      endtime      = MPI_WTIME()
      tottime      = 0.0D0
      tottime_save = 0.0D0
      tottime      = endtime - starttime
C     TIMING FOR PARTS 2.1 - 2.2 
      tottime_save = tottime
      WALLTID      = SECTID(tottime) 
      WRITE(LUWRT,9350) WALLTID
C
      IF( TIMING_par )THEN
C
        xixidletime_save = 0.0D0
        xixidletime_save = xixidletime
        WALLTID = SECTID(xixidletime)
C
C       print idle time
C
        WRITE(LUWRT,'(/A,1X,A)') 
     &  ' accumulated idle time in part 2                 :',WALLTID
        xpercent = (xixidletime_save/tottime_save) * 100
        WRITE(LUWRT,'(A,F14.9,A/)') 
     &  ' ratio (idle time)/(time part 2) =',xpercent,' %'
C
      END IF
C
C
      IF(IROOTHOMING.EQ.1) THEN
C
C      Reorder roots so the NROOT with the largest overlap with
C      the original roots become the first
C
C      Norm of wavefunction in previous space
       DO IVEC = 1, NVEC
         XJEP(IVEC) = INPROD(AVEC(1+(IVEC-1)*NROOT),
     &                AVEC(1+(IVEC-1)*NROOT),NROOT)
       END DO
       WRITE(LUWRT,*)
     & ' Norm of projections to previous vector space '
       CALL WRTMATMN(XJEP,1,NVEC,1,NVEC,LUWRT)
C      My sorter arranges in increasing order, multiply with minus 1
C      so the eigenvectors with largest overlap comes out first
       ONEM = -1.0D0
       CALL DSCAL(NVEC,ONEM,XJEP,1)
       CALL SORLOW(XJEP,XJEP(1+NVEC),IXJEP,NVEC,NVEC,NSORT,IPRT)
       IF(NSORT.LT.NVEC) THEN
         WRITE(LUWRT,*) ' Warning : Some elements lost in sorting '
         WRITE(LUWRT,*) ' NVEC,NSORT = ', NSORT,NVEC
       END IF
       IF(IPRT.GE.3) THEN
         WRITE(LUWRT,*) ' New roots choosen as vectors '
         CALL IWRTMAMN(IXJEP,1,NROOT,1,NROOT,LUWRT)
       END IF
C      Reorder
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         CALL DCOPY(NVEC,AVEC(1+(IOLD-1)*NVEC),1,
     &              XJEP(1+(INEW-1)*NVEC),1)
       END DO
       CALL DCOPY(NROOT*NVEC,XJEP,1,AVEC,1)
       DO INEW = 1, NVEC
         IOLD = IXJEP(INEW)
         XJEP(INEW*(INEW+1)/2) = WORK(IOLD*(IOLD+1)/2)
       END DO
       DO INEW = 1, NVEC
         WORK(INEW*(INEW+1)/2) = XJEP(INEW*(INEW+1)/2)
       END DO
C
       IF(IPRT.GE.3) THEN
         WRITE(LUWRT,*) ' Reordered WORK and AVEC arrays '
         CALL PRSYM(WORK,NVEC)
         CALL WRTMATMN(AVEC,NVEC,NVEC,NVEC,NVEC,LUWRT)
       END IF
C
      END IF
C     ^ End of root homing procedure
C
      DO IROOT = 1, NROOT
        EIG(ITER,IROOT) = WORK(KAPROJ-1+IROOT*(IROOT+1)/2)
      END DO
C
csk   IPRT = 5
C
      IF(IPRT .GE. 3 ) THEN
        WRITE(LUWRT,'(A,I4)') ' Eigenvalues of iteration ..', ITER
        WRITE(LUWRT,'(5F18.13)')
     &  ( EIG(ITER,IROOT)+EIGSHF,IROOT=1,NROOT)
        WRITE(LUWRT,'(A)') ' Norm of Residuals (Previous it) '
        WRITE(LUWRT,'(5F18.13)')
     &  ( RNRM(ITER-1,IROOT),IROOT=1,NROOT)
      END IF
C
      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUWRT,*) ' Projected matrix and eigen pairs '
        CALL PRSYM(APROJ,NVEC)
        WRITE(LUWRT,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMATMN(AVEC,NVEC,NROOT,NVEC,NROOT,LUWRT)
      END IF
*
*
*     check timing
      timer3 = 0.0D0
      starttimer = MPI_WTIME()
*
*===========================================================
*                         PART 3                           =
*                                                          =
*      perhaps reset or assemble converged eigenvectors    =
*                                                          =
*===========================================================
*
*     case 1 : Only NROOT vectors can be stored
*              save current eigenvector approximations
*     case 2 : Atleast 2*NROOT eigenvectors can be saved
*              Current eigenvactor approximations+
*              vectors allowing generation of previous approxs.
*
      CALL IZERO(LU5LIST,IALL_LU5)
*
*     c vectors to ILU1
*
      DO IROOT = 1, NROOT
*
        CALL P3_B_PAR_RL_LUCI1(VEC1,VEC2,AVEC,
     &                         LU1LIST,LU3LIST,LU5LIST,
     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                         MY_LU1_OFF,MY_LU3_OFF,MY_LU5_OFF,
     &                         NVEC,NROOT,IROOT,ILU1,ILU3,ILU5)
*
      END DO
*
*     update WORK array to get correct scaling factor
*
      CALL DZERO(WORK,NROOT)
*
*     no scaling, we should already work in a normalized basis
*
      CALL SETVEC(WORK,1.0D0,NROOT)
      CALL IZERO(LU1LIST,IALL_LU1)
*
      DO IROOT = 1, NROOT
*
        CALL COPVCD_PP_B_RL(VEC1,LU5LIST,LU1LIST,
     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                      MY_LU5_OFF,MY_LU1_OFF,IROOT,ILU5,ILU1)
*
      END DO
*
      CALL IZERO(LU5LIST,IALL_LU5)
*
*     corresponding sigma vectors to ILU2
*
      DO IROOT = 1, NROOT
*
        CALL P3_B_PAR_RL_LUCI1(VEC1,VEC2,AVEC,
     &                         LU2LIST,LU4LIST,LU5LIST,
     &                         NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                         MY_LU2_OFF,MY_LU4_OFF,MY_LU5_OFF,
     &                         NVEC,NROOT,IROOT,ILU2,ILU4,ILU5)
*
      END DO
*
*     update WORK array to get correct scaling factor
*
      CALL DZERO(WORK,NROOT)
*
*     no scaling, we should already work in a normalized basis
*
      CALL SETVEC(WORK,1.0D0,NROOT)
      CALL IZERO(LU2LIST,IALL_LU2)
*
      DO IROOT = 1, NROOT
*
        CALL COPVCD_PP_B_RL(VEC1,LU5LIST,LU2LIST,
     &                      NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                      MY_LU5_OFF,MY_LU2_OFF,IROOT,ILU5,ILU2)
*
      END DO
*
      CALL IZERO(LU5LIST,IALL_LU5)
*
      NNVEC = NROOT
*
      IF( 3*NROOT .le. MAXVEC )THEN
*
*       Orthogonalize the
*       last set of correction vectors to the current
*       eigenvectors on LU1, and save on LU2
*       Overlap with root approximations
*       Start of last set of trial vectors
*
        ISTART = NVEC-NROOT-IADD+1
*
        CALL IZERO(LU5LIST,IALL_LU5)
        CALL IZERO(LU7LIST,IALL_LU7)
*
        !****************!
        DO JVEC = 1, IADD
        !****************!
*
*
*         Orthogonalize to vectors on LU1 and to trial vectors on LU5
*
          CALL ORTHG_VEC_BATCH_LUCI1(VEC1,VEC2,WORK,tmp_space1,
     &                               tmp_space2,NROOT,JVEC,ISTART,
     &                               NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                               MY_LU3_OFF,MY_LU1_OFF,MY_LU5_OFF,
     &                               MY_LU7_OFF,LU3LIST,LU1LIST,LU5LIST,
     &                               LU7LIST,ILU3,ILU1,ILU5,ILU7)
*
*
          WORK(NROOT+JVEC+(JVEC-1)*2*NROOT) = 1.0D0
*
*         current vector to work on has been written to ILU7
*
          CALL ORTHG_VEC_BATCH_LUCI2(VEC1,VEC2,WORK(1+(JVEC-1)*2*NROOT),
     &                               WORK(NROOT+1+(JVEC-1)*2*NROOT),
     &                               NROOT,JVEC,
     &                               NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                               MY_LU1_OFF,MY_LU7_OFF,MY_LU5_OFF,
     &                               MY_LU5_OFF,MY_LU3_OFF,
     &                               WORK((JVEC-1)*2*NROOT+1),
     &                               NROOT+JVEC,IADD,
     &                               LU1LIST,LU7LIST,LU5LIST,
     &                               LU5LIST,LU3LIST,
     &                               ILU1,ILU7,ILU5,ILU5,ILU3)
*
        !*****!
        END DO
        !*****!
*       end of loop over orthogonalized directions -- IADD
*
        CALL IZERO(LU5LIST,IALL_LU5)
        CALL IZERO(LU7LIST,IALL_LU7)
*
*       sigma vectors corresponding to orthogonalized directions
*
        !****************!
        DO JVEC = 1, IADD
        !****************!
*
*
          FACT = WORK(NROOT+JVEC+(JVEC-1)*2*NROOT)
*
*         orthogonalize sigma vectors
*
          CALL ORTHG_VEC_BATCH_LUCI3(VEC1,VEC2,WORK(1+(JVEC-1)*2*NROOT),
     &                               WORK(NROOT+1+(JVEC-1)*2*NROOT),
     &                               FACT,NROOT,JVEC,
     &                               NBATV,LBATV,LEBATV,I1BATV,IBATV,
     &                               MY_LU4_OFF,MY_LU2_OFF,
     &                               MY_LU4_OFF,MY_LU4_OFF,
     &                               LU4LIST,LU2LIST,LU4LIST,LU4LIST,
     &                               ISTART,ILU4,ILU2,ILU4,ILU4)
*
        !*****!
        END DO
        !*****!
*       End of loop over orthogonalized directions
*
        NNVEC = NROOT + IADD
*
      END IF
*     ^ end if more than NROOT vectors to be reset
*
      NVEC = NNVEC
*
*===========================================================
*                                                          =
*     new elements for subspace Hamiltonian                =
*                                                          =
*===========================================================
*    
      ISTRED     = 0
      IMUSTRED   = 0
      ILOOP      = 0
      call dzero(tmp_space0,maxvec*(maxvec+1)/2)
*
      DO IVEC = 1, NVEC
*
        IF( IVEC .gt. NROOT ) ILOOP = ILOOP + 1
*
        CALL CALC_SUBSPACE_H_LUCI(VEC1,VEC2,tmp_space0,IMUSTRED,IVEC,
     &                            NROOT,
     &                            ISTRED,NBATV,LBATV,LEBATV,I1BATV,
     &                            IBATV,ILOOP,MY_LU1_OFF,MY_LU3_OFF,
     &                            MY_LU2_OFF,MY_LU4_OFF,
     &                            LU1LIST,LU3LIST,LU2LIST,LU4LIST,
     &                            ILU1,ILU3,ILU2,ILU4)
*
      END DO
*
*     communicate new subspace hamiltonian ...
*
      call dzero(aproj(istred),imustred)
      starttime = MPI_WTIME()
      call mpi_allreduce(tmp_space0,aproj(istred),imustred,
     &                   mpi_real8,mpi_sum,mpi_comm_world,ierr)
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      IF( TIMING_par ) WRITE(luwrt,9460) WALLTID
*
*     finish timing
*
      timer3 = timer3 + MPI_WTIME() - starttimer
      WALLTSTEP = SECTID(timer3)
      WRITE(luwrt,9600) WALLTSTEP
*
*     timing for this iteration
*
      WALLITR2 = MPI_WTIME()
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(luwrt,9300) WALLTID
*
*     end of resetting business
*
C
C     save current solution vectors after each 6th iteration
      IF( CHECKPOINT_LUCIX .and. (ITER_CHECKP.ge.ICHPARAM))THEN
C
C       reset ITER_CHECKP
        ITER_CHECKP = 0
C
C       copy c-vectors from nodes and master back to the master
        CALL REWINE(61,-1)
        call copy_vector_mpi_2_sequential_io(ILU1,61,VEC1,
     &                                       MY_LU1_OFF,LU1LIST,
     &                                       NBLOCKDN,IBLOCKL,
     &                                       MPI_COMM_WORLD,
     &                                       NUM_BLOCKS2,nroot,1)
        IF(luci_myproc .eq. luci_master) CALL REWINE(61,-1)
      END IF
C
      IF( ITER .LE. MAXIT .AND. .NOT. CONVER) GOTO 1000
C
 1001 CONTINUE
C
C     ( End of loop over iterations )
C
C      let's synchronize all processes!
C
 1003  DO 1004 I = 1, 2
        CALL MPI_BARRIER(ICOMM, IERR)
        TIME1 = MPI_Wtime()
        IF( IPRT .gt. 100 ) WRITE(LUWRT,*)
     &     'I AM IN SYNC PROC., LUCI_MYPROC = ', LUCI_MYPROC
        CALL MPI_BARRIER(ICOMM, IERR)
        TOTAL_TIME = MPI_Wtime() - TIME1
 1004  CONTINUE
C
       IF( IPRT .gt. 100) 
     &   WRITE(luwrt,*) 'SYNC. PROCESS HAS USED', TOTAL_TIME, 'SEC.'
C
      IF( .NOT. CONVER ) THEN
C        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(luwrt,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
C        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(luwrt,1180) ITER
 1180    FORMAT(1X,' Convergence was obtained in ',I3,' iterations')
      END IF
C
      do iroot = 1, nroot
        write(luwrt,'(/1x,a)')    '-------------------'
        write(luwrt,'( 1x,a,i6)') 'root number =',IROOT
        write(luwrt,'(1x,a/)')    '-------------------'
        do i = 1, iter
          write(luwrt,'( 1x,a9,i4,f20.10,4x,1p,e10.3)') 
     & 'iteration',i,eig(i,iroot)+eigshf,rnrm(i,iroot)
        end do
      end do

      write(LUWRT,'(/1x,a)') '**********************************'//
     &           '******************************'
      write(LUWRT,*) '  iter   root         energy          '//
     &           'residual     convergence'
      write(LUWRT,'(1x,a )') '**********************************'//
     &           '******************************'
      DO 1601 IROOT = 1, NROOT
         FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
         if(rtcnv(iroot)) then 
           convergence = ' converged'
         else
           convergence = 'NOT converged'
         end if
         WRITE(LUWRT,'(i7,i7,f20.10,4x,1p,E10.3,4x,a)')
     &         ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &         convergence
 1601 CONTINUE
C
      deallocate(rtcnv)
      deallocate(tmp_space2)
      deallocate(tmp_space1)
      deallocate(tmp_space0)

      IF(IROOTHOMING.EQ.1) THEN
        deallocate( xjep)
        deallocate(ixjep)
      END IF
C
      RETURN
 1030 FORMAT(3X,7F15.8,/,(3X,7F15.8))
 1120 FORMAT(3X,I3,7F15.8,/,(6X,7F15.8))
 9300 FORMAT(' >>>  WALL TIME FOR CURRENT ITERATION            : ',A)
 9410 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH NROOT: ',A)
 9420 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
 9430 FORMAT(' >>>  WALL TIME FOR REDUCING APROJ OF LENGTH     :
     & ',A,I4)
 9440 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
 9450 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
 9460 FORMAT(' >>>  WALL TIME FOR REDUCING APROJ OF LENGTH     : ',A)
 9400 FORMAT(/' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A)
 9401 FORMAT(' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A)
 9600 FORMAT(' >>>  WALL TIME IN STEP 3 OF CURRENT ITERATION   : ',A)
 9250 FORMAT(' >>>  WALL TIME IN STEP 1 OF CURRENT ITERATION   : ',A)
 9350 FORMAT(/' >>>  WALL TIME FOR PART 2.1 - 2.2               : ',A)
 9777 FORMAT(/' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A)
 9778 FORMAT(' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A)
 9888 FORMAT(/' >>>  WALL TIME FOR INITIAL ITERATION            : ',A/)
      END
#else /* defined VAR_MPI*/
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE PAR_DUMMY
C     dummy routine for normal compilation.
      END 
#endif /* defined VAR_MPI*/
