C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
      ZERO = 0.0D0
      CALL SETVEC(SCALFAC,ZERO,NBLOCK)
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
      SUBROUTINE GASDIAS_PAR(NAEL,IASTR,NBEL,IBSTR,
     &           NORB,DIAG,NSMST,H,XA,XB,SCR,RJ,RK,
     &           NSSOA,NSSOB,LUDIA,ECORE,
     &           PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,I12,
     &           IBLTP,NBLOCK,IBLKFO,NPARBLOCK)
*
* Calculate determinant diagonal
* Turbo-ras version
*
* Driven by IBLKFO, May 97
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, July 1995, GAS version
*
* I12 = 1 => only one-body part
*     = 2 =>      one+two-body part
*
* Parallel adaption in January 2007, Stefan Knecht
*
      IMPLICIT REAL*8           (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
C     REAL * 8  INPROD
*.General input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION H(NORB)
*. Specific input
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK),NPARBLOCK(NBLOCK)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION DIAG(*)
*
      NTEST =  00
      NTEST = MAX(NTEST,IPRNT)
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0
      ELSE
         XADD = 0.0D0
      END IF
*

      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMT_LU(H,1,NORB,1,NORB)
        WRITE(6,*) ' Core energy ', ECORE
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMT_LU(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMT_LU(RK,NORB,NORB,NTOOB,NTOOB)
        END IF
*
        WRITE(6,*) ' TTSS for Blocks '
        DO IBLOCK = 1, NBLOCK
          WRITE(6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
        END DO
*
        WRITE(6,*) ' I12 = ',I12
      END IF
*
*  Diagonal elements according to Handys formulae
*   (corrected for error)
*
*   DIAG(IDET) = HII*(NIA+NIB)
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIA*NJA
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIB*NJB
*              +         J(I,J) * NIA*NJB
*
*. K goes to J - K
      IF(I12.EQ.2)
     &CALL VECSUM(RK,RK,RJ,-1.0D0,+1.0D0,NTOOB **2)
      IDET = 0
      ITDET = 0
      IF(LUDIA.NE.0) CALL REWINO(LUDIA)
*
      DO 100 IBLK = 1, NBLOCK
*
*. write only blocks to disc that are necessary for this node
        IF(NPARBLOCK(IBLK).NE.LUCI_MYPROC) GOTO 100
*
        IATP = IBLKFO(1,IBLK)
        IBTP = IBLKFO(2,IBLK)
        IASM = IBLKFO(3,IBLK)
        IBSM = IBLKFO(4,IBLK)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
        IOFF =  1
        DO IA = 1, NSSOA(IASM,IATP)
          EAA = 0.0D0
          DO IEL = 1, NAEL
            IAEL = IASTR(IEL,IA)
            EAA = EAA + H(IAEL)
            IF(I12.EQ.2) THEN
              DO JEL = 1, NAEL
                EAA =   EAA + 0.5D0*RK(IASTR(JEL,IA),IAEL )
              END DO
            END IF
          END DO
          RJKAA(IA-IOFF+1) = EAA
        END DO
*. Obtain beta strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
        IBSTRT = 1
        IBSTOP =  NSSOB(IBSM,IBTP)
        DO IB = IBSTRT,IBSTOP
          IBREL = IB - IBSTRT + 1
*
*. Terms depending only on IB
*
          HB = 0.0D0
          RJBB = 0.0D0
          CALL SETVEC(XB,0.0D0,NORB)
*
          DO IEL = 1, NBEL
            IBEL = IBSTR(IEL,IB)
            HB = HB + H(IBEL )
*
            IF(I12.EQ.2) THEN
              DO JEL = 1, NBEL
                RJBB = RJBB + RK(IBSTR(JEL,IB),IBEL )
              END DO
*
              DO IORB = 1, NORB
                XB(IORB) = XB(IORB) + RJ(IORB,IBEL)
              END DO
            END IF
          END DO
          EB = HB + 0.5D0*RJBB + ECORE
*
          IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1
          END IF
          IASTOP = NSSOA(IASM,IATP)
*
          DO IA = IASTRT,IASTOP
            IDET = IDET + 1
            ITDET = ITDET + 1
            X = EB + RJKAA(IA-IOFF+1)
            DO IEL = 1, NAEL
              X = X +XB(IASTR(IEL,IA))
            END DO
            DIAG(IDET) = X
            IF(IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
          END DO
*         ^ End of loop over alpha strings|
        END DO
*       ^ End of loop over betastrings
*. Yet a RAS block of the diagonal has been constructed
        IF(ICISTR.GE.2) THEN
          IF(NTEST.GE.100) THEN
            if(myproc.eq.master)then
            write(6,*) ' number of diagonal elements to disc ',IDET
            CALL WRTMT_LU(DIAG,1,IDET,1,IDET)
            endif
          END IF
          CALL ITODS(IDET,1,-1,LUDIA)
          CALL TODSC_LUCI(DIAG,IDET,-1,LUDIA)
          IDET = 0
        END IF
 100  CONTINUE
*        ^ End of loop over blocks

      IF(NTEST.GE.5) WRITE(6,*)
     &' Number of diagonal elements generated ',ITDET
CSK      WRITE(6,*) ' Number of diagonal elements generated by node',ITDET,myproc
*
      IF(NTEST .GE.100 .AND.ICISTR.LE.1 ) THEN
        WRITE(6,*) ' CIDIAGONAL '
        CALL WRTMT_LU(DIAG(1),1,IDET,1,IDET)
      END IF
*
      IF ( ICISTR.GE.2 ) CALL ITODS(-1,1,-1,LUDIA)
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE GASDIAT_PAR(DIAG,LUDIA,ECORE,ICISTR,I12,
     &                      IBLTP,NBLOCK,IBLKFO,NPARBLOCK)
*
* CI diagonal in SD basis for state with symmetry ISM in internal
* space ISPC
*
* GAS version, Winter of 95
*
* Driven by table of TTS blocks, May97
*
      IMPLICIT REAL*8(A-H,O-Z)
* =====
*.Input
* =====
*
*
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
*
      DIMENSION IBLTP(*)
      DIMENSION IBLKFO(8,NBLOCK),NPARBLOCK(NBLOCK)
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
*
      CALL QENTER('CIDIA')
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRDIA)
*
** Specifications of internal space
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' GASDIAT_PAR speaking'
        WRITE(6,*) ' ===================='
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
      END IF
*
**. Local memory
*
      IDUM = 0
      CALL MEMMAN(IDUM,  IDUM,    'MARK  ',IDUM,'GASDIA')
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')
      CALL MEMMAN(KLXA  ,NACOB,   'ADDL  ',2,'KLXA  ')
      CALL MEMMAN(KLXB  ,NACOB,   'ADDL  ',2,'KLXB  ')
      CALL MEMMAN(KLSCR ,2*NACOB, 'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KLH1D ,NACOB,   'ADDL  ',2,'KLH1D ')
*. Space for blocks of strings
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
*
      MAXA = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      CALL GT1DIA(WORK(KLH1D))
      CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
      IF( LUDIA .GT. 0 ) CALL REWINO(LUDIA)
      CALL GASDIAS_PAR(NAEL,WORK(KLASTR),NBEL,WORK(KLBSTR),
     &     NACOB,DIAG,NSMST,WORK(KLH1D),
     &     WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &     WORK(KLK),WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &     LUDIA,ECORE,PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &     WORK(KLRJKA),I12,IBLTP,NBLOCK,IBLKFO,NPARBLOCK)
*.Flush local memory
      CALL MEMMAN(IDUM,  IDUM,    'FLUSM ',IDUM,'GASDIA')
      CALL QEXIT('CIDIA')
*
C?    Call Abend2( ' Stefan forced me to stop after GASDIA ' )
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
      SUBROUTINE H0LNSL_PAR(PHP,PHQ,QHQ,NP1DM,NP2DM,NQDM,
     &           X,RHS,S,SCR,NTESTG)
*
* Matrix H0 of the form
*
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*             *    *                  * *
*             ***************************
*
* Solve the set of equations
*
*     ( H0+S ) X = RHS

*
* =========================
* Jeppe Olsen , May 1 1990
* =========================
*
* Modified to allow solution by conjugate gradient, March 1993
* =====
* Input
* =====
* PHP : The matrix in the P1+P2 space, given in lower
*       Triangular form
* PHQ : PHQ block of matrix
* QHQ : Diagonal approximation in Q-Q space
* NP1DM : Dimension of P1 space
* NP2DM : Dimension of P2 space
* NQDM  : Dimension of Q space
* RHS   : Right hand side of equations
*
* ======
* Output
* ======
* X : solution to linear equations
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      LOGICAL CONVER
* Input
      DIMENSION PHP(*),PHQ(*),QHQ(*),RHS(*)
* Output
      DIMENSION X(*)
* Scratch
      DIMENSION SCR(*), ERROR(20+1)
*.SCR Should atleast be dimensioned 2 *(NP1DM+NP2DM)** 2 + 2 NPQDM
      DOUBLE PRECISION INPROD
      COMMON/SHFT/SHIFT
*
      EXTERNAL HPQTVM
*.
* The Q-space can be partitioned into the P -space
* to give the effective linear equation
*
* (PHP+S - PHQ  (QHQ+S)**-1 QHP ) XP = RHSP - HPQ(QHQ+S)-1 RHSQ
*
* This leads to a simple iterative scheme
*
      CALL QENTER('H0LNS')
      NTESTL =  00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST .GE. 5 ) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' H0LNSL speaking '
        WRITE(6,*) ' =============== '
      END IF
*
      NPDM = NP1DM + NP2DM
      NPQDM = NPDM + NQDM
      IROUTE = 2
*
      IF( IROUTE.EQ.1. OR. IROUTE. EQ.3 ) THEN
*. Solve by partitioning theory
*. A bit of memory
*
      KLFREE = 1
*. Space for two local P-P matrix
      KLPP1 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*

      KLPP2 = KLFREE
      KLFREE = KLFREE + NPDM ** 2
*. Two vectors in space
      KLV1 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
      KLV2 = KLFREE
      KLFREE = KLFREE + NPDM + NQDM
* =========================
*  RHSP - HPQ(QHQ+S)-1 RHSQ
* =========================
*          DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
      CALL DIAVC3_LUCI(SCR(KLV1),RHS(1+NPDM),QHQ,S,NQDM,XDUMMY)
      CALL MATML4(SCR(KLV2),PHQ,SCR(KLV1),NP1DM,1,NP1DM,NQDM,NQDM,1,0)
      CALL VECSUM(SCR(KLV1),RHS,SCR(KLV2),1.0D0,-1.0D0,NP1DM)
      CALL COPVEC(RHS(1+NP1DM),SCR(KLV1+NP1DM),NP2DM)
* ===============================
* (PHP+S - PHQ  (QHQ+S)**-1 QHP )
* ===============================
C          XDIXT2(XDX,X,DIA,NXRDM,NXCDM,SHIFT,SCR)
      CALL XDIXT2(SCR(KLPP1),PHQ,QHQ,NP1DM,NQDM,S,SCR(KLV2))
C                TRIPAK(AUTPAK,APAK,IWAY,MATDIM,NDIM,SIGNTP)
          CALL SETVEC(SCR(KLPP2),0.0D0,NPDM*(NPDM+1)/2)
          SIGNTP = 1.0
          CALL TRIPAK_LUCI(SCR(KLPP1),SCR(KLPP2),1,NP1DM,NP1DM,SIGNTP)
          CALL VECSUM(SCR(KLPP1),SCR(KLPP2),PHP,-1.0D0,1.0D0,
     &                NPDM*(NPDM+1)/2)
C                ADDDIA(A,FACTOR,NDIM,IPACK)
          CALL ADDDIA(SCR(KLPP1),S,NPDM,1)
*. Pack to full matrix
          CALL TRIPAK_LUCI(SCR(KLPP2),SCR(KLPP1),2,NPDM,NPDM,SIGNTP)
          IF(NTEST.GE.5) THEN
            WRITE(6,*) ' Partitioned matrix '
            CALL WRTMT_LU(SCR(KLPP2),NPDM,NPDM,NPDM,NPDM)
          END IF
*.Solve p equations by inverting and multiplying
C               INVMAT(A,B,MATDIM,NDIM)
           CALL INVMAT(SCR(KLPP2),SCR(KLPP1),NPDM,NPDM)
C            MATVCB(MATRIX,VECIN,VECOUT,MATDIM,NDIM,ITRNSP)
           CALL MATVCB(SCR(KLPP2),SCR(KLV1),X,NPDM,NPDM,0)
*. q part of solution
* ==================================
* XQ = (QHQ+SHIFT)**-1 (RHS Q - QHP XP)
* ==================================

         CALL MATML4(SCR(KLV1),PHQ,X,
     &        NQDM,1,NP1DM,NQDM,NP1DM,1,1)
         CALL VECSUM(SCR(KLV2),RHS(NPDM+1),SCR(KLV1),1.0D0,
     &               -1.0D0,NQDM)
*
C             DIAVC3_LUCI(VECOUT,VECIN,DIAG,SHIFT,NDIM,VDSV)
         CALL DIAVC3_LUCI(X(NPDM+1),SCR(KLV2),QHQ,
     &               S,NQDM,XDUMMY)
*
         IF(NTEST.GE.2) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMT_LU(X,1,NPQDM,1,NPQDM)
         END IF
      END IF
*
      IF (IROUTE. EQ. 2 .OR. IROUTE .EQ. 3 ) THEN
*. Use preconditioned conjugate gradient
        LU1 = 34
        LU2 = 35
        LU3 = 36
        LUDIA = 37
*
        KLV1 = 1
        KLFREE = KLV1 + NPQDM
        KLV2 = KLFREE
        KLFREE = KLFREE + NPQDM
*. Diagonal
        CALL XTRCDI(PHP,SCR(KLV1),NPDM ,1)
        CALL COPVEC(QHQ,SCR(KLV1+NPDM),NQDM)
        CALL REWINE(LUDIA,-1)
        CALL TODSC_LUCI(SCR(KLV1),NPQDM,-1,LUDIA)
*. Initial Guess
        CALL REWINE(LU1,-1)
        CALL SETVEC(SCR(KLV1),0.0D0,NPQDM)
        CALL TODSC_LUCI(SCR(KLV1),NPQDM,-1,LU1)
*. Right hand side
        CALL REWINE(LU2,-1)
        CALL TODSC_LUCI(RHS,NPQDM,-1,LU2)
*
        MAXIT = 20
        CONVER = .FALSE.
        TEST = 1.0D-9 * SQRT(INPROD(RHS,RHS,NPQDM))
        SHIFT = S
        ILNPRT = MAX(NTEST-10,0)
        CALL MINGCG(HPQTVM,LU1,LU2,LU3,LUDIA,SCR(KLV1),SCR(KLV2),
     &              MAXIT,CONVER,TEST,S,ERROR,NPQDM,0,ILNPRT)
        CALL REWINE(LU1,-1)
        CALL FRMDSC_LUCI(SCR(KLV1),NPQDM,-1,LU1,IMZERO,IAMPACK)
        CALL COPVEC(SCR(KLV1),X,NPQDM)
*
         IF(NTEST.GE.50) THEN
           WRITE(6,*) ' Solution to linear equations '
           CALL WRTMT_LU(X,1,NPQDM,1,NPQDM)
         END IF
*
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE H0M1TD_PAR(LUOUT,LUDIA,LUIN,LBLK,NPQDM,IPNTR,
     &                     H0,SHIFT,WORK,XH0PSX,
     &                     NP1,NP2,NQ,VEC1,VEC2,NTESTG)
*
* Calculate inverted general preconditioner matrix times vector
*
* Disc version
*
*  Vecut=  (H0 + shift )-1 Vecin
*
*  LUOUT       LUDIA        LUIN
*
*  and XH0PSX = X(T) (H0 + shift )** - 1 X
*
* Where H0 consists of a diagonal on file LUDIA
* and a block matrix of the form
*
*              P1    P2        Q
*             ***************************
*             *    *     *              *
*         P1  * Ex *  Ex *   Ex         *    Ex : exact H matrix
*             ***************************         is used in this block
*         P2  *    *     *              *
*             * Ex *  Ex *     Diag     *    Diag : Diagonal
*             ************              *           appriximation used
*             *    *      *             *
*             *    *        *           *
*             * Ex *  Diag    *         *
*         Q   *    *            *       *
*             *    *              *     *
*             *    *                *   *
*
* Note : The diagonal elements on LUDIA corresponding to
*        elements in the subspace are neglected,
*        i.e. their elements can have arbitrary value
*        without affecting the results
*
* The block matrix is defined by
* ==============================
*  NPQDM  : Total dimension of PQ subspace
*  NP1,NP2,NQ : Dimensions of the three subspaces
*  IPNTR(I) : Scatter array, gives adress of subblock element
*             I in full matrix
*             IPNTR gives first all elements in P1,
*             the all elements in P2,an finally all elements in Q
*  H0       : contains PHP,PHQ and QHQ in this order
*
* Jeppe Olsen , September 1993
*
*
*
* =====
* Input
* =====
*
* LUOUT : File to contain output vector
* LUDIA : File Containing diagonal of H0
* LUIN  : File Containing input vector
* LBLK : Defines format of files
* NPQDM,H0,NP1,NP2,NQ,IPNTR : Defines PQ subspace, see above
* SHIFT : constant ADDED to diagonal
* WORK : Scratch array , at least 2*(NP1DM+NP2DM) ** 2 + 4 NPQDM
*
* ======
* Output
* ======
*
* LUOUT : contains output vector, not rewinded
* XH0PSX  = X(T)(H0+SHIFT)**(-1)X
*
* =======
* Scratch
* =======
*
* VEC1,VEC2 : Must each be able to hold largest segment of vector
*
* ==========
* Externals: GATVEC,DIAVC2,SCAVEC,SBINTV,WRTMAT
* ==========
*
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
*
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION IPNTR(*),H0(*)
      DIMENSION WORK(*)
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.10)
     & write(6,*) ' H0M1TD_PAR , NPQDM = ', NPQDM
      IF(NPQDM.NE.0)THEN
       WRITE(6,*)'Problem in H0M1TD_PAR: NPQDM is not 0:',NPQDM
       Call Abend2('Problem in H0M1TD_PAR during parallel execution
     & detected')
      ENDIF
*
      KLFREE = 1
      KLV1 = KLFREE
      KLFREE = KLV1 + NPQDM
*
      KLV2 = KLFREE
      KLFREE = KLV2 + NPQDM
*
      KLSCR = KLFREE
*
      IF(NPQDM.NE.0) THEN
*. Obtain subspace components of input vector
C            GATVCD(LU,LBLK,NGAT,IGAT,XGAT,SEGMNT,NTESTG)
        IZERO = 0
        CALL GATVCD_PAR(LUIN,LBLK,NPQDM,IPNTR,WORK(KLV1),VEC1,
     &                 IZERO)
*. Solve linear equations in subspace
         KLPHP = 1
         KLPHQ = KLPHP + (NP1+NP2) *(NP1+NP2+1)/2
         KLQHQ = KLPHQ + NP1 * NQ
*
         CALL H0LNSL_PAR(H0(KLPHP),H0(KLPHQ),H0(KLQHQ),NP1,NP2,NQ,
     &                  WORK(KLV2),WORK(KLV1),SHIFT,WORK(KLSCR),
     &                  NTEST )
      END IF
*
*. Calculate inverse diagonal and scatter results from subspace,
*. Write to file LUOUT
C     DMTVDS(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
C    &                  ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
CSK      WRITE(6,*) 'I call DMTVDS_PAR with LUIN,LUOUT',LUIN,LUOUT,myproc
      IF(CSCREEN) THEN
C
        CALL DMTVDS_PAR2(VEC1,VEC2,LUDIA,LUIN,LUOUT,SHIFT,1,1,
     &                  IPNTR,WORK(KLV2),NPQDM,LBLK,XH0PSX)
C
      ELSE
        CALL DMTVDS_PAR(VEC1,VEC2,LUDIA,LUIN,LUOUT,SHIFT,1,1,
     &                 IPNTR,WORK(KLV2),NPQDM,LBLK,XH0PSX)
      END IF
*
      IF(NTEST.GT. 100 ) THEN
        WRITE(6,*) ' Output vector from H0M1TD_PAR '
        WRITE(6,*) ' ============================= '
*. Note : works only if result vector is first file on LUOUT
C            WRTVCD(SEGMNT,LU,IREW,LBLK)
        CALL WRTVCD(VEC1,LUOUT,1,LBLK)
        WRITE(6,*) ' Overlap between input and output vector',
     &               XH0PSX
      END IF
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
     &                      SCRRED,SCRRED2,RCCTOS,LU1LIST,LU2LIST,
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
!     not in use yet
!     use par_mcci_io
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      INTEGER STATUS_ARRAY(MPI_STATUS_SIZE,LUCI_NMPROC-1), IERR
      INTEGER IREQ(LUCI_NMPROC-1)
#include "parluci.h"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION RNRM(MAXIT,NROOT),EIG(MAXIT,NROOT)
      DIMENSION APROJ(*),AVEC(*),WORK(*)
      DIMENSION H0(*),IPNTR(1)
      DIMENSION H0SCR(*), SCRRED(*), SCRRED2(*)
      DIMENSION LU1LIST(*), LU2LIST(*), LU3LIST(*), LU7LIST(*)
      DIMENSION LU4LIST(*), LU5LIST(*), LU6LIST(*), LUCLIST(*)
      DIMENSION LBATV(*), LEBATV(*), I1BATV(*), IBATV(8,*)
*
* Dimensioning required of local vectors
*     APROJ  : MAXVEC*(MAXVEC+1)/2
*     AVEC   : MAXVEC ** 2
*     WORK   : MAXVEC*(MAXVEC+1)/2
*     H0SCR  : 2*(NP1+NP2) ** 2 +  4 * (NP1+NP2+NQ)
*     SCRRED : MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
*     SCRRED2: MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
*
*     IROOTHOMING : Do roothoming, i.e. select the
*     eigenvectors in iteration n+1 as the approximations
*     with largest overlap with the previous space
*
      DIMENSION FINEIG(1)
      LOGICAL CONVER,RTCNV(100)
      REAL*8 INPRDD, INPROD 
      DOUBLE PRECISION tottime,endtime,starttime
      CHARACTER SECTID*12, CPUTID*12, WALLTID*12
*     characters used for timing
      CHARACTER WALLTSTEP*12,  WPART22*12
*. Notice XJEP is also used for ROOTHOMING, should be allocated
* outside (for roothoming :dim = 3*MAXVEC )
      DIMENSION XJEP(10000)
      INTEGER   IXJEP(10000)
      INTEGER RCCTOS(*), IGROUPLIST(LUCI_NMPROC)
      DIMENSION IBLOCKL(*),NBLOCKDN(*), IPROCLIST(LUCI_NMPROC)
      INTEGER NDATATYPE, NREDTYPE, IONE, NZERO
      REAL*8 REDSCRVAR
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      IONE  = 1
      NZERO = 0
*
      TIMING = .FALSE.
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

      WRITE(LUCIWT,'(7X,A)') '                             '
      WRITE(LUCIWT,'(7X,A)') 
     &'  ***************************************************'
      WRITE(LUCIWT,'(7X,A)') 
     &'     entering MICDV6_par (parallel solver routine)   '
      WRITE(LUCIWT,'(7X,A)') 
     &'  ***************************************************'
      WRITE(LUCIWT,'(7X,A)') '                             '
*
C
C     transfer root information to common block in parluci.h
C
      NROOT_INFO = NROOT
C
*
*     initialize LZERO_SCRRED
*
      LZERO_SCRRED = 0
      LZERO_SCRRED = MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2 )
*
      IPICO = 0
      IPRT = 000
       
      IOLSTM = 1
      IF(IPRT.GT.1.AND.IOLSTM.NE.0)
     &WRITE(LUCIWT,*) ' Inverse iteration modified Davidson '
      IF(IPRT.GT.1.AND.IOLSTM.EQ.0)
     &WRITE(LUCIWT,*) ' Normal Davidson method '
      IF( MAXVEC .LT. 2 * NROOT ) THEN
        WRITE(6,*) ' Sorry MICDV6 wounded , MAXVEC .LT. 2*NROOT '
        WRITE(6,*) ' NROOT, MAXVEC  :',NROOT,MAXVEC
        WRITE(6,*) ' Raise MXCIV to be at least 2 * Nroot '
        WRITE(6,*) ' Enforced stop on MICDV6 '
        Call Abend1( 20 )
      END IF
*
      IF(IROOTHOMING.EQ.1) THEN
        WRITE(LUWRTOUT,*) ' Root homing performed '
      END IF
      KAPROJ = 1
      KFREE = KAPROJ+ MAXVEC*(MAXVEC+1)/2
      TEST = 1.0D-4
      CONVER = .FALSE.
C
C     ===================
C      Initial iteration
C     ===================
C
C     start timing of initial iteration
      WALLITR1 = MPI_WTIME()

      CALL DZERO(APROJ, MAXVEC*(MAXVEC+1)/2 )
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
csk     WRITE(LUCIWT,*) 'LUCLIST for the 1st time:'
csk     CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,LUCIWT)
csk     WRITE(LUCIWT,*) 'LU1LIST for the 1st time:'
csk     CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,LUCIWT)
*
!       call mcci_cp_vcd_batch(ilu1,iluc,vec1,
!    &                         nbatv,lbatv,lebatv,i1batv,
!    &                         ibatv,my_lu1_off,my_luc_off,lu1list,
!    &                         luclist,ivec-1)
        CALL COPVCD_PP_CC_B(ILU1,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
     &                      IBATV,MY_LU1_OFF,MY_LUC_OFF,LU1LIST,
     &                      LUCLIST,IBLOCKL,IVEC-1)
*
csk     WRITE(LUCIWT,*) 'LUCLIST for the 2nd time:'
csk     CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,LUCIWT)
csk     WRITE(LUCIWT,*) 'LU1LIST for the 2nd time:'
csk     CALL IWRTMAMN(LU1LIST,1,IALL_LU1,1,IALL_LU1,LUCIWT)
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
        WRITE(LUCIWT,9777) WALLTID
        WRITE(LUCIWT,*) '   '
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
      WRITE(LUCIWT,9888) WALLTID
      WRITE(LUCIWT,*) '   '
*
*     synchronize MPI_COMM_WORLD
*
      IREDL = 0
      IREDL = NINVEC*(NINVEC-1)/2 + NINVEC
      CALL DZERO(SCRRED,IREDL)
      CALL REDVEC_REL(APROJ,SCRRED,IREDL,2,MPI_SUM,MPI_COMM_WORLD,-1)
      CALL DCOPY(IREDL,SCRRED,1,APROJ,1)
*
*
csk   IPRT = 5
*
      IF( IPRT .GE.5 ) THEN
        WRITE(6,*) ' INITIAL PROJECTED MATRIX  '
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
        WRITE(LUCIWT,*) ' Eigenvalues of initial iteration (with
     &  shift)'
        WRITE(LUCIWT,'(5F18.13)')
     &  ( EIG(1,IROOT)+EIGSHF,IROOT=1,NROOT)
      END IF
      IF( IPRT  .GE. 5 ) THEN
        WRITE(LUCIWT,*) ' Initial set of eigen values (no shift) '
        CALL WRTMATMN(EIG(1,1),1,NROOT,MAXIT,NROOT,LUCIWT)
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
      write(LUCIWT,*)
      write(LUCIWT,'(A21,3X,I3)') ' Info from iteration ',ITER
      write(LUCIWT,*) '_______________________'
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
     &                         SCRRED,ILU2,ILU1,ILU5)
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
     &                           MY_LU5_OFF,SCRRED,NVEC,IADD,
     &                           ILU1,ILU6,ILU3,ILU5)
*
        END IF
*       ^ converged?
  100 CONTINUE
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      IF( TIMING )
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
      CALL DZERO(SCRRED,LZERO_SCRRED)
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
!       call mcci_cp_vcd_batch(ilu3,iluc,vec1,
!    &                         nbatv,lbatv,lebatv,i1batv,
!    &                         ibatv,my_lu3_off,my_luc_off,lu3list,
!    &                         luclist,nvec+ivec-1-nroot)
        CALL COPVCD_PP_CC_B(ILU3,ILUC,VEC1,NBATV,LBATV,LEBATV,I1BATV,
     &                      IBATV,MY_LU3_OFF,MY_LUC_OFF,LU3LIST,
     &                      LUCLIST,IBLOCKL,NVEC+IVEC-1-NROOT)
*
csk     WRITE(LUWRT,*) 'LUCLIST for the 2nd time:'
csk     CALL IWRTMAMN(LUCLIST,1,IALL_LUC,1,IALL_LUC,LUWRT)
csk     WRITE(LUWRT,*) 'LU3LIST for the 2nd time:'
csk     CALL IWRTMAMN(LU3LIST,1,IALL_LU3,1,IALL_LU3,LUWRT)
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
        WRITE(LUWRT,9400) WALLTID
        WRITE(LUWRT,*) '   '
*
*       augment projected matrix using batch structure of CI vector(s)
*
        CALL INPROD_B_PAR_RL_LUCI2(ILU4,ILU1,ILU3,VEC1,VEC2,SCRRED,
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
      CALL DZERO(APROJ(ISTRED),IMUSTRED)
      CALL REDVEC_REL(SCRRED(ISTRED),APROJ(ISTRED),IMUSTRED,2,
     &                MPI_SUM,MPI_COMM_WORLD,-1)
      xixidletime = xixidletime - xidletime + MPI_WTIME()
*
csk   write(LUWRTOUT,*)'APROJ is'
csk   CALL WRTMATMN(APROJ,1,MAXVEC*(MAXVEC+1)/2,1,
csk  &              MAXVEC*(MAXVEC+1)/2,LUWRT)

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
      WALLTID = SECTID(tottime) 
      WRITE(LUWRT,9350) WALLTID
C
      IF( TIMING )THEN
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
csk     CALL PRSYM(APROJ,NVEC)
        WRITE(LUWRT,'(2X,E13.7)') (EIG(ITER,IROOT),IROOT = 1, NROOT)
        CALL WRTMATMN(AVEC,NVEC,NROOT,NVEC,NROOT,LUWRT)
      END IF
      IPRT = 0
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
          CALL ORTHG_VEC_BATCH_LUCI1(VEC1,VEC2,WORK,SCRRED,
     &                               SCRRED2,NROOT,JVEC,ISTART,
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
      ISTRED   = 0
      IMUSTRED = 0
      ILOOP    = 0
      ISETSCRRED = 0
      ISETSCRRED = MAX( ISETSCRRED, MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
      CALL DZERO(SCRRED,ISETSCRRED)
*
      DO IVEC = 1, NVEC
*
        IF( IVEC .gt. NROOT ) ILOOP = ILOOP + 1
*
        CALL CALC_SUBSPACE_H_LUCI(VEC1,VEC2,SCRRED,IMUSTRED,IVEC,NROOT,
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
      CALL DZERO(APROJ(ISTRED),IMUSTRED)
      starttime = MPI_WTIME()
      CALL REDVEC(SCRRED,APROJ(ISTRED),IMUSTRED,2,MPI_SUM,
     &            MPI_COMM_WORLD,-1)
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      IF( TIMING )
     &WRITE(LUWRTOUT,9460)WALLTID
*
*     finish timing
*
      timer3 = timer3 + MPI_WTIME() - starttimer
      WALLTSTEP = SECTID(timer3)
      WRITE(LUWRTOUT,9600) WALLTSTEP
*
*     timing for this iteration
*
      WALLITR2 = MPI_WTIME()
      WALLTID = SECTID(WALLITR2-WALLITR1)
      WRITE(LUWRTOUT,9300) WALLTID
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
        DO JXROOT = 1, NROOT
          CALL REWINE(97,-1)
          CALL COPVCD_PAR_BDRIV_REL(ILU1,97,VEC1,NBLOCKDN,NUM_BLOCKS2,
     &                              IBLOCKL,MPI_COMM_WORLD,-1,JXROOT,
     &                              LU1LIST,MY_LU1_OFF,1)
          IF(MYPROC .eq. MASTER) THEN
             CALL REWINE(97,-1)
             CALL COPVCD(97,61,VEC1,0,-1)
          END IF
        END DO
        IF(MYPROC .eq. MASTER) CALL REWINE(61,-1)
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
     &   WRITE(LUWRTOUT,*) 'SYNC. PROCESS HAS USED', TOTAL_TIME, 'SEC.'
C
      IF( .NOT. CONVER ) THEN
C        CONVERGENCE WAS NOT OBTAINED
         IF(IPRT .GE. 2 )
     &   WRITE(LUWRTOUT,1170) MAXIT
 1170    FORMAT('0  Convergence was not obtained in ',I3,' iterations')
      ELSE
C        CONVERGENCE WAS OBTAINED
         ITER = ITER - 1
         IF (IPRT .GE. 2 )
     &   WRITE(LUWRTOUT,1180) ITER
 1180    FORMAT(1H0,' Convergence was obtained in ',I3,' iterations')
        END IF
C
      IF ( IPRT .GT. 0 ) THEN
        CALL REWINE(LU1,LBLK)
        DO 1600 IROOT = 1, NROOT
          write(LUWRTOUT,*)
          write(LUWRTOUT,*) '------------------------'
          write(LUWRTOUT,*) 'Root number  ',IROOT
          write(LUWRTOUT,*) '------------------------'
          DO I=1,ITER
            WRITE(LUWRTOUT,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
          end do
          WRITE(LUWRTOUT,*)
          WRITE(LUWRTOUT,'(A,I3)')
     &  ' Information about convergence for root... ' ,IROOT
          WRITE(LUWRTOUT,*)
     &    '============================================'
          WRITE(LUWRTOUT,*)
          FINEIG(IROOT) = EIG(ITER,IROOT)
          WRITE(LUWRTOUT,1190) FINEIG(IROOT)+EIGSHF
 1190     FORMAT(' The final approximation to eigenvalue ',F18.10)
          IF(IPRT.GE.400) THEN
            WRITE(LUWRTOUT,1200)
 1200       FORMAT(1H0,'The final approximation to eigenvector')
            CALL WRTVCD(VEC1,LU1,0,LBLK)
          END IF
 1340     FORMAT(1H ,6X,I4,8X,F20.13,2X,E12.5)
 1600   CONTINUE
      ELSE
        write(LUWRTOUT,*)
 1310   FORMAT
     &  (1H0,' Iteration point        Eigenvalue         Residual ')
        DO IROOT = 1, NROOT
          write(LUWRT,*) '------------------------'
          write(LUWRT,*) 'Root number  ',IROOT
          write(LUWRT,*) '------------------------'
          DO I=1,ITER
            WRITE(LUWRT,1340) I,EIG(I,IROOT)+EIGSHF,RNRM(I,IROOT)
          end do
        end do
        write(LUWRT,*)
        write(LUWRT,*) '**********************************'//
     &             '**************************'
        write(LUWRT,*) '   Iter  Root       Energy        '//
     &             'RESIDUAL     RESRATIO '
        write(LUWRT,*) '**********************************'//
     &             '**************************'
        write(LUWRT,*)
        DO 1601 IROOT = 1, NROOT
           FINEIG(IROOT) = EIG(ITER,IROOT)+EIGSHF
           WRITE(LUWRT,'(3X,I3,2X,I3,3X,F18.10,2X,E10.3,3X,E10.3)')
     &               ITER,IROOT,FINEIG(IROOT),RNRM(ITER,IROOT),
     &               RNRM(1,IROOT)/RNRM(ITER,IROOT)
           write(LUWRT,'(A,F18.10)') 
     &               ' Final energy ',FINEIG(IROOT)
 1601   CONTINUE
      END IF
C
      CALL FLSHFO(LUWRT)
C
      Call Add_Info('E_CI',FINEIG,NROOT)
C
CSK      WRITE(LUCIWT,'(A)' ) ' stop in MICDV6_PAR'
CSK      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CSK      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CSK      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CSK      CALL Abend2('Stefan forced me to stop in MICDV6_PAR')
*
      RETURN
 1030 FORMAT(1H0,2X,7F15.8,/,(1H ,2X,7F15.8))
 1120 FORMAT(1H0,2X,I3,7F15.8,/,(1H ,5X,7F15.8))
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
 9400 FORMAT(' >>>  WALL TIME FOR SIGMA VECTOR CALL            : ',A)
 9600 FORMAT(' >>>  WALL TIME IN STEP 3 OF CURRENT ITERATION   : ',A)
 9250 FORMAT(' >>>  WALL TIME IN STEP 1 OF CURRENT ITERATION   : ',A)
 9350 FORMAT(' >>>  WALL TIME FOR PART 2.1 - 2.2               : ',A)
 9777 FORMAT(' >>>  WALL TIME FOR INITIAL SIGMA VECTOR CALL    : ',A)
 9888 FORMAT(' >>>  WALL TIME FOR INITIAL ITERATION            : ',A)
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DIREC_OPTDRV(VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &                        NROOTS,MAXVEC,NVEC,APROJ,AVEC,WORK_SP,
     &                        LBLK,IBKINF,NBLOCKDN,IGROUPLIST,IPROCLIST,
     &                        BSCRFAC,IADD,SCRRED,LUWRTOUT,
     &                        ICOLOR)
C                                                                       
C     driver routine for the generation of a new best direction for 
C     the c-/s-vectors.
C                                                                      
C     file structure is adapted to batch scheme used also in MV7_PAR
C     parallel version: MPI-FILE I/O is used.                            
C
C     Last revision: S. Knecht           - March  2007
C                   
************************************************************************
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
#include "parluci.h"
#include "infpar.h"
#include "mpif.h"
#include "mxpdim.inc"
C
C     Definition of c and sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
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
C        Used: NSMOB
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
      DIMENSION IBKINF(NUM_BLOCKS2), BSCRFAC(NUM_BLOCKS2) 
      DIMENSION NBLOCKDN(NUM_BLOCKS2)
C     MY_VEC_IOFF == absolute length of 1 vec on a given cpu
C     IOFF_PROC   == max length of 1 vec for all cpus
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF, IOFF_PROC
C     MY_IOFF_125 == absolute offset for files ILU1,ILU2,ILU5 for a
C     given cpu
C     MY_IOFF_347 == absolute offset for files ILU3,ILU4,ILU7 for a
C     given cpu
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_125, MY_IOFF_347
      INTEGER ICOPY, IZERO, NBLK_ACT, NACT_MAX
      INTEGER MY_INT_IOFF_125, MY_INT_IOFF_347
      CHARACTER*12 SECTID, WALLTIME
      LOGICAL CHUNK_FILE
      IZERO = 0
      ICOPY = 0
C
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'D_OPTD')
C
      IATP = 1
      IBTP = 2
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
C     Offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
C
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
C     Arrays giving allowed type combinations
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KSIOIO),IPRDIA)
C
      KSVST = 1
C     Arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
C     Arrays for partitioning of sigma
      NTTS = MXNTTS
      CALL MEMMAN(KLSLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLSLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLSI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLSIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
C     initialize bscrfac and assume all cpu-blocks to be active
C     bscrfac is then used as determing array in the partition routine
      ZERO = 0.0D0
      ONE  = 1.0D0
      CALL DZERO(BSCRFAC,NUM_BLOCKS2)
      CALL SET_ACTIVE(BSCRFAC,NBLOCKDN,ONE,NUM_BLOCKS2)
C     partition the CI-vactor in batches
      ITTSS_ORD = 2
      IEXCLBLK  = 1
      IF (IODENSEPAR.eq.2.and.IAMSURE.eq.1) THEN
C         
        CALL PART_CIV_PAR2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &                    WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &                    WORK(KSIOIO),ISMOST(1,ISSM),
     &                    NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                    WORK(KLSI1BT),WORK(KLSIBT),0,ITTSS_ORD,
     &                    BSCRFAC,NUM_BLOCKS2,IEXCLBLK)
C
      ELSE 
      
        CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &                 WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &                 WORK(KSIOIO),ISMOST(1,ISSM),
     &                 NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                 WORK(KLSI1BT),WORK(KLSIBT),0,ITTSS_ORD)
      END IF
C
C
      IF( IAMSURE .eq. 1 ) THEN
C
C       more setup to follow ...
C
        MY_VEC_IOFF = 0
        IOFF_PROC   = 0
        NBLK_ACT    = 0
        NACT_MAX    = 0
C       I know that this is not smart but it is for the moment the 
C       easiest. -->  TODO: change LUCI_NMPROC --> NEWCOMM_PROC properly.
        DO IPROC = 1, LUCI_NMPROC
          DO IBLK = 1, NUM_BLOCKS2
            IF(NBLOCKDN(IBLK) .eq. IPROC-1 ) THEN
              MY_VEC_IOFF = MY_VEC_IOFF + IBKINF(IBLK)
              NBLK_ACT = NBLK_ACT + 1
            END IF
          END DO
          IOFF_PROC = MAX(IOFF_PROC,MY_VEC_IOFF)
          NACT_MAX = MAX(NACT_MAX,NBLK_ACT)
          MY_VEC_IOFF = 0
          NBLK_ACT = 0
        END DO 
C
        MY_VEC_IOFF = 0
        NBLK_ACT = 0
C
C       length of each vector + number of active blocks on a given cpu
C
        DO IBLK = 1, NUM_BLOCKS2
          IF( NBLOCKDN(IBLK) .eq. LUCI_MYPROC ) THEN 
            MY_VEC_IOFF = MY_VEC_IOFF + IBKINF(IBLK)
            NBLK_ACT = NBLK_ACT + 1
          END IF
        END DO
C
C       check number of copies on LU3/LU4
C
        IF(MAX(NVEC-NROOT,NROOTS).LE.NVEC-NROOTS)THEN
               ICOPY = MAX(NVEC-NROOTS,NROOTS)
        ELSE
               ICOPY = MIN(NVEC-NROOTS,NROOTS)
        END IF
C       check ICOPY
        IF( ICOPY .lt. 1 ) CALL Abend2('Problem in DIREC_OPTDRV since 
     & the number of vectors to copy (ICOPY) is lower than 1')
C       copy sequential files in MPI-IO readable format
        CALL REWINE(LU1,-1)
        CALL REWINE(LU2,-1)
        CALL REWINE(LU3,-1)
        CALL REWINE(LU4,-1)
C
CSK        CALL MEMMAN(KINFOACTIVE,NBLK_ACT,'ADDL  ',1,'ICPBLA')
CSK        CALL ISETVC(WORK(KINFOACTIVE),IZERO,NBLK_ACT)
C
C       calculate correct offset for each cpu (in a group)
C
        MY_IOFF_125      = IOFF_PROC * NROOTS * MYNEW_ID
        MY_INT_IOFF_125  = NACT_MAX  * NROOTS * MYNEW_ID
        MY_IOFF_347      = IOFF_PROC * ICOPY * MYNEW_ID
        MY_INT_IOFF_347  = NACT_MAX  * ICOPY * MYNEW_ID
C
        startcopy = MPI_WTIME()
C
        IF( CHUNK_FILE ) THEN
C
C       allocation for work-control arrays
C
          CALL MEMMAN(KCONTRL1,NBLK_ACT*NROOTS,'ADDL  ',1,'ICTRL1')
          CALL MEMMAN(KCONTRL2,NBLK_ACT*NROOTS,'ADDL  ',1,'ICTRL2')
          CALL MEMMAN(KCONTRL5,NBLK_ACT*NROOTS,'ADDL  ',1,'ICTRL5')
          CALL MEMMAN(KCONTRL3,NBLK_ACT*ICOPY,'ADDL  ',1,'ICTRL3')
          CALL MEMMAN(KCONTRL4,NBLK_ACT*ICOPY,'ADDL  ',1,'ICTRL4')
          CALL MEMMAN(KCONTRL7,NBLK_ACT*ICOPY,'ADDL  ',1,'ICTRL7')
          CALL ISETVC(WORK(KCONTRL1),0,NBLK_ACT*NROOTS)
          CALL ISETVC(WORK(KCONTRL2),0,NBLK_ACT*NROOTS)
          CALL ISETVC(WORK(KCONTRL5),0,NBLK_ACT*NROOTS)
          CALL ISETVC(WORK(KCONTRL3),0,NBLK_ACT*ICOPY)
          CALL ISETVC(WORK(KCONTRL4),0,NBLK_ACT*ICOPY)
          CALL ISETVC(WORK(KCONTRL7),0,NBLK_ACT*ICOPY)
C
          CALL COPY_STOMPI_DRV2(LU1,LU2,LU3,LU4,VEC1,NBLK_ACT,
     &                          NBATCH,
     &                          WORK(KLSLBT),WORK(KLSLEBT),
     &                          WORK(KLSI1BT),WORK(KLSIBT),
     &                          NROOTS,ICOPY,MY_IOFF_125,
     &                          MY_IOFF_347,MY_VEC_IOFF,
     &                          MY_INT_IOFF_125,
     &                          MY_INT_IOFF_347,
     &                          WORK(KCONTRL1),WORK(KCONTRL2),
     &                          WORK(KCONTRL3),WORK(KCONTRL4))
C
        ELSE
C
          CALL COPY_STOMPI_DRV1(LU1,LU2,LU3,LU4,VEC1,NBLK_ACT,
     &                          NBATCH,
     &                          WORK(KLSLBT),WORK(KLSLEBT),
     &                          WORK(KLSI1BT),WORK(KLSIBT),
     &                          NROOTS,ICOPY,MY_IOFF_125,
     &                          MY_IOFF_347,MY_VEC_IOFF)
C
        END IF
C
        finishcopy = MPI_WTIME()
        IPRINTSTAT = 1
        IF ( IPRINTSTAT .eq. 1 ) then
          WALLTIME = SECTID( finishcopy - startcopy )
          WRITE(LUWRTOUT,9450) WALLTIME
        END IF
C?        IF( LUCI_MYPROC .eq. LUCI_MASTER ) THEN
C?         write(6,*) 'these are my integer arrays!',myproc
C?         CALL IWRTMA(WORK(KCONTRL1),1,NBLK_ACT*NROOTS,1,NBLK_ACT*NROOTS)
C?         CALL IWRTMA(WORK(KCONTRL2),1,NBLK_ACT*NROOTS,1,NBLK_ACT*NROOTS)
C?         CALL IWRTMA(WORK(KCONTRL3),1,NBLK_ACT*ICOPY,1,NBLK_ACT*ICOPY)
C?         CALL IWRTMA(WORK(KCONTRL4),1,NBLK_ACT*ICOPY,1,NBLK_ACT*ICOPY)
C?        END IF
C
        IALLOCLEN = 0
        IALLOCLEN = MAX(MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
        CALL MEMMAN(KSCRRED2,IALLOCLEN,'ADDL  ',2,'SCRRD2')
        CALL SETVEC(WORK(KSCRRED2),0.0D0,IALLOCLEN)
C
C       setup is done, now we can call the proper working horse
C
        IF( CHUNK_FILE ) THEN
C
        CALL DIREC_OPTWRK_3(VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &                      NROOTS,MAXVEC,NVEC,APROJ,AVEC,WORK_SP,
     &                      LBLK,IBKINF,NBLOCKDN,IGROUPLIST,IPROCLIST,
     &                      BSCRFAC,IADD,SCRRED,
     &                      WORK(KSCRRED2),
     &                      LUWRTOUT,
     &                      NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                      WORK(KLSI1BT),WORK(KLSIBT),
     &                      MY_VEC_IOFF,MY_IOFF_125,MY_IOFF_347,
     &                      ICOPY,MY_INT_IOFF_125,MY_INT_IOFF_347,
     &                      NBLK_ACT,
     &                      WORK(KCONTRL1),WORK(KCONTRL2),
     &                      WORK(KCONTRL5),WORK(KCONTRL3),
     &                      WORK(KCONTRL4),WORK(KCONTRL7))
C
        ELSE
C
        CALL DIREC_OPTWRK_2(VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &                      NROOTS,MAXVEC,NVEC,APROJ,AVEC,WORK_SP,
     &                      LBLK,IBKINF,NBLOCKDN,IGROUPLIST,IPROCLIST,
     &                      BSCRFAC,IADD,SCRRED,
     &                      WORK(KSCRRED2),
     &                      LUWRTOUT,
     &                      NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                      WORK(KLSI1BT),WORK(KLSIBT),
     &                      MY_VEC_IOFF,MY_IOFF_125,MY_IOFF_347,
     &                      ICOPY,NBLK_ACT)
C
        END IF
C
C       ... all MPI-File I/O parts finished
C
      ELSE
C
        IF(MAX(NVEC-NROOT,NROOTS).LE.NVEC-NROOTS)THEN
               ICOPY = MAX(NVEC-NROOTS,NROOTS)
        ELSE
               ICOPY = MIN(NVEC-NROOTS,NROOTS)
        END IF
C
        CALL DIREC_OPTWRK(VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &                    NROOTS,MAXVEC,NVEC,APROJ,AVEC,WORK_SP,
     &                    LBLK,IBKINF,NBLOCKDN,IGROUPLIST,IPROCLIST,
     &                    BSCRFAC,IADD,SCRRED,LUWRTOUT,
     &                    NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                    WORK(KLSI1BT),WORK(KLSIBT))
C
C
      END IF
C
C     ... optimization finished ... back to the main DAVIDSON routine
C
C     Eliminate local memory
      IDUM = 0
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'D_OPTD')
C
 9450 FORMAT(' >>>  WALL TIME FOR COPY ROUTINES SQ --> MPI     :
     & ',A)
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DIREC_OPTWRK(VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &                        NROOT,MAXVEC,NVEC,APROJ,AVEC,WORK_SP,
     &                        LBLK,IBKINF,NBLOCKDN,IGROUPLIST,IPROCLIST,
     &                        BSCRFAC,IADD,SCRRED,LUWRTOUT,
     &                        NBATCH,LBATS,LEBATS,I1BATS,IBATS)
C                                                                       
C     working routine for the generation of a new best direction for 
C     the c-/s-vectors.        
C                                                                      
C     file structure is adapted to batch scheme used also in MV7_PAR
C     no MPI-FILE I/O is used inside this routine.                            
C
C     Last revision: S. Knecht           - March  2007
C                   
************************************************************************
#include "implicit.h"
#include "parluci.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      DIMENSION VEC1(*), VEC2(*)
C     batch structure of the CI vector
      DIMENSION LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
      DIMENSION IBKINF(NUM_BLOCKS2), NBLOCKDN(NUM_BLOCKS2)
      DIMENSION IPROCLIST(LUCI_NMPROC), IGROUPLIST(LUCI_NMPROC)
      INTEGER ISTRED, IMUSTRED
      INTEGER NNVEC
C
C     Dimensioning required of local vectors
C     APROJ     : MAXVEC*(MAXVEC+1)/2
C     AVEC      : MAXVEC ** 2
C     WORK_SP   : MAXVEC*(MAXVEC+1)/2
C     SCRRED    : MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
      DIMENSION APROJ(*),AVEC(*),WORK_SP(*), SCRRED(*)
      REAL*8 INPRDD, INPROD
      DOUBLE PRECISION tottime,endtime,starttime
      CHARACTER SECTID*12, WALLTID*12, WALLTIME

C
      ISTRED   = 0
      IMUSTRED = 0
      NNVEC    = 0
C
C     3.1 : Current wave function approximations, collect on LU7
C
      CALL REWINE( LU7,LBLK)
      DO IROOT = 1, NROOT
C       From LU1 to LU5
        CALL REWINO(LU1)
        CALL MVCSMD_PAR(LU1,AVEC((IROOT-1)*NVEC+1),
     &  LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of c reset, part 1'
C
C       and add LU3 to LU5
C
        CALL REWINO(LU3)
        ONE = 1.0D0
        CALL MVCSMD2_PAR(LU3,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &  LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of c reset, part 2'
*
        SCALE  = 1.0D0
        WORK_SP(IROOT) = SCALE
        CALL REWINO(LU5)
        CALL COPVCD(LU5,LU7,VEC1,0,LBLK)
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of c reset, part 3'
      END DO
C     Transfer C vectors from LU7 to LU1
      CALL REWINE( LU7,LBLK)
      CALL REWINE( LU1,LBLK)
      DO IVEC = 1,NROOT
        CALL COPVCD(LU7,LU1,VEC1,0,LBLK)
      END DO
      IF(IPRT.GE.1000)
     &write(6,*) ' end of c reset, part 4'
      IF(IPRT.GE.1000) THEN
        WRITE(6,*) ' Reset C-vectors on LU1 (myproc)',myproc
        CALL REWINO(LU1)
        DO IVEC = 1, NROOT
           CALL WRTVCD(VEC1,LU1,0,LBLK)
        END DO
      END IF
C     3.2 : corresponding sigma vectors
      CALL REWINE( LU7,LBLK)
      DO IROOT = 1, NROOT
C       From LU2 to LU5
        CALL REWINO(LU2)
        CALL MVCSMD_PAR(LU2,AVEC((IROOT-1)*NVEC+1),
     &                  LU5,LU6,VEC1,VEC2,NROOT,1,LBLK)
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of s reset, part 1'
C       and add LU4 to LU5
        CALL REWINO(LU4)
        CALL MVCSMD2_PAR(LU4,AVEC((IROOT-1)*NVEC+NROOT+1),ONE,
     &                   LU5,LU6,VEC1,VEC2,NVEC-NROOT,1,LBLK)
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of s reset, part 2'
        SCALE  = WORK_SP(IROOT)
C       scale LU5 => LU7
        CALL REWINO(LU5)
        CALL COPVCD(LU5,LU7,VEC1,0,LBLK)
        IF(IPRT.GE.1000)
     &  write(6,*) ' end of s reset, part 3'
      END DO
C     Transfer sigma  vectors from LU7 to LU2
      CALL REWINE( LU7,LBLK)
      CALL REWINE( LU2,LBLK)
      DO IVEC = 1,NROOT
        CALL COPVCD(LU7,LU2,VEC1,0,LBLK)
      END DO
      IF(IPRT.GE.1000)
     &  write(6,*) ' end of s reset, part 4'
C
C     next step
C
      NNVEC = NROOT
C
      IF(3*NROOT.LE.MAXVEC) THEN
C
C       Orthogonalize the
C       last set of correction vectors to the current
C       eigenvectors on LU1, and save on LU2
C       Overlap with root approximations
C       Start of last set of trial vectors
C
        ISTART = NVEC-NROOT-IADD+1
        CALL SKPVCD(LU3,ISTART-1,VEC1,1,LBLK)
C
        CALL REWINO(LU5)
        DO JVEC = 1, IADD
C         Orthogonalize to vectors on LU1
          CALL REWINO(LU7)
          CALL COPVCD(LU3,LU7,VEC1,0,LBLK)
C
          CALL REWINO(LU1)
          IMUSTRED = 0
          DO IROOT = 1, NROOT
            CALL REWINO(LU7)
            WORK_SP(IROOT+(JVEC-1)*2*NROOT) = 
     &     -INPRDD(VEC1,VEC2,LU7,LU1,0,LBLK)
            IF(IROOT.EQ.1) ISTRED = IROOT+(JVEC-1)*2*NROOT
            IMUSTRED = IMUSTRED + 1
          END DO
          starttime = MPI_WTIME()
          CALL REDVEC(WORK_SP(ISTRED),SCRRED,IMUSTRED,2,
     &                MPI_SUM,MPI_COMM_WORLD,-1)
          endtime = MPI_WTIME()
          WALLTID = SECTID(endtime-starttime)
          WRITE(LUWRTOUT,9440) WALLTID,IMUSTRED
          CALL COPVEC(SCRRED,WORK_SP(ISTRED),IMUSTRED)
C                 
          IF(IPRT.GE.1000)
     &    write(6,*) ' c, part1 finito'
C         and to trial vectors on LU5
          CALL REWINO(LU5)
          IMUSTRED = 0
          DO KVEC = 1, JVEC-1
            CALL REWINO(LU7)
            WORK_SP(NROOT+KVEC+(JVEC-1)*2*NROOT) = 
     &     -INPRDD(VEC1,VEC2,LU7,LU5,0,LBLK)
            IF(KVEC.EQ.1) ISTRED = NROOT+KVEC+(JVEC-1)*2*NROOT
            IMUSTRED = IMUSTRED + 1
          END DO
          IF(IMUSTRED.GT.0)THEN
            starttime = MPI_WTIME()
            CALL REDVEC(WORK_SP(ISTRED),SCRRED,IMUSTRED,2,
     &                  MPI_SUM,MPI_COMM_WORLD,-1)
            endtime = MPI_WTIME()
            WALLTID = SECTID(endtime-starttime)
            WRITE(LUWRTOUT,9450)WALLTID,IMUSTRED
            CALL COPVEC(SCRRED,WORK_SP(ISTRED),IMUSTRED)
          ENDIF
          WORK_SP(NROOT+JVEC+(JVEC-1)*2*NROOT) = 1.0D0
          IF(IPRT.GE.1000)
     &    write(6,*) ' c, part2 finito'
C
          ONE = 1.0D0
          CALL REWINO(LU1)
          CALL MVCSMD2_PAR(LU1,WORK_SP(1+(JVEC-1)*2*NROOT),ONE ,
     &                     LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
C
          ONE = 1.0D0
          CALL MVCSMD2_PAR(LU5,WORK_SP(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &                     LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
          IF(IPRT.GE.1000) THEN
            write(6,*) ' c, part4 finito'
            write(6,*) ' Vector after sec ort'
            CALL WRTVCD(VEC1,LU7,1,LBLK)
          END IF
C
          REDSCR = INPRDD(VEC1,VEC1,LU7,LU7,1,LBLK)
          CALL REDVEC(REDSCR,FACTOR,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
          SCALE = 1.0D0/SQRT(FACTOR)
          WRITE(6,*) ' SCALE in part 3.2 = ', SCALE
          CALL SCALVE(WORK_SP((JVEC-1)*2*NROOT+1),SCALE,
     &                NROOT+JVEC)
          CALL REWINO(LU7)
          CALL SCLVCD(LU7,LU5,SCALE,VEC1,0,LBLK)
          IF(IPRT.GE.1000)
     &      write(6,*) ' c, part6 finito'
        END DO
C       /\ End of loop over orthogonalized directions
C       Transfer modified directions to LU3
        CALL REWINO(LU3)
        CALL REWINO(LU5)
        DO JVEC =1, IADD
          CALL COPVCD(LU5,LU3,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000) THEN
          write(6,*) ' c, part7 finito'
          WRITE(6,*) ' Additional trial vectors on LU3 myproc',myproc
          CALL REWINO(LU3)
          DO JVEC = 1, IADD
            CALL WRTVCD(VEC1,LU3,0,LBLK)
          END DO
        END IF
C       Sigma vectors corresponding to orthogonalized directions
        CALL SKPVCD(LU4,ISTART-1,VEC1,1,LBLK)
        CALL REWINO(LU5)
        DO JVEC = 1, IADD
           CALL REWINO(LU7)
           CALL COPVCD(LU4,LU7,VEC1,0,LBLK)
C
           FACT = WORK_SP(NROOT+JVEC+(JVEC-1)*2*NROOT)
           CALL MVCSMD2_PAR(LU2,WORK_SP(1+(JVEC-1)*2*NROOT),FACT,
     &                      LU7,LU6,VEC1,VEC2,NROOT,1,LBLK)
C
           IF(IPRT.GE.1000)
     &      write(6,*) ' s, part 1 finito '
           ONE = 1.0D0
           CALL MVCSMD2_PAR(LU5,WORK_SP(NROOT+1+(JVEC-1)*2*NROOT),ONE,
     &                      LU7,LU6,VEC1,VEC2,JVEC-1,1,LBLK)
           IF(IPRT.GE.1000)
     &     write(6,*) ' s, part 2 finito '
           CALL REWINO(LU7)
           CALL COPVCD(LU7,LU5,VEC1,0,LBLK)
           IF(IPRT.GE.1000)
     &     write(6,*) ' s, part 3 finito '
        END DO
C       /\ End of loop over orthogonalized directions
C       Copy back to LU4
        CALL REWINO(LU4)
        CALL REWINO(LU5)
        DO JVEC = 1, IADD
           CALL COPVCD(LU5,LU4,VEC1,0,LBLK)
        END DO
        IF(IPRT.GE.1000)
     &  write(LUWRTOUT,*) ' s, part 4 finito '
C 
        NNVEC = NROOT + IADD
      END IF
C     ^ End if more than NROOT vectors to be reset
C
C
      IF(IPRT.GE.1000) THEN
        WRITE(LUWRTOUT,*)' Additional sigma vectors on LU4',myproc
        CALL REWINO(LU4)
        DO JVEC = 1, IADD
          CALL WRTVCD(VEC1,LU4,0,LBLK)
        END DO
      END IF
C
      NVEC = NNVEC
C
      IMUSTRED = 0
C     New subspace
C     Calculate subspace Hamiltonian from actual vectors
C     on disc
      IF(IPRT.GE.1000) write(6,*) ' Subspace hamiltonian'
        CALL REWINO(LU1)
        CALL REWINO(LU3)
C
      DO IVEC = 1, NVEC
C
        CALL REWINO(LU5)
        IF(IVEC.LE.NROOT) THEN
          CALL COPVCD(LU1,LU5,VEC1,0,LBLK)
        ELSE
          CALL COPVCD(LU3,LU5,VEC1,0,LBLK)
        END IF
C
        CALL REWINO(LU2)
        DO JVEC = 1, MIN(IVEC,NROOT)
          CALL REWINO(LU5)
          IJ = IVEC*(IVEC-1)/2+JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU2,0,LBLK)
          IF(IVEC.EQ.1.AND.JVEC.EQ.1) ISTRED = IJ
          IMUSTRED = IMUSTRED + 1
        END DO
        CALL REWINO(LU4)
        DO JVEC = NROOT+1,IVEC
          CALL REWINO(LU5)
          IJ = IVEC*(IVEC-1)/2+JVEC
          APROJ(IJ) = INPRDD(VEC1,VEC2,LU5,LU4,0,LBLK)
          IMUSTRED = IMUSTRED + 1
        END DO
      END DO
      CALL SETVEC(SCRRED,0.0D0,IMUSTRED)
      starttime = MPI_WTIME()
      CALL REDVEC(APROJ(ISTRED),SCRRED,IMUSTRED,2,MPI_SUM,
     &            MPI_COMM_WORLD,-1)
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      WRITE(LUWRTOUT,9460)WALLTID,IMUSTRED
      CALL COPVEC(SCRRED,APROJ(ISTRED),IMUSTRED)
C
 9440 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
 9450 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
 9460 FORMAT(' >>>  WALL TIME FOR REDUCING APROJ OF LENGTH     :
     & ',A,I4)
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DIREC_OPTWRK_2(VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &                          NROOT,MAXVEC,NVEC,APROJ,AVEC,WORK_SP,
     &                          LBLK,IBKINF,NBLOCKDN,IGROUPLIST,
     &                          IPROCLIST,BSCRFAC,IADD,
     &                          SCRRED,SCRRED2,LUWRTOUT,
     &                          NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                          MY_VEC_IOFF,MY_IOFF_125,MY_IOFF_347,
     &                          NOCOPIES,NBLK_ACT)

C                                                                       
C     working routine for the generation of a new best direction for 
C     the c-/s-vectors.        
C                                                                      
C     read in of blocks is adapted to batch scheme used also in MV7_PAR,
C     MPI-FILE I/O is used.
C
C     Last revision: S. Knecht           - March  2007
C                   
************************************************************************
#include "implicit.h"
#include "parluci.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      DIMENSION VEC1(*), VEC2(*)
      DIMENSION BSCRFAC(*)
C     batch structure of the CI vector
      DIMENSION LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
      DIMENSION IBKINF(NUM_BLOCKS2), NBLOCKDN(NUM_BLOCKS2)
      DIMENSION IPROCLIST(LUCI_NMPROC), IGROUPLIST(LUCI_NMPROC)
      INTEGER ISTRED, IMUSTRED, ISETSCRRED
      INTEGER NNVEC
C
C     Dimensioning required of local vectors
C     APROJ     : MAXVEC*(MAXVEC+1)/2
C     AVEC      : MAXVEC ** 2
C     WORK_SP   : MAXVEC*(MAXVEC+1)/2
C     SCRRED    : MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
C     SCRRED2   : MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
      DIMENSION APROJ(*),AVEC(*),WORK_SP(*), SCRRED(*), SCRRED2(*)
      REAL*8 INPRDD, INPROD
      DOUBLE PRECISION tottime,endtime,starttime
      CHARACTER SECTID*12, WALLTID*12, WALLTIME*12
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_125,MY_IOFF_347
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_START_ILU1
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_START_ILU2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_START_ILU3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_START_ILU4
C
      ISTRED     = 0
      IMUSTRED   = 0
      NNVEC      = 0
      IOFF_SCR   = 0
C
         call mpi_barrier(mpi_comm_world,ierr)
         DO II = 1, LUCI_NMPROC
         call mpi_barrier(mpi_comm_world,ierr)
         call mpi_barrier(mpi_comm_world,ierr)
           IF (LUCI_MYPROC.eq.II-1) THEN 
              WRITE(6,*) 'THIS IS MY BATCH LENGTH:',LUCI_MYPROC
           DO ISBATCH = 1, NBATCH
      WRITE(6,*) 'Length for batch ISBATCH',ISBATCH,'=',LEBATS(ISBATCH)
           END DO
           END IF
         call mpi_barrier(mpi_comm_world,ierr)
         call mpi_barrier(mpi_comm_world,ierr)
         END DO
         IF(LUCI_MYPROC.eq.LUCI_MASTER) THEN
         IOFF_SCR = MY_IOFF_125
         DO ISR = 1, NROOT
         DO ISBATCH = 1, NBATCH
          call setvec(VEC1,0.0D0,LEBATS(ISBATCH))
          call mpi_file_read_at(ILU1,IOFF_SCR,VEC1,LEBATS(ISBATCH),
     &                       MPI_REAL8,ISTAT,IERR)
          write(6,*) 'this is my batch isbatch for lu1',isbatch,IOFF_SCR
          write(6,*) 'this is  for lu1',isbatch,IOFF_SCR
          call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
          IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH) 
         END DO
         END DO
         ENDIF
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_barrier(mpi_comm_world,ierr)
      IOFF_SCR = MY_IOFF_347
      IF(myproc.eq.master) then
      DO ISR = 1, NOCOPIES
      DO ISBATCH = 1, NBATCH
       call mpi_file_read_at(ilu3,IOFF_SCR,VEC1,LEBATS(ISBATCH),
     &                       MPI_REAL8,ISTAT,IERR)
       write(6,*) 'this is my batch isbatch for lu3',isbatch,IOFF_SCR
       call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
       IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
      END DO
      END DO
      end if
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_barrier(mpi_comm_world,ierr)
C
C     3.1 : Current wave function approximations, collect on ILU5
C
      DO IROOT = 1, NROOT
        CALL MVCSMD_BATCH_PAR(VEC1,VEC2,AVEC((IROOT-1)*NVEC+1),
     &                        AVEC((IROOT-1)*NVEC+NROOT+1),NROOT,
     &                        NVEC-NROOT,
     &                        NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                        MY_VEC_IOFF,MY_IOFF_125,MY_IOFF_347,
     &                        MY_IOFF_125,NBLK_ACT,LU1,LU3,LU5,IROOT-1)
*
        WORK_SP(IROOT) = 1.0D0
      END DO
C     Transfer C vectors from ILU5 to ILU1
      DO IROOT = 1, NROOT
        CALL COPVCD_BATCH_PAR(VEC1,NBATCH,LBATS,LEBATS,
     &                        I1BATS,IBATS,MY_VEC_IOFF,
     &                        MY_IOFF_125,MY_IOFF_125,
     &                        NBLK_ACT,LU5,LU1,IROOT-1)
      END DO
C
C
C      IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
         IOFF_SCR = MY_IOFF_125
         WRITE(6,*) 'VECTOR(S) on ILU5',LUCI_MYPROC
         DO ISR = 1, NROOT
         DO ISBATCH = 1, NBATCH
          call mpi_file_read_at(ilu5,IOFF_SCR,VEC1,LEBATS(ISBATCH),
     &                       MPI_REAL8,ISTAT,IERR)
          write(6,*) 'this is my batch isbatch for lu5',isbatch,myproc
          call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
          IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
         END DO
         END DO

C      END IF
C
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_barrier(mpi_comm_world,ierr)
         IF( LUCI_MYPROC.EQ.LUCI_MASTER) THEN
         IOFF_SCR = MY_IOFF_125
         WRITE(6,*) 'reset VECTOR(S) on ILU1',LUCI_MYPROC
         DO ISR = 1, NROOT
         DO ISBATCH = 1, NBATCH
          CALL SETVEC(VEC1,0.0D0,LEBATS(ISBATCH))
          call mpi_file_read_at(ilu1,IOFF_SCR,VEC1,LEBATS(ISBATCH),
     &                       MPI_REAL8,ISTAT,IERR)
          write(6,*) 'this is my batch isbatch for lu1',isbatch,myproc
          write(6,*) 'this is the corresponding offset',IOFF_SCR
          call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
          IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
         END DO
         END DO
         END IF
         call mpi_barrier(mpi_comm_world,ierr)
         call mpi_barrier(mpi_comm_world,ierr)
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_barrier(mpi_comm_world,ierr)
         IF( LUCI_MYPROC.EQ.LUCI_MASTER) THEN
         IOFF_SCR = MY_IOFF_125
         WRITE(6,*) 'reset VECTOR(S) on ILU5',LUCI_MYPROC
         DO ISR = 1, NROOT
         DO ISBATCH = 1, NBATCH
          CALL SETVEC(VEC1,0.0D0,LEBATS(ISBATCH))
          call mpi_file_read_at(ilu1,IOFF_SCR,VEC1,LEBATS(ISBATCH),
     &                       MPI_REAL8,ISTAT,IERR)
          write(6,*) 'this is my batch isbatch for lu5',isbatch,myproc
          write(6,*) 'this is the corresponding offset',IOFF_SCR
          call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
          IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
         END DO
         END DO
         END IF
         call mpi_barrier(mpi_comm_world,ierr)
         call mpi_barrier(mpi_comm_world,ierr)
CCC      call abend2('stefan wants me to stop here')
C                  
C     3.2 : corresponding sigma vectors
C
      DO IROOT = 1, NROOT
        CALL MVCSMD_BATCH_PAR(VEC1,VEC2,AVEC((IROOT-1)*NVEC+1),
     &                        AVEC((IROOT-1)*NVEC+NROOT+1),NROOT,
     &                        NVEC-NROOT,
     &                        NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                        MY_VEC_IOFF,MY_IOFF_125,MY_IOFF_347,
     &                        MY_IOFF_125,NBLK_ACT,LU2,LU4,LU5,IROOT-1)
C
      END DO
C     Transfer sigma vectors from ILU5 to ILU2
      DO IROOT = 1, NROOT
        CALL COPVCD_BATCH_PAR(VEC1,NBATCH,LBATS,LEBATS,
     &                        I1BATS,IBATS,MY_VEC_IOFF,
     &                        MY_IOFF_125,MY_IOFF_125,
     &                        NBLK_ACT,LU5,LU2,IROOT-1)
      END DO
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_barrier(mpi_comm_world,ierr)
      IF( LUCI_MYPROC.EQ.LUCI_MASTER ) THEN
      IOFF_SCR = MY_IOFF_125
      WRITE(6,*) 'VECTOR(S) on ILU2',LUCI_MYPROC
      DO ISR = 1, NROOT
      DO ISBATCH = 1, NBATCH
       call mpi_file_read_at(ilu2,IOFF_SCR,VEC1,LEBATS(ISBATCH),
     &                       MPI_REAL8,ISTAT,IERR)
       write(6,*) 'this is my batch isbatch for ilu2',isbatch,myproc
       call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
       IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
      END DO
      END DO
      END IF
      call mpi_barrier(mpi_comm_world,ierr)
      call mpi_barrier(mpi_comm_world,ierr)
CCC   write(6,*) 'bye bye'
CCC   stop
C
C     first part ready, go on ...
C
      NNVEC = NROOT
C
      IF(3*NROOT.LE.MAXVEC) THEN
C
C       Orthogonalize the
C       last set of correction vectors to the current
C       eigenvectors on LU1, and save on LU2
C       Overlap with root approximations
C       Start of last set of trial vectors
C
        ISTART = NVEC-NROOT-IADD+1
CSK        WRITE(6,*) 'I have to skip ISTART-1 vectors on ILU3',
CSK     &              ISTART-1,LUCI_MYPROC
        !****************!
        DO JVEC = 1, IADD
        !****************!
C
C         calculate offset w.r.t ISTART
C
          MY_IOFF_START_ILU3     = 0
          MY_IOFF_START_ILU3 = MY_IOFF_347 +
     &                       ( MY_VEC_IOFF * (ISTART - 1) ) +
     &                       ( MY_VEC_IOFF * (JVEC   - 1) )
C
C         Orthogonalize to vectors on LU1 and to trial vectors on LU5
C
          CALL ORTHG_VEC_BATCH_PAR(VEC1,VEC2,WORK_SP,SCRRED,SCRRED2,
     &                             NROOT,JVEC-1,
     &                             NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                             MY_VEC_IOFF,
     &                             MY_IOFF_START_ILU3,MY_IOFF_125,
     &                             MY_IOFF_125,MY_IOFF_347,
     &                             NBLK_ACT,LU3,LU1,LU5,LU7,LUWRTOUT)
C
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      IOFF_SCR = MY_IOFF_347 +  (MY_VEC_IOFF * (JVEC -1))
CSK      WRITE(6,*) 'VECTOR on ILU7',LUCI_MYPROC
CSK      DO ISBATCH = 1, NBATCH
CSK       call mpi_file_read_at(ilu7,IOFF_SCR,VEC1,LEBATS(ISBATCH),
CSK     &                       MPI_REAL8,ISTAT,IERR)
CSK       write(6,*) 'this is my batch isbatch for lu7',isbatch,myproc
CSK       call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
CSK       IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
CSK      END DO
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
C
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      IF( JVEC-1 .gt. 0) THEN
CSK      IOFF_SCR = MY_IOFF_125 +  (MY_VEC_IOFF * (JVEC - 2))
CSK      WRITE(6,*) 'VECTOR on ILU5',LUCI_MYPROC
CSK      DO ISBATCH = 1, NBATCH
CSK       call mpi_file_read_at(ilu5,IOFF_SCR,VEC1,LEBATS(ISBATCH),
CSK     &                       MPI_REAL8,ISTAT,IERR)
CSK       write(6,*) 'this is my batch isbatch for lu5',isbatch,myproc
CSK       call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
CSK       IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
CSK      END DO
CSK      END IF
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
C
          WORK_SP(NROOT+JVEC+(JVEC-1)*2*NROOT) = 1.0D0
C
CSK          IF(JVEC-1.GT.0) THEN
CSK            WRITE(6,*) 'WORK_SP(NROOT+1+(JVEC-1)*2*NROOT) = ',
CSK     &      WORK_SP(NROOT+1+(JVEC-1)*2*NROOT)
CSK          WRITE(6,*) 'WORK_SP(1+(JVEC-1)*2*NROOT)) = '
CSK          CALL WRTMT_LU(WORK_SP(1+(JVEC-1)*2*NROOT),1,NROOT,1,NROOT)
CSK          END IF
CSK          WRITE(6,*) 'WORK_SP(NROOT+1+(JVEC-1)*2*NROOT) = '
C
C         current vector to work on has been written to ILU7
C
          CALL ORTHG_VEC_BATCH_PAR2(VEC1,VEC2,
     &                              WORK_SP(1+(JVEC-1)*2*NROOT),
     &                              WORK_SP(NROOT+1+(JVEC-1)*2*NROOT),
     &                              NROOT,JVEC-1,
     &                              NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                              MY_VEC_IOFF,MY_IOFF_125,MY_IOFF_347,
     &                              MY_IOFF_125,MY_IOFF_125,MY_IOFF_347,
     &                              NBLK_ACT,LU1,LU7,LU5,
     &                              LU5,LU3,LUWRTOUT,
     &                              REDSCR,WORK_SP((JVEC-1)*2*NROOT+1),
     &                              NROOT+JVEC,IADD)
C
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      IF( JVEC-1 .eq. 0) THEN
CSK      IOFF_SCR = MY_IOFF_125 +  (MY_VEC_IOFF * (JVEC - 1))
CSK      WRITE(6,*) 'VECTOR that has been written to ILU5',LUCI_MYPROC
CSK      DO ISBATCH = 1, NBATCH
CSK       call mpi_file_read_at(ilu5,IOFF_SCR,VEC1,LEBATS(ISBATCH),
CSK     &                       MPI_REAL8,ISTAT,IERR)
CSK       write(6,*) 'this is my written batch isbatch for lu5',isbatch,
CSK     &             myproc
CSK       call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
CSK       IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
CSK      END DO
CSK      END IF
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
        !*****!
        END DO
        !*****!
C       end of loop over orthogonalized directions -- IADD
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      IOFF_SCR = MY_IOFF_347
CSK      DO JVEC = 1, IADD
CSK        DO ISBATCH = 1, NBATCH
CSK          call mpi_file_read_at(ilu3,IOFF_SCR,VEC1,LEBATS(ISBATCH),
CSK     &                         MPI_REAL8,ISTAT,IERR)
CSK          write(6,*) 'this is my batch isbatch for lu3 and offset',
CSK     &                isbatch,IOFF_SCR,myproc
CSK          call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
CSK          IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
CSK        END DO
CSK      END DO
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
C
C       Sigma vectors corresponding to orthogonalized directions
C
CSK        WRITE(6,*) 'I have to skip ISTART-1 vectors on ILU4',
CSK     &              ISTART-1,LUCI_MYPROC
        !****************!
        DO JVEC = 1, IADD
        !****************!
C
C         calculate offset w.r.t ISTART
C
          MY_IOFF_START_ILU4     = 0
          MY_IOFF_START_ILU4 = MY_IOFF_347 +
     &                       ( MY_VEC_IOFF * (ISTART - 1) ) +
     &                       ( MY_VEC_IOFF * (JVEC   - 1) )
C
          FACT = WORK_SP(NROOT+JVEC+(JVEC-1)*2*NROOT)
C?         WRITE(LUWRTOUT,*) 'FACT for LU4 is =', FACT
C
C         orthogonalize sigma vectors
C
          CALL ORTHG_VEC_BATCH_PAR3(VEC1,VEC2,
     &                              WORK_SP(1+(JVEC-1)*2*NROOT),FACT,
     &                              WORK_SP(NROOT+1+(JVEC-1)*2*NROOT),
     &                              NROOT,JVEC-1,
     &                              NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                              MY_VEC_IOFF,MY_IOFF_START_ILU4,
     &                              MY_IOFF_125,MY_IOFF_347,MY_IOFF_347,
     &                              NBLK_ACT,LU4,LU2,LU4,LU4,
     &                              LUWRTOUT,IADD)
C
C
        !*****!
        END DO
        !*****!
C       End of loop over orthogonalized directions
C
        NNVEC = NROOT + IADD
C
      END IF
C    /\ End if more than NROOT vectors to be reset
C
C
CSK      call mpi_barrier(mpi_comm_world,ierr)
CSK      call mpi_barrier(mpi_comm_world,ierr)
C
         WRITE(LUWRTOUT,*)' Additional sigma vectors on LU4',myproc
         call mpi_barrier(mpi_comm_world,ierr)
         call mpi_barrier(mpi_comm_world,ierr)
         IOFF_SCR = MY_IOFF_347
         DO JVEC = 1, IADD
           DO ISBATCH = 1, NBATCH
            call mpi_file_read_at(ilu4,IOFF_SCR,VEC1,LEBATS(ISBATCH),
     &                         MPI_REAL8,ISTAT,IERR)
             write(6,*) 'this is my batch isbatch for lu4 and offset',
     &                isbatch,IOFF_SCR,myproc
             call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
             IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
           END DO
         END DO
         call mpi_barrier(mpi_comm_world,ierr)
         call mpi_barrier(mpi_comm_world,ierr)
C
      NVEC = NNVEC
C
C     New subspace
C     Calculate subspace Hamiltonian from actual vectors
C     on disc
C
C     set offsets to 0
      MY_IOFF_START_ILU1     = 0
      MY_IOFF_START_ILU3     = 0
      ISTRED   = 0
      IMUSTRED = 0
      ILOOP    = 0
      ISETSCRRED = 0
      ISETSCRRED = MAX( ISETSCRRED, MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
CSK      IF(LUCI_MYPROC.EQ.LUCI_MASTER) WRITE(6,*)'THIS IS ISETSCRRED',ISETSCRRED
      CALL DZERO(SCRRED,ISETSCRRED)
C
CSK      write(LUWRTOUT,*) ' Subspace hamiltonian'
C
      DO IVEC = 1, NVEC
C
C       calculate offset for the corresponding file
        IF( IVEC .le. NROOT )THEN
C         start with ILU1
          MY_IOFF_START_ILU1 = MY_IOFF_125 +
     &                       ( MY_VEC_IOFF * (IVEC - 1) )
        ELSE
C       ... ILU3
          ILOOP = ILOOP + 1
C
          MY_IOFF_START_ILU3 = MY_IOFF_347 + 
     &                       ( MY_VEC_IOFF * ( ILOOP - 1) )
        END IF
C
CSK        IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
CSK          WRITE(6,*) 'IVEC,IMUSTRED',IVEC,IMUSTRED
CSK        END IF
C
        CALL CALC_SUBSPACE_H(VEC1,VEC2,SCRRED,IMUSTRED,IVEC,NROOT,
     &                       ISTRED,
     &                       NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                       MY_VEC_IOFF,
     &                       MY_IOFF_START_ILU1,MY_IOFF_START_ILU3,
     &                       MY_IOFF_125,MY_IOFF_347,
     &                       NBLK_ACT,LU1,LU3,LU2,LU4)
C
CSK      IF(LUCI_MYPROC.EQ.LUCI_MASTER) WRITE(6,*)'THIS IS ISTRED',ISTRED
      END DO
C
C     communicate new subspace hamiltonian ...
C
      CALL DZERO(APROJ(ISTRED),IMUSTRED)
      starttime = MPI_WTIME()
      CALL REDVEC(SCRRED,APROJ(ISTRED),IMUSTRED,2,MPI_SUM,
     &            MPI_COMM_WORLD,-1)
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      WRITE(LUWRTOUT,9460)WALLTID,IMUSTRED
CCC      call mpi_barrier(mpi_comm_world,ierr)
CCC      call mpi_barrier(mpi_comm_world,ierr)
CCC   call mpi_barrier(mpi_comm_world,ierr)
CCC   call mpi_barrier(mpi_comm_world,ierr)
CCC   WRITE(6,*)'NVEC,NROOT,IMUSTRED,MYP..',NVEC,NROOT,IMUSTRED,LUCI_MYPROC
CCC   IF( myproc .eq. master ) then
CCC     WRITE(6,*)'THIS IS MY NEW SUBSPACE HAMILTONIAN'
CCC     CALL WRTMT_LU(APROJ(ISTRED),1,IMUSTRED,1,IMUSTRED)
CCC   end if
CCC   call mpi_barrier(mpi_comm_world,ierr)
CCC      call mpi_barrier(mpi_comm_world,ierr)
CCC      write(luwrtout,*) 'Stefan forced me to stop'
CCC      call Abend2('i want to test')
C        ... now copy back!
         startcopy = MPI_WTIME()
C        LU1 and LU2 first..
         CALL REWINO(LU1)
         CALL REWINO(LU2)
         DO IROOT = 1, NROOT
C
C          find correct offset!
           MY_IOFF_START_ILU1 = MY_IOFF_125 +
     &                       ( MY_VEC_IOFF * (IROOT - 1) )
           MY_IOFF_START_ILU2 = MY_IOFF_125 +
     &                       ( MY_VEC_IOFF * (IROOT - 1) )
C
           CALL COPVC_MPITOS_BDRV(LU1,LU1,VEC1,NBLK_ACT,
     &                            NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                            MY_IOFF_START_ILU1)
C
           CALL COPVC_MPITOS_BDRV(LU2,LU2,VEC1,NBLK_ACT,
     &                            NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                            MY_IOFF_START_ILU2)
C
         END DO
C
         CALL REWINO(LU3)
         CALL REWINO(LU4)
C        ICOPY vectors on LU3 and LU4
C
         DO IVEC = 1, NOCOPIES
C
           MY_IOFF_START_ILU3 = MY_IOFF_347 +
     &                       ( MY_VEC_IOFF * (IVEC - 1) )
           MY_IOFF_START_ILU4 = MY_IOFF_347 +
     &                       ( MY_VEC_IOFF * (IVEC - 1) )
C
           CALL COPVC_MPITOS_BDRV(LU3,LU3,VEC1,NBLK_ACT,
     &                            NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                            MY_IOFF_START_ILU3)
C
           CALL COPVC_MPITOS_BDRV(LU4,LU4,VEC1,NBLK_ACT,
     &                            NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                            MY_IOFF_START_ILU4)
C
         END DO
         finishcopy = MPI_WTIME()
         IF ( IPRINTSTAT .eq. 1 ) then
           WALLTIME = SECTID(finishcopy-startcopy)
           WRITE(LUWRTOUT,9470) WALLTIME
         END IF
C
 9440 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
 9450 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
 9460 FORMAT(' >>>  WALL TIME FOR REDUCING APROJ OF LENGTH     :
     & ',A,I4)
 9470 FORMAT(' >>>  WALL TIME FOR COPY ROUTINES MPI --> SQ     :
     & ',A)
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DIREC_OPTWRK_3(VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &                          NROOT,MAXVEC,NVEC,APROJ,AVEC,WORK_SP,
     &                          LBLK,IBKINF,NBLOCKDN,IGROUPLIST,
     &                          IPROCLIST,BSCRFAC,IADD,
     &                          SCRRED,SCRRED2,LUWRTOUT,
     &                          NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                          MY_VEC_IOFF,MY_IOFF_125,MY_IOFF_347,
     &                          NOCOPIES,
     &                          MY_INT_IOFF_125,MY_INT_IOFF_347,
     &                          NBLK_ACT,
     &                          IVCOFF_1,IVCOFF_2,IVCOFF_5,
     &                          IVCOFF_3,IVCOFF_4,IVCOFF_7)

C                                                                       
C     working routine for the generation of a new best direction for 
C     the c-/s-vectors. Adopted for use of INTEGER arrays instead of a 
C     MPI-INTEGER file.
C                                                                      
C     read in of blocks is adapted to batch scheme used also in MV7_PAR,
C     MPI-FILE I/O is used.
C
C     Last revision: S. Knecht           - April  2007
C                   
************************************************************************
#include "implicit.h"
#include "parluci.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      DIMENSION VEC1(*), VEC2(*)
      DIMENSION BSCRFAC(*)
C     batch structure of the CI vector
      DIMENSION LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
      DIMENSION IBKINF(NUM_BLOCKS2), NBLOCKDN(NUM_BLOCKS2)
      DIMENSION IPROCLIST(LUCI_NMPROC), IGROUPLIST(LUCI_NMPROC)
      DIMENSION IVCOFF_1(*),IVCOFF_2(*),IVCOFF_5(*)
      DIMENSION IVCOFF_3(*),IVCOFF_4(*),IVCOFF_7(*)
      INTEGER ISTRED, IMUSTRED, ISETSCRRED
      INTEGER NNVEC
C
C     Dimensioning required of local vectors
C     APROJ     : MAXVEC*(MAXVEC+1)/2
C     AVEC      : MAXVEC ** 2
C     WORK_SP   : MAXVEC*(MAXVEC+1)/2
C     SCRRED    : MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
C     SCRRED2   : MAX( MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
C
      DIMENSION APROJ(*),AVEC(*),WORK_SP(*), SCRRED(*), SCRRED2(*)
      DOUBLE PRECISION tottime,endtime,starttime
      CHARACTER SECTID*12, WALLTID*12, WALLTIME*12
      INTEGER(KIND=MPI_OFFSET_KIND) MY_VEC_IOFF
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_125,MY_IOFF_347
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_START_ILU1
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_START_ILU2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_START_ILU3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_START_ILU4
C
      ISTRED     = 0
      IMUSTRED   = 0
      NNVEC      = 0
      IOFF_SCR   = 0
C
C     3.1 : Current wave function approximations, collect on ILU5
C
      DO IROOT = 1, NROOT
C        
        CALL MVCSMD_BATCH_PAR_USE_IA(VEC1,VEC2,AVEC((IROOT-1)*NVEC+1),
     &                               AVEC((IROOT-1)*NVEC+NROOT+1),NROOT,
     &                               NVEC-NROOT,
     &                               NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                               MY_VEC_IOFF,MY_IOFF_125,
     &                               MY_IOFF_347,MY_IOFF_125,
     &                               NBLK_ACT,LU1,LU3,LU5,
     &                               IROOT-1,IVCOFF_1,IVCOFF_3,IVCOFF_5)
C
        WORK_SP(IROOT) = 1.0D0
      END DO
C     Transfer C vectors from ILU5 to ILU1
      DO IROOT = 1, NROOT
        CALL COPVCD_BATCH_PAR2(VEC1,NBATCH,LBATS,LEBATS,
     &                         I1BATS,IBATS,MY_VEC_IOFF,
     &                         MY_IOFF_125,MY_IOFF_125,
     &                         NBLK_ACT,LU5,LU1,
     &                         IROOT-1,IVCOFF_5,IVCOFF_1)
      END DO
C                  
C     3.2 : corresponding sigma vectors
C
      DO IROOT = 1, NROOT
C
        CALL MVCSMD_BATCH_PAR_USE_IA(VEC1,VEC2,AVEC((IROOT-1)*NVEC+1),
     &                               AVEC((IROOT-1)*NVEC+NROOT+1),NROOT,
     &                               NVEC-NROOT,
     &                               NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                               MY_VEC_IOFF,MY_IOFF_125,
     &                               MY_IOFF_347,MY_IOFF_125,
     &                               NBLK_ACT,LU2,LU4,
     &                               LU5,IROOT-1,
     &                               IVCOFF_2,IVCOFF_4,IVCOFF_5)
C
      END DO
C     Transfer sigma vectors from ILU5 to ILU2
      DO IROOT = 1, NROOT
        CALL COPVCD_BATCH_PAR2(VEC1,NBATCH,LBATS,LEBATS,
     &                         I1BATS,IBATS,MY_VEC_IOFF,
     &                         MY_IOFF_125,MY_IOFF_125,
     &                         NBLK_ACT,LU5,LU2,
     &                         IROOT-1,IVCOFF_5,IVCOFF_2)
      END DO
C
C     first part ready, go on ...
C
C
CCC   call mpi_barrier(mpi_comm_world,ierr)
CCC   call mpi_barrier(mpi_comm_world,ierr)
CCC   IF( LUCI_MYPROC.EQ.LUCI_MASTER +2 ) THEN
CCC   IOFF_SCR = MY_IOFF_125
CCC   WRITE(6,*) 'VECTOR(S) on ILU2',LUCI_MYPROC
CCC   DO ISR = 1, NROOT
CCC   DO ISBATCH = 1, NBATCH
CCC    call mpi_file_read_at(ilu2,IOFF_SCR,VEC1,LEBATS(ISBATCH),
CCC  &                       MPI_REAL8,ISTAT,IERR)
CCC    write(6,*) 'this is my batch isbatch for ilu2',isbatch,myproc
CCC    call WRTMT_LU(VEC1,1,LEBATS(ISBATCH),1,LEBATS(ISBATCH))
CCC    IOFF_SCR = IOFF_SCR + LEBATS(ISBATCH)
CCC   END DO
CCC   END DO
CCC   END IF
CCC   call mpi_barrier(mpi_comm_world,ierr)
CCC   call mpi_barrier(mpi_comm_world,ierr)
CCC   write(6,*) 'bye bye'
CCC   stop
C
      NNVEC = NROOT
C
      IF(3*NROOT.LE.MAXVEC) THEN
C
C       Orthogonalize the
C       last set of correction vectors to the current
C       eigenvectors on LU1, and save on LU2
C       Overlap with root approximations
C       Start of last set of trial vectors
C
        ISTART = NVEC-NROOT-IADD+1
CSK        WRITE(6,*) 'I have to skip ISTART-1 vectors on ILU3',
CSK     &              ISTART-1,LUCI_MYPROC
        !****************!
        DO JVEC = 1, IADD
        !****************!
C
C         calculate offset w.r.t ISTART
C
          MY_IOFF_START_ILU3     = 0
          MY_IOFF_START_INT_ILU3 = 0
          MY_IOFF_START_ILU3 = MY_IOFF_347 +
     &                       ( MY_VEC_IOFF * (ISTART - 1) ) +
     &                       ( MY_VEC_IOFF * (JVEC   - 1) )
C
          MY_IOFF_START_INT_ILU3 = ( NBLK_ACT * (ISTART - 1) )
     &                           + ( NBLK_ACT * (JVEC   - 1) )
     &                           + 1
C
C         Orthogonalize to vectors on LU1 and to trial vectors on LU5
C
          CALL ORTHG_VEC_BATCH_PAR_USE_IA(VEC1,VEC2,WORK_SP,SCRRED,
     &                                    SCRRED2,NROOT,JVEC-1,
     &                                    NBATCH,LBATS,LEBATS,
     &                                    I1BATS,IBATS,
     &                                    MY_VEC_IOFF,
     &                                    MY_IOFF_START_ILU3,
     &                                    MY_IOFF_125,MY_IOFF_125,
     &                                    MY_IOFF_347,NBLK_ACT,
     &                                    LU3,LU1,LU5,LU7,LUWRTOUT,
     &                                    MY_IOFF_START_INT_ILU3,
     &                                    IVCOFF_3,IVCOFF_1,IVCOFF_5,
     &                                    IVCOFF_7)
C
C
          WORK_SP(NROOT+JVEC+(JVEC-1)*2*NROOT) = 1.0D0
C
C         current vector to work on has been written to ILU7
C
          CALL ORTHG_VEC_BATCH_PAR2_USE_IA(VEC1,VEC2,
     &                                     WORK_SP(1+(JVEC-1)*2*NROOT),
     &                                WORK_SP(NROOT+1+(JVEC-1)*2*NROOT),
     &                                     NROOT,JVEC-1,NBATCH,LBATS,
     &                                     LEBATS,I1BATS,IBATS,
     &                                     MY_VEC_IOFF,MY_IOFF_125,
     &                                     MY_IOFF_347,MY_IOFF_125,
     &                                     MY_IOFF_125,MY_IOFF_347,
     &                                     NBLK_ACT,
     &                                     LU1,LU7,LU5,LU5,LU3,LUWRTOUT,
     &                                     REDSCR,
     &                                     WORK_SP((JVEC-1)*2*NROOT+1),
     &                                     NROOT+JVEC,IADD,
     &                                     IVCOFF_1,IVCOFF_7,IVCOFF_5,
     &                                     IVCOFF_5,IVCOFF_3)
C
        !*****!
        END DO
        !*****!
C       end of loop over orthogonalized directions -- IADD
C
C       Sigma vectors corresponding to orthogonalized directions
C
CSK        WRITE(6,*) 'I have to skip ISTART-1 vectors on ILU4',
CSK     &              ISTART-1,LUCI_MYPROC
        !****************!
        DO JVEC = 1, IADD
        !****************!
C
C         calculate offset w.r.t ISTART
C
          MY_IOFF_START_ILU4     = 0
          MY_IOFF_START_INT_ILU4 = 0
          MY_IOFF_START_ILU4 = MY_IOFF_347 +
     &                       ( MY_VEC_IOFF * (ISTART - 1) ) +
     &                       ( MY_VEC_IOFF * (JVEC   - 1) )
C
          MY_IOFF_START_INT_ILU4 = 1 +
     &                           ( NBLK_ACT * (ISTART - 1) ) +
     &                           ( NBLK_ACT * (JVEC   - 1) )
C
          FACT = WORK_SP(NROOT+JVEC+(JVEC-1)*2*NROOT)
C?         WRITE(LUWRTOUT,*) 'FACT for LU4 is =', FACT
C
C         orthogonalize sigma vectors
C
          CALL ORTHG_VEC_BATCH_PAR3_USE_IA(VEC1,VEC2,
     &                                     WORK_SP(1+(JVEC-1)*2*NROOT),
     &                                     FACT,
     &                                WORK_SP(NROOT+1+(JVEC-1)*2*NROOT),
     &                                     NROOT,JVEC-1,
     &                                     NBATCH,LBATS,LEBATS,
     &                                     I1BATS,IBATS,
     &                                     MY_VEC_IOFF,
     &                                     MY_IOFF_START_ILU4,
     &                                     MY_IOFF_125,MY_IOFF_347,
     &                                     MY_IOFF_347,
     &                                     NBLK_ACT,
     &                                     LU4,LU2,LU4,LU4,LUWRTOUT,
     &                                     IADD,MY_IOFF_START_INT_ILU4,
     &                                     IVCOFF_4,IVCOFF_2,
     &                                     IVCOFF_4,IVCOFF_4)
C
        !*****!
        END DO
        !*****!
C       End of loop over orthogonalized directions
C
        NNVEC = NROOT + IADD
C
      END IF
C    /\ End if more than NROOT vectors to be reset
C
C
      NVEC = NNVEC
C
C     New subspace
C     Calculate subspace Hamiltonian from actual vectors
C     on disc
C
C     set offsets to 0
      MY_IOFF_START_ILU1     = 0
      MY_IOFF_START_INT_ILU1 = 0
      MY_IOFF_START_ILU3     = 0
      MY_IOFF_START_INT_ILU3 = 0
      ISTRED   = 0
      IMUSTRED = 0
      ILOOP    = 0
      ISETSCRRED = 0
      ISETSCRRED = MAX( ISETSCRRED, MAXVEC*(MAXVEC+1)/2 , MAXVEC ** 2)
CSK      IF(LUCI_MYPROC.EQ.LUCI_MASTER) WRITE(6,*)'THIS IS ISETSCRRED',ISETSCRRED
      CALL DZERO(SCRRED,ISETSCRRED)
C
CSK      write(LUWRTOUT,*) ' Subspace hamiltonian'
C
      DO IVEC = 1, NVEC
C
C       calculate offset for the corresponding file
        IF( IVEC .le. NROOT )THEN
C         start with ILU1
          MY_IOFF_START_ILU1 = MY_IOFF_125 +
     &                       ( MY_VEC_IOFF * (IVEC - 1) )
          MY_IOFF_START_INT_ILU1 = 1 +
     &                          ( NBLK_ACT * (IVEC - 1) )
        ELSE
C       ... ILU3
          ILOOP = ILOOP + 1
C
          MY_IOFF_START_ILU3 = MY_IOFF_347 + 
     &                       ( MY_VEC_IOFF * ( ILOOP - 1) )
          MY_IOFF_START_INT_ILU3 = 1 + 
     &                          ( NBLK_ACT * ( ILOOP - 1) )
        END IF
C
CSK        IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
CSK          WRITE(6,*) 'IVEC,IMUSTRED',IVEC,IMUSTRED
CSK        END IF
C
        CALL CALC_SUBSPACE_H_USE_IA(VEC1,VEC2,SCRRED,IMUSTRED,
     &                              IVEC,NROOT,ISTRED,NBATCH,LBATS,
     &                              LEBATS,I1BATS,IBATS,
     &                              MY_VEC_IOFF,
     &                              MY_IOFF_START_ILU1,
     &                              MY_IOFF_START_ILU3,
     &                              MY_IOFF_125,MY_IOFF_347,
     &                              NBLK_ACT,
     &                              LU1,LU3,LU2,LU4,
     &                              MY_IOFF_START_INT_ILU1,
     &                              MY_IOFF_START_INT_ILU3,
     &                              IVCOFF_1,IVCOFF_3,
     &                              IVCOFF_2,IVCOFF_4,IVCOFF_5)
C
CSK      IF(LUCI_MYPROC.EQ.LUCI_MASTER) WRITE(6,*)'THIS IS ISTRED',ISTRED
      END DO
C
C     communicate new subspace hamiltonian ...
C
      CALL DZERO(APROJ(ISTRED),IMUSTRED)
      starttime = MPI_WTIME()
      CALL REDVEC(SCRRED,APROJ(ISTRED),IMUSTRED,2,MPI_SUM,
     &            MPI_COMM_WORLD,-1)
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
      WRITE(LUWRTOUT,9460)WALLTID,IMUSTRED
CCC   IF( myproc .eq. master ) then
CCC     WRITE(6,*)'THIS IS MY NEW SUBSPACE HAMILTONIAN'
CCC     CALL WRTMT_LU(APROJ(ISTRED),1,IMUSTRED,1,IMUSTRED)
CCC   end if
CCC   call mpi_barrier(mpi_comm_world,ierr)
CCC      call mpi_barrier(mpi_comm_world,ierr)
CCC      write(luwrtout,*) 'Stefan forced me to stop'
CCC      call Abend2('i want to test')
C     ... now copy back!
      startcopy = MPI_WTIME()
C     LU1 and LU2 first..
      CALL REWINO(LU1)
      CALL REWINO(LU2)
      DO IROOT = 1, NROOT
C
C       find correct offset!
        MY_IOFF_START_ILU1 = MY_IOFF_125 +
     &                    ( MY_VEC_IOFF * (IROOT - 1) )
        MY_IOFF_START_INT_ILU1 = 1 +
     &                       ( NBLK_ACT * (IROOT - 1) )
        MY_IOFF_START_ILU2 = MY_IOFF_125 +
     &                    ( MY_VEC_IOFF * (IROOT - 1) )
        MY_IOFF_START_INT_ILU2 = 1 +
     &                       ( NBLK_ACT * (IROOT - 1) )
C
        CALL COPVC_MPITOS_BDRV2(LU1,LU1,VEC1,NBLK_ACT,
     &                          NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                          MY_IOFF_START_ILU1,
     &                          MY_IOFF_START_INT_ILU1,IVCOFF_1)
C
        CALL COPVC_MPITOS_BDRV2(LU2,LU2,VEC1,NBLK_ACT,
     &                          NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                          MY_IOFF_START_ILU2,
     &                          MY_IOFF_START_INT_ILU2,IVCOFF_2)
C
      END DO
C
      CALL REWINO(LU3)
      CALL REWINO(LU4)
C     ICOPY vectors on LU3 and LU4
C
      DO IVEC = 1, NOCOPIES
C
        MY_IOFF_START_ILU3 = MY_IOFF_347 +
     &                    ( MY_VEC_IOFF * (IVEC - 1) )
        MY_IOFF_START_INT_ILU3 = 1 +
     &                       ( NBLK_ACT * (IVEC - 1) )
        MY_IOFF_START_ILU4 = MY_IOFF_347 +
     &                    ( MY_VEC_IOFF * (IVEC - 1) )
        MY_IOFF_START_INT_ILU4 = 1 +
     &                       ( NBLK_ACT * (IVEC - 1) )
C
        CALL COPVC_MPITOS_BDRV2(LU3,LU3,VEC1,NBLK_ACT,
     &                          NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                          MY_IOFF_START_ILU3,
     &                          MY_IOFF_START_INT_ILU3,IVCOFF_3)
C
        CALL COPVC_MPITOS_BDRV2(LU4,LU4,VEC1,NBLK_ACT,
     &                          NBATCH,LBATS,LEBATS,I1BATS,IBATS,
     &                          MY_IOFF_START_ILU4,
     &                          MY_IOFF_START_INT_ILU4,IVCOFF_4)
C
      END DO
      finishcopy = MPI_WTIME()
      IF ( IPRINTSTAT .eq. 1 ) THEN
        WALLTIME = SECTID(finishcopy-startcopy)
        WRITE(LUWRTOUT,9470) WALLTIME
      END IF
C
 9440 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
 9450 FORMAT(' >>>  WALL TIME FOR REDUCING WORK OF LENGTH      :
     & ',A,I4)
 9460 FORMAT(' >>>  WALL TIME FOR REDUCING APROJ OF LENGTH     :
     & ',A,I4)
 9470 FORMAT(' >>>  WALL TIME FOR COPY ROUTINES MPI --> SQ     :
     & ',A)
C
      END
C
#else
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
