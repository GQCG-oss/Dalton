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

***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ADD_TTS_BLK_TO_VECTOR
     &           (NBLOCKI,IBLOCKI,IOFFI,VECI,
     &            NBLOCKO,IBLOCKO,IOFFO,VECO,FACTOR)
* A vector VECI containing NBLOCKI blocks defined by IBLOCKI
* are given. Add those blocks to a vector VECO,
* defined by NBLOCKO,IBLOCKO,IOFFO,VECO as
*
*  VECO = VECO + FACTOR*VECI
*
* Jeppe Olsen,  November 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBLOCKI(8,*)
      INTEGER IBLOCKO(8,*)
      DIMENSION VECI(*)
*. output
      DIMENSION VECO(*)
*
      DO JBLOCKI = IOFFI,IOFFI+NBLOCKI-1
        DO JBLOCKO = IOFFO, IOFFO + NBLOCKO-1
          IF(IBLOCKI(1,JBLOCKI).EQ. IBLOCKO(1,JBLOCKO).AND.
     &    IBLOCKI(2,JBLOCKI).EQ. IBLOCKO(2,JBLOCKO)   .AND.
     &    IBLOCKI(3,JBLOCKI).EQ. IBLOCKO(3,JBLOCKO)   .AND.
     &    IBLOCKI(4,JBLOCKI).EQ. IBLOCKO(4,JBLOCKO)          ) THEN
            IFROM = IBLOCKI(6,JBLOCKI)-IBLOCKI(6,IOFFI)+1
            ITO   = IBLOCKO(6,JBLOCKO)-IBLOCKO(6,IOFFO)+1
            NELMNT = IBLOCKI(8,JBLOCKI)
            ONE = 1.0D0
            CALL VECSUM(VECO(ITO),VECO(ITO),VECI(IFROM),
     &                  ONE,FACTOR,NELMNT)
          END IF
        END DO
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ADDBLKV(VEC,VECA,FACTOR,NBLKA,IBLKA,IBLOCK,IOFF)
*
* A blocked vector VEC is given
* add factor * VECADD to this
* mapping of VECADD blocks to VEC is given by IBLKA
*
* Jeppe Olsen, August 95
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(*),VECA(*)
*. Base and length of blocks in complete vector
      INTEGER IBLOCK(8,*)
*. Blocks to be added
      INTEGER IBLKA(NBLKA)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input to ADDBLKV '
        WRITE(6,*) ' Factor ', FACTOR
        WRITE(6,*) ' blocks to be updated '
        CALL IWRTMA(IBLKA,1,NBLKA,1,NBLKA)
      END IF
*
      IOFFI = 1
      MAXE = 0
      DO JBLK = 1, NBLKA
C       IOFFO  = IBLOCK(6,IBLKA(JBLK)-1+IOFF)
C       NELMNT = IBLOCK(8,IBLKA(JBLK)-1+IOFF)
        IOFFO  = IBLOCK(6,IBLKA(JBLK))
        NELMNT = IBLOCK(8,IBLKA(JBLK))
        MAXE = MAX(MAXE,IOFFO+NELMNT-1)
        ONE = 1.0D0
        CALL VECSUM(VEC(IOFFO),VEC(IOFFO),VECA(IOFFI),
     &              ONE,FACTOR,NELMNT)
        IOFFI = IOFFI + NELMNT
       END DO
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' output from ADDBLKV '
       CALL WRTMT_LU(VEC,1,MAXE,1,MAXE)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE BVEC(B,IBSM,LUC,LUB,VEC1,VEC2)
*
* Construct B defined as
*
*       B = B !0> (no projection)
*
* Where B is a one-electron operator with symmetry IBSM
* and !0> is assumed stored on LUC
*
* Jeppe Olsen , Jan. 1  1990
*               Sept97 : Modified to LUCIA
*
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION B(*),VEC1(*),VEC2(*)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "multd2h.inc"
#include "oper.inc"
#include "orbinp.inc"
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*. Ensure that symmetries are consistent
      ISSMP = MULTD2H(IBSM,ICSM)
      IF(ISSMP .NE. ISSM ) THEN
        WRITE(6,*) ' Inconsistent symmetries in BVEC'
        WRITE(6,*) ' ICSM, ISSM, IBSM ', ICSM,ISSM,IBSM
        Call Abend2(       ' Inconsistent symmetries in BVEC' )
      END IF
*. Get integrals in place
      CALL SWAPVE(B,WORK(KINT1),NTOOB**2)
*. Just one-electron operator
      I12 = 1
*. No additional approximations
      IAPR = 0
      IPERTOP = 0
*. And : DO IT
      CALL MV7(VEC1,VEC2,LUC,LUB)
*. Restore order
      CALL SWAPVE(B,WORK(KINT1),NTOOB**2)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE CLASS_PROD(VEC1,VEC2,NOCTPA,NOCTPB,
     &                      IBLOCK_OFF,NBLOCK,IBLOCK,NOCCLS,IOCCLS,
     &                      CLSVEC)
*
* Two vectors in blocked form are given.
* Find contributions to each occupation class
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      REAL*8 INPROD
*. Input
      DIMENSION VEC1(*),VEC2(*)
      INTEGER IBLOCK(8,*)
*. Input/output
      DIMENSION CLSVEC(*)

*
      IOFF = 1
      DO JBLOCK = IBLOCK_OFF,IBLOCK_OFF-1+NBLOCK
*
        JATP = IBLOCK(1,JBLOCK)
        JBTP = IBLOCK(2,JBLOCK)
        JASM = IBLOCK(3,JBLOCK)
        NELMNT = IBLOCK(8,JBLOCK)
C       NELMNT = NOOS(JATP,JBTP,JASM)
*. Corresponding occupation class
        CALL SPGSPG_TO_CLASS(JATP,JBTP,JOCCLS,NOCCLS,IOCCLS)
C?      WRITE(6,*)
C?   &  ' CLASS_PROD : JBLOCK, IOFF NELMNT ' ,JBLOCK,  IOFF, NELMNT
        IF(JOCCLS.EQ.0) THEN
          WRITE(6,*) ' JOCCLS = 0 returned from  SPGSPG_TO_CLASS'
        ELSE
          XTERM = INPROD(VEC1(IOFF),VEC2(IOFF),NELMNT)
C?        WRITE(6,*) ' CLASS_PROD : XTERM = ', XTERM
          CLSVEC(JOCCLS) = CLSVEC(JOCCLS) + XTERM
        END IF
        IOFF = IOFF + NELMNT
      END DO
*
      NTEST = 0
      IF(NTEST.GT.0) THEN
         WRITE(6,*) ' Updated CLSVEC '
         CALL WRTMT_LU(CLSVEC,1,NOCCLS,1,NOCCLS)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE EXPCIV(ISM,ISPCIN,LUIN,
     &                  ISPCUT,LUUT,LBLK,
     &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
*
* Expand CI vector in CI space ISPCIN to CI vector in ISPCUT
* Input vector is supposed to be on LUIN
* Output vector will be placed on unit LUUT
*. If ICOPY .ne. 0 the output vectors will be copied to LUIN
*
* Storage form is defined by ICISTR
*
* Jeppe Olsen, February 1994
* GAS version August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLIABI, KLIABU, KLBLIN, KLBLUT, KLTTSII, KLTTSNI, 
     &          KLTTSIU, KLTTSNU, KLBLI, KLBLU
!               for addressing of WORK
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "strbas.inc"
#include "stinf.inc"
#include "csm.inc"
#include "cgas.inc"
#include "gasstr.inc"

*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXPCIV')
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
*
*. Allowed combinations of strings types for input and output
*. spaces
*
      CALL MEMMAN(KLIABI,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABI')
      CALL MEMMAN(KLIABU,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABU')
      CALL IAIBCM_GAS(LCMBSPC(ISPCIN),ICMBSPC(1,ISPCIN),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KLIABI),NTEST)
*
      CALL IAIBCM_GAS(LCMBSPC(ISPCUT),ICMBSPC(1,ISPCUT),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KLIABU),NTEST)
*
* type of each symmetry block ( full, lower diagonal, absent )
*
        CALL MEMMAN(KLBLIN,NSMST,'ADDL  ',1,'KLBLIN')
        CALL MEMMAN(KLBLUT,NSMST,'ADDL  ',1,'KLBLUT')
        CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLIN),IDUMMY)
        CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLUT),IDUMMY)

*
*. Number of dets etc per TTS block
*

        CALL MEMMAN(KLTTSII,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSI')
        CALL MEMMAN(KLTTSNI,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSN')
        CALL ZOOS(ISMOST(1,ISM),WORK(KLBLIN),NSMST,WORK(KLIABI),
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &            NOCTPA,NOCTPB,IDC,WORK(KLTTSII),WORK(KLTTSNI),
     &            NCOMBI,0)
*
        CALL MEMMAN(KLTTSIU,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSI')
        CALL MEMMAN(KLTTSNU,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSN')
        CALL ZOOS(ISMOST(1,ISM),WORK(KLBLUT),NSMST,WORK(KLIABU),
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &            NOCTPA,NOCTPB,IDC,WORK(KLTTSIU),WORK(KLTTSNU),
     &            NCOMBU,0)
*
*. Allocate memory for blocks of input and output space in
*  accordance with ICISTR
*
      IF(ICISTR.EQ.1) THEN
        LENGTHI = NCOMBI
        LENGTHU = NCOMBU
      ELSE IF (ICISTR .EQ. 2 ) THEN
        LENGTHI = MXSB
        LENGTHU = MXSB
      ELSE IF (ICISTR.EQ. 3 ) THEN
        LENGTHI = MXSOOB
        LENGTHU = MXSOOB
      END IF
*
      CALL MEMMAN(KLBLI,LENGTHI,'ADDL  ',2,'KLBLI ')
      CALL MEMMAN(KLBLU,LENGTHU,'ADDL  ',2,'KLBLU ')

*
*. and now : Let another subroutine complete the taks
*
      CALL REWINE(LUIN,-1)
      CALL REWINE(LUUT,-1)
*
*. Print for testing initial vectors out
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Initial vectors in EXPCIV '
*
        DO IROOT = 1, NROOT
        WRITE(6,*) ' Root number ', IROOT
*
          IF(ICISTR.EQ.1) THEN
            CALL WRTMT_LU(WORK(KLBLI),1,NCOMBUT,1,NCOMBUT)
          ELSE
            CALL WRTVCD(WORK(KLBLI),LUIN,0,-1)
          END IF
        END DO
        CALL REWINE(LUIN,-1)
      END IF
*
      DO IROOT = 1, NROOT
*. Input vector should be first vector on file so
        IF(IROOT.EQ.1) THEN
          LLUIN = LUIN
        ELSE
          IF(ICISTR.EQ.1) THEN
            CALL REWINE(LUSCR,-1)
            CALL FRMDSC_LUCI(WORK(KLBLI),NCOMBI,-1,LUIN,IMZERO,IAMPACK)
            CALL  TODSC_LUCI(WORK(KLBLI),NCOMBI,-1,LUSCR)
            CALL REWINE(LUSCR,-1)
            LLUIN = LUSCR
          ELSE
*. With the elegance of an elephant
            CALL REWINE(LUSCR,-1)
            CALL REWINE(LUIN,-1)
            DO JROOT = 1, IROOT
              CALL REWINE(LUSCR,-1)
              CALL COPVCD(LUIN,LUSCR,WORK(KLBLI),0,-1)
            END DO
            CALL REWINE(LUSCR,-1)
            LLUIN = LUSCR
          END IF
        END IF
*. Expcivs may need the IAMPACK parameter ( in case it must write
*  a zero block before any blocks have been read in.
*  Use IDIAG to decide
       IF(IDIAG.EQ.1) THEN
         IAMPACK = 0
       ELSE
         IAMPACK = 1
       END IF
       WRITE(6,*) ' IAMPACK in EXPCIV ', IAMPACK
*
        ITTSS_ORD = 2
        CALL EXPCIVS(LLUIN,WORK(KLBLI),NCOMBI,WORK(KLTTSNI),
     &       WORK(KLTTSII),
     &       WORK(KLIABI),NOCTPA,NOCTPB,WORK(KLBLIN),
     &       LUUT,WORK(KLBLU),NCOMBU,WORK(KLTTSNU),
     &       WORK(KLTTSIU),WORK(KLIABU),
     &       WORK(KLBLUT),
     &       ICISTR,IDC,NSMST,
     &       LBLK,IAMPACK,ITTSS_ORD)
*
      END DO
*
      IF(ICOPY.NE.0) THEN
*. Copy expanded vectors to LUIN
        CALL REWINE(LUIN,-1)
        CALL REWINE(LUUT,-1)
        DO IROOT = 1, NROOT
          IF(ICISTR.EQ.1) THEN
            CALL FRMDSC_LUCI(WORK(KLBLU),NCOMBU,-1,LUUT,IMZERO,IAMPACK)
            CALL  TODSC_LUCI(WORK(KLBLU),NCOMBU,-1,LUIN)
          ELSE
            CALL COPVCD(LUUT,LUIN,WORK(KLBLU),0,-1)
          END IF
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output  vectors in EXPCIV '
*
        CALL REWINE(LUUT,-1)
        DO IROOT = 1, NROOT
        WRITE(6,*) ' Root number ', IROOT
*
          IF(ICISTR.EQ.1) THEN
            CALL WRTMT_LU(WORK(KLBLU),1,NCOMBUT,1,NCOMBUT)
          ELSE
            CALL WRTVCD(WORK(KLBLU),LUUT,0,-1)
          END IF
        END DO
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EXPCV2')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE EXPCIVS(LUIN,VECIN,NCOMBIN,NTTSIN,ITTSIN,IABIN,
     &                   NOCTPA,NOCTPB,IBLTPIN,
     &                   LUUT,VECUT,NCOMBUT,NTTSUT,ITTSUT,IABUT,
     &                   IBLTPUT,
     &                   ICISTR,IDC,NSMST,LBLK,IAMPACKED_IN)
*
* Obtain those part of vector in cispace UT ,
* that can be obtained from terms in cispace IN
*
* Input vector on LUIN, Output vector in LUUT
* Output vector is supposed on start of vector
*
* LUIN is assumed to be single vector file,
* so rewinding will place vector on start of vector
*
* Both files are assumed on start of vector
*
* Jeppe Olsen, February 1994
*
      IMPLICIT REAL*8 (A-H,O-Z)
*. Input
      INTEGER IABIN(NOCTPA,NOCTPB),IABUT(NOCTPA,NOCTPB)
      INTEGER IBLTPIN(NSMST),IBLTPUT(NSMST)
      INTEGER NTTSIN(NOCTPA,NOCTPB,NSMST),ITTSIN(NOCTPA,NOCTPB,NSMST)
      INTEGER NTTSUT(NOCTPA,NOCTPB,NSMST),ITTSUT(NOCTPA,NOCTPB,NSMST)
      DIMENSION VECIN(*)
C     DIMENSION IATPUI(*),IBTPUI(*)
*. Output
      DIMENSION VECUT(*)
*
C     WRITE(6,*) ' EXPCIVS in action'
C     WRITE(6,*) '        Number of input  parameters',NCOMBIN
C     WRITE(6,*) '        Number of output parameters',NCOMBUT
*
      IF(ICISTR.EQ.1) THEN
        CALL FRMDSC_LUCI(VECIN,NCOMBIN,-1,LUIN,IMZERO,IAMPACK)
C       WRITE(6,*) ' Input vector '
C       CALL WRTMT_LU(VECIN,1,NCOMBIN,1,NCOMBIN)
      END IF
*
*. Loop over TTS blocks of output vector
*
*. Changed for new ordering of blocks
      IATPIN = 1
      IBTPIN = 1
      IASMIN = 0
*
      IATPUT = 1
      IBTPUT = 1
      IASMUT = 0
*
 1000 CONTINUE
*. Next output block
        CALL NXTBLK(IATPUT,IBTPUT,IASMUT,NOCTPA,NOCTPB,NSMST,
     &              IBLTPUT,IDC,NONEWUT,IABUT,ITTSS_ORD)
*. Corresponding input TTS block
        JATPIN = IATPUT
        JBTPIN = IBTPUT
        IF(IABIN(JATPIN,JBTPIN).EQ.0) THEN
          IZERO = 1
        ELSE
          IZERO = 0
        END IF
*
        NELMNT = NTTSUT(IATPUT,IBTPUT,IASMUT)
        IF(ICISTR.EQ.1) THEN
          IF(IZERO.EQ.0)
     &    IOFFIN = ITTSIN(JATPIN,JBTPIN,IASMUT)
          IOFFUT = ITTSUT(IATPUT,IBTPUT,IASMUT)
          IF(IZERO.EQ.0) THEN
            CALL COPVEC(VECIN(IOFFIN),VECUT(IOFFUT),NELMNT)
          ELSE
            ZERO = 0.0D0
            CALL SETVEC(VECUT(IOFFUT),ZERO,NELMNT)
          END IF
        END IF
*
        IF(ICISTR.NE.1.AND.NONEWUT.EQ.0) THEN
          IF(IZERO.EQ.1) THEN
            CALL ITODS(NELMNT,1,-1,LUUT)
            CALL ZERORC(-1,LUUT,IAMPACKED_IN)
          ELSE
*. Obtain input block
  999      CONTINUE
           CALL NXTBLK(IATPIN,IBTPIN,IASMIN,NOCTPA,NOCTPB,NSMST,
     &                 IBLTPIN,IDC,NONEWIN,IABIN,ITTSS_ORD)
           IF(NONEWIN.NE.0) THEN
             CALL REWINE(LUIN,-1)
*. Changed for new ordering
             IATPIN = 1
             IBTPIN = 1
             IASMIN = 0
             GOTO 999
           ELSE IF (NONEWIN.EQ. 0 ) THEN
             CALL IFRMDS(LENGTH,1,-1,LUIN)
             CALL FRMDSC_LUCI(VECIN,LENGTH,-1,LUIN,IMZERO,IAMPACK)
             IF(IATPIN.EQ.JATPIN.AND.IBTPIN.EQ.JBTPIN.AND.
     &          IASMIN.EQ.IASMUT) THEN
*. Correct block, save it
                  CALL ITODS(LENGTH,1,-1,LUUT)
                  IF(IMZERO.EQ.1) THEN
                    CALL ZERORC(-1,LUUT,IAMPACKED_IN)
                  ELSE
                    IF(IAMPACK.EQ.0) THEN
                      CALL TODSC_LUCI(VECIN,LENGTH,-1,LUUT)
                    ELSE
                      CALL TODSCP(VECIN,LENGTH,-1,LUUT)
                    END IF
                  END IF
             ELSE
               GOTO 999
             END IF
           END IF
          END IF
        END IF
      IF(NONEWUT.EQ.0) GOTO 1000
*. End of loop over output blocks
      IF(ICISTR.EQ.1) THEN
        CALL TODSC_LUCI(VECUT,NCOMBUT,-1,LUUT)
      ELSE
        CALL ITODS(-1,1,-1,LUUT)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' EXPCIVS Speaking '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' ============ '
        WRITE(6,*) ' Input Vector '
        WRITE(6,*) ' ============ '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMT_LU(VECIN,1,NCOMBIN,1,NCOMBIN)
        ELSE
          CALL WRTVCD(VECIN,LUIN,1,LBLK)
        END IF
        WRITE(6,*)
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' Output Vector '
        WRITE(6,*) ' =============== '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMT_LU(VECUT,1,NCOMBUT,1,NCOMBUT)
        ELSE
          CALL WRTVCD(VECUT,LUUT,1,LBLK)
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE EXTRROW(INMAT,IROW,NROW,NCOL,IOUTVEC)
*
* Extract row IROW from integer matrix INMAT
*
* Jeppe Olsen, Winter 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION INMAT(NROW,NCOL)
      DIMENSION IOUTVEC(NCOL)
*
      DO ICOL = 1, NCOL
        IOUTVEC(ICOL) = INMAT(IROW,ICOL)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output vector from EXTRROW '
        WRITE(6,*) ' Extracted ROW ', IROW
        CALL IWRTMA(IOUTVEC,1,NCOL,1,NCOL)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE EXTRROW2(INMAT,IROW,ICOLOFF,NROW,NCOL,IOUTVEC)
*
* Extract row IROW from integer matrix INMAT, starting from column
* ICOLOFF
*
* Jeppe Olsen, Winter 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION INMAT(NROW,ICOLOFF-1+NCOL)
      DIMENSION IOUTVEC(NCOL)
*
      DO ICOL = ICOLOFF, ICOLOFF - 1 +  NCOL
        IOUTVEC(ICOL-ICOLOFF+1) = INMAT(IROW,ICOL)
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output vector from EXTRROW '
        WRITE(6,*) ' Extracted ROW ', IROW
        CALL IWRTMA(IOUTVEC,1,NCOL,1,NCOL)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GASANA(C,NBLOCK,IBLOCK,IBLTP,LUC,ICISTR)
*
*
*
* Analyze CI vector
*
* Jeppe Olsen, August 1995
*
* Driven By IBLOCK, May 1997
*           String occupations added, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KNCPMT, KWCPMT, KLASTR, KLBSTR, KLOCCLS
!               for addressing of WORK
* =====
*.Input
* =====
*
#include "priunit.h"
#include "mxpdim.inc"
#include "orbinp.inc"
#include "strbas.inc"
#include "cicisp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
*
      DIMENSION IBLOCK(8,NBLOCK),IBLTP(*)
      CALL QENTER('ANACI')
      CALL MEMMAN(KLOFF,DUMMY,'MARK  ',DUMMY,'GASANA')
*
** Specifications of internal space
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
* Type of alpha and beta strings
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      IF(NTEST.GE.10) THEN
        WRITE(lupri,*) ' ================'
        WRITE(lupri,*) ' GASANA speaking '
        WRITE(lupri,*) ' ================'
        WRITE(lupri,*) ' IATP IBTP NAEL NBEL '
        WRITE(lupri,*)   IATP,IBTP,NAEL,NBEL
        WRITE(lupri,*) ' NOCTPA NOCTPB ', NOCTPA,NOCTPB
      END IF
*
**. Info on block structure of space
*
*
      IF( ICISTR .NE. 1 ) CALL REWINE(LUC,-1)
*. Number of terms to be printed
      IF (IPRNCIV.EQ.0) THEN
        THRES  = 1.0D0
        MAXTRM =  100
      ELSE IF (IPRNCIV.EQ.1) THEN
        THRES  = 5.0D-2
        MAXTRM =  250
      ELSE IF (IPRNCIV.EQ.2) THEN
        THRES  = 5.0D-3
        MAXTRM = 1000
      END IF
*
      IUSLAB = 0
*. Number of occupation classes
      IWAY = 1
      NEL = NAEL + NBEL
      CALL OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),
     &           0,0)
*. and then the occupation classes
      CALL MEMMAN(KLOCCLS,NGAS*NOCCLS,'ADDL  ',1,'KLOCCL')
      IWAY = 2
      CALL OCCLS(IWAY,NOCCLS,WORK(KLOCCLS),NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),
     &           0,0)
*
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
*
      LENGTH = NOCCLS*10
      CALL MEMMAN(KNCPMT,LENGTH,'ADDL  ',1,'KNCPMT')
      CALL MEMMAN(KWCPMT,LENGTH,'ADDL  ',2,'KWCPMT')
*
*. Occupation of strings of given sym and supergroup
      CALL GASANAS(C,LUC,WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &            NOCTPA,NOCTPB,MXPNGAS,IOCTPA,IOCTPB,
     &            NBLOCK,IBLOCK,
     &            THRES,MAXTRM,NAEL,NBEL,
     &            WORK(KLASTR),WORK(KLBSTR),
     &            IBLTP,NSMST,IUSLAB,
     &            IDUMMY,WORK(KNCPMT),WORK(KWCPMT),NELFSPGP,
     &            NOCCLS,NGAS,WORK(KLOCCLS),ICISTR,NTOOB,IPRNCIV)

      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'GASANA')
      CALL QEXIT('ANACI')
*
C?    Call Abend2( ' Jeppe forced me to stop after INTANA ' )
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GASANAS(C,LUC,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                   MXPNGAS,IOCTPA,IOCTPB,NBLOCK,IBLOCK,
     &                   THRES,MAXTRM,NAEL,NBEL,
     &                   IASTR,IBSTR,IBLTP,NSMST,IUSLAB,
     &                   IOBLAB,NCPMT,WCPMT,NELFSPGP,NOCCLS,NGAS,
     &                   IOCCLS,ICISTR,NORB,IPRNCIV)
*
* Analyze CI vector :
*
*      1) Print atmost MAXTRM  combinations with coefficients
*         larger than THRES
*         Currently the corresponding dets are not GIVEN !!
*
*      2) Number of coefficients in given range
*
*      3) Number of coefficients in given range for given
*         occupation class
*
* Jeppe Olsen , Jan. 1989 ,
*               Aug 1995 : GAS version
*               May 1997 : BLOCK driven
*

*. If IUSLAB  differs from zero Character*6 array IOBLAB is used to identify
*  Orbitals
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
      DIMENSION C(*)
*. General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION IBLTP(*)
      CHARACTER*6 IOBLAB(*)
      DIMENSION NELFSPGP(MXPNGAS,*)
      DIMENSION IOCCLS(NGAS,NOCCLS)
*. Specific input
      DIMENSION IBLOCK(8,NBLOCK)
*. Output
      DIMENSION NCPMT(10,NOCCLS)
      DIMENSION WCPMT(10,NOCCLS)
*
*
      IF(IUSLAB.NE.0) THEN
       WRITE(lupri,*)
       WRITE(lupri,*)
     & ' Labels for orbitals are of the type n l ml starting with n = 1'
       WRITE(lupri,*)
     & ' so the user should not be  alarmed by labels like 1 f+3 '
       WRITE(lupri,*)
      END IF

C     WRITE(6,*) 'C(1) = ',C(1)
      MINPRT = 0
      ITRM = 0
      IDET = 0
      IIDET = 0
      ILOOP = 0
      NCIVAR = 0
      IF(THRES .LT. 0.0D0 ) THRES = ABS(THRES)
      CNORM = 0.0D0
2001  CONTINUE
      IF( ICISTR .GE. 2 ) CALL REWINE(LUC,-1)
      IIDET = 0
      ILOOP = ILOOP + 1
      IF ( ILOOP  .EQ. 1 ) THEN
        XMAX = 1.0D0
        XMIN = 1.0D0/SQRT(10.0D0)
      ELSE
        XMAX = XMIN
        XMIN = XMIN/SQRT(10.0D0)
      END IF
      IF(XMIN .LT. THRES  ) XMIN =  THRES
      IF(IPRNCIV.GE.2) THEN
*. Print in one shot
       XMAX = 3006.1956
       XMIN = 0.0D0
      END IF
      IDET     = 0
      ipot_max = 5
C
      WRITE(lupri,'(/A,E10.4,A,E10.4)')
     &'   Printout of coefficients in interval  ',XMIN,' to ',XMAX
      WRITE(lupri,'(A )')
     &'   ========================================================='
*
      iprint_hline = 0
C?    WRITE(6,*) ' GASANAS : NBLOCK = ',NBLOCK
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
*. Obtain Beta  strings of sym IBSM and type IBTP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                           NORB,0,IDUM,IDUM)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
*
        IMZERO = 0
        IF( ICISTR.GE.2 ) THEN
*. Read in a Type-Type-symmetry block
          CALL IFRMDS(IDET,1,-1,LUC)
          CALL FRMDSC_LUCI(C,IDET,-1,LUC,IMZERO,IAMPACK)
C?        WRITE(6,*) ' Number of elements readin ',IDET
          IDET = 0
        END IF
        IF(IMZERO.NE.1) THEN
*
        IBBAS = 1
        IABAS = 1
*

        DO  IB = IBBAS,IBBAS+NIB-1
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB - IBBAS + IABAS
          ELSE
            MINIA = IABAS
          END IF
          DO  IA = MINIA,IABAS+NIA-1
*
            IF(ILOOP .EQ. 1 ) NCIVAR = NCIVAR + 1
            IDET = IDET + 1
            IF( XMAX .GE. ABS(C(IDET)) .AND.
     &      ABS(C(IDET)).GT. XMIN ) THEN
              ITRM = ITRM + 1
              IIDET = IIDET + 1
              IF( ITRM .LE. MAXTRM ) THEN
                CNORM = CNORM + C(IDET) ** 2
                if(iprint_hline.gt.0)then
                  WRITE(lupri,'(A )')
     &            '                  =================== '
                  iprint_hline = iprint_hline + 1
                end if
*
                WRITE(lupri,'(A,I8,A,E14.8)')
     &          '   Coefficient of combination ',IDET,' is ',
     &          C(IDET)
                WRITE(lupri,'(A)')
     &          '   Corresponding alpha - and beta string '
                IF(IUSLAB.EQ.0) THEN
                  WRITE(lupri,'(4X,10I4)')
     &            (IASTR(IEL,IA),IEL = 1, NAEL )
                  WRITE(lupri,'(4X,10I4)')
     &            (IBSTR(IEL,IB),IEL = 1, NBEL )
                ELSE
                  WRITE(lupri,'(4X,10(1X,A6))')
     &            (IOBLAB(IASTR(IEL,IA)),IEL = 1, NAEL )
                  WRITE(lupri,'(4X,10(1X,A6))')
     &            (IOBLAB(IBSTR(IEL,IB)),IEL = 1, NBEL )
                END IF
              END IF
            END IF
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over beta strings
        END IF
*       ^ End of if statement for nonvanishing blocks
        END DO
*       ^ End of loop over blocks
       IF(IIDET .EQ. 0 ) WRITE(lupri,*) '    ( no coefficients )'
       IF( XMIN .GT. THRES .AND. ILOOP .LE. 30 ) GOTO 2001
*
       WRITE(lupri,'(A,E15.8)')
     & '   Norm of printed CI vector .. ', CNORM
*
*.Size of CI coefficients
*
*
      IDET = 0
      IF(ICISTR .GE. 2 ) CALL REWINE(LUC,-1)
      CALL ISETVC(NCPMT,0    ,10*NOCCLS)
      CALL SETVEC(WCPMT,0.0D0,10*NOCCLS)
      DO JBLOCK = 1, NBLOCK
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
*
        IF(IBLTP(IASM).EQ.2) THEN
          IRESTR = 1
        ELSE
          IRESTR = 0
        END IF
*. Occupation class corresponding to given occupation
        JOCCLS = 0
        DO JJOCCLS = 1, NOCCLS
          IM_THE_ONE = 1
          DO IGAS = 1, NGAS
            IF(NELFSPGP(IGAS,IATP-1+IOCTPA)+
     &         NELFSPGP(IGAS,IBTP-1+IOCTPB).NE.IOCCLS(IGAS,JJOCCLS))
     &         IM_THE_ONE = 0
          END DO
          IF(IM_THE_ONE .EQ. 1 ) JOCCLS = JJOCCLS
        END DO
*
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
*
        IMZERO = 0
        IF( ICISTR.GE.2 ) THEN
*. Read in a Type-Type-symmetry block
          CALL IFRMDS(IDET,1,-1,LUC)
          CALL FRMDSC_LUCI(C,IDET,-1,LUC,IMZERO,IAMPACK)
          IDET = 0
        END IF
        IF(IMZERO.EQ.0) THEN

        DO IB = IBBAS,IBBAS+NIB-1
          IF(IRESTR.EQ.1.AND.IATP.EQ.IBTP) THEN
            MINIA = IB - IBBAS + IABAS
          ELSE
            MINIA = IABAS
          END IF
          DO IA = MINIA,IABAS+NIA-1
            IDET = IDET + 1
            DO IPOT = 1, ipot_max
              IF(10.0D0 ** (-IPOT+1).GE.ABS(C(IDET)).AND.
     &           ABS(C(IDET)).GT. 10.0D0 ** ( - IPOT )) THEN
                 NCPMT(IPOT,JOCCLS)= NCPMT(IPOT,JOCCLS)+ 1
                 WCPMT(IPOT,JOCCLS)= WCPMT(IPOT,JOCCLS)+ C(IDET) ** 2
              END IF
            END DO
*           ^ End of loop over powers of ten
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over beta strings
        END IF
*       ^ End of test for novanishing blocks
      END DO
*     ^ End of loop over blocks
*
      WRITE(lupri,'(/A)') '     Magnitude of CI coefficients '
      WRITE(lupri,'(A/)') '    =============================='
      WACC = 0.0D0
      NACC = 0
      DO 300 IPOT = 1, ipot_max
        W = 0.0D0
        N = 0
        DO 290 JOCCLS = 1, NOCCLS
            N = N + NCPMT(IPOT,JOCCLS)
            W = W + WCPMT(IPOT,JOCCLS)
  290   CONTINUE
        WACC = WACC + W
        NACC = NACC + N
        IF (N .GT. 0) WRITE(lupri,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &  '   10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC
  300 CONTINUE
*
      WRITE(lupri,'(/4x,a,i9/)') 
     & 'Number of coefficients less than 10^{-11} is ',NCIVAR-NACC
*
      IF(NOCCLS.NE.1) THEN
      WRITE(lupri,'(A)')
     & '    Magnitude of CI coefficients for each excitation level '
      WRITE(lupri,'(A)')
     & '   ========================================================='
      DO 400 JOCCLS = 1, NOCCLS
          N = 0
          DO 380 IPOT = 1, ipot_max
            N = N + NCPMT(IPOT,JOCCLS)
  380     CONTINUE
          IF(N .NE. 0 ) THEN
            WRITE(lupri,'(/7X,A,15I3)')' Occupation of active sets :',
     &      (IOCCLS(IGAS,JOCCLS),IGAS=1, NGAS)
            WRITE(lupri,'(A,I9)')
     &      '        Number of coefficients larger than 10^-{11}:', N
            WRITE(lupri,*)
            WACC = 0.0D0
            DO 370 IPOT = 1, ipot_max
              N =  NCPMT(IPOT,JOCCLS)
              W =  WCPMT(IPOT,JOCCLS)
              WACC = WACC + W
              IF (N .GT. 0)
     &           WRITE(lupri,'(A,I2,A,I2,3X,I9,1X,E15.8,3X,E15.8)')
     &           '   10-',IPOT,' TO 10-',(IPOT-1),N,W,WACC
  370       CONTINUE
          END IF
  400 CONTINUE
*
*. Total weight and number of dets per excitation level
*
      WRITE(lupri,'(/A)')
     & '    Total weight and number of SD''s (> 10^{-11}):  '
      WRITE(lupri,'(A/A/A)')
     &  '  =====================================================',
     &  '         N        Weight    Acc. Weight   Occupation',
     &  '  ====================================================='
      WACC = 0.0D0
      DO 500 JOCCLS = 1, NOCCLS
          N = 0
          W = 0.0D0
          DO 480 IPOT = 1, ipot_max
            N = N + NCPMT(IPOT,JOCCLS)
            W = W + WCPMT(IPOT,JOCCLS)
  480     CONTINUE
          WACC = WACC + W
          IF(N .NE. 0 ) THEN
            WRITE(lupri,'(I11,F14.8,F14.8,2X,16I3)')
     &      N,W,WACC,(IOCCLS(IGAS,JOCCLS),IGAS=1,NGAS)
          END IF
  500 CONTINUE
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_BLOCKS_FROM_DISC
     &           (LU,NBLOCK,IOFF,IBLOCK,NBLOCKT,C,IREW)
*
* Obtain blocks IOFF - IOFF + NBLOCK from file LU
*
* Jeppe Olsen, January 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(*)
      DIMENSION IBLOCK(8,*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "csm.inc"
#include "stinf.inc"
*
      PARAMETER (IATP = 1, IBTP = 2)
*
#ifdef NEW_OR_OLDis2
*. This is "old" version, activate with #define NEW_OR_OLDis2
*. Local Scratch
*
      PARAMETER (MXPNBLK = 33000)
      DIMENSION ISCR(3*MXPNBLK)
#endif
*
C?    WRITE(6,*) ' ENTERING GET_BLOCKS_FROM_DISC'
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
CM    IREW = 1
      ISCAL = 0
#ifdef NEW_OR_OLDis2
      IF(NBLOCKT.GT.MXPNBLK) THEN
         WRITE(6,*) ' Increase parameter MXPNBLK in'
         WRITE(6,*) ' GET_BLOCKS_FROM_DISC '
         WRITE(6,*) ' Current and required values ',
     &     MXPNBLK,NBLOCKT
         Call Abend2( 'GET_BLOCKS_FROM_DISC' )
      END IF
      CALL GET_TTS_BATCH(C,NBLOCK,IBLOCK(1,IOFF),NBLOCKT,IBLOCK,
     &                   NOCTPA,NOCTPB,NSMST,
     &                   WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &                   IDC,LU,ISCR,IREW,ISCAL)
#else
      CALL GET_TTS_BATCHN(C,NBLOCK,IBLOCK(1,IOFF),NBLOCKT,IBLOCK,
     &                    NOCTPA,NOCTPB,NSMST,
     &                    WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &                    IDC,LU,IREW,ISCAL)
#endif
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_DIAG_DOM(AIN,AOUT,NDIM,IREO)
*
* A square matrix AIN is given. Permute columns to obtain
* Largest elements on diagonal
*
* Jeppe Olsen, Feb.98
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      DIMENSION AIN(NDIM,NDIM)
*.Output
      DIMENSION AOUT(NDIM,NDIM)
      INTEGER IREO(*)
*
*. IREO : Old => New
      IZERO = 0
      CALL ISETVC(IREO,IZERO,NDIM)
      DO ICOL = 1, NDIM
*. Column with largest values at row ICOL
        XVAL = 0.0D0
        IICOL = 0
        DO JCOL = 1, NDIM
          IF(ABS(AIN(ICOL,JCOL)).GT.XVAL.AND.IREO(JCOL).EQ.0) THEN
            IICOL = JCOL
            XVAL = ABS(AIN(ICOL,JCOL))
          END IF
        END DO
        IF(IICOL.GT.0) THEN
          CALL COPVEC(AIN(1,IICOL),AOUT(1,ICOL),NDIM)
          IREO(IICOL) = ICOL
        ELSE
          CALL COPVEC(AIN(1,ICOL),AOUT(1,ICOL),NDIM)
          IREO(ICOL) = ICOL
        END IF
      END DO
*
      NTEST = 10
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
        WRITE(6,*) ' OUTPUT FROM GET_DIAG_DOM '
        WRITE(6,*) ' ======================== '
        WRITE(6,*)
        WRITE(6,*) ' IREO array : old => new order '
        CALL IWRTMA(IREO,1,NDIM,1,NDIM)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reordered matrix '
        CALL WRTMT_LU(AOUT,NDIM,NDIM,NDIM,NDIM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_NMLW(VEC,N,M,LW,N_MAX,M_MAX,
     &           LUCORR,LUREF,IREAD,LUOUT,IREW,LBLK,LUOUTEFF)
*
* Obtain correction vector ! N,M,LW> obtained in GNDBTFREQ
* and save on file LUOUT
* If IREAD = 0, the file is positioned
* for reading the file, but no transfer is realized
* If IREW ne 0, LUOUT is rewinded
*
*. Jeppe in Pisa, July 17, 97
*                 Sept. 97, Changed to disc, ICISTR = 2,3
*
*
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)' GET_NMLW : Vector in Demand, N,M,LW :', N,M,LW
      END IF
*
      LUOUTEFF = LUCORR
      IF(IREW.NE.0) CALL REWINE(LUOUT,-1)
*
      IF(N.EQ.0.AND.M.EQ.0) THEN
*. You asked for the reference vector and that will be
        LUOUTEFF = LUREF
        CALL REWINE(LUREF,-1)
        IVEC = 1
        IF(IREAD.NE.0 )  THEN
          CALL COPVCDP(LUREF,LUOUT,VEC,1,LBLK)
        END IF
      ELSE
*. A real correction vector
*. Position to start of requested vector
        CALL REWINE(LUCORR,-1)
        IVEC = 0
        DO NP = 0, N
          IF(NP.EQ.N) THEN
           MXMP = M
          ELSE
           MXMP = M_MAX
          END IF
          DO  MP = 0, MXMP
           IF(NP.EQ.N.AND.MP.EQ.M)THEN
             MXLWP = LW-1
           ELSE
             MXLWP = NP
           END IF
           DO LWP = -NP, MXLWP,2
C?           WRITE(6,*) ' Next vector NP,MP,LWP',NP,MP,LWP
*. Skip
             IF(.NOT.(NP.EQ.0.AND.MP.EQ.0)) THEN
C?            WRITE(6,*) ' SKPVCD will be called '
              CALL SKPVCD(LUCORR,1,VEC,0,LBLK)
              IVEC = IVEC + 1
C?            write(6,*) ' Vector skipped ', IVEC
             END IF
           END DO
          END DO
        END DO
*. And then we are ready to read the vector of interest
        IF(IREAD.NE.0) THEN
          CALL COPVCDP(LUCORR,LUOUT,VEC,0,LBLK)
          IVEC = IVEC + 1
        END IF
      END IF
*
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' Vector to be read '
        WRITE(6,*) ' N M LW =',N,M,LW
        WRITE(6,*) ' obtained as vector ',IVEC
      END IF
*
      IF(NTEST.GE.100) THEN
       IF(IREAD.NE.0) THEN
         IF(IREW.NE.0) THEN
           WRITE(6,*) ' Vector read in '
           CALL WRTVCD(VEC,LUOUT,1,LBLK)
         ELSE
           WRITE(6,*) ' Vector is not first vector'
           WRITE(6,*) ' so it will not be printed '
         END IF
       ELSE
         WRITE(6,*) 'vector positioned at end of vector',IVEC
       END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_TTS_BATCH(CTTS,NBLOCK,IBLOCK,NBLOCKC,IBLOCKC,
     &                  NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                  IDC,LUC,ISCR,IREW,ISCALE)
*
* Read in a batch of C blocks from file LUC.
*
* The complete file is defined by NBLOCKC,IBLOCKC,
* and the blocks of the actual batch is defined by NBLOCK,IBLOCK.
* Vector packing is defined by IDC
*
* Should be initialized with rewind on LUC i.e. IREW = 1
*
*. Feb 96 : Modified to accomodate packed files
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CTTS(*),NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IBLOCKC(8,NBLOCKC),IBLOCK(8,NBLOCK)
*
      COMMON/HIDLUC/IBLK,IATP,IBTP,IASM,IBSM
*. Local scratch : should atleast be 3*NBLOCK
      DIMENSION ISCR(*)
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to GET_TTS_BATCH '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' Number of blocks in batch ', NBLOCK
        WRITE(6,*)
        WRITE(6,*) ' IATP  IBTP  IASM  IBSM '
        WRITE(6,*) ' ====================== '
        DO JBLOCK = 1, NBLOCK
          WRITE(6,'(4I4)') (IBLOCK(II,JBLOCK),II = 1, 4 )
        END DO
      END IF
*

*
      IF(IREW.EQ.1) THEN
        REWIND (LUC)
        IBLK = 1
        IATP = IBLOCKC(1,1)
        IBTP = IBLOCKC(2,1)
        IASM = IBLOCKC(3,1)
        IBSM = IBLOCKC(4,1)
      END IF
*
*. Map batch blocks to global blocks
*
      DO JBLOCK = 1, NBLOCK
        DO JBLOCKC = 1, NBLOCKC
          IF(IBLOCK(1,JBLOCK).EQ.IBLOCKC(1,JBLOCKC) .AND.
     &       IBLOCK(2,JBLOCK).EQ.IBLOCKC(2,JBLOCKC) .AND.
     &       IBLOCK(3,JBLOCK).EQ.IBLOCKC(3,JBLOCKC) .AND.
     &       IBLOCK(4,JBLOCK).EQ.IBLOCKC(4,JBLOCKC)       ) THEN
           ISCR(JBLOCK) = JBLOCKC
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Blocks mapped to global order '
        CALL IWRTMA(ISCR,1,NBLOCK,1,NBLOCK)
      END IF
*. Order
C     ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
*
* ORDER A STRING OF INTEGERS TO ASCENDING ORDER
*
* IINST : INPUT STRING
* IOUTST : OUTPUT STRING
* NELMNT : NUMBER OF INTEGERS
* INO : Mapping array from new to old order
      KLORIG = 1
      KLORDER = KLORIG + NBLOCK
      KLNTOO = KLORDER + NBLOCK
      CALL ORDINT(ISCR(KLORIG),ISCR(KLORDER),NBLOCK,ISCR(KLNTOO),NTEST)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Blocks ordered and original place '
        CALL IWRTMA(ISCR(KLORDER),1,NBLOCK,1,NBLOCK)
        CALL IWRTMA(ISCR(KLNTOO),1,NBLOCK,1,NBLOCK)
      END IF
*
*. Ready to loop over blocks
*
C     IOFFO = 1
      DO JJBLK = 1, NBLOCK
        JBLK = ISCR(KLORDER-1+JJBLK)
        IOFFO   = IBLOCK(6,ISCR(KLNTOO-1+JJBLK))
        IBSM = IBLOCKC(4,1)
*. Loop to this block
        IPACK = 1
        IMZERO= 0
        DO KBLK = IBLK,JBLK-1
          CALL IFRMDS(LBL,1,-1,LUC)
          CALL SKPRCD2(LBL,-1,LUC)
C         CALL IFRMDS(LBL,1,-1,LUC)
C         IF(IPACK.EQ.1) CALL IFRMDS(IMZERO,1,-1,LUC)
C         IF(IPACK.EQ.0.OR.IMZERO.EQ.0) THEN
C           READ(LUC)
C         END IF
        END DO
*. and : READIN !!!
        CALL IFRMDS(LBL,1,-1,LUC)
        CALL FRMDSC_LUCI(CTTS(IOFFO),LBL,-1,LUC,IMZERO,IAMPACK)
C?      WRITE(6,*) ' GET_TTS, JJBLK,IOFFO LBL ', JJBLK,IOFFO,LBL
C       IOFFO = IOFFO + LBL
        IBLK = JBLK+1
      END DO
*
*. Rescale from combination form to determinant
*
C     ISCALE = 0
      IF(IDC.EQ.2.AND.ISCALE.NE.0) THEN
        CALL SCDTTS(CTTS,IBLOCK,NBLOCK,
     &  NSMST,NICTOA,NOCTPB,NSASO,NSBSO,
     &  IDC,2,NTEST)
C     SCDTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,IDC,IWAY,IPRNT)
      END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Batch of blocks read in '
         WRITE(6,*) ' ========================'
         WRITE(6,*)
         CALL WRTTTS(CTTS,IBLOCK,NBLOCK,NSMST,NOCTPA,NOCTPB,
     &               NSASO,NSBSO,IDC)
      END IF
C     WRTTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,ISC)

      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_TTS_BLK_IN_VECTOR
     &           (NBLOCKI,IBLOCKI,IOFFI,VECI,
     &            NBLOCKO,IBLOCKO,IOFFO,
     &            NBLOCKIO,LBLOCKIO,IBLOCKIO,VECIO )
* A vector VECI containing NBLOCKI blocks defined by IBLOCKI
* are given. Obtain those blocks in VECI that corresponds to
* IBLOCKO(1,IOFFO) - IBLOCKO(1,IOFFO+NBLOCKO-1)
*
*. The number of common blocks is NBLOCKIO
*. The number of elements in these blocks is LBLOCKIO
*. Save the common blocks in VECIO and the corresponding block info in
*  IBLOCKIO
*
*
* Jeppe Olsen,  November 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBLOCKI(8,*)
      INTEGER IBLOCKO(8,*)
      DIMENSION VECI(*)
*. output
      INTEGER IBLOCKIO(8,*)
      DIMENSION VECIO(*)
*. Output
*
      IOFFP = 1
      IOFF = 1
      NBLOCKIO = 0
      LBLOCKIO = 0
      DO JBLOCKI = IOFFI,IOFFI+NBLOCKI-1
        DO JBLOCKO = IOFFO, IOFFO + NBLOCKO-1
          IF(IBLOCKI(1,JBLOCKI).EQ. IBLOCKO(1,JBLOCKO).AND.
     &    IBLOCKI(2,JBLOCKI).EQ. IBLOCKO(2,JBLOCKO)   .AND.
     &    IBLOCKI(3,JBLOCKI).EQ. IBLOCKO(3,JBLOCKO)   .AND.
     &    IBLOCKI(4,JBLOCKI).EQ. IBLOCKO(4,JBLOCKO)          ) THEN
            NBLOCKIO = NBLOCKIO + 1
            CALL ICOPVE(IBLOCKI(1,JBLOCKI),IBLOCKIO(1,NBLOCKIO),8)
            IBLOCKIO(5,NBLOCKIO) = IOFF
            IBLOCKIO(6,NBLOCKIO) = IOFFP
            CALL COPVEC(VECI(IBLOCKI(6,JBLOCKI)),VECIO(IOFFP),
     &                  IBLOCKI(8,JBLOCKI))
            IOFF = IOFF + IBLOCKI(7,JBLOCKI)
            IOFFP= IOFFP+ IBLOCKI(8,JBLOCKI)
          END IF
        END DO
      END DO
      LBLOCKIO = IOFFP - 1
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
     &   ' Number of blocks obtained in GET_TTS_FROM_VECTOR ',NBLOCKIO
        WRITE(6,*)
     &   ' Number of elements obtained ',LBLOCKIO
        WRITE(6,*) ' Blocks obtained : '
        CALL IWRTMA(IBLOCKIO,8,NBLOCKIO,8,NBLOCKIO)
        WRITE(6,*) ' corresponding vector '
        NELMNTIO = IOFFP - 1
        CALL WRTMT_LU(VECIO,1,NELMNTIO,1,NELMNTIO)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0TV2(VEC1,VEC2,LUC,LUHC,LU0,LUSCR1,E0,LBLK)
*
* Multiply vector in LUC with H0 where H0 is defined as
*
* H0 = (1-|0><0|) H apr  (1-|0><0>) + E0 |0><0>
*
* Where H apr is defined by call to MV7 ( with IPERTOP = 1)
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
      REAL*8 INPRDD
*. For communicating with sigma routine
#include "oper.inc"
#include "mxpdim.inc"
#include "crun.inc"
*
C     WRITE(6,*) ' H0TV2 : LUC, LUHC, LU0, LUSCR1 ',
C    &                     LUC, LUHC, LU0, LUSCR1
*. Overlap <C|0>
      IF(LU0.GT.0) THEN
        SC0 = INPRDD(VEC1,VEC2,LUC,LU0,1,LBLK)
*. C -  <C|0> |0> on LUSCR1
        FAC1 = 1.0D0
        FAC2 = -SC0
        CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUC,LU0,LUSCR1,1,LBLK)
      ELSE
        CALL COPVCD(LUC,LUSCR1,VEC1,1,LBLK)
      END IF
*. Multiply with H apr, result on LUHC
*. Place for trouble in the future.
      IF(IH0INSPC(1).NE.4) THEN
         IPERTOP = 1
      END IF
      WRITE(6,*) ' MV7 will be called in a few NANOSECONDS'
      CALL MV7(VEC1,VEC2,LUSCR1,LUHC)
*. Orthogonalize LUHC to LU0
      IF(LU0.NE.0) THEN
        SSIGMA0 = INPRDD(VEC1,VEC2,LUHC,LU0,1,LBLK)
        FAC1 = 1.0D0
        FAC2 = -SSIGMA0
        CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUHC,LU0,LUSCR1,1,LBLK)
*. and add E0 <C|0> |0>
        FAC1 = 1.0D0
        FAC2 = E0 * SC0
        CALL VECSMD(VEC1,VEC2,FAC1,FAC2,LUSCR1,LU0,LUHC,1,LBLK)
      ELSE
CSEPT29 CALL COPVCD(LUSCR1,LUHC,VEC1,1,LBLK)
      END IF
*.
      NTEST = 2
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' results from H0TV2 '
        WRITE(6,*) ' ==================='
        write(6,*) ' SC0, SSIGMA0 ', SC0,SSIGMA0
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Input vector '
        CALL WRTVCD(VEC1,LUC,1,LBLK)
        WRITE(6,*)
        WRITE(6,*) ' Output vector '
        CALL WRTVCD(VEC1,LUHC,1,LBLK)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE H0TVM(VEC1,VEC2,LLUC,LLUHC)
*
* Outer routine for zero order operator + shift times vector
*
*. Input  vector : on LLUC
*. Output fector : on LLUHC
*
* Jeppe Olsen, February 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "clunit.inc"
C
*. Transfer of zero order energy
      COMMON/CENOT/E0
*. Transfer of shift
      COMMON/CSHIFT/SHIFT,IPROJ
*. Default block parameter
      LBLK = -1
*.  Zero order vector is assumed on LUC
      IF(IPROJ.EQ.0) THEN
       LU0 = 0
      ELSE IF (IPROJ.EQ.1) THEN
       LU0 = LUC
      ELSE
       WRITE(6,*)  ' H0TVM, Unknown IPROJ = ', IPROJ
       Call Abend2( ' H0TVM, Unknown IPROJ  ' )
      END IF
      LUSCR1 = LUSC40
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*)
        WRITE(6,*) '============== '
        WRITE(6,*) ' H0TVM entered '
        WRITE(6,*) '============== '
        WRITE(6,*)
        WRITE(6,*) ' LLUC LLUHC LU0 and LUSCR1 ',
     &               LLUC,LLUHC,LU0,LUSCR1
        WRITE(6,*) ' E0 , Shift : ', E0 , SHIFT
      END IF
*. A scratch file not used by linear solver in SIMPRT : LUSCR1
      IF(SHIFT.EQ.0.0D0) THEN
        CALL H0TV2(VEC1,VEC2,LLUC,LLUHC,LU0,LUSCR1,E0,LBLK)
      ELSE
*. H0TV on LUSCR1
        CALL H0TV2(VEC1,VEC2,LLUC,LUSCR1,LU0,LLUHC,E0,LBLK)
*. Add shift and save on LLUHC
        ONE = 1.0D0
        CALL VECSMD(VEC1,VEC2,ONE,SHIFT,LUSCR1,LLUC,LLUHC,1,LBLK)
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input and output vectors from H0TVM '
        CALL WRTVCD(VEC1,LLUC,1,LBLK)
        WRITE(6,*)
        CALL WRTVCD(VEC1,LLUHC,1,LBLK)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE HPQTV(NP1,NP2,NQ,PHP,PHQ,QHQ,VECIN,VECUT)
*
* Multiply P1P2Q preconditioner with a vector
* Jeppe Olsen , July 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General Input
      DIMENSION PHP(*),PHQ(NP1,NQ),QHQ(NQ)
*. PHP is in lower triangular form
*. Specific input
      DIMENSION VECIN(*)
*.Output
      DIMENSION VECUT(*)
*
      NP = NP1+NP2
      NPQ = NP + NQ
*
      call dzero(vecut,npq)
*. PHQ * VECIN
      CALL MATVCC(PHQ,VECIN(1+NP),VECUT,NP1,NQ,0)
*. PHP * VECIN
      IJ = 0
      DO 60 I = 1, NP
      DO 50 J = 1, I
        IJ = IJ + 1
        VECUT(I) = VECUT(I)+PHP(IJ)*VECIN(J)
        VECUT(J) = VECUT(J)+PHP(IJ)*VECIN(I)
   50 CONTINUE
      VECUT(I) = VECUT(I)-PHP(IJ)*VECIN(I)
   60 CONTINUE
* QHP * VECIN
      CALL MATVCC(PHQ,VECIN,VECUT(1+NP),NP1,NQ,1)
*. QHQ * VECIN
      DO 100 I = 1, NQ
        VECUT(I+NP) = VECUT(I+NP)+QHQ(I)*VECIN(NP+I)
  100 CONTINUE
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input and output vectors from HPQTV '
        CALL WRTMT_LU(VECIN,1,NPQ,1,NPQ)
        CALL WRTMT_LU(VECUT,1,NPQ,1,NPQ)
      END IF
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE HPQTVM(VECIN,VECUT,IDUM,JDUM)
*
* Outer routine for Preconditioner times vector
*
*. LUCIA version
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "glbbas.inc"
#include "crun.inc"
C     COMMON/BIGGY/WORK(MXPWRD)
#include "wrkspc.inc"
*
      COMMON/SHFT/SHIFT
*
      NPDM = MXP1+MXP2
      NQDM = MXQ
      KLPHP = KH0
      KLPHQ = KH0 + NPDM*(NPDM+1)/2
      KLQHQ = KLPHQ + MXP1*MXQ
*
C          HPQTV(NP1,NP2,NQ,PHP,PHQ,QHQ,VECIN,VECUT)
      CALL HPQTV(MXP1,MXP2,MXQ,
     &           WORK(KLPHP),WORK(KLPHQ),WORK(KLQHQ),
     &           VECIN,VECUT)
*
      IF(SHIFT.NE.0.0D0) THEN
       CALL VECSUM(VECUT,VECUT,VECIN,1.0D0,SHIFT,NPDM+NQDM)
      END IF
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,*) ' Input and output vectors from HPQTVM'
        CALL WRTMT_LU(VECIN,1,NPDM+NQDM,1,NPDM+NQDM)
        CALL WRTMT_LU(VECUT,1,NPDM+NQDM,1,NPDM+NQDM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE LULU(A,L,U,NDIM)
C
C LU DECOMPOSITION OF MATRIX A
C
C     A = L * U
C
C WHERE L IS A LOWER TRIANGULAR MATRIX WITH A
C UNIT DIAGONAL AND U IS AN UPPER DIAGONAL
C
C L AND U ARE STORED AS ONE DIMENSIONAL ARRAYS
C
C   L(I,J) = L(I*(I-1)/2 + J ) ( I .GE. J )
C
C   U(I,J) = U(J*(J-1)/2 + I ) ( J .GE. I )
C
C THIS ADRESSING SCHEMES SUPPORTS VECTORIZATION OVER COLUMNS
C FOR L AND  OVER ROWS FOR U .
C
C
C NO PIVOTING IS DONE HERE , SO THE SCHEME GOES :
C
C     LOOP OVER R=1, NDIM
C        LOOP OVER J = R, NDIM
C          U(R,J) = A(R,J) - SUM(K=1,R-1) L(R,K) * U(K,J)
C        END OF LOOP OVER J
C
C        LOOP OVER I = R+1, NDIM
C          L(I,R) = (A(I,R) - SUM(K=1,R-1)L(I,K) * U(K,R) ) /U(R,R)
C        END OF LOOP OVER I
C     END OF LOOP OVER R
C
C JEPPE OLSEN , OCTOBER 1988
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(NDIM,NDIM)
      DOUBLE PRECISION  L(*),U(*)
      REAL * 8  INPROD
      INTEGER R
C
C
      DO 1000 R = 1, NDIM
C
        DO 100 J = R, NDIM
         U(J*(J-1)/2 + R ) = A(R,J) -
     &   INPROD(L(R*(R-1)/2+1),U(J*(J-1)/2+1),R-1)
  100   CONTINUE
C
        XFACI = 1.0D0/ U(R*(R+1)/2)
        L(R*(R+1)/2 ) = 1.0D0
        DO 200 I = R+1, NDIM
          L(I*(I-1)/2 + R) = (A(I,R) -
     &   INPROD(L(I*(I-1)/2+1),U(R*(R-1)/2+1),R-1) ) * XFACI
  200  CONTINUE
C
 1000 CONTINUE
C
      NTEST = 0
      IF ( NTEST .NE. 0 ) THEN
         WRITE(6,*) ' L MATRIX '
         CALL PRSYM(L,NDIM)
         WRITE(6,*) ' U MATRIX ( TRANSPOSED ) '
         CALL PRSYM(U,NDIM)
      END IF
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MULT_BLOC_MAT(C,A,B,NBLOCK,LCROW,LCCOL,
     &                         LAROW,LACOL,LBROW,LBCOL,ITRNSP)
*
* Multiply two blocked matrices
*
* ITRNSP = 0 => C = A * B
* ITRNSP = 1 => C = A(T) * B
* ITRNSP = 2 => C = A * B(T)
*
* Jeppe Olsen
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION A(*),B(*)
      INTEGER LCROW(NBLOCK),LCCOL(NBLOCK)
      INTEGER LAROW(NBLOCK),LACOL(NBLOCK)
      INTEGER LBROW(NBLOCK),LBCOL(NBLOCK)
*. Output
      DIMENSION C(*)
*
      DO IBLOCK = 1, NBLOCK
       IF(IBLOCK.EQ.1)  THEN
         IOFFA = 1
         IOFFB = 1
         IOFFC = 1
       ELSE
         IOFFA = IOFFA + LACOL(IBLOCK-1)*LACOL(IBLOCK-1)
         IOFFB = IOFFB + LBCOL(IBLOCK-1)*LBCOL(IBLOCK-1)
         IOFFC = IOFFC + LCCOL(IBLOCK-1)*LCCOL(IBLOCK-1)
       END IF
*
       ZERO = 0.0D0
       ONE = 1.0D0
       CALL MATML7(C(IOFFC),A(IOFFA),B(IOFFB),
     &             LCROW(IBLOCK),LCCOL(IBLOCK),
     &             LAROW(IBLOCK),LACOL(IBLOCK),
     &             LBROW(IBLOCK),LBCOL(IBLOCK),
     &             ZERO,ONE,ITRNSP)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' output matrix from MULT_BLOC_MAT '
        WRITE(6,*) ' ================================='
        WRITE(6,*)
        CALL APRBLM2(C,LCROW,LCCOL,NBLOCK,0)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE PAMTMT(X,T,WORK,NORB)
*
* GENERATE PER AKE'S T MATRIX FROM AN
* ORBITAL ROTATION MATRIX X
*
* T IS OBTAINED AS A STRICTLY LOWER TRIANGULAR
* MATRIX TL AND AN UPPER TRIANGULAR MATRIX TU
*
*         TL = 1 - L
*         TU = U ** -1
*
* WHERE L AND U ARISES FROM A LU DECOMPOSITION OF
* X :
*         X = L * U
* WITH L BEING A LOWER TRIANGULAR MATRIX WITH UNIT ON THE
* DIAGONAL AND U IS AN UPPER TRIANGULAR MATRIX
*
* JEPPE OLSEN OCTOBER 1988
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NORB,NORB),T(NORB,NORB)
      DIMENSION WORK(*)
* DIMENSION OF WORK : NORB ** 2 + NORB*(NORB+1) / 2
*
      NTEST = 00
      IF(NTEST.GE.2) THEN
        WRITE(6,*) ' Wellcome to PAMTMT '
        WRITE(6,*) ' =================='
        WRITE(6,*)
      END IF
*. Allocate local memory
      KLFREE = 1
C     KLL = KFLREE
      KLL = KLFREE
      KLFREE = KLL + NORB*(NORB+1)/2
      KLU = KLFREE
      KLFREE = KLU + NORB ** 2
*.LU factorize X
      CALL LULU(X,WORK(KLL),WORK(KLU),NORB)
*.Expand U to full matrix
      CALL SETVEC(T,0.0D0,NORB ** 2 )
      DO 10 I = 1,NORB
      DO 10 J = I,NORB
        T(I,J) = WORK(KLU-1+J*(J-1)/2+I)
   10 CONTINUE
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' MATRIX TO BE INVERTED '
        CALL WRTMT_LU(T,NORB,NORB,NORB,NORB)
      END IF
*.Invert U
      CALL INVMAT(T,WORK(KLU),NORB,NORB)
      IF ( NTEST .GE. 100 ) THEN
        WRITE(6,*) ' INVERTED MATRIX '
        CALL WRTMT_LU(T,NORB,NORB,NORB,NORB)
      END IF
*.Subtract L
      DO 20 I = 1, NORB
      DO 20 J = 1,I-1
       T(I,J)= - WORK(KLL-1+I*(I-1)/2+J)
   20 CONTINUE
*
      IF( NTEST .GE. 2 ) THEN
        WRITE(6,*) ' INPUT X MATRIX '
        CALL WRTMT_LU(X,NORB,NORB,NORB,NORB)
        WRITE(6,*) ' T MATRIX '
        CALL WRTMT_LU(T,NORB,NORB,NORB,NORB)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE PR_MATML_STAT
*
* Print statistics from matrix multiplier
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/MATMLST/XNFLOP,XNCALL,XLCROW,XLCCOL,XLCROWCOL,TMULT
*
      WRITE(6,*)
      WRITE(6,*) ' Information about MATML7 calls : '
      WRITE(6,*) ' ================================'
      WRITE(6,*)
      WRITE(6,*) ' Number of calls ', XNCALL
      WRITE(6,*) ' Number of flops executed ', XNFLOP
      WRITE(6,*) ' Average row length of C ',
     &           XLCROWCOL/XLCCOL
      WRITE(6,*) ' Average column  length of C ',
     &           XLCROWCOL/XLCROW
      WRITE(6,*) ' Average number of operations per element of  C ',
     &           XNFLOP/XLCROWCOL
      WRITE(6,*) ' Average number of operations per per call ',
     &           XNFLOP/XNCALL
      WRITE(6,*) ' Number of seconds spent in MATML7', TMULT
*
      WRITE(6,*) ' Average MFLOPS ',XNFLOP/TMULT/1000000.0D0
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE PRHONE(H,NFUNC,IHSM,NSM,IPACK)
*
*. Print one-electron integrals with symmetry IHSM
*. If IPACK = 1, then diagonal blocks are assumed packed
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION H(*), NFUNC(*)
#include "multd2h.inc"
*
      IOFF = 1
      DO IRSM = 1, NSM
        ICSM = MULTD2H(IHSM,IRSM)
        NR = NFUNC(IRSM)
        NC = NFUNC(ICSM)
        WRITE(6,'(A,2I3)')
     &  ' Blocks with row- and column-symmetry',IRSM,ICSM
        WRITE(6,'(A)')
     &  ' =========================================== '
        IF(IRSM.EQ.ICSM.AND.IPACK.EQ.1) THEN
          CALL PRSYM(H(IOFF),NR)
          IOFF = IOFF + NR*(NR+1)/2
        ELSE IF(IRSM.LT.ICSM .AND. IPACK.EQ.1) THEN
*. Upper block, neglected
        ELSE
          CALL WRTMT_LU(H(IOFF),NR,NC,NR,NC)
          IOFF = IOFF + NR*NC
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE RDVEC(LU,NSYM,NBAS,NORB,CMO,OCC,LOCC,TITLE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION NBAS(NSYM),NORB(NSYM),CMO(*),OCC(*)
      CHARACTER*(*) TITLE
      CHARACTER LINE*80,FMT*40
      LOGICAL SET
      FMT='(4E18.12)'
      REWIND (LU)
      KCMO  = 0
      NDIV  = 4
      TITLE = ' '
      SET   = .FALSE.
      DO 100 ISYM=1,NSYM
         DO 110 IORB=1,NORB(ISYM)
            DO 111 IBAS=1,NBAS(ISYM),NDIV
112            READ(LU,'(A80)',END=999,ERR=999) LINE
               IF(LINE(1:1).EQ.'*') THEN
                  IF(.NOT. SET) THEN
                     TITLE=LINE
                     SET=.TRUE.
                  END IF
                  GOTO 112
               END IF
               READ(LINE,FMT)
     &             (CMO(I+KCMO),I=IBAS,MIN(IBAS+3,NBAS(ISYM)))
111         CONTINUE
            KCMO=KCMO+NBAS(ISYM)
110      CONTINUE
100   CONTINUE
      IF(LOCC.EQ.0) RETURN
      KOCC=0
      DO 200 ISYM=1,NSYM
         DO 210 IORB=1,NORB(ISYM),NDIV
212         READ(LU,'(A80)',END=999,ERR=999) LINE
            IF(LINE(1:1).EQ.'*') THEN
               IF(.NOT. SET) THEN
                  TITLE=LINE
                  SET=.TRUE.
               END IF
               GOTO 212
            END IF
            READ(LINE,FMT) (OCC(I+KOCC),I=IORB,MIN(IORB+3,NORB(ISYM)))
210      CONTINUE
         KOCC=KOCC+NORB(ISYM)
200   CONTINUE
      RETURN
999   CONTINUE
      WRITE(*,*) '* ERROR IN RDVEC WHILE READING VECTOR SOURCE FILE *'
      Call Abend1( 20 )
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE REDBLK(NBLOCKI,IBLOCKI,IBBLOCKI,IOCOC,
     &                  NOCTPA,NOCTPB,
     &                  NBLOCKO,IBLOCKO,IMAP)
*
* a set of TTS blocks are given in IBLOCKI starting from
* block IBBLOCKI,
*
* Obtain those blocks that are allowed in space defined by IOCOC
*
* Jeppe Olsen , August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IBLOCKI(8,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
*. Output
      INTEGER IBLOCKO(8,*)
      INTEGER IMAP(*)
*
      IOFF = 1
      IOFFP = 1
      NBLOCKO = 0
      DO JBLOCKI = IBBLOCKI,IBBLOCKI+NBLOCKI-1
        IATP = IBLOCKI(1,JBLOCKI)
        IBTP = IBLOCKI(2,JBLOCKI)
        IF(IOCOC(IATP,IBTP).GT.0) THEN
          NBLOCKO = NBLOCKO + 1
          IMAP(NBLOCKO) = JBLOCKI
          CALL ICOPVE(IBLOCKI(1,JBLOCKI),IBLOCKO(1,NBLOCKO),8)
          IBLOCKO(5,NBLOCKO) = IOFF
          IBLOCKO(6,NBLOCKO) = IOFFP
          IOFF = IOFF + IBLOCKO(7,NBLOCKO)
          IOFFP= IOFFP+ IBLOCKO(8,NBLOCKO)
        END IF
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of blocks obtained in REDBLK ', NBLOCKO
        WRITE(6,*) ' output => input block map '
        CALL IWRTMA(IMAP,1,NBLOCKO,1,NBLOCKO)
        WRITE(6,*) ' Blocks obtained : '
        CALL IWRTMA(IBLOCKO,8,NBLOCKO,8,NBLOCKO)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE SETVCD(LUIN,LUOUT,SEGMNT,VALUE,IREW,LBLK)
*
* Construct a vector with the same structure as LUIN
* with values VALUE
*
* LBLK DEFINES STRUCTURE OF FILE
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
*
      IF( IREW .NE. 0 ) THEN
        CALL REWINE( LUIN ,LBLK)
        CALL REWINE( LUOUT ,LBLK)
      END IF

*
* Loop over blocks
*
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUOUT)
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC_LUCI(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          CALL SETVEC(SEGMNT,VALUE,LBL)
          CALL todsc_luci(SEGMNT,LBL,KBLK,LUOUT)
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE TYPE_TO_SYM_REO_MAT(XIN,XOUT)
*
*. a matrix XIN is given as NTOOB X NTOOB matrix in type form
*
* If ISYM.eq.1 matrix is assumed to be packed - lower half as usual
*.
*. Reorder to symmetry-ordered and -blocked matrix to give XOUT
*
*. Matrix is assumed to exclude inactive orbitals !!
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
#include "orbinp.inc"
*.
      DIMENSION XIN(NTOOB,NTOOB)
      DIMENSION XOUT(*)
*
      NTEST = 0
*
      DO ISMOB = 1, NSMOB
*. IOBOFF : Offset for active orbitals in symmetry order
        IF(ISMOB.EQ.1) THEN
          IOBOFF = NINOBS(1)+1
          IMTOFF = 1
        ELSE
          IOBOFF =
     &    IOBOFF + NTOOBS(ISMOB-1)-NINOBS(ISMOB-1)+NINOBS(ISMOB)
          IMTOFF = IMTOFF + NACOBS(ISMOB-1)**2
        END IF
        LOB = NACOBS(ISMOB)
*
*. Extract symmetry block of matrix
*
        DO IOB = IOBOFF,IOBOFF + LOB-1
           DO JOB = IOBOFF, IOBOFF + LOB-1
*. Corresponding type indeces
             IOBP = IREOST(IOB)
             JOBP = IREOST(JOB)
             XOUT(IMTOFF-1+(JOB-IOBOFF)*LOB+IOB-IOBOFF+1)
     &     = XIN(IOBP,JOBP)
           END DO
        END DO
*
      END DO
*. (End of loop over orbital symmetries )
      IF(NTEST.GE.10 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Symmetry ordered matrix '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        CALL APRBLM2(XOUT,NACOBS,NACOBS,NSMOB,0)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSMDP(VEC1,VEC2,FAC1,FAC2, LU1,LU2,LU3,IREW,LBLK)
C
C DISC VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
*
* Packed version, May 1996
C
C LBLK DEFINES STRUCTURE OF FILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
C
      IF(IREW .NE. 0 ) THEN
        CALL REWINE( LU1,LBLK)
        CALL REWINE( LU2,LBLK)
        CALL REWINE( LU3,LBLK)
      END IF
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE IF(LBLK .EQ. 0 ) THEN
          READ(LU1) NBL1
          READ(LU2) NBL2
          WRITE(LU3) NBL1
        ELSE IF (LBLK .LT. 0 ) THEN
          CALL IFRMDS( NBL1,1,-1,LU1)
          CALL IFRMDS( NBL2,1,-1,LU2)
          CALL ITODS ( NBL1,1,-1,LU3)
        END IF
        IF( NBL1 .NE. NBL2 ) THEN
        WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &  NBL1,NBL2
        Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
      END IF
C
      IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
        NO_ZEROING = 1
        CALL FRMDSC2(VEC1,NBL1,KBLK,LU1,IMZERO1,IAMPACK,NO_ZEROING)
        CALL FRMDSC2(VEC2,NBL1,KBLK,LU2,IMZERO2,IAMPACK,NO_ZEROING)
        IF( NBL1 .GT. 0 ) THEN
          IF(IMZERO1.EQ.1.AND.IMZERO2.EQ.1) THEN
*. Simple zero record
            CALL ZERORC(NBL1,LU3,IAMPACK)
          ELSE
*. Nonvanishing record
            ZERO = 0.0D0
            IF(IMZERO1.EQ.1) THEN
              CALL VECSUM(VEC1,VEC1,VEC2,ZERO,FAC2,NBL1)
            ELSE IF(IMZERO2.EQ.1) THEN
              CALL VECSUM(VEC1,VEC1,VEC2,FAC1,ZERO,NBL1)
            ELSE
              CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
            END IF
            CALL TODSCP(VEC1,NBL1,KBLK,LU3)
          END IF
        ELSE IF (NBL1.EQ.0) THEN
          CALL TODSCP(VEC1,NBL1,KBLK,LU3)
        END IF
      END IF
C
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
C                WRSVCD(LU1,LBLK,VEC1,ISCR1(IBASE),SCR1,NDEG,NDIM,
C    &           LUDIA)
      SUBROUTINE WRSVCD(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,JPACK)
*
* Write scattered vector to disc.
*.Vector is always written in packed form
*
      IMPLICIT REAL*8(A-H,O-Z)
      
*.Input
      DIMENSION IPLAC(*),VAL(*)
*.Scratch
      DIMENSION VEC1(*)
*
      IF(LBLK.GT.0) THEN
*. Write the vector without markers in one block
         CALL SETVEC(VEC1,0.0D0,NDIM)
         DO 100 IEFF = 1, NSCAT
           VEC1(IPLAC(IEFF)) = VAL(IEFF)
  100    CONTINUE
         CALL TODSC_LUCI(VEC1,NDIM,-1,LU)
      ELSE
*. Write the vector with the block format of file LUFORM
        CALL REWINE(LUFORM,-1)
        IBOT = 1
*. Loop over records
 1000   CONTINUE
*.Length
        CALL IFRMDS(LBL,1,-1,LUFORM)
        CALL ITODS (LBL,1,-1,LU)
        IF(LBL.GE.0) THEN
C?        write(6,*) ' IBOT, LBL  ',IBOT, LBL
          CALL SETVEC(VEC1,0.0D0,LBL)
          DO 200 IEFF = 1, NSCAT
            IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
     &      VEC1(IPLAC(IEFF)-IBOT+1) = VAL(IEFF)
*.CSK
  200     CONTINUE
          IF(JPACK.EQ.1) THEN
            CALL TODSCP(VEC1,LBL,-1,LU)
          ELSE
            CALL TODSC_LUCI(VEC1,LBL,-1,LU)
          END IF
*.Skip record on LUFORM
          CALL FRMDSC_LUCI(VEC1,LBL,-1,LUFORM,IMZERO,IAMPACK)
          IBOT = IBOT + LBL
          GOTO 1000
        END IF
*
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE WRTBLKN(VEC,NBLOCK,LBLOCK)
*
* Write the NBLOCK blocks of VEC
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC(*)
      DIMENSION LBLOCK(NBLOCK)
*
      IOFF = 1
      DO IBLOCK = 1, NBLOCK
        IF( LBLOCK(IBLOCK).GT.0) THEN
          WRITE(6,*) ' Block : ', IBLOCK
          WRITE(6,*) ' ==================='
          WRITE(6,*)
          WRITE(6,*) ' Length : ', LBLOCK(IBLOCK)
          WRITE(6,*)
          CALL WRITVE(VEC(IOFF),LBLOCK(IBLOCK))
          IOFF = IOFF + LBLOCK(IBLOCK)
          WRITE(6,*)
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE WRTRS2(VECTOR,ISMOST,ICBLTP,IOCOC,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,NSMST)
*
* Write RAS vector . Storage form is defined by ICBLTP
*
      IMPLICIT REAL*8           (A-H,O-Z)
*
      LOGICAL DIAGBL
*
      DIMENSION VECTOR(*)
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION NSASO(NSMST,* ),NSBSO(NSMST,* )
      DIMENSION ICBLTP(*),ISMOST(*)
*
*
      IBASE = 1
      DO 1000 IASM = 1, NSMST
        IBSM = ISMOST(IASM)
        IF(IBSM.EQ.0.OR.ICBLTP(IASM).EQ.0) GOTO 1000
*
        DO 900 IATP = 1, NOCTPA
          IF(ICBLTP(IASM).EQ.2) THEN
            IBTPMX = IATP
          ELSE
            IBTPMX = NOCTPB
          END IF
          NAST = NSASO(IASM,IATP)
*
          DO 800 IBTP = 1 , IBTPMX
            IF(IOCOC(IATP,IBTP) .EQ. 0 ) GOTO 800
            NBST = NSBSO(IBSM,IBTP)
            IF(ICBLTP(IASM).EQ.2.AND.IATP.EQ.IBTP ) THEN
* Diagonal block
              NELMNT = NAST*(NAST+1)/2
              IF(NELMNT.NE.0) THEN
                WRITE(6,'(A,3I3)')
     &          '  Iasm iatp ibtp : ', IASM,IATP,IBTP
                WRITE(6,'(A)')
     &          '  ============================'
                CALL PRSM2(VECTOR(IBASE),NAST)
                IBASE = IBASE + NELMNT
              END IF
            ELSE
              NELMNT = NAST*NBST
              IF(NELMNT.NE.0) THEN
                WRITE(6,'(A,3I3)')
     &          '  Iasm iatp ibtp : ', IASM,IATP,IBTP
                WRITE(6,'(A)')
     &          '  ============================'
                CALL WRTMT_LU(VECTOR(IBASE),NAST,NBST,NAST,NBST)
                IBASE = IBASE + NELMNT
              END IF
            END IF
  800     CONTINUE
  900   CONTINUE
 1000 CONTINUE
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE WRTTTS2(BLOCKS,IBLOCK,NBLOCK,IOFF,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,ISC)
*
* Print a batch of TTS blocks as given by IBLOCK.
* The blocks starts from block IOFF
*
*
* ISC = 1 : In slater determinant form
* ISC = 2 : In Combination        form
*
*. Jeppe Olsen, August 1995
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*. General input
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKS(*)
      INTEGER IBLOCK(8,*)
*
*
      WRITE(6,*) '  Batch of blocks '
      WRITE(6,*) ' ================= '
      WRITE(6,*)
      WRITE(6,'(A,I4)') ' Number of blocks in batch ', NBLOCK
*
      DO JBLOCK = IOFF,IOFF - 1 + NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
*
        IF (ISC.EQ.1 ) THEN
          JOFF = IBLOCK(5,JBLOCK)
        ELSE
          JOFF = IBLOCK(6,JBLOCK)
        END IF
*
*. Is this block diagonal
        IF(ISC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
C?      write(6,*) ' iatp ibtp iasm ibsm nia nib ',
C?   &  iatp,ibtp,iasm,ibsm,nia,nib
*
        IF(IPACK.EQ.1) THEN
          NELMNT = NIA*(NIA+1)/2
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            CALL PRSM2(BLOCKS(JOFF) ,NIA)
          END IF
        ELSE
          NELMNT = NIA*NIB
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            CALL WRTMT_LU(BLOCKS(JOFF) ,NIA,NIB,NIA,NIB)
          END IF
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE WRTTTSC(BLOCKS,IBLOCK,NBLOCK,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  NSASO,NSBSO,ISC,I0CHK,I0BLK)
*
* Print a batch of TTS blocks as given by IBLOCK
*
*
* ISC = 1 : In slater determinant form
* ISC = 2 : In Combination        form
*
*. Jeppe Olsen, from WRTTTS, July 97
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*. General input
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
*.
      DIMENSION BLOCKS(*)
      INTEGER IBLOCK(8,NBLOCK)
      INTEGER I0BLK(NBLOCK)
*
*
      WRITE(6,*) '  Batch of blocks '
      WRITE(6,*) ' ================= '
      WRITE(6,*)
      WRITE(6,'(A,I4)') ' Number of blocks in batch ', NBLOCK
*
      DO JBLOCK = 1, NBLOCK
*
        IATP = IBLOCK(1, JBLOCK)
        IBTP = IBLOCK(2, JBLOCK)
        IASM = IBLOCK(3, JBLOCK)
        IBSM = IBLOCK(4, JBLOCK)
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
*
        IF (ISC.EQ.1 ) THEN
          IOFF = IBLOCK(5,JBLOCK)
        ELSE
          IOFF = IBLOCK(6,JBLOCK)
        END IF
*
*. Is this block diagonal
        IF(ISC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
        NIA = NSASO(IASM,IATP)
        NIB = NSBSO(IBSM,IBTP)
C?      write(6,*) ' iatp ibtp iasm ibsm nia nib ',
C?   &  iatp,ibtp,iasm,ibsm,nia,nib
*
        IF(I0CHK.NE.0) THEN
          IF(I0BLK(JBLOCK).NE.0) THEN
            I_AM_ZERO = 1
          ELSE
            I_AM_ZERO = 0
          END IF
        ELSE
          I_AM_ZERO = 0
        END IF
        IF(IPACK.EQ.1) THEN
          NELMNT = NIA*(NIA+1)/2
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            IF(I_AM_ZERO.EQ.0) THEN
              CALL PRSM2(BLOCKS(IOFF) ,NIA)
            ELSE
              WRITE(6,*) ' Vanishing block'
            END IF
          END IF
        ELSE
          NELMNT = NIA*NIB
          IF(NELMNT.NE.0) THEN
            WRITE(6,'(A,3I3)')
     &      '  Iasm iatp ibtp : ', IASM,IATP,IBTP
            WRITE(6,'(A)')
     &      '  ============================'
            IF(I_AM_ZERO.EQ.0) THEN
              CALL WRTMT_LU(BLOCKS(IOFF) ,NIA,NIB,NIA,NIB)
            ELSE
              WRITE(6,*) ' Vanishing block'
            END IF
          END IF
        END IF
*
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE XDIAXT(XDX,X,DIA,NDIM,SCR)
*
* Obtain XDX = X * DIA * X(Transposed)
* where DIA is an diagonal matrix stored as a vector
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XDX(NDIM,NDIM)
      DIMENSION X(NDIM,NDIM),DIA(NDIM)
      DIMENSION SCR(NDIM,NDIM)
*
* DIA * X(transposed)
      DO 100 I=1,NDIM
        CALL COPVEC(X(1,I),SCR(1,I),NDIM)
        CALL SCALVE(SCR(1,I),DIA(I),NDIM)
  100 CONTINUE
* X * DIA * X(Transposed)
      CALL MATML4(XDX,X,SCR,NDIM,NDIM,NDIM,NDIM,NDIM,NDIM,2)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE XDIXT2(XDX,X,DIA,NXRDM,NXCDM,SHIFT,SCR)
*
* Obtain XDX = X * (DIA+Shift)-1 * X(Transposed)
* where DIA is an diagonal matrix stored as a vector
*
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION XDX(NXRDM,NXRDM)
      DIMENSION X(NXRDM,NXCDM),DIA(NXCDM)
      DIMENSION SCR(NXCDM)
*
      CALL SETVEC(XDX,0.0D0,NXRDM ** 2 )
      THRES = 1.0D-9
      DO 100 J=1,NXRDM
*.Scr(k) = X(J,K)/(SHIFT+DIA(K)
        DO 10 K = 1, NXCDM
          IF(ABS(SHIFT+DIA(K)) .GT. THRES ) THEN
            SCR(K) = X(J,K)/(SHIFT+DIA(K))
          ELSE
            SCR(K) = X(J,K)/THRES
          END IF
   10   CONTINUE
*
        DO 20 K = 1, NXCDM
          FACTOR = SCR(K)
          CALL VECSUM(XDX(1,J),XDX(1,J),X(1,K),1.0D0,FACTOR,NXRDM)
   20   CONTINUE
*
  100 CONTINUE
*
      NTEST = 00
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' x (Dia + shift ) - 1 X(T) from XDIXT2 '
        CALL WRTMT_LU(XDX,NXRDM,NXRDM,NXRDM,NXRDM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ZAP_BLOCK_VEC(LUIN,LBLK,IBLKS_A,SEGMNT,LUSCR)
*
* Zap blocks in vector in file LUIN for which IBLKS_A is zero
*
* Vector is initially constructed on LUSCR, and is copied back to
* LUIN after use
*
* Note : Files are always rewinded
* Packed version
*
* LBLK defines file type
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
      INTEGER IBLKS_A(*)
*
      CALL REWINE(LUIN  ,LBLK)
      CALL REWINE(LUSCR ,LBLK)
*
      IBLK = 0
*. Loop over blocks
C?      write(6,*) ' ZAP_BLOCK_VEC :  LBLK = ', LBLK
 1000 CONTINUE
        IBLK = IBLK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
          WRITE(LUSCR) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
          CALL ITODS (LBL,1,-1,LUSCR)
        END IF
*
        IF( LBL .GE. 0 ) THEN
*
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
*
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &    NO_ZEROING)
*
* Not packed
*
          IF(IAMPACK.EQ.0) THEN
            IF(IBLKS_A(IBLK).EQ.0) THEN
              ZERO = 0.0D0
              CALL SETVEC(SEGMNT,ZERO,LBL)
            END IF
            CALL TODSC_LUCI(SEGMNT,LBL,KBLK,LUSCR)
          END IF
*
* Packed
*
          IF(IAMPACK.EQ.1) THEN
            IF(IMZERO.EQ.1.OR.IBLKS_A(IBLK).EQ.0) THEN
              CALL ZERORC(LBL,LUSCR,IAMPACK)
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUSCR)
            END IF
          END IF
*
        END IF
*       ^ End if LBL .GE. 0


      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
*
*. And then copy back to LUIN
C          COPVCDP(LUIN,LUOUT,SEGMNT,IREW,LBLK)
      IREW = 1
      CALL COPVCD(LUSCR,LUIN,SEGMNT,IREW,LBLK)
*
      RETURN
      END
