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
      SUBROUTINE DIA0TRM_GAS(ITASK,LUIN,LUOUT,VEC,VEC2,
     &                       FACTOR)
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 2)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 1)
*
* Note : Opposite to normal DIATRM_GAS !!!
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLJ, KLK, KLSCR2, KLXA, KLXB, KLSCR, KLH1D
      INTEGER*8 KLASTR, KLBSTR, KLRJKA, KSVST
      INTEGER*8 KLVIOIO, KLVBLTP, KLVLBT, KLVLEBT, KLVI1BT, KLVIBT
      ! for addressing of WORK
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
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "oper.inc"
#include "crun.inc"
#include "glbbas.inc"
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC

*
      DIMENSION VEC(*)
*
      CALL QENTER('DIA0T')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIA0TR')
*
      ISM = ISSM
      ISPC = ISSPC
      WRITE(6,*) ' DIA0TRM : ISSM ISSPC :', ISSM,ISSPC
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRDIA)
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
        WRITE(6,*) ' ================'
        WRITE(6,*) ' DIA0TRM speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        write(6,*) ' ISPC = ', ISPC
        WRITE(6,*) ' LUIN,LUOUT ', LUIN,LUOUT
        WRITE(6,*) ' IPART = ', IPART
        WRITE(6,*) ' FACTOR = ', FACTOR
      END IF
*
*
*. Perturbation operator
*
      IF(IPART.EQ.1) THEN
*. Moller-Plesset partitioning
        I12 = 1
        IPERTOP = 1
      ELSE IF(IPART.EQ.2) THEN
*. Epstein-Nesbet Partitioning
       I12 = 2
       IPERTOP = 0
      END IF
*
      IF(NTEST.GE.10)WRITE(6,*) ' I12, IPERTOP',I12,IPERTOP

*. A bit of scratch
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')
      CALL MEMMAN(KLXA  ,NACOB,   'ADDL  ',2,'KLXA  ')
      CALL MEMMAN(KLXB  ,NACOB,   'ADDL  ',2,'KLXB  ')
      CALL MEMMAN(KLSCR ,2*NACOB, 'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KLH1D ,NACOB,   'ADDL  ',2,'KLH1D ')
*. Space for blocks of strings
C     WRITE(6,*) ' MXNSTR in DIATERM', MXNSTR
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NAEL,'ADDL  ',1,'KLBSTR')
      MAXA = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')
*. Diagonal of one-body integrals and coulomb and exchange integrals
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFIO),WORK(KINT1O),NINT1)
      CALL GT1DIA(WORK(KLH1D))
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFIO),WORK(KINT1O),NINT1)
      WRITE(6,*) ' DIA0TRM_GAS : IPERTOP ', IPERTOP
*
      IF(I12.EQ.2)
     &CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
*. Interchange
      IF(ITASK.EQ.1) THEN
        JTASK = 2
      ELSE
        JTASK = 1
      END IF
*
*
*. Iblock driven, so just set up ....
*
      NTTS = MXNTTS
      NOOS = NOCTPA*NOCTPB*NSMCI
      CALL MEMMAN(KLVIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'VIOIO ')
      CALL MEMMAN(KLVBLTP,NSMST,'ADDL  ',2,'VBLTP ')
*
      CALL IAIBCM(ISPC,WORK(KLVIOIO))
      KSVST = 1
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLVBLTP),WORK(KSVST))
*
*. Space for partitioning of vectors
      NTTS = MXNTTS
      CALL MEMMAN(KLVLBT ,NTTS  ,'ADDL  ',1,'VLBT  ')
      CALL MEMMAN(KLVLEBT ,NTTS ,'ADDL  ',1,'VLEBT ')
      CALL MEMMAN(KLVI1BT,NTTS  ,'ADDL  ',1,'VI1BT ')
      CALL MEMMAN(KLVIBT ,8*NTTS,'ADDL  ',1,'VIBT  ')
*
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KLVBLTP),WORK(KNSTSO(IATP)),
     &               WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,
     &               NSMST,LBLOCK,WORK(KLVIOIO),ISMOST(1,ISM),
     &               NSBATCH,WORK(KLVLBT),WORK(KLVLEBT),
     &               WORK(KLVI1BT),WORK(KLVIBT),1,ITTSS_ORD)
C    PART_CIV2(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
C    &                  NSMST,MXLNG,IOCOC,ISMOST,
C    &                  NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP,
C    &                  ITTSS_ORD)
*
      NBLOCKT = IFRMR(WORK(KLVLBT),1,1)
C     WRITE(6,*) ' NBLOCKT = ', NBLOCKT

      ECORES = 0.0D0
      CALL DIATERMS_GAS(NAEL,WORK(KLASTR),NBEL,WORK(KLBSTR),
     &             NACOB,VEC,NSMST,WORK(KLH1D),
     &             IDC,WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &             WORK(KLK),WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &             ECORES,LUIN,LUOUT,
     &             IPRDIA,NTOOB,
     &             WORK(KLRJKA),
     &             I12,WORK(KLVIBT),NBLOCKT,JTASK,FACTOR,0,0)

      IF(NTEST.GE.100.AND.LUOUT.EQ.0) THEN
        WRITE(6,*)  ' output vector from DIA0TRM '
        CALL WRTTTS(VEC,WORK(KLVIBT),NBLOCKT,
     &              NSMST,NOCTPA,NOCTPB,
     &              WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),IDC)
      END IF
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'DIA0TR')
      CALL QEXIT('DIA0T')
*
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIABLK(IASM,IATP,IBSM,IBTP,IFULL,DIAG)
*
* Obtain diagonal block of determinant block
* IASM,IATP,IBSM,IBTP
*
*
*. Type of operator is taken from OPER
*. Type of CI space is taken from CANDS
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLJ, KLK, KLSCR2, KLXA, KLXB, KLSCR, KLH1D
      INTEGER*8 KLSMOS, KLBLTP, KLIOIO, KLSVST, KLRJKA
      ! for addressing of WORK
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cprnt.inc"
#include "crun.inc"
#include "oper.inc"
*
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
C     COMMON/OPER/I12,IPERTOP
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
*
      CALL QENTER('DIABL')
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRDIA)
*
** Specifications of internal space
*
      ISM = ICSM
      ISPC = ICSPC
*
      IAGP = IASTFI(ISPC)
      IBGP = IBSTFI(ISPC)
      NAEL = NELEC(IAGP)
      NBEL = NELEC(IBGP)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' DIABLK speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IAGP IBGP NAEL NBEL '
        WRITE(6,*)   IAGP,IBGP,NAEL,NBEL
        WRITE(6,*) ' I12,IPERTOP',I12,IPERTOP
      END IF

      NOCTPA = NOCTYP(IAGP)
      NOCTPB = NOCTYP(IBGP)
      MNRS1C = MNR1IC(ISPC)
      MXRS3C = MXR3IC(ISPC)
      IF(NTEST.GE.10) THEN
        write(6,*) ' NOCTPA NOCTPB MNRS1C,MXRS3C'
        write(6,*)   NOCTPA,NOCTPB,MNRS1C,MXRS3C
      END IF
*
**. Local memory
*
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',IDUM,'DIABLK')
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')
      CALL MEMMAN(KLXA  ,NACOB,   'ADDL  ',2,'KLXA  ')
      CALL MEMMAN(KLXB  ,NACOB,   'ADDL  ',2,'KLXB  ')
      CALL MEMMAN(KLSCR ,2*NACOB, 'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KLH1D ,NACOB,   'ADDL  ',2,'KLH1D ')
      CALL MEMMAN(KLSMOS,NSMST,   'ADDL  ',2,'KLSMOS')
      CALL MEMMAN(KLBLTP,NSMST,   'ADDL  ',2,'KLSMOS')
      CALL MEMMAN(KLIOIO,NOCTPA*NOCTPB,   'ADDL  ',2,'KLIOIO')
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KLSVST,NSMST,   'ADDL  ',2,'KLSVST')
      ELSE
        KLSVST = 1
      END IF
      MAXA = IMNMX(WORK(KNSTSO(IAGP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')
*
** Info on block structure of internal state
*
      IF(IDC.EQ.3.OR.IDC.EQ.4)
     &CALL SIGVST(WORK(KLSVST),NSMST)
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLTP),WORK(KLSVST))
      CALL QUIT( 'In DIABLK: IAIBCM_GAS call argument list wrong!')
      CALL IAIBCM_GAS(MNRS1C,MXRS3C,
     &                NOCTPA,NOCTPB,
     &     WORK(KEL1(IAGP)),WORK(KEL3(IAGP)),
     &     WORK(KEL1(IBGP)),WORK(KEL3(IBGP)),WORK(KLIOIO),IPRDIA)
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
      CALL GT1DIA(WORK(KLH1D))
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
      ECORES = 0.0D0
      IF(I12.EQ.2)
     &CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
*
      CALL DIABLKS(NAEL,WORK(KOCSTR(IAGP)),NBEL,WORK(KOCSTR(IBGP)),
     &             NACOB,NSMST,WORK(KLH1D),
     &             ISMOST(1,ISM),WORK(KLBLTP),
     &             WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &             WORK(KLK),WORK(KNSTSO(IAGP)),WORK(KNSTSO(IBGP)),
     &             WORK(KLIOIO),NOCTPA,NOCTPB,WORK(KISTSO(IAGP)),
     &             WORK(KISTSO(IBGP)),ECORES,
     &             PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &             WORK(KLRJKA),I12,
     &             IASM,IATP,IBSM,IBTP,IFULL,DIAG)
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'DIABLK')
      CALL QEXIT('DIABL')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIABLKS(NAEL,IASTR,NBEL,IBSTR,
     &                  NORB,NSMST,H,
     &                  ISMOST,IBLTP,XA,XB,SCR,RJ,RK,
     &                  NSSOA,NSSOB,IOCOC,NOCTPA,NOCTPB,
     &                  ISSOA,ISSOB,ECORE,
     &                  PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,I12,
     &                  IASM,IATP,IBSM,IBTP,IFULL,DIAG)
*
*
* Obtain specific block of diagonal
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, February 1994 , obtained from CIDIA4
*
* I12 = 1 => only one-body part
*     = 2 =>      one+two-body part
*
      IMPLICIT REAL*8           (A-H,O-Z)
      REAL * 8  INPROD
*.General input
      DIMENSION NSSOA(NOCTPA,*),NSSOB(NOCTPB,* )
      DIMENSION ISSOA(NOCTPA,*),ISSOB(NOCTPB,*)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION H(NORB)
*. Specific input
      DIMENSION IOCOC(NOCTPA,NOCTPB)
      DIMENSION ISMOST(*),IBLTP(*)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION DIAG(*)
*
      NTESTL = 000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' DIABLKS in action '
        WRITE(6,*) ' ================= '
        WRITE(6,*)
        WRITE(6,*) ' IATP IBTP IASM IBSM '
        WRITE(6,*)   IATP,IBTP,IASM,IBSM
      END IF
*
      IF( NTEST .GE. 2000 ) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMT_LU(H,1,NORB,1,NORB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMT_LU(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMT_LU(RK,NORB,NORB,NTOOB,NTOOB)
        END IF
      END IF
*
**3 Diagonal elements according to Handys formulae
*   (corrected for error)
*
*   DIAG(IDET) = HII*(NIA+NIB)
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIA*NJA
*              + 0.5 * ( J(I,J)-K(I,J) ) * NIB*NJB
*              +         J(I,J) * NIA*NJB
*
*. K goes to J - K
      ONE = 1.0D0
      ONEM = -1.0D0
      IF(I12.EQ.2)
     &CALL VECSUM(RK,RK,RJ,ONEM,ONE,NTOOB **2)
*
        IF(IFULL.EQ.0.AND.IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*
          IF(IREST1.EQ.1) THEN
            MXBTP = IATP
          ELSE
            MXBTP = NOCTPB
          END IF
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
          IOFF =  ISSOA(IATP,IASM)
          DO IA = IOFF,IOFF+NSSOA(IATP,IASM)-1
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
*
          IBSTRT = ISSOB(IBTP,IBSM)
          IBSTOP = IBSTRT + NSSOB(IBTP,IBSM)-1
*
          IDET = 0
          DO 899 IB = IBSTRT,IBSTOP
            IBREL = IB - IBSTRT + 1
*
*. Array for terms depending only upon IB
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
                DO  JEL = 1, NBEL
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
              IASTRT = ISSOA(IATP,IASM) - 1 + IBREL
            ELSE
              IASTRT = ISSOA(IATP,IASM)
            END IF
            IASTOP = ISSOA(IATP,IASM) + NSSOA(IATP,IASM) - 1
            DO 800 IA = IASTRT,IASTOP
              IDET = IDET + 1
              X = EB + RJKAA(IA-IOFF+1)
              DO 890 IEL = 1, NAEL
                X = X +XB(IASTR(IEL,IA))
  890         CONTINUE
              DIAG(IDET) = X
  800       CONTINUE
  899     CONTINUE
*
          IF(Ntest.ge.1000) THEN
             WRITE(6,*) ' Next batch of diagonal elements '
             CALL WRTMT_LU(DIAG,1,IDET,1,IDET)
          END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIATERM2_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,
     &                       JPERT,J12,JDC)
* = DIATERM_GAS, just J12 added !
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
*
* For the NBLOCKS givem in IBLOCK starting from BLOCK IOFF
*
* If JPERT.NE.0, the perturbation operator as defined by IPART is used.
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLJ, KLK, KLSCR2, KLXA, KLXB, KLSCR, KLH1D, KLASTR,
     &          KLBSTR, KLRJKA
!               for addressing of WORK
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
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "glbbas.inc"
#include "oper.inc"
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      INTEGER IBLOCK(8,*)
*
      DIMENSION VEC(*)
*
      CALL QENTER('DIATR')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIATRM')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
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
*
C     IF(JPERT.EQ.0) THEN
*. Use full Hamiltonian
C       I12 = 2
C       IPERTOP = 0
C     ELSE
*. Use perturbation operator
C       IF(IPART.EQ.1) THEN
*. Moller-Plesset partitioning
C         I12 = 1
C         IPERTOP = 1
C       ELSE IF(IPART.EQ.2) THEN
*. Epstein-Nesbet Partitioning
C         I12 = 2
C         IPERTOP = 0
C       END IF
C     END IF

      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ========================='
        WRITE(6,*) '   DIATERM2_GAS speaking '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        WRITE(6,*) ' JPERT,IPART,J12,IPERTOP',JPERT,IPART,J12,IPERTOP
      END IF
*. A bit of scracth
  
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')
      CALL MEMMAN(KLXA  ,NACOB,   'ADDL  ',2,'KLXA  ')
      CALL MEMMAN(KLXB  ,NACOB,   'ADDL  ',2,'KLXB  ')
      CALL MEMMAN(KLSCR ,2*NACOB, 'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KLH1D ,NACOB,   'ADDL  ',2,'KLH1D ')
*. Space for blocks of strings
C     WRITE(6,*) ' MXNSTR in DIATERM', MXNSTR
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NAEL,'ADDL  ',1,'KLBSTR')
      MAXA = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')
*. Diagonal of one-body integrals and coulomb and exchange integrals
*. Integrals assumed in place so :
C!    IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
      CALL GT1DIA(WORK(KLH1D))
C!    IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
      IF(J12.EQ.2)
     &CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
*. Core energy not included
      ECOREP = 0.0D0
      CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
*
      CALL DIATERMS_GAS(NAEL,WORK(KLASTR),NBEL,WORK(KLBSTR),
     &             NACOB,VEC,NSMST,WORK(KLH1D),
     &             JDC,WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &             WORK(KLK),WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &             ECOREP,0,0,
     &             IPRDIA,NTOOB,
     &             WORK(KLRJKA),
     &             J12,IBLOCK(1,IOFF),NBLOCK,ITASK,FACTOR,0,0)
C    &                  IBLOCK,NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'DIATRM')
      CALL QEXIT('DIATR')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' output vector from DIATRM '
        CALL WRTTTS(VEC,IBLOCK(1,IOFF),NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,
     &              WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),IDC)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIATERM_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,
     &                       JPERT,I0CHK,I0BLK)
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
*
* For the NBLOCKS givem in IBLOCK starting from BLOCK IOFF
*
* If JPERT.NE.0, the perturbation operator as defined by IPART is used.
*
* IF ICHBLKS = 1, entries in IZBLKS are checked for zero blocks
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLJ, KLK, KLSCR2, KLXA, KLXB, KLSCR, KLH1D, KLASTR,
     &          KLBSTR, KLRJKA
!               for addressing of WORK
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
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "glbbas.inc"
#include "oper.inc"
*
C     COMMON/CECORE/ECORE,ECORE_ORIG
#include "cecore.inc"
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      INTEGER IBLOCK(8,*)
*
      DIMENSION VEC(*)
*
      CALL QENTER('DIATR')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIATRM')
*
      NTEST = 00
      NTEST = MAX(NTEST,IPRDIA)
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
*
      IF(JPERT.EQ.0) THEN
*. Use full Hamiltonian
        I12 = 2
        IPERTOP = 0
      ELSE
*. Use perturbation operator
        IF(IPART.EQ.1) THEN
*. Moller-Plesset partitioning
          I12 = 1
          IPERTOP = 1
        ELSE IF(IPART.EQ.2) THEN
*. Epstein-Nesbet Partitioning
          I12 = 2
          IPERTOP = 0
        END IF
      END IF

      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' DIATERM speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        WRITE(6,*) ' JPERT,IPART,I12,IPERTOP',JPERT,IPART,I12,IPERTOP
      END IF
*. A bit of scratch
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')
      CALL MEMMAN(KLXA  ,NACOB,   'ADDL  ',2,'KLXA  ')
      CALL MEMMAN(KLXB  ,NACOB,   'ADDL  ',2,'KLXB  ')
      CALL MEMMAN(KLSCR ,2*NACOB, 'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KLH1D ,NACOB,   'ADDL  ',2,'KLH1D ')
*. Space for blocks of strings
C     WRITE(6,*) ' MXNSTR in DIATERM', MXNSTR
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NAEL,'ADDL  ',1,'KLBSTR')
      MAXA = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')
*. Diagonal of one-body integrals and coulomb and exchange integrals
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1O),NINT1)
      CALL GT1DIA(WORK(KLH1D))
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Diagonal of 1-el ints '
        CALL WRTMT_LU(WORK(KLH1D),1,NTOOB,1,NTOOB)
      END IF
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1O),NINT1)
      IF(I12.EQ.2)
     &CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
*. Core energy not included
      ECOREP = 0.0D0
*
      SHIFT = ECORE_ORIG-ECORE
      FACTORX = FACTOR + SHIFT
      IF(NTEST.GE.10 ) THEN
        WRITE(6,*) ' SHIFT and FACTORX', SHIFT,FACTORX
        WRITE(6,*) ' ECORE_ORIG, ECORE', ECORE_ORIG,ECORE
      END IF
*
      CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
      CALL DIATERMS_GAS(NAEL,WORK(KLASTR),NBEL,WORK(KLBSTR),
     &             NACOB,VEC,NSMST,WORK(KLH1D),
     &             IDC,WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &             WORK(KLK),WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &             ECOREP,0,0,
     &             IPRDIA,NTOOB,
     &             WORK(KLRJKA),
     &             I12,IBLOCK(1,IOFF),NBLOCK,ITASK,FACTORX,I0CHK,I0BLK)
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'DIATRM')
      CALL QEXIT('DIATR')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' output vector from DIATRM '
        CALL WRTTTSC(VEC,IBLOCK(1,IOFF),NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,
     &              WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),IDC,
     &              I0CHK,I0BLK)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIATERMS_GAS(NAEL,IASTR,NBEL,IBSTR,
     &                  NORB,VEC,NSMST,H,
     &                  IDC,XA,XB,SCR,RJ,RK,
     &                  NSSOA,NSSOB,
     &                  ECORE,LUIN,LUOUT,
     &                  IPRNT,NTOOB,RJKAA,I12,
     &                  IBLOCK,NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)
*
* Terms from diagonal to specific blocks
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
*
* Calculate determinant diagonal
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
      IMPLICIT REAL*8           (A-H,O-Z)
*.General input
      DIMENSION NSSOA(NSMST,*), NSSOB(NSMST,*)
      DIMENSION H(NORB)
      DIMENSION IBLOCK(8,*)
*.
      INTEGER I0BLK(*)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
*. Output
      DIMENSION VEC (*)
*
      NTEST =  00
      NTEST = MAX(NTEST,IPRNT)
C?    WRITE(6,*) ' NTEST = ',NTEST
*
      IF(LUIN.GT.0) REWIND LUIN
      IF(LUOUT.GT.0) REWIND LUOUT

      IF( NTEST .GE. 20 ) THEN
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' DIATERMS_GAS in action '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' LUIN,LUOUT = ', LUIN,LUOUT
        WRITE(6,*) ' NBLOCK =', NBLOCK
        WRITE(6,*) ' I0CHK = ', I0CHK
      END IF
*
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Diagonal one electron integrals'
        CALL WRTMT_LU(H,1,NORB,1,NORB)
        IF(I12.EQ.2) THEN
          WRITE(6,*) ' Coulomb and exchange integrals '
          CALL WRTMT_LU(RJ,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*)
          CALL WRTMT_LU(RK,NORB,NORB,NTOOB,NTOOB)
          WRITE(6,*) ' I12 and ITASK = ', I12,ITASK
        END IF
      WRITE(6,*) ' FACTOR = ',FACTOR
      END IF
*
**3 Diagonal elements according to Handys formulae
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
*
      ITDET = 0
      IDET = 0
      DO JBLOCK = 1, NBLOCK
        IF(IBLOCK(1,JBLOCK).GT.0) THEN
        IATP = IBLOCK(1,JBLOCK)
        IBTP = IBLOCK(2,JBLOCK)
        IASM = IBLOCK(3,JBLOCK)
        IBSM = IBLOCK(4,JBLOCK)
        IOFF = IBLOCK(6,JBLOCK)
        IF(NTEST.GE.20) THEN
         WRITE(6,*) ' Block in action : IATP IBTP IASM IBSM ',
     &               IATP,IBTP,IASM,IBSM
        END IF
*
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP) THEN
          IPACK = 1
        ELSE
          IPACK = 0
        END IF
*
*
*. Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
*                           0.5*SUM(I,J) ( J(I,J) - K(I,J))*N(I)*N(J)
*
*. Obtain alpha strings of sym IASM and type IATP
        IDUM = 0
        CALL GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,
     &                           NORB,0,IDUM,IDUM)
        IF(NTEST.GE.1000) THEN
          write(6,*) ' After GETSTR for A strings '
          WRITE(6,*) ' alpha strings obtained '
          NAST = NSSOA(IASM,IATP)
          CALL IWRTMA(IASTR,NAEL,NAST,NAEL,NAST)
        END IF
*
        IOFF =  1
        NIA = NSSOA(IASM,IATP)
        DO IA = 1 ,NSSOA(IASM,IATP)
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
*. Obtain alpha strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,
     &                         NORB,0,IDUM,IDUM)
        NIB =  NSSOB(IBSM,IBTP)
*
        IMZERO=0
        IF(LUIN.GT.0) THEN
          CALL IFRMDS(LDET,1,-1,LUIN)
          IDET = 0
          CALL FRMDSC_LUCI(VEC(1),LDET,-1,LUIN,IMZERO,IAMPACK)
        END IF
*
        IF(I0CHK.EQ.1) THEN
          IMZERO = I0BLK(JBLOCK)
          IF(IMZERO.EQ.1) THEN
*.Update offset to next block
            IF(IPACK.EQ.1.AND.IATP.EQ.IBTP) THEN
              IDET = IDET + NIA*(NIA+1)/2
            ELSE
              IDET = IDET + NIA*NIB
            END IF
          END IF
        END IF
C?      WRITE(6,*) ' DIATERMS_GAS : I0CHK,JBLOCK IMZERO',
C?   &  I0CHK,JBLOCK,IMZERO
*
        IF(IMZERO.NE.1) THEN
*. Calculate ...
*
        DO IB = 1 ,NIB
*
*. Terms depending only on IB
*
          HB = 0.0D0
          RJBB = 0.0D0
          CALL SETVEC(XB,0.0D0,NORB)

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
          IF(IPACK.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1
          END IF
*
          IASTOP = NSSOA(IASM,IATP)
          DO IA = IASTRT,IASTOP
            IDET = IDET + 1
            ITDET = ITDET + 1
            X = EB + RJKAA(IA-IOFF+1)
            DO IEL = 1, NAEL
              X = X +XB(IASTR(IEL,IA))
            END DO
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
            IF(ITASK.EQ.1) THEN
              IF(ABS(X+FACTOR) .GT. 1.0D-10) THEN
                VEC(IDET) = VEC(IDET)/(X+FACTOR)
              ELSE
                VEC(IDET) = 0.0D0
              END IF
            ELSE
              VEC(IDET) = VEC(IDET)*(X+FACTOR)
            END IF
C?         write(6,*) ' IDET,X,VEC(IDET) ', IDET,X,VEC(IDET)
          END DO
        END DO
        END IF
*
        IF(LUOUT.GT.0) THEN
          CALL ITODS(LDET,1,-1,LUOUT)
          CALL TODSC_LUCI(VEC,LDET,-1,LUOUT)
C?        WRITE(6,*) ' Number of elements transferred to DISC ',
C?   &    LDET
          IDET = 0
        END IF
*
      END IF
      END DO
*
      IF(LUOUT.GT.0) THEN
       IONEM = -1
       CALL ITODS(IONEM,1,-1,LUOUT)
      END IF
*
C?    WRITE(6,*) ' Mission DIATERMS finished '
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE GASDIAS(NAEL,IASTR,NBEL,IBSTR,
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
* last parallel revision: March 2008
*
      IMPLICIT REAL*8           (A-H,O-Z)
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) IDIA_OFFSET
#endif
#include "parluci.h"
C     REAL * 8  INPROD
*.General input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION H(NORB)
*. Specific input
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*), NPARBLOCK(*)
*. Output
      DIMENSION DIAG(*)
*
      NTEST =  0000
      NTEST = MAX(NTEST,IPRNT)
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0
      ELSE
         XADD = 0.0D0
      END IF
      IDIA_OFFSET = 0
*     proper offset
      IDIA_OFFSET = IDIA_OFFSET + MY_DIA_OFF
*
      IF( NTEST .GE. 20 ) THEN
        WRITE(LUWRT,*) ' Diagonal one electron integrals'
        CALL WRTMATMN(H,1,NORB,1,NORB,LUWRT)
        WRITE(LUWRT,*) ' Core energy ', ECORE
        IF(I12.EQ.2) THEN
          WRITE(LUWRT,*) ' Coulomb and exchange integrals '
          CALL WRTMATMN(RJ,NORB,NORB,NTOOB,NTOOB,LUWRT)
          WRITE(LUWRT,*)
          CALL WRTMATMN(RK,NORB,NORB,NTOOB,NTOOB,LUWRT)
        END IF
*
        WRITE(LUWRT,*) ' TTSS for Blocks '
        DO IBLOCK = 1, NBLOCK
          WRITE(LUWRT,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
        END DO
*
        WRITE(luwrt,*) ' I12 = ',I12
        if(LUCI_NMPROC .gt. 1)then 
          WRITE(LUWRT,*) ' NPARBLOCK:'
          CALL IWRTMA(NPARBLOCK,1,NBLOCK,1,NBLOCK)
        end if
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
      DO IBLK = 1, NBLOCK

         if(LUCI_NMPROC .gt. 1)then
           IF( NPARBLOCK( IBLK ) .ne. LUCI_MYPROC ) GOTO 2425
         end if
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
          CALL DZERO(XB,NORB)
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
          IF(NTEST.GE.20) THEN
            if(IDET.gt.0)then
             write(LUWRT,*) ' number of diagonal elements to disc ',IDET
             CALL WRTMT_LU(DIAG,1,IDET,1,IDET)
            end if
          END IF
          IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
             CALL MPI_FILE_WRITE_AT(IDIA,IDIA_OFFSET,DIAG,IDET,
     &                           MPI_REAL8,my_STATUS,IERR)
C            new offset
             IDIA_OFFSET = IDIA_OFFSET + IDET
#else
             call quit('LUCI_NMPROC .gt. 1 but not VAR_MPI')
#endif
          ELSE
             CALL ITODS(IDET,1,-1,LUDIA)
             CALL TODSC_LUCI(DIAG,IDET,-1,LUDIA)
          END IF
          IDET = 0
        END IF
2425    CONTINUE
      END DO
*     ^ End of loop over blocks

      IF(NTEST.GE.5) WRITE(LUWRT,*)
     &' Number of diagonal elements generated ',ITDET
*
      IF(NTEST .GE.20 .AND.ICISTR.LE.1 )THEN
        WRITE(LUWRT,*) ' CIDIAGONAL '
        CALL WRTMATMN(DIAG(1),1,IDET,1,IDET,LUWRT)
      END IF
*
      IF ( ICISTR.GE.2 ) CALL ITODS(-1,1,-1,LUDIA)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE GASDIAT(DIAG,LUDIA,ECORE,ICISTR,I12,
     &                   IBLTP,NBLOCK,IBLKFO,NPARBLOCK)
*
* CI diagonal in SD basis for state with symmetry ISM in internal
* space ISPC
*
* GAS version, Winter of 95
*
* Driven by table of TTS blocks, May97
*
* parallel adaption by S. Knecht - March 2008
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLJ, KLK, KLSCR2, KLXA, KLXB, KLSCR, KLH1D, KLASTR,
     &          KLBSTR, KLRJKA
!               for addressing of WORK
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
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
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
C
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif
#include "parluci.h"

*
      DIMENSION IBLTP(*)
      DIMENSION IBLKFO(8,NBLOCK), NPARBLOCK(*)
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
*
      CALL QENTER('GASDIAT')
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
        WRITE(luwrt,*) ' ================'
        WRITE(luwrt,*) ' GASDIA speaking '
        WRITE(luwrt,*) ' ================'
        WRITE(luwrt,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(luwrt,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(luwrt,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
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
CE-JULY29-99      CALL MEMMAN(KLBSTR,MXNSTR*NAEL,'ADDL  ',1,'KLBSTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
*
C     IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
C       CALL MEMMAN(KLSVST,NSMST,   'ADDL  ',2,'KLSVST')
C     ELSE
C       KLSVST = 1
C     END IF
      MAXA = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      CALL GT1DIA(WORK(KLH1D))
      CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
      IF( LUDIA .GT. 0 ) CALL REWINO(LUDIA)
      CALL GASDIAS(NAEL,WORK(KLASTR),NBEL,WORK(KLBSTR),
     &     NACOB,DIAG,NSMST,WORK(KLH1D),
     &     WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &     WORK(KLK),WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &     LUDIA,ECORE,PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &     WORK(KLRJKA),I12,IBLTP,NBLOCK,IBLKFO,NPARBLOCK)
*.Flush local memory
      CALL MEMMAN(IDUM,  IDUM,    'FLUSM ',IDUM,'GASDIA')
      CALL QEXIT('GASDIAT')
!     call quit('*** stefan forced me to stop after diagonal. ***')
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_DIAG_BLOC_MAT(A,ADIAG,NBLOCK,LBLOCK,ISYM)
*
* Obtain diagonal elements from symmetry blocked matrix
*
*
* ISYM = 1 => Input and output are     triangular packed
*      else=> Input and Output are not triangular packed
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A,H,O-Z)
*. Input
      DIMENSION A(*)
      INTEGER LBLOCK(*)
*. Output
      DIMENSION ADIAG(*)
*
      DO IBLOCK = 1, NBLOCK
        IF(IBLOCK.EQ.1) THEN
          IOFF = 1
          LOFF = 1
        ELSE
          IF(ISYM.EQ.1) THEN
            IOFF = IOFF + LBLOCK(IBLOCK-1)*(LBLOCK(IBLOCK-1)+1)/2
          ELSE
            IOFF = IOFF + LBLOCK(IBLOCK-1)** 2
          END IF
          LOFF = LOFF + LBLOCK(IBLOCK-1)
        END IF
*
        L = LBLOCK(IBLOCK)
        CALL COPDIA(A(IOFF),ADIAG(LOFF),L,ISYM)
C            COPDIA(A,VEC,NDIM,IPACK)
      END DO
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        NDIM = IELSUM(LBLOCK,NBLOCK)
        WRITE(luwrt,*) ' output matrix GET_DIAG_BLOC_MAT '
        WRITE(luwrt,*) ' ==============================='
        CALL WRTMT_LU(ADIAG,1,NDIM,1,NDIM)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GT1DIA(H1DIA)
*
* Obtain diagonal of one electron matrix over active
* orbitals
*
*. Dec 97 : obtained from KINT1O
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"

*.GLobal pointers
C     COMMON/GLBBAS/KINT1,KINT2,KPINT1,KPINT2,KLSM1,KLSM2,KRHO1
#include "glbbas.inc"

#include "lucinp.inc"
#include "orbinp.inc"
*
CINA  CALL GT1DIS(H1DIA,IREOTS(1+NINOB),WORK(KPINT1),WORK(KINT1),
CINA &            ISMFTO,IBSO,NACOB)
      CALL GT1DIS(H1DIA,IREOTS(1),WORK(KPINT1),WORK(KINT1O),
     &            ISMFTO,IBSO,NACOB)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GT1DIS(H1DIA,IREOTS,IPNT,H,ISMFTO,IBSO,NACOB)
*
* diagonal of one electron integrals over active orbitals
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
*.Input
      INTEGER IREOTS(*),IPNT(*),ISMFTO(*),IBSO(*)
      DIMENSION H(*)
*.Output
      DIMENSION H1DIA(*)
*
      DO 100 IIOB = 1, NACOB
        IOB = IREOTS(IIOB)
        ISM = ISMFTO(IIOB)
        IOBREL = IOB-IBSO(ISM)+1
C?      WRITE(6,*) ' IIOB IOB ISM IOBREL '
C?      WRITE(6,*)   IIOB,IOB,ISM,IOBREL
        H1DIA(IIOB) = H(IPNT(ISM)-1+IOBREL*(IOBREL+1)/2)
  100 CONTINUE
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(lupri,*) ' Diagonal one electron integrals '
        CALL WRTMT_LU(H1DIA,1,NACOB,1,NACOB)
      END IF
*
      RETURN
      END
