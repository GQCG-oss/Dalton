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
      SUBROUTINE ADVICE_SIGMA(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,LADVICE)
*
* Advice Sigma routine about best route to take
*
* ITERM : Term  to be studied :
*         =1 alpha-beta term
*         ....... ( to be continued )
*
* LADVICE : ADVICE given ( short, an integer !!)
*
* For ITERM = 1 :
*           LADVICE = 1 : Business as usual, no transpose of matrix
*                         (resolution on alpha strings, direct exc 
*                         on beta)
*           LADVICE = 2 = Transpose matrices
*                         (resolution on beta strings, direct exc on 
*                         alpha)
*
* Jeppe Olsen, Tirstrup Airport, Jan 12, 98
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "gasstr.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "crun.inc"
*. Specific input
      INTEGER IAOCC(*),IBOCC(*),JAOCC(*),JBOCC(*)
*. Local Scratch
       DIMENSION ITP(16),JTP(16),KTP(16),LTP(16)
*
      NTEST = 00
      IF(ITERM.EQ.1) THEN
*.
*. sigma(i,Ka,Ib) = sum(i,kl)<Ib!Eb_kl!Jb>(ij!kl)C(j,Ka,Jb)
*
* Number of ops : Number of sx(kl) N_i*N_j_dimension of C(j,Ka,Jb)
*.No absolute calc of flops is made, only a relative measure
*
* Single excitations connecting the two types
*
C            SXTYP2_GAS(NSXTYP,ITP,JTP,NGAS,ILTP,IRTP,IPHGAS)
        CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOCC,JAOCC,IPHGAS)
        CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOCC,JBOCC,IPHGAS)
C?      WRITE(6,*) 'NIJTYP, NKLTYP', NIJTYP,NKLTYP
*. P-h modifications ( I cannot predict these at the moment
        IF(NIJTYP.GE.1.AND.NKLTYP.GE.1) THEN
*
        IF((IPHGAS(ITP(1)).EQ.2.AND.IPHGAS(JTP(1)).EQ.2).OR.
     &     (IPHGAS(KTP(1)).EQ.2.AND.IPHGAS(LTP(1)).EQ.2)     ) THEN
           IPHMODI = 1
         ELSE
           IPHMODI = 0
         END IF
        ELSE
           IPHMODI = 0
        END IF

*
        IF(IPHMODI.EQ.1.OR.NIJTYP.NE.1.OR.NKLTYP.NE.1
     &     .OR.IADVICE.EQ.0) THEN
*. Several connections, i.e. the alpha or the beta blocks are identical,
*. or ph modifications
*. just continue
          LADVICE = 1
        ELSE
* =========================================
*.. Index for flops along C(j,Ka,Jb) route
* =========================================
*.Dim of C(j,Ka,Jb) relative to C(Ja,Jb)
*. going from Ja to  Ka reduces occ by one elec, changes dim by n/(N-n+1)
          XNJOB = FLOAT(NOBPT(JTP(1)))
          XNJEL = FLOAT(JAOCC(JTP(1)))
          XCJKAJB = XNJOB*XNJEL/(XNJOB-XNJEL+1)
*. Number of kl excitations per beta string :
          XNKLSX = FLOAT((NOBPT(KTP(1))-JBOCC(KTP(1)))*JBOCC(LTP(1)))
*. Number of ops (relative to dim of C)
          XNIOB = FLOAT(NOBPT(ITP(1)))
          XFLOPA = XCJKAJB*XNKLSX*XNIOB
* =========================================
*.. Index for flops along C(l,Ja,Kb) route
* =========================================
*.Dim of C(l,Ja,Kb) relative to C(Ja,Jb)
          XNLOB = FLOAT(NOBPT(LTP(1)))
          XNLEL = FLOAT(JBOCC(LTP(1)))
          XCLJAKB = XNLOB*XNLEL/(XNLOB-XNLEL+1)
*. Number of ij excitations per alpha string :
          XNIJSX = FLOAT((NOBPT(ITP(1))-JAOCC(ITP(1)))*JAOCC(JTP(1)))
*. Number of ops (relative to dim of C)
          XNKOB = FLOAT(NOBPT(KTP(1)))
          XFLOPB = XCLJAKB*XNIJSX*XNKOB
*. Switch to second route if atleast 20 percent less work
          IF(XFLOPB.LE.0.8*XFLOPA) THEN
            LADVICE = 2
          ELSE
            LADVICE = 1
          END IF
*. Well, an additional consideration :
* If the C block involes the smallest allowed number of elecs in hole space,
* and the annihilation is in hole space
* then we do the annihilation in the space with the smallest number of
* hole electrons.
          LHOLEA =0
          LHOLEB =0
          DO IGAS = 1, NGAS
            IF(IPHGAS(IGAS).EQ.2) THEN
              LHOLEA = LHOLEA + JAOCC(IGAS)
              LHOLEB = LHOLEB + JBOCC(IGAS)
            END IF
          END DO
*
          IF(LHOLEA+LHOLEB.EQ.MNHL.AND.
     &       (IPHGAS(JTP(1)).EQ.2.OR.IPHGAS(LTP(1)).EQ.2))  THEN
*
             if (IPHGAS(JTP(1)).eq.2) then
               KHOLEA = LHOLEA - 1
               KHOLEB = LHOLEB
             else
               KHOLEA = LHOLEA
               KHOLEB = LHOLEB - 1
             end if
*
             IF(KHOLEA.EQ.KHOLEB) THEN
               LLADVICE = LADVICE
             ELSE IF(KHOLEA.LT.KHOLEB) THEN
               LLADVICE= 1
             ELSE
               LLADVICE = 2
             END IF
             IF(NTEST.GE.100.AND.LADVICE.NE.LLADVICE) THEN
               WRITE(6,*) ' Advice changed by hole considetions'
               WRITE(6,*) ' LADVICE, LLADVICE', LADVICE,LLADVICE
             END IF
             LADVICE = LLADVICE
          END IF
*
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' ADVICE active '
            WRITE(6,*) ' IAOCC JAOCC IBOCC JBOCC'
            CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
            CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
            WRITE(6,*) ' ITP JTP KTP LTP ',ITP(1),JTP(1),KTP(1),LTP(1)
            WRITE(6,*) ' XFLOPA,XFLOPB', XFLOPA,XFLOPB
            WRITE(6,*) ' ADVICE given : ', LADVICE
          END IF
        END IF
*       ^ End if several types/ph modi
      END IF
*     ^ End if ITERM test ( type of excitation)
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ALG_ROUTERX(ILOC,IROC,NOP,IOP_TYP,IOP_AC,IOP_REO,
     &                       SIGN)
*
* Decide route for calculating <ILOC! operator string !ROC>
* with smallest amount of operations/storage
*
* The operator string contains NOP elementary operators
* defined by
*
* IOP_TYP : Orbital type of each operator
* IOP_AC  : Creation/annihilation operator
*
*. Output :
* ==========
*
* IOP_REO : New to old order
* sign    : sign of full rank operator
*
* Method : use IPHGAS , and move creation of holes
*                       and annihilation of particles to the left
*
* Jeppe Olsen, October 1997
* version  of : dec 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "cgas.inc"
*. Specific input, I*OC : Contains number of electrons in each space
*  (not used !!!)
       INTEGER ILOC(NGAS),IROC(NGAS)
       INTEGER IOP_TYP(*), IOP_AC(*)
*. Output, new to old order :
       INTEGER IOP_REO(*)
*
*. The loops are assumed to be of the resolution type, so
*. we will insert a resolution of the identity after NOP/2
*  operators(we assume that NOP is even ) :
*
* <ILOC! IOPL !K><K! IOPR ! IROC >
*=<ILOC! IOPL !K><IROC! (IOPR)+ ! K >
*
*
* and the question is now how to order the operators.
* The operation count can in general be mininized by
*. According to above : Should operators stand left or
*                       right of the resolution K
*
      NTEST = 00
*. First time around : not very general
      NMOVE = 0
      IF(NOP.EQ.2) THEN
*
* One-Electron operator
*
        IF(IOP_AC(1).EQ.2.AND.IOP_AC(2).EQ.1) THEN
* a+1 a 2 : Move around only if both spaces are hole spaces
          NMOVE = 0
          IF(IPHGAS(IOP_TYP(1)).EQ.2) NMOVE = 1
          IF(IPHGAS(IOP_TYP(2)).EQ.2) NMOVE = NMOVE + 1
*
          IF(NTEST.GE.100) THEN
            WRITE(6,*) '  IOP_TYP(1), IOP_TYP(2), NMOVE ',
     &                    IOP_TYP(1), IOP_TYP(2), NMOVE
          END IF
        ELSE
          WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
          WRITE(6,*) ' NOP  = ', NOP
          WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
          Call Abend2( ' ALG_ROUTER : Path not implemented ' )
        END IF
*
        IF(NMOVE.EQ.2) THEN
*. Move the operators around in expansion
          IOP_REO(1) = 2
          IOP_REO(2) = 1
          SIGN = -1.0D0
        ELSE
*. No reorganization adviced
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          SIGN = 1.0D0
        END IF
      ELSE IF(NOP.EQ.4) THEN
*
* Two-electron operator
        IF(IOP_AC(1).EQ.2.AND.IOP_AC(2).EQ.2.AND.
     &     IOP_AC(3).EQ.1.AND.IOP_AC(4).EQ.1) THEN
* a+1 a+2 a3 a4
          IF(IPHGAS(IOP_TYP(1)).EQ.2) THEN
            IMOVE1 = 1
          ELSE
            IMOVE1 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(2)).EQ.2) THEN
            IMOVE2 = 1
          ELSE
            IMOVE2 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(3)).EQ.2) THEN
            IMOVE3 = 1
          ELSE
            IMOVE3 = 0
          END IF
*
          IF(IPHGAS(IOP_TYP(4)).EQ.2) THEN
            IMOVE4 = 1
          ELSE
            IMOVE4 = 0
          END IF
        ELSE
*. Not yet implemented
          WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
          WRITE(6,*) ' NOP  = ', NOP
          WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
          Call Abend2( ' ALG_ROUTER : Path not implemented ' )
        END IF
*. Number of left and right operators that would like to be moved
        IMOVEL = IMOVE1 + IMOVE2
        IMOVER = IMOVE3 + IMOVE4
        IMOVE = MIN(IMOVEL,IMOVER)
        IF(IMOVE.EQ.2) THEN
*. Well, we are all moving
          IOP_REO(1) = 3
          IOP_REO(2) = 4
          IOP_REO(3) = 1
          IOP_REO(4) = 2
          SIGN = 1.0D0
        ELSE IF(IMOVE.EQ.0) THEN
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          IOP_REO(3) = 3
          IOP_REO(4) = 4
          SIGN = 1.0D0
        ELSE IF(IMOVEL.EQ.1.AND.IMOVER.EQ.1) THEN
* a+ a a+ a
* 1 : Position the two operators to be moved as operators 2, 3
          IOP_REO(1) = 1
          IOP_REO(2) = 2
          IOP_REO(3) = 3
          IOP_REO(4) = 4
          SIGN = 1.0D0
          IF(IMOVE1.EQ.1) THEN
            IOP_REO(1) = 2
            IOP_REO(2) = 1
            SIGN = -SIGN
          END IF
          IF(IMOVE4.EQ.1) THEN
           IOP_REO(3) = 4
           IOP_REO(4) = 3
           SIGN = - SIGN
          END IF
* 2 : and interchange operators 2 and 3
          ISAVE = IOP_REO(2)
          IOP_REO(2) = IOP_REO(3)
          IOP_REO(3) = ISAVE
          SIGN = - SIGN
        ELSE IF (IMOVEL.EQ.2.AND.IMOVER.EQ.1) THEN
* Final operator is a a+ a+ a
          IF(IMOVE4.EQ.1) THEN
            IOP_REO(3) = 4
            IOP_REO(4) = 3
            SIGN = -1.0D0
          ELSE
            IOP_REO(3) = 3
            IOP_REO(4) = 4
            SIGN = 1.0D0
          END IF
          IOP_REO(1) = IOP_REO(3)
          IOP_REO(2) = 1
          IOP_REO(3) = 2
          SIGN = SIGN
        ELSE IF( IMOVEL.EQ.1 .AND. IMOVER .EQ. 2 ) THEN
* Final operator is a+ a a a+
          IF( IMOVE1.EQ.1) THEN
            IOP_REO(1) = 2
            IOP_REO(2) = 1
            SIGN = -1.0D0
          ELSE
            IOP_REO(1) = 1
            IOP_REO(2) = 2
            SIGN = 1.0D0
          END IF
          IOP_REO(4) = IOP_REO(2)
          IOP_REO(2) = 3
          IOP_REO(3) = 4
          SIGN = SIGN
        END IF
*       ^ End of number of pairs to be moved around
      ELSE
*. Not one- or two- electron operator
        WRITE(6,*) ' ALG_ROUTER : Path not yet implemented : '
        WRITE(6,*) ' NOP  = ', NOP
        WRITE(6,*) ' IOP_AC = ', (IOP_AC(II),II=1, NOP)
        Call Abend2( ' ALG_ROUTER : Path not implemented ' )
      END IF
*

      IF(NOP.EQ.4.AND.NTEST.GE.100) THEN
        WRITE(6,*) ' Information from ALG_ROUTER '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*)
        WRITE(6,*)
     &  ' Input : anni- or crea- operators + orbital types'
        CALL IWRTMA(IOP_AC,1,NOP,1,NOP)
        CALL IWRTMA(IOP_TYP,1,NOP,1,NOP)
        WRITE(6,*)
        WRITE(6,*) ' Suggested output order '
        CALL IWRTMA(IOP_REO,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE CKAJJB_PA(CKAJJB,CKAJJBPA,IWAY,NKA,NJ,NJB,
     &                     IREO,NJB_PAS,NJB_ACT,NSMST,SIGN,
     &                     I_ADD_COPY)
*
* Reform C(Ka,j,Jb) between usual and active/passive division of Jb :
*
* C(Ka,j,Jb) <=> C(Ka,Jb_pa,j,Jb_ac)
*
* IWAY : =1 => normal to passive/active form
* IWAY : =2 => passive/active to normal form

* Jeppe Olsen, March 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CKAJJB(*)
      DIMENSION NJB_PAS(*),NJB_ACT(*),IREO(*)
*. Output
      DIMENSION CKAJJBPA(*)
CT    CALL QENTER('CKA_PA')
*
      NTEST = 000
*
      IJB = 0
      DO JB_ACT_SM = 1, NSMST
*
        IF(JB_ACT_SM.EQ.1) THEN
         IOFF = 1
         IOFF_PA = 1
        ELSE
         IOFF = IOFF + NJB_PAS(JB_ACT_SM-1)*NJB_ACT(JB_ACT_SM-1)
         IOFF_PA = IOFF_PA
     &           + NJB_ACT(JB_ACT_SM-1)*NJ*NJB_PAS(JB_ACT_SM-1)*NKA
        END IF
*
        LJB_ACT = NJB_ACT(JB_ACT_SM)
        LJB_PAS = NJB_PAS(JB_ACT_SM)
        IB_OUT_ADD = LJB_PAS*NKA
* C(Ka,j,Jb) <=> C(Ka,Jb_pa,j,Jb_ac)
        DO JB_ACT = 1, LJB_ACT
          DO JB_PAS = 1, LJB_PAS
            IJB =  IJB + 1
            IJB_IN = IREO(IJB)
            IB_IN = (IJB_IN-1)*NJ*NKA-NKA
            IB_OUT = (JB_ACT-1)*NJ*LJB_PAS*NKA
     &             + (JB_PAS-1)           *NKA + IOFF_PA - 1
     &             - IB_OUT_ADD
            DO J = 1, NJ
*. Address of C(1,j,Jb)
C             IB_IN = (IJB_IN-1)*NJ*NKA+(J-1)*NKA
              IB_IN = IB_IN + NKA
*. Address of C(1,Jb_pa,j,Jb_ac)
              IB_OUT = IB_OUT + IB_OUT_ADD
C             IB_OUT = (JB_ACT-1)*NJ*LJB_PAS*NKA
C    &               + (J     -1)   *LJB_PAS*NKA
C    &               + (JB_PAS-1)           *NKA + IOFF_PA - 1
C             WRITE(6,*) 'IB_IN, IB_OUT', IB_IN, IB_OUT
C             WRITE(6,*) ' IJB_IN, NJA,NKA ',  IJB_IN, NJA,NKA
              IF(IWAY.EQ.1) THEN
*. Normal => passive/active form
                DO KA = 1, NKA
                  CKAJJBPA(IB_OUT+KA) = CKAJJB(IB_IN+KA)
                END DO
              ELSE IF(IWAY.EQ.2) THEN
C               IF(I_ADD_COPY.EQ.1) THEN
*. Passive/Active => Normal form
                  DO KA = 1, NKA
                    CKAJJB(IB_IN+KA) =  CKAJJB(IB_IN+KA) +
     &              SIGN*CKAJJBPA(IB_OUT+KA)
                  END DO
C               ELSE IF (I_ADD_COPY.EQ.2) THEN
C                 DO KA = 1, NKA
C                   CKAJJB(IB_IN+KA) =
C    &              SIGN*CKAJJBPA(IB_OUT+KA)
C                 END DO
C               END IF
*               ^ End of ADD_COPY switch
              END IF
*             ^ End of switch passive/active <=  => normal
            END DO
*           ^ End of loop over J-orbitals
          END DO
*         ^ End of loop over JB_pa
        END DO
*       ^ End of loop over JB_ac
      END DO
*     ^ End of loop over symmetry of active beta strings
*
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' =================='
       WRITE(6,*) ' CKAJJB_PA speaking '
       WRITE(6,*) ' =================='
       WRITE(6,*)
       IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' Normal to passive/active form '
       ELSE
         WRITE(6,*) ' passive/active  to normal form '
       END IF
*
       WRITE(6,*) ' C(ka,j,Jb) as C(kaj,Jb) '
       NKAJ = NKA*NJ
       CALL WRTMT_LU(CKAJJB,NKAJ,NJB,NKAJ,NJB)
*
       WRITE(6,*) ' C(Ka,Jb_pa,j,Jb_ac) as blocks C(KaJB_paj,Jb_ac)'
       DO JB_ACT_SM = 1, NSMST
         WRITE(6,*) ' Symmetry of Jb_act ', JB_ACT_SM
         IF(JB_ACT_SM.EQ.1) THEN
          IOFF_PA = 1
         ELSE
          IOFF_PA = IOFF_PA
     &            + NJB_ACT(JB_ACT_SM-1)*NJ*NJB_PAS(JB_ACT_SM-1)*NKA
         END IF
         NROW = NKA*NJB_PAS(JB_ACT_SM)*NJ
         NCOL = NJB_ACT(JB_ACT_SM)
*
C        WRITE(6,*)  ' NJB_ACT(JB_ACT_SM)', NJB_ACT(JB_ACT_SM)
C        WRITE(6,*)  ' NJB_PAS(JB_ACT_SM)', NJB_PAS(JB_ACT_SM)
C        WRITE(6,*)  ' NJ and NKA ', NJ,NKA
C        WRITE(6,*) ' Offset IOFF_PA', IOFF_PA
         CALL WRTMT_LU(CKAJJBPA(IOFF_PA),NROW,NCOL,NROW,NCOL)
       END DO
      END IF
*
CT    CALL QEXIT('CKA_PA')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIHDJ_LUCI(IASTR,IBSTR,NIDET,
     &                 JASTR,JBSTR,NJDET,
     &                 NAEL,NBEL,
     &IWORK,LWORK,NORB,ONEBOD,HAMIL,ISYM,NINOB,ECORE,ICOMBI,PSIGN,
     & IASTRM,IBSTRM,JASTRM,JBSTRM,
     & IGENSG,IASGN,IBSGN,JASGN,JBSGN,LIA,LIB,NDIF0,NDIF1,NDIF2,
     & IPRT)
*
* A set of left hand side determinants defined by string numbers
* IASTR and IBSTR and a set of right hand side determinants
* defined by JASTR and JBSTR are given.
*
* Obtain Hamiltonian matrix  < IA IB ! H ! JA JB >
*
* If Icombi .NE. 0 Spin combinations are assumed  for alpha and
* beta strings with different orbital configurations
*   1/SQRT(2) * ( !I1A I2B! + PSIGN * !I2A I1B! )
*
* If ISYM .EQ. 0 FULL Hamiltonian is constructed
* If ISYM .NE. 0 LOWER half of hamiltonian is constructed
*
* JEPPE OLSEN JANUARY 1989
*
*. Modifed to work with string numbers instead of strings
*. March 93
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION IASTR(*),IBSTR(*)
      DIMENSION JASTR(*),JBSTR(*)
      DIMENSION IASTRM(NAEL,*),IBSTRM(NBEL,*)
      DIMENSION JASTRM(NAEL,*),JBSTRM(NBEL,*)
      DIMENSION IASGN(*),IBSGN(*),JASGN(*),JBSGN(*)
*
      DIMENSION IWORK(*), HAMIL(*), ONEBOD(NORB,NORB)
      DIMENSION LIA(NAEL),LIB(NBEL)
*
      NTESTL = 0
      NTEST = MAX(NTESTL,IPRT)
      IF( NTEST .GE. 5 ) THEN
         WRITE(6,*) ' =========================='
         WRITE(6,*) ' DIHDJ reporting to service'
         WRITE(6,*) ' =========================='
      END IF
*
      IF(NTEST.GE. 200) THEN
        WRITE(6,*) ' Input alpha strings,IASTR,JASTR '
        CALL IWRTMA(IASTR,1,NIDET,1,NIDET)
        CALL IWRTMA(JASTR,1,NJDET,1,NJDET)
        WRITE(6,*) ' Input beta  strings,IBSTR,JBSTR '
        CALL IWRTMA(IBSTR,1,NIDET,1,NIDET)
        CALL IWRTMA(JBSTR,1,NJDET,1,NJDET)
      END IF
*
*. Scratch space : 4 vectors of length NORB
      KLFREE = 1
      KLIAE  = KLFREE
      KLFREE = KLIAE + NORB
      KLIBE  = KLFREE
      KLFREE = KLIBE + NORB
*
      KLJAE = KLFREE
      KLFREE = KLJAE + NORB
      KLJBE = KLFREE
      KLFREE = KLJBE + NORB
*
      IF( ISYM .EQ. 0 ) THEN
        LHAMIL = NIDET*NJDET
      ELSE
        LHAMIL = NIDET*(NIDET+1) / 2
      END IF
      CALL SETVEC(HAMIL,0.0D0,LHAMIL)
*
      NTERMS= 0
      NDIF0 = 0
      NDIF1 = 0
      NDIF2 = 0
*
*. Loop over J determinants
*
      DO 1000 JDET = 1,NJDET
        IF( NTEST .GE. 100 ) WRITE(6,*) '  ****** 1000 JDET ', JDET
* Expand JDET
        JASTAC =JASTR(JDET)
        JBSTAC =JBSTR(JDET)
*
        IF(IGENSG .GT. 0 ) THEN
         JXSGN = JASGN(JASTAC)*JBSGN(JBSTAC)
        ELSE
         JXSGN = 1
        END IF
*
        CALL ISETVC(IWORK(KLJAE),0,NORB)
        CALL ISETVC(IWORK(KLJBE),0,NORB)
        DO 40 IAEL = 1, NAEL
          IWORK(KLJAE-1+JASTRM(IAEL,JASTAC) ) = 1
   40   CONTINUE
*
        DO 50 IBEL = 1, NBEL
          IWORK(KLJBE-1+JBSTRM(IBEL,JBSTAC) ) = 1
   50   CONTINUE
*
        IF( ICOMBI .NE. 0 ) THEN
          IF(JASTAC .EQ. JBSTAC) THEN
             JAEQJB = 1
          ELSE
             JAEQJB = 0
          END IF
        END IF
*
*
        IF( NTEST .GE. 100 ) THEN
          WRITE(6,*) ' LOOP 1000 JDET =  ',JDET
          WRITE(6,*) ' JASTAC AND JBSTAC ', JASTAC,JBSTAC
          WRITE(6,*) ' Expanded ALPHA and BETA string '
          CALL IWRTMA(IWORK(KLJAE),1,NORB,1,NORB)
          CALL IWRTMA(IWORK(KLJBE),1,NORB,1,NORB)
        END IF
*
        IF( ISYM .EQ. 0 ) THEN
          MINI = 1
        ELSE
          MINI = JDET
        END IF
*
* Loop over I determinants
*
        DO 900 IDET = MINI, NIDET
          IF(NTEST.GE. 100 )
     &    WRITE(6,*) '   LOOP 900 IDET .... ',IDET
          IASTAC = IASTR(IDET)
          IBSTAC = IBSTR(IDET)
*
          IF(IGENSG .GT. 0 ) THEN
           IXSGN = IASGN(IASTAC)*IBSGN(IBSTAC)
          ELSE
           IXSGN = 1
          END IF
*
          IF(IASTAC.EQ.IBSTAC) THEN
            IAEQIB = 1
          ELSE
            IAEQIB = 0
          END IF

*
          IF(ICOMBI.EQ.1 .AND. IAEQIB+JAEQJB.EQ.0 ) THEN
              NLOOP = 2
          ELSE
              NLOOP = 1
          END IF
C
          DO 899 ILOOP = 1, NLOOP
           NTERMS = NTERMS + 1
* For second part of spin combinations strings should be swopped
           IF(ILOOP.EQ.1) THEN
             CALL ICOPVE(IASTRM(1,IASTAC),LIA,NAEL)
             CALL ICOPVE(IBSTRM(1,IBSTAC),LIB,NBEL)
           ELSE IF (ILOOP.EQ.2) THEN
             CALL ICOPVE(IASTRM(1,IASTAC),LIB,NAEL)
             CALL ICOPVE(IBSTRM(1,IBSTAC),LIA,NBEL)
           END IF
*
* ==============================
*. Number of orbital differences
* ==============================
*
           NACM = 0
           DO 61 IAEL = 1, NAEL
             NACM = NACM + IWORK(KLJAE-1+LIA(IAEL))
   61      CONTINUE
C
           NBCM = 0
           DO 62 IBEL = 1, NBEL
             NBCM = NBCM + IWORK(KLJBE-1+LIB(IBEL))
   62      CONTINUE
C
           NADIF = NAEL-NACM
           NBDIF = NBEL-NBCM
           IF( NTEST .GE. 100 ) THEN
             WRITE(6,*) '  LOOP 900 IDET ',IDET
             WRITE(6,*) ' Comparison : NADIF ,NBDIF ', NADIF,NBDIF
           END IF
*
           IF(NADIF+NBDIF .GT. 2 ) GOTO 898
*. Factor for combinations
           IF( ICOMBI .EQ. 0 ) THEN
             CONST = 1.0D0
           ELSE
             IF((JAEQJB +IAEQIB) .EQ.2 ) THEN
               CONST = 1.0D0
             ELSE IF( (JAEQJB+IAEQIB) .EQ. 1 ) THEN
               CONST = 1.0D0/SQRT(2.0D0)*(1.0D0+PSIGN)
              ELSE IF( (JAEQJB+IAEQIB) .EQ. 0 ) THEN
               IF( ILOOP .EQ. 1)  THEN
                 CONST = 1.0D0
               ELSE
                 CONST = PSIGN
               END IF
             END IF
           END IF
*. External sign factor
           IF(IXSGN*JXSGN .EQ. -1 ) CONST = - CONST
           IF(NTEST.GE.100) WRITE(6,*) ' CONST ', CONST
*
* ==================================================
*.. Find differing orbitals and sign for permutation
* ==================================================
*
* Expand idet
           CALL ISETVC(IWORK(KLIAE),0,NORB)
           CALL ISETVC(IWORK(KLIBE),0,NORB)
*
           DO 42 IAEL = 1, NAEL
             IWORK(KLIAE-1+LIA(IAEL)) = 1
   42      CONTINUE
*
           DO 52 IBEL = 1, NBEL
             IWORK(KLIBE-1+LIB(IBEL) ) = 1
   52      CONTINUE
*
*. One pair of differing alpha electrons
*
           IF(NADIF .EQ. 1 ) THEN
             DO 120 IAEL = 1,NAEL
               IF(IWORK(KLJAE-1+LIA(IAEL)).EQ.0) THEN
                 IA = LIA(IAEL)
                 IEL1 = IAEL
                 GOTO 121
               END IF
  120        CONTINUE
  121        CONTINUE
*
             DO 130 JAEL = 1,NAEL
               IF(IWORK(KLIAE-1+JASTRM(JAEL,JASTAC)).EQ.0) THEN
                 JA = JASTRM(JAEL,JASTAC)
                 JEL1 = JAEL
                 GOTO 131
                END IF
  130        CONTINUE
  131        CONTINUE
             SIGNA = (-1)**(JEL1+IEL1)
             IF(NTEST.GE.100)
     &       WRITE(6,*) ' IA JA SIGNA... ',IA,JA,SIGNA
           END IF
*
*. One pair of differing beta electrons
*
           IF(NBDIF .EQ. 1 ) THEN
             DO 220 IBEL = 1,NBEL
               IF(IWORK(KLJBE-1+LIB(IBEL) ).EQ.0) THEN
                 IB = LIB(IBEL)
                 IEL1 = IBEL
                 GOTO 221
                END IF
  220        CONTINUE
  221        CONTINUE
C
             DO 230 JBEL = 1,NBEL
               IF(IWORK(KLIBE-1+JBSTRM(JBEL,JBSTAC)).EQ.0) THEN
                 JB = JBSTRM(JBEL,JBSTAC)
                 JEL1 = JBEL
                 GOTO 231
                END IF
  230        CONTINUE
  231        CONTINUE
             SIGNB = (-1)**(JEL1+IEL1)
             IF(NTEST.GE.100)
     &       WRITE(6,*) ' IB JB SIGNB... ',IB,JB,SIGNB
           END IF
*
*. Two pairs of differing alpha electrons
*
           IF(NADIF .EQ. 2 ) THEN
             IDIFF = 0
             DO 320 IAEL = 1,NAEL
               IF(IWORK(KLJAE-1+LIA(IAEL)          ).EQ.0) THEN
                 IF( IDIFF .EQ. 0 ) THEN
                   IDIFF = 1
                   I1 = LIA(IAEL)
                   IPERM = IAEL
                 ELSE
                   I2 = LIA(IAEL)
                   IPERM = IAEL + IPERM
                   GOTO 321
                 END IF
               END IF
  320        CONTINUE
  321        CONTINUE
*
             JDIFF = 0
             DO 330 JAEL = 1,NAEL
               IF(IWORK(KLIAE-1+JASTRM(JAEL,JASTAC)).EQ.0) THEN
                 IF( JDIFF .EQ. 0 ) THEN
                   JDIFF = 1
                   J1 = JASTRM(JAEL,JASTAC)
                   JPERM = JAEL
                 ELSE
                   J2 = JASTRM(JAEL,JASTAC)
                   JPERM = JAEL + JPERM
                   GOTO 331
                 END IF
               END IF
  330        CONTINUE
  331        CONTINUE
             SIGN = (-1)**(IPERM+JPERM)
           END IF
*
*. Two pairs of differing beta electrons
*
           IF(NBDIF .EQ. 2 ) THEN
             IDIFF = 0
             DO 420 IBEL = 1,NBEL
               IF(IWORK(KLJBE-1+LIB(IBEL)          ).EQ.0) THEN
                 IF( IDIFF .EQ. 0 ) THEN
                   IDIFF = 1
                   I1 = LIB(IBEL)
                   IPERM = IBEL
                 ELSE
                   I2 = LIB(IBEL)
                   IPERM = IBEL + IPERM
                   GOTO 421
                  END IF
               END IF
  420        CONTINUE
  421        CONTINUE
*
             JDIFF = 0
             DO 430 JBEL = 1,NBEL
               IF(IWORK(KLIBE-1+JBSTRM(JBEL,JBSTAC)).EQ.0) THEN
                 IF( JDIFF .EQ. 0 ) THEN
                   JDIFF = 1
                   J1 = JBSTRM(JBEL,JBSTAC)
                   JPERM = JBEL
                 ELSE
                   J2 = JBSTRM(JBEL,JBSTAC)
                   JPERM = JBEL + JPERM
                   GOTO 431
                 END IF
               END IF
  430        CONTINUE
  431        CONTINUE
             SIGN = (-1)**(IPERM+JPERM)
           END IF
*
* =======================
* Value of matrix element
* =======================
*
        IF( NADIF .EQ. 2 .OR. NBDIF .EQ. 2 ) THEN
* 2 differences in alpha or beta strings
          NDIF2 = NDIF2 + 1
* SIGN * (I1 J1 ! I2 J2 ) - ( I1 J2 ! I2 J1 )
          XVAL = SIGN*( GTIJKL(I1,J1,I2,J2)-GTIJKL(I1,J2,I2,J1) )
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 1 ) THEN
*. 1 difference in alpha strings and one difference in beta string
          NDIF2 = NDIF2 + 1
* SIGN * (IA JA ! IB JB )
          XVAL = SIGNA*SIGNB* GTIJKL(IA,JA,IB,JB)
        ELSE IF( NADIF .EQ. 1 .AND. NBDIF .EQ. 0 .OR.
* 1 differences in alpha or beta strings
     &           NADIF .EQ. 0 .AND. NBDIF .EQ. 1 )THEN
          NDIF1 = NDIF1 + 1
* SIGN *
*(  H(I1 J1 ) +
*  (SUM OVER ORBITALS OF BOTH      SPIN TYPES  ( I1 J1 ! JORB JORB )
* -(SUM OVER ORBITALS OF DIFFERING SPIN TYPE   ( I1 JORB ! JORB J1 ) )
          IF( NADIF .EQ. 1 ) THEN
            I1 = IA
            J1 = JA
            SIGN = SIGNA
          ELSE
            I1 = IB
            J1 = JB
            SIGN = SIGNB
          END IF
*
          XVAL = GETH1I(I1,J1)
          DO 520 JAEL = 1, NAEL
            JORB = JASTRM(JAEL,JASTAC)
            XVAL = XVAL + GTIJKL(I1,J1,JORB,JORB)
C?     write(6,*) ' I1 J1 JORB JORB ', I1, J1, JORB, JORB
C?     write(6,*) ' Integral : ', GTIJKL(I1,J1,JORB,JORB)
C?     write(6,*) 'XVAL ', xval
  520     CONTINUE
          DO 521 JBEL = 1, NBEL
            JORB = JBSTRM(JBEL,JBSTAC)
            XVAL = XVAL + GTIJKL(I1,J1,JORB,JORB)
C?     write(6,*) ' I1 J1 JORB JORB ', I1, J1, JORB, JORB
C?     write(6,*) ' Integral : ', GTIJKL(I1,J1,JORB,JORB)
C?     write(6,*) 'XVAL ', xval
  521     CONTINUE
          IF( NADIF .EQ. 1 ) THEN
            DO 522 JAEL = 1, NAEL
              JORB = JASTRM(JAEL,JASTAC)
              XVAL = XVAL - GTIJKL(I1,JORB,JORB,J1)
C?     write(6,*) ' I1  JORB JORB J1', I1,  JORB, JORB, J1
C?     write(6,*) ' Integral : ', GTIJKL(I1,JORB,JORB,J1)
C?     write(6,*) 'XVAL ', xval
  522       CONTINUE
          ELSE
            DO 523 JBEL = 1, NBEL
              JORB = JBSTRM(JBEL,JBSTAC)
              XVAL = XVAL - GTIJKL(I1,JORB,JORB,J1)
C?     write(6,*) ' I1  JORB JORB J1', I1,  JORB, JORB, J1
C?     write(6,*) ' Integral : ', GTIJKL(I1,JORB,JORB,J1)
C?     write(6,*) 'XVAL ', xval
  523       CONTINUE
          END IF
          XVAL = XVAL * SIGN
        ELSE IF( NADIF .EQ. 0 .AND. NBDIF .EQ. 0 ) THEN
*. Diagonal elements
          NDIF0 = NDIF0 + 1
C SUM(I,J OF JDET) H(I,J) + (I I ! J J ) - (I J ! J I )
C
          XVAL = ECORE
          DO 650 IAB = 1, 2
            IF(IAB .EQ. 1 ) THEN
              NIABEL = NAEL
            ELSE
              NIABEL = NBEL
            END IF
            DO 640 JAB = 1, 2
              IF(JAB .EQ. 1 ) THEN
                NJABEL = NAEL
              ELSE
                NJABEL = NBEL
              END IF
              DO 630 IEL = 1, NIABEL
                IF( IAB .EQ. 1 ) THEN
                  IORB = LIA(IEL)
                ELSE
                  IORB = LIB(IEL)
                END IF
                IF(IAB .EQ. JAB ) XVAL = XVAL + GETH1I(IORB,IORB)
C?     write(6,*) ' XVAL ', XVAL
                DO 620 JEL = 1, NJABEL
                  IF( JAB .EQ. 1 ) THEN
                    JORB = LIA(JEL)
                  ELSE
                    JORB = LIB(JEL)
                  END IF
                  XVAL = XVAL + 0.5D0*GTIJKL(IORB,IORB,JORB,JORB)
*. test
C?        FAC = GTIJKL(IORB,IORB,JORB,JORB)
C?        write(6,*) ' IORB JORB (IORB IORB ! JORB JORB )'
C?        WRITE(6,*) IORB,JORB,FAC
C?     write(6,*) ' XVAL ', XVAL

                  IF( IAB . EQ. JAB )
     &            XVAL = XVAL - 0.5D0*GTIJKL(IORB,JORB,JORB,IORB)
*. test
          FAC = GTIJKL(IORB,JORB,JORB,IORB)
C?        write(6,*) ' IORB JORB (IORB JORB ! JORB IORB )'
C?        WRITE(6,*) IORB,JORB,FAC
C?     write(6,*) ' XVAL ', XVAL
  620           CONTINUE
  630         CONTINUE
  640       CONTINUE
  650     CONTINUE
        END IF

        IF( NTEST .GE. 100 )
     &  WRITE(6,*) ' CONST XVAL  ', CONST ,XVAL
        IF( ISYM .EQ. 0 ) THEN
          HAMIL((JDET-1)*NIDET+IDET) =
     &    HAMIL((JDET-1)*NIDET+IDET) + CONST * XVAL
        ELSE
          HAMIL((IDET-1)*IDET/2 + JDET ) =
     &    HAMIL((IDET-1)*IDET/2 + JDET ) + CONST * XVAL
        END IF
  898 CONTINUE
  899 CONTINUE
  900 CONTINUE
 1000 CONTINUE

      IF( IPRT .GT. 0 ) THEN
      WRITE(6,*)
     &'  Number of elements differing by 0 excitation.. ',NDIF0

      WRITE(6,*)
     &'  Number of elements differing by 1 excitation.. ',NDIF1

      WRITE(6,*)
     &'  Number of elements differing by 2 excitation.. ',NDIF2

      WRITE(6,*)
     &'  Number of vanishing elments                    ',
     &   NTERMS - NDIF0 - NDIF1 - NDIF2
      END IF
      IF( IPRT .GE. 10) THEN
        WRITE(6,*) '  Hamiltonian matrix '
        IF( ISYM .EQ. 0 ) THEN
          CALL WRTMT_LU(HAMIL,NIDET,NJDET,NIDET,NJDET)
        ELSE
          CALL PRSYM(HAMIL,NIDET)
        END IF
      END IF
C
C!    Call Abend2( ' ENFORCED STOP AT END OF DIHDJ ' )
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MV7(C,HC,LUC,LUHC
#ifdef VAR_MPI
     &               ,LUCLIST,LUHCLIST,IBLOCKL,NPARBLOCK,RCCTOS,
     &                IGROUPLIST,IPROCLIST
#endif
     &            )
*
* Outer routine for sigma vector generation
* GAS version !!!!
*
* Written in terms of RASG3/SBLOCK, May 1997
* Works only for ICISTR .ne. 1, i.e. pure memory version
* has been eliminated.
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KSIOIO, KSVST, KLSLBT, KLSLEBT, KLSI1BT, KLSIBT, KSBLTP,
     &          KISCLFAC, KSCALLOC
!               for addressing of WORK
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
*./ORBINP/ : NACOB used
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
C
#ifdef VAR_MPI
      DIMENSION LUCLIST(*),LUHCLIST(*),IBLOCKL(*),NPARBLOCK(*)
      INTEGER RCCTOS(*),IGROUPLIST(*),IPROCLIST(*)
#endif
C
      CALL QENTER('MV7  ')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'MV7   ')

*
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*) ' MV7 does not work for ICISTR = 1'
        WRITE(6,*) ' SWITCH to ICISTR = 2,3 or program'
        Call Abend2(' MV7 does not work for ICISTR = 1')
      END IF
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*     offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*     arrays giving allowed type combinations
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &                IGSOCCX,NOCTPA,NOCTPB,
     &                ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),NELFGP,
     &                MXPNGAS,NGAS,WORK(KSIOIO),IPRDIA)
*     arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
         
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
*     arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',1,'SBLTP ')
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
*     arrays for partitioning of sigma
      NTTS = MXNTTS
      CALL MEMMAN(KLSLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLSLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLSI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLSIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
*     batches  of s vector
      ITTSS_ORD = 2

      IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
         CALL PART_CIV_PAR1(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &                   WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,
     &                   LBLOCK_LUCI,
     &                   WORK(KSIOIO),ISMOST(1,ISSM),
     &                   NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &                   WORK(KLSI1BT),WORK(KLSIBT),0,ITTSS_ORD,
     &                   NPARBLOCK,NUM_BLOCKS2)
#endif
      ELSE ! (LUCI_NMPROC .GT. 1) THEN
         CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &               WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LBLOCK,
     &               WORK(KSIOIO),ISMOST(1,ISSM),
     &               NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &               WORK(KLSI1BT),WORK(KLSIBT),0,ITTSS_ORD)
      END IF ! (LUCI_NMPROC .GT. 1) THEN
C
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
      IF(ICISTR.EQ.1) THEN
       LLUC = 0
       LLUHC = 0
      ELSE
       LLUC = LUC
       LLUHC = LUHC
      END IF
C
C     scaling factor array
C
      CALL MEMMAN(KISCLFAC,NUM_BLOCKS2,'ADDL  ',2,'SSCLFC')
      CALL DZERO(WORK(KISCLFAC),NUM_BLOCKS2)
C
C     ISCLFAC_GROUP array for c-blocks:
C     0 = block is not active, > 0: block is active
C
      CALL MEMMAN(KSCALLOC,NUM_BLOCKS2,'ADDL  ',1,'KITCLL')
      
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*) 'ICISTR = 1 route currently not active'
        Call Abend2('MV7 : ICISTR = 1 route currently not active')
      ELSE
        CALL RASSG3(C,HC,NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &       WORK(KLSI1BT),WORK(KLSIBT),LLUC,LLUHC,
     &       WORK(KISCLFAC),WORK(KSCALLOC)
#ifdef VAR_MPI
     &               ,LUCLIST,LUHCLIST,IBLOCKL,NPARBLOCK,RCCTOS,
     &                IGROUPLIST,IPROCLIST
#endif
     &             )
C            RASSG3(CB,SB,LBATS,LEBATS,I1BATS,IBATS,
C    &       LUC,LUHC)
      END IF
*. Eliminate local memory
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'MV7   ')
*
      CALL QEXIT('MV7  ')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MXRESC(IAB,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                  NSMST,NSTFSMSPGP,
     &                  MXPNSMST,
     &                  NSMOB,MXPTOB,NTPOB,NTSOB,NTESTG,MXPKA,
     &                  NEL1234,
     &                  MXCJ,MXCIJA,
     &                  MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
*
* Find largest dimension of matrix C(Ka,Ib,J)
* Find largest dimension of matrix C(ij,Ka,Ib)
* Find largest dimension of matrix C(ij,Ia,Kb)
* Find largest dimension of matrix C(ij,Ka,Kb)
*
* Largest block of single excitations MXSXBL

*. Input
* IAB :allowed combination of alpha and beta supergroups
* IOCPTA : Number of first active alpha supergroup
* IOCPTB : Number of first active beta  supergroup
* NOCTPA : Number of active alpha supergroups
* NOCTPB : Number of active alpha supergroups

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAB(NOCTPA,NOCTPB)
      DIMENSION NSTFSMSPGP(MXPNSMST,*)
      DIMENSION NTSOB(MXPTOB,NSMOB)
      DIMENSION NEL1234(MXPTOB,*)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.100) WRITE(6,*) ' MXRESC : MXPKA ', MXPKA
*
* matrix C(j,Ka,Ib)
*
*. Note : Only done for alpha-strings, problems when transposing
*         constructing C(J,Ia,Kb)
      MXCJ = 0
      MXADKBLK = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(j,Ka,Ib)'
        write(6,*)
      end if
      DO 100 IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO 200 IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXB = 0
            DO 210 ISM = 1, NSMST
              MXB =MAX(MXB,NSTFSMSPGP(ISM,IBTPABS))
  210       CONTINUE
            IF(NTEST.GE.100) WRITE(6,*) ' MXB = ', MXB
            DO 300 IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO 310 KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
  310           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                MXKAO = MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKA= MXPKA
                MXSOB = 0
                DO 320 ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
  320           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKAO)
                LCJBLK = MXSOB*MXKA*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
  300       CONTINUE
          END IF
  200   CONTINUE
  100 CONTINUE
*
* matrix C(j,Ia,Kb)
*
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(j,Ia,Kb)'
        write(6,*)
      end if
*
      DO IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXA = 0
            DO ISM = 1, NSMST
              MXA =MAX(MXA,NSTFSMSPGP(ISM,IATPABS))
            END DO
            IF(NTEST.GE.100) WRITE(6,*) ' MXA = ', MXA
            DO IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IBTPABS,1,IOBTP,1,KBTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KBTP ',IOBTP,KBTP
              IF(KBTP.GT.0) THEN
                MXKB = 0
                DO KSM = 1, NSMST
                  MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                MXKBO = MXKB
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKB= MXPKA
                MXSOB = 0
                DO ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKBO)
                LCJBLK = MXSOB*MXKB*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXADKBLK,MXCJ ', MXADKBLK,MXCJ
      END IF
*
* matrix C(ij,Ka,Ib)
* both Ka and Ib blocked
*
      MXCIJA = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ka,Ib)'
        write(6,*)
      end if
*
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA -1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIB = 0
            DO  ISM = 1, NSMST
              MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
            END DO
            IF(MXIB.GT.MXPKA) MXIB = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,K1ATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1ATP ',IOBTP,K1ATP
              IF(K1ATP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1ATP,1,JOBTP,1,KATP)
                  IF(KATP.GT.0) THEN
                    MXKA = 0
                    DO KSM = 1, NSMST
                      MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                    IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &              MXKA= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKA*MXIB
                    MXCIJA = MAX(MXCIJA,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJA ', MXCIJA
      END IF
*
*
* matrix C(ij,Ia,Kb)
* both Ka and Ib blocked
*
      IF(NTEST.GE.100) WRITE(6,*) ' MXCIJB under development '
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ia,Kb)'
        write(6,*)
      end if
*
      MXCIJB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIA = 0
            DO  ISM = 1, NSMST
              MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS ))
            END DO
            IF(MXIA.GT.MXPKA) MXIA = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IBTP
              CALL NEWTYP(IBTPABS,1,IOBTP,1,K1BTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1BTP ',IOBTP,K1BTP
              IF(K1BTP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1BTP,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXIA
                    MXCIJB = MAX(MXCIJB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC : MXCIJB ', MXCIJB
      END IF
*
*
* matrix C(ij,Ka,Kb)
* both Ka and Kb blocked
*
      MXCIJAB = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Doing C(ij,Ka,Kb)'
        write(6,*)
      end if
*
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            DO  IOBTP = 1, NTPOB
*. type of Ka string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) MXKA= MXPKA
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from IBTP
                  CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXKA
                    MXCIJAB = MAX(MXCIJAB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
*
* Largest block of single excitations :
* Strings of given type and sym, orbitals of given type and sym
*
* Largest block of creations : a+i !kstring> where K string is
* obtained as single annihilations
      MXSXBL = 0
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Largest block of creations a+i |kstring>'
        write(6,*) 'alpha strings'
        write(6,*)
      end if
*
*. For alpha strings :
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        MXIA = 0
        DO  ISM = 1, NSMST
          MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP(IATPABS,1,JOBTP,1,KATP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KATP ',JOBTP,KATP
          IF(KATP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IATPABS))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP(KATP,2,IOBTP,1,JATP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
*. For beta  strings :
      if (NTEST.ge.100) then
        write(6,*)
        write(6,*) 'Largest block of creations a+i |kstring>'
        write(6,*) 'beta strings'
        write(6,*)
      end if
*
      DO  IBTP = 1, NOCTPB
        IBTPABS = IBTP + IOCTPB - 1
        MXIB = 0
        DO  ISM = 1, NSMST
          MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KBTP ',JOBTP,KBTP
          IF(KBTP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IBTP))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP(KBTP,2,IOBTP,1,JBTP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.00) THEN
        WRITE(6,*) 'MXRESC: MXSXBL : ', MXSXBL
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MXRESCPH(IAB,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                  NSMST,NSTFSMSPGP,
     &                  MXPNSMST,
     &                  NSMOB,MXPTOB,NTPOB,NTSOB,NTESTG,MXPKA,
     &                  NEL1234,
     &                  MXCJ,MXCIJA,
     &                  MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,IPHGAS,
     &                  NHLFSPGP,MNHL,IADVICE)
*
* Find largest dimension of matrix C(Ka,Ib,J)
* Find largest dimension of matrix C(ij,Ka,Ib)
* Find largest dimension of matrix C(ij,Ia,Kb)
* Find largest dimension of matrix C(ij,Ka,Kb)
*
* Particle hole version : hole electrons added, particle elec removed
*
* Largest block of single excitations MXSXBL

*. Input
* IAB :allowed combination of alpha and beta supergroups
* IOCPTA : Number of first active alpha supergroup
* IOCPTB : Number of first active beta  supergroup
* NOCTPA : Number of active alpha supergroups
* NOCTPB : Number of active alpha supergroups
*
* Version of Jan 98 : IPHGAS added

      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IAB(NOCTPA,NOCTPB)
      DIMENSION NSTFSMSPGP(MXPNSMST,*)
      DIMENSION NTSOB(MXPTOB,NSMOB)
      DIMENSION NEL1234(MXPTOB,*)
      DIMENSION IPHGAS(*)
      INTEGER NHLFSPGP(*)
*
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.100) WRITE(6,*) ' MXRESC : MXPKA ', MXPKA
*
* matrix C(j,Ka,Ib)
*
*. Note : Only done for alpha-strings, problems when transposing
*         constructing C(J,Ia,Kb)
      MXCJ = 0
      MXADKBLK = 0
      DO IAORC= 1,2
      DO 100 IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO 200 IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXB = 0
            DO 210 ISM = 1, NSMST
              MXB =MAX(MXB,NSTFSMSPGP(ISM,IBTPABS))
  210       CONTINUE
            IF(NTEST.GE.100) WRITE(6,*) ' MXB = ', MXB
            DO 300 IOBTP = 1, NTPOB
*. No K strings obtained from creation in particle space
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1) GOTO 300
*. type of K string obtained
              CALL NEWTYP(IATPABS,IAORC,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
*. addi constraint to avoid calc with long columns and few rows
*. Works only in connection with active advice routine !
              IF(KATP.GT.0) THEN
                IF(IAORC.EQ.1.AND.IADVICE.EQ.1.AND.
     &          NHLFSPGP(IBTPABS)+NHLFSPGP(KATP).LT.MNHL.AND.
     &          NHLFSPGP(IATPABS).GT.(NHLFSPGP(IBTPABS)+1)) THEN
C                 WRITE(6,*) ' N-1 hole space eliminated '
C                 WRITE(6,*) ' IOBTP,IBTPABS,KATP',
C    &            IOBTP,IBTPABS,KATP
                  KATP = 0
                END IF
              END IF
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO 310 KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
  310           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                MXKAO = MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &          MXKA= MXPKA
                MXSOB = 0
                DO 320 ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
  320           CONTINUE
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKAO)
                LCJBLK = MXSOB*MXKA*MXB
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
  300       CONTINUE
          END IF
  200   CONTINUE
  100 CONTINUE
      END DO
*     ^ End of anni/crea map
*
* matrix C(j,Ia,Kb)
*
      DO IAORC = 1, 2
      DO IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA-1
        DO IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            IF(NTEST.GE.100)
     &      WRITE(6,*) ' allowed IATP,IBTP', IATP,IBTP
            MXA = 0
            DO ISM = 1, NSMST
              MXA =MAX(MXA,NSTFSMSPGP(ISM,IATPABS))
            END DO
            IF(NTEST.GE.100) WRITE(6,*) ' MXA = ', MXA
            DO IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1) GOTO 2812
              CALL NEWTYP(IBTPABS,IAORC,IOBTP,1,KBTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KBTP ',IOBTP,KBTP
              IF(KBTP.GT.0) THEN
                IF(IAORC.EQ.1.AND.IADVICE.EQ.1.AND.
     &          NHLFSPGP(IATPABS)+NHLFSPGP(KBTP).LT.MNHL.AND.
     &          NHLFSPGP(IBTPABS).GT.NHLFSPGP(IATPABS)+1) THEN
C                 WRITE(6,*) ' N-1 hole space eliminated '
C                 WRITE(6,*) ' IOBTP,IATPABS,KBTP',
C    &            IOBTP,IATPABS,KBTP
                  KBTP = 0
                END IF
              END IF
              IF(KBTP.GT.0) THEN
                MXKB = 0
                DO KSM = 1, NSMST
                  MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                MXKBO = MXKB
                IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &          MXKB= MXPKA
                MXSOB = 0
                DO ISMOB = 1, NSMOB
                  MXSOB = MAX(MXSOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXSOB = ', MXSOB
*
                MXADKBLK = MAX(MXADKBLK,MXSOB*MXKBO)
CE JULY-30-99   LCJBLK = MXSOB*MXKB*MXB
                LCJBLK = MXSOB*MXKB*MXA
                MXCJ = MAX(MXCJ,LCJBLK)
*
              END IF
 2812       CONTINUE
            END DO
          END IF
        END DO
      END DO
      END DO
*     ^ End of loop over creation/annihilation
      IF(NTEST.GT.100) THEN
        WRITE(6,*) 'MXRESC : MXADKBLK,MXCJ ', MXADKBLK,MXCJ
      END IF
*
* matrix C(ij,Ka,Ib)
* both Ka and Ib blocked
*
      MXCIJA = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA -1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1

          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIB = 0
            DO  ISM = 1, NSMST
              MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
            END DO
            IF(MXIB.GT.MXPKA) MXIB = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
            DO IAORC = 1, 2
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,IAORC,IOBTP,1,K1ATP)
*. No N+1 mappings for particle spaces
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1) K1ATP = 0
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1ATP ',IOBTP,K1ATP
              IF(K1ATP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JAORC = 1, 2
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1ATP,JAORC,JOBTP,1,KATP)
                  IF(JAORC.EQ.2.AND.IPHGAS(JOBTP).EQ.1) KATP = 0
                  IF(KATP.GT.0) THEN
                    MXKA = 0
                    DO KSM = 1, NSMST
                      MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                    IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA)
     &              MXKA= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKA*MXIB
                    MXCIJA = MAX(MXCIJA,LBLK)
                  END IF
                END DO
                END DO
*               ^ End of loop over JOBTP, JAORC
              END IF
            END DO
            END DO
*           ^ End of loop over IOBTP, IAORC
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) 'MXRESC : MXCIJA ', MXCIJA
      END IF
*
*
* matrix C(ij,Ia,kb)
* both Ka and Ib blocked
*
      MXCIJB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            MXIA = 0
            DO  ISM = 1, NSMST
              MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS ))
            END DO
            IF(MXIA.GT.MXPKA) MXIA = MXPKA
            IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
            DO IAORC = 1, 2
            DO  IOBTP = 1, NTPOB
*. type of K string obtained by removing one elec of type IOPBTP from IBTP
              CALL NEWTYP(IBTPABS,IAORC,IOBTP,1,K1BTP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP K1BTP ',IOBTP,K1BTP
              IF(IAORC.EQ.2.AND.IPHGAS(IOBTP).EQ.1)K1BTP = 0
              IF(K1BTP.GT.0) THEN
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JAORC = 1, 2
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from K1ATP
                  CALL NEWTYP(K1BTP,JAORC,JOBTP,1,KBTP)
                  IF(JAORC.EQ.2.AND.IPHGAS(JOBTP).EQ.1) KBTP = 0
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXIA
                    MXCIJB = MAX(MXCIJB,LBLK)
                  END IF
                END DO
                END DO
*               ^ End of loop over JOBTP,JAORC
              END IF
            END DO
            END DO
*           ^ End of loop over IOBTP,IAORC
          END IF
        END DO
      END DO
*
      IF(NTEST.GT.10) THEN
        WRITE(6,*) 'MXRESC : MXCIJB ', MXCIJB
      END IF
*
*
* matrix C(ij,Ka,kb)
* both Ka and Kb blocked
*
      MXCIJAB = 0
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        DO  IBTP = 1, NOCTPB
          IBTPABS = IBTP + IOCTPB - 1
          IF(IAB(IATP,IBTP).NE.0) THEN
            DO  IOBTP = 1, NTPOB
*. type of Ka string obtained by removing one elec of type IOPBTP from IATP
              CALL NEWTYP(IATPABS,1,IOBTP,1,KATP)
              IF(NTEST.GE.100)
     &        WRITE(6,*) ' IOBTP KATP ',IOBTP,KATP
              IF(KATP.GT.0) THEN
                MXKA = 0
                DO KSM = 1, NSMST
                  MXKA = MAX(MXKA,NSTFSMSPGP(KSM,KATP))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXKA = ',MXKA
                IF(MXPKA .GT. 0 .AND. MXKA .GT. MXPKA) MXKA= MXPKA
                MXISOB = 0
                DO ISMOB = 1, NSMOB
                  MXISOB = MAX(MXISOB,NTSOB(IOBTP,ISMOB))
                END DO
                IF(NTEST.GE.100) WRITE(6,*) ' MXISOB = ', MXISOB
                DO JOBTP = 1, NTPOB
*  type of K string obtained by removing one elec of type JOPBTP from IBTP
                  CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
                  IF(KBTP.GT.0) THEN
                    MXKB = 0
                    DO KSM = 1, NSMST
                      MXKB = MAX(MXKB,NSTFSMSPGP(KSM,KBTP))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXKB = ',MXKB
                    IF(MXPKA .GT. 0 .AND. MXKB .GT. MXPKA)
     &              MXKB= MXPKA
                    MXJSOB = 0
                    DO JSMOB = 1, NSMOB
                      MXJSOB = MAX(MXJSOB,NTSOB(JOBTP,JSMOB))
                    END DO
                    IF(NTEST.GE.100) WRITE(6,*) ' MXJSOB = ', MXJSOB
*
                    LBLK = MXISOB*MXJSOB*MXKB*MXKA
                    MXCIJAB = MAX(MXCIJAB,LBLK)
                  END IF
                END DO
              END IF
            END DO
          END IF
        END DO
      END DO
*
*
* Largest block of single excitations :
* Strings of given type and sym, orbitals of given type and sym
*
* Largest block of creations : a+i !kstring> where K string is
* obtained as single annihilations
      MXSXBL = 0
*. For alpha strings :
      DO  IATP = 1, NOCTPA
        IATPABS = IATP + IOCTPA - 1
        MXIA = 0
        DO  ISM = 1, NSMST
          MXIA = MAX(MXIA,NSTFSMSPGP(ISM,IATPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIA = ', MXIA
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP(IATPABS,1,JOBTP,1,KATP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KATP ',JOBTP,KATP
          IF(KATP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IATPABS))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP(KATP,2,IOBTP,1,JATP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
*
*. For beta  strings :
      DO  IBTP = 1, NOCTPB
        IBTPABS = IBTP + IOCTPB - 1
        MXIB = 0
        DO  ISM = 1, NSMST
          MXIB = MAX(MXIB,NSTFSMSPGP(ISM,IBTPABS))
        END DO
        IF(NTEST.GE.100) WRITE(6,*) ' MXIB = ', MXIB
*. Orbitals to be removed
        DO  JOBTP = 1, NTPOB
*. Is this removal allowed ??
          CALL NEWTYP(IBTPABS,1,JOBTP,1,KBTP)
          IF(NTEST.GE.100)
     &    WRITE(6,*) ' JOBTP KBTP ',JOBTP,KBTP
          IF(KBTP.GT.0) THEN
*. Number of possible choices of J orbitals
            MXJOB = 0
            DO JSMOB = 1, NSMOB
               MXJOB = MAX(MXJOB,NTSOB(JOBTP,JSMOB))
            END DO
            MXJOB = MIN(MXJOB,NEL1234(JOBTP,IBTP))
            IF(NTEST.GE.100) WRITE(6,*) ' MXJOB = ', MXJOB
*. Then  : add an electron
            DO IOBTP = 1, NTPOB
*  Allowed ?
              CALL NEWTYP(KBTP,2,IOBTP,1,JBTP)
              IF(JATP.GT.0) THEN
                MXIOB = 0
                DO ISMOB = 1, NSMOB
                  MXIOB = MAX(MXIOB,NTSOB(IOBTP,ISMOB))
                END DO
*
                MXSXBL = MAX(MXSXBL,MXIOB*MXJOB*MXIA)
              END IF
            END DO
          END IF
        END DO
      END DO
      IF(NTEST.GT.10) THEN
        WRITE(6,*) 'MXRESC: MXSXBL : ', MXSXBL
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE RASSG3(CB,SB,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &                  LUC,LUHC,SSCLFAC,ISCLFAC_GROUP
#ifdef VAR_MPI
     &                  ,LUCLIST,LUHCLIST,IBLOCKL,NPARBLOCK,RCCTOS,
     &                   IGROUPLIST,IPROCLIST
#endif
     &                 )
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen   Winter of 1991
*               May 1997 : Connected to SBLOCK
*
* =====
* Input
* =====
*
      IMPLICIT REAL*8(A-H,O-Z)
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      DIMENSION LUCLIST(*), IBLOCKL(*), NPARBLOCK(*), IGROUPLIST(*)
      INTEGER IPROCLIST(*), RCCTOS(*)
      DIMENSION ISCLFAC_GROUP(*), LUHCLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR1, IOFF_SCR2
      CHARACTER SECTID*12,WALLTID*12,WALLTID2*12
      LOGICAL WRITE_SB
#endif
#include "parluci.h"
*     batches of sigma
      INTEGER LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
*     scratch
      DIMENSION SB(*),CB(*),SSCLFAC(*)
      NZERO = 0
      IONE  = 1
      xreadtimebi = 0.0D0
      xcomputesi  = 0.0D0
      ISI_CALC_BL = 0
      IBI_MULT_BL = 0
*
      CALL QENTER('RASSG')

      NTEST     = 0000 ! normal print level
!     NTEST     = 9999 ! debug
      NPTESTVAR = NTEST

      IF(NTEST.GE.20) THEN
        WRITE(luwrt,*) ' ================='
        WRITE(luwrt,*) ' RASSG3 speaking :'
        WRITE(luwrt,*) ' ================='
        WRITE(luwrt,*) ' RASSG3 : NBATS = ',NBATS
      END IF
*
      IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
*
*     make an update of list ISCLFAC_GROUP using LUCLIST...
*
      BLOCKTIME = 0.0D0
*
*     only real part needed
*
      IRILP = 1
*
      CALL UPDATE_LUC_LIST2(ISCLFAC_GROUP,LUCLIST,RCCTOS,CB,
     &                      NPARBLOCK,IBLOCKL,IGROUPLIST,
     &                      IPROCLIST,IRILP,BLOCKTIME)
C
C     debug printing of 'complete' C vector on node/global file
C
      IF( NPTESTVAR .ge. 20 ) THEN
        IOFFSET_SCRATCH = 0
        WRITE(LUWRT,*) 'left vector = C-vector'
        DO IBLK = 1, NUM_BLOCKS2
C
           LENGTH_BLK = IBLOCKL(IBLK)
C
           IF( ISCLFAC_GROUP( IBLK ) .ne. 0 ) THEN
             CALL MPI_FILE_READ_AT(LUC,IOFFSET_SCRATCH,CB,LENGTH_BLK,
     &                             MPI_REAL8,my_STATUS,IERR)
             WRITE(LUWRT,'(A,1X,I12,1X,A,1X,I6,1X,I6)') 'Read-in at',
     &             IOFFSET_SCRATCH,'for block', IBLK, LENGTH_BLK
             CALL WRTMATMN(CB,1,LENGTH_BLK,1,LENGTH_BLK,LUWRT)
           END IF
           IOFFSET_SCRATCH = IOFFSET_SCRATCH + LENGTH_BLK
          END DO
      END IF
      IOFFSET_SCRATCH = 0
C
C     initialize writing of vector file
C
      IOFF_SCR1 = 0
      IOFF_SCR2 = 0
C     file offset
      IOFF_SCR1 = MY_LU4_OFF
     &          + MY_VEC2_IOFF * ( JVEC_SF )
      INT_IOFF1  = 0
      INT_IOFF2  = 0
C     block list offset
      INT_IOFF1  = 1 + MY_ACT_BLK2 * ( JVEC_SF )
C
C?    WRITE(LUCIWT,*) 'INT_IOFF1 in RASSG3',INT_IOFF1
C
C     include in C-vector batches only C-blocks which are 
C     needed to compute a s-block for a given CPU
C
      CALL UPDATE_GEN_LIST(ISCLFAC_GROUP,RCCTOS,NUM_BLOCKS2)
      CALL IFACTOSFAC(ISCLFAC_GROUP,SSCLFAC,NUM_BLOCKS2)
C
C     debug info
C
      IF( NPTESTVAR .ge. 20 ) THEN
         WRITE(LUWRT,*) ' ISCLFAC_GROUP '
         CALL IWRTMAMN(ISCLFAC_GROUP,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUWRT)
         WRITE(LUWRT,*) ' RCCTOS '
         CALL IWRTMAMN(RCCTOS,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUWRT)
         WRITE(LUWRT,*) ' SSCLFAC '
         CALL WRTMATMN(SSCLFAC,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUWRT)
      END IF
C
*=====================================================
#endif
      ELSE ! (LUCI_NMPROC .GT. 1) THEN
         IF(LUHC.GT.0) CALL REWINE(LUHC,-1)
C
C        get scaling factor - 0.0D0 or 1.0D0
C
         CALL FIND_ACTIVE_BLOCKS(LUC,-1,SSCLFAC,CB)
      END IF ! (LUCI_NMPROC .GT. 1) THEN
C
C     loop over batches over sigma blocks
C
      DO JBATS = 1, NBATS
C
#ifdef VAR_MPI
       IF (LUCI_NMPROC .GT. 1) THEN
*
*       start various timings...
*
        starttime  = mpi_wtime()
        sbatchtime = 0.0D0
        writetime  = 0.0D0
C
C       do we need to calculate a single s-block?
C
        ICOMPUTE = 0
*
        DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
           IF( NPARBLOCK( ISBLK ) .eq. LUCI_MYPROC ) THEN
             ICOMPUTE = 1
           END IF
        END DO
csk     WRITE(LUWRT,*) ' ICOMPUTE is',ICOMPUTE
        IF( ICOMPUTE .eq. 0 ) THEN
           starttimer = MPI_WTIME()
           GOTO 60
        ENDIF
       END IF ! (LUCI_NMPROC .GT. 1) THEN
*
#endif /* if def VAR_MPI */
C       zero SB
        LS = LEBATS(JBATS)
!       WRITE(LUWRT,*) ' LS - sbtach length =',LS
        CALL DZERO(SB,LS)
C
        CALL SBLOCK(LBATS(JBATS),IBATS(1,I1BATS(JBATS)),1,
     &       CB,SB,LUC,0,0,0,0,0,SSCLFAC,ISCLFAC_GROUP
#ifdef VAR_MPI
     &      ,LUCLIST,NPARBLOCK,IBLOCKL
#endif
     &                             )
C  
C       transfer S block to permanent storage
C
 60     CONTINUE
C
#ifdef VAR_MPI
      IF (LUCI_NMPROC .GT. 1) THEN
        starttimer = MPI_WTIME()
        IF( ICOMPUTE .eq. 0 ) GOTO 70
      END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif /* if def VAR_MPI */
C
        if (NTEST.ge.100) then
            WRITE(luwrt,*) ' **********************************'
            WRITE(luwrt,*) ' *                                *'
            WRITE(luwrt,*) ' * Array containing final sbatch  *'
            WRITE(luwrt,*) ' *                                *'
            WRITE(luwrt,*) ' **********************************'
            CALL WRTMATMN(SB,1,LS,1,LS,LUCIWT)
        end if
C
        DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
           IOFF = IBATS(6,ISBLK)
           LEN  = IBATS(8,ISBLK)
         IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
*          check again if block is needed
           IF( NPARBLOCK(ISBLK) .ne. LUCI_MYPROC ) GOTO 65
*
*          new offset
*
           IOFF_SCR1 = IOFF_SCR1 + IOFF_SCR2
           INT_IOFF1 = INT_IOFF1 + INT_IOFF2
*          check norm of the block
           XXX = 0.0D0
           XXX = DDOT(LEN,SB(IOFF),1,SB(IOFF),1)
           WRITE_SB = .FALSE.
           IF( XXX .eq. 0.0D0 ) THEN
              LUHCLIST(INT_IOFF1) = NZERO
           ELSE
*             mark block as nonzero
              LUHCLIST(INT_IOFF1) = IONE
              WRITE_SB = .TRUE.
           END IF
*          transfer block to disc
           IF( WRITE_SB ) THEN
#ifdef LUCI_DEBUG
              WRITE(luwrt,*) 'I transfer to file block ISBLK',
     &                        ISBLK
              CALL WRTMATMN(SB(IOFF),1,LEN,1,LEN,luwrt)
              WRITE(luwrt,*) 'offset,length',IOFF_SCR1,LEN
#endif
              CALL MPI_FILE_WRITE_AT(LUHC,IOFF_SCR1,SB(IOFF),LEN,
     &                               MPI_REAL8,my_STATUS,IERR)
              ISI_CALC_BL = ISI_CALC_BL + 1
           END IF
*          keep track of correct offset
           IOFF_SCR2 = LEN
           INT_IOFF2 = 1
 65        CONTINUE
C
#endif
          ELSE ! (LUCI_NMPROC .GT. 1) THEN
           CALL ITODS(LEN,1,-1,LUHC)
           CALL TODSC_LUCI(SB(IOFF),LEN,-1,LUHC)
          END IF ! (LUCI_NMPROC .GT. 1) THEN
        END DO
 70     CONTINUE
#ifdef VAR_MPI
        IF (LUCI_NMPROC .GT. 1) THEN
          IF( ICOMPUTE .ne. 0 ) THEN
            writetime = writetime + MPI_WTIME() - starttimer
            sbatchtime = sbatchtime + MPI_WTIME() - starttime
CSK too much information?                WALLTID  = SECTID(writetime)
            WALLTID2 = SECTID(sbatchtime)
CSK too much information?                 WRITE(LUWRT,50000) WALLTID
            IF( TIMING_par ) THEN
              WRITE(LUWRT,60001) WALLTID2
            END IF
CSK         WRITE(LUWRT,60000) WALLTID2, JBATABS
          END IF
        END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif /* if def VAR_MPI */
      END DO
*
C
      IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
      WALLTID2 = SECTID(BLOCKTIME)
      IF( TIMING_par ) THEN
        WRITE(LUWRT,80000) WALLTID2
      END IF
*
50000 FORMAT(' >>>  WALL TIME FOR WRITING SB TO DISK           :
     & ',A)
60001 FORMAT(' >>>  WALL TIME IN TOTAL FOR CURRENT S-BATCH     :
     & ',A)
80000 FORMAT(' >>>  WALL TIME FOR C-COEFFICIENT EXCHANGE       :
     & ',A)
*
C
      IF( TIMING_par ) THEN
C
C       print statistics
C
        WRITE(LUWRT,'(/A)')
     &  '               H x b_i contraction run statistics  '
        WRITE(LUWRT,'(A/)')
     &  '              ____________________________________ '
        WRITE(LUWRT,'(2X,A,1X,I9)')
     &  ' number of s_i blocks calculated          : ',ISI_CALC_BL
        WRITE(LUWRT,'(2X,A,1X,I9)')
     &  ' total number of b_i read from disk       : ',IBI_MULT_BL
        WALLTID = SECTID(xreadtimebi)
        WRITE(LUWRT,'(2X,A,1X,A)')
     &  ' read time for b_i blocks from disk       : ', WALLTID
        WALLTID = SECTID(xcomputesi)
        WRITE(LUWRT,'(2X,A,1X,A/)')
     &  ' pure matrix x vector multiplication time : ', WALLTID
C
      END IF
C
#endif
      ELSE ! (LUCI_NMPROC .GT. 1) THEN
C     write EOF mark
      CALL ITODS(-1,1,-1,LUHC)
C
      IF(NTEST.GE.100) THEN
        WRITE(LUCIWT,*) ' Final S-vector on disc'
        CALL WRTVCD(SB,LUHC,1,-1)
      END IF
      END IF ! (LUCI_NMPROC .GT. 1) THEN
*
      CALL QEXIT('RASSG')
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE RSBB1E(ISCSM,ISCTP,ISBCTP,ICCSM,ICCTP,ICBCTP,
     &                  IGRP,NROW,
     &                  NGAS,ISEL,ICEL,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,
     &                  MXPNGAS,NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,MOC,
     &                  NSCOL,SCLFAC,
     &                  IUSE_PH,IPHGAS,NTESTG)
*
* One electron excitations on column strings
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ISBCTP : Base for sigma column types
* ICCSM,ICCTP : Symmetry and type of C     columns
* ICBCTP : Base for C     column types
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active sets
* ISEL : Occupation in each active set for sigma block
* ICEL : Occupation in each active set for C     block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* MOC  : Use MOC method ( instead of N-1 resolution method )
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for one electron integrals
*
* Jeppe Olsen, Winter of 1991
*              IUSE_PH added winter of 97
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
      REAL*8 INPROD
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER IPHGAS(NGAS)
*.Specific Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),H(*)
      DIMENSION I2(*),XI2S(*)
*.Local arrays ( assume MPNGAS = 16 ) !!!
      DIMENSION ITP(16),JTP(16)
      DIMENSION ISGRP(16),ICGRP(16)
*. For transposing integral block
      PARAMETER(MXPTSOB = 35)
      DIMENSION HSCR(MXPTSOB*MXPTSOB)
*
      DIMENSION IJ_REO(2),IJ_DIM(2),IJ_SM(2),IJ_TP(2),IJ_AC(2)
      DIMENSION ISCR(2)
      CALL QENTER('RS1   ')
* Type of single excitations that connects the two column strings
C     MOC = 1
      NTESTL = 000
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500)THEN
        WRITE(lupri,*)
        WRITE(lupri,*) ' ======================= '
        WRITE(lupri,*) ' Information from RSBB1E '
        WRITE(lupri,*) ' ======================= '
        WRITE(lupri,*)
        WRITE(lupri,*) ' RSBB1E : MOC,IUSE_PH ', MOC,IUSE_PH
        WRITE(lupri,*) ' ISEL : '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(lupri,*) ' ICEL : '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*. Obtain groups
C     GET_SPGP_INF(ISPGP,ITP,IGRP)
      CALL GET_SPGP_INF(ICCTP,IGRP,ICGRP)
      CALL GET_SPGP_INF(ISCTP,IGRP,ISGRP)
*
      IFRST = 1
      JFRST = 1
*. Types of single excitations that connect ISEL and ICEL
      CALL SXTYP2_GAS(NSXTP,ITP,JTP,NGAS,ISEL,ICEL,IPHGAS)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        IF(NTEST.GE.2000)
     &  write(lupri,*) ' ITYP JTYP ', ITYP,JTYP
*. Is this combination of types allowed
        IJ_ACT = I_SX_ACT(ITYP,JTYP)
        IF(IJ_ACT.EQ.0) GOTO 900
*. Hvilken vej skal vi valge,
        NOP = 2
        IJ_AC(1) = 2
        IJ_AC(2) = 1
        IJ_TP(1) = ITYP
        IJ_TP(2) = JTYP
        IF(IUSE_PH.EQ.1) THEN
          CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TP,IJ_AC,IJ_REO,SIGNIJ)
        ELSE
          IJ_REO(1) = 1
          IJ_REO(2) = 2
          SIGNIJ = 1.0D0
        END IF
*
        ISCR(1) = IJ_AC(1)
        ISCR(2) = IJ_AC(2)
        IJ_AC(1) = ISCR(IJ_REO(1))
        IJ_AC(2) = ISCR(IJ_REO(2))
*
        ISCR(1) = ITYP
        ISCR(2) = JTYP
        IJ_TP(1) = ISCR(IJ_REO(1))
        IJ_TP(2) = ISCR(IJ_REO(2))
*
        DO 800 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
*. New intermediate strings will be accessed so
          KFRST = 1
          IF(JSM.EQ.0) GOTO 800
          IF(NTEST.GE.2000)
     &    write(lupri,*) ' ISM JSM ', ISM,JSM
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
*. Reorder
*
          ISCR(1) = ISM
          ISCR(2) = JSM
          IJ_SM(1) = ISCR(IJ_REO(1))
          IJ_SM(2) = ISCR(IJ_REO(2))
*
          ISCR(1) = NIORB
          ISCR(2) = NJORB
          IJ_DIM(1) = ISCR(IJ_REO(1))
          IJ_DIM(2) = ISCR(IJ_REO(2))
*
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*. Fetch integrals : For CI-transformations using RSBB1E
*. most of the blocks vanishes
*.Obtain one electron integrals (ISM,ITP,JSM,JTP) transposed
           IF(IJ_REO(1).EQ.1) THEN
*. obtain integrals h(j,i)
             CALL GETH1(HSCR,IJ_SM(1),IJ_TP(1),IJ_SM(2),IJ_TP(2))
             CALL TRPMAT(HSCR,IJ_DIM(1),IJ_DIM(2),H)
           ELSE
*. Obtain integrals h(i,j)
             CALL GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
           END IF
           XNORM = INPROD(H,H,IJ_DIM(1)*IJ_DIM(2))
           IF(XNORM.EQ.0) GOTO 800
          IF(MOC.EQ.0) THEN
*
*
* ======================================================================
*.                   Use N-1 resolution method
* ======================================================================
*
*
*. Obtain annihilation/creation maps for all K strings
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          SCLFACS = SIGNIJ*SCLFAC
          IF(NTEST.GE.1000)
     &    WRITE(lupri,*) ' IJ_SM,IJ_TP,IJ_AC',IJ_SM(2),IJ_TP(2),IJ_AC(2)
          CALL ADAST_GAS(IJ_SM(2),IJ_TP(2),NGAS,ICGRP,ICCSM,
     &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFACS,IJ_AC(1))
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          ONE = 1.0D0
          CALL ADAST_GAS(IJ_SM(1),IJ_TP(1),NGAS,ISGRP,ISCSM,
     &         I2,XI2S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJ_AC(1))
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I2,XI2S,IJ_DIM(1),
     &             NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF
*. Loop over partitionings of the row strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = MIN(KTOP + MAXK,NKAEFF)
              IF(KTOP.EQ.NKAEFF) THEN
                KEND = 1
              ELSE
                KEND = 0
              END IF
              LKABTC = KTOP - KBOT +1
*. This is the place to start over partitioning of I strings
              DO 701 IPART = 1, NIPART
                IBOT = (IPART-1)*MAXI+1
                ITOP = MIN(IBOT+MAXI-1,NROW)
                NIBTC = ITOP - IBOT + 1
* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
                DO JJORB = 1,IJ_DIM(2)
                  ICGOFF = 1 + (JJORB-1)*LKABTC*NIBTC
                  CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                 LKABTC,I1(KBOT+(JJORB-1)*NKASTR),
     &                 XI1S(KBOT+(JJORB-1)*NKASTR) )
                END DO
*.Obtain one electron integrals (ISM,ITP,JSM,JTP) transposed
C               CALL GETH1(HSCR,IJ_SM(1),IJ_TP(1),IJ_SM(2),IJ_TP(2))
C               CALL TRPMAT(HSCR,IJ_DIM(1),IJ_DIM(2),H)
*. Problems when HOLE switches blocks around ?
C               CALL GETH1(H,IJ_SM(2),IJ_TP(2),IJ_SM(1),IJ_TP(1))
                IF(NTEST.GE.1000) THEN
                  WRITE(lupri,*) ' RSBB1E H BLOCK '
                CALL WRTMT_LU(H,IJ_DIM(2),IJ_DIM(1),IJ_DIM(2),IJ_DIM(1))
                END IF
*.Sscr(I,K,i) = CSCR(I,K,j)*h(j,i)
                NIK = NIBTC*LKABTC
                FACTORC = 0.0D0
                FACTORAB = 1.0D0
                IF(NTEST.GE.2000) THEN
                  WRITE(lupri,*) ' CSCR array,NIK X NJORB array '
                  CALL WRTMT_LU(CSCR,NIK,IJ_DIM(2),NIK,IJ_DIM(2))
                END IF
                CALL MATML7(SSCR,CSCR,H,
     &               NIK,IJ_DIM(1),NIK,IJ_DIM(2),IJ_DIM(2),IJ_DIM(1),
     &               FACTORC,FACTORAB,0)
                IF(NTEST.GE.2000) THEN
                  WRITE(lupri,*) ' SSCR array,NIK X NIORB array '
                  CALL WRTMT_LU(SSCR,NIK,IJ_DIM(1),NIK,IJ_DIM(1))
                END IF
*.S(I,a+ K) =  S(I, a+ K) + sgn*Sscr(I,K,i)
                DO IIORB = 1,IJ_DIM(1)
                  ISBOFF = 1+(IIORB-1)*LKABTC*NIBTC
                  CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                 LKABTC,I2(KBOT+(IIORB-1)*NKASTR),
     &                 XI2S(KBOT+(IIORB-1)*NKASTR))
                END DO
*
  701       CONTINUE
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
*. End of loop over I partitioninigs
          END IF
*.(End of algorithm switch)
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
      CALL QEXIT('RS1  ')
C!    WRITE(6,*) ' Enforced stop in RSBB1E '
C!    Call Abend2( ' Enforced stop in RSBB1E ' )
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE RSBB2A(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NSCOL,
     &                  NGAS,ISOC,ICOC,
     &                  SB,CB,
     &                  ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &                  RIKSX,RKJSX,MXSXBL,IMOC,SCLFAC,NTESTG,
     &                  NSEL2E,ISEL2E,IUSE_PH,IPHGAS,MXTSOB,SCR)
*
* two electron excitations on column strings
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NSCOL : Number of columns in S block
* ISEL1(3) : Number of electrons in RAS1(3) for S block
* ICEL1(3) : Number of electrons in RAS1(3) for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* DXSTST : Sym of sx,!st> => sym of sx !st>
* STSTDX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* SXDXSX : Symmetry of SX1,SX1*SX2 => symmetry of SX2
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals, strings,
*                           single excitations, double excitations
* MAXI   : Largest number of 'spectator strings' treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
*
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),DXSTST(NSMDX,NSMST)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER IPHGAS(NGAS)
*
      INTEGER ISEL2E(*)
*.Input
      DIMENSION CB(*)
      INTEGER ISOC(NGAS),ICOC(NGAS)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),XINT(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
      DIMENSION RIKSX(MXSXBL,4),RKJSX(MXSXBL,4)
      dimension SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
*.Local arrays
      DIMENSION ITP(256),JTP(256),KTP(256),LTP(256)
      INTEGER I4_DIM(4),I4_SM(4),I4_TP(4),I4_REO(4),ISCR(4)
      INTEGER I4_AC(4)
*
      INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
CTF   PARAMETER(MXTSOB=40)
CTF   COMMON/SOMESCR/SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
*
#include "comjep.inc"
*
      DIMENSION IACAR(2),ITPAR(2)
      CALL QENTER('RS2A')
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.1000) THEN
        WRITE(lupri,*) ' ================'
        WRITE(lupri,*) ' RSBB2A speaking '
        WRITE(lupri,*) ' ================'
        WRITE(lupri,*) ' MXSXBL = ', MXSXBL
        WRITE(lupri,*) ' RSBB2A : IMOC,IUSE_PH ', IMOC, IUSE_PH
        WRITE(lupri,*) ' ISOC and ICOC : '
        CALL IWRTMA(ISOC,1,NGAS,1,NGAS)
        CALL IWRTMA(ICOC,1,NGAS,1,NGAS)
      END IF
      IFRST = 1
      JFRST = 1
*
*.Types of DX that connects the two strings
*
      IDXSM = STSTDX(ISCSM,ICCSM)
      IF(IDXSM.EQ.0)  GOTO 2001
*. Connecting double excitations
      CALL DXTYP2_GAS(NDXTYP,ITP,JTP,KTP,LTP,NGAS,ISOC,ICOC,IPHGAS)
      DO 2000 IDXTYP = 1, NDXTYP
        ITYP = ITP(IDXTYP)
        JTYP = JTP(IDXTYP)
        KTYP = KTP(IDXTYP)
        LTYP = LTP(IDXTYP)
*. Is this combination of types allowed
         IJKL_ACT = I_DX_ACT(ITYP,KTYP,LTYP,JTYP)
         IF(IJKL_ACT.EQ.0) GOTO 2000
*
C?      write(6,*) ' test inserted in RSBB2AN'
C?      NPTOT = 0
C?      IF(ITYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(JTYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(KTYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(LTYP.EQ.3) NPTOT = NPTOT + 1
C?      IF(NPTOT.EQ.3) GOTO 2000
*
        ITYP_ORIG = ITYP
        JTYP_ORIG = JTYP
        KTYP_ORIG = KTYP
        LTYP_ORIG = LTYP
*
        NIJKL1 = 0
        IF(ITYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(JTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(KTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(LTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(NIJKL1.EQ.0) CALL QENTER('BB2A0')
        IF(NIJKL1.EQ.1) CALL QENTER('BB2A1')
        IF(NIJKL1.EQ.2) CALL QENTER('BB2A2')
        IF(NIJKL1.EQ.3) CALL QENTER('BB2A3')
        IF(NIJKL1.EQ.4) CALL QENTER('BB2A4')
*. Optimal ordering of operators
        I4_AC(1) = 2
        I4_AC(2) = 2
        I4_AC(3) = 1
        I4_AC(4) = 1
        I4_TP(1) = ITYP
        I4_TP(2) = KTYP
        I4_TP(3) = LTYP
        I4_TP(4) = JTYP
        IF(IUSE_PH.EQ.1) THEN
          NOP = 4
          CALL ALG_ROUTERX(ISOC,JSOC,NOP,I4_TP,I4_AC,I4_REO,SIGN4)
        ELSE
          DO IJKL = 1, 4
            I4_REO(IJKL) = IJKL
          END DO
          SIGN4 = 1.0D0
        END IF
*. Type of operators : TP and AC
        DO IJKL = 1, 4
C         ISCR( I4_REO(IJKL) ) = I4_TP(IJKL)
          ISCR(IJKL) = I4_TP( I4_REO(IJKL) )
        END DO
        DO IJKL = 1, 4
          I4_TP(IJKL) = ISCR(IJKL)
        END DO
        ITYP = I4_TP(1)
        KTYP = I4_TP(2)
        LTYP = I4_TP(3)
        JTYP = I4_TP(4)
        DO IJKL = 1, 4
          ISCR(IJKL) = I4_AC( I4_REO(IJKL) )
        END DO
        DO IJKL = 1, 4
          I4_AC(IJKL) = ISCR(IJKL)
        END DO
        IF(NTEST.GE.100) THEN
          WRITE(lupri,*) ' I4_AC, IT_TP  defined '
          WRITE(lupri,*) ' I4_AC, I4_TP '
          CALL IWRTMA(I4_AC,1,4,1,4)
          CALL IWRTMA(I4_TP,1,4,1,4)
        END IF
*
*      ==================================
        IF(I4_AC(1).EQ.I4_AC(2) ) THEN
*      ==================================
*
*. a+ a+ a a or a a a+ a+
*. Integrals included ??
          IF(NSEL2E.NE.0) THEN
            IAMOKAY=0
            IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP) THEN
              DO JSEL2E = 1, NSEL2E
                IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY=1
              END DO
            END IF
            IF(IAMOKAY.EQ.0) GOTO 2000
          END IF
*. Largest possible number of orbital pairs
          MI = 0
          MJ = 0
          MK = 0
          ML = 0
          DO IOBSM = 1, NSMST
            MI = MAX(MI,NOBPTS(ITYP,IOBSM))
            MJ = MAX(MJ,NOBPTS(JTYP,IOBSM))
            MK = MAX(MK,NOBPTS(KTYP,IOBSM))
            ML = MAX(ML,NOBPTS(LTYP,IOBSM))
          END DO
          MXPAIR = MAX(MI*MK,MJ*ML)
*. Largest posssible
*. Symmetry of allowed Double excitation,loop over excitations
          DO 1950 IKOBSM = 1, NSMOB
            JLOBSM = SXDXSX(IKOBSM,IDXSM)
            IF(NTEST.GE.100)
     &      WRITE(lupri,*) ' IKOBSM,JLOBSM', IKOBSM,JLOBSM
            IF(JLOBSM.EQ.0) GOTO 1950
*. types + symmetries defined => K strings are defined
            KFRST = 1
*
*. Number of batchs of symmetry pairs IK
*
            LENGTH = 0
            NIKBT = 0
            NBLK = 0
            NBLKT = 0
            DO ISM = 1, NSMOB
              KSM = ADSXA(ISM,IKOBSM)
              NI = NOBPTS(ITYP,ISM)
              NK = NOBPTS(KTYP,KSM)
              IF(NTEST.GE.100) write(lupri,*) ' NI, NK' , NI,NK
*
              IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                NIK = NI*(NI+1)/2
              ELSE IF(ITYP.GT.KTYP.OR.(ITYP.EQ.KTYP.AND.ISM.GT.KSM))THEN
                NIK = NI*NK
              ELSE
                NIK = 0
              END IF
              IF(NIK.NE.0) THEN
                NBLKT = NBLKT + 1
                IF(LENGTH+NIK .GT. MXPAIR) THEN
*. The present batch is complete
                  NIKBT = NIKBT+1
                  IKBT(1,NIKBT) = NBLKT - NBLK
                  IKBT(2,NIKBT) = NBLK
                  IKBT(3,NIKBT) = LENGTH
                  LENGTH = 0
                  NBLK = 0
                END IF
                NBLK = NBLK + 1
                LENGTH = LENGTH + NIK
                IKSMBT(1,NBLKT) = ISM
                IKSMBT(2,NBLKT) = KSM
              END IF
            END DO
*. The last batch
            IF(NBLK.NE.0) THEN
              NIKBT = NIKBT+1
              IKBT(1,NIKBT) = NBLKT - NBLK + 1
              IKBT(2,NIKBT) = NBLK
              IKBT(3,NIKBT) = LENGTH
            END IF

*.
            IF(NTEST.GE.2000) THEN
              WRITE(lupri,*) ' ITYP, KTYP, IKOBSM,  NIKBT = ',
     &                     ITYP, KTYP, IKOBSM,  NIKBT
              WRITE(lupri,*) ' IKBT : Offset, number, length '
              DO JIKBT = 1, NIKBT
                WRITE(lupri,'(3i3)') (IKBT(II,JIKBT), II = 1, 3)
              END DO
              WRITE(lupri,*) ' IKSMBT '
              CALL IWRTMA(IKSMBT,2,NBLKT,2,8)
            END IF
*
*. Number of batchs of symmetry pairs JL
*
            LENGTH = 0
            NJLBT = 0
            NBLK = 0
            NBLKT = 0
            DO JSM = 1, NSMOB
              LSM = ADSXA(JSM,JLOBSM)
              NJ = NOBPTS(JTYP,JSM)
              NL = NOBPTS(LTYP,LSM)
*
              IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                NJL = NJ*(NJ+1)/2
              ELSE IF(JTYP.GT.LTYP.OR.(JTYP.EQ.LTYP.AND.JSM.GT.LSM))THEN
                NJL = NJ*NL
              ELSE
                NJL = 0
              END IF
              IF(NJL.NE.0) THEN
                NBLKT = NBLKT + 1
                IF(LENGTH+NJL .GT. MXPAIR) THEN
*. The present batch is complete
                  NJLBT = NJLBT+1
                  JLBT(1,NJLBT) = NBLKT - NBLK
                  JLBT(2,NJLBT) = NBLK
                  JLBT(3,NJLBT) = LENGTH
                  LENGTH = 0
                  NBLK = 0
                END IF
                NBLK = NBLK + 1
                LENGTH = LENGTH + NJL
                JLSMBT(1,NBLKT) = JSM
                JLSMBT(2,NBLKT) = LSM
              END IF
            END DO
*. The last batch
            IF(NBLK.NE.0) THEN
              NJLBT = NJLBT+1
              JLBT(1,NJLBT) = NBLKT - NBLK + 1
              JLBT(2,NJLBT) = NBLK
              JLBT(3,NJLBT) = LENGTH
            END IF
*.
            IF(NTEST.GE.2000) THEN
              WRITE(lupri,*) ' JTYP, LTYP, JLOBSM,  NJLBT = ',
     &                     JTYP, LTYP, JLOBSM,  NJLBT
              WRITE(lupri,*) ' JLBT : Offset, number, length '
              DO JJLBT = 1, NJLBT
                WRITE(lupri,'(3i3)') (JLBT(II,JJLBT), II = 1, 3)
              END DO
              WRITE(lupri,*) ' JLSMBT '
              CALL IWRTMA(JLSMBT,2,NBLKT,2,8)
            END IF
*
*. Loop over batches of IK strings
            DO 1940 IKBTC = 1, NIKBT
              IF(NTEST.GE.1000) WRITE(lupri,*) ' IKBTC = ', IKBTC
*. Loop over batches of JL strings
              DO 1930 JLBTC = 1, NJLBT
                IFIRST = 1
*. Loop over batches of I strings
                NPART = NROW/MAXI
                IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
                IF(NTEST.GE.2000)
     &          write(lupri,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
                DO 1801 IPART = 1, NPART
                  IBOT = 1+(IPART-1)*MAXI
                  ITOP = MIN(IBOT+MAXI-1,NROW)
                  NIBTC = ITOP-IBOT+1
*.Loop over batches of intermediate strings
                  KBOT = 1- MAXK
                  KTOP = 0
 1800             CONTINUE
                    KBOT = KBOT + MAXK
                    KTOP = KTOP + MAXK
*
                    IONE = 1
                    JLBOFF = 1
                    NJLT = JLBT(3,JLBTC)
                    DO JLPAIR = 1, JLBT(2,JLBTC)
                      JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                      LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
                      NJ = NOBPTS(JTYP,JSM)
                      NL = NOBPTS(LTYP,LSM)
                      IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                        NJL = NJ*(NJ+1)/2
                        JLSM = 1
                      ELSE
                        NJL = NJ * NL
                        JLSM = 0
                      END IF

*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Obtain all double excitations from this group of K strings
CT                    CALL QENTER('ADADS')
                      I12 = 1
                      K12 = 1
                      IONE = 1
C?       write(6,*) ' Before ADAADAST '
*. Creation / annihilation maps , conjugated of above
                      IF(I4_AC(4).EQ.1) THEN
                        JAC = 2
                      ELSE
                        JAC = 1
                      END IF
                      IF(I4_AC(3).EQ.1) THEN
                        LAC = 2
                      ELSE
                        LAC = 1
                      END IF
                      CALL ADAADAST_GAS(IONE,JSM,JTYP,NJ,JAC,
     &                                  IONE,LSM,LTYP,NL,LAC,
     &                            ICCTP,ICCSM,IGRP,
     &                            KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                            JFRST,KFRST,I12,K12,SCLFAC)
                      JFRST = 0
                      KFRST = 0
*
CT                    CALL QEXIT('ADADS')
                      IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                      J = 0
                      L = 1
*
CT                    CALL QENTER('MATCG')
                      DO  IJL = 1, NJL
                        CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                        I1JL = (L-1)*NJ+J
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                        JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                        IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                          ZERO = 0.0D0
                          CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                        ELSE
                          CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,
     &                              NKBTC,I1(1,I1JL),XI1S(1,I1JL))
                        END IF
                      END DO
CT                    CALL QEXIT ('MATCG')
*
                      JLBOFF = JLBOFF + NJL
                    END DO
*
*. ( End of loop over jlpair in batch )
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
*.Obtain two electron integrals xint(ik,jl) = (ij!kl)-(il!kj)
                    IF(IFIRST.EQ.1) THEN
                      IXCHNG = 1
* Obtain integrals in ik batch
                      NIKT = IKBT(3,IKBTC)
                      NJLT = JLBT(3,JLBTC)
                      JLOFF = 1
                      DO JLPAIR = 1, JLBT(2,JLBTC)
                      IKOFF = 1
                      DO IKPAIR = 1, IKBT(2,IKBTC)
*
                        ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                        KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                        JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                        LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
*
                        IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                          IKSM = 1
                          NIK =
     &                    NOBPTS(ITYP,ISM)*(NOBPTS(ITYP,ISM)+1)/2
                        ELSE
                          IKSM = 0
                          NIK =
     &                    NOBPTS(ITYP,ISM)*NOBPTS(KTYP,KSM)
                        END IF
*
                        IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                          JLSM = 1
                          NJL =
     &                    NOBPTS(JTYP,JSM)*(NOBPTS(JTYP,JSM)+1)/2
                        ELSE
                          JLSM = 0
                          NJL =
     &                    NOBPTS(JTYP,JSM)*NOBPTS(LTYP,LSM)
                        END IF
* ================================================================
*. Required form of integrals : Coulomb - Exchange of just Coulomb
* ================================================================
                        ICOUL = 0
*. Use coulomb - exchange
                        IXCHNG = 1
*. fetch integrals
                        CALL LGETINT(SCR,ITYP,ISM,JTYP,JSM,KTYP,KSM,
     &                              LTYP,LSM,IXCHNG,IKSM,JLSM,ICOUL)
                        DO JL = 1, NJL
                          CALL COPVEC(SCR((JL-1)*NIK+1),
     &                         XINT((JLOFF-1+JL-1)*NIKT+IKOFF),NIK)
                        END DO
                        IKOFF = IKOFF + NIK
                      END DO
                      JLOFF = JLOFF + NJL
                      END DO
                    END IF
*                   ^ End if integrals should be fetched
                    IFIRST = 0
*.and now ,to the work
                    LIKB = NIBTC*NKBTC
                    IF(NTEST.GE.3000) THEN
                     WRITE(lupri,*) ' Integral block '
                     CALL WRTMT_LU(XINT,NIKT,NJLT,NIKT,NJLT)
                    END IF
                    IF(NTEST.GE.3000) THEN
                      WRITE(lupri,*) ' CSCR matrix '
                      CALL WRTMT_LU(CSCR,LIKB,NJLT,LIKB,NJLT)
                    END IF
*
C?                  MXACIJO = MXACIJ
                    MXACIJ = MAX(MXACIJ,LIKB*NJLT,LIKB*NIKT)
C?                  IF(MXACIJ.GT.MXACIJO) THEN
C?                    write(6,*) ' New max MXACIJ = ', MXACIJ
C?                    write(6,*) ' ISCTP,ICCTP', ISCTP,ICCTP
C?                    WRITE(6,*) ' ITYP,JTYP,KTYP,LTYP',
C?   &                             ITYP,JTYP,KTYP,LTYP
C?                    WRITE(6,*)'NIJT, NJLT, NIBTC NKBTC',
C?   &                           NIJT, NJLT,NIBTC,NKBTC
C?                  END IF
*
                    FACTORC = 0.0D0
                    FACTORAB = 1.0D0
                    CALL MATML7(SSCR,CSCR,XINT,
     &                          LIKB,NIKT,LIKB,NJLT,NIKT,NJLT,
     &                          FACTORC,FACTORAB,2)
                    IF(NTEST.GE.3000) THEN
                      WRITE(lupri,*) ' SSCR matrix '
                      CALL WRTMT_LU(SSCR,LIKB,NIKT,LIKB,NIKT)
                    END IF
* ============================
* Loop over ik and scatter out
* ============================
*. Generate double excitations from K strings
*. I strings connected with K strings in batch <I!a+i a+k!K)
                    I12 = 2
*
                    IONE = 1
                    IKBOFF = 1
                    DO IKPAIR = 1, IKBT(2,IKBTC)
                      ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                      KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                      NI = NOBPTS(ITYP,ISM)
                      NK = NOBPTS(KTYP,KSM)
                      IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                        NIK = NI*(NI+1)/2
                        IKSM = 1
                      ELSE
                        NIK = NI * NK
                        IKSM = 0
                      END IF
CT                    CALL QENTER('ADADS')
                      IF(IFRST.EQ.1) KFRST = 1
                      ONE = 1.0D0
*
                      IAC = I4_AC(1)
                      KAC = I4_AC(2)
*
                      CALL ADAADAST_GAS(IONE,ISM,ITYP,NI,IAC,
     &                                  IONE,KSM,KTYP,NK,KAC,
     &                                ISCTP,ISCSM,IGRP,
     &                                KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                                IFRST,KFRST,I12,K12,ONE          )
*
                      IFRST = 0
                      KFRST = 0
CT                    CALL QEXIT ('ADADS')
*
CT                    CALL QENTER('MATCS')
                      I = 0
                      K = 1
                      DO IK = 1, NIK
                        CALL NXTIJ(I,K,NI,NK,IKSM,NONEW)
                        IKOFF = (K-1)*NI + I
                        ISBOFF = 1+(IKBOFF-1+IK-1)*NIBTC*NKBTC
                        IF(IKSM.EQ.1.AND.I.EQ.k) THEN
* a+ i a+i gives trivially zero
                        ELSE
                          CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                                NKBTC,I1(1,IKOFF),XI1S(1,IKOFF))
                        END IF
                      END DO
CT                    CALL QEXIT ('MATCS')
                      IKBOFF = IKBOFF + NIK
*
                    END DO
*                   ^ End of loop over IKPAIRS in batch
*
                  IF(KEND.EQ.0) GOTO 1800
*.                ^ End of loop over partitionings of resolution strings
 1801           CONTINUE
*               ^ End of loop over partionings of I strings
 1930         CONTINUE
*             ^ End of loop over batches of JL
 1940       CONTINUE
*           ^ End of loop over batches of IK
 1950     CONTINUE
*         ^ End of loop over IKOBSM
*
*
*      ==============================================
        ELSE IF(.NOT.( I4_AC(1).EQ. I4_AC(2)) ) THEN
*      ==============================================
*
*
* Three types of operators :
* a+ a  a+ a
* a+ a  a  a+
* a  a+ a+ a
*
* There first end up with
* -a+ i ak a+l aj X2(ik,jl)
*
* Number two and three end up with
* -a i a k a l aj XC(ik,jl)  ( In coulomb form)
*
*. Integrals included ??
          JLSM = 0
          IKSM = 0
          IF(NSEL2E.NE.0) THEN
            IAMOKAY=0
            IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP) THEN
              DO JSEL2E = 1, NSEL2E
                IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY=1
              END DO
            END IF
            IF(IAMOKAY.EQ.0) GOTO 2000
          END IF
*. Symmetry of allowed Double excitation,loop over excitations
          DO 2950 IKOBSM = 1, NSMOB
            JLOBSM = SXDXSX(IKOBSM,IDXSM)
            IF(JLOBSM.EQ.0) GOTO 2950
*. types + symmetries defined => K strings are defined
            KFRST = 1
            K2FRST = 1
            DO ISM = 1, NSMOB
              KSM = ADSXA(ISM,IKOBSM)
              DO JSM = 1, NSMOB
                LSM = ADSXA(JSM,JLOBSM)
                IF(NTEST.GE.2000) WRITE(lupri,*) ' ISM KSM LSM JSM',
     &          ISM,KSM,LSM,JSM
                ISCR(I4_REO(1)) = ISM
                ISCR(I4_REO(2)) = KSM
                ISCR(I4_REO(3)) = LSM
                ISCR(I4_REO(4)) = JSM
*
                ISM_ORIG = ISCR(1)
                KSM_ORIG = ISCR(2)
                LSM_ORIG = ISCR(3)
                JSM_ORIG = ISCR(4)
*
C           DO ISM_ORIG = 1, NSMOB
C             KSM_ORIG = ADSXA(ISM_ORIG,IKOBSM)
C             DO JSM_ORIG = 1, NSMOB
C               LSM_ORIG = ADSXA(JSM_ORIG,JLOBSM)
*
C               ISCR(1) = ISM_ORIG
C               ISCR(2) = KSM_ORIG
C               ISCR(3) = LSM_ORIG
C               ISCR(4) = JSM_ORIG
*
C               ISM = ISCR(I4_REO(1))
C               KSM = ISCR(I4_REO(2))
C               LSM = ISCR(I4_REO(3))
C               JSM = ISCR(I4_REO(4))
*
                NI = NOBPTS(ITYP,ISM)
                NJ = NOBPTS(JTYP,JSM)
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
                NIK = NI*NK
                NJL = NJ*NL
                IF(NIK.EQ.0.OR.NJL.EQ.0) GOTO 2803
*
                ITPSM_ORIG = (ITYP_ORIG-1)*NSMOB + ISM_ORIG
                JTPSM_ORIG = (JTYP_ORIG-1)*NSMOB + JSM_ORIG
                KTPSM_ORIG = (KTYP_ORIG-1)*NSMOB + KSM_ORIG
                LTPSM_ORIG = (LTYP_ORIG-1)*NSMOB + LSM_ORIG
*
                IF(ITPSM_ORIG.GE.KTPSM_ORIG.AND.
     &             JTPSM_ORIG.GE.LTPSM_ORIG) THEN
*
                IFIRST = 1
*. Loop over batches of I strings
                NPART = NROW/MAXI
                IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
                IF(NTEST.GE.2000)
     &          write(lupri,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
                DO 2801 IPART = 1, NPART
                  IBOT = 1+(IPART-1)*MAXI
                  ITOP = MIN(IBOT+MAXI-1,NROW)
                  NIBTC = ITOP-IBOT+1
*.Loop over batches of intermediate strings
                  KBOT = 1- MAXK
                  KTOP = 0
 2800             CONTINUE
                    KBOT = KBOT + MAXK
                    KTOP = KTOP + MAXK
*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Obtain all double excitations from this group of K strings
CT                  CALL QENTER('ADADS')
                    I12 = 1
                    K12 = 1
                    IONE = 1
*. Creation / annihilation maps , conjugated of above
                    IF(I4_AC(4).EQ.1) THEN
                      JAC = 2
                    ELSE
                      JAC = 1
                    END IF
                    IF(I4_AC(3).EQ.1) THEN
                      LAC = 2
                    ELSE
                      LAC = 1
                    END IF
C                   KFRST = 1
                    CALL ADAADAST_GAS(IONE,JSM,JTYP,NJ,JAC,
     &                                IONE,LSM,LTYP,NL,LAC,
     &                          ICCTP,ICCSM,IGRP,
     &                          KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                          JFRST,KFRST,I12,K12,SCLFAC)
                    JFRST = 0
                    KFRST = 0
*
CT                  CALL QEXIT('ADADS')
                    IF(NKBTC.EQ.0) GOTO 2801
*. Loop over jl in TS classes and gather
CT                  CALL QENTER('MATCG')
                    J = 0
                    L = 1
                    DO  IJL = 1, NJL
                      CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                      I1JL = (L-1)*NJ+J
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                      JLOFF = (IJL-1)*NKBTC*NIBTC+1
                      CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,
     &                         NKBTC,I1(1,I1JL),XI1S(1,I1JL))
                    END DO
CT                  CALL QEXIT ('MATCG')
*
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
*.Obtain two electron integrals as xint(ik,jl) = (ij!kl)-(il!kj)
                    IKSM = 0
                    JLSM = 0
                    IF(IFIRST.EQ.1) THEN
                      IF(I4_AC(1).EQ.I4_AC(3)) THEN
* a+ a a+ a
                        ICOUL = 2
                      ELSE
* a+ a a a+ or a+ a a a+
                        ICOUL = 1
                      END IF
*. Use coulomb - exchange or just coulomb integrals ?
                      IF(ITPSM_ORIG.EQ.KTPSM_ORIG
     &                .AND.JTPSM_ORIG.EQ.LTPSM_ORIG)THEN
*. No use of exchange
                        IXCHNG = 0
                        FACX = -0.5D0
                      ELSE IF(ITPSM_ORIG.NE.KTPSM_ORIG
     &                .OR.JTPSM_ORIG.NE.LTPSM_ORIG) THEN
*. Exchange used, combines two terms
                        IXCHNG = 1
                        FACX = -0.5D0
                      END IF
                      IF(ITPSM_ORIG.NE.KTPSM_ORIG
     &                .AND.JTPSM_ORIG.NE.LTPSM_ORIG)THEN
*. Exchange used, combines four terms
                        IXCHNG = 1
                        FACX = -1.0D0
                      END IF
           IF( NTEST.GE.1000) WRITE(lupri,*)
     &   ' ITPSM_ORIG,KTPSM_ORIG,JTPSM_ORIG,LTPSM_ORIG,FACX',
     &     ITPSM_ORIG,KTPSM_ORIG,JTPSM_ORIG,LTPSM_ORIG,FACX
*. fetch integrals
* we want the operator in the form a+i ak a+l aj ((ij!lk)-(ik!lj))
                      IF(ICOUL.EQ.2) THEN
                      CALL LGETINT(XINT,ITYP,ISM,JTYP,JSM,LTYP,LSM,
     &                            KTYP,KSM,IXCHNG,IKSM,JLSM,ICOUL)
                      ELSE IF (ICOUL.EQ.1) THEN
                        CALL LGETINT(XINT,ITYP,ISM,KTYP,KSM,JTYP,JSM,
     &                              LTYP,LSM,IXCHNG,IKSM,JLSM,ICOUL)
                      END IF
*
                    END IF
*                   ^ End if integrals should be fetched
                    IFIRST = 0
*.and now ,to the work
                    LIKB = NIBTC*NKBTC
                    IF(NTEST.GE.3000) THEN
                     WRITE(lupri,*) ' Integral block '
                     CALL WRTMT_LU(XINT,NIK,NJL,NIK,NJL)
                    END IF
                    IF(NTEST.GE.3000) THEN
                      WRITE(lupri,*) ' CSCR matrix '
                      CALL WRTMT_LU(CSCR,LIKB,NJL,LIKB,NJL)
                    END IF
*
C?                  MXACIJO = MXACIJ
                    MXACIJ = MAX(MXACIJ,LIKB*NJL,LIKB*NIK)
C?                  IF(MXACIJ.GT.MXACIJO) THEN
C?                    write(6,*) ' New max MXACIJ = ', MXACIJ
C?                    write(6,*) ' ISCTP,ICCTP', ISCTP,ICCTP
C?                    WRITE(6,*) ' ITYP,JTYP,KTYP,LTYP',
C?   &                             ITYP,JTYP,KTYP,LTYP
C?                    WRITE(6,*)'NIJ NJL NIBTC NKBTC',
C?   &                           NIJ,NJL,NIBTC,NKBTC
C?                  END IF
*
                    FACTORC = 0.0D0
                    FACTORAB = FACX
                    CALL MATML7(SSCR,CSCR,XINT,
     &                          LIKB,NIK,LIKB,NJL,NIK,NJL,
     &                          FACTORC,FACTORAB,2)
                    IF(NTEST.GE.3000) THEN
                      WRITE(lupri,*) ' SSCR matrix '
                      CALL WRTMT_LU(SSCR,LIKB,NIK,LIKB,NIK)
                    END IF
* ============================
* Loop over ik and scatter out
* ============================
*. Generate double excitations from K strings
*. I strings connected with K strings in batch <I!a+i a+k!K)
                    I12 = 2
*
                    IONE = 1
CT                  CALL QENTER('ADADS')
                    IF(IFRST.EQ.1) KFRST = 1
                    ONE = 1.0D0
*
                    IAC = I4_AC(1)
                    KAC = I4_AC(2)
*
C                   KFRST = 1
                    CALL ADAADAST_GAS(IONE,ISM,ITYP,NI,IAC,
     &                                IONE,KSM,KTYP,NK,KAC,
     &                              ISCTP,ISCSM,IGRP,
     &                              KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                              IFRST,KFRST,I12,K12,ONE          )
*
                    IFRST = 0
                    KFRST = 0
CT                  CALL QEXIT ('ADADS')
*
CT                  CALL QENTER('MATCS')
                    I = 0
                    K = 1
                    DO IK = 1, NIK
                      CALL NXTIJ(I,K,NI,NK,IKSM,NONEW)
                      IKOFF = (K-1)*NI + I
                      ISBOFF = 1+(IK-1)*NIBTC*NKBTC
                      CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                     NKBTC,I1(1,IKOFF),XI1S(1,IKOFF))
                    END DO
C                   write(6,*) ' first element of updated SB', SB(1)
CT                  CALL QEXIT ('MATCS')
*
                  IF(KEND.EQ.0) GOTO 2800
*. End of loop over partitionings of resolution strings
 2801           CONTINUE
*               ^ End of loop over batches of I strings
              END IF
*             ^ End of if I. ge. K, J.ge. L
 2803         CONTINUE
              END DO
*             ^ End of loop over KSM
            END DO
*           ^ End of loop over ISM
 2950     CONTINUE
        END IF
*       ^ End of a+ a+ a a/a a a+ a+ versus a+ a a+ a switch

        IF(NIJKL1.EQ.0) CALL QEXIT ('BB2A0')
        IF(NIJKL1.EQ.1) CALL QEXIT ('BB2A1')
        IF(NIJKL1.EQ.2) CALL QEXIT ('BB2A2')
        IF(NIJKL1.EQ.3) CALL QEXIT ('BB2A3')
        IF(NIJKL1.EQ.4) CALL QEXIT ('BB2A4')
 2000 CONTINUE
*
 2001 CONTINUE
      CALL QEXIT('RS2A ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                   JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                   IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                   SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                   NOBPTS,       MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                   XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                   CJRES,SIRES,SCLFAC,NTESTG,
     &                   NSEL2E,ISEL2E,IUSE_PH,IPHGAS,
     &                   MXTSOB,H,XIJILS)
*
* Combined alpha-beta double excitation
* contribution from given C block to given S block
*. If IUSAB only half the terms are constructed
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha  strings in sigma
* IBSM,IBTP : Symmetry and type of beta   strings in sigma
* JASM,JATP : Symmetry and type of alpha  strings in C
* JBSM,JBTP : Symmetry and type of beta   strings in C
* NIA,NIB : Number of alpha-(beta-) strings in sigma
* NJA,NJB : Number of alpha-(beta-) strings in C
* IAGRP : String group of alpha strings
* IBGRP : String group of beta strings
* IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
* IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
* JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
* JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
* CB   : Input C block
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
*       single excitations
* MAXK   : Largest number of inner resolution strings treated at simult.
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT  : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + CJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks
* October 96 : New routines for accessing annihilation information
*             Cleaned and shaved, only IROUTE = 3 option active
* October   97 : allowing for N-1/N+1 switch
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
      REAL*8 INPROD
*
      INTEGER ISEL2E(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*)
      DIMENSION CJRES(*),SIRES(*)
*
CTF   PARAMETER(MXTSOB = 40)
CTF   DIMENSION H(MXTSOB*MXTSOB)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
      DIMENSION IOP_TYP(2),IOP_AC(2),IOP_REO(2)
*
      DIMENSION IJ_TYP(2),IJ_DIM(2),IJ_REO(2),IJ_AC(2),IJ_SYM(2)
      DIMENSION KL_TYP(2),KL_DIM(2),KL_REO(2),KL_AC(2),KL_SYM(2)
*
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
*. Arrays for reorganization
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
C    &          SIGNREO,NADOP,NADDEL,IADDEL,ADSIGN)
*
#include "comjep.inc"
      CALL QENTER('RS2B ')
*
      NTESTL = 0
      NTEST = MAX(NTESTG,NTESTL)
*
      IF(NTEST.GE.500) THEN
        WRITE(lupri,*) ' =============== '
        WRITE(lupri,*) ' RSBB2BN speaking '
        WRITE(lupri,*) ' =============== '
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Groups defining each supergroup
      CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
      CALL GET_SPGP_INF(JATP,IAGRP,JASPGP)
      CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
      CALL GET_SPGP_INF(JBTP,IBGRP,JBSPGP)
*
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(lupri,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(lupri,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
      CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
*
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Should N-1 or N+1 projection be used for alpha strings
          IJ_TYP(1) = ITYP
          IJ_TYP(2) = JTYP
          IJ_AC(1)  = 2
          IJ_AC(2) =  1
          NOP = 2
          IF(IUSE_PH.EQ.1) THEN
            CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TYP,IJ_AC,IJ_REO,
     &           SIGNIJ)
          ELSE
*. Enforced a+ a
            IJ_REO(1) = 1
            IJ_REO(2) = 2
            SIGNIJ = 1.0D0
          END IF
*. Two choices here :
*  1 : <Ia!a+ ia!Ka><Ja!a+ ja!Ka> ( good old creation mapping)
*  2 :-<Ia!a  ja!Ka><Ja!a  ia!Ka>  + delta(i,j)
C?        WRITE(6,*) ' RSBB2BN : IOP_REO : ', (IOP_REO(II),II=1,2)
          IF(IJ_REO(1).EQ.1.AND.IJ_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
            IJAC = 2
            SIGNIJ2 = SCLFAC
*
            IJ_DIM(1) = NI
            IJ_DIM(2) = NJ
            IJ_SYM(1) = ISM
            IJ_SYM(2) = JSM
            IJ_TYP(1) = ITYP
            IJ_TYP(2) = JTYP
*
            NOP1   = NI
            IOP1SM = ISM
            IOP1TP = ITYP
            NOP2   = NJ
            IOP2SM = JSM
            IOP2TP = JTYP
          ELSE
*. Terra Nova, annihilation map
            IJAC = 1
            SIGNIJ2 = -SCLFAC
*
            IJ_DIM(1) = NJ
            IJ_DIM(2) = NI
            IJ_SYM(1) = JSM
            IJ_SYM(2) = ISM
            IJ_TYP(1) = JTYP
            IJ_TYP(2) = ITYP
*
            NOP1   = NJ
            IOP1SM = JSM
            IOP1TP = JTYP
            NOP2   = NI
            IOP2SM = ISM
            IOP2TP = ITYP
          END IF
*
*. Generate creation- or annihilation- mappings for all Ka strings
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          CALL ADAST_GAS(IJ_SYM(2),IJ_TYP(2),NGAS,JASPGP,JASM,
     &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SIGNIJ2,IJAC)
C         CALL ADAST_GAS(JSM,JTYP,JATP,JASM,IAGRP,
C    &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFACS,IJ_AC)
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          CALL ADAST_GAS(IJ_SYM(1),IJ_TYP(1),NGAS,IASPGP,IASM,
     &         I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJAC)
C         CALL ADAST_GAS(ISM,ITYP,NGAS,IASPGP,IASM,
C    &         I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJ_AC)
*. Compress list to common nonvanishing elements
          IDOCOMP = 0
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I3,XI3S,IJ_DIM(1),
     &                        NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF

*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
*
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, IJ_DIM(2)
              CALL GET_CKAJJB(CB,IJ_DIM(2),NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
            IF(NTEST.GE.500) THEN
              WRITE(lupri,*) ' Updated CJRES as C(Kaj,Jb)'
              CALL WRTMT_LU(CJRES,NKASTR*NJ,NJB,NKASTR*NJ,NJB)
            END IF
*
            MXACJ=MAX(MXACJ,NIB*LKABTC*IJ_DIM(1),NJB*LKABTC*IJ_DIM(2))
            CALL SETVEC(SIRES,ZERO,NIB*LKABTC*IJ_DIM(1))
            FACS = 1.0D0
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
*. Allowed double excitation ?
              IJKL_ACT = I_DX_ACT(ITYP,KTYP,LTYP,JTYP)
              IF(IJKL_ACT.EQ.0) GOTO 2000
              IF(NTEST.GE.100) THEN
                WRITE(lupri,*) ' KTYP, LTYP', KTYP, LTYP
              END IF
*. Should this group of excitations be included
              IF(NSEL2E.NE.0) THEN
               IAMOKAY=0
               IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP)THEN
                 DO JSEL2E = 1, NSEL2E
                   IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY = 1
                 END DO
               END IF
               IF(IAMOKAY.EQ.0) GOTO 2000
              END IF
*
              KL_TYP(1) = KTYP
              KL_TYP(2) = LTYP
              KL_AC(1)  = 2
              KL_AC(2) =  1
              NOP = 2
              IF(IUSE_PH.EQ.1) THEN
                CALL ALG_ROUTERX(IBOC,JBOC,NOP,KL_TYP,KL_AC,KL_REO,
     &               SIGNKL)
              ELSE
*. Enforced a+ a
                KL_REO(1) = 1
                KL_REO(2) = 2
                SIGNKL = 1.0D0
              END IF
*
              DO 1930 KSM = 1, NSMOB
                IFIRST = 1
                LSM = ADSXA(KSM,KLSM)
                IF(NTEST.GE.100) THEN
                  WRITE(lupri,*) ' KSM, LSM', KSM, LSM
                END IF
                IF(LSM.EQ.0) GOTO 1930
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*
                IF(KL_REO(1).EQ.1.AND.KL_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
                  KLAC = 2
                  KL_DIM(1) = NK
                  KL_DIM(2) = NL
                  KL_SYM(1) = KSM
                  KL_SYM(2) = LSM
                  KL_TYP(1) = KTYP
                  KL_TYP(2) = LTYP
                ELSE
*. Terra Nova, annihilation map
                  KLAC = 1
                  KL_DIM(1) = NL
                  KL_DIM(2) = NK
                  KL_SYM(1) = LSM
                  KL_SYM(2) = KSM
                  KL_TYP(1) = LTYP
                  KL_TYP(2) = KTYP
                END IF
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
*. currently we are using creation mappings for kl
*. (Modify to use ADAST later )
                CALL ADAST_GAS(KL_SYM(2),KL_TYP(2),NGAS,JBSPGP,JBSM,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,SIGNKL,KLAC)
C               CALL ADSTN_GAS(LSM,LTYP,JBTP,JBSM,IBGRP,
C    &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADAST_GAS(KL_SYM(1),KL_TYP(1),NGAS,IBSPGP,IBSM,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,KLAC)
C               CALL ADSTN_GAS(KSM,KTYP,IBTP,IBSM,IBGRP,
C    &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*
* Fetch Integrals as (iop2 iop1 |  k l )
*
                IXCHNG = 0
                ICOUL = 1
                CALL LGETINT(XINT,IJ_TYP(2),IJ_SYM(2),
     &               IJ_TYP(1),IJ_SYM(1),
     &               KL_TYP(1),KL_SYM(1),KL_TYP(2),KL_SYM(2),IXCHNG,
     &               0,0,ICOUL)
C               CALL LGETINT(XINT,JTYP,JSM,ITYP,ISM,KTYP,KSM,
C    &                     LTYP,LSM,IXCHNG,0,0,ICOUL)
*
* S(Ka,i,Ib) = sum(j,k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)
*
                IJKL_DIM = IJ_DIM(1)*IJ_DIM(2)*KL_DIM(1)*KL_DIM(2)
                IF(INPROD(XINT,XINT,IJKL_DIM).NE.0.0D0) THEN
                IROUTE = 3
                CALL SKICKJ(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,XINT,IJ_DIM(1),IJ_DIM(2),
     &               KL_DIM(1),KL_DIM(2),
     &               NKBSTR,I4,XI4S,I2,XI2S,IKORD,
     &               FACS,IROUTE,MXTSOB,XIJILS)
                END IF
*
                IF(NTEST.GE.500) THEN
                  WRITE(lupri,*) ' Updated Sires as S(Kai,Ib)'
                  CALL WRTMT_LU(SIRES,LKABTC*NI,NIB,LKABTC*NI,NIB)
                END IF
*
 1930         CONTINUE
*             ^ End of loop over KSM
 2000       CONTINUE
*           ^ End of loop over KLTYP
*
*. Scatter out from s(Ka,Ib,i)
*
            IF(NTEST.GE.1000) THEN
              WRITE(lupri,*) ' S(Ka,Ib,i) as S(Ka,Ibi)'
              CALL WRTMT_LU(SIRES,LKABTC,NIB*IJ_DIM(1),LKABTC,IJ_DIM(1))
            END IF
*
            DO II = 1, IJ_DIM(1)
              CALL ADD_SKAIIB(SB,IJ_DIM(1),NIA,SIRES,LKABTC,NIB,II,
     &             I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
 1801     CONTINUE
*.        ^End of loop over partitioning of alpha strings
 1940   CONTINUE
*       ^ End of loop over ISM
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      CALL QEXIT('RS2B ')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE RSBB2BVN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                   JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                   IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                   SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                   NOBPTS,       MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                   XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                   CJRES,SIRES,SCLFAC,NTESTG,
     &                   NSEL2E,ISEL2E,IUSE_PH,IPHGAS,
     &                   SIRESPA,CJRESPA,MXTSOB,H,XIJILS)
*
* Combined alpha-beta double excitation
* contribution from given C block to given S block
*. If IUSAB only half the terms are constructed
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha  strings in sigma
* IBSM,IBTP : Symmetry and type of beta   strings in sigma
* JASM,JATP : Symmetry and type of alpha  strings in C
* JBSM,JBTP : Symmetry and type of beta   strings in C
* NIA,NIB : Number of alpha-(beta-) strings in sigma
* NJA,NJB : Number of alpha-(beta-) strings in C
* IAGRP : String group of alpha strings
* IBGRP : String group of beta strings
* IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
* IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
* JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
* JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
* CB   : Input C block
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
*       single excitations
* MAXK   : Largest number of inner resolution strings treated at simult.
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT  : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + CJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks
* October 96 : New routines for accessing annihilation information
*             Cleaned and shaved, only IROUTE = 3 option active
* October   97 : allowing for N-1/N+1 switch
* March 98   : Allows for splitting of strings into active and passive groups
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
*
      INTEGER ISEL2E(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*)
      DIMENSION CJRES(*),SIRES(*)
*
CTF   PARAMETER(MXTSOB = 40)
CTF   DIMENSION H(MXTSOB*MXTSOB)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
      DIMENSION IOP_TYP(2),IOP_AC(2),IOP_REO(2)
*
      DIMENSION IJ_TYP(2),IJ_DIM(2),IJ_REO(2),IJ_AC(2),IJ_SYM(2)
      DIMENSION KL_TYP(2),KL_DIM(2),KL_REO(2),KL_AC(2),KL_SYM(2)
*
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
*. Arrays for reorganization
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
*. Arrays for active/passive division
      INTEGER IACIA(20),IPAIA(20),IACIB(20),IPAIB(20)
      INTEGER IACJA(20),IPAJA(20),IACJB(20),IPAJB(20)
*
C     IBSPGP,NGAS,IBSM,NIB,2,KL_TYP,NIBAC_S,NIBPA_S,IREO_IB,IBREO_IB)
      INTEGER NIBAC_S(8),NIBPA_S(8),IBREO_IB(8)
      INTEGER NJBAC_S(8),NJBPA_S(8),IBREO_JB(8)
*. Dimension ^ : Number of string symmetries
*. The following must be moved outside !
      INTEGER IREO_IB(20000),IREO_JB(20000)
*.    ^ Dimension  : Largest number of strings of given sym
C     DIMENSION SIRESPA(10000),CJRESPA(10000)
      DIMENSION SIRESPA(*),CJRESPA(*)
*     ^ Dimension : Same as SIRES and CJRES

#include "comjep.inc"
      CALL QENTER('RS2B ')
      NTESTL = 0
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ================= '
        WRITE(6,*) ' RSBB2BVN speaking '
        WRITE(6,*) ' ================= '
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Use passive/active splitting ?
      IUSE_PA = 1
*. Groups defining each supergroup
      CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
      CALL GET_SPGP_INF(JATP,IAGRP,JASPGP)
      CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
      CALL GET_SPGP_INF(JBTP,IBGRP,JBSPGP)
*
      IF(IUSE_PA .EQ. 0 ) THEN
*. No use of active/passive division
        NPAIA = 0
        NPAJA = 0
        NPAIB = 0
        NPAJB = 0
*
        NACIA = NGAS
        NACJA = NGAS
        NACIB = NGAS
        NACJB = NGAS
*
        DO IGAS = 1, NGAS
         IACIA(IGAS) = IASPGP(IGAS)
         IACIB(IGAS) = IBSPGP(IGAS)
         IACJA(IGAS) = JASPGP(IGAS)
         IACJB(IGAS) = JBSPGP(IGAS)
        END DO
      END IF
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
      CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
*
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Should N-1 or N+1 projection be used for alpha strings
          IJ_TYP(1) = ITYP
          IJ_TYP(2) = JTYP
          IJ_AC(1)  = 2
          IJ_AC(2) =  1
          NOP = 2
          IF(IUSE_PH.EQ.1) THEN
            CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TYP,IJ_AC,IJ_REO,
     &           SIGNIJ)
          ELSE
*. Enforced a+ a
            IJ_REO(1) = 1
            IJ_REO(2) = 2
            SIGNIJ = 1.0D0
          END IF
*. Two choices here :
*  1 : <Ia!a+ ia!Ka><Ja!a+ ja!Ka> ( good old creation mapping)
*  2 :-<Ia!a  ja!Ka><Ja!a  ia!Ka>  + delta(i,j)
          IF(IJ_REO(1).EQ.1.AND.IJ_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
            IJAC = 2
            SIGNIJ2 = SCLFAC
*
            IJ_DIM(1) = NI
            IJ_DIM(2) = NJ
            IJ_SYM(1) = ISM
            IJ_SYM(2) = JSM
            IJ_TYP(1) = ITYP
            IJ_TYP(2) = JTYP
*
            NOP1   = NI
            IOP1SM = ISM
            IOP1TP = ITYP
            NOP2   = NJ
            IOP2SM = JSM
            IOP2TP = JTYP
          ELSE
*. Terra Nova, annihilation map
            IJAC = 1
            SIGNIJ2 = -SCLFAC
*
            IJ_DIM(1) = NJ
            IJ_DIM(2) = NI
            IJ_SYM(1) = JSM
            IJ_SYM(2) = ISM
            IJ_TYP(1) = JTYP
            IJ_TYP(2) = ITYP
*
            NOP1   = NJ
            IOP1SM = JSM
            IOP1TP = JTYP
            NOP2   = NI
            IOP2SM = ISM
            IOP2TP = ITYP
          END IF
*
*. Generate creation- or annihilation- mappings for all Ka strings
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          CALL ADAST_GAS(IJ_SYM(2),IJ_TYP(2),NGAS,JASPGP,JASM,
     &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SIGNIJ2,IJAC)
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          CALL ADAST_GAS(IJ_SYM(1),IJ_TYP(1),NGAS,IASPGP,IASM,
     &         I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJAC)
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I3,XI3S,IJ_DIM(1),
     &                        NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF

*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
*
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, IJ_DIM(2)
              CALL GET_CKAJJB(CB,IJ_DIM(2),NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*
            MXACJ=MAX(MXACJ,NIB*LKABTC*IJ_DIM(1),NJB*LKABTC*IJ_DIM(2))
C           IF(IUSE_PA.EQ.0) THEN
              CALL SETVEC(SIRES,ZERO,NIB*LKABTC*IJ_DIM(1))
C           END IF
C           IF(IUSE_PA.EQ.1) THEN
C             CALL SETVEC(SIRESPA,ZERO,NIB*LKABTC*IJ_DIM(1))
C           END IF
            FACS = 1.0D0
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
              IF(NTEST.GE.100) THEN
                WRITE(6,*) ' KTYP, LTYP', KTYP, LTYP
              END IF
*. Should this group of excitations be included
              IF(NSEL2E.NE.0) THEN
               IAMOKAY=0
               IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP)THEN
                 DO JSEL2E = 1, NSEL2E
                   IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY = 1
                 END DO
               END IF
               IF(IAMOKAY.EQ.0) GOTO 2000
              END IF
*
              KL_TYP(1) = KTYP
              KL_TYP(2) = LTYP
              KL_AC(1)  = 2
              KL_AC(2) =  1
              NOP = 2
              IF(IUSE_PH.EQ.1) THEN
                CALL ALG_ROUTERX(IBOC,JBOC,NOP,KL_TYP,KL_AC,KL_REO,
     &               SIGNKL)
              ELSE
*. Enforced a+ a
                KL_REO(1) = 1
                KL_REO(2) = 2
                SIGNKL = 1.0D0
              END IF
                IF(IUSE_PA.GT.0) THEN
*. Split IB strings into active/passive part
                  CALL REO_STR_SPGRP3(IBSPGP,NGAS,IBSM,NIB,2,KL_TYP,
     &                 NACIB,IACIB,NIBAC_S,NIBPA_S,IBREO_IB,IREO_IB,
     &                 SIGNPAI)
*. Split JB strings into active/passive part
                  CALL REO_STR_SPGRP3(JBSPGP,NGAS,JBSM,NJB,2,KL_TYP,
     &                 NACJB,IACJB,NJBAC_S,NJBPA_S,IBREO_JB,IREO_JB,
     &                 SIGNPAJ)
*. Reorganize C(Ka,j,Jb) to C(Ka,Jb_pa,j,Jb_pa)
                  XDUM = 0.0D0
                  IDUM = 1
                  CALL CKAJJB_PA(CJRES,CJRESPA,1,LKABTC,IJ_DIM(2),NJB,
     &                           IREO_JB,NJBPA_S,NJBAC_S,NSMST,XDUM,
     &                           IDUM)
                  CALL SETVEC(SIRESPA,ZERO,NIB*LKABTC*IJ_DIM(1))
                  I_ADD_COPY = 2
                END IF
*
              DO 1930 KSM = 1, NSMOB
*
                IFIRST = 1
                LSM = ADSXA(KSM,KLSM)
                IF(NTEST.GE.100) THEN
                  WRITE(6,*) ' KSM, LSM', KSM, LSM
                END IF
                IF(LSM.EQ.0) GOTO 1930
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*
                IF(KL_REO(1).EQ.1.AND.KL_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
                  KLAC = 2
                  KL_DIM(1) = NK
                  KL_DIM(2) = NL
                  KL_SYM(1) = KSM
                  KL_SYM(2) = LSM
                  KL_TYP(1) = KTYP
                  KL_TYP(2) = LTYP
                ELSE
*. Terra Nova, annihilation map
                  KLAC = 1
                  KL_DIM(1) = NL
                  KL_DIM(2) = NK
                  KL_SYM(1) = LSM
                  KL_SYM(2) = KSM
                  KL_TYP(1) = LTYP
                  KL_TYP(2) = KTYP
                END IF
*
*
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Loop over symmetries of active strings
                IF(IUSE_PA.EQ.0) THEN
                  JBSMMN = JBSM
                  JBSMMX = JBSM
                ELSE
                  JBSMMN = 1
                  JBSMMX = NSMST
                END IF
                DO JBSM_AC = JBSMMN,JBSMMX
                  IF(NTEST.GE.1000) WRITE(6,*) 'JBSM_AC=',JBSM_AC
                  IF(IUSE_PA.EQ.1) THEN
*. Symmetry of active part of I_b string
                    DO IIBSM_AC = 1, NSMST
                      IF(STSTSX(IIBSM_AC,JBSM_AC).EQ.KLSM)
     &                IBSM_AC = IIBSM_AC
                    END DO
*. Block in C(Ka,Jb_pa,j,Jb_ac)
                    IF(JBSM_AC.EQ.1) THEN
                      ICK_OFF = 1
                    ELSE
                      ICK_OFF = ICK_OFF
     &              + LKABTC*IJ_DIM(2)*
     &                NJBPA_S(JBSM_AC-1)*NJBAC_S(JBSM_AC-1)
                    END IF
*. Offset in S(Ka,Ib_pa,i,Ib_ac)
                    ISK_OFF = 1
                    DO IIBSM_AC = 1, IBSM_AC-1
                      ISK_OFF = ISK_OFF
     &              + LKABTC*IJ_DIM(1)
     &              *NIBPA_S(IIBSM_AC)*NIBAC_S(IIBSM_AC)
                    END DO
                    IF(NTEST.GE.1000) THEN
                      WRITE(6,*) ' JBSM_AC,IBSM_AC', JBSM_AC,IBSM_AC
                      WRITE(6,*) ' ISK_OFF = ', ISK_OFF
                    END IF
                  IF(NIBPA_S(IBSM_AC)*NIBAC_S(IBSM_AC).EQ.0)GOTO 1912
                  IF(NJBPA_S(JBSM_AC)*NJBAC_S(JBSM_AC).EQ.0)GOTO 1912
*
                  END IF

                IF(IUSE_PA.EQ.0) THEN
                CALL ADAST_GAS(KL_SYM(2),KL_TYP(2),NGAS,JBSPGP,JBSM,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,SIGNKL,KLAC)
                ELSE
                CALL ADAST_GAS(KL_SYM(2),KL_TYP(2),NACJB,IACJB,JBSM_AC,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,SIGNKL,KLAC)
                END IF
                IF(NTEST.GE.2000) WRITE(6,*) ' NKBSTR = ', NKBSTR
                IF(NKBSTR.EQ.0) GOTO 1912
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                IF(IUSE_PA.EQ.0) THEN
                CALL ADAST_GAS(KL_SYM(1),KL_TYP(1),NGAS,IBSPGP,IBSM,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,KLAC)
                ELSE
                CALL ADAST_GAS(KL_SYM(1),KL_TYP(1),NACIB,IACIB,IBSM_AC,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,KLAC)
                END IF
*
* Fetch Integrals as (iop2 iop1 |  k l )
*
                IXCHNG = 0
                ICOUL = 1
                CALL LGETINT(XINT,IJ_TYP(2),IJ_SYM(2),
     &               IJ_TYP(1),IJ_SYM(1),
     &               KL_TYP(1),KL_SYM(1),KL_TYP(2),KL_SYM(2),IXCHNG,
     &               0,0,ICOUL)
*
* S(Ka,j,Ib) = sum(k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)
*
                IROUTE = 3
                IF(IUSE_PA.EQ.0) THEN
                  CALL SKICKJ(SIRES,CJRES,LKABTC,NIB,NJB,
     &                 NKBSTR,XINT,IJ_DIM(1),IJ_DIM(2),
     &                 KL_DIM(1),KL_DIM(2),
     &                 NKBSTR,I4,XI4S,I2,XI2S,IKORD,
     &                 FACS,IROUTE,MXTSOB,XIJILS)
                ELSE
                  CALL SKICKJ(SIRESPA(ISK_OFF),CJRESPA(ICK_OFF),
     &                 LKABTC*NJBPA_S(JBSM_AC),
     &                 NIBAC_S(IBSM_AC),NJBAC_S(JBSM_AC),
     &                 NKBSTR,XINT,IJ_DIM(1),IJ_DIM(2),
     &                 KL_DIM(1),KL_DIM(2),
     &                 NKBSTR,I4,XI4S,I2,XI2S,IKORD,
     &                 FACS,IROUTE,MXTSOB,XIJILS)
                END IF
*               ^ End of switch IUSE_PA
 1912         CONTINUE
              END DO
*             ^ End of loop over JBSM_AC
 1930         CONTINUE
*             ^ End of loop over KSM
              IF(IUSE_PA.EQ.1) THEN
                SIGNPA = SIGNPAI*SIGNPAJ
C?              WRITE(6,*) SIGNPA
                CALL CKAJJB_PA(SIRES,SIRESPA,2,LKABTC,IJ_DIM(1),
     &                      NIB,IREO_IB,
     &                      NIBPA_S,NIBAC_S,NSMST,SIGNPA,I_ADD_COPY)
                I_ADD_COPY = 1
              END IF
 2000       CONTINUE
*           ^ End of loop over KLTYP
*. Scatter out from s(Ka,i,Ib)
*
            IF(NTEST.GE.1000) THEN
              WRITE(6,*) ' S(Ka,i,Ib) as S(Kai,Ibi)'
              CALL WRTMT_LU(SIRES,LKABTC*IJ_DIM(1),NIB,
     &                          LKABTC*IJ_DIM(1),NIB)
            END IF
*
            DO II = 1, IJ_DIM(1)
              CALL ADD_SKAIIB(SB,IJ_DIM(1),NIA,SIRES,LKABTC,NIB,II,
     &             I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
 1801     CONTINUE
*.        ^End of loop over partitioning of alpha strings
 1940   CONTINUE
*       ^ End of loop over ISM
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      CALL QEXIT('RS2B ')
      RETURN
      END
***********************************************************************
      SUBROUTINE RSSBCB2(IASM,IATP,IOCPTA,IBSM,IBTP,IOCTPB,
     &                  JASM,JATP,JBSM,JBTP,NGAS,
     &                  IAOC,IBOC,JAOC,JBOC,NAEL,NBEL,
     &                  IJAGRP,IJBGRP,SB,CB,JDOH2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &                  NSMOB,NSMST,NSMSX,NSMDX,
     &                  NIA,NIB,NJA,NJB,MXPOBS,IDC,PS,
     &                  ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &                  MXSXBL,MOCAA,MOCAB,IAPR,IPRNT,IHAPR,
     &                  IPTSPC,JPTSPC,IJOP,NNSEL2E,ISEL2E,SCLFAC,
     &                  IUSE_PH,IPHGAS,I_RES_AB,IUSE_PA,CJPA,SIPA,
     &                  MXTSOB,SOMESCR,SOMEH,XIJILS,luwrt)
*
* Contributions to sigma block (iasm iatp, ibsm ibtp ) from
* C block (jasm jatp , jbsm, jbtp)
*
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha strings in sigma
* IBSM,IBTP : Symmetry and type of beta  strings in sigma
* JASM,JATP : Symmetry and type of alpha strings in C
* JBSM,JBTP : Symmetry and type of beta  strings in C
* NGAS      : Number of active spaces in calculation
* IAOC,IBOC : Number of electrons in each AS for sigma supergroups
* JAOC,JBOC : Number of electrons in each AS for C     supergroups
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* IDOH2 : = 0 => no two electron operator
* IDOH2 : = 1 =>    two electron operator
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
*
* IHAPR : .ne. 0 implies thatt the exact Hamiltonian shoulf not be uses
* In the case IPTSPC and JPTSPC defined the perturbation spaces
* a nonvanishing perturbation is allowed inside each subspace.
* The actual type of approximate Hamiltonian in each subspace is defined by
* IHFORM
* NNSEL2E : Only selected 2e terms will be included
* ISEL2E : orbital spaces in which 2e terms are included
*          (Currently : all indeces identical )
*
* ======
* Output
* ======
* SB : fresh sigma block
*
* =======
* Scratch
* =======
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* C2 : Must hold largest STT block of sigma or C
*
* XINT : Scratch space for integrals.
*
* Jeppe Olsen , Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
*. Output
      DIMENSION CB(*),SB(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*)
      DIMENSION C2(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION IBLOCK(8)
      DIMENSION IPHGAS(*)
*. For H(apr)
      DIMENSION ISEL2E(*)
*
      NTEST = 0000
      NTEST = MAX(NTEST,IPRNT)

*
      IF(NTEST.GE.200) THEN
        WRITE(luwrt,*) ' ==============================='
        WRITE(luwrt,*) ' RSSBCB2 :  C block (transposed)'
        WRITE(luwrt,*) ' ================================'
        CALL WRTMAT(CB,NJB,NJA,NJB,NJA)
        WRITE(luwrt,*) ' ======================================='
        WRITE(luwrt,*) ' RSSBCB2 : Initial  S block(transposed) '
        WRITE(luwrt,*) ' ======================================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
        WRITE(luwrt,*) ' Overall scalefactor ',SCLFAC
        WRITE(luwrt,*) ' IHAPR,JDOH2 = ', IHAPR,JDOH2
        WRITE(luwrt,*) ' IUSE_PH,I_RES_AB = ', IUSE_PH,I_RES_AB
      END IF
*
      IF(NTEST.GE.500) THEN
        WRITE(luwrt,*) ' IAOC and IBOC '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(luwrt,*) ' JAOC and JBOC  : '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(luwrt,*) ' IASM IATP JASM JATP ', IASM,IATP,JASM,JATP
        WRITE(luwrt,*) ' IBSM IBTP JBSM JBTP ', IBSM,IBTP,JBSM,JBTP
        WRITE(luwrt,*) ' NAEL NBEL ', NAEL, NBEL
      END IF
* Should the corresponding Hamiltonian matrix block be
* calculated exactly or approximately
      IF(IHAPR.NE.0) THEN
        CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &               IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
     &               IMZERO,IDIAG)
        IF(NTEST.GE. 20)
     &  WRITE(luwrt,*) ' RSSBCBN : ', NNSEL2E, ISEL2E(1)
        NSEL2E = NNSEL2E
        IF(IMZERO.NE.0) GOTO 9999
      ELSE
*. Operator specified by input
        IAPRLEV =-1
        IDOH2 = JDOH2
        IDIAG = 0
        NSEL2E = 0
      END IF
      IF(NTEST.GE. 20)
     &WRITE(luwrt,*) ' IHAPR, IDIAG IDOH2 ' , IHAPR,IDIAG, IDOH2
*
*
      IF(IDC.EQ.2.AND.IATP.EQ.IBTP.AND.IASM.EQ.IBSM .AND.
     &   I_RES_AB.EQ.0.AND.JASM.EQ.JBSM.AND.JATP.EQ.JBTP) THEN
*. Diagonal sigma block, use alpha-beta symmetry to reduce computations.
        IUSEAB = 1
      ELSE
        IUSEAB = 0
      END IF
*
      IF(IDIAG.EQ.0) THEN
*
* Calculate block exactly
*
      IF(I_RES_AB.NE.1.AND.IUSEAB.EQ.0.
     &   AND.IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
* Sigma beta beta contribution
* =============================
*
* Sigma aa(IA,IB) = sum(i.gt.k,j.gt.l)<IB!Eb(ij)Eb(kl)!JB>
*                 * ((ij!kl)-(il!kj)) C(IA,JB)
*                 + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)
*.One electron part
          CALL TRPMT3(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)

        IF(NBEL.GE.0) THEN
           IF(NTEST.GE.500) THEN
             WRITE(luwrt,*) ' SB before RSBB1E'
             call wrtmat(sb,nia,nib,nia,nib)
           END IF
          IF(NTEST.GE.101)
     &    WRITE(luwrt,*) ' I am going to call RSBB1E'
          CALL RSBB1E(IBSM,IBTP,IOCTPB,JBSM,JBTP,IOCTPB,
     &                IJBGRP,NIA,NGAS,IBOC,JBOC,
     &                SB,CB,ADSXA,SXSTST,STSTSX,
     &                MXPNGAS,NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                NSMOB,NSMST,NSMSX,MXPOBS,MOCAA,
     &                NIB,SCLFAC,IUSE_PH,IPHGAS,NTEST)
           IF(NTEST.GE.500) THEN
             WRITE(luwrt,*) ' SB after RSBB1E'
             call wrtmat(sb,nib,nia,nib,nia)
           END IF
        END IF
        IF(IDOH2.NE.0.AND.NBEL.GE.0) THEN
*. Two electron part
          IF(NTEST.GE.101)
     &    WRITE(luwrt,*) ' I am going to call RSBB2A'
          CALL RSBB2A(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NIB,
     &                NGAS,IBOC,JBOC,SB,CB,
     &                ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGAS,
     &                NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &                CJRES,SIRES,MXSXBL,MOCAA,SCLFAC,NTEST,
     &                0,0,IUSE_PH,IPHGAS,MXTSOB,SOMESCR)
           IF(NTEST.GE.500) THEN
             WRITE(luwrt,*) ' SB after RSBB2a'
             call wrtmat(sb,nib,nia,nib,nia)
           END IF
        END IF
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
      END IF
*
* =============================
* Sigma alpha beta contribution
* =============================
*
      IF(IDOH2.NE.0.AND.NAEL.GE.0.AND.NBEL.GE.0) THEN
        IF(NTEST.GE.101)
     &  WRITE(luwrt,*) ' I am going to call RSBB2B'
        IIITRNS = 1
        IF(IIITRNS.EQ.1) THEN
*. Call advice routine
C     ADVICE_SIGMA(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,LADVICE)
           CALL ADVICE_SIGMA(IAOC,IBOC,JAOC,JBOC,1,LADVICE)
*. LADVICE = 2 => implies transpose
           IF(LADVICE.EQ.2) THEN
             JJJTRNS = 1
           ELSE
             JJJTRNS = 0
           END IF
        ELSE
           JJJTRNS = 0
        END IF
*
C       WRITE(luwrt,*) ' IUSE_PA = ', IUSE_PA
*
        IF (JJJTRNS.EQ.0) THEN
          IF( IUSE_PA.EQ.0 ) THEN
          CALL RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &         JASM,JATP,JBSM,JBTP,NJA,NJB,
     &         IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &         SB,CB,ADSXA,STSTSX,MXPNGAS,
     &         NOBPTS,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &         XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &         IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &         MXTSOB,SOMEH,XIJILS)
          ELSE IF ( IUSE_PA.EQ.1 ) THEN
C         WRITE(luwrt,*) ' RSBB2BVN will be called '
          CALL RSBB2BVN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &         JASM,JATP,JBSM,JBTP,NJA,NJB,
     &         IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &         SB,CB,ADSXA,STSTSX,MXPNGAS,
     &         NOBPTS,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &         XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &         IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &         CJPA,SIPA,MXTSOB,SOMEH,XIJILS)
          END IF
*
         ELSE IF ( JJJTRNS.EQ.1) THEN
*. well lets give the transpose routine some more practice : Transpose back
          CALL TRPMT3(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
*
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
C         WRITE(luwrt,*) ' RSSBCB2 : Transpose path choosen'
*
          IF(IUSE_PA.EQ.0) THEN
*. No division into active/passive
          CALL RSBB2BN(IBSM,IBTP,IASM,IATP,NIB,NIA,
     &         JBSM,JBTP,JASM,JATP,NJB,NJA,
     &         IJBGRP,IJAGRP,NGAS,IBOC,IAOC,JBOC,JAOC,
     &         SB,CB,ADSXA,STSTSX,MXPNGAS,
     &         NOBPTS,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &         XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &         IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &         MXTSOB,SOMEH,XIJILS)
          ELSE
*. Divide into active/passive
          CALL RSBB2BVN(IBSM,IBTP,IASM,IATP,NIB,NIA,
     &         JBSM,JBTP,JASM,JATP,NJB,NJA,
     &         IJBGRP,IJAGRP,NGAS,IBOC,IAOC,JBOC,JAOC,
     &         SB,CB,ADSXA,STSTSX,MXPNGAS,
     &         NOBPTS,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &         XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &         IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS,
     &         CJPA,SIPA,MXTSOB,SOMEH,XIJILS)
           END IF

*. Transpose ( To compensate later transposition )
          CALL TRPMT3(SB,NIA,NIB,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJA,NJB,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
        END IF
           IF(NTEST.GE.101) THEN
             WRITE(luwrt,*) ' SB after RSBB2B, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.500) THEN
             WRITE(luwrt,*) ' SB after RSBB2b'
             call wrtmat(sb,nib,nia,nib,nia)
           END IF
      END IF
*
* =============================
* Sigma alpha alpha contribution
* =============================
*
      IF(I_RES_AB.NE.-1.AND.
     &   NAEL.GE.0.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* alpha single excitation
*
        IF(NTEST.GE.101)
     &  WRITE(luwrt,*) ' I am going to call RSBB1E (last time )'
        CALL RSBB1E(IASM,IATP,IOCTPA,JASM,JATP,IOCTPA,
     &              IJAGRP,NIB,
     &              NGAS,IAOC,JAOC,
     &              SB,CB,
     &              ADSXA,SXSTST,STSTSX,
     &              MXPNGAS,NOBPTS,IOBPTS,
     &              ITSOB,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &              NSMOB,NSMST,NSMSX,MXPOBS,MOCAA,
     &              NIA,SCLFAC,IUSE_PH,IPHGAS,NTEST)
           IF(NTEST.GE.101) THEN
             WRITE(luwrt,*) ' SB transposed after RSBB1, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.500) THEN
             WRITE(luwrt,*) ' SB transposed  after RSBB1E'
             call wrtmat(SB,nib,nia,nib,nia)
           END IF
*
* alpha double excitation
*
        IF(IDOH2.NE.0.AND.NAEL.GE.0) THEN
          IF(NTEST.GE.101)
     &    WRITE(luwrt,*) ' I am going to call RSBB2A (last time )'
          CALL RSBB2A(IASM,IATP,JASM,JATP,IJAGRP,NIB,NIA,
     &                NGAS,IAOC,JAOC,
     &                SB,CB,
     &                ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGAS,
     &                NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &                CJRES,SIRES,MXSXBL,MOCAA,SCLFAC,NTEST,
     &                0,0,IUSE_PH,IPHGAS,MXTSOB,SOMESCR)
        END IF
*
           IF(NTEST.GE.101) THEN
             WRITE(luwrt,*)' SB transposed after RSBB2A, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
        IF(NTEST.GE.500) THEN
          WRITE(luwrt,*) ' SB after RSBB2A'
          call wrtmat(sb,nia,nib,nia,nib)
        END IF
      END IF
*
      ELSE IF (IDIAG.EQ.1) THEN
*
*. Diagonal approxiation
*
       IBLOCK(1) = IATP
       IBLOCK(2) = IBTP
       IBLOCK(3) = IASM
       IBLOCK(4) = IBSM
       IBLOCK(5) = 1
       IBLOCK(6) = 1
       IF(IDOH2.EQ.0) THEN
         I12 = 1
       ELSE
         I12 = 2
       END IF
C?     WRITE(luwrt,*) ' IDOH2, I12 ', IDOH2,I12
       ITASK = 2
       FACTOR = 0.0D0
*. Input is in det basis
       IIDC = 1
*. Well, we are working with transposed matrices so
       CALL TRPMT3(CB,NIB,NIA,C2)
C?     WRITE(luwrt,*) ' DIATERM2_GAS will be called '
       CALL DIATERM2_GAS(FACTOR,ITASK,C2,1,IBLOCK,1,0,I12,IIDC)
       CALL TRPMT3(C2,NIA,NIB,CB)
       ONE = 1.0D0
       IF(IUSEAB.EQ.0) THEN
         FACTOR = 1.0D0*SCLFAC
       ELSE
         FACTOR = 0.5D0*SCLFAC
       END IF
       CALL VECSUM(SB,SB,CB,ONE,FACTOR,NIA*NIB)
      END IF
*
 9999 CONTINUE
      IF(IHAPR.NE.0) THEN
*. Clean up
        CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &               IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
     &               IMZERO,IDIAG)
      END IF
*
      IF(NTEST.GE.200) THEN
        WRITE(luwrt,*) ' ==================================='
        WRITE(luwrt,*) ' RSSBCB : Final S block (transposed)'
        WRITE(luwrt,*) ' ==================================='
        CALL WRTMAT(SB,NIB,NIA,NIB,NIA)
      END IF
      NTESTO = NTEST
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
     &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,
     &                  SSCLFAC,ISCLFAC_GROUP
#ifdef VAR_MPI
     &                 ,LUCLIST,NPARBLOCK,IBLOCKL
#endif
     &                  )
*
* Generate a set of sigma blocks,
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order
* If LUCBLK.GT.0, the blocks of C corresponding to IBLOCK
* are stored on LUCBLK
*
* ICONSPA,ICONSPB  added October 1996
* ICBAT_RES, ICBAT_INI, IBBAT_END added august 1997
*
* If ICBAT_RES .eq.1 then it as assumed that only
* Cbatches ICBAT_INI to ICBAT_END are stored on  LUC
*
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KSIOIO, KSVST, KLSLBT, KLSLEBT, KLSI1BT, KLSIBT, KSBLTP,
     &          KISCLFAC, KSCALLOC, KLLBT, KLLEBT, KLI1BT, KLIBT, 
     &          KLSCLFAC, KI1, KXI1S, KI2, KXI2S, KI3, KXI3S, KI4,KXI4S,
     &          KC2, KCJRES, KSIRES, KSSCR, KCSCR, KCJPA, KSIPA, KSOMEH,
     &          KSOMESCR, KCONSPA, KCONSPB, KSTSTS, KSTSTD, KINSCR, 
     &          KCIOIO, KCBLTP
      INTEGER*8 KLOCSTR, KLREO, KLZ, KLZSCR
!               for addressing of WORK
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
*./ORBINP/ : NACOB used
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
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      dimension XIJILS(MXTSOB)
      INTEGER MLTHCHCK, ITEMPP
      DIMENSION SSCLFAC(*),ISCLFAC_GROUP(*)
*
*     set local mark memory
* 
      IDUM = 0
      CALL MEMMAN(KDUM,IDUM,'MARK  ',IDUM,'SBLOCK')
*
CTF
*     scratch space for rsbb2a, rsbb2bn, rsbb2bvn, skickj
CTF
      LSOMESQ = MXTSOB*MXTSOB
      call memman(KSOMEH,LSOMESQ,'ADDS',2,'SOMSCH')
      call memman(KSOMESCR,LSOMESQ**2,'ADDS',2,'SOMSCR')
CTF
CTF   dimension SOMESCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
CTF   dimension SOMEH(MXTSOB*MXTSOB)
      NTEST = 0
      IF(NTEST.GE.5)
     &WRITE(6,*) ' SBLOCK : ISSPC,ICSPC ', ISSPC,ICSPC
C?    WRITE(6,*) ' LUC in SBLOCK ', LUC
      IF(LUCBLK.GT.0) THEN
        CALL REWINE(LUCBLK,-1)
      END IF
*
*     info for this internal space
*     type of alpha and beta strings
      IATP = 1
      IBTP = 2
*     alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*     alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*     number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*     offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
*     connection matrices for supergroups
*
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
C     SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB),IPRCIX)
*     string sym, string sym => sx sym
*     string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*     largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*     largest number of strings of given symmetry and type
      MAXA = 0
      MAXA0 = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXA = MAX(MAXA,MAXA0)
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
*
      MAXB = 0
      MAXB0 = IMNMX(WORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MAXB = MAX(MAXB,MAXB0)
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.3 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*     largest number of resolution strings and spectator strings
*     that can be treated simultaneously
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
*     largest active orbital block belonging to given type and symmetry
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
        DO IOBSM = 1, NSMOB
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        END DO
      END DO
C?    WRITE(6,*) ' MXTSOB_L = ', MXTSOB_L
      MAXIJ = MXTSOB_L ** 2
*     local scratch arrays for blocks of C and sigma
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' ICISTR,LBLOCK ',ICISTR,LBLOCK
*     SCRATCH space for integrals
*     SCRATCH space for integrals
*     4 index integral block with four indeces belonging OS class
      INTSCR = MXTSOB_L ** 4
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. sigma needed for MXRESC
C          IAIBCM(ICISPC,IAIB)
      CALL IAIBCM(ISSPC,WORK(KSIOIO))
      CALL IAIBCM(ICSPC,WORK(KCIOIO))
*. Arrays giving block type
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      CALL MXRESCPH(WORK(KSIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSTFSMSPGP,MXPNSMST,
     &            NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &            NELFSPGP,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK,
     &            IPHGAS,NHLFSPGP,MNHL,IADVICE)
      IF(IPRCIX.GE.3) THEN
        WRITE(6,*) 'SBLOCK : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
         WRITE(6,*) 'SBLOCK : MXADKBLK ', MXADKBLK
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB)
      IF(IPRCIX.GE.3)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*
      IF(IPRCIX.GE.3)  WRITE(6,*) ' LSCR2 = ', LSCR2
C
      LSCR12 = MAX(LBLOCK,2*LSCR2)
csk?  LSCR12 = MIN(LBLOCK,2*LSCR2)
C     I assume memory was allocated for blocks, so
      KC2 = KVEC3
*
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*     space for active/passive blocks
      KCJPA = KC2 + 2*LSCR2
      KSIPA = KC2 + 3*LSCR2
*
*     vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*     I1 and Xi1s must also be able to hold largest st block
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0)
      CALL MEMMAN(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
*
      CALL MEMMAN(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
*
      CALL MEMMAN(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
*
      CALL MEMMAN(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*     some TTS arrays
      NOOS = NOCTPA*NOCTPB*NSMCI
      NTTS = MXNTTS
C
C     for partitioning of vector
C
      CALL MEMMAN(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      CALL MEMMAN(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      CALL MEMMAN(KLI1BT,NTTS  ,'ADDL  ',1,'I1BTC ')
      CALL MEMMAN(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
C
C     for scaling for each TTS block
C
      CALL MEMMAN(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')
C
C     space for four blocks of string occupations and arrays of
C     reordering arrays
C     
      LZSCR = (MAX(NAEL,NBEL)+3)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)+2) * NOCOB
      DO I1234 = 1, 1
        CALL MEMMAN(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      END DO

      DO I1234 = 1, 2
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*     Place perturbation integrals over one body integrals
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*
*
      CALL SBLOCKS(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KC2),
     &     WORK(KCIOIO),ISMOST(1,ICSM),WORK(KCBLTP),
     &     NACOB,WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &     NAEL,IATP,NBEL,IBTP,
     &     IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &     NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &     ITSOB,MAXIJ,MAXK,MAXI,INSCR,LBLOCK,
     &     LBLOCK,WORK(KINSCR),WORK(KCSCR),WORK(KSSCR),
     &     SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &     ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &     MXTSOB,WORK(KSOMESCR),WORK(KSOMEH),XIJILS,
     &     WORK(KI1),WORK(KXI1S),
     &     WORK(KI2),WORK(KXI2S),IDOH2,MXPOBS,WORK(KSVST),
     &     PSSIGN,IPRDIA,LUC,ICJKAIB,WORK(KCJRES),
     &     WORK(KSIRES),WORK(KI3),WORK(KXI3S),
     &     WORK(KI4),WORK(KXI4S),MXSXBL,
     &     MOCAA,MOCAB,IAPR,
     &     WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &     IRESTRICT,WORK(KCONSPA),WORK(KCONSPB),WORK(KLSCLFAC),
     &     LUCBLK,IPERTOP,IH0INSPC,WORK(KLH0SPC),
     &     ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,I_RES_AB,
     &     IUSE_PA,WORK(KCJPA),WORK(KSIPA),ISCLFAC_GROUP,SSCLFAC
#ifdef VAR_MPI
     &    ,LUCLIST,NPARBLOCK,IBLOCKL
#endif
     &                  )
*
      IF(IDC.EQ.2) THEN
*       reform
        CALL RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*       scale
        CALL SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK,NSMST,NOCTPA,NOCTPB,
     &              WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &              IDC,1,NTEST)
      END IF
*
      IF(LUCBLK.GT.0) THEN
        CALL ITODS(-1,1,-1,LUCBLK)
      END IF
*     eliminate local memory
      IDUM = 0
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'SBLOCK')
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE SBLOCKS(NSBLOCK,ISBLOCK,CB,SB,C2,
     &           ICOCOC,ICSMOS,ICBLTP,NACOB,NSSOA,NSSOB,
     &           NAEL,IAGRP,NBEL,IBGRP,
     &           IOCTPA,IOCTPB, NOCTPA,NOCTPB,
     &           NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
     &           ITSOB,MAXIJ,MAXK,MAXI,LI,LC,LS,
     &           XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &           SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &           MXTSOB,SOMESCR,SOMEH,XIJILS,
     &           I1,XI1S,I2,XI2S,IDOH2,MXPOBS,ISTRFL,
     &           PS,IPRNT,LUC,ICJKAIB,CJRES,SIRES,I3,XI3S,
     &           I4,XI4S,MXSXBL,MOCAA,MOCAB,IAPR,
     &           LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,IRESTRICT,
     &           ICONSPA,ICONSPB,SCLFAC,
     &           LUCBLK,IPERTOP,IH0INSPC,IH0SPC,
     &           ICBAT_RES,ICBAT_INI,ICBAT_END,IUSE_PH,IPHGAS,
     &           I_RES_AB,IUSE_PA,CJPA,SIPA,ISCLFAC_GROUP,SSCLFAC
#ifdef VAR_MPI
     &          ,LUCLIST,NPARBLOCK,IBLOCKL
#endif
     &           )
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen , Winter of 1991
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included
*   ISBLOCK(1,*) : alpha type of block
*   ISBLOCK(2,*) : beta type of block
*   ISBLOCK(3,*) : sym of alpha in block
*   ISBLOCK(4,*) : Offset of block
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
* ICJKAIB = 1 => construct C(Ka,Jb,j) and S(Ka,IB,i) as intermediate terms
*         = 0 => do not construct the above montioned matrices
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* If IRESTRICT.NE. 0 THEN we are after :
* sigma(iblk) = summa(jblk.le.iblk) (2-delta(iblk,jblk))/2
*                                                 * <Iblk!H!Jblk>C(Jblk)
      IMPLICIT REAL*8(A-H,O-Z)
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*. Specific input
      INTEGER ISBLOCK(8,*)
*.General input
      INTEGER ICOCOC(NOCTPA,NOCTPB)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA(NSMST ,*), NSSOB(NSMST ,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST), ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA,NOCTPA), ICONSPB(NOCTPB,NOCTPB)
*     scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),I3(*),XI1S(*),XI2S(*),XI3S(*)
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*     zero order Hamiltonian
      INTEGER IH0SPC(NOCTPA,NOCTPB)
      INTEGER IH0INSPC(*)
*     resolution matrices
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION SCLFAC(*),SSCLFAC(*)
#ifdef VAR_MPI
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCR
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR_C1, IOFF_SCR_C_EXPL
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR_C2
      DIMENSION LUCLIST(*),NPARBLOCK(*),IBLOCKL(*)
#endif
*     new scratch
      DOUBLE PRECISION CPUSIGXXX1,CPUSIGXXX2,WALLSIGXXX1,WALLSIGXXX2
      DOUBLE PRECISION CPUSIGXXXX1,CPUSIGXXXX2
      DOUBLE PRECISION WALLSIGXXXX1,WALLSIGXXXX2
      DOUBLE PRECISION CPUSIGXX1,CPUSIGXX2,WALLSIGXX1,WALLSIGXX2
      DOUBLE PRECISION CPUSIGX1,CPUSIGX2,WALLSIGX1,WALLSIGX2
      CHARACTER SECTID*12, CPUTIDXXX*12, WALLTIDXXX*12
      CHARACTER CPUTIDXXXX*12, WALLTIDXXXX*12
      CHARACTER CPUTIDX*12, WALLTIDX*12
      CHARACTER CPUTIDXX*12, WALLTIDXX*12
      DIMENSION ISCLFAC_GROUP(*)
#include "comjep.inc"
*
      CALL QENTER('SBLOC')
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.10) THEN
        WRITE(luwrt,*) ' ================='
        WRITE(luwrt,*) ' SBLOCKS speaking :'
        WRITE(luwrt,*) ' ================='
        WRITE(luwrt,*)
        WRITE(luwrt,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(luwrt,*) ' TTSS for each ACTIVE sigma block'
          DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0)
     &      WRITE(luwrt,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
          END DO
          WRITE(luwrt,*) ' IDC PS IPERTOP', IDC,PS,IPERTOP
          WRITE(luwrt,*) ' IDOH2 = ',IDOH2
          WRITE(luwrt,*) ' IRES_AB=',IRES_AB
      END IF
C     copy scaling factors from restore array to working array
      CALL DCOPY(NUM_BLOCKS2,SSCLFAC,1,SCLFAC,1)
      IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
C     initialize offsets in C-vector file
      IOFFSET_SCR     = 0
      IOFF_SCR_C1     = 0
      IOFF_SCR_C2     = 0
      IOFF_SCR_C_EXPL = 0
      INT_IOFF_C1     = 1
      INT_IOFF_C2     = 0
C
C     find batches of C-blocks; restrict to 'active' blocks
C     for a given CPU
C
      ITTSS_ORD = 2
      CALL PART_CIV_PAR2(IDC,ICBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,NSMST,
     &                   LBLOCK_LUCI,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,
     &                   LECBLOCK,I1CBLOCK,ICBLOCK,0,ITTSS_ORD,
     &                   SCLFAC,NUM_BLOCKS2,0)
#endif
      ELSE ! (LUCI_NMPROC .GT. 1) THEN
      IF(NTEST.GE.50) THEN
        WRITE(LUCIWT,*) ' Initial C vector '
        CALL WRTVCD(CB,LUC,1,-1)
      END IF
C
C     find batches of C-blocks
C
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,ICBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,NSMST,
     &              LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &              I1CBLOCK,ICBLOCK,0,ITTSS_ORD)
C
C     rewind C block file
C
      REWIND LUC
C
      END IF ! (LUCI_NMPROC .GT. 1) THEN
C     set max. excitations allowed by operator IDOH2
      IF(IDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE
        MXEXC = 1
      END IF
C?    WRITE(6,*) ' IDOH2, MXEXC',IDOH2,MXEXC
      IF(ICBAT_RES.EQ.1) THEN
        WRITE(6,*) ' Restricted set of C batches '
        WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
        JCBAT_INI = ICBAT_INI
        JCBAT_END = ICBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = NCBATCH
      END IF
C
C     loop over batches over C blocks
C     -------------------------------
C
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END
C
C       read C blocks into core
C
        ICOFF = 1
        NJBLOCK = LCBLOCK(JCBATCH)

#ifdef VAR_MPI
        IF (LUCI_NMPROC .GT. 1) THEN
C
C       C block non-zero?
C
        ICOMPUTE_C  = 0
C
        DO ICBLK = I1CBLOCK(JCBATCH), I1CBLOCK(JCBATCH)-1+NJBLOCK
          IF( ISCLFAC_GROUP( ICBLK ) .gt. 0 ) THEN
             ICOMPUTE_C = 1
csk       WRITE(LUCIWT,*) ' active C block found',ICBLK
          END IF
        END DO ! loop over blocks in this C-batch
        IF( ICOMPUTE_C .eq. 0 ) THEN
          DO ICBLK = I1CBLOCK(JCBATCH), I1CBLOCK(JCBATCH)-1+NJBLOCK
            LBL             = IBLOCKL( ICBLK )
            IOFF_SCR_C_EXPL = LBL
            IOFF_SCR_C1     = IOFF_SCR_C1 + IOFF_SCR_C_EXPL
            SCLFAC( ICBLK ) = 0.0D0
          END DO
          GOTO 20000
        END IF
        END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif
        DO JJCBLOCK = 1, NJBLOCK
          JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
          INTERACT = 0
!!!#ifdef VAR_MPI  -- for Stefan to check/remove later !!! 2-dec-2010
!!!          IF( ISCLFAC_GROUP( JBLOCK ) .gt. 0 ) THEN
!!!#else
          IF( SCLFAC(JBLOCK) .ne. 0.0D0 ) THEN
!!!#endif
csk       WRITE(LUCIWT,*) ' Nonvanishing block ', JBLOCK
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
*           loop over sigma blocks in batch
            DO JSBLOCK = 1, NSBLOCK
            IDENT = 0
C
C           parallel CI: ISBLOCK(1,JSBLOCK) == 0 for 
C           blocks which are not assigned to a given CPU
C
            IF( ISBLOCK(1,JSBLOCK) .gt. 0 ) THEN
              IATP = ISBLOCK(1,JSBLOCK)
              IBTP = ISBLOCK(2,JSBLOCK)
              IASM = ISBLOCK(3,JSBLOCK)
              IBSM = ISBLOCK(4,JSBLOCK)
*             are the two blocks connected by allowed excitation
              IF(MXEXC.EQ.1) THEN
                 IF((ICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (ICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((ICONSPA(IATP,LLATP).LE.1.AND.
     &               ICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (ICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (ICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
*             are they identical ?
              IDENT = 0
              IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.
     &           IBSM.EQ.JBSM.AND.IBTP.EQ.JBTP) IDENT = 1
*
            END IF
C           ^ IF( ISBLOCK(1,JSBLOCK) .gt. 0 ) THEN
            END DO
C           ^ s-block loop
          END DO
*         ^ IPERM loop
          END IF
*         ^ Checking was only done for nonvanishing blocks
*
          ISCALE = 0
        IF (LUCI_NMPROC .GT. 1) THEN
#ifdef VAR_MPI
C         length of block
          LBL         = IBLOCKL(JBLOCK)
          IOFF_SCR_C2 = LBL
          IF( INTERACT .eq. 1 ) THEN
C
C            file offset
C
csk          WRITE(LUWRT,*) 'CB read in IOFFSET',IOFF_SCR_C1
             xxxreadtime = MPI_WTIME()
             CALL MPI_FILE_READ_AT(LUC,IOFF_SCR_C1,CB(JOFF),
     &                             LBL,MPI_REAL8,my_STATUS,IERR)
csk          CALL WRTMATMN(CB(JOFF),1,LBL,1,LBL,LUWRT)
             xreadtimebi = xreadtimebi - xxxreadtime + MPI_WTIME()
             IBI_MULT_BL = IBI_MULT_BL + 1
             SCLFAC(JBLOCK) = 1.0D0
          ELSE
             SCLFAC(JBLOCK) = 0.0D0
          END IF
          IOFF_SCR_C1 = IOFF_SCR_C1 + IOFF_SCR_C2
csk       WRITE(LUWRT,*) ' set new IOFF_SCR_C2 for JBLOCK   ',
csk  &                     IOFF_SCR_C2, JBLOCK
#endif
        ELSE ! (LUCI_NMPROC .GT. 1) THEN
          IF(INTERACT.EQ.1) THEN
            CALL GSTTBL(C,CB(JOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                  NOCTPA,NOCTPB,NSSOA,NSSOB,PS,ICOOSC,IDC,
     &                  PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK))
csk      write(luciwt,*)'read in done for JBLOCK',JBLOCK
csk      CALL WRTMATMN(CB(JOFF),NSSOB(JBSM,JBTP),
csk  &                 NSSOA(JASM,JATP),NSSOB(JBSM,JBTP),
csk  &                 NSSOA(JASM,JATP),LUCIWT)
          ELSE
            CALL IFRMDS(LBL,-1,1,LUC)
            CALL SKPRCD2(LBL,-1,LUC)
            SCLFAC(JBLOCK) = 0.0D0
          END IF
        END IF ! (LUCI_NMPROC .GT. 1) THEN
C
          IF(IDENT.EQ.1 .AND. LUCBLK.GT.0)THEN
            write(6,*) ' block will be copied to LUCLBLK'
            NJA = NSSOA(JASM,JATP)
            NJB = NSSOB(JBSM,JBTP)
            IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM)THEN
*             block is unpacked, pack first
              IF(SCLFAC(JBLOCK).EQ.0.0D0) THEN
                CALL DZERO(C2,NJA*(NJA+1)/2)
              ELSE
                CALL TRIPK2(CB(ICOFF),C2,1,NJA,NJA,PS)
              END IF
              CALL TODSCP(C2,NJA*(NJA+1)/2,-1,LUCBLK)
            ELSE
              IF(SCLFAC(ICBLK).EQ.0.0D0) THEN
                CALL SETVEC(CB(ICOFF),ZERO,NJA*NJB)
              END IF
              CALL TODSCP(CB(ICOFF),NJA*NJB,-1,LUCBLK)
            END IF
          END IF
*
          IF(NTEST.GE.100) THEN
            IF(INTERACT.EQ.1) THEN
              WRITE(luwrt,*) ' TTSS for C block read in  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            ELSE
              WRITE(luwrt,*) ' TTSS for C block skipped  '
              CALL IWRTMA(ICBLOCK(1,JBLOCK),4,1,4,1)
            END IF
          END IF
*
        END DO
C       ^ loop over blocks in C-batch
C
C       loop again over blocks in C-batch -- contraction run
C
C                          sigma_n = Hb_n
C                          --------------
C
        DO 9000 ICBLK = I1CBLOCK(JCBATCH), I1CBLOCK(JCBATCH)-1+NJBLOCK
C
!!!#ifdef VAR_MPI
!!!            IF( ISCLFAC_GROUP( ICBLK ) .gt. 0 ) THEN
!!!#else
            IF( SCLFAC( ICBLK ) .ne. 0.0D0 ) THEN
!!!#endif
              JATP = ICBLOCK(1,ICBLK)
              JBTP = ICBLOCK(2,ICBLK)
              JASM = ICBLOCK(3,ICBLK)
              JBSM = ICBLOCK(4,ICBLK)
              ICOFF = ICBLOCK(5,ICBLK)
              NJA = NSSOA(JASM,JATP)
              NJB = NSSOB(JBSM,JBTP)
C             other symmetry blocks that can be obtained from this block
              CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                    LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
C             start with transposed block
              DO 8765 IPERM = NPERM,1, -1
                LLASM = LASM(IPERM)
                LLBSM = LBSM(IPERM)
                LLATP = LATP(IPERM)
                LLBTP = LBTP(IPERM)
                NLLA = NSSOA(LLASM,LLATP)
                NLLB = NSSOB(LLBSM,LLBTP)
C               the routines assumes on input that the blocks are transposed, 
C               so, initial tour, IPERM = 1 corresponds always to no transpose, 
C               so transpose!
                IF(IPERM.EQ.1) THEN
                  IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM) THEN
C                   diagonal blocks, Transposing corresponds to scaling
                    IF(PS.EQ.-1.0D0) THEN
                      CALL DSCAL(NJA*NJB,PS,CB(ICOFF),1)
                    END IF
                  ELSE
C                   offdiagonal blocks, explicit transposing
                    CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
                    CALL DCOPY(NJA*NJB,C2,1,CB(ICOFF),1)
                    if (LSGN(IPERM).eq.-1) then
                      write(6,*) 'LROW,LCOL undef. in sblocks'
                      write(6,*) 'and LSGN(IPERM).eq.-1'
                      Call Abend2('Quitting.')
CTF                   CALL SCALVE(CB(ICOFF),-1.0D0,LROW*LCOL)
                    end if
                  END IF
                END IF
C
C               loop over sigma blocks in this batch
C
C               NOTE: in parallel execution ISBLOCK(1,ISBLK) is 
C                     set to 0 if block is not included in block 
C                     distribution list NPARBLOCK for a given CPU.
C
                DO 10000 ISBLK = 1, NSBLOCK
                IF( ISBLOCK(1,ISBLK) .GT. 0 ) THEN
                  IATP = ISBLOCK(1,ISBLK)
                  IBTP = ISBLOCK(2,ISBLK)
                  IASM = ISBLOCK(3,ISBLK)
                  IBSM = ISBLOCK(4,ISBLK)
                  ISOFF = ISBLOCK(5,ISBLK)
                  NIA = NSSOA(IASM,IATP)
                  NIB = NSSOB(IBSM,IBTP)
C
                  IF(NIA*NIB.EQ.0) GOTO 10000
                  IF(IRESTRICT.EQ.1.AND.
     &               (JASM.GT.IASM.OR.
     &               JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
     &               JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
     &            GOTO 10000
C             are the two blocks connected by allowed excitation
              INTERACT = 0
              IF(MXEXC.EQ.1) THEN
                 IF((ICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (ICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((ICONSPA(IATP,LLATP).LE.1.AND.
     &               ICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (ICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (ICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
C
                 IF(INTERACT.EQ.0) GOTO 10000

                 IF(NTEST.GE.100) THEN
                   WRITE(luwrt,*) ' Next s block in batch : '
                   write(luwrt,*) ' ISBLK IASM IBSM IATP IBTP'
                   write(luwrt,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                 END IF
C
                IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &            ((LLBSM.GT.LLASM).OR.
     &            (LLASM.EQ.LLBSM).AND.(LLBTP.GT.LLATP)))
     &            GOTO 8764
C
                IF(NTEST.GE.60) THEN
                  WRITE(LUCIWT,*) ' RSSBCB will be called for '
                  WRITE(LUCIWT,*) ' Sigma block : '
                  WRITE(LUCIWT,*) ' ISOFF ', ISOFF
                  WRITE(LUCIWT,*) ' ISBLK IASM IBSM IATP IBTP'
                  WRITE(LUCIWT,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  WRITE(LUCIWT,*) ' C     block : '
                  WRITE(LUCIWT,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                  WRITE(LUCIWT,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                  WRITE(LUCIWT,*) ' ICOFF ', ICOFF
                  WRITE(LUCIWT,*) ' Overall scale',SCLFAC(ICBLK)
                END IF
*
                IF(IRESTRICT.EQ.1.AND.
     &             ((IASM.EQ.LLASM.AND.IBSM.EQ.LLBSM.AND.
     &               IATP.EQ.LLATP.AND.IBTP.EQ.LLBTP     ) .OR.
     &              (IDC.EQ.2.AND.
     &               IASM.EQ.LLBSM.AND.IBSM.EQ.LLASM.AND.
     &               IATP.EQ.LLBTP.AND.IBTP.EQ.LLATP     )     ))THEN
                     XFAC = 0.5D0*SCLFAC(ICBLK)
                ELSE
                     XFAC = SCLFAC(ICBLK)
                END IF
C               form of operator in action
                IF(IPERTOP.NE.0) THEN
C                  not exact Hamiltonian in use
                   IPTSPC = IH0SPC(IATP,IBTP)
                   JPTSPC = IH0SPC(JATP,JBTP)
                   IJOP   = IH0INSPC(IPTSPC)
*
C?                 WRITE(6,*) ' SBLOCKS:IPTSPC,JPTSPC,IJOP',
C?   &                         IPTSPC,JPTSPC,IJOP
*
                ELSE
                   IPTSPC = -1
                   JPTSPC = -1
                   IJOP   = -1
                END IF
*
#ifdef VAR_MPI
                IF (LUCI_NMPROC .GT. 1) THEN
                   xxxcomputtime = MPI_WTIME()
                END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif
*
                CALL RSSBCB2(IASM,IATP,IOCTPA,
     &               IBSM,IBTP,IOCTPB,
     &               LLASM,LLATP,LLBSM,LLBTP,NGAS,
     &               NELFSPGP(1,IATP+IOCTPA-1),
     &               NELFSPGP(1,IBTP+IOCTPB-1),
     &               NELFSPGP(1,LLATP+IOCTPA-1),
     &               NELFSPGP(1,LLBTP+IOCTPB-1),
     &               NAEL,NBEL,
     &               IAGRP,IBGRP,
     &               SB(ISOFF),CB(ICOFF),IDOH2,
     &               ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &               NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &               SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &               NSMOB,NSMST,NSMSX,NSMDX,
     &               NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
     &               ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &               MXSXBL,MOCAA,MOCAB,IAPR,
     &               IPRNT,IPERTOP,IPTSPC,JPTSPC,IJOP,0,IDUM,XFAC,
     &               IUSE_PH,IPHGAS,I_RES_AB,IUSE_PA,CJPA,SIPA,
     &               MXTSOB,SOMESCR,SOMEH,XIJILS,luwrt)
#ifdef VAR_MPI
                IF (LUCI_NMPROC .GT. 1) THEN
                   xcomputesi = xcomputesi - xxxcomputtime + MPI_WTIME()
                END IF ! (LUCI_NMPROC .GT. 1) THEN
#endif
C
 8764         CONTINUE
              END IF
10000         CONTINUE
C             ^ end of loop over sigma blocks
 8765         CONTINUE
            END IF
C           ^ non-vanishing C-block
 9000       CONTINUE
C           ^ end of loop over C blocks in Batch
20000   CONTINUE
C       ^ end of loop over batches of C blocks
C
C     order
C
      DO  ISBLK = 1 , NSBLOCK
      IF(ISBLOCK(1,ISBLK).GT.0) THEN
        IATP = ISBLOCK(1,ISBLK)
        IBTP = ISBLOCK(2,ISBLK)
        IASM = ISBLOCK(3,ISBLK)
        IBSM = ISBLOCK(4,ISBLK)
        ISOFF  = ISBLOCK(5,ISBLK)
        ISOFFP = ISBLOCK(6,ISBLK)
        NIA = NSSOA(IASM,IATP)
        NIB = NSSOB(IBSM,IBTP)
        IF(ICJKAIB.NE.0) THEN
C
C          tranposed sigma block was obtained, 
C          transpose to obtain correct block
C          ----------------------------------
C
           CALL TRPMT3(SB(ISOFF),NSSOB(IBSM,IBTP),
     &                 NSSOA(IASM,IATP),C2)
           CALL DCOPY(NSSOA(IASM,IATP)*NSSOB(IBSM,IBTP),C2,1,
     &                SB(ISOFF),1)
        END IF
        IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP)
     &  CALL TRPAD3(SB(ISOFF),PS,NSSOA(IASM,IATP))
*
      END IF
      END DO
*
      IF(NTEST.GE.50) THEN
        WRITE(luwrt,*) ' output blocks from SBLOCKS '
        CALL WRTTTS(SB,ISBLOCK,NSBLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSSOA,NSSOB,1)
      END IF
C
      CALL QEXIT('SBLOC')
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE SKICKJ(SKII,CKJJ,NKA,NIB,NJB,NKB,XIJKL,
     &                  NI,NJ,NK,NL,MAXK,
     &                  KBIB,XKBIB,KBJB,XKBJB,IKORD,
     &                  FACS,IROUTE,MXTSOB,XIJILS)
*
*
* Calculate S(Ka,Ib,i) = FACS*S(Ka,Ib,i)
*          +SUM(j,k,l,Kb) <Ib!a+ kb!Kb><Kb!a lb !Jb>*(ij!kl)*C(Ka,Jb,j)
*
*
*
* Jeppe Olsen, Spring of 94
*
* : Note : Route 1 has retired, March 97
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION CKJJ(NKA,*)
      DIMENSION XIJKL(*)
      DIMENSION KBIB(MAXK,*),XKBIB(MAXK,*)
      DIMENSION KBJB(MAXK,*),XKBJB(MAXK,*)
*. Input and output
      DIMENSION SKII(NKA,*)
*. Scratch
CTF   PARAMETER(MXTSOB=40)
C     COMMON/SOMESCR/SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
C     DIMENSION IBOFF(MXTSOB*MXTSOB),JBOFF(MXTSOB*MXTSOB)
C     DIMENSION SCR(NI*NL*NJ),IBOFF(NI),JBOFF(NL*NJ)
      DIMENSION XIJILS(MXTSOB)
      IF(NI.GT.MXTSOB.OR.NJ.GT.MXTSOB.OR.NK.GT.MXTSOB
     &   .OR.NL.GT.MXTSOB) THEN
         WRITE(6,*) ' SKICKJ : Too many orbs : NI > MXTSOB '
         WRITE(6,*) ' NI, MXTSOB ',MAX(NI,NJ,NK,NL),MXTSOB
         Call Abend2( ' Redim MXTSOB in SKICKJ' )
      END IF
*
C     CALL QENTER('SKICK')
      IF(IROUTE.EQ.3) THEN
* S(Ka,i,Ib) = S(Ka,i,Ib) + sum(j) (ji!kl) C(Ka,j,Jb)
        DO KB = 1, NKB
*. Number of nonvanishing connections from KB
         LL = 0
         KK = 0
         DO L = 1, NL
           IF(KBJB(KB,L).NE.0) LL = LL + 1
         END DO
         DO K = 1, NK
           IF(KBIB(KB,K).NE.0) KK = KK + 1
         END DO
*
         IF(KK.NE.0.AND.LL.NE.0) THEN
           DO K = 1, NK
             IB = KBIB(KB,K)
             IF(IB.NE.0) THEN
               SGNK = XKBIB(KB,K)
               DO L = 1, NL
                 JB = KBJB(KB,L)
                 IF(JB.NE.0) THEN
                   SGNL = XKBJB(KB,L)
                   FACTOR = SGNK*SGNL
*. We have now a IB and Jb string, let's do it
                   ISOFF = (IB-1)*NI*NKA + 1
                   ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
                   IMAX = NI
*
                   IF(IKORD.NE.0) THEN
*. Restrict so (ij) .le. (kl)
                     IMAX  = K
                     JKINTOF = INTOF + (K-1)*NJ
C                    CALL COPVEC(XIJKL(JKINTOF),XIJILS,NJ)
                     DO J = L,NL
                       XIJILS(J) = XIJKL(JKINTOF-1+J)
                     END DO
                     XIJKL(JKINTOF-1+L) = 0.5D0*XIJKL(JKINTOF-1+L)
                     DO J = L+1, NL
                      XIJKL(JKINTOF-1+J) = 0.0D0
                     END DO
                   END IF
C                  ONE = 1.0D0
                   CALL MATML7(SKII(ISOFF,1),  CKJJ(ICOFF,1),
     &                         XIJKL(INTOF),NKA,IMAX,NKA,NJ,
     &                         NJ,IMAX,FACS,FACTOR ,0)
                   IF(IKORD.NE.0) THEN
                      DO J = L,NL
                        XIJKL(JKINTOF-1+J) =  XIJILS(J)
                      END DO
C                    CALL COPVEC(XIJILS,XIJKL(JKINTOF),NJ)
                   END IF
*
                 END IF
               END DO
*
             END IF
           END DO
         END IF
       END DO
*. (end over loop over Kb strings )
      ELSE IF(IROUTE.EQ.2) THEN
* S(I,Ka,Ib) = S(I,Ka,Ib) + sum(j) (ij!kl) C(j,Ka,Jb)
        DO KB = 1, NKB
*. Number of nonvanishing connections from KB
         LL = 0
         KK = 0
         DO L = 1, NL
           IF(KBJB(KB,L).NE.0) LL = LL + 1
         END DO
         DO K = 1, NK
           IF(KBIB(KB,K).NE.0) KK = KK + 1
         END DO
*
         IF(KK.NE.0.AND.LL.NE.0) THEN
           DO K = 1, NK
             IB = KBIB(KB,K)
             IF(IB.NE.0) THEN
               SGNK = XKBIB(KB,K)
               DO L = 1, NL
                 JB = KBJB(KB,L)
                 IF(JB.NE.0) THEN
                   SGNL = XKBJB(KB,L)
                   FACTOR = SGNK*SGNL
*. We have now a IB and Jb string, let's do it
                   ISOFF = (IB-1)*NI*NKA + 1
                   ICOFF = (JB-1)*NJ*NKA + 1
                   INTOF = ((L-1)*NK + K - 1 )*NI*NJ + 1
*
                   JMAX = NJ
                   IF(IKORD.NE.0) THEN
*. Restrict so (ji) .le. (kl)
                     JMAX  = K
                     IKINTOF = INTOF + (K-1)*NI
                     CALL COPVEC(XIJKL(IKINTOF),XIJILS,NI)
                     XIJKL(IKINTOF-1+L) = 0.5D0*XIJKL(IKINTOF-1+L)
                     DO I = L+1, NL
                      XIJKL(IKINTOF-1+I) = 0.0D0
                     END DO
                   END IF
*
C                  ONE = 1.0D0
                   CALL MATML7(SKII(ISOFF,1), XIJKL(INTOF),
     &                         CKJJ(ICOFF,1),NI,NKA,NI,NJ,
     &                         NJ,NKA,FACS,FACTOR,0)
*
                 IF(IKORD.NE.0) THEN
                   CALL COPVEC(XIJILS,XIJKL(IKINTOF),NI)
                 END IF
*
                 END IF
               END DO
             END IF
           END DO
         END IF
       END DO
*. (end over loop over Kb strings )
*

      ELSE IF (IROUTE.EQ.1) THEN
         WRITE(6,*) ' Sorry route 1 has retired, March 1997'
         Call Abend2( 'SKICKJ:Invalid route=1' )
C     DO 1000 KB = 1, NKB
*. Number of nonvanishing a+lb !Kb>
C       LL = 0
C       DO L = 1, NL
C         IF(KBJB(KB,L).NE.0) LL = LL + 1
C       END DO
*
C       IKEFF = 0
C       DO 900 K = 1, NK
C         IB = KBIB(KB,K)
C         IF(IB.EQ.0) GOTO 900
C         SGNK = XKBIB(KB,K)
*
C         IF(IKORD.EQ.0) THEN
C            LI = NI
C            IMIN = 1
C         ELSE
C            LI = NI-K+1
C            IMIN = K
C         END IF
*
C         DO 700 I = IMIN, NI
C           IKEFF = IKEFF + 1
C           IOFF = (IKEFF-1)*NJ*LL
*. Offset for S(1,IB,i)
C           IBOFF(IKEFF)  = (I-1)*NIB+IB
C           LEFF = 0
C           DO 800 L = 1, NL
C             JB = KBJB(KB,L)
C             IF(JB.EQ.0) GOTO 800
C             LEFF = LEFF + 1
C             SGNL = XKBJB(KB,L)
C             IF(IKORD.EQ.1.AND.I.EQ.K)THEN
C                FACTOR = 0.5D0*SGNK*SGNL
C             ELSE
C                FACTOR =       SGNK*SGNL
C             END IF
C             JL0 = (LEFF-1)*NJ
C             JLIK0 = (K-1)*NJ*NL*NI
C    &              + (I-1)*NJ*NL
C    &              + (L-1)*NJ
C             DO 600 J = 1, NJ
C               JL = JL0 + J
*. Offsets for C(1,JB,j)
C               JBOFF(JL) = (J-1)*NJB + JB
*. integral * signs in SCR(jl,ik)
*. Integrals are stored as (j l i k )
C               SCR((IKEFF-1)*NJ*LL+JL) = FACTOR*XIJKL(JLIK)
C               SCR(IOFF+JL) = FACTOR*XIJKL(JLIK0+J)
C 600         CONTINUE
C 800       CONTINUE
C 700     CONTINUE
C 900   CONTINUE
*
C       CALL GSAXPY(SKII,CKJJ,SCR,IKEFF,NJ*LL,NKA,IBOFF,JBOFF)
C1000 CONTINUE
      END IF
*. End of IROUTE branchning
*
C     CALL QEXIT('SKICK')
      END

!**********************************************************************
      FUNCTION I_DX_ACT(IGAS,KGAS,LGAS,JGAS)
!**********************************************************************
*
* Is double excitation a+ igas a+ kgas a lgas a jgas active
*
      IMPLICIT REAL*8(A-H,O-Z)

*. Normal CI, extensions for perturbation etc can be inseted here
      I_DX_ACT = 1
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(I_DX_ACT.EQ.1) THEN
          WRITE(6,*)
     &    ' allowed excitation a+(igas) a (jgas) a+(kgas) a (lgas) for'
     &    ,IGAS,JGAS,KGAS,LGAS
        ELSE IF(I_DX_ACT.EQ.0) THEN
          WRITE(6,*)
     &    ' excluded excitation a+(igas) a (jgas) a+(kgas) a (lgas) for'
     &    ,IGAS,JGAS,KGAS,LGAS
        END IF
      END IF
      END

!**********************************************************************
      FUNCTION I_SX_ACT(IGAS,JGAS)
!**********************************************************************
*
* Is single excitation a+ igas a jgas active
*
* Hiding cc restrictions etc
*
      IMPLICIT REAL*8(A-H,O-Z)
*
*. Normal CI, extensions for perturbation etc can be inseted here
      I_SX_ACT = 1
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(I_SX_ACT.EQ.1) THEN
          WRITE(6,*)
     &    ' allowed excitation a+ (igas) a (jgas) for igas,jgas=',
     &    IGAS,JGAS
        ELSE IF(I_SX_ACT.EQ.0) THEN
          WRITE(6,*)
     &    ' exluded excitation a+ (igas) a (jgas) for igas,jgas=',
     &    IGAS,JGAS
        END IF
      END IF
      END
!**********************************************************************
      SUBROUTINE MATML7(C,A,B,NCROW,NCCOL,NAROW,NACOL,
     &                  NBROW,NBCOL,FACTORC,FACTORAB,ITRNSP )
C
C MULTIPLY A AND B TO GIVE C
C
C     C =  FACTORC*C + FACTORAB* A * B             FOR ITRNSP = 0
C
C     C =  FACTORC*C + FACTORAB* A(TRANSPOSED) * B FOR ITRNSP = 1
C
C     C =  FACTORC*C + FACTORAB* A * B(TRANSPOSED) FOR ITRNSP = 2
C
C... JEPPE OLSEN,
C
*. Notice : If the summation index has dimension zero nothing
*           is performed
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION A(NAROW,NACOL),B(NBROW,NBCOL)
      DIMENSION C(NCROW,NCCOL)
      COMMON/MATMLST/XNFLOP,XNCALL,XLCROW,XLCCOL,XLCROWCOL,TMULT
*
#ifdef LUCI_DEBUG
      NTEST = 00
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*) ' Initial A, B and C MATRIX FROM MATML7 ' 
      WRITE(6,*)      ' NCROW NCCOL NAROW NACOL NBROW NBCOL ' 
      WRITE(6,'(6I6)')  NCROW,NCCOL,NAROW,NACOL,NBROW,NBCOL
      WRITE(6,*) ' ITRNSP : ', ITRNSP
        IF (NCROW .lt. 0 .or. NCCOL .lt. 0 .or. NAROW .lt. 0 .or.
     &      NACOL .lt. 0 .or. NBROW .lt. 0 .or. NBCOL .lt. 0) THEN
           WRITE(6,'(/A)') 'FATAL ERROR in MATML7, neagtive dimensions'
           call our_own_traceback
        END IF
      CALL WRTMT_LU(A,NAROW,NACOL,NAROW,NACOL)
      CALL WRTMT_LU(B,NBROW,NBCOL,NBROW,NBCOL)
      CALL WRTMT_LU(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
#endif
*. Statistics
      XNCALL = XNCALL + 1
      XLCROW = XLCROW + NCROW
      XLCCOL = XLCCOL + NCCOL
      XLCROWCOL = XLCROWCOL + NCROW*NCCOL
      T_INI = Second()
*
      IF(ITRNSP.EQ.1) THEN
        XNFLOP = XNFLOP + 2*NCROW*NCCOL*NAROW
      ELSE
        XNFLOP = XNFLOP + 2*NCROW*NCCOL*NACOL
      END IF
C
C     do matrix-matrix multiplication using DGEMM
C
C     input
      LDA = MAX(1,NAROW)
      LDB = MAX(1,NBROW)
C     output
      LDC = MAX(1,NCROW)
      IF(ITRNSP.EQ.0) THEN
        CALL DGEMM('N','N',NAROW,NBCOL,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
      ELSE IF (ITRNSP.EQ.1) THEN
        CALL DGEMM('T','N',NACOL,NBCOL,NAROW,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
      ELSE IF(ITRNSP.EQ.2) THEN
        CALL DGEMM('N','T',NAROW,NBROW,NACOL,FACTORAB,A,LDA,
     &                 B,LDB,FACTORC,C,LDC)
      END IF
C
#ifdef LUCI_DEBUG
      IF ( NTEST .NE. 0 ) THEN
      WRITE(6,*) ' Final C MATRIX FROM MATML7 ' 
      CALL WRTMT_LU(C,NCROW,NCCOL,NCROW,NCCOL)
      END IF
#endif
C
      T_END = Second()
      TMULT = TMULT + T_END - T_INI 
C
      END
