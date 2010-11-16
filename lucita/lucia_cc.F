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
*fordeck lucia_cc

* single- and double- amplitudes
*
* =================
* Singles (E_{IJ})
* =================
*
* form 1 : Only symmetry correct terms included
*          organization :
*          Loop over symmetry of I
*          Find symmetry of J
*           Loop over type of I (a+)
*            Loop over type of J ( a)
*             IF (a+i aj is allowed ) :
*              Loop over J
*               Loop over I
*               End of loop over I
*              End of loop over J
*             End if ( allowed types )
*            End of loop over type of J
*           End of loop over type of I
*          End of loop over symmetry of I
*
* Form 2 : in complete matrix form IJ = (J-1)*NORB + I,
*          where I and J are the absolute orbital numbers
*          and the orbitals are type ordered
*
*
* ===================
* Doubles E(IJ) E(KL)
* ===================
*
* Form 1 : Nonredundant form without singlet/triplet splitting
*
* Loop over symmetry of I
*  Loop over symmetry of J
*   Loop over symmetry of K
*     => Symmetry of L
*    Loop over type of I
*     Loop over type of J
*      Loop over type of K
*       Loop over type of L
*        If (type is allowed ) then
*         Loop over L
*          Loop over K
*           Loop over J
*            Loop over I
*            End of loop over I
*           End of Loop over J
*          End of Loop over K
*         End of Loop over L
*        End if (types are allowed )
*       End of loop over types of L
*      End of loop over types of K
*     End of loop over type over J
*    End of loop over types of I
*
* Form 2 : Nonredundant form with singlet/triplet splitting
*
*   A singlet-singlet coupled
*     Loop structure as above
*   B triplet-triplet coupled
*     Loop structure as above
*
* Form 3 : redundant form in symmetric matrix (IJ,KL),IJ.GT.KL
*          and the orbitals are type ordered
*          ( this form will probably disappear later but ... )

*
      SUBROUTINE BIO_TO_STANDARD(V_BIO,V_STANDARD,IWAY)
*
* Transform s vector between standard format and biorthgonal
*
* Jeppe Olsen, Sept. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "lucinp.inc"
#include "orbinp.inc"
*. Specific input/output
      DIMENSION V_BIO(*),V_STANDARD(*)
*. For the 2s He case ( STANDARD => BIO)
      V_BIO(1) = 0.5*V_STANDARD(1)
      V_BIO(2) = 0.25*V_STANDARD(2)
*
      RETURN
      END
C     CALL GET_DIAG_BLMAT(WORK(KFI),WORK(KFDIA),NSMOB,NTOOB,1)
      SUBROUTINE GET_DIAG_BLMAT(A,DIAG,NBLK,LBLK,ISYM)
*
* Extract diagonal from blocked matrix
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION A(*)
      INTEGER LBLK(*)
*. Output
      DIMENSION DIAG(*)
*
      IOFF1 = 1
      IOFF2 = 1
C?    WRITE(6,*) '  GET_DIA..  NBLK = ', NBLK
      DO IBLK = 1, NBLK
        L = LBLK(IBLK)
C?      WRITE(6,*) ' IBLK and L ', IBLK,L
        CALL COPDIA(A(IOFF2),DIAG(IOFF1),L,ISYM)
*
        IOFF1 = IOFF1 + L
        IF(ISYM.EQ.1) THEN
          IOFF2 = IOFF2 + L*(L+1)/2
        ELSE
          IOFF2 = IOFF2 + L ** 2
        END IF
*
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Diagonal elements of blocked matrix '
        LEN = IOFF1 - 1
        CALL WRTMAT(DIAG,1,LEN,1,LEN)
      END IF
*
      RETURN
      END
      SUBROUTINE GET_H0_DIAG(DIAG,I_SYM,F,ISCALE)
*
* Obtain SCF H0 matrix  :
*
*           <HF![F,E]!HF>
*           <HF![F,EE]!HF>
* If Iscale .ne. 0, the diagonal is scaled as
*
* diag(a,i,b,j) => 1/2*(2-delta(ai,bj)) diag(a,i,b,j)
*
* Jeppe Olsen, Sept. 98
*
* Blocks are packed column wise
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "lucinp.inc"
*. Input : Diagonal elements of F
      DIMENSION F(*)
*. Output
      DIMENSION DIAG(*)
*
* One-electron part
*
      IJ = 0
      DO ISM =1, NSMOB
        JSM = MULTD2H(ISM,I_SYM)
        DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
          IJACT = I_SX_CCACT(IGAS,JGAS)
          IF(IJACT.EQ.1) THEN
*. Offsets
            I_OFF = IOBPTS(IGAS,ISM)
            J_OFF = IOBPTS(JGAS,JSM)
*. Numbers
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
*
            DO J = J_OFF,J_OFF+NJ-1
            DO I = I_OFF,I_OFF+NI-1
              IJ = IJ + 1
              DIAG(IJ) = F(I) - F(J)
C?            write(6,*) ' I J IJ DIAG(IJ) ', I,J,IJ,DIAG(IJ)
            END DO
            END DO
          END IF
        END DO
        END DO
      END DO
*
* Two-electron part
*
      IJKL = IJ
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,I_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
C?        WRITE(6,*) ' IGAS,JGAS,KGAS,LGAS,IJKL_ACT',
C?   &                 IGAS,JGAS,KGAS,LGAS,IJKL_ACT
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0
            ELSE
              IJBL_GT_KLBL = -1
            END IF
          END IF
C?        WRITE(6,*) ' IJBL_GT_KLBL' ,  IJBL_GT_KLBL
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            I_OFF = IOBPTS(IGAS,ISM)
*
            NJ = NOBPTS(JGAS,JSM)
            J_OFF = IOBPTS(JGAS,JSM)
*
            NK = NOBPTS(KGAS,KSM)
            K_OFF = IOBPTS(KGAS,KSM)
*
            NL = NOBPTS(LGAS,LSM)
            L_OFF = IOBPTS(LGAS,LSM)
*
            DO L = L_OFF,L_OFF+NL-1
            DO K = K_OFF,K_OFF+NK-1
C
            IF(IJBL_GT_KLBL .EQ. 0 ) THEN
             J_MIN = L
            ELSE
             J_MIN = J_OFF
            END IF
            DO J = J_MIN,J_OFF+NJ-1
*
            IF( IJBL_GT_KLBL .EQ. 1 ) THEN
              I_MIN = I_OFF
            ELSE IF ( IJBL_GT_KLBL .EQ. 0 ) THEN
              IF(J.GT.L) THEN
                I_MIN = I_OFF
              ELSE
                I_MIN = K
              END IF
            END IF
            DO I = I_MIN,I_OFF+NI-1
C
*
              IJKL = IJKL + 1
C?            WRITE(6,*) 'I,J,K,L,IJKL',
C?   &                    I,J,K,L,IJKL
*. I and K corresponds to creation ops, J and L to annihilation ops
              DIAG(IJKL) = F(I) + F(K) - F(J) - F(L)
*. Scale
              IF(ISCALE.EQ.1) THEN
                IF(I.EQ.K.AND.J.EQ.L) DIAG(IJKL) = 0.5D0*DIAG(IJKL)
              END IF
            END DO
            END DO
            END DO
            END DO
*           ^ End of loop over orbitals over given TS
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        IF(ISCALE.EQ.0) THEN
          WRITE(6,*) ' Unscaled <HF![F,t]!HF> matrix '
        ELSE
          WRITE(6,*) '   scaled <HF![F,t]!HF> matrix '
        END IF
        CALL WRT_CC_VEC(DIAG,LU)
      END IF
*
      RETURN
      END
      SUBROUTINE OPTIM_CC(CC_AMP,VEC1,VEC2)
*
* Optimize CC wave function
*
* Initial version,Sept 98
*
* In this initial version a simple sequence of Fock type iterations
* are used
*
* Jeppe Olsen
*
* Input :
* ========
*
* CC_AMP : Initial set of CC amplitiudes
*
* Output
* ========
*
* CC_AMP : Final set of amplitudes
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "cecore.inc"
#include "clunit.inc"
*. Scratch : 2 vector blocks for CI
      DIMENSION VEC1(*),VEC2(*)
*. Input and output
      DIMENSION CC_AMP(*)
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      NTEST = 05
*
      IDUMMY = 0
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ',IDUMMY,'CC_OPT')
*. Coupled cluster flag
      ICC_EXC = 1
*
*. (NSXE and NDXEE are stored in /CRUN/
*. Allocate space for additional vectors
      WRITE(6,*) ' CC_OPTIM : NSXE NDXEE N_CC_AMP', NSXE,NDXEE,N_CC_AMP
      LEN = N_CC_AMP
      CALL MEMMAN(KCC4,2*LEN,'ADDL  ',2,'CC4  ')
      CALL MEMMAN(KCC5,2*LEN,'ADDL  ',2,'CC5  ')
      CALL MEMMAN(KCC6,2*LEN,'ADDL  ',2,'CC6  ')
      CALL MEMMAN(KCC7,2*LEN,'ADDL  ',2,'CC7  ')
      CALL MEMMAN(KDIA,2*LEN,'ADDL  ',2,'CCDIA')
      CALL MEMMAN(KFDIA,NTOOB,'ADDL  ',2,'FDIA ')
*
* Obtain Fock matrix
*
*. One-body density matrix
      LBLK = -1
      CALL COPVCD(LUC,LUSC1,VEC1,1,LBLK)
      CALL DENSI2(1,WORK(KRHO1),WORK(KRHO2),VEC1,VEC2,
     &             LUSC1,LUC,EXPS2)
*. And Fock matrix
      CALL COPVEC(WORK(KINT1O),WORK(KFI),NINT1)
      CALL FIFAM(WORK(KFI))
      CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
      ECORE_H = 0.0D0
*. Extract diagonal
C      CALL GT1DIS(H1DIA,IREOTS(1),WORK(KPINT1),WORK(KINT1O),
C    &             ISMFTO,IBSO,NACOB)
      CALL GT1DIS(WORK(KFDIA),IREOTS,WORK(KPINT1),
     &            WORK(KFI),ISMFTO,IBSO,NACOB)
C     CALL GET_DIAG_BLMAT(WORK(KFI),WORK(KFDIA),NSMOB,NTOOBS,1)
*. Construct scaled diagonal
      ISCALE = 1
      CALL GET_H0_DIAG(WORK(KDIA),1,WORK(KFDIA),ISCALE)
*. (Maxit is obtained from CRUN)
      DO ITER = 1, MAXIT
        IF(NTEST.GE.5) THEN
           WRITE(6,*)
           WRITE(6,*) ' =================================='
           WRITE(6,*) ' Information from iteration ', ITER
           WRITE(6,*) ' =================================='
           WRITE(6,*)
        END IF
*. Calculate CC_vector function for current set of parameters
        IBIO = 1
        CALL CC_VEC_FNC(CC_AMP,WORK(KCC4),ECC1,VEC1,VEC2,IBIO)
        ECC = ECC1 + ECORE
        XNORM = SQRT(INPROD(WORK(KCC4),WORK(KCC4),N_CC_AMP))
*
        IF(NTEST.GE.5) WRITE(6,'(A,I5,2F25.12)')
     &  ' It, Energy, Vecnorm',ITER,ECC,XNORM
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' vector function : '
          CALL WRT_CC_VEC(WORK(KCC4),LU)
* (LU is not active)
        END IF
*. Transfer to Biorthogonal rep.( Done in CC_VEC_FNC)
C     BIO_TO_STANDARD(V_BIO,V_STANDARD,IWAY)
C       CALL BIO_TO_STANDARD(WORK(KCC5),WORK(KCC4),2)
C       CALL COPVEC(WORK(KCC5),WORK(KCC4),N_CC_AMP)
*. find change of CC coefficients as - Scaled(diag)-1 * CC_VEC_FNC
C            DIAVC2(VECOUT,VECIN,DIAG,SHIFT,NDIM)
        ZERO = 0.0D0
        CALL DIAVC2(WORK(KCC5),WORK(KCC4),WORK(KDIA),ZERO,N_CC_AMP)
        ONEM = -1.0D0
        CALL SCALVE(WORK(KCC5),ONEM,N_CC_AMP)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Correction  vector '
          CALL WRT_CC_VEC(WORK(KCC5),LU)
* (LU is not active)
        END IF
*. And add to current set of amplitudes
        ONE = 1.0D0
        CALL VECSUM(CC_AMP,CC_AMP,WORK(KCC5),ONE,ONE,N_CC_AMP)
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' Updated amplitudes '
          CALL WRT_CC_VEC(CC_AMP,LU)
* (LU is not active)
        END IF
*. Dump to LU_CCAMP
        CALL REWINO(LU_CCAMP)
        DO I = 1, N_CC_AMP
          WRITE(LU_CCAMP,'(E25.15)') CC_AMP(I)
        END DO
*
      END DO
*     ^ End of loop over iterations
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Final set of CC_amplitudes '
        CALL WRT_CC_VEC(CC_AMP,LU)
      END IF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'CC_OPT')
      RETURN
      END
      SUBROUTINE WRT_CC_VEC(CC,LU)
*
* Print vector of CC amplitudes
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Generel input
#include "mxpdim.inc"
#include "crun.inc"
      DIMENSION CC(*)
*
      WRITE(6,*)
      WRITE(6,*) ' ====================== '
      WRITE(6,*) ' Single excitation part '
      WRITE(6,*) ' ====================== '
      WRITE(6,*)
      CALL WRT_SX(CC,1)
*
      WRITE(6,*)
      WRITE(6,*) ' ====================== '
      WRITE(6,*) ' Double excitation part '
      WRITE(6,*) ' ====================== '
      WRITE(6,*)
      CALL WRT_DX1(CC(1+NSXE),1)
*
      RETURN
      END
      SUBROUTINE CC_VEC_FNC(CC_AMP,CC_VEC,E_CC,VEC1,VEC2,IBIO)
*
* Calculate energy and CC vector function for
* a set of CC amplitudes defined by CC_AMP
*
* Jeppe Olsen, Summer of 98
*
* Initial very dirty and slow version, allows however open-shell cc
* as well as multireference cc expansions
*
* The CC vector function reads
*
*     <\my! Exp(-T) H Exp (T) !0>
*
* where my is one of the following
*
*  Ibio = 0 :
*
* <0!E(ij)
* <0!e(ijkl)
* <K!
*
* IBIO = 1 :
*
* 1/2<0!E(ij)
* 1/3<0!e(ijkl)+1/6<0!e(ijkl)
* <K!
*
* With IBIO = 1, we thus transform to the biorthonormal basis

* Where E(ij) and e(ijkl) are the non-redundant orbital operators,
* and |K> is an orthogonal complement space for the reference space
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      DIMENSION CC_AMP(*)
*. Output
      DIMENSION CC_VEC(*)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*.
#include "mxpdim.inc"
#include "cc_exc.inc"
#include "clunit.inc"
#include "glbbas.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "cprnt.inc"
*
      REAL*8 INPRDD
*
      LBLK = -1
*
* ==========================
* .1 : Exp(T) !0> (on LUHC)
* ==========================
*
* Zero order state is assumed on LUC, results will be stored on LUHC
* CC amplitudes are assumed located in KCC1
       MX_TERM = 100
       ICC_EXC = 1
       THRES_E = 1.0D-12
       CALL EXPT_REF(LUC,LUHC,LUSC1,LUSC2,LUSC3,THRES_E,MX_TERM,
     &             VEC1,VEC2)
*
* ==========================
*. 2 H Exp(T) !0> (On LUSC1)
* ==========================
*
       ICC_EXC = 0
C?     WRITE(6,*) ' Input file to MV7 '
C?     CALL WRTVCD(VEC1,LUHC,1,LBLK)
       CALL MV7(VEC1,VEC2,LUHC,LUSC1)
*
* =================================
* 3 Exp(-T) H Exp(T) !0> (On LUSC2)
* =================================
*
       ONEM = -1.0D0
       ICC_EXC = 1
C?     WRITE(6,*) ' N_CC_AMP = ', N_CC_AMP
       CALL SCALVE(WORK(KCC1),ONEM,N_CC_AMP)
       CALL EXPT_REF(LUSC1,LUSC2,LUSC3,LUSC34,LUSC35,THRES_E,MX_TERM,
     &             VEC1,VEC2)
       CALL SCALVE(WORK(KCC1),ONEM,N_CC_AMP)
*
* ========================================================
* 4 <0!E Exp(-T) H Exp(T) !0>, <0!EE Exp(-T) H Exp(T) !0>
* ========================================================
*
*. Densities with !0> and  Exp(-T) H Exp(T) !0>
*           DENSI2(I12,RHO1,RHO2,L,R,LUL,LUR,EXPS2)
       IPRDEN_ORIG = IPRDEN
       IPRDEN = 0
       CALL DENSI2(2,WORK(KRHO1),WORK(KRHO2),VEC1,VEC2,
     &             LUSC2,LUC,EXPS2)
       IPRDEN = IPRDEN_ORIG
*
* ================================
* 5 Reformat to CC form
* ================================
*
      CALL REF_SX(WORK(KRHO1),CC_VEC(1),2,1,1,IBIO)
C          REF_SX(CIN,COUT,INFRM,IOUTFRM,ISX_SYM)
      CALL REF_DX_EXP_COMP(WORK(KRHO2),CC_VEC(1+NSXE),1,1,IBIO)
C          REF_DX_EXP_COMP(CEXP,CCOM,IDX_SYM,IWAY)
* ===============
* 6 : CC energy
* ===============
*
      LBLK = -1
      E_CC = INPRDD(VEC1,VEC2,LUC,LUSC2,1,LBLK)
C?    WRITE(6,'(A,F25.12)') ' Coupled Cluster Energy ', E_CC
*
*. Mission over, print
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' CC-vector function '
        WRITE(6,*) ' ==================='
        WRITE(6,*)
        CALL WRT_CC_VEC(CC_VEC,6)
      END IF
*.
      RETURN
      END

      SUBROUTINE GET_SX_BLK(HBLK,H,IGAS,ISM,JGAS,JSM)
*
* Fetch block of one-electron excitations from H
*
* Jeppe Olsen, August 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "orbinp.inc"
*. Input
      DIMENSION H(*)
*. Output
      DIMENSION HBLK(*)
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
*
      IJ_ACT = I_SX_ACT(IGAS,JGAS)
      IF(IJ_ACT.EQ.0) THEN
*. Just another empty block
       ZERO = 0.0D0
       CALL SETVEC(HBLK,ZERO,NI*NJ)
      ELSE
*. Block assumed total symmetric
       ISX_SYM = 1
       CALL I_OFF_SX(IOFF,IGAS,ISM,JGAS,JSM,ISX_SYM)
       CALL COPVEC(H(IOFF),HBLK,NI*NJ)
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,'(A,A, 4I3)')
     &   ' Block of single excitations with sym and type ',
     &   ' (ism,itp,jsm,jtp)', ISM,IGAS,JSM,JGAS
         CALL WRTMAT(HBLK,NI,NJ,NI,NJ)
      END IF
*
      RETURN
      END
      SUBROUTINE I_OFF_SX(IOFF,IIGAS,IISM,JJGAS,JJSM,ISX_SYM)
*
* Offset for single excitation block
*
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "multd2h.inc"
*
      IJ_OFF = 1
      DO ISM =1, NSMOB
        JSM = MULTD2H(ISM,ISX_SYM)
        DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
C                 I_SX_CCACT(IGAS,JGAS)
          IJACT = I_SX_CCACT(IGAS,JGAS)
          IF(IJACT.EQ.1) THEN
*
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
*
            IF(ISM.EQ.IISM.AND.IGAS.EQ.IIGAS.AND.
     &         JSM.EQ.JJSM.AND.JGAS.EQ.JJGAS     ) THEN
               IOFF = IJ_OFF
            END IF
            IJ_OFF = IJ_OFF + NI*NJ
          END IF
        END DO
        END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Offset for sx block (gas and sym for i and j )',
     &             IISM,IIGAS,JJSM,JJGAS ,' is ', IOFF
      END IF
*
      RETURN
      END
      SUBROUTINE GET_DX_BLK(IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,
     &                      C,CBLK,IEXP,IXCHNG,IKLJ,IKSM,JLSM,SCR,
     &                      IJ_TRNSP)
*
* Fetch block of double excitation coefficients
*
* IXCNG : Obtain coulom - exchange
* IKLJ  : = 1 => Obtain in dirac form <IK!LJ>
*         = 0 => Obtain in Coulomb form (IJ!KL)
*
* IJ_TRNSP = 1 : Block of interest is C(ji,kl) in form CBLK(IJ,KL)
*
* If the integrals are exported in dirac form
* there is the additional possibility of supplying integrals in
* symmetric forms
* IKSM : use i.ge.k
* IKSM : use j.ge.l
*
* Jeppe Olsen, August 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION C(*)
*.output
      DIMENSION CBLK(*)
*. Scratch
      DIMENSION SCR(*)
*. General information
#include "mxpdim.inc"
#include "orbinp.inc"
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
      NK = NOBPTS(KGAS,KSM)
      NL = NOBPTS(LGAS,LSM)
*. Double excitations assumed symmetrix
      IDX_SYM = 1
*
*. Obtain C_{ijkl} in usual (IJ!KL) form
*
      CALL FETCH_DX_BLK(IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,C,CBLK,
     &                  IDX_SYM,IEXP,IJ_TRNSP)
      IF(IXCHNG.EQ.1) THEN
*. Obtain C_{iljk}
        CALL FETCH_DX_BLK(IGAS,ISM,LGAS,LSM,KGAS,KSM,JGAS,JSM,C,SCR ,
     &                    IDX_SYM,IEXP,IJ_TRNSP)
*. C(I,J,K,L) = C_{ijkl} - C_{ilkj}
        CALL EXCHN_TO_COUL(CBLK,SCR,NI,NJ,NK,NL,1)
C            EXCHN_TO_COUL(C_COUL,C_EXCHN,NI,NJ,NK,NL,ISC)
      END IF
*. Into Dirac form ?
      IF(IKLJ.EQ.1) THEN
        CALL REO_DXBLK_MUL_TO_DIR(NI,NJ,NK,NL,IKSM,JLSM,CBLK,SCR,1)
      END IF
*
      NTEST = 00
      IF(ISM+JSM+KSM+LSM.EQ.16) THEN
        NTEST = 00
      ELSE
        NTEST = 0
      END IF
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*)
        WRITE(6,*) ' ============================='
        WRITE(6,*) ' Output block from GET_DX_BLK '
        WRITE(6,*) ' ============================='
        WRITE(6,*)
*
        IF(IKLJ.EQ.0) THEN
          WRITE(6,*) ' ISM JSM JSM LSM ', ISM,JSM,KSM,LSM
          WRITE(6,*) ' Coulomb form C(IJ,KL) '
          WRITE(6,*)
          CALL WRTMAT(CBLK,NI*NJ,NK*NL,NI*NJ,NK*NL)
        ELSE IF (IKLJ .EQ. 1 ) THEN
          WRITE(6,*) ' Exchange form C(IK,JL) '
          WRITE(6,*)
          IF(IKSM.EQ.1) THEN
            NIK = NI*(NI+1)/2
          ELSE
            NIK = NI*NK
          END IF
          IF(JLSM.EQ.1) THEN
            NJL = NJ*(NJ+1)/2
          ELSE
            NJL = NJ*NL
          END IF
          CALL WRTMAT(CBLK,NIK,NJL,NIK,NKL)
        END IF
*
      END IF
*
      RETURN
      END
      SUBROUTINE REO_DXBLK_MUL_TO_DIR(NI,NJ,NK,NL,IKSM,JLSM,DIN,DOUT,
     &                                ICOPY)
*
* Reorganize block of 4-electron terms from Mullikan to Dirac form
*
* If Icopy .eq. 1, then output block is copied over inputblock
* Jeppe Olsen, August 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. input D(I,J,K,L)
      DIMENSION DIN(NI,NJ,NK,NL)
*. output D(IK,JL)
      DIMENSION DOUT(*)
*
      IF(IKSM.EQ.0) THEN
        NIK = NI*NK
      ELSE
        NIK = NI*(NI+1)/2
      END IF
*
      IF(JLSM.EQ.0) THEN
        NJL = NJ*NL
      ELSE
        NJL = NJ*(NJ+1)/2
      END IF
*
      DO L = 1, NL
       IF(JLSM.EQ.1) THEN
         JMIN = L
       ELSE
         JMIN = 1
       END IF
       DO J = JMIN, NJ
        IF(JLSM.EQ.0) THEN
         JL = (L-1)*NJ + J
        ELSE
C        JL = J*(J-1)/2+L
         JL = (L-1)*NJ + J - L*(L-1)/2
        END IF
        DO K = 1, NK
         IF(IKSM.EQ.1) THEN
          IMIN = K
         ELSE
          IMIN = 1
         END IF
         DO I = IMIN,NI
          IF(IKSM.EQ.0) THEN
           IK = (K-1)*NI+I
          ELSE
C          IK = I*(I-1)/2+K
           IK = (K-1)*NI + I - K*(K-1)/2
          END IF
*
          DOUT((JL-1)*NIK+IK) = DIN(I,J,K,L)
*
         END DO
        END DO
       END DO
      END DO
*
      IF(ICOPY.EQ.1) THEN
       CALL COPVEC(DOUT,DIN,NIK*NJL)
      END IF
*
C     WRITE(6,*) ' NI NJ NK NL NIK NJL', NI,NJ,NK,NL,NIK,NJL
      NTEST = 00
      IF(NTEST.GE.100) THEN
*
        WRITE(6,*) ' Double coefficients in DIRAC format D(IK,JL)'
        WRITE(6,*)
        CALL WRTMAT(DOUT,NIK,NJL,NIK,NJL)
      END IF
*
      RETURN
      END
      SUBROUTINE EXCHN_TO_COUL(C_COUL,C_EXCHN,NI,NJ,NK,NL,ISC)
*
* Form exchange to coulomb format for block of doubles coefs
*
* ISC = 1 C_COUL(I,J,K,L) =  C_COUL(I,J,K,L) - C_EXCHN(I,L,K,J)
* ISC = 2 C_COUL(I,J,K,L) =  C_EXCHN(I,L,K,J)
*
* Jeppe Olsen, August 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION C_COUL(NI,NJ,NK,NL)
*. Output
      DIMENSION C_EXCHN(NI,NL,NK,NJ)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The two input matrices to EXCHN.... '
        CALL WRTMAT(C_COUL ,NI*NJ,NK*NL,NI*NJ,NK*NL)
        CALL WRTMAT(C_EXCHN,NI*NL,NK*NJ,NI*NL,NK*NJ)
      END IF
*
      DO I = 1, NI
       DO J = 1, NJ
        DO K = 1, NK
         DO L = 1, NL
           IF(ISC.EQ.2) THEN
             C_COUL(I,J,K,L) = C_EXCHN(I,L,K,J)
           ELSE
             C_COUL(I,J,K,L) =  C_COUL(I,J,K,L) - C_EXCHN(I,L,K,J)
           END IF
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from COUL_TO_EXCHN '
        WRITE(6,*)
        CALL WRTMAT(C_COUL,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*
      RETURN
      END
      SUBROUTINE FETCH_DX_BLK(IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,
     &                        C,CBLK,IDX_SYM,IEXP,IJ_TRNSP)
*
* Fetch block with given type and sym from list of double excitations
*
* If block is packed and IEXP = 1, then the block is expanded
*
* IF IJ_TRNSP = 1, then we are after the doubles block C(ji,kl) stored
* as Cblk(ij,kl). Introduced to comply with call from RSBB2BN
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. input
      DIMENSION C(*)
*. Output
      DIMENSION CBLK(*)
*. General input
#include "mxpdim.inc"
#include "orbinp.inc"
*. Local scratch (dirty, but simple)
      DIMENSION SCR(200,200)
*
C?    WRITE(6,*) ' First 5 elements of C in FETCH_DX.. '
C?    CALL WRTMAT(C,1,5,1,5)
      IF(IJ_TRNSP.EQ.1) THEN
        L = IGAS
        IGAS = JGAS
        JGAS = L
        L = ISM
        ISM = JSM
        JSM = L
      END IF
*. Is combination of types allowed
      IJKL_ACT = I_DX_ACT(IGAS,KGAS,LGAS,JGAS)
*. Complete or packed block
      IF(IGAS.EQ.KGAS.AND.ISM.EQ.KSM .AND.
     &   JGAS.EQ.LGAS.AND.JSM.EQ.LSM       ) THEN
        ISYM = 1
      ELSE
        ISYM = 0
      END IF
*
      NI = NOBPTS(IGAS,ISM)
      NJ = NOBPTS(JGAS,JSM)
      NK = NOBPTS(KGAS,KSM)
      NL = NOBPTS(LGAS,LSM)
      IF(ISYM.EQ.0) THEN
        LEN = NI*NJ*NK*NL
      ELSE
        LEN = (NI*NJ+1)*NI*NJ/2
      END IF
      IF(IJKL_ACT.EQ.0) THEN
*. Zero block
        ZERO = 0.0D0
        CALL SETVEC(CBLK,ZERO,LEN)
      ELSE
*. Obtain offset for block
        CALL I_OFF_DX(IOFF,ITRNSP,
     &         IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM,IDX_SYM)
        IF(IEXP.EQ.0.OR.ISYM.EQ.0) THEN
          IF(ITRNSP.EQ.0) THEN
            CALL COPVEC(C(IOFF),CBLK(1),LEN)
          ELSE IF(ITRNSP.EQ.1 ) THEN
            CALL TRPMAT(C(IOFF),NK*NL,NI*NJ,CBLK(1) )
          END IF
        ELSE
          SIGN = 1.0D0
          CALL TRIPK2(CBLK(1),C(IOFF),2,NI*NJ,NI*NJ,SIGN)
        END IF
      END IF
*
      IF(IJ_TRNSP.EQ.1) THEN
* C(JI,KL) obtained, transpose to C(IJ,KL)
       DO KL = 1, NK*NL
         KLOFF = 1 + (KL-1)*NI*NJ
         CALL TRPMAT(CBLK(KLOFF),NI,NJ,SCR)
         CALL COPVEC(SCR,CBLK(KLOFF),NI*NJ)
        END DO
*. Clean up
        L = IGAS
        IGAS = JGAS
        JGAS = L
        L = ISM
        ISM = JSM
        JSM = L
      END IF
*
      NTEST = 00
      IF(ISM+JSM+KSM+LSM.EQ.16) THEN
        NTEST = 00
      ELSE
        NTEST = 0
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Block of double-coefficients '
        IF(IJ_TRNSP.EQ.0) THEN
          WRITE(6,*) ' Type and symmetry of indeces (I,J,K,L) ' ,
     &                 IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
        ELSE IF(IJ_TRNSP.EQ.1) THEN
          WRITE(6,*) ' Type and symmetry of indeces (J,I,K,L) ' ,
     &                 JGAS,JSM,IGAS,ISM,KGAS,KSM,LGAS,LSM
         WRITE(6,*) ' (Transposed in indeces I and J ) '
        END IF
        CALL WRTMAT(CBLK,NI*NJ,NK*NL,NI*NJ,NK*NL)
      END IF
*
      RETURN
      END
      SUBROUTINE I_OFF_DX(IOFF,ITRNSP,
     &           IXGAS,IXSM,JXGAS,JXSM,KXGAS,KXSM,LXGAS,LXSM,IDX_SYM)
*
* Obtain offset for  block of double excitations.
*
* Coefficients are stored in blocks C(ai,bj), (ij).ge.(kl)
*
* If ijkl corresponds to (ij).lt.(kl), the ITRNSP flag
* is set, and the offset to the corresponding klij block is returned
* Blocks with (kl).gt. (ij) are referred to the corresponding
* klij block, and the transpose flag is set.
*
*
*. General input
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "lucinp.inc"
*
      IIGAS = IXGAS
      IISM  = IXSM
*
      JJGAS = JXGAS
      JJSM  = JXSM
*
      KKGAS = KXGAS
      KKSM  = KXSM
*
      LLGAS = LXGAS
      LLSM  = LXSM
*. C(ij,kl) or C(kl,ij) ( Transposed)
      ITRNSP = 0
      IJSM_INDEX = (IISM-1)*NSMOB + JJSM
      KLSM_INDEX = (KKSM-1)*NSMOB + LLSM
      IJGAS_INDEX = (IIGAS-1)*NGAS+ JJGAS
      KLGAS_INDEX = (KKGAS-1)*NGAS+ LLGAS
*. Transpose if
      IF(IJSM_INDEX.LT.KLSM_INDEX .OR.
     &   IJSM_INDEX.EQ.KLSM_INDEX.AND.IJGAS_INDEX.LT.KLGAS_INDEX) THEN
        ITRNSP = 1
*
        IIGAS = KXGAS
        IISM  = KXSM
*
        JJGAS = LXGAS
        JJSM  = LXSM
*
        KKGAS = IXGAS
        KKSM  = IXSM
*
        LLGAS = JXGAS
        LLSM  = JXSM
*
      END IF
*
      IJKL_OFF = 1
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
C?     WRITE(6,*) ' ISM JSM KSM LSM', ISM,JSM,KSM,LSM
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,LGAS,JGAS)
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0
            ELSE
              IJBL_GT_KLBL = -1
            END IF
          END IF
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            I_OFF = IOBPTS(IGAS,ISM)
*
            NJ = NOBPTS(JGAS,JSM)
            J_OFF = IOBPTS(JGAS,JSM)
*
            NK = NOBPTS(KGAS,KSM)
            K_OFF = IOBPTS(KGAS,KSM)
*
            NL = NOBPTS(LGAS,LSM)
            L_OFF = IOBPTS(LGAS,LSM)
*
            IF(IIGAS.EQ.IGAS.AND.IISM.EQ.ISM.AND.
     &         JJGAS.EQ.JGAS.AND.JJSM.EQ.JSM.AND.
     &         KKGAS.EQ.KGAS.AND.KKSM.EQ.KSM.AND.
     &         LLGAS.EQ.LGAS.AND.LLSM.EQ.LSM ) THEN
               IOFF = IJKL_OFF
            END IF
*
            IF(IJBL_GT_KLBL.EQ.1) THEN
              IJKL_OFF = IJKL_OFF + NI*NJ*NK*NL
            ELSE IF (IJBL_GT_KLBL .EQ. 0) THEN
              IJKL_OFF = IJKL_OFF + NI*NJ*(NI*NJ+1)/2
            END IF
C?          WRITE(6,*) ' IGAS,JGAS,KGAS,LGAS,IJKL_OFF',
C?   &                   IGAS,JGAS,KGAS,LGAS,IJKL_OFF
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Offset for block of double excitations'
        WRITE(6,*) ' ======================================'
        WRITE(6,*)
        WRITE(6,*) ' Type and Sym for orbitals I,J,K,L ',
     &               IIGAS,IISM,JJGAS,JJSM,KKGAS,KKSM,LLGAS,LLSM
        WRITE(6,*) ' IDX_SYM and IOFF = ', IDX_SYM,IOFF
      END IF
*
      RETURN
      END
      SUBROUTINE RENRM_DX(C,IWAY,IDX_SYM)
*
* Switch normalization of double excitations corresponding
* to switching between restricted and unrestricted summation
*
* Restricted   summation : T2 = sum(ij.ge.kl) CR_{ijkl} e_{ijkl}
* Unrestricted summation : T2 = sum(ij    kl) CU_{ijkl} e_{ijkl}
*
* Relation between restricted and unrestricted summation
*
* CU(ijkl) = (1+delta((ij),kl))/2 CR(ijkl)
*
* IWAY = 1 restricted to unrestricted
* IWAY = 2 unrestricted to restricted

* Jeppe Olsen, Summer of 98
*
*. Input and output
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "lucinp.inc"
*. Specific input/output
      DIMENSION C(*)
*
      IF(IWAY.EQ.1) THEN
        FACTOR = 0.5D0
      ELSE
        FACTOR = 2.0D0
      END IF
      FACTORI = 1.0D0/FACTOR
*
      IJKL_OFF = 1
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0
            ELSE
              IJBL_GT_KLBL = -1
            END IF
          END IF
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
            NK = NOBPTS(KGAS,KSM)
            NL = NOBPTS(LGAS,LSM)
*
            WRITE(6,*) ' REN.. ACTIVE IGAS,JGAS,KGAS,LGAS,IJKL_OFF',
     &      IGAS,JGAS,KGAS,LGAS,IJKL_OFF
            IF(IJBL_GT_KLBL.EQ.1) THEN
              CALL SCALVE(C(IJKL_OFF),FACTOR,NI*NJ*NK*NL)
              IJKL_OFF = IJKL_OFF + NI*NJ*NK*NL
            ELSE IF (IJBL_GT_KLBL .EQ. 0) THEN
              CALL SCALVE(C(IJKL_OFF),FACTOR,(NI*NJ+1)*NI*NJ/2)
              CALL SCLDIA(C(IJKL_OFF),FACTORI,NI*NJ,1)
C                  SCLDIA(A,FACTOR,NDIM,IPACK)
              IJKL_OFF = IJKL_OFF + NI*NJ*(NI*NJ+1)/2
            END IF
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*)
       WRITE(6,*) ' ================================='
       WRITE(6,*) ' Renormalized double coefficients '
       WRITE(6,*) ' ================================='
       WRITE(6,*)
       WRITE(6,*) ' Symmetry = ', IDX_SYM
       CALL WRT_DX1(C,IDX_SYM)
      END IF
*
      RETURN
      END
      FUNCTION I_DX_ACT(IGAS,KGAS,LGAS,JGAS)
*
* Is double excitation a+ igas a+ kgas a lgas a jgas active
*
* Hiding cc restrictions etc
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "cc_exc.inc"
*
      IF(ICC_EXC.EQ.1) THEN
*. Check for coupled cluster restrictions
        I_DX_ACT = I_DX_CCACT(IGAS,KGAS,LGAS,JGAS)
      ELSE
*. Normal CI, extensions for perturbation etc can be inseted here
        I_DX_ACT = 1
      END IF
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
*
      RETURN
      END
      FUNCTION I_SX_ACT(IGAS,JGAS)
*
* Is single excitation a+ igas a jgas active
*
* Hiding cc restrictions etc
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "cc_exc.inc"
*
      IF(ICC_EXC.EQ.1) THEN
*. Check for coupled cluster restrictions
        I_SX_ACT = I_SX_CCACT(IGAS,JGAS)
      ELSE
*. Normal CI, extensions for perturbation etc can be inseted here
        I_SX_ACT = 1
      END IF
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
*
      RETURN
      END
      SUBROUTINE WRT_DX1(C,IDX_SYM)
*
*
* Print list of double excitations in compressed form
* - without singlet-singlet, triplet-triplet separation
* Jeppe Olsen, summer of 98
*

      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "lucinp.inc"
*. Input
      DIMENSION C(*)
*
C     WRITE(6,*)
C     WRITE(6,*) ' ================================================'
C     WRITE(6,*) ' List of double excitations without S/T splitting'
C     WRITE(6,*) ' ================================================'
C     WRITE(6,*)
      WRITE(6,*)
     & '  (Blocks for E(IJ) E(KL) written as matrices C(IJ,KL) )'
      IJKL_OFF = 1
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0
            ELSE
              IJBL_GT_KLBL = -1
            END IF
          END IF
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            I_OFF = IOBPTS(IGAS,ISM)
*
            NJ = NOBPTS(JGAS,JSM)
            J_OFF = IOBPTS(JGAS,JSM)
*
            NK = NOBPTS(KGAS,KSM)
            K_OFF = IOBPTS(KGAS,KSM)
*
            NL = NOBPTS(LGAS,LSM)
            L_OFF = IOBPTS(LGAS,LSM)
*
            IF(NI*NJ*NK*NL.GT.0) THEN
              WRITE(6,'(A,8I3)')
     &        ' Orbital indeces I,J,K,L (type and sym)',
     &          IGAS,ISM,JGAS,JSM,KGAS,KSM,LGAS,LSM
              IF(IJBL_GT_KLBL.EQ.1) THEN
                CALL WRTMAT(C(IJKL_OFF),NI*NJ,NK*NL,NI*NJ,NL*NL)
                IJKL_OFF = IJKL_OFF + NI*NJ*NK*NL
              ELSE IF (IJBL_GT_KLBL .EQ. 0) THEN
                CALL PRSM2(C(IJKL_OFF),NI*NJ)
                IJKL_OFF = IJKL_OFF + NI*NJ*(NI*NJ+1)/2
              END IF
            END IF
*           ^ End if nonvanishing block
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      RETURN
      END
      SUBROUTINE REF_DX_EXP_COMP(CEXP,CCOM,IDX_SYM,IWAY,IBIO)
*
* Reform double excitations between expanded form and
* compressed form - without singlet-triplet separation
*
* IWAY = 1 : Expanded to compressed form
* IWAY = 2 : Compressed to expanded form
*
* Jeppe Olsen, summer of 98
*
*. Modified to column wise packing
*
* Diagonal elements are divided with two.
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "lucinp.inc"
*. Input and output
      DIMENSION CEXP(*),CCOM(*)
*
      FACTOR1 = 1.0D0/3.0D0
      FACTOR2 = 1.0D0/6.0D0
      HALF = 0.5D0
*
C?    WRITE(6,*) ' IBIO = ', IBIO
*
      IJKL_COM = 0
      DO ISM = 1, NSMOB
      DO JSM = 1, NSMOB
      DO KSM = 1, ISM
       IJSM = MULTD2H(ISM,JSM)
       IJKSM = MULTD2H(IJSM,KSM)
       LSM   = MULTD2H(IJKSM,IDX_SYM)
       IF(ISM.GT.KSM.OR.(ISM.EQ.KSM.AND.JSM.GT.LSM)) THEN
         IJSM_GT_KLSM = 1
       ELSE IF( ISM.EQ.KSM.AND.JSM.EQ.LSM) THEN
         IJSM_GT_KLSM = 0
       ELSE
         IJSM_GT_KLSM = -1
       END IF
C?     WRITE(6,*) ' ISM,JSM,KSM,LSM,IJSM_GT_KLSM',
C?   &              ISM,JSM,KSM,LSM,IJSM_GT_KLSM
       IF( IJSM_GT_KLSM.GE.0) THEN
         DO IGAS = 1, NGAS
         DO JGAS = 1, NGAS
         DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
          IJKL_ACT = I_DX_CCACT(IGAS,KGAS,JGAS,LGAS)
C?        WRITE(6,*) ' IGAS,JGAS,KGAS,LGAS,IJKL_ACT',
C?   &                 IGAS,JGAS,KGAS,LGAS,IJKL_ACT
*. Check of block fulfills (IJ.GE.KL)
          IF( IJSM_GT_KLSM .EQ. 1 ) THEN
            IJBL_GT_KLBL = 1
          ELSE IF ( IJSM_GT_KLSM .EQ. 0 ) THEN
            IF(IGAS.GT.KGAS.OR.(IGAS.EQ.KGAS.AND.JGAS.GT.LGAS)) THEN
              IJBL_GT_KLBL = 1
            ELSE IF(IGAS.EQ.KGAS.AND.JGAS.EQ.LGAS) THEN
              IJBL_GT_KLBL = 0
            ELSE
              IJBL_GT_KLBL = -1
            END IF
          END IF
C?        WRITE(6,*) ' IJBL_GT_KLBL' ,  IJBL_GT_KLBL
          IF(IJKL_ACT.EQ.1 .AND. IJBL_GT_KLBL.GE.0 ) THEN
*
            NI = NOBPTS(IGAS,ISM)
            I_OFF = IOBPTS(IGAS,ISM)
*
            NJ = NOBPTS(JGAS,JSM)
            J_OFF = IOBPTS(JGAS,JSM)
*
            NK = NOBPTS(KGAS,KSM)
            K_OFF = IOBPTS(KGAS,KSM)
*
            NL = NOBPTS(LGAS,LSM)
            L_OFF = IOBPTS(LGAS,LSM)
*
            DO L = L_OFF,L_OFF+NL-1
            DO K = K_OFF,K_OFF+NK-1
C
            IF(IJBL_GT_KLBL .EQ. 0 ) THEN
             J_MIN = L
            ELSE
             J_MIN = J_OFF
            END IF
            DO J = J_MIN,J_OFF+NJ-1
*
            IF( IJBL_GT_KLBL .EQ. 1 ) THEN
              I_MIN = I_OFF
            ELSE IF ( IJBL_GT_KLBL .EQ. 0 ) THEN
              IF(J.GT.L) THEN
                I_MIN = I_OFF
              ELSE
                I_MIN = K
              END IF
            END IF
            DO I = I_MIN,I_OFF+NI-1
*
              IJ = (J-1)*NTOOB+I
              KL = (L-1)*NTOOB+K
*
              IF(IBIO.EQ.1) THEN
               IL = (L-1)*NTOOB + I
               KJ = (J-1)*NTOOB + K
               ILKJ_EXP = MAX(IL,KJ)*(MAX(IL,KJ)-1)/2 + MIN(IL,KJ)
              END IF
              IJKL_EXP = IJ*(IJ-1)/2+KL
              IJKL_EXP = MAX(IJ,KL)*(MAX(IJ,KL)-1)/2+MIN(IJ,KL)
              IJKL_COM = IJKL_COM + 1
C?            WRITE(6,*) 'I,J,K,L,IJ,KL,IJKL_EXP,IJKL_COM',
C?   &                    I,J,K,L,IJ,KL,IJKL_EXP,IJKL_COM
              IF(IWAY.EQ.1) THEN
*. expanded to compressed
                IF(IBIO.EQ.0) THEN
                  CCOM(IJKL_COM) = CEXP(IJKL_EXP)
                  IF(IJ.EQ.KL) CCOM(IJKL_COM) = HALF*CCOM(IJKL_COM)
                ELSE
                  CCOM(IJKL_COM) = FACTOR1*CEXP(IJKL_EXP) +
     &                             FACTOR2*CEXP(ILKJ_EXP)
                  IF(IJ.EQ.KL) CCOM(IJKL_COM) = HALF*CCOM(IJKL_COM)
                END IF
              ELSE
                CEXP(IJKL_EXP) = CCOM(IJKL_COM)
              END IF
*
            END DO
            END DO
            END DO
            END DO
*           ^ End of loop over orbitals over given TS
          END IF
*         ^ End if allowed block
         END DO
         END DO
         END DO
         END DO
*        ^ End of loop over gasspaces
       END IF
*      ^ End if IJ_SM .GT. KL_SM
      END DO
      END DO
      END DO
*     ^ End of loop over orbital symmetries
*
      RETURN
      END
      SUBROUTINE WRT_SX(C,ISX_SYM)
*
* Write list of single excitations given in  compact form
*
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "multd2h.inc"
*. Input
      DIMENSION C(*)
*
      WRITE(6,*)
      WRITE(6,*) ' List of single excitations in compressed form '
      WRITE(6,*) ' ============================================= '
      WRITE(6,*)
      IJ_OFF = 1
      DO ISM =1, NSMOB
        JSM = MULTD2H(ISM,ISX_SYM)
        DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
C                 I_SX_CCACT(IGAS,JGAS)
          IJACT = I_SX_CCACT(IGAS,JGAS)
          IF(IJACT.EQ.1) THEN
*
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
*
            IF(NI*NJ.GT.0) THEN
              WRITE(6,'(A,A, 4I3)')
     &        ' Block of single excitations with sym and type ',
     &        ' (ism,itp,jsm,jtp)', ISM,IGAS,JSM,JGAS
              CALL WRTMAT(C(IJ_OFF),NI,NJ,NI,NJ)
              IJ_OFF = IJ_OFF + NI*NJ
            END IF
*           ^ End if nonvanishing block
          END IF
*         ^ End if active block
        END DO
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE REF_SX(CIN,COUT,INFRM,IOUTFRM,ISX_SYM,IBIO)
*
* Reform single excitations
*
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "multd2h.inc"
*. Input and output
      DIMENSION CIN(*),COUT(*)
*
      IF(INFRM.EQ.1.AND.IOUTFRM.EQ.2) THEN
*. restricted to unrestricted
        IWAY = 1
      ELSE
        IWAY = 2
      END IF
*
      IF(IBIO.EQ.0) THEN
        FACTOR = 1.0D0
      ELSE
        FACTOR = 0.5D0
      END IF
*
      IJ_COMP = 0
      DO ISM =1, NSMOB
        JSM = MULTD2H(ISM,ISX_SYM)
        DO IGAS = 1, NGAS
        DO JGAS = 1, NGAS
C                 I_SX_CCACT(IGAS,JGAS)
          IJACT = I_SX_CCACT(IGAS,JGAS)
          IF(IJACT.EQ.1) THEN
*. Offsets
            I_OFF = IOBPTS(IGAS,ISM)
            J_OFF = IOBPTS(JGAS,JSM)
*. Numbers
            NI = NOBPTS(IGAS,ISM)
            NJ = NOBPTS(JGAS,JSM)
*
            DO J = J_OFF,J_OFF+NJ-1
            DO I = I_OFF,I_OFF+NI-1
              IJ_COMP = IJ_COMP + 1
              IF(IWAY.EQ.1) THEN
               COUT((J-1)*NTOOB + I) = CIN(IJ_COMP)*FACTOR
              ELSE
               COUT(IJ_COMP) = CIN((J-1)*NTOOB + I)*FACTOR
              END IF
            END DO
            END DO
          END IF
        END DO
        END DO
      END DO
*
      RETURN
      END
      SUBROUTINE SD_TO_EXC(IREFA,IREFB,IEXCA,IEXCB,NAEL,NBEL,
     &                     NAEXC,IAEXC,NBEXC,IBEXC)
*
* Transfer between slater determinant representation and
* excitation representation
*
* Jeppe Olsen, summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IREFA(NAEL),IEXCA(NAEL)
      INTEGER IREFB(NBEL),IEXCB(NBEL)
*.Output
      INTEGER IAEXC(*),IBEXC(*)
*. For alpha string
*. COding interruptus
      RETURN
      END
      SUBROUTINE ST_TO_EXC(IREF,IEXC,NEL,NEXC,IEXC_OP)
*
* Find excitation that Ex | IREF > = | IEXC >
*
*. Jeppe Olsen, Summer of 98
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      INTEGER IREF(NEL),IEXC(NEL)
*. Output
      INTEGER IEXC_OP(*)
*. Scratch
#include "mxpdim.inc"
      INTEGER IORB(MXPORB),IEL(MXPORB)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reference string '
        CALL IWRTMA(IREF,1,NEL,1,NEL)
        WRITE(6,*) ' Excited string '
        CALL IWRTMA(IEXC,1,NEL,1,NEL)
      END IF
*
*. Strings are assumed in ascending order
*. Find orbital occuring only in one string.
      JEEL = 1
      JDIF = 1
      DO JREL = 1, NEL
       IF(IREF(JREL).GT.IEXC(JEEL)) THEN
*. IEXC(JEEL) is only in IEXC
        IORB(JDIF) = IEXC(JEEL)
        IEL (JDIF) = JEEL
        JEEL = JEEL + 1
        JDIF = JDIF + 1
       ELSE IF(IREF(JREL).LT.IEXC(JEEL) ) THEN
*. IREF(JREL) is only in IREF
        IORB(JDIF) = -IREF(JREL)
        IEL(JDIF) = JREL
        JDIF = JDIF + 1
       END IF
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of differences detected ', NDIF
        WRITE(6,*) ' IORB and IEL arrays '
        CALL IWRTMA(IORB,1,NDIF,1,NDIF)
        CALL IWRTMA(IEL ,1,NDIF,1,NDIF)
      END IF
*. We will write the excitation operator in descending order
* !IEXC> = +/-a+(IEXC(1)) ... a+(IEXC(NDIF))
*             a (IEXC(NDIF+1)) ... a (IEXC(2*NDIF)) ! IREF>
*  with IEXC(I).GT.IEXC(I+1),I=1, NEXC-1,..., NEXC+1, ...., 2*NEXC
*
       ICREA = 0
       IANNI = 0
       DO JEL = 1, 2*NEXC
         IF(IORB(JEL).GT.0) THEN
           ICREA = ICREA + 1
           IEXC_OP(NEXC+1 -ICREA) = IORB(JEL)
         ELSE IF(IORB(JEL).LT.0) THEN
           IANNI = IANNI + 1
           IEXC_OP(2*NEXC+1-IANNI) = ABS(IORB(JEL))
         END IF
       END DO
*. Then, only the sign is missing
*. 1 : sign for moving orbitals to be annihilated outside
*      (and ordering them in ascending order outside )
       ISIGN = 1
       IANNI = 0
       DO JEL = 1, 2*NEXC
         IF(IORB(JEL).LT.0) THEN
           ISIGN = ISIGN*(-1)**(IEL(JEL)-1-IANNI)
           IANNI = IANNI + 1
         END IF
       END DO
*. 2 : Sign for adding electrons to annihilated list,
*      done from list of creation operators in
*.1 : Sign for changing moving all differing electrons in
       NEXC = NDIF/2
*. We have now tabulated all differences
*
* Coding iterruptus
      RETURN
      END
      SUBROUTINE CC_AC_SPACES
*
* Divide orbital spaces ( IHPVGAS ) into
*  Hole spaces     : Only annihilation allowed
*  particle spaces : Only creation allowed
*  valence  spaces : Both annihilation and creation allowed
*
* Division based upon occupation in first CI spaces
* Used for Coupled Cluster Calculations
*
* Jeppe Olsen, Summer of 98 ( not much of an summer !)
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "strinp.inc"
*
      NEL_REF = NELEC(1) + NELEC(2)
C     WRITE(6,*) ' NELEC(1), NELEC(2) ', NELEC(1),NELEC(2)
*
      DO IGAS = 1, NGAS
*.. Occupation in CI space 1 ;
*. hole space : doubly occupied =>2
*. part space : un     occupied
*. val  space : Various occupations
*
*
       IF(IGAS.EQ.1) THEN
         NEL_MAX = 2*NGSOBT(IGAS)
       ELSE
         NEL_MAX = NEL_MAX + 2*NGSOBT(IGAS)
       END IF
*
       IF(IGSOCCX(IGAS,1,1) .EQ. NEL_MAX  .AND.
     &    IGSOCCX(IGAS,2,1) .EQ. NEL_MAX       ) THEN
*. hole space
          IHPVGAS(IGAS) = 1
       ELSE IF(IGAS.GT.1.AND.IGSOCCX(IGAS-1,1,1) .EQ. NEL_REF ) THEN
*. Particle space
          IHPVGAS(IGAS) = 2
       ELSE
*. Valence space
          IHPVGAS(IGAS) = 3
       END IF
*
      END DO
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CC division of orbitals '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' Hole =>1, Part=>2, Val=>3 '
        WRITE(6,*)
        CALL IWRTMA(IHPVGAS,1,NGAS,1,NGAS)
      END IF
*
      RETURN
      END
      FUNCTION I_SX_CCACT(IGAS,JGAS)
*
*  Is excitation a+Igas a Jgas active
*
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CC_EXC/ICC_EXC
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
*
      IF(ICC_EXC .EQ. 0 ) THEN
        IACT = 1
      ELSE IF (ICC_EXC .EQ. 1 ) THEN
        IF(IHPVGAS(IGAS) .GE.2 .AND. IHPVGAS(JGAS) .NE. 2
     &     .AND..NOT.(IHPVGAS(IGAS).EQ.3.AND.IHPVGAS(JGAS).EQ.3) ) THEN
          IF(I_IAD(IGAS).EQ.2.AND.I_IAD(JGAS).EQ.2) THEN
            IACT = 1
          ELSE
            IACT = 0
          END IF
        ELSE
          IACT = 0
        END IF
      END IF
*
      I_SX_CCACT = IACT
*
      RETURN
      END
      FUNCTION I_DX_CCACT(IGAS,KGAS,LGAS,JGAS)
*
*  Is excitation a+Igas a+Kgas a Lgas a Jgas active
*
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/CC_EXC/ICC_EXC
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
*
      IF(ICC_EXC .EQ. 0 ) THEN
        IACT = 1
      ELSE IF (ICC_EXC .EQ. 1 ) THEN
*. Not allowed to excite into inactive, annihilate from
*  secondary, furthermore not all four indeces can be active
        IACT = 0
        IF(I_IAD(IGAS).EQ.2.AND.I_IAD(JGAS).EQ.2.AND.
     &     I_IAD(KGAS).EQ.2.AND.I_IAD(LGAS).EQ.2     ) THEN
          IACT = 1
          IF(IHPVGAS(IGAS).EQ.1.OR.IHPVGAS(KGAS).EQ.1.OR.
     &       IHPVGAS(JGAS).EQ.2.OR.IHPVGAS(LGAS).EQ.2) THEN
             IACT = 0
          END IF
*. 4 indeces in valence space not allowed
          IF(IHPVGAS(IGAS).EQ.3.AND.IHPVGAS(JGAS).EQ.3.AND.
     &       IHPVGAS(KGAS).EQ.3.AND.IHPVGAS(LGAS).EQ.3) THEN
            IACT = 0
          END IF
        END IF
      END IF
*
      I_DX_CCACT = IACT
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*) ' IGAS(c) JGAS(a) KGAS(c) LGAS(a) and IACT',
     &   IGAS,JGAS,KGAS,LGAS,IACT
      END IF
*
      RETURN
      END
      SUBROUTINE FIND_N_CC_AMP(IEXSYM,NSXA,NSXB,NDXAA,NDXBB,NDXAB,
     &                    NSXE,NDXEE)
*
* Number of coupled cluster amplitudes
*
* Jeppe Olsen, Summer of 98
*
* Two sets of amplitudes : (Spin-) Restricted and unrestricted.
*. Currently only unrestricted set is implemented
*
* Restricted
* T = sum(ai) C(ai)E(ai) + 1/2 sum(ai>=bj) C(aibj) E(ai) E(bj)

* Unrestricted
* T = sum(ai) Ca(ai) Ea(ai) + sum(ai) Cb(ai)Eb(ai)
*   + sum(a>b,i>j) Caa(abij) Ea(ai)Ea(bj)
*   + sum(a>b,i>j) Cbb(abij) Eb(ai)Eb(bj)
*   + sum(ab,ij )  Cab(abij) Ea(ai)Eb(bj)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
      COMMON/CC_EXC/ICC_EXC
*. ICC_EXC should be set to one before call
*
*. Single excitations
*
      NSX = 0
      DO IGAS = 1, NGAS
       DO JGAS = 1, NGAS
*. Is this excitation allowed
        IACT = I_SX_ACT(IGAS,JGAS)
        IF(IACT.EQ.1) THEN
          DO ISYM = 1, NSMST
            JSYM = MULTD2H(ISYM,IEXSYM)
            NSX = NSX + NGSOB(ISYM,IGAS)*NGSOB(JSYM,JGAS)
          END DO
        END IF
       END DO
      END DO
      NSXA = NSX
      NSXB = NSX
      NSXE = NSX
*
* Double excitations
*
      NDXAA = 0
      NDXAB = 0
      NDXEE = 0
      DO IGAS = 1, NGAS
       DO JGAS = 1, NGAS
        DO KGAS = 1, NGAS
         DO LGAS = 1, NGAS
           IACT = I_DX_CCACT(IGAS,KGAS,LGAS,JGAS)
           IF(IACT.EQ.1) THEN
            DO ISYM = 1, NSMST
             DO JSYM = 1, NSMST
              DO KSYM = 1, NSMST
*
               IJSYM = MULTD2H(ISYM,JSYM)
               IJKSYM = MULTD2H(IJSYM,KSYM)
               LSYM = MULTD2H(IJKSYM,IEXSYM)
*
               I_INDEX = (ISYM-1)*NGAS + IGAS
               J_INDEX = (JSYM-1)*NGAS + JGAS
               K_INDEX = (KSYM-1)*NGAS + KGAS
               L_INDEX = (LSYM-1)*NGAS + LGAS
*
               NI = NGSOB(ISYM,IGAS)
               NJ = NGSOB(JSYM,JGAS)
               NK = NGSOB(KSYM,KGAS)
               NL = NGSOB(LSYM,LGAS)
*. Alpha-alpha and beta-beta excitations
               NIK = 0
               NJL = 0
               IF(I_INDEX.GT.K_INDEX) THEN
                 NIK =  NI*NK
               ELSE IF(I_INDEX.EQ.K_INDEX) THEN
                 NIK =  NK*(NK-1)/2
               END IF
               IF(J_INDEX.GT.L_INDEX) THEN
                 NJL = NJ*NL
               ELSE IF(J_INDEX.EQ.L_INDEX) THEN
                 NJL = NJ*(NJ-1)/2
               END IF
               NDXAA = NDXAA + NIK*NJL
*. Alpha-beta excitations
               NDXAB = NDXAB + NI*NJ*NK*NL
*. Restricted EE excitations
               IF(I_INDEX.GT.K_INDEX.OR.
     &            I_INDEX.EQ.K_INDEX.AND.J_INDEX.GT.L_INDEX) THEN
*. Allowed sym/types, no restrictions
                   NDXEE = NDXEE + NI*NJ*NK*NL
               ELSE IF(I_INDEX.EQ.K_INDEX.AND.J_INDEX.EQ.L_INDEX) THEN
*. Allowed sym/types, restriced sum
                   NDXEE = NDXEE + (NI*NJ+1)*NI*NJ/2
               END IF
              END DO
             END DO
            END DO
          END IF
         END DO
        END DO
       END DO
      END DO
      NDXBB = NDXAA
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' NSXA,NSXB = ', NSXA,NSXB
       WRITE(6,*) ' NDXAA,NDXBB,NDXAB = ', NDXAA,NDXBB,NDXAB
       WRITE(6,*) ' NDXEE = ', NDXEE
      END IF
*
      RETURN
      END
*
      SUBROUTINE INI_CC_AMP(CC,IFORM)
*
* Initialize Coupled Cluster amplitudes
*
*. Method of initialization depends upon IFORM
*
* IFORM = 1 => Set to zero
* IFORM = 2 => Read in from LU_CCAMP
*
* Jeppe Olsen, Summer of 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "crun.inc"
*(    ^contains number of single and double excitations and N_CC_AMP)
#include "clunit.inc"
*. Amplitudes to be defined
      DIMENSION CC(*)
*
      IF(IFORM.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(CC,ZERO,N_CC_AMP)
        WRITE(6,*) ' Initial set of amplitudes set to zero '
*
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
C       WRITE(6,*) ' Playing around in INI_CC_AMP'
*
      ELSE IF(IFORM.EQ.2) THEN
        WRITE(6,*) ' Reading in CC amplitudes from ', LU_CCAMP
        CALL REWINO(LU_CCAMP)
        DO I = 1, N_CC_AMP
          READ(LU_CCAMP,*) CC(I)
        END DO
      ELSE
        WRITE(6,*) ' Unknown parameter in INI_CC_AMP ',IFORM
        Call Abend2(       ' Unknown parameter in INI_CC_AMP ' )
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Initial set of amplitudes '
        WRITE(6,*) ' ========================= '
        WRITE(6,*)
C            WRT_CC_VEC(CC,LU)
        CALL WRT_CC_VEC(CC,6)
      END IF
*
      RETURN
      END
      SUBROUTINE LUCIA_CC(ISM,ISPC,IPRNT,ECC,II_RESTRT_CC,I_TRANS_WF)
*
* Coupled Cluster calculations with LUCIA
*
* Jeppe Olsen, March 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPRDD
#include "mxpdim.inc"
#include "crun.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "cands.inc"
#include "clunit.inc"
#include "cecore.inc"
      COMMON/CC_EXC/ICC_EXC
*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'LUC_CC')
      WRITE(6,*)
      WRITE(6,*) '===================='
      WRITE(6,*) 'LUCIA_CC in control '
      WRITE(6,*) '===================='
      WRITE(6,*)
      WRITE(6,*) ' Symmetry and space ', ISM,ISPC
*.Transfer to CANDS
      ICSM = ISM
      ISSM = ISM
      ICSPC = ISPC
      ISSPC = ISPC
*. Complete operator
      I12 = 2

*. Coupled cluster flag
      ICC_EXC = 1

*. Divide orbital space into inactive, active, secondary
      CALL CC_AC_SPACES
*. Number of single excitation and double excition amplitudes
      IEXSYM = 1
      CALL FIND_N_CC_AMP
     &(IEXSYM,NSXA,NSXB,NDXAA,NDXBB,NDXAB,NNSXE,NNDXEE)
*
      NSXE = NNSXE
      NDXEE = NNDXEE
*. (NSXE and NDXEE are stored in /CRUN/
*. Allocate space for 3 CC vectors
* : Three type of operators
*  1 : single excitations
*  2 : Double excitations
*  3 : Determinants orthogonal to reference (I guess i am going multiref cc)
* The amplitudes will be stored as above
      NNDET = 0
      N_CC_AMP = NSXE + NDXEE + NNDET
      WRITE(6,*) ' NSXE NDXEE N_CC_AMP', NSXE,NDXEE,N_CC_AMP
      LEN = N_CC_AMP
      CALL MEMMAN(KCC1,2*LEN,'ADDL  ',2,'CC1  ')
      CALL MEMMAN(KCC2,2*LEN,'ADDL  ',2,'CC2  ')
      CALL MEMMAN(KCC3,2*LEN,'ADDL  ',2,'CC3  ')
*. Scratch space for CI - behind the curtain
       CALL GET_3BLKS(KVEC1,KVEC2,KVEC3)
*. Initialize CC amplitudes in CC1
C          INI_CC_AMP(CC_AMP,IFORM)
      IF(II_RESTRT_CC.EQ.0) THEN
       IFORM = 1
      ELSE
       IFORM = 2
      END IF
      CALL INI_CC_AMP(WORK(KCC1),IFORM)
*. Find CC wave function
      CALL OPTIM_CC(WORK(KCC1),WORK(KVEC1),WORK(KVEC2))
*. Calculate CC vector function for initial set of amplitudes
*. And store in WORK(KCC3)
      IBIO = 1
C     CALL CC_VEC_FNC(WORK(KCC1),WORK(KCC3),ECC,WORK(KVEC1),WORK(KVEC2),
C    &                IBIO)
*
* =======================
* Variational energy etc
* =======================
*
*
*
* Exp(t) !ref> = sum(n=1,large) 1/n! T^n !ref>
*
      MX_TERM = 100
      ICC_EXC = 1
      XCONV = 1.0D-14
      CALL EXPT_REF(LUC,LUSC1,LUHC,LUSC2,LUSC3,XCONV,MX_TERM,
     &             WORK(KVEC1),WORK(KVEC2))
*
*
*. H Exp T !ref>
*
      ICC_EXC = 0
      CALL MV7(WORK(KVEC1),WORK(KVEC2),LUSC1,LUHC)
*. E = <ref! exp(T)+ H exp(T)!ref>/<ref! exp(T)+exp(T)!ref>
      LBLK = -1
      CHC = INPRDD(WORK(KVEC1),WORK(KVEC2),LUSC1,LUHC  ,1,LBLK)
      CC  = INPRDD(WORK(KVEC1),WORK(KVEC2),LUSC1,LUSC1 ,1,LBLK)
      WRITE(6,*)
      WRITE(6,*) ' Energy as coupled cluster expectation value : '
      WRITE(6,*) ' ============================================='
      WRITE(6,*)
      WRITE(6,'(5X,A,E25.12)')
     &'  <ref! exp(T)+ H exp(T)!ref> (- Ecore)= ', CHC
      WRITE(6,'(5X,A,E25.12)')
     &'  <ref! exp(T)+   exp(T)!ref> = ', CC
      WRITE(6,'(5X,A,E25.12)')
     &' Expectation value coupled cluster energy   = ', CHC/CC + ECORE
      WRITE(6,*)
*
      IF(I_TRANS_WF.EQ.1) THEN
        WRITE(6,*) ' cc wf is normalized and transferred to LUC '
        XNORM = SQRT(CC)
        FACTOR = 1.0D0/XNORM
C            SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
        CALL SCLVCD(LUSC1,LUC,FACTOR,WORK(KVEC1),1,-1)
      END IF

      WRITE(6,*) ' Returning from LUCIA_CC '

      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'LUC_CC')
*
      RETURN
      END
C       EXPT_REF(LUC,LUHC,LUSC1,LUSC2,THRES_E,MX_TERM,
C    &             WORK(KVEC1),WORK(KVEC2))
      SUBROUTINE EXPT_REF(LUC,LUHC,LUSC1,LUSC2,LUSC3,THRES_C,MX_TERM,
     &                    VEC1,VEC2)
*
* Obtain Exp (T) !ref> by Taylor expansion of exponential
*
* Jeppe Olsen, March 1998
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPRDD
*
      DIMENSION VEC1(*),VEC2(*)
      COMMON/CINT_CC/INT_CC
*
      LBLK = -1
*
C?    MX_TERM = 10
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF'
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF'
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF'
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF'
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF'
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF'
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF'
C?    WRITE(6,*) ' Note MX_TERM reduced to 10  in EXPT_REF'
*
      NTEST = 001
      IF(NTEST.GE.1) THEN
       WRITE(6,*)
       WRITE(6,*) '==================='
       WRITE(6,*) 'EXPT_REF in action '
       WRITE(6,*) '==================='
       WRITE(6,*)
      END IF
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' LUC,LUHC,LUSC1,LUSC2',LUC,LUHC,LUSC1,LUSC2
       WRITE(6,*) ' Initial vector on LUC '
       CALL WRTVCD(VEC1,LUC,1,LBLK)
      END IF
* Tell integral fetcher to fetch cc amplitudes, not integrals
      INT_CC = 1
*. Loop over orders of expansion
      N = 0
      XFACN = 1.0D0
*
      CALL COPVCD(LUC,LUSC1,VEC1,1,LBLK)
      CALL COPVCD(LUC,LUSC3,VEC1,1,LBLK)
*
 1000 CONTINUE
       N = N+1
       IF(NTEST.GE.5) THEN
         WRITE(6,*) ' Info for N = ', N
       END IF
*. T^N  times vector on LUSC1
C?     WRITE(6,*) ' Input vector to MV7 '
C?     CALL WRTVCD(VEC1,LUSC1,1,LBLK)
       CALL MV7(VEC1,VEC2,LUSC1,LUHC)
       CALL COPVCD(LUHC,LUSC1,VEC1,1,LBLK)
       IF(NTEST.GE.500) THEN
         WRITE(6,*) ' T**(N) |0> '
         WRITE(6,*) ' ==========='
         CALL WRTVCD(VEC1,LUSC1,1,LBLK)
       END IF
*. Norm of this correction term
       LBLK = -1
       XNORM2 = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
       XFACN = XFACN/FLOAT(N)
       XNORM = SQRT(XNORM2)/XFACN
       IF(NTEST.GE.5) THEN
         WRITE(6,*) ' Norm of correction ', XNORM
       END IF
*. Update output file with 1/N! T^N !ref>
       ONE = 1.0D0
       CALL VECSMD(VEC1,VEC2,ONE,XFACN,LUSC3,LUSC1,LUSC2,1,LBLK)
       CALL COPVCD(LUSC2,LUSC3,VEC1,1,LBLK)
*. Take another turn ?
      IF(XNORM.GT. THRES_C.AND. N .LT. MX_TERM) GOTO 1000
*
*. Result on LUHC
      CALL COPVCD(LUSC3,LUHC,VEC1,1,LBLK)
      IF(NTEST.GT.0) THEN
        WRITE(6,*) ' Convergence obtained in ', N, ' iterations'
        WRITE(6,*) ' Norm of last correction ', XNORM
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============'
        WRITE(6,*) ' Exp T |ref> '
        WRITE(6,*) ' ============'
        WRITE(6,*)
        CALL WRTVCD(VEC1,LUHC,1,LBLK)
      END IF
*
      RETURN
      END








