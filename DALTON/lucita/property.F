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
      SUBROUTINE GENH1(XLAMBDA_X)
*
* Construct the general one-electron operator
*
* H = XLAMBDA*H(NORMAL) + (1-XLAMBDA)FIFA
*
*
* Where H(Normal) is the normal one-electron operator
* and FIFA is the sum of the inactive and active Fock matrices 
* used in CASPTN theory
*
* The correct one-electron density is assumed in place
*
* Jeppe Olsen, March 1996
*
*. Note : Correct Lambda is transferred through CRUN as of Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "oper.inc"
#include "crun.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
*. Construct FIFA in WORK(KFI)
*
      IF(IPART.NE.3) THEN
*. Normal M-P Partitioning
        WRITE(6,*) 'Normal MP partitioning'
        WRITE(6,*) 'Normal MP partitioning'
        WRITE(6,*) 'Normal MP partitioning'
        CALL COPVEC(WORK(KINT1O),WORK(KFI),NINT1)
        CALL FIFAM(WORK(KFI))
        CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
        ECORE_H = 0.0D0
        IF(IUSE_PH.EQ.1) THEN
          CALL FI(WORK(KFI),ECORE_H,0)
          ECORE = ECORE_ORIG + ECORE_H 
        END IF
*. And write to disc
        LU18 = 18
        REWIND (LU18)
        WRITE(6,*) ' H0 written to disc '
        CALL TODSC_LUCI(WORK(KFI),NINT1,-1,LU18)
      ELSE IF(IPART.EQ.3) THEN
        WRITE(6,*) ' Zero-order Hamiltonian readin '
        WRITE(6,*) ' Zero-order Hamiltonian readin '
        WRITE(6,*) ' Zero-order Hamiltonian readin '
        WRITE(6,*) ' Zero-order Hamiltonian readin '
        LU18 = 18
        REWIND (LU18)
        CALL FRMDSC_LUCI(WORK(KFI),NINT1,-1,LU18,IMZERO,IAMPACK)
      END IF
*
*. And obtain modified operator in KINT1 ( No return !! )
*
      FAC2 = 1.0D0 - XLAMBDA
      CALL VECSUM(WORK(KINT1),WORK(KINT1),WORK(KFI),XLAMBDA,FAC2,
     &            NINT1)
*
      CALL COPVEC(WORK(KINT1),WORK(KINT1O),NINT1)
      NTEST = 100
      IF(NTEST.GE.100) THEN
       WRITE(6,*) 'Modified matrix as delivered by GENH1 '
       WRITE(6,*) '======================================'
       CALL APRBLM2(WORK(KINT1),NTOOBS,NTOOBS,NSMOB,1)
      END IF
*
      RETURN
      END 
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE CI_RESPONS(LU1,LU2,LU3,LU4,LU5,LU6,LU7,LUC,LUDIA,
     &                      VEC1,VEC2,ENOT,CCALC)
*
* Master routine for LUCIA calculation of properties
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KINTOP1, KINTOP2, KLSCR, KINTAV, KLENM, KINTAV2
!               for addressing of WORK
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "lucinp.inc"
#include "glbbas.inc"
#include "orbinp.inc"
#include "multd2h.inc"
#include "crun.inc"
#include "cprnt.inc"
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*)
*. Local scratch
       INTEGER IAVE_SYM(20)
*
      NTEST = 1000
*
      IDUM = -1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'RESPON')
      IRFSM = ICSM
*
C?    write(6,*) ' ENOT IN CIRESP ', ENOT
*
      WRITE(6,*) ' Number of CI response calculations ', NRESP
*
      DO ICALC= 1, NRESP
       WRITE(6,*)
       WRITE(6,'(A)')
     & '  ********************************************** '
       WRITE(6,'(A,I3)')
     & '  Information about response calculation ..',ICALC
       WRITE(6,'(A)')
     & '  ********************************************** '
       WRITE(6,*)
*
*
**  General procedure for for double perturbation theory
*                 through arbitrary order
*

          WRITE(6,*)
          WRITE(6,*) ' ================================================'
          WRITE(6,*) ' General expansion in two external perturbations'
          WRITE(6,*) ' ================================================'
          WRITE(6,*)
          WRITE(6,*)
          WRITE(6,*)
     &    ' Operators used for expectation values ( A-operators) '
          WRITE(6,*)
     &    ' ====================================================='
          WRITE(6,*)
          DO IAVE = 1, N_AVE_OP
            WRITE(6,'(1H ,5X,A)') AVE_OP(IAVE)
          END DO
          WRITE(6,*)
          WRITE(6,*) ' Perturbations : '
          WRITE(6,*) ' =============== '
          WRITE(6,*)
          WRITE(6,*) '     Operator 1 : ', RESP_OP(1,ICALC)
          WRITE(6,*) '     Operator 2 : ', RESP_OP(2,ICALC)
          WRITE(6,*)
     &    '     Max order in operator 1 : ', MAXORD_OP(1,ICALC)
          WRITE(6,*)
     &    '     Max order in operator 2 : ', MAXORD_OP(2,ICALC)
          IF(RESP_W(ICALC).EQ.0.0D0) THEN
            WRITE(6,*) '     Static perturbations'
          ELSE
            WRITE(6,*)
     &      '    frequency of second operator  ',RESP_W(ICALC)
          END IF
          WRITE(6,*)
*. Space for one_electron integrals
          CALL MEMMAN(KINTOP1,  NTOOB**2,'ADDL  ',2,'INTOP1')
          CALL MEMMAN(KINTOP2,  NTOOB**2,'ADDL  ',2,'INTOP2')
          CALL MEMMAN(KLSCR  ,5*NTOOB**2,'ADDL  ',2,'RSPSCR')
          CALL MEMMAN(KINTAV,N_AVE_OP*NTOOB **2,'ADDL  ',2,'INTAVE')
*. Energy corrections
          NCORR = (MAXORD_OP(1,ICALC)+1)*(MAXORD_OP(2,ICALC)+1)
          CALL MEMMAN(KLENM,NCORR,'ADDL  ',2,'ENM   ')
          
*. Symmetry of perturbations
C     SYM_FOR_OP(OPER,IXYX,ISYM)
*. W : Op1, V : Op 2
          CALL SYM_FOR_OP(RESP_OP(1,ICALC),IXYZSYM,IWSM)
          CALL SYM_FOR_OP(RESP_OP(2,ICALC),IXYZSYM,IVSM)
*. Symmetry of average operators
          DO I_AVE_OP = 1, N_AVE_OP
            CALL SYM_FOR_OP(AVE_OP(I_AVE_OP),IXYZSYM,
     &                      IAVE_SYM(I_AVE_OP))
          END DO
*. Perturbations times reference
* W !0>
          IWRFSM = MULTD2H(IWSM,IRFSM)
* V !0>
          IVRFSM = MULTD2H(IVSM,IRFSM)
* VW !0>
          IWVRFSM = MULTD2H(IWSM,IVRFSM)
*. Obtain property integrals
          CALL GET_PROPINT(WORK(KINTOP1),IWSM,RESP_OP(1,ICALC),
     &                     WORK(KLSCR),NTOOBS,NTOOBS,NSMOB,1 )
          CALL GET_PROPINT(WORK(KINTOP2),IVSM,RESP_OP(2,ICALC),
     &                     WORK(KLSCR),NTOOBS,NTOOBS,NSMOB,1 )
          DO I_AVE_OP = 1, N_AVE_OP
            KINTAV2 = KINTAV + (I_AVE_OP-1)*NTOOB**2
            CALL GET_PROPINT(WORK(KINTAV2),IAVE_SYM(I_AVE_OP),
     &           AVE_OP(I_AVE_OP),WORK(KLSCR),NTOOBS,NTOOBS,NSMOB,1 )
          END DO
*
          MAXW = MAXORD_OP(1,ICALC)
          MAXV = MAXORD_OP(2,ICALC)
C?        WRITE(6,*) ' MAXW, MAXV ', MAXW,MAXV
*
          FREQ = RESP_W(ICALC)
          NTOOB2 = NTOOB ** 2
          IF(FREQ.EQ.0.0D0) THEN
*. Static perturbation
            ZERO = 0.0D0
            CALL GNDBPTFREQ(WORK(KINTOP2),WORK(KINTOP1),IRFSM,
     &           ENOT,LUC,
     &           IVRFSM,MAXV,IWRFSM,MAXW,IWVRFSM,IVSM,IWSM,
     &           LUN,VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &           LUDIA,ZERO,WORK(KLENM),
     &           N_AVE_OP,AVE_OP,WORK(KINTAV),NTOOB2,IAVE_SYM,
     &           IPRRSP)
          ELSE
*. Frequency dependent perturbation
            CALL GNDBPTFREQ(WORK(KINTOP2),WORK(KINTOP1),IRFSM,
     &           ENOT,LUC,
     &           IVRFSM,MAXV,IWRFSM,MAXW,IWVRFSM,IVSM,IWSM,
     &           LUN,VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &           LUDIA,FREQ,WORK(KLENM),
     &           N_AVE_OP,AVE_OP,WORK(KINTAV),NTOOB2,IAVE_SYM,
     &           IPRRSP)
        END IF
*       ^ End of switch between static/dynamic perturbation
      END DO
*.    ^ End of loop over calculations
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'RESPON')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_PROP_RANK(PROPER,IRANK)
*
* Obtain rank for property PROPERTY (gives as CHAR*6 as usual)
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*6 PROPER
*
      IF(PROPER.EQ.'DIPOLE') THEN
        IRANK = 1
      ELSE IF(PROPER.EQ.'THETA ' .OR.
     &        PROPER.EQ.'QUADRU' .OR.
     &        PROPER.EQ.'SECMOM' .OR.
     &        PROPER(1:3).EQ.'EFG' ) THEN
        IRANK  = 2
      ELSE
        WRITE(6,'(A,A)') ' Unknown operator ',PROPER
        IRANK = -1
      END IF
*
      NTEST = 0
      IF(NTEST.GE.5) THEN
        WRITE(6,'(A,A,I3)') ' Property and rank : ', PROPER,IRANK
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GNDBPTFREQ(V,W,IRFSM,ERF,LURF,
     &           IVRFSM,MAXV,IWRFSM,MAXW,IDRFSM,IVSM,IWSM,
     &           LUN,VEC1,VEC2,LU1,LU2,LU3,LU4,LU5,LU6,LU7,
     &           LUDIA,FREQ,ENM,
     &           N_AVE_OP,AVE_OP,AVINT,LAVINT,IAVE_SYM,IPRNT)
*
* Solve the double perturbation expansion arising from a
* static perturbation(V) and a frequency dependent perturbation(W).
*
* Susceptibilities corresponding to expectation values of
* the N_AVE_OP operators defined by AVINT are also constructed
*
* The perturbations or the squared perturbations are assumed
* to be symmetric
*
* Jeppe Olsen, adapted from LUCAS codes
* The wave function corrections reads
*
* d(n,m,lw) = -(e[2] - lw I)-1 *
*      ( Vw[2] d(n-1,m,(l-1)w)
*       +Vo[2] d(n,m-1,kw)
*       +V-w[2]d(n-1,m,(l+1)w)
*       -sum(ld=1,n;ls=1,m-1;ll)Vo[1]t  d(ld,ls,lw)d(n-ld,m-ls-1,(l-ll)w)
*       -sum(ld=1,n-1;ls=1,m;ll)Vw[1]t  d(ld,ls,lw)d(n-ld-1,m-ls,(l-ll-1)w)
*       -sum(ld=1,n.1;ls=1,m;ll)V-w[1]t d(ld,ls,lw)d(n-ld-1,m-ls,(l-ll+1)w)
*
* Where
* n : is the order in the frequency dependent expansion
* m : is the order in the static expansion
* w : is the frequency (FREQ)
*
*
*
* Input
* =====
* V : one-electron integrals defining  V
* W : one-electron integrals defining  W
* IRFSM : symmetry of reference vector
* LURF : file containing reference vector
* IVRFSM : symmetry of V times reference vector
* IWRFSM : symmetry of W times reference vector
* IDRFSM : symmetry of VW times reference vector
* LUN : file number for file to contain perturbation vectors
* ERF   : Reference energy
* MAXV : Order in V through which the equations should be solved
* MAXW : Order in W through which the equations should be solved
* VEC1,VEC2 : Scratch vectors , should be able to hold a complete vector
* LU1, LU2, LU3, LU4, LU5, LU6, LU7  : scratch files
*
* N_AVE_OP : Number of operators for which expectation values are expanded
*            (A- operators)
* AVE_OP   : Labels of above operators
* AVINT    : Integrals of A operators AVINT(*,IAV) contains integrals
*            for A operator number IAV
* LAVINT   : Row length of AVINT
* IAVE_SYM : Symmetry of A operators
*
* FREQ
*
* Output
* ======
* LUN : contains the MAXORD perturbation vectors
*       Correction nm is stored as record (n-1)*MAXV + m
* ENM : Contains the energy corrections
*
* Internal links
* ===============
* The actual solutions of linear eqs and the multiplicartion
* of V times a vector is realized through calls to BVEC and HMWITV
* /CANDS/ and /OPER( is used to communicate to these routines
*
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL * 8 INPRDD
      DIMENSION V(*),VEC1(*),VEC2(*)
      DIMENSION ENM(MAXV+1,MAXW+1)
*
      DIMENSION ISMVW(2,2)
      PARAMETER(MAXORD = 10)
      DIMENSION V01D(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
      DIMENSION VW1D(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
*. Inner products <d!d>
      DIMENSION DD(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
*. Inner products <d!W!d>
      DIMENSION DWD(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
*. And susceptibilities
      DIMENSION CHI(0:MAXORD,0:MAXORD,-MAXORD:MAXORD)
*. A-operators
      CHARACTER*8 AVE_OP(N_AVE_OP)
      DIMENSION  AVINT(LAVINT,N_AVE_OP),IAVE_SYM(N_AVE_OP)
C    &           N_AVE_OP,AVE_OP,AVINT,LAVINT,IAVE_SYM)
*  V ** n W ** m !ref> : Symmetry is  ISMVW(MOD(n,2)+1,MOD(m,2)+1)
      COMMON/CANDS/ICSM,ISSM,ICPSC,ISSPC
*. Local scratch
      DIMENSION SCR(10*MAXORD)
#include "oper.inc"
#include "multd2h.inc"
*. Transfer of shift
      COMMON/CSHIFT/SHIFT,IPROJ
*. Transfer of zero order energy
      COMMON/CENOT/E0
*. Still testing after all these years
      NTESTL =   1
      WRITE(6,*) ' Input Print flag ', IPRNT
      NTEST = MAX(NTESTL,IPRNT)
*. redefine LUN
      LUN = 41
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' LU1, LU2, LU3, LU4, LU5, LU6, LU7 '
        WRITE(6,*)   LU1, LU2, LU3, LU4, LU5, LU6, LU7
        WRITE(6,*) ' LUN LURF ', LUN,LURF
      END IF
C     WRITE(6,*) ' FREQ = ', FREQ
*
      E0 = 0.0D0
      LBLK = -1
*
      NNNVEC = 0
      ISMVW(1,1) = IRFSM
      ISMVW(2,1) = IVRFSM
      ISMVW(1,2) = IWRFSM
      ISMVW(2,2) = IDRFSM
*
      IF(NTEST.GE.5) THEN
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) '*************************************************'
      WRITE(6,*) '*                                               *'
      WRITE(6,*) '* Welcome to GNDBPTFREQ                         *'
      WRITE(6,*) '* General double perturbation calculations      *'
      WRITE(6,*) '* Involving a static and a dynamic perturbation *'
      WRITE(6,*) '*                                               *'
      WRITE(6,*) '*                      Best wishes from         *'
      WRITE(6,*) '*                      Antonio and Jeppe        *'
      WRITE(6,*) '*                                               *'
      WRITE(6,*) '*************************************************'
      WRITE(6,*)
      WRITE(6,*)
      END IF
C     WRITE(6,*)  ' LUN and LURF ',LUN,LURF

*.  Real perturbation
      IRC = 1
*. For a imaginary pert set IRC to 2
*
*
*. The vectors P V0!0> and P Vw!0> ( P = 1 - !0><0!)
*  will be used a number of times in the
*. following. Construct these guys and store on LU5
*
*
      CALL REWINE(LU5,-1)
*
* V !0> as first vector on LU5
*
      ICSM = IRFSM
      ISSM = IVRFSM
      IPERTOP = 0
      IAPR = 0
C          BVEC(B,IBSM,LUC,LUB,VEC1,VEC2)
      CALL BVEC(V,IVSM,LURF,LU5,VEC1,VEC2)
      IF(IVSM.NE.1) THEN
       V00 = 0.0D0
      ELSE
       V00 = INPRDD(VEC1,VEC2,LU5,LURF,1,LBLK)
       V1NORM = INPRDD(VEC1,VEC2,LU5,LU5,1,LBLK)
       IF(NTEST.GE.10) THEN
         WRITE(6,*) ' V1NORM = ', V1NORM
         WRITE(6,*) ' V00 = ', V00
       END IF
       ONE = 1.0D0
       FACTOR = -V00
       CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU5,LURF,LU4,1,LBLK)
       CALL COPVCDP(LU4,LU5,VEC1,1,LBLK)
      END IF
*
* W !0> as second vector on LU5
*
C?    WRITE(6,*) ' IWSM = ', IWSM
      ISSM  = IWRFSM
      CALL BVEC(W,IWSM,LURF,LU4,VEC1,VEC2)
      IF(IWSM.NE.1) THEN
       W00 = 0.0D0
       CALL REWINE(LU4,-1)
       CALL COPVCDP(LU4,LU5,VEC1,0,LBLK)
      ELSE
       W00 = INPRDD(VEC1,VEC2,LU4,LURF,1,LBLK)
C?     WRITE(6,*) ' W00 = ', W00
       ONE = 1.0D0
       FACTOR = -W00
       CALL REWINE(LU4,-1)
       CALL REWINE(LURF,-1)
       CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU4,LURF,LU5,0,LBLK)
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The two vectors on LU5 '
        CALL REWINE(LU5,-1)
        DO IVEC = 1, 2
           WRITE(6,*) ' Vector ', IVEC
           CALL WRTVCD(VEC1,LU5,0,LBLK)
        END DO
      END IF

*
*. A note on the organization of the correction vectors on disc :
*  The vectors are stored as Loop over freq dep op, Loop over static op,
*  Loop over freq. components

      DO N = 0, MAXV
*
       IF(N.EQ.0) THEN
        MMIN = 1
       ELSE
        MMIN = 0
       END IF
*
       DO M = MMIN, MAXW
*. Symmetry of this correction  vector
         MNSM  = ISMVW(MOD(N,2)+1,MOD(M,2)+1)
*. Frequency components
         DO LW = -N,N,2
*
          WRITE(6,*)
          WRITE(6,*) ' Correction vector will be solved for '
          WRITE(6,*)
          WRITE(6,*) '    Order in dynamic perturbation ', N
          WRITE(6,*) '    Order in static  perturbation ', M
          WRITE(6,*) '    Frequency component           ', LW
          WRITE(6,*) '    Symmetry of correction vector ', MNSM
          WRITE(6,*)
*
          ISSM = ISMVW(MOD(N,2)+1,MOD(M,2)+1)

* =========================================================================
* 1 : construct (Vw-Vw00)!N-1,M,LW-1>,(Vw-Vw00)!N-1,M,LW+1>,(W-W00)!N,M-1,LW> ( and save on LU4)
* =========================================================================
          CALL REWINE(LU4,-1)
*. Components !N-1,M,LW-1>, !N-1,M,LW+1> corresponds to allowed freq.
          IF(N.GE.1.AND.(ABS(LW-1).LE.N-1)) THEN
            IMIN = 1
          ELSE
            IMIN = 0
          END IF
*
          IF(N.GE.1.AND.(ABS(LW+1).LE.N-1)) THEN
            IPLUS= 1
          ELSE
            IPLUS= 0
          END IF
C?        WRITE(6,*) ' N LW IPLUS',N,LW,IPLUS
          ICSM = ISMVW(MOD(N-1,2)+1,MOD(M,2)+1)
*
* (Vw - <0!Vw!0>) !N-1,M,LW-1>( on LU4)
*
          IF(IMIN.EQ.1) THEN
            IF(NTEST.GE.100) WRITE(6,*) ' IMIN active '
* : Place !0(N-1,M,LW-1)> on LU6
            CALL GET_NMLW(VEC1,N-1,M,LW-1,MAXV,MAXW,LUN,LURF,1,
     &                    LU6,1,LBLK,LUOUTEFF)
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Correction vector !0(N-1,M,LW-1)> read '
            END IF
            IF(NTEST.GE.100) THEN
              CALL WRTVCD(VEC1,LU6,1,LBLK)
            END IF
            I12 = 1
            CALL BVEC(V,IVSM,LU6,LU7,VEC1,VEC2)
*. - <0!Vw!0> !N-1,M,LW-1>
C           IF(V00.NE.0.0D0 .AND. N+M.NE.1) THEN
            IF(V00.NE.0.0D0 ) THEN
              ONE = 1.0D0
              FACTOR = -V00
              CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU7,LU6,LU4,1,LBLK)
            ELSE
              CALL COPVCDP(LU7,LU4,VEC1,1,LBLK)
            END IF
*
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' V times Correction vector on LU4 '
              CALL WRTVCD(VEC1,LU4,1,LBLK)
            END IF
          END IF
*
* Vw!N-1,M,LW+1> on LU4
*
          IF(IPLUS.EQ.1) THEN
* : Place !0(N-1,M,LW+1)> in VEC1
            IF(NTEST.GE.100) WRITE(6,*) ' IPLUS active '
            CALL GET_NMLW(VEC1,N-1,M,LW+1,MAXV,MAXW,LUN,LURF,1,
     &                    LU6,1,LBLK,LUOUTEFF)
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' Correction vector !0(N-1,M,LW+1)> read '
            END IF
            IF(NTEST.GE.100) THEN
              CALL WRTVCD(VEC1,LU6,1,LBLK)
            END IF
            I12 = 1
            CALL BVEC(V,IVSM,LU6,LU7,VEC1,VEC2)
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' V times Correction vector'
              CALL WRTVCD(VEC2,LU7,1,LBLK)
            END IF
*. - <0!Vw!0> !N-1,M,LW+1>
C           IF(V00.NE.0.0D0 .AND. N+M .NE. 1) THEN
            IF(V00.NE.0.0D0 ) THEN
              ONE = 1.0D0
              FACTOR = -V00
              CALL REWINE(LU7,-1)
              CALL REWINE(LU6,-1)
              CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU7,LU6,LU4,0,LBLK)
            ELSE
              CALL REWINE(LU7,-1)
              CALL COPVCDP(LU7,LU4,VEC1,0,LBLK)
            END IF
          END IF
*
*  W!N,M-1,LW> on LU4
*
          IF(M.GT.0.AND.ABS(LW).LE.N) THEN
            INOT = 1
          ELSE
            INOT = 0
          END IF
          IF(INOT.EQ.1) THEN
            IF(NTEST.GE.100) WRITE(6,*) ' INOT active '
            ICSM = ISMVW(MOD(N,2)+1,MOD(M-1,2)+1)
*. !N,M-1,LW> on LU6
            CALL GET_NMLW(VEC1,N,M-1,LW,MAXV,MAXW,LUN,LURF,1,
     &                    LU6,1,LBLK,LUOUTEFF)
            IF(NTEST.GE.100) THEN
             WRITE(6,*) ' !n,m-1> read in'
            END IF
            IF(NTEST.GE.100) THEN
              CALL WRTVCD(VEC1,LU6,1,LBLK)
            END IF
C?          RHSNORM = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
C?          write(6,*) ' norm of input vector to BVEC ',RHSNORM
            I12 = 1
            CALL REWINE(LU7,-1)
            CALL BVEC(W,IWSM,LU6,LU7,VEC1,VEC2)
C?          RHSNORM = INPRDD(VEC1,VEC2,LU7,LU7,1,LBLK)
C?          WRITE(6,*) ' Norm of RHS after BVEC ', RHSNORM

            IF(NTEST.GE.10000) THEN
              WRITE(6,*) ' W times Correction vector on LU7 '
              CALL WRTVCD(VEC2,LU7,1,LBLK)
            END IF
*. - <0!W!0> !N,M-1,LW>
C           IF(W00.NE.0.0D0 .AND. N+M .NE. 1) THEN
            IF(W00.NE.0.0D0 ) THEN
              ONE = 1.0D0
              FACTOR = -W00
              CALL REWINE(LU7,-1)
              CALL REWINE(LU6,-1)
              CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU7,LU6,LU4,0,LBLK)
            ELSE
C?            write(6,*) ' route 1 '
              CALL REWINE(LU7,-1)
              CALL COPVCDP(LU7,LU4,VEC1,0,LBLK)
            END IF
            RHSNORM = INPRDD(VEC1,VEC2,LU7,LU7,1,LBLK)
C?          WRITE(6,*) ' Norm of vector on LU7  ', RHSNORM
          END IF
          IF(NTEST.GE.100) THEN
            write(6,*) ' v,w times vectors calculated '
          END IF
*. Add the above vectors and save on LU6
          NVEC = IMIN + IPLUS + INOT
C?        write(6,*) ' imin iplus, inot', imin,iplus,inot
C?          RHSNORM = INPRDD(VEC1,VEC2,LU4,LU4,1,LBLK)
C?          WRITE(6,*) ' Norm of first vector on lu4  ', RHSNORM
          DO IVEC = 1, NVEC
            SCR(IVEC) = 1.0D0
          END DO
          ZERO = 0.0D0
          CALL MVCSMD2(LU4,SCR,ZERO,LU6,LU7,VEC1,VEC2,
     &                 NVEC,1,LBLK)
          IF(NTEST.GE.100) THEN
           WRITE(6,*) ' rhs before sum V01D* D '
           CALL WRTVCD(VEC1,LU6,1,LBLK)
          END IF
C?        RHSNORM = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
C?        WRITE(6,*) ' Norm of RHS before sum V01D*D ', RHSNORM
*
*       -sum(ld=1,n;ls=1,m-1;llw)Vo[1]t d(ld,ls,llw)d(n-ld,m-ls-1,lw-llw)
*
*. d(ld,ls,lw) must be of the same symmety as Vo=W
*. d(n-ld,m-ls-1,(k-l)w) will then have the same symmetry as d(n,m,*)
*
* Add up on LU6, using LU7 as scratch
          ID1SM  = IWRFSM
          ID2SM = MNSM
*
C?        write(6,*) ' TERM 1'
*. Avoid occurance of zero order vector
          DO LD = 0, N
            IF(LD.EQ.0) THEN
              LSMIN = 1
            ELSE
              LSMIN = 0
            END IF
            IF(LD.EQ.N) THEN
              LSMAX = M-2
            ELSE
              LSMAX = M-1
            END IF
            DO LS = LSMIN, LSMAX
              DO LLW = -LD,LD, 2
*. within frequency bounds ?
              IF(ABS(LW-LLW).LE.N-LD) THEN
                IF(NTEST.GE.5) THEN
                  WRITE(6,*) ' Within frequency bounds for LD,LS,LLW=',
     &            LD,LS,LLW
                END IF
                ID1PSM  = ISMVW(MOD(LD,2)+1,MOD(LS,2)+1)
                ID2PSM  = ISMVW(MOD(N-LD,2)+1,MOD(M-LS-1,2)+1)
                IF(ID1PSM.EQ.ID1SM.AND.ID2PSM.EQ.ID2SM) THEN
                  IF(NTEST.GE.5) THEN
                    WRITE(6,*)
     &              ' term (Vo[1]t d(ld,ls,llw))d(n-ld,m-ls-1,lw-llw)',
     &              ' with right symmetry'
                    WRITE(6,*)
     &              'LD LS',LD,LS
                  END IF
*. Position LUN for reading D(n-ld,m-ls-1,(k-l)w)
                  CALL GET_NMLW(VEC2,N-LD,M-LS-1,LW-LLW,
     &                 MAXV,MAXW,LUN,LURF,0,LU7,1,LBLK,LUOUTEFF)
CW                IF(LD.EQ.N-LD.AND.LS.EQ.M-LS-1.AND.LLW.EQ.LW-LLW)
CW   &            THEN
CW                  COEF = -2.0D0*V01D(LD,LS,LLW)
CW                ELSE
                    COEF = -V01D(LD,LS,LLW)
CW                END IF
C?                write(6,*) ' COEF = ',coef
                  ONE = 1.0D0
                  CALL REWINE(LU6,-1)
                  CALL REWINE(LU7,-1)
                  CALL VECSMDP(VEC1,VEC2,ONE,COEF,LU6,LUOUTEFF,
     &            LU7,0,LBLK)
                  CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
                END IF
              END IF
              END DO
            END DO
          END DO
*
*       -sum(ld=1,n-1;ls=1,m;llw)Vw[1]t d(ld,ls,lw)d(n-ld-1,m-ls,lw-llw-1)
*
C?        write(6,*) ' TERM 2 '
          ID1SM  = IVRFSM
          DO LD = 0, N-1
            IF(LD.EQ.0) THEN
              LSMIN = 1
            ELSE
              LSMIN = 0
            END IF
            IF(LD.EQ.N-1) THEN
              LSMAX = M-1
            ELSE
              LSMAX = M
            END IF
            DO LS = LSMIN, LSMAX
              DO LLW = -LD,LD,2
*. within frequency bounds ?
              IF(ABS(LW-LLW-1).LE.N-LD-1) THEN
                IF(NTEST.GE.5) THEN
                  WRITE(6,*) ' Within frequency bounds for LD,LS,LLW=',
     &            LD,LS,LLW
                END IF
                ID1PSM  = ISMVW(MOD(LD,2)+1,MOD(LS,2)+1)
                ID2PSM  = ISMVW(MOD(N-LD-1,2)+1,MOD(M-LS,2)+1)
                IF(ID1PSM.EQ.ID1SM.AND.ID2PSM.EQ.ID2SM) THEN
                  IF(NTEST.GE.5) THEN
                    WRITE(6,*)
     &              ' term (Vw[1]t d(ld,ls,lw))d(n-ld-1,m-ls,lw-llw-1)',
     &              ' with right symmetry'
                    WRITE(6,*)
     &              'LD LS',LD,LS
                  END IF
*. Position LUN at start of D(n-ld-1,m-ls,(k-l)w)
                  CALL GET_NMLW(VEC2,N-LD-1,M-LS,LW-LLW-1,
     &                 MAXV,MAXW,LUN,LURF,0,LU7,0,LBLK,LUOUTEFF)
CW                IF(LD.EQ.N-LD-1.AND.LS.EQ.M-LS.AND.LLW.EQ.LW-LLW-1)
CW   &            THEN
CW                  COEF = -2.0D0*VW1D(LD,LS,LLW)
CW                ELSE
                    COEF = -VW1D(LD,LS,LLW)
CW                END IF
C?                write(6,*) ' COEF = ',coef
                  ONE = 1.0D0
                  CALL REWINE(LU6,-1)
                  CALL REWINE(LU7,-1)
                  CALL VECSMDP(VEC1,VEC2,ONE,COEF,LU6,
     &                 LUOUTEFF,LU7,0,LBLK)
                  CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
                END IF
              END IF
              END DO
            END DO
          END DO
*
*       -sum(ld=1,n-1;ls=1,m;ll)(V-w[1]t d(ld,ls,lw))d(n-ld-1,m-ls,(l-ll+1)w)
*
C?        write(6,*) ' TERM 3 '
          DO LD = 0, N-1
            IF(LD.EQ.0) THEN
              LSMIN = 1
            ELSE
              LSMIN = 0
            END IF
            IF(LD.EQ.N-1) THEN
              LSMAX = M-1
            ELSE
              LSMAX = M
            END IF
            DO LS = LSMIN, LSMAX
              DO LLW = -LD,LD,2
*. within frequency bounds ?
              IF(ABS(LW-LLW+1).LE.N-LD-1) THEN
                IF(NTEST.GE.5) THEN
                  WRITE(6,*) ' Within frequency bounds for LD,LS,LLW=',
     &            LD,LS,LLW
                END IF
                ID1PSM  = ISMVW(MOD(LD,2)+1,MOD(LS,2)+1)
                ID2PSM  = ISMVW(MOD(N-LD-1,2)+1,MOD(M-LS,2)+1)
                IF(ID1PSM.EQ.ID1SM.AND.ID2PSM.EQ.ID2SM) THEN
                  IF(NTEST.GE.5) THEN
                    WRITE(6,*)
     &              'term (V-w[1]t d(ld,ls,lw))d(n-ld-1,m-ls,lw-llw+1)'
     &              ,' with right symmetry'
                    WRITE(6,*)
     &              'LD LS',LD,LS
                  END IF
*. Position LUN for reading D(n-ld-1,m-ls,(k-l)w)
                  CALL GET_NMLW(VEC2,N-LD-1,M-LS,LW-LLW+1,
     &                 MAXV,MAXW,LUN,LURF,0,LU7,0,LBLK,LUOUTEFF)
*. We assume that the perturbation is real,
                  IF(IRC.EQ.1) THEN
                    COEF = -VW1D(LD,LS,LLW)
                  ELSE
                    COEF =  VW1D(LD,LS,LLW)
                  END IF
CW                IF(LD.EQ.N-LD-1.AND.LS.EQ.M-LS.AND.LLW.EQ.LW-LLW+1)
CW   &            THEN
CW                  COEF = 2.0D0*COEF
CW                END IF
C?                write(6,*) ' COEF = ',coef
                  ONE = 1.0D0
                  CALL REWINE(LU6,-1)
                  CALL REWINE(LU7,-1)
                  CALL VECSMDP(VEC1,VEC2,ONE,COEF,LU6,
     &                 LUOUTEFF,LU7,0,LBLK)
                  CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
                END IF
              END IF
              END DO
            END DO
          END DO
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' After the 3 terms '
            CALL WRTVCD(VEC1,LU6,1,LBLK)
          END IF
*. We have now assambled the right hand side vector
*. Well there was a minus in front of (e[2]-wI)**-1
          ONEM = -1.0D0
          CALL SCLVCD(LU6,LU7,ONEM,VEC1,1,LBLK)
CSEPT29   CALL COPVCDP(LU7,LU6,VEC ,1,LBLK)
          CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
          IF(NTEST.GE.5) THEN
             WRITE(6,*) ' Assembling right hand site vector finished'
          END IF
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Right hand side vector before projection '
            CALL WRTVCD(VEC1,LU6,1,LBLK)
          END IF
*. Symmetry of right hand side is also ISSM
          ICSM = ISSM
C?        RHSNORM = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
C?        WRITE(6,*) ' Norm of RHS 1 ', RHSNORM
*. project !0> component out, save result on LU6
          IF(ICSM.EQ.IRFSM) THEN
            OVLAP = INPRDD(VEC1,VEC2,LURF,LU6,1,LBLK)
C?          write(6,*) ' ovlap ', OVLAP
            FACTOR = -OVLAP
            ONE = 1.0D0
            CALL VECSMDP(VEC1,VEC2,ONE,FACTOR,LU6,LURF,LU7,1,LBLK)
            CALL COPVCDP(LU7,LU6,VEC1,1,LBLK)
          END IF
C?        RHSNORM = INPRDD(VEC1,VEC2,LU6,LU6,1,LBLK)
C?        WRITE(6,*) ' Norm of RHS ', RHSNORM
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Right hand side vector after projection '
            CALL WRTVCD(VEC1,LU6,1,LBLK)
          END IF
*.save in LU1 and solve linear eq .
          CALL REWINE(LU1,-1)
          CALL COPVCDP(LU6,LU1,VEC1,1,LBLK)
*. Should !0> be printed out in each micro
          IF(ICSM.EQ.IRFSM) THEN
            LUPROJ = LURF
            IPROJ = 1
          ELSE
            LUPROJ = 0
            IPROJ = 0
          END IF
*. Change core energy to include -lw *FREQ
          SHIFT = -ERF + LW*FREQ
C?        WRITE(6,*) ' Shift in GNDB.. ', SHIFT
          I12 = 2
C?        WRITE(6,*) ' LUPROJ and SHIFT', LUPROJ,SHIFT
          CALL HINTV(LU1,LU2,SHIFT,SHIFT,VEC1,VEC2,LBLK,LUPROJ)
*. Solution to equations are hiding on LU2, transfer to LUN
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' NEW Correction vector'
            CALL WRTVCD(VEC1,LU2,1,LBLK)
          END IF
*. Add to our collection of correction vectors
          CALL GET_NMLW(VEC1,N,M,LW,
     &         MAXV,MAXW,LUN,LURF,0,LU7,0,LBLK,LUOUTEFF)
          CALL REWINE(LU2,-1)
          CALL COPVCDP(LU2,LUN,VEC1,0,LBLK)
*. Test overlap with reference
          IF(ICSM.EQ.IRFSM) THEN
            OVLAP = INPRDD(VEC1,VEC2,LURF,LU2,1,LBLK)
C?          write(6,*) ' overlap between reference and correction ',
C?   &      IORD,OVLAP
          END IF
*. Obtain the V[1]t D(n,m,lw) arrays
*
          CALL REWINE(LU5,-1)
*
C?        write(6,*) ' MNSM IVRFSM IWRFSM ', MNSM,IVRFSM,IWRFSM
          IF(MNSM.EQ.IVRFSM) THEN
            CALL REWINE(LU2,-1)
            VW1D(N,M,LW) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
          ELSE
            CALL SKPVCD(LU5,1,VEC1,0,LBLK)
            VW1D(N,M,LW) = 0.0D0
          END IF
          IF(NTEST.GE.5)
     &    write(6,*) '  VW1D(N,M,LW) =  ', VW1D(N,M,LW)
*
          IF(MNSM.EQ.IWRFSM) THEN
            CALL REWINE(LU2,-1)
            V01D(N,M,LW) = INPRDD(VEC1,VEC2,LU2,LU5,0,LBLK)
          ELSE
            CALL SKPVCD(LU5,1,VEC1,0,LBLK)
            V01D(N,M,LW) = 0.0D0
          END IF
          IF(NTEST.GE.5)
     &    write(6,*) ' V01D(N,M,LW) =  ', V01D(N,M,LW)
*
         END DO
*        ^ End of loop over frequency components
        END DO
*       ^ End of loop over W
       END DO
*       ^ End of loop over V
*
       IF(NTEST.GE.5) THEN
*
       WRITE(6,*)
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' V0[1]t D(ld,ls,lw) : '
       WRITE(6,*) ' ==================== '
       WRITE(6,*)
       DO LD = 0, MAXV
         IF(LD.EQ.0) THEN
           LSMIN = 1
         ELSE
           LSMIN = 0
         END IF
         DO LS = LSMIN, MAXW
           DO LW = -LD,LD,2
             WRITE(6,'(A,3I3,A,E18.10)')
     &       ' ld, ls, lw = ',LD,LS,LW,'   ,   ',V01D(LD,LS,LW)
           END DO
         END DO
       END DO
*
       WRITE(6,*)
       WRITE(6,*) ' ==================== '
       WRITE(6,*) ' Vw[1]t D(ld,ls,lw) : '
       WRITE(6,*) ' ==================== '
       WRITE(6,*)
       DO LD = 0, MAXV
         IF(LD.EQ.0) THEN
           LSMIN = 1
         ELSE
           LSMIN = 0
         END IF
         DO LS = LSMIN, MAXW
           DO LW = -LD,LD,2
             WRITE(6,'(A,3I3,A,E18.10)')
     &       ' ld, ls, lw = ',LD,LS,LW,'   ,   ',Vw1D(LD,LS,LW)
           END DO
         END DO
       END DO
*
       END IF
*      ^ End of print section
*
*. And then the susceptibilities
*
* The n'th ,m'th order susceptibility with w the operator to be studied
*  is <<w;v,v...,w,w..>>f1,f2,f3,...fn,0,0,0,
* where f1,f2,f3 ... fn equals +/- freq. It depends only upon the total
* frequency  ftot= sum(i=1,n)fi and can be written as
*
* Khi(ftot,f,f,f,f,f...,0,0,0,...)w,v,v,v,....,w,w,...
*
*. The susceptibility is the term on the Taylor expansion of
* <0!w!0> with the right frequencies and perturbations
* and by insisting on symmetry on all indeces after the
*
*. Inner products <d!d>
*
*. Currently we only asemble up to MAXV,MAXW
      DO LD = 0,MAXV
        DO LS = 0,MAXW
          DO LW = -LD,LD,2
            IF(NTEST.GE.3) THEN
              WRITE(6,*)
              write(6,*) ' Terms to DD for LD LS LW =',LD,LS,LW
              WRITE(6,*)
            END IF
            DD(LD,LS,LW) = 0.0D0
*. Obtain <d|d>(ld,ls,lw) =
*  sum(ldl,lsl,lwl)dt(ldl,lsl,-lwl)d(ld-ldl,ls-lsl,lw-lwl)
            DO LDL = 0, LD
              DO LSL = 0, LS
                LDR = LD - LDL
                LSR = LS - LSL
*. Identical symmetries
                ILSM = ISMVW(MOD(LDL,2)+1,MOD(LSL,2)+1)
                IRSM = ISMVW(MOD(LDR,2)+1,MOD(LSR,2)+1)
                IF(LDR.LE.MAXV.AND.LSR.LE.MAXW.AND.ILSM.EQ.IRSM) THEN
                 DO LWL = -LDL,LDL,2
                  LWR  = LW -LWL
*. Other frequency component allowed ?
                  IF(ABS(LWR).LE.LDR) THEN
*. Fetch D(LDL,LSL,-LWL) and save on LU7
                    CALL GET_NMLW(VEC2,LDL,LSL,-LWL,
     &                   MAXV,MAXW,LUN,LURF,1,LU7,1,LBLK,LUOUTEFF)
*. Position at start of D(LDR,LSR,LWR)
                    CALL GET_NMLW(VEC1,LDR,LSR, LWR,
     &                   MAXV,MAXW,LUN,LURF,0,LU7,0,LBLK,LUOUTEFF)
                    CALL REWINE(LU7,-1)
                    CONT = INPRDD(VEC1,VEC2,LUOUTEFF,LU7,0,LBLK)
*
                    IF(NTEST.GE.3) THEN
                      WRITE(6,*) ' D(ldl,lsl,-lwl)D(ldr,lsr,lwr) for '
                      WRITE(6,*) '  ldl lsl lwl ldr lsr lwr :',
     &                LDL,LSL,LWL,LDR,LSR,LWR
                      WRITE(6,*) '    is ',CONT
                    END IF
                    DD(LD,LS,LW) = DD(LD,LS,LW)+CONT
                  END IF
                 END DO
                END IF
*               ^ End of check of correct symmetries
              END DO
*             ^ End of loop over LSL
            END DO
*           ^ End of loop over LSR
          IF(NTEST.GE.3) WRITE(6,*) ' DD(LD,LS,LW)= ',DD(LS,LS,LW)
          END DO
*         ^ End of loop over LW
        END DO
*       ^ End of loop over LS
      END DO
*     ^ End of loop over LD
      WRITE(6,*)
      WRITE(6,*) ' ==================== '
      WRITE(6,*) '    DD(ld,ls,lw) : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      DO LD = 0, MAXV
        IF(LD.EQ.0) THEN
          LSMIN = 1
        ELSE
          LSMIN = 0
        END IF
        DO LS = LSMIN, MAXW
          DO LW = -LD,LD,2
            WRITE(6,'(A,3I3,A,E18.10)')
     &      ' ld, ls, lw = ',LD,LS,LW,'   ,   ',DD(LD,LS,LW)
          END DO
        END DO
      END DO
*. We will now calculate elements <d!A!d> and the
*  corresponding response functions.
* This involves a loop over A- operators so :
      DO I_AVE_OP = 1, N_AVE_OP
        WRITE(6,*)
        WRITE(6,'(A)')
     &  ' *****************************************'
        WRITE(6,'(A,A)')
     &  ' Info for A - operator ', AVE_OP(I_AVE_OP)
        WRITE(6,'(A)')
     &  ' *****************************************'
        WRITE(6,*)
*
        IASM = IAVE_SYM(I_AVE_OP)
C?      WRITE(6,*) ' IASM = ', IASM
*
*. Inner products <d!A!d>
*
* Notation A => W
      DO LD = 0,MAXV
        DO LS = 0,MAXW

          DO LW = -LD,LD,2
            DWD(LD,LS,LW) = 0.0D0
*. Obtain <d|w!d>(ld,ls,lw) =
*  sum(ldl,lsl,lwl)dt(ldl,lsl,-lwl)d(ld-ldl,ls-lsl,lw-lwl)
            IF(NTEST.GE.5) THEN
              WRITE(6,*)
     &        ' DAD under construction for LD,LS,LW=',LD,LS,LW
            END IF
*
            DO LDL = 0, LD
              DO LSL = 0, LS
                LDR = LD - LDL
                LSR = LS - LSL
*. Identical symmetries
                ILSM = ISMVW(MOD(LDL,2)+1,MOD(LSL,2)+1)
                IRSM = ISMVW(MOD(LDR,2)+1,MOD(LSR,2)+1)
*W*right
C               IWRSM = ISMVW(MOD(LDR,2)+1,MOD(LSR+1,2)+1)
                IWRSM = MULTD2H(IRSM,IASM)
C?              WRITE(6,*) 'IWRSM = ', IWRSM
                IF(LDR.LE.MAXV.AND.LSR.LE.MAXW.AND.ILSM.EQ.IWRSM) THEN
                 DO LWL = -LDL,LDL,2
                  LWR  = LW -LWL
*. Other frequency component allowed ?
                  IF(ABS(LWR).LE.LDR) THEN
*. Fetch D(LDR,LSR,LWR) on LU7
                    CALL GET_NMLW(VEC1,LDR,LSR, LWR,
     &                   MAXV,MAXW,LUN,LURF,1,LU7,1,LBLK,LUOUTEFF)
*. A time D(LDR,LSR,LWR) on LU6
                    ICSM = ISMVW(MOD(LDR,2)+1,MOD(LSR,2)+1)
C                   ISSM= ISMVW(MOD(LDR,2)+1,MOD(LSR+1,2)+1)
                    ISSM = MULTD2H(ICSM,IASM)
C?                  WRITE(6,*) ' ISSM = ', ISSM
*
                    I12 = 1
                    CALL BVEC(AVINT(1,I_AVE_OP),IASM,LU7,LU6,
     &                        VEC1,VEC2)
C                        BVEC(B,IBSM,LUC,LUB,VEC1,VEC2)
*. Fetch D(LDL,LSL,-LWL) on LU7
                    CALL GET_NMLW(VEC1,LDL,LSL,-LWL,
     &                   MAXV,MAXW,LUN,LURF,1,LU7,1,LBLK,LUOUTEFF)
                    CONT = INPRDD(VEC1,VEC2,LU7,LU6,1,LBLK)
                    IF(NTEST.GE.3) THEN
                     WRITE(6,*) ' Contribution to DAD : '
                     WRITE(6,*) '    LDL LSL LWL',LDL,LSL,LWL
                     WRITE(6,*) '    LDR LSR LWR',LDR,LSR,LWR
                     WRITE(6,*) '    CONT = ',CONT
                    END IF
                    DWD(LD,LS,LW) = DWD(LD,LS,LW)+CONT
                  END IF
                 END DO
                END IF
*               ^ End of check of correct symmetries
              END DO
*             ^ End of loop over LSL
            END DO
*           ^ End of loop over LSR
          END DO
*         ^ Enf of loop over LW
        END DO
*       ^ End of loop over LS
      END DO
*     ^ End of loop over LD
      WRITE(6,*)
      WRITE(6,*) ' ==================== '
      WRITE(6,*) '    DAD(ld,ls,lw) : '
      WRITE(6,*) ' ==================== '
      WRITE(6,*)
      DO LD = 0, MAXV
        IF(LD.EQ.0) THEN
          LSMIN = 1
        ELSE
          LSMIN = 0
        END IF
        DO LS = LSMIN, MAXW
          DO LW = -LD,LD,2
            WRITE(6,'(A,3I3,A,E18.10)')
     &      ' ld, ls, lw = ',LD,LS,LW,'   ,   ',DWD(LD,LS,LW)
          END DO
        END DO
      END DO
*
* **********************************
* . And then : The susceptibilities
* **********************************
*
* The khi's are to obtained as expansion of (<0!A!0>/<0!0>)
*
* It is slightly inconvenient to expand 1/(<0!0>) to arbitary order.
*. Instead we multiply both sides by <0!0>
*  Khi*<0!0> = <0!A!0> and expand this
*
      DO LD = 0, MAXV
       DO LS = 0, MAXW
         DO LW = -LD,LD,2
         IF(NTEST.GE.3) THEN
           WRITE(6,*) ' Susceptibility for LD,LS,LW =',LD,LS,LW
           WRITE(6,*) ' ========================================='
           WRITE(6,*)
         END IF
*. Obtain chi(ld,ls,-lw) =
*<0!A!0>(ld,ls,lw) -
*sum(ldp,lsp,lwp,ldp+lsp.ge.0)chi(ldp,lsp,-lwp)<0!0>(ld-ldp,ls-lsp,lw-lwp)
         CHI(LD,LS,-LW) = DWD(LD,LS,LW)
         IF(NTEST.GE.3)THEN
            WRITE(6,*)
     &      ' Initial term <0!A!0>(ld,ls,lw) ',CHI(LD,LS,-LW)
            WRITE(6,'(A)')
     &      '  Terms -chi(ldp,lsp,-lwp) <0!0>(ld-ldp,ls-lsp,lw-lwp) : '
         END IF
         DO LDP = 0, LD
*
          IF(LDP.EQ.LD) THEN
           LSPMAX = LS-1
          ELSE
           LSPMAX = LS
          END IF
*
          DO LSP = 0,LSPMAX
           DO LWP = -LDP,LDP
*. Allowed frequency splitting
            IF(ABS(LW-LWP).LE.(LD-LDP)) THEN
              TERM =- CHI(LDP,LSP,-LWP)*DD(LD-LDP,LS-LSP,LW-LWP)
              CHI(LD,LS,-LW) = CHI(LD,LS,-LW) + TERM
              IF(NTEST.GE.3) THEN
                WRITE(6,*)
     &          '   ldp lsp -lwp = ', LDP,LSP,-LWP ,' is ', term
              END IF
            END IF
           END DO
          END DO
*         ^ End of LSP,LSW loops
         END DO
*.       ^ End of LDP loop, end thereby end of terms for given susceptibility
         IF(NTEST.GE.3) THEN
         WRITE(6,*) ' Unscaled Susceptibility for LD,LS,-LW = ',
     &              LD,LS,-LW,  ' is ',   CHI(LD,LS,-LW)
         END IF
*
         END DO
*        ^ End of loop over LW
       END DO
*      ^ End of loop over LS
      END DO
*     ^ End of loop over LD
*
*
*. What we have calculated now is the term in a powerseries expansion
*  of given orders and frequency. Three scalings must be performed
* 1 : susceptibilities are taylor coefficients, multiply with ( ld+ls)!
* 2 : Susceptibilities like <<A;Vd,Vs>> and <<A;Vs,Vd>> are identical
*     and only their sum was obtained. Divide by number of different
*     components to obtain individual susceptibilities
*     but separate.
* 3 : There can be several frequence components underlying a givrn total
*     frequence, f.ex, second order zero freq corresponds to +w,-w
*     and -w +w. These two should be considered separate susceptibilities.
       WRITE(6,*)
       WRITE(6,*) ' ================================= '
       WRITE(6,*) ' Response functions Chi(ld,ls,lw) : '
       WRITE(6,*) ' ================================= '
       WRITE(6,*)
       WRITE(6,*)
       DO LD = 0, MAXV
         IF(LD.EQ.0) THEN
           LSMIN = 1
         ELSE
           LSMIN = 0
         END IF
         DO LS = LSMIN, MAXW
           DO LW = -LD,LD,2
*
             FAC1 = IFAC(LD+LS)
*. Number of ld, ls combinations
             FAC2 = IBION(LD+LS,LD)
*. Number of components with frequency up
             NUP = (LD+LW)/2
             NCOMB = IBION(LD,NUP)
             FAC = FAC1/FLOAT(NCOMB)/FAC2
C?           WRITE(6,*) ' FAC1, FAC2, NCOMB ',FAC1,FAC2, NCOMB
             CHI(LD,LS,LW) = FAC * CHI(LD,LS,LW)
*
             WRITE(6,'(A,3I3,A,E18.10)')
     &       ' ld, ls, lw = ',LD,LS,LW,'   ,   ',CHI(LD,LS,LW)
           END DO
         END DO
       END DO
*
       END DO
*      ^ End of loop over A-operators
*
      WRITE(6,*) ' Returning from GNDBPTFREQ '
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      FUNCTION NDIM_1EL_MAT(IHSM,NRPSM,NCPSM,NSM,IPACK)
*
* Length of one-electron matrix with symmetry IHSM
*
* Jeppe Olsen, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION NRPSM(*),NCPSM(*)
*
      LENGTH = 0
      DO IRSM = 1, NSM
        CALL SYMCOM(2,5,IRSM,ICSM,IHSM)
C       WRITE(6,*) ' IHSM,IRSM,ICSM', IHSM,IRSM,ICSM
        IF(IPACK.EQ.0.OR.(IPACK.EQ.1.AND.IRSM.GT.ICSM)) THEN
          LENGTH = LENGTH + NRPSM(IRSM)*NCPSM(ICSM)
        ELSE IF(IPACK.EQ.1.AND.IRSM.EQ.ICSM) THEN
          LENGTH = LENGTH + NRPSM(IRSM)*(NRPSM(IRSM)+1)/2
        END IF
      END DO
*
      NDIM_1EL_MAT = LENGTH
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' 2 Dim array , sym and length : ', IHSM,LENGTH
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE ONE_EL_PROP
*
* Calculate one-electron properties
* One-electron density is assumed in WORK(KRHO1)
*
* Jeppe Olsen, June 1997
*              Updated Feb. 98 ( Natural orbital analysis added )
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLOCC, KLHDIA, KLHNAT, KLXNAT, KLSCR, KLRHO1S, 
     &          KLCMO, KLHONE
!               for addressing of WORK
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "multd2h.inc"
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
*
      REAL*8 INPROD
      CHARACTER*1 XYZ(3)
      CHARACTER*8 LABEL
      CHARACTER*6 LABEL2
      DATA XYZ/'X','Y','Z'/
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'ONE_EL')
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================'
        WRITE(6,*) ' Welcome to ONE_EL_PROP'
        WRITE(6,*) ' ======================'
        WRITE(6,*)
      END IF
*. A bit of local memory
      LHONE = NTOOB * NTOOB
      WRITE(6,*) ' Max size of one-electron operator',LHONE
      CALL MEMMAN(KLHONE,LHONE,'ADDL  ',2,'KLHONE')
      CALL MEMMAN(KLCMO,LHONE,'ADDL  ',2,'KLCMO ')
      CALL MEMMAN(KLRHO1S,LHONE,'ADDL  ',2,'LRHO1S')
      CALL MEMMAN(KLSCR,4*LHONE,'ADDL  ',2,'KLCMO ')
*. Natural orbital expansion
      CALL MEMMAN(KLXNAT,LHONE,'ADDL  ',2,'XNAT  ')
*. Integrals in natural orbital basis
      CALL MEMMAN(KLHNAT,LHONE,'ADDL  ',2,'HNAT  ')
*. Diagonal of integrals in natorb basis
      CALL MEMMAN(KLHDIA,NTOOB,'ADDL  ',2,'HDIA  ')
*. Occupation numbers
      CALL MEMMAN(KLOCC,NTOOB,'ADDL  ',2,'OCCNUM')
*
*
*. Assumed symmetry for one-electron density-
*. for making  change to general transitions densities later
*
      IRHO1SM = 1
*. Extract symmetry blocks from complete one-electron density
C     REORHO1(RHO1I,RHO1O,IRHO1SM)
      CALL REORHO1(WORK(KRHO1),WORK(KLRHO1S),IRHO1SM)
*. Number of elements in symmetry blocks of integrals and density
      LRHO1S = 0
      DO ISM = 1, NSMOB
        JSM = MULTD2H(ISM,IRHO1SM)
        LRHO1S = LRHO1S + NTOOBS(ISM)*NTOOBS(JSM)
      END DO
      WRITE(6,*) ' Number of elements in symmetry blocks ',
     &LRHO1S
*. Natural orbitals
      CALL LNATORB(WORK(KRHO1),NSMOB,NTOOBS,NACOBS,NINOBS,
     &            IREOST,WORK(KLXNAT),
     &            WORK(KLHNAT),WORK(KLOCC),NACOB,
     &            WORK(KLSCR),IPRDEN)
      WRITE(6,*) ' memchk 1 '
      CALL LMEMCHK('after LNATORB')
*
      DO IPROP =1, NPROP
        WRITE(6,'(A,A)') ' Property to be calculated',
     &                     PROPER(IPROP)
*. one- or two-dimensional tensor ?
        LABEL2 = PROPER(IPROP)
C            GET_PROP_RANK(PROPER,IRANK)
        WRITE(6,'(A,A)') ' LABEL2 =', LABEL2
        CALL GET_PROP_RANK(LABEL2,NRANK)
C       IF(LABEL2.EQ.'DIPOLE') THEN
C         NRANK = 1
C       ELSE IF(LABEL2.EQ.'THETA ' .OR.
C    &          LABEL2.EQ.'QUADRU' .OR.
C    &          LABEL2(1:3).EQ.'EFG' ) THEN
C         NRANK  = 2
C       ELSE
        IF(NRANK.EQ.-1) THEN
         WRITE(6,'(A,A)') ' Unknown operator ',PROPER(IPROP)
         NRANK = 0
        END IF
        WRITE(6,*) ' Rank of operator ', NRANK
*.
        IF(NRANK.EQ.1) THEN
          DO ICOMP = 1, 3
            IF(IXYZSYM(ICOMP).EQ.IRHO1SM) THEN
              WRITE(6,*) ' right symmetry for component',ICOMP
*. Label of integrals on file- currently DALTON FORM !!
              IF(PROPER(IPROP).EQ.'DIPOLE') THEN
                LABEL =XYZ(ICOMP)//'DIPLEN '
              END IF
              WRITE(6,'(A,A)') ' Label ', LABEL
*. Obtain one-electron integrals
              CALL GET_PROPINT(WORK(KLHONE),IRHO1SM,LABEL,
     &             WORK(KLSCR),NTOOBS,NTOOBS,NSMOB,0)
*. and then : Expectation value
              EXPVAL = INPROD(WORK(KLRHO1S),WORK(KLHONE),LRHO1S)
                   WRITE(6,'(A,A,E22.14)')
     &             ' Expectation value of ',LABEL,  EXPVAL
*. Analysis in terms of natural orbitals
*. Transform  integrals to nat orb basis
              CALL TRAN_SYM_BLOC_MAT2(WORK(KLHONE),WORK(KLXNAT),
     &             NSMOB,NTOOBS,WORK(KLHNAT),WORK(KLSCR),0)
*. Extract diagonal integrals
              CALL GET_DIAG_BLOC_MAT(WORK(KLHNAT),WORK(KLHDIA),
     &             NSMOB,NTOOBS,0)
              CALL PROP_NATORB(WORK(KLHDIA),WORK(KLOCC),NTOOBS,NSMOB)
            END IF
          END DO
        ELSE IF(NRANK.EQ.2) THEN
          DO ICOMP = 1,3
            DO JCOMP = 1,ICOMP
              IF(MULTD2H(IXYZSYM(ICOMP),IXYZSYM(JCOMP))
     &        .EQ.IRHO1SM) THEN
                WRITE(6,*) ' Right symmetry for components',
     &          ICOMP,JCOMP
*
C               IF(LABEL2.EQ.'THETA ') THEN
*. Buckinghams traceless quadrupole moment
C                  LABEL=XYZ(JCOMP)//XYZ(ICOMP)//'THETA'
C               ELSE IF(LABEL2.EQ.'QUADRU') THEN
C                 LABEL=XYZ(JCOMP)//XYZ(ICOMP)//'QUADRU'
C               ELSE IF (LABEL2(1:3).EQ.'EFG' ) THEN
*. Again DALTON format, we assume a efg component for a specific
* nuclei so the label has the form EFGabc
                LABEL=XYZ(JCOMP)//XYZ(ICOMP)//LABEL2
C               END IF
*. Obtain one-electron integrals
                WRITE(6,*) ' memchk 2'
                CALL LMEMCHK('ONE_EL_PROP before GET_PROPINT')
                CALL GET_PROPINT(WORK(KLHONE),IRHO1SM,LABEL,
     &               WORK(KLSCR),NTOOBS,NTOOBS,NSMOB,0)
*. and then : Expectation value
                EXPVAL = INPROD(WORK(KLRHO1S),WORK(KLHONE),LRHO1S)
                   WRITE(6,'(A,A,E22.14)')
     &             ' Expectation value of ',LABEL,  EXPVAL
*. Analysis in terms of natural orbitals
*. Transform  integrals to nat orb basis
                CALL TRAN_SYM_BLOC_MAT2(WORK(KLHONE),WORK(KLXNAT),
     &               NSMOB,NTOOBS,WORK(KLHNAT),WORK(KLSCR),0)
                WRITE(6,*) ' memchk 3'
                CALL LMEMCHK('ONE_EL_PROP after TRAN_SYM_BLOC_MAT2')
*. Extract diagonal integrals
                CALL GET_DIAG_BLOC_MAT(WORK(KLHNAT),WORK(KLHDIA),
     &               NSMOB,NTOOBS,0)
                CALL PROP_NATORB(WORK(KLHDIA),WORK(KLOCC),
     &               NTOOBS,NSMOB)
                WRITE(6,*) ' memchk 4'
                CALL LMEMCHK('ONE_EL_PROP after PROP_NATORB')
              END IF
            END DO
          END DO
        END IF
*
      END DO
*. and then : Expectation value
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'ONE_EL')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
*
      SUBROUTINE PROP_NATORB(HDIAG,OCCNUM,NTOOBS,NSMOB)
*
* Analyze property in natural orbital representation
*
* Jeppe Olsen, Feb 98
*
* Input :
* ========
*
* HDIAG : Diagonal of 1-e integrals over nat orbs symmetry order
* OCCNUM: Occupation numbers
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION HDIAG(*),OCCNUM(*)
*
      DIMENSION NTOOBS(*)
*
      SUM = 0.0D0
      WRITE(6,'(A)')
     & ' Orb   Sym    Integral   Occ number   Contribution    Sum '
      WRITE(6,'(A)')
     & ' ========================================================='
      DO ISYM = 1, NSMOB
      WRITE(6,*)
        IF(ISYM.EQ.1) THEN
          IOFF = 1
        ELSE
          IOFF = IOFF + NTOOBS(ISYM-1)
        END IF
        DO IORB = 1, NTOOBS(ISYM)
          IORBEF = IORB+IOFF-1
          CONT = OCCNUM(IORBEF)*HDIAG(IORBEF)
          SUM = SUM + CONT
          WRITE(6,'(1H ,I4, I4,4E13.5)')
     &    IORB,ISYM,HDIAG(IORBEF),OCCNUM(IORBEF),CONT,SUM
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
      SUBROUTINE PROP_PERT(LU0,LUN,N,ISM,ISPC)
*
* Perturbation expansion of one-electron properties
*
* It is assumed that this calculation is preceded  by
* a call to the perturbation routine to obtain the
* wave function corrections to the neutral state.
*
* Input
*       LUN : File containing wave function corrections
*       LU0 : File containing reference wave funcrtion
*         N : Max order of expansion
*      ISM : Symmetry of reference state
*      ISPC : Space of referencestate
*
*
* Jeppe Olsen, April 98 ( on the train for once )
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*8 KLSI, KLSIJ, KLSCR4, KLFSCR, KLDEN2P, KLDEN2, KLDEN1,
     &          KLDEN1P, KLVEC2, KLVEC1, IOFF, JOFF, III
!               for addressing of WORK
      REAL *8 INPRDD
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "lucinp.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*. Local scratch
      PARAMETER(MXNORD = 100)
      DIMENSION SIJ(MXNORD*(MXNORD+1)/2)
      DIMENSION SI(MXNORD)
*
      NTEST = 5
*
      WRITE(6,*)
      WRITE(6,*) ' ============================ '
      WRITE(6,*) ' PROP_PERT is now in CONTROL '
      WRITE(6,*) ' ============================ '
      WRITE(6,*)
      IF(IRELAX.EQ.0) THEN
        WRITE(6,*) ' Property evaluated as expectation value'
      ELSE
        WRITE(6,*) ' Property evaluated as derivative '
      END IF

* a bit on files :
* LUSC36 is LUN.
* Two additional scratch files to be used are  LUSC1 and LUSC2
*
      LBLK = -1
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'PROPPT')
*
*
*     ========================
* 1 : Local memory allocation
*     ========================
*
*. Allocate space for two vector chunks
      CALL MEMMAN(KLVEC1,LBLOCK,'ADDL  ',2,'VEC1  ')
      CALL MEMMAN(KLVEC2,LBLOCK,'ADDL  ',2,'VEC2  ')
* space for one-body Density matrices through order n
      NMAT = N+1
      LENGTH = NMAT * NTOOB ** 2
      CALL MEMMAN(KLDEN1,LENGTH,'ADDL  ',2,'DENN1 ')
*. And an extra set of density matrices
      CALL MEMMAN(KLDEN1P,LENGTH,'ADDL  ',2,'DENN1P')
*. Two-electron densities
      IF(IRELAX.EQ.0) THEN
        KLDEN2 = KLDEN1
      ELSE
        LENGTH = NMAT * NTOOB ** 2 * (NTOOB**2 + 1)/2
        CALL MEMMAN(KLDEN2,LENGTH,'ADDL  ',2,'DENN2 ')
*. And 2e- dens for normalized wf
        CALL MEMMAN(KLDEN2P,LENGTH,'ADDL  ',2,'DENN2P')
      END IF
*. A scratch matrix ( not a nice thing to say about a matrix )
      LENGTH =  2*NINT1
      CALL MEMMAN(KLFSCR,LENGTH,'ADDL  ',2,'FSCR  ')
      CALL MEMMAN(KLSCR4,LENGTH,'ADDL  ',2,'SCR4  ')
*. S(i,j) matrix for overlap of corrections
      CALL MEMMAN(KLSIJ,(N+1)**2,'ADDL  ',2,'KLSIJ ')
      CALL MEMMAN(KLSI , N+1    ,'ADDL  ',2,'KLSI  ')
*
* =========================================================================
*.1 :  overlap of correction vectors ( intermediate normalization  assumed )
* =========================================================================
*
*. Sij(i,j) = <i!j>
      CALL REWINE(LUN,-1)
      DO I = 1, N
*. LUN is positioned at end of vector I-1, copy vector I to LUSC1
        CALL REWINE(LUSC1,-1)
        CALL COPVCD(LUN,LUSC1,WORK(KLVEC1),0,LBLK)
        WRITE(6,*) ' Memcheck I '
        CALL LMEMCHK('PROP_PERT I')
*
        CALL REWINE(LUN,-1)
        DO J = 1, I
          IJ = I*(I-1)/2 + j
          CALL REWINE(LUSC1,-1)
          WORK(KLSIJ-1+IJ) =
     &    INPRDD (WORK(KLVEC1),WORK(KLVEC2),LUSC1,LUN,0,LBLK)
        END DO
        WRITE(6,*) ' Memcheck II '
        CALL LMEMCHK('PROP_PERT II')
      END DO
* SI(i) = sum(j=1,i-1)S(j,i-j)
      DO I = 1, N
        X = 0.0D0
        DO J = 1, I-1
          IMJ = I-J
          IJ = MAX(J,IMJ)*(MAX(J,IMJ)-1)/2+MIN(J,IMJ)
          X = X + WORK(KLSIJ-1+IJ)
        END DO
        WORK(KLSI-1+I) = X
      END DO
*
      IF(NTEST.GE.5) THEN
        WRITE(6,*) ' The S(i,j) Matrix '
        WRITE(6,*) ' ================= '
        CALL PRSYM(WORK(KLSIJ),N)
        WRITE(6,*)
        WRITE(6,*) ' The S(i) array '
        WRITE(6,*) ' ================= '
        CALL WRTMT_LU(WORK(KLSI),N,1,N,1)
      END IF
        WRITE(6,*) ' Memcheck III '
        CALL LMEMCHK('PROP_PERT III')
*
* ===============================================
* 2 : Construct density matrices through order N
* ===============================================
*
      IF(IRELAX.EQ.0) THEN
       ILRHO2 = 0
       LRHO2 = 0
      ELSE
       ILRHO2 = 1
       LRHO2 = NTOOB**2*(NTOOB**2 + 1)/ 2
      END IF
*
      LRHO1 = NTOOB**2
*
*. No print in density matrices
      IPRDEN_SAVE = IPRDEN
      IPRDEN = 0
      I12_SAVE = I12
      IF(IRELAX.EQ.0) THEN
        I12 = 1
      ELSE
        I12 = 2
      END IF
      DO K = 0, N
        CALL PERTDN(K,LU0,LUN,ISM,ISPC,WORK(KLVEC1),WORK(KLVEC2),
     &       WORK(KLDEN1+(K-0)*LRHO1),
     &       WORK(KLDEN2+(K-0)*LRHO2),LUSC1,LUSC2)
      END DO
      IPRDEN = IPRDEN_SAVE
      I12 = I12_SAVE
C     PERTDN
C    &(N,LU0,LUN,ISM,ISPC,VEC1,VEC2,RHO1N,RHO2N,LUSC1,LUSC2,I12)
*
* Change the densities so the correspond to order expansion of
* normalized wf
* Rho'(n) = Rho(n) - sum(j=1,n) Si(j)Rho'(n-j)
*
      ONE = 1.0D0
      DO I = 0, N
        CALL COPVEC(WORK(KLDEN1 +(I-0)*LRHO1),
     &              WORK(KLDEN1P+(I-0)*LRHO1),LRHO1)
        IF(ILRHO2.EQ.1) THEN
          CALL COPVEC(WORK(KLDEN2 +(I-0)*LRHO2),
     &                WORK(KLDEN2P+(I-0)*LRHO2),LRHO2)
        END IF
*
        DO J = 1, I
          FACTOR = -WORK(KLSI-1+J)
          IOFF = KLDEN1P+(I-0)*LRHO1
          JOFF = KLDEN1P+(I-J-0)*LRHO1
          CALL VECSUM(WORK(IOFF),WORK(IOFF),WORK(JOFF),
     &                ONE,FACTOR,LRHO1)
*
          IF(ILRHO2.EQ.1) THEN
            IOFF = KLDEN2P+(I-0)*LRHO2
            JOFF = KLDEN2P+(I-J-0)*LRHO2
            CALL VECSUM(WORK(IOFF),WORK(IOFF),WORK(JOFF),
     &                  ONE,FACTOR,LRHO2)
          END IF
*
        END DO
*
C?      WRITE(6,*) ' Density correction for NORMALIZED wf '
C?      CALL WRTMT_LU(WORK(KLDEN1P+I*LRHO1),NTOOB,NTOOB,NTOOB,NTOOB)
      END DO
*HER ER JEG
*
*. Properties for each order
*
      DO IORD = 0, N
        WRITE(6,*)
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' Information for order ', IORD
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        III = KLDEN1P + (IORD-0)*LRHO1
        CALL COPVEC(WORK(III),WORK(KRHO1),LRHO1)
        CALL ONE_EL_PROP
      END DO
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'PROPPT')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE LTRAPRP
*
* Calculate one-electron transition properties
*
* Jeppe Olsen, June 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER*8 KLRHOT, KLSCR, KLRHO1S, KLCMO, KLHONE, KLRHOTP 
!               for addressing of WORK
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "multd2h.inc"
#include "glbbas.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "cgas.inc"
#include "cstate.inc"
#include "clunit.inc"
*
      REAL*8 INPROD
      CHARACTER*1 XYZ(3)
      CHARACTER*8 LABEL
      CHARACTER*6 LABEL2
      DATA XYZ/'X','Y','Z'/
*. Common block for communicating with sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'TRAPRP')
*
      NCALC = NROOT*NEXCSTATE
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Welcome to TRAPRP'
        WRITE(6,*)
        WRITE(6,*) ' Number of transition densities ',NCALC
      END IF
*. A bit of local memory
      LHONE = NTOOB * NTOOB
      WRITE(6,*) ' Max size of one-electron operator',LHONE
      CALL MEMMAN(KLHONE,LHONE,'ADDL  ',2,'KLHONE')
      CALL MEMMAN(KLCMO,LHONE,'ADDL  ',2,'KLCMO ')
      CALL MEMMAN(KLRHO1S,LHONE,'ADDL  ',2,'LRHO1S')
      CALL MEMMAN(KLSCR,4*LHONE,'ADDL  ',2,'KLCMO ')
      CALL MEMMAN(KLRHOT,NCALC*LHONE,'ADDL  ',2,'RHOT  ')
*
* 1 : Obtain transition densities
*
*. Right states : current ci vectors
*. Left states  : states on LUEXC
*
*. Fill /CANDS/
      ICSM = IREFSM
      ISSM = IEXCSYM
*. The space is assumed to be the final space
      ICSPC = NCMBSPC
      ISSPC = NCMBSPC
*. and then the transition densities
* will be stored in KLRHOT as Loop over left , Loop over right
C     TRADEN(I12,RHO1,RHO2,NL,NR,LUL,LUR)
C     I12 = 1
C     ILR = 0
C     DO IL = 1, NEXCSTATE
C       DO IR = 1, NROOT
C         ILR = ILR + 1
C         KLRHOTP = KLRHOT + (ILR-1)*LHONE
          CALL TRADEN(I12,WORK(KLRHOT),RHO2,NEXCSTATE,
     &                NROOT,LUEXC,LUC)
C       END DO
C     END DO
*
      WRITE(6,*) ' First Transition density '
      CALL WRTMT_LU(WORK(KLRHOT),NTOOB,NTOOB,NTOOB,NTOOB)
*. Symmetry of transition densities - and therefore of operators
      CALL SYMCOM(3,1,IREFSM,IEXCSYM,ITRASYM)
      WRITE(6,*) ' Symmetry of transition densities',ITRASYM
*
      IRHO1SM = ITRASYM
      ILR = 0
      DO IL = 1, NEXCSTATE
        DO IR = 1, NROOT
          WRITE(6,*) ' info for IL, IR =', IL,IR
*. Extract symmetry blocks from complete one-electron density
          ILR = ILR + 1
          KLRHOTP = KLRHOT + (ILR-1)*LHONE
C              REORHO1(RHO1I,RHO1O,IRHO1SM)
          CALL REORHO1(WORK(KLRHOTP),WORK(KLRHO1S),IRHO1SM)
*. Number of elements in symmetry blocks of integrals and density
          LRHO1S = 0
          DO ISM = 1, NSMOB
            JSM = MULTD2H(ISM,IRHO1SM)
            LRHO1S = LRHO1S + NTOOBS(ISM)*NTOOBS(JSM)
          END DO
          WRITE(6,*) ' Number of elements in symmetry blocks ',
     &    LRHO1S
*
          DO IPROP =1, NPROP
            WRITE(6,'(A,A)') ' Property to be calculated',
     &                     PROPER(IPROP)
*. one- or two-dimensional tensor ?
            LABEL2 = PROPER(IPROP)
            WRITE(6,'(A,A)') ' LABEL2 =', LABEL2
C           GET_PROP_RANK(PROPER,IRANK)
            CALL GET_PROP_RANK(LABEL2,NRANK)
C           IF(LABEL2.EQ.'DIPOLE') THEN
C             NRANK = 1
C           ELSE IF(LABEL2.EQ.'THETA ' .OR.
C    &              LABEL2.EQ.'QUADRU' .OR.
C    &              LABEL2(1:3).EQ.'EFG' ) THEN
C             NRANK  = 2
C           ELSE
            IF(NRANK.EQ.-1) THEN
             WRITE(6,'(A,A)') ' Unknown operator ',PROPER(IPROP)
             NRANK = -1
            END IF
            WRITE(6,*) ' Rank of operator ', NRANK
*.
            IF(NRANK.EQ.1) THEN
              DO ICOMP = 1, 3
                IF(IXYZSYM(ICOMP).EQ.IRHO1SM) THEN
                  WRITE(6,*) ' right symmetry for component',ICOMP
*. Label of integrals on file- currently DALTON FORM !!
                  IF(PROPER(IPROP).EQ.'DIPOLE') THEN
                    LABEL =XYZ(ICOMP)//'DIPLEN '
                  END IF
                  WRITE(6,'(A,A)') ' Label ', LABEL
*. Obtain one-electron integrals
                   CALL GET_PROPINT(WORK(KLHONE),IRHO1SM,LABEL,
     &                  WORK(KLSCR),NTOOBS,NTOOBS,NSMOB,0)
*. and then : Expectation value
                   EXPVAL = INPROD(WORK(KLRHO1S),WORK(KLHONE),LRHO1S)
                   WRITE(6,'(A,A,E22.14)')
     &             ' Expectation value of ',LABEL,  EXPVAL
                 END IF
               END DO
            ELSE IF(NRANK.EQ.2) THEN
              DO ICOMP = 1,3
                DO JCOMP = 1,ICOMP
                  IF(MULTD2H(IXYZSYM(ICOMP),IXYZSYM(JCOMP))
     &            .EQ.IRHO1SM) THEN
                    WRITE(6,*) ' Right symmetry for components',
     &              ICOMP,JCOMP
*
C                   IF(LABEL2.EQ.'THETA ') THEN
*. Buckinghams traceless quadrupole moment
C                      LABEL=XYZ(JCOMP)//XYZ(ICOMP)//'THETA'
C                   ELSE IF(LABEL2.EQ.'QUADRU') THEN
C                     LABEL=XYZ(JCOMP)//XYZ(ICOMP)//'QUADRU'
C                   ELSE IF (LABEL2(1:3).EQ.'EFG' ) THEN
C
                    LABEL=XYZ(JCOMP)//XYZ(ICOMP)//LABEL2
C                   END IF
*. Obtain one-electron integrals
                   CALL GET_PROPINT(WORK(KLHONE),IRHO1SM,LABEL,
     &                  WORK(KLSCR),NTOOBS,NTOOBS,NSMOB,0)
*. and then : Expectation value
                   EXPVAL =
     &             INPROD(WORK(KLRHO1S),WORK(KLHONE),LRHO1S)
                   WRITE(6,'(A,A,E22.14)')
     &             ' Expectation value of ',LABEL,  EXPVAL
                  END IF
                END DO
              END DO
            END IF
*
          END DO
*
        END DO
      END DO
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'TRAPRP')
      RETURN
      END
