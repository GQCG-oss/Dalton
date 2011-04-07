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
*---------------------------------------------------------------------*
c/* Deck CC_HMAT */
*=====================================================================*
       SUBROUTINE CC_HMAT( LISTL,  ! inp: Zeta vector list
     &                     LISTA,  ! inp: A amplitude list
     &                     LISTB,  ! inp: B amplitude list
     &                     LISTC,  ! inp: C amplitude list
     &                     LISTD,  ! inp: D amplitude list
     &                     NHTRAN, ! inp: nb. of H matrix transf.
     &                     MXVEC,  ! inp: max. nb. of dot products 
     &                     IHTRAN, ! inp: index array of H mat. transf.
     &                     IHDOTS, ! inp: index array of dot products
     &                     HCON,   ! out: array for dot products
     &                     WORK,   ! inp: work space
     &                     LWORK,  ! scr: length of work space
     &                     IOPTRES)! inp: output option 
*---------------------------------------------------------------------*
*
*    Purpose: calculation of linear transformations of
*             3 amplitude vectors, T^A, T^B and T^C with the 
*             H matrix (third partial derivative of the lagrangian
*             with respect to T) for index combinations for Lambda,
*             T^A, T^B, T^C passed in IHTRAN.
* 
*    epsilon_mu = < Lambda | [[[[H,T^A],T^B],T^C],tau_mu] | HF >
*
*    return of the result vectors:
*
*           IOPTRES = 0 :  note used (but compare with CC_BMAT to see
*                          for which purpose it is reserved)
*
*           IOPTRES = 1 :  the vectors are kept and returned in WORK
*                          if possible, start addresses returned in
*                          IHTRAN(5,*). N.B.: if WORK is not large
*                          enough, CC_HMAT will stop!! 
*
*           IOPTRES = 3 :  each result vector is written to its own
*                          file by a call to CC_WRRSP, LISTD is used
*                          as list type and IHTRAN(5,*) as list index
*                          NOTE that IHTRAN(5,*) is in this case input!
*
*           IOPTRES = 4 :  each result vector is added to a vector on
*                          file by a call to CC_WARSP, LISTD is used
*                          as list type and IHTRAN(5,*) as list index
*                          NOTE that IHTRAN(5,*) is in this case input!
*
*           IOPTRES = 5 :  the result vectors are dotted on a array
*                          of vectors, the type of the arrays given
*                          by LISTD and the indeces from IHDOTS
*                          the result of the dot products is returned
*                          in the HCON array
*
*
*    symmetries/variables:
*              
*           EPSI1,  ISYRES : H matrix transformation
*           CTR2,   ISYCTR : response lagrangian multipliers 
*           T1AMPA, ISYMTA : A response amplitudes 
*           T1AMPB, ISYMTB : B response amplitudes
*           T1AMPC, ISYMTC : C response amplitudes
*           T1AMPD, ISYMTD : D response amplitudes
*
*     uses approximately V^2 O^2 + O^4 + O^3 work space
*            CC2:  (ia|jb), CTR2
*            CCSD: (ia|jb), CTR2
* 
*     Written by Christof Haettig, Februar 1998.
*     CC3 noddy version, April 2002, Christof Haettig
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccnoddy.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      CHARACTER MSGDBG*(17)
      PARAMETER (MSGDBG='[debug] CC_HMAT> ')

#if defined (SYS_CRAY)
      REAL ZERO, ONE, TWO, THREE
#else
      DOUBLE PRECISION ZERO, ONE, TWO, THREE
#endif
      PARAMETER (ZERO = 0.0d0, ONE = 1.0d0, TWO = 2.0d0, THREE = 3.0d0)

      INTEGER KDUM
      PARAMETER (KDUM   = +99 999 999)  ! dummy address
      INTEGER ISYOVOV
      PARAMETER (ISYOVOV = 1)


      CHARACTER*(*) LISTL, LISTA, LISTB, LISTC, LISTD
      INTEGER LWORK, IOPTRES, MXVEC, NHTRAN
      INTEGER IHTRAN(5,NHTRAN)
      INTEGER IHDOTS(MXVEC,NHTRAN)

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL HCON(MXVEC,NHTRAN)
      REAL DDOT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION HCON(MXVEC,NHTRAN)
      DOUBLE PRECISION DDOT
#endif

      CHARACTER MODEL*(10), MODELW*(10)
      LOGICAL LLSAME
      INTEGER IZETAV, ITAMPA, ITAMPB, ITAMPC, ITAMPD, ITRAN, IFILE
      INTEGER KT1AMPA, KT1AMPB, KT1AMPC, KT1AMPD, KEPSI1, KSTART
      INTEGER ISYRES, ISYCTR, ISYMTA, ISYMTB, ISYMTC, ISYABCD
      INTEGER KXIAJB, KCTR2, KEND1, LEND1, IOPT, IOPTW, KEND0, LEND0,
     &        KEPSI2, ILEN, IDBL

      INTEGER ILSTSYM

  
*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*

* check coupled cluster model:
      IF (CCS) THEN
         MODELW = 'CCS       '
         IOPTW  = 1
      ELSE IF (CC2) THEN
         MODELW = 'CC2       '
         IOPTW  = 1
      ELSE IF (CCSD) THEN
         MODELW = 'CCSD      '
         IOPTW  = 1
      ELSE IF (CC3) THEN
         MODELW = 'CC3       '
         IOPTW  = 3
      ELSE 
        CALL QUIT('Unknown CC model in H matrix transformation')
      END IF

      IF ( .not. (CCS .or. CC2 .or. CCSD .or. CC3) ) THEN
        WRITE(LUPRI,'(/a)') ' CC_HMAT called for a Coupled Cluster '
     &          //'method not implemented in CC_HMAT...'
        CALL QUIT('Unknown CC method in CC_HMAT.')
      END IF

* check list combination:
      IF (LISTA(1:1).NE.'R' 
     &    .OR. LISTB(1:1).NE.'R' 
     &    .OR. LISTC(1:1).NE.'R' 
     &    .OR. .NOT.(LISTD(1:1).EQ.'R' .OR. LISTD(1:1).EQ.'X')
     &    .OR. LISTL(1:1).NE.'L'                         ) THEN
        WRITE(LUPRI,'(2A,/A)') 
     &    ' In CC_HMAT LISTA, LISTB, LISTC, LISTD', 
     &    ' should refer to t-amplitude vectors and LISTL to a ',
     &    ' multiplier vector.'
        WRITE(LUPRI,'(2A)')
     &    ' LISTL: ',LISTL(1:3),
     &    ' LISTA: ',LISTA(1:3),
     &    ' LISTB: ',LISTB(1:3),
     &    ' LISTC: ',LISTC(1:3),
     &    ' LISTD: ',LISTD(1:3)
        CALL QUIT('Strange LIST combination in CC_HMAT.')
      END IF

* check IOPTRES and initialize output array HCON, if used:
      IF (IOPTRES.EQ.1  .OR.  IOPTRES.EQ.3  .OR.  IOPTRES.EQ.4 ) THEN
        CONTINUE
      ELSE IF (IOPTRES.EQ.5) THEN
        IF (MXVEC*NHTRAN .GT. 0) CALL DZERO(HCON,MXVEC*NHTRAN)
      ELSE
        CALL QUIT('Illegal value of IOPTRES in CC_HMAT.')
      END IF
  
* in debug mode print length of work space:
      IF (LOCDBG) THEN
        WRITE(LUPRI,'(/1x,a,i15)') 'work space in CC_HMAT:',LWORK
        Call FLSHFO(LUPRI)
      END IF

* flush print unit
      Call FLSHFO(LUPRI)

* initialize start address on the work space:
      KSTART = 1

*=====================================================================*
* start loop over all requeste H matrix transformations:
*=====================================================================*
      DO ITRAN = 1, NHTRAN

        IZETAV = IHTRAN(1,ITRAN)
        ITAMPA = IHTRAN(2,ITRAN)
        ITAMPB = IHTRAN(3,ITRAN)
        ITAMPC = IHTRAN(4,ITRAN)
        ITAMPD = IHTRAN(5,ITRAN)
        IFILE  = ITAMPD

* set & check symmetries:
        ISYCTR = ILSTSYM(LISTL,IZETAV)
        ISYMTA = ILSTSYM(LISTA,ITAMPA)
        ISYMTB = ILSTSYM(LISTB,ITAMPB)
        ISYMTC = ILSTSYM(LISTC,ITAMPC)

        ISYRES = MULD2H(MULD2H(ISYMTA,ISYMTB),MULD2H(ISYCTR,ISYMTC))

* allocated & initialize single excitation part of result vector EPSI1:
        KEPSI1 = KSTART
        KEND0  = KEPSI1 + NT1AM(ISYRES)
        LEND0  = LWORK - KEND0

        IF (LEND0 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_HMAT. (a)')
        END IF

        Call DZERO (WORK(KEPSI1), NT1AM(ISYRES))
      
*---------------------------------------------------------------------*
* calculate H matrix transformation (for CCS the result is zero!):
*---------------------------------------------------------------------*
       IF (.NOT. CCS) THEN

*---------------------------------------------------------------------*
* initialize pointer for work space & read single parts of the vectors
*---------------------------------------------------------------------*
        KT1AMPA = KEND0
        KT1AMPB = KT1AMPA + NT1AM(ISYMTA)
        KT1AMPC = KT1AMPB + NT1AM(ISYMTB)
        KCTR2   = KT1AMPC + NT1AM(ISYMTC)
        KXIAJB  = KCTR2   + NT2AM(ISYCTR)
        KEND1   = KXIAJB  + NT2AM(ISYOVOV)
        LEND1   = LWORK - KEND1

        IF (LEND1 .LT. 0) THEN
          CALL QUIT('Insufficient work space in CC_HMAT. (0)')
        END IF

* * read packed (ia|jb) integrals:
        Call CCG_RDIAJB(WORK(KXIAJB),NT2AM(ISYOVOV))

* read packed multipliers:
        IOPT = 2
        Call CC_RDRSP(LISTL,IZETAV,ISYCTR,IOPT,MODEL,
     &                WORK(KDUM),WORK(KCTR2))

* read A response amplitudes:
        IOPT = 1
        Call CC_RDRSP(LISTA,ITAMPA,ISYMTA,IOPT,MODEL,
     &                WORK(KT1AMPA),WORK(KDUM))

* read B response amplitudes:
        IOPT = 1
        Call CC_RDRSP(LISTB,ITAMPB,ISYMTB,IOPT,MODEL,
     &                WORK(KT1AMPB),WORK(KDUM))

* read C response amplitudes:
        IOPT = 1
        Call CC_RDRSP(LISTC,ITAMPC,ISYMTC,IOPT,MODEL,
     &                WORK(KT1AMPC),WORK(KDUM))

        IF (LOCDBG) THEN
          Call AROUND('debug_CC_HMAT> response T1 amplitudes B:')
          WRITE (LUPRI,*) 'List, Index, Symmetry:',LISTA, ITAMPA, ISYMTA
          Call CC_PRP(WORK(KT1AMPA),WORK(KDUM),ISYMTA,1,0)

          Call AROUND('debug_CC_HMAT> response T1 amplitudes C:')
          WRITE (LUPRI,*) 'List, Index, Symmetry:',LISTB, ITAMPB, ISYMTB
          Call CC_PRP(WORK(KT1AMPB),WORK(KDUM),ISYMTB,1,0)

          Call AROUND('debug_CC_HMAT> response T1 amplitudes D:')
          WRITE (LUPRI,*) 'List, Index, Symmetry:',LISTC, ITAMPC, ISYMTC
          Call CC_PRP(WORK(KT1AMPC),WORK(KDUM),ISYMTC,1,0)

          Call FLSHFO(LUPRI)
        END IF

*---------------------------------------------------------------------*
* calculate C and D terms for the cyclic permutations of B, C, D:
* (permutation of the first two vectors is treated inside CC_HCD)
*---------------------------------------------------------------------*
        CALL CC_HCD( WORK(KXIAJB),  ISYOVOV,!inp: (ia|jb) integrals
     &               WORK(KCTR2),   ISYCTR, !inp: double of Zeta vector
     &               WORK(KT1AMPA), ISYMTA, !inp: T1 amplitudes for B
     &               WORK(KT1AMPB), ISYMTB, !inp: T1 amplitudes for C
     &               WORK(KT1AMPC), ISYMTC, !inp: T1 amplitudes for D
     &               WORK(KEPSI1),  ISYRES, !out: Epsilon result vector
     &               WORK(KEND1),   LEND1  )!wrk: work space    

        CALL CC_HCD( WORK(KXIAJB),  ISYOVOV,!inp: (ia|jb) integrals
     &               WORK(KCTR2),   ISYCTR, !inp: double of Zeta vector
     &               WORK(KT1AMPB), ISYMTB, !inp: T1 amplitudes for C
     &               WORK(KT1AMPC), ISYMTC, !inp: T1 amplitudes for D
     &               WORK(KT1AMPA), ISYMTA, !inp: T1 amplitudes for B
     &               WORK(KEPSI1),  ISYRES, !out: Epsilon result vector
     &               WORK(KEND1),   LEND1  )!wrk: work space    
 
        CALL CC_HCD( WORK(KXIAJB),  ISYOVOV,!inp: (ia|jb) integrals
     &               WORK(KCTR2),   ISYCTR, !inp: double of Zeta vector
     &               WORK(KT1AMPC), ISYMTC, !inp: T1 amplitudes for D
     &               WORK(KT1AMPA), ISYMTA, !inp: T1 amplitudes for B
     &               WORK(KT1AMPB), ISYMTB, !inp: T1 amplitudes for C
     &               WORK(KEPSI1),  ISYRES, !out: Epsilon result vector
     &               WORK(KEND1),   LEND1  )!wrk: work space    

       END IF ! (.NOT. CCS) 

*---------------------------------------------------------------------*
* add CC3 contributions:
*---------------------------------------------------------------------*
        IF (CCSDT) THEN
          ! allocate and initialize doubles result vector:
          KEPSI1 = KSTART
          KEPSI2 = KEPSI1 + NT1AM(ISYRES)
          KEND0  = KEPSI2 + NT2AM(ISYRES)
          LEND0  = LWORK - KEND0

          IF (LEND0 .LT. 0) THEN
            CALL QUIT('Insufficient work space in CC_HMAT. (ccsdt)')
          END IF

          Call DZERO (WORK(KEPSI2), NT2AM(ISYRES))

          IF (NODDY_HMAT) THEN

             CALL CCSDT_HMAT_NODDY(LISTL,IZETAV,LISTA,ITAMPA,
     &                          LISTB,ITAMPB,LISTC,ITAMPC,
     &                          WORK(KEPSI1),WORK(KEPSI2),
     &                          WORK(KEND0),LEND0)

          ELSE 

             CALL CC3_HMAT(LISTL,IZETAV,LISTA,ITAMPA,
     &                     LISTB,ITAMPB,LISTC,ITAMPC,
     &                     WORK(KEPSI1),WORK(KEPSI2),ISYRES,
     &                     WORK(KEND0),LEND0)

          END IF

        ELSE
          KEPSI2 = KSTART
C         ! to not get undefined address in calls below /hjaaj Aug 07
        END IF
*=====================================================================*
* store result vectors:
*=====================================================================*
        IF ( IOPTRES .EQ. 1 ) THEN

*         output returned in work space: just increase start address

          IHTRAN(5,ITRAN) = KSTART
          KSTART = KSTART + NT1AM(ISYRES)

        ELSE IF (IOPTRES .EQ. 3) THEN

*         write to file using CC_WRRSP, pass LISTD as list and
*         IFILE as index

          CALL CC_WRRSP(LISTD,IFILE,ISYRES,IOPTW,MODELW,WORK(KDUM),
     &                  WORK(KEPSI1),WORK(KEPSI2),WORK(KEND0),LEND0)

        ELSE IF (IOPTRES .EQ. 4) THEN

*         add to vector on file using CC_WRRSP, pass LISTD as list
*         and IFILE as index

          CALL CC_WARSP(LISTD,IFILE,ISYRES,IOPTW,MODELW,WORK(KDUM),
     &                  WORK(KEPSI1),WORK(KEPSI2),WORK(KEND0),LEND0)

        ELSE IF (IOPTRES .EQ. 5) THEN

*         calculate list of dot products

          IF (IOPTW.GE.2) CALL CCLR_DIASCL(WORK(KEPSI2),TWO,ISYRES)
          CALL CCDOTRSP(IHDOTS,HCON,IOPTW,LISTD,ITRAN,NHTRAN,MXVEC,
     &                  WORK(KEPSI1),WORK(KEPSI2),ISYRES,
     &                  WORK(KEND0),LEND0)

        ELSE
          CALL QUIT('Illegal value for IOPTRES in CC_HMAT.')
        END IF

*---------------------------------------------------------------------*
* print debug information:
*---------------------------------------------------------------------*
        IF (LOCDBG) THEN
          WRITE (LUPRI,*) MSGDBG, ' Lambda vector     : ',LISTL,IZETAV
          WRITE (LUPRI,*) MSGDBG, ' A amplitude vector: ',LISTA,ITAMPA
          WRITE (LUPRI,*) MSGDBG, ' B amplitude vector: ',LISTB,ITAMPB
          WRITE (LUPRI,*) MSGDBG, ' C amplitude vector: ',LISTC,ITAMPC
          WRITE (LUPRI,*) MSGDBG, ' epsilon vector:'
          ILEN = NT1AM(ISYRES)
          IDBL = 0
          IF (IOPTW.GE.2) THEN
            ILEN = ILEN + NT2AM(ISYRES)
            IDBL = 1
          END IF
          Call CC_PRP(WORK(KEPSI1),WORK(KEPSI2),ISYRES,1,IDBL)
          WRITE (LUPRI,*) MSGDBG, ' Norm^2 of EPSI1:',
     &       DSQRT(DDOT(ILEN,WORK(KEPSI1),1,WORK(KEPSI1),1))
          Call FLSHFO(LUPRI)
        END IF

*---------------------------------------------------------------------*
      END DO ! ITRAN 

      RETURN
      END
*=====================================================================*
*                END OF SUBROUTINE CC_HMAT
*=====================================================================*

*---------------------------------------------------------------------*
c/* Deck CC_HCD */
*=====================================================================*
       SUBROUTINE CC_HCD( XIAJB,  ISYOVOV, ! inp: (ia|jb) integrals
     &                    CTR2,   ISYCTR,  ! inp: double of Zeta vector
     &                    T1AMPB, ISYMTB,  ! inp: T1 amplitudes for B
     &                    T1AMPC, ISYMTC,  ! inp: T1 amplitudes for C
     &                    T1AMPD, ISYMTD,  ! inp: T1 amplitudes for D
     &                    EPSI1,  ISYRES,  ! out: Epsilon result vector
     &                    WORK,   LWORK   )! wrk: work space    
*---------------------------------------------------------------------*
*
*    Purpose: calculate C and D term contributions for H matrix
*             for a particular permutation of the amplitude vectors
*
*    Written by Christof Haettig, Februar 1998.
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>

* local parameters:
      CHARACTER MSGDBG*(16)
      PARAMETER (MSGDBG='[debug] CC_HCD> ')
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYMTB, ISYMTC, ISYMTD, ISYCTR, ISYRES, ISYOVOV, LWORK
  
#if defined (SYS_CRAY)
      REAL T1AMPB(*)      ! dimension (NT1AM(ISYMTB))
      REAL T1AMPC(*)      ! dimension (NT1AM(ISYMTC))
      REAL T1AMPD(*)      ! dimension (NT1AM(ISYMTD))
      REAL EPSI1(*)       ! dimension (NT1AM(ISYRES))
      REAL XIAJB(*)       ! dimension (NT2AM(ISYOVOV))
      REAL CTR2(*)        ! dimension (NT2AM(ISYCTR))
      REAL WORK(LWORK)
      REAL DDOT
#else
      DOUBLE PRECISION T1AMPB(*)      ! dimension (NT1AM(ISYMTB))
      DOUBLE PRECISION T1AMPC(*)      ! dimension (NT1AM(ISYMTC))
      DOUBLE PRECISION T1AMPD(*)      ! dimension (NT1AM(ISYMTD))
      DOUBLE PRECISION EPSI1(*)       ! dimension (NT1AM(ISYRES))
      DOUBLE PRECISION XIAJB(*)       ! dimension (NT2AM(ISYOVOV))
      DOUBLE PRECISION CTR2(*)        ! dimension (NT2AM(ISYCTR))
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION DDOT
#endif

      INTEGER ISYMTBC, ISYMA, ISYMLJK, ISYC4O, ISYX4O, ISYMI
      INTEGER KC4O, KX4O, KEND1, LEND1, NAI, KOFF, KXLJKA, KCLJKA, IOPT


*---------------------------------------------------------------------*
* begin:
*---------------------------------------------------------------------*
      ISYMTBC = MULD2H(ISYMTB,ISYMTC)

*---------------------------------------------------------------------*
* D term:
*---------------------------------------------------------------------*

* calculate double one-index transformed Zeta vector:
      ISYC4O  = MULD2H(ISYCTR,ISYMTBC)

      KC4O  = 1
      KEND1 = KC4O  + N3ORHF(ISYC4O)
      LEND1 = LWORK - KEND1

      IF ( LEND1 .LT. 0 ) THEN
        CALL QUIT('Insufficient work space in CC_HCD. (1)')
      END IF

      IOPT = 3
      Call CCG_4O(WORK(KC4O),ISYC4O, CTR2,ISYCTR,
     &            T1AMPB,ISYMTB,T1AMPC,ISYMTC,
     &            WORK(KEND1),LEND1,IOPT)


* start loop over virtual index a:
      DO ISYMA = 1, NSYM

        ISYMI   = MULD2H(ISYRES,ISYMA)
        ISYMLJK = MULD2H(ISYOVOV,MULD2H(ISYMTD,ISYMA))

        KXLJKA  = KEND1
        If ( LEND1 .LT. NMAIJK(ISYMLJK) ) THEN
          CALL QUIT('Insufficient work space in CC_HCD. (2)')
        END IF

      DO A = 1, NVIR(ISYMA)

* calculate batch of (l j^D | k a) integrals with fixed a:
        IOPT = 1
        Call CCG_TRANS2(WORK(KXLJKA),ISYMLJK,XIAJB,ISYOVOV,
     &                  T1AMPD,ISYMTD,A,ISYMA,IOPT)

* contract (l j^D|k a) with [ Zeta(l^C j|k^B i) + Zeta(l^B j|k^C i) ]:
        DO I = 1, NRHF(ISYMI)
          NAI  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
          KOFF = I3ORHF(ISYMLJK,ISYMI) + NMAIJK(ISYMLJK)*(I-1)
          EPSI1(NAI) = EPSI1(NAI) + 
     &       DDOT(NMAIJK(ISYMLJK),WORK(KC4O+KOFF),1,WORK(KXLJKA),1)
        END DO

      END DO ! A
      END DO ! ISYMA
      

*---------------------------------------------------------------------*
* D term:
*---------------------------------------------------------------------*
      ISYX4O  = MULD2H(ISYOVOV,ISYMTBC)

      KX4O  = 1
      KEND1 = KX4O  + N3ORHF(ISYX4O)
      LEND1 = LWORK - KEND1

      IF ( LEND1 .LT. 0 ) THEN
        CALL QUIT('Insufficient work space in CC_HCD. (3)')
      END IF

* calculate double one-index transformed integrals:
      IOPT = 3
      Call CCG_4O(WORK(KX4O),ISYX4O, XIAJB,ISYOVOV,
     &            T1AMPB,ISYMTB,T1AMPC,ISYMTC,
     &            WORK(KEND1),LEND1,IOPT)

C     WRITE (LUPRI,*) 'KX4O:'
C     WRITE (LUPRI,'(3X,I5,3X,D25.14)'), 
C    &  (J+1,WORK(KX4O+J),J=0,N3ORHF(ISYX4O)-1)

* start loop over virtual index a:
      DO ISYMA = 1, NSYM

        ISYMI   = MULD2H(ISYRES,ISYMA)
        ISYMLJK = MULD2H(ISYCTR,MULD2H(ISYMTD,ISYMA))

        KCLJKA  = KEND1
        If ( LEND1 .LT. NMAIJK(ISYMLJK) ) THEN
          CALL QUIT('Insufficient work space in CC_HCD. (4)')
        END IF

      DO A = 1, NVIR(ISYMA)

        CALL DCOPY(NMAIJK(ISYMLJK),999.99d0,0,WORK(KCLJKA),1)
* calculate batch of one-index transf. Zeta(l j^D | k a) with fixed a:
        IOPT = 1
        Call CCG_TRANS2(WORK(KCLJKA),ISYMLJK,CTR2,ISYCTR,
     &                  T1AMPD,ISYMTD,A,ISYMA,IOPT)

* contract Zeta(l j^D|k a) with [ (l^C j|k^B i) + (l^B j|k^C i) ]:
        DO I = 1, NRHF(ISYMI)
          NAI  = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
          KOFF = I3ORHF(ISYMLJK,ISYMI) + NMAIJK(ISYMLJK)*(I-1)
          EPSI1(NAI) = EPSI1(NAI) + 
     &       DDOT(NMAIJK(ISYMLJK),WORK(KX4O+KOFF),1,WORK(KCLJKA),1)

C         WRITE (LUPRI,*) 'KX4O               LJKA:'
C         WRITE (LUPRI,'(3X,I5,3X,2D25.14)'), 
C    &    (J+1,WORK(KX4O+KOFF+J),WORK(KCLJKA+J),J=0,NMAIJK(ISYMLJK)-1)
        END DO

      END DO ! A
      END DO ! ISYMA

      RETURN
      END
*=====================================================================*
*                   END OF SUBROUTINE CC_HCD                          *
*=====================================================================*

*=====================================================================*
      SUBROUTINE CC_HTST(WORK,LWORK)
*=====================================================================*
C
C perform finite difference test for H matrix:
C
C----------------------------------------------------------------------
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccqrinf.h>

      INTEGER MXHTRAN, MXVEC
      PARAMETER ( MXHTRAN = 5, MXVEC = 1)

      CHARACTER*(3) LISTB, LISTC, LISTD, LISTL, LISTA
      CHARACTER*(10) MODEL
      INTEGER IDLSTB, IDLSTC, IDLSTD, IDLSTL
      INTEGER KT1AMPB, KT2AMPB, KEND, LEND, LWORK
      INTEGER ISYMB, ISYMC, ISYMD, ISYML, ISYRES
      INTEGER KRESLT1, KRESLT2, IOPT, KEPSI1, KEPSI1A, IOPTRES
      INTEGER IHTRAN(5,MXHTRAN), NHTRAN
      INTEGER IHDOTS(MXVEC,MXHTRAN)

#if defined (SYS_CRAY)
      REAL WORK(LWORK)
      REAL HCON(MXVEC,MXHTRAN)
      REAL DDOT
#else
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION HCON(MXVEC,MXHTRAN)
      DOUBLE PRECISION DDOT
#endif

      INTEGER ILSTSYM
      INTEGER IR1TAMP

*---------------------------------------------------------------------*
* call G matrix transformation
*---------------------------------------------------------------------*
      LISTL = 'L0'
      LISTA = 'R1'
      LISTB = 'R1'
      LISTC = 'R1'
      LISTD = 'R1'

      IDLSTL = 0
      IDLSTB = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDLSTC = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,1)
      IDLSTD = IR1TAMP('XDIPLEN ',.FALSE.,0.0D0,1)

      IHTRAN(1,1) = IDLSTL
      IHTRAN(2,1) = IDLSTB
      IHTRAN(3,1) = IDLSTC
      IHTRAN(4,1) = IDLSTD
      NHTRAN = 1

      ISYML  = ILSTSYM(LISTL,IDLSTL)
      ISYMB  = ILSTSYM(LISTB,IDLSTB)
      ISYMC  = ILSTSYM(LISTB,IDLSTC)
      ISYMD  = ILSTSYM(LISTD,IDLSTD)
      ISYRES = MULD2H(MULD2H(ISYML,ISYMD),MULD2H(ISYMB,ISYMC))

      KEPSI1  = 1
      KEPSI1A = KEPSI1  + NT1AM(ISYRES)
      KEND    = KEPSI1A + NT1AM(ISYRES)
      LEND    = LWORK - KEND

C
C old H matrix routine:
C
C     WRITE (LUPRI,*) 'CC_HTST: CALL NOW CCCR_K'
C     CALL CCCR_K(LISTL,IDLSTL,LISTB,IDLSTB,LISTC,IDLSTC,LISTD,IDLSTD,
C    &            WORK(KEPSI1),ISYRES,WORK(KEND),LEND)
C     WRITE (LUPRI,*) 'CC_HTST: RETURNED FROM CCCR_K'

      WRITE (LUPRI,*) 'CC_HTST: CALL NOW CC_HMAT'
      IOPTRES = 1
      CALL CC_HMAT(LISTL,LISTB,LISTC,LISTD,LISTA,NHTRAN,MXVEC,
     &             IHTRAN,IHDOTS,HCON,
     &             WORK(KEPSI1A),LWORK-KEPSI1A,IOPTRES)
      WRITE (LUPRI,*) 'CC_HTST: RETURNED FROM CC_HMAT'

      KT1AMPB = KEND
      KT2AMPB = KT1AMPB + NT1AM(ISYMB)
      KRESLT1 = KT2AMPB + NT2AM(ISYMB)
      KRESLT2 = KRESLT1 + NT1AM(ISYRES)
      KEND    = KRESLT2 + NT2AM(ISYRES)
      LEND    = LWORK - KEND
  
      IF (LEND .LT. 0) THEN
        CALL QUIT('Insufficient work space in CC_HTST.')
      END IF

      IOPT = 3
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KT1AMPB),WORK(KT2AMPB))

* zero doubles of B and/or C vector:
C     CALL DZERO(WORK(KT2AMPB),NT2AM(ISYMB))
       
      CALL CC_FDH(NT1AM(ISYMB),NT2AM(ISYMB),
     >              LISTC,IDLSTC,LISTD,IDLSTD,
     >              WORK(KT1AMPB), WORK(KRESLT1),
     >              WORK(KEND), LEND)

      IF (CCS) CALL DZERO(WORK(KRESLT2),NT2AM(ISYRES))

      IF (.TRUE.) THEN
        WRITE (LUPRI,*) 'finite difference Epsilon vector:'
        Call CC_PRP(WORK(KRESLT1),WORK(KRESLT2),ISYRES,1,1)
        WRITE (LUPRI,*) 'analytical Epsilon vector (CCCR_K):'
        Call CC_PRP(WORK(KEPSI1),WORK,ISYRES,1,0)
        WRITE (LUPRI,*) 'analytical Epsilon vector (CC_HMAT):'
        Call CC_PRP(WORK(KEPSI1A),WORK,ISYRES,1,0)
      END IF

      Call DAXPY(NT1AM(ISYRES),-1.0d0,WORK(KRESLT1),1,WORK(KEPSI1),1)

C
C compare with result of old H matrix routine:
C
C     WRITE (LUPRI,*) 'Norm of difference between analytical Epsilon '
C    >           // 'vector (CCCR_K) and the numerical result:'
C     WRITE (LUPRI,*) 'singles excitation part:',
C    > DSQRT(DDOT(NT1AM(ISYRES),WORK(KEPSI1),1,WORK(KEPSI1),1))
C     WRITE (LUPRI,*) 'double excitation part: ',
C    > DSQRT(DDOT(NT2AM(ISYRES),WORK(KRESLT2),1,WORK(KRESLT2),1))


C
C compare with result of new H matrix routine:
C
      Call DAXPY(NT1AM(ISYRES),-1.0d0,WORK(KRESLT1),1,WORK(KEPSI1A),1)

      WRITE (LUPRI,*) 'Norm of difference between analytical Epsilon '
     >           // 'vector (CC_HMAT) and the numerical result:'
      WRITE (LUPRI,*) 'singles excitation part:',
     > DSQRT(DDOT(NT1AM(ISYRES),WORK(KEPSI1A),1,WORK(KEPSI1A),1))
      WRITE (LUPRI,*) 'double excitation part: ',
     > DSQRT(DDOT(NT2AM(ISYRES),WORK(KRESLT2),1,WORK(KRESLT2),1))
      CALL FLSHFO(LUPRI)

      RETURN
      END 
C----------------------------------------------------------------------
      SUBROUTINE CC_FDH(NC1VEC,NC2VEC,LISTB,ITAMPB,LISTC,ITAMPC,
     &                  TZAM,RESULT,WORK,LWORK)
C
C----------------------------------------------------------------------
C     Test routine for calculating the CC K*t*t*t vector by
C     finite difference on the H-matrix transformation.
C     C.Haettig and Ove Christiansen oktober 1996, februar 1997
C
C----------------------------------------------------------------------
C
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxorb.h>
#include <iratdef.h>
#include <ccorb.h>
#include <aovec.h>
#include <ccsdinp.h>
#include <cclr.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <leinf.h>
C
      DIMENSION WORK(LWORK),ITADR(2)
      PARAMETER (XHALF = 0.5D00,XMTWO = -2.0D00, DELTA = 1.0D-07)
      PARAMETER (ONE = 1.0D00, ZERO =  0.0D00 )
      CHARACTER*10 MODEL
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      IF (CCR12) CALL QUIT('Finite-difference on H-matrix for CCR12 '//
     &                     'not adapted')
C
      IF (IPRINT.GT.5) THEN
         CALL AROUND( 'IN CC_FDH  : MAKING FINITE DIFF. CC H-Matrix')
      ENDIF
C
C----------------------------
C     Work space allocations.
C----------------------------
C
      ISYMTR     = 1
      ISYMOP     = 1
C
      NTAMP      = NT1AM(ISYMTR) + NT2AM(ISYMTR)
      NTAMP2     = NTAMP*(NC1VEC + NC2VEC )
      KF         = 1
      KRHO1      = KF       + NTAMP2
      KRHO2      = KRHO1    + NT1AMX
      KC1AM      = KRHO2    + MAX(NT2AMX,NT2AM(ISYMTR))
      KC2AM      = KC1AM    + NT1AM(ISYMTR)
      KEND1      = KC2AM 
     *           + MAX(NT2AMX,NT2AM(ISYMTR),NT2AO(ISYMTR),
     *                 2*NT2ORT(ISYMTR))
      LWRK1      = LWORK    - KEND1
C
      KRHO1D     = KEND1
      KRHO2D     = KRHO1D   + NT1AMX
      KEND2      = KRHO2D     
     *           + MAX(NT2AMX,NT2AM(ISYMTR),NT2AO(ISYMTR),
     *                 2*NT2ORT(ISYMTR))
      LWRK2      = LWORK      - KEND1
C
      IF (IPRINT .GT. 100 ) THEN
         WRITE(LUPRI,*) ' IN CC_FDH: KF      =  ',KF     
         WRITE(LUPRI,*) ' IN CC_FDH: KRHO1   =  ',KRHO1
         WRITE(LUPRI,*) ' IN CC_FDH: KRHO2   =  ',KRHO2
         WRITE(LUPRI,*) ' IN CC_FDH: KC1AM   =  ',KC1AM
         WRITE(LUPRI,*) ' IN CC_FDH: KC2AM   =  ',KC2AM
         WRITE(LUPRI,*) ' IN CC_FDH: KRHO1D  =  ',KRHO1D
         WRITE(LUPRI,*) ' IN CC_FDH: KRHO2D  =  ',KRHO2D
         WRITE(LUPRI,*) ' IN CC_FDH: KEND2   =  ',KEND2
         WRITE(LUPRI,*) ' IN CC_FDH: LWRK2   =  ',LWRK2
      ENDIF
      IF (LWRK2.LT.0 ) THEN
         WRITE(LUPRI,*) 'Too little work space in CC_FDH '
         WRITE(LUPRI,*) 'AVAILABLE: LWORK   =  ',LWORK
         WRITE(LUPRI,*) 'NEEDED (AT LEAST)  =  ',KEND2
         CALL QUIT('TOO LITTLE WORKSPACE IN CC_FDH ')
      ENDIF
      KF2   = KF      + NC1VEC*NTAMP
C
C---------------------
C     Initializations.
C---------------------
C
      CALL DZERO(WORK(KC1AM),NT1AMX)
      CALL DZERO(WORK(KC2AM),NT2AMX)
      CALL DZERO(WORK(KF),NTAMP2)
      IF (ABS(DELTA) .GT. 1.0D-15 ) THEN 
         DELTAI = 1.0D00/DELTA
      ELSE
         DELTAI = 1
      ENDIF
      X11 = 0.0D00
      X12 = 0.0D00
      X21 = 0.0D00
      X22 = 0.0D00
      XNJ = 0.0D00
C
C------------------------------------------------
C     Read the CC reference amplitudes From disk.
C------------------------------------------------
C
      IOPT = 3
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KC1AM),WORK(KC2AM))
C
C----------------------------------------------
C     Save the CC reference amplitudes on disk.
C----------------------------------------------
C
      LUTAM = -1
      CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND(LUTAM)
      WRITE(LUTAM) (WORK(KC1AM + I -1 ), I = 1, NT1AMX)
      WRITE(LUTAM) (WORK(KC2AM + I -1 ), I = 1, NT2AMX)
      CALL GPCLOSE(LUTAM,'KEEP')
C
      IF (IPRINT.GT.125) THEN
         RHO1N = DDOT(NT1AMX,WORK(KC1AM),1,WORK(KC1AM),1)
         RHO2N = DDOT(NT2AMX,WORK(KC2AM),1,WORK(KC2AM),1)
         WRITE(LUPRI,*) 'Norm of T1AM: ',RHO1N
         WRITE(LUPRI,*) 'Norm of T2AM: ',RHO2N
         CALL CC_PRP(WORK(KC1AM),WORK(KC2AM),1,1,1)
      ENDIF
      RSPIM = .TRUE.
C
C------------------------------------
C     Calculate reference G*T*T vector.
C------------------------------------
C
      CALL QUIT('CC_FDH has to be fixed because of new G matrix.')
C     CALL CCQR_G('L0',0,LISTB,ITAMPB,LISTC,ITAMPC,WORK(KRHO1D),
C    &            ISYMTR,WORK(KRHO2D),LWORK-KRHO2D)

C
C-------------------------
C     Zero out components.
C-------------------------
C
      IF (LCOR .OR. LSEC) THEN
C
         CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
C
      ENDIF
C
      IF (IPRINT.GT.2) THEN
         RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
         RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
         WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'ref'
         WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'ref'
      ENDIF
      IF (IPRINT.GT.125) THEN
         CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
      ENDIF
C
      CALL DCOPY(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1),1)
      CALL DCOPY(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2),1)
C
C=============================================
C     Calculate H-matrix by finite difference.
C=============================================
C
      DO 100 I = 1, NC1VEC
           WRITE (LUPRI,*) 'singles index:',I
C
C----------------------------------------
C        Add finite displadement to t and 
C        calculate new intermediates.
C----------------------------------------
C
         LUTAM = -1
         CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
         READ(LUTAM) (WORK(KC1AM + J -1 ) , J = 1, NT1AMX)
         READ(LUTAM) (WORK(KC2AM + J -1 ) , J = 1, NT2AMX)
         CALL GPCLOSE(LUTAM,'KEEP')
C
         TI   = SECOND()
         WORK(KC1AM +I -1) = WORK(KC1AM +I -1 ) + DELTA
         IF (LCOR .OR. LSEC) THEN
            CALL CC_CORE(WORK(KC1AM),WORK(KC2AM),ISYMTR)
         ENDIF
C
         IOPT = 3
         CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *                 WORK(KC2AM),WORK(KEND2),LWRK2)
C
         RSPIM = .TRUE.
         CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *               WORK(KC2AM),WORK(KEND2),LWRK2,'XXX') 
C
C-----------------------------------------
C        Get the CC response vector again.
C-----------------------------------------
C
C        CALL DCOPY(NTAMP,TXAM,1,WORK(KC1AM),1)
C
C---------------------------------------
C        For Test zero part of L vector.
C---------------------------------------
C     
C        IF ( L1TST ) THEN
C           CALL DZERO(WORK(KC2AM),NT2AMX)
C        ENDIF
C        IF ( L2TST ) THEN
C           CALL DZERO(WORK(KC1AM),NT1AMX)
C        ENDIF
C
C------------------
C        Transform.
C------------------
C
         CALL QUIT('CC_FDH has to be fixed because of new G matrix.')
C        CALL CCQR_G('L0',0,LISTB,ITAMPB,LISTC,ITAMPC,WORK(KRHO1D),
C    &               ISYMTR,WORK(KRHO2D),LWORK-KRHO2D)
C
         IF (LCOR .OR. LSEC) THEN
            CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
         ENDIF
C
         IF (IPRINT.GT.2) THEN
            RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
            RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
            WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'ai=',I
            WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'ai=',I
         ENDIF
         IF (IPRINT.GT.125) THEN
            CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
         ENDIF
         CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
         CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
         CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
         CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)
         CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *              WORK(KF+NTAMP*(I-1)),1)
         CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *             WORK(KF+NTAMP*(I-1)+NT1AMX),1)
         X11 = X11 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
         X21 = X21 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
C
         TI   = SECOND() - TI
         IF (IPRINT.GT.5 ) THEN
            WRITE(LUPRI,*) '  '
            WRITE(LUPRI,*) 'FDH ROW NR. ',I,' DONE IN ',TI,' SEC.'
         ENDIF
C
 100  CONTINUE
C
C----------------------------------------------------------------
C     Loop over T2 amplitudes. Take care of diagonal t2 elements
C     is in a different convention in the energy code.
C     Factor 1/2 from right , and factor 2 from left.
C----------------------------------------------------------------
C
      DO 200 NAI = 1, NT1AMX
        DO 300 NBJ = 1, NAI
         I = INDEX(NAI,NBJ)
C
         IF (I.LE.NC2VEC) THEN

           WRITE (LUPRI,*) 'doubles index:',I
C
C--------------------------------------------
C          Add finite displacement to t and
C          calculate new intermediates.
C-------------------------------------------
C
           LUTAM = -1
           CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',
     *                 IDUMMY,.FALSE.)
           READ(LUTAM) (WORK(KC1AM + J -1 ) , J = 1, NT1AMX)
           READ(LUTAM) (WORK(KC2AM + J -1 ) , J = 1, NT2AMX)
           CALL GPCLOSE(LUTAM,'KEEP')
C
           TI   = SECOND()
           DELT = DELTA
           IF (NAI.EQ.NBJ) DELT = 2*DELTA
           WORK(KC2AM + I -1) = WORK(KC2AM+I -1) + DELT
           IF (LCOR .OR. LSEC) THEN
             CALL CC_CORE(WORK(KC1AM),WORK(KC2AM),ISYMTR)
           ENDIF
C
           IOPT = 3
           CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *                   WORK(KC2AM),WORK(KEND2),LWRK2)
C
           RSPIM = .TRUE.
           CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     *                 WORK(KC2AM),WORK(KEND2),LWRK2,'XXX') 
C
C-------------------------------------------
C          Get the CC response vector again.
C-------------------------------------------
C
C          CALL DCOPY(NTAMP,TXAM,1,WORK(KC1AM),1)
C
C-----------------------------------------
C          For Test zero part of L vector.
C-----------------------------------------
C     
C          IF ( L1TST ) THEN
C             CALL DZERO(WORK(KC2AM),NT2AMX)
C          ENDIF
C          IF ( L2TST ) THEN
C             CALL DZERO(WORK(KC1AM),NT1AMX)
C          ENDIF
C
C          RHO1N = DDOT(NT1AMX,WORK(KC1AM),1,WORK(KC1AM),1)
C          RHO2N = DDOT(NT2AMX,WORK(KC2AM),1,WORK(KC2AM),1)
C          IF ( DEBUG ) THEN
C             WRITE(LUPRI,*) 'Norm of L1AM-inp: ',RHO1N
C             WRITE(LUPRI,*) 'Norm of L2AM-inp: ',RHO2N
C          ENDIF
C
C--------------------
C          Transform.
C--------------------
C
           CALL QUIT('CC_FDH has to be fixed because of new G matrix.')
C          CALL CCQR_G('L0',0,LISTB,ITAMPB,LISTC,ITAMPC,WORK(KRHO1D),
C    &                 ISYMTR,WORK(KRHO2D),LWORK-KRHO2D)
C
           IF (LCOR .OR. LSEC) THEN
              CALL CC_CORE(WORK(KRHO1D),WORK(KRHO2D),ISYMTR)
           ENDIF
C
           IF (IPRINT.GT.2) THEN
             RHO1N = DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
             RHO2N = DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
             WRITE(LUPRI,*) 'Norm of RHO1: ',RHO1N,'aibj=',I
             WRITE(LUPRI,*) 'Norm of RHO2: ',RHO2N,'aibj=',I
           ENDIF
           IF (IPRINT.GT.125) THEN
            CALL CC_PRP(WORK(KRHO1D),WORK(KRHO2D),1,1,1)
           ENDIF
C
           CALL DAXPY(NT1AMX,-1.0D00,WORK(KRHO1),1,WORK(KRHO1D),1)
           CALL DAXPY(NT2AMX,-1.0D00,WORK(KRHO2),1,WORK(KRHO2D),1)
           CALL DSCAL(NT1AMX,DELTAI,WORK(KRHO1D),1)
           CALL DSCAL(NT2AMX,DELTAI,WORK(KRHO2D),1)
           CALL DCOPY(NT1AMX,WORK(KRHO1D),1,
     *              WORK(KF2+NTAMP*(I-1)),1)
           CALL DCOPY(NT2AMX,WORK(KRHO2D),1,
     *              WORK(KF2+NTAMP*(I-1)+NT1AMX),1)
C
           X12 = X12 + DDOT(NT1AMX,WORK(KRHO1D),1,WORK(KRHO1D),1)
           X22 = X22 + DDOT(NT2AMX,WORK(KRHO2D),1,WORK(KRHO2D),1)
           TI   = SECOND() - TI
           IF (IPRINT.GT.5 ) THEN
              WRITE(LUPRI,*) '  '
              WRITE(LUPRI,*) 'FDG ROW NR. ',I+NT1AMX,
     *                  ' DONE IN ',TI,' SEC.'
           ENDIF
C
         ENDIF
C
 300    CONTINUE
 200  CONTINUE
C
      WRITE(LUPRI,*) '    '
      WRITE(LUPRI,*) '**  FINITE DIFF WITH DELTA ',DELTA, '**'
      WRITE(LUPRI,*) '    '
      IF ((IPRINT .GT. 4)) THEN
         CALL AROUND( 'FINITE DIFF. CC K*Tx*Ty-Matrix - 11 & 21 PART ')
         CALL OUTPUT(WORK(KF),1,NTAMP,1,NC1VEC,NTAMP,NC1VEC,1,LUPRI)
         CALL AROUND( 'FINITE DIFF. CC K*Tx*Ty-Matrix - 12 & 22 PART ')
         CALL OUTPUT(WORK(KF+NTAMP*NC1VEC),1,NTAMP,1,NC2VEC,
     *               NTAMP,NC2VEC,1,LUPRI)
      ENDIF
      IF (.TRUE.) THEN
         XNJ = X11 + X12 + X21 + X22
         WRITE(LUPRI,*) '  '
         WRITE(LUPRI,*) ' NORM OF FIN. DIFF. K*tx*ty-Matrix.', SQRT(XNJ)
         WRITE(LUPRI,*) '  '
         WRITE (LUPRI,*) ' NORM OF 11 PART OF FD. K*tx*ty-mat.: ',
     &        SQRT(X11)
         WRITE (LUPRI,*) ' NORM OF 21 PART OF FD. K*tx*ty-mat.: ',
     &        SQRT(X21)
         WRITE (LUPRI,*) ' NORM OF 12 PART OF FD. K*tx*ty-mat.: ',
     &        SQRT(X12)
         WRITE (LUPRI,*) ' NORM OF 22 PART OF FD. K*tx*ty-mat.: ',
     &        SQRT(X22)
      ENDIF
C
C--------------------------------------
C     Calculate Matrix times Tz vector.
C--------------------------------------
C
      CALL DGEMV('N',NTAMP,NTAMP,ONE,WORK(KF),NTAMP,TZAM,1,
     *           ZERO,RESULT,1)
C
C-------------------------------------------------
C     Restore the CC reference amplitudes on disk.
C-------------------------------------------------
C
      LUTAM = -1
      CALL GPOPEN(LUTAM,'TAM_SAV','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *            .FALSE.)
      REWIND(LUTAM)
      READ(LUTAM) (WORK(KC1AM + I -1 ) , I = 1, NT1AMX)
      READ(LUTAM) (WORK(KC2AM + I -1 ) , I = 1, NT2AMX)
      CALL GPCLOSE(LUTAM,'DELETE')
C
      IOPT = 3
      CALL CC_WRRSP('R0',0,1,IOPT,MODEL,DUMMY,WORK(KC1AM),
     *              WORK(KC2AM),WORK(KEND2),LWRK2)
C
      RSPIM = .TRUE.
      CALL CCRHSN(WORK(KRHO1D),WORK(KRHO2D),WORK(KC1AM),
     &            WORK(KC2AM),
     *            WORK(KEND2),LWRK2,'XXX') 
C
      IF (IPRINT .GT. 10) THEN
         CALL AROUND(' END OF CC_FDH ')
      ENDIF
C
      RETURN
      END
*=====================================================================*
*=====================================================================*
      SUBROUTINE CCSDT_HMAT_NODDY(LISTL,IDLSTL,LISTB,IDLSTB,
     &                            LISTC,IDLSTC,LISTD,IDLSTD,
     &                            OMEGA1,OMEGA2,WORK,LWORK)
*---------------------------------------------------------------------*
*
*    Purpose: compute triples contribution to G matrix transformation
*
*  (H T^B T^C T^D)^eff_1,2 = (H T^B T^C T^D)_1,2(CCSD) 
*                             + (H T^B T^C)_1,2(L_3)
*        
*     Written by Christof Haettig, April 2002 
*     based on CCSDT_GMAT_NODDY
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdinp.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccfield.h>
#include <ccorb.h>
#include <dummy.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.false.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      CHARACTER*3 LISTL, LISTB, LISTC, LISTD
      INTEGER LWORK, IDLSTL, IDLSTB, IDLSTC, IDLSTD

#if defined (SYS_CRAY)
      REAL WORK(LWORK), TWO, FREQL, DDOT
      REAL OMEGA1(NT1AMX), OMEGA2(NT2AMX)
#else
      DOUBLE PRECISION WORK(LWORK), TWO, FREQL, DDOT
      DOUBLE PRECISION OMEGA1(NT1AMX), OMEGA2(NT2AMX)
#endif
      PARAMETER ( TWO = 2.0D0 )

      CHARACTER*10 MODEL
      INTEGER KXINT, KXIAJB, KYIAJB, KT1AMP0, KLAMP0, KLAMH0, KEND1,
     &        KFOCK0, LWRK1, KLAMPB, KLAMHB, KLAMPC, KLAMHC, KT1AMPB,
     &        KLAMPD, KLAMHD, KT1AMPD, KT1AMPC, 
     &        KINT1T0, KINT2T0, 
     &        KINT1SBCD, KINT2SBCD,
     &        KBDIOVVO, KBDIOOVV, KBDIOOOO, KBDIVVVV,
     &        KCDIOVVO, KCDIOOVV, KCDIOOOO, KCDIVVVV,
     &        KBCIOVVO, KBCIOOVV, KBCIOOOO, KBCIVVVV,
     &        KOME1, KOME2, KL1AM, KL2AM, KL3AM, KT2AM, KFOCKD,
     &        KSCR1, KDUM
      INTEGER ISYMD, ILLL, IDEL, ISYDIS, IJ, NIJ, LUSIFC, INDEX,
     &        ISYMC, ISYMB, LUFOCK, KEND2, LWRK2, IOPT, ILSTSYM, ISYML

      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J 

      KDUM = 1
      CALL QENTER('CCSDT_HMAT_NODDY')

      IF (DIRECT) CALL QUIT('CCSDT_HMAT_NODDY: DIRECT NOT IMPLEMENTED')

*---------------------------------------------------------------------*
*     Memory allocation:
*---------------------------------------------------------------------*
      KSCR1   = 1
      KFOCKD  = KSCR1  + NT1AMX
      KEND1   = KFOCKD + NORBT

      KFOCK0  = KEND1
      KEND1   = KFOCK0  + NORBT*NORBT

      KT1AMP0 = KEND1
      KLAMP0  = KT1AMP0 + NT1AMX
      KLAMH0  = KLAMP0  + NLAMDT
      KEND1   = KLAMH0  + NLAMDT

      KT1AMPB = KEND1
      KLAMPB  = KT1AMPB + NT1AMX
      KLAMHB  = KLAMPB  + NLAMDT
      KEND1   = KLAMHB  + NLAMDT

      KT1AMPC = KEND1
      KLAMPC  = KT1AMPC + NT1AMX
      KLAMHC  = KLAMPC  + NLAMDT
      KEND1   = KLAMHC  + NLAMDT

      KT1AMPD = KEND1
      KLAMPD  = KT1AMPD + NT1AMX
      KLAMHD  = KLAMPD  + NLAMDT
      KEND1   = KLAMHD  + NLAMDT

      KXIAJB  = KEND1
      KYIAJB  = KXIAJB  + NT1AMX*NT1AMX
      KEND1   = KYIAJB  + NT1AMX*NT1AMX

      KINT1T0 = KEND1
      KINT2T0 = KINT1T0 + NT1AMX*NVIRT*NVIRT
      KEND1   = KINT2T0 + NRHFT*NRHFT*NT1AMX

      KBDIOVVO = KEND1
      KBDIOOVV = KBDIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KBDIOOOO = KBDIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KBDIVVVV = KBDIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1    = KBDIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KCDIOVVO = KEND1
      KCDIOOVV = KCDIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KCDIOOOO = KCDIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KCDIVVVV = KCDIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1    = KCDIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KBCIOVVO = KEND1
      KBCIOOVV = KBCIOVVO + NRHFT*NVIRT*NVIRT*NRHFT
      KBCIOOOO = KBCIOOVV + NRHFT*NVIRT*NVIRT*NRHFT
      KBCIVVVV = KBCIOOOO + NRHFT*NRHFT*NRHFT*NRHFT
      KEND1    = KBCIVVVV + NVIRT*NVIRT*NVIRT*NVIRT

      KINT1SBCD = KEND1
      KINT2SBCD = KINT1SBCD + NT1AMX*NVIRT*NVIRT
      KEND1     = KINT2SBCD + NRHFT*NRHFT*NT1AMX

      KOME1   = KEND1
      KOME2   = KOME1  + NT1AMX
      KEND1   = KOME2  + NT1AMX*NT1AMX

      KL1AM   = KEND1
      KL2AM   = KL1AM  + NT1AMX
      KL3AM   = KL2AM  + NT1AMX*NT1AMX
      KEND1   = KL3AM + NT1AMX*NT1AMX*NT1AMX

      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in CCSDT_HMAT_NODDY')
      ENDIF

*---------------------------------------------------------------------*
*     Read SCF orbital energies from file:
*---------------------------------------------------------------------*
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND LUSIFC
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD+I-1), I=1,NORBT)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     Get zeroth-order Lambda matrices:
*---------------------------------------------------------------------*
      IOPT   = 1
      Call CC_RDRSP('R0',0,ISYM0,IOPT,MODEL,WORK(KT1AMP0),WORK(KDUM))

      Call LAMMAT(WORK(KLAMP0),WORK(KLAMH0),WORK(KT1AMP0),
     &            WORK(KEND1),LWRK1)

*---------------------------------------------------------------------*
*     Get response  Lambda matrices:
*---------------------------------------------------------------------*
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 1
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KT1AMPB),WORK(KDUM))

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPB),
     &                 WORK(KLAMH0),WORK(KLAMHB),WORK(KT1AMPB),ISYMB)

      ISYMC = ILSTSYM(LISTC,IDLSTC)
      IOPT  = 1
      Call CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &              WORK(KT1AMPC),WORK(KDUM))

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPC),
     &                 WORK(KLAMH0),WORK(KLAMHC),WORK(KT1AMPC),ISYMC)

      ISYMD = ILSTSYM(LISTD,IDLSTD)
      IOPT  = 1
      Call CC_RDRSP(LISTD,IDLSTD,ISYMD,IOPT,MODEL,
     &              WORK(KT1AMPD),WORK(KDUM))

      CALL CCLR_LAMTRA(WORK(KLAMP0),WORK(KLAMPD),
     &                 WORK(KLAMH0),WORK(KLAMHD),WORK(KT1AMPD),ISYMD)

*---------------------------------------------------------------------*
*     read zeroth-order AO Fock matrix from file: 
*---------------------------------------------------------------------*
      LUFOCK = -1
      CALL GPOPEN(LUFOCK,'CC_FCKH','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUFOCK)
      READ(LUFOCK) (WORK(KFOCK0-1+I),I=1,N2BST(ISYM0))
      CALL GPCLOSE(LUFOCK,'KEEP')

      CALL CC_FCKMO(WORK(KFOCK0),WORK(KLAMP0),WORK(KLAMH0),
     &              WORK(KEND1),LWRK1,ISYM0,ISYM0,ISYM0)      

*---------------------------------------------------------------------*
*     Compute integrals needed for the following contributions:
*---------------------------------------------------------------------*

      CALL DZERO(WORK(KBDIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBDIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBDIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KBDIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KCDIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KCDIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KCDIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KCDIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KBCIOVVO),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBCIOOVV),NRHFT*NVIRT*NVIRT*NRHFT)
      CALL DZERO(WORK(KBCIOOOO),NRHFT*NRHFT*NRHFT*NRHFT)
      CALL DZERO(WORK(KBCIVVVV),NVIRT*NVIRT*NVIRT*NVIRT)

      CALL DZERO(WORK(KXIAJB),NT1AMX*NT1AMX)
      CALL DZERO(WORK(KYIAJB),NT1AMX*NT1AMX)

      CALL DZERO(WORK(KINT1T0),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2T0),NRHFT*NRHFT*NT1AMX)

      CALL DZERO(WORK(KINT1SBCD),NT1AMX*NVIRT*NVIRT)
      CALL DZERO(WORK(KINT2SBCD),NRHFT*NRHFT*NT1AMX)

      DO ISYMD = 1, NSYM
         DO ILLL = 1,NBAS(ISYMD)
            IDEL   = IBAS(ISYMD) + ILLL
            ISYDIS = MULD2H(ISYMD,ISYMOP)
 
C           ----------------------------
C           Work space allocation no. 2.
C           ----------------------------
            KXINT  = KEND1
            KEND2  = KXINT + NDISAO(ISYDIS)
            LWRK2  = LWORK - KEND2
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
               CALL QUIT('Insufficient space in CCSDT_HMAT_NODDY')
            ENDIF
 
C           ---------------------------
C           Read in batch of integrals.
C           ---------------------------
            CALL CCRDAO(WORK(KXINT),IDEL,1,WORK(KEND2),LWRK2,
     *                  WORK(KDUM),DIRECT)
 
C           ----------------------------------
C           Calculate integrals needed in CC3:
C           ----------------------------------
            CALL CCSDT_TRAN1(WORK(KINT1T0),WORK(KINT2T0),
     &                       WORK(KLAMP0),WORK(KLAMH0),
     &                       WORK(KXINT),IDEL)

            CALL CC3_TRAN2(WORK(KXIAJB),WORK(KYIAJB),
     &                     WORK(KLAMP0),WORK(KLAMH0),
     &                     WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barB K-barC|B-barD D)
            ! XINT2S = XINT2S + (C-barB K-barC|L J-barD)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBCD),WORK(KINT2SBCD),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHC),
     &                         WORK(KLAMPD),WORK(KLAMHD),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barB K-barD|B-barC D)
            ! XINT2S = XINT2S + (C-barB K-barD|L J-barC)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBCD),WORK(KINT2SBCD),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHD),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barC K-barB|B-barD D)
            ! XINT2S = XINT2S + (C-barC K-barB|L J-barD)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBCD),WORK(KINT2SBCD),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHB),
     &                         WORK(KLAMPD),WORK(KLAMHD),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barD K-barB|B-barC D)
            ! XINT2S = XINT2S + (C-barD K-barB|L J-barC)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBCD),WORK(KINT2SBCD),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPD),WORK(KLAMHB),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barC K-barD|B-barB D)
            ! XINT2S = XINT2S + (C-barC K-barD|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBCD),WORK(KINT2SBCD),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHD),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

            ! XINT1S = XINT1S + (C-barD K-barC|B-barB D)
            ! XINT2S = XINT2S + (C-barD K-barC|L J-barB)
            CALL CCSDT_TRAN3_R(WORK(KINT1SBCD),WORK(KINT2SBCD),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPD),WORK(KLAMHC),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

C           ----------------------------------------------
C             (OV|VO)-BD (OO|VV)-BD (OO|OO)-BD (VV|VV)-BD
C           ----------------------------------------------
            CALL CCFOP_TRAN1_R(WORK(KBDIOVVO),WORK(KBDIOOVV),
     &                         WORK(KBDIOOOO),WORK(KBDIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KLAMPD),WORK(KLAMHD),
     &                         WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(KBDIOVVO),WORK(KBDIOOVV),
     &                         WORK(KBDIOOOO),WORK(KBDIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPD),WORK(KLAMHD),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

C           ----------------------------------------------
C             (OV|VO)-CD (OO|VV)-CD (OO|OO)-CD (VV|VV)-CD
C           ----------------------------------------------
            CALL CCFOP_TRAN1_R(WORK(KCDIOVVO),WORK(KCDIOOVV),
     &                         WORK(KCDIOOOO),WORK(KCDIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KLAMPD),WORK(KLAMHD),
     &                         WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(KCDIOVVO),WORK(KCDIOOVV),
     &                         WORK(KCDIOOOO),WORK(KCDIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPD),WORK(KLAMHD),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

C           ----------------------------------------------
C           (OV|VO)-BC  (OO|VV)-BC  (OO|OO)-BC  (VV|VV)-BC
C           ----------------------------------------------
            CALL CCFOP_TRAN1_R(WORK(KBCIOVVO),WORK(KBCIOOVV),
     &                         WORK(KBCIOOOO),WORK(KBCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KXINT),IDEL)

            CALL CCFOP_TRAN1_R(WORK(KBCIOVVO),WORK(KBCIOOVV),
     &                         WORK(KBCIOOOO),WORK(KBCIVVVV),
     &                         WORK(KLAMP0),WORK(KLAMH0),
     &                         WORK(KLAMPC),WORK(KLAMHC),
     &                         WORK(KLAMPB),WORK(KLAMHB),
     &                         WORK(KXINT),IDEL)

         END DO   
      END DO  

*---------------------------------------------------------------------*
*     Compute L^0_3 multipliers:
*---------------------------------------------------------------------*
      ISYML = ILSTSYM(LISTL,IDLSTL)

      CALL DZERO(WORK(KL3AM),NT1AMX*NT1AMX*NT1AMX)

      IF (LISTL(1:3).EQ.'L0 ') THEN

        IF (LWRK1 .LT. NT2AMX) 
     *     CALL QUIT('Insufficient space in CCSDT_HMAT_NODDY')

        ! read left amplitudes from file and square up the doubles part
        IOPT   = 3
        Call CC_RDRSP(LISTL,IDLSTL,ISYML,IOPT,MODEL,
     *                WORK(KL1AM),WORK(KEND1))
        CALL CC_T2SQ(WORK(KEND1),WORK(KL2AM),ISYML)

        CALL CCSDT_L03AM(WORK(KL3AM),WORK(KINT1T0),WORK(KINT2T0),
     *                  WORK(KXIAJB),WORK(KFOCK0),WORK(KL1AM),
     *                  WORK(KL2AM),WORK(KSCR1),WORK(KFOCKD),
     *                  DUMMY,DUMMY)

      ELSE IF (LISTL(1:3).EQ.'L1 ' .OR. LISTL(1:3).EQ.'LE ' .OR.
     &         LISTL(1:3).EQ.'M1 ' .OR. LISTL(1:3).EQ.'N2 ' .OR.
     &         LISTL(1:3).EQ.'E0 '                              ) THEN

        CALL CCSDT_TBAR31_NODDY(WORK(KL3AM),FREQL,LISTL,IDLSTL,
     &                        WORK(KLAMP0),WORK(KLAMH0),
     &                        WORK(KFOCK0),WORK(KFOCKD),WORK(KSCR1),
     &                        WORK(KXIAJB),WORK(KINT1T0),WORK(KINT2T0),
     &                        WORK(KEND1),LWRK1)

      ELSE

        CALL QUIT('CCSDT_HMAT_NODDY> LISTL NOT AVAILABLE:'//LISTL)
      
      END IF

      CALL DSCAL(NT1AMX*NT1AMX*NT1AMX,-1.0D0,WORK(KL3AM),1)

      IF (LOCDBG) THEN
        WRITE(LUPRI,*) 'NORM^2(L3AM):',
     *    DDOT(NT1AMX*NT1AMX*NT1AMX,WORK(KL3AM),1,WORK(KL3AM),1)
      END IF

*---------------------------------------------------------------------*
*     Compute contribution from  <L_3|[[H^BC, T^D_2],\tau_nu_1]|HF>
*                          and   <L_3|[[H^BD, T^C_2],\tau_nu_1]|HF>
*                          and   <L_3|[[H^CD, T^B_2],\tau_nu_1]|HF>
*---------------------------------------------------------------------*
      KT2AM  = KEND1
      KEND1  = KT2AM + NT1AMX*NT1AMX
      LWRK1  = LWORK  - KEND1
      IF (LWRK1 .LT. NT2AMX) THEN
         CALL QUIT('Insufficient space in CCSDT_HMAT_NODDY')
      ENDIF

      CALL DZERO(WORK(KOME1),NT1AMX)

      ! read T^D doubles amplitudes from file and square up 
      ISYMD = ILSTSYM(LISTD,IDLSTD)
      IOPT  = 2
      Call CC_RDRSP(LISTD,IDLSTD,ISYMD,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      Call CCLR_DIASCL(WORK(KEND1),TWO,ISYMD)
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMD)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KBCIOOOO),WORK(KBCIOVVO),
     &                         WORK(KBCIOOVV),WORK(KBCIVVVV),
     &                         WORK(KT2AM))


      ! read T^C doubles amplitudes from file and square up 
      ISYMC = ILSTSYM(LISTC,IDLSTC)
      IOPT  = 2
      Call CC_RDRSP(LISTC,IDLSTC,ISYMC,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      Call CCLR_DIASCL(WORK(KEND1),TWO,ISYMC)
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMC)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KBDIOOOO),WORK(KBDIOVVO),
     &                         WORK(KBDIOOVV),WORK(KBDIVVVV),
     &                         WORK(KT2AM))


      ! read T^B doubles amplitudes from file and square up 
      ISYMB = ILSTSYM(LISTB,IDLSTB)
      IOPT  = 2
      Call CC_RDRSP(LISTB,IDLSTB,ISYMB,IOPT,MODEL,
     &              WORK(KDUM),WORK(KEND1))
      Call CCLR_DIASCL(WORK(KEND1),TWO,ISYMB)
      CALL CC_T2SQ(WORK(KEND1),WORK(KT2AM),ISYMB)

      CALL CC3_L3_OMEGA1_NODDY(WORK(KOME1),WORK(KL3AM),
     &                         WORK(KCDIOOOO),WORK(KCDIOVVO),
     &                         WORK(KCDIOOVV),WORK(KCDIVVVV),
     &                         WORK(KT2AM))

      DO I = 1,NT1AMX
         OMEGA1(I) = OMEGA1(I) + WORK(KOME1+I-1)
      END DO

*---------------------------------------------------------------------*
*     Compute contribution from  <L_3|[H^BCD,\tau_nu_2]|HF>
*---------------------------------------------------------------------*
      CALL DZERO(WORK(KOME2),NT1AMX*NT1AMX)

      CALL CC3_L3_OMEGA2_NODDY(WORK(KOME2),WORK(KL3AM),
     *                         WORK(KINT1SBCD),WORK(KINT2SBCD))

      DO I = 1,NT1AMX
         DO J = 1,I
            IJ = NT1AMX*(I-1) + J
            NIJ = INDEX(I,J)
            OMEGA2(NIJ) = OMEGA2(NIJ) + WORK(KOME2+IJ-1)
         END DO
      END DO

      CALL QEXIT('CCSDT_HMAT_NODDY')

      RETURN
      END 
*---------------------------------------------------------------------*
*              END OF SUBROUTINE CCSDT_HMAT_NODDY                     *
*---------------------------------------------------------------------*
