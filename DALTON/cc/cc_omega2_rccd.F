!  /* Deck cc_omega2_rccd */

      SUBROUTINE cc_omega2_rccd(model,omega1,omega2,work,lwork)
!
!     Purpose: generate all contributions to the Omega_aibj function
!     in RCCD, DRCCD and RTCCD to test that everything works.
!     It involves construction of (ia|jb) (from file), (ia|bj) and (ij|ab)
!     omega1 and omega2 are output. Omega1=0
!     Omega function expression for RCCD as derived by Thomas B. Pedersen
!     Omega function expression for RTCCD from triplet contributions
!     from Szabo and Ostlund, JCP ...
!     Sonia & FRAN, nov 2009, RTCCD by Sonia, Nov 2010
!
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, THRDEM
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, FOUR=4.0d0)
      PARAMETER(HALF = 0.5D0)
      PARAMETER (THRDEM = 1.0D-12)
      CHARACTER MODEL*10
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION WORK(LWORK)
      DIMENSION OMEGA1(*),OMEGA2(*)
      integer del, delI, idel
      LOGICAL LOCDBG
      parameter (LOCDBG = .false.)
#include "ccorb.h"
#include "ccisao.h"
#include "r12int.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccinftap.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eritap.h"
#include "cclr.h"
#include "ccfro.h"
#include "ccfield.h"
#include "dftcom.h"
!#include "oepopt.h"
!#include "ccandy.h"
      logical FOCKMAT
      parameter (FOCKMAT=.true.)
      LOGICAL ETRAN,FCKCON
!
      CALL QENTER('CC_OMEGA2_RCCD')
!
C      CALL HEADER('Noddy code for Omega_2 in cc_RCCD',-1)
      if (LOCDBG) then
        write(lupri,*)'OMEGA2 FOR RCCD: NONHF=', NONHF
        write(lupri,*)'OMEGA2 FOR RCCD: NHFFIELD=',NHFFIELD
        write(lupri,*)'OMEGA2 FOR RCCD: NFIELD=', NFIELD
      end if
!
!------------------------------------------------------------------
!     NO SYMMETRY
!------------------------------------------------------------------
!
      ISYMTR = 1
      ISYMOP = 1
!
! Set Omega1=0, we don't need it.
!
      CALL DZERO(OMEGA1,NT1AMX)
      CALL DZERO(OMEGA2,NT2AMX)
!
!------------------------------------------------------------------
!     ALLOCATE FOR INDIVIDUAL CONTRIBUTIONS to Omega2
!------------------------------------------------------------------

      KOMEGAF = 1
      KOMEGAE = KOMEGAF + NT2AMX
      KOMEGAD = KOMEGAE + NT2SQ(1)
      KOMEGAC = KOMEGAD + NT2SQ(1)
      KOMEGA2 = KOMEGAC + NT2SQ(1)
      KFOCKD  = KOMEGA2 + NT2SQ(1)
      KOMEpk  = KFOCKD  + NORBTS
      Komepk2 = KOMEpk +  NT2AMX
      KEND1   = KOMEpk2 + NT2AMX
      LWRK1   = LWORK - KEND1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 1.st allocation '//
     &                'in CC_RCCD_NODDY')
      ENDIF
      call dzero(WORK(KOMEGA2),Nt2sq(1))
      call dzero(WORK(KOMEGAF),NT2AMX)
      call dzero(WORK(KOMEPK),NT2AMX)
      call dzero(WORK(KOMEPK2),NT2AMX)
!-------------------------------------
!     Read canonical orbital energies.
!-------------------------------------
!
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUSIFC)
!
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD + I - 1), I = 1,NORBTS)
!
      CALL GPCLOSE(LUSIFC,'KEEP')
!
!----------------------------------------------------------------
!     Change symmetry ordering of the canonical orbital energies.
!----------------------------------------------------------------
!
      IF (FROIMP)
     *    CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
      CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
      
      if (locdbg) then
        do i=1,NORBTS
          write(lupri,*) 'Epsilon_',i,' = ', WORK(KFOCKD+i-1)
        end do
      end if
!
!-----------------------------
! F term is: i
!   L_aibj  for RCCD
!   g_ajbi  for RTCCD
!   2g_aibj for dRCCD
!-----------------------------
!
      KIAJB  = KEND1
      KLIAJB = KIAJB  + NT2AMX
      KEND1  = KLIAJB + NT2AMX
      LWRK1  = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 1.st allocation '//
     &                'in CC_RCCD_NODDY')
      ENDIF
      call dzero(WORK(KIAJB),NT2AMX)
      REWIND(LUIAJB)
      READ(LUIAJB) (WORK(KIAJB + I - 1), I = 1,NT2AMX)
!      
      IF (IPRINT .GT.45) THEN
        write(lupri,*) ' IAJB integrals (pck)'
        CALL CC_PRP(WORK(KEND1),WORK(KIAJB),1,0,1)
      ENDIF
!
      call DCOPY(NT2AMX,WORK(KIAJB),1,WORK(KLIAJB),1)
!
      IF (RCCD) then
         !L_iajb, i.e. 2C-E'
         IOPT = 1
         CALL CCSD_TCMEPK(WORK(KLIAJB),1.0d0,1,IOPT)
      else if (RTCCD) then
         IOPT = 2
         CALL CCSD_TCMEPK(WORK(KLIAJB),1.0d0,1,IOPT)
         !call DSCAL(NT2AMX,-1.0d0,WORK(KLIAJB),1)
      else
         !generate 2*g_iajb
         call DSCAL(NT2AMX,2.0d0,WORK(KLIAJB),1)
      ENDIF
!
      IF ((IPRINT.GT.45).or.LOCDBG) THEN
        write(lupri,*) 'L_IAJB integrals (pck) in RCCD'
        write(lupri,*) '2*G_IAJB integrals (pck) in DRCCD'
        write(lupri,*) 'G_IBJA integrals (pck) in RTCCD'
        CALL CC_PRP(WORK(KEND1),WORK(KLIAJB),1,0,1)
      END IF
!
      if ((RCCD).or.(DRCCD)) then 
        CALL DAXPY(NT2AMX,1.0d0,WORK(KLIAJB),1,WORK(KOMEGAF),1)
      else if (RTCCD) then
        !SONIA FIXME change sign because of Szabo-Ostlund
        !CALL DAXPY(NT2AMX,-1.0d0,WORK(KLIAJB),1,WORK(KOMEGAF),1)
        CALL DAXPY(NT2AMX,1.0d0,WORK(KLIAJB),1,WORK(KOMEGAF),1)
      end if
!      
      IF (IPRINT.GT.45) THEN
        write(lupri,*)'Omega of RCCD methods: F contrib (pck)'
        CALL CC_PRP(WORK(KEND1),WORK(KOMEGAF),1,0,1)
      ENDIF
!
!--------------------------------------------------------------------
!     Read from file the T2 amplitudes 
!--------------------------------------------------------------------
!
      KT2AM   = KEND1
      KT2SQ   = KT2AM + NT2AMX
      KEND1   = KT2SQ + NT2SQ(1)
      LWRK1   = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 2nd allocation '//
     &                'in CC_OMEGA2_RCCD')
      ENDIF
      IOPT = 2
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KEND1),WORK(KT2AM))
      CALL DZERO(WORK(KT2SQ),NT2SQ(1))
      CALL CC_T2SQ(WORK(KT2AM),WORK(KT2SQ),1)
!
      KT1AM = KEND1
      KEND1 = KT1AM + NT1AMX
      KLAMDP = KEND1
      KLAMDH = KLAMDP + NLAMDT
      KEND1 = KLAMDH + NLAMDT
      LWRK1  = LWORK  - KEND1
!
      CALL DZERO(WORK(KT1AM),NT1AMX)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *               LWRK1)
!-------------------------------------
! Compute E term
!-------------------------------------

      call RCCD_E(work(KFockD),WORK(KT2SQ),WORK(KOMEGAE),
     &                WORK(KEND1),LWRK1,WORK(KT2SQ),WORK(komepk2))
!     
      IF ((IPRINT.GT.45).or.(LOCDBG)) THEN
        WRITE(LUPRI,*) 'cc_omega2_rccd: the E term'
        call cc_prsq(work(kend1),WORK(KOMEGAE),1,0,1) 
      ENDIF
      iopt=0
      call cc_t2pk(WORK(komepk),WORK(KOMEGAE),1,iopt)
      
      if (iprint.gt.45) then
         WRITE(LUPRI,*) 'The E term packed'
         call cc_prp(work(kend1),WORK(KOMEPK),1,0,1) 
         WRITE(LUPRI,*) 'The FF contribution to E term packed'
         call cc_prp(work(kend1),WORK(KOMEPK2),1,0,1) 
      end if
!--------------------------------------------------------------------
!     Recover the full MO coefficient matrix 
!--------------------------------------------------------------------

      IF ((RCCD).or.(RTCCD)) THEN
        !compute the exchange intergrals needed. Requires MO mat
        KCMO  = KEND1
        KEND1 = KCMO  + NLAMDS
        LWRK1 = LWORK - KEND1
!        
        IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation '//
     &                'in RCCD_NODDY_tb')
        ENDIF
!        
        CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)

!--------------------------------------------------------------------
!     Allocate for integrals
!--------------------------------------------------------------------
        KIOOVV = KEND1
        KEND1  = KIOOVV + NRHFT*NRHFT*NVIRT*NVIRT
        KIaick = KEND1
        KEND1  = KIaick + NT2SQ(1)
        LWRK1  = LWORK - KEND1
!        
        IF (LWRK1 .LT. 0) THEN
           WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
           CALL QUIT('Insufficient memory for allocation '//
     &                'in cc_omega2_rccd')
        ENDIF
!
        call dzero(WORK(KIOOVV),NRHFT*NRHFT*NVIRT*NVIRT)
        call dzero(WORK(KIaick),NT2SQ(1))

!------------------------------------------------------------------
! Start loop on integral distributions to build (OO|VV)
!------------------------------------------------------------------

        KENDS2 = KEND1
        LWRKS2 = LWRK1
!
        IF (DIRECT) THEN
         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
           KCCFB1 = KEND1
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK1  = LWORK  - KEND1
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                 KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                 WORK(KEND1),LWRK1,IPRERI)
           KEND1 = KFREE
           LWRK1 = LFREE
         END IF
         NTOSYM = 1
        ELSE
         NTOSYM = NSYM
        ENDIF
!
        KENDSV = KEND1
        LWRKSV = LWRK1
!
        ICDEL1 = 0

        DO 100 ISYMD1 = 1,NTOSYM
!
         IF (DIRECT) THEN
            IF (HERDIR) THEN
              NTOT = MAXSHL
            ELSE
              NTOT = MXCALL
            END IF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
!
         DO 110 ILLL = 1,NTOT
!
!---------------------------------------------
!           If direct calculate the integrals.
!---------------------------------------------
!
            IF (DIRECT) THEN
!
               KEND1 = KENDSV
               LWRK1 = LWRKSV
!
               IF (HERDIR) THEN
                 CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                       IPRERI)
               ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                       WORK(KODCL1),WORK(KODCL2),
     *                       WORK(KODBC1),WORK(KODBC2),
     *                       WORK(KRDBC1),WORK(KRDBC2),
     *                       WORK(KODPP1),WORK(KODPP2),
     *                       WORK(KRDPP1),WORK(KRDPP2),
     *                       WORK(KCCFB1),WORK(KINDXB),
     *                       WORK(KEND1), LWRK1,IPRERI)
               END IF
!
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'CC_OMEGA2_RCCD')
               END IF
!
            ELSE
               NUMDIS = 1
               KRECNR = KENDSV
            ENDIF
!
!-----------------------------------------------------
!           Loop over number of distributions in disk.
!-----------------------------------------------------
!
            DO 120 IDEL2 = 1,NUMDIS
!
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
!CN                  ISYMD = ISAO(IDEL)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
!                  WRITE(LUPRI,*)'DIRECT CASE ISYMD',ISYMD
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
!                  WRITE(LUPRI,*)'NON-DIRECT CASE ISYMD',ISYMD
               ENDIF
!
!----------------------------------------
!              Work space allocation two.
!----------------------------------------
!
               ISYDIS = MULD2H(ISYMD,ISYMOP)
!
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
!
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in CC_OMEGA2_RCCD')
               ENDIF
!--------------------------------------------
!              Read AO integral distribution 
!              alpha>beta,gamma (^delta) da file.
!--------------------------------------------
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
!     
!-------------------------------------------
!                 Make oovv
!-------------------------------------------
!
!               IF(LOCDBG) THEN
!                  !Stupid and SLOW way, storage is OOVV
!                  CALL MAKE_OOVV(WORK(KIOOVV),WORK(KCMO),
!     *            WORK(KXINT),IDEL)
!               ELSE
                 !Faster way, Storage is VOOV
                  CALL RCCD_2O2V(WORK(KXINT),1, WORK(KIOOVV),
     *                       WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                       WORK(KLAMDP),1,WORK(KLAMDH),1,1,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
                                      
!               ENDIF
  120       CONTINUE
  110     CONTINUE
  100   CONTINUE

        if (IPRINT.GE.45) then
!          if (LOCDBG) then
!           write(lupri,*)'----- Integrals OOVV ----------------'
!           CALL OUTPUT(WORK(KIOOVV),1,NRHFT*NRHFT,1,NVIRT*NVIRT,
!     &                      NRHFT*NRHFT,NVIRT*NVIRT,1,LUPRI)
!          else
           write(lupri,*)'----- Integrals VOOV----------------'
           CALL OUTPUT(WORK(KIOOVV),1,NT1AMX,1,NT1AMX,
     &                      NT1AMX,NT1AMX,1,LUPRI)
!         ENDIF
        end if
!-------------------------
! Resort I_kiac to I_ai,ck
!-------------------------
!        if (locdbg) then
!          CALL Resort_I_kiac(WORK(KIaick),WORK(KIOOVV))
!          IF (IPRINT.GT.45) THEN
!             CALL AROUND('resort I_kiac -> I_aick ????')
!             CALL OUTPUT(WORK(KIaick),1,NT1AMX,1,NT1AMX,
!     &                      NT1AMX,NT1AMX,1,LUPRI)
!          ENDIF
!        ELSE  
!---------------------------------------
! Resort I_ikac (sorted aikc) to I_ai,ck
!---------------------------------------

          CALL Resort2_I_aikc(WORK(KIaick),WORK(KIOOVV))

          IF (IPRINT.GT.45) THEN
          CALL AROUND('resort2 I_ikac(aikc) -> I_aick')
          CALL OUTPUT(WORK(KIaick),1,NT1AMX,1,NT1AMX,
     &                      NT1AMX,NT1AMX,1,LUPRI)
          ENDIF
!        ENDIF  
      END IF

!------------------
! Omega_D like term
!------------------
      ! first linear term
      KT2SQ = KEND1            !squared amplitudes
      KLXSQ = KT2SQ + NT2SQ(1) !squared integrals
      KDint = KLXSQ + NT2SQ(1)
      KEND1 = KDint + NT2SQ(1)
      LWRK1 = LWORK-KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for D allocation '//
     &                'in CC_omega2_rccd')
      ENDIF
      CALL DZERO(WORK(KDint),NT2SQ(1))
      CALL DZERO(WORK(KT2SQ),NT2SQ(1))
      CALL DZERO(WORK(KLXSQ),NT2SQ(1))

      CALL CC_T2SQ(WORK(KT2AM),WORK(KT2SQ),1)
      
      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'Squared amplitudes'
        call CC_PRSQ(WORK(KEND1),work(KT2SQ),1,0,1)
      ENDIF
      !
      if ((RCCD).or.(DRCCD)) then
         CALL CC_T2SQ(WORK(KIAJB),WORK(KLXSQ),1)
      end if
      !
      IF ((IPRINT.GT.45).or.LOCDBG) THEN
        WRITE(LUPRI,*) 'Squared integrals (=0 for RTCCD)'
        call CC_PRSQ(WORK(KEND1),work(KLXSQ),1,0,1)
      ENDIF
      !generate L_ai,kc = 2g_ai,ck-g_ac,ki (kiac resorted ai,ck)
      !generate L_ck,jb = 2g_ck,bj-g_cb,jk (jkcb resorted ck,bj)
      IF ((RCCD).or.(RTCCD)) THEN
         !ADD EXCHANGE PART TO INTEGRALS: I=g(C)-0.5g(ex)
C         WRITE(LUPRI,*)'ADDING X PART, SCALE BY HFXFAC'
         CALL DAXPY(NT2SQ(1),-HFXFAC*0.5d0,WORK(KIaick),1,WORK(KLXSQ),1)
      END IF
      !L=2*I
      call DSCAL(NT2SQ(1),2.0d0,WORK(KLXSQ),1)

      KOFF1 = IT2SQ(1,1) + KT2SQ
      KOFF2 = IT2SQ(1,1) + KLXSQ
      KOFF3 = IT2SQ(1,1) + KDint

      NTOTAI = MAX(NT1AM(1),1)
      NTOTCK = MAX(NT1AM(1),1)
      NTOTDL = MAX(NT1AM(1),1)

      CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,WORK(KOFF1),NTOTAI,WORK(KOFF2),NTOTCK,
     &            ONE,WORK(KOFF3),NTOTAI)

      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'Linear term nr 1 '
        call CC_PRSQ(WORK(KEND1),work(KDInt),1,0,1)
      ENDIF

      KOFF1 = IT2SQ(1,1) + KT2SQ
      KOFF2 = IT2SQ(1,1) + KLXSQ
      KOFF3 = IT2SQ(1,1) + KDINT

      NTOTAI = MAX(NT1AM(1),1)
      NTOTCK = MAX(NT1AM(1),1)
      NTOTBJ = MAX(NT1AM(1),1)

      CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,WORK(KOFF2),NTOTAI,WORK(KOFF1),NTOTCK,
     &            ONE,WORK(KOFF3),NTOTAI)
C     
      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'After Addition of second Linear term'
        call CC_PRSQ(WORK(KEND1),work(KDint),1,0,1)
      ENDIF

      !the Quadratic term: square up g_ajbi

      CALL CC_T2SQ(WORK(KLIAJB),WORK(KLXSQ),1)
      if (RTCCD) then
         !SONIA FIXME: change sign because of Szabo-Ostlund
         call DSCAL(NT2AMX,-1.0d0,WORK(KLIAJB),1)
      end if

      KOFF1 = IT2AM(1,1) + KT2SQ
      KOFF2 = IT2SQ(1,1) + KLXSQ
      KOFF3 = IT2SQ(1,1) + KOMEGAD
      KOFF4 = IT2SQ(1,1) + KDint

      NTOTAI = MAX(NT1AM(1),1)
      NTOTCK = MAX(NT1AM(1),1)
      NTOTDL = MAX(NT1AM(1),1)

      CALL DZERO(WORK(KOMEGAD),NT2SQ(1))

      CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,WORK(KOFF2),NTOTck,WORK(KOFF1),NTOTdl,
     &            ONE,WORK(KOFF3),NTOTdl)

!     WRITE(LUPRI,*) 'Second multiplication and add to linear terms'

      CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,WORK(KOFF1),NTOTai,WORK(KOFF3),NTOTck,
     &            ONE,WORK(KOFF4),NTOTai)

      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'OMEGAD (squared)'
        call CC_PRSQ(WORK(KEND1),work(KDint),1,0,1)
      ENDIF

      !Pack contribution into Omega2 for solver
      IOPT=0
      CALL CC_T2PK(OMEGA2,WORK(KDint),1,IOPT)
C      
      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'OMEGA2 (after D) (packed)'
        call CC_PRP(WORK(KEND1),OMEGA2(1),1,0,1)
      ENDIF

!Add E term: it is taken twice!!!
      if (.not.(FOCKMAT)) then
         CALL DAXPY(NT2AMX,2.0d0,WORK(KOMEPK),1,OMEGA2,1)
C      
        IF(IPRINT.GT.45) THEN
          WRITE(LUPRI,*) 'OMEGA2 + E (packed)'
          call CC_PRP(WORK(KEND1),OMEGA2(1),1,0,1)
        ENDIF
        if ((nfield.gt.0).and.NONHF) then
           CALL DAXPY(NT2AMX,2.0d0,WORK(KOMEPK2),1,OMEGA2,1)
          !if (iprint.ge.45) then
             WRITE(LUPRI,*) 'OMEGA2 + E_field (packed)'
             call CC_PRP(WORK(KEND1),OMEGA2(1),1,0,1)
          !end if
        end if
       end if

!Add F term
      CALL DAXPY(NT2AMX,1.0d0,WORK(KOMEGAF),1,OMEGA2,1)

      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'OMEGA2 +F (packed)'
        call CC_PRP(WORK(KEND1),OMEGA2(1),1,0,1)
      ENDIF

!Scale the diagonal????
!      WRITE(LUPRI,*) 'Scaling the diagonal'
!      CALL CCLR_DIASCL(OMEGA2,HALF,1)

!
!If FOCKMAT compute E term via the HF Fock matrix + field
!
      if (fockmat) then
         kfckHF = kend1
         kcmof  = kfckHF + N2BST(1)
         kfield = kcmof  + NLAMDS
         kt2am  = kfield + N2BST(1)
         KT2SQ  = KT2AM + NT2AMX
         KEND1  = KT2SQ + NT2SQ(1)
         lwrk1  = lwork - kend1
         CALL DZERO(WORK(KCMOF),NLAMDS)
         CALL CMO_ALL(WORK(KCMOF),WORK(KEND1),LWRK1)
         call DZERO(WORK(KField),N2BST(1))
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KEND1),WORK(KT2AM))
C         CALL AROUND( 'In CC_RHS_RCCD: The T2 for E term other way' )
C         call cc_prp(work(kend1),WORK(KT2AM),1,0,1) 
         CALL DZERO(WORK(KT2SQ),NT2SQ(1))
         CALL CC_T2SQ(WORK(KT2AM),WORK(KT2SQ),1)
!     This AO Fock matrix is constructed from the CMO transformed density
         LUFCK=-1
         CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
!
         REWIND(LUFCK)
         READ(LUFCK)(WORK(KFCKHF + I-1),I = 1,N2BST(1))
         CALL GPCLOSE(LUFCK,'KEEP' )
!
         IF (IPRINT .GT. 10) THEN
            CALL AROUND( 'HARTREE-FOCK Fock AO matrix' )
            CALL CC_PRFCKAO(WORK(KFCKHF),1)
         ENDIF
         if ((nfield.gt.0).and.(nonhf)) then
           CALL AROUND( 'CC_RHS_RCCD: ADD Finite FIELD TO FOCK!')
            DO IFI = 1, NFIELD
               FF =  EFIELD(IFI)
               CALL CC_ONEP(WORK(KFIELD),WORK(KEND1),LWRK1,FF,1,
     &              LFIELD(IFI))
            END DO
            CALL DAXPY(N2BST(1),1.0d0,WORK(KFIELD),1,WORK(KFCKHF),1)
            CALL CC_FCKMO(WORK(KFIELD),WORK(KCMOF),WORK(KCMOF),
     *                 WORK(KEND1),LWRK1,1,1,1)
            IF (IPRINT .GT. 10) THEN
              CALL AROUND( 'CC_RHS_RCCD: FF Fock MO matrix' )
              CALL CC_PRFCKMO(WORK(KFIELD),1)
            ENDIF
         end if
         !trasform to MO
         CALL CC_FCKMO(WORK(KFCKHF),WORK(KCMOF),WORK(KCMOF),
     *                 WORK(KEND1),LWRK1,1,1,1)
         IF (IPRINT .GT. 10) THEN
            CALL AROUND( 'In CC_RHS_RCCD: HF Fock MO matrix' )
            CALL CC_PRFCKMO(WORK(KFCKHF),1)
         ENDIF

         KEMAT1 = kend1
         KEMAT2 = KEMAT1 + NMATAB(1)
         KEpack = KEMAT2 + NMATIJ(1)
         kend1  = kEpack + nt2amx
         lwrk1  = lwork-kend1
         CALL DZERO(WORK(KEMAT1),NMATAB(1))
         CALL DZERO(WORK(KEMAT2),NMATIJ(1))
         ETRAN  = .FALSE.
         FCKCON = .TRUE.
         CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),WORK(KCMOF),
     *                   WORK(KFCKHF),WORK(KEND1),LWRK1,FCKCON,
     *                   ETRAN,1)
C        WRITE (LUPRI,*) 'EMAT2:'
C        WRITE (LUPRI,'(5f12.6)') (WORK(KEMAT2-1+I),I=1,NMATIJ(1))
C        WRITE (LUPRI,*) 'EMAT1:'
C        WRITE (LUPRI,'(2f12.6)') (WORK(KEMAT1-1+I),I=1,NMATAB(1))

         CALL DZERO(WORK(KEpack),Nt2amx)
         CALL CCRHS_E(work(KEpack),WORK(KT2SQ),WORK(KEMAT1),
     &                 WORK(KEMAT2),WORK(KEND1),LWRK1,1,1)

C         CALL AROUND( 'In CC_RHS_RTCCD: The E term the other way' )
C         call cc_prp(work(kend1),WORK(KEpacK),1,0,1) 

         CALL DAXPY(NT2AMX,2.0d0,WORK(KEpacK),1,OMEGA2,1)
         !CALL DAXPY(NT2AMX,1.0d0,WORK(KEpacK),1,OMEGA2,1)
      end if
      CALL QEXIT('CC_OMEGA2_RCCD')
      RETURN
      END

!-------------------------------------------------------
!  /* Deck make_oovv */
      SUBROUTINE make_oovv(XINT,CMO,AOINT,IDEL)
!
! Generate g_oovv (XINT) from AOINT (a>b,g; delta)
!
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
      DIMENSION XINT(NRHFT,NRHFT,NVIRT,NVIRT)
      DIMENSION CMO(NBAST,NORBT)
      DIMENSION AOINT(NNBAST,NBAST)
      LOGICAL LOCDBG
!
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
!
      LOCDBG = .FALSE.

!
!----------------------------------------
!     Calculate integrals :
!----------------------------------------
!
      DO 100 G = 1,NBAST
         DO 110 IB = 1,NBAST
            DO 120 A = 1,NBAST
               NAB = INDEX(A,IB)
!
               if (aoint(nab,g) .eq. 0.0d0) goto 120
!
               DO ND = 1,NVIRT
                  DO NC = 1,NVIRT
                     DO NL = 1,NRHFT
                        DO NK = 1,NRHFT
!
                              XINT(NK,NL,NC,ND) = XINT(NK,NL,NC,ND)
     *                  + AOINT(NAB,G)*CMO(A,NK)*CMO(IB,NL)
     *                         *CMO(G,NRHFT+NC)*CMO(IDEL,NRHFT+ND)
!
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
!

      RETURN
      END
!-------------------------------------------------------
      SUBROUTINE Resort_I_kiac(Xaick,XIoovv)
!
! Resort I_kiac to I_ai,ck
!
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
      DIMENSION Xaick(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XIoovv(NRHFT,NRHFT,NVIRT,NVIRT)
      LOGICAL LOCDBG
      LOCDBG = .FALSE.
!
!----------------------------------------
!     Resort integrals :
!----------------------------------------
!
      DO c = 1,NVIRT
         DO a= 1,NVIRT
            DO i = 1,NRHFT
               DO k = 1,NRHFT
                  Xaick(A,I,C,K) = XIOOVV(K,I,A,C) 
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
!==============================================================
      SUBROUTINE RCCD_E(Epsi,T2SQ,OmegaE,Work,Lwork,T2,Omepk)
!
! E-like terms in RPA energy, including finite field contributions
! NB: T2 is also SQUARED
! Sonia
!
#include "implicit.h"
#include "priunit.h"
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccfield.h"
      DIMENSION Epsi(NORBTS),WORK(LWORK)
      DIMENSION T2SQ(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION OMEGAE(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION OMEpk(*),T2(*)
      LOGICAL LOCDBG,ETRAN,FCKCON
      LOCDBG = .FALSE.

!      write(lupri,*) 'In RCCD_E, NONHF =', NONHF,'NFIELD=', NFIELD
!----------------------------------------
!     OmegaE_aibj=T_ai,bj*epsilon_b - T_ai,bj*epsilon_j 
!                +epsilon_a*t_ai,bj - epsilon_i*t_ai,bj
!----------------------------------------
!
      do j=1,nrhft
         do b=1,nvirt
            do i = 1,nrhft
               do a=1,nvirt
                  OMEGAE(a,i,b,j) = !OMEGAE(a,i,b,j) +
     &             t2sq(a,i,b,j)*Epsi(nrhft+b)  
     &           - t2sq(a,i,b,j)*Epsi(j)  
     &           + t2sq(a,i,b,j)*Epsi(nrhft+a) 
     &           - t2sq(a,i,b,j)*Epsi(i)  
               end do
            end do
         end do
      end do

      if ((nfield.gt.0).and.(nonhf)) then
         kfock = 1
         kcmof = KFOCK  + N2BST(1)
         kend1 = kcmof  + NLAMDS
         lwrk1 = lwork-kend1
         CALL DZERO(WORK(KCMOF),NLAMDS)
         CALL CMO_ALL(WORK(KCMOF),WORK(KEND1),LWRK1)
         CALL DZERO(WORK(KFOCK),N2BST(1))
         DO IFI = 1, NFIELD
            FF =  EFIELD(IFI)
            CALL CC_ONEP(WORK(KFOCK),WORK(KEND1),LWRK1,FF,1,LFIELD(IFI))
         END DO
         CALL CC_FCKMO(WORK(KFOCK),WORK(KCMOF),WORK(KCMOF),WORK(KEND1),
     *                 LWRK1,1,1,1)
         call CC_PRFCKMO(work(kfock),1)
         KEMAT1 = kend1
         KEMAT2 = KEMAT1 + NEMAT1(1)
         Kend1  = KEMAT2 + NMATIJ(1)
         lwrk1  = lwork-kend1
         CALL DZERO(WORK(KEMAT1),NEMAT1(1))
         CALL DZERO(WORK(KEMAT2),NMATIJ(1))

         ETRAN  = .FALSE.
         FCKCON = .TRUE.
         CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),WORK(KCMOF),
     *                   WORK(KFOCK),WORK(KEND1),LWRK1,FCKCON,
     *                   ETRAN,1)
         CALL CCRHS_E(OMEPK,T2,WORK(KEMAT1),WORK(KEMAT2),
     *                WORK(KEND1),LWRK1,1,1)

      end if
      RETURN
      END
!==============================================================
      SUBROUTINE Resort2_I_aikc(Xaick,XIoovv)
!
! Resort I_ijac (stored aij,c) to I_ai,ck
! Sonia, Jan 2010
!
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
      DIMENSION Xaick(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XIoovv(NVIRT,NRHFT,NRHFT,NVIRT)
      LOGICAL LOCDBG
      LOCDBG = .FALSE.
!----------------------------------------
!     Resort integrals :
!     Rewrite using a DAXPY operation!!!!
!----------------------------------------
      DO c = 1,NVIRT
         DO k= 1,NRHFT
            DO i = 1,NRHFT
               DO a = 1,NVIRT
                  Xaick(A,I,C,K) = XIOOVV(A,I,K,C) 
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
C-------------------------------------------------------
!      SUBROUTINE Resort3_I_aikc(Xaick,XIoovv)
!C
!C Resort I_ijac (stored aij,c) to I_ai,ck
!C
!#include "implicit.h"
!#include "priunit.h"
!#include "inforb.h"
!#include "ccsdinp.h"
!#include "ccsdsym.h"
!      DIMENSION Xaick(NT1AMX,NVIRT,NRHFT)
!      DIMENSION XIoovv(NT1AMX,NRHFT,NVIRT)
!      LOGICAL LOCDBG
!      LOCDBG = .TRUE.
!C
!C----------------------------------------
!C     Resort integrals :
!C----------------------------------------
!C
!      DO c = 1,NVIRT
!         DO k= 1,NRHFT
!            KOFF1 = 
!            KOFF2 =
!            call DCOPY(NT1AMX,XIOOVV(KOFF1),1,WORK(KOFF2),1)
!            DO i = 1,NRHFT
!               DO a = 1,NVIRT
!                  Xaick(A,I,C,K) = XIOOVV(A,I,K,C) 
!               ENDDO
!            ENDDO
!         ENDDO
!      ENDDO
!      end if
!
!      RETURN
!      END

!  /* Deck rccd_2o2v */
      SUBROUTINE RCCD_2O2V(AOINT,ISYMAO,XOOVV,
     *                    XLAMP0,ISYMLP0,XLAMH0,ISYMLH0,
     *                    XLAMP1,ISYMLP1,XLAMH1,ISYMLH1,
     *                    ISYINT,WORK,LWORK,IDEL,ISYMD)
!
!     Written by Fran, based on routine by K. Hald, 2010
!     Calculate integrals that are needed for the RCCD right hand side.
!     VVOO (O=occ. V=vir.) are needed.
!     (k^Lambda_p0 l^Lambda_h1 | c^Lambda_p1 d^Lambda_h0)
!
      IMPLICIT NONE
!
      INTEGER ISYMAO, ISYMLP0, ISYMLH0, ISYMLP1, ISYMLH1
      INTEGER ISYINT, LWORK, IDEL, ISYMD
      INTEGER ISYABJ, ISYTMP, ISYAIJ, KVOO, KEND0, LWRK0, KEND1, LWRK1
      INTEGER KEND2, LWRK2, ISYMJ, ISALBE, ISYMAI, KSCR1, KSCR2
      INTEGER ISYMI, ISYMAL, ISYMBE, ISYMA, KOFF1, KOFF2, KOFF3
      INTEGER NTOTAL, NTOTA, NTOTB, ISYMB
      INTEGER ISYABG, KOOV, ISYMG, ISYMIJ, KINT, NBASAL, NBASBE
      INTEGER NRHFI, NIJ, NBASG, ISAIJD
!
#if defined (SYS_CRAY)
      REAL AOINT(*),  XOOVV(*)
      REAL XLAMP0(*), XLAMH0(*)
      REAL XLAMP1(*), XLAMH1(*)
      REAL WORK(LWORK), ZERO, ONE
#else
      DOUBLE PRECISION AOINT(*),  XOOVV(*)
      DOUBLE PRECISION XLAMP0(*), XLAMH0(*)
      DOUBLE PRECISION XLAMP1(*), XLAMH1(*)
      DOUBLE PRECISION WORK(LWORK), ZERO, ONE
#endif
C
#include "priunit.h"
#include "ccsdsym.h"
#include "ccorb.h"
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
      CALL QENTER('RCCD_2O2V')

!=========================================
!     Calculate the integrals g_{oovv}
!=========================================
!
      ISYABG = MULD2H(ISYMAO,ISYMD)
!
      ISAIJD = MULD2H(ISYINT,ISYMLH0)
      ISYAIJ = MULD2H(ISAIJD,ISYMD)
!
      KOOV  = 1
      KEND1  = KOOV  + NCKI(ISYAIJ)
      LWRK1  = LWORK - KEND1
!
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Out of memory in RCCD_2O2V (g_{oovv})')
      ENDIF
!
      CALL DZERO(WORK(KOOV),NCKI(ISYAIJ))
!
      DO ISYMG = 1, NSYM
         ISYMA  = MULD2H(ISYMG,ISYMLP1)
         ISALBE = MULD2H(ISYABG,ISYMG)
         ISYMIJ = MULD2H(ISYAIJ,ISYMA)
         ISYTMP = MULD2H(ISYMIJ,ISYMLH1)
C
         KINT   = KEND1
         KSCR1  = KINT  + NMATIJ(ISYMIJ)*NBAS(ISYMG)
         KSCR2  = KSCR1 + N2BST(ISALBE)
         KEND2  = KSCR2 + NT1AO(ISYTMP)
         LWRK2  = LWORK - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT('Out of memory in RCCD_2O2V (2)')
         ENDIF
C
         DO G = 1, NBAS(ISYMG)
C
            KOFF1 = IDSAOG(ISYMG,ISYMD) + NNBST(ISALBE)*(G-1) + 1
            CALL CCSD_SYMSQ(AOINT(KOFF1),ISALBE,WORK(KSCR1))
C
            DO ISYMJ = 1,NSYM
C
               ISYMBE = MULD2H(ISYMJ,ISYMLH1)
               ISYMAL = MULD2H(ISYMBE,ISALBE)
               ISYMI  = MULD2H(ISYMAL,ISYMLP0)
C
               KOFF1 = KSCR1 
     *               + IAODIS(ISYMAL,ISYMBE)
               KOFF2 = IGLMRH(ISYMBE,ISYMJ) + 1
               KOFF3 = KSCR2
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASBE = MAX(NBAS(ISYMBE),1)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMJ),NBAS(ISYMBE),
     *                    ONE,WORK(KOFF1),NBASAL,XLAMH1(KOFF2),NBASBE,
     *                    ZERO,WORK(KOFF3),NBASAL)
C
               KOFF1 = IGLMRH(ISYMAL,ISYMI) + 1
               KOFF2 = KSCR2
               KOFF3 = KINT 
     *               + NMATIJ(ISYMIJ)*(G - 1)
     *               + IMATIJ(ISYMI,ISYMJ)
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NRHFI  = MAX(NRHF(ISYMI),1)
C
               CALL DGEMM('T','N',NRHF(ISYMI),NRHF(ISYMJ),NBAS(ISYMAL),
     *                    ONE,XLAMP0(KOFF1),NBASAL,WORK(KOFF2),NBASAL,
     *                    ZERO,WORK(KOFF3),NRHFI)
C
            ENDDO
C
         ENDDO
C
         KOFF2 = IGLMVI(ISYMG,ISYMA)  + 1
         KOFF3 = KOOV 
     *         + IMAIJA(ISYMIJ,ISYMA)
C
         NIJ    = MAX(NMATIJ(ISYMIJ),1)
         NBASG  = MAX(NBAS(ISYMG),1)
C
         CALL DGEMM('N','N',NMATIJ(ISYMIJ),NVIR(ISYMA),NBAS(ISYMG),
     *              ONE,WORK(KINT),NIJ,XLAMP1(KOFF2),NBASG,
     *              ONE,WORK(KOFF3),NIJ)
C
      ENDDO
!
!--------------------------------------------------------------
!     Do the final contraction of delta and store in XOOVV
!     Final storage is actually O1,O2,V1,V2 --> V1O1C2,V2
!--------------------------------------------------------------
!
      ISYMB = MULD2H(ISYMD,ISYMLH0)
C
      DO B = 1, NVIR(ISYMB)
C
         KOFF1 = IGLMVI(ISYMD,ISYMB)
     *         + NBAS(ISYMD)*(B-1)
     *         + (IDEL - ibas(isymd))
C
         DO ISYMA = 1, NSYM
            ISYMIJ = MULD2H(ISYAIJ,ISYMA)
            DO ISYMI = 1, NSYM
               ISYMJ  = MULD2H(ISYMIJ,ISYMI)
               ISYMAI = MULD2H(ISYMA,ISYMI)
               DO A = 1, NVIR(ISYMA)
                  DO I = 1, NRHF(ISYMI)
C
                     KOFF2 = KOOV - 1
     *                     + IMAIJA(ISYMIJ,ISYMA)
     *                     + NMATIJ(ISYMIJ)*(A-1)
     *                     + IMATIJ(ISYMI,ISYMJ)
     *                     + I
                     KOFF3 = IT2SP(ISYAIJ,ISYMB)
     *                     + NCKI(ISYAIJ)*(B-1)
     *                     + ICKI(ISYMAI,ISYMJ)
     *                     + IT1AM(ISYMA,ISYMI)
     *                     + NVIR(ISYMA)*(I-1)
     *                     + A
C
                     CALL DAXPY(NRHF(ISYMJ),XLAMH0(KOFF1),
     *                          WORK(KOFF2),NRHF(ISYMI),
     *                          XOOVV(KOFF3),NT1AM(ISYMAI))
                  ENDDO
               ENDDO
C
            ENDDO
         ENDDO
      ENDDO
!
!--------------------------
!     End.
!--------------------------
!
      CALL QEXIT('RCCD_2O2V')
!
      RETURN
      END
!  /* Deck cc_omega2_rtccd */

      SUBROUTINE cc_omega2_rtccd(model,omega1,omega2,work,lwork)
!
!     Purpose: generate all contributions to the Omega_aibj function
!     in RTCCD 
!     It involves construction of (ia|jb) (from file), (ia|bj) and (ij|ab)
!     omega1 and omega2 are output. Omega1=0
!     'Triplet' Omega function expression as given by Szabo-Ostlund
!     (I call it closed-shell RCCD!)
!     Sonia & FRAN, nov 2009
!
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxash.h"
#include "maxorb.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, THRDEM
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, FOUR=4.0d0)
      PARAMETER(HALF = 0.5D0)
      PARAMETER (THRDEM = 1.0D-12)
      CHARACTER MODEL*10
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION WORK(LWORK)
      DIMENSION OMEGA1(*),OMEGA2(*)
      integer del, delI, idel
      LOGICAL LOCDBG
      parameter (LOCDBG = .false.)
#include "ccorb.h"
#include "ccisao.h"
#include "r12int.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "ccinftap.h"
#include "ccsdio.h"
#include "distcl.h"
#include "cbieri.h"
#include "eritap.h"
#include "cclr.h"
#include "ccfro.h"
#include "ccfield.h"
#include "dftcom.h"
      logical FOCKMAT
      parameter (FOCKMAT=.false.)
      LOGICAL ETRAN,FCKCON
!
      CALL QENTER('CC_OMEGA2_RTCCD')
!
      CALL HEADER('Noddy code for Omega_2 in cc_RTCCD',-1)
      write(lupri,*)'OMEGA2 FOR RTCCD: NONHF=', NONHF
      write(lupri,*)'OMEGA2 FOR RTCCD: NHFFIELD=',NHFFIELD
      write(lupri,*)'OMEGA2 FOR RTCCD: NFIELD=', NFIELD
!
!------------------------------------------------------------------
!     NO SYMMETRY
!------------------------------------------------------------------
!
      ISYMTR = 1
      ISYMOP = 1
!
! Set Omega1=0, we don't need it.
!
      CALL DZERO(OMEGA1,NT1AMX)
      CALL DZERO(OMEGA2,NT2AMX)
!
!------------------------------------------------------------------
!     ALLOCATE FOR INDIVIDUAL CONTRIBUTIONS to Omega2
!------------------------------------------------------------------

      KOMEGAF = 1
      KOMEGAE = KOMEGAF + NT2AMX
      KOMEGAD = KOMEGAE + NT2SQ(1)
      KOMEGAC = KOMEGAD + NT2SQ(1)
      KOMEGA2 = KOMEGAC + NT2SQ(1)
      KFOCKD  = KOMEGA2 + NT2SQ(1)
      KOMEpk  = KFOCKD  + NORBTS
      Komepk2 = KOMEpk +  NT2AMX
      KEND1   = KOMEpk2 + NT2AMX
      LWRK1   = LWORK - KEND1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 1.st allocation '//
     &                'in cc_omega2_rtccd')
      ENDIF
      call dzero(WORK(KOMEGA2),Nt2sq(1))
      call dzero(WORK(KOMEGAF),NT2AMX)
      call dzero(WORK(KOMEPK),NT2AMX)
      call dzero(WORK(KOMEPK2),NT2AMX)
!-------------------------------------
!     Read canonical orbital energies.
!-------------------------------------
!
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUSIFC)
!
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD + I - 1), I = 1,NORBTS)
!
      CALL GPCLOSE(LUSIFC,'KEEP')
!
!----------------------------------------------------------------
!     Change symmetry ordering of the canonical orbital energies.
!----------------------------------------------------------------
!
      IF (FROIMP)
     *    CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
      CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)
      
      if (locdbg) then
        do i=1,NORBTS
          write(lupri,*) 'Epsilon_',i,' = ', WORK(KFOCKD+i-1)
        end do
      end if
!
!-----------------------------
! F term  is L_aibj
!-----------------------------
!
      KIAJB  = KEND1
      KLIAJB = KIAJB  + NT2AMX
      KEND1  = KLIAJB + NT2AMX
      LWRK1  = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 1.st allocation '//
     &                'in CC_RCCD_NODDY')
      ENDIF
      call dzero(WORK(KIAJB),NT2AMX)
      REWIND(LUIAJB)
      READ(LUIAJB) (WORK(KIAJB + I - 1), I = 1,NT2AMX)
!      
      IF (IPRINT .GT.45) THEN
        write(lupri,*) ' IAJB integrals (pck)'
        CALL CC_PRP(WORK(KEND1),WORK(KIAJB),1,0,1)
      ENDIF
!
      call DCOPY(NT2AMX,WORK(KIAJB),1,WORK(KLIAJB),1)
!
      IF (RTCCD) then
         !-g_jaib, i.e. 2C-E'
         write(lupri,*)'call tcme with coulomb=0'
         IOPT = 2
         CALL CCSD_TCMEPK(WORK(KLIAJB),1.0d0,1,IOPT)
      else
         !generate 2*g_iajb
         call DSCAL(NT2AMX,2.0d0,WORK(KLIAJB),1)
      ENDIF
!
      IF(IPRINT.GT.45) THEN
        write(lupri,*) '-g_JAIB integrals (pck) in RTCCD'
        write(lupri,*) '2*G_IAJB integrals (pck) in DRCCD'
        CALL CC_PRP(WORK(KEND1),WORK(KLIAJB),1,0,1)
      END IF
!
      CALL DAXPY(NT2AMX,1.0d0,WORK(KLIAJB),1,WORK(KOMEGAF),1)
!      
      IF(IPRINT.GT.45) THEN
        write(lupri,*) 'NODDY CCRHS: F contribution (pck)'
        CALL CC_PRP(WORK(KEND1),WORK(KOMEGAF),1,0,1)
      ENDIF
!
!--------------------------------------------------------------------
!     Read from file the T2 amplitudes 
!--------------------------------------------------------------------
!
      KT2AM   = KEND1
      KT2SQ   = KT2AM + NT2AMX
      KEND1   = KT2SQ + NT2SQ(1)
      LWRK1   = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for 2nd allocation '//
     &                'in RCCD_NODDY')
      ENDIF
      IOPT = 2
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KEND1),WORK(KT2AM))
      CALL DZERO(WORK(KT2SQ),NT2SQ(1))
      CALL CC_T2SQ(WORK(KT2AM),WORK(KT2SQ),1)
!
      KT1AM = KEND1
      KEND1 = KT1AM + NT1AMX
      KLAMDP = KEND1
      KLAMDH = KLAMDP + NLAMDT
      KEND1 = KLAMDH + NLAMDT
      LWRK1  = LWORK  - KEND1
!
      CALL DZERO(WORK(KT1AM),NT1AMX)
      CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),WORK(KEND1),
     *               LWRK1)
!-------------------------------------
! Compute E term
!-------------------------------------

      call RCCD_E(work(KFockD),WORK(KT2SQ),WORK(KOMEGAE),
     &                WORK(KEND1),LWRK1,WORK(KT2SQ),WORK(komepk2))
!     
      IF ((IPRINT.GT.45).or.(LOCDBG)) THEN
        WRITE(LUPRI,*) 'RTCCD The E term'
        call cc_prsq(work(kend1),WORK(KOMEGAE),1,0,1) 
      ENDIF
      iopt=0
      call cc_t2pk(WORK(komepk),WORK(KOMEGAE),1,iopt)
      
      if (iprint.gt.45) then
         WRITE(LUPRI,*) 'RTCCD The E term packed'
         call cc_prp(work(kend1),WORK(KOMEPK),1,0,1) 
         WRITE(LUPRI,*) 'RTCCD The FF contribution to E term packed'
         call cc_prp(work(kend1),WORK(KOMEPK2),1,0,1) 
      end if
!--------------------------------------------------------------------
!     Recover the full MO coefficient matrix 
!--------------------------------------------------------------------

      IF (RTCCD) THEN
        !compute the exchange intergrals needed. Requires MO mat
        KCMO  = KEND1
        KEND1 = KCMO  + NLAMDS
        LWRK1 = LWORK - KEND1
!        
        IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation '//
     &                'in RCCD_NODDY_tb')
        ENDIF
!        
        CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)

!--------------------------------------------------------------------
!     Allocate for integrals
!--------------------------------------------------------------------
        KIOOVV = KEND1
        KEND1  = KIOOVV + NRHFT*NRHFT*NVIRT*NVIRT
        KIaick = KEND1
        KEND1  = KIaick + NT2SQ(1)
        LWRK1  = LWORK - KEND1
!        
        IF (LWRK1 .LT. 0) THEN
           WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
           CALL QUIT('Insufficient memory for allocation '//
     &                'in RCCD_NODDY_tb')
        ENDIF
!
        call dzero(WORK(KIOOVV),NRHFT*NRHFT*NVIRT*NVIRT)
        call dzero(WORK(KIaick),NT2SQ(1))

!------------------------------------------------------------------
! Start loop on integral distributions to build (OO|VV)
!------------------------------------------------------------------

        KENDS2 = KEND1
        LWRKS2 = LWRK1
!
        IF (DIRECT) THEN
         IF (HERDIR) THEN
           CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
           KCCFB1 = KEND1
           KINDXB = KCCFB1 + MXPRIM*MXCONT
           KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
           LWRK1  = LWORK  - KEND1
           CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     *                 KODPP1,KODPP2,KRDPP1,KRDPP2,
     *                 KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     *                 WORK(KEND1),LWRK1,IPRERI)
           KEND1 = KFREE
           LWRK1 = LFREE
         END IF
         NTOSYM = 1
        ELSE
         NTOSYM = NSYM
        ENDIF
!
        KENDSV = KEND1
        LWRKSV = LWRK1
!
        ICDEL1 = 0

        DO 100 ISYMD1 = 1,NTOSYM
!
         IF (DIRECT) THEN
            IF (HERDIR) THEN
              NTOT = MAXSHL
            ELSE
              NTOT = MXCALL
            END IF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
!
         DO 110 ILLL = 1,NTOT
!
!---------------------------------------------
!           If direct calculate the integrals.
!---------------------------------------------
!
            IF (DIRECT) THEN
!
               KEND1 = KENDSV
               LWRK1 = LWRKSV
!
               IF (HERDIR) THEN
                 CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                       IPRERI)
               ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     *                       WORK(KODCL1),WORK(KODCL2),
     *                       WORK(KODBC1),WORK(KODBC2),
     *                       WORK(KRDBC1),WORK(KRDBC2),
     *                       WORK(KODPP1),WORK(KODPP2),
     *                       WORK(KRDPP1),WORK(KRDPP2),
     *                       WORK(KCCFB1),WORK(KINDXB),
     *                       WORK(KEND1), LWRK1,IPRERI)
               END IF
!
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'CC_OMEGA2_RTCCD')
               END IF
!
            ELSE
               NUMDIS = 1
               KRECNR = KENDSV
            ENDIF
!
!-----------------------------------------------------
!           Loop over number of distributions in disk.
!-----------------------------------------------------
!
            DO 120 IDEL2 = 1,NUMDIS
!
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
!CN                  ISYMD = ISAO(IDEL)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
!                  WRITE(LUPRI,*)'DIRECT CASE ISYMD',ISYMD
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
!                  WRITE(LUPRI,*)'NON-DIRECT CASE ISYMD',ISYMD
               ENDIF
!
!----------------------------------------
!              Work space allocation two.
!----------------------------------------
!
               ISYDIS = MULD2H(ISYMD,ISYMOP)
!
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
!
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in CC_OMEGA2_RTCCD')
               ENDIF
!--------------------------------------------
!              Read AO integral distribution 
!              alpha>beta,gamma (^delta) da file.
!--------------------------------------------
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
!     
!-------------------------------------------
!                 Make oovv
!-------------------------------------------
!
!               IF(LOCDBG) THEN
!                  !Stupid and SLOW way, storage is OOVV
!                  CALL MAKE_OOVV(WORK(KIOOVV),WORK(KCMO),
!     *            WORK(KXINT),IDEL)
!               ELSE
                 !Faster way, Storage is VOOV
                  CALL RCCD_2O2V(WORK(KXINT),1, WORK(KIOOVV),
     *                       WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                       WORK(KLAMDP),1,WORK(KLAMDH),1,1,
     *                       WORK(KEND2),LWRK2,IDEL,ISYMD)
                                      
!               ENDIF
  120       CONTINUE
  110     CONTINUE
  100   CONTINUE

        if (IPRINT.GE.45) then
!          if (LOCDBG) then
!           write(lupri,*)'----- Integrals OOVV ----------------'
!           CALL OUTPUT(WORK(KIOOVV),1,NRHFT*NRHFT,1,NVIRT*NVIRT,
!     &                      NRHFT*NRHFT,NVIRT*NVIRT,1,LUPRI)
!          else
           write(lupri,*)'----- Integrals VOOV----------------'
           CALL OUTPUT(WORK(KIOOVV),1,NT1AMX,1,NT1AMX,
     &                      NT1AMX,NT1AMX,1,LUPRI)
!         ENDIF
        end if
!-------------------------
! Resort I_kiac to I_ai,ck
!-------------------------
!        if (locdbg) then
!          CALL Resort_I_kiac(WORK(KIaick),WORK(KIOOVV))
!          IF (IPRINT.GT.45) THEN
!             CALL AROUND('resort I_kiac -> I_aick ????')
!             CALL OUTPUT(WORK(KIaick),1,NT1AMX,1,NT1AMX,
!     &                      NT1AMX,NT1AMX,1,LUPRI)
!          ENDIF
!        ELSE  
!---------------------------------------
! Resort I_ikac (sorted aikc) to I_ai,ck
!---------------------------------------

          CALL Resort2_I_aikc(WORK(KIaick),WORK(KIOOVV))

          IF (IPRINT.GT.45) THEN
          CALL AROUND('resort2 I_ikac(aikc) -> I_aick')
          CALL OUTPUT(WORK(KIaick),1,NT1AMX,1,NT1AMX,
     &                      NT1AMX,NT1AMX,1,LUPRI)
          ENDIF
!        ENDIF  
      END IF

!------------------
! Omega_D like term
!------------------
      ! first linear term
      KT2SQ = KEND1            !squared amplitudes
      KLXSQ = KT2SQ + NT2SQ(1) !squared integrals
      KDint = KLXSQ + NT2SQ(1)
      KEND1 = KDint + NT2SQ(1)
      LWRK1 = LWORK-KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for D allocation '//
     &                'in CC_omega2_rccd')
      ENDIF
      CALL DZERO(WORK(KDint),NT2SQ(1))
      CALL DZERO(WORK(KT2SQ),NT2SQ(1))
      CALL DZERO(WORK(KLXSQ),NT2SQ(1))

      CALL CC_T2SQ(WORK(KT2AM),WORK(KT2SQ),1)
      
      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'Squared amplitudes'
        call CC_PRSQ(WORK(KEND1),work(KT2SQ),1,0,1)
      ENDIF
      !
      !CALL CC_T2SQ(WORK(KIAJB),WORK(KLXSQ),1)
      !
      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'Squared integrals'
        call CC_PRSQ(WORK(KEND1),work(KLXSQ),1,0,1)
      ENDIF
      !generate L_ai,kc = 2g_ai,ck-g_ac,ki (kiac resorted ai,ck)
      !generate L_ck,jb = 2g_ck,bj-g_cb,jk (jkcb resorted ck,bj)
      IF (RTCCD) THEN
         !ADD EXCHANGE PART TO INTEGRALS: I=g(C)-0.5g(ex)
C         WRITE(LUPRI,*)'ADDING X PART, SCALE BY HFXFAC'
         CALL DAXPY(NT2SQ(1),-HFXFAC*0.5d0,WORK(KIaick),1,WORK(KLXSQ),1)
      END IF
      !L=2*I
      call DSCAL(NT2SQ(1),2.0d0,WORK(KLXSQ),1)

      KOFF1 = IT2SQ(1,1) + KT2SQ
      KOFF2 = IT2SQ(1,1) + KLXSQ
      KOFF3 = IT2SQ(1,1) + KDint

      NTOTAI = MAX(NT1AM(1),1)
      NTOTCK = MAX(NT1AM(1),1)
      NTOTDL = MAX(NT1AM(1),1)

      CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,WORK(KOFF1),NTOTAI,WORK(KOFF2),NTOTCK,
     &            ONE,WORK(KOFF3),NTOTAI)

      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'RTCCD Linear term nr 1 '
        call CC_PRSQ(WORK(KEND1),work(KDInt),1,0,1)
      ENDIF

      KOFF1 = IT2SQ(1,1) + KT2SQ
      KOFF2 = IT2SQ(1,1) + KLXSQ
      KOFF3 = IT2SQ(1,1) + KDINT

      NTOTAI = MAX(NT1AM(1),1)
      NTOTCK = MAX(NT1AM(1),1)
      NTOTBJ = MAX(NT1AM(1),1)

      CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,WORK(KOFF2),NTOTAI,WORK(KOFF1),NTOTCK,
     &            ONE,WORK(KOFF3),NTOTAI)
C     
      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'After Addition of second Linear term'
        call CC_PRSQ(WORK(KEND1),work(KDint),1,0,1)
      ENDIF

      !the Quadratic term: square up L_iajb

      CALL CC_T2SQ(WORK(KLIAJB),WORK(KLXSQ),1)

      KOFF1 = IT2AM(1,1) + KT2SQ
      KOFF2 = IT2SQ(1,1) + KLXSQ
      KOFF3 = IT2SQ(1,1) + KOMEGAD
      KOFF4 = IT2SQ(1,1) + KDint

      NTOTAI = MAX(NT1AM(1),1)
      NTOTCK = MAX(NT1AM(1),1)
      NTOTDL = MAX(NT1AM(1),1)

      CALL DZERO(WORK(KOMEGAD),NT2SQ(1))

      CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,WORK(KOFF2),NTOTck,WORK(KOFF1),NTOTdl,
     &            ONE,WORK(KOFF3),NTOTdl)

!     WRITE(LUPRI,*) 'Second multiplication and add to linear terms'

      CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,WORK(KOFF1),NTOTai,WORK(KOFF3),NTOTck,
     &            ONE,WORK(KOFF4),NTOTai)

      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'RTCCD OMEGAD (squared)'
        call CC_PRSQ(WORK(KEND1),work(KDint),1,0,1)
      ENDIF

      !Pack contribution into Omega2 for solver
      IOPT=0
      CALL CC_T2PK(OMEGA2,WORK(KDint),1,IOPT)
C      
      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'OMEGA2 (after D) (packed)'
        call CC_PRP(WORK(KEND1),OMEGA2(1),1,0,1)
      ENDIF

!Add E term: it is taken twice!!!
      if (.not.(FOCKMAT)) then
         CALL DAXPY(NT2AMX,2.0d0,WORK(KOMEPK),1,OMEGA2,1)
C      
        IF(IPRINT.GT.45) THEN
          WRITE(LUPRI,*) 'OMEGA2 + E (packed)'
          call CC_PRP(WORK(KEND1),OMEGA2(1),1,0,1)
        ENDIF
        if ((nfield.gt.0).and.NONHF) then
           CALL DAXPY(NT2AMX,2.0d0,WORK(KOMEPK2),1,OMEGA2,1)
          !if (iprint.ge.45) then
             WRITE(LUPRI,*) 'OMEGA2 + E_field (packed)'
             call CC_PRP(WORK(KEND1),OMEGA2(1),1,0,1)
          !end if
        end if
       end if

!Add F term
      CALL DAXPY(NT2AMX,1.0d0,WORK(KOMEGAF),1,OMEGA2,1)

      IF(IPRINT.GT.45) THEN
        WRITE(LUPRI,*) 'OMEGA2 +F (packed)'
        call CC_PRP(WORK(KEND1),OMEGA2(1),1,0,1)
      ENDIF

!Scale the diagonal????
!      WRITE(LUPRI,*) 'Scaling the diagonal'
!      CALL CCLR_DIASCL(OMEGA2,HALF,1)

!
!If FOCKMAT compute E term via the HF Fock matrix + field
!
      if (fockmat) then
         kfckHF = kend1
         kcmof  = kfckHF + N2BST(1)
         kfield = kcmof  + NLAMDS
         kt2am  = kfield + N2BST(1)
         KT2SQ  = KT2AM + NT2AMX
         KEND1  = KT2SQ + NT2SQ(1)
         lwrk1  = lwork - kend1
         CALL DZERO(WORK(KCMOF),NLAMDS)
         CALL CMO_ALL(WORK(KCMOF),WORK(KEND1),LWRK1)
         call DZERO(WORK(KField),N2BST(1))
         IOPT = 2
         CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KEND1),WORK(KT2AM))
         CALL AROUND( 'In CC_RHS_RCCD: The T2 for E term other way' )
         call cc_prp(work(kend1),WORK(KT2AM),1,0,1) 
         CALL DZERO(WORK(KT2SQ),NT2SQ(1))
         CALL CC_T2SQ(WORK(KT2AM),WORK(KT2SQ),1)
!     This AO Fock matrix is constructed from the CMO transformed density
         LUFCK=-1
         CALL GPOPEN(LUFCK,'CC_FCKREF','UNKNOWN',' ','UNFORMATTED',
     *               IDUMMY,.FALSE.)
!
         REWIND(LUFCK)
         READ(LUFCK)(WORK(KFCKHF + I-1),I = 1,N2BST(1))
         CALL GPCLOSE(LUFCK,'KEEP' )
!
         IF (IPRINT .GT. 10) THEN
            CALL AROUND( 'HARTREE-FOCK Fock AO matrix' )
            CALL CC_PRFCKAO(WORK(KFCKHF),1)
         ENDIF
         if ((nfield.gt.0).and.(nonhf)) then
           CALL AROUND( 'CC_RHS_RCCD: ADD Finite FIELD TO FOCK!')
            DO IFI = 1, NFIELD
               FF =  EFIELD(IFI)
               CALL CC_ONEP(WORK(KFIELD),WORK(KEND1),LWRK1,FF,1,
     &              LFIELD(IFI))
            END DO
            CALL DAXPY(N2BST(1),1.0d0,WORK(KFIELD),1,WORK(KFCKHF),1)
            CALL CC_FCKMO(WORK(KFIELD),WORK(KCMOF),WORK(KCMOF),
     *                 WORK(KEND1),LWRK1,1,1,1)
            IF (IPRINT .GT. 10) THEN
              CALL AROUND( 'CC_RHS_RCCD: FF Fock MO matrix' )
              CALL CC_PRFCKMO(WORK(KFIELD),1)
            ENDIF
         end if
         !trasform to MO
         CALL CC_FCKMO(WORK(KFCKHF),WORK(KCMOF),WORK(KCMOF),
     *                 WORK(KEND1),LWRK1,1,1,1)
         IF (IPRINT .GT. 10) THEN
            CALL AROUND( 'In CC_RHS_RCCD: HF Fock MO matrix' )
            CALL CC_PRFCKMO(WORK(KFCKHF),1)
         ENDIF

         KEMAT1 = kend1
         KEMAT2 = KEMAT1 + NMATAB(1)
         KEpack = KEMAT2 + NMATIJ(1)
         kend1  = kEpack + nt2amx
         lwrk1  = lwork-kend1
         CALL DZERO(WORK(KEMAT1),NMATAB(1))
         CALL DZERO(WORK(KEMAT2),NMATIJ(1))
         ETRAN  = .FALSE.
         FCKCON = .TRUE.
         CALL CCRHS_EFCK(WORK(KEMAT1),WORK(KEMAT2),WORK(KCMOF),
     *                   WORK(KFCKHF),WORK(KEND1),LWRK1,FCKCON,
     *                   ETRAN,1)
        WRITE (LUPRI,*) 'EMAT2:'
        WRITE (LUPRI,'(5f12.6)') (WORK(KEMAT2-1+I),I=1,NMATIJ(1))
        WRITE (LUPRI,*) 'EMAT1:'
        WRITE (LUPRI,'(2f12.6)') (WORK(KEMAT1-1+I),I=1,NMATAB(1))

         CALL DZERO(WORK(KEpack),Nt2amx)
         CALL CCRHS_E(work(KEpack),WORK(KT2SQ),WORK(KEMAT1),
     &                 WORK(KEMAT2),WORK(KEND1),LWRK1,1,1)

C         CALL AROUND( 'In CC_RHS_RTCCD: The E term the other way' )
C         call cc_prp(work(kend1),WORK(KEpacK),1,0,1) 

         CALL DAXPY(NT2AMX,2.0d0,WORK(KEpacK),1,OMEGA2,1)
         !CALL DAXPY(NT2AMX,1.0d0,WORK(KEpacK),1,OMEGA2,1)
      end if
      CALL QEXIT('CC_OMEGA2_RTCCD')
      RETURN
      END
