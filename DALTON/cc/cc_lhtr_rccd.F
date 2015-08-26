C  /* Deck cc_lhtr_rccd */
      SUBROUTINE CC_LHTR_RCCD(ECURR,
     *                   FRHO1,LUFR1,FRHO2,LUFR2,
     *                   FC1AM,LUFC1,FC2AM,LUFC2,
     *                   RHO1,RHO2,CTR1,CTR2,WORK,LWORK,
     *                   NSIMTR,IVEC,ITR,LRHO1)
C
C     Written by Sonia and Fran, Nov-Dec 2009, based on CC_LHTR
C     Purpose:
C     Calculate left hand side transformation of a trial vector
C     in a noddy fashion for RCCD=RPA and DRCCD.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "ccorb.h"
#include "ccisao.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsections.h"
#include "ccfield.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccinftap.h"
#include "distcl.h"
#include "cbieri.h"
#include "cclr.h"
#include "eritap.h"
#include "ccslvinf.h"
#include "qm3.h"
#include "ccnoddy.h"
#include "r12int.h"
C
      CHARACTER*10 MODEL
      CHARACTER*8  FRHO2, FRHO1, FC2AM, FC1AM, FRHO12, FC12AM
      CHARACTER*1  LR
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOURTH = 0.25D0, THREE = 3.0D0)
      DIMENSION INDEXA(MXCORB_CC)
      DIMENSION RHO1(LRHO1,NSIMTR), CTR1(LRHO1,NSIMTR)
      DIMENSION RHO2(*), CTR2(*)
      DIMENSION WORK(LWORK)
C
      LOGICAL LOCDBG, FCKCON, ETRAN
      PARAMETER (LOCDBG = .false.)
      logical FOCKMAT
      parameter (FOCKMAT=.true.)

C
      CALL QENTER('CC_LHTR_RCCD')
C
C-----------------------------
C     Work-space allocation 1.
C-----------------------------
C
      ISYMOP = 1
      ISYRES = MULD2H(1,1)
      KT1AM  = 1
      KT2AM  = KT1AM + NT1AMX
      KT2SQ  = KT2AM + NT2AMX 
      KIAJB  = KT2SQ + NT2SQ(isymop) 
      KLIAJB = KIAJB + NT2AMX
      KI2SQ  = KLIAJB + NT2AMX
      KL2SQ  = KI2SQ + NT2SQ(isymop)
      KEND1  = KL2SQ + NT2SQ(isymop)
      LWRK1  = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation '//
     &                'in CC_LHTR_RCCD')
      ENDIF
C
C-------------------------------------------
C     Read zero'th order cluster amplitudes.
C-------------------------------------------
C
      IOPT = 2
      CALL CC_RDRSP('R0',0,1,IOPT,MODEL,WORK(KT1AM),WORK(KT2AM))
      CALL CC_T2SQ(WORK(KT2AM),WORK(KT2SQ),1)
      CALL DZERO(WORK(KT1AM),NT1AMX)
      if (iprint.ge.45) then
         write(lupri,*)' The squared T amplitudes'
         call cc_prsq(work(kend1),WORK(KT2SQ),1,0,1)
      end if
C
C-------------------------------------------------------------------
C     Read integrals ( ma | nb )=( am | bn ) from disc.
C     They are stored am,bn
C-------------------------------------------------------------------
C
      REWIND(LUIAJB)
      READ(LUIAJB) (WORK(KIAJB + J - 1), J = 1,NT2AM(ISYMOP))
      call DCOPY(NT2AMX,WORK(KIAJB),1,WORK(KLIAJB),1)
      if (RCCD) then
         IOPT = 1
         !generate L_iajb
         WRITE(LUPRI,*)'GENERATING 2C-E, MAKE 2C-FRAC*E'
         CALL CCSD_TCMEPK(WORK(KLIAJB),1.0d0,1,IOPT)
      else
         !generate 2*g_iajb
         call DSCAL(NT2AMX,TWO,WORK(KLIAJB),1)
      end if
      !square IAJB integrals for later use
      CALL CC_T2SQ(WORK(KIAJB),WORK(KI2SQ),1)
      if ((iprint.ge.45).or.LOCDBG) then
         write(lupri,*) 'CC_LHTR_RCCD: g_IAJB integrals (pck)'
         CALL CC_PRP(WORK(KEND1),WORK(KIAJB),1,0,1)
         call flshfo(lupri)
         write(lupri,*) 'CC_LHTR_RCCD: g_IAJB integrals (sqr)'
         call cc_prsq(work(kend1),WORK(KI2SQ),1,0,1)

         write(lupri,*) 'CC_LHTR_RCCD: L_IAJB integrals (pck)'
         CALL CC_PRP(WORK(KEND1),WORK(KLIAJB),1,0,1)
      end if
      !square L_iajb
      CALL CC_T2SQ(WORK(KLIAJB),WORK(KL2SQ),1)
      if (iprint.ge.45) then
         write(lupri,*) 'CC_LHTR_RCCD: L_IAJB integrals (sqr)'
         call cc_prsq(work(kend1),WORK(KL2SQ),1,0,1)
         call flshfo(lupri)
      end if

C-------------------------------------------
C     Intermediates for B-mat contributions
C-------------------------------------------

      KB_bj = KEND1
      KB_ai = KB_bj + NT2SQ(1)
      KEND1 = KB_ai + NT2SQ(1)
      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0 )
     &     CALL QUIT('INSUFFICIENT WORK SPACE IN CC_LHTR_RCCD')

      CALL DZERO(WORK(KB_bj),NT2SQ(1))
      CALL DZERO(WORK(KB_ai),NT2SQ(1))

      NTOTAI = MAX(NT1AM(1),1)
      NTOTBJ = MAX(NT1AM(1),1)
      NTOTCK = MAX(NT1AM(1),1)
      NTOTDL = MAX(NT1AM(1),1)

      !Primo intermedio B_bj,ck = sum_dl L_bj,dl*t_dl,ck =
      !                 = B_ck,bj = sum_dl t_ck,dl L_dl,bj
      CALL DGEMM('N','N',NTOTCK,NTOTBJ,NTOTDL,
     &            4.0d0,WORK(KT2SQ),NTOTCK,WORK(KL2SQ),NTOTDL,
     &            ONE,WORK(KB_bj),NTOTCK)
      if (iprint.ge.45) then
         write(lupri,*) 'CC_LHTR_RCCD: B_bj,ck=B_ck,bj intermediate'
         call cc_prsq(work(kend1),WORK(KB_bj),1,0,1)
      end if

!      !Secondo intermedio B_dl,ai=sum_ck tamp_dl,ck L_ck,ai
!      CALL DGEMM('N','N',NTOTDL,NTOTCK,NTOTAI,
!     &           4.0d0,WORK(KT2SQ),NTOTDL,WORK(KL2SQ),NTOTCK,
!     &            ONE,WORK(KB_AI),NTOTAI)
!      write(lupri,*) 'NODDY LHTR: B_dl,ai intermediate'
!      call cc_prsq(work(kend1),WORK(KB_ai),1,0,1)

      !Secondo intermedio B_ai,dl=sum_ck tamp_dl,ck L_ck,ai
      !                          =sum_ck L_ai,ck tamp_ck,dl
      CALL DGEMM('N','N',NTOTAI,NTOTDL,NTOTCK,
     &            4.0d0,WORK(KL2SQ),NTOTAI,WORK(KT2SQ),NTOTCK,
     &            ONE,WORK(KB_ai),NTOTAI)
      if (iprint.ge.45) then
         write(lupri,*) 'CC_LHTR_RCCD: B_ai,dl intermediate'
         call cc_prsq(work(kend1),WORK(KB_ai),1,0,1)
      end if

C----------------------------------------------------------
C     Initialize result-arrays and zero out single-vectors 
C----------------------------------------------------------
      NRHO2 = NT2AM(1)
      CALL DZERO(RHO1(1,1),NT1AM(ISYRES)*NSIMTR)
      CALL DZERO(RHO2,NRHO2)
      !set T1AM = 0, TBAR1=0
      CALL DZERO(WORK(KT1AM),NT1AMX)
      CALL DZERO(CTR1(1,1),NT1AM(ISYMTR)*NSIMTR)
C---------------------------------------
C Prepare for tbar E-term like contributions
C---------------------------------------
      KFOCKD = KEND1
      KEND1  = KFOCKD  + NORBTS
      LWRK1  = LWORK - KEND1
      IF (LWRK1 .LT. 0 )
     &     CALL QUIT('INSUFF. WORK SPACE IN CC_LHTR_RCCD (E prep)')
C-------------------------------------
C     Read canonical orbital energies.
C-------------------------------------
      LUSIFC = -1
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND (LUSIFC)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ (LUSIFC)
      READ (LUSIFC) (WORK(KFOCKD + I - 1), I = 1,NORBTS)
      CALL GPCLOSE(LUSIFC,'KEEP')
C
C----------------------------------------------------------------
C     Change symmetry ordering of the canonical orbital energies.
C----------------------------------------------------------------
C
      IF (FROIMP)
     *    CALL CCSD_DELFRO(WORK(KFOCKD),WORK(KEND1),LWRK1)
      CALL FOCK_REORDER(WORK(KFOCKD),WORK(KEND1),LWRK1)

      if ((iprint.ge.45).or.locdbg) then
        do i=1,NORBTS
         write(lupri,*) 'Epsilon_',i,' = ', WORK(KFOCKD+i-1)
        end do
      end if
C
C--------------------------------------------------------------------
C     Recover the full MO coefficient matrix for I_oovv 
C     for alternative E term. Otherwise only used for RCCD
C--------------------------------------------------------------------
C
      IF (RCCD) THEN
        KCMO  = KEND1
        KEND1 = KCMO  + NLAMDS
        LWRK1 = LWORK - KEND1
        IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation '//
     &                'in CC_LHTR_RCCD')
        ENDIF
        CALL DZERO(WORK(KCMO),NLAMDS)
        CALL CMO_ALL(WORK(KCMO),WORK(KEND1),LWRK1)
      END IF
C
C--------------------------------------------------------------------
C     If RCCD:
C     Allocate for OOVV integrals
C--------------------------------------------------------------------
C
      IF (RCCD) THEN
        KLAMDP = KEND1
        KLAMDH = KLAMDP + NLAMDT
        KEND1  = KLAMDH + NLAMDT
        LWRK1  = LWORK  - KEND1
        IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insufficient memory for allocation '//
     &                'in CC_LHTR_RCCD')
        ENDIF
C
        CALL DZERO(WORK(KLAMDP),NLAMDT)
        CALL DZERO(WORK(KLAMDH),NLAMDT)
        CALL LAMMAT(WORK(KLAMDP),WORK(KLAMDH),WORK(KT1AM),
     &       WORK(KEND1),LWRK1)

        KIOOVV = KEND1
        KIaick = KIOOVV + NRHFT*NRHFT*NVIRT*NVIRT
        KEND1  = KIaick + NT2SQ(1)
        LWRK1  = LWORK - KEND1
        IF (LWRK1 .LT. 0) THEN
           WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
           CALL QUIT('Insufficient memory for allocation '//
     &                'in CC_LHTR_RCCD')
        ENDIF
        call dzero(WORK(KIOOVV),NRHFT*NRHFT*NVIRT*NVIRT)
        call dzero(WORK(KIaick),NT2SQ(1))
C------------------------------------------------------------------
C Start loop on integral distributions to build (OO|VV)
C------------------------------------------------------------------
        KENDS2 = KEND1
        LWRKS2 = LWRK1
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
        KENDSV = KEND1
        LWRKSV = LWRK1
C
        ICDEL1 = 0
        DO ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
              NTOT = MAXSHL
            ELSE
              NTOT = MXCALL
            END IF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
         DO ILLL = 1,NTOT
C---------------------------------------------
C           If direct calculate the integrals.
C---------------------------------------------
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               IF (HERDIR) THEN
                 CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                       IPRERI)
               ELSE
                 CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                       WORK(KODCL1),WORK(KODCL2),
     &                       WORK(KODBC1),WORK(KODBC2),
     &                       WORK(KRDBC1),WORK(KRDBC2),
     &                       WORK(KODPP1),WORK(KODPP2),
     &                       WORK(KRDPP1),WORK(KRDPP2),
     &                       WORK(KCCFB1),WORK(KINDXB),
     &                       WORK(KEND1), LWRK1,IPRERI)
               END IF
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'CC_LHTR_RCCD')
               END IF
            ELSE
               NUMDIS = 1
               KRECNR = KENDSV
            ENDIF
C-----------------------------------------------------
C           Loop over number of distributions in disk.
C-----------------------------------------------------
            DO IDEL2 = 1,NUMDIS
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  IF (NOAUXB) THEN
                     IDUM = 1
                     CALL IJKAUX(IDEL,IDUM,IDUM,IDUM)
                  END IF
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C----------------------------------------
C              Work space allocation two.
C----------------------------------------
               ISYDIS = MULD2H(ISYMD,ISYMOP)
               KXINT  = KEND1
               KEND2  = KXINT + NDISAO(ISYDIS)
               LWRK2  = LWORK - KEND2
               IF (LWRK2 .LT. 0) THEN
                  WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
                  CALL QUIT('Insufficient memory for integrals '//
     &                      'in CC_LHTR_RCCD')
               ENDIF
C------------------------------------------------
C              Read AO integral distribution 
C              alpha>beta,gamma (^delta) da file.
C------------------------------------------------
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C-------------------------------------------
C              Make oovv, they come out sorted v1o1o2,v2
C-------------------------------------------
               CALL RCCD_2O2V(WORK(KXINT),1, WORK(KIOOVV),
     *                        WORK(KLAMDP),1,WORK(KLAMDH),1,
     *                        WORK(KLAMDP),1,WORK(KLAMDH),1,1,
     *                        WORK(KEND2),LWRK2,IDEL,ISYMD)

            END DO
         END DO
        END DO
C
C Done computing additional 2-electron integrals
C
      if ((iprint.ge.45).or.LOCDBG) then
        write(lupri,*)'----- Integrals OOVV ----------------'
        CALL OUTPUT(WORK(KIOOVV),1,NCKI(1),1,NVIR(1),
     &                    NCKI(1),NVIR(1),1,LUPRI)
      end if
C-----------------------------------------------------------------
C Resort I_kiac=I_ikac(aikc) ->  I_ai,ck 
C-----------------------------------------------------------------
        CALL Resort2_I_aikc(WORK(KIaick),WORK(KIOOVV))
        if (iprint.ge.45) then
           CALL AROUND('CC_LHTR_RCCD: resort I_kiac -> I_aick')
           CALL OUTPUT(WORK(KIaick),1,NT1AMX,1,NT1AMX,
     &                    NT1AMX,NT1AMX,1,LUPRI)
        end if

      ELSE

        KENDS2 = KEND1
        LWRKS2 = LWRK1

      END IF !RCCD
      !----------------------------------------
      !Prepare integrals for A-mat contribution
      !----------------------------------------
      
      if (iprint.ge.45) then
        WRITE(LUPRI,*) 'Squared IAJB integrals once more'
        call CC_PRSQ(WORK(KEND1),work(KI2SQ),1,0,1)
      end if
      !generate L_ai,kc = 2g_ai,ck-g_ac,ki (kiac resorted ai,ck)
      !generate L_ck,jb = 2g_ck,bj-g_cb,jk (jkcb resorted ck,bj)
      IF (RCCD) THEN
         !ADD EXCHANGE CONTRIBUTION
         WRITE(LUPRI,*)'TAKE CARE! ADD HFXFAC*EXCONT'
         CALL DAXPY(NT2SQ(1),-0.5d0,WORK(KIaick),1,WORK(KI2SQ),1)
         call DSCAL(NT2SQ(1),2.0d0,WORK(KI2SQ),1)
         if (iprint.ge.45) then
            WRITE(LUPRI,*) 'L_aijb integrals'
            call CC_PRSQ(WORK(KEND1),work(KI2SQ),1,0,1)
         end if
      else
         call DSCAL(NT2SQ(1),2.0d0,WORK(KI2SQ),1)
         if (iprint.ge.45) then
            WRITE(LUPRI,*) 'DRCCD: 2*g_aijb integrals'
            call CC_PRSQ(WORK(KEND1),work(KI2SQ),1,0,1)
         end if
      END IF
C
C-----------------------------
C     Loop over trial vectors.
C-----------------------------
C
      DO 125 IV = 1,NSIMTR
C
         IF (.NOT. MINSCR) THEN
C
C-------------------------------------------------------
C           Read CTR2 from disc into RHO2 and square up into CTR2.
C-------------------------------------------------------
C
            CALL CC_RVEC(LUFC2,FC2AM,NT2AM(ISYMTR),NT2AM(ISYMTR),
     *                   IVEC+IV-1,RHO2)
C
            CALL CC_T2SQ(RHO2,CTR2,ISYMTR)
C
            IF ( LOCDBG ) THEN
               RHO1N = DDOT(NT1AM(ISYMTR),CTR1(1,IV),1,CTR1(1,IV),1)
               RHO2N = DDOT(NT2SQ(ISYMTR),CTR2,1,CTR2,1)
               WRITE(LUPRI,1) 'Norm of CTR1 -Read in after loop:  ',
     &              RHO1N
               WRITE(LUPRI,1) 'Norm of CTR2 -Read in after loop:  ',
     &              RHO2N
            ENDIF
C
C
C----------------------------------------
C           Read result vector from disc.
C----------------------------------------
C
            CALL CC_RVEC(LUFR2,FRHO2,NRHO2,NRHO2,IV+ITR-1,RHO2)

         ENDIF
C
         IF ( LOCDBG ) THEN
            RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
            RHO2N = DDOT(NRHO2,RHO2,1,RHO2,1)
            WRITE(LUPRI,1) 'Norm of RHO1 loop over vect. 1.   ', RHO1N
            WRITE(LUPRI,1) 'Norm of RHO2 loop over vect. 1.   ', RHO2N
         ENDIF
C
         !KRHO2 = KENDS2
         KRHO2 = KEND1
         KEND1 = KRHO2 + NT2SQ(1)
         LWRK1 = LWORK - KEND1
         CALL DZERO(WORK(KRHO2),NT2SQ(1))
C
C-------------------------------------------------------------------
C  Compute for tbar E-term like contributions 
C-------------------------------------------------------------------
C
         CALL RCCD_E(work(KFockD),CTR2,WORK(KRHO2),
     &                WORK(KEND1),LWRK1,CTR2,WORK(KEND1))
         CALL DSCAL(NT2SQ(1),2.0d0,WORK(KRHO2),1)

         if (iprint.ge.45) then
            WRITE(LUPRI,*) 'The E-term like contribution (2x)'
            call cc_prsq(work(kend1),WORK(KRHO2),1,0,1)
            call flshfo(lupri)
         end if
C
C----------------------------------------------------------
C  Contract CTR2 with B-Mat intermediates and add to result 
C  CTR2 is SQUARED!!!
C  names of lengths are wrong, but all equally long!
C----------------------------------------------------------
C
         NTOTAI = MAX(NT1AM(1),1)
         NTOTCK = MAX(NT1AM(1),1)
         NTOTDL = MAX(NT1AM(1),1)
         NTOTBJ = MAX(NT1AM(1),1)

         CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            1.0d0,CTR2(1),NTOTAI,WORK(KB_bj),NTOTCK,
     &            ONE,WORK(KRHO2),NTOTAI)

         CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            1.0d0,WORK(KB_ai),NTOTAI,CTR2(1),NTOTDL,
     &            ONE,WORK(KRHO2),NTOTAI)

!-------------------------------------------
! COMPUTE A-mat CONTRIBUTIONS tbar*L_aijb
! L_aijb (or 2*G_iajb for DRCCD) is in KI2SQ 
!-------------------------------------------

         NTOTAI = MAX(NT1AM(1),1)
         NTOTCK = MAX(NT1AM(1),1)
         NTOTBJ = MAX(NT1AM(1),1)

         CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,CTR2(1),NTOTAI,WORK(KI2SQ),NTOTBJ,
     &            ONE,WORK(KRHO2),NTOTAI)

         if (iprint.ge.45) then
            WRITE(LUPRI,*) 
     &   'CC_LHTR_RCCD: After Addition of first Linear term'
            call CC_PRSQ(WORK(KEND1),work(KRHO2),1,0,1)
         end if
         NTOTAI = MAX(NT1AM(1),1)
         NTOTCK = MAX(NT1AM(1),1)
         NTOTBJ = MAX(NT1AM(1),1)

         CALL DGEMM('N','N',NT1AM(1),NT1AM(1),NT1AM(1),
     &            2.0d0,WORK(KI2SQ),NTOTAI,CTR2,NTOTCK,
     &            ONE,WORK(KRHO2),NTOTAI)
         IF (iprint.ge.45) THEN
           WRITE(LUPRI,*) 
     &     'CC_LHTR_RCCD: After Addition of second Linear term'
           CALL CC_PRSQ(WORK(KEND1),work(KRHO2),1,0,1)
         END IF
C
C--------------------------------------------------------------------
C     Print out result vectors - zero out single vectors 
C--------------------------------------------------------------------
C
         !Pack final vector into RHO2
         IOPT=0
         !SONIA e FRAN: SCALE THE RHO VECTOR BY 2????
         !CALL DSCAL(NT2SQ(1),0.5d0,WORK(KRHO2),1)
         CALL CC_T2PK(RHO2,WORK(KRHO2),1,IOPT)
         CALL DZERO(RHO1(1,IV),NT1AM(ISYRES))
C
         IF ((IPRINT .GT. 50).or.LOCDBG) THEN
            CALL AROUND('Transformed vectors out of CC_LHTR_RCCD')
            CALL CC_PRP(RHO1(1,IV),RHO2,ISYRES,1,1)
         ENDIF
C
         IF ( LOCDBG ) THEN
           RHO1N = DDOT(NT1AM(ISYRES),RHO1(1,IV),1,RHO1(1,IV),1)
           RHO2N = DDOT(NT2AM(ISYRES),RHO2,1,RHO2,1)
           WRITE(LUPRI,1) 'Norm of RHO1 coming out CC_LHTR_RCCD:', RHO1N
           WRITE(LUPRI,1) 'Norm of RHO2 coming out CC_LHTR_RCCD:', RHO2N
         ENDIF

C----------------------------------
C     Write out transformed vector.
C     Write out in all cases.
C----------------------------------

         CALL CC_WVEC(LUFR1,FRHO1,NT1AM(ISYMTR),NT1AM(ISYMTR),
     *                IV + IVEC -1,RHO1(1,IV))
C
         CALL CC_WVEC(LUFR2,FRHO2,NRHO2,NT2AM(ISYMTR),
     *                   IV + ITR -1,RHO2)
C
  125 CONTINUE
C
C-------------------------------
C     Write out program timings.
C-------------------------------
C
   1  FORMAT(1x,A35,1X,E20.10)
C
      CALL QEXIT('CC_LHTR_RCCD')
      RETURN
      END
