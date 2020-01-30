!--------------------------------------------------------------------
      SUBROUTINE CC_LANCZOS_LR(LABEL,J,ER,EI,R,L,F,U,V)
!--------------------------------------------------------------------
!     Given a set of Lanzcos T eigenvalues and eigenvector this
!     calculates the linear response function for 
!     given frequencies and gammas
!     J: Chain Length 
!     ER: Eigenvalues, real part 
!     EI: Eigenvalues, imag part 
!     R: Eigenvectors right (stored as cols)
!     L: Eigenvectors Left  (stored as cols, e.i. eigvectors of T^T)
!     F: F-matrix "trans" (ikke double q-transformed) 
!     U, V are normalization factors
!
!     FIXME: todo list:
!     1) rethink handling of complex roots
!     2) make loops more efficient (vectors maybe?)
!     3) "weight" thing
!
!     Sonia Coriani, 2010-2012
!--------------------------------------------------------------------
      IMPLICIT NONE
#include "priunit.h"
#include "cclrlancz.h"
#include "codata.h"
     
      INTEGER J,IFREQ,NFREQ
#if defined (SYS_CRAY) 
      REAL U,V
      REAL FSTART,FSTOP,FSTEP
      REAL ER(J),EI(J),R(J,J),L(J,J),F(J,J)
      REAL ZERO, ONE, TWO, PDIFFI,PSUMI,PSUMK,XGAMMA
      REAL LR_OM_RE,LR_OM_RE_T,LR_OM_RE_F
      REAL LR_OM_IM,LR_OM_IM_T,LR_OM_IM_F
      REAL DPI, DPK,FACTOR,FACTOR2, FREQ, DMI
#else 
      DOUBLE PRECISION U,V
      DOUBLE PRECISION FSTART,FSTOP,FSTEP
      DOUBLE PRECISION ER(J),EI(J),R(J,J),L(J,J),F(J,J)
      DOUBLE PRECISION ZERO, ONE, TWO, PDIFFI,PSUMI,PSUMK,XGAMMA
      DOUBLE PRECISION LR_OM_RE,LR_OM_RE_T,LR_OM_RE_F
      DOUBLE PRECISION LR_OM_IM,LR_OM_IM_T,LR_OM_IM_F
      DOUBLE PRECISION DPI, DPK,FACTOR,FACTOR2, FREQ, DMI
      DOUBLE PRECISION XGAM_PLUS,XGAM_MINUS,FACTORR,FACTORI
      Double precision add_to_real, add_to_im
      Double precision tot_add_real, tot_add_im
      double precision weight(j,2)
#endif 
      LOGICAL NOIMAG
      INTEGER I, K, IDAMP, IDUMMY
!
      PARAMETER (ZERO = 0.0D0, ONE  = 1.0D0, TWO = 2.0D0)
      !Specifics of the file containing the Lanczos abs spectrum 
      CHARACTER*7 FNSPECTRUM
      PARAMETER (FNSPECTRUM='CCSPECX')
      INTEGER LUSPEC
      !
      CHARACTER*(8) LABEL
      LOGICAL LOCDBG
      PARAMETER (LOCDBG=.false.)
!
      LOGICAL SKIPNEXT 
      SKIPNEXT=.FALSE.
!
! Set these using input data
!
      FSTART = FREQ_RANGE(1)
      FSTOP  = FREQ_RANGE(2)
      FSTEP  = FREQ_RANGE(3)
      XGAMMA = DAMPING(1)
      LR_OM_RE_T=zero
      LR_OM_IM_T=zero
      !CALL DZERO(weight,j*2)
      
      NFREQ=1+INT((FSTOP-FSTART)/FSTEP)
      WRITE (LUPRI,*) "1. FREQ", FSTART," NFREQ:",NFREQ, " STEP:", FSTEP
!
! Open unit with results
!      
      LUSPEC=-1
      CALL GPOPEN(LUSPEC,FNSPECTRUM,'NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

      !initialize to zero
      add_to_real = zero
      add_to_im   = zero
      tot_add_real = zero
      tot_add_im   = zero
C
C     Loop over chain 
C
      DO IDAMP=1,NDAMP
         XGAMMA = DAMPING(IDAMP)

         WRITE (LUPRI,*) " FOR OPERATOR LABEL = ", LABEL
         WRITE (LUPRI,'(1X,A27,F12.8,A3,I10)') 
     *               " Response fun with gamma = ",XGAMMA, ' J=',J
         WRITE(LUPRI,*)  " Frequency (eV)  & Real  (au)& Imaginary(au)",
     &                   "&  sigma^I (au)"
         WRITE(LUPRI,'(X,"---------------------------")')
!
!Also on file LUSPEC
!
         WRITE (LUSPEC,*) " FOR OPERATOR LABEL = ", LABEL
         WRITE (LUSPEC,'(1X,A27,F12.8,A3,I10)')
     *               " Response fun with gamma = ",XGAMMA, ' J=',J
         WRITE(LUSPEC,*) " Frequency (eV)  & Real  (au)& Imaginary(au)",
     &                   "&  sigma^I (au)"
         WRITE(LUSPEC,'(X,"---------------------------")')

         DO IFREQ=1,NFREQ

             FREQ = FSTART+FSTEP*(IFREQ-1)
             LR_OM_RE=ZERO
             LR_OM_IM=ZERO
             LR_OM_RE_F=ZERO
             LR_OM_IM_F=ZERO
!
!reinitialize inside freq loop
!
             add_to_real = zero
             add_to_im   = zero
             tot_add_real = zero
             tot_add_im   = zero

             DO I=1,J
               IF (SKIPNEXT) THEN
                  SKIPNEXT = .FALSE.
               ELSE
                !check on imaginary eigenvalue = 0
                !IF (ABS(EI(I)).EQ.ZERO)  THEN
                IF (EI(I).EQ.ZERO)  THEN
                   FACTOR=U*V*R(1,I)*L(1,I)
C
C                  No imaginary eigenvector, calc eta*t conts 
C
                   PDIFFI = FREQ-ER(I)
                   PSUMI  = FREQ+ER(I)
                   !denominators
                   DMI = PDIFFI*PDIFFI + XGAMMA*XGAMMA
                   DPI = PSUMI*PSUMI + XGAMMA*XGAMMA
                   LR_OM_RE = LR_OM_RE + FACTOR*(PDIFFI/DMI-PSUMI/DPI)
                   LR_OM_IM = LR_OM_IM - FACTOR*XGAMMA*(ONE/DMI-ONE/DPI)
                   !weight(i,1) = ER(I)
                   !weight(i,2) = R(1,I)*L(1,I)
C
C                  Include F-part cont 
C
                   FACTOR = V*V*L(1,I)
                   DO K=1,J
                      FACTOR2 = FACTOR*L(1,K)*F(I,K)
                      if (abs(ei(k)).eq.zero) then
                         PSUMK  = FREQ+ER(K)
                         !denominator
                         DPK = PSUMK*PSUMK + XGAMMA*XGAMMA
                         LR_OM_RE_F =  LR_OM_RE_F + FACTOR2*
     &                   (PSUMK*PDIFFI-XGAMMA*XGAMMA)/(DMI*DPK)
                         LR_OM_IM_F =  LR_OM_IM_F - FACTOR2*
     &                   XGAMMA*(PSUMK+PDIFFI)/(DMI*DPK)
                      else
                         !write(lupri,*)'LANCZOS, do not know, index=',k
                         !SHOULD DECIDE HOW TO HANDLE THIS
!                        FIXME 1)
                      end if
                   ENDDO

                ELSE

                   FACTORR=U*V*(R(1,I)*L(1,I)-R(1,I+1)*L(1,I+1))
                   FACTORI=U*V*(R(1,I+1)*L(1,I)+R(1,I)*L(1,I+1))
!
!               Include also states with imaginary excitation energy 
!
                   PDIFFI = FREQ-ER(I)
                   PSUMI  = FREQ+ER(I)
                   xgam_minus = xgamma - EI(I)
                   xgam_plus  = xgamma + EI(I)
                   !denominators
                   DMI = PDIFFI*PDIFFI + XGAM_minus*XGAM_minus
                   DPI = PSUMI*PSUMI+ XGAM_plus*XGAM_plus
                   !
                   !real response function (eta part)
                   !
                   add_to_real = FACTORR*(PDIFFI/DMI-PSUMI/DPI)
     &                        + FACTORI*(xgam_minus/DMI-xgam_plus/DPI)
                   tot_add_real = tot_add_real + add_to_real
                   LR_OM_RE = LR_OM_RE + FACTORR*(PDIFFI/DMI-PSUMI/DPI)
     &                        + FACTORI*(xgam_minus/DMI-xgam_plus/DPI)
                   !
                   !imaginary response function (eta part)
                   !
                   add_to_im = FACTORR*(XGAM_minus/DMI-
     &                        xgam_plus/DPI) + FACTORI*(PDIFFI/DMI
     &                        - PSUMI/DPI)
                   tot_add_im = tot_add_im + add_to_im
                   LR_OM_IM = LR_OM_IM - FACTORR*(XGAM_minus/DMI-
     &                        xgam_plus/DPI) + FACTORI*(PDIFFI/DMI
     &                        - PSUMI/DPI)
C
C                  Include F-part cont 
C
!                   FIXME 1)
!                   write(lupri,*)'WARNING: CMPLX ROOT, SKIP F PART'
!                   FACTOR = V*V*L(1,I)
!                   DO K=1,J
!                      FACTOR2 = FACTOR*L(1,K)*F(I,K)
!                      PSUMK  = FREQ+ER(K)
!                      !denominator
!                      DPK = PSUMK*PSUMK + XGAMMA*XGAMMA
!                      LR_OM_RE_F =  LR_OM_RE_F + FACTOR2*
!     &                (PSUMK*PDIFFI-XGAMMA*XGAMMA)/(DMI*DPK)
!                      LR_OM_IM_F =  LR_OM_IM_F - FACTOR2*
!     &                XGAMMA*(PSUMK+PDIFFI)/(DMI*DPK)
!                   ENDDO
                    SKIPNEXT=.true.
                ENDIF
              end if !skipnext
C
            ENDDO
           
            if (locdbg) then
              WRITE (LUPRI,*)'tot_add_real', tot_add_real
              WRITE (LUPRI,*)'tot_add_im', tot_add_im
              WRITE (LUPRI,*)'Freq; F contrib Re', FREQ, LR_OM_RE_F
              WRITE (LUPRI,*)'Freq; F contrib Im', FREQ, LR_OM_IM_F
            end if

            LR_OM_RE_T = LR_OM_RE - LR_OM_RE_F
            LR_OM_IM_T = LR_OM_IM - LR_OM_IM_F

       WRITE (LUSPEC,'(6F16.8)') 
     &       FREQ*XTEV, LR_OM_RE_T, LR_OM_IM_T, -LR_OM_IM_T*FREQ, 
     &       LR_OM_RE,LR_OM_IM 
 
         IF (IFREQ.LE.10) then
            !dump first 10 values on output file as well
       WRITE (LUPRI,'(6F16.8)')
     &       FREQ*XTEV, LR_OM_RE_T, LR_OM_IM_T, -LR_OM_IM_T*FREQ,
     &       LR_OM_RE,LR_OM_IM
         end if

         ENDDO
         WRITE(LUPRI,'(X,"--- remaining values on LUSPEC file ---")')
         WRITE(LUSPEC,'(X,"---------------------------")')
      ENDDO

      CALL GPCLOSE(LUSPEC,'KEEP')

      RETURN
      END
!------------------------------------------------------------
      SUBROUTINE CC_BIORTOCOMP(J,ER,EI,EVR,EVL,WORK,LWORK) 
!------------------------------------------------------------
!     From a DGEEV output, take eigenvectors  with unit norm and 
!     make them biorthogonal so that <R|R>=1,<L|R>=1, by 
!     scaling the left eigenvectors with 1/<L|R>.
!     NB: Take special care with complex that are assumed
!     to come with eigenvectors order as 
!     Right eigenvector for w_k_r + iw_k_i is R(k) + i R(k+1)
!     Right eigenvector for w_k_r - iw_k-I is R(k) - i R(k+1)
!     and similar for left eigenvectors 
!     Ove Christiansen & Sonia Coriani, November 2010
!
!     J = chain length, equal to length of each eigenvector
!         and to nr of eigenvalues
!     ER = real parts of eigenvalues
!     EI = imaginary parts of eigenvalues
!     EVR = Right eigenvectors (DGEEV solutions)
!     EVL = left eigenvectors (DGEEV solutions)
C--------------------------------------------------------------------
      IMPLICIT NONE
#include "priunit.h"
      INTEGER J,LWORK
#if defined (SYS_CRAY) 
      REAL ER(J),EI(J),EVR(J,J),EVL(J,J),WORK(LWORK)
      REAL XNORM,THRESH
      REAL ONE,ZERO,DDOT
#else 
      DOUBLE PRECISION ER(J),EI(J),EVR(J,J),EVL(J,J),WORK(LWORK)
      DOUBLE PRECISION XNORM,THRESH
      DOUBLE PRECISION ONE,ZERO,DDOT
#endif 
      PARAMETER (ONE=1.0d0,ZERO=0.0d0,THRESH=1.0d-16)
      INTEGER K, IDAMP, KLRO, KLIO
      LOGICAL SKIPNEXT 
      
      SKIPNEXT=.FALSE.
      !write(lupri,*)'I AM INSIDE BIORTO COMP'

      KLRO = 1        !L real Output
      KLIO = KLRO + J !L imaginary Output
      IF (KLIO+J .GT. LWORK)
     &   CALL QUIT("Too little space in CC_BIORTOCOMP")

      DO K=1,J
         IF (SKIPNEXT) THEN
           SKIPNEXT = .FALSE. 
         ELSE 
            IF (EI(K).EQ.ZERO) THEN 
!===================================
!              Real eigenvalue 
!===================================
               XNORM = DDOT(J,EVL(1,K),1,EVR(1,K),1) 
               !write(lupri,*)'Biortho of RL', xnorm
               IF (ABS(XNORM).LE. THRESH) THEN 
                 CALL QUIT("Error in CC_BIORTOCOMP") 
               ENDIF 
               XNORM = ONE/XNORM 
               !write(lupri,*)'Biortho of RL', xnorm
               CALL DSCAL(J,XNORM,EVL(1,K),1) 
               !XNORM = DDOT(J,EVL(1,K),1,EVR(1,K),1) 
               !write(lupri,*)'Biortho of RL dopo',xnorm,' K=', K, 'J=',J
            ELSE 
!=============================
!              complex pair 
!=============================
               write(lupri,*)'WARNING: COMPLEX PAIR alternative norm'
               call flshfo(lupri)
               CALL CC_BIORTOCOMP2(J,EVR(1,K),EVR(1,K+1),
     *                             EVL(1,K),EVL(1,K+1),
     *                             WORK(KLRO),WORK(KLIO)) 
      !SUBROUTINE CC_BIORTOCOMP2(N,EVRR,EVRI,EVLR,EVLI,EVLRO,EVLIO)
               call flshfo(lupri)
               CALL DCOPY(J,WORK(KLRO),1,EVL(1,K),1)
               CALL DCOPY(J,WORK(KLIO),1,EVL(1,K+1),1)

               SKIPNEXT =.TRUE. 
            ENDIF
         ENDIF
      ENDDO 

      RETURN 
      END 
!--------
      SUBROUTINE CC_BIORTOCOMP1(N,EVRR,EVRI,EVLR,EVLI,EVLRO,EVLIO)
C--------------------------------------------------------------------
C     OC,03/11/2010
C     Take eigenvectors in with unit norm. 
C     Create new left vectors eigenvectors such that they are biorthog
C     to right eigenvectors. 
C
C     EVRR: Right eigenvector, Real part. 
C     EVRI: Right eigenvector, Imaginary part. 
C     EVLR: Eigenvector, left,  real part. 
C     EVLI: Eigenvector, left,  imag. part. 
C     EVLRO: Eigenvectors, left,  real part. OUTPUT
C     EVLIO: Eigenvectors, left,  imag. part. OUTPUT
C     N : length of vector
C     (we pass one eigenvector at a time)
C--------------------------------------------------------------------
C
      IMPLICIT NONE
#include "priunit.h"
      INTEGER N 
#if defined (SYS_CRAY) 
      REAL EVRR(N),EVRI(N),EVLR(N),EVLI(N),EVLRO(N),EVLIO(N) 
      REAL CR,CI,CNSQ,DDOT
#else 
      !input
      DOUBLE PRECISION EVRR(N),EVRI(N),EVLR(N),EVLI(N)
      !output
      Double precision EVLRO(N),EVLIO(N)
      !local
      DOUBLE PRECISION CR,CI,CNSQ,DDOT
#endif 

      write(lupri,*)'I AM INSIDE BIORTCOMP1'
      call flshfo(lupri)
      write(lupri,*)'<rr|rr>=',ddot(N,EVRR,1,EVRR,1)
      write(lupri,*)'<ri|ri>=',ddot(N,EVRI,1,EVRI,1)
      write(lupri,*)'<lr|lr>=',ddot(N,EVLR,1,EVLR,1)
      write(lupri,*)'<li|li>=',ddot(N,EVLI,1,EVLI,1)
C
      write(lupri,*)'<lr|Rr>=',ddot(N,EVLR,1,EVRR,1)
      write(lupri,*)'<li|Ri>=',ddot(N,EVLI,1,EVRI,1)
      write(lupri,*)'<lr|Ri>=',ddot(N,EVLR,1,EVRI,1)
      write(lupri,*)'<li|Rr>=',ddot(N,EVLI,1,EVRR,1)
      !<ReR|ReL>-<ImR|ImL>
      CR = DDOT(N,EVLR,1,EVRR,1)-DDOT(N,EVLI,1,EVRI,1)
      write(lupri,*)'CR=', CR
      !<ReR|ImL>+<ImR|ReL>
      CI = DDOT(N,EVLI,1,EVRR,1)+DDOT(N,EVLR,1,EVRI,1)
      write(lupri,*)'CI=', CI
      CNSQ=CR*CR+CI*CI
      write(lupri,*)'CNSQ', CNSQ
      write(lupri,*)'CI/CNSQ', CI/CNSQ, 'CR/CNSQ', Cr/Cnsq
!---
!      !<ReL|ReI>+<ImL|ImR>
!      CRp = DDOT(N,EVLR,1,EVRR,1)+DDOT(N,EVLI,1,EVRI,1)
!      !<ReL|ImR>-<ImL|ReR>
!      CIm = DDOT(N,EVLR,1,EVRI,1)-DDOT(N,EVLI,1,EVRR,1)
!      CNmSQ=CRp*CRp+CIm*CIm
C
      CALL DZERO(EVLRO,N)
      CALL DZERO(EVLIO,N)
C
      !write(lupri,*)'CR/CNSQ=',CNSQ
      !call flshfo(lupri)
      CALL DCOPY(N,EVLR,1,EVLRO,1)
      CALL DCOPY(N,EVLI,1,EVLIO,1)
      CALL DSCAL(N,(CR/CNSQ),EVLRO,1) 
      CALL DSCAL(N,(CR/CNSQ),EVLIO,1) 
C
      CALL DAXPY(N,(CI/CNSQ),EVLI,1,EVLRO,1)
      CALL DAXPY(N,(-CI/CNSQ),EVLR,1,EVLIO,1)
!
      WRITE(LUPRI,*)"This should be one ", 
     *      DDOT(N,EVLRO,1,EVRR,1)-DDOT(N,EVLIO,1,EVRI,1) 
      call flshfo(lupri)
      WRITE(LUPRI,*)"This should be zero ", 
     *      DDOT(N,EVLRO,1,EVRI,1)+DDOT(N,EVLIO,1,EVRR,1) 
      call flshfo(lupri)

      WRITE(LUPRI,*)"<new LrO|Rr>", DDOT(N,EVLRO,1,EVRR,1)
      WRITE(LUPRI,*)"<new LiO|Ri>", DDOT(N,EVLiO,1,EVRi,1)
      WRITE(LUPRI,*)"<new LrO|Ri>", DDOT(N,EVLRO,1,EVRi,1)
      WRITE(LUPRI,*)"<new LiO|Rr>", DDOT(N,EVLiO,1,EVRR,1)
      call flshfo(lupri)

      RETURN 
      END
!--------
      SUBROUTINE CC_BIORTOCOMP2(N,EVRR,EVRI,EVLR,EVLI,EVLRO,EVLIO)
C--------------------------------------------------------------------
C     SC,09/12/2010
C     Biorthonormalization of complex left eigenvectors
C     Create new left vectors eigenvectors such that they are biorthog
C     to right eigenvectors. Based on BIORT1. 
C     Complex input left eigenv to dot on R is taken as <Lr|-i<Li|
C     Complex left eigenv in output is taken as <Lr_o|+i<Li_o|
C
C     EVRR: Right eigenvector, Real part. 
C     EVRI: Right eigenvector, Imaginary part. 
C     EVLR: Eigenvector, left,  real part. 
C     EVLI: Eigenvector, left,  imag. part. 
C     EVLRO: Eigenvectors, left,  real part. OUTPUT
C     EVLIO: Eigenvectors, left,  imag. part. OUTPUT
C     N : length of vector
C     (we pass one eigenvector at a time)
C--------------------------------------------------------------------
C
      IMPLICIT NONE
#include "priunit.h"
      INTEGER N 
#if defined (SYS_CRAY) 
      REAL EVRR(N),EVRI(N),EVLR(N),EVLI(N),EVLRO(N),EVLIO(N) 
      REAL CR,CI,CNSQ,DDOT
#else 
      !input
      DOUBLE PRECISION EVRR(N),EVRI(N),EVLR(N),EVLI(N)
      !output
      Double precision EVLRO(N),EVLIO(N)
      !local
      DOUBLE PRECISION CR,CI,CNSQ,DDOT
#endif 

      write(lupri,*)'INSIDE BIORTCOMP1'
      call flshfo(lupri)
      write(lupri,*)'<rr|rr>=',ddot(N,EVRR,1,EVRR,1)
      write(lupri,*)'<ri|ri>=',ddot(N,EVRI,1,EVRI,1)
      write(lupri,*)'<lr|lr>=',ddot(N,EVLR,1,EVLR,1)
      write(lupri,*)'<li|li>=',ddot(N,EVLI,1,EVLI,1)
C
      write(lupri,*)'<lr|Rr>=',ddot(N,EVLR,1,EVRR,1)
      write(lupri,*)'<li|Ri>=',ddot(N,EVLI,1,EVRI,1)
      write(lupri,*)'<lr|Ri>=',ddot(N,EVLR,1,EVRI,1)
      write(lupri,*)'<li|Rr>=',ddot(N,EVLI,1,EVRR,1)
      !<ReL|ReR>+<ImL|ImR>
      CR = DDOT(N,EVLR,1,EVRR,1)+DDOT(N,EVLI,1,EVRI,1)
      write(lupri,*)'CR=', CR
      !<ReL|ImR>-<ImL|ReR>
      CI = DDOT(N,EVLR,1,EVRI,1)-DDOT(N,EVLI,1,EVRR,1)
      write(lupri,*)'CI=', CI
      CNSQ=CR*CR+CI*CI
      write(lupri,*)'CNSQ', CNSQ
      write(lupri,*)'CI/CNSQ', CI/CNSQ, 'CR/CNSQ', Cr/Cnsq
C
      CALL DZERO(EVLRO,N)
      CALL DZERO(EVLIO,N)
C
      !call flshfo(lupri)
      CALL DCOPY(N,EVLR,1,EVLRO,1)
      CALL DCOPY(N,EVLI,1,EVLIO,1)
      CALL DSCAL(N,(CR/CNSQ),EVLRO,1) 
      CALL DSCAL(N,(-CR/CNSQ),EVLIO,1) 
C
      CALL DAXPY(N,(-CI/CNSQ),EVLI,1,EVLRO,1)
      CALL DAXPY(N,(-CI/CNSQ),EVLR,1,EVLIO,1)
!
      WRITE(LUPRI,*)"This should be one ", 
     *      DDOT(N,EVLRO,1,EVRR,1)-DDOT(N,EVLIO,1,EVRI,1) 
      call flshfo(lupri)
      WRITE(LUPRI,*)"This should be zero ", 
     *      DDOT(N,EVLRO,1,EVRI,1)+DDOT(N,EVLIO,1,EVRR,1) 
      call flshfo(lupri)

      WRITE(LUPRI,*)"<new LrO|Rr>", DDOT(N,EVLRO,1,EVRR,1)
      WRITE(LUPRI,*)"<new LiO|Ri>", DDOT(N,EVLiO,1,EVRi,1)
      WRITE(LUPRI,*)"<new LrO|Ri>", DDOT(N,EVLRO,1,EVRi,1)
      WRITE(LUPRI,*)"<new LiO|Rr>", DDOT(N,EVLiO,1,EVRR,1)
      call flshfo(lupri)

      RETURN 
      END
!--------
      SUBROUTINE CC_LANCZOS_Ftrans(ioptrw,Lchain,Lspace,IFQTRAN,XVEC,
     &                            EVR,Ftrans,WORK,LWORK,ER,EI,ISYHOP)
C-------------------------------------------------------------------
C    SC,03/11/2010
C    Read in the FQ vectors and generate the trans^F matrix 
!    Lchain is length of Chain
!    Lspace is length of exc. space
!    Xvec is QR
!    EVR is eigenvector
!    ER is (real part of) eigenvalue
!    EI is (imaginary part of) eigenvalue
C--------------------------------------------------------------------
C
      IMPLICIT NONE
#include "priunit.h"
#include "ccsdsym.h"
      INTEGER Lchain,Lspace,IFQTRAN(3,*),LWORK,ioptrw 
#if defined (SYS_CRAY) 
      REAL XVEC(Lspace,Lchain),EVR(Lchain,Lchain),ER(Lchain), EI(Lchain)
      REAL Ftrans(Lchain,Lchain), FQX(Lchain,Lchain)
      REAL DDOT, WORK(LWORK)
#else 
      DOUBLE PRECISION XVEC(Lspace,Lchain),EVR(Lchain,Lchain),ER(Lchain)
      DOUBLE PRECISION EI(Lchain)
      DOUBLE PRECISION Ftrans(Lchain,Lchain), FQX(Lchain,Lchain)
      DOUBLE PRECISION DDOT, WORK(LWORK)
      DOUBLE PRECISION ONE,ZERO
#endif 
      INTEGER kFq1,kFq2,kend,lend,isyhop,ivec, n2vec
      CHARACTER*(10) MODEL
      PARAMETER (ONE=1.0d0,ZERO=0.0d0)


!      CALL AROUND('trans^F print input stuff')
!      DO J = 1, Lchain
!      WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
!     *   ' - R EIGENVALUE NO.',J, ' - EIGENVALUE',ER(J),
!     *   ' - R EIGENVECTOR :'
!         WRITE (LUPRI,'(5F15.8)') (EVR(I,J),I=1,Lchain)
!      WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
!     *   ' - EIGENVALUE NO.',J, ' - EIGENVALUE',ER(J),
!     *   ' - QR EIGENVECTOR :'
!         WRITE (LUPRI,'(5F15.8)') (XVEC(I,J),I=1,Lspace)
!      END DO
      n2vec = 1
      if (ioptrw.eq.1) n2vec=0

      CALL DZERO(FQX,Lchain*LCHAIN)
      kFq1 = 1
      !kfq2 = kFq1 + NT1AMX 
      !kend = kFq2 + NT2AMX
      kfq2 = kFq1 + NT1AM(isyhop) 
      kend = kFq2 + NT2AM(isyhop)
      lend = lwork - kend
      if (lend.lt.0) call quit('Insufficient memory for trans^F')
      !isyhop=1
      do i = 1, Lchain
         ivec = ifqtran(3,i)
         CALL CC_RDRSP('FQ ',IVEC,ISYHOP,IOPTRW,MODEL,
     &                  WORK(KFQ1),WORK(KFQ2))
         
!        CALL AROUND('FQ transformed vector ')
!         CALL CC_PRP(WORK(KFQ1),WORK(KFQ2),ISYHOP,1,N2VEC)

         do j = 1, Lchain
            if (EI(j).ne.zero) then
               !write(lupri,*) 'WARNING: FQX of cmplx, EI(j)=', j, EI(j)
            end if
            FQX(j,i) = ddot(Lspace,WORK(KFQ1),1,XVEC(1,j),1)
         end do
      end do
!      CALL AROUND('--- The FQQR matrix ---')
!      CALL OUTPUT(FQX,1,Lchain,1,Lchain,Lchain,Lchain,1,LUPRI)


      CALL DGEMM('N','N',Lchain,Lchain,Lchain,
     *            ONE,FQX,Lchain,
     *            EVR,Lchain,ZERO,
     *            Ftrans,Lchain)
!      CALL AROUND('--- The trans^F matrix ---')
!      CALL OUTPUT(Ftrans,1,Lchain,1,Lchain,Lchain,Lchain,1,LUPRI)

      RETURN 
      END

C  /* Deck cc_pram_lanczos */
      SUBROUTINE CC_PRAM_lanczos(CAM,PT1,ISYMTR,LGRS,work,lwork)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     Based on Ove's cc_pram
C
C     Purpose: Writes out vector:
C              %T1 and %T2 and ||T1||/||T2||
C     makes an analysis of contributions from given 
C     occupied orbital
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
C
      dimension work(lwork)
      PARAMETER (TWO = 2.0D0, THPRT = 1.0D-9, HUNDRED = 100.0D0)
      PARAMETER (PT1THR = 75.0D0)
C
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
#include "priunit.h"
Cholesky
#include "maxorb.h"
#include "ccdeco.h"
C
      LOGICAL lxray
      LOGICAL CCSEFF
Cholesky
C
C
      LOGICAL LGRS
      DIMENSION CAM(*)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
Cholesky
      CCSEFF = CCS .OR. (CHOINT.AND.CC2)
Cholesky
C
C------------------------
C     Add up the vectors.
C------------------------
C
      C1NOSQ = 0.0D0
      C2NOSQ = 0.0D0
      KC1 = 1
      KC2 = KC1 + NT1AM(ISYMTR)
      !<t1|t1>
      C1NOSQ = DDOT(NT1AM(ISYMTR),CAM(KC1),1,CAM(KC1),1)
Chol  IF (.NOT. CCS) C2NOSQ = DDOT(NT2AM(ISYMTR),CAM(KC2),1,CAM(KC2),1)
      IF (.NOT. CCSEFF)
      !<t2|t2>
     &   C2NOSQ = DDOT(NT2AM(ISYMTR),CAM(KC2),1,CAM(KC2),1)
      CNOSQ  = C1NOSQ + C2NOSQ
C
      IF (.NOT. (CCSEFF.OR.CCP2)) THEN
C
         WRITE(LUPRI,'(//10X,A)')
     *     'CC_PRAM:Overall Contribution of the Different Components'
         WRITE(LUPRI,'(10X,A//)')
     *     '--------------------------------------------------------'
         WRITE(LUPRI,'(/10X,A,10X,F10.4,A)')
     *     'Single Excitation Contribution <t1|t1>/<t|t>*100: ',
     *     (C1NOSQ/CNOSQ)*HUNDRED,' %'
         WRITE(LUPRI,'(/10X,A,10X,F10.4,A)')
     *     'Double Excitation Contribution <t2|t2>/<t|t>*100: ',
     *     (C2NOSQ/CNOSQ)*HUNDRED,' %'
         WRITE(LUPRI,'(/10X,A,10X,F10.4)')
     *     '||T1||/||T2||=sqrt(<t1|t1>/<t2|t2>)   : ',
     *      SQRT(C1NOSQ/C2NOSQ)
         IF (LGRS) WRITE(LUPRI,'(/10X,A,10X,F10.4)')
     *     'Tau1 diagnostic                : ',
     *      SQRT(C1NOSQ/(TWO*DFLOAT(NRHFT)))
         PT1 = (C1NOSQ/CNOSQ)*HUNDRED
      ELSE
         PT1 = HUNDRED
      ENDIF
      WRITE(LUPRI,'(/10X,A,10X,F10.4)')
     *  'Norm of Total Amplitude Vector : ',SQRT(CNOSQ)
C
      CALL FLSHFO(LUPRI)
C
C----------------------------------------------
C     Initialize threshold etc from Printlevel.
C----------------------------------------------
C
      NL = MAX(1,2*IPRINT)
C
      CNOSQ = MAX(CNOSQ,THPRT)
C
      THR1 = SQRT(C1NOSQ/CNOSQ)/NL
      THR2 = SQRT(C2NOSQ/CNOSQ)/NL
      THR1 = MAX(THR1,1.0D-02)
      THR2 = MAX(THR2,1.0D-02)
      SUMOFP = 0.0D00
C
      IF (DEBUG) THR1 = 0.0D0
C
C-------------------------------------
C     Loop until a few is Printed out.
C-------------------------------------
C
C
C---------------------------------------
C     Loop through One excitation part.
C---------------------------------------
C
      WRITE(LUPRI,'(//A)')
     *     ' +=============================================='
     *    //'===============================+'
      WRITE(LUPRI,'(1X,A)')
     *     '| symmetry|  orbital index  |   Excitation Numbers'
     *     //'             |   Amplitude  |'
      WRITE(LUPRI,'(1X,A)')
     *     '|  Index  |   a   b   i   j |      NAI      NBJ |'
     *     //'     NAIBJ    |              |'
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      ISYMAI = MULD2H(ISYMTR,ISYMOP)
C
  1   CONTINUE
      N1 = 0
C
      DO 100 ISYMA = 1,NSYM
C
         ISYMI = MULD2H(ISYMAI,ISYMA)
C
         DO 110 I = 1,NRHF(ISYMI)
C
            MI = IORB(ISYMI) + I
C
            DO 120 A=1,NVIR(ISYMA)
C
               NAI = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
C
               MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
               IF (ABS(CAM(NAI)) .GE. THR1 ) THEN
C
                  WRITE(LUPRI,9990) ISYMA,ISYMI,A,I,NAI,CAM(NAI)
C
                  N1 = N1 + 1
                  SUMOFP = SUMOFP + CAM(NAI)*CAM(NAI)
C
               ENDIF
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      IF ((N1.LT.1).AND.(SQRT(C1NOSQ/CNOSQ).GT.1.0D-3)) THEN
         THR1 = THR1/5.0D0
         GOTO 1
      ENDIF

      !skod test with no symmetry
!      maxlength=0
!      do isym=1,nsym
!         maxnocc=max(maxnocc,nrhf(isym))
!      end do

      IF (PT1.ge.pt1thr) then
        
!         maxocc = 0
!         do isym=1,nsym
!            maxocc = max(maxocc,nrhf(isym))
!         end do

         kstart = 1
         kend   = nrhft
         lend   = lwork - kend
C
         call dzero(work(kstart),nrhft)
!      DO I = 1,NRHFt
!         MI = IORB(1) + I
!         NA = NVIRT*(I-1) + 1
!         work(mi) = ddot(nvirt,CAM(NA),1,CAM(NA),1)
!         write(lupri,*) 'iocc=', mi, ' sum_a=', work(mi)
!      end do
         WRITE (LUPRI,*)
     &   'Important occupied orbital contributions (normalized)'
         WRITE (LUPRI,*)
     &   'i, mi, contribution  contrib/t1norm**2'

         DO ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            DO I  = 1,NRHF(ISYMI)
               MI = IORB(ISYMI) + I
               NA = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + 1
               work(mi) = 
     &         ddot(nvir(isyma),CAM(NA),1,CAM(NA),1)
               write(lupri,*) i, mi, work(mi), work(mi)/C1NOSQ
            end do
         end do
      !call FINDMAXN(X,IXLEN,IMAX,NMAX)
C
        CALL FLSHFO(LUPRI)
      end if
C
C--------------------------------------------
C     Loop through Double excitation vector.
C     If not ccs or ccp2
C--------------------------------------------
C
      IF (.NOT. ( CCSEFF .OR. CCP2 )) THEN
C
      WRITE(LUPRI,'(A)')
     *     ' +----------------------------------------------'
     *    //'-------------------------------+'
C

 2    CONTINUE
      N2 = 0
C
      DO 200 ISYMAI = 1,NSYM
C
         ISYMBJ = MULD2H(ISYMAI,ISYMTR)
C
         DO 210 ISYMJ = 1,NSYM
C
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
C
            DO 220 ISYMI = 1,NSYM
C
               ISYMA = MULD2H(ISYMI,ISYMAI)
C
               DO 230 J = 1,NRHF(ISYMJ)
C
                  MJ = IORB(ISYMJ) + J
C
                  DO 240 B = 1,NVIR(ISYMB)
C
                     NBJ = IT1AM(ISYMB,ISYMJ)
     *                   + NVIR(ISYMB)*(J - 1) + B
C
                     MB = IORB(ISYMB) + NRHF(ISYMB) + B
C
                     DO 250 I = 1,NRHF(ISYMI)
C
                        MI = IORB(ISYMI) + I
C
                        DO 260 A = 1,NVIR(ISYMA)
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
C
                           MA = IORB(ISYMA) + NRHF(ISYMA) +  A
C
                           IF (((ISYMAI.EQ.ISYMBJ).AND.
     *                         (NAI .LT. NBJ)).OR.(ISYMAI.LT.ISYMBJ))
     *                          GOTO 260
C
                           IF (ISYMAI.EQ.ISYMBJ) THEN
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                             + INDEX(NAI,NBJ)
                           ELSE
                              NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                            + NT1AM(ISYMAI)*(NBJ-1) + NAI
                           ENDIF
C
                           KAIBJ = NAIBJ + NT1AM(ISYMTR)
                           IF (ABS(CAM(KAIBJ)) .GT. THR2 ) THEN
C
                              WRITE(LUPRI,9991) ISYMA,ISYMB,ISYMI,ISYMJ,
     *                                      A,B,I,J,NAI,NBJ,NAIBJ,
     *                                      CAM(KAIBJ)
                              N2 = N2 + 1
C
                              SUMOFP = SUMOFP + CAM(KAIBJ)*CAM(KAIBJ)
C
                           ENDIF
C
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
  230          CONTINUE
  220       CONTINUE
  210    CONTINUE
  200 CONTINUE
C
      IF ((N2.LT.1).AND.(SQRT(C2NOSQ/CNOSQ).GT.1.0D-3)) THEN
         THR2 = THR2/5D00
         GOTO 2
      ENDIF
C
      ENDIF
C
      WRITE(LUPRI,'(A)')
     *     ' +=============================================='
     *    //'===============================+'
C
      WRITE(LUPRI,'(//10X,A,8X,F10.4)')
     *     'Norm of Printed Amplitude Vector : ',SQRT(SUMOFP)
      WRITE(LUPRI,'(//10X,A43,1X,F9.6)')
     *     'Printed all single excitations greater than',THR1
      IF (.NOT. (CCSEFF.OR.CCP2)) THEN
         WRITE(LUPRI,'(//10X,A43,1X,F9.6)')
     *     'Printed all double excitations greater than',THR2
      ENDIF
C
 9990 FORMAT(1X,'| ',I1,3X,I1,2X,' | ',I3,5X,I3,4X,' | ',I8,9x,
     *       ' | ',12x,' | ',1x, F10.6,'  |')
 9991 FORMAT(1X,'| ',I1,1X,I1,1X,I1,1X,I1,' | ',
     *       I3,1X,I3,1X,I3,1X,I3,' | ',
     *       I8,1x,I8,' | ',I12,' | ',1x,F10.6,'  |')
C
      RETURN
      END
C=================================================================
      SUBROUTINE CC_build_Rfull(LUMAT,FMAT,
     &                          Lchain,Lspace,XVEC,EVR,
     &                          ER,EI,
     &                          WORK,LWORK)
C--------------------------------------------------------------------
C    SC,10/03/2011: generate pseudo eigenvectors in full space by:
C    - reading in the Q (or PT) matrix
!    - ddotting with Lanczos eigenvector
!    Lchain is length of Chain
!    Lspace is length of full excitation space
!    Xvec is R (or L) eigenvector in full space - OUTPUT
!    EVR is Lanczos eigenvector - INPUT
!    ER is lanczos eigenvalue (real part), EI is lanczos eigenvalue (im part)
!    LUMAT/FMAT are the unit number and file name of the Q (PT) matrix
C--------------------------------------------------------------------
C
      IMPLICIT NONE
#include "priunit.h"
#include "ccsdsym.h"
      INTEGER Lchain,Lspace,LWORK,LUMAT 
      CHARACTER FMAT*(*)
#if defined (SYS_CRAY) 
      REAL XVEC(Lspace*Lchain),EVR(Lchain*Lchain)
      !REAL XVEC(Lspace,Lchain),EVR(Lchain,Lchain)
      REAL ER(Lchain), EI(Lchain)
      REAL DDOT, WORK(LWORK)
#else 
      !DOUBLE PRECISION XVEC(Lspace,Lchain),EVR(Lchain,Lchain)
      DOUBLE PRECISION XVEC(Lspace*Lchain),EVR(Lchain*Lchain)
      DOUBLE PRECISION ER(Lchain), EI(Lchain)
      DOUBLE PRECISION DDOT, WORK(LWORK)
      DOUBLE PRECISION ONE,ZERO
#endif 
      INTEGER kend,lend,isyhop,ivec, n2vec
      INTEGER KQMAT, IOFF
      CHARACTER*(10) MODEL
      PARAMETER (ONE=1.0d0,ZERO=0.0d0)


      KQMAT = 1 
      KEND  = kQMAT + Lspace*Lchain
      LEND  = LWORK - KEND
      IF (LEND .LT. 0) THEN
         WRITE(LUPRI,*) 'Needed:', KEND, 'Available:', LWORK
         CALL QUIT('Insufficient work space in CC_build_Rfull')
      END IF

      CALL DZERO(WORK(KQMAT),Lspace*Lchain)

      !Read in the entire Q matrix 
      ioff = 0
      do k=1,Lchain
         CALL CC_RVEC(LUMAT,FMAT,Lspace,Lspace,K,WORK(KQMAT+ioff))
         ioff = ioff + Lspace
      enddo
!      write(lupri,*)'---- '
!      write(lupri,*)' The Q (or P) matrix '
!      call OUTPUT (WORK(KQMAT),1,Lspace,1,Lchain,Lspace,Lchain,
!     *                   1,LUPRI)
!      write(lupri,*)'---- '

      !Compute R_full = Qmat*R_lanczos_mat (Qmat dim N*J, R_l_mat dim J*J)
      CALL DGEMM('N','N',Lspace,Lchain,Lchain,
     *            ONE,WORK(KQMAT),Lspace,
     *            EVR,Lchain,ZERO,
     *            XVEC,Lspace)

      write(lupri,*)'Un-normalized R_full matrix computed'
      call flshfo(lupri)

      RETURN 
      END
!--------------------------------------------------------------------
      SUBROUTINE CC_LANCZOS_SUMRULES(LABEL,J,ER,EI,R,L,F,U,V)
C--------------------------------------------------------------------
C     Given a set of Lanzcos T eigenvalues and eigenvector this
C     calculates the sum rules 
C     J: Chain Length 
C     ER: Eigenvalues, real part 
C     EI: Eigenvalues, imag part 
C     R: Eigenvectors right (stored as cols)
C     L: Eigenvectors Left  (stored as cols, e.i. eigvectors of T^T)
C     F: F-matrix "trans" (not double q-transformed) 
C     U, V are normalization factors
C--------------------------------------------------------------------
      IMPLICIT NONE
#include "priunit.h"
#include "cclrlancz.h"
#include "codata.h"
     
      INTEGER J,IFREQ,NFREQ
      DOUBLE PRECISION U,V
      DOUBLE PRECISION ER(J),EI(J),R(J,J),L(J,J),F(J,J)
      DOUBLE PRECISION ZERO, ONE, TWO, PSUMK
      DOUBLE PRECISION sum_s0,sum_s0_f,sum_s0_t(-6:2)
      DOUBLE PRECISION sum_l0,sum_l0_f,sum_l0_t(-6:2)
      DOUBLE PRECISION mean_exc_energy(-6:2)
      DOUBLE PRECISION FACTOR,FACTOR2

      LOGICAL NOIMAG
      INTEGER I, K, IDAMP, N
!
      PARAMETER (ZERO = 0.0D0, ONE  = 1.0D0, TWO = 2.0D0)
      CHARACTER*(8) LABEL
!
      LOGICAL SKIPNEXT 
      SKIPNEXT=.FALSE.
C
C     Loop over chain 
C
      WRITE (LUPRI,*) " FOR OPERATOR LABEL = ", LABEL

      !WRITE (LUPRI,*) 'check FREQ=', FREQ
      CALL DZERO(SUM_S0_T,9)
      CALL DZERO(SUM_L0_T,9)
      CALL DZERO(mean_exc_energy,9)

      DO N=-6,2
      SUM_S0=ZERO
      SUM_L0=ZERO
      SUM_S0_F=ZERO
      SUM_L0_F=ZERO
      DO I=1,J
         IF (SKIPNEXT) THEN
            SKIPNEXT = .FALSE.
         ELSE
            !check on imaginary eigenvalue = 0
            IF (EI(I).EQ.ZERO)  THEN
               FACTOR=U*V*R(1,I)*L(1,I)
C              No imaginary eigenvector, calc eta*t conts 
               SUM_S0 = SUM_S0 + FACTOR*ER(I)**(N+1)
               SUM_L0 = SUM_L0 + FACTOR*(ER(I)**(N+1))*DLOG(ER(I))
C
C              Include F-part cont 
C
               FACTOR = V*V*L(1,I)*ER(I)**(N+1)
               DO K=1,J
                  FACTOR2 = FACTOR*L(1,K)*F(I,K)
                  if (abs(ei(k)).eq.zero) then
                     !denominator
                     PSUMK  = ER(I)+ER(K)
                     SUM_S0_F = SUM_S0_F - FACTOR2/PSUMK
                     SUM_L0_F = SUM_L0_F - FACTOR2*DLOG(ER(I))/PSUMK
                  else
                     write(lupri,*)'LANCZOS, do not know, index=',k
                  end if
                ENDDO

                ELSE
                    SKIPNEXT=.true.
                ENDIF
              end if !skipnext
C
      ENDDO

      SUM_S0_T(N) = (SUM_S0 + SUM_s0_F)*TWO
      SUM_L0_T(N) = (SUM_L0 + SUM_l0_F)*TWO
      END DO
      DO N=-6,2
         if (SUM_S0_T(N).eq.zero) then
           mean_exc_energy(N) = zero
         else
           mean_exc_energy(N) = xtev*dexp(SUM_L0_T(N)/SUM_S0_T(N))
         end if
      END DO

      CALL HEADER('Oscillator strength sum rules',30)
      WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,A,A,/)')
     &   'S(N) Sum Rules : Dipole Length Approximation in a.u.',
     &   'N',LABEL(1:1),LABEL(1:1),' - component'
         WRITE (LUPRI,'(9(5X,I3,(4X,G13.6)/))')
     &         (N,SUM_S0_T(N),N=-6,2)
      WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,A,A,/)')
     &   'L(N) Sum Rules : Dipole Length Approximation in a.u.',
     &   'N',LABEL(1:1),LABEL(1:1),' - component'
         WRITE (LUPRI,'(9(5X,I3,(4X,G13.6)/))')
     &         (N,SUM_L0_T(N),N=-6,2)
      WRITE (LUPRI,'(//,14X,A,/,6X,A,5X,A,A,A,/)')
     &   'I(N) Sum Rules : Dipole Length Approximation in eV ',
     &   'N',LABEL(1:1),LABEL(1:1),' - component'
         WRITE (LUPRI,'(9(5X,I3,(4X,G13.6)/))')
     &         (N,mean_exc_energy(N),N=-6,2)

      RETURN
      END
!--------------------------------------------------------------------
      SUBROUTINE CC_LANCZOS_OSCSTR(LABEL,J,ER,EI,R,L,F,U,V)
C--------------------------------------------------------------------
C     Given a set of Lanzcos T eigenvalues and eigenvector this
C     calculates the oscillator strengths
C     J: Chain Length 
C     ER: Eigenvalues, real part 
C     EI: Eigenvalues, imag part 
C     R: Eigenvectors right (stored as cols)
C     L: Eigenvectors Left  (stored as cols, e.i. eigvectors of T^T)
C     F: F-matrix "trans" (not double q-transformed) 
C     U, V are normalization factors
C--------------------------------------------------------------------
      IMPLICIT NONE
#include "priunit.h"
#include "cclrlancz.h"
#include "codata.h"
     
      INTEGER J,IFREQ,NFREQ
      DOUBLE PRECISION U,V
      DOUBLE PRECISION ER(J),EI(J),R(J,J),L(J,J),F(J,J)
      DOUBLE PRECISION F0I(J)
      DOUBLE PRECISION ZERO, ONE, TWO, PSUMK, THREE
      DOUBLE PRECISION sum_s0,sum_s0_f
      DOUBLE PRECISION FACTOR,FACTOR2

      LOGICAL NOIMAG
      INTEGER I, K, IDAMP, N
!
      PARAMETER (ZERO = 0.0D0, ONE  = 1.0D0, TWO = 2.0D0,THREE=3.0d0)
      CHARACTER*7 FNEXCIL
      PARAMETER (FNEXCIL='CCEXCIL')
      INTEGER IDUMMY, LUEXCIL

      CHARACTER*(8) LABEL
!
      LOGICAL SKIPNEXT 
      SKIPNEXT=.FALSE.
C
C     Loop over chain 
C
      WRITE (LUPRI,*) " FOR OPERATOR LABEL = ", LABEL
      LUEXCIL=-1
      CALL GPOPEN(LUEXCIL,FNEXCIL,'NEW',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
!
      CALL DZERO(F0I,J)

      DO I=1,J
        IF (SKIPNEXT) THEN
          SKIPNEXT = .FALSE.
        ELSE
          !check on imaginary eigenvalue = 0
          IF (EI(I).EQ.ZERO)  THEN
            FACTOR=U*V*R(1,I)*L(1,I)
C           No imaginary eigenvector, calc eta*t conts 
            F0I(I) = FACTOR*ER(I)
C           Include F-part cont 
            FACTOR = V*V*L(1,I)*ER(I)
            SUM_S0_F = ZERO
            DO K=1,J
               FACTOR2 = FACTOR*L(1,K)*F(I,K)
               if (ei(k).eq.zero) then
                  !denominator
                  PSUMK  = ER(I)+ER(K)
                  SUM_S0_F = SUM_S0_F - FACTOR2/PSUMK
               else
                  write(lupri,*)'LANCZOS, do not know, index=',k
               end if
             ENDDO
             F0I(I) = (F0I(I) + SUM_S0_F)*TWO/THREE

          ELSE
             SKIPNEXT=.true.
          ENDIF
        end if !skipnext
C
      ENDDO

      !CALL HEADER('Oscillator strengths',30)
      WRITE (LUEXCIL,'(14X,A,/,6X,A,8X,A,A,A)')
     &   'Oscillator strengths: Dipole Length Approximation in a.u.',
     &   'I',LABEL(1:1),LABEL(1:1),' - component'
!         WRITE (LUPRI,'(J(5X,I3,(4X,G13.6)/))')
!     &         (I,ER(I),F0I(I),I=1,J)
      WRITE (LUEXCIL,*) 'T MATRIX EIGENVALUES (eV) FOR JCHAIN = ',J
      WRITE (LUEXCIL,*)
     & 'INDEX & REAL (eV) & IMAG (eV) & REAL (au) & IMAG (au) & OSCST'
      do i=1,j
         write(luexcil,'(X,I4,5F16.8)') i, ER(i)*XTEV,
     &                  EI(I)*XTEV, ER(I), EI(I), F0I(i)
      end do

      CALL GPCLOSE(LUEXCIL,'KEEP')
      RETURN
      END
!----------------
C=================================================================
      SUBROUTINE CC_build_Rfull2(LUMAT,FMAT,
     &                          nexci,Lchain,Lspace,
     &                          XVEC,EVR,ER,EI,
     &                          WORK,LWORK)
C--------------------------------------------------------------------
C    SC,05/2011: generate pseudo eigenvectors in full space by:
C    - reading in the Q (or PT) matrix
!    - multiplying with Lanczos eigenvector
!    Nexci is the nr of vectors in full space that will be generated
!    Lchain is length of Chain
!    Lspace is length of full excitation space
!    Xvec is R (or L) eigenvector matrix in full space - OUTPUT
!    EVR is Lanczos eigenvectors - INPUT
!    ER is lanczos eigenvalue (real part), EI is lanczos eigenvalue (im part)
!    LUMAT/FMAT are the unit number and file name of the Q (PT) matrix
C--------------------------------------------------------------------
      IMPLICIT NONE
#include "priunit.h"
#include "ccsdsym.h"
      INTEGER Lchain,Lspace,LWORK,LUMAT,nexci
      CHARACTER FMAT*(*)
#if defined (SYS_CRAY) 
      REAL XVEC(Lspace*nexci),EVR(Lchain*nexci)
      REAL ER(nexci), EI(nexci)
      REAL DDOT, WORK(LWORK)
#else 
      DOUBLE PRECISION XVEC(Lspace*nexci),EVR(Lchain*nexci)
      DOUBLE PRECISION ER(nexci), EI(nexci)
      DOUBLE PRECISION DDOT, WORK(LWORK)
      DOUBLE PRECISION ONE,ZERO
#endif 
      INTEGER kend,lend,isyhop,ivec, n2vec
      INTEGER KQMAT, IOFF
      CHARACTER*(10) MODEL
      PARAMETER (ONE=1.0d0,ZERO=0.0d0)
      KQMAT = 1
      KEND  = kQMAT + Lspace*Lchain
      LEND  = LWORK - KEND
      IF (LEND .LT. 0) THEN
         WRITE(LUPRI,*) 'Needed:', KEND, 'Available:', LWORK
         CALL QUIT('Insufficient work space in CC_build_Rfull')
      END IF

      CALL DZERO(WORK(KQMAT),Lspace*Lchain)

      !Read in the entire Q matrix 
      ioff = 0
      do k=1,Lchain
         CALL CC_RVEC(LUMAT,FMAT,Lspace,Lspace,K,WORK(KQMAT+ioff))
         ioff = ioff + Lspace
      enddo
!      write(lupri,*)'---- '
!      write(lupri,*)' The Q (or P) matrix '
!      call OUTPUT (WORK(KQMAT),1,Lspace,1,Lchain,Lspace,Lchain,
!     *                   1,LUPRI)
!
!      write(lupri,*)'---- '

      !Compute R_full = Qmat*R_lanczos_mat (Qmat dim N*J, R_l_mat dim J*nexci)
      CALL DGEMM('N','N',Lspace,nexci,Lchain,
     *            ONE,WORK(KQMAT),Lspace,
     *            EVR,Lchain,ZERO,
     *            XVEC,Lspace)

      write(lupri,*)'Un-normalized R/L_full matrix block computed'
      call flshfo(lupri)

      RETURN
      END

C=================================================================
      SUBROUTINE CC_check_resid(ECURR,iside,nexc,Lspace,isyhop,
     &                          XVEC,ER,EI,
     &                          WORK,LWORK,APROXR12,N2VEC,
     &                          tstomega)
C--------------------------------------------------------------------
C    Check the residual on the pseudo eigenvectors in full space by:
C    - Performing AR_i-w_iR
!    - computing norm of resulting vector
!    Nexci is the nr of vectors in full space that will be generated
!    Lspace is length of full excitation space
!    Xvec is the matrix of R (or L) eigenvectors
!    ER is lanczos eigenvalue (real part), EI is lanczos eigenvalue (im part)
C--------------------------------------------------------------------
      IMPLICIT NONE
#include "priunit.h"
#include "ccsdsym.h"
#include "dummy.h"
      INTEGER Lspace,LWORK,nexc,iside
      DOUBLE PRECISION XVEC(Lspace,nexc)
      DOUBLE PRECISION ER(nexc), EI(nexc)
      DOUBLE PRECISION DDOT, WORK(LWORK)
      DOUBLE PRECISION ONE,ZERO,xnorm,ECURR
      DOUBLE PRECISION omega1,Rdag_R,Rdag_A_R,xnorm2

      INTEGER katrr,kend,lend,isyhop,ivec,n2vec
      INTEGER IOFF
      logical tstomega
      CHARACTER*(10) MODEL
      CHARACTER*(3) APROXR12
      PARAMETER (ONE=1.0d0,ZERO=0.0d0)

      KATRR = 1
      KEND  = KATRR + Lspace
      LEND  = LWORK - KEND
      IF (LEND .LT. 0) THEN
         WRITE(LUPRI,*) 'Needed:', KEND, 'Available:', LWORK
         CALL QUIT('Insufficient work space in CC_check_resid')
      END IF

      !Transform one excitation vector at a time.
      !ioff = 0
      do k=1,nexc
         call dcopy(Lspace,XVEC(1,k),1,work(katrr),1)
         CALL cc_atrr(ECURR,ISYHOP,ISIDE,WORK(Katrr),Lend,.false.,
     &                dummy,APROXR12,.false.)
!      SUBROUTINE CC_ATRR(ECURR,ISYMV,ISIDE,WORK,LWORK,LCONT,CONT,
!     &                   APROXR12)

         !WRITE(LUPRI,*) 'AR vector'
         !CALL CC_PRP(WORK(KAtrr),WORK(KAtrr+nt1amx),1,1,N2VEC)
         !call dcopy(Lspace,XVEC(1,k),1,work(kend),1)
         !call dscal(Lspace,ER(k),work(kend),1)
         !WRITE(LUPRI,*) 'wR vector'
         !CALL CC_PRP(WORK(Kend),WORK(Kend+nt1amx),1,1,N2VEC)
         !CALL DAXPY(Lspace,-1.0d0,WORK(Kend),1,WORK(Katrr),1)
         !WRITE(LUPRI,*) 'AR-wR vector'
         !CALL CC_PRP(WORK(KAtrr),WORK(KAtrr+nt1amx),1,1,N2VEC)
         CALL DAXPY(Lspace,-ER(k),XVEC(1,k),1,WORK(Katrr),1)
         xnorm=dsqrt(ddot(Lspace,WORK(KAtrr),1,WORK(KAtrr),1))
         write(lupri,*)'Norm residual eival=', er(k), ' is =', xnorm
         !ioff = ioff + Lspace
      enddo
      !
      ! Now try and compute omega'=<R|A|R>/<R|R> and Res=AR-w'R
      !
      if (tstomega) then
      do k=1,nexc
         call dcopy(Lspace,XVEC(1,k),1,work(katrr),1)
         CALL cc_atrr(ECURR,ISYHOP,ISIDE,WORK(Katrr),Lend,.false.,
     &                dummy,APROXR12,.false.)
         Rdag_A_R=ddot(Lspace,XVEC(1,k),1,WORK(KAtrr),1)
         Rdag_R=ddot(Lspace,XVEC(1,k),1,XVEC(1,k),1)
         omega1=Rdag_A_R/Rdag_R
         !WRITE(LUPRI,*) 'AR vector'
         !CALL CC_PRP(WORK(KAtrr),WORK(KAtrr+nt1amx),1,1,N2VEC)
         !call dcopy(Lspace,XVEC(1,k),1,work(kend),1)
         !call dscal(Lspace,ER(k),work(kend),1)
         !WRITE(LUPRI,*) 'wR vector'
         !CALL CC_PRP(WORK(Kend),WORK(Kend+nt1amx),1,1,N2VEC)
         !CALL DAXPY(Lspace,-1.0d0,WORK(Kend),1,WORK(Katrr),1)
         !WRITE(LUPRI,*) 'AR-wR vector'
         !CALL CC_PRP(WORK(KAtrr),WORK(KAtrr+nt1amx),1,1,N2VEC)
         CALL DAXPY(Lspace,-omega1,XVEC(1,k),1,WORK(Katrr),1)
         xnorm2=dsqrt(ddot(Lspace,WORK(KAtrr),1,WORK(KAtrr),1))
         write(lupri,*)'Norm residual omega1=', omega1,' is =',xnorm2
         write(lupri,*)'Diff |omega1-ER(k)|=', abs(omega1-ER(k))
         !ioff = ioff + Lspace
      enddo
      end if
      RETURN
      END

C  /* Deck cc_normalize */
      SUBROUTINE cc_normalize(LENGTH,N,ZR,ZL,ER,EI,ISYHOP,IOPT)
C
C  normalize pseudo RE
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION ZR(LENGTH,N), ZL(LENGTH,N)
      DIMENSION ER(N), EI(N)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
C
C     Normalize RE eigenvectors 
C
      if (iopt.eq.1) then
      DO I = 1,N
        ZNRM = DDOT(LENGTH,ZR(1,I),1,ZR(1,I),1)
        !write(lupri,*)'Norm RrRr', ZNRM
        ZNRM = D1 / SQRT(ZNRM)
        !IMAX = IDAMAX(N,ZR(1,I),1)
        !IF (ZR(IMAX,I) .LT. D0) ZNRM = -ZNRM
        CALL DSCAL(LENGTH,ZNRM,ZR(1,I),1)
        !
!        WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
!     *  ' - XR EIGENVALUE NO.',I, ' - EIGENVALUE',ER(I),
!     *  ' - XR EIGENVECTOR :'
!        CALL CC_PRAM(ZR(1,I),PT1,ISYHOP,.FALSE.)
        !
      END DO
      end if
C
C     binormalize LE eigenvectors 
C
      if (iopt.eq.2) then
      DO I = 1,N
        ZNRM = DDOT(LENGTH,ZL(1,I),1,ZR(1,I),1)
        write(lupri,*)'Norm LrRr of vector nr:',I,' is=', ZNRM
        ZNRM = D1 / (ZNRM)
        !IMAX = IDAMAX(N,ZR(1,I),1)
        !IF (ZR(IMAX,I) .LT. D0) ZNRM = -ZNRM
        CALL DSCAL(LENGTH,ZNRM,ZL(1,I),1)
        !
        ZNRM = DDOT(LENGTH,ZL(1,I),1,ZR(1,I),1)
        !write(lupri,*)'Check Norm LrRr of vector nr:',I,' is=', ZNRM
!        WRITE (LUPRI,'(/A,I4,A,1P,D15.6/A)')
!     &        '- LPT EIGENVALUE NO.',I,' - EIGENVALUE',
!     &                        ER(I),
!     &        '- LPT EIGENVECTOR :'
!        CALL CC_PRAM(ZL(1,I),PT1,ISYHOP,.FALSE.)

      END DO
      end if
      RETURN
      END
C=================================================================
      SUBROUTINE CC_dump_onfile(length,nex,XVEC,isyhop,lista,
     &                          work,lwork,model,ioptrw)
C--------------------------------------------------------------------
!    Dump on file, according to lista type
C--------------------------------------------------------------------
      IMPLICIT NONE
#include "priunit.h"
#include "ccsdsym.h"
#include "dummy.h"
      INTEGER length,nex,isyhop,lwork,ioptrw
#if defined (SYS_CRAY) 
      REAL XVEC(length,nex),work(lwork)
      REAL DDOT
#else 
      DOUBLE PRECISION XVEC(length,nex), work(lwork)
      DOUBLE PRECISION DDOT
      DOUBLE PRECISION ONE,ZERO
#endif 
      INTEGER kend,lend,ivec,n2vec
      INTEGER IOFF
      CHARACTER*(10) MODEL
      CHARACTER*(3) LISTA
      PARAMETER (ONE=1.0d0,ZERO=0.0d0)

      !Dump
      kend=1
      lend=lwork-kend
      IF (LEND .LT. 0) THEN
         WRITE(LUPRI,*) 'Needed:', KEND, 'Available:', LWORK
         CALL QUIT('Insufficient work space in CC_dump_onfile')
      END IF
      do I=1,nex
         CALL CC_WRRSP(lista,I,ISYHOP,IOPTRW,MODEL,WORK(KEND),
     *                 XVEC(1,I),XVEC(1+nt1am(isyhop),I),
     *                 WORK(kend),Lend)
      enddo
      RETURN
      END
C  /* Deck cc_resid_lanczos */
      SUBROUTINE cc_resid_lanczos(beta,LENGTHQ,N_RE,Qi,ZR,ER,EI)
C
C  Check residual according to Lanczos recipe
C  N_RE is the number of Lanczos eigenvectors, each of them long N_RE
C  (squared matrix). We use the last element of each Lanczos eigenvector.
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION ZR(N_RE,N_RE), Qi(lengthq)
      DIMENSION ER(N_RE), EI(N_RE)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
C
      xqnorm = ddot(lengthq,Qi,1,Qi,1)
      xqnorm = sqrt(xqnorm)
      write(lupri,*)'norm of Q_i+1 is ', xqnorm
      write(lupri,*)'beta_jchain=', beta
      

      DO I = 1,N_RE
        residREi = beta*xqnorm*abs(ZR(N_RE,I))
        write(lupri,*)'Residual norm Lanczos wi:',ER(i),' =', 
     &                 residREi
      END DO
C
      RETURN
      END

!----------------------
C  /* Deck ccbiort */
      SUBROUTINE CCBIORT(LUQMAT,FQMAT,LUPMAT,FPMAT,ichain,
     &                   NPREV,NCCVAR,
     &                   QNEW,PNEW,ISYMTR,THRLDP,WRK,LWRK)
C 
C Modified version of CCORT to biorthogonalize p^t and q
C vectors of Lanczos chain
C two-sided modified Gram-Schmidt (TSMGS) process
C Notice we assume to enter the routine with just 1 Q_new and 1 P_new
C
C Purpose:
C  Bi-orthogonalize new p and q vectors against all previous p 
C  and q vectors and among themselves, and renormalize.
C  (Orthogonalization is performed twice if round-off is large,
C   if larger than THRRND).
C
C Input:
C  LUQMAT - file with previous Q vectors with name FQMAT
C  LUPMAT - file with previous P vectors with name FPMAT
C  QNEW, PNEW - non-orthogonal q and p new vectors
C  NBPREV - number of previous q and p vectors on FQMAT/FPMAT
C  THRLDP - threshold for linear dependence
C
C Output:
C  QNEW and PNEW on output are biorthogonal to previous ones
C
C Scratch:
C  require WRK of length NCCVAR
C
#include "implicit.h"
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
      CHARACTER FQMAT*(*),FPMAT*(*)
      DIMENSION QNEW(NCCVAR),PNEW(NCCVAR), WRK(NCCVAR)
      LOGICAL LOCDBG
      PARAMETER (THRRND=1.D-4, THRTT=1.D-6, THRZER=1.0D-35)
      PARAMETER (THRTST=1.0D-12)
      PARAMETER (D1 = 1.0D0)
      PARAMETER (LOCDBG=.false.)
      
      INTEGER SVEC, RVEC, LVEC
C
C
      IF (LOCDBG) THEN
         WRITE (LUPRI,'(//A/)') ' --- Debug output from CCORT ---'
      END IF
      KEND1  = 1
      LWRK1  = LWRK - KEND1
      IF (LWRK1.LT.0) THEN
        WRITE(LUPRI,*) 'LWRK, KEND1: ',LWRK, KEND1
        CALL QUIT('Insufficient memory in CCORT')
      END IF
C
C     Check bi-norm of input Q and P vectors 
C     Issue warning if norm .le. THRZER
C
      TTPQ = DDOT(NCCVAR,PNEW(1),1,QNEW(1),1)
!      WRITE (LUPRI,*) ' Initial overlap <Pnew|Qnew>:'
!     &                 ,TTPQ, 'for ichain=', ichain
      !THRZER=1.0D-35
      IF (TTPQ.LT.THRZER) THEN
         WRITE(LUPRI,*) 'CCBIORTHO: WARNING!! <pi|qi> <= 0'
      END IF
C
C work space allocation
C
      KQPRE  = KEND1
      KPPRE  = KQPRE + NCCVAR
      KEND1  = KPPRE + NCCVAR
      IROUND = 0
      ITURN  = 0

 1500 CONTINUE
      ITURN  = ITURN + 1
C
C     Orthogonalize new q and p vectors against previous ones
C
      DO K = 1, NPREV
         CALL CC_RVEC(LUQMAT,FQMAT,nccvar,
     *               nccvar,K,WRK(KQPRE))
         CALL CC_RVEC(LUPMAT,FPMAT,nccvar,
     *               nccvar,K,WRK(KPPRE))

         !<Pold|Qnew>
         TTQ   = DDOT(NCCVAR,WRK(KPPRE),1,QNEW(1),1)
         !<Pnew|Qold>
         TTP   = DDOT(NCCVAR,PNEW(1),1,WRK(KQPRE),1)
         !<Pold|Qold>
         TTDQP = DDOT(NCCVAR,WRK(KPPRE),1,WRK(KQPRE),1)
         IF (LOCDBG) THEN
           WRITE (LUPRI,*)'start <P(',K,')|Qnew(',ichain,')> =', TTQ
           WRITE (LUPRI,*)'start <Pnew(',ichain,')| Q(',K,')> =', TTP
           WRITE (LUPRI,*)'start <P(',K,')|Q(',K,')>=', TTDQP
         END IF
         !Qnew = |Qnew> - |Qold>*<Pold|Qnew>
         CALL DAXPY(NCCVAR,-TTQ,WRK(KQPRE),1,QNEW(1),1)
         !<Pnew| = <Pnew| - <Pnew|Qold><Pold|
         CALL DAXPY(NCCVAR,-TTP,WRK(KPPRE),1,PNEW(1),1)
C ortho test
         TTQ   = DDOT(NCCVAR,WRK(KPPRE),1,QNEW(1),1)
         TTP   = DDOT(NCCVAR,PNEW(1),1,WRK(KQPRE),1)
         TTDQP = DDOT(NCCVAR,WRK(KPPRE),1,WRK(KQPRE),1)
         if (LOCDBG) then
          WRITE(LUPRI,*)'BIORT ortho: <P(',k,')|Qnew(',ichain,')>=', TTQ
          WRITE(LUPRI,*)'BIORT ortho: <Pnew(',ichain,')|Q(',k,')>=', TTP
          WRITE(LUPRI,*)'Biort ortho: <P(',K,')|Q(',K,')>=', TTDQP
         end if
      END DO
C
C     Bi-normalize new vectors 
C     (and maybe check if vectors are linear dependent?) 
C
      TTN   = DDOT(NCCVAR,PNEW(1),1,QNEW(1),1)
      WRITE(LUPRI,*) 'BIORT check ortho: <Pnew|Qnew>(',ichain,')=', TTN
      TTqnn  = DDOT(NCCVAR,QNEW(1),1,QNEW(1),1)
      WRITE(LUPRI,*) 'BIORT check ortho: <Qnew|Qnew>(',ichain,')=',TTqnn
!      IF (abs(TTN) .LT. THRRND) IROUND = IROUND+1
!      SQRTTN = SQRT(TTN)
!      IF (SQRTTN.LT.THRTT) THEN
!         CALL DSCAL(NCCVAR,(D1/THRTT),PNEW(1),1)
!         TT = DDOT(NCCVAR,PNEW(1),1,QNEW(1),1)
!         SQRTT = SQRT(TT)
!      END IF
      CALL DSCAL(NCCVAR,(D1/TTN),PNEW(1),1)
      !factor=sign(ttn)*(D1/sqrt((abs(ttn)))
      !CALL DSCAL(NCCVAR,factor,PNEW(1),1)
      !CALL DSCAL(NCCVAR,factor,QNEW(1),1)
C
      !IF (IROUND.GT.0  .AND. ITURN.EQ.1 ) GO TO 1500
      IF (ITURN.EQ.1 ) GO TO 1500

C
C *** End of subroutine CCBiORT
C
      IF (LOCDBG) THEN
        WRITE (LUPRI,'(/A//)')' --- End debug output from CCBiORT ---'
      END IF

      RETURN
      END
*=============================================
C  /* Deck ccbiortPQ */
      SUBROUTINE CCBIORTPQ(LUQMAT,FQMAT,LUPMAT,FPMAT,ichain,
     &                   NPREV,NCCVAR,
     &                   QNEW,PNEW,ISYMTR,THRLDP,WRK,LWRK)
C 
C Modified version of CCORT to biorthogonalize p^t and q
C vectors of Lanczos chain
C two-sided modified Gram-Schmidt (TSMGS) process
C Notice we assume to enter the routine with just 1 Q_new and 1 P_new
C NORMALIZE BY SHARING PQ NORM ON BOTH P AND Q
C
C Purpose:
C  Bi-orthogonalize new p and q vectors against all previous p 
C  and q vectors and among themselves, and renormalize.
C  (Orthogonalization is performed twice if round-off is large,
C   if larger than THRRND).
C
C Input:
C  LUQMAT - file with previous Q vectors with name FQMAT
C  LUPMAT - file with previous P vectors with name FPMAT
C  QNEW, PNEW - non-orthogonal q and p new vectors
C  NBPREV - number of previous q and p vectors on FQMAT/FPMAT
C  THRLDP - threshold for linear dependence
C
C Output:
C  QNEW and PNEW on output are biorthogonal to previous ones
C
C Scratch:
C  require WRK of length NCCVAR
C
#include "implicit.h"
#include "priunit.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
C
      CHARACTER FQMAT*(*),FPMAT*(*)
      DIMENSION QNEW(NCCVAR),PNEW(NCCVAR), WRK(NCCVAR)
      PARAMETER (THRRND=1.D-4, THRTT=1.D-6, THRZER=1.0D-35)
      PARAMETER (THRTST=1.0D-12)
      PARAMETER (D1 = 1.0D0)
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .false.)
      
      INTEGER SVEC, RVEC, LVEC
C
C
      IF (LOCDBG) THEN
         WRITE (LUPRI,'(//A/)') ' --- Debug output from CCORT ---'
      END IF
      KEND1  = 1
      LWRK1  = LWRK - KEND1
      IF (LWRK1.LT.0) THEN
        WRITE(LUPRI,*) 'LWRK, KEND1: ',LWRK, KEND1
        CALL QUIT('Insufficient memory in CCORT')
      END IF
C
C     Check bi-norm of input Q and P vectors 
C     Issue warning if norm .le. THRZER
C
      TTPQ = DDOT(NCCVAR,PNEW(1),1,QNEW(1),1)
!      WRITE (LUPRI,*) ' Initial overlap <Pnew|Qnew>:'
!     &                 ,TTPQ, 'for ichain=', ichain
      !THRZER=1.0D-35
      IF (TTPQ.LT.THRZER) THEN
         WRITE(LUPRI,*) 'CCBIORTHO: WARNING!! <pi|qi> <= 0'
      END IF
C
C work space allocation
C
      KQPRE  = KEND1
      KPPRE  = KQPRE + NCCVAR
      KEND1  = KPPRE + NCCVAR
      IROUND = 0
      ITURN  = 0

 1500 CONTINUE
      ITURN  = ITURN + 1
C
C     Orthogonalize new q and p vectors against previous ones
C
      DO K = 1, NPREV
         CALL CC_RVEC(LUQMAT,FQMAT,nccvar,
     *               nccvar,K,WRK(KQPRE))
         CALL CC_RVEC(LUPMAT,FPMAT,nccvar,
     *               nccvar,K,WRK(KPPRE))

         !<Pold|Qnew>
         TTQ   = DDOT(NCCVAR,WRK(KPPRE),1,QNEW(1),1)
         !<Pnew|Qold>
         TTP   = DDOT(NCCVAR,PNEW(1),1,WRK(KQPRE),1)
         !<Pold|Qold>
         TTDQP = DDOT(NCCVAR,WRK(KPPRE),1,WRK(KQPRE),1)
         IF (LOCDBG) THEN
           WRITE (LUPRI,*)'start <P(',K,')|Qnew(',ichain,')> =', TTQ
           WRITE (LUPRI,*)'start <Pnew(',ichain,')| Q(',K,')> =', TTP
           WRITE (LUPRI,*)'start <P(',K,')|Q(',K,')>=', TTDQP
         END IF
         !Qnew = |Qnew> - |Qold>*<Pold|Qnew>
         CALL DAXPY(NCCVAR,-TTQ,WRK(KQPRE),1,QNEW(1),1)
         !<Pnew| = <Pnew| - <Pnew|Qold><Pold|
         CALL DAXPY(NCCVAR,-TTP,WRK(KPPRE),1,PNEW(1),1)
C ortho test
         TTQ   = DDOT(NCCVAR,WRK(KPPRE),1,QNEW(1),1)
         TTP   = DDOT(NCCVAR,PNEW(1),1,WRK(KQPRE),1)
         TTDQP = DDOT(NCCVAR,WRK(KPPRE),1,WRK(KQPRE),1)
         if (LOCDBG) then
          WRITE(LUPRI,*)'BIORT ortho: <P(',k,')|Qnew(',ichain,')>=', TTQ
          WRITE(LUPRI,*)'BIORT ortho: <Pnew(',ichain,')|Q(',k,')>=', TTP
          WRITE(LUPRI,*)'Biort ortho: <P(',K,')|Q(',K,')>=', TTDQP
         end if
      END DO
C
C     Bi-normalize new vectors 
C     (and maybe check if vectors are linear dependent?) 
C
      TTN   = DDOT(NCCVAR,PNEW(1),1,QNEW(1),1)
      WRITE(LUPRI,*) 'BiORTPQ check <Pnew|Qnew>(',ichain,')=', TTN
      TTqnn  = DDOT(NCCVAR,QNEW(1),1,QNEW(1),1)
      WRITE(LUPRI,*) 'BiORTPQ check <Qnew|Qnew>(',ichain,')=',TTqnn
!      IF (abs(TTN) .LT. THRRND) IROUND = IROUND+1
!      SQRTTN = SQRT(ABS(TTN))
!      IF (SQRTTN.LT.THRTT) THEN
!         CALL DSCAL(NCCVAR,(D1/THRTT),PNEW(1),1)
!         TT = DDOT(NCCVAR,PNEW(1),1,QNEW(1),1)
!         SQRTT = SQRT(TT)
!      END IF
!      CALL DSCAL(NCCVAR,(D1/TTN),PNEW(1),1)
!
      !distribute norm between the two vectors
      xn=D1/sqrt((abs(ttn)))
      factor=sign(xn,ttn)
      write(lupri,*)'BiOrthPQ: factor 1/sqrt(|<p|q>|)=',factor
      CALL DSCAL(NCCVAR,factor,PNEW(1),1)
      CALL DSCAL(NCCVAR,factor,QNEW(1),1)
      write(lupri,*)'BiOrthPQ 2: <p|q>=', 
     &               DDOT(NCCVAR,PNEW(1),1,QNEW(1),1)
      write(lupri,*)'BiOrthPQ 2: <q|q>=', 
     &               DDOT(NCCVAR,QNEW(1),1,QNEW(1),1)
C
      !IF (IROUND.GT.0  .AND. ITURN.EQ.1 ) GO TO 1500
      IF (ITURN.EQ.1 ) GO TO 1500

C
C *** End of subroutine CCBiORTPQ
C
      IF (LOCDBG) THEN
        WRITE (LUPRI,'(/A//)')' --- End debug output from CCBiORTPQ ---'
      END IF

      RETURN
      END
*=============================================
C  /* Deck rgord2 */
      SUBROUTINE RGORD2(NM,N,WR,WI,ZR,ZL,LSORT)
C
C  15-Aug-1989 Hans Joergen Aa. Jensen
C  modified by Ove Christiansen 20-dec 1994
C  to normalize correct!
C  modified by Sonia Coriani 28-oct 2010 to
C  to sort both left and right
C
C  Reorder (no normalize!) eigenpairs from DGEEV
C          (they are already normalized in DGEEV)
C  LSORT = .FALSE. --> ascending order
C  LSORT = .TRUE.  --> descending order
#include "implicit.h"

      LOGICAL LSORT
      DIMENSION WR(N), WI(N), ZR(NM,N), ZL(NM,N)
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0 )
      DO 100 I = 1,N-1
         ILOW = I
         XLOW = WR(I)
         DO 50 J = I+1,N
          IF (.NOT.LSORT) THEN
            !descending order
            IF (WR(J) .LT. XLOW) THEN
               ILOW = J
               XLOW = WR(J)
            END IF
          ELSE
            IF (WR(J) .GT. XLOW) THEN
               ILOW = J
               XLOW = WR(J)
            END IF
          END IF
   50    CONTINUE
         IF (ILOW .NE. I) THEN
            WR(ILOW) = WR(I)
            WR(I)    = XLOW
            XLOW     = WI(ILOW)
            WI(ILOW) = WI(I)
            WI(I)    = XLOW
            CALL DSWAP(N,ZR(1,I),1,ZR(1,ILOW),1)
            CALL DSWAP(N,ZL(1,I),1,ZL(1,ILOW),1)
         END IF
  100 CONTINUE
C
C     Normalize eigenvectors (should already be normalized actually)
C
!      DO 200 I = 1,N
!        ZNRM = DDOT(N,ZR(1,I),1,ZR(1,I),1)
!        write(lupri,*)'Norm RrRr', ZNRM
!        ZNRM = D1 / SQRT(ZNRM)
!        IMAX = IDAMAX(N,ZR(1,I),1)
!        IF (ZR(IMAX,I) .LT. D0) ZNRM = -ZNRM
!        CALL DSCAL(N,ZNRM,ZR(1,I),1)
!        !
!        ZNRM = DDOT(N,ZL(1,I),1,ZL(1,I),1)
!        write(lupri,*)'Norm LrLr', ZNRM
!        ZNRM = D1 / SQRT(ZNRM)
!        IMAX = IDAMAX(N,ZL(1,I),1)
!        IF (ZL(IMAX,I) .LT. D0) ZNRM = -ZNRM
!        CALL DSCAL(N,ZNRM,ZL(1,I),1)
! 200 CONTINUE
      RETURN
      END
!
