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
!  #define SIRFCK_DEBUG
C  /* Deck sirfck */
      SUBROUTINE SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRECT,
     &                  WRK,LWORK)
C
C     March 1997 - tsaue Screening
C     921201-hjaaj; partly based on GETFCK
C     written by Henrik Koch and Trygve Helgaker 22-November-1991.
C
C PURPOSE : Driver routine for the calculation of the two-electron
C           part of the fock matrices from AO integrals.
C           We assume the densities and fock matrices are full
C           squares and without symmtry reduction.
C           The subroutine offers the possibility of calculating
C           the fock matrices directly or by reading the integrals
C           from disk.
C
C output:
C     FMAT   :  NDMAT Fock matrices with 2-el. cont. added
C  input:
C     FMAT   :  NDMAT Fock matrices
C     DMAT   :  NDMAT density matrices
C     NDMAT  :  number of density matrices and Fock matrices
C     ISYMDM :  point group symmetry of DMAT and FMAT
C     IFCTYP = +/-XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C       + sign: alpha + beta matrix (singlet)
C       - sign: alpha - beta matrix (triplet)
C
C     DIRECT = .TRUE.   Direct calculation using hermit.
C              .FALSE.  Read integrals from AOTWOINT.
Cholesky       Overriden by the CHOINT flag
C
#include <implicit.h>
      LOGICAL   DIRECT, DOERG, DOATR, LSRHYBR
      LOGICAL   DO_SPINDNS, DFT_SPINDNS_SAVE
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,*)
      DIMENSION ISYMDM(NDMAT), IFCTYP(NDMAT), MBAS(8)
      DIMENSION WRK(*), IFCTMP_SR(3)
      DIMENSION DFCOEFF(NBAST) !density-fitting coefficients
C
#include <dummy.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D4 = 4.0D0, 
     &           DMP5 = -0.5D0)
C
C Used from common blocks
C   GNRINF: SRINTS, CHIVAL
C   INFORB: NSYM, NBAS(:), N2BASX
C   INFTAP: LUINTA
C   INFPRI: IPRFCK
C   CBIREA: ATOMDF, NODDYDF
C
#include <maxorb.h>
#include <priunit.h>
#include <gnrinf.h>
#include <infinp.h>
#include <inforb.h>
#include <inftap.h>
#include <infpri.h>
#include <dftcom.h>
#include <dfterg.h>
#include <cbirea.h>
C
Cholesky
#include <ccdeco.h>
#include <sirftim.h>
#include <choscf.h>
Cholesky
C
#include <ibtfun.h>
C
      LOGICAL SFIRST, PREHF, IOPTVP, DO_CHOINT
      SAVE SFIRST
      DATA SFIRST  /.TRUE./
C
C
#include <memint.h>
      CALL QENTER('SIRFCK')
Cholesky
      DTIME  = SECOND()
      ITSLDA = .FALSE.
      ICOUN3 = 0
Cholesky
C
      iprsrdft = iprfck
      IF (NDMAT.EQ.0) GO TO 9999
C
C     Find out if we have spin-density terms
C     Find out if we shall and can use Cholesky
C
      DO_SPINDNS = .FALSE.
      DO_CHOINT = CHOINT .and. CHIVAL .lt. 0.0D0 ! no Cholesky for short-range integrals (yet)
      DO I = 1,NDMAT
         IF (IFCTYP(I) .LT. 0) DO_SPINDNS = .TRUE.
         IF (IFCTYP(I) .ne. 13 .and. IFCTYP(I) .ne. 12 .and.
     &       IFCTYP(I) .ne.  3 .and. IFCTYP(I) .ne.  2 .and.
     /       IFCTYP(I) .ne. 22 .and. IFCTYP(I) .ne. 11)
     &       DO_CHOINT = .FALSE.
      END DO
      IF (DO_SPINDNS .AND. .NOT. ADDSRI) THEN
         CALL QUIT('SIRFCK can only do spin dens. for srDFT')
      END IF
      IF (CHOINT .AND. .NOT. DO_CHOINT) THEN
         NINFO = NINFO + 1
         WRITE(LUPRI,'(3(/A),//A)')
     &   ' INFO: Cholesky decomposed two-electron integrals can',
     &   ' only be applied for IFCTYP equal to 2-3, 11-13 or 22.',
     &   ' You may have to implement the request yourself !',
     &   ' Proceeding with standard AO-integral direct '//
     &     'Fock matrix construction.'
      END IF
C
#ifdef SIRFCK_DEBUG
      iprsrdft = max(1,iprsrdft)  ! development code, always print
      write(lupri,*)'SIRFCK: entering'
      write(lupri,*)'SIRFCK: NDMAT=', NDMAT
      write(lupri,*)'SIRFCK: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
      write(lupri,*)'SIRFCK: DO_SPINDNS,DFT_SPINDNS = ',
     & DO_SPINDNS,DFT_SPINDNS
#endif
C
      IF (DO_SPINDNS) THEN
C        DC, DV, DS or DXC, DXV, DXS
         NSR_DMAT = NDMAT / 3
         NDCandDV = NSR_DMAT * 2
         IF (NSR_DMAT*3 .ne. NDMAT) THEN
C        ... check if programming was OK
            write(lupri,*)'Inconsistency(1), NSR_DMAT*3 ne NDMAT'
            write(lupri,*)'SIRFCK: NDMAT,NSR_DMAT=', NDMAT,NSR_DMAT
            write(lupri,*)'SIRFCK: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
            CALL QUIT('Inconsistency(1), NSR_DMAT*3 ne NDMAT')
         END IF
         IF (NSR_DMAT .gt. 1) THEN
C        off-sets below assume DXC1, DXC2, DXV1, DXV2
C        but the correct order is DXC1, DXV1, DXC2, DXV2
            call quit('hjaaj: not correct off-sets in sirfck')
C        REMEMBER ALSO zero the spin-density FMAT, because that is not done
C        in the calculation of the long-range contribution
         END IF
      ELSE IF (NASHT .EQ. 0 .OR. NDMAT .EQ. 1 .OR. .NOT. ADDSRI) THEN
C        ... NDMAT typically .eq. 1 for CI or CI-srDFT
C        DC or DXC (closed shell SCF) / Dtot
         NSR_DMAT = NDMAT
         NDCandDV = NDMAT
CFRAN JUST A TEST !! DEBUG REMOVE AGAIN WHEN LOCALSPIN HESSIAN IS PROPERLY IMPLEMENTED TODO !!!!!!!!
      ELSE IF(DFT_SPINDNS) THEN
C        DC, DV, DS or DXC, DXV, DXS
         NSR_DMAT = NDMAT / 3
         NDCandDV = NSR_DMAT * 2
         IF (NSR_DMAT*3 .ne. NDMAT) THEN
C        ... check if programming was OK
            write(lupri,*)'Inconsistency(1a), NSR_DMAT*3 ne NDMAT'
            write(lupri,*)'SIRFCK: NDMAT,NSR_DMAT=', NDMAT,NSR_DMAT
            write(lupri,*)'SIRFCK: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
            CALL QUIT('Inconsistency(1a), NSR_DMAT*3 ne NDMAT')
         END IF
CFRAN JUST A TEST END
      ELSE
#ifdef SIRFCK_DEBUG
CFRAN 17/02/2010
      WRITE(LUPRI,*) 'I do DC(and DV) or DXC(and DXV) '
#endif
C        DC, DV or DXC, DXV
         NSR_DMAT = NDMAT / 2
         NDCandDV = NDMAT
         IF (NSR_DMAT*2 .ne. NDMAT) THEN
            write(lupri,*)'Inconsistency(2), NSR_DMAT*2 ne NDMAT'
CFRAN JUST TO UNDERSTAND END
            write(lupri,*)'SIRFCK: NDMAT,NSR_DMAT=', NDMAT,NSR_DMAT
            write(lupri,*)'SIRFCK: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
            CALL QUIT('Inconsistency(2), NSR_DMAT*2 ne NDMAT')
         END IF
      END IF

C     ...PREHF: Find out if this is a HF calculation preceding
C        a Multiconfigurational SR-DFT calculation.
      PREHF = DOSCF.AND.(DOMCSRDFT.OR.DOCISRDFT).AND.
     &                  (.NOT.DOHFSRDFT).AND.(MCTYPE.EQ.0)
      LSRHYBR = SRHYBR .AND. .NOT.DOHFSRDFT .AND. .NOT.PREHF
C
#ifdef SIRFCK_DEBUG
CFRAN 17/02/2010
      WRITE(LUPRI,*) 'PREHF   = ', PREHF
      WRITE(LUPRI,*) 'LSRHYBR = ', LSRHYBR
CFRAN 17/02/2010 END
#endif
      HFXFACSAVE = HFXFAC
      CHIVALSAVE = CHIVAL
      IF (DO_CHOINT) THEN
C
C        *************** Cholesky section *************** 
C
C        ... Outermost loop over density matrices in order
C            to maximize work space for batches of Cholesky
C            vectors. We may want to change this to avoid
C            reading through the Cholesky file(s) NDMAT times.
C
         DO I = 1,NDMAT
C
           IFCSAV = IFCTYP(I)
           IF (IFCTYP(I) .EQ. 11)    THEN       ! No exchange
              IFCTYP(I) = 13
              ITSLDA = .TRUE.
           END IF
C
           JSYMDM = ISYMDM(I)
C
           IF (IFCTYP(I) .EQ. 13) THEN
               IOPTVP = .TRUE.                ! Coulomb included
C           ... singlet, (tot.) symmetric density matrices
               CALL CHO_FCKDS1(FMAT(1,I),DMAT(1,I),WRK,LWORK,
     &                         ISYMDM(I),IOPTVP)
C
           ELSE IF (IFCTYP(I) .EQ. 12) THEN
               IOPTVP = .FALSE.
C           ... singlet, (tot.) symmetric density matrices (exchange only)
               CALL CHO_FCKDS1(FMAT(1,I),DMAT(1,I),WRK,LWORK,
     &                         ISYMDM(I),IOPTVP)
C
           ELSE IF (IFCTYP(I) .EQ.  3) THEN
               IOPTVP = .TRUE.
C           ...  singlet, non-symmetric (general) density matrices
               ICOUN3 = ICOUN3 + 1
               CALL CHO_FCKDS3(FMAT(1,I),DMAT(1,I),WRK,LWORK,
     &                         ISYMDM(I),IOPTVP)
C
           ELSE IF ((IFCTYP(I) .EQ. 2) .OR. (IFCTYP(I) .EQ.  22)) THEN
               IOPTVP = .FALSE.
C           ...  singlet, non-symmetric density matrices (exchange only)
               ICOUN3 = ICOUN3 + 1
               CALL CHO_FCKDS3(FMAT(1,I),DMAT(1,I),WRK,LWORK,
     &                         ISYMDM(I),IOPTVP)
C
           ELSE IF (IFCTYP(I) .EQ.  0) THEN
               GO TO 201
           ELSE
               GO TO 8000
           END IF
C
           ITSLDA    = .FALSE.
           IFCTYP(I) = IFCSAV
C
  201     CONTINUE
         ENDDO
      ELSE IF (DIRECT) THEN ! not Cholesky, AO-direct?
C
C
C     ****************** AO direct Fock matrices ***************
C
C
         IF (SFIRST) THEN
            IPRHER = 2*IPRFCK
            CALL SETHER(IPRHER,.FALSE.,WRK,LWORK)
C           CALL SETHER(IPRHER,NEWGEO,WORK,LWORK)
            SFIRST = .FALSE.
         END IF
         DO I = 1,NDMAT
            ISYMDM(I) = ISYMDM(I) - 1
         END DO
         IPRHER = IPRFCK/100
C
C        Special code for short range/long range integrals
C        ... never short-range integrals in first call of HERFCK
C        ... if normal HF then normal 2-el ints (CHIVAL=0)
         SRINTS = .FALSE.
         IF (ADDSRI .OR. PREHF) THEN
C        ... full exchange for long range part and for PREHF
            HFXFAC = D1
            DO I = 1,NDCandDV
               J = MOD(IFCTYP(I),10)
               IF (IBTAND(J,2) .eq. 0) IFCTYP(I) = IFCTYP(I) + 2
C              ... make sure the "do exchange" bit is set in IFCTYP
            END DO
         END IF
         IF (PREHF) THEN
C           this is e.g. PREHF, i.e. standard HF
            CHIVAL = -1.0D0
         END IF
#ifdef SIRFCK_DEBUG
         write(lupri,*) 'sirfck: before calling herfck ...'
         write(lupri,*) 'SIRFCK: ADDSRI =', ADDSRI
         write(lupri,*) 'SIRFCK: calling HERFCK with CHIVAL ',CHIVAL
         write(lupri,*) 'SIRFCK: calling HERFCK with HFXFAC ',HFXFAC
         write(lupri,*) 'SIRFCK: ISYMDM ',(ISYMDM(I),I=1,NDMAT)
         write(lupri,*) 'SIRFCK: IFCTYP ',(IFCTYP(I),I=1,NDMAT)
!        write(lupri,*) 'SIRFCK: first DMAT'
!        call output(DMAT,1,nbast,1,nbast,nbast,nbast,-1,lupri)
#endif
C
       IF (CHIVAL .ne. D0) THEN
C        ... No long-range integrals if CHIVAL .eq. 0
         IF (NSYM.GT.1 .AND. (DODFT.OR.DOHFSRDFT.OR.DOMCSRDFT)) THEN
C           HERFCK destroys DMAT when called with symmetry
            KDMAT = 1
            KLAST = KDMAT + NDCandDV*N2BASX
            LWRK  = LWORK - KLAST + 1
            IF (KLAST.GT.LWORK) CALL STOPIT('SIRFCK','DMAT',KLAST,LWORK)
            CALL DCOPY(NDCandDV*N2BASX,DMAT,1,WRK(KDMAT),1)
            if (atomdf .AND. CHIVAL .LT. D0) then
C             hjaaj mar 2006: dens_fit not for short range integrals (yet)
              call dens_fit(fmat,wrk(kdmat),NDCandDV,ISYMDM,IFCTYP,
     &                    IPRHER,WRK(KLAST),LWRK)
            else
              CALL HERFCK(FMAT,WRK(KDMAT),NDCandDV,ISYMDM,IFCTYP,
     &                    IPRHER,WRK(KLAST),LWRK)
            end if
         ELSE
            if (atomdf .AND. CHIVAL .LT. D0) then
C             hjaaj mar 2006: dens_fit not for short range integrals (yet)
              call dens_fit(fmat,dmat,NDCandDV,ISYMDM,IFCTYP,IPRHER,
     &                      wrk,lwork)
            else
              CALL HERFCK(FMAT,DMAT,NDCandDV,ISYMDM,IFCTYP,IPRHER,
     &                      WRK,LWORK)
            end if
         END IF
       END IF ! (CHIVAL .ne. D0)
#ifdef SIRFCK_DEBUG
!        write(lupri,*) 'SIRFCK: first FMAT'
!        call output(FMAT,1,nbast,1,nbast,nbast,nbast,-1,lupri)
#endif
         HFXFAC = HFXFACSAVE
         CHIVAL = CHIVALSAVE
C
C        Add Short-range integrals (HF/MC-SRDFT hybrid model)
C
#ifdef SIRFCK_DEBUG
         write(lupri,*) 'add direct sr integrals ? ADDSRI = ', ADDSRI 
#endif 
         IF (ADDSRI .AND. (DOHFSRDFT.OR.DOMCSRDFT)) THEN
            IF (NASHT.GT.0 .AND. (.NOT.DOCISRDFT))THEN
C              ... get total DMAT = DC + DV if NASHT.gt.0
C                  for Hartree SR term
               IF (DO_SPINDNS) THEN
                  JFMAT = 3
               ELSE
                  JFMAT = 2
               END IF
               DO I = 1,NSR_DMAT
!!! new code !!!
C
                  CALL DAXPY(N2BASX,D1,DMAT(1,2*I),1,DMAT(1,2*I-1),1)
C                 ... add DV(I) to DC(I)
!
! Manu: the core one-index transformed and transition density matrices
!       now contains the total ones 

!!!!! inital code !!!
!                 CALL DAXPY(N2BASX,D1,DMAT(1,NSR_DMAT+I),1,DMAT(1,I),1)
C                 ... add DV(I) to DC(I)
!!!!!!!!!!!!!!!!!!!!!
               END DO
            ELSE
C              ... total DMAT = DC when NASHT.eq.0 or CISRDFT
               JFMAT   = 1
               IF (DO_SPINDNS) THEN
                  CALL QUIT('srDFT spin density for closed shell')
               END IF
            END IF
C
            IF (HFXFAC .NE. D0 .and. DO_SPINDNS) THEN
               CALL QUIT('"HFEXCH" not implemented yet'//
     &         ' for non-singlet MC-srDFT')
            END IF
C
            SRINTS = .TRUE.
C           HERFCK will zero FMAT with long range term, thus :
            KFSRAO = 1
            KLAST = KFSRAO + NSR_DMAT*N2BASX
            LWRK  = LWORK - KLAST + 1
            IF (HFXFAC .EQ. D0) THEN
C              remove any "do exchange" bits
               DO I = 1,NSR_DMAT
                  J = MOD(IFCTYP(I),10)
                  IF (IBTAND(J,2) .ne. 0) IFCTYP(I) = IFCTYP(I) - 2
               END DO
            END IF
            IF (KLAST.GT.LWORK) CALL STOPIT('SIRFCK','FMAT',KLAST,LWORK)
#ifdef SIRFCK_DEBUG
            write(lupri,*) 'SIRFCK: calling HERFCK with CHIVAL ',CHIVAL
            write(lupri,*) 'SIRFCK: calling HERFCK with HFXFAC ',HFXFAC
            write(lupri,*) 'SIRFCK: ISYMDM ',(ISYMDM(I),I=1,NDMAT)
            write(lupri,*) 'SIRFCK: IFCTYP ',(IFCTYP(I),I=1,NDMAT)
         do i = 1,ndmat+1
            write(lupri,*) 'SIRFCK: DMAT no.',i,' out of',ndmat
            call output(DMAT(1,i),1,nbast,1,nbast,nbast,nbast,-1,lupri)
         end do
#endif
            IF (NSYM.GT.1) THEN

!!!! start new code for symmetry !!!
                KDMAT = KLAST
                KLAST = KDMAT + NSR_DMAT*N2BASX
                LWRK  = LWORK - KLAST + 1
                IF (KLAST.GT.LWORK)
     &             CALL STOPIT('SIRFCK sr direct','DMAT',KLAST,LWORK)
C
                DO ISRMAT = 1, NSR_DMAT
                   IDMAT = JFMAT*(ISRMAT-1) +1
                   KIDMAT = KDMAT +(ISRMAT-1)*N2BASX 
                   CALL DCOPY(N2BASX,DMAT(1,IDMAT),1,WRK(KIDMAT),1)
!                  CALL DAXPY(N2BASX,D1,DMAT(1,IDMAT+1),1,
!    &                              WRK(KIDMAT),1)
                END DO 
C
                if (atomdf .AND. CHIVAL .LT. D0) then
C                  hjaaj mar 2006: dens_fit not for short range integrals (yet)
                   call dens_fit(wrk(KFSRAO),wrk(kdmat),NSR_DMAT,
     &                         ISYMDM,IFCTYP,IPRHER,WRK(KLAST),LWRK)
                else
                   CALL HERFCK(WRK(KFSRAO),WRK(KDMAT),NSR_DMAT,
     &                       ISYMDM,IFCTYP,IPRHER,WRK(KLAST),LWRK)
                end if

!!!! end new code for symmetry !!!

!!!!! start initial code !!!!!!!
C              HERFCK destroys DMAT when called with symmetry
C              We need DMAT for the "CALL SRDFT" below
!*             KDMAT = KLAST
!*             KLAST = KDMAT + NSR_DMAT*N2BASX
!*             LWRK  = LWORK - KLAST + 1
!*             IF (KLAST.GT.LWORK)
!*   &         CALL STOPIT('SIRFCK sr direct','DMAT',KLAST,LWORK)
!*             CALL DCOPY(NSR_DMAT*N2BASX,DMAT,1,WRK(KDMAT),1)
!*             if (atomdf .AND. CHIVAL .LT. D0) then
C               hjaaj mar 2006: dens_fit not for short range integrals (yet)
!*               call dens_fit(wrk(KFSRAO),wrk(kdmat),NSR_DMAT,
!*   &                         ISYMDM,IFCTYP,IPRHER,WRK(KLAST),LWRK)
!*             else
!*               CALL HERFCK(WRK(KFSRAO),WRK(KDMAT),NSR_DMAT,
!*   &                       ISYMDM,IFCTYP,IPRHER,WRK(KLAST),LWRK)
!*             end if
!!!!! end initial code !!!!!!!
            ELSE
C
C Manu - march 2011: allocation of additional memory 
C                    for total one-index transformed 
C                    and transition density matrices 
C                    in order to calculate the corresponding 
C                    srH Kernel contributions to the hessian. 
C                    The initial implementation, that used 
C                    directly DMAT, did not work when the two
C                    matrices were non-zero (which happens in the linear
C                    response when the trial vector has both an orbital AND a
C                    configuration part) because they were not ordered
C                    one after the other (DTVAO "stands" in between)  
C                         
C
#ifdef SIRFCK_DEBUG
                write(lupri,*) 'JFMAT =', JFMAT
#endif 
                   KDMAT = KLAST
                   KLAST = KDMAT + NSR_DMAT*N2BASX
                   LWRK  = LWORK - KLAST + 1
                   IF (KLAST.GT.LWORK)
     &             CALL STOPIT('SIRFCK sr direct','DMAT',KLAST,LWORK)
C
                   DO ISRMAT = 1, NSR_DMAT
                      IDMAT = JFMAT*(ISRMAT-1) +1
                      KIDMAT = KDMAT +(ISRMAT-1)*N2BASX 
                      CALL DCOPY(N2BASX,DMAT(1,IDMAT),1,WRK(KIDMAT),1)
!                     CALL DAXPY(N2BASX,D1,DMAT(1,IDMAT+1),1,
!    &                              WRK(KIDMAT),1)
                   END DO 

C
               if (atomdf .AND. CHIVAL .LT. D0) then
C               hjaaj mar 2006: dens_fit not for short range integrals (yet)
                 call dens_fit(wrk(KFSRAO),dmat,NSR_DMAT,
     &                         ISYMDM,IFCTYP,IPRHER,wrk(KLAST),LWRK)
               else
#ifdef SIRFCK_DEBUG
                   write(lupri,*) 'calling herfck with no symmetry!'     
                   write(lupri,*) 'NSR_DMAT = ', NSR_DMAT
#endif
C
C 
                       CALL HERFCK(WRK(KFSRAO),WRK(KDMAT),NSR_DMAT,
     &                       ISYMDM,IFCTYP,IPRHER,WRK(KLAST),LWRK)

! initial code (DMAT was used instead of WRK(KDMAT))
!
!*                     CALL HERFCK(WRK(KFSRAO),DMAT,NSR_DMAT,
!*   &                       ISYMDM,IFCTYP,IPRHER,WRK(KLAST),LWRK)
               end if
            END IF
C
            EJSR = DP5*DDOT(N2BASX,DMAT(1,1),1,WRK(KFSRAO),1)
            IF (NASHT.GT.0) THEN
             EJKVSR =-DP5*
     &                DDOT(N2BASX,DMAT(1,NSR_DMAT+1),1,WRK(KFSRAO),1)
            ELSE
             EJKVSR = D0
            END IF
 
            IF (JFMAT .EQ. 1) THEN
               CALL DAXPY(NSR_DMAT*N2BASX,D1,WRK(KFSRAO),1,FMAT,1)
            ELSE
C
               DO I = 1,NSR_DMAT
                  JFSRAO = KFSRAO + (I-1)*N2BASX
                  JF = JFMAT*(I-1)+ 1
                  CALL DAXPY(N2BASX,D1,WRK(JFSRAO),1,FMAT(1,JF),1)
               END DO
            END IF
            SRINTS = .FALSE.
            IF (NASHT.GT.0) THEN
               DO I = 1,NSR_DMAT
!!!! new code !!!!!!
                 CALL DAXPY(N2BASX,-D1,DMAT(1,2*I),1,DMAT(1,2*I-1),1)
C                 ... restore DC(I) = DTOT(I) -  DV(I)
!
! Manu: the core one-index transformed and transition density matrices
!       are now restored. 
!
!!!! initial code !!!
C                CALL DAXPY(N2BASX,-D1,DMAT(1,NSR_DMAT+I),1,DMAT(1,I),1)
C                 ... restore DC(I) = DTOT(I) -  DV(I)
!!!!!!!!!!!!!!!!!!!!!!!
               END DO
            END IF
         ENDIF
         DO I = 1,NDMAT
            ISYMDM(I) = ISYMDM(I) + 1
         END DO
#ifdef SIRFCK_DEBUG
           do i = 1,nsr_dmat
              write(lupri,*)
     &        'SIRFCK: direct short-range contribution to FMAT no.', i
           call output(WRK(KFSRAO+n2basx*(i-1)),1,nbast,1,nbast,
     &                 nbast,nbast,-1,lupri)
           end do
           write(lupri,*) 'SIRFCK: FMAT (lr + sr)'
           call output(FMAT,1,nbast,1,nbast,nbast,nbast,-1,lupri)
#endif
C
      ELSE ! not Cholesky and not AO-direct
C
C     ****************** "conventional" AO Fock matrices ***************
C                        (from AO integrals or supermatrix integrals on file)
C
         IF (HFXMU .NE. D0) CALL QUIT(
     &       'SIRFCK2: nonzero HFXMU requires .DIRECT calculation.') 
C
C        ... if SRHYBR (DFT-SRDFT hybrid) only include Coulomb terms 
C            (core is treated with BLYP -> hence no exchange here!)
         IF(LSRHYBR) THEN
            IFCTMP_SR(1) = IFCTYP(1)
            IFCTYP(1) = 11
C           ... FC,lr = JC,lr = J,lr * DC for SR-hybrid
            IF (NDMAT .EQ. 2) THEN
               IFCTMP_SR(2) = IFCTYP(2)
               IFCTYP(2) = 13
            END IF
C           ... FV,lr = G,lr * DV
            IF (NDMAT .GT. 2) CALL QUIT('NDMAT.gt.2 for SRHYBR & lr')
         ELSE
         END IF
         IF (ADDSRI .OR. PREHF) THEN
C        ... full exchange for long range part and for PREHF
            HFXFAC = D1
            DO I = 1,NDCandDV
               J = MOD(IFCTYP(I),10)
               IF (IBTAND(J,2) .eq. 0) IFCTYP(I) = IFCTYP(I) + 2
C              ... make sure the "do exchange" bit is set in IFCTYP
            END DO
         END IF
#ifdef SIRFCK_DEBUG
         write(lupri,*)'SIRFCK: calling SIRFCK2 "AOTWOINT" with'
         write(lupri,*)'SIRFCK: NDMAT=',NDMAT
         write(lupri,*)'SIRFCK: IFCTYP=',(IFCTYP(I),I=1,NDMAT)
!        write(lupri,*)'SIRFCK: first DMAT'
!        call output(DMAT,1,nbast,1,nbast,nbast,nbast,-1,lupri)
#endif
         IF (CHIVAL .ne. D0) THEN
C        ... no long-range integrals if chival .eq. 0
            CALL SIRFCK2(LUINTA,'AOTWOINT',
     &                   FMAT,DMAT,NDCandDV,ISYMDM,IFCTYP,WRK,LWORK)
          END IF ! (CHIVAL .ne. D0)
         if (lsrhybr) then
            do i = 1,ndmat
               ifctyp(i) = ifctmp_sr(i)
            end do
         end if
         HFXFAC = HFXFACSAVE
C
C        Add Short-range integrals to FC part of FMAT (SRDFT hybrid models)
C
         IF (ADDSRI) THEN
            IF (HFXFAC .NE. D0 .and. DO_SPINDNS) THEN
               CALL QUIT('"HFEXCH" not implemented yet'//
     &         ' for non-singlet non-direct MC-srDFT')
            END IF
            IF (NASHT.GT.0 .AND. (.NOT.DOCISRDFT))THEN
C              ... get total DMAT = DC + DV if NASHT.gt.0
C                  for Hartree SR term
               DO I = 1,NSR_DMAT
                  CALL DAXPY(N2BASX,D1,DMAT(1,NSR_DMAT+I),1,DMAT(1,I),1)
C                 ... add DV(I) to DC(I)
               END DO
            END IF
            IF (NSR_DMAT.GT.1) CALL QUIT('NSR_DMAT .GT. 1 for ADDSRI')
C
            IFCTMP_SR(1) = IFCTYP(1)
            IF (HFXFAC.NE.D0) THEN
C           ... for normal/presceding HF and for SRDFT with exact exch.
C               (add SR coul. and exch. to LR to give full 2-el part)
                IFCTYP(1) = 13
            ELSE 
C           ... add only Coul
C               (only SR Hartree term for SRDFT hybrids)
                IFCTYP(1) = 11
            END IF
            KFSRAO = 1
            KLAST  = KFSRAO + NSR_DMAT*N2BASX
            LWRK   = LWORK  - KLAST + 1
            IF (HFXFAC .EQ. D0) THEN
C              remove any "do exchange" bits
               DO I = 1,NSR_DMAT
                  J = MOD(IFCTYP(I),10)
                  IF (IBTAND(J,2) .ne. 0) IFCTYP(I) = IFCTYP(I) - 2
               END DO
            END IF
#ifdef SIRFCK_DEBUG
            write(lupri,*)'SIRFCK: calling SIRFCK2 "AOSR2INT" with'
            write(lupri,*)'SIRFCK: NSR_DMAT,NDMAT=',NSR_DMAT,NDMAT
            write(lupri,*)'SIRFCK: IFCTYP=',(IFCTYP(I),I=1,NDMAT)
!           write(lupri,*)'SIRFCK: first DMAT'
!           call output(DMAT,1,nbast,1,nbast,nbast,nbast,-1,lupri)
#endif
            CALL DZERO(WRK(KFSRAO),NSR_DMAT*N2BASX)
            CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFSRAO),
     &                   DMAT,NSR_DMAT,ISYMDM,IFCTYP,WRK(KLAST),LWRK)
            CALL DAXPY(N2BASX,D1,WRK(KFSRAO),1,FMAT,1)
            EJSR = DP5*DDOT(N2BASX,DMAT(1,1),1,WRK(KFSRAO),1)
            IF (NASHT.GT.0) THEN
               EJKVSR =-DP5*
     &            DDOT(N2BASX,DMAT(1,NSR_DMAT+1),1,WRK(KFSRAO),1)
            ELSE
               EJKVSR = D0
            END IF
            if (iprsrdft .gt. 0 .and.
     &         nsr_dmat.eq.1 .and. HFXFAC .NE. 0.0D0) then
               CALL DZERO(WRK(KFSRAO),N2BASX)
               IFCTYP_x = 12 ! exchange only
               CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFSRAO),
     &                   DMAT,NSR_DMAT,ISYMDM,IFCTYP_x,WRK(KLAST),LWRK)
               EKSR = DP5*DDOT(N2BASX,DMAT(1,1),1,WRK(KFSRAO),1)
               EJSR = EJSR - EKSR
            else
               EKSR = D0
            end if
            IF (NASHT.GT.0) THEN
               DO I = 1,NSR_DMAT
                  CALL DAXPY(N2BASX,-D1,DMAT(1,NSR_DMAT+I),1,
     &                                  DMAT(1,I)         ,1)
C                 ... restore DC(I) = DTOT(I) -  DV(I)
               END DO
            END IF
         END IF
      END IF
      HFXFAC = HFXFACSAVE
C
#ifdef OLD_DFT_CODE
!hjaaj: can probably be deleted, but want to be sure about srdft first. March 2006.
C
C     Add DFT contribution
C
      IF (DODFT) THEN
         DOERG = DFTADD
         DOATR = .NOT.DFTADD
         IF(DOATR)THEN
            J = NDMAT + 1
C           DTOTAO is put in DMAT(1,NDMAT+1) in TR1FCK.
C           We now need to copy DXAO to DMAT(1,NDMAT+2) for EXCMAT
            CALL DCOPY(N2BASX,DMAT,1,DMAT(1,J+1),1)
            NDMTMP = 2
            KSRAO  = 1
            KLAST  = KSRAO  + N2BASX
            LWRK  = LWORK - KLAST + 1
            CALL DZERO(WRK(KSRAO),N2BASX)
            CALL DFTEXC(DMAT(1,J),NDMTMP,WRK(KSRAO),1,
     &           .FALSE.,DOERG,.FALSE.,.FALSE.,DOATR,
     &           .FALSE.,DUMMY,.FALSE.,WRK(KLAST),LWRK,IPRFCK)
C           CALL HEADER('AO Fock matrix contribution in DFTATR /SIRFCK',-1)
C           CALL OUTPUT(WRK(KSRAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            CALL DAXPY(N2BASX,D1,WRK(KSRAO),1,FMAT,1)
         ELSE
            J = 1
            NDMTMP = 1
            CALL DFTEXC(DMAT(1,J),NDMTMP,FMAT,1,.FALSE.,DOERG,.FALSE.,
     &           .FALSE.,DOATR,.FALSE.,DUMMY,.FALSE.,WRK,LWORK,IPRFCK)
         ENDIF
C        ... DFTADD tells if add normal DFT Vxc or second order Vxc
         IF (NDMAT.GT.1) CALL QUIT('DFTADD and NDMAT.gt.1')
      END IF
#endif
C
C     ***************** ADD SR-DFT contribution if MCSCF-DFT or if HF-SRDFT
C
      IF (ADDSRI.AND.(DOHFSRDFT.OR.DOMCSRDFT).AND.(.NOT.PREHF)) THEN
#ifndef MOD_SRDFT
         call quit('srdft not included in this version')
#else
#ifdef SIRFCK_DEBUG
         do i = 1,ndmat+1
            write(lupri,*) 'SIRFCK srdft: DMAT no.',i,' out of',ndmat+1
            call output(DMAT(1,i),1,nbast,1,nbast,nbast,nbast,-1,lupri)
         end do
#endif
         IF (LSRHYBR .AND. DO_SPINDNS) THEN
            call quit('hybrid srDFT only for singlet MC-srDFT')
         END IF
         DOERG = DFTADD
         DOATR = .NOT.DFTADD
C        ... DFTADD tells if add normal DFT Vxc or second order Vxc
         IF (DO_SPINDNS) THEN
            NSRAO = 2
            NDOFF = 3
C           3 because DC, DV, DS or DXC, DXV, DXS
         ELSE IF (NASHT.GT.0 .AND. (.NOT.DOCISRDFT))THEN
C           ... get total DMAT = DC + DV if NASHT.gt.0
C               for xc SR term
            NSRAO = 1
            NDOFF = 2
C           2 because DC, DV or DXC, DXV
         ELSE
            NSRAO = 1
            NDOFF = 1
C           1 because DC or DXC
         END IF
         KSRAO = 1
         IF (LSRHYBR) NSRAO = NSRAO + 1
C        ... one more EXCMAT from SRDFT when LSRHYBR
         KLAST = KSRAO + NSRAO*N2BASX
         IF (DOERG .AND. NSR_DMAT .GT. 1) THEN
            WRITE(LUPRI,*) 'SRFCK.SRDFT only for one SRMAT when DOERG',
     &         NSR_DMAT
            CALL QUIT('SRFCK.SRDFT only for one SRMAT when DOERG')
         END IF
         IF (LSRHYBR) THEN
            IF (NDMAT .NE. 2) THEN
               write(lupri,*) 'ERROR LSRHYBR: NDMAT.ne.2 but is ',NDMAT
               CALL QUIT('NDMAT .ne. 2 for LSRHYBR')
            END IF
            NDOFF  = 2
         END IF
         LWRK  = LWORK - KLAST + 1
         IF (DOATR) THEN
            JDAO = NDMAT + 1
C           ... DTOTAO is in DMAT(1,NDMAT+1) from TR1FCK or RSPOLI
C           For SRHYBR we also need DVAO in DMAT(1,NDMAT+2) from TR1FCK.:
            IF (LSRHYBR) THEN
               JDXAO = JDAO + 2
            ELSE
               JDXAO = JDAO + 1
            END IF
            IF (DO_SPINDNS) THEN
C              DT, DXT, DS, DXS
               JDSAO  = JDXAO + 1
               JDXSAO = JDSAO + 1
            ELSE
C              DT, DXT
            END IF
         END IF
C
C        loop over NSRMAT output matrices (later do the loop inside SRDFT for
C        efficiency!!!!!) /hjaaj June 2009
C
         DO JSRMAT = 1, NSR_DMAT
          JDMAT = NDOFF*(JSRMAT-1) + 1
          IF (DOATR) THEN
#ifdef SIRFCK_DEBUG
              write(lupri,*) 'DOATR in sirfck: ndoff =', ndoff
              write(lupri,*) 'JDMAT = ', JDMAT
              write(lupri,*) 'JDXAO = ', JDXAO
              WRITE(LUPRI,*) 'DMAT(1,JDMAT)'
              call output(DMAT(1,JDMAT),1,NBAST,1,NBAST,NBAST,
     &                    NBAST,-1,lupri)
              WRITE(LUPRI,*) 'DMAT(1,JDMAT+1)'
              call output(DMAT(1,JDMAT+1),1,NBAST,1,NBAST,NBAST,
     &                NBAST,-1,lupri)
#endif
C           we need the symmetric part of DXTOTAO in DMAT(1,JDXAO) for SRDFT :
            IF (NDOFF .EQ. 1) THEN
C              closed shell SCF case
               CALL DGETSI(NBAST,DMAT(1,JDMAT),DMAT(1,JDXAO))
            ELSE
               CALL DCOPY(N2BASX,DMAT(1,JDMAT),1,DMAT(1,JDXAO),1)
               CALL DAXPY(N2BASX,D1,DMAT(1,JDMAT+1),1,
     &                              DMAT(1,JDXAO),1)
C              ... DXTOTAO = DXCAO + DXVAO
               CALL DGETSI(NBAST,DMAT(1,JDXAO),DMAT(1,JDXAO))
C              ... extract symmetric component (matrix is assumed
C              symmetric in SRDFT, and it will in general not have any
C              matrix symmetry e.g. when SIRFCK is called from RSPOLI)
               IF (DO_SPINDNS) THEN
                  CALL DGETSI(NBAST,DMAT(1,JDMAT+2),DMAT(1,JDXSAO))
C              ... extract symmetric component of DXSAO
               END IF
               IF (LSRHYBR) THEN
                  CALL DGETSI(NBAST,DMAT(1,JDMAT+1),DMAT(1,JDXAO+1))
C                 ... copy symmetric part of DXVAO
               END IF
            END IF
          ELSE
            JDAO = JDMAT
            IF (NDOFF.NE.1)
     &      CALL DAXPY(N2BASX,D1,DMAT(1,JDMAT+1),1,DMAT(1,JDAO),1)
C           ... get total DMAT = DC + DV if NASHT .gt. 0
          END IF
C         SRDFT(NDMAT,EXCMAT,DMAT,EDFTY,DOERG,DOGRD,DOATR,WORK,LWORK,IPRINT)
#ifdef SIRFCK_DEBUG
          WRITE(LUPRI,*) 'SIRFCK: calling SRDFT, DOERG',DOERG  !JT
          WRITE(LUPRI,*) 'DMAT(1,JDAO), JDAO =',JDAO
          call output(DMAT(1,JDAO),1,NBAST,1,NBAST,NBAST,NBAST,-1,lupri)
          WRITE(LUPRI,*) 'DMAT(1,JDAO+1), JDAO+1 =',JDAO+1
          call output(DMAT(1,JDAO+1),1,NBAST,1,NBAST,NBAST,
     &                NBAST,-1,lupri)
          WRITE(LUPRI,*) 'SIRFCK: DFT_SPINDNS = ', DFT_SPINDNS, 
     &                   'BEFORE SRDFT'
#endif
          CALL DZERO(WRK(KSRAO),NSRAO*N2BASX)
          DFT_SPINDNS_SAVE = DFT_SPINDNS
          DFT_SPINDNS = DO_SPINDNS 
          CALL SRDFT(1,WRK(KSRAO),DMAT(1,JDAO),EDFT,
     &              DOERG,.FALSE.,DOATR,WRK(KLAST),LWRK,IPRSRDFT)
          DFT_SPINDNS = DFT_SPINDNS_SAVE
C         IF (DOATR)
C    &      ESRLTR = ESRLTR + DDOT(N2BASX,WRK(KSRAO),1,DMAT(1,JDAO),1)
C           ... for MCSRDFT csf sigma-vectors : - 2 ESRLTR cref(i)
C               NOT NEEDED because we project out CREF component in sigma vectors
C
          CALL DAXPY(N2BASX,D1,WRK(KSRAO),1,FMAT(1,JDMAT),1)
          IF (DO_SPINDNS)
     &    CALL DAXPY(N2BASX,D1,WRK(KSRAO+N2BASX),1,FMAT(1,JDMAT+2),1)

          IF (LSRHYBR) CALL DCOPY(N2BASX,
     &        WRK(KSRAO+N2BASX),1,FMAT(1,JDMAT+1),1)
C
         END DO  ! JSRMAT = 1, NSR_DMAT
C
C
         IF (DOERG) THEN
            IF (NDOFF.NE.1)
     &         CALL DAXPY(N2BASX,-D1,DMAT(1,2),1,DMAT(1,1),1)
C        ... restore DC = DTOT - DV
            ECSR = -DP5*DDOT(N2BASX,DMAT(1,1),1,WRK(KSRAO),1)
            IF (NASHT.GT.0) THEN
               EVSR = -DDOT(N2BASX,DMAT(1,2),1,WRK(KSRAO),1)
            ELSE
               EVSR = D0
            END IF
            ESRDFTY = EDFT + ECSR + EVSR + EJKVSR
C
C           EDFT = Exsr + Ecsr
C           ECSR = -0.5*Tr(Vxcsr Dc)
C           EVSR = -Tr(Vxcsr Dv)
C           EJKVSR = -0.5*Tr( (Jsr+HFXFAC*Ksr) Dv)
C           note: ECSR + EVSR + EJKVSR corrections are added
C                 to get right contributions from MCSCF code.
C
C           For print, short-range Coulomb and Exchange energies
C           EJSR = 0.5*Tr( (Jsr Dtot)
C           EKSR = 0.5*Tr( (Ksr Dtot) * HFXFAC
C
C           IF (iprsrdft .GT. 0) then
               WRITE (LUPRI,'()')
               WRITE (LUPRI,'(A,F20.10)')
     &         '   Exsr + Ecsr                    ',EDFT,
     &         ' - 0.5 Tr(Vxcsr Dcore)            ',ECSR,
     &         ' -     Tr(Vxcsr Dval)             ',EVSR,
     &         ' - 0.5 Tr( (Jsr+HFXFAC*Ksr) Dval )',EJKVSR,
     &         ' = ESRDFTY                        ',ESRDFTY,
     &         '   EJSR                           ',EJSR,
     &         '   EKSR                           ',EKSR
C           END IF
         END IF ! (DOERG)
      END IF ! (ADDSRI.AND.(DOHFSRDFT.OR.DOMCSRDFT).AND.(.NOT.PREHF))
C
      IF (NODDYDF) THEN
         IF (IFCTYP(1).NE.11) THEN
            WRITE(LUPRI,*)'IFCTYP(1) = ',IFCTYP(1)
            WRITE(LUPRI,*)'DenFit noddy code works only for Coulomb'
            CALL QUIT('Ooops... Probably trying to calculate exchange '
     &                //'with density-fitting noddy code.')
         END IF
         IF (NDMAT.GT.1) THEN
            WRITE(LUPRI,*)'NDMAT = ',NDMAT
            WRITE(LUPRI,*)'DenFit noddy code works only for NDMAT=1'
            CALL QUIT('Illegal NDMAT value in DenFit noddy code')
         END IF
C
         call dzero(FMAT(1,1),N2BASX)
C
         !debugging tool... Simulates real J_ab (no DenFit)
*        call jab_noddy(nbast,dmat(1,1),fmat(1,1))
C
         !calculate density-fitting coefficeints...
         CALL DENFIT_COEFFS_NODDY(NBAST,NBAST,DMAT(1,1),DFCOEFF)
c
         WRITE (LUPRI,'(//,1X,A,I3)') ' Noddy fitting coeffs'
         CALL OUTPUT(DFCOEFF,1,NBAST,1,1,NBAST,1,1,LUPRI)
c
         !...and use them to build the (fitted) Fock matrix
         CALL JAB_DENFIT_NODDY(NBAST,NBAST,DFCOEFF,FMAT(1,1))
#endif  /* MOD_SRDFT */
      END IF
C
#if defined (VAR_OLDCODE)
--- keep print option here ??? 921201-hjaaj
C
C------------------------------------------------------
C     Write out densities and associated fock matrices.
C------------------------------------------------------
C
      IF (IPRFCK.GT.2) THEN
         CALL HEADER('Density and Fock matrices in SIRFCK',-1)
         DO 300 I = 1, NDMAT
            WRITE (LUPRI,'(//A,I3)') ' Density matrix No.',I
            CALL OUTPUT(DMAT(1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            WRITE (LUPRI,'(//A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(1,I),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
  300    CONTINUE
      END IF
#endif
C
Cholesky
      DTIME = SECOND() - DTIME
      TSIRF = TSIRF    + DTIME
      IF (ICOUN3 .EQ. NDMAT) THEN
         WRITE(LUPRI,'(/,A,I3,A,F10.2,A,/)')
     &   '>>> SIRFCK: ',NDMAT,
     &   ' non-symmetric Fock matrices calculated in ',DTIME,' seconds'
      ENDIF
Cholesky
C
 9999 CALL QEXIT('SIRFCK')
      RETURN
C
C
 8000 CONTINUE
      WRITE (LUPRI,'(/A/,(8I10))')
     &   ' FATAL ERROR SIRFCK: Unknown IFCTYP; dump of IFCTYP:',
     &   (IFCTYP(I),I=1,NDMAT)
C
      CALL QTRACE(LUPRI)
      CALL QUIT('SIRFCK ERROR: Unknown IFCTYP')
C
      END
C  /* Deck sirfck2 */
      SUBROUTINE SIRFCK2(LU2INT,FN2INT,
     &                   FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,WRK,LWORK)
C
C     Extracted non-Direct code from SIRFCK into SIRFCK2,
C     with file unit and name in parameter list.
C     This means that this can also be used e.g. for
C     the short range Hartree terms in the MC-DFT model.
C                \jkp+hjaaj
C
C output:
C     FMAT   :  NDMAT Fock matrices with 2-el. cont. added
C  input:
C     FMAT   :  NDMAT Fock matrices
C     DMAT   :  NDMAT density matrices
C     NDMAT  :  number of density matrices and Fock matrices
C     ISYMDM :  symmetry of DMAT and FMAT
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 0 no contribution !
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
#include <implicit.h>
      CHARACTER*(*) FN2INT
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,NDMAT)
      DIMENSION ISYMDM(NDMAT), IFCTYP(NDMAT), MBAS(8)
      DIMENSION WRK(LWORK)
C
#include <iratdef.h>
#include <dummy.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D4 = 4.0D0, DMP5 = -0.5D0)
C
C Used from common blocks
C   INFORB: NSYM, NBAS(:), N2BASX
C   INFPRI: IPRFCK
C
#include <maxorb.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infpri.h>
C
C
#include <dftcom.h>
#include <memint.h>
      CALL QENTER('SIRFCK2')
#ifdef SIRFCK_DEBUG
      write(lupri,*)'SIRFCK2: NDMAT=', NDMAT
      write(lupri,*)'SIRFCK2: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
#endif
C
C     Transfer constant factors to DMAT
C     IFCTYP 11,12,13,22; symmetrization or antisymmetrization:
C          factor 2 because factor 0.5 in FCKEN1 and FCKEN2 in order
C          to preserve previous content of FMAT
C     IFCTYP 11,13: factor  2 (see note in FCKD11/FCKD13)
C     IFCTYP 12,22: factor  0.5 from exchange integral
C     IFCTYP 02   : factor -0.5 from exchange integral
      DO I = 1,NDMAT
         IF (IFCTYP(I) .EQ. 11 .OR. IFCTYP(I) .EQ. 13) THEN
            CALL DSCAL(N2BASX,D4,DMAT(1,I),1)
         ELSE IF (IFCTYP(I) .EQ. 2) THEN
            IF (HFXFAC .EQ. D0) THEN
C              This exchange-only matrix is zero
               IFCTYP(I) = 0
            ELSE
               CALL DSCAL(N2BASX,HFXFAC*DMP5,DMAT(1,I),1)
            END IF
         END IF
      END DO
C
C
      IF (LU2INT .GT. 0) CALL GPCLOSE(LU2INT,'KEEP')
      CALL GPOPEN(LU2INT,FN2INT,'OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
C
      CALL REWSPL(LU2INT)
      CALL MOLLAB('BASINFO ',LU2INT,LUPRI)
      READ (LU2INT) MSYM, MBAS, LBUF, NIBUF, NBITS
      IF (MSYM .NE. NSYM) THEN
         CALL QUIT('SIRFCK2: NSYM from AOxxxINT .ne. NSYM in common')
      END IF
#if !defined (VAR_INT64)
      IF (NBITS .NE. NIBUF*8) THEN
         CALL QUIT('SIRFCK2: NBITS from AOxxxINT .ne. NIBUF*8')
      END IF
#endif
      DO I = 1,MSYM
         IF (NBAS(I) .NE. MBAS(I)) GO TO 110
      END DO
      GO TO 120
  110    WRITE (LUPRI,'(/A/3A,8I5)')
     &      ' SIRFCK2 FATAL ERROR:',
     &      ' NBAS from ',FN2INT,' :',(MBAS(I),I=1,NSYM)
         WRITE (LUPRI,'(A,8I5)')
     &      ' NBAS from /INFORB/ :',(NBAS(I),I=1,NSYM)
         CALL QUIT('SIRFCK2: NBAS(:) from AOxxxINT not consistent')
  120 CONTINUE
C
      LENINT = IRAT*LBUF + NIBUF*LBUF + 1
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('INTE',KINT ,LENINT,WRK,KFREE,LFREE)
      KIINT = KINT + LBUF
      CALL MEMGET('INTE',KIIN4,4*LBUF,WRK,KFREE,LFREE)
C
      CALL MOLLAB('BASTWOEL',LU2INT,LUPRI)
      JBUF = 608
C
    1 CONTINUE
C
#if defined (VAR_SPLITFILES)         
         CALL READSI(LU2INT,LENINT,WRK(KINT),JBUF)
#else
         CALL READI(LU2INT,LENINT,WRK(KINT))
#endif
         CALL AOLAB4(WRK(KIINT),LBUF,NIBUF,NBITS,WRK(KIIN4),NINT)
         IF (NINT .EQ. 0)  GOTO 1
         IF (NINT .EQ. -1) GOTO 2
C
         DO 200 I = 1,NDMAT
            IF (IFCTYP(I) .EQ. 13) THEN
C           ... e.g. singlet, symmetric density matrix
              CALL FCKD13(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 22
     &          .OR. IFCTYP(I) .EQ. 23      
     &          .OR. IFCTYP(I) .EQ. 12) THEN
C           ... Symmetric or antisymmetric, Only Exchange
C               (type 23: Coulomb term is zero because antisymmetric)
C           ... e.g. singlet, antisymmetric density matrix OR
C           ... or   triplet,     symmetric density matrix OR
C           ... or   triplet, antisymmetric density matrix
              CALL FCKD12(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 11) THEN
C           ... Symmetric, Only Coulomb
              CALL FCKD11(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 3) THEN
C           ... Non-symmetric, Coulomb+Exchange
C           ... e.g. singlet, non-symmetric (general) density matrix
              CALL FCKD03(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 2) THEN
C           ... Non-symmetric, Only Exchange
C           ... e.g. triplet, non-symmetric (general) density matrix
              CALL FCKD02(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 1) THEN
C           ... Non-symmetric, Only Coulomb
              CALL FCKD01(FMAT(1,I),DMAT(1,I),WRK(KINT),WRK(KIIN4),NINT)
            ELSE IF (IFCTYP(I) .EQ. 0 .OR. IFCTYP(I) .EQ. 21) THEN
C              00:  nothing
C              21:  antisymmetric, only Coulomb gives a zero matrix
               GOTO 200
            ELSE
              GO TO 8000
            END IF
  200    CONTINUE
C
      GOTO 1
C
    2 CONTINUE
      CALL GPCLOSE(LU2INT,'KEEP')
C
C     Finish symmetric and antisymmetric Fock matrices
C     Restore original density matrices
C
      DO 300 I = 1,NDMAT
         IF (IFCTYP(1) .EQ. 11 .OR. IFCTYP(I) .EQ. 12 .OR.
     &       IFCTYP(I) .EQ. 13) THEN
            CALL FCKEN1(FMAT(1,I),NBAST)
         ELSE IF (IFCTYP(I) .EQ. 22 .OR. IFCTYP(I) .EQ. 23) THEN
            CALL FCKEN2(FMAT(1,I),NBAST)
         END IF
C
         IF (IFCTYP(I) .EQ. 11 .OR. IFCTYP(I) .EQ. 13) THEN
            DFAC = D1 / D4
            CALL DSCAL(N2BASX,DFAC,DMAT(1,I),1)
         ELSE IF (IFCTYP(I) .EQ. 2) THEN
            IF (HFXFAC .NE. D0) THEN
               DFAC = D1 / (DMP5*HFXFAC)
               CALL DSCAL(N2BASX,DFAC,DMAT(1,I),1)
            END IF
         END IF
C
  300 CONTINUE
C
      CALL QEXIT('SIRFCK2')
      RETURN
C
 8000 CONTINUE
      WRITE (LUPRI,'(/A/,(8I10))')
     &   ' FATAL ERROR SIRFCK2: Unknown IFCTYP; dump of IFCTYP:',
     &   (IFCTYP(I),I=1,NDMAT)
      CALL QTRACE(LUPRI)
      CALL QUIT('SIRFCK2 ERROR: Unknown IFCTYP')
      END
C  /* Deck fckd13 */
      SUBROUTINE FCKD13(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C Hans Joergen Aa. Jensen October 1994: general routine for
C symmetric singlet Fock matrix from symmetric singlet density matrix
C (based on FCKD03)
C
C This subroutine adds (derivative) two-electron integrals to
C Fock matrices. The Fock matrices are assumed
C to be square matrices in full dimension without symmetry reduction
C in size. Remember to zero out the Fock matrices before starting
C to accumulate.
C
C tsaue Feb 1997: indices permuted for exchange
C HJAaJ October 1994: NOTE:
C     Symmetry is implicitly taken into account by the fact that
C     only non-zero integrals are considered.  The symmetry of DMAT
C     is not taken into account.
C DFT modifications T. Helgaker
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (DP25 = 0.25D0)
      INTEGER P, Q, R, S
#include <inforb.h>
#include <dftcom.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      EFAC = DP25*HFXFAC
C     ... 0.25 because DMAT is multiplied by extra 2 in SIRFCK
C         in order to get full Coulomb contribution:
C         DINT*(DMAT(R,S) + DMAT(S,R)) = 2*DINT*DMAT(R,S)
C         Compare general code in FCKD03.
      DO INT = 1, LENGTH
C        factor "4" is multiplied on DMAT in SIRFCK
C        DINT = 4*BUF(INT)
         DINT =   BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DINT / 2
         IF (R.EQ.S)              DINT = DINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DINT / 2
         EINT = EFAC*DINT
         FMAT(P,Q) = FMAT(P,Q) + DINT*DMAT(R,S)
         FMAT(P,R) = FMAT(P,R) - EINT*DMAT(Q,S)
         FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(P,S)
         FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(P,R)
         FMAT(P,S) = FMAT(P,S) - EINT*DMAT(Q,R)
         FMAT(R,S) = FMAT(R,S) + DINT*DMAT(P,Q)
C        ... DMAT(i,j) = DMAT(j,i) has been used above
      END DO
      RETURN
      END
C  /* Deck fckd12 */
      SUBROUTINE FCKD12(FMAT,XMAT,BUF,IBUF,LENGTH)
C
C     Hans Joergen Aa. Jensen and Rika Kobayashi 24-Apr-1992
C     (based on FCKDIR by Henrik Koch and Trygve Helgaker)
C     970304-tsaue: index permuted
C     941011-hjaaj: renamed from LFKDIR to FCKD12
C                   more than halved the work by using
C                   the antisymmetry also in the exchange part
C     DFT modifications T. Helgaker
C
C     This subroutine adds two-electron integrals multiplied
C     by antisymmetric effective density matrices to Fock matrices
C     (e.g. from London reorthonormalization one-index transformations).
C     The Fock matrices are assumed to be square matrices
C     in full dimension without symmetry reduction in size.
C     Remember to zero out the fock matrices before starting
C     to accumulate.
C
C     The matrix is finished by antisymmetrization in FCKES2.
C
C     FMAT(P,R)  = - 0.5 SUM(q,s) (p q | r s) X(s,q)
C                = + 0.5 SUM(q,s) (p q | r s) X(q,s)
C     Note that X(q,s) = -X(s,q)
C
#include <implicit.h>
#include <priunit.h>
      INTEGER P, Q, R, S
#include <inforb.h>
#include <dftcom.h>
      DIMENSION FMAT(NBAST,NBAST), XMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      EFAC = HFXFAC
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = EFAC*BUF(INT)
         IF (P.EQ.Q)              GINT = GINT / 2
         IF (R.EQ.S)              GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
         FMAT(P,R) = FMAT(P,R) - GINT*XMAT(Q,S)
         FMAT(P,S) = FMAT(P,S) - GINT*XMAT(Q,R)
         FMAT(Q,R) = FMAT(Q,R) - GINT*XMAT(P,S)
         FMAT(Q,S) = FMAT(Q,S) - GINT*XMAT(P,R)
CFCKD03 MAERK 941013: bemaerk 0.5 faktor mangler sammenlignet med FCKD03
CFCKD03               and XMAT transposed using XMAT(i,j) = -XMAT(j,i)
CFCKD03  FMAT(P,R) = FMAT(P,R) - EINT*DMAT(S,Q)
CFCKD03  FMAT(P,S) = FMAT(P,S) - EINT*DMAT(R,Q)
CFCKD03  FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(S,P)
CFCKD03  FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckd11 */
      SUBROUTINE FCKD11(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C Hans Joergen Aa. Jensen March 2003: add to FMAT for
C IX = 1 (symmetric DMAT), IY = 1 (only Coulomb)
C (based on FCKD13)
C
C The Fock matrices are assumed to be square matrices
C in full dimension without symmetry reduction in size.
C Remember to zero out the fock matrices before starting
C to accumulate.
C
C HJAaJ October 1994: NOTE:
C     Symmetry is implicitly taken into account by the fact that
C     only non-zero integrals are considered.  The symmetry of DMAT
C     is not taken into account.
C
#include <implicit.h>
#include <priunit.h>
      INTEGER P, Q, R, S
#include <inforb.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      DO INT = 1, LENGTH
C        factor "4" is multiplied on DMAT in SIRFCK
C        DINT = 4*BUF(INT)
         DINT =   BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              DINT = DINT / 2
         IF (R.EQ.S)              DINT = DINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) DINT = DINT / 2
         FMAT(P,Q) = FMAT(P,Q) + DINT*DMAT(R,S)
         FMAT(R,S) = FMAT(R,S) + DINT*DMAT(P,Q)
      END DO
      RETURN
      END
C  /* Deck fckd03 */
      SUBROUTINE FCKD03(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     Henrik Koch and Trygve Helgaker 18-NOV-1991.
C     970303-tsaue : index permutation
C     941011-hjaaj: renamed from FCKDI1 to FCKD03
C     DFT modifications T. Helgaker
C
C     This subroutine adds derivative two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include <implicit.h>
#include <priunit.h>
      INTEGER P, Q, R, S
#include <inforb.h>
#include <dftcom.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
#include <ibtfun.h>
C
      EFAC = HFXFAC / 2
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = BUF(INT)
         IF (P.EQ.Q)              GINT = GINT / 2
         IF (R.EQ.S)              GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
C           coulomb:
         FADD = GINT*(DMAT(R,S)+DMAT(S,R))
         FMAT(P,Q) = FMAT(P,Q) + FADD
         FMAT(Q,P) = FMAT(Q,P) + FADD
         FADD = GINT*(DMAT(P,Q)+DMAT(Q,P))
         FMAT(R,S) = FMAT(R,S) + FADD
         FMAT(S,R) = FMAT(S,R) + FADD
C           exchange:
         GINT = EFAC*GINT
         FMAT(P,R) = FMAT(P,R) - GINT*DMAT(Q,S)
         FMAT(R,P) = FMAT(R,P) - GINT*DMAT(S,Q)
         FMAT(P,S) = FMAT(P,S) - GINT*DMAT(Q,R)
         FMAT(S,P) = FMAT(S,P) - GINT*DMAT(R,Q)
         FMAT(Q,R) = FMAT(Q,R) - GINT*DMAT(P,S)
         FMAT(R,Q) = FMAT(R,Q) - GINT*DMAT(S,P)
         FMAT(Q,S) = FMAT(Q,S) - GINT*DMAT(P,R)
         FMAT(S,Q) = FMAT(S,Q) - GINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckd02 */
      SUBROUTINE FCKD02(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     Hans Joergen Aa. Jensen 11-Oct-1994
C     General routine for triplet Fock matrices, based on FCKD03.
C     970304 - tsaue: Index permutation
C
C     This subroutine adds two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include <implicit.h>
#include <priunit.h>
      INTEGER P, Q, R, S
#include <inforb.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
#include <ibtfun.h>
C
      DO 100 INT = 1, LENGTH
C        factor "-0.5" is multiplied on DMAT in SIRFCK
C        factor HFXFAC is multiplied on DMAT in SIRFCK
C        GINT = - ( HFXFAC * BUF(INT) ) / 2
         GINT = BUF(INT)
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         IF (P.EQ.Q)              GINT = GINT / 2
         IF (R.EQ.S)              GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
         FMAT(P,R) = FMAT(P,R) + GINT*DMAT(Q,S)
         FMAT(Q,R) = FMAT(Q,R) + GINT*DMAT(P,S)
         FMAT(R,P) = FMAT(R,P) + GINT*DMAT(S,Q)
         FMAT(S,P) = FMAT(S,P) + GINT*DMAT(R,Q)
         FMAT(P,S) = FMAT(P,S) + GINT*DMAT(Q,R)
         FMAT(Q,S) = FMAT(Q,S) + GINT*DMAT(P,R)
         FMAT(R,Q) = FMAT(R,Q) + GINT*DMAT(S,P)
         FMAT(S,Q) = FMAT(S,Q) + GINT*DMAT(R,P)
  100 CONTINUE
      RETURN
      END
C  /* Deck fckd01 */
      SUBROUTINE FCKD01(FMAT,DMAT,BUF,IBUF,LENGTH)
C
C     HJAaJ Mar 2003: based on FCKD03
C
C     The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include <implicit.h>
      INTEGER P, Q, R, S
#include <inforb.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LENGTH), IBUF(4,LENGTH)
C
      DO 100 INT = 1, LENGTH
         P    = IBUF(1,INT)
         Q    = IBUF(2,INT)
         R    = IBUF(3,INT)
         S    = IBUF(4,INT)
         GINT = BUF(INT)
         IF (P.EQ.Q)              GINT = GINT / 2
         IF (R.EQ.S)              GINT = GINT / 2
         IF (P.EQ.R .AND. S.EQ.Q) GINT = GINT / 2
C           coulomb:
         FADD = GINT*(DMAT(R,S)+DMAT(S,R))
         FMAT(P,Q) = FMAT(P,Q) + FADD
         FMAT(Q,P) = FMAT(Q,P) + FADD
         FADD = GINT*(DMAT(P,Q)+DMAT(Q,P))
         FMAT(R,S) = FMAT(R,S) + FADD
         FMAT(S,R) = FMAT(S,R) + FADD
  100 CONTINUE
      RETURN
      END
C  /* Deck fcken1 */
      SUBROUTINE FCKEN1(FMAT,NBAST)
C
C     Hans Joergen Aa. Jensen October 1994.
C
C     This subroutine completes symmetric Fock matrices
C     accumulated by FCKD13 or FCKD02 by symmetrization.
C
C     Note: 0.5 factor is required in order to preserve
C     symmetric content of FMAT from before SIRFCK call.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION FMAT(NBAST,NBAST)
      INTEGER   P, Q
C
      DO P = 1,NBAST
         DO Q = 1,P
            FPQ = ( FMAT(P,Q) + FMAT(Q,P) ) / 2
            FMAT(P,Q) = FPQ
            FMAT(Q,P) = FPQ
         END DO
      END DO
      RETURN
      END
C  /* Deck fcken2 */
      SUBROUTINE FCKEN2(FMAT,NBAST)
C
C     Hans Joergen Aa. Jensen October 1994.
C
C     This subroutine completes antisymmetric Fock matrices
C     accumulated by FCKD12 by antisymmetrization.
C
C     Note: 0.5 factor is required in order to preserve
C     symmetric content of FMAT from before SIRFCK call.
C
#include <implicit.h>
#include <priunit.h>
      DIMENSION FMAT(NBAST,NBAST)
      INTEGER   P, Q
C
      DO P = 1,NBAST
         DO Q = 1,P
            FPQ = ( FMAT(P,Q) - FMAT(Q,P) ) / 2
            FMAT(P,Q) =  FPQ
            FMAT(Q,P) = -FPQ
         END DO
      END DO
      RETURN
      END
C
C  /* Deck cho_fckds1 */
      SUBROUTINE CHO_FCKDS1(FOCK,DEN,WORK,LWORK,ISYDEN,IOPTVP)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Calculate two-electron integral contribution to
C        FOCK matrix using Cholesky decomposition of
C        (1) two-electron integrals, AND
C        if CMO not available on disk,
C        (2) density matrix DEN (decomposed locally).
C
C     Notes:
C        (a) The density matrix *must* be tot. sym. (ISYDEN = 1).
C        (b) FOCK and DEN are assumed stored as full square
C            without symmetry reduction.
C
#include <implicit.h>
      DIMENSION FOCK(*), DEN(*), WORK(LWORK)
#include <infpri.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <choscf.h>
#include <priunit.h>
#include <ccdeco.h>

      INTEGER NVECD(8), IVECD(8,8)

      PARAMETER (HALF = 0.50D0, ONE = 1.00D0)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_FCKDS1')

      LOGICAL LOCDBG, TINFO
      PARAMETER (LOCDBG = .FALSE., TINFO = .FALSE.)

      LOGICAL IOPTVP
C
C
C     Start timing.
C     -------------

      TIMT = SECOND()

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Entry')
         CALL OUTPUT(DEN,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL AROUND(SECNAM//': Fock Matrix on Entry')
         CALL OUTPUT(FOCK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF

C     Check ISYDEN.
C     -------------

      IF (ISYDEN .NE. 1) THEN
         WRITE(LUPRI,'(//,A,A)')
     &   SECNAM,': FATAL ERROR: Density *must* be tot. symmetric !'
         WRITE(LUPRI,'(A,I10,/)')
     &   'Symmetry of DEN on input: ',ISYDEN
         CALL QUIT('SYMMETRY MISMATCH IN '//SECNAM)
      ENDIF

C     Allocation.
C     -----------

      KFOCK = 1
      KDEN  = KFOCK + N2BST(ISYDEN)
      KEND1 = KDEN  + N2BST(ISYDEN)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,8X,A)')
     &   'Insufficient memory in ',SECNAM,':',
     &   'Allocation: Symmetry pack'
         WRITE(LUPRI,'(8X,A,I10,/,8X,A,I10,/)')
     &   'Need      : ',KEND1-1,
     &   'Available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Symmetry pack DEN and scale.
C     ----------------------------

      TIMPK = SECOND()
      CALL CHO_DENPK(DEN,WORK(KDEN),ISYDEN)
      CALL DSCAL(N2BST(ISYDEN),HALF,WORK(KDEN),1)
      TIMPK = SECOND() - TIMPK

C     Cholesky decompose symmetry blocks of DEN if necessary;
C     else read occupied part of CMO from disk.
C     Set up index array for the CMO 'vectors'.
C     -------------------------------------------------------

      TIMDC = SECOND()

      IF (NEWSCF) CCMODSK = .FALSE.
      IF (.NOT. CCMODSK) THEN

         CALL DCOPY(N2BST(ISYDEN),WORK(KDEN),1,WORK(KFOCK),1)

         KVECD  = KEND1
         ICOUNT = 0

         DO ISYM = 1,NSYM

            IVECD(ISYM,ISYM) = ICOUNT

            KEND2 = KVECD + IVECD(ISYM,ISYM)
            LWRK2 = LWORK - KEND2 + 1
            IF (LWRK2 .LT. 0) THEN
               WRITE(LUPRI,'(//,5X,A,A,A,/,8X,A,I1)')
     &         'Insufficient memory in ',SECNAM,':',
     &         'Allocation: Density decomposition, sym. block: ',ISYM
               WRITE(LUPRI,'(8X,A,I10,/,8X,A,I10,/)')
     &         'Need      : ',KEND2-1,
     &         'Available : ',LWORK
               CALL QUIT('Insufficient memory in '//SECNAM)
            ENDIF

            KOFFD = KFOCK + IAODIS(ISYM,ISYM)
            CALL SIR_CHODECO(WORK(KOFFD),NBAS(ISYM),WORK(KEND2),LWRK2,
     &                       THRDC,NVECD(ISYM))

            ICOUNT = ICOUNT + NBAS(ISYM)*NVECD(ISYM)

         ENDDO

         KEND2 = KVECD + ICOUNT
         LWRK2 = LWORK - KEND2

      ELSE

         ICOUNT = 0
         DO ISYM = 1,NSYM
            IVECD(ISYM,ISYM) = ICOUNT
            NVECD(ISYM) = NDVCS(ISYM)
            ICOUNT = ICOUNT + NBAS(ISYM)*NVECD(ISYM)
         ENDDO
         LEN  = ICOUNT

         KVECD = KEND1
         KEND2 = KVECD + LEN
         LWRK2 = LWORK - KEND2 + 1

         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A,/,8X,A)')
     &      'Insufficient memory in ',SECNAM,':',
     &      'Allocation: occ. CMO read'
            WRITE(LUPRI,'(8X,A,I10,/,8X,A,I10,/)')
     &      'Need      : ',KEND2-1,
     &      'Available : ',LWORK
            CALL QUIT('Insufficient memory in '//SECNAM)
         ENDIF

         CALL WOPEN2(LUCCMO,FCCMO,64,0)
         IADR = 1
         CALL GETWA2(LUCCMO,FCCMO,WORK(KVECD),IADR,LEN)
         CALL WCLOSE2(LUCCMO,FCCMO,'KEEP')

      ENDIF

      TIMDC = SECOND() - TIMDC

C     Debug: Print number of vectors.
C     Debug: Check decomposition/CMO read.
C     ------------------------------------

      IF (LOCDBG) THEN
         IF (.NOT. CCMODSK) THEN
            WRITE(LUPRI,'(/,2X,A,A,/)')
     &      SECNAM,': Cholesky decomposition of density:'
         ELSE
            WRITE(LUPRI,'(/,2X,A,A,/)')
     &      SECNAM,': occ. CMO read:'
         ENDIF
         WRITE(LUPRI,'(2X,A,/,2X,A)')
     &   'Sym.  No. vectors  Offset',
     &   '-------------------------'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(3X,I1,4X,I6,7X,I6)')
     &      ISYM,NVECD(ISYM),IVECD(ISYM,ISYM)
         ENDDO
         WRITE(LUPRI,'(2X,A,/)')
     &   '-------------------------'
      ENDIF

C     Calculate contributions.
C     ------------------------

      TIMCL = SECOND()
      CALL DZERO(WORK(KFOCK),N2BST(ISYDEN))
      CALL SCF_CHOFCK1(WORK(KDEN),WORK(KVECD),WORK(KFOCK),
     &                 WORK(KEND2),LWRK2,NVECD,IVECD,TINFO,IOPTVP)
      TIMCL = SECOND() - TIMCL

C     Add result to FOCK matrix.
C     --------------------------

      TIMAD = SECOND()
      DO ISYMB = 1,NSYM
         ISYMA = MULD2H(ISYMB,ISYDEN)
         DO IB = 1,NBAS(ISYMB)
            B = IBAS(ISYMB) + IB
            KOFF1 = KFOCK + IAODIS(ISYMA,ISYMB)
     &            + NBAS(ISYMA)*(IB - 1)
            KOFF2 = NBAST*(B - 1) + IBAS(ISYMA) + 1
            CALL DAXPY(NBAS(ISYMA),ONE,WORK(KOFF1),1,FOCK(KOFF2),1)
         ENDDO
      ENDDO
      TIMAD = SECOND() - TIMAD

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Exit')
         CALL OUTPUT(DEN,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL AROUND(SECNAM//': Fock Matrix on Exit')
         CALL OUTPUT(FOCK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF

C     Print timing.
C     -------------

      IF (TINFO .OR. LOCDBG) THEN
         TIMT = SECOND() - TIMT
         WRITE(LUPRI,'(/,5X,A,A,A)')
     &   ' - exiting ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sym. packing density : ',TIMPK,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for decomposing density  : ',TIMDC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating contr.   : ',TIMCL,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for adding to Fock matr. : ',TIMAD,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                         : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck scf_chofck1 */
      SUBROUTINE SCF_CHOFCK1(DEN,CMO,FOCK,WORK,LWORK,NVECD,IVECD,
     &                       TINFO,IOPTVP)
C
C     Written by Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Calculate two-electron integral contributions to the AO Fock
C        matrix based on Cholesky decomposed integrals.
C        Result is added into Fock matrix which is assumed allocated
C        as symmetry packed (!) full square.
C        The density matrix is used for the Coulomb part
C        and the MO coefficient matrix (or the Cholesky vector
C        representing the density matrix) for the exchange part.
C
C     Version 1:
C        The density and the MO coefficients are assumed
C        TOTAL SYMMETRIC, as appropriate for an SCF ground state
C        closed-shell optimization.
C
C     NB: NO SPARSITY IS USED WHATSOEVER !!!!!
C
#include <implicit.h>
      DIMENSION DEN(*), CMO(*), FOCK(*), WORK(LWORK)
      INTEGER IVECD(8,8), NVECD(8)
      LOGICAL TINFO
#include <iratdef.h>
#include <infpri.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <sirftim.h>
#include <priunit.h>
Casm
#include <dftcom.h>
#include <choscf.h>
Casm

      INTEGER IHOLE(8,8), NACHH(8), IALPHK(8,8), NALPHK(8)

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (ISTORE = 2)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'SCF_CHOFCK1')

      LOGICAL IOPTVP


C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMO = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIMX = ZERO

C     Read reduce index array.
C     ------------------------

      KIND1 = 1
      CALL CC_GETIND1(WORK(KIND1),LWORK,LIND1)
      KEND0 = KIND1 + LIND1
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND0-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Entry')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IAODIS(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            CALL OUTPUT(DEN(KOFF),1,NUMB,1,NUMB,NUMB,NUMB,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': MO Coefficient Matrix on Entry')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IVECD(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            NUMK = NVECD(ISYM)
            CALL OUTPUT(CMO(KOFF),1,NUMB,1,NUMK,NUMB,NUMK,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': Fock Matrix on Entry')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IAODIS(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            CALL OUTPUT(FOCK(KOFF),1,NUMB,1,NUMB,NUMB,NUMB,1,LUPRI)
         ENDDO
      ENDIF

C     Set up local index arrays.
C     --------------------------

      DO ISYCHO = 1,NSYM

         ICOUNT = 0

         DO ISYMK = 1,NSYM

            ISYMA = MULD2H(ISYMK,ISYCHO)

            IALPHK(ISYMA,ISYMK) = ICOUNT

            ICOUNT = ICOUNT + NBAS(ISYMA)*NVECD(ISYMK)

         ENDDO

         NALPHK(ISYCHO) = ICOUNT

      ENDDO

C     Start of Cholesky loop.
C     -----------------------

      DO ISYCHO = 1,NSYM

         IF (N2BST(ISYCHO)  .LE. 0)  GO TO 999
         IF (NUMCHO(ISYCHO) .LE. 0)  GO TO 999
         IF (NALPHK(ISYCHO) .LE. 0)  GO TO 999

C        Allocation: Cholesky.
C        ---------------------

         LREAD = MEMRD(1,ISYCHO,ISTORE)
         LSCR1 = MAX(LREAD,NALPHK(ISYCHO))

         KCHOL = KEND0
         KSCR1 = KCHOL + N2BST(ISYCHO)
         KEND1 = KSCR1 + LSCR1
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - allocation: Cholesky'
            WRITE(LUPRI,'(5X,A,I10,/,A,I10,/)')
     &      'Need     : ',KEND1,
     &      'Available: ',LWORK
            CALL QUIT('INSUFFICIENT MEMORY IN '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MINMEM = NALPHK(ISYCHO)
         NVEC   = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - Allocation: Batch'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Need per vector    : ',MINMEM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available for batch: ',LWRK1
            CALL QUIT('INSUFFICIENT MEMORY IN '//SECNAM)
         ENDIF

         KHOLE = KEND1
         KEND2 = KHOLE + NVEC*NALPHK(ISYCHO)
         LWRK2 = LWORK - KEND2 + 1

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

         IF (LOCDBG) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'BATCH INFO FROM ',SECNAM,':'
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Cholesky symmetry                : ',ISYCHO
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Total number of vectors          : ',NUMCHO(ISYCHO)
            WRITE(LUPRI,'(8X,A,I10)')
     &      'AO vector dimension              : ',N2BST(ISYCHO)
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Total memory                     : ',LWORK
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Memory used for AO vectors       : ',KEND1-1
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Memory left for (alpha,k) vectors: ',LWRK1
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Min. memory for (alpha,k) vectors: ',MINMEM
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Number of batches                : ',NBATCH
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Number of vectors per batch      : ',NVEC
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Memory remainder                 : ',LWRK2
         ENDIF

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

C           Set up index arrays for HOLE.
C           -----------------------------

            ICOUNH = 0
            DO ISYMK = 1,NSYM

               ISYMA = MULD2H(ISYMK,ISYCHO)

               IHOLE(ISYMA,ISYMK) = ICOUNH

               NACHH(ISYMA) = NBAS(ISYMA)*NUMV

               ICOUNH = ICOUNH + NACHH(ISYMA)*NVECD(ISYMK)

            ENDDO

            JCHOL1 = NVEC*(IBATCH - 1) + 1

            DO JVEC = 1,NUMV

               JCHOL = JCHOL1 + JVEC - 1

C              Read vector.
C              ------------

               DTIME = SECOND()
               CALL CHO_READN(WORK(KCHOL),JCHOL,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,ISTORE,WORK(KSCR1),LSCR1)
c              xnorm = ddot(nnbst(isycho),work(kchol),1,work(kchol),1)
c              write(lupri,'(a,i4,a,i2,a,d12.4)') 'Vector',jchol,
c    &               ' of symmetry',isycho,'. Norm: ',xnorm 
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

C              Coulomb part.
C              -------------

               IF ((ISYCHO .EQ. 1) .AND. IOPTVP) THEN
                  DTIME = SECOND()
                  CONST = TWO*DDOT(N2BST(ISYCHO),DEN,1,WORK(KCHOL),1)
                  CALL DAXPY(N2BST(ISYCHO),CONST,WORK(KCHOL),1,
     &                       FOCK,1)
                  DTIME = SECOND() - DTIME
                  TIMC  = TIMC     + DTIME
               ENDIF

Casm
C              No exchange for LDA. Thus, no transformation
C
               IF (ITSLDA) GOTO 788
Casm

C              Hole transformation.
C              --------------------

               DTIME = SECOND()
               DO ISYMD = 1,NSYM

                  ISYMA = MULD2H(ISYMD,ISYCHO)
                  ISYMK = ISYMD

                  NA    = NBAS(ISYMA)
                  NTOTA = MAX(NA,1)
                  ND    = NBAS(ISYMD)
                  NTOTD = MAX(ND,1)
                  NK    = NVECD(ISYMK)

                  KOFFC = KCHOL + IAODIS(ISYMD,ISYMA)
                  KOFFH = IVECD(ISYMD,ISYMK) + 1
                  KOFFS = KSCR1 + IALPHK(ISYMA,ISYMK)

                  CALL DGEMM('T','N',NA,NK,ND,ONE,
     &                       WORK(KOFFC),NTOTD,CMO(KOFFH),NTOTD,
     &                       ZERO,WORK(KOFFS),NTOTA)

               ENDDO
               DTIME = SECOND() - DTIME
               TIMO  = TIMO     + DTIME

C              Copy into HOLE array.
C              ---------------------

               DTIME = SECOND()
               DO ISYMK = 1,NSYM

                  ISYMA = MULD2H(ISYMK,ISYCHO)

                  DO K = 1,NVECD(ISYMK)

                     KOFFS = KSCR1 + IALPHK(ISYMA,ISYMK)
     &                     + NBAS(ISYMA)*(K - 1)
                     KOFFH = KHOLE + IHOLE(ISYMA,ISYMK)
     &                     + NACHH(ISYMA)*(K - 1)
     &                     + NBAS(ISYMA)*(JVEC - 1)

                     CALL DCOPY(NBAS(ISYMA),WORK(KOFFS),1,
     &                          WORK(KOFFH),1)

                  ENDDO

               ENDDO

  788          CONTINUE

               DTIME = SECOND() - DTIME
               TIMS  = TIMS     + DTIME
               
            ENDDO
Casm
C           Do not include exchange for DFT/LDA
C
            IF (ITSLDA) GOTO 789
Casm

C           Exchange part.
C           --------------

            DTIME = SECOND()
            DO ISYMK = 1,NSYM

               ISYMA = MULD2H(ISYMK,ISYCHO)
               ISYMB = ISYMA

               NA    = NBAS(ISYMA)
               NTOTA = MAX(NA,1)
               NB    = NBAS(ISYMB)
               NTOTB = MAX(NB,1)
               NVK   = NUMV*NVECD(ISYMK)

               KOFFF = IAODIS(ISYMA,ISYMB) + 1
               KOFFH = KHOLE + IHOLE(ISYMA,ISYMK)
Casm
C              Very dirty hack to solve weird bug when linking
C              against some ATLAS libraries.
C              See the end of cc_cholesky.F for more info
Casm
               CALL OTHER_DGEMM('N','T',NA,NB,NVK,
     &                    HFXFAC*XMONE,WORK(KOFFH),NTOTA,WORK(KOFFH),
     &                    NTOTB,ONE,FOCK(KOFFF),NTOTA)

            ENDDO
            DTIME = SECOND() - DTIME
            TIMX  = TIMX     + DTIME
 
  789       CONTINUE

         ENDDO

  999    CONTINUE

      ENDDO

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Exit')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IAODIS(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            CALL OUTPUT(DEN(KOFF),1,NUMB,1,NUMB,NUMB,NUMB,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': MO Coefficient Matrix on Exit')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IVECD(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            NUMK = NVECD(ISYM)
            CALL OUTPUT(CMO(KOFF),1,NUMB,1,NUMK,NUMB,NUMK,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': Fock Matrix on Exit')
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(/,20X,A,I1)')
     &      'Symmetry Block: ',ISYM
            KOFF = IAODIS(ISYM,ISYM) + 1
            NUMB = NBAS(ISYM)
            CALL OUTPUT(FOCK(KOFF),1,NUMB,1,NUMB,NUMB,NUMB,1,LUPRI)
         ENDDO
      ENDIF

C     Print timing.
C     -------------

      CSIRF = CSIRF + TIMC
      XSIRF = XSIRF + TIMX
      RSIRF = RSIRF + TIMR
      OSIRF = OSIRF + TIMO
      SSIRF = SSIRF + TIMS

      IF (TINFO .OR. LOCDBG) THEN
c      IF (.true.) THEN
         TIMT = SECOND() - TIMT
         WRITE(LUPRI,'(/,5X,A,A,A)')
     &   ' - exiting ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O and squaring     : ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for Coulomb  part        : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for partial MO transform.: ',TIMO,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sorting (copy)       : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for exchange part        : ',TIMX,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                         : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck cho_denpk */
      SUBROUTINE CHO_DENPK(DEN,DENPK,ISYDEN)
C
C     Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Symmetry pack DEN.
C
#include <implicit.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>

      DIMENSION DEN(*), DENPK(*)

      DO ISYMB = 1,NSYM
         ISYMA = MULD2H(ISYMB,ISYDEN)
         DO IB = 1,NBAS(ISYMB)
            B = IBAS(ISYMB) + IB
            DO IA = 1,NBAS(ISYMA)
               A = IBAS(ISYMA) + IA
               KABPK = IAODIS(ISYMA,ISYMB) + NBAS(ISYMA)*(IB - 1) + IA
               KABFL = NBAST*(B - 1) + A
               DENPK(KABPK) = DEN(KABFL)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END
C  /* Deck sir_chodeco */
      SUBROUTINE SIR_CHODECO(XMAT,NDIM,WORK,LWORK,THRS,NVECS)
C
C     Henrik Koch   10-may-2002
C
#include <implicit.h>
#include <infpri.h>
#include <priunit.h>
C
      DIMENSION XMAT(NDIM,NDIM),WORK(NDIM,*)
C
      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'SIR_CHODECO')
      PARAMETER (TOL = 5.00D-13)
C
c     DO I = 1,NDIM
c        WRITE(LUPRI,*) 'Diagonal :',I,XMAT(I,I)
c     ENDDO
C
      IF (LWORK .LT. NDIM*NDIM) STOP 'Not enough memory'
C
      DO I = 1,NDIM
C
         XMAX = -1.0D10
         DO J = 1,NDIM
            IF (ABS(XMAT(J,J)) .GT. XMAX) THEN
               XMAX = XMAT(J,J)
               IMAX = J
            ENDIF
         ENDDO
C 
         IF (XMAX .LT. 0.0D0) THEN
            WRITE(LUPRI,'(A,A,A)')
     &      '*** ',SECNAM,': NOTICE ***'
            WRITE(LUPRI,'(A,1P,D30.15)')
     &      '    Negative diagonal encountered: ',XMAX
            IF (ABS(XMAX) .GT. TOL) THEN
               WRITE(LUPRI,'(A,/)')
     &         '    - unable to continue!'
               CALL QUIT('ERROR IN '//SECNAM)
            ELSE
               WRITE(LUPRI,'(A)')
     &         '    - program continues nevertheless....'
            ENDIF
         ENDIF
C
         IF (ABS(XMAX) .LT. THRS) THEN
            NVECS= I-1
            RETURN
         ENDIF
C
c        WRITE(LUPRI,*)  'Largest diagonal  : ',NDIM,I,IMAX,XMAX
C
         DO J = 1,NDIM
            WORK(J,I) = XMAT(J,IMAX)/SQRT(XMAX)
         ENDDO
C
         DO J = 1,NDIM
            DO K = 1,NDIM
               XMAT(K,J) = XMAT(K,J) - WORK(K,I)*WORK(J,I)
            ENDDO
         ENDDO
C
         DO J = 1,NDIM
            XMAT(J,IMAX) = 0.0D0
            XMAT(IMAX,J) = 0.0D0
         ENDDO
C
      ENDDO
C
      RETURN
      END
C  /* Deck cho_fckds3 */
      SUBROUTINE CHO_FCKDS3(FOCK,DEN,WORK,LWORK,ISYDEN,IOPTVP)
C
C     Thomas Bondo Pedersen, July 2002.
C
C     Purpose:
C        Calculate two-electron integral contribution to
C        FOCK matrix using Cholesky decomposition of
C        the two-electron integrals.
C
C     Notes:
C        (a) FOCK and DEN are assumed stored as full square
C            without symmetry reduction.
C        (b) This is a pilot version: density driven but
C            with use of neither sparsity of the Cholesky vectors
C            nor screening.
C
#include <implicit.h>
      DIMENSION FOCK(*), DEN(*), WORK(LWORK)
#include <infpri.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccisao.h>
#include <ccsdsym.h>
#include <priunit.h>

      PARAMETER (HALF = 0.50D0, ONE = 1.00D0)

      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'CHO_FCKDS3')

      LOGICAL LOCDBG, TINFO, IOPTVP
      PARAMETER (LOCDBG = .FALSE., TINFO = .FALSE.)

C     Start timing.
C     -------------

      TIMT = SECOND()

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Entry')
         CALL OUTPUT(DEN,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL AROUND(SECNAM//': Fock Matrix on Entry')
         CALL OUTPUT(FOCK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF

C     Allocation.
C     -----------

      KFOCK = 1
      KDEN  = KFOCK + N2BST(ISYDEN)
      KEND1 = KDEN  + N2BST(ISYDEN)
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A,/,8X,A)')
     &   'Insufficient memory in ',SECNAM,':',
     &   'Allocation: Symmetry pack'
         WRITE(LUPRI,'(8X,A,I10,/,8X,A,I10,/)')
     &   'Need      : ',KEND1-1,
     &   'Available : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Symmetry pack DEN and scale.
C     ----------------------------

      TIMPK = SECOND()
      CALL CHO_DENPK(DEN,WORK(KDEN),ISYDEN)
      CALL DSCAL(N2BST(ISYDEN),HALF,WORK(KDEN),1)
      TIMPK = SECOND() - TIMPK

C     Calculate contributions.
C     ------------------------

      TIMCL = SECOND()
      CALL DZERO(WORK(KFOCK),N2BST(ISYDEN))
      CALL SCF_CHOFCK3(WORK(KDEN),WORK(KFOCK),WORK(KEND1),LWRK1,
     &                 ISYDEN,TINFO,IOPTVP)
      TIMCL = SECOND() - TIMCL

C     Add result to FOCK matrix.
C     --------------------------

      TIMAD = SECOND()
      DO ISYMB = 1,NSYM
         ISYMA = MULD2H(ISYMB,ISYDEN)
         DO IB = 1,NBAS(ISYMB)
            B = IBAS(ISYMB) + IB
            KOFF1 = KFOCK + IAODIS(ISYMA,ISYMB)
     &            + NBAS(ISYMA)*(IB - 1)
            KOFF2 = NBAST*(B - 1) + IBAS(ISYMA) + 1
            CALL DAXPY(NBAS(ISYMA),ONE,WORK(KOFF1),1,FOCK(KOFF2),1)
         ENDDO
      ENDDO
      TIMAD = SECOND() - TIMAD

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Exit')
         CALL OUTPUT(DEN,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         CALL AROUND(SECNAM//': Fock Matrix on Exit')
         CALL OUTPUT(FOCK,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF

C     Print timing.
C     -------------

      IF (TINFO .OR. LOCDBG) THEN
         TIMT = SECOND() - TIMT
         WRITE(LUPRI,'(/,5X,A,A,A)')
     &   ' - exiting ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sym. packing density : ',TIMPK,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for calculating contr.   : ',TIMCL,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for adding to Fock matrix: ',TIMAD,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                         : ',TIMT,' seconds'
      ENDIF

      RETURN
      END
C  /* Deck scf_chofck3 */
      SUBROUTINE SCF_CHOFCK3(DEN,FOCK,WORK,LWORK,ISYDEN,TINFO,IOPTVP)
C
C     Written by Thomas Bondo Pedersen, May 2002.
C
C     Purpose:
C        Calculate two-electron integral contributions to the AO Fock
C        matrix based on Cholesky decomposed integrals.
C        Result is added into Fock matrix which is assumed allocated
C        as symmetry packed (!) full square.
C        The density matrix is used for both Coulomb and exchange part.
C
C     Formula used:
C        F(ab) = F(ab) + 2 * Sum(gd) D(dg) (ab|dg)
C                      -     Sum(gd) D(dg) (ad|bg)
C
C     NB: NO SPARSITY IS USED WHATSOEVER !!!!!
C
#include <implicit.h>
      DIMENSION DEN(*), FOCK(*), WORK(LWORK)
      LOGICAL TINFO
#include <iratdef.h>
#include <infpri.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccdeco.h>
#include <priunit.h>
#include <dftcom.h>

      INTEGER IOFFC(8)

      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'SCF_CHOFCK3')

      LOGICAL IOPTVP

C     Initialize timings.
C     -------------------

      TIMT = SECOND()
      TIMR = ZERO
      TIMS = ZERO
      TIMC = ZERO
      TIMX = ZERO

C     Read reduce index array.
C     ------------------------

      KIND1 = 1
      CALL CC_GETIND1(WORK(KIND1),LWORK,LIND1)
      KEND0 = KIND1 + LIND1
      LWRK0 = LWORK - KEND0 + 1

      IF (LWRK0 .LT. 0) THEN
         WRITE(LUPRI,'(//,5X,A,A,A)')
     &   'Insufficient memory in ',SECNAM,' - allocation: index'
         WRITE(LUPRI,'(5X,A,I10,/,5X,A,I10,/)')
     &   'Need (more than): ',KEND0-1,
     &   'Available       : ',LWORK
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Entry')
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYDEN)
            WRITE(LUPRI,'(/,15X,A,I1,1X,I1)')
     &      'Symmetry Block (alpha,beta): ',ISYMA,ISYMB
            KOFF = IAODIS(ISYMA,ISYMB) + 1
            NUMB = NBAS(ISYMB)
            NUMA = NBAS(ISYMA)
            CALL OUTPUT(DEN(KOFF),1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': Fock Matrix on Entry')
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYDEN)
            WRITE(LUPRI,'(/,15X,A,I1,1X,I1)')
     &      'Symmetry Block (alpha,beta): ',ISYMA,ISYMB
            KOFF = IAODIS(ISYMA,ISYMB) + 1
            NUMB = NBAS(ISYMB)
            NUMA = NBAS(ISYMA)
            CALL OUTPUT(FOCK(KOFF),1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
         ENDDO
      ENDIF

C     Start Cholesky loop.
C     --------------------

      DO ISYCHO = 1,NSYM

         IF (N2BST(ISYCHO)  .LE. 0) GO TO 999
         IF (NUMCHO(ISYCHO) .LE. 0) GO TO 999

         ISYMX = MULD2H(ISYCHO,ISYDEN)

C        Allocation: Cholesky.
C        ---------------------

         LREAD = 2*NNBST(ISYCHO) + (NNBST(ISYCHO) - 1)/2 + 2

         KCSCR = KEND0
         KREAD = KCSCR + N2BST(ISYCHO)
         KEND1 = KREAD + LREAD
         LWRK1 = LWORK - KEND1 + 1

         IF (LWRK1 .LT. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - allocation: Cholesky'
            WRITE(LUPRI,'(5X,A,I10,/,A,I10,/)')
     &      'Need     : ',KEND1-1,
     &      'Available: ',LWORK
            CALL QUIT('INSUFFICIENT MEMORY IN '//SECNAM)
         ENDIF

C        Set up batch.
C        -------------

         MAXAG = -1
         DO ISYMG = 1,NSYM
            ISYMA = MULD2H(ISYMG,ISYMX)
            MAXAG = MAX(MAXAG,NBAS(ISYMA)*NBAS(ISYMG))
         ENDDO

         MINMEM = N2BST(ISYCHO) + MAXAG
         NVEC   = MIN(LWRK1/MINMEM,NUMCHO(ISYCHO))

         IF (NVEC .LE. 0) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'Insufficient memory in ',SECNAM,
     &      ' - Allocation: Batch'
            WRITE(LUPRI,'(5X,A,I10)')
     &      'Need per vector    : ',MINMEM
            WRITE(LUPRI,'(5X,A,I10,/)')
     &      'Available for batch: ',LWRK1
            CALL QUIT('INSUFFICIENT MEMORY IN '//SECNAM)
         ENDIF

         NBATCH = (NUMCHO(ISYCHO) - 1)/NVEC + 1

         IF (LOCDBG) THEN
            WRITE(LUPRI,'(//,5X,A,A,A)')
     &      'BATCH INFO FROM ',SECNAM,':'
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Cholesky symmetry                : ',ISYCHO
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Total number of vectors          : ',NUMCHO(ISYCHO)
            WRITE(LUPRI,'(8X,A,I10)')
     &      'AO Cholesky vector dimension     : ',N2BST(ISYCHO)
            WRITE(LUPRI,'(8X,A,I10)')
     &      'AO intermediate vector dimension : ',MAXAG
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Total memory                     : ',LWORK
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Min. memory needed for batch     : ',MINMEM
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Number of batches                : ',NBATCH
            WRITE(LUPRI,'(8X,A,I10)')
     &      'Number of vectors per batch      : ',NVEC
         ENDIF

         DO IBATCH = 1,NBATCH

            NUMV = NVEC
            IF (IBATCH .EQ. NBATCH) THEN
               NUMV = NUMCHO(ISYCHO) - NVEC*(NBATCH - 1)
            ENDIF

            JCHOL1 = NVEC*(IBATCH - 1) + 1

C           Complete allocation.
C           --------------------

            KCHOL = KEND1
            KXMAT = KCHOL + N2BST(ISYCHO)*NUMV

C           Set up index array IOFFC.
C           -------------------------

            ICOUNT = 0
            DO ISYMG = 1,NSYM
               ISYMB = MULD2H(ISYMG,ISYCHO)
               IOFFC(ISYMG) = ICOUNT
               ICOUNT = ICOUNT + NBAS(ISYMB)*NUMV*NBAS(ISYMG)
            ENDDO

            DO JVEC = 1,NUMV

               JCHOL = JCHOL1 + JVEC - 1

C              Read vector.
C              ------------

               DTIME = SECOND()
               IOPT  = 2
               CALL CHO_READN(WORK(KCSCR),JCHOL,1,WORK(KIND1),IDUM2,
     &                        ISYCHO,IOPT,WORK(KREAD),LREAD)
               DTIME = SECOND() - DTIME
               TIMR  = TIMR     + DTIME

C              Coulomb part.
C              -------------

               IF (IOPTVP .AND. (ISYCHO .EQ. ISYDEN)) THEN
                  DTIME = SECOND()
                  CONST = TWO*DDOT(N2BST(ISYDEN),DEN,1,WORK(KCSCR),1)
                  CALL DAXPY(N2BST(ISYDEN),CONST,WORK(KCSCR),1,
     &                       FOCK,1)
                  DTIME = SECOND() - DTIME
                  TIMC  = TIMC     + DTIME
               ENDIF

C              Reorder: L(bg,J) -> L(b,J,g).
C              -----------------------------

               DTIME = SECOND()
               DO ISYMG = 1,NSYM

                  ISYMB = MULD2H(ISYMG,ISYCHO)

                  IF (NBAS(ISYMB) .GT. 0) THEN

                     DO G = 1,NBAS(ISYMG)

                        KOFF1 = KCSCR + IAODIS(ISYMB,ISYMG)
     &                        + NBAS(ISYMB)*(G - 1)
                        KOFF2 = KCHOL + IOFFC(ISYMG)
     &                        + NBAS(ISYMB)*NUMV*(G - 1)
     &                        + NBAS(ISYMB)*(JVEC - 1)

                        CALL DCOPY(NBAS(ISYMB),WORK(KOFF1),1,
     &                                         WORK(KOFF2),1)

                     ENDDO

                  ENDIF

               ENDDO
               DTIME = SECOND() - DTIME
               TIMS  = TIMS     + DTIME

            ENDDO

            DTIME = SECOND()
            DO ISYMG = 1,NSYM

               ISYMB = MULD2H(ISYMG,ISYCHO)
               ISYMD = MULD2H(ISYMG,ISYDEN)
               ISYMA = MULD2H(ISYMB,ISYDEN)

C              Calculate intermediate:
C              X(aJ,g) = Sum(d) L(aJ,d) * D(g,d).
C              ----------------------------------

               NTOAJ = MAX(NBAS(ISYMA)*NUMV,1)
               NTOTD = MAX(NBAS(ISYMD),1)

               KOFF1 = KCHOL + IOFFC(ISYMD)
               KOFF2 = IAODIS(ISYMD,ISYMG) + 1

               CALL DGEMM('N','N',
     &                    NBAS(ISYMA)*NUMV,NBAS(ISYMG),NBAS(ISYMD),
     &                    ONE,WORK(KOFF1),NTOAJ,DEN(KOFF2),NTOTD,
     &                    ZERO,WORK(KXMAT),NTOAJ)

C              Calculate exchange contribution:
C              F(a,b) = F(a,b) - Sum(g) X(a,Jg) * L(b,Jg).
C              -------------------------------------------

               NTOTA = MAX(NBAS(ISYMA),1)
               NTOTB = MAX(NBAS(ISYMB),1)

               KOFF3 = KCHOL + IOFFC(ISYMG)
               KOFF4 = IAODIS(ISYMA,ISYMB) + 1

               CALL DGEMM('N','T',
     &                    NBAS(ISYMA),NBAS(ISYMB),NUMV*NBAS(ISYMG),
     &                    HFXFAC*XMONE,WORK(KXMAT),NTOTA,WORK(KOFF3),
     &                    NTOTB,ONE,FOCK(KOFF4),NTOTA)

            ENDDO
            DTIME = SECOND() - DTIME
            TIMX  = TIMX     + DTIME

         ENDDO

  999    CONTINUE

      ENDDO

C     Debug: Print matrices.
C     ----------------------

      IF (LOCDBG) THEN
         CALL AROUND(SECNAM//': Density Matrix on Exit')
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYDEN)
            WRITE(LUPRI,'(/,15X,A,I1,1X,I1)')
     &      'Symmetry Block (alpha,beta): ',ISYMA,ISYMB
            KOFF = IAODIS(ISYMA,ISYMB) + 1
            NUMB = NBAS(ISYMB)
            NUMA = NBAS(ISYMA)
            CALL OUTPUT(DEN(KOFF),1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
         ENDDO
         CALL AROUND(SECNAM//': Fock Matrix on Exit')
         DO ISYMB = 1,NSYM
            ISYMA = MULD2H(ISYMB,ISYDEN)
            WRITE(LUPRI,'(/,15X,A,I1,1X,I1)')
     &      'Symmetry Block (alpha,beta): ',ISYMA,ISYMB
            KOFF = IAODIS(ISYMA,ISYMB) + 1
            NUMB = NBAS(ISYMB)
            NUMA = NBAS(ISYMA)
            CALL OUTPUT(FOCK(KOFF),1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
         ENDDO
      ENDIF

C     Print timing.
C     -------------

      IF (TINFO .OR. LOCDBG) THEN
         TIMT = SECOND() - TIMT
         WRITE(LUPRI,'(/,5X,A,A,A)')
     &   ' - exiting ',SECNAM,':'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for I/O and squaring     : ',TIMR,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for Coulomb  part        : ',TIMC,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for sorting (copy)       : ',TIMS,' seconds'
         WRITE(LUPRI,'(5X,A,F10.2,A)')
     &   'Time used for exchange part        : ',TIMX,' seconds'
         WRITE(LUPRI,'(5X,A)')
     &   '-------------------------------------------------------'
         WRITE(LUPRI,'(5X,A,F10.2,A,/)')
     &   'Total time                         : ',TIMT,' seconds'
      ENDIF

      RETURN
      END

