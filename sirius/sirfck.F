C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
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
C     DIRECT = .TRUE.   Direct calculation using hermit.
C              .FALSE.  Read integrals from AOTWOINT.
#include <implicit.h>
      LOGICAL   DIRECT, DOERG, DOATR, LSRHYBR
      DIMENSION FMAT(N2BASX,NDMAT), DMAT(N2BASX,*)
      DIMENSION ISYMDM(NDMAT), IFCTYP(NDMAT), MBAS(8)
      DIMENSION WRK(*), IFCTMP(10)
      DIMENSION DFCOEFF(NBAST) !denisty-fitting coefficients
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
      LOGICAL SFIRST,PREHF
      SAVE SFIRST
      DATA SFIRST  /.TRUE./
C
C
#include <memint.h>
      IF (NDMAT.EQ.0) RETURN
      CALL QENTER('SIRFCK')
!      print *,'SIRFCK: entering'
!      print *,'SIRFCK: NDMAT=', NDMAT
!      print *,'SIRFCK: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
C     ...PREHF: Find out if this is a HF calculation preceeding
C        a Multiconfigurational SR-DFT calculation.
      PREHF = DOSCF.AND.(DOMCSRDFT.OR.DOCISRDFT).AND.
     &                  (.NOT.DOHFSRDFT).AND.(MCTYPE.EQ.0)
      LSRHYBR = SRHYBR .AND. .NOT.DOHFSRDFT .AND. .NOT.PREHF
C
      IF (DIRECT) THEN
         IF (SFIRST) THEN
            IPRHER = 2*IPRFCK
            CALL SETHER(IPRHER,.FALSE.,WRK,LWORK)
C           CALL SETHER(IPRHER,NEWGEO,WORK,LWORK)
            SFIRST = .FALSE.
         END IF
         DO I = 1,NDMAT
            ISYMDM(I) = ISYMDM(I) - 1
         END DO
C
C        Special code for short range/long range integrals
C        ... never short range integrals in first call of HERFCK
C        ... if normal HF then normal 2-el ints (CHIVAL=0)
         SRINTS = .FALSE.
         CHIVALSAVE = CHIVAL
         HFXFACSAVE = HFXFAC
         IF (ADDSRI .OR. PREHF) HFXFAC = D1
C        ... full exchange for long range part
         IF (ADDSRI .AND. .NOT.(DOHFSRDFT.OR.DOMCSRDFT)) THEN
C           this is e.g. PREHF, i.e. standard HF
            CHIVAL = D0
         END IF
C
         IF (NSYM.GT.1 .AND. (DODFT.OR.DOHFSRDFT.OR.DOMCSRDFT)) THEN
C           HERFCK destroys DMAT when called with symmetry
            KDMAT = 1
            KLAST = KDMAT + NDMAT*N2BASX
            LWRK  = LWORK - KLAST + 1
            IF (KLAST.GT.LWORK) CALL STOPIT('SIRFCK','DMAT',KLAST,LWORK)
            CALL DCOPY(NDMAT*N2BASX,DMAT,1,WRK(KDMAT),1)
            if (atomdf .AND. CHIVAL .LE. D0) then
C             hjaaj mar 2006: dens_fit not for short range integrals (yet)
              call dens_fit(fmat,wrk(kdmat),ndmat,ISYMDM,IFCTYP,IPRFCK,
     &                    WRK(KLAST),LWRK)
            else
              CALL HERFCK(FMAT,WRK(KDMAT),NDMAT,ISYMDM,IFCTYP,IPRFCK,
     &                    WRK(KLAST),LWRK)
            end if
         ELSE
            if (atomdf .AND. CHIVAL .LE. D0) then
C             hjaaj mar 2006: dens_fit not for short range integrals (yet)
              call dens_fit(fmat,dmat,ndmat,ISYMDM,IFCTYP,IPRFCK,
     &                      wrk,lwork)
            else
              CALL HERFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,IPRFCK,
     &                      WRK,LWORK)
            end if
         END IF
         CHIVAL = CHIVALSAVE
         IF (.NOT. PREHF) HFXFAC = HFXFACSAVE
C
C        Add Short-range integrals (HFSRDFT hybrid model)
C
         IF (ADDSRI .AND. (DOHFSRDFT.OR.DOMCSRDFT)) THEN
            IF (DOMCSRDFT)
     &      CALL QUIT('DOMCSRDFT not implemented DIRECT')
C           ... we need to add density matrices as below for not DIRECT
C               in order to get this working
            IF (NDMAT .GT. 10)
     &         CALL QUIT('NDMAT .gt. 10 in SIRFCK, increase dimension!')
            SRINTS = .TRUE.
            IF (HFXFAC .NE. 0.0D0) THEN
C           ... Use exact SR-exchange for the HF-DFT hybrid
C               (i.e both exchange and coulomb : IFCTYP = 13)
                IFCVAL = 13
            ELSE
C           ... Exchange in functional ->  add only Coul (IFCTYP = 11)
                IFCVAL = 11
            ENDIF
            DO I = 1,NDMAT
               IFCTMP(I) = IFCVAL
            END DO
!            print *,'SIRFCK: after IFCTMP(I) = IFCVAL'       !JT
!            print *,'SIRFCK: NDMAT=',NDMAT                   !JT
!            print *,'SIRFCK: IFCTMP=',(IFCTMP(I),I=1,NDMAT)  !JT
C           HERFCK zeroes FMAT with long range term, thus :
            KFMAT = 1
            KLAST = KFMAT + NDMAT*N2BASX
            LWRK  = LWORK - KLAST + 1
            IF (KLAST.GT.LWORK) CALL STOPIT('SIRFCK','FMAT',KLAST,LWORK)
            CALL HERFCK(WRK(KFMAT),DMAT,NDMAT,ISYMDM,IFCTMP,IPRFCK,
     &                  WRK(KLAST),LWRK)
            CALL DAXPY(NDMAT*N2BASX,D1,WRK(KFMAT),1,FMAT,1)
            SRINTS = .FALSE.
         ENDIF
         DO I = 1,NDMAT
            ISYMDM(I) = ISYMDM(I) + 1
         END DO
         HFXFAC = HFXFACSAVE
      ELSE
         IF (HFXMU .NE. D0) CALL QUIT(
     &       'SIRFCK2: nonzero HFXMU requires .DIRECT calculation.') 
C     ... not DIRECT
C
C        ... if SRHYBR (DFT-SRDFT hybrid) only include Coulomb terms 
C            (core is treated with BLYP -> hence no exchange here!)
         IF(LSRHYBR) THEN
            IFCTMP(1) = 11
C           ... FC,lr = JC,lr = J,lr * DC for SR-hybrid
            IF (NDMAT .EQ. 2) IFCTMP(2) = 13
C           ... FV,lr = G,lr * DV
            IF (NDMAT .GT. 2) CALL QUIT('NDMAT.gt.2 for SRHYBR & lr')
         ELSE
!JT           IFCTMP(1) = IFCTYP(1)
!JT           IFCTMP(2) = IFCTYP(2)
           DO I = 1,NDMAT          !JT
            IFCTMP(I) = IFCTYP(I)  !JT
           ENDDO                   !JT
         END IF
         HFXFACSAVE = HFXFAC
         IF (ADDSRI .OR. PREHF) HFXFAC = D1
!         print *,'SIRFCK: calling SIRFCK2 with'
!         print *,'SIRFCK: NDMAT=',NDMAT
!         print *,'SIRFCK: IFCTMP=',(IFCTMP(I),I=1,NDMAT)
         if (atomdf .AND. CHIVAL .LE. D0) then
C             hjaaj mar 2006: dens_fit not for short range integrals (yet)
           call dens_fit(fmat,dmat,ndmat,ISYMDM,IFCTMP,IPRFCK,wrk,lwork)
         else
           CALL SIRFCK2(LUINTA,'AOTWOINT',
     &                FMAT,DMAT,NDMAT,ISYMDM,IFCTMP,WRK,LWORK)
         end if
C
C        Add Short-range integrals to FC part of FMAT (SRDFT hybrid models)
C
         IF (ADDSRI) THEN
            IF (.NOT.PREHF) HFXFAC = HFXFACSAVE
             IF (NASHT.GT.0 .AND. (.NOT.DOCISRDFT))THEN
C              ... get total DMAT = DC + DV if NASHT.gt.0
C                  for Hartree SR term
               NSRDMAT = NDMAT/2
               IF (2*NSRDMAT .NE. NDMAT) THEN
                  WRITE(LUPRI,'(3X,A,3I4)')
     &           'SIRFCK : NSRDMAT, NDMAT, NASHT = ',NSRDMAT,NDMAT,NASHT
                  CALL QUIT('2*NSRDMAT.ne.NDMAT for '//
     &                      'ADDSRI and NASHT.gt.0')
               ENDIF
               DO I = 1,NSRDMAT
                  CALL DAXPY(N2BASX,D1,DMAT(1,NSRDMAT+I),1,DMAT(1,I),1)
C                 ... add DV(I) to DC(I)
               END DO
            ELSE
C              ... total DMAT = DC when NASHT.eq.0 or CISRDFT
               NSRDMAT = NDMAT
            END IF
C
            IF(HFXFAC.NE.D0) THEN
C           ... for normal/presceding HF and for SRDFT with exact exch.
C               (add SR coul. and exch. to LR to give full 2-el part)
                IFCVAL = 13
            ELSE 
C           ... add only Coul
C               (only SR Hartree term for SRDFT hybrids)
                IFCVAL = 11
            END IF
            IF (NSRDMAT .GT. 10)
     &         CALL QUIT('NSRDMAT.gt.10 in SIRFCK, increase dimension!')
            DO I = 1,NSRDMAT
               IFCTMP(I) = IFCVAL
            END DO
            IF (NSRDMAT.GT.1) CALL QUIT('NSRDMAT .GT. 1 for ADDSRI')
            KSRAO = 1
            KLAST = KSRAO  + NSRDMAT*N2BASX
            LWRK  = LWORK - KLAST + 1
            CALL DZERO(WRK(KSRAO),NSRDMAT*N2BASX)
            CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KSRAO),
     &                   DMAT,NSRDMAT,ISYMDM,IFCTMP,WRK(KLAST),LWRK)
            CALL DAXPY(N2BASX,D1,WRK(KSRAO),1,FMAT,1)
            EJSR = DDOT(N2BASX,DMAT(1,1),1,WRK(KSRAO),1)
            IF (NASHT.GT.0) THEN
               EJVSR =-DP5*DDOT(N2BASX,DMAT(1,NSRDMAT+1),1,WRK(KSRAO),1)
            ELSE
               EJVSR = D0
            END IF
            if (HFXFAC .NE. 0.0D0) then
               nsrdmat = 1
               IFCTMP(1) = 12
               CALL DZERO(WRK(KSRAO),NSRDMAT*N2BASX)
               CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KSRAO),
     &                   DMAT,NSRDMAT,ISYMDM,IFCTMP,WRK(KLAST),LWRK)
               EKSR = DDOT(N2BASX,DMAT(1,1),1,WRK(KSRAO),1)
               EJSR = EJSR - EKSR
            else
               EKSR = D0
            end if
            IF (NASHT.GT.0) THEN
               DO I = 1,NSRDMAT
                  CALL DAXPY(N2BASX,-D1,DMAT(1,NSRDMAT+I),1,DMAT(1,I),1)
C                 ... restore DC(I) = DTOT(I) -  DV(I)
               END DO
            END IF
         END IF
         HFXFAC = HFXFACSAVE
      END IF
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
C     ADD SR-DFT contribution if MCSCF-DFT or if HF-SRDFT
C
      IF (ADDSRI.AND.(DOHFSRDFT.OR.DOMCSRDFT).AND.(.NOT.PREHF)) THEN
         IF (NDMAT .GT. 2)
     &   CALL QUIT('SR-DFT not implemented for NDMAT.gt.2')
C        we need to do as above for ADDSRI !!!!!!!
         DOERG = DFTADD
         DOATR = .NOT.DFTADD
C        ... DFTADD tells if add normal DFT Vxc or second order Vxc
         IF (LSRHYBR) THEN
            NSRMAT = 2
         ELSE
            NSRMAT = 1
         END IF
         KSRAO  = 1
         KLAST  = KSRAO  + NSRMAT*N2BASX
         LWRK  = LWORK - KLAST + 1
         IF (DOATR) THEN
            JDAO = NDMAT + 1
C           ... DTOTAO is in DMAT(1,NDMAT+1) from TR1FCK.
C           For SRHYBR we also need DVAO in DMAT(1,NDMAT+2) from TR1FCK.:
            IF (LSRHYBR) THEN
               JDXAO = JDAO + 2
            ELSE
               JDXAO = JDAO + 1
            END IF
C
C           we need DXTOTAO in DMAT(1,JDXAO) for SRDFT :
            CALL DCOPY(N2BASX,DMAT(1,1),1,DMAT(1,JDXAO),1)
            IF (NDMAT.EQ.2) THEN
               CALL DAXPY(N2BASX,D1,DMAT(1,2),1,DMAT(1,JDXAO),1)
C              ... DXTOTAO = DXCAO + DXVAO
               IF (LSRHYBR)
     &         CALL DCOPY(N2BASX,DMAT(1,2),1,DMAT(1,JDXAO+1),1)
C              ... copy DXVAO
            END IF
         ELSE
            JDAO = 1
            IF (NDMAT.EQ.2)
     &      CALL DAXPY(N2BASX,D1,DMAT(1,2),1,DMAT(1,1),1)
C           ... get total DMAT = DC + DV if NDMAT .gt. 1
         END IF
         CALL DZERO(WRK(KSRAO),NSRMAT*N2BASX)
C        SRDFT(EXCMAT,DMAT,EDFTY,DOERG,DOATR,WORK,LWORK,IPRINT)
!         WRITE(LUPRI,*) 'SIRFCK: calling SRDFT'  !JT
         CALL SRDFT(WRK(KSRAO),DMAT(1,JDAO),EDFT,DOERG,DOATR,
     &              WRK(KLAST),LWRK,IPRFCK)
C        IF (DOATR)
C    &      ESRLTR = ESRLTR + DDOT(N2BASX,WRK(KSRAO),1,DMAT(1,JDAO),1)
C           ... for MCSRDFT csf sigma-vectors : - 2 ESRLTR cref(i)
C               NOT NEEDED because we project out CREF component in sigma vectors
         IF (NDMAT.EQ.2 .AND. .NOT.DOATR)
     &      CALL DAXPY(N2BASX,-D1,DMAT(1,2),1,DMAT(1,1),1)
C        ... restore DC = DTOT - DV
         CALL DAXPY(N2BASX,D1,WRK(KSRAO),1,FMAT,1)
         IF (LSRHYBR) CALL DCOPY(N2BASX,
     &      WRK(KSRAO+N2BASX),1,FMAT(1,NDMAT+1),1)
         IF (DOERG) THEN
            ECSR = -DP5*DDOT(N2BASX,DMAT(1,1),1,WRK(KSRAO),1)
            IF (NASHT.GT.0) THEN
               EVSR = -DDOT(N2BASX,DMAT(1,2),1,WRK(KSRAO),1)
            ELSE
               EVSR = D0
            END IF
            ESRDFTY = EDFT + ECSR + EVSR + EJVSR
C           IF (IPRFCK .GT. 0)
            WRITE (LUPRI,'(A,F25.10)')
     &         '   EDFT    ',EDFT,
     &         ' + ECSR    ',ECSR,
     &         ' + EVSR    ',EVSR,
     &         ' + EJVSR   ',EJVSR,
     &         ' = ESRDFTY ',ESRDFTY,
     &         '   EJSR    ',EJSR,
     &         '   EKSR    ',EKSR
         END IF
      ENDIF
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
      CALL QEXIT('SIRFCK')
      RETURN
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
!
!      print *,'SIRFCK2: NDMAT=', NDMAT
!      print *,'SIRFCK2: IFCTYP=', (IFCTYP(I),I=1,NDMAT)
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
#if !defined (SYS_CRAY) && !defined (SYS_T3D) && !defined (SYS_T90)
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
