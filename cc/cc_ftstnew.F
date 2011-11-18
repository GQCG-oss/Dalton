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
*=====================================================================*
      SUBROUTINE CC_FTSTNEW(WORK,LWORK,APROXR12)
      IMPLICIT NONE
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccfro.h>
#include <dummy.h>

      CHARACTER*3 APROXR12,LISTL,LISTR
      CHARACTER*6 FILFMAT
      CHARACTER*8 CDUMMY
      CHARACTER*10 MODEL,MODELW

      LOGICAL LOCDBG
      PARAMETER(LOCDBG = .FALSE.)
      INTEGER NFTRAN,IFTRAN(3,1),IOPTRES,KRHOANA1,KRHOANA2,
     &        KEND1,KRHOANA12,KRHODIF1,KRHODIF2,KRHODIF12,LWORK,
     &        LWRK1,IDLSTL,IDLSTR,ISYM,IR1TAMP,LENALL,LUFMAT,
     &        iopt,kt0amp0,kt12,lenmod,ISYCTR,ISYAMP,ISYRES,
     &        isymkl,isymij

#if defined (SYS_CRAY) 
      REAL WORK(LWORK), ERRNRM, DDOT
#else
      DOUBLE PRECISION WORK(LWORK), ERRNRM, DDOT
#endif
     
      !external functions:
      INTEGER ILSTSYM

      WRITE(LUPRI,*) 'Entered CC_FTSTNEW!'

      ISYM   = 1

      LISTL  = 'L0 '
      IDLSTL = 0

      LISTR  = 'R1 '
      IDLSTR = IR1TAMP('ZDIPLEN ',.FALSE.,0.0D0,ISYM)
      NFTRAN = 1

      IFTRAN(1,1) = IDLSTL
      IFTRAN(2,1) = IDLSTR

      IOPTRES = 0
      FILFMAT = 'CC_FMAT'

C     determine symmetry of result vector: 
      ISYCTR = ILSTSYM(LISTL,IDLSTL)
      ISYAMP = ILSTSYM(LISTR,IDLSTR)
      ISYRES = MULD2H(ISYAMP,ISYCTR)
      if (locdbg) then
        write(lupri,*) 'ISYCTR, ISYAMP, ISYRES:', isyctr, isyamp, isyres
      end if

C Models:
      if (CCS) then
        modelw = 'CCS       '
      else if (CC2) then
        modelw = 'CC2       '
      else if (CCSD) then
        modelw = 'CCSD      '
      else
        call quit('unknown CC method in CC_FTSTNEW!')
      end if
      CALL CCSD_MODEL(MODELW,LENMOD,24,MODELW,10,APROXR12)
C     write(lupri,*) 'MODEL in CC_FTSTNEW: ',MODELW
     
C     Generate R12-parts for Lagrange-multipliers and trial-vector
C     from ground state amplitudes (which are total-symmetric)
C     (just to have some values...)
C     if (CCR12) then
C       memory allocation
C       kt0amp0 = 1
C       kt12 = kt0amp0 + 2*nallai(1) 
C       kend1  = kt12 + ngamma(1)
C       lwrk1  = lwork - kend1
C       if (lwrk1 .lt. 0) then
C         call quit('Insufficient memory in CC_FTSTNEW for R12!')
C       end if
C
C       check symmetry:
C       if (ISYCTR .NE. 1) then
C         call quit('Symmetry mismatch while generating 
C    &               Lagrange-multipliers')
C       else if (ISYAMP .NE. 1) then 
C         call quit('Symmetry mismatch while generating response
C    &               amplitudes')
C       end if
C       
C       iopt=32
C       call cc_rdrsp('R0 ',0,1,iopt,model,dummy,work(kt12))
C
C       if (locdbg) then
C         write(lupri,*) 'Norm^2(Ground state R12-ampl.)',
C    &      ddot(ngamma(1),work(kt12),1,work(kt12),1) 
C         call cc_prpr12(work(kt12),1,1,.TRUE.)
C       end if
C
C
C       do i=1, ngamma(1)
C         work(kt12 + i - 1) = 100.0d0 * 
C    &           work(kt12 + i - 1)*dabs(work(kt12+i-1))
C       end do
C
C       if (locdbg) then
C         write(lupri,*) 'Norm^2(R12-Zeta)',
C    &      ddot(ngamma(1),work(kt12),1,work(kt12),1)
C         call cc_prpr12(work(kt12),1,1,.TRUE.) 
C       end if
C
C       iopt=32
C       call cc_wrrsp(listl,idlstl,1,iopt,modelw,work(kt0amp0),dummy,
C    &                work(kt12),work(kend1),lwrk1)
C       do i=1, ngamma(1)
C         work(kt12 + i - 1) = 100.0d0 * 
C    &           work(kt12 + i - 1)*dabs(work(kt12+i-1))
C       end do
C
C       if (locdbg) then
C         write(lupri,*) 'Norm^2(R12 response ampl.)',
C    &      ddot(ngamma(1),work(kt12),1,work(kt12),1) 
C         call cc_prpr12(work(kt12),1,1,.TRUE.)
C       end if
C
C
C       iopt=32
C       call cc_wrrsp(listr,idlstr,1,iopt,modelw,work(kt0amp0),dummy,
C    &                work(kt12),work(kend1),lwrk1)
C     end if 
CCN
C     !TEST: Zero ALL singles amplit., resp. vecs. and lagr. mult.:
C     CALL DZERO(WORK,NT1AM(1))
C     IOPT = 1
C     CALL CC_WRRSP('R0 ',0,1,IOPT,MODELW,WORK,WORK,DUMMY,
C    &              WORK(1+NT1AM(1)),LWORK-NT1AM(1)-1)
C     IOPT = 1
C     CALL CC_WRRSP(LISTR,IDLSTR,1,IOPT,MODELW,WORK,WORK,DUMMY,
C    &              WORK(1+NT1AM(1)),LWORK-NT1AM(1)-1)
C     IOPT = 1
C     CALL CC_WRRSP(LISTL,IDLSTL,1,IOPT,MODELW,WORK,WORK,DUMMY,
C    &              WORK(1+NT1AM(1)),LWORK-NT1AM(1)-1)
CCN
      CALL CC_FMATRIX(IFTRAN,NFTRAN,LISTL,LISTR,IOPTRES,FILFMAT,
     &                IDUMMY,DUMMY,0,WORK,LWORK)

      KRHOANA1  = IFTRAN(3,1)
      KRHOANA2  = KRHOANA1 + NT1AM(ISYRES)
      KEND1     = KRHOANA2 + NT2AM(ISYRES)
      IF (CCR12) THEN
        KRHOANA12 = KEND1
        KEND1     = KRHOANA12 + NTR12AM(ISYRES)
      END IF  
      LWRK1 = LWORK - KEND1
      IF (LWRK1.LT.0) CALL QUIT('Insufficient memory in CC_FTSTNEW')

      LENALL = NT1AM(ISYRES) + NT2AM(ISYRES)
      IF (CCS)   LENALL = NT1AM(ISYRES)
      IF (CCR12) LENALL = LENALL + NTR12AM(ISYRES)

      CALL DZERO(WORK(KRHOANA1),LENALL)

      LUFMAT = -1
      CALL WOPEN2(LUFMAT,FILFMAT,64,0)
      CALL GETWA2(LUFMAT,FILFMAT,WORK(KRHOANA1),1,LENALL) 
      CALL WCLOSE2(LUFMAT,FILFMAT,'DELETE')

      IF (NSYM .EQ. 1) then
        KRHODIF1 = KEND1
        KRHODIF2 = KRHODIF1 + NT1AMX
        KEND1    = KRHODIF2 + NT2AMX
        IF (CCR12) THEN
          KRHODIF12 = KEND1
          KEND1     = KRHODIF12 + NTR12AM(1)
        END IF

        LWRK1 = LWORK - KEND1
        IF (LWRK1.LT.0) CALL QUIT('Insufficient memory in CC_FTSTNEW')
        CALL DZERO(WORK(KRHODIF1),LENALL)

        CALL CC_FDF(WORK(KRHODIF1),WORK(KRHODIF2),WORK(KRHODIF12),
     &              LISTL,IDLSTL,LISTR,IDLSTR,
     &              WORK(KEND1),LWRK1,APROXR12)
      ELSE
        write (lupri,*) 'CC_FDF must run in C1 symmetry!'
        write (lupri,*) 'Will ignore finite differences....'
      END IF

      N = 1
      IF (CCS) N = 0

      IF (NSYM .EQ. 1) THEN
        WRITE(LUPRI,'(/A)') 'Finite difference Result for FMATRIX:'
        CALL CC_PRP(WORK(KRHODIF1),WORK(KRHODIF2),1,1,N)
        IF (CCR12) THEN
           CALL CC_PRPR12(WORK(KRHODIF12),1,1,.TRUE.)
        END IF
      END IF
      WRITE(LUPRI,'(//A)') 'Analytical Result for FMATRIX:'
      CALL CC_PRP(WORK(KRHOANA1),WORK(KRHOANA2),ISYRES,1,N)
      IF (CCR12) THEN
         CALL CC_PRPR12(WORK(KRHOANA12),ISYRES,1,.TRUE.)
         WRITE(LUPRI,*) 'Norm^2 of FMATRIX: ',
     &         DDOT(NT1AM(ISYRES),WORK(KRHOANA1),1,WORK(KRHOANA1),1) +
     &         DDOT(NT2AM(ISYRES),WORK(KRHOANA2),1,WORK(KRHOANA2),1) +
     &         DDOT(NTR12AM(ISYRES),WORK(KRHOANA12),1,WORK(KRHOANA12),1)
      END IF

      IF (NSYM .EQ. 1) THEN
        CALL DAXPY(NT1AMX,-1.0D0,WORK(KRHODIF1),1,WORK(KRHOANA1),1)
        CALL DAXPY(NT2AMX,-1.0D0,WORK(KRHODIF2),1,WORK(KRHOANA2),1)
        ERRNRM = DDOT(NT1AMX,WORK(KRHOANA1),1,WORK(KRHOANA1),1) +
     &           DDOT(NT2AMX,WORK(KRHOANA2),1,WORK(KRHOANA2),1)
        IF (CCR12) THEN
          CALL DAXPY(NTR12AM(1),-1.0D0,WORK(KRHODIF12),1,
     &               WORK(KRHOANA12),1)
          ERRNRM = ERRNRM + 
     &             DDOT(NTR12AM(1),WORK(KRHOANA12),1,WORK(KRHOANA12),1)
        END IF  
     
        WRITE(LUPRI,'(//A)') 'Norm of Difference between Results:'
        WRITE(LUPRI,*) 'Norm:',ERRNRM
        WRITE(LUPRI,*) 'Singles part: ',
     &                 DDOT(NT1AMX,WORK(KRHOANA1),1,WORK(KRHOANA1),1)
        WRITE(LUPRI,*) 'Doubles part: ',
     &                 DDOT(NT2AMX,WORK(KRHOANA2),1,WORK(KRHOANA2),1)
        IF (CCR12) WRITE(LUPRI,*) 'R12-Doubles part: ',
     &             DDOT(NTR12AM(1),WORK(KRHOANA12),1,WORK(KRHOANA12),1)
        CALL CC_PRP(WORK(KRHOANA1),WORK(KRHOANA2),1,1,N)
        IF (CCR12) THEN
           CALL CC_PRPR12(WORK(KRHOANA12),1,1,.TRUE.)
        END IF
      END IF

      WRITE(LUPRI,*) 'LEAVING CC_FTSTNEW!'
      RETURN
      END
*=====================================================================*
      SUBROUTINE CC_FDF(RHODIF1,RHODIF2,RHODIF12,
     &                  LISTL,IDLSTL,LISTR,IDLSTR,
     &                  WORK,LWORK,APROXR12)
C
C---------------------------------------------------------------------*
C      Christian Neiss, Christof Haettig, october 2004
C      
C      test routine for F-matrix. calculates contractions of the 
C      F-matrix by numerical differentiation of the Jacobian.
C---------------------------------------------------------------------*
C
      IMPLICIT NONE  
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccfro.h>
#include <dummy.h>
#include <r12int.h>
#include <ccr12int.h>

* local parameters:
      LOGICAL LOCDBG, LTESTV
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER LWORK
#if defined (SYS_CRAY) 
      REAL WORK(LWORK) 
      REAL RHODIF1(*), RHODIF2(*), RHODIF12(*)
      REAL ECURR, FAC, DDOT
      REAL HALF, ZERO, ONE, TWO, DELTA, DELINV
#else 
      DOUBLE PRECISION WORK(LWORK) 
      DOUBLE PRECISION RHODIF1(*), RHODIF2(*), RHODIF12(*)
      DOUBLE PRECISION ECURR, FAC, DDOT
      DOUBLE PRECISION HALF, ZERO, ONE, TWO, DELTA, DELINV
#endif 

      CHARACTER*(3) LISTR, LISTL, APROXR12
      CHARACTER*(10) MODEL, MODELW
      INTEGER IDLSTL, IDLSTR, ISYM, IOPT, ISIDE
      INTEGER KT1AMPSAV, KT2AMPSAV, KT12AMPSAV, KT1AMP0, KT2AMP0,
     &        KT12AMP0, KT1AMPA, KT2AMPA, KT12AMPA, KOMEGA1, KOMEGA2
      INTEGER KZETA1, KZETA2, KZETA12, KT0AMP0, ISYM0, ISYMR, ISYML
      INTEGER KRHO1, KRHO2, KRHO12, KETA1, KETA2, KETA12, LUNIT
      INTEGER IDUM, LDUM, LENMOD
      INTEGER KEND1, LWRK1, KEND2, LWRK2
      INTEGER KLAMDP, KLAMDH, KVABKL, KVDIFF, NVCDKL, KVCDKL
      
      ! external functions:      
      INTEGER IR1TAMP
  

      PARAMETER (HALF=0.5D0, ZERO=0.0D0, ONE=1.0D0, TWO=2.0D0)
      
C---------------------------------------------------------------------*
C displacement parameter for finite difference:
C---------------------------------------------------------------------*
      PARAMETER (DELTA = 1.0D-07, DELINV = 1.0D+07)
C      PARAMETER (DELTA = 1.0D-02, DELINV = 1.0D+02)

      
      CALL QENTER('CC_FDF')

C---------------------------------------------------------------------*
C Initialisation
C---------------------------------------------------------------------*

      LTESTV = .FALSE..AND.CCR12.AND..NOT.(CCS.OR.CC2)

C calc. dimension:
      NVCDKL = 0
      DO I = 1, NSYM
        NVCDKL = NVCDKL + NVIR(1)*NRHFB(1)*(NVIR(1)*NRHFB(1)+1)/2
      END DO

C only C1-symmetry:
      if (nsym.gt.1) then
        write (lupri,*) 'CC_FDF must run in C1 symmetry!'
        write (lupri,*) 'Will ignore finite differences....'
        goto 1000 
      end if
      isym  = 1
      isym0 = 1
      isyml = 1
      isymr = 1

C Models:
      modelw = 'UNKNOWN   '
      if (CCS) modelw  = 'CCS       '
      if (CC2) modelw  = 'CC2       '
      if (CCSD) modelw = 'CCSD      '
      
      CALL CCSD_MODEL(MODELW,LENMOD,24,MODELW,10,APROXR12)

C memory allocation:
      kt1ampsav = 1
      kt2ampsav = kt1ampsav + nt1amx
      kend1     = kt2ampsav + nt2amx
      if (CCR12) then
        kt12ampsav = kend1
        kend1      = kt12ampsav + ntr12am(isym)
      end if

      kt1ampa = kend1
      kt2ampa = kt1ampa + nt1amx
      kend1   = kt2ampa + nt2amx
      if (CCR12) then
        kt12ampa = kend1
        kend1    = kt12ampa + ntr12am(isym)
C     dirty hack; otherwise problems with g77-compiler!
      else
        kt12ampa = lwork-1
      end if

      keta1  = kend1
      keta2  = keta1 + nt1amx
      kend1  = keta2 + nt2amx
      if (CCR12) then
        keta12 = kend1
        kend1  = keta12 + ntr12am(1)
      end if

      kzeta1 = kend1
      kzeta2 = kzeta1 + nt1amx
      kend1  = kzeta2 + nt2amx
      if (CCR12) then
        kzeta12 = kend1
        kend1   = kzeta12 + ntr12am(isym)
      end if

      kt0amp0 = kend1 
      kt1amp0 = kt0amp0 + 2*nallai(isym)
      komega1 = kt1amp0 + nt1amx
      komega2 = komega1 + nt1amx
      kt2amp0 = komega2 + max(nt2amx,2*nt2ort(1),nt2ao(1))
      kend1   = kt2amp0 + max(nt2sq(1),nt2r12(1))
      if (CCR12) then
        kt12amp0 = kend1
        kend1    = kt12amp0 + ntr12am(isym)
      end if

      if (ltestv) then
        kvdiff = kend1
        kvcdkl = kvdiff + nvcdkl
        kend1  = kvcdkl + nvcdkl
      end if

      krho1  = kend1
      krho2  = krho1 + nt1amx
      kend1  = krho2 + nt2amx
      if (CCR12) then
        krho12 = kend1
        kend1  = krho12 + ntr12am(isym)
      end if
      lwrk1 = lwork - kend1

C test if work space is enough:
      if (lwrk1 .lt. 0) then
        write(lupri,*) 'Not enough work memory in cc_FDF:'
        write(lupri,*) 'Available: LWORK = ',lwork
        write(lupri,*) 'Needed: ', kend1
        call quit('Too little work space in CC_FDF!')
      end if

C zero out result arrays:
      call dzero(rhodif1,nt1amx)
      call dzero(rhodif2,nt2amx)
      if (CCR12) call dzero(rhodif12,ntr12am(isym))

      if (ltestv) then
        call dzero(work(kvcdkl),nvcdkl)
        call dzero(work(kvdiff),nvcdkl)
      end if

C---------------------------------------------------------------------*
C Read in Lagrange multipliers
C---------------------------------------------------------------------*
      iopt=3
      call cc_rdrsp(listl,idlstl,isyml,iopt,model,
     &              work(kzeta1),work(kzeta2))
      if (locdbg) then
        write(lupri,*) 'norm of zeta1 at 0:',
     &    ddot(nt1amx,work(kzeta1),1,work(kzeta1),1)
      end if
      if (CCR12) then
       iopt=32
       call cc_rdrsp(listl,idlstl,isyml,iopt,model,
     &               dummy,work(kzeta12))
      end if
      
C---------------------------------------------------------------------*
C Read in cluster amplitudes
C---------------------------------------------------------------------*
      iopt=3
      call cc_rdrsp('R0 ',0,isym0,iopt,model,work(kt1ampsav),
     &               work(kt2ampsav))
      if (CCR12) then
        iopt=32
        call cc_rdrsp('R0 ',0,isym0,iopt,model,
     &                dummy,work(kt12ampsav))
      end if

C---------------------------------------------------------------------*
C Read in trial vector
C---------------------------------------------------------------------*
      iopt=3
      call cc_rdrsp(listr,idlstr,isymr,iopt,model,work(kt1ampa),
     &              work(kt2ampa))
C     trial vector diagonal elements must be scaled by 2 due to 
C     different format than ground state amplitudes:
C     (CCRHSN expects ground state format)
      call cclr_diascl(work(kt2ampa),two,1)
      if (CCR12) then
        iopt=32
        call cc_rdrsp(listr,idlstr,isymr,iopt,model,dummy,
     &                 work(kt12ampa))
C       trial vector diagonal elements must be scaled by KETSCL due to
C       different format than ground state amplitudes:
C       (CCRHSN expects ground state format)
        call cclr_diasclr12(work(kt12ampa),ketscl,1)
        if (locdbg) then
          write(lupri,*) 'T12AMPA in CC_FDF:'
          call outpak(work(kt12ampa),nmatki(1),1,lupri)
        end if
      end if
      if (locdbg) then 
        write(lupri,*) 'norm of zeta1 at 1:',
     &    ddot(nt1amx,work(kzeta1),1,work(kzeta1),1)
      end if

C---------------------------------------------------------------------*
C Loop for +/- displacement DELTA
C---------------------------------------------------------------------*
      do n = 1, 2

C---------------------------------------------------------------------*
C   add/subtract finite displacement Delta*Trialvector to t 
C   and recalculate intermediates
C---------------------------------------------------------------------*
        call dcopy(nt1amx,work(kt1ampsav),1,work(kt1amp0),1)
        call dcopy(nt2amx,work(kt2ampsav),1,work(kt2amp0),1)
        if (CCR12) then
          call dcopy(ntr12am(isym),work(kt12ampsav),1,work(kt12amp0),1)
        end if       
        
        if (n.eq.1) fac=one
        if (n.eq.2) fac=-one

        call daxpy(nt1amx,fac*delta,work(kt1ampa),1,work(kt1amp0),1)
        call daxpy(nt2amx,fac*delta,work(kt2ampa),1,work(kt2amp0),1)
        if (CCR12) then
          call daxpy(ntr12am(isym),fac*delta,work(kt12ampa),1,
     &               work(kt12amp0),1)
        end if

        iopt=3
        call cc_wrrsp('R0 ',0,isym0,iopt,modelw,work(kt0amp0),
     &              work(kt1amp0),work(kt2amp0),work(kend1),lwrk1)
        if (CCR12) then
          iopt = 32
          call cc_wrrsp('R0 ',0,isym0,iopt,modelw,work(kt0amp0),
     &                  dummy,work(kt12amp0),work(kend1),lwrk1)
        end if

        if (CCR12.AND..NOT.(CCS.OR.CC2)) then
          klamdp = kend1
          klamdh = klamdp + nlamdt
          kvabkl = klamdh + nlamdt
          kend2  = kvabkl + nvabkl(1)
          lwrk2  = lwork - kend2
          if (lwrk2.lt.0) then
            call quit('Insufficient work space for Vabkl in CC_FDF')
          end if

          call lammat(work(klamdp),work(klamdh),work(kt1amp0),
     &                work(kend2),lwrk2)

          lunit = -1
          call gpopen(lunit,fvabkl,'OLD',' ','UNFORMATTED',idummy,
     &                .FALSE.)
          read(lunit) (work(kvabkl+i-1), i=1, nvabkl(1))
          call gpclose(lunit,'KEEP')
       
          iopt = 1 
          call cc_r12mkvirt(work(kvabkl),work(klamdp),1,work(klamdp),1,
     &                      'R12VCTDTKL',iopt,work(kend2),lwrk2)
        end if

        rspim=.TRUE.
        call ccrhsn(work(komega1),work(komega2),
     &              work(kt1amp0),work(kt2amp0),
     &              work(kend1),lwrk1,aproxr12)

C---------------------------------------------------------------------*
C   Do left hand transformation with Lagrange multipliers 
C---------------------------------------------------------------------*
        if (locdbg) then 
          write(lupri,*) 'norm^2 of zeta1 at 2:',
     &      ddot(nt1amx,work(kzeta1),1,work(kzeta1),1)
        end if
        call dcopy(nt1amx,work(kzeta1),1,work(krho1),1)
        call dcopy(nt2amx,work(kzeta2),1,work(krho2),1)
        if (CCR12) then
          call dcopy(ntr12am(isym),work(kzeta12),1,work(krho12),1)
        end if

        ecurr=0.0D0 
        iside=-1
      
        call cc_atrr(ecurr,isym,iside,work(krho1),lwrk1,.FALSE.,dummy,
     &               aproxr12)

        call cc_eta(work(keta1),work(kend1),lwrk1)

C       we do not need to add doubles- and r12-contributions of eta
C       to rho, since these quantities are independent on amplitudes:
        call daxpy(nt1amx,one,work(keta1),1,work(krho1),1)
        
        if (locdbg) then
          write(lupri,*) 'norm of zeta1 at 3:',
     &      ddot(nt1amx,work(kzeta1),1,work(kzeta1),1)
        end if
C---------------------------------------------------------------------*
C   divide by 2*delta
C---------------------------------------------------------------------*
        call daxpy(nt1amx,fac*half*delinv,work(krho1),1,rhodif1,1)
        call daxpy(nt2amx,fac*half*delinv,work(krho2),1,rhodif2,1)
        if (CCR12) then
          call daxpy(ntr12am(isym),fac*half*delinv,work(krho12),1,
     &                rhodif12,1)
        end if

        if (ltestv) then
          lunit = -1
          call gpopen(lunit,'R12VCTDTKL','OLD',' ','UNFORMATTED',idummy,
     &                .FALSE.)
          read(lunit) (work(kvcdkl+i-1), i=1, nvcdkl)
          call gpclose(lunit,'KEEP')

          call daxpy(nvcdkl,fac*half*delinv,work(kvcdkl),1,
     &               work(kvdiff),1)

          if (n.eq.2) then
            call around('Finite Diff. Result of Vabkl in CC_FDF')
            call outpkb(work(kvdiff),nvir(1)*nrhfb(1),1,1,lupri)
          end if 
        end if

C---------------------------------------------------------------------*
C End of Loop over +/- DELTA
C---------------------------------------------------------------------*
      end do

C---------------------------------------------------------------------*
C restore original amplitudes
C---------------------------------------------------------------------*
      call dcopy(nt1amx,work(kt1ampsav),1,work(kt1amp0),1)
      call dcopy(nt2amx,work(kt2ampsav),1,work(kt2amp0),1)
      if (CCR12) then
        call dcopy(ntr12am(isym),work(kt12ampsav),1,work(kt12amp0),1)
      end if

      iopt=3
      call cc_wrrsp('R0 ',0,isym0,iopt,modelw,work(kt0amp0),
     &              work(kt1amp0),work(kt2amp0),work(kend1),lwrk1)
      if (CCR12) then
        iopt = 32
        call cc_wrrsp('R0 ',0,isym0,iopt,modelw,work(kt0amp0),
     &                dummy,work(kt12amp0),work(kend1),lwrk1)
      end if

C---------------------------------------------------------------------*
C restore original intermediates
C---------------------------------------------------------------------*
      if (CCR12.AND..NOT.(CCS.OR.CC2)) then
        klamdp = kend1
        klamdh = klamdp + nlamdt
        kvabkl = klamdh + nlamdt
        kend2  = kvabkl + nvabkl(1)
        lwrk2  = lwork - kend2
        if (lwrk2.lt.0) then
          call quit('Insufficient work space for Vabkl in CC_FDF')
        end if

        call lammat(work(klamdp),work(klamdh),work(kt1amp0),
     &              work(kend2),lwrk2)

        lunit = -1
        call gpopen(lunit,fvabkl,'OLD',' ','UNFORMATTED',idummy,
     &              .FALSE.)
        read(lunit) (work(kvabkl+i-1), i=1, nvabkl(1))
        call gpclose(lunit,'KEEP')

        iopt = 1
        call cc_r12mkvirt(work(kvabkl),work(klamdp),1,work(klamdp),1,
     &                    'R12VCTDTKL',iopt,work(kend2),lwrk2)
      end if

      rspim=.TRUE.
      call ccrhsn(work(komega1),work(komega2),
     &            work(kt1amp0),work(kt2amp0),
     &            work(kend1),lwrk1,aproxr12)

1000  continue
      CALL QEXIT('CC_FDF')

      RETURN
      END 
*=====================================================================*
