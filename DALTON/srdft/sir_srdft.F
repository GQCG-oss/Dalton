!
!...   Copyright (c) 2014 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2014 (2014), see http://daltonprogram.org"
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
C File: sir_srdft.F
C
C Special rouitnes for srDFT in sirius/
C Responsible author: Hans Joergen Aa. Jensen.


C*****************************************************************************
      SUBROUTINE SRDFT_DENS(DAO,UDV,PV,CMO,WRK,KFREE,LFREE)
C*****************************************************************************
C      Hans Joergen Aa. Jensen, 2-Jun-2009
C
C      Purpose: calculate DCAO and DVAO matrices for srDFT
C               (done in a new routine for future flexibility
C                for e.g. also  on-top pair density)
C
C      July 2017 hjaaj: If (DFT_SPINDNS) then UDV and DVAO are doubled,
C      containing both "DV" active charge density matrix
C      and "DS" active spin density matrix. This is done
C      inside GTDMSO, thus we do not need to check for DFT_SPINDNS here.
C
C*****************************************************************************
#include "implicit.h"
C
#include "inforb.h"
      DIMENSION DAO(N2BASX,*), UDV(NASHT,NASHT), PV(*), CMO(*), WRK(*)
C
#include "dftcom.h"
C
      CALL QENTER('SRDFT_DENS')
      IF (DFT_SPINDNS) THEN
         call quit('DFT_SPINDNS not implemented in SRDFT_DENS yet')
      END IF
      IF (LFREE .LT. N2BASX) CALL QUIT('Insufficient WRK dim')
      CALL GTDMSO(UDV,CMO,DAO(1,1),DAO(1,2),WRK(KFREE))
      IF (NASHT .GT. 0) THEN
          CALL DAXPY(N2BASX,1.0D0,DAO(1,2),1,DAO(1,1),1)
      END IF
      CALL QEXIT ('SRDFT_DENS')
      RETURN
      END
C
      SUBROUTINE FCK_MAKE_FCEFF_LOCALSPIN_MODEL1(FC,FS,DV,WRK,LWRK)
C
C     Feb/Oct 2010 Hans Joergen Aa. Jensen
C
!>    DS   = f(DV) DV
!>    DS'  = f'(DV) DV + f(DV) DV'
!>    Model 1:
!>      f(DV)  = (2 - DV) DV
!>      f'(DV) = - DV' DV + (2 - DV) DV' = 2 DV' - DV' DV - DV DV'
C
C
#include "implicit.h"
#include "priunit.h"
C
#include "inforb.h"
C
      DIMENSION FC(*), FS(*), DV(*), WRK(LWRK)
C
      KFRSAV= 1
      KFREE = KFRSAV
      LFREE = LWRK
      CALL MEMGET('REAL',KUFC,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUFS,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDV,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP,N2ORBX,WRK,KFREE,LFREE)
C
C     Unpack FS and DV
C
      IF (NSYM.GT.1) THEN
         CALL PKSYM1(WRK(KTMP),FS,NORB,NSYM,-1)
         CALL DSPTGE(NORBT,WRK(KTMP),WRK(KUFS))
      ELSE
         CALL DSPTGE(NORBT,FS,WRK(KUFS))
      END IF
      CALL DSPTGE(NASHT,DV,WRK(KUDV))

#ifdef VAR_DEBUG
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: FS packed'
      call outpkb(fs,norb,nsym,-1,lupri)
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: FS unpacked'
      call output(wrk(kufs),1,norbt,1,norbt,norbt,norbt,-1,lupri)
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: DV packed'
      call outpAK(DV,NASHT,-1,lupri)
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: DS packed'
      call outpAK(DV(1+NNASHX),NASHT,-1,lupri)
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: DV unpacked'
      call output(wrk(kudv),1,nasht,1,nasht,nasht,nasht,-1,lupri)
#endif
C
C     calculate sum_r ( Dtot(p,r) FS(r,q) ) into tmp(p,q)
C     and add tmp(transposed) to tmp
      CALL FCK_MAKE_DxFS(WRK(KUFS),WRK(KUDV),WRK(KTMP))
#ifdef VAR_DEBUG
      write (lupri,*) ' - D FS - FS D unpacked'
      call output(wrk(ktmp),1,norbt,1,norbt,norbt,norbt,-1,lupri)
#endif
C
      CALL DAXPY(N2ORBX,2.0D0,WRK(KUFS),1,WRK(KTMP),1)
#ifdef VAR_DEBUG
      write (lupri,*) ' 2 FS - D FS - FS D unpacked'
      call output(wrk(ktmp),1,norbt,1,norbt,norbt,norbt,-1,lupri)
#endif
C
      CALL FCK_ESRDFTY_LOCALSPIN(WRK(KTMP),WRK(KUDV))
C     pack the correction matrix
      CALL DGETSP(NORBT,WRK(KTMP),WRK(KUFS))
      IF (NSYM.GT.1) CALL PKSYM1(WRK(KUFS),WRK(KUFS),NORB,NSYM,2)
#ifdef VAR_DEBUG
      write (lupri,*) ' 2 FS - D FS - FS D packed'
      call outpkb(wrk(kufs),norb,nsym,-1,lupri)
#endif
C     Add the effective FS (spin-density Fock matrix) contribution
C     to the FC matrix for the effective DFT_LOCALSPIN FC matrix.
      CALL DAXPY(NNORBT,1.0D0,WRK(KUFS),1,FC,1)
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: FC packed'
      call outpkb(fc,norb,nsym,-1,lupri)
      CALL MEMREL('FCK_MAKE_FCEFF_LOCALSPIN',WRK,KFRSAV,KFRSAV,
     &            KFREE,LFREE)
C      CALL MEMREL('FCK_..._LOCALSPIN',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END


C*****************************************************************************
      SUBROUTINE FCK_MAKE_FCEFF_LOCALSPIN(FC,FS,DV,WRK,LWRK)
C*****************************************************************************
C
C     Nov 2011 HJJ+MNP                     
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "inforb.h"
C
      DIMENSION FC(*), FS(*), DV(*), WRK(LWRK)
C
      KFRSAV= 1
      KFREE = KFRSAV
      LFREE = LWRK
      CALL MEMGET('REAL',KUFC,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUFS,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDV,N2ASHX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDT,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMP,N2ORBX,WRK,KFREE,LFREE)
C
C     Unpack FS and DV
C
      IF (NSYM.GT.1) THEN
         CALL PKSYM1(WRK(KTMP),FS,NORB,NSYM,-1)
         CALL DSPTGE(NORBT,WRK(KTMP),WRK(KUFS))
      ELSE
         CALL DSPTGE(NORBT,FS,WRK(KUFS))
      END IF
      CALL DSPTGE(NASHT,DV,WRK(KUDV))

#ifdef VAR_DEBUG
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: FS packed'
      call outpkb(fs,norb,nsym,-1,lupri)
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: FS unpacked'
      call output(wrk(kufs),1,norbt,1,norbt,norbt,norbt,-1,lupri)
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: DV packed'
      call outpAK(DV,NASHT,-1,lupri)
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: DS packed'
      call outpAK(DV(1+NNASHX),NASHT,-1,lupri)
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: DV unpacked'
      call output(wrk(kudv),1,nasht,1,nasht,nasht,nasht,-1,lupri)
#endif
!     Calculate Dtot:
      call FCK_MAKE_UDT(wrk(kudt),wrk(kudv))
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: DT' 
      call output(wrk(kudt),1,norbt,1,norbt,norbt,norbt,-1,lupri)
!     Calculate P = (2 - Dtot)
      CALL MEMGET('REAL',KUDP ,N2ORBX,WRK,KFREE,LFREE)
      CALL FCK_MAKE_UDP(wrk(kudt),wrk(kudp))
#ifdef VAR_DEBUG
      write (lupri,*) ' FCK_MAKE_FCEFF_LOCALSPIN: P = (2-DT)' 
      call output(wrk(kudp),1,norbt,1,norbt,norbt,norbt,-1,lupri)
#endif
C     Calculate the four terms in the effective FS matrix
!     First term: D'*D*P = DT*P**2*FS 
      CALL MEMGET('REAL',KUDPI,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDP2FS,N2ORBX,WRK,KFREE,LFREE)
!     P**2 = PI
      CALL DGEMM('N','N',norbt,norbt,norbt,1.0D0,
     &      WRK(KUDP),norbt,
     &      WRK(KUDP),norbt, 0.0D0,
     &      WRK(KUDPI),norbt)
!     P**2 * FS = P2FS
      CALL DGEMM('N','N',norbt,norbt,norbt,1.0D0,
     &      WRK(KUDPI),norbt,
     &      WRK(KUFS),norbt, 0.0D0,
     &      WRK(KUDP2FS),norbt)
!     DT * P2FS
      CALL DGEMM('N','N',norbt,norbt,norbt,1.0D0,
     &      WRK(KUDT),norbt,
     &      WRK(KUDP2FS),norbt, 0.0D0,
     &      WRK(KTMP),norbt)
!     end first term
!     second term = D*D'*P**2 = P**2*FS*DT
!     FS * DT = FSDT
      CALL MEMGET('REAL',KUFSDT,N2ORBX,WRK,KFREE,LFREE)
      CALL DGEMM('N','N',norbt,norbt,norbt,1.0D0,
     &      WRK(KUFS),norbt,
     &      WRK(KUDT),norbt, 0.0D0,
     &      WRK(KUFSDT),norbt)
!     PI * FSDT = P**2 * FS * DT
      CALL DGEMM('N','N',norbt,norbt,norbt,1.0D0,
     &      WRK(KUDPI),norbt,
     &      WRK(KUFSDT),norbt, 1.0D0,
     &      WRK(KTMP),norbt)
!     end second term
!     Third term = - D**2 * P * D' = - FS * DT**2 * P
!     DT**2 =DT2
      CALL MEMGET('REAL',KUDT2,N2ORBX,WRK,KFREE,LFREE)
      CALL MEMGET('REAL',KUDT2P,N2ORBX,WRK,KFREE,LFREE)
      CALL DGEMM('N','N',norbt,norbt,norbt,1.0D0,
     &      WRK(KUDT),norbt,
     &      WRK(KUDT),norbt, 0.0D0,
     &      WRK(KUDT2),norbt)
!     DT2 * P = DT2P
      CALL DGEMM('N','N',norbt,norbt,norbt,1.0D0,
     &      WRK(KUDT2),norbt,
     &      WRK(KUDP),norbt, 0.0D0,
     &      WRK(KUDT2P),norbt)
!     - FS * DT2P
      CALL DGEMM('N','N',norbt,norbt,norbt,-1.0D0,
     &      WRK(KUFS),norbt,
     &      WRK(KUDT2P),norbt, 1.0D0,
     &      WRK(KTMP),norbt)
!     end third term
!     Fourth term = - D**2 * D' * P = - P * FS * DT**2
!     FS * DT**2 = FSDT2
      CALL MEMGET('REAL',KUFSDT2,N2ORBX,WRK,KFREE,LFREE)
      CALL DGEMM('N','N',norbt,norbt,norbt,1.0D0,
     &      WRK(KUFS),norbt,
     &      WRK(KUDT2),norbt, 0.0D0,
     &      WRK(KUFSDT2),norbt)
!     - P * FSDT2
      CALL DGEMM('N','N',norbt,norbt,norbt,-1.0D0,
     &      WRK(KUDP),norbt,
     &      WRK(KUFSDT2),norbt, 1.0D0,
     &      WRK(KTMP),norbt)
!     end fourth term
C
      CALL FCK_ESRDFTY_LOCALSPIN(WRK(KTMP),WRK(KUDV))
C     pack the correction matrix
      CALL DGETSP(NORBT,WRK(KTMP),WRK(KUFS))
      IF (NSYM.GT.1) CALL PKSYM1(WRK(KUFS),WRK(KUFS),NORB,NSYM,2)
#ifdef VAR_DEBUG
      write (lupri,*) 'Corrected FS packed' 
      call outpkb(wrk(kufs),norb,nsym,-1,lupri)
#endif
C     Add the effective FS (spin-density Fock matrix) contribution
C     to the FC matrix for the effective DFT_LOCALSPIN FC matrix.
      CALL DAXPY(NNORBT,1.0D0,WRK(KUFS),1,FC,1)
      CALL MEMREL('FCK_MAKE_FCEFF_LOCALSPIN',WRK,KFRSAV,KFRSAV,
     &            KFREE,LFREE)
C      CALL MEMREL('FCK_..._LOCALSPIN',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C
      SUBROUTINE FCK_MAKE_UDP(UDT,TMP)
C
C     Nov 2011 mnp   
C     Utility routine called by FCK_MAKE_FCEFF_LOCALSPIN
C     to setup P =  2-DV  
C
#include "implicit.h"
#include "priunit.h"
C
#include "inforb.h"
      dimension udv(norbt,norbt), tmp(norbt,norbt)

      call dzero(tmp,n2orbx)
      do i = 1, norbt
         tmp(i,i) = 2.0d0
      end do
      call daxpy(n2orbx,-1.0d0,udt,1,tmp,1)
      return
      end


C*****************************************************************************
      SUBROUTINE FCK_MAKE_UDT(UDT,UDV)
C*****************************************************************************
C
C     Nov 2011 hjaaj
C     Utility routine to calculate UDT(NORBT,NORBT)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "inforb.h"
      dimension UDT(norbt,norbt), udv(nasht,nasht)
C
      call dzero(UDT,n2orbx)
      do isym = 1,nsym
         iorbi = iorb(isym)
         do ir = iorbi+1,iorbi+nish(isym)
            UDT(IR,IR) = 2.0D0
         end do
         iashi    = iash(isym)
         iash_off = iorbi + nish(isym)
         do irx = 1,nash(isym)
            do ipx = 1,nash(isym)
               UDT(iash_off+ipx,iash_off+irx) =
     &              UDV(iashi+ipx,iashi+irx)
            end do ! ipx
         end do ! irx
      end do ! isym
      RETURN
      END


C*****************************************************************************
      SUBROUTINE FCK_MAKE_DxFS(UFS,UDV,TMP)
C*****************************************************************************
C
C     Feb 2010 hjaaj
C     Utility routine called by FCK_MAKE_FCEFF_LOCALSPIN
C     to calculate "- sum_r ( Dtot(p,r) FS(r,q) )" into tmp(p,q)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "inforb.h"
      dimension UFS(norbt,norbt), udv(nasht,nasht), tmp(norbt,norbt)
C
      call dzero(tmp,n2orbx)
      do isym = 1,nsym
         iorbi = iorb(isym)
         iashi = iash(isym)
         do iq = iorbi+1,iorbi+norb(isym)
            do ir = iorbi+1,iorbi+nish(isym)
               TMP(IR,IQ) = -2.0D0*UFS(IR,IQ)
            end do ! ir
            iash_off = iorbi + nish(isym)
            do irx = iashi+1,iashi+nash(isym)
               FSRQ = UFS(iash_off+irx,iq)
               do ipx = iashi+1,iashi+nash(isym)
                  TMP(iash_off+ipx,IQ) = TMP(iash_off+ipx,IQ)
     &               - UDV(ipx,irx) * FSRQ
               end do ! ipx
            end do ! irx
         end do ! iq
      end do ! isym
C calculate TMP = TMP + TMP(transposed)
      do iq = 1,norbt
         do ip = 1,iq
            X = TMP(IP,IQ) + TMP(IQ,IP)
            TMP(IP,IQ) = X
            TMP(IQ,IP) = X
         end do
      end do
      RETURN
      END
      SUBROUTINE FCK_ESRDFTY_LOCALSPIN(UFC_LS,UDV)
C
C     Feb 2010 hjaaj
C     Utility routine called by FCK_MAKE_FCEFF_LOCALSPIN
C     to calculate contribution to ESRDFTY from
C     local spin contribution to EVSR:
C           EVSR = -Tr(Vxcsr Dv)
C     The "regular" contributions are calculated in SIRFCK,
C     but the local spin part was not (and could not be) calculated in
C     SIRFCK.
C
C     OBS! No local spin correction to this term
C           ECSR = -0.5*Tr(Vxcsr Dc)
C     because the local spin FC matrix wasn't added to Vxcsr when
C     EMCMY was calculated in FCKMAO.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxash.h"
C
#include "inforb.h"
#include "infind.h"
#include "dfterg.h"
      dimension UFC_LS(norbt,norbt), udv(nasht,nasht)
!     write (lupri,*) 'UFC_LS matrix, MO basis'
!     call output(ufc_ls,1,norbt,1,norbt,norbt,norbt,1,lupri)
!     write (lupri,*) 'UDV matrix, MO basis'
!     call output(udv,1,nasht,1,nasht,nasht,nasht,1,lupri)
!     write (lupri,*) 'ICH array'
!     write (lupri,'(10I5)') (ich(j),j=1,norbt)

      EVSR_LS = 0.0D0
      DO j = 1,norbt
         ichj = ich(j)
         if (ichj .gt. 0) then
            do i = 1,j-1
               if (ich(i) .gt. 0) then
                  EVSR_LS = EVSR_LS - 2.0D0*UDV(ich(i),ichj)*UFC_LS(i,j)
               end if
            end do
            EVSR_LS = EVSR_LS - UDV(ichj,ichj)*UFC_LS(j,j)
         end if
      end do

      ESRDFTY = ESRDFTY + EVSR_LS
      WRITE (LUPRI,'(/A)') '* Local spin corrections:'
      WRITE (LUPRI,'(A,F20.10)')
     &         ' -     Tr(Vxcsr_ls Dval)          ',EVSR_LS,
     &         '   ESRDFTY revised                ',ESRDFTY

      RETURN
      END
