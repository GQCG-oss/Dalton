
! Define some functions that can be called from 
! the old dft code without f90 complications

subroutine xcfun_select(conf_string)
  use xcfun_dalton
  character(*) conf_string
  integer :: i
  i = xcfun_select_by_name(conf_string)
end subroutine xcfun_select

! this requires that xcfun has been set up with xc_eval_set,
! use the prepare subroutines
SUBROUTINE xcfun_ksmlda(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO, &
     RHOA,GRADA,COORD,WGHT,FKSM)
!C
!C     P. Salek and T. Helgaker oct 2003
!C
  use xcfun_dalton
#include <implicit.h>
  PARAMETER (D2 = 2.0D0)
#include <inforb.h>
  DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),&
       RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8),&
       FKSM(*)
  EXTERNAL DFTENE
  COMMON /DKSMPRIV/ ENERGY_DFTKSM
  
  DIMENSION VXC(NBLEN), VX(5), DENS(1,NBLEN), OUT(2,NBLEN)
  PARAMETER (DUMMY = 0D0)
  DO IPNT = 1, NBLEN
     DENS(1,IPNT) = RHOA(IPNT) ! CHECK: xcfun wants total density, correct??
  ENDDO
  print *,'block starts, points = ',nblen
  DO IPNT = 1, NBLEN
     CALL DFTPTF0(RHOA(IPNT),DUMMY,WGHT(IPNT),VX)
     VXC(IPNT) = 2D0*VX(1)
     ENERGY_DFTKSM = ENERGY_DFTKSM &
             + DFTENE(RHOA(IPNT),DUMMY)*WGHT(IPNT)
     print *,'rhoa',RHOA(IPNT)
     print *,'old dft vx',VX
     print *,'old dft en',DFTENE(RHOA(IPNT),DUMMY)*WGHT(IPNT)
  END DO

  iorder = 1
  ires = xc_eval_setup(funobj,XC_N,XC_PARTIAL_DERIVATIVES,iorder)
  if (ires.ne.0) then
     print *,'a xc_eval_setup failed, this should not happen',ires
     stop
  endif
  call xc_eval(funobj,nblen,dens,out)
  DO IPNT = 1, NBLEN
!     ENERGY_DFTKSM = ENERGY_DFTKSM + WGHT(IPNT)*OUT(1,IPNT)
!     VXC(IPNT) = WGHT(IPNT)*OUT(2,IPNT)
     print *,'ulf weights = ',WGHT(IPNT)
     print *,'ulf VXC(IPNT) = ',WGHT(IPNT)*OUT(2,IPNT)
     print *,'ulf ene = ',WGHT(IPNT)*OUT(1,IPNT)
  END DO
  CALL DISTLDAB(NBLEN,NBLCNT,NBLOCKS,LDAIB,1,NBLEN, &
       VXC,GAO,FKSM)
  RETURN
END SUBROUTINE xcfun_ksmlda

#if 0
SUBROUTINE xcfun_ksmgga(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO, &
     RHOA,GRADA,COORD,WGHT,FKSM)
!C
!C     P. Salek and T. Helgaker oct 2003
!C
  use xcfun_dalton
#include <implicit.h>
  PARAMETER (D2 = 2.0D0)
#include <inforb.h>
  DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN),&
       RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8),&
       FKSM(*)
  EXTERNAL DFTENE
  COMMON /DKSMPRIV/ ENERGY_DFTKSM
  
  DIMENSION VXC(2,NBLEN), VX(5), DENS(2,NBLEN), OUT(3,NBLEN)
  PARAMETER (DUMMY = 0D0)
  DO IPNT = 1, NBLEN
     DENS(1,IPNT) = RHOA(IPNT) ! xcfun wants total density
     DENS(2,IPNT) = & ! gradient squared, different from Salek routines
          1.0D0*(GRADA(1,IPNT)**2 + &
          GRADA(2,IPNT)**2 + &
          GRADA(3,IPNT)**2)
  ENDDO
  iorder = 1
  ires = xc_eval_setup(funobj,XC_N_GNN,XC_PARTIAL_DERIVATIVES,iorder)
  if (ires.ne.0) then
     print *,'b xc_eval_setup failed, this shoulad not happen',ires
     stop
  endif
!#if 0
  DO IPNT = 1, NBLEN
     GRD = SQRT(GRADA(1,IPNT)**2+GRADA(2,IPNT)**2+GRADA(3,IPNT)**2)
     CALL DFTPTF0(RHOA(IPNT),GRD,WGHT(IPNT),VX)
     ENERGY_DFTKSM = ENERGY_DFTKSM &
          + DFTENE(RHOA(IPNT),GRD)*WGHT(IPNT)
     VXC(1,IPNT) = VX(1)
     IF(GRD.GT.1D-40) THEN
        VXC(2,IPNT) = VX(2)/(0.5*GRD)
     ELSE
        VXC(2,IPNT) = 0D0
     END IF
  END DO
!#endif
#if 0
  call xc_eval(funobj,nblen,dens,out)
  DO IPNT = 1, NBLEN
     ENERGY_DFTKSM = ENERGY_DFTKSM + OUT(1,IPNT)*WGHT(IPNT)
!     print *,'n, gnn, w',dens(1,ipnt),dens(2,ipnt), wght(ipnt)
!     print *,'lda contrib (xcfun, old)',OUT(2,IPNT)*WGHT(IPNT),VXC(1,IPNT)
!     print *,'gga contrib (xcfun, old)',OUT(3,IPNT)*WGHT(IPNT)*4.0D0 &
!          ,VXC(2,IPNT)
     VXC(1,IPNT) = OUT(2,IPNT)*WGHT(IPNT)
     VXC(2,IPNT) = OUT(3,IPNT)*WGHT(IPNT)*4.0D0
! todo: gradient contribution
  END DO
#endif
  CALL DISTGGAB(NBLEN,NBLCNT,NBLOCKS,LDAIB,1,NBLEN, &
       VXC,GAO,FKSM)
  RETURN
END SUBROUTINE xcfun_ksmgga
#endif

SUBROUTINE xcfun_KSMGGA(NBLEN,NBLCNT,NBLOCKS,LDAIB,GAO, &
     RHOA,GRADA,COORD,WGHT,FKSM)
  !C     P. Salek and T. Helgaker oct 2003
  use xcfun_dalton
#include <implicit.h>
#include <inforb.h>
  DIMENSION GAO(NBLEN,NBAST,*), COORD(3,NBLEN),WGHT(NBLEN), &
       RHOA(NBLEN), GRADA(3,NBLEN),NBLCNT(8),NBLOCKS(2,LDAIB,8), &
       FKSM(*)
  EXTERNAL DFTENE
  
  DIMENSION VXC(2,NBLEN),VX(5), DENS(2,NBLEN), OUT(3,NBLEN)
  COMMON /DKSMPRIV/ ENERGY_DFTKSM
  !C
  !C     Exchange-correlation contribution to Kohn-Sham matrix
  !C
  DO IPNT = 1, NBLEN
     DENS(1,IPNT) = RHOA(IPNT) ! xcfun wants total density
     DENS(2,IPNT) = & ! gradient squared, different from Salek routines
          1.0D0*(GRADA(1,IPNT)**2 + &
          GRADA(2,IPNT)**2 + &
          GRADA(3,IPNT)**2)
  ENDDO
  iorder = 1
  ires = xc_eval_setup(funobj,XC_N_GNN,XC_PARTIAL_DERIVATIVES,iorder)
  if (ires.ne.0) then
     print *,'b xc_eval_setup failed, this shoulad not happen',ires
     stop
  endif
#if 0
  DO IPNT = 1, NBLEN
     GRD = SQRT(GRADA(1,IPNT)**2+GRADA(2,IPNT)**2+GRADA(3,IPNT)**2)
     CALL DFTPTF0(RHOA(IPNT),GRD,WGHT(IPNT),VX)
     ENERGY_DFTKSM = ENERGY_DFTKSM &
          + DFTENE(RHOA(IPNT),GRD)*WGHT(IPNT)
     VXC(1,IPNT) = VX(1)
     IF(GRD.GT.1D-40) THEN
        VXC(2,IPNT) = VX(2)/(0.5*GRD)
     ELSE
        VXC(2,IPNT) = 0D0
     END IF
  END DO
#endif
  call xc_eval(funobj,nblen,dens,out)
  DO IPNT = 1, NBLEN
     ENERGY_DFTKSM = ENERGY_DFTKSM + OUT(1,IPNT)*WGHT(IPNT)
!     print *,'n, gnn, w',dens(1,ipnt),dens(2,ipnt), wght(ipnt)
!     print *,'lda contrib (xcfun, old)',OUT(2,IPNT)*WGHT(IPNT),VXC(1,IPNT)
!     print *,'gga contrib (xcfun, old)',OUT(3,IPNT)*WGHT(IPNT)*4.0D0 &
!          ,VXC(2,IPNT)
!     GRD = SQRT(GRADA(1,IPNT)**2+GRADA(2,IPNT)**2+GRADA(3,IPNT)**2)
!     print *,'ene ', DFTENE(RHOA(IPNT),GRD), OUT(1,IPNT)
     VXC(1,IPNT) = OUT(2,IPNT)*WGHT(IPNT)
     VXC(2,IPNT) = OUT(3,IPNT)*WGHT(IPNT)*4.0D0
  END DO
  CALL DISTGGAB(NBLEN,NBLCNT,NBLOCKS,LDAIB,1,NBLEN, &
       VXC,GAO,GRADA,FKSM)
END SUBROUTINE xcfun_KSMGGA
