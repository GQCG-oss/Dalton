!
! FILE: wrkrsp.h
!
! Contains info for current "work" in response solver
! while infrsp.h contains more permanent information.
!
! Many values in wrkrsp.h are set in RSPVAR in "rsp/rspmai.F";
! RSPVAR must be called each time the operator symmetry KSYMOP or
! the operator frequency AFREQ change in the RESPONS module in "rsp/".
!
      REAL*8          AFREQ
      INTEGER         KSYMOP, KZVAR,  KZYVAR, KZCONF, KZYCON, KZWOPT,   &
     &                KZYWOP, KZRED , KZYRED, KEXCNV, KEXSIM, KEXSTV,   &
     &                KLRSTV, JEXSIM, KOFFTY, KCONV,  KSYMST
      LOGICAL         RESTPP, ABCHK,  ABSYM,  RESTLR, RESTC6, DETERM
      COMMON /WRKRSP/ AFREQ,                                            &
     &                KSYMOP, KZVAR,  KZYVAR, KZCONF, KZYCON, KZWOPT,   &
     &                KZYWOP, KZRED , KZYRED, KEXCNV, KEXSIM, KEXSTV,   &
     &                KLRSTV, JEXSIM, KOFFTY, KCONV,  KSYMST,           &
     &                RESTPP, ABCHK,  ABSYM,  RESTLR, RESTC6, DETERM
! -- end of wrkrsp.h --
