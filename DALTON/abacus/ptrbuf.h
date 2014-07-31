! FILE: ptrbuf.h -- only used in DALTON/abacus/abaptr.F
!                -- depends on: maxorb.h
      INTEGER         MAXADR, MAXCHN
      PARAMETER       (MAXADR = 200, MAXCHN = 200)

      INTEGER         LASTAD, IADR
      INTEGER         MX1BUF, L1BUF,  MEMS,   MEMT,  MXABUF,
     &                LABUF,  LDAMAX, MX2BUF, L2BUF, NCHAIN, NBLOCK
      COMMON /PTRBUF/
     &                LASTAD(MAXADR), IADR(MAXADR),
     &                MX1BUF, L1BUF,  MEMS,   MEMT,  MXABUF,
     &                LABUF,  LDAMAX, MX2BUF, L2BUF, NCHAIN, NBLOCK
! end of ptrbuf.h
