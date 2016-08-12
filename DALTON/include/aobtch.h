! FILE: aobtch.h
! Used for eri and densfit

      LOGICAL ACTVBT

      INTEGER AOBTCHlast

      COMMON /AOBTCH/ EXPBT(MXPRIM),                                     ! real*8
     &                CORXBT(MXSHEL), CORYBT(MXSHEL), CORZBT(MXSHEL),    ! integer
     &                NHKTBT(MXSHEL), KCKTBT(MXSHEL), KHKTBT(MXSHEL),
     &                NPRFBT(MXSHEL), NCTFBT(MXSHEL), ISTBBT(MXSHEL),
     &                MULTBT(MXSHEL), NCNTBT(MXSHEL), KCLSBT(MXSHEL),
     &                KNDXBT(MXSHEL), KCMTBT(MXSHEL),
     &                KEXPBT(MXSHEL), KCCFBT(MXSHEL),
     &                KAOSRT(MXSHEL), NORBBT(MXSHEL),
     &                NAOBCH, NBASE, NGAB, IORBRP(0:7),
     &                MAXQN, KQNBT(MXQN), NQNBT(MXQN),
     &                ACTVBT(MXSHEL,4),                                  ! logical
     &   AOBTCHlast !  Very important:
      !  Always keep AOBTCHlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
! -- end of aobtch.h
