C
C$Id: aobtch.h,v 1.1.1.1 2001-02-08 13:33:25 hjj Exp $
C
      LOGICAL ACTVBT
      COMMON /AOBTCH/ EXPBT(MXPRIM),  
     &                CORXBT(MXSHEL), CORYBT(MXSHEL), CORZBT(MXSHEL),
     &                NHKTBT(MXSHEL), KCKTBT(MXSHEL), KHKTBT(MXSHEL),
     &                NPRFBT(MXSHEL), NCTFBT(MXSHEL), ISTBBT(MXSHEL),
     &                MULTBT(MXSHEL), NCNTBT(MXSHEL),
     &                KNDXBT(MXSHEL),
     &                KEXPBT(MXSHEL), KCCFBT(MXSHEL),
     &                KAOSRT(MXSHEL),
     &                NORBBT(MXSHEL), ACTVBT(MXSHEL,4),
     &                NAOBCH, NBASIS, NGAB, IORBRP(0:7),
     &                MAXQN, KQNBT(MXQN), NQNBT(MXQN), NODD
