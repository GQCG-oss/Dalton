C
C$Id: aobtch.h,v 1.2 2001-02-12 18:17:55 vebjornb Exp $
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
     &                NAOBCH, NBASE, NGAB, IORBRP(0:7),
     &                MAXQN, KQNBT(MXQN), NQNBT(MXQN)
