      LOGICAL ACTVBT

      INTEGER AOBTCHLAST

      COMMON /AOBTCH/ EXPBT(MXPRIM),                                    &
     &                CORXBT(MXSHEL), CORYBT(MXSHEL), CORZBT(MXSHEL),   &
     &                NHKTBT(MXSHEL), KCKTBT(MXSHEL), KHKTBT(MXSHEL),   &
     &                NPRFBT(MXSHEL), NCTFBT(MXSHEL), ISTBBT(MXSHEL),   &
     &                MULTBT(MXSHEL), NCNTBT(MXSHEL), KCLSBT(MXSHEL),   &
     &                KNDXBT(MXSHEL), KCMTBT(MXSHEL),                   &
     &                KEXPBT(MXSHEL), KCCFBT(MXSHEL),                   &
     &                KAOSRT(MXSHEL),                                   &
     &                NORBBT(MXSHEL), ACTVBT(MXSHEL,4),                 &
     &                NAOBCH, NBASE, NGAB, IORBRP(0:7),                 & 
     &                MAXQN, KQNBT(MXQN), NQNBT(MXQN)

      COMMON /AOBTCH/ AOBTCHLAST
      !   Very important !!!
      !  Always keep ERIMEMlast as the last variable in the common block. 
      !  See GETBYTESPAN(firstvar, <name>last, SizeInBytes) for explanation.
