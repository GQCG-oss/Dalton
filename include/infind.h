      PARAMETER (LIROW  = 3*MXCORB + 3)
      PARAMETER (JTFRO  = 1, JTINAC = 2, JTACT  = 3, JTSEC  = 4)
      PARAMETER (JTFRFR =-1, JTINFR =-2, JTACFR =-3, JTSEFR =-4,
     &           JTININ = 1, JTACIN = 2, JTACAC = 3, JTSEIN = 4,
     &           JTSEAC = 5, JTSESE = 6)
      COMMON /INFIND/ IROW(LIROW),  ISMO(MXCORB),ISAO(MXCORB),
     &                ISW(MXCORB),  ISX(MXCORB),
     &                ICH(MXCORB),  LOC(MXCORB), IOBTYP(MXCORB),
     &                NSM(MAXASH),  IACTYP(MAXASH),
     &                ISSMO(MXCORB),ISSORD(MXCORB)
