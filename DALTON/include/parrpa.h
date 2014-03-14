C     PARRPA = Parallel Random Phase Approximation
C
C     LPARRP Is the number of variables in /parrpa/. Update this parameter 
C     if you change the number of variables in the common block. 
C     It's imperative for synchronization of slaves with the 
C     master in RPA and SOPPA calculations.
      integer, parameter :: LPARRP=42
C
C
      integer :: NUMDIS, KODCL1, KODCL2, KODBC1, 
     &           KODBC2, KRDBC1, KRDBC2, KODPP1,
     &           KODPP2, KRDPP1, KRDPP2, KCCFB1,
     &           KINDXB, NTOT,
     &           KCMO  , KEND1 , LWORK1,
     &           LTR1E , LTR1D , LRES1E, LRES1D,
     &           LFOCK , LDENS , LBTR1E, LBTR1D, 
     &           LBTJ1E, LBTJ1D, KTR1E , KTR1D ,
     &           KRES1E, KRES1D, KFOCK , KDENS ,
     &           KBTR1E, KBTR1D, KBTJ1E, KBTJ1D,
     &           KEND2,  LWORK2, KENDSV, LWORKSV,
     &           ISYMPAR
C
C
      common /parrpa/
     & NTOT  , NUMDIS, KODCL1, KODCL2, KODBC1, LWORKSV,
     & KODBC2, KRDBC1, KRDBC2, KODPP1, KODPP2, KRDPP1, 
     & KRDPP2, KCCFB1, KINDXB, KCMO  , KEND1 , LWORK1,
     & LTR1E , LTR1D , LRES1E, LRES1D, LFOCK , LDENS , 
     & LBTR1E, LBTR1D, LBTJ1E, LBTJ1D, KTR1E , KTR1D ,
     & KRES1E, KRES1D, KFOCK , KDENS , KBTR1E, KBTR1D, 
     & KBTJ1E, KBTJ1D, KEND2 , LWORK2, KENDSV, ISYMPAR
C
      save /parrpa/
