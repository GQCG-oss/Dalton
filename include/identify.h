C
C$Id: identify.h,v 1.1.1.1 2001-02-08 13:33:27 hjj Exp $
C
      CHARACTER*72 STARS,RELEAS
#if defined (VAR_MOTECC)
      CHARACTER*12 VERSN
#else
      CHARACTER*3  VERSN
#endif
      DATA STARS(1:36)  /'************************************'/
      DATA STARS(37:72) /'************************************'/
      DATA RELEAS(1:36) /'*SIRIUS* a direct, restricted step, '/
      DATA RELEAS(37:72)/'second order MCSCF program  *Apr 96*'/
#if defined (VAR_MOTECC)
      DATA VERSN        /'METECC-94.n08'/
#else
      DATA VERSN        /'n08'/
#endif
