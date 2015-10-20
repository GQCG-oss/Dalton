C
C     Parallellization is sensitive to changes herein.
C
      INTEGER   MTAG0, MTAG1, MTAG2, MTAG3, MTAG4, MTAG5, MTAG6
      INTEGER   MTAG7, MTAG8
      INTEGER   MTAG51, MTAG52, MTAG53, MTAG54
      INTEGER   MTAG61, MTAG62, MTAG63, MTAG64
      PARAMETER (MTAG0 =  9, MTAG1 = 10, MTAG2 = 20, MTAG3 = 30,
     &     MTAG4 = 40, MTAG5 = 50, MTAG6 = 60, 
     &     MTAG7 = 70, MTAG8 = 80,
     &     MTAG51 = 51, MTAG52 = 52, MTAG53 = 53, MTAG54 = 54,
     &     MTAG61 = 61, MTAG62 = 62, MTAG63 = 63, MTAG64 = 64)
C
C     PCM-specific communication tags
C
      INTEGER    MPTAG1, MPTAG2
      PARAMETER (MPTAG1 = 1, MPTAG2 = 2)
