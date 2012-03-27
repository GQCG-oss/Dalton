      INTEGER JCHAIN, IPRLRLCZ, NDAMP
      INTEGER JCHAINOLD,JCHAINNEW
      LOGICAL ABS_RANGE, LCHAINADD
      LOGICAL ABSANALYZE, DUMP_EIGFIL, DUMP_ALLFIL
      LOGICAL REDMEML,Debug_sym
      LOGICAL SUM_RULES
      CHARACTER LABELO*(8)
#if defined (SYS_CRAY)
      REAL DAMPING(20), FREQ_RANGE(3), EIG_RANGE(2)
#else
      DOUBLE PRECISION DAMPING(20),FREQ_RANGE(3),EIG_RANGE(2)
#endif
      COMMON /CCLRLANCZ/ ABS_RANGE, 
     &                 LCHAINADD,
     &                 JCHAIN, IPRLRLCZ,NDAMP,
     &                 JCHAINOLD,JCHAINNEW,
     &                 LABELO,
     &                 ABSANALYZE,DUMP_EIGFIL,DUMP_ALLFIL,
     &                 SUM_RULES,
     &                 REDMEML,Debug_sym,
     &                 DAMPING, FREQ_RANGE, EIG_RANGE
