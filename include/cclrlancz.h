!  FILE: cclrlancz.h
!
      REAL*8  DAMPING(20), FREQ_RANGE(3), EIG_RANGE(2)
      INTEGER JCHAIN, IPRLRLCZ, NDAMP
      INTEGER JCHAINOLD,JCHAINNEW
      LOGICAL ABS_RANGE, LCHAINADD
      LOGICAL ABSANALYZE, DUMP_EIGFIL, DUMP_ALLFIL
      LOGICAL REDMEML,Debug_sym
      LOGICAL SUM_RULES
      CHARACTER LABELO*(8)
      COMMON /CCLRLANCZ/
     &                 DAMPING, FREQ_RANGE, EIG_RANGE,
     &                 JCHAIN, IPRLRLCZ,NDAMP,
     &                 JCHAINOLD,JCHAINNEW,
     &                 ABS_RANGE, LCHAINADD,
     &                 ABSANALYZE,DUMP_EIGFIL,DUMP_ALLFIL,
     &                 REDMEML,Debug_sym,
     &                 SUM_RULES,
     &                 LABELO
!  -- end of cclrlancz.h --
