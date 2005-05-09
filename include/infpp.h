      PARAMETER ( MAXLPP = MAXLBL )
c     LCMEXC(ISYM,IEXC) - whether transition moments involving IEXC:th
c     state in ISYM symmetry should be computed.
      LOGICAL LPPOP, EXMOM, EXMTES, LCMEXC
      CHARACTER*8 LBLPP
      COMMON /INFPP/  THCPP,  IPRPP,  MAXITP, IPREXM,  NPPCNV(8),
     *                NPPSIM(8), NPPSTV(8), NGPPP(8), EXMOM ,
     *                LPPOP(MAXLPP),EXMTES,LCMEXC(MAXLPP,8)
      COMMON /CHRPP/  LBLPP(8,MAXLPP)
