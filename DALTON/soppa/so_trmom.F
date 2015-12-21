C
C  /* Deck so_trmom */
      SUBROUTINE SO_TRMOM(MODEL,ISYMTR,NLBTOT,NEXCI,
     &                    T2AM,LT2AM,DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,
CSPAS:23/5-11: second and third moment sum rules
C    &                    TRLEN,TRVEL,TRLON,TRMAG,BSRLON,EXENG,
C    &                    WORK,LWORK)
     &                    TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,TRLON,TRMAG,
     &                    BSRLON,EXENG,WORK,LWORK)
CKeinSPASmehr
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate RPA, RPA(D), SOPPA and SOPPA(CCSD) transition 
C              moments.
C
      use so_info, only: FN_RDENSE, FN_RDENSD
C
#include "implicit.h"
#include "priunit.h"
C
#include "cbiexc.h"
CSPAS:15/03-2006: merge with Dalton-2.0  
C#include "cbilrs.h"
CKeinSPASmehr
#include "ccsdsym.h"
#include "ccorb.h"
#include "soppinf.h"
#include "mxcent.h"
#include "nuclei.h"
C
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION TRLEN(3,NSYM,MXNEXI), TRVEL(3,NSYM,MXNEXI)
CSPAS:23/5-11: second and third moment sum rules
      DIMENSION TQLEN(3,3,NSYM,MXNEXI), TQVEL(3,3,NSYM,MXNEXI)
      DIMENSION TTLEN(10,NSYM,MXNEXI)
CKeinSPASmehr
      DIMENSION TRLON(3,NSYM,MXNEXI), TRMAG(3,NSYM,MXNEXI)
      DIMENSION BSRLON(3,NSYM,MXNEXI), EXENG(NSYM,MXNEXI)
      DIMENSION T2AM(LT2AM), WORK(LWORK)
C
      PARAMETER( ONE=1.0D0 )
C      
      CHARACTER*8 LABEL, PDENS_LABEL
      CHARACTER*5 MODEL
      CHARACTER*8 RTNLBL(2)
      LOGICAL  IMAGPROP,DOUBLES

C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TRMOM')
C
C------------------------
C     Do we need doubles?
C------------------------
C
      DOUBLES = MODEL.EQ.'AOSOP'.OR.MODEL.EQ.'AOSOC'.OR.
     &          MODEL.EQ.'DCRPA'
C
      WRITE(PDENS_LABEL,'(A7,I1)') 'EXCITA ', ISYMTR
      LURDENSE = -1
      LURDENSD = -1
      CALL GPOPEN(LURDENSE,FN_RDENSE,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.) 
      CALL GPOPEN(LURDENSD,FN_RDENSD,'UNKNOWN',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.) 
C      
C-------------------------------
C     Loop over property labels.
C-------------------------------
C
      DO 200 IPRLBL = 1, NLBTOT
C
         LABEL  = LABAPP(IPRLBL)
         KSYMOP = LABSYM(IPRLBL)
C
C-----------------------------------------------------------------
C        Check if integrals are second moment of charge integrals.
C-----------------------------------------------------------------
C
CSPAS:23/5-11: second and third moment sum rules
C        IF ( LABEL(3:8) .EQ. 'SECMOM' ) GOTO 200
CKeinSPASmehr
C
C---------------------------------------------------------------------
C       If KSYMOP (the operator symmetry) equals ISYMTR (the reference
C       state symmetry times excited state symmmety) calculate the
C       transition moment.
C---------------------------------------------------------------------
C
         IF (KSYMOP .EQ. ISYMTR) THEN
C
C-------------------------------------
C           Allocation of work space.
C-------------------------------------
C
            LTRMOM  = NEXCI
C
            KTRMOM  = 1 
C
C----------------------------------
C     Space for MO property matrix.
C----------------------------------
C
            LPRP1 = N2BST(ISYMTR)
            KPRP1 = KTRMOM + LTRMOM
            KEND1 = KPRP1 + LPRP1
            LWORK1  = LWORK  - KEND1
C

C     And perturbed density matrix            
            LPDENSIJ = NIJDEN(ISYMTR)
            IF (MODEL.EQ.'AORPA') THEN
               LPDENSAB = 0
            ELSE
               LPDENSAB = NABDEN(ISYMTR)
            ENDIF
            LPDENSAI = NAIDEN(ISYMTR)
         
            KPDENSIJ = KEND1 
            KPDENSAB = KPDENSIJ + LPDENSIJ
            KPDENSAI = KPDENSAB + LPDENSAB
            KEND2    = KPDENSAI + LPDENSAI
            LWORK2   = LWORK - KEND2

            CALL SO_MEMMAX ('SO_TRMOM',LWORK2)
            IF (LWORK2 .LT. 0)  CALL STOPIT('SO_TRMOM',' ',KEND2,LWORK)
C
C----------------------------------------------------------------------
C           Transform the property integrals to MO basis 
C----------------------------------------------------------------------
C
            CALL SO_ONEPMO(WORK(KPRP1),LPRP1,LABEL,ISYMTR,
     &                     RTNLBL,WORK(KEND1),LWORK1)
C
            IMAGPROP = RTNLBL(2).EQ.'ANTISYMM' 
C
C           Factor for D-ex part of gradient
            DFACTOR = -ONE
            IF (IMAGPROP) DFACTOR = ONE
C
C-----------------------------------------------------------
C           Loop over excitations of the specified symmetry.
C-----------------------------------------------------------
C
            DO 100 IEXCI = 1,NEXCI
C
C-------------------------------------------------------------
C              Calculate excitation part of transition moment.
C-------------------------------------------------------------
C
               EEXCIT=EXENG(ISYMTR,IEXCI)
               CALL SO_FNDPOS(PDENS_LABEL,EEXCIT,LURDENSE)
               LPDENSTOT = LPDENSIJ + LPDENSAB + LPDENSAI
               IF (MODEL.EQ.'AORPA') THEN
                  CALL SO_REAVE(WORK(KPDENSAI),LPDENSAI,ISYMTR,
     &                          PDENS_LABEL,EEXCIT,LURDENSE)
                  CALL DZERO(WORK(KPDENSIJ),LPDENSIJ)
               ELSE
                  CALL SO_REAVE(WORK(KPDENSIJ),LPDENSTOT,ISYMTR,
     &                          PDENS_LABEL,EEXCIT,LURDENSE) 
               ENDIF
               CALL SO_PROPMO(ISYMTR,CONTE,
     &                        MODEL.NE.'AORPA',IMAGPROP,
     &                        WORK(KPRP1),LPRP1,
     &                        WORK(KPDENSIJ),LPDENSIJ,
     &                        WORK(KPDENSAB),LPDENSAB,
     &                        WORK(KPDENSAI),LPDENSAI)
C
C----------------------------------------------------------------
C              Calculate de-excitation part of transition moment.
C----------------------------------------------------------------
C
               CALL SO_FNDPOS(PDENS_LABEL,EEXCIT,LURDENSD)
               LPDENSTOT = LPDENSIJ + LPDENSAB + LPDENSAI
               IF (MODEL.EQ.'AORPA') THEN
                  CALL SO_REAVE(WORK(KPDENSAI),LPDENSAI,ISYMTR,
     &                          PDENS_LABEL,EEXCIT,LURDENSD)
                  CALL DZERO(WORK(KPDENSIJ),LPDENSIJ)
               ELSE
                  CALL SO_REAVE(WORK(KPDENSIJ),LPDENSTOT,ISYMTR,
     &                          PDENS_LABEL,EEXCIT,LURDENSD) 
               ENDIF
               CALL SO_PROPMO(ISYMTR,CONTD,
     &                        MODEL.NE.'AORPA',IMAGPROP,
     &                        WORK(KPRP1),LPRP1,
     &                        WORK(KPDENSIJ),LPDENSIJ,
     &                        WORK(KPDENSAB),LPDENSAB,
     &                        WORK(KPDENSAI),LPDENSAI)
C
C-------------------------------------------------------------------
C              Calculate total transiton moment and write to output.
C-------------------------------------------------------------------
C
               WORK(KTRMOM+IEXCI-1) = CONTE + DFACTOR * CONTD
               IF ( IPRSOP .GE. 4) THEN
C
                  WRITE(LUPRI,'(/,1X,A,A,I3,A)')
     &            LABEL,' transition moment for state',IEXCI,':'
                  WRITE(LUPRI,'(2(F10.5,A),F10.5)')
     &            CONTE,' + ',CONTD*DFACTOR,
     &            ' = ',WORK(KTRMOM+IEXCI-1)
C
               END IF
C
C---------------------------------------------------
C              Collect transition moments in arrays.
C---------------------------------------------------
C
               CALL SO_COLLECT_TM(ISYMTR,LABEL,IEXCI,
     &                            WORK(KTRMOM+IEXCI-1),TRLEN,TRVEL,
CSPAS:23/5-11: second and third moment sum rules
     &                            TQLEN,TQVEL,TTLEN,
CKeinSPASmehr
     &                            TRLON,TRMAG,BSRLON,EXENG)
C
  100       CONTINUE
C
         END IF
C
  200 CONTINUE

      CALL GPCLOSE(LURDENSE,'KEEP')
      CALL GPCLOSE(LURDENSD,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TRMOM')
C
      RETURN
C
 9001 FORMAT(1X,'================================================',
     &       '===================')
 9002 FORMAT(1X,'------------------------------------------------',
     &       '-------------------')
C
      END
