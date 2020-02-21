C
C  /* Deck so_trmom */
      SUBROUTINE SO_TRMOM_1(MODEL,LABEL,KSYMOP,ISYMTR,NEXCI,
     &                    T2AM,LT2AM,DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                    TRLEN,TRVEL,TQLEN,TQVEL,TTLEN,TRLON,TRMAG,
     &                    BSRLON,EXENG,TRGOS,WORK,LWORK,
     &                    LURDENSE,LURDENSD)
CClark:7/1/2016
C     This routine is to read integrals and process them
C     for certain label in SOPPA program.
C
C     PURPOSE: Calculate RPA, RPA(D), SOPPA and SOPPA(CCSD) transition
C              moments.
CClark:endC
      use so_info, only: FN_RDENSE, FN_RDENSD, so_full_name
C
#include "implicit.h"
#include "priunit.h"
#include "cbiexc.h"
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
      LOGICAL  IMAGPROP
      DOUBLE PRECISION :: DLABEL
C
CClark:7/1/2016
      DIMENSION   TRGOS(3,MXNEXI)
CClark:end
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TRMOM_1')
      WRITE(PDENS_LABEL,'(A7,I1)') 'EXCITA ', ISYMTR
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

            IF ( IPRSOP .GE. 1) THEN

                WRITE(LUPRI,'(/,1X,4A,I3)')
     &          LABEL,' ',TRIM(SO_FULL_NAME(MODEL)),
     &          ' transition moments for excitation symmetry',
     &          ISYMTR
                WRITE(LUPRI,9001)
                WRITE(LUPRI,'(A)')
     &          ' Exci.    1p-1h        1h-1p        Total'
                WRITE(LUPRI,9002)
            ENDIF
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
CRF          Note - Label is now number within irrep!
               DLABEL = DBLE(IEXCI)
               LPDENSTOT = LPDENSIJ + LPDENSAB + LPDENSAI
               IF (MODEL.EQ.'AORPA') THEN
                  CALL SO_REAVE(WORK(KPDENSAI),LPDENSAI,ISYMTR,
     &                          PDENS_LABEL,DLABEL,LURDENSE)
                  CALL DZERO(WORK(KPDENSIJ),LPDENSIJ)
               ELSE
                  CALL SO_REAVE(WORK(KPDENSIJ),LPDENSTOT,ISYMTR,
     &                          PDENS_LABEL,DLABEL,LURDENSE)
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
               LPDENSTOT = LPDENSIJ + LPDENSAB + LPDENSAI
               IF (MODEL.EQ.'AORPA') THEN
                  CALL SO_REAVE(WORK(KPDENSAI),LPDENSAI,ISYMTR,
     &                          PDENS_LABEL,DLABEL,LURDENSD)
                  CALL DZERO(WORK(KPDENSIJ),LPDENSIJ)
               ELSE
                  CALL SO_REAVE(WORK(KPDENSIJ),LPDENSTOT,ISYMTR,
     &                          PDENS_LABEL,DLABEL,LURDENSD)
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
               IF ( IPRSOP .GE. 1) THEN
C
                  WRITE(LUPRI,'(1X,I3,2X,2(F10.5,A),F10.5)') IEXCI,
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
CClark:7/1/2016 add TRGOS to the argument list
     &                            TRLON,TRMAG,TRGOS,BSRLON,EXENG)
CClark:end
C

  100       CONTINUE
C
C---------------------------------------------------------
C              Write tabel of transition moments to output.
C---------------------------------------------------------
C
C
             IF ( IPRSOP .GE. 1) WRITE(LUPRI,9001)
C
         END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_TRMOM_1')
C
      RETURN
C
 9001 FORMAT(1X,'================================================',
     &       '===================')
 9002 FORMAT(1X,'------------------------------------------------',
     &       '-------------------')
C
      END
