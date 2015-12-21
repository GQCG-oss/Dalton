C
C  /* Deck so_polar */
      SUBROUTINE SO_POLAR(MODEL,ISYMTR,IDIP,LABEL1,T2AM,LT2AM,
     &                    DENSIJ,LDENSIJ,
     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,POLDD,POLDQ,
     &                    POLDL,POLDA,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Stephan P. A. Sauer: 5.12.2003
C     Rasmus Faber: Nov. 2015 --> Implemented pert. dens. approach.
C         
C     PURPOSE: Calculates the frequency dependent linear response 
C              properties from the perturbed density matrices and 
C              appropriate property integrals with the atomic 
C              integral direct SOPPA program.
C
      use so_info, only: fn_rdens, sop_stat_trh
C
#include "implicit.h"
#include "priunit.h"
C
#include "soppinf.h"
#include "cbilnr.h"
#include "cbiexc.h"
#include "absorp.h"
#include "ccsdsym.h"
c#include "infdim.h"
c#include "inforb.h"
c#include "iratdef.h"
c#include "maxaqn.h"
c#include "maxorb.h"
c#include "mxcent.h"
c#include "nuclei.h"
c#include "symmet.h"
c#include "wrkrsp.h"
C
C
      DIMENSION POLDD(2,3,3,NFRVAL), POLDQ(2,3,3,3,NFRVAL)
      DIMENSION POLDL(2,3,3,NFRVAL), POLDA(2,3,3,NFRVAL)
      DIMENSION DENSIJ(LDENSIJ), DENSAB(LDENSAB), DENSAI(LDENSAI)
      DIMENSION T2AM(LT2AM)
      DIMENSION WORK(LWORK)
C
      CHARACTER*8 LABEL1, LABEL2
      CHARACTER*5 MODEL
      CHARACTER*8 RTNLBL(2)
      DIMENSION   SNDPRP(2)
C
      LOGICAL   HASDOUBLES, IMAGPROP
      LOGICAL   THIS_STATIC
      PARAMETER (DP5=0.5D0)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_POLAR')
C
C------------------------------
C     Allocation of work space.
C------------------------------
C
      HASDOUBLES = MODEL .NE. 'AORPA'
      KEND1 = 1
C      IF (.NOT.HASDOUBLES) THEN
CC
C         LGPVC1  = NT1AM(ISYMTR) 
CCRF         LGPVEC  = 2 * LGPVC1
C         LGPVEC  = LGPVC1
C         LSLV1   = NT1AM(ISYMTR) 
CCRF         LSLV1D  = NT1AM(ISYMTR) 
CC
C         KGPVEC  = 1
C         KGPVC1  = KGPVEC
CC         KGPVC1D = KGPVC1 + LGPVC1
C         KSLV1   = KGPVEC + LGPVEC 
C         KSLV1E  = KSLV1
C         KSLV1D  = KSLV1E + LSLV1
C         KEND1   = KSLV1D + LSLV1
CC 
C         LGPVC2  = 1
C         KGPVC2  = KEND1         
CC
C      ELSE ! IF (AOSOP .OR. AOSOC) THEN
CC
CCRF         LGPVC1 = 2*NT1AM(ISYMTR) 
C         LGPVC1 = NT1AM(ISYMTR) 
CCPFP   2009-Jun
CCRF         LGPVC2 = 2*NT2AM(ISYMTR) 
C         LGPVC2 = NT2AM(ISYMTR) 
CC         LGPVC2 = N2P2HOP(ISYMTR) 
CCend-PFP
CCRf         LGPVEC = 2 * (LGPVC1 + LGPVC2)
C         LGPVEC = (LGPVC1 + LGPVC2)
C         LSLV1  = NT1AM(ISYMTR)
CCPFP   2009-Jun
C         LSLV2  = NT2AM(ISYMTR)
CC         LSLV2  = N2P2HOP(ISYMTR)
CCend-PFP
C         KGPVC1  = 1
C         KGPVC2  = KGPVC1 + LGPVC1
C         KSLV1   = KGPVC2 + LGPVC2
C         KSLV1E  = KSLV1
C         KSLV1D  = KSLV1  + LSLV1
C         KSLV2   = KSLV1D + LSLV1 
C         KSLV2E  = KSLV2 
C         KSLV2D  = KSLV2  + LSLV2 
C         KEND1   = KSLV2D + LSLV2 
CC     
C      END IF
C
C----------------------------------
C     Space for MO property matrix.
C----------------------------------
C
      LPRP1 = N2BST(ISYMTR)
      KPRP1 = KEND1
      KEND1 = KPRP1 + LPRP1
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_POLAR.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('SO_POLAR.1',' ',KEND1,LWORK)
C
C--------------------------------------
C     Open files with solution vectors.
C--------------------------------------
C
C      LUSV1E = -1 
C      CALL GPOPEN(LUSV1E,FNSV1E,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
C      LUSV1D = -1
C      CALL GPOPEN(LUSV1D,FNSV1D,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
CC
C      IF ( HASDOUBLES ) THEN
CC
C       LUSV2E = -1
C       CALL GPOPEN(LUSV2E,FNSV2E,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
C       LUSV2D = -1
C       CALL GPOPEN(LUSV2D,FNSV2D,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
CC
C      ENDIF
C 
C     Open file with perturbed densities 
C  
      LURDENS = -1
      CALL GPOPEN(LURDENS, FN_RDENS, 'OLD',' ','UNFORMATTED',
     &            IDUMMY,.FALSE.)
C     
C=============================================
C     Loop over the second property operators.
C=============================================
C     
cspas      DO 200 IPRLBL = 1, NLBTOT
      DO 200 IPRLBL = 1, 3
C
C---------------------------------------------------     
C        Find label and symmetry of second operator.
C---------------------------------------------------     
C     
         LABEL2 = LABAPP(IPRLBL)
         KSYM   = LABSYM(IPRLBL)
C
C--------------------------------------------------------     
C        If symmetry of first operator equals symmetry of
C        second operator, that is if ISYMTR = KSYM, then.
C--------------------------------------------------------
C     
         IF (KSYM .EQ. ISYMTR) THEN
C
C----------------------------------------------------------------------
C           Calculate the gradient property vectors and write to output
C----------------------------------------------------------------------
C
C            IF ( MODEL .EQ. 'AORPA' ) THEN
C
C               CALL RP_GETGP(WORK(KGPVEC),LGPVEC,LABEL2,ISYMTR,
C     &                       WORK(KEND1),LWORK1)
C
C            ELSE IF ( (MODEL .EQ. 'AOSOP') .OR.
C     &                (MODEL .EQ. 'AOSOC')       ) THEN
C
cLig <> since I have changed the SO_GETGP function you have to call with
cLig    LGPVC1*2 and LGPVC2*2
CRF  I changed that back !!!!
CRF  Also we should be able to call the same GETGP for RPA
C            CALL SO_GETGP(WORK(KGPVC1),LGPVC1,WORK(KGPVC2),LGPVC2,
C     &                    LABEL2,ISYMTR,IMAGPROP,MODEL,
C     &                    T2AM,LT2AM,DENSIJ,LDENSIJ,
C     &                    DENSAB,LDENSAB,DENSAI,LDENSAI,
C     &                    WORK(KEND1),LWORK1)
C
C--------------------------------------------------
C           Get the property integrals in MO basis.
C--------------------------------------------------
C            
            CALL SO_ONEPMO(WORK(KPRP1),LPRP1,LABEL2,ISYMTR,
     &                     RTNLBL,WORK(KEND1),LWORK1)
C
            IMAGPROP = RTNLBL(2).EQ.'ANTISYMM' 

            DFACTOR = -1.0D0
            IF (IMAGPROP) DFACTOR = 1.0D0

            

C
C            END IF
C
C            IF ( IPRSOP .GE. 5) THEN
C               IF ( MODEL .EQ. 'AORPA' ) THEN
C                  CALL AROUND('In SO_POLAR:  RPA gradient property '//
C     &                         'vector.'//LABEL2)
C                  CALL OUTPUT(WORK(KGPVEC),1,LGPVC1,1,1,LGPVC1,1,1,
C     &                        LUPRI)
C               ELSE IF ( MODEL .EQ. 'AOSOP' ) THEN
C                  CALL AROUND('In SO_POLAR:  SOPPA gradient property '//
C     &                         'vector.'//LABEL2)
C                  CALL OUTPUT(WORK(KGPVC1),1,LGPVC1,1,1,LGPVC1,1,1,
C     &                        LUPRI)
C                  CALL OUTPUT(WORK(KGPVC2),1,LGPVC2,1,1,LGPVC2,1,1,
C     &                        LUPRI)
C               ELSE IF ( MODEL .EQ. 'AOSOC' ) THEN
C                  CALL AROUND('In SO_POLAR:  SOPPA(CCSD) gradient '//
C     &                         'property vector.'//LABEL2)
C                  CALL OUTPUT(WORK(KGPVC1),1,LGPVC1,1,1,LGPVC1,1,1,
C     &                        LUPRI)
C                  CALL OUTPUT(WORK(KGPVC2),1,LGPVC2,1,1,LGPVC2,1,1,
C     &                        LUPRI)
C               END IF
C            END IF
C
C===============================================
C           Form second order properties SNDPRP.
C===============================================
C     
            IF (.NOT. ABSORP) THEN
C     
               DO 100 IFRVAL = 1, NFRVAL
C                   THIS_STATIC = ABS(FRVAL(IFRVAL)).LT.SOP_STAT_TRH
                  THIS_STATIC = .FALSE.
C
C----------------------------------------------------------
C           Find and read the appropriate solution vectors.
C----------------------------------------------------------
C
C                  CALL SO_FNDPOS(LABEL1,FRVAL(IFRVAL),LUSV1E)
C                  CALL SO_FNDPOS(LABEL1,FRVAL(IFRVAL),LUSV1D)
C                  CALL SO_REAVE(WORK(KSLV1E),LSLV1,ISYMTR,LABEL1,
C     &                          FRVAL(IFRVAL),LUSV1E)
C                  CALL SO_REAVE(WORK(KSLV1D),LSLV1,ISYMTR,LABEL1,
C     &                          FRVAL(IFRVAL),LUSV1D)
CC
C                  IF ( HASDOUBLES ) THEN
CC
CC
C                    CALL SO_FNDPOS(LABEL1,FRVAL(IFRVAL),LUSV2E)
C                    CALL SO_FNDPOS(LABEL1,FRVAL(IFRVAL),LUSV2D)
CC
C                    CALL SO_REAVE(WORK(KSLV2E),LSLV2,ISYMTR,LABEL1,
C     &                            FRVAL(IFRVAL),LUSV2E)
C                    CALL SO_REAVE(WORK(KSLV2D),LSLV2,ISYMTR,LABEL1,
C     &                            FRVAL(IFRVAL),LUSV2D)
C                  ENDIF
CC
C                  IF (IPRSOP .GT. 2) THEN
C                     WRITE(LUPRI,'(A,I4,3A,F15.8)') 
C     &                   'Solution Vector no. ',IFRVAL,
C     &                   ' for operator ',LABEL1,
C     &                   ' and frequency ',FRVAL(IFRVAL)
C                     CALL OUTPUT(WORK(KSLV1),1,LSLV1,1,2,LSLV1,
C     &                           2,1,LUPRI)
C                  ENDIF
C
C---------------------------------------------------------------------
C           Calculate second order properties SNDPRP.
C---------------------------------------------------------------------
C
C           Singles contribution                  
C                  PRPOLD1 = DDOT(LGPVC1,WORK(KGPVC1),1,WORK(KSLV1E),1)
C     &                      +DFACTOR *
C     &                       DDOT(LGPVC1,WORK(KGPVC1),1,WORK(KSLV1D),1)
C
C                  IF ( MODEL .EQ. 'AORPA' ) THEN
C
C                     SNDPRP(1)=DDOT(LGPVEC,WORK(KGPVEC),1,WORK(KSLV1),1)
C
C                  IF ( HASDOUBLES ) THEN
C           Doubles contribution
C                     PRPOLD2 = 
C     &                      + DDOT(LGPVC2,WORK(KGPVC2),1,WORK(KSLV2E),1)
C     &                      + DFACTOR *
C     &                        DDOT(LGPVC2,WORK(KGPVC2),1,WORK(KSLV2D),1)
C
C                  ENDIF
C
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
                  IF (LWORK2 .LT. 0) WRITE(LUPRI,*) 'Insufficient '//
     &                       'memory, new so_polar'
C
C
C        Test --> calculate perturbed density instead
C
C                  IF(.FALSE.)THEN
C                  IF (THIS_STATIC) THEN
C                     FACTOR = 2.0D0
C                  ELSE 
C                     FACTOR = 1.0D0
C                  END IF
CC
CC           Initialize perturbed densities
C                  CALL DZERO(WORK(KPDENSIJ),LPDENSIJ+LPDENSAB+LPDENSAI)
C
C                  CALL SO_PERTD1(MODEL,ISYMTR,FACTOR,
C     &                           WORK(KPDENSIJ),LPDENSIJ,
C     &                           WORK(KPDENSAB),LPDENSAB,
C     &                           WORK(KPDENSAI),LPDENSAI,
C     &                           WORK(KSLV1E),LSLV1,
C     &                           DENSIJ,LDENSIJ,
C     &                           DENSAB,LDENSAB,DENSAI,LDENSAI,
C     %                           WORK(KEND2),LWORK2)
C                  IF (HASDOUBLES) THEN
C                     CALL SO_PERTD2(ISYMTR,FACTOR,
C     &                              WORK(KPDENSIJ),LPDENSIJ,
C     &                              WORK(KPDENSAB),LPDENSAB,
C     &                              T2AM,LT2AM,
C     &                              WORK(KSLV2E),LSLV2)
C                  ENDIF
C
C                  IF (.NOT.THIS_STATIC) THEN !NOT STATIC
C               
C                     IF (IMAGPROP) THEN  
C                        FACTOR =  1.0D0
C                     ELSE 
C                        FACTOR = -1.0D0
C                     END IF
C
C                     CALL SO_PERTD1(MODEL,ISYMTR,FACTOR,
C     &                              WORK(KPDENSIJ),LPDENSIJ,
C     &                              WORK(KPDENSAB),LPDENSAB,
C     &                              WORK(KPDENSAI),LPDENSAI,
C     &                              WORK(KSLV1D),LSLV1,
C     &                              DENSIJ,LDENSIJ,
C     &                              DENSAB,LDENSAB,DENSAI,LDENSAI,
C     %                              WORK(KEND2),LWORK2)
C                     IF (HASDOUBLES) THEN
C                        CALL SO_PERTD2(ISYMTR,FACTOR,
C     &                                 WORK(KPDENSIJ),LPDENSIJ,
C     &                                 WORK(KPDENSAB),LPDENSAB,
C     &                                 T2AM,LT2AM,
C     &                                 WORK(KSLV2D),LSLV2)
C                     ENDIF
C     
C
C                  END IF
C                 
C-----------------------------------------------------------
C           Get the perturbed density matrix from file
C-----------------------------------------------------------
C                  
                  CALL SO_FNDPOS(LABEL1,FRVAL(IFRVAL),LURDENS)
                  LPDENSTOT = LPDENSIJ + LPDENSAB + LPDENSAI
                  IF (MODEL.EQ.'AORPA') THEN
                     CALL SO_REAVE(WORK(KPDENSAI),LPDENSAI,ISYMTR,
     &                             LABEL1,FRVAL(IFRVAL),LURDENS)
                     CALL DZERO(WORK(KPDENSIJ),LPDENSIJ)
                  ELSE
                     CALL SO_REAVE(WORK(KPDENSIJ),LPDENSTOT,ISYMTR,
     &                             LABEL1,FRVAL(IFRVAL),LURDENS) 
                  ENDIF

C                  WRITE(LUPRI,'(a)') 'Perturbed density of '//LABEL1//
C     &                               'read'
C                  WRITE(LUPRI,'(5f14.8)') (WORK(KPDENSIJ-1+I),
C     &                           I=1,LPDENSIJ)
C                  WRITE(LUPRI,'(5f14.8)') (WORK(KPDENSAB-1+I),
C     &                           I=1,LPDENSAB)
C                  WRITE(LUPRI,'(5f14.8)') (WORK(KPDENSAI-1+I),
C     &                           I=1,LPDENSAI)
C
C---------------------------------------------------------------------
C           Calculate second order properties SNDPRP.
C---------------------------------------------------------------------
C
                  CALL SO_PROPMO(ISYMTR,SNDPRP(1),
     &                           MODEL.NE.'AORPA',IMAGPROP,
     &                           WORK(KPRP1),LPRP1,
     &                           WORK(KPDENSIJ),LPDENSIJ,
     &                           WORK(KPDENSAB),LPDENSAB,
     &                           WORK(KPDENSAI),LPDENSAI)
C
                  IF (IPRSOP .GT. 4) THEN
                     WRITE (LUPRI,'(/,A,F15.8)')
     &                   ' Frequency = ',FRVAL(IFRVAL)
                     WRITE (LUPRI,'(4A,F15.8)')
     &                   ' Second order property for ',
     &                   LABEL2,LABEL1,' = ',SNDPRP(1)
C                     WRITE (LUPRI,'(A,F15.8)') 'Old singles: ',
C     &                         PRPOLD1
C                     IF (HASDOUBLES) 
C     &                      WRITE (LUPRI,'(A,F15.8)') 'Old doubles: ',
C     &                                PRPOLD2
C                     WRITE (LUPRI,'(A,F15.8)') 'New doubles: ',
C     &                         SNDPRP(1) - PRPOLD1

                  ENDIF
C     
C---------------------------------------------------------------------
C                 Write properties into the various property matrices.
C---------------------------------------------------------------------
C     
                  IF (LABEL2(2:7).EQ.'DIPLEN') THEN
C
                     IF (LABEL2(1:1).EQ.'X') THEN
                        POLDD(1,IDIP,1,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:1).EQ.'Y') THEN
                        POLDD(1,IDIP,2,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:1).EQ.'Z') THEN
                        POLDD(1,IDIP,3,IFRVAL) = SNDPRP(1)
                     END IF
C
                  ELSE IF (LABEL2(3:8).EQ.'THETA ') THEN
C     
                     IF (LABEL2(1:2).EQ.'XX') THEN
                        POLDQ(1,IDIP,1,1,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:2).EQ.'XY') THEN
                        POLDQ(1,IDIP,1,2,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:2).EQ.'XY') THEN
                        POLDQ(1,IDIP,2,1,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:2).EQ.'XZ') THEN
                        POLDQ(1,IDIP,1,3,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:2).EQ.'XZ') THEN
                        POLDQ(1,IDIP,3,1,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:2).EQ.'YY') THEN
                        POLDQ(1,IDIP,2,2,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:2).EQ.'YZ') THEN
                        POLDQ(1,IDIP,2,3,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:2).EQ.'YZ') THEN
                        POLDQ(1,IDIP,3,2,IFRVAL) = SNDPRP(1)
                     END IF
                     IF (LABEL2(1:2).EQ.'ZZ') THEN
                        POLDQ(1,IDIP,3,3,IFRVAL) = SNDPRP(1)
                     END IF
C     
                  ELSE IF (LABEL2(2:7).EQ.'LONMAG') THEN
C     
C----------------------------------------------------------------
C                 Multiply with minus the Bohr-magneton (-0.5) to
C                 create the magnetic dipole operator from the 
C                 angular momentum operator.
C----------------------------------------------------------------
C     
                     IF (LABEL2(1:1).EQ.'X') THEN
                        POLDL(1,IDIP,1,IFRVAL) = -DP5*SNDPRP(1)
                     END IF
                     IF (LABEL2(1:1).EQ.'Y') THEN
                        POLDL(1,IDIP,2,IFRVAL) = -DP5*SNDPRP(1)
                     END IF
                     IF (LABEL2(1:1).EQ.'Z') THEN
                        POLDL(1,IDIP,3,IFRVAL) = -DP5*SNDPRP(1)
                     END IF
C     
                  ELSE IF (LABEL2(2:7).EQ.'ANGMOM') THEN
C     
                     IF (LABEL2(1:1).EQ.'X') THEN
                        POLDA(1,IDIP,1,IFRVAL) = -DP5*SNDPRP(1)
                     END IF
                     IF (LABEL2(1:1).EQ.'Y') THEN
                        POLDA(1,IDIP,2,IFRVAL) = -DP5*SNDPRP(1)
                     END IF
                     IF (LABEL2(1:1).EQ.'Z') THEN
                        POLDA(1,IDIP,3,IFRVAL) = -DP5*SNDPRP(1)
                     END IF
C     
                  END IF
C
  100          CONTINUE
C
            END IF
C
         END IF
C
  200 CONTINUE
C
C--------------------------------------------
C           Close files with solution vectors.
C--------------------------------------------
C
C      CALL GPCLOSE(LUSV1E,'KEEP')
C      CALL GPCLOSE(LUSV1D,'KEEP')
CC
C      IF ( HASDOUBLES ) THEN
CC
C         CALL GPCLOSE(LUSV2E,'KEEP')
C         CALL GPCLOSE(LUSV2D,'KEEP')
CC
C      ENDIF
C
      CALL GPCLOSE(LURDENS,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_POLAR')
C
      RETURN
      END
