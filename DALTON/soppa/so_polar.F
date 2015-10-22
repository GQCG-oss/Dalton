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
C
C     PURPOSE: Calculates the frequency dependent linear response 
C              properties from the solution and residual vectors and 
C              appropriate property gradient vectors with the atomic 
C              integral direct SOPPA program.
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
      DIMENSION   SNDPRP(2)
C
      LOGICAL   HASDOUBLES
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
      IF (.NOT.HASDOUBLES) THEN
C
         LGPVC1  = NT1AM(ISYMTR) 
         LGPVEC  = 2 * LGPVC1
         LSLV1E  = NT1AM(ISYMTR) 
         LSLV1D  = NT1AM(ISYMTR) 
C
         KGPVEC  = 1
         KGPVC1  = KGPVEC
         KGPVC1D = KGPVC1 + LGPVC1
         KSLV1   = KGPVEC + LGPVEC 
         KSLV1E  = KSLV1
         KSLV1D  = KSLV1E + LSLV1E
         KEND1   = KSLV1D + LSLV1D 
C
      ELSE ! IF (AOSOP .OR. AOSOC) THEN
C
         LGPVC1 = 2*NT1AM(ISYMTR) 
CPFP   2009-Jun
         LGPVC2 = 2*NT2AM(ISYMTR) 
C         LGPVC2 = N2P2HOP(ISYMTR) 
Cend-PFP
         LGPVEC = 2 * (LGPVC1 + LGPVC2)
         LSLV1  = NT1AM(ISYMTR)
CPFP   2009-Jun
         LSLV2  = NT2AM(ISYMTR)
C         LSLV2  = N2P2HOP(ISYMTR)
Cend-PFP
         KGPVC1  = 1
         KGPVC2  = KGPVC1 + LGPVC1
         KSLV1   = KGPVC2 + LGPVC2
         KSLV1E  = KSLV1
         KSLV1D  = KSLV1  + LSLV1
         KSLV2   = KSLV1D + LSLV1 
         KSLV2E  = KSLV2 
         KSLV2D  = KSLV2  + LSLV2 
         KEND1   = KSLV2D + LSLV2 
C     
      END IF
C
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('SO_POLAR.1',LWORK1)
      IF (LWORK1 .LT.0) CALL STOPIT('SO_POLAR.1',' ',KEND1,LWORK)
C
C--------------------------------------
C     Open files with solution vectors.
C--------------------------------------
C
      LUSV1E = -1 
      CALL GPOPEN(LUSV1E,FNSV1E,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
      LUSV1D = -1
      CALL GPOPEN(LUSV1D,FNSV1D,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
C
      IF ( HASDOUBLES ) THEN
C
       LUSV2E = -1
       CALL GPOPEN(LUSV2E,FNSV2E,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
       LUSV2D = -1
       CALL GPOPEN(LUSV2D,FNSV2D,'OLD',' ','UNFORMATTED',IDUMMY,.FALSE.)
C
      ENDIF

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
            IF ( MODEL .EQ. 'AORPA' ) THEN
C
               CALL RP_GETGP(WORK(KGPVEC),LGPVEC,LABEL2,ISYMTR,
     &                       WORK(KEND1),LWORK1)
C
            ELSE IF ( (MODEL .EQ. 'AOSOP') .OR.
     &                (MODEL .EQ. 'AOSOC')       ) THEN
C
cLig <> since I have changed the SO_GETGP function you have to call with
cLig    LGPVC1*2 and LGPVC2*2
               CALL SO_GETGP(WORK(KGPVC1),LGPVC1,WORK(KGPVC2),LGPVC2,
     &                       LABEL2,ISYMTR,T2AM,LT2AM,DENSIJ,LDENSIJ,
     &                       DENSAB,LDENSAB,DENSAI,LDENSAI,
     &                       WORK(KEND1),LWORK1)
C
            END IF
C
            IF ( IPRSOP .GE. 5) THEN
               IF ( MODEL .EQ. 'AORPA' ) THEN
                  CALL AROUND('In SO_POLAR:  RPA gradient property '//
     &                         'vector.'//LABEL2)
                  CALL OUTPUT(WORK(KGPVEC),1,LGPVC1,1,2,LGPVC1,2,1,
     &                        LUPRI)
               ELSE IF ( MODEL .EQ. 'AOSOP' ) THEN
                  CALL AROUND('In SO_POLAR:  SOPPA gradient property '//
     &                         'vector.'//LABEL2)
                  CALL OUTPUT(WORK(KGPVC1),1,LGPVC1/2,1,2,LGPVC1/2,2,1,
     &                        LUPRI)
                  CALL OUTPUT(WORK(KGPVC2),1,LGPVC2/2,1,2,LGPVC2/2,2,1,
     &                        LUPRI)
               ELSE IF ( MODEL .EQ. 'AOSOC' ) THEN
                  CALL AROUND('In SO_POLAR:  SOPPA(CCSD) gradient '//
     &                         'property vector.'//LABEL2)
                  CALL OUTPUT(WORK(KGPVC1),1,LGPVC1,1,2,LGPVC1,2,1,
     &                        LUPRI)
                  CALL OUTPUT(WORK(KGPVC2),1,LGPVC2,1,2,LGPVC2,2,1,
     &                        LUPRI)
               END IF
            END IF
C
C===============================================
C           Form second order properties SNDPRP.
C===============================================
C     
            IF (.NOT. ABSORP) THEN
C     
               DO 100 IFRVAL = 1, NFRVAL
C
C----------------------------------------------------------
C           Find and read the appropriate solution vectors.
C----------------------------------------------------------
C
                  CALL SO_FNDPOS(LABEL1,FRVAL(IFRVAL),LUSV1E)
                  CALL SO_FNDPOS(LABEL1,FRVAL(IFRVAL),LUSV1D)
C
                  IF ( HASDOUBLES ) THEN
C
                    CALL SO_REAVE(WORK(KSLV1E),LSLV1,ISYMTR,LABEL1,
     &                            FRVAL(IFRVAL),LUSV1E)
                    CALL SO_REAVE(WORK(KSLV1D),LSLV1,ISYMTR,LABEL1,
     &                            FRVAL(IFRVAL),LUSV1D)
C
                    CALL SO_FNDPOS(LABEL1,FRVAL(IFRVAL),LUSV2E)
                    CALL SO_FNDPOS(LABEL1,FRVAL(IFRVAL),LUSV2D)
C
                    CALL SO_REAVE(WORK(KSLV2E),LSLV2,ISYMTR,LABEL1,
     &                            FRVAL(IFRVAL),LUSV2E)
                    CALL SO_REAVE(WORK(KSLV2D),LSLV2,ISYMTR,LABEL1,
     &                            FRVAL(IFRVAL),LUSV2D)
C
                  ELSE
C
                    CALL SO_REAVE(WORK(KSLV1E),LSLV1E,ISYMTR,LABEL1,
     &                            FRVAL(IFRVAL),LUSV1E)
                    CALL SO_REAVE(WORK(KSLV1D),LSLV1D,ISYMTR,LABEL1,
     &                            FRVAL(IFRVAL),LUSV1D)
C
                  ENDIF
C
                  IF (IPRSOP .GT. 3) THEN
                     WRITE(LUPRI,'(A,I4,3A,F15.8)') 
     &                   'Solution Vector no. ',IFRVAL,
     &                   ' for operator ',LABEL1,
     &                   ' and frequency ',FREQ
                     CALL OUTPUT(WORK(KSLV1),1,LSLV1E,1,2,LSLV1E,
     &                           2,1,LUPRI)
                  ENDIF
C
C---------------------------------------------------------------------
C           Calculate second order properties SNDPRP.
C---------------------------------------------------------------------
C
                  IF ( MODEL .EQ. 'AORPA' ) THEN
C
                     SNDPRP(1)=DDOT(LGPVEC,WORK(KGPVEC),1,WORK(KSLV1),1)
C
                  ELSE IF ( (MODEL .EQ. 'AOSOP') .OR. 
     &                      (MODEL .EQ. 'AOSOC')      ) THEN
C
                     SNDPRP(1)=DDOT(LGPVC1,WORK(KGPVC1),1,WORK(KSLV1),1)
C
                     SNDPRP(1)=SNDPRP(1)
     &                        +DDOT(LGPVC2,WORK(KGPVC2),1,WORK(KSLV2),1)
C
                  ENDIF
C
                  IF (IPRSOP .GT. 2) THEN
                     WRITE (LUPRI,'(/,A,F15.8)')
     &                   ' Frequency = ',FRVAL(IFRVAL)
                     WRITE (LUPRI,'(4A,F15.8)')
     &                   ' Second order property for ',
     &                   LABEL2,LABEL1,' = ',SNDPRP(1)
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
      CALL GPCLOSE(LUSV1E,'KEEP')
      CALL GPCLOSE(LUSV1D,'KEEP')
C
      IF ( HASDOUBLES ) THEN
C
         CALL GPCLOSE(LUSV2E,'KEEP')
         CALL GPCLOSE(LUSV2D,'KEEP')
C
      ENDIF
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
