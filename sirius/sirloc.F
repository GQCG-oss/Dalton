!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
 
 
C  /* Deck sirloc */
      SUBROUTINE SIRLOC(CMO,WORK,LWORK)
C     Written by Claire C.M. Samson (University of Karlsruhe, 28 April 2003).
 
#include <implicit.h>
      PARAMETER (THRES = 1.0D-10)
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h> 
#include <infvar.h>
#include <r12int.h>
      DIMENSION CMO(*), WORK(LWORK)
      LOGICAL ANTIS
 
      CALL TITLER('BOYS LOCALIZATION','*',125)
      CALL NUMLOC(NOC,NBA,NB1,NB2)
      CALL FROLOC(NCO)
C
      WRITE(LUPRI,*) 'Number of frozen orbitals    :',NCO
      WRITE(LUPRI,*) 'Number of localized orbitals :',NOC-NCO
      WRITE(LUPRI,*) 'Number of molecular orbitals :',NB1
      WRITE(LUPRI,*) 'Number of auxilliary orbitals:',NB2
      WRITE(LUPRI,*) 'Number of basis functions    :',NBA
C
C Parameters are set for valence orbitals only
 
      NOC = NOC - NCO
      NTR  = ((NOC-1)*NOC)/2
      NN   = NBA*NBA
      NOC3 = NOC*3
      NTR3 = NTR*3
      NCMO = NCMOT - NCO*NBA
      NKB   = NBA*NOC
      NKCL  = NOC*NOC
      NNBA1 = NBA*(NBA+1)/2
 
C     *****   memory allocation  *****
 
      KINI   = 1
      KCMO   = KINI
      KDIPX  = KCMO   + NCMOT + NBA*(NOC + NCO)
      KDIPY  = KDIPX  + NNBASX
      KDIPZ  = KDIPY  + NNBASX
      KILIFQ = KDIPZ  + NNBASX
      KIKY   = KILIFQ + NOC
      KMAXI0  = KIKY   + NBA
      KMINI0  = KMAXI0  + NN
      KRI    = KMINI0  + NN
      KRIJ   = KRI    + NN  
      KCL    = KRIJ   + NTR3
      KB     = KCL    + NKCL
      KQPIX  = KB     + NKB 
      KQPJX  = KQPIX  + NOC
      KPAO   = KQPJX  + NOC
      KEVC   = KPAO   + NB1*NBA
      KMOC   = KEVC   + NB1*NBA
      KOVER1 = KMOC   + (NOC+NCO)*NBA
      KOVER2 = KOVER1 + NNBAST
      KEND   = KOVER2 + NN
      LFREE  = LWORK  - KEND
 
      KCMOC  = KCMO  + NCO*NBA
 
C molecular coefficients and dipole moment integrals
 
      JRDMO = 9
      CALL READMO(WORK(KCMO),JRDMO,WORK(KEND),LFREE)
C     CALL AROUND('Canonical molecular orbitals')
C     CALL PRINTORB(NB1,NBA,WORK(KCMO))

      CALL RDPROP('XDIPLEN ',WORK(KDIPX),ANTIS)
      CALL RDPROP('YDIPLEN ',WORK(KDIPY),ANTIS)
      CALL RDPROP('ZDIPLEN ',WORK(KDIPZ),ANTIS)
 
      CALL BOYLOC(WORK(KCMOC),WORK(KDIPX),WORK(KILIFQ),WORK(KIKY),
     &            WORK(KMAXI0),WORK(KMINI0),WORK(KRI),WORK(KRIJ),
     &            WORK(KCL),WORK(KB),WORK(KQPIX),WORK(KQPJX),NCMO,
     &            NNBASX,NTR,NOC,NBA,WORK(KEND),LFREE)

Cccms CALL DCOPY(NCMOT,WORK(KCMO),1,CMO,1)
 
C Projected atomic orbitals
c 
c     NOC = NOC + NCO
c     CALL RDONEL('OVERLAP ',.TRUE.,WORK(KOVER1),NNBAST)
c     CALL PROAO(WORK(KPAO),WORK(KEVC),WORK(KCMO),WORK(KMOC),
c    *           WORK(KOVER1),WORK(KOVER2),NOC,NBA,NB1,NNBA1)
c
c     KSHI = KCMO + NBA*NB1
c     LSHI = NOC*NBA
c     CALL DZERO(WORK(KSHI),LSHI)
c
Cccms KSHI = KCMO + NBA*NOC
C     LSHI = NB1*NBA
C     CALL DZERO(WORK(KSHI),LSHI)
C     DO I = 1, NB1
C        KSHIP = KPAO + NB1*(I-1)
C        KSHIC = KCMO + NBA*(NOC+I-1)
C        CALL DCOPY(NB1,WORK(KSHIP),1,WORK(KSHIC),1)
C     END DO
C
C     CALL AROUND('Localized orbitals')
C     CALL PRINTORB(NB1,NBA,WORK(KCMO))
C
C *** IF R12EIN COMBINED WITH LOCAL MOS 
 
CWMK  IF (R12EOR) THEN
 
C *** Memory allocation

         KFCAO  = KEND
         KFVAO  = KFCAO  + N2BASX
         KFCMO  = KFVAO  + N2BASX
         KDCAO  = KFCMO  + N2BASX
         KDVAO  = KDCAO  + N2BASX  
         KFPAO  = KDVAO  + N2BASX
         KVEC   = KFPAO  + NNBAST
         KWRK   = KVEC   + NN
         LWRK   = LFREE  - KWRK 
 
C *** Tranform AO fock matrix in local MO fock matrix
 
         CALL FCKDEN(.TRUE.,.FALSE.,WORK(KDCAO),WORK(KDVAO),
     &                                 WORK(KCMO),DUMMY,WORK(KWRK),LWRK)
         CALL FCKMAO(.TRUE.,EMCMY,WORK(KFCAO),WORK(KFVAO),
     &               WORK(KDCAO),WORK(KDVAO),DV_DUMMY,WORK(KCMO),
     &               WORK(KWRK),LWRK)
         CALL DCOPY(N2BASX,WORK(KFCAO),1,WORK(KDCAO),1)
         CALL DGETSP(NBAST,WORK(KDCAO),WORK(KFCAO))
         CALL UTHU(WORK(KFCAO),WORK(KFCMO),WORK(KCMO),
     &             WORK(KDCAO),NBA,NB1)
C ccms diag pao
C         NTFMO = NOBL*(NOBL+1)/2
C         NB1T  = NB1 * (NB1+1)/2 
C      DO I = 1, NB1
C         KSHIP = KCMO + NBA*(NOC+I-1)
C         KSHIC = KVEC + NB1*(I-1)
C         CALL DCOPY(NB1,WORK(KSHIP),1,WORK(KSHIC),1)
C      END DO
C         CALL DIAPAO(WORK(KFCMO),WORK(KFPAO),WORK(KVEC),
C     *                                           NTFMO,NB1T,NB1,NOC)
C      DO I = 1, NB1
C         KSHIP = KCMO + NBA*(NOC+I-1)
C         KSHIC = KVEC + NB1*(I-1)
C         CALL DCOPY(NB1,WORK(KSHIC),1,WORK(KSHIP),1)
C      END DO
C      CALL DCOPY(NBA*NOC,CMO,1,WORK(KCMO),1)

      CALL WRLOCMO(WORK(KCMO),WORK(KFCMO),NOC,NCO,NBA,NB1)

C        WRITE(LUPRI,'(A)') 'SUBROUTINE FCKMAT'
C        CALL FCKMAT(.TRUE.,DUMMY,WORK(KCMO),EMY,WORK(KFCAO),DUMMY,
C    *               WORK(KWRK),LWRK)
C        CALL OUTPAK(WORK(KFCAO),NBAST,1,LUPRI)
C        CALL WRLOCMO(WORK(KCMO),WORK(KFCAO),NOC,NCO,NBA,NB1)
c
c     CALL AROUND('Molecular orbitals and Canonical PAOs')
c     CALL PRINTORB(NB1+NOC,NBA,WORK(KCMO))
c
CWMK  END IF
 
      CALL TITLER('END BOYS LOCALIZATION','*',103)
      RETURN
      END

C  /* Deck diapao  */
      SUBROUTINE DIAPAO(FCMO,FPAO,VEC,NTFMO,NB1T,NB1,NOC)
C     Written by Claire C.M. Samson (University of Karlsruhe, 28 April 2003).
#include <implicit.h>
#include <priunit.h>
      DIMENSION FCMO(NTFMO),FPAO(NB1T),VEC(NB1,NB1)

      NOBL = NOC + NB1
      IJ = 0
      KL = 0
      DO I = 1, NOBL
         DO J = 1, I
            IJ = IJ + 1
            IF (J.GT.NOC) THEN
               KL = KL +1
               FPAO(KL)=FCMO(IJ)
            END IF
         END DO
      END DO
      CALL JACOBI(FPAO,VEC,NB1,NB1)
      IJ = 0
      KL = 0
      DO I = 1, NOBL
         DO J = 1, I
            IJ = IJ + 1
            IF (J.GT.NOC) THEN
               KL = KL +1
               FCMO(IJ)=FPAO(KL)
            END IF
         END DO
      END DO

      END

C  /* Deck numloc  */
      SUBROUTINE NUMLOC(NOC,NBA,NB1,NB2)
C     Written by Claire C.M. Samson (University of Karlsruhe, 28 April 2003).
C     Revised by hjaaj Aug 2004: set NOC = NISHT instead of NOC = NOCCT

#include <implicit.h>
#include <inforb.h>
#include <r12int.h>

C INFORMATION ON CANONICAL OCCUPIED ORBITALS
C Symmetry is switched off
 
      NOC = 0
      NBA = 0
      NB1 = 0
      NB2 = 0
      DO ISYM = 1,NSYM
         NOC  = NOC + NISH(ISYM)
C        ... localize all doubly occupied orbitals
C            but not the active orbitals. /hjaaj aug 2004
         NBA  = NBA + NBAS(ISYM)
         NB1  = NB1 + NORB1(ISYM)
         NB2  = NB2 + NORB2(ISYM)
      END DO
 
      RETURN
      END

C  /* Deck froloc  */
      SUBROUTINE FROLOC(NCO)
C     Written by Claire C.M. Samson (University of Karlsruhe, 28 April 2003).
C
C     Revised by hjaaj Aug 2004 to use NFRO(ISYM) from inforb.h instead
C     of NRHFFR(ISYM) from ccorb.h which is not defined yet in Sirius!
C     NOTE that NRHFFR(isym) frozen orbitals from CC from ccorb.h not
C     necessarily are the first, thus also for that reason it cannot be used.
C     This will also work for MCSCF.
C
C     Revised by wmk Aug 2005 to use LOCFRO(ISYM).
C
#include <implicit.h>
#include <inforb.h>
#include <r12int.h>

C     LGLO = .TRUE.
C     hjaaj aug 2004: is in ccorb.h but is never used in cc/, thus
C     I just disabled the definition.
      NCO = 0
      DO ISYM = 1, NSYM
cwmk     NCO = NCO + NFRO(ISYM)
         NCO = NCO + LOCFRO(ISYM)
      END DO
         
      RETURN
      END

C  /* Deck boyloc */ 
      SUBROUTINE BOYLOC(CMO,DIPXYZ,ILIFQ,IKY,MAXI0,MINI0,RI,RIJ,CL,B,
     &                  QPIX,QPJX,NCMO,NDIP,NTR,NOC,NBA,WORK,LWORK)
C Purpose :
C     To localize molecular orbitals with Boys procedure
C     Adapted from the GAMESS-UK code
C
C     Written by Claire C.M. Samson (University of Utrecht, 30 November 2002).
 
#include <implicit.h>
       PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     &            D1M = -1.0D0, D2M = -2.0D0, DP4 = 0.25D0,
     &            THR10 = 1.0D-10, THR8 = 1.0D-8, THR6= 1.0D-6,
     &            THR8M = -1.0D-8)
       DIMENSION WORK(LWORK),CMO(NCMO),DIPXYZ(NDIP,3),ILIFQ(NOC),
     &           IKY(NBA),MAXI0(NBA,NBA),MINI0(NBA,NBA),RI(NBA,NBA),
     &           RIJ(NTR,3),CL(NOC,NOC),QPIX(NOC),QPJX(NOC),B(NBA,NOC)
#include <priunit.h>
 
C PRINT THE DIPOLE MOMENT MATRICES
 
      IPRT = 0
      IF (IPRT .GE. 1) THEN
         CALL AROUND('DIPX')
         CALL OUTPAK(DIPXYZ(1,1),NBA,1,LUPRI) 
         CALL AROUND('DIPY')
         CALL OUTPAK(DIPXYZ(1,2),NBA,1,LUPRI)
         CALL AROUND('DIPZ')
         CALL OUTPAK(DIPXYZ(1,3),NBA,1,LUPRI)
      END IF
 
C END OF THE PRINT
 
      IF (NOC.EQ.0) GO TO 500
 
C Tools to read the triangular density matrix
 
      DO I= 1, NOC
         ILIFQ(I)=(I-1)*NBA
      END DO
      DO I= 1, NBA
         IKY(I)=((I-1)*I)/2
         DO J= 1, NBA
            IF (J.le.I) THEN
               MAXI0(I,J)=I
               MINI0(I,J)=J
            ELSE
               MAXI0(I,J)=J
               MINI0(I,J)=I
            END IF
         END DO
      END DO
 
 
C     Calculate dipole moment coordinates for the molecular orbitals I,J
C     (I/COOR/J) ; COOR=X,Y,Z
C     RI are diagonal elements
C     RIJ are the off diagonal
  
      DO IC = 1, 3
         MM=0
         DO I = 1, NOC
            DO J = 1, I   
               IF (I.NE.J) MM = MM + 1
               SUM = 0.0D0
               DO K = 1, NBA        
                  KI = K + ILIFQ(I) 
                  CCMO = CMO(KI)
                  DO L = 1, NBA
                     LJ = L + ILIFQ(J)
                     KL = IKY(MAXI0(K,L)) + MINI0(K,L)
                     SUM = SUM + CCMO*CMO(LJ)*DIPXYZ(KL,IC)
                  END DO
               END DO
               IF (I.EQ.J) THEN
                  RI(I,IC) = SUM
               ELSE
                  RIJ(MM,IC) = SUM
               END IF
            END DO
         END DO
      END DO
 
C Now do the rotations
C ---------------------
C Initialize the array CL

      DO I = 1 , NOC
         DO J = 1 , I
            CL(I,J) = D0     
            CL(J,I) = D0     
         END DO  
         CL(I,I) = D1     
      END DO  
      ITER = 0
      SHIFT = DATAN(D1)      
 100  CHANGE = D0
      ITER = ITER + 1

C     2X2 unitary transformation
C         psi prime(i)=   cos(t)*psi(i)+sin(t)*psi(j) 
C         psi prime (j)= -sin(t)*psi(i)+cos(t)*psi(j).
C     BOYS LOCALIZATION : maximize the sum of the squares of molecular
C     dipole moment integrals.

      DO 150 I = 1 , NOC 
         IM = I - 1
         JM = 1
         IJM = IM*(IM-1)/2 + 1
         RM = D0
         TM = D0
         SM = D0
         CM = D1
         DO 120 J = 1 , NOC 
            IF (I.LT.J) THEN
               IJ = (J-1)*(J-2)/2 + I
            ELSE IF (I.EQ.J) THEN
               GO TO 120
            ELSE
               IJ = IM*(IM-1)/2 + J
            END IF
            T = D0
            TX = D0
            DO KK = 1 , 3
               T = T + D4*RIJ(IJ,KK)**2 - RI(I,KK)**2 - RI(J,KK)**2
     &               + D2*RI(I,KK)*RI(J,KK)
               TX = TX + RIJ(IJ,KK)*(RI(J,KK)-RI(I,KK))
            END DO   
            IF ((DABS(T).LE.THR10).AND.(DABS(TX).LE.THR10))
     &          GO TO 120
            TX = D4*TX
            T = DATAN2(TX,T)
            T = T*DP4
            SIGN = D1
            IF (T.GT.D0) SIGN = D1M
            T = T + SIGN*SHIFT
            ITIM = 0
 110        S = DSIN(T)
            ITIM = ITIM + 1
            CO = DCOS(T)
            RIN = D0
            DO KK = 1 , 3
               QPI = CO*CO*RI(I,KK) + S*S*RI(J,KK) +
     &               D2*CO*S*RIJ(IJ,KK)
               QPJ = CO*CO*RI(J,KK) + S*S*RI(I,KK) +
     &               D2M*CO*S*RIJ(IJ,KK)
               RIN = RIN + QPI*QPI + QPJ*QPJ - RI(I,KK)**2 - RI(J,KK)**2
            END DO   
            TTEST = DABS(T) - SHIFT
            IF ((DABS(T).GT.THR8).AND.(DABS(TTEST).GT.THR8)) THEN
               IF (RIN.LT.THR8M) THEN
                  IF (ITIM.LE.1) THEN
                     SIGN = D1
                     IF (T.GT.D0) SIGN = D1M     
                     T = T + SHIFT*SIGN
                     GO TO 110
                  ELSE
                     GO TO 160
                  END IF
               END IF
            END IF
            IF (RIN.GT.RM) THEN
               IJM = IJ
               RM = RIN
               JM = J
               SM = S
               CM = CO
               TM = T
            END IF
 120     CONTINUE
         T = TM
         RIN = RM
         S = SM
         CO = CM
         J = JM
         IJ = IJM
         CHANGE = CHANGE + T*T
         DO KK = 1 , 3
            QPI = CO*CO*RI(I,KK) + S*S*RI(J,KK) + D2*CO*S*RIJ(IJ,KK)
            QPJ = CO*CO*RI(J,KK) + S*S*RI(I,KK) + D2M*CO*S*RIJ(IJ,KK)
            QPIJ = (CO*CO-S*S)*RIJ(IJ,KK) + CO*S*(RI(J,KK)-RI(I,KK))
            DO 130 K = 1 , NOC 
               IF (I.LT.K) THEN
                  IK = (K-1)*(K-2)/2 + I
               ELSE IF (I.EQ.K) THEN
                  GO TO 130
               ELSE
                  IK = (I-1)*(I-2)/2 + K
               END IF
               IF (J.LT.K) THEN
                  JK = (K-1)*(K-2)/2 + J
               ELSE IF (J.EQ.K) THEN
                  GO TO 130
               ELSE
                  JK = (J-1)*(J-2)/2 + K
               END IF
               QPIX(K) = CO*RIJ(IK,KK) + S*RIJ(JK,KK)
               QPJX(K) = CO*RIJ(JK,KK) - S*RIJ(IK,KK)
 130        CONTINUE
            DO 140 K = 1 , NOC 
               IF (I.LT.K) THEN
                  IK = (K-1)*(K-2)/2 + I
               ELSE IF (I.EQ.K) THEN
                  GO TO 140
               ELSE
                  IK = (I-1)*(I-2)/2 + K
               END IF
               IF (J.LT.K) THEN
                  JK = (K-1)*(K-2)/2 + J
               ELSE IF (J.EQ.K) THEN
                  GO TO 140
               ELSE
                  JK = (J-1)*(J-2)/2 + K
               END IF
               RIJ(IK,KK) = QPIX(K)
               RIJ(JK,KK) = QPJX(K)
 140        CONTINUE
            RIN = RIN + QPI + QPJ - RI(I,KK) - RI(J,KK)
            RI(I,KK) = QPI
            RI(J,KK) = QPJ
            RIJ(IJ,KK) = QPIJ
         END DO   
         DO K = 1,NOC 
            C1 = CO*CL(K,I)+S*CL(K,J)
            C2 = -S*CL(K,I)+CO*CL(K,J)
            CL(K,I) = C1
            CL(K,J) = C2
         END DO   
 150  CONTINUE

C      if convergence has not been reached start another series
C      of two center rotations.

      CHANGE = DSQRT(D2*CHANGE/(NOC*(NOC-1)))
      IF (ITER.LE.75) THEN
         IF (CHANGE.GE.THR10) GO TO 100
      END IF
 160  WRITE(LUPRI,2200) ITER
      IF (ITER.GE.75 .OR. CHANGE.GT.THR6) THEN
      WRITE(LUPRI,2300) 
      END IF

C *** load localized molecular orbitals in CMO array

      DO I = 1 , NOC
         CALL DCOPY(NOC,CL(1,I),1,B(1,I),1)
      END DO   

      DO I = 1 , NOC 
C 
C        (WK/UniKA/07-11-2003)
C------> CALL VCLR(RI(1,1),1,NBA)
         CALL DZERO(RI,NBA)
C
         DO J = 1 , NOC  
            IIJ = ILIFQ(J)
            CALL DAXPY(NBA,B(J,I),CMO(IIJ+1),1,RI(1,1),1)
         END DO    
         CALL DCOPY(NBA,RI(1,1),1,B(1,I),1)
      END DO   

      DO I = 1 , NOC
         II = ILIFQ(I)
         CALL DCOPY(NBA,B(1,I),1,CMO(II+1),1)
      END DO   

C End of the Localization

 500  CONTINUE
      RETURN
 2200 FORMAT(/9x,'** LOCALIZATION CONVERGED AFTER',i3,' ITERATIONS **'/)
 2300 FORMAT(/'** LOCALIZATION HAS BEEN UNSUCESSFUL **')

      END

C  /* Deck proao  */
      SUBROUTINE PROAO(PAO,EVC,CMO,OMOC,OVLAT,OVLA,NOC,NBA,NB1,NNBA)
C     Written by Claire C.M. Samson (University of Karlsruhe, 28 April 2003).
#include <implicit.h>
      PARAMETER ( D0 = 0.0D0, D1 = 1.0D0, D1M = -1.0D0 )
      DIMENSION CMO(NBA,NOC),OMOC(NB1,NOC)
      DIMENSION OVLAT(NNBA),OVLA(NB1,NB1)
      DIMENSION EVC(NB1,NB1),PAO(NB1,NB1)

      DO J = 1, NOC
         DO I = 1, NB1
            OMOC(I,J)=CMO(I,J)
         END DO
      END DO

      CALL DZERO(PAO,NB1*NB1)
      IJ = 0
      DO I = 1, NB1
         PAO(I,I) = D1
         DO J = 1, I
            IJ = IJ + 1
            OVLA(I,J) = OVLAT(IJ)
            OVLA(J,I) = OVLAT(IJ)
         END DO
      END DO

      CALL DGEMM('N','T',NB1,NB1,NOC,D1,OMOC,NB1,OMOC,NB1,D0,EVC,NB1)
      CALL DGEMM('N','N',NB1,NB1,NB1,D1M,EVC,NB1,OVLA,NB1,D1,PAO,NB1)
      DO I=1,NB1
       CALL DGEMV('N',NB1,NB1,D1,OVLA,NB1,PAO(1,I),1,D0,EVC,1)
       XNO=DDOT(NB1,PAO(1,I),1,EVC,1)
       XNO=D1/DSQRT(XNO)
       CALL DSCAL(NB1,XNO,PAO(1,I),1)
      ENDDO
C test
C      CALL DGEMM('T','N',NB1,NB1,NB1,D1,PAO,NB1,OVLA,NB1,D0,EVC,NB1)
C      CALL DGEMM('N','N',NB1,NOC,NB1,D1,EVC,NB1,OMOC,NB1,D0,OMOC,NB1)
      
      RETURN
      END

C  /* Deck wrlocmo  */
      SUBROUTINE WRLOCMO(CMO,FKMAT,NOC,NCO,NBA,NB1)
C     Written by Claire C.M. Samson (University of Karlsruhe, 28 April 2003).
#include <implicit.h>
#include <priunit.h>
      DIMENSION FKMAT(*), CMO(*)

      IJ = 0
      OPEN(99,FILE='FLOCA')
      DO I = 1, NB1+NOC
         DO J = 1, I
            IJ = IJ + 1
            IF (J.GT.NCO) WRITE(99,'(D30.20)') FKMAT(IJ)
         END DO
      END DO
      CLOSE(99)
C
      OPEN(99,FILE='LOCMO',FORM='FORMATTED')
      WRITE(99,'(D30.20)') (CMO(IJ), IJ = 1, NBA * (NB1+NOC))
      CLOSE(99)
      RETURN
      END
C  /* Deck printorb  */
      SUBROUTINE PRINTORB(NENDI,NBASI,ARRA)
C     Written by Claire C.M. Samson (University of Karlsruhe, 28 April 2003).
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <maxash.h>
#include <infinp.h>
#include <inforb.h>
      DIMENSION ARRA(*)
C
      ISTBAS = 0
      IF (NENDI.EQ.0) GO TO 20
      ICMOI = ICMO(1)
      ISTORB = IORB(1)
      IEND=0
   10 IST =IEND+1
      ISTMO=IEND*NBASI+ICMOI
      IEND=IEND+7
      IF(IEND.GT.NENDI) IEND=NENDI
      IEMO=NBASI*(IEND-1)+ICMOI
      WRITE(LUPRI,3100) (I,I=IST,IEND)
      DO I=1,NBASI
         JSMO=ISTMO+I
         JEMO=IEMO+I
         WRITE(LUPRI,3200) I,CENT(I+ISTBAS),TYPE(I+ISTBAS),
     *                  (ARRA(J),J=JSMO,JEMO,NBASI)
      END DO
      IF (IEND.NE.NENDI) GO TO 10
   20 CONTINUE
      RETURN
 3100 FORMAT(/' Orbital  ',7I9)
 3200 FORMAT(1X,I3,2X,2A4,7F9.4)
      END
