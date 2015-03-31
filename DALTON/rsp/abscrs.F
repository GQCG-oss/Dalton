!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C ============================================================================
C     newabs.F: Main author Joanna Kauczor
C     This implementation is published in:
C
C     Reference C:
C        J. Kauczor, P. Jorgensen, and P. Norman, 
C        "On the efficiency of algorithms for solving Hartree-Fock and Kohn-Sham
C        response equations",
C        J. Chem. Theory Comput. 7 (2011) 1610.
C ============================================================================
C
      SUBROUTINE ABS_GAMMA_SETUP
C                                                                                             
C ======================================================================                      
C                                                                                             
C  Purpose : determine what is to be calculated for the requested CR                          
C                                                                                             
C   Author : Tobias Fahleson, 2014                                                            
C                                                                                             
C ======================================================================                      
C                                                                                             
#include "implicit.h"
#include "priunit.h"
#include "inforb.h"
#include "abslrs.h"
#include "abscrs.h"
C                                                                                             
      PARAMETER (THRZERO = 1.0D-6)
      LOGICAL DOHYP
      CHARACTER*8 ALAB,BLAB,CLAB,DLAB
C                                                                                             
      ABS_NCRF = 0
      ABS_CRF_NEQTO = 0
C                                                                                             
C     Form the three lists (acting as starting points for the CR                              
C     calculations) ABS_CRFLAB(ABS_NCRF,4), ABS_CRFSYM(ABS_NCRF,4),                           
C     ABS_CRFFRQ(ABS_NCRF,4), containing operator labels, symmetries of                       
C     said operators, and frequencies, for each ABS_NCRF number of CRF.                       
C     =================================================================                       
C                                                                                             
C     For every possible frequency combination                                                
      DO IBFR=1,ABS_NFREQ_GAMMA_B
       DO ICFR=1,ABS_NFREQ_GAMMA_C
        DO IDFR=1,ABS_NFREQ_GAMMA_D
           FREQB = ABS_FREQ_GAMMA_B(IBFR)
           FREQC = ABS_FREQ_GAMMA_C(ICFR)
           FREQD = ABS_FREQ_GAMMA_D(IDFR)
           FREQA = -(FREQB+FREQC+FREQD)
C
C          Set the specific symmetries for operators A to D                                   
           DO ISYMA=1,NSYM
            DO ISYMB=1,NSYM
             DO ISYMC=1,NSYM
                ISYMD = MULD2H(ISYMC, MULD2H(ISYMA, ISYMB))
C                                                                                             
C               If the current symmetries represents a nonzero number of                      
C               operators                                                                     
                IF (ABS_NOPER(ISYMA) .GT. 0 .AND. ABS_NOPER(ISYMB) .GT.
     &             0 .AND. ABS_NOPER(ISYMC) .GT. 0 .AND.
     &             ABS_NOPER(ISYMD).GT.0) THEN
C                                                                                             
C                  For every possible operator combination                                    
                   DO IAOP=1,ABS_NOPER(ISYMA)
                    DO IBOP=1,ABS_NOPER(ISYMB)
                     DO ICOP=1,ABS_NOPER(ISYMC)
                      DO IDOP=1,ABS_NOPER(ISYMD)
C                                                                                             
C                        Operator labels are set depending on symmetry                        
                         ALAB = ABS_LABOP(IAOP, ISYMA)
                         BLAB = ABS_LABOP(IBOP, ISYMB)
                         CLAB = ABS_LABOP(ICOP, ISYMC)
                         DLAB = ABS_LABOP(IDOP, ISYMD)
C                                                                                             
C                        Check if this operator combination at certain                        
C                        frequencies and with certain symmetries should                       
C                        be included in the CR calculation                                    
                         CALL ABS_CRCHK(DOHYP,
     &                                  ALAB,BLAB,CLAB,DLAB,
     &                                  ISYMA,ISYMB,ISYMC,ISYMD,
     &                                  FREQA,FREQB,FREQC,FREQD)
C                                                                                             
C                        If the operator combination checks out, append                       
C                        the labels, symmetries, and frequencies, to the                      
C                        lists                                                                
                         IF (DOHYP) THEN
                            ABS_NCRF = ABS_NCRF + 1
                            ABS_CRFLAB(ABS_NCRF,1) = ALAB
                            ABS_CRFLAB(ABS_NCRF,2) = BLAB
                            ABS_CRFLAB(ABS_NCRF,3) = CLAB
                            ABS_CRFLAB(ABS_NCRF,4) = DLAB
                            ABS_CRFSYM(ABS_NCRF,1) = ISYMA
                            ABS_CRFSYM(ABS_NCRF,2) = ISYMB
                            ABS_CRFSYM(ABS_NCRF,3) = ISYMC
                            ABS_CRFSYM(ABS_NCRF,4) = ISYMD
                            ABS_CRFFRQ(ABS_NCRF,1) = FREQA
                            ABS_CRFFRQ(ABS_NCRF,2) = FREQB
                            ABS_CRFFRQ(ABS_NCRF,3) = FREQC
                            ABS_CRFFRQ(ABS_NCRF,4) = FREQD
                         END IF
C                                                                                             
                      END DO
                     END DO
                    END DO
                   END DO
                END IF
             END DO
            END DO
           END DO
 100    CONTINUE
        END DO
       END DO
      END DO
C                                                                                             
C     Sort list of LR frequencies while removing                                              
C     overlapping frequencies separated by less than                                          
C     THRZERO                                                                                 
C     ==============================================                                          
C                                                                                             
      CALL GPDSRT(ABS_NFREQ_ALPHA, ABS_FREQ_ALPHA, THRZERO)
C                                                                                             
C     Print what has been determined by the gamma setup                                       
C     =================================================                                       
C                                                                                             
      CALL AROUND('Setup of Second-Order Hyperpolarizability '//
     &            'Calculation')
C                                                                                             
      WRITE (LUPRI,'(2(/A),I4,A)')
     & ' This calculations requires the solution of linear response',
     & ' equations for electric dipole operators at ',ABS_NFREQ_GAMMA_B,
     & ' frequencies:'
      WRITE(LUPRI,'(/A,5(4F12.8,/,9X))')
     &        ' LR FREQ:',(ABS_FREQ_GAMMA_B(I),I=1,ABS_NFREQ_GAMMA_B)
      WRITE (LUPRI,'(/A,I4,A)')
     & ' and the evaluation of ',ABS_NCRF,
     & ' cubic response functions:'
      WRITE(LUPRI,'(/2A,/A3,8A12,/2A)')
     & '--------------------------------------------------',
     & '-------------------------------------------------',
     & ' No','A-oper','B-oper','C-oper','D-oper',
     & 'A-freq','B-freq','C-freq','D-freq',
     & '--------------------------------------------------',
     & '-------------------------------------------------'
      DO ICRF=1,ABS_NCRF
         WRITE(LUPRI,'(I3,4A12,4F12.8)') ICRF,
     &        (ABS_CRFLAB(ICRF,I), I=1,4),
     &        (ABS_CRFFRQ(ICRF,I), I=1,4)
      END DO
      WRITE(LUPRI,'(2A)')
     & '--------------------------------------------------',
     & '-------------------------------------------------'
C                                                                                             
C     End of ABS_GAMMA_SETUP                                                                  
C                                                                                             
      RETURN
      END
C                                                                                             
C                                                                                             
C
      SUBROUTINE ABS_CRCHK(DOHYP,
     &                     ALAB,BLAB,CLAB,DLAB,
     &                     ISYMA,ISYMB,ISYMC,ISYMD,
     &                     FREQA,FREQB,FREQC,FREQD)
C                                                                                             
C ======================================================================                      
C                                                                                             
C  Purpose : Determine if the supplied combination of operators and                           
C            corresponding symmetries and frequencies should be included                      
C            in the CR list                                                                   
C                                                                                             
C   Author : Tobias Fahleson, 2014                                                            
C                                                                                             
C ======================================================================                      
C                                                                                             
#include "implicit.h"
#include "priunit.h"
#include "abslrs.h"
#include "abscrs.h"
C                                                                                             
      PARAMETER (THRZERO = 1.0D-6)
      LOGICAL NEWFRQ, DOHYP, DOHYPLAB, DOHYPFREQ, DONXY
      DIMENSION FREQ(4)
      CHARACTER*8 ALAB, BLAB, CLAB, DLAB, LAB(4)

C     Nonoverlapping IDRI elements required for orientational avarage---                      
C     not considering symmetry of the system which might reduce the list                      
C     even further
      CHARACTER*4 IDRI_ELEMENTS(15)/'XXXX','XXYY','XYXY','XXZZ','XZXZ',
     &                              'YYXX','YXYX','YYYY','YYZZ','YZYZ',
     &                              'ZZXX','ZXZX','ZZYY','ZYZY','ZZZZ'/
C                                                                                             
C     Nonoverlapping required NXY vectors for IDRI evaluation where the                       
C     first row is for (omega, +/-omega) and the second row is only                           
C     for (omega, -omega)                                                                     
      CHARACTER*2 IDRI_NXY_ELEMENTS(9)/'XX','XY','XZ','YY','YZ','ZZ',
     &                                 'YX','ZX','ZY'/
C                                                                                             
      DOHYP     = .FALSE.
      DOHYPLAB  = .FALSE.
      DOHYPFREQ = .FALSE.
      DONXY     = .TRUE.
C                                                                                             
C     Skip the current frequency combination if IDRI has been                                 
C     requested and the frequency requiremens are not met                                     
      IF (ABSLRS_IDRI .AND. .NOT.
     &     (FREQB .EQ. -FREQC .AND. FREQB .EQ. FREQD)) THEN
         DOHYPFREQ = .FALSE.
      END IF
C                                                                                             
C     If IDRI has been requested, check for current element in the list                       
C     of relevant IDRI elements                                                               
C     =================================================================                       
C                                                                                             
      IIDRI_ELEMENT = 0
      IF (ABSLRS_IDRI) THEN
         DO IIE=1,15
            IF ( ALAB(1:1) .EQ. IDRI_ELEMENTS(IIE)(1:1) .AND.
     &           BLAB(1:1) .EQ. IDRI_ELEMENTS(IIE)(2:2) .AND.
     &           CLAB(1:1) .EQ. IDRI_ELEMENTS(IIE)(3:3) .AND.
     &           DLAB(1:1) .EQ. IDRI_ELEMENTS(IIE)(4:4) ) THEN
               DOHYPLAB = .TRUE.
C                                                                                             
C              If an IDRI element matched, save its index for                                 
C              later reference                                                                
               IIDRI_ELEMENT = IIE
            END IF
         END DO
      END IF
C                                                                                             
C     If GAMMA has been requested, check if only tensor-average                               
C     elements, all elements, or specific elements, are wanted                                
C     =========================================================                               
C                                                                                             
C     Allow only elements contributing to average GAMMA                                       
      IF (ABSLRS_GAMMA .AND. .NOT. ABS_GAMMA_ALLELE
     &                 .AND. .NOT. ABS_GAMMA_ELEMENT) THEN
         IF ((ALAB.EQ.BLAB .OR. ALAB.EQ.CLAB .OR. ALAB.EQ.DLAB) .AND.
     &       (BLAB.EQ.ALAB .OR. BLAB.EQ.CLAB .OR. BLAB.EQ.DLAB) .AND.
     &       (CLAB.EQ.ALAB .OR. CLAB.EQ.BLAB .OR. CLAB.EQ.DLAB) .AND.
     &       (DLAB.EQ.ALAB .OR. DLAB.EQ.BLAB .OR. DLAB.EQ.CLAB)) THEN
            DOHYPLAB = .TRUE.
         END IF
C                                                                                             
C     Allow all GAMMA elements                                                                
      ELSEIF (ABSLRS_GAMMA .AND. ABS_GAMMA_ALLELE) THEN
         DOHYPLAB = .TRUE.
C                                                                                             
C     Allow only specified GAMMA elements                                                     
      ELSEIF (ABSLRS_GAMMA .AND. ABS_GAMMA_ELEMENT) THEN
         DO IELEMENT=1,ABS_NELEMENTS_GAMMA
           IF (ALAB(1:1) .EQ. ABS_ELEMENTS_GAMMA(IELEMENT)(1:1)  .AND.
     &         BLAB(1:1) .EQ. ABS_ELEMENTS_GAMMA(IELEMENT)(2:2)  .AND.
     &         CLAB(1:1) .EQ. ABS_ELEMENTS_GAMMA(IELEMENT)(3:3)  .AND.
     &         DLAB(1:1) .EQ. ABS_ELEMENTS_GAMMA(IELEMENT)(4:4)) THEN
               DOHYPLAB = .TRUE.
           END IF
         END DO
      END IF
C                                                                                             
C     If one-to-one frequencies has been requested for GAMMA or IDRI                          
C     through '.OTOFGA' keyword, check if current frequency combinations                      
C     match                                                                                   
C     ==================================================================                      
C                                                                                             
      IF ((ABSLRS_GAMMA .OR. ABSLRS_IDRI) .AND. ABS_OTOF_GAMMA) THEN
         IF (ABS_NFREQ_GAMMA_B .NE. ABS_NFREQ_GAMMA_C  .AND.
     &       ABS_NFREQ_GAMMA_B .NE. ABS_NFREQ_GAMMA_D) THEN
             CALL QUIT('One-to-one frequencies requested for                                  
     &       cubic response, while number of frequencies not matched')
         END IF
         DO IOTOFREQ = 1,ABS_NFREQ_GAMMA_B
            IF (FREQB .EQ. ABS_FREQ_GAMMA_B(IOTOFREQ)  .AND.
     &          FREQC .EQ. ABS_FREQ_GAMMA_C(IOTOFREQ)  .AND.
     &          FREQD .EQ. ABS_FREQ_GAMMA_D(IOTOFREQ)) THEN
                DOHYPFREQ = .TRUE.
            END IF
         END DO
      ELSEIF ((ABSLRS_GAMMA .OR. ABSLRS_IDRI) .AND.
     &         .NOT. ABS_OTOF_GAMMA) THEN
         DOHYPFREQ = .TRUE.
      END IF
C                                                                                             
C     Check if labels and frequencies made it through the                                     
C     test above for GAMMA.                                                                   
C     ===================================================                                     
C                                                                                             
      IF (DOHYPLAB .AND. DOHYPFREQ) THEN
         DOHYP = .TRUE.
      ELSE
         DOHYP = .FALSE.
         RETURN
      END IF
C                                                                                             
C     Check if an equivalent CRF (permutational symmetry) is already                          
C     in the lists                                                                            
C     ==============================================================                          
C
      IF (ABS_DAMP .EQ. 0.0D0) THEN
C                                                                                             
C        Overall permutational symmetry                                                       
         JCTR=4
         KCTR=1
         LCTR=1
         MCTR=1
      ELSE
C                                                                                             
C        Only intrinsic permutational symmetry, so operator A has to                          
C        match the first operator in the list of response functions                           
         JCTR=1
         KCTR=2
         LCTR=2
         MCTR=2
      END IF
C                                                                                             
      DO ICRF = 1,ABS_NCRF
       DO J = 1,JCTR
        DO K = KCTR,4
         DO L = LCTR,4
          DO M = MCTR,4
C                                                                                             
             IF ( J.NE.K .AND. J.NE.L .AND. J.NE.M .AND.
     &            (K.NE.L .AND. K.NE.M .AND.
     &            (L.NE.M)) ) THEN
C                                                                                             
C               Go through the already accepted CRF elements and see                          
C               if the current element is identical to some other                             
C               element by permutational symmetry                                             
                IF ( ALAB  .EQ. ABS_CRFLAB(ICRF,J) .AND.
     &               BLAB  .EQ. ABS_CRFLAB(ICRF,K) .AND.
     &               CLAB  .EQ. ABS_CRFLAB(ICRF,L) .AND.
     &               DLAB  .EQ. ABS_CRFLAB(ICRF,M) .AND.
     &               ISYMA .EQ. ABS_CRFSYM(ICRF,J) .AND.
     &               ISYMB .EQ. ABS_CRFSYM(ICRF,K) .AND.
     &               ISYMC .EQ. ABS_CRFSYM(ICRF,L) .AND.
     &               ISYMD .EQ. ABS_CRFSYM(ICRF,M) .AND.
     &               ABS(FREQA - ABS_CRFFRQ(ICRF,J)) .LT. THRZERO .AND.
     &               ABS(FREQB - ABS_CRFFRQ(ICRF,K)) .LT. THRZERO .AND.
     &               ABS(FREQC - ABS_CRFFRQ(ICRF,L)) .LT. THRZERO .AND.
     &               ABS(FREQD - ABS_CRFFRQ(ICRF,M)) .LT. THRZERO) THEN
                   DOHYP = .FALSE.
C                                                                                             
C                  If so, store the current element's information and                         
C                  link it to the matching ICRF number
                   ABS_CRF_NEQTO = ABS_CRF_NEQTO + 1
                   ABS_CRF_EQTO(ABS_CRF_NEQTO) = ICRF
                   ABS_CRFLAB_EQTO(ABS_CRF_NEQTO,1) = ALAB
                   ABS_CRFLAB_EQTO(ABS_CRF_NEQTO,2) = BLAB
                   ABS_CRFLAB_EQTO(ABS_CRF_NEQTO,3) = CLAB
                   ABS_CRFLAB_EQTO(ABS_CRF_NEQTO,4) = DLAB
                   ABS_CRFSYM_EQTO(ABS_CRF_NEQTO,1) = ISYMA
                   ABS_CRFSYM_EQTO(ABS_CRF_NEQTO,2) = ISYMB
                   ABS_CRFSYM_EQTO(ABS_CRF_NEQTO,3) = ISYMC
                   ABS_CRFSYM_EQTO(ABS_CRF_NEQTO,4) = ISYMD
                   ABS_CRFFRQ_EQTO(ABS_CRF_NEQTO,1) = FREQA
                   ABS_CRFFRQ_EQTO(ABS_CRF_NEQTO,2) = FREQB
                   ABS_CRFFRQ_EQTO(ABS_CRF_NEQTO,3) = FREQC
                   ABS_CRFFRQ_EQTO(ABS_CRF_NEQTO,4) = FREQD
                   RETURN
                END IF
C                                                                                             
             END IF
C                                                                                             
          END DO
         END DO
        END DO
       END DO
      END DO
C      END IF                                                                                 
C                                                                                             
C     Check if this CRF will inflict new LR solver frequencies                                
C     ========================================================                                
C                                                                                             
      IF (DOHYP) THEN
C                                                                                             
         FREQ(1) = FREQA
         FREQ(2) = FREQB
         FREQ(3) = FREQC
         FREQ(4) = FREQD
C                                                                                             
C        For every frequency A to D in current CRF                                            
         DO I = 1,4
            NEWFRQ = .TRUE.
C                                                                                             
C           For every LR solver frequency                                                     
            DO IFR = 1,ABS_NFREQ_ALPHA
C                                                                                             
C              Check if the difference between the CR                                         
C              frequency and some LR frequency is negligible, i.e. the                        
C              LR list already has it or not                                                  
               IF (ABS(ABS_FREQ_ALPHA(IFR)-ABS(FREQ(I))) .LT. THRZERO)
     &            THEN
                  NEWFRQ = .FALSE.
               END IF
            END DO
C                                                                                             
C           If the the CR frequency wasn't found in the LR list, append                       
C           it                                                                                
            IF (NEWFRQ) THEN
               IF (ABS_NFREQ_ALPHA .GE. NMXFREQ) THEN
                  WRITE(LUPRI,'(2(/A),I4,A,/A)')
     & ' The specified calculation requires more than the allowed',
     & ' number of frequencies in the LR solver NMXFREQ=',NMXFREQ,'.',
     & ' The program will stop.'
                  CALL QUIT('Too many frequencies in LR solver.')
               END IF
               ABS_NFREQ_ALPHA = ABS_NFREQ_ALPHA + 1
               ABS_FREQ_ALPHA(ABS_NFREQ_ALPHA) = ABS(FREQ(I))
            END IF
         END DO
      END IF
C                                                                                             
C     If GAMMA has been requested, add all necessary NXY vectors to                           
C     the fold                                                                                
C     =============================================================                           
C                                                                                             
      IF (ABSLRS_GAMMA .AND. DOHYP) THEN
C                                                                                             
C        Check for operator combinations B and C                                              
C        =======================================                                              
C                                                                                             
C        Check if this NXY is already in the list                                             
         IF (ABS_NXYVEC .GT. 0) THEN
            DONXY = .TRUE.
            DO ICHK = 1,ABS_NXYVEC
               IF ( BLAB  .EQ. ABS_XYLAB(ICHK,1) .AND.
     &              CLAB  .EQ. ABS_XYLAB(ICHK,2) .AND.
     &              ISYMB .EQ. ABS_XYSYM(ICHK,1) .AND.
     &              ISYMC .EQ. ABS_XYSYM(ICHK,2) .AND.
     &              ABS(FREQB - ABS_XYFREQ(ICHK,1)) .LT. THRZERO .AND.
     &              ABS(FREQC - ABS_XYFREQ(ICHK,2)) .LT. THRZERO )
     &              THEN
                  DONXY = .FALSE.
               END IF
            END DO
         END IF
C                                                                                             
C        If current NXY is unique so far, add it to the list
         IF (DONXY) THEN
            ABS_NXYVEC = ABS_NXYVEC + 1
            ABS_XYFREQ(ABS_NXYVEC,1) = FREQB
            ABS_XYFREQ(ABS_NXYVEC,2) = FREQC
            ABS_XYLAB( ABS_NXYVEC,1) = BLAB
            ABS_XYLAB( ABS_NXYVEC,2) = CLAB
            ABS_XYSYM( ABS_NXYVEC,1) = ISYMB
            ABS_XYSYM( ABS_NXYVEC,2) = ISYMC
         END IF
C                                                                                             
C        Check for operator combinations B and D                                              
C        =======================================                                              
C                                                                                             
C        Check if this NXY is already in the list                                             
         IF (ABS_NXYVEC .GT. 0) THEN
            DONXY = .TRUE.
            DO ICHK = 1,ABS_NXYVEC
               IF ( BLAB  .EQ. ABS_XYLAB(ICHK,1) .AND.
     &              DLAB  .EQ. ABS_XYLAB(ICHK,2) .AND.
     &              ISYMB .EQ. ABS_XYSYM(ICHK,1) .AND.
     &              ISYMD .EQ. ABS_XYSYM(ICHK,2) .AND.
     &              ABS(FREQB - ABS_XYFREQ(ICHK,1)) .LT. THRZERO .AND.
     &              ABS(FREQD - ABS_XYFREQ(ICHK,2)) .LT. THRZERO )
     &              THEN
                  DONXY = .FALSE.
               END IF
            END DO
         END IF
C                                                                                             
C        If current NXY is unique so far, add it to the list                                  
         IF (DONXY) THEN
            ABS_NXYVEC = ABS_NXYVEC + 1
            ABS_XYFREQ(ABS_NXYVEC,1) = FREQB
            ABS_XYFREQ(ABS_NXYVEC,2) = FREQD
            ABS_XYLAB( ABS_NXYVEC,1) = BLAB
            ABS_XYLAB( ABS_NXYVEC,2) = DLAB
            ABS_XYSYM( ABS_NXYVEC,1) = ISYMB
            ABS_XYSYM( ABS_NXYVEC,2) = ISYMD
         END IF
C                                                                                             
C        Check for operator combinations C and D                                              
C        =======================================                                              
C                                                                                             
C        Check if this NXY is already in the list
         IF (ABS_NXYVEC .GT. 0) THEN
            DONXY = .TRUE.
            DO ICHK = 1,ABS_NXYVEC
               IF ( CLAB  .EQ. ABS_XYLAB(ICHK,1) .AND.
     &              DLAB  .EQ. ABS_XYLAB(ICHK,2) .AND.
     &              ISYMC .EQ. ABS_XYSYM(ICHK,1) .AND.
     &              ISYMD .EQ. ABS_XYSYM(ICHK,2) .AND.
     &              ABS(FREQC - ABS_XYFREQ(ICHK,1)) .LT. THRZERO .AND.
     &              ABS(FREQD - ABS_XYFREQ(ICHK,2)) .LT. THRZERO )
     &              THEN
                  DONXY = .FALSE.
               END IF
            END DO
         END IF
C                                                                                             
C        If current NXY is unique so far, add it to the list                                  
         IF (DONXY) THEN
            ABS_NXYVEC = ABS_NXYVEC + 1
            ABS_XYFREQ(ABS_NXYVEC,1) = FREQC
            ABS_XYFREQ(ABS_NXYVEC,2) = FREQD
            ABS_XYLAB( ABS_NXYVEC,1) = CLAB
            ABS_XYLAB( ABS_NXYVEC,2) = DLAB
            ABS_XYSYM( ABS_NXYVEC,1) = ISYMC
            ABS_XYSYM( ABS_NXYVEC,2) = ISYMD
         END IF
      END IF
C                                                                                             
C     Check all the elements in the table IDRI_NXY_ELEMENTS,                                  
C     containing necessary NXY type vectors, and see if some match                            
C     the current CRF operator combination of BLAB, CLAB, and DLAB                            
C     ============================================================                            
C                                                                                             
      IF (ABSLRS_IDRI .AND. DOHYP) THEN
C                                                                                             
C     Check the (omega,-omega) case for operator combinations B and C                         
C     ===============================================================                         
C
      DO I = 1,9
         IF ( IDRI_ELEMENTS(IIDRI_ELEMENT)(2:3) .EQ.
     &        IDRI_NXY_ELEMENTS(I) ) THEN
C                                                                                             
C           Check if this NXY is already in the list                                          
            IF (ABS_NXYVEC .GT. 0) THEN
               DO ICHK = 1,ABS_NXYVEC
                  IF ( BLAB  .EQ. ABS_XYLAB(ICHK,1) .AND.
     &                 CLAB  .EQ. ABS_XYLAB(ICHK,2) .AND.
     &                 ISYMB .EQ. ABS_XYSYM(ICHK,1) .AND.
     &                 ISYMC .EQ. ABS_XYSYM(ICHK,2) .AND.
     &            ABS( FREQB - ABS_XYFREQ(ICHK,1)).LT.THRZERO.AND.
     &            ABS(-FREQB - ABS_XYFREQ(ICHK,2)).LT.THRZERO )
     &            THEN
                     DONXY = .FALSE.
                  END IF
               END DO
            END IF
C                                                                                             
C           If current NXY is unique so far, add it to the list
            IF (DONXY) THEN
               ABS_NXYVEC = ABS_NXYVEC + 1
               ABS_XYFREQ(ABS_NXYVEC,1) =  FREQB
               ABS_XYFREQ(ABS_NXYVEC,2) = -FREQB
               ABS_XYLAB( ABS_NXYVEC,1) =  BLAB
               ABS_XYLAB( ABS_NXYVEC,2) =  CLAB
               ABS_XYSYM( ABS_NXYVEC,1) =  ISYMB
               ABS_XYSYM( ABS_NXYVEC,2) =  ISYMC
            END IF
         END IF
C                                                                                             
C        Reset DONXY for next turn in the loop                                                
         DONXY = .TRUE.
      END DO
C                                                                                             
C     Check the (omega,omega) case for operator combinations B and C                          
C     ==============================================================                          
C
      DO I = 1,6
         IF ( IDRI_ELEMENTS(IIDRI_ELEMENT)(2:3) .EQ.
     &        IDRI_NXY_ELEMENTS(I) ) THEN
C                                                                                             
C           Check if this NXY is already in the list                                          
            IF (ABS_NXYVEC .GT. 0) THEN
               DO ICHK = 1,ABS_NXYVEC
                  IF ( BLAB  .EQ. ABS_XYLAB(ICHK,1) .AND.
     &                 CLAB  .EQ. ABS_XYLAB(ICHK,2) .AND.
     &                 ISYMB .EQ. ABS_XYSYM(ICHK,1) .AND.
     &                 ISYMC .EQ. ABS_XYSYM(ICHK,2) .AND.
     &            ABS( FREQB - ABS_XYFREQ(ICHK,1)).LT.THRZERO.AND.
     &            ABS( FREQB - ABS_XYFREQ(ICHK,2)).LT.THRZERO )
     &            THEN
                     DONXY = .FALSE.
                  END IF
               END DO
            END IF
C                                                                                             
C           If current NXY is unique so far, add it to the list                               
            IF (DONXY) THEN
               ABS_NXYVEC = ABS_NXYVEC + 1
               ABS_XYFREQ(ABS_NXYVEC,1) =  FREQB
               ABS_XYFREQ(ABS_NXYVEC,2) =  FREQB
               ABS_XYLAB( ABS_NXYVEC,1) =  BLAB
               ABS_XYLAB( ABS_NXYVEC,2) =  CLAB
               ABS_XYSYM( ABS_NXYVEC,1) =  ISYMB
               ABS_XYSYM( ABS_NXYVEC,2) =  ISYMC
            END IF
         END IF
C                                                                                             
C        Reset DONXY for next turn in the loop                                                
         DONXY = .TRUE.
      END DO
C                                                                                             
C     Check the (omega,-omega) case for operator combinations C and D                         
C     ===============================================================                         
C
      DO I = 1,9
         IF ( IDRI_ELEMENTS(IIDRI_ELEMENT)(3:4) .EQ.
     &        IDRI_NXY_ELEMENTS(I) .AND. .NOT.
     &        IDRI_ELEMENTS(IIDRI_ELEMENT)(2:3) .EQ.
     &        IDRI_ELEMENTS(IIDRI_ELEMENT)(3:4)) THEN
C                                                                                             
C           Check if this NXY is already in the list                                          
            IF (ABS_NXYVEC .GT. 0) THEN
               DO ICHK = 1,ABS_NXYVEC
                  IF ( CLAB  .EQ. ABS_XYLAB(ICHK,1) .AND.
     &                 DLAB  .EQ. ABS_XYLAB(ICHK,2) .AND.
     &                 ISYMC .EQ. ABS_XYSYM(ICHK,1) .AND.
     &                 ISYMD .EQ. ABS_XYSYM(ICHK,2) .AND.
     &            ABS( FREQB - ABS_XYFREQ(ICHK,1)).LT.THRZERO.AND.
     &            ABS(-FREQB - ABS_XYFREQ(ICHK,2)).LT.THRZERO )
     &            THEN
                     DONXY = .FALSE.
                  END IF
               END DO
            END IF
C                                                                                             
C           If current NXY is unique so far, add it to the list                               
            IF (DONXY) THEN
               ABS_NXYVEC = ABS_NXYVEC + 1
               ABS_XYFREQ(ABS_NXYVEC,1) =  FREQB
               ABS_XYFREQ(ABS_NXYVEC,2) = -FREQB
               ABS_XYLAB( ABS_NXYVEC,1) =  CLAB
               ABS_XYLAB( ABS_NXYVEC,2) =  DLAB
               ABS_XYSYM( ABS_NXYVEC,1) =  ISYMC
               ABS_XYSYM( ABS_NXYVEC,2) =  ISYMD
            END IF
         END IF
C                                                                                             
C        Reset DONXY for next turn in the loop                                                
         DONXY = .TRUE.
      END DO
C                                                                                             
C     Check the (omega,omega) case for operator combinations C and D                          
C     ==============================================================                          
C
      DO I = 1,6
         IF ( IDRI_ELEMENTS(IIDRI_ELEMENT)(3:4) .EQ.
     &        IDRI_NXY_ELEMENTS(I) .AND. .NOT.
     &        IDRI_ELEMENTS(IIDRI_ELEMENT)(2:3) .EQ.
     &        IDRI_ELEMENTS(IIDRI_ELEMENT)(3:4)) THEN
C                                                                                             
C           Check if this NXY is already in the list                                          
            IF (ABS_NXYVEC .GT. 0) THEN
               DO ICHK = 1,ABS_NXYVEC
                  IF ( CLAB  .EQ. ABS_XYLAB(ICHK,1) .AND.
     &                 DLAB  .EQ. ABS_XYLAB(ICHK,2) .AND.
     &                 ISYMC .EQ. ABS_XYSYM(ICHK,1) .AND.
     &                 ISYMD .EQ. ABS_XYSYM(ICHK,2) .AND.
     &            ABS( FREQB - ABS_XYFREQ(ICHK,1)).LT.THRZERO.AND.
     &            ABS( FREQB - ABS_XYFREQ(ICHK,2)).LT.THRZERO )
     &            THEN
                     DONXY = .FALSE.
                  END IF
               END DO
            END IF
C                                                                                             
C           If current NXY is unique so far, add it to the list                               
            IF (DONXY) THEN
               ABS_NXYVEC = ABS_NXYVEC + 1
               ABS_XYFREQ(ABS_NXYVEC,1) =  FREQB
               ABS_XYFREQ(ABS_NXYVEC,2) =  FREQB
               ABS_XYLAB( ABS_NXYVEC,1) =  CLAB
               ABS_XYLAB( ABS_NXYVEC,2) =  DLAB
               ABS_XYSYM( ABS_NXYVEC,1) =  ISYMC
               ABS_XYSYM( ABS_NXYVEC,2) =  ISYMD
            END IF
         END IF
C                                                                                             
C        Reset DONXY for the next turn in the loop                                            
         DONXY = .TRUE.
      END DO
C                                                                                             
      END IF
C                                                                                             
C     If we are breaching the maximum number of allowed CRF, call quit                        
C     ================================================================                        
C
      IF (ABS_NCRF.GE.MXCRF .AND. DOHYP) THEN
         WRITE(LUPRI,'(2(/A),I4,A,/A)')
     & ' The specified calculation requires more than the allowed',
     & ' number of cubic response functions MXCRF=',MXCRF,'.',
     & ' The program will stop'
         CALL QUIT('Too many cubic response functions specified.')
      END IF
C                                                                                             
C     End of ABS_CRCHK                                                                        
C                                                                                             
      RETURN
      END
C                                                                                             
C                                                                                             
C
      SUBROUTINE ABS_NXY(CMO, UDV, PV, XINDX, MJWOP, FC, FCAC, FV,
     &                   H2AC, RESLRF, WRK, LWRK)
C                                                                                             
C ======================================================================                      
C                                                                                             
C  Purpose : Form the right-hand sides XY[1](omega_a,omega_b) required                        
C            for the CRF operator and frequency combinations set up by                        
C            ABS_GAMMA_SETUP                                                                  
C                                                                                             
C   Author : Tobias Fahleson, 2014                                                            
C                                                                                             
C ======================================================================                      
C                                                                                             
#include "implicit.h"
#include "inforb.h"
#include "priunit.h"
#include "absorp.h"
#include "abslrs.h"
#include "abscrs.h"
#include "inftap.h"
#include "rspprp.h"
#include "inflr.h"
#include "qrinf.h"
#include "infdim.h"
#include "infvar.h"
#include "wrkrsp.h"
C
      PARAMETER ( D0=0.0D0, D1=1.0D0, DM1=-1.0D0 )
      DIMENSION WRK(LWRK),CMO(*),UDV(*),PV(*),XINDX(*),MJWOP(*)
      DIMENSION FC(*), FCAC(*), FV(*), H2AC(*), RESLRF(*)
      CHARACTER*8 XLAB, YLAB, BLANK
      PARAMETER (BLANK='        ')
      LOGICAL FOUND, CONV, ISOLATE
      PARAMETER ( IDUMMY = - 9999999 )
C                                                                                             
      ISOLATE = .TRUE.
C                                                                                             
      KFREE = 1
      LFREE = LWRK
C                                                                                             
C     Main loop that loops over all ABS_NXYVEC operator and frequency                         
C     combinations set up by ABS_GAMMA_SETUP                                                  
C     ===============================================================                         
C
C ---
C PRINT FLAG
      WRITE(LUPRI,*) 'Starting loop over ABS_NXYVEC'
C ---                                                                                             
      DO INXY = 1,ABS_NXYVEC
C                                                                                             
C     Set the labels, frequencies, and symmetries of operators X and Y                        
      XLAB   = ABS_XYLAB(INXY,1)
      YLAB   = ABS_XYLAB(INXY,2)
      FREQX  = ABS_XYFREQ(INXY,1)
      FREQY  = ABS_XYFREQ(INXY,2)
      ISYMX  = ABS_XYSYM(INXY,1)
      ISYMY  = ABS_XYSYM(INXY,2)
      ISYMXY = MULD2H(ISYMX, ISYMY)
C                                                                                             
C     Set the length of the vectors Nx, Ny, and Nyx
      KZYVX  = MZYVAR(ISYMX)
      KZYVY  = MZYVAR(ISYMY)
      KZYVXY = MZYVAR(ISYMXY)
C                                                                                             
C     Only singlets                                                                           
      ISPINX  = 0
      ISPINY  = 0
      ISPINXY = 0
C                                                                                             
C     Allocate memory for current INXY vector and register the work-                          
C     array referens in ABS_XYMEMREF. In addition, allocate memory for                        
C     the three vectors T[3]NxNy, X[1]Ny, and Y[1]Nx, and also for Nx                         
C     and Ny                                                                                  
C     ================================================================                        
C                                                                                             
      ABS_XYMEMREF(INXY) = KFREE
      KVECT              = ABS_XYMEMREF(INXY) + 4*KZVAR
      KGRAD              = KVECT              + 4*KZVAR
      KVECX              = KGRAD              + 4*KZVAR
      KVECY              = KVECX              + 4*KZVAR
      KFREE              = KVECY              + 4*KZVAR
      LFREE              = LWRK - KFREE
C                                                                                             
C     Read in Nx                                                                              
C     ==========                                                                              
C                                                                                             
      CALL READ_XVEC2(LUABSVECS,2*KZYVX,WRK(KVECX),XLAB,BLANK,ISYMX,
     &               ABS(FREQX),D0,ABS_THCLR,FOUND,CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',XLAB,
     &           ' with frequency ',FREQX, ' and symmetry',
     &           ISYMX,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',XLAB,
     &           ' with frequency ',FREQX, ' and symmetry',
     &           ISYMX,' not converged on file LUABSVECS'
         END IF
      END IF
      IF (FREQX .LT. D0) THEN
         CALL DSWAP(KZYVX/2,WRK(KVECX),1,WRK(KVECX+KZYVX/2),1)
         CALL DSWAP(KZYVX/2,WRK(KVECX+KZYVX),1,
     &              WRK(KVECX+KZYVX+KZYVX/2),1)
         CALL DSCAL(KZYVX,DM1,WRK(KVECX),1)
      END IF
C                                                                                             
C     Read in Ny                                                                              
C     ==========                                                                              
C
      CALL READ_XVEC2(LUABSVECS,2*KZYVY,WRK(KVECY),YLAB,BLANK,ISYMY,
     &               ABS(FREQY), D0, ABS_THCLR, FOUND, CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F7.5,A,I3,/A)') ' Response label ',YLAB,
     &           ' with frequency ',FREQY, ' and symmetry',
     &           ISYMY,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F7.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',YLAB,
     &           ' with frequency ',FREQY, ' and symmetry',
     &           ISYMY,' not converged on file LUABSVECS'
         END IF
      END IF
      IF (FREQY .LT. D0) THEN
         CALL DSWAP(KZYVY/2,WRK(KVECY),1,WRK(KVECY+KZYVY/2),1)
         CALL DSWAP(KZYVY/2,WRK(KVECY+KZYVY),1,
     &              WRK(KVECY+KZYVY+KZYVY/2),1)
         CALL DSCAL(KZYVY,DM1,WRK(KVECY),1)
      END IF
C                                                                                             
C     Calculate T[3] Nx Ny with permutations done in T3DRV (R[3] term                         
C     absent here until further notice due to it being zero in                                
C     single-determinant processes (HF, DFT, ...))                                            
C     ===============================================================                         
C
      CALL DZERO(WRK(KVECT),4*KZVAR)
C                                                                                             
      CALL T3DRV(1,ISYMXY,ISYMX,ISYMY,WRK(KVECX),WRK(KVECY), .FALSE.,
     &           WRK(KVECX),-FREQX,-FREQY,XINDX,UDV,PV,MJWOP,
     &           WRK(KFREE),LFREE,CMO,FC,FV)
      CALL DAXPY(KZYVXY, D1, WRK(KFREE), 1, WRK(KVECT), 1)
C                                                                                             
      CALL T3DRV(1,ISYMXY,ISYMX,ISYMY,WRK(KVECX),WRK(KVECY+KZYVY),
     &           .FALSE.,WRK(KVECX),-FREQX,-FREQY,XINDX,UDV,PV,MJWOP,
     &           WRK(KFREE),LFREE,CMO,FC,FV)
      CALL DAXPY(KZYVXY, D1, WRK(KFREE), 1, WRK(KVECT+KZYVXY), 1)
C                                                                                             
      CALL T3DRV(1,ISYMXY,ISYMX,ISYMY,WRK(KVECX+KZYVX),WRK(KVECY),
     &           .FALSE.,WRK(KVECX),-FREQX,-FREQY,XINDX,UDV,PV,MJWOP,
     &           WRK(KFREE),LFREE,CMO,FC,FV)
      CALL DAXPY(KZYVXY, D1, WRK(KFREE), 1, WRK(KVECT+KZYVXY), 1)
C                                                                                             
      CALL T3DRV(1,ISYMXY,ISYMX,ISYMY,WRK(KVECX+KZYVX),WRK(KVECY+KZYVY),
     &           .FALSE.,WRK(KVECX),-FREQX,-FREQY,XINDX,UDV,PV,MJWOP,
     &           WRK(KFREE),LFREE,CMO,FC,FV)
      CALL DAXPY(KZYVXY, DM1, WRK(KFREE), 1, WRK(KVECT), 1)
C                                                                                             
C     Calculate X[2]Ny and Y[2]Nx                                                             
C     ===========================                                                             
C
      CALL DZERO(WRK(KGRAD),4*KZVAR)
      CALL X2INIT(1, KZYVXY, KZYVY, ISYMXY, ISPINXY, ISYMY, ISPINY,
     &            WRK(KVECX), WRK(KVECY), WRK(KGRAD), XINDX, UDV, PV,
     &            XLAB, ISYMX, ISPINX, CMO, MJWOP, WRK(KFREE), LFREE)
      CALL DAXPY(KZYVXY, DM1, WRK(KGRAD), 1, WRK(KVECT), 1)
C                                                                                             
      CALL DZERO(WRK(KGRAD),4*KZVAR)
      CALL X2INIT(1, KZYVXY, KZYVY, ISYMXY, ISPINXY, ISYMY, ISPINY,
     &           WRK(KVECX), WRK(KVECY+KZYVY), WRK(KGRAD), XINDX, UDV,
     &           PV, XLAB, ISYMX, ISPINX, CMO, MJWOP, WRK(KFREE), LFREE)
      CALL DAXPY(KZYVXY, DM1, WRK(KGRAD), 1, WRK(KVECT+KZYVXY), 1)
C                                                                                             
      CALL DZERO(WRK(KGRAD),4*KZVAR)
      CALL X2INIT(1, KZYVXY, KZYVX, ISYMXY, ISPINXY, ISYMX, ISPINX,
     &            WRK(KVECX), WRK(KVECX), WRK(KGRAD), XINDX, UDV, PV,
     &            YLAB, ISYMY, ISPINY, CMO, MJWOP, WRK(KFREE), LFREE)
      CALL DAXPY(KZYVXY, DM1, WRK(KGRAD), 1, WRK(KVECT),1)
C                                                                                             
      CALL DZERO(WRK(KGRAD),4*KZVAR)
      CALL X2INIT(1, KZYVXY, KZYVX, ISYMXY, ISPINXY, ISYMX, ISPINX,
     &           WRK(KVECX), WRK(KVECX+KZYVX), WRK(KGRAD), XINDX, UDV,
     &           PV, YLAB, ISYMY, ISPINY, CMO, MJWOP, WRK(KFREE), LFREE)
      CALL DAXPY(KZYVXY, DM1, WRK(KGRAD), 1, WRK(KVECT+KZYVXY),1)
C                                                                                             
C     Transfer the results from temporary memory WRK(KVECT) to a                              
C     permanent save at WRK(ABS_XYMEMREF(INXY))                                               
C     ==========================================================                              
C                                                                                             
      CALL DZERO(WRK(ABS_XYMEMREF(INXY)),4*KZVAR)
      CALL DAXPY(2*KZYVXY, D1, WRK(KVECT), 1, WRK(ABS_XYMEMREF(INXY)),1)
C                                                                                             
C     Nullify temporary memory, and rewind memory to just after                               
C     the latest NXY vector                                                                   
C     =========================================================                               
C
      CALL DZERO(WRK(KVECT),4*KZVAR)
      CALL DZERO(WRK(KGRAD),4*KZVAR)
      CALL DZERO(WRK(KVECX),4*KZVAR)
      CALL DZERO(WRK(KVECY),4*KZVAR)
      KFREE = KVECT
C                                                                                             
      END DO
C                                                                                             
C     Allocate memory for the NXY solution, which will be written                             
C     to file                                                                                 
      KXSOL = KFREE
      KFREE = KXSOL + 4*KZVAR
      LFREE = LWRK  - KFREE
C                                                                                             
C     Logical in ABS_CTL that enables the double-indexed solver mode                          
      ABS_IDRI_NXY = .TRUE.
C                                                                                             
C     Putting NFREQS and NGD to 1 so that ABS_CTL is called once for                          
C     every ABS_NXYVEC                                                                        
      NFREQS = 1
      NGD = 1
C                                                                                             
C     Call the LR solver for every ABS_NXYVEC                                                 
C     =======================================                                                 
C                                                                                             
      DO INXY = 1,ABS_NXYVEC
C                                                                                             
C        Open temporary storage file for ABS_CTL                                              
         LUTEMP   = -1
         CALL GPOPEN(LUTEMP,'ABS_NXY_TEMP','NEW',' ',' ',
     &        IDUMMY,.FALSE.)
         WRITE(LUTEMP) 'EOFLABEL'
C                                                                                             
C        Set labels, frequencies, and symmetries
         XLAB   = ABS_XYLAB(INXY,1)
         YLAB   = ABS_XYLAB(INXY,2)
         FREQX  = ABS_XYFREQ(INXY,1)
         FREQY  = ABS_XYFREQ(INXY,2)
         ISYMX  = ABS_XYSYM(INXY,1)
         ISYMY  = ABS_XYSYM(INXY,2)
         ISYMXY = MULD2H(ISYMX, ISYMY)
C                                                                                             
C        Form the sum of frequencies for term (omega1 + omega2)S[2]                           
         FREQS = FREQX + FREQY
C                                                                                             
C        Determine the norm of the r.h.s.                                                     
         VNORM = DNRM2(4*KZVAR,WRK(ABS_XYMEMREF(INXY)),1)
C                                                                                             
         CALL DZERO(WRK(KXSOL),4*KZVAR)
C                                                                                             
C        Check if norm of vector is sufficiently big for solver                               
         IF (VNORM .GT. 1.0D-6) THEN
C                                                                                             
C        Open temporary auxiliary files                                                       
         LUSB = -1
         LUAB = -1
         LUSS = -1
         LUAS = -1
         CALL GPOPEN(LUSB,'ABS_SB','NEW',' ',' ',IDUMMY,.FALSE.)
         CALL GPOPEN(LUAB,'ABS_AB','NEW',' ',' ',IDUMMY,.FALSE.)
         CALL GPOPEN(LUSS,'ABS_SS','NEW',' ',' ',IDUMMY,.FALSE.)
         CALL GPOPEN(LUAS,'ABS_AS','NEW',' ',' ',IDUMMY,.FALSE.)
C                                                                                             
         CALL ABS_CTL('XYLAB   ',KZVAR,WRK(ABS_XYMEMREF(INXY)),
     &        NGD, WRK(KXSOL), FREQS, NFREQS, LUTEMP,
     &        MJWOP, CMO, UDV, FC, FCAC, FV, PV, XINDX,
     &        RESLRF, WRK(KFREE),LFREE)
C
         CALL GPCLOSE(LUSB,'DELETE')
         CALL GPCLOSE(LUAB,'DELETE')
         CALL GPCLOSE(LUSS,'DELETE')
         CALL GPCLOSE(LUAS,'DELETE')
C                                                                                             
C        Transfer the solution vector from a temporary to a                                   
C        permanent storage                                                                    
C        ==================================================                                   
C                                                                                             
C        Read in the vector
         CALL DZERO(WRK(KXSOL),4*KZVAR)
         CALL READ_XVEC2(LUTEMP, 4*KZVAR, WRK(KXSOL), 'XYLAB   ', BLANK,
     &                   ISYMXY, FREQS, D0, ABS_THCLR, FOUND, CONV)
         IF (.NOT. (FOUND .AND. CONV)) THEN
            IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/4A,2(F7.5),A,2(I3),/A)')
     &           ' Response label ','XYLAB',BLANK,' with frequencies ',
     &           FREQX,FREQY,' and symmetries',ISYMX,ISYMY,
     &           ' not found on file LUTEMP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/4A,2(F7.5),/A,2(I3),A)') ' @WARNING>>>>'//
     &           ' Response label ','XYLAB',BLANK,
     &           ' with frequencies ',FREQX,FREQY,
     &           ' and symmetries', ISYMX, ISYMY,
     &           ' not converged on file LUABSVECS'
            END IF
         END IF
         IF (FREQS .LT. D0) THEN
            CALL DSWAP(KZVAR,WRK(KXSOL),1,WRK(KXSOL+KZVAR),1)
            CALL DSWAP(KZVAR,WRK(KXSOL+2*KZVAR),1,
     &           WRK(KXSOL+2*KZVAR+KZVAR),1)
            CALL DSCAL(2*KZVAR,DM1,WRK(KXSOL),1)
         END IF
C                                                                                             
         END IF
C                                                                                             
C        Write the vector to file                                                             
         CALL WRITE_XVEC2(LUABSVECS,4*KZVAR,WRK(KXSOL),XLAB,YLAB,
     &                    FREQX,FREQY,ABS_THCLR)
C                                                                                             
         CALL GPCLOSE(LUTEMP, 'DELETE')
C                                                                                             
      END DO
C                                                                                             
C     End of ABS_NXY                                                                          
C                                                                                             
      RETURN
      END
C                                                                                             
C                                                                                             
C
