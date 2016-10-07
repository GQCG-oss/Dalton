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
C     abscrs.F: Main author Tobias Fahleson
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
C                  Check if it already has been registered as a doublet
C                  (this unintentional double counting may occur for 
C                  e.g. highly symmetric systems)
                   DO INEQTO=1,ABS_CRF_NEQTO
                   IF ( ALAB  .EQ. ABS_CRFLAB_EQTO(INEQTO,1) .AND.
     &                  BLAB  .EQ. ABS_CRFLAB_EQTO(INEQTO,2) .AND.
     &                  CLAB  .EQ. ABS_CRFLAB_EQTO(INEQTO,3) .AND.
     &                  DLAB  .EQ. ABS_CRFLAB_EQTO(INEQTO,4) .AND.
     &                  ISYMA .EQ. ABS_CRFSYM_EQTO(INEQTO,1) .AND.
     &                  ISYMB .EQ. ABS_CRFSYM_EQTO(INEQTO,2) .AND.
     &                  ISYMC .EQ. ABS_CRFSYM_EQTO(INEQTO,3) .AND.
     &                  ISYMD .EQ. ABS_CRFSYM_EQTO(INEQTO,4) .AND.
     &                  FREQA .EQ. ABS_CRFFRQ_EQTO(INEQTO,1) .AND.
     &                  FREQB .EQ. ABS_CRFFRQ_EQTO(INEQTO,2) .AND.
     &                  FREQC .EQ. ABS_CRFFRQ_EQTO(INEQTO,3) .AND.
     &                  FREQD .EQ. ABS_CRFFRQ_EQTO(INEQTO,4) ) THEN
                      RETURN
                   END IF
                   END DO
C
C                  Check main list if it identically appears there
                   DO INCRF=1,ABS_NCRF
                   IF ( ALAB  .EQ. ABS_CRFLAB(INCRF,1) .AND.
     &                  BLAB  .EQ. ABS_CRFLAB(INCRF,2) .AND.
     &                  CLAB  .EQ. ABS_CRFLAB(INCRF,3) .AND.
     &                  DLAB  .EQ. ABS_CRFLAB(INCRF,4) .AND.
     &                  ISYMA .EQ. ABS_CRFSYM(INCRF,1) .AND.
     &                  ISYMB .EQ. ABS_CRFSYM(INCRF,2) .AND.
     &                  ISYMC .EQ. ABS_CRFSYM(INCRF,3) .AND.
     &                  ISYMD .EQ. ABS_CRFSYM(INCRF,4) .AND.
     &                  FREQA .EQ. ABS_CRFFRQ(INCRF,1) .AND.
     &                  FREQB .EQ. ABS_CRFFRQ(INCRF,2) .AND.
     &                  FREQC .EQ. ABS_CRFFRQ(INCRF,3) .AND.
     &                  FREQD .EQ. ABS_CRFFRQ(INCRF,4) ) THEN
                      RETURN
                   END IF
                   END DO
C
C                  If not, register current element's information and
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
     &                   H2AC, FOCK, RESLRF, WRK, LWRK)
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
      DIMENSION FC(*), FCAC(*), FV(*), H2AC(*), FOCK(*), RESLRF(*)
      CHARACTER*8 XLAB, YLAB, BLANK
      PARAMETER (BLANK='        ')
      LOGICAL FOUND, CONV, ISOLATE
      PARAMETER ( IDUMMY = - 9999999 )
C                                                                                             
      ISOLATE = .FALSE.
C                                                                                             
      KFREE = 1
      LFREE = LWRK
C
C     Determine MZYVAR for all symmetries NSYM
C     ========================================
C
      DO ISYM=1,NSYM
         KSYMOP = ISYM
         CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KFREE),LFREE)
         CALL SETZY(MJWOP)
      END DO
C                                                                                             
C     Main loop that loops over all ABS_NXYVEC operator and frequency                         
C     combinations set up by ABS_GAMMA_SETUP                                                  
C     ===============================================================                         
C
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
      KVECT              = ABS_XYMEMREF(INXY) + 2*KZYVXY
      KGRAD              = KVECT              + 2*KZYVXY
      KVECX              = KGRAD              + 2*KZYVXY
      KVECY              = KVECX              + 2*KZYVX
      KFREE              = KVECY              + 2*KZYVY
      LFREE              = LWRK - KFREE
C
C    If by symmetry the vectors Nx or Ny are zero, set Nxy to zero and
C    skip the calculation
C    =================================================================
C
      IF ((KZYVX .EQ. 0) .OR. (KZYVY .EQ. 0)) THEN
         CALL DZERO(WRK(ABS_XYMEMREF(INXY)),2*KZYVXY)
      ELSE
C                                                                                             
C     Read in Nx                                                                              
C     ==========                                                                              
C                                                                                             
      CALL READ_XVEC2(LUABSVECS,2*KZYVX,WRK(KVECX),XLAB,BLANK,ISYMX,
     &               ABS(FREQX),D0,ABS_THCLR,FOUND,CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F8.5,A,I3,/A)') ' Response label ',XLAB,
     &           ' with frequency ',FREQX, ' and symmetry',
     &           ISYMX,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F8.5,/A,I3,A)') ' @WARNING>>>>'//
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
            WRITE (LUPRI,'(/3A,F8.5,A,I3,/A)') ' Response label ',YLAB,
     &           ' with frequency ',FREQY, ' and symmetry',
     &           ISYMY,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F8.5,/A,I3,A)') ' @WARNING>>>>'//
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
      CALL DZERO(WRK(KVECT),2*KZYVXY)
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
      CALL DZERO(WRK(KGRAD),2*KZYVXY)
      CALL X2INIT(1, KZYVXY, KZYVY, ISYMXY, ISPINXY, ISYMY, ISPINY,
     &            WRK(KVECX), WRK(KVECY), WRK(KGRAD), XINDX, UDV, PV,
     &            XLAB, ISYMX, ISPINX, CMO, MJWOP, WRK(KFREE), LFREE)
      CALL DAXPY(KZYVXY, DM1, WRK(KGRAD), 1, WRK(KVECT), 1)
C                                                                                             
      CALL DZERO(WRK(KGRAD),2*KZYVXY)
      CALL X2INIT(1, KZYVXY, KZYVY, ISYMXY, ISPINXY, ISYMY, ISPINY,
     &           WRK(KVECX), WRK(KVECY+KZYVY), WRK(KGRAD), XINDX, UDV,
     &           PV, XLAB, ISYMX, ISPINX, CMO, MJWOP, WRK(KFREE), LFREE)
      CALL DAXPY(KZYVXY, DM1, WRK(KGRAD), 1, WRK(KVECT+KZYVXY), 1)
C                                                                                             
      CALL DZERO(WRK(KGRAD),2*KZYVXY)
      CALL X2INIT(1, KZYVXY, KZYVX, ISYMXY, ISPINXY, ISYMX, ISPINX,
     &            WRK(KVECX), WRK(KVECX), WRK(KGRAD), XINDX, UDV, PV,
     &            YLAB, ISYMY, ISPINY, CMO, MJWOP, WRK(KFREE), LFREE)
      CALL DAXPY(KZYVXY, DM1, WRK(KGRAD), 1, WRK(KVECT),1)
C                                                                                             
      CALL DZERO(WRK(KGRAD),2*KZYVXY)
      CALL X2INIT(1, KZYVXY, KZYVX, ISYMXY, ISPINXY, ISYMX, ISPINX,
     &           WRK(KVECX), WRK(KVECX+KZYVX), WRK(KGRAD), XINDX, UDV,
     &           PV, YLAB, ISYMY, ISPINY, CMO, MJWOP, WRK(KFREE), LFREE)
      CALL DAXPY(KZYVXY, DM1, WRK(KGRAD), 1, WRK(KVECT+KZYVXY),1)
C                                                                                             
C     Transfer the results from temporary memory WRK(KVECT) to a                              
C     permanent save at WRK(ABS_XYMEMREF(INXY))                                               
C     ==========================================================                              
C                                                                                             
      CALL DZERO(WRK(ABS_XYMEMREF(INXY)),2*KZYVXY)
      CALL DAXPY(2*KZYVXY, D1, WRK(KVECT), 1, WRK(ABS_XYMEMREF(INXY)),1)
C                                                                                             
C     Nullify temporary memory, and rewind memory to just after                               
C     the latest NXY vector                                                                   
C     =========================================================                               
C
      CALL DZERO(WRK(KVECT),2*KZYVXY)
      CALL DZERO(WRK(KGRAD),2*KZYVXY)
      CALL DZERO(WRK(KVECX),2*KZYVX)
      CALL DZERO(WRK(KVECY),2*KZYVY)
      KFREE = KVECT
C
      END IF
C                                                                                             
      END DO
C                                                                                             
C     Allocate memory for the NXY solution, which will be written                             
C     to file (this was moved inside the loop)
C      KXSOL = KFREE
C      KFREE = KXSOL + 4*KZVAR
C      LFREE = LWRK  - KFREE
C                                                                                             
C     Logical in ABS_CTL that enables the double-indexed solver mode
C     (this was originally introduced in accordance with multiple r.h.s.
C      but that is put on ice)
      ABS_IDRI_XY = .TRUE.
C                                                                                             
C     Putting NFREQS and NGD to 1 so that ABS_CTL is called once for                          
C     every ABS_NXYVEC                                                                        
      NFREQS = 1
C
C     Setting NGD to 1 would mean only 1 r.h.s., but multiple r.h.s. is
C     put on ice
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
C        Determine various response variables for symmetry ISYMXY
         KSYMOP = ISYMXY
         abs_gradsym = isymxy
         CALL RSPVAR(UDV,FOCK,FC,FV,FCAC,H2AC,XINDX,WRK(KFREE),LFREE)
         CALL SETZY(MJWOP)
C
C        Set the length of vector Nxy
         KZYVXY = MZYVAR(ISYMXY)
C                                                                                             
C        Form the sum of frequencies for term (omega1 + omega2)S[2]                           
         FREQS = FREQX + FREQY
C
C        Allocate memory for the NXY solution, which will be written
C        to file
         KXSOL = KFREE
         KFREE = KXSOL + 2*KZYVXY
         LFREE = LWRK  - KFREE
C                                                                                             
C        Determine the norm of the r.h.s.                                                     
         VNORM = DNRM2(2*KZYVXY,WRK(ABS_XYMEMREF(INXY)),1)
C                                                                                             
         CALL DZERO(WRK(KXSOL),2*KZYVXY)
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
C        Solve the linear response equation
         CALL ABS_CTL('XYLAB   ', KZYVXY/2, WRK(ABS_XYMEMREF(INXY)),
     &        NGD, WRK(KXSOL), FREQS, NFREQS, LUTEMP,
     &        MJWOP, CMO, UDV, FC, FCAC, FV, PV, XINDX,
     &        RESLRF, WRK(KFREE), LFREE)
C
C        Delete temporary files
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
         CALL DZERO(WRK(KXSOL),2*KZYVXY)
         CALL READ_XVEC2(LUTEMP, 2*KZYVXY, WRK(KXSOL), 'XYLAB   ',BLANK,
     &                   ISYMXY, FREQS, D0, ABS_THCLR, FOUND, CONV)
         IF (.NOT. (FOUND .AND. CONV)) THEN
            IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/4A,2(F8.5),A,2(I3),/A)')
     &           ' Response label ','XYLAB',BLANK,' with frequencies ',
     &           FREQX,FREQY,' and symmetries',ISYMX,ISYMY,
     &           ' not found on file LUTEMP'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/4A,2(F8.5),/A,2(I3),A)') ' @WARNING>>>>'//
     &           ' Response label ','XYLAB',BLANK,
     &           ' with frequencies ',FREQX,FREQY,
     &           ' and symmetries', ISYMX, ISYMY,
     &           ' not converged on file LUABSVECS'
            END IF
         END IF
         IF (FREQS .LT. D0) THEN
            CALL DSWAP(KZYVXY/2,WRK(KXSOL),1,WRK(KXSOL+KZYVXY/2),1)
            CALL DSWAP(KZYVXY/2,WRK(KXSOL+KZYVXY),1,
     &           WRK(KXSOL+KZYVXY+KZYVXY/2),1)
            CALL DSCAL(KZYVXY,DM1,WRK(KXSOL),1)
         END IF
C                                                                                             
         END IF
C                                                                                             
C        Write the vector to file                                                             
         CALL WRITE_XVEC2(LUABSVECS,2*KZYVXY,WRK(KXSOL),XLAB,YLAB,
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
      SUBROUTINE ABS_CRF(CMO, UDV, PV, FOCK, FC, FV, FCAC, H2AC, XINDX,  
     &                   MJWOP, WRK, LWRK)
C
C ======================================================================
C
C  Purpose : Compute the cubic response functions and print the results
C
C   Author : Tobias Fahleson, 2014
C
C ======================================================================
C
#include "implicit.h"
#include "priunit.h"
#include "absorp.h"
#include "wrkrsp.h"
#include "infvar.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "symmet.h"
#include "cbilrs.h"
#include "infrsp.h"
#include "inforb.h"
#include "qrinf.h"
#include "infdim.h"
#include "abslrs.h"
#include "abscrs.h"
#include "thrzer.h"
C
      PARAMETER ( D0=0.0D0, DM1=-1.0D0 )
C
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*),FV(*),FCAC(*),H2AC(*)
      DIMENSION XINDX(*), MJWOP(*), WRK(LWRK)
C
      CHARACTER*8 ALAB, BLAB, CLAB, DLAB
      DIMENSION HYPVAL(2), TMPVAL(2), VAL(4), IBCDCOMPLEX(2,2,2)
      LOGICAL ISOLATE
      ISOLATE = .TRUE.
C
C     Spin multiplicity is always 1, i.e. only singlets 
      ISPINA  = 0
      ISPINB  = 0
      ISPINC  = 0
      ISPIND  = 0
      ISPINBC = 0
      ISPINBD = 0
      ISPINCD = 0
C
C     Allocating appropriate memory for response vectors
      MZYVMX = 2*NVARMA
      KVECA  = 1
      KVECB  = KVECA  + 2*MZYVMX
      KVECC  = KVECB  + 2*MZYVMX
      KVECD  = KVECC  + 2*MZYVMX
      KVECBC = KVECD  + 2*MZYVMX
      KVECBD = KVECBC + 2*MZYVMX
      KVECCD = KVECBD + 2*MZYVMX
      KGRAD  = KVECCD + 2*MZYVMX
      KFREE  = KGRAD  + 2*MZYVMX
      LFREE  = LWRK   - KFREE
C
C     Opening the file where response vectors are saved
      LURSPRES = -1
      CALL GPOPEN(LURSPRES,'RESULTS.RSP','UNKNOWN','SEQUENTIAL',
     &     'FORMATTED',0,.FALSE.)
C
C     Main loop that loops over all the CRF set up by ABS_GAMMA_SETUP
C     ===============================================================
C
      DO ICRF = 1,ABS_NCRF
C
C     Real and imaginary values of the hyperpolarizability
      HYPVAL(1) = D0
      HYPVAL(2) = D0
C
C     Labels, symmetries, and frequencies, of operators A to D
      ALAB   = ABS_CRFLAB(ICRF,1)
      BLAB   = ABS_CRFLAB(ICRF,2)
      CLAB   = ABS_CRFLAB(ICRF,3)
      DLAB   = ABS_CRFLAB(ICRF,4)
      ISYMA  = ABS_CRFSYM(ICRF,1)
      ISYMB  = ABS_CRFSYM(ICRF,2)
      ISYMC  = ABS_CRFSYM(ICRF,3)
      ISYMD  = ABS_CRFSYM(ICRF,4)
      ISYMBC = MULD2H(ISYMB,ISYMC)
      ISYMBD = MULD2H(ISYMB,ISYMD)
      ISYMCD = MULD2H(ISYMC,ISYMD)
      FREQA  = ABS_CRFFRQ(ICRF,1)
      FREQB  = ABS_CRFFRQ(ICRF,2)
      FREQC  = ABS_CRFFRQ(ICRF,3)
      FREQD  = ABS_CRFFRQ(ICRF,4)
C
C     Write details of current CRF to output
C     ======================================
C
      WRITE(LUPRI,'(/,2(A,I4),A,4(/A,A10,I4,F10.6))')
     &        ' Cubic response function no',ICRF,' of',ABS_NCRF,'.',
     &        ' A operator, symmetry, frequency: ', ALAB, ISYMA, FREQA,
     &        ' B operator, symmetry, frequency: ', BLAB, ISYMB, FREQB,
     &        ' C operator, symmetry, frequency: ', CLAB, ISYMC, FREQC,
     &        ' D operator, symmetry, frequency: ', DLAB, ISYMD, FREQD
C     
C     Read in requested response vectors
C     ==================================
C
      CALL ABS_CRRDVE(ISYMA,       ISYMB,       ISYMC,       ISYMD,
     &                ISYMBC,      ISYMBD,      ISYMCD,
     &                ALAB,        BLAB,        CLAB,        DLAB,
     &                FREQA,       FREQB,       FREQC,       FREQD,
     &                KZYVA,       KZYVB,       KZYVC,       KZYVD,
     &                KZYVBC,      KZYVBD,      KZYVCD,
     &                WRK(KVECA),  WRK(KVECB),  WRK(KVECC),  WRK(KVECD),
     &                WRK(KVECBC), WRK(KVECBD), WRK(KVECCD))
C
C     Check if some of the single-indexed vectors are equal or zero
C     =============================================================
C
      CALL ABS_ABCDCHK_INIT(IAEQ,  IBCDCOMPLEX,
     &     FREQA,      FREQB,      FREQC,      FREQD,
     &     ISYMA,      ISYMB,      ISYMC,      ISYMD,
     &     KZYVA,      KZYVB,      KZYVC,      KZYVD,
     &     WRK(KVECA), WRK(KVECB), WRK(KVECC), WRK(KVECD))
C
C     Print individual contributions if requested
      IF (IPRABSLRS .GT. 0) THEN
      WRITE(LUPRI,'(//A13,2A36,/A13,2A36,/A13,4A18,/A)')
     &   'Contribution','Term','Accumulated',
     &   ' ','----------------------','----------------------',
     &   ' ','Real','Imag','Real','Imag',
     &   ' ---------------------------------------------------------'//
     &   '---------------------------'
      END IF
C
C     1. Calculate Na T[3] Nb Ncd type terms (three permutations)
C     ===========================================================
C
      TMPVAL(1) = D0
      TMPVAL(2) = D0
C
C     If WRK(KVECA) is zero, add zero T[3] contribution
      IF (IAEQ .NE. 0) THEN
C
C     This is the A B CD permutation
      CALL T3COMPLEX(KZYVA,      KZYVB,      KZYVCD,
     &               ISYMA,      ISYMB,      ISYMCD,
     &               WRK(KVECA), WRK(KVECB), WRK(KVECCD),
     &               VAL(1),     FREQB,      FREQC+FREQD,
     &               ABS_DAMP, XINDX, UDV, PV, MJWOP, CMO, FC, FV,
     &               WRK(KFREE), LFREE)
      VAL(1) = -VAL(1)
      VAL(2) = -VAL(2)
      TMPVAL(1) = VAL(1)
      TMPVAL(2) = VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
C     If the real and imaginary parts of all single-indexed vectors are
C     equal, add twice the value of the first permutation
      IF (IBCDCOMPLEX(1,1,1).EQ.24 .AND. IBCDCOMPLEX(2,2,2).EQ.24) THEN
         TMPVAL(1) = TMPVAL(1) + 2*VAL(1)
         TMPVAL(2) = TMPVAL(2) + 2*VAL(2)
         HYPVAL(1) = HYPVAL(1) + 2*VAL(1)
         HYPVAL(2) = HYPVAL(2) + 2*VAL(2)
C
C     (The following cases have been commented out due to incomp-
C      atibility with IBCDCOMPLEX as of 27/5 2015) 
C
C     If vectors WRK(KVECB) and WRK(KVECC) are equal, add the value of
C     the first permutation and calculate the A D BC permutation
C      ELSE IF (IBCDEQ.EQ.2) THEN
C         TMPVAL(1) = TMPVAL(1) + VAL(1)
C         TMPVAL(2) = TMPVAL(2) + VAL(2)
C         HYPVAL(1) = HYPVAL(1) + VAL(1)
C         HYPVAL(2) = HYPVAL(2) + VAL(2)
C         CALL T3COMPLEX(KZYVA,      KZYVD,      KZYVBC,
C     &                  ISYMA,      ISYMD,      ISYMBC,
C     &                  WRK(KVECA), WRK(KVECD), WRK(KVECBC),
C     &                  VAL(1),     FREQD,      FREQB+FREQC,
C     &                  ABS_DAMP, XINDX, UDV, PV, MJWOP, CMO, FC, FV,
C     &                  WRK(KFREE), LFREE)
C         VAL(1) = -VAL(1)
C         VAL(2) = -VAL(2)
C         TMPVAL(1) = TMPVAL(1) + VAL(1)
C         TMPVAL(2) = TMPVAL(2) + VAL(2)
C         HYPVAL(1) = HYPVAL(1) + VAL(1)
C         HYPVAL(2) = HYPVAL(2) + VAL(2)
C
C     If vectors WRK(KVECB) and WRK(KVECD) are equal, add the value of
C     the first permutation and calculate the A C BD permutation
C      ELSE IF (IBCDEQ.EQ.3) THEN
C         TMPVAL(1) = TMPVAL(1) + VAL(1)
C         TMPVAL(2) = TMPVAL(2) + VAL(2)
C         HYPVAL(1) = HYPVAL(1) + VAL(1)
C         HYPVAL(2) = HYPVAL(2) + VAL(2)
C         CALL T3COMPLEX(KZYVA,      KZYVC,      KZYVBD,
C     &                  ISYMA,      ISYMC,      ISYMBD,
C     &                  WRK(KVECA), WRK(KVECC), WRK(KVECBD),
C     &                  VAL(1),     FREQC,      FREQB+FREQD,
C     &                  ABS_DAMP, XINDX, UDV, PV, MJWOP, CMO, FC, FV,
C     &                  WRK(KFREE), LFREE)
C         VAL(1) = -VAL(1)
C         VAL(2) = -VAL(2)
C         TMPVAL(1) = TMPVAL(1) + VAL(1)
C         TMPVAL(2) = TMPVAL(2) + VAL(2)
C         HYPVAL(1) = HYPVAL(1) + VAL(1)
C         HYPVAL(2) = HYPVAL(2) + VAL(2)
C
C     If vectors WRK(KVECC) and WRK(KVECD) are equal, calculate the
C     A C BD permutation and add twice its value
C      ELSE IF (IBCDEQ.EQ.4) THEN
C         CALL T3COMPLEX(KZYVA,      KZYVC,      KZYVBD,
C     &                  ISYMA,      ISYMC,      ISYMBD,
C     &                  WRK(KVECA), WRK(KVECC), WRK(KVECBD),
C     &                  VAL(1),     FREQC,      FREQB+FREQD,
C     &                  ABS_DAMP, XINDX, UDV, PV, MJWOP, CMO, FC, FV,
C     &                  WRK(KFREE), LFREE)
C         VAL(1) = -VAL(1)
C         VAL(2) = -VAL(2)
C         TMPVAL(1) = TMPVAL(1) + 2*VAL(1)
C         TMPVAL(2) = TMPVAL(2) + 2*VAL(2)
C         HYPVAL(1) = HYPVAL(1) + 2*VAL(1)
C         HYPVAL(2) = HYPVAL(2) + 2*VAL(2)
C
C     If no single-indexed vectors are equal, add the remaining two
C     permutations A C BD and A D BC
      ELSE
         CALL T3COMPLEX(KZYVA,      KZYVC,      KZYVBD,
     &                  ISYMA,      ISYMC,      ISYMBD,
     &                  WRK(KVECA), WRK(KVECC), WRK(KVECBD),
     &                  VAL(1),     FREQC,      FREQB+FREQD,
     &                  ABS_DAMP, XINDX, UDV, PV, MJWOP, CMO, FC, FV,
     &                  WRK(KFREE), LFREE)
         VAL(1) = -VAL(1)
         VAL(2) = -VAL(2)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
C
         CALL T3COMPLEX(KZYVA,      KZYVD,      KZYVBC,
     &                  ISYMA,      ISYMD,      ISYMBC,
     &                  WRK(KVECA), WRK(KVECD), WRK(KVECBC),
     &                  VAL(1),     FREQD,      FREQB+FREQC,
     &                  ABS_DAMP, XINDX, UDV, PV, MJWOP, CMO, FC, FV,
     &                  WRK(KFREE), LFREE)
         VAL(1) = -VAL(1)
         VAL(2) = -VAL(2)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
      END IF
C
      IF (IPRABSLRS .GT. 0) THEN
      WRITE(LUPRI,'(A13,4F18.8)') 'Na T[3] NbNcd',
     &        TMPVAL(1),TMPVAL(2),HYPVAL(1),HYPVAL(2)
      END IF
C
C     2. Calculate Na T[4] Nb Nc Nd type terms (all permutations
C     done in a single call)
C     ==========================================================
C
      CALL T4COMPLEX(IAEQ,       IBCDCOMPLEX,
     &               ISYMA,      ISYMB,      ISYMC,      ISYMD,
     &               KZYVA,      KZYVB,      KZYVC,      KZYVD,
     &               WRK(KVECA), WRK(KVECB), WRK(KVECC), WRK(KVECD),
     &                           FREQB,      FREQC,      FREQD,
     &               ABS_DAMP, XINDX, UDV, PV, MJWOP, CMO, FC,
     &               WRK(KFREE), LFREE, VAL(1))
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL R4DRV(IAEQ,       IBCDCOMPLEX,
     &              ISYMA,      ISYMB,       ISYMC,      ISYMD,
     &              KZYVA,      KZYVB,       KZYVC,      KZYVD,
     &              WRK(KVECA), WRK(KVECB),  WRK(KVECC), WRK(KVECD),
     &              ABS_DAMP, XINDX, UDV, MJWOP,
     &              WRK(KFREE), LFREE, VAL(3))
      END IF
C
      TMPVAL(1) = VAL(1) + VAL(3)
      TMPVAL(2) = VAL(2) + VAL(4)
      HYPVAL(1) = HYPVAL(1) + VAL(1) + VAL(3)
      HYPVAL(2) = HYPVAL(2) + VAL(2) + VAL(4)
C
      IF (IPRABSLRS .GT. 0) THEN
      WRITE(LUPRI,'(A13,4F18.8)') 'Na T[4] NbNcNd',
     &        TMPVAL(1),TMPVAL(2),HYPVAL(1),HYPVAL(2)
      END IF
C
C     3. Calculate Na B[2] Ncd type terms (three permutations)
C     ========================================================
C
C     This is the A B CD permutation
      CALL DZERO(WRK(KGRAD),2*MZYVMX)
      CALL X2INIT(1, KZYVA, KZYVCD, ISYMA, ISPINA, ISYMCD, ISPINCD,
     &            WRK(KVECA), WRK(KVECCD), WRK(KGRAD), XINDX, UDV, PV,
     &            BLAB, ISYMB, ISPINB, CMO, MJWOP, WRK(KFREE), LFREE)
      VAL(1) = DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA), 1)
      VAL(2) = DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA + KZYVA), 1)
      TMPVAL(1) = VAL(1)
      TMPVAL(2) = VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL X2INIT(1, KZYVA, KZYVCD, ISYMA, ISPINA, ISYMCD, ISPINCD,
     &            WRK(KVECA), WRK(KVECCD + KZYVCD), WRK(KGRAD), XINDX,
     &            UDV, PV, BLAB, ISYMB, ISPINB, CMO, MJWOP,
     &            WRK(KFREE), LFREE)
         VAL(1) = -DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA + KZYVA), 1)
         VAL(2) =  DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA), 1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
C     This is the A C BD permutation
      CALL DZERO(WRK(KGRAD),2*MZYVMX)
      CALL X2INIT(1, KZYVA, KZYVBD, ISYMA, ISPINA, ISYMBD, ISPINBD,
     &            WRK(KVECA), WRK(KVECBD), WRK(KGRAD), XINDX, UDV, PV,
     &            CLAB, ISYMC, ISPINC, CMO, MJWOP, WRK(KFREE), LFREE)
      VAL(1) = DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA), 1)
      VAL(2) = DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA + KZYVA), 1)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL X2INIT(1, KZYVA, KZYVBD, ISYMA, ISPINA, ISYMBD, ISPINBD,
     &            WRK(KVECA), WRK(KVECBD + KZYVBD), WRK(KGRAD), XINDX,
     &            UDV, PV, CLAB, ISYMC, ISPINC, CMO, MJWOP,
     &            WRK(KFREE), LFREE)
         VAL(1) = -DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA + KZYVA), 1)
         VAL(2) =  DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA), 1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
C     This is the A D BC permutation
      CALL DZERO(WRK(KGRAD),2*MZYVMX)
      CALL X2INIT(1, KZYVA, KZYVBC, ISYMA, ISPINA, ISYMBC, ISPINBC,
     &            WRK(KVECA), WRK(KVECBC), WRK(KGRAD), XINDX, UDV, PV,
     &            DLAB, ISYMD, ISPIND, CMO, MJWOP, WRK(KFREE), LFREE)
      VAL(1) = DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA), 1)
      VAL(2) = DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA + KZYVA), 1)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL X2INIT(1, KZYVA, KZYVBC, ISYMA, ISPINA, ISYMBC, ISPINBC,
     &            WRK(KVECA), WRK(KVECBC + KZYVBC), WRK(KGRAD), XINDX,
     &            UDV, PV, DLAB, ISYMD, ISPIND, CMO, MJWOP,
     &            WRK(KFREE), LFREE)
         VAL(1) = -DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA + KZYVA), 1)
         VAL(2) =  DDOT(KZYVA, WRK(KGRAD), 1, WRK(KVECA), 1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
      IF (IPRABSLRS .GT. 0) THEN
      WRITE(LUPRI,'(A13,4F18.8)') 'Na B[2] Ncd',
     &        TMPVAL(1),TMPVAL(2),HYPVAL(1),HYPVAL(2)
      END IF
C
C     4. Calculate Na B[3] Nc Nd type terms (two of six permutations in
C     each call)
C     =================================================================
C
C     This is the A B C D permutation
      CALL X3COMPLEX(KZYVA,      KZYVC,      KZYVD,
     &               ISYMA,      ISYMC,      ISYMD,
     &               WRK(KVECA), WRK(KVECC), WRK(KVECD),
     &               VAL(1), BLAB, ISYMB, XINDX, UDV, CMO, MJWOP,
     &               WRK(KFREE), LFREE)
      VAL(1) = -VAL(1)
      VAL(2) = -VAL(2)
      TMPVAL(1) = VAL(1)
      TMPVAL(2) = VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
C     This is the A C B D permutation
      CALL X3COMPLEX(KZYVA,      KZYVB,      KZYVD,
     &               ISYMA,      ISYMB,      ISYMD,
     &               WRK(KVECA), WRK(KVECB), WRK(KVECD),
     &               VAL(1), CLAB, ISYMC, XINDX, UDV, CMO, MJWOP,
     &               WRK(KFREE), LFREE)
      VAL(1) = -VAL(1)
      VAL(2) = -VAL(2)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
C     This is the A D B C permutation
      CALL X3COMPLEX(KZYVA,      KZYVB,      KZYVC,
     &               ISYMA,      ISYMB,      ISYMC,
     &               WRK(KVECA), WRK(KVECB), WRK(KVECC),
     &               VAL(1), DLAB, ISYMD, XINDX, UDV, CMO, MJWOP,
     &               WRK(KFREE), LFREE)
      VAL(1) = -VAL(1)
      VAL(2) = -VAL(2)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (IPRABSLRS .GT. 0) THEN
      WRITE(LUPRI,'(A13,4F18.8)') 'Na B[3] NcNd',
     &        TMPVAL(1),TMPVAL(2),HYPVAL(1),HYPVAL(2)
      END IF
C
C     5. Calculate Nb A[2] Ncd type terms (six permutations)
C     ======================================================
C
C     This is the B A CD permutation
      CALL DZERO(WRK(KGRAD),2*MZYVMX)
      CALL A2INIT(1,KZYVB,KZYVCD,ISYMB,ISPINB,ISYMCD,ISPINCD,1,
     &            WRK(KVECCD),WRK(KGRAD),
     &            XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &            WRK(KFREE),LFREE)
      VAL(1) = DDOT(KZYVB,WRK(KGRAD),1,WRK(KVECB),1)
      VAL(2) = DDOT(KZYVB,WRK(KGRAD),1,WRK(KVECB+KZYVB),1)
      TMPVAL(1) = VAL(1)
      TMPVAL(2) = VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL A2INIT(1,KZYVB,KZYVCD,ISYMB,ISPINB,ISYMCD,ISPINCD,1,
     &               WRK(KVECCD+KZYVCD),WRK(KGRAD),
     &               XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &               WRK(KFREE),LFREE)
         VAL(1) = -DDOT(KZYVB,WRK(KGRAD),1,WRK(KVECB+KZYVB),1)
         VAL(2) =  DDOT(KZYVB,WRK(KGRAD),1,WRK(KVECB),1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
C     This is the CD A B permutation
      CALL DZERO(WRK(KGRAD),2*MZYVMX)
      CALL A2INIT(1,KZYVCD,KZYVB,ISYMCD,ISPINCD,ISYMB,ISPINB,1,
     &            WRK(KVECB),WRK(KGRAD),
     &            XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &            WRK(KFREE),LFREE)
      VAL(1) = DDOT(KZYVCD,WRK(KGRAD),1,WRK(KVECCD),1)
      VAL(2) = DDOT(KZYVCD,WRK(KGRAD),1,WRK(KVECCD+KZYVCD),1)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL A2INIT(1,KZYVCD,KZYVB,ISYMCD,ISPINCD,ISYMB,ISPINB,1,
     &               WRK(KVECB+KZYVB),WRK(KGRAD),
     &               XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &               WRK(KFREE),LFREE)
         VAL(1) = -DDOT(KZYVCD,WRK(KGRAD),1,WRK(KVECCD+KZYVCD),1)
         VAL(2) =  DDOT(KZYVCD,WRK(KGRAD),1,WRK(KVECCD),1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
C     This is the C A BD permutation
      CALL DZERO(WRK(KGRAD),2*MZYVMX)
      CALL A2INIT(1,KZYVC,KZYVBD,ISYMC,ISPINC,ISYMBD,ISPINBD,1,
     &            WRK(KVECBD),WRK(KGRAD),
     &            XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &            WRK(KFREE),LFREE)
      VAL(1) = DDOT(KZYVC,WRK(KGRAD),1,WRK(KVECC),1)
      VAL(2) = DDOT(KZYVC,WRK(KGRAD),1,WRK(KVECC+KZYVC),1)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL A2INIT(1,KZYVC,KZYVBD,ISYMC,ISPINC,ISYMBD,ISPINBD,1,
     &               WRK(KVECBD+KZYVBD),WRK(KGRAD),
     &               XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &               WRK(KFREE),LFREE)
         VAL(1) = -DDOT(KZYVC,WRK(KGRAD),1,WRK(KVECC+KZYVC),1)
         VAL(2) =  DDOT(KZYVC,WRK(KGRAD),1,WRK(KVECC),1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
C     This is the BD A C permutation
      CALL DZERO(WRK(KGRAD),2*MZYVMX)
      CALL A2INIT(1,KZYVBD,KZYVC,ISYMBD,ISPINBD,ISYMC,ISPINC,1,
     &            WRK(KVECC),WRK(KGRAD),
     &            XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &            WRK(KFREE),LFREE)
      VAL(1) = DDOT(KZYVBD,WRK(KGRAD),1,WRK(KVECBD),1)
      VAL(2) = DDOT(KZYVBD,WRK(KGRAD),1,WRK(KVECBD+KZYVBD),1)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL A2INIT(1,KZYVBD,KZYVC,ISYMBD,ISPINBD,ISYMC,ISPINC,1,
     &               WRK(KVECC+KZYVC),WRK(KGRAD),
     &               XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &               WRK(KFREE),LFREE)
         VAL(1) = -DDOT(KZYVBD,WRK(KGRAD),1,WRK(KVECBD+KZYVBD),1)
         VAL(2) =  DDOT(KZYVBD,WRK(KGRAD),1,WRK(KVECBD),1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
C     This is the D A BC permutation
      CALL DZERO(WRK(KGRAD),2*MZYVMX)
      CALL A2INIT(1,KZYVD,KZYVBC,ISYMD,ISPIND,ISYMBC,ISPINBC,1,
     &            WRK(KVECBC),WRK(KGRAD),
     &            XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &            WRK(KFREE),LFREE)
      VAL(1) = DDOT(KZYVD,WRK(KGRAD),1,WRK(KVECD),1)
      VAL(2) = DDOT(KZYVD,WRK(KGRAD),1,WRK(KVECD+KZYVD),1)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL A2INIT(1,KZYVD,KZYVBC,ISYMD,ISPIND,ISYMBC,ISPINBC,1,
     &               WRK(KVECBC+KZYVBC),WRK(KGRAD),
     &               XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &               WRK(KFREE),LFREE)
         VAL(1) = -DDOT(KZYVD,WRK(KGRAD),1,WRK(KVECD+KZYVD),1)
         VAL(2) =  DDOT(KZYVD,WRK(KGRAD),1,WRK(KVECD),1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
C     This is the BC A D permutation
      CALL DZERO(WRK(KGRAD),2*MZYVMX)
      CALL A2INIT(1,KZYVBC,KZYVD,ISYMBC,ISPINBC,ISYMD,ISPIND,1,
     &            WRK(KVECD),WRK(KGRAD),
     &            XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &            WRK(KFREE),LFREE)
      VAL(1) = DDOT(KZYVBC,WRK(KGRAD),1,WRK(KVECBC),1)
      VAL(2) = DDOT(KZYVBC,WRK(KGRAD),1,WRK(KVECBC+KZYVBC),1)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (ABS_DAMP .GT. D0) THEN
         CALL DZERO(WRK(KGRAD),2*MZYVMX)
         CALL A2INIT(1,KZYVBC,KZYVD,ISYMBC,ISPINBC,ISYMD,ISPIND,1,
     &               WRK(KVECD+KZYVD),WRK(KGRAD),
     &               XINDX,UDV,PV,ALAB,ISYMA,ISPINA,CMO,MJWOP,
     &               WRK(KFREE),LFREE)
         VAL(1) = -DDOT(KZYVBC,WRK(KGRAD),1,WRK(KVECBC+KZYVBC),1)
         VAL(2) =  DDOT(KZYVBC,WRK(KGRAD),1,WRK(KVECBC),1)
         TMPVAL(1) = TMPVAL(1) + VAL(1)
         TMPVAL(2) = TMPVAL(2) + VAL(2)
         HYPVAL(1) = HYPVAL(1) + VAL(1)
         HYPVAL(2) = HYPVAL(2) + VAL(2)
      END IF
C
      IF (IPRABSLRS .GT. 0) THEN
      WRITE(LUPRI,'(A13,4F18.8)') 'Nb A[2] Ncd',
     &        TMPVAL(1),TMPVAL(2),HYPVAL(1),HYPVAL(2)
      END IF
C
C     6. Calculate Nb A[3] Nc Nd type terms (two of six permutations in
C     each call)
C     =================================================================
C
C     This is the B A C D permutation
      CALL A3COMPLEX(KZYVB,      KZYVC,      KZYVD,
     &               ISYMB,      ISYMC,      ISYMD,
     &               WRK(KVECB), WRK(KVECC), WRK(KVECD),
     &               VAL(1),     ALAB,       ISYMA,
     &               XINDX, UDV, CMO, MJWOP,
     &               WRK(KFREE), LFREE)
      VAL(1) = -VAL(1)
      VAL(2) = -VAL(2)
      TMPVAL(1) = VAL(1)
      TMPVAL(2) = VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
C     This is the C A B D permutation
      CALL A3COMPLEX(KZYVC,      KZYVB,      KZYVD,
     &               ISYMC,      ISYMB,      ISYMD,
     &               WRK(KVECC), WRK(KVECB), WRK(KVECD),
     &               VAL(1),     ALAB,       ISYMA,
     &               XINDX, UDV, CMO, MJWOP,
     &               WRK(KFREE), LFREE)
      VAL(1) = -VAL(1)
      VAL(2) = -VAL(2)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
C     This is the D A B C permutation
      CALL A3COMPLEX(KZYVD,      KZYVB,      KZYVC,
     &               ISYMD,      ISYMB,      ISYMC,
     &               WRK(KVECD), WRK(KVECB), WRK(KVECC),
     &               VAL(1),     ALAB,       ISYMA,
     &               XINDX, UDV, CMO, MJWOP,
     &               WRK(KFREE), LFREE)
      VAL(1) = -VAL(1)
      VAL(2) = -VAL(2)
      TMPVAL(1) = TMPVAL(1) + VAL(1)
      TMPVAL(2) = TMPVAL(2) + VAL(2)
      HYPVAL(1) = HYPVAL(1) + VAL(1)
      HYPVAL(2) = HYPVAL(2) + VAL(2)
C
      IF (IPRABSLRS .GT. 0) THEN
      WRITE(LUPRI,'(A13,4F18.8)') 'Nb A[3] NcNd',
     &        TMPVAL(1),TMPVAL(2),HYPVAL(1),HYPVAL(2)
      END IF
C
C     Print essential information on a single line
C     ============================================
C
      WRITE(LUPRI,'(/A4,7A1,A2,3F11.6,2F30.10)') 
     & '@ <<',ALAB(1:1),';',BLAB(1:1),',',CLAB(1:1),',',DLAB(1:1),'>>',
     & FREQB,FREQC,FREQD,HYPVAL(1),HYPVAL(2)
C
C     Save value of response function as <<A;B,C,D>>
C     ==============================================
C
      RES_CRF(ICRF,1) = HYPVAL(1)
      RES_CRF(ICRF,2) = HYPVAL(2)
C
C     Check if current element has any equals by permutational symmetry
C     =================================================================
C
      DO IEQTO=1,ABS_CRF_NEQTO
         IF (ABS_CRF_EQTO(IEQTO) .EQ. ICRF) THEN
      WRITE(LUPRI,'(/A10,7A1,A2,A16,3F11.6,/A26,7A1,A2,A16,3F11.6)')
     & 'Element <<',ALAB(1:1),';',BLAB(1:1),',',CLAB(1:1),',',DLAB(1:1),
     & '>>','with frequencies',FREQB,FREQC,FREQD,
     & 'is identical to element <<',
     & ABS_CRFLAB_EQTO(IEQTO,1)(1:1),';',
     & ABS_CRFLAB_EQTO(IEQTO,2)(1:1),',',
     & ABS_CRFLAB_EQTO(IEQTO,3)(1:1),',',
     & ABS_CRFLAB_EQTO(IEQTO,4)(1:1),'>>',
     & 'with frequencies',ABS_CRFFRQ_EQTO(IEQTO,2),
     & ABS_CRFFRQ_EQTO(IEQTO,3),ABS_CRFFRQ_EQTO(IEQTO,4)
C
       WRITE(LUPRI,'(/A4,7A1,A2,3F11.6,2F26.6)')
     & '@ <<',ABS_CRFLAB_EQTO(IEQTO,1)(1:1),';',
     &        ABS_CRFLAB_EQTO(IEQTO,2)(1:1),',',
     &        ABS_CRFLAB_EQTO(IEQTO,3)(1:1),',',
     &        ABS_CRFLAB_EQTO(IEQTO,4)(1:1),'>>',
     &        ABS_CRFFRQ_EQTO(IEQTO,2),ABS_CRFFRQ_EQTO(IEQTO,3),
     &        ABS_CRFFRQ_EQTO(IEQTO,4),HYPVAL(1),HYPVAL(2)
         END IF
      END DO
C
C     End of loop over ABS_NCRF
      END DO
C
C     Close the file where response vectors are saved while keeping it
      CALL GPCLOSE(LURSPRES,'KEEP')
C
C     End of ABS_CRF
C
      RETURN
      END
C
C
C
      SUBROUTINE ABS_BCDCHK(IBCDEQ,
     &                      FREQB, FREQC, FREQD,
     &                      ISYMB, ISYMC, ISYMD,
     &                      KZYVB, KZYVC, KZYVD,
     &                      VECB,  VECC,  VECD)
C
C ======================================================================
C
C  Purpose : Check if some of the vectors VECB, VECC, and VECD are equal
C            or zero. IBCDEQ indicates the
C            result as in
C            0    at least one vector is equal to zero
C            1    all vectors are unequal
C            2    VECB = VECC
C            3    VECB = VECD
C            4    VECC = VECD
C            24   all vectors are equal
C
C   Author : Tobias Fahleson, 2014
C 
C ======================================================================
C
#include "implicit.h"
C
      PARAMETER ( D0=0.0D0, THRSML=1.0D-9, THRZER=1.0D-14 )
      DIMENSION VECB(*), VECC(*), VECD(*)
C
C     Assume that all vectors are unequal
      IBCDEQ = 1
C
C     Compare vectors VECB, VECC, and VECD
C     ====================================
C
      IF ( ISYMB.EQ.ISYMC .AND. FREQB.EQ.FREQC ) THEN
         VMIN = D0
         DO I = 1,2*KZYVB
            VMIN = VMIN + ABS(VECB(I) - VECC(I))
         END DO
         IF (VMIN .LE. THRZER) THEN
            IBCDEQ = 2
         END IF
      END IF
C
      IF ( ISYMB.EQ.ISYMD .AND. FREQB.EQ.FREQD ) THEN
         VMIN = D0
         DO I = 1,2*KZYVB
            VMIN = VMIN + ABS(VECB(I) - VECD(I))
         END DO
         IF (VMIN .LE. THRZER) THEN
            IBCDEQ = IBCDEQ * 3
         END IF
      END IF
C
      IF ( ISYMC.EQ.ISYMD .AND. FREQC.EQ.FREQD ) THEN
         VMIN = D0
         DO I = 1,2*KZYVC
            VMIN = VMIN + ABS(VECC(I) - VECD(I))
         END DO
         IF (VMIN .LE. THRZER) THEN
            IBCDEQ = IBCDEQ * 4
         END IF
      END IF
C
C     Check if some vector is equal to zero
C     =====================================
C
      VNORM = DNRM2(2*KZYVB,VECB,1)
      IF (VNORM .LT. THRSML) THEN
         IBCDEQ = 0
      END IF
C
      VNORM = DNRM2(2*KZYVC,VECC,1)
      IF (VNORM .LT. THRSML) THEN
         IBCDEQ = 0
      END IF
C
      VNORM = DNRM2(2*KZYVD,VECD,1)
      IF (VNORM .LT. THRSML) THEN
         IBCDEQ = 0
      END IF
C
C     End of ABS_BCDCHK
C
      RETURN
      END
C
C
C

      SUBROUTINE T4COMPLEX(IAEQ,  IBCDCOMPLEX,
     &                     ISYMA, ISYMB, ISYMC,  ISYMD,
     &                     KZYVA, KZYVB, KZYVC,  KZYVD,
     &                     VECA,  VECB,  VECC,   VECD,
     &                            FREQB, FREQC,  FREQD,
     &                     FACTOR, XINDX, UDV, PV, MJWOP, CMO, FC,
     &                     WRK, LWRK, RESULT)
C
C ======================================================================
C
C  Purpose : Calculate T[4] times 3 complex vectors
C
C  Author  : Tobias Fahleson, 2014
C
C ======================================================================
C
#include "implicit.h"
#include "priunit.h"
#include "infvar.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inforb.h"
#include "infdim.h"
#include "abslrs.h"
C
      PARAMETER ( D0=0.0D0, D1=1.0D0, DM1=-1.0D0 )
      DIMENSION VECA(2*KZYVA),VECB(2*KZYVB),VECC(2*KZYVC),VECD(2*KZYVD)
      DIMENSION WRK(LWRK),XINDX(*),UDV(*),PV(*),MJWOP(*),CMO(*),FC(*)
      DIMENSION RESULT(2), TMPVAL(2), IBCDCOMPLEX(2,2,2)
C
      KFREE = 1
      LFREE = LWRK - KFREE
C
      RESULT(1) = D0
      RESULT(2) = D0
C
C     If VECA is zero, return zero contribution
C     =========================================
C
      IF (IAEQ .EQ. 0) RETURN
C
C     Calculate four real and four imaginary terms of T[4] NbNcNd and
C     multiply with Na
C     ===============================================================
C
C     Four real terms
C     ===============
C
C     Na T[4] Re{Nb}Re{Nc}Re{Nd}
      IF (IBCDCOMPLEX(1,1,1) .NE. 0) THEN
      CALL T4DRV(IBCDCOMPLEX(1,1,1),
     &           ISYMA,   ISYMB,   ISYMC,   ISYMD,
     &           VECA(1), VECB(1), VECC(1), VECD(1),
     &                   -FREQB,  -FREQC,  -FREQD,
     &           XINDX,UDV,PV,MJWOP,WRK(KFREE),LFREE,CMO,FC)
      TMPVAL(1) = + DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      TMPVAL(2) = + DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      RESULT(1) = RESULT(1) + DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      RESULT(2) = RESULT(2) + DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      END IF
C
      IF (ABS_DAMP .GT. D0) THEN
C
C     Na T[4] (-)Re{Nb}Im{Nc}Im{Nd}
      IF (IBCDCOMPLEX(1,2,2) .NE. 0) THEN
      CALL T4DRV(IBCDCOMPLEX(1,2,2),
     &           ISYMA,   ISYMB,   ISYMC,         ISYMD,
     &           VECA(1), VECB(1), VECC(1+KZYVC), VECD(1+KZYVD),
     &                   -FREQB,  -FREQC,        -FREQD,
     &           XINDX,UDV,PV,MJWOP,WRK(KFREE),LFREE,CMO,FC)
      TMPVAL(1) = - DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      TMPVAL(2) = - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      RESULT(1) = RESULT(1) - DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      RESULT(2) = RESULT(2) - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      END IF
C
C     Na T[4] (-)Im{Nb}Re{Nc}Im{Nd}
      IF (IBCDCOMPLEX(2,1,2) .NE. 0) THEN
      CALL T4DRV(IBCDCOMPLEX(2,1,2),
     &           ISYMA,   ISYMB,         ISYMC,   ISYMD,
     &           VECA(1), VECB(1+KZYVB), VECC(1), VECD(1+KZYVD),
     &                   -FREQB,        -FREQC,  -FREQD,
     &           XINDX,UDV,PV,MJWOP,WRK(KFREE),LFREE,CMO,FC)
      TMPVAL(1) = - DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      TMPVAL(2) = - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      RESULT(1) = RESULT(1) - DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      RESULT(2) = RESULT(2) - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      END IF
C
C     Na T[4] (-)Im{Nb}Im{Nc}Re{Nd}
      IF (IBCDCOMPLEX(2,2,1) .NE. 0) THEN
      CALL T4DRV(IBCDCOMPLEX(2,2,1),
     &           ISYMA,   ISYMB,         ISYMC,         ISYMD,
     &           VECA(1), VECB(1+KZYVB), VECC(1+KZYVC), VECD(1),
     &                   -FREQB,        -FREQC,        -FREQD,
     &           XINDX,UDV,PV,MJWOP,WRK(KFREE),LFREE,CMO,FC)
      TMPVAL(1) = - DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      TMPVAL(2) = - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      RESULT(1) = RESULT(1) - DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      RESULT(2) = RESULT(2) - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      END IF
C
C     Four imaginary terms
C     ===================
C
C     Na T[4] (+i)Re{Nb}Re{Nc}Im{Nd}
      IF (IBCDCOMPLEX(1,1,2) .NE. 0) THEN
      CALL T4DRV(IBCDCOMPLEX(1,1,2),
     &           ISYMA,   ISYMB,   ISYMC,   ISYMD,
     &           VECA(1), VECB(1), VECC(1), VECD(1+KZYVD),
     &                   -FREQB,  -FREQC,  -FREQD,
     &           XINDX,UDV,PV,MJWOP,WRK(KFREE),LFREE,CMO,FC)
      TMPVAL(1) = - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      TMPVAL(2) = + DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      RESULT(1) = RESULT(1) - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      RESULT(2) = RESULT(2) + DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      END IF
C
C     Na T[4] (+i)Re{Nb}Im{Nc}Re{Nd}
      IF (IBCDCOMPLEX(1,2,1) .NE. 0) THEN
      CALL T4DRV(IBCDCOMPLEX(1,2,1),
     &           ISYMA,   ISYMB,   ISYMC,         ISYMD,
     &           VECA(1), VECB(1), VECC(1+KZYVC), VECD(1),
     &                   -FREQB,  -FREQC,        -FREQD,
     &           XINDX,UDV,PV,MJWOP,WRK(KFREE),LFREE,CMO,FC)
      TMPVAL(1) = - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      TMPVAL(2) = + DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      RESULT(1) = RESULT(1) - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      RESULT(2) = RESULT(2) + DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      END IF
C
C     Na T[4] (+i)Im{Nb}Re{Nc}Re{Nd}
      IF (IBCDCOMPLEX(2,1,1) .NE. 0) THEN
      CALL T4DRV(IBCDCOMPLEX(2,1,1),
     &           ISYMA,   ISYMB,         ISYMC,   ISYMD,
     &           VECA(1), VECB(1+KZYVB), VECC(1), VECD(1),
     &                   -FREQB,        -FREQC,  -FREQD,
     &           XINDX,UDV,PV,MJWOP,WRK(KFREE),LFREE,CMO,FC)
      TMPVAL(1) = - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      TMPVAL(2) = + DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      RESULT(1) = RESULT(1) - DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      RESULT(2) = RESULT(2) + DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      END IF
C
C     Na T[4] (-i)Im{Nb}Im{Nc}Im{Nd}
      IF (IBCDCOMPLEX(2,2,2) .NE. 0) THEN
      CALL T4DRV(IBCDCOMPLEX(2,2,2),
     &           ISYMA,   ISYMB,         ISYMC,         ISYMD,
     &           VECA(1), VECB(1+KZYVB), VECC(1+KZYVC), VECD(1+KZYVD),
     &                   -FREQB,        -FREQC,        -FREQD,
     &           XINDX,UDV,PV,MJWOP,WRK(KFREE),LFREE,CMO,FC)
      TMPVAL(1) = + DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      TMPVAL(2) = - DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      RESULT(1) = RESULT(1) + DDOT(KZYVA,VECA(1+KZYVA),1,WRK(KFREE),1)
      RESULT(2) = RESULT(2) - DDOT(KZYVA,VECA(1),1,WRK(KFREE),1)
      END IF
C
      END IF
C
C     End of T4COMPLEX
C
      RETURN
      END
C
C
C
      SUBROUTINE R4DRV(IAEQ,   IBCDCOMPLEX,
     &                 ISYMA,  ISYMB, ISYMC,  ISYMD,
     &                 KZYVA,  KZYVB, KZYVC,  KZYVD,
     &                 VECA,   VECB,  VECC,   VECD,
     &                 FACTOR, XINDX, UDV,    MJWOP,
     &                 WRK,    LWRK,  RESULT)
C
C ======================================================================
C
C  Purpose : Calculate the R[4] contribution as Nb(bar) S[4] Na(bar)NcNd
C            (see Eq. (80) in JCP 123, 194103 (2005), Norman et al.)
C
C  Author  : Tobias Fahleson, 2014
C
C ======================================================================
C
#include "implicit.h"
#include "priunit.h"
#include "infvar.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inforb.h"
#include "infdim.h"
C
      PARAMETER ( D0=0.0D0, D1=1.0D0 )
      DIMENSION WRK(LWRK),XINDX(*),MJWOP(*),UDV(*)
      DIMENSION VECA(2*KZYVA),VECB(2*KZYVB),VECC(2*KZYVC),VECD(2*KZYVD)
      DIMENSION RESULT(2), IBCDCOMPLEX(2,2,2)
      LOGICAL ISOLATE
      ISOLATE = .TRUE.
C
      RESULT(1) = D0
      RESULT(2) = D0
C
      KVECABAR = 1
      KVECBBAR = KVECABAR + 2*KZYVA
      KVECCBAR = KVECBBAR + 2*KZYVB
      KVECDBAR = KVECCBAR + 2*KZYVC
      KS4      = KVECDBAR + 2*KZYVD
C
C     If any vector is zero, return zero contribution
C     ===============================================
C
      IF (IAEQ .EQ. 0 .OR.
     & IBCDCOMPLEX(1,1,1).EQ.0 .OR. IBCDCOMPLEX(1,1,2).EQ.0 .OR.
     & IBCDCOMPLEX(1,2,1).EQ.0 .OR. IBCDCOMPLEX(1,2,2).EQ.0 .OR.
     & IBCDCOMPLEX(2,1,1).EQ.0 .OR. IBCDCOMPLEX(2,1,2).EQ.0 .OR.
     & IBCDCOMPLEX(2,2,1).EQ.0 .OR. IBCDCOMPLEX(2,2,2).EQ.0) RETURN
C
C     Form the bar-vectors
C     ====================
C
      CALL DCOPY(2*KZYVA,VECA(1),1,WRK(KVECABAR),1)
      CALL DCOPY(2*KZYVB,VECB(1),1,WRK(KVECBBAR),1)
      CALL DCOPY(2*KZYVC,VECC(1),1,WRK(KVECCBAR),1)
      CALL DCOPY(2*KZYVD,VECD(1),1,WRK(KVECDBAR),1)
C
      CALL DSWAP(KZYVA/2,WRK(KVECABAR),1,
     &                   WRK(KVECABAR+KZYVA/2),1)
      CALL DSWAP(KZYVA/2,WRK(KVECABAR+KZYVA),1,
     &                   WRK(KVECABAR+KZYVA+KZYVA/2),1)
C
      CALL DSWAP(KZYVB/2,WRK(KVECBBAR),1,
     &                   WRK(KVECBBAR+KZYVB/2),1)
      CALL DSWAP(KZYVB/2,WRK(KVECBBAR+KZYVB),1,
     &                   WRK(KVECBBAR+KZYVB+KZYVB/2),1)
C
      CALL DSWAP(KZYVC/2,WRK(KVECCBAR),1,
     &                   WRK(KVECCBAR+KZYVC/2),1)
      CALL DSWAP(KZYVC/2,WRK(KVECCBAR+KZYVC),1,
     &                   WRK(KVECCBAR+KZYVC+KZYVC/2),1)
C
      CALL DSWAP(KZYVD/2,WRK(KVECDBAR),1,
     &                   WRK(KVECDBAR+KZYVD/2),1)
      CALL DSWAP(KZYVD/2,WRK(KVECDBAR+KZYVD),1,
     &                   WRK(KVECDBAR+KZYVD+KZYVD/2),1)
C
C     Calculate Na(bar) S[4] Nb(bar)NcNd (six permutations, two in each call)
C     =======================================================================
C
      KFREE = KS4  + 2*KZYVB
      LFREE = LWRK - KFREE
C
C     Nb(bar) S[4] [Na(bar)NcNd + Na(bar)NdNc]
      CALL DZERO(WRK(KS4),2*KZYVB)
      CALL S4COMPLEX(KZYVB,    KZYVA,         KZYVC, KZYVD,
     &               ISYMB,    ISYMA,         ISYMC, ISYMD,
     &               WRK(KS4), WRK(KVECABAR), VECC,  VECD,
     &               XINDX, UDV, MJWOP, WRK(KFREE), LFREE)
      RESULT(1)=RESULT(1) - DDOT(KZYVB,WRK(KVECBBAR),1,WRK(KS4+KZYVB),1)
     &                    - DDOT(KZYVB,WRK(KVECBBAR+KZYVB),1,WRK(KS4),1)
      RESULT(2)=RESULT(2) - DDOT(KZYVB,WRK(KVECBBAR+KZYVB),1,
     &                                 WRK(KS4+KZYVB),1)
     &                    + DDOT(KZYVB,WRK(KVECBBAR),1,WRK(KS4),1)
C
      KFREE = KS4  + 2*KZYVC
      LFREE = LWRK - KFREE
C
C     Nc(bar) S[4] [Na(bar)NbNd + Na(bar)NdNb]
      CALL DZERO(WRK(KS4),2*KZYVC)
      CALL S4COMPLEX(KZYVC,    KZYVA,         KZYVB, KZYVD,
     &               ISYMC,    ISYMA,         ISYMB, ISYMD,
     &               WRK(KS4), WRK(KVECABAR), VECB,  VECD,
     &               XINDX, UDV, MJWOP, WRK(KFREE), LFREE)
      RESULT(1)=RESULT(1) - DDOT(KZYVC,WRK(KVECCBAR),1,WRK(KS4+KZYVC),1)
     &                    - DDOT(KZYVC,WRK(KVECCBAR+KZYVC),1,WRK(KS4),1)
      RESULT(2)=RESULT(2) - DDOT(KZYVC,WRK(KVECCBAR+KZYVC),1,
     &                                 WRK(KS4+KZYVC),1)
     &                    + DDOT(KZYVC,WRK(KVECCBAR),1,WRK(KS4),1)
C
      KFREE = KS4  + 2*KZYVD
      LFREE = LWRK - KFREE
C
C     Nd(bar) S[4] [Na(bar)NbNc + Na(bar)NcNb]
      CALL DZERO(WRK(KS4),2*KZYVD)
      CALL S4COMPLEX(KZYVD,    KZYVA,         KZYVB, KZYVC,
     &               ISYMD,    ISYMA,         ISYMB, ISYMC,
     &               WRK(KS4), WRK(KVECABAR), VECB,  VECC,
     &               XINDX, UDV, MJWOP, WRK(KFREE), LFREE)
      RESULT(1)=RESULT(1) - DDOT(KZYVD,WRK(KVECDBAR),1,WRK(KS4+KZYVD),1)
     &                    - DDOT(KZYVD,WRK(KVECDBAR+KZYVD),1,WRK(KS4),1)
      RESULT(2)=RESULT(2) - DDOT(KZYVD,WRK(KVECDBAR+KZYVD),1,
     &                                 WRK(KS4+KZYVD),1)
     &                    + DDOT(KZYVD,WRK(KVECDBAR),1,WRK(KS4),1)
C
C     Finish by multiplying with minus the damping parameter: -FACTOR
C     ===============================================================
C
      RESULT(1) = -FACTOR*RESULT(1)
      RESULT(2) = -FACTOR*RESULT(2)
C
C     End of R4DRV
C
      RETURN
      END
C
C
C
      SUBROUTINE S4COMPLEX(KZYV1,  KZYV2, KZYV3, KZYV4,
     &                     ISYM1,  ISYM2, ISYM3, ISYM4,
     &                     RESKS4, VEC2,  VEC3,  VEC4,
     &                     XINDX,  UDV,   MJWOP, WRK,   LWRK)
C
C ======================================================================
C
C  Purpose : Calculate S[4] times three complex vectors
C
C  Author  : Tobias Fahleson, 2014
C
C ======================================================================
C
#include "implicit.h"
#include "priunit.h"
#include "infvar.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inforb.h"
#include "infdim.h"
C
      PARAMETER ( D0=0.0D0, D1=1.0D0, DM1=-1.0D0 )
      DIMENSION RESKS4(2*KZYV1)
      DIMENSION VEC2(2*KZYV2), VEC3(2*KZYV3), VEC4(2*KZYV4)
      DIMENSION WRK(LWRK), XINDX(*), MJWOP(*), UDV(*)
C
      KS4   = 1
      KFREE = KS4  + KZYV1
      LFREE = LWRK - KFREE
C
C     Calculate four real and four imaginary terms of S[4] N1N2N3 and
C     add to RESKS4
C     ===============================================================
C
C     Four real terms
C     ===============
C
C     S[4] Re{N1}Re{N2}Re{N3}
      CALL DZERO(WRK(KS4),KZYV1)
      CALL S4INIT(KZYV1, KZYV2,   KZYV3,   KZYV4,
     &            ISYM1, ISYM2,   ISYM3,   ISYM4,
     &                   VEC2(1), VEC3(1), VEC4(1),
     &            WRK(KS4),XINDX,UDV,MJWOP,WRK(KFREE),LFREE)
      CALL DAXPY(KZYV1, D1, WRK(KS4), 1, RESKS4(1), 1)
C
C     S[4] (-)Re{N1}Im{N2}Im{N3}
      CALL DZERO(WRK(KS4),KZYV1)
      CALL S4INIT(KZYV1, KZYV2,   KZYV3,         KZYV4,
     &            ISYM1, ISYM2,   ISYM3,         ISYM4,
     &                   VEC2(1), VEC3(1+KZYV3), VEC4(1+KZYV4),
     &            WRK(KS4),XINDX,UDV,MJWOP,WRK(KFREE),LFREE)
      CALL DAXPY(KZYV1, DM1, WRK(KS4), 1, RESKS4(1), 1)
C
C     S[4] (-)Im{N1}Re{N2}Im{N3}
      CALL DZERO(WRK(KS4),KZYV1)
      CALL S4INIT(KZYV1, KZYV2,         KZYV3,   KZYV4,
     &            ISYM1, ISYM2,         ISYM3,   ISYM4,
     &                   VEC2(1+KZYV2), VEC3(1), VEC4(1+KZYV4),
     &            WRK(KS4),XINDX,UDV,MJWOP,WRK(KFREE),LFREE)
      CALL DAXPY(KZYV1, DM1, WRK(KS4), 1, RESKS4(1), 1)
C
C     S[4] (-)Im{N1}Im{N2}Re{N3}
      CALL DZERO(WRK(KS4),KZYV1)
      CALL S4INIT(KZYV1, KZYV2,         KZYV3,         KZYV4,
     &            ISYM1, ISYM2,         ISYM3,         ISYM4,
     &                   VEC2(1+KZYV2), VEC3(1+KZYV3), VEC4(1),
     &            WRK(KS4),XINDX,UDV,MJWOP,WRK(KFREE),LFREE)
      CALL DAXPY(KZYV1, DM1, WRK(KS4), 1, RESKS4(1), 1)
C
C     Four imaginary terms
C     ====================
C
C     S[4] Re{N1}Re{N2}Im{N3}
      CALL DZERO(WRK(KS4),KZYV1)
      CALL S4INIT(KZYV1, KZYV2,   KZYV3,   KZYV4,
     &            ISYM1, ISYM2,   ISYM3,   ISYM4,
     &                   VEC2(1), VEC3(1), VEC4(1+KZYV4),
     &            WRK(KS4),XINDX,UDV,MJWOP,WRK(KFREE),LFREE)
      CALL DAXPY(KZYV1, D1, WRK(KS4), 1, RESKS4(1+KZYV1), 1)
C
C     S[4] Re{N1}Im{N2}Re{N3}
      CALL DZERO(WRK(KS4),KZYV1)
      CALL S4INIT(KZYV1, KZYV2,   KZYV3,         KZYV4,
     &            ISYM1, ISYM2,   ISYM3,         ISYM4,
     &                   VEC2(1), VEC3(1+KZYV3), VEC4(1),
     &            WRK(KS4),XINDX,UDV,MJWOP,WRK(KFREE),LFREE)
      CALL DAXPY(KZYV1, D1, WRK(KS4), 1, RESKS4(1+KZYV1), 1)
C
C     S[4] Im{N1}Re{N2}Re{N3}
      CALL DZERO(WRK(KS4),KZYV1)
      CALL S4INIT(KZYV1, KZYV2,         KZYV3,   KZYV4,
     &            ISYM1, ISYM2,         ISYM3,   ISYM4,
     &                   VEC2(1+KZYV2), VEC3(1), VEC4(1),
     &            WRK(KS4),XINDX,UDV,MJWOP,WRK(KFREE),LFREE)
      CALL DAXPY(KZYV1, D1, WRK(KS4), 1, RESKS4(1+KZYV1), 1)
C
C     S[4] (-)Im{N1}Im{N2}Im{N3}
      CALL DZERO(WRK(KS4),KZYV1)
      CALL S4INIT(KZYV1, KZYV2,         KZYV3,   KZYV4,
     &            ISYM1, ISYM2,         ISYM3,   ISYM4,
     &                   VEC2(1+KZYV2), VEC3(1+KZYV3), VEC4(1+KZYV4),
     &            WRK(KS4),XINDX,UDV,MJWOP,WRK(KFREE),LFREE)
      CALL DAXPY(KZYV1, DM1, WRK(KS4), 1, RESKS4(1+KZYV1), 1)
C
C     End of S4COMPLEX
C
      RETURN
      END
C
C
C
      SUBROUTINE T3COMPLEX(KZYV1,  KZYV2, KZYV3,
     &                     ISYM1,  ISYM2, ISYM3,
     &                     VEC1,   VEC2,  VEC3,
     &                     RESKT3, FREQ2, FREQ3,
     &                     FACTOR, XINDX, UDV, PV, MJWOP, CMO, FC, FV,
     &                     WRK, LWRK)
C
C ======================================================================
C
C  Purpose : Calculate N1 (T[3] - i*gamma*R[3]) N2N3 with
C            complex vectors
C
C  Author  : Tobias Fahleson, 2014
C
C ======================================================================
C
#include "implicit.h"
#include "priunit.h"
#include "infvar.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inforb.h"
#include "infdim.h"
#include "abslrs.h"
C
      PARAMETER ( D0=0.0D0, D1=1.0D0, DM1=-1.0D0 )
      DIMENSION VEC1(2*KZYV1), VEC2(2*KZYV2), VEC3(2*KZYV3)
      DIMENSION WRK(LWRK), RESKT3(2)
      DIMENSION XINDX(*), UDV(*), PV(*), MJWOP(*), CMO(*), FC(*), FV(*)
C
      KT3   = 1
      KFREE = KT3  + KZYV1
      LFREE = LWRK - KFREE
C
      RESKT3(1) = D0
      RESKT3(2) = D0
C
C     Calculate T[3] N2N3 and multiply with complex vector N1
C     =======================================================
C
C     (Re{N1} + iIm{N1}) T[3] Re{N2}Re{N3}
      CALL DZERO(WRK(KT3),KZYV1)
      CALL T3DRV(1,ISYM1,ISYM2,ISYM3,VEC2(1),VEC3(1),
     &           .FALSE.,VEC1,-FREQ2,-FREQ3,XINDX,UDV,PV,MJWOP,
     &           WRK(KT3),LFREE,CMO,FC,FV)
      RESKT3(1) = RESKT3(1) + DDOT(KZYV1,VEC1(1),1,WRK(KT3),1)
      RESKT3(2) = RESKT3(2) + DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KT3),1)
C
      IF (ABS_DAMP .GT. D0) THEN
C
C     (Re{N1} + iIm{N1}) T[3] (-)Im{N2}Im{N3}
      CALL DZERO(WRK(KT3),KZYV1)
      CALL T3DRV(1,ISYM1,ISYM2,ISYM3,VEC2(1+KZYV2),VEC3(1+KZYV3),
     &           .FALSE.,VEC1,-FREQ2,-FREQ3,XINDX,UDV,PV,MJWOP,
     &           WRK(KT3),LFREE,CMO,FC,FV)
      RESKT3(1) = RESKT3(1) + DDOT(KZYV1,VEC1(1),1,WRK(KT3),1)
      RESKT3(2) = RESKT3(2) - DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KT3),1)
C
C     (Re{N1} + iIm{N1}) T[3] i*Im{N2}Re{N3}
      CALL DZERO(WRK(KT3),KZYV1)
      CALL T3DRV(1,ISYM1,ISYM2,ISYM3,VEC2(1+KZYV2),VEC3(1),
     &           .FALSE.,VEC1,-FREQ2,-FREQ3,XINDX,UDV,PV,MJWOP,
     &           WRK(KT3),LFREE,CMO,FC,FV)
      RESKT3(1) = RESKT3(1) - DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KT3),1)
      RESKT3(2) = RESKT3(2) + DDOT(KZYV1,VEC1(1),1,WRK(KT3),1)
C
C     (Re{N1} + iIm{N1}) T[3] i*Re{N2}Im{N3}
      CALL DZERO(WRK(KT3),KZYV1)
      CALL T3DRV(1,ISYM1,ISYM2,ISYM3,VEC2(1),VEC3(1+KZYV3),
     &           .FALSE.,VEC1,-FREQ2,-FREQ3,XINDX,UDV,PV,MJWOP,
     &           WRK(KT3),LFREE,CMO,FC,FV)
      RESKT3(1) = RESKT3(1) - DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KT3),1)
      RESKT3(2) = RESKT3(2) + DDOT(KZYV1,VEC1(1),1,WRK(KT3),1)
C
      END IF
C
C     End of T3COMPLEX
C
      RETURN
      END
C
C
C
      SUBROUTINE X3COMPLEX(KZYV1,  KZYV2,   KZYV3,
     &                     ISYM1,  ISYM2,   ISYM3,
     &                     VEC1,   VEC2,    VEC3,
     &                     RESKX3, GRADLAB, ISYMGRAD,
     &                     XINDX, UDV, CMO, MJWOP,
     &                     WRK, LWRK)
C
C ======================================================================
C
C  Purpose : Calculate N1 X[3] N2N3 with complex vectors
C
C  Author  : Tobias Fahleson, 2014
C
C ======================================================================
C
#include "implicit.h"
#include "priunit.h"
#include "infvar.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inforb.h"
#include "infdim.h"
#include "abslrs.h"
C
      PARAMETER ( D0=0.0D0, D1=1.0D0, DM1=-1.0D0 )
      DIMENSION VEC1(2*KZYV1), VEC2(2*KZYV2), VEC3(2*KZYV3)
      DIMENSION WRK(LWRK), RESKX3(2)
      DIMENSION XINDX(*), UDV(*),  CMO(*), MJWOP(*)
      CHARACTER*8 GRADLAB
C
      KGRAD = 1
      KFREE = KGRAD + KZYV1
      LFREE = LWRK  - KFREE
C
      RESKX3(1) = D0
      RESKX3(2) = D0
C
C     Calculate X[3] N2N3 and multiply with complex vector N1
C     =======================================================
C
C     (Re{N1} + iIm{N1}) X[3] Re{N2}Re{N3}
      CALL DZERO(WRK(KGRAD),KZYV1)
      CALL X3INIT(KZYV1, KZYV2, KZYV3, ISYM1, ISYM2, ISYM3,
     &            GRADLAB, ISYMGRAD, VEC2(1), VEC3(1),
     &            WRK(KGRAD), XINDX, UDV, CMO, MJWOP, WRK(KFREE), LFREE)
      RESKX3(1) = RESKX3(1) + DDOT(KZYV1,VEC1(1),1,WRK(KGRAD),1)
      RESKX3(2) = RESKX3(2) + DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KGRAD),1)
C
      IF (ABS_DAMP .GT. D0) THEN
C
C     (Re{N1} + iIm{N1}) X[3] (-)Im{N2}Im{N3}
      CALL DZERO(WRK(KGRAD),KZYV1)
      CALL X3INIT(KZYV1, KZYV2, KZYV3, ISYM1, ISYM2, ISYM3,
     &            GRADLAB, ISYMGRAD, VEC2(1+KZYV2), VEC3(1+KZYV3), 
     &            WRK(KGRAD), XINDX, UDV, CMO, MJWOP, WRK(KFREE), LFREE)
      RESKX3(1) = RESKX3(1) - DDOT(KZYV1,VEC1(1),1,WRK(KGRAD),1)
      RESKX3(2) = RESKX3(2) - DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KGRAD),1)
C
C     (Re{N1} + iIm{N1}) X[3] i*Im{N2}Re{N3}
      CALL DZERO(WRK(KGRAD),KZYV1)
      CALL X3INIT(KZYV1, KZYV2, KZYV3, ISYM1, ISYM2, ISYM3,
     &            GRADLAB, ISYMGRAD, VEC2(1+KZYV2), VEC3(1),
     &            WRK(KGRAD), XINDX, UDV, CMO, MJWOP, WRK(KFREE), LFREE)
      RESKX3(1) = RESKX3(1) - DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KGRAD),1)
      RESKX3(2) = RESKX3(2) + DDOT(KZYV1,VEC1(1),1,WRK(KGRAD),1)
C
C     (Re{N1} + iIm{N1}) X[3] i*Re{N2}Im{N3}
      CALL DZERO(WRK(KGRAD),KZYV1)
      CALL X3INIT(KZYV1, KZYV2, KZYV3, ISYM1, ISYM2, ISYM3,
     &            GRADLAB, ISYMGRAD, VEC2(1), VEC3(1+KZYV3),
     &            WRK(KGRAD), XINDX, UDV, CMO, MJWOP, WRK(KFREE), LFREE)
      RESKX3(1) = RESKX3(1) - DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KGRAD),1)
      RESKX3(2) = RESKX3(2) + DDOT(KZYV1,VEC1(1),1,WRK(KGRAD),1)
C
      END IF
C
C     End of X3COMPLEX
C
      RETURN
      END
C
C
C
      SUBROUTINE A3COMPLEX(KZYV1,  KZYV2,   KZYV3,
     &                     ISYM1,  ISYM2,   ISYM3,
     &                     VEC1,   VEC2,    VEC3,
     &                     RESKA3, GRADLAB, ISYMGRAD,
     &                     XINDX, UDV, CMO, MJWOP,
     &                     WRK, LWRK)
C
C ======================================================================
C
C  Purpose : Calculate N1 A[3] N2N3 with complex vectors
C
C  Author  : Tobias Fahleson, 2014
C
C ======================================================================
C
#include "implicit.h"
#include "priunit.h"
#include "infvar.h"
#include "wrkrsp.h"
#include "infrsp.h"
#include "inforb.h"
#include "infdim.h"
#include "abslrs.h"
C
      PARAMETER ( D0=0.0D0, D1=1.0D0, DM1=-1.0D0 )
      DIMENSION VEC1(2*KZYV1), VEC2(2*KZYV2), VEC3(2*KZYV3)
      DIMENSION WRK(LWRK), RESKA3(2)
      DIMENSION XINDX(*), UDV(*), CMO(*), MJWOP(*)
      CHARACTER*8 GRADLAB
C
      KGRAD = 1
      KFREE = KGRAD + KZYV1
      LFREE = LWRK  - KFREE
C
      RESKA3(1) = D0
      RESKA3(2) = D0
C
C     Calculate A[3] N2N3 and multiply with complex vector N1
C     =======================================================
C
C     (Re{N1} + iIm{N1}) A[3] Re{N2}Re{N3}
      CALL A3INIT(KZYV1, KZYV2, KZYV3, ISYM1, ISYM2, ISYM3,
     &            GRADLAB, ISYMGRAD, VEC2(1), VEC3(1),
     &            WRK(KGRAD), XINDX, UDV, CMO, MJWOP, WRK(KFREE), LFREE)
      RESKA3(1) = RESKA3(1) + DDOT(KZYV1,VEC1(1),1,WRK(KGRAD),1)
      RESKA3(2) = RESKA3(2) + DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KGRAD),1)
C
      IF (ABS_DAMP .GT. D0) THEN
C
C     (Re{N1} + iIm{N1}) A[3] (-)Im{N2}Im{N3}
      CALL A3INIT(KZYV1, KZYV2, KZYV3, ISYM1, ISYM2, ISYM3,
     &            GRADLAB, ISYMGRAD, VEC2(1+KZYV2), VEC3(1+KZYV3),
     &            WRK(KGRAD), XINDX, UDV, CMO, MJWOP, WRK(KFREE), LFREE)
      RESKA3(1) = RESKA3(1) - DDOT(KZYV1,VEC1(1),1,WRK(KGRAD),1)
      RESKA3(2) = RESKA3(2) - DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KGRAD),1)
C
C     (Re{N1} + iIm{N1}) A[3] i*Im{N2}Re{N3}
      CALL A3INIT(KZYV1, KZYV2, KZYV3, ISYM1, ISYM2, ISYM3,
     &            GRADLAB, ISYMGRAD, VEC2(1+KZYV2), VEC3(1),
     &            WRK(KGRAD), XINDX, UDV, CMO, MJWOP, WRK(KFREE), LFREE)
      RESKA3(1) = RESKA3(1) - DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KGRAD),1)
      RESKA3(2) = RESKA3(2) + DDOT(KZYV1,VEC1(1),1,WRK(KGRAD),1)
C
C     (Re{N1} + iIm{N1}) A[3] i*Re{N2}Im{N3}
      CALL A3INIT(KZYV1, KZYV2, KZYV3, ISYM1, ISYM2, ISYM3,
     &            GRADLAB, ISYMGRAD, VEC2(1), VEC3(1+KZYV3),
     &            WRK(KGRAD), XINDX, UDV, CMO, MJWOP, WRK(KFREE), LFREE)
      RESKA3(1) = RESKA3(1) - DDOT(KZYV1,VEC1(1+KZYV1),1,WRK(KGRAD),1)
      RESKA3(2) = RESKA3(2) + DDOT(KZYV1,VEC1(1),1,WRK(KGRAD),1)
C
      END IF
C
C     End of A3COMPLEX
C
      RETURN
      END
C
C
C
      SUBROUTINE ABS_CRRDVE(ISYMA,  ISYMB,  ISYMC,  ISYMD,
     &                      ISYMBC, ISYMBD, ISYMCD,
     &                      ALAB,   BLAB,   CLAB,   DLAB,
     &                      FREQA,  FREQB,  FREQC,  FREQD,
     &                      KZYVA,  KZYVB,  KZYVC,  KZYVD,
     &                      KZYVBC, KZYVBD, KZYVCD,
     &                      VECA,   VECB,   VECC,   VECD,
     &                      VECBC,  VECBD,  VECCD)
C
C ======================================================================
C
C  Purpose : Read Nx and Nxy vectors from file and return them
C
C   Author : Tobias Fahleson, 2014
C
C ======================================================================
C
#include "implicit.h"
#include "priunit.h"
#include "absorp.h"
#include "abslrs.h"
#include "abscrs.h"
#include "inftap.h"
#include "rspprp.h"
#include "inflr.h"
#include "qrinf.h"
#include "wrkrsp.h"
C
      LOGICAL FOUND, CONV
      CHARACTER*8 ALAB, BLAB, CLAB, DLAB, BLANK
      PARAMETER ( D0=0.0D0, DM1=-1.0D0 )
      PARAMETER (BLANK='        ')
      DIMENSION VECA(*), VECB(*), VECC(*), VECD(*)
      DIMENSION VECBC(*), VECBD(*), VECCD(*)
C
C     Set the length of the vectors Na to Nd, Nbc, Nbd, and Ncd
      KZYVA  = MZYVAR(ISYMA)
      KZYVB  = MZYVAR(ISYMB)
      KZYVC  = MZYVAR(ISYMC)
      KZYVD  = MZYVAR(ISYMD)
      KZYVBC = MZYVAR(ISYMBC)
      KZYVBD = MZYVAR(ISYMBD)
      KZYVCD = MZYVAR(ISYMCD)
C
C     Read in Na
C     ==========
C
c      CALL READ_XVEC(LUABSVECS, 2*KZYVA, VECA, ALAB, ISYMA,
c     &               ABS(FREQA), ABS_THCLR, FOUND, CONV)
      CALL READ_XVEC2(LUABSVECS, 2*KZYVA, VECA, ALAB, BLANK, ISYMA,
     &               ABS(FREQA), D0, ABS_THCLR, FOUND, CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F8.5,A,I3,/A)') ' Response label ',ALAB,
     &           ' with frequency ',FREQA, ' and symmetry',
     &           ISYMA,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F8.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',ALAB,
     &           ' with frequency ',FREQA, ' and symmetry',
     &           ISYMA,' not converged on file LUABSVECS'
         END IF
      END IF
      IF (FREQA .GT. D0) THEN
         CALL DSWAP(KZYVA/2,VECA,1,VECA(1+KZYVA/2),1)
         CALL DSWAP(KZYVA/2,VECA(1+KZYVA),1,VECA(1+KZYVA+KZYVA/2),1)
         CALL DSCAL(KZYVA,DM1,VECA,1)
      END IF
C
C     Read in Nb
C     ==========
C
c      CALL READ_XVEC(LUABSVECS,2*KZYVB,VECB,BLAB,ISYMB,
c     &              ABS(FREQB),ABS_THCLR,FOUND,CONV)
      CALL READ_XVEC2(LUABSVECS,2*KZYVB,VECB,BLAB,BLANK,ISYMB,
     &              ABS(FREQB),D0,ABS_THCLR,FOUND,CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F8.5,A,I3,/A)') ' Response label ',BLAB,
     &           ' with frequency ',FREQB, ' and symmetry',
     &           ISYMB,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F8.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',BLAB,
     &           ' with frequency ',FREQB, ' and symmetry',
     &           ISYMB,' not converged on file LUABSVECS'
         END IF
      END IF
      IF (FREQB .LT. D0) THEN
         CALL DSWAP(KZYVB/2,VECB,1,VECB(1+KZYVB/2),1)
         CALL DSWAP(KZYVB/2,VECB(1+KZYVB),1,VECB(1+KZYVB+KZYVB/2),1)
         CALL DSCAL(KZYVB,DM1,VECB,1)
      END IF
C
C     Read in Nc
C     ==========
C
c      CALL READ_XVEC(LUABSVECS,2*KZYVC,VECC,CLAB,ISYMC,
c     &              ABS(FREQC),ABS_THCLR,FOUND,CONV)
      CALL READ_XVEC2(LUABSVECS,2*KZYVC,VECC,CLAB,BLANK,ISYMC,
     &              ABS(FREQC),D0,ABS_THCLR,FOUND,CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F8.5,A,I3,/A)') ' Response label ',CLAB,
     &           ' with frequency ',FREQC, ' and symmetry',
     &           ISYMC,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F8.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',CLAB,
     &           ' with frequency ',FREQC, ' and symmetry',
     &           ISYMC,' not converged on file LUABSVECS'
         END IF
      END IF
      IF (FREQC .LT. D0) THEN
         CALL DSWAP(KZYVC/2,VECC,1,VECC(1+KZYVC/2),1)
         CALL DSWAP(KZYVC/2,VECC(1+KZYVC),1,VECC(1+KZYVC+KZYVC/2),1)
         CALL DSCAL(KZYVC,DM1,VECC,1)
      END IF
C
C     Read in Nd
C     ==========
C
c      CALL READ_XVEC(LUABSVECS,2*KZYVD,VECD,DLAB,ISYMD,
c     &              ABS(FREQD),ABS_THCLR,FOUND,CONV)
      CALL READ_XVEC2(LUABSVECS,2*KZYVD,VECD,DLAB,BLANK,ISYMD,
     &              ABS(FREQD),D0,ABS_THCLR,FOUND,CONV)
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/3A,F8.5,A,I3,/A)') ' Response label ',DLAB,
     &           ' with frequency ',FREQD, ' and symmetry',
     &           ISYMD,' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/3A,F8.5,/A,I3,A)') ' @WARNING>>>>'//
     &           ' Response label ',DLAB,
     &           ' with frequency ',FREQD, ' and symmetry',
     &           ISYMD,' not converged on file LUABSVECS'
         END IF
      END IF
      IF (FREQD .LT. D0) THEN
         CALL DSWAP(KZYVD/2,VECD,1,VECD(1+KZYVD/2),1)
         CALL DSWAP(KZYVD/2,VECD(1+KZYVD),1,VECD(1+KZYVD+KZYVD/2),1)
         CALL DSCAL(KZYVD,DM1,VECD,1)
      END IF
C
C     Read in Nbc
C     ===========
C
C     Reading vector according to input
      CALL READ_XVEC2(LUABSVECS,2*KZYVBC,VECBC,BLAB,CLAB,ISYMBC,
     &                FREQB, FREQC, ABS_THCLR,FOUND,CONV)
C
      IF (.NOT. FOUND) THEN
C
C     If not found, attempt to find vector with interchanged
C     operators and frequencies
      CALL READ_XVEC2(LUABSVECS,2*KZYVBC,VECBC,CLAB,BLAB,ISYMBC,
     &                FREQC, FREQB, ABS_THCLR,FOUND,CONV)
      END IF
C
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/4A,2(F10.5),A,2(I3),/A)')
     &            ' Response label ',BLAB, CLAB,' with frequencies ',
     &            FREQB,FREQC,' and symmetries',ISYMB,ISYMC,
     &            ' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/4A,2(F10.5),/A,2(I3),A)') ' @WARNING>>>>'//
     &           ' Response label ',BLAB,CLAB,
     &           ' with frequencies ',FREQB,FREQC,
     &           ' and symmetries', ISYMB, ISYMC,
     &           ' not converged on file LUABSVECS'
         END IF
      END IF
      IF ((FREQB + FREQC) .LT. D0) THEN
         CALL DSWAP(KZYVBC/2,VECBC,1,VECBC(1+KZYVBC/2),1)
         CALL DSWAP(KZYVBC/2,VECBC(1+KZYVBC),
     &              1,VECBC(1+KZYVBC+KZYVBC/2),1)
         CALL DSCAL(KZYVBC,DM1,VECBC,1)
      END IF
C
C     Read in Nbd
C     ===========
C
C     Reading vector according to input
      CALL READ_XVEC2(LUABSVECS,2*KZYVBD,VECBD,BLAB,DLAB,ISYMBD,
     &                FREQB, FREQD, ABS_THCLR,FOUND,CONV)
C
      IF (.NOT. FOUND) THEN
C
C     If not found, attempt to find vector with interchanged
C     operators and frequencies
      CALL READ_XVEC2(LUABSVECS,2*KZYVBD,VECBD,DLAB,BLAB,ISYMBD,
     &                FREQD, FREQB, ABS_THCLR,FOUND,CONV)
      END IF
C
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/4A,2(F10.5),A,2(I3),/A)')
     &            ' Response label ',BLAB, DLAB,' with frequencies ',
     &            FREQB,FREQD,' and symmetries',ISYMB,ISYMD,
     &            ' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/4A,2(F10.5),/A,2(I3),A)') ' @WARNING>>>>'//
     &           ' Response label ',BLAB,DLAB,
     &           ' with frequencies ',FREQB,FREQD,
     &           ' and symmetries', ISYMB, ISYMD,
     &           ' not converged on file LUABSVECS'
         END IF
      END IF
      IF ((FREQB + FREQD) .LT. D0) THEN
         CALL DSWAP(KZYVBD/2,VECBD,1,VECBD(1+KZYVBD/2),1)
         CALL DSWAP(KZYVBD/2,VECBD(1+KZYVBD),
     &              1,VECBD(1+KZYVBD+KZYVBD/2),1)
         CALL DSCAL(KZYVBD,DM1,VECBD,1)
      END IF
C
C     Read in Ncd
C     ===========
C
C     Reading vector according to input
      CALL READ_XVEC2(LUABSVECS,2*KZYVCD,VECCD,CLAB,DLAB,ISYMCD,
     &                FREQC, FREQD, ABS_THCLR,FOUND,CONV)
C
      IF (.NOT. FOUND) THEN
C
C     If not found, attempt to find vector with interchanged
C     operators and frequencies
      CALL READ_XVEC2(LUABSVECS,2*KZYVCD,VECCD,DLAB,CLAB,ISYMCD,
     &                FREQD, FREQC, ABS_THCLR,FOUND,CONV)   
      END IF
C
      IF (.NOT. (FOUND .AND. CONV)) THEN
         IF (.NOT. FOUND) THEN
            WRITE (LUPRI,'(/4A,2(F10.5),A,2(I3),/A)')
     &            ' Response label ',CLAB, DLAB,' with frequencies ',
     &            FREQC,FREQD,' and symmetries',ISYMC,ISYMD,
     &            ' not found on file LUABSVECS'
            CALL QUIT('Response vector not found on file')
         ELSE
            WRITE (LUPRI,'(/4A,2(F10.5),/A,2(I3),A)') ' @WARNING>>>>'//
     &           ' Response label ',CLAB,DLAB,
     &           ' with frequencies ',FREQC,FREQD,
     &           ' and symmetries', ISYMC, ISYMD,
     &           ' not converged on file LUABSVECS'
         END IF
      END IF
      IF ((FREQC + FREQD) .LT. D0) THEN
         CALL DSWAP(KZYVCD/2,VECCD,1,VECCD(1+KZYVCD/2),1)
         CALL DSWAP(KZYVCD/2,VECCD(1+KZYVCD),
     &              1,VECCD(1+KZYVCD+KZYVCD/2),1)
         CALL DSCAL(KZYVCD,DM1,VECCD,1)
      END IF
C
C     End of ABS_CRRDVE
C
      RETURN
      END
C
C
C
      SUBROUTINE ABS_ABCDCHK_INIT(IAEQ, IBCDCOMPLEX,
     &      FREQA, FREQB, FREQC, FREQD,
     &      ISYMA, ISYMB, ISYMC, ISYMD,
     &      KZYVA, KZYVB, KZYVC, KZYVD,
     &      VECA,  VECB,  VECC,  VECD)
C
C ======================================================================
C
C  Purpose : Call ABS_BCDCHK for all 8 terms included in complex
C            multiplication of vectors VECB,C,D. Also check if VECA is
C            zero.
C
C   Author : Tobias Fahleson, 2015
C
C ======================================================================
C
#include "implicit.h"
C
      PARAMETER ( THRSML=1.0D-9 )
      DIMENSION IBCDCOMPLEX(2,2,2),VECA(*),VECB(*),VECC(*),VECD(*)
C
C     Check if VECA is zero
C     =====================
C
      VNORM = DNRM2(2*KZYVA,VECA,1)
      IF (VNORM .LT. THRSML) THEN
         IAEQ = 0
      ELSE
         IAEQ = 1
      END IF
C
C     Check vectors VECB,C,D for (in)equality
C     =======================================
C
C     Re{Nb}Re{Nc}Re{Nd}
      CALL ABS_BCDCHK(IBCDCOMPLEX(1,1,1),
     &     FREQB,   FREQC,   FREQD,
     &     ISYMB,   ISYMC,   ISYMD,
     &     KZYVB,   KZYVC,   KZYVD,
     &     VECB(1), VECC(1), VECD(1))
C
C     Re{Nb}Re{Nc}Im{Nd}
      CALL ABS_BCDCHK(IBCDCOMPLEX(1,1,2),
     &     FREQB,   FREQC,   FREQD,
     &     ISYMB,   ISYMC,   ISYMD,
     &     KZYVB,   KZYVC,   KZYVD,
     &     VECB(1), VECC(1), VECD(1+KZYVD))
C
C     Re{Nb}Im{Nc}Re{Nd}
      CALL ABS_BCDCHK(IBCDCOMPLEX(1,2,1),
     &     FREQB,   FREQC,         FREQD,
     &     ISYMB,   ISYMC,         ISYMD,
     &     KZYVB,   KZYVC,         KZYVD,
     &     VECB(1), VECC(1+KZYVC), VECD(1))
C
C     Re{Nb}Im{Nc}Im{Nd}
      CALL ABS_BCDCHK(IBCDCOMPLEX(1,2,2),
     &     FREQB,   FREQC,         FREQD,
     &     ISYMB,   ISYMC,         ISYMD,
     &     KZYVB,   KZYVC,         KZYVD,
     &     VECB(1), VECC(1+KZYVC), VECD(1+KZYVD))
C
C     Im{Nb}Re{Nc}Re{Nd}
      CALL ABS_BCDCHK(IBCDCOMPLEX(2,1,1),
     &     FREQB,         FREQC,   FREQD,
     &     ISYMB,         ISYMC,   ISYMD,
     &     KZYVB,         KZYVC,   KZYVD,
     &     VECB(1+KZYVB), VECC(1), VECD(1))
C
C     Im{Nb}Re{Nc}Im{Nd}
      CALL ABS_BCDCHK(IBCDCOMPLEX(2,1,2),
     &     FREQB,         FREQC,   FREQD,
     &     ISYMB,         ISYMC,   ISYMD,
     &     KZYVB,         KZYVC,   KZYVD,
     &     VECB(1+KZYVB), VECC(1), VECD(1+KZYVD))
C
C     Im{Nb}Im{Nc}Re{Nd}
      CALL ABS_BCDCHK(IBCDCOMPLEX(2,2,1),
     &     FREQB,         FREQC,         FREQD,
     &     ISYMB,         ISYMC,         ISYMD,
     &     KZYVB,         KZYVC,         KZYVD,
     &     VECB(1+KZYVB), VECC(1+KZYVC), VECD(1))
C
C     Im{Nb}Im{Nc}Im{Nd}
      CALL ABS_BCDCHK(IBCDCOMPLEX(2,2,2),
     &     FREQB,         FREQC,         FREQD,
     &     ISYMB,         ISYMC,         ISYMD,
     &     KZYVB,         KZYVC,         KZYVD,
     &     VECB(1+KZYVB), VECC(1+KZYVC), VECD(1+KZYVD))
C
C     End of ABS_BCDCHK_INIT
C
      RETURN
      END
C
C
C
