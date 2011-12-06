C  /* Deck nocc_t2larg */
      SUBROUTINE NOCC_T2LARG(T2AM,ISYM,LAI,ISYM1,LBJ,ISYM2,LAIBJ,T2MX)
C
C     Thomas Bondo Pedersen, August 2001.
C
C     Purpose: Find the largest (absolute value) element and its location
C              in the packed doubles array T2AM of symmetry ISYM.
C
#include <implicit.h>
      INTEGER AI, BJ, AIBJ
#include <ccorb.h>
#include <ccsdsym.h>
C
      DIMENSION T2AM(*)
C
      PARAMETER (ZERO = 0.00D0)
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      T2MX  = ZERO
      LAI   = 0
      ISYM1 = 0
      LBJ   = 0
      ISYM2 = 0
      LAIBJ = 0
C
      IF (ISYM .EQ. 1) THEN
C
         DO ISYMBJ = 1,NSYM
            ISYMAI = ISYMBJ
            DO BJ = 1,NT1AM(ISYMBJ)
               DO AI = 1,BJ
                  AIBJ = IT2AM(ISYMAI,ISYMBJ) + INDEX(AI,BJ)
                  IF (DABS(T2AM(AIBJ)) .GT. DABS(T2MX)) THEN
                     T2MX  = T2AM(AIBJ)
                     LAI   = AI
                     ISYM1 = ISYMAI
                     LBJ   = BJ
                     ISYM2 = ISYMBJ
                     LAIBJ = AIBJ
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
C
      ELSE
C
         DO ISYMBJ = 1,NSYM
            ISYMAI = MULD2H(ISYMBJ,ISYM)
            IF (ISYMAI .LT. ISYMBJ) THEN
               DO BJ = 1,NT1AM(ISYMBJ)
                  DO AI = 1,NT1AM(ISYMAI)
                     AIBJ = IT2AM(ISYMAI,ISYMBJ)
     &                    + NT1AM(ISYMAI)*(BJ - 1)
     &                    + AI
                     IF (DABS(T2AM(AIBJ)) .GT. DABS(T2MX)) THEN
                        T2MX  = T2AM(AIBJ)
                        LAI   = AI
                        ISYM1 = ISYMAI
                        LBJ   = BJ
                        ISYM2 = ISYMBJ
                        LAIBJ = AIBJ
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
C
      ENDIF
C
      RETURN
      END
C  /* Deck nocc_prt */
      SUBROUTINE NOCC_PRT(XMAT,ISYM,TYP)
C
C     Thomas Bondo Pedersen, November 2000.
C
C     Purpose: Print XMAT to output.
C
C     ISYM is the symmetry of XMAT.
C
C     TYP identifies the matrix:
C
C        TYP(1:4) = 'IJ  ' : Print occ.-occ. matrix.
C
C        TYP(1:4) = 'AI  ' : Print vir.-occ. matrix.
C
C        TYP(1:4) = 'IA  ' : Print occ.-vir. matrix.
C
C        TYP(1:4) = 'AB  ' : Print vir.-vir. matrix.
C
C        TYP(1:4) = 'AIBJ' : Print squared doubles array.
C
C        TYP(1:4) = 'SIDO' : Print full singles and squared doubles array.
C
C        TYP(1:4) = 'ALBE' : Print squared (alfa,beta) matrix.
C
C        TYP(1:4) = 'LMDA' : Print Lambda type matrix.
C
C        TYP(1:4) = 'AIJK' : Print (ai,jk) matrix.
C
C     NOTE: No header is written - must be done in calling routine!
C           (Except for 'LMDA' option where occupied and virtual blocks
C            are printed seperately.)
C
#include <implicit.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <noccsym.h>
#include <priunit.h>
C
      DIMENSION XMAT(*)
C
      CHARACTER*(*) TYP
C
      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'NOCC_PRT')
C
C============================
C     Print according to TYP.
C============================
C
      IF (TYP(1:4) .EQ. 'IJ  ') THEN
C
         DO ISYMJ = 1,NSYM
C
            ISYMI = MULD2H(ISYMJ,ISYM)
C
            NI = NRHF(ISYMI)
            NJ = NRHF(ISYMJ)
C
            KOFF = IMATIJ(ISYMI,ISYMJ) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of I: ',ISYMI,' Symmetry of J: ',ISYMJ
C
            IF ((NI.EQ.0) .OR. (NJ.EQ.0)) THEN
               WRITE(LUPRI,'(10X,A,/)')
     &         'Empty block'
            ELSE
               CALL OUTPUT(XMAT(KOFF),1,NI,1,NJ,NI,NJ,1,LUPRI)
            ENDIF
C
         ENDDO
C
      ELSE IF (TYP(1:4) .EQ. 'AI  ') THEN
C
         DO ISYMI = 1,NSYM
C
            ISYMA = MULD2H(ISYMI,ISYM)
C
            NA = NVIR(ISYMA)
            NI = NRHF(ISYMI)
C
            KOFF = IT1AM(ISYMA,ISYMI) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of A: ',ISYMA,' Symmetry of I: ',ISYMI
C
            IF ((NA.EQ.0) .OR. (NI.EQ.0)) THEN
               WRITE(LUPRI,'(10X,A,/)')
     &         'Empty block'
            ELSE
               CALL OUTPUT(XMAT(KOFF),1,NA,1,NI,NA,NI,1,LUPRI)
            ENDIF
C
         ENDDO
C
      ELSE IF (TYP(1:4) .EQ. 'IA  ') THEN
C
         DO ISYMA = 1,NSYM
C
            ISYMI = MULD2H(ISYMA,ISYM)
C
            NA = NVIR(ISYMA)
            NI = NRHF(ISYMI)
C
            KOFF = IT1AMT(ISYMI,ISYMA) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of I: ',ISYMI,' Symmetry of A: ',ISYMA
C
            IF ((NI.EQ.0) .OR. (NA.EQ.0)) THEN
               WRITE(LUPRI,'(10X,A,/)')
     &         'Empty block'
            ELSE
               CALL OUTPUT(XMAT(KOFF),1,NI,1,NA,NI,NA,1,LUPRI)
            ENDIF
C
         ENDDO
C
      ELSE IF (TYP(1:4) .EQ. 'AB  ') THEN
C
         DO ISYMA = 1,NSYM
C
            ISYMB = MULD2H(ISYMA,ISYM)
C
            NA = NVIR(ISYMA)
            NB = NVIR(ISYMB)
C
            KOFF = IMATAB(ISYMA,ISYMB) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of A: ',ISYMA,' Symmetry of B: ',ISYMB
C
            IF ((NA.EQ.0) .OR. (NB.EQ.0)) THEN
               WRITE(LUPRI,'(10X,A,/)')
     &         'Empty block'
            ELSE
               CALL OUTPUT(XMAT(KOFF),1,NA,1,NB,NA,NB,1,LUPRI)
            ENDIF
C
         ENDDO
C
      ELSE IF (TYP(1:4) .EQ. 'AIBJ') THEN
C
         CALL CC_PRSQ(DUMMY,XMAT,ISYM,0,1)
C
      ELSE IF (TYP(1:4) .EQ. 'SIDO') THEN
C
         KOFF1 = 1
         KOFF2 = KOFF1 + NT1AM(ISYM)
C
         CALL CC_PRSQ(XMAT(KOFF1),XMAT(KOFF2),ISYM,1,1)
C
      ELSE IF (TYP(1:4) .EQ. 'ALBE') THEN
C
         DO ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYM)
C
            NA = NBAS(ISYMA)
            NB = NBAS(ISYMB)
C
            KOFF = IAODIS(ISYMA,ISYMB) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of Alpha: ',ISYMA,' Symmetry of Beta: ',ISYMB
C
            IF ((NA.EQ.0) .OR. (NB.EQ.0)) THEN
               WRITE(LUPRI,'(10X,A,/)')
     &         'Empty block'
            ELSE
               CALL OUTPUT(XMAT(KOFF),1,NA,1,NB,NA,NB,1,LUPRI)
            ENDIF
C
         ENDDO
C
      ELSE IF (TYP(1:4) .EQ. 'LMDA') THEN
C
         CALL HEADER('Occupied Part',-1)
C
         DO ISYMMO = 1,NSYM
C
            ISYMAO = MULD2H(ISYMMO,ISYM)
C
            NALPHA = NBAS(ISYMAO)
            NI     = NRHF(ISYMMO)
C
            KOFF = IGLMRH(ISYMAO,ISYMMO) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of Alpha: ',ISYMAO,' Symmetry of I: ',ISYMMO
C
            IF ((NALPHA.EQ.0) .OR. (NI.EQ.0)) THEN
               WRITE(LUPRI,'(10X,A,/)')
     &         'Empty block'
            ELSE
               CALL OUTPUT(XMAT(KOFF),1,NALPHA,1,NI,NALPHA,NI,1,LUPRI)
            ENDIF
C
         ENDDO
C
         CALL HEADER('Virtual Part',-1)
C
         DO ISYMMO = 1,NSYM
C
            ISYMAO = MULD2H(ISYMMO,ISYM)
C
            NALPHA = NBAS(ISYMAO)
            NA     = NVIR(ISYMMO)
C
            KOFF = IGLMVI(ISYMAO,ISYMMO) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of Alpha: ',ISYMAO,' Symmetry of A: ',ISYMMO
C
            IF ((NALPHA.EQ.0) .OR. (NA.EQ.0)) THEN
               WRITE(LUPRI,'(10X,A,/)')
     &         'Empty block'
            ELSE
               CALL OUTPUT(XMAT(KOFF),1,NALPHA,1,NA,NALPHA,NA,1,LUPRI)
            ENDIF
C
         ENDDO
C
      ELSE IF (TYP(1:4) .EQ. 'AIJK') THEN
C
         DO ISYMJK = 1,NSYM
C
            ISYMAI = MULD2H(ISYMJK,ISYM)
C
            NAI = NT1AM(ISYMAI)
            NJK = NMATIJ(ISYMJK)
C
            KOFF = ICAIJK(ISYMAI,ISYMJK) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of AI: ',ISYMAI,' Symmetry of JK: ',ISYMJK
C
            IF ((NAI.EQ.0) .OR. (NJK.EQ.0)) THEN
               WRITE(LUPRI,'(10X,A,/)')
     &         'Empty block'
            ELSE
               CALL OUTPUT(XMAT(KOFF),1,NAI,1,NJK,NAI,NJK,1,LUPRI)
            ENDIF
C
         ENDDO
C
      ELSE
C
C----------------------------
C        Unknown TYP: Return.
C----------------------------
C
         WRITE(LUPRI,'(/,1X,A,A,A4,/)')
     &   SECNAM,': Returning immediately; unknown TYP: ',
     &   TYP(1:4)
C
         RETURN
C
      ENDIF
C
      RETURN
      END
C  /* Deck nocc_prtpk */
      SUBROUTINE NOCC_PRTPK(XMAT,ISYM,TYP)
C
C     Thomas Bondo Pedersen, November 2000.
C
C     Purpose: Print packed XMAT to output.
C
C     ISYM is the symmetry of XMAT.
C
C     TYP identifies the matrix:
C
C        TYP(1:4) = 'AIBJ' : Print packed doubles array (ai .LE. bj).
C
C        TYP(1:4) = 'SIDO' : Print full singles and packed doubles array.
C
C        TYP(1:4) = 'IJKL' : Print packed (ij,kl) matrix  (ij .LE. kl).
C
C        TYP(1:4) = 'ALBE' : Print packed (alfa,beta) matrix (alfa .LE. beta).
C
C     NOTE: No header is written - must be done in calling routine!
C
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <ccsdsym.h>
#include <ccorb.h>
#include <symsq.h>
C
      DIMENSION XMAT(*)
C
      CHARACTER*(*) TYP
C
      CHARACTER*10 SECNAM
      PARAMETER (SECNAM = 'NOCC_PRTPK')
C
C============================
C     Print according to TYP.
C============================
C
      IF (TYP(1:4) .EQ. 'AIBJ') THEN
C
C-----------------------------------
C        Print packed doubles array.
C-----------------------------------
C
         CALL CC_PRP(DUMMY,XMAT,ISYM,0,1)
C
      ELSE IF (TYP(1:4) .EQ. 'SIDO') THEN
C
C--------------------------------------------------
C        Print full singles + packed doubles array.
C--------------------------------------------------
C
         KOFF = NT1AM(ISYM) + 1
C
         CALL CC_PRP(XMAT(1),XMAT(KOFF),ISYM,1,1)
C
      ELSE IF (TYP(1:4) .EQ. 'IJKL') THEN
C
C----------------------------------
C        Print packed ij,kl matrix.
C----------------------------------
C
         DO ISYMKL = 1,NSYM
C
            ISYMIJ = MULD2H(ISYMKL,ISYM)
C
            NIJ = NMATIJ(ISYMIJ)
            NKL = NMATIJ(ISYMKL)
C
            KOFF = IGAMMA(ISYMIJ,ISYMKL) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of IJ: ',ISYMIJ,' Symmetry of KL: ',ISYMKL
C
            IF (ISYMIJ .EQ. ISYMKL) THEN
C
               IF ((NIJ.EQ.0) .OR. (NKL.EQ.0)) THEN
                  WRITE(LUPRI,'(10X,A,/)')
     &            'Empty block'
               ELSE
                  CALL OUTPAK(XMAT(KOFF),NIJ,1,LUPRI)
               ENDIF
C
            ELSE IF (ISYMIJ .LT. ISYMKL) THEN
C
               IF ((NIJ.EQ.0) .OR. (NKL.EQ.0)) THEN
                  WRITE(LUPRI,'(10X,A,/)')
     &            'Empty block'
               ELSE
                  CALL OUTPUT(XMAT(KOFF),1,NIJ,1,NKL,NIJ,NKL,1,LUPRI)
               ENDIF
C
            ENDIF
C
         ENDDO
C
      ELSE IF (TYP(1:4) .EQ. 'ALBE') THEN
C
C----------------------------------------
C        Print packed (alfa,beta) matrix.
C----------------------------------------
C
         DO ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYM)
C
            NA = NBAS(ISYMA)
            NB = NBAS(ISYMB)
C
            KOFF = IAODPK(ISYMA,ISYMB) + 1
C
            WRITE(LUPRI,'(/,5X,A,I1,A,I1,/)')
     &      'Symmetry of Alpha: ',ISYMA,' Symmetry of Beta: ',ISYMB
C
            IF (ISYMA .EQ. ISYMB) THEN
C
               IF (NA .EQ. 0) THEN
                  WRITE(LUPRI,'(10X,A,/)')
     &            'Empty block'
               ELSE
                  CALL OUTPAK(XMAT(KOFF),NA,1,LUPRI)
               ENDIF
C
            ELSE IF (ISYMA .LT. ISYMB) THEN
C
               IF ((NA.EQ.0) .OR. (NB.EQ.0)) THEN
                  WRITE(LUPRI,'(10X,A,/)')
     &            'Empty block'
               ELSE
                  CALL OUTPUT(XMAT(KOFF),1,NA,1,NB,NA,NB,1,LUPRI)
               ENDIF
C
            ENDIF
C
         ENDDO
C
      ELSE
C
C----------------------------
C        Unknown TYP: Return.
C----------------------------
C
         WRITE(LUPRI,'(/,1X,A,A,A4,/)')
     &   SECNAM,': Returning immediately; unknown TYP: ',
     &   TYP(1:4)
C
         RETURN
C
      ENDIF
C
      RETURN
      END
