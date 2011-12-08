C
C  /* Deck rp_eres */
      SUBROUTINE RP_ERES(NOLDTR, NNEWTR, FOCKD, LFOCKD, ISYMTR, 
     &                   WORK,   LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, October 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Driver routine for making a linear transformation of
C              a trialvector with the SOPPA hessian matricx E[2]. 
C              The trial vector consists of two parts TR1E and TR1D. 
C              E refers to excitations and D to 
C              de-excitations. 1 refer to the one-particle part and
C              2 to the two-particle part. The linear transformed 
C              trialvector is refered to as the resultvector and is
C              kept in four corresponding arrays. For the linear
C              transformation with E[2] the result vector is in RES1E,
C              RES1D.
C              The linear transformation is driven over atomic orbitals,
C              and E[2] is not constructed explicitly.
C
#include <implicit.h>
#include <priunit.h>
C
#include <maxorb.h>
#include <maxash.h>
CSPAS:12/11-2003: we need MXCOOR for eribuf.h
#include <mxcent.h>
CKeinSPASmehr
CSPAS:10/11-2003: is merged with maxorb.h
C#include <mxorb.h>
CKeinSPASmehr
#include <aovec.h>
#include <iratdef.h>
#include <eribuf.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION INDEXA(MXCORB)
      DIMENSION FOCKD(LFOCKD)
      DIMENSION WORK(LWORK)
C
#include <ccorb.h>
#include <infind.h>
#include <blocks.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <ccsdio.h>
#include <distcl.h>
#include <cbieri.h>
#include <soppinf.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('RP_ERES')
C
C------------------------------------------------------------------
C     Determine the symmetri of the result vector from the symmetry
C     of the trial vector ISYMTR, and the opperator symmtry ISYMOP.
C------------------------------------------------------------------
C
      ISYRES  = MULD2H(ISYMOP,ISYMTR)
C
C---------------------------------
C     Work space allocation no. 1.
C---------------------------------
C
      LCMO    = NLAMDT
C
      KCMO    = 1
      KEND1   = KCMO   + LCMO  
      LWORK1  = LWORK  - KEND1
C
      CALL SO_MEMMAX ('RP_ERES.1',LWORK1)
      IF (LWORK1 .LT. 0) CALL STOPIT('RP_ERES.1',' ',KEND1,LWORK)
C
C-------------------------------------------------------
C     Get the matrix which contains the MO coefficients.
C-------------------------------------------------------
C
      DTIME      = SECOND()
      CALL SO_GETMO(WORK(KCMO),LCMO,WORK(KEND1),LWORK1)
      DTIME      = SECOND()   - DTIME
      SOTIME(1)  = SOTIME(1) + DTIME
C
C---------------------------------
C     Work space allocation no. 2.
C---------------------------------
C
      LTR1E   = NT1AM(ISYMTR)
      LTR1D   = NT1AM(ISYMTR)
      LRES1E  = NT1AM(ISYMTR)
      LRES1D  = NT1AM(ISYMTR)
      LFOCK   = N2BST(ISYRES)
      LDENS   = N2BST(ISYMTR)
      LBTR1E  = NT1AO(ISYMTR)
      LBTR1D  = NT1AO(ISYMTR)
      LBTJ1E  = NMATAV(ISYMTR)
      LBTJ1D  = NMATAV(ISYMTR)
C
      KTR1E   = KEND1
      KTR1D   = KTR1E   + LTR1E
      KRES1E  = KTR1D   + LTR1D
      KRES1D  = KRES1E  + LRES1E
      KFOCK   = KRES1D  + LRES1D
      KDENS   = KFOCK   + LFOCK
      KBTR1E  = KDENS   + LDENS
      KBTR1D  = KBTR1E  + LBTR1E
      KBTJ1E  = KBTR1D  + LBTR1D
      KBTJ1D  = KBTJ1E  + LBTJ1E
      KEND2   = KBTJ1D  + LBTJ1D
      LWORK2  = LWORK   - KEND2
C
      CALL SO_MEMMAX ('RP_ERES.2',LWORK2)
      IF (LWORK2 .LT. 0) CALL STOPIT('RP_ERES.2',' ',KEND2,LWORK)
C
C----------------------------------------------
C     Open files with trial and result vectors.
C----------------------------------------------
C
      CALL SO_OPEN(LUTR1E,FNTR1E,LTR1E)
      CALL SO_OPEN(LUTR1D,FNTR1D,LTR1D)
      CALL SO_OPEN(LURS1E,FNRS1E,LRES1E)
      CALL SO_OPEN(LURS1D,FNRS1D,LRES1D)
C
      IF ( IPRSOP. GE. 7 ) THEN
C
C------------------------------------------
C        Write new trial vectors to output.
C------------------------------------------
C
         DO 50 INEWTR = 1,NNEWTR
C
C----------------------------------------------------
C           Determine pointer to INEWTR trial vector.
C----------------------------------------------------
C
            INEW = NOLDTR + INEWTR
C
            CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
            CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
C
            WRITE(LUPRI,'(/,I3,A)') 
     &            INEWTR,'. new trial vector in RP_ERES'
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KTR1E+I-1),WORK(KTR1D+I-1),I=1,LTR1E)
C
   50    CONTINUE
C
      END IF
C
C================================================
C     Loop over number of excitations considered.
C================================================
C
      DO 100 INEWTR = 1,NNEWTR
C
C-------------------------------------------------
C        Determine pointer to INEWTR trial vector.
C-------------------------------------------------
C
         INEW = NOLDTR + INEWTR
C
C----------------------------------------
C        Initialize RES1E, RES1D and FOCK
C----------------------------------------
C
         CALL DZERO(WORK(KRES1E),LRES1E)
         CALL DZERO(WORK(KRES1D),LRES1D)
C
         CALL DZERO(WORK(KFOCK),LFOCK)
C
C--------------------------
C        Read trial vector.
C--------------------------
C
         CALL SO_READ(WORK(KTR1E),LTR1E,LUTR1E,FNTR1E,INEW)
         CALL SO_READ(WORK(KTR1D),LTR1D,LUTR1D,FNTR1D,INEW)
C
C---------------------------------------------------
C        Calculate RPA-density matrices in AO basis.
C---------------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_AODENS(WORK(KDENS),LDENS,WORK(KCMO),LCMO,
     &                  WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,ISYMTR,
     &                  WORK(KEND2),LWORK2)
         DTIME     = SECOND()  - DTIME
         SOTIME(6) = SOTIME(6) + DTIME
C
C--------------------------------------------
C        Backtransformation of trial vectors.
C--------------------------------------------
C
         DTIME     = SECOND()
         CALL SO_BCKTR(WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,WORK(KBTR1E),
     &                 LBTR1E,WORK(KBTR1D),LBTR1D,WORK(KBTJ1E),LBTJ1E,
     &                 WORK(KBTJ1D),LBTJ1D,WORK(KCMO),LCMO,ISYMTR)
         DTIME     = SECOND()  - DTIME
         SOTIME(7) = SOTIME(7) + DTIME
C
C=======================================================
C        Start the loop over distributions of integrals.
C=======================================================
C
         IF (DIRECT) THEN
            NTOSYM = 1
            DTIME     = SECOND()
C           CALL HERDI1(WORK(KEND2),LWORK2,IPRINT)
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND2,
     &                  WORK(KEND2),LWORK2,IPRERI)
            KEND2  = KFREE
            LWORK2 = LFREE
            DTIME     = SECOND()  - DTIME
            SOTIME(8) = SOTIME(8) + DTIME
         ELSE
            NTOSYM = NSYM
         ENDIF
C
         KENDSV  = KEND2
         LWORKSV = LWORK2
C
         ICDEL1 = 0
C
         DO 210 ISYMD1 = 1,NTOSYM
C
            IF (DIRECT) THEN
               NTOT = MXCALL
            ELSE
               NTOT = NBAS(ISYMD1)
            ENDIF
C
            DO 220 ILLL = 1,NTOT
C
C------------------------------------------------
C              If direct calculate the integrals.
C------------------------------------------------
C
               IF (DIRECT) THEN
C
                  KEND2  = KENDSV
                  LWORK2 = LWORKSV
C
                  DTIME     = SECOND()
C                 CALL HERDI2(WORK(KEND2),LWORK2,INDEXA,ILLL,NUMDIS,
C    &                             IPRINT)
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,
     &                        WORK(KODCL1),WORK(KODCL2),
     &                        WORK(KODBC1),WORK(KODBC2),
     &                        WORK(KRDBC1),WORK(KRDBC2),
     &                        WORK(KODPP1),WORK(KODPP2),
     &                        WORK(KRDPP1),WORK(KRDPP2),
     &                        WORK(KEND2), LWORK2,IPRERI)
                  DTIME     = SECOND()  - DTIME
                  SOTIME(9) = SOTIME(9) + DTIME
C
                  LRECNR = ( (NBUF -1) / IRAT ) + 1
C
                  KRECNR  = KEND2
                  KEND2   = KRECNR + LRECNR
                  LWORK2  = LWORK  - KEND2
C
                  CALL SO_MEMMAX ('RP_ERES.2B',LWORK2)
                  IF (LWORK2 .LT. 0) 
     &               CALL STOPIT('RP_ERES.2B',' ',KEND2,LWORK)
C
               ELSE
                  NUMDIS = 1
               ENDIF
C
C--------------------------------------------------------
C              Loop over number of distributions in disk.
C--------------------------------------------------------
C
               DO 230 IDEL2 = 1,NUMDIS
C
                  IF (DIRECT) THEN
                     IDEL  = INDEXA(IDEL2)
                     ISYMD = ISAO(IDEL)
                  ELSE
                     IDEL  = IBAS(ISYMD1) + ILLL
                     ISYMD = ISYMD1
                  ENDIF
C
                  ISYDIS = MULD2H(ISYMD,ISYMOP)
C
                  IT2DEL(IDEL) = ICDEL1
                  ICDEL1       = ICDEL1 + NT2BCD(ISYDIS)
C
C---------------------------------------------
C                 Work space allocation no. 3.
C---------------------------------------------
C
                  LXINT  = NDISAO(ISYDIS)
C
                  KXINT   = KEND2
                  KEND3   = KXINT  + LXINT
                  LWORK3  = LWORK  - KEND3
C
                  CALL SO_MEMMAX ('RP_ERES.3',LWORK3)
                  IF (LWORK3 .LT. 0) 
     &               CALL STOPIT('RP_ERES.3',' ',KEND3,LWORK)
C
C--------------------------------------------
C                 Read in batch of integrals.
C--------------------------------------------
C
                  DTIME      = SECOND()
CSPAS:14/11-2003:I am not sure whether it should be IDEL2 here or not,
Cbecause it is a new argument
                  CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND3),LWORK3,
     &                        WORK(KRECNR),DIRECT)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(10) = SOTIME(10) + DTIME
C
C---------------------------------------------
C                 Work space allocation no. 4.
C---------------------------------------------
C
                  LDSRHF = NDSRHF(ISYMD)
C
                  KDSRHF  = KEND3
                  KEND4   = KDSRHF + LDSRHF
                  LWORK4  = LWORK  - KEND4
C
                  CALL SO_MEMMAX ('RP_ERES.4',LWORK4)
                  IF (LWORK4 .LT. 0) 
     &               CALL STOPIT('RP_ERES.4',' ',KEND4,LWORK)
C
C----------------------------------------------------------------
C                 Transform one index in the integral batch to an
C                 occupied index.
C----------------------------------------------------------------
C
                  DTIME      = SECOND()
                  ISYMLP = 1
                  CALL CCTRBT(WORK(KXINT),WORK(KDSRHF),WORK(KCMO),
     &                        ISYMLP,WORK(KEND4),LWORK4,ISYDIS)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(13) = SOTIME(13) + DTIME
C
C----------------------------------------------
C                 Calculate the AO-Fock matrix.
C----------------------------------------------
C
                  DTIME      = SECOND()
                  CALL CC_AOFOCK(WORK(KXINT),WORK(KDENS),WORK(KFOCK),
     &                           WORK(KEND4),LWORK4,IDEL,ISYMD,.FALSE.,
     &                           DUMMY,ISYMTR)
                  DTIME      = SECOND()   - DTIME
                  SOTIME(11) = SOTIME(11) + DTIME
C
  230          CONTINUE
C
  220       CONTINUE
C
  210    CONTINUE
C
C====================================================
C        End of loop over distributions of integrals.
C====================================================
C
C---------------------------------------------
C        Transform AO Fock matrix to MO basis.
C---------------------------------------------
C
         DTIME      = SECOND()
         CALL CC_FCKMO(WORK(KFOCK),WORK(KCMO),WORK(KCMO),
     &                    WORK(KEND2),LWORK2,ISYRES,1,1)
         DTIME      = SECOND()   - DTIME
         SOTIME(24) = SOTIME(24) + DTIME
C
C------------------------------------------------------------------
C        Calculate and add the RPA two-particle parts to the result 
C        vectors.
C------------------------------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_TWOFOCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,
     &                   WORK(KFOCK),LFOCK,ISYRES)
         DTIME      = SECOND()   - DTIME
         SOTIME(25) = SOTIME(25) + DTIME
C
C------------------------------------------------------------------
C        Calculate and add the RPA one-particle parts to the result
C        vectors.
C------------------------------------------------------------------
C
         DTIME      = SECOND()
         CALL SO_ONEFOCK(WORK(KRES1E),LRES1E,WORK(KRES1D),LRES1D,FOCKD,
     &                   LFOCKD,WORK(KTR1E),LTR1E,WORK(KTR1D),LTR1D,
     &                   ISYRES,ISYMTR)
         DTIME      = SECOND()   - DTIME
         SOTIME(26) = SOTIME(26) + DTIME
C
C----------------------------------------
C        Write new resultvectors to file.
C----------------------------------------
C
         CALL SO_WRITE(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
         CALL SO_WRITE(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,INEW)
C
  100 CONTINUE
C
C==================================
C     End of loop over excitations.
C==================================
C
      IF ( IPRSOP. GE. 7 ) THEN
C
C------------------------------------------
C        Write new resultvectors to output.
C------------------------------------------
C
         DO 400 INEWTR = 1,NNEWTR
C
C----------------------------------------------------
C           Determine pointer to INEWTR trial vector.
C----------------------------------------------------
C
            INEW = NOLDTR + INEWTR
C
            WRITE(LUPRI,'(/,I3,A)') INEWTR,
     &         '. new E[2] linear transformed'//
     &         ' trial vector'
C
            CALL SO_READ(WORK(KRES1E),LRES1E,LURS1E,FNRS1E,INEW)
            CALL SO_READ(WORK(KRES1D),LRES1D,LURS1D,FNRS1D,INEW)
C
            WRITE(LUPRI,'(I8,1X,F14.8,5X,F14.8)') 
     &           (I,WORK(KRES1E+I-1),WORK(KRES1D+I-1),I=1,LRES1E)
C
  400    CONTINUE
C
      END IF
C
C-----------------
C     Close files.
C-----------------
C
      CALL SO_CLOSE(LUTR1E,FNTR1E,'KEEP')
      CALL SO_CLOSE(LUTR1D,FNTR1D,'KEEP')
      CALL SO_CLOSE(LURS1E,FNRS1E,'KEEP')
      CALL SO_CLOSE(LURS1D,FNRS1D,'KEEP')
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('RP_ERES')
C
      RETURN
      END
