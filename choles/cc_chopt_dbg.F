C  /* Deck cc_chopt_dbg */
      SUBROUTINE CC_CHOPT_DBG(WORK,LWORK)
C
C     JLC, BFR, TBP, HK, and AS, Oct. 2002.
C
C     Purpose: Driver routine for debugging Cholesky denominator
C              CCSD(T). After setting up the necessary data arrays,
C              control is handed over to CCHO_DBGDRV.
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <priunit.h>
C
      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CC_CHOPT_DBG')
C
C------------------
C     Print header.
C------------------
C
      CALL AROUND('Output from '//SECNAM)
C
C-------------------
C     Symmetry test.
C-------------------
C
      IF (NSYM .NE. 1) THEN
         WRITE(LUPRI,'(5X,A,A,A,/,5X,A,I10,/)')
     &   'ERROR: Symmetry not implemented in ',SECNAM,':',
     &   'NSYM = ',NSYM
         CALL QUIT('Symmetry not implemented in '//SECNAM)
      ENDIF
C
C----------------------
C     Allocation o.int.
C----------------------
C
      KIAJK = 1
      KEND1 = KIAJK + NRHFT*NRHFT*NRHFT*NVIRT
      LWRK1 = LWORK - KEND1 + 1

      write(6,*) 'Check NTRAOC: ',NTRAOC(1),NRHFT*NRHFT*NRHFT*NVIRT
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//': o.int')
      ENDIF
C
C------------------------------------------------------
C     Get integrals (ia|jk) sorted as (ijka) [in core],
C                   (ai|bc) sorted as (aicb) [on disk].
C------------------------------------------------------
C
      CALL CCHO_TRAINT(WORK(KEND1),LWRK1,WORK(KIAJK))

      XNRM = DSQRT(DDOT(NTRAOC(1),WORK(KIAJK),1,WORK(KIAJK),1))
      write(6,*) 'Norm of (ia|jk): ',XNRM
      write(6,*) 'after traint...'
      call flshfo(6)
C
C----------------------
C     Allocation: rest.
C----------------------
C
      KT1AM = KEND1
      KT2AM = KT1AM + NVIRT*NRHFT
      KIAJB = KT2AM + NVIRT*NRHFT*NVIRT*NRHFT
      KAIBC = KIAJB + NVIRT*NRHFT*NVIRT*NRHFT
      KFDIA = KAIBC + NVIRT*NRHFT*NVIRT*NVIRT
      KEND2 = KFDIA + NORBTS
      LWRK2 = LWORK - KEND2 + 1
C
      KSCR  = KAIBC
      KENDS = KSCR  + NT2AM(1)
      LWRKS = LWORK - KENDS + 1
C
      LWRKT = MIN(LWRK2,LWRKS)
      IF (LWRKT .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM//': rest')
      ENDIF
C
C--------------------------------
C     Read amplitudes and square.
C--------------------------------
C
      OPEN(UNIT=53,STATUS='OLD',FORM='UNFORMATTED',FILE='CCSD_TAM')
      REWIND (UNIT=53)
      READ(53) (WORK(KT1AM+I-1), I = 1,NT1AM(1))
      READ(53) (WORK(KSCR+I-1),  I = 1,NT2AM(1))
      CLOSE(53,STATUS='KEEP')
      CALL CC_T2SQ(WORK(KSCR),WORK(KT2AM),1)

      T1NRM = DSQRT(DDOT(NT1AMX,WORK(KT1AM),1,WORK(KT1AM),1))
      T2NRM = DSQRT(DDOT(NT2AMX,WORK(KSCR),1,WORK(KSCR),1))
      write(6,*) 'Norm of t1am: ',T1NRM
      write(6,*) 'Norm of t2am: ',T2NRM

C
C------------------------------------------------
C     Read (ia|jb) [sorted as (aibj)] and square.
C------------------------------------------------
C
      CALL CCG_RDIAJB(WORK(KSCR),NT2AM(1))
      CALL CC_T2SQ(WORK(KSCR),WORK(KIAJB),1)

      XNRM = DSQRT(DDOT(NT2AMX,WORK(KSCR),1,WORK(KSCR),1))
      write(6,*) 'Norm of (ia|jb): ',XNRM

      write(6,*) 'after rdiajb...'
      call flshfo(6)
C
C-------------------------------------
C     Read (ai|bc) [sorted as (aicb)].
C-------------------------------------
C
      CALL CCHO_RDINT(61,'CHO_VI1',WORK(KAIBC),1,NVIRT,1)

      LEN  = NVIRT*NRHFT*NVIRT*NVIRT
      XNRM = DSQRT(DDOT(LEN,WORK(KAIBC),1,WORK(KAIBC),1))
      write(6,*) 'Norm of (ai|bc): ',XNRM
      write(6,*) 'after rdint...'
      call flshfo(6)
C
C----------------------------------------------------
C     Read orbital energies; delete frozen; re-order.
C----------------------------------------------------
C
      OPEN(53,STATUS='OLD',FORM='UNFORMATTED',FILE='SIRIFC')
      REWIND (UNIT=53)
      CALL MOLLAB('TRCCINT ',53,LUPRI)
      READ (53)
      READ (53) (WORK(KFDIA+I-1), I=1,NORBTS)
      CLOSE(53,STATUS='KEEP')
      IF (FROIMP .OR. FROEXP) THEN
         CALL CCSD_DELFRO(WORK(KFDIA),WORK(KEND2),LWRK2)
      ENDIF
      CALL FOCK_REORDER(WORK(KFDIA),WORK(KEND2),LWRK2)
C
C----------------------
C     Get the job done.
C----------------------
C

      write(6,*) 'calling dbg..'
      call flshfo(6)

      KFOCC = KFDIA
      KFVIR = KFOCC + NRHFT
      CALL CCHO_DBGDRV(WORK(KIAJK),WORK(KT1AM),WORK(KT2AM),WORK(KIAJB),
     &                 WORK(KAIBC),WORK(KFOCC),WORK(KFVIR),
     &                 WORK(KEND2),LWRK2)
C
C------------------------
C     Print exit message.
C------------------------
C
      CALL AROUND('End of '//SECNAM)
C
      RETURN
      END
C  /* Deck ccho_dbgdrv */
      SUBROUTINE CCHO_DBGDRV(XIAJK,T1AM,T2AM,XIAJB,XAIBC,FOCKO,FOCKV,
     &                       WORK,LWORK)
C
C     JLC, BFR, TBP, HK, and AS, October 2002.
C
C     Purpose: Debug code for calculating CCSD(T) energy correction.
C
#include <implicit.h>
      DIMENSION XIAJK(*),T1AM(*),T2AM(*),XIAJB(*),XAIBC(*)
      DIMENSION FOCKO(*),FOCKV(*)
      DIMENSION WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CCHO_DBGDRV')
C
      PARAMETER (TWO = 2.00D0)

      write(6,*) 'Entering CCHO_DBGDRV'
      LIAJK  = NRHFT*NRHFT*NRHFT*NVIRT
      LAI    = NVIRT*NRHFT
      LAIBC  = NVIRT*NRHFT*NVIRT*NVIRT
      XNIAJK = DSQRT(DDOT(LIAJK,XIAJK,1,XIAJK,1))
      XNT1AM = DSQRT(DDOT(LAI,T1AM,1,T1AM,1))
      XNT2AM = 0.0D0
      DO NBJ = 1,LAI
         DO NAI = 1,NBJ
            KOFF   = LAI*(NBJ - 1) + NAI
            XNT2AM = XNT2AM + T2AM(KOFF)*T2AM(KOFF)
         ENDDO
      ENDDO
      XNT2AM = DSQRT(XNT2AM)
      XNIAJB = 0.0D0
      DO NBJ = 1,LAI
         DO NAI = 1,NBJ
            KOFF   = LAI*(NBJ - 1) + NAI
            XNIAJB = XNIAJB + XIAJB(KOFF)*XIAJB(KOFF)
         ENDDO
      ENDDO
      XNIAJB = DSQRT(XNIAJB)
      XNAIBC = DSQRT(DDOT(LAIBC,XAIBC,1,XAIBC,1))
      write(6,*) 'Norm of (ia|jk): ',XNIAJK
      write(6,*) 'Norm of t1am   : ',XNT1AM
      write(6,*) 'Norm of t2am   : ',XNT2AM
      write(6,*) 'Norm of (ia|jb): ',XNIAJB
      write(6,*) 'Norm of (ai|bc): ',XNAIBC
C
C---------------------------------------
C     Set up orbital energy differences.
C---------------------------------------
C
      KDIFF = 1
      KEND1 = KDIFF + NVIRT*NRHFT
      LWRK1 = LWORK - KEND1 + 1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient memory in '//SECNAM)
      ENDIF
C
      DO I = 1,NRHFT
         DO A = 1,NVIRT
            KAI = KDIFF + NVIRT*(I - 1) + A - 1
            WORK(KAI) = FOCKV(A) - FOCKO(I)
         ENDDO
      ENDDO

      write(6,*) 'Occupied orbital energies: '
      write(6,'(8F10.6)') (FOCKO(I), I=1,NRHFT)
      write(6,*) 'Virtual orbital energies: '
      write(6,'(8F10.6)') (FOCKV(I), I=1,NVIRT)
      write(6,*) 'Differences:'
      call output(WORK(KDIFF),1,NVIRT,1,NRHFT,NVIRT,NRHFT,1,6)
      call chkint(XIAJK,XAIBC,NVIRT,NRHFT)

C
C------------
C     A term.
C------------
C
      CALL CCHO_ATERM_DBG(T2AM,XAIBC,WORK(KDIFF),NVIRT,NRHFT,E4AT)

      write(6,*) ' after a...'
      call flshfo(6)
C
C------------
C     B term.
C------------
C
      CALL CCHO_BTERM_DBG(FOCKO,FOCKV,T2AM,XAIBC,NVIRT,NRHFT,E4BT)

      write(6,*) ' after b...'
      call flshfo(6)
C
C------------
C     C term.
C------------
C
C      CALL CCHO_CTERM_DBG(T2AM,XIAJK,WORK(KDIFF),NVIRT,NRHFT,E4C)
      CALL CCHO_CTERM_DBG(FOCKO,FOCKV,T2AM,XIAJK,NVIRT,NRHFT,E4C)

      write(6,*) ' after c...'
      call flshfo(6)
C
C------------
C     D term.
C------------
C
      CALL CCHO_DTERM_DBG(T2AM,XIAJK,XAIBC,WORK(KDIFF),NVIRT,NRHFT,E4D)

      write(6,*) ' after d...'
      call flshfo(6)
C
C------------
C     E term.
C------------
C
      CALL CCHO_ETERM_DBG(FOCKO,FOCKV,T2AM,XIAJK,XIAJB,XAIBC,NRHFT,
     &                    NVIRT,WORK(KEND1),LWRK1,E4E1,E4E2)

      write(6,*) ' after e...'
      call flshfo(6)
C
C------------
C     F term.
C------------
C
      CALL CCHO_FTERM_DBG(FOCKO,FOCKV,T2AM,XIAJK,XIAJB,XAIBC,NRHFT,
     &                    NVIRT,WORK(KEND1),LWRK1,E4F1,E4F2)

      write(6,*) ' after f...'
      call flshfo(6)
C
C------------
C     G term.
C------------
C
      CALL CCHO_GTERM_DBG(FOCKO,FOCKV,T1AM,T2AM,XIAJK,XIAJB,XAIBC,NRHFT,
     &                    NVIRT,WORK(KEND1),LWRK1,E5G1,E5G)

      write(6,*) ' after g...'
      call flshfo(6)
C
C------------
C     H term.
C------------
C
      CALL CCHO_HTERM_DBG(FOCKO,FOCKV,T1AM,T2AM,XIAJK,XIAJB,XAIBC,NRHFT,
     &                    NVIRT,WORK(KEND1),LWRK1,E5H1,E5H)

      write(6,*) ' after h...'
      call flshfo(6)
C
C------------
C     I term.
C------------
C
      CALL CCHO_ITERM_DBG(T2AM,XAIBC,WORK(KDIFF),NVIRT,NRHFT,E4I)

      write(6,*) ' after i...'
      call flshfo(6)
C
C------------
C     J term.
C------------
C
      CALL CCHO_JTERM_DBG(T2AM,XIAJK,XAIBC,WORK(KDIFF),NVIRT,NRHFT,E4J)

      write(6,*) ' after j...'
      call flshfo(6)

      write(6,*) 'End of CCHO_DBGDRV'
      LIAJK  = NRHFT*NRHFT*NRHFT*NVIRT
      LAI    = NVIRT*NRHFT
      LAIBC  = NVIRT*NRHFT*NVIRT*NVIRT
      XNIAJK = DSQRT(DDOT(LIAJK,XIAJK,1,XIAJK,1))
      XNT1AM = DSQRT(DDOT(LAI,T1AM,1,T1AM,1))
      XNT2AM = 0.0D0
      DO NBJ = 1,LAI
         DO NAI = 1,NBJ
            KOFF   = LAI*(NBJ - 1) + NAI
            XNT2AM = XNT2AM + T2AM(KOFF)*T2AM(KOFF)
         ENDDO
      ENDDO
      XNT2AM = DSQRT(XNT2AM)
      XNIAJB = 0.0D0
      DO NBJ = 1,LAI
         DO NAI = 1,NBJ
            KOFF   = LAI*(NBJ - 1) + NAI
            XNIAJB = XNIAJB + XIAJB(KOFF)*XIAJB(KOFF)
         ENDDO
      ENDDO
      XNIAJB = DSQRT(XNIAJB)
      XNAIBC = DSQRT(DDOT(LAIBC,XAIBC,1,XAIBC,1))
      write(6,*) 'Norm of (ia|jk): ',XNIAJK
      write(6,*) 'Norm of t1am   : ',XNT1AM
      write(6,*) 'Norm of t2am   : ',XNT2AM
      write(6,*) 'Norm of (ia|jb): ',XNIAJB
      write(6,*) 'Norm of (ai|bc): ',XNAIBC
C
C-----------
C     Print.
C-----------
C
      CALL HEADER('CCSD(T) Debug: Energy Corrections',-1)
C
      E4TH = E4AT + E4BT + E4C  + E4D + E4E1 + E4E2
     &     + E4F1 + E4F2 + E4I  + E4J
      E5TH = E5G1 + E5G  + E5H  + E5H1
      ETOT = E4TH + E5TH
C
      WRITE(LUPRI,'(10X,A,D20.10)') 'A   : ',E4AT*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'B   : ',E4BT*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'C   : ',E4C*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'D   : ',E4D*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'E1  : ',E4E1*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'E2  : ',E4E2*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'F1  : ',E4F1*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'F2  : ',E4F2*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'G   : ',E5G*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'G1  : ',E5G1*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'H   : ',E5H*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'H1  : ',E5H1*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'I   : ',E4I*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'J   : ',E4J*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'E4th: ',E4TH*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'E5th: ',E5TH*TWO
      WRITE(LUPRI,'(10X,A,D20.10)') 'ETOT: ',ETOT*TWO
C
      RETURN
      END
C  /* Deck ccho_aterm_dbg */
      SUBROUTINE CCHO_ATERM_DBG(T2AM,XAIBC,DIFF,NVIRT,NRHFT,E4AT)
#include <implicit.h>
      DIMENSION T2AM(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION DIFF(NVIRT,NRHFT)
C
      INTEGER A,B,C,D,E
      PARAMETER (ZERO = 0.00D0, TWO = 2.00D0)
C
      E4AT = ZERO
C
      DO A = 1,NVIRT
         DO B = 1,NVIRT
            DO C = 1,NVIRT
               DO D = 1,NVIRT
                  DO E = 1,NVIRT

                     DO I = 1,NRHFT
                        DO J = 1,NRHFT
                           DO L = 1,NRHFT

         E4AT = E4AT
     &        - (T2AM(C,I,A,J) * (TWO*T2AM(C,I,B,J)-T2AM(C,J,B,I)))
     &         *(XAIBC(D,L,A,E)* (TWO*XAIBC(D,L,E,B)-XAIBC(E,L,D,B)))
     &         /(DIFF(C,I) + DIFF(D,J) + DIFF(E,L))

                           ENDDO
                        ENDDO
                     ENDDO

                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      write(6,*) 'E4AT: ',E4AT
C
      RETURN
      END
C  /* Deck ccho_bterm_dbg */
      SUBROUTINE CCHO_BTERM_DBG(FO,FV,T2AM,XAIBC,NVIRT,NRHFT,E4BT)
#include <implicit.h>
      DIMENSION T2AM(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION FO(NRHFT), FV(NVIRT)
C
      INTEGER A,B,C,D,E
      PARAMETER (ZERO = 0.00D0, TWO = 2.00D0)
C
      E4BT = ZERO
C
      DO A = 1,NVIRT
         DO B = 1,NVIRT
            DO C = 1,NVIRT
               DO D = 1,NVIRT

                  DO E = 1,NVIRT
                     DO I = 1,NRHFT
                        DO J = 1,NRHFT
                           DO M = 1,NRHFT

         E4BT = E4BT
     &        - (T2AM(B,J,C,I) * (TWO*T2AM(D,J,A,I)-T2AM(D,I,A,J)))
     &         *(XAIBC(E,M,C,A)* (TWO*XAIBC(E,M,B,D)-XAIBC(B,M,E,D))
     &          -XAIBC(A,M,C,E)* XAIBC(E,M,B,D))
     &         /( FV(A) + FV(B) + FV(E)
     &          - FO(I) - FO(J) - FO(M) )

         E4BT = E4BT
     &        + (T2AM(B,J,C,I) * (TWO*T2AM(D,J,A,I)-T2AM(D,I,A,J)))
     &         *(XAIBC(D,M,C,E)* XAIBC(B,M,E,A))
     &         /( FV(D) + FV(B) + FV(E)
     &          - FO(I) - FO(J) - FO(M) )

                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO

      write(6,*) 'E4BT without DIFF: ',E4BT
C
      RETURN
      END
C  /* Deck ccho_cterm_dbg */
      SUBROUTINE CCHO_CTERM_DBG(FOCKO,FOCKV,T2AM,XIAJK,NVIRT,NRHFT,E4CT)
#include <implicit.h>
      DIMENSION T2AM(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XIAJK(NRHFT,NRHFT,NRHFT,NVIRT)
      DIMENSION FOCKO(NRHFT),FOCKV(NVIRT)
C
      INTEGER A,B,C,D,E
      PARAMETER (ZERO = 0.00D0, TWO = 2.00D0)
C
      E4CT = ZERO
C
      DO I = 1,NRHFT
         DO J = 1,NRHFT
            DO K = 1,NRHFT
               DO L = 1,NRHFT

                  DO M = 1,NRHFT
                     DO A = 1,NVIRT
                        DO B = 1,NVIRT
                           DO C = 1,NVIRT

         E4CT = E4CT
     &        - (T2AM(A,K,B,I) * (TWO*T2AM(A,K,B,J)-T2AM(A,J,B,K)))
     &         *(XIAJK(L,I,M,C)* (TWO*XIAJK(L,M,J,C)-XIAJK(M,L,J,C)))
     &         /( FOCKV(A) + FOCKV(B) + FOCKV(C)
     &          - FOCKO(M) - FOCKO(L) - FOCKO(K) )

                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO

                  DO M = 1,NRHFT
                     DO A = 1,NVIRT
                        DO B = 1,NVIRT
                           DO E = 1,NVIRT

         E4CT = E4CT
     &        - (T2AM(A,I,B,J) * (TWO*T2AM(A,K,B,L) - T2AM(A,L,B,K)))
     &         *(XIAJK(M,J,L,E)* (TWO*XIAJK(M,I,K,E)-XIAJK(I,M,K,E))
     &          -XIAJK(L,J,M,E)* XIAJK(M,I,K,E))
     &         /( FOCKV(A) + FOCKV(B) + FOCKV(E)
     &          - FOCKO(I) - FOCKO(L) - FOCKO(M) )

         E4CT = E4CT
     &        + (T2AM(A,I,B,J) * (TWO*T2AM(A,K,B,L) - T2AM(A,L,B,K)))
     &         *(XIAJK(K,J,M,E)* XIAJK(I,M,L,E))
     &         /( FOCKV(A) + FOCKV(B) + FOCKV(E)
     &          - FOCKO(I) - FOCKO(K) - FOCKO(M) )

                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO

      write(6,*) 'E4CT without DIFF(,): ',E4CT
C
      RETURN
      END
C  /* Deck ccho_dterm_dbg */
      SUBROUTINE CCHO_DTERM_DBG(T2AM,XIAJK,XAIBC,DIFF,NVIRT,NRHFT,E4DT)
#include <implicit.h>
      DIMENSION T2AM(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XIAJK(NRHFT,NRHFT,NRHFT,NVIRT)
      DIMENSION XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION DIFF(NVIRT,NRHFT)
C
      INTEGER A,B,C,D
      PARAMETER (ZERO = 0.00D0, TWO = 2.00D0)
C
      E4DT = ZERO
C
      DO A = 1,NVIRT
         DO B = 1,NVIRT
            DO C = 1,NVIRT
               DO D = 1,NVIRT

                  DO I = 1,NRHFT
                     DO J = 1,NRHFT
                        DO K = 1,NRHFT
                           DO L = 1,NRHFT

         E4DT = E4DT
     &        + ((TWO*T2AM(B,I,C,J)-T2AM(B,J,C,I))*XAIBC(C,J,B,A))
     &         *((TWO*T2AM(A,K,D,L)-T2AM(A,L,D,K))*XIAJK(L,I,K,D))
     &         /(DIFF(B,J) + DIFF(C,K) + DIFF(D,L))

                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      E4DT = TWO*E4DT

      write(6,*) 'E4DT: ',E4DT
C
      RETURN
      END
C  /* Deck ccho_iterm_dbg */
      SUBROUTINE CCHO_ITERM_DBG(T2AM,XAIBC,DIFF,NVIRT,NRHFT,E4IT)
#include <implicit.h>
      DIMENSION T2AM(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION DIFF(NVIRT,NRHFT)
C
      INTEGER A,B,C,D,E
      PARAMETER (ZERO = 0.00D0, TWO = 2.00D0)
C
      E4IT = ZERO
C
      DO A = 1,NVIRT
         DO B = 1,NVIRT
            DO C = 1,NVIRT
               DO D = 1,NVIRT
                  DO E = 1,NVIRT

                     DO I = 1,NRHFT
                        DO J = 1,NRHFT
                           DO K = 1,NRHFT

         E4IT = E4IT
     &        - (T2AM(A,I,C,K) * (TWO*T2AM(B,J,C,K)-T2AM(B,K,C,J)))
     &         *(XAIBC(D,J,A,E)*(TWO*XAIBC(E,I,D,B)-XAIBC(D,I,E,B)))
     &         /(DIFF(D,I) + DIFF(E,J) + DIFF(C,K))

         E4IT = E4IT
     &        + (T2AM(A,K,C,I)*(TWO*T2AM(B,J,C,K)-T2AM(B,K,C,J)))
     &         *(XAIBC(D,J,A,E)*XAIBC(E,I,D,B))
     &         /(DIFF(D,I) + DIFF(E,J) + DIFF(C,K))

         E4IT = E4IT
     &        + (T2AM(A,K,C,I)*(TWO*T2AM(B,K,C,J)-T2AM(B,J,C,K)))
     &         *(XAIBC(E,J,A,D)*XAIBC(E,I,D,B))
     &         /(DIFF(D,I) + DIFF(E,J) + DIFF(C,K))

                           ENDDO
                        ENDDO
                     ENDDO

                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      write(6,*) 'E4IT: ',E4IT
C
      RETURN
      END
C  /* Deck ccho_jterm_dbg */
      SUBROUTINE CCHO_JTERM_DBG(T2AM,XIAJK,XAIBC,DIFF,NVIRT,NRHFT,E4JT)
#include <implicit.h>
      DIMENSION T2AM(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XIAJK(NRHFT,NRHFT,NRHFT,NVIRT)
      DIMENSION XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION DIFF(NVIRT,NRHFT)
C
      INTEGER A,B,C,D,E
      PARAMETER (ZERO = 0.00D0, TWO = 2.00D0)
C
      E4JT = ZERO
      E4JI = ZERO
C
      DO A = 1,NVIRT
         DO I = 1,NRHFT
            DO B = 1,NVIRT
               DO J = 1,NRHFT

                  DO C = 1,NVIRT
                     DO K = 1,NRHFT
                        DO L = 1,NRHFT
                           DO M = 1,NRHFT

         E4JT = E4JT
     &        - (TWO*T2AM(A,I,C,K) - T2AM(A,K,C,I))
     &         *(TWO*T2AM(B,J,C,K) - T2AM(B,K,C,J))
     &         *(XIAJK(L,I,M,B)*XIAJK(M,L,J,A))
     &         /(DIFF(A,M) + DIFF(B,L) + DIFF(C,K))

         E4JT = E4JT
     &        + (T2AM(A,I,C,K)*(TWO*T2AM(B,J,C,K)-T2AM(B,K,C,J))
     &          +T2AM(A,K,C,I)*(TWO*T2AM(B,K,C,J)-T2AM(B,J,C,K)))
     &         *(XIAJK(L,I,M,B)*XIAJK(L,M,J,A))
     &         /(DIFF(A,M) + DIFF(B,L) + DIFF(C,K))

                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO

                  DO C = 1,NVIRT
                     DO K = 1,NRHFT
                        DO D = 1,NVIRT
                            DO L = 1,NRHFT

         E4JI = E4JI
     &        + (TWO*T2AM(A,I,C,K) - T2AM(A,K,C,I))
     &         *(TWO*T2AM(B,J,C,K) - T2AM(B,K,C,J))
     &         *(XIAJK(J,I,L,D)*XAIBC(A,L,D,B))
     &         /(DIFF(A,J) + DIFF(C,K) + DIFF(D,L))

         E4JI = E4JI
     &        + (T2AM(A,I,C,K)*(TWO*T2AM(B,J,C,K)-T2AM(B,K,C,J))
     &          +T2AM(A,K,C,I)*(TWO*T2AM(B,K,C,J)-T2AM(B,J,C,K)))
     &         *(XIAJK(L,I,J,D)*(TWO*XAIBC(D,L,A,B)-XAIBC(A,L,D,B))
     &          -XIAJK(J,I,L,D)*XAIBC(D,L,A,B))
     &         /(DIFF(A,J) + DIFF(C,K) + DIFF(D,L))

                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      E4JT = E4JT + TWO*E4JI

      write(6,*) 'E4JT: ',E4JT
C
      RETURN
      END
C  /* ETERM */
C  /* Deck ccho_eterm_dbg */
      SUBROUTINE CCHO_ETERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                          NVIRT,WORK,LWORK,E1,E2)
C
      IMPLICIT NONE
C
      INTEGER NRHFT,NVIRT,LWORK
      DOUBLE PRECISION FOCC(NRHFT),FVIR(NVIRT),
     &               T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT), 
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION E1,E2
C     
      CALL CCHO_E1TERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                     NVIRT,WORK,LWORK,E1)
      CALL CCHO_E2TERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                     NVIRT,WORK,LWORK,E2)
C     
      END
C
C  /* Deck ccho_e1term_dbg */
      SUBROUTINE CCHO_E1TERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                           NVIRT,WORK,LWORK,E_E1)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      INTEGER I,J,K,L,M,A,B,C,D,E
C      
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      NV = NVIRT
      NO = NRHFT
C      
      E_E11 = ZERO
C      
      DO A = 1,NV
         DO I = 1,NO
            DO B = 1,NV
               DO C = 1,NV
                  DO E = 1,NV
                     DO M = 1,NO
                        DO D = 1,NV
                           DO L = 1,NO
         E_E11 = E_E11
     &         + ( -XAIBC(E,M,B,A)*( TWO*T2AM(E,M,C,I)-T2AM(E,I,C,M) )
     &           - XAIBC(B,M,E,A)*( TWO*T2AM(E,I,C,M)-T2AM(E,M,C,I) ) )
     &         * ( -T2AM(A,L,D,I)*XAIBC(D,L,C,B)
     &           + T2AM(A,I,D,L)*( TWO*XAIBC(D,L,C,B)-XAIBC(B,L,C,D) ) )
     &         / ( FVIR(E) + FVIR(B) + FVIR(D)
     &           - FOCC(M) - FOCC(L) - FOCC(I) )
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      E_E12 = ZERO
C     
      DO A = 1,NV
         DO I = 1,NO
            DO B = 1,NV
               DO C = 1,NV
                  DO E = 1,NV
                     DO M = 1,NO
                        DO D = 1,NV
                           DO L = 1,NO
         E_E12 = E_E12
     &         + ( TWO*XAIBC(E,M,B,A)-XAIBC(B,M,E,A) )
     &         * ( TWO*T2AM(E,M,C,I)-T2AM(E,I,C,M) )  
     &         * ( -T2AM(A,L,D,I)*XAIBC(B,L,C,D) )
     &         / ( FVIR(E) + FVIR(B) + FVIR(D)
     &           - FOCC(M) - FOCC(L) - FOCC(I) )
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      E_E1 = E_E11 + E_E12
C     
      END
C  /* Deck ccho_e2term_dbg */
      SUBROUTINE CCHO_E2TERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                           NVIRT,WORK,LWORK,E_E2)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      INTEGER I,J,K,L,M,A,B,C,D,E
C      
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      NV = NVIRT
      NO = NRHFT
C      
      E_E21 = ZERO
C      
      DO A = 1,NV
         DO I = 1,NO
            DO B = 1,NV
               DO C = 1,NV
                  DO J = 1,NO
                     DO K = 1,NO
                        DO D = 1,NV
                           DO L = 1,NO
         E_E21 = E_E21
     &         + XIJKA(J,I,K,B)*( TWO*T2AM(A,K,C,J)-T2AM(A,J,C,K) ) 
     &         * ( -T2AM(A,L,D,I)*XAIBC(D,L,C,B)
     &           + T2AM(A,I,D,L)*( TWO*XAIBC(D,L,C,B)-XAIBC(B,L,C,D) ) )
     &         / ( FVIR(B) + FVIR(D) + FVIR(A)
     &           - FOCC(J) - FOCC(K) - FOCC(L) )
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      E_E22 = ZERO
C     
      DO A = 1,NV
         DO I = 1,NO
            DO B = 1,NV
               DO C = 1,NV
                  DO J = 1,NO
                     DO K = 1,NO
                        DO D = 1,NV
                           DO L = 1,NO
         E_E22 = E_E22
     &         + ( XIJKA(K,I,J,B)*( TWO*T2AM(A,K,C,J)-T2AM(A,J,C,K) ) )
     &         * ( -T2AM(A,L,D,I)*XAIBC(B,L,C,D) )
     &         / ( FVIR(B) + FVIR(D) + FVIR(A)
     &           - FOCC(J) - FOCC(K) - FOCC(L) )
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      E_E2 = TWO*E_E21 + TWO*E_E22
C     
      END
C  /* FTERM */
C  /* Deck ccho_fterm_dbg */
      SUBROUTINE CCHO_FTERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                          NVIRT,WORK,LWORK,EF1,EF2)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      DOUBLE PRECISION EF1,EF2
C     
      CALL CCHO_F1TERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                     NVIRT,WORK,LWORK,EF1)
      CALL CCHO_F2TERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                     NVIRT,WORK,LWORK,EF2)
C     
      END
C      
C  /* Deck ccho_f1term_dbg */
      SUBROUTINE CCHO_F1TERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                           NVIRT,WORK,LWORK,E_F1)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      INTEGER I,J,K,L,M,A,B,C,D,E
C      
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C      
      E_F11 = ZERO
C      
      DO A = 1,NVIRT
         DO I = 1,NRHFT
            DO J = 1,NRHFT
               DO K = 1,NRHFT
                  DO E = 1,NVIRT
                     DO M = 1,NRHFT
                        DO D = 1,NVIRT
                           DO L = 1,NRHFT
         E_F11 = E_F11
     &         + ( -XIJKA(M,I,J,E)*( TWO*T2AM(E,M,A,K)-T2AM(E,K,A,M) )
     &           -  XIJKA(J,I,M,E)*( TWO*T2AM(E,K,A,M)-T2AM(E,M,A,K) ) )
     &         * ( - T2AM(A,L,D,I)*XIJKA(L,J,K,D)
     &           + T2AM(A,I,D,L)*( TWO*XIJKA(L,J,K,D)-XIJKA(J,L,K,D) ) )
     &         / ( FVIR(E) + FVIR(D) + FVIR(A)
     &           - FOCC(M) - FOCC(J) - FOCC(L) )
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      E_F12 = ZERO
C
      DO A = 1,NVIRT
         DO I = 1,NRHFT
            DO J = 1,NRHFT
               DO K = 1,NRHFT
                  DO E = 1,NVIRT
                     DO M = 1,NRHFT
                        DO D = 1,NVIRT
                           DO L = 1,NRHFT
         E_F12 = E_F12
     &         + ( TWO*XIJKA(M,I,J,E)-XIJKA(J,I,M,E) )
     &         * ( TWO*T2AM(E,M,A,K)-T2AM(E,K,A,M) ) 
     &         * ( - T2AM(A,L,D,I)*XIJKA(J,L,K,D) )
     &         / ( FVIR(E) + FVIR(D) + FVIR(A)
     &            - FOCC(M) - FOCC(J) - FOCC(L) )   
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      E_F1 = E_F11 + E_F12

      write(6,*) 'E_F1: ',E_F1
C      
      END
C  /* Deck ccho_f2term_dbg */
      SUBROUTINE CCHO_F2TERM_DBG(FOCC,FVIR,T2AM,XIJKA,XAIBJ,XAIBC,NRHFT,
     &                           NVIRT,WORK,LWORK,E_F2)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      INTEGER I,J,K,L,M,A,B,C,D,E
C      
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C      
      E_F21 = ZERO
C      
      DO A = 1,NVIRT
         DO I = 1,NRHFT
            DO J = 1,NRHFT
               DO K = 1,NRHFT
                  DO B = 1,NVIRT
                     DO C = 1,NVIRT
                        DO D = 1,NVIRT
                           DO L = 1,NRHFT
         E_F21 = E_F21
     &         + XAIBC(B,J,C,A)*( TWO*T2AM(B,K,C,I)-T2AM(B,I,C,K) ) 
     &         * ( -T2AM(A,L,D,I)*XIJKA(L,J,K,D)
     &           + T2AM(A,I,D,L)*( TWO*XIJKA(L,J,K,D)-XIJKA(J,L,K,D) ) )
     &         / ( FVIR(B) + FVIR(C) + FVIR(D)
     &            -FOCC(J) - FOCC(L) - FOCC(I) )
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      E_F22 = ZERO
C     
      DO A = 1,NVIRT
         DO I = 1,NRHFT
            DO J = 1,NRHFT
               DO K = 1,NRHFT
                  DO B = 1,NVIRT
                     DO C = 1,NVIRT
                        DO D = 1,NVIRT
                           DO L = 1,NRHFT
         E_F22 = E_F22
     &         + XAIBC(C,J,B,A)*( TWO*T2AM(B,K,C,I)-T2AM(B,I,C,K) ) 
     &         * ( -T2AM(A,L,D,I)*XIJKA(J,L,K,D) )
     &         / ( FVIR(B) + FVIR(C) + FVIR(D)
     &            -FOCC(J) - FOCC(L) - FOCC(I) )
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      E_F2 = TWO*(E_F21 + E_F22)

      write(6,*) 'E_F2: ',E_F2
C      
      END
C  /* GTERM */
C  /* Deck ccho_gterm_dbg */
      SUBROUTINE CCHO_GTERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                          NRHFT,NVIRT,WORK,LWORK,EG1,EG2)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T1AM(NVIRT,NRHFT),
     &               T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      DOUBLE PRECISION EG1, EG2
C
      write(6,*) 'In GTERM:'
      write(6,*) 'Occupied orbital energies:'
      write(6,'(8F10.6)') (FOCC(I),I=1,NRHFT)
      write(6,*) 'Virtual  orbital energies:'
      write(6,'(8F10.6)') (FVIR(I),I=1,NVIRT)
C     
      CALL CCHO_G1TERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                     NRHFT,NVIRT,WORK,LWORK,EG1)
      CALL CCHO_G2TERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                     NRHFT,NVIRT,WORK,LWORK,EG2)
C
C Further debugging of X and R intermediates:
C
c     KRMAT = 1
c     KXMAT = KRMAT + NVIRT*NRHFT*NRHFT*NRHFT
c     KEND  = KXMAT + NVIRT*NRHFT*NRHFT*NRHFT
c     LWRK  = LWORK - KEND + 1
c     IF (LWRK.LT.0) CALL QUIT('Insufficient memory in CCHO_GTERM_DBG')
c     CALL CCHO_RXMAT_DBG(T1AM,T2AM,XAIBJ,XAIBC,NVIRT,NRHFT,
c    &                    WORK(KRMAT),WORK(KXMAT))
C
      END
C  /* Deck ccho_rxmat_dbg */
      SUBROUTINE CCHO_RXMAT_DBG(T1AM,T2AM,XAIBJ,XAIBC,NVIRT,NRHFT,
     &                          RMAT,XMAT)
C
#include <implicit.h>
      DIMENSION T1AM(NVIRT,NRHFT)
      DIMENSION T2AM(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION RMAT(NVIRT,NRHFT,NRHFT,NRHFT)
      DIMENSION XMAT(NVIRT,NRHFT,NRHFT,NRHFT)

      INTEGER   A, B, C, D

      PARAMETER (TWO = 2.00D0)


      LENMAT = NVIRT*NRHFT*NRHFT*NRHFT
      CALL DZERO(XMAT,LENMAT)
      CALL DZERO(RMAT,LENMAT)

C     Calculate X(b,j,i,k).
C     ---------------------

      DO K = 1,NRHFT
         DO I = 1,NRHFT
            DO J = 1,NRHFT
               DO B = 1,NVIRT

                  DO A = 1,NVIRT
                     DO D = 1,NVIRT

                        XMAT(B,J,I,K) = XMAT(B,J,I,K)
     &                  + XAIBC(D,J,B,A)
     &                   *(TWO*XAIBJ(A,I,D,K) - XAIBJ(A,K,D,I))

                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     Calculate R(b,i,k,j).
C     ---------------------

      DO J = 1,NRHFT
         DO K = 1,NRHFT
            DO I = 1,NRHFT
               DO B = 1,NVIRT

                  DO C = 1,NVIRT

                     RMAT(B,I,K,J) = RMAT(B,I,K,J)
     &               + T2AM(B,I,C,K)*T1AM(C,J)

                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     Get the product of X and R.
C     ---------------------------

      XR = 0.00D0

      DO J = 1,NRHFT
         DO K = 1,NRHFT
            DO I = 1,NRHFT
               DO B = 1,NVIRT

                  XR = XR + XMAT(B,J,I,K)*RMAT(B,I,K,J)

               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     Print norms.
C     ------------

      XNORM = DNORM2(LENMAT,XMAT,1)
      RNORM = DNORM2(LENMAT,RMAT,1)

      WRITE(6,*)
      WRITE(6,*) '   CCHO_RXMAT_DBG:'
      WRITE(6,*) '   ==============='
      WRITE(6,*) '   NOTICE: Unscaled quantitites!'
      WRITE(6,*) '   Norm of XMAT: ',XNORM
      WRITE(6,*) '   Norm of RMAT: ',RNORM
      WRITE(6,*) '   Product X*R : ',XR
      WRITE(6,*)
      WRITE(6,*) '   XMAT(bj,ik):'
      NBJ = NVIRT*NRHFT
      NIK = NRHFT*NRHFT
      CALL OUTPUT(XMAT,1,NBJ,1,NIK,NBJ,NIK,1,6)
      WRITE(6,*)
      WRITE(6,*) '   RMAT(bik,j):'
      NBIK = NVIRT*NRHFT*NRHFT
      NJ   = NRHFT
      CALL OUTPUT(RMAT,1,NBIK,1,NJ,NBIK,NJ,1,6)
      WRITE(6,*)

      RETURN
      END
C  /* Deck ccho_g1term_dbg */
      SUBROUTINE CCHO_G1TERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                           NRHFT,NVIRT,WORK,LWORK,EG1)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T1AM(NVIRT,NRHFT),
     &               T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      INTEGER I,J,K,L,M,A,B,C,D,E
C
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      NV = NVIRT
      NO = NRHFT
C     
      EG1 = ZERO
C     
      DO A = 1,NV
         DO B = 1,NV
            DO C = 1,NV
               DO I = 1,NO
                  DO J = 1,NO
                     DO D = 1,NV
                        DO L = 1,NO
         EG1 = EG1
     &       - T2AM(A,I,C,J)*( TWO*XAIBJ(C,J,B,I)-XAIBJ(C,I,B,J) )
     &       * ( TWO*XAIBC(D,L,A,B)-XAIBC(B,L,A,D) )*T1AM(D,L) 
     &       / ( FVIR(C) + FVIR(D) + FVIR(B)
     &         - FOCC(I) - FOCC(J) - FOCC(L) )

                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      END
C  /* Deck ccho_g2term_dbg */
      SUBROUTINE CCHO_G2TERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                           NRHFT,NVIRT,WORK,LWORK,EG2)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T1AM(NVIRT,NRHFT),
     &               T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      INTEGER I,J,K,L,M,A,B,C,D,E
C      
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      NV = NVIRT
      NO = NRHFT
C     
      EG21 = ZERO
C     
      DO A = 1,NV
         DO I = 1,NO
            DO B = 1,NV
               DO J = 1,NO
                  DO C = 1,NV
                     DO K = 1,NO
                        DO D = 1,NV
          EG21 = EG21
     &         - ( TWO*T2AM(A,I,B,J)-T2AM(A,J,B,I) )*T1AM(B,J)
     &         * ( XAIBC(C,K,A,D)*(TWO*XAIBJ(C,K,D,I)-XAIBJ(C,I,D,K)) )
     &         / ( FVIR(B) + FVIR(C) + FVIR(D)
     &           - FOCC(J) - FOCC(I) - FOCC(K) )

                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      EG22 = ZERO
C     
      DO A = 1,NV
         DO I = 1,NO
            DO J = 1,NO
               DO K = 1,NO
                  DO B = 1,NV
                     DO C = 1,NV
                        DO D = 1,NV
          EG22 = EG22
     &         + T2AM(A,I,B,J)*T1AM(B,K)
     &         * XAIBC(C,K,A,D)*( TWO*XAIBJ(D,I,C,J)-XAIBJ(D,J,C,I) )
     &         / ( FVIR(B) + FVIR(C) + FVIR(D)
     &           - FOCC(J) - FOCC(I) - FOCC(K) )
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      EG23 = ZERO
C     
      DO A = 1,NV
         DO I = 1,NO
            DO B = 1,NV
               DO J = 1,NO
                  DO C = 1,NV
                     DO K = 1,NO
                        DO D = 1,NV
          EG23 = EG23
     &         - ( TWO*T2AM(C,K,A,I)-T2AM(C,I,A,K)   )
     &         * ( TWO*XAIBJ(C,K,B,J)-XAIBJ(C,J,B,K) ) 
     &         * XAIBC(B,J,A,D)*T1AM(D,I) 
     &         / ( FVIR(C) + FVIR(B) + FVIR(D)
     &           - FOCC(K) - FOCC(I) - FOCC(J) )
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      EG24 = ZERO
C     
      DO A = 1,NV
         DO I = 1,NO
            DO B = 1,NV
               DO J = 1,NO
                  DO D = 1,NV
                     DO C = 1,NV
                        DO K = 1,NO
          EG24 = EG24
     &         + XAIBC(D,J,A,B)*T1AM(D,I) 
     &         * ( T2AM(A,I,C,K)
     &           * ( TWO*XAIBJ(C,K,B,J)-XAIBJ(B,K,C,J) )
     &           + T2AM(A,K,C,I)
     &           * ( TWO*XAIBJ(B,K,C,J)-XAIBJ(C,K,B,J) )
     &           )
     &         / ( FVIR(B) + FVIR(D) + FVIR(C)
     &           - FOCC(J) - FOCC(K) - FOCC(I) )
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      EG2 = EG21 + EG22 + EG23 + EG24

c     write(6,*)
c     write(6,*) '   End of GTERM_DBG:'
c     write(6,*) '   EG21: ',EG21,' (Y*S)'
c     write(6,*) '   EG22: ',EG22,' (X*R)'
c     write(6,*) '   EG23: ',EG23,' (W*Q)'
c     write(6,*) '   EG24: ',EG24,' (V*P)'
c     write(6,*)
C     
      END
C  /* HTERM */
C   /* Deck ccho_hterm_dbg */
      SUBROUTINE CCHO_HTERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                          NRHFT,NVIRT,WORK,LWORK,EH1,EH2)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T1AM(NVIRT,NRHFT),
     &               T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      DOUBLE PRECISION EH1, EH2
C     
      CALL CCHO_H1TERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                     NRHFT,NVIRT,WORK,LWORK,EH1)
      CALL CCHO_H2TERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                     NRHFT,NVIRT,WORK,LWORK,EH2)
C
c     WRITE(6,*)
c     WRITE(6,*) '   HTERM_DBG:'
c     WRITE(6,*) '   - calculating V,W,P,Q matrices:'
c     WRITE(6,*)
c     KVMAT = 1
c     KWMAT = KVMAT + NVIRT*NRHFT*NRHFT*NRHFT
c     KPMAT = KWMAT + NVIRT*NRHFT*NRHFT*NRHFT
c     KQMAT = KPMAT + NVIRT*NRHFT*NRHFT*NRHFT
c     KEND  = KQMAT + NVIRT*NRHFT*NRHFT*NRHFT
c     LWRK  = LWORK - KEND + 1
c     IF (LWRK. LT. 0) THEN
c        WRITE(6,*) '   - insufficient memory for test. Skipping...'
c        WRITE(6,*)
c     ELSE
c        CALL H_VWPQ(T1AM,T2AM,XIJKA,XAIBJ,WORK(KVMAT),WORK(KWMAT),
c    &               WORK(KPMAT),WORK(KQMAT),NVIRT,NRHFT)
c     ENDIF
C     
      END
C  /* Deck h_vw_2 */
      SUBROUTINE H_VW_2(XK1,XK2,XM1,XM2,VMAT,WMAT,NVIRT,NRHFT)
C
#include <implicit.h>
      DIMENSION XK1(NVIRT,NRHFT,NRHFT,NRHFT)
      DIMENSION XK2(NVIRT,NRHFT,NRHFT,NRHFT)
      DIMENSION XM1(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XM2(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION VMAT(NVIRT,NRHFT,NRHFT,NRHFT)
      DIMENSION WMAT(NVIRT,NRHFT,NRHFT,NRHFT)

      INTEGER A,B,D

      PARAMETER (TWO = 2.00D0)

      NO2   = NRHFT*NRHFT
      NOV   = NVIRT*NRHFT
      NO3V  = NO2*NOV
      NO2V2 = NOV*NOV

      WRITE(6,*)
      WRITE(6,*) '   Entering H_VW_2:'
      WRITE(6,*) '   ================'
      WRITE(6,*) '   Norm of K1: ',DNORM2(NO3V,XK1,1)
      WRITE(6,*) '   Norm of K2: ',DNORM2(NO3V,XK2,1)
      WRITE(6,*) '   Norm of M1: ',DNORM2(NO2V2,XM1,1)
      WRITE(6,*) '   Norm of M2: ',DNORM2(NO2V2,XM2,1)
      WRITE(6,*)

C     VMAT Coul.
C     ----------

      CALL DZERO(VMAT,NO3V)
      DO K = 1,NRHFT
         DO J = 1,NRHFT
            DO I = 1,NRHFT
               DO A = 1,NVIRT

                  DO L = 1,NRHFT
                     DO D = 1,NVIRT

                        VMAT(A,I,J,K) = VMAT(A,I,J,K)
     &                  + XM1(D,L,A,I)*XK1(D,L,J,K)

                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*) '   Norm of VMAT after Coul.: ',
     &               DNORM2(NO3V,VMAT,1)

C     WMAT Coul.
C     ----------

      CALL DCOPY(NO3V,VMAT,1,WMAT,1)
      WRITE(6,*) '   Norm of WMAT after Coul.: ',
     &               DNORM2(NO3V,WMAT,1)

C     VMAT Exch.
C     ----------

      CALL DSCAL(NO3V,TWO,VMAT,1)
      DO K = 1,NRHFT
         DO J = 1,NRHFT
            DO I = 1,NRHFT
               DO A = 1,NVIRT

                  DO L = 1,NRHFT
                     DO D = 1,NVIRT

                        VMAT(A,I,J,K) = VMAT(A,I,J,K)
     &                  - XM1(D,L,A,I)*XK2(D,L,J,K)

                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*) '   Norm of VMAT after Exch.: ',
     &               DNORM2(NO3V,VMAT,1)

C     WMAT Exch.
C     ----------

      DO K = 1,NRHFT
         DO J = 1,NRHFT
            DO I = 1,NRHFT
               DO A = 1,NVIRT

                  DO L = 1,NRHFT
                     DO D = 1,NVIRT

                        WMAT(A,I,J,K) = WMAT(A,I,J,K)
     &                  + XM2(D,L,A,I)*XK2(D,L,J,K)

                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*) '   Norm of WMAT after Exch.: ',
     &               DNORM2(NO3V,WMAT,1)
      WRITE(6,*)

      NAI = NVIRT*NRHFT
      NJK = NRHFT*NRHFT

      WRITE(6,*) '   VMAT(ai,jk):'
      WRITE(6,*) '   Norm of V: ',DNORM2(NAI*NJK,VMAT,1)
      CALL OUTPUT(VMAT,1,NAI,1,NJK,NAI,NJK,1,6)
      WRITE(6,*)

      WRITE(6,*) '   WMAT(ai,jk):'
      WRITE(6,*) '   Norm of W: ',DNORM2(NAI*NJK,WMAT,1)
      CALL OUTPUT(WMAT,1,NAI,1,NJK,NAI,NJK,1,6)
      WRITE(6,*)

      WRITE(6,*)
      WRITE(6,*) '   Exiting H_VW_2:'
      WRITE(6,*) '   ==============='
      WRITE(6,*) '   Norm of K1: ',DNORM2(NO3V,XK1,1)
      WRITE(6,*) '   Norm of K2: ',DNORM2(NO3V,XK2,1)
      WRITE(6,*) '   Norm of M1: ',DNORM2(NO2V2,XM1,1)
      WRITE(6,*) '   Norm of M2: ',DNORM2(NO2V2,XM2,1)
      WRITE(6,*)

      RETURN
      END
C  /* Deck h_vwpq */
      SUBROUTINE H_VWPQ(T1AM,T2AM,XIJKA,XAIBJ,VMAT,WMAT,PMAT,QMAT,
     &                  NVIRT,NRHFT)
C
#include <implicit.h>
      DIMENSION T1AM(NVIRT,NRHFT)
      DIMENSION T2AM(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION XIJKA(NRHFT,NRHFT,NRHFT,NVIRT)
      DIMENSION XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT)
      DIMENSION VMAT(NVIRT,NRHFT,NRHFT,NRHFT)
      DIMENSION WMAT(NVIRT,NRHFT,NRHFT,NRHFT)
      DIMENSION PMAT(NVIRT,NRHFT,NRHFT,NRHFT)
      DIMENSION QMAT(NVIRT,NRHFT,NRHFT,NRHFT)
C
      INTEGER A,B,D
C
      PARAMETER (TWO = 2.00D0)

      NO2   = NRHFT*NRHFT
      NOV   = NVIRT*NRHFT
      NO3V  = NO2*NOV
      NO2V2 = NOV*NOV

      XT1NRM = DNORM2(NOV,T1AM,1)
      XT2NRM = DNORM2(NO2V2,T2AM,1)
      XNIJKA = DNORM2(NO3V,XIJKA,1)
      XNAIBJ = 0.00D0
      DO J = 1,NRHFT
         DO B = 1,NVIRT
            DO I = 1,NRHFT
               DO A = 1,NVIRT
                  XNAIBJ = XNAIBJ
     &            + (TWO*XAIBJ(A,I,B,J)-XAIBJ(A,J,B,I))**2
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      XNAIBJ = DSQRT(XNAIBJ)

      WRITE(6,*) '   H_VWPQ: Norm of T1AM   : ',XT1NRM
      WRITE(6,*) '   H_VWPQ: Norm of T2AM   : ',XT2NRM
      WRITE(6,*) '   H_VWPQ: Norm of (ia|jk): ',XNIJKA
      WRITE(6,*) '   H_VWPQ: Norm of L(iajb): ',XNAIBJ
      WRITE(6,*)

C     VMAT Coul.
C     ----------

      CALL DZERO(VMAT,NO3V)
      DO K = 1,NRHFT
         DO J = 1,NRHFT
            DO I = 1,NRHFT
               DO A = 1,NVIRT

                  DO L = 1,NRHFT
                     DO D = 1,NVIRT

                        VMAT(A,I,J,K) = VMAT(A,I,J,K)
     &                  + XIJKA(L,K,J,D)
     &                   *(TWO*XAIBJ(D,L,A,I) - XAIBJ(D,I,A,L))

                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*) '   Norm of VMAT after Coul.: ',
     &               DNORM2(NO3V,VMAT,1)

C     WMAT Coul.
C     ----------

      CALL DCOPY(NO3V,VMAT,1,WMAT,1)
      WRITE(6,*) '   Norm of WMAT after Coul.: ',
     &               DNORM2(NO3V,WMAT,1)

C     VMAT Exch.
C     ----------

      CALL DSCAL(NO3V,TWO,VMAT,1)
      DO K = 1,NRHFT
         DO J = 1,NRHFT
            DO I = 1,NRHFT
               DO A = 1,NVIRT

                  DO L = 1,NRHFT
                     DO D = 1,NVIRT

                        VMAT(A,I,J,K) = VMAT(A,I,J,K)
     &                  - XIJKA(J,K,L,D)
     &                   *(TWO*XAIBJ(D,L,A,I) - XAIBJ(D,I,A,L))

                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*) '   Norm of VMAT after Exch.: ',
     &               DNORM2(NO3V,VMAT,1)

C     WMAT Exch.
C     ----------

      DO K = 1,NRHFT
         DO J = 1,NRHFT
            DO I = 1,NRHFT
               DO A = 1,NVIRT

                  DO L = 1,NRHFT
                     DO D = 1,NVIRT

                        WMAT(A,I,J,K) = WMAT(A,I,J,K)
     &                  + XIJKA(J,K,L,D)
     &                   *(TWO*XAIBJ(D,I,A,L) - XAIBJ(D,L,A,I))

                     ENDDO
                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*) '   Norm of WMAT after Exch.: ',
     &               DNORM2(NO3V,WMAT,1)
      WRITE(6,*)

C     PMAT
C     ----

      CALL DZERO(PMAT,NO3V)
      DO J = 1,NRHFT
         DO K = 1,NRHFT
            DO I = 1,NRHFT
               DO A = 1,NVIRT

                  DO B = 1,NVIRT

                     PMAT(A,I,K,J) = PMAT(A,I,K,J)
     &               + T2AM(A,I,B,K)*T1AM(B,J)

                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO

C     QMAT
C     ----

      CALL DZERO(QMAT,NO3V)
      DO J = 1,NRHFT
         DO K = 1,NRHFT
            DO I = 1,NRHFT
               DO A = 1,NVIRT

                  DO B = 1,NVIRT

                     QMAT(A,I,K,J) = QMAT(A,I,K,J)
     &               - T2AM(A,K,B,I)*T1AM(B,J)

                  ENDDO

               ENDDO
            ENDDO
         ENDDO
      ENDDO

      NAI = NVIRT*NRHFT
      NJK = NRHFT*NRHFT

      WRITE(6,*) '   VMAT(ai,jk):'
      WRITE(6,*) '   Norm of V: ',DNORM2(NAI*NJK,VMAT,1)
      CALL OUTPUT(VMAT,1,NAI,1,NJK,NAI,NJK,1,6)
      WRITE(6,*)

      WRITE(6,*) '   WMAT(ai,jk):'
      WRITE(6,*) '   Norm of W: ',DNORM2(NAI*NJK,WMAT,1)
      CALL OUTPUT(WMAT,1,NAI,1,NJK,NAI,NJK,1,6)
      WRITE(6,*)

      NAIK = NAI*NRHFT
      NJ   = NRHFT

      WRITE(6,*) '   PMAT(aik,j):'
      WRITE(6,*) '   Norm of P: ',DNORM2(NAIK*NJ,PMAT,1)
      CALL OUTPUT(PMAT,1,NAIK,1,NJ,NAIK,NJ,1,6)
      WRITE(6,*)

      WRITE(6,*) '   QMAT(aik,j):'
      WRITE(6,*) '   Norm of Q: ',DNORM2(NAIK*NJ,QMAT,1)
      CALL OUTPUT(QMAT,1,NAIK,1,NJ,NAIK,NJ,1,6)
      WRITE(6,*)

      RETURN
      END
C  /* Deck ccho_h1term_dbg */
      SUBROUTINE CCHO_H1TERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                           NRHFT,NVIRT,WORK,LWORK,EH1)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T1AM(NVIRT,NRHFT),
     &               T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      INTEGER I,J,K,L,M,A,B,C,D,E
C      
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      NV = NVIRT
      NO = NRHFT
C     
      EH1 = ZERO
C     
      DO A = 1,NV
         DO I = 1,NO
            DO B = 1,NV
               DO J = 1,NO
                  DO C = 1,NV
                     DO K = 1,NO
                        DO L = 1,NO
          EH1 = EH1 
     &        + ( TWO*T2AM(A,I,B,J)-T2AM(A,J,B,I) )*T1AM(B,J)
     &        * XIJKA(K,I,L,C)*( TWO*XAIBJ(C,K,A,L)-XAIBJ(A,K,C,L) )
     &        / ( FVIR(A) + FVIR(B) + FVIR(C)
     &          - FOCC(J) - FOCC(K) - FOCC(L) )
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      END
C  /* Deck ccho_h2term_dbg */
      SUBROUTINE CCHO_H2TERM_DBG(FOCC,FVIR,T1AM,T2AM,XIJKA,XAIBJ,XAIBC,
     &                           NRHFT,NVIRT,WORK,LWORK,EH2)
C
C
#include <implicit.h>
C
      DIMENSION FOCC(NRHFT),FVIR(NVIRT),T1AM(NVIRT,NRHFT),
     &               T2AM(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XIJKA(NRHFT,NRHFT,NRHFT,NVIRT),
     &               XAIBJ(NVIRT,NRHFT,NVIRT,NRHFT),
     &               XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
      DIMENSION WORK(LWORK)
C     
      INTEGER I,J,K,L,M,A,B,C,D,E
C      
      PARAMETER (XMONE = -1.0D0, ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      NV = NVIRT
      NO = NRHFT
C     
      EH21 = ZERO
C     
      DO I = 1,NO
         DO J = 1,NO
            DO A = 1,NV
               DO B = 1,NV
                  DO K = 1,NO
                     DO D = 1,NV
                        DO L = 1,NO
          EH21 = EH21
     &         + T2AM(A,I,B,K)*( TWO*XAIBJ(B,K,A,J)-XAIBJ(A,K,B,J) )
     &         * ( TWO*XIJKA(L,I,J,D)-XIJKA(J,I,L,D) )*T1AM(D,L)
     &         / ( FVIR(A) + FVIR(B) + FVIR(D)
     &           - FOCC(K) - FOCC(L) - FOCC(J) )
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      EH22 = ZERO
C     
      DO I = 1,NO
         DO J = 1,NO
            DO K = 1,NO
               DO L = 1,NO
                  DO A = 1,NV
                     DO B = 1,NV
                        DO C = 1,NV
          EH22 = EH22
     &         - T2AM(A,I,B,K)*( TWO*XAIBJ(A,J,B,L)-XAIBJ(B,J,A,L) )
     &         * XIJKA(L,I,J,C)*T1AM(C,K)
     &         / ( FVIR(A) + FVIR(B) + FVIR(C)
     &           - FOCC(K) - FOCC(L) - FOCC(J) )
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      EH23 = ZERO
C     
      DO A = 1,NV
         DO I = 1,NO
            DO J = 1,NO
               DO K = 1,NO
                  DO B = 1,NV
                     DO D = 1,NV
                        DO L = 1,NO
          EH23 = EH23
     &         + T2AM(A,I,B,J)*T1AM(B,K)
     &         * ( TWO*XIJKA(L,J,K,D)-XIJKA(K,J,L,D) )
     &         * ( TWO*XAIBJ(D,L,A,I)-XAIBJ(A,L,D,I) )
     &         / ( FVIR(A) + FVIR(B) + FVIR(D)
     &           - FOCC(I) - FOCC(L) - FOCC(K) )
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      EH24 = ZERO
C     
      DO A = 1,NV
         DO I = 1,NO
            DO J = 1,NO
               DO K = 1,NO
                  DO B = 1,NV
                     DO D = 1,NV
                        DO L = 1,NO
          EH24 = EH24
     &         - T2AM(A,J,B,I)*T1AM(B,K)
     &         * (
     &            XIJKA(L,J,K,D)*( TWO*XAIBJ(A,I,D,L)-XAIBJ(A,L,D,I) )
     &           + XIJKA(K,J,L,D)*( TWO*XAIBJ(A,L,D,I)-XAIBJ(A,I,D,L) )
     &           )
     &         / ( FVIR(A) + FVIR(B) + FVIR(D)
     &           - FOCC(I) - FOCC(L) - FOCC(K) )
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C     
      EH2 = EH21 + EH22 + EH23 + EH24
C
      write(6,*)
      write(6,*)
      write(6,*) '   H TERMS:'
      write(6,*) '   ========'
      write(6,*) '   EVP = ',EH23
      write(6,*) '   EWQ = ',EH24
      write(6,*) '   EXR = ',EH22
      write(6,*) '   EYS = ',EH21
      write(6,*)
      write(6,*)
C     
      END
C  /* Deck chkint */
      SUBROUTINE CHKINT(XIAJK,XAIBC,NVIRT,NRHFT)
C
#include <implicit.h>
      DIMENSION XIAJK(NRHFT,NRHFT,NRHFT,NVIRT)
      DIMENSION XAIBC(NVIRT,NRHFT,NVIRT,NVIRT)
C
      INTEGER A,B,C
C
      ICOUNO = 0
      ICOUNV = 0
C
      DO A = 1,NVIRT
         DO K = 1,NRHFT
            DO J = K+1,NRHFT
               DO I = 1,NRHFT
                  TEST = XIAJK(I,J,K,A) - XIAJK(I,K,J,A)
                  IF (DABS(TEST) .GT. 1.0D-14) THEN
                     WRITE(6,*) 'For I,J,K,A = ',I,J,K,A,':'
                     WRITE(6,*) 'int,diff: ',XIAJK(I,J,K,A),TEST
                     ICOUNO = ICOUNO + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      IF (ICOUNO .GT. 0) WRITE(6,*)
C
      DO B = 1,NVIRT
         DO C = B+1,NVIRT
            DO I = 1,NRHFT
               DO A = 1,NVIRT
                  TEST = XAIBC(A,I,C,B) - XAIBC(A,I,B,C)
                  IF (DABS(TEST) .GT. 1.0D-14) THEN
                     WRITE(6,*) 'For A,I,B,C = ',A,I,B,C,':'
                     WRITE(6,*) 'int,diff: ',XAIBC(A,I,C,B),TEST
                     ICOUNV = ICOUNV + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      IF (ICOUNV .GT. 0) WRITE(6,*)
C
      WRITE(6,*) 'Occupied non-symmetries: ',ICOUNO
      WRITE(6,*) 'Virtual  non-symmetries: ',ICOUNV
      WRITE(6,*)
C
      END
C  /* Deck chk_vec1 */
      SUBROUTINE CHK_VEC1(FOCKD,CHOELE,NUMCHO,ICHO,
     &                    CHOL,CHOO,CHOV,ISYMMN,IA1,NUMIA,ISYMA,
     &                    J1,NUMIJ,ISYMJ,WORK,LWORK)
C
#include <implicit.h>
      DIMENSION FOCKD(*), CHOELE(*)
      DIMENSION CHOL(*), CHOO(*), CHOV(*)
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
      INTEGER IOFF2(8)

      write(6,*) 'CHK_VEC1: LWORK = ',LWORK
      write(6,*) 'CHK_VEC1: IA1,NUMIA,ISYMA: ',IA1,NUMIA,ISYMA
      write(6,*) 'CHK_VEC1: J1,NUMIJ,ISYMJ : ',J1,NUMIJ,ISYMJ
      write(6,*) 'CHK_VEC1: ISYMMN         : ',ISYMMN

      ICOUNT = 0
      DO ISYM = 1,NSYM
         IOFF2(ISYM) = ICOUNT
         ICOUNT = ICOUNT + NT1AM(ISYM)
      ENDDO
      NTOT1 = ICOUNT

      KCHOO = 1
      KCHOV = KCHOO + NTOT1*NUMIJ
      KEND1 = KCHOV + NTOT1*NUMIA
      LWRK1 = LWORK - KEND1 + 1

      IF (LWRK1 .LT. 0) THEN
         WRITE(6,*) 'Insufficient memory for Cholesky test in CHK_VEC1'
         WRITE(6,*) 'Available: ',LWORK
         WRITE(6,*) 'Need     : ',KEND1-1
         CALL QUIT('Insufficient memory in CHK_VEC1')
      ENDIF

      CALL CCHO_DECHO(FOCKD,CHOELE,NUMCHO,ICHO,WORK(KCHOO),WORK(KCHOV),
     &                IA1,NUMIA,ISYMA,J1,NUMIJ,ISYMJ)

      ICOUNV = 0
      DO A = 1,NUMIA
         DO ISYMCK = 1,NSYM
            DO ISYMK = 1,NSYM
               ISYMC = MULD2H(ISYMK,ISYMCK)
               DO K = 1,NRHF(ISYMK)
                  DO C = 1,NVIR(ISYMC)
                     KOFF1 = IOFF2(ISYMCK)
     &                     + NT1AM(ISYMCK)*(A - 1)
     &                     + IT1AM(ISYMC,ISYMK)
     &                     + NVIR(ISYMC)*(K - 1)
     &                     + C
                     KOFF2 = KCHOV + KOFF1 - 1
                     TEST  = WORK(KOFF2) - CHOV(KOFF1)
                     IF (DABS(TEST) .GT. 1.0D-15) THEN
                        ICOUNV = ICOUNV + 1
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*)
      WRITE(6,*) '   CHK_VEC1: Virtual Part'
      WRITE(6,*) '   ======================'
      WRITE(6,*) '  ',ICOUNV,' errors out of ',NTOT1*NUMIA
      WRITE(6,*) '   IA1,NUMIA,ISYMA: ',IA1,NUMIA,ISYMA

      IF (ICOUNV .GT. 0) THEN
         WRITE(6,*) '   NOTICE: Replacing virtual part with',
     &              ' decho vectors!'
         CALL DCOPY(NTOT1*NUMIA,WORK(KCHOV),1,CHOV,1)
      ENDIF

      ICOUNO = 0
      DO J = 1,NUMIJ
         DO ISYMCK = 1,NSYM
            DO ISYMK = 1,NSYM
               ISYMC = MULD2H(ISYMK,ISYMCK)
               DO K = 1,NRHF(ISYMK)
                  DO C = 1,NVIR(ISYMC)
                     KOFF1 = IOFF2(ISYMCK)
     &                     + NT1AM(ISYMCK)*(J - 1)
     &                     + IT1AM(ISYMC,ISYMK)
     &                     + NVIR(ISYMC)*(K - 1)
     &                     + C
                     KOFF2 = KCHOO + KOFF1 - 1
                     TEST  = WORK(KOFF2) - CHOO(KOFF1)
                     IF (DABS(TEST) .GT. 1.0D-15) THEN
                        ICOUNO = ICOUNO + 1
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*)
      WRITE(6,*) '   CHK_VEC1: Occupied Part'
      WRITE(6,*) '   ======================='
      WRITE(6,*) '  ',ICOUNO,' errors out of ',NTOT1*NUMIJ
      WRITE(6,*) '   J1,NUMIJ,ISYMJ: ',J1,NUMIJ,ISYMJ

      IF (ICOUNO .GT. 0) THEN
         WRITE(6,*) '   NOTICE: Replacing occupied part with',
     &              ' decho vectors!'
         CALL DCOPY(NTOT1*NUMIJ,WORK(KCHOO),1,CHOO,1)
      ENDIF

      WRITE(6,*)

      RETURN
      END
C  /* Deck chk_vec2 */
      SUBROUTINE CHK_VEC2(FOCKD,CHOELE,NUMCHO,ICHO,CHOL,OCCHO,VICHO,
     &                    ISYDJA,IB1,NUMIB,ISYMB,II1,NUMII,ISYMI,WORK,
     &                    LWORK)
C
#include <implicit.h>
      DIMENSION FOCKD(*), CHOELE(*)
      DIMENSION CHOL(*), OCCHO(*), VICHO(*)
      DIMENSION WORK(LWORK)
#include <ccorb.h>
#include <ccsdsym.h>
      INTEGER IOFT1(8)

      WRITE(6,*)
      WRITE(6,*) 'Entering CHK_VEC2'
      WRITE(6,*) 'ICHO,NUMCHO    : ',ICHO,NUMCHO
      WRITE(6,*) 'IB1,NUMIB,ISYMB: ',IB1,NUMIB,ISYMB
      WRITE(6,*) 'II1,NUMII,ISYMI: ',II1,NUMII,ISYMI

C     Set up IOFT1.
C     -------------

      ICOUNT = 0
      DO ISYM = 1,NSYM
         IOFT1(ISYM) = ICOUNT
         ICOUNT = ICOUNT + NT1AM(ISYM)
      ENDDO
      NTOT1 = ICOUNT

C     Check consistency of CHOL and VICHO.
C     ------------------------------------

      ICOUNE = 0
      ISYMDJ = MULD2H(ISYMB,ISYDJA)
      DO B = 1,NUMIB
         IB = IB1 + B - 1
         DO NDJ = 1,NT1AM(ISYMDJ)
            KOFFC = ICKATR(ISYMDJ,ISYMB)
     &            + NT1AM(ISYMDJ)*(IB - 1) + NDJ
            KOFFV = NTOT1*(B - 1)
     &            + IOFT1(ISYMDJ) + NDJ
            TEST  = CHOL(KOFFC) - VICHO(KOFFV)
            IF (DABS(TEST) .GT. 1.0D-15) THEN
               WRITE(6,*) 'CHK_VEC2: Error: ISYMDJ,DJ,B: ',
     &                    ISYMDJ,NDJ,B,' CHOL - VICHO: ',TEST
               ICOUNE = ICOUNE + 1
            ENDIF
         ENDDO
      ENDDO
      WRITE(6,*) 'CHK_VEC2: ',ICOUNE,' errors comparing CHOL and VICHO'

C     Check OCCHO.
C     ------------

      MAXCKI = -1
      DO ISYM = 1,NSYM
         MAXCKI = MAX(MAXCKI,NCKI(ISYM))
      ENDDO

      KCHOT = 1
      KCHOO = KCHOT + NTOT1*NUMII
      KEND  = KCHOO + MAXCKI
      LWRK  = LWORK - KEND + 1

      IF (LWRK .LT. 0) CALL QUIT(' crappy test in CHK_VEC2')

      DUMMY = 0.0D0
      IDUM1 = -1
      IDUM2 = 0
      DO ISYCKI = 1,NSYM
         CALL CCHO_DECHO3(FOCKD,CHOELE,NUMCHO,ICHO,WORK(KCHOO),DUMMY,
     &                    ISYCKI,IDUM1,1,IDUM2)
         ISYMCK = MULD2H(ISYCKI,ISYMI)
         DO I = 1,NUMII
            II = II1 + I - 1
            KOFF1 = KCHOO + ICKI(ISYMCK,ISYMI) + NT1AM(ISYMCK)*(II - 1)
            KOFF2 = KCHOT + NTOT1*(I - 1) + IOFT1(ISYMCK)
            CALL DCOPY(NT1AM(ISYMCK),WORK(KOFF1),1,WORK(KOFF2),1)
         ENDDO
      ENDDO

      ICOUNE = 0
      DO I = 1,NUMII
         DO ISYMCK = 1,NSYM
            DO NCK = 1,NT1AM(ISYMCK)
               KORIG = NTOT1*(I - 1) + IOFT1(ISYMCK) + NCK
               KOFF1 = KCHOT + KORIG - 1
               TEST  = OCCHO(KORIG) - WORK(KOFF1)
               IF (DABS(TEST) .GT. 1.0D-15) THEN
                  WRITE(6,*) 'CHK_VEC2: Error: ISYMCK,CK,I: ',
     &                       ISYMCK,NCK,II1+I-1,
     &                       ' OCCHO - DECHO3: ',TEST
                  ICOUNE = ICOUNE + 1                  
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*) 
     & 'CHK_VEC2: ',ICOUNE,' errors comparing OCCHO and DECHO3'

C     Check VICHO.
C     ------------

      KCHOO = 1
      KCHOV = KCHOO + NTOT1*NUMII
      KEND  = KCHOV + NTOT1*NUMIB
      LWRK  = LWORK - KEND + 1

      IF (LWRK .LT. 0) CALL QUIT(' crappy test in CHK_VEC2 (virt.)')

      CALL CCHO_DECHO(FOCKD,CHOELE,NUMCHO,ICHO,WORK(KCHOO),WORK(KCHOV),
     &                IB1,NUMIB,ISYMB,II1,NUMII,ISYMI)

      ICOUNE = 0
      DO B = 1,NUMIB
         DO ISYMDJ = 1,NSYM
            DO NDJ = 1,NT1AM(ISYMDJ)
               IADR = NTOT1*(B - 1) + IOFT1(ISYMDJ) + NDJ
               KOFF = KCHOV + IADR - 1
               TEST = VICHO(IADR) - WORK(KOFF)
               IF (DABS(TEST) .GT. 1.0D-15) THEN
                  WRITE(6,*) 'CHK_VEC2: Error: ISYMDJ,DJ,B: ',
     &                       ISYMDJ,NDJ,IB1+B-1,
     &                       ' VICHO - DECHO : ',TEST
                  ICOUNE = ICOUNE + 1                  
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*)
     & 'CHK_VEC2: ',ICOUNE,' errors comparing VICHO and DECHO'

      ICOUNE = 0
      DO I = 1,NUMII
         DO ISYMCK = 1,NSYM
            DO NCK = 1,NT1AM(ISYMCK)
               KORIG = NTOT1*(I - 1) + IOFT1(ISYMCK) + NCK
               KOFF1 = KCHOO + KORIG - 1
               TEST  = OCCHO(KORIG) - WORK(KOFF1)
               IF (DABS(TEST) .GT. 1.0D-15) THEN
                  WRITE(6,*) 'CHK_VEC2: Error: ISYMCK,CK,I: ',
     &                       ISYMCK,NCK,II1+I-1,
     &                       ' OCCHO - DECHO : ',TEST
                  ICOUNE = ICOUNE + 1                  
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*) 
     & 'CHK_VEC2: ',ICOUNE,' errors comparing OCCHO and DECHO'
      WRITE(6,*)


      RETURN
      END
C  /* Deck chk_gt */
      SUBROUTINE CHK_GT(T1,T4,II1,NUMII,ISYMI,IB1,NUMIB,ISYMB,
     &                  IOFF5,IOFF6)
C
#include <implicit.h>
      DIMENSION T1(*), T4(*)
      INTEGER IOFF5(8), IOFF6(8)
#include <ccorb.h>
#include <ccsdsym.h>

      ISYMCK = MULD2H(ISYMB,ISYMI)
      ICOUNE = 0
      NTST   = 0

      DO I = 1,NUMII
         DO B = 1,NUMIB
            DO ISYMK = 1,NSYM
               ISYMC = MULD2H(ISYMK,ISYMCK)
               DO K = 1,NRHF(ISYMK)
                  DO C = 1,NVIR(ISYMC)
                     KOFF1 = IOFF5(ISYMB)
     &                     + NT1AM(ISYMCK)*NUMIB*(I - 1)
     &                     + NT1AM(ISYMCK)*(B - 1)
     &                     + IT1AM(ISYMC,ISYMK)
     &                     + NVIR(ISYMC)*(K - 1)
     &                     + C
                     KOFF4 = IOFF6(ISYMC)
     &                     + NUMIB*NUMII*NRHF(ISYMK)*(C - 1)
     &                     + NUMIB*NUMII*(K - 1)
     &                     + NUMIB*(I - 1)
     &                     + B
                     TEST  = T1(KOFF1) - T4(KOFF4)
                     NTST  = NTST + 1
                     IF (DABS(TEST) .GT. 1.0D-15) THEN
                        WRITE(6,*) 'CHK_GT: ISYMC,ISYMK,ISYMB,ISYMI: ',
     &                                      ISYMC,ISYMK,ISYMB,ISYMI
                        WRITE(6,*) 'CHK_GT: C,K,B,I: ',
     &                                      C,K,IB1+B-1,II1+I-1
                        WRITE(6,*) 'CHK_GT: Diff.  : ',TEST
                        ICOUNE = ICOUNE + 1
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      NDIM = NT1AM(ISYMCK)*NUMIB*NUMII
      WRITE(6,*) 'CHK_GT: ',ICOUNE,' errors comparing T1 and T4',
     &           ' out of ',NTST,' tested (dim = ',NDIM,')'
      WRITE(6,*)

      RETURN
      END
C  /* Deck chk_hk12 */
      SUBROUTINE CHK_HK12(XK1,XK2,IK1,NUMIK,ISYMK,IOFK12)
C
#include <implicit.h>
      DIMENSION XK1(*), XK2(*)
      INTEGER IOFK12(8)
#include <ccorb.h>
#include <ccsdsym.h>

      ICOUNE = 0
      NTST   = 0
      ISYDLJ = ISYMK
      DO ISYMJ = 1,NSYM
         ISYMDL = MULD2H(ISYMJ,ISYDLJ)
         DO K = 1,NUMIK
            DO J = 1,NRHF(ISYMJ)
               DO ISYML = 1,NSYM
                  ISYMD = MULD2H(ISYML,ISYMDL)
                  ISYMDJ = MULD2H(ISYMD,ISYMJ)
                  DO L = 1,NRHF(ISYML)
                     DO D = 1,NVIR(ISYMD)
                        KDLJK = IOFK12(ISYMJ)
     &                        + NT1AM(ISYMDL)*NRHF(ISYMJ)*(K - 1)
     &                        + NT1AM(ISYMDL)*(J - 1)
     &                        + IT1AM(ISYMD,ISYML)
     &                        + NVIR(ISYMD)*(L - 1) + D
                        KDJLK = IOFK12(ISYML)
     &                        + NT1AM(ISYMDJ)*NRHF(ISYML)*(K - 1)
     &                        + NT1AM(ISYMDJ)*(L - 1)
     &                        + IT1AM(ISYMD,ISYMJ)
     &                        + NVIR(ISYMD)*(J - 1) + D
                        NTST  = NTST + 1
                        TEST  = XK1(KDLJK) - XK2(KDJLK)
                        IF (DABS(TEST) .GT. 1.0D-15) THEN
      write(6,*) 'CHK_HK12: ISYMD,ISYML,ISYMJ,ISYMK: ',
     &                         ISYMD,ISYML,ISYMJ,ISYMK
      write(6,*) 'CHK_HK12: D,L,J,K: ',D,L,J,IK1+K-1
      write(6,*) 'CHK_HK12: Diff.  : ',TEST
                           ICOUNE = ICOUNE + 1
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      ndim = NCKI(ISYMK)*NUMIK
      write(6,*) 'CHK_HK12: ',ICOUNE,' interchange errors out of ',
     &           NTST,' tested (dim = ',ndim,')'

      RETURN
      END
