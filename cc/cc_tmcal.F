*---------------------------------------------------------------------*
c /* deck cc_tmcal */
*=====================================================================*
       SUBROUTINE CC_TMCAL(WRK,LWRK) 
*---------------------------------------------------------------------*
*
*    Purpose: Third moment calculations
*    
*    Written by: P.Joergensen and C.Haettig 1997
*    Clean up  : S. Coriani 2000
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#include <priunit_dec.h>
#include <ccorb_dec.h>
#include <ccsdsym_dec.h>
#else
#include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <cctm.h>
#include <cctminf.h>
#include <ccrspprp.h>
#include <ccexci.h>
#include <ccroper.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB, LABELC,
     *            LABELD, LABELE, LABELF 
      CHARACTER*10 MODFIL
      INTEGER ISYMB, ISYMC, ISYMA, ISYMD, ISYME, ISYMF,  ISYMABC
      INTEGER IFREQ, INUM, IOPER, IDX, IOFFST, LWRK, IOPT
      INTEGER K1VEC1, K1VEC2, K2VEC1, K2VEC2, IM11
      INTEGER IX3AC0F, IX3DF0F, IO3AC0F, IO3ACF0, IO3DF0F, IO3DFF0
      INTEGER NCCVAR1, NCCVAR2

#if defined (SYS_CRAY)
      REAL HALF, FREQEX, FREQB, FREQC, EIGV, WRK(LWRK)
      REAL SMLM, SMCLM, SMRM, SMCRM
      REAL X1, X2, Y1, Y2
      REAL DDOT, ZERO
#else
      DOUBLE PRECISION HALF, FREQEX, FREQB, FREQC, EIGV, WRK(LWRK) 
      DOUBLE PRECISION SMLM, SMCLM, SMRM, SMCRM
      DOUBLE PRECISION X1, X2, Y1, Y2
      DOUBLE PRECISION DDOT, ZERO
#endif

      PARAMETER ( HALF = 0.5D00, ZERO = 0.0D00 )

* external functions:
      INTEGER IRHSR3
      INTEGER ILRMAMP
      INTEGER ICHI3

* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_TMIND> NTMOP = ',NTMOPER
      END IF

*---------------------------------------------------------------------*
* find list entries for the required response vectors
* and excitation vectors:
*---------------------------------------------------------------------*

      DO IOPER = 1, NTMOPER
        LABELA = LBLOPR(IATMOP(IOPER))
        LABELB = LBLOPR(IBTMOP(IOPER))
        LABELC = LBLOPR(ICTMOP(IOPER))
        LABELD = LBLOPR(IDTMOP(IOPER))
        LABELE = LBLOPR(IETMOP(IOPER))
        LABELF = LBLOPR(IFTMOP(IOPER))

        ISYMA  = ISYOPR(IATMOP(IOPER))
        ISYMB  = ISYOPR(IBTMOP(IOPER))
        ISYMC  = ISYOPR(ICTMOP(IOPER))
        ISYMD  = ISYOPR(IDTMOP(IOPER))
        ISYME  = ISYOPR(IETMOP(IOPER))
        ISYMF  = ISYOPR(IFTMOP(IOPER))

        WRITE (LUPRI,*) 'LABELA:',LABELA
        WRITE (LUPRI,*) 'LABELB:',LABELB
        WRITE (LUPRI,*) 'LABELC:',LABELC
        WRITE (LUPRI,*) 'LABELD:',LABELD
        WRITE (LUPRI,*) 'LABELE:',LABELE
        WRITE (LUPRI,*) 'LABELF:',LABELF

        WRITE (LUPRI,*) 'ISYMA:', ISYMA
        WRITE (LUPRI,*) 'ISYMB:', ISYMB
        WRITE (LUPRI,*) 'ISYMC:', ISYMC
        WRITE (LUPRI,*) 'ISYMD:', ISYMD
        WRITE (LUPRI,*) 'ISYME:', ISYME
        WRITE (LUPRI,*) 'ISYMF:', ISYMF

        CALL FLSHFO(LUPRI)
        
        ISYMABC = MULD2H(MULD2H(ISYMA,ISYMB),ISYMC)
        IF (ISYMABC .EQ. MULD2H(ISYMD,MULD2H(ISYME,ISYMF))) THEN
      
          NCCVAR1 = NT1AM(ISYMABC)
          NCCVAR2 = NT2AM(ISYMABC)
          K1VEC1   = 1
          K1VEC2   = K1VEC1 + NCCVAR1
          K2VEC1   = K1VEC2 + NCCVAR2
          K2VEC2   = K2VEC1 + NCCVAR1
     
          DO I = 1, NTMSELX(ISYMABC)  
              IFREQ  = ITMSELX( MULD2H(ISYMA,ISYMB) ) + I
              FREQEX  = EXTMFR(IFREQ)
              FREQB  = BTMFR(IFREQ)
              FREQC  = CTMFR(IFREQ)
              IF (LOCDBG) THEN
                WRITE (LUPRI,*) 'CC_TMCAL> put on the list:',
     &            LABELA,'(',FREQEX,'),  ', LABELB,'(',FREQB ,'),  ',
     &            LABELC,'(',FREQC, '),  ', LABELD,'(',FREQEX,'),  ',
     &            LABELE,'(',FREQB, '),  ', LABELF,'(',FREQC ,'),  '
              END IF
*           request third order chi vectors:

           IX3AC0F = ICHI3(LABELA,-FREQEX+FREQB+FREQC,ISYMA,
     &                  LABELB,-FREQB,ISYMB,LABELC,-FREQC,ISYMC)
           IX3DF0F = ICHI3(LABELD,-FREQEX+FREQB+FREQC,ISYMD,
     &                  LABELE,-FREQB,ISYME,LABELF,-FREQC,ISYMF)


*           request third order rhs vectors

           IO3AC0F = IRHSR3(LABELA,-FREQEX+FREQB+FREQC,ISYMA,
     &                  LABELB,-FREQB,ISYMB,LABELC,-FREQC,ISYMC)
           IO3ACF0 = IRHSR3(LABELA,+FREQEX-FREQB-FREQC,ISYMA,
     &                  LABELB,+FREQB,ISYMB,LABELC,+FREQC,ISYMC)
           IO3DF0F = IRHSR3(LABELD,-FREQEX+FREQB+FREQC,ISYMD,
     &                  LABELE,-FREQB,ISYME,LABELF,-FREQC,ISYMF)
           IO3DFF0 = IRHSR3(LABELD,+FREQEX-FREQB-FREQC,ISYMD,
     &                  LABELE,+FREQB,ISYME,LABELF,+FREQC,ISYMF)

*           request M vectors for different excitation energies

             IOFFST = ISYOFE(ISYMABC) +  ITMSEL(IFREQ,2)
             EIGV   = EIGVAL(IOFFST)
             IM11   = ILRMAMP(IOFFST,EIGV,ISYMABC)

*            calculate left  moment ABC-w contribution

             IOPT = 3
             CALL CC_RDRSP('X3',IX3AC0F,ISYMABC,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             X1 = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K1VEC1),1)
             IF (.NOT.CCS) THEN
               X2 = DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K1VEC2),1) 
             ELSE
               X2 = ZERO
             END IF
             WRITE (LUPRI,*) ' norm^2 of X3 vector:',X1,X2

             IOPT = 3
             CALL CC_RDRSP('RE',IOFFST,ISYMABC,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             y1 = DDOT(NCCVAR1,WRK(K2VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               Y2 = DDOT(NCCVAR2,WRK(K2VEC2),1,WRK(K2VEC2),1) 
             ELSE
               Y2 = ZERO
             END IF
             WRITE (LUPRI,*) ' norm^2 of RE vector:',Y1,Y2

             SMCLM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMCLM=SMCLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF
             WRITE (LUPRI,*)' SMCLM1',SMCLM


             IOPT = 3
             CALL CC_RDRSP('M1',IM11,ISYMABC,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             X1 = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K1VEC1),1)
             IF (.NOT.CCS) THEN
               X2 = DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K1VEC2),1) 
             ELSE
               X2 = ZERO
             END IF
             WRITE (LUPRI,*) 'NORM^2 OF M1:',X1,X2,X1+X2

             IOPT = 3
             CALL CC_RDRSP('O3',IO3AC0F,ISYMABC,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             y1 = DDOT(NCCVAR1,WRK(K2VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               y2 = DDOT(NCCVAR2,WRK(K2VEC2),1,WRK(K2VEC2),1) 
             ELSE
               y2 = ZERO
             END IF
             WRITE (LUPRI,*) ' norm^2 of O3 vector:',y1,y2,y1+y2

             CALL CCLR_DIASCL(WRK(K2VEC2),HALF,ISYMABC)

             SMCLM = SMCLM + DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMCLM=SMCLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF

             WRITE (LUPRI,*) ' smclm2',smclm
*            calculate right moment DEF+w contribution

             IOPT = 3
             CALL CC_RDRSP('LE',IOFFST,ISYMABC,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             IOPT = 3
             CALL CC_RDRSP('O3',IO3DFF0,ISYMABC,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),HALF,ISYMABC)

             SMCRM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMCRM=SMCRM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF

*            calculate left  moment DEF-w contribution

             IOPT = 3
             CALL CC_RDRSP('X3',IX3DF0F,ISYMABC,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))

             IOPT = 3
             CALL CC_RDRSP('RE',IOFFST,ISYMABC,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             SMLM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMLM = SMLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF

             IOPT = 3
             CALL CC_RDRSP('M1',IM11,ISYMABC,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))

             IOPT = 3
             CALL CC_RDRSP('O3',IO3DF0F,ISYMABC,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),HALF,ISYMABC)

             SMLM = SMLM + DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMLM = SMLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF

*            calculate right moment AB+w contribution 

             IOPT = 3
             CALL CC_RDRSP('LE',IOFFST,ISYMABC,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             IOPT = 3
             CALL CC_RDRSP('O3',IO3ACF0,ISYMABC,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),HALF,ISYMABC)

             SMRM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1) 
             IF (.NOT.CCS) THEN
               SMRM = SMRM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF

             WRITE (LUPRI,*) ' ablmom',smclm,' cdrmom',smcrm
             WRITE (LUPRI,*) ' cdlmom',smlm,' abrmom',smrm
             WRITE (LUPRI,*) ' w moment product', smclm*smcrm
             WRITE (LUPRI,*) ' -w moment product', smlm*smrm
             WRITE (LUPRI,*) ' tot moment product',
     &            0.5d0*(smclm*smcrm+smlm*smrm)

          END DO
        END IF

      END DO


      RETURN
      END
*=====================================================================*
*---------------------------------------------------------------------*

       SUBROUTINE CC_TMSORT
*---------------------------------------------------------------------*
*
*    Purpose: sort the selected states for which third moment 
*             calculation is carried. if no selected states are
*             chosen use all states specified in the excitation
*             energy calculation is used
*
*=====================================================================*

#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <cclr_dec.h>
#else
#include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <cctm.h>
#include <cctminf.h>
#include <ccexci.h>
#include <cclr.h>


* local parameters:

      INTEGER ISYM, IST, ISEL, I, ISAVE, JSEL, J, IOFF 
      INTEGER ISYMSV, ISTSV, JSTSV, ISTATE 
C

#if defined (SYS_CRAY)
      REAL D3, BTMFRSV, CTMFRSV
#else
      DOUBLE PRECISION D3, BTMFRSV, CTMFRSV
#endif
      PARAMETER ( D3 = 3.0d00 )

C
C sort the selected states for which second order transition
C matrix elements are calculated
C
      DO 50 ISYM = 1,NSYM
         NTMSELX(ISYM) = 0
 50   CONTINUE
C
      IF ( SELTMST ) THEN
C
C sort list according to symmetry
C
         ITMSELX(1) = 0
         DO 100 ISYM = 1,NSYM
            IST = ITMSELX(ISYM) + 1
            DO 200 I = IST,NTMSEL
               IF ( ITMSEL(I,1).EQ.ISYM) THEN
                  NTMSELX(ISYM) = NTMSELX(ISYM) + 1
               ELSE
                  DO 300 J = I+1,NTMSEL
                     IF ( ITMSEL(J,1).EQ.ISYM) THEN
                        ISYMSV = ITMSEL(J,1)          
                        ISTSV  = ITMSEL(J,2)          
                        BTMFRSV = BTMFR(J)
                        CTMFRSV = CTMFR(J)
                        ITMSEL(J,1) = ITMSEL(I,1)
                        ITMSEL(J,2) = ITMSEL(I,2)
                        BTMFR(J)   = BTMFR(I)
                        CTMFR(J)   = CTMFR(I)
                        ITMSEL(I,1) = ISYMSV
                        ITMSEL(I,2) = ISTSV
                        BTMFR(I)   = BTMFRSV
                        CTMFR(I)   = CTMFRSV
                        NTMSELX(ISYM) = NTMSELX(ISYM) + 1
                        GO TO 200
                     END IF
 300              CONTINUE
               END IF
 200        CONTINUE
            IF ( ISYM .LT. NSYM ) THEN
               ITMSELX(ISYM+1) = ITMSELX(ISYM) + NTMSELX(ISYM)
            END IF
            WRITE (LUPRI,*) 'SORT:',ITMSELX(ISYM),NTMSELX(ISYM),IST
 100     CONTINUE
         WRITE (LUPRI,*) ' after sort of  symmetry '
         WRITE (LUPRI,*) 'ntmsel',ntmsel
         do 210 i = 1,ntmsel
         WRITE (LUPRI,*) ' itmsel(i,1),itmsel(i,2),i'
         WRITE (LUPRI,*)  itmsel(i,1),itmsel(i,2),i
 210     continue     
         do 211 i = 1,nsym
         WRITE (LUPRI,*) ' itmselx(i),ntmselx(i),i'
         WRITE (LUPRI,*) itmselx(i),ntmselx(i),i
 211     continue
C
C sort list according to state number
C 
         DO 400 ISYM = 1,NSYM
            IOFF = ITMSELX(ISYM)
            DO 500 ISEL = 1,NTMSELX(ISYM)
               I = IOFF + ISEL
               ISTSV = ITMSEL(I,2)
               ISAVE  = I
               DO 600 JSEL = ISEL+1,NTMSELX(ISYM)
                  J = IOFF + JSEL
                  JSTSV = ITMSEL(J,2)
                  IF ( JSTSV.LT. ISTSV ) THEN 
                     ISTSV  = JSTSV
                     ISAVE  = J
                  END IF
 600           CONTINUE
               IF ( I.NE.ISAVE ) THEN
                  ISYMSV = ITMSEL(ISAVE,1)          
                  ISTSV  = ITMSEL(ISAVE,2)          
                  BTMFRSV = BTMFR(ISAVE)
                  CTMFRSV = CTMFR(ISAVE)
                  ITMSEL(ISAVE,1) = ITMSEL(I,1)
                  ITMSEL(ISAVE,2) = ITMSEL(I,2)
                  BTMFR(ISAVE)   = BTMFR(I)
                  CTMFR(ISAVE)   = CTMFR(I)
                  ITMSEL(I,1) = ISYMSV
                  ITMSEL(I,2) = ISTSV
                  BTMFR(I)   = BTMFRSV
                  CTMFR(I)   = CTMFRSV
               END IF
 500        CONTINUE
 400     CONTINUE
         WRITE (LUPRI,*) ' after sort of both symmetry and state'
         WRITE (LUPRI,*) 'ntmsel',ntmsel
         do 212 i = 1,ntmsel
         WRITE (LUPRI,*) ' itmsel(i,1),itmsel(i,2),i'
         WRITE (LUPRI,*) itmsel(i,1),itmsel(i,2),i
 212     continue     
         do 213 i = 1,nsym
         WRITE (LUPRI,*) ' itmselx(i),ntmselx(i),i'
         WRITE (LUPRI,*) itmselx(i),ntmselx(i),i
 213     continue
C
C if .HALFFR not specified find frequencies for AOPERATOR
C
         DO 550 ISYM = 1,NSYM
            IOFF = ITMSELX(ISYM)
            WRITE (LUPRI,*) 'isym,ioff'
            WRITE (LUPRI,*) isym,ioff
            DO 560 I = 1,NTMSELX(ISYM)
               ISTSV = ITMSEL(IOFF+I,2)
         WRITE (LUPRI,*) 'istsv,ioff,isym,i'
         WRITE (LUPRI,*) istsv,ioff,isym,i
         WRITE (LUPRI,*) ' isyofe(isym)'
         WRITE (LUPRI,*) isyofe(isym)
         WRITE (LUPRI,*) ' eigval(1)'
         call flshfo(LUPRI)
         WRITE (LUPRI,*) eigval(1)
         call flshfo(LUPRI)
         WRITE (LUPRI,*) ' eigval(isyofe(isym)+istsv)'
         call flshfo(LUPRI)
         WRITE (LUPRI,*) eigval(isyofe(isym)+istsv)
         call flshfo(LUPRI)
               EXTMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTSV)
         call flshfo(LUPRI)
              WRITE (LUPRI,*) ' EXTMFR(IOFF+I) '
         call flshfo(LUPRI)
              WRITE (LUPRI,*) EXTMFR(IOFF+I) 
 560        CONTINUE
         call flshfo(LUPRI)
         WRITE (LUPRI,*) ' isym loop slut',isym
 550     CONTINUE
      END IF
C
C if selected states not specified for second moment calculations
C then carry out calculations for all specified excited states
C and use frequencies that are half the excitation energy
C
      IF ( .NOT. SELTMST ) THEN
         call flshfo(LUPRI)
c        WRITE(LUPRI,*) ' jeg er 1'
         call flshfo(LUPRI)
         ITMSELX(1) = 0 
         NTMSEL = 0
         DO 700 ISYM = 1,NSYM
            DO 750 I = 1,NCCEXCI(ISYM,1)
               NTMSEL = NTMSEL + 1
               ITMSEL(NTMSEL,1) = ISYM
               ITMSEL(NTMSEL,2) = I
               NTMSELX(ISYM)    = NTMSELX(ISYM) + 1
 750        CONTINUE
         call flshfo(LUPRI)
c        WRITE(LUPRI,*) ' jeg er 2'
            ITMSELX(ISYM+1) = ITMSELX(ISYM) + NTMSELX(ISYM)
 700     CONTINUE
         THIRDFR = .TRUE.
         call flshfo(LUPRI)
c        WRITE(LUPRI,*) ' jeg er 3'
      END IF
C 
C
      IF (THIRDFR) THEN
         call flshfo(LUPRI)
c        WRITE(LUPRI,*) ' jeg er 4'
         DO 800  ISYM = 1,NSYM
            IOFF = ITMSELX(ISYM)
            DO 850 I = 1,NTMSELX(ISYM)
               ISTATE = ITMSEL(IOFF+I,2)
               BTMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTATE)/ D3
               CTMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTATE)/ D3
               EXTMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTATE)
 850        CONTINUE
         call flshfo(LUPRI)
c        WRITE(LUPRI,*) ' jeg er 5'
 800     CONTINUE
         call flshfo(LUPRI)
c        WRITE(LUPRI,*) ' jeg er 6'
      END IF 
c        call flshfo(LUPRI)
c      WRITE(LUPRI,*) ' ' leaving sort'
c        do i = 1,ntmsel
c        WRITE(LUPRI,*) ' itmsel(i,1),itmsel(i,2),extmfr(i),i'
c        call flshfo(LUPRI)
c        WRITE(LUPRI,*) itmsel(i,1),itmsel(i,2),extmfr(i),i
c        end do     
c        call flshfo(LUPRI)
                
      RETURN
      END 
*=====================================================================*
