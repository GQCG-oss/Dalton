*---------------------------------------------------------------------*
c /* deck cc_smcal */
*=====================================================================*
       SUBROUTINE CC_SMCAL(WRK,LWRK) 
*---------------------------------------------------------------------*
*
*    Purpose: Second moment calculations
*
*    Written by: Poul Joergensen and Christof Haettig  1997
*    Clean up  : Sonia Coriani 2000
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
#include <ccsm.h>
#include <ccsminf.h>
#include <ccrspprp.h>
#include <ccexci.h>
#include <ccroper.h>

* local parameters:
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

* variables:
      CHARACTER*8 LABELA, LABELB, LABELC, LABELD
      CHARACTER*10 MODFIL
      INTEGER ISYMB, ISYMC, ISYMA, ISYMD, ISYMAB
      INTEGER IFREQ, INUM, IOPER, IDX, IOFFST, LWRK, IM11, IOPT
      INTEGER IX2AB0F, IX2CD0F, IO2AB0F, IO2ABF0, IO2CD0F, IO2CDF0
      INTEGER K1VEC1, K1VEC2, K2VEC1, K2VEC2, NCCVAR1, NCCVAR2 

#if defined (SYS_CRAY)
      REAL WRK(LWRK)
      REAL ZERO, FACTOR, FREQEX, FREQB, EIGV
      REAL SMLM, SMCLM, SMRM, SMCRM
      REAL X1, X2, Y1, Y2, TESTP1, TESTP2
      REAL DDOT
#else
      DOUBLE PRECISION WRK(LWRK)
      DOUBLE PRECISION ZERO, FACTOR, FREQEX, FREQB, EIGV
      DOUBLE PRECISION SMLM, SMCLM, SMRM, SMCRM
      DOUBLE PRECISION X1, X2, Y1, Y2, TESTP1, TESTP2
      DOUBLE PRECISION DDOT
#endif

      PARAMETER ( ZERO = 0.0D00, FACTOR = 0.5D00 )


* external functions:
      INTEGER IRHSR2
      INTEGER ILRMAMP
      INTEGER ICHI2
* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./


      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_SMCAL> NSMOP = ',NSMOPER
      END IF

*---------------------------------------------------------------------*
* find list entries for the required response vectors
* and excitation vectors:
*---------------------------------------------------------------------*

      DO IOPER = 1, NSMOPER
        LABELA = LBLOPR(IASMOP(IOPER))
        LABELB = LBLOPR(IBSMOP(IOPER))
        LABELC = LBLOPR(ICSMOP(IOPER))
        LABELD = LBLOPR(IDSMOP(IOPER))

        ISYMA  = ISYOPR(IASMOP(IOPER))
        ISYMB  = ISYOPR(IBSMOP(IOPER))
        ISYMC  = ISYOPR(ICSMOP(IOPER))
        ISYMD  = ISYOPR(IDSMOP(IOPER))

        WRITE (LUPRI,*) 'LABELA:',LABELA
        WRITE (LUPRI,*) 'LABELB:',LABELB
        WRITE (LUPRI,*) 'LABELC:',LABELC
        WRITE (LUPRI,*) 'LABELD:',LABELD

        WRITE (LUPRI,*) 'ISYMA:', ISYMA
        WRITE (LUPRI,*) 'ISYMB:', ISYMB
        WRITE (LUPRI,*) 'ISYMC:', ISYMC
        WRITE (LUPRI,*) 'ISYMD:', ISYMD

        CALL FLSHFO(LUPRI)
        
        IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
      
          ISYMAB = MULD2H(ISYMA,ISYMB)
          NCCVAR1 = NT1AM(ISYMAB)
          NCCVAR2 = NT2AM(ISYMAB)
          K1VEC1  = 1
          K1VEC2  = K1VEC1 + NCCVAR1
          K2VEC1  = K1VEC2 + NCCVAR2
          K2VEC2  = K2VEC1 + NCCVAR1
     
          DO I = 1, NSMSELX(ISYMAB)  
             IFREQ  = ISMSELX( ISYMAB ) + I
             FREQEX = EXSMFR(IFREQ)
             FREQB  = BSMFR(IFREQ)
             IF (LOCDBG) THEN
                WRITE (LUPRI,*) 'CC_SMCAL> put on the list:',
     &          LABELA,'(',FREQEX,'),  ', LABELB,'(',FREQB,'),  ',
     &          LABELC,'(',FREQEX,'),  ', LABELD,'(',FREQB,'),  '
             END IF

*           request second order chi vectors:
 
        IX2AB0F = ICHI2(LABELA,.FALSE.,-FREQEX+FREQB,ISYMA,
     &                  LABELB,.FALSE.,-FREQB,ISYMB)
C       IX2CD0F = ICHI2(LABELC,.FALSE.,-FREQEX-FREQB,ISYMC,
C    &                  LABELD,.FALSE.,+FREQB,ISYMD)
        IX2CD0F = ICHI2(LABELC,.FALSE.,-FREQEX+FREQB,ISYMC,
     &                  LABELD,.FALSE.,-FREQB,ISYMD)

*           request second-order rhs vectors

        IO2AB0F =IRHSR2(LABELA,.FALSE.,-FREQEX+FREQB,ISYMA,
     &                  LABELB,.FALSE.,-FREQB,ISYMB)
        IO2ABF0 =IRHSR2(LABELA,.FALSE.,+FREQEX-FREQB,ISYMA,
     &                  LABELB,.FALSE.,+FREQB,ISYMB)
        IO2CD0F =IRHSR2(LABELC,.FALSE.,-FREQEX+FREQB,ISYMC,
     &                  LABELD,.FALSE.,-FREQB,ISYMD)
        IO2CDF0 =IRHSR2(LABELC,.FALSE.,+FREQEX-FREQB,ISYMC,
     &                  LABELD,.FALSE.,+FREQB,ISYMD)

*           request M vectors for different excitation energies

             IOFFST = ISYOFE(ISYMAB) +  ISMSEL(IFREQ,2)
             EIGV   = EIGVAL(IOFFST)
             IM11   = ILRMAMP(IOFFST,EIGV,ISYMAB)

*            calculate left  moment AB-w contribution

             IOPT = 3
             CALL CC_RDRSP('X2',IX2AB0F,ISYMAB,IOPT,MODFIL,
     &                     WRK(K1VEC1),WRK(K1VEC2))
             X1 = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K1VEC1),1)
             IF (.NOT.CCS) THEN
               X2 = DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K1VEC2),1) 
             ELSE
               X2 = ZERO
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) ' NORM^2 of X2 vector:',X1,X2,X1+X2

             IOPT = 3
             CALL CC_RDRSP('RE',IOFFST,ISYMAB,IOPT,MODFIL,
     &                      WRK(K2VEC1),WRK(K2VEC2))
             Y1 = DDOT(NCCVAR1,WRK(K2VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               Y2 = DDOT(NCCVAR2,WRK(K2VEC2),1,WRK(K2VEC2),1) 
             ELSE
               Y2 = ZERO
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) ' NORM^2 of RE vector:',Y1,Y2,Y1+Y2

             SMCLM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMCLM=SMCLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) ' SMCLM:',SMCLM

             IOPT = 3
             CALL CC_RDRSP('M1',IM11,ISYMAB,IOPT,MODFIL,
     &                      WRK(K1VEC1),WRK(K1VEC2))
             X1 = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K1VEC1),1)
             IF (.NOT.CCS) THEN
               X2 = DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K1VEC2),1) 
             ELSE
               X2 = ZERO
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) ' NORM^2 of M1 vector:',X1,X2,X1+X2

             IOPT = 3
             CALL CC_RDRSP('O2',IO2AB0F,ISYMAB,IOPT,MODFIL,
     &                     WRK(K2VEC1),WRK(K2VEC2))
             Y1 = DDOT(NCCVAR1,WRK(K2VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               Y2 = DDOT(NCCVAR2,WRK(K2VEC2),1,WRK(K2VEC2),1) 
             ELSE
               Y2 = ZERO
             END IF
             WRITE (LUPRI,*) ' norm of O2 vector:',y1,y2,y1+y2

             CALL CCLR_DIASCL(WRK(K2VEC2),FACTOR,ISYMAB)

cs
             testp1 = ddot(nccvar1,wrk(k1vec1),1,wrk(k2vec1),1)
             if (.not.ccs) then
               testp2 = ddot(nccvar2,wrk(k1vec2),1,wrk(k2vec2),1)
             else
               testp2 = 0.0d0
             end if
             WRITE (LUPRI,*) ' product Mf*CSIab :', testp1 + testp2
cs

             SMCLM = SMCLM + DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMCLM=SMCLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF
             WRITE (LUPRI,*) ' smclm:',smclm

*            calculate right moment CD+w contribution

             IOPT = 3
             CALL CC_RDRSP('LE',IOFFST,ISYMAB,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             IOPT = 3
             CALL CC_RDRSP('O2',IO2CDF0,ISYMAB,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),FACTOR,ISYMAB)

             SMCRM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMCRM=SMCRM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF

*            calculate left  moment CD-w contribution

             IOPT = 3
             CALL CC_RDRSP('X2',IX2CD0F,ISYMAB,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))

             IOPT = 3
             CALL CC_RDRSP('RE',IOFFST,ISYMAB,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             SMLM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMLM = SMLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF

             IOPT = 3
             CALL CC_RDRSP('M1',IM11,ISYMAB,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))

             IOPT = 3
             CALL CC_RDRSP('O2',IO2CD0F,ISYMAB,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),FACTOR,ISYMAB)

             SMLM = SMLM + DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               SMLM = SMLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF

*            calculate right moment AB+w contribution 

             IOPT = 3
             CALL CC_RDRSP('LE',IOFFST,ISYMAB,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             IOPT = 3
             CALL CC_RDRSP('O2',IO2ABF0,ISYMAB,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),FACTOR,ISYMAB)

             SMRM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1) 
             IF (.NOT.CCS) THEN
               SMRM = SMRM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF


             WRITE(LUPRI,'(/A,/I4,3X,I5,3X,F16.6)')
     *       ' SYM   STATE NUM   EXCITATION ENERGY',
     *        ISYMAB, ISMSEL(IFREQ,2), EIGV
             WRITE(LUPRI,'(/A,/8A)')
     *       'OPERATOR LABELS. '
     *       ,' A: ', LABELA,' B: ',LABELB,' C: ',LABELC,' D: ',LABELD
             WRITE(LUPRI,'(/A,F12.6)')
     *       ' LASER FREQUENCES   W: ', FREQB
             WRITE(LUPRI,'(/A,/A)')
     *       ' TRANSITION STRENGTH FOR SECOND ORDER MOMENT ' ,
     *       ' TSM = 0.5*( ABLMOM * CDRMOM + CDLMOM * ABRMOM ) '

             WRITE (LUPRI,*) ' TSM: ',0.5D0*(smclm*smcrm+smlm*smrm)
             WRITE (LUPRI,*) ' ABLMOM',smclm,' CDRMOM',smcrm
             WRITE (LUPRI,*) ' CDLMOM',smlm,' ABRMOM',smrm
             WRITE (LUPRI,*) ' w moment product  (ablmom*cdlmom): ',
     &            smclm*smcrm
             WRITE (LUPRI,*) ' -w moment product (cdlmom*abrmom): ', 
     &            smlm*smrm

          END DO
        END IF

      END DO


      RETURN
      END
*=====================================================================*
*---------------------------------------------------------------------*
       SUBROUTINE CC_SMSORT
*---------------------------------------------------------------------*
*
*    Purpose: sort the selected states for which second moment 
*             calculation is carried. if no selected states are
*             chosen use all states specified in excitation
*
*    Written by: Poul Joergensen and Christof Haettig  1997
*    Clean up  : Sonia Coriani 2000                            
*=====================================================================*

#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#include <priunit_dec.h>
#include <ccorb_dec.h>
#include <cclr_dec.h>
#else
#include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsm.h>
#include <ccsminf.h>
#include <ccexci.h>
#include <cclr.h>


* local parameters:

      INTEGER ISYM, IST, ISEL, I, ISAVE, JSEL, J, IOFF 
      INTEGER ISYMSV, ISTSV, JSTSV, ISTATE 
      LOGICAL LOCDBG
C

#if defined (SYS_CRAY)
      REAL D2, BSMFRSV
#else
      DOUBLE PRECISION D2, BSMFRSV
#endif
      PARAMETER ( D2 = 2.0d00 )
      PARAMETER ( LOCDBG = .FALSE. )

C
C sort the selected states for which second order transition
C matrix elements are calculated
C
      DO 50 ISYM = 1,NSYM
         NSMSELX(ISYM) = 0
 50   CONTINUE
C
      IF ( SELSMST ) THEN
C
C sort list according to symmetry
C
         ISMSELX(1) = 0
         DO 100 ISYM = 1,NSYM
            IST = ISMSELX(ISYM) + 1
            DO 200 I = IST,NSMSEL
               IF ( ISMSEL(I,1).EQ.ISYM) THEN
                  NSMSELX(ISYM) = NSMSELX(ISYM) + 1
               ELSE
                  DO 300 J = I+1,NSMSEL
                     IF ( ISMSEL(J,1).EQ.ISYM) THEN
                        ISYMSV = ISMSEL(J,1)          
                        ISTSV  = ISMSEL(J,2)          
                        BSMFRSV = BSMFR(J)
                        ISMSEL(J,1) = ISMSEL(I,1)
                        ISMSEL(J,2) = ISMSEL(I,2)
                        BSMFR(J)   = BSMFR(I)
                        ISMSEL(I,1) = ISYMSV
                        ISMSEL(I,2) = ISTSV
                        BSMFR(I)   = BSMFRSV
                        NSMSELX(ISYM) = NSMSELX(ISYM) + 1
                        GO TO 200
                     END IF
 300              CONTINUE
               END IF
 200        CONTINUE
            IF ( ISYM .LT. NSYM ) THEN
               ISMSELX(ISYM+1) = ISMSELX(ISYM) + NSMSELX(ISYM)
            END IF
            IF (LOCDBG)
     &         WRITE (LUPRI,*) 
     &               'CC_SMSORT:',ISMSELX(ISYM),NSMSELX(ISYM),IST
 100     CONTINUE
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) ' after sort of  symmetry '
            WRITE(LUPRI,*) 'nsmsel',nsmsel
            do 210 i = 1,nsmsel
               WRITE(LUPRI,*) ' ismsel(i,1),ismsel(i,2),i',
     *                      ismsel(i,1),ismsel(i,2),i
 210        continue     
            do 211 i = 1,nsym
               WRITE(LUPRI,*) ' ismselx(i),nsmselx(i),i',
     *                      ismselx(i),nsmselx(i),i
 211        continue
         END IF
C
C sort list according to state number
C 
         DO 400 ISYM = 1,NSYM
            IOFF = ISMSELX(ISYM)
            DO 500 ISEL = 1,NSMSELX(ISYM)
               I = IOFF + ISEL
               ISTSV = ISMSEL(I,2)
               ISAVE  = I
               DO 600 JSEL = ISEL+1,NSMSELX(ISYM)
                  J = IOFF + JSEL
                  JSTSV = ISMSEL(J,2)
                  IF ( JSTSV.LT. ISTSV ) THEN 
                     ISTSV  = JSTSV
                     ISAVE  = J
                  END IF
 600           CONTINUE
               IF ( I.NE.ISAVE ) THEN
                  ISYMSV = ISMSEL(ISAVE,1)          
                  ISTSV  = ISMSEL(ISAVE,2)          
                  BSMFRSV = BSMFR(ISAVE)
                  ISMSEL(ISAVE,1) = ISMSEL(I,1)
                  ISMSEL(ISAVE,2) = ISMSEL(I,2)
                  BSMFR(ISAVE)   = BSMFR(I)
                  ISMSEL(I,1) = ISYMSV
                  ISMSEL(I,2) = ISTSV
                  BSMFR(I)   = BSMFRSV
               END IF
 500        CONTINUE
 400     CONTINUE
         IF (LOCDBG) THEN
           WRITE(LUPRI,*) ' after sort of both symmetry and state'
           WRITE(LUPRI,*) 'nsmsel',nsmsel
           do 212 i = 1,nsmsel
              WRITE(LUPRI,*) ' ismsel(i,1),ismsel(i,2),i',
     *                     ismsel(i,1),ismsel(i,2),i
 212       continue     
           do 213 i = 1,nsym
             WRITE(LUPRI,*) ' ismselx(i),nsmselx(i),i',
     *                    ismselx(i),nsmselx(i),i
 213       continue
         END IF
C
C if .HALFFR not specified find frequencies for AOPERATOR
C
         DO 550 ISYM = 1,NSYM
            IOFF = ISMSELX(ISYM)

            IF (LOCDBG) 
     &        WRITE (LUPRI,*) 'CC_SMSORT: ISYM,IOFF: ',ISYM,IOFF

            DO 560 I = 1,NSMSELX(ISYM)
               ISTSV = ISMSEL(IOFF+I,2)
               IF (LOCDBG) THEN
                 WRITE(LUPRI,*) 'istsv,ioff,isym,i', istsv,ioff,isym,i
                 WRITE(LUPRI,*) ' isyofe(isym)', isyofe(isym)
                 WRITE(LUPRI,*) ' eigval(1)', eigval(1)
                 WRITE(LUPRI,*) ' eigval(isyofe(isym)+istsv)', 
     &                        eigval(isyofe(isym)+istsv)
                 CALL FLSHFO(LUPRI)
                END IF

                EXSMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTSV)

                IF (LOCDBG) print*,'EXSMFR(IOFF+I) ', EXSMFR(IOFF+I)

 560        CONTINUE
 550     CONTINUE
      END IF
C
C if selected states not specified for second moment calculations
C then carry out calculations for all specified excited states
C and use frequencies that are half the excitation energy
C
      IF ( .NOT. SELSMST ) THEN
         ISMSELX(1) = 0 
         NSMSEL = 0
         DO 700 ISYM = 1,NSYM
            DO 750 I = 1,NCCEXCI(ISYM,1)
               NSMSEL = NSMSEL + 1
               ISMSEL(NSMSEL,1) = ISYM
               ISMSEL(NSMSEL,2) = I
               NSMSELX(ISYM)    = NSMSELX(ISYM) + 1
 750        CONTINUE
            ISMSELX(ISYM+1) = ISMSELX(ISYM) + NSMSELX(ISYM)
 700     CONTINUE
         HALFFR = .TRUE.
      END IF
C 
C
      IF (HALFFR) THEN
         DO 800  ISYM = 1,NSYM
            IOFF = ISMSELX(ISYM)
            DO 850 I = 1,NSMSELX(ISYM)
               ISTATE         = ISMSEL(IOFF+I,2)
               BSMFR(IOFF+I)  = EIGVAL(ISYOFE(ISYM)+ISTATE)/ D2
               EXSMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTATE)
 850        CONTINUE
 800     CONTINUE
      END IF 
      IF (LOCDBG) THEN
        CALL FLSHFO(LUPRI)
        WRITE(LUPRI,*) ' leaving sort'
        do 215 i = 1,nsmsel
           WRITE(LUPRI,*)' ismsel(i,1),ismsel(i,2),exsmfr(i),i',
     *                 ismsel(i,1),ismsel(i,2),exsmfr(i),i
           CALL FLSHFO(LUPRI)
215     CONTINUE     
        CALL FLSHFO(LUPRI)
      END IF
C
      RETURN
      END 
*=====================================================================*
