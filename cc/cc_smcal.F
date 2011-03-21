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
#include <single.h>
C
*---------------------------------------------------------------------*
c /* deck cc_smcal */
*=====================================================================*
       SUBROUTINE CC_SMCAL(WRK,LWRK) 
*---------------------------------------------------------------------*
*
*    Purpose: Second moment calculations
*
*    Written by: Poul Joergensen and Christof Haettig  1997
*    Clean up and new output: Sonia Coriani 2000,2001
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
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
      CHARACTER MODFIL*10, MODPRI*5
      INTEGER ISYMB, ISYMC, ISYMA, ISYMD, ISYMAB
      INTEGER IFREQ, INUM, IOPER, IDX, IOFFST, LWRK, IM11, IOPT
      INTEGER IX2AB0F, IX2CD0F, IO2AB0F, IO2ABF0, IO2CD0F, IO2CDF0
      INTEGER K1VEC1, K1VEC2, K2VEC1, K2VEC2, NCCVAR1, NCCVAR2 
      INTEGER KEND0,KTWOPH

#if defined (SYS_CRAY)
      REAL WRK(LWRK)
      REAL ZERO, FACTOR, FREQEX, FREQB, EIGV, DP5
      REAL ABLM, ABRM, CDLM, CDRM, TSM
      REAL X1, X2, Y1, Y2, TESTP1, TESTP2
      REAL DDOT
#else
      DOUBLE PRECISION WRK(LWRK)
      DOUBLE PRECISION ZERO, FACTOR, FREQEX, FREQB, EIGV, DP5
      DOUBLE PRECISION ABLM, ABRM, CDLM, CDRM, TSM
      DOUBLE PRECISION X1, X2, Y1, Y2, TESTP1, TESTP2
      DOUBLE PRECISION DDOT
#endif

      PARAMETER ( ZERO = 0.0D00, FACTOR = 0.5D00, DP5 = 0.5D00)

* external functions:
      INTEGER IRHSR2
      INTEGER ILRMAMP
      INTEGER ICHI2
* data:
      LOGICAL FIRSTCALL
      SAVE    FIRSTCALL
      DATA    FIRSTCALL /.TRUE./
*---------------------------------------------------------------------*
* print header for second order moments section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<< OUTPUT FROM COUPLED CLUSTER Q',
     &                                  'UADRATIC RESPONSE >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<< CALCULATION OF TWO PHOTON TRANS',
     &                                    'ITION STRENGTHS >>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************'

*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_SMCAL> NSMOP = ',NSMOPER
      END IF
*---------------------------------------------------------------------*
* set MODFIL, MODPRI, IOPT for calls to CC_RDRSP and print out
*---------------------------------------------------------------------*
      IF (CCS) THEN
         MODFIL = 'CCS       '
         MODPRI  = 'CCS  '
         IOPT = 1
      ELSE IF (CC2) THEN
         MODFIL = 'CC2       '
         MODPRI  = 'CC2  '
         IOPT = 3
      ELSE IF (CCSD) THEN
         MODFIL = 'CCSD      '
         MODPRI  = 'CCSD '
         IOPT = 3
      ELSE
         CALL QUIT('Unknown coupled cluster model in CC_TMCAL')
      END IF
*---------------------------------------------------------------------*
*     set length of result matrix (in fieri, Sonia)
*---------------------------------------------------------------------*
*      NTWOPH = NSMOPER * 
*      KTWOPH = 1
*      KEND0  = KTWOPH + 
*      LWRK0  = LWRK - KEND0
*---------------------------------------------------------------------*
* find list entries for the required response vectors
* and excitation vectors:
*---------------------------------------------------------------------*

      KEND0 = 1

* +-------------------------------------------------------------------------------------------+
* | SYM | STATE#| EXCIT.EN. |    A    |    B    |    C    |    D    | w (au) | S^of_AB,CD (w) |
* +-------------------------------------------------------------------------------------------+
*      WRITE(LUPRI,'(/1x,(70("-")))')
*      WRITE(LUPRI,'(3(a))')
*     & '| SYM | STATE#| EXCIT.EN. |    A    |    B    |',
*     &                            '    C    |    D    |',
*     &                            ' w (au) | S^of_AB,CD (w) |'

      DO IOPER = 1, NSMOPER
        LABELA = LBLOPR(IASMOP(IOPER))
        LABELB = LBLOPR(IBSMOP(IOPER))
        LABELC = LBLOPR(ICSMOP(IOPER))
        LABELD = LBLOPR(IDSMOP(IOPER))

        ISYMA  = ISYOPR(IASMOP(IOPER))
        ISYMB  = ISYOPR(IBSMOP(IOPER))
        ISYMC  = ISYOPR(ICSMOP(IOPER))
        ISYMD  = ISYOPR(IDSMOP(IOPER))

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'LABELA:',LABELA
           WRITE (LUPRI,*) 'LABELB:',LABELB
           WRITE (LUPRI,*) 'LABELC:',LABELC
           WRITE (LUPRI,*) 'LABELD:',LABELD

           WRITE (LUPRI,*) 'ISYMA:', ISYMA
           WRITE (LUPRI,*) 'ISYMB:', ISYMB
           WRITE (LUPRI,*) 'ISYMC:', ISYMC
           WRITE (LUPRI,*) 'ISYMD:', ISYMD

           CALL FLSHFO(LUPRI)
        END IF
        
        IF (MULD2H(ISYMA,ISYMB).EQ.MULD2H(ISYMC,ISYMD)) THEN
      
          ISYMAB = MULD2H(ISYMA,ISYMB)
          NCCVAR1 = NT1AM(ISYMAB)
          NCCVAR2 = NT2AM(ISYMAB)
          K1VEC1  = KEND0
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

*       request second order chi vectors (for left moments):
 
        IX2AB0F = ICHI2(LABELA,.FALSE.,-FREQEX+FREQB,ISYMA,
     &                  LABELB,.FALSE.,-FREQB,ISYMB)
        IX2CD0F = ICHI2(LABELC,.FALSE.,-FREQEX+FREQB,ISYMC,
     &                  LABELD,.FALSE.,-FREQB,ISYMD)
C
C       IX2CD0F = ICHI2(LABELC,.FALSE.,-FREQEX-FREQB,ISYMC,
C    &                  LABELD,.FALSE.,+FREQB,ISYMD)
C

*       request second-order rhs vectors (for both left and right mom):

        IO2AB0F =IRHSR2(LABELA,.FALSE.,-FREQEX+FREQB,ISYMA,
     &                  LABELB,.FALSE.,-FREQB,ISYMB)
        IO2ABF0 =IRHSR2(LABELA,.FALSE.,+FREQEX-FREQB,ISYMA,
     &                  LABELB,.FALSE.,+FREQB,ISYMB)
        IO2CD0F =IRHSR2(LABELC,.FALSE.,-FREQEX+FREQB,ISYMC,
     &                  LABELD,.FALSE.,-FREQB,ISYMD)
        IO2CDF0 =IRHSR2(LABELC,.FALSE.,+FREQEX-FREQB,ISYMC,
     &                  LABELD,.FALSE.,+FREQB,ISYMD)

*       request M vectors for different excitation energies (for left mom)

             IOFFST = ISYOFE(ISYMAB) +  ISMSEL(IFREQ,2)
             EIGV   = EIGVAL(IOFFST)
             IM11   = ILRMAMP(IOFFST,EIGV,ISYMAB)

*--------------------------------------------------------------------*
*            calculate left  moment M^AB_of(-w) contribution
*--------------------------------------------------------------------*
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

             ABLM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               ABLM=ABLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) ' ABLM',ABLM

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

             CALL CC_RDRSP('O2',IO2AB0F,ISYMAB,IOPT,MODFIL,
     &                     WRK(K2VEC1),WRK(K2VEC2))
             Y1 = DDOT(NCCVAR1,WRK(K2VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               Y2 = DDOT(NCCVAR2,WRK(K2VEC2),1,WRK(K2VEC2),1) 
             ELSE
               Y2 = ZERO
             END IF
             IF (LOCDBG) 
     &          WRITE (LUPRI,*) ' norm of O2 vector:',y1,y2,y1+y2

             CALL CCLR_DIASCL(WRK(K2VEC2),FACTOR,ISYMAB)

             ABLM = ABLM + DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               ABLM=ABLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) 'Left moment M^AB_of (-w) :', ABLM
*--------------------------------------------------------------------*
*            calculate right moment M^CD_fo(+w) contribution
*--------------------------------------------------------------------*
             CALL CC_RDRSP('LE',IOFFST,ISYMAB,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             CALL CC_RDRSP('O2',IO2CDF0,ISYMAB,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),FACTOR,ISYMAB)

             CDRM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               CDRM=CDRM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF
             IF (LOCDBG) 
     &           WRITE (LUPRI,*) 'Right moment M^CD_fo(w) :', CDRM
*--------------------------------------------------------------------*
*            calculate left  moment M^CD_of(-w) contribution
*--------------------------------------------------------------------*
             CALL CC_RDRSP('X2',IX2CD0F,ISYMAB,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))

             CALL CC_RDRSP('RE',IOFFST,ISYMAB,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CDLM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               CDLM = CDLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF

             CALL CC_RDRSP('M1',IM11,ISYMAB,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))

             CALL CC_RDRSP('O2',IO2CD0F,ISYMAB,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),FACTOR,ISYMAB)

             CDLM = CDLM + DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               CDLM = CDLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF
             IF (LOCDBG) 
     &          WRITE (LUPRI,*) 'Left moment M^CD_of (-w) :', CDLM
*--------------------------------------------------------------------*
*            calculate right moment M^AB_fo(+w) contribution 
*--------------------------------------------------------------------*
             CALL CC_RDRSP('LE',IOFFST,ISYMAB,IOPT,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             CALL CC_RDRSP('O2',IO2ABF0,ISYMAB,IOPT,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),FACTOR,ISYMAB)

             ABRM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1) 
             IF (.NOT.CCS) THEN
               ABRM = ABRM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF
             IF (LOCDBG) 
     &          WRITE (LUPRI,*) 'Right moment M^AB_fo (w) :', ABRM
*--------------------------------------------------------------------*
*             Total strength
*--------------------------------------------------------------------*

             TSM = DP5*(ABLM*CDRM+CDLM*ABRM)

*--------------------------------------------------------------------*
*      Write output
*---------------------------------------------------------------------*
         WRITE(LUPRI,'(1x,65("-"),/1x,a,f10.6,a,i1,a,i1)')
     &  'For trans. to |f(',EIGV,')>, state nr.',ISMSEL(IFREQ,2),
     &                              ' of symm. ',ISYMAB
         WRITE(LUPRI,'(2(/1x,a5,a,a1,i1,a1,a6,a,a1,i1,a1))')
     &     ' A : ',LABELA,'(',ISYMA,')', '; B : ',LABELB,'(',ISYMB,')',
     &     ' C : ',LABELC,'(',ISYMC,')', '; D : ',LABELD,'(',ISYMD,')'
         WRITE(LUPRI,'(/1x,a,f16.6)')
     &     ' Laser frequency w (au) = ', FREQB
         WRITE(LUPRI,'(2(/1x,a,f15.9,1x,a,f15.9))')
     &       ' M^AB_of(-w): ', ablm,' M^CD_fo(w): ', cdrm,
     &       ' M^CD_of(-w): ', cdlm,' M^AB_fo(w): ', abrm
         WRITE(LUPRI,'(2(/1x,a,f15.9))')
     &       ' M^AB_of(-w) x M^CD_fo(w)   = ', ablm*cdrm,
     &       '[M^CD_of(-w) x M^AB_fo(w)]* = ', cdlm*abrm
         WRITE(LUPRI,'(/1x,a5,a,/,1x,a5,a,f10.6,a,f15.9,/1x,65("-"))')
     &     MODPRI,'Transition strength for Second Order Moment: ',
     &     MODPRI,'S^of_AB,CD(',FREQB,') = ',TSM  
*---------------------------------------------------------------------*
*       WRITE(LUPRI,
*     &'(/1x,a2,i1,a3,i5,a3,f8.4,a3,a,a3,a,a3,a,a3,a,f8.7,a3,f16.8,a2)')
*     &'| ',ISYMAB,' | ',ISMSEL(IFREQ,2),' | ', EIGV, ' | ', 
*     &     LABELA,' | ',LABELB,' | ',LABELC,' | ',LABELD,' | ',
*     &      FREQB,' | ',TSM,' |'

          END DO
        END IF
      END DO
*
* It should look like this ......... Sonia
* +-------------------------------------------------------------------------------------------+
* | SYM | STATE#| EXCIT.EN. |    A    |    B    |    C    |    D    | w (au) | S^of_AB,CD (w) |
* +-------------------------------------------------------------------------------------------+
* |  1  |  345  | 0000.0000 | XDIPLEN | XDIPLEN | XDIPLEN | XDIPLEN | 0.0000 | 0000.000000000 |
*
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
      REAL D2, DP5, BSMFRSV
#else
      DOUBLE PRECISION D2, DP5, BSMFRSV
#endif
      PARAMETER ( D2  = 2.0D00, DP5 = 0.5D00 )
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

                IF (LOCDBG) 
     &            WRITE(LUPRI,*) 'EXSMFR(IOFF+I) ', EXSMFR(IOFF+I)

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
               BSMFR(IOFF+I)  = EIGVAL(ISYOFE(ISYM)+ISTATE)* DP5
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
