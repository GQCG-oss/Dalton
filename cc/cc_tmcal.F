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
c /* deck cc_tmcal */
*=====================================================================*
       SUBROUTINE CC_TMCAL(WRK,LWRK) 
*---------------------------------------------------------------------*
*
*    Purpose: Third moment calculations
*    
*    Written by: P.Joergensen and C.Haettig 1997
*    Clean up and new output style: S. Coriani 2001
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE  
#else
# include <implicit.h>
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
      CHARACTER*8 LABELA, LABELB, LABELC,LABELD, LABELE, LABELF 
      CHARACTER MODFIL*10, MODPRI*5
      INTEGER ISYMB, ISYMC, ISYMA, ISYMD, ISYME, ISYMF, ISYMABC
      INTEGER IFREQ, INUM, IOPER, IDX, IOFFST, LWRK, IOPTRD
      INTEGER K1VEC1, K1VEC2, K2VEC1, K2VEC2, IM11
      INTEGER IX3AC0F, IX3DF0F, IO3AC0F, IO3ACF0, IO3DF0F, IO3DFF0
      INTEGER NCCVAR1, NCCVAR2

#if defined (SYS_CRAY)
      REAL HALF, FREQEX, FREQB, FREQC, EIGV, WRK(LWRK)
      REAL SMLM, SMCLM, SMRM, SMCRM
      REAL ABCLM,DEFLM,ABCRM,DEFRM
      REAL X1, X2, Y1, Y2
      REAL DDOT, ZERO
#else
      DOUBLE PRECISION HALF, FREQEX, FREQB, FREQC, EIGV, WRK(LWRK) 
      DOUBLE PRECISION SMLM, SMCLM, SMRM, SMCRM
      DOUBLE PRECISION ABCLM,DEFLM,ABCRM,DEFRM,THREEPH
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
*---------------------------------------------------------------------*
* print header for third order moments section
*---------------------------------------------------------------------*
      WRITE (LUPRI,'(7(/1X,2A),/)')
     & '************************************',
     &                               '*******************************',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<< OUTPUT FROM COUPLED CLUSTER C',
     &                                  'UBIC RESPONSE >>>>>>>>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '*<<<<<<<< CALCULATION OF THREE PHOTON TRANS',
     &                                      'ITION STRENGTHS >>>>>>>*',
     & '*                                   ',
     &                               '                              *',
     & '************************************',
     &                               '*******************************'

*---------------------------------------------------------------------*
* print debug info
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'DEBUG_CC_TMIND> NTMOP = ',NTMOPER
      END IF
*---------------------------------------------------------------------*
* set MODFIL, MODPRI, IOPTRD for calls to CC_RDRSP and print out
*---------------------------------------------------------------------*
      IF (CCS) THEN
         MODFIL = 'CCS       '
         MODPRI  = 'CCS  '
         IOPTRD = 1
      ELSE IF (CC2) THEN
         MODFIL = 'CC2       '
         MODPRI  = 'CC2  '
         IOPTRD = 3
      ELSE IF (CCSD) THEN
         MODFIL = 'CCSD      '
         MODPRI  = 'CCSD '
         IOPTRD = 3
      ELSE
         CALL QUIT('Unknown coupled cluster model in CC_TMCAL')
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

        IF (LOCDBG) THEN
           WRITE (LUPRI,*) 'LABELA:',LABELA, ' ISYMA:', ISYMA
           WRITE (LUPRI,*) 'LABELB:',LABELB, ' ISYMB:', ISYMB
           WRITE (LUPRI,*) 'LABELC:',LABELC, ' ISYMC:', ISYMC
           WRITE (LUPRI,*) 'LABELD:',LABELD, ' ISYMD:', ISYMD
           WRITE (LUPRI,*) 'LABELE:',LABELE, ' ISYME:', ISYME
           WRITE (LUPRI,*) 'LABELF:',LABELF, ' ISYMF:', ISYMF
           CALL FLSHFO(LUPRI)
        END IF
        
        ISYMABC = MULD2H(MULD2H(ISYMA,ISYMB),ISYMC)
        IF (ISYMABC .EQ. MULD2H(ISYMD,MULD2H(ISYME,ISYMF))) THEN
      
          NCCVAR1 = NT1AM(ISYMABC)
          NCCVAR2 = NT2AM(ISYMABC)
          K1VEC1   = 1
          K1VEC2   = K1VEC1 + NCCVAR1
          K2VEC1   = K1VEC2 + NCCVAR2
          K2VEC2   = K2VEC1 + NCCVAR1
     
          DO I = 1, NTMSELX(ISYMABC)  
C bug fix
c             IFREQ  = ITMSELX( MULD2H(ISYMA,ISYMB) ) + I
              IFREQ  = ITMSELX(ISYMABC) + I
              FREQEX  = EXTMFR(IFREQ)
              FREQB  = BTMFR(IFREQ)
              FREQC  = CTMFR(IFREQ)
              IF (LOCDBG) THEN
                WRITE (LUPRI,*) 'CC_TMCAL> put on the list:',
     &            LABELA,'(',FREQEX,'),  ', LABELB,'(',FREQB ,'),  ',
     &            LABELC,'(',FREQC, '),  ', LABELD,'(',FREQEX,'),  ',
     &            LABELE,'(',FREQB, '),  ', LABELF,'(',FREQC ,'),  '
              END IF

*    request third order chi vectors:

           IX3AC0F = ICHI3(LABELA,-FREQEX+FREQB+FREQC,ISYMA,
     &                  LABELB,-FREQB,ISYMB,LABELC,-FREQC,ISYMC)
           IX3DF0F = ICHI3(LABELD,-FREQEX+FREQB+FREQC,ISYMD,
     &                  LABELE,-FREQB,ISYME,LABELF,-FREQC,ISYMF)

*    request third order rhs vectors

           IO3AC0F = IRHSR3(LABELA,-FREQEX+FREQB+FREQC,ISYMA,
     &                  LABELB,-FREQB,ISYMB,LABELC,-FREQC,ISYMC)
           IO3ACF0 = IRHSR3(LABELA,+FREQEX-FREQB-FREQC,ISYMA,
     &                  LABELB,+FREQB,ISYMB,LABELC,+FREQC,ISYMC)
           IO3DF0F = IRHSR3(LABELD,-FREQEX+FREQB+FREQC,ISYMD,
     &                  LABELE,-FREQB,ISYME,LABELF,-FREQC,ISYMF)
           IO3DFF0 = IRHSR3(LABELD,+FREQEX-FREQB-FREQC,ISYMD,
     &                  LABELE,+FREQB,ISYME,LABELF,+FREQC,ISYMF)

*    request M vectors for different excitation energies

             IOFFST = ISYOFE(ISYMABC) +  ITMSEL(IFREQ,2)
             EIGV   = EIGVAL(IOFFST)
             IM11   = ILRMAMP(IOFFST,EIGV,ISYMABC)
*--------------------------------------------------------------*
*            calculate left  moment M_of^ABC(-w1,-w2) contrib.
*            previously called SMCLM
*--------------------------------------------------------------*

             CALL CC_RDRSP('X3',IX3AC0F,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             X1 = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K1VEC1),1)
             IF (.NOT.CCS) THEN
               X2 = DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K1VEC2),1) 
             ELSE
               X2 = ZERO
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) ' norm^2 of X3 vector:',X1,X2

             CALL CC_RDRSP('RE',IOFFST,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             Y1 = DDOT(NCCVAR1,WRK(K2VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               Y2 = DDOT(NCCVAR2,WRK(K2VEC2),1,WRK(K2VEC2),1) 
             ELSE
               Y2 = ZERO
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) ' norm^2 of RE vector:',Y1,Y2

             ABCLM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               ABCLM=ABCLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF

             CALL CC_RDRSP('M1',IM11,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             X1 = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K1VEC1),1)
             IF (.NOT.CCS) THEN
               X2 = DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K1VEC2),1) 
             ELSE
               X2 = ZERO
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) 'Norm^2 of M1:',X1,X2,X1+X2

             CALL CC_RDRSP('O3',IO3AC0F,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             Y1 = DDOT(NCCVAR1,WRK(K2VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               Y2 = DDOT(NCCVAR2,WRK(K2VEC2),1,WRK(K2VEC2),1) 
             ELSE
               Y2 = ZERO
             END IF
             IF (LOCDBG)
     &          WRITE (LUPRI,*) ' Norm^2 of O3 vector:',y1,y2,y1+y2

             CALL CCLR_DIASCL(WRK(K2VEC2),HALF,ISYMABC)

             ABCLM = ABCLM + DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               ABCLM=ABCLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF

*--------------------------------------------------------------*
*            calculate right moment M_fo^DEF(w1,w2) contribution
*            previously called SMCRM
*--------------------------------------------------------------*

             CALL CC_RDRSP('LE',IOFFST,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             CALL CC_RDRSP('O3',IO3DFF0,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),HALF,ISYMABC)

             DEFRM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               DEFRM=DEFRM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1) 
             END IF

*--------------------------------------------------------------*
*            calculate left moment M_of^DEF(-w1,-w2) contrib.
*            (previously SMLM)
*--------------------------------------------------------------*

             CALL CC_RDRSP('X3',IX3DF0F,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))

             CALL CC_RDRSP('RE',IOFFST,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             DEFLM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               DEFLM = DEFLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF

             CALL CC_RDRSP('M1',IM11,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))

             CALL CC_RDRSP('O3',IO3DF0F,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),HALF,ISYMABC)

             DEFLM = DEFLM + DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1)
             IF (.NOT.CCS) THEN
               DEFLM = DEFLM + DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF

*--------------------------------------------------------------*
*            calculate right moment M_fo^ABC(w1,w2) contribution 
*            (previously SMRM)
*--------------------------------------------------------------*

             CALL CC_RDRSP('LE',IOFFST,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K1VEC1),WRK(K1VEC2))
             CALL CC_RDRSP('O3',IO3ACF0,ISYMABC,IOPTRD,MODFIL,
     *                     WRK(K2VEC1),WRK(K2VEC2))
             CALL CCLR_DIASCL(WRK(K2VEC2),HALF,ISYMABC)

             ABCRM = DDOT(NCCVAR1,WRK(K1VEC1),1,WRK(K2VEC1),1) 
             IF (.NOT.CCS) THEN
               ABCRM = ABCRM+DDOT(NCCVAR2,WRK(K1VEC2),1,WRK(K2VEC2),1)
             END IF
*--------------------------------------------------------------*
*            Final three-photon transition strength:
*--------------------------------------------------------------*
              
             THREEPH = HALF*(ABCLM*DEFRM+DEFLM*ABCRM)

*--------------------------------------------------------------*
* Write results on output
*--------------------------------------------------------------*
         WRITE(LUPRI,'(65("-"),/1x,a,f10.6,a,i1,a,i1)')
     &  'For trans. to |f(',EIGV,')>, state nr. ',ITMSEL(IFREQ,2),
     &                              ' of symm. ',ISYMABC
         WRITE(LUPRI,'(/,3(1x,a5,a,a1,i1,a1))')
     &     ' A:  ',LABELA,'(',ISYMA,')', '; B: ',LABELB,'(',ISYMB,')',
     &     '; C: ',LABELC,'(',ISYMC,')' 
         WRITE(LUPRI,'(3(1x,a5,a,a1,i1,a1))')
     &     ' D:  ',LABELD,'(',ISYMD,')', '; E: ',LABELE,'(',ISYME,')',
     &     '; F: ',LABELF,'(',ISYMF,')' 
         WRITE(LUPRI,'(1x,a,f10.6,a,f10.6)')
     &     ' Laser frequencies (au): w1 = ', FREQB, '; w2 = ', FREQC
C         IF (LOCDBG) THEN
         WRITE(LUPRI,'(2(/1x,a,f15.9,1x,a,f15.9))')
     & ' M^ABC_of(-w1,-w2): ',ABCLM,' M^DEF_fo(w1,w2): ',DEFRM,
     & ' M^DEF_of(-w1,-w2): ',DEFLM,' M^ABC_fo(w1,w2): ',ABCRM
         WRITE(LUPRI,'(2(1x,a,f15.9,/))')
     & ' M^ABC_of(-w1,-w2) x M^DEF_fo(w1,w2)   = ', abclm*defrm,
     & '[M^DEF_of(-w1,-w2) x M^ABC_fo(w1,w2)]* = ', deflm*abcrm
C         END IF
        WRITE(LUPRI,'(1x,a5,a,/,1x,a5,a,f10.6,a1,f10.6,a,f15.9)')
     &   MODPRI,'Transition strength for Third Order Moment: ',
     &   MODPRI,'S^of_ABC,DEF(',FREQB,',',FREQC,') = ', THREEPH
        WRITE(LUPRI,'(65("-"))')

*--------------------------------------------------------------*

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
*    P. Joergensen, C. Haettig 1997
*    Clean up, new output. Sonia 2001
*=====================================================================*

#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
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
      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

#if defined (SYS_CRAY)
      REAL D3, BTMFRSV, CTMFRSV
#else
      DOUBLE PRECISION D3, BTMFRSV, CTMFRSV
#endif
      PARAMETER ( D3 = 3.0D00 )

C
C sort the selected states for which third order transition
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
            IF (LOCDBG)
     &      WRITE (LUPRI,*) 'SORT:',ITMSELX(ISYM),NTMSELX(ISYM),IST
 100     CONTINUE
         IF (LOCDBG) THEN
           WRITE (LUPRI,*) ' after sort of  symmetry '
           WRITE (LUPRI,*) 'ntmsel',ntmsel
           do 210 i = 1,ntmsel
             WRITE (LUPRI,*) ' itmsel(i,1),itmsel(i,2),i'
             WRITE (LUPRI,*)  itmsel(i,1),itmsel(i,2),i
 210       continue     
           do 211 i = 1,nsym
             WRITE (LUPRI,*) ' itmselx(i),ntmselx(i),i'
             WRITE (LUPRI,*) itmselx(i),ntmselx(i),i
 211       continue
         END IF
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
         IF (LOCDBG) THEN
           WRITE (LUPRI,*) ' after sort of both symmetry and state'
           WRITE (LUPRI,*) 'ntmsel',ntmsel
           do 212 i = 1,ntmsel
             WRITE (LUPRI,*) ' itmsel(i,1),itmsel(i,2),i'
             WRITE (LUPRI,*) itmsel(i,1),itmsel(i,2),i
 212       continue     
           do 213 i = 1,nsym
             WRITE (LUPRI,*) ' itmselx(i),ntmselx(i),i'
             WRITE (LUPRI,*) itmselx(i),ntmselx(i),i
 213       continue
         END IF
C
C if .HALFFR not specified find frequencies for AOPERATOR
C
         DO 550 ISYM = 1,NSYM
            IOFF = ITMSELX(ISYM)
            WRITE (LUPRI,*) 'isym, ioff', isym, ioff
            DO 560 I = 1,NTMSELX(ISYM)
               ISTSV = ITMSEL(IOFF+I,2)
               EXTMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTSV)

               IF (LOCDBG) THEN
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
                  WRITE (LUPRI,*) ' EXTMFR(IOFF+I) '
                  call flshfo(LUPRI)
                  WRITE (LUPRI,*) EXTMFR(IOFF+I) 
                  call flshfo(LUPRI)
               END IF
 560        CONTINUE
         IF (LOCDBG) THEN
            WRITE (LUPRI,*) ' isym loop slut',isym
            call flshfo(LUPRI)
         END IF
 550     CONTINUE
      END IF
C
C if selected states not specified for second moment calculations
C then carry out calculations for all specified excited states
C and use frequencies that are half the excitation energy
C
      IF ( .NOT. SELTMST ) THEN
         ITMSELX(1) = 0 
         NTMSEL = 0
         DO 700 ISYM = 1,NSYM
            DO 750 I = 1,NCCEXCI(ISYM,1)
               NTMSEL = NTMSEL + 1
               ITMSEL(NTMSEL,1) = ISYM
               ITMSEL(NTMSEL,2) = I
               NTMSELX(ISYM)    = NTMSELX(ISYM) + 1
 750        CONTINUE
            ITMSELX(ISYM+1) = ITMSELX(ISYM) + NTMSELX(ISYM)
 700     CONTINUE
         THIRDFR = .TRUE.
      END IF
C 
C
      IF (THIRDFR) THEN
         DO 800  ISYM = 1,NSYM
            IOFF = ITMSELX(ISYM)
            DO 850 I = 1,NTMSELX(ISYM)
               ISTATE = ITMSEL(IOFF+I,2)
               BTMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTATE)/ D3
               CTMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTATE)/ D3
               EXTMFR(IOFF+I) = EIGVAL(ISYOFE(ISYM)+ISTATE)
 850        CONTINUE
 800     CONTINUE
      END IF 
         IF (LOCDBG) THEN
           WRITE(LUPRI,*) ' leaving sort'
           do i = 1,ntmsel
           WRITE(LUPRI,*) ' itmsel(i,1),itmsel(i,2),extmfr(i),i'
           call flshfo(LUPRI)
           WRITE(LUPRI,*) itmsel(i,1),itmsel(i,2),extmfr(i),i
           end do     
           call flshfo(LUPRI)
         END IF
                
      RETURN
      END 
*=====================================================================*
