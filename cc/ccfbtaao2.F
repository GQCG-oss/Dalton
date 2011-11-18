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
*=====================================================================*
      SUBROUTINE CCFBTAAO2(X0INT, ISY0DIS, X1INT, ISY1DIS,  
     &                     XMGD1, XMGD2, 
     &                     BF0RHF, BF1RHF,
     &                     XLAMDP0, XLAMDH0, 
     &                     XLAMDPQ, XLAMDHQ, 
     &                     XLAMDPA, XLAMDHA, 
     &                     XLAMPQA, XLAMHQA, 
     &                     FAIM, FQAIM,
     &                     FILBFZA,  LUBFZA,  IADRBZA, IADRBZ1, 
     &                     FILBFZQA, LUBFZQA, IADRBZQA,IADRBZ2, 
     &                     FNBFDE0, LUBFDE0, IADRE0,
     &                     FNBFDE1, LUBFDE1, IADRE1,
     &                     FILPQ0,  LUPQ0,   IADRPQ0,  
     &                     FILPQ1,  LUPQ1,   IADRPQ1,  
     &                     FILPQA0,  LUPQA0,   IADRPQA0,  
     &                     FILPQA1,  LUPQA1,   IADRPQA1,  
     &                     LRELAX, LTWOEL,  LZERO, LNEWTA, LNEWZ,
     &                     LX1ISQ,
     &                     ISYCTR, ISYMTA,ISYHOP,IREAL,
     &                     IDEL,   WORK,    LWORK)
*---------------------------------------------------------------------*
*
* Purpose:
*
*     Precalculate some intermediates for F^bT^a vector  which depend 
*     on both T^A and Zeta vectors (and IOPER) and require the 2-electrons
*     AO integrals
*
*     LRELAX, LTWOEL not carried thru yet
*     zeroth-order intermediates are only computed for LZERO ?
*     Ta-Zeta dependent only calculated for LNEWTA or LNEWZ
*
*     Input: 
*     X0INT, X1INT   = (al bet ga; del) integrals
*     BF0RHF, BF1RHF = presorted (al be, k;del) integrals (as for BF)
*     XMGD1 = effective density for BZ(A) (symmetry specified inside)
*     XMGD2 = effective density for BZ(QA)(symmetry specified inside)
*
*     Sonia Coriani February 1999
*     Version: 17/11-1999
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE           
#else
#  include <implicit.h>
#endif     
#include <priunit.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
#include <maxorb.h>
#include <ccorb.h>
#include <ccisao.h>

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      LOGICAL LRELAX, LTWOEL, LZERO, LNEWTA, LNEWZ, LX1ISQ
      CHARACTER*(*) FILBFZA, FILBFZQA
      CHARACTER*(*) FILPQ0, FILPQ1, FILPQA0,FILPQA1
      CHARACTER*(*) FNBFDE0, FNBFDE1
      INTEGER ISYCTR, ISYMTA, ISY0DIS, ISY1DIS, IDEL, LWORK
      INTEGER ISYHOP, IREAL
      INTEGER LUBFZA,  IADRBZA(MXCORB)
      INTEGER LUBFZQA, IADRBZQA(MXCORB)
      INTEGER LUPQ0, IADRPQ0(MXCORB)
      INTEGER LUPQ1, IADRPQ1(MXCORB)
      INTEGER LUPQA0, IADRPQA0(MXCORB)
      INTEGER LUPQA1, IADRPQA1(MXCORB)
      INTEGER LUBFDE0,IADRE0(MXCORB)
      INTEGER LUBFDE1,IADRE1(MXCORB)
      
#if defined (SYS_CRAY)
      REAL XLAMDP0(*), XLAMDH0(*), XLAMDPQ(*), XLAMDHQ(*)
      REAL XLAMDPA(*), XLAMDHA(*), XLAMPQA(*), XLAMHQA(*)
      REAL FAIM(*), FQAIM(*)
      REAL XMGD1(*), XMGD2(*), BF0RHF(*), BF1RHF(*)
      REAL X0INT(*), X1INT(*), WORK(LWORK), DUMMY
#else
      DOUBLE PRECISION XLAMDP0(*), XLAMDH0(*), XLAMDPQ(*), XLAMDHQ(*)
      DOUBLE PRECISION XLAMDPA(*), XLAMDHA(*), XLAMPQA(*), XLAMHQA(*)
      DOUBLE PRECISION FAIM(*), FQAIM(*)
      DOUBLE PRECISION XMGD1(*), XMGD2(*), BF0RHF(*), BF1RHF(*)
      DOUBLE PRECISION X0INT(*), X1INT(*), WORK(LWORK), DUMMY
#endif
      LOGICAL LSQRUP
      INTEGER ISYDEL, ISYZTA, NMGD, IADRSTR, KEND1, LWRK1
      INTEGER KPINT0, KQINT0, KPINT1, KQINT1, IOPT
      INTEGER KPAIN0, KQAIN0, KPAIN1, KQAIN1, IADR
      INTEGER ISYRES, ISYBF0, ISYBF1, IDUMMY
      INTEGER ISYZT0, ISYZT1, ISZTA0, ISZTA1, ISYHTA
      INTEGER LEN00, LEN01, LENA0, LENA1, LWRK2, LDSRHF, LENMGD
      INTEGER ISYHZE, ISYABK, ISQABK, ISMGD0, ISMGD1
      INTEGER KDSRHF, KMGD, KEND2, IGZ
      INTEGER IADRBZ1, IADRBZ2
#if defined (SYS_CRAY)
      REAL SIGN, ONE
#else
      DOUBLE PRECISION SIGN, ONE
#endif
      PARAMETER (ONE = 1.0D0)

*---------------------------------------------------------------------*
*     nothing to do for CCS
*---------------------------------------------------------------------*

      ISYDEL = ISAO(IDEL)
      D      = IDEL - IBAS(ISYDEL)

      ISYHTA = MULD2H(ISYMTA,ISYHOP)
      ISYZTA = MULD2H(ISYCTR,ISYMTA)    !sym of Dens(A)
      ISYRES = MULD2H(ISYZTA,ISYHOP)    !sym of Dens(QA) (and result)
      !sym presorted g(0)_(al be,k) for given delta
      ISYBF0 = MULD2H(ISYDEL,ISYM0)
      !sym presorted g(1)_(al be,k) for given delta
      ISYBF1 = MULD2H(ISYDEL,ISYHOP)
* 
*  Note that  ISYBF0 = ISY0DIS*ISYM0 = (ISYDEL*ISYM0)*ISYM0
*             ISYBF1 = ISY1DIS*ISYM0 = (ISYDEL*ISYHOP)*ISYM0 
*                                    = ISY0DIS*ISYHOP
*
*---------------------------------------------------------------------*
*     for CCSD calculate the BZeta intermediates:
*---------------------------------------------------------------------*
      IF (.NOT.CC2) THEN
C
cs as far as I can see LRELAX is the only condition
C
         IF (LRELAX) THEN
            
*        -------------------------------------------------------
*        calculate the BZ(QA) intermediate
*        from the XMGD2 density (QA)  and presorted I(0)+I(1)
*        Result is written on file directly (for each delta)
*        CHECK what IADR is (start address?)
*        multiply by -1 for lao? (pass sign common block)
*        -------------------------------------------------------

            CALL CC_BFIF1(BF0RHF,ISYBF0,XMGD1,ISYZTA,
     &                    BF1RHF,ISYBF1,XMGD2,ISYRES,
     &                    LUBFZQA,FILBFZQA,IADRBZQA,IADRBZ2,IDEL,
     &                                    WORK,LWORK)

            IF (LNEWTA.OR.LNEWZ) THEN     
*           ------------------------------------------------------
*           calculate also BZ(A) intermediate:
*           ------------------------------------------------------

               CALL CC_BFIF(BF0RHF,ISYBF0,XMGD1,ISYZTA,
     &                      LUBFZA,FILBFZA,IADRBZA,IADRBZ1,IDEL,
     &                                  WORK,LWORK)
           END IF
         END IF

      END IF

*---------------------------------------------------------------------*
*     for CC2 calculate the G intermediates:
*---------------------------------------------------------------------*
*      IF (CC2) THEN
*
*        ---------------------------------------------------------
*        ...to be implemented...
*        ---------------------------------------------------------
*
*      END IF
*---------------------------------------------------------------------*
*     calculate the FAIM intemediates for the 21F term
*---------------------------------------------------------------------*
      IF (.NOT.CC2) THEN 
*---------------------------------------------------------------------*
* calculate the GZeta-breve intermediates for the (E1')' term
* The result is added to the FAIM intermediates:
* GZ-breve to F(A) and GZ-bar-breve to F(QA)
*---------------------------------------------------------------------*
*
* ALL WHAT FOLLOWS FOR LRELAX = TRUE!!!!!!!!!!!
*

*
* GZ-breve to F(A), GZ-barbreve to F(QA)
*
         ISYHZE = MULD2H(ISYCTR,ISYHOP)
         ISYABK = MULD2H(ISY0DIS,ISYMTA)        !al-bet-brevK integral distrib
         ISMGD0 = MULD2H(ISYDEL,ISYCTR)         !eff dens d (for given delta)
         ISQABK = MULD2H(ISY0DIS,ISYHTA)        !alfa-beta-Kbrevbar integral distrib
         ISMGD1 = MULD2H(ISYDEL,ISYHZE)         !eff dens d (for given delta)
*
         IF ((LTWOEL).AND.(LX1ISQ)) THEN
            LDSRHF = MAX(NDSRHF(ISYABK),NDSRHFSQ(ISQABK))
            LSQRUP = .TRUE.
         ELSE
            LDSRHF = MAX(NDSRHF(ISYABK),NDSRHF(ISQABK))
            LSQRUP = .FALSE.
         END IF
         LENMGD = MAX(NT2BGD(ISMGD0),NT2BGD(ISMGD1))

         KDSRHF = 1
         KMGD   = KDSRHF + LDSRHF
         KEND2  = KMGD   + LENMGD
         LWRK2  = LWORK  - KEND2
*
         IF (LWRK2 .LT. 0) THEN
            CALL QUIT(
     &       'Insufficient memory in CCFBTAAO2. (21F and E1_2)')
         END IF
* 
* obs: MO transform 3. index (here called beta) with HOLE matrices
*
* 1.st zero order integral contribution to GZeta-bar
*
         IOPT = 0
                                                    !algam,k^bar(hole)
         CALL CCTRBT2(X0INT,WORK(KDSRHF),XLAMHQA,ISYHTA, WORK(KEND2),
     &                      LWRK2, ISY0DIS,IOPT,.FALSE.,LSQRUP,ONE)
*
* Derivative integral contribution to GZeta-bar (add to previous)
*
         IF (LTWOEL) THEN                       !Husk multiply by -1 for LAOs
            IOPT = 1
            SIGN = DBLE(IREAL)
                                                           !(1)algam,k 
            CALL CCTRBT2(X1INT,WORK(KDSRHF),XLAMDHA,ISYMTA, WORK(KEND2),
     &                        LWRK2, ISY1DIS,IOPT,LX1ISQ,.FALSE.,SIGN)

         END IF

         IADR = IADRE0(IDEL)                       ! d_gami,k density
C
         CALL GETWA2(LUBFDE0,FNBFDE0,WORK(KMGD),IADR,NT2BGD(ISMGD0))
C
         IOPT = 0
         CALL CC_GIM1(WORK(KDSRHF),ISQABK,WORK(KMGD),ISMGD0, FQAIM,
     &                               IOPT,LX1ISQ,WORK(KEND2),LWRK2)
*
* 2.nd zero order integral contribution to both GZeta-bar and GZeta
*
                                                      !algam,k^(hole)
         CALL CCTRBT(X0INT,WORK(KDSRHF),XLAMDHA,ISYMTA,WORK(KEND2),
     &                                               LWRK2,ISY0DIS)
c
         IF (LNEWTA.OR.LNEWZ) THEN
            IOPT = 0
                                                                !GZ-breve
            CALL CC_GIM(WORK(KDSRHF),ISYABK,WORK(KMGD),ISMGD0,FAIM,IOPT,
     &                                                WORK(KEND2),LWRK2)
         END IF

         IADR = IADRE1(IDEL)                             !get d-bar density
         CALL GETWA2(LUBFDE1,FNBFDE1,WORK(KMGD),IADR,NT2BGD(ISMGD1))

         IOPT = 0
         CALL CC_GIM1(WORK(KDSRHF),ISYABK,WORK(KMGD),ISMGD1,FQAIM,
     &                             IOPT, .FALSE.,WORK(KEND2),LWRK2)

*---------------------------------------------------------------------*
*     calculate the calF (FAIM,FQAIM) intermediates for the 21F term 
*---------------------------------------------------------------------*
         
         ISYZT0 = MULD2H(ISYDEL,ISYCTR)    !Symm of P_ci,j (^delta)
         ISYZT1 = MULD2H(ISYZT0,ISYHOP)    !Symm of Pbar_ci,j (^delta)
         ISZTA0 = MULD2H(ISYDEL,ISYZTA)    !Symm of calP_ci,j (^delta)
         ISZTA1 = MULD2H(ISZTA0,ISYHOP)    !Symm of calPbar_ci,j (^delta)

         LEN00   = NT2BCD(ISYZT0)          !P^del
         LEN01   = NT2BCD(ISYZT1)          !P-bar^del
         LENA0   = NT2BCD(ISZTA0)          !calP^del
         LENA1   = NT2BCD(ISZTA1)          !calP-bar^del

c         WRITE (LUPRI,*) 'LEN00   =', LEN00   
c         WRITE (LUPRI,*) 'LEN01   =', LEN01   
c         WRITE (LUPRI,*) 'LENA0   =', LENA0   
c         WRITE (LUPRI,*) 'LENA1   =', LENA1   


         KPINT0 = 1
         KQINT0 = KPINT0 + LEN00
         KPINT1 = KQINT0 + LEN00
         KQINT1 = KPINT1 + LEN01
         KPAIN0 = KQINT1 + LEN01 
         KQAIN0 = KPAIN0 + LENA0
         KPAIN1 = KQAIN0 + LENA0
         KQAIN1 = KPAIN1 + LENA1
         KEND1  = KQAIN1 + LENA1
         LWRK1  = LWORK  - KEND1

c         WRITE (6,*) 'FBTAAO2: Need:',KEND1,'Available:',LWORK
         IF (LWRK1 .LT. 0) THEN
           CALL QUIT('Insufficient work space in CC_FBTA. (4c)')
         END IF
c
c Read in memory the various P^del, Q^del, etc for the
c given delta from files
c
         IADR = IADRPQ0(IDEL)
         CALL GETWA2(LUPQ0,FILPQ0,WORK(KPINT0),IADR,LEN00)
         IADR = IADRPQ0(IDEL) + LEN00
         CALL GETWA2(LUPQ0,FILPQ0,WORK(KQINT0),IADR,LEN00)

         IADR = IADRPQ1(IDEL)
         CALL GETWA2(LUPQ1,FILPQ1,WORK(KPINT1),IADR,LEN01)
         IADR = IADRPQ1(IDEL) + LEN01
         CALL GETWA2(LUPQ1,FILPQ1,WORK(KQINT1),IADR,LEN01)
 
         IADR = IADRPQA0(IDEL)
         CALL GETWA2(LUPQA0,FILPQA0,WORK(KPAIN0),IADR,LENA0)
         IADR = IADRPQA0(IDEL) + LENA0
         CALL GETWA2(LUPQA0,FILPQA0,WORK(KQAIN0),IADR,LENA0)

         IADR = IADRPQA1(IDEL)
         CALL GETWA2(LUPQA1,FILPQA1,WORK(KPAIN1),IADR,LENA1)
         IADR = IADRPQA1(IDEL) + LENA1
         CALL GETWA2(LUPQA1,FILPQA1,WORK(KQAIN1),IADR,LENA1)
c
c  Calculate the F(A) intermediate (depends only on T^A and Zeta)
c  no barred contributions. 
c
C         IF (LNEWTA.OR.LNEWZ) THEN
            IOPT = 2
            CALL CC_21I2(FAIM,X0INT,ISY0DIS,DUMMY,IDUMMY,
     *                WORK(KPINT0),WORK(KQINT0),ISYZT0,
     *                WORK(KPAIN0),WORK(KQAIN0),ISZTA0,
     *                XLAMDP0,XLAMDH0,ISYM0,XLAMDPA,ISYMTA,
     *                WORK(KEND1),LWRK1,IOPT,
     *                .TRUE.,.FALSE.,.FALSE.)

C         END IF
c
c  Calculate the F(QA) intermediate (always, depends on everything)
c
         IOPT = 2
         CALL CC_21I3(FQAIM,X0INT,ISY0DIS,X1INT,ISY1DIS,
     &                WORK(KPINT0),WORK(KQINT0),ISYZT0,
     &                WORK(KPINT1),WORK(KQINT1),ISYZT1,
     &                WORK(KPAIN0),WORK(KQAIN0),ISZTA0,
     &                WORK(KPAIN1),WORK(KQAIN1),ISZTA1,
     &                XLAMDP0,XLAMDH0,ISYM0,XLAMDPQ,ISYHOP,
     &                XLAMDPA,ISYMTA,XLAMPQA,ISYHTA,
     &                WORK(KEND1),LWRK1,IOPT,
     &                .TRUE.,LTWOEL,LX1ISQ)
c        CALL CC_21I3(RHO1, XINT0, ISYDIS0, XINT1, ISYDIS1,
c    *                PINT0,QINT0,ISYPQ0,PINT1,QINT1,ISYPQ1,
c    *                PAINT0,QAINT0,ISYPQA0,PAINT1,QAINT1,ISYPQA1,
c    *                XLAMP0,XLAMH0, ISYML0,XLAMPQ,ISYMLQ,
c    *                XLAMPA,ISYMLA,XLAMPQA,ISYMLQA,
c    *                WORK,LWORK,IOPT,LRHOAO,LDERIV,LXI1SQ)


      END IF
*---------------------------------------------------------------------*
*     That's it; return:
*---------------------------------------------------------------------*
      RETURN
      END 
*=====================================================================*
*                   END OF SUBROUTINE CCFBTAAO2
*=====================================================================*
