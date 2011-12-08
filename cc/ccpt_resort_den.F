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
C  /* Deck resort_daibc */
      SUBROUTINE RESORT_DAIBC(IOPT,ISYDEN,FACT,LU1,FN1,
     &                              LU2,FN2,WORK,LWORK)
C
C     S. Coriani, spring 2002. 
C     Based on CCFOP_SINT by K. Hald
C
C     Purpose : Resort some of the (T) densities with
C               3 virtual indices and 1 occupied.
C     IOPT = 1 (default): resort biac->aibc
C     IOPT = 2, first resort bica->biac, then biac->aibc
C     FACT = controls sign of resorted density
C
      IMPLICIT NONE
C
      INTEGER LU1, LU2, IOPT
      CHARACTER*9 FN1, FN2
      INTEGER ISYDEN, LWORK
      INTEGER ISYMD, ISYCKB, ISYMB, ISYMC, ISYMK, ISYMCK, ISYCKD, ISYMDK
      INTEGER KTRVI1, KTRVI2, KEND1, LWRK1, IOFF1, IOFF2, ISYMBK
      INTEGER KOFF1, KOFF2, KOFF3
C
#if defined (SYS_CRAY)
      REAL WORK(LWORK), HALF, ONE, TWO, FACT
#else
      DOUBLE PRECISION WORK(LWORK), HALF, ONE, TWO, FACT
#endif
C
#include <priunit.h>
#include <ccinftap.h>
#include <ccsdsym.h>
#include <ccorb.h>
C
      PARAMETER (HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      CALL QENTER('RESORT_DAIBC')
C
C-------------------------------
C     Resort from ckbd -> ckdb
!     'resort_daibc: resort bica->biac'
C-------------------------------
C
      if (iopt.eq.2) then

         do isymd = 1, nsym
            isyckb = muld2h(isyden,isymd)
            do isymb = 1, nsym
               isymck = muld2h(isyckb,isymb)
               isyckd = muld2h(isyden,isymb)
C
               ktrvi1 = 1
               kend1  = ktrvi1 + nt1am(isymck)
               lwrk1  = lwork - kend1
C
               if (nt1am(isymck) .gt. 0) then
                  do d = 1, nvir(isymd)
                     do b = 1, nvir(isymb)
C
                        ioff1 = ickbd(isyckb,isymd) 
     &                        + nckatr(isyckb)*(d - 1)
     &                        + ickatr(isymck,isymb)
     &                        + nt1am(isymck)*(b-1)
     &                        + 1
C
                        ioff2 = ickbd(isyckd,isymb) 
     &                        + nckatr(isyckd)*(b - 1)
     &                        + ickatr(isymck,isymd)
     &                        + nt1am(isymck)*(d-1)
     &                        + 1
C
                        call getwa2(lu1,fn1,work(ktrvi1),ioff1,
     &                              nt1am(isymck))
C
                        call putwa2(lu2,fn2,work(ktrvi1),ioff2,
     &                              nt1am(isymck))
C
                     enddo
                  enddo
               endif
            enddo
         enddo
      end if
C
C------------------------------------------------------
C     Sort ckbd to bkcd  (substitute on the same file)
!     'resort_daibc: resort biac->aibc (ciba->bica)'
C------------------------------------------------------
C
      do isymd = 1, nsym
         isyckb = muld2h(isyden,isymd)
         if (nckatr(isyckb) .gt. 0) then
C
            ktrvi1 = 1
            ktrvi2 = ktrvi1 + nckatr(isyckb)
            kend1  = ktrvi2 + nckatr(isyckb)
            lwrk1  = lwork  - kend1
C
            do d = 1, nvir(isymd)
               ioff1 = ickbd(isyckb,isymd) 
     &               + nckatr(isyckb)*(d - 1)
     &               + 1

                call getwa2(lu2,fn2,work(ktrvi1),ioff1,
     &                      nckatr(isyckb))

                !if (fact.ne.one) then
                   call dscal(NCKATR(isyckb),FACT,WORK(ktrvi1),1)
                !end if
C
                do isymb = 1, nsym
                   isymck = muld2h(isyckb,isymb)
                   do isymc = 1, nsym
                      isymk  = muld2h(isymck,isymc)
                      isymbk = muld2h(isymb,isymk)
                      do b = 1, nvir(isymb)
                         do c = 1, nvir(isymc)

                            koff1 = ktrvi1 - 1
     &                            + ickatr(isymck,isymb)
     &                            + nt1am(isymck)*(b-1)
     &                            + it1am(isymc,isymk)
     &                            + c

                            koff2 = ktrvi2 - 1
     &                            + ickatr(isymbk,isymc)
     &                            + nt1am(isymbk)*(c-1)
     &                            + it1am(isymb,isymk)
     &                            + b
C
                            call dcopy(nrhf(isymk),
     &                                 work(koff1),nvir(isymc),
     &                                 work(koff2),nvir(isymb))
C
                         enddo
                     enddo
                   enddo
                enddo
C
                call putwa2(lu2,fn2,work(ktrvi2),ioff1,
     &                      nckatr(isyckb))
C
            enddo
C
         endif
      enddo
C
C-------------
C     End
C-------------
C
      CALL QEXIT('RESORT_DAIBC')
C
      RETURN
C
    1 FORMAT(7X,'Time used in',2X,A12,F12.2,' seconds')
C
      END
