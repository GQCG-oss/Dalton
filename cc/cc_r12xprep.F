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
C  /* Deck r12qv */
      SUBROUTINE R12QV(FROS,IFEL,QQ,S0,X1,Y1,Z1,X2,Y2,Z2,UU,ZZ,
     *                 WORK,LWORK,NBAST,NORBT,NOCCT,
     *                 NVIRT,FEL)
C
C     This subroutine computes the integrals (ak|r12**2|jl).
C
#include <implicit.h>
#include <priunit.h>
#include <ccfro.h>
      DIMENSION QQ(*),X1(*),Y1(*),Z1(*),X2(*),Y2(*),Z2(*),
     *          UU(*),ZZ(*),WORK(LWORK),S0(*)
      LOGICAL FROS,FOUND,FEL,IFEL
      KEND = 1 + NBAST*NBAST
      LWRK = LWORK - KEND
      CALL RDPROP('CM000000',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,S0,WORK(KEND),LWRK,NBAST,NORBT)
      CALL RDPROP('CM010000',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,X1,WORK(KEND),LWRK,NBAST,NORBT)
      CALL RDPROP('CM000100',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,Y1,WORK(KEND),LWRK,NBAST,NORBT)
      CALL RDPROP('CM000001',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,Z1,WORK(KEND),LWRK,NBAST,NORBT)
      CALL RDPROP('CM020000',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,X2,WORK(KEND),LWRK,NBAST,NORBT)
      CALL RDPROP('CM000200',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,Y2,WORK(KEND),LWRK,NBAST,NORBT)
      CALL RDPROP('CM000002',WORK,FOUND)
      CALL UTHZ(WORK,UU,ZZ,Z2,WORK(KEND),LWRK,NBAST,NORBT)
      IF (FROS .EQV. .TRUE.) THEN
         IF (FEL) THEN
            CALL MAKEF(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
         ELSEIF (IFEL) THEN
            CALL MAKF(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
         ELSE 
            CALL MAKEQF(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
         ENDIF
      ELSE 
         IF (IFEL .AND. FEL) THEN
            CALL MAKEQM1(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT)
         ELSEIF (FEL) THEN
            CALL MAKEV(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
         ELSEIF (IFEL) THEN
            CALL MAKV(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
         ELSE 
            CALL MAKEQV(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
         ENDIF
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
C  /* Deck makeqv */
      SUBROUTINE MAKEQV(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
C
C     This subroutine computes the integrals (ak|r12**2|jl).
C
#include <implicit.h>
      PARAMETER (D2 = 2D0)
#include <priunit.h>
      DIMENSION X1(NORBT,NORBT,2),Y1(NORBT,NORBT,2),
     *          Z1(NORBT,NORBT,2),X2(NORBT,NORBT,2),
     *          Y2(NORBT,NORBT,2),Z2(NORBT,NORBT,2),
     *          S0(NORBT,NORBT,2),
     *          QQ(NVIRT,NOCCT,NOCCT,NOCCT)
      DO 100 I = 1, NVIRT 
        DO 200 J = 1, NOCCT
          DO 300 K = 1, NOCCT 
            DO 400 L = 1,NOCCT
              QQ(I,J,K,L) =
     *        (X2(I+NOCCT,L,2) + Y2(I+NOCCT,L,2) + Z2(I+NOCCT,L,2))
     *         * S0(J,K,1) +
     *        (X2(J,K,1) + Y2(J,K,1) +
     *        Z2(J,K,1)) * S0(I+NOCCT,L,2) -
     *        D2 * ( X1(I+NOCCT,L,2) * X1(J,K,1) +
     *               Y1(I+NOCCT,L,2) * Y1(J,K,1) +
     *               Z1(I+NOCCT,L,2) * Z1(J,K,1) )
 400        CONTINUE
 300      CONTINUE
 200    CONTINUE
 100  CONTINUE
      RETURN
      END
*---------------------------------------------------------------------------
c-----------------------------------------------------------------------
C  /* Deck makeqf */
      SUBROUTINE MAKEQF(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
C
C     This subroutine computes the integrals (ak|r12**2|jl).
C
#include <implicit.h>
      PARAMETER (D2 = 2D0)
#include <priunit.h>
      DIMENSION X1(NORBT,NORBT,2),Y1(NORBT,NORBT,2),
     *          Z1(NORBT,NORBT,2),X2(NORBT,NORBT,2),
     *          Y2(NORBT,NORBT,2),Z2(NORBT,NORBT,2),
     *          S0(NORBT,NORBT,2),
     *          QQ(NVIRT,NOCCT,NOCCT,NOCCT)
      DO 100 I = 1, NVIRT
        DO 200 J = 1, NOCCT
          DO 300 K = 1, NOCCT
            DO 400 L = 1,NOCCT
              QQ(I,J,K,L) =
     *        (X2(I,L+NVIRT,2) + Y2(I,L+NVIRT,2) + Z2(I,L+NVIRT,2))
     *         * S0(J+NVIRT,K+NVIRT,1) +
     *        (X2(J+NVIRT,K+NVIRT,1) + Y2(J+NVIRT,K+NVIRT,1) +
     *        Z2(J+NVIRT,K+NVIRT,1)) * S0(I,L+NVIRT,2) -
     *        D2 * ( X1(I,L+NVIRT,2) * X1(J+NVIRT,K+NVIRT,1) +
     *               Y1(I,L+NVIRT,2) * Y1(J+NVIRT,K+NVIRT,1) +
     *               Z1(I,L+NVIRT,2) * Z1(J+NVIRT,K+NVIRT,1) )
 400        CONTINUE
 300      CONTINUE
 200    CONTINUE
 100  CONTINUE
      RETURN
      END
*---------------------------------------------------------------------------

c-----------------------------------------------------------------------
C  /* Deck makeqvi */
      SUBROUTINE MAKEQVI(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
C
C     This subroutine computes the integrals (ak|r12**2|jl).
C
#include <implicit.h>
      PARAMETER (D2 = 2D0)
#include <priunit.h>
      DIMENSION X1(NORBT,NORBT,2),Y1(NORBT,NORBT,2),
     *          Z1(NORBT,NORBT,2),X2(NORBT,NORBT,2),
     *          Y2(NORBT,NORBT,2),Z2(NORBT,NORBT,2),
     *          S0(NORBT,NORBT,2),
     *          QQ(NVIRT,NOCCT,NOCCT,NOCCT)
      DO 100 I = 1, NVIRT
        DO 200 J = 2, NOCCT
          DO 300 K = 2, NOCCT
            DO 400 L = 2,NOCCT
              QQ(I,J,K,L) =
     *        (X2(I+NOCCT+1,L,2) + Y2(I+NOCCT+1,L,2)+ Z2(I+NOCCT+1,L,2))
     *         * S0(J,K,1) +
     *        (X2(J,K,1) + Y2(J,K,1) +
     *        Z2(J,K,1)) * S0(I+NOCCT+1,L,2) -
     *        D2 * ( X1(I+NOCCT+1,L,2) * X1(J,K,1) +
     *               Y1(I+NOCCT+1,L,2) * Y1(J,K,1) +
     *               Z1(I+NOCCT+1,L,2) * Z1(J,K,1) )
 400        CONTINUE
 300      CONTINUE
 200    CONTINUE
 100  CONTINUE
      RETURN
      END
*---------------------------------------------------------------------------
C  /* Deck makev */
      SUBROUTINE MAKEV(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
C
C     This subroutine computes the integrals (ak|r12**2|jl).
C
#include <implicit.h>
      PARAMETER (D2 = 2D0)
#include <priunit.h>
      DIMENSION X1(NORBT,NORBT,2),Y1(NORBT,NORBT,2),
     *          Z1(NORBT,NORBT,2),X2(NORBT,NORBT,2),
     *          Y2(NORBT,NORBT,2),Z2(NORBT,NORBT,2),
     *          S0(NORBT,NORBT,2),
     *          QQ(NVIRT,NOCCT,NOCCT,NOCCT)
      DO I = 1,NOCCT
        DO J = 1,NVIRT
          DO K = 1,NOCCT 
            DO L = 1,NOCCT
              QQ(J,I,L,K) =
     *        (X2(I,L,2) + Y2(I,L,2) + Z2(I,L,2))
     *         * S0(J+NOCCT,K,1) +
     *        (X2(J+NOCCT,K,1) + Y2(J+NOCCT,K,1) +
     *        Z2(J+NOCCT,K,1)) * S0(I,L,2) -
     *        D2 * ( X1(I,L,2) * X1(J+NOCCT,K,1) +
     *               Y1(I,L,2) * Y1(J+NOCCT,K,1) +
     *               Z1(I,L,2) * Z1(J+NOCCT,K,1) )
            ENDDO   
          ENDDO    
        ENDDO   
      ENDDO   
      RETURN
      END
*---------------------------------------------------------------------------
*---------------------------------------------------------------------------
C  /* Deck makef */
      SUBROUTINE MAKEF(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
C
C     This subroutine computes the integrals (ak|r12**2|jl).
C
#include <implicit.h>
      PARAMETER (D2 = 2D0)
#include <priunit.h>
      DIMENSION X1(NORBT,NORBT,2),Y1(NORBT,NORBT,2),
     *          Z1(NORBT,NORBT,2),X2(NORBT,NORBT,2),
     *          Y2(NORBT,NORBT,2),Z2(NORBT,NORBT,2),
     *          S0(NORBT,NORBT,2),
     *          QQ(NVIRT,NOCCT,NOCCT,NOCCT)
      DO I = 1,NOCCT
        DO J = 1,NVIRT
          DO K = 1,NOCCT
            DO L = 1,NOCCT
              QQ(J,I,L,K) =
     *        (X2(I+NVIRT,L+NVIRT,2) + Y2(I+NVIRT,L+NVIRT,2) +
     &         Z2(I+NVIRT,L+NVIRT,2))
     *         * S0(J,K+NVIRT,1) +
     *        (X2(J,K+NVIRT,1) + Y2(J,K+NVIRT,1) +
     *        Z2(J,K+NVIRT,1)) * S0(I+NVIRT,L+NVIRT,2) -
     *        D2 * ( X1(I+NVIRT,L+NVIRT,2) * X1(J,K+NVIRT,1) +
     *               Y1(I+NVIRT,L+NVIRT,2) * Y1(J,K+NVIRT,1) +
     *               Z1(I+NVIRT,L+NVIRT,2) * Z1(J,K+NVIRT,1) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
*---------------------------------------------------------------------------

*---------------------------------------------------------------------------
C  /* Deck makv */
      SUBROUTINE MAKV(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
C
C     This subroutine computes the integrals (ak|r12**2|jl).
C
#include <implicit.h>
      PARAMETER (D2 = 2D0)
#include <priunit.h>
      DIMENSION X1(NORBT,NORBT,2),Y1(NORBT,NORBT,2),
     *          Z1(NORBT,NORBT,2),X2(NORBT,NORBT,2),
     *          Y2(NORBT,NORBT,2),Z2(NORBT,NORBT,2),
     *          S0(NORBT,NORBT,2),
     *          QQ(NVIRT,NOCCT,NOCCT,NOCCT)
      DO I = 1,NOCCT
        DO J = 1,NOCCT
          DO K = 1,NVIRT
            DO L = 1,NOCCT
              QQ(K,L,J,I) =
     *        (X2(I,K+NOCCT,2) + Y2(I,K+NOCCT,2) + Z2(I,K+NOCCT,2))
     *         * S0(J,L,1) +
     *        (X2(J,L,1) + Y2(J,L,1) +
     *        Z2(J,L,1)) * S0(I,K+NOCCT,2) -
     *        D2 * ( X1(I,K+NOCCT,2) * X1(J,L,1) +
     *               Y1(I,K+NOCCT,2) * Y1(J,L,1) +
     *               Z1(I,K+NOCCT,2) * Z1(J,L,1) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END
*-------------------------------------------------------------------
*=====================================================================*
*---------------------------------------------------------------------------
C  /* Deck makf */
      SUBROUTINE MAKF(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT,NVIRT)
C
C     This subroutine computes the integrals (ak|r12**2|jl).
C
#include <implicit.h>
      PARAMETER (D2 = 2D0)
#include <priunit.h>
      DIMENSION X1(NORBT,NORBT,2),Y1(NORBT,NORBT,2),
     *          Z1(NORBT,NORBT,2),X2(NORBT,NORBT,2),
     *          Y2(NORBT,NORBT,2),Z2(NORBT,NORBT,2),
     *          S0(NORBT,NORBT,2),
     *          QQ(NVIRT,NOCCT,NOCCT,NOCCT)
      DO I = 1,NOCCT
        DO J = 1,NOCCT
          DO K = 1,NVIRT
            DO L = 1,NOCCT
              QQ(K,L,J,I) =
     *        (X2(I+NVIRT,K,2) + Y2(I+NVIRT,K,2) + Z2(I+NVIRT,K,2))
     *         * S0(J+NVIRT,L+NVIRT,1) +
     *        (X2(J+NVIRT,L+NVIRT,1) + Y2(J+NVIRT,L+NVIRT,1) +
     *        Z2(J+NVIRT,L+NVIRT,1)) * S0(I+NVIRT,K,2) -
     *        D2 * ( X1(I+NVIRT,K,2) * X1(J+NVIRT,L+NVIRT,1) +
     *               Y1(I+NVIRT,K,2) * Y1(J+NVIRT,L+NVIRT,1) +
     *               Z1(I+NVIRT,K,2) * Z1(J+NVIRT,L+NVIRT,1) )
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      RETURN
      END
*-------------------------------------------------------------------
*=====================================================================*

      subroutine cc_r12xsort(sajkl,blajkl,nvir0)
c---------------------------------------------------------------------

      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccfro.h>
      integer isymkl,isymij,isyml,isymk,isymaj
      integer isymj,isyma,idxkl,isym,nvir0(8)
      integer kt,idx,nvir0t,koffj,koffa,koffk,koffl
      integer icoun1,kl,aj,ivir0(8),ajkl,nrhfst
#if defined (SYS_CRAY)
      real blajkl(*),sajkl(*)
#else
      double precision blajkl(*),sajkl(*)
#endif
c
      call qenter('cc_r12xsort')
c
      nvir0t = 0
      nrhfst = 0
      do isym = 1, nsym
         nvir0t      = nvir0t + nvir0(isym)
         nrhfst      = nrhfst + nrhfs(isym)
      end do

      icoun1 = 0
      do isym = 1,nsym
         ivir0(isym) = icoun1
         icoun1 = icoun1 + nvir0(isym)
      enddo   

      
      idx = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhfs(isymk)
              koffk = irhfs(isymk)+k
             do l = 1, nrhfs(isyml)
                koffl = irhfs(isyml)+l
                idxkl=imatij(isymk,isyml)+nrhfs(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    do j = 1, nrhfs(isymj)
                       koffj = irhfs(isymj)+j
                       do a = 1, nvir0(isyma)
                          idx = idx + 1
                          koffa = ivir0(isyma)+a
                          kl = (koffk-1)*nrhfst + koffl
                          aj = (koffj-1)*nvir0t + koffa
                          ajkl = aj +(kl-1)*nvir0t*nrhfst
                          blajkl(idx) = sajkl(ajkl)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
         enddo
      enddo

      call qexit('cc_r12xsort')
      end

*---------------------------------------------------------------------------
*=====================================================================*
      subroutine cc_r12xsort1(sajkl,blajkl)
c---------------------------------------------------------------------

      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccfro.h>
      integer isymkl,isymij,isyml,isymk,isymaj
      integer isymj,isyma,idxkl,isym,nvir0(8)
      integer kt,idx,nvir0t,koffj,koffa,koffk,koffl
      integer icoun1,kl,aj,ivir0(8),ajkl,nrhfst
#if defined (SYS_CRAY)
      real blajkl(*),sajkl(*)
#else
      double precision blajkl(*),sajkl(*)
#endif
c
      call qenter('cc_r12xsort1')
c
      nvir0t = 0
      nrhfst = 0
      do isym = 1, nsym
         nvir0(isym) = norb1(isym) - nrhfs(isym) 
         nvir0t      = nvir0t + nvir0(isym)
         nrhfst      = nrhfst + nrhfs(isym)
      end do

      icoun1 = 0
      do isym = 1,nsym
         ivir0(isym) = icoun1
         icoun1 = icoun1 + nvir0(isym)
      enddo   

      
      idx = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhfs(isymk)
              koffk = irhfs(isymk)+k
             do l = 1, nrhfs(isyml)
                koffl = irhfs(isyml)+l
                idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    do j = 1, nrhfs(isymj)
                       koffj = irhfs(isymj)+j
                       do a = 1, nrhfs(isyma)
                          koffa = irhfs(isyma)+a
                          idx = idx + 1
                          kl = (koffl-1)*nrhfst + koffk
                          aj = (koffj-1)*nrhfst + koffa
                          ajkl = aj +(kl-1)*nrhfst*nrhfst
                          blajkl(idx) = sajkl(ajkl)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
         enddo
      enddo

      call qexit('cc_r12xsort1')
      end

*---------------------------------------------------------------------------
*=====================================================================*
C  /* Deck utz */
      SUBROUTINE UTZ(AA,UU,ZZ,BB,WORK,LWORK,NBAST,NOCCT)
C
C              B(1) = Z(Transpose) * A * Z
C              B(2) =  A * Z
C
C     On input, AA is the upper triangle of the symmetric
C     matrix A. B is returned through the array BB as
C     square matrix.
C
c     Modified UTHZ
#include <implicit.h>
#include <priunit.h>
      DIMENSION AA(*),UU(NBAST,NOCCT),ZZ(NBAST,NOCCT),
     *          WORK(NBAST,NOCCT),BB(NOCCT,NOCCT,2)
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
      IF (NBAST*NOCCT .GT. LWORK)
     *      CALL STOPIT('UTZ',' ',LWORK,NBAST*NOCCT)
      CALL DZERO(WORK,NBAST*NOCCT)
      DO  K = 1, NBAST
         DO  J = 1, NOCCT
            ZZKJ = ZZ(K,J)
            DO  I = 1, NBAST
               IK = INDEX(I,K)
               WORK(I,J) = WORK(I,J) + AA(IK) * ZZKJ
            ENDDO   
         ENDDO   
      ENDDO   
      DO I = 1, NOCCT
         DO J = 1, NOCCT
            BB(I,J,1) = DDOT(NBAST,ZZ(1,I),1,WORK(1,J),1)
         ENDDO   
      ENDDO   
      DO I = 1, NBAST
         DO J = 1, NOCCT
             BB(I,J,2) = WORK(I,J)
         ENDDO
      ENDDO
      RETURN
      END
*=====================================================================*
c---------------------------------------------------------------------------
c--------------------------------------------------------------------------
C  /* Deck makeqm1 */
      SUBROUTINE MAKEQM1(S0,X1,Y1,Z1,X2,Y2,Z2,QQ,NOCCT,NORBT)
C
C     This subroutine computes the integrals (i'k|r12**2|jl).
C
C     Written by Wim Klopper (University of Karlsruhe, 31 October 2002).
C
#include <implicit.h>
      PARAMETER (D2 = 2D0)
#include <priunit.h>
      DIMENSION X1(NORBT,NORBT,2),Y1(NORBT,NORBT,2),
     *          Z1(NORBT,NORBT,2),X2(NORBT,NORBT,2),
     *          Y2(NORBT,NORBT,2),Z2(NORBT,NORBT,2),
     *          S0(NORBT,NORBT,2),
     *          QQ(NORBT,NOCCT,NOCCT,NOCCT)
      DO 100 I = 1,NORBT 
        DO 200 J = 1, NOCCT
          DO 300 K = 1, NOCCT
            DO 400 L = 1,NOCCT
              QQ(I,J,K,L) =
     *        (X2(I,L,2) + Y2(I,L,2) + Z2(I,L,2))
     *         * S0(J,K,1) +
     *        (X2(J,K,1) + Y2(J,K,1) +
     *        Z2(J,K,1)) * S0(I,L,2) -
     *        D2 * ( X1(I,L,2) * X1(J,K,1) +
     *               Y1(I,L,2) * Y1(J,K,1) +
     *               Z1(I,L,2) * Z1(J,K,1) )
 400        CONTINUE
 300      CONTINUE
 200    CONTINUE
 100  CONTINUE
      RETURN
      END
c---------------------------------------------------------------------------
*=====================================================================*
*=====================================================================*
      subroutine cc_r12xsort2(sajkl,blajkl,nvir0)
c---------------------------------------------------------------------

      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
      integer isymkl,isymij,isyml,isymk,isymaj
      integer isymj,isyma,idxkl,isym
      integer nvir0(8),kt,idx,koffj,koffa,koffk,koffl
      integer nvir0t,icoun1,kl,aj,iorb1(8),ajkl,ivir0(8)
#if defined (SYS_CRAY)
      real blajkl(*),sajkl(*)
#else
      double precision blajkl(*),sajkl(*)
#endif
c
      call qenter('cc_r12xsort2')

c
      nvir0t = 0
      do isym = 1, nsym
         nvir0t      = nvir0t + nvir0(isym)
      end do

      icoun1 = 0
      do isym = 1,nsym
         ivir0(isym) = icoun1
         icoun1 = icoun1 + nvir0(isym)
      enddo


      idx = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhf(isymk)
              koffk = irhf(isymk)+k
             do l = 1, nrhf(isyml)
                koffl = irhf(isyml)+l
                idxkl=imatij(isymk,isyml)+nrhfs(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    do j = 1, nrhf(isymj)
                       koffj = irhf(isymj)+j
                       do a = 1, nvir0(isyma)
                          idx = idx + 1
                          koffa = ivir0(isyma)+a
                          kl = (koffk-1)*nrhft + koffl
                          aj = (koffj-1)*nvir0t + koffa
                          ajkl = aj +(kl-1)*nvir0t*nrhft
                          blajkl(idx) = sajkl(ajkl)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
         enddo
      enddo


c

      call qexit('cc_r12xsort2')
      end

*---------------------------------------------------------------------------
*-------------------------------------------------------------------
*=====================================================================*
      subroutine cc_r12xsort3(sajkl,blajkl,nvir0)
c---------------------------------------------------------------------

      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccfro.h>
      integer isymkl,isymij,isyml,isymk,isymaj
      integer isymj,isyma,idxkl,isym,nvir0(8)
      integer kt,idx,nvir0t,koffj,koffa,koffk,koffl
      integer icoun1,kl,aj,ivir0(8),ajkl
#if defined (SYS_CRAY)
      real blajkl(*),sajkl(*)
#else
      double precision blajkl(*),sajkl(*)
#endif
c
      call qenter('cc_r12xsort3')
c
      nvir0t = 0
      do isym = 1, nsym
         nvir0t      = nvir0t + nvir0(isym)
      end do

      icoun1 = 0
      do isym = 1,nsym
         ivir0(isym) = icoun1
         icoun1 = icoun1 + nvir0(isym)
      enddo


      idx = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhf(isymk)
              koffk = irhf(isymk)+k
             do l = 1, nrhf(isyml)
                koffl = irhf(isyml)+l
                idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    do j = 1, nrhf(isymj)
                       koffj = irhf(isymj)+j
                       do a = 1, nvir0(isyma)
                          idx = idx + 1
                          koffa = ivir0(isyma)+a
                          kl = (koffk-1)*nrhft + koffj
                          aj = (koffl-1)*nvir0t + koffa
                          ajkl = aj +(kl-1)*nvir0t*nrhft
                          blajkl(idx) = sajkl(ajkl)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
         enddo
      enddo

      call qexit('cc_r12xsort3')
      end

*=====================================================================*
      subroutine cc_r12xsort4(sajkl,blajkl)
c---------------------------------------------------------------------

      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccfro.h>
      integer isymkl,isymij,isyml,isymk,isymaj
      integer irhffr(8),isymj,isyma,idxkl,isym,nvir0(8)
      integer kt,idx,nvir0t,koffj,koffa,koffk,koffl
      integer nrhffrt,nrhfst,icoun1,kl,aj,ivir0(8),ajkl
#if defined (SYS_CRAY)
      real blajkl(*),sajkl(*)
#else
      double precision blajkl(*),sajkl(*)
#endif
c
      call qenter('cc_r12xsort4')
c
      nvir0t = 0
      nrhfst= 0
      nrhffrt= 0
      do isym = 1, nsym
         nvir0(isym) = norb1(isym) - nrhfs(isym)
         nvir0t      = nvir0t + nvir0(isym)
         nrhfst      = nrhfst + nrhfs(isym)
         nrhffrt      = nrhffrt + nrhffr(isym)
      end do

      icoun1 = 0
      do isym = 1,nsym
         irhffr(isym) = icoun1
         icoun1 = icoun1 + nrhffr(isym)
      enddo


      idx = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhfs(isymk)
              koffk = irhfs(isymk)+k
             do l = 1, nrhfs(isyml)
                koffl = irhfs(isyml)+l
                idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    do j = 1, nrhfs(isymj)
                       koffj = irhfs(isymj)+j
                       do a = 1, nrhfs(isyma)
                          idx = idx + 1
                          koffa = irhfs(isyma)+a
                          kl = (koffk-1)*nrhfst + koffl
                          aj = (koffj-1)*nrhfst + koffa
                          ajkl = aj +(kl-1)*nrhfst*nrhfst
                          blajkl(idx) = sajkl(ajkl)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
         enddo
      enddo
      call qexit('cc_r12xsort4')
      end

*=====================================================================*
      subroutine cc_r12fsym(imatmunu,sajkl,blajkl)
c---------------------------------------------------------------------

      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
      integer isymkl,isymij,isyml,isymk,isymaj
      integer isymj,isyma,idxkl,isym
      integer kt,idx,koffj,koffa,koffk,koffl
      integer ibasi(8),icoun1,kl,aj,iorb1(8),ajkl
      integer imatmunu(8,8)
#if defined (SYS_CRAY)
      real blajkl(*),sajkl(*)
#else
      double precision blajkl(*),sajkl(*)
#endif
c
      call qenter('cc_r12fsym')
c
      nbast = 0
      do isym = 1, nsym
         nbast      = nbast + nbas(isym)
      end do

      icoun1 = 0
      do isym = 1,nsym
         ibasi(isym) = icoun1
         icoun1 = icoun1 + nbas(isym)
      enddo

      idx = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nbas(isymk)
               koffk = ibasi(isymk)+k
               do l = 1, nbas(isyml)
                  koffl = ibasi(isyml)+l
                  idxkl=imatmunu(isyml,isymk)+nbas(isyml)*(k-1)+l
                  idx = idx + 1
                  kl = (koffk-1)*nbast + koffl
                  blajkl(kl) = sajkl(idx)
               enddo
            enddo
         enddo
      enddo


      call qexit('cc_r12fsym')
      end

*---------------------------------------------------------------------------
*=====================================================================*
      subroutine cc_r12mksfr(V2AM,F2AM,WORK,LWRK1)
c---------------------------------------------------------------------
c     Make s_{kl}^{Ij*}
c---------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccfro.h>
#include <ccsdinp.h>
#if defined (SYS_CRAY)
      real zero,one
      real work(*),v2am(*),f2am(*)
#else
      double precision zero,one
      double precision work(*),v2am(*),f2am(*)
#endif
      logical ldum,test,fel
      parameter (test = .false.,zero = 0.0d0, one = 1.0d0)
      integer isymai,ISYMIP,ISYMI,ISYMP,ISYMK,NPK,KGX,npi,nki
      integer npiki,kug,kcmo,kend1,LUSIRG,index,LWRK1,isymik
      integer isym,it1vm(8,8),it2vm(8,8)
      integer nrhffrt,isyma,lunit,IGX,IUG,nb1fro(8)
      integer ISMRS0,ISMRX1,ISMRY1,ISMRZ1,ISMRX2,ISMRY2,ISMRZ2
      integer nrhfst,ntotgu,kemo,lu43,kugz,kcmoz,LQQ,KQF,IEND
      integer kend2,norb1t,ismo0,ismo2,ismo1,ismou,isymj
      integer iv2fro(8,8),iv1fro(8,8),icoun1,icoun2,icoun3,icoun4
      integer kdmo,i2frofr(8,8),icoun5
      integer kemoz,kdmox,kdmoz,isymkl,isymaj,isyml

      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J

      call qenter('cc_r12mksfr')

         do isymk = 1,nsym
            nv1fro(isymk) = 0
            do isymj = 1,nsym
               isymi  = muld2h(isymk,isymj)
               nv1fro(isymk) = nv1fro(isymk) + nrhffr(isymj)*
     &                          (norb1(isymi)-nrhffr(isymi))
            enddo
         enddo

 
         do isymk = 1, nsym
            icoun1 = 0
            icoun2 = 0
            icoun3 = 0
            icoun4 = 0
            icoun5 = 0 
            do isymj = 1,nsym
               isymi  = muld2h(isymk,isymj)
               iv2fro(isymi,isymj) = icoun3
               iv1fro(isymi,isymj) = icoun4
               i2frofr(isymi,isymj) = icoun5
               icoun3 = icoun3 + nt1fro(isymi)*nv1fro(isymj)
               icoun4 = icoun4 + (nrhffr(isymi))
     &                   *(norb1(isymj)-nrhffr(isymj))
               icoun5 = icoun5 + nt1fro(isymi)*nfrofr(isymj)
            enddo
         enddo


      DO ISYMAI = 1,NSYM
         Nb1FRO(ISYMAI) = 0
         DO ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            Nb1FRO(ISYMAI) = Nb1FRO(ISYMAI) + NRHFFR(ISYMA)*NBAS(ISYMI)
         ENDDO   
      ENDDO   

      NRHFFRT = 0
      NRHFST = 0
      NORB1T = 0
      DO ISYM = 1,NSYM
         NRHFFRT = NRHFFRT  +  NRHFFR(ISYM)
         NRHFST = NRHFST  +  NRHFS(ISYM)
         NORB1T = NORB1T  +  NORB1(ISYM)
      ENDDO


      kgx   = 1
      kcmo  = kgx + nt1fro(1) 
      kug   = kcmo+ nlamds   
      kemo  = kug + nb1fro(1)
      kdmo  = kemo+ nlamds 
      kend1 = kdmo+ nlamds 

      call dzero(work(kemo),nlamds)
      lu43 = -1
      call gpopen(lu43,'GUMAT.1','UNKNOWN',' ','UNFORMATTED',
     &                      IDUMMY,LDUM)
      REWIND(LU43)
      READ(LU43) NTOTGU
      READ(LU43) (work(kemo+I-1), I = 1, NTOTGU)
      CALL GPCLOSE(LU43,'KEEP')





      kcmoz = kcmo
      kugz  = kug

      call dzero(work(kcmo),nlamds)
      LUSIRG = -1
      CALL GPOPEN(LUSIRG,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,LDUM)
      REWIND LUSIRG

      CALL MOLLAB(LABEL,LUSIRG,LUPRI)
      READ (LUSIRG)
C
      READ (LUSIRG)
      READ (LUSIRG) (WORK(kcmo+i-1), I=1,NLAMDS)
C
      CALL GPCLOSE(LUSIRG,'KEEP')


      call dzero(work(kGX),nt1fro(1))
      call dzero(work(kug),nt1fro(1))

      DO ISYMIP = 1, NSYM
         ISYMIK = ISYMIP
         DO ISYMI = 1, NSYM
            ISYMP = MULD2H(ISYMI,ISYMIP)
            ISYMK = ISYMP
            DO I = 1, NRHF(ISYMI)
               NPK = IT1FRO(ISYMP,ISYMK)
               DO K = 1, NRHFFR(ISYMK)
                  DO P = 1, NVIR(ISYMP)
                     NPK = NPK + 1
                     NPI = IT1AM(ISYMP,ISYMI)
     *                   + NVIR(ISYMP)*(I-1) + P
                     NKI = IT1AM(ISYMK,ISYMI)
     *                   + NVIR(ISYMK)*(I-1) + K
                     NPIKI = IT2AM(ISYMIP,ISYMIK)
     *                     + INDEX(NPI,NKI)
                     work(kGX+NPK-1) = work(kGX+NPK-1) + V2AM(NPIKI)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO ISYMIP = 1, NSYM
         ISYMIK = ISYMIP
         DO ISYMI = 1, NSYM
            ISYMP = MULD2H(ISYMI,ISYMIP)
            ISYMK = ISYMP
            DO I = 1, NRHFFR(ISYMI)
               NPK = IT1FRO(ISYMP,ISYMK)
               DO K = 1, NRHFFR(ISYMK)
                  DO P = 1, NVIR(ISYMP)
                     NPK = NPK + 1
                     NKI = IFROFR(ISYMI,ISYMK)
     *                   + NRHFFR(ISYMI)*(K-1) + I
                     NPI = IT1FRO(ISYMP,ISYMI)
     *                   + NVIR(ISYMP)*(I-1) + P
                     NPIKI = i2frofr(ISYMIP,ISYMIK)
     *                     + NT1FRO(ISYMIP)*(NKI-1)
     *                     + NPI
                     work(kGX+NPK-1) = work(kGX+NPK-1) + F2AM(NPIKI)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      DO ISYM = 1, NSYM
         NPK = IT1fro(ISYM,ISYM) + NORB1(ISYM) + 1
         DO K = 1, NRHFFR(ISYM)
            CALL DZERO(work(kGX+NPK-1),NORB2(ISYM))
            NPK = NPK + NVIR(ISYM)
         ENDDO
      ENDDO   


      IGX = 1
      IUG = 1
      DO ISYM = 1, NSYM
         KCMO = KCMO + NRHFS(ISYM)*NBAS(ISYM)
         CALL MPAB(WORK(KCMO),NBAS(ISYM),NVIR(ISYM),
     *                         NBAS(ISYM),NVIR(ISYM),
     *             work(kGX+IGX-1),NVIR(ISYM),NRHFFR(ISYM),
     *                     NVIR(ISYM),NRHFFR(ISYM),
     *             work(kUG+IUG-1),NBAS(ISYM),NRHFFR(ISYM))
         IF (NORB1(ISYM)*NBAS(ISYM) .NE. 0 .AND. IPRINT .GT.3) THEN
            WRITE(LUPRI,'(/A,I2)') ' U(NRHFFR) x X/ISYM =',ISYM
            CALL OUTPUT(work(kUG+IUG-1),1,NBAS(ISYM),1,NRHFFR(ISYM),
     *                  NBAS(ISYM),NRHFFR(ISYM),1,LUPRI)
            WRITE(LUPRI,'(/A,I2)') ' G(NRHFFR) x X/ISYM =',ISYM
            CALL OUTPUT(work(kGX+IGX-1),1,NBAS(ISYM),1,NRHFFR(ISYM),
     *                  NBAS(ISYM),NRHFFR(ISYM),1,LUPRI)
         ENDIF 
         KCMO  = KCMO  + NBAS(ISYM)*NVIRS(ISYM)
         IGX = IGX + NVIR(ISYM)*NRHFFR(ISYM)
         IUG = IUG + NBAS(ISYM)*NRHFFR(ISYM)
      ENDDO

      kdmoz = kdmo
      kdmox = kdmo
      call dzero(work(kdmo),nlamds)
      do isym =1, nsym
         do i = 1,nrhffr(isym)
            call dcopy(nbas(isym),work(kugz),1,work(kdmo),1)
            kdmo   = kdmo   + nbas(isym)
            kugz   = kugz   + nbas(isym)
         enddo
         kdmo = kdmo + nrhf(isym)*nbas(isym)   
      enddo
      kemoz = kemo
      do isym =1, nsym
         kdmox = kdmox + nrhffr(isym)*nbas(isym)
         do i = 1,nrhf(isym)
            call dcopy(nbas(isym),work(kemoz),1,work(kdmox),1)
            kdmox   = kdmox   + nbas(isym)
            kemoz   = kemoz   + nbas(isym)
         enddo
         kemoz = kemoz + nbas(isym)*(norb1(isym)-nrhfs(isym))
      enddo
      
      

      ISMO0 = 1 + nt1fro(1) 
      ISMOU = 1
      ISMO1 = KEND1 + NLAMDS
      ISMO2 = ISMO1 + NORB1T*NBAST
      KEND2 = ISMO2 + NLAMDS
      CALL DZERO(WORK(ISMO1),NLAMDS)
      CALL DZERO(WORK(ISMO2),NLAMDS)
      DO ISYM=1,NSYM
         ISMO0 = ISMO0 + NRHFS(ISYM)*NBAS(ISYM) 
         DO I = 1, NRHFFR(ISYM)
           CALL DCOPY(NBAS(ISYM),WORK(ISMO0),1,WORK(ISMO1),1)
           CALL DCOPY(NBAS(ISYM),work(kug+ISMOU-1),1,WORK(ISMO2),1)
           ISMO1 = ISMO1 + NBAST
           ISMO2 = ISMO2 + NBAST
           ISMO0 = ISMO0 + NBAS(ISYM)
           ISMOU = ISMOU + NBAS(ISYM)
         ENDDO   
         ISMO0 = ISMO0 + 
     &           (NORB2(ISYM)+NORB1(ISYM)-NRHFFR(ISYM))*NBAS(ISYM)
         ISMO1 = ISMO1 + NBAS(ISYM)
         ISMO2 = ISMO2 + NBAS(ISYM)
      ENDDO   
      ISMO0 = 1 + nt1fro(1)
      ISMOU = 1
      ISMO1 = KEND1 + NLAMDS+NRHFFRT*NBAST
      ISMO2 = ISMO1 + NORB1T*NBAST
      KEND2 = ISMO2 + NLAMDS
      DO ISYM=1,NSYM
         ISMO0 = ISMO0 + (NRHFS(ISYM)+NRHFFR(ISYM))*NBAS(ISYM)
         DO I = 1, NRHF(ISYM)
           CALL DCOPY(NBAS(ISYM),WORK(ISMO0),1,WORK(ISMO1),1)
           CALL DCOPY(NBAS(ISYM),work(kemo+ISMOU-1),1,WORK(ISMO2),1)
           ISMO1 = ISMO1 + NBAST
           ISMO2 = ISMO2 + NBAST
           ISMO0 = ISMO0 + NBAS(ISYM)
           ISMOU = ISMOU + NBAS(ISYM)
         ENDDO
         ISMO0 = ISMO0 +
     &           (NORB2(ISYM)+NORB1(ISYM)-NRHFS(ISYM))*NBAS(ISYM)
         ISMOU=ISMOU + (NORB1(ISYM)-NRHFS(ISYM))*NBAS(ISYM)
         ISMO1 = ISMO1 + NBAS(ISYM)
         ISMO2 = ISMO2 + NBAS(ISYM)
      ENDDO


      LUNIT = -1
      CALL GPOPEN(LUNIT,'EXFRO','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUM)
      WRITE(LUNIT,'(4E30.20)') (work(kUG+J-1), J = 1,NB1FRO(1)) 
      CALL GPCLOSE(LUNIT,'KEEP')


      ISMO1   = KEND1 + NLAMDS
      ISMO2   = ISMO1 + NORB1T*NBAST 
      ISMRS0  = KEND2
      ISMRX1  = ISMRS0  + 2*NRHFST*NRHFST
      ISMRY1  = ISMRX1  + 2*NRHFST*NRHFST
      ISMRZ1  = ISMRY1  + 2*NRHFST*NRHFST
      ISMRX2  = ISMRZ1  + 2*NRHFST*NRHFST
      ISMRY2  = ISMRX2  + 2*NRHFST*NRHFST
      ISMRZ2  = ISMRY2  + 2*NRHFST*NRHFST
      KQF     = ISMRZ2  + 2*NRHFST*NRHFST
      IEND    = KQF     + NRHFFRT*NRHFT**3
      LQQ   = LWRK1 - IEND

      

      CALL DZERO(WORK(ISMRS0),14*NRHFST*NRHFST)
      CALL DZERO(WORK(KQF),NRHFFRT*NRHFT**3)
      CALL R12QV(.TRUE.,.FALSE.,WORK(KQF),WORK(ISMRS0),WORK(ISMRX1),
     *           WORK(ISMRY1),
     *           WORK(ISMRZ1),
     *           WORK(ISMRX2),WORK(ISMRY2),
     *           WORK(ISMRZ2),
     *           WORK(ISMO2),WORK(ISMO1),
     *           WORK(IEND),LQQ,NBAST,NRHFST,NRHFT,NRHFFRT,.FALSE.)
      LUNIT = -1
      CALL GPOPEN(LUNIT,'AUXQSF12','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUM)
      WRITE(LUNIT,'(4E30.20)') (WORK(KQF+J-1), J = 1, NRHFFRT*NRHFT**3)
      CALL GPCLOSE(LUNIT,'KEEP')


      CALL DZERO(WORK(ISMRS0),14*NRHFST*NRHFST)
      CALL DZERO(WORK(KQF),NRHFFRT*NRHFT**3)
      CALL R12QV(.TRUE.,.FALSE.,WORK(KQF),WORK(ISMRS0),WORK(ISMRX1),
     *           WORK(ISMRY1),
     *           WORK(ISMRZ1),
     *           WORK(ISMRX2),WORK(ISMRY2),
     *           WORK(ISMRZ2),
     *           WORK(ISMO2),WORK(ISMO1),
     *           WORK(IEND),LQQ,NBAST,NRHFST,NRHFT,NRHFFRT,.TRUE.)
      LUNIT = -1
      CALL GPOPEN(LUNIT,'AUXQSFJ12','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUM)
      WRITE(LUNIT,'(4E30.20)') (WORK(KQF+J-1), J = 1, NRHFFRT*NRHFT**3)
      CALL GPCLOSE(LUNIT,'KEEP')

      CALL DZERO(WORK(ISMRS0),14*NRHFST*NRHFST)
      CALL DZERO(WORK(KQF),NRHFFRT*NRHFT**3)
      CALL R12QV(.TRUE.,.TRUE.,WORK(KQF),WORK(ISMRS0),WORK(ISMRX1),
     *           WORK(ISMRY1),
     *           WORK(ISMRZ1),
     *           WORK(ISMRX2),WORK(ISMRY2),
     *           WORK(ISMRZ2),
     *           WORK(ISMO2),WORK(ISMO1),
     *           WORK(IEND),LQQ,NBAST,NRHFST,NRHFT,NRHFFRT,.FALSE.)
      LUNIT = -1
      CALL GPOPEN(LUNIT,'AUXQSFI12','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUM)
      WRITE(LUNIT,'(4E30.20)') (WORK(KQF+J-1), J = 1, NRHFFRT*NRHFT**3)
      CALL GPCLOSE(LUNIT,'KEEP')

      


      call qexit('cc_r12mksfr')
      end
*=====================================================================*

      subroutine cc_r12xsort5(sajkl,blajkl,nvir0)
c---------------------------------------------------------------------

      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccfro.h>
      integer isymkl,isymij,isyml,isymk,isymaj
      integer isymj,isyma,idxkl,isym,nvir0(8)
      integer kt,idx,nvir0t,koffj,koffa,koffk,koffl
      integer icoun1,kl,aj,ivir0(8),ajkl,nrhfst
#if defined (SYS_CRAY)
      real blajkl(*),sajkl(*)
#else
      double precision blajkl(*),sajkl(*)
#endif
c
      call qenter('cc_r12xsort5')
c
      nvir0t = 0
      nrhfst = 0
      do isym = 1, nsym
         nvir0t      = nvir0t + nvir0(isym)
         nrhfst      = nrhfst + nrhfs(isym)
      end do

      icoun1 = 0
      do isym = 1,nsym
         ivir0(isym) = icoun1
         icoun1 = icoun1 + nvir0(isym)
      enddo


      idx = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhf(isymk)
              koffk = irhf(isymk)+k
             do l = 1, nrhf(isyml)
                koffl = irhf(isyml)+l
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    do j = 1, nrhf(isymj)
                       koffj = irhf(isymj)+j
                       do a = 1, nvir0(isyma)
                          idx = idx + 1
                          koffa = ivir0(isyma)+a
                          kl = (koffk-1)*nrhft + koffl
                          aj = (koffj-1)*nvir0t + koffa
                          ajkl = aj +(kl-1)*nvir0t*nrhft
                          blajkl(idx) = sajkl(ajkl)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
         enddo
      enddo

      call qexit('cc_r12xsort5')
      end

*---------------------------------------------------------------------------
      subroutine cc_r12xsort6(sajkl,blajkl,nvir0)
c---------------------------------------------------------------------

      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <ccfro.h>
      integer isymkl,isymij,isyml,isymk,isymaj
      integer isymj,isyma,idxkl,isym,nvir0(8)
      integer kt,idx,nvir0t,koffj,koffa,koffk,koffl
      integer icoun1,kl,aj,ivir0(8),ajkl
#if defined (SYS_CRAY)
      real blajkl(*),sajkl(*)
#else
      double precision blajkl(*),sajkl(*)
#endif
c
      call qenter('cc_r12xsort6')
c
      nvir0t = 0
      do isym = 1, nsym
         nvir0t      = nvir0t + nvir0(isym)
      end do

      icoun1 = 0
      do isym = 1,nsym
         ivir0(isym) = icoun1
         icoun1 = icoun1 + nvir0(isym)
      enddo


      idx = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhf(isymk)
              koffk = irhf(isymk)+k
             do l = 1, nrhf(isyml)
                koffl = irhf(isyml)+l
                idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    do j = 1, nrhf(isymj)
                       koffj = irhf(isymj)+j
                       do a = 1, nvir0(isyma)
                          idx = idx + 1
                          koffa = ivir0(isyma)+a
                          kl = (koffk-1)*nrhft + koffj
                          aj = (koffl-1)*nvir0t + koffa
                          ajkl = aj +(kl-1)*nvir0t*nrhft
                          blajkl(idx) = sajkl(ajkl)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
         enddo
      enddo

      call qexit('cc_r12xsort6')
      end


*=====================================================================*
