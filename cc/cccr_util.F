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
c/* Deck CC_GATHEROO */
*=====================================================================*
      SUBROUTINE CC_GATHEROO(XMO,XOO,ISYINT)
*---------------------------------------------------------------------*
*     Purpose: gather occupied/occupied integrals from one-electron
*              MO integral matrix and store them according to IMATIJ
*
*     Christof Haettig, January 1997
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYINT, KOFF1, KOFF2, ISYMJ, ISYMI
      
#if defined (SYS_CRAY)
      REAL XMO(*)
      REAL XOO(*)      ! dimension (NMATIJ(ISYINT))
#else
      DOUBLE PRECISION XMO(*)
      DOUBLE PRECISION XOO(*)      ! dimension (NMATIJ(ISYINT))
#endif

      DO ISYMJ = 1, NSYM
        ISYMI = MULD2H(ISYMJ,ISYINT)
        DO J = 1, NRHF(ISYMJ)
          KOFF1 = IFCRHF(ISYMI,ISYMJ) + NORB(ISYMI)*(J-1) + 1
          KOFF2 = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J-1) + 1
          CALL DCOPY(NRHF(ISYMI),XMO(KOFF1),1,XOO(KOFF2),1)
        END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CC_GATHEROV */
*=====================================================================*
      SUBROUTINE CC_GATHEROV(XMO,XOV,ISYINT)
*---------------------------------------------------------------------*
*     Purpose: gather occupied/virtual integrals from one-electron
*              MO integral matrix and store them according to IT1AM
*
*     Christof Haettig, January 1997
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYINT, KOFF1, KOFF2, ISYMJ, ISYMB
      
#if defined (SYS_CRAY)
      REAL XMO(*)
      REAL XOO(*)     ! dimension (NT1AM(ISYINT))
#else
      DOUBLE PRECISION XMO(*)
      DOUBLE PRECISION XOV(*)     ! dimension (NT1AM(ISYINT))
#endif

      DO ISYMJ = 1, NSYM
        ISYMB = MULD2H(ISYMJ,ISYINT)
        DO J = 1, NRHF(ISYMJ)
          KOFF1 = IFCVIR(ISYMJ,ISYMB) + J 
          KOFF2 = IT1AM(ISYMB,ISYMJ)  + NVIR(ISYMB)*(J-1) + 1
          CALL DCOPY(NVIR(ISYMB),XMO(KOFF1),NORB(ISYMJ),XOV(KOFF2),1)
        END DO
      END DO

      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CC_GATHERVV */
*=====================================================================*
      SUBROUTINE CC_GATHERVV(XMO,XVV,ISYINT)
*---------------------------------------------------------------------*
*     Purpose: gather occupied/occupied integrals from one-electron
*              MO integral matrix and store them according to IMATAB
*
*     Christof Haettig, January 1997
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYINT, KOFF1, KOFF2, ISYMB, ISYMC
      
#if defined (SYS_CRAY)
      REAL XMO(*)
      REAL XOO(*)  ! dimension (NMATAB(ISYINT))
#else
      DOUBLE PRECISION XMO(*)
      DOUBLE PRECISION XVV(*)  ! dimension (NMATAB(ISYINT))
#endif

      DO ISYMC = 1, NSYM
        ISYMB = MULD2H(ISYMC,ISYINT)
        DO C = 1, NVIR(ISYMC)
          KOFF1 = IFCVIR(ISYMB,ISYMC) + NORB(ISYMB)*(C-1)+NRHF(ISYMB)+1
          KOFF2 = IMATAB(ISYMB,ISYMC) + NVIR(ISYMB)*(C-1)+1
          CALL DCOPY(NVIR(ISYMB),XMO(KOFF1),1,XVV(KOFF2),1)
        END DO
      END DO


      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CC_XBAR */
*=====================================================================*
      SUBROUTINE CC_XBAR(XBAR,XLIAJB,ISYOVOV,T2AMP,ISYAMP,WORK,LWORK)
*---------------------------------------------------------------------*
*     Purpose: calculate XBAR intermediate needed for EMAT2
*              the XBAR intermediate in similar to the X intermediate
*              in the left transformation, but the two indeces 
*              transposed and the Zeta vector exchanged by L(ia|jb). 
*
*     symmetries:   ISYOVOV -- XLIAJB
*                   ISYAMP  -- T2AMP
*
*     Christof Haettig, January 1997
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYMX, ISYMI, ISYMJ, MAXJ, NIJ, NJI
      INTEGER ISYOVOV, ISYAMP, LWORK
      
#if defined (SYS_CRAY)
      REAL XBAR(*)      ! dimension (NMATIJ(MULD2H(ISYOVOV,ISYAMP)))
      REAL XLIAJB(*)    ! dimension (NT2SQ(ISYOVOV))
      REAL T2AMP(*)     ! dimension (NT2AM(ISYAMP))
      REAL WORK(LWORK)
      REAL SWAP
#else
      DOUBLE PRECISION XBAR(*)   ! dim. (NMATIJ(MULD2H(ISYOVOV,ISYAMP)))
      DOUBLE PRECISION XLIAJB(*) ! dim. (NT2SQ(ISYOVOV))
      DOUBLE PRECISION T2AMP(*)  ! dim. (NT2AM(ISYAMP))
      DOUBLE PRECISION WORK(LWORK)
      DOUBLE PRECISION SWAP
#endif


* call CC_XI for the expansive part:
      CALL CC_XI(XBAR,XLIAJB,ISYOVOV,T2AMP,ISYAMP,WORK,LWORK)

* transpose the result:
      ISYMX = MULD2H(ISYOVOV,ISYAMP)
      DO ISYMI = 1, NSYM
        ISYMJ = MULD2H(ISYMI,ISYMX)
        IF (ISYMJ .LE. ISYMI) THEN
          DO I = 1, NRHF(ISYMI)
            MAXJ =  NRHF(ISYMJ)
            IF (ISYMJ .EQ. ISYMI) MAXJ = I-1
          DO J = 1, MAXJ
            NIJ = IMATIJ(ISYMI,ISYMJ) + NRHF(ISYMI)*(J-1) + I
            NJI = IMATIJ(ISYMJ,ISYMI) + NRHF(ISYMJ)*(I-1) + J
            SWAP      = XBAR(NIJ)
            XBAR(NIJ) = XBAR(NJI)
            XBAR(NJI) = SWAP
          END DO
          END DO
        END IF
      END DO


      RETURN
      END
*---------------------------------------------------------------------*
c/* Deck CC_YBAR */
*=====================================================================*
      SUBROUTINE CC_YBAR(YBAR,XLIAJB,ISYOVOV,T2AMP,ISYAMP,WORK,LWORK)
*---------------------------------------------------------------------*
*     Purpose: calculate YBAR intermediate needed for EMAT1
*              the YBAR intermediate in similar to the X intermediate
*              in the left transformation, but has the opposite sign
*              and the Zeta is vector exchanged by L(ia|jb). 
*
*     symmetries:   ISYOVOV -- XLIAJB
*                   ISYAMP  -- T2AMP
*
*     Christof Haettig, January 1997
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
# include <implicit.h>
#endif
#include <priunit.h>
#include <ccsdsym.h>
#include <ccorb.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYMY, LEN
      INTEGER ISYOVOV, ISYAMP, LWORK
      
#if defined (SYS_CRAY)
      REAL YBAR(*)       ! dimension (NMATAB(MULD2H(ISYOVOV,ISYAMP)))
      REAL XLIAJB(*)     ! dimension (NT2SQ(ISYOVOV))
      REAL T2AMP(*)      ! dimension (NT2AM(ISYAMP))
      REAL WORK(LWORK)
#else
      DOUBLE PRECISION YBAR(*)   ! dim. (NMATAB(MULD2H(ISYOVOV,ISYAMP)))
      DOUBLE PRECISION XLIAJB(*) ! dim. (NT2SQ(ISYOVOV))
      DOUBLE PRECISION T2AMP(*)  ! dim. (NT2AM(ISYAMP))
      DOUBLE PRECISION WORK(LWORK)
#endif


* call CC_YI for the expansive part:
      CALL CC_YI(YBAR,XLIAJB,ISYOVOV,T2AMP,ISYAMP,WORK,LWORK)

* invert the sign of the result:
      ISYMY = MULD2H(ISYOVOV,ISYAMP)
      LEN   = NMATAB(ISYMY)
      CALL DSCAL(LEN,-1.0d0,YBAR,1)


      RETURN
      END
*---------------------------------------------------------------------*
c /* deck iccset1 */
*=====================================================================*
       INTEGER FUNCTION ICCSET1(IARRAY,LIST,IDLST,NARR,MAXARR,LAPPEND)
*---------------------------------------------------------------------*
*
*    Purpose: set up and maintain a nonredundant array of response
*             vectors characterized by LIST and IDLST
*
*             IARRAY - array of vectors
*             LIST   - vector type, given as a list name 'R0','R1',...
*             IDLST  - index of the vector in LIST
*        
*             NARR   - used length of IARRAY
*             MAXARR - maximum dimension of IARRAY
*
*    Written by Christof Haettig, january 1997.
*    PL1 added march 2000, Sonia
*
*=====================================================================*
      IMPLICIT NONE  

#include <priunit.h>
#include "cciccset.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      LOGICAL LAPPEND
      CHARACTER*(*) LIST
      INTEGER NARR,MAXARR, IDLST, ILIST, I, INDEX, LEN
      INTEGER IARRAY(2,MAXARR)

      LOGICAL LTST, LEND
      INTEGER IT,IV,IDX

* statement  functions:
      LTST(IDX,IT,IV) = IARRAY(1,IDX).EQ.IV .AND. IARRAY(2,IDX).EQ.IT
      LEND(IDX) = IDX.GT.NARR

*---------------------------------------------------------------------*
* if LOCDBG echo input:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'LIST, IDLST:',LIST, IDLST
        WRITE (LUPRI,*) 'NARR,MAXARR:',NARR,MAXARR
        WRITE (LUPRI,*) 'LAPPEND:',LAPPEND
        WRITE (LUPRI,*) 'IARRAY:'
        WRITE (LUPRI,'((/5X,2I5))')  
     &        ((IARRAY(I,INDEX),I=1,2),INDEX=1,NARR)
      END IF

*---------------------------------------------------------------------*
* maintain list of nonredundant vectors:
*---------------------------------------------------------------------*

* convert LIST to an integer:
      IF (      LIST(1:1).EQ.'R' .OR. LIST(1:1).EQ.'O'
     &     .OR. LIST(1:1).EQ.'L' .OR. LIST(1:1).EQ.'X'
     &     .OR. LIST(1:1).EQ.'N' .OR. LIST(1:1).EQ.'M' 
     &     .OR. LIST(1:2).EQ.'D0'                      ) THEN
         LEN = 2
      ELSE IF  (LIST(1:1).EQ.'E' .OR. LIST(1:1).EQ.'C' 
     &                           .OR. LIST(1:1).EQ.'P') THEN
         LEN = 3
      END IF
      ILIST = 0
      DO I = 1, MAXTAB
       IF (VTABLE(I)(1:LEN).EQ.LIST(1:LEN)) ILIST = I
      END DO
      IF (ILIST .EQ.0) CALL QUIT('Unknown list in ICCSET1.')

      INDEX = 1
      DO WHILE ( .NOT. (LTST(INDEX,ILIST,IDLST) .OR. LEND(INDEX)) )
        INDEX = INDEX + 1
      END DO

      IF (INDEX.GT.MAXARR) THEN
        CALL QUIT('ERROR> list overflow in CCCR_SET1.')
      END IF

      IF (.NOT.LTST(INDEX,ILIST,IDLST) .OR. LEND(INDEX)) THEN ! append:
        IF (LAPPEND) THEN
          IARRAY(1,INDEX) = IDLST
          IARRAY(2,INDEX) = ILIST
          NARR = NARR + 1
        ELSE
          CALL QUIT('ICCSET1 failed: requested element was '//
     &          'not registered.')
        END IF
      END IF

      ICCSET1 = INDEX

      RETURN 
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE ICCSET1                              *
*---------------------------------------------------------------------*
c /* deck iccset2 */
*=====================================================================*
      INTEGER FUNCTION ICCSET2(IARRAY,LISTA,IDLSTA,LISTB,IDLSTB,
     &                         NARR,MAXARR,LAPPEND)
*---------------------------------------------------------------------*
*
*    Purpose: set up and maintain a nonredundant array of pairs of
*             response vectors characterized by LIST and IDLST
*
*            IARRAY   - array of vectors
*            LISTA/B  - vector type, given as a list name 'R0','R1',...
*            IDLSTA/B - index of the vector in LIST
*       
*            MAXARR - maximum dimension of IARRAY
*
*    Written by Christof Haettig, january 1997.
*    PL1 vectors, Sonia 2000
*
*=====================================================================*
      IMPLICIT NONE  
#include <priunit.h>

#include "cciccset.h"

      LOGICAL LAPPEND
      CHARACTER*(*) LISTA, LISTB
      INTEGER NARR,MAXARR, IDLSTA, IDLSTB, I
      INTEGER IARRAY(4,MAXARR)

      LOGICAL LTST, LEND
      INTEGER ITA,ITB,IVA,IVB,IDX,ILISTA,ILISTB, INDEX, LENA, LENB

* statement  functions:
      LTST(IDX,ITA,IVA,ITB,IVB) = 
     &  ( IARRAY(1,IDX).EQ.IVA .AND. IARRAY(2,IDX).EQ.ITA .AND.
     &    IARRAY(3,IDX).EQ.IVB .AND. IARRAY(4,IDX).EQ.ITB      ) .OR.
     &  ( IARRAY(3,IDX).EQ.IVA .AND. IARRAY(4,IDX).EQ.ITA .AND.
     &    IARRAY(1,IDX).EQ.IVB .AND. IARRAY(2,IDX).EQ.ITB      ) 

      LEND(IDX) = IDX.GT.NARR

*---------------------------------------------------------------------*
* maintain list of nonredundant vectors:
*---------------------------------------------------------------------*

* convert LISTA & LISTB to an integer:
      IF (      LISTA(1:1).EQ.'R' .OR. LISTA(1:1).EQ.'O'
     &     .OR. LISTA(1:1).EQ.'L' .OR. LISTA(1:1).EQ.'X'
     &     .OR. LISTA(1:1).EQ.'N' .OR. LISTA(1:1).EQ.'M' 
     &     .OR. LISTA(1:2).EQ.'D0'                       ) THEN
         LENA = 2
      ELSE IF  (LISTA(1:1).EQ.'E' .OR. LISTA(1:1).EQ.'C' 
     &     .OR. LISTA(1:1).EQ.'P'                        ) THEN
         LENA = 3
      END IF
      IF (      LISTB(1:1).EQ.'R' .OR. LISTB(1:1).EQ.'O'
     &     .OR. LISTB(1:1).EQ.'L' .OR. LISTB(1:1).EQ.'X'
     &     .OR. LISTB(1:1).EQ.'N' .OR. LISTB(1:1).EQ.'M' 
     &     .OR. LISTB(1:2).EQ.'D0'                       ) THEN
         LENB = 2
      ELSE IF  (LISTB(1:1).EQ.'E' .OR. LISTB(1:1).EQ.'C' 
     &     .OR. LISTB(1:1).EQ.'P'                        ) THEN
         LENB = 3
      END IF
      ILISTA = 0
      ILISTB = 0
      DO I = 1, MAXTAB
       IF (VTABLE(I)(1:LENA).EQ.LISTA(1:LENA)) ILISTA = I
       IF (VTABLE(I)(1:LENB).EQ.LISTB(1:LENB)) ILISTB = I
      END DO
      IF (ILISTA .EQ.0) CALL QUIT('Unknown list in ICCSET2.')
      IF (ILISTB .EQ.0) CALL QUIT('Unknown list in ICCSET2.')


      INDEX = 1
      DO WHILE ( .NOT. (LTST(INDEX,ILISTA,IDLSTA,ILISTB,IDLSTB) 
     &                  .OR. LEND(INDEX)) )
        INDEX = INDEX + 1
      END DO

      IF (INDEX.GT.MAXARR) THEN
        CALL QUIT('ERROR> list overflow in CCCR_SET2.')
      END IF

      IF (.NOT.LTST(INDEX,ILISTA,IDLSTA,ILISTB,IDLSTB)
     &    .OR.LEND(INDEX)) THEN 
        IF (LAPPEND) THEN
          IARRAY(1,INDEX) = IDLSTA
          IARRAY(2,INDEX) = ILISTA
          IARRAY(3,INDEX) = IDLSTB
          IARRAY(4,INDEX) = ILISTB
          NARR = NARR + 1
        ELSE
          WRITE (LUPRI,*) 'error in ICCSET2:'
          WRITE (LUPRI,*) 'requested pair was not registered.'
          WRITE (LUPRI,*) IDLSTA, ILISTA, IDLSTB, ILISTB
          WRITE (LUPRI,*) 'LIST:'
          WRITE(*,'((/5X,4I5))') ((IARRAY(I,IDX),I=1,4),IDX=1,NARR)
          CALL QUIT(
     *         'ICCSET2 failed: requested pair was not registered.')
        END IF
      END IF

      ICCSET2 = INDEX

      RETURN 
      END 

*---------------------------------------------------------------------*
*              END OF SUBROUTINE ICCSET2                              *
*---------------------------------------------------------------------*
