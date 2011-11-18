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
C  /* Deck cc_den2_pt */
      SUBROUTINE CC_DEN2_PT(D2IJG,D2AIG,D2IAG,D2ABG,
     &                      XCMO,ISYCMO,WORK,LWORK,
     &                      IDEL,ISYMD)
C
C     Written by Sonia Coriani, 03/02-2002, based on CC_DEN2
C
C     Version: 1.0
C
C     Purpose: Directs the calculation of the 2 electron CC-(T) density
C              d(pq,gam;del) for a given del (IDEL). The 4 blocks pq
C              of the result are returned in D2IJG through D2ABG!
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      INTEGER ISYCMO
      DIMENSION D2IJG(*), D2AIG(*), D2IAG(*), D2ABG(*)
      DIMENSION XCMO(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
C Locally define density file names for easier debug
C
      CHARACTER FNDSIJKD*9, FNDSIAJD*9, FNDSAIJD*9, FNDSABID*9
      CHARACTER FNDSAIBD*9, FNDSIABD*9
      PARAMETER (FNDSIJKD = 'PT_DSIJKD', FNDSIAJD = 'PT_DSIAJD')
      PARAMETER (FNDSAIJD = 'PT_DSAIJD', FNDSABID = 'PT_DSABID')
      PARAMETER (FNDSAIBD = 'PT_DSAIBD', FNDSIABD = 'PT_DSIABD')
C
      LOGICAL LOCDBG
      PARAMETER(LOCDBG = .false.)
C
      CALL QENTER('CC_DEN2_PT')
C
C-------------------------------
C     set some symmetries:
C-------------------------------
C
      ISYD2H = MULD2H(ISYMD,ISYCMO)
      ISYDEN = MULD2H(ISYD2H,ISYCMO) 
      D = IDEL - IBAS(ISYMD)
C
C---------------------------------------------------------
C     Open files for backtransformed/symmetrized densities
C---------------------------------------------------------
C
      LUIJKDS  = -1
      LUABIDS  = -1
C
      LUIAJDE  = -1
      LUAIJDE  = -1
      LUAIBDE  = -1
      LUIABDE  = -1
C
C     d^s_{ijkdelta}
      CALL WOPEN2(LUIJKDS,FNDSIJKD,64,0)
C     d^s_{abidelta}
      CALL WOPEN2(LUABIDS,FNDSABID,64,0)
C     d^s_{iajdelta}
      CALL WOPEN2(LUIAJDE,FNDSIAJD,64,0)
C     d^s_{aijdelta}
      CALL WOPEN2(LUAIJDE,FNDSAIJD,64,0)
C     d^s_{aibdelta}
      CALL WOPEN2(LUAIBDE,FNDSAIBD,64,0)
C     d^s_{iabdelta}
      CALL WOPEN2(LUIABDE,FNDSIABD,64,0)
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KD2IJK = 1
      KEND1  = KD2IJK + NMAIJK(ISYD2H)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insuff. space for 1st alloc. in CC_DEN2_PT')
      ENDIF
C
C-------------------------------------------------------------
C     Read-in the symmetrized (occ.occ,occ;del) density block.
C-------------------------------------------------------------
C
      IOFF = I3ODEL(ISYD2H,ISYMD) +  NMAIJK(ISYD2H)*(D-1) + 1
      CALL GETWA2(LUIJKDS,FNDSIJKD,WORK(KD2IJK),
     &                             IOFF,NMAIJK(ISYD2H))
      if (locdbg) then
         xtest = ddot(NMAIJK(ISYD2H),WORK(KD2IJK),1,WORK(KD2IJK),1)
         write(lupri,*)'DEN2_PT: norm (i.j,k;del)', xtest, 
     &                                           isymd, idel, d
      end if
C
C------------------------------------------------------------------
C     Backtransform third occupied index to AO and store in result.
C------------------------------------------------------------------
C
      if (.true.) then
      ICON = 1
      CALL CCD2_PQAO(D2IJG,ISYDEN,WORK(KD2IJK),ISYD2H,
     *               XCMO,ISYCMO,ICON)
      end if
C
C-------------------------------
C     Work space allocation two.
C-------------------------------
C
      KD2ABI = 1
      KEND2  = KD2ABI + NCKASR(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insuff. space for 2nd alloc. in CC_DEN2_PT')
      ENDIF
C
C------------------------------------------------------
C     Read-in the symmetrized (vir.vir,occ;del) density
C------------------------------------------------------
C

      IOFF = ICDKAO(ISYD2H,ISYMD) + NCKASR(ISYD2H)*(D-1)+1
      CALL GETWA2(LUABIDS,FNDSABID,WORK(KD2ABI),
     &                                IOFF,NCKASR(ISYD2H))
      if (locdbg) then
      xtest = ddot(NCKASR(ISYD2H),WORK(KD2ABI),1,WORK(KD2ABI),1)
      write(lupri,*) 'DEN2_PT: norm (a.b,i;del)', xtest, 
     &                                           isymd, idel, d
      end if
C
C---------------------------------------
C     Backtransform third occupied index 
C---------------------------------------
C
      ICON = 3
      CALL CCD2_PQAO(D2ABG,ISYDEN,WORK(KD2ABI),ISYD2H,
     *               XCMO,ISYCMO,ICON)
C
C---------------------------------------------------------------
C     Read-in (occ.vir,occ;del) density [stored vir.occ,occ;del]
C     Read-in (occ.vir,vir;del) density [stored vir.occ,vir;del]
C---------------------------------------------------------------
C
      KD2IAJ = 1
      KD2IAB = KD2IAJ + NT2BCD(ISYD2H)
      KEND2  = KD2IAB + NCKATR(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insuff. space for 3rd alloc. in CC_DEN2_PT')
      ENDIF

      IOFF = ICKID(ISYD2H,ISYMD) + NT2BCD(ISYD2H)*(D-1) + 1 
      CALL GETWA2(LUIAJDE,FNDSIAJD,WORK(KD2IAJ),IOFF,NT2BCD(ISYD2H))

      if (locdbg) then
      xtest = ddot(NT2BCD(ISYD2H),WORK(KD2IAJ),1,WORK(KD2IAJ),1)
      write(lupri,*) 'DEN2_PT: norm (i.a,j;del)', xtest, 
     &                                           isymd, idel, d
      end if

      IOFF = ICKDAO(ISYD2H,ISYMD) + NCKATR(ISYD2H)*(D-1) + 1 
      CALL GETWA2(LUIABDE,FNDSIABD,WORK(KD2IAB),IOFF,NCKATR(ISYD2H))

      if (locdbg) then
      xtest = ddot(NCKATR(ISYD2H),WORK(KD2IAB),1,WORK(KD2IAB),1)
      write(lupri,*) 'DEN2_PT: norm (i.a,b;del)', xtest, 
     &                                           isymd, idel, d
      end if
C
C---------------------------------------------------------------------
C     Backtransform occupied/virtual index to AO and store in results.
C---------------------------------------------------------------------
C
      ICON = 2
      CALL CCD2_PQAO(D2IAG,ISYDEN,WORK(KD2IAJ),ISYD2H,
     *               XCMO,ISYCMO,ICON)

      ICON = 5
      CALL CCD2_PQAO(D2IAG,ISYDEN,WORK(KD2IAB),ISYD2H,
     *               XCMO,ISYCMO,ICON)
C
C--------------------------------------------------
C     Read-in (vir.occ,vir;del) & (vir.occ,occ;del).
C--------------------------------------------------
C
      KD2AIB = 1
      KD2AIJ = KD2AIB + NCKATR(ISYD2H)
      KEND2  = KD2AIJ + NT2BCD(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insufficient space for fifth allocation in CC_DEN2')
      ENDIF
C
      IOFF = ICKDAO(ISYD2H,ISYMD) + NCKATR(ISYD2H)*(D-1) + 1
      CALL GETWA2(LUAIBDE,FNDSAIBD,WORK(KD2AIB),IOFF,NCKATR(ISYD2H))

      if (locdbg) then
      xtest = ddot(NCKATR(ISYD2H),WORK(KD2AIB),1,WORK(KD2AIB),1)
      write(lupri,*) 'DEN2_PT: norm (a.i,b;del)', xtest, 
     &                                            isymd, idel, d
      end if

      IOFF = ICKID(ISYD2H,ISYMD) + NT2BCD(ISYD2H)*(D-1) + 1
      CALL GETWA2(LUAIJDE,FNDSAIJD,WORK(KD2AIJ),IOFF,NT2BCD(ISYD2H))

      if (locdbg) then
      xtest = ddot(NT2BCD(ISYD2H),WORK(KD2AIJ),1,WORK(KD2AIJ),1)
      write(lupri,*) 'DEN2_PT: norm (a.i,j;del)', xtest, 
     &                                            isymd, idel, d
      end if

C-------------------------------------------------------
C     Backtransform third index of the (vir.occ,vir;del)
C     and (vir.occ,occ;del) blocks to AO-basis and 
C     store in result.
C-------------------------------------------------------
C
      ICON = 5
      CALL CCD2_PQAO(D2AIG,ISYDEN,WORK(KD2AIB),ISYD2H,
     *               XCMO,ISYCMO,ICON)

      ICON = 2
      CALL CCD2_PQAO(D2AIG,ISYDEN,WORK(KD2AIJ),ISYD2H,
     *               XCMO,ISYCMO,ICON)
C
C---------------------------------------------------
C     Close files
C---------------------------------------------------
C
      CALL WCLOSE2(LUIJKDS,FNDSIJKD,'KEEP')
      CALL WCLOSE2(LUABIDS,FNDSABID,'KEEP')
      CALL WCLOSE2(LUIAJDE,FNDSIAJD,'KEEP')
      CALL WCLOSE2(LUAIJDE,FNDSAIJD,'KEEP')
      CALL WCLOSE2(LUAIBDE,FNDSAIBD,'KEEP')
      CALL WCLOSE2(LUIABDE,FNDSIABD,'KEEP')
C
      CALL QEXIT('CC_DEN2_PT')
      RETURN
      END
C------------------------------------------------------
C  /* Deck cc_den2_ptanti */
      SUBROUTINE CC_DEN2_PTanti(D2IJG,D2AIG,D2IAG,D2ABG,
     &                      XCMO,ISYCMO,WORK,LWORK,
     &                      IDEL,ISYMD)
C
C     Written by Sonia Coriani, 03/02-2002, based on CC_DEN2
C
C     Version: 1.0
C
C     Purpose: Directs the calculation of the 2 electron CC-(T) density
C              d(pq,gam;del) for a given del (IDEL). The 4 blocks pq
C              of the result are returned in D2IJG through D2ABG!
C
C     Special antisymmetrizing version for orbital-orbital Breit-Pauli
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      INTEGER ISYCMO
      DIMENSION D2IJG(*), D2AIG(*), D2IAG(*), D2ABG(*)
      DIMENSION XCMO(*), WORK(LWORK)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <ccsdinp.h>
C
C Locally define density file names
C
      CHARACTER FNDAIJKD*9,FNDAABID*9
      PARAMETER (FNDAIJKD = 'PT_DAIJKD', FNDAABID = 'PT_DAABID')

      CHARACTER FNDSAIBD*9, FNDSIABD*9, FNDSIAJD*9, FNDSAIJD*9
      PARAMETER (FNDSAIJD = 'PT_DSAIJD', FNDSIAJD = 'PT_DSIAJD')
      PARAMETER (FNDSAIBD = 'PT_DSAIBD', FNDSIABD = 'PT_DSIABD')
C
      LOGICAL LOCDBG
      PARAMETER(LOCDBG = .false.)
C
      CALL QENTER('CC_DEN2_PTanti')
C
C-------------------------------
C     set some symmetries:
C-------------------------------
C
      ISYD2H = MULD2H(ISYMD,ISYCMO)
      ISYDEN = MULD2H(ISYD2H,ISYCMO) 
      D = IDEL - IBAS(ISYMD)
C
C---------------------------------------------------------
C     Open files for backtransformed/symmetrized densities
C---------------------------------------------------------
C
      LUIJKDA  = -1
      LUABIDA  = -1
C
      LUIAJDE  = -1
      LUAIJDE  = -1
      LUAIBDE  = -1
      LUIABDE  = -1
C
C     d^s_{ijkdelta}
      CALL WOPEN2(LUIJKDA,FNDAIJKD,64,0)
C     d^s_{abidelta}
      CALL WOPEN2(LUABIDA,FNDAABID,64,0)

C     d^s_{iajdelta}
      CALL WOPEN2(LUIAJDE,FNDSIAJD,64,0)
C     d^s_{aijdelta}
      CALL WOPEN2(LUAIJDE,FNDSAIJD,64,0)
C     d^s_{aibdelta}
      CALL WOPEN2(LUAIBDE,FNDSAIBD,64,0)
C     d^s_{iabdelta}
      CALL WOPEN2(LUIABDE,FNDSIABD,64,0)
C
C-------------------------------
C     Work space allocation one.
C-------------------------------
C
      KD2IJK = 1
      KEND1  = KD2IJK + NMAIJK(ISYD2H)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND1
         CALL QUIT('Insuff. space for 1st alloc. in CC_DEN2_PT')
      ENDIF
C
C-------------------------------------------------------------
C     Read-in the symmetrized (occ.occ,occ;del) density block.
C-------------------------------------------------------------
C
      IOFF = I3ODEL(ISYD2H,ISYMD) +  NMAIJK(ISYD2H)*(D-1) + 1
      CALL GETWA2(LUIJKDA,FNDAIJKD,WORK(KD2IJK),
     &                             IOFF,NMAIJK(ISYD2H))
      if (locdbg) then
         !call dzero(WORK(KD2IJK),NMAIJK(ISYD2H))
         xtest = ddot(NMAIJK(ISYD2H),WORK(KD2IJK),1,WORK(KD2IJK),1)
         write(lupri,*)'DEN2_PT: norm (i.j,k;del)', xtest, 
     &                                           isymd, idel, d
      end if
C
C------------------------------------------------------------------
C     Backtransform third occupied index to AO and store in result.
C------------------------------------------------------------------
C
      if (.true.) then
      ICON = 1
      CALL CCD2_PQAO(D2IJG,ISYDEN,WORK(KD2IJK),ISYD2H,
     *               XCMO,ISYCMO,ICON)
      end if
C
C-------------------------------
C     Work space allocation two.
C-------------------------------
C
      KD2ABI = 1
      KEND2  = KD2ABI + NCKASR(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insuff. space for 2nd alloc. in CC_DEN2_PT')
      ENDIF
C
C------------------------------------------------------
C     Read-in the symmetrized (vir.vir,occ;del) density
C------------------------------------------------------
C

      IOFF = ICDKAO(ISYD2H,ISYMD) + NCKASR(ISYD2H)*(D-1)+1
      CALL GETWA2(LUABIDA,FNDAABID,WORK(KD2ABI),
     &                                IOFF,NCKASR(ISYD2H))
      if (locdbg) then
      !call dzero(WORK(KD2ABI),NCKASR(ISYD2H))
      xtest = ddot(NCKASR(ISYD2H),WORK(KD2ABI),1,WORK(KD2ABI),1)
      write(lupri,*) 'DEN2_PT: norm (a.b,i;del)', xtest, 
     &                                           isymd, idel, d
      end if
C
C---------------------------------------
C     Backtransform third occupied index 
C---------------------------------------
C
      ICON = 3
      CALL CCD2_PQAO(D2ABG,ISYDEN,WORK(KD2ABI),ISYD2H,
     *               XCMO,ISYCMO,ICON)
C
C---------------------------------------------------------------
C     Read-in (occ.vir,occ;del) density [stored vir.occ,occ;del]
C     Read-in (occ.vir,vir;del) density [stored vir.occ,vir;del]
C---------------------------------------------------------------
C
      KD2IAJ = 1
      KD2IAB = KD2IAJ + NT2BCD(ISYD2H)
      KEND2  = KD2IAB + NCKATR(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insuff. space for 3rd alloc. in CC_DEN2_PT')
      ENDIF

      IOFF = ICKID(ISYD2H,ISYMD) + NT2BCD(ISYD2H)*(D-1) + 1 
      CALL GETWA2(LUIAJDE,FNDSIAJD,WORK(KD2IAJ),IOFF,NT2BCD(ISYD2H))

      if (locdbg) then
      xtest = ddot(NT2BCD(ISYD2H),WORK(KD2IAJ),1,WORK(KD2IAJ),1)
      write(lupri,*) 'DEN2_PT: norm (i.a,j;del)', xtest, 
     &                                           isymd, idel, d
      end if

      IOFF = ICKDAO(ISYD2H,ISYMD) + NCKATR(ISYD2H)*(D-1) + 1 
      CALL GETWA2(LUIABDE,FNDSIABD,WORK(KD2IAB),IOFF,NCKATR(ISYD2H))

      if (locdbg) then
      xtest = ddot(NCKATR(ISYD2H),WORK(KD2IAB),1,WORK(KD2IAB),1)
      write(lupri,*) 'DEN2_PT: norm (i.a,b;del)', xtest, 
     &                                           isymd, idel, d
      end if
C
C---------------------------------------------------------------------
C     Backtransform occupied/virtual index to AO and store in results.
C---------------------------------------------------------------------
C
      ICON = 2
      CALL CCD2_PQAO(D2IAG,ISYDEN,WORK(KD2IAJ),ISYD2H,
     *               XCMO,ISYCMO,ICON)

      ICON = 5
      CALL CCD2_PQAO(D2IAG,ISYDEN,WORK(KD2IAB),ISYD2H,
     *               XCMO,ISYCMO,ICON)
C
C--------------------------------------------------
C     Read-in (vir.occ,vir;del) & (vir.occ,occ;del).
C--------------------------------------------------
C
      KD2AIB = 1
      KD2AIJ = KD2AIB + NCKATR(ISYD2H)
      KEND2  = KD2AIJ + NT2BCD(ISYD2H)
      LWRK2  = LWORK  - KEND2
C
      IF (LWRK2 .LT. 0) THEN
         WRITE(LUPRI,*) 'Available:', LWORK, 'Needed:', KEND2
         CALL QUIT('Insufficient space for fifth allocation in CC_DEN2')
      ENDIF
C
      IOFF = ICKDAO(ISYD2H,ISYMD) + NCKATR(ISYD2H)*(D-1) + 1
      CALL GETWA2(LUAIBDE,FNDSAIBD,WORK(KD2AIB),IOFF,NCKATR(ISYD2H))

      if (locdbg) then
      xtest = ddot(NCKATR(ISYD2H),WORK(KD2AIB),1,WORK(KD2AIB),1)
      write(lupri,*) 'DEN2_PT: norm (a.i,b;del)', xtest, 
     &                                            isymd, idel, d
      end if

      IOFF = ICKID(ISYD2H,ISYMD) + NT2BCD(ISYD2H)*(D-1) + 1
      CALL GETWA2(LUAIJDE,FNDSAIJD,WORK(KD2AIJ),IOFF,NT2BCD(ISYD2H))

      if (locdbg) then
      xtest = ddot(NT2BCD(ISYD2H),WORK(KD2AIJ),1,WORK(KD2AIJ),1)
      write(lupri,*) 'DEN2_PT: norm (a.i,j;del)', xtest, 
     &                                            isymd, idel, d
      end if

C-------------------------------------------------------
C     Backtransform third index of the (vir.occ,vir;del)
C     and (vir.occ,occ;del) blocks to AO-basis and 
C     store in result.
C-------------------------------------------------------
C
      ICON = 5
      CALL CCD2_PQAO(D2AIG,ISYDEN,WORK(KD2AIB),ISYD2H,
     *               XCMO,ISYCMO,ICON)

      ICON = 2
      CALL CCD2_PQAO(D2AIG,ISYDEN,WORK(KD2AIJ),ISYD2H,
     *               XCMO,ISYCMO,ICON)
C
C---------------------------------------------------
C     Close files
C---------------------------------------------------
C
      CALL WCLOSE2(LUIJKDA,FNDAIJKD,'KEEP')
      CALL WCLOSE2(LUABIDA,FNDAABID,'KEEP')

      CALL WCLOSE2(LUIAJDE,FNDSIAJD,'KEEP')
      CALL WCLOSE2(LUAIJDE,FNDSAIJD,'KEEP')
      CALL WCLOSE2(LUAIBDE,FNDSAIBD,'KEEP')
      CALL WCLOSE2(LUIABDE,FNDSIABD,'KEEP')
C
      CALL QEXIT('CC_DEN2_PTanti')
      RETURN
      END
