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
C  /* Deck rspfxd */
      SUBROUTINE RSPFXD (FI,FA,QA,QB,H2AX,
     *                  IDAX,INTSYM, ISYMDN,DEN1,DEN2,
     *                  PVDEN,H2,H2X,WRK,KFREE,LFREE,
     *                  LCON,LORB,IPRFX,IGRSPI,ISPIN1,ISPIN2,ISPIN3,
     *                  IKLVL,ZYMAT1,ISYM1,ZYMAT2,ISYM2,ZYMAT3,ISYM3)
C
C
C 15-Nov-1991 hjaaj (pulled out of RSP2GR)
C
C     This routine computes the FX matrices, that is, the Fock
C     type matrices FI, FA, QA, and QB with transformed ("X")
C     integrals.  In addition, If LCON true then the H2AX matrix
C     with transformed integrals is extracted for the CI routines.
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, DH = 0.50D0 )
C
C Used from common blocks:
C  ?
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infinp.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infrsp.h>
#include <infpri.h>
#include <orbtypdef.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infspi.h>
#include <cbgetdis.h>
#include <infdis.h>
C
 
      DIMENSION FI(NORBT,NORBT),FA(NORBT,NORBT)
      DIMENSION QA(NORBT,NASHDI),QB(NORBT,NASHDI)
      DIMENSION H2(NORBT,NORBT), H2X(NORBT,NORBT)
      DIMENSION H2AX(N2ASHX*N2ASHX,*)
      DIMENSION DEN1(NASHDI,NASHDI), DEN2(*)
      DIMENSION PVDEN(NASHDI,NASHDI)
      DIMENSION ZYMAT1(NORBT,NORBT),ZYMAT2(NORBT,NORBT)
      DIMENSION ZYMAT3(NORBT,NORBT)
      DIMENSION WRK(*)
      DIMENSION NEEDED(6)
      DATA NEEDED /1,1,1,1,1,0/
C
      LOGICAL LCON, LORB
C
C     Put up the structure for reading the (transformed) integrals:
C     IEND < 0 flags the completeness of the distributions
C
      CALL QENTER('RSPFXD')
      IF (DOHFSRDFT .OR. DOMCSRDFT) THEN
         CALL QUIT('srDFT not implemented yet in RSPFXD')
      END IF
      KWRK = KFREE
      LWRK = LFREE
      ID1 = 1
      ID2 = N2ASHX*N2ASHX + 1
      CALL DZERO(H2,N2ORBX)
      IF (IKLVL.GE.2) NEEDED(6) = 1
      IF (IPRFX.GT.200) THEN
         WRITE(LUPRI,'(/A)')   '>>>>>> ENTERING RSPFXD <<<<<<'
         WRITE(LUPRI,'(A,I5)') 'IDAX =',IDAX
         WRITE(LUPRI,'(A,I5)') 'INTSYM =',INTSYM
         WRITE(LUPRI,'(A,I5)') 'ISYMDN =',ISYMDN
         WRITE(LUPRI,'(A,L1)') 'LCON =',LCON
         WRITE(LUPRI,'(A,L2)') 'LORB =',LORB
         WRITE(LUPRI,'(A,I5)') 'IPRFX =',IPRFX
         WRITE(LUPRI,'(A,I5)') 'IGRSPI =',IGRSPI
         WRITE(LUPRI,'(A,I5)') 'ISPIN1 =',ISPIN1
         WRITE(LUPRI,'(A,I5)') 'ISPIN2 =',ISPIN2
         WRITE(LUPRI,'(A,I5)') 'ISPIN3 =',ISPIN3
         WRITE(LUPRI,'(A,I5)') 'IKLVL =',IKLVL
         WRITE(LUPRI,'(A,I5)') 'ISYM1 =',ISYM1
         WRITE(LUPRI,'(A,I5)') 'ISYM2 =',ISYM2
         WRITE(LUPRI,'(A,I5)') 'ISYM3 =',ISYM3
         WRITE(LUPRI,'(/A)') 'One-electron FI'
         CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(/A)') 'One-electron density'
         CALL OUTPUT(DEN1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A)') 'Two-electron density'
         CALL PRIAC2(DEN2,NASHT,LUPRI)
         IF (TRPFLG) THEN
            WRITE(LUPRI,'(/A)') 'Two-electron density, triplet'
            CALL PRIAC2(DEN2(ID2),NASHT,LUPRI)
         END IF
      END IF
      IDIST = 0
      IPERCD = 1
C
      CALL MEMGET('REAL',KFIRST,1,WRK,KWRK,LWRK)
      IF (IKLVL.EQ.0) THEN
         CALL MEMGET('REAL',KFI,N2ORBX,WRK,KWRK,LWRK)
         CALL DZERO(WRK(KFI),N2ORBX)
      ELSE
         CALL MEMGET('REAL',KFI,0     ,WRK,KWRK,LWRK)
      END IF
      IF (IKLVL.EQ.3) THEN
         CALL MEMGET('REAL',KH2,N2ORBX,WRK,KWRK,LWRK)
      ELSE
         CALL MEMGET('REAL',KH2,0     ,WRK,KWRK,LWRK)
      END IF
      LSCR = 0
      IF (NASHT.GT.0) THEN
         CALL MEMGET('REAL',KH2XAC,N2ASHX,WRK,KWRK,LWRK)
         CALL DZERO(WRK(KH2XAC),N2ASHX)
         IF (LORB) LSCR = MAX(NORBT,N2ASHX)
      ELSE
         CALL MEMGET('REAL',KH2XAC,0,WRK,KWRK,LWRK)
      END IF
      CALL MEMGET('REAL',KSCR,LSCR,WRK,KWRK,LWRK)
C
1000  CALL NXTH2M(IC,ID,H2,NEEDED,WRK,KWRK,LWRK,IDIST)
      IF ( IDIST .LT. 0) GO TO 2000
C
C
C     Determine type of distribution
C
      ICTYP = IOBTYP(IC)
      IDTYP = IOBTYP(ID)
      ICDTYP = IDBTYP(ICTYP,IDTYP)
C
C     Determine symmetry of distribution
C
      ICSYM = ISMO( IC )
      IDSYM = ISMO( ID )
      ICDSYM = MULD2H(ICSYM,IDSYM)
C
C
C
      IF ( IPRFX .GT. 300) THEN
         WRITE(LUPRI,'(//A)') 'Characteristics of distribution:'
         WRITE(LUPRI,'(A)')   '================================'
         WRITE(LUPRI,'(A,2I5)')'IC ID   = ', IC,ID
         WRITE(LUPRI,'(3A)')   'TYP     = ', COBTYP(ICTYP),COBTYP(IDTYP)
         WRITE(LUPRI,'(A,2I5)')'Symmetry= ', ICSYM, IDSYM
         WRITE(LUPRI,'(A)')   '================================'
C
         CALL HEADER('Two-electron integrals from this distribution',3)
         CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
      IF (IDAX .EQ. 0) THEN
C
C Regular MO-integrals
C
         IF (IKLVL.EQ.0) THEN
            IPDA = IPDENS(IGRSPI,0)
            IPDB = IPDENS(0,IGRSPI)
            CALL EXTRAC(H2,WRK(KH2XAC))
            CALL FQDIS1(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                  WRK(KFI),FA,QA,QB,H2,WRK(KH2XAC),H2AX,
     *                  DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                  LORB,LCON,
     *                  IGRSPI,ISPIN1,ISPIN2,IPERCD)
         ELSE IF (IKLVL.EQ.1) THEN
            IH2SYM = MULD2H(ICDSYM,INTSYM)
C
C Transform left (~| )
C
            CALL DZERO(H2X,N2ORBX)
            CALL OITH1(ISYM1,ZYMAT1,H2,H2X,IH2SYM)
            INTSY1 = MULD2H(INTSYM,ISYM1)
            IF (IPRFX.GT.300) THEN
                WRITE(LUPRI,'(/A,I5)')'Left transform H2 of symmetry',
     *              IH2SYM
                CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A,I5)') 'with ZYMAT1 of symmetry',ISYM1
                CALL OUTPUT(ZYMAT1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A)') 'to H2X '
                CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
               IPDA = IPDENS(MULSP(IGRSPI,ISPIN1),ISPIN2)
               IPDB = IPDENS(ISPIN1,MULSP(IGRSPI,ISPIN2))
               CALL EXTRAC(H2X,WRK(KH2XAC))
               CALL FQDIS1(INTSY1,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                     FI,FA,QA,QB,H2X,WRK(KH2XAC),H2AX,
     *                     DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                     LORB,LCON,
     *                     IGRSPI,ISPIN1,ISPIN2,IPERCD)
         ELSE IF (IKLVL.EQ.2) THEN
            IH2SYM = MULD2H(ICDSYM,INTSYM)
C
C Transform left (~| )
C
            CALL DZERO(H2X,N2ORBX)
            CALL OITH1(ISYM1,ZYMAT1,H2,H2X,IH2SYM)
            INTSY1 = MULD2H(INTSYM,ISYM1)
            IF (IPRFX.GT.300) THEN
                WRITE(LUPRI,'(/A,I5)')'Left transform H2 of symmetry',
     *              IH2SYM
                CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A,I5)') 'with ZYMAT1 of symmetry',ISYM1
                CALL OUTPUT(ZYMAT1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A)') 'to H2X '
                CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
C
C Transform left (~~| )
C
               IH2SY1 = MULD2H(ICDSYM,INTSY1)
               INTSY2 = MULD2H(INTSY1,ISYM2)
               IF (ICDTYP.NE.JTSESE) THEN
                  CALL DZERO(H2,N2ORBX)
                  CALL OITH1(ISYM2,ZYMAT2,H2X,H2,IH2SY1)
                  IF (IPRFX.GT.300) THEN
                     WRITE(LUPRI,'(/A,I5)')
     *                   'Left transform H2X of symmetry',IH2SY1
                      CALL OUTPUT(H2X,1,NORBT,1,NORBT,
     *                   NORBT,NORBT,1,LUPRI)
                      WRITE(LUPRI,'(/A,I5)')
     *                   'with ZYMAT2 of symmetry',ISYM2
                      CALL OUTPUT(ZYMAT2,1,NORBT,1,NORBT,
     *                   NORBT,NORBT,1,LUPRI)
                      WRITE(LUPRI,'(/A)') 'to H2'
                      CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                            LUPRI)
                  END IF
                  ISP12=MULSP(ISPIN1,ISPIN2)
                  IF (ISP12.NE.IGRSPI) CALL QUIT('RSPFXD:SPIN ERROR')
                  IPDA = IPDENS(MULSP(IGRSPI,ISP12),0)
                  IPDB = IPDENS(ISP12,IGRSPI)
                  CALL EXTRAC(H2,WRK(KH2XAC))
                  CALL FQDIS1(INTSY2,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                        FI,FA,QA,QB,H2,WRK(KH2XAC),H2AX,
     *                        DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                        LORB,LCON,
     *                        IGRSPI,ISP12,0,IPERCD)
               END IF
C
C Transform right (~|~)
C
               IPDA = IPDENS(MULSP(IGRSPI,ISPIN1),ISPIN2)
               IPDB = IPDENS(ISPIN1,MULSP(IGRSPI,ISPIN2))
               CALL EXTRAC(H2X,WRK(KH2XAC))
               CALL FQDIS2(INTSY1,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                     FI,FA,QA,QB,H2X,WRK(KH2XAC),H2AX(1,2),
     *                     ZYMAT2,ISYM2,
     *                     DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                     LORB,LCON,
     *                     IGRSPI,ISPIN1,ISPIN2,IPERCD,WRK(KSCR))
         ELSE IF (IKLVL.EQ.3) THEN
C
C Transform left (~~| ) with k1,k2
C
            IH2SYM = MULD2H(ICDSYM,INTSYM)
C
C Transform left (~| ) with k1
C
            CALL DZERO(WRK(KH2),N2ORBX)
            CALL OITH1(ISYM1,ZYMAT1,H2,WRK(KH2),IH2SYM)
            INTSY1 = MULD2H(INTSYM,ISYM1)
            IF (IPRFX.GT.300) THEN
                WRITE(LUPRI,'(/A,I5)')'Left transform H2 of symmetry',
     *              IH2SYM
                CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A,I5)') 'with ZYMAT1 of symmetry',ISYM1
                CALL OUTPUT(ZYMAT1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A)') 'to H2~ '
                CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                      LUPRI)
            END IF
            IH2SY1 = MULD2H(ICDSYM,INTSY1)
C
C Transform left (~~| ) with k2
C
            CALL DZERO(H2X,N2ORBX)
            CALL OITH1(ISYM2,ZYMAT2,WRK(KH2),H2X,IH2SY1)
            INTSY2 = MULD2H(INTSY1,ISYM2)
                  IF (IPRFX.GT.300) THEN
                     WRITE(LUPRI,'(/A,I5)')
     *                   'Left transform H2~ of symmetry',IH2SY1
                      CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,
     *                   NORBT,NORBT,1,LUPRI)
                      WRITE(LUPRI,'(/A,I5)')
     *                   'with ZYMAT2 of symmetry',ISYM2
                      CALL OUTPUT(ZYMAT2,1,NORBT,1,NORBT,
     *                   NORBT,NORBT,1,LUPRI)
                      WRITE(LUPRI,'(/A)') 'to H2X'
                      CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,
     *                            1,LUPRI)
                  END IF
C
C Transform left (~~~| ) with k3
C
            IF (ICDTYP.NE.JTSESE) THEN
               IH2SY2 = MULD2H(ICDSYM,INTSY2)
               CALL DZERO(WRK(KH2),N2ORBX)
               CALL OITH1(ISYM3,ZYMAT3,H2X,WRK(KH2),IH2SY2)
               INTSY3=MULD2H(INTSY2,ISYM3)
               IF (IPRFX.GT.300) THEN
                  WRITE(LUPRI,'(/A,I5)')
     *                'Left transform H2~~ of symmetry',IH2SY2
                   CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,
     *                NORBT,NORBT,1,LUPRI)
                   WRITE(LUPRI,'(/A,I5)')
     *                'with ZYMAT3 of symmetry',ISYM3
                   CALL OUTPUT(ZYMAT2,1,NORBT,1,NORBT,
     *                NORBT,NORBT,1,LUPRI)
                   WRITE(LUPRI,'(/A)') 'to H2~~~'
                   CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,NORBT,NORBT,
     *                         1,LUPRI)
               END IF
               ISP12=MULSP(ISPIN1,ISPIN2)
               ISP123=MULSP(ISP12,ISPIN3)
               IF (ISP123.NE.IGRSPI) CALL QUIT('RSPFXD:SPIN ERROR')
               IPDA = IPDENS(MULSP(IGRSPI,ISP123),0)
               IPDB = IPDENS(ISP123,IGRSPI)
               CALL EXTRAC(WRK(KH2),WRK(KH2XAC))
               CALL FQDIS1(INTSY3,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                     FI,FA,QA,QB,WRK(KH2),WRK(KH2XAC),H2AX,
     *                     DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                     LORB,LCON,
     *                     IGRSPI,ISP123,0,IPERCD)
            END IF
C
C Transform right (~~|~) k1k2,k3
C
            ISP12=MULSP(ISPIN1,ISPIN2)
            IPDA = IPDENS(MULSP(IGRSPI,ISP12),ISPIN3)
            IPDB = IPDENS(ISP12,MULSP(IGRSPI,ISPIN3))
            CALL EXTRAC(H2X,WRK(KH2XAC))
            CALL FQDIS2(INTSY2,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                  FI,FA,QA,QB,H2X,WRK(KH2XAC),H2AX(1,2),
     *                  ZYMAT3,ISYM3,
     *                  DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                  LORB,LCON,
     *                  IGRSPI,ISP12 ,ISPIN3,IPERCD,WRK(KSCR))
C
C Transform (~~|~) k1k3,k2
C
C
C Transform left (~| ) with k1
C
            CALL DZERO(WRK(KH2),N2ORBX)
            CALL OITH1(ISYM1,ZYMAT1,H2,WRK(KH2),IH2SYM)
            INTSY1 = MULD2H(INTSYM,ISYM1)
            IF (IPRFX.GT.300) THEN
                WRITE(LUPRI,'(/A,I5)')'Left transform H2 of symmetry',
     *              IH2SYM
                CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A,I5)') 'with ZYMAT1 of symmetry',ISYM1
                CALL OUTPUT(ZYMAT1,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A)') 'to H2~ '
                CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                      LUPRI)
            END IF
            IH2SY1 = MULD2H(ICDSYM,INTSY1)
C
C Transform left (~~| ) with k3
C
            CALL DZERO(H2X,N2ORBX)
            CALL OITH1(ISYM3,ZYMAT3,WRK(KH2),H2X,IH2SY1)
            INTSY2 = MULD2H(INTSY1,ISYM3)
                  IF (IPRFX.GT.300) THEN
                     WRITE(LUPRI,'(/A,I5)')
     *                   'Left transform H2~ of symmetry',IH2SY1
                      CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,
     *                   NORBT,NORBT,1,LUPRI)
                      WRITE(LUPRI,'(/A,I5)')
     *                   'with ZYMAT3 of symmetry',ISYM3
                      CALL OUTPUT(ZYMAT3,1,NORBT,1,NORBT,
     *                   NORBT,NORBT,1,LUPRI)
                      WRITE(LUPRI,'(/A)') 'to H2X'
                      CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,
     *                            1,LUPRI)
                  END IF
C
C Transform right (~~|~) k2
C
            ISP13=MULSP(ISPIN1,ISPIN3)
            IPDA = IPDENS(MULSP(IGRSPI,ISP13),ISPIN2)
            IPDB = IPDENS(ISP13,MULSP(IGRSPI,ISPIN2))
            CALL EXTRAC(H2X,WRK(KH2XAC))
            CALL FQDIS2(INTSY2,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                  FI,FA,QA,QB,H2X,WRK(KH2XAC),H2AX(1,3),
     *                  ZYMAT2,ISYM2,
     *                  DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                  LORB,LCON,
     *                  IGRSPI,ISP13 ,ISPIN2,IPERCD,WRK(KSCR))
C
C Transform (~~|~) k2k3,k1
C
C
C Transform left (~| ) with k2
C
            CALL DZERO(WRK(KH2),N2ORBX)
            CALL OITH1(ISYM2,ZYMAT2,H2,WRK(KH2),IH2SYM)
            INTSY1 = MULD2H(INTSYM,ISYM2)
            IF (IPRFX.GT.300) THEN
                WRITE(LUPRI,'(/A,I5)')'Left transform H2 of symmetry',
     *              IH2SYM
                CALL OUTPUT(H2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A,I5)') 'with ZYMAT2 of symmetry',ISYM2
                CALL OUTPUT(ZYMAT2,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
                WRITE(LUPRI,'(/A)') 'to H2~ '
                CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,NORBT,NORBT,1,
     &                      LUPRI)
            END IF
            IH2SY1 = MULD2H(ICDSYM,INTSY1)
C
C Transform left (~~| ) with k3
C
            CALL DZERO(H2X,N2ORBX)
            CALL OITH1(ISYM3,ZYMAT3,WRK(KH2),H2X,IH2SY1)
            INTSY2 = MULD2H(INTSY1,ISYM3)
                  IF (IPRFX.GT.300) THEN
                     WRITE(LUPRI,'(/A,I5)')
     *                   'Left transform H2~ of symmetry',IH2SY1
                      CALL OUTPUT(WRK(KH2),1,NORBT,1,NORBT,
     *                   NORBT,NORBT,1,LUPRI)
                      WRITE(LUPRI,'(/A,I5)')
     *                   'with ZYMAT3 of symmetry',ISYM3
                      CALL OUTPUT(ZYMAT3,1,NORBT,1,NORBT,
     *                   NORBT,NORBT,1,LUPRI)
                      WRITE(LUPRI,'(/A)') 'to H2X'
                      CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,
     *                            1,LUPRI)
                  END IF
C
C Transform right (~~|~) k1
C
            ISP23=MULSP(ISPIN2,ISPIN3)
            IPDA = IPDENS(MULSP(IGRSPI,ISP23),ISPIN1)
            IPDB = IPDENS(ISP23,MULSP(IGRSPI,ISPIN1))
            CALL EXTRAC(H2X,WRK(KH2XAC))
            CALL FQDIS2(INTSY2,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                  FI,FA,QA,QB,H2X,WRK(KH2XAC),H2AX(1,4),
     *                  ZYMAT1,ISYM1,
     *                  DEN1,DEN2(IPDA),DEN2(IPDB),ISYMDN,PVDEN,
     *                  LORB,LCON,
     *                  IGRSPI,ISP23 ,ISPIN1,IPERCD,WRK(KSCR))
C           END IF
         ELSE
            WRITE(LUPRI,'(/A)') ' WRONG VALUE OF IKLVL IN RSPFXD'
            WRITE(LUPRI,'(/A,I5)') ' IKLVL =',IKLVL
            CALL QUIT(' WRONG VALUE OF IKLVL IN RSPFXD')
         END IF
      ELSE
         CALL QUIT(
     &        'RSPFXD not implemented for pre-transformed integrals')
      END IF
C
C all for now
C
C
         IF( IPRFX .GT. 300 ) THEN
            WRITE(LUPRI,'(/A)') ' PARTIAL MATRICES IN  RSPFXD'
            WRITE(LUPRI,'(A)')  ' ==========================='
            WRITE(LUPRI,'(//A)') ' Inactive fock matrix'
            CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Active fock matrix'
            CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Qa matrix'
            CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(//A)') ' Qb matrix'
            CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            WRITE(LUPRI,'(/A)') ' Active two-electron integrals'
            CALL PRIAC2(H2AX,NASHT,LUPRI)
            IF (IKLVL.EQ.2) THEN
               WRITE(LUPRI,'(/A)') ' Active two-electron integrals'
               CALL PRIAC2(H2AX(1,2),NASHT,LUPRI)
            END IF
         END IF
C
C
C     We have now completed processing this load, so get the next
C
      GO TO 1000
2000  CONTINUE
C
C Account for 1/2 in definition of two-electron integrals
C
C and for the fact that we calculated 2 times FA
C
      IF (IKLVL.EQ.0) THEN
         CALL DAXPY(N2ORBX,DH,WRK(KFI),1,FI,1)
         IF (NASHT.GT.0) CALL DSCAL(N2ORBX,DH*DH,FA,1)
         CALL DSCAL(NORBT*NASHT,DH,QA,1)
         CALL DSCAL(NORBT*NASHT,DH,QB,1)
      ELSE
         IF (NASHT.GT.0) CALL DSCAL(NORBT*NORBT,DH,FA,1)
      END IF
C
C Set DISTYP for the cases we may have ci-gradients
C
      IF (IKLVL.EQ.0) THEN
         INFDIS(1) = 9
         INFDIS(2) = 10
      END IF
      IF (IKLVL.EQ.1) THEN
         INFDIS(1) = 15
         INFDIS(2) = 16
      END IF
      IF (IKLVL.EQ.2) THEN
         INFDIS(1) = 11
         INFDIS(2) = 12
      END IF
      IF (IKLVL.EQ.3) THEN
         INFDIS(1) = 27
         INFDIS(2) = 28
      END IF
      IF (IPRFX.GT.200) THEN
         WRITE(LUPRI,'(/A)') ' Completed matrices in RSPFXD'
         WRITE(LUPRI,'(/A)') ' Inactive Fock matrix'
         CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' Active Fock matrix'
         CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' QA matrix'
         CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' QB matrix'
         CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' Active two-electron integrals'
         CALL PRIAC2(H2AX,NASHT,LUPRI)
         IF (IKLVL.EQ.2) THEN
            WRITE(LUPRI,'(/A)') ' Active two-electron integrals'
            CALL PRIAC2(H2AX(1,2),NASHT,LUPRI)
         END IF
      END IF
c
      CALL MEMREL('RSPFXD',WRK,KFIRST,KFIRST,KWRK,LWRK)
      CALL QEXIT('RSPFXD')
      RETURN
      END
C  /* Deck fqdis1 */
      SUBROUTINE FQDIS1(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                  FI,FA,QA,QB,H2X,H2XAC,H2ACAC,
     *                  DEN1,DEN2A,DEN2B,ISYMDN,PVDEN,
     *                  LORB,LCON,
     *                  IGRSPI,ISPIN1,ISPIN2,IPERCD)
C
C Olav Vahtras
C Dec 9 1991
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infrsp.h>
#include <infpri.h>
#include <orbtypdef.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infspi.h>
      LOGICAL LORB,LCON
      DIMENSION FI(*),FA(*),QA(*),QB(*),H2X(*),H2XAC(*),H2ACAC(*)
      DIMENSION DEN1(*),DEN2A(*),DEN2B(*),PVDEN(*)
C
      PARAMETER(DM1 = -1.0D0)
C
      CALL QENTER('FQDIS1')
      CALL FIXAD1(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *            FI,H2X,ISPIN1,ISPIN2)
      IF (NASHT.GT.0 .AND. LORB) THEN
         CALL FAXAD1(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *               DEN1,ISYMDN,FA,H2X,IGRSPI,ISPIN1,ISPIN2,H2XAC)
         CALL QXAD1(INTSYM,IC,ID,ICDTYP,
     *              ICDSYM,ICSYM,IDSYM,QA,QB, H2X,
     *              DEN2A,DEN2B,ISYMDN,PVDEN,H2XAC)
      END IF
      IF (LCON) THEN
         CALL DISAC1(IC,ID,H2XAC,H2ACAC)
      END IF
      IF (IC.NE.ID .AND. IPERCD.NE.0) THEN
         IF (IPERCD.EQ.1) THEN
C  nothing
         ELSEIF (IPERCD.EQ.-1) THEN
            CALL DSCAL(N2ORBX,DM1,H2X,1)
            CALL DSCAL(N2ASHX,DM1,H2XAC,1)
         ELSE
            WRITE(LUPRI,'(/A)') ' FQDIS1: WRONG VALUE OF IPERCD',IPERCD
            WRITE(LUPRI,'(/A)') ' ALLOWED IPERCD = 1 0 -1'
            CALL QUIT('FQDIS1: WRONG VALUE OF IPERCD')
         END IF
         CALL FIXAD1(INTSYM,ID,IC,ICDTYP,ICDSYM,IDSYM,ICSYM,
     *               FI,H2X,ISPIN1,ISPIN2)
         IF (NASHT.GT.0 .AND. LORB) THEN
            CALL FAXAD1(INTSYM,ID,IC,ICDTYP,ICDSYM,IDSYM,ICSYM,
     *                  DEN1,ISYMDN,FA,H2X,IGRSPI,ISPIN1,ISPIN2,
     *                  H2XAC)
            CALL QXAD1(INTSYM,ID,IC,ICDTYP,
     *                 ICDSYM,IDSYM,ICSYM,QA,QB, H2X,
     *                 DEN2A,DEN2B,ISYMDN,PVDEN,H2XAC)
         END IF
         IF (LCON) THEN
            CALL DISAC1(ID,IC,H2XAC,H2ACAC)
         END IF
         IF (IPERCD.EQ.-1) THEN
            CALL DSCAL(N2ORBX,DM1,H2X,1)
            CALL DSCAL(N2ASHX,DM1,H2XAC,1)
         END IF
      END IF
      CALL QEXIT('FQDIS1')
      RETURN
      END
C  /* Deck fqdis2 */
      SUBROUTINE FQDIS2(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                  FI,FA,QA,QB,H2X,H2XAC,H2ACAC,
     *                  ZYMAT,IZYSYM,
     *                  DEN1,DEN2A,DEN2B,ISYMDN,PVDEN,
     *                  LORB,LCON,
     *                  IGRSPI,ISPIN1,ISPIN2,IPERCD,SCR)
C
C Olav Vahtras
C Dec 9 1991
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infind.h>
#include <infdim.h>
#include <infrsp.h>
#include <infpri.h>
#include <orbtypdef.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infspi.h>
      LOGICAL LORB,LCON
      DIMENSION FI(*),FA(*),QA(*),QB(*)
      DIMENSION H2X(*),H2XAC(*),H2ACAC(*),ZYMAT(*)
      DIMENSION DEN1(*),DEN2A(*),DEN2B(*),PVDEN(*)
      DIMENSION SCR(*)
C
      PARAMETER(DM1 = -1.0D0)
C
      CALL QENTER('FQDIS2')
      ICTYP = IOBTYP(IC)
      IDTYP = IOBTYP(ID)
      CALL FIXAD2(INTSYM,IC,ID,ICDTYP,ICTYP,IDTYP,
     *            ICDSYM,ICSYM,IDSYM,
     *            FI,H2X,ISPIN1,ISPIN2,ZYMAT,IZYSYM)
      IF (NASHT.GT.0 .AND. LORB) THEN
         CALL FAXAD2(INTSYM,IC,ID,ICDTYP,ICTYP,IDTYP,
     *               ICDSYM,ICSYM,IDSYM,
     *               DEN1,ISYMDN,FA,H2X,IGRSPI,ISPIN1,ISPIN2,
     *               ZYMAT,IZYSYM,
     *               H2XAC,SCR)
         CALL QXAD2(INTSYM,IC,ID,ICDTYP,
     *              ICDSYM,ICSYM,IDSYM,QA,QB, H2X,DEN2A,DEN2B,
     *              ISYMDN,PVDEN,
     *              H2XAC,SCR,ZYMAT,IZYSYM)
      END IF
      IF (LCON) THEN
         CALL DISAC2(IC,ID,H2XAC,H2ACAC,ZYMAT,IZYSYM)
      END IF
      IF (IC.NE.ID .AND. IPERCD.NE.0) THEN
         IF (IPERCD.EQ.1) THEN
C  nothing
         ELSEIF (IPERCD.EQ.-1) THEN
            CALL DSCAL(N2ORBX,DM1,H2X,1)
            CALL DSCAL(N2ASHX,DM1,H2XAC,1)
         ELSE
            WRITE(LUPRI,'(/A)') ' FQDIS2: WRONG VALUE OF IPERCD',IPERCD
            WRITE(LUPRI,'(/A)') ' ALLOWED IPERCD = 1 0 -1'
            CALL QUIT('FQDIS2: WRONG VALUE OF IPERCD')
         END IF
         CALL FIXAD2(INTSYM,ID,IC,ICDTYP,IDTYP,ICTYP,
     *               ICDSYM,IDSYM,ICSYM,
     *               FI,H2X,ISPIN1,ISPIN2,ZYMAT,IZYSYM)
         IF (NASHT.GT.0 .AND. LORB) THEN
            CALL FAXAD2(INTSYM,ID,IC,ICDTYP,IDTYP,ICTYP,
     *                  ICDSYM,IDSYM,ICSYM,
     *                  DEN1,ISYMDN,FA,H2X,IGRSPI,ISPIN1,ISPIN2,
     *                  ZYMAT,IZYSYM,
     *                  H2XAC,SCR)
            CALL QXAD2(INTSYM,ID,IC,ICDTYP,
     *                  ICDSYM,IDSYM,ICSYM,QA,QB, H2X,DEN2A,DEN2B,
     *                  ISYMDN,PVDEN,
     *                  H2XAC,SCR,ZYMAT,IZYSYM)
         END IF
         IF (LCON) THEN
            CALL DISAC2(ID,IC,H2XAC,H2ACAC,ZYMAT,IZYSYM)
         END IF
         IF (IPERCD.EQ.-1) THEN
            CALL DSCAL(N2ORBX,DM1,H2X,1)
            CALL DSCAL(N2ASHX,DM1,H2XAC,1)
         END IF
      END IF
      CALL QEXIT('FQDIS2')
      RETURN
      END
C  /* Deck fixad1 */
      SUBROUTINE FIXAD1(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                 FI,H2X,ISPIN1,ISPIN2)
C
C
C     Olav Vahtras
C
C     Add contribution presently in H2X  to the inactive Fock matrix.
C
C     F(p,q) = (1+S1)(kk|pq) + (1+S2)(pq|kk) + (pk|kq) + (kq|pq)
C
C     Convention used in this routine: Distribution CD is the matrix
C                                      HCD(A,B)=(AB|CD)
C
C Case 1a)
C        F(p,q) = F(p,q) + (1+S1) (j j|p q)
C        F(C,D) = F(C,D) + (1+S1) (j j|C D)
C
C Case 1b)
C        F(p,q) = F(p,q) + (1+S2) (p q|j j)   C=D=inactive
C        F(p,q) = F(p,q) + (1+S2) (p q|C D)
C
C Case 2)
C        F(p,q) = F(p,q) -  (pj|jq)      C inactive
C        F(p,D) = F(p,D) -  (pC|CD)
C
C Case 3)
C        F(p,q) = F(p,q) - (j q|p j)     D inactive
C        F(C,q) = F(C,q) - (D q|C D)
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, D2= 2.0D0, DM1 = -1.0D0 )
C
#include <dftcom.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infdim.h>
#include <infind.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION FI(NORBT,NORBT), H2X(NORBT,NORBT)
      CALL QENTER('FIXAD1')
      IF (ICDTYP.EQ.JTSESE) GOTO 999
C
      IF (IPRRSP.GT.300) THEN
         WRITE(LUPRI,'(/A)') 'ENTRE FIXAD1'
         WRITE(LUPRI,'(A,I5)') 'INTSYM =',INTSYM
         WRITE(LUPRI,'(A,2I5)') 'IC ID =',IC,ID
         WRITE(LUPRI,'(A,I5)') 'ICDTYP =',ICDTYP
         WRITE(LUPRI,'(A,I5)') 'ICDSYM =',ICDSYM
         WRITE(LUPRI,'(A,I5)') 'ICSYM =',ICSYM
         WRITE(LUPRI,'(A,I5)') 'IDSYM =',IDSYM
         WRITE(LUPRI,'(A,I5)') 'ISPIN1 =',ISPIN1
         WRITE(LUPRI,'(A,I5)') 'ISPIN2 =',ISPIN2
         WRITE(LUPRI,'(A)') 'FI'
         CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'H2X'
         CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
      ICTYP = IOBTYP(IC)
      IDTYP = IOBTYP(ID)
C
C
C     Process case 1a)  : all distributions
C
C        F(p,q) = F(p,q) + (1+S1) (j j|p q)
C        F(C,D) = F(C,D) + (1+S1) (j j|C D)
C
      IF (ISPIN1.NE.1) THEN
         DO 10 I=1,NISHT
            IX=ISX(I)
            FI(IC,ID) = FI(IC,ID) + D2*H2X(IX,IX)
10       CONTINUE
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FIXAD1: FI after case 1b)'
            CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
      END IF
C
C     Process case 1b)  inactive - inactive distributions coulomb part
C
C        F(p,q) = F(p,q) + (1+S2) (p q|j j)
C        F(p,q) = F(p,q) + (1+S2) (p q|C D)
C
      IF (( IC .EQ. ID) .AND. (ICDTYP .EQ. JTININ) .AND. (ISPIN2.NE.1))
     *THEN
         CALL DAXPY(NORBT*NORBT,D2,H2X,1,FI,1)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FIXAD1: FI after case 1a)'
            CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
      END IF
C
C     Now process all exchange contributions
C
C     Process case 2)   :  inactive - general distribution
C
C        F(p,q) = F(p,q) -  (pj|jq)
C        F(p,D) = F(p,D) -  (pC|CD)
C
      IF (ICTYP.EQ.JTINAC .AND. HFXFAC.NE.D0) THEN
         IPSYM  = MULD2H(IDSYM,INTSYM)
         IF (IDTYP.EQ.JTSEC) THEN
            NORBP = NOCC(IPSYM)
         ELSE
            NORBP = NORB(IPSYM)
         END IF
         IPST = IORB(IPSYM) + 1
         IF(NORBP .GT. 0) THEN
            CALL DAXPY(NORBP,-HFXFAC,H2X(IPST,IC),1,
     *                               FI(IPST,ID),1)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FIXAD1: FI after case 3)'
               CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
      END IF
C
C     Process case 3)    : general - inactive distribution
C
C        F(p,q) = F(p,q) - (j q|p j)
C        F(C,q) = F(C,q) - (D q|C D)
C
      IF (IDTYP.EQ.JTINAC .AND. HFXFAC.NE.D0) THEN
         IQSYM   = MULD2H(ICSYM,INTSYM)
         IF (ICTYP.EQ.JTSEC) THEN
            NORBQ  =  NOCC(IQSYM)
         ELSE
            NORBQ  =  NORB(IQSYM)
         END IF
         IQST  =   IORB(IQSYM) + 1
         IF( NORBQ .GT. 0) THEN
            CALL DAXPY(NORBQ,-HFXFAC,H2X(ID,IQST),NORBT,
     *                               FI(IC,IQST),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FIXAD1: FI after case 2)'
               CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
      END IF
999   CALL QEXIT('FIXAD1')
      RETURN
      END
C  /* Deck fixad2 */
      SUBROUTINE FIXAD2(INTSYM,IC,ID,ICDTYP,ICTYP,IDTYP,
     *                 ICDSYM,ICSYM,IDSYM,
     *                 FI,H2X,ISPIN1,ISPIN2,ZYMAT,IZYSYM)
C
C     Olav Vahtras
C     Nov 22 1991
C
C     Purpose: One-index transform the right-hand side of (AB|CD)
C     i.e. (AB|CD) -> (AB|C~D)
C     and distribute in the inactive Fock matrix
C     Convention used in this routine: Distribution CD is the matrix
C                                      HCD(A,B)=(AB|CD)
C
C     Expression:
C
C     F(p,q) = (1+S1)(kk|p~q) + (1+S2)(pq|k~k) - (pk|k~q) - (kq|p~k)
C
C     1) (kk|p~q) = ZYMAT(p,r)(kk|rq) - (kk|pr)ZYMAT(r,q)
C
C     Conditions: S1 = +, ICDSYM = INTSYM
C
C     1a) - 2(kk|pr)ZYMAT(r,q) -> F(p,q)
C         - 2(kk|CD)ZYMAT(D,q) -> F(C,q)
C
C         Restrictions: IQSYM=IFCKSY X ICSYM
C
C     1b) 2*ZYMAT(p,r)(kk|rq) -> F(p,q)
C         2*ZYMAT(p,C)(kk|CD) -> F(p,D)
C
C          Restrictions: IPSYM=IFCKSY X IDSYM
C
C     2)  (pq|k~k) = ZYMAT(k,r)(pq|rk) - (pq|rk)ZYMAT(k,r)
C
C         Conditions: S2 = +
C                     ICDTYP not inactive-inactive or secondary-secondary
C
C     2a) 2*ZYMAT(k,r)(pq|rk) -> F(p,q)
C         2*ZYMAT(D,C)(pq|CD) -> F(p,q)             D Inactive
C
C     2b) -2*(pq|rk)ZYMAT(k,r) -> F(p,q)
C         -2*(pq|CD)ZYMAT(D,C) -> F(p,q)            C Inactive
C
C     3) -(pk|k~q) = (pk|kr)ZYMAT(r,q) - ZYMAT(k,r)(pk|rq)
C
C     3a) (pk|kr)ZYMAT(r,q) -> F(p,q)
C         (pC|CD)ZYMAT(D,q) -> F(p,q)               C Inactive
C
C         Restrictions: IPSYM=IDSYM X INTSYM
C                       IQSYM=IDSYM X IZYSYM
C
C     3b) -ZYMAT(k,r)(pk|rq) -> F(p,q)
C         -ZYMAT(k,C)(pk|CD) -> F(p,D)             C not inactive
C
C         Restrictions: IPSYM=IDSYM X IFCKSY
C                       IKSYM=ICSYM X IZYSYM
C
C     4) -(kq|p~k) = (kq|pr)ZYMAT(r,k) - ZYMAT(p,r)(kq|rk)
C
C     4a) (kq|pr)ZYMAT(r,k) -> F(p,q)
C         (kq|CD)ZYMAT(D,k) -> F(C,q)              D not inactive
C
C         Restrictions: IQSYM = ICSYM X IFCKSY
C                       IKSYM = IDSYM X IZYSYM
C
C     4b) -ZYMAT(p,r)(kq|rk) -> F(p,q)
C         -ZYMAT(p,C)(Dq|CD) -> F(p,q)             D inactive
C
C         Restrictions: IPSYM = ICSYM X IZYSYM
C                       IQSYM = ICSYM X INTSYM
C
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0, DM1 = -1.0D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infdim.h>
#include <infind.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
#include <dftcom.h>
C
      DIMENSION FI(NORBT,NORBT), H2X(NORBT,NORBT)
      DIMENSION ZYMAT(NORBT,NORBT)
      CALL QENTER('FIXAD2')
C
      IFCKSY = MULD2H(INTSYM,IZYSYM)
      IF (IPRRSP.GT.300) THEN
         WRITE(LUPRI,'(/A)') 'ENTRE FIXAD2'
         WRITE(LUPRI,'(A,I5)') 'INTSYM =',INTSYM
         WRITE(LUPRI,'(A,2I5)') 'IC ID =',IC,ID
         WRITE(LUPRI,'(A,I5)') 'ICDTYP =',ICDTYP
         WRITE(LUPRI,'(A,I5)') 'ICDSYM =',ICDSYM
         WRITE(LUPRI,'(A,I5)') 'ICSYM =',ICSYM
         WRITE(LUPRI,'(A,I5)') 'IDSYM =',IDSYM
         WRITE(LUPRI,'(A,I5)') 'ISPIN1 =',ISPIN1
         WRITE(LUPRI,'(A,I5)') 'ISPIN2 =',ISPIN2
         WRITE(LUPRI,'(A,I5)') 'IZYSYM =',IZYSYM
         WRITE(LUPRI,'(A)') 'FI'
         CALL OUTPUT(FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'H2X'
         CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'ZYMAT'
         CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
      END IF
C
C Case 1)
C
      IF (ISPIN1.EQ.0 .AND. ICDSYM.EQ.INTSYM) THEN
C
C Evaluate inactive trace
C
         DTRACE = D0
         DO 10 I=1,NISHT
            IX=ISX(I)
            DTRACE = DTRACE + H2X(IX,IX)
10       CONTINUE
         DTRACE = D2*DTRACE
C
C Case 1a)
C     1a) - 2(kk|pr)ZYMAT(r,q) -> F(p,q)
C         - 2(kk|CD)ZYMAT(D,q) -> F(C,q)
C
         IQSYM = MULD2H(IFCKSY,ICSYM)
         NORBQ = NORB(IQSYM)
        IF (NORBQ .GT. 0) THEN
         IQST = IORB(IQSYM) + 1
         CALL DAXPY(NORBQ,-DTRACE,ZYMAT(ID,IQST),NORBT,
     *                               FI(IC,IQST),NORBT)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FIXAD2:FI after case 1a'
            CALL OUTPUT (FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
        END IF
C
C Case 1b)
C         2*ZYMAT(p,r)(kk|rq) -> F(p,q)
C         2*ZYMAT(p,C)(kk|CD) -> F(p,D)
C
         IPSYM=MULD2H(IFCKSY,IDSYM)
         NORBP=NORB(IPSYM)
        IF (NORBP .GT. 0) THEN
         IPST=IORB(IPSYM)+1
         CALL DAXPY(NORBP,DTRACE,ZYMAT(IPST,IC),1,FI(IPST,ID),1)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FIXAD2:FI after case 1b'
            CALL OUTPUT (FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
        END IF
      END IF
C
C Case 2)
C
      IF (ISPIN2.EQ.0 .AND. ICDSYM.EQ.IZYSYM
     *   .AND. ICDTYP.NE.JTININ .AND. ICDTYP.NE.JTSESE) THEN
C
C Case 2a)
C         2*ZYMAT(k,r)(pq|rk) -> F(p,q)
C         2*ZYMAT(D,C)(pq|CD) -> F(p,q)             D Inactive
C
         ZYDC=2*ZYMAT(ID,IC)
         IF (IDTYP.EQ.JTINAC) THEN
            CALL DAXPY(N2ORBX,ZYDC,H2X,1,FI,1)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FIXAD2:FI after case 2a'
               CALL OUTPUT (FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
C
C Case 2b)
C         -2*(pq|rk)ZYMAT(k,r) -> F(p,q)
C         -2*(pq|CD)ZYMAT(D,C) -> F(p,q)            C Inactive
C
         IF (ICTYP.EQ.JTINAC) THEN
            CALL DAXPY(N2ORBX,-ZYDC,H2X,1,FI,1)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FIXAD2:FI after case 2b'
               CALL OUTPUT (FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
      END IF
C
C Case 3)
C
      IF (HFXFAC.NE.D0) THEN
      IF (ICTYP.EQ.JTINAC) THEN
C
C Case 3a)
C         (pk|kr)ZYMAT(r,q) -> F(p,q)
C         (pC|CD)ZYMAT(D,q) -> F(p,q)               C Inactive
C
         IPSYM=MULD2H(INTSYM,IDSYM)
         IQSYM=MULD2H(IZYSYM,IDSYM)
         NORBP=NORB(IPSYM)
         NORBQ=NORB(IQSYM)
        IF (NORBP .GT. 0 .AND. NORBQ .GT. 0) THEN
         IPST=IORB(IPSYM)+1
         IQST=IORB(IQSYM)+1
         CALL DGEMM('N','N',NORBP,NORBQ,1,HFXFAC,
     &              H2X(IPST,IC),NORBT,
     &              ZYMAT(ID,IQST),NORBT,1.D0,
     &              FI(IPST,IQST),NORBT)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FIXAD2:FI after case 3a'
            CALL OUTPUT (FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
        END IF
      ELSE
C
C Case 3b)
C         -ZYMAT(k,r)(pk|rq) -> F(p,q)
C         -ZYMAT(k,C)(pk|CD) -> F(p,D)             C not inactive
C
         IPSYM=MULD2H(IFCKSY,IDSYM)
         IKSYM=MULD2H(IZYSYM,ICSYM)
         NORBP=NORB(IPSYM)
         NISHK=NISH(IKSYM)
        IF (NORBP .GT. 0 .AND. NISHK .GT. 0) THEN
         IPST=IORB(IPSYM)+1
         IKST=IORB(IKSYM)+1
         CALL DGEMM('N','N',NORBP,1,NISHK,-HFXFAC,
     &              H2X(IPST,IKST),NORBT,
     &              ZYMAT(IKST,IC),NORBT,1.D0,
     &              FI(IPST,ID),NORBT)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FIXAD2:FI after case 3b'
            CALL OUTPUT (FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
        END IF
      END IF
C
C Case 4)
C
      IF (IDTYP.EQ.JTINAC) THEN
C
C Case 4b)
C         -ZYMAT(p,r)(kq|rk) -> F(p,q)
C         -ZYMAT(p,C)(Dq|CD) -> F(p,q)             D inactive
C
         IPSYM=MULD2H(IZYSYM,ICSYM)
         IQSYM=MULD2H(INTSYM,ICSYM)
         NORBP=NORB(IPSYM)
         NORBQ=NORB(IQSYM)
        IF (NORBP .GT. 0 .AND. NORBQ .GT. 0) THEN
         IPST=IORB(IPSYM)+1
         IQST=IORB(IQSYM)+1
         CALL DGEMM('N','N',NORBP,NORBQ,1,-HFXFAC,
     &              ZYMAT(IPST,IC),NORBT,
     &              H2X(ID,IQST),NORBT,1.D0,
     &              FI(IPST,IQST),NORBT)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FIXAD2:FI after case 4b'
            CALL OUTPUT (FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
        END IF
C
C Case 4a)
C         (kq|pr)ZYMAT(r,k) -> F(p,q)
C         (kq|CD)ZYMAT(D,k) -> F(C,q)              D not inactive
C
      ELSE
         IKSYM=MULD2H(IZYSYM,IDSYM)
         IQSYM=MULD2H(IFCKSY,ICSYM)
         NISHK=NISH(IKSYM)
         NORBQ=NORB(IQSYM)
        IF (NISHK .GT. 0 .AND. NORBQ .GT. 0) THEN
         IKST=IORB(IKSYM)+1
         IQST=IORB(IQSYM)+1
         CALL DGEMM('N','N',1,NORBQ,NISHK,HFXFAC,
     &              ZYMAT(ID,IKST),NORBT,
     &              H2X(IKST,IQST),NORBT,1.D0,
     &              FI(IC,IQST),NORBT)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FIXAD2:FI after case 4a'
            CALL OUTPUT (FI,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
        END IF
      END IF
      END IF ! if (hfxfac.ne.0.0d0)
      CALL QEXIT('FIXAD2')
      RETURN
      END
C  /* Deck faxad1 */
      SUBROUTINE FAXAD1(INTSYM,IC,ID,ICDTYP,ICDSYM,ICSYM,IDSYM,
     *                  DEN1,IDENSY,FA,H2X,IGRSPI,ISPIN1,ISPIN2,H2XAC)
C
C     Add contribution presently in H2X  to the active Fock matrix.
C
C
C F(p,q) = [(1+S1Sa)(pq|xy) + (1+S2Sa)(xy|pq) - (py|xq) - (xq|py)]D(x,y)
C
C     Convention used in this routine: Distribution CD is the matrix
C                                      HCD(A,B)=(AB|CD)
C Case 1a)
C        F(p,q) = F(p,q) + (1+S1Sa) (p q|x y) D(x,y)   C,D active
C        F(p,q) = F(p,q) + (1+S1Sa) (p q|C D) D(C,D)
C
C Case 1b)
C        F(p,q) = F(p,q) + (1+S2Sa) (x y|p q) D(x,y)
C        F(C,D) = F(C,D) + (1+S2Sa) (x y|C D) D(x,y)
C
C Case 2)
C        F(p,q) = F(p,q) -  (py|xq) D(x,y)      C active
C        F(p,D) = F(p,D) -  (py|CD) D(C,y)
C
C Case 3)
C        F(p,q) = F(p,q) - (x q|p y) D(x,y)     D active
C        F(C,q) = F(C,q) - (x q|C D) D(x,D)
C
#include <implicit.h>
C
      PARAMETER ( D0 = 0.0D0, D2= 2.0D0, DM1 = -1.0D0 , DMH = -.5D0)
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <wrkrsp.h>
#include <infrsp.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infdim.h>
#include <infind.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION DEN1(NASHT,NASHT),FA(NORBT,NORBT), H2X(NORBT,NORBT)
      DIMENSION H2XAC(NASHT,NASHT)
C     LOGICAL NATORB
C
      IF (ICDTYP.EQ.JTSESE) RETURN
      CALL QENTER('FAXAD1')
      IFCKSY = MULD2H(INTSYM,IDENSY)
      ICTYP = IOBTYP(IC)
      IDTYP = IOBTYP(ID)
      ICW = ISW( IC )
      IDW = ISW( ID )
      NCW = ICW - NISHT
      NDW = IDW - NISHT
      ISP1A = MULSP(ISPIN1,IGRSPI)
      ISP2A = MULSP(ISPIN2,IGRSPI)
      IF (IPRRSP.GT.300) THEN
         WRITE(LUPRI,'(/A)') 'ENTRE FAXAD1'
         WRITE(LUPRI,'(A,I5)') 'INTSYM =',INTSYM
         WRITE(LUPRI,'(A,I5)') 'IDENSY =',IDENSY
         WRITE(LUPRI,'(A,2I5)') 'IC ID =',IC,ID
         WRITE(LUPRI,'(A,I5)') 'ICDTYP =',ICDTYP
         WRITE(LUPRI,'(A,I5)') 'ICDSYM =',ICDSYM
         WRITE(LUPRI,'(A,I5)') 'ICSYM =',ICSYM
         WRITE(LUPRI,'(A,I5)') 'IDSYM =',IDSYM
         WRITE(LUPRI,'(A,I5)') 'IGRSPI =',IGRSPI
         WRITE(LUPRI,'(A,I5)') 'ISPIN1 =',ISPIN1
         WRITE(LUPRI,'(A,I5)') 'ISPIN2 =',ISPIN2
C        WRITE(LUPRI,'(A,L5)') 'NATORB =',NATORB
         WRITE(LUPRI,'(A)') 'FA'
         CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'DEN1'
         CALL OUTPUT(DEN1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'H2X'
         CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'H2XAC'
         CALL OUTPUT(H2XAC,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
C
C     Process case 1)  active - active distributions coulomb part
C
C Case 1a)
C        F(p,q) = F(p,q) + (1+S1Sa) (p q|x y) D(x,y)   C,D active
C        F(p,q) = F(p,q) + (1+S1Sa) (p q|C D) D(C,D)
C
      IF ((ICDTYP .EQ. JTACAC) .AND. (ISP1A.NE.1)) THEN
         DENFAC = D2*DEN1(NCW,NDW)
         IF (DENFAC .NE. D0) THEN
            CALL DAXPY(NORBT*NORBT,DENFAC,H2X,1,FA,1)
         END IF
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FAXAD1: FA after case 1a)'
            CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
      END IF
C
C Case 1b)
C        F(p,q) = F(p,q) + (1+S2Sa) (x y|p q) D(x,y)
C        F(C,D) = F(C,D) + (1+S2Sa) (x y|C D) D(x,y)
C
      IF (ISP2A.NE.1 .AND. ICDSYM.EQ.IFCKSY .AND.ICDTYP.NE.JTININ) THEN
         FA(IC,ID) = FA(IC,ID) + D2*DDOT(N2ASHX,H2XAC,1,DEN1,1)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FAXAD1: FA after case 1b)'
            CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
      END IF

C
C     Now process all exchange contributions
C
C Case 2)
C        F(p,q) = F(p,q) -  (py|xq) D(x,y)      C active
C        F(p,D) = F(p,D) -  (py|CD) D(C,y)
C
C
      IF (ICTYP.EQ.JTACT) THEN
         IPSYM  = MULD2H(IDSYM,IFCKSY)
         IYSYM = MULD2H(ICSYM,IDENSY)
         IF (IDTYP.EQ.JTSEC) THEN
            NORBP = NOCC(IPSYM)
         ELSE
            NORBP = NORB(IPSYM)
         END IF
         NASHY = NASH(IYSYM)
         IF(NORBP.GT.0 .AND. NASHY.GT.0) THEN
            IPST = IORB(IPSYM) + 1
            IYST = IORB(IYSYM) + NISH(IYSYM) + 1
            NYW = IASH(IYSYM) + 1
            CALL DGEMM('N','T',NORBP,1,NASHY,-1.D0,
     &                 H2X(IPST,IYST),NORBT,
     &                 DEN1(NCW,NYW),NASHT,1.D0,
     &                 FA(IPST,ID),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FAXAD1: FA after case 2)'
               CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
      END IF
C
C Case 3)
C        F(p,q) = F(p,q) - (x q|p y) D(x,y)     D active
C        F(C,q) = F(C,q) - (x q|C D) D(x,D)
C
      IF (IDTYP.EQ.JTACT) THEN
         IQSYM = MULD2H(ICSYM,IFCKSY)
         IXSYM = MULD2H(IDSYM,IDENSY)
         IF (ICTYP.EQ.JTSEC) THEN
            NORBQ =  NOCC(IQSYM)
         ELSE
            NORBQ =  NORB(IQSYM)
         END IF
         NASHX = NASH(IXSYM)
         IF (NORBQ.GT.0 .AND. NASHX.GT.0) THEN
            IQST = IORB(IQSYM) + 1
            IXST = IORB(IXSYM) + NISH(IXSYM) + 1
            NXW = IASH(IXSYM) + 1
            CALL DGEMM('T','N',1,NORBQ,NASHX,-1.D0,
     &                 DEN1(NXW,NDW),NASHT,
     &                 H2X(IXST,IQST),NORBT,1.D0,
     &                 FA(IC,IQST),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FAXAD1: FA after case 3)'
               CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
      END IF
      CALL QEXIT('FAXAD1')
      RETURN
      END
C  /* Deck faxad2 */
      SUBROUTINE FAXAD2(INTSYM,IC,ID,ICDTYP,ICTYP,IDTYP,
     *                 ICDSYM,ICSYM,IDSYM,
     *                 DEN1,IDENSY,FA,H2X,IGRSPI,ISPIN1,ISPIN2,
     *                 ZYMAT,IZYSYM, H2XAC,WRK)
C
C     Olav Vahtras
C
C     Purpose: One-index transform the right-hand side of (AB|CD)
C     i.e. (AB|CD) -> (AB|C~D)
C     and distribute in the inactive Fock matrix
C     Convention used in this routine: Distribution CD is the matrix
C                                      HCD(A,B)=(AB|CD)
C
C     Expression:
C
C F(p,q) = [(1+S2Sa)(xy|p~q)+(1+S1Sa)(pq|x~y)-(px|y~q)-(xq|p~y)] D(x,y)
C
C     1) (xy|p~q) = ZYMAT(p,r)(xy|rq) - (xy|pr)ZYMAT(r,q)
C
C     Conditions: S2Sa = +, ICDSYM =
C
C     1a) 2*ZYMAT(p,r)(xy|rq)D(x,y) -> F(p,q)
C         2*ZYMAT(p,C)(xy|CD)D(x,y) -> F(p,D)
C
C          Restrictions: IPSYM=IFCKSY X IDSYM
C
C     1b) - 2(xy|pr)ZYMAT(r,q) D(x,y) -> F(p,q)
C         - 2(xy|CD)ZYMAT(D,q) D(x,y) -> F(C,q)
C
C         Restrictions: IQSYM=IFCKSY X ICSYM
C
C     2)  (pq|x~y) = ZYMAT(x,r)(pq|ry) - (pq|xr)ZYMAT(r,y)
C
C         Conditions: S1Sa = +
C                     ICDTYP active-general
C
C     2a) 2*ZYMAT(x,r)(pq|ry) D(x,y) -> F(p,q)
C         2*ZYMAT(x,C)(pq|CD) D(x,D) -> F(p,q)             D Active
C
C     2b) -2*(pq|xr)ZYMAT(r,y) D(x,y) -> F(p,q)
C         -2*(pq|CD)ZYMAT(D,y) D(C,y) -> F(p,q)            C Active
C
C     3) -(py|x~q) = (py|xr)ZYMAT(r,q) - ZYMAT(x,r)(py|rq)
C
C     3a) (py|xr)ZYMAT(r,q) D(x,y) -> F(p,q)
C         (py|CD)ZYMAT(D,q) D(C,y) -> F(p,q)               C Active
C
C         Restrictions: IPSYM=IDENSY X IDSYM X INTSYM
C                       IQSYM=IDSYM X IZYSYM
C
C     3b) -ZYMAT(x,r)(py|rq) D(x,y) -> F(p,q)
C         -ZYMAT(x,C)(py|CD) D(x,y) -> F(p,D)
C
C         Restrictions: IPSYM=IDSYM X IFCKSY
C                       IXSYM=ICSYM X IZYSYM
C
C     4) -(xq|p~y) = (xq|pr)ZYMAT(r,y) - ZYMAT(p,r)(xq|ry)
C
C     4a) (xq|pr)ZYMAT(r,y) D(x,y) -> F(p,q)
C         (xq|CD)ZYMAT(D,y) D(x,y) -> F(C,q)
C
C         Restrictions: IQSYM = ICSYM X IFCKSY
C                       IYSYM = IDSYM X IZYSYM
C
C     4b) -ZYMAT(p,r)(xq|ry) D(x,y) -> F(p,q)
C         -ZYMAT(p,C)(xq|CD) D(x,D) -> F(p,q)             D active
C
C         Restrictions: IPSYM = ICSYM X IZYSYM
C                       IQSYM = ICSYM X INTSYM
C
C
#include <implicit.h>
C
      PARAMETER ( D2= 2.0D0, DM1 = -1.0D0 )
C
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infrsp.h>
#include <infdim.h>
#include <infind.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
C
      DIMENSION FA(NORBT,NORBT), H2X(NORBT,NORBT)
      DIMENSION ZYMAT(NORBT,NORBT)
      DIMENSION DEN1(NASHT,NASHT)
      DIMENSION H2XAC(NASHT,NASHT)
      DIMENSION WRK(*)
      CALL QENTER('FAXAD2')
C
      IF (IPRRSP.GT.300) THEN
         WRITE(LUPRI,'(/A)') 'ENTRE FAXAD2'
         WRITE(LUPRI,'(A,I5)') 'INTSYM =',INTSYM
         WRITE(LUPRI,'(A,2I5)') 'IC ID =',IC,ID
         WRITE(LUPRI,'(A,I5)') 'ICDTYP =',ICDTYP
         WRITE(LUPRI,'(A,I5)') 'ICDSYM =',ICDSYM
         WRITE(LUPRI,'(A,I5)') 'ICSYM =',ICSYM
         WRITE(LUPRI,'(A,I5)') 'IDSYM =',IDSYM
         WRITE(LUPRI,'(A,I5)') 'IGRSPI =',IGRSPI
         WRITE(LUPRI,'(A,I5)') 'ISPIN1 =',ISPIN1
         WRITE(LUPRI,'(A,I5)') 'ISPIN2 =',ISPIN2
         WRITE(LUPRI,'(A,I5)') 'IZYSYM =',IZYSYM
         WRITE(LUPRI,'(A)') 'FA'
         CALL OUTPUT(FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'H2X'
         CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'ZYMAT'
         CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'DEN1'
         CALL OUTPUT(DEN1,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
      IDZY = MULD2H(IDENSY,IZYSYM)
      IFCKSY=MULD2H(INTSYM,IDZY)
      ISP1A = MULSP(ISPIN1,IGRSPI)
      ISP2A = MULSP(ISPIN2,IGRSPI)
      ICW = ISW( IC )
      IDW = ISW( ID )
      NCW = ICW - NISHT
      NDW = IDW - NISHT
C
C Case 1)
C        (xy|p~q) = ZYMAT(p,r)(xy|rq) - (xy|pr)ZYMAT(r,q)
C
C
C Case 1a)
C         2*ZYMAT(p,r)(xy|rq)D(x,y) -> F(p,q)
C         2*ZYMAT(p,C)(xy|CD)D(x,y) -> F(p,D)
C
      IF (ISP2A.EQ.0 .AND. ICDSYM.EQ.MULD2H(IZYSYM,IFCKSY) )THEN
         DTRACE=D2*DDOT(N2ASHX,H2XAC,1,DEN1,1)
         IPSYM=MULD2H(IFCKSY,IDSYM)
         IF (ICTYP.EQ.JTSEC .OR. IDTYP.EQ.JTSEC) THEN
            NORBP=NOCC(IPSYM)
         ELSE
            NORBP=NORB(IPSYM)
         END IF
         IF (NORBP.GT.0) THEN
            IPST=IORB(IPSYM)+1
            CALL DAXPY(NORBP,DTRACE,ZYMAT(IPST,IC),1,FA(IPST,ID),1)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FAXAD2:FA after case 1a'
               CALL OUTPUT (FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
C
C Case 1b
C         - 2(xy|pr)ZYMAT(r,q) D(x,y) -> F(p,q)
C         - 2(xy|CD)ZYMAT(D,q) D(x,y) -> F(C,q)
C
C
         IQSYM = MULD2H(IFCKSY,ICSYM)
         IF (ICTYP.EQ.JTSEC .OR. IDTYP.EQ.JTSEC) THEN
            NORBQ = NOCC(IQSYM)
         ELSE
            NORBQ = NORB(IQSYM)
         END IF
         IF (NORBQ.GT.0) THEN
            IQST = IORB(IQSYM) + 1
            CALL DAXPY(NORBQ,-DTRACE,ZYMAT(ID,IQST),NORBT,
     *                                  FA(IC,IQST),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FAXAD2:FA after case 1b'
               CALL OUTPUT (FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
           END IF
         END IF
      END IF
C
C Case 2)
C
      IF (ISP1A.EQ.0 .AND. ICDSYM.EQ.MULD2H(IFCKSY,INTSYM)) THEN
C
C Case 2a)
C     2a) 2*ZYMAT(x,r)(pq|ry) D(x,y) -> F(p,q)
C         2*ZYMAT(x,C)(pq|CD) D(x,D) -> F(p,q)             D Active
C
         IF (IDTYP.EQ.JTACT ) THEN
            IXSYM=MULD2H(IZYSYM,ICSYM)
            NASHX=NASH(IXSYM)
            IF (NASHX.GT.0) THEN
               IXST = IORB(IXSYM) + NISH(IXSYM) + 1
               NXW = IASH(IXSYM) + 1
               ZYDEN = D2*DDOT(NASHX,ZYMAT(IXST,IC),1,DEN1(NXW,NDW),1)
               CALL DAXPY(N2ORBX,ZYDEN,H2X,1,FA,1)
               IF (IPRRSP.GT.300) THEN
                  WRITE(LUPRI,'(/A)') 'FAXAD2:FA after case 2a'
                  CALL OUTPUT (FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
            END IF
         END IF
C
C Case 2b)
C     2b) -2*(pq|xr)ZYMAT(r,y) D(x,y) -> F(p,q)
C         -2*(pq|CD)ZYMAT(D,y) D(C,y) -> F(p,q)            C Active
C
         IF (ICTYP.EQ.JTACT) THEN
            IYSYM = MULD2H(IZYSYM,IDSYM)
            NASHY = NASH(IYSYM)
            IF (NASHY.GT.0) THEN
               IYST = IORB(IYSYM) + NISH(IYSYM) + 1
               NYW = IASH(IYSYM) + 1
               ZYDEN = D2*DDOT(NASHY,ZYMAT(ID,IYST),NORBT,
     *                               DEN1(NCW,NYW),NASHT)
               CALL DAXPY(N2ORBX,-ZYDEN,H2X,1,FA,1)
               IF (IPRRSP.GT.300) THEN
                  WRITE(LUPRI,'(/A)') 'FAXAD2:FA after case 2b'
                  CALL OUTPUT (FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
               END IF
            END IF
         END IF
      END IF
C
C Case 3)
C        -(py|x~q) = (py|xr)ZYMAT(r,q) - ZYMAT(x,r)(py|rq)
C
C     3a) (py|xr)ZYMAT(r,q) D(x,y) -> F(p,q)
C         (py|CD)ZYMAT(D,q) D(C,y) -> F(p,q)               C Active
C
C
      IF (ICTYP.EQ.JTACT) THEN
         IQSYM=MULD2H(IZYSYM,IDSYM)
         IPSYM=MULD2H(IFCKSY,IQSYM)
         IYSYM=MULD2H(IDENSY,ICSYM)
         NORBP=NORB(IPSYM)
         IF (IDTYP.EQ.JTSEC) THEN
            NORBQ=NOCC(IQSYM)
         ELSE
            NORBQ=NORB(IQSYM)
         END IF
         NASHY=NASH(IYSYM)
         IF (NORBP.GT.0 .AND. NORBQ.GT.0 .AND.NASHY.GT.0) THEN
            IPST=IORB(IPSYM)+1
            IQST=IORB(IQSYM)+1
            IYST=IORB(IYSYM)+NISH(IYSYM)+1
            NYW = IASH(IYSYM)+1
            CALL DGEMM('N','T',NORBP,1,NASHY,1.D0,
     &                 H2X(IPST,IYST),NORBT,
     &                 DEN1(NCW,NYW),NASHT,0.D0,
     &                 WRK,NORBP)
            CALL DGEMM('N','N',NORBP,NORBQ,1,1.D0,
     &                 WRK,NORBP,
     &                 ZYMAT(ID,IQST),NORBT,1.D0,
     &                 FA(IPST,IQST),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FAXAD2:FA after case 3a'
               CALL OUTPUT (FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
      END IF
C
C Case 3b)
C         -ZYMAT(x,r)(py|rq) D(x,y) -> F(p,q)
C         -ZYMAT(x,C)(py|CD) D(x,y) -> F(p,D)
C
C
      IPSYM=MULD2H(IFCKSY,IDSYM)
      IXSYM=MULD2H(IZYSYM,ICSYM)
      IYSYM=MULD2H(IDENSY,IXSYM)
      IF (IDTYP.EQ.JTSEC) THEN
         NORBP=NOCC(IPSYM)
      ELSE
         NORBP=NORB(IPSYM)
      END IF
      NASHX=NASH(IXSYM)
      NASHY=NASH(IYSYM)
      IF (NORBP.GT.0 .AND. NASHX.GT.0 .AND. NASHY.GT.0) THEN
         IPST = IORB(IPSYM) + 1
         IXST = IORB(IXSYM) + NISH(IXSYM) + 1
         IYST = IORB(IYSYM) + NISH(IYSYM) + 1
         NXW = IASH(IXSYM) + 1
         NYW = IASH(IYSYM) + 1
         CALL DGEMM('T','N',1,NASHY,NASHX,1.D0,
     &              ZYMAT(IXST,IC),NORBT,
     &              DEN1(NXW,NYW),NASHT,0.D0,
     &              WRK,1)
         CALL DGEMM('N','N',NORBP,1,NASHY,-1.D0,
     &              H2X(IPST,IYST),NORBT,
     &              WRK,NASHY,1.D0,
     &              FA(IPST,ID),NORBT)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FAXAD2:FA after case 3b'
            CALL OUTPUT (FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
      END IF
C
C Case 4)
C
C
C Case 4b)
C         -ZYMAT(p,r)(xq|ry) D(x,y) -> F(p,q)
C         -ZYMAT(p,C)(xq|CD) D(x,D) -> F(p,q)             D active
C
      IF (IDTYP.EQ.JTACT) THEN
         IPSYM=MULD2H(IZYSYM,ICSYM)
         IQSYM=MULD2H(IFCKSY,IPSYM)
         IXSYM=MULD2H(IDENSY,IDSYM)
         IF (ICTYP.EQ.JTSEC) THEN
            NORBP=NOCC(IPSYM)
         ELSE
            NORBP=NORB(IPSYM)
         END IF
         NORBQ=NORB(IQSYM)
         NASHX=NASH(IXSYM)
         IF (NORBP.GT.0 .AND. NORBQ.GT.0 .AND. NASHX.GT.0) THEN
            IPST = IORB(IPSYM) + 1
            IQST = IORB(IQSYM) + 1
            IXST = IORB(IXSYM) + NISH(IXSYM) + 1
            NXW = IASH(IXSYM) + 1
            CALL DGEMM('T','N',1,NORBQ,NASHX,1.D0,
     &                 DEN1(NXW,NDW),NASHT,
     &                 H2X(IXST,IQST),NORBT,0.D0,
     &                 WRK,1)
            CALL DGEMM('N','N',NORBP,NORBQ,1,-1.D0,
     &                 ZYMAT(IPST,IC),NORBT,
     &                 WRK,1,1.D0,
     &                 FA(IPST,IQST),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') 'FAXAD2:FA after case 4b'
               CALL OUTPUT (FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
            END IF
         END IF
      END IF
C
C Case 4a)
C         (xq|pr)ZYMAT(r,y) D(x,y) -> F(p,q)
C         (xq|CD)ZYMAT(D,y) D(x,y) -> F(C,q)
C
      IYSYM=MULD2H(IZYSYM,IDSYM)
      IXSYM=MULD2H(IDENSY,IYSYM)
      IQSYM=MULD2H(IFCKSY,ICSYM)
      NASHX=NASH(IXSYM)
      NASHY=NASH(IYSYM)
      IF (ICTYP.EQ.JTSEC) THEN
         NORBQ=NOCC(IQSYM)
      ELSE
         NORBQ=NORB(IQSYM)
      END IF
      IF (NORBQ.GT.0 .AND. NASHX.GT.0 .AND. NASHY.GT.0) THEN
         IXST=IORB(IXSYM)+NISH(IXSYM)+1
         IYST=IORB(IYSYM)+NISH(IYSYM)+1
         NXW=IASH(IXSYM)+1
         NYW=IASH(IYSYM)+1
         IQST=IORB(IQSYM)+1
         CALL DGEMM('N','T',1,NASHX,NASHY,1.D0,
     &              ZYMAT(ID,IYST),NORBT,
     &              DEN1(NXW,NYW),NASHT,0.D0,
     &              WRK,1)
         CALL DGEMM('N','N',1,NORBQ,NASHX,1.D0,
     &              WRK,1,
     &              H2X(IXST,IQST),NORBT,1.D0,
     &              FA(IC,IQST),NORBT)
         IF (IPRRSP.GT.300) THEN
            WRITE(LUPRI,'(/A)') 'FAXAD2:FA after case 4a'
            CALL OUTPUT (FA,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         END IF
      END IF
      CALL QEXIT('FAXAD2')
      RETURN
      END
C  /* Deck qxad1 */
      SUBROUTINE QXAD1(INTSYM,IC,ID,ICDTYP,
     *                ICDSYM,ICSYM,IDSYM,QA,QB, H2X,
     *                DEN2A,DEN2B,IDENSY,PVDEN,H2XAC)
C
C
C Olav Vahtras
C Dec 9 1991
C
C     Process contributions from the Q matrices. Formulas:
C
C     QA(p,q) = QA(p,q)
C              + (xp|yz)d2a(xq;yz) + (xy|zp)d2b(xy;zq)
C
C     QB(p,q) = QB(p,q)
C              + (px|yz)d2a(qx;yz) + (xy|pz)d2b(xy;qz)
C
C     where d2a=<e(SA*S1,S2> and d2b = <e(S1,SA*S2)>
C     The sum over x and v is taken by reading in all active contribu-
C     tions.
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <infind.h>
#include <wrkrsp.h>
#include <infpri.h>
#include "rspprp.h"
#include <infhyp.h>
C
      DIMENSION QA(NORBT,NASHDI),QB(NORBT,NASHDI)
      DIMENSION DEN2A(NASHT,NASHT,NASHT,NASHT)
      DIMENSION DEN2B(NASHT,NASHT,NASHT,NASHT)
      DIMENSION PVDEN(NASHDI,NASHDI)
      DIMENSION H2X(NORBT,NORBT),  H2XAC(NASHT,NASHT)
C
      IF (ICDTYP.EQ.JTSESE .OR. NASHT.LE.1)  RETURN
C
      CALL QENTER('QXAD1')
      IFCKSY = MULD2H(INTSYM,IDENSY)
      ICTYP = IOBTYP(IC)
      IDTYP = IOBTYP(ID)
      IF (IPRRSP.GT.300) THEN
         WRITE(LUPRI,'(/A)') 'ENTER QXAD1'
         WRITE(LUPRI,'(A,I5)') 'INTSYM =',INTSYM
         WRITE(LUPRI,'(A,2I5)') 'IC ID =',IC,ID
         WRITE(LUPRI,'(A,I5)') 'ICDTYP =',ICDTYP
         WRITE(LUPRI,'(A,I5)') 'ICDSYM =',ICDSYM
         WRITE(LUPRI,'(A,I5)') 'ICSYM =',ICSYM
         WRITE(LUPRI,'(A,I5)') 'IDSYM =',IDSYM
         WRITE(LUPRI,'(A)') 'QA'
         CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'QB'
         CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'DEN2A'
         CALL OUTPUT(DEN2A,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
         WRITE(LUPRI,'(A)') 'DEN2B'
         CALL OUTPUT(DEN2B,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
         WRITE(LUPRI,'(A,I5)') 'IDENSY =',IDENSY
         WRITE(LUPRI,'(/A)') 'H2X'
         CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(/A)') 'H2XAC'
         CALL OUTPUT(H2XAC,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
C
C Case 1)
C
C        (xp|yz)d2a(xq;yz)  -> QA(p,q)
C        (xp|CD)d2a(xq;CD)  -> QA(p,q)
C
C
C
      ICW = ISW( IC )
      IDW = ISW( ID )
      NCW = ICW - NISHT
      NDW = IDW - NISHT
      IF (ICDTYP.EQ.JTACAC) THEN
C
         DO 2000 IPSYM = 1,NSYM
            IQSYM = MULD2H(IFCKSY,IPSYM)
            IABSYM = MULD2H(INTSYM,ICDSYM)
            IXSYM = MULD2H(IPSYM,IABSYM)
            NORBP = NORB(IPSYM)
            NASHQ = NASH(IQSYM)
            NASHX = NASH(IXSYM)
            IF (NORBP.GT.0 .AND. NASHQ.GT.0 .AND. NASHX.GT.0) THEN
               IPST = IORB(IPSYM) + 1
               NQW = IASH(IQSYM) + 1
               NXW = IASH(IXSYM) + 1
               IXST = ISX(NXW + NISHT)
C
C        (xp|yz)d2a(xq;yz)  -> QA(p,q)
C        (xp|CD)d2a(xq;CD)  -> QA(p,q)
C
               CALL DGEMM('T','N',NORBP,NASHQ,NASHX,1.D0,
     &                    H2X(IXST,IPST),NORBT,
     &                    DEN2A(NXW,NQW,NCW,NDW),NASHT,1.D0,
     &                    QA(IPST,NQW),NORBT)
               IF (IPRRSP.GT.300) THEN
                  WRITE(LUPRI,'(/A,I2)')
     *                     ' QXAD1: QA after Case 1: IPSYM =', IPSYM
                  CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
               END IF
C
C Case 2)
C            (px|yz)d2a(qx;yz) -> QB(p,q)
C            (px|CD)d2a(qx;CD) -> QB(p,q)
C
               CALL DGEMM('N','T',NORBP,NASHQ,NASHX,1.D0,
     &                    H2X(IPST,IXST),NORBT,
     &                    DEN2A(NQW,NXW,NCW,NDW),NASHT,1.D0,
     &                    QB(IPST,NQW),NORBT)
               IF (IPRRSP.GT.300) THEN
                  WRITE(LUPRI,'(/A,I2)')
     *                     ' QXAD1: QB after Case 2: IPSYM =', IPSYM
                  CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
               END IF
            END IF
2000     CONTINUE
      END IF
C
C
C Case 3)
C                (xy|zp)d2b(xy;zq) -> QA(p,q)
C                (xy|CD)d2b(xy;Cq) -> QA(D,q)
C
      IF (ICTYP.EQ.JTACT) THEN
         IQSYM=MULD2H(IFCKSY,IDSYM)
         NASHQ = NASH(IQSYM)
         IF (NASHQ.GT.0) THEN
            NQW = IASH(IQSYM) + 1
            CALL DGEMM('N','N',1,NASHQ,N2ASHX,1.D0,
     &                 H2XAC,1,
     &                 DEN2B(1,1,NCW,NQW),NASHT*N2ASHX,1.D0,
     &                 QA(ID,NQW),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)')
     *            ' QXAD1: QA after Case 3'
               CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            END IF
         END IF
      END IF
C
C Case 4)
C                (xy|pz)d2b(xy;qz) -> QB(p,q)
C                (xy|CD)d2b(xy;qD) -> QB(C,q)
C
      IF (IDTYP.EQ.JTACT) THEN
         IQSYM=MULD2H(IFCKSY,ICSYM)
         NASHQ=NASH(IQSYM)
         IF (NASHQ.GT.0) THEN
            NQW = IASH(IQSYM) + 1
            CALL DGEMM('N','N',1,NASHQ,N2ASHX,1.D0,
     &                 H2XAC,1,
     &                 DEN2B(1,1,NQW,NDW),N2ASHX,1.D0,
     &                 QB(IC,NQW),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)')
     *            ' QXAD1: QB after Case 4'
               CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            END IF
         END IF
      END IF
      CALL QEXIT('QXAD1')
      RETURN
      END
C  /* Deck qxad2 */
      SUBROUTINE QXAD2(INTSYM,IC,ID,ICDTYP,
     *                ICDSYM,ICSYM,IDSYM,QA,QB, H2X,DEN2A,DEN2B,
     *                IDENSY,PVDEN, H2XAC,SCR,ZYMAT,IZYSYM)
C
C
C Olav Vahtras
C Dec 9 1991
C
C     Process contributions from the Q matrices. Formulas:
C
C     QA(p,q) = QA(p,q)
C              + (xp|y~z)d2a(xq;yz) + (xy|z~p)d2b(xy;zq)
C
C     QB(p,q) = QB(p,q)
C              + (px|y~z)d2a(qx;yz) + (xy|p~z)d2b(xy;qz)
C
C Case 1)
C     (xp|y~z) = ZYMAT(y,r)(xp|rz) - (xp|yr)ZYMAT(r,z)
C
C     1a)  ZYMAT(y,r)(xp|rz)d2a(xq;yz) -> QA(p,q)
C          ZYMAT(y,C)(xp|CD)d2a(xq;yD) -> QA(p,q)    D active
C
C     1b) -(xp|yr)ZYMAT(r,z)d2a(xq;yz) -> QA(p,q)
C         -(xp|CD)ZYMAT(D,z)d2a(xq;Cz) -> QA(p,q)    C active
C
C
C Case 2)
C     (px|y~z) = ZYMAT(y,r)(px|rz) - (px|yr)ZYMAT(r,z)
C
C     2a)  ZYMAT(y,r)(px|rz)d2a(qx;yz) -> QB(p,q)
C          ZYMAT(y,C)(px|CD)d2a(qx;yD) -> QB(p,q)    D active
C
C     2b) -(px|yr)ZYMAT(r,z)d2a(qx;yz) -> QB(p,q)
C         -(px|CD)ZYMAT(D,z)d2a(qx;Cz) -> QB(p,q)    C active
C
C Case 3)
C    (xy|z~p) = ZYMAT(z,r)(xy|rp) - (xy|zr)ZYMAT(r,p)
C
C     3a)  ZYMAT(z,r)(xy|rp)d2b(xy;zq) -> QA(p,q)
C          ZYMAT(z,C)(xy|CD)d2b(xy;zq) -> QA(D,q)
C
C     3b) -(xy|zr)ZYMAT(r,p)d2b(xy;zq) -> QA(p,q)
C         -(xy|CD)ZYMAT(D,p)d2b(xy;Cq) -> QA(p,q)      C active
C
C Case 4)
C   (xy|p~z) = ZYMAT(p,r)(xy|rz) - (xy|pr) ZYMAT(r,z)
C
C     4a)  ZYMAT(p,r)(xy|rz)d2b(xy;qz) -> QB(p,q)
C          ZYMAT(p,C)(xy|CD)d2b(xy;qD) -> QB(p,q)      D active
C
C     4b) -(xy|pr)ZYMAT(r,z)d2b(xy;qz) -> QB(p,q)
C         -(xy|CD)ZYMAT(D,z)d2b(xy;qz) -> QB(C,q)
C
C     where d2a=<e(SA*S1,S2> and d2 = <e(S1,SA*S2)>
C     The sum over x and v is taken by reading in all active contribu-
C     tions.
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <infrsp.h>
#include <inforb.h>
#include <infdim.h>
#include <infind.h>
#include <wrkrsp.h>
#include <infpri.h>
#include "rspprp.h"
#include <infhyp.h>
C
      DIMENSION QA(NORBT,NASHDI),QB(NORBT,NASHDI)
      DIMENSION DEN2A(NASHT,NASHT,NASHT,NASHT)
      DIMENSION DEN2B(NASHT,NASHT,NASHT,NASHT)
      DIMENSION PVDEN(NASHDI,NASHDI)
      DIMENSION H2X(NORBT,NORBT), H2XAC(NASHT,NASHT)
C Allocated for temporary space max(n2ashx,norbt)
      DIMENSION SCR(*)
      DIMENSION ZYMAT(NORBT,NORBT)
C
      IF (NASHT.LE.1) RETURN
C
      CALL QENTER('QXAD2')
      ICTYP = IOBTYP(IC)
      IDTYP = IOBTYP(ID)
      IF (IPRRSP.GT.300) THEN
         WRITE(LUPRI,'(/A)') 'ENTER QXAD2'
         WRITE(LUPRI,'(A,I5)') 'INTSYM =',INTSYM
         WRITE(LUPRI,'(A,I5)') 'IZYSYM =',IZYSYM
         WRITE(LUPRI,'(A,2I5)') 'IC ID =',IC,ID
         WRITE(LUPRI,'(A,I5)') 'ICDTYP =',ICDTYP
         WRITE(LUPRI,'(A,I5)') 'ICDSYM =',ICDSYM
         WRITE(LUPRI,'(A,I5)') 'ICSYM =',ICSYM
         WRITE(LUPRI,'(A,I5)') 'IDSYM =',IDSYM
         WRITE(LUPRI,'(A,I5)') 'IDENSY =',IDENSY
         WRITE(LUPRI,'(/A)') 'QA'
         CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A)') 'QB'
         CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(A)') 'DEN2A'
         CALL OUTPUT(DEN2A,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
         WRITE(LUPRI,'(A)') 'DEN2B'
         CALL OUTPUT(DEN2B,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
         WRITE(LUPRI,'(/A)') 'ZYMAT'
         CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(/A)') 'H2X'
         CALL OUTPUT(H2X,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(/A)') 'H2XAC'
         CALL OUTPUT(H2XAC,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
      END IF
      ICW = ISW( IC )
      IDW = ISW( ID )
      NCW = ICW - NISHT
      NDW = IDW - NISHT
      IFCKSY = MULD2H(INTSYM,MULD2H(IDENSY,IZYSYM))
C
C
C     1a)  ZYMAT(y,r)(xp|rz)d2a(xq;yz) -> QA(p,q)
C          ZYMAT(y,C)(xp|CD)d2a(xq;yD) -> QA(p,q)    D active
C
C     2a)  ZYMAT(y,r)(px|rz)d2a(qx;yz) -> QB(p,q)
C          ZYMAT(y,C)(px|CD)d2a(qx;yD) -> QB(p,q)    D active
C
C
      IF (IDTYP.EQ.JTACT) THEN
         DO 1000 IPSYM = 1,NSYM
            IQSYM = MULD2H(IFCKSY,IPSYM)
            IXSYM = MULD2H(IPSYM,MULD2H(INTSYM,ICDSYM))
            IYSYM = MULD2H(IZYSYM,ICSYM)
            NORBP = NORB(IPSYM)
            NASHQ = NASH(IQSYM)
            NASHX = NASH(IXSYM)
            NASHY = NASH(IYSYM)
            IF (NORBP.GT.0 .AND. NASHQ.GT.0 .AND.
     *         NASHX.GT.0 .AND. NASHY.GT.0) THEN
               IPST = IORB(IPSYM) + 1
               NQW = IASH(IQSYM) + 1
               NXW = IASH(IXSYM) + 1
               IXST = ISX(NXW+NISHT)
               DO 1500 IYW = 1,NASHY
                  NYW = IASH(IYSYM) + IYW
                  IY = ISX(NYW + NISHT)
                  ZYYC = ZYMAT(IY,IC)
C
C     1a)  ZYMAT(y,r)(xp|rz)d2a(xq;yz) -> QA(p,q)
C          ZYMAT(y,C)(xp|CD)d2a(xq;yD) -> QA(p,q)    D active
C
                  CALL DGEMM('T','N',NORBP,NASHQ,NASHX,ZYYC,
     &                       H2X(IXST,IPST),NORBT,
     &                       DEN2A(NXW,NQW,NYW,NDW),NASHT,1.D0,
     &                       QA(IPST,NQW),NORBT)
                  IF (IPRRSP.GT.300) THEN
                     WRITE(LUPRI,'(/A,I2,A,I2,A,I2)')
     *                  ' QXAD2: QA after Case 1a, IYW =', IYW,
     *                  '(',NASHY,') ;IPSYM = ',IPSYM
                     CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
                  END IF
C
C     2a)  ZYMAT(y,r)(px|rz)d2a(qx;yz) -> QB(p,q)
C          ZYMAT(y,C)(px|CD)d2a(qx;yD) -> QB(p,q)    D active
C
                  CALL DGEMM('N','T',NORBP,NASHQ,NASHX,ZYYC,
     &                       H2X(IPST,IXST),NORBT,
     &                       DEN2A(NQW,NXW,NYW,NDW),NASHT,1.D0,
     &                       QB(IPST,NQW),NORBT)
                  IF (IPRRSP.GT.300) THEN
                     WRITE(LUPRI,'(/A,I2,A,I2,A,I2)')
     *                  ' QXAD2: QB after Case 2a, IYW =',IYW,
     *                  '(',NASHY,'); IPSYM = ',IPSYM
                     CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
                  END IF
1500           CONTINUE
            END IF
1000     CONTINUE
      END IF


C
C     1b) -(xp|yr)ZYMAT(r,z)d2a(xq;yz) -> QA(p,q)
C         -(xp|CD)ZYMAT(D,z)d2a(xq;Cz) -> QA(p,q)    C active
C
C     2b) -(px|yr)ZYMAT(r,z)d2a(qx;yz) -> QB(p,q)
C         -(px|CD)ZYMAT(D,z)d2a(qx;Cz) -> QB(p,q)    C active
C
      IF (ICTYP.EQ.JTACT) THEN
         DO 2000 IPSYM=1,NSYM
            IQSYM = MULD2H(IFCKSY,IPSYM)
            IZSYM = MULD2H(IZYSYM,IDSYM)
            IXSYM = MULD2H(IPSYM,MULD2H(INTSYM,ICDSYM))
            NORBP = NORB(IPSYM)
            NASHQ = NASH(IQSYM)
            NASHX = NASH(IXSYM)
            NASHZ = NASH(IZSYM)
            IF (NORBP.GT.0 .AND. NASHQ.GT.0 .AND.
     *         NASHX.GT.0 .AND. NASHZ.GT.0) THEN
               IPST = IORB(IPSYM) + 1
               NQW = IASH(IQSYM) + 1
               NXW = IASH(IXSYM) + 1
               IXST = ISX(NXW+NISHT)
               DO 2500 IZW = 1,NASHZ
                  NZW = IZW + IASH(IZSYM)
                  IZ = ISX(NZW + NISHT)
                  ZYDZ = ZYMAT(ID,IZ)
C
C     1b) -(xp|yr)ZYMAT(r,z)d2a(xq;yz) -> QA(p,q)
C         -(xp|CD)ZYMAT(D,z)d2a(xq;Cz) -> QA(p,q)    C active
C
                  CALL DGEMM('T','N',NORBP,NASHQ,NASHX,-ZYDZ,
     &                       H2X(IXST,IPST),NORBT,
     &                       DEN2A(NXW,NQW,NCW,NZW),NASHT,1.D0,
     &                       QA(IPST,NQW),NORBT)
                  IF (IPRRSP.GT.300) THEN
                     WRITE(LUPRI,'(/A,I2)')
     *                  ' QXAD2: QA after Case 1b: IZW =',IZW
                     CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
                  END IF
C
C     2b) -(px|yr)ZYMAT(r,z)d2a(qx;yz) -> QB(p,q)
C         -(px|CD)ZYMAT(D,z)d2a(qx;Cz) -> QB(p,q)    C active
C
                  CALL DGEMM('N','T',NORBP,NASHQ,NASHX,-ZYDZ,
     &                       H2X(IPST,IXST),NORBT,
     &                       DEN2A(NQW,NXW,NCW,NZW),NASHT,1.D0,
     &                       QB(IPST,NQW),NORBT)
                  IF (IPRRSP.GT.300) THEN
                     WRITE(LUPRI,'(/A,I2)')
     *                  ' QXAD2: QB after Case 2b: IZW =',IZW
                     CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
                  END IF
2500           CONTINUE
            END IF
2000     CONTINUE
      END IF

C
C
C     3a)  ZYMAT(z,r)(xy|rp)d2b(xy;zq) -> QA(p,q)
C          ZYMAT(z,C)(xy|CD)d2b(xy;zq) -> QA(D,q)
C
      IZSYM = MULD2H(IZYSYM,ICSYM)
      IQSYM = MULD2H(IFCKSY,IDSYM)
      NASHZ = NASH(IZSYM)
      NASHQ = NASH(IQSYM)
      IF (NASHZ.GT.0 .AND. NASHQ.GT.0) THEN
         NZW = IASH(IZSYM) + 1
         IZ = ISX(NZW + NISHT)
         DO 3000 IQW = 1,NASHQ
            NQW = IQW + IASH(IQSYM)
            CALL DGEMM('N','N',1,NASHZ,N2ASHX,1.D0,
     &                 H2XAC,1,
     &                 DEN2B(1,1,NZW,NQW),N2ASHX,0.D0,
     &                 SCR,1)
            CALL DGEMM('N','N',1,1,NASHZ,1.D0,
     &                 SCR,1,
     &                 ZYMAT(IZ,IC),NORBT,1.D0,
     &                 QA(ID,NQW),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A,I2,A,I1,A)')
     *            ' QXAD2: QA after Case 3a: IQW =',IQW,' (',NASHQ,')'
               CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            END IF
3000     CONTINUE
      END IF
C
C     4b) -(xy|pr)ZYMAT(r,z)d2b(xy;qz) -> QB(p,q)
C         -(xy|CD)ZYMAT(D,z)d2b(xy;qz) -> QB(C,q)
C
      IZSYM = MULD2H(IZYSYM,IDSYM)
      IQSYM = MULD2H(IFCKSY,ICSYM)
      NASHZ = NASH(IZSYM)
      NASHQ = NASH(IQSYM)
      IF (NASHZ.GT.0 .AND. NASHQ.GT.0) THEN
         NQW = IASH(IQSYM) + 1
         DO 4000 IZW = 1,NASHZ
            NZW = IASH(IZSYM) + IZW
            IZ = ISX(NZW + NISHT)
            ZYDZ = ZYMAT(ID,IZ)
            CALL DGEMM('N','N',1,NASHQ,N2ASHX,1.D0,
     &                 H2XAC,1,
     &                 DEN2B(1,1,NQW,NZW),N2ASHX,0.D0,
     &                 SCR,1)
            CALL DAXPY(NASHQ,-ZYDZ,SCR,1,QB(IC,NQW),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A,I2,A,I1,A)')
     *            ' QXAD2: QB after Case 4b: IZW =',IZW,'(',NASHZ,')'
               CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            END IF
4000     CONTINUE
      END IF
C
C     3b) -(xy|zr)ZYMAT(r,p)d2b(xy;zq) -> QA(p,q)
C         -(xy|CD)ZYMAT(D,p)d2b(xy;Cq) -> QA(p,q)      C active
C
      IF (ICTYP.EQ.JTACT) THEN
         IPSYM = MULD2H(IZYSYM,IDSYM)
         IQSYM = MULD2H(IFCKSY,IPSYM)
         IF (IDTYP.EQ.JTSEC) THEN
            NORBP = NOCC(IPSYM)
         ELSE
            NORBP = NORB(IPSYM)
         END IF
         NASHQ = NASH(IQSYM)
         IF (NORBP.GT.0 .AND. NASHQ.GT.0) THEN
            IPST = IORB(IPSYM) + 1
            NQW = IASH(IQSYM) + 1
            CALL DGEMM('N','N',1,NASHQ,N2ASHX,1.D0,
     &                 H2XAC,1,
     &                 DEN2B(1,1,NCW,NQW),N2ASHX*NASHT,0.D0,
     &                 SCR,1)
            CALL DGEMM('T','N',NORBP,NASHQ,1,-1.D0,
     &                 ZYMAT(ID,IPST),NORBT,
     &                 SCR,1,1.D0,
     &                 QA(IPST,NQW),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') ' QXAD2: QA after Case 3b'
               CALL OUTPUT(QA,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            END IF
         END IF
      END IF
C
C     4a)  ZYMAT(p,r)(xy|rz)d2b(xy;qz) -> QB(p,q)
C          ZYMAT(p,C)(xy|CD)d2b(xy;qD) -> QB(p,q)      D active
C
      IF (IDTYP.EQ.JTACT) THEN
         IPSYM = MULD2H(IZYSYM,ICSYM)
         IQSYM = MULD2H(IFCKSY,IPSYM)
         IF (ICTYP.EQ.JTSEC) THEN
            NORBP = NOCC(IPSYM)
         ELSE
            NORBP = NORB(IPSYM)
         END IF
         NASHQ = NASH(IQSYM)
         IF (NORBP.GT.0 .AND. NASHQ.GT.0) THEN
            IPST = IORB(IPSYM) + 1
            NQW = IASH(IQSYM) + 1
            CALL DGEMM('N','N',1,NASHQ,N2ASHX,1.D0,
     &                 H2XAC,1,
     &                 DEN2B(1,1,NQW,NDW),N2ASHX,0.D0,
     &                 SCR,1)
            CALL DGEMM('N','N',NORBP,NASHQ,1,1.D0,
     &                 ZYMAT(IPST,IC),NORBT,
     &                 SCR,1,1.D0,
     &                 QB(IPST,NQW),NORBT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') ' QXAD2: QB after Case 4a'
               CALL OUTPUT(QB,1,NORBT,1,NASHT,NORBT,NASHT,1,LUPRI)
            END IF
         END IF
      END IF
      CALL QEXIT('QXAD2')
      RETURN
      END
C  /* Deck extrac */
      SUBROUTINE EXTRAC(H,HAC)
C
C Extract the active part of H into HAC
C
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
      DIMENSION H(NORBT,NORBT), HAC(NASHT,NASHT)
      DO 110 JW = 1, NASHT
         JX = ISX(NISHT + JW)
         DO 111 IW = 1, NASHT
            IX = ISX(NISHT + IW)
            HAC(IW,JW) = H(IX,JX)
111      CONTINUE
110   CONTINUE
      RETURN
      END
C  /* Deck disac1 */
      SUBROUTINE DISAC1(IC,ID,H2XAC,H2ACAC)
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infind.h>
      DIMENSION H2XAC(N2ASHX)
      DIMENSION H2ACAC(N2ASHX,NASHT,NASHT)
      PARAMETER (D1 = 1.0D0)
      IF (IOBTYP(IC).NE.JTACT .OR. IOBTYP(ID).NE.JTACT) RETURN
      NCW = ISW(IC) - NISHT
      NDW = ISW(ID) - NISHT
      CALL DAXPY(N2ASHX,D1,H2XAC,1,H2ACAC(1,NCW,NDW),1)
      RETURN
      END
C  /* Deck disac2 */
      SUBROUTINE DISAC2(IC,ID,H2XAC,H2ACAC,ZYMAT,IZYSYM)
#include <implicit.h>
#include <maxorb.h>
#include <maxash.h>
#include <priunit.h>
#include <inforb.h>
#include <infrsp.h>
#include <infind.h>
#include "rspprp.h"
#include <infhyp.h>
#include <infpri.h>
      DIMENSION H2XAC(N2ASHX)
      DIMENSION H2ACAC(N2ASHX,NASHT,NASHT)
      DIMENSION ZYMAT(NORBT,NORBT)
      ICTYP = IOBTYP(IC)
      IDTYP = IOBTYP(ID)
      ICSYM = ISMO(IC)
      IDSYM = ISMO(ID)
      IF (ICTYP.NE.JTACT .AND. IDTYP.NE.JTACT) RETURN
      IF (IPRRSP.GT.300) THEN
         WRITE(LUPRI,'(/A)') 'ENTER DISAC2'
         WRITE(LUPRI,'(A,2I3)') 'IC ID',IC,ID
         WRITE(LUPRI,'(/A)') ' H2XAC'
         CALL OUTPUT(H2XAC,1,NASHT,1,NASHT,NASHT,NASHT,1,LUPRI)
         WRITE(LUPRI,'(/A)') ' H2ACAC'
         CALL OUTPUT(H2ACAC,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,LUPRI)
         WRITE(LUPRI,'(/A)') 'ZYMAT'
         CALL OUTPUT(ZYMAT,1,NORBT,1,NORBT,NORBT,NORBT,1,LUPRI)
         WRITE(LUPRI,'(/A,I3)') 'IZYSYM',IZYSYM
      END IF
      NCW = ISW(IC) - NISHT
      NDW = ISW(ID) - NISHT
C
C (xy|z~w) = ZYMAT(z,r)(xy|rw) - (xy|zr)ZYMAT(r,w)
C
C Case 1)    ZYMAT(z,C)(xy|CD) -> (xy|z~D)        D active
C
C Case 2)  - (xy|CD)ZYMAT(D,w) -> (xy|C~w)        C active
C
C
C
C Case 1)    ZYMAT(z,C)(xy|CD) -> (xy|z~D)        D active
C
      IF (IDTYP.EQ.JTACT) THEN
         IZSYM = MULD2H(IZYSYM,ICSYM)
         NASHZ = NASH(IZSYM)
         IF (NASHZ.GT.0) THEN
            NZW = IASH(IZSYM) + 1
            IZ = ISX(NZW + NISHT)
            CALL DGEMM('N','N',N2ASHX,NASHZ,1,1.D0,
     &                 H2XAC,N2ASHX,
     &                 ZYMAT(IZ,IC),1,1.D0,
     &                 H2ACAC(1,NZW,NDW),N2ASHX)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') ' H2ACAC after case 1'
              CALL OUTPUT(H2ACAC,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,
     &                    LUPRI)
            END IF
         END IF
      END IF
C
C Case 2)  - (xy|CD)ZYMAT(D,w) -> (xy|C~w)        C active
C
      IF (ICTYP.EQ.JTACT) THEN
         IWSYM = MULD2H(IDSYM,IZYSYM)
         NASHW = NASH(IWSYM)
         IF (NASHW.GT.0) THEN
            NWW = IASH(IWSYM) + 1
            IW = ISX(NWW + NISHT)
            CALL DGEMM('N','N',N2ASHX,NASHW,1,-1.D0,
     &                 H2XAC,N2ASHX,
     &                 ZYMAT(ID,IW),NORBT,1.D0,
     &                 H2ACAC(1,NCW,NWW),N2ASHX*NASHT)
            IF (IPRRSP.GT.300) THEN
               WRITE(LUPRI,'(/A)') ' H2ACAC after case 2'
              CALL OUTPUT(H2ACAC,1,N2ASHX,1,N2ASHX,N2ASHX,N2ASHX,1,
     &                    LUPRI)
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck priac2 */
      SUBROUTINE PRIAC2(AC2,N,LU)
#include <implicit.h>
      DIMENSION AC2(N*N,N,N)
      DO 10 I=1,N
         DO 20 J=1,I
            WRITE(LU,'(/A,2I3,A)') ' Distribution (* *,',I,J,')'
            CALL OUTPUT(AC2(1,I,J),1,N,1,N,N,N,1,LU)
            IF (I.NE.J) THEN
               WRITE(LU,'(/A,2I3,A)') ' Distribution (* *,',J,I,')'
               CALL OUTPUT(AC2(1,J,I),1,N,1,N,N,N,1,LU)
            END IF
20       CONTINUE
10    CONTINUE
      RETURN
      END
