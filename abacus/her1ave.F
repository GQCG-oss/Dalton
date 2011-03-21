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
      SUBROUTINE AVENA1(ADER,NATOMC,SECDER,NCENTC,JCENTC,MAXCMP,JSYMC,
     &                  SIGNC,DSHELL,HESSNA)
C
C     tuh 1984
C     modified for symmetry Jul 5 1988 tuh
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.00 D00, D1 = 1.00 D00, D2 = 2.00D00)
      LOGICAL   SECDER, DOAX, DOAY, DOAZ, DOCX, DOCY, DOCZ
      DIMENSION ADER(KCKTAB,NATOMC,*), DSHELL(KHKTAB,*), NCENTC(*),
     &          JCENTC(*), SIGNC(3,*), JSYMC(*), HESSNA(MXCOOR,MXCOOR)
C
#include "onecom.h"
#include "ader.h"
#include "symmet.h"
#include "dorps.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
#include "ibtfun.h"
C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
C
C     Run over atoms in the operator
C
      DO 200 IATOM = 1,NATOMC
         ICENTC = NCENTC(IATOM)
         IF (ICENTC .NE. ICENTA) THEN
C
C           ********************************************
C           ***** Multiply densities and integrals *****
C           ********************************************
C
            DA0000 = D0
            DA000X = D0
            DA000Y = D0
            DA000Z = D0
            DA00XX = D0
            DA00XY = D0
            DA00XZ = D0
            DA00YY = D0
            DA00YZ = D0
            DA00ZZ = D0
            DO 300 ICOMP = 1, KHKTAB
               DENS = DSHELL(ICOMP,1)
               DA0000 = DA0000 + DENS*ADER(ICOMP,IATOM,IA0000)
               DA000X = DA000X + DENS*ADER(ICOMP,IATOM,IA000X)
               DA000Y = DA000Y + DENS*ADER(ICOMP,IATOM,IA000Y)
               DA000Z = DA000Z + DENS*ADER(ICOMP,IATOM,IA000Z)
               IF (SECDER) THEN
                  DA00XX = DA00XX + DENS*ADER(ICOMP,IATOM,IA00XX)
                  DA00XY = DA00XY + DENS*ADER(ICOMP,IATOM,IA00XY)
                  DA00XZ = DA00XZ + DENS*ADER(ICOMP,IATOM,IA00XZ)
                  DA00YY = DA00YY + DENS*ADER(ICOMP,IATOM,IA00YY)
                  DA00YZ = DA00YZ + DENS*ADER(ICOMP,IATOM,IA00YZ)
                  DA00ZZ = DA00ZZ + DENS*ADER(ICOMP,IATOM,IA00ZZ)
               END IF
  300       CONTINUE
C
            KCENTC = JCENTC(IATOM)
            NCX    = 3*KCENTC - 2
            NCY    = 3*KCENTC - 1
            NCZ    = 3*KCENTC
            SCX    = SIGNC(1,IATOM)
            SCY    = SIGNC(2,IATOM)
            SCZ    = SIGNC(3,IATOM)
C
C           ***** Energy *****
C
            ENERNA = ENERNA + DA0000
C
C           ********************
C           ***** Gradient *****
C           ********************
C
            IF (DOREPS(0)) THEN
               IAX  = IPTCNT(NAX,0,1)
               IAY  = IPTCNT(NAY,0,1)
               IAZ  = IPTCNT(NAZ,0,1)
               ICX  = IPTCNT(NCX,0,1)
               ICY  = IPTCNT(NCY,0,1)
               ICZ  = IPTCNT(NCZ,0,1)
               DOAX = IAX .NE. 0
               DOAY = IAY .NE. 0
               DOAZ = IAZ .NE. 0
               DOCX = ICX .NE. 0
               DOCY = ICY .NE. 0
               DOCZ = ICZ .NE. 0
C
C              A nuclear-attraction gradient elements:
C
               IF (DOAX) GRADNA(IAX) = GRADNA(IAX) - DA000X
               IF (DOAY) GRADNA(IAY) = GRADNA(IAY) - DA000Y
               IF (DOAZ) GRADNA(IAZ) = GRADNA(IAZ) - DA000Z
C
C              C nuclear-attraction gradient elements:
C
               IF (DOCX) GRADNA(ICX) = GRADNA(ICX) + SCX*DA000X
               IF (DOCY) GRADNA(ICY) = GRADNA(ICY) + SCY*DA000Y
               IF (DOCZ) GRADNA(ICZ) = GRADNA(ICZ) + SCZ*DA000Z
            END IF
C
C           *******************
C           ***** Hessian *****
C           *******************
C
            IF (SECDER) THEN
               ISYMPC = JSYMC(IATOM)
               FAC = D1
               IF (NCENTA .EQ. KCENTC) FAC = D2
C
C              Run over irreps
C
               DO 400 IREP = 0, MAXREP
               IF (DOREPS(IREP)) THEN
                  CHIC = PT(IBTAND(ISYMPC,IREP))
                  CSCX = CHIC*SCX
                  CSCY = CHIC*SCY
                  CSCZ = CHIC*SCZ
                  IAX  = IPTCNT(NAX,IREP,1)
                  IAY  = IPTCNT(NAY,IREP,1)
                  IAZ  = IPTCNT(NAZ,IREP,1)
                  ICX  = IPTCNT(NCX,IREP,1)
                  ICY  = IPTCNT(NCY,IREP,1)
                  ICZ  = IPTCNT(NCZ,IREP,1)
                  DOAX = IAX .NE. 0
                  DOAY = IAY .NE. 0
                  DOAZ = IAZ .NE. 0
                  DOCX = ICX .NE. 0
                  DOCY = ICY .NE. 0
                  DOCZ = ICZ .NE. 0
C
C                 A-A nuclear-attraction Hessian elements:
C
                  IF (DOAX)          HESSNA(IAX,IAX) =
     *                               HESSNA(IAX,IAX) + DA00XX
                  IF (DOAX.AND.DOAY) HESSNA(IAX,IAY) =
     *                               HESSNA(IAX,IAY) + DA00XY
                  IF (DOAX.AND.DOAZ) HESSNA(IAX,IAZ) =
     *                               HESSNA(IAX,IAZ) + DA00XZ
                  IF (DOAY)          HESSNA(IAY,IAY) =
     *                               HESSNA(IAY,IAY) + DA00YY
                  IF (DOAY.AND.DOAZ) HESSNA(IAY,IAZ) =
     *                               HESSNA(IAY,IAZ) + DA00YZ
                  IF (DOAZ)          HESSNA(IAZ,IAZ) =
     *                               HESSNA(IAZ,IAZ) + DA00ZZ
C
C                 A-C nuclear-attraction Hessian elements:
C
                  IF (DOAX.AND.DOCX) HESSNA(IAX,ICX) =
     *                               HESSNA(IAX,ICX) - FAC*CSCX*DA00XX
                  IF (DOAX.AND.DOCY) HESSNA(IAX,ICY) =
     *                               HESSNA(IAX,ICY) - CSCY*DA00XY
                  IF (DOAX.AND.DOCZ) HESSNA(IAX,ICZ) =
     *                               HESSNA(IAX,ICZ) - CSCZ*DA00XZ
                  IF (DOAY.AND.DOCX) HESSNA(IAY,ICX) =
     *                               HESSNA(IAY,ICX) - CSCX*DA00XY
                  IF (DOAY.AND.DOCY) HESSNA(IAY,ICY) =
     *                               HESSNA(IAY,ICY) - FAC*CSCY*DA00YY
                  IF (DOAY.AND.DOCZ) HESSNA(IAY,ICZ) =
     *                               HESSNA(IAY,ICZ) - CSCZ*DA00YZ
                  IF (DOAZ.AND.DOCX) HESSNA(IAZ,ICX) =
     *                               HESSNA(IAZ,ICX) - CSCX*DA00XZ
                  IF (DOAZ.AND.DOCY) HESSNA(IAZ,ICY) =
     *                               HESSNA(IAZ,ICY) - CSCY*DA00YZ
                  IF (DOAZ.AND.DOCZ) HESSNA(IAZ,ICZ) =
     *                               HESSNA(IAZ,ICZ) - FAC*CSCZ*DA00ZZ
C
C                 C-C nuclear-attraction Hessian elements:
C
                  IF (DOCX)          HESSNA(ICX,ICX) =
     *                               HESSNA(ICX,ICX) + DA00XX
                  IF (DOCX.AND.DOCY) HESSNA(ICX,ICY) =
     *                               HESSNA(ICX,ICY) + SCX*SCY*DA00XY
                  IF (DOCX.AND.DOCZ) HESSNA(ICX,ICZ) =
     *                               HESSNA(ICX,ICZ) + SCX*SCZ*DA00XZ
                  IF (DOCY)          HESSNA(ICY,ICY) =
     *                               HESSNA(ICY,ICY) + DA00YY
                  IF (DOCY.AND.DOCZ) HESSNA(ICY,ICZ) =
     *                               HESSNA(ICY,ICZ) + SCY*SCZ*DA00YZ
                  IF (DOCZ)          HESSNA(ICZ,ICZ) =
     *                               HESSNA(ICZ,ICZ) + DA00ZZ
               END IF
  400          CONTINUE
            END IF
         END IF
  200 CONTINUE
      RETURN
      END
C  /* Deck avekfs */
      SUBROUTINE AVEKFS(STDER0,STDER1,STDER2,ISYMOP,MAXCMP,SECDER,
     &                  DSHELL,FSHELL,HESSKE,HESFS2)
C
C     tuh 1984
C
C     modified for symmetry tuh Jul 5 1988
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.00 D00)
      LOGICAL SECDER, DOAX, DOAY, DOAZ, DOBX, DOBY, DOBZ
      DIMENSION DSHELL(KHKTAB), FSHELL(KHKTAB), STDER0(KCKTAB,2),
     &          STDER1(KCKTAB,3,2), STDER2(KCKTAB,6,2),
     &          HESSKE(MXCOOR,MXCOOR), HESFS2(MXCOOR,MXCOOR)
#include <onecom.h>
#include <symmet.h>
#include <dorps.h>
#include <energy.h>
#include <ibtfun.h>
C
C     Multiply densities and Fock elements with integrals
C
      DERT0  = D0
      DERTX  = D0
      DERTY  = D0
      DERTZ  = D0
      DERSX  = D0
      DERSY  = D0
      DERSZ  = D0
      IF (SECDER) THEN
         DERTXX = D0
         DERTXY = D0
         DERTXZ = D0
         DERTYY = D0
         DERTYZ = D0
         DERTZZ = D0
         DERSXX = D0
         DERSXY = D0
         DERSXZ = D0
         DERSYY = D0
         DERSYZ = D0
         DERSZZ = D0
      END IF
      DO 200 ICOMP = 1, KHKTAB
         DENS   = DSHELL(ICOMP)
         FOCK   = FSHELL(ICOMP)
         DERSX  = DERSX  - FOCK*STDER1(ICOMP,1,1)
         DERSY  = DERSY  - FOCK*STDER1(ICOMP,2,1)
         DERSZ  = DERSZ  - FOCK*STDER1(ICOMP,3,1)
         DERT0  = DERT0  + DENS*STDER0(ICOMP,2)
         DERTX  = DERTX  + DENS*STDER1(ICOMP,1,2)
         DERTY  = DERTY  + DENS*STDER1(ICOMP,2,2)
         DERTZ  = DERTZ  + DENS*STDER1(ICOMP,3,2)
         IF (SECDER) THEN
            DERSXX = DERSXX - FOCK*STDER2(ICOMP,1,1)
            DERSXY = DERSXY - FOCK*STDER2(ICOMP,2,1)
            DERSXZ = DERSXZ - FOCK*STDER2(ICOMP,3,1)
            DERSYY = DERSYY - FOCK*STDER2(ICOMP,4,1)
            DERSYZ = DERSYZ - FOCK*STDER2(ICOMP,5,1)
            DERSZZ = DERSZZ - FOCK*STDER2(ICOMP,6,1)
            DERTXX = DERTXX + DENS*STDER2(ICOMP,1,2)
            DERTXY = DERTXY + DENS*STDER2(ICOMP,2,2)
            DERTXZ = DERTXZ + DENS*STDER2(ICOMP,3,2)
            DERTYY = DERTYY + DENS*STDER2(ICOMP,4,2)
            DERTYZ = DERTYZ + DENS*STDER2(ICOMP,5,2)
            DERTZZ = DERTZZ + DENS*STDER2(ICOMP,6,2)
         END IF
  200 CONTINUE
C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
      NBX = 3*NCENTB - 2
      NBY = 3*NCENTB - 1
      NBZ = 3*NCENTB
C
C     Undifferentiated kinetic energy
C     Note: One-center terms are not included!
C
      ENERKE = ENERKE + DERT0
C
C     ********************
C     ***** Gradient *****
C     ********************
C
      IF (DOREPS(0)) THEN
         IAX  = IPTCNT(NAX,0,1)
         IAY  = IPTCNT(NAY,0,1)
         IAZ  = IPTCNT(NAZ,0,1)
         IBX  = IPTCNT(NBX,0,1)
         IBY  = IPTCNT(NBY,0,1)
         IBZ  = IPTCNT(NBZ,0,1)
         DOAX = IAX .NE. 0
         DOAY = IAY .NE. 0
         DOAZ = IAZ .NE. 0
         DOBX = IBX .NE. 0
         DOBY = IBY .NE. 0
         DOBZ = IBZ .NE. 0
C
C        A gradient elements:
C
         IF (DOAX) THEN
            GRADKE(IAX)   = GRADKE(IAX)   + DERTX
            GRADFS(IAX)   = GRADFS(IAX)   + DERSX
         END IF
         IF (DOAY) THEN
            GRADKE(IAY)   = GRADKE(IAY)   + DERTY
            GRADFS(IAY)   = GRADFS(IAY)   + DERSY
         END IF
         IF (DOAZ) THEN
            GRADKE(IAZ)   = GRADKE(IAZ)   + DERTZ
            GRADFS(IAZ)   = GRADFS(IAZ)   + DERSZ
         END IF
C
C        B gradient elements:
C
         IF (DOBX) THEN
            GRADKE(IBX)   = GRADKE(IBX)   - SIGNBX*DERTX
            GRADFS(IBX)   = GRADFS(IBX)   - SIGNBX*DERSX
         END IF
         IF (DOBY) THEN
            GRADKE(IBY)   = GRADKE(IBY)   - SIGNBY*DERTY
            GRADFS(IBY)   = GRADFS(IBY)   - SIGNBY*DERSY
         END IF
         IF (DOBZ) THEN
            GRADKE(IBZ)   = GRADKE(IBZ)   - SIGNBZ*DERTZ
            GRADFS(IBZ)   = GRADFS(IBZ)   - SIGNBZ*DERSZ
         END IF
      END IF
C
C     *******************
C     ***** Hessian *****
C     *******************
C
      IF (SECDER) THEN
         DO 300 IREP = 0, MAXREP
         IF (DOREPS(IREP)) THEN
            CHI = PT(IBTAND(ISYMOP,IREP))
            CSBX = CHI*SIGNBX
            CSBY = CHI*SIGNBY
            CSBZ = CHI*SIGNBZ
            IF (NCENTA .EQ. NCENTB) THEN
               CASBX = CSBX + CSBX
               CASBY = CSBY + CSBY
               CASBZ = CSBZ + CSBZ
            ELSE
               CASBX = CSBX
               CASBY = CSBY
               CASBZ = CSBZ
            END IF
            IAX  = IPTCNT(NAX,IREP,1)
            IAY  = IPTCNT(NAY,IREP,1)
            IAZ  = IPTCNT(NAZ,IREP,1)
            IBX  = IPTCNT(NBX,IREP,1)
            IBY  = IPTCNT(NBY,IREP,1)
            IBZ  = IPTCNT(NBZ,IREP,1)
            DOAX = IAX .NE. 0
            DOAY = IAY .NE. 0
            DOAZ = IAZ .NE. 0
            DOBX = IBX .NE. 0
            DOBY = IBY .NE. 0
            DOBZ = IBZ .NE. 0
C
C           A-A and A-B Hessian elements:
C
C           First is Ax
C
            IF (DOAX) THEN
                  HESSKE(IAX,IAX) = HESSKE(IAX,IAX) + DERTXX
                  HESFS2(IAX,IAX) = HESFS2(IAX,IAX) + DERSXX
               IF (DOAY) THEN
                  HESSKE(IAX,IAY) = HESSKE(IAX,IAY) + DERTXY
                  HESFS2(IAX,IAY) = HESFS2(IAX,IAY) + DERSXY
               END IF
               IF (DOAZ) THEN
                  HESSKE(IAX,IAZ) = HESSKE(IAX,IAZ) + DERTXZ
                  HESFS2(IAX,IAZ) = HESFS2(IAX,IAZ) + DERSXZ
               END IF
               IF (DOBX) THEN
                  HESSKE(IAX,IBX) = HESSKE(IAX,IBX) - CASBX*DERTXX
                  HESFS2(IAX,IBX) = HESFS2(IAX,IBX) - CASBX*DERSXX
               END IF
               IF (DOBY) THEN
                  HESSKE(IAX,IBY) = HESSKE(IAX,IBY) - CSBY*DERTXY
                  HESFS2(IAX,IBY) = HESFS2(IAX,IBY) - CSBY*DERSXY
               END IF
               IF (DOBZ) THEN
                  HESSKE(IAX,IBZ) = HESSKE(IAX,IBZ) - CSBZ*DERTXZ
                  HESFS2(IAX,IBZ) = HESFS2(IAX,IBZ) - CSBZ*DERSXZ
               END IF
            END IF
C
C           First is Ay
C
            IF (DOAY) THEN
                  HESSKE(IAY,IAY) = HESSKE(IAY,IAY) + DERTYY
                  HESFS2(IAY,IAY) = HESFS2(IAY,IAY) + DERSYY
               IF (DOAZ) THEN
                  HESSKE(IAY,IAZ) = HESSKE(IAY,IAZ) + DERTYZ
                  HESFS2(IAY,IAZ) = HESFS2(IAY,IAZ) + DERSYZ
               END IF
               IF (DOBX) THEN
                  HESSKE(IAY,IBX) = HESSKE(IAY,IBX) - CSBX*DERTXY
                  HESFS2(IAY,IBX) = HESFS2(IAY,IBX) - CSBX*DERSXY
               END IF
               IF (DOBY) THEN
                  HESSKE(IAY,IBY) = HESSKE(IAY,IBY) - CASBY*DERTYY
                  HESFS2(IAY,IBY) = HESFS2(IAY,IBY) - CASBY*DERSYY
               END IF
               IF (DOBZ) THEN
                  HESSKE(IAY,IBZ) = HESSKE(IAY,IBZ) - CSBZ*DERTYZ
                  HESFS2(IAY,IBZ) = HESFS2(IAY,IBZ) - CSBZ*DERSYZ
               END IF
            END IF
C
C           First is Az
C
            IF (DOAZ) THEN
                  HESSKE(IAZ,IAZ) = HESSKE(IAZ,IAZ) + DERTZZ
                  HESFS2(IAZ,IAZ) = HESFS2(IAZ,IAZ) + DERSZZ
               IF (DOBX) THEN
                  HESSKE(IAZ,IBX) = HESSKE(IAZ,IBX) - CSBX*DERTXZ
                  HESFS2(IAZ,IBX) = HESFS2(IAZ,IBX) - CSBX*DERSXZ
               END IF
               IF (DOBY) THEN
                  HESSKE(IAZ,IBY) = HESSKE(IAZ,IBY) - CSBY*DERTYZ
                  HESFS2(IAZ,IBY) = HESFS2(IAZ,IBY) - CSBY*DERSYZ
               END IF
               IF (DOBZ) THEN
                  HESSKE(IAZ,IBZ) = HESSKE(IAZ,IBZ) - CASBZ*DERTZZ
                  HESFS2(IAZ,IBZ) = HESFS2(IAZ,IBZ) - CASBZ*DERSZZ
               END IF
            END IF
C
C           B-B Hessian elements:
C
C           First is Bx
C
            IF (DOBX) THEN
                  HESSKE(IBX,IBX) = HESSKE(IBX,IBX) + DERTXX
                  HESFS2(IBX,IBX) = HESFS2(IBX,IBX) + DERSXX
               IF (DOBY) THEN
                  HESSKE(IBX,IBY) = HESSKE(IBX,IBY) + CSBX*CSBY*DERTXY
                  HESFS2(IBX,IBY) = HESFS2(IBX,IBY) + CSBX*CSBY*DERSXY
               END IF
               IF (DOBZ) THEN
                  HESSKE(IBX,IBZ) = HESSKE(IBX,IBZ) + CSBX*CSBZ*DERTXZ
                  HESFS2(IBX,IBZ) = HESFS2(IBX,IBZ) + CSBX*CSBZ*DERSXZ
               END IF
            END IF
C
C           First is By
C
            IF (DOBY) THEN
                  HESSKE(IBY,IBY) = HESSKE(IBY,IBY) + DERTYY
                  HESFS2(IBY,IBY) = HESFS2(IBY,IBY) + DERSYY
               IF (DOBZ) THEN
                  HESSKE(IBY,IBZ) = HESSKE(IBY,IBZ) + CSBY*CSBZ*DERTYZ
                  HESFS2(IBY,IBZ) = HESFS2(IBY,IBZ) + CSBY*CSBZ*DERSYZ
               END IF
            END IF
C
C           First is Bz
C
            IF (DOBZ) THEN
                  HESSKE(IBZ,IBZ) = HESSKE(IBZ,IBZ) + DERTZZ
                  HESFS2(IBZ,IBZ) = HESFS2(IBZ,IBZ) + DERSZZ
            END IF
         END IF
  300    CONTINUE
      END IF
      RETURN
      END
C  /* Deck avena2 */
      SUBROUTINE AVENA2(ADER,NATOMC,ISYMOP,SECDER,NCENTC,MAXCMP,JCENTC,
     &                  JSYMC,SIGNC,DSHELL,HESSNA)
C
C     tuh 1984
C     modified for symmetry tuh Jul 5 1988
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.00 D00, D1 = 1.00 D00, D2 = 2.00D00)
      LOGICAL SECDER, DOAX, DOAY, DOAZ, DOBX, DOBY, DOBZ,
     &        DOCX, DOCY, DOCZ
      DIMENSION ADER(KCKTAB,NATOMC,*), DSHELL(KHKTAB,*), NCENTC(*),
     &          JCENTC(*), SIGNC(3,*), JSYMC(*), HESSNA(MXCOOR,MXCOOR)
C
#include "onecom.h"
#include "ader.h"
#include "symmet.h"
#include "dorps.h"
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
#include "ibtfun.h"
C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
      NBX = 3*NCENTB - 2
      NBY = 3*NCENTB - 1
      NBZ = 3*NCENTB
      FAB = D1
      IF (NCENTA .EQ. NCENTB) FAB = D2
C
C    Loop over atoms in operator
C
      DO 200 IATOM = 1,NATOMC
         ICENTC = NCENTC(IATOM)
         ISYMPC = JSYMC(IATOM)
         KCENTC = JCENTC(IATOM)
         NCX    = 3*KCENTC - 2
         NCY    = 3*KCENTC - 1
         NCZ    = 3*KCENTC
         SCX    = SIGNC(1,IATOM)
         SCY    = SIGNC(2,IATOM)
         SCZ    = SIGNC(3,IATOM)
C
C        Multiply densities and integrals
C
         DA0000 = D0
         DA0X00 = D0
         DA0Y00 = D0
         DA0Z00 = D0
         DA000X = D0
         DA000Y = D0
         DA000Z = D0
         IF (SECDER) THEN
            DAXX00 = D0
            DAXY00 = D0
            DAXZ00 = D0
            DAYY00 = D0
            DAYZ00 = D0
            DAZZ00 = D0
            DA00XX = D0
            DA00XY = D0
            DA00XZ = D0
            DA00YY = D0
            DA00YZ = D0
            DA00ZZ = D0
            DA0X0X = D0
            DA0X0Y = D0
            DA0X0Z = D0
            DA0Y0X = D0
            DA0Y0Y = D0
            DA0Y0Z = D0
            DA0Z0X = D0
            DA0Z0Y = D0
            DA0Z0Z = D0
         END IF
         DO 300 ICOMP = 1, KHKTAB
            DENS = DSHELL(ICOMP,1)
            DA0000 = DA0000 + DENS*ADER(ICOMP,IATOM,IA0000)
            DA0X00 = DA0X00 + DENS*ADER(ICOMP,IATOM,IA0X00)
            DA0Y00 = DA0Y00 + DENS*ADER(ICOMP,IATOM,IA0Y00)
            DA0Z00 = DA0Z00 + DENS*ADER(ICOMP,IATOM,IA0Z00)
            DA000X = DA000X + DENS*ADER(ICOMP,IATOM,IA000X)
            DA000Y = DA000Y + DENS*ADER(ICOMP,IATOM,IA000Y)
            DA000Z = DA000Z + DENS*ADER(ICOMP,IATOM,IA000Z)
            IF (SECDER) THEN
               DAXX00 = DAXX00 + DENS*ADER(ICOMP,IATOM,IAXX00)
               DAXY00 = DAXY00 + DENS*ADER(ICOMP,IATOM,IAXY00)
               DAXZ00 = DAXZ00 + DENS*ADER(ICOMP,IATOM,IAXZ00)
               DAYY00 = DAYY00 + DENS*ADER(ICOMP,IATOM,IAYY00)
               DAYZ00 = DAYZ00 + DENS*ADER(ICOMP,IATOM,IAYZ00)
               DAZZ00 = DAZZ00 + DENS*ADER(ICOMP,IATOM,IAZZ00)
               DA00XX = DA00XX + DENS*ADER(ICOMP,IATOM,IA00XX)
               DA00XY = DA00XY + DENS*ADER(ICOMP,IATOM,IA00XY)
               DA00XZ = DA00XZ + DENS*ADER(ICOMP,IATOM,IA00XZ)
               DA00YY = DA00YY + DENS*ADER(ICOMP,IATOM,IA00YY)
               DA00YZ = DA00YZ + DENS*ADER(ICOMP,IATOM,IA00YZ)
               DA00ZZ = DA00ZZ + DENS*ADER(ICOMP,IATOM,IA00ZZ)
               DA0X0X = DA0X0X + DENS*ADER(ICOMP,IATOM,IA0X0X)
               DA0X0Y = DA0X0Y + DENS*ADER(ICOMP,IATOM,IA0X0Y)
               DA0X0Z = DA0X0Z + DENS*ADER(ICOMP,IATOM,IA0X0Z)
               DA0Y0X = DA0Y0X + DENS*ADER(ICOMP,IATOM,IA0Y0X)
               DA0Y0Y = DA0Y0Y + DENS*ADER(ICOMP,IATOM,IA0Y0Y)
               DA0Y0Z = DA0Y0Z + DENS*ADER(ICOMP,IATOM,IA0Y0Z)
               DA0Z0X = DA0Z0X + DENS*ADER(ICOMP,IATOM,IA0Z0X)
               DA0Z0Y = DA0Z0Y + DENS*ADER(ICOMP,IATOM,IA0Z0Y)
               DA0Z0Z = DA0Z0Z + DENS*ADER(ICOMP,IATOM,IA0Z0Z)
            END IF
  300    CONTINUE
C
C        Undifferentiated Nuclear Attraction Energy
C
         ENERNA = ENERNA + DA0000
C
C        ***** two-center case: C = A *****
C
         IF (ICENTC .EQ. ICENTA) THEN
            IF (DOREPS(0)) THEN
               TERMX  = DA0X00 + SCX*DA000X
               TERMY  = DA0Y00 + SCY*DA000Y
               TERMZ  = DA0Z00 + SCZ*DA000Z
               IAX    = IPTCNT(NAX,0,1)
               IAY    = IPTCNT(NAY,0,1)
               IAZ    = IPTCNT(NAZ,0,1)
               IBX    = IPTCNT(NBX,0,1)
               IBY    = IPTCNT(NBY,0,1)
               IBZ    = IPTCNT(NBZ,0,1)
               DOAX   = IAX .NE. 0
               DOAY   = IAY .NE. 0
               DOAZ   = IAZ .NE. 0
               DOBX   = IBX .NE. 0
               DOBY   = IBY .NE. 0
               DOBZ   = IBZ .NE. 0
C
C              A nuclear-attraction gradient elements
C
               IF (DOAX) GRADNA(IAX) = GRADNA(IAX) + TERMX
               IF (DOAY) GRADNA(IAY) = GRADNA(IAY) + TERMY
               IF (DOAZ) GRADNA(IAZ) = GRADNA(IAZ) + TERMZ
C
C              B nuclear-attraction gradient elements
C
               IF (DOBX) GRADNA(IBX) = GRADNA(IBX) - SIGNBX*TERMX
               IF (DOBY) GRADNA(IBY) = GRADNA(IBY) - SIGNBY*TERMY
               IF (DOBZ) GRADNA(IBZ) = GRADNA(IBZ) - SIGNBZ*TERMZ
            END IF
C
C           Second Derivatives:
C
            IF (SECDER) THEN
               TERMXX = DAXX00 + DA0X0X + DA0X0X + DA00XX
               TERMXY = DAXY00 + DA0X0Y + DA0Y0X + DA00XY
               TERMXZ = DAXZ00 + DA0X0Z + DA0Z0X + DA00XZ
               TERMYY = DAYY00 + DA0Y0Y + DA0Y0Y + DA00YY
               TERMYZ = DAYZ00 + DA0Y0Z + DA0Z0Y + DA00YZ
               TERMZZ = DAZZ00 + DA0Z0Z + DA0Z0Z + DA00ZZ
               DO 400 IREP = 0, MAXREP
               IF (DOREPS(IREP)) THEN
                  CHIB = PT(IBTAND(ISYMOP,IREP))
                  CSBX = CHIB*SIGNBX
                  CSBY = CHIB*SIGNBY
                  CSBZ = CHIB*SIGNBZ
                  IAX  = IPTCNT(NAX,IREP,1)
                  IAY  = IPTCNT(NAY,IREP,1)
                  IAZ  = IPTCNT(NAZ,IREP,1)
                  IBX  = IPTCNT(NBX,IREP,1)
                  IBY  = IPTCNT(NBY,IREP,1)
                  IBZ  = IPTCNT(NBZ,IREP,1)
                  DOAX = IAX .NE. 0
                  DOAY = IAY .NE. 0
                  DOAZ = IAZ .NE. 0
                  DOBX = IBX .NE. 0
                  DOBY = IBY .NE. 0
                  DOBZ = IBZ .NE. 0
C
C
C                 A-A nuclear-attraction Hessian elements
C
                  IF (DOAX)
     *               HESSNA(IAX,IAX) = HESSNA(IAX,IAX) + TERMXX
                  IF (DOAX.AND.DOAY)
     *               HESSNA(IAX,IAY) = HESSNA(IAX,IAY) + TERMXY
                  IF (DOAX.AND.DOAZ)
     *               HESSNA(IAX,IAZ) = HESSNA(IAX,IAZ) + TERMXZ
                  IF (DOAY)
     *               HESSNA(IAY,IAY) = HESSNA(IAY,IAY) + TERMYY
                  IF (DOAY.AND.DOAZ)
     *               HESSNA(IAY,IAZ) = HESSNA(IAY,IAZ) + TERMYZ
                  IF (DOAZ)
     *               HESSNA(IAZ,IAZ) = HESSNA(IAZ,IAZ) + TERMZZ
C
C                 B-B nuclear-attraction Hessian elements
C
                  IF (DOBX)
     *               HESSNA(IBX,IBX) = HESSNA(IBX,IBX) + TERMXX
                  IF (DOBX.AND.DOBY)
     *               HESSNA(IBX,IBY) = HESSNA(IBX,IBY)
     &                                   + SIGNBX*SIGNBY*TERMXY
                  IF (DOBX.AND.DOBZ)
     *               HESSNA(IBX,IBZ) = HESSNA(IBX,IBZ)
     &                                   + SIGNBX*SIGNBZ*TERMXZ
                  IF (DOBY)
     *               HESSNA(IBY,IBY) = HESSNA(IBY,IBY) + TERMYY
                  IF (DOBY.AND.DOBZ)
     *               HESSNA(IBY,IBZ) = HESSNA(IBY,IBZ)
     &                                   + SIGNBY*SIGNBZ*TERMYZ
                  IF (DOBZ.AND.DOBZ)
     *               HESSNA(IBZ,IBZ) = HESSNA(IBZ,IBZ) + TERMZZ
C
C                 A-B nuclear-attraction Hessian elements
C
                  IF (DOAX.AND.DOBX)
     *               HESSNA(IAX,IBX) = HESSNA(IAX,IBX) - FAB*CSBX*TERMXX
                  IF (DOAX.AND.DOBY)
     *               HESSNA(IAX,IBY) = HESSNA(IAX,IBY) - CSBY*TERMXY
                  IF (DOAX.AND.DOBZ)
     *               HESSNA(IAX,IBZ) = HESSNA(IAX,IBZ) - CSBZ*TERMXZ
                  IF (DOAY.AND.DOBX)
     *               HESSNA(IAY,IBX) = HESSNA(IAY,IBX) - CSBX*TERMXY
                  IF (DOAY.AND.DOBY)
     *               HESSNA(IAY,IBY) = HESSNA(IAY,IBY) - FAB*CSBY*TERMYY
                  IF (DOAY.AND.DOBZ)
     *               HESSNA(IAY,IBZ) = HESSNA(IAY,IBZ) - CSBZ*TERMYZ
                  IF (DOAZ.AND.DOBX)
     *               HESSNA(IAZ,IBX) = HESSNA(IAZ,IBX) - CSBX*TERMXZ
                  IF (DOAZ.AND.DOBY)
     *               HESSNA(IAZ,IBY) = HESSNA(IAZ,IBY) - CSBY*TERMYZ
                  IF (DOAZ.AND.DOBZ)
     *               HESSNA(IAZ,IBZ) = HESSNA(IAZ,IBZ) - FAB*CSBZ*TERMZZ
               END IF
  400          CONTINUE
            END IF
C
C        ***** two-center case: C = B *****
C
         ELSE IF (ICENTC .EQ. ICENTB) THEN
            IF (DOREPS(0)) THEN
               IAX  = IPTCNT(NAX,0,1)
               IAY  = IPTCNT(NAY,0,1)
               IAZ  = IPTCNT(NAZ,0,1)
               IBX  = IPTCNT(NBX,0,1)
               IBY  = IPTCNT(NBY,0,1)
               IBZ  = IPTCNT(NBZ,0,1)
               DOAX = IAX .NE. 0
               DOAY = IAY .NE. 0
               DOAZ = IAZ .NE. 0
               DOBX = IBX .NE. 0
               DOBY = IBY .NE. 0
               DOBZ = IBZ .NE. 0
C
C              A nuclear-attraction gradient elements
C
               IF (DOAX) GRADNA(IAX) = GRADNA(IAX) + DA0X00
               IF (DOAY) GRADNA(IAY) = GRADNA(IAY) + DA0Y00
               IF (DOAZ) GRADNA(IAZ) = GRADNA(IAZ) + DA0Z00
C
C              B nuclear-attraction gradient elements
C
               IF (DOBX) GRADNA(IBX) = GRADNA(IBX) - SIGNBX*DA0X00
               IF (DOBY) GRADNA(IBY) = GRADNA(IBY) - SIGNBY*DA0Y00
               IF (DOBZ) GRADNA(IBZ) = GRADNA(IBZ) - SIGNBZ*DA0Z00
            END IF
C
C           Second Derivatives:
C
            IF (SECDER) THEN
               DO 500 IREP = 0, MAXREP
               IF (DOREPS(IREP)) THEN
                  CHIB = PT(IBTAND(ISYMOP,IREP))
                  CSBX = CHIB*SIGNBX
                  CSBY = CHIB*SIGNBY
                  CSBZ = CHIB*SIGNBZ
                  IAX  = IPTCNT(NAX,IREP,1)
                  IAY  = IPTCNT(NAY,IREP,1)
                  IAZ  = IPTCNT(NAZ,IREP,1)
                  IBX  = IPTCNT(NBX,IREP,1)
                  IBY  = IPTCNT(NBY,IREP,1)
                  IBZ  = IPTCNT(NBZ,IREP,1)
                  DOAX = IAX .NE. 0
                  DOAY = IAY .NE. 0
                  DOAZ = IAZ .NE. 0
                  DOBX = IBX .NE. 0
                  DOBY = IBY .NE. 0
                  DOBZ = IBZ .NE. 0
C
C                 A-A nuclear-attraction Hessian elements
C
                  IF (DOAX)
     *               HESSNA(IAX,IAX) = HESSNA(IAX,IAX) + DAXX00
                  IF (DOAX.AND.DOAY)
     *               HESSNA(IAX,IAY) = HESSNA(IAX,IAY) + DAXY00
                  IF (DOAX.AND.DOAZ)
     *               HESSNA(IAX,IAZ) = HESSNA(IAX,IAZ) + DAXZ00
                  IF (DOAY)
     *               HESSNA(IAY,IAY) = HESSNA(IAY,IAY) + DAYY00
                  IF (DOAY.AND.DOAZ)
     *               HESSNA(IAY,IAZ) = HESSNA(IAY,IAZ) + DAYZ00
                  IF (DOAZ)
     *               HESSNA(IAZ,IAZ) = HESSNA(IAZ,IAZ) + DAZZ00
C
C                 B-B nuclear-attraction Hessian elements
C
                  IF (DOBX)
     *               HESSNA(IBX,IBX) = HESSNA(IBX,IBX) + DAXX00
                  IF (DOBX.AND.DOBY)
     *               HESSNA(IBX,IBY) = HESSNA(IBX,IBY)
     &                                   + SIGNBX*SIGNBY*DAXY00
                  IF (DOBX.AND.DOBZ)
     *               HESSNA(IBX,IBZ) = HESSNA(IBX,IBZ)
     &                                   + SIGNBX*SIGNBZ*DAXZ00
                  IF (DOBY)
     *               HESSNA(IBY,IBY) = HESSNA(IBY,IBY) + DAYY00
                  IF (DOBY.AND.DOBZ)
     *               HESSNA(IBY,IBZ) = HESSNA(IBY,IBZ)
     &                                   + SIGNBY*SIGNBZ*DAYZ00
                  IF (DOBZ)
     *               HESSNA(IBZ,IBZ) = HESSNA(IBZ,IBZ) + DAZZ00
C
C                 A-B nuclear-attraction Hessian elements
C
                  IF (DOAX.AND.DOBX)
     *               HESSNA(IAX,IBX) = HESSNA(IAX,IBX) - FAB*CSBX*DAXX00
                  IF (DOAX.AND.DOBY)
     *               HESSNA(IAX,IBY) = HESSNA(IAX,IBY) - CSBY*DAXY00
                  IF (DOAX.AND.DOBZ)
     *               HESSNA(IAX,IBZ) = HESSNA(IAX,IBZ) - CSBZ*DAXZ00
                  IF (DOAY.AND.DOBX)
     *               HESSNA(IAY,IBX) = HESSNA(IAY,IBX) - CSBX*DAXY00
                  IF (DOAY.AND.DOBY)
     *               HESSNA(IAY,IBY) = HESSNA(IAY,IBY) - FAB*CSBY*DAYY00
                  IF (DOAY.AND.DOBZ)
     *               HESSNA(IAY,IBZ) = HESSNA(IAY,IBZ) - CSBZ*DAYZ00
                  IF (DOAZ.AND.DOBX)
     *               HESSNA(IAZ,IBX) = HESSNA(IAZ,IBX) - CSBX*DAXZ00
                  IF (DOAZ.AND.DOBY)
     *               HESSNA(IAZ,IBY) = HESSNA(IAZ,IBY) - CSBY*DAYZ00
                  IF (DOAZ.AND.DOBZ)
     *               HESSNA(IAZ,IBZ) = HESSNA(IAZ,IBZ) - FAB*CSBZ*DAZZ00
               END IF
  500          CONTINUE
            END IF
C
C        ***** three-center case *****
C
         ELSE
            IF (DOREPS(0)) THEN
               IAX  = IPTCNT(NAX,0,1)
               IAY  = IPTCNT(NAY,0,1)
               IAZ  = IPTCNT(NAZ,0,1)
               IBX  = IPTCNT(NBX,0,1)
               IBY  = IPTCNT(NBY,0,1)
               IBZ  = IPTCNT(NBZ,0,1)
               ICX  = IPTCNT(NCX,0,1)
               ICY  = IPTCNT(NCY,0,1)
               ICZ  = IPTCNT(NCZ,0,1)
               DOAX = IAX .NE. 0
               DOAY = IAY .NE. 0
               DOAZ = IAZ .NE. 0
               DOBX = IBX .NE. 0
               DOBY = IBY .NE. 0
               DOBZ = IBZ .NE. 0
               DOCX = ICX .NE. 0
               DOCY = ICY .NE. 0
               DOCZ = ICZ .NE. 0
C
C              A nuclear-attraction gradient elements:
C
               IF (DOAX) GRADNA(IAX) = GRADNA(IAX) + DA0X00
               IF (DOAY) GRADNA(IAY) = GRADNA(IAY) + DA0Y00
               IF (DOAZ) GRADNA(IAZ) = GRADNA(IAZ) + DA0Z00
C
C              B nuclear-attraction gradient elements:
C
               IF (DOBX) GRADNA(IBX) = GRADNA(IBX)
     *                               - SIGNBX*(DA0X00 + DA000X)
               IF (DOBY) GRADNA(IBY) = GRADNA(IBY)
     *                               - SIGNBY*(DA0Y00 + DA000Y)
               IF (DOBZ) GRADNA(IBZ) = GRADNA(IBZ)
     *                               - SIGNBZ*(DA0Z00 + DA000Z)
C
C              C nuclear-attraction gradient elements:
C
               IF (DOCX) GRADNA(ICX) = GRADNA(ICX) + SCX*DA000X
               IF (DOCY) GRADNA(ICY) = GRADNA(ICY) + SCY*DA000Y
               IF (DOCZ) GRADNA(ICZ) = GRADNA(ICZ) + SCZ*DA000Z
            END IF
C
C           Second Derivatives:
C
            IF (SECDER) THEN
               FAC = D1
               IF (NCENTA .EQ. KCENTC) FAC = D2
               FBC = D1
               IF (NCENTB .EQ. KCENTC) FBC = D2
               DO 600 IREP = 0, MAXREP
               IF (DOREPS(IREP)) THEN
                  CHIB = PT(IBTAND(ISYMOP,IREP))
                  CHIC = PT(IBTAND(ISYMPC,IREP))
                  CSBX = CHIB*SIGNBX
                  CSBY = CHIB*SIGNBY
                  CSBZ = CHIB*SIGNBZ
                  CSCX = CHIC*SCX
                  CSCY = CHIC*SCY
                  CSCZ = CHIC*SCZ
                  IAX  = IPTCNT(NAX,IREP,1)
                  IAY  = IPTCNT(NAY,IREP,1)
                  IAZ  = IPTCNT(NAZ,IREP,1)
                  IBX  = IPTCNT(NBX,IREP,1)
                  IBY  = IPTCNT(NBY,IREP,1)
                  IBZ  = IPTCNT(NBZ,IREP,1)
                  ICX  = IPTCNT(NCX,IREP,1)
                  ICY  = IPTCNT(NCY,IREP,1)
                  ICZ  = IPTCNT(NCZ,IREP,1)
                  DOAX = IAX .NE. 0
                  DOAY = IAY .NE. 0
                  DOAZ = IAZ .NE. 0
                  DOBX = IBX .NE. 0
                  DOBY = IBY .NE. 0
                  DOBZ = IBZ .NE. 0
                  DOCX = ICX .NE. 0
                  DOCY = ICY .NE. 0
                  DOCZ = ICZ .NE. 0
C
C                 A-A nuclear-attraction Hessian elements:
C
                  IF (DOAX)
     *               HESSNA(IAX,IAX) = HESSNA(IAX,IAX) + DAXX00
                  IF (DOAX.AND.DOAY)
     *               HESSNA(IAX,IAY) = HESSNA(IAX,IAY) + DAXY00
                  IF (DOAX.AND.DOAZ)
     *               HESSNA(IAX,IAZ) = HESSNA(IAX,IAZ) + DAXZ00
                  IF (DOAY)
     *               HESSNA(IAY,IAY) = HESSNA(IAY,IAY) + DAYY00
                  IF (DOAY.AND.DOAZ)
     *               HESSNA(IAY,IAZ) = HESSNA(IAY,IAZ) + DAYZ00
                  IF (DOAZ)
     *               HESSNA(IAZ,IAZ) = HESSNA(IAZ,IAZ) + DAZZ00
C
C                 A-B nuclear-attraction Hessian elements:
C
                  IF (DOAX.AND.DOBX)
     *               HESSNA(IAX,IBX) = HESSNA(IAX,IBX)
     *                               - FAB*CSBX*(DAXX00 + DA0X0X)
                  IF (DOAX.AND.DOBY)
     *               HESSNA(IAX,IBY) = HESSNA(IAX,IBY)
     *                               - CSBY*(DAXY00 + DA0X0Y)
                  IF (DOAX.AND.DOBZ)
     *               HESSNA(IAX,IBZ) = HESSNA(IAX,IBZ)
     *                               - CSBZ*(DAXZ00 + DA0X0Z)
                  IF (DOAY.AND.DOBX)
     *               HESSNA(IAY,IBX) = HESSNA(IAY,IBX)
     *                               - CSBX*(DAXY00 + DA0Y0X)
                  IF (DOAY.AND.DOBY)
     *               HESSNA(IAY,IBY) = HESSNA(IAY,IBY)
     *                               - FAB*CSBY*(DAYY00 + DA0Y0Y)
                  IF (DOAY.AND.DOBZ)
     *               HESSNA(IAY,IBZ) = HESSNA(IAY,IBZ)
     *                               - CSBZ*(DAYZ00 + DA0Y0Z)
                  IF (DOAZ.AND.DOBX)
     *               HESSNA(IAZ,IBX) = HESSNA(IAZ,IBX)
     *                               - CSBX*(DAXZ00 + DA0Z0X)
                  IF (DOAZ.AND.DOBY)
     *               HESSNA(IAZ,IBY) = HESSNA(IAZ,IBY)
     *                               - CSBY*(DAYZ00 + DA0Z0Y)
                  IF (DOAZ.AND.DOBZ)
     *               HESSNA(IAZ,IBZ) = HESSNA(IAZ,IBZ)
     *                               - FAB*CSBZ*(DAZZ00 + DA0Z0Z)
C
C                 A-C nuclear-attraction Hessian elements:
C
                  IF (DOAX.AND.DOCX)
     *               HESSNA(IAX,ICX) = HESSNA(IAX,ICX) + FAC*CSCX*DA0X0X
                  IF (DOAX.AND.DOCY)
     *               HESSNA(IAX,ICY) = HESSNA(IAX,ICY) + CSCY*DA0X0Y
                  IF (DOAX.AND.DOCZ)
     *               HESSNA(IAX,ICZ) = HESSNA(IAX,ICZ) + CSCZ*DA0X0Z
                  IF (DOAY.AND.DOCX)
     *               HESSNA(IAY,ICX) = HESSNA(IAY,ICX) + CSCX*DA0Y0X
                  IF (DOAY.AND.DOCY)
     *               HESSNA(IAY,ICY) = HESSNA(IAY,ICY) + FAC*CSCY*DA0Y0Y
                  IF (DOAY.AND.DOCZ)
     *               HESSNA(IAY,ICZ) = HESSNA(IAY,ICZ) + CSCZ*DA0Y0Z
                  IF (DOAZ.AND.DOCX)
     *               HESSNA(IAZ,ICX) = HESSNA(IAZ,ICX) + CSCX*DA0Z0X
                  IF (DOAZ.AND.DOCY)
     *               HESSNA(IAZ,ICY) = HESSNA(IAZ,ICY) + CSCY*DA0Z0Y
                  IF (DOAZ.AND.DOCZ)
     *               HESSNA(IAZ,ICZ) = HESSNA(IAZ,ICZ) + FAC*CSCZ*DA0Z0Z
C
C                 B-B nuclear-attraction Hessian elements:
C
                  IF (DOBX)
     *               HESSNA(IBX,IBX) = HESSNA(IBX,IBX) +
     *               (DAXX00 + DA00XX + DA0X0X + DA0X0X)
                  IF (DOBX.AND.DOBY)
     *               HESSNA(IBX,IBY) = HESSNA(IBX,IBY)
     &                                   + SIGNBX*SIGNBY*
     *               (DAXY00 + DA00XY + DA0X0Y + DA0Y0X)
                  IF (DOBX.AND.DOBZ)
     *               HESSNA(IBX,IBZ) = HESSNA(IBX,IBZ)
     &                                   + SIGNBX*SIGNBZ*
     *               (DAXZ00 + DA00XZ + DA0X0Z + DA0Z0X)
                  IF (DOBY)
     *               HESSNA(IBY,IBY) = HESSNA(IBY,IBY) +
     *               (DAYY00 + DA00YY + DA0Y0Y + DA0Y0Y)
                  IF (DOBY.AND.DOBZ)
     *               HESSNA(IBY,IBZ) = HESSNA(IBY,IBZ)
     &                                   + SIGNBY*SIGNBZ*
     *               (DAYZ00 + DA00YZ + DA0Y0Z + DA0Z0Y)
                  IF (DOBZ)
     *               HESSNA(IBZ,IBZ) = HESSNA(IBZ,IBZ) +
     *               (DAZZ00 + DA00ZZ + DA0Z0Z + DA0Z0Z)
C
C                 B-C nuclear-attraction Hessian elements:
C
                  IF (DOBX.AND.DOCX) HESSNA(IBX,ICX) = HESSNA(IBX,ICX)
     *                  - FBC*CSBX*CSCX*(DA0X0X + DA00XX)
                  IF (DOBX.AND.DOCY) HESSNA(IBX,ICY) = HESSNA(IBX,ICY)
     *                  - CSBX*CSCY*(DA0X0Y + DA00XY)
                  IF (DOBX.AND.DOCZ) HESSNA(IBX,ICZ) = HESSNA(IBX,ICZ)
     *                  - CSBX*CSCZ*(DA0X0Z + DA00XZ)
                  IF (DOBY.AND.DOCX) HESSNA(IBY,ICX) = HESSNA(IBY,ICX)
     *                  - CSBY*CSCX*(DA0Y0X + DA00XY)
                  IF (DOBY.AND.DOCY) HESSNA(IBY,ICY) = HESSNA(IBY,ICY)
     *                  - FBC*CSBY*CSCY*(DA0Y0Y + DA00YY)
                  IF (DOBY.AND.DOCZ) HESSNA(IBY,ICZ) = HESSNA(IBY,ICZ)
     *                  - CSBY*CSCZ*(DA0Y0Z + DA00YZ)
                  IF (DOBZ.AND.DOCX) HESSNA(IBZ,ICX) = HESSNA(IBZ,ICX)
     *                  - CSBZ*CSCX*(DA0Z0X + DA00XZ)
                  IF (DOBZ.AND.DOCY) HESSNA(IBZ,ICY) = HESSNA(IBZ,ICY)
     *                  - CSBZ*CSCY*(DA0Z0Y + DA00YZ)
                  IF (DOBZ.AND.DOCZ) HESSNA(IBZ,ICZ) = HESSNA(IBZ,ICZ)
     *                  - FBC*CSBZ*CSCZ*(DA0Z0Z + DA00ZZ)
C
C                 C-C nuclear-attraction Hessian elements:
C
                  IF (DOCX)
     *               HESSNA(ICX,ICX) = HESSNA(ICX,ICX) + DA00XX
                  IF (DOCX.AND.DOCY)
     *               HESSNA(ICX,ICY) = HESSNA(ICX,ICY) + SCX*SCY*DA00XY
                  IF (DOCX.AND.DOCZ)
     *               HESSNA(ICX,ICZ) = HESSNA(ICX,ICZ) + SCX*SCZ*DA00XZ
                  IF (DOCY)
     *               HESSNA(ICY,ICY) = HESSNA(ICY,ICY) + DA00YY
                  IF (DOCY.AND.DOCZ)
     *               HESSNA(ICY,ICZ) = HESSNA(ICY,ICZ) + SCY*SCZ*DA00YZ
                  IF (DOCZ)
     *               HESSNA(ICZ,ICZ) = HESSNA(ICZ,ICZ) + DA00ZZ
               END IF
  600          CONTINUE
            END IF
         END IF
  200 CONTINUE
      RETURN
      END
C  /* Deck avedip */
      SUBROUTINE AVEDIP(SINT0,DINT1,ISYMOP,DSHELL,MAXCMP)
C
C     tuh 1985
C     symmetry 081288 tuh
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (D0 = 0.00 D00)
      DIMENSION DSHELL(KHKTAB), SINT0(KCKTAB), DINT1(KCKTAB,3,3)
#include "onecom.h"
#include "symmet.h"
#include "dorps.h"
#include "dipole.h"
C
#include "ibtfun.h"
C
         NAX = 3*NCENTA - 2
         NAY = 3*NCENTA - 1
         NAZ = 3*NCENTA
         NBX = 3*NCENTB - 2
         NBY = 3*NCENTB - 1
         NBZ = 3*NCENTB
         IX  = IPTAX(1,1)
         IY  = IPTAX(2,1)
         IZ  = IPTAX(3,1)
C
C        Expectation value of undifferentiated elements
C
         SAVR0  = D0
#if defined (SYS_HAL)
#include <scalar.h>
#endif
         DO 200 ICOMP = 1, KHKTAB
            DENS  = DSHELL(ICOMP)
            SAVR0 = SAVR0 + DENS*SINT0(ICOMP)
 200     CONTINUE
         IF (ONECEN) THEN
            DO 300 IREP = 0, MAXREP
            IF (DOREPS(IREP)) THEN
               IAX = IPTCNT(NAX,IREP,1)
               IAY = IPTCNT(NAY,IREP,1)
               IAZ = IPTCNT(NAZ,IREP,1)
               IF (ISYMAX(1,1) .EQ. IREP .AND. IAX .GT. 0) THEN
                  DDIPE(IX,IAX) = DDIPE(IX,IAX) - SAVR0
               END IF
               IF (ISYMAX(2,1) .EQ. IREP .AND. IAY .GT. 0) THEN
                  DDIPE(IY,IAY) = DDIPE(IY,IAY) - SAVR0
               END IF
               IF (ISYMAX(3,1) .EQ. IREP .AND. IAZ .GT. 0) THEN
                  DDIPE(IZ,IAZ) = DDIPE(IZ,IAZ) - SAVR0
               END IF
            END IF
  300       CONTINUE
         ELSE
            XAVRX  = D0
            XAVRY  = D0
            XAVRZ  = D0
            YAVRX  = D0
            YAVRY  = D0
            YAVRZ  = D0
            ZAVRX  = D0
            ZAVRY  = D0
            ZAVRZ  = D0
            DO 400 ICOMP = 1, KHKTAB
               DENS  = DSHELL(ICOMP)
               XAVRX = XAVRX - DENS*DINT1(ICOMP,1,1)
               XAVRY = XAVRY - DENS*DINT1(ICOMP,1,2)
               XAVRZ = XAVRZ - DENS*DINT1(ICOMP,1,3)
               YAVRX = YAVRX - DENS*DINT1(ICOMP,2,1)
               YAVRY = YAVRY - DENS*DINT1(ICOMP,2,2)
               YAVRZ = YAVRZ - DENS*DINT1(ICOMP,2,3)
               ZAVRX = ZAVRX - DENS*DINT1(ICOMP,3,1)
               ZAVRY = ZAVRY - DENS*DINT1(ICOMP,3,2)
               ZAVRZ = ZAVRZ - DENS*DINT1(ICOMP,3,3)
  400       CONTINUE
            DO 500 IREP = 0, MAXREP
            IF (DOREPS(IREP)) THEN
               CHI  = PT(IBTAND(ISYMOP,IREP))
               CSBX = CHI*SIGNBX
               CSBY = CHI*SIGNBY
               CSBZ = CHI*SIGNBZ
               IAX = IPTCNT(NAX,IREP,1)
               IAY = IPTCNT(NAY,IREP,1)
               IAZ = IPTCNT(NAZ,IREP,1)
               IBX = IPTCNT(NBX,IREP,1)
               IBY = IPTCNT(NBY,IREP,1)
               IBZ = IPTCNT(NBZ,IREP,1)
               IF (ISYMAX(1,1) .EQ. IREP) THEN
                  IF (IAX.GT.0)
     *               DDIPE(IX,IAX) = DDIPE(IX,IAX) + XAVRX
                  IF (IAY.GT.0)
     *               DDIPE(IX,IAY) = DDIPE(IX,IAY) + XAVRY
                  IF (IAZ.GT.0)
     *               DDIPE(IX,IAZ) = DDIPE(IX,IAZ) + XAVRZ
                  IF (IBX.GT.0)
     *               DDIPE(IX,IBX) = DDIPE(IX,IBX) - CSBX*XAVRX - SAVR0
                  IF (IBY.GT.0)
     *               DDIPE(IX,IBY) = DDIPE(IX,IBY) - CSBY*XAVRY
                  IF (IBZ.GT.0)
     *               DDIPE(IX,IBZ) = DDIPE(IX,IBZ) - CSBZ*XAVRZ
               END IF
               IF (ISYMAX(2,1) .EQ. IREP) THEN
                  IF (IAX.GT.0)
     *               DDIPE(IY,IAX) = DDIPE(IY,IAX) + YAVRX
                  IF (IAY.GT.0)
     *               DDIPE(IY,IAY) = DDIPE(IY,IAY) + YAVRY
                  IF (IAZ.GT.0)
     *               DDIPE(IY,IAZ) = DDIPE(IY,IAZ) + YAVRZ
                  IF (IBX.GT.0)
     *               DDIPE(IY,IBX) = DDIPE(IY,IBX) - CSBX*YAVRX
                  IF (IBY.GT.0)
     *               DDIPE(IY,IBY) = DDIPE(IY,IBY) - CSBY*YAVRY - SAVR0
                  IF (IBZ.GT.0)
     *               DDIPE(IY,IBZ) = DDIPE(IY,IBZ) - CSBZ*YAVRZ
               END IF
               IF (ISYMAX(3,1) .EQ. IREP) THEN
                  IF (IAX.GT.0)
     *               DDIPE(IZ,IAX) = DDIPE(IZ,IAX) + ZAVRX
                  IF (IAY.GT.0)
     *               DDIPE(IZ,IAY) = DDIPE(IZ,IAY) + ZAVRY
                  IF (IAZ.GT.0)
     *               DDIPE(IZ,IAZ) = DDIPE(IZ,IAZ) + ZAVRZ
                  IF (IBX.GT.0)
     *               DDIPE(IZ,IBX) = DDIPE(IZ,IBX) - CSBX*ZAVRX
                  IF (IBY.GT.0)
     *               DDIPE(IZ,IBY) = DDIPE(IZ,IBY) - CSBY*ZAVRY
                  IF (IBZ.GT.0)
     *               DDIPE(IZ,IBZ) = DDIPE(IZ,IBZ) - CSBZ*ZAVRZ - SAVR0
               END IF
            END IF
  500       CONTINUE
         END IF
      RETURN
      END
C  /* Deck aveexp */
      SUBROUTINE AVEEXP(EXPDER,DSHELL,FSHELL)
C
C     tuh & fj feb 2002
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.00 D00)
      DIMENSION DSHELL(KHKTAB), FSHELL(KHKTAB), 
     &          EXPDER(KCKTAB,2,2,NUCA,NUCB)
#include <onecom.h>
#include <symmet.h>
#include <expopt.h>
C
      DO I = 1, NUCA
      DO J = 1, NUCB
         DERA  = D0
         DERB  = D0
         DO K = 1, KHKTAB
            DENS  = DSHELL(K)
            FOCK  = FSHELL(K)
            DERA = DERA - FOCK*EXPDER(K,1,1,I,J) 
     &                  + DENS*EXPDER(K,1,2,I,J)
            DERB = DERB - FOCK*EXPDER(K,2,1,I,J) 
     &                  + DENS*EXPDER(K,2,2,I,J)
         END DO
         ALPGRD(JSTA+I) = ALPGRD(JSTA+I) + DERA
         ALPGRD(JSTB+J) = ALPGRD(JSTB+J) + DERB
      END DO
      END DO
      RETURN
      END
C  /* Deck aveqdp */
      SUBROUTINE AVEQDP(SINT1,QINT1,ISYMOP,DSHELL,MAXCMP)
C
C     K.Ruud, spring 2003
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.00 D00, D2 = 2.0D0)
      DIMENSION DSHELL(KHKTAB), SINT1(KCKTAB,3), QINT1(KCKTAB,3,6)
#include <onecom.h>
#include <symmet.h>
#include <dorps.h>
#include <cbiqpg.h>
#include <difsec.h>
C
#include <ibtfun.h>
C
         NAX = 3*NCENTA - 2
         NAY = 3*NCENTA - 1
         NAZ = 3*NCENTA
         NBX = 3*NCENTB - 2
         NBY = 3*NCENTB - 1
         NBZ = 3*NCENTB
         IX  = IPTAX(1,1)
         IY  = IPTAX(2,1)
         IZ  = IPTAX(3,1)
C
C        Expectation value of undifferentiated elements
C
         SAVRX  = D0
         SAVRY  = D0
         SAVRZ  = D0
#if defined (SYS_CONVEX) || defined (SYS_HAL)
#include <scalar.h>
#endif
         DO 200 ICOMP = 1, KHKTAB
            DENS  = DSHELL(ICOMP)
            SAVRX = SAVRX + DENS*SINT1(ICOMP,1)
            SAVRY = SAVRY + DENS*SINT1(ICOMP,2)
            SAVRZ = SAVRZ + DENS*SINT1(ICOMP,3)
 200     CONTINUE
         IF (ONECEN) THEN
            DO 300 IREP = 0, MAXREP
            IF (DOREPS(IREP)) THEN
               IAX = IPTCNT(NAX,IREP,1)
               IAY = IPTCNT(NAY,IREP,1)
               IAZ = IPTCNT(NAZ,IREP,1)
               IF (IREP .EQ. 0) THEN
                  IF (IAX .GT. 0)
     &                 DSECE(IX,IX,IAX) = DSECE(IX,IX,IAX) - D2*SAVRX
                  IF (IAY .GT. 0)
     &                 DSECE(IY,IY,IAY) = DSECE(IY,IY,IAY) - D2*SAVRY
                  IF (IAZ .GT. 0)
     &                 DSECE(IZ,IZ,IAZ) = DSECE(IZ,IZ,IAZ) - D2*SAVRZ
               END IF
               IF (ISYMAX(3,2) .EQ. IREP .AND. IAX .GT. 0) THEN
                  DSECE(IX,IY,IAX) = DSECE(IX,IY,IAX) - SAVRY
                  DSECE(IY,IX,IAX) = DSECE(IY,IX,IAX) - SAVRY
               END IF
               IF (ISYMAX(3,2) .EQ. IREP .AND. IAY .GT. 0) THEN
                  DSECE(IX,IY,IAY) = DSECE(IX,IY,IAY) - SAVRX
                  DSECE(IY,IX,IAY) = DSECE(IY,IX,IAY) - SAVRX
               END IF
               IF (ISYMAX(2,2) .EQ. IREP .AND. IAX .GT. 0) THEN
                  DSECE(IX,IZ,IAX) = DSECE(IX,IZ,IAX) - SAVRZ
                  DSECE(IZ,IX,IAX) = DSECE(IZ,IX,IAX) - SAVRZ
               END IF
               IF (ISYMAX(2,2) .EQ. IREP .AND. IAZ .GT. 0) THEN
                  DSECE(IX,IZ,IAZ) = DSECE(IX,IZ,IAZ) - SAVRX
                  DSECE(IZ,IX,IAZ) = DSECE(IZ,IX,IAZ) - SAVRX
               END IF
               IF (ISYMAX(1,2) .EQ. IREP .AND. IAY .GT. 0) THEN
                  DSECE(IY,IZ,IAY) = DSECE(IY,IZ,IAY) - SAVRZ
                  DSECE(IZ,IY,IAY) = DSECE(IZ,IY,IAY) - SAVRZ
               END IF
               IF (ISYMAX(1,2) .EQ. IREP .AND. IAZ .GT. 0) THEN
                  DSECE(IY,IZ,IAZ) = DSECE(IY,IZ,IAZ) - SAVRY
                  DSECE(IZ,IY,IAZ) = DSECE(IZ,IY,IAZ) - SAVRY
               END IF
            END IF
  300       CONTINUE
         ELSE
            XAVRXX  = D0
            XAVRXY  = D0
            XAVRXZ  = D0
            XAVRYY  = D0
            XAVRYZ  = D0
            XAVRZZ  = D0
            YAVRXX  = D0
            YAVRXY  = D0
            YAVRXZ  = D0
            YAVRYY  = D0
            YAVRYZ  = D0
            YAVRZZ  = D0
            ZAVRXX  = D0
            ZAVRXY  = D0
            ZAVRXZ  = D0
            ZAVRYY  = D0
            ZAVRYZ  = D0
            ZAVRZZ  = D0
            DO 400 ICOMP = 1, KHKTAB
               DENS  = DSHELL(ICOMP)
               XAVRXX = XAVRXX - DENS*QINT1(ICOMP,1,1)
               XAVRXY = XAVRXY - DENS*QINT1(ICOMP,1,2)
               XAVRXZ = XAVRXZ - DENS*QINT1(ICOMP,1,3)
               XAVRYY = XAVRYY - DENS*QINT1(ICOMP,1,4)
               XAVRYZ = XAVRYZ - DENS*QINT1(ICOMP,1,5)
               XAVRZZ = XAVRZZ - DENS*QINT1(ICOMP,1,6)
               YAVRXX = YAVRXX - DENS*QINT1(ICOMP,2,1)
               YAVRXY = YAVRXY - DENS*QINT1(ICOMP,2,2)
               YAVRXZ = YAVRXZ - DENS*QINT1(ICOMP,2,3)
               YAVRYY = YAVRYY - DENS*QINT1(ICOMP,2,4)
               YAVRYZ = YAVRYZ - DENS*QINT1(ICOMP,2,5)
               YAVRZZ = YAVRZZ - DENS*QINT1(ICOMP,2,6)
               ZAVRXX = ZAVRXX - DENS*QINT1(ICOMP,3,1)
               ZAVRXY = ZAVRXY - DENS*QINT1(ICOMP,3,2)
               ZAVRXZ = ZAVRXZ - DENS*QINT1(ICOMP,3,3)
               ZAVRYY = ZAVRYY - DENS*QINT1(ICOMP,3,4)
               ZAVRYZ = ZAVRYZ - DENS*QINT1(ICOMP,3,5)
               ZAVRZZ = ZAVRZZ - DENS*QINT1(ICOMP,3,6)
  400       CONTINUE
            DO 500 IREP = 0, MAXREP
            IF (DOREPS(IREP)) THEN
               CHI  = PT(IBTAND(ISYMOP,IREP))
               CSBX = CHI*SIGNBX
               CSBY = CHI*SIGNBY
               CSBZ = CHI*SIGNBZ
               IAX = IPTCNT(NAX,IREP,1)
               IAY = IPTCNT(NAY,IREP,1)
               IAZ = IPTCNT(NAZ,IREP,1)
               IBX = IPTCNT(NBX,IREP,1)
               IBY = IPTCNT(NBY,IREP,1)
               IBZ = IPTCNT(NBZ,IREP,1)
               IF (IREP .EQ. 0) THEN
                  IF (IAX.GT.0) THEN
                     DSECE(IX,IX,IAX) = DSECE(IX,IX,IAX) + XAVRXX
                     DSECE(IY,IY,IAX) = DSECE(IY,IY,IAX) + XAVRYY
                     DSECE(IZ,IZ,IAX) = DSECE(IZ,IZ,IAX) + XAVRZZ
                  END IF
                  IF (IAY.GT.0) THEN
                     DSECE(IX,IX,IAY) = DSECE(IX,IX,IAY) + YAVRXX
                     DSECE(IY,IY,IAY) = DSECE(IY,IY,IAY) + YAVRYY
                     DSECE(IZ,IZ,IAY) = DSECE(IZ,IZ,IAY) + YAVRZZ
                  END IF
                  IF (IAZ.GT.0) THEN
                     DSECE(IX,IX,IAZ) = DSECE(IX,IX,IAZ) + ZAVRXX
                     DSECE(IY,IY,IAZ) = DSECE(IY,IY,IAZ) + ZAVRYY
                     DSECE(IZ,IZ,IAZ) = DSECE(IZ,IZ,IAZ) + ZAVRZZ
                  END IF
                  IF (IBX.GT.0) THEN
                     DSECE(IX,IX,IBX) = DSECE(IX,IX,IBX) - CSBX*XAVRXX
     &                                - D2*SAVRX
                     DSECE(IY,IY,IBX) = DSECE(IY,IY,IBX) - CSBX*XAVRYY
                     DSECE(IZ,IZ,IBX) = DSECE(IZ,IZ,IBX) - CSBX*XAVRZZ
                  END IF
                  IF (IBY.GT.0) THEN
                     DSECE(IX,IX,IBY) = DSECE(IX,IX,IBY) - CSBY*YAVRXX
                     DSECE(IY,IY,IBY) = DSECE(IY,IY,IBY) - CSBY*YAVRYY
     &                                - D2*SAVRY
                     DSECE(IZ,IZ,IBY) = DSECE(IZ,IZ,IBY) - CSBY*YAVRZZ
                  END IF
                  IF (IBZ.GT.0) THEN
                     DSECE(IX,IX,IBZ) = DSECE(IX,IX,IBZ) - CSBZ*ZAVRXX
                     DSECE(IY,IY,IBZ) = DSECE(IY,IY,IBZ) - CSBZ*ZAVRYY
                     DSECE(IZ,IZ,IBZ) = DSECE(IZ,IZ,IBZ) - CSBZ*ZAVRZZ
     &                                - D2*SAVRZ
                  END IF
               END IF
               IF (ISYMAX(3,2) .EQ. IREP) THEN
                  IF (IAX .GT. 0) THEN
                      DSECE(IX,IY,IAX) = DSECE(IX,IY,IAX) + XAVRXY
                      DSECE(IY,IX,IAX) = DSECE(IY,IX,IAX) + XAVRXY
                  END IF
                  IF (IAY .GT. 0) THEN
                      DSECE(IX,IY,IAY) = DSECE(IX,IY,IAY) + YAVRXY
                      DSECE(IY,IX,IAY) = DSECE(IY,IX,IAY) + YAVRXY
                  END IF
                  IF (IAZ .GT. 0) THEN
                      DSECE(IX,IY,IAZ) = DSECE(IX,IY,IAZ) + ZAVRXY
                      DSECE(IY,IX,IAZ) = DSECE(IY,IX,IAZ) + ZAVRXY
                  END IF
                  IF (IBX .GT. 0) THEN
                      DSECE(IX,IY,IBX) = DSECE(IX,IY,IBX) - CSBX*XAVRXY
     *                                 - SAVRY
                      DSECE(IY,IX,IBX) = DSECE(IY,IX,IBX) - CSBX*XAVRXY
     *                                 - SAVRY
                  END IF
                  IF (IBY .GT. 0) THEN
                      DSECE(IX,IY,IBY) = DSECE(IX,IY,IBY) - CSBY*YAVRXY
     *                                 - SAVRX
                      DSECE(IY,IX,IBY) = DSECE(IY,IX,IBY) - CSBY*YAVRXY
     *                                 - SAVRX
                  END IF
                  IF (IBZ .GT. 0) THEN
                      DSECE(IX,IY,IBZ) = DSECE(IX,IY,IBZ) - CSBZ*ZAVRXY
                      DSECE(IY,IX,IBZ) = DSECE(IY,IX,IBZ) - CSBZ*ZAVRXY
                  END IF
               END IF
               IF (ISYMAX(2,2) .EQ. IREP) THEN
                  IF (IAX .GT. 0) THEN
                      DSECE(IX,IZ,IAX) = DSECE(IX,IZ,IAX) + XAVRXZ
                      DSECE(IZ,IX,IAX) = DSECE(IZ,IX,IAX) + XAVRXZ
                  END IF
                  IF (IAY .GT. 0) THEN
                      DSECE(IX,IZ,IAY) = DSECE(IX,IZ,IAY) + YAVRXZ
                      DSECE(IZ,IX,IAY) = DSECE(IZ,IX,IAY) + YAVRXZ
                  END IF
                  IF (IAZ .GT. 0) THEN
                      DSECE(IX,IZ,IAZ) = DSECE(IX,IZ,IAZ) + ZAVRXZ
                      DSECE(IZ,IX,IAZ) = DSECE(IZ,IX,IAZ) + ZAVRXZ
                  END IF
                  IF (IBX .GT. 0) THEN
                      DSECE(IX,IZ,IBX) = DSECE(IX,IZ,IBX) - CSBX*XAVRXZ
     *                                 - SAVRZ
                      DSECE(IZ,IX,IBX) = DSECE(IZ,IX,IBX) - CSBX*XAVRXZ
     *                                 - SAVRZ
                  END IF
                  IF (IBY .GT. 0) THEN
                      DSECE(IX,IZ,IBY) = DSECE(IX,IZ,IBY) - CSBY*YAVRXZ
                      DSECE(IZ,IX,IBY) = DSECE(IZ,IX,IBY) - CSBY*YAVRXZ
                  END IF
                  IF (IBZ .GT. 0) THEN
                      DSECE(IX,IZ,IBZ) = DSECE(IX,IZ,IBZ) - CSBZ*ZAVRXZ
     *                                 - SAVRX
                      DSECE(IZ,IX,IBZ) = DSECE(IZ,IX,IBZ) - CSBZ*ZAVRXZ
     *                                 - SAVRX
                  END IF
               END IF
               IF (ISYMAX(1,2) .EQ. IREP) THEN
                  IF (IAX .GT. 0) THEN
                      DSECE(IY,IZ,IAX) = DSECE(IY,IZ,IAX) + XAVRYZ
                      DSECE(IZ,IY,IAX) = DSECE(IZ,IY,IAX) + XAVRYZ
                  END IF
                  IF (IAY .GT. 0) THEN
                      DSECE(IY,IZ,IAY) = DSECE(IY,IZ,IAY) + YAVRYZ
                      DSECE(IZ,IY,IAY) = DSECE(IZ,IY,IAY) + YAVRYZ
                  END IF
                  IF (IAZ .GT. 0) THEN
                      DSECE(IY,IZ,IAZ) = DSECE(IY,IZ,IAZ) + ZAVRYZ
                      DSECE(IZ,IY,IAZ) = DSECE(IZ,IY,IAZ) + ZAVRYZ
                  END IF
                  IF (IBX .GT. 0) THEN
                      DSECE(IY,IZ,IBX) = DSECE(IY,IZ,IBX) - CSBX*XAVRYZ
                      DSECE(IZ,IY,IBX) = DSECE(IZ,IY,IBX) - CSBX*XAVRYZ
                  END IF
                  IF (IBY .GT. 0) THEN
                      DSECE(IY,IZ,IBY) = DSECE(IY,IZ,IBY) - CSBY*YAVRYZ
     *                                 - SAVRZ
                      DSECE(IZ,IY,IBY) = DSECE(IZ,IY,IBY) - CSBY*YAVRYZ
     *                                 - SAVRZ
                  END IF
                  IF (IBZ .GT. 0) THEN
                     DSECE(IY,IZ,IBZ) = DSECE(IY,IZ,IBZ) - CSBZ*ZAVRYZ
     *                                 - SAVRY
                     DSECE(IZ,IY,IBZ) = DSECE(IZ,IY,IBZ) - CSBZ*ZAVRYZ
     *                                 - SAVRY
                  END IF
               END IF
            END IF
  500       CONTINUE
         END IF
      RETURN
      END
C  /* Deck avepcm */
      SUBROUTINE AVEPCM(PCMINT,MAXCMP,DSHELL,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <maxorb.h>
      PARAMETER (D0 = 0.00 D00, D1 = 1.00 D00, D2 = 2.00D00)
      PARAMETER (DP5 = 0.5 D00)
C
C Used from common blocks:
C  ONECOM : ONECEN, ?
C  CBISOL : NCNTCV
C
#include <onecom.h>
#include <symmet.h>
#include <pcmdef.h>
#include <pcm.h>
#include <pcmlog.h>
#include <taysol.h>
C
      LOGICAL DOAX, DOAY, DOAZ, DOBX, DOBY, DOBZ
      DIMENSION DSHELL(MAXCMP)
      DIMENSION PCMINT(KCKTAB,7,NTSIRR)
#include <ibtfun.h>
C
      NAX = 3*NCENTA - 2
      NAY = 3*NCENTA - 1
      NAZ = 3*NCENTA
      NBX = 3*NCENTB - 2
      NBY = 3*NCENTB - 1
      NBZ = 3*NCENTB
C
C        Multiply densities and integrals
C
      DA0000 = D0
      DA0X00 = D0
      DA0Y00 = D0
      DA0Z00 = D0
      DA000X = D0
      DA000Y = D0
      DA000Z = D0
      DO 300 ICOMP = 1, MAXCMP
         DT0000 = D0
         DT0X00 = D0
         DT0Y00 = D0
         DT0Z00 = D0
         DT000X = D0
         DT000Y = D0
         DT000Z = D0
         DO ISYM = 0, MAXREP
            MTS = ISYM * NTSIRR
            DO 200 ITS = 1, NTSIRR
               QT = QSE(ITS) + QSN(ITS)
               I = ITS + MTS
               DT0000 = DT0000 + PCMINT(ICOMP,1,I)*QT
               IF(.NOT. ONECEN) THEN
                  DT0X00 = DT0X00 + PCMINT(ICOMP,2,I)*QT
                  DT0Y00 = DT0Y00 + PCMINT(ICOMP,3,I)*QT
                  DT0Z00 = DT0Z00 + PCMINT(ICOMP,4,I)*QT
               END IF 
               DT000X = DT000X + PCMINT(ICOMP,5,I)*QT
               DT000Y = DT000Y + PCMINT(ICOMP,6,I)*QT
               DT000Z = DT000Z + PCMINT(ICOMP,7,I)*QT
 200        CONTINUE
         ENDDO
C
         DENS = DSHELL(ICOMP)
C
         DA0000 = DA0000 + DT0000 * DENS
         IF(.NOT. ONECEN) THEN
            DA0X00 = DA0X00 + DT0X00 * DENS
            DA0Y00 = DA0Y00 + DT0Y00 * DENS
            DA0Z00 = DA0Z00 + DT0Z00 * DENS
            DA000X = DA000X + SIGNBX * DT000X * DENS
            DA000Y = DA000Y + SIGNBY * DT000Y * DENS
            DA000Z = DA000Z + SIGNBZ * DT000Z * DENS
         ELSE
            DA000X = DA000X + DT000X * DENS
            DA000Y = DA000Y + DT000Y * DENS
            DA000Z = DA000Z + DT000Z * DENS
         END IF
  300 CONTINUE
      IAX  = IPTCNT(NAX,0,1)
      IAY  = IPTCNT(NAY,0,1)
      IAZ  = IPTCNT(NAZ,0,1)
      IBX  = IPTCNT(NBX,0,1)
      IBY  = IPTCNT(NBY,0,1)
      IBZ  = IPTCNT(NBZ,0,1)
      DOAX = IAX .NE. 0 .AND. .NOT. ONECEN
      DOAY = IAY .NE. 0 .AND. .NOT. ONECEN
      DOAZ = IAZ .NE. 0 .AND. .NOT. ONECEN
      DOBX = IBX .NE. 0
      DOBY = IBY .NE. 0
      DOBZ = IBZ .NE. 0
C
C           Undifferentiated expectation value
C           of SUM(pq) Dpq*((-2 Sum(lm) gl Tlm) Rlm) * (0.5)
C                      DENS*sum(lm) (FCM(lm) * R(lm)) * DP5
C
      IF (IPRINT .GT. 4) THEN
         WRITE (LUPRI,*) 'old ESOLTT  ',ESOLTT
         WRITE (LUPRI,*) 'contribution',DP5*DA0000
         WRITE (LUPRI,*) 'new ESOLTT  ',ESOLTT + DP5*DA0000
      END IF
      ESOLTT = ESOLTT + DA0000 * DP5
C
C           <u(A)/Rlm(C)/u(B)>
C           Consider u(A) derivatives;
C           SUM(pq) Dpq(tpq)a =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)a):
C
      IF (DOAX) GSOLTT(IAX) = GSOLTT(IAX) + DA0X00
      IF (DOAY) GSOLTT(IAY) = GSOLTT(IAY) + DA0Y00
      IF (DOAZ) GSOLTT(IAZ) = GSOLTT(IAZ) + DA0Z00
C
C           <u(A)/Rlm(C)/u(B)>
C           Consider u(B) derivatives;
C           SUM(pq) Dpq(tpq)b =
C             SUM(pq) Dpq*(-2Sum(lm) gl Tlm (Rlm)b):
C
      IF (DOBX) GSOLTT(IBX) = GSOLTT(IBX) + DA000X
      IF (DOBY) GSOLTT(IBY) = GSOLTT(IBY) + DA000Y
      IF (DOBZ) GSOLTT(IBZ) = GSOLTT(IBZ) + DA000Z
      RETURN
      END
