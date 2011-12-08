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
C
C FILE : hergroup.F
C
C Purpose: set up point group symmetry information
C
#include <single.h>
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck syminp */
      SUBROUTINE SYMINP(NSYMOP,KASYM,IFXYZ,CLASS)
C*****************************************************************************
C
C     This subroutine sets up point group symmetry and looks at
C     the behaviour of principal axes and rotations under point group
C     symmetry
C
C     tsaue - 940825 - major revision: tagged on SYMGRP + polish
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "symmet.h"
#include "ccom.h"
#include "cbirea.h"
#include "chrxyz.h"
      CHARACTER*1 KASYM(3,3)
      CHARACTER*(*) CLASS
      DIMENSION IFXYZ(3),IS(0:7),IGEN(3)
      DATA IS/0,1,1,2,1,2,2,3/
#include "ibtfun.h"
C
C     Initialization
C     ==============
C
C     PT is parity of a bitstring:
C       1 for an even number of ones: 000,011,110,101
C      -1 for an odd  number of ones: 001,010,100,111
C
      PT(0) =  1.0D0
      PT(1) = -1.0D0
      PT(2) = -1.0D0
      PT(3) =  1.0D0
      PT(4) = -1.0D0
      PT(5) =  1.0D0
      PT(6) =  1.0D0
      PT(7) = -1.0D0
      DO 5 I = 1,3
         IFXYZ(I)    = 0
         ISYMAX(I,1) = 0
         IGEN(I)     = 0
    5 CONTINUE
C
C     Determine:
C     IGEN(I)   - basic operations
C     ISYMAX(I) - behavior of principal axes under basic operations
C     ===============================================================
C
      IAXIS = 0
      MAXREP = 2**NSYMOP - 1
      IF (NSYMOP.GT.0.AND.IPREAD .GT. 5) THEN
         CALL HEADER('Symmetry Operations',1)
         WRITE (LUPRI,'(A,I2)') '  Symmetry operations:',NSYMOP
      END IF
      DO 100 J = 1,NSYMOP
        DO 110 I = 1,3
        IF (KASYM(I,J).NE.' ') THEN
          K = ICHAR(KASYM(I,J)) - ICHAR('W')
          IGEN(J)     = IGEN(J)     + 2**(K-1)
          ISYMAX(K,1) = ISYMAX(K,1) + 2**(J-1)
        END IF
  110   CONTINUE
        IAXIS = IBTOR(IAXIS,IGEN(J))
  100 CONTINUE
C
C     Determine IFXYZ
C     ===============
C     Do we really need it ????
C
      IND = 0
      DO 140 I = 1,NSYMOP
        IND = IBTOR(IND,IGEN(I))
  140 CONTINUE
      DO 141 I = 1,3
        IFXYZ(I) = IBTAND(IBTSHR(IND,(I-1)),1)
  141 CONTINUE
C
C     Determine:
C     ISYMAX(I,2) - behaviour of principal rotations under basic operations
C     =====================================================================
C
      ISYMAX(1,2) = IBTXOR(ISYMAX(2,1),ISYMAX(3,1))
      ISYMAX(2,2) = IBTXOR(ISYMAX(3,1),ISYMAX(1,1))
      ISYMAX(3,2) = IBTXOR(ISYMAX(1,1),ISYMAX(2,1))
C
C     Determine:
C     IPTAX   - coordinate axis: pointer analogous to IPTSYM
C     NAXREP  - number of coordinate axis in each symmetry
C     ========================================================
C
      CALL IZERO(IPTAX,6)
      CALL IZERO(IPTXYZ(1,0,1),48)
      DO 200 ITYPE = 1,2
        IPTAXI = 0
        DO 205 IREP = 0, MAXREP
           NAXIS = 0
           DO 210 ICOOR = 1, 3
              IF (IBTXOR(IREP,ISYMAX(ICOOR,ITYPE)) .EQ. 0) THEN
                 NAXIS  = NAXIS + 1
                 IPTAXI = IPTAXI + 1
                 IPTAX(ICOOR,ITYPE) = IPTAXI
                 IPTXYZ(NAXIS,IREP,ITYPE) = ICOOR
              END IF
  210      CONTINUE
           NAXREP(IREP,ITYPE) = NAXIS
  205   CONTINUE
  200 CONTINUE
      IF (IPREAD .GT. 5) THEN
         WRITE (LUPRI,'(A,3I5)') '  IPTAX(*,1)  ',
     &                           (IPTAX(I,1),I=1,3)
         WRITE (LUPRI,'(A,3I5)') '  IPTAX(*,2)  ',
     &                           (IPTAX(I,2),I=1,3)
         WRITE (LUPRI,'(A,8I5)') '  NAXREP(*,1) ',
     &                           (NAXREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(A,8I5)') '  NAXREP(*,2) ',
     &                           (NAXREP(I,2),I=0,MAXREP)
      END IF
C
C     Determine group & properties
C
      CALL SYMGRP(IGEN,NSYMOP,CLASS)
C
C     For DIRAC: double group symmetry
C
      IF (DIRAC) CALL DBLGRP
C
C     Determine MULT(I) and FMULT(I) - multiplicity of center
C     =======================================================
C
      DO 101 I = 0,7
         MULT(I)  = 2**MAX(0,NSYMOP-IS(I))
         FMULT(I) = MULT(I)
  101 CONTINUE
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck symgrp */
      SUBROUTINE SYMGRP(IGEN,NSYMOP,CLASS)
C*****************************************************************************
C
C     Given the group generators, this subroutine will identify
C     Abelian subgroup and set up group multiplication table,
C     character table and direct product table. Irreps are identified.
C
C     tsaue - august 1994
C     Sep 24 1996 - tsaue : Included GROUPS,JSOP and IPAR in COMMON SYMMET.
C                           Modified activation of LSYMOP
C                           GROUPS and SYMOP are initialized in BLOCK DATA
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "cbirea.h"
#ifdef PRG_DIRAC
      PARAMETER (LUCME = -1)
#else
#include <inftap.h>
#endif
      LOGICAL LSYMOP(0:7)
      DIMENSION IGEN(3),IROTS(3),IREFL(3),JPAR(0:7)
      CHARACTER*(*) CLASS
      DATA IROTS/3,5,6/
      DATA IREFL/4,2,1/
      DATA JPAR/ 1,-1,-1, 1,-1, 1, 1,-1/
#include "ibtfun.h"

#ifdef MOD_ECP
#include "argoscom.h"
#endif

      DO 5 I = 0,7
        LSYMOP(I) = .FALSE.
    5 CONTINUE
C
C     Activate all symmetry operations of the group
C     =============================================
C
      LSYMOP(0) = .TRUE.
      JSOP(0)   = 0
      IPAR(0)   = 1
      DO 10 I   = 1,MAXREP
        I0 = IBTAND(1,I)          *IGEN(1)
        I1 = IBTAND(1,IBTSHR(I,1))*IGEN(2)
        I2 = IBTAND(1,IBTSHR(I,2))*IGEN(3)
        IND = IBTXOR(IBTXOR(I0,I1),I2)
        LSYMOP(IND) = .TRUE.
        IPAR(I)     = JPAR(IND)
   10 CONTINUE
C
C     List group operations in preferred order,
C       that is construct pointer LSOP
C     =========================================
C
C     Identity
C
      IND = 0
      JSOP(IND) = 0
C
C     Rotations
C
      NROTS = 0
      DO 40 I = 1,3
         IF(LSYMOP(IROTS(I))) THEN
            IND         = IND + 1
            JSOP(IND) = IROTS(I)
            NROTS       = NROTS + 1
         ENDIF
   40 CONTINUE
C
C     Inversion
C
      NINVC = 0
      IF(LSYMOP(7)) THEN
        IND         = IND + 1
        JSOP(IND) = 7
        NINVC       = 1
      ENDIF
C
C     Reflections
C
      NREFL = 0
      DO 50 I = 1,3
      IF(LSYMOP(IREFL(I))) THEN
        IND         = IND + 1
        JSOP(IND) = IREFL(I)
        NREFL       = NREFL + 1
      ENDIF
   50 CONTINUE

      IF(IND.NE.MAXREP) CALL QUIT('SYMGRP: IND.NE.MAXREP!')
C
C     Classify group
C     ==============
C     tsaue - Here I have devised a highly empirical formula, but
C             it works !!!
C
      IGROUP = MIN(7,NINT((4*NROTS+8*NINVC+6*NREFL)/3.0))
      GROUP  = GROUPS(IGROUP)

c     export it to ARGOS
c     ==================
#ifdef MOD_ECP
      ndptag = ndptags(IGROUP)
#endif
C
C     Generate character table
C     ========================
C
      DO 60 I = 0,MAXREP
        IXVAL(0,I) = 1
        DO 70 J = 1,NSYMOP
          IXVAL(IGEN(J),I) = NINT(PT(IBTAND(IBTSHR(I,(J-1)),1)))
          DO 80 K = 1,(J-1)
            IND      = IBTXOR(IGEN(J),IGEN(K))
            IXVAL(IND,I)  = IXVAL(IGEN(J),I)*IXVAL(IGEN(K),I)
            DO 90 L = 1,(K-1)
              IXVAL(IBTXOR(IND,IGEN(L)),I)
     &           = IXVAL(IND,I)*IXVAL(IGEN(L),I)
   90       CONTINUE
   80     CONTINUE
   70   CONTINUE
   60 CONTINUE
C
C     Classify irrep
C     ==============
C
      DO 100 I = 0,MAXREP
        REP(I) = 'A  '
        IPOS = 2
C
C       Rotational symmetry
C
        IF(NROTS.EQ.3) THEN
          IND = (1-IXVAL(JSOP(1),I))+(1-IXVAL(JSOP(2),I))/2
          IF(IND.NE.0) THEN
            REP(I)(1:1) = 'B'
            REP(I)(2:2) = CHAR(ICHAR('0')+IND)
            IPOS = 3
          ENDIF
        ELSEIF(NROTS.EQ.1) THEN
          IF(IXVAL(JSOP(1),I).EQ.-1) REP(I)(1:1) = 'B'
          IF(NREFL.EQ.2) THEN
            IF(IBTAND(IBTSHR(JSOP(1),1),1).EQ.1) THEN
              IND = 2
            ELSE
              IND = 3
            ENDIF
            IF(IXVAL(JSOP(IND),I).EQ.1) THEN
              REP(I)(2:2) = '1'
            ELSE
              REP(I)(2:2) = '2'
            ENDIF
          ENDIF
        ELSEIF(NREFL.EQ.1) THEN
C
C       Mirror symmetry
C
          IF(IXVAL(JSOP(1),I).EQ.1) THEN
            REP(I)(2:2) = ''''
          ELSEIF(IXVAL(JSOP(1),I).EQ.-1) THEN
            REP(I)(2:2) = '"'
          ENDIF
        ENDIF
C
C       Inversion symmetry
C
        IF(NINVC.EQ.1) THEN
          IND = NROTS+1
          IF(IXVAL(JSOP(IND),I).EQ.1) THEN
            REP(I)(IPOS:IPOS) = 'g'
          ELSE
            REP(I)(IPOS:IPOS) = 'u'
          ENDIF
        ENDIF
#ifdef MOD_ECP
c           
c       For ARGOS ECP
c       ============½=
        itypag(i) = REP(I)
#endif            
  100 CONTINUE
C
C     Output section
C     ==============
C
      IF(IPREAD.GT.0) THEN
        CALL HEADER('SYMGRP: Point group information',-1)
C
C       Group name
C
        IF (CLASS(1:3) .EQ. 'N/A') THEN
           WRITE(LUPRI,'(A,A3)') 'Point group: ',GROUP
           IF (LUCME.GT.0)
     &     WRITE(LUCME,'(A,A3)') 'Point group: ',GROUP
        ELSE
           WRITE(LUPRI,'(A,A)')  'Full point group is: ',CLASS
           WRITE(LUPRI,'(A,A3)') 'Represented as:      ',GROUP
           IF (LUCME.GT.0) THEN
              WRITE(LUCME,'(A,A)')  'Full point group is: ',CLASS
              WRITE(LUCME,'(A,A3)') 'Represented as:      ',GROUP
           END IF
        END IF
C
C       Group generators
C
        IF (NSYMOP.GT.0) THEN
          WRITE(LUPRI,'(/3X,A/)')
     &        '* The point group was generated by:'
          DO I = 1,NSYMOP
            IF    (SYMOP(IGEN(I))(1:1).EQ.'C') THEN
              WRITE(LUPRI,'(6X,3A)')
     &        'Rotation about the ',SYMOP(IGEN(I))(3:3),'-axis'
            ELSEIF(SYMOP(IGEN(I))(1:1).EQ.'O') THEN
              WRITE(LUPRI,'(6X,3A)')
     &        'Reflection in the ',SYMOP(IGEN(I))(2:3),'-plane'
            ELSE
              WRITE(LUPRI,'(6X,A)') 'Inversion center'
            END IF
          END DO
C
C         Group multiplication table
C
          WRITE(LUPRI,'(/3X,A/)') '* Group multiplication table'
          WRITE(LUPRI,'(8X,A1,8(1X,A3,1X))')'|',(SYMOP(JSOP(I)),
     &        I = 0,MAXREP)
          WRITE(LUPRI,'(3X,A6,8A5)') '-----+',('-----',I = 0,MAXREP)
          DO I = 0,MAXREP
            WRITE(LUPRI,'(4X,A3,1X,A1,8(1X,A3,1X))')
     &        SYMOP(JSOP(I)),'|',
     &        (SYMOP(IBTXOR(JSOP(I),JSOP(J))),J = 0,MAXREP)
          END DO
C
C         Character table
C
          WRITE(LUPRI,'(/3X,A/)') '* Character table'
          WRITE(LUPRI,'(8X,A1,8(1X,A3,1X))') '|',(SYMOP(JSOP(J)),
     &               J = 0,MAXREP)
          WRITE(LUPRI,'(3X,A6,8A5)') '-----+',('-----',I = 0,MAXREP)
          DO I = 0,MAXREP
            WRITE(LUPRI,'(4X,A3,1X,A1,8(1X,I3,1X))')
     &      REP(I),'|',(IXVAL(JSOP(J),I),J=0,MAXREP)
          END DO
C
C         Direct product table
C
          WRITE(LUPRI,'(/3X,A/)') '* Direct product table'
          WRITE(LUPRI,'(8X,A1,8(1X,A3,1X))')'|',(REP(I),I = 0,MAXREP)
          WRITE(LUPRI,'(3X,A6,8A5)') '-----+',('-----',I = 0,MAXREP)
          DO I = 0,MAXREP
            WRITE(LUPRI,'(3X,1X,A3,1X,A1,8(1X,A3,1X))')
     &        REP(I),'|',(REP(IBTXOR(I,J)),J = 0,MAXREP)
          END DO
        ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dblgrp */
      SUBROUTINE DBLGRP
C*****************************************************************************
C
C     This routine will analyze the fermion irreps in the molecular symmetry
C     and, if possible, set up a transformation within AO - basis to a
C     symmetry-adapted spinor basis. details are given below.
C
C     Written by T.Saue - october 1994 - odense
C     Sep 24 1996 - tsaue: Updated version
C     Febr.2007, M.Ilias - added operator type 20
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(DM1 = -1.0D0,D1 = 1.0D0)
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
C
#include "cbirea.h"
#include "symmet.h"
#include "pgroup.h"
#include "dgroup.h"
#include "qalgebra.h"
C
      LOGICAL LBUF(8)
      EXTERNAL BDQALGEBRA
C
#include "ibtfun.h"
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
      CALL IZERO(IQTOPQ,32)
      CALL IZERO(IPQTOQ,32)
C
C     Determine number of boson and fermion irreps
C     ============================================
C
      NBSYM = MAXREP + 1
      NFSYM = 1
      IF(NINVC.EQ.1) NFSYM = 2
      IF(NFSYM.EQ.1) THEN
        FREP(1) = 'E1 '
      ELSE
        FREP(1) = 'E1g'
        FREP(2) = 'E1u'
      ENDIF
C
C     Determine NZ - parameter
C     ========================
C       NZ = 4  quaternionic group        no    totally symmetric rotations
C       NZ = 2  complex group             one   totally symmetric rotation
C       NZ = 1  real group                three totally symmetric rotations
C
      NZ = 1
      DO 30 I = 1,3
        IF(ISYMAX(I,2).EQ.0) THEN
          NZ   = NZ + 1
          ITOT = I
        ENDIF
   30 CONTINUE
C
C     Align rotations and quaternion irreps,
C     for complex groups this is important....
C
C     The three quaternion units are equivalent and can
C     thus be interchanged. Normally one makes the mapping
C       (x,y,z) --> (k,j,i)
C     and this is the default, including complex groups with
C     a totally symmetric z-rotation
C
C     Variables IXQ, IYQ and IZQ gives the mapping from coordinates
C     to quaternion units 1,2, and 3, whereas IQ1, IQ2 and IQ3 gives
C     the reverse mapping.
C
      IQ1 = 3
      IQ2 = 2
      IQ3 = 1
      IZQ = 2
      IYQ = 3
      IXQ = 4
C.....complex groups
      IF(NZ.EQ.2) THEN
        IF(ITOT.EQ.1) THEN
C.........totally symmetric x-rotation
C         (x,y,z) --> (i,k,j)
          IQ1 = 1
          IQ2 = 3
          IQ3 = 2
          IXQ = 2
          IZQ = 3
          IYQ = 4
        ELSEIF(ITOT.EQ.2) THEN
C.........totally symmetric y-rotation
C         (x,y,z) --> (j,i,k)
          IQ1 = 2
          IQ2 = 1
          IQ3 = 3
          IYQ = 2
          IXQ = 3
          IZQ = 4
        ENDIF
      ENDIF
C
C
C     Set up direct product table for fermion ircops
C     ==============================================
C
C     IF there is a totally symmetric irrep, it should be the
C     symmetry of Im(p_1,q_1)......
C
C     Re(p_1,q_1)
      IFDIRP(1,1) = 0
C     Im(p_1,q_1)
      IFDIRP(2,1) = ISYMAX(IQ1,2)
C     Re(p_1,Q_1)
      IFDIRP(3,1) = ISYMAX(IQ2,2)
C     Im(p_1,Q_1)
      IFDIRP(4,1) = ISYMAX(IQ3,2)
C     Re(p_1,q_2)
      IFDIRP(1,2) = IXYZ
C     Im(p_1,q_2)
      IFDIRP(2,2) = ISYMAX(IQ1,1)
C     Re(p_1,Q_2)
      IFDIRP(3,2) = ISYMAX(IQ2,1)
C     Im(p_1,Q_2)
      IFDIRP(4,2) = ISYMAX(IQ3,1)
C
C
C     Determine distribution of boson irreps in spinors
C     =================================================
C
      DO 10 IM = 1,4
        JBTOF(IFDIRP(IM,1),1) = 1
        JSPINR(IM,1,1) = IFDIRP(IM,1)
        JBTOF(IFDIRP(IM,2),2) = 1
        JSPINR(IM,2,1) = IFDIRP(IM,2)
   10 CONTINUE
      IF(NFSYM.EQ.2) THEN
        DO 20 IM = 1,4
          JBTOF(IFDIRP(IM,2),1) = 2
          JSPINR(IM,1,2) = IFDIRP(IM,2)
          JBTOF(IFDIRP(IM,1),2) = 2
          JSPINR(IM,2,2) = IFDIRP(IM,1)
   20   CONTINUE
      ENDIF
C
C     Determine boson symmetries connected by a given
C     fermion irrep
C
      IF(NFSYM.EQ.1) THEN
        DO ISYM = 1,NBSYM
          JFSYM(ISYM,1) = ISYM
          JFSYM(ISYM,2) = ISYM
        ENDDO
      ELSEIF(NFSYM.EQ.2) THEN
        ISYMA = 0
        ISYMB = 0
        DO ISYM = 1,NBSYM
          IF(JBTOF(ISYM-1,1).EQ.1) THEN
            ISYMA = ISYMA + 1
            JFSYM(ISYMA,1) = ISYM
          ELSE
            ISYMB = ISYMB + 1
            JFSYM(ISYMB,2) = ISYM
          ENDIF
        ENDDO
      ENDIF
C
C     Matrix packing
C     ==============
C
      JMROI(1) = 1
      IQDEF(1) = 1
      IQDEF(2) = 2
      IQDEF(3) = 3
      IQDEF(4) = 4
C
C     Quaternion groups
C
      IF(NZ.EQ.4) THEN
        JMROI (2) = 1
        JMROI (3) = 1
        JMROI (4) = 1
        IF(NFSYM.EQ.1) THEN
C       C1:
          DO I = 1,4
            IPQTOQ(I,   0) = I
          ENDDO
        ELSE
C       Ci:
          DO I = 1,4
            IPQTOQ(I,   0) = I
            IPQTOQ(I,IXYZ) = I
          ENDDO
        ENDIF
C
C     Complex groups:
C
      ELSEIF(NZ.EQ.2) THEN
C
C       Totally symmetric Z-rotation
C
          JMROI (2) = 1
          JMROI (3) = 4
          JMROI (4) = 4
C         ..non-totally symmetric rotation
          IRNT           = ISYMAX(IQ2,2)
          IPQTOQ(1,0   ) = 1
          IPQTOQ(2,0   ) = 2
          IPQTOQ(1,IRNT) = 4
          IPQTOQ(2,IRNT) = 3
        IF(NFSYM.EQ.2) THEN
          IRNU = IBTXOR(IRNT,IXYZ)
          DO IZ = 1,NZ
            IPQTOQ(IZ,IXYZ) = IPQTOQ(IZ,0)
            IPQTOQ(IZ,IRNU) = IPQTOQ(IZ,IRNT)
          ENDDO
        ENDIF
      ELSE
        JMROI(2)  = 2
        JMROI(3)  = 3
        JMROI(4)  = 4
        IF(NFSYM.EQ.1) THEN
          IPQTOQ(1,0) = 1
          DO I = 2,4
            JMROI(I) = I
            IPQTOQ(1,ISYMAX(5-I,2)) = I
          ENDDO
        ELSE
          IPQTOQ(1,   0) = 1
          IPQTOQ(1,IXYZ) = 1
          DO I = 2,4
            JMROI(I) = I
            IPQTOQ(1,ISYMAX(5-I,2)) = I
            IPQTOQ(1,ISYMAX(5-I,1)) = I
          ENDDO
        ENDIF
      ENDIF
C
C     Assign IQTOPQ
C
      DO IREP = 0,MAXREP
        DO IPQ = 1,NZ
          IQ = IPQTOQ(IPQ,IREP)
          IQTOPQ(IQ,IREP) = IPQ
        ENDDO
      ENDDO
C
C     Assign JQBAS - gives quaternion vector for given
C     component and bosonirrep
C     ================================================
C
      DO 60 IFR = 1,NFSYM
        DO 70 IC = 1,2
          DO 80 IM = 1,4
            JQBAS(JSPINR(IM,IC,IFR),IC) = JMROI(IM)
   80     CONTINUE
   70   CONTINUE
   60 CONTINUE
C
C     Give signs of quaternion phases: 1,i,-j,k
C
      IQPH(1,1) =  1
      IQPH(2,1) =  1
      IQPH(3,1) =  1
      IQPH(4,1) =  1
C
C     ...and their Hermitian conjugates: 1,-i, j,-k
C
      IQPH(1,2) =  1
      IQPH(2,2) = -1
      IQPH(3,2) = -1
      IQPH(4,2) = -1
C
C     ********************************************
C     *****   Define 4-component operators   *****
C     ********************************************
C
C     M = I_4,A_Z,A_Y,A_X, G_5,S_Z,S_Y,S_X
C
C
C     I_4       -->  I_2
C     ==================
C
      JM4REP(0) =  0
      JM4POS(0) =  1
      JM4FAS(0) =  1
C
C     i(alpha_z)   --> i(sigma_x)
C     ===========================
C
      JM4REP(1) =  ISYMAX(3,1)
      JM4POS(1) =  IZQ
      JM4FAS(1) =  1
C
C     i(alpha_y)   --> j(sigma_x)
C     ===========================
C
      JM4REP(2) =  ISYMAX(2,1)
      JM4POS(2) =  IYQ
      JM4FAS(2) =  1
C
C     i(alpha_x)   --> k(sigma_x)
C     ===========================
C
      JM4REP(3) =  ISYMAX(1,1)
      JM4POS(3) =  IXQ
      JM4FAS(3) =  1
C
C     | 0 I |       | 0 1 |
C     | I 0 |  -->  | 1 0 |
C     =====================
C
      JM4REP(4) =  IXYZ
      JM4POS(4) =  1
      JM4FAS(4) =  1
C
C     i(sigma_z)   --> i(I_2)
C     =======================
C
      JM4REP(5) =  ISYMAX(3,2)
      JM4POS(5) =  IZQ
      JM4FAS(5) =  1
C
C     i(sigma_y)   --> j(I_2)
C     =======================
C
      JM4REP(6) =  ISYMAX(2,2)
      JM4POS(6) =  IYQ
      JM4FAS(6) =  1
C
C     i(sigma_x)   --> k(I_2)
C     =======================
C
      JM4REP(7) =  ISYMAX(1,2)
      JM4POS(7) =  IXQ
      JM4FAS(7) =  1
C
C     ***************************************
C     ***** Definition of full operator *****
C     ***************************************
C
C
C     1. P             * scalar operator
C     ==================================
C
      MCMP(1)   = 1
      JM4 (1,1) = 0
      JCOM(1,1) = 1
C
C     2. i[alpha_x]P    * x-component of alpha times scalar operator
C     =============================================================
C
      MCMP(2)   =  1
      JM4 (1,2) =  3
      JCOM(1,2) =  1
C
C     3. i[alpha_y]P    * y-component of alpha times scalar operator
C     =============================================================
C
      MCMP(3)   =  1
      JM4 (1,3) =  2
      JCOM(1,3) =  1
C
C     4. i[alpha_z]P    * z-component of alpha times scalar operator
C     =============================================================
C
      MCMP(4)   =  1
      JM4 (1,4) =  1
      JCOM(1,4) =  1
C
C     5. i[Alpha x P]_x * vector product of alpha and vector operator,
C                         x-component
C     ================================================================
C
      MCMP(5)   =  2
      JM4 (1,5) =  2
      JM4 (2,5) =  1
      JCOM(1,5) =  1
      JCOM(2,5) = -1
C
C     6. i[Alpha x P]_y * vector product of alpha and vector operator,
C                         y-component
C     ================================================================
C
      MCMP(6)   =  2
      JM4 (1,6) =  1
      JM4 (2,6) =  3
      JCOM(1,6) =  1
      JCOM(2,6) = -1
C
C     7. i[Alpha x P]_z * vector product of alpha and vector operator,
C                         z-component
C     ================================================================
C
      MCMP(7)   =  2
      JM4 (1,7) =  3
      JM4 (2,7) =  2
      JCOM(1,7) =  1
      JCOM(2,7) = -1
C
C     8. iA.P           * dot-product of alpha and vector operator
C     ===========================================================
C
      MCMP(8)   =  3
      JM4 (1,8) =  3
      JM4 (2,8) =  2
      JM4 (3,8) =  1
      JCOM(1,8) =  1
      JCOM(2,8) =  1
      JCOM(3,8) =  1
C
C     9. gamma5 P       * gamma5 times scalar operator
C     ===========================================================
C
      MCMP(9)   =  1
      JM4 (1,9) =  4
      JCOM(1,9) =  1
C
C     10. i[Sigma_x]P   * x-component of sigma times scalar operator
C     ==============================================================
C
      MCMP(10)  =  1
      JM4 (1,10)=  7
      JCOM(1,10)=  1
C
C     11. i[Sigma_y]P   * y-component of sigma times scalar operator
C     ==============================================================
C
      MCMP(11)  =  1
      JM4 (1,11)=  6
      JCOM(1,11)=  1
C
C     12. i[Sigma_z]P   * z-component of sigma times scalar operator
C     ==============================================================
C
      MCMP(12)  =  1
      JM4 (1,12)=  5
      JCOM(1,12)=  1
C
C     13. i[betaSig_x]P * x-component of beta sigma times scalar operator
C     ==================================================================
C
      MCMP(13)  =  1
      JM4 (1,13)=  7
      JCOM(1,13)=  1
C
C     14. i[betaSig_y]P * y-component of beta sigma times scalar operator
C     ==================================================================
C
      MCMP(14)  =  1
      JM4 (1,14)=  6
      JCOM(1,14)=  1
C
C     15. i[betaSig_z]P * z-component of beta sigma times scalar operator
C     ==================================================================
C
      MCMP(15)  =  1
      JM4 (1,15)=  5
      JCOM(1,15)=  1
C
C     16. i[betaalp_x]P * x-component of beta alpha times scalar operator
C     ==================================================================
C
      MCMP(16)  =  1
      JM4 (1,16)=  3
      JCOM(1,16)=  1
C
C     17. i[betaalp_y]P * y-component of beta alpha times scalar operator
C     ==================================================================
C
      MCMP(17)  =  1
      JM4 (1,17)=  2
      JCOM(1,17)=  1
C
C     18. i[betaalp_z]P * z-component of beta alpha times scalar operator
C     ===================================================================
C
      MCMP(18)  =  1
      JM4 (1,18)=  1
      JCOM(1,18)=  1
C
C
C     19. beta         * scalar operator
C     ==================================
C
      MCMP(19)  = 1
      JM4 (1,19)= 0
      JCOM(1,19)= 1
C       
C         
C     20. iS.P         * dot-product of sigma and vector operator
C     ==================================
C
      MCMP(20)   =  3
      JM4 (1,20) =  3
      JM4 (2,20) =  2
      JM4 (3,20) =  1
      JCOM(1,20) =  1
      JCOM(2,20) =  1
      JCOM(3,20) =  1
C
C     Matrix symmetry
C     ===============
C
C     The matrix of an operator that is symmetric or
C     antisymmetric under time reversal has the following
C     structure:
C
C       A   B          A^{dagger} =   hA
C     -tB* tA*         B^T        = -thB
C
C     The matrix symmetry of matrices A and B
C     can be summarized as follows:
C
C     IH   ITIM   AR  AI  BR  BI
C      1     1     1   2   2   2
C     -1     1     2   1   1   1
C      1    -1     1   2   1   1
C     -1    -1     2   1   2   2
C
C     Only ITIM=1 is used (ITIM=-1 is transferred to IH=-1
C     by extracting an imaginary i), and the matrix symmetry
C     of matrices A and B is stored in IHQMAT(4,IH):
C
C     H-
      IHQMAT(1,-1) = 2
      IHQMAT(2,-1) = 1
      IHQMAT(3,-1) = 1
      IHQMAT(4,-1) = 1
C     H0
      IHQMAT(1, 0) = 0
      IHQMAT(2, 0) = 0
      IHQMAT(3, 0) = 0
      IHQMAT(4, 0) = 0
C     H+
      IHQMAT(1, 1) = 1
      IHQMAT(2, 1) = 2
      IHQMAT(3, 1) = 2
      IHQMAT(4, 1) = 2
C
C     IRQMAT gives the irrep of a given component of
C     quaternion matrix:
C
      DO IZ = 1,4
        DO IBRP = 0,MAXREP
          IRQMAT(IZ,IBRP) = IBTXOR(IFDIRP(IZ,1),IBRP)
        ENDDO
      ENDDO
C
C     Output section
C     ==============
C
      IF(IPREAD.GE.1) THEN
        CALL TITLER('Output from DBLGRP','*',103)
        IF    (NFSYM.EQ.1) THEN
          WRITE(LUPRI,'(3X,A,2X,A3)')
     &      '* One fermion irrep: ',FREP(1)
        ELSEIF(NFSYM.EQ.2) THEN
          WRITE(LUPRI,'(3X,A,2(2X,A3))')
     &      '* Two fermion irreps:',FREP(1),FREP(2)
        ENDIF
        IF    (NZ.EQ.4) THEN
          WRITE(LUPRI,'(3X,A)') '* Quaternionic group. NZ = 4'
        ELSEIF(NZ.EQ.2) THEN
          WRITE(LUPRI,'(3X,A)') '* Complex group. NZ = 2'
        ELSEIF(NZ.EQ.1) THEN
          WRITE(LUPRI,'(3X,A)') '* Real group. NZ = 1'
        ENDIF
        WRITE(LUPRI,'(3X,A)') '* Direct product decomposition:'
        IND = 1
        DO 110 I1 = 1,NFSYM
          DO 120 I2 = I1,NFSYM
            IND = IND + 1
            J = MOD(IND,2)+1
            WRITE(LUPRI,'(10X,2(A3,A),3(A3,A),A3)')
     &        FREP(I2),' x ',FREP(I1),' : ',
     &        REP(IFDIRP(1,J)),' + ',
     &        REP(IFDIRP(2,J)),' + ',
     &        REP(IFDIRP(3,J)),' + ',
     &        REP(IFDIRP(4,J))
  120     CONTINUE
  110   CONTINUE
        CALL HEADER('Spinor structure',-1)
        WRITE(LUPRI,'(/2(3X,A,I2,9X)/)')
     &     ('* Fermion irrep no.:',I,I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &     ('La','|',REP(JSPINR(1,1,I)),
     &            '(',JQBAS(JSPINR(1,1,I),1),')',
     &               REP(JSPINR(2,1,I)),
     &            '(',JQBAS(JSPINR(2,1,I),1),')','|',I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &     ('Sa','|',REP(JSPINR(1,2,I)),
     &            '(',JQBAS(JSPINR(1,2,I),2),')',
     &               REP(JSPINR(2,2,I)),
     &            '(',JQBAS(JSPINR(2,2,I),2),')','|',I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
C        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3),2X,A1,10X))')
     &     ('Lb','|',REP(JSPINR(3,1,I)),
     &            '(',JQBAS(JSPINR(3,1,I),1),')',
     &               REP(JSPINR(4,1,I)),
     &            '(',JQBAS(JSPINR(4,1,I),1),')','|',I=1,NFSYM)
        WRITE(LUPRI,'(2(6X,A2,2X,A1,2(2X,A3,A1,I1,A1),2X,A1,10X))')
     &   ('Sb','|',REP(JSPINR(3,2,I)),
     &            '(',JQBAS(JSPINR(3,2,I),2),')',
     &             REP(JSPINR(4,2,I)),
     &            '(',JQBAS(JSPINR(4,2,I),2),')','|',I=1,NFSYM)
        CALL HEADER('Quaternion symmetries',-1)
        WRITE(LUPRI,'(4X,A3,2X,A4)') 'Rep','T(+)'
        CALL PRSYMB(LUPRI,'-',29,4)
        DO IBRP = 0,NBSYM-1
          WRITE(LUPRI,'(4X,A3,4(2X,A1))')
     &       REP(IBRP),(QUNIT(IPQTOQ(IZ,IBRP)),IZ=1,NZ)
        ENDDO
      ENDIF
      RETURN
C
      END
!  -- end of hergroup.F --
