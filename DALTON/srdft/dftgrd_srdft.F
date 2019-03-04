!
!...   Copyright (c) 2014 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2014 (2014), see http://daltonprogram.org"
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


C*****************************************************************************
      SUBROUTINE MAKE_DFTGRID(WORK,LWORK,NTOTAL,IPRINT,TEST)
C*****************************************************************************
C                                                                      
C     Interface routine for grid points                                
C     T. Helgaker                                                      
C                                                                      
C*****************************************************************************

#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "infpar.h"
#include "nuclei.h"
#include "dftcom.h"
#include "symmet.h"
#include "ccom.h"


      PARAMETER (D0=0.0D0)
      LOGICAL TEST
      INTEGER NAT, NUM, ZAN
      COMMON /INFOA/ NAT,NUM,ZAN(MXCENT),C(3,MXCENT)

#include "dftwrk.h"

      DIMENSION WORK(LWORK)

      CALL TIMER('START ',TIMSTR,TIMEND)

#include "memint.h"

C
C     ********************************
C     ***** Set up common blocks *****
C     ********************************
C
C     INFOA
C     =====
CTROND: Modified to handle symmetry dependent atoms as well
      NAT = 0
      NUM = NBASIS
      DO ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         IF(NAMN(ICENT)(1:2).NE.'Gh' .AND.
     &      NAMN(ICENT)(1:2).NE.'IP' .AND.
     &      .NOT. NOORBT(ICENT)) THEN
         IF(IZATOM(ICENT).GT.0) THEN ! omit ghost centers
           DO ISYMOP = 0, MAXOPR         
           IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
             NAT = NAT + 1
             ZAN(NAT) = IZATOM(ICENT)
             C(1,NAT) = PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
             C(2,NAT) = PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
             C(3,NAT) = PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
           ENDIF
           ENDDO
         ENDIF
         ENDIF
      ENDDO
C
C     DFTWRK
C     ======
C
      RNDMAX = 20.0D0
      RADERR = RADINT
      NSCHEME = ANGINT
      LEBMIN = ANGMIN

C
C     *******************************************
C     ***** Prepare arguments for CONSTRUCT *****
C     *******************************************
C
CTROND      IBIG  = 100000
      IBIG  = 1000000
      IVECL = 100
      NNAT  = NAT*(NAT-1)/2
      N2AT = NAT*NAT 

      CALL MEMGET('REAL',KSP   ,4*NUM    ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KX8   ,IBIG     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KY8   ,IBIG     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KZ8   ,IBIG     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KW8   ,IBIG     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPSMU ,IVECL*NAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRJ   ,IVECL*NAT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KACCUM,IVECL    ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXMUJN,IVECL    ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXMUJN,IVECL    ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXMUJ2,IVECL    ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRIJ  ,NNAT     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAIJ  ,NNAT     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRBC  ,N2AT     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRVEC ,3*IVECL  ,WORK,KFREE,LFREE)
C..added for the MOLCAS scheme:
      NDIM=2*NHTYP*NUCIND
      CALL MEMGET('INTE',KNUCO,NDIM  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAA  ,2*NDIM,WORK,KFREE,LFREE)
C
C     Open unit LUQUAD
C     ================
C     (When MPI parallel, then GPOPEN adds a node number
C      as e.g. DALTON.QUAD.n0001 for mynum=1)
C
      LUQUAD = -1
      CALL GPOPEN(LUQUAD,'DALTON.QUAD',
     &     'UNKNOWN','SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUQUAD
C     
C
C     Calculate abscissas and weights
C     ===============================
C
      CALL CONSTRUC(WORK(KSP),WORK(KX8),WORK(KY8),WORK(KZ8),WORK(KW8),
     &              WORK(KPSMU),WORK(KRJ),WORK(KACCUM),
     &              WORK(KXMUJN),WORK(KXMUJ2),
     &              WORK(KRIJ),WORK(KAIJ),IVECL,IBIG,
     &              WORK(KRBC),WORK(KRVEC),
     &              WORK(KNUCO),WORK(KAA),LUQUAD,NTOTAL,IPRINT)
c      IF (TEST) CALL DFTEST(WORK,LWORK,IBIG,LUQUAD,IPRINT)
C
C     Close unit LUQUAD
C     =================
C
      CALL GPCLOSE(LUQUAD,'KEEP')

      CALL TIMER('MAKE_DFTGRID',TIMSTR,TIMEND)

      CALL MEMREL('MAKE_DFTGRID',WORK,1,1,KFREE,LFREE)

      RETURN
      END


C*****************************************************************************
      SUBROUTINE CONSTRUC(WRKSP,X8,Y8,Z8,WT8,PSMU,RJ,ACCUM,XMUIJN,
     &                    XMUIJ2,RIJ,AIJ,IVECL,IBIG,RBC,RVEC,
     &                    NUCORB,AA,LUQUAD,NTOTAL,IPRINT)
C*****************************************************************************
C                                                                      
C      Creates a quadrature grid and writes it to disk                 
C      Original version by C.W.Murray. Rewritten by                    
C      A.M.Lee and D.J.Tozer.                                          
C                                                                      
C      Adapted for Dalton by T. Helgaker                               
C      Adapted for DIRAC by T.Saue                                     
C                                                                      
C*****************************************************************************

#include "implicit.h"
#include "priunit.h"

      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,D3=3.0D0)

      LOGICAL SWITCH
      DIMENSION X8(IBIG),Y8(IBIG),Z8(IBIG),WT8(IBIG)

      CHARACTER SPDCAR*1

#include "dftwrk.h"
#include "mxcent.h"

      INTEGER NAT, NUM, ZAN
      COMMON /INFOA/ NAT,NUM,ZAN(MXCENT),C(3,MXCENT)

#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
#include "ccom.h"

      DIMENSION WRKSP(4*NUM),
     &          ACCUM(IVECL),XMUIJN(IVECL),XMUIJ2(IVECL),
     &          PSMU(IVECL,NAT),RJ(IVECL,NAT)
      DIMENSION NUCORB(NHTYP,2,NUCIND),AA(2,NHTYP,2,NUCIND)
C
C     CARE: Some equivalencing through argument list!
C
      DIMENSION RIJ(NAT*(NAT-1)/2),AIJ(NAT*(NAT-1)/2)
C
C     Bragg contains Bragg-Slater radii for atoms. 
C     Following Becke, the hydrogen radius is 0.35. 
C     For the noble gases, the values have been guessed.
C
      DIMENSION RBC(NAT,NAT)
      DIMENSION BRAGG(0:103)
      DIMENSION RVEC(IVECL,3)
      DIMENSION RADWT(2000), RADNDE(2000)
      DIMENSION ANGWT(9000), XANG(9000), YANG(9000), ZANG(9000)
      DIMENSION XPASC(20)

C Trond Saue:
C     The below data gives atomic radii in Angstroms and stems from table I of 
C     J.C.Slater: "Atomic Radii in Crystals"
C     J.Chem.Phys. 41(1964) 3199-3204
C     Values for elements marked with an asterisk has been
C     guessed/interpolated
C
      DATA BRAGG/ 0.75D0,
C       H      He*   
     &  0.35D0,0.35D0,
C       Li     Be     B      C      N      O      F      Ne*
     &  1.45D0,1.05D0,0.85D0,0.70D0,0.65D0,0.60D0,0.50D0,0.45D0,
C       Na     Mg     Al     Si     P      S      Cl     Ar*
     &  1.80D0,1.50D0,1.25D0,1.10D0,1.00D0,1.00D0,1.00D0,1.00D0,
C       K      Ca     Sc     Ti     V      Cr     Mn     Fe     Co 
     &  2.20D0,1.80D0,1.60D0,1.40D0,1.35D0,1.40D0,1.40D0,1.40D0,1.35D0,
C       Ni     Cu     Zn     Ga     Ge     As     Se     Br     Kr*
     &  1.35D0,1.35D0,1.35D0,1.30D0,1.25D0,1.15D0,1.15D0,1.15D0,1.10D0,
C       Rb     Sr     Y      Zr     Nb     Mo     Tc     Ru     Rh
     &  2.35D0,2.00D0,1.80D0,1.55D0,1.45D0,1.45D0,1.35D0,1.30D0,1.35D0,
C       Pd     Ag     Cd     In     Sn     Sb     Te     I      Xe*
     &  1.40D0,1.60D0,1.55D0,1.55D0,1.45D0,1.45D0,1.40D0,1.40D0,1.40D0,
C       Cs     Ba     La    
     &  2.60D0,2.15D0,1.95D0,
C       Ce     Pr     Nd     Pm     Sm     Eu     Gd
     &  1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,1.80D0,
C       Tb     Dy     Ho     Er     Tm     Yb     Lu
     &  1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,
C       Hf     Ta     W      Re     Os     Ir     Pt     Au     Hg
     &  1.55D0,1.45D0,1.35D0,1.30D0,1.30D0,1.35D0,1.35D0,1.35D0,1.50D0,
C       Tl     Pb*    Bi     Po     At*    Rn*
     &  1.90D0,1.75D0,1.60D0,1.90D0,1.50D0,1.50D0,
C       Fr*    Ra     Ac     
     &  2.15D0,2.15D0,1.95D0,
CTROND rad(U): 1.75 --> 1.37D0
C       Th     Pa     U      Np     Pu     Am     Cm*     
     &  1.80D0,1.80D0,1.37D0,1.75D0,1.75D0,1.75D0,1.75D0,
CTROND       Th     Pa     U      Np     Pu     Am     Cm*     
CTROND     &  1.80D0,1.80D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,
C       Bk*    Cf*    Es*    Fm*    Md*    No*    Lw*
     &  1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0/ 

C
C     USE INPUTTED INFORMATION IN COMMON BLOCK DFTWRK TO DEFINE THE TYPE
C     OF QUADRATURE
C
C     Experimental value of NTRANS
C
      NTRANS=10
C
C     Form Pascals triangle in XPASC for fuzzy Voronoi polyhedra code
C
      ISIGN=-1
      DO 5 I=NTRANS,1,-1
         ISIGN=-ISIGN
         XPASC(I+1) = ISIGN*FACULT(NTRANS)/(FACULT(I)*FACULT(NTRANS-I))
         XPASC(1)   = 1.0D0
5     CONTINUE

      APASC=0.0D0
      DO 6 I=1,NTRANS+1
         XPASC(I)=XPASC(I)/DBLE(2*I-1)
         APASC=APASC+XPASC(I)
6     CONTINUE
      APASC=0.5D0/APASC

      ITEMP=0
      CALL DZERO(RBC,NAT*NAT)
      DO 7 INA=1,NAT
      DO 7 JNA=1,INA-1
          ITEMP=ITEMP+1
          RIJ(ITEMP) = 1.0D0/(sqrt((C(1,INA)-C(1,JNA))**2
     &                             +(C(2,INA)-C(2,JNA))**2
     &                             +(C(3,INA)-C(3,JNA))**2))
          RBC(INA,JNA) = 1.0D0/RIJ(ITEMP)
          RBC(JNA,INA) = 1.0D0/RIJ(ITEMP)
          CHI=BRAGG(ZAN(INA))/BRAGG(ZAN(JNA))
          TEMP=(CHI-1.0D0)/(CHI+1.0D0)
          AIJ(ITEMP)=TEMP/(TEMP*TEMP-1.0D0)
          IF (AIJ(ITEMP).GT.0.5D0) AIJ(ITEMP)=0.5D0
          IF (AIJ(ITEMP).LT.-0.5D0) AIJ(ITEMP)=-0.5D0
7     CONTINUE

C
C     Get information about basis sets
C
      CALL NUCBAS(NUCORB,AA,IPRINT)
C
C     Loop over atoms
C     ===============
C 
      MXOPR = 0 
      NTOTAL = 0
      NALL = 0
      BOHR   = 0.529177249D0
      NATOM = 1
CTROND: define scaling factor for Bragg: set to one !!!!
      TFAC = 1.0D0
      DO 15 IATOM = 1, NUCIND
         IF(NOORBT(IATOM)) GOTO 15
         ! Skip if no orbitals.
         IF(NAMN(IATOM)(1:2).EQ.'Gh') GOTO 15
         IF(NAMN(IATOM)(1:2).EQ.'IP') GOTO 15
         ! Skip if atom name starts with either Gh or IP.
         ! Note: IP is intended for ultra-diffuse orbitals, e.g. 1.d-6, for
         ! ionization potential calculations, for which grid generation will fail.
         ! Anyway, the necessary grid points for points with significant
         ! density will be available from the other grid points.
         ! /hjaaj-2016
         IF(IZATOM(IATOM) .LE. 0) GOTO 15
      DO 10 ISYMOP = 0, MXOPR         
      IF (IAND(ISYMOP,ISTBNU(IATOM)).EQ.0) THEN 
         NDEG  = MULT(ISTBNU(IATOM))
         WRITE(LUPRI,'(A,A4)') '** Atom ',NAMN(IATOM)
         MULA  = ISTBNU(IATOM)
         NDEG  = MULT(MULA)
         IPT   = 0
         INDEX = 0
         NDEX  = 1
         NTHIS = 0
         NTHAT = 0
C
C     Radial quadrature 
C     =================
C     
C     As proposed by 
C     Roland Lindh, Per-Aake Malmqvist and Laura Gagliardi
C

C        Grid spacing....
         WRITE(LUPRI,'(A)') '* Grid spacing'
         H = DUMMY
         DO LL = 1,NHTYP
           L = LL-1
           NBAS=NUCORB(LL,1,IATOM)+NUCORB(LL,2,IATOM)
           IF(NBAS.GT.0) THEN
             HTMP = GRID_DISERR(RADERR,L)
             H = MIN(H,HTMP)
             IF(IPRINT.GE.3) THEN
               WRITE(LUPRI,'(3X,A1,A,F6.3)')
     &        SPDCAR(L),'-orbitals --> ',HTMP
             ENDIF
           ENDIF
         ENDDO
         EPH = EXP(H)
         WRITE(LUPRI,'(A,F6.3)') ' Value chosen: ',H
C...     Inner grid point
         AH = D2*MAX(AA(1,1,1,IATOM),AA(1,1,2,IATOM))
         WRITE(LUPRI,*) 'AH = ',AH
         WRITE(LUPRI,*) 'RADERR = ',RADERR
         RL = ((1.9D0+LOG(RADERR))/D3)-(LOG(AH)/D2)
         RL = EXP(RL)
         WRITE(LUPRI,'(A,1P,E12.5)') '* Inner grid point:',RL
C...     Outer point
         WRITE(LUPRI,'(A)') '* Outer point:'
         RH = D0
         DO LL = 1,NHTYP
           L = LL-1
           AL=DUMMY
           IF(NUCORB(LL,1,IATOM).GT.0) AL=AA(2,LL,1,IATOM)
           IF(NUCORB(LL,2,IATOM).GT.0) AL=MIN(AL,AA(2,LL,2,IATOM))
           IF(AL.LT.DUMMY) THEN
             AL = AL+AL
             RHTMP = GRID_OUTERR(AL,L,RADERR)           
             RH=MAX(RH,RHTMP)
             IF(IPRINT.GE.3) THEN
               WRITE(LUPRI,'(3X,A1,A,F6.3)')
     &            SPDCAR(L),'-orbitals --> ',RHTMP
             ENDIF
           ENDIF
         ENDDO
         WRITE(LUPRI,'(A,F9.3)')  ' Value chosen: ',RH
         GRDC = RL/(EPH-D1)
         WRITE(LUPRI,'(A,1P,E12.5)') ' Constant c:',GRDC
         NR = NINT(LOG(D1+(RH/GRDC))/H)
         WRITE(LUPRI,'(A,I6)') ' Number of points:',NR
         RADNDE(1) = RL
         RADWT(1)  = (RL+GRDC)*RL*RL*H
         DO IR = 2,NR
           RADNDE(IR) = (RADNDE(IR-1)+GRDC)*EPH-GRDC
           RADWT(IR) = (RADNDE(IR)+GRDC)*RADNDE(IR)*RADNDE(IR)*H
         ENDDO
C 
C        loop over radial points
C        =======================
C
         SWITCH = .FALSE.
         RBRAGG = BRAGG(ZAN(NATOM))/(5.0D0*BOHR)
         DO 20 IR = 1, NR
            RWT   = RADWT(IR)
            RNODE = RADNDE(IR)
COLAV            IF (RNODE.GT.RNDMAX) GOTO 20
C 
C           angular abscissas and weights for this radial point
C           ===================================================
C
            IF (.NOT.SWITCH) THEN
               CALL SLEB(XANG,YANG,ZANG,ANGWT,NSCHEME,NANG,NATOM,RNODE,
     &                   RBRAGG)
               IF (RNODE.GT.RBRAGG) SWITCH = .TRUE.
            END IF
            IF (NR*NANG.GT.IBIG) THEN
               WRITE(LUPRI,*) 'Have:',IBIG,' need:',NR*NANG,NANG
               CALL QUIT('Storage error in CONSTRUCT')
            END IF
            IF(IPRINT.GE.3) THEN
              WRITE(LUPRI,'(A,I5,F16.8,I5,I5,5X,L1)') 
     &          '* Grid:',IR,RNODE,NANG,NSCHEME,SWITCH
            ENDIF
C
C           loop over angular points
C           ========================
C
            DO 30 IQ = 1, NANG
               IPT   = IPT + 1
               INDEX = INDEX + 1
               X8(INDEX)  = XANG(IQ)*RNODE+C(1,NATOM)
               Y8(INDEX)  = YANG(IQ)*RNODE+C(2,NATOM)
               Z8(INDEX)  = ZANG(IQ)*RNODE+C(3,NATOM)
               WT8(INDEX) = RWT*ANGWT(IQ)
C
C              Becke partitioning
C              ==================
C
               IF (IPT.EQ.IVECL) THEN
                  CALL BWGHT(RJ,PSMU,RIJ,AIJ,APASC,XPASC,
     &                       X8(NDEX),Y8(NDEX),Z8(NDEX),WT8(NDEX),
     &                       ACCUM,XMUIJN,XMUIJ2,IVECL,IPT,NTRANS,
     &                       NATOM)
                  NDEX   = NDEX  + IPT
CVT               NTOTAL = NTOTAL+ IPT
                  NTHIS  = NTHIS + IPT
                  INDEX  = INDEX - IVECL + IPT
                  IPT    = 0
               END IF
30          CONTINUE
20       CONTINUE
C
C        Becke partitioning
C        ==================
C
         IF (IPT.NE.0) THEN
            JPT = IPT
            CALL BWGHT(RJ,PSMU,RIJ,AIJ,APASC,XPASC,
     &                 X8(NDEX),Y8(NDEX),Z8(NDEX),WT8(NDEX),ACCUM,
     &                 XMUIJN,XMUIJ2,IVECL,IPT,NTRANS,NATOM)
            NDEX   = NDEX   + IPT
CVT         NTOTAL = NTOTAL + IPT
            NTHIS  = NTHIS  + IPT
            INDEX  = INDEX  - JPT + IPT
            IPT = 0
         END IF
         IF(NDEG.GT.0) THEN
C          Scale weights with multiplicity of center
           FAC = NDEG 
           CALL DSCAL(NTHIS,FAC,WT8,1)
         ENDIF
         NTHAT = NTHIS
C
C        Write to file LUQUAD
C        ====================
C
         CALL WRTQUA_OLD(X8,Y8,Z8,WT8,LUQUAD,NTHIS)
         NTOTAL = NTOTAL + NTHIS
         NALL = NALL + NTHAT
         NATOM = NATOM + NDEG
         WRITE(LUPRI,'(A4,I3,3I8)') NAMN(IATOM),NDEG,NTHIS,NR,NANG

      END IF
10    CONTINUE
15    CONTINUE

      NLAST = -1
      CALL WRTQUA_OLD(X8,Y8,Z8,WT8,LUQUAD,NLAST)

      WRITE (LUPRI,'(/,2X,A,I10,A,F5.1,A/)') 
     &     ' Number of grid points in quadrature:',NTOTAL,
     &     ' (',FLOAT(100*NTOTAL)/FLOAT(NALL),'%)'
C
C     Print file
C
      IF (IPRINT .GT. 20) CALL PRIQUA(X8,Y8,Z8,WT8,LUQUAD)

      RETURN
      END


C*****************************************************************************
      SUBROUTINE WRTQUA_OLD(X,Y,Z,WT,LUQUAD,NPOINT)
C*****************************************************************************
C
C     T. Helgaker
C
C*****************************************************************************
#include "implicit.h" 
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      WRITE(LUQUAD) NPOINT
      IF (NPOINT .GT. 0) WRITE(LUQUAD) X,Y,Z,WT
      RETURN
      END


C*****************************************************************************
      SUBROUTINE REAQUA_OLD(X,Y,Z,WT,LUQUAD,NPOINT)
C*****************************************************************************
C
C     T. Helgaker
C
C*****************************************************************************
#include "implicit.h" 
#include "priunit.h"
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      READ(LUQUAD) X,Y,Z,WT
      RETURN
      END


C*****************************************************************************
      SUBROUTINE PRIQUA(X,Y,Z,WT,LUQUAD)
C*****************************************************************************
C
C     T. Helgaker
C
C*****************************************************************************
#include "implicit.h" 
      DIMENSION X(*),Y(*),Z(*),WT(*)
      REWIND LUQUAD
  100 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT .LT. 0) GO TO 200
      CALL PRIQU2(X,Y,Z,WT,LUQUAD,NPOINT)
      GO TO 100
  200 CONTINUE
      RETURN
      END


C*****************************************************************************
      SUBROUTINE PRIQU2(X,Y,Z,WT,LUQUAD,NPOINT)
C*****************************************************************************
C
C     T. Helgaker
C
C*****************************************************************************
#include "implicit.h" 
#include "priunit.h"
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      READ(LUQUAD) X,Y,Z,WT
      DO 100 I = 1, NPOINT
         WRITE (LUPRI,'(2X,4F20.12)') X(I),Y(I),Z(I),WT(I)
  100 CONTINUE
      RETURN
      END


C*****************************************************************************
      integer pure FUNCTION NLEB(NSCHEME)
C*****************************************************************************
      integer, intent(in) :: NSCHEME
      integer :: NTOT
      NTOT = 5810
      
      IF (NSCHEME.LE.131) NTOT=5810
      IF (NSCHEME.LE.125) NTOT=5294
      IF (NSCHEME.LE.119) NTOT=4802
      IF (NSCHEME.LE.113) NTOT=4334
      IF (NSCHEME.LE.107) NTOT=3890
      IF (NSCHEME.LE.101) NTOT=3470
      IF (NSCHEME.LE.95) NTOT=3074
      IF (NSCHEME.LE.89) NTOT=2702
      IF (NSCHEME.LE.83) NTOT=2354
      IF (NSCHEME.LE.77) NTOT=2030
      IF (NSCHEME.LE.71) NTOT=1730
      IF (NSCHEME.LE.65) NTOT=1454
      IF (NSCHEME.LE.59) NTOT=1202
      IF (NSCHEME.LE.53) NTOT=974
      IF (NSCHEME.LE.47) NTOT=770
      IF (NSCHEME.LE.41) NTOT=590
      IF (NSCHEME.LE.35) NTOT=434
      IF (NSCHEME.LE.31) NTOT=350
      IF (NSCHEME.LE.29) NTOT=302
      IF (NSCHEME.LE.27) NTOT=266
      IF (NSCHEME.LE.25) NTOT=230
      IF (NSCHEME.LE.23) NTOT=194
      IF (NSCHEME.LE.21) NTOT=170
      IF (NSCHEME.LE.19) NTOT=146
      IF (NSCHEME.LE.17) NTOT=110
      IF (NSCHEME.LE.15) NTOT=86
      IF (NSCHEME.LE.13) NTOT=74
      IF (NSCHEME.LE.11) NTOT=50
      IF (NSCHEME.LE.9) NTOT=38

      NLEB = NTOT

      RETURN
      END


C*****************************************************************************
      integer pure FUNCTION NLEBSET(NTOT)
C*****************************************************************************
      integer, intent(in) :: NTOT 
      integer :: NSCHEME
      IF (NTOT.GE.5810) THEN
         NSCHEME = 131
      ELSEIF (NTOT.GE.5294) THEN
         NSCHEME = 125
      ELSEIF (NTOT.GE.4802) THEN
         NSCHEME = 119
      ELSEIF (NTOT.GE.4334) THEN
         NSCHEME = 113
      ELSEIF (NTOT.GE.3890) THEN
         NSCHEME = 107
      ELSEIF (NTOT.GE.3470) THEN
         NSCHEME = 101
      ELSEIF (NTOT.GE.3074) THEN
         NSCHEME = 95
      ELSEIF (NTOT.GE.2702) THEN
         NSCHEME = 89
      ELSEIF (NTOT.GE.2354) THEN
         NSCHEME = 83
      ELSEIF (NTOT.GE.1730) THEN
         NSCHEME = 71
      ELSEIF (NTOT.GE.1453) THEN
         NSCHEME = 65
      ELSEIF (NTOT.GE.1201) THEN
         NSCHEME = 59
      ELSEIF (NTOT.GE.973) THEN
         NSCHEME = 53
      ELSEIF (NTOT.GE.769) THEN
         NSCHEME=47
      ELSEIF (NTOT.GE.589) THEN
         NSCHEME=41
      ELSEIF (NTOT.GE.433) THEN
         NSCHEME=35
      ELSEIF (NTOT.GE.349) THEN
         NSCHEME=31
      ELSEIF (NTOT.GE.301) THEN
         NSCHEME=29
      ELSEIF (NTOT.GE.265) THEN
         NSCHEME=27
      ELSEIF (NTOT.GE.229) THEN
         NSCHEME=25
      ELSEIF (NTOT.GE.193) THEN
         NSCHEME=23
      ELSEIF (NTOT.GE.169) THEN
         NSCHEME=21
      ELSEIF (NTOT.GE.145) THEN
         NSCHEME=19
      ELSEIF (NTOT.GE.109) THEN
         NSCHEME=17
      ELSEIF (NTOT.GE.85) THEN
         NSCHEME=15
      ELSEIF (NTOT.GE.73) THEN
         NSCHEME=13
      ELSEIF (NTOT.GE.49) THEN
         NSCHEME=11
      ELSE   
         NSCHEME=9
      ENDIF

      NLEBSET=NSCHEME
 
      RETURN
      END


C*****************************************************************************
      SUBROUTINE SLEB(X,Y,Z,WT,NSCHEME,NCOUNT,NATOM,
     &RNODE,RBRAGG)
C*****************************************************************************
#include "implicit.h" 
#include "priunit.h" 
#include "dftcom.h"
#include "pi.h"

      PARAMETER (MAXNUM=1455)
      DIMENSION X(MAXNUM),Y(MAXNUM),Z(MAXNUM),WT(MAXNUM)

C
C     Lebedev schemes programmed by Cristoph van Wüllen
C     and used by kind permission.
C
C     For principal references, see file Lebedev-Laikov.F
C 
      NTOT=NLEB(NSCHEME)

C     Chose the Lebedev grid to use as a function of the radius (rnode)
C     closer to the nuclei scale the number of angular points as
      NOLD=NSCHEME
      IF(NSCHEME.GT.LEBMIN.AND.RNODE.LT.RBRAGG.AND..NOT.NOPRUN) THEN
         IANG=NINT(NTOT*RNODE/RBRAGG)
         NTOT=MIN(NTOT,IANG)
         NSCHEME=NLEBSET(NTOT)
         NTOT=NLEB(NSCHEME)
C
C     If the pruned value of NSCHEME is less than LEBMIN, we set 
C     NSCHEME = LEBMIN, map NSCHEME to NTOT and then we map NTOT
C     back to a Lebedev scheme that exists on the list.
C
         IF(NSCHEME.LT.LEBMIN) THEN
            NSCHEME = LEBMIN
            NTOT = NLEB(NSCHEME)
            NSCHEME = NLEBSET(NTOT)
         END IF 
      ENDIF

      DO K=1,NTOT
         X(K)=0.0D0
         Y(K)=0.0D0
         Z(K)=0.0D0
         WT(K)=0.0D0
      ENDDO

      IF(NSCHEME.GT.131) THEN
         CALL LD5810(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.131.AND.NSCHEME.GT.125)THEN
         CALL LD5810(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.125.AND.NSCHEME.GT.119)THEN
         CALL LD5294(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.119.AND.NSCHEME.GT.113)THEN
         CALL LD4802(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.113.AND.NSCHEME.GT.107)THEN
         CALL LD4334(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.107.AND.NSCHEME.GT.101)THEN
         CALL LD3890(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.101.AND.NSCHEME.GT.95)THEN
         CALL LD3470(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.95.AND.NSCHEME.GT.89)THEN
         CALL LD3074(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.89.AND.NSCHEME.GT.83)THEN
         CALL LD2702(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.83.AND.NSCHEME.GT.77)THEN
         CALL LD2354(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.77.AND.NSCHEME.GT.71)THEN
         CALL LD2030(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.71.AND.NSCHEME.GT.65)THEN
         CALL LD1730(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.65.AND.NSCHEME.GT.59) THEN
         CALL LD1454(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.59.AND.NSCHEME.GT.53) THEN
         CALL LD1202(X,Y,Z,WT,NCOUNT) 
      ELSEIF(NSCHEME.LE.53.AND.NSCHEME.GT.47) THEN
         CALL LD0974(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.47.AND.NSCHEME.GT.41) THEN
         CALL LD0770(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.41.AND.NSCHEME.GT.35) THEN
         CALL LD0590(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.35.AND.NSCHEME.GT.31) THEN
         CALL LD0434(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.31.AND.NSCHEME.GT.29) THEN
        CALL LD0350(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.29.AND.NSCHEME.GT.27) THEN
        CALL LD0302(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.27.AND.NSCHEME.GT.25) THEN
        CALL LD0266(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.25.AND.NSCHEME.GT.23) THEN
        CALL LD0230(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.23.AND.NSCHEME.GT.21) THEN
        CALL LD0194(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.21.AND.NSCHEME.GT.19) THEN
        CALL LD0170(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.19.AND.NSCHEME.GT.17) THEN
        CALL LD0146(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.17.AND.NSCHEME.GT.15) THEN
        CALL LD0110(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.15.AND.NSCHEME.GT.13) THEN
        CALL LD0086(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.13.AND.NSCHEME.GT.11) THEN
        CALL LD0074(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.11.AND.NSCHEME.GT.9) THEN
        CALL LD0050(X,Y,Z,WT,NCOUNT)
      ELSEIF(NSCHEME.LE.9) THEN
        CALL LD0038(X,Y,Z,WT,NCOUNT)
      END IF

      DO J=1,NCOUNT
         WT(J)=WT(J)*4.0D0*PI
      ENDDO

      IF (NTOT.NE.NCOUNT) THEN
         WRITE(LUPRI,*)'Want',NTOT,', but get',ncount
         CALL QUIT('Error 1 in subroutine_SLEB, get Aaron.')
      ENDIF
       
C     Check values of nodes, and weights; calculate absolute errors.
C     WTERROR returns a value of 0 if any of the weights is zero.
C     EABS is the sum of all absolute deviations of point from
C     the surface of the unit sphere. It is then divided by NCOUNT.

      EABS=0.0D0
      WTERROR=1.0D0
      DO I=1,NCOUNT
         XDEV=abs(1-sqrt(X(I)**2 + Y(I)**2 + Z(I)**2))
         EABS=EABS+XDEV
         IF (XDEV.GT.1D-8) THEN
            WRITE(LUPRI,*) 'NODE ERROR(1) AT I',I,' deviation ',xdev
            CALL QUIT('Error 2 in subroutine_SLEB, get Aaron.')
         ENDIF
         IF (abs(WT(I)).LT.1D-5) THEN
            WTERROR=0
            WRITE(LUPRI,*) 'BWGHT ERROR(1) AT I',I,' weight ',wt(I)
            CALL QUIT('Error 3 in subroutine_SLEB, get Aaron.')
         ENDIF
      ENDDO
      EABS=EABS/NCOUNT
      IF (EABS.LT.1D-15) EABS=1.0D-15

C     Check values of nodes, and weights; calculate absolute errors.
C     WTERROR returns a value of 0 if any of the weights is zero.
C     EABS is the sum of all absolute deviations of point from
C     the surface of the unit sphere. It is then divided by NCOUNT.

      EABS=0.0D0
      WTERROR=0.0D0
      DO I=1,NCOUNT
         XDEV=abs(1.0D0 - sqrt(X(I)**2 + Y(I)**2 + Z(I)**2))
         EABS=EABS+XDEV
         WTERROR=WTERROR+WT(I)
         IF (XDEV.GT.1D-8) THEN
            WRITE(LUPRI,*) 'NODE ERROR(2) AT I',I,' deviation ',xdev
            CALL QUIT('NODE ERROR')
         ENDIF
      ENDDO
      XDEV=abs(WTERROR/4.0D0/PI-1.0D0)
      IF (XDEV.GT.1D-9) THEN
         WRITE(LUPRI,*) 'SUM OF WEIGHTS NOT EQUAL TO 1',XDEV
         CALL QUIT('SUM OF WEIGHTS IN LEBEDEV NOT EQUAL TO 1')
      ENDIF
      EABS=EABS/NCOUNT
      IF (EABS.LT.1D-15) EABS=1.0D-15
      
      NSCHEME=NOLD

      RETURN
      END


C*****************************************************************************
      SUBROUTINE BWGHT(RJ,PSMU,RIJ,AIJ,APASC,XPASC,X,Y,Z,WT,
     &                 ACCUM,XMUIJN,XMUIJ2,IVECL,IPT,NTRANS,NATOM)
C*****************************************************************************
C                                                                      
C      Written by C. W. Murray                                         
C      BWGHT calculates the weights associated with                    
C      the Becke partitioning amongst the atoms.                       
C                                                                      
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"

      PARAMETER (MAXBFN=1000)
      DIMENSION RJ(IVECL,NAT),
     &          PSMU(IVECL,NAT),
     &          RIJ(NAT*(NAT-1)/2),
     &          AIJ(NAT*(NAT-1)/2),
     &          X(IVECL),Y(IVECL),Z(IVECL),WT(IVECL),
     &          XMUIJN(IVECL),XMUIJ2(IVECL),ACCUM(IVECL),XPASC(20)

      INTEGER NAT, NUM, ZAN
      COMMON /INFOA/ NAT,NUM,ZAN(MXCENT),C(3,MXCENT)
C
C     WORK OUT THE WEIGHT FUNCTION BY BECKE PARTITIONING
C
      DO 100 INA=1,NAT
      DO 100 M=1,IPT
         RJ(M,INA)=sqrt((C(1,INA)-X(M))**2
     &                  +(C(2,INA)-Y(M))**2
     &                  +(C(3,INA)-Z(M))**2)
         PSMU(M,INA)=1.0D0
  100 CONTINUE
C
      ITEMP=0
      DO 200 INA=1,NAT
      DO 200 JNA=1,INA-1
         ITEMP=ITEMP+1
         DO 210 M=1,IPT
            XMU=(RJ(M,INA)-RJ(M,JNA))*RIJ(ITEMP)
            XMUIJ=XMU+AIJ(ITEMP)*(1-XMU*XMU)
            XMUIJ2(M)=XMUIJ*XMUIJ
            XMUIJN(M)=XMUIJ
            ACCUM(M)=0.0D0
  210    CONTINUE
         DO 220 I=1,NTRANS+1
         DO 220 M=1,IPT
            ACCUM(M)=ACCUM(M)+XPASC(I)*XMUIJN(M)
            XMUIJN(M)=XMUIJN(M)*XMUIJ2(M)
  220    CONTINUE
         DO 230 M=1,IPT
            PSMU(M,INA)=PSMU(M,INA)*(0.5D0-APASC*ACCUM(M))
            PSMU(M,JNA)=PSMU(M,JNA)*(0.5D0+APASC*ACCUM(M))
  230    CONTINUE
  200 CONTINUE
C
      CALL DZERO(ACCUM,IPT)
      DO 300 INA=1,NAT
      DO 300 M=1,IPT
         ACCUM(M)=ACCUM(M)+PSMU(M,INA)
  300 CONTINUE
C
C     Contract ACCUM into PSMU
C
      DO 400 INA=1,NAT
      DO 400 M=1,IPT
         PSMU(M,INA)=PSMU(M,INA)/ACCUM(M)
  400 CONTINUE
C
      DO 500 M=1,IPT
         WT(M)=WT(M)*PSMU(M,NATOM)
  500 CONTINUE

      RETURN
      END


C*****************************************************************************
      SUBROUTINE CMPRSQ(X,Y,Z,WT,NPOINT,KPOINT)                  
C*****************************************************************************
C                                                                      
C     T. Helgaker                                                      
C                                                                      
C*****************************************************************************
#include "implicit.h"

      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)

      J = 0
      DO I = 1, NPOINT
         NDEG = MLTPNT(X(I),Y(I),Z(I))
         IF (NDEG.GT.0) THEN
            J = J + 1 
            X(J) = X(I)
            Y(J) = Y(I)
            Z(J) = Z(I)
            WT(J) = NDEG*WT(I) 
         END IF
      END DO
      KPOINT = J

      RETURN
      END


C*****************************************************************************
      FUNCTION MLTPNT(PX,PY,PZ)
C*****************************************************************************
C
C     Symmetry multiplicity of a point in space
C
C     T. Helgaker  Feb 01
C
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
      PARAMETER (A0 = 1.0D-8)
      LOGICAL LBTAX
#include "symmet.h"
      LBTAX(I,J) = IAND(2**(I-1),ISYMAX(J,1)) .GT. 0
C
      IF (MAXREP.EQ.0) NSYMOP = 0 
      IF (MAXREP.EQ.1) NSYMOP = 1 
      IF (MAXREP.EQ.3) NSYMOP = 2 
      IF (MAXREP.EQ.7) NSYMOP = 3 
      ISTAB = 0
      DO I = 1, NSYMOP
         IF(LBTAX(I,1) .AND. ABS(PX).GT.A0) GOTO 100
         IF(LBTAX(I,2) .AND. ABS(PY).GT.A0) GOTO 100
         IF(LBTAX(I,3) .AND. ABS(PZ).GT.A0) GOTO 100
            ISTAB = ISTAB + 2**(I-1)
 100     CONTINUE 
      END DO
      MLTPNT = MULT(ISTAB)
      DO I = 1, NSYMOP
         IF(LBTAX(I,1) .AND. PX.LT.-A0) MLTPNT = 0
         IF(LBTAX(I,2) .AND. PY.LT.-A0) MLTPNT = 0
         IF(LBTAX(I,3) .AND. PZ.LT.-A0) MLTPNT = 0
      END DO

      RETURN
      END
