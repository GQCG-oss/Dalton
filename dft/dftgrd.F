C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
#include <single.h>
C
C  /* Deck dftgrd */
      SUBROUTINE DFTGRD(WORK,LWORK,IPRINT,TEST)
C
C     Interface routine for grid points
C     T. Helgaker 
C
#include <implicit.h>
#include <priunit.h>
C
#include <dummy.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <infpar.h>
#include <nuclei.h>
#include <dftcom.h>
#include <symmet.h>
C
      PARAMETER (MAXCEN=100)
      LOGICAL TEST, SPCDIV
      COMMON /INFOA/ NAT,NUM,ZAN(MAXCEN),C(3,MAXCEN)
      LOGICAL KHNSHM,DFTGRA,GRDWT
      COMMON /DFTCWM/ DFTGRA,GRDWT,RNDMAX,IQUAD,NRHX,NROX,NRSX,NTMAXX,
     &                NPHIX,KHNSHM
      COMMON /DFTOZER/ NRTYPF
      COMMON /DFTAML/ NSCHEME
C
      CHARACTER*16 QUADNAME
      DIMENSION WORK(LWORK)
#include <chrnos.h>
#include <ibtfun.h>
C
      CALL TIMER('START ',TIMSTR,TIMEND)
      QUADNAME = '                '
C
C     ********************************
C     ***** Set up common blocks *****
C     ********************************
C
C     INFOA
C     =====
      NAT = 0
      NUM = NBASIS
      DO ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         ZCNT   = CHARGE(ICENT)
         IF(NAMN(ICENT)(1:2).NE.'Gh') THEN
           DO ISYMOP = 0, MAXOPR         
           IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
             NAT = NAT + 1
             ZAN(NAT) = ZCNT
             C(1,NAT) = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
             C(2,NAT) = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
             C(3,NAT) = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
           ENDIF
           ENDDO
         ENDIF
      ENDDO
C
C     DFTCWM
C     ======
C
      SPCDIV = IQDTYP .GT. 0 .AND. MAXOPR .GT. 0
      spcdiv = .false.
      IQDTYP = ABS(IQDTYP)
      IF (IQDTYP.EQ.1) IQUAD = 2
      IF (IQDTYP.EQ.2) IQUAD = 4
      IF (IQDTYP.EQ.3) IQUAD = 5
      IF (IQDTYP.EQ.4) IQUAD = 1
      IF (IQDTYP.EQ.5) IQUAD = 7
C
C     IQUAD = 2     very coarse
C     IQUAD = 4     coarse
C     IQUAD = 5     medium
C     IQUAD = 1     dense
C     IQUAD = 7     very dense
C
      IF (IQUAD.NE.1) THEN
         NRHX   = 12345
         NROX   = 12345
         NRSX   = 12345
         NTMAXX = 12345
         NPHIX  = 12345
      END IF
      KHNSHM = .TRUE.
      DFTGRA = .FALSE.
      GRDWT  = .FALSE.
      RNDMAX = 20.0D0
C     RNDMAX = 100.0D0
C
C     DFTOZER
C     =======
C
C     NRTYPF = 2
      NRTYPF = 1
C
C     DFTAML
C     ======
C
      NSCHEME = 15
C
C     *******************************************
C     ***** Prepare arguments for CONSTRUCT *****
C     *******************************************
C
      IBIG  = 200000
      IVECL = 200
C
      KSP    = 1
      KX8    = KSP    + 4*NUM
      KY8    = KX8    + IBIG
      KZ8    = KY8    + IBIG
      KW8    = KZ8    + IBIG
      KPSMU  = KW8    + IBIG 
      KRJ    = KPSMU  + IVECL*NAT
      KACCUM = KRJ    + IVECL*NAT
      KXMUJN = KACCUM + IVECL
      KXMUJ2 = KXMUJN + IVECL
      KRIJ   = KXMUJ2 + IVECL
      KAIJ   = KRIJ   + NAT*(NAT-1)/2
      KRBC   = KAIJ   + NAT*(NAT-1)/2
      KRVEC  = KRBC   + NAT*NAT
      KLAST  = KRVEC  + 3*IVECL
      IF (KLAST .GT. LWORK) CALL STOPIT('DFTGRA','CONSTRUC',KLAST,LWORK)
C
C     Open unit LUQUAD
C     ================
C
      LUQUAD = -1
#if defined (VAR_MPI)
      QUADNAME = 'DALTON.QUAD.'//chrnos(mynum/1000)
     &     //chrnos((mynum-(mynum/1000)*1000)/100)
     &     //chrnos((mynum-(mynum/100)*100)/10)
     &     //chrnos(mynum-(mynum/10)*10)
      LUQUAD = -1
      CALL GPOPEN(LUQUAD,QUADNAME,
     &     'UNKNOWN','SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUQUAD
C     
#else 
      CALL GPOPEN(LUQUAD,'DALTON.QUAD',
     &     'UNKNOWN','SEQUENTIAL','UNFORMATTED',IDUMMY,.FALSE.)
      REWIND LUQUAD
#endif
C
C     Calculate abscissas and weights
C     ===============================
C
      CALL CONSTRUC(WORK(KSP),WORK(KX8),WORK(KY8),WORK(KZ8),WORK(KW8),
     &              WORK(KPSMU),WORK(KRJ),WORK(KACCUM),
     &              WORK(KXMUJN),WORK(KXMUJ2),
     &              WORK(KRIJ),WORK(KAIJ),IVECL,IBIG,
     &              WORK(KRBC),WORK(KRVEC),LUQUAD,SPCDIV,IPRINT)
C
      IF (TEST) CALL DFTTST(WORK,LWORK,IBIG,LUQUAD,IPRINT)
C
C     Close unit LUQUAD
C     =================
C
      CALL GPCLOSE(LUQUAD,'KEEP')
C
C      CALL TIMER('DFTGRD',TIMSTR,TIMEND)
C
      RETURN
      END
C  /* Deck construc */
      SUBROUTINE CONSTRUC(WRKSP,X8,Y8,Z8,WT8,PSMU,RJ,ACCUM,XMUIJN,
     &                    XMUIJ2,RIJ,AIJ,IVECL,IBIG,RBC,RVEC,LUQUAD,
     &                    SPCDIV,IPRINT)
C-------------------------------------------------------------
C
C      Creates a quadrature grid and writes it to disk
C      Original version by C.W.Murray. Rewritten by 
C      A.M.Lee and D.J.Tozer.
C
C      Adapted for Dalton by T. Helgaker (with SPCDIV added)
C
C-------------------------------------------------------------
#include <implicit.h> 
#include <priunit.h> 
C
      LOGICAL SWITCH
      DIMENSION X8(IBIG),Y8(IBIG),Z8(IBIG),WT8(IBIG)
C
      PARAMETER (MAXCEN=100)
      COMMON /INFOA/ NAT,NUM,ZAN(MAXCEN),C(3,MAXCEN)
      LOGICAL KHNSHM,DFTGRA,GRDWT,SPCDIV
      COMMON /DFTCWM/ DFTGRA,GRDWT,RNDMAX,IQUAD,NRHX,NROX,NRSX,NTMAXX,
     &                NPHIX,KHNSHM
      COMMON /DFTOZER/ NRTYPF
      COMMON /DFTAML/ NSCHEME
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <nuclei.h>
#include <symmet.h>
C
      DIMENSION WRKSP(4*NUM),
     &          ACCUM(IVECL),XMUIJN(IVECL),XMUIJ2(IVECL),
     &          PSMU(IVECL,NAT),RJ(IVECL,NAT)
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
COLAV      DIMENSION BRAGG(0:54)
      DIMENSION BRAGG(0:103)	
      DIMENSION RVEC(IVECL,3)
      DIMENSION RADWT(2000), RADNDE(2000)
      DIMENSION ANGWT(9000), XANG(9000), YANG(9000), ZANG(9000)
      DIMENSION NRAR(6)
      DIMENSION EXPAZ(50,3)
      DIMENSION XPASC(20)
C
C     Data statements
C
      DATA (EXPAZ(I,1),I=1,10)/7.5D0,12.5D0,5.0D0,7.5D0,10.0D0,10.0D0,
     &                               12.5D0,15.0D0,17.5D0,20.0D0/
C
C   (I,1) I=13...17 values guessed...density innaccurate
      DATA (EXPAZ(I,1),I=11,18)/12.5D0,12.5D0,12.5D0,12.5D0,12.5D0
     1,12.5D0,12.5D0,12.5D0/
      DATA (EXPAZ(I,2),I=1,10)/12.5D0,12.5D0,7.5D0,12.5D0,15.0D0,17.5D0
     1,20.0D0,22.5D0,25.0D0,30.0D0/
      DATA (EXPAZ(I,2),I=11,18)/17.5D0,17.5D0,17.5D0,15.0D0,20.0D0
     1,25.0D0,17.5D0,25.0D0/
      DATA (EXPAZ(I,3),I=1,10)/17.5D0,22.5D0,10.0D0,20.0D0,20.0D0
     1,20.0D0,20.0D0,22.5D0,25.0D0,27.5D0/
      DATA (EXPAZ(I,3),I=11,18)/27.5D0,27.5D0,-1.0D0,25.0D0,30.0D0
     1,30.0D0,32.5D0,35.0D0/
C
      DATA BRAGG/ 0.75D0,
     1            0.35D0,0.35D0,
     11.45D0,1.05D0,0.85D0,0.70D0,0.65D0,0.60D0,0.50D0,0.45D0,
     11.80D0,1.50D0,1.25D0,1.10D0,1.00D0,1.00D0,1.00D0,1.00D0,
     12.20D0,1.80D0,
     11.60D0,1.40D0,1.35D0,1.40D0,1.40D0,
     11.40D0,1.35D0,1.35D0,1.35D0,1.35D0,
     1              1.30D0,1.25D0,1.15D0,1.15D0,1.15D0,1.15D0,
     2        18*1.3D0,
C
COLAV For Z=55-103 the data are taken from DIRAC's dftgrd.F:
     &2.60D0,2.15D0,1.95D0,
     &1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,1.80D0,
     &1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,
     &1.55D0,1.45D0,1.35D0,1.30D0,1.30D0,1.35D0,1.35D0,1.35D0,1.50D0,
     &1.90D0,1.75D0,1.60D0,1.90D0,1.50D0,1.50D0,
     &2.15D0,2.15D0,1.95D0,
     &1.80D0,1.80D0,1.37D0,1.75D0,1.75D0,1.75D0,1.75D0,
     &1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0/
C
#include <sdpre.h>
#include <ibtfun.h>
C
C
      IF (DFTGRA) CALL DZERO(WRKSP(NUM+1),NUM*3)
C
C     USE INPUTTED INFORMATION IN COMMON BLOCK DFTCWM TO DEFINE THE TYPE
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
C
      APASC=0.0D0
      DO 6 I=1,NTRANS+1
         XPASC(I)=XPASC(I)/DBLE(2*I-1)
         APASC=APASC+XPASC(I)
6     CONTINUE
      APASC=0.5D0/APASC
C
      ITEMP=0
      CALL DZERO(RBC,NAT*NAT)
      DO 7 INA=1,NAT
      DO 7 JNA=1,INA-1
          ITEMP=ITEMP+1
          RIJ(ITEMP) = 1.0D0/(DSQRT((C(1,INA)-C(1,JNA))**2
     &                             +(C(2,INA)-C(2,JNA))**2
     &                             +(C(3,INA)-C(3,JNA))**2))
          RBC(INA,JNA) = 1.0D0/RIJ(ITEMP)
          RBC(JNA,INA) = 1.0D0/RIJ(ITEMP)
          CHI=BRAGG(NINT(ZAN(INA)))/BRAGG(NINT(ZAN(JNA)))
          TEMP=(CHI-1.0D0)/(CHI+1.0D0)
          AIJ(ITEMP)=TEMP/(TEMP*TEMP-1.0D0)
          IF (AIJ(ITEMP).GT.0.5D0) AIJ(ITEMP)=0.5D0
          IF (AIJ(ITEMP).LT.-0.5D0) AIJ(ITEMP)=-0.5D0
7     CONTINUE
C
C     NRTYP1 determines the coarseness of the quadrature grid
C
      NRTYP1 = IQUAD
      CALL ASSIG(NRH,NRO,NRS,NRAR,NTMAX,NPMAX,NSCHEME,RHOTES,
     &           NRTYP1,NRTYPF,ILIM1,ILIM2)
C
C     Loop over atoms
C     ===============
C 
      IF (SPCDIV) THEN
         MXOPR = MAXOPR
      ELSE
         MXOPR = 0
      END IF
C
      NTOTAL = 0
      NALL = 0
      BOHR   = 0.529177249D0
      NATOM = 1
      DO 10 IATOM = 1, NUCIND
      DO 10 ISYMOP = 0, MXOPR         
      IF (IBTAND(ISYMOP,ISTBNU(IATOM)).EQ.0) THEN 
         IF (SPCDIV) THEN
            NDEG = 1
         ELSE
            NDEG  = MULT(ISTBNU(IATOM))
         END IF
C
         IPT   = 0
         INDEX = 0
         NDEX  = 1
         NTHIS = 0
         NTHAT = 0
C
C        Euler-Maclaurin radial quadrature
C        =================================
C
         NR    = NRWHAT(ZAN(NATOM),NRAR)
         BQEXP = BRAGG(NINT(ZAN(NATOM)))/BOHR
         DO IR = 1, NR
            XNODE = DBLE(IR)/DBLE(NR+1)
            RNODE = BQEXP*(XNODE/(1.0D0-XNODE))**2
            RWT   = 2.0D0*BQEXP**3*XNODE**5/(1.0D0-XNODE)**7
            RWT   = RWT/DBLE(NR+1)
            RADNDE(IR) = RNODE
            RADWT(IR)  = RWT
         END DO
C 
C        loop over radial points
C        =======================
C
         SWITCH = .FALSE.
         RBRAGG = BRAGG(NINT(ZAN(NATOM)))/(5.0D0*BOHR)
         DO 20 IR = 1, NR
            RWT   = RADWT(IR)
            RNODE = RADNDE(IR)
            IF (RNODE.GT.RNDMAX) GOTO 20
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
               IF(GRDWT)THEN
                  RVEC(IPT,1) = X8(INDEX)
                  RVEC(IPT,2) = Y8(INDEX)
                  RVEC(IPT,3) = Z8(INDEX)
               ENDIF
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
CTROND
         IF(.NOT.SPCDIV .AND. NDEG.GT.0) THEN
C          Scale weights with multiplicity of center
           FAC = SDPRE(NDEG) 
           CALL DSCAL(NTHIS,FAC,WT8,1)
         ENDIF
         NTHAT = NTHIS
         IF (SPCDIV) THEN
            CALL CMPRSQ(X8,Y8,Z8,WT8,NTHIS,KPOINT)
            NTHIS = KPOINT
         END IF
C
C        Write to file LUQUAD
C        ====================
C
         CALL WRTQUA(X8(1),Y8(1),Z8(1),WT8(1),LUQUAD,NTHIS)
         NTOTAL = NTOTAL + NTHIS
         NALL = NALL + NTHAT
         NATOM = NATOM + NDEG
         WRITE(LUPRI,'(A4,I3,3I8)') NAMN(IATOM),NDEG,NTHIS,NR,NANG
C
      END IF
10    CONTINUE
C
      NLAST = -1
      CALL WRTQUA(X8(1),Y8(1),Z8(1),WT8(1),LUQUAD,NLAST)
C
      WRITE (LUPRI,'(/,2X,A,I10,A,F5.1,A/)') 
     &  ' Number of grid points in quadrature:',NTOTAL,
     &  ' (',FLOAT(100*NTOTAL)/FLOAT(NALL),'%)'
C
C     Print file
C
      IF (IPRINT .GT. 20) CALL PRIQUA(X8,Y8,Z8,WT8,LUQUAD)
C
      RETURN
      END
C  /* Deck cmprsq */
      SUBROUTINE CMPRSQ(X,Y,Z,WT,NPOINT,KPOINT)
C
C     T. Helgaker
C
#include <implicit.h> 
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
#include <sdpre.h>
      J = 0
      DO I = 1, NPOINT
         NDEG = MLTPNT(X(I),Y(I),Z(I))
         IF (NDEG.GT.0) THEN
            J = J + 1 
            X(J) = X(I)
            Y(J) = Y(I)
            Z(J) = Z(I)
            WT(J) = SDPRE(NDEG)*WT(I) 
         END IF
      END DO
      KPOINT = J
      RETURN
      END
C  /* Deck wrtqua */
      SUBROUTINE WRTQUA(X,Y,Z,WT,LUQUAD,NPOINT)
C
C     T. Helgaker
C
#include <implicit.h> 
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      WRITE(LUQUAD) NPOINT
      IF (NPOINT .GT. 0) WRITE(LUQUAD) X,Y,Z,WT
      RETURN
      END
C  /* Deck reaqua */
      SUBROUTINE REAQUA(X,Y,Z,WT,LUQUAD,NPOINT)
C
C     T. Helgaker
C
#include <implicit.h> 
#include <priunit.h>
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      READ(LUQUAD) X,Y,Z,WT
      RETURN
      END
C  /* Deck priqua */
      SUBROUTINE PRIQUA(X,Y,Z,WT,LUQUAD)
C
C     T. Helgaker
C
#include <implicit.h> 
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
C  /* Deck priqu2 */
      SUBROUTINE PRIQU2(X,Y,Z,WT,LUQUAD,NPOINT)
C
C     T. Helgaker
C
#include <implicit.h> 
#include <priunit.h>
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      READ(LUQUAD) X,Y,Z,WT
      DO 100 I = 1, NPOINT
         WRITE (LUPRI,'(2X,4F20.12)') X(I),Y(I),Z(I),WT(I)
  100 CONTINUE
      RETURN
      END
C  /* Deck assig */
      SUBROUTINE ASSIG(NRH,NRO,NRS,NRAR,NTMAX,NPMAX,ISCHEME,
     &                 RHOTES,NRTYP1,NRTYPF,ILIM1,ILIM2)

#include <implicit.h> 
#include <priunit.h> 
C
      LOGICAL KHNSHM,DFTGRA,GRDWT
      COMMON /DFTCWM/ DFTGRA,GRDWT,RNDMAX,IQUAD,NRHX,NROX,NRSX,NTMAXX,
     &                NPHIX,KHNSHM
C
      COMMON /DFTAML/ NSCHEME
      DIMENSION NRAR(6)
C 
C    CWM, AML and DJT
C
      IF (IQUAD.GT.0) THEN
          ILIM1=IQUAD
          ILIM2=IQUAD
      ELSE
          ILIM1=1
          ILIM2=5
      ENDIF
      NRTYP1=IQUAD

      IF (.NOT.KHNSHM) NRTYP1=NRTYPF

      IF (NRTYP1.EQ.6) THEN
       NRO=NROX
       NRH=NRHX
       NRS=NRSX
       NRAR(1)=NRH
       NRAR(2)=NRO
       NRAR(3)=NRS
       NRAR(4)=NRS+10
       NRAR(5)=NRS+20
       NRAR(6)=NRS+30
       NTMAX=NTMAXX
       NPMAX=NPHIX
       RHOTES=1.0D-12
C      ISCHEME is either assigned or uses default value for LEBEDEV angular
C      quadratures
      ELSE
      IF (NRTYP1.EQ.1) THEN
       NRO=64
       NRH=60
       NRS=72
       NRAR(1)=60
       NRAR(2)=64
       NRAR(3)=72
       NRAR(4)=96
       NRAR(5)=110
       NRAR(6)=120
       NTMAX=24
       NPMAX=48
       RHOTES=1.0D-10
      ENDIF
       IF(NRTYP1.EQ.7)THEN
       NRO=96
       NRH=72
       NRS=108
       NRAR(1)=72
       NRAR(2)=96
       NRAR(3)=108
       NRAR(4)=120
       NRAR(5)=136
       NRAR(6)=160
       NTMAX=54
       NPMAX=108
       RHOTES=1d-12
       ENDIF
      IF (NRTYP1.EQ.2.OR.NRTYP1.EQ.3) THEN
       NRO=25
       NRH=20
       NRS=30
       NRAR(1)=20
       NRAR(2)=25
       NRAR(3)=30
       NRAR(4)=36
       NRAR(5)=42
       NRAR(6)=48
       NTMAX=8
       NPMAX=16
       RHOTES=1.0D-6
      ENDIF
      IF (NRTYP1.EQ.4) THEN
       NRAR(1)=30
       NRAR(2)=36
       NRAR(3)=42
       NRAR(4)=50
       NRAR(5)=60
       NRAR(6)=70
       NRO=36 
       NRH=30
       NRS=42
       NTMAX=12
       NPMAX=24
       RHOTES=1.0D-7
      ENDIF
      IF (NRTYP1.EQ.5) THEN
       NRO=54
       NRH=46
       NRS=60
       NRAR(1)=46
       NRAR(2)=54
       NRAR(3)=60
       NRAR(4)=68
       NRAR(5)=76 
       NRAR(6)=84 
       NTMAX=16 
       NPMAX=32 
       RHOTES=1.0D-8 
       ENDIF
      ENDIF
      ISCHEME=NLEBSET(NTMAX)
      RETURN
      END
C  /* Deck nrwhat */
      FUNCTION NRWHAT(ZAN,NRAR)
#include <implicit.h> 
      INTEGER NRWHAT , NRAR(6)
      IF(ZAN.LE.2.0D0)NR=NRAR(1)
      IF(ZAN.GT.2.0D0.AND.ZAN.LE.10.0D0)NR=NRAR(2)
      IF(ZAN.GT.10.0D0.AND.ZAN.LE.18.1D0)NR=NRAR(3)
      IF(ZAN.GT.18.0D0.AND.ZAN.LE.36.0D0)NR=NRAR(4)
      IF(ZAN.GT.36.0D0.AND.ZAN.LE.54.0D0)NR=NRAR(5)
      IF(ZAN.GE.55.0D0)NR=NRAR(6)
     
      IF(ZAN.EQ.0D0)NR=NRAR(2)

      NRWHAT=NR
      RETURN
      END
C  /* Deck nleb */
      FUNCTION NLEB(NSCHEME)
#include <implicit.h> 
C
      IF (NSCHEME.EQ.9) NTOT=38
      IF (NSCHEME.EQ.11) NTOT=50
      IF (NSCHEME.EQ.13) NTOT=78
      IF (NSCHEME.EQ.15) NTOT=86
      IF (NSCHEME.EQ.17) NTOT=116
      IF (NSCHEME.EQ.19) NTOT=146
      IF (NSCHEME.EQ.23) NTOT=194
      IF (NSCHEME.EQ.29) NTOT=302
      IF (NSCHEME.EQ.41) NTOT=590
      IF (NSCHEME.EQ.47) NTOT=770
      IF (NSCHEME.EQ.53) NTOT=974
C    
      NLEB=NTOT
      RETURN
      END
C  /* Deck nlebset2 */
      FUNCTION NLEBSET2(NTOT)
#include <implicit.h> 
C
       IF (NTOT.GE.973) THEN
      NSCHEME=53
       ELSEIF (NTOT.GE.769) THEN
      NSCHEME=47
       ELSEIF (NTOT.GE.589) THEN
      NSCHEME=41
       ELSEIF (NTOT.GE.301) THEN
      NSCHEME=29
       ELSEIF (NTOT.GE.193) THEN
      NSCHEME=23
       ELSEIF (NTOT.GE.145) THEN
      NSCHEME=19
       ELSEIF (NTOT.GE.115) THEN
      NSCHEME=17
       ELSEIF (NTOT.GE.85) THEN
      NSCHEME=15
       ELSEIF (NTOT.GE.77) THEN
      NSCHEME=13
       ELSEIF (NTOT.GE.49) THEN
      NSCHEME=11
       ELSE
      NSCHEME=9
       ENDIF
      NLEBSET2=NSCHEME
      RETURN
      END
C  /* Deck nlebset */
      FUNCTION NLEBSET(NTMAX)
#include <implicit.h> 
      IF (NTMAX.LE.5) THEN
        ISCHEME=9
      ELSEIF (NTMAX.LE.6) THEN
        ISCHEME=11
      ELSEIF (NTMAX.LE.7) THEN
        ISCHEME=13
      ELSEIF (NTMAX.LE.8) THEN
        ISCHEME=15
      ELSEIF (NTMAX.LE.9) THEN
        ISCHEME=17
      ELSEIF (NTMAX.LE.10) THEN
        ISCHEME=19
      ELSEIF (NTMAX.LE.12) THEN
        ISCHEME=23
      ELSEIF (NTMAX.LE.17) THEN
        ISCHEME=29
      ELSEIF (NTMAX.LE.35) THEN
        ISCHEME=41
      ELSEIF (NTMAX.LE.52) THEN
        ISCHEME=47
      ELSE
        ISCHEME=53
      ENDIF
      NLEBSET=ISCHEME

      RETURN
      END
C  /* Deck sleb */
      SUBROUTINE SLEB(XNEW,YNEW,ZNEW,WTNEW,NSCHEME,NNEQV,NATOM,
     &RNODE,RBRAGG)
#include <implicit.h> 
#include <priunit.h> 
      PARAMETER (MAXNUM=975,MAXCLS=32)
      DIMENSION WTNEW(MAXNUM),XNEW(MAXNUM),YNEW(MAXNUM),ZNEW(MAXNUM),
     &WT(MAXNUM),X(MAXNUM),Y(MAXNUM),Z(MAXNUM),
     &II(48)
      DIMENSION XB(12),XC(12),XD(12),XM(12),XL(12),
     &P(12),Q(12),R(12),U(12),W(12)
      DIMENSION IUNIQUE(MAXCLS,48),IDEGEN(MAXCLS,48),IMAP(MAXCLS,48)
      DIMENSION INEQVC(60)
      PARAMETER (MAXCEN=100)
      COMMON /INFOA/ NAT,NUM,ZAN(MAXCEN),C(3,MAXCEN)
      LOGICAL NEGWT,EQUIV,TEQUIV
C*******************************************************************
C     Lebedev scheme programmed by Aaron M. Lee, Cambridge 1993
C     Principal references:
C     (1) Lebedev, V.I., Zh. Vychisl. Mat. mat. Fiz., 1976, 16, 293.
C     (2) Lebedev, V.I., Zh. Vychisl. Mat. mat. Fiz., 1975, 15, 48.
C     (3) Lebedev, V.I., Sib. Mat. Zh., 1975, 18(1), 132-142.
C     Scheme designed to exactly integrate functions which are linear
C     combinations of spherical harmonics up to a given order,NSCHEME.

C     All symmetry by Aaron M. Lee
C*******************************************************************
 
      PI=DACOS(-1.0D0)

      NEGWT=.FALSE.
      NTOT=NLEB(NSCHEME)

C     Chose the Lebedev grid to use as a function of the radius (rnode)
C     closer to the nuclei scale the number of angular points as

       NOLD=NSCHEME
      IF (NSCHEME.GT.9.AND.RNODE.LT.RBRAGG) THEN
       IANG=NINT(NTOT*RNODE/RBRAGG)
       NTOT=MIN(NTOT,IANG)

       NSCHEME=NLEBSET2(NTOT)
       NTOT=NLEB(NSCHEME)

      ENDIF

      DO K=1,NTOT
       X(K)=0.0D0
       Y(K)=0.0D0
       Z(K)=0.0D0
       WT(K)=0.0D0
       XNEW(K)=0.0D0
       YNEW(K)=0.0D0
       ZNEW(K)=0.0D0
       WTNEW(K)=0.0D0
      ENDDO
       

      
C     Some values redefined because of poor precision, namely
C     where the resulting nodes lie off the unit sphere. If values
C     were reported, these are still preserved in the code.
C     Scheme 13.1(neg weight) replaced by 13.2
C     Scheme 17.1(useless) replaced by 17.2
C     Schemes 25 and 27 negative weight quadratures and therefore unstable.

      IF(NSCHEME.EQ.9)THEN
      A1=1.0D0/105.0D0
      A2=0.0D0
      A3=9.0D0/280.0D0
      XC(1)=1.0D0/35.0D0
      P(1)=0.888073833977D0
      Q(1)=0.459700843381D0
      N1=0 
      N2=1 
      N3=0
      NTOT=14+N1*24+N2*24+N3*48
      ENDIF     

      IF(NSCHEME.EQ.11)THEN
      A1=4.0D0/315.0D0
      A2=64.0D0/2835.0D0
      A3=27.0D0/1280.0D0
      XB(1)=11.0D0*11.0D0*11.0D0*11.0D0/725760.0D0
      XL(1)=0.301511344578D0
      XM(1)=0.904534033733D0
      N1=1 
      N2=0 
      N3=0
      NTOT=26+N1*24+N2*24+N3*48
      ENDIF    

      IF(NSCHEME.EQ.13)THEN
      A1=0.0138665921047D0
      A2=0.0D0
      A3=0.0D0
      XB(1)=0.0130509318626D0
      XB(2)=0.0132064232231D0
      XC(1)=0.0119426635549D0
      XM(1)=0.914152532416D0
      XM(2)=0.359236381200D0
      XL(1)=0.286640146767D0
      XL(2)=0.659905001656D0
      P(1)=0.841991943785D0
      Q(1)=0.539490098706D0
      N1=2 
      N2=1 
      N3=0 
      NTOT=6+N1*24+N2*24+N3*48
      ENDIF     

      IF(NSCHEME.EQ.15)THEN
      A1=0.0115440115441D0
      A2=0.0D0
      A3=0.0119439090859D0
      XB(1)=0.0111105557106D0
      XB(2)=0.0118765012945D0
      XC(1)=0.0118123037469D0
      XL(1)=0.369602846454D0
      XL(2)=0.694354006603D0
      XM(1)=0.852518311701D0
      XM(2)=0.189063552885D0
      P(1)=0.927330657151D0
      Q(1)=0.374243039090D0
      N1=2 
      N2=1 
      N3=0 
      NTOT=14+N1*24+N2*24+N3*48
      ENDIF      

      IF(NSCHEME.EQ.17)THEN
      A1=0.0D0
      A2=0.00200918797730D0
      A3=0.00988550016044D0
      XB(1)=0.00844068048232D0
      XB(2)=0.00987390742389D0
      XB(3)=0.00935732169000D0
      XC(1)=0.00969499636166D0
      XM(1)=0.973314565209D0
      XM(2)=0.840255982384D0
      XM(3)=0.238807866929D0
      XL(1)=0.162263300152D0
      XL(2)=0.383386152638D0
      XL(3)=0.686647945709D0
      P(1)=0.878158910604D0
      Q(1)=0.478369028812D0
      N1=3 
      N2=1 
      N3=0
      NTOT=20+N1*24+N2*24+N3*48
      ENDIF 


      IF(NSCHEME.EQ.19)THEN
      A1=5.99631368862D-4
      A2=7.37299971862D-3
      A3=7.21051536014D-3
      XB(1)=7.57439415905D-3
      XB(2)=6.75382948631D-3
      XB(3)=7.11635549312D-3
      XD(1)=6.99108735330D-3
      XL(1)=0.157467667204D0
      XL(2)=0.417496122797D0
      XL(3)=0.676441040011D0
      XM(1)=0.974888643677D0
      XM(2)=0.807089818360D0
      XM(3)=0.291298882210D0
      R(1)=0.882270011260D0
      U(1)=0.140355381171D0
      W(1)=0.449332832327D0
      N1=3  
      N2=0  
      N3=1
      NTOT=26+N1*24+N2*24+N3*48
      ENDIF     

      IF(NSCHEME.EQ.23)THEN
      A1=1.78234044724D-3
      A2=5.71690594998D-3
      A3=5.57338317884D-3
      XB(1)=5.51877146727D-3
      XB(2)=5.15823771181D-3
      XB(3)=5.60870408259D-3
      XB(4)=4.10677702817D-3
      XC(1)=5.05184606462D-3
      XD(1)=5.53024891623D-3
      XL(1)=0.444693317871D0
      XL(2)=0.289246562758D0
      XL(3)=0.671297344270D0
      XL(4)=0.129933544765D0
      XM(1)=0.777493219315D0
      XM(2)=0.912509096867D0
      XM(3)=0.314196994183D0
      XM(4)=0.982972302707D0
      P(1)=0.938319218138D0
      Q(1)=0.345770219761D0
      R(1)=0.836036015482D0
      U(1)=0.159041710538D0
      W(1)=0.525118572443D0
      N1=4 
      N2=1 
      N3=1
      NTOT=26+N1*24+N2*24+N3*48
      ENDIF     

      IF(NEGWT.AND.NSCHEME.EQ.25)THEN
      A1=-5.52263991974D-2
      A2=0.0D0
      A3=4.45027460745D-3
      XB(1)=3.97640801805D-3
      XB(2)=4.40140065038D-3
      XB(3)=4.49684106792D-3
      XB(4)=5.04915345048D-3
      XB(5)=1.72454435055D-2
      XC(1)=5.19806986406D-3
      XC(2)=4.23108309536D-3
      XD(1)=4.69572097257D-3
      XL(1)=0.698190665845D0
      XL(2)=0.658740524346D0
      XL(3)=0.449204468740D0
      XL(4)=0.252041949021D0
      XL(5)=4.03854405009D-2
      XM(1)=0.158302205463D0
      XM(2)=0.363485684957D0
      XM(3)=0.772289253148D0
      XM(4)=0.934317778846D0
      XM(5)=0.998367683968D0
      P(1)=0.935022745881D0
      P(2)=0.812913653173D0
      Q(1)=0.354587739052D0
      Q(2)=0.582384230972D0
      R(1)=0.486466535887D0
      U(1)=0.843636521069D0
      W(1)=0.227218180900D0
      N1=5 
      N2=2 
      N3=1 
      NTOT=14+N1*24+N2*24+N3*48
      ENDIF     

      IF(NEGWT.AND.NSCHEME.EQ.27)THEN
      A1=-1.31376912733D-3
      A2=-2.52272870489D-3
      A3=4.18685388170D-3
      XB(1)=5.31516797782D-3
      XB(2)=4.25613135143D-3
      XB(3)=4.11248239441D-3
      XB(4)=3.59558489976D-3
      XB(5)=4.04714237709D-3
      XC(1)=4.22958270065D-3
      XD(1)=4.07146759383D-3
      XD(2)=4.08091422578D-3
      XL(1)=0.703937339159D0
      XL(2)=0.662033866370D0
      XL(3)=0.464744872642D0
      XL(4)=0.327742065497D0
      XL(5)=0.101252624857D0  
      XM(1)=0.0945750764036D0
      XM(2)=0.351315128565D0
      XM(3)=0.753673939251D0
      XM(4)=0.886098344997D0
      XM(5)=0.986948074629D0
      XL(5)=DSQRT(0.50D0*(1.0D0-XM(5)*XM(5)))
      P(1)=0.850650808352D0
      Q(1)=0.525731112119D0
      R(1)=0.819343388819D0
      U(1)=0.524493924092D0
      W(1)=0.231479015871D0
      R(2)=0.939227929750D0
      U(2)=0.323348454269D0
      W(2)=0.115311196765D0
      N1=5 
      N2=1 
      N3=2 
      NTOT=26+N1*24+N2*24+N3*48
      ENDIF     

      IF(NSCHEME.EQ.29)THEN
      A1=8.54591172878D-4
      A2=0.0D0
      A3=3.59911928502D-3
      XB(1)=3.65004580768D-3
      XB(2)=3.60482260142D-3
      XB(3)=3.57672966173D-3
      XB(4)=3.44978842429D-3
      XB(5)=3.10895312238D-3
      XB(6)=2.35210141366D-3
      XC(1)=3.60082093222D-3
      XC(2)=2.98234496317D-3
      XD(1)=3.57154055427D-3
      XD(2)=3.39231220501D-3
      XL(1)=0.701176641609D0
      XL(2)=0.656632941022D0
      XL(3)=0.472905413258D0
      XL(4)=0.351564034558D0
      XL(5)=0.221964523631D0
      XL(6)=0.0961830852303D0
      XM(1)=0.129238672710D0
      XM(2)=0.371034178385D0
      XM(3)=0.743452042987D0
      XM(4)=0.867643624544D0
      XM(5)=0.949454317226D0
      XM(6)=0.990705621379D0
      P(1)=0.820326419828D0
      P(2)=0.964408914879D0
      Q(1)=0.571895589188D0
      Q(2)=0.264415288706D0
      R(1)=0.251003475177D0
      R(2)=0.902442529533D0
      U(1)=0.800072749407D0
      U(2)=0.412772408317D0
      W(1)=0.544867737258D0
      W(2)=0.123354853258D0
      N1=6 
      N2=2 
      N3=2 
      NTOT=14+N1*24+N2*24+N3*48
      ENDIF      


      IF(NSCHEME.EQ.41)THEN
      A1=0.3095121295d-3
      A2=0.0D0
      A3=0.1852379698d-2

      XB(1)=0.9764331164d-3
      XB(2)=0.1384737234d-2
      XB(3)=0.1617210647d-2
      XB(4)=0.1749564657d-2
      XB(5)=0.1818471778d-2
      XB(6)=0.1846715956d-2
      XB(7)=0.1852028828d-2
      XB(8)=0.1858812585d-2
      XB(9)=0.1871790639d-2

      XC(1)=0.1300321685d-2
      XC(2)=0.1705153996d-2
      XC(3)=0.1857161196d-2

      XD(1)=0.1555213603d-2
      XD(2)=0.1802239128d-2
      XD(3)=0.1849830560d-2
      XD(4)=0.1713904507d-2
      XD(5)=0.1802658934d-2
      XD(6)=0.1842866472d-2

      XL(1)=0.6095034115d-1
      XL(2)=0.1459036449d0
      XL(3)=0.2384736701d0
      XL(4)=0.3317920736d0
      XL(5)=0.4215761784d0
      XL(6)=0.5044419707d0
      XL(7)=0.6372546939d0
      XL(8)=0.6807744066d0
      XL(9)=0.7040954938d0

      XM(1)=0.9962781297d0
      XM(2)=0.9784805837d0
      XM(3)=0.9414141582d0
      XM(4)=0.8830787279d0
      XM(5)=0.8028368773d0
      XM(6)=0.7007685753d0
      XM(7)=0.4333738687d0
      XM(8)=0.2703560883d0
      XM(9)=0.9219040707d-1

      P(1)=0.9850133350d0
      P(2)=0.9180452877d0
      P(3)=0.7911019296d0

      Q(1)=0.1724782009d0
      Q(2)=0.3964755348d0
      Q(3)=0.6116843442d0

      U(1)=0.8213021581d-1
      U(2)=0.8999205842d-1
      U(3)=0.1816640840d0
      U(4)=0.1720795225d0
      U(5)=0.2634716655d0
      U(6)=0.3518280927d0

      R(1)=0.2778673190d0
      R(2)=0.5033564271d0
      R(3)=0.5984126497d0
      R(4)=0.3791035407d0
      R(5)=0.4742392842d0
      R(6)=0.5610263808d0

      W(1)=0.9571020743d0
      W(2)=0.8593798558d0
      W(3)=0.7803207424d0
      W(4)=0.9092134750d0
      W(5)=0.8400474883d0
      W(6)=0.7493106119d0

      N1=9 
      N2=3 
      N3=6 
      NTOT=14+N1*24+N2*24+N3*48
      ENDIF     

      IF(NSCHEME.EQ.47)THEN
      A1=0.2192942090D-3
      A2=0.1436433617D-2
      A3=0.1421940344D-2

      XB(1)=0.6798123510D-3
      XB(2)=0.9913184235D-3
      XB(3)=0.1180207833D-2
      XB(4)=0.1296599602D-2
      XB(5)=0.1365871427D-2
      XB(6)=0.1402988604D-2
      XB(7)=0.1418645563D-2
      XB(8)=0.1421376741D-2
      XB(9)=0.1423996475D-2
      XB(10)=0.1431554042D-2

      XC(1)=0.9254401499D-3
      XC(2)=0.1250239995D-2
      XC(3)=0.1394365843D-2

      XD(1)=0.1127089094D-2
      XD(2)=0.1345753761D-2
      XD(3)=0.1424957283D-2
      XD(4)=0.1261523341D-2
      XD(5)=0.1392547106D-2
      XD(6)=0.1418761677D-2
      XD(7)=0.1338366684D-2
      XD(8)=0.1393700862D-2
      XD(9)=0.1415914757D-2

      XL(1)=0.5087204410D-1
      XL(2)=0.1228198790D0
      XL(3)=0.2026890814D0
      XL(4)=0.2847745156D0
      XL(5)=0.3656719078D0
      XL(6)=0.4428264886D0
      XL(7)=0.5140619627D0
      XL(8)=0.6306401219D0
      XL(9)=0.6716883332D0
      XL(10)=0.6979792685D0

      XM(1)=0.9974086776D0
      XM(2)=0.9847997535D0
      XM(3)=0.9580366759D0
      XM(4)=0.9153179504D0
      XM(5)=0.8559019286D0
      XM(6)=0.7796213195D0
      XM(7)=0.6866444472D0
      XM(8)=0.4523119203D0
      XM(9)=0.3125213050D0 
      XM(10)=0.1601558034D0 

      P(1)=0.9894775374D0
      P(2)=0.9407768787D0
      P(3)=0.8457493051D0

      Q(1)=0.1446865674D0
      Q(2)=0.3390263475D0
      Q(3)=0.5335804651D0

      U(1)=0.6944024393D-1
      U(2)=0.2269004109D0
      U(3)=0.8025574608D-1
      U(4)=0.1467999527D0
      U(5)=0.1571507769D0
      U(6)=0.2365702993D0
      U(7)=0.7714815844D-1
      U(8)=0.3062936666D0
      U(9)=0.3822477379D0

      R(1)=0.2355187894D0
      R(2)=0.4102182474D0
      R(3)=0.6214302417D0
      R(4)=0.3245284345D0
      R(5)=0.5224482189D0
      R(6)=0.6017546634D0
      R(7)=0.4346575516D0
      R(8)=0.4908826589D0
      R(9)=0.5648768149D0

      W(1)=0.9693858634D0
      W(2)=0.8833103605D0
      W(3)=0.7793481057D0
      W(4)=0.9344148270D0
      W(5)=0.8380641334D0
      W(6)=0.7628406246D0
      W(7)=0.8972853361D0
      W(8)=0.8156092232D0
      W(9)=0.7313007936D0

      N1=10
      N2=3
      N3=9
      NTOT=26+N1*24+N2*24+N3*48
      ENDIF     

      IF(NSCHEME.EQ.53)THEN
      A1=0.1438294190D-3
      A2=0.0D0               
      A3=0.1125772288D-2

      XB(1)=0.4948029342D-3
      XB(2)=0.7357990108D-3
      XB(3)=0.8889132771D-3
      XB(4)=0.9888347838D-3
      XB(5)=0.1053299681D-2
      XB(6)=0.1092778807D-2
      XB(7)=0.1114389394D-2
      XB(8)=0.1123724788D-2
      XB(9)=0.1125239325D-2
      XB(10)=0.1126153271D-2
      XB(11)=0.1130286931D-2
      XB(12)=0.1134986534D-2

      XC(1)=0.6823367927D-3
      XC(2)=0.9454158160D-3
      XC(3)=0.1074429975D-2
      XC(4)=0.1129300086D-2

      XD(1)=0.8436884500D-3
      XD(2)=0.1075255720D-2
      XD(3)=0.1108577236D-2
      XD(4)=0.9566475323D-3
      XD(5)=0.1080663250D-2
      XD(6)=0.1126797131D-2
      XD(7)=0.1022568715D-2
      XD(8)=0.1108960267D-2
      XD(9)=0.1122790653D-2
      XD(10)=0.1032401847D-2
      XD(11)=0.1107249382D-2
      XD(12)=0.1121780048D-2

      XL(1)=0.4292963545D-1
      XL(2)=0.1051426854D0
      XL(3)=0.1750024867D0
      XL(4)=0.2477653379D0
      XL(5)=0.3206567123D0
      XL(6)=0.3916520749D0
      XL(7)=0.4590825874D0
      XL(8)=0.5214563888D0
      XL(9)=0.6253170244D0
      XL(10)=0.6637926744D0
      XL(11)=0.6910410398D0
      XL(12)=0.7052907007D0

      XM(1)=0.9981553450D0
      XM(2)=0.9888832243D0
      XM(3)=0.9688902204D0
      XM(4)=0.9366027304D0
      XM(5)=0.8912679426D0
      XM(6)=0.8325967237D0
      XM(7)=0.7605829053D0
      XM(8)=0.6754009691D0
      XM(9)=0.4668589056D0 
      XM(10)=0.3446136545D0 
      XM(11)=0.2119541518D0 
      XM(12)=0.7162440144D-1

      P(1)=0.9923235654D0
      P(2)=0.9557815124D0
      P(3)=0.8827859807D0
      P(4)=0.7737784472D0

      Q(1)=0.1236686762D0
      Q(2)=0.2940777114D0
      Q(3)=0.4697753849D0
      Q(4)=0.6334563241D0

      U(1)=0.5974048614D-1
      U(2)=0.1375760408D0
      U(3)=0.3391016526D0
      U(4)=0.1271675191D0
      U(5)=0.2693120740D0
      U(6)=0.1419786452D0
      U(7)=0.6709284600D-1
      U(8)=0.7057738183D-1
      U(9)=0.2783888477D0
      U(10)=0.1979578938D0
      U(11)=0.2087307061D0
      U(12)=0.4055122137D0

      R(1)=0.2029128752D0
      R(2)=0.4602621942D0
      R(3)=0.5030673999D0
      R(4)=0.2817606422D0
      R(5)=0.4331561291D0
      R(6)=0.6256167328D0
      R(7)=0.3798395291D0
      R(8)=0.5517505421D0
      R(9)=0.6029619156D0
      R(10)=0.3589606329D0
      R(11)=0.5348666438D0
      R(12)=0.5674997546D0

      W(1)=0.9773727228D0
      W(2)=0.8770584618D0
      W(3)=0.7949422999D0
      W(4)=0.9510201693D0
      W(5)=0.8601434616D0
      W(6)=0.7671021862D0
      W(7)=0.9226161107D0
      W(8)=0.8310175524D0
      W(9)=0.7476206108D0
      W(10)=0.9121183784D0
      W(11)=0.8187485362D0
      W(12)=0.7165918454D0

      N1=12
      N2=4
      N3=12
      NTOT=14+N1*24+N2*24+N3*48

      ENDIF

C **********************************************************************
C     The following code, for points of index={1,26}, is general for all
C     schemes given appropriate values of {Ak}.
C **********************************************************************

C     Vertices of the octahedron
      NCOUNT=1 
      IF (A1.NE.0) THEN
      X(NCOUNT)=0.0D0
      Y(NCOUNT)=0.0D0
      Z(NCOUNT)=1.0D0

      X(NCOUNT+1)=0.0D0
      Y(NCOUNT+1)=0.0D0
      Z(NCOUNT+1)=-1.0D0

      X(NCOUNT+2)=0.0D0 
      Y(NCOUNT+2)=1.0D0
      Z(NCOUNT+2)=0.0D0

      X(NCOUNT+3)=0.0D0 
      Y(NCOUNT+3)=-1.0D0
      Z(NCOUNT+3)=0.0D0

      X(NCOUNT+4)=1.0D0
      Y(NCOUNT+4)=0.0D0
      Z(NCOUNT+4)=0.0D0

      X(NCOUNT+5)=-1.0D0
      Y(NCOUNT+5)=0.0D0
      Z(NCOUNT+5)=0.0D0
      DO IVERT=0,5
       WT(IVERT+NCOUNT)=A1
      ENDDO
      NCOUNT=NCOUNT+6
      ENDIF
C     Midpoints of the edges of the octahedron 
      IF (A2.NE.0) THEN
      S2=1.0D0/DSQRT(2.0D0)
      X(NCOUNT)=S2
      Y(NCOUNT)=S2
      Z(NCOUNT)=0.0D0

      X(NCOUNT+1)=-S2
      Y(NCOUNT+1)=S2
      Z(NCOUNT+1)=0.0D0

      X(NCOUNT+2)=-S2
      Y(NCOUNT+2)=-S2
      Z(NCOUNT+2)=0.0D0

      X(NCOUNT+3)=S2
      Y(NCOUNT+3)=-S2
      Z(NCOUNT+3)=0.0D0

      X(NCOUNT+4)=S2
      Y(NCOUNT+4)=0.0D0
      Z(NCOUNT+4)=S2

      X(NCOUNT+5)=-S2
      Y(NCOUNT+5)=0.0D0
      Z(NCOUNT+5)=S2

      X(NCOUNT+6)=-S2
      Y(NCOUNT+6)=0.0D0
      Z(NCOUNT+6)=-S2

      X(NCOUNT+7)=S2
      Y(NCOUNT+7)=0.0D0
      Z(NCOUNT+7)=-S2

      X(NCOUNT+8)=0.0D0
      Y(NCOUNT+8)=S2
      Z(NCOUNT+8)=S2

      X(NCOUNT+9)=0.0D0
      Y(NCOUNT+9)=-S2
      Z(NCOUNT+9)=S2

      X(NCOUNT+10)=0.0D0
      Y(NCOUNT+10)=-S2
      Z(NCOUNT+10)=-S2

      X(NCOUNT+11)=0.0D0
      Y(NCOUNT+11)=S2
      Z(NCOUNT+11)=-S2

      DO IMID=0,11
       I=IMID+NCOUNT
       WT(I)=A2
      ENDDO
      NCOUNT=NCOUNT+12
      ENDIF

C      The midpoints of the faces
      IF (A3.NE.0) THEN
      S3=1.0D0/DSQRT(3.0D0)
      J1=-1
      DO I1=0,1
       J1=-J1
       J2=-1
       DO I2=0,1
        J2=-J2
        J3=-1
        DO I3=0,1
         J3=-J3
         X(NCOUNT)=S3*DBLE(J1)
         Y(NCOUNT)=S3*DBLE(J2)
         Z(NCOUNT)=S3*DBLE(J3)
         WT(NCOUNT)=A3
         NCOUNT=NCOUNT+1
        ENDDO
       ENDDO
      ENDDO
      ENDIF

C########################################################################
C     on the bisectors of the faces
      DO K=0,N1-1
       XLR=XL(K+1)
       XMR=XM(K+1)
       J1=-1
       DO I1=0,1   
        J1=-J1
        J2=-1
        DO I2=0,1
         J2=-J2
         J3=-1
         DO I3=0,1
          J3=-J3
          KK=I1+2*I2+4*I3

          X(NCOUNT+KK)=DBLE(J1)*XLR
          Y(NCOUNT+KK)=DBLE(J2)*XLR
          Z(NCOUNT+KK)=DBLE(J3)*XMR

          X(NCOUNT+8+KK)=DBLE(J1)*XMR
          Y(NCOUNT+8+KK)=DBLE(J2)*XLR
          Z(NCOUNT+8+KK)=DBLE(J3)*XLR

          X(NCOUNT+16+KK)=DBLE(J1)*XLR
          Y(NCOUNT+16+KK)=DBLE(J2)*XMR
          Z(NCOUNT+16+KK)=DBLE(J3)*XLR

          WT(NCOUNT+KK)=XB(K+1)
          WT(NCOUNT+8+KK)=XB(K+1)
          WT(NCOUNT+16+KK)=XB(K+1)
         ENDDO
        ENDDO
       ENDDO
       NCOUNT=NCOUNT+24
      ENDDO
C     on the edges of the faces
      DO K=0,N2-1
      RP=P(K+1)
      RQ=Q(K+1)
       J1=-1
       DO I1=0,1   
        J1=-J1
        J2=-1
        DO I2=0,1
         J2=-J2
         J3=-1
         DO I3=0,1
          J3=-J3
          IHH=I1+2*I2
C1
          X(NCOUNT+IHH)=RP*DBLE(J1)
          Y(NCOUNT+IHH)=RQ*DBLE(J2)
          Z(NCOUNT+IHH)=0.0D0
C2
          X(NCOUNT+4+IHH)=RP*DBLE(J1)
          Y(NCOUNT+4+IHH)=0.0D0
          Z(NCOUNT+4+IHH)=RQ*DBLE(J2)
C3
          X(NCOUNT+8+IHH)=0.0D0
          Y(NCOUNT+8+IHH)=RP*DBLE(J1)
          Z(NCOUNT+8+IHH)=RQ*DBLE(J2)
C4
          X(NCOUNT+12+IHH)=RQ*DBLE(J1)
          Y(NCOUNT+12+IHH)=RP*DBLE(J2)
          Z(NCOUNT+12+IHH)=0.0D0
C5
          X(NCOUNT+16+IHH)=RQ*DBLE(J1)
          Y(NCOUNT+16+IHH)=0.0D0
          Z(NCOUNT+16+IHH)=RP*DBLE(J2)
C6
          X(NCOUNT+20+IHH)=0.0D0
          Y(NCOUNT+20+IHH)=RQ*DBLE(J1)
          Z(NCOUNT+20+IHH)=RP*DBLE(J2)

          WT(NCOUNT+IHH)=XC(K+1)
          WT(NCOUNT+IHH+4)=XC(K+1)
          WT(NCOUNT+IHH+8)=XC(K+1)
          WT(NCOUNT+IHH+12)=XC(K+1)
          WT(NCOUNT+IHH+16)=XC(K+1)
          WT(NCOUNT+IHH+20)=XC(K+1)
         ENDDO
        ENDDO
       ENDDO
      NCOUNT=NCOUNT+24
      ENDDO
      DO K=0,N3-1
      RR=R(K+1)
      RU=U(K+1)
      RW=W(K+1)
       J1=-1
       DO I1=0,1   
        J1=-J1
        J2=-1
        DO I2=0,1
         J2=-J2
         J3=-1
         DO I3=0,1
          J3=-J3
          KK=I1+2*I2+4*I3
          X(NCOUNT+KK)=DBLE(J1)*RR
          Y(NCOUNT+KK)=DBLE(J2)*RU
          Z(NCOUNT+KK)=DBLE(J3)*RW

          X(NCOUNT+8+KK)=DBLE(J1)*RR
          Y(NCOUNT+8+KK)=DBLE(J2)*RW
          Z(NCOUNT+8+KK)=DBLE(J3)*RU
 
          X(NCOUNT+16+KK)=DBLE(J1)*RU
          Y(NCOUNT+16+KK)=DBLE(J2)*RR
          Z(NCOUNT+16+KK)=DBLE(J3)*RW
 
          X(NCOUNT+24+KK)=DBLE(J1)*RU
          Y(NCOUNT+24+KK)=DBLE(J2)*RW
          Z(NCOUNT+24+KK)=DBLE(J3)*RR

          X(NCOUNT+32+KK)=DBLE(J1)*RW
          Y(NCOUNT+32+KK)=DBLE(J2)*RU
          Z(NCOUNT+32+KK)=DBLE(J3)*RR

          X(NCOUNT+40+KK)=DBLE(J1)*RW
          Y(NCOUNT+40+KK)=DBLE(J2)*RR
          Z(NCOUNT+40+KK)=DBLE(J3)*RU
 
          WT(NCOUNT+KK)=XD(K+1)
          WT(NCOUNT+KK+8)=XD(K+1)
          WT(NCOUNT+KK+16)=XD(K+1)
          WT(NCOUNT+KK+24)=XD(K+1)
          WT(NCOUNT+KK+32)=XD(K+1)
          WT(NCOUNT+KK+40)=XD(K+1)
         ENDDO
        ENDDO
       ENDDO
      NCOUNT=NCOUNT+48
      ENDDO
      NCOUNT=NCOUNT-1



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
       XDEV=DABS(1-DSQRT(x(I)**2+y(I)**2+z(I)**2))
       EABS=EABS+XDEV
       IF (XDEV.GT.1D-8) THEN
       WRITE(LUPRI,*) 'NODE ERROR(1) AT I',I,' deviation ',xdev
       CALL QUIT('Error 2 in subroutine_SLEB, get Aaron.')
       ENDIF
       IF (DABS(WT(I)).LT.1D-5) THEN
        WTERROR=0
        WRITE(LUPRI,*) 'BWGHT ERROR(1) AT I',I,' weight ',wt(I)
        CALL QUIT('Error 3 in subroutine_SLEB, get Aaron.')
       ENDIF
      ENDDO
      EABS=EABS/NCOUNT
      IF (EABS.LT.1D-15) EABS=1.0D-15


C***********************************************************************
C     Testing procedure on real and imaginary harmonics:l=[0,29],m=[0,l]
C     Do integration by quadrature in subroutine HARMONICS(l,m)
C     CALL HARMONICS



C*****INITIALISING all necessary values***************************

C     Check to see if the group posseses a threefold axis,
C     (excluding T,Td, O and Oh). i.e. use the subgroup of the
C     octahedron with the highest symmetry.

      DO J=1,NCOUNT
         XNEW(J)=X(J)
         YNEW(J)=Y(J)
         ZNEW(J)=Z(J)
         WTNEW(J)=WT(J)*4.0D0*PI
      ENDDO
      NNEQV=NCOUNT

C     Check values of nodes, and weights; calculate absolute errors.
C     WTERROR returns a value of 0 if any of the weights is zero.
C     EABS is the sum of all absolute deviations of point from
C     the surface of the unit sphere. It is then divided by NCOUNT.

      EABS=0.0D0
      WTERROR=0.0D0
      DO I=1,NNEQV
       XDEV=DABS(1-DSQRT(xnew(I)**2+ynew(I)**2+znew(I)**2))
       EABS=EABS+XDEV
       WTERROR=WTERROR+WTNEW(I)
       IF (XDEV.GT.1D-8) THEN
        WRITE(LUPRI,*) 'NODE ERROR(2) AT I',I,' deviation ',xdev
        CALL QUIT('NODE ERROR')
       ENDIF
      ENDDO
      XDEV=DABS(WTERROR/4.0D0/PI-1.0D0)
      IF (XDEV.GT.1D-9) THEN
       WRITE(LUPRI,*) 'SUM OF WEIGHTS NOT EQUAL TO 1',XDEV
       CALL QUIT('SUM OF WEIGHTS IN LEBEDEV NOT EQUAL TO 1')
      ENDIF
      EABS=EABS/NCOUNT
      IF (EABS.LT.1D-15) THEN
       EABS=1.0D-15
       ENDIF

      NSCHEME=NOLD
      RETURN
      END
C  /* Deck bweght */
      SUBROUTINE BWGHT(RJ,PSMU,RIJ,AIJ,APASC,XPASC,X,Y,Z,WT,
     &                 ACCUM,XMUIJN,XMUIJ2,IVECL,IPT,NTRANS,NATOM)
C
C      Written by C. W. Murray
C      BWGHT calculates the weights associated with
C      the Becke partitioning amongst the atoms.
C
#include <implicit.h> 
      PARAMETER (MAXBFN=1000)
      PARAMETER (MAXCEN=100)
      DIMENSION RJ(IVECL,NAT),
     &          PSMU(IVECL,NAT),
     &          RIJ(NAT*(NAT-1)/2),
     &          AIJ(NAT*(NAT-1)/2),
     &          X(IVECL),Y(IVECL),Z(IVECL),WT(IVECL),
     &          XMUIJN(IVECL),XMUIJ2(IVECL),ACCUM(IVECL),XPASC(20)
C
      COMMON /INFOA/ NAT,NUM,ZAN(MAXCEN),C(3,MAXCEN)
C
C     WORK OUT THE WEIGHT FUNCTION BY BECKE PARTITIONING
C
      DO 100 INA=1,NAT
      DO 100 M=1,IPT
         RJ(M,INA)=DSQRT((C(1,INA)-X(M))**2
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
C /* Deck mltpnt */
      FUNCTION MLTPNT(PX,PY,PZ)
C
C     Symmetry multiplicity of a point in space
C
C     T. Helgaker  Feb 01
C
#include <implicit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <maxaqn.h>
      PARAMETER (A0 = 1.0D-8)
      LOGICAL LBTAX
#include <symmet.h>
#include <ibtfun.h>
      LBTAX(I,J) = IBTAND(2**(I-1),ISYMAX(J,1)) .GT. 0
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
