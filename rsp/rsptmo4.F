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
      SUBROUTINE CRTMO(VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD,
     *                 RESVEC,CMO,UDV,PV,FOCK,FC,FV,
     *                 XINDX,MJWOP,WRK,LWRK)
C
#include <implicit.h>
#include <iratdef.h>
#include <dummy.h>
C
C PURPOSE:
C CALCULATION OF THIRD ORDER TRANSITION MOMENTS
C
      LOGICAL DOMOM, DIPLEN, ATEST
C
      CHARACTER*8 ALAB,BLAB,CLAB,DLAB
C
      DIMENSION VECA(*),VECB(*),VECC(*),VECD(*)
      DIMENSION VECBC(*),VECCD(*),VECBD(*)
      DIMENSION RESVEC(*)
      DIMENSION CMO(*),UDV(*),PV(*),FOCK(*),FC(*)
      DIMENSION XINDX(*),WRK(*)
C
      PARAMETER ( D0 = 0.0D0, ZEROTHR = 1.0D-10, IMXEXC = 20 )
      DIMENSION RESTOM(3,3,3,IMXEXC,8), RESFRE(IMXEXC,8)
      DIMENSION IEQTO(5)
C
#include <priunit.h>
#include <infrsp.h>
#include <maxorb.h>
#include <infvar.h>
#include <inforb.h>
#include <infpri.h>
#include <infspi.h>
#include <wrkrsp.h>
#include <tstjep.h>
#include <infhso.h>
#include <qrinf.h>
      DIMENSION MJWOP(2,MAXWOP,8)
#include <rspprp.h>
#include <indcr.h>
#include <infcr.h>
#include <inftmo.h>
#include <codata.h>
C
      CALL QENTER('CRTMO')
      ATEST = .FALSE.
C
C     We write the results to a complementary output file. This will then
C     both serve as a file for getting a summary of the results, but more
C     importantly, it will serve as a way of avoiding already calculated
C     gamma-components in case of a crashed calculation. Check for calculated
C     components are done in BCDCHK.
C
      LURSPRES = -1
      CALL GPOPEN(LURSPRES,'RESULTS.RSP','UNKNOWN',' ','FORMATTED',
     &            IDUMMY,.FALSE.)
      WRITE(LUPRI,'(//A,A)')
     *' <<<<< CALCULATING CONTRIBUTIONS TO THIRD ORDER',
     *' TRANSITION MOMENT >>>>>'
C
      IF (THREEPHOTON) THEN
        CALL DZERO(RESTOM,3*3*3*IMXEXC*8)
        CALL DZERO(RESFRE,IMXEXC*8)
      END IF
C
      DO 200 ISYMD = 1,NSYM
      DO 300 ISYMC = 1,NSYM
      DO 400 ISYMB = 1,NSYM
C
      ISYMDX = MULD2H(IREFSY,ISYMD)
      ISYMA = MULD2H(ISYMD,MULD2H(ISYMC,ISYMB))
      IF ( (NTMCNV(ISYMD).GT.0) .AND. (NCTMOP(ISYMC).GT.0) .AND.
     *     (NBTMOP(ISYMB).GT.0) .AND. (NATMOP(ISYMA).GT.0) ) THEN
C
      DO 500 ID   = 1,NTMCNV(ISYMD)
C
      DO 600 ICOP = 1,NCTMOP(ISYMC)
      DO 650 ICFR = 1,NCTMFR
C
      DO 700 IBOP = 1,NBTMOP(ISYMB)
      DO 750 IBFR = 1,NBTMFR
C
C     If three photon absorption or harmonic generation calculation 
C     only certain frequencies are considered.
C
      IF ( THREEPHOTON .AND.
     &     ((ABS(EXCIT2(ISYMD,ID)-3*BTMFR(IBFR)).GT.ZEROTHR) .OR.
     &      (ABS(EXCIT2(ISYMD,ID)-3*CTMFR(ICFR)).GT.ZEROTHR)) ) GOTO 750
C
      IF (CTMOHG) THEN
         IF (IBFR.NE.ICFR) GO TO 750
      END IF
C
      DO 800 IAOP = 1,NATMOP(ISYMA)
C
C     Initialize variables.
C     Check if an equivalent moment calculation already has been done,
C     DOMOM indicates the result.
C     Read response vectors and eigen vectors from disk.
C     Check if some of the response vectors are equal or zero,
C     IBCDEQ indicates the result
C
      CALL BCDCHK(DOMOM,IBCDEQ,LURSPRES,DIPLEN,DUMMY,
     *            ISYMA,ISYMB,ISYMC,ISYMD,ISYMBC,ISYMBD,ISYMCD,
     *            ALAB,BLAB,CLAB,DLAB,IAOP,IBOP,ICOP,0,
     *            IBFR,ICFR,ID,FREQA,FREQB,FREQC,FREQD,
     *            KZYVA,KZYVB,KZYVC,KZYVD,KZYVBC,KZYVBD,KZYVCD,
     *            VECA,VECB,VECC,VECD,VECBC,VECBD,VECCD,IEQTO)
C
      IF (THREEPHOTON .AND. (.NOT. DOMOM)) THEN
        CALL DIPLAB(ALAB,IA)
        CALL DIPLAB(BLAB,IB)
        CALL DIPLAB(CLAB,IC)
        RESFRE(ID,ISYMD)=EXCIT2(ISYMD,ID)
        RESTOM(IA,IB,IC,ID,ISYMD)=
     &  RESTOM(IEQTO(1),IEQTO(2),IEQTO(3),IEQTO(4),IEQTO(5))
      END IF
      IF (.NOT.DOMOM) GOTO 800
C
C    Initialize third order transition moment
C
      TMOM = 0
C
      IF (IPRRSP.GT.0) WRITE(LUPRI,'(//A15,2A20,/A)')
     *   'Contribution','Term','Accumulated',
     *   ' ------------------------------------------------------'
C
C
C     Calculate Na T[4] Nb Nc Nd
C
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL T4DRV(IBCDEQ,ISYMA,ISYMB,ISYMC,ISYMD,VECA,VECB,VECC,VECD,
     *           -FREQB,-FREQC,-FREQD,XINDX,UDV,PV,MJWOP,
     *           WRK,LWRK,CMO,FC)
      VAL = -DDOT(KZYVA,WRK,1,VECA,1)
      TMOM = TMOM + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A17,F18.8,F20.8)')' Na T[4] Nb Nc Nd',VAL,TMOM
C
      IF (IPRRSP.GT.5) CALL TIMER('T4DRV ',TIMSTR,TIMEND)
C
C
C     Calculate Na T[3] Nb Ncd type terms (three permutations)
C
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL T3DRV(1,ISYMA,ISYMB,ISYMCD,VECB,VECCD,ATEST,VECA,
     *           -FREQB,-FREQC-FREQD,XINDX,UDV,PV,MJWOP,
     &           WRK,LWRK,CMO,FC,FV)
      VAL = DDOT(KZYVA,WRK,1,VECA,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      IF (IBCDEQ.EQ.2) THEN
         TMPVAL = TMPVAL + VAL
         TMOM = TMOM + VAL
         CALL T3DRV(1,ISYMA,ISYMD,ISYMBC,VECD,VECBC,ATEST,VECA,
     *              -FREQD,-FREQB-FREQC,XINDX,UDV,PV,MJWOP,
     &              WRK,LWRK,CMO,FC,FV)
         VAL = DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         TMOM = TMOM + VAL
      ELSE
         CALL T3DRV(1,ISYMA,ISYMC,ISYMBD,VECC,VECBD,ATEST,VECA,
     *              -FREQC,-FREQB-FREQD,XINDX,UDV,PV,MJWOP,
     *              WRK,LWRK,CMO,FC,FV)
         VAL = DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         TMOM = TMOM + VAL
         CALL T3DRV(1,ISYMA,ISYMD,ISYMBC,VECD,VECBC,ATEST,VECA,
     *              -FREQD,-FREQB-FREQC,XINDX,UDV,PV,MJWOP,
     *              WRK,LWRK,CMO,FC,FV)
         VAL = DDOT(KZYVA,WRK,1,VECA,1)
         TMPVAL = TMPVAL + VAL
         TMOM = TMOM + VAL
      END IF
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)')' Na T[3] Nx Nyz',TMPVAL,TMOM
C
      IF (IPRRSP.GT.5) CALL TIMER('T3DRV ',TIMSTR,TIMEND)

C
C     Calculate Na B[3] Nc Nd type terms 
C     (two of four permutations in each call)
C
C
      IF (IPRRSP.GT.5) CALL TIMER('START ',TIMSTR,TIMEND)
C
      CALL X3INIT(KZYVA,KZYVC,KZYVD,ISYMA,ISYMC,ISYMD,BLAB,
     *            ISYMB,VECC,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      CALL X3INIT(KZYVA,KZYVB,KZYVD,ISYMA,ISYMB,ISYMD,CLAB,
     *            ISYMC,VECB,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVA,RESVEC,1,VECA,1)
      TMOM = TMOM + VAL
      TMPVAL = TMPVAL + VAL
C
      IF (IPRRSP.GT.0)
     *  WRITE(LUPRI,'(A15,3F20.8)')' Na X[3] Ny Nz ',TMPVAL,TMOM
C
C
C     Calculate Nb A[3] Nc Nd type terms
C     (two of six permutations in each call)
C
C
      CALL A3INIT(KZYVB,KZYVC,KZYVD,ISYMB,ISYMC,ISYMD,ALAB,
     *            ISYMA,VECC,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVB,RESVEC,1,VECB,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      CALL A3INIT(KZYVC,KZYVB,KZYVD,ISYMC,ISYMB,ISYMD,ALAB,
     *            ISYMA,VECB,VECD,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVC,RESVEC,1,VECC,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A3INIT(KZYVD,KZYVB,KZYVC,ISYMD,ISYMB,ISYMC,ALAB,
     *            ISYMA,VECB,VECC,RESVEC,XINDX,UDV,CMO,MJWOP,WRK,LWRK)
      VAL = DDOT(KZYVD,RESVEC,1,VECD,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      IF (IPRRSP.GT.0)
     *  WRITE(LUPRI,'(A15,2F20.8)')' Nx A[3] Ny Nz ',TMPVAL,TMOM
C
C
C     Calculate Na B[2] Ncd type terms (two permutations)
C
C
      CALL X2INIT(1,KZYVA,KZYVCD,ISYMA,ISPINA,ISYMCD,0,WRK(1),VECCD,
     *            RESVEC,XINDX,UDV,PV,BLAB,ISYMB,0,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      CALL X2INIT(1,KZYVA,KZYVBD,ISYMA,ISPINA,ISYMBD,0,1,VECBD,
     *            RESVEC,XINDX,UDV,PV,CLAB,ISYMC,ISPINC,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVA,RESVEC,1,VECA,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)') ' Na X[2] Nyz   ',TMPVAL,TMOM
C
C
C     Calculate Nb A[2] Ncd type terms (six permutations)
C
C
      CALL A2INIT(1,KZYVB,KZYVCD,ISYMB,ISPINB,ISYMCD,0,1,VECCD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVB,RESVEC,1,VECB,1)
      TMPVAL = VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVCD,KZYVB,ISYMCD,0,ISYMB,ISPINB,1,VECB,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVCD,RESVEC,1,VECCD,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVC,KZYVBD,ISYMC,ISPINC,ISYMBD,0,1,VECBD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVC,RESVEC,1,VECC,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVBD,KZYVC,ISYMBD,0,ISYMC,ISPINC,1,VECC,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVBD,RESVEC,1,VECBD,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVD,KZYVBC,ISYMD,ISPIND,ISYMBC,0,1,VECBC,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVD,RESVEC,1,VECD,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      CALL A2INIT(1,KZYVBC,KZYVD,ISYMBC,0,ISYMD,ISPIND,1,VECD,
     *            RESVEC,XINDX,UDV,PV,ALAB,ISYMA,ISPINA,
     *            CMO,MJWOP,WRK,LWRK)
      VAL = -DDOT(KZYVBC,RESVEC,1,VECBC,1)
      TMPVAL = TMPVAL + VAL
      TMOM = TMOM + VAL
C
      IF (IPRRSP.GT.0)
     * WRITE(LUPRI,'(A15,2F20.8)') ' Nx A[2] Nyz   ',TMPVAL,TMOM
C
      IF (IPRRSP.GT.5) CALL TIMER('OTHERS',TIMSTR,TIMEND)
C
      WRITE(LUPRI,'(/A,3(/A,A10,I4,F10.6))')
     *     '@ Third order transition moment in a.u. for',
     *     '@ A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *     '@ B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB,
     *     '@ C operator, symmetry, frequency: ',CLAB,ISYMC,-FREQC
      WRITE(LUPRI,'(/A,2I4,F10.6//A,F20.8/)')
     *     '@ State no., symmetry, excitation energy:',ID,ISYMD,FREQD,
     &     '@ < 0 | ABC | f >  = ', TMOM
C
C     Write out to result file
C
      WRITE(LURSPRES,'(/A,3(/A,A10,I4,F10.6))')
     *     ' Third order transition moment in a.u. for',
     *     ' A operator, symmetry, frequency: ',ALAB,ISYMA,FREQA,
     *     ' B operator, symmetry, frequency: ',BLAB,ISYMB,-FREQB,
     *     ' C operator, symmetry, frequency: ',CLAB,ISYMC,-FREQC
      WRITE(LURSPRES,'(/A,2I4,F10.6)')
     *     ' State no., symmetry, excitation energy:',ID,ISYMD,FREQD
      WRITE(LURSPRES,'(/A,F20.8)') ' < 0 | ABC | f >  = ', TMOM
C
      IF (THREEPHOTON) THEN
        CALL DIPLAB(ALAB,IA)
        CALL DIPLAB(BLAB,IB)
        CALL DIPLAB(CLAB,IC)
        RESTOM(IA,IB,IC,ID,ISYMD)=TMOM
        RESFRE(ID,ISYMD)=EXCIT2(ISYMD,ID)
      END IF
C
 800  CONTINUE
 750  CONTINUE
 700  CONTINUE
 650  CONTINUE
 600  CONTINUE
 500  CONTINUE
C
      END IF
C
 400  CONTINUE
 300  CONTINUE
 200  CONTINUE
      CALL GPCLOSE(LURSPRES,'KEEP')
C
C
      IF (THREEPHOTON) THEN
         CALL PRINT_THREEPHOTON(RESTOM,RESFRE)   
      END IF
C
C    End of subroutine CRTMO
C
      CALL QEXIT('CRTMO')
      RETURN
      END
      SUBROUTINE PRINT_THREEPHOTON(RESTOM,RESFRE)
C
#include <implicit.h>
#include <priunit.h>
#include <inforb.h>
#include <rspprp.h>
#include <inftmo.h>
#include <codata.h>
C
      PARAMETER ( D0 = 0.0D0, ZERO = 1.0D-10, IMXEXC = 20 )
      DIMENSION RESTOM(3,3,3,IMXEXC,8), RESFRE(IMXEXC,8)
      CHARACTER*8 ALAB,BLAB,CLAB,DLAB
      CHARACTER*1 ITINDX,JTINDX
C
      CALL TITLER('FINAL RESULTS FROM THREE-PHOTON CALCULATION',
     &       '*',112)
      WRITE(LUPRI,'(A64,4(/A64))')
     & ' The three-photon absorption strength for an average molecular',
     & ' orientation is computed according to formulas given by       ',
     & ' W.M. McClain in J. Chem. Phys. 57:2264, 1972. The absorption ',
     & ' depends on the light polarization. A monochromatic light     ',
     & ' source is assumed.                                           '
      WRITE(LUPRI,'(2(/A64))')
     & ' All results are presented in atomic units, except the        ',
     & ' excitation energy which is given in eV.                      '
      WRITE(LUPRI,'(/,3(/A52))')
     &        '+----------------------------------+',
     &        '| Three-photon transition tensor T |',
     &        '+----------------------------------+'
      WRITE(LUPRI,'(2A)')
     & ' ----------------------------------------------------',
     & '--------------'
C
      DO ISYMD=1,NSYM
      DO ID=1,NTMCNV(ISYMD)
C
        DF=D0
        DG=D0
        WRITE(LUPRI,'(A12,I2,A9,I2,A11,F8.3,A3)') '   Symmetry ',ISYMD,
     &  '   State ',ID,'   Energy: ',RESFRE(ID,ISYMD)*XTEV,' eV'
        WRITE(LUPRI,'(A1)')
        DO I=1,3
        DO J=1,3
            CALL LABDIP(ITINDX,I)
            CALL LABDIP(JTINDX,J)
          WRITE(LUPRI,'(3(A5,2A1,A3,F12.2))')
     &   '  T_{',ITINDX,JTINDX,'x}=',RESTOM(I,J,1,ID,ISYMD),
     &   '  T_{',ITINDX,JTINDX,'y}=',RESTOM(I,J,2,ID,ISYMD),
     &   '  T_{',ITINDX,JTINDX,'z}=',RESTOM(I,J,3,ID,ISYMD)
        END DO
        END DO
      WRITE(LUPRI,'(2A)')
     & ' ----------------------------------------------------',
     & '--------------'
      END DO
      END DO
      WRITE(LUPRI,'(/,6(/A60))')
     &        ' Transition probabilities (a.u.)         ',
     &        '--------------------------------------------------',
     &        ' D  =  (3*Df + 2*Dg)/35, Linear   polarization',
     &        ' D  = (-3*Df + 5*Dg)/35, Circular polarization',
     &        ' Df = sum(i,j,k){ T_iij * T_kkj }             ',
     &        ' Dg = sum(i,j,k){ T_ijk * T_ijk }             '
      WRITE(LUPRI,'(3(/A53))')
     &        '       Polarization ratio      ',
     &        '-------------------------------',
     &        '    R  = (-3*Df+5*Dg)/(3*Df+2*Dg)  '
      WRITE(LUPRI,'(/,3(/A56))')
     &        '+-------------------------------------+',
     &        '| Three-photon transition probability |',
     &        '+-------------------------------------+'
      WRITE(LUPRI,'(A62,A)')
     &'-----------------------------------------------------------',
     &'-----------'
      WRITE(LUPRI,'(A6,A4,A8,A14,3A11,A8/A62,A)') 'Sym','No',
     &        'Energy','Polarization','Df','Dg','D','R',
     &'-----------------------------------------------------------',
     &'-----------'
      DO ISYMD=1,NSYM
      DO ID = 1,NTMCNV(ISYMD)
         IF (ID.LE.IMXEXC) THEN
            DF=0.0D0
            DG=0.0D0
            DO I=1,3
            DO J=1,3
            DO K=1,3
               DF=DF+RESTOM(I,I,J,ID,ISYMD)*RESTOM(K,K,J,ID,ISYMD)
               DG=DG+RESTOM(I,J,K,ID,ISYMD)**2
            END DO
            END DO
            END DO
            R=(-3*DF+5*DG)/(3*DF+2*DG)
            D=(3*DF+2*DG)/35
            WRITE(LUPRI,'(A2,2I4,F8.2,A14,3E11.3,F8.2)')' ',ISYMD,ID,
     &           RESFRE(ID,ISYMD)*XTEV,'Linear     ',DF,DG,D,R
            D=(-3*DF+5*DG)/35
            WRITE(LUPRI,'(A2,2I4,F8.2,A14,3E11.3,F8.2)')' ',ISYMD,ID,
     &              RESFRE(ID,ISYMD)*XTEV,'Circular   ',DF,DG,D,R
         ELSE
            WRITE(LUPRI,'(/A)')
     &           'WARNING!! NOT ENOUGH SPACE TO IN RESULT MATRIX'
         END IF
      END DO
      END DO
      WRITE(LUPRI,'(A62,A)')
     &'-----------------------------------------------------------',
     &'-----------'
      RETURN
      END
      SUBROUTINE LABDIP(LAB,I)
C
#include <implicit.h>
C
      CHARACTER*1 LAB
C
C Map integer numbers to dipole operators
C 1 => x, 2 => y, and 3 => z
C
      IF (I.EQ.1) LAB(1:1)='x'
      IF (I.EQ.2) LAB(1:1)='y'
      IF (I.EQ.3) LAB(1:1)='z'
C
      RETURN
      END
