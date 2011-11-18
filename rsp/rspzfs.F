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
C  /* Deck zfsdrv */
      SUBROUTINE ZFSDRV(WORK,LWORK)
C
C CALCULATE AVERAGE VALUE OF PROPERTIES
C
C
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <dummy.h>
C
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 , D8 = 8.0D0 )
      PARAMETER (DP5=0.5D0, D1P5=1.5D0)
      PARAMETER ( D1INF = 0.99999D0 , CKMXPR = 1.0D-6 )
      DIMENSION ISYMDM(2), IFCTYP(2)
      CHARACTER SPD(7)
      DATA SPD/'S','P','D','F','G','H','I'/
C
#include <maxorb.h>
#include <maxash.h>
#include <inforb.h>
#include <infrsp.h>
#include <infind.h>
#include <wrkrsp.h>
#include "rspprp.h"
#include "infave.h"
#include <infpri.h>
#include <inftap.h>
#include <iratdef.h>
#include <aovec.h>
#include <codata.h>
#include <gfac.h>
#include <priunit.h>
#include <blocks.h>
#include <zfs.h>
#include "infesr.h"
#include <infinp.h>
C
      LOGICAL ANTI, PANTI, NODPTR, NODV, NOPV, NOCONT, TTIME, 
     &   RETUR, NOBLK, DEBUG, DIA2SO, ZFS2EL
      CHARACTER*5 STRING
      IF (.NOT.ZFSCAL) RETURN

      CALL QENTER('ZFSDRV')
C
      STRING="     "
      IPRINT = IPRESR
      NODPTR = IPRINT.GT.10
      NOBLK = .FALSE.
      KFREE = 1
      LFREE = LWORK
      CALL HEADER("Output from ZFSDRV",0)

C
C     *******************************************************
C     ***** Set up COMMON /BLOCKS/ for PSORT and TWOINT *****
C     *******************************************************
C
      CALL MEMGET('INTE',KJSTRS,MXSHEL*MXAOVC*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNPRIM,MXSHEL*MXAOVC*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCONT,MXSHEL*MXAOVC*2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIORBS,MXSHEL*MXAOVC  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJORBS,MXSHEL*MXAOVC  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KKORBS,MXSHEL*MXAOVC  ,WORK,KFREE,LFREE)
C
      IPRPAO=0
      CALL PAOVEC(WORK(KJSTRS),WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            WORK(KJORBS),WORK(KKORBS),0,NOBLK,IPRPAO)
C
      CALL MEMREL('AVEDIA:PAOVEC',WORK,KJORBS,KJORBS,KFREE,LFREE)
C
      IF (HSROHF) THEN
         NODV = .FALSE.
         NOPV = .TRUE.
         NDMAT = 2
      ELSE
         NODV = .TRUE.
         NOPV = NASHT .LE. 1
         NDMAT = 0
      END IF
C
      IF (.NOT.NOPV) THEN
C
C     Transform two-electron density matrix to SO basis
C
         ANTI   = .FALSE.
         PANTI  = .FALSE.
         DIA2SO = .FALSE.
         ZFS2EL = .TRUE.
         KDTAO  = KFREE
         JPRINT = IPRINT
         CALL GPOPEN(LUPAO,'ZFSPAO','NEW',' ',' ',IDUMMY,.FALSE.)
         CALL PTRAN(NODPTR,WORK(KFREE),LFREE,JPRINT,ANTI,PANTI,DIA2SO,
     &      ZFS2EL)
         CALL PSORG(WORK(KFREE),WORK(KFREE),LFREE,WORK(KNCONT),JPRINT,
     &   ANTI,PANTI)
C
      ELSE IF (HSROHF) THEN
         CALL MEMGET('REAL',KDMAT,NNASHX,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KDTAO,2*N2BASX,WORK,KFREE,LFREE)
         KDVAO=KDTAO+N2BASX
         DO I=1,NASHT
            II= KDMAT + I*(I+1)/2 - 1
            WORK(II) = D1
         END DO
         CALL DZERO(WORK(KDTAO),N2BASX)
         CALL GETDMT(WORK(KDTAO),NDMAT,WORK(KFREE),LFREE,.TRUE.,.FALSE.,
     &      .TRUE.,IPRINT)
      ELSE
         CALL QUIT(
     &      'Zero-field splitting requires at least two open shells')
      END IF
C
C     Call HERMIT to evaluate expectation value. A lot of these variables 
C     may be of interest to control through a input routine.
C
      ISYMDM(1) = 1
      ISYMDM(2) = 1
      IFCTYP(1) = 14
      IFCTYP(2) = 14
      ITYPE     = 11
      MAXDIF    = 2
      JATOM     = 0
      NOCONT    = .FALSE.
      TTIME     = .FALSE.
      RETUR     = .FALSE.
      IPRNTA    = 0
      IPRNTB    = 0
      IPRNTC    = 0
      IPRNTD    = 0
      I2TYP     = 0
      CALL TWOINT(WORK(KFREE),LFREE,WORK(KFREE),WORK(KFREE),WORK(KDTAO),
     &            NDMAT,IREPDM,IFCTYP,DUMMY,IDUMMY,IDUMMY,1,ITYPE,
     &            MAXDIF,JATOM,NODV,NOPV,NOCONT,TTIME,IPRINT,IPRNTA,
     &            IPRNTB,IPRNTC,IPRNTD,RETUR,IDUMMY,I2TYP,WORK(KJSTRS),
     &            WORK(KNPRIM),WORK(KNCONT),WORK(KIORBS),
     &            IDUMMY,IDUMMY,DUMMY,DUMMY,DUMMY,
     &            DUMMY,.FALSE.)
C
      IF (.NOT. NOPV)  CALL GPCLOSE(LUPAO,'DELETE')
      CALL MEMREL('ZFSDRV',WORK,1,1,KFREE,LFREE)
C
C Analyze and print
C
      CALL ZFSANA(WORK,LWORK,IPRINT)
      CALL QEXIT('ZFSDRV')
      RETURN
      END
C  /* Deck zfsana */
      SUBROUTINE ZFSANA(WORK,LWORK,IPRINT)
#include <implicit.h>
      DIMENSION WORK(LWORK)
#include <maxorb.h>
#include <infinp.h>
      INTEGER S
      CALL QENTER('ZFSANA')
C
C Assuming determinant basis and high spin projection
C
      MULT=ISPIN
      MULT2=MULT*MULT


      KFREE=1
      LFREE=LWORK
      CALL MEMGET('COMP',KHZFS,MULT2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRZFS,MULT2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCG,MULT2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KHEIG,MULT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KRWORK,3*MULT-2,WORK,KFREE,LFREE)
      CALL ZFSAN1(MULT,WORK(KHZFS),WORK(KRZFS),WORK(KCG),WORK(KHEIG),
     &   WORK(KRWORK),WORK,KFREE,LFREE,IPRINT)
      CALL MEMREL('ZFSANA',WORK,1,1,KFREE,LFREE)

      CALL QEXIT('ZFSANA')
      END
C  /* Deck zfsan1 */
      SUBROUTINE ZFSAN1(MULT,HZFS,RZFS,CG,HEIG,RWORK,WORK,KFREE,LFREE,
     &   IPRINT)
#include <implicit.h>
      DOUBLE COMPLEX HZFS(MULT,MULT)
      DIMENSION RZFS(MULT,MULT)
      DIMENSION CG(MULT,MULT)
      DIMENSION HEIG(MULT)
      DIMENSION WORK(*), RWORK(*)
#include <priunit.h>
#include <zfs.h>
#include <codata.h>
#include <gfac.h>
      DOUBLE COMPLEX U(-2:2)
      DIMENSION RU(-2:2)
      DIMENSION ZFSEIG(3)
      DOUBLE COMPLEX I
      PARAMETER ( D1 = 1.0D0 , D2 = 2.0D0 , D6=6.0D0, D8 = 8.0D0 )
      PARAMETER (DP5=0.5D0, D1P5=1.5D0, D0=0.0D0)
      LOGICAL TESTZ 
      DATA TESTZ /.FALSE./
      
      CALL QENTER('ZFSAN1')

CBS   I=(D0,D1)
      I=(0D0,1D0)
      S=DBLE(MULT-1)/2

      IF (IPRINT.GT.5) THEN
         CALL HEADER('2-electron field gradient tensor',0)
         CALL OUTPUT(ZFS,1,3,1,3,3,3,1,LUPRI)
      END IF
C
C    Quintet density needs proper scaling
C
      DENFAC=D1/(3*S*S - S*(S+1))
C
C  Traceless part
C
      TR3=(ZFS(1,1)+ZFS(2,2)+ZFS(3,3))/3
      DO K=1,3
         ZFS(K,K)=ZFS(K,K)-TR3
      END DO
C
C  All in cm-1
C
      ZFSFAC=-DENFAC*XTKAYS*(GFAC/2)**2*ALPHAC**2
      CALL DSCAL(9,ZFSFAC,ZFS,1)
      CALL HEADER('Trace-less zero field splitting tensor (cm-1)',0)
      CALL OUTPUT(ZFS,1,3,1,3,3,3,1,LUPRI)
C
C  We now have a pure second rank tensor, get the spherical representation
C
      U2R=(ZFS(1,1)-ZFS(2,2))/2
      U2I=(ZFS(1,2)+ZFS(2,1))/2
      U1R=(ZFS(1,3)+ZFS(3,1))/2
      U1I=(ZFS(2,3)+ZFS(3,2))/2
      U( 2) = U2R + I*U2I
      U(-2) = DCONJG(U(2))
      U( 1) =-U1R - I*U1I
      U(-1) =-DCONJG(U(1))
      U( 0) = (2*ZFS(3,3)-ZFS(1,1)-ZFS(2,2))/SQRT(D6)
      IF (IPRINT.GT.5) THEN
         CALL HEADER('Spherical zero field splitting tensor (cm-1)',0)
         CALL COUTPUT(U,1,5,1,1,5,1,1,LUPRI)
      END IF
C
C Tabulate Clebsh Gordan CG_S2S(M,N) = <S2NM-N|SM>
C

      CALL DZERO(CG,MULT*MULT)
      DO IM=1,MULT
         RM = -S-1+IM
         DO IN=1,MULT
            RN = -S-1+IN
            MN=(IM-IN)
            IF (MN.EQ.2) THEN
               CG(IM,IN)=SQRT(
     &            (3*(S+RM-1)*(S+RM)*(S-RM+1)*(S-RM+2))/
     &            ((2*S-1)*2*S*(S+1)*(2*S+3))
     &            )
            ELSE IF (MN.EQ.1) THEN
               CG(IM,IN)=(1-2*RM)
     &           * SQRT(
     &               DBLE(3*(S-RM+1)*(S+RM)) /
     &               ((2*S-1)*S*(2*S+2)*(2*S+3))
     &               )
            ELSE IF (MN .EQ. 0) THEN
               CG(IM,IN)=DBLE(3*RM*RM-S*(S+1))/
     &            SQRT( DBLE((2*S-1)*S*(S+1)*(2*S+3)) )
            ELSE IF (MN.EQ.-1) THEN
               CG(IM,IN)=(2*RM+1)
     &           * SQRT(
     &               DBLE(3*(S-RM)*(S+RM+1))/
     &               ((2*S-1)*S*(2*S+2)*(2*S+3))
     &               )
            ELSE IF (MN.EQ.-2) THEN
               CG(IM,IN)=SQRT(
     &            DBLE(3*(S-RM-1)*(S-RM)*(S+RM+1)*(S+RM+2))/
     &            ((2*S-1)*S*(2*S+2)*(2*S+3))
     &            )
            END IF
         END DO
      END DO
      IF (IPRINT.GT.5) THEN
         CALL HEADER('Clebsh Gordan coefficients',0)
         CALL OUTPUT(CG,1,MULT,1,MULT,MULT,MULT,1,LUPRI)
      END IF
C
C Reduced matrix element of rank 2 evaluated as
C
C     TSS=<SS|T(2,0)|SS>/<S2S0|SS>
C
      TSS=SQRT( (2*S-1) * S*(S+1) * (2*S+3)/6 )
C
C  Build hamiltonian matrix
C
C     H(M,N)=<SM|T|SN>U(M-N)* = <S||T||S><S2NM-N|SM>U(M-N)*

      DO IM=1,MULT
         DO IN=1,MULT
            HZFS(IM,IN) = TSS*CG(IM,IN)*DCONJG(U(IM-IN))
         END DO
      END DO
      IF (IPRINT.GT.5) THEN
         CALL HEADER('ZFS Hamiltonian (spherical basis)',0)
         CALL COUTPUT(HZFS,1,MULT,1,MULT,MULT,MULT,1,LUPRI)
      END IF
C
C  Diagonalize HZFS for energy splittings
C
      CALL ZHEEV('N','U',MULT,HZFS,MULT,HEIG,WORK(KFREE),LFREE,
     &   RWORK,INFO)
      IF (INFO.LT.0) THEN
         WRITE(LUERR)'ZFSDRV:CSYEV:INFO=',INFO
         CALL QUIT('ZFSDRV:LAPACK DIAGONALIZATION CSYEV FAILED')
      END IF
      CALL HEADER('ZFS energy eigenvalues (cm-1)',0)
      CALL OUTPUT(HEIG,1,MULT,1,1,MULT,1,1,LUPRI)
C
C Conventions for triplet states
C
      IF (MULT.EQ.3) THEN
C
C Z the level of greatest splitting, X mid level
C
         DZ1=ABS(HEIG(1))
         DZ2=ABS(HEIG(2))
         DZ3=ABS(HEIG(3))
         IF (DZ1.GE.DZ3) THEN
            IZ=1
            IX=2
            IY=3
         ELSE
            IZ=3
            IX=2
            IY=1
         END IF
         DX=HEIG(IX)
         DY=HEIG(IY)
         DZ=HEIG(IZ)
C
C  D and E values for triplet reference states
C
         D=-D1P5*DZ
         IF (D.GT.0) THEN
            E=ABS(DX-DY)/2
         ELSE
            E=-ABS(DX-DY)/2
         END IF
         CALL HEADER('Zero field splitting paramters',0)
         WRITE(LUPRI,'(A,F10.6,A,F10.2,A)')'@ZFS parameter D = ',
     &   D, ' cm-1 = ',D*XTHZ*1D-6/XTKAYS,' MHz'
         WRITE(LUPRI,'(A,F10.6,A,F10.2,A)')'@ZFS parameter E = ',
     &   E, ' cm-1 = ',E*XTHZ*1D-6/XTKAYS,' MHz'
C
C Check convention (ref Langhoff...)
C
         IF (ABS(D).LT.3*ABS(E) .OR. D*E.LT.0) THEN
            WRITE(LUPRI,*) 'WARNING:ZFS Principle axis test failed'
            NWARN=NWARN+1
         END IF
C
C Get principal axes by diagonalizing D
C
         CALL DSYEV('V','U',MULT,ZFS,3,ZFSEIG,WORK(KFREE),LFREE,
     &         INFO)
         IF (INFO.LT.0) THEN
            WRITE(LUERR)'ZFSDRV:DSYEV:INFO=',INFO
            CALL QUIT('ZFSDRV:LAPACK DIAGONALIZATION DSYEV FAILED')
         END IF
         CALL HEADER('ZFS tensor eigenvalues and cosines',0)
         CM2MHZ=XTHZ*1D-6/XTKAYS
         WRITE(LUPRI,'(2A16,A24)')'cm-1','MHz','    direction cosines'
         DO K=1,3
            WRITE(LUPRI,'(2F16.6,3F8.4)') ZFSEIG(K),ZFSEIG(K)*CM2MHZ,
     &         (ZFS(J,K),J=1,3) 
         END DO
         IF (TESTZ) THEN
C
C Diagonalize D first -> U, H will be real
C

            RU(2)=(ZFSEIG(1)-ZFSEIG(2))/2
            RU(-2)=RU(2)
            RU(1)=D0
            RU(-1)=D0
            RU( 0) = (2*ZFSEIG(3)-ZFSEIG(1)-ZFSEIG(2))/SQRT(D6)
C
C  Build hamiltonian matrix
C
C     H(M,N)=<SM|T|SN>U(M-N)* = <S||T||S><S2NM-N|SM>U(M-N)*

            DO IM=1,MULT
               DO IN=1,MULT
                  RZFS(IM,IN) = TSS*CG(IM,IN)*RU(IM-IN)
               END DO
            END DO
            CALL HEADER('RZFS Hamiltonian (spherical basis)',0)
            CALL OUTPUT(RZFS,1,MULT,1,MULT,MULT,MULT,1,LUPRI)
            CALL DSYEV('N','U',MULT,RZFS,MULT,ZFSEIG,WORK(KFREE),LFREE,
     &         INFO)
            IF (INFO.LT.0) THEN
               WRITE(LUERR)'ZFSDRV:DSYEV:INFO=',INFO
               CALL QUIT('ZFSDRV:LAPACK DIAGONALIZATION DSYEV FAILED')
            END IF
            CALL HEADER('ZFS energy eigenvalues (cm-1)',0)
            CALL OUTPUT(ZFSEIG,1,MULT,1,1,MULT,1,1,LUPRI)
         END IF
      END IF
C
C Clean up and exit
C
      CALL QEXIT('ZFSAN1')
      END
