#define LDA
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRDFT(EXCMAT,DMAT,EDFTY,DOERG,DOATR,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Mar. 2003
C
C     Purpose : Calculate SR DFT-energy and potential contribution for 
C               SR-DFT hybrid methods.
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <dftinf.h>
#include <dftcom.h>
C
      LOGICAL DOGGA,DOLND,DOGRD,DOERG,DOATR
      DIMENSION EXCMAT(*),DMAT(*),WORK(LWORK)
C
!      WRITE(LUPRI,*) 'SRDFT: entering'  !JT

      CALL QENTER('SRDFT')
      CALL TIMER('START ',TIMSTR,TIMEND)
C     =================================
C     Initializations :
C     Number of AOs and their addresses
C     =================================
      DFTHR0 = 1.0D-8 
      DFTHRL = 1.0D-10
      DFTHRI = 1.0D-13
      DOGRD  = .FALSE.
      DOLND  = .FALSE.

C     ... jkp : For gradient correction to SR-LDA we will
C               be calling The Becke and LYP functionals.
      DOGGA  = DOSRXGGA.OR.DOSRCGGA.OR.DOSRBCK.OR.DOSRGGA2.OR.
     >         DOSRXWIB.OR.DOSRCWIB.OR.DOSRXPBEHSE.OR.DOSRXPBETCS.OR.
     >         DOSRCPBETCS.OR.DOSRXRI.OR.DOSRCRI ! Manu 01-02-2006 
c

      NDER  = 0
      IF (DOGRD) NDER = 1
      IF (SRHYBR.AND..NOT.DOHFSRDFT) THEN
         NDER = NDER + 2
C        ... we use DFTPOT for SRHYBR
      ELSE IF (DOGGA) THEN
         NDER = NDER + 2
C        ... we always use DFTPOT in this version
CC       NDER = NDER + 1
CC       IF (DFTPOT) NDER = NDER + 1
      END IF
C
C


      IF (NDER.EQ.0) NTYPSO =  1
      IF (NDER.EQ.1) NTYPSO =  4
      IF (NDER.EQ.2) NTYPSO = 10
      IF (NDER.GT.2) CALL QUIT('NDER.gt.2 not implemented')
      NSO0 = 1
      NSO1 = 2
      NSO2 = 5
      IF (DOLND) THEN
         NTYPSO = NTYPSO + 3
         NSOB   = NTYPSO - 2
         IF (DOGGA) THEN
            NTYPSO = NTYPSO + 9
            NSOB1  = NTYPSO - 8
         END IF
      END IF
      KSO0 = (NSO0-1)*NBAST + 1
      KSO1 = (NSO1-1)*NBAST + 1
      KSO2 = (NSO2-1)*NBAST + 1
      KSOB = (NSOB-1)*NBAST + 1
      KSOB1= (NSOB1-1)*NBAST + 1
C     =============================
C     Allocations
C     =============================
      NBUF   = 1000000
      KX     = 1
      KY     = KX    + NBUF
      KZ     = KY    + NBUF
      KW     = KZ    + NBUF
      KGSO   = KW    + NBUF
      KCNT   = KGSO  + NTYPSO*NBAST
      KDGASR = KCNT  + NBAST
      IF (SRHYBR.AND. .NOT.DOHFSRDFT) THEN
         KDGAT = KDGASR + NBAST
      ELSE
         KDGAT = KDGASR
      END IF
      KWRK  = KDGAT  + NBAST
      LWRK  = LWORK + 1 - KWRK 
      IF (KWRK .GT. LWORK) CALL ERRWRK('SRDFT',-KWRK,LWORK)

C JT 15-02-05 beg
      IF(ISJT) THEN
      CALL SRDFT1JT(EXCMAT,EDFTY,DMAT,DOLND,WORK(KX),WORK(KY),
     &            WORK(KZ),WORK(KW),NBUF,WORK(KGSO),WORK(KCNT),
     &            WORK(KDGASR),WORK(KDGAT),DOGGA,DOERG,DOATR,
     &            WORK(KWRK),LWRK,IPRINT)
C JT 15-02-05 end
      ELSE  !JT 15-02-05
      CALL SRDFT1(EXCMAT,EDFTY,DMAT,DOLND,WORK(KX),WORK(KY),
     &            WORK(KZ),WORK(KW),NBUF,WORK(KGSO),WORK(KCNT),
     &            WORK(KDGASR),WORK(KDGAT),DOGGA,DOERG,DOATR,
     &            WORK(KWRK),LWRK,IPRINT)
      ENDIF  !JT 15-02-05

C
      CALL TIMER('SRDFT ',TIMSTR,TIMEND)
      CALL QEXIT('SRDFT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRDFT1(EXCMAT,EDFTY,DMAT,DOLND,CORX,CORY,CORZ,
     &                  WEIGHT,NBUF,GSO,NCNT,DMAGAOSR,DMAGAOT,DOGGA,
     &                  DOERG,DOATR,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Mar. 2003
C
C     Purpose : Calculate DFT-energy and potential contribution for 
C               SR-DFT hybrid methods.
C
C*****************************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pi.h>
#include <dummy.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0, D1 = 1.0D0, D2 = 2.0D0)
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infpar.h>
#include <nuclei.h>
#include <dftinf.h>
#include <dftcom.h>
#include <gnrinf.h>
C
      LOGICAL   DOLND,DOGGA,DOERG,DOATR, LDOGGA,LSRHYBR,RTEST,ldebug
      DIMENSION EXCMAT(NBAST,NBAST,*),DMAT(N2BASX,*),CORX(NBUF),
     &          CORY(NBUF),CORZ(NBUF),WEIGHT(NBUF),GSO(NBAST*NTYPSO),
     &          NCNT(NBAST),DMAGAOSR(NBAST),DMAGAOT(NBAST),WORK(LWORK),
     &          RHGT(3),RHGSR(3),ZKSR(3,3),XKSR(3,3)
      CHARACTER*16 QUADNAME
C
#include <chrnos.h>
#include <sdpre.h>
C
!      WRITE(LUPRI,*) 'SRDFT1: entering'  !JT
!      WRITE(LUPRI,*) 'SRDFT1: DOGGA=',DOGGA   !JT
!      WRITE(LUPRI,*) 'SRDFT1: ERFEXP=',ERFEXP !JT


      ldebug = .false.
      LSRHYBR = SRHYBR .AND. .NOT.DOHFSRDFT
      IF (LSRHYBR) THEN
         NXCMAT = 2
      ELSE
         NXCMAT = 1
      END IF
      IPNTstart = 8168
c     IPNTend   = 10210
      IPNTend   = 0
#ifdef jkp_debug
      WRITE(LUPRI,'(//A//)') ' <<< Total EXCMAT on entry >>>'
      CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      WRITE(LUPRI,'(//A//)') ' <<< Total density matrix >>>'
      CALL OUTPUT(DMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      IF (LSRHYBR) THEN
      WRITE(LUPRI,'(//A//)') ' <<< Valence density matrix >>>'
      CALL OUTPUT(DMAT(1,2),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      IF (DOATR) THEN
      WRITE(LUPRI,'(//A//)') ' <<< Total DTRMAT >>>'
      CALL OUTPUT(DMAT(1,3),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      WRITE(LUPRI,'(//A//)') ' <<< Valence DTRMAT >>>'
      CALL OUTPUT(DMAT(1,4),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
      END IF
#endif
      QUADNAME='                '
C     =============================
C     Calculate grid
C     =============================
C
      IF (.NOT.GRDONE) THEN
         CALL DFTGRD(WORK,LWORK,NTOTAL,1,.FALSE.)
         CALL CONDFT
         GRDONE = .TRUE.
      END IF
C
      LDOGGA = DOGGA .OR. LSRHYBR
C
      ELCTRN = D0
      ELCTSR = D0
      ELCTRSKIP = D0
C     ... Initialize energies to zero. They should come out 
C         zero if NULL functional is defined.
      EX     = D0
      EC     = D0
      EXC    = D0
      ETX    = D0
      ETC    = D0
      VXC    = D0
      VTX    = D0
      VTC    = D0
      VXBT   = D0
      VXBV   = D0
      VTOT   = D0
C
C     For SRHYBR:
C
      EXSR   = D0
      EXTOT  = D0
      ECSR   = D0
      ECTOT  = D0
C
C     HJAaJ: THRELCTRN makes sure that error in
C            no. of electrons from grid integration less than DFTHR0
C            (NTOTAL is the total number of grid points)
C
c1    THRELCTRN = 1.0D-6 / NTOTAL
      THRELCTOT = 1.0D-7
      THRELCTRN = THRELCTOT/100
c     THRDRHO   = DFTHRI
      THRDRHO   = 1.0D-10
C Debug values: include everything
c#ifdef DFT_NOSCREEN
      THRELCTRN = D0
      THRDRHO   = D0
      DFTHR0    = D0
      DFTHRI    = D0
      DFTHRL    = D0
c#endif
C
C     HJAaJ: skip eval. of neighbor points with same RHO
C
      RHOOLD = D0
      NSKIP1 = 0
      NSKIP2 = 0
C
      LUQUAD = -1
#if defined (VAR_MPI) 
C             
C     Make quadname: Can take 10000 procs 
C
      QUADNAME = 'DALTON.QUAD.'//chrnos(mynum/1000)
     &     //chrnos((mynum-(mynum/1000)*1000)/100)
     &     //chrnos((mynum-(mynum/100)*100)/10)
     &     //chrnos(mynum-(mynum/10)*10)
C
      CALL GPOPEN(LUQUAD,QUADNAME,'OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,LDUMMY)   
#else
      CALL GPOPEN(LUQUAD,'DALTON.QUAD','OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
#endif
      NPNTS = 0
  200 CONTINUE
      READ(LUQUAD) NPOINT
      kpoint = npoint/10
      IF (NPOINT.GT.0) THEN
D         jpnt = 0
         NPNTS = NPNTS + NPOINT
         CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
#if defined (VAR_MPI)
         DO 300 IPNT = 1+MYNUM, NPOINT, NODTOT+1
#else
         DO 300 IPNT = 1, NPOINT
#endif
            IF (IPRINT .GT. 100) THEN
               WRITE (LUPRI,'(2X,I6,4F12.6)') 
     &            IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
            END IF
C
            WGHT = WEIGHT(IPNT)
C
C           AOs
C           ===
C
            THRINT = DFTHRI/WGHT
            CALL GETSOS(GSO,NCNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &                  WORK,LWORK,NBAST,DOLND,LDOGGA,THRINT,IPRINT)
C           =============================
C           Density
C           =============================
            IF (LSRHYBR) THEN
               CALL GETRHO_OLD(DMAT(1,2),GSO(KSO0),RHOSR,RHOSR13,
     &                     DMAGAOSR,THRINT,IPRINT)
               CALL GETRHO_OLD(DMAT(1,1),GSO(KSO0),RHOTOT,RHOT13,
     &                     DMAGAOT,THRINT,IPRINT)
            ELSE
               CALL GETRHO_OLD(DMAT(1,1),GSO(KSO0),RHOSR,RHOSR13,
     &                     DMAGAOSR,THRINT,IPRINT)
               RHOTOT = RHOSR
               RHOT13 = RHOSR13
               CALL DCOPY(NBAST,DMAGAOSR,1,DMAGAOT,1)
            END IF
            DELCTRN = WGHT*RHOTOT
            DELCTSR = WGHT*RHOSR
Chj            IF (RHOTOT.LE.DFTHR0) THEN
            IF (ABS(DELCTRN).LE.THRELCTRN) THEN
               NSKIP1 = NSKIP1 + 1
               ELCTRSKIP = ELCTRSKIP + DELCTRN
               THRELCTRN = (THRELCTOT-ELCTRSKIP)/100
            ELSE
C
               RHOTST = RHOTOT
               IF (LDOGGA) THEN
C             ... For the gradient correction to SR-LDA we need
C                 the gradient of the density for the energy and
C                 the hessian of the density for the potential.
C
C                 Gradient of density
C                 ===================
                  CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMAGAOSR,
     &                         1,D0,RHGSR,1)
                  RHOSRGRD = SQRT(RHGSR(1)**2 + RHGSR(2)**2 
     &                     + RHGSR(3)**2)
                  IF (LSRHYBR) THEN
                     CALL DGEMV('T',NBAST,3,D2,GSO(KSO1),NBAST,DMAGAOT,
     &                          1,D0,RHGT,1)
                     RHOTGRD = SQRT(RHGT(1)**2 + RHGT(2)**2  
     &                       + RHGT(3)**2)
                  ELSE
                     RHGT(1) = RHGSR(1)
                     RHGT(2) = RHGSR(2)
                     RHGT(3) = RHGSR(3)
                     RHOTGRD = RHOSRGRD
                  ENDIF
C
C                 Hessian of density
C                 ===================
C
                  IF (LSRHYBR) THEN
                     CALL DFTRHH(DMAT(1,2),DMAGAOSR,GSO(KSO0),GSO(KSO1),
     &                           GSO(KSO2),RHGSR,RHOSRLAP,RHOSRGHG)
                     CALL DFTRHH(DMAT,DMAGAOT,GSO(KSO0),GSO(KSO1),
     &                           GSO(KSO2),RHGT,RHOTLAP,RHOTGHG)
                  ELSE
                     CALL DFTRHH(DMAT,DMAGAOSR,GSO(KSO0),GSO(KSO1),
     &                           GSO(KSO2),RHGSR,RHOSRLAP,RHOSRGHG)
                     RHOTLAP = RHOSRLAP
                     RHOTGHG = RHOSRGHG
                  END IF
                  RHOTST = RHOTST + RHOTGRD + RHOTLAP + RHOTGHG
               ENDIF
C              =============================
C              Number of electrons
C              =============================
               ELCTRN = ELCTRN + DELCTRN
               ELCTSR = ELCTSR + DELCTSR
C              =============================
C              ETX,ETC: Energy
C              VTX,VTC: Exchange-correlation potential
C              ==============================
               IF (DOERG) THEN
c              RTEST = (ABS(RHOTST-RHOOLD) .GT. THRDRHO)
c              IF (RTEST) THEN
               IF (.true.) THEN
C              ... hjaaj 10. July 2003: we could make this better by using that
C                  EXC(RHO) = EXC(RHOOLD+DRHO) ~ EXC(RHOOLD) + VXC(RHOOLD)*DRHO
C                  = ETX+ETC + (VTX+VTC)*DRHO (as EXC and VXC are calculated for RHOOLD)
                  RHOOLD = RHOTST
C                 ... only reevaluate if new RHOTOT /hjaaj July 2003
c
c                 -------------------------------------------------------
c                 Exchange Correlation energy and potential
c                 -------------------------------------------------------
c
                  if (IPNT.GE.IPNTstart.AND.IPNT.LE.IPNTend) then
                     ldebug=.true.
                  else
                     ldebug=.false.
                  end if
                  CALL SRDFTEXC(XKSR,ZKSR,RHOSR,RHOSR13,RHOSRGRD,
     &                          RHOSRLAP,RHOSRGHG,CHIVAL,ERFEXP,
     &                         .TRUE.,VXB,ldebug)
               ELSE
C                  we can use ETX, ETC, VTX, VTC from last IPNT
                   NSKIP2 = NSKIP2 + 1
               END IF
               EX   = EX + WGHT*XKSR(1,1)
               EC   = EC + WGHT*ZKSR(1,1)
               VXC  = WGHT*(XKSR(2,1) + ZKSR(2,1))
c
c             ... SRHYBR (RHOTOT) contribution
c              IF (RTEST) THEN
               IF (.true.) THEN
                  IF (LSRHYBR) THEN
                     VXB   = -WGHT*VXB
                     EXSR  = EXSR + WGHT*XKSR(1,2)
                     ECSR  = ECSR + WGHT*ZKSR(1,2)
                     VXC   = VXC - WGHT*(XKSR(2,2)+ZKSR(2,2))
                     CALL SRDFTEXC(XKSR,ZKSR,RHOTOT,RHOT13,RHOTGRD,
     &                             RHOTLAP,RHOTGHG,CHIVAL,ERFEXP,
     &                            .FALSE.,VXB,ldebug)
                     EXTOT = EXTOT + WGHT*XKSR(1,2)
                     ECTOT = ECTOT + WGHT*ZKSR(1,2)
                     VXCT  = WGHT*(XKSR(2,2)+ZKSR(2,2))
                     VXBT  = WGHT*VXB
                  ENDIF
               ENDIF
#ifdef debug
C               if (mod(IPNT,kpoint) .eq. 0) then
               if (IPNT.GE.IPNTSTART.AND.IPNT.LE.IPNTEND) then
                  write (lupri,*) 'IPNT,ex,ec,vx,vc',IPNT,ex,ec,vx,vc
                  write (lupri,*) 'wght,xksr,zksr',wght,xksr(1),zksr(1)
                  write (lupri,*) 
     &            'rhotot,rhot13,rhotgrd,rhotlap,rhotghg : ',
     &            rhotot,rhot13,rhotgrd,rhotlap,rhotghg
                  write (lupri,*) 
     &            'rhosr,rhosr13,rhosrgrd,rhosrlap,rhosrghg : ',
     &            rhosr,rhosr13,rhosrgrd,rhosrlap,rhosrghg
                  write (lupri,*) 'nskip1,nskip2',nskip1,nskip2
                  write (lupri,*) 'excmat(1,1,1),gso(1) before: ',
     &                             excmat(1,1,1),gso(1)
               end if
#endif
C              =====================================================
C              Exchange-correlation contribution to Kohn-Sham matrix
C              =====================================================
C -hjaaj: implement screening>             IF (ABS(VXC) .GT. SCRTHR)
C
C  --- NOTE ! This uses FROMVX.
!      WRITE(LUPRI,*) 'SRDFT1JT: IPNT=',IPNT,' VXC=',VXC !JT

               CALL DFTKSM(EXCMAT,GSO(KSO0),GSO(KSO1),RHGSR,VXC,VXB,
     &                     LDOGGA,.true.,DFTHRL)
               IF (LSRHYBR) THEN
c             ... SRHYBR (RHOTOT) contribution
                  CALL DFTKSM(EXCMAT(1,1,2),GSO(KSO0),GSO(KSO1),RHGT,
     &                        VXCT,VXBT,LDOGGA,.true.,DFTHRL)
               END IF
C
c#ifdef debug
c               if (mod(IPNT,kpoint) .eq. 0) then
               IF (IPNT.GE.IPNTSTART.AND.IPNT.LE.IPNTEND) then
                  write (lupri,*) 'excmat(1,1,1),gso(1) after: ',
     &                             excmat(1,1,1),gso(1)
               endif
c#endif
               ENDIF
C
C              Hessian transformation
C              ======================
C
               IF (DOATR) THEN
                  IF (LSRHYBR) THEN
C                    matrix 1 & 2 are DTAO and DVAO
C                    matrix 3 & 4 are DXTAO and DXVAO
                     JDXAO = 4
                     CALL SRDFTLTR(JWOPSY,DMAT(1,JDXAO),EXCMAT,
     &                             WGHT, GSO(KSO0),GSO(KSO1),CHIVAL,
     &                             ERFEXP,.TRUE.,RHOSR,RHOSR13,
     &                             RHOSRGRD,RHOSRLAP, RHOSRGHG,RHGSR,
     &                             DMAGAOSR)
                     JDXAO = 3
                     CALL SRDFTLTR(JWOPSY,DMAT(1,JDXAO),EXCMAT(1,1,2),
     &                             WGHT, GSO(KSO0),GSO(KSO1),CHIVAL,
     &                             ERFEXP,.FALSE.,RHOTOT,RHOT13,RHOTGRD,
     &                             RHOTLAP, RHOTGHG,RHGT,DMAGAOT)
                  ELSE
C                    matrix 1 is DTAO
C                    matrix 2 is DXTAO
                     JDXAO = 2
                     CALL SRDFTLTR(JWOPSY,DMAT(1,JDXAO),EXCMAT,
     &                             WGHT, GSO(KSO0),GSO(KSO1),CHIVAL,
     &                             ERFEXP,.TRUE.,RHOTOT,RHOT13,
     &                             RHOTGRD,RHOTLAP,RHOTGHG,RHGT,DMAGAOT)
                  END IF
               END IF
            ENDIF
C
  300    CONTINUE
C
         GO TO 200
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 200
      END IF
C
      CALL GPCLOSE(LUQUAD,'KEEP')
#if defined  (VAR_MPI)
      call quit('MPI code not correct, finish coding ...')
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
C
      CALL MPI_REDUCE(ELCTRN,EELCTRN,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR)
      CALL MPI_REDUCE(EC,EEC,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR)
      CALL MPI_REDUCE(EX,EEX,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR)
      CALL MPI_REDUCE(EXCMAT,WORK,NBAST*NBAST*NXCMAT,
     &     MPI_DOUBLE_PRECISION,MPI_SUM,
     %     0, MPI_COMM_WORLD, IERR)
C     ... Collect on master
      IF (MYNUM .EQ. 0) then
         ELCTRN = EELCTRN
         EC     = EEC
         EX     = EEX
         CALL DCOPY(NBAST*NBAST,WORK,1,EXCMAT,1)
      ENDIF
C     ... Release nodes
      IF (MYNUM.NE.0) RETURN
#endif
C
C
      IF (DOERG) THEN
         EDFTY = EX + EC
         IF (LSRHYBR) THEN
            EDFTY = EDFTY + EXTOT - EXSR + ECTOT - ECSR
         END IF
         FAC = DP5
      ELSE IF (DOATR) THEN
         FAC = D1
      ELSE
         CALL QUIT('both DOERG and DOATR are false')
      END IF
      DO K = 1, NXCMAT
      DO I = 1, NBAST
      DO J = 1, I - 1
         AVERAG = FAC*(EXCMAT(I,J,K) + EXCMAT(J,I,K))
         EXCMAT(I,J,K) = AVERAG
         EXCMAT(J,I,K) = AVERAG
      END DO
      END DO
      END DO
C
C     Print section
C
!     IF (IPRINT .GE. 3) THEN
      IF (DOERG) THEN
         ELCTOT = ELCTRN + ELCTRSKIP
         WRITE (LUPRI,'(3(/1X,A,F18.10))')
     &'No. of active electrons from numerical integration  :',ELCTRN,
     &'No. of active electrons skipped in num. integration :',ELCTRSKIP,
     &'Sum of active and skipped electrons                 :',ELCTOT
         IF (LSRHYBR) THEN
            WRITE (LUPRI,'(1X,A,F18.10)')
     & 'No. of valence electrons  (SRHYBR)                  :',ELCTSR
         WRITE (LUPRI,'(/1X,A,F18.10)')
     &'Exchange energy (SRHYBR: SR valence part)           :',EX
         WRITE (LUPRI,'(1X,A,F18.10)')
     &'Correlation energy  (SRHYBR: SR valence part)       :',EC,
     &'Exchange energy (SRHYBR: BLYP valence part)         :',EXSR,
     &'Correlation energy  (SRHYBR: BLYP valence part)     :',ECSR,
     &'Exchange energy (SRHYBR: BLYP total)                :',EXTOT,
     &'Correlation energy  (SRHYBR: BLYP total)            :',ECTOT
            EX = EX + EXTOT - EXSR
            EC = EC + ECTOT - ECSR
         END IF
         WRITE (LUPRI,'(/1X,A,F18.10)')
     &'Exchange energy                                     :',EX
         WRITE (LUPRI,'(1X,A,F18.10)')
     &'Correlation energy                                  :',EC
         WRITE (LUPRI,'(1X,A,F18.10)')
     &'DFT exchange-correlation energy                     :',EX+EC
         IF (IPRINT.GE.5) THEN
            WRITE(LUPRI,'(//A//)') ' <<< Kohn-Sham matrix >>>'
            CALL OUTPUT(EXCMAT,1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
            IF (LSRHYBR) THEN
            WRITE(LUPRI,'(//A//)')
     &        ' <<< Valence SR - valence Kohn-Sham matrix >>>'
            CALL OUTPUT(EXCMAT(1,1,2),1,NBAST,1,NBAST,NBAST,NBAST,
     &        1,LUPRI)
            END IF
         END IF
      END IF
!     END IF
C
      RETURN 
      END
C /* Deck srdftexc */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SRDFTEXC(XKSR,ZKSR,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG,
     &                    CHIVAL,ERFEXP,LSRLDA,VXB,ldebug)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Dec. 2003
C
C     Purpose : Driver for Short-range DFT exchange and correlation
C               contributions to energe and potential.
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxorb.h>
#include <inforb.h>
#include <infinp.h>
#include <dftinf.h>
#include <dftcom.h>
      LOGICAL ERFEXP,LSRLDA,ldebug
      DIMENSION XKSR(3,2), ZKSR(3,2), ZMU(3)
c
      VB    = 0.0D0
      VL    = 0.0D0
      VXB   = 0.0D0
      CALL DZERO(XKSR,6)
      CALL DZERO(ZKSR,6)
      CALL DZERO(ZMU ,3)
c
c     -------------------------------------------------------
c     Exchange energy and potential
c     -------------------------------------------------------
c
      IF (DOSRXLDA.OR.(SRHYBR.AND.LSRLDA)) THEN
          CALL VXSRLDA(XKSR,RHO,RHO13,CHIVAL,1,ERFEXP)
      END IF
      IF (DOSRXGGA) THEN
          CALL VXSRGGA(XKSR,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG,
     &                 CHIVAL,ERFEXP)
      END IF 
      IF (DOSRGGA2) THEN
          CALL VXSRGGA2(XKSR,RHO,RHO13,RHOGRD,RHOLAP,
     &                  RHOGHG,CHIVAL,ERFEXP)
      END IF
      IF (DOSRXWIB) THEN
          CALL EBCK(XKSR(1,1),RHO,RHO13,RHOGRD)
          CALL VBCK(XKSR(2,1),RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
          IF (IWINT.EQ.1) THEN
             CALL VXSRLDA(ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
             CALL VSRWINT(XKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE IF (IWINT.EQ.2.OR.IWINT.EQ.3) THEN
             CALL VXSRLDA(ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
             CALL EDRC(ED,RHO,RHO13)
             CALL VDRC(VD,RHO13)
             XKSR(1,1) = XKSR(1,1) + ED
             XKSR(2,1) = XKSR(2,1) + VD
             CALL VSRWINT2(XKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE IF (IWINT.GT.10) THEN
             CALL EDRC(ED,RHO,RHO13)
             CALL VDRC(VD,RHO13)
             XKSR(1,1) = XKSR(1,1) + ED
             XKSR(2,1) = XKSR(2,1) + VD
             CALL VSRWINT2(XKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE
            CALL QUIT ('Wrong dft functional interpolation scheme'//
     >                 'in DFTEXC (IWINT).')
          END IF
      END IF
      IF (DOSRBCK ) THEN
          CALL EDRC(ED,RHO,RHO13)
          CALL EBCK(EB,RHO,RHO13,RHOGRD)
          XKSR(1,1) = ED + EB
          CALL VDRC(VD,RHO13)
          CALL VBCK(VB,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
          XKSR(2,1) = VD + VB
      ENDIF
      IF (SRHYBR) THEN
          CALL EDRC(ED,RHO,RHO13)
#ifndef LDA
          CALL EBCK(EB,RHO,RHO13,RHOGRD)
#else
          EB  = 0.0D0
#endif
          XKSR(1,2) = ED + EB
          CALL VDRC(VD, RHO13)
          CALL GBCK(VB,VXB,RHO,RHO13,RHOGRD)
#ifdef LDA
          VB   = 0.0D0
          VXB  = 0.0D0
#endif
          XKSR(2,2) = VD  + VB
      END IF
c
c     -------------------------------------------------------
c     Correlation energy and potential
c     -------------------------------------------------------
c
      IF ((DOSRCLDA.OR.DOSRMULO).OR.(SRHYBR.AND.LSRLDA)) THEN
          CALL VCSRLDA(ZKSR,RHO,RHO13,CHIVAL,DOSRMULO,1,ERFEXP)
      END IF
      IF (DOSRCGGA) THEN
          CALL VCSRGGA(ZKSR,RHO,RHO13,RHOGRD,RHOLAP,
     &                 CHIVAL,ERFEXP)
      END IF
      IF (DOSRCWIB) THEN
          CALL ELYP(ZKSR(1,1),RHO,RHO13,RHOGRD)
          CALL VLYP(ZKSR(2,1),RHO,RHO13,RHOGRD,RHOLAP)
          IF (IWINT.EQ.1) THEN
             CALL VCSRLDA(ZMU,RHO,RHO13,CHIVAL,DOSRMULO,1,ERFEXP)
             CALL WVWN(VLDA,RHO,RHO13,ELDA,.TRUE.,.TRUE.)
             ZKSR(1,1) = ZKSR(1,1) - ELDA*RHO
             ZKSR(2,1) = ZKSR(2,1) - VLDA
             CALL VSRWINT(ZKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE IF (IWINT.GE.2.OR.IWINT.EQ.3) THEN
             CALL VCSRLDA(ZMU,RHO,RHO13,CHIVAL,DOSRMULO,1,ERFEXP)
             CALL VSRWINT2(ZKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE IF (IWINT.GT.10) THEN
             CALL VSRWINT2(ZKSR,ZMU,RHO,RHO13,CHIVAL,1,ERFEXP)
          ELSE
             CALL QUIT ('Wrong dft functional interpolation scheme'//
     >                  'in DFTEXC (IWINT).')
          END IF
      END IF
      IF (DOSRLYPT) THEN
          CALL ELYP(ZKSR,RHO,RHO13,RHOGRD)
          CALL GLYP(ZKSR(2,1),VL2,RHO,RHO13,RHOGRD)
          VXBT = VL2
      ENDIF
      IF (SRHYBR) THEN
          CALL ELYP(EL,RHO,RHO13,RHOGRD)
#ifdef LDA
          CALL WVWN(VL,RHO,RHO13,EL,.TRUE.,.TRUE.)
          VL2 = 0.0D0
          EL  = EL*RHO
#endif
          ZKSR(1,2) = EL
#ifndef LDA
          CALL GLYP(VL,VL2,RHO,RHO13,RHOGRD)
#endif
          VXB  = VXB  + VL2
          ZKSR(2,2) = VL
      ENDIF
c
      RETURN
      END
#ifdef JKP_DEBUG
C /* Deck mcrho */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MCRHO(CMO,INDXCI,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Jan. 2003
C
C     Purpose : Generate active density from multiconfigurationel
C               wavefunction on a grid.
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
#include <dftinf.h>
#include <dftcom.h>
C
      LOGICAL DOGGA,DOBCK,DOLYP,DOLND,DOGRD,DODCAO
      DIMENSION CMO(*),WORK(LWORK)
C
      CALL TIMER('START ',TIMSTR,TIMEND)
C     =================================
C     Initializations :
C     Number of AOs and their addresses
C     =================================
      DOLND = .FALSE.
      DOBCK = .FALSE.
      DOLYP = .FALSE.
      DOGRD = .FALSE.
      DOGGA = DOBCK .OR. DOLYP
c
      NDER  = 0
      IF (DOGRD) NDER = 1
      IF (DOGGA) THEN
         NDER = NDER + 1
         IF (DFTPOT) NDER = NDER + 1
      END IF
C
      IF (NDER.EQ.0) NTYPSO =  1
      IF (NDER.EQ.1) NTYPSO =  4
      IF (NDER.EQ.2) NTYPSO = 10
      NSO0 = 1
      NSO1 = 2
      NSO2 = 5
      IF (DOLND) THEN
         NTYPSO = NTYPSO + 3
         NSOB   = NTYPSO - 2
         IF (DOGGA) THEN
            NTYPSO = NTYPSO + 9
            NSOB1  = NTYPSO - 8
         END IF
      END IF
      KSO0 = (NSO0-1)*NBAST + 1
      KSO1 = (NSO1-1)*NBAST + 1
      KSO2 = (NSO2-1)*NBAST + 1
      KSOB = (NSOB-1)*NBAST + 1
      KSOB1= (NSOB1-1)*NBAST + 1
C     =============================
C     Allocations
C     =============================
      DODCAO = (NISHT.GT.0)
      NBUF   = 1000000
      KCREF  = 1
      KX     = KCREF + NCONF
      KY     = KX    + NBUF
      KZ     = KY    + NBUF
      KW     = KZ    + NBUF
      KGSO   = KW    + NBUF
      KCNT   = KGSO  + NTYPSO*NBAST
      KDGA   = KCNT  + NBAST
      KDV    = KDGA  + NBAST
      KDVAO  = KDV   + NNASHX
      KDCAO  = KDVAO + N2BASX
      IF (DODCAO) THEN
          KWRK  = KDCAO  + N2BASX
      ELSE
          KWRK  = KDCAO
      ENDIF
      LWRK  = LWORK + 1 - KWRK 
      CALL MCRHO1(CMO,WORK(KDV),DOLND,WORK(KX),WORK(KY),
     &            WORK(KZ),WORK(KW),NBUF,WORK(KGSO),WORK(KCNT),
     &            WORK(KDGA),DOGGA,WORK(KCREF),INDXCI,DODCAO,
     &            WORK(KDCAO),WORK(KDVAO),WORK(KWRK),LWRK,IPRINT)
C
      CALL TIMER('MCRHO',TIMSTR,TIMEND)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE MCRHO1(CMO,DV,DOLND,CORX,CORY,CORZ,WEIGHT,NBUF,GSO,
     &                  NCNT,DMAGAO,DOGGA,CREF,INDXCI,DODCAO,DCAO,DVAO,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Written by Jesper Kielberg Pedersen, Jan. 2003
C
C     Purpose : Generate density from multiconfigurationel
C               wavefunction on a grid.
C
C*****************************************************************************
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pi.h>
#include <dummy.h>
      PARAMETER (D0 = 0.0D0, DP5 = 0.5D0)
C
#include <maxorb.h>
#include <infinp.h>
#include <inforb.h>
#include <infvar.h>
#include <inftap.h>
#include <infpar.h>
#include <nuclei.h>
#include <dftinf.h>
#include <dftcom.h>
C
      LOGICAL   DOLND,DOGGA,DODCAO
      DIMENSION CMO(*),DV(NNASHX),CORX(NBUF),CORY(NBUF),
     &          CORZ(NBUF),WEIGHT(NBUF),GSO(NBAST*NTYPSO),
     &          NCNT(NBAST),DMAGAO(NBAST),CREF(NCONF),
     &          DCAO(NBAST,NBAST),DVAO(NBAST,NBAST),WORK(LWORK)
      CHARACTER*16 QUADNAME
C
#include <chrnos.h>
#include <sdpre.h>
C
C     =============================
C     Open file for dumping density
C     =============================
      CALL GPOPEN(LUDDUMP,'DENSITY.GRID','UNKNOWN',
     &                ' ','FORMATTED',IDUMMY,.FALSE.) 
      QUADNAME='                '
C
C     =============================
C     Retrieve CREF
C     =============================
C     
      REWIND LUIT1
      CALL MOLLAB('STARTVEC ',LUIT1,LUERR)
      DO 100 I = 1, (ISTACI-1)
         READ (LUIT1)
  100 CONTINUE
      CALL READT(LUIT1,NCONF,CREF)
C     
C     =============================
C     Calculate DV
C     =============================
C
      CALL MAKDV(CREF,DV,INDXCI,WORK,LWORK)
C     
C     =============================
C     Calculate inactive and active
C     density matrises in AO-basis
C     =============================
C
      CALL FCKDEN(DODCAO,.TRUE.,DCAO,DVAO,CMO,DV,WORK,LWORK)
C
C     =============================
C     Calculate grid
C     =============================
C
      IF (.NOT.GRDONE) THEN
         CALL DFTGRD(WORK,LWORK,NTOTAL,1,.FALSE.)
         CALL CONDFT
         GRDONE = .TRUE.
      END IF
C
      CELCTRN = D0
      VELCTRN = D0
C
      LUQUAD = -1
#if defined (VAR_MPI)
C            
C     Make quadname: Can take 10000 procs 
C     
      QUADNAME = 'DALTON.QUAD.'//chrnos(mynum/1000)
     &     //chrnos((mynum-(mynum/1000)*1000)/100)
     &     //chrnos((mynum-(mynum/100)*100)/10)
     &     //chrnos(mynum-(mynum/10)*10)
C     
      CALL GPOPEN(LUQUAD,QUADNAME,'OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,LDUMMY)      
#else
      CALL GPOPEN(LUQUAD,'DALTON.QUAD','OLD','SEQUENTIAL',
     &     'UNFORMATTED',IDUMMY,.FALSE.)
#endif
      NPNTS = 0
  200 CONTINUE
      READ(LUQUAD) NPOINT
      IF (NPOINT.GT.0) THEN
         NPNTS = NPNTS + NPOINT
         CALL REAQUA_OLD(CORX,CORY,CORZ,WEIGHT,LUQUAD,NPOINT)
#if defined (VAR_MPI)
         DO 300 IPNT = 1+MYNUM, NPOINT, NODTOT+1 
#else 
         DO 300 IPNT = 1, NPOINT
#endif
            IF (IPRINT .GT. 100) THEN
               WRITE (LUPRI,'(2X,I6,4F12.6)') 
     &            IPNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),WEIGHT(IPNT)
            END IF
C
            WGHT = WEIGHT(IPNT)
C
C           AOs
C           ===
C
            THRINT = DFTHRI/WGHT
            CALL GETSOS(GSO,NCNT,CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &                  WORK,LWORK,NBAST,DOLND,DOGGA,THRINT,IPRINT)
C
C           Density & Number of electrons
C           =============================
C
C           ... inactive part
            IF (DODCAO) THEN
               CALL GETRHO_OLD(DCAO,GSO,RHO,RHO13,DMAGAO,THRINT,IPRINT)
               CELCTRN = CELCTRN + WGHT*RHO
            ENDIF
C           ... active part
            CALL GETRHO_OLD(DVAO,GSO,RHO,RHO13,DMAGAO,THRINT,IPRINT)
            VELCTRN = VELCTRN + WGHT*RHO
C
C           Dump active density to file
C           =============================
C
            WRITE(LUDDUMP,'(4F25.13)') CORX(IPNT),CORY(IPNT),
     &           CORZ(IPNT),RHO
C
  300    CONTINUE
C
         GO TO 200
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 200
      END IF
C
      CALL GPCLOSE(LUQUAD,'KEEP')
C
C     Test on the number of electrons
C
      ELCTRX = FLOAT(2*NRHFT)
      ELCTRN = CELCTRN + VELCTRN
      ERROR  = ELCTRN - ELCTRX
C
C     Print section
C
      WRITE (LUPRI,'(3(/2X,A,F14.6))')
     &' Nr. of inactive electrons from numerical integration :',CELCTRN,
     &' Nr. of active electrons from numerical integration   :',VELCTRN,
     &' Total nr. of electrons from numerical integration    :',ELCTRN
      WRITE (LUPRI,'(2(/2X,A,F14.6))')
     &' Number of electrons from orbial occupations          :',ELCTRX,
     &' Error in the number of electrons                     :',ERROR
      IF (ABS(ERROR) .GT. DFTELS) THEN
         WRITE (LUPRI,'(/2X,A,F14.6,/2X,A)')
     &' Error larger than DFTELS (set input)                 :',DFTELS,
     &   ' Calculation aborted.'
         CALL QUIT
     &    ('Wrong number of electrons in DFTDRV. Calculation aborted.')
      END IF
C
      IF (IPRINT .GT. 5) THEN 
         WRITE (LUPRI,'(/2X,A,F14.7,6X,D8.2,I14)')
     &      ' Number of electrons/abscissas:  ', 
     &        ELCTRN,ELCTRN-ELCTRX,NPNTS
      END IF
      RETURN 
      END
#endif  /* JKP_DEBUG */
C /* Deck srdftltr */
      SUBROUTINE SRDFTLTR(KSYMOP,DTRMAT,EXCMAT,WGHT,GAO,GAO1,CHIVAL,
     &                    ERFEXP,LSRLDA,RHO,RHO13,RHOGRD,RHOLAP,RHOGHG,
     &                    RHG,DTGAO)
C
C     Jesper K. Pedersen. Nov. 2003
C     Purpose : Generate the 1-index transformed dft Hessian
C               needed for short-range DFT hybrids.
C     Based on DFTLTR used in regular DFT.
C
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <pi.h>
C
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, DP5 = 0.5D0,
     >           D13=1.0D0/3.0D0,
     >          RSFAC = 0.62035 04908 99400 08660 D0)
C               RSFAC = (D3/(D4*PI))**THIRD
      LOGICAL INVGGA,DOGGA,ERFEXP,LSRLDA
C
#include <maxorb.h>
#include <inforb.h>
#include <infinp.h>
#include <nuclei.h>
#include <dftinf.h>
#include <dftcom.h>
C
      DIMENSION DTRMAT(NBAST,NBAST),GAO(NBAST),GAO1(NBAST,3),
     &          EXCMAT(NBAST,NBAST),RHG(3),DTGAO(NBAST),
     &          DTGTAO(NBAST),XKSR(3,3),ZKSR(3,3),ZMU(3),BT3(3),BV3(3)
      DIMENSION B3(3)
C
      FRZTMP = 0D0
      FZZTMP = 0D0
C
C     Check for valid functional
C
      INVGGA = DOSRXGGA.OR.DOSRBCK.OR.DOSRGGA2.OR.DOSRCGGA.OR.DOSRLYPT
      DOGGA  = DOSRXWIB.OR.DOSRCWIB.OR.SRHYBR
C
      CALL DGEMV('N',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
      B0 = DDOT(NBAST,DTGAO,1,GAO,1)
C
      IF (INVGGA) THEN
         CALL QUIT('SRDFTLTR called with GGA type functional. '//
     &             'This is not implemented!')
      END IF
C
C     ---------------------------------------------------------
C     SRLDA,MULOCAL,SRHYBR or Interpolated GGA/SRLDA functional
C     ---------------------------------------------------------
C
      IF (DOGGA) THEN
C        B3 = GAO1'*DTGAO
         CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D0,B3,1)
C        DTGAO= DTRMAT'*GAO
         CALL DGEMV('T',NBAST,NBAST,D1,DTRMAT,NBAST,GAO,1,D0,DTGAO,1)
C        B3 = B3 + GAO1'*DTGAO
         CALL DGEMV('T',NBAST,3,D1,GAO1,NBAST,DTGAO,1,D1,B3,1)
         BMAX = MAX(DABS(B0),DABS(B3(1)),DABS(B3(2)),DABS(B3(3)))
C
         RS   = RSFAC/RHO13
         X = CHIVAL*RS*COPFAC
         DERFCX = DERFC(X)
         GAUX=(D2/SQRTPI)*X*EXP(-D13*X**2)
      ELSE
         BMAX = DABS(B0)
         CALL DZERO(B3,3)
      END IF
C
      IF (BMAX.GT.DFTHRL) THEN 
C
C        *****************************
C        ** Calculate exchange part **
C        *****************************
C
         IF (.NOT.DOSRMULO) THEN
c
            IF (LSRLDA) THEN
c        ...  Dont do SRLDA if SRDFTLTR is called with LSRLDA false
c             (regular Kohm-Sham part with total density for SRHYBR).
              CALL VXSRLDA(XKSR,RHO,RHO13,CHIVAL,2,ERFEXP)
            END IF
c
c           Below are various corrections to the SRLDA energy/potential
c
            IF (DOSRXWIB) THEN
c         ... Special case for Weighted Interpolation functionals
              CALL DCOPY(3,XKSR,1,ZMU,1)
              CALL DZERO(XKSR,3)
              CALL EBCK(XKSR(1,1),RHO,RHO13,RHOGRD)
              CALL VBCK(XKSR(2,1),RHO,RHO13,RHOGRD,RHOLAP,RHOGHG)
              CALL V1BCK(FR0BCK,FZ0BCK,FRRBCK,FRZBCK,
     &                   FZZBCK,RHO,RHOGRD)
              XKSR(3,1) = FRRBCK
              IF (IWINT.EQ.1) THEN
                 CALL VSRWINT(XKSR,ZMU,RHO,RHO13,CHIVAL,2,ERFEXP)
              ELSE IF (IWINT.EQ.2.OR.IWINT.EQ.3.OR.IWINT.GT.10) THEN
                 IF (IWINT.GT.10) CALL DZERO(ZMU,3)
                 CALL EDRC(ED,RHO,RHO13)
                 CALL VDRC(VD,RHO13)
                 CALL V1DRC(VDRC0,VDRC1,RHO,RHO13)
                 XKSR(1,1) = XKSR(1,1) + ED
                 XKSR(2,1) = XKSR(2,1) + VD
                 XKSR(3,1) = XKSR(3,1) + VDRC1
                 CALL VSRWINT2(XKSR,ZMU,RHO,RHO13,CHIVAL,2,ERFEXP)
              ELSE
                 write(lupri,*) 'IWINT : ',IWINT
                 CALL QUIT ('Wrong dft functional interpolation '//
     &                      'scheme in SRDFTLTR (IWINT).')
              END IF
c             --- Remember mixed terms from Becke
              IF (IWINT.EQ.3) THEN
                 FRZTMP = FRZBCK*GAUX
                 FZZTMP = FZZBCK*GAUX
              ENDIF
              FRZBCK = FRZBCK*DERFCX + FRZTMP
              FZZBCK = FZZBCK*DERFCX + FZZTMP
c
            ELSE IF (SRHYBR) THEN
c          ... Special case for SRDFT in DFT embedding.
              CALL V1DRC(VDRC0,FRRDRC,RHO,RHO13)
#ifndef LDA
              CALL V1BCK(FR0BCK,FZ0BCK,FRRBCK,FRZBCK,FZZBCK,
     &                   RHO,RHOGRD)
              IF (LSRLDA) THEN
                 XKSR(3,1) =  XKSR(3,1) - FRRDRC - FRRBCK
                 FR0BCK    = -FR0BCK
                 FZ0BCK    = -FZ0BCK
                 FRRBCK    = -FRRBCK
                 FRZBCK    = -FRZBCK
                 FZZBCK    = -FZZBCK
              ELSE
                  XKSR(3,1) = FRRDRC + FRRBCK
              END IF
#else
              FR0BCK = 0D0
              FZ0BCK = 0D0
              FRRBCK = 0D0
              FRZBCK = 0D0
              FZZBCK = 0D0 
#endif
            ELSE
              FR0BCK = 0D0
              FZ0BCK = 0D0
              FRRBCK = 0D0
              FRZBCK = 0D0
              FZZBCK = 0D0
            ENDIF
         ENDIF
C
C        ********************************
C        ** Calculate correlation part **
C        ********************************
C
         IF (LSRLDA) THEN
c     ... Dont do SRLDA if SRDFTLTR is called with LSRHYBR true
c          (regular Kohm-Sham part with total density for SRHYBR).
           CALL VCSRLDA(ZKSR,RHO,RHO13,CHIVAL,DOSRMULO,2,ERFEXP)
         END IF
c
c     ... Below are various corrections to the SRLDA energy/potential
c
         IF (DOSRCWIB) THEN
c        ... Special case for Weighted Interpolation functionals
           CALL DCOPY(3,ZKSR,1,ZMU,1)
           CALL DZERO(XKSR,3)
c          -- lyp
           CALL ELYP(ZKSR(1,1),RHO,RHO13,RHOGRD)
           CALL VLYP(ZKSR(2,1),RHO,RHO13,RHOGRD,RHOLAP)
           RHOA = DP5*RHO
           RHGA = (DP5*RHOGRD)**2
           CALL GLYPCO(DF1000,DF0100,DF0010,DF0001,
     &                 DF00001,RHO,RHO13,RHOGRD,.TRUE.)
           CALL VTLYP (DF2000,DF0200,DF1100,DF1010,
     &                 DF0101,DF1001,DF0110,DF10001,
     &                 DF01001,RHOA,RHOA,RHGA,RHGA,RHGA)
           FZ0LYP = DP5*(DF0010 + DF00001)*RHOGRD
           FRRLYP = DP5*(DF2000 + DF1100)
           FRZLYP = DP5*(DF1010 + DF1001+DF10001)*RHOGRD
           FZZLYP = FZ0LYP/RHOGRD 
           ZKSR(3,1) = FRRLYP
c
           IF (IWINT.EQ.1) THEN
              CALL WVWN(VLDA,RHO,RHO13,ELDA,.TRUE.,.TRUE.)
              CALL V1VWN(VVWN0,FRRVWN,RHO,RHO13)
              ZKSR(1,1) = ZKSR(1,1) - ELDA*RHO
              ZKSR(2,1) = ZKSR(2,1) - VLDA
              ZKSR(3,1) = ZKSR(3,1) - FRRVWN
              CALL VSRWINT(ZKSR,ZMU,RHO,RHO13,CHIVAL,2,ERFEXP)
           ELSE IF (IWINT.EQ.2.OR.IWINT.EQ.3 .OR. IWINT.GT.10) THEN
              IF (IWINT.GT.10) CALL DZERO(ZMU,3)
              CALL VSRWINT2(ZKSR,ZMU,RHO,RHO13,CHIVAL,2,ERFEXP)
           ELSE
              CALL QUIT ('Wrong dft functional interpolation '//
     &       'scheme in correlation part of SRDFTLTR (IWINT).')
           END IF
c          --- Remember mixed terms from lyp
           IF (IWINT.EQ.3) THEN
              FRZTMP = FRZLYP*GAUX
              FZZTMP = FZZLYP*GAUX
           ENDIF
           FRZLYP = FRZLYP*DERFCX + FRZTMP
           FZZLYP = FZZLYP*DERFCX + FZZTMP
c
         ELSE IF (SRHYBR) THEN
c        ... Special case for SRDFT in DFT embedding.
            RHOA = DP5*RHO
            RHGA = (DP5*RHOGRD)**2
c
#ifdef LDA
            CALL V1VWN(VVWN0,FRRLYP,RHO,RHO13)
            FZ0LYP = 0.0D0
            FRZLYP = 0.0D0
            FZZLYP = 0.0D0
#else
            CALL GLYPCO(DF1000,DF0100,DF0010,DF0001,
     &                  DF00001,RHO,RHO13,RHOGRD,.TRUE.)
            CALL VTLYP (DF2000,DF0200,DF1100,DF1010,
     &                  DF0101,DF1001,DF0110,DF10001,
     &                  DF01001,RHOA,RHOA,RHGA,RHGA,RHGA)
            FZ0LYP = DP5*(DF0010 + DF00001)*RHOGRD
            FRRLYP = DP5*(DF2000 + DF1100)
            FRZLYP = DP5*(DF1010 + DF1001+DF10001)*RHOGRD
            FZZLYP = FZ0LYP/RHOGRD
c
#endif
            IF (LSRLDA) THEN
               ZKSR(3,1) =  ZKSR(3,1) - FRRLYP
               FZ0LYP = -FZ0LYP
               FRRLYP = -FRRLYP
               FRZLYP = -FRZLYP
               FZZLYP = -FZZLYP
            ELSE
               ZKSR(3,1) =  FRRLYP
            END IF
C
         ELSE
            FRRVWN = 0D0
            FZ0LYP = 0D0
            FRRLYP = 0D0
            FRZLYP = 0D0
            FZZLYP = 0D0
         END IF
C
C        ************************
C        ** Total contribution **
C        ************************
C
         ZNV     = D1/RHOGRD
#ifdef LDA
         ZNV     = 0.d0
#endif
         FRR = WGHT*(ZKSR(3,1) + XKSR(3,1))
         FZ0 = ZNV*WGHT*(FZ0BCK + FZ0LYP)
         FRZ = WGHT*(FRZBCK + FRZLYP)
         FZZ = WGHT*(FZZBCK + FZZLYP)
C
         RX = ZNV*RHG(1)
         RY = ZNV*RHG(2)
         RZ = ZNV*RHG(3)
C
C        ***************************
C        ** Linear transformation **
C        ***************************
C
         IF (.NOT.SRHYBR .OR. DOHFSRDFT) THEN
            IF (DOGGA) THEN
               BR = B3(1)*RX + B3(2)*RY + B3(3)*RZ
               FAC0 = FRR*B0 + FRZ*BR
               FACR = FRZ*B0 + FZZ*BR
            ELSE
               FAC0 = FRR*B0
            END IF
            IF (NSYM.EQ.1) THEN
               DO I = 1, NBAST
                  G0 = GAO(I)
                  IF (DOGGA) THEN
                     GX = GAO1(I,1)
                     GY = GAO1(I,2)
                     GZ = GAO1(I,3)
                  END IF
                  DO J = 1, I 
                     EXCMAT(J,I) = EXCMAT(J,I)+FAC0*G0*GAO(J)
                     IF (DOGGA) THEN
                        AX = GX*GAO(J) + G0*GAO1(J,1)
                        AY = GY*GAO(J) + G0*GAO1(J,2)
                        AZ = GZ*GAO(J) + G0*GAO1(J,3)
                        AR = AX*RX + AY*RY + AZ*RZ
                        AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) - AR*BR
                        EXCMAT(J,I) = EXCMAT(J,I)+FACR*AR+FZ0*AB
                     END IF
                  END DO
               END DO
            ELSE
               DO ISYM = 1, NSYM
                  ISTR = IBAS(ISYM) + 1
                  IEND = IBAS(ISYM) + NBAS(ISYM)
                  JSYM = MULD2H(ISYM,KSYMOP)
                  IF (ISYM.GE.JSYM) THEN
                     JSTR = IBAS(JSYM) + 1
                     JEND = IBAS(JSYM) + NBAS(JSYM)
                     DO I = ISTR, IEND
                        G0 = GAO(I)
                        IF (DOGGA) THEN
                           GX = GAO1(I,1)
                           GY = GAO1(I,2)
                           GZ = GAO1(I,3)
                        END IF
                        DO J = JSTR, MIN(I,JEND) 
                           EXCMAT(J,I) = EXCMAT(J,I) + FAC0*G0*GAO(J)
                           IF (DOGGA) THEN
                              AX = GX*GAO(J) + G0*GAO1(J,1)
                              AY = GY*GAO(J) + G0*GAO1(J,2)
                              AZ = GZ*GAO(J) + G0*GAO1(J,3)
                              AR = AX*RX + AY*RY + AZ*RZ
                              AB = AX*B3(1) + AY*B3(2) + AZ*B3(3) -
     &                             AR*BR
                              EXCMAT(J,I) = EXCMAT(J,I)
     &                                     + FACR*AR + FZ0*AB
                           END IF
                        END DO
                     END DO
                  END IF
               END DO
            END IF
         END IF
      END IF
c
      RETURN
      END
C  /* Deck srfmat */
      SUBROUTINE SRFMAT(FMAT,CMO,DVREF,EJCSR,EJVSR,EDSR,EDFT,
     >                  WRK,LWRK,IPRINT)
C
C     20-Mar-2003 Jesper K. Pedersen
C
C     Driver for Fock matrix construction of short-range terms
C     needed for CI-DFT (J_ee,E_xc)
C
#include <implicit.h>
      DIMENSION FMAT(*), CMO(*), DVREF(*), WRK(*)
      DIMENSION ISYMDM(2), IFCTYP(2)
      LOGICAL   LSRHYBR
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0, DM1 = -1.0D0)
      SAVE ISYMDM, IFCTYP
      DATA ISYMDM/1,1/, IFCTYP/999,999/
C
#include <dummy.h>
#include <maxorb.h>
#include <mxcent.h>
#include <priunit.h>
#include <gnrinf.h>
#include <infvar.h>
#include <inftap.h>
#include <inforb.h>
#include <infinp.h>
#include <infpri.h>
#include <dftcom.h>
C
      CALL QENTER('SRFMAT')
C
      LSRHYBR = SRHYBR .AND. .NOT.DOHFSRDFT
      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWRK
      CALL MEMGET('REAL',KDCAO,2*N2BASX,WRK,KFREE,LFREE)
      KDVAO = KDCAO + N2BASX
      IF (LSRHYBR) THEN
         NF = 3
      ELSE
         NF = 2
      END IF
      CALL MEMGET('REAL',KFCAO,NF*N2BASX,WRK,KFREE,LFREE)
      KFVAO = KFCAO + N2BASX
      KFHAO = KFVAO + N2BASX
C
      CALL FCKDEN(.TRUE.,.TRUE.,
     >            WRK(KDCAO),WRK(KDVAO),CMO,DVREF,WRK(KFREE),LFREE)
C
      IF (IPRINT.GT.5) THEN
         WRITE(LUPRI,'(//5X,A//)') 
     >   '<<<<  Valence AO density matrix >>>>'
         CALL OUTPUT(WRK(KDVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
C
C
C ... Get J^c,sr adn J^v,sr
C
      IF(HFXFAC .NE. D0) THEN
C        Coulomb + Exchange
C     ... Warning : assuming DIM(IFCTYP) = 2
         IFCTYP(1) = 13
         IFCTYP(2) = 13
      ELSE
C        only Coulomb (Hartree) term
         IFCTYP(1) = 11
         IFCTYP(2) = 11
      ENDIF
      CALL DZERO(WRK(KFCAO),NF*N2BASX)
      CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KFCAO),WRK(KDCAO),2,
     >             ISYMDM,IFCTYP,WRK(KFREE),LFREE)
      EJCSR =  DP5*DDOT(N2BASX,WRK(KDCAO),1,WRK(KFCAO),1)
      EJVSR = -DP5*DDOT(N2BASX,WRK(KDVAO),1,WRK(KFVAO),1)
      IF (IPRINT.GT.5) THEN
         WRITE(LUPRI,'(//5X,A//)') 
     >   '<<<<  Short-range core Hartree matrix >>>>'
         CALL OUTPUT(WRK(KFCAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
         WRITE(LUPRI,'(//5X,A//)') 
     >   '<<<<  Short-range valence Hartree matrix >>>>'
         CALL OUTPUT(WRK(KFVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
C
c ... Add up in KFCAO
C
      CALL DAXPY(N2BASX,D1,WRK(KFVAO),1,WRK(KFCAO),1)
c
c ... Get V^sr in WRK(KFVAO) and E_xc^sr
c     using total density matrix (stored in KDCAO)
      CALL DAXPY(N2BASX,D1,WRK(KDVAO),1,WRK(KDCAO),1)
c
      CALL DZERO(WRK(KFVAO),N2BASX)
C
C     SRDFT(EXCMAT,DMAT,EDFTY,DOERG,DOATR,WORK,LWORK,IPRINT)
      CALL SRDFT(WRK(KFVAO),WRK(KDCAO),EDFT,.TRUE.,.FALSE.,
     >           WRK(KFREE),LFREE,IPRFCK)
      IF (IPRINT.GT.5) THEN
         WRITE(LUPRI,'(//5X,A//)') 
     >   '<<<<  Short-range Exchange-Correlation Fock-matrix >>>>'
         CALL OUTPUT(WRK(KFVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
      ENDIF
      EDSR = -DDOT(N2BASX,WRK(KDVAO),1,WRK(KFVAO),1)
      IF (LSRHYBR) THEN
         write (lupri,*) 'WARNING, srdft-hybr exc matrix not used'
C TODO
      END IF
c
c ... Add the Kohn-Sham Fock matrix
c
      CALL DAXPY(N2BASX,D1,WRK(KFVAO),1,WRK(KFCAO),1)
c
      IF (DOSRBCK) THEN
c     ... Get K^lr using total density and subtract it!
          IFCTYP(1) = 12
C         ... only exchange
          CALL DZERO(WRK(KFVAO),N2BASX)
          CALL SIRFCK2(LUINTA,'AOTWOINT',WRK(KFVAO),WRK(KDCAO),1,
     >                 ISYMDM,IFCTYP,WRK(KFREE),LFREE)
          EKVLR =     DDOT(N2BASX,WRK(KDVAO),1,WRK(KFVAO),1)
          CALL DAXPY(N2BASX,DM1,WRK(KFVAO),1,WRK(KFCAO),1)
          EKTLR = DP5*DDOT(N2BASX,WRK(KDCAO),1,WRK(KFVAO),1)
          EDSR = EDSR + EKVLR
          EDFT = EDFT - EKTLR
      ENDIF
C
C     ===============================================
C      Transform short-range Fock-matrix to MO basis
C     ===============================================
C
      CALL DCOPY(N2BASX,WRK(KFCAO),1,WRK(KDVAO),1)
      CALL DGETSP(NBAST,WRK(KDVAO),WRK(KFCAO))
      IF (NSYM.GT.1) CALL PKSYM1(WRK(KFCAO),WRK(KFCAO),NBAS,NSYM,2)
      CALL UTHUB(WRK(KFCAO),FMAT,CMO,WRK(KDVAO),NSYM,NBAS,NORB)
C
      IF (P6FLAG(14)) THEN
        WRITE(LUPRI,1200)
        CALL OUTPKB(FMAT,NORB,NSYM,1,LUPRI)
      END IF
 1200 FORMAT(/' Short-range Fock matrix (MO-basis)')
C
      CALL MEMREL('SRFMAT',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('SRFMAT')
      RETURN
      END
C  /* Deck srdftso */
      SUBROUTINE SRDFTSO(DCVAL,DV,FTC_SR,SOVEC)
C
C  Copyright Dec-2003 Hans Joergen Aa. Jensen
C  for MCSRDFT (based on SOLGO code).
C
C  Purpose: Add FTC_SR orbital sigma vector contribution
C           to SOVEC.
C
#include <implicit.h>
      DIMENSION DV(*), FTC_SR(NORBT,*), SOVEC(*)
C
C  Used from common blocks:
C    INFORB: NISHT,NASHT
C    INFVAR: NWOPH,JWOP(2,*)
C    INFIND: IROW(*),IOBTYP(*),ISX(*)
C
#include <maxash.h>
#include <maxorb.h>
#include <inforb.h>
#include <infvar.h>
#include <infind.h>
#include <dftcom.h>
C
      PARAMETER ( D0 = 0.0D0, D2 = 2.0D0 )
C
      CALL QENTER('SRDFTSO ')
      DCFAC = D2 * DCVAL
C     ... DCVAL is the value of DC = <0| Eii |B>,
C         = 2 <0|0> = 2 normally,
C         = 2 <0|B> = 0 for transition density matrix
      DO 300 IG = 1,NWOPH
         K     = JWOP(1,IG)
         L     = JWOP(2,IG)
         ITYPK = IOBTYP(K)
         ITYPL = IOBTYP(L)
         IF (ITYPK .EQ. JTINAC) THEN
C           first index inactive
            SOVEC(IG) = SOVEC(IG) + DCFAC * FTC_SR(K,L)
         ELSE
C           first index active
            NK   = ICH(K)
            TEMP = D0
            DO 100 NX = 1,NASHT
               IF (NX .LE. NK) THEN
                  DVKX  = DV(IROW(NK)+NX)
               ELSE
                  DVKX  = DV(IROW(NX)+NK)
               END IF
               IX = ISX(NISHT + NX)
               TEMP = TEMP + DVKX*FTC_SR(IX,L)
  100       CONTINUE
            SOVEC(IG) = SOVEC(IG) + D2 * TEMP
         END IF
         IF (ITYPL .EQ. JTACT) THEN
C           second index active
            NL   = ICH(L)
            TEMP = D0
            DO 200 NX = 1,NASHT
               IF (NX .LE. NL) THEN
                  DVLX  = DV(IROW(NL)+NX)
               ELSE
                  DVLX  = DV(IROW(NX)+NL)
               END IF
               IX = ISX(NISHT + NX)
               TEMP = TEMP + DVLX*FTC_SR(IX,K)
  200       CONTINUE
            SOVEC(IG) = SOVEC(IG) - D2 * TEMP
         END IF
  300 CONTINUE
C
      CALL QEXIT('SRDFTSO ')
      RETURN
C     end of srdftso.
      END
