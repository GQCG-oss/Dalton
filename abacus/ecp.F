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
      SUBROUTINE ECPINP(Q,QEFF,IATOM,IPRINT,ECPNAME)
C
C Purpose:
C     Find the atom corresponding to charge Q in the ecp pseudo
C     basis file and read in information. 
C     Set QEFF to be the effective charge of the atom. 
C     Store the position IATOM in the list of all atoms to later
C     identify the position of ECP atoms.
C
#include <implicit.h>
#include <priunit.h>
#include <gnrinf.h>
#include <maxorb.h>
#include <infpar.h>
#include <ecpinf.h>
C
      CHARACTER*80 STRING
      CHARACTER LAB
C...  added by Bin Gao, for effective core potentials, Dec. 30, 2010
      character*(*) ECPNAME
C...
C
C     Slaves can return here.
C
      IF (SLAVE) RETURN
C
C     Initialize ECP specifics.
C
      IF (.NOT. ECP) THEN
         CALL CORTAB
         ECP = .TRUE.
      END IF
C
C...  modified by Bin Gao, Dec. 30, 2010
C...      STRING = 'ecp-sdd-pseudo '
      STRING = "ecp_data/"//trim(ECPNAME)
      LEN_STRING = LNBLNK(STRING)
      LU = -1
      CALL FIND_BASFIL(STRING(1:LEN_STRING),LU,IPRINT)
C
C     Search for label identfying atom
C
      INTQINP = NINT(Q)
      REWIND(LU)
 20   CONTINUE
      READ(LU,'(A80)', ERR=900,END=910)
     &     STRING
      READ (STRING, '(A1)') LAB
      IF (LAB.EQ.'A' .OR. LAB.EQ.'a') THEN
         READ (STRING, '(A1,I4)') LAB,INTQ
         IF (INTQ.EQ.INTQINP) THEN
            GOTO 100
         ELSE
            GOTO 20
         END IF
      ELSE
         GOTO 20
      END IF
C     
C     Atomic label found
C
 100  CONTINUE
C
      NTYECP = NTYECP + 1
      IF (NTYECP.GT.MXECP) THEN
         WRITE(LUPRI,*) ' NUMBER OF ECP ATOM TYPES EXCEEDS MXECP=',MXECP
         CALL QUIT('Too many ECP atom types')
      END IF
      IF (NTYECP.EQ.1) THEN
         IOFF = 0
      ELSE 
         IOFF = NKCRU(LCR(NTYECP-1)+1,NTYECP-1)
      END IF
C
C     Read in ECP parameters
C
      READ(LU,*,ERR=900,END=900)
      READ(LU,'(2I4)',ERR=900,END=900) LMAX,NCORE
C
      QEFF = Q - NCORE
      LCR(NTYECP) = LMAX
C
C     Local term of the ECP
C     
      READ(LU,*,ERR=900,END=900) NEXP
      NKCRL(1,NTYECP) =    1 + IOFF
      NKCRU(1,NTYECP) = NEXP + IOFF
      DO I = NKCRL(1,NTYECP),NKCRU(1,NTYECP)
         READ(LU,'(I2,2(F20.12))',ERR=900,END=900) 
     &        NCR(I),ZCR(I),CCR(I)
      END DO
C
C     Nonlocal term of the ECP
C
      DO I = 1,LMAX
         READ(LU,*,ERR=900,END=900) NEXP
         NKCRL(I+1,NTYECP) = NKCRU(I,NTYECP) + 1
         NKCRU(I+1,NTYECP) = NKCRU(I,NTYECP) + NEXP
         DO J = NKCRL(I+1,NTYECP),NKCRU(I+1,NTYECP)
            READ(LU,'(I2,2(F20.12))',ERR=900,END=900) 
     &           NCR(J),ZCR(J),CCR(J)
         END DO
      END DO
C
      CALL GPCLOSE(LU,'KEEP')
C
C     Print
C
      IF (IPRINT.GE.0) THEN
         WRITE(LUPRI,'(5X,A,I2,A,/5X,A,2(/5X,A,F5.1),2(/5X,A,I5))')
     &        'Pseudo potential parameters for ECP type ',NTYECP,'.',
     &        '-----------------------------------------------',
     &        'Atomic charge  :', Q,
     &        'Effectiv charge:', QEFF,
     &        'Core electrons :', NCORE,
     &        'Max angular mom:', LMAX-1 
      END IF
C
      IF (IPRINT.GE.2) THEN
         WRITE(LUPRI,'(2(/A12,A4,2A20))')
     &        ' Type       ',' L ','Exponent       ','Coefficient    ',
     &        ' --------   ','---','---------------','---------------'
         DO I = 1,LMAX+1
            DO J = NKCRL(I,NTYECP),NKCRU(I,NTYECP)
               IF (I.EQ.1) THEN
                  WRITE(LUPRI,'(A12,A3,2F20.12)')
     &                 'Local      ',' ',ZCR(J),CCR(J)
               ELSE
                  WRITE(LUPRI,'(A12,I3,2F20.12)')
     &                 'Nonlocal   ',I-2,ZCR(J),CCR(J)
               END IF
            END DO
         END DO
      END IF
C
      RETURN
C
C     Error messages 
C
 900  CONTINUE
      WRITE (LUPRI,'(/A)')
     &   'ECPINP: Error when reading from basis file ECP_pseudo'
      CALL QUIT('Read error in ECPINP reading basis file ECP_pseudo')
 910  CONTINUE
      WRITE (LUPRI,'(/I3,A)') INTQINP,
     &     ' is an unsupported element in basis file ECP_pseudo'
      CALL QUIT('Unsupported element in ECPINP')
C
      END
      SUBROUTINE ECPINT(X2,Y2,Z2,X1,Y1,Z1,EXP2,EXP1,L2,L1,C2,C1,ADER,
     &                  WORK,LWORK,IPRINT)
C
C Purpose: Calculatate ECP integrals
C
#include <implicit.h>
#include <priunit.h>
#include <ecpinf.h>
#include <mxcent.h>
#include <nuclei.h>
#include <pi.h>
#include <maxorb.h>
#include <maxaqn.h>
#include <symmet.h>
C
      COMMON /TRQ88/  KCRS,LCRU
      COMMON /TRQ77/  XI,YI,ZI,XJ,YJ,ZJ,XC,YC,ZC
C
      PARAMETER ( D4 = 4.0D0, D1 = 1D0 )
      DIMENSION WORK(*),ADER(L1*(L1+3)/2+1,L2*(L2+3)/2+1)
C
C     Used to scatter GOUT into ADER, the reason for the difference
C     in ordering is that PSEUDO routines order according to ZYX whereas
C     Dalton uses XYZ. Maximum g-functions at the moment.
C
      PARAMETER (MAX_L = 4)
      DIMENSION IORDER(35),IOFF(0:MAX_L)
      DATA IORDER/1,
     &            1,2,3,
     &            1,4,6,2,3,5,
     &            1,7,10,2,3,4,8,6,9,5,
     &            1,11,15,2,3,7,12,10,14,4,6,13,5,8,9/
      DATA IOFF/0,1,4,10,20/  
C
#include <ibtfun.h>
C
      IF (L1 .GT. MAX_L .OR. L2 .GT. MAX_L) THEN
         write(lupri,'(//A,3I8)')
     &   'ERROR ECPINT: L1 or L2 .gt. MAX_L',L1,L2,MAX_L
         CALL QUIT('ECPINT: L1 or L2 .gt. MAX_L')
      END IF
C
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KGOUT,((MXANG*MXANG+3*MXANG+2)/2)**2,
     &            WORK,KFREE,LFREE)
C
C     Transfer coordinates of primitives to common block.
C
      XI = X1
      YI = Y1
      ZI = Z1
      XJ = X2
      YJ = Y2
      ZJ = Z2
C
C     Loop over number of all ECP centers
C
      CALL DZERO(WORK(KGOUT),((MXANG*MXANG+3*MXANG+2)/2)**2)
C
      DO ITYP = 1,NTYECP
         DO INONT = 1,NECP(ITYP)
            IATOM = INDECP(ITYP)+INONT-1
            KCRS  = ITYP
            LCRU  = LCR(ITYP)
            XCI   = CORD(1,IATOM) 
            YCI   = CORD(2,IATOM) 
            ZCI   = CORD(3,IATOM) 
            DO ISYMOP = 0, MAXOPR
               IF (IBTAND(ISYMOP,ISTBNU(IATOM)) .EQ. 0) THEN
                  XC = PT(IBTAND(ISYMAX(1,1),ISYMOP))*XCI
                  YC = PT(IBTAND(ISYMAX(2,1),ISYMOP))*YCI
                  ZC = PT(IBTAND(ISYMAX(3,1),ISYMOP))*ZCI
                  CALL PSEUD1(WORK(KGOUT),L1+1,L2+1,EXP1,EXP2)
                  CALL PSEUD2(WORK(KGOUT),L1+1,L2+1,EXP1,EXP2)
               END IF
            END DO
         END DO
      END DO
C
      CALL DSCAL((L1*(L1+3)/2+1)*(L2*(L2+3)/2+1),D4*PI,WORK(KGOUT),1)
C
      IF (IPRINT.GE.10) THEN
         WRITE(LUPRI,'(2(/A))') 
     &        ' ADER in ECPINT before ecp integrals added',
     &        '===========================================' 
         DO I = 1,L1*(L1+3)/2+1
            WRITE(LUPRI,'(I3,A4,15F12.8)') I,' ',
     &      (ADER(I,J), J=1,L2*(L2+3)/2+1)
         END DO
         WRITE(LUPRI,'(2(/A,F12.8))')
     &        ' Contraction coefficent 1 = ', C1,
     &        ' Contraction coefficent 2 = ', C2
      END IF
C
      IJ = 0
      DO I = 1,L1*(L1+3)/2+1
         DO J = 1,L2*(L2+3)/2+1
            ADER(IORDER(IOFF(L1)+I),IORDER(IOFF(L2)+J)) =
     &           ADER(IORDER(IOFF(L1)+I),IORDER(IOFF(L2)+J)) +
     &           C1*C2*WORK(KGOUT+IJ)
            IJ = IJ + 1
         END DO
      END DO
C
C     Print
C
      IF (IPRINT.GE.10) THEN
         WRITE(LUPRI,'(2(/A))') 
     &        ' ADER in ECPINT after ecp integrals added',
     &        '==========================================' 
         DO I = 1,L1*(L1+3)/2+1
            WRITE(LUPRI,'(I3,A4,15F12.8)') I,' ',
     &      (ADER(I,J), J=1,L2*(L2+3)/2+1)
         END DO
      END IF
C
      RETURN
      END
      subroutine pseud1(gout,lit,ljt,ai,aj)
C ***************************************************************************
C These subroutines were written by Dr. Heinz-Juergen Flad, Stuttgart 
C ***************************************************************************
C
c pseudopotential parameters
c
c Vpp = sum_i ( r^(ncr(i)-2) ccr(i) exp(-zcr(i)*r^2) )      (local part)
c     + sum_l (sum_i ( r^(ncr(i)-2) ccr(i) exp(-zcr(i)*r^2) ) |l><l| ) 
c                                                        (nonlocal part)
c ncr(i) exponents n for r**(n-2)
c zcr(i) exponents
c ccr(i) coefficients
c lcr(k) maximal l value + 1 of the pseudopotential
c nkcrl(l,k) lower index of the pseudopotential i for projection operator l
c nkcru(l,k) upper index of the pseudopotential i for projection operator l
c
c kcrs specifies the pseudopotential
c lcru lcr(kcrs)
c nkcrl(1,kcrs) to nkcrl(1,kcrs), local part of the pseudopotential
c nkcrl(l+1,kcrs) to nkcrl(l+1,kcrs), nonlocal parts of the pseudopotential
c
c xi.., xj.., xc... coordinates of the basis functions i,j and of the 
c pseudopotential c
c ai, aj exponents of the basis functions in pseud1, pseud2
c lit, ljt l value + 1 of the basis functions in pseud1, pseud2
c
c Integrals between two shells are contained in gout.
c
c*******************************************************
c     ----- computes type 1 core potential integrals -----
      implicit double precision (a-h,o-z)
#include <ecpinf.h>      
c*
      dimension llt(5,2),crda(5,3),crdb(5,3),gout(225)
      data llt /1,2,5,11,21,1,4,10,20,35/
CMM   common /trq99/ lcr(50),ncr(1000),nkcrl(6,50),nkcru(6,50),
CMM  a               zcr(1000),ccr(1000)
      common /trq88/ kcrs,lcru
      common /trq77/ xi,yi,zi,xj,yj,zj,xc,yc,zc
c     first center, second center, nucleus
c*    common /dfac/ dfac(23)
      common /dfac/ dfac(29)
c*
c     logical esf, esfc, igueq1, jgueq1
c     common /parma/ cutoff, tol, mblu, mccu, mconu, mcu, mcxu,
c    1   mrcru, mstu, msu, nblu, ng, nlist, ns, nst, ntape
c     common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
c    1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ijx(225),ijy(225),ijz(225),ic,
c    2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
c    3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
c     common /ccali / crda(4,3),crdb(4,3),xka,yka,zka,ca,xkb,ykb,zkb,
c    1   cb,tai,taj,aa,taa,aarr1,aarr2,xk,yk,zk,fctr2,ltot1,kcrs,lcru
c*    common/qstore/q(7,7),alpha,rk,t
      common/qstore/q(9,9),alpha,rk,t
c     dimension ccr(2), gout(2), ncr(2), nkcrl(6,2), nkcru(6,2), zcr(2)
c*    dimension ang(7,7), qsum(7,7), xab(7), yab(7), zab(7)
      dimension ang(9,9), qsum(9,9), xab(9), yab(9), zab(9)
      dimension nap(35), lap(35), map(35)
      data nap /0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,
     1           4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/,
     2     lap /0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,
     3           0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/,
     4     map /0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,
     5           0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
      data a0 /0.0d0/
c*
      tol=25.32838
      fctr2=1.d0
      itl=llt(lit,1)
      itu=llt(lit,2)
      jtl=llt(ljt,1)
      jtu=llt(ljt,2)
      aa=ai+aj
      ltot1=lit+ljt-1
      do 120 i=1,3
      crda(1,i)=1.d0
  120 crdb(1,i)=1.d0
      xka=xc-xi
      yka=yc-yi
      zka=zc-zi
      ca=dsqrt(xka*xka+yka*yka+zka*zka)
      if(lit.eq.1) go to 220
      crda(2,1)=xka
      crda(2,2)=yka
      crda(2,3)=zka
      if(lit.eq.2) go to 220
      do 210 i=1,3
      do 210 j=3,lit
  210 crda(j,i)=crda(2,i)*crda(j-1,i)
  220 xkb=xc-xj
      ykb=yc-yj
      zkb=zc-zj
      cb=dsqrt(xkb*xkb+ykb*ykb+zkb*zkb)
      if(ljt.eq.1) go to 240
      crdb(2,1)=xkb
      crdb(2,2)=ykb
      crdb(2,3)=zkb
      if(ljt.eq.2) go to 240
      do 230 i=1,3
      do 230 j=3,ljt
  230 crdb(j,i)=crdb(2,i)*crdb(j-1,i)
  240 continue
      xij=0.5d0*(xi+xj)
      yij=0.5d0*(yi+yj)
      zij=0.5d0*(zi+zj)
      xijm=0.5d0*(xi-xj)
      yijm=0.5d0*(yi-yj)
      zijm=0.5d0*(zi-zj)
      rr=(xi-xj)**2+(yi-yj)**2+(zi-zj)**2
      aaa=(ai-aj)/aa
      xk=xij+aaa*xijm-xc
      yk=yij+aaa*yijm-yc
      zk=zij+aaa*zijm-zc
      aarr1=(ai*aj/aa)*rr
      taa=aa+aa
c*
      rp2=xk*xk+yk*yk+zk*zk
      if(rp2.ne.a0) go to 10
      rp=a0
      arp2=a0
      alpt=a0
      rk=a0
      lamu=1
      go to 20
   10 rp=dsqrt(rp2)
      xk=xk/rp
      yk=yk/rp
      zk=zk/rp
      arp2=aa*rp2
      alpt=aa*arp2
      rk=taa*rp
      lamu=ltot1
c     ----- compute radial integrals and sum over potential terms -----
c* 20 call dzero(qsum,(7*lamu))
   20 call dzero(qsum,(9*lamu))
      kcrl=nkcrl(1,kcrs)
      kcru=nkcru(1,kcrs)
      do 40 kcr=kcrl,kcru
      alpha=aa+zcr(kcr)
c     ----- exponential factor from q functions included in dum -----
      dum=aarr1+zcr(kcr)*arp2/alpha
      if(dum.gt.tol) go to 40
      prd=fctr2*ccr(kcr)*dexp(-dum)
      t=alpt/alpha
      call recur1(ncr(kcr),ltot1)
      do 30 lam=1,lamu
      nhi=ltot1-mod(ltot1-lam,2)
      do 30 n=lam,nhi,2
   30 qsum(n,lam)=qsum(n,lam)+prd*q(n,lam)
   40 continue
      ijt=0
      do 90 it=itl,itu
      na1=nap(it)+1
      la1=lap(it)+1
      ma1=map(it)+1
      do 90 jt=jtl,jtu
      ijt=ijt+1
      s=a0
      nb1=nap(jt)+1
      lb1=lap(jt)+1
      mb1=map(jt)+1
c     ----- compute angular integrals -----
      call facab(na1,nb1,crda(1,1),crdb(1,1),xab)
      call facab(la1,lb1,crda(1,2),crdb(1,2),yab)
      call facab(ma1,mb1,crda(1,3),crdb(1,3),zab)
      call ang1(na1+nb1-1,la1+lb1-1,ma1+mb1-1,lamu,xab,yab,zab,
     1   xk,yk,zk,ang)
c     ----- combine angular and radial integrals -----
      do 80 lam=1,lamu
      nhi=ltot1-mod(ltot1-lam,2)
      do 80 n=lam,nhi,2
   80 s=s+ang(n,lam)*qsum(n,lam)
   90 gout(ijt)=gout(ijt)+s
      return
      end
c
      subroutine pseud2(gout,lit,ljt,ai,aj)
c     ----- COMPUTES TYPE 2 CORE POTENTIAL INTEGRALS -----
      implicit double precision (a-h,o-z)
#include <ecpinf.h>
      dimension llt(5,2),crda(5,3),crdb(5,3),gout(225)
      data llt /1,2,5,11,21,1,4,10,20,35/
CMM   common /trq99/ lcr(50),ncr(1000),nkcrl(6,50),nkcru(6,50),
CMM  a               zcr(1000),ccr(1000)
      common /trq88/ kcrs,lcru
      common /trq77/ xi,yi,zi,xj,yj,zj,xc,yc,zc
c*    common /dfac/ dfac(23)
      common /dfac/ dfac(29)
c     include "common/mxang"
c     parameter (mxx=((mxang+1)*(mxang+2)*(mxang+3))/6)
c     parameter (mxy=(mxang+1)*(mxang+1), mxang1=mxang+1)
c     parameter (mxz=((mxang+1)*(mxang+2)/2)**2)
c     logical es, esf, esfc, igueq1, jgueq1
c     common /parma/ cutoff, tol, mblu, mccu, mconu, mcu, mcxu,
c    1   mrcru, mstu, msu, nblu, ng, nlist, ns, nst, ntape
c*    COMMON /PARMI/ MBLU, MCCU, MCONU, MCU, MCXU,
c*   1   MRCRU, MSTU, MSU, NBLU, NG, NS, NST
c     common /one/dxij,dyij,dzij,fnfct,rr,xij,xijm,yij,yijm,zij,zijm,
c    1  ibl1,ibl2,icxi1,icxi2,ij,ijsf,ijx(mxz),ijy(mxz),ijz(mxz),ic,
c    2  icons,igu,ircru,is,isf,itl,itu,jc,jcons,jgu,jrcru,js,jsf,jtl,
c    3  jtu,lit,ljt,nblt1,nc2,nc1,nop,ntij1,ntij2,esf,esfc,igueq1,jgueq1
c*    COMMON /CCALI/ CRDA(4,3),CRDB(4,3),XKA,YKA,ZKA,CA,XKB,YKB,ZKB,
c*   1   CB,TAI,TAJ,AA,TAA,AARR1,AARR2,XK,YK,ZK,FCTR2,LTOT1,KCRS,LCRU
c     common /ccali/ crda(5,3),crdb(5,3),xka,yka,zka,ca,xkb,ykb,zkb,
c    1   cb,tai,taj,aa,taa,aarr1,aarr2,xk,yk,zk,fctr2,ltot1,kcrs,lcru
c     COMMON /CALLIN/ CRDA(4,3),CRDB(4,3),XKA,YKA,ZKA,CA,XKB,YKB,ZKB,
c    1   CB,TAI,TAJ,AA,TAA,AARR1,AARR2,XK,YK,ZK,FCTR2,LTOT1,KCRS,LCRU
c     dimension ccr(2),gout(2),ncr(2),nkcrl(6,2),nkcru(6,2),zcr(2)
c*    DIMENSION ANGA(4,7,7), ANGB(4,7,7), QSUM(7,7,7)
      dimension anga(5,9,9), angb(5,9,9), qsum(9,9,9)
c*    DATA A0 /0.0D0/
      data a0,eps1,a1,a4,a50 /0.0d0,1.0d-15,1.0d0,4.0d0,50.0d0/
c*
      tol=25.32838
      fctr2=1.d0
      itl=llt(lit,1)
      itu=llt(lit,2)
      jtl=llt(ljt,1)
      jtu=llt(ljt,2)
      tai=ai+ai
      taj=aj+aj
      aa=ai+aj
      ltot1=lit+ljt-1
      do 120 i=1,3
      crda(1,i)=1.d0
  120 crdb(1,i)=1.d0
      xka=xc-xi
      yka=yc-yi
      zka=zc-zi
      ca=dsqrt(xka*xka+yka*yka+zka*zka)
      if(lit.eq.1) go to 220
      crda(2,1)=xka
      crda(2,2)=yka
      crda(2,3)=zka
      if(lit.eq.2) go to 220
      do 210 i=1,3
      do 210 j=3,lit
  210 crda(j,i)=crda(2,i)*crda(j-1,i)
  220 xkb=xc-xj
      ykb=yc-yj
      zkb=zc-zj
      cb=dsqrt(xkb*xkb+ykb*ykb+zkb*zkb)
      if(ljt.eq.1) go to 240
      crdb(2,1)=xkb
      crdb(2,2)=ykb
      crdb(2,3)=zkb
      if(ljt.eq.2) go to 240
      do 230 i=1,3
      do 230 j=3,ljt
  230 crdb(j,i)=crdb(2,i)*crdb(j-1,i)
  240 continue
      aarr2=(ai*aj/aa)*(ca-cb)**2
CJV 
C      IF (LIT.NE.2 .OR. LJT.NE.2) GOTO 999
C      WRITE(6,*) 'initialized variables in PSEUD2/DEMON' 
C      WRITE(6,'(A,2F12.5,4I3)') ' tol,fctr2,itl,itu,jtl,jtu ', 
C     &     tol,fctr2,itl,itu,jtl,jtu 
C      WRITE(6,'(A,3F12.5,I3)') ' ai,aj,aa,ltot1 ',tai/2.0,taj/2.0,aa,ltot1 
C      WRITE(6,'(A,4F12.5)') ' xka,yka,zka,ca ',xka,yka,zka,ca 
C      DO iiii = 1,5,1 
C         WRITE(6,'(A,3F12.5)') ' crda(j, ) ', 
C     &        (crda(iiii,jjjj), jjjj = 1,3,1) 
C      ENDDO 
C      WRITE(6,'(A,4F12.5)') ' xkb,ykb,zkb,cb ',xkb,ykb,zkb,cb  
C      DO iiii = 1,5,1  
C         WRITE(6,'(A,3F12.5)') ' crdb(j, ) ', 
C     &        (crdb(iiii,jjjj), jjjj = 1,3,1)  
C      ENDDO  
C      WRITE(6,'(A,F12.5)') ' aarr2 ',aarr2 
C 999  CONTINUE
c*
      if(ca.ne.a0) go to 12
      rka=a0
      lmau=1
      go to 16
   12 xka=-xka/ca
      yka=-yka/ca
      zka=-zka/ca
      rka=tai*ca
      lmau=lcru+(lit-1)
   16 if(cb.ne.a0) go to 20
      rkb=a0
      lmbu=1
      go to 24
   20 xkb=-xkb/cb
      ykb=-ykb/cb
      zkb=-zkb/cb
      rkb=taj*cb
      lmbu=lcru+(ljt-1)
   24 llo=1
      lhi=lcru
      inc=1
      if((ca.ne.a0).or.(cb.ne.a0)) go to 28
      lhi=min(lcru,lit,ljt)
      llo=mod((lit-1),2)+1
      if(llo.ne.mod((ljt-1),2)+1.or.llo.gt.lhi) return
      go to 36
   28 if(ca.ne.a0) go to 32
      lhi=min(lcru,lit)
      llo=mod((lit-1),2)+1
      if(llo.gt.lhi) return
      go to 36
   32 if(cb.ne.a0) go to 40
      lhi=min(lcru,ljt)
      llo=mod((ljt-1),2)+1
      if(llo.gt.lhi) return
   36 inc=2
   40 do 88 l=llo,lhi,inc
      mhi=l+l-1
      lmalo=max(l-(lit-1),1)
      lmahi=min(lmau,l+(lit-1))
      lmblo=max(l-(ljt-1),1)
      lmbhi=min(lmbu,l+(ljt-1))
c     ----- COMPUTE RADIAL INTEGRALS -----
c*    CALL DZERO(QSUM,(49*LMAHI))
      call dzero(qsum,(81*lmahi))
      kcrl=nkcrl(l+1,kcrs)
      kcru=nkcru(l+1,kcrs)
      do 64 kcr=kcrl,kcru
      npi=ncr(kcr)
      alpha=aa+zcr(kcr)
      rc=(rka+rkb)/(alpha+alpha)
      arc2=alpha*rc*rc
      dum=aarr2+zcr(kcr)*arc2/aa
      if(dum.gt.tol) go to 64
      if(arc2.lt.a50) go to 60
c     ----- use pts and wts method -----
      prd=fctr2*ccr(kcr)*dexp(-dum)
c     ----- estimate radial integrals and compare to threshold -----
      qlim=dabs(prd)/(dmax1(a1,(rc+rc)*rka)*dmax1(a1,(rc+rc)*rkb))*
     1   dsqrt(a4*(tai+tai)**lit*(taj+taj)**ljt*dsqrt(tai*taj)/alpha)
      if(rc.ge.ca) go to 44
      nlim=npi
      qlim=qlim*ca**(lit-1)
      go to 48
   44 nlim=npi+(lit-1)
   48 if(rc.ge.cb) go to 52
      qlim=qlim*cb**(ljt-1)
      go to 56
   52 nlim=nlim+(ljt-1)
   56 if(qlim*rc**nlim.lt.eps1) go to 64
      call ptwt(npi,l,lit,ljt,ltot1,lmalo,lmahi,lmblo,lmbhi,
     1   alpha,rc,rka,rkb,arc2,prd,qsum)
      go to 64
c     ----- use partially asymptotic method -----
   60 prd_init = fctr2*ccr(kcr)
      call qpasy(npi,l,lit,ljt,ltot1,lmalo,lmahi,lmblo,lmbhi,
     1   alpha,rka,rkb,prd_init,dum+arc2,qsum)
   64 continue
c*    CALL RAD2(CCR,KCRL,KCRU,L,LMAHI,LMALO,LMBHI,LMBLO,NCR,QSUM,
c*   1   RKA,RKB,ZCR)
      ijt=0
      do 88 it=itl,itu
      call ang2(it,crda,l,lmalo,lmahi,xka,yka,zka,anga)
      do 88 jt=jtl,jtu
      ijt=ijt+1
      s=a0
      call ang2(jt,crdb,l,lmblo,lmbhi,xkb,ykb,zkb,angb)
c     ----- COMBINE ANGULAR AND RADIAL INTEGRALS -----
      do 84 lama=lmalo,lmahi
      ldifa1=abs(l-lama)+1
      nlmau=lit-mod(lit-ldifa1,2)
      do 84 lamb=lmblo,lmbhi
      ldifb=abs(l-lamb)
      nlmbu=(ljt-1)-mod((ljt-1)-ldifb,2)
      nlo=ldifa1+ldifb
      nhi=ltot1-mod(lit-ldifa1,2)-mod((ljt-1)-ldifb,2)
      do 84 n=nlo,nhi,2
      nlmalo=max(ldifa1,n-nlmbu)
      nlmahi=min(nlmau,n-ldifb)
      angp=a0
      do 80 m=1,mhi
      do 80 nlma=nlmalo,nlmahi
   80 angp=angp+anga(nlma,m,lama)*angb((n+1)-nlma,m,lamb)
   84 s=s+angp*qsum(n,lamb,lama)
   88 gout(ijt)=gout(ijt)+s
CJV
C      if (lit.ne.2 .or. ljt.ne.2) goto 998
C      WRITE(6,997)(GOUT(I),I=1,9) 
C 997  FORMAT(' gout ', 3F16.9) 
C 998  continue
      return
      end
c
      subroutine ang1(nanb,lalb,mamb,lamu,xab,yab,zab,xk,yk,zk,ang)
c     ----- computes angular integrals for type 1 integrals -----
      implicit double precision (a-h,o-z)
c*    common /zlmtab/ zlm(130)
      common /zlmtab/ zlm(377)
      common /dfac/ dfac(29)
c*    common /dfac/ dfac(23)
c*    dimension ang(7,2), xab(2), yab(2), zab(2)
chj   dimension ang(9,2), xab(2), yab(2), zab(2)
      dimension ang(9,*), xab(*), yab(*), zab(*)
c*    dimension lf(7),lmf(49),lml(49),lmx(130),lmy(130),lmz(130)
      dimension lf(9),lmf(81),lml(81),lmx(377),lmy(377),lmz(377)
c*    data lf /0,1,4,9,16,25,36/
      data lf /0,1,4,9,16,25,36,49,64/
c*    data lmf/1,2,3,4,5,7,8,9,10,12,14,16,18,19,21,23,25,28,30,32,
c*   1   34,38,40,42,44,47,50,53,56,58,62,66,70,72,75,78,81,85,88,91,
c*   2   94,100,104,108,112,118,121,124,127/
c*    data lml/1,2,3,4,6,7,8,9,11,13,15,17,18,20,22,24,27,29,31,33,
c*   1   37,39,41,43,46,49,52,55,57,61,65,69,71,74,77,80,84,87,90,93,
c*   2   99,103,107,111,117,120,123,126,130/
c*    data lmx/0,0,0,1,0,0,0,1,1,2,0,0,0,0,0,1,1,1,2,0,2,0,3,1,0,0,0,
c*   1   0,0,1,1,1,1,2,0,2,0,2,0,3,1,3,1,4,2,0,0,0,0,0,0,0,1,1,1,1,1,
c*   2   2,0,2,0,2,0,2,0,3,1,3,1,3,1,4,2,0,4,2,0,5,3,1,0,0,0,0,0,0,0,
c*   3   1,1,1,1,1,1,2,0,2,0,2,0,2,0,2,0,3,1,3,1,3,1,3,1,4,2,0,4,2,0,
c*   4   4,2,0,5,3,1,5,3,1,6,4,2,0/
c*    data lmy/0,0,1,0,0,0,1,0,1,0,2,0,0,1,1,0,0,1,0,2,1,3,0,2,0,0,0,
c*   1   1,1,0,0,1,1,0,2,0,2,1,3,0,2,1,3,0,2,4,0,0,0,1,1,1,0,0,0,1,1,
c*   2   0,2,0,2,1,3,1,3,0,2,0,2,1,3,0,2,4,1,3,5,0,2,4,0,0,0,0,1,1,1,
c*   3   0,0,0,1,1,1,0,2,0,2,0,2,1,3,1,3,0,2,0,2,1,3,1,3,0,2,4,0,2,4,
c*   4   1,3,5,0,2,4,1,3,5,0,2,4,6/
c*    data lmz/0,1,0,0,2,0,1,1,0,0,0,3,1,2,0,2,0,1,1,1,0,0,0,0,4,2,0,
c*   1   3,1,3,1,2,0,2,2,0,0,1,1,1,1,0,0,0,0,0,5,3,1,4,2,0,4,2,0,3,1,
c*   2   3,3,1,1,2,2,0,0,2,2,0,0,1,1,1,1,1,0,0,0,0,0,0,6,4,2,0,5,3,1,
c*   3   5,3,1,4,2,0,4,4,2,2,0,0,3,3,1,1,3,3,1,1,2,2,0,0,2,2,2,0,0,0,
c*   4   1,1,1,1,1,1,0,0,0,0,0,0,0/
      data lmf/1,2,3,4,5,7,8,9,10,12,14,16,18,19,21,23,25,28,30,32,
     1   34,38,40,42,44,47,50,53,56,58,62,66,70,72,75,78,81,85,88,91,
     2   94,100,104,108,112,118,121,124,127,131,135,139,143,146,152,
     3   159,166,171,180,189,198,204,214,224,234,239,243,247,251,259,
     4   266,273,280,292,301,310,319,333,343,353,363/
      data lml/1,2,3,4,6,7,8,9,11,13,15,17,18,20,22,24,27,29,31,33,
     1   37,39,41,43,46,49,52,55,57,61,65,69,71,74,77,80,84,87,90,93,
     2   99,103,107,111,117,120,123,126,130,134,138,142,145,151,158,
     3   165,170,179,188,197,203,213,223,233,238,242,246,250,258,265,
     4   272,279,291,300,309,318,332,342,352,362,377/
      data lmx/0,0,0,1,0,0,0,1,1,2,0,0,0,0,0,1,1,1,2,0,2,0,3,1,0,0,0,
     1   0,0,1,1,1,1,2,0,2,0,2,0,3,1,3,1,4,2,0,0,0,0,0,0,0,1,1,1,1,1,
     2   2,0,2,0,2,0,2,0,3,1,3,1,3,1,4,2,0,4,2,0,5,3,1,0,0,0,0,0,0,0,
     3   1,1,1,1,1,1,2,0,2,0,2,0,2,0,2,0,3,1,3,1,3,1,3,1,4,2,0,4,2,0,
     4   4,2,0,5,3,1,5,3,1,6,4,2,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,
     5   2,0,0,0,0,0,0,0,0,0,0,1,3,1,1,1,3,3,1,1,1,1,1,2,2,2,4,4,0,0,
     6   0,0,0,0,0,0,0,0,0,0,0,1,3,5,1,1,1,3,3,5,1,1,1,1,1,1,2,2,2,4,
     7   4,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3,5,7,1,1,1,3,3,5,0,0,0,0,
     8   0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,
     9   1,1,3,3,3,1,1,1,1,1,1,1,2,2,2,2,4,4,4,0,0,0,0,0,0,0,0,0,0,0,
     9   0,0,0,1,1,1,1,3,3,3,5,5,1,1,1,1,1,1,1,1,1,0,2,4,6,2,2,2,4,4,
     9   6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,3,3,3,5,5,7,1,1,1,1,1,
     9   1,1,1,1,1,0,2,4,6,8,2,2,2,4,4,6,0,0,0,0/
      data lmy/0,0,1,0,0,0,1,0,1,0,2,0,0,1,1,0,0,1,0,2,1,3,0,2,0,0,0,
     1   1,1,0,0,1,1,0,2,0,2,1,3,0,2,1,3,0,2,4,0,0,0,1,1,1,0,0,0,1,1,
     2   0,2,0,2,1,3,1,3,0,2,0,2,1,3,0,2,4,1,3,5,0,2,4,0,0,0,0,1,1,1,
     3   0,0,0,1,1,1,0,2,0,2,0,2,1,3,1,3,0,2,0,2,1,3,1,3,0,2,4,0,2,4,
     4   1,3,5,0,2,4,1,3,5,0,2,4,6,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,0,0,
     5   0,2,2,2,1,3,1,1,1,3,3,0,0,0,0,0,0,0,1,3,1,3,1,0,0,0,0,0,0,0,
     6   0,0,1,3,5,1,1,1,3,3,5,0,0,0,0,0,0,0,0,0,1,1,1,3,3,5,0,0,0,0,
     7   0,0,0,0,0,0,1,3,5,7,1,1,1,3,3,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     8   0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,2,2,2,2,1,1,1,1,3,3,3,0,0,
     9   0,0,0,0,0,1,1,1,1,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,3,3,
     9   3,5,5,0,0,0,0,0,0,0,0,0,1,3,5,1,1,1,3,3,5,0,0,0,0,0,0,0,0,0,
     9   0,0,0,0,0,1,1,1,1,3,3,3,5,5,7,0,0,0,0,0,0,0,0,0,0,1,3,5,7,1,
     9   1,1,3,3,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data lmz/0,1,0,0,2,0,1,1,0,0,0,3,1,2,0,2,0,1,1,1,0,0,0,0,4,2,0,
     1   3,1,3,1,2,0,2,2,0,0,1,1,1,1,0,0,0,0,0,5,3,1,4,2,0,4,2,0,3,1,
     2   3,3,1,1,2,2,0,0,2,2,0,0,1,1,1,1,1,0,0,0,0,0,0,6,4,2,0,5,3,1,
     3   5,3,1,4,2,0,4,4,2,2,0,0,3,3,1,1,3,3,1,1,2,2,0,0,2,2,2,0,0,0,
     4   1,1,1,1,1,1,0,0,0,0,0,0,0,1,3,5,7,0,2,4,6,0,2,4,6,1,3,5,1,3,
     5   5,1,3,5,0,0,2,4,6,2,4,0,0,2,4,6,2,4,1,1,3,3,5,1,3,5,1,3,1,3,
     6   5,7,0,0,0,2,4,6,2,4,2,0,0,0,2,4,6,2,4,2,1,3,5,1,3,1,1,3,5,1,
     7   3,1,1,3,5,7,0,0,0,0,2,4,6,2,4,2,0,0,0,0,2,4,6,2,4,2,0,2,4,6,
     8   8,1,3,5,7,1,3,5,7,0,2,4,6,0,2,4,6,0,2,4,6,1,3,5,7,1,3,5,1,3,
     9   5,7,1,3,5,0,2,4,6,0,2,4,0,2,4,6,0,2,4,0,2,4,6,8,1,3,5,7,1,3,
     9   5,1,3,1,3,5,7,1,3,5,1,3,0,0,0,2,4,6,2,4,2,0,0,0,0,2,4,6,2,4,
     9   2,2,4,6,8,1,3,5,7,1,3,5,1,3,1,1,3,5,7,1,3,5,1,3,1,0,0,0,0,2,
     9   4,6,2,4,2,0,0,0,0,0,2,4,6,2,4,2,2,4,6,8/
      data a0 /0.0d0/, a1 /1.0d0/
c*    call dzero(ang,7*lamu)
      call dzero(ang,9*lamu)
      do 96 n=1,nanb
      if(xab(n).eq.a0) go to 96
      do 94 l=1,lalb
      if(yab(l).eq.a0) go to 94
      do 92 m=1,mamb
      if(zab(m).eq.a0) go to 92
      nlm=((n-2)+l)+m
      lamlo=mod(nlm-1,2)+1
      lamhi=min0(nlm,lamu)
      if(lamlo.gt.lamhi) go to 92
      do 90 lam=lamlo,lamhi,2
      l2=lam+lam-1
      angt=a0
      loc=lf(lam)
      do 80 mu1=1,l2
      istart=lmf(loc+mu1)
      if(mod(n,2).eq.mod(lmx(istart),2).or.
     1   mod(l,2).eq.mod(lmy(istart),2).or.
     2   mod(m,2).eq.mod(lmz(istart),2)) go to 80
      pre=a0
      aint=a0
      iend=lml(loc+mu1)
      do 70 i=istart,iend
      indx=lmx(i)
      indy=lmy(i)
      indz=lmz(i)
      if(indx.ne.0) go to 10
      xkp=a1
      go to 20
   10 xkp=xk**indx
   20 if(indy.ne.0) go to 30
      ykp=a1
      go to 40
   30 ykp=yk**indy
   40 if(indz.ne.0) go to 50
      zkp=a1
      go to 60
   50 zkp=zk**indz
   60 pre=pre+zlm(i)*xkp*ykp*zkp
   70 aint=aint+zlm(i)*dfac(n+indx)*dfac(l+indy)*dfac(m+indz)/
     1   dfac((n+indx)+(l+indy)+(m+indz))
      angt=angt+pre*aint
   80 continue
   90 ang(nlm,lam)=ang(nlm,lam)+((xab(n)*yab(l))*zab(m))*angt
   92 continue
   94 continue
   96 continue
      return
      end
      subroutine ang2(it,crda,l,lmlo,lmhi,xk,yk,zk,ang)
c     ----- computes angular integrals for type 2 integrals -----
c**   ge#ndert durch h. stoll, oktober 1987
      implicit double precision (a-h,o-z)
      common /zlmtab/ zlm(377)
c*    common /zlmtab/ zlm(130)
      common /dfac/ dfac(29)
c*    common /dfac/ dfac(23)
c*    dimension ang(4,7,7), crda(4,3)
      dimension ang(5,9,9), crda(5,3)
c*    dimension lf(7),lmf(49),lml(49),lmx(130),lmy(130),lmz(130)
      dimension lf(9),lmf(81),lml(81),lmx(377),lmy(377),lmz(377)
c*    dimension binom(10), ind(4), nap(35), lap(35), map(35)
      dimension binom(15), ind(5), nap(35), lap(35), map(35)
c*    data lf /0,1,4,9,16,25,36/
      data lf /0,1,4,9,16,25,36,49,64/
c*    data lmf/1,2,3,4,5,7,8,9,10,12,14,16,18,19,21,23,25,28,30,32,
c*   1   34,38,40,42,44,47,50,53,56,58,62,66,70,72,75,78,81,85,88,91,
c*   2   94,100,104,108,112,118,121,124,127/
c*    data lml/1,2,3,4,6,7,8,9,11,13,15,17,18,20,22,24,27,29,31,33,
c*   1   37,39,41,43,46,49,52,55,57,61,65,69,71,74,77,80,84,87,90,93,
c*   2   99,103,107,111,117,120,123,126,130/
c*    data lmx/0,0,0,1,0,0,0,1,1,2,0,0,0,0,0,1,1,1,2,0,2,0,3,1,0,0,0,
c*   1   0,0,1,1,1,1,2,0,2,0,2,0,3,1,3,1,4,2,0,0,0,0,0,0,0,1,1,1,1,1,
c*   2   2,0,2,0,2,0,2,0,3,1,3,1,3,1,4,2,0,4,2,0,5,3,1,0,0,0,0,0,0,0,
c*   3   1,1,1,1,1,1,2,0,2,0,2,0,2,0,2,0,3,1,3,1,3,1,3,1,4,2,0,4,2,0,
c*   4   4,2,0,5,3,1,5,3,1,6,4,2,0/
c*    data lmy/0,0,1,0,0,0,1,0,1,0,2,0,0,1,1,0,0,1,0,2,1,3,0,2,0,0,0,
c*   1   1,1,0,0,1,1,0,2,0,2,1,3,0,2,1,3,0,2,4,0,0,0,1,1,1,0,0,0,1,1,
c*   2   0,2,0,2,1,3,1,3,0,2,0,2,1,3,0,2,4,1,3,5,0,2,4,0,0,0,0,1,1,1,
c*   3   0,0,0,1,1,1,0,2,0,2,0,2,1,3,1,3,0,2,0,2,1,3,1,3,0,2,4,0,2,4,
c*   4   1,3,5,0,2,4,1,3,5,0,2,4,6/
c*    data lmz/0,1,0,0,2,0,1,1,0,0,0,3,1,2,0,2,0,1,1,1,0,0,0,0,4,2,0,
c*   1   3,1,3,1,2,0,2,2,0,0,1,1,1,1,0,0,0,0,0,5,3,1,4,2,0,4,2,0,3,1,
c*   2   3,3,1,1,2,2,0,0,2,2,0,0,1,1,1,1,1,0,0,0,0,0,0,6,4,2,0,5,3,1,
c*   3   5,3,1,4,2,0,4,4,2,2,0,0,3,3,1,1,3,3,1,1,2,2,0,0,2,2,2,0,0,0,
c*   4   1,1,1,1,1,1,0,0,0,0,0,0,0/
      data lmf/1,2,3,4,5,7,8,9,10,12,14,16,18,19,21,23,25,28,30,32,
     1   34,38,40,42,44,47,50,53,56,58,62,66,70,72,75,78,81,85,88,91,
     2   94,100,104,108,112,118,121,124,127,131,135,139,143,146,152,
     3   159,166,171,180,189,198,204,214,224,234,239,243,247,251,259,
     4   266,273,280,292,301,310,319,333,343,353,363/
      data lml/1,2,3,4,6,7,8,9,11,13,15,17,18,20,22,24,27,29,31,33,
     1   37,39,41,43,46,49,52,55,57,61,65,69,71,74,77,80,84,87,90,93,
     2   99,103,107,111,117,120,123,126,130,134,138,142,145,151,158,
     3   165,170,179,188,197,203,213,223,233,238,242,246,250,258,265,
     4   272,279,291,300,309,318,332,342,352,362,377/
      data lmx/0,0,0,1,0,0,0,1,1,2,0,0,0,0,0,1,1,1,2,0,2,0,3,1,0,0,0,
     1   0,0,1,1,1,1,2,0,2,0,2,0,3,1,3,1,4,2,0,0,0,0,0,0,0,1,1,1,1,1,
     2   2,0,2,0,2,0,2,0,3,1,3,1,3,1,4,2,0,4,2,0,5,3,1,0,0,0,0,0,0,0,
     3   1,1,1,1,1,1,2,0,2,0,2,0,2,0,2,0,3,1,3,1,3,1,3,1,4,2,0,4,2,0,
     4   4,2,0,5,3,1,5,3,1,6,4,2,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,
     5   2,0,0,0,0,0,0,0,0,0,0,1,3,1,1,1,3,3,1,1,1,1,1,2,2,2,4,4,0,0,
     6   0,0,0,0,0,0,0,0,0,0,0,1,3,5,1,1,1,3,3,5,1,1,1,1,1,1,2,2,2,4,
     7   4,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,3,5,7,1,1,1,3,3,5,0,0,0,0,
     8   0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0,1,1,
     9   1,1,3,3,3,1,1,1,1,1,1,1,2,2,2,2,4,4,4,0,0,0,0,0,0,0,0,0,0,0,
     9   0,0,0,1,1,1,1,3,3,3,5,5,1,1,1,1,1,1,1,1,1,0,2,4,6,2,2,2,4,4,
     9   6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,3,3,3,5,5,7,1,1,1,1,1,
     9   1,1,1,1,1,0,2,4,6,8,2,2,2,4,4,6,0,0,0,0/
      data lmy/0,0,1,0,0,0,1,0,1,0,2,0,0,1,1,0,0,1,0,2,1,3,0,2,0,0,0,
     1   1,1,0,0,1,1,0,2,0,2,1,3,0,2,1,3,0,2,4,0,0,0,1,1,1,0,0,0,1,1,
     2   0,2,0,2,1,3,1,3,0,2,0,2,1,3,0,2,4,1,3,5,0,2,4,0,0,0,0,1,1,1,
     3   0,0,0,1,1,1,0,2,0,2,0,2,1,3,1,3,0,2,0,2,1,3,1,3,0,2,4,0,2,4,
     4   1,3,5,0,2,4,1,3,5,0,2,4,6,0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,0,0,
     5   0,2,2,2,1,3,1,1,1,3,3,0,0,0,0,0,0,0,1,3,1,3,1,0,0,0,0,0,0,0,
     6   0,0,1,3,5,1,1,1,3,3,5,0,0,0,0,0,0,0,0,0,1,1,1,3,3,5,0,0,0,0,
     7   0,0,0,0,0,0,1,3,5,7,1,1,1,3,3,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     8   0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,0,2,2,2,2,1,1,1,1,3,3,3,0,0,
     9   0,0,0,0,0,1,1,1,1,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,3,3,
     9   3,5,5,0,0,0,0,0,0,0,0,0,1,3,5,1,1,1,3,3,5,0,0,0,0,0,0,0,0,0,
     9   0,0,0,0,0,1,1,1,1,3,3,3,5,5,7,0,0,0,0,0,0,0,0,0,0,1,3,5,7,1,
     9   1,1,3,3,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      data lmz/0,1,0,0,2,0,1,1,0,0,0,3,1,2,0,2,0,1,1,1,0,0,0,0,4,2,0,
     1   3,1,3,1,2,0,2,2,0,0,1,1,1,1,0,0,0,0,0,5,3,1,4,2,0,4,2,0,3,1,
     2   3,3,1,1,2,2,0,0,2,2,0,0,1,1,1,1,1,0,0,0,0,0,0,6,4,2,0,5,3,1,
     3   5,3,1,4,2,0,4,4,2,2,0,0,3,3,1,1,3,3,1,1,2,2,0,0,2,2,2,0,0,0,
     4   1,1,1,1,1,1,0,0,0,0,0,0,0,1,3,5,7,0,2,4,6,0,2,4,6,1,3,5,1,3,
     5   5,1,3,5,0,0,2,4,6,2,4,0,0,2,4,6,2,4,1,1,3,3,5,1,3,5,1,3,1,3,
     6   5,7,0,0,0,2,4,6,2,4,2,0,0,0,2,4,6,2,4,2,1,3,5,1,3,1,1,3,5,1,
     7   3,1,1,3,5,7,0,0,0,0,2,4,6,2,4,2,0,0,0,0,2,4,6,2,4,2,0,2,4,6,
     8   8,1,3,5,7,1,3,5,7,0,2,4,6,0,2,4,6,0,2,4,6,1,3,5,7,1,3,5,1,3,
     9   5,7,1,3,5,0,2,4,6,0,2,4,0,2,4,6,0,2,4,0,2,4,6,8,1,3,5,7,1,3,
     9   5,1,3,1,3,5,7,1,3,5,1,3,0,0,0,2,4,6,2,4,2,0,0,0,0,2,4,6,2,4,
     9   2,2,4,6,8,1,3,5,7,1,3,5,1,3,1,1,3,5,7,1,3,5,1,3,1,0,0,0,0,2,
     9   4,6,2,4,2,0,0,0,0,0,2,4,6,2,4,2,2,4,6,8/
      data nap /0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,
     1           4,0,0,3,3,1,0,1,0,2,2,0,2,1,1/,
     2     lap /0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,
     3           0,4,0,1,0,3,3,0,1,2,0,2,1,2,1/,
     4     map /0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1,
     5           0,0,4,0,1,0,1,3,3,0,2,2,1,1,2/
c*    data binom/4*1.0d0,2.0d0,2*1.0d0,2*3.0d0,1.0d0/
c*    data ind/0,1,3,6/
      data binom/4*1.0d0,2.0d0,2*1.0d0,2*3.0d0,1.0d0,1.0d0,4.0d0,
     a           6.0d0,4.0d0,1.0d0/
      data ind/0,1,3,6,10/
      data a0 /0.0d0/, a1 /1.0d0/
c*    call dzero(ang,28*lmhi)
      call dzero(ang,45*lmhi)
      na1=nap(it)+1
      la1=lap(it)+1
      ma1=map(it)+1
      naind=ind(na1)
      laind=ind(la1)
      maind=ind(ma1)
      loc1=lf(l)
      mhi=l+l-1
      do 85 ia=1,na1
      pab1=binom(naind+ia)*crda((na1+1)-ia,1)
      if(pab1.eq.a0) go to 85
      do 80 ib=1,la1
      pab2=pab1*binom(laind+ib)*crda((la1+1)-ib,2)
      if(pab2.eq.a0) go to 80
      do 75 ic=1,ma1
      pab3=pab2*binom(maind+ic)*crda((ma1+1)-ic,3)
      if(pab3.eq.a0) go to 75
      n=((ia-3)+ib)+ic
c**   lamlo=max0(l-n,lmlo+mod(l+n+lmlo,2))
c**   lamhi=min0(l+n,lmhi-mod(l+n+lmhi,2))
      lamlo=max0(l-n,lmlo)
      lamhi=min0(l+n,lmhi)
      if(mod(lamhi-lamlo,2).ne.0) lamlo=lamlo+1
c**
      if(lamlo.gt.lamhi) go to 75
      do 20 m=1,mhi
      mstart=lmf(loc1+m)
      mend=lml(loc1+m)
      do 20 lam=lamlo,lamhi,2
      l2=lam+lam-1
      angt=a0
      loc2=lf(lam)
      do 30 mu=1,l2
      istart=lmf(loc2+mu)
      if(mod(ia+lmx(mstart)+lmx(istart),2).ne.1.or.
     1   mod(ib+lmy(mstart)+lmy(istart),2).ne.1.or.
     2   mod(ic+lmz(mstart)+lmz(istart),2).ne.1) go to 30
      pre=a0
      iend=lml(loc2+mu)
      aint=a0
      do 40 i=istart,iend
      indx=lmx(i)
      indy=lmy(i)
      indz=lmz(i)
      if(indx.ne.0)go to 4
      xkp=a1
      go to 5
 4    xkp=xk**indx
 5    if(indy.ne.0)go to 6
      ykp=a1
      go to 7
 6    ykp=yk**indy
 7    if(indz.ne.0)go to 8
      zkp=a1
      go to 9
 8    zkp=zk**indz
 9    pre=pre+zlm(i)*xkp*ykp*zkp
      do 40 j=mstart,mend
      mndx=lmx(j)
      mndy=lmy(j)
      mndz=lmz(j)
   40 aint=aint+zlm(i)*zlm(j)*dfac(ia+indx+mndx)*dfac(ib+indy+mndy)*
     1   dfac(ic+indz+mndz)/dfac(ia+indx+mndx+ib+indy+mndy+ic+indz+mndz)
      angt=angt+pre*aint
   30 continue
c**20 ang(n+1,m,lam)=angt*pab3
   20 ang(n+1,m,lam)=ang(n+1,m,lam)+angt*pab3
c**
   75 continue
   80 continue
   85 continue
      return
      end
      subroutine facab(na1,nb1,crda,crdb,xab)
      implicit double precision (a-h,o-z)
c*    dimension crda(2), crdb(2), bic(10), ind(4), xab(2)
      dimension crda(*), crdb(*), bic(15), ind(5), xab(*)
c*    data bic/4*1.0d0,2.0d0,2*1.0d0,2*3.0d0,1.0d0/
      data bic/4*1.0d0,2.0d0,2*1.0d0,2*3.0d0,1.0d0,1.0d0,4.0d0,6.0d0,
     a         4.0d0,1.0d0/
c*    data ind/0,1,3,6/
      data ind/0,1,3,6,10/
      call dzero(xab,na1+nb1-1)
      naind=ind(na1)
      nbind=ind(nb1)
      do 10 ia1=1,na1
      do 10 ib1=1,nb1
   10 xab((ia1-1)+ib1)=xab((ia1-1)+ib1)+(bic(naind+ia1)*
     1   crda((na1+1)-ia1))*bic(nbind+ib1)*crdb((nb1+1)-ib1)
      return
      end
      subroutine ptwt(npi,l,lit,ljt,ltot1,lmalo,lmahi,lmblo,lmbhi,
     1   alpha,rc,rka,rkb,arc2,prd,qsum)
c     ----- computes q(n,la,lb) using points and weights method -----
c     ----- scaled by exp(-arc2) to prevent overflows           -----
      implicit double precision (a-h,o-z)
c*    dimension qsum(7,7,2)
c*    dimension abess(20,7), bbess(20,7), ptpow(20,7), w(35), z(35)
      dimension qsum(9,9,9)
      dimension abess(20,9), bbess(20,9), ptpow(20,9), w(35), z(35)
      data              z /-.20201828704561d+01,-.95857246461382d+00,
     + .00000000000000d+00, .95857246461382d+00, .20201828704561d+01,
     +                     -.34361591188377d+01,-.25327316742328d+01,
     +-.17566836492999d+01,-.10366108297895d+01,-.34290132722370d+00,
     + .34290132722370d+00, .10366108297895d+01, .17566836492999d+01,
     + .25327316742328d+01, .34361591188377d+01,
     +                     -.53874808900112d+01,-.46036824495507d+01,
     +-.39447640401156d+01,-.33478545673832d+01,-.27888060584281d+01,
     +-.22549740020893d+01,-.17385377121166d+01,-.12340762153953d+01,
     +-.73747372854539d+00,-.24534070830090d+00, .24534070830090d+00,
     + .73747372854539d+00, .12340762153953d+01, .17385377121166d+01,
     + .22549740020893d+01, .27888060584281d+01, .33478545673832d+01,
     + .39447640401156d+01, .46036824495507d+01, .53874808900112d+01/
      data              w / .19953242059046d-01, .39361932315224d+00,
     + .94530872048294d+00, .39361932315224d+00, .19953242059046d-01,
     +                      .76404328552326d-05, .13436457467812d-02,
     + .33874394455481d-01, .24013861108231d+00, .61086263373533d+00,
     + .61086263373533d+00, .24013861108231d+00, .33874394455481d-01,
     + .13436457467812d-02, .76404328552326d-05,
     +                      .22293936455342d-12, .43993409922732d-09,
     + .10860693707693d-06, .78025564785321d-05, .22833863601635d-03,
     + .32437733422379d-02, .24810520887464d-01, .10901720602002d+00,
     + .28667550536283d+00, .46224366960061d+00, .46224366960061d+00,
     + .28667550536283d+00, .10901720602002d+00, .24810520887464d-01,
     + .32437733422379d-02, .22833863601635d-03, .78025564785321d-05,
     + .10860693707693d-06, .43993409922732d-09, .22293936455342d-12/
      data a500 /500.0d0/, a50000 /50000.0d0/
      if(arc2.gt.a50000) go to 20
      if(arc2.gt.a500) go to 10
      npt=20
      idif=15
      go to 30
   10 npt=10
      idif=5
      go to 30
   20 npt=5
      idif=0
   30 sqalp=dsqrt(alpha)
      prd=prd/sqalp
      do 90 i=1,npt
      pt=rc+z(i+idif)/sqalp
      do 50 lama=lmalo,lmahi
   50 abess(i,lama)=bess(rka*pt,lama-1)
      do 60 lamb=lmblo,lmbhi
   60 bbess(i,lamb)=bess(rkb*pt,lamb-1)
      if(npi.gt.0) go to 65
      ptpow(i,1)=prd
      go to 68
   65 ptpow(i,1)=prd*pt**npi
   68 if(ltot1.eq.1) go to 72
      do 70 n=2,ltot1
   70 ptpow(i,n)=pt*ptpow(i,n-1)
   72 do 90 lama=lmalo,lmahi
      ldifa1=iabs(l-lama)+1
      do 90 lamb=lmblo,lmbhi
      ldifb=iabs(l-lamb)
      nlo=ldifa1+ldifb
      nhi=(ltot1-mod(lit-ldifa1,2))-mod((ljt-1)-ldifb,2)
      do 90 n=nlo,nhi,2
   90 qsum(n,lamb,lama)=qsum(n,lamb,lama)+((w(i+idif)*abess(i,lama))*
     1   bbess(i,lamb))*ptpow(i,n)
      return
      end
      subroutine qpasy(npi,l,lit,ljt,ltot1,lmalo,lmahi,lmblo,lmbhi,
     1   alp,xka1,xkb1,prd,dum,qsum)
c     ----- partially asymptotic form for q(n,la1,lb1) -----
c     ----- scaled by exp(-xkb*xkb/(4*alpha)) to prevent overflows
      implicit double precision (a-h,o-z)
c*    common/qstore/dum1(49),alpha,xk,t
      common/qstore/dum1(81),alpha,xk,t
      common /dfac/ dfac(29)
c*    common /dfac/ dfac(23)
c*    dimension qsum(7,7,2)
c*    dimension fctr(7), sum(7), term(7)
      dimension qsum(9,9,9)
      dimension fctr(9), sum(9), term(9)
      data a0,accrcy,a1s4,a1s2,a1 /0.0d0,1.0d-13,0.25d0,0.5d0,1.0d0/
c     ----- first set up xkb as larger -----
      save
      if(xka1.gt.xkb1) go to 10
      xka=xka1
      xkb=xkb1
      go to 12
   10 xka=xkb1
      xkb=xka1
c     ----- set up parameters for qcomp using xkb -----
   12 alpha=a1
      sqalp=dsqrt(alp)
      xk=xkb/sqalp
      t=a1s4*xk*xk
      prd=prd*dexp(-(dum-t))
      tk=xka*xka/(alp+alp)
      do 90 lama=lmalo,lmahi
      ldifa1=iabs(l-lama)+1
      if(xka1.gt.xkb1) go to 14
      la=lama-1
      go to 16
   14 lb=lama-1
   16 do 90 lamb=lmblo,lmbhi
      ldifb=iabs(l-lamb)
      nlo=ldifa1+ldifb
      nhi=(ltot1-mod(lit-ldifa1,2))-mod((ljt-1)-ldifb,2)
      if(xka1.gt.xkb1) go to 18
      lb=lamb-1
      go to 20
   18 la=lamb-1
c     ----- run power series using xka, obtaining initial    -----
c     ----- q(n,l) values from qcomp, then recurring upwards -----
c     ----- j=0 term in sum -----
   20 qold2=qcomp(npi+nlo-1+la,lb)/dfac(la+la+3)
      fctr(nlo)=a1
      sum(nlo)=qold2
      if(nlo.eq.nhi.and.tk.eq.a0) go to 60
c     ----- j=1 term in sum -----
      nprime=npi+nlo+la+1
      qold1=qcomp(nprime,lb)/dfac(la+la+3)
      if(nlo.ne.nhi) fctr(nlo+2)=fctr(nlo)
      f1=(la+la+3)
      fctr(nlo)=tk/f1
      term(nlo)=fctr(nlo)*qold1
      sum(nlo)=sum(nlo)+term(nlo)
      if(nlo.ne.nhi) go to 22
      qold2=fctr(nlo)*qold2
      qold1=term(nlo)
      go to 24
   22 nlo2=nlo+2
      sum(nlo2)=qold1
      if(nlo2.eq.nhi.and.tk.eq.a0) go to 60
   24 j=1
c     ----- increment j for next term -----
   30 j=j+1
      nprime=nprime+2
      f1=(nprime+nprime-5)
      f2=((lb-nprime+4)*(lb+nprime-3))
      qnew=(t+a1s2*f1)*qold1+a1s4*f2*qold2
      nlojj=nlo+j+j
      if(nlo.eq.nhi) go to 40
      nhitmp=min0(nlojj,nhi)
      do 38 n=nlo2,nhitmp,2
      nrev=nhitmp+nlo2-n
   38 fctr(nrev)=fctr(nrev-2)
   40 f1=(j*(la+la+j+j+1))
      fctr(nlo)=tk/f1
      if(nlojj.gt.nhi) go to 44
      nhitmp=nlojj-2
      term(nlojj)=qnew
      sum(nlojj)=term(nlojj)
      do 42 n=nlo,nhitmp,2
      nrev=nhitmp+nlo-n
      term(nrev)=fctr(nrev)*term(nrev+2)
   42 sum(nrev)=sum(nrev)+term(nrev)
      if(nlojj.eq.nhi.and.tk.eq.a0) go to 60
      qold2=qold1
      qold1=qnew
      go to 30
   44 qold2=fctr(nhi)*qold1
      qold1=fctr(nhi)*qnew
      term(nhi)=qold1
      sum(nhi)=sum(nhi)+term(nhi)
      if(nlo.eq.nhi) go to 47
      nhitmp=nhi-2
      do 46 n=nlo,nhitmp,2
      nrev=nhitmp+nlo-n
      term(nrev)=fctr(nrev)*term(nrev+2)
   46 sum(nrev)=sum(nrev)+term(nrev)
   47 do 48 n=nlo,nhi,2
   48 if(term(n).gt.accrcy*sum(n)) go to 30
   60 if(la.ne.0) go to 62
      prefac=prd/sqalp**(npi+nlo+la)
      go to 64
   62 prefac=prd*xka**la/sqalp**(npi+nlo+la)
   64 do 66 n=nlo,nhi,2
      qsum(n,lamb,lama)=qsum(n,lamb,lama)+prefac*sum(n)
   66 prefac=prefac/alp
   90 continue
      return
      end
      subroutine recur1(npi,ltot1)
c     ----- controls computation of q(n,l) by recurrence -----
      implicit double precision (a-h,o-z)
      common /pifac/ pi, sqpi2
c*    common/qstore/q(7,7),alpha,xk,t
      common/qstore/q(9,9),alpha,xk,t
      data a0, a1s2, a1, a2, a3 /0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0/
      np1l=mod(npi,2)+1
      np1u=npi+1
      if(xk.ne.a0) go to 24
c     ----- xk=0 -----
      f1sta=a1/alpha
      if(np1l.eq.2) go to 12
      q(1,1)=a1s2*dsqrt(pi/alpha)
      fista=a1s2*f1sta
      go to 14
   12 q(1,1)=a1s2*f1sta
      fista=f1sta
   14 if(np1l.eq.np1u) go to 18
      np1l=np1l+2
      do 16 np1=np1l,np1u,2
      q(1,1)=fista*q(1,1)
   16 fista=fista+f1sta
   18 if(ltot1.le.2) return
      do 20 lp1=3,ltot1,2
      q(lp1,1)=fista*q(lp1-2,1)
   20 fista=fista+f1sta
      return
   24 if(ltot1-2) 28,26,30
c     ----- ltot1 1 or 2 -----
   26 q(2,2)=qcomp(npi+1,1)
   28 q(1,1)=qcomp(npi,0)
      return
   30 talph=alpha+alpha
      if(npi.gt.0) go to 32
      n2=3
      go to 34
   32 n2=1
   34 if(np1l.eq.2) go to 48
c     ----- npi even -----
      q(n2,1)=qcomp(2,0)
      lu=ltot1-n2
      if(lu.eq.0) go to 38
c     ----- use (38a) -----
      do 36 l=1,lu
   36 q(n2+l,l+1)=(xk/talph)*q(n2+l-1,l)
      if(npi-2) 38,68,60
c     ----- npi=0 -----
   38 if(t.gt.a3) go to 42
c     ----- recur down using (38c) etc. -----
      q(ltot1,ltot1)=qcomp(ltot1-1,ltot1-1)
      f2lp1=(ltot1+ltot1-3)
      q(ltot1-1,ltot1-1)=xk*(q(ltot1,ltot1-2)-q(ltot1,ltot1))/f2lp1
      do 40 lr=3,ltot1
      l=ltot1-lr
      f2lp1=f2lp1-a2
   40 q(l+1,l+1)=(talph*q(l+3,l+1)-xk*q(l+2,l+2))/f2lp1
      if(ltot1-4) 74,74,68
c     ----- recur up using (38b) etc. -----
   42 q(1,1)=qcomp(0,0)
      f2lm1=a1
      do 44 lp2=3,ltot1
      q(lp2-1,lp2-1)=(talph*q(lp2,lp2-2)-f2lm1*q(lp2-2,lp2-2))/xk
   44 f2lm1=f2lm1+a2
      q(ltot1,ltot1)=q(ltot1,ltot1-2)-f2lm1*q(ltot1-1,ltot1-1)/xk
      if(ltot1-4) 74,74,68
c     ----- npi odd -----
   48 if(t.gt.a3) go to 52
c     ----- recur down using (38f) -----
      q(ltot1,ltot1)=qcomp(ltot1,ltot1-1)
      q(ltot1-1,ltot1-1)=qcomp(ltot1-1,ltot1-2)
      f2lp2=(ltot1+ltot1-2)
      do 50 lr=3,ltot1
      l=ltot1-lr
      f2lp2=f2lp2-a2
   50 q(l+1,l+1)=(talph*q(l+3,l+3)-(xk-(f2lp2+a1)*(talph/xk))*
     1   q(l+2,l+2))/f2lp2
      go to 56
c     ----- recur up using (38e) -----
   52 q(1,1)=qcomp(1,0)
      q(2,2)=qcomp(2,1)
      f2lm2=a0
      do 54 lp1=3,ltot1
      f2lm2=f2lm2+a2
   54 q(lp1,lp1)=(f2lm2*q(lp1-2,lp1-2)+(xk-(f2lm2+a1)*(talph/xk))*
     1   q(lp1-1,lp1-1))/talph
   56 if(npi.eq.1) go to 68
c     ----- npi greater than 2.  use (38d) etc. -----
   60 np1l=mod(npi-3,2)+4
      fnm3=(np1l-6)
      do 64 np1=np1l,np1u,2
      fnm3=fnm3+a2
      fnplm1=fnm3
      do 62 lp2=2,ltot1
      fnplm1=fnplm1+a2
   62 q(lp2-1,lp2-1)=(fnplm1*q(lp2-1,lp2-1)+xk*q(lp2,lp2))/talph
   64 q(ltot1,ltot1)=((fnm3+a1)*q(ltot1,ltot1)+xk*
     1   q(ltot1-1,ltot1-1))/talph
c     ----- use (38d) -----
   68 nbnd=(ltot1-n2)/2
      lp1mx=(ltot1-n2)+1
      fnm3=(npi+n2-4)
      do 70 ibnd=1,nbnd
      fnm3=fnm3+a2
      fnplm1=fnm3
      n2=n2+2
      lp1mx=lp1mx-2
      do 70 lp1=1,lp1mx
      fnplm1=fnplm1+a2
   70 q(n2+lp1-1,lp1)=(fnplm1*q(n2+lp1-3,lp1)+xk*q(n2+lp1-2,lp1+1))/
     1   talph
   74 return
      end
      function qcomp(n,l)
c     ----- computes q(n,l)                                -----
c     ----- scaled by exp(-t) to prevent overflows         -----
c     ----- arguments are alpha, xk, and t=xk**2/(4*alpha) -----
c     ----- no restriction on the magnitude of t           -----
c     ----- increase dfac array to raise n, l restrictions -----
      implicit double precision (a-h,o-z)
      common /pifac/ pi, sqpi2
      common /dfac/ dfac(29)
c*    common /dfac/ dfac(23)
c*    common/qstore/dum1(49),alpha,xk,t
      common/qstore/dum1(81),alpha,xk,t
      dimension tmin(9)
      data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,
     1   18.0d0,15.0d0/
      data am1,   a0,   accpow, accasy, a1,   a2,   a4
     1    /-1.0d0,0.0d0,1.0d-14,1.0d-10,1.0d0,2.0d0,4.0d0/
      if(mod(n+l,2).ne.0.or.n.le.l) go to 30
c     ----- use alternating series (n+l.le.22.and.l.le.10) -----
      if(l.eq.0) xkp=a1
      if(l.ne.0) xkp=(xk/(alpha+alpha))**l
      prefac=sqpi2*xkp*dfac(n+l+1)/
     1   ((alpha+alpha)**((n-l)/2)*dsqrt(alpha+alpha)*dfac(l+l+3))
      num=l-n+2
      xden=(l+l+3)
      term=a1
      sum=term
      xc=am1
   10 if(num.eq.0) go to 20
      fnum=num
      term=term*fnum*t/(xden*xc)
      xc=xc+am1
      sum=sum+term
      num=num+2
      xden=xden+a2
      go to 10
   20 qcomp=prefac*sum
      return
   30 if(t.lt.tmin(min0(n,8)+1)) go to 60
c     ----- use asymptotic series (arbitrary n,l) -----
      xkp=(xk/(alpha+alpha))**(n-2)
      prefac=xkp*sqpi2/((alpha+alpha)*dsqrt(alpha+alpha))
      sum=a1
      term=a1
      fac1=(l-n+2)
      fac2=(1-l-n)
      xc=a1
   40 term=term*fac1*fac2/(a4*xc*t)
      if(term.eq.a0) go to 50
      sum=sum+term
      if(dabs(term/sum).lt.accasy) go to 50
      fac1=fac1+a2
      fac2=fac2+a2
      xc=xc+a1
      go to 40
   50 qcomp=prefac*sum
      return
c     ----- use power series (n+l.le.22.and.l.le.10) -----
   60 if(l.eq.0) xkp=a1
      if(l.ne.0) xkp=(xk/(alpha+alpha))**l
      prefac=dexp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
      if(mod(n+l,2).eq.0) prefac=prefac*sqpi2/dsqrt(alpha+alpha)
      xnum=(l+n-1)
      xden=(l+l+1)
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=a0
   70 xnum=xnum+a2
      xden=xden+a2
      xj=xj+a1
      term=term*t*xnum/(xj*xden)
      sum=sum+term
      if((term/sum).gt.accpow) go to 70
      qcomp=prefac*sum
      return
      end
      function bess(z,l)
      implicit double precision (a-h,o-z)
      common /dfac/ dfac(29)
c*    common /dfac/ dfac(23)
      common/fact/fac(17),fprod(9,9)
      data am1,a0,accrcy,a1s2,a1,a5,a16p1
     1   /-1.0d0,0.0d0,5.0d-14,0.5d0,1.0d0,5.0d0,16.1d0/
      if(z.gt.a5) go to 50
      if(z.eq.a0) go to 40
      if(z.lt.a0) go to 35
      zp=a1s2*z*z
      term=(z**l)/dfac(l+l+3)
      bess=term
      j=0
 5    j=j+1
      fjlj1=(j*(l+l+j+j+1))
      term=term*zp/fjlj1
      bess=bess+term
      if(dabs(term/bess).gt.accrcy) go to 5
      bess=bess*dexp(-z)
      go to 100
 35   bess=a0
      go to 100
 40   if(l.ne.0)go to 45
      bess=a1
      go to 100
 45   bess=a0
      go to 100
 50   if(z.gt.a16p1) go to 60
      rp=a0
      rm=a0
      tzp=z+z
      tzm=-tzp
      l1=l+1
      do 55 k1=1,l1
      k=k1-1
      rp=rp+fprod(k1,l1)/tzp**k
      rm=rm+fprod(k1,l1)/tzm**k
 55   continue
      bess=(rm-(am1**l)*rp*dexp(tzm))/tzp
      go to 100
 60   rm=a0
      tzm=-z-z
      l1=l+1
      do 65 k1=1,l1
      k=k1-1
      rm=rm+fprod(k1,l1)/tzm**k
 65   continue
      bess=rm/(-tzm)
 100  return
      end
      subroutine cortab
      implicit double precision (a-h,o-z)
      common /pifac/ pi, sqpi2
      common /dfac/ dfac(29)
c*    common /dfac/ dfac(23)
c*    common/fact/fac(13),fprod(7,7)
      common/fact/fac(17),fprod(9,9)
c*    common /zlmtab/ zlm(130)
      common /zlmtab/ zlm(377)
      data a0, a1s2, a1, a2, a3 /0.0d0, 0.5d0, 1.0d0, 2.0d0, 3.0d0/
c     ----- this value of pi may have to be changed if the     -----
c     ----- compiler will not accept extra significant figures -----
      pi=3.141592653589793238462643383279502884197169399d0
      sqpi2=sqrt(a1s2*pi)
      rt3=sqrt(a3)
      rt6=sqrt(6.0d0)
      rt5s2=sqrt(2.5d0)
      rt3s2=sqrt(1.5d0)
      dfac(1)=a1
      dfac(2)=a1
      fi=a0
c*    do 10 i=1,21
      do 10 i=1,27
      fi=fi+a1
   10 dfac(i+2)=fi*dfac(i)
      fac(1)=a1
      fi=a0
c*    do 20 i=1,12
      do 20 i=1,16
      fi=fi+a1
   20 fac(i+1)=fi*fac(i)
c*    do 30 l1=1,7
      do 30 l1=1,9
      do 30 k1=1,l1
   30 fprod(k1,l1)=fac(l1+(k1-1))/(fac(k1)*fac(l1-(k1-1)))
c     ----- real spherical harmonics in the form of   -----
c     ----- linear combinations of cartesian products -----
c     ----- l=0 -----
      zlm(1)=a1
c     ----- l=1 -----
      zlm(2)=rt3
      zlm(3)=zlm(2)
      zlm(4)=zlm(2)
c     ----- l=2 -----
      zlm(5)=a3*(a1s2*sqrt(5.0d0))
      zlm(6)=-(a1s2*sqrt(5.0d0))
      zlm(7)=sqrt(15.0d0)
      zlm(8)=zlm(7)
      zlm(9)=zlm(7)
      zlm(10)=a1s2*zlm(7)
      zlm(11)=-zlm(10)
c     ----- l=3 -----
      zlm(12)=5.0d0*(a1s2*sqrt(7.0d0))
      zlm(13)=-a3*(a1s2*sqrt(7.0d0))
      zlm(14)=5.0d0*(sqrt(21.0d0/8.0d0))
      zlm(15)=-(sqrt(21.0d0/8.0d0))
      zlm(16)=zlm(14)
      zlm(17)=zlm(15)
      zlm(18)=sqrt(105.0d0)
      zlm(19)=a1s2*zlm(18)
      zlm(20)=-zlm(19)
      zlm(21)=a3*(sqrt(35.0d0/8.0d0))
      zlm(22)=-(sqrt(35.0d0/8.0d0))
      zlm(23)=-zlm(22)
      zlm(24)=-zlm(21)
c     ----- l=4 -----
      zlm(25)=105.0d0/8.0d0
      zlm(26)=-45.0d0/4.0d0
      zlm(27)=9.0d0/8.0d0
      zlm(28)=21.0d0*(sqrt(5.0d0/8.0d0))
      zlm(29)=-9.0d0*(sqrt(5.0d0/8.0d0))
      zlm(30)=zlm(28)
      zlm(31)=zlm(29)
      zlm(32)=7.0d0*zlm(5)
      zlm(33)=-zlm(5)
      zlm(34)=a1s2*zlm(32)
      zlm(35)=-zlm(34)
      zlm(36)=a1s2*zlm(33)
      zlm(37)=-zlm(36)
      zlm(38)=a3*zlm(21)
      zlm(39)=zlm(24)
      zlm(40)=zlm(21)
      zlm(41)=-zlm(38)
      zlm(42)=4.0d0*(sqrt(315.0d0)/8.0d0)
      zlm(43)=-zlm(42)
      zlm(44)=(sqrt(315.0d0)/8.0d0)
      zlm(45)=-6.0d0*zlm(44)
      zlm(46)=zlm(44)
c     ----- l=5 -----
      zlm(47)=63.0d0*(sqrt(11.0d0)/8.0d0)
      zlm(48)=-70.0d0*(sqrt(11.0d0)/8.0d0)
      zlm(49)=15.0d0*(sqrt(11.0d0)/8.0d0)
      zlm(50)=21.0d0*(sqrt(165.0d0)/8.0d0)
      zlm(51)=-14.0d0*(sqrt(165.0d0)/8.0d0)
      zlm(52)=(sqrt(165.0d0)/8.0d0)
      zlm(53)=zlm(50)
      zlm(54)=zlm(51)
      zlm(55)=zlm(52)
      zlm(56)=6.0d0*(sqrt(1155.0d0)/4.0d0)
      zlm(57)=-a2*(sqrt(1155.0d0)/4.0d0)
      zlm(58)=a3*(sqrt(1155.0d0)/4.0d0)
      zlm(59)=-zlm(58)
      zlm(60)=-(sqrt(1155.0d0)/4.0d0)
      zlm(61)=(sqrt(1155.0d0)/4.0d0)
      zlm(62)=27.0d0*(sqrt(385.0d0/128.0d0))
      zlm(63)=-9.0d0*(sqrt(385.0d0/128.0d0))
      zlm(64)=-a3*(sqrt(385.0d0/128.0d0))
      zlm(65)=(sqrt(385.0d0/128.0d0))
      zlm(66)=-zlm(63)
      zlm(67)=-zlm(62)
      zlm(68)=-zlm(65)
      zlm(69)=-zlm(64)
      zlm(70)=4.0d0*(sqrt(3465.0d0)/8.0d0)
      zlm(71)=-zlm(70)
      zlm(72)=(sqrt(3465.0d0)/8.0d0)
      zlm(73)=-6.0d0*zlm(72)
      zlm(74)=zlm(72)
      zlm(75)=5.0d0*(sqrt(693.0d0/128.0d0))
      zlm(76)=-10.0d0*(sqrt(693.0d0/128.0d0))
      zlm(77)=(sqrt(693.0d0/128.0d0))
      zlm(78)=zlm(77)
      zlm(79)=zlm(76)
      zlm(80)=zlm(75)
c     ----- l=6 -----
      zlm(81)=231.0d0*(sqrt(13.0d0)/16.0d0)
      zlm(82)=-315.0d0*(sqrt(13.0d0)/16.0d0)
      zlm(83)=105.0d0*(sqrt(13.0d0)/16.0d0)
      zlm(84)=-5.0d0*(sqrt(13.0d0)/16.0d0)
      zlm(85)=33.0d0*(sqrt(273.0d0)/8.0d0)
      zlm(86)=-30.0d0*(sqrt(273.0d0)/8.0d0)
      zlm(87)=5.0d0*(sqrt(273.0d0)/8.0d0)
      zlm(88)=zlm(85)
      zlm(89)=zlm(86)
      zlm(90)=zlm(87)
      zlm(91)=66.0d0*(sqrt(1365.0d0/512.0d0))
      zlm(92)=-36.0d0*(sqrt(1365.0d0/512.0d0))
      zlm(93)=a2*(sqrt(1365.0d0/512.0d0))
      zlm(94)=33.0d0*(sqrt(1365.0d0/512.0d0))
      zlm(95)=-zlm(94)
      zlm(96)=-18.0d0*(sqrt(1365.0d0/512.0d0))
      zlm(97)=-zlm(96)
      zlm(98)=(sqrt(1365.0d0/512.0d0))
      zlm(99)=-zlm(98)
      zlm(100)=zlm(91)
      zlm(101)=-22.0d0*zlm(98)
      zlm(102)=zlm(96)
      zlm(103)=6.0d0*zlm(98)
      zlm(104)=-zlm(101)
      zlm(105)=-zlm(91)
      zlm(106)=-zlm(103)
      zlm(107)=zlm(97)
      zlm(108)=44.0d0*(sqrt(819.0d0)/16.0d0)
      zlm(109)=-zlm(108)
      zlm(110)=-4.0d0*(sqrt(819.0d0)/16.0d0)
      zlm(111)=-zlm(110)
      zlm(112)=11.0d0*(sqrt(819.0d0)/16.0d0)
      zlm(113)=-66.0d0*(sqrt(819.0d0)/16.0d0)
      zlm(114)=zlm(112)
      zlm(115)=-(sqrt(819.0d0)/16.0d0)
      zlm(116)=6.0d0*(sqrt(819.0d0)/16.0d0)
      zlm(117)=zlm(115)
      zlm(118)=5.0d0*(sqrt(9009.0d0/128.0d0))
      zlm(119)=-10.0d0*(sqrt(9009.0d0/128.0d0))
      zlm(120)=(sqrt(9009.0d0/128.0d0))
      zlm(121)=zlm(120)
      zlm(122)=zlm(119)
      zlm(123)=zlm(118)
      zlm(124)=6.0d0*(sqrt(3003.0d0/512.0d0))
      zlm(125)=-20.0d0*(sqrt(3003.0d0/512.0d0))
      zlm(126)=zlm(124)
      zlm(127)=(sqrt(3003.0d0/512.0d0))
      zlm(128)=-15.0d0*zlm(127)
      zlm(129)=-zlm(128)
      zlm(130)=-zlm(127)
c     ----- l=7 -----
c   ** 0
      zlm(131)=-35.d0*sqrt(15.d0)/16.d0
      zlm(132)=-9.d0*zlm(131)
      zlm(133)=-693.d0*sqrt(15.d0)/16.d0
      zlm(134)=429.d0*sqrt(15.d0)/16.d0
c  ** -1
      zlm(135)=-35.d0*sqrt(15.d0/7.d0)/32.d0
      zlm(136)=-27.d0*zlm(135)
      zlm(137)=99.d0*zlm(135)
      zlm(138)=3003.d0*sqrt(15.d0/7.d0)/32.d0
c  ** 1
      zlm(139)=zlm(135)
      zlm(140)=zlm(136)
      zlm(141)=zlm(137)
      zlm(142)=zlm(138)
c  ** -2
      zlm(143)=315.d0*sqrt(5.d0/14.d0)/8.d0
      zlm(144)=-2310.d0*sqrt(5.d0/14.d0)/8.d0
      zlm(145)=3003.d0*sqrt(5.d0/14.d0)/8.d0
c  **  2
      zlm(146)=zlm(143)/2.d0
      zlm(147)=zlm(144)/2.d0
      zlm(148)=zlm(145)/2.d0
      zlm(149)=-zlm(143)/2.d0
      zlm(150)=-zlm(144)/2.d0
      zlm(151)=-zlm(145)/2.d0
c  ** -3
      zlm(152)=189.d0*sqrt(5.d0/7.d0)/32.d0
      zlm(153)=-252.d0*sqrt(5.d0/7.d0)/32.d0
      zlm(154)=-23.d0*zlm(152)
      zlm(155)=13167.d0*sqrt(5.d0/7.d0)/32.d0
      zlm(156)=-9009.d0*sqrt(5.d0/7.d0)/32.d0
      zlm(157)=-22.d0*zlm(153)
      zlm(158)=-12012.d0*sqrt(5.d0/7.d0)/32.d0
c  **  3
      zlm(159)=-zlm(152)
      zlm(160)=-zlm(153)
      zlm(161)=-zlm(154)
      zlm(162)=-zlm(155)
      zlm(163)=-zlm(156)
      zlm(164)=-zlm(157)
      zlm(165)=-zlm(158)
c  ** -4
      zlm(166)=-693.d0*sqrt(5.d0/77.d0)/4.d0
      zlm(167)=-2.d0*zlm(166)
      zlm(168)=924.d0*sqrt(5.d0/77.d0)
      zlm(169)=-3003.d0*sqrt(5.d0/77.d0)/2.d0
      zlm(170)=zlm(169)/2.d0
c   ** 4
      zlm(171)=693.d0*sqrt(5.d0/77.d0)/2.d0
      zlm(172)=-1848.d0*sqrt(5.d0/77.d0)
      zlm(173)=3003.d0*sqrt(5.d0/77.d0)/2.d0
      zlm(174)=-zlm(171)
      zlm(175)=zlm(173)
      zlm(176)=-693.d0*sqrt(5.d0/77.d0)/16.d0
      zlm(177)=4389.d0*sqrt(5.d0/77.d0)/16.d0
      zlm(178)=-6699.d0*sqrt(5.d0/77.d0)/16.d0
      zlm(179)=3003.d0*sqrt(5.d0/77.d0)/16.d0
c  **  -5
      zlm(180)=-1155.d0*sqrt(5.d0/77.d0)/32.d0
      zlm(181)=-4.d0*zlm(180)
      zlm(182)=-3696.d0*sqrt(5.d0/77.d0)/32.d0
      zlm(183)=-15.d0*zlm(180)
      zlm(184)=27.d0*zlm(180)
      zlm(185)=-13.d0*zlm(180)
      zlm(186)=56.d0*zlm(180)
      zlm(187)=-52.d0*zlm(180)
      zlm(188)=-13.d0*zlm(182)
c  **  5
      zlm(189)=zlm(180)
      zlm(190)=zlm(181)
      zlm(191)=zlm(182)
      zlm(192)=zlm(183)
      zlm(193)=zlm(184)
      zlm(194)=zlm(185)
      zlm(195)=zlm(186)
      zlm(196)=zlm(187)
      zlm(197)=zlm(188)
c  **  -6
      zlm(198)=9009.d0*sqrt(5.d0/2002.d0)/8.d0
      zlm(199)=-2.d0*zlm(198)
      zlm(200)=zlm(198)
      zlm(201)=-6006.d0*sqrt(5.d0/2002.d0)
      zlm(202)=-zlm(201)
      zlm(203)=-zlm(201)
c  **  6
      zlm(204)=27027.d0*sqrt(5.d0/2002.d0)/8.d0
      zlm(205)=-2.d0*zlm(204)
      zlm(206)=zlm(204)
      zlm(207)=-9009.d0*sqrt(5.d0/2002.d0)
      zlm(208)=-zlm(207)
      zlm(209)=6006.d0*sqrt(5.d0/2002.d0)
      zlm(210)=-3003.d0*sqrt(5.d0/2002.d0)/16.d0
      zlm(211)=-3.d0*zlm(210)
      zlm(212)=3.d0*zlm(210)
      zlm(213)=-zlm(210)
c  ** -7
      zlm(214)=3003.d0*sqrt(5.d0/143.d0)/32.d0
      zlm(215)=-8.d0*zlm(214)
      zlm(216)=16.d0*zlm(214)
      zlm(217)=-858.d0*sqrt(5.d0/143.d0)
      zlm(218)=-3.d0*zlm(214)
      zlm(219)=3.d0*zlm(214)
      zlm(220)=-zlm(214)
      zlm(221)=zlm(216)
      zlm(222)=zlm(215)
      zlm(223)=-zlm(216)
c  **  7
      zlm(224)=-zlm(214)
      zlm(225)=-zlm(215)
      zlm(226)=-zlm(216)
      zlm(227)=-zlm(217)
      zlm(228)=-zlm(218)
      zlm(229)=-zlm(219)
      zlm(230)=-zlm(220)
      zlm(231)=-zlm(221)
      zlm(232)=-zlm(222)
      zlm(233)=-zlm(223)
c     ----- l=8 -----
c  **  0
      zlm(234)=35.d0*sqrt(17.d0)/128.d0
      zlm(235)=-36.d0*zlm(234)
      zlm(236)=198.d0*zlm(234)
      zlm(237)=-3003.d0*sqrt(17.d0)/32.d0
      zlm(238)=6435.d0*sqrt(17.d0)/128.d0
c  ** -1
      zlm(239)=-105.d0*sqrt(17.d0)/32.d0
      zlm(240)=-11.d0*zlm(239)
      zlm(241)=-3003.d0*sqrt(17.d0)/32.d0
      zlm(242)=2145.d0*sqrt(17.d0)/32.d0
c  **  1
      zlm(243)=zlm(239)
      zlm(244)=zlm(240)
      zlm(245)=zlm(241)
      zlm(246)=zlm(242)
c  ** -2
      zlm(247)=-105.d0*sqrt(17.d0/70.d0)/16.d0
      zlm(248)=-33.d0*zlm(247)
      zlm(249)=143.d0*zlm(247)
      zlm(250)=-zlm(249)
c  **  2
      zlm(251)=-105.d0*sqrt(17.d0/70.d0)/32.d0
      zlm(252)=-33.d0*zlm(251)
      zlm(253)=143.d0*zlm(251)
      zlm(254)=-zlm(253)
      zlm(255)=-zlm(251)
      zlm(256)=-zlm(252)
      zlm(257)=-zlm(253)
      zlm(258)=-zlm(254)
c  ** -3
      zlm(259)=10395.d0*sqrt(17.d0/1155.d0)/32.d0
      zlm(260)=-100485.d0*sqrt(17.d0/1155.d0)/32.d0
      zlm(261)=225225.d0*sqrt(17.d0/1155.d0)/32.d0
      zlm(262)=-13.d0*zlm(259)
      zlm(263)=-3465.d0*sqrt(17.d0/1155.d0)/8.d0
      zlm(264)=15015.d0*sqrt(17.d0/1155.d0)/4.d0
      zlm(265)=13.d0*zlm(263)
c  **  3
      zlm(266)=-zlm(259)
      zlm(267)=-zlm(260)
      zlm(268)=-zlm(261)
      zlm(269)=-zlm(262)
      zlm(270)=-zlm(263)
      zlm(271)=-zlm(264)
      zlm(272)=-zlm(265)
c  ** -4
      zlm(273)=231.d0*sqrt(17.d0/77.d0)/16.d0
      zlm(274)=-27.d0*zlm(273)
      zlm(275)=91.d0*zlm(273)
      zlm(276)=-65.d0*zlm(273)
      zlm(277)=-2.d0*zlm(273)
      zlm(278)=52.d0*zlm(273)
      zlm(279)=-130.d0*zlm(273)
c  **  4
      zlm(280)=-231.d0*sqrt(17.d0/77.d0)/8.d0
      zlm(281)=-27.d0*zlm(280)
      zlm(282)=91.d0*zlm(280)
      zlm(283)=-65.d0*zlm(280)
      zlm(284)=-zlm(280)
      zlm(285)=26.d0*zlm(280)
      zlm(286)=-65.d0*zlm(280)
      zlm(287)=231.d0*sqrt(17.d0/77.d0)/64.d0
      zlm(288)=-28.d0*zlm(287)
      zlm(289)=118.d0*zlm(287)
      zlm(290)=-156.d0*zlm(287)
      zlm(291)=65.d0*zlm(287)
c  ** -5
      zlm(292)=-15015.d0*sqrt(17.d0/1001.d0)/32.d0
      zlm(293)=-7.d0*zlm(292)
      zlm(294)=11.d0*zlm(292)
      zlm(295)=-5.d0*zlm(292)
      zlm(296)=-4.d0*zlm(292)
      zlm(297)=24.d0*zlm(292)
      zlm(298)=-20.d0*zlm(292)
      zlm(299)=-3003.d0*sqrt(17.d0/1001.d0)/2.d0
      zlm(300)=-5.d0*zlm(299)
c  **  5
      zlm(301)=zlm(292)
      zlm(302)=zlm(293)
      zlm(303)=zlm(294)
      zlm(304)=zlm(295)
      zlm(305)=zlm(296)
      zlm(306)=zlm(297)
      zlm(307)=zlm(298)
      zlm(308)=zlm(299)
      zlm(309)=zlm(300)
c  ** -6
      zlm(310)=-1287.d0*sqrt(17.d0/858.d0)/16.d0
      zlm(311)=429.d0*sqrt(17.d0/858.d0)
      zlm(312)=-zlm(311)
      zlm(313)=-17.d0*zlm(310)
      zlm(314)=31.d0*zlm(310)
      zlm(315)=-15.d0*zlm(310)
      zlm(316)=-16.d0*zlm(311)
      zlm(317)=15.d0*zlm(311)
      zlm(318)=zlm(317)
c  **  6
      zlm(319)=429.d0*sqrt(17.d0/858.d0)/32.d0
      zlm(320)=-18.d0*zlm(319)
      zlm(321)=48.d0*zlm(319)
      zlm(322)=-32.d0*zlm(319)
      zlm(323)=306.d0*zlm(319)
      zlm(324)=-558.d0*zlm(319)
      zlm(325)=270.d0*zlm(319)
      zlm(326)=-768.d0*zlm(319)
      zlm(327)=720.d0*zlm(319)
      zlm(328)=480.d0*zlm(319)
      zlm(329)=-18.d0*zlm(319)
      zlm(330)=48.d0*zlm(319)
      zlm(331)=-46.d0*zlm(319)
      zlm(332)=15.d0*zlm(319)
c  ** -7
      zlm(333)=15015.d0*sqrt(17.d0/715.d0)/32.d0
      zlm(334)=-3.d0*zlm(333)
      zlm(335)=3.d0*zlm(333)
      zlm(336)=-zlm(333)
      zlm(337)=-8.d0*zlm(333)
      zlm(338)=16.d0*zlm(333)
      zlm(339)=-8.d0*zlm(333)
      zlm(340)=16.d0*zlm(333)
      zlm(341)=-16.d0*zlm(333)
      zlm(342)=-4290.d0*sqrt(17.d0/715.d0)
c  **  7
      zlm(343)=-zlm(333)
      zlm(344)=-zlm(334)
      zlm(345)=-zlm(335)
      zlm(346)=-zlm(336)
      zlm(347)=-zlm(337)
      zlm(348)=-zlm(338)
      zlm(349)=-zlm(339)
      zlm(350)=-zlm(340)
      zlm(351)=-zlm(341)
      zlm(352)=-zlm(342)
c  ** -8
      zlm(353)=2145.d0*sqrt(17.d0/715.d0)/16.d0
      zlm(354)=-10.d0*zlm(353)
      zlm(355)=24.d0*zlm(353)
      zlm(356)=-16.d0*zlm(353)
      zlm(357)=-3.d0*zlm(353)
      zlm(358)=3.d0*zlm(353)
      zlm(359)=-zlm(353)
      zlm(360)=20.d0*zlm(353)
      zlm(361)=-10.d0*zlm(353)
      zlm(362)=-24.d0*zlm(353)
c  **  8
      zlm(363)=2145.d0*sqrt(17.d0/715.d0)/128.d0
      zlm(364)=-2145.d0*sqrt(17.d0/715.d0)/4.d0
      zlm(365)=-5.d0*zlm(364)
      zlm(366)=8.d0*zlm(364)
      zlm(367)=-4.d0*zlm(364)
      zlm(368)=-3.d0*zlm(364)
      zlm(369)=3.d0*zlm(364)
      zlm(370)=-zlm(364)
      zlm(371)=10.d0*zlm(364)
      zlm(372)=-5.d0*zlm(364)
      zlm(373)=-8.d0*zlm(364)
      zlm(374)=-4.d0*zlm(363)
      zlm(375)=6.d0*zlm(363)
      zlm(376)=-4.d0*zlm(363)
      zlm(377)=zlm(363)
      return
      end
