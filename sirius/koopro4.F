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
      subroutine koopro4(WRK,LWRK)
C-----------------------------------------------------------------------
C
C     This program calculates the second order correction to the energy
C     using the n-electron valence state perturbation thoery (NEVPT)
C     in the partially contracted (PC) ad strongly contracted (SC)
C     variants using the Dyall's Hamiltonian for the definition of the
C     zero order energies. All the equations are based on the spin free
C     formulation of NEVPT.
C
C     Relevant papers:
C     1) C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J-P Malrieu,
C        J. Chem. Phys., 114, 10252, (2001)
C        Introduction of NEVPT
C     2) C. Angeli, R. Cimiraglia, J-P Malrieu,
C        Chem. Phys. Lett., 350, 297-305, (2001)
C        Introduction of the formalism based on the "Koopmans' operators"
C     3) C. Angeli, R. Cimiraglia, J-P Malrieu,
C        J. Chem. Phys., 117, 10525-10533, (2002)
C        Spin free version of NEVPT
C
C
C-----------------------------------------------------------------------
#include <single.h>
#include <implicit.h>
#include <maxorb.h>
#include <mxpdim.h>
#include <priunit.h>
#include <inforb.h>
#include <infopt.h>
#include <infinp.h>
#include <infvar.h>
#include <infdim.h>
#include <inftap.h>
#include <infpt2.h>
#include <strnum.h>
#include <dummy.h>
#include <csfbas.h>

! ** The nevpt2 code was not compiled unless PGF90, GFORTRAN or IFORT compiler.
! ** Now code is compiled for all compilers, if your compiler cannot handle
! ** the code then do NOT #define GO_ON_COMP in dalcip.F and koopro4.F (below) /Jan. 2011

#define GO_ON_COMP

#ifdef GO_ON_COMP
      REAL*8  koopaa,koopbb,koopab,koopeaa,koopebb,koopeab
      logical*1 zseg
      COMMON /CIP/NORBFE,NOCB,NCF,ZSEG
      common /actspace/mact,nact2,nact3,nact4
#include "nevpt2.h"
c      real*4 tarray(2)
      logical zcapos,zion
      integer*1 ioccfe
      integer*2 icomp,iconf,isegno
      allocatable coopipa(:,:),coopeaa(:,:)
      allocatable dal(:,:),koopaa(:,:,:,:),koopeaa(:,:,:,:)
      allocatable koopbb(:,:,:,:),koopab(:,:)
      allocatable daaa(:,:,:,:,:,:)
      allocatable taa(:,:,:,:)
      allocatable amat(:,:,:,:,:,:),bmat(:,:,:,:),cmat(:,:,:,:),dmat(:,:
     $     )
      allocatable atmat(:,:,:,:,:,:),btmat(:,:,:,:),ctmat(:,:,:,:)
     $     ,dtmat(:,:)
      allocatable ro4(:,:,:,:,:)
      allocatable atwo(:,:,:,:)
      allocatable ne(:),nd(:),ispinfe(:),iorbfe(:),itsym(:)
      allocatable cvec(:)
      allocatable trou(:),part(:)
      allocatable f(:),aj(:),ak(:),num(:),indic(:),jndic(:),lndic(:)
      allocatable ioccfe(:,:)
      allocatable icomp(:),iconf(:),isegno(:),zcapos(:)

      INTEGER*2 ispinfe,IORBFE
      INTEGER*2 NE,TROU,PART
      logical *1 HFEX,ZGEL
      dimension WRK(LWRK),its(8,8)
      character*8 RTNLBL(2)
      integer*1 flip
      allocatable istr(:),jstr(:),onel(:),zgel(:)
      allocatable in2d(:),id2n(:)
      allocatable detvec(:)
      allocatable iocca(:),ioccb(:),flip(:)
C--- Cele  Debug
c      character*1 nexc(13)
c      allocatable sppm(:),orbfe(:)
c      integer orbfe
c      character*1 sppm
C--- Cele  Debug
c     indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      CALL QENTER('koopro4')
      KFRSAV=1
      KFREE=KFRSAV
      LFREE=LWRK
C--- Cele  Debug
c       nexc(1)='F'
c       nexc(2)='M'
c       nexc(3)='D'
c       nexc(4)='T'
c       nexc(5)='Q'
c       nexc(6)='P'
c       nexc(7)='H'
c       nexc(8)='E'
c       nexc(9)='O'
c       nexc(10)='N'
c       nexc(11)='X'
c       nexc(12)='U'
c       nexc(13)='Z'
c       allocate(sppm(2*norbt+1))
c       allocate(orbfe(2*norbt+1))
c       do i=1,norbt
c       sppm(i)='-'
c       orbfe(i)=i
c       sppm(i+norbt)='+'
c       orbfe(i+norbt)=i  
c       enddo
c       sppm(norbt+norbt+1)='+'
c       orbfe(norbt+norbt+1)=norbt+1
C--- Cele  Debug

      if(DOCI.AND.ISTNEVCI.eq.0)then
         write(lupri,*)' '
         write(lupri,*)'<<<<<<<<<<<===========>>>>>>>>>>>>'
         write(lupri,*)'ATTENTION ATTENTION ATTENTION ATTENTION'
         write(lupri,'(a,$)'
     $        )'Number of state is required with .CI keyword'
         write(lupri,*)
         write(lupri,*)'Please restart'
         CALL QUIT('NEVPT2 error: .CI but .STATE not defined')
      endif
      nwall=0
      norbfe=norbt
      t=second()
      write (lupri,'(///A)')' ***************************************'//
     * '*******************************'
      IF (DoMCSRDFT) THEN
         write (lupri,'(A)')' *              SC and PC NEVPT2-srDFT '//
     *   'Calculation                    *'
      ELSE
         write (lupri,'(A)')' *                    SC and PC NEVPT2 '//
     *   'Calculation                    *'
      END IF
      write (lupri,'(A///)')' ***************************************'//
     * '*******************************'
      call flshfo(lupri)

       do i=1,8
        do j=1,8
        its(i,j)=muld2h(i,j)
       enddo
       enddo

C
C     First of all we prepare a vector (in2d) which converts the OM indices
C     from the NEVPT order to the DALTON order 
C     
      call nwadd(nwall,2*norbt*4+norbt)
      allocate (in2d(norbt))
      allocate (id2n(norbt))
      allocate (zgel(norbt))

C
C     Getting the coefficients of the determinants
C
      call nwadd(nwall,ncdets*8)
      allocate (cvec(ncdets))
      rewind (luit1)
      CALL MOLLB2 ('STARTVEC',RTNLBL,LUIT1,LUERR)
      if (.not.doci) then
         ioldstate=istate
         istate=1
      endif
      IF (ISTATE .GT. 0) THEN
C            ... CI with ISTACI .gt. 0 or MC calculation
            DO 102 I = 1,(ISTNEVCI-1)
               READ (LUIT1)
 102        CONTINUE
            ILOW = I
            IHGH = I
      ELSE
C        else MCSCF with only reference vector saved (NRSAVE)
            ILOW = ABS(I)
            IHGH = ABS(I)
      END IF
      if (.not.doci) then
         istate=ioldstate
      endif
      if(.NOT.FLAG(27)) then
C
C     Getting the informations on the determinants in the MR wavefunction
C
         CALL MEMGET('REAL',KCINDX,LCINDX,WRK,KFREE,LFREE)
         CALL MEMGET('REAL',KDETVEC,NCDETS,WRK,KFREE,LFREE)
         CALL GETCIX(WRK(KCINDX),LSYM,LSYM,WRK(KFREE),LFREE,0)
         CALL READT(LUIT1,NCONF,CVEC)
         call csdtvc(cvec,WRK(KDETVEC),1,WRK(KCINDX+KDTOC-1)
     $        ,WRK(KCINDX+KICTS(1)-1),lsym,1,0)
         CALL MEMREL('after csdtvc',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      else
         CALL READT(LUIT1,NCDETS,CVEC)
      endif
C     CALL WRTMAT(CVEC,1,NCDETS,1,NCDETS,0)

C
C     Getting the one-particle density matrix
C     (Not used for the moment!)
C
C Manu April 2009: now it is used in the NEVPT2-srDFT :-)

      
      do i=1,norbt
       zgel(i)=.false. ! is this orbital frozen ?
      enddo
C
C     in2d(jorb_nevpt2) = jorb_dalton : from nevpt2 to dalton orbital order
C     nevpt2: loop type, symmetry; dalton: loop symmetry, type
C
      indom=0
C     Core
      do is=1,nsym
       do i=1,nish(is)
       indom=indom+1
       in2d(indom)=iorb(is)+i
       if (i.le.NFRNEVPT2(is)) zgel(indom)=.true.
       enddo
      enddo
C     Active
      do is=1,nsym
       do i=1,nash(is)
       indom=indom+1
       in2d(indom)=iorb(is)+nish(is)+i
       enddo
      enddo
C     Virtual
      do is=1,nsym
       do i=1,norb(is)-nish(is)-nash(is)
       indom=indom+1
       in2d(indom)=iorb(is)+nish(is)+nash(is)+i
       enddo
      enddo

C
C     The active orbitals are ordered following their occupation number
C
C--- Cele Next Work
C
C      Questo pezzo e' da rivedere: l'ordinamento viene fatto bene,
C      ma con il nuovo ordinamento la sub rofour non fornisce la buona 
C      matrice densita' : approfondire.
C
C      do ia=1,nasht
C       do ja=ia+1,nasht
C       if (dal(ia,ia).lt.dal(ja,ja)) then
C        dum=dal(ja,ja)
C        dal(ja,ja)=dal(ia,ia)
C        dal(ia,ia)=dum
C        idum=in2d(nisht+ja)
C        in2d(nisht+ja)=in2d(nisht+ia)
C        in2d(nisht+ia)=idum
C        idum=in2da(ia)
C        in2da(ia)=in2da(ja)
C        in2da(ja)=idum
C       endif
C       enddo
C      enddo
C--- Cele Next Work
C
C     Defining the inverse of in2d (id2n) and in2da (id2na)
C
      do i=1,norbt
       id2n(in2d(i))=i
      enddo


      thr=thrnevpt
      NCF=NCDETS

      lfermi=nactel/2 
      zion=mod(nactel,2).eq.1
      if(zion)lfermi=(nactel+1)/2

C
C     Compute the number of holes/particles
C
      ntothp=0
      LUN98=-1000
      CALL GPOPEN(LUN98,'DETSTR',' ',' ','UNFORMATTED',IDUMMY,.FALSE.)
      LUN31=-13300
      CALL GPOPEN(LUN31,'SCR31','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *           .FALSE.)
      call nwadd(nwall,2*nael*4)
      allocate (istr(NAEL))
      allocate (jstr(NBEL))
      do id=1,ncf
         read(LUN98) dum,(istr(i),i=1,nael),(jstr(i),i=1,nbel)
         do ia=1,nael
c            istr(ia)=istr(ia)
            if (istr(ia).gt.lfermi) ntothp=ntothp+1 
         enddo
         do ib=1,nbel
c            jstr(ib)=jstr(ib)
            if (jstr(ib).gt.lfermi) ntothp=ntothp+1 
         enddo
C--- Cele Next Work
C       write (lupri,*) 
C       write (lupri,*) 'Determ umber',id
C       write (lupri,*) 
C       write (lupri,*) 'Original'
C       write (LUPRI,'(20I3)') (istr(i),i=1,nael)
C       write (LUPRI,'(20I3)') (jstr(i),i=1,nbel)
C       write (lupri,*) 'New order'
C       write (LUPRI,'(20I3)') (id2na(istr(i)),i=1,nael)
C       write (LUPRI,'(20I3)') (id2na(jstr(i)),i=1,nbel)
C       write (lupri,*) 'ntothp',ntothp
C--- Cele Next Work
      enddo ! do id=1,ncf
      if(zion)ntothp=ntothp+ncf


      call nwdel(nwall,2*nael*4)
      deallocate (istr)
      deallocate (jstr)
      call nwadd(nwall,(ncf+1000)*(2+4+2+2+2+4))
      allocate(ne(ncf+1000))
      allocate(nd(ncf+1000))
      allocate(icomp(ncf+1000))
      allocate(iconf(ncf+1000))
      allocate(isegno(ncf+1000))
      allocate(zcapos(ncf+1000))
      call nwadd(nwall,2*(ntothp+1000)*2)
      allocate(trou(ntothp+1000))
      allocate(part(ntothp+1000))
c      call nwadd(nwall,2*(nbel*ncf)*4)
c      allocate (istr(nael*ncf))
c      allocate (jstr(nbel*ncf))
      call dalcip(nael,nbel,ncf,cvec,ne,trou,part,nd,norbt,
     *            nish,nash,nasht,lfermi,nisht,nsym,LUN98,ispin)
      tdal=second()
      write(lupri,'(a,f8.2,a)')
     &   ' Elapsed time after dalcip ',tdal-t,' sec.'
      call flshfo(lupri)
      CALL GPCLOSE(LUN98,'KEEP')
      call nwdel(nwall,2*(nbel*ncf)*4)
c      deallocate (istr)
c      deallocate (jstr)

C--- Cele debug
C
C     Debug printing of determinants in the hole/part formalism
C
c      do id=1,ncf
c         write (LUPRI,111) cvec(id),nexc(ne(id)+1),(orbfe(trou(nd(id)+i)
c     $        ),sppm(trou(nd(id)+i)),i=1,ne(id)),(orbfe(part(nd(id)+i))
c     $        ,sppm(part(nd(id)+i)),i=1,ne(id))
c      enddo
c 111  format (F12.8,A,1X,30(I3,a))
c      deallocate (sppm)
c      deallocate (orbfe)
C--- Cele debug

      iditrpa=ntothp
C
C     Getting itsym
C
      call nwadd(nwall,(norbt+norbt+1)*4)
      allocate(itsym(norbt+norbt+1))
      io=0
      do is=1,nsym
       do i=1,norb(is)
        io=io+1
        itsym(id2n(io))=is
        itsym(id2n(io)+norbt)=is
       enddo
      enddo

      call nwadd(nwall,2*(norbt+norbt+1)*2)
      allocate(ispinfe(norbt+norbt+1))
      allocate(iorbfe(norbt+norbt+1))
      call nwadd(nwall,2*nnorbx*8)
      allocate(f(nnorbx))
      allocate(aj(nnorbx))
      NOCB=nisht+lfermi

      if (DOCI) then
        REWIND LUIT1
        CALL MOLLAB('ENERGIES',LUIT1,LUERR)
        READ (LUIT1) (escf,I=1,ISTNEVCI)
      else
        escf=emcscf
      endif

      write (lupri,'(/1x,a,i8)')      'NORB   =',norbt
      write (lupri,'(1x,a,i8)')       'NOCB   =',NOCB
      write (lupri,'(1x,a,i8)')       'NCF    =',ncf
      write (lupri,'(1x,a,i8)')       'NTOTHP =',ntothp
      write (lupri,'(1x,a,8I4)')      'FROZEN =',(NFRNEVPT2(i),i=1,nsym)
      if(doci)then
         write (lupri,'(1x,a,i8)')    'STATE  =',ISTNEVCI
         write (lupri,'(1x,a,f16.8)') 'E-CASCI=',escf
      else
         write (lupri,'(1x,a,f16.8)') 'E-MCSCF=',escf
      endif
      write (lupri,'(1x,a,1p,d16.2)') 'THR    =',thr
      call flshfo(lupri)
      mact=nasht
      nact2=nasht**2
      nact3=nasht**3
      nact4=nasht**4
      nele=nactel           
c      if(zion)nele=nele-1
      ndim2=nasht*(nasht+1)/2
      ndim3=ndim2*(ndim2+1)/2
      if(nasht.gt.max_act)then
         write(lupri,*)
     &     'no. of active orbitals greater than',max_act,
     &     ': not supported'
         CALL QUIT('nevpt2 error: nact greater than "max_act" ')
      endif
      do i=1,nasht
        lindi(i)=(i-1)*i*(i+1)*(i+2)/24
        lindj(i)=(i-1)*i*(i+1)/6
        lindk(i)=i*(i-1)/2
      enddo
      do l=1,nasht
       do k=1,nasht
        do j=1,nasht
         do i=1,nasht
          call getnorm(i,j,k,l,norm)
          normind(i,j,k,l)=norm
         enddo
        enddo
       enddo
      enddo
      nwords=nasht*(nasht+1)*(nasht+2)*(nasht+3)/24
      call nwadd(nwall,nasht**6*8)
      allocate(daaa(1:nasht,1:nasht,1:nasht,1:nasht,1:nasht,1:nasht))
      call nwadd(nwall,nasht**2*8)
      allocate(dal(1:nasht,1:nasht))
      call nwadd(nwall,nasht**4*8)
      allocate(taa(1:nasht,1:nasht,1:nasht,1:nasht))
      call nwadd(nwall,nasht**6*8)
      allocate(amat(1:nasht,1:nasht,1:nasht,1:nasht,1:nasht,1:nasht))
      call nwadd(nwall,nasht**4*8)
      allocate(atwo(1:nasht,1:nasht,1:nasht,1:nasht))
      call nwadd(nwall,nwords*nasht**4*8)
      allocate(ro4(1:nwords,1:nasht,1:nasht,1:nasht,1:nasht))
      do  i=1,norbt
         ispinfe(i)=0
         ispinfe(i+norbt)=1
         iorbfe(i)=i
         iorbfe(i+norbt)=i
      enddo
      if(zion)then
         iorbfe(norbt+norbt+1)=norbt+1
         ispinfe(norbt+norbt+1)=1
      endif
C
C Manu comment: the one-electron part of the Hamiltonian is built here !
C
      call nwadd(nwall,NNORBX*8)
      ALLOCATE (onel(NNORBX))
C      
C ***********************************************************
C Manu b
      IF (DOMCSRDFT) THEN
C        Build active density matrix wrk(dv) with ro1 ... 

         CALL MEMGET('REAL',KDV,NNASHX,WRK,KFREE,LFREE)
         call dzero(dal,n2ashx)
         call ro1(cvec,dal,nisht,nasht,thr,trou
     $        ,part,ne,nd,ispinfe
     $        ,iorbfe,itsym)
         CALL DGETSP(NASHT,dal,wrk(kdv)) 
      ELSE
         CALL MEMGET('REAL',KDV,0,WRK,KFREE,LFREE)
      END IF
      call readint(onel,atwo,aj,wrk(kdv),
     &     wrk(kfree),lfree,id2n,in2d,lupri)
C Manu e      

      CALL MEMREL('DV',WRK,KFRSAV,KFRSAV,KFREE,LFREE)
      
!     write (lupri,'(/A)') 'Done with readint'
!     call flshfo(lupri)

      call dzero(ro4,nwords*nasht**4)
      call dzero(daaa,nasht**6)
      call dzero(taa,nasht**4)
      call dzero(dal,nasht**2)

      write (lupri,'(/A,I5)') ' Spin multiplicity is',ispin
c      if(ispin.eq.1.or.ispin.eq.2.or.ispin.eq.3)then
      allocate(ioccfe(norbt,1))
      call esclass(nd,ne,trou,part,icomp,iconf,isegno,ncf,NOCB,norbt
     $     ,ispin,ioccfe)
      deallocate(ioccfe)
c      else
c         do i=1,ncf
c            icomp(i)=1
c         enddo
c      endif
      ncapos=0
      inext=1
      do i=1,ncf
         if(i.eq.inext)then
            zcapos(i)=.true.
            inext=i+icomp(i)
         else
            zcapos(i)=.false.
         endif
         if(zcapos(i))ncapos=ncapos+1
      enddo
      write (lupri,'(a,i7,2a)') ' There are',ncapos,' different orbital'
     $     ,' occupation patterns'
      call flshfo(lupri)
      call nwadd(nwall,norbt*ncapos*2)
      allocate(ioccfe(1:norbt,1:ncapos))
      icapos=0
      do i=1,ncf
         if(zcapos(i))then
            icapos=icapos+1
            call giveocc(ioccfe(1,icapos),i,nd,ne,trou,part,NOCB,norbt)
         endif
      enddo
!     write(lupri,'(/A)) ' Density matrices calculation'
!     call flshfo(lupri)
      allocate(iocca(norbfe+1))
      allocate(ioccb(norbfe+1))
      allocate(flip(ncf+1000))
      if(nele.ge.4)then
         call rofour(cvec,ro4,nisht,nasht,thr,icomp,iconf,ioccfe,zcapos,
     $        ncapos,trou,part,ne,nd,ispinfe,iorbfe,itsym,iocca,ioccb
     $        ,flip)
         call bro3(ro4,daaa,nasht,nele)
         call bro2(daaa,taa,nasht,nele)
         call bro1(taa,dal,nasht,nele)
      elseif(nele.eq.3)then
         call ro3(cvec,daaa,nisht,nasht,thr,icomp,iconf,ioccfe,zcapos
     $        ,ncapos,trou,part,ne,nd,ispinfe,iorbfe,itsym,iocca,ioccb
     $        ,flip)
         call bro2(daaa,taa,nasht,nele)
         call bro1(taa,dal,nasht,nele)
      elseif(nele.eq.2)then
         call ro2(cvec,taa,nisht,nasht,thr,trou,part,ne,nd,
     *        ispinfe,iorbfe,itsym,iocca,ioccb,flip)
         call bro1(taa,dal,nasht,nele)
      elseif(nele.eq.1)then
         write(lupri,*) 'just one active electron: are you kidding?'
         call ro1(cvec,dal,nisht,nasht,thr,trou,part,ne,nd,
     *        ispinfe,iorbfe,itsym)
      endif
      deallocate(iocca)
      deallocate(ioccb)
      deallocate(flip)
!     write(lupri,'(/A)') ' Done with density matrices calculation'
!     call flshfo(lupri)

      call nwdel(nwall,ncdets*8)
      deallocate (cvec)
      call nwdel(nwall,(ncf+1000)*(2+4+2+2+2+1))
      deallocate (ne)
      deallocate (nd)
      deallocate (trou)
      deallocate (part)
      deallocate (icomp)
      deallocate (iconf)
      deallocate (isegno)
      deallocate (zcapos)
      call nwdel(nwall,2*(norbt+norbt+1)*2)
      deallocate (ispinfe)
      deallocate (iorbfe)
      call nwdel(nwall,norbt*ncapos*2)
      deallocate (ioccfe)
      t2=second()
      write(lupri,'(a,f10.2,a)')
     &   ' Elapsed time after density matrices',t2-tdal,' seconds'
      n=nasht
      nstop=(n+n*n)/2
c---calcolo Koopmans esteso doppie ionizzazioni---------------------
      call nwadd(nwall,2*nasht**2*8)
      allocate(coopipa(1:nasht,1:nasht))
      allocate(coopeaa(1:nasht,1:nasht))
      call nwadd(nwall,2*nasht**4*8)
      allocate(koopaa(1:nasht,1:nasht,1:nasht,1:nasht))
      allocate(koopbb(1:nasht,1:nasht,1:nasht,1:nasht))
      call nwadd(nwall,nasht**2*8)
      allocate(koopab(1:nasht,1:nasht))
      call nwadd(nwall,nasht**4*8)
      allocate(koopeaa(1:nasht,1:nasht,1:nasht,1:nasht))

      call koopE(atwo,coopipa,coopeaa,taa,dal,nasht,nisht,f,aj)
c         write (lupri,*) 'Koopipa'
c         do i=1,6
c          write (lupri,*) i,(coopipa(i,j),j=1,6)
c         enddo
c         write (lupri,*) 'Koopeaa'
c         do i=1,6
c          write (lupri,*) i,(coopeaa(i,j),j=1,6)
c         enddo
      call wr1(dal,coopeaa,coopipa,nasht,lun31)
      call nwdel(nwall,2*nasht**2*8)
      deallocate (coopipa)
      deallocate (coopeaa)
      call wr2(koopaa,taa,dal,nasht,lun31)
!     write(lupri,*)
!     write(lupri,*) 'Done with one part/hole Koopmans'' matrices'
      call koop2E(atwo,daaa,taa,dal,koopaa,koopeaa,nasht,nisht
     $     ,f,aj,itsym,its)
!     write(lupri,*)
!     write(lupri,*) 'Done with two part/hole Koopmans'' matrices'
      t3=second()
      write(lupri,'(a,f10.2,a)')
     &   ' KoopE and KoopE2 have taken        ',t3-t2,' seconds'
      call wr3(koopeaa,taa,koopaa,nasht,lun31)
      call koopman0pE(atwo,daaa,taa,dal,koopaa,koopbb,koopab,
     $ nasht,nisht,f,aj,itsym,its)
!     write(lupri,*)
!     write(lupri,*) 'Done with excitation Koopmans'' matrix'
      t4=second()
      write(lupri,'(a,f10.2,a)')
     &   ' Koopman0pE has taken               ',t4-t3,' seconds'
      fmax=0.d0
      do j=1,nasht
         do i=1,nasht
            fprov=abs(koopab(i,j))
            if(fprov.gt.fmax)fmax=fprov
         enddo
      enddo
      write(lupri,'(/a,1p,d10.3)') ' Maximum value of F matrix is:',fmax
      call wr4(koopaa,koopbb,koopab,nasht,lun31)
c-----scrittura matrici A e D-----------------------------------

      call nwdel(nwall,2*nasht**4*8)
      deallocate(koopaa)
      deallocate(koopbb)
      call nwdel(nwall,nasht**2*8)
      deallocate(koopab)
      call nwdel(nwall,nasht**4*8)
      deallocate(koopeaa)

      call nwadd(nwall,nasht**6*8)
      allocate(atmat(1:nasht,1:nasht,1:nasht,1:nasht,1:nasht,1:nasht))
      call nwadd(nwall,4*nasht**4*8)
      allocate(bmat(1:nasht,1:nasht,1:nasht,1:nasht))
      allocate(cmat(1:nasht,1:nasht,1:nasht,1:nasht))
      allocate(btmat(1:nasht,1:nasht,1:nasht,1:nasht))
      allocate(ctmat(1:nasht,1:nasht,1:nasht,1:nasht))
      call nwadd(nwall,2*nasht**2*8)
      allocate(dmat(1:nasht,1:nasht))
      allocate(dtmat(1:nasht,1:nasht))

      call dzero(amat,nasht**6)
      call dzero(atmat,nasht**6)
      call bamat(amat,atmat,atwo,ro4,daaa,taa,dal,nasht,nisht,f,aj,itsym
     $     ,its,norbt)

      call nwdel(nwall,nwords*nasht**4*8)
      deallocate (ro4)
      t6=second()
      write(lupri,'(a,f10.2,a)')' BAMAT has taken ',t6-t4,
     * ' seconds'
      call bbmat(atwo,bmat,btmat,daaa,taa,dal,nasht,nisht,f,aj)
      t7=second()
      write(lupri,'(a,f10.2,a)')' BBMAT has taken ',t7-t6,' seconds'
      call bcmat(atwo,cmat,ctmat,daaa,taa,dal,nasht,nisht,f,aj)
      t8=second()
      write(lupri,'(a,f10.2,a)')' BCMAT has taken ',t8-t7,' seconds'
      call fdiff(bmat,cmat,btmat,ctmat,nasht)
      call bdmat(atwo,dmat,dtmat,taa,dal,nasht,nisht,f,aj)
      t9=second()
      write(lupri,'(a,f10.2,a)')' BDMAT has taken ',t9-t8,' seconds'
      call wr5(daaa,amat,atmat,bmat,btmat,cmat,ctmat,dmat,dtmat,nasht
     $     ,lun31)

      write(lupri,'(/a,f8.2,a)')' Total elapsed time in Koopro4',t9-t,
     * ' seconds'
      call flshfo(lupri)

      call nwdel(nwall,nasht**6*8)
      deallocate (daaa)
      call nwdel(nwall,nasht**2*8)
      deallocate (dal)
      call nwdel(nwall,nasht**4*8)
      deallocate (taa)
      call nwdel(nwall,nasht**6*8)
      deallocate (amat)
      call nwdel(nwall,nasht**4*8)
      deallocate (atwo)
      call nwdel(nwall,(norbt+norbt+1)*4)
      deallocate (f)
      call nwdel(nwall,nasht**6*8)
      deallocate(atmat)
      call nwdel(nwall,4*nasht**4*8)
      deallocate(bmat)
      deallocate(cmat)
      deallocate(btmat)
      deallocate(ctmat)
      call nwdel(nwall,2*nasht**2*8)
      deallocate(dmat)
      deallocate(dtmat)

      call dypc(escf,onel,id2n,in2d,WRK,LWRK,
     $           norbt,nisht,nasht,its,itsym,NOCB,
     $           LUPRI,LUN98,LUN31,zgel)

      CALL GPCLOSE(LUN31,'DELETE')
      call nwdel(nwall,2*norbt*4+norbt)
      deallocate (id2n)
      deallocate (in2d)
      deallocate (zgel)
      call nwdel(nwall,(norbt+norbt+1)*4)
      deallocate (itsym)
      call nwdel(nwall,(norbt+norbt+1)*4)
      deallocate (aj)
      call nwdel(nwall,NNORBX*8)
      deallocate (onel)
c----that's all folks
      CALL QEXIT('koopro4')
      return
      end
c-------------------------------------------------------------------
      subroutine rofour(c,amat,nisht,nasht,thr,icomp,iconf,ioccfe,
     $     zcapos,ncaptot,trou,part,ne,nd,ispinfe,iorbfe,itsym,iocca
     $     ,ioccb,flip)
      IMPLICIT REAL*8  (A-H,O-Y), LOGICAL*1 (Z)
#include <priunit.h>
      dimension amat(nwords,nasht,nasht,nasht,nasht)
#include "nevpt2.h"
      dimension c(*),zcapos(*)
      integer*2 icomp(*),iconf(*)
      integer*1 ioccfe(norb,ncaptot),flip
      dimension iocca(*),ioccb(*),flip(*)
      COMMON /CIP/NORB,NOCB,NCF,ZSEG
      dimension ne(*),nd(*),ispinfe(*),iorbfe(*),itsym(*)
      INTEGER*2 ispinfe,iorbfe
      INTEGER*2 NE,TROU,PART
      integer*8 nhits,ijcal
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      integer op,pp,a,b,op1,op2,op3,pp1,pp2,pp3
      integer segno,segno0,segno1,segno2,segno3,segno4,segno5
      nocc=NOCB-nisht
      do i=1,ncf
       flip(i)=0
      enddo
      nfours=0
      ntrips=0
      ndoubs=0
      nsings=0
      nhits=0
c     nhitsd=0
      mnext=1
      iconfi=1
      mcapos=0
      znocalc=.false.
      mfin=0
      ijcal=0
      call esame(flip,ne,nd,trou,part,ispinfe,iorbfe,ncf)
      iwr=0
      do  mcap=1,ncaptot
         minit=mfin+1
         mfin=minit+icomp(minit)-1
         if (mod(mcap,100).eq.0)then
            iwr=iwr+1
            imd=mod(iwr,8)
            if(imd.eq.1)then
               write(lupri,'(a,i8,a,$)')'mcap=',mcap,','
               call flshfo(lupri)
            elseif(imd.gt.0)then
               write(lupri,'(i8,a,$)')mcap,','
               call flshfo(lupri)
            elseif(imd.eq.0)then
               write(lupri,'(i8)')mcap
               call flshfo(lupri)
            endif
         endif
         nec=ne(minit)
         ndm=nd(minit)
         nfin=0
         do ncap = 1,mcap      !do 2
            ninit=nfin+1
            nfin=ninit+icomp(ninit)-1
            nedif=abs(ne(minit)-ne(ninit))
            znocalc=nedif.gt.4
c           if(znocalc)nhitsd=nhitsd+1
            if(.not.znocalc)znocalc=ndiff(ioccfe(nisht+1,mcap)
     $           ,ioccfe(nisht+1,ncap),nasht).gt.8
            if(znocalc)nhits=nhits+1
            if(znocalc)goto 2
            do m=minit,mfin
               nec=ne(m)
               ndm=nd(m)
               call giveocc2(m,iocca,ioccb,nasht,nisht,nocca,noccb,nocc
     $              ,nec,ndm,trou,part,iorbfe,ispinfe) 
               if(abs(c(m)).lt.thr) goto 3
               if(flip(m).eq.2)goto 3
               if(mcap.eq.ncap)then
                  nbeg=minit
                  nend=m-1
c                  call case0
c     write(lupri,*)'case 0 m=',m
c---  contributo a daaaa
                  if(flip(m).eq.1)then
                     valcm=2.d0*c(m)
                  else
                     valcm=c(m)
                  endif
                  cij=valcm*c(m)
                  val=cij
                  do i=1,nocca
                     ia=iocca(i)
                     do j=i+1,nocca
                        ic=iocca(j)
                        do k=j+1,nocca
                           ie=iocca(k)
                           do l=k+1,nocca
                              ig=iocca(l)
                              call ordsame(ia,ic,ie,ig,ia,ic,ie,ig,val
     $                             ,amat,nasht,m,m)
                           enddo
                        enddo
                     enddo
                  enddo
c---  contributo a dbbbb
                  do i=1,noccb
                     ia=ioccb(i)
                     do j=i+1,noccb
                        ic=ioccb(j)
                        do k=j+1,noccb
                           ie=ioccb(k)
                           do l=k+1,noccb
                              ig=ioccb(l)
                              call ordsame(ia,ic,ie,ig,ia,ic,ie,ig,val
     $                             ,amat,nasht,m,m)
                           enddo
                        enddo
                     enddo
                  enddo
c---  contributi a daaab,daaba,dabaa,dbaaa
                  do i=1,nocca
                     ia=iocca(i)
                     do j=i+1,nocca
                        ic=iocca(j)
                        do k=j+1,nocca
                           ie=iocca(k)
                           do l=1,noccb
                              ig=ioccb(l)
                              call ord31(ia,ic,ie,ig,ia,ic,ie,ig,val
     $                             ,amat,nasht,m,m)
                           enddo
                        enddo
                     enddo
                  enddo
c---  contributi a dbbba,dbbab,dbabb,dabbb
                  do i=1,noccb
                     ia=ioccb(i)
                     do j=i+1,noccb
                        ic=ioccb(j)
                        do k=j+1,noccb
                           ie=ioccb(k)
                           do l=1,nocca
                              ig=iocca(l)
                              call ord31(ia,ic,ie,ig,ia,ic,ie,ig,val
     $                             ,amat,nasht,m,m)
                           enddo
                        enddo
                     enddo
                  enddo
c---  contributo a daabb,dbbaa,dabab,dbaba,dabba,dbaab
                  do i=1,nocca
                     ia=iocca(i)
                     do j=i+1,nocca
                        ic=iocca(j)
                        do k=1,noccb
                           ie=ioccb(k)
                           do l=k+1,noccb
                              ig=ioccb(l)
                              call ord22(ia,ic,ie,ig,ia,ic,ie,ig,val
     $                             ,amat,nasht,m,m)
                           enddo
                        enddo
                     enddo
                  enddo
               else
                  nbeg=ninit
                  nend=nfin
               endif
               do n=nbeg,nend
                  if(abs(c(n)).lt.thr) goto 4
                  ij=jdifgen(ne,nd,trou,part,m,n,nv,nw,4,zseg)
                  ijcal=ijcal+1
                  if(ij.eq.0)nbdif=5
                  if(nbdif.gt.4)goto 4
                  if(zseg)then
                     segno=-1
                  else
                     segno=1
                  endif
                  if(nbdif.eq.1)then
                     nsings=nsings+1
                     n1=iorbfe(nv(1))
                     ns1=ispinfe(nv(1))
                     n2=iorbfe(nw(1))
                     ns2=ispinfe(nw(1))
c                     call case10
                     n1=n1-nisht
                     n2=n2-nisht
                     if(flip(m).eq.1.and.n.ne.m-1)then
                        valcm=2.d0*c(m)
                     else
                        valcm=c(m)
                     endif
                     cij=valcm*c(n)*segno
                     val=cij
                     if(ns1.eq.1.and.ns2.eq.1)then
c---  a sinistra il piu' grande
                        if(n1.ge.n2)then
                           nl=n1
                           nr=n2
                        else
                           nl=n2
                           nr=n1
                        endif
c---  contributo a daaaa
                        do i=1,nocca
                           ia=iocca(i)
                           do j=i+1,nocca
                              ic=iocca(j)
                              do k=j+1,nocca
                                 ie=iocca(k)
                                 if(ia.ne.n1.and.ic.ne.n1.and.ie.ne.n1
     $                                .and.ia.ne.n2.and.ic.ne.n2.and.ie
     $                                .ne.n2)then
                                    call ordsame(ia,ic,ie,nl,ia,ic,ie,nr
     $                                   ,val,amat,nasht,m,n)
                                 endif
                              enddo
                           enddo
                        enddo
c---  contributo a dbbba
                        do i=1,noccb
                           ia=ioccb(i)
                           do j=i+1,noccb
                              ic=ioccb(j)
                              do k=j+1,noccb
                                 ie=ioccb(k)
c     write (6,*) 'calling ord31'
                                 call ord31(ia,ic,ie,nl,ia,ic,ie,nr,val
     $                                ,amat,nasht,m,n)
                              enddo
                           enddo
                        enddo
c---  contributi a daaab,daaba,dabaa,dbaaa
                        do i=1,nocca
                           ia=iocca(i)
                           do j=i+1,nocca
                              ic=iocca(j)
                              if(ia.ne.n1.and.ic.ne.n1.and.ia.ne.n2.and
     $                             .ic.ne.n2)then
                                 do k=1,noccb
                                    ie=ioccb(k)
c     write (6,*) 'calling ord31 aaab'
c     write (6,*) ia,ic,ie,nl,nr
                                    call ord31(ia,ic,nl,ie,ia,ic,nr,ie
     $                                   ,val,amat,nasht,m,n)
                                 enddo
                              endif
                           enddo
                        enddo
c---  contributi a daabb,dbbaa,dabab,dbaba,dabba,dbaab
                        do i=1,nocca
                           ia=iocca(i)
                           do j=1,noccb
                              ic=ioccb(j)
                              do k=j+1,noccb
                                 ie=ioccb(k)
                                 if(ia.ne.n1.and.ia.ne.n2)then
c     write (6,*) 'calling ord22'
                                    call ord22(ia,nl,ic,ie,ia,nr,ic,ie
     $                                   ,val,amat,nasht,m,n)
                                 endif
                              enddo
                           enddo
                        enddo
                     elseif(ns1.eq.0.and.ns2.eq.0)then
c---  a sinistra il piu' grande
                        if(n1.ge.n2)then
                           nl=n1
                           nr=n2
                        else
                           nl=n2
                           nr=n1
                        endif
c---  contributo a dbbbb
                        do i=1,noccb
                           ia=ioccb(i)
                           do j=i+1,noccb
                              ic=ioccb(j)
                              do k=j+1,noccb
                                 ie=ioccb(k)
                                 if(ia.ne.n1.and.ic.ne.n1.and.ie.ne.n1
     $                                .and.ia.ne.n2.and.ic.ne.n2.and.ie
     $                                .ne.n2)then
c     write (6,*) 'calling ordsame'
                                    call ordsame(ia,ic,ie,nl,ia,ic,ie,nr
     $                                   ,val,amat,nasht,m,n)
                                 endif
                              enddo
                           enddo
                        enddo
c---  contributo a daaab
                        do i=1,nocca
                           ia=iocca(i)
                           do j=i+1,nocca
                              ic=iocca(j)
                              do k=j+1,nocca
                                 ie=iocca(k)
c     write (6,*) 'calling ord31'
                                 call ord31(ia,ic,ie,nl,ia,ic,ie,nr,val
     $                                ,amat,nasht,m,n)
                              enddo
                           enddo
                        enddo
c---  contributi a dbbba,dbbab,dbabb,dabbb
                        do i=1,noccb
                           ia=ioccb(i)
                           do j=i+1,noccb
                              ic=ioccb(j)
                              if(ia.ne.n1.and.ic.ne.n1.and.ia.ne.n2.and
     $                             .ic.ne.n2)then
                                 do k=1,nocca
                                    ie=iocca(k)
c     write (6,*) 'calling ord31'
                                    call ord31(ia,ic,nl,ie,ia,ic,nr,ie
     $                                   ,val,amat,nasht,m,n)
                                 enddo
                              endif
                           enddo
                        enddo
c---  contributi a dbbaa,daabb,dbaba,dabab,dbaab,dabba
                        do i=1,noccb
                           ia=ioccb(i)
                           do j=1,nocca
                              ic=iocca(j)
                              do k=j+1,nocca
                                 ie=iocca(k)
                                 if(ia.ne.n1.and.ia.ne.n2)then
c     write (6,*) 'calling ord22'
                                    call ord22(ic,ie,ia,nl,ic,ie,ia,nr
     $                                   ,val,amat,nasht,m,n)
                                 endif
                              enddo
                           enddo
                        enddo
                     endif
c     goto 2
                  elseif(nbdif.eq.2)then
                     ndoubs=ndoubs+1
                     call ord2(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                     n1=iorbfe(nv(1))
                     ns1=ispinfe(nv(1))
                     n2=iorbfe(nv(2))
                     ns2=ispinfe(nv(2))
                     n3=iorbfe(nw(1))
                     ns3=ispinfe(nw(1))
                     n4=iorbfe(nw(2))
                     ns4=ispinfe(nw(2))
c     goto 20
c                     call case20
                     n1=n1-nisht
                     n2=n2-nisht
                     n3=n3-nisht
                     n4=n4-nisht
                     if(flip(m).eq.1.and.n.ne.m-1)then
                        valcm=2.d0*c(m)
                     else
                        valcm=c(m)
                     endif
                     cij=valcm*c(n)*segno
                     val=cij
c---  n1l,n2l a sinistra 
                     if((n1+n2).ge.(n3+n4))then
                        n1l=n1
                        n2l=n2
                        n1r=n3
                        n2r=n4
                     else
                        n1l=n3
                        n2l=n4
                        n1r=n1
                        n2r=n2
                     endif
                     if(ns1.eq.1.and.ns2.eq.1)then
c---  contributo a daaaa
                        do i=1,nocca
                           ia=iocca(i)
                           do j=i+1,nocca
                              ic=iocca(j)
                              if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and
     $                             .ia.ne.n4.and.ic.ne.n1.and.ic.ne.n2
     $                             .and.ic.ne.n3.and.ic.ne.n4)then
                                 call ordsame(ia,ic,n1l,n2l,ia,ic,n1r
     $                                ,n2r,val,amat,nasht,m,n)
                              endif
                           enddo
                        enddo
c---  contributo a daaab,daaba,dabaa,dbaaa
                        do i=1,nocca
                           ia=iocca(i)
                           if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia
     $                          .ne.n4)then
                              do j=1,noccb
                                 ic=ioccb(j)
                                 call ord31(ia,n1l,n2l,ic,ia,n1r,n2r,ic
     $                                ,val,amat,nasht,m,n)
                              enddo
                           endif
                        enddo
c---  contributo a dbbaa,daabb,dabab,dbaba,dabba,dbaab
                        do i=1,noccb
                           ia=ioccb(i)
                           do j=i+1,noccb
                              ic=ioccb(j)
                              call ord22(n1l,n2l,ia,ic,n1r,n2r,ia,ic,val
     $                             ,amat,nasht,m,n)
                           enddo
                        enddo
                     elseif(ns1.eq.0.and.ns2.eq.0)then
c---  contributo a dbbbb
                        do i=1,noccb
                           ia=ioccb(i)
                           do j=i+1,noccb
                              ic=ioccb(j)
                              if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and
     $                             .ia.ne.n4.and.ic.ne.n1.and.ic.ne.n2
     $                             .and.ic.ne.n3.and.ic.ne.n4)then
                                 call ordsame(ia,ic,n1l,n2l,ia,ic,n1r
     $                                ,n2r,val,amat,nasht,m,n)
                              endif
                           enddo
                        enddo
c---  contributo a dbbba,dbbab,dbabb,dabbb
                        do i=1,noccb
                           ia=ioccb(i)
                           if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia
     $                          .ne.n4)then
                              do j=1,nocca
                                 ic=iocca(j)
                                 call ord31(ia,n1l,n2l,ic,ia,n1r,n2r,ic
     $                                ,val,amat,nasht,m,n)
                              enddo
                           endif
                        enddo
c---  contributo a daabb,dbbaa,dbaba,dabab,dbaab,dabba
                        do i=1,nocca
                           ia=iocca(i)
                           do j=i+1,nocca
                              ic=iocca(j)
                              call ord22(ia,ic,n1l,n2l,ia,ic,n1r,n2r,val
     $                             ,amat,nasht,m,n)
                           enddo
                        enddo
                     elseif(ns1.eq.1.and.ns2.eq.0)then
c---  contributi a daaab,daaba,dabaa,dbaaa
                        do i=1,nocca
                           ia=iocca(i)
                           do j=i+1,nocca
                              ic=iocca(j)
                              if(ia.ne.n1.and.ia.ne.n3.and.ic.ne.n1.and
     $                             .ic.ne.n3)then
c---- daaab
                                 call ord31(ia,ic,n1l,n2l,ia,ic,n1r,n2r
     $                                ,val,amat,nasht,m,n)
                              endif
                           enddo
                        enddo
c---  contributi a dbbba,dbbab,dbabb,dabbb
                        do i=1,noccb
                           ia=ioccb(i)
                           do j=i+1,noccb
                              ic=ioccb(j)
                              if(ia.ne.n2.and.ia.ne.n4.and.ic.ne.n2.and
     $                             .ic.ne.n4)then
c---- dbbba
                                 call ord31(ia,ic,n2l,n1l,ia,ic,n2r,n1r
     $                                ,val,amat,nasht,m,n)
                              endif
                           enddo
                        enddo
c---  contributi a daabb,dbbaa,dabab,dbaba,dabba,dbaab
                        do i=1,nocca
                           ia=iocca(i)
                           do j=1,noccb
                              ic=ioccb(j)
                              if(ia.ne.n1.and.ia.ne.n3.and.ic.ne.n2.and
     $                             .ic.ne.n4)then
c---- daabb
                                 call ord22(ia,n1l,ic,n2l,ia,n1r,ic,n2r
     $                                ,val,amat,nasht,m,n)
                              endif
                           enddo
                        enddo    
                     endif
c     goto 2
                  elseif(nbdif.eq.3)then
                     ntrips=ntrips+1
                     call ord3(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                     n1=iorbfe(nv(1))
                     ns1=ispinfe(nv(1))
                     n2=iorbfe(nv(2))
                     ns2=ispinfe(nv(2))
                     n3=iorbfe(nv(3))
                     ns3=ispinfe(nv(3))
                     n4=iorbfe(nw(1))
                     ns4=ispinfe(nw(1))
                     n5=iorbfe(nw(2))
                     ns5=ispinfe(nw(2))
                     n6=iorbfe(nw(3))
                     ns6=ispinfe(nw(3))
c     goto 30
c                     call case30
                     n1=n1-nisht
                     n2=n2-nisht
                     n3=n3-nisht
                     n4=n4-nisht
                     n5=n5-nisht
                     n6=n6-nisht
                     if(flip(m).eq.1.and.n.ne.m-1)then
                        valcm=2.d0*c(m)
                     else
                        valcm=c(m)
                     endif
                     cij=valcm*c(n)*segno
                     val=cij
c---  n1l,n2l,n3l a sinistra 
                     if((n1+n2+n3).ge.(n4+n5+n6))then
                        n1l=n1
                        n2l=n2
                        n3l=n3
                        n1r=n4
                        n2r=n5
                        n3r=n6
                     else
                        n1l=n4
                        n2l=n5
                        n3l=n6
                        n1r=n1
                        n2r=n2
                        n3r=n3
                     endif
c     casi aaa, bbb, aab, bba
                     if(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1)then
c---  contributo daaaa
                        do i=1,nocca
                           ia=iocca(i)
                           if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia
     $                          .ne.n4.and.ia.ne.n5.and.ia.ne.n6)then
                              call ordsame(ia,n1l,n2l,n3l,ia,n1r,n2r,n3r
     $                             ,val,amat,nasht,m,n)
                           endif
                        enddo
c---  contributi a daaab,daaba,dabaa,dbaaa
                        do i=1,noccb
                           ia=ioccb(i)
                           call ord31(n1l,n2l,n3l,ia,n1r,n2r,n3r,ia,val
     $                          ,amat,nasht,m,n)
                        enddo
                     elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0)then
c---  contributo dbbbb
                        do i=1,noccb
                           ia=ioccb(i)
                           if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia
     $                          .ne.n4.and.ia.ne.n5.and.ia.ne.n6)then
                              call ordsame(ia,n1l,n2l,n3l,ia,n1r,n2r,n3r
     $                             ,val,amat,nasht,m,n)
                           endif
                        enddo
c---  contributi a dbbba,dbbab,dbabb,dabbb
                        do i=1,nocca
                           ia=iocca(i)
                           call ord31(n1l,n2l,n3l,ia,n1r,n2r,n3r,ia,val
     $                          ,amat,nasht,m,n)
                        enddo
                     elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0)then
c---  contributi a daaab,daaba,dabaa,dbaaa
                        do i=1,nocca
                           ia=iocca(i)
                           if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n4.and.ia
     $                          .ne.n5)then
c---- daaab
                              call ord31(ia,n1l,n2l,n3l,ia,n1r,n2r,n3r
     $                             ,val,amat,nasht,m,n)
                           endif
                        enddo
c---  contributi a daabb,dbbaa,dabab,dbaba,dabba,dbaab
                        do i=1,noccb
                           ia=ioccb(i)
                           if(ia.ne.n3.and.ia.ne.n6)then
c---- daabb
                              call ord22(n1l,n2l,ia,n3l,n1r,n2r,ia,n3r
     $                             ,val,amat,nasht,m,n)
                           endif
                        enddo
                     elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0)then
c---  contributi a dbbba,dbbab,dbabb,dabbb
                        do i=1,noccb
                           ia=ioccb(i)
                           if(ia.ne.n2.and.ia.ne.n3.and.ia.ne.n5.and.ia
     $                          .ne.n6)then
c---- dbbba
                              call ord31(ia,n2l,n3l,n1l,ia,n2r,n3r,n1r
     $                             ,val,amat,nasht,m,n)
                           endif
                        enddo
c---  contributi a dbbaa,daabb,dbaba,dabab,dbaab,dabba
                        do i=1,nocca
                           ia=iocca(i)
                           if(ia.ne.n1.and.ia.ne.n4)then
c---- dbbaa
                              call ord22(n1l,ia,n2l,n3l,n1r,ia,n2r,n3r
     $                             ,val,amat,nasht,m,n)
                           endif
                        enddo
                     endif
c     goto 2
                  else          !nbdif.eq.4
c---  renzo debug
                     nfours=nfours+1
c---- 
                     call ord4(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                     n1=iorbfe(nv(1))
                     ns1=ispinfe(nv(1))
                     n2=iorbfe(nv(2))
                     ns2=ispinfe(nv(2))
                     n3=iorbfe(nv(3))
                     ns3=ispinfe(nv(3))
                     n4=iorbfe(nv(4))
                     ns4=ispinfe(nv(4))
                     n5=iorbfe(nw(1))
                     ns5=ispinfe(nw(1))
                     n6=iorbfe(nw(2))
                     ns6=ispinfe(nw(2))
                     n7=iorbfe(nw(3))
                     ns7=ispinfe(nw(3))
                     n8=iorbfe(nw(4))
                     ns8=ispinfe(nw(4))
c                     goto 40
c                  call case40
                     n1=n1-nisht
                     n2=n2-nisht
                     n3=n3-nisht
                     n4=n4-nisht
                     n5=n5-nisht
                     n6=n6-nisht
                     n7=n7-nisht
                     n8=n8-nisht
                     if(flip(m).eq.1.and.n.ne.m-1)then
                        valcm=2.d0*c(m)
                     else
                        valcm=c(m)
                     endif
                     cij=valcm*c(n)*segno
                     val=cij
c---  n1l,n2l,n3l a sinistra 
                     if((n1+n2+n3+n4).ge.(n5+n6+n7+n8))then
                        n1l=n1
                        n2l=n2
                        n3l=n3
                        n4l=n4
                        n1r=n5
                        n2r=n6
                        n3r=n7
                        n4r=n8
                     else
                        n1l=n5
                        n2l=n6
                        n3l=n7
                        n4l=n8
                        n1r=n1
                        n2r=n2
                        n3r=n3
                        n4r=n4
                     endif
                     if(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1.and.ns4.eq.1
     $                    )then
c     contributo a daaaa
                        call ordsame(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val
     $                       ,amat,nasht,m,n)
                     elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0.and.ns4
     $                       .eq.0)then
c     contributo a dbbbb
                        call ordsame(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val
     $                       ,amat,nasht,m,n)
                     elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1.and.ns4
     $                       .eq.0)then
c     contributo a daaab
                        call ord31(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val
     $                       ,amat,nasht,m,n)
                     elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0.and.ns4
     $                       .eq.0)then
c     contributo a dbbba
                        call ord31(n2l,n3l,n4l,n1l,n2r,n3r,n4r,n1r,val
     $                       ,amat,nasht,m,n)
                     elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0.and.ns4
     $                       .eq.0)then
c     contributo a daabb
                        call ord22(n1l,n2l,n3l,n4l,n1r,n2r,n3r,n4r,val
     $                       ,amat,nasht,m,n)
                     endif
                  endif
 4             continue
               enddo
 3          continue
            enddo
 2       continue
         enddo
 1    continue
      enddo
c---renzo debug
      write(lupri,*)
      write(lupri,*)'number of single    differences=',nsings
      write(lupri,*)'number of double    differences=',ndoubs
      write(lupri,*)'number of triple    differences=',ntrips
      write(lupri,*)'number of quadruple differences=',nfours
      write(lupri,*)'number of jdifgen calls         ',ijcal
c---contributi diagonali
      return
      end
c--------------------------------------------------------------------
      subroutine ro3(c,daaa,nisht,nasht,thr,icomp,iconf,ioccfe,zcapos
     $,ncaptot,trou,part,ne,nd,ispinfe,iorbfe,itsym,iocca,ioccb,flip)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1 (Z)
#include <priunit.h>
c      include 'rotwo.par'
      allocatable dbbb(:,:,:,:,:,:),daab(:,:,:,:,:,:),dbba(:,:,:,:,:,:)
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht)
      dimension c(*),zcapos(*)
      COMMON /CIP/NORB,NOCB,NCF,ZSEG
      INTEGER*2 ispinfe,iorbfe
      INTEGER*2 NE,TROU,PART
      integer*2 icomp(*),iconf(*)
      integer*1 ioccfe(norb,ncaptot)
      dimension ne(*),nd(*),ispinfe(*),iorbfe(*),itsym(*)
      integer segno
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      integer*1 flip
      dimension iocca(*),ioccb(*),flip(*)
      nocc=NOCB-nisht
      write(lupri,*)'subr. ro3'
      write(lupri,*)'nocc,NOCB',nocc,NOCB
      write(lupri,*)'nasht=',nasht,' nocc=',nocc,' ncf=',ncf

      allocate(dbbb(1:nasht,1:nasht,1:nasht,1:nasht,1:nasht,1:nasht))
      allocate(daab(1:nasht,1:nasht,1:nasht,1:nasht,1:nasht,1:nasht))
      allocate(dbba(1:nasht,1:nasht,1:nasht,1:nasht,1:nasht,1:nasht))
c---azzeramento
      do i=1,nasht
         do j=1,nasht
            do k=1,nasht
               do l=1,nasht
                  do m=1,nasht
                     do n=1,nasht
                        daaa(i,j,k,l,m,n)=0.d0
                        dbbb(i,j,k,l,m,n)=0.d0
                        daab(i,j,k,l,m,n)=0.d0
                        dbba(i,j,k,l,m,n)=0.d0
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      do i=1,ncf
       flip(i)=0
      enddo
      mnext=1
      iconfi=1
      mcapos=0
      znocalc=.false.
      mfin=0
      ijcal=0
      call esame(flip,ne,nd,trou,part,ispinfe,iorbfe,ncf)
      iwr=0
      do  mcap=1,ncaptot
         minit=mfin+1
         mfin=minit+icomp(minit)-1
         if (mod(mcap,100).eq.0)then
            iwr=iwr+1
            imd=mod(iwr,8)
            if(imd.eq.1)then
               write(lupri,'(a,i8,a,$)')'mcap=',mcap,','
               call flush(lupri)
            elseif(imd.gt.0)then
               write(lupri,'(i8,a,$)')mcap,','
               call flush(lupri)
            elseif(imd.eq.0)then
               write(lupri,'(i8)')mcap
               call flush(lupri)
            endif
         endif
         nec=ne(minit)
         ndm=nd(minit)
         nfin=0
         do  ncap=1,mcap        !do 2
            ninit=nfin+1
            nfin=ninit+icomp(ninit)-1
            nedif=abs(ne(minit)-ne(ninit))
            znocalc=nedif.gt.3
c           if(znocalc)nhitsd=nhitsd+1
            if(.not.znocalc)znocalc=ndiff(ioccfe(nisht+1,mcap)
     $           ,ioccfe(nisht+1,ncap),nasht).gt.6
            if(znocalc)goto 2
            do m=minit,mfin
               nec=ne(m)
               ndm=nd(m)
               call giveocc2(m,iocca,ioccb,nasht,nisht,nocca,noccb,nocc
     $              ,nec,ndm,trou,part,iorbfe,ispinfe) 
               if(abs(c(m)).lt.thr) goto 3
               if(flip(m).eq.2)goto 3
               if(mcap.eq.ncap)then
                  nbeg=minit
                  nend=m-1
c     call case0
c---- renzo 
                  if(flip(m).eq.1)then
                     valcm=2.d0*c(m)
                  else
                     valcm=c(m)
                  endif
                  cij=valcm*c(m)
                  val=cij
c     costruzione di daaa
                  do i=1,nocca
                     ia=iocca(i)
                     do j=i+1,nocca
                        ja=iocca(j)
                        do k=j+1,nocca
                           ka=iocca(k)
                           daaa(ia,ja,ka,ia,ja,ka)=daaa(ia,ja,ka,ia,ja
     $                          ,ka)+val
                           call perm3(daaa,ia,ja,ka,ia,ja,ka,nasht)
                        enddo
                     enddo
                  enddo
c     costruzione di dbbb            
                  do i=1,noccb
                     ib=ioccb(i)
                     do j=i+1,noccb
                        jb=ioccb(j)
                        do k=j+1,noccb
                           kb=ioccb(k)
                           dbbb(ib,jb,kb,ib,jb,kb)=dbbb(ib,jb,kb,ib,jb
     $                          ,kb)+val
                           call perm3(dbbb,ib,jb,kb,ib,jb,kb,nasht)
                        enddo
                     enddo
                  enddo
c     costruzione di daab
                  do i=1,nocca
                     ia=iocca(i)
                     do j=i+1,nocca
                        ja=iocca(j)
                        do k=1,noccb
                           kb=ioccb(k)
                           daab(ia,ja,kb,ia,ja,kb)=daab(ia,ja,kb,ia,ja
     $                          ,kb)+val
                           call perm2(daab,ia,ja,kb,ia,ja,kb,nasht)
                        enddo
                     enddo
                  enddo
c     costruzione di dbba            
                  do i=1,noccb
                     ib=ioccb(i)
                     do j=i+1,noccb
                        jb=ioccb(j)
                        do k=1,nocca
                           ka=iocca(k)
                           dbba(ib,jb,ka,ib,jb,ka)=dbba(ib,jb,ka,ib,jb
     $                          ,ka)+val
                           call perm2(dbba,ib,jb,ka,ib,jb,ka,nasht)
                        enddo
                     enddo
                  enddo
               else
                  nbeg=ninit
                  nend=nfin
               endif
               do n=nbeg,nend
                  if(abs(c(n)).lt.thr) goto 4
                  ij=jdifgen(ne,nd,trou,part,m,n,nv,nw,3,zseg)
                  if(ij.eq.0)nbdif=5
                  if(nbdif.gt.3)goto 4
                  if(zseg)then
                     segno=-1
                  else
                     segno=1
                  endif
                  if(nbdif.eq.1)then
                     n1=iorbfe(nv(1))
                     ns1=ispinfe(nv(1))
                     n2=iorbfe(nw(1))
                     ns2=ispinfe(nw(1))
                     n1=n1-nisht
                     n2=n2-nisht
                     if(flip(m).eq.1.and.n.ne.m-1)then
                        valcm=2.d0*c(m)
                     else
                        valcm=c(m)
                     endif
                     cij=valcm*c(n)*segno
                     val=cij
c     cij=c(m)*c(n)*segno
                     if(ns1.eq.1.and.ns2.eq.1)then
c---  contributo a daaa
                        do i=1,nocca
                           ia=iocca(i)
                           do j=i+1,nocca
                              ja=iocca(j)
                              if(ia.ne.n1.and.ja.ne.n1.and.ia.ne.n2.and
     $                             .ja.ne.n2)then
                                 daaa(n1,ia,ja,n2,ia,ja)=daaa(n1,ia,ja
     $                                ,n2,ia,ja)+val
                                 daaa(n2,ia,ja,n1,ia,ja)=daaa(n2,ia,ja
     $                                ,n1,ia,ja)+val
                                 call perm3(daaa,n1,ia,ja,n2,ia,ja,nasht
     $                                )
                                 call perm3(daaa,n2,ia,ja,n1,ia,ja,nasht
     $                                )
                              endif
                           enddo
                        enddo
c---- contributo a dbba
                        do i=1,noccb
                           ib=ioccb(i)
                           do j=i+1,noccb
                              jb=ioccb(j)
                              dbba(ib,jb,n1,ib,jb,n2)=dbba(ib,jb,n1,ib
     $                             ,jb,n2)+val
                              call perm2(dbba,ib,jb,n1,ib,jb,n2,nasht)
                           enddo
                        enddo
c---- contributo a daab
                        do i=1,nocca
                           ia=iocca(i)
                           if(ia.ne.n1.and.ia.ne.n2)then
                              do j=1,noccb
                                 jb=ioccb(j)
                                 daab(n1,ia,jb,n2,ia,jb)=daab(n1,ia,jb
     $                                ,n2,ia,jb)+val
                                 call perm2(daab,n1,ia,jb,n2,ia,jb,nasht
     $                                )
                              enddo
                           endif
                        enddo
                     elseif(ns1.eq.0.and.ns2.eq.0)then
c---- contributo a dbbb
                        do i=1,noccb
                           ib=ioccb(i)
                           do j=i+1,noccb
                              jb=ioccb(j)
                              if(ib.ne.n1.and.jb.ne.n1.and.ib.ne.n2.and
     $                             .jb.ne.n2)then
                                 dbbb(n1,ib,jb,n2,ib,jb)=dbbb(n1,ib,jb
     $                                ,n2,ib,jb)+val
                                 dbbb(n2,ib,jb,n1,ib,jb)=dbbb(n2,ib,jb
     $                                ,n1,ib,jb)+val
                                 call perm3(dbbb,n1,ib,jb,n2,ib,jb,nasht
     $                                )
                                 call perm3(dbbb,n2,ib,jb,n1,ib,jb,nasht
     $                                )
                              endif
                           enddo
                        enddo
c---  contrinuto a daab
                        do i=1,nocca
                           ia=iocca(i)
                           do j=i+1,nocca
                              ja=iocca(j)
                              daab(ia,ja,n1,ia,ja,n2)=daab(ia,ja,n1,ia
     $                             ,ja,n2)+val
                              call perm2(daab,ia,ja,n1,ia,ja,n2,nasht)
                           enddo
                        enddo
c---- contributo a dbba
                        do i=1,noccb
                           ib=ioccb(i)
                           if(ib.ne.n1.and.ib.ne.n2)then
                              do j=1,nocca
                                 ja=iocca(j)
                                 dbba(n1,ib,ja,n2,ib,ja)=dbba(n1,ib,ja
     $                                ,n2,ib,ja)+val
                                 call perm2(dbba,n1,ib,ja,n2,ib,ja,nasht
     $                                )
                              enddo
                           endif
                        enddo
                     endif
                  elseif(nbdif.eq.2)then
                     call ord2(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                     n1=iorbfe(nv(1))
                     ns1=ispinfe(nv(1))
                     n2=iorbfe(nv(2))
                     ns2=ispinfe(nv(2))
                     n3=iorbfe(nw(1))
                     ns3=ispinfe(nw(1))
                     n4=iorbfe(nw(2))
                     ns4=ispinfe(nw(2))
                     n1=n1-nisht
                     n2=n2-nisht
                     n3=n3-nisht
                     n4=n4-nisht
                     if(flip(m).eq.1.and.n.ne.m-1)then
                        valcm=2.d0*c(m)
                     else
                        valcm=c(m)
                     endif
                     cij=valcm*c(n)*segno
                     val=cij
                     if(ns1.eq.1.and.ns2.eq.1)then
                        do i=1,nocca
                           ia=iocca(i)
                           if(ia.ne.n1.and.ia.ne.n2.and.ia.ne.n3.and.ia
     $                          .ne.n4)then
                              daaa(n1,n2,ia,n3,n4,ia)=daaa(n1,n2,ia,n3
     $                             ,n4,ia)+val
                              daaa(n3,n4,ia,n1,n2,ia)=daaa(n3,n4,ia,n1
     $                             ,n2,ia)+val
                              call perm3(daaa,n1,n2,ia,n3,n4,ia,nasht)
                              call perm3(daaa,n3,n4,ia,n1,n2,ia,nasht)
                           endif
                        enddo
                        do i=1,noccb
                           ib=ioccb(i)
                           daab(n1,n2,ib,n3,n4,ib)=daab(n1,n2,ib,n3,n4
     $                          ,ib)+val
                           call perm2(daab,n1,n2,ib,n3,n4,ib,nasht)
                        enddo
                        
                     elseif(ns1.eq.0.and.ns2.eq.0)then
                        do i=1,noccb
                           ib=ioccb(i)
                           if(ib.ne.n1.and.ib.ne.n2.and.ib.ne.n3.and.ib
     $                          .ne.n4)then
                              dbbb(n1,n2,ib,n3,n4,ib)=dbbb(n1,n2,ib,n3
     $                             ,n4,ib)+val
                              dbbb(n3,n4,ib,n1,n2,ib)=dbbb(n3,n4,ib,n1
     $                             ,n2,ib)+val
                              call perm3(dbbb,n1,n2,ib,n3,n4,ib,nasht)
                              call perm3(dbbb,n3,n4,ib,n1,n2,ib,nasht)
                           endif
                        enddo
                        do i=1,nocca
                           ia=iocca(i)
                           dbba(n1,n2,ia,n3,n4,ia)=dbba(n1,n2,ia,n3,n4
     $                          ,ia)+val
                           call perm2(dbba,n1,n2,ia,n3,n4,ia,nasht)
                        enddo
                     elseif(ns1.eq.1.and.ns2.eq.0)then
                        do i=1,nocca
                           ia=iocca(i)
                           if(ia.ne.n1.and.ia.ne.n3)daab(n1,ia,n2,n3,ia
     $                          ,n4)=daab(n1,ia,n2,n3,ia,n4)+val
                           call perm2(daab,n1,ia,n2,n3,ia,n4,nasht)
                        enddo
                        do i=1,noccb
                           ib=ioccb(i)
                           if(ib.ne.n2.and.ib.ne.n4)dbba(n2,ib,n1,n4,ib
     $                          ,n3)=dbba(n2,ib,n1,n4,ib,n3)+val
                           call perm2(dbba,n2,ib,n1,n4,ib,n3,nasht)
                        enddo
                     endif
                  else   !nbdif is 3
                     call ord3(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                     n1=iorbfe(nv(1))
                     ns1=ispinfe(nv(1))
                     n2=iorbfe(nv(2))
                     ns2=ispinfe(nv(2))
                     n3=iorbfe(nv(3))
                     ns3=ispinfe(nv(3))
                     n4=iorbfe(nw(1))
                     ns4=ispinfe(nw(1))
                     n5=iorbfe(nw(2))
                     ns5=ispinfe(nw(2))
                     n6=iorbfe(nw(3))
                     ns6=ispinfe(nw(3))
                     n1=n1-nisht
                     n2=n2-nisht
                     n3=n3-nisht
                     n4=n4-nisht
                     n5=n5-nisht
                     n6=n6-nisht
                     if(flip(m).eq.1.and.n.ne.m-1)then
                        valcm=2.d0*c(m)
                     else
                        valcm=c(m)
                     endif
                     cij=valcm*c(n)*segno
                     val=cij
c     cij=c(m)*c(n)*segno
c     casi aaa, bbb, aab, bba
                     if(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.1)then
                        daaa(n1,n2,n3,n4,n5,n6)=daaa(n1,n2,n3,n4,n5,n6)
     $                       +val
                        daaa(n4,n5,n6,n1,n2,n3)=daaa(n4,n5,n6,n1,n2,n3)
     $                       +val
                        call perm3(daaa,n1,n2,n3,n4,n5,n6,nasht)
                        call perm3(daaa,n4,n5,n6,n1,n2,n3,nasht)
                     elseif(ns1.eq.0.and.ns2.eq.0.and.ns3.eq.0)then
                        dbbb(n1,n2,n3,n4,n5,n6)=dbbb(n1,n2,n3,n4,n5,n6)
     $                       +val
                        dbbb(n4,n5,n6,n1,n2,n3)=dbbb(n4,n5,n6,n1,n2,n3)
     $                       +val
                        call perm3(dbbb,n1,n2,n3,n4,n5,n6,nasht)
                        call perm3(dbbb,n4,n5,n6,n1,n2,n3,nasht)
                     elseif(ns1.eq.1.and.ns2.eq.1.and.ns3.eq.0)then
                        daab(n1,n2,n3,n4,n5,n6)=daab(n1,n2,n3,n4,n5,n6)
     $                       +val
                        call perm2(daab,n1,n2,n3,n4,n5,n6,nasht)
                     elseif(ns1.eq.1.and.ns2.eq.0.and.ns3.eq.0)then
                        dbba(n2,n3,n1,n5,n6,n4)=dbba(n2,n3,n1,n5,n6,n4)
     $                       +val
                        call perm2(dbba,n2,n3,n1,n5,n6,n4,nasht)
                     endif
c---- renzo            if(abs(c(n)).lt.thr)goto 2
                  endif
 4          continue
            enddo
 3          continue
            enddo
 2       continue
         enddo
 1    continue
      enddo
c--costruzione matrice ro3 spinless
      do i=1,nasht
         do j=1,nasht
            do k=1,nasht
               do ip=1,nasht
                  do jp=1,nasht
                     do kp=1,nasht
                        daaa(i,j,k,ip,jp,kp)=daaa(i,j,k,ip,jp,kp)+
     $                       dbbb(i,j,k,ip,jp,kp)+daab(i,j,k,ip,jp,kp)
     $                       +dbba(i,j,k,ip,jp,kp)+daab(i,k,j,ip,kp,jp)
     $                       +daab(j,k,i,jp,kp,ip)+dbba(i,k,j,ip,kp,jp)
     $                       +dbba(k,j,i,kp,jp,ip)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      deallocate(dbbb)
      deallocate(daab)
      deallocate(dbba)
      return
      end
c-------------------------------------------------------
      subroutine ro2(c,d2alal,nisht,nasht,thr,trou,part,ne,nd,
     $               ispinfe,iorbfe,itsym,iocca,ioccb,flip)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1 (Z)
#include <priunit.h>
      dimension d2alal(nasht,nasht,nasht,nasht)
      allocatable d2bebe(:,:,:,:),d2albe(:,:,:,:)
      dimension c(*)
      COMMON /CIP/NORB,NOCB,NCF,ZSEG
      INTEGER*2 ispinfe,iorbfe
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispinfe(*),iorbfe(*),itsym(*)
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      logical z1,z2,z3,z4,z5
      integer segno
      integer*1 flip
      dimension iocca(*),ioccb(*),flip(*)
      external index
      parameter (pt5=0.5d0)
c     indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c common e dimensionamenti mettere dopo
      allocate(d2bebe(1:nasht,1:nasht,1:nasht,1:nasht))
      allocate(d2albe(1:nasht,1:nasht,1:nasht,1:nasht))
      n=nasht
      nocc=NOCB-nisht
      nbig=(n*n+n)/2
      nstop=nbig
      nbig=(nbig*nbig+nbig)/2
      do i=1,n
         do j=1,n
            do k=1,n
               do l=1,n
                  d2alal(i,j,k,l)=0.d0
                  d2bebe(i,j,k,l)=0.d0
                  d2albe(i,j,k,l)=0.d0
               enddo
            enddo
         enddo
      enddo
      write(lupri,*)'subr. ro2'
      write(lupri,*)'n=',n,' nocc=',nocc,' ncf=',ncf
      do 1 m=1,ncf
         if (mod(m,100).eq.0)print *,'m=',m
         call flshfo(lupri)
         if(abs(c(m)).lt.thr)goto 1
         ndi=nd(m)
         nec=ne(m)
         call giveocc2(m,iocca,ioccb,nasht,nisht,nocca,noccb,nocc,nec
     $        ,ndi,trou,part,iorbfe,ispinfe)
         do 2 n=1,m
            if(abs(c(n)).lt.thr)goto 2
            if(n.eq.m)goto 3
            ndj=nd(n)
            ij=jdifgen(ne,nd,trou,part,m,n,nv,nw,2,zseg)
            if(ij.eq.0)nbdif=4
            if(nbdif.gt.2)goto 2
               if(zseg)then
                  segno=-1
               else
                  segno=1
               endif
               if(nbdif.eq.1)then
                  n1=iorbfe(nv(1))
                  ns1=ispinfe(nv(1))
                  n2=iorbfe(nw(1))
                  ns2=ispinfe(nw(1))
                  goto 10
               else
                  call ord2(nv,nw,zseg,segno) !prima gli alfa, poi i beta
                  n1=iorbfe(nv(1))
                  ns1=ispinfe(nv(1))
                  n2=iorbfe(nv(2))
                  ns2=ispinfe(nv(2))
                  n3=iorbfe(nw(1))
                  ns3=ispinfe(nw(1))
                  n4=iorbfe(nw(2))
                  ns4=ispinfe(nw(2))
c                  print '(a,4i3,l2)','due diff: ',n1,n3,n2,n4,zseg
                  goto 20
               endif
c            goto(10,20),nbdif
 3          continue
c     nessuna differenza
            ci2=c(m)*c(m)
            do i=1,nocca
               ia=iocca(i)
               do j=i+1,nocca
                  ja=iocca(j)
                  d2alal(ia,ja,ia,ja)=d2alal(ia,ja,ia,ja)+ci2
                  call permbis(d2alal,ia,ja,ia,ja,nasht)
               enddo
            enddo
            do i=1,noccb
               ib=ioccb(i)
               do j=i+1,noccb
                  jb=ioccb(j)
                  d2bebe(ib,jb,ib,jb)=d2bebe(ib,jb,ib,jb)+ci2
                  call permbis(d2bebe,ib,jb,ib,jb,nasht)
               enddo
            enddo
            do i=1,nocca
               ia=iocca(i)
               do j=1,noccb
                  jb=ioccb(j)
                  d2albe(ia,jb,ia,jb)=d2albe(ia,jb,ia,jb)+ci2
               enddo
            enddo
               goto 2
 10            continue
c     una differenza
               n1=n1-nisht
               n2=n2-nisht
               cij=c(m)*c(n)*segno
               if(ns1.eq.1.and.ns2.eq.1)then
                  do i=1,nocca
                     ia=iocca(i)
                     if(ia.ne.n1.and.ia.ne.n2)then
                        d2alal(n1,ia,n2,ia)=d2alal(n1,ia,n2,ia)+cij
                        call permbis(d2alal,n1,ia,n2,ia,nasht)
                     endif
                  enddo
                  do i=1,noccb
                     ib=ioccb(i)
                     d2albe(n1,ib,n2,ib)=d2albe(n1,ib,n2,ib)+cij
                     d2albe(n2,ib,n1,ib)=d2albe(n2,ib,n1,ib)+cij
                  enddo
               else
                  do i=1,noccb
                     ib=ioccb(i)
                     if(ib.ne.n1.and.ib.ne.n2)then
                        d2bebe(n1,ib,n2,ib)=d2bebe(n1,ib,n2,ib)+cij
                        call permbis(d2bebe,n1,ib,n2,ib,nasht)
                     endif
                  enddo
                  do i=1,nocca
                     ia=iocca(i)
                     d2albe(ia,n1,ia,n2)=d2albe(ia,n1,ia,n2)+cij
                     d2albe(ia,n2,ia,n1)=d2albe(ia,n2,ia,n1)+cij
                  enddo
               endif
               goto 2
c     due differenze
 20            continue
               n1=n1-nisht
               n2=n2-nisht
               n3=n3-nisht
               n4=n4-nisht
               cij=c(m)*c(n)*segno
c               print '(a,i1,a,i1,a,i1,a,i1)','Due diff.: ns1=',ns1,
c     $              ' ns2=',ns2,' ns3=',ns3,' ns4=',ns4
               if(ns1.eq.1.and.ns2.eq.1)then
                  d2alal(n1,n2,n3,n4)=d2alal(n1,n2,n3,n4)+cij
                  call permbis(d2alal,n1,n2,n3,n4,nasht)
               elseif(ns1.eq.0.and.ns2.eq.0)then
                  d2bebe(n1,n2,n3,n4)=d2bebe(n1,n2,n3,n4)+cij
                  call permbis(d2bebe,n1,n2,n3,n4,nasht)
               elseif(ns1.eq.1.and.ns2.eq.0)then
                  d2albe(n1,n2,n3,n4)=d2albe(n1,n2,n3,n4)+cij
                  d2albe(n3,n4,n1,n2)=d2albe(n3,n4,n1,n2)+cij
               endif
 2          continue
 1       continue
c--costruzione ro2 spinless
         do i=1,nasht
            do j=1,nasht
               do ip=1,nasht
                  do jp=1,nasht
                     d2alal(i,j,ip,jp)=d2alal(i,j,ip,jp)+
     $                    d2albe(i,j,ip,jp)+d2albe(j,i,jp,ip)
     $                    +d2bebe(i,j,ip,jp)
                  enddo
               enddo
            enddo
         enddo
         deallocate(d2bebe)
         deallocate(d2albe)
         return
         end
c-----------------------------------------------------
      subroutine ro1(c,d1,nisht,nasht,thr,trou,part,ne,nd,ispinfe,iorbfe
     $     ,itsym)

      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1 (Z)

#include <priunit.h>
      dimension d1(nasht,nasht)
      integer segno
      dimension c(*)
      COMMON /CIP/NORB,NOCB,NCF,ZSEG
      INTEGER*2 ispinfe,iorbfe
      INTEGER*2 NE,TROU,PART
      dimension ne(*),nd(*),ispinfe(*),iorbfe(*),itsym(*)
      dimension trou(*),part(*)
      common /diff/ n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      dimension nv(4),nw(4)
      n=nasht
      nocc=NOCB-nisht
C     note: NOCB = nisht + lfermie
      nn=(n*n+n)/2
      write(lupri,*)'subr. ro1'
      write(lupri,*)'n=',n,' nocc=',nocc,' ncf=',ncf
      nac=nocc
      do j=1,nac
         d1(j,j)=2.d0
      enddo
      do 1 i=1,ncf
         if(abs(c(i)).lt.thr)goto 1
         do 2 j=1,i
            if(abs(c(j)).lt.thr)goto 2
            if(i.eq.j)goto 3
            ij=jdifgen(ne,nd,trou,part,i,j,nv,nw,1,zseg)
            if(ij.eq.0)goto 2
            if(zseg)then
               segno=-1
            else
               segno=1
            endif
            n1=iorbfe(nv(1))
            ns1=ispinfe(nv(1))
            n2=iorbfe(nw(1))
            ns2=ispinfe(nw(1))
            n2=n2-nisht
            n1=n1-nisht
            d1(n1,n2)=d1(n1,n2)+c(i)*c(j)*segno
            goto 2
 3          nec=ne(i)
            if(nec.eq.0)goto 2
            cmi2=c(i)*c(i)
            ndi=nd(i)
            do 30 l=1,nec
               lt=trou(l+ndi)
               lp=part(l+ndi)
               lt=iorbfe(lt)-nisht
               lp=iorbfe(lp)-nisht
               d1(lt,lt)=d1(lt,lt)-cmi2
               d1(lp,lp)=d1(lp,lp)+cmi2
 30         continue
 2       continue
 1    continue
      return
      end
c-------------------------------------------------
      subroutine ord2(nv,nw,zseg,segno)
c--prima gli alfa poi i beta
#include <priunit.h>
      logical*1 zseg
      integer segno
      dimension nv(*),nw(*)
      if(nv(2).gt.nv(1))then
         n=nv(1)
         nv(1)=nv(2)
         nv(2)=n
         zseg=.not.zseg
      endif
      if(nw(2).gt.nw(1))then
         n=nw(1)
         nw(1)=nw(2)
         nw(2)=n
         zseg=.not.zseg
      endif
      if(zseg)then
         segno=-1
      else
         segno=1
      endif
      return
      end
c----------------------------------------------------------------
      subroutine ord3(nv,nw,zseg,segno)
c--prima gli alfa poi i beta
#include <priunit.h>
      logical*1 zseg
      integer segno
      dimension nv(*),nw(*)
      if(nv(2).gt.nv(1))then
         n=nv(1)
         nv(1)=nv(2)
         nv(2)=n
         zseg=.not.zseg
      endif
      if(nv(3).gt.nv(1))then
         n=nv(1)
         nv(1)=nv(3)
         nv(3)=n
         zseg=.not.zseg
      endif
      if(nv(3).gt.nv(2))then
         n=nv(2)
         nv(2)=nv(3)
         nv(3)=n
         zseg=.not.zseg
      endif
      if(nw(2).gt.nw(1))then
         n=nw(1)
         nw(1)=nw(2)
         nw(2)=n
         zseg=.not.zseg
      endif
      if(nw(3).gt.nw(1))then
         n=nw(1)
         nw(1)=nw(3)
         nw(3)=n
         zseg=.not.zseg
      endif
      if(nw(3).gt.nw(2))then
         n=nw(2)
         nw(2)=nw(3)
         nw(3)=n
         zseg=.not.zseg
      endif
      if(zseg)then
         segno=-1
      else
         segno=1
      endif
      return
      end
c----------------------------------------------------------------
      subroutine ord4(nv,nw,zseg,segno)
c--prima gli alfa poi i beta
#include <priunit.h>
      logical*1 zseg
      integer segno
      dimension nv(*),nw(*)
      do i=1,4
         imax=nv(i)
         do j=i+1,4
            if(nv(j).gt.imax)then
               n=imax
               imax=nv(j)
               nv(j)=n
               zseg=.not.zseg
            endif
            nv(i)=imax
         enddo
      enddo
      do i=1,4
         imax=nw(i)
         do j=i+1,4
            if(nw(j).gt.imax)then
               n=imax
               imax=nw(j)
               nw(j)=n
               zseg=.not.zseg
            endif
            nw(i)=imax
         enddo
      enddo
      if(zseg)then
         segno=-1
      else
         segno=1
      endif
      return
      end
c----------------------------------------------------------------
      subroutine giveocc2(m,iocca,ioccb,nasht,nisht,nocca,noccb,nocc,nec
     $     ,ndm,trou,part,iorbfe,ispinfe) 
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1 (Z)
#include <priunit.h>
      dimension iocca(*),ioccb(*),trou(*),part(*),iorbfe(*),ispinfe(*)
      integer*2 trou,part,iorbfe,ispinfe
      do i=1,nasht
         if(i.le.nocc)then
            iocca(i)=i
            ioccb(i)=i
         else
            iocca(i)=0
            ioccb(i)=0
         endif
      enddo
      do i=1,nec
         it=trou(ndm+i)
         ip=part(ndm+i)
         its=ispinfe(it)
         ips=ispinfe(ip)
         it=iorbfe(it)-nisht
         ip=iorbfe(ip)-nisht
         if(its.eq.1)then
            iocca(it)=0
         else
            ioccb(it)=0
         endif
         if(ips.eq.1)then
            iocca(ip)=ip
         else
            ioccb(ip)=ip
         endif
      enddo
      nocca=0
      noccb=0
      do i=1,nasht
         if(iocca(i).ne.0)then
            nocca=nocca+1
            iocca(nocca)=iocca(i)
         endif
         if(ioccb(i).ne.0)then
            noccb=noccb+1
            ioccb(noccb)=ioccb(i)
         endif
      enddo
      return
      end
c--------------------------------------------------------------------
      subroutine perm3(d,i,j,k,l,m,n,nasht)

      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1 (Z)
#include <priunit.h>
      dimension d(nasht,nasht,nasht,nasht,nasht,nasht)
c--block 1
      a=d(i,j,k,l,m,n)
      d(j,i,k,l,m,n)=-a
      d(i,k,j,l,m,n)=-a
      d(k,j,i,l,m,n)=-a
      d(k,i,j,l,m,n)=a
      d(j,k,i,l,m,n)=a
c--block 2
      d(i,j,k,m,l,n)=-a
      d(j,i,k,m,l,n)=a
      d(i,k,j,m,l,n)=a
      d(k,j,i,m,l,n)=a
      d(k,i,j,m,l,n)=-a
      d(j,k,i,m,l,n)=-a
c--block 3
      d(i,j,k,l,n,m)=-a
      d(j,i,k,l,n,m)=a
      d(i,k,j,l,n,m)=a
      d(k,j,i,l,n,m)=a
      d(k,i,j,l,n,m)=-a
      d(j,k,i,l,n,m)=-a
c--block 4
      d(i,j,k,n,m,l)=-a
      d(j,i,k,n,m,l)=a
      d(i,k,j,n,m,l)=a
      d(k,j,i,n,m,l)=a
      d(k,i,j,n,m,l)=-a
      d(j,k,i,n,m,l)=-a
c--block 5
      d(i,j,k,n,l,m)=a
      d(j,i,k,n,l,m)=-a
      d(i,k,j,n,l,m)=-a
      d(k,j,i,n,l,m)=-a
      d(k,i,j,n,l,m)=a
      d(j,k,i,n,l,m)=a
c--block 6
      d(i,j,k,m,n,l)=a
      d(j,i,k,m,n,l)=-a
      d(i,k,j,m,n,l)=-a
      d(k,j,i,m,n,l)=-a
      d(k,i,j,m,n,l)=a
      d(j,k,i,m,n,l)=a
      return
      end
c-----------------------------------------------------------
      subroutine perm2(d,i,j,k,l,m,n,nasht)
#include <implicit.h>
#include <priunit.h>
      dimension d(nasht,nasht,nasht,nasht,nasht,nasht)
c--block 1
      a=d(i,j,k,l,m,n)
      d(j,i,k,l,m,n)=-a
      d(l,m,n,i,j,k)=a
      d(l,m,n,j,i,k)=-a
c--block 2
      d(i,j,k,m,l,n)=-a
      d(j,i,k,m,l,n)=a
      d(m,l,n,i,j,k)=-a
      d(m,l,n,j,i,k)=a
      return
      end
c--------------------------------------------------------------------
      subroutine permbis(d,i,j,k,l,nasht)
#include <implicit.h>
#include <priunit.h>
      dimension d(nasht,nasht,nasht,nasht)
c--block 1
      a=d(i,j,k,l)
      d(j,i,k,l)=-a
      d(k,l,i,j)=a
      d(k,l,j,i)=-a
c--block 2
      d(i,j,l,k)=-a
      d(j,i,l,k)=a
      d(l,k,i,j)=-a
      d(l,k,j,i)=a
      return
      end
c--------------------------------------------------------------------
      subroutine permut123(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
#include <implicit.h>
#include <priunit.h>
      integer segno
      lp=l
      if(ivolte.eq.1)then
c---permutazioni pari
         ip=i
         jp=j
         kp=k
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=k
         kp=i
         segno=1
      elseif(ivolte.eq.3)then
         ip=k
         jp=i
         kp=j
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.4)then
         ip=j
         jp=i
         kp=k
         segno=-1
      elseif(ivolte.eq.5)then
         ip=k
         jp=j
         kp=i
         segno=-1
      elseif(ivolte.eq.6)then
         ip=i
         jp=k
         kp=j
         segno=-1
      endif
      return
      end
c----------------------------------------------------------
      subroutine permut124(i,j,l,k,ip,jp,lp,kp,segno,ivolte)
#include <implicit.h>
#include <priunit.h>
      integer segno
      lp=l
      if(ivolte.eq.1)then
c---permutazioni pari
         ip=i
         jp=j
         kp=k
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=k
         kp=i
         segno=1
      elseif(ivolte.eq.3)then
         ip=k
         jp=i
         kp=j
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.4)then
         ip=j
         jp=i
         kp=k
         segno=-1
      elseif(ivolte.eq.5)then
         ip=k
         jp=j
         kp=i
         segno=-1
      elseif(ivolte.eq.6)then
         ip=i
         jp=k
         kp=j
         segno=-1
      endif
      return
      end
c--------------------------------------------------------------------
      subroutine permut134(i,l,j,k,ip,lp,jp,kp,segno,ivolte)
#include <implicit.h>
#include <priunit.h>
      integer segno
      lp=l
      if(ivolte.eq.1)then
c---permutazioni pari
         ip=i
         jp=j
         kp=k
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=k
         kp=i
         segno=1
      elseif(ivolte.eq.3)then
         ip=k
         jp=i
         kp=j
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.4)then
         ip=j
         jp=i
         kp=k
         segno=-1
      elseif(ivolte.eq.5)then
         ip=k
         jp=j
         kp=i
         segno=-1
      elseif(ivolte.eq.6)then
         ip=i
         jp=k
         kp=j
         segno=-1
      endif
      return
      end
c--------------------------------------------------------------------
      subroutine permut234(l,i,j,k,lp,ip,jp,kp,segno,ivolte)
#include <implicit.h>
#include <priunit.h>
      integer segno
      lp=l
      if(ivolte.eq.1)then
c---permutazioni pari
         ip=i
         jp=j
         kp=k
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=k
         kp=i
         segno=1
      elseif(ivolte.eq.3)then
         ip=k
         jp=i
         kp=j
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.4)then
         ip=j
         jp=i
         kp=k
         segno=-1
      elseif(ivolte.eq.5)then
         ip=k
         jp=j
         kp=i
         segno=-1
      elseif(ivolte.eq.6)then
         ip=i
         jp=k
         kp=j
         segno=-1
      endif
      return
      end
c--------------------------------------------------------------------
      subroutine permut1234(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
#include <implicit.h>
#include <priunit.h>
      integer segno
c--permutazioni pari
      if(ivolte.eq.1)then
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
         ip=j
         jp=i
         kp=l
         lp=k
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.3)then
         ip=j
         jp=i
         kp=k
         lp=l
         segno=-1
      elseif(ivolte.eq.4)then
         ip=i
         jp=j
         kp=l
         lp=k
         segno=-1
      endif
      return
      end
c--------------------------------------------------------------------
      subroutine permut1324(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
#include <implicit.h>
#include <priunit.h>
      integer segno
c--permutazioni pari
      if(ivolte.eq.1)then
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
         ip=k
         jp=l
         kp=i
         lp=j
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.3)then
         ip=k
         jp=j
         kp=i
         lp=l
         segno=-1
      elseif(ivolte.eq.4)then
         ip=i
         jp=l
         kp=k
         lp=j
         segno=-1
      endif
      return
      end
c--------------------------------------------------------------------
      subroutine permut1423(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
#include <implicit.h>
#include <priunit.h>
      integer segno
c--permutazioni pari
      if(ivolte.eq.1)then
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
         ip=l
         jp=k
         kp=j
         lp=i
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.3)then
         ip=l
         jp=j
         kp=k
         lp=i
         segno=-1
      elseif(ivolte.eq.4)then
         ip=i
         jp=k
         kp=j
         lp=l
         segno=-1
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine permut4(i,j,k,l,ip,jp,kp,lp,segno,ivolte)
#include <implicit.h>
#include <priunit.h>
      integer segno
      if(ivolte.eq.1)then
c--permutazioni pari      
         ip=i
         jp=j
         kp=k
         lp=l
         segno=1
      elseif(ivolte.eq.2)then
         ip=k
         jp=i
         kp=j
         lp=l
         segno=1
      elseif(ivolte.eq.3)then
         ip=l
         jp=i
         kp=k
         lp=j
         segno=1
      elseif(ivolte.eq.4)then
         ip=j
         jp=k
         kp=i
         lp=l
         segno=1
      elseif(ivolte.eq.5)then
         ip=j
         jp=l
         kp=k
         lp=i
         segno=1
      elseif(ivolte.eq.6)then
         ip=j
         jp=i
         kp=l
         lp=k
         segno=1
      elseif(ivolte.eq.7)then
         ip=l
         jp=j
         kp=i
         lp=k
         segno=1
      elseif(ivolte.eq.8)then
         ip=k
         jp=l
         kp=i
         lp=j
         segno=1
      elseif(ivolte.eq.9)then
         ip=k
         jp=j
         kp=l
         lp=i
         segno=1
      elseif(ivolte.eq.10)then
         ip=l
         jp=k
         kp=j
         lp=i
         segno=1
      elseif(ivolte.eq.11)then
         ip=i
         jp=l
         kp=j
         lp=k
         segno=1
      elseif(ivolte.eq.12)then
         ip=i
         jp=k
         kp=l
         lp=j
         segno=1
c---permutazioni dispari
      elseif(ivolte.eq.13)then
         ip=j
         jp=i
         kp=k
         lp=l
         segno=-1
      elseif(ivolte.eq.14)then
         ip=k
         jp=j
         kp=i
         lp=l
         segno=-1
      elseif(ivolte.eq.15)then
         ip=l
         jp=j
         kp=k
         lp=i
         segno=-1
      elseif(ivolte.eq.16)then
         ip=i
         jp=k
         kp=j
         lp=l
         segno=-1
      elseif(ivolte.eq.17)then
         ip=i
         jp=l
         kp=k
         lp=j
         segno=-1
      elseif(ivolte.eq.18)then
         ip=i
         jp=j
         kp=l
         lp=k
         segno=-1
      elseif(ivolte.eq.19)then
         ip=l
         jp=i
         kp=j
         lp=k
         segno=-1
      elseif(ivolte.eq.20)then
         ip=k
         jp=l
         kp=j
         lp=i
         segno=-1
      elseif(ivolte.eq.21)then
         ip=k
         jp=i
         kp=l
         lp=j
         segno=-1
      elseif(ivolte.eq.22)then
         ip=l
         jp=k
         kp=i
         lp=j
         segno=-1
      elseif(ivolte.eq.23)then
         ip=j
         jp=l
         kp=i
         lp=k
         segno=-1
      elseif(ivolte.eq.24)then
         ip=j
         jp=k
         kp=l
         lp=i
         segno=-1
      endif
      return
      end
c---------------------------------------------------------------------
      FUNCTION Jdifgen(ne,nd,trou,part,II,JJ,nv1,nv2,nma,zsig)
C     
C     SSP DE RECHERCHE DES ELEMENTS D'INTERACTION
C     ENTRE DEUX DETERMINANTS, NUMEROTES II ET JJ
C     
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1 (Z)
#include <priunit.h>
      parameter (nmax=4)
      dimension ne(*),nd(*),trou(*),part(*)
      INTEGER*2 NE,TROU,PART
      DIMENSION ITA(99),ITB(99),NV1(nmax),NV2(nmax),NATU1(nmax)
     $     ,NATU2(nmax)
      common/diff/n1,n2,n3,n4,ns1,ns2,ns3,ns4,nbdif
      KR=II
      KL=JJ
C     
C     LE DETERMINANT - 1 - EST CHOISI COMME LE
C     DETERMINANT AYANT LE PLUS D'O.M. EXCITEES
C     
      IF (NE(KR) .lt. NE(KL)) then
        KKKK=KR
        KR=KL
        KL=KKKK
        zshift=.true.
      else
        zshift=.false.
      end if
      NE1=NE(KR)
C     TEST SUR LE FONDAMENTAL
      IF (NE1.EQ.0) GOTO 44
      NE2=NE(KL)
C     NBDIF : EST LE NOMBRE DE SPIN-ORBITALES
C     DIFFERENTES
      NBDIF=NE1-NE2
      IF (NBDIF.GT.nma) GOTO 24
C     CONSTRUCTION POUR LE DETERMINANT -1- ET -2-
C     DE IT ET NS, QUI CONTIENNENT LES NUMEROS
C     DES ORBITALES ET LEURS SPINS
C     
      ND1=ND(KR)
      ND2=ND(KL)
      DO 15 J=1,NE1
         ITA(J)=TROU(J+ND1)
         ITA(J+NE1)=PART(J+ND1)
 15   CONTINUE
C     
C     ZSIG :  EST LA PARITE DU  NOMBRE DE CROISEMENTS
C     DANS LE DIAGRAMME D INTERACTION
      ZSIG=.FALSE.
C     SI NE2=0 LE DETERMINANT -2- EST LE FONDAMENTAL
      IF (NE2.EQ.0) GOTO 125
 121  CONTINUE
      DO 25 J= 1,NE2
         ITB(J)=TROU(J+ND2)
 25   ITB(J+NE2)=PART(J+ND2)
    4 DO 8 K=1,2
         K1=(K-1)*NE1
         K2=(K-1)*NE2
         DO 10I=1,NE2
            NIT1=ITB(K2+I)
            J1=1+K1
            J2=NE1+K1
            DO 12 J=J1,J2
               NJT2=ITA(J)
               IF (NIT1.NE.NJT2) GOTO 12
               IF (I.EQ.(J-K1)) GOTO 10
               ITA(J)=ITA(I+K1)
               ITA(I+K1)=NJT2
               ZSIG=.NOT.ZSIG
               GOTO 10
 12         CONTINUE
            NBDIF=NBDIF+1
            IF (NBDIF.GT.nma) GOTO 24
            NV1(NBDIF)=NIT1
            NATU1(NBDIF)=K
            NATU2(NBDIF)=I
 10         CONTINUE
 8       CONTINUE
 125     IF (NBDIF.LE.0) GOTO 44
C     
C     NAR : EXCITATION SUPLEMENTAIRE DU DETERMINANT -1-
C     PAR RAPPORT AU DETERMINANT -2-
 26      NAR=NE1-NE2
         IF (NAR.LE.0) GOTO 28
 30      DO 32 I=1,NAR
            NET=NE1+1-I
            NV1(I)=ITA(NET)
            NV2(I)=ITA(NET+NE1)
 32      CONTINUE
 28      NBAR=NBDIF-NAR
         IF (NBAR.LE.0) GOTO 34
 36      NAR1=NAR+1
         DO 38 I=NAR1,NBDIF
            K=NATU1(I)
            NI=NATU2(I)
            NU=NV1(I)
            NV2(I)=ITA(NI+NE1*(K-1))
            IF (K.EQ.2) GOTO 38
 40         NV1(I)=NV2(I)
            NV2(I)=NU
            ZSIG=.NOT.ZSIG
 38      CONTINUE
 34      CONTINUE
C     CALCUL DE L ELEMENT DE MATRICE
 2       continue
         if(zshift)then
            do i=1,nbdif
               ndum=nv1(i)
               nv1(i)=nv2(i)
               nv2(i)=ndum
            enddo
         endif
         jdifgen=1
         RETURN
 24      Jdifgen=0
         RETURN
 44      WRITE (lupri,1006) II,JJ,NE1,NE2
 1006    FORMAT (//5X,'IDENTITE DANS HNTD'/' II,JJ,NE1,NE2 =',4I4)
         WRITE (lupri,1007) (ITA(K),K=1,NE1)
         WRITE (lupri,1008) (ITA(K+NE1),K=1,NE1)
         WRITE (lupri,1007) (ITB(K),K=1,NE2)
         WRITE (lupri,1008) (ITB(K+NE2),K=1,NE2)
 1007    FORMAT (' TROU',8I4)
 1008    FORMAT (' PART',8I4)
         STOP
         END
c------------------------------------------------------------------
      integer function index(i,j,k,l)
#include <implicit.h>
      if(i.ge.k)then
         ik=(i*i-i)/2+k
      else
         ik=(k*k-k)/2+i
      endif
      if(j.ge.l)then
         jl=(j*j-j)/2+l
      else
         jl=(l*l-l)/2+j
      endif
      if(ik.ge.jl)then
         index=(ik*ik-ik)/2+jl
      else
         index=(jl*jl-jl)/2+ik
      endif
      return
      end
c------------------------------------------------------------------
      subroutine accum(amat,ip,jp,kp,lp,mp,np,op,pp,val,nasht)
#include <implicit.h>
#include <priunit.h>
#include "nevpt2.h"
      dimension amat(nwords,nasht,nasht,nasht,nasht)
      integer op,pp
      lindice(i,j,k,l)=lindi(i)+lindj(j)+lindk(k)+l

      ind=lindice(ip,jp,kp,lp)
      amat(ind,mp,np,op,pp)=amat(ind,mp,np,op,pp)+val
      return
      end
c-----------------------------------------------------------------
      SUBROUTINE MATOUT(N,N2,A,NN)
#include <implicit.h>
#include <priunit.h>
c     CHARACTER*3 SYMB(NN)
c     integer symb(*)
      DIMENSION A(NN,*)
C 
   20 FORMAT(//,8X,11(5X,I3,3X))
   30 FORMAT(1X,I4,4X,11(F11.6))
   35 FORMAT(8X,11(5X,I3,3X))
   37 FORMAT(8X,11(4X,A4,3X))
   36 FORMAT(/)
C 
c      DO 40 M=1,N2,11
      DO 40 M=1,N2,6
c      K=M+10
      K=M+5
      IF(K.LE.N2) GOTO 10
      K=N2
   10 WRITE(lupri,20) (J,J=M,K)
      WRITE(lupri,36)
c     WRITE(lupri,35) (ME(I),I=M,K)
c     WRITE(lupri,37) (SYMB(I),I=M,K)
      WRITE(lupri,36)
      DO 40 I=1,N
      WRITE(lupri,30) I,(A(I,J),J=M,K)
 40   CONTINUE
      call flshfo(lupri)
      RETURN
      END
c-------------------------------------------------------------------
      subroutine koop2E(atwo,daaa,taa,dal,koopaa,koopeaa,nasht
     $     ,nisht,f,aj,itsym,its)
#include <implicit.h>
#include <priunit.h>
c---F e' la matrice dei monoelettronici!
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht),
     $ taa(nasht,nasht,nasht,nasht),dal(nasht,nasht)
      dimension atwo(nasht,nasht,nasht,nasht)
      dimension koopaa(nasht,nasht,nasht,nasht),
     $ koopeaa(nasht,nasht,nasht,nasht)
      REAL*8  koopaa,koopbb,koopab,koopeaa,koopebb,koopeab
      
      dimension itsym(*),its(8,8)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,apc,bpc
      dimension w(100),work(1000)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c--building heff'
c--in AJ is heff
      do a=1,nasht
         ac=a+nisht
         do b=a,nasht
            bc=b+nisht
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nasht
               cc=c+nisht
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-dum*0.5d0
            enddo
         enddo
      enddo
c-------------------
      do a=1,nasht
         ac=a+nisht
	 isa=itsym(ac)
         do b=1,nasht
            bc=b+nisht
	    isb=itsym(bc)
	    isab=its(isa,isb)
            do ap=1,nasht
               apc=ap+nisht
	       isd=itsym(apc)
               do bp=1,nasht
                  bpc=bp+nisht
		  isc=itsym(bpc)
		  iscd=its(isc,isd)
                  koopaa(ap,bp,a,b)=0.d0
		  if (its(isab,iscd).ne.1) goto 1
                  do d=1,nasht
                     dc=d+nisht
                     jad=indice(ac,dc)
                     jbd=indice(bc,dc)
                     koopaa(ap,bp,a,b)=koopaa(ap,bp,a,b)-f(jad)*taa(ap
     $                    ,bp,d,b)-f(jbd)*taa(ap,bp,a,d)
                     do ie=1,nasht
                        do ief=1,nasht
                           dint=atwo(d,ie,a,ief)
                           dum=daaa(ap,bp,d,ief,b,ie)
                           if(d.eq.ief)dum=dum+0.5d0*taa(ap,bp,ie,b)
                           if(d.eq.b)dum=dum+0.5d0*taa(ap,bp,ief,ie)
                           koopaa(ap,bp,a,b)=koopaa(ap,bp,a,b)-dint*dum
                           dint=atwo(d,ie,b,ief)
                           dum=daaa(ap,bp,d,a,ief,ie)
                           if(d.eq.a)dum=dum+0.5d0*taa(ap,bp,ie,ief)
                           if(d.eq.ief)dum=dum+0.5d0*taa(ap,bp,a,ie)
                           koopaa(ap,bp,a,b)=koopaa(ap,bp,a,b)-dint*dum
                        enddo
                     enddo
                  enddo
 1             continue
               enddo
            enddo
         enddo
      enddo
      do a=1,nasht
         ac=a+nisht
	 isa=itsym(ac)
         do b=1,nasht
            bc=b+nisht
	    isb=itsym(bc)
	    isab=its(isa,isb)
            do ap=1,nasht
               apc=ap+nisht
	       isd=itsym(apc)
               do bp=1,nasht
                  bpc=bp+nisht
		  isc=itsym(bpc)
		  iscd=its(isc,isd)
                  koopeaa(ap,bp,a,b)=0.d0
		  if (its(isab,iscd).ne.1) goto 2
                  do c=1,nasht
                     cc=c+nisht
                     jac=indice(ac,cc)
                     jbc=indice(bc,cc)
                     dum1=ro2t(ap,bp,c,b,taa,dal,nasht)
                     dum2=ro2t(ap,bp,a,c,taa,dal,nasht)
                     koopeaa(ap,bp,a,b)=koopeaa(ap,bp,a,b)+f(jac)*dum1
     $                    +f(jbc)*dum2
                     do d=1,nasht
                        do ie=1,nasht
                           dint=atwo(c,ie,d,a)
                           dum=-ro3t(ap,bp,ie,d,b,c,daaa,taa,dal,nasht)
                           if(c.eq.ie)dum=dum+2.d0*ro2t(ap,bp,d,b,taa
     $                          ,dal,nasht)
                           if(d.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,c,b,taa
     $                          ,dal,nasht)
                           if(b.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,d,c,taa
     $                          ,dal,nasht)
                           koopeaa(ap,bp,a,b)=koopeaa(ap,bp,a,b)+dint
     $                          *dum
                           dint=atwo(c,ie,d,b)
                           dum=-ro3t(ap,bp,ie,a,d,c,daaa,taa,dal,nasht)
                           if(c.eq.ie)dum=dum+2.d0*ro2t(ap,bp,a,d,taa
     $                          ,dal,nasht)
                           if(a.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,c,d,taa
     $                          ,dal,nasht)
                           if(d.eq.ie)dum=dum-0.5d0*ro2t(ap,bp,a,c,taa
     $                          ,dal,nasht)
                           koopeaa(ap,bp,a,b)=koopeaa(ap,bp,a,b)+dint
     $                          *dum
                        enddo
                     enddo
                  enddo
 2             continue
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine koopman0pE(atwo,daaa,taa,dal,koopaa,koopbb,koopab
     $     ,nasht,nisht,f,aj,itsym,its)
#include <implicit.h>
#include <priunit.h>
c---F e' la matrice dei monoelettronici!
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht),taa(nasht
     $     ,nasht,nasht,nasht),atwo(nasht,nasht,nasht,nasht),dal(nasht
     $     ,nasht)
      dimension koopaa(nasht,nasht,nasht,nasht),koopbb(nasht,nasht,nasht
     $     ,nasht),koopab(nasht,nasht)
      REAL*8  koopaa,koopbb,koopab
      dimension itsym(*),its(8,8)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,apc,bpc
      dimension w(100),work(1000)
      parameter (two=2.d0)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      do a=1,nasht
         ac=a+nisht
         do b=a,nasht
            bc=b+nisht
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nasht
               cc=c+nisht
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-0.5d0*dum
            enddo
         enddo
      enddo
      do a=1,nasht
         ac=a+nisht
         do b=1,nasht
            bc=b+nisht
            do ap=1,nasht
               apc=ap+nisht
               do bp=1,nasht
                  bpc=bp+nisht
                  koopaa(bp,ap,a,b)=0.d0
                  do c=1,nasht
                     cc=c+nisht
                     ica=indice(cc,ac)
                     ibc=indice(bc,cc)
                     koopaa(bp,ap,a,b)=koopaa(bp,ap,a,b)+2.d0*(f(ica)
     $                    *ee2(bp,ap,c,b,taa,dal,nasht)-f(ibc)
     $                    *ee2(bp,ap,a,c,taa,dal,nasht))
                     do d=1,nasht
                        dc=d+nisht
                        do ie=1,nasht
                           iec=ie+nisht
                           dum1=eee2(bp,ap,c,ie,d,b,daaa
     $                          ,taa,dal,nasht)+eee2(bp,ap,d,b
     $                          ,c,ie,daaa,taa,dal,nasht)
                           dum2=eee2(bp,ap,a,ie,c,d,daaa,taa,dal,nasht)
     $                          +eee2(bp,ap,c,d,a,ie,daaa,taa,dal,nasht)
                           koopaa(bp,ap,a,b)=koopaa(bp,ap,a,b)+atwo(c
     $                          ,ie,d,a)*dum1-atwo(b,ie,c,d)*dum2
                        enddo
                     enddo
                  enddo
                  koopbb(bp,ap,a,b)=0.d0
                  do c=1,nasht
                     cc=c+nisht
                     ica=indice(cc,ac)
                     ibc=indice(bc,cc)
                     dum1=-ee2(a,ap,bp,c,taa,dal,nasht)
                     if(a.eq.ap)dum1=dum1+two*dal(bp,c)
                     if(bp.eq.ap)dum1=dum1+dal(a,c)
                     dum2=-ee2(c,ap,bp,b,taa,dal,nasht)
                     if(ap.eq.c)dum2=dum2+two*dal(bp,b)
                     if(bp.eq.ap)dum2=dum2+dal(c,b)
                     koopbb(bp,ap,a,b)=koopbb(bp,ap,a,b)-f(ibc)*dum1
     $                    +f(ica)*dum2
                     do d=1,nasht
                        dc=d+nisht
                        do ie=1,nasht
                           iec=ie+nisht
                           dum1=-eee2(c,ie,a,ap,bp,d,daaa,taa,dal,nasht)
     $                          -eee2(a,ap,bp,d,c,ie,daaa,taa,dal,nasht)
                           if(ap.eq.a)dum1=dum1+two*(ee2(c,ie,bp,d,taa
     $                          ,dal,nasht)+ee2(bp,d,c,ie,taa
     $                          ,dal,nasht))
                           if(bp.eq.ap)dum1=dum1+ee2(c,ie,a,d,taa,dal
     $                          ,nasht)+ee2(a,d,c,ie,taa,dal,nasht)
                           if(ap.eq.c)then
                              dum1=dum1-ee2(a,ie,bp,d,taa,dal,nasht)
                              if(ie.eq.a)dum1=dum1+two*dal(bp,d)
                           endif
                           if(bp.eq.ie)then
                              dum1=dum1+ee2(a,ap,c,d,taa,dal,nasht)
                              if(ap.eq.a)dum1=dum1-two*dal(c,d)
                           endif
                           dum2=-eee2(c,ie,d,ap,bp,b,daaa
     $                          ,taa,dal,nasht)-eee2(d,ap,bp,b
     $                          ,c,ie,daaa,taa,dal,nasht)
                           if(ap.eq.d)dum2=dum2+two*(ee2(c,ie,bp,b,taa
     $                          ,dal,nasht)+ee2(bp,b,c,ie,taa
     $                          ,dal,nasht))
                           if(bp.eq.ap)dum2=dum2+ee2(c,ie,d,b,taa,
     $                          dal,nasht)+ee2(d,b,c,ie,taa,dal,nasht)
                           if(ap.eq.c)then
                              dum2=dum2-ee2(d,ie,bp,b,taa,dal,nasht)
                              if(d.eq.ie)dum2=dum2+two*dal(bp,b)
                           endif
                           if(bp.eq.ie)then
                              dum2=dum2+ee2(d,ap,c,b,taa,dal,nasht)
                              if(d.eq.ap)dum2=dum2-two*dal(c,b)
                           endif
                           koopbb(bp,ap,a,b)=koopbb(bp,ap,a,b)-0.5d0*
     $                          atwo(c,ie,b,d)*dum1+0.5d0*atwo(c,ie,d
     $                          ,a)*dum2
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            koopab(a,b)=0.d0
            do c=1,nasht
               cc=c+nisht
               ica=indice(cc,ac)
               ibc=indice(bc,cc)
               dum1=dal(c,b)
               dum2=dal(a,c)
               koopab(a,b)=koopab(a,b)+two*(f(ica)*dum1
     $              -f(ibc)*dum2)
               do d=1,nasht
                  dc=d+nisht
                  do ie=1,nasht
                     iec=ie+nisht
                     dum1=ee2(c,ie,d,b,taa,dal,nasht)+ee2(d,b,c
     $                    ,ie,taa,dal,nasht)
                     dum2=ee2(a,ie,c,d,taa,dal,nasht)+ee2(c,d,a
     $                    ,ie,taa,dal,nasht)
                     koopab(a,b)=koopab(a,b)+atwo(c,ie,d
     $                    ,a)*dum1-atwo(b,ie,c,d)*dum2
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine bamat(amat,atmat,atwo,ro4,daaa,taa,dal,nasht,nisht,f,aj
     $     ,itsym,its,norb)
#include <implicit.h>
#include <priunit.h>
c---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht),taa(nasht
     $     ,nasht,nasht,nasht),dal(nasht,nasht)
#include "nevpt2.h"
      dimension ro4(nwords,nasht,nasht,nasht,nasht)
      dimension amat(nasht,nasht,nasht,nasht,nasht,nasht)
      dimension atmat(nasht,nasht,nasht,nasht,nasht,nasht)
      dimension atwo(nasht,nasht,nasht,nasht)
      dimension itsym(*),its(8,8)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      allocatable indice(:,:)
      allocate(indice(norb,norb))
c-------------------
      do i=1,norb
       do j=1,norb
       indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
       enddo
      enddo
c--building heff'
c--in AJ is heff
      do a=1,nasht
         ac=a+nisht
         do b=a,nasht
            bc=b+nisht
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nasht
               cc=c+nisht
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-dum*0.5d0
            enddo
         enddo
      enddo
c-------------------
      do a=1,nasht
         ac=a+nisht
	 isa=itsym(ac)
         do b=1,nasht
            bc=b+nisht
	    isb=itsym(bc)
	    isab=its(isa,isb)
            do c=1,nasht
               cc=c+nisht
	       isc=itsym(cc)
	       isabc=its(isab,isc)
               isbc=its(isb,isc)
               isac=its(isa,isc)
               do ap=1,nasht
                  apc=ap+nisht
		  isap=itsym(apc)
                  do bp=1,nasht
                     bpc=bp+nisht
		     isbp=itsym(bpc)
		     isabp=its(isap,isbp)
                     do cp=1,nasht
                        cpc=cp+nisht
			iscp=itsym(cpc)
			isabcp=its(isabp,iscp)
			isatot=its(isabc,isabcp)
			if (isatot.ne.1) then
			goto 1
			endif
                        do d=1,nasht
                           dc=d+nisht
			   isd=itsym(dc)
                           ida=indice(dc,ac)
                           icd=indice(cc,dc)
                           ibd=indice(bc,dc)
                           amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b
     $                          ,c)+f(ida)*eee2(cp,ap,bp,b,d,c,daaa,taa
     $                          ,dal,nasht)-f(icd)*eee2(cp,ap,bp,b,a,d
     $                          ,daaa,taa,dal,nasht)-f(ibd)*eee2(cp,ap
     $                          ,bp,d,a,c,daaa,taa,dal,nasht)
                           atmat(ap,bp,cp,a,b,c)=atmat(ap,bp,cp,a,b
     $                          ,c)+f(ida)*eee2t(cp,ap,bp,b,d,c,daaa,taa
     $                          ,dal,nasht)-f(icd)*eee2t(cp,ap,bp,b,a,d
     $                          ,daaa,taa,dal,nasht)+f(ibd)*eee2t(cp,ap
     $                          ,bp,d,a,c,daaa,taa,dal,nasht)
                           do ie=1,nasht
                              iec=ie+nisht
			      ise=itsym(iec)
                              dint1=atwo(c,ie,d,a)
                              dint2=atwo(d,ie,ie,a)
                              dint3=atwo(c,d,d,ie)
                              dint4=atwo(b,d,d,ie)
                              dum1=eee2(cp,ap,bp,b,d,ie,daaa,taa,dal
     $                             ,nasht)
                              dum2=eee2(cp,ap,bp,b,d,c,daaa,taa,dal
     $                             ,nasht)
                              dum3=eee2(cp,ap,bp,b,a,ie,daaa,taa,dal
     $                             ,nasht)
                              dum4=eee2(cp,ap,bp,ie,a,c,daaa,taa,dal
     $                             ,nasht)
                              dum1t=-dum1
                              dum2t=-dum2
                              dum3t=-dum3
                              if(bp.eq.b)then
                                 dum1t=dum1t+2.d0*ee2(cp,ap,d
     $                                ,ie,taa,dal,nasht)
                                 dum2t=dum2t+2.d0*ee2(cp,ap,d,c
     $                                ,taa,dal,nasht)
                                 dum3t=dum3t+2.d0*ee2(cp,ap,a,ie
     $                                ,taa,dal,nasht)
                              endif
                              atmat(ap,b,cp,a,bp,c)=atmat(ap,b,cp,a,bp,c
     $                             )+dum1t*dint1-0.5d0*(dum2t*dint2
     $                             +dum3t*dint3)
                              dum4t=eee2t(cp,ap,bp,ie,a,c,daaa,taa,dal
     $                             ,nasht)
                              amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b,c)+
     $                             dum1*dint1-0.5d0*(dum2*dint2+dum3
     $                             *dint3-dum4*dint4)
                              atmat(ap,bp,cp,a,b,c)=atmat(ap,bp,cp,a,b,c
     $                             )+0.5d0*dum4t*dint4
                              do ief=1,nasht
                                 iefc=ief+nisht
				 isf=itsym(iefc)
				 isfd=its(isf,isd)
				 isfde=its(isfd,ise)
                                 istot=its(isabcp,isfde)
				 if (its(isfde,isa).eq.1.and.its(istot
     $                                ,isbc).eq.1) then
                                 dum1=eeee(cp,ap,bp,b,d,ief,ie,c,daaa
     $                                   ,taa,dal,ro4,nasht)
                                 dum1t=-dum1
                                 if(bp.eq.b)dum1t=dum1t+2.d0*eee2(cp,ap
     $                                ,d,ief,ie,c,daaa,taa,dal,nasht)
                                 dint1=atwo(d,ief,ie,a)
                                 amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b
     $                                ,c)+dum1*dint1
                                 atmat(ap,b,cp,a,bp,c)=atmat(ap,b,cp,a
     $                                ,bp,c)+dum1t*dint1
				 endif
				 if (its(isfde,isc).eq.1.and.its(istot
     $                                ,isab).eq.1) then
                                 dum2=eeee(cp,ap,bp,b,d,ief,a,ie,daaa
     $                                   ,taa,dal,ro4,nasht)
                                 dum2t=-dum2
                                 if(bp.eq.b)dum2t=dum2t+2.d0*eee2(cp,ap
     $                                ,d,ief,a,ie,daaa,taa,dal,nasht)
                                 dint2=atwo(d,ief,c,ie)
                                 amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b
     $                                ,c)-dum2*dint2
                                 atmat(ap,b,cp,a,bp,c)=atmat(ap,b,cp,a
     $                                ,bp,c)-dum2t*dint2
				 endif
				 if (its(isfde,isb).eq.1.and.its(istot
     $                                ,isac).eq.1) then
                                 dum3=eeee(cp,ap,bp,ie,d,ief,a,c,daaa
     $                                ,taa,dal,ro4,nasht)
                                 dum3t=eeeet(cp,ap,bp,ie,d,ief,a,c,daaa
     $                                ,taa,dal,ro4,nasht)
                                 dint3=atwo(d,ief,b,ie)
                                 amat(ap,bp,cp,a,b,c)=amat(ap,bp,cp,a,b
     $                                ,c)-dum3*dint3
                                 atmat(ap,bp,cp,a,b,c)=atmat(ap,bp,cp,a
     $                                ,b,c)+dum3t*dint3
                                 endif
                              enddo
                           enddo
                        enddo
   1                 continue
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      deallocate(indice)
      return
      end
c-------------------------------------------------------------------
      subroutine bbmat(atwo,bmat,btmat,daaa,taa,dal,nasht,nisht,f,aj)
#include <implicit.h>
#include <priunit.h>
c---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht),taa(nasht
     $     ,nasht,nasht,nasht),dal(nasht,nasht)
      dimension bmat(nasht,nasht,nasht,nasht)
      dimension btmat(nasht,nasht,nasht,nasht)
      dimension atwo(nasht,nasht,nasht,nasht)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c---building matrix heff''
      do a=1,nasht
         ac=a+nisht
         do b=a,nasht
            bc=b+nisht
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nasht
               cc=c+nisht
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-dum
            enddo
         enddo
      enddo
c-------------------------------------
      do ap=1,nasht
         apc=ap+nisht
         do bp=1,nasht
            bpc=bp+nisht
            do cp=1,nasht
               cpc=cp+nisht
               do a=1,nasht
                  ac=a+nisht
                  bmat(ap,bp,cp,a)=0.d0
                  btmat(ap,bp,cp,a)=0.d0
                  do c=1,nasht
                     cc=c+nisht
                     iac=indice(cc,ac)
                     bmat(ap,bp,cp,a)=bmat(ap,bp,cp,a)-f(iac)
     $                    *ee2(cp,ap,bp,c,taa,dal,nasht)
                     btmat(ap,bp,cp,a)=btmat(ap,bp,cp,a)+aj(iac)
     $                    *ee2tr(cp,ap,bp,c,taa,dal,nasht)
                     do ie=1,nasht
                        iec=ie+nisht
                        do ief=1,nasht
                           iefc=ief+nisht
                           dum=eee2(cp,ap,bp,ie,c,ief,daaa,taa
     $                          ,dal,nasht)
                           dumt=eee2t(cp,ap,bp,ie,c,ief,daaa,taa
     $                          ,dal,nasht)
                           dint=atwo(a,ie,c,ief)
                           bmat(ap,bp,cp,a)=bmat(ap,bp,cp,a)-dum
     $                          *dint
                           btmat(ap,bp,cp,a)=btmat(ap,bp,cp,a)+dumt
     $                          *dint
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine bcmat(atwo,cmat,ctmat,daaa,taa,dal,nasht,nisht,f,aj)
#include <implicit.h>
#include <priunit.h>
c---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht),taa(nasht
     $     ,nasht,nasht,nasht),dal(nasht,nasht)
      dimension cmat(nasht,nasht,nasht,nasht)
      dimension ctmat(nasht,nasht,nasht,nasht)
      dimension atwo(nasht,nasht,nasht,nasht)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c--building heff'
      do a=1,nasht
         ac=a+nisht
         do b=a,nasht
            bc=b+nisht
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nasht
               cc=c+nisht
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-dum*0.5d0
            enddo
         enddo
      enddo
c-------------------
      do ap=1,nasht
         apc=ap+nisht
         do a=1,nasht
            ac=a+nisht
            do b=1,nasht
               bc=b+nisht
               do c=1,nasht
                  cc=c+nisht
                  cmat(ap,a,b,c)=0.d0
                  ctmat(ap,a,b,c)=0.d0
                  do d=1,nasht
                     dc=d+nisht
                     ida=indice(dc,ac)
                     icd=indice(cc,dc)
                     ibd=indice(bc,dc)
                     cmat(ap,a,b,c)=cmat(ap,a,b,c)+f(ida)*ee2(ap,b,d,c
     $                    ,taa,dal,nasht)-f(icd)*ee2(ap,b,a,d
     $                    ,taa,dal,nasht)-f(ibd)*ee2(ap,d,a,c
     $                    ,taa,dal,nasht)
                     ctmat(ap,a,b,c)=ctmat(ap,a,b,c)+f(ida)*ee2t(ap,b,d
     $                    ,c,taa,dal,nasht)-f(icd)*ee2t(ap,b,a,d,taa,dal
     $                    ,nasht)+f(ibd)*ee2t(ap,d,a,c,taa,dal,nasht)
                     do ie=1,nasht
                        iec=ie+nisht
                        dum1=ee2(ap,b,d,ie,taa,dal,nasht)
                        dum2=ee2(ap,b,d,c,taa,dal,nasht)
                        dum3=ee2(ap,b,a,ie,taa,dal,nasht)
                        dum4=ee2(ap,ie,a,c,taa,dal,nasht)
                        dum1t=ee2t(ap,b,d,ie,taa,dal,nasht)
                        dum2t=ee2t(ap,b,d,c,taa,dal,nasht)
                        dum3t=ee2t(ap,b,a,ie,taa,dal,nasht)
                        dum4t=ee2t(ap,ie,a,c,taa,dal,nasht)
                        dint1=atwo(c,ie,d,a)
                        dint2=atwo(d,ie,ie,a)
                        dint3=atwo(c,d,d,ie)
                        dint4=atwo(b,d,d,ie)
                        cmat(ap,a,b,c)=cmat(ap,a,b,c)+dum1*dint1-0.5d0
     $                       *(dum2*dint2+dum3*dint3-dum4*dint4)
                        ctmat(ap,a,b,c)=ctmat(ap,a,b,c)+dum1t*dint1-0
     $                       .5d0*(dum2t*dint2+dum3t*dint3-dum4t*dint4)
                        do ief=1,nasht
                           iefc=ief+nisht
                           dum1=eee2(ap,b,d,ief,ie,c,daaa,taa
     $                          ,dal,nasht)
                           dum1t=eee2tl(ap,b,d,ief,ie,c,daaa,taa
     $                          ,dal,nasht)
                           dint1=atwo(d,ief,ie,a)
                           dum2=eee2(ap,b,d,ief,a,ie,daaa,taa
     $                          ,dal,nasht)
                           dum2t=eee2tl(ap,b,d,ief,a,ie,daaa,taa
     $                          ,dal,nasht)
                           dint2=atwo(d,ief,c,ie)
                           dum3=eee2(ap,ie,d,ief,a,c,daaa,taa
     $                          ,dal,nasht)
                           dum3t=eee2tl(ap,ie,d,ief,a,c,daaa,taa
     $                          ,dal,nasht)
                           dint3=atwo(d,ief,b,ie)
                           cmat(ap,a,b,c)=cmat(ap,a,b,c)+dum1*dint1-dum2
     $                          *dint2-dum3*dint3
                           ctmat(ap,a,b,c)=ctmat(ap,a,b,c)+dum1t*dint1
     $                          -dum2t*dint2+dum3t*dint3
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      subroutine bdmat(atwo,dmat,dtmat,taa,dal,nasht,nisht,f,aj)
#include <implicit.h>
#include <priunit.h>
c---F e' la matrice dei monoelettronici con modif. Dyall e altra modif!
      dimension taa(nasht,nasht,nasht,nasht),dal(nasht,nasht)
      dimension dmat(nasht,nasht)
      dimension dtmat(nasht,nasht)
      dimension atwo(nasht,nasht,nasht,nasht)
      dimension f(*),aj(*)
      integer a,b,c,d,ac,bc,cc,dc,ap,bp,cp,apc,bpc,cpc
      parameter (two=2.d0)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c---building matrix heff''
      do a=1,nasht
         ac=a+nisht
         do b=a,nasht
            bc=b+nisht
            jab=indice(ac,bc)
            f(jab)=aj(jab)
            do c=1,nasht
               cc=c+nisht
               dum=atwo(a,c,c,b)
               f(jab)=f(jab)-dum
            enddo
         enddo
      enddo
c---------------------------
      do ap=1,nasht
         apc=ap+nisht
         do a=1,nasht
            ac=a+nisht
            dmat(ap,a)=0.d0
            dtmat(ap,a)=0.d0
            do c=1,nasht
               cc=c+nisht
               iac=indice(cc,ac)
               dmat(ap,a)=dmat(ap,a)-f(iac)*dal(ap,c)
               dum=-dal(ap,c)
               if(ap.eq.c)dum=dum+2.d0
               dtmat(ap,a)=dtmat(ap,a)+aj(iac)*dum
               do ie=1,nasht
                  iec=ie+nisht
                  do ief=1,nasht
                     iefc=ief+nisht
                     dum=ee2(ap,ie,c,ief,taa,dal,nasht)
                     dumt=ee2t(ap,ie,c,ief,taa,dal,nasht)
                     dint=atwo(a,ie,c,ief)
                     dmat(ap,a)=dmat(ap,a)-dum*dint
                     dtmat(ap,a)=dtmat(ap,a)+dumt*dint
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------------
      function eeee(a,b,c,d,e,f,g,h,d3,d2,d1,ro4,nasht)
#include <implicit.h>
      integer a,b,c,d,e,f,g,h
      dimension d3(nasht,nasht,nasht,nasht,nasht,nasht),d2(nasht,nasht
     $     ,nasht,nasht),d1(nasht,nasht)
#include "nevpt2.h"
      dimension ro4(nwords,nasht,nasht,nasht,nasht)
c----four-particle density matrix dealt with here
c--a,c,e,g    b,d,f,h
      lindice(i,j,k,l)=lindi(i)+lindj(j)+lindk(k)+l
      call ord8(a,c,e,g,b,d,f,h,i1,i2,i3,i4,i5,i6,i7,i8)
      ind=lindice(i1,i2,i3,i4)
      eeee=ro4(ind,i5,i6,i7,i8)
      if(b.eq.c)eeee=eeee+d3(a,e,g,d,f,h)
      if(b.eq.e)eeee=eeee+d3(a,c,g,f,d,h)
      if(b.eq.g)eeee=eeee+d3(a,c,e,h,d,f)
      if(d.eq.e)eeee=eeee+d3(a,g,c,b,h,f)
      if(d.eq.g)eeee=eeee+d3(a,c,e,b,h,f)
      if(f.eq.g)eeee=eeee+d3(a,e,c,b,h,d)
      if(b.eq.c)then
         if(d.eq.e)eeee=eeee+d2(a,g,f,h)
         if(d.eq.g)eeee=eeee+d2(a,e,h,f)
         if(f.eq.g)eeee=eeee+d2(a,e,d,h)
      endif
      if(b.eq.e.and.d.eq.g)eeee=eeee+d2(a,c,f,h)
      if(b.eq.e.and.f.eq.g)eeee=eeee+d2(a,c,h,d)
      if(b.eq.g.and.d.eq.e)eeee=eeee+d2(a,c,h,f)
      if(d.eq.e.and.f.eq.g)eeee=eeee+d2(a,c,b,h)
      if(b.eq.c.and.d.eq.e.and.f.eq.g)eeee=eeee+d1(a,h)
      return
      end
c-------------------------------------------------------------------
      function eeeet(a,b,c,d,e,f,g,h,d3,d2,d1,ro4,nasht)
#include <implicit.h>
      integer a,b,c,d,e,f,g,h
      dimension d3(nasht,nasht,nasht,nasht,nasht,nasht),d2(nasht,nasht
     $     ,nasht,nasht),d1(nasht,nasht)
#include "nevpt2.h"
      dimension ro4(nwords,nasht,nasht,nasht,nasht)
      eeeet=-eeee(a,b,d,c,e,f,g,h,d3,d2,d1,ro4,nasht)
      if(c.eq.d)eeeet=eeeet+2.d0*eee2(a,b,e,f,g,h,d3,d2,d1,nasht)
      return
      end
c-------------------------------------------------------------------
      function eee2(a,b,c,d,e,f,daaa,taa,dal,nasht)
#include <implicit.h>
      integer a,b,c,d,e,f
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht),taa(nasht
     $     ,nasht,nasht,nasht),dal(nasht,nasht)
c--the same as eee but with the spinless density matrices as input
      eee2=daaa(a,c,e,b,d,f)
      if(e.eq.d)eee2=eee2+taa(a,c,b,f)
      if(e.eq.b)eee2=eee2+taa(a,c,f,d)
      if(b.eq.c)eee2=eee2+taa(a,e,d,f)
      if(b.eq.c.and.e.eq.d)eee2=eee2+dal(a,f)
      return
      end
c-------------------------------------------------------------------
      function eee2t(a,b,c,d,e,f,daaa,taa,dal,nasht)
#include <implicit.h>
      integer a,b,c,d,e,f
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht),taa(nasht
     $     ,nasht,nasht,nasht),dal(nasht,nasht)
      eee2t=-eee2(a,b,d,c,e,f,daaa,taa,dal,nasht)
      if(c.eq.d)eee2t=eee2t+2.d0*ee2(a,b,e,f,taa,dal,nasht)
      return
      end
c-------------------------------------------------------------------
      function eee2tl(a,b,c,d,e,f,daaa,taa,dal,nasht)
#include <implicit.h>
      integer a,b,c,d,e,f
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht),taa(nasht
     $     ,nasht,nasht,nasht),dal(nasht,nasht)
      eee2tl=-eee2(b,a,c,d,e,f,daaa,taa,dal,nasht)
      if(a.eq.b)eee2tl=eee2tl+2.d0*ee2(c,d,e,f,taa,dal,nasht)
      return
      end
c-------------------------------------------------------------------
      function ee2(a,b,c,d,taa,dal,nasht)
#include <implicit.h>
      integer a,b,c,d
      dimension taa(nasht,nasht,nasht,nasht),dal(nasht,nasht)
      ee2=taa(a,c,b,d)
      if(b.eq.c)ee2=ee2+dal(a,d)
      return
      end
c-------------------------------------------------------------------
      function ee2t(a,b,c,d,taa,dal,nasht)
#include <implicit.h>
      integer a,b,c,d
      dimension taa(nasht,nasht,nasht,nasht),dal(nasht,nasht)
      ee2t=-ee2(b,a,c,d,taa,dal,nasht)
      if(a.eq.b)ee2t=ee2t+2.d0*dal(c,d)
      return
      end
c-------------------------------------------------------------------
      function ee2tr(a,b,c,d,taa,dal,nasht)
#include <implicit.h>
      integer a,b,c,d
      dimension taa(nasht,nasht,nasht,nasht),dal(nasht,nasht)
      ee2tr=-ee2(a,b,d,c,taa,dal,nasht)
      if(c.eq.d)ee2tr=ee2tr+2.d0*dal(a,b)
      return
      end
c----------------------------------------------------------------------
      function daaat(i,j,k,l,m,n,daaa,daa,da,nasht)
#include <implicit.h>
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht),daa(nasht
     $     ,nasht,nasht,nasht),da(nasht,nasht)
      v=-daaa(i,j,k,l,m,n)
      if(k.eq.n)v=v+daat(i,j,l,m,daa,da,nasht)
      if(j.eq.n)v=v-daat(i,k,l,m,daa,da,nasht)
      if(i.eq.n)v=v+daat(j,k,l,m,daa,da,nasht)
      if(m.eq.k)then
         v=v-daa(n,l,j,i)
         if(i.eq.l)v=v+da(n,j)
         if(j.eq.l)v=v-da(n,i)
      endif
      if(m.eq.j)then
         v=v+daa(n,l,k,i)
         if(k.eq.l)v=v+da(n,i)
         if(i.eq.l)v=v-da(n,k)
      endif
      if(m.eq.i)then
         v=v-daa(n,l,k,j)
         if(k.eq.l)v=v-da(n,j)
         if(j.eq.l)v=v+da(n,k)
      endif
      if(k.eq.l)v=v+daa(i,j,m,n)
      if(l.eq.j)v=v-daa(i,k,m,n)
      if(l.eq.i)v=v+daa(j,k,m,n)
      daaat=v
      return
      end
c----------------------------------------------------------------------
      function dbbbt(i,j,k,l,m,n,dbbb,dbb,db,nasht)
#include <implicit.h>
      dimension dbbb(nasht,nasht,nasht,nasht,nasht,nasht),dbb(nasht
     $     ,nasht,nasht,nasht),db(nasht,nasht)
      v=-dbbb(i,j,k,l,m,n)
      if(k.eq.n)v=v+dbbt(i,j,l,m,dbb,db,nasht)
      if(j.eq.n)v=v-dbbt(i,k,l,m,dbb,db,nasht)
      if(i.eq.n)v=v+dbbt(j,k,l,m,dbb,db,nasht)
      if(m.eq.k)then
         v=v-dbb(n,l,j,i)
         if(i.eq.l)v=v+db(n,j)
         if(j.eq.l)v=v-db(n,i)
      endif
      if(m.eq.j)then
         v=v+dbb(n,l,k,i)
         if(k.eq.l)v=v+db(n,i)
         if(i.eq.l)v=v-db(n,k)
      endif
      if(m.eq.i)then
         v=v-dbb(n,l,k,j)
         if(k.eq.l)v=v-db(n,j)
         if(j.eq.l)v=v+db(n,k)
      endif
      if(k.eq.l)v=v+dbb(i,j,m,n)
      if(l.eq.j)v=v-dbb(i,k,m,n)
      if(l.eq.i)v=v+dbb(j,k,m,n)
      dbbbt=v
      return
      end
c----------------------------------------------------------------------
      function daabt(i,j,k,l,m,n,daab,daa,dab,da,db,nasht)
#include <implicit.h>
      dimension daab(nasht,nasht,nasht,nasht,nasht,nasht),dab(nasht
     $     ,nasht,nasht,nasht),daa(nasht,nasht,nasht,nasht),da(nasht
     $     ,nasht),db(nasht,nasht)
      v=-daab(i,j,k,l,m,n)
      if(k.eq.n)v=v+daat(i,j,l,m,daa,da,nasht)
      if(m.eq.j)v=v+dab(l,n,i,k)
      if(m.eq.i)v=v-dab(l,n,j,k)
      if(l.eq.j)v=v-dab(i,k,m,n)
      if(l.eq.i)v=v+dab(j,k,m,n)
      if(m.eq.j.and.i.eq.l)v=v-db(n,k)
      if(m.eq.i.and.j.eq.l)v=v+db(n,k)
      daabt=v
      return
      end
c----------------------------------------------------------------------
      function dbbat(i,j,k,l,m,n,dbba,dbb,dab,da,db,nasht)
#include <implicit.h>
      dimension dbba(nasht,nasht,nasht,nasht,nasht,nasht),dab(nasht
     $     ,nasht,nasht,nasht),dbb(nasht,nasht,nasht,nasht),da(nasht
     $     ,nasht),db(nasht,nasht)
      v=-dbba(i,j,k,l,m,n)
      if(k.eq.n)v=v+dbbt(i,j,l,m,dbb,db,nasht)
      if(m.eq.j)v=v+dab(n,l,k,i)
      if(m.eq.i)v=v-dab(n,l,k,j)
      if(l.eq.j)v=v-dab(k,i,n,m)
      if(l.eq.i)v=v+dab(k,j,n,m)
      if(m.eq.j.and.i.eq.l)v=v-da(n,k)
      if(m.eq.i.and.j.eq.l)v=v+da(n,k)
      dbbat=v
      return
      end
c--------------------------------------------------------------------------
      function daat(a,b,c,d,daa,da,nasht)
#include <implicit.h>
      integer a,b,c,d
      dimension daa(nasht,nasht,nasht,nasht),da(nasht,nasht)
      v=daa(a,b,c,d)
      if(b.eq.d)v=v+dat(c,a,da,nasht)
      if(a.eq.d)v=v-dat(c,b,da,nasht)
      if(b.eq.c)v=v+da(a,d)
      if(a.eq.c)v=v-da(b,d)
      daat=v
      return
      end
c--------------------------------------------------------------------------
      function dbbt(a,b,c,d,dbb,db,nasht)
#include <implicit.h>
      integer a,b,c,d
      dimension dbb(nasht,nasht,nasht,nasht),db(nasht,nasht)
      v=dbb(a,b,c,d)
      if(b.eq.d)v=v+dbt(c,a,db,nasht)
      if(a.eq.d)v=v-dbt(c,b,db,nasht)
      if(b.eq.c)v=v+db(a,d)
      if(a.eq.c)v=v-db(b,d)
      dbbt=v
      return
      end
c--------------------------------------------------------------------------
      function dabt(a,b,c,d,dab,da,db,nasht)
#include <implicit.h>
      integer a,b,c,d
      dimension dab(nasht,nasht,nasht,nasht),da(nasht,nasht),db(nasht
     $     ,nasht)
      v=dab(a,b,c,d)
      if(b.eq.d)v=v+dat(c,a,da,nasht)
      if(a.eq.c)v=v-db(b,d)
      dabt=v
      return
      end
c------------------------------------------------------------------------
      function dat(a,b,da,nasht)
#include <implicit.h>
      integer a,b
      dimension da(nasht,nasht)
      v=-da(a,b)
      if(a.eq.b)v=v+1.0d0
      dat=v
      return
      end
c------------------------------------------------------------------------
      function dbt(a,b,db,nasht)
#include <implicit.h>
      integer a,b
      dimension db(nasht,nasht)
      v=-db(a,b)
      if(a.eq.b)v=v+1.0d0
      dbt=v
      return
      end
c------------------------------------------------------------------------
      subroutine koopE(atwo,koopipa,koopeaa,taa,dal,nasht,nisht,f
     $     ,aj)
#include <implicit.h>
#include <priunit.h>
      dimension taa(nasht,nasht,nasht,nasht),dal(nasht,nasht)
      dimension koopipa(nasht,nasht),koopeaa(nasht,nasht)
      dimension atwo(nasht,nasht,nasht,nasht)
      REAL*8  koopipa,koopipb,koopeaa,koopeab
      dimension f(*),aj(*)
      integer a,b,ac,bc,c,ap,apc,cc,d,dc
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
C
C     Calcolo le matrici di Koopmans
C     Caso -1: matrici di Koopmans per la ionizzazione attiva
C
C     La parte relativa alle singole differenze attive
C     e alle doppie differenze, una attiva e una inattiva
C     viene calcolata mediante rho_ab
C
c----building Dyall's h eff----------
      do a=1,nasht
cele
       ac=a+nisht
         do b=a,nasht
cele
          bc=b+nisht
            jab=indice(ac,bc)
            f(jab)=aj(jab)
         enddo
      enddo
      do a=1,nasht
         ac=a+nisht
         do ap=1,nasht
            apc=ap+nisht
            koopipa(ap,a)=0.d0
            do c=1,nasht
               cc=c+nisht
               irind=indice(ac,cc)
               dum=f(irind)
               koopipa(ap,a)=koopipa(ap,a)-dum*dal(ap,c)
c	       write (lupri,*) 'mono',a,ap,c,dum,dal(ap,c)
               do d=1,nasht
                  do ie=1,nasht
                     koopipa(ap,a)=koopipa(ap,a)-atwo(c,ie,a,d)*taa(ap,c
     $                    ,d,ie)
c      write (lupri,*) 'bi',ap,a,c,d,ie,atwo(c,ie,a,d)*taa(ap,c,d,ie)
c     $          ,atwo(c,ie,a,d),taa(ap,c,d,ie)
                  enddo
               enddo
            enddo
            koopipa(ap,a)=2.d0*koopipa(ap,a)
         enddo
      enddo
      do ap=1,nasht
         do a=1,nasht
            ind=indice(ap+nisht,a+nisht)
            dum=2.d0*f(ind)
            do d=1,nasht
               do ie=1,nasht
                  dum=dum+(2.d0*atwo(d,ie,ap,a)-atwo(ap,ie,d,a))*dal(d
     $                 ,ie)
               enddo
            enddo
            koopeaa(ap,a)=koopipa(ap,a)+2.d0*dum
         enddo
      enddo
      return
      end
c------------------------------------------------------
      subroutine ord31(ia,ic,ie,ig,ia2,ic2,ie2,ig2,val,amat,nasht,
     * m,n)
#include <implicit.h>
      integer segno,segno2,op,pp,op2,pp2,m,n,mp2,np2,mp,np
      integer ia,ic,ie,ig,ia2,ic2,ie2,ig2,nasht,ivolte,jvolte
      integer i1,i2,i3,i4,i5,i6,i7,i8,norm1,norm2,nwords,normind
      integer lindi,lindj,lindk
      DIMENSION amat(nwords,nasht,nasht,nasht,nasht)
#include "nevpt2.h"

      norm1=normind(ia,ic,ie,ig)
      norm2=normind(ia2,ic2,ie2,ig2)
      if (norm1.ge.norm2) then
      i1=ia
      i2=ic
      i3=ie
      i4=ig
      i5=ia2
      i6=ic2
      i7=ie2
      i8=ig2
      else
      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig
      endif

      do ivolte=1,6  !aaab
      call permut123(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6   
        call permut123(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,6  !aaba
      call permut124(i1,i2,i4,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6   
        call permut124(i5,i6,i8,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,6  !abaa
      call permut134(i1,i4,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
c      write (6,'(''ord31 3'',8I3)') mp,np,op,pp
        do jvolte=1,6   
        call permut134(i5,i8,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,6  !baaa
      call permut234(i4,i1,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
c      write (6,'(''ord31 4'',8I3)') mp,np,op,pp
        do jvolte=1,6   
        call permut234(i8,i5,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      if (norm1.ne.norm2) return
      if (m.eq.n) return

      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig

      do ivolte=1,6  !aaab
      call permut123(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6   
        call permut123(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,6  !aaba
      call permut124(i1,i2,i4,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6   
        call permut124(i5,i6,i8,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,6  !abaa
      call permut134(i1,i4,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6   
        call permut134(i5,i8,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,6  !baaa
      call permut234(i4,i1,i2,i3,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,6   
        call permut234(i8,i5,i6,i7,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      return
      end
c-----------------------------------------------------------
      subroutine ord22(ia,ic,ie,ig,ia2,ic2,ie2,ig2,val,amat,nasht,
     * m,n)
#include <implicit.h>
#include "nevpt2.h"
      DIMENSION amat(nwords,nasht,nasht,nasht,nasht)
      integer segno,segno2,op,pp,op2,pp2,m,n

      norm1=normind(ia,ic,ie,ig)
      norm2=normind(ia2,ic2,ie2,ig2)
      if (norm1.ge.norm2) then
      i1=ia
      i2=ic
      i3=ie
      i4=ig
      i5=ia2
      i6=ic2
      i7=ie2
      i8=ig2
      else
      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig
      endif

      do ivolte=1,4  !aabb
      call permut1234(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 1'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !abab
      call permut1324(i1,i3,i2,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 2'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i5,i7,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !abba
      call permut1423(i1,i3,i4,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 3'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i5,i7,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !baab
      call permut1423(i3,i1,i2,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 4'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i7,i5,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !baba
      call permut1324(i3,i1,i4,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 5'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i7,i5,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !bbaa
      call permut1234(i3,i4,i1,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 6'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i7,i8,i5,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      if (norm1.ne.norm2) return
      if (m.eq.n) return

      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig

c      write (6,'(''ord22 in'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8

      do ivolte=1,4  !aabb
      call permut1234(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 1'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !abab
      call permut1324(i1,i3,i2,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 2'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i5,i7,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !abba
      call permut1423(i1,i3,i4,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 3'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i5,i7,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !baab
      call permut1423(i3,i1,i2,i4,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 4'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1423(i7,i5,i6,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !baba
      call permut1324(i3,i1,i4,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 5'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1324(i7,i5,i8,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      do ivolte=1,4  !bbaa
      call permut1234(i3,i4,i1,i2,mp,np,op,pp,segno,ivolte)
c      write (6,'(''ord22 6'',8I3)') mp,np,op,pp
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,4
        call permut1234(i7,i8,i5,i6,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
        enddo
      endif
      enddo

      return
      end
c------------------------------------------------------
      subroutine ordsame(ia,ic,ie,ig,ia2,ic2,ie2,ig2,val,amat,nasht,
     * m,n)
#include <implicit.h>
#include "nevpt2.h"
      DIMENSION amat(nwords,nasht,nasht,nasht,nasht)
      integer segno,segno2,op,pp,op2,pp2,n,m

c      call getnorm(ia,ic,ie,ig,norm1)
c      call getnorm(ia2,ic2,ie2,ig2,norm2)
      norm1=normind(ia,ic,ie,ig)
      norm2=normind(ia2,ic2,ie2,ig2)
c      write (6,*) 'calling ordsame'
      if (norm1.ge.norm2) then
      i1=ia
      i2=ic
      i3=ie
      i4=ig
      i5=ia2
      i6=ic2
      i7=ie2
      i8=ig2
      else
      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig
      endif

      do ivolte=1,24
      call permut4(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,24
        call permut4(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
	enddo
      endif
      enddo

      if (norm1.ne.norm2) return
      if (m.eq.n) return

      i1=ia2
      i2=ic2
      i3=ie2
      i4=ig2
      i5=ia
      i6=ic
      i7=ie
      i8=ig

      do ivolte=1,24
      call permut4(i1,i2,i3,i4,mp,np,op,pp,segno,ivolte)
      if (mp.ge.np.and.np.ge.op.and.op.ge.pp) then
        do jvolte=1,24
        call permut4(i5,i6,i7,i8,mp2,np2,op2,pp2,segno2,jvolte)
        vals=val*segno*segno2
        call accum(amat,mp,np,op,pp,mp2,np2,op2,pp2,vals,nasht)
	enddo
      endif
      enddo


      return
      end
c-----------------------------------------------------
      subroutine ord9(ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2,
     * i1,i2,i3,i4,i5,i6,i7,i8)
#include <implicit.h>
c---  ordina in maniera decrescente ia1, ic1, ie1, ig1 in i1 i2 i3 i4
c---  ed fa la stessa perm anche su ia2, ic2, ie2, ig2 in i5 i6 i7 i8
      dimension nv1(4),nv2(4)
      common /actspace/nasht,nact2,nact3,nact4
#include "nevpt2.h"
c      write (6,'(''ord8 in'',8I3)') ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2

c      call getnorm(ia1,ic1,ie1,ig1,norm1)
c      call getnorm(ia2,ic2,ie2,ig2,norm2)
      norm1=normind(ia1,ic1,ie1,ig1)
      norm2=normind(ia2,ic2,ie2,ig2)


      if (norm1.ge.norm2) then
      nv1(1)=ia1
      nv1(2)=ic1
      nv1(3)=ie1
      nv1(4)=ig1
      nv2(1)=ia2
      nv2(2)=ic2
      nv2(3)=ie2
      nv2(4)=ig2
      else
      nv1(1)=ia2
      nv1(2)=ic2
      nv1(3)=ie2
      nv1(4)=ig2
      nv2(1)=ia1
      nv2(2)=ic1
      nv2(3)=ie1
      nv2(4)=ig1
      endif
      do i=1,4
         imax=nv1(i)
        do j=i+1,4
            if(nv1(j).gt.imax)then
               n=imax
               imax=nv1(j)
               nv1(j)=n
	       n=nv2(i)
	       nv2(i)=nv2(j)
	       nv2(j)=n
            endif
            nv1(i)=imax
         enddo
      enddo
      i1=nv1(1)
      i2=nv1(2)
      i3=nv1(3)
      i4=nv1(4)
      i5=nv2(1)
      i6=nv2(2)
      i7=nv2(3)
      i8=nv2(4)
c      write (6,'(''ord8 ou'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8         
      return
      end
c-----------------------------------------------------
      subroutine ord8(ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2,
     * i1,i2,i3,i4,i5,i6,i7,i8)
#include <implicit.h>
c---  ordina in maniera decrescente ia1, ic1, ie1, ig1 in i1 i2 i3 i4
c---  ed fa la stessa perm anche su ia2, ic2, ie2, ig2 in i5 i6 i7 i8
      dimension nv1(4),nv2(4)
      common /actspace/nasht,nact2,nact3,nact4
#include "nevpt2.h"
c      write (6,'(''ord8 in'',8I3)') ia1,ic1,ie1,ig1,ia2,ic2,ie2,ig2

c      call getnorm(ia1,ic1,ie1,ig1,norm1)
c      call getnorm(ia2,ic2,ie2,ig2,norm2)
      norm1=normind(ia1,ic1,ie1,ig1)
      norm2=normind(ia2,ic2,ie2,ig2)


      if (norm1.ge.norm2) then
      nv1(1)=ia1
      nv1(2)=ic1
      nv1(3)=ie1
      nv1(4)=ig1
      nv2(1)=ia2
      nv2(2)=ic2
      nv2(3)=ie2
      nv2(4)=ig2
      else
      nv1(1)=ia2
      nv1(2)=ic2
      nv1(3)=ie2
      nv1(4)=ig2
      nv2(1)=ia1
      nv2(2)=ic1
      nv2(3)=ie1
      nv2(4)=ig1
      endif
      do i=1,4
         imax=nv1(i)
	 jmax=i
         do j=i+1,4
            if(nv1(j).gt.imax)then
               imax=nv1(j)
               jmax=j  
            endif
         enddo
         if (jmax.ne.i) then
         ndum=nv1(i)
         nv1(i)=nv1(jmax)
         nv1(jmax)=ndum
         ndum=nv2(i)
         nv2(i)=nv2(jmax)
         nv2(jmax)=ndum
         endif
c      write (6,'(''ord8 ou'',i6,3x,8I3)') i,nv1(1),nv1(2),nv1(3),nv1(4),
c     *                                      nv2(1),nv2(2),nv2(3),nv2(4)
      enddo
      i1=nv1(1)
      i2=nv1(2)
      i3=nv1(3)
      i4=nv1(4)
      i5=nv2(1)
      i6=nv2(2)
      i7=nv2(3)
      i8=nv2(4)
c      write (6,'(''ord8 ou'',8I3)') i1,i2,i3,i4,i5,i6,i7,i8         
      return
      end
c-------------------------------------------------------------------
      function ro3t(a,b,c,ap,bp,cp,d3,d2,d1,nasht)
#include <implicit.h>
      integer a,b,c,ap,bp,cp
      dimension d3(nasht,nasht,nasht,nasht,nasht,nasht),d2(nasht,nasht
     $     ,nasht,nasht),d1(nasht,nasht)
      ro3t=-d3(ap,bp,cp,a,b,c)
      if(ap.eq.a)ro3t=ro3t+2.d0*d2(bp,cp,b,c)
      if(bp.eq.b)ro3t=ro3t+2.d0*d2(ap,cp,a,c)
      if(cp.eq.c)ro3t=ro3t+2.d0*d2(ap,bp,a,b)
      if(bp.eq.c)ro3t=ro3t-d2(ap,cp,a,b)
      if(ap.eq.c)ro3t=ro3t-d2(cp,bp,a,b)
      if(cp.eq.b)ro3t=ro3t-d2(bp,ap,c,a)
      if(ap.eq.b)ro3t=ro3t-d2(bp,cp,a,c)
      if(cp.eq.a)ro3t=ro3t-d2(ap,bp,c,b)
      if(bp.eq.a)ro3t=ro3t-d2(ap,cp,b,c)
      if(cp.eq.c.and.bp.eq.b)ro3t=ro3t-4.d0*d1(ap,a)
      if(cp.eq.b.and.bp.eq.c)ro3t=ro3t+2.d0*d1(ap,a)
      if(cp.eq.c.and.ap.eq.a)ro3t=ro3t-4.d0*d1(bp,b)
      if(cp.eq.a.and.ap.eq.c)ro3t=ro3t+2.d0*d1(bp,b)
      if(ap.eq.a.and.bp.eq.b)ro3t=ro3t-4.d0*d1(cp,c)
      if(ap.eq.b.and.bp.eq.a)ro3t=ro3t+2.d0*d1(cp,c)
      if(cp.eq.c.and.ap.eq.b)ro3t=ro3t+2.d0*d1(bp,a)
      if(cp.eq.b.and.ap.eq.c)ro3t=ro3t-d1(bp,a)
      if(cp.eq.c.and.bp.eq.a)ro3t=ro3t+2.d0*d1(ap,b)
      if(cp.eq.a.and.bp.eq.c)ro3t=ro3t-d1(ap,b)
      if(bp.eq.b.and.ap.eq.c)ro3t=ro3t+2.d0*d1(cp,a)
      if(bp.eq.c.and.ap.eq.b)ro3t=ro3t-d1(cp,a)
      if(bp.eq.b.and.cp.eq.a)ro3t=ro3t+2.d0*d1(ap,c)
      if(bp.eq.a.and.cp.eq.b)ro3t=ro3t-d1(ap,c)
      if(ap.eq.a.and.cp.eq.b)ro3t=ro3t+2.d0*d1(bp,c)
      if(ap.eq.b.and.cp.eq.a)ro3t=ro3t-d1(bp,c)
      if(ap.eq.a.and.bp.eq.c)ro3t=ro3t+2.d0*d1(cp,b)
      if(ap.eq.c.and.bp.eq.a)ro3t=ro3t-d1(cp,b)
      if(ap.eq.a.and.bp.eq.b.and.cp.eq.c)ro3t=ro3t+8.d0
      if(ap.eq.b.and.bp.eq.a.and.cp.eq.c)ro3t=ro3t-4.d0
      if(ap.eq.a.and.bp.eq.c.and.cp.eq.b)ro3t=ro3t-4.d0
      if(ap.eq.c.and.bp.eq.b.and.cp.eq.a)ro3t=ro3t-4.d0
      if(ap.eq.b.and.bp.eq.c.and.cp.eq.a)ro3t=ro3t+2.d0
      if(ap.eq.c.and.bp.eq.a.and.cp.eq.b)ro3t=ro3t+2.d0
      return
      end
c-------------------------------------------------------------------
      function ro2t(ap,bp,a,b,d2,d1,nasht)
#include <implicit.h>
      integer a,b,c,ap,bp,cp
      dimension d2(nasht,nasht,nasht,nasht),d1(nasht,nasht)
      ro2t=d2(a,b,ap,bp)
      if(ap.eq.a)ro2t=ro2t-2.d0*d1(b,bp)
      if(bp.eq.b)ro2t=ro2t-2.d0*d1(a,ap)
      if(ap.eq.b)ro2t=ro2t+d1(a,bp)
      if(bp.eq.a)ro2t=ro2t+d1(b,ap)
      if(ap.eq.a.and.bp.eq.b)ro2t=ro2t+4.d0
      if(ap.eq.b.and.bp.eq.a)ro2t=ro2t-2.d0
      return
      end
c-----------------------------------------------------------------
      subroutine bro3(ro4,daaa,nasht,nele)
#include <implicit.h>
#include <priunit.h>
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht)
#include "nevpt2.h"
      dimension ro4(nwords,nasht,nasht,nasht,nasht)
      lindice(i,j,k,l)=lindi(i)+lindj(j)+lindk(k)+l
cc---controllo su contrazione della ro4
      do i=1,nasht
         do j=1,nasht
            do k=1,nasht
               do ii=1,nasht
                  do jj=1,nasht
                     do kk=1,nasht
                        dum4=0.d0
                        do ll=1,nasht
			call ord8(i,j,k,ll,ii,jj,kk,ll,
     $                            i1,i2,i3,i4,i5,i6,i7,i8)
			ind=lindice(i1,i2,i3,i4)
                        dum4=dum4+ro4(ind,i5,i6,i7,i8)
                        enddo
                        daaa(i,j,k,ii,jj,kk)=dum4/(nele-3)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------
      subroutine bro2(daaa,taa,nasht,nele)
#include <implicit.h>
#include <priunit.h>
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht)
      dimension taa(nasht,nasht,nasht,nasht)
      do i=1,nasht
         do j=1,nasht
            do ii=1,nasht
               do jj=1,nasht
                  dum4=0.d0
                  do ll=1,nasht
                     dum4=dum4+daaa(i,j,ll,ii,jj,ll)
                  enddo
                  taa(i,j,ii,jj)=dum4/(nele-2)
               enddo
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------
      subroutine bro1(taa,dal,nasht,nele)
#include <implicit.h>
      dimension taa(nasht,nasht,nasht,nasht),dal(nasht,nasht)
      do i=1,nasht
         do j=1,nasht
            dum4=0.d0
            do ll=1,nasht
               dum4=dum4+taa(i,ll,j,ll)
            enddo
            dal(i,j)=dum4/(nele-1)
         enddo
      enddo
      return
      end
c------------------------------------------------------------------------
      subroutine wr1(dal,coopeaa,coopipa,nasht,lun31)
#include <implicit.h>
      dimension dal(nasht,nasht),coopeaa(nasht,nasht)
     $     ,coopipa(nasht,nasht)
      write(lun31) nasht
      write(lun31)((dal(k,l),k=1,nasht),l=1,nasht)
      write(lun31)((coopeaa(k,l),k=1,nasht),l=1,nasht)
      write(lun31)((coopipa(k,l),k=1,nasht),l=1,nasht)
      return
      end
c------------------------------------------------------------------------
      subroutine wr2(koopaa,taa,dal,nasht,lun31)
#include <implicit.h>
#include <priunit.h>
      REAL*8  KOOPAA
      DIMENSION dal(nasht,nasht),koopaa(nasht,nasht,nasht,nasht
     $     ),taa(nasht,nasht,nasht,nasht)
      do i=1,nasht
         do j=1,nasht
            do k=1,nasht
               do l=1,nasht
                  adum=ro2t(i,j,k,l,taa,dal,nasht)
                  koopaa(i,j,k,l)=adum
               enddo
            enddo
         enddo
      enddo
      do i=1,nasht
         do j=1,nasht
            write(lun31)((koopaa(i,j,k,l),k=1,nasht),l=1,nasht)
         enddo
      enddo
      return
      end
c------------------------------------------------------------------------
      subroutine wr3(koopeaa,taa,koopaa,nasht,lun31)
#include <implicit.h>
#include <priunit.h>
      REAL*8  koopaa(nasht,nasht,nasht,nasht),
     &           taa(nasht,nasht,nasht,nasht),
     &       koopeaa(nasht,nasht,nasht,nasht)
      do i=1,nasht
         do j=1,nasht
            write(lun31)((koopeaa(i,j,k,l),k=1,nasht),l=1,nasht)
         enddo
      enddo
c-------scrittura matrici densita a due particelle
      do i=1,nasht
         do j=1,nasht
            write(lun31)((taa(i,j,k,l),k=1,nasht),l=1,nasht)
         enddo
      enddo
c-------scrittura matrici koopman doppia ionizzazione
      do i=1,nasht
         do j=1,nasht
            write(lun31)((koopaa(i,j,k,l),k=1,nasht),l=1,nasht)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine wr4(koopaa,koopbb,koopab,nasht,lun31)
#include <implicit.h>
#include <priunit.h>
      REAL*8  KOOPAA, KOOPBB, KOOPAB
      DIMENSION koopaa(nasht,nasht,nasht,nasht
     $     ),koopbb(nasht,nasht,nasht,nasht),koopab(nasht,nasht)
      do i=1,nasht
         do j=1,nasht
            write(lun31)((koopaa(i,j,k,l),k=1,nasht),l=1,nasht)
            write(lun31)((koopbb(i,j,k,l),k=1,nasht),l=1,nasht)
         enddo
      enddo
      write(lun31)((koopab(k,l),k=1,nasht),l=1,nasht)
      return
      end
c------------------------------------------------------------------
      subroutine fdiff(bmat,cmat,btmat,ctmat,nasht)
#include <implicit.h>
#include <priunit.h>
      dimension bmat(nasht,nasht,nasht,nasht),cmat(nasht,nasht,nasht
     $     ,nasht)
      dimension btmat(nasht,nasht,nasht,nasht),ctmat(nasht,nasht,nasht
     $     ,nasht)
      bigdiff=0.d0
      do i =1,nasht
         do j=1,nasht
            do ii=1,nasht
               do jj=1,nasht
                  val=bmat(i,j,ii,jj)
                  val2=cmat(jj,i,j,ii)
                  pdiff=abs(val-val2)
                  if(pdiff.gt.bigdiff)bigdiff=pdiff
               enddo
            enddo
         enddo
      enddo
      write (lupri,'(/a,1P,d10.2)')
     & ' Largest difference between bmat  and cmat  =',bigdiff
      bigdiff=0.d0
      do i =1,nasht
         do j=1,nasht
            do ii=1,nasht
               do jj=1,nasht
                  val=btmat(i,j,ii,jj)
                  val2=ctmat(jj,i,j,ii)
                  pdiff=abs(val-val2)
                  if(pdiff.gt.bigdiff)bigdiff=pdiff
               enddo
            enddo
         enddo
      enddo
      write (lupri,'(/a,1P,d10.2)')
     & ' Largest difference between btmat and ctmat =',bigdiff
      call flshfo(lupri)
      return
      end
c------------------------------------------------------------------------------
      subroutine wr5(daaa,amat,atmat,bmat,btmat,cmat,ctmat,dmat,dtmat
     $     ,nasht,lun31)
#include <implicit.h>
#include <priunit.h>
      dimension bmat(nasht,nasht,nasht,nasht),cmat(nasht,nasht,nasht
     $     ,nasht)
      dimension amat(nasht,nasht,nasht,nasht,nasht,nasht),dmat(nasht
     $     ,nasht)
      dimension btmat(nasht,nasht,nasht,nasht),ctmat(nasht,nasht,nasht
     $     ,nasht)
      dimension atmat(nasht,nasht,nasht,nasht,nasht,nasht),dtmat(nasht
     $     ,nasht)
      dimension daaa(nasht,nasht,nasht,nasht,nasht,nasht)
c-----scrittura matrice spinless a 3 part e matrice amat
      do i=1,nasht
         do j=1,nasht
            do k=1,nasht
               do l=1,nasht
                  write(lun31)((daaa(i,j,k,l,ip,jp),ip=1,nasht),jp=1
     $                 ,nasht)
                  write(lun31)((amat(i,j,k,l,ip,jp),ip=1,nasht),jp=1
     $                 ,nasht)
               enddo
            enddo
         enddo
      enddo
c----scrittura matrici bmat e cmat
      do i=1,nasht
         do j=1,nasht
            write(lun31)((bmat(i,j,k,l),k=1,nasht),l=1,nasht)
            write(lun31)((cmat(i,j,k,l),k=1,nasht),l=1,nasht)
         enddo
      enddo
c----scrittura matrice dmat
      write(lun31)((dmat(i,j),i=1,nasht),j=1,nasht)
c-----scrittura matrice atmat
      do i=1,nasht
         do j=1,nasht
            do k=1,nasht
               do l=1,nasht
                  write(lun31)((atmat(i,j,k,l,ip,jp),ip=1,nasht),jp=1
     $                 ,nasht)
               enddo
            enddo
         enddo
      enddo
c----scrittura matrici btmat e ctmat
      do i=1,nasht
         do j=1,nasht
            write(lun31)((btmat(i,j,k,l),k=1,nasht),l=1,nasht)
            write(lun31)((ctmat(i,j,k,l),k=1,nasht),l=1,nasht)
         enddo
      enddo
c----scrittura matrice dtmat
      write(lun31)((dtmat(i,j),i=1,nasht),j=1,nasht)
      return
      end
c-----------------------------------------------------------
      subroutine giveocc(ioccfe,m,nd,ne,trou,part,nocc,norb)
#include <implicit.h>
#include <priunit.h>
      integer*2 ne(*),trou(*),part(*)
      integer*1 ioccfe(*)
      integer nd(*)
      do i=1,nocc
         ioccfe(i)=2
      enddo
      do i=nocc+1,norb
         ioccfe(i)=0
      enddo
      nstart=nd(m)
      do i=1,ne(m)
         ib=trou(nstart+i)
         ip=part(nstart+i)
         if(ib.gt.norb)ib=ib-norb
         if(ip.gt.norb)ip=ip-norb
         ioccfe(ib)=ioccfe(ib)-1
	 if (ip.le.norb) then
         ioccfe(ip)=ioccfe(ip)+1
	 endif
      enddo
      return
      end
c------------------------------------------------------------------
      integer function ndiff(ioccfe,jocc,norb)
#include <implicit.h>
      integer*1 ioccfe(*),jocc(*)
      ndiff=0
      do i=1,norb
         ia=ioccfe(i)-jocc(i)
         if(ia.ne.0)ndiff=ndiff+abs(ia)
         if(ndiff.gt.8)return
      enddo
      return
      end
c--------------------------------------------------------
      subroutine esclass(nd,ne,trou,part,icomp,iconf,isegno,ncf,nocc
     $     ,norb,ispin,iocc)
#include <implicit.h>
#include <priunit.h>
      integer*2 ne,trou,part,icomp,iconf,isegno
      dimension nd(*),ne(*),trou(*),part(*),icomp(*),iconf(*),isegno(*)
      equivalence (ip1,ib1),(ip2,ib2)
      integer*1 iocc(*)
      isz2=ispin-1
      ifirst=1
      ncapos=0
      do idet=1,ncf
         if(idet.eq.ifirst)then
c            print*,'calling giveocc with idet=',idet
c            call flush(6)
            call giveocc(iocc,idet,nd,ne,trou,part,nocc,norb)
            nsing=0
            do i=1,norb
               if(iocc(i).eq.1)nsing=nsing+1
            enddo
            nha=(nsing-isz2)/2
            nde=noverk(nsing,nha)
            ilast=ifirst+nde-1
            icomp(idet)=nde
            ncapos=ncapos+1
         endif
         if(idet.gt.ifirst.and.idet.le.ilast)then
            icomp(idet)=icomp(ifirst)
            goto 999
         endif
 999     if(idet.eq.ilast)ifirst=ilast+1
      enddo
c      print*,'esclass found ',ncapos,' capostipiti'
c      call flshfo(lupri)
      return
      end
c------------------------------------------------------------------------
      SUBROUTINE ESAME(flip,ne,nd,trou,part,ispinfe,iorbfe,ncf)

      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1 (Z)

#include <priunit.h>
      integer*2 ne,trou,part,ispinfe,iorbfe
      integer*1 flip
      dimension flip(*),ne(*),nd(*),trou(*),part(*),ispinfe(*),iorbfe(*)
      DIMENSION IT(12),IP(12)
      DO2M=1,NCF
      IF(M.EQ.1)GOTO5
      IF(FLIP(M-1).EQ.2)GOTO22
    5 NEC=NE(M)
      IF(NEC.EQ.0)GOTO3
      NAD=ND(M)
      ISP=0
      DO10I=1,NEC
      IT(I)=TROU(NAD+I)
      IP(I)=PART(NAD+I)
      IF(ispinfe(IT(I)).EQ.0)ISP=ISP-1
      IF(ispinfe(IT(I)).EQ.1)ISP=ISP+1
      IF(ispinfe(IP(I)).EQ.0)ISP=ISP+1
      IF(ispinfe(IP(I)).EQ.1)ISP=ISP-1
   10 CONTINUE
      IF(ISP.NE.0)RETURN
      IF((NEC/2)*2.NE.NEC)GOTO21
      DO20I=1,NEC
      ZREV1=.TRUE.
      ZREV2=.TRUE.
      IIT=iorbfe(IT(I))
      IIP=iorbfe(IP(I))
      DO30J=1,NEC
      IF(J.EQ.I)GOTO30
      JJT=iorbfe(IT(J))
      JJP=iorbfe(IP(J))
      IF(JJT.EQ.IIT)ZREV1=.FALSE.
      IF(JJP.EQ.IIP)ZREV2=.FALSE.
   30 CONTINUE
      ZREV=ZREV1.OR.ZREV2
      IF(ZREV)GOTO21
   20 CONTINUE
      FLIP(M)=0
      GOTO2
   21 FLIP(M)=2
      GOTO2
   22 FLIP(M)=1
      GOTO2
    3 FLIP(M)=0
    2 CONTINUE
      RETURN
      END
c-------------------------------------------------------------------

C
C Manu april 2009: DV (active density matrix in the MO basis) added as argument so that we can calculate
C the density of the MCSCF wave function in case we need it (for
C computing the srHxc potential in a NEVPT2-srDFT calculation for example)
C
               
      subroutine readint(onel,atwo,aj,dv,wrk,lfree,id2n,in2d,lupri)
C      
C      
#include <implicit.h>
#include <maxorb.h>
#include <inforb.h>
C Manu: for NEVPT2-srDFT
#include <inftap.h>
#include <infinp.h>
C#include <priunit.h>
      DIMENSION NEEDTP(6)
c      DATA NEEDTP /0,0,1,0,0,0/
      allocatable h2m(:,:)
      ALLOCATABLE CMO(:)
      ALLOCATABLE H1AO(:)
      ALLOCATABLE H1MO(:)
      dimension id2n(*),in2d(*),aj(*),wrk(*),onel(*)
      dimension atwo(nasht,nasht,nasht,nasht)
      integer a,b,ac,bc,acd,bcd
C Manu: needed for the NEVPT2-srDFT ...      
C
      DIMENSION DV(*)
      DOUBLE PRECISION EDFT 
      INTEGER LFREE
      PARAMETER (D1 = 1.0D0) 
C
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

C
      call dzero(aj,nnorbx)
      call dzero(atwo,n2ashx*n2ashx)
      do i=1,6
       needtp(i)=0
      enddo
      needtp(3)=1
      KFRSAV = 1
      KFREE = KFRSAV
      IDIST = 0
      allocate (h2m(norbt,norbt))
c      write (lupri,*) 'Starting reading the integrals'
c      call flshfo(lupri)
  100 CALL NXTH2M(IC,ID,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
c      write (lupri,*) 'read',ic,id
c      call flshfo(lupri)
      IF (IDIST .LT. 0) GO TO 800
C     ... if idist .lt. 0 then no more distributions
C
      icn=id2n(ic)
      if (icn.le.nisht.or.icn.gt.(nisht+nasht)) then
       write (lupri,*) 'Error IC not active'
       write (lupri,*) 'IC,ICN',ic,icn
       call quit('Error IC not active')
      endif 
      k=icn-nisht
      idn=id2n(id)
      if (idn.le.nisht.or.idn.gt.(nisht+nasht)) then
       write (lupri,*) 'Error ID not active'
       write (lupri,*) 'ID,IDN',id,idn
       call quit('Error ID not active')
      endif 
      l=idn-nisht
c      write (lupri,*) 'tranformed to nevpt',k,l


c---building matrix atwo
      do i=1,nasht
         do j=1,nasht
           id=in2d(i+nisht)
           jd=in2d(j+nisht)
c           write (lupri,*) 'nevpt index',i,j
c           write (lupri,*) 'dalt  index',id,jd
           atwo(i,j,k,l)=h2m(id,jd)
           atwo(i,j,l,k)=h2m(id,jd)
         enddo
      enddo

c----building the two el part of Dyall's h eff in aj (J contr)----------
      jab=indice(icn,idn)
      aj(jab)=0.d0
      do j=1,nisht
        jd=in2d(j)
        aj(jab)=aj(jab)+2.d0*h2m(jd,jd)
c	write (lupri,*) 'j',icn,idn,j,2.d0*h2m(jd,jd)
      enddo
       
C
      GO TO 100

  800 continue

c----Reading again the integrals: we need K for Heff 
      needtp(3)=0
      needtp(2)=1
      KFREE=1
      IDIST = 0
  200 CALL NXTH2M(IC,ID,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 900
C     ... if idist .lt. 0 then no more distributions
C
      icn=id2n(ic)
      idn=id2n(id)


c----building the two el part of Dyall's h eff in aj (K contr.)----------
      if (icn.le.idn) then ! idn is active and icn is core
        ac=idn
        icord=ic
      else  ! icn is active and idn is core
        ac=icn
        icord=id
      endif
      jn=id2n(icord)
      do bc=nisht+1,nisht+nasht
      jab=indice(bc,ac)
        bcd=in2d(bc)
        dum=h2m(bcd,icord)
        if (ac.ne.bc) dum=dum*.5d0
        aj(jab)=aj(jab)-dum
c	write (lupri,*) 'k',bc,ac,jn,dum
      enddo
       
C
      GO TO 200

  900 continue

c      do i=1,nasht
c         do j=1,nasht
c           do k=1,nasht
c              do l=1,nasht
c              write (lupri,*) i,j,k,l,atwo(i,j,k,l)
c              enddo
c            enddo
c         enddo
c      enddo
C
C     Generate 1-electron Hamiltonian matrix in MO basis
C
      ALLOCATE (CMO(NCMOT))
      ALLOCATE (H1AO(NNBAST))
      ALLOCATE (H1MO(NNORBT))
C
C     We read from LUIT1 using NEWORB label
C
      JRDMO =9
      CALL READMO(CMO,JRDMO,WRK(KFREE),LFREE)
c      CALL PRORB(CMO,.FALSE.,LUPRI) ! Write MO sym packed
C
C
      CALL SIRH1(H1AO,WRK(KFREE),LFREE)
C
c     CALL OUTPKB(H1AO,NBAS,NSYM,1,LUPRI) ! Write H1AO sym packed 

C
C************ NEVPT2-srDFT -- start ***************
C
C March 2009 Emmanuel Fromager

C Purpose: add the frozen srHXC potential calculated
C for the MCSCF-srDFT density to the one-electron part of the 
C Hamiltonian. 
C

      IF (DOMCSRDFT) THEN
#ifndef MOD_SRDFT
         call quit('srdft not implemented in this version')
#else

C        Allocation of memory for ...       

C        full (inactive+active) DM in the AO basis             
         CALL MEMGET('REAL',KDAO,N2BASX,WRK,KFREE,LFREE) 
C        sr H or xc potential in the AO basis          
         CALL MEMGET('REAL',KVSRAO,N2BASX,WRK,KFREE,LFREE) 
C        active part of the DM in the AO basis         
         CALL MEMGET('REAL',KDVAO,N2BASX,WRK,KFREE,LFREE) 

C Manu b
!        write(lupri,*) 'print DV before calling FCKDEN'
!        call outpak(dv,nasht,1,lupri)
C Manu e

C        Backtransform the inactive and active MCSCF-srDFT density matrix to the AO
C        basis       

         CALL FCKDEN(.TRUE.,.TRUE.,WRK(KDAO),WRK(KDVAO),CMO,DV,
     &               WRK(KFREE),LFREE)

C Manu b
         
!        write(lupri,*) 'print DVAO just after calling FCKDEN'
!        CALL OUTPUT(WRK(KDVAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
!        write(lupri,*) 'print DCAO just after calling FCKDEN'
!        CALL OUTPUT(WRK(KDAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
C Manu e

C        make total density matrix (in the AO basis) in WORK(KDAO)

         CALL DAXPY(N2BASX,D1,WRK(KDVAO),1,WRK(KDAO),1)
C

C        Release memory allocated for DVAO 
C
         CALL MEMREL('Releasing mem. used for DVAO',
     &             WRK,KFRSAV,KDVAO,KFREE,LFREE)  

C        calculate the srH potential (for the MCSCF-srDFT density) 

         CALL DZERO(WRK(KVSRAO),N2BASX)

         IFCTMP = 11 ! sr Hartree term only
         NSRDMAT = 1 ! only one density matrix (the MCSCF-srDFT one)
         ISYMDM = 1  ! DMAT of symmetry 1 

         CALL SIRFCK2(LUSRINT,'AOSR2INT',WRK(KVSRAO),
     &                WRK(KDAO),NSRDMAT,ISYMDM,IFCTMP,WRK(KFREE),LFREE)


C manu b
!        write(lupri,*) 'print VSRAO after calling sirfck2'
!        CALL OUTPUT(WRK(KVSRAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
!        write(lupri,*) 'print H1AO after calling sirfck2'
!        CALL OUTPKB(H1AO,NBAS,NSYM,1,LUPRI) ! Write H1AO sym packed 
       
C manu e

C        Add the srH potential to the one-electron part of the Hamiltonian

         CALL ADD_SU_SP(WRK(KVSRAO),H1AO)
C manu b
!        write(lupri,*) 'print H1AO after 1st call to ADD_SU_SP'
!        CALL OUTPKB(H1AO,NBAS,NSYM,1,LUPRI) ! Write H1AO sym packed 
       
C manu e

C        Calculate the srxc potential for the MCSCF-srDFT density  

         CALL DZERO(WRK(KVSRAO),N2BASX)

         IPRFCK = 0

         CALL SRDFT(1,WRK(KVSRAO),WRK(KDAO),EDFT,.TRUE.,.FALSE.,.FALSE.,
     &              WRK(KFREE),LFREE,IPRFCK)

C Comment Manu: we do not care about the energy contribution EDFT. It
C is just that it is calculated with the srxc potential ...

C manu b
!        write(lupri,*) 'print VSRAO after calling srdft'
!        CALL OUTPUT(WRK(KVSRAO),1,NBAST,1,NBAST,NBAST,NBAST,1,LUPRI)
!        write(lupri,*) 'print H1AO after calling srdft'
!        CALL OUTPKB(H1AO,NBAS,NSYM,1,LUPRI) ! Write H1AO sym packed 
       
C manu e

C        Add the srxc potential to the one-electron part of the Hamiltonian

         CALL ADD_SU_SP(WRK(KVSRAO),H1AO)

C        Release memory

         CALL MEMREL('Releasing mem. used for DAO and VSRAO',
     &             WRK,KFRSAV,KDAO,KFREE,LFREE)  

C
C************ NEVPT2-srDFT -- end ***************
C
#endif
      END IF 
         
      CALL UTHUB(H1AO,H1MO,CMO,WRK(KFREE),NSYM,NBAS,NORB)
      CALL PKSYM1(onel,H1MO,NORB,NSYM,-1) ! Unpack H1MO in low tri onel
c     CALL OUTPKB(H1MO,NORB,NSYM,1,LUPRI)
!     write(lupri,*) 'print onel (unpacked H1MO)'
c     CALL OUTPAK(onel,NORBT,1,LUPRI)

c----building Dyall's h eff in aj----------
      do a=1,nasht
         ac=a+nisht
         acd=in2d(ac)
         do b=a,nasht
            bc=b+nisht
            bcd=in2d(bc)
            jab=indice(ac,bc)
            jabd=indice(acd,bcd)
c           write (lupri,*) ac,bc,onel(jabd),acd,bcd
c	    call flshfo(lupri)
            aj(jab)=aj(jab)+onel(jabd)
c           write (lupri,*) 'one',ac,bc,aj(jab)
c	    call flshfo(lupri)
         enddo
      enddo

      deallocate (CMO)
      deallocate (H1AO)
      deallocate (H1MO)

      return
      end

c-----------------------------------------------------
C  /* Deck ADD_SU_SP */
      SUBROUTINE ADD_SU_SP(VSRAO,H1AO)
C
C April-2009 Emmanuel Fromager 
C
C This subroutine adds the symmetric UNPACKED matrix VSRAO to the
C symmetry PACKED H1AO matrix 
C
#include <implicit.h>
      DIMENSION VSRAO(NBAST,NBAST),H1AO(*)   
C 
#include <inforb.h>
#include <priunit.h>
C 

         IH1AO = 0 
         IVSRAO = 0
         DO 2500 ISYM = 1,NSYM
            NBASI  = NBAS(ISYM)
C ***** initial code - assuming packed wrt column ... ******
C           DO 2600 J = 1,NBASI 
C              DO 2700 I = J,NBASI
C***********************************************************     
            DO 2600 I = 1,NBASI 
               DO 2700 J = 1,I
C Manu b
C        write(lupri,*) 'debug printing in ADD_SU_SP ...'
C        write(lupri,*) 'H1AO(',
C    &   IH1AO+(I-J+1)+(J-1)*NBASI-((J-1)*(J-2)/2),') =',
C    &   H1AO(IH1AO+(I-J+1)+(J-1)*NBASI-((J-1)*(J-2)/2))
C        write(lupri,*) 'VSRAO(', IVSRAO+I,',',IVSRAO+J,
C    &   ')=', VSRAO(IVSRAO+I,IVSRAO+J) 
C Manu e
C ***** initial code - assuming packed wrt column ... ******
C                 H1AO(IH1AO+(I-J+1)+(J-1)*NBASI-((J-1)*(J-2)/2)) =
C    &            H1AO(IH1AO+(I-J+1)+(J-1)*NBASI-((J-1)*(J-2)/2)) +
C    &            VSRAO(IVSRAO+I,IVSRAO+J)    
C***********************************************************     
                  H1AO(IH1AO+J+(I*(I-1)/2)) =
     &            H1AO(IH1AO+J+(I*(I-1)/2)) +
     &            VSRAO(IVSRAO+I,IVSRAO+J)    
 2700          CONTINUE
 2600       CONTINUE
            IVSRAO = IVSRAO + NBASI
            IH1AO  = IH1AO + ((NBASI*(NBASI+1))/2)
 2500    CONTINUE
C
      RETURN
      END
C
c-----------------------------------------------------
      subroutine getnorm(ia,ic,ie,ig,norm)
#include <implicit.h>
#include <priunit.h>
      dimension nv(4)
      integer ia,ic,ie,ig,norm,nv,i,j,jmax,imax,ndum

c      norm=ia+ic+ie+ig
c      return

      nv(1)=ia
      nv(2)=ic
      nv(3)=ie
      nv(4)=ig
      

c      write (6,'(''norm in'',8I3)') ia,ic,ie,ig

      do i=1,4
       imax=nv(i)
       jmax=i
        do j=i+1,4
	if(nv(j).gt.imax) then
	imax=nv(j)
	jmax=j
	endif
	enddo
	if(jmax.ne.i) then
	ndum=nv(i)
	nv(i)=nv(jmax)
	nv(jmax)=ndum
	endif
      enddo
c      write (6,'(''norm ou'',8I3)') nv(1),nv(2),nv(3),nv(4)

      norm=(nv(1)-1)*nv(1)*(nv(1)+1)*(nv(1)+2)/24+
     *     (nv(2)-1)*nv(2)*(nv(2)+1)/6+
     *     nv(3)*(nv(3)-1)/2+nv(4)

      return
      end
c-----------------------------------------------------
      subroutine nwadd(nwall,nadd)
#include <single.h>
#include <implicit.h>
#include <priunit.h>
      nwall=nwall+nadd
c      write (lupri,*) 'Alloc  :',nadd/8,nwall/8
c      write (lupri,*) 'Alloc  :',nadd,nwall
c      call flshfo(lupri)
      return
      end
c-----------------------------------------------------
      subroutine nwdel(nwall,ndel)
#include <single.h>
#include <implicit.h>
#include <priunit.h>
      nwall=nwall-ndel
c      write (lupri,*) 'Dealloc:',ndel/8,nwall/8
c      write (lupri,*) 'Dealloc:',ndel,nwall
c      call flshfo(lupri)
      return
      end
      function vald2r8(a,i,j,n)
#include <implicit.h>
      dimension a(n,n)
      vald2r8 = a(i,j)
      return
      end
#endif
! --- end of koopro4.F ---
