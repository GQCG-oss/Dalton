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
      subroutine dalcip(nela,nelb,ndet,c,ne,trou,part,nd,nm,
     *           inact,iact,nact,lfermi,nisht,nsym,LUN98,mult)
#include <implicit.h>
#include <priunit.h>

! ** Before 2011, the nevpt2 code was not compiled unless PGF90, GFORTRAN or IFORT compiler.
! ** Now code is compiled for all compilers, if your compiler cannot handle
! ** the code then do NOT #define GO_ON_COMP in dalcip.F (below) and koopro4.F /Jan. 2011

! #if !(defined(VAR_IFORT)||defined(VAR_PGF90)||defined(VAR_GFORTRAN))
!       return
!       end
! #else
#define GO_ON_COMP
! #endif      
#ifdef GO_ON_COMP
      dimension c(*),ne(*),trou(*),part(*),nd(*)
c      dimension c(*),ast(nela,*),bst(nelb,*),ne(*),trou(*),part(*),nd(*)
c      integer*1 ast,bst
      integer*2 ne,trou,part,ispin(600),iorb(600)
      integer astr(1600),bstr(1600),iecci(0:100)
      character*1 zwspin(600),zplus,zmoins
      dimension inact(8),iact(8)
c     integer tri
c     tri(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      zplus='+'
      zmoins='-'
      ngelo=0
      ina=0
      nocc=lfermi+nisht
      do i=1,nsym
         ina=ina+inact(i)
      enddo
      norb=nm
      nd(1)=0
      do  i=1,norb
         ispin(i)=0
         zwspin(i)=zmoins
         ispin(i+norb)=1
         zwspin(i+norb)=zplus
         iorb(i)=i
         iorb(i+norb)=i
      enddo
c--renzo
      if(mod((nela-nelb),2).ne.0)then !per un doppietto, quadrupletto,...
         if(nela.ne.nelb)then
            ispin(norb+norb+1)=1
            zwspin(norb+norb+1)=zplus
            iorb(norb+norb+1)=norb+1
         endif
      endif
      rewind LUN98
      do idet=1,ndet
         read(LUN98)c(idet),(astr(i),i=1,nela),(bstr(i),i=1,nelb)
c remco
c         if (dabs(c(idet)).gt.0.10) then
c            write(2,'(F10.6,5x,14I3)')c(idet),(astr(jj),jj=1,nela)
c            write(2,'(15x,14I3)')(bstr(jj),jj=1,nelb)
c         endif 
c remco
c        if(nela.eq.nelb)then
         if(mult.eq.1)then
            call bupa(astr,bstr,norb,nela,nelb,ne,trou,part,nd,ina,
     $           idet,isign)
c        elseif (mod((nela-nelb),2).eq.1)then !doublets, quartets,...
         elseif (mod(mult,2).eq.0)then !doublets, quartets,...
            call bupa2(astr,bstr,norb,nela,nelb,ne,trou,part,nd,ina,
     $           idet,isign,lfermi)
         else                   !triplets, quintets,...
            call bupa3(astr,bstr,norb,nela,nelb,ne,trou,part,nd,ina,
     $           idet,isign)
         endif
         if(isign.eq.-1)c(idet)=-c(idet)
         nec=ne(idet)
         if(idet.lt.ndet)nd(idet+1)=nd(idet)+nec
        enddo
      call ascend(nd,ne,trou,part,c,ndet,nm,iecci,nexmax,nocc,norb
     $     ,mult,iorb,ispin)
c     if(nela.eq.nelb)then
c      if(mult.eq.1)then
c         call orderfe3(nexmax,iecci,ndet,c,nd,ne,trou,part,nm,ina,ast,
c     $        bst,nela)
c      endif
      return
      end
c---------------------------------------------------------------
      subroutine bupa2(astr,bstr,norb,nela,nelb,ne,trou,part,nd,
     $     ina,idet,isign,nacto)
      dimension astr(*),bstr(*),trou(*),part(*),ne(*),nd(*)
      dimension c(1600)
      integer*2 trou,part,ne,istri(20),jstri(20)
      integer astr,bstr,c
      logical zhole
      ip=0
      irest=0
      ni=nd(idet)
      do i=1,20
         istri(i)=0
         jstri(i)=0
      enddo
      do i=1,nelb
         if(bstr(i).gt.nacto)then
            ip=ip+1
            part(ni+ip)=bstr(i)+ina
         else
            irest=irest+1
            c(irest)=bstr(i)
         endif
      enddo
      it=0
      do i=1,nacto
         zhole=.true.
         do j=1,irest
            if(c(j).eq.i)then
               zhole=.false.
               goto 1
            endif
         enddo
 1       continue
         if(zhole)then
            it=it+1
            trou(ni+it)=i+ina
         endif
      enddo
      irest=0
      do i=1,nela
         if(astr(i).gt.nacto)then
            ip=ip+1
            part(ni+ip)=astr(i)+ina+norb
         else
            irest=irest+1
            c(irest)=astr(i)
         endif
      enddo
      if(nelb.lt.nela)then
         ip=ip+1
         part(ni+ip)=norb+norb+1 !practically infinity
      endif
      do i=1,nacto
         zhole=.true.
         do j=1,irest
            if(c(j).eq.i)then
               zhole=.false.
               goto 2
            endif
         enddo
 2       continue
         if(zhole)then
            it=it+1
            trou(ni+it)=i+ina+norb
         endif
      enddo
      nec=it
      ne(idet)=nec
      do i=1,nacto
         istri(i)=i  !vacuum state definition
         jstri(i)=i  !vacuum state definition
      enddo
      imore=0
      do i=1,nec
         it=trou(ni+i)
         ip=part(ni+i)
         if(it.gt.norb)then
            if(ip.gt.norb)then
               istri(it-norb-ina)=ip-norb-ina
            else
               istri(it-norb-ina)=0
            endif
         else
            if(ip.gt.norb)then
               jstri(it-ina)=0
               imore=imore+1
               istri(nacto+imore)=ip-norb-ina !for triplets
            else
               jstri(it-ina)=ip-ina
            endif
         endif
      enddo
      isign=1
      do i=1,nacto+imore   !this is for triplets (and for doublets a la cipsi)
         do j=i+1,nacto+imore !same remark a s above
            if(istri(i).gt.istri(j))isign=-isign
         enddo
      enddo
      do i=1,nacto
         do j=i+1,nacto
            if(jstri(i).gt.jstri(j))isign=-isign
         enddo
      enddo
      return
      end
c**********************************************
      subroutine bupa3(astr,bstr,norb,nela,nelb,ne,trou,part,nd,
     $     ina,idet,isign)
c--routine for triplets (high spin)
      dimension astr(*),bstr(*),trou(*),part(*),ne(*),nd(*)
      dimension c(1600)
      integer*2 trou,part,ne,istri(20),jstri(20)
      integer astr,bstr,c
      logical zhole
      ip=0
      irest=0
      ni=nd(idet)
      nel=(nela+nelb)/2
      do i=1,nelb
         if(bstr(i).gt.nel)then
            ip=ip+1
            part(ni+ip)=bstr(i)+ina
         else
            irest=irest+1
            c(irest)=bstr(i)
         endif
      enddo
      it=0
      do i=1,nel
         zhole=.true.
         do j=1,irest
            if(c(j).eq.i)then
               zhole=.false.
               goto 1
            endif
         enddo
 1       continue
         if(zhole)then
            it=it+1
            trou(ni+it)=i+ina
         endif
      enddo
      irest=0
      do i=1,nela
         if(astr(i).gt.nel)then
            ip=ip+1
            part(ni+ip)=astr(i)+ina+norb
         else
            irest=irest+1
            c(irest)=astr(i)
         endif
      enddo
      do i=1,nel
         zhole=.true.
         do j=1,irest
            if(c(j).eq.i)then
               zhole=.false.
               goto 2
            endif
         enddo
 2       continue
         if(zhole)then
            it=it+1
            trou(ni+it)=i+ina+norb
         endif
      enddo
      nec=it
      ne(idet)=nec
      do i=1,nel
         istri(i)=i  !vacuum state definition
         jstri(i)=i  !vacuum state definition
      enddo
      imore=0
      do i=1,nec
         it=trou(ni+i)
         ip=part(ni+i)
         if(it.gt.norb)then
            if(ip.gt.norb)then
               istri(it-norb-ina)=ip-norb-ina
            else
               istri(it-norb-ina)=0
            endif
         else
            if(ip.gt.norb)then
               jstri(it-ina)=0
               imore=imore+1
               istri(nel+imore)=ip-norb-ina !for triplets
            else
               jstri(it-ina)=ip-ina
            endif
         endif
      enddo
      isign=1
      do i=1,nel+imore   !this is for triplets (and for doublets a la cipsi)
         do j=i+1,nel+imore !same remark a s above
            if(istri(i).gt.istri(j))isign=-isign
         enddo
      enddo
      do i=1,nel
         do j=i+1,nel
            if(jstri(i).gt.jstri(j))isign=-isign
         enddo
      enddo
      return
      end
c**************************************************************
      subroutine bupa(astr,bstr,norb,nela,nelb,ne,trou,part,nd,
     $     ina,idet,isign)
      dimension astr(*),bstr(*),trou(*),part(*),ne(*),nd(*)
      dimension c(1600)
      integer*2 trou,part,ne,istri(20),jstri(20)
      integer astr,bstr,c
      logical zhole
      ip=0
      irest=0
      ni=nd(idet)
      do i=1,nelb
         if(bstr(i).gt.nela)then
            ip=ip+1
            part(ni+ip)=bstr(i)+ina
         else
            irest=irest+1
            c(irest)=bstr(i)
         endif
      enddo
      it=0
      do i=1,nela
         zhole=.true.
         do j=1,irest
            if(c(j).eq.i)then
               zhole=.false.
               goto 1
            endif
         enddo
 1       continue
         if(zhole)then
            it=it+1
            trou(ni+it)=i+ina
         endif
      enddo
      irest=0
      do i=1,nela
         if(astr(i).gt.nela)then
            ip=ip+1
            part(ni+ip)=astr(i)+ina+norb
         else
            irest=irest+1
            c(irest)=astr(i)
         endif
      enddo
      do i=1,nela
         zhole=.true.
         do j=1,irest
            if(c(j).eq.i)then
               zhole=.false.
               goto 2
            endif
         enddo
 2       continue
         if(zhole)then
            it=it+1
            trou(ni+it)=i+ina+norb
         endif
      enddo
      nec=it
      ne(idet)=nec
      do i=1,nela
         istri(i)=i  !vacuum state definition
         jstri(i)=i  !vacuum state definition
      enddo
      do i=1,nec
         it=trou(ni+i)
         ip=part(ni+i)
         if(it.gt.norb)then
            istri(it-norb-ina)=ip-norb-ina
         else
            jstri(it-ina)=ip-ina
         endif
      enddo
      isign=1
      do i=1,nela
         do j=i+1,nela
            if(istri(i).gt.istri(j))isign=-isign
         enddo
      enddo
      do i=1,nela
         do j=i+1,nela
            if(jstri(i).gt.jstri(j))isign=-isign
         enddo
      enddo
      return
      end
c****************************************************************
      subroutine orderfe(ifirst,ilast,c,ne,nd,trou,part,npos)
#include <implicit.h>      
      DIMENSION c(*)
      integer*2 ne(*),trou(*),part(*)
      dimension nd(*),npos(*)
      do idet=ifirst,ilast
         nec=ne(idet)
         inu=idet-ifirst+1
         if(npos(inu).eq.inu) goto 1
         do jdet=idet+1,ilast
            jnu=jdet-ifirst+1
            if(npos(jnu).eq.inu)then
               call swapbp(nec,idet,jdet,nd,trou,part)
               cdum=c(idet)
               c(idet)=c(jdet)
               c(jdet)=cdum
               ndum=npos(inu)
               npos(inu)=npos(jnu)
               npos(jnu)=ndum
            endif
         enddo
 1    continue
      enddo
      return
      end
c----------------------------------------------------------
      subroutine giveoccij(iocc,m,nd,ne,trou,part,nocc,norb)
      integer*2 ne(*),trou(*),part(*)
      integer*1 iocc(*)
      integer nd(*)
      do i=1,nocc
         iocc(i)=2
      enddo
      do i=nocc+1,norb
         iocc(i)=0
      enddo
      nstart=nd(m)
      do i=1,ne(m)
         ib=trou(nstart+i)
         ip=part(nstart+i)
         if(ib.gt.norb)ib=ib-norb
         if(ip.gt.norb)ip=ip-norb
         iocc(ib)=iocc(ib)-1
         if (ip.le.norb) then
         iocc(ip)=iocc(ip)+1
         endif
      enddo
      return
      end
C------------------------------------------------------------------
      integer function ndiffij(iocc,jocc,norb)
      integer*1 iocc(*),jocc(*)
      ndiffij=0
      do i=1,norb
         ia=iocc(i)-jocc(i)
         if(ia.ne.0)then
            ndiffij=1
            return
         endif
      enddo
      return
      end
c*************************************************************
         subroutine swapbp(nec,idet,jdet,nd,trou,part)
         dimension nd(*),trou(*),part(*)
         dimension troud(100),partd(100)
         integer*2 trou,part,troud,partd
         ndi=nd(idet)
         do i=1,nec
            troud(i)=trou(ndi+i)
            partd(i)=part(ndi+i)
         enddo
         ndj=nd(jdet)
         do i=1,nec
            trou(ndi+i)=trou(ndj+i)
            part(ndi+i)=part(ndj+i)
            trou(ndj+i)=troud(i)
            part(ndj+i)=partd(i)
         enddo
         return
         end
c******************************************************************
       subroutine givepos(iel,n,k,npos)
c
c      n=numero di orbitali singoli occupati
c      k=numero di elettroni spaiati alpha
c      iel(i),i=1,k = posizione dell'i-esimo elettrone 
c                     nelle n caselle
c
       integer iel
       dimension iel(*)
       logical*1 zsz0

       zsz0=n.eq.k*2

       if (zsz0) then
          ntot=1
          do i=k+1,n
             ntot=ntot*i
          enddo
          do i=1,k
             ntot=ntot/i
          enddo
       endif

       if (k.eq.1) then
          npos=iel(1)
          goto 110
       endif

       npos=1
       kpos=0
       kposold=0
       do ijk=1,k-1
          kpos=iel(ijk)-ijk-kposold
       if (kpos.eq.0) goto 111
          nogg=k-ijk
          do kl=ijk+kposold,iel(ijk)-1
             ncasel=n-kl
             ndistr=1
             do l=nogg+1,ncasel
                ndistr=ndistr*l
             enddo
             do l=1,ncasel-nogg
                ndistr=ndistr/l
             enddo
             npos=npos+ndistr
          enddo
          kposold=kposold+kpos
 111   continue
       enddo
       npos=npos+iel(k)-iel(k-1)-1
       
       
 110   if (zsz0) then
          if (npos*2.le.ntot) then
             npos=npos*2-1
          else
             npos=(ntot-npos)*2+2	 
          endif
       endif
       
       return
       end
c---------------------------------------------------------------
      subroutine ascend(nd,ne,trou,part,c,ndet,nm,iecci,nexmax,nocc,norb
     $     ,mult,iorb,ispin)
#include <implicit.h>
      DIMENSION c(*)
      dimension nd(*),ne(*),trou(*),part(*),iorb(*),ispin(*)
      integer*2 ne,trou,part,ibu(100),ipa(100),iorb,ispin
      integer iecci(0:100)
      integer*1 o
      integer*1 occu
      allocatable occu(:,:),iel(:),npos(:)
      nexmax=0
      do i=0,100
         iecci(i)=0
      enddo
      do i=1,ndet
         if(ne(i).gt.nexmax) nexmax=ne(i)
      enddo
      LUN15=-13300
      CALL GPOPEN(LUN15,'scr15','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     *           .FALSE.)
      do iec=0,nexmax
         do idet=1,ndet
            nec=ne(idet)
            if(nec.eq.iec)then
               iecci(iec)=iecci(iec)+1
               ndi=nd(idet)
               write(lun15)c(idet),nec,(trou(i),part(i),i=ndi+1,ndi+nec)
            endif
         enddo
      enddo
      rewind lun15
      nd(1)=0
      idet1=1
      allocate(occu(norb,ndet))
      allocate(iel(norb))
c      allocate(npos(norb))   !wrong should be maximum of a deter space
      allocate(npos(1000))    !covers 12-open shell
      do iec=0,nexmax
       do idet=idet1,idet1+iecci(iec)-1
        read(lun15)c(idet),nec,(ibu(i),ipa(i),i=1,nec)
        ne(idet)=nec
        if(idet.lt.ndet)nd(idet+1)=nd(idet)+nec
        ndi=nd(idet)
        do i=1,nec
         trou(ndi+i)=ibu(i)
         part(ndi+i)=ipa(i)
        enddo
        call giveoccij(occu(1,idet),idet,nd,ne,trou,part,nocc,norb)
       enddo
       ifirst=idet1
       do idet=idet1,idet1+iecci(iec)-1
        nec=ne(idet)
        inext=idet+1
        if(idet.eq.ifirst)then
           nsing=0
           do i=1,norb
            if(occu(i,idet).eq.1)nsing=nsing+1
           enddo
           isz2=mult-1
           nha=(nsing-isz2)/2
           nde=noverk(nsing,nha)
           ilast=ifirst+nde-1
           ifound=1
        endif
        if(idet.gt.ifirst.and.idet.le.ilast)goto 999
        do jdet=idet+1,idet1+iecci(iec)-1
         if(ndiff(occu(1,idet),occu(1,jdet),norb).eq.0)then
            ifound=ifound+1
            if(jdet.ne.inext)then
               call swapbp(nec,inext,jdet,nd,trou
     $              ,part)
               cdum=c(inext)
               c(inext)=c(jdet)
               c(jdet)=cdum
               do io=1,norb
                o=occu(io,inext)
                occu(io,inext)=occu(io,jdet)
                occu(io,jdet)=o
               enddo
            endif
            inext=inext+1
         endif
         if(ifound.eq.nde)goto 999
        enddo
 999    if(idet.eq.ilast)then
           do kdet=ifirst,ilast
            isingle=0
            indalf=0
            do k=1,norb
             if(occu(k,kdet).eq.1)then
                isingle=isingle+1
                ndk=nd(kdet)
                do i=1,ne(kdet)
                 ispb=ispin(trou(i+ndk))
                 iorbb=iorb(trou(i+ndk))
                 ispp=ispin(part(i+ndk))
                 iorbp=iorb(part(i+ndk))
                 if((iorbp.eq.k.and.ispp.eq.1).or.
     $              (iorbb.eq.k.and.ispb.eq.0))then
                    indalf=indalf+1
                    iel(indalf)=isingle
                    if(indalf.gt.norb)then
                       print*,'ERROR: indalf .gt. norb',indalf,norb
                       CALL QUIT('givepos: indalf .gt. norb')
                    endif
                 endif
                enddo
             endif
            enddo
            if((kdet-ifirst+1).gt.1000)then
               print*,'givepos: arg of npos too large'
               CALL QUIT('givepos: arg of npos too large')
            endif
            if (indalf.gt.0) then
               call givepos(iel,isingle,indalf,npos(kdet-ifirst+1))
            else
               npos(kdet-ifirst+1) = -999 999 999
            end if
           enddo
           call orderfe(ifirst,ilast,c,ne,nd,trou,part,npos)
           ifirst=ilast+1
        endif
       enddo
       idet1=idet1+iecci(iec)
      enddo
      CALL GPCLOSE(LUN15,'DELETE')
c     close(15,status='DELETE')
      return
      end
c****************************************************************
      integer function noverk(n,k)
      if(k.eq.0)then
         noverk=1
         return
      elseif(k.eq.1)then
         noverk=n
         return
      else
         num=n
         iden=1
         do i=2,k
            num=num*(n-i+1)
            iden=iden*i
         enddo
         noverk=num/iden
         return
      endif
      end
#endif
