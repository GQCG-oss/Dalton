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
*=====================================================================*
      subroutine cc_r12mkxajkl(kvajkl,cmo,work,lwork,fel)
c---------------------------------------------------------------------
c     purpose: calculate X_{aj}^{kl}, X = V,T necessary for the
c     calculation of MP2-R12 contributions to the right hand side
c     of Z-Vector Equations (cc_r12grad)
c---------------------------------------------------------------------

      implicit none
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      logical ldum,test,fel
      parameter (test = .false.,zero = 0.0d0, one = 1.0d0)

      integer  lwork,luvajkl,lwrk1,lwrk2,isymai
      integer isymkl,isymij,isymmu,isyml,isymk,isymak,isymaj,icount7,
     &     ntota,isymv,isymmuj,koffl,koffv,kvmujkl,
     &     isymj,isyma,idxkl,xajkl,iglmrv(8,8)
      integer kend1,kend2,koffc,isymi,ntklaj,yajkl
      integer icou2,isymmua,ntotmu,nr1bas(8),ir1bas(8,8),icount2
      integer noffro,ntota1,noffset,isym,nvir0(8)
      integer idxij,idxji,kcmoa,idum,ntotgu,ntotj,lu43
      integer nklaj(8),nr1vir(8),iglmro(8,8),icou3
      character*10 fnelena
      data fnelena/'DDR12xxxxx'/
      real*8 cmo(*),work(*),kvajkl(*)
c
      call qenter('cc_r12mkxajkl')
      fnelena = fnelena(1:5)//fnvajkl(6:10)

      if (fnvajkl .eq. 'CCR12QAJKL' .and. fel .eqv. .true.) then
         lu43 = -43
         call dzero(cmo,nlamds)
         if (R12HYB) then
             call gpopen(lu43,'GUMAT.1','UNKNOWN',' ','UNFORMATTED',
     &                      IDUM,LDUM)
         else
             call quit('Sorry, calculation of R12 corrections to'// 
     &                  ' first order properties not implemented'//
     &                  ' with .NO HYB')
         endif 
         REWIND(LU43)
         READ(LU43) NTOTGU
         READ(LU43) (cmo(I), I = 1, NTOTGU)
         CALL GPCLOSE(LU43,'KEEP')
      endif
c
      if (fnvajkl .eq. 'CCR12QAJKL' .and. fel .eqv. .false.) then
         kend1 = 1
         lwrk1 = lwork - kend1
         lu43 = -43
         call dzero(cmo,nlamds)
         call gpopen(lu43,'SIRIFC','OLD',' ','UNFORMATTED',
     &                      IDUMMY,.FALSE.)
         REWIND(LU43)
         CALL MOLLAB(LABEL,LU43,LUPRI)
         READ(LU43)
         READ(LU43)
         READ(LU43) (cmo(I), I = 1, NLAMDS)
         CALL GPCLOSE(LU43,'KEEP')
         CALL CMO_REORDER(CMO,WORK(KEND1),LWRK1)
 
      endif
c

      do isymak = 1, nsym
        nr1bas(isymak) = 0
        nr1vir(isymak) = 0
        do isymk = 1, nsym
           isyma = muld2h(isymak,isymk)
           nvir0(isyma)    = norb1(isyma) - nrhfs(isyma)  
           nr1bas(isymak)  = nr1bas(isymak) +mbas1(isyma)*nrhf(isymk)
           nr1vir(isymak)  = nr1vir(isymak) +nvir0(isyma)*nrhf(isymk)
        end do
      end do

      do isymak = 1, nsym
         nvajkl(isymak) = 0
         icount7 = 0
         icount2 = 0
         do isymk = 1, nsym
            isyma = muld2h(isymak,isymk)
            ivajkl(isyma,isymk) = icount7
            nvajkl(isymak) = nvajkl(isymak)+ nr1bas(isyma)*nmatij(isymk)
            icount7 = icount7 + nr1bas(isyma)*nmatij(isymk)
            ir1bas(isyma,isymk)  = icount2
            icount2 = icount2 + nrhf(isymk)*mbas1(isyma)
         end do
      end do

      ntklaj = 0
      do isymkl = 1,nsym
         isymaj = isymkl
         nklaj(isymaj) = 0
         do isymj = 1,nsym
            isyma = muld2h(isymaj,isymj)
            nklaj(isymaj) = nklaj(isymaj) + nmatij(isymj)*nr1vir(isyma)
            ntklaj = ntklaj + nmatij(isymkl)*nr1vir(isymj)
         enddo
      enddo

      do isymmua = 1,nsym
         icou2 = 0
         icou3 = 0
         do isyma = 1,nsym
            isymmu = muld2h(isymmua,isyma)
            iglmrv(isymmu,isyma) = icou2
            iglmro(isymmu,isyma) = icou3
            icou2 = icou2 + nbas(isymmu)*(nvir(isyma)-nrhffr(isyma))
            icou3 = icou3 + nbas(isymmu)*(norb1(isyma)-nrhffr(isyma))
         enddo
      enddo


      noffro = 0
      do isym = 1,nsym
         noffro = noffro + nbas(isym)*nrhf(isym)
      enddo


      kend1 = 1
      yajkl = kend1
      kcmoa =  yajkl + ntklaj
      kend2  = kcmoa + nlamds
      lwrk2  = lwork - kend2

      if (lwrk2 .lt.0) then
         call quit('Insufficient work space in cc_r12mkxajkl')
      end if
c
c
      call dzero(work(yajkl),ntklaj)
      xajkl = yajkl
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymmuj = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  noffset = 0
                  do isymmu =1, nsym
                     isyma = muld2h(1,isymmu)
                     isymj  = muld2h(isymmuj,isymmu)
                     nvir0(isyma) = norb1(isyma) - nrhfs(isyma)
                     if (fel) then 
                        koffc = 1 + iglmro(isymmu,isyma) 
     &                          +nrhf(isymmu)*nbas(isymmu) 
                     else
                        noffset = noffset+ nbas(isymmu)*nrhffr(isymmu)
                        koffc = 1 + nrhf(isymmu)*nbas(isymmu)
     &                        +iglmrv(isymmu,isyma)+noffset
     &                        +noffro
                     endif
                     kvmujkl = 1+ivajkl(isymmuj,isymkl)+
     &                         nr1bas(isymmuj)*(idxkl-1)
     &                         +ir1bas(isymmu,isymj)
                     ntotmu = max(1,mbas1(isymmu))
                     ntota1 = max(1,nbas(isymmu))
                     ntota = max(1,nvir0(isyma))
                     ntotj = max(1,nrhf(isymj))
c
                     call dgemm('T','N',nvir0(isyma),
     &                 nrhf(isymj),mbas1(isymmu),one,cmo(koffc),
     &                 ntota1,kvajkl(kvmujkl),ntotmu,zero,
     &                 work(xajkl),ntota)

                     xajkl = xajkl + nvir0(isyma)*nrhf(isymj)
                  end do
               end do
            end do
         end do
      end do

      if (r12cbs) call dscal(nvajkl(1),1.00d0,work(yajkl),1)


      luvajkl = -1
      call gpopen(luvajkl,fnelena,'unknown',' ','unformatted',
     &               idummy,.false.)
      rewind(luvajkl)
      write(luvajkl) (work(yajkl+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'KEEP')


      call qexit('cc_r12mkxajkl')
      end
C=====================================================================*
C  /* Deck cc_r12grad */
      subroutine cc_r12grad(ipdd,orb,v2am,work,lwork)
C     Calculate MP2-R12 contributions to the right hand side of Z-Vector
C     Equations
      implicit none
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>
#include <cbirea.h>
#include <r12int.h>
#include <priunit.h>
#include <ccfro.h>
      logical locdbg
      logical lauxd,ldum
      parameter (locdbg = .false.)
      character*10 CCR12YAJIJ
      integer kend1,kccr12,kend2,lwrk1,lwork,nrhftria,isymij
      integer lu43,idum,n2,lwrk2,ipdd
      integer krsing,krtrip,kend3,lwrk3,nr1bas(8)
      integer isymkl,isyml,isymk,idxij,idxji,idxklij,idxklji
      integer isymai,isymv,isymi,isymj,idxkl,kctil,ntklaj
      integer kbajkl,kbijal,kvajkl,kvijal,isymaj,isyma,luvajkl
      integer lcvai,koffc,nvir0(8),isymmu,ntota,ntotmu,kcvai
      integer ir1bas(8,8),isymmuj
      integer dklij,ntotk,lcvia,kvmujkl,ntoti
      integer koffd,bdajij,lccr12
      integer dklijz,kajij,kcvia,ntotj,imataj(8,8),icou3
      integer nr1orb(8),nr1xbas(8),nrgkl,n2bst1(8),ir1orb(8,8)
      integer ir2bas(8,8),ir2xbas(8,8),irgkl,iaodis1(8,8),nr2bas(8)
      integer ir1xbas(8,8),nijkl(8)
      integer koffe,kctil2,idxijkl,lctil,ncvai
      integer res,kscr1

#if defined (SYS_CRAY)
      real work(lwork),one,two,zero
      real orb(*),v2am(*)
#else
      double precision work(lwork),one,zero,two,orb(*),v2am(*)
#endif
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)


      call qenter('cc_r12grad')


      nrhftria = nrhft*(nrhft+1)/2
      n2 = nrhftria * nrhftria

      kccr12  = 1
      krsing  = kccr12  + ngamsq(1)
      krtrip  = krsing  + n2
      kctil   = krtrip  + n2
      kctil2  = kctil   + ngamsq(1)
      kend1   = kctil2  + ngamsq(1)
      lwrk1   = lwork   - kend1

      call dzero(work(kccr12),ngamsq(1))
      call dzero(work(krsing),n2)
      call dzero(work(krtrip),n2)
      call dzero(work(kctil),ngamsq(1))
      call dzero(work(kctil2),ngamsq(1))

      isymv = 1

      if (lwrk1 .lt. 0) then
         call quit('Insufficient work space in cc_r12grad')
      END IF

c     ----------------------------------------
c     read R12 amplitudes
c     ----------------------------------------


      lu43 = -43
      call gpopen(lu43,'CCR12_C','unknown',' ','formatted',
     &                   idummy,.false.)
      read(lu43,'(4e30.20)') (work(krsing+i), i = 0, n2-1 )
      read(lu43,'(4e30.20)') (work(krtrip+i), i = 0, n2-1 )
      call gpclose(lu43,'KEEP')

      call dscal(nrhftria*nrhftria,0.5d0,work(krsing),1)
      call dscal(nrhftria*nrhftria,0.5d0,work(krtrip),1)

      call cc_r12vpack(work(kccr12),isymv,work(krsing),work(krtrip),
     &     nrhf,nrhfa,nijkl)

c     ---------------------------------
C     get c_tilde = 2*c_ij^kl - c_ji^kl
c     ---------------------------------

      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  do isymi = 1, nsym
                     isymj = muld2h(isymij,isymi)
                     do j = 1, nrhf(isymj)
                        do i = 1, nrhf(isymi)
                           idxij = imatij(isymi,isymj)+
     &                          nrhf(isymi)*(j-1) + i
                           idxji = imatij(isymj,isymi)+
     &                          nrhf(isymj)*(i-1) + j
                           idxklij =  igamsq(isymkl,isymij)+
     &                          nmatij(isymkl)*(idxij-1)+idxkl
                           idxklji =igamsq(isymkl,isymij)+
     &                          nmatij(isymkl)*(idxji-1)+idxkl
                           work(kctil-1+idxklij) =
     &                          two*work(kccr12-1+idxklij) -
     &                          work(kccr12-1+idxklji)
                        end do
                     end do

                  end do
               end do
            end do
         end do
      end do

      if (locdbg) then
         write(lupri,*) 'resorted R12 amplitudes (ctilde):'
         do isymkl = 1, nsym
            isymij = muld2h(isymkl,1)
            write(lupri,*) 'symmetry block ',isymkl,isymij
            call output(work(kctil+igamsq(isymkl,isymij)),
     &           1,nmatij(isymkl),1,nmatij(isymij),
     &            nmatij(isymkl),nmatij(isymij),1,lupri)
         end do
      endif 

c     ---------------------------------
C     Transpone c_tilde
c     ---------------------------------

      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  do isymi = 1, nsym
                     isymj = muld2h(isymij,isymi)
                     do j = 1, nrhf(isymj)
                        do i = 1, nrhf(isymi)
                           idxij = imatij(isymi,isymj)+
     &                          nrhf(isymi)*(j-1) + i
                           idxji = imatij(isymj,isymi)+
     &                          nrhf(isymj)*(i-1) + j
                           idxklij =  igamsq(isymij,isymkl)+
     &                          nmatij(isymkl)*(idxij-1)+idxkl
                           idxklji =igamsq(isymij,isymkl)+
     &                          nmatij(isymkl)*(idxji-1)+idxkl
                           idxijkl = igamsq(isymkl,isymij)+
     &                          nmatij(isymij)*(idxkl-1) +idxij
                           work(kctil2-1+idxklij) =
     &                          work(kctil-1+idxijkl)
                        end do
                     end do

                  end do
               end do
            end do
         end do
      end do


c
C---------------------------------------------
C     Calculate some offsets and dimensions
C---------------------------------------------
      ntklaj = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  do isyma =1, nsym
                     isymj  = muld2h(isymaj,isyma)
                     isymi  = isyma
                     ntklaj = ntklaj + nrhf(isymj) *
     &                        (norb1(isyma) - nrhfs(isyma))
                  enddo
               enddo
            enddo
         enddo
      enddo


      icou3 = 0
      do isymai = 1,nsym
         icou3 = 0
         do isyma = 1,nsym
            isymi = muld2h(isymai,isyma)
            imataj(isyma,isymi)= icou3
            icou3 = icou3+nrhf(isymi)*(norb1(isyma) - nrhfs(isyma))
         enddo
      enddo

      do isymaj = 1,nsym
         isymij = isymaj
        ncvai = 0
         do isyma = 1,nsym
            isymj = muld2h(isymaj,isyma)
            isymi = muld2h(isymij,isymj)
            ncvai = ncvai + (norb1(isyma) - nrhfs(isyma))*nrhf(isymi)
         enddo
      enddo


      kbajkl  = kend1
      kbijal  = kbajkl  + ntklaj
      kvajkl  = kbijal  + ntklaj
      kvijal  = kvajkl  + ntklaj
      kcvai   = kvijal  + ntklaj
      kcvia   = kcvai   + ncvai
      kajij   = kcvia   + ncvai
      dklij   = kajij   + ncvai
      res     = dklij   + ngamsq(1)
      kend2   = res     + ncvai
      lwrk2   = lwork   - kend2

c     ----------------------------------------------------
c     Read V_{aj}^{kl},V_{kl}^{aj},B_{aj}^{kl},B_{kl}^{aj}
c               and transform it with Ctilde
c     -----------------------------------------------------

      if (lwrk2 .lt.0) then
         call quit('Insufficient work space in cc_r12grad')
      end if

      call dzero(work(kbajkl),ntklaj)
      call dzero(work(kbijal),ntklaj)
      luvajkl = -1
      call gpopen(luvajkl,'DDR12BAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kbajkl+i-1), i = 1,ntklaj)
      if (ipdd .eq. 5) then
          call gpclose(luvajkl,'DELETE')
      else
          call gpclose(luvajkl,'KEEP')
      endif


      call gpopen(luvajkl,'DDR12BIJAL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kbijal+i-1), i = 1,ntklaj)
      if (ipdd .eq. 5) then
          call gpclose(luvajkl,'DELETE')
      else
          call gpclose(luvajkl,'KEEP')
      endif
c
      call daxpy(ntklaj,one,work(kbajkl),1,work(kbijal),1)

      call dzero(work(kvajkl),ntklaj)
      call dzero(work(kvijal),ntklaj)
      call gpopen(luvajkl,'DDR12VAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kvajkl+i-1), i = 1,ntklaj)
      if (ipdd .eq. 5) then
          call gpclose(luvajkl,'DELETE')
      else
          call gpclose(luvajkl,'KEEP')
      endif

      call gpopen(luvajkl,'DDR12VIJAL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kvijal+i-1), i = 1,ntklaj)
      if (ipdd .eq. 5) then
          call gpclose(luvajkl,'DELETE')
      else
          call gpclose(luvajkl,'KEEP')
      endif



c     Calculate d_{kl}^{ij} = c_{kl}^{mn}*ctilde_{ij}^{mn}

      call dzero(work(dklij),ngamsq(1))
      do isymkl = 1,nsym
         isymij  = isymkl
         lccr12 = kccr12 + igamsq(isymij,isymkl)
         lctil  = kctil + igamsq(isymij,isymkl)
         dklijz = dklij + igamsq(isymij,isymkl)
         ntoti = max(1,nmatij(isymij))
         ntotk = max(1,nmatij(isymkl))
         call dgemm('T','N',nmatij(isymij),nmatij(isymkl),
     &              nmatij(isymkl),one,work(lccr12),
     &              ntotk,work(lctil),ntoti,zero,
     &              work(dklijz),ntotk)
      end do


      call dzero(work(kcvai),ncvai)
      lcvai = kcvai
      call dzero(work(kcvia),ncvai)
      lcvia = kcvia
      call dzero(work(kajij),ncvai)
      bdajij = kajij

c     Calculate V_{aj}^{kl}*ctilde_{kl}^{ij}


      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhf(isymk)
             do l = 1, nrhf(isyml)
                idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    isymi  = muld2h(isymij,isymj)
                    nvir0(isyma)  = norb1(isyma) - nrhfs(isyma)
                    koffc  = kctil + igamsq(isymij,isymkl)
     &                 + nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    ntotmu = max(1,nvir0(isyma))
                    ntotj = max(1,nrhf(isymj))
                    ntoti = max(1,nrhf(isymi))
                    lcvai  = kcvai + imataj(isyma,isymi)
                    call dgemm('N','T',nvir0(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kvajkl),
     &                    ntotmu,work(koffc),ntoti,one,
     &                    work(lcvai),ntotmu)

                    kvajkl = kvajkl + nvir0(isyma)*nrhf(isymj)

c     Calculate V_{kl}^{aj}*ctilde_{ij}^{kl}

                    lcvia  = kcvia + imataj(isyma,isymi)
                    koffe = kctil2+ igamsq(isymkl,isymij)
     &              +  nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    call dgemm('N','T',nvir0(isyma),nrhf(isymi),
     &                         nrhf(isymj),one,work(kvijal),
     &                         ntotmu,work(koffe),ntoti,one,
     &                         work(lcvia),ntotmu)
                    kvijal = kvijal + nvir0(isyma)*nrhf(isymj)
c
c     Calculate B_{aj}^{kl}*d_{kl}^{ij}
c
                    koffd = dklij + igamsq(isymij,isymkl)
     &              +  nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    bdajij  = kajij + imataj(isyma,isymi)
                    call dgemm('N','T',nvir0(isyma),nrhf(isymi),
     &                         nrhf(isymj),one,work(kbijal),
     &                         ntotmu,work(koffd),ntoti,one,
     &                         work(bdajij),ntotmu)
                    kbijal  = kbijal + nvir0(isyma)*nrhf(isymj)
                end do
              end do
            end do
         end do
      end do
      call dscal(ncvai,two,work(kcvai),1)
      call dscal(ncvai,two,work(kcvia),1)
      call daxpy(ncvai,one,work(kcvai),1,work(kcvia),1)
      call daxpy(ncvai,one,work(kcvia),1,work(kajij),1)
      call dscal(ncvai,two,work(kajij),1)
      if (.not. lmulbs) call dscal(ncvai,-one,work(kajij),1)

      if (locdbg) then
         write(lupri,*) 'Matrixelemente MP2-R12 Gradient',ipdd
         do isymaj = 1,nsym
            do isyma = 1,nsym
               isymi =isyma
               write(lupri,*) 'symmetry block',imataj(isyma,isymi)
               call output(work(kajij+imataj(isyma,isymi)),
     &              1,nvir0(isyma),1,nrhf(isymi),
     &              nvir0(isyma),nrhf(isymi),1,lupri)

            enddo
         enddo
      end if

      if (ipdd .eq. 2) then
        luvajkl = -1
        call gpopen(luvajkl,'CCR12YAJIJ','unknown',' ','unformatted',
     &               idummy,.false.)
        write(luvajkl) (work(kajij+i-1), i = 1,ncvai)
        call gpclose(luvajkl,'KEEP')
      else if (ipdd .eq. 3) then 
        call dzero(work(res),ncvai)
        call r12gradap(ipdd,work(res),orb,work(kend2),v2am,lwrk2)
        call dscal(ncvai,2.0D0,work(res),1)
        call daxpy(ncvai,one,work(kajij),1,work(res),1)
        luvajkl = -1
        call gpopen(luvajkl,'CCR12ZAJIJ','unknown',' ','unformatted',
     &               idummy,.false.)
        write(luvajkl) (work(res+i-1), i = 1,ncvai)
        call gpclose(luvajkl,'KEEP')
      else if (ipdd .eq. 5) then
        call dzero(work(res),ncvai)
        call r12gradap(ipdd,work(res),orb,work(kend2),v2am,lwrk2)
        call dscal(ncvai,2.00D0,work(res),1)
        call daxpy(ncvai,one,work(kajij),1,work(res),1)
        luvajkl = -1
        call gpopen(luvajkl,'CCR12XAJIJ','unknown',' ','unformatted',
     &               idummy,.false.)
        write(luvajkl) (work(res+i-1), i = 1,ncvai)
        call gpclose(luvajkl,'KEEP')
 
      endif
c
c
      call qexit('cc_r12grad')
      end
c------------------------------------------------------------------------
C=====================================================================*
C  /* Deck r12gradap */
      subroutine r12gradap(ipdd,res,orb,work,v2am,lwork)
      implicit none
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>
#include <r12int.h>
#include <priunit.h>
#include <ccfro.h>
#include <ccsdinp.h>
      integer kccr12,krsing,krtrip,kctil,kend1,lwrk1,nrhftria,n2 
      integer lwork,lu43,idxkl,isymv,isymi,isymj,isymk,idxij
      integer isyml,isymkl,isymij,idxji,idxklij,idxijkl,idxklji
      integer kxajkl,lunit,luvajkl,ntklaj,isyma,isymaj,kend2
      integer nmatan(8),isyman,isymn,ntotij,ntotan,nvir0(8)
      integer lctil,klijan,ntota,ntoti,ntotn,lcvai,isymai
      integer kijan,dklij,dklijz,ntotk,koffc,lccr12,koffd
      integer idxijlk,idxlk,isymin,kbajkl,ncvai,imatan(8,8),icou3 
      integer lwrk2,lctil1,klij,kcotil,kdotil,koccr,koctil
      integer idxlkji,idxlkij,idxjilk,krklaj,krajkl,ntaj
      integer kxsing,kxtrip,kxxr12,klijmn,lxxr12,ij,al,idx
      integer isym,nvir0t,blajkl,klijz,krajklz,klijanz,kxajklz
      integer ipdd,kxtil,koffk,koffl,koffj,koffi,isymmu,isymmuj
      integer icou4,imataj(8,8),blajklf,kl,aj,ajkl,nrhfst
      integer imatmn(8,8),icou6,nmatkl1(8),imatkl1(8,8),icou5
      integer imatka(8,8),icou7,nmatak(8),nijkl(8)
#if defined (SYS_CRAY)
      real res(*),work(*),orb(*),v2am(*),one,two,zero
#else
      double precision work(*),orb(*),v2am(*),one,zero,two
      double precision res(*)
#endif
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)

      call qenter('r12gradap')

      isymv = 1

      nrhftria = nrhft*(nrhft+1)/2
      n2 = nrhftria * nrhftria

      nrhfst = 0
      do isymi = 1,nsym
         nrhfst = nrhfst + nrhfs(isymi)
      enddo

      kccr12  = 1
      krsing  = kccr12  + ngamsq(1)
      krtrip  = krsing  + n2
      kctil   = krtrip  + n2
      kend1   = kctil  + ngamsq(1)
      lwrk1   = lwork   - kend1

      call dzero(work(krsing),2*n2)
      call dzero(work(kccr12),ngamsq(1))
      call dzero(work(kctil),ngamsq(1))

      if (lwrk1 .lt. 0) then
         call quit('Insufficient work space in r12gradap')
      end if

c     ----------------------------------------
c     read R12 amplitudes
c     ----------------------------------------
      lu43 = -43
      call gpopen(lu43,'CCR12_C','unknown',' ','formatted',
     &                   idummy,.false.)
      read(lu43,'(4e30.20)') (work(krsing+i), i = 0, n2-1 )
      read(lu43,'(4e30.20)') (work(krtrip+i), i = 0, n2-1 )
      call gpclose(lu43,'KEEP')


      call dscal(nrhftria*nrhftria,0.5d0,work(krsing),1)
      call dscal(nrhftria*nrhftria,0.5d0,work(krtrip),1)

      call cc_r12vpack(work(kccr12),isymv,work(krsing),work(krtrip),
     &     nrhf,nrhfa,nijkl)

c     ---------------------------------
C     get c_tilde = 2*c_ij^kl - c_ji^kl
c     ---------------------------------

      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  do isymi = 1, nsym
                     isymj = muld2h(isymij,isymi)
                     do j = 1, nrhf(isymj)
                        do i = 1, nrhf(isymi)
                           idxij = imatij(isymi,isymj)+
     &                          nrhf(isymi)*(j-1) + i
                           idxji = imatij(isymj,isymi)+
     &                          nrhf(isymj)*(i-1) + j
                           idxklij =  igamsq(isymij,isymkl)+
     &                          nmatij(isymkl)*(idxij-1)+idxkl
                           idxklji =igamsq(isymkl,isymij)+
     &                          nmatij(isymkl)*(idxji-1)+idxkl
                           idxijkl = igamsq(isymkl,isymij)+
     &                          nmatij(isymij)*(idxkl-1) +idxij
                           work(kctil-1+idxklij) =
     &                          two*work(kccr12-1+idxklij) -
     &                          work(kccr12-1+idxklji)
                        end do
                     end do

                  end do
               end do
            end do
         end do
      end do


      
      ntklaj = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  do isyma =1, nsym
                     isymj  = muld2h(isymaj,isyma)
                     isymi  = isyma
                     ntklaj = ntklaj + nrhf(isymj) *
     &                        (norb1(isyma) - nrhfs(isyma))
                  enddo
               enddo
            enddo
         enddo
      enddo

      do isyman = 1,nsym
         icou3 = 0
         icou4 = 0
         icou5 = 0
         do isyma = 1,nsym
            isymn = muld2h(isyman,isyma)
            nvir0(isyma)  = norb1(isyma) - nrhfs(isyma)
            icou3 = icou3 + nvir0(isyma)*nrhfs(isymn)
            icou4 = icou4 + nrhfs(isyma)*nrhfs(isymn)
            icou5 = icou5 + nvir0(isyma)*nrhf(isymn)
          enddo
          nmatan(isyman) = icou3
          nmatkl1(isyman) = icou4
          nmatak(isyman) = icou5
      enddo

      icou3 = 0
      do isyman = 1,nsym
         icou3 = 0
         icou4 = 0
         icou5 = 0
         icou6 = 0
         icou7 = 0
         do isyma = 1,nsym
            isymn = muld2h(isyman,isyma)
            imatan(isyma,isymn)= icou3
            imataj(isymn,isyma)= icou4
            imatkl1(isyma,isymn)= icou5
            imatmn(isymn,isyma)= icou6
            imatka(isymn,isyma)= icou7
            icou3 = icou3+nrhf(isymn)*(norb1(isyma) - nrhfs(isyma))
            icou4 = icou4+nrhfs(isymn)*(norb1(isyma) - nrhfs(isyma))
            icou5 = icou5+nmatan(isymn)* nmatkl1(isyma)
            icou6 = icou6+nrhfs(isyma)*nrhfs(isymn)
            icou7 = icou7+nmatak(isyma)*nmatij(isymn)
         enddo
      enddo
      do isyman = 1,nsym
         do isyma = 1,nsym
            isymn = muld2h(isyma,isyman)
         enddo
      enddo

      do isymaj = 1,nsym
         isymij = isymaj
         ncvai = 0
         do isyma = 1,nsym
            isymj = muld2h(isymaj,isyma)
            isymi = muld2h(isymij,isymj)
            ncvai = ncvai + (norb1(isyma) - nrhfs(isyma))*nrhf(isymi)
         enddo
      enddo

      ntaj = 0
      do isymkl = 1,nsym
         isymij = isymkl
         do isymj = 1,nsym
            isyma = muld2h(isymij,isymj)
             ntaj = ntaj + nrhf(isymj)
     &            *(norb1(isyma) - nrhfs(isyma))
         enddo
      enddo

      nvir0t = 0
      do isym = 1, nsym
         nvir0(isym) = norb1(isym) - nrhfs(isym)
         nvir0t      = nvir0t + nvir0(isym)
      end do


      kxajkl  = kend1 
      klijan  = kxajkl + ntklaj
      dklij   = klijan + ncvai 
      koctil  = dklij  + ngamsq(1)
      kdotil  = koctil + ngamsq(1)
      klij    = kdotil + ngamsq(1)
      krajkl  = klij   + ngamsq(1)
      kxxr12  = krajkl + nvir0t*nrhfst**3 
      kxsing  = kxxr12 + ngamsq(1) 
      kxtrip  = kxsing + n2
      klijmn  = kxtrip + n2
      blajkl  = klijmn + ngamsq(1) 
      kxtil   = blajkl + nvir0t*nrhfst**3 
      blajklf = kxtil  + ngamsq(1)
      kend2   = blajklf + ntklaj
      lwrk2   = lwork  - kend2      




      call dzero(work(koctil),ngamsq(1))  
       do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  idxlk=imatij(isyml,isymk)+nrhf(isyml)*(k-1)+l
                  do isymi = 1, nsym
                     isymj = muld2h(isymij,isymi)
                     do j = 1, nrhf(isymj)
                        koffj = irhf(isymj)+j
                        do i = 1, nrhf(isymi)
                           ij = (i-1)*nrhf(isymj)+j
                           koffi = irhf(isymi)+i
                           idxij = imatij(isymi,isymj)+
     &                          nrhf(isymi)*(j-1) + i
                           idxji = imatij(isymj,isymi)+
     &                          nrhf(isymj)*(i-1) + j
                           idxklij =  igamsq(isymij,isymkl)+
     &                          nmatij(isymkl)*(idxij-1)+idxkl
                           idxklji =igamsq(isymkl,isymij)+
     &                          nmatij(isymkl)*(idxji-1)+idxkl
                           idxijkl = igamsq(isymkl,isymij)+
     &                          nmatij(isymij)*(idxkl-1) +idxij
                           idxijlk = igamsq(isymij,isymkl)+
     &                          nmatij(isymij)*(idxlk-1) +idxij
                           work(koctil-1+idxijlk) =
     &                           work(koctil-1+idxijlk) -two*
     &                          work(kctil-1+idxijlk) * orb(koffi)
                           work(koctil-1+idxijlk) =
     &                          work(koctil-1+idxijlk) -two*
     &                          work(kctil-1+idxijlk) * orb(koffj)
                        end do
                     end do

                  end do
               end do
            end do
         end do
      end do

c     Calculate d_{kl}^{ij} = c_{kl}^{mn}*ctilde_{ij}^{mn}

      call dzero(work(dklij),ngamsq(1))
      call dzero(work(klij),ngamsq(1))
      do isymkl = 1,nsym
         isymij  = isymkl
         lccr12 = kccr12 + igamsq(isymij,isymkl)
         lctil  = koctil + igamsq(isymij,isymkl)
         lctil1  = kctil + igamsq(isymij,isymkl)
         dklijz = dklij + igamsq(isymij,isymkl)
         klijz = klij + igamsq(isymij,isymkl)
         ntoti = max(1,nmatij(isymij))
         ntotk = max(1,nmatij(isymkl))
         call dgemm('T','N',nmatij(isymij),nmatij(isymkl),
     &              nmatij(isymkl),one,work(lccr12),
     &              ntotk,work(lctil),ntoti,zero,
     &              work(dklijz),ntotk)
         call dgemm('T','N',nmatij(isymij),nmatij(isymkl),
     &              nmatij(isymkl),one,work(lccr12),
     &              ntotk,work(lctil1),ntoti,zero,
     &              work(klijz),ntotk)
      end do



      call dzero(work(kdotil),ngamsq(1))
      do isymkl = 1, nsym
        isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               koffk = irhf(isymk)+k
               do l = 1, nrhf(isyml)
                  koffl = irhf(isyml)+l
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  idxlk=imatij(isyml,isymk)+nrhf(isyml)*(k-1)+l
                  do isymi = 1, nsym
                     isymj = muld2h(isymij,isymi)
                     do j = 1, nrhf(isymj)
                        koffj = irhf(isymj)+j
                        do i = 1, nrhf(isymi)
                           koffi = irhf(isymi)+i
                           idxij = imatij(isymi,isymj)+
     &                          nrhf(isymi)*(j-1) + i
                           idxji = imatij(isymj,isymi)+
     &                          nrhf(isymj)*(i-1) + j
                           idxklij =  igamsq(isymij,isymkl)+
     &                          nmatij(isymkl)*(idxij-1)+idxkl
                           idxklji =igamsq(isymkl,isymij)+
     &                          nmatij(isymkl)*(idxji-1)+idxkl
                           idxijkl = igamsq(isymkl,isymij)+
     &                          nmatij(isymij)*(idxkl-1) +idxij
                           idxijlk = igamsq(isymkl,isymij)+
     &                          nmatij(isymij)*(idxlk-1) +idxij
                           work(kdotil-1+idxijlk) =
     &                          work(kdotil-1+idxijlk) +
     &                          work(klij-1+idxijlk) * orb(koffi)
                           work(kdotil-1+idxijlk) =
     &                          work(kdotil-1+idxijlk) +
     &                          work(klij-1+idxijlk) * orb(koffj)
                           work(kdotil-1+idxijlk) =
     &                          work(kdotil-1+idxijlk) +
     &                          work(klij-1+idxijlk) * orb(koffk)
                           work(kdotil-1+idxijlk) =
     &                          work(kdotil-1+idxijlk) +
     &                          work(klij-1+idxijlk) * orb(koffl)
                        end do
                     end do

                  end do
               end do
            end do
         end do
      end do

      call daxpy(ngamsq(1),one,work(dklij),1,work(kdotil),1)


      luvajkl = -1
      call dzero(work(kxajkl),ntklaj)
      call gpopen(luvajkl,'DDR12XAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kxajkl+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'KEEP')



      call dzero(work(krajkl),nvir0t*nrhfst**3)
      luvajkl = -1
      call gpopen(luvajkl,'AUXQA12','unknown',' ','formatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(krajkl+i-1), i=1,nvir0t*nrhfst**3)
      call gpclose(luvajkl,'KEEP')

      call dzero(work(blajkl),ntklaj)
      call cc_r12xsort(work(krajkl),work(blajkl),nvir0)

      if (froimp) then
         call dzero(work(blajklf),ntklaj)
         idx = 0
         do isymkl = 1, nsym
            isymaj  = isymkl
            isymij  = isymkl
            do isyml =1, nsym
               isymk = muld2h(isymkl,isyml)
               do k = nrhffr(isymk)+1, nrhfs(isymk)
                  do l = nrhffr(isyml)+1, nrhfs(isyml)
                     do isyma =1, nsym
                        isymj  = muld2h(isymaj,isyma)
                        nvir0(isyma)  = norb1(isyma) - nrhfs(isyma)
                        do j = nrhffr(isymj)+1, nrhfs(isymj)
                           do a = 1, nvir0(isyma)
                              idx = idx + 1
                              kl =imatmn(isymk,isyml)+ 
     &                            (k-1)*nrhfs(isyml) + l
                              aj = imataj(isymj,isyma)+
     &                            (j-1)*nvir0(isyma) + a
                              ajkl = imatkl1(isymaj,isymkl)+
     &                          aj +(kl-1)*nmatan(isymaj)
                              work(blajklf+idx-1) = work(blajkl+ajkl-1)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif 
      

      if (froimp) then
         call daxpy(ntklaj,one,work(blajklf),1,work(kxajkl),1)
      else 
         call daxpy(ntklaj,one,work(blajkl),1,work(kxajkl),1)
      endif

      call dzero(work(klijan),ncvai)
      do isymkl = 1, nsym
         isyman  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhf(isymk)
             do l = 1, nrhf(isyml)
                idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymn  = muld2h(isyman,isyma)
                    isymi  = muld2h(isymij,isymn)
                    nvir0(isyma)  = norb1(isyma) - nrhfs(isyma)
                    ntota = max(1,nvir0(isyma))
                    ntotn = max(1,nrhf(isymn))
                    ntoti = max(1,nrhf(isymi))
                    koffd = kdotil+ igamsq(isymij,isymkl)
     &              +  nmatij(isymij)*(idxkl-1)+imatij(isymi,isymn)
                    klijanz = klijan + imatan(isyma,isymi)
                    call dgemm('N','T',nvir0(isyma),nrhf(isymi),
     &                         nrhf(isymn),one,work(kxajkl),
     &                         ntota,work(koffd),ntoti,one,
     &                         work(klijanz),ntota)
                    kxajkl = kxajkl + nvir0(isyma)*nrhf(isymn)
                    blajklf = blajklf + nvir0(isyma)*nrhf(isymn)

                 enddo
              enddo
           enddo
         enddo
      enddo  
      


      call r12gradaxd(ipdd,work(kccr12),work(kctil),work(klij),
     &                work(kend2),
     &                v2am,res,lwrk2)
      call daxpy(ncvai,one,work(klijan),1,res(1),1)


      call qexit('r12gradap')
      end
c--------------------------------------------------------------------------
C=====================================================================*
C  /* Deck r12gradaxd */
      subroutine r12gradaxd(ipdd,kccr12,kctil,dklij,work,v2am,res,
     &                      lwork)
      implicit none
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>
#include <r12int.h>
#include <priunit.h>
#include <ccfro.h>
      character*8 CCR12YIJ
      integer krsing,krtrip,kend1,lwrk1,nrhftria,n2
      integer lwork,lu43,idxkl,isymv,isymi,isymj,isymk,idxij
      integer isyml,isymkl,isymij,idxji,idxklij,idxijkl,idxklji
      integer isymia,luvajkl,ntklaj,isymaj,kend2
      integer nmatan(8),isyman,isyma,isymn,ntotij,ntotan,nvir0(8)
      integer lctil,klijan,ntota,ntoti,ntotn,lcvai,isymai
      integer kijan,ntotk,koffc,lccr12,koffd
      integer idxijlk,idxlk,isymin,kbajkl,ncvai,imatan(8,8),icou3
      integer lwrk2,lctil1,kcotil,kdotil,koccr,koctil
 
      integer kxijkl,kaklij,nak,isymik,nai,naiki,nkj,naikj,isymjk
      integer ntotmu,dklijz,eklij,eklijz,lxxr12,kxxr12,kxsing,kxtrip
      integer kiakj,nia,njk,niajk,idx,ntotj,kgiakj,kgaijk,kaxdia 
      integer fklij,fklijz,kxcijkl,kxcijklz,ntklij,ntaj
      integer naj,nki,najki,nji,naijk,nakji,kxcijklt,nka
      integer rklij,koffh,dijkl,hklij,kaxcia,kaxdai,kiajk,kgaikj
      integer dlkmnz,dlkmn,kctil1,kcxijklt,kcxijkl,krcc12,lunit
      integer kcxklmnz,kcxmnklz,kcxmnkl,kcxklmn,kxxr12n,idxjikl,idxjilk
      integer aijkl,hijkl,hklijz,lxxr12n,dijklz,eijkl,eijklz
      integer njaik,nik,nja,najik,nikja,aklij,kctil2,aijklz,lrcc12,bijkl
      integer isymja,isymka,ntoto,isymoj,isymo,idxlkij,kres,nikaj
      integer kiakjz,aj,ajkl,kl,isymkj,kiajkz,rklijz,index
      integer ipdd,resb,kxtil,kresz,icoun4,nmataj(8),ncvij,nijkl(8)
#if defined (SYS_CRAY)
      real work(*),dklij(*),v2am(*),kctil(*),kccr12(*),one,two,zero
      real res(*) 
#else
      double precision work(*),dklij(*),v2am(*),kctil(*),one,zero,two
      double precision kccr12(*),res(*)
#endif
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)

      index(i,j) = max(i,j)*(max(i,j)-3)/2 +i+j

      call qenter('r12gradaxd')

C     Calculate some offsets

      isymv = 1

      nrhftria = nrhft*(nrhft+1)/2
      n2 = nrhftria * nrhftria



      do isymaj = 1,nsym
         isymij = isymaj
         ncvij = 0
         do isyma = 1,nsym
            isymj = muld2h(isymaj,isyma)
            isymi = muld2h(isymij,isymj)
            ncvij = ncvij +  nrhf(isyma)*nrhf(isymi)
         enddo
      enddo
   
       


      do isymaj = 1,nsym
         isymij = isymaj
        ntaj = 0
         do isyma = 1,nsym
            isymj = muld2h(isymaj,isyma)
            isymi = muld2h(isymij,isymj)
            ntaj = ntaj + (norb1(isyma) - nrhfs(isyma))*nrhf(isymi)
         enddo
      enddo


      ntklaj = 0
      do isymkl = 1,nsym
         isymaj = isymkl
         do isymj = 1,nsym
            isyma = muld2h(isymaj,isymj)
            ntklaj = ntklaj + nmatij(isymkl)*nrhf(isymj)
     &               *(norb1(isyma) - nrhfs(isyma))
         enddo
      enddo



      kxxr12  = 1
      kxsing  = kxxr12  + ngamsq(1)
      kxtrip  = kxsing  + n2 
      eklij   = kxtrip  + n2
      kgiakj  = eklij   + ncvij 
      kgaijk  = kgiakj  + ntklaj 
      kxcijkl = kgaijk  + ntklaj 
      fklij   = kxcijkl + ngamsq(1) 
      kiajk   = fklij   + ncvij 
      kgaikj  = kiajk   + ntklaj
      rklij   = kgaikj  + ntklaj 
      kcxijkl = rklij   + ncvij 
      eijkl   = kcxijkl + ngamsq(1) 
      kxtil   = eijkl   + ncvij 
      resb    = kxtil   + ngamsq(1)
      kend2   = resb    + ntaj 
      lwrk2   = lwork   - kend2 
      
C     Read matrix X

      lu43 = -43
      call gpopen(lu43,'CCR12_XP','unknown',' ','formatted',
     &                   idummy,.false.)
      read(lu43,*) i
      read(lu43,'(4e30.20)') (work(kxsing+i), i = 0, n2-1 )
      read(lu43,'(4e30.20)') (work(kxtrip+i), i = 0, n2-1 )
      if (ipdd.eq.3) then 
          call gpclose(lu43,'KEEP')
      else 
          call gpclose(lu43,'DELETE')
      endif
      call dscal(n2,2d0,work(kxsing),1)
      call dscal(n2,2d0/3d0,work(kxtrip),1)

      call cc_r12vpack(work(kxxr12),isymv,work(kxsing),work(kxtrip),
     &     nrhf,nrhfa,nijkl)


C     Calculate XC_{ij}^{kl} = X_{ij}^{mn} * Ctilde_{mn}^{kl}   

      call dzero(work(kxcijkl),ngamsq(1))
      do isymkl = 1,nsym
         isymij  = isymkl
         lxxr12   = kxxr12 + igamsq(isymij,isymkl)
         lctil    = 1 + igamsq(isymij,isymkl)
         kxcijklz = kxcijkl + igamsq(isymij,isymkl)
         ntoti = max(1,nmatij(isymij))
         ntotk = max(1,nmatij(isymkl))
         call dgemm('N','T',nmatij(isymij),nmatij(isymkl),
     &              nmatij(isymkl),one,work(lxxr12),
     &              ntotk,kctil(lctil),ntoti,zero,
     &              work(kxcijklz),ntotk)
      end do
c     ---------------------------------
C     Transpone CX_{ij}^{kl}
c     ---------------------------------

      call dzero(work(kcxijkl),ngamsq(1))
      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  do isymi = 1, nsym
                     isymj = muld2h(isymij,isymi)
                     do j = 1, nrhf(isymj)
                        do i = 1, nrhf(isymi)
                           idxij = imatij(isymi,isymj)+
     &                          nrhf(isymi)*(j-1) + i
                           idxklij =  igamsq(isymkl,isymij)+
     &                          nmatij(isymkl)*(idxij-1)+idxkl
                           idxijkl = igamsq(isymij,isymkl)+
     &                          nmatij(isymij)*(idxkl-1) +idxij
                           work(kcxijkl-1+idxijkl) =
     &                         work(kxcijkl-1+idxklij)
                        end do
                     end do

                  end do
               end do
            end do
         end do
      end do



      call dzero(work(eklij),ncvij)
      call dzero(work(eijkl),ncvij)
      call dzero(work(fklij),ncvij)
      do isymkl = 1, nsym
         isymoj  = isymkl
         isymij = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhf(isymk)
             do l = 1, nrhf(isyml)
                idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                 do isymo =1, nsym
                    isymj  = muld2h(isymoj,isymo)
                    isymi  = muld2h(1,isymo) 

C     Calculate X_{ol}^{ij} * D_{ij}^{kl}

                    koffd  = 1 + igamsq(isymij,isymkl)
     &                 + nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    lxxr12 = kxxr12 + igamsq(isymij,isymkl)
     &                 + nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    ntotj = max(1,nrhf(isymj))
                    ntoti = max(1,nrhf(isymi))
                    ntoto = max(1,nrhf(isymo)) 
                    eklijz = eklij + imatij(isymo,isymi) 
                    call dgemm('N','T',nrhf(isymo),nrhf(isymi),
     &                    nrhf(isymj),one,work(lxxr12),
     &                    ntoto,dklij(koffd),ntoti,one,
     &                    work(eklijz),ntoto)
                    eijklz = eijkl + imatij(isymi,isymo)
                    call dgemm('N','T',nrhf(isymi),nrhf(isymo),
     &                    nrhf(isymj),one,dklij(koffd),
     &                    ntoti,work(lxxr12),ntoto,one,
     &                    work(eijklz),ntoti)

C     Calculate   F_{}^{} = C_{ij}^{kl} * XC_{ij}^{ol} 

                    lccr12 = 1+ igamsq(isymij,isymkl)
     &                 + nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    kcxijklt = kcxijkl + igamsq(isymij,isymkl)
     &                 + nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    fklijz = fklij + imatij(isymo,isymi) 
                    call dgemm('N','T',nrhf(isymi),nrhf(isymo),
     &                    nrhf(isymj),one,kccr12(lccr12),
     &                    ntoti,work(kcxijklt),ntoti,one,
     &                    work(fklijz),ntoti)

                 enddo
              enddo
            enddo
         enddo
      end do

      call dzero(work(rklij),ncvij)
      call daxpy(ncvij,one,work(eklij),1,work(rklij),1)
      call daxpy(ncvij,one,work(eijkl),1,work(rklij),1)
      call daxpy(ncvij,-two,work(fklij),1,work(rklij),1)
      call dscal(ncvij,two,work(rklij),1)
      lunit = -1
      if (ipdd .eq. 3) then 
         call gpopen(lunit,'CCR12YIJ','unknown',' ','unformatted',
     &               idummy,.false.)
      else if (ipdd .eq. 5) then
         call gpopen(lunit,'CCR12YIJB','unknown',' ','unformatted',
     &               idummy,.false.)
      endif
      write(lunit)(work(rklij+i-1),i=1,ncvij)
      call gpclose(lunit,'KEEP') 


      call dzero(work(kgiakj),ntklaj)
      call dzero(work(kgaijk),ntklaj)
      call dzero(work(kgaikj),ntklaj)
      call dzero(work(kiajk),ntklaj)

      idx = 0
      do isymka = 1,nsym
         isymij = isymka
         do isyma = 1, nsym
            nvir0(isyma) = norb1(isyma) - nrhfs(isyma)
            isymk = muld2h(isymka,isyma)
            do isymj = 1,nsym
               isymi = muld2h(isymij,isymj)
               do i = 1,nrhf(isymi)
                  do j = 1,nrhf(isymj)
                     isymia = muld2h(isymi,isyma)
                     isymkj = muld2h(isymk,isymj)
                     isymja = muld2h(isymj,isyma)
                     isymik = muld2h(isymi,isymk)
                     do k = 1, nrhf(isymk)
                        do a = 1,nvir0(isyma)
                          idx = idx + 1 
                          nai = it1am(isyma,isymi)
     *                        + nvir(isyma)*(i-1) 
     &                         + a + nrhfs(isyma)
                          njk = it1am(isymj,isymk)
     *                         + nvir(isymj)*(k-1) + j+nrhffr(isymj)
                          naijk = it2am(isymia,isymkj)
     *                           + index(nai,njk)
                          nak = it1am(isyma,isymk)
     *                         + nvir(isyma)*(k-1) 
     &                         + a + nrhfs(isyma)
                          nji = it1am(isymj,isymi)
     *                          + nvir(isymj)*(i-1) + j+nrhffr(isymj)
                          nakji = it2am(isymka,isymij)
     *                            + index(nak,nji)
                          naj = it1am(isyma,isymj)
     *                          + nvir(isyma)*(j-1) 
     &                          + a + nrhfs(isyma)
                          nki = it1am(isymk,isymi)
     *                          + nvir(isymk)*(i-1) + k +nrhffr(isymk)
                          nik = it1am(isymi,isymk)
     *                          + nvir(isymi)*(k-1) + i
                          najki = it2am(isymja,isymik)
     *                            + index(naj,nki)
                          work(kgiakj+idx-1) = v2am(naijk)
                          work(kgaijk+idx-1) = v2am(nakji)
                          work(kgaikj+idx-1) = v2am(najki)
                          work(kiajk+idx-1)  = 4.0d0*work(kgaijk+idx-1)
     &                        - work(kgaikj+idx-1) - work(kgiakj+idx-1)
                      enddo
                    enddo
                  enddo
               enddo   
            enddo 
         enddo
      enddo


   
      do isymi = 1,nsym  
          nvir0(isymi) = norb1(isymi)-nrhfs(isymi)
          nmataj(isymi) = nvir0(isymi)*nrhf(isymi)
      enddo



      call dzero(res(1),ntaj)
      kres = 1
      do isymij = 1,nsym
         isymkl = isymij
         isymaj = isymij
         do isyma =1, nsym
            isymj  = muld2h(isymaj,isyma)
            isymi  = muld2h(1,isyma)
            nvir0(isyma) = norb1(isyma)-nrhfs(isyma)
            ntotk = max(1,nvir0(isyma)*nrhf(isymj)) 
            ntoti = max(1,nmatij(isymkl)) 
            call dgemm('N','N',nrhf(isymj)*nvir0(isyma),1,
     &                 nmatij(isymkl),one,work(kiajk),ntotk,
     &                 work(rklij),ntoti,zero,res(kres),ntotk) 
            kres = kres + nvir0(isyma)*nrhf(isymj)
            kiajk = kiajk + nrhf(isymj)*nvir0(isyma)
     &                      * nmatij(isymkl) 
         enddo
      enddo 

      call dscal(ntaj,0.5d0,res(1),1) 


      if (ipdd .eq. 5) then
         call dscal(ntaj,1.0d0,res(1),1) 
         call dzero(work(resb),ntaj)
         call r12gradb(work(resb),dklij,work(kend2),ntklaj,ntaj,lwrk2)
         
         call daxpy(ntaj,one,work(resb),1,res(1),1)
      endif
      call qexit('r12gradaxd')
      end
c--------------------------------------------------------------------------
C=====================================================================*
C  /* Deck r12gradb */
      subroutine r12gradb(res,dklij,work,ntklaj,ncvai,lwork)
      implicit none
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>
#include <r12int.h>
#include <priunit.h>
#include <ccfro.h>
      logical locdbg
      parameter (locdbg = .false.) 
      integer kqajkl,kqijal,kuijal,krajkl,kend1
      integer isym,luvajkl,nvir0t,lwork,lwrk1,ntklaj
      integer icou3,imataj(8,8),koffd,ksajkl,nvir0(8),ncvai
      integer kcqia,kcqai,lcqia,lcqai,isymai,isyma,isymi
      integer isymaj,isymkl,isymij,isyml,isymk,isymj,idxkl,ntotj
      integer isymv,idxlk,idxklij,idxlkij,idxij,idxji
      integer idxjikl,idxijkl,idxijlk,koff,ntota,ntoti
      integer kuajkl,ktajkl,ktijal,kcuia,kcuai,lcuai,lcuia
      integer kssijal,ksrajkl,kstijal,idxklji,koff1,kdjikl
      integer kuijals,kqijals,kssajkl

      integer kccr12,krsing,krtrip,kctil,lccr12,lctil
      integer nrhftria,n2,lu43,ntotk
#if defined (SYS_CRAY)
      real dklij(*),res(*),work(*),one,two,zero
#else
      double precision dklij(*),res(*),work(*),one,zero,two
#endif
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)

      call qenter('r12gradb')

      isymv = 1
      
      nrhftria = nrhft*(nrhft+1)/2
      n2 = nrhftria * nrhftria


      nvir0t = 0
      do isym = 1, nsym
         nvir0(isym) = norb1(isym) - nrhfs(isym)
         nvir0t      = nvir0t + nvir0(isym)
      end do

      icou3 = 0
      do isymai = 1,nsym
         icou3 = 0
         do isyma = 1,nsym
            isymi = muld2h(isymai,isyma)
            imataj(isyma,isymi)= icou3
            icou3 = icou3+nrhf(isymi)*(norb1(isyma) - nrhfs(isyma))
         enddo
      enddo

  
      kqajkl = 1
      kqijal = kqajkl + ntklaj
      kuijal = kqijal + ntklaj
      krajkl = kuijal + ntklaj
      ksajkl = krajkl + nvir0t*nrhft**3 
      kcqai  = ksajkl + nvir0t*nrhft**3
      kcqia  = kcqai  + ncvai
      ktijal = kcqia  + ncvai
      kcuai  = ktijal + nvir0t*nrhft**3 
      kuajkl = kcuai  + ncvai
      kcuia  = kuajkl + ntklaj
      kstijal= kcuia  + ncvai
      ksrajkl= kstijal+ ntklaj
      kssijal= ksrajkl+ ntklaj 
      kssajkl= kssijal+ ntklaj
      kuijals= kssajkl+ ntklaj 
      kqijals= kuijals+ ntklaj
      kend1  = kqijals+ ntklaj 
      lwrk1  = lwork  - kend1


      luvajkl = -1
      call dzero(work(kqajkl),ntklaj)
      call gpopen(luvajkl,'DDR12QAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kqajkl+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'DELETE')


      luvajkl = -1
      call dzero(work(kqijal),ntklaj)
      call gpopen(luvajkl,'DDR12QIJAL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kqijal+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'DELETE')
  
      call dzero(work(kqijals),ntklaj)

      luvajkl = -1
      call dzero(work(kuijal),ntklaj)
      call gpopen(luvajkl,'DDR12UIJAL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kuijal+i-1), i = 1,ntklaj)
          call gpclose(luvajkl,'DELETE')

      call dzero(work(kuijals),ntklaj)

  
      luvajkl = -1
      call dzero(work(kuajkl),ntklaj)
      call gpopen(luvajkl,'DDR12UAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kuajkl+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'DELETE')

      call dzero(work(krajkl),nvir0t*nrhft**3)
      luvajkl = -1
      call gpopen(luvajkl,'AUXQSA12','unknown',' ','formatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(krajkl+i-1), i=1,nvir0t*nrhft**3)
      call gpclose(luvajkl,'DELETE')

      call dzero(work(ksrajkl),ntklaj)  
      call cc_r12xsort2(work(krajkl),work(ksrajkl),nvir0)


      call dzero(work(ksajkl),nvir0t*nrhft**3)
      luvajkl = -1
      call gpopen(luvajkl,'AUXQSAJ12','unknown',' ','formatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(ksajkl+i-1), i=1,nvir0t*nrhft**3)
      call gpclose(luvajkl,'DELETE')

      call dzero(work(kssijal),ntklaj)  
      call cc_r12xsort2(work(ksajkl),work(kssijal),nvir0) 

      call dzero(work(kssajkl),ntklaj)  
      call cc_r12xsort3(work(ksajkl),work(kssajkl),nvir0) 



      call dzero(work(ktijal),nvir0t*nrhft**3)
      luvajkl = -1
      call gpopen(luvajkl,'AUXQSAI12','unknown',' ','formatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(ktijal+i-1), i=1,nvir0t*nrhft**3)
      call gpclose(luvajkl,'DELETE')

      call dzero(work(kstijal),ntklaj)  
      call cc_r12xsort3(work(ktijal),work(kstijal),nvir0)


      call daxpy(ntklaj,one,work(kssijal),1,work(kqijal),1)
      call daxpy(ntklaj,one,work(kssajkl),1,work(kuijal),1)
      call daxpy(ntklaj,one,work(kstijal),1,work(kuajkl),1)
      call daxpy(ntklaj,one,work(ksrajkl),1,work(kqajkl),1)



      call dzero(work(kcqai),ncvai)
      lcqai = kcqai
      call dzero(work(kcqia),ncvai)
      lcqia = kcqia
      call dzero(work(kcuai),ncvai)
      lcuai = kcuai
      call dzero(work(kcuia),ncvai)
      lcuia = kcuia

c     Calculate Q_{aj}^{kl}*D_{kl}^{ij}


      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
           do k = 1, nrhf(isymk)
             do l = 1, nrhf(isyml)
                idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    isymi  = muld2h(isymij,isymj)
                    nvir0(isyma)  = norb1(isyma) - nrhfs(isyma)
                    koffd  = 1 + igamsq(isymij,isymkl)
     &                 + nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    ntota = max(1,nvir0(isyma))
                    ntotj = max(1,nrhf(isymj))
                    ntoti = max(1,nrhf(isymi))
                    lcqai  = kcqai + imataj(isyma,isymi)
                    call dgemm('N','T',nvir0(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kqajkl),
     &                    ntota,dklij(koffd),ntoti,one,
     &                    work(lcqai),ntota)

                    kqajkl = kqajkl + nvir0(isyma)*nrhf(isymj)
                    ksrajkl = ksrajkl + nvir0(isyma)*nrhf(isymj)

                    ntota = max(1,nvir0(isyma))
                    ntotj = max(1,nrhf(isymj))
                    ntoti = max(1,nrhf(isymi))
                    lcqia  = kcqia + imataj(isyma,isymi)
                    call dgemm('N','T',nvir0(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kqijal),
     &                    ntota,dklij(koffd),ntoti,one,
     &                    work(lcqia),ntota)

                    kqijal = kqijal + nvir0(isyma)*nrhf(isymj)

                    lcuai  = kcuai + imataj(isyma,isymi)
                    call dgemm('N','T',nvir0(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kuijal),
     &                    ntota,dklij(koffd),ntoti,one,
     &                    work(lcuai),ntota)

                    kuijal = kuijal + nvir0(isyma)*nrhf(isymj)
                    kssijal = kssijal + nvir0(isyma)*nrhf(isymj)
                    kssajkl = kssajkl + nvir0(isyma)*nrhf(isymj)

                    lcuia  = kcuia + imataj(isyma,isymi)
                    call dgemm('N','T',nvir0(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kuajkl),
     &                    ntota,dklij(koffd),ntoti,one,
     &                    work(lcuia),ntota)

                    kuajkl  = kuajkl + nvir0(isyma)*nrhf(isymj)
                    kstijal = kstijal + nvir0(isyma)*nrhf(isymj)
                    ktijal  = ktijal + nvir0(isyma)*nrhf(isymj)

                 enddo
              enddo
           enddo
         enddo
      enddo 


      call r12gradbd(res,dklij,work(kend1),lwrk1)

      call daxpy(ncvai,one,work(kcqai),1,res,1)
      call daxpy(ncvai,one,work(kcuai),1,res,1)
      call daxpy(ncvai,one,work(kcqia),1,res,1)
      call daxpy(ncvai,one,work(kcuia),1,res,1)

      if (locdbg) then
         write(lupri,*) 'Ergebnisse r12gradb'
         do isymaj = 1,nsym
            do isyma = 1,nsym
               isymi =isyma
               write(lupri,*) 'symmetry block',isyma
               call output(res(1+imataj(isyma,isymi)),
     &                 1,nvir0(isyma),1,nrhf(isymi),
     &                 nvir0(isyma),nrhf(isymi),1,lupri)
            enddo
         enddo
       endif


      call qexit('r12gradb')
      end
c------------------------------------------------------------------------
C=====================================================================*
C  /* Deck r12gradbd */
      subroutine r12gradbd(res,dklij,work,lwork)
      implicit none
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>
#include <r12int.h>
#include <priunit.h>
#include <ccfro.h>
#include <ccsdinp.h>
      integer lxdmmu,nr1bast,isymak,isymk,isyma,nr1bas(8)
      integer lwork,kxmujkl,kxdmmu,kend1,nmum,luvajkl
      integer isymkl,isymij,isymi,isymj,idxkl,koffd,isymaj
      integer isymai,icou3,imatmuj(8,8),ntota,ntoti,ntotj,isyml
      integer smujkl,koffc,lusifc,lwrk1,kcmo,ntota1,kxdnumu
      integer kxpjkl,ntklpj,isymmua,iglmrv(8,8),isymmu
      integer kxajkl,ntotmu,isymmuj,kxpjklz,icount2,icou2,icount7
      integer nr1xbas(8),koff,idxji,idxij,kzmunu,kxdpmu,ir1bas(8,8)
      integer nrhftria,kxsing,kxtrip,n2,isymv,lu43,nxajkl(8),kcmoa
      integer isym,noffset,kspjkl,nvir0t,nvir0(8),ntklaj,krajkl
      integer imbas1(8),kdens,lxdnumu,icou4,imatmunu(8,8),lxdpmu
      integer nrhfst,ksqijklz,ksqjkl,ksajkl,sajkl,icoun1,norb1t
      integer imataj(8,8),idx,kl,aj,ajkl,kspjklf,ksajklf
      integer imatmn(8,8),icou6,nmatkl1(8),imatkl1(8,8),icou5
      integer nmatan(8),icou7
      integer imatoi(8,8),icou8,isyman,isymn
      integer iglmrv1(8,8),noffset1,icou9,nrhffrt
      integer ksijklf,ij,ijkl,ksijkl,sijklf
#if defined (SYS_CRAY)
      real res(*),work(*),dklij(*),zero,one,two
#else
      double precision res(*),work(*),dklij(*),zero,one,two
#endif
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)

      call qenter('r12gradbd')

      isymv = 1

      nrhftria = nrhft*(nrhft+1)/2
      n2 = nrhftria * nrhftria


      nvir0t = 0
      norb1t = 0
      nrhfst = 0
      nrhffrt = 0
      do isymi = 1, nsym
         nvir0(isymi) = norb1(isymi) - nrhfs(isymi)
         nvir0t      = nvir0t + nvir0(isymi)
         norb1t      = norb1t + norb1(isymi)
         nrhfst      = nrhfst + nrhfs(isymi)
         nrhffrt     = nrhffrt+ nrhffr(isymi)
      end do


      do isymmua = 1,nsym
         icou2 = 0
         icou3 = 0
         do isyma = 1,nsym
            isymmu = muld2h(isymmua,isyma)
            iglmrv(isymmu,isyma) = icou2
            iglmrv1(isymmu,isyma) = icou3
            icou2 = icou2 + nbas(isymmu)*(nvir(isyma))
            icou3 = icou3 + nbas(isymmu)*(nvir(isyma)-nrhffr(isyma))
         enddo
      enddo

      do isyman = 1,nsym
         icou3 = 0
         icou4 = 0
         icou7 = 0
         nmatkl1(isyman) = 0
         do isyma = 1,nsym
            isymn = muld2h(isyman,isyma)
            nvir0(isyma)  = norb1(isyma) - nrhfs(isyma)
            icou3 = icou3 + nvir0(isyma)*nrhfs(isymn)
            icou4 = icou4 + nrhfs(isyma)*nrhfs(isymn)
          enddo
          nmatan(isyman) = icou3
          nmatkl1(isyman) = icou4
      enddo

      do isyman = 1,nsym
         icou4 = 0
         icou5 = 0
         icou6 = 0
         icou8 = 0
         do isyma = 1,nsym
            isymn = muld2h(isyman,isyma)
            imataj(isyma,isymn)= icou4
            imatkl1(isyma,isymn)= icou5
            imatmn(isyma,isymn)= icou6
            imatoi(isyma,isymn)= icou8
            icou4 = icou4+nrhfs(isymn)*(norb1(isyma) - nrhfs(isyma))
            icou5 = icou5+nmatan(isymn)* nmatkl1(isyma)
            icou6 = icou6+nrhfs(isyma)*nrhfs(isymn)
            icou8 = icou8+nmatkl1(isyma)*nmatkl1(isymn)
         enddo
      enddo



      do isymak = 1, nsym
        nr1bas(isymak) = 0
        nr1xbas(isymak) = 0
        do isymk = 1, nsym
           isyma = muld2h(isymak,isymk)
           nr1bas(isymak)  = nr1bas(isymak) +mbas1(isyma)*nrhf(isymk)
           nr1xbas(isymak) = nr1xbas(isymak)+mbas2(isyma)*nrhf(isymk)
        end do
      end do
 
      nr1bast = 0
      do isyma = 1,nsym
         nr1bast = nr1bast + nr1bas(isyma)
      enddo

      ntklpj = 0
      do isymkl = 1,nsym
         isymaj = isymkl
         do isymj = 1,nsym
            isyma = muld2h(isymaj,isymj)
            ntklpj = ntklpj + nmatij(isymkl)*nrhf(isymj)*norb1(isyma)
         enddo
      enddo
      


      do isymak = 1, nsym
         nvajkl(isymak) = 0
         icount7        = 0
         icount2 = 0
         do isymk = 1, nsym
            isyma = muld2h(isymak,isymk)
            ivajkl(isyma,isymk) = icount7
            nvajkl(isymak) = nvajkl(isymak)+ nr1bas(isyma)*nmatij(isymk)
            icount7 = icount7 + nr1bas(isyma)*nmatij(isymk)
            ir1bas(isyma,isymk)  = icount2
            icount2 = icount2 + nrhf(isymk)*mbas1(isyma)
         end do
      end do

      icou3 = 0
      do isymai = 1,nsym
         icou3 = 0
         icou4 = 0
         do isyma = 1,nsym
            isymi = muld2h(isyma,isymai)
            imatmuj(isyma,isymi)= icou3
            icou3 = icou3+nrhf(isymi)*norb1(isyma)
            imatmunu(isymi,isyma) = icou4
            icou4 = icou4+mbas1(isymi)*mbas1(isyma) 
         enddo
      enddo

      ntklaj = 0
      do isymkl = 1, nsym
         isymaj  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  do isyma =1, nsym
                     isymj  = muld2h(isymaj,isyma)
                     isymi  = isyma
                     ntklaj = ntklaj + nrhf(isymj) *
     &                        (norb1(isyma) - nrhfs(isyma))
                  enddo
               enddo
            enddo
         enddo
      enddo



      kxdmmu  = 1
      kxdnumu = kxdmmu  + nr1bast 
      kcmo    = kxdnumu + 10*mbas1t*mbas1t 
      smujkl  = kcmo    + nlamds
      kxpjkl  = smujkl  + nrhfst**4
      kxdpmu  = kxpjkl  + ntklpj
      kzmunu  = kxdpmu  + mbas1t*nrhft  
      kcmoa   = kzmunu  + n2basx
      kxmujkl = kcmoa   + nlamds
      kspjkl  = kxmujkl + nvajkl(1)
      kdens   = kspjkl  + ntklpj   
      sajkl   = kdens   + n2basx
      ksajkl  = sajkl   + nvir0t*nrhfst**3
      ksqjkl  = ksajkl  + nvir0t*nrhfst**3 
      ksajklf = ksqjkl  + ntklpj
      kspjklf = ksajklf + nvir0t*nrhfst**3 
      ksijklf = kspjklf + nrhfst**4 
      ksijkl  = ksijklf + nrhfst*nrhfst**3
      kend1   = ksijkl  + nrhfst*nrhfst**3 
      lwrk1   = lwork   - kend1

      
      call dzero(work(kxmujkl),nvajkl(1))
      luvajkl = -1
      call gpopen(luvajkl,'CCR12QAJKL','unknown',' ','unformatted',
     &               idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kxmujkl+i-1), i = 1,nvajkl(1))
      call gpclose(luvajkl,'KEEP')
c
      call dzero(work(smujkl),nrhfst**4)
      luvajkl = -1

      call gpopen(luvajkl,'AUXQ12','unknown',' ','formatted',
     &               idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(smujkl+i-1), i = 1,
     &                        nrhfst**4)
      call gpclose(luvajkl,'KEEP')
c
      luvajkl = -1
      call gpopen(luvajkl,'AUXQA12','unknown',' ','formatted',
     &               idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(sajkl+i-1), i = 1,
     &                        nvir0t*nrhfst**3)
      call gpclose(luvajkl,'DELETE')



      call dzero(work(kspjkl),ntklpj)
      call dzero(work(ksijkl),nrhfst**4)
      call cc_r12xsort1(work(smujkl),work(kspjkl))
      call cc_r12xsort(work(sajkl),work(ksajkl),nvir0)
      call cc_r12xsort4(work(smujkl),work(ksijkl))

      if (froimp) then
         call dzero(work(ksijklf),nrhfst**4)
         idx = 0
         do isymkl = 1, nsym
            isymij  = isymkl
            do isyml =1, nsym
               isymk = muld2h(isymkl,isyml)
              do k = nrhffr(isymk)+1, nrhfs(isymk)
                 do l = nrhffr(isyml)+1, nrhfs(isyml)
                    idxkl=imatij(isymk,isyml)+nrhfs(isymk)*(l-1)+k
                    do isymi =1, nsym
                       isymj  = muld2h(isymij,isymi)
                       do j = nrhffr(isymj)+1, nrhfs(isymj)
                          do i = 1, nrhffr(isymi)
                            idx = idx + 1
                             kl = imatmn(isymk,isyml)+
     &                             (l-1)*nrhfs(isymk) + k
                             ij = imatmn(isymi,isymj)+
     &                             (j-1)*nrhfs(isymi) + i
                             ijkl = imatoi(isymij,isymkl)+
     &                             (kl-1)*nmatkl1(isymij) + ij
                             work(ksijklf+idx-1) = work(ksijkl+ijkl-1)
                         enddo
                       enddo
                    enddo
                 enddo
              enddo
            enddo
         enddo



         call dzero(work(ksajklf),ntklaj)
         idx = 0
         do isymkl = 1, nsym
            isymaj  = isymkl
            isymij  = isymkl
            do isyml =1, nsym
               isymk = muld2h(isymkl,isyml)
              do k = nrhffr(isymk)+1, nrhfs(isymk)
                do l = nrhffr(isyml)+1, nrhfs(isyml)
                   idxkl=imatij(isymk,isyml)+nrhfs(isymk)*(l-1)+k
                    do isyma =1, nsym
                       isymj  = muld2h(isymaj,isyma)
                       nvir0(isyma)  = norb1(isyma) - nrhfs(isyma)
                       do j = nrhffr(isymj)+1, nrhfs(isymj)
                          do a = 1, nvir0(isyma)
                             idx = idx + 1
                             kl =imatmn(isyml,isymk)+
     &                            (k-1)*nrhfs(isyml) + l
                              aj = imataj(isyma,isymj)+
     &                            (j-1)*nvir0(isyma) + a
                              ajkl = imatkl1(isymaj,isymkl)+
     &                          aj +(kl-1)*nmatan(isymaj)
                             work(ksajklf+idx-1) = work(ksajkl+ajkl-1)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
            enddo
         enddo

        call dzero(work(kspjklf),nrhfst**4)
         idx = 0
         do isymkl = 1, nsym
            isymaj  = isymkl
            isymij  = isymkl
            do isyml =1, nsym
               isymk = muld2h(isymkl,isyml)
              do k = nrhffr(isymk)+1, nrhfs(isymk)
                do l = nrhffr(isyml)+1, nrhfs(isyml)
                    do isyma =1, nsym
                       isymj  = muld2h(isymaj,isyma)
                       nvir0(isyma)  = norb1(isyma) - nrhfs(isyma)
                       do j = nrhffr(isymj)+1, nrhfs(isymj)
                          do a = nrhffr(isyma)+1, nrhfs(isyma)
                             idx = idx + 1
                             kl =imatmn(isyml,isymk)+
     &                            (k-1)*nrhfs(isyml) + l
                             aj = imatmn(isyma,isymj)+
     &                            (j-1)*nrhfs(isyma) + a
                             ajkl = imatoi(isymaj,isymkl)+
     &                          aj +(kl-1)*nmatkl1(isymaj)
                             work(kspjklf+idx-1) = work(kspjkl+ajkl-1)
                          enddo
                       enddo
                    enddo
                 enddo
              enddo
            enddo
         enddo

      endif

c     Read MO-coefficients

      call dzero(work(kcmo),nlamds)
      lusifc = -1
      call gpopen(lusifc,'SIRIFC','OLD',' ','UNFORMATTED',
     &            idummy,.false.)
      rewind(lusifc)
      call mollab(label,lusifc,lupri)
      read(lusifc)
      read(lusifc)
      read(lusifc) (work(kcmo+i-1),i=1,nlamds)
      call gpclose(lusifc,'KEEP')

!     reorder to occupied, symmetry; virtual, symmetry
!             from symmetry, (occupied,virtual)

      call cmo_reorder(work(kcmo),work(kend1),lwrk1)

      noffset  = 0
      noffset1 = 0
      do isym = 1,nsym
         noffset  = noffset  + nbas(isym)*nrhf(isym)
         noffset1 = noffset1 + nbas(isym)*nrhfs(isym)
      enddo

c     Transform into orbital basis

      call dzero(work(kxpjkl),ntklpj)
      kxpjklz = kxpjkl 
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymmuj = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  do isymmu =1, nsym
                     isyma = muld2h(1,isymmu)
                     isymj  = muld2h(isymmuj,isymmu)
                     koffc = kcmo + noffset+ 
     &                      iglmrv(isymmu,isyma)
                     kxajkl = kxmujkl+ivajkl(isymmuj,isymkl)+
     &                         nr1bas(isymmuj)*(idxkl-1)
     &                         +ir1bas(isymmu,isymj)
                     ntotmu = max(1,mbas1(isymmu))
                     ntota1 = max(1,nbas(isymmu))
                     ntota = max(1,norb1(isyma))
                     ntotj = max(1,nrhf(isymj))
c
                     call dgemm('T','N',norb1(isyma),
     &                 nrhf(isymj),mbas1(isymmu),one,work(koffc),
     &                 ntota1,work(kxajkl),ntotmu,zero,
     &                 work(kxpjkl),ntota)
                     kxpjkl = kxpjkl + norb1(isyma)*nrhf(isymj)
                  end do
               end do
            end do
         end do
      end do


      if (froimp) then
         ksajkl = ksajklf
         kspjkl = kspjklf
         ksijkl = ksijklf
      endif 
      ksqijklz = ksqjkl
      do isymkl = 1, nsym
         isymaj  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
               do l = 1, nrhf(isyml)
                  idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                  do isym =1, nsym
                     isyma = muld2h(1,isym)
                     isymj  = muld2h(isymkl,isym)
                     do i = 1,nrhf(isymj)
                        call dcopy(nrhffr(isym),work(ksijklf),1,
     &                             work(ksqjkl),1)
                        ksqjkl = ksqjkl + nrhffr(isym)
                        ksijklf = ksijklf+ nrhffr(isym)
                        call dcopy(nrhf(isym),work(kspjkl),1,
     &                             work(ksqjkl),1)
                        ksqjkl = ksqjkl + nrhf(isym)
                        kspjkl = kspjkl + nrhf(isym)
                        call dcopy(nvir0(isym),work(ksajkl),1,
     &                             work(ksqjkl),1)
                        ksqjkl = ksqjkl + nvir0(isym)
                        ksajkl = ksajkl + nvir0(isym)
                      enddo
                  enddo
               enddo
             enddo
          enddo
      enddo


      call daxpy(ntklpj,one,work(ksqijklz),1,work(kxpjklz),1)
c               
      
c--------------------------------------------------------------------
c     Calculate X_{pj}^{kl}*D_{kl}^{ij}
c---------------------------------------------------------------------

      call dzero(work(kxdmmu),nr1bast)
      do isymkl = 1, nsym
         isymaj  = isymkl
         isymij  = isymkl
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            do k = 1, nrhf(isymk)
              do l = 1, nrhf(isyml)
                 idxkl=imatij(isymk,isyml)+nrhf(isymk)*(l-1)+k
                 do isyma =1, nsym
                    isymj  = muld2h(isymaj,isyma)
                    isymi  = muld2h(isymij,isymj)
                    koffd  = 1 + igamsq(isymij,isymkl)
     &                 + nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    ntota = max(1,norb1(isyma))
                    ntotj = max(1,nrhf(isymj))
                    ntoti = max(1,nrhf(isymi))
                    ntota1= max(1,nbas(isyma)) 
                    lxdmmu= kxdmmu + imatmuj(isyma,isymi)
                    call dgemm('N','T',norb1(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kxpjklz),
     &                    ntota,dklij(koffd),ntoti,one,
     &                    work(lxdmmu),ntota)

                    kxpjklz = kxpjklz + norb1(isyma)*nrhf(isymj)
                    ksqijklz= ksqijklz+ norb1(isyma)*nrhf(isymj)
                    ksijkl  = ksijkl+ nrhffr(isyma)*nrhf(isymj)
                 enddo
              enddo
           enddo
         enddo
      enddo



    
      call dzero(work(kxdpmu),mbas1t*nrhft)
      call dzero(work(kxdnumu),10*mbas1t**2)
      do isymaj = 1, nsym
         isymkl  = isymaj
         isymij  = isymaj
         noffset1 = 0
         do isyma =1, nsym
            isymmu = muld2h(1,isyma)
            isymj  = muld2h(isymaj,isyma)
            isymi  = muld2h(isymij,isymj)
            ntota = max(1,mbas1(isymi))
            ntotj = max(1,norb1(isyma))
            ntoti = max(1,mbas1(isymi))
            ntota1= max(1,nbas(isymi))
            koffc = kcmo + noffset +
     &              iglmrv(isyma,isymi)
            lxdmmu = kxdmmu + imatmuj(isyma,isymi)
            lxdpmu = kxdpmu + imatmuj(isyma,isymi)
            lxdnumu = kxdnumu + imatmunu(isymi,isymi)

            call dgemm('N','N',mbas1(isymi),
     &           nrhf(isyma),norb1(isyma),one,
     &           work(koffc),
     &           ntota1,work(lxdmmu),ntotj,zero,
     &           work(lxdpmu),ntota)
          
            noffset1 = noffset1 + nbas(isymi)*nrhffr(isymi)
            koffd = kcmo + noffset + noffset1 + 
     &              iglmrv1(isyma,isymi)

            call dgemm('N','T',mbas1(isymi),
     &           mbas1(isymi),nrhf(isyma),one,work(koffd),
     &           ntota1,work(lxdpmu),ntoti,zero,
     &           work(lxdnumu),ntota)

         enddo
      enddo



      call dzero(work(kzmunu),n2basx)
      call dzero(work(kdens),n2basx) 
      koff = 0
      do isymij = 1,nsym
         isymj = isymij
         isymi = isymij
         koff = imatmunu(isymi,isymj)
         do j = 1,mbas1(isymj)
            do i = 1,mbas1(isymi)
               koff = koff+1
               idxij = iaodis(isymi,isymj)+nbas(isymi)*
     &                 (j-1)+i
               work(kzmunu+idxij-1) = work(kxdnumu+koff-1)
             enddo
          enddo
      enddo


      call cc_r12fsym(imatmunu,work(kzmunu),work(kdens))

       
      call r12exdens(res,work(kcmo),work(kdens),work(kend1),lwrk1)

 
      call qexit('r12gradbd')
      end
c------------------------------------------------------------------------
C=====================================================================*
C  /* Deck r12exdens*/
      subroutine r12exdens(res,cmo,dens,work,lwork)
      implicit none
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>
#include <r12int.h>
#include <priunit.h>
#include <ccfro.h>
#include <inftap.h>
#include <ccsdinp.h>
ccc#include <inforb.h>
      integer kcmo,lwork,kfrsav,kfree,lfree,kfvao,kfcao,kdcao,kdvao,kdv
      integer pq,kcao,isym,ioff,ioff1,kexia,isymi,NNORB1
      integer kctao,koff,index,ISYMAI,isyma,nnbast
      integer idxkl,idxlk,idxij,idxji,idxklij,idxijkl,idxijlk
      integer idxjikl,isymj,isymij,isymkl,isymk,isyml,idxklji
      integer koffc,kexiaz,noffset,nvir0(8),kextia,kfctao
      integer icoun1,imatpq(8,8),icoun2,imatsq(8,8),ij,ji
      integer kajfr,isymaj,luvajkl,ncfai,kexfa,kexfaz
      logical onlyfc
c
      integer icoun3,imatfq(8,8),kxij,kxijij,ncvij,nvirf(8)
#include <maxorb.h>
#include <infinp.h>
      real*8 res(*),cmo(*),emcmy,zero,one,two,work(*),dens(*),densx
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)

      index(i,j) = max(i,j)*(max(i,j)-3)/2 +i+j

      call qenter('r12exdens')


      nrhft = 0
      do isymi = 1, nsym
         nvir0(isymi) = norb1(isymi) - nrhfs(isymi)
         nvirf(isymi) = nvir(isymi) - nrhffr(isymi)
         nrhft        = nrhft + nrhf(isymi)
      end do

      do isymaj = 1,nsym
         isymij = isymaj
         ncfai = 0
         do isyma = 1,nsym
            isymj = muld2h(isymaj,isyma)
            isymi = muld2h(isymij,isymj)
            ncfai = ncfai + (norb1(isyma) - nrhfs(isyma))*nrhffr(isymi)
         enddo
      enddo


      do isymi = 1,nsym
         icoun1 = 0
         icoun2 = 0
         icoun3 = 0
         do isymj = 1,nsym
            isymk = muld2h(isymj,isymi)
            imatpq(isymk,isymj) = icoun1
            imatsq(isymk,isymj) = icoun2
            imatfq(isymk,isymj) = icoun3
            icoun1 = icoun1 + nvir(isymk)*(nvir(isymj)+1)/2
            icoun2 = icoun2 + nrhf(isymk)*nvir0(isymj)
            icoun3 = icoun3 + nrhffr(isymk)*nvir0(isymj)
         enddo    
      enddo

      do isymaj = 1,nsym
         isymij = isymaj
         ncvij = 0
         do isyma = 1,nsym
            isymj = muld2h(isymaj,isyma)
            isymi = muld2h(isymij,isymj)
            ncvij = ncvij +  nrhffr(isyma)*nrhf(isymi)
         enddo
      enddo


      

      NBAST  = 0
      NNBAST = 0
      DO I = 1,NSYM
         NBAST  = NBAST  + NBAS(I)
         NNBAST = NNBAST + NBAS(I)*(NBAS(I)+1)/2
      ENDDO   
      
      noffset = 0
      do isym = 1,nsym
         noffset = noffset+nbas(isym)*nrhf(isym)
      enddo

      kfrsav = 1
      KFREE  = 1
      LFREE  = lwork
      CALL MEMGET('REAL',kfcao,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kfvao,     0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kdcao,N2BASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kdvao,     0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kdv  ,     0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kcmo, nlamds,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kexia,nnbast,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kexfa,nnbast,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kajfr,ncfai ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kxijij,ncvij,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',kxij  ,ncvij,WORK,KFREE,LFREE)

      call erisfk(dens,nbast,1)
      call dscal(n2basx,0.50d0,dens,1)
      

      DIRFCK = .TRUE.
      R12TRA = .TRUE.
      R12INT = .FALSE.
      V12INT = .TRUE.
      NORXR = NORXR .OR. R12HYB
      mbsmax = 4

      LUSUPM = -1
      call dzero(work(kfcao),n2basx)
      onlyfc = .TRUE.
      call fckmao(onlyfc,emcmy,work(kfcao),work(kfvao),
     &     dens,work(kdvao),work(kdv),cmo(1+noffset),work(kfree),lfree)
      
!     Note: dens is destroyed in fckmao because of dirfck, but that doesn't
!     matter as it is not used any more /Jan 2011

c
      WRITE(LUPRI,'(/A,I1,A)')
     &   ' Computation of exchange matrix done [',MBSMAX,']'
      

      R12INT = .TRUE.
      V12INT = .FALSE.
      call dzero(WORK(KDCAO),N2BASX)
      call dzero(WORK(kexia),NNBAST)
      CALL DCOPY(N2BASX,WORK(KFCAO),1,WORK(KDCAO),1)
c N2BASX = NBAST*NBAST!
      CALL DGETSP(NBAST,WORK(KDCAO),WORK(KFCAO))
      IF (NSYM .GT. 1) CALL PKSYM1(WORK(KFCAO),WORK(KFCAO),NBAS,NSYM,2)
      CALL DZERO(WORK(kexia),NNBAST)
      kexiaz = kexia
      koffc= 1
      CALL UTHUB(WORK(KFCAO),WORK(kexia),
     *           CMO(1+noffset),WORK(kfree),NSYM,nbas,nvir)


      koff = 0
      do isym = 1,nsym
         do p = 1+nrhfs(isym),norb1(isym)
            do q = 1+nrhffr(isym),nrhfs(isym)
               koff = imatsq(isym,isym)+nvir0(isym)*(q-nrhffr(isym)-1)
     &                        +p-nrhfs(isym)
               pq = imatpq(isym,isym)+index(p,q)
               res(koff) = 4.00d0*work(kexia+pq-1) 
            enddo
         enddo
      enddo

      if (froimp) then
      
          koff = 0
          do isym = 1,nsym
             do p = 1+nrhfs(isym),norb1(isym)
                do q = 1,nrhffr(isym)
                   koff = imatfq(isym,isym)+nvir0(isym)
     &                       *(q-1)+p - nrhfs(isym)
                   pq = imatpq(isym,isym)+index(p,q)
                   work(kajfr+koff-1) = 8.00d0*work(kexia+pq-1) 
                enddo
            enddo
          enddo

          luvajkl = -1
          call gpopen(luvajkl,'CCR12YAIFR','unknown',' ','unformatted',
     &                idummy,.false.)
          write(luvajkl) (work(kajfr+i-1), i = 1,ncfai)
          call gpclose(luvajkl,'KEEP')

      
      endif 


      CALL MEMREL('r12exdens',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      call qexit('r12exdens')
      end
c--------------------------------------------------------------------------
