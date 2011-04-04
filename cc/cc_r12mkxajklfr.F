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
      subroutine cc_r12mkxijkl(kvajkl,cmo,work,lwork,fel)
c---------------------------------------------------------------------
c     purpose: calculate X_{Ij}^{kl}, X = V,T necessary for the
c     calculation of MP2-R12 frozen core contributions to the
c     MP2-R12 gradient
c---------------------------------------------------------------------

      implicit none
#include <priunit.h>
#include <dummy.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
      logical ldum
#if defined (SYS_CRAY)
      real zero,one
#else
      double precision zero,one
#endif
      logical fel,test
      parameter (test = .true.,zero = 0.0d0, one = 1.0d0)

      integer  lwork,luvajkl,lwrk1,lwrk2,isymai
      integer isymkl,isymij,isymmu,isyml,isymk,isymak,isymaj,icount7,
     &     ntota,isymv,isymmuj,koffl,koffv,kvmujkl,
     &     isymj,isyma,idxkl,idum,xajkl,iglmrv(8,8)
      integer kend1,kend2,koffc,isymi,ntklaj,yajkl
      integer icou2,isymmua,ntotmu,nr1bas(8),ir1bas(8,8),icount2
      integer iglmrf(8,8),noffro(8),ntota1,isym
      integer icou3,iglmro(8,8),noffset,nb1fro(8),lu43,ntotj
      character*10 fnelena
      data fnelena/'FFR12xxxxx'/
#if defined (SYS_CRAY)
      real cmo(*),work(*),kvajkl(*)
#else
      double precision cmo(*),work(*),kvajkl(*)
#endif
c
      call qenter('cc_r12mkxijkl')
      fnelena = fnelena(1:5)//fnvajkl(6:10)
c
c
      DO ISYMAI = 1,NSYM
         Nb1FRO(ISYMAI) = 0
         DO ISYMI = 1,NSYM
            ISYMA = MULD2H(ISYMAI,ISYMI)
            Nb1FRO(ISYMAI) = Nb1FRO(ISYMAI) + NRHFFR(ISYMA)*NBAS(ISYMI)
         ENDDO
      ENDDO


      if (fnvajkl .eq. 'CCR12QAJKL' .and. fel .eqv. .true.) then
         lu43 = -43
         call dzero(cmo,nlamds)
         CALL GPOPEN(LU43,'EXFRO','UNKNOWN','SEQUENTIAL',
     &                   'FORMATTED',IDUMMY,LDUM)
         REWIND(LU43)
         READ(LU43,'(4E30.20)') (cmo(I), I = 1,NB1FRO(1))
         CALL GPCLOSE(LU43,'DELETE')
      endif


      do isymak = 1, nsym
        nr1bas(isymak) = 0
        do isymk = 1, nsym
           isyma = muld2h(isymak,isymk)
           nr1bas(isymak)  = nr1bas(isymak) +mbas1(isyma)*nrhf(isymk)
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
         do isymj = 1,nsym
            isyma = muld2h(isymaj,isymj)
            ntklaj = ntklaj + nmatij(isymkl)*nrhf(isymj)
     &               * nrhffr(isyma)
         enddo
      enddo

      do isymmua = 1,nsym
         icou2 = 0
         icou3 = 0
         do isyma = 1,nsym
            isymmu = muld2h(isymmua,isyma)
            iglmrf(isymmu,isyma) = icou2
            iglmro(isymmu,isyma) = icou3
            icou2 = icou2 + nbas(isymmu)*(nvir(isyma))
            icou3 = icou3 + nbas(isymmu)*nrhffr(isyma)
         enddo
      enddo

      noffset = 0
      do isym = 1,nsym
         noffset = noffset + nbas(isym)*nrhf(isym)
      enddo


      kend1 = 1 
      yajkl = kend1
      kend2 =  yajkl + ntklaj
      lwrk2  = lwork - kend2

      if (lwrk2 .lt.0) then
         call quit('Insufficient work space in cc_r12mkxijkl')
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
                  do isymmu =1, nsym
                     isyma = muld2h(1,isymmu)
                     isymj  = muld2h(isymmuj,isymmu)
                     if (fnvajkl .eq. 'CCR12QAJKL') then
                        koffc   = 1 + iglmro(isymmu,isyma)  
                     else
                        koffc   = 1 + noffset
     &                            + iglmrf(isymmu,isyma)
                     endif  
                     kvmujkl = 1+ivajkl(isymmuj,isymkl)+
     &                         nr1bas(isymmuj)*(idxkl-1)
     &                         +ir1bas(isymmu,isymj)
                     ntotmu = max(1,mbas1(isymmu))
                     ntota1 = max(1,nbas(isymmu))
                     ntota = max(1,nrhffr(isyma))
                     ntotj = max(1,nrhf(isymj))
c
                     call dgemm('T','N',nrhffr(isyma),
     &                 nrhf(isymj),mbas1(isymmu),one,cmo(koffc),
     &                 ntota1,kvajkl(kvmujkl),ntotmu,zero,
     &                 work(xajkl),ntota)
                     xajkl = xajkl + nrhffr(isyma)*nrhf(isymj)
                  end do
               end do
            end do
         end do
      end do


      luvajkl = -1
      call gpopen(luvajkl,fnelena,'unknown',' ','unformatted',
     &               idummy,.false.)
      rewind(luvajkl)
      write(luvajkl) (work(yajkl+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'KEEP')


      call qexit('cc_r12mkxijkl')
      end
C=====================================================================*
C  /* Deck cc_r12gradfr */
      subroutine cc_r12gradfr(ipdd,orb,v2am,work,lwork)
C     Calculate MP2-R12 contributions for frozen core gradients 
      implicit none
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>
#include <cbirea.h>
#include <r12int.h>
#include <priunit.h>
#include <ccfro.h>
      logical locdbg,ldum
      parameter (locdbg = .false.)
      character*10 CCR12YIJIJ
      integer kend1,kccr12,kend2,lwrk1,lwork,nrhftria,isymij
      integer lu43,idum,n2,lwrk2
      integer krsing,krtrip,kend3,lwrk3,nr1bas(8)
      integer isymkl,isyml,isymk,idxij,idxji,idxklij,idxklji
      integer isymai,isymv,isymi,isymj,idxkl,kctil,ntklaj
      integer kbajkl,kbijal,kvajkl,kvijal,isymaj,isyma,luvajkl
      integer lcvai,koffc,isymmu,ntota,ntotmu,kcvai
      integer ir1bas(8,8),isymmuj
      integer dklij,ntotk,lcvia,kvmujkl,ntoti
      integer koffd,bdajij,lccr12
      integer dklijz,kajij,kcvia,ntotj,imataj(8,8),icou3
      integer nr1orb(8),nr1xbas(8),nrgkl,n2bst1(8),ir1orb(8,8)
      integer ir2bas(8,8),ir2xbas(8,8),irgkl,iaodis1(8,8),nr2bas(8)
      integer kxij,ir1xbas(8,8),nijkl(8)
      integer ipdd,res,koffe,kctil2,idxijkl,lctil,ncvai

#if defined (SYS_CRAY)
      real work(*),one,two,zero
      real orb(*),v2am(*)
#else
      double precision work(*),one,zero,two
      double precision orb(*),v2am(*)
#endif
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)

      call qenter('cc_r12gradfr')

      nrhftria = nrhft*(nrhft+1)/2
      n2 = nrhftria * nrhftria

      kccr12  = 1
      krsing  = kccr12  + ngamsq(1)
      krtrip  = krsing  + n2
      kctil   = krtrip  + n2
      kctil2  = kctil   + ngamsq(1)
      kend1   = kctil2  + ngamsq(1)
      lwrk1   = lwork   - kend1

      call dzero(work(krsing),2*n2)
      call dzero(work(kccr12),ngamsq(1))
      call dzero(work(kctil),ngamsq(1))

      isymv = 1

      if (lwrk1 .lt. 0) then
         call quit('Insufficient work space in cc_r12gradfr')
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

      if (locdbg) then
         write(lupri,*) 'resorted R12 amplitudesfr (ctilde):'
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
     &                        ( nrhffr(isyma))
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
            icou3 = icou3+nrhf(isymi)*nrhffr(isyma)
         enddo
      enddo

      do isymaj = 1,nsym
         isymij = isymaj
         ncvai = 0
         do isyma = 1,nsym
            isymj = muld2h(isymaj,isyma)
            isymi = muld2h(isymij,isymj)
            ncvai = ncvai +  nrhffr(isyma)*nrhf(isymi)
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
      kxij    = res     + ncvai 
      kend2   = kxij    + ncvai 
      lwrk2   = lwork   - kend2

c     ----------------------------------------------------
c     Read V_{aj}^{kl},V_{kl}^{aj},B_{aj}^{kl},B_{kl}^{aj}
c               and transform it with Ctilde
c     -----------------------------------------------------

      if (lwrk2 .lt.0) then
         call quit('Insufficient work space in cc_r12mkxajkl')
      end if

      call dzero(work(kbajkl),ntklaj)
      call dzero(work(kbijal),ntklaj)
      luvajkl = -1
      call gpopen(luvajkl,'FFR12BAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kbajkl+i-1), i = 1,ntklaj)
      if (ipdd .eq. 5) then 
          call gpclose(luvajkl,'DELETE')
      else
          call gpclose(luvajkl,'KEEP')
      endif

      call gpopen(luvajkl,'FFR12BIJAL','unknown',' ','unformatted',
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
      call gpopen(luvajkl,'FFR12VAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kvajkl+i-1), i = 1,ntklaj)
      if (ipdd .eq. 5) then 
          call gpclose(luvajkl,'DELETE')
      else
          call gpclose(luvajkl,'KEEP')
      endif

      call gpopen(luvajkl,'FFR12VIJAL','unknown',' ','unformatted',
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
                    koffc  = kctil + igamsq(isymij,isymkl)
     &                 + nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    ntotmu = max(1,nrhffr(isyma))
                    ntotj = max(1,nrhf(isymj))
                    ntoti = max(1,nrhf(isymi))
                    lcvai  = kcvai + imataj(isyma,isymi)
                    call dgemm('N','T',nrhffr(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kvajkl),
     &                    ntotmu,work(koffc),ntoti,one,
     &                    work(lcvai),ntotmu)

                    kvajkl = kvajkl + nrhffr(isyma)*nrhf(isymj)

c     Calculate V_{kl}^{aj}*ctilde_{ij}^{kl}

                    lcvia  = kcvia + imataj(isyma,isymi)
                    koffe = kctil2+ igamsq(isymkl,isymij)
     &              +  nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    call dgemm('N','T',nrhffr(isyma),nrhf(isymi),
     &                         nrhf(isymj),one,work(kvijal),
     &                         ntotmu,work(koffe),ntoti,one,
     &                         work(lcvia),ntotmu)
                    kvijal = kvijal + nrhffr(isyma)*nrhf(isymj)
c
c     Calculate B_{aj}^{kl}*d_{kl}^{ij}
c
                    koffd = dklij + igamsq(isymij,isymkl)
     &              +  nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    bdajij  = kajij + imataj(isyma,isymi)
                    call dgemm('N','T',nrhffr(isyma),nrhf(isymi),
     &                         nrhf(isymj),one,work(kbijal),
     &                         ntotmu,work(koffd),ntoti,one,
     &                         work(bdajij),ntotmu)
                    kbijal  = kbijal + nrhffr(isyma)*nrhf(isymj)
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
         write(lupri,*) 'Matrixelemente MP2-R12 Gradient(frozen)'
         do isymaj = 1,nsym
            do isyma = 1,nsym
               isymi =isyma
               write(lupri,*) 'symmetry block',isyma
               call output(work(kajij+imataj(isyma,isymi)),
     &              1,nrhf(isyma),1,nrhffr(isymi),
     &              nrhf(isyma),nrhffr(isymi),1,lupri)

            enddo
         enddo
      end if

      if (ipdd .eq. 2) then
         luvajkl = -1
         call gpopen(luvajkl,'CCR12YIJIJ','unknown',' ','unformatted',
     &               idummy,.false.)
         write(luvajkl) (work(kajij+i-1), i = 1,ncvai)
         call gpclose(luvajkl,'KEEP')
      elseif (ipdd .eq. 3) then
         call dzero(work(res),ncvai)
         call r12gradapfr(ipdd,work(res),orb,work(kend2),v2am,lwrk2)
         call dscal(ncvai,2.0D0,work(res),1)
         call daxpy(ncvai,one,work(kajij),1,work(res),1)
         luvajkl = -1
         call gpopen(luvajkl,'CCR12ZIJIJ','unknown',' ','unformatted',
     &                idummy,.false.)
         write(luvajkl) (work(res+i-1), i = 1,ncvai)
         call gpclose(luvajkl,'KEEP')
      else if (ipdd .eq. 5) then
        call dzero(work(res),ncvai)
        call r12gradapfr(ipdd,work(res),orb,work(kend2),v2am,lwrk2)
        call dscal(ncvai,2.00D0,work(res),1)
        call daxpy(ncvai,one,work(kajij),1,work(res),1)
        luvajkl = -1
        call gpopen(luvajkl,'CCR12XIJIJ','unknown',' ','unformatted',
     &               idummy,.false.)
        write(luvajkl) (work(res+i-1), i = 1,ncvai)
        call gpclose(luvajkl,'KEEP')
      endif 
c
c
      call qexit('cc_r12gradfr')
      end
c------------------------------------------------------------------------
C=====================================================================*
C  /* Deck r12gradapfr */
      subroutine r12gradapfr(ipdd,res,orb,work,v2am,lwork)
      implicit none
#include <ccsdsym.h>
#include <ccorb.h>
#include <dummy.h>
#include <r12int.h>
#include <priunit.h>
#include <ccfro.h>
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
      integer blajklz,aj,kl,ijkl,nrhfst,nrhffrt
      integer imataj(8,8),imatka(8,8),nmatkl1(8),nmatak(8)
      integer icou4,icou5,icou6,icou7,imatmn(8,8),imatkl1(8,8)
      integer nijkl(8)
#if defined (SYS_CRAY)
      real res(*),work(*),orb(*),v2am(*),one,two,zero
#else
      double precision work(*),orb(*),v2am(*),one,zero,two
      double precision res(*)
#endif
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)

      call qenter('r12gradapfr')

      isymv = 1

      nrhftria = nrhft*(nrhft+1)/2
      n2 = nrhftria * nrhftria

      do isyman = 1,nsym
         nmatan(isyman) = 0
         nmatkl1(isyman) = 0
         nmatak(isyman) = 0
         icou3 = 0
         icou4 = 0
         icou5 = 0
         do isyma = 1,nsym
            isymn = muld2h(isyman,isyma)
            icou4 = icou4 + nrhfs(isyma)*nrhfs(isymn)
            icou5 = icou5 + nrhffr(isyma)*nrhf(isymn)
            nmatan(isyman) = nmatan(isyman) + nrhffr(isyma)*nrhfs(isymn)
            nmatak(isyman) = nmatak(isyman) + nrhffr(isyma)*nrhf(isymn)
            nmatkl1(isyman) = nmatkl1(isyman)+nrhfs(isyma)*nrhfs(isymn)
          enddo
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
            imataj(isymn,isyma)= icou4
            imatkl1(isyma,isymn)= icou5
            imatmn(isyma,isymn)= icou6
            imatka(isymn,isyma)= icou7
            icou4 = icou4+nrhfs(isymn)*nrhffr(isyma)
            icou5 = icou5+nmatkl1(isymn)* nmatkl1(isyma)
            icou6 = icou6+nrhfs(isyma)*nrhfs(isymn)
            icou7 = icou7+nmatak(isyma)*nmatij(isymn)
         enddo
      enddo
      do isyman = 1,nsym
         do isyma = 1,nsym
            isymn = muld2h(isyma,isyman)
         enddo
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
     &                        (nrhffr(isyma))
                  enddo
               enddo
            enddo
         enddo
      enddo

      icou3 = 0
      do isyman = 1,nsym
         icou3 = 0
         do isyma = 1,nsym
            isymn = muld2h(isyman,isyma)
            imatan(isyma,isymn)= icou3
            icou3 = icou3+nrhf(isymn)*(nrhffr(isyma))
         enddo
      enddo

      do isymaj = 1,nsym
         isymij = isymaj
         ncvai = 0
         do isyma = 1,nsym
            isymj = muld2h(isymaj,isyma)
            isymi = muld2h(isymij,isymj)
            ncvai = ncvai + (nrhffr(isyma))*nrhf(isymi)
         enddo
      enddo

      ntaj = 0
      do isymkl = 1,nsym
         isymij = isymkl
         do isymj = 1,nsym
            isyma = muld2h(isymij,isymj)
             ntaj = ntaj + nrhf(isymj)
     &            *(nrhffr(isyma))
         enddo
      enddo

      nrhffrt = 0
      nrhfst = 0
      do isym = 1, nsym
         nrhfst     = nrhfst+ nrhfs(isym)
         nrhffrt     = nrhffrt+ nrhffr(isym)
      end do


      kxajkl  = kend1 
      klijan  = kxajkl + ntklaj
      dklij   = klijan + ncvai 
      koctil  = dklij  + ngamsq(1)
      kdotil  = koctil + ngamsq(1)
      klij    = kdotil + ngamsq(1)
      krajkl  = klij   + ngamsq(1)
      kxxr12  = krajkl + nrhfst*nrhfst**3 
      kxsing  = kxxr12 + ngamsq(1) 
      kxtrip  = kxsing + n2
      klijmn  = kxtrip + n2
      blajkl  = klijmn + ngamsq(1) 
      kxtil   = blajkl + nrhfst*nrhfst**3 
      blajklz = kxtil  + ngamsq(1)
      kend2   = blajklz+ ntklaj 
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
c



      luvajkl = -1
      call dzero(work(kxajkl),ntklaj)
      call gpopen(luvajkl,'FFR12XAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kxajkl+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'KEEP')



      call dzero(work(krajkl),nrhfst*nrhfst**3)
      luvajkl = -1
      call gpopen(luvajkl,'AUXQ12','unknown',' ','formatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(krajkl+i-1), i=1,nrhfst*nrhfst**3)
      call gpclose(luvajkl,'KEEP')



      call dzero(work(blajkl),nrhfst**4)
      call cc_r12xsort4(work(krajkl),work(blajkl))

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
     &                          (l-1)*nrhfs(isymk) + k
                          ij = imatmn(isymi,isymj)+
     &                          (j-1)*nrhfs(isymi) + i
                          ijkl = imatkl1(isymij,isymkl)+
     &                          ij +(kl-1)*nmatkl1(isymij)
                          work(blajklz+idx-1) = work(blajkl+ijkl-1)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
         enddo
      enddo


      call daxpy(ntklaj,one,work(blajklz),1,work(kxajkl),1)

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
                    ntota = max(1,nrhffr(isyma))
                    ntotn = max(1,nrhf(isymn))
                    ntoti = max(1,nrhf(isymi))
                    koffd = kdotil + igamsq(isymij,isymkl)
     &              +  nmatij(isymij)*(idxkl-1)+imatij(isymi,isymn)
                    klijanz = klijan + imatan(isyma,isymi)
                    call dgemm('N','T',nrhffr(isyma),nrhf(isymi),
     &                         nrhf(isymn),one,work(kxajkl),
     &                         ntota,work(koffd),ntoti,one,
     &                         work(klijanz),ntota)

                    kxajkl = kxajkl + nrhffr(isyma)*nrhf(isymn)
                    blajklz = blajklz + nrhffr(isyma)*nrhf(isymn)
                 enddo
              enddo
           enddo
         enddo
      enddo  


      call r12gradaxdfr(ipdd,work(kccr12),work(kctil),work(klij),
     &                work(kend2),v2am,res,lwrk2)
      call daxpy(ncvai,one,work(klijan),1,res(1),1)


      call qexit('r12gradapfr')
      end
c--------------------------------------------------------------------------
C=====================================================================*
C  /* Deck r12gradaxd */
      subroutine r12gradaxdfr(ipdd,kccr12,kctil,dklij,work,v2am,res,
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
      integer ipdd,resb,kxtil,kresz,icoun4,nmataj(8),ncvij
      integer nt1fm(8),nij,LUCJDK,kintfr,icoun1
      integer iv1fro(8,8),it1vm(8,8),icoun3,iv2fro(8,8)
      integer ir2fro(8,8)
#if defined (SYS_CRAY)
      real work(*),dklij(*),v2am(*),kctil(*),kccr12(*),one,two,zero
      real res(*) 
#else
      double precision work(*),dklij(*),v2am(*),kctil(*),one,zero,two
      double precision kccr12(*),res(*)
#endif
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)

      index(i,j) = max(i,j)*(max(i,j)-3)/2 +i+j

      call qenter('r12gradaxdfr')

C     Calculate some offsets

      isymv = 1

      nrhftria = nrhft*(nrhft+1)/2
      n2 = nrhftria * nrhftria

      do isymai = 1,nsym
         icoun1 = 0
         do isymi = 1,nsym
            isyma = muld2h(isymai,isymi)
            ir2fro(isyma,isymi) = icoun1
            icoun1 = nmatij(isyma)*nf1fro(isymi)
         enddo
      enddo



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
            ntaj = ntaj + nrhffr(isyma)*nrhf(isymi)
         enddo
      enddo


      ntklaj = 0
      do isymkl = 1,nsym
         isymaj = isymkl
         do isymj = 1,nsym
            isyma = muld2h(isymaj,isymj)
            ntklaj = ntklaj + nmatij(isymkl)*nrhf(isymj)
     &               *(nrhffr(isyma))
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
      kintfr  = resb    + ntaj 
      kend2   = kintfr  + nt2fro(1) 
      lwrk2   = lwork   - kend2 
      

      call dzero(work(rklij),ncvij)

      lunit = -1
      if (ipdd .eq. 3) then 
         call gpopen(lunit,'CCR12YIJ','unknown',' ','unformatted',
     &               idummy,.false.)
      else if (ipdd .eq. 5) then
         call gpopen(lunit,'CCR12YIJB','unknown',' ','unformatted',
     &               idummy,.false.)
      endif
      read(lunit)(work(rklij+i-1),i=1,ncvij)
      call gpclose(lunit,'KEEP') 


      call dzero(work(kintfr),nt2fro(1))
      LUCJDK = -1
      CALL GPOPEN(LUCJDK,'INCJDK','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUCJDK)
      READ(LUCJDK) (WORK(KINTFR+I-1), I = 1,NT2FRO(1))
      CALL GPCLOSE(LUCJDK,'KEEP')


      call dzero(work(kgiakj),ntklaj)
      call dzero(work(kgaijk),ntklaj)
      call dzero(work(kgaikj),ntklaj)
      call dzero(work(kiajk),ntklaj)

      idx = 0
      do isymka = 1,nsym
         isymij = isymka
         do isyma = 1, nsym
            isymk = muld2h(isymka,isyma)
            do isymj = 1,nsym
               isymi = muld2h(isymij,isymj)
               do i = 1,nrhf(isymi)
                  do j = 1,nrhf(isymj)
                     isymia = muld2h(isymi,isyma)
                     isymkj = muld2h(isymk,isymj)
                     isymja = muld2h(isymj,isyma)
                     isymik = muld2h(isymi,isymk)
                     isymaj = muld2h(isyma,isymj)
                     do k =  1, nrhf(isymk)
                        do a = 1,nrhffr(isyma)
                          idx = idx + 1 
                          nai = it1fro(isymi,isyma)
     *                        + nvir(isymi)*(a-1) + i+nrhffr(isymi)
                          njk = it1am(isymk,isymj)
     *                        + nvir(isymk)*(j-1) + k+nrhffr(isymk)
                          naijk = it2fro(isymia,isymkj)
     *                        + nt1fro(isymia)*(njk-1)+nai
                          nak = it1fro(isymk,isyma)
     *                         + nvir(isymk)*(a-1) + k+nrhffr(isymk)
                          nji = it1am(isymi,isymj)
     *                          + nvir(isymi)*(j-1) + i+nrhffr(isymi)
                          nakji = it2fro(isymka,isymij)
     *                            + nt1fro(isymka)*(nji-1)+nak 
                          nik = it1am(isymi,isymk)
     *                          + nvir(isymi)*(k-1) + i + nrhffr(isymi)
                          naj = it1fro(isymj,isyma)
     *                          + nvir(isymj)*(a-1) + j  +nrhffr(isymj)
                          najki = it2fro(isymaj,isymik)
     *                          + nt1fro(isymja)*(nik-1)+naj
                          work(kgiakj+idx-1) = work(kintfr+naijk-1)
                          work(kgaijk+idx-1) = work(kintfr+nakji-1)
                          work(kgaikj+idx-1) = work(kintfr+najki-1)
                          work(kiajk+idx-1)  = 4.0d0*work(kgaijk+idx-1)
     &                   - work(kgaikj+idx-1) - work(kgiakj+idx-1)
                      enddo
                    enddo
                  enddo
               enddo   
            enddo 
         enddo
      enddo




      call dzero(res(1),ntaj)
      kres = 1
      do isymij = 1,nsym
         isymkl = isymij
         isymaj = isymij
         do isyma =1, nsym
            isymj  = muld2h(isymaj,isyma)
            isymi  = muld2h(1,isyma)
            ntotk = max(1,nrhffr(isyma)*nrhf(isymj)) 
            ntoti = max(1,nmatij(isymkl)) 
            call dgemm('N','N',nrhf(isymj)*nrhffr(isyma),1,
     &                 nmatij(isymkl),one,work(kiajk),ntotk,
     &                 work(rklij),ntoti,zero,res(kres),ntotk) 
            kres = kres + nrhffr(isyma)*nrhf(isymj)
            kiajk = kiajk + nrhf(isymj)*nrhffr(isyma)
     &                      * nmatij(isymkl) 
         enddo
      enddo 

      call dscal(ntaj,0.5d0,res(1),1) 



      if (ipdd .eq. 5) then
         call dscal(ntaj,1.0d0,res(1),1) 
         call dzero(work(resb),ntaj)
         call r12gradbfr(work(resb),dklij,work(kend2),ntklaj,ntaj,lwrk2)
         call dscal(ntaj,1.0d0,work(resb),1) 
         call daxpy(ntaj,one,work(resb),1,res(1),1)
      endif
      call qexit('r12gradaxdfr')
      end
c--------------------------------------------------------------------------
C=====================================================================*
C  /* Deck r12gradbfr */
      subroutine r12gradbfr(res,dklij,work,ntklaj,ncvai,lwork)
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
      integer kcmo,dens,nrhffrt,kuijals,kqijals,kssajkl

      integer kccr12,krsing,krtrip,kctil,lccr12,lctil
      integer lusifc,nrhftria,n2,lu43,ntotk
#if defined (SYS_CRAY)
      real dklij(*),res(*),work(*),one,two,zero
#else
      double precision dklij(*),res(*),work(*),one,zero,two
#endif
      parameter (zero = 0.0d0,one = 1.0d0, two = 2.0d0)

      call qenter('r12gradbfr')

      isymv = 1
      
      nrhftria = nrhft*(nrhft+1)/2
      n2 = nrhftria * nrhftria


      nrhffrt = 0
      do isym = 1, nsym
         nrhffrt = nrhffrt + nrhffr(isym)
      end do

      icou3 = 0
      do isymai = 1,nsym
         icou3 = 0
         do isyma = 1,nsym
            isymi = muld2h(isymai,isyma)
            imataj(isyma,isymi)= icou3
            icou3 = icou3+nrhf(isymi)*(nrhffr(isyma))
         enddo
      enddo

  
      kqajkl = 1
      kqijal = kqajkl + ntklaj
      kuijal = kqijal + ntklaj
      krajkl = kuijal + ntklaj
      ksajkl = krajkl + nrhffrt*nrhft**3 
      kcqai  = ksajkl + nrhffrt*nrhft**3
      kcqia  = kcqai  + ncvai
      ktijal = kcqia  + ncvai
      kcuai  = ktijal + nrhffrt*nrhft**3 
      kuajkl = kcuai  + ncvai
      kcuia  = kuajkl + ntklaj
      kstijal= kcuia  + ncvai
      ksrajkl= kstijal+ ntklaj
      kssijal= ksrajkl+ ntklaj 
      kssajkl= kssijal+ ntklaj
      kuijals= kssajkl+ ntklaj 
      kqijals= kuijals+ ntklaj
      dens   = kqijals+ ntklaj 
      kcmo   = dens   + nbast**2
      kend1  = kcmo   + nlamds
      lwrk1  = lwork  - kend1




      luvajkl = -1
      call dzero(work(kqajkl),ntklaj)
      call gpopen(luvajkl,'FFR12QAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kqajkl+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'DELETE')


      luvajkl = -1
      call dzero(work(kqijal),ntklaj)
      call gpopen(luvajkl,'FFR12QIJAL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kqijal+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'DELETE')
  
      call dzero(work(kqijals),ntklaj)
      call cc_r12xsort5(work(kqijal),work(kqijals),nrhffr)

      luvajkl = -1
      call dzero(work(kuijal),ntklaj)
      call gpopen(luvajkl,'FFR12UIJAL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kuijal+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'DELETE')

      call dzero(work(kuijals),ntklaj)
      call cc_r12xsort3(work(kuijal),work(kuijals),nrhffr)

  
      luvajkl = -1
      call dzero(work(kuajkl),ntklaj)
      call gpopen(luvajkl,'FFR12UAJKL','unknown',' ','unformatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl) (work(kuajkl+i-1), i = 1,ntklaj)
      call gpclose(luvajkl,'DELETE')

      call dzero(work(krajkl),nrhffrt*nrhft**3)
      luvajkl = -1
      call gpopen(luvajkl,'AUXQSF12','unknown',' ','formatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(krajkl+i-1), i=1,nrhffrt*nrhft**3)
      call gpclose(luvajkl,'KEEP')

      call dzero(work(ksrajkl),ntklaj)  
      call cc_r12xsort5(work(krajkl),work(ksrajkl),nrhffr)


      call dzero(work(ksajkl),nrhffrt*nrhft**3)
      luvajkl = -1
      call gpopen(luvajkl,'AUXQSFJ12','unknown',' ','formatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(ksajkl+i-1), i=1,nrhffrt*nrhft**3)
      call gpclose(luvajkl,'KEEP')

      call dzero(work(kssijal),ntklaj)  
      call cc_r12xsort5(work(ksajkl),work(kssijal),nrhffr) 

      call dzero(work(kssajkl),ntklaj)  
      call cc_r12xsort3(work(ksajkl),work(kssajkl),nrhffr) 


      call dzero(work(ktijal),nrhffrt*nrhft**3)
      luvajkl = -1
      call gpopen(luvajkl,'AUXQSFI12','unknown',' ','formatted',
     &                   idummy,.false.)
      rewind(luvajkl)
      read(luvajkl,'(4E30.20)') (work(ktijal+i-1), i=1,nrhffrt*nrhft**3)
      call gpclose(luvajkl,'KEEP')
 
      call dzero(work(kstijal),ntklaj)  
      call cc_r12xsort3(work(ktijal),work(kstijal),nrhffr)




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
                    koffd  = 1 + igamsq(isymij,isymkl)
     &                 + nmatij(isymij)*(idxkl-1)+imatij(isymi,isymj)
                    ntota = max(1,nrhffr(isyma))
                    ntotj = max(1,nrhf(isymj))
                    ntoti = max(1,nrhf(isymi))
                    lcqai  = kcqai + imataj(isyma,isymi)
                    call dgemm('N','T',nrhffr(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kqajkl),
     &                    ntota,dklij(koffd),ntoti,one,
     &                    work(lcqai),ntota)
                    kqajkl = kqajkl + nrhffr(isyma)*nrhf(isymj)
                    ksrajkl = ksrajkl + nrhffr(isyma)*nrhf(isymj)

                    ntota = max(1,nrhffr(isyma))
                    ntotj = max(1,nrhf(isymj))
                    ntoti = max(1,nrhf(isymi))
                    lcqia  = kcqia + imataj(isyma,isymi)
                    call dgemm('N','T',nrhffr(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kqijal),
     &                    ntota,dklij(koffd),ntoti,one,
     &                    work(lcqia),ntota)

                    kqijal = kqijal + nrhffr(isyma)*nrhf(isymj)

                    lcuai  = kcuai + imataj(isyma,isymi)
                    call dgemm('N','T',nrhffr(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kuijal),
     &                    ntota,dklij(koffd),ntoti,one,
     &                    work(lcuai),ntota)

                    kuijal = kuijal + nrhffr(isyma)*nrhf(isymj)
                    kssijal = kssijal + nrhffr(isyma)*nrhf(isymj)
                    kssajkl = kssajkl + nrhffr(isyma)*nrhf(isymj)

                    lcuia  = kcuia + imataj(isyma,isymi)
                    call dgemm('N','T',nrhffr(isyma),nrhf(isymi),
     &                    nrhf(isymj),one,work(kuajkl),
     &                    ntota,dklij(koffd),ntoti,one,
     &                    work(lcuia),ntota)

                    kuajkl  = kuajkl + nrhffr(isyma)*nrhf(isymj)
                    kstijal = kstijal + nrhffr(isyma)*nrhf(isymj)
                    ktijal  = ktijal + nrhffr(isyma)*nrhf(isymj)

                 enddo
              enddo
           enddo
         enddo
      enddo 




      call daxpy(ncvai,one,work(kcqai),1,res,1)
      call daxpy(ncvai,one,work(kcuai),1,res,1)
      call daxpy(ncvai,one,work(kcqia),1,res,1)
      call daxpy(ncvai,one,work(kcuia),1,res,1)

      if (locdbg) then
         write(lupri,*) 'Ergebnisse r12gradbfr'
         do isymaj = 1,nsym
            do isyma = 1,nsym
               isymi =isyma
               write(lupri,*) 'symmetry block',isyma
               call output(res(1+imataj(isyma,isymi)),
     &                 1,nrhffr(isyma),1,nrhf(isymi),
     &                 nrhffr(isyma),nrhf(isymi),1,lupri)
            enddo
         enddo
      endif


      call qexit('r12gradbfr')
      end
C=====================================================================*
