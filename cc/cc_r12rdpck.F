*=====================================================================*
      subroutine ccr12pck(vpck,isymv,sing,trip)
c----------------------------------------------------------------------
c   purpose: pack R12 amplitudes stored as singlet and triplet matrices 
c            (without exploiting symmetry) to a symmetry packed 
c            triangular matrix using ngamma,igamma and nmatij,imatij
c
c   H. Fliegl, C. Haettig spring 2003         
c----------------------------------------------------------------------
      implicit none
#include <ccorb.h>
#include <ccsdsym.h>

      integer isymv,isymi,isymj,isymk,isyml,isymij
      integer kt,lt, it,jt,isymki,isymlj,idxki,idxlj,idxkj,idxli,klij
      integer idxkilj,idxljki,idxkjli,idxlikj,idum,ij,isymkl,kl
      integer isymli,isymkj, nrhftria

#if defined (SYS_CRAY)
      real sing(*),trip(*),vpck(*),ffkl,ff,half
#else
      double precision sing(*),trip(*),vpck(*),ffkl,ff,half
#endif
      parameter (half = 0.5D0)
      integer index 
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j  
      
      call qenter('ccr12pck')

      nrhftria = nrhft*(nrhft+1)/2

      call dzero(vpck,ngamma(isymv))
 
      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            
            do k = 1, nrhf(isymk)
               kt = irhf(isymk) + k
               do l = 1, nrhf(isyml)
                  lt = irhf(isyml) + l
                  
                  if (lt.eq.kt) then
                     ffkl = sqrt(2d0)
                  else
                     ffkl = 1d0
                  end if
 
                  if (lt.le.kt) then
                     kl = index(kt,lt)

                     do isymi =1, nsym
                        isymj = muld2h(isymij,isymi)
                        isymki = muld2h(isymk,isymi)
                        isymlj = muld2h(isyml,isymj)
                        isymkj = muld2h(isymk,isymj)
                        isymli = muld2h(isyml,isymi)

                        do j = 1, nrhf(isymj)
                           jt = irhf(isymj)+j
                           do i = 1, nrhf(isymi)
                              it = irhf(isymi)+i
                              
                              if (jt.eq.it) then
                                 ff = ffkl * sqrt(2d0)
                              else
                                 ff = ffkl
                              end if

                              if (jt.le.it) then
                                 ij = index(it,jt)
                                 idxki = imatij(isymk,isymi)+
     &                                nrhf(isymk)*(i-1)+k
                                 idxlj = imatij(isyml,isymj)+
     &                                nrhf(isyml)*(j-1)+l
                                 idxkj = imatij(isymk,isymj)+
     &                                nrhf(isymk)*(j-1)+k
                                 idxli =imatij(isyml,isymi)+
     &                                nrhf(isyml)*(i-1)+l
                                 
                                 klij = nrhftria*(ij-1)+kl
                                 
                                 if      (isymki.eq.isymlj) then
                                    idxkilj = igamma(isymki,isymlj)+
     &                                   index(idxlj,idxki)
                                    vpck(idxkilj)=ff*half*(sing(klij)
     &                                   +       trip(klij))
                                 else if (isymki.lt.isymlj) then
                                    idxkilj = igamma(isymki,isymlj)+
     &                                   nmatij(isymki)*(idxlj-1)+idxki
                                    vpck(idxkilj)=ff*half*(sing(klij)
     &                                   +       trip(klij))
                                 else if (isymki.gt.isymlj) then
                                    idxljki =igamma(isymlj,isymki)+
     &                                   nmatij(isymlj)*(idxki-1)+idxlj
                                    vpck(idxljki)=ff*half*(sing(klij)
     &                                   +       trip(klij))
                                 end if

                                 if (isymli.eq.isymkj) then
                                    idxkjli =igamma(isymkj,isymli)+
     &                                   index(idxli,idxkj)
                                    vpck(idxkjli)=ff*half*(sing(klij)
     &                                   -       trip(klij))
                                 else if (isymli.lt.isymkj) then
                                    idxlikj =igamma(isymli,isymkj)+
     &                                   nmatij(isymli)*(idxkj-1)+idxli
                                    vpck(idxlikj)=ff*half*(sing(klij)
     &                                   -       trip(klij))
                                 else if (isymli.gt.isymkj) then
                                    idxkjli =igamma(isymkj,isymli)+
     &                                   nmatij(isymkj)*(idxli-1)+idxkj
                                    vpck(idxkjli)=ff*half*(sing(klij)
     &                                   -       trip(klij))
                                 end if                        
                              end if
                              
                           end do
                        end do
                        
                     end do
                     
                  end if
                  
               end do
            end do
            
         end do
      end do                
      
      call qexit('ccr12pck')
      return

      end 
*=====================================================================*
      subroutine ccr12unpck(vpck,isymv,sing,trip)
c----------------------------------------------------------------------
c   purpose: unpack R12 amplitudes stored as a symmetry packed lower
c            triangular matrix (using ngamma,igamma and nmatij,imatij)
c            to singlet and triplet matrices (without exploiting sym.)
c
c   H. Fliegl, C. Haettig spring 2003
c----------------------------------------------------------------------
      implicit none
#include <ccorb.h>
#include <ccsdsym.h>

      integer isymv,isymi,isymj,isymk,isyml,isymij
      integer kt,lt, it,jt,isymki,isymlj,idxki,idxlj,idxkj,idxli,klij
      integer idxkilj,idxljki,idxkjli,idxlikj,idum,ij,isymkl,kl
      integer isymli,isymkj, nrhftria
#if defined (SYS_CRAY)
      real sing(*),trip(*),vpck(*),ffkl,ff,vklij,vlkij
#else
      double precision sing(*),trip(*),vpck(*),ffkl,ff,vklij,vlkij
#endif

      integer index 
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j  

      call qenter('ccr12unpck')

      nrhftria = nrhft*(nrhft+1)/2
      
      call dzero(sing,nrhftria*nrhftria)
      call dzero(trip,nrhftria*nrhftria)
   
      do isymkl = 1, nsym
        isymij = muld2h(isymkl,isymv)
        do isyml =1, nsym
          isymk = muld2h(isymkl,isyml)
          
          do k = 1, nrhf(isymk)
             kt = irhf(isymk) + k
             do l = 1, nrhf(isyml)
                lt = irhf(isyml) + l
                
                if (lt.eq.kt) then
                   ffkl = 1.0d0/sqrt(2.0d0) 
                else
                   ffkl = 1.0d0
                end if
 
                if (lt.le.kt) then
                   kl = index(kt,lt)

                   do isymi =1, nsym
                      isymj = muld2h(isymij,isymi)
                      isymki = muld2h(isymk,isymi)
                      isymlj = muld2h(isyml,isymj)
                      isymkj = muld2h(isymk,isymj)
                      isymli = muld2h(isyml,isymi)

                      do j = 1, nrhf(isymj)
                         jt = irhf(isymj)+j
                         do i = 1, nrhf(isymi)
                            it = irhf(isymi)+i
                            
                            if (jt.eq.it) then
                               ff = ffkl / sqrt(2.0d0)
                            else
                               ff = ffkl
                            end if

                            if (jt.le.it) then
                               ij = index(it,jt)
                               idxki = imatij(isymk,isymi)+
     &                              nrhf(isymk)*(i-1)+k
                               idxlj = imatij(isyml,isymj)+
     &                              nrhf(isyml)*(j-1)+l
                               idxkj = imatij(isymk,isymj)+
     &                              nrhf(isymk)*(j-1)+k
                               idxli =imatij(isyml,isymi)+
     &                              nrhf(isyml)*(i-1)+l
                               
                               if      (isymki.eq.isymlj) then
                                  idxkilj = igamma(isymki,isymlj)+
     &                                 index(idxlj,idxki)
                                  vklij = vpck(idxkilj)
                               else if (isymki.lt.isymlj) then
                                  idxkilj = igamma(isymki,isymlj)+
     &                                 nmatij(isymki)*(idxlj-1)+idxki
                                  vklij = vpck(idxkilj)
                               else if (isymki.gt.isymlj) then
                                  idxljki =igamma(isymlj,isymki)+
     &                                 nmatij(isymlj)*(idxki-1)+idxlj
                                  vklij = vpck(idxljki)
                               end if

                               if (isymli.eq.isymkj) then
                                  idxkjli =igamma(isymkj,isymli)+
     &                                 index(idxli,idxkj)
                                  vlkij = vpck(idxkjli)
                               else if (isymli.lt.isymkj) then
                                  idxlikj =igamma(isymli,isymkj)+
     &                                 nmatij(isymli)*(idxkj-1)+idxli
                                  vlkij = vpck(idxlikj)
                               else if (isymli.gt.isymkj) then
                                  idxkjli =igamma(isymkj,isymli)+
     &                                 nmatij(isymkj)*(idxli-1)+idxkj
                                  vlkij = vpck(idxkjli)
                               end if                        
                      
                               klij = nrhftria*(ij-1)+kl
                               
                               sing(klij)=ff*(vklij+vlkij)
                               trip(klij)=ff*(vklij-vlkij)

                            end if
                            
                         end do
                      end do
                      
                   end do
                   
                end if
                
             end do
          end do
          
        end do
      end do                

      call qexit('ccr12unpck')
      end 
*=====================================================================*
      subroutine cclr_diasclr12(tampr12,fac,isymt)
c---------------------------------------------------------------------
c     purpose: scale the diagonal elements of R12 amplitudes with fac
c
c       tampr12 : R12 amplitudes stored as symmetry packed triangular
c                 matrix using ngamma,igamma and nmatij,imatij
c       fac     : scale factor
c       isymt   : symmetry of R12 amplitudes
c
c     H. Fliegl, C. Haettig spring 2003 
c---------------------------------------------------------------------
      implicit none
#include <ccorb.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      real tampr12(*),fac
#else
      double precision tampr12(*),fac
#endif
      integer isymt,isymlj,idxlj,idxljlj
 
      if (isymt.eq.1) then
         do isymlj = 1, nsym
            do idxlj = 1,nmatij(isymlj)
               idxljlj = igamma(isymlj,isymlj) + idxlj*(idxlj+1)/2
               tampr12(idxljlj) = tampr12(idxljlj)*fac
            end do
         end do
      endif
 
      return
      end
*======================================================================*
      subroutine cc_r12generaltf(xint,xtf,idelta,isymd,lam1,isyml1,
     &                           lam2,isyml2,lidxtf1,ioffl1,ioffl2,
     &                           nmat1,nr1,ixint,nxint,ixtf,nxtf,ixtf1,
     &                           igamxtf,work,lwork)
c----------------------------------------------------------------------
c     purpose: two index transformation of R12 and Coulomb intergrals
c              with different transformation matrices
c
c              r^ab_kl --> r^MtNt_kl
c              g^ab_MN --> g^itjt_MN
c
c              xint: input integral
c              ixint, nxint: offset and dimension for xint
c              xtf : transformed integral
c              ixtf, nxtf, igamxtf: offsets and dimension for xtf
c              ixtf1: offset for one index transformed result
c              lam1: transformation matrix for first index
c              lam2: transformation matrix for second index
c              ioffl1: offset for lam1
c              ioffl2: offset for lam2
c              nmat1, nr1: dimensions for first transformation
c
c     H. Fliegl, C. Haettig, winter 2004
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
#include <r12int.h>
#include <dummy.h>

      logical lidxtf1,locdbg
      parameter (locdbg = .false.)
      integer idelta,isymd,isyml1,isyml2,ioffl1(8,8),ioffl2(8,8),
     &        nmat1(8),nr1(8),ixint(8,8),nxint(8),ixtf(8,8),nxtf(8),
     &        ixtf1(8,8),igamxtf(8,8),lwork,isyma,isymkl,isymm,nkl,nm,
     &        koff1,koff2,ntotkl,ntota,ntotm,krtf,isymn,isymmn,idxmn,
     &        idxdn 
#if defined (SYS_CRAY)
      real xint(*),xtf(*),lam1(*),lam2(*),work(*),ddot,one,zero
#else
      double precision xint(*),xtf(*),lam1(*),lam2(*),work(*),ddot,
     &                 one,zero
#endif
      parameter(one = 1.0d0, zero = 0.0d0)
      call qenter('generaltf')

c     transform first index
      do isyma = 1, nsym
        isymkl = muld2h(isyma,isymd)
        isymm  = muld2h(isyml1,isyma)

        nkl = nmat1(isymkl)
        nm  = nr1(isymm)
        if (lwork.lt.nkl*nm) then
          call quit('Insufficient work space in generaltf')
        end if

        koff1  = 1 + ixint(isyma,isymkl)
        koff2  = 1 + ioffl1(isyma,isymm)

        ntotkl = max(1,nkl)
        ntota  = max(1,nbas(isyma))
        ntotm  = max(1,nm)

          call dgemm('T','N',nm,nkl,nbas(isyma),
     &             one,lam1(koff2),ntota,
     &             xint(koff1),ntota,
     &             zero,work,ntotm)

        if (locdbg) then
          write(lupri,*)'after first transformation in tfs:'
          write(lupri,*)'isymd,idelta,isyma,isymkl:',
     &                  isymd,idelta,isyma,isymkl
          write(lupri,*)'norm^2 after first tf',
     &      ddot(nmat1(isymkl)*nr1(isymm),work,1,work,1)
          write(lupri,*)'R^d_alpha,kl:'
          call output(xint(koff1),1,nbas(isyma),1,nmat1(isymkl),
     &                nbas(isyma),nmat1(isymkl),1,lupri)
          write(lupri,*)'work:'
          call output(work,1,nr1(isymm),1,nmat1(isymkl),
     &                nr1(isymm),nmat1(isymkl),1,lupri)
        end if
c
        if (lidxtf1) then
          krtf = 1 + ixtf1(isymkl,isymm)
          call dcopy(nkl*nm,work,1,xtf(krtf),1)
        else
c         transform second index
          isymn  = muld2h(isyml2,isymd)
          isymmn = muld2h(isymm,isymn) !isymmn = isymkl

          do m = 1, nr1(isymm)
            do n = 1, nr1(isymn)
              idxmn = ixtf(isymm,isymn) +
     &                nr1(isymm)*(n-1) + m
              krtf  = igamxtf(isymmn,isymkl) + idxmn
              idxdn = ioffl2(isymd,isymn) +
     &                nbas(isymd)*(n-1)+idelta

              call daxpy(nmat1(isymkl),lam2(idxdn),
     &                   work(m),nr1(isymm),
     &                   xtf(krtf),nxtf(isymmn))
            end do
          end do
c
        end if
      end do
      
      call qexit('generaltf')
      end
*======================================================================*
       subroutine ccr12pck2(vpck,isymv,lproj,vunpck,ltrans)
C----------------------------------------------------------------------
C  purpose: transform array of four occupied indices, e.g. V_(kl)^(ij)
C           where indices (kl) and (ij) are packed, into V'_(kl)^(ij) =
C           V_(kl)^(ij) (+ V_(lk)^(ji) if lproj=.TRUE.)
C           V' is stored in triangular format with the indices (ki) and 
C           (lj) packed.
C           Dimensions/Offsets: V:  ngamsq(isymv),igamsq(isym1,isym2)
C                               V': ngamma(isymv),igamma(isym1,isym2) 
C           During transformation it is (certainly) not assumed that
C           V_(kl)^(ij) = V_(lk)^(ji).
C 
C           isymv    symmetry of matrix elements
C           lproj    flag whether to apply projection operator or not
C           ltrans   flag wheter lower indices (kl) are leading 
C                    (-> .FALSE.), or upper indices (ij) are leading
C                    (-> .TRUE.) in vunpck
C
C  C. Neiss autumn 2004
C----------------------------------------------------------------------

      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      integer isymv,isymi,isymj,isymk,isyml
      integer isymij,isymkl
      integer idxij,idxji,idxkl,idxlk,idxklij,idxlkji
      integer isymki,isymlj,isymli,isymkj
      integer idxki,idxlj,idxli,idxkj
      integer idxkilj,idxljki,idxkjli,idxlikj

      logical ltrans
      logical lproj, locdbg
      PARAMETER (locdbg = .FALSE.)      

#if defined (SYS_CRAY)
      real precision vpck(*), vunpck(*)
#else
      double precision vpck(*), vunpck(*)
#endif
      integer index 
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j  
      
      call qenter('ccr12pck2')

      if (locdbg) write(lupri,*) 'Entered CCR12PCK2'

      call dzero(vpck,ngamma(isymv))

      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            
            do l = 1, nrhf(isyml)
               do k = 1, nrhf(isymk)
                    idxkl = imatij(isymk,isyml) + nrhf(isymk)*(l-1) + k
                    idxlk = imatij(isyml,isymk) + nrhf(isyml)*(k-1) + l
                    do isymi =1, nsym
                       isymj = muld2h(isymij,isymi)
                       isymki = muld2h(isymk,isymi)
                       isymlj = muld2h(isyml,isymj)
                       isymkj = muld2h(isymk,isymj)
                       isymli = muld2h(isyml,isymi)

                       do j = 1, nrhf(isymj)
                          do i = 1, nrhf(isymi)
                             idxij = imatij(isymi,isymj)+
     &                            nrhf(isymi)*(j-1)+i
                             idxji = imatij(isymj,isymi)+
     &                            nrhf(isymj)*(i-1)+j
                             idxki = imatij(isymk,isymi)+
     &                            nrhf(isymk)*(i-1)+k
                             idxlj = imatij(isyml,isymj)+
     &                            nrhf(isyml)*(j-1)+l
                             if (.not. ltrans) then
                               idxklij = igamsq(isymkl,isymij)+
     &                              nmatij(isymkl)*(idxij-1)+idxkl
                               idxlkji = igamsq(isymkl,isymij)+
     &                              nmatij(isymkl)*(idxji-1)+idxlk
                             else 
                               idxklij = igamsq(isymij,isymkl)+
     &                              nmatij(isymij)*(idxkl-1)+idxij
                               idxlkji = igamsq(isymij,isymkl)+
     &                              nmatij(isymij)*(idxlk-1)+idxji
                             end if
                           
                             if (isymki.eq.isymlj) then
                               if (idxki .le. idxlj) then 
                                 idxkilj = igamma(isymki,isymlj)+
     &                                index(idxki,idxlj)
                                 if (lproj) then
                                   vpck(idxkilj)=vunpck(idxklij)+
     &                                           vunpck(idxlkji)
                                 else 
                                   vpck(idxkilj)=vunpck(idxklij)
                                 end if
                               end if
                             else if (isymki.lt.isymlj) then
                               idxkilj = igamma(isymki,isymlj)+
     &                              nmatij(isymki)*(idxlj-1)+idxki
                               if (lproj) then
                                 vpck(idxkilj)=vunpck(idxklij)+
     &                                         vunpck(idxlkji)
                               else
                                 vpck(idxkilj)=vunpck(idxklij)
                               end if
C                             else if (isymki.gt.isymlj) then
C                               idxljki =igamma(isymlj,isymki)+
C     &                              nmatij(isymlj)*(idxki-1)+idxlj
C                               if (lproj) then
C                                 vpck(idxljki)=vunpck(idxklij)+
C     &                                         vunpck(idxlkji)
C                               else
C                                 vpck(idxljki)=vunpck(idxklij)
C                               end if

       if (locdbg) then
       write(lupri,*) 'idxij, idxkl, idxklij, idxki, idxlj, idxkilj: ',
     &                 idxij, idxkl, idxklij, idxki, idxlj, idxkilj
       if (lproj) then
       write(lupri,*) 'idxji, idxlk, idxlkji: ',
     &                 idxji, idxlk, idxlkji 
       write(lupri,*) 'vunpck(idxklij),vunpck(idxlkji),vpck(idxkilj): ',
     &                 vunpck(idxklij), vunpck(idxlkji), vpck(idxkilj)
       end if
       end if

C                               if      (isymli.eq.isymkj) then
C                                  idxkjli = igamma(isymkj,isymli)+
C     &                                 index(idxli,idxkj)
C                                  if (lproj) then
C                                    vpck(idxkjli)=vunpck(idxijkl)+
C     &                                 vunpck(idxjilk)
C                                  else
C                                    vpck(idxkjli)=vunpck(idxijkl)
C                                  end if
C                               else if (isymli.lt.isymkj) then
C                                  idxlikj = igamma(isymkj,isymli)+
C     &                                 nmatij(isymli)*(idxkj-1)+idxli
C                                  if (lproj) then
C                                    vpck(idxlikj)=vunpck(idxijkl)+
C     &                                 vunpck(idxjilk)
C                                  else
C                                    vpck(idxlikj)=vunpck(idxijkl)
C                                  end if
C                               else if (isymli.gt.isymkj) then
C                                  idxkjli =igamma(isymkj,isymli)+
C     &                                 nmatij(isymkj)*(idxli-1)+idxkj
C                                  if (lproj) then
C                                    vpck(idxkjli)=vunpck(idxijkl)+
C     &                                 vunpck(idxjilk)
C                                  else
C                                    vpck(idxkjli)=vunpck(idxijkl)
C                                  end if
C                               end if

                             end if

                          end do
                       end do
                    end do
               end do
            end do
         end do
      end do                
      
      if (locdbg) write(lupri,*) 'Leaving CCR12PCK2'
      call qexit('ccr12pck2')
      return

      end 
*======================================================================

*======================================================================
       subroutine ccr12unpck2(vpck,isymv,vunpck,ltrans)
C----------------------------------------------------------------------
C  purpose: do the inverse of ccr12pck2 (with lproj=.F.): 
C           transform a symmetry packed triangular matrix
C           with indices (ki) and (lj) packed into a
C           square matrix with (kl) and (ij) packed
C 
C           isymv    symmetry of matrix elements
C           ltrans   flag wheter lower indices (kl) are leading
C                    (-> .FALSE.), or upper indices (ij) are leading
C                    (-> .TRUE.) in vunpck
C
C  C. Neiss   04. Feb. 2005
C----------------------------------------------------------------------

      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      integer isymv,isymi,isymj,isymk,isyml
      integer isymij,isymkl
      integer idxij,idxji,idxkl,idxlk,idxklij,idxlkji
      integer isymki,isymlj,isymli,isymkj
      integer idxki,idxlj,idxli,idxkj
      integer idxkilj,idxljki,idxkjli,idxlikj
      logical ltrans

      logical locdbg
      PARAMETER (locdbg = .FALSE.)      

#if defined (SYS_CRAY)
      real precision vpck(*), vunpck(*)
#else
      double precision vpck(*), vunpck(*)
#endif
      integer index 
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j  
      
      call qenter('ccr12unpck2')

      if (locdbg) write(lupri,*) 'Entered CCR12UNPCK2'

      call dzero(vunpck,ngamsq(isymv))

      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            
            do l = 1, nrhf(isyml)
               do k = 1, nrhf(isymk)
                 idxkl = imatij(isymk,isyml) + nrhf(isymk)*(l-1) + k
                 idxlk = imatij(isyml,isymk) + nrhf(isyml)*(k-1) + l
                 do isymi =1, nsym
                    isymj = muld2h(isymij,isymi)
                    isymki = muld2h(isymk,isymi)
                    isymlj = muld2h(isyml,isymj)
                    isymkj = muld2h(isymk,isymj)
                    isymli = muld2h(isyml,isymi)

                    do j = 1, nrhf(isymj)
                       do i = 1, nrhf(isymi)
                          idxij = imatij(isymi,isymj)+
     &                         nrhf(isymi)*(j-1)+i
                          idxji = imatij(isymj,isymi)+
     &                         nrhf(isymj)*(i-1)+j
                          idxki = imatij(isymk,isymi)+
     &                         nrhf(isymk)*(i-1)+k
                          idxlj = imatij(isyml,isymj)+
     &                         nrhf(isyml)*(j-1)+l
                          if (.not. ltrans) then
                            idxklij = igamsq(isymkl,isymij)+
     &                           nmatij(isymkl)*(idxij-1)+idxkl
                            idxlkji = igamsq(isymkl,isymij)+
     &                           nmatij(isymkl)*(idxji-1)+idxlk
                          else 
                            idxklij = igamsq(isymij,isymkl)+
     &                           nmatij(isymij)*(idxkl-1)+idxij
                            idxlkji = igamsq(isymij,isymkl)+
     &                           nmatij(isymij)*(idxlk-1)+idxji
                          end if
                        
                          if (isymki.eq.isymlj) then
                            if (idxki .le. idxlj) then
                              idxkilj = igamma(isymki,isymlj)+
     &                              index(idxki,idxlj)
                              vunpck(idxklij)=vpck(idxkilj)
                              vunpck(idxlkji)=vpck(idxkilj)
                            end if
                          else if (isymki.lt.isymlj) then
                            idxkilj = igamma(isymki,isymlj)+
     &                            nmatij(isymki)*(idxlj-1)+idxki
                            vunpck(idxklij)=vpck(idxkilj)
                          else if (isymki.gt.isymlj) then
                            idxljki =igamma(isymlj,isymki)+
     &                            nmatij(isymlj)*(idxki-1)+idxlj
                            vunpck(idxklij)=vpck(idxljki)
                          end if

       if (locdbg) then
       write(lupri,*) 'idxij, idxkl, idxklij, idxki, idxlj, idxkilj: ',
     &                 idxij, idxkl, idxklij, idxki, idxlj, idxkilj
       end if
C                               if      (isymli.eq.isymkj) then
C                                  idxkjli = igamma(isymkj,isymli)+
C     &                                 index(idxli,idxkj)
C                                  vunpck(idxklij)=vpck(idxkjli)
C                               else if (isymli.lt.isymkj) then
C                                  idxlikj = igamma(isymkj,isymli)+
C     &                                 nmatij(isymli)*(idxkj-1)+idxli
C                                  vunpck(idxijkl)=vpck(idxlikj)
C                               else if (isymli.gt.isymkj) then
C                                  idxkjli =igamma(isymkj,isymli)+
C     &                                 nmatij(isymkj)*(idxli-1)+idxkj
C                                  vunpck(idxijkl)=vpck(idxkjli)
C                               end if

                       end do
                    end do
                 end do
               end do
            end do
         end do
      end do                
      
      if (locdbg) write(lupri,*) 'Leaving CCR12UNPCK2'
      call qexit('ccr12unpck2')
      return

      end 
*======================================================================

*=====================================================================*
      subroutine cclr_nondiasclr12(tampr12,fac,isymt)
c---------------------------------------------------------------------
c     purpose: scale the non-diagonal elements of R12 amplitudes 
c              (triangular matrix) with fac
c
c       tampr12 : R12 amplitudes stored as symmetry packed triangular
c                 matrix using ngamma,igamma and nmatij,imatij
c       fac     : scale factor
c       isymt   : symmetry of R12 amplitudes
c
c     C. Neiss spring 2005 
c---------------------------------------------------------------------
      implicit none
#include <ccorb.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      real tampr12(*),fac
#else
      double precision tampr12(*),fac
#endif
      integer isymt,isymki,isymlj,idxki,idxlj,idxkilj
      integer index

      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
 
      if (isymt.eq.1) then
         do isymki = 1, nsym
           isymlj  = isymki
           do idxki = 1, nmatij(isymki)
             do idxlj = 1, nmatij(isymlj)
               if (idxki.ne.idxlj) then
                 idxkilj = igamma(isymki,isymlj) + index(idxki,idxlj)
                 tampr12(idxkilj) = tampr12(idxkilj)*fac
               end if
             end do
           end do
         end do
      else
         call dscal(ngamma(isymt),fac,tampr12,1) 
      endif
 
      return
      end
*======================================================================*

*=====================================================================*
      subroutine cclr_diasclsqr12(matrix,fac,isym)
c---------------------------------------------------------------------
c     purpose: scale diagonal elements of square matrix with four
c              "occupied" indices  with fac
c              pure test routine!
c
c     C. Neiss spring 2005
c---------------------------------------------------------------------
      implicit none
#include <ccorb.h>
#include <ccsdsym.h>

#if defined (SYS_CRAY)
      real matrix(*),fac
#else
      double precision matrix(*),fac
#endif
      integer isym,isymkl,idxkl,idxklkl
      integer index

      if (isym.eq.1) then
        do isymkl = 1, nsym
          do idxkl = 1, nmatij(isymkl)
            idxklkl = igamsq(isymkl,isymkl) + 
     &                nmatij(isymkl)*(idxkl-1) + idxkl
            matrix(idxklkl) = matrix(idxklkl)*fac
          end do
        end do
      endif

      return
      end
*======================================================================*

