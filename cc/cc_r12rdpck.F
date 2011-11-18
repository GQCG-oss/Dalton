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
      subroutine ccr12pck(vpck,isymv,sing,trip,nr12orb,nrhforb,nkilj)
c----------------------------------------------------------------------
c   purpose: pack R12 four-index quantities stored as singlet and 
c            triplet matrices (without exploiting symmetry)
c            to a symmetry packed triangular matrix
c
c   H. Fliegl, C. Haettig spring 2003 
c
c   modified by C. Neiss summer 2005:
c   nr12orb  input: number of r12-orbitals (index k, l)
c   nrhforb  input: number of active occ. orbitals (index i, j)
c   nkilj    output-dimension
c   dimensions are calculated from scratch, since this routine can
c   be called from the MP2-R12 OR the CC-R12 code, which have 
c   different definitions for "nrhf", thus we cannot use the
c   dimensions from the common block in general.        
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      integer isymv,isymi,isymj,isymk,isyml,isymij
      integer kt,lt, it,jt,isymki,isymlj,idxki,idxlj,idxkj,idxli,klij
      integer idxkilj,idxljki,idxkjli,idxlikj,idum,ij,isymkl,kl
      integer isymli,isymkj, nrhftria,isym,icoun1,icoun2
      integer nr12orb(8),nrhforb(8),irhforb(8),ir12orb(8),
     &        nkj(8),ikj(8,8),nkilj(8),ikilj(8,8),nr12t

      real*8 sing(*),trip(*),vpck(*),ffkl,ff,half

      logical locdbg
      parameter (locdbg = .false.)
      parameter (half = 0.5D0)
      integer index 
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j  
      
      call qenter('ccr12pck')
C
      nr12t  = 0
      icoun1 = 0
      icoun2 = 0
      do isym = 1, nsym
        irhforb(isym) = icoun1
        ir12orb(isym) = icoun2
        icoun1 = icoun1 + nrhforb(isym)
        icoun2 = icoun2 + nr12orb(isym)
        nr12t  = nr12t  + nr12orb(isym)
      end do

      do isym = 1, nsym
        icoun1 = 0
        do isymj = 1, nsym
          isymi = muld2h(isymj,isym)
          ikj(isymi,isymj) = icoun1
          icoun1 = icoun1 + nr12orb(isymi)*nrhforb(isymj)
        end do
        nkj(isym) = icoun1
      end do

      do isym = 1, nsym
        icoun1 = 0
        do isymlj = 1, nsym
          isymki = muld2h(isymlj,isym)
          ikilj(isymki,isymlj) = icoun1
          ikilj(isymlj,isymki) = icoun1
          if (isymlj.gt.isymki) then
            icoun1 = icoun1 + nkj(isymki)*nkj(isymlj)
          else if (isymlj.eq.isymki) then
            icoun1 = icoun1 + nkj(isymki)*(nkj(isymlj)+1)/2
          end if
        end do
        nkilj(isym) = icoun1
      end do
C
      nrhftria = nr12t*(nr12t+1)/2
      call dzero(vpck,nkilj(isymv))
 
      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            
            do k = 1, nr12orb(isymk)
               kt = ir12orb(isymk) + k
               do l = 1, nr12orb(isyml)
                  lt = ir12orb(isyml) + l
                  
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

                        do j = 1, nrhforb(isymj)
                           jt = irhforb(isymj)+j
                           do i = 1, nrhforb(isymi)
                              it = irhforb(isymi)+i
                              
                              if (jt.eq.it) then
                                 ff = ffkl * sqrt(2d0)
                              else
                                 ff = ffkl
                              end if

                              if (jt.le.it) then
                                 ij = index(it,jt)
                                 idxki = ikj(isymk,isymi)+
     &                                nr12orb(isymk)*(i-1)+k
                                 idxlj = ikj(isyml,isymj)+
     &                                nr12orb(isyml)*(j-1)+l
                                 idxkj = ikj(isymk,isymj)+
     &                                nr12orb(isymk)*(j-1)+k
                                 idxli =ikj(isyml,isymi)+
     &                                nr12orb(isyml)*(i-1)+l
                                 
                                 klij = nrhftria*(ij-1)+kl
                                 
                                 if      (isymki.eq.isymlj) then
                                    idxkilj = ikilj(isymki,isymlj)+
     &                                   index(idxlj,idxki)
                                    vpck(idxkilj)=ff*half*(sing(klij)
     &                                   +       trip(klij))
                                 else if (isymki.lt.isymlj) then
                                    idxkilj = ikilj(isymki,isymlj)+
     &                                   nkj(isymki)*(idxlj-1)+idxki
                                    vpck(idxkilj)=ff*half*(sing(klij)
     &                                   +       trip(klij))
                                 else if (isymki.gt.isymlj) then
                                    idxljki =ikilj(isymlj,isymki)+
     &                                   nkj(isymlj)*(idxki-1)+idxlj
                                    vpck(idxljki)=ff*half*(sing(klij)
     &                                   +       trip(klij))
                                 end if

                                 if (isymli.eq.isymkj) then
                                    idxkjli =ikilj(isymkj,isymli)+
     &                                   index(idxli,idxkj)
                                    vpck(idxkjli)=ff*half*(sing(klij)
     &                                   -       trip(klij))
                                 else if (isymli.lt.isymkj) then
                                    idxlikj =ikilj(isymli,isymkj)+
     &                                   nkj(isymli)*(idxkj-1)+idxli
                                    vpck(idxlikj)=ff*half*(sing(klij)
     &                                   -       trip(klij))
                                 else if (isymli.gt.isymkj) then
                                    idxkjli =ikilj(isymkj,isymli)+
     &                                   nkj(isymkj)*(idxli-1)+idxkj
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
      
      if (locdbg) then
        write(lupri,*) 'Result in CCR12PCK:'
C        call cc_prpr12(vpck,isymv,1,.false.)
        DO ISYMLJ = 1,NSYM
            ISYMKI = MULD2H(ISYMLJ,ISYMV)
            IF (ISYMKI.EQ.ISYMLJ) THEN
               WRITE(LUPRI,*) 'Symmetry block number(ki,lj): ',
     &              ISYMKI,ISYMLJ
               IF (NMATKI(ISYMKI).EQ.0) THEN
                  WRITE(LUPRI,*) 'This symmetry is empty'
               ELSE
                  CALL OUTPAK(vpck(ikilj(isymki,isymlj)+1),
     &                        nkj(ISYMKI),1,LUPRI)
               END IF
               WRITE(LUPRI,*) ' '
            ELSE IF (ISYMLJ.GT.ISYMKI) THEN
               WRITE(LUPRI,*) 'Symmetry block number(ki,lj): ',
     &              ISYMKI,ISYMLJ
               IF (nkj(ISYMKI).EQ.0 .OR. nkj(ISYMLJ).EQ.0) THEN
                  WRITE(LUPRI,*) 'This symmetry is empty'
               ELSE
                  CALL OUTPUT(vpck(ikilj(isymki,isymlj)+1),
     &                        1,nkj(ISYMKI),1,nkj(ISYMLJ),
     *                        nkj(ISYMKI),nkj(ISYMLJ),1,LUPRI)
               END IF
               WRITE(LUPRI,*) ' '
            ENDIF
         END DO
      end if

      call qexit('ccr12pck')
      return

      end 
*=====================================================================*
      subroutine ccr12unpck(vpck,isymv,sing,trip,nr12orb,nrhforb)
c----------------------------------------------------------------------
c   purpose: unpack R12 four-index quantities stored as a symmetry 
c            packed lower triangular matrix 
c            to singlet and triplet matrices
c
c   H. Fliegl, C. Haettig spring 2003
c
c   modified by C. Neiss summer 2005
c----------------------------------------------------------------------
      implicit none
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      integer isymv,isymi,isymj,isymk,isyml,isymij
      integer kt,lt, it,jt,isymki,isymlj,idxki,idxlj,idxkj,idxli,klij
      integer idxkilj,idxljki,idxkjli,idxlikj,idum,ij,isymkl,kl
      integer isymli,isymkj, nrhftria
      integer icoun1,icoun2,isym
      integer nrhforb(8),irhforb(8),nkj(8),ikj(8,8),ikilj(8,8)
      integer nr12orb(8),ir12orb(8),nr12t
      logical locdbg 
      parameter (locdbg = .false.)

      real*8 sing(*),trip(*),vpck(*),ffkl,ff,vklij,vlkij

      integer index 
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j  

      call qenter('ccr12unpck')
C
      nr12t  = 0
      icoun1 = 0
      icoun2 = 0
      do isym = 1, nsym
        irhforb(isym) = icoun1
        ir12orb(isym) = icoun2
        icoun1 = icoun1 + nrhforb(isym)
        icoun2 = icoun2 + nr12orb(isym)
        nr12t  = nr12t  + nr12orb(isym)
      end do
      do isym = 1, nsym
        icoun1 = 0
        do isymj = 1, nsym
          isymi = muld2h(isymj,isym)
          ikj(isymi,isymj) = icoun1
          icoun1 = icoun1 + nr12orb(isymi)*nrhforb(isymj)
        end do
        nkj(isym) = icoun1
      end do
      do isym = 1, nsym
        icoun1 = 0
        do isymlj = 1, nsym
          isymki = muld2h(isymlj,isym)
          ikilj(isymki,isymlj) = icoun1
          ikilj(isymlj,isymki) = icoun1
          if (isymlj.gt.isymki) then
            icoun1 = icoun1 + nkj(isymki)*nkj(isymlj)
          else if (isymlj.eq.isymki) then
            icoun1 = icoun1 + nkj(isymki)*(nkj(isymlj)+1)/2
          end if
        end do
      end do
C
      nrhftria = nr12t*(nr12t+1)/2
      call dzero(sing,nrhftria*nrhftria)
      call dzero(trip,nrhftria*nrhftria)

      do isymkl = 1, nsym
        isymij = muld2h(isymkl,isymv)
        do isyml =1, nsym
          isymk = muld2h(isymkl,isyml)
          
          do k = 1, nr12orb(isymk)
             kt = ir12orb(isymk) + k
             do l = 1, nr12orb(isyml)
                lt = ir12orb(isyml) + l
                
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

                      do j = 1, nrhforb(isymj)
                         jt = irhforb(isymj)+j
                         do i = 1, nrhforb(isymi)
                            it = irhforb(isymi)+i
                            
                            if (jt.eq.it) then
                               ff = ffkl / sqrt(2.0d0)
                            else
                               ff = ffkl
                            end if

                            if (jt.le.it) then
                               ij = index(it,jt)
                               idxki = ikj(isymk,isymi)+
     &                              nr12orb(isymk)*(i-1)+k
                               idxlj = ikj(isyml,isymj)+
     &                              nr12orb(isyml)*(j-1)+l
                               idxkj = ikj(isymk,isymj)+
     &                              nr12orb(isymk)*(j-1)+k
                               idxli =ikj(isyml,isymi)+
     &                              nr12orb(isyml)*(i-1)+l
                               
                               if (isymki.eq.isymlj) then
                                  idxkilj = ikilj(isymki,isymlj)+
     &                                 index(idxlj,idxki)
                                  vklij = vpck(idxkilj)
                               else if (isymki.lt.isymlj) then
                                  idxkilj = ikilj(isymki,isymlj)+
     &                                 nkj(isymki)*(idxlj-1)+idxki
                                  vklij = vpck(idxkilj)
                               else if (isymki.gt.isymlj) then
                                  idxljki = ikilj(isymlj,isymki)+
     &                                 nkj(isymlj)*(idxki-1)+idxlj
                                  vklij = vpck(idxljki)
                               end if

                               if (isymli.eq.isymkj) then
                                  idxkjli = ikilj(isymkj,isymli)+
     &                                 index(idxli,idxkj)
                                  vlkij = vpck(idxkjli)
                               else if (isymli.lt.isymkj) then
                                  idxlikj = ikilj(isymli,isymkj)+
     &                                 nkj(isymli)*(idxkj-1)+idxli
                                  vlkij = vpck(idxlikj)
                               else if (isymli.gt.isymkj) then
                                  idxkjli = ikilj(isymkj,isymli)+
     &                                 nkj(isymkj)*(idxli-1)+idxkj
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

      if (locdbg) then
        write(lupri,*) 'Result in CCR12UNPCK:'
        write(lupri,*) 'Singlet:'
        CALL OUTPUT(sing,1,nrhftria,1,nrhftria,nrhftria,nrhftria,
     &              1,LUPRI)
        write(lupri,*) 'Triplet:'
        CALL OUTPUT(trip,1,nrhftria,1,nrhftria,nrhftria,nrhftria,
     &              1,LUPRI)
      end if

      call qexit('ccr12unpck')
      end 
*=====================================================================*
      subroutine cclr_diasclr12(tampr12,fac,isymt)
c---------------------------------------------------------------------
c     purpose: scale the diagonal elements of R12 amplitudes with fac
c
c       tampr12 : R12 amplitudes stored as symmetry packed triangular
c                 matrix using ntr12am,itr12am and nmatki,imatki
c       fac     : scale factor
c       isymt   : symmetry of R12 amplitudes
c
c     H. Fliegl, C. Haettig spring 2003 
c---------------------------------------------------------------------
      implicit none
#include <ccorb.h>
#include <ccsdsym.h>

      real*8 tampr12(*),fac

      integer isymt,isymlj,idxlj,idxljlj
 
      if (isymt.eq.1) then
         do isymlj = 1, nsym
            do idxlj = 1,nmatki(isymlj)
               idxljlj = itr12am(isymlj,isymlj) + idxlj*(idxlj+1)/2
               tampr12(idxljlj) = tampr12(idxljlj)*fac
            end do
         end do
      endif
 
      return
      end
*======================================================================*
      subroutine cc_r12generaltf(xint,xtf,idelta,isymd,lam1,isyml1,
     &                           lam2,isyml2,lidxtf1,ioffl1,ioffl2,
     &                           nmat1,nr1,nb1,ixint,ixtf,nxtf,
     &                           ixtf1,igamxtf,work,lwork)
c----------------------------------------------------------------------
c     purpose: two index transformation of R12 and Coulomb integrals
c              with different transformation matrices
c
c              r^ab_kl --> r^MtNt_kl
c              g^ab_MN --> g^itjt_MN
c
c              xint: input integral
c              ixint: offset for xint
c              xtf : transformed integral
c              ixtf, nxtf, igamxtf: offsets and dimension for xtf
c              ixtf1: offset for one index transformed result
c              lam1: transformation matrix for first index
c              lam2: transformation matrix for second index
c              ioffl1: offset for lam1
c              ioffl2: offset for lam2
c              nr1, nmat1, nb1: dimensions for first transformation
c
c     lidxtf1 = .T.: transform only first index
c
c     result is ADDED to xtf!       
c
c     H. Fliegl, C. Haettig, winter 2004
c     modified by C. Neiss, 2005/2006
c----------------------------------------------------------------------
      implicit none
#include <priunit.h>
#include <ccorb.h>
#include <dummy.h>

      logical lidxtf1,locdbg
      parameter (locdbg = .false.)
      integer idelta,isymd,isyml1,isyml2,ioffl1(8,8),ioffl2(8,8),
     &        nmat1(8),nr1(8),ixint(8,8),ixtf(8,8),nxtf(8),
     &        ixtf1(8,8),igamxtf(8,8),lwork,isyma,isymkl,isymm,nkl,nm,
     &        koff1,koff2,ntotkl,ntota,ntotm,krtf,isymn,isymmn,idxmn,
     &        idxdn,nb1(8),na,m,n 

      real*8 xint(*),xtf(*),lam1(*),lam2(*),work(*),ddot,one,zero

      parameter(one = 1.0d0, zero = 0.0d0)
      call qenter('generaltf')

c     transform first index
      do isyma = 1, nsym
        isymkl = muld2h(isyma,isymd)
        isymm  = muld2h(isyml1,isyma)

        nkl = nmat1(isymkl)
        nm  = nr1(isymm)
        na  = nb1(isyma)
        if (lwork.lt.nkl*nm) then
          call quit('Insufficient work space in generaltf')
        end if

        koff1  = 1 + ixint(isyma,isymkl)
        koff2  = 1 + ioffl1(isyma,isymm)

        ntotkl = max(1,nkl)
        ntota  = max(1,na)
        ntotm  = max(1,nm)

        call dgemm('T','N',nm,nkl,na,
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
          call output(xint(koff1),1,nb1(isyma),1,nmat1(isymkl),
     &                nb1(isyma),nmat1(isymkl),1,lupri)
          write(lupri,*)'work:'
          call output(work,1,nr1(isymm),1,nmat1(isymkl),
     &                nr1(isymm),nmat1(isymkl),1,lupri)
        end if
c
        if (lidxtf1) then
          krtf = 1 + ixtf1(isymkl,isymm)
c         call dcopy(nkl*nm,work,1,xtf(krtf),1)
          call daxpy(nkl*nm,one,work,1,xtf(krtf),1)
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
     &                nb1(isymd)*(n-1)+idelta

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
       subroutine ccr12pck2(vpck,isymv,lproj,vunpck,trans,iopt)
C----------------------------------------------------------------------
C  purpose: transform array of four occupied indices, e.g. V_(kl)^(ij)
C           where indices (kl) and (ij) are packed, into V'_(kl)^(ij) =
C           V_(kl)^(ij) (+ V_(lk)^(ji) if lproj=.TRUE.)
C           V' is stored in triangular format with the indices (ki) and 
C           (lj) packed.
C           Dimensions/Offsets: V:  ntr12sq(isymv),itr12sq(isym1,isym2)
C                               V': ntr12am(isymv),itr12am(isym1,isym2) 
C           During transformation it is (certainly) not assumed that
C           V_(kl)^(ij) = V_(lk)^(ji).
C 
C           isymv    symmetry of matrix elements
C           lproj    flag whether to apply projection operator or not
C           trans    flag wheter lower indices (kl) are leading 
C                    (-> 'N'), or upper indices (ij) are leading
C                    (-> 'T') in vunpck
C           iopt     = 0: indices k, l, i, j are ALL occupied indices
C                    = 1: indices i, j are occ., k, l are R12-indices 
C                    = 2: indices k, l, i, j are ALL R12-indices
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
      integer iopt
      integer nocc(8),nr12(8),nij(8),nkl(8),nki(8),nkilj(8),nklij(8),
     &        iij(8,8),ikl(8,8),iki(8,8),ikilj(8,8),iklij(8,8),
     &        iijkl(8,8)
      character*1 trans
      logical lproj, locdbg
      PARAMETER (locdbg = .FALSE.)      

      real*8 vpck(*), vunpck(*)

      integer index 
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j  
      
      call qenter('ccr12pck2')

      if (locdbg) write(lupri,*) 'Entered CCR12PCK2'

      if ((trans.ne.'T').and.(trans.ne.'N')) 
     &   call quit('Illegal value for "TRANS" in CCR12PCK2')
c
      if (iopt.eq.0) then
        call icopy(8,nrhfa,1,nocc,1)
        call icopy(8,nrhfa,1,nr12,1)
        call icopy(8,nmatij,1,nij,1)
        call icopy(8,nmatij,1,nkl,1)
        call icopy(8,nmatij,1,nki,1)
        call icopy(64,imatij,1,iij,1)
        call icopy(64,imatij,1,ikl,1)
        call icopy(64,imatij,1,iki,1)
        call icopy(8,ngamma,1,nkilj,1)
        call icopy(64,igamma,1,ikilj,1)
        call icopy(64,igamsq,1,iklij,1)
        call icopy(64,igamsq,1,iijkl,1)
      else if (iopt.eq.1) then
        call icopy(8,nrhfa,1,nocc,1)
        call icopy(8,nrhfb,1,nr12,1)
        call icopy(8,nmatij,1,nij,1)
        call icopy(8,nmatkl,1,nkl,1)
        call icopy(8,nmatki,1,nki,1)
        call icopy(64,imatij,1,iij,1)
        call icopy(64,imatkl,1,ikl,1)
        call icopy(64,imatki,1,iki,1)
        call icopy(8,ntr12am,1,nkilj,1)
        call icopy(64,itr12am,1,ikilj,1)
        call icopy(64,itr12sq,1,iklij,1)
        call icopy(64,itr12sqt,1,iijkl,1)
      else if (iopt.eq.2) then
        call icopy(8,nrhfb,1,nocc,1)
        call icopy(8,nrhfb,1,nr12,1)
        call icopy(8,nmatkl,1,nij,1)
        call icopy(8,nmatkl,1,nkl,1)
        call icopy(8,nmatkl,1,nki,1)
        call icopy(64,imatkl,1,iij,1)
        call icopy(64,imatkl,1,ikl,1)
        call icopy(64,imatkl,1,iki,1)
        call icopy(8,nr12r12p,1,nkilj,1)
        call icopy(64,ir12r12p,1,ikilj,1)
        call icopy(64,ir12r12sq,1,iklij,1)
        call icopy(64,ir12r12sq,1,iijkl,1)
      else 
        call quit('Unknown IOPT in CCR12PCK2')
      end if
C
      call dzero(vpck,nkilj(isymv))
C
      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            
            do l = 1, nr12(isyml)
               do k = 1, nr12(isymk)
                    idxkl = ikl(isymk,isyml) + nr12(isymk)*(l-1) + k
                    idxlk = ikl(isyml,isymk) + nr12(isyml)*(k-1) + l
                    do isymi =1, nsym
                       isymj = muld2h(isymij,isymi)
                       isymki = muld2h(isymk,isymi)
                       isymlj = muld2h(isyml,isymj)
                       isymkj = muld2h(isymk,isymj)
                       isymli = muld2h(isyml,isymi)

                       do j = 1, nocc(isymj)
                          do i = 1, nocc(isymi)
                             idxij = iij(isymi,isymj)+
     &                            nocc(isymi)*(j-1)+i
                             idxji = iij(isymj,isymi)+
     &                            nocc(isymj)*(i-1)+j
                             idxki = iki(isymk,isymi)+
     &                            nr12(isymk)*(i-1)+k
                             idxlj = iki(isyml,isymj)+
     &                            nr12(isyml)*(j-1)+l
                             if (trans.eq.'N') then
                               idxklij = iklij(isymkl,isymij)+
     &                              nkl(isymkl)*(idxij-1)+idxkl
                               idxlkji = iklij(isymkl,isymij)+
     &                              nkl(isymkl)*(idxji-1)+idxlk
                             else if (trans.eq.'T') then
                               idxklij = iijkl(isymij,isymkl)+
     &                              nij(isymij)*(idxkl-1)+idxij
                               idxlkji = iijkl(isymij,isymkl)+
     &                              nij(isymij)*(idxlk-1)+idxji
                             end if
                           
                             if (isymki.eq.isymlj) then
                               if (idxki .le. idxlj) then 
                                 idxkilj = ikilj(isymki,isymlj)+
     &                                index(idxki,idxlj)
                                 if (lproj) then
                                   vpck(idxkilj)=vunpck(idxklij)+
     &                                           vunpck(idxlkji)
                                 else 
                                   vpck(idxkilj)=vunpck(idxklij)
                                 end if
                               end if
                             else if (isymki.lt.isymlj) then
                               idxkilj = ikilj(isymki,isymlj)+
     &                              nki(isymki)*(idxlj-1)+idxki
                               if (lproj) then
                                 vpck(idxkilj)=vunpck(idxklij)+
     &                                         vunpck(idxlkji)
                               else
                                 vpck(idxkilj)=vunpck(idxklij)
                               end if
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
       subroutine ccr12unpck2(vpck,isymv,vunpck,trans,iopt)
C----------------------------------------------------------------------
C  purpose: do the inverse of ccr12pck2 (with lproj=.F.): 
C           transform a symmetry packed triangular matrix
C           with indices (ki) and (lj) packed into a
C           square matrix with (kl) and (ij) packed
C 
C           isymv   symmetry of matrix elements
C           trans   flag wheter lower indices (kl) are leading
C                    (-> 'N'), or upper indices (ij) are leading
C                    (-> 'T') in vunpck
C           iopt    = 0: indices k, l, i, j are ALL occ. indices
C                   = 1: indices i, j are occ., k, l are R12-indices
C                   = 2: indices k, l, i, j are ALL R12-indices
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
      integer iopt
      integer nocc(8),nr12(8),nij(8),nkl(8),nki(8),nkilj(8),nklij(8),
     &        iij(8,8),ikl(8,8),iki(8,8),ikilj(8,8),iklij(8,8),
     &        iijkl(8,8)
      character*1 trans

      logical locdbg
      PARAMETER (locdbg = .false.)      

      real*8 vpck(*), vunpck(*)

      integer index 
      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j  
      
      call qenter('ccr12unpck2')

      if (locdbg) write(lupri,*) 'Entered CCR12UNPCK2'

      if ((trans.ne.'T').and.(trans.ne.'N')) then
         write(lupri,*) 'TRANS = ',trans
         call quit('Illegal value for "TRANS" in CCR12UNPCK2')
      end if
c
      if (iopt.eq.0) then
        call icopy(8,nrhfa,1,nocc,1)
        call icopy(8,nrhfa,1,nr12,1)
        call icopy(8,nmatij,1,nij,1)
        call icopy(8,nmatij,1,nkl,1)
        call icopy(8,nmatij,1,nki,1)
        call icopy(64,imatij,1,iij,1)
        call icopy(64,imatij,1,ikl,1)
        call icopy(64,imatij,1,iki,1)
        call icopy(8,ngamsq,1,nklij,1)
        call icopy(64,igamma,1,ikilj,1)
        call icopy(64,igamsq,1,iklij,1)
        call icopy(64,igamsq,1,iijkl,1)
      else if (iopt.eq.1) then
        call icopy(8,nrhfa,1,nocc,1)
        call icopy(8,nrhfb,1,nr12,1)
        call icopy(8,nmatij,1,nij,1)
        call icopy(8,nmatkl,1,nkl,1)
        call icopy(8,nmatki,1,nki,1)
        call icopy(64,imatij,1,iij,1)
        call icopy(64,imatkl,1,ikl,1)
        call icopy(64,imatki,1,iki,1)
        call icopy(8,ntr12sq,1,nklij,1)
        call icopy(64,itr12am,1,ikilj,1)
        call icopy(64,itr12sq,1,iklij,1)
        call icopy(64,itr12sqt,1,iijkl,1)
      else if (iopt.eq.2) then
        call icopy(8,nrhfb,1,nocc,1)
        call icopy(8,nrhfb,1,nr12,1)
        call icopy(8,nmatkl,1,nij,1)
        call icopy(8,nmatkl,1,nkl,1)
        call icopy(8,nmatkl,1,nki,1)
        call icopy(64,imatkl,1,iij,1)
        call icopy(64,imatkl,1,ikl,1)
        call icopy(64,imatkl,1,iki,1)
        call icopy(8,nr12r12sq,1,nklij,1)
        call icopy(64,ir12r12p,1,ikilj,1)
        call icopy(64,ir12r12sq,1,iklij,1)
        call icopy(64,ir12r12sq,1,iijkl,1)
      else
        call quit('Unknown IOPT in CCR12UNPCK2')
      end if
C
      call dzero(vunpck,nklij(isymv))

      do isymkl = 1, nsym
         isymij = muld2h(isymkl,isymv)
         do isyml =1, nsym
            isymk = muld2h(isymkl,isyml)
            
            do l = 1, nr12(isyml)
               do k = 1, nr12(isymk)
                 idxkl = ikl(isymk,isyml) + nr12(isymk)*(l-1) + k
                 idxlk = ikl(isyml,isymk) + nr12(isyml)*(k-1) + l
                 do isymi =1, nsym
                    isymj = muld2h(isymij,isymi)
                    isymki = muld2h(isymk,isymi)
                    isymlj = muld2h(isyml,isymj)
                    isymkj = muld2h(isymk,isymj)
                    isymli = muld2h(isyml,isymi)

                    do j = 1, nocc(isymj)
                       do i = 1, nocc(isymi)
                          idxij = iij(isymi,isymj)+
     &                         nocc(isymi)*(j-1)+i
                          idxji = iij(isymj,isymi)+
     &                         nocc(isymj)*(i-1)+j
                          idxki = iki(isymk,isymi)+
     &                         nr12(isymk)*(i-1)+k
                          idxlj = iki(isyml,isymj)+
     &                         nr12(isyml)*(j-1)+l
                          if (trans.eq.'N') then
                            idxklij = iklij(isymkl,isymij)+
     &                           nkl(isymkl)*(idxij-1)+idxkl
                            idxlkji = iklij(isymkl,isymij)+
     &                           nkl(isymkl)*(idxji-1)+idxlk
                          else if (trans.eq.'T') then
                            idxklij = iijkl(isymij,isymkl)+
     &                           nij(isymij)*(idxkl-1)+idxij
                            idxlkji = iijkl(isymij,isymkl)+
     &                           nij(isymij)*(idxlk-1)+idxji
                          end if
                        
                          if (isymki.eq.isymlj) then
                            if (idxki .le. idxlj) then
                              idxkilj = ikilj(isymki,isymlj)+
     &                              index(idxki,idxlj)
                              vunpck(idxklij)=vpck(idxkilj)
                              vunpck(idxlkji)=vpck(idxkilj)
                            end if
                          else if (isymki.lt.isymlj) then
                            idxkilj = ikilj(isymki,isymlj)+
     &                            nki(isymki)*(idxlj-1)+idxki
                            vunpck(idxklij)=vpck(idxkilj)
                          else if (isymki.gt.isymlj) then
                            idxljki =ikilj(isymlj,isymki)+
     &                            nki(isymlj)*(idxki-1)+idxlj
                            vunpck(idxklij)=vpck(idxljki)
                          end if

C      if (locdbg) then
C      write(lupri,*) 'idxij, idxkl, idxklij, idxki, idxlj, idxkilj: ',
C    &                 idxij, idxkl, idxklij, idxki, idxlj, idxkilj
C      end if
C                               if      (isymli.eq.isymkj) then
C                                  idxkjli = itr12am(isymkj,isymli)+
C     &                                 index(idxli,idxkj)
C                                  vunpck(idxklij)=vpck(idxkjli)
C                               else if (isymli.lt.isymkj) then
C                                  idxlikj = itr12am(isymkj,isymli)+
C     &                                 nmatki(isymli)*(idxkj-1)+idxli
C                                  vunpck(idxijkl)=vpck(idxlikj)
C                               else if (isymli.gt.isymkj) then
C                                  idxkjli =itr12am(isymkj,isymli)+
C     &                                 nmatki(isymkj)*(idxli-1)+idxkj
C                                  vunpck(idxijkl)=vpck(idxkjli)
C                               end if

                       end do
                    end do
                 end do
               end do
            end do
         end do
      end do                
      
      if (locdbg) then
        if (iopt.eq.1) then
          write(lupri,*) 'Result in CCR12UNPCK2:'
          call cc_prsqr12(vunpck,isymv,trans,1,.false.)
        end if
        write(lupri,*) 'Leaving CCR12UNPCK2'
      end if
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
c                 matrix using ntr12am,itr12am and nmatki,imatki
c       fac     : scale factor
c       isymt   : symmetry of R12 amplitudes
c
c     C. Neiss spring 2005 
c---------------------------------------------------------------------
      implicit none
#include <ccorb.h>
#include <ccsdsym.h>

      real*8 tampr12(*),fac

      integer isymt,isymki,isymlj,idxki,idxlj,idxkilj
      integer index

      index(i,j) = max(i,j)*(max(i,j)-3)/2 + i + j
 
      if (isymt.eq.1) then
         do isymki = 1, nsym
           isymlj  = isymki
           do idxki = 1, nmatki(isymki)
             do idxlj = 1, nmatki(isymlj)
               if (idxki.ne.idxlj) then
                 idxkilj = itr12am(isymki,isymlj) + index(idxki,idxlj)
                 tampr12(idxkilj) = tampr12(idxkilj)*fac
               end if
             end do
           end do
         end do
      else
         call dscal(ntr12am(isymt),fac,tampr12,1) 
      endif
 
      return
      end
*======================================================================*


