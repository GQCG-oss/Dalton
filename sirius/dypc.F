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
      SUBROUTINE dypc(escf,onel,id2n,in2d,WRK,LWRK,
     $           norb,ncore,nact,its,itsym,nocb,
     $           LUPRI,LUN98,LUN31,zgel)
C-----------------------------------------------------------------------
C
C     This program calculates the second order correction to the energy
C     using the n-electron valence state perturbation theory (NEVPT)
C     in the partially contracted (PC) ad strongly contracted (SC)
C     variants using the Dyall's Hamiltonian for the definition of the
C     zero order energies. All equations are based on the spin free
C     formulation of NEVPT.
C
C     Relevant papers:
C     1) C. Angeli, R. Cimiraglia, J-P Malrieu, 
C        Chem. Phys. Lett.,  317, 472-480, (2000)
C        Definition of the inactive (core+virtual) orbital energies.
C     2) C. Angeli, R. Cimiraglia, S. Evangelisti, T. Leininger, J-P Malrieu,
C        J. Chem. Phys., 114, 10252, (2001)
C        Introduction of NEVPT
C     3) C. Angeli, R. Cimiraglia, J-P Malrieu,
C        Chem. Phys. Lett., 350, 297-305, (2001)
C        Introduction of the formalism based on the "Koopmans' operators"
C     4) C. Angeli, R. Cimiraglia, J-P Malrieu,
C        J. Chem. Phys., 117, 10525-10533, (2002)
C        Spin free version of NEVPT
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      COMMON /NEVPT/ e(1),e2pc(1),e2sc(1),psisc(1),psipc(1),
     *               e2pcc(8),psipcc(8),e2scc(8),psiscc(8),
     *               ntotdet,ntotconf,ntotcap ,zro
      dimension id2n(*),in2d(*),onel(*)
      dimension WRK(LWRK)
      dimension itsym(*),its(8,8),zgel(*)
      allocatable heffd(:),fmpb(:),igelo(:)
C
C     Generalized Koopmans matrices
C 
C     gkp1a = generalized Koopmans matrix for alpha electronic affinities
C     gkm1a = generalized Koopmans matrix for alpha ionization potentials
C     gkp2aa= generalized Koopmans matrix for alpha double EA 
C     ......
C     gk0pa = A matrix used for the 0' class (see report)
C     gk0pd = D matrix used for the 0' class (see report)
C     gk0pf = F matrix used for the 0' class (see report)
C
      allocatable gkp1a(:,:),gkm1a(:,:)
      allocatable gkp2aa(:,:,:,:)
      allocatable gkm2aa(:,:,:,:)
      allocatable gk0pa(:,:,:,:),gk0pd(:,:,:,:),gk0pf(:,:)
      allocatable ro3(:,:,:,:,:,:),amat(:,:,:,:,:,:),bmat(:,:,:,:),
     *            cmat(:,:,:,:),dmat(:,:)
      allocatable atmat(:,:,:,:,:,:),btmat(:,:,:,:),
     *            ctmat(:,:,:,:),dtmat(:,:)
      allocatable rho(:,:),w(:),work(:)
C
C     Density matrices SPIN LESS
C
C     da  = one-particle density matrix
C     daa = two-particle density matrix
C     daat= two-hole density matrix
C
      allocatable da(:,:)
      allocatable daa(:,:,:,:)
      allocatable daat(:,:,:,:)
c     indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

      CALL QENTER('DYPC')

      write(lupri,'(/A//)') ' STARTING THE PERTURBATION SUMMATION'
 
      nwall=0
      ntotdet=0
      ntotconf=0
      ntotcap=0
      zro=.false.
      t=second()
c      CALL OFILES
      rewind lun31
      ZTR=.TRUE.
      ZFA=.FALSE.

      do i=1,8
         e2scc(i)=0.d0
         psiscc(i)=0.d0
      enddo
      E2SC(1)=0.d0
      E2PC(1)=0.d0
      call nwadd(nwall,norb*(8+1+4+1))
      allocate(fmpb(norb))

      e(1)=escf

      call nwadd(nwall,nact**2*8)
      allocate(da(1:nact,1:nact))
      read(lun31)
      read(lun31)((da(k,l),k=1,nact),l=1,nact)

C
C     Let us compute the inactive orbitals energies
C
      call nwadd(nwall,(norb*(norb+1)/2)*8)
      allocate(heffd(norb*(norb+1)/2))
      call getorbeps(WRK,LWRK,onel,heffd,fmpb,da,norb,nact,
     *               ncore,id2n,in2d)
c
c     Debug printing of the frozen core orbitals
c
       allocate (igelo(norb))
       ngel=0
       do i=1,norb
       if (zgel(i)) then
        ngel=ngel+1
        igelo(ngel)=i
       endif
       enddo
      if (ngel.gt.0) write (lupri,271) ngel,(igelo(i),i=1,ngel)
  271 format(/i5,' FROZEN CORE ORBITALS:'/10(5x,18i4/))
      write (lupri,'(/)')
      call flshfo(lupri)

cele
      LUINTD = -1000
      CALL GPOPEN(LUINTD,'DRCINT','UNKNOWN',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
cele

c
c
c     V^0 term (double inactive)
c
      ncoat=ncore+nact
      nvirt=norb-ncoat     
      if(ncore.gt.0.and.nvirt.gt.0)then
         write (lupri,*)'V(0) class'
         call v0(luintd,norb,ncore,nact,in2d,id2n,itsym,its,wrk,lwrk,
     *        fmpb,zgel,lupri,nwall)
         write (lupri,125)'E2SC =',e2scc(1)
         write (lupri,125)'E2PC =',e2pcc(1)
         write (lupri,125)'PSISC=',psiscc(1)
         write (lupri,125)'PSIPC=',psipcc(1)
         t0=second()
         write (lupri,'(a,f8.2,a)')' Elapsed time ',t0-t,' sec.'
         write (lupri,*)
      endif
c
c     V^(+1) term
c
  125 format(1x,a,f12.8)
      if(ncore.gt.0.and.nvirt.gt.0)then
         write (lupri,*)'V(+1) class'
         call nwadd(nwall,nact**2*8)
         allocate(gkp1a(1:nact,1:nact))
         read(lun31)((gkp1a(k,l),k=1,nact),l=1,nact)
         call v1k(luintd,da,gkp1a,norb,ncore,nact,in2d,id2n,itsym,its,
     *          wrk,lwrk,fmpb,zgel,lupri,nwall)
         call nwdel(nwall,nact**2*8)
         deallocate (gkp1a)
         write (lupri,125)'E2SC =',e2scc(2)
         write (lupri,125)'E2PC =',e2pcc(2)
         write (lupri,125)'PSISC=',psiscc(2)
         write (lupri,125)'PSIPC=',psipcc(2)
         tp1=second()
         write (lupri,'(a,f8.2,a)')' Elapsed time ',tp1-t0,' sec.'
         write (lupri,*)
      else
         read(lun31)
      endif

      CALL GPCLOSE(LUINTD,'KEEP')
c
c     V^(-1) term
c
      if(ncore.gt.0.and.nvirt.gt.0)then
         write (lupri,*)'V(-1) class'
         call nwadd(nwall,nact**2*8)
         allocate(gkm1a(1:nact,1:nact))
         read(lun31)((gkm1a(k,l),k=1,nact),l=1,nact)
         call vm1k(da,gkm1a,norb,ncore,nact,in2d,id2n,itsym,its,
     *        wrk,lwrk,fmpb,zgel,lupri,nwall)
         call nwdel(nwall,nact**2*8)
         deallocate (gkm1a)
         write (lupri,125)'E2SC =',e2scc(3)
         write (lupri,125)'E2PC =',e2pcc(3)
         write (lupri,125)'PSISC=',psiscc(3)
         write (lupri,125)'PSIPC=',psipcc(3)
         tm1=second()
         write (lupri,'(a,f8.2,a)')' Elapsed time ',tm1-tp1,' sec.'
         write (lupri,*)
      else
         read(lun31)
      endif
c
c     V^(+2) term
c
      if(ncore.gt.0)then
         write(lupri,*)'V(+2) class'
         call nwadd(nwall,2*nact**4*8)
         allocate(gkp2aa(1:nact,1:nact,1:nact,1:nact))
         allocate(daat(1:nact,1:nact,1:nact,1:nact))
         do i=1,nact
            do j=1,nact
               read(lun31)((daat(i,j,k,l),k=1,nact),l=1,nact)
            enddo
         enddo
         do i=1,nact
            do j=1,nact
               read(lun31)((gkp2aa(i,j,k,l),k=1,nact),l=1,nact)
            enddo
         enddo
         call v2mod(gkp2aa,daat,norb,ncore,nact,in2d,id2n,itsym,its,
     *        wrk,lwrk,fmpb,zgel,lupri,nwall)
         call nwdel(nwall,2*nact**4*8)
         deallocate (gkp2aa)
         deallocate (daat)
         write(lupri,125)'E2SC =',e2scc(4)
         write(lupri,125)'E2PC =',e2pcc(4)
         write(lupri,125)'PSISC=',psiscc(4)
         write(lupri,125)'PSIPC=',psipcc(4)
         tp2=second()
         write (lupri,'(a,f8.2,a)')' Elapsed time ',tp2-tm1,' sec.'
         write(lupri,*)
      else
         do i=1,nact
            do j=1,nact
               read(lun31)
               read(lun31)
            enddo
         enddo
      endif
c     
c     V^(-2) term
c
      allocate(daa(1:nact,1:nact,1:nact,1:nact))
      do i=1,nact
      do j=1,nact
      read(lun31)((daa(i,j,k,l),k=1,nact),l=1,nact)
      enddo
      enddo
      if(nvirt.gt.0)then
         write (lupri,*)'V(-2) class'
         call nwadd(nwall,2*nact**4*8)
         allocate(gkm2aa(1:nact,1:nact,1:nact,1:nact))
         do i=1,nact
            do j=1,nact
               read(lun31)((gkm2aa(i,j,k,l),k=1,nact),l=1,nact)
            enddo
         enddo
         call vm2mod(daa,gkm2aa,norb,ncore,nact,in2d,id2n,itsym,its,
     *        wrk,lwrk,fmpb,zgel,lupri,nwall)
         call nwdel(nwall,nact**4*8)
         deallocate (gkm2aa)
         write (lupri,125)'E2SC= ',e2scc(5)
         write (lupri,125)'E2PC= ',e2pcc(5)
         write (lupri,125)'PSISC=',psiscc(5)
         write (lupri,125)'PSIPC=',psipcc(5)
         tm2=second()
         write (lupri,'(a,f8.2,a)')' Elapsed time ',tm2-tp2,' sec.'
         write (lupri,*)
      else
         do i=1,nact
            do j=1,nact
               read(lun31)
            enddo
         enddo
      endif
c     
c     V^(0p) term
c
      if(ncore.gt.0.and.nvirt.gt.0)then
         write(lupri,*)'V(0'') class'
         call nwadd(nwall,2*nact**4*8)
         allocate(gk0pa(1:nact,1:nact,1:nact,1:nact))
         allocate(gk0pd(1:nact,1:nact,1:nact,1:nact))
         call nwadd(nwall,nact**2*8)
         allocate(gk0pf(1:nact,1:nact))
         do i=1,nact
            do j=1,nact
               read(lun31)((gk0pa(i,j,k,l),k=1,nact),l=1,nact)
               read(lun31)((gk0pd(i,j,k,l),k=1,nact),l=1,nact)
            enddo
         enddo
         read(lun31)((gk0pf(k,l),k=1,nact),l=1,nact)
         call v0pmod(da,daa,gk0pa,gk0pd,gk0pf,heffd,norb,ncore
     $        ,nact,in2d,id2n,itsym,its,wrk,lwrk,fmpb,zgel,lupri,nwall)
         call nwdel(nwall,2*nact**4*8)
         deallocate (gk0pa)
         deallocate (gk0pd)
         call nwdel(nwall,nact**2*8)
         deallocate (gk0pf)
         write(lupri,125)'E2SC =',e2scc(8)
         write(lupri,125)'E2PC =',e2pcc(8)
         write(lupri,125)'PSISC=',psiscc(8)
         write(lupri,125)'PSIPC=',psipcc(8)
         t0p=second()
         write (lupri,'(a,f8.2,a)')' Elapsed time ',t0p-tm2,' sec.'
         write (lupri,*)
         call flshfo(lupri)
      else
         do i=1,nact
            do j=1,nact
               read(lun31)
               read(lun31)
            enddo
         enddo
         read(lun31)
      endif
c
c     V^(-1p) term
c


      call nwadd(nwall,2*nact**6*8)
      allocate(ro3(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
      allocate(amat(1:nact,1:nact,1:nact,1:nact,1:nact,1:nact))
      call nwadd(nwall,2*nact**4*8)
      allocate(bmat(1:nact,1:nact,1:nact,1:nact))
      allocate(cmat(1:nact,1:nact,1:nact,1:nact))
      call nwadd(nwall,nact**2*8)
      allocate(dmat(1:nact,1:nact))
c-----Reading of the 3 part spinless matrix and amat matrix
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  read(lun31)((ro3(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
                  read(lun31)((amat(i,j,k,l,ip,jp),ip=1,nact),jp=1,nact)
               enddo
            enddo
         enddo
      enddo
c----Reading of the bmat, cmat, btmat and ctmat matrices
      do i=1,nact
         do j=1,nact
            read(lun31)((bmat(i,j,k,l),k=1,nact),l=1,nact)
            read(lun31)((cmat(i,j,k,l),k=1,nact),l=1,nact)
         enddo
      enddo
c----Reading of the dmat and dtmat matrices
      read(lun31)((dmat(i,j),i=1,nact),j=1,nact)
      if(nvirt.gt.0)then
         write (lupri,125)'V(-1)'' class'
         call vmupE(ro3,amat,bmat,cmat,dmat,daa,da,heffd,norb,ncore,nact
     $        ,in2d,id2n,itsym,its,wrk,lwrk,fmpb,zgel,lupri,nwall)
         write (lupri,125)'E2SC =',e2scc(7)
         write (lupri,125)'E2PC =',e2pcc(7)
         write (lupri,125)'PSISC=',psiscc(7)
         write (lupri,125)'PSIPC=',psipcc(7)
         tmp=second()
         write (lupri,'(a,f8.2,a)')' Elapsed time ',tmp-t0p,' sec.'
         write (lupri,125)
      endif

c
c     V^(+1p) term
c

      if(ncore.gt.0)then
         write (lupri,*)'V(+1)'' class'
c-----Reading of the matrices for the V+1' term
         do i=1,nact
            do j=1,nact
               do k=1,nact
                  do l=1,nact
                     read(lun31)((amat(i,j,k,l,ip,jp),ip=1,nact),jp=1
     $                    ,nact)
                  enddo
               enddo
            enddo
         enddo
c---- Reading of the matrices bmat and cmat matrices
         do i=1,nact
            do j=1,nact
               read(lun31)((bmat(i,j,k,l),k=1,nact),l=1,nact)
               read(lun31)((cmat(i,j,k,l),k=1,nact),l=1,nact)
            enddo
         enddo
c---- Reading of the dmat matrix
         read(lun31)((dmat(i,j),i=1,nact),j=1,nact)
         call vpupE(ro3,amat,bmat,cmat,dmat,daa,da,heffd,norb,ncore,nact
     $        ,in2d,id2n,itsym,its,wrk,lwrk,fmpb,zgel,lupri,nwall)
         write (lupri,125)'E2SC =',e2scc(6)
         write (lupri,125)'E2PC =',e2pcc(6)
         write (lupri,125)'PSISC=',psiscc(6)
         write (lupri,125)'PSIPC=',psipcc(6)
         tpp=second()
         write (lupri,'(a,f8.2,a)')' Elapsed time ',tpp-tmp,' sec.'
         write (lupri,*)
c
      endif
      call nwdel(nwall,(norb*(norb+1)/2)*8)
      deallocate (heffd)
      write (lupri,*) ' Strongly Contracted and Partially Contracted ',
     * 'NEVPT2 results for each class'
      write (lupri,*) 
      write (lupri,'(10x,a)')
     $      '   Norm  SC       Energy SC      Norm  PC       Energy PC '
      write (lupri,123) ' (0) ',psiscc(1),e2scc(1),psipcc(1),e2pcc(1)
      write (lupri,123) '(+1) ',psiscc(2),e2scc(2),psipcc(2),e2pcc(2)
      write (lupri,123) '(-1) ',psiscc(3),e2scc(3),psipcc(3),e2pcc(3)
      write (lupri,123) '(+2) ',psiscc(4),e2scc(4),psipcc(4),e2pcc(4)
      write (lupri,123) '(-2) ',psiscc(5),e2scc(5),psipcc(5),e2pcc(5)
      write (lupri,123) "(+1)'",psiscc(6),e2scc(6),psipcc(6),e2pcc(6)
      write (lupri,123) "(-1)'",psiscc(7),e2scc(7),psipcc(7),e2pcc(7)
      write (lupri,123) " (0)'",psiscc(8),e2scc(8),psipcc(8),e2pcc(8)
      write (lupri,'(2x,65a1)')('-',i=1,65)
      write (lupri,123) 'Total',psisc(1),e2sc(1),psipc(1),e2pc(1)
  123 format (2x,A,4f15.10)
      
      deallocate(ro3)
      deallocate(amat)
      deallocate(bmat)
      deallocate(cmat)
      deallocate(dmat)
      deallocate(daa)
cro      if(zro)then
cro         allocate(rho(1:norb,1:norb))
cro         call dzero(rho,norb**2)
cro         do i=1,nact
cro            do j=1,nact
cro               rho(ncore+i,ncore+j)=da(i,j)
cro            enddo
cro         enddo
cro         do i=1,ncore
cro            rho(i,i)=2.d0
cro         enddo
cro         rewind 27
cro         read(27)((rho(i,j),i=1,ncore),j=ncore+nact+1,norb)
cro         do i=1,ncore
cro            do j=ncore+nact+1,norb
cro               rho(j,i)=rho(i,j)
cro            enddo
cro         enddo
cro         read(27)((rho(i,j),i=ncore+1,ncore+nact),j=ncore+nact+1,norb)
cro         do i=ncore+1,ncore+nact
cro            do j=ncore+nact+1,norb
cro               rho(j,i)=rho(i,j)
cro            enddo
cro         enddo
cro         read(27)((rho(i,j),i=1,ncore),j=ncore+1,ncore+nact)
cro         do i=1,ncore
cro            do j=ncore+1,ncore+nact
cro               rho(j,i)=rho(i,j)
cro            enddo
cro         enddo
cro         if(zverbose)then
cro            print*,'Density matrix 0+1'
cro            call matout(norb,norb,rho,norb,1.d0)
cro         endif
cro         print*,'calculation of natural occupation numbers'
cro         allocate(w(norb))
cro         allocate(work(3*norb))
croc         call dsyev('V','U',norb,rho,norb,w,work,3*norb,info)
cro         print*,'core occupation numbers:'
cro         print '(6f12.8)',(w(i),i=norb-ncore+1,norb)
cro         print*,'active occupation numbers:'
cro         print '(6f12.8)',(w(i),i=norb-ncore-nact+1,norb-ncore)
cro         print*,'virtual occupation numbers:'
cro         print '(6f12.8)',(w(i),i=1,norb-ncore-nact)
cro      endif
      WRITE (lupri,'(/)')
      WRITE (lupri,2000)
      WRITE (lupri,1060)
      WRITE (lupri,2000)
      WRITE (lupri,'(/)')
c      WRITE (lupri,1105) ESCF
      WRITE (lupri,1070)
      E2SC(1)=E2SC(1)+E(1)
      E2PC(1)=E2PC(1)+E(1)
      WRITE (lupri,1100) 1,E(1),E2SC(1),E2PC(1)
      write (lupri,'(///a,f8.2,a/)')' Total elapsed time in dypc',
     *  tpp-t,' sec.'
 1060 FORMAT (5X,'TOTAL ENERGY CONTRIBUTIONS')
 1070 FORMAT (' STATE',3X,'DIAGONALISATION        PERTURBATION SC-D',
     * '        PERTURBATION PC-D')
 1100 FORMAT (I5,F17.8,3F24.8)
c 1105 FORMAT (' REFERENCE: ',F16.8,//)
 2000 FORMAT (5X,32('*'))
      CALL QEXIT('DYPC')
      return
      END
c----------------------------------------------------------------------
      subroutine v0(luintd,norb,ncore,nact,in2d,id2n,itsym,its,wrk,lwrk,
     * fmpb,zgel,lupri,nwall)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension in2d(*),id2n(*),itsym(*),its(8,8),needtp(6)
      COMMON /NEVPT/ e(1),e2pc(1),e2sc(1),psisc(1),psipc(1),
     *               e2pcc(8),psipcc(8),e2scc(8),psiscc(8),
     *               ntotdet,ntotconf,ntotcap,zro
      dimension WRK(LWRK)
      allocatable h2m(:,:),tint(:,:)
      dimension fmpb(*),zgel(*)
cele
#include <dummy.h>
cele

      
      ncoat=ncore+nact
      nvirt=norb-ncoat     


      do i=1,6
       needtp(i)=0
      enddo


celec------ two-electron part with Mulliken distr (ir,ii|**)
c------ two-electron part with Dirac distr (ii,ij|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
cele      needtp(4)=1
      needtp(1)=1
cele      icount=0
      call nwadd(nwall,norb*norb*8)
      allocate(h2m(norb,norb))
      call nwadd(nwall,ncore*ncore*nvirt*8)
cele      allocate(tint(ncore,ncore,nvirt))
      allocate(tint(nvirt,nvirt))
      call dzero(tint,nvirt*nvirt)
cele  900 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
  900 CALL NXTH2D(ICD,IDD,H2M,NEEDTP,LUINTD,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1000
cele
      ii=id2n(icd)
      ij=id2n(idd)
cele
      if (zgel(ii)) goto 900
      if (zgel(ij)) goto 900
cele      ic=id2n(icd)
cele      id=id2n(idd)
Cele      if (ic.gt.id) then
Cele       ir=ic
Cele       ii=id
Cele      else
Cele       ir=id
Cele       ii=ic
Cele      endif
cele
cele
cele      icount=icount+1
cele      do ij=1,ncore
cele       do is=ncore+nact+1,norb
cele        ijd=in2d(ij)
cele        isd=in2d(is)
cele        tint(ii,ij,is-ncoat)=h2m(ijd,isd)
cele       enddo
cele      enddo
      do ir=ncore+nact+1,norb
       do is=ir,norb
        ird=in2d(ir)
        isd=in2d(is)
        tint(ir-ncoat,is-ncoat)=h2m(ird,isd)
        tint(is-ncoat,ir-ncoat)=h2m(isd,ird)
       enddo
      enddo
cele      if (icount.lt.ncore) then
cele       goto 900
cele      else
cele       icount=0
cele       goto 11
cele      endif

cele  11  irsy=itsym(ir)
cele      do 800 ii=1,ncore
cele      if (zgel(ii)) goto 800
cele      isy=itsym(ii)
cele      do 700 ij=ii,ncore
cele      if (zgel(ij)) goto 700
cele      jsy=itsym(ij)
cele      ijsy=its(isy,jsy)
cele      irijsy=its(ijsy,irsy)
 
cele
      isy=itsym(ii)
      jsy=itsym(ij)
      ijsy=its(isy,jsy)
      do ir=ncore+nact+1,norb
         irsy=itsym(ir)
         irijsy=its(ijsy,irsy)
cele

      do is=ir,norb
         issy=itsym(is)
         irsijsy=its(irijsy,issy)
         if(irsijsy.ne.1)goto 500
      

         zrs=ir.eq.is
         zij=ii.eq.ij
cele      ab1=tint(ij,ii,is-ncoat)
cele      ab2=tint(ii,ij,is-ncoat)
         ab1=tint(ir-ncoat,is-ncoat)
         ab2=tint(is-ncoat,ir-ncoat)
cele
      
         if (zrs.and.zij) dint=ab1
         if (zrs.neqv.zij) dint=sqrt(2.d0)*ab1
         if (.not.(zrs.or.zij)) dint=2.d0*sqrt(ab1**2+ab2**2-ab1*ab2)
         co=-dint/(fmpb(ir)+fmpb(is)-fmpb(ii)-fmpb(ij))
         e2sc(1)=e2sc(1)+dint*co
         psisc(1)=psisc(1)+co*co
         e2scc(1)=e2scc(1)+dint*co
         psiscc(1)=psiscc(1)+co*co
  500 continue
      enddo

cele
      enddo
cele

cele  700 enddo

cele  800 enddo


      e2pcc(1)=e2scc(1)
      psipcc(1)=psiscc(1)
      e2pc(1)=e2sc(1)
      psipc(1)=psisc(1)

      goto 900
 1000 continue

      call nwdel(nwall,norb*norb*8)
      deallocate(h2m)
      call nwdel(nwall,ncore*ncore*nvirt*8)
      deallocate(tint)

      return
      end
c----------------------------------------------------------------
      subroutine v1k(luintd,da,gkp1a,norb,ncore,nact,in2d,id2n,itsym,
     *               its,wrk,lwrk,fmpb,zgel,lupri,nwall)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension in2d(*),id2n(*),itsym(*),its(8,8),needtp(6)
      dimension da(nact,nact),gkp1a(nact,nact)
      COMMON /NEVPT/ e(1),e2pc(1),e2sc(1),psisc(1),psipc(1),
     *               e2pcc(8),psipcc(8),e2scc(8),psiscc(8),
     *               ntotdet,ntotconf,ntotcap,zro
      dimension WRK(LWRK)
      dimension fmpb(*),zgel(*)
      allocatable h2m(:,:),tintar(:,:),tintra(:,:)
      allocatable s(:,:),coef(:,:),w(:),work(:),caux(:,:),vct(:,:)
#include <dummy.h>
      
      ncoat=ncore+nact
      nvirt=norb-ncoat

c---- Initializating for partially contracted NEV-PT ------------------
      tini=second()
      call nwadd(nwall,(nact**2+4*nact)*8)
      allocate(s(1:nact,1:nact))
      allocate(coef(1:nact,1:nact))
      allocate(vct(1:nact,1:nact))
      allocate(w(nact))
      allocate(work(3*nact))
      call fill2(s,da,nact)
      call fill(coef,gkp1a,nact**2)
      call halve(coef,nact**2)
      call rs(nact,nact,s,w,1,vct,work,work(nact+1),ierr)  !eispack
c      call dsyev('V','U',nact,s,nact,w,work,3*nact,info)  !lapack
c      call dsygv(1,'V','U',nact,coef,nact,s,nact,w,work,3*nact,info)
c      if(info.ne.0)then
c         print*,'something wrong in diagonalization: sub v1k'
c         stop 'v1k'
c      endif
c      thr=1.d-8                 !perhaps too strict
      thr=1.d-6
      ndep=0
      do i=1,nact
       if(w(i).lt.0.d0)then
        ndep=ndep+1
        write (lupri,'(a,i3,a,f12.8)')'negative eigenvalue: eps(',i,')='
     $           ,w(i)
        elseif(w(i).lt.thr)then
            ndep=ndep+1
        endif
      enddo
      if (ndep.gt.0) then
      write (lupri,*)'Linear dependencies analysis:1'
      write (lupri,'(a,i5,a)')' There are',ndep,' linear dependencies'
      endif
      nnew=nact-ndep
c--   costruzione della matrice di trasformazione
      call nwadd(nwall,(nact*nnew)*8)
      allocate (caux(1:nact,1:nnew))
      call caux0p(caux,vct,w,nact,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact,nnew,1.d0,coef,nact,caux,nact,0.d0,s
     $     ,nact)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact,1.d0,caux,nact,s,nact,0.d0
     $     ,coef,nact)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
c      call dsyev('V','U',nnew,coef,nact,w,work,3*nact,info)
      call rs(nact,nnew,coef,w,1,vct,work,work(nact+1),ierr)  !eispack
c--   back transformation
c--   moltiplico CAUX(nact,nnew) per coef(nnew,nnew)
c--   e metto in S(nact,nnew)
      call dgemm('N','N',nact,nnew,nnew,1.d0,caux,nact,vct,nact,0
     $     .d0,s,nact)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact,nnew)
      call nwdel(nwall,(nact*nnew)*8)
      deallocate(caux)
      deallocate(vct)
c---costruzione matrice S
      do ia=1,nact
         do mu=1,nnew
            s(ia,mu)=0.d0
            do iap=1,nact
               if(iap.ne.ia)then
                  ro=-da(iap,ia)
               else
                  ro=2.d0-da(iap,ia)
               endif
               s(ia,mu)=s(ia,mu)+ro*coef(iap,mu)
            enddo
         enddo
      enddo


      do i=1,6
       needtp(i)=0
      enddo

c------ two-electron part with Dirac distr (ii,ij|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
      needtp(1)=1
      call nwadd(nwall,norb*norb*8)
      allocate(h2m(norb,norb))
      call nwadd(nwall,2*nact*nvirt*8)
      allocate(tintar(nact,nvirt))
      allocate(tintra(nvirt,nact))
      call dzero(tintar,nact*nvirt)
      call dzero(tintra,nact*nvirt)
  900 CALL NXTH2D(ICD,IDD,H2M,NEEDTP,LUINTD,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1000
      ii=id2n(icd)
      ij=id2n(idd)
      if (zgel(ii)) goto 900
      if (zgel(ij)) goto 900
      do ir=ncore+nact+1,norb
       do ia=ncore+1,ncore+nact
        ird=in2d(ir)
        iad=in2d(ia)
        tintar(ia-ncore,ir-ncoat)=h2m(iad,ird)
        tintra(ir-ncoat,ia-ncore)=h2m(ird,iad)
       enddo
      enddo
c
C --------- Strongly Contracted -----------------
      isy=itsym(ii)
      jsy=itsym(ij)
      ijsy=its(isy,jsy)

      do is=ncore+nact+1,norb
       issy=itsym(is)
       isijsy=its(ijsy,issy)

c       write (lupri,*) 'Eccitazione con gli indici',ii,ij,' -> ',is
c       call flshfo(lupri)

       dnormk=0.d0
       epsimk=0.d0
       do ia=1,nact
        iaa=ia+ncore
        do iap=1,nact
         iaap=iap+ncore
         contj=tintra(is-ncoat,ia)
         contjp=tintra(is-ncoat,iap)
         if (ij.eq.ii) then
          cint=0.5d0*contj*contjp
         else
          contk=tintar(ia,is-ncoat)
          contkp=tintar(iap,is-ncoat)
          cint=contj*contjp+contk*contkp-
     *       0.5d0*contj*contkp-0.5d0*contk*contjp
         endif
         if (ia.eq.iap) then
          ro=2.d0-da(ia,iap)
         else
          ro=-da(ia,iap)
         endif
         dnormk=dnormk+2.d0*ro*cint
         epsimk=epsimk+gkp1a(ia,iap)*cint
        enddo
       enddo

       if(abs(dnormk).lt.1.d-15)goto 600    
       if(dnormk.lt.0.d0)then
         write(lupri,*)'Warning v1k: dnormk=',dnormk
         dnormk=0.d0
       endif
       if (dnormk.eq.0.d0) goto 600
       epsimk=epsimk/dnormk
       den=fmpb(ii)+fmpb(ij)-fmpb(is)-epsimk
       co=dnormk/den
       con=co/den
       e2sc(1)=e2sc(1)+co
       psisc(1)=psisc(1)+con
       e2scc(2)=e2scc(2)+co
       psiscc(2)=psiscc(2)+con

  600 continue
      enddo

c  700 enddo

c  800 enddo

C --------- Partially Contracted -----------------
      call e2p1cont(ii,ij,tintar,tintra,e2,psi2,coef,s,w,da,fmpb,
     *              nact,ncore,norb,nnew,zgel)
      e2pcc(2)=e2pcc(2)+e2
      psipcc(2)=psipcc(2)+psi2
      e2pc(1)=e2pc(1)+e2
      psipc(1)=psipc(1)+psi2

      goto 900
 1000 continue

      call nwdel(nwall,(nact**2+4*nact)*8)
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      call nwdel(nwall,norb*norb*8)
      deallocate(h2m)
      call nwdel(nwall,2*nact*nvirt*8)
      deallocate(tintar)
      deallocate(tintra)

      return
      end
c----------------------------------------------------------------
      subroutine fill2 (s,da,nact)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension s(nact,nact),da(nact,nact)
      do i=1,nact
         do j=1,nact
            if(i.ne.j)then
               s(i,j)=-da(i,j)
            else
               s(i,i)=2.d0-da(i,i)
            endif
         enddo
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine vm1k(da,gkm1a,norb,ncore,nact,in2d,id2n,itsym,its,
     *          wrk,lwrk,fmpb,zgel,lupri,nwall)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension in2d(*),id2n(*),itsym(*),its(8,8),needtp(6)
      dimension da(nact,nact),gkm1a(nact,nact)
      COMMON /NEVPT/ e(1),e2pc(1),e2sc(1),psisc(1),psipc(1),
     *               e2pcc(8),psipcc(8),e2scc(8),psiscc(8),
     *               ntotdet,ntotconf,ntotcap,zro
      dimension WRK(LWRK)
      dimension fmpb(*),zgel(*)
      allocatable h2m(:,:),tint(:,:,:,:)
      allocatable s(:,:),coef(:,:),w(:),work(:),caux(:,:),vct(:,:)
      
      ncoat=ncore+nact
      nvirt=norb-ncoat

c---- Initializating for partially contracted NEV-PT ------------------
      tini=second()
      call nwadd(nwall,(nact**2+4*nact)*8)
      allocate(s(1:nact,1:nact))
      allocate(coef(1:nact,1:nact))
      allocate(vct(1:nact,1:nact))
      allocate(w(nact))
      allocate(work(3*nact))
      call fill(s,da,nact**2)
      call fill(coef,gkm1a,nact**2)
      call halve(coef,nact**2)
      call rs(nact,nact,s,w,1,vct,work,work(nact+1),ierr)  !eispack
c      call dsyev('V','U',nact,s,nact,w,work,3*nact,info)
c      call dsygv(1,'V','U',nact,coef,nact,s,nact,w,work,3*nact,info)
c      if(info.ne.0)then
c         print*,'something wrong in diagonalization: sub vm1k'
c         stop 'vm1k'
c      endif
c      thr=1.d-8                 !perhaps too strict
      thr=1.d-6
      ndep=0
      do i=1,nact
       if(w(i).lt.0.d0)then
        ndep=ndep+1
        write (lupri,'(a,i3,a,f12.8)')'negative eigenvalue: eps(',i,')='
     $           ,w(i)
        elseif(w(i).lt.thr)then
            ndep=ndep+1
        endif
      enddo
      if (ndep.gt.0) then
      write (lupri,*)'Linear dependencies analysis:2'
      write (lupri,'(a,i5,a)')' There are',ndep,' linear dependencies'
      endif
      nnew=nact-ndep
c--   costruzione della matrice di trasformazione
      call nwadd(nwall,(nact*nnew)*8)
      allocate (caux(1:nact,1:nnew))
      call caux0p(caux,vct,w,nact,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact,nnew,1.d0,coef,nact,caux,nact,0.d0,s
     $     ,nact)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact,1.d0,caux,nact,s,nact,0.d0
     $     ,coef,nact)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
c      call dsyev('V','U',nnew,coef,nact,w,work,3*nact,info)
      call rs(nact,nnew,coef,w,1,vct,work,work(nact+1),ierr)  !eispack
c--   back transformation
c--   moltiplico CAUX(nact,nnew) per coef(nnew,nnew)
c--   e metto in S(nact,nnew)
      call dgemm('N','N',nact,nnew,nnew,1.d0,caux,nact,vct,nact,0
     $     .d0,s,nact)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact,nnew)
      call nwdel(nwall,(nact*nnew)*8)
      deallocate(caux)
      deallocate(vct)

      do i=1,6
       needtp(i)=0
      enddo

c------ two-electron part with mulliken distr (ir,ii|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
      needtp(4)=1
      L_tint=nvirt*ncore*nvirt*nact
      call nwadd(nwall,norb*norb*8)
      allocate(h2m(norb,norb))
      call nwadd(nwall,L_tint*8)
      allocate(tint(nvirt,ncore,nvirt,nact))
      call dzero(tint,L_tint)
  900 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1000
      ic=id2n(icd)
      id=id2n(idd)
      if (ic.gt.id) then
       ir=ic
       ii=id
      else
       ir=id
       ii=ic
      endif
      irc=ir-ncoat
      do ia=1,nact
       iad=in2d(ia+ncore)
       do is=1,nvirt
        isd=in2d(is+ncoat)
        tint(irc,ii,is,ia)=h2m(isd,iad)
       enddo
      enddo
      goto 900
 1000 continue


c
C --------- Strongly Contracted -----------------
      do ii=1,ncore ! do 800
      if (zgel(ii)) goto 800
      isy=itsym(ii)

      do irc=1,nvirt ! do 700
      ir=irc+ncoat
      irsy=itsym(ir)
      irisy=its(isy,irsy)

      do isc=irc,nvirt ! do 600
      is=isc+ncoat
      issy=itsym(is)
      irsisy=its(irisy,issy)

      irisy=its(isy,irsy)

      irsisy=its(irisy,issy)

      dnormk=0.d0
      epsimk=0.d0
      do ia=1,nact
       iaa=ia+ncore
       contj=tint(irc,ii,isc,ia)
       contk=tint(isc,ii,irc,ia)
       do iap=1,nact
        iaap=iap+ncore
        contjp=tint(irc,ii,isc,iap)
       if (ir.eq.is) then
        cint=.5d0*contj*contjp
       else
        contkp=tint(isc,ii,irc,iap)
        cint=contj*contjp+contk*contkp-
     *       0.5d0*contj*contkp-0.5d0*contk*contjp
       endif
        dnormk=dnormk+2.d0*da(ia,iap)*cint
        epsimk=epsimk+gkm1a(ia,iap)*cint
       enddo ! iap
      enddo ! ia

      if (abs(dnormk).lt.1.d-15) goto 600
      if(dnormk.lt.0.d0)then
         write(lupri,*)'Warning vm1k: dnormk=',dnormk
         dnormk=0.d0
      endif

      if (dnormk.ne.0.d0) then
        epsimk=-epsimk/dnormk
        den=epsimk+fmpb(ii)-fmpb(ir)-fmpb(is)
        co=dnormk/den
        con=co/den
        e2sc(1)=e2sc(1)+co
        psisc(1)=psisc(1)+con
        e2scc(3)=e2scc(3)+co
        psiscc(3)=psiscc(3)+con
      endif



  600 continue
      enddo 

  700 continue
      enddo

  800 continue
      enddo
c----Partially contracted NEV-PT---------------------------------
      call e2m1cont(tint,e2,psi2,coef,s,w,da,fmpb,nact,ncore,norb,nnew
     $     ,zgel)
      e2pcc(3)=e2
      psipcc(3)=psi2
      e2pc(1)=e2pc(1)+e2
      psipc(1)=psipc(1)+psi2


      call nwdel(nwall,(nact**2+4*nact)*8)
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      call nwdel(nwall,norb*norb*8)
      deallocate(h2m)
      call nwdel(nwall,ncore*nvirt*nact*8)
      deallocate(tint)
      return
      end
c-------------------------------------------------------------
      subroutine fill(a,b,n)
#include <implicit.h>
      DIMENSION a(*),b(*)
      do i=1,n
         a(i)=b(i)
      enddo
      return
      end
c----------------------------------------------
      subroutine halve(a,n)
#include <implicit.h>
      DIMENSION a(*)
      do i=1,n
         a(i)=0.5d0*a(i)
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2m1cont(tint,e2,psi2,coef,s,w,da,eps,nact,ncore,
     *                    norb,nnew,zgel)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension w(*),eps(*),coef(nact,nact),s(nact,nact),da(nact,nact)
      dimension zgel(*),tint(norb-ncore-nact,ncore,norb-ncore-nact,nact)
      integer a,ap
      ncoat=ncore+nact
c---costruzione matrice S
      do a=1,nact
         do mu=1,nnew
            s(a,mu)=0.d0
            do ap=1,nact
               s(a,mu)=s(a,mu)+da(ap,a)*coef(ap,mu)
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do i=1,ncore
        if (zgel(i)) goto 1
         do ir=ncore+nact+1,norb
         irc=ir-ncoat
            do is=ir,norb
            isc=is-ncoat
               do mu=1,nnew
                  rsmui=0.d0
                  rsmuip=0.d0
                  do a=1,nact
                     if(ir.ne.is)then
                        aint=tint(irc,i,isc,a)+tint(isc,i,irc,a)
                        bint=tint(irc,i,isc,a)-tint(isc,i,irc,a)
                        rsmui=rsmui+aint*s(a,mu)
                        rsmuip=rsmuip+bint*s(a,mu)
                     else
                        aint=tint(irc,i,irc,a)
                        rsmui=rsmui+aint*s(a,mu)
                     endif
                  enddo
                  deno=w(mu)+eps(ir)+eps(is)-eps(i)
                  if(ir.ne.is)then
                     contri=(0.5d0*rsmui**2+1.5d0*rsmuip**2)/deno
                     e2=e2-contri
                     psi2=psi2+contri/deno
                  else
                     contri=rsmui**2/deno
                     e2=e2-contri
                     psi2=psi2+contri/deno
                  endif
               enddo
            enddo
         enddo
 1    continue
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2m2cont(tint,e2,psi2,coef,s,w,daa,eps,nact,ncore,
     *                    norb,nnew,zgel)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension w(*),eps(*),coef(nact**2,nact**2),s(nact**2,nact**2)
     $     ,daa(nact,nact,nact,nact),
     *      tint(norb-nact-ncore,nact,norb-nact-ncore,nact)
      dimension zgel(*)
      integer a,ap,b,bp,ab,abp
      ncoat=ncore+nact
c--caso r#s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            do mu=1,nnew
               s(ab,mu)=0.d0
               abp=0
               do ap=1,nact
                  do bp=1,nact
                     abp=abp+1
                     s(ab,mu)=s(ab,mu)+daa(ap,bp,a,b)*coef(abp,mu)
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do ir=ncore+nact+1,norb
      irc=ir-ncoat
         do is=ir+1,norb
         isc=is-ncoat
            do mu=1,nnew
               rsmu=0.d0
               ab=0
               do a=1,nact
                  do b=1,nact
                     ab=ab+1
                     rsmu=rsmu+tint(irc,b,isc,a)*s(ab,mu)
                  enddo
               enddo
               deno=w(mu)+eps(ir)+eps(is)
               contri=rsmu**2/deno
               e2=e2-contri
               psi2=psi2+contri/deno
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2m2contp(tint,e2,psi2,coef,s,w,daa,eps,nact,ncore,
     *                     norb,nnew,zgel)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension w(*),eps(*),coef(nact*(nact+1)/2,nact*(nact+1)/2),s(nact
     $     *(nact+1)/2,nact*(nact+1)/2),daa(nact,nact,nact,nact),
     *     tint(norb-nact-ncore,nact,norb-nact-ncore,nact)
      dimension zgel(*)
      integer a,ap,b,bp,ab,abp
      ncoat=ncore+nact
c--caso r=s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=a,nact
            ab=ab+1
            do mu=1,nnew
               s(ab,mu)=0.d0
               abp=0
               do ap=1,nact
                  do bp=ap,nact
                     abp=abp+1
                     s(ab,mu)=s(ab,mu)+0.5d0*(daa(ap,bp,a,b)+daa(ap
     $                    ,bp,b,a))*coef(abp,mu)
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do ir=ncore+nact+1,norb
      irc=ir-ncoat
         do mu=1,nnew
            rsmu=0.d0
            ab=0
            do a=1,nact
               ab=ab+1
                  rsmu=rsmu+tint(irc,a,irc,a)*s(ab,mu)
               do b=a+1,nact
                  ab=ab+1
                  rsmu=rsmu+tint(irc,b,irc,a)*s(ab,mu)*2.d0
               enddo
            enddo
            deno=w(mu)+eps(ir)+eps(ir)
            contri=rsmu**2/deno
            e2=e2-contri
            psi2=psi2+contri/deno
         enddo
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2p2cont(tint,e2,psi2,coef,s,w,daa,eps,nact,ncore,
     *                    norb,nnew,zgel)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension w(*),eps(*),coef(nact**2,nact**2),s(nact**2,nact**2)
     $     ,daa(nact,nact,nact,nact),zgel(*),
     *      tint(nact,ncore,nact,ncore)
      integer a,ap,b,bp,ab,abp
c--caso r#s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            do mu=1,nnew
               s(ab,mu)=0.d0
               abp=0
               do ap=1,nact
                  do bp=1,nact
                     abp=abp+1
                     s(ab,mu)=s(ab,mu)+daa(ap,bp,a,b)*coef(abp,mu)
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do i=1,ncore
        if (zgel(i)) goto 1
         do j=i+1,ncore
           if (zgel(j)) goto 2
            do mu=1,nnew
               aijmu=0.d0
               ab=0
               do a=1,nact
                  do b=1,nact
                     ab=ab+1
                     aijmu=aijmu+tint(b,i,a,j)*s(ab,mu)
                  enddo
               enddo
               deno=w(mu)-eps(i)-eps(j)
               contri=aijmu**2/deno
               e2=e2-contri
               psi2=psi2+contri/deno
            enddo
 2       continue
         enddo
 1    continue
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine e2p2contp(tint,e2,psi2,coef,s,w,daa,eps,nact,ncore,
     *                     norb,nnew,zgel)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension w(*),eps(*),coef(nact*(nact+1)/2,nact*(nact+1)/2),s(nact
     $    *(nact+1)/2,nact*(nact+1)/2),daa(nact,nact,nact,nact),zgel(*),
     *     tint(nact,ncore,nact,ncore)
      integer a,ap,b,bp,ab,abp
c--caso r=s
c---costruzione matrice S
      ab=0
      do a=1,nact
         do b=a,nact
            ab=ab+1
            do mu=1,nnew
               s(ab,mu)=0.d0
               abp=0
               do ap=1,nact
                  do bp=ap,nact
                     abp=abp+1
                     s(ab,mu)=s(ab,mu)+0.5d0*(daa(ap,bp,a,b)+daa(ap
     $                    ,bp,b,a))*coef(abp,mu)
                  enddo
               enddo
            enddo
         enddo
      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do i=1,ncore
        if (zgel(i)) goto 1
         do mu=1,nnew
            aijmu=0.d0
            ab=0
            do a=1,nact
               ab=ab+1
              aijmu=aijmu+tint(a,i,a,i)*s(ab,mu)
               do b=a+1,nact
                  ab=ab+1
                  aijmu=aijmu+tint(b,i,a,i)*s(ab,mu)*2.d0
               enddo
            enddo
            deno=w(mu)-eps(i)-eps(i)
            contri=aijmu**2/deno
            e2=e2-contri
            psi2=psi2+contri/deno
         enddo
 1    continue
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine e2p1cont(ii,ij,tintar,tintra,e2,psi2,coef,s,w,da,eps,
     *                    nact,ncore,norb,nnew,zgel)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension w(*),eps(*),coef(nact,nact),s(nact,nact),da(nact,nact)
     $     ,zgel(*),tintar(nact,norb-nact-ncore),
     $      tintra(norb-nact-ncore,nact)
      integer a,ap
      logical*1 zcond
c---costruzione matrice S
c      do a=1,nact
c         do mu=1,nnew
c            s(a,mu)=0.d0
c            do ap=1,nact
c               if(ap.ne.a)then
c                  ro=-da(ap,a)
c               else
c                  ro=2.d0-da(ap,a)
c               endif
c               s(a,mu)=s(a,mu)+ro*coef(ap,mu)
c            enddo
c         enddo
c      enddo
c--------------------------------
      e2=0.d0
      psi2=0.d0
      do ir=ncore+nact+1,norb
       irc=ir-ncore-nact
        i=ii
        j=ij
c         do i=1,ncore
c          if (zgel(i)) goto 2
c            do j=i,ncore
c              if (zgel(j)) goto 3
               do mu=1,nnew
                  rjimu=0.d0
                  rjimup=0.d0
                  do a=1,nact
                     if(i.ne.j)then
                        aint=tintar(a,irc)+tintra(irc,a)
                        bint=tintar(a,irc)-tintra(irc,a)
                        rjimu=rjimu+aint*s(a,mu)
                        rjimup=rjimup+bint*s(a,mu)
                     else
                        aint=tintar(a,irc)
                        rjimu=rjimu+aint*s(a,mu)
                     endif
                  enddo
                  deno=w(mu)+eps(ir)-eps(i)-eps(j)
                  if(i.ne.j)then
                     contri=(0.5d0*rjimu**2+1.5d0*rjimup**2)/deno
                     e2=e2-contri
                     psi2=psi2+contri/deno
                  else
                     contri=rjimu**2/deno
                     e2=e2-contri
                     psi2=psi2+contri/deno
                  endif
               enddo
c 3          enddo
c 2       enddo
      enddo
      return
      end
c------------------------------------------------------
      subroutine v2mod(gkp2aa,daat,norb,ncore,nact,in2d,id2n,itsym,its,
     *          wrk,lwrk,fmpb,zgel,lupri,nwall)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension in2d(*),id2n(*),itsym(*),its(8,8),needtp(6)
c      dimension da(nact,nact),gkm1a(nact,nact)
      COMMON /NEVPT/ e(1),e2pc(1),e2sc(1),psisc(1),psipc(1),
     *               e2pcc(8),psipcc(8),e2scc(8),psiscc(8),
     *               ntotdet,ntotconf,ntotcap,zro
      dimension WRK(LWRK)
      dimension fmpb(*),zgel(*)
      dimension gkp2aa(nact,nact,nact,nact)
      dimension daat(nact,nact,nact,nact)
      allocatable h2m(:,:),tint(:,:,:,:)
      allocatable coef(:,:),s(:,:),w(:),work(:),caux(:,:),vct(:,:)
      
      ncoat=ncore+nact
      nvirt=norb-ncoat

      do i=1,6
       needtp(i)=0
      enddo

c------ two-electron part with mulliken distr (ia,ii|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
      needtp(2)=1
      L_tint = nact*ncore*nact*ncore
      call nwadd(nwall,norb*norb*8)
      allocate(h2m(norb,norb))
      call nwadd(nwall,L_tint*8)
      allocate(tint(nact,ncore,nact,ncore))
      call dzero(tint,L_tint)
  900 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1000
      ic=id2n(icd)
      id=id2n(idd)
      if(ic.gt.id) then
      ia=ic
      ii=id
      else
      ia=id
      ii=ic
      endif
      iac=ia-ncore
      do ib=1,nact 
      ibd=in2d(ib+ncore)
       do ij=1,ncore
        ijd=in2d(ij)
        tint(iac,ii,ib,ij)=h2m(ibd,ijd)
       enddo
      enddo
      goto 900
 1000 continue
      call nwdel(nwall,norb*norb*8)
      deallocate(h2m)

c
C --------- Strongly Contracted -----------------

      do ir=1,ncore ! do 1800
      if (zgel(ir)) goto 1800
      irsy=itsym(ir)

      do is=ir,ncore ! do 1700
      if (zgel(is)) goto 1700
      issy=itsym(is)
      irssy=its(irsy,issy)

c
c     Si inizia il ciclo sugli attivi 
c
      
      dn=0.d0
      dk=0.d0

      do ic=ncore+1,ncore+nact
         ica=ic-ncore
         do id=ncore+1,ncore+nact
            ida=id-ncore
            icsy=itsym(ic)
            idsy=itsym(id)
            icdsy=its(icsy,idsy)
            isimm=its(icdsy,irssy)
            if (isimm.ne.1) goto 1121
            do icp=ncore+1,ncore+nact
               icpa=icp-ncore
               do idp=ncore+1,ncore+nact
                  idpa=idp-ncore
                  icpsy=itsym(icp)
                  idpsy=itsym(idp)
                  icpdpsy=its(icpsy,idpsy)
                  isimm=its(icpdpsy,irssy)
                  if (isimm.ne.1) goto 1120
                  dum=daat(icpa,idpa,ica,ida)
                  dumk=gkp2aa(icpa,idpa,ica,ida)
                  dint1=tint(icpa,is,idpa,ir)
                  dint2=tint(ica,is,ida,ir)
                cint=dint1*dint2
                if (ir.eq.is) then
                cint=cint*0.25d0
                dum=dum+daat(icpa,idpa,ida,ica)
                dumk=dumk+gkp2aa(idpa,icpa,ica,ida)
                endif
                  dn=dn+cint*dum
                  dk=dk+cint*dumk
 1120          continue
               enddo
            enddo
 1121    continue
         enddo
      enddo
      if(abs(dn).lt.1.d-15) goto 1700
      if(dn.lt.0.d0)then
         write(lupri,*)'Warning v2mod: dn=',dn
         dn=0.0d0
      endif
      if(dn.gt.0.d0)then
         epsimk=dk/dn                     
         den=fmpb(ir)+fmpb(is)-epsimk
         co=dn/den
         con=co/den
         e2sc(1)=e2sc(1)+co
         psisc(1)=psisc(1)+con
         e2scc(4)=e2scc(4)+co
         psiscc(4)=psiscc(4)+con
      endif
      
 1700 continue
      enddo
 1800 continue
      enddo

c---Partially contracted NEV-PT-----------------------------------
      tini=second()
      call nwadd(nwall,(nact**4+4*nact**2)*8)
      allocate(s(1:nact**2,1:nact**2))
      allocate(coef(1:nact**2,1:nact**2))
      allocate(vct(1:nact**2,1:nact**2))
      allocate(w(nact**2))
      allocate(work(3*nact**2))
      call fillcm2(s,daat,nact)
      call fillcm2(coef,gkp2aa,nact)
      nnew=nact**2
      nact3=nact**2
c      call dsyev('V','U',nact3,s,nact3,w,work,3*nact3,info)
      call rs(nact3,nact3,s,w,1,vct,work,work(nact3+1),ierr)  !eispack
c      thr=1.d-8                 !perhaps too strict
      thr=1.d-6
      ndep=0
      do i=1,nact3
       if(w(i).lt.0.d0)then
        ndep=ndep+1
        write (lupri,'(a,i3,a,f12.8)')'negative eigenvalue: eps(',i,')='
     $           ,w(i)
        elseif(w(i).lt.thr)then
            ndep=ndep+1
        endif
      enddo
      if (ndep.gt.0) then
      write (lupri,*)'Linear dependencies analysis:3'
      write (lupri,'(a,i5,a)')' There are',ndep,' linear dependencies'
      endif
      nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
      call nwadd(nwall,(nact3*nnew)*8)
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,vct,w,nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
c      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      call rs(nact3,nnew,coef,w,1,vct,work,work(nact3+1),ierr)  !eispack
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,vct,nact3,0
     $     .d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact3,nnew)
      call nwdel(nwall,(nact3*nnew)*8)
      deallocate(caux)
      deallocate(vct)
      call e2p2cont(tint,e2,psi2,coef,s,w,daat,fmpb,nact,ncore,
     *              norb,nnew,zgel)
      e2pcc(4)=e2
      psipcc(4)=psi2
      e2pc(1)=e2pc(1)+e2
      psipc(1)=psipc(1)+psi2
      nact2=nact*(nact+1)/2
      call nwdel(nwall,(nact**4+4*nact**2)*8)
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      call nwadd(nwall,(nact**4+4*nact**2)*8)
      allocate(s(1:nact2,1:nact2))
      allocate(coef(1:nact2,1:nact2))
      allocate(vct(1:nact2,1:nact2))
      allocate(w(1:nact2))
      allocate(work(1:3*nact2))
      call fillcm2p(s,daat,nact)
      call fillcm2p(coef,gkp2aa,nact)
      nact3=nact2
      nnew=nact2
c      call dsyev('V','U',nact3,s,nact3,w,work,3*nact3,info)
      call rs(nact3,nact3,s,w,1,vct,work,work(nact3+1),ierr)  !eispack
c      thr=1.d-8                 !perhaps too strict
      thr=1.d-6
      ndep=0
      do i=1,nact3
       if(w(i).lt.0.d0)then
         ndep=ndep+1
        write (lupri,'(a,i3,a,f12.8)')'negative eigenvalue: eps(',i,')='
     $           ,w(i)
       elseif(w(i).lt.thr)then
         ndep=ndep+1
       endif
      enddo
      if (ndep.gt.0) then
      write (lupri,*)'Linear dependencies analysis:4'
      write (lupri,'(a,i5,a)')' There are',ndep,' linear dependencies'
      endif
      nnew=nact3-ndep
c--   costruzione della matrice di trasformazione
      call nwadd(nwall,(nact3*nnew)*8)
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,vct,w,nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
c      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      call rs(nact3,nnew,coef,w,1,vct,work,work(nact3+1),ierr)  !eispack
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,vct,nact3,0
     $     .d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact3,nnew)
      call nwdel(nwall,(nact3*nnew)*8)
      deallocate(caux)
      deallocate(vct)
      call e2p2contp(tint,e2,psi2,coef,s,w,daat,fmpb,nact,ncore,
     *               norb,nnew,zgel)
      e2pcc(4)=e2pcc(4)+e2
      psipcc(4)=psi2+psipcc(4)
      e2pc(1)=e2pc(1)+e2
      psipc(1)=psipc(1)+psi2
      call nwdel(nwall,(nact**4+4*nact**2)*8)
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      call nwdel(nwall,ncore*nact*nvirt*8)
      deallocate(tint)
      return
      end
cc----------------------------------------------------------------
      subroutine vm2mod(daa,gkm2aa,norb,ncore,nact,in2d,id2n,itsym,its,
     *          wrk,lwrk,fmpb,zgel,lupri,nwall)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension in2d(*),id2n(*),itsym(*),its(8,8),needtp(6)
c      dimension da(nact,nact),gkm1a(nact,nact)
      COMMON /NEVPT/ e(1),e2pc(1),e2sc(1),psisc(1),psipc(1),
     *               e2pcc(8),psipcc(8),e2scc(8),psiscc(8),
     *               ntotdet,ntotconf,ntotcap,zro
      dimension WRK(LWRK)
      dimension fmpb(*),zgel(*)
      dimension daa(nact,nact,nact,nact)
      dimension gkm2aa(nact,nact,nact,nact)
      allocatable h2m(:,:),tint(:,:,:,:)
      allocatable coef(:,:),s(:,:),w(:),work(:),caux(:,:),vct(:,:)
      
      ncoat=ncore+nact
      nvirt=norb-ncoat

      do i=1,6
       needtp(i)=0
      enddo

c------ two-electron part with mulliken distr (ir,ia|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
      needtp(5)=1
      L_tint = nvirt*nact*nvirt*nact
      call nwadd(nwall,norb*norb*8)
      allocate(h2m(norb,norb))
      call nwadd(nwall,L_tint*8)
      allocate(tint(nvirt,nact,nvirt,nact))
      call dzero(tint,L_tint)
  900 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1000
      ic=id2n(icd)
      id=id2n(idd)
      if(ic.gt.id) then
      ir=ic
      ia=id
      else
      ir=id
      ia=ic
      endif
      irc=ir-ncoat
      iac=ia-ncore
      do is=1,nvirt
      isd=in2d(is+ncoat)
       do ib=1,nact
        ibd=in2d(ib+ncore)
        tint(irc,iac,is,ib)=h2m(isd,ibd)
       enddo
      enddo
      goto 900
 1000 continue
      call nwdel(nwall,norb*norb*8)
      deallocate(h2m)
C --------- Strongly Contracted -----------------

      do ir=ncore+nact+1,norb ! do 1800
      irc=ir-ncoat
      irsy=itsym(ir)

      do is=ir,norb ! do 1700
      isc=is-ncoat
      issy=itsym(is)
      irssy=its(irsy,issy)

c
c     Si inizia il ciclo sugli attivi 
c
      
      dn=0.d0
      dk=0.d0

      do ic=ncore+1,ncore+nact
         ica=ic-ncore
         do id=ncore+1,ncore+nact
            ida=id-ncore
            icsy=itsym(ic)
            idsy=itsym(id)
            icdsy=its(icsy,idsy)
            isimm=its(icdsy,irssy)
            if (isimm.ne.1) goto 1121
            do icp=ncore+1,ncore+nact
               icpa=icp-ncore
               do idp=ncore+1,ncore+nact
                  idpa=idp-ncore
                  icpsy=itsym(icp)
                  idpsy=itsym(idp)
                  icpdpsy=its(icpsy,idpsy)
                  isimm=its(icpdpsy,irssy)
                  if (isimm.ne.1) goto 1120
                  dum=daa(icpa,idpa,ica,ida)
                  dumk=gkm2aa(icpa,idpa,ica,ida)
                  dint1=tint(isc,icpa,irc,idpa)
                  dint2=tint(isc,ica,irc,ida)
                cint=dint1*dint2
                if (ir.eq.is) then
                cint=cint*0.25d0
                dum=dum+daa(icpa,idpa,ida,ica)
                dumk=dumk+gkm2aa(idpa,icpa,ica,ida)
                endif
                  dn=dn+dum*cint
                  dk=dk+dumk*cint 
 1120          continue
               enddo
            enddo
 1121    continue
         enddo
      enddo
      if(abs(dn).lt.1.d-15)goto 1700
      if(dn.lt.0.d0)then
         write(lupri,*)'Warning vm2mod: dn=',dn
         dn=0.0d0
      endif
      if(dn.gt.0.d0)then
         epsimk=-dk/dn                     
         den=epsimk-fmpb(ir)-fmpb(is)
         co=dn/den
         con=co/den
         e2sc(1)=e2sc(1)+co
         psisc(1)=psisc(1)+con
         e2scc(5)=e2scc(5)+co
         psiscc(5)=psiscc(5)+con
      endif
      
 1700 continue
      enddo
 1800 continue
      enddo

c---Partially contracted NEV-PT-----------------------------------
      tini=second()
      call nwadd(nwall,(nact**4+4*nact**2)*8)
      allocate(s(1:nact**2,1:nact**2))
      allocate(coef(1:nact**2,1:nact**2))
      allocate(vct(1:nact**2,1:nact**2))
      allocate(w(nact**2))
      allocate(work(3*nact**2))
      call fillcm2(s,daa,nact)
      call fillcm2(coef,gkm2aa,nact)
      nnew=nact**2
      nact3=nact**2
c      call dsyev('V','U',nact3,s,nact3,w,work,3*nact3,info)
      call rs(nact3,nact3,s,w,1,vct,work,work(nact3+1),ierr)  !eispack
c      thr=1.d-8                 !perhaps too strict
      thr=1.d-6
      ndep=0
      do i=1,nact3
       if(w(i).lt.0.d0)then
         ndep=ndep+1
         write(lupri,'(a,i3,a,f12.8)')'negative eigenvalue: eps(',i,')='
     $        ,w(i)
       elseif(w(i).lt.thr)then
          ndep=ndep+1
       endif
      enddo
      if (ndep.gt.0) then
      write (lupri,*)'Linear dependencies analysis:5'
      write (lupri,'(a,i5,a)')' There are',ndep,' linear dependencies'
      endif
      nnew=nact3-ndep
      if(nnew.le.0)then
         deallocate(s)
         deallocate(coef)
         deallocate(vct)
         deallocate(w)
         deallocate(work)
         goto 123
      endif
c--   costruzione della matrice di trasformazione
      call nwadd(nwall,(nact3*nnew)*8)
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,vct,w,nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
c      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      call rs(nact3,nnew,coef,w,1,vct,work,work(nact3+1),ierr)  !eispack
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,vct,nact3,0
     $     .d0,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact3,nnew)
      call nwdel(nwall,(nact3*nnew)*8)
      deallocate(caux)
      deallocate(vct)
      call e2m2cont(tint,e2,psi2,coef,s,w,daa,fmpb,nact,ncore,
     *              norb,nnew,zgel)
      e2pcc(5)=e2
      psipcc(5)=psi2
      e2pc(1)=e2pc(1)+e2
      psipc(1)=psipc(1)+psi2
      call nwdel(nwall,(nact**4+4*nact**2)*8)
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
 123  nact2=nact*(nact+1)/2
      call nwadd(nwall,(nact**4+4*nact**2)*8)
      allocate(s(1:nact2,1:nact2))
      allocate(coef(1:nact2,1:nact2))
      allocate(vct(1:nact2,1:nact2))
      allocate(w(1:nact2))
      allocate(work(1:3*nact2))
      call fillcm2p(s,daa,nact)
      call fillcm2p(coef,gkm2aa,nact)
      nact3=nact2
      nnew=nact2
c      call dsyev('V','U',nact3,s,nact3,w,work,3*nact3,info)
      call rs(nact3,nact3,s,w,1,vct,work,work(nact3+1),ierr)  !eispack
c      thr=1.d-8                 !perhaps too strict
      thr=1.d-6
      ndep=0
      do i=1,nact3
       if(w(i).lt.0.d0)then
         ndep=ndep+1
         write(lupri,'(a,i3,a,f12.8)')'negative eigenvalue: eps(',i,')='
     $       ,w(i)
       elseif(w(i).lt.thr)then
          ndep=ndep+1
       endif
      enddo
      if (ndep.gt.0) then
      write (lupri,*)'Linear dependencies analysis:6'
      write (lupri,'(a,i5,a)')' There are',ndep,' linear dependencies'
      endif
      nnew=nact3-ndep
      if(nnew.le.0)then
         deallocate(s)
         deallocate(coef)
         deallocate(vct)
         deallocate(w)
         deallocate(work)
         return
      endif
c--   costruzione della matrice di trasformazione
      call nwadd(nwall,(nact3*nnew)*8)
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,vct,w,nact3,nnew,thr)
c---  trasformazione della matrice K
c---  moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---  moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
c--   diagonalizzo la matrice coef(hamiltoniana trasformata)
c      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      call rs(nact3,nnew,coef,w,1,vct,work,work(nact3+1),ierr)  !eispack
c--   back transformation
c--   moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c--   e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,vct,nact3,0.d0
     $     ,s,nact3)
c--   metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact3,nnew)
      call nwdel(nwall,(nact3*nnew)*8)
      deallocate(caux)
      deallocate(vct)
      call e2m2contp(tint,e2,psi2,coef,s,w,daa,fmpb,nact,ncore,
     *               norb,nnew,zgel)
      e2pcc(5)=e2+e2pcc(5)
      psipcc(5)=psipcc(5)+psi2
      e2pc(1)=e2pc(1)+e2
      psipc(1)=psipc(1)+psi2
      call nwdel(nwall,(nact**4+4*nact**2)*8)
      deallocate(s)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      call nwdel(nwall,nact*nact*nvirt*8)
      deallocate(tint)
      return
      end
c----------------------------------------------------------------
      subroutine fillcm2(s,daa,nact)
#include <implicit.h>
      DIMENSION s(nact**2,nact**2),daa(nact,nact,nact,nact)
c---caso r#s
      ij=0
      do i=1,nact
         do j=1,nact
            ij=ij+1
            kl=0
            do k=1,nact
               do l=1,nact
                  kl=kl+1
                  s(ij,kl)=daa(i,j,k,l)
               enddo
            enddo
         enddo
      enddo
      return
      end
c-------------------------------------------------------------
      subroutine fillcm2p(s,daa,nact)
#include <implicit.h>
      DIMENSION s(nact*(nact+1)/2,nact*(nact+1)/2),
     &          daa(nact,nact,nact,nact)
c---caso r=s
      ij=0
      do i=1,nact
         do j=i,nact
            ij=ij+1
            kl=0
            do k=1,nact
               do l=k,nact
                  kl=kl+1
                  s(ij,kl)=daa(i,j,k,l)+daa(i,j,l,k)
               enddo
            enddo
         enddo
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine fillc0p(s,coef,ga,gd,daa,da,nact)
#include <implicit.h>
      DIMENSION s(2*nact**2,2*nact**2),daa(nact,nact,nact,nact)
      DIMENSION coef(2*nact**2,2*nact**2),da(nact,nact)
      DIMENSION ga(nact,nact,nact,nact),gd(nact,nact,nact,nact)
      integer a,b,ap,bp,ab,abp
      nact2=nact**2
c---parte nord-ovest
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            abp=0
            do ap=1,nact
               do bp=1,nact
                  abp=abp+1
                  coef(ab,abp)=ga(b,a,ap,bp)
                  s(ab,abp)=2.d0*ee2(bp,ap,a,b,daa,da,nact)
               enddo
            enddo
         enddo
      enddo
c---parte nord-est
      ab=0
      do a=1,nact
         do b=1,nact
            ab=ab+1
            abp=nact2
            do ap=1,nact
               do bp=1,nact
                  abp=abp+1
                  coef(ab,abp)=-0.5d0*ga(b,a,ap,bp)
                  s(ab,abp)=-ee2(bp,ap,a,b,daa,da,nact)
               enddo
            enddo
         enddo
      enddo
c---parte sud-ovest
      ab=nact2
      do a=1,nact
         do b=1,nact
            ab=ab+1
            abp=0
            do ap=1,nact
               do bp=1,nact
                  abp=abp+1
                  coef(ab,abp)=-0.5d0*ga(b,a,ap,bp)
                  s(ab,abp)=-ee2(bp,ap,a,b,daa,da,nact)
               enddo
            enddo
         enddo
      enddo
c---parte sud-est
      ab=nact2
      do a=1,nact
         do b=1,nact
            ab=ab+1
            abp=nact2
            do ap=1,nact
               do bp=1,nact
                  abp=abp+1
                  coef(ab,abp)=gd(b,a,ap,bp)
                  s(ab,abp)=-daa(bp,a,b,ap)
                  if(a.eq.ap)s(ab,abp)=s(ab,abp)+2.d0*da(bp,b)
               enddo
            enddo
         enddo
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine fillcm1p(s,coef,amat,ro3,daa,da,nact,nact3)
#include <implicit.h>
      DIMENSION s(nact3,nact3),daa(nact,nact,nact,nact)
      DIMENSION coef(nact3,nact3),da(nact,nact)
      DIMENSION ro3(nact,nact,nact,nact,nact,nact)
      DIMENSION amat(nact,nact,nact,nact,nact,nact)
      integer a,b,c,ap,bp,cp,abc,abcp
      abc=0
      do a=1,nact
         do b=1,nact
            do c=1,nact
               abc=abc+1
               abcp=0
               do ap=1,nact
                  do bp=1,nact
                     do cp=1,nact
                        abcp=abcp+1
                        coef(abc,abcp)=amat(ap,bp,cp,a,b,c)
                        s(abc,abcp)=eee2(cp,ap,bp,b,a,c,ro3,daa,da,nact)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine fillcp1p(s,coef,amat,ro3,daa,da,nact,nact3)
#include <implicit.h>
      DIMENSION s(nact3,nact3),daa(nact,nact,nact,nact)
      DIMENSION coef(nact3,nact3),da(nact,nact)
      DIMENSION ro3(nact,nact,nact,nact,nact,nact)
      DIMENSION amat(nact,nact,nact,nact,nact,nact)
      integer a,b,c,ap,bp,cp,abc,abcp
      abc=0
      do a=1,nact
         do b=1,nact
            do c=1,nact
               abc=abc+1
               abcp=0
               do ap=1,nact
                  do bp=1,nact
                     do cp=1,nact
                        abcp=abcp+1
                        coef(abc,abcp)=amat(ap,bp,cp,a,b,c)
                        s(abc,abcp)=eee2t(cp,ap,bp,b,a,c,ro3,daa,da,nact
     $                       )
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      return
      end
c---------------------------------------------------------------
      subroutine v0pmod(da,daa,gk0pa,gk0pd,gk0pf,heffd,norb,ncore,
     *         nact,in2d,id2n,itsym,its,wrk,lwrk,fmpb,zgel,lupri,nwall)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension in2d(*),id2n(*),itsym(*),its(8,8),needtp(6)
      COMMON /NEVPT/ e(1),e2pc(1),e2sc(1),psisc(1),psipc(1),
     *               e2pcc(8),psipcc(8),e2scc(8),psiscc(8),
     *               ntotdet,ntotconf,ntotcap,zro
      dimension WRK(LWRK)
      dimension fmpb(*),zgel(*)
      allocatable h2m(:,:),tintc(:,:,:,:),tinta(:,:,:,:)
      dimension heffd(*)
      dimension gk0pa(nact,nact,nact,nact),gk0pd(nact,nact,nact,nact),
     * gk0pf(nact,nact)
      dimension daa(nact,nact,nact,nact)
      dimension da(nact,nact)
      allocatable s(:,:),coef(:,:),w(:),work(:),caux(:,:),sp(:,:),spp(:
     $     ,:),vct(:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      
      ncoat=ncore+nact
      nvirt=norb-ncoat

      do i=1,6
       needtp(i)=0
      enddo
      tini=second()

c------ two-electron part with mulliken distr (ir,ia|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
      needtp(5)=1
      L_tinta = nvirt*nact*ncore*nact
      call nwadd(nwall,norb*norb*8)
      allocate(h2m(norb,norb))
      call nwadd(nwall,L_tinta*8)
      allocate(tinta(nvirt,nact,ncore,nact))
      CALL DZERO(tinta,L_tinta)
  900 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1000
      ic=id2n(icd)
      id=id2n(idd)
      if(ic.gt.id) then
      ir=ic
      ia=id
      else
      ir=id
      ia=ic
      endif
      irc=ir-ncoat
      iac=ia-ncore
      do ii=1,ncore
      iid=in2d(ii)
       do ib=1,nact
        ibd=in2d(ib+ncore)
        tinta(irc,iac,ii,ib)=h2m(iid,ibd)
c       write (lupri,*) 'stoa',irc,iac,ii,ib,h2m(iid,ibd)
       enddo
      enddo
      goto 900
 1000 continue
      do i=1,6
       needtp(i)=0
      enddo
c------ two-electron part with mulliken distr (ir,ii|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
      needtp(4)=1
      L_tintc = nvirt*ncore*nact*nact
      call nwadd(nwall,L_tintc*8)
      allocate(tintc(nvirt,ncore,nact,nact))
      call dzero(tintc,L_tintc)
  910 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1010
      ic=id2n(icd)
      id=id2n(idd)
      if(ic.gt.id) then
      ir=ic
      ii=id
      else
      ir=id
      ii=ic
      endif
      irc=ir-ncoat
      do ia=1,nact
      iad=in2d(ia+ncore)
       do ib=1,nact
        ibd=in2d(ib+ncore)
        tintc(irc,ii,ia,ib)=h2m(iad,ibd)
       enddo
      enddo
      goto 910
 1010 continue
      call nwdel(nwall,norb*norb*8)
      deallocate(h2m)
C --------- Strongly Contracted -----------------

      do ii=1,ncore ! do 1800
      if (zgel(ii)) goto 1800
      iis=itsym(ii)

      do ir=ncore+nact+1,norb ! do 1700
      irc=ir-ncoat
      irs=itsym(ir)
      iirs=its(iis,irs)

c
c     Si inizia il ciclo sugli attivi 
c
      
      dn=0.d0
      dkoo=0.d0

      do ia=ncore+1,ncore+nact
         iaa=ia-ncore
         ias=itsym(ia)
         do ib=ncore+1,ncore+nact
            iba=ib-ncore
            ibs=itsym(ib)
            iabs=its(ias,ibs)
            if (its(iabs,iirs).ne.1) goto 1121
            do iap=ncore+1,ncore+nact
               iapa=iap-ncore
               iaps=itsym(iap)
               do ibp=ncore+1,ncore+nact
                  ibpa=ibp-ncore
                  ibps=itsym(ibp)
                  iabps=its(iaps,ibps)
                  if (its(iabps,iirs).ne.1) goto 1120


                  djp=tintc(irc,ii,iapa,ibpa)
                  dj=tintc(irc,ii,iaa,iba)
                  dkp=tinta(irc,ibpa,ii,iapa)
                  dk=tinta(irc,iba,ii,iaa)

                  dint1=2.d0*djp*dj-djp*dk-dkp*dj
                  dint2=dkp*dk

                  dum1=ee2(ibpa,iapa,iaa,iba,daa,da,nact)
                  dum2=ee2tr(ibpa,iba,iapa,iaa,daa,da,nact)
                  if(iaa.eq.iba) dum2=dum2+da(ibpa,iapa)

                  dn=dn+dint1*dum1+dint2*dum2

                  dumk1=gk0pa(ibpa,iapa,iaa,iba)
                  dumk2=gk0pd(ibpa,iapa,iaa,iba)
                  dint1=djp*dj-0.5d0*djp*dk-0.5d0*dj*dkp
                  dint2=dkp*dk
                  dkoo=dkoo+dint1*dumk1+dint2*dumk2

 1120          continue
               enddo ! ibp
            enddo ! iap
 1121    continue
         enddo ! ib
      enddo ! ia

      dummy=0.d0
      dummyk=0.d0
      do ia=ncore+1,ncore+nact
         iaa=ia-ncore
         do ib=ncore+1,ncore+nact
            iba=ib-ncore
            dint=2.d0*tintc(irc,ii,iaa,iba)-tinta(irc,iba,ii,iaa)
            dum=da(iba,iaa)
            dummy=dummy+dum*dint
            dint=tintc(irc,ii,iaa,iba)-0.5d0*tinta(irc,iba,ii,iaa)
            dumk=gk0pf(iaa,iba)
            dummyk=dummyk+dumk*dint
         enddo
      enddo

      iir=indice(ir,ii)
      dn=dn+2.d0*heffd(iir)*dummy+2.d0*heffd(iir)**2
      dkoo=dkoo+heffd(iir)*dummyk
      
      if(abs(dn).lt.1.d-15) goto 1700
      if(dn.lt.0.d0)then
         write(lupri,*)'Warning v0pmod: dn=',dn
         dn=0.d0
      endif
      if(dn.gt.0.d0)then
         epsimk=dkoo/dn
         den=fmpb(ii)-fmpb(ir)-epsimk
         co=dn/den
         con=co/den
         e2sc(1)=e2sc(1)+co
         psisc(1)=psisc(1)+con
         e2scc(8)=e2scc(8)+co
         psiscc(8)=psiscc(8)+con
      endif
      
 1700 continue
      enddo
 1800 continue
      enddo
      tsc=second()
      write (lupri,'(a,f8.2,a)')' Elapsed time for SC',tsc-tini,' sec.'
c---Partially contracted NEV-PT      
      nact2=2*nact**2
      call nwadd(nwall,(2*nact2**2+4*nact2)*8)
      allocate(s(1:nact2,1:nact2))
      allocate(coef(1:nact2,1:nact2))
      allocate(vct(1:nact2,1:nact2))
      allocate(w(nact2))
      allocate(work(3*nact2))
      call fillc0p(s,coef,gk0pa,gk0pd,daa,da,nact)
      call rs(nact2,nact2,s,w,1,vct,work,work(nact2+1),info)  !eispack
c      call dsyev('V','U',nact2,s,nact2,w,work,3*nact2,info)
      if(info.ne.0)then
         write(lupri,*)'something wrong in diagonalization: sub v0pmod'
         CALL QUIT('There is a problem with diagonalization in v0pmod')
      endif
      ndep=0
c      thr=1.d-8
      thr=1.d-6
      do i=1,nact2
       if(w(i).lt.0.d0)then
        ndep=ndep+1
c        write(lupri,'(a,i2,a,f12.8)')' Negative eigenvalue: eps(',i,')='
c     $         ,w(ioffw+i-1)
       elseif(w(i).lt.thr)then
          ndep=ndep+1
       endif
      enddo
      if (ndep.gt.0) then
      write (lupri,*)'Linear dependencies analysis:7'
      write (lupri,'(a,i5,a)')' There are',ndep,' linear dependencies'
      endif
      nnew=nact2-ndep
c--costruzione della matrice di trasformazione
      call nwadd(nwall,(nact2*nnew)*8)
      allocate (caux(1:nact2,1:nnew))
      call caux0p(caux,vct,w,nact2,nnew,thr)
c---trasformazione della matrice K
c---moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact2,nnew,1.d0,coef,nact2,caux,nact2,0.d0,s
     $     ,nact2)
c---moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact2,1.d0,caux,nact2,s,nact2,0.d0
     $     ,coef,nact2)
c--diagonalizzo la matrice coef(hamiltoniana trasformata)
c      call dsyev('V','U',nnew,coef,nact2,w,work,3*nact2,info)
      call rs(nact2,nnew,coef,w,1,vct,work,work(nact2+1),info)  !eispack
      if(info.ne.0)then
         write(lupri,*)'Something wrong in diagonalization: sub v0pmod'
         CALL QUIT('There is a problem with diagonalization in v0pmod')
      endif
      call flshfo(lupri)
c--back transformation
c--moltiplico CAUX(nact2,nnew) per coef(nnew,nnew)
c-- e metto in S(nact2,nnew)
      call dgemm('N','N',nact2,nnew,nnew,1.d0,caux,nact2,vct,nact2,0.d0
     $     ,s,nact2)
c--metto in coef i coefficienti trasformati (for clarity sake)
      call copia(coef,s,nact2,nnew)
c--inizio il calcolo perturbativo
      call nwdel(nwall,(nact2**2)*8)
      deallocate(s)
      call nwdel(nwall,(nact2*nnew)*8)
      deallocate(caux)
      deallocate(vct)
      call nwadd(nwall,3*(nact**2*nnew)*8)
      allocate(s(1:nact**2,1:nnew))
      allocate(sp(1:nact**2,1:nnew))
      allocate(spp(1:nact**2,1:nnew))
      call e20pcont(e2,psi2,heffd,coef,s,sp,spp,w,daa,da,fmpb,nact,ncore
     $     ,norb,nnew,zro,zgel,tintc,tinta)
      e2pcc(8)=e2
      psipcc(8)=psi2
      e2pc(1)=e2pc(1)+e2
      psipc(1)=psipc(1)+psi2
      call nwdel(nwall,3*(nact**2*nnew)*8)
      deallocate(s)
      deallocate(sp)
      deallocate(spp)
      call nwadd(nwall,(nact2**2+4*nact2)*8)
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      call nwdel(nwall,nvirt*nact*ncore*nact*8)
      deallocate(tinta)
      call nwdel(nwall,nvirt*nact*ncore*nact*8)
      deallocate(tintc)
      tpc=second()
      write (lupri,'(a,f8.2,a)')' Elapsed time for PC',tpc-tsc,' sec.'
      return
      end
c----------------------------------------------------------------
      subroutine caux0p(caux,s,w,nact2,nnew,thr)
#include <implicit.h>
      dimension caux(nact2,nnew),s(nact2,nact2),w(*)
      jnew=0
      do j=1,nact2
         if(w(j).ge.thr)then
            jnew=jnew+1
            val=sqrt(w(j))
            do i=1,nact2
               caux(i,jnew)=s(i,j)/val
            enddo
         endif
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine copia(coef,s,nact2,nnew)
#include <implicit.h>
      dimension coef(nact2,nact2),s(nact2,nact2)
      do mu=1,nnew
         do i=1,nact2
            coef(i,mu)=s(i,mu)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------
      subroutine e20pcont(e2,psi2,heffd,coef,s,sp,spp,w,daa,da,fmpb,nact
     $     ,ncore,norb,nnew,zro,zgel,tintc,tinta)
#include <implicit.h>
      logical*1 zro,zgel
c      allocatable rhoir(:,:),vaux(:)
      dimension coef(2*nact**2,2*nact**2),s(nact**2,nnew),
     $     sp(nact**2,nnew),spp(nact**2,nnew),w(*),fmpb(*),
     $     daa(nact,nact,nact,nact),da(nact,nact),heffd(*),zgel(*)
      dimension tintc(norb-ncore-nact,ncore,nact,nact)
      dimension tinta(norb-ncore-nact,nact,ncore,nact)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      nn=nact**2
      nvirt=norb-ncore-nact
      ncoat=ncore+nact
cro      if(zro)then
cro         open(27,file='RHOFILE',form='UNFORMATTED',status='UNKNOWN')
cro         rewind 27
cro         allocate(rhoir(1:ncore,1:nvirt))
cro         allocate(vaux(1:nnew))
cro         call dzero(rhoir,ncore*nvirt)
cro         call dzero(vaux,nnew)
cro      endif
      call dzero(s,nn*nnew)
      call dzero(sp,nn*nnew)
      call dzero(spp,nn*nnew)
c--costruiamo le matrici ausiliarie S(ab,mu),S'(ab,mu) e S''(ab,mu)
c--e il vettore vaux(mu) se si vuole la rho perturbata al primo ordine
      iab=0
      do ia=1,nact
         do ib=1,nact
            iab=iab+1
            iabp=0
            do iap=1,nact
               do ibp=1,nact
                  iabp=iabp+1
                  val=ee2(ibp,iap,ia,ib,daa,da,nact)
                  val2=-daa(ibp,ia,ib,iap)
                  if(iap.eq.ia)val2=val2+2.d0*da(ibp,ib)
                  do mu=1,nnew
                     s(iab,mu)=s(iab,mu)+val*coef(iabp,mu)
                     sp(iab,mu)=sp(iab,mu)+val*coef(iabp+nn,mu)
                     spp(iab,mu)=spp(iab,mu)+val2*coef(iabp+nn,mu)
cro                     if(zro.and.iabp.eq.1)vaux(mu)=vaux(mu)+da(ia,ib)*
cro     $                    (2.d0*coef(iab,mu)-coef(iab+nn,mu))
                  enddo
               enddo
            enddo
         enddo
      enddo
      e2=0.d0
      psi2=0.d0
      do i=1,ncore              !loop sui core
         if (zgel(i)) goto 1
         do ir=ncore+nact+1,norb !loop sui virtuali
         irc=ir-ncoat
            iir=indice(ir,i)
            hir=heffd(iir)
            do mu=1,nnew
               rimu=0.d0
               rimup=0.d0
               iab=0
               do ia=1,nact
                  do ib=1,nact
                     iab=iab+1
                     raib=tintc(irc,i,ia,ib)
                     rabi=tinta(irc,ib,i,ia)
                     rimu=rimu+(2.d0*raib-rabi)*s(iab,mu)+2.d0*da(ib,ia)
     $                    *hir*coef(iab,mu)
                     rimup=rimup-raib*sp(iab,mu)+rabi*spp(iab,mu)-hir
     $                    *da(ib,ia)*coef(iab+nn,mu)
                  enddo
               enddo
               co=rimu+rimup
               deno=w(mu)+fmpb(ir)-fmpb(i)
               cirmu=co/deno
               cont=co*cirmu
               e2=e2-cont
               psi2=psi2+cirmu**2
c               if(zro)then
c                  ivirt=ir-ncore-nact
c                  rhoir(i,ivirt)=rhoir(i,ivirt)-cirmu*vaux(mu)
c               endif
            enddo
c            print*,'i,ivirt,rho',i,ivirt,rhoir(i,ivirt)
         enddo
 1    continue
      enddo
cro      if(zro)write(27)((rhoir(i,ivirt),i=1,ncore),ivirt=1,nvirt)
cro      if(zro)deallocate(rhoir)
      return
      end
c----------------------------------------------------------------
      subroutine e21pcont(e2,psi2,heffd,coef,s,sp,w,ro3,daa,da,fmpb,nact
     $     ,ncore,norb,nnew,zro,zgel,tint)
#include <implicit.h>
      logical*1 zro,zgel
      dimension coef(nact**3,nact**3),s(nnew,nact**3),
     $     sp(nnew,nact**3),w(*),fmpb(*),ro3(nact,nact,nact,nact,nact
     $     ,nact),daa(nact,nact,nact,nact),da(nact,nact),heffd(*),zgel(
     $     *)
      dimension tint(ncore,nact,nact,nact)
      DIMENSION tarray(2)
      logical*1 zcond
c      allocatable rhoid(:,:),caux(:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      t3=second()
      nact3=nact**3
c     ro      if(zro)then
c     ro         allocate(rhoid(1:ncore,1:nact))
c     ro         allocate(caux(1:nact,1:nnew))
c     ro         call dzero(rhoid,nact*ncore)
c     ro         call dzero(caux,nact*nnew)
c     ro      endif
      call dzero(s,nact3*nnew)
      call dzero(sp,nact3*nnew)
      iabc=0
      do ia=1,nact
       do ib=1,nact
        do ic=1,nact
         zcond=ib.eq.1.and.ic.eq.1
         iabc=iabc+1
         iabcp=0
         do iap=1,nact
          do ibp=1,nact
           do icp=1,nact
            iabcp=iabcp+1
            val=eee2t(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact)
            if(zcond)val2=ee2tr(icp,iap,ibp,ia,daa,da,nact)
c     do mu=1,nnew
c     s(mu,iabc)=s(mu,iabc)+val*coef(mu,iabcp)
c     if(zcond)sp(mu,ia)=sp(mu,ia)+val2*coef(mu
c     $                          ,iabcp)
c     enddo
            call daxpy(nnew,val,coef(1,iabcp),1,s(1,iabc),1)
            if(zcond)call daxpy(nnew,val2,coef(1,iabcp),1 ,sp(1,ia),1)
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo
c     ro      if(zro)then
c     ro         do id=1,nact
c     ro            do mu=1,nnew
c     ro               iabc=0
c     ro               do ia=1,nact
c     ro                  do ib=1,nact
c     ro                     do ic=1,nact
c     ro                        iabc=iabc+1
c     ro                        val=ee2tl(id,ib,ia,ic,daa,da,nact)
c     ro                        caux(id,mu)=caux(id,mu)+coef(mu,iabc)
C     *val
c     ro                     enddo
c     ro                  enddo
c     ro               enddo
c     ro            enddo
c     ro         enddo
c     ro      endif
c-------------------------------------------------------
      t4=second()
      e2=0.d0
      psi2=0.d0
      do i=1,ncore              !loop sugli orbitali di core
       if (zgel(i)) goto 1 
       do mu=1,nnew
        aimu=0.d0
        iabc=0
        do ia=1,nact
         ira=indice(i,ia+ncore)
         hia=heffd(ira)
         aimu=aimu+hia*sp(mu,ia)
         do ib=1,nact
          do ic=1,nact
           iabc=iabc+1
           baic=tint(i,ib,ia,ic)
           aimu=aimu+baic*s(mu,iabc)
          enddo
         enddo
        enddo
        deno=w(mu)-fmpb(i)
        co=aimu/deno
        cont=aimu*co
        e2=e2-cont
        psi2=psi2+co**2
c       if(zro)then
c       do id=1,nact
c        rhoid(i,id)=rhoid(i,id)-co*caux(id,mu)
c       enddo
c       endif
       enddo
 1    continue
      enddo
      t5=second()
c     if(zro)write(27)((rhoid(i,id),i=1,ncore),id=1,nact)
c     if(zro)deallocate(rhoid)
      return
      end
c------------------------------------------------------------------
      subroutine e2m1pcont(e2,psi2,heffd,coef,s,sp,w,ro3,daa,da,fmpb
     $     ,nact,ncore,norb,nnew,zro,zgel,tint)
#include <implicit.h>
      logical*1 zro,zgel
      dimension coef(nact**3,nact**3),s(nnew,nact**3),
     $     sp(nnew,nact**3),w(*),fmpb(*),ro3(nact,nact,nact,nact,nact
     $     ,nact),daa(nact,nact,nact,nact),da(nact,nact),heffd(*),zgel(
     $     *)
      dimension tint(norb-nact-ncore,nact,nact,nact)
c      allocatable rhodr(:,:),caux(:,:)
      DIMENSION tarray(2)
      logical*1 zcond
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      ncoat=nact+ncore
      t3=second()
      nact3=nact**3
      nvirt=norb-ncore-nact
cro      if(zro)then
cro         allocate(rhodr(1:nact,1:nvirt))
cro         allocate(caux(1:nact,1:nnew))
cro         call dzero(rhodr,nact*nvirt)
cro         call dzero(caux,nact*nnew)
cro      endif
      call dzero(s,nact3*nnew)
      call dzero(sp,nact3*nnew)
      iabc=0
      do ia=1,nact
         do ib=1,nact
            do ic=1,nact
               zcond=ib.eq.1.and.ic.eq.1
               iabc=iabc+1
               iabcp=0
               do iap=1,nact
                  do ibp=1,nact
                     do icp=1,nact
                        iabcp=iabcp+1
                        val=eee2(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact)
                        if(zcond)val2=ee2(icp,iap,ibp,ia,daa,da,nact)
c$$$                        do mu=1,nnew
c$$$                           s(mu,iabc)=s(mu,iabc)+val*coef(mu,iabcp)
c$$$                           if(zcond)sp(mu,ia)=sp(mu,ia)+val2*coef(mu
c$$$     $                          ,iabcp)
c$$$                        enddo
                        call daxpy(nnew,val,coef(1,iabcp),1,s(1,iabc),1)
                        if(zcond)call daxpy(nnew,val2,coef(1,iabcp),1
     $                       ,sp(1,ia),1)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
cro      if(zro)then
cro         do id=1,nact
cro            do mu=1,nnew
cro               iabc=0
cro               do ia=1,nact
cro                  do ib=1,nact
cro                     do ic=1,nact
cro                        iabc=iabc+1
cro                        val=ee2(id,ib,ia,ic,daa,da,nact)
cro                        caux(id,mu)=caux(id,mu)+coef(mu,iabc)*val
cro                     enddo
cro                  enddo
cro               enddo
cro            enddo
cro         enddo
cro      endif
c-------------------------------------------------------
      t4=second()
      e2=0.d0
      psi2=0.d0
      do ir=ncore+nact+1,norb   !loop sui virtuali
      irc=ir-ncoat
         do mu=1,nnew
            rmu=0.d0
            iabc=0
            do ia=1,nact
               ira=indice(ir,ia+ncore)
               hra=heffd(ira)
               rmu=rmu+hra*sp(mu,ia)
               do ib=1,nact
                  do ic=1,nact
                     iabc=iabc+1
                     rabc=tint(irc,ib,ia,ic)
                     rmu=rmu+rabc*s(mu,iabc)
                  enddo
               enddo
            enddo
            deno=w(mu)+fmpb(ir)
            co=rmu/deno
            cont=rmu*co
            e2=e2-cont
            psi2=psi2+co**2
c            if(zro)then
c               ivirt=ir-nact-ncore
c               do id=1,nact
c                  rhodr(id,ivirt)=rhodr(id,ivirt)-co*caux(id,mu)
c               enddo
c            endif
         enddo
c         print*,'ir,rho',ir,(rhodr(id,ivirt),id=1,nact)
      enddo
      t5=second()
c      if(zro)write(27)((rhodr(id,ivirt),id=1,nact),ivirt=1,nvirt)
c      if(zro)deallocate(rhodr)
c      write(lupri,'(a,f8.2,a)')'time for pert. calcn. ',t5-t4,' sec.'
      return
      end
c----------------------------------------------------------------
      subroutine vpupE(ro3,amat,bmat,cmat,dmat,daa,da,heffd,norb,
     * ncore,nact,in2d,id2n,itsym,its,wrk,lwrk,fmpb,zgel,lupri,nwall)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension in2d(*),id2n(*),itsym(*),its(8,8),needtp(6)
      COMMON /NEVPT/ e(1),e2pc(1),e2sc(1),psisc(1),psipc(1),
     *               e2pcc(8),psipcc(8),e2scc(8),psiscc(8),
     *               ntotdet,ntotconf,ntotcap,zro
      dimension WRK(LWRK)
      dimension fmpb(*),zgel(*)
      allocatable h2m(:,:),tint(:,:,:,:)
      dimension heffd(*)
      dimension ro3(nact,nact,nact,nact,nact,nact),
     *          amat(nact,nact,nact,nact,nact,nact),
     *          bmat(nact,nact,nact,nact),
     *          cmat(nact,nact,nact,nact),
     *          daa(nact,nact,nact,nact),
     *          dmat(nact,nact),
     *          da(nact,nact)
      allocatable s(:,:),coef(:,:),sp(:,:),caux(:,:),w(:),work(:),
     $     vct(:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      
      tini=second()

      ncoat=ncore+nact
      nvirt=norb-ncoat

      do i=1,6
       needtp(i)=0
      enddo

c------ two-electron part with mulliken distr (ia,ii|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
      needtp(2)=1
      L_tint = ncore*nact*nact*nact
      call nwadd(nwall,norb*norb*8)
      allocate(h2m(norb,norb))
      call nwadd(nwall,L_tint*8)
      allocate(tint(ncore,nact,nact,nact))
      call dzero(tint,L_tint)
  900 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1000
      ic=id2n(icd)
      id=id2n(idd)
      if(ic.gt.id) then
      ia=ic
      ii=id
      else
      ia=id
      ii=ic
      endif
      iac=ia-ncore
      do ic=1,nact 
      icd=in2d(ic+ncore)
       do ib=1,nact
        ibd=in2d(ib+ncore)
        tint(ii,iac,ib,ic)=h2m(ibd,icd)
       enddo
      enddo
      goto 900
 1000 continue
      call nwdel(nwall,norb*norb*8)
      deallocate(h2m)

C --------- Strongly Contracted -----------------

      do ii=1,ncore ! do 800
      if (zgel(ii)) goto 800

      dumm=0.d0
      dumk=0.d0
c      dumk1=0.d0
c      dumk2=0.d0
c      dumk3=0.d0

      do ia=1,nact
       iac=ia+ncore
       iia=indice(ii,iac)
       do iap=1,nact
        iapc=iap+ncore
        iiap=indice(ii,iapc)
        rodum=-da(iap,ia)
        if (iap.eq.ia) rodum=rodum+2.d0
        dumm=dumm+heffd(iiap)*rodum*heffd(iia)
        dumk=dumk+heffd(iiap)*dmat(iap,ia)*heffd(iia)
c        dumk1=dumk1+heffd(iiap)*dmat(iap,ia)*heffd(iia)
        do ibp=1,nact
         ibpc=ibp+ncore
         do icp=1,nact
          icpc=icp+ncore
          dint=tint(ii,ibp,iap,icp)
          dumm=dumm+2.d0*dint*ee2tr(icp,iap,ibp,ia,daa,da,nact)
     *             *heffd(iia)
          dumk=dumk+dint*bmat(iap,ibp,icp,ia)*heffd(iia)
          dumk=dumk+dint*cmat(ia,iap,ibp,icp)*heffd(iia)
c          dumk2=dumk2+dint*bmat(iap,ibp,icp,ia)*heffd(iia)
c          dumk2=dumk2+dint*cmat(ia,iap,ibp,icp)*heffd(iia)
          do ib=1,nact
           ibc=ib+ncore
           do ic=1,nact
            icc=ic+ncore
            dint2=tint(ii,ib,ia,ic)
            dumm=dumm+dint*eee2t(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact)
     *           *dint2
            dumk=dumk+dint*amat(iap,ibp,icp,ia,ib,ic)*dint2
c            dumk3=dumk3+dint*amat(iap,ibp,icp,ia,ib,ic)*dint2
            enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      dn=dumm
      dk=dumk
      if(abs(dn).lt.1.d-15) goto 800
      if(dn.lt.0.d0)then
         write(lupri,*)'Warning v1pupe: dn=',dn
         dn=0.0d0
      endif
      if(dn.gt.0.d0)then
         epsimk=dk/dn                     
         den=fmpb(ii)-epsimk
         co=dn/den
         con=co/den
         e2sc(1)=e2sc(1)+co
         psisc(1)=psisc(1)+con
         e2scc(6)=e2scc(6)+co
         psiscc(6)=psiscc(6)+con
      endif

  800 continue
      enddo
      tsc=second()
      write (lupri,'(a,f8.2,a)')' Elapsed time for SC',tsc-tini,' sec.'
c--Partially contracted NEV-PT class (+1')
      nact3=nact**3
      call nwadd(nwall,(2*nact3**2+4*nact3*8))
      allocate(s(1:nact3,1:nact3))
      allocate(coef(1:nact3,1:nact3))
      allocate(vct(1:nact3,1:nact3))
      allocate(w(nact3))
      allocate(work(3*nact3))
      call fillcp1p(s,coef,amat,ro3,daa,da,nact,nact3)
      call rs(nact3,nact3,s,w,1,vct,work,work(nact3+1),info)  !eispack
c      call dsyev('V','U',nact3,s,nact3,w,work,3*nact3,info)
      if(info.ne.0)then
         write(lupri,*)'something wrong in diagonalization: sub vpupe'
         CALL QUIT('There is a problem with diagonalization in vpupe')
      endif
      t=second()
      write(lupri,'(a,f8.2,a)')' Time for S diag. of S',
     *                       t-tsc,' sec.'
c      thr=1.d-8  !perhaps too strict
      thr=1.d-6  !perhaps too strict
c      thr=1.d-5  !perhaps too large
      ndep=0
      do i=1,nact3
       if(w(i).lt.0.d0)then
        ndep=ndep+1
c        write(lupri,'(a,i3,a,f12.8)')'negative eigenvalue: eps(',i,')='
c     $      ,w(ioffw+i-1)
       elseif(w(i).lt.thr)then
        ndep=ndep+1
       endif
      enddo
      if (ndep.gt.0) then
      write (lupri,*)'Linear dependencies analysis:8'
      write (lupri,'(a,i5,a)')' There are',ndep,' linear dependencies'
      endif
      nnew=nact3-ndep
c--costruzione della matrice di trasformazione
      call nwadd(nwall,(nact3*nnew*8))
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,vct,w,nact3,nnew,thr)
c---trasformazione della matrice K
c---moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
      t2=second()
      write(lupri,'(a,f8.2,a)')' Time for matrix transf.',t2-t,' sec.'
c--diagonalizzo la matrice coef(hamiltoniana trasformata)
c      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      call rs(nact3,nnew,coef,w,1,vct,work,work(nact3+1),info)  !eispack
      if(info.ne.0)then
         write(lupri,*)'Something wrong in diagonalization: sub vpupe'
         CALL QUIT('There is a problem with diagonalization in vpupe')
      endif
      t3=second()
      write(lupri,'(a,f8.2,a)')' Time for H diag.',t3-t2,' sec.'
c--back transformation
c--moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c-- e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,vct,nact3,0.d0
     $     ,s,nact3)
c--metto in coef i coefficienti trasformati (for clarity sake)
      do mu=1,nnew
         do i=1,nact3
            coef(mu,i)=s(i,mu)
         enddo
      enddo
c--inizio il calcolo perturbativo
c--costruiamo le matrici ausiliarie S(abc,mu) e S'(a,mu)
      call nwdel(nwall,(nact3**2))
      deallocate(s)
      call nwdel(nwall,(nact3*nnew*8))
      deallocate(caux)
      deallocate(vct)
      call nwadd(nwall,2*(nact3*nnew*8))
      allocate(s(1:nnew,1:nact3))
      allocate(sp(1:nnew,1:nact3))
      call e21pcont(e2,psi2,heffd,coef,s,sp,w,ro3,daa,da,fmpb,nact,ncore
     $     ,norb,nnew,zro,zgel,tint)
      e2pcc(6)=e2
      psipcc(6)=psi2
      e2pc(1)=e2pc(1)+e2
      psipc(1)=psipc(1)+psi2
      call nwdel(nwall,2*(nact3*nnew*8))
      deallocate(s)
      deallocate(sp)
      call nwdel(nwall,(nact3**2+4*nact3*8))
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      call nwdel(nwall,ncore*nact*nact*nact*8)
      deallocate(tint)
      tpc=second()
      write (lupri,'(a,f8.2,a)')' Elapsed time for PC',tpc-tsc,' sec.'
      return
      end
c----------------------------------------------------------------
      subroutine vmupE(ro3,amat,bmat,cmat,dmat,daa,da,heffd,norb,ncore,
     *          nact,in2d,id2n,itsym,its,wrk,lwrk,fmpb,zgel,lupri,nwall)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      dimension in2d(*),id2n(*),itsym(*),its(8,8),needtp(6)
      COMMON /NEVPT/ e(1),e2pc(1),e2sc(1),psisc(1),psipc(1),
     *               e2pcc(8),psipcc(8),e2scc(8),psiscc(8),
     *               ntotdet,ntotconf,ntotcap,zro
      dimension WRK(LWRK)
      dimension fmpb(*),zgel(*)
      allocatable h2m(:,:),tint(:,:,:,:)
      dimension ro3(nact,nact,nact,nact,nact,nact),
     *          amat(nact,nact,nact,nact,nact,nact),
     *          bmat(nact,nact,nact,nact),
     *          cmat(nact,nact,nact,nact),
     *          daa(nact,nact,nact,nact),
     *          dmat(nact,nact),
     *          da(nact,nact)
      dimension heffd(*)
      allocatable s(:,:),coef(:,:),w(:),work(:),caux(:,:),sp(:,:),
     $     vct(:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      
      tini=second()

      ncoat=ncore+nact
      nvirt=norb-ncoat

      do i=1,6
       needtp(i)=0
      enddo

c------ two-electron part with mulliken distr (ir,ia|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
      needtp(5)=1
      L_tint = nvirt*nact*nact*nact
      call nwadd(nwall,norb*norb*8)
      allocate(h2m(norb,norb))
      call nwadd(nwall,L_tint*8)
      allocate(tint(nvirt,nact,nact,nact))
      call dzero(tint,L_tint)
  900 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1000
      ic=id2n(icd)
      id=id2n(idd)
      if(ic.gt.id) then
      ir=ic
      ia=id
      else
      ir=id
      ia=ic
      endif
      irc=ir-ncoat
      iac=ia-ncore
      do ic=1,nact 
      icd=in2d(ic+ncore)
       do ib=1,nact
        ibd=in2d(ib+ncore)
        tint(irc,iac,ib,ic)=h2m(ibd,icd)
       enddo
      enddo
      goto 900
 1000 continue
      call nwdel(nwall,norb*norb*8)
      deallocate(h2m)

C --------- Strongly Contracted -----------------

      do ir=ncore+nact+1,norb ! do 800
      irc=ir-ncoat
      irsy=itsym(ir)

c      write (lupri,*) 
c      write (lupri,*) 'Eccitazione con gli indici',' -> ',ir


      dumm=0.d0
c      dumm1=0.d0
c      dumm2=0.d0
c      dumm3=0.d0
      dumk=0.d0
c      dumk1=0.d0
c      dumk2=0.d0
c      dumk3=0.d0

      do ia=1,nact
       iac=ia+ncore
       ira=indice(ir,iac)
       do iap=1,nact
        iapc=iap+ncore
        irap=indice(ir,iapc)
        dumm=dumm+heffd(irap)*da(iap,ia)*heffd(ira)
c        dumm1=dumm1+heffd(irap)*da(iap,ia)*heffd(ira)
        dumk=dumk+heffd(irap)*dmat(iap,ia)*heffd(ira)
c        dumk1=dumk1+heffd(irap)*dmat(iap,ia)*heffd(ira)
        do ibp=1,nact
         ibpc=ibp+ncore
         do icp=1,nact
          icpc=icp+ncore
          dint=tint(irc,ibp,iap,icp)
          dumm=dumm+2.d0*dint*ee2(icp,iap,ibp,ia,daa,da,nact)
     *             *heffd(ira)
c          dumm2=dumm2+2.d0*dint*ee2(icp,iap,ibp,ia,daa,da,nact)
c     *             *heffd(ira)
          dumk=dumk+dint*bmat(iap,ibp,icp,ia)*heffd(ira)
          dumk=dumk+dint*cmat(ia,iap,ibp,icp)*heffd(ira)
c          dumk2=dumk2+dint*bmat(iap,ibp,icp,ia)*heffd(ira)
c          dumk2=dumk2+dint*cmat(ia,iap,ibp,icp)*heffd(ira)
          do ib=1,nact
           ibc=ib+ncore
           do ic=1,nact
            icc=ic+ncore
            dint2=tint(irc,ib,ia,ic)
c            dumm3=dumm3+dint*eee2(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact)
c     *           *dint2
            dumm=dumm+dint*eee2(icp,iap,ibp,ib,ia,ic,ro3,daa,da,nact)
     *           *dint2
            dumk=dumk+dint*amat(iap,ibp,icp,ia,ib,ic)*dint2
c            dumk3=dumk3+dint*amat(iap,ibp,icp,ia,ib,ic)*dint2
            enddo
          enddo
         enddo
        enddo
       enddo
      enddo
      dn=dumm
      dk=dumk
      if(abs(dn).lt.1.d-15) goto 800
      if(dn.lt.0.d0)then
         write(lupri,*)'Warning vmupe: dn=',dn
         dn=0.0d0
      endif
      if(dn.ne.0.d0)then
         epsimk=-dk/dn                     
         den=epsimk-fmpb(ir)
         co=dn/den
         con=co/den
         e2sc(1)=e2sc(1)+co
         psisc(1)=psisc(1)+con
         e2scc(7)=e2scc(7)+co
         psiscc(7)=psiscc(7)+con
      endif

  800 continue
      enddo
      tsc=second()
      write (lupri,'(a,f8.2,a)')' Elapsed time for SC',tsc-tini,' sec.'
c--Partially contracted NEV-PT class (-1')
      nact3=nact**3
      call nwadd(nwall,(2*nact3**2+4*nact3*8))
      allocate(s(1:nact3,1:nact3))
      allocate(coef(1:nact3,1:nact3))
      allocate(vct(1:nact3,1:nact3))
      allocate(w(nact3))
      allocate(work(3*nact3))
      call fillcm1p(s,coef,amat,ro3,daa,da,nact,nact3)
c      call matout(nact3,nact3,s,nact3)   !renzo debug
      call rs(nact3,nact3,s,w,1,vct,work,work(nact3+1),info)  !eispack
c      call dsyev('V','U',nact3,s(i_s),nact3,w(i_w),work(i_work),3*nact3
c     $,info)  !lapack
      if(info.ne.0)then
         write(lupri,*)'something wrong in diagonalization: sub vmupe'
         CALL QUIT('There is a problem with diagonalization in vmupe')
      endif
      t=second()
      write(lupri,'(a,f8.2,a)')' Time for S diag.',t-tsc,' sec.'
c      thr=1.d-8  !perhaps too strict
      thr=1.d-6
c      thr=1.d-5  !perhaps too large
      ndep=0
      do i=1,nact3
c      write(lupri,'(f12.8)')w(i)    !renzo debug
       if(w(i).lt.0.d0)then
         ndep=ndep+1
c        write(lupri,'(a,i3,a,f12.8)')' Negative eigenvalue: eps(',i,')='
c     $        ,w(ioffw+i-1)
       elseif(w(i).lt.thr)then
         ndep=ndep+1
       endif
      enddo
      if (ndep.gt.0) then
      write (lupri,*)'Linear dependencies analysis:9'
      write (lupri,'(a,i5,a)')' There are',ndep,' linear dependencies'
      endif
      nnew=nact3-ndep
      if(nnew.le.0)then
         deallocate(s)
         deallocate(coef)
         deallocate(vct)
         deallocate(w)
         deallocate(work)
         return
      endif
         
c--costruzione della matrice di trasformazione
      call nwadd(nwall,(nact3*nnew*8))
      allocate (caux(1:nact3,1:nnew))
      call caux0p(caux,vct,w,nact3,nnew,thr)
c---trasformazione della matrice K
c---moltiplico K(coef) per CAUX e metto il risultato in S
      call dsymm('L','U',nact3,nnew,1.d0,coef,nact3,caux,nact3,0.d0,s
     $     ,nact3)
c---moltiplico CAUX(tilde) per S e metto il risultato in coef
      call dgemm('T','N',nnew,nnew,nact3,1.d0,caux,nact3,s,nact3,0.d0
     $     ,coef,nact3)
      t2=second()
      write(lupri,'(a,f8.2,a)')' Time for matrix transf.',t2-t,' sec.'
c--diagonalizzo la matrice coef(hamiltoniana trasformata)
c      call dsyev('V','U',nnew,coef,nact3,w,work,3*nact3,info)
      call rs(nact3,nnew,coef,w,1,vct,work,work(nact3+1),info)  !eispack
      if(info.ne.0)then
         write(lupri,*)'something wrong in diagonalization: sub vmupe'
         CALL QUIT('There is a problem with diagonalization in vmupe')
      endif
      t3=second()
      write(lupri,'(a,f8.2,a)')' Time for H diag. ',t3-t2,' sec.'
c--back transformation
c--moltiplico CAUX(nact3,nnew) per coef(nnew,nnew)
c-- e metto in S(nact3,nnew)
      call dgemm('N','N',nact3,nnew,nnew,1.d0,caux,nact3,vct,nact3,0.d0
     $     ,s,nact3)
c--metto in coef i coefficienti trasformati (for clarity sake)
c      call copia(coef,s,nact3,nnew)
      do mu=1,nnew
         do i=1,nact3
            coef(mu,i)=s(i,mu)
         enddo
      enddo
c--inizio il calcolo perturbativo
c--costruiamo le matrici ausiliarie S(abc,mu) e S'(a,mu)
      call nwdel(nwall,(nact3**2))
      deallocate(s)
      call nwdel(nwall,(nact3*nnew*8))
      deallocate(caux)
      deallocate(vct)
      call nwadd(nwall,2*(nact3*nnew*8))
      allocate(s(1:nnew,1:nact3))
      allocate(sp(1:nnew,1:nact3))
      call e2m1pcont(e2,psi2,heffd,coef,s,sp,w,ro3,daa,da,fmpb,nact
     $     ,ncore,norb,nnew,zro,zgel,tint)
      e2pcc(7)=e2
      psipcc(7)=psi2
      e2pc(1)=e2pc(1)+e2
      psipc(1)=psipc(1)+psi2
      call nwdel(nwall,2*(nact3*nnew*8))
      deallocate(s)
      deallocate(sp)
      call nwdel(nwall,(nact3**2+4*nact3*8))
      deallocate(coef)
      deallocate(w)
      deallocate(work)
      call nwdel(nwall,nvirt*nact*nact*nact*8)
      deallocate(tint)
      tpc=second()
      write (lupri,'(a,f8.2,a)')' Elapsed time for PC',tpc-tsc,' sec.'
      return
      end
c-------------------------------------------------------------------
      function ee2tl(a,b,c,d,taa,dal,nact)
      IMPLICIT REAL*8 (A-H,O-Y), LOGICAL*1(Z)
      integer a,b,c,d
      dimension taa(nact,nact,nact,nact),dal(nact,nact)
      ee2tl=-ee2(b,a,c,d,taa,dal,nact)
      if(a.eq.b)ee2tl=ee2tl+2.d0*dal(c,d)
      return
      end
c----------------------------------------------------
      subroutine getorbeps(WRK,LWRK,onel,heffd,fmpb,da,norb,nact,ncore,
     *           id2n,in2d)
#include <implicit.h>
#include <priunit.h>
      dimension onel(*),fmpb(*),da(nact,*),id2n(*),in2d(*),heffd(*)
      dimension wrk(lwrk),needtp(6)
      allocatable h2m(:,:)
      allocatable fock(:,:)
      indice(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)

Cele-debug
C
C     The Fock matrix is computed in order to verify that
C     the core and virtual orbitals are canonical
C
      allocate (fock(norb,norb))
Cele-debug
       
      do i=1,6
       needtp(i)=0
      enddo
C
C     One-electron contributions to heffd
C

      do i=1,norb
       id=in2d(i)
       do j=i,norb
        jd=in2d(j)
        ij=indice(i,j)
        ijd=indice(id,jd)
        heffd(ij)=onel(ijd)
        fock(i,j)=onel(ijd)
        fock(j,i)=onel(ijd)
       enddo
      enddo


c---The orbital energies are computed following the formula:
c--eps(i)=h(i,i)+Sum(j=1,ncore)(2*J(i,j)-K(i,j))+
c--+Sum(a,b=1,nact)(<ia|ib>-0.5*<ia|bi>)*da(b,a)
c--and analogously for eps(r)

C
C     One-electron contributions to orbital energies
C
c------ Core
      do i=1,ncore
        id=in2d(i)
        indi=indice(id,id)
        fmpb(i)=onel(indi)
      enddo
c------ Virt
      do ir=ncore+nact+1,norb
        ird=in2d(ir)
        indi=indice(ird,ird)
        fmpb(ir)=onel(indi)
      enddo

      allocate (h2m(norb,norb))
c------ two-electron part with mulliken distr (ij|**)
      KFREE=1
      LFREE=LWRK
      IDIST = 0
      needtp(1)=1
  100 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 200
      ic=id2n(icd)
      id=id2n(idd)

c----- (tj,ju) part of heff (t=iar,u=i)
      if (ic.gt.id) then
        ii=ic
        ij=id
      else
        ii=id
        ij=ic
      endif
      iid=in2d(ii)
      ijd=in2d(ij)
      do it=1,norb 
        itd=in2d(it)
        jtj=indice(it,ij)
        heffd(jtj)=heffd(jtj)-h2m(itd,iid)
        if (ii.ne.ij) then
          jti=indice(it,ii)
          heffd(jti)=heffd(jti)-h2m(itd,ijd)
        endif
      enddo

      if (ic.eq.id) then

c----- (tu,jj) part of heff (tu=ri,ra,ai)
      do it=1,norb
       itd=in2d(it)
       do iu=it,norb
        iud=in2d(iu)
        itu=indice(it,iu)
        heffd(itu)=heffd(itu)+2.d0*h2m(itd,iud)
       enddo
      enddo

c----- Core eps
        indi=indice(ic,ic)
         do j=1,ncore
          jd=in2d(j)
          fmpb(ic)=fmpb(ic)+2.d0*h2m(jd,jd)
          fock(ic,ic)=fock(ic,ic)+2.d0*h2m(jd,jd)
          if (ic.ne.j) then
           fock(ic,j)=fock(ic,j)-h2m(icd,jd)
          endif
          if (j.eq.ic) then
          fmpb(ic)=fmpb(ic)-h2m(jd,jd)
          fock(ic,ic)=fock(ic,ic)-h2m(jd,jd)
          endif
         enddo
         do ia=1,nact
            do ib=1,nact
             iad=in2d(ia+ncore)
             ibd=in2d(ib+ncore)
             fmpb(ic)=fmpb(ic)+h2m(iad,ibd)*da(ib,ia)
             fock(ic,ic)=fock(ic,ic)+h2m(iad,ibd)*da(ib,ia)
            enddo
         enddo
c----- Virt eps
         do ir=ncore+nact+1,norb
          ird=in2d(ir)
          fmpb(ir)=fmpb(ir)+2.d0*h2m(ird,ird)
          fock(ir,ir)=fock(ir,ir)+2.d0*h2m(ird,ird)
         do is=ir+1,norb
          isd=in2d(is)
          fock(ir,is)=fock(ir,is)+2.d0*h2m(ird,isd)
          fock(is,ir)=fock(is,ir)+2.d0*h2m(ird,isd)
         enddo
         enddo 

      else  ! if (ic.eq.id) then

c----- Core eps
         fmpb(ic)=fmpb(ic)-h2m(idd,icd)
         fmpb(id)=fmpb(id)-h2m(idd,icd)
         do j=1,ncore
          jd=in2d(j)
          fock(ic,id)=fock(ic,id)+2.d0*h2m(jd,jd)
          fock(id,ic)=fock(id,ic)+2.d0*h2m(jd,jd)
          fock(ic,j)=fock(ic,j)-h2m(idd,jd)
          fock(id,j)=fock(id,j)-h2m(icd,jd)
         enddo
         do ia=1,nact
            do ib=1,nact
             iad=in2d(ia+ncore)
             ibd=in2d(ib+ncore)
             fock(ic,id)=fock(ic,id)+h2m(iad,ibd)*da(ib,ia)
             fock(id,ic)=fock(id,ic)+h2m(iad,ibd)*da(ib,ia)
            enddo
         enddo

      endif ! if (ic.eq.id) then

      goto 100
  200 continue
c------ two-electron part with mulliken distr (ia|**)
      KFREE=1
      IDIST = 0
      needtp(1)=0
      needtp(2)=1
  300 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 400
      ic=id2n(icd)
      id=id2n(idd)
      if (ic.gt.id) then
        ii=id
        ia=ic-ncore
      else
        ii=ic
        ia=id-ncore
      endif
c----- Core
      iid=in2d(ii)
      do ib=1,nact
        ibd=in2d(ib+ncore)
        fmpb(ii)=fmpb(ii)-0.5*h2m(iid,ibd)*da(ib,ia)
         do ij=1,ncore
          ijd=in2d(ij)
          fock(ii,ij)=fock(ii,ij)-0.5*h2m(ijd,ibd)*da(ib,ia)
         enddo
      enddo
      goto 300
  400 continue
c------ two-electron part with mulliken distr (ia,ib|**)
      KFREE=1
      IDIST = 0
      needtp(2)=0
      needtp(3)=1
  500 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 600
      ic=id2n(icd)-ncore
      id=id2n(idd)-ncore
c----- Virt
      do ir=nact+ncore+1,norb
        ird=in2d(ir)
        dum=h2m(ird,ird)*da(ic,id)
        if (ic.ne.id) dum=dum*2.d0
        fmpb(ir)=fmpb(ir)+dum
        fock(ir,ir)=fock(ir,ir)+dum                        
        do is=ir+1,norb
         isd=in2d(is)
         dum=h2m(ird,isd)*da(ic,id)
         if (ic.ne.id) dum=dum*2.d0
         fock(ir,is)=fock(ir,is)+dum
         fock(is,ir)=fock(is,ir)+dum
        enddo
      enddo
      goto 500
  600 continue
c------ two-electron part with mulliken distr (ir,ia|**)
      KFREE=1
      IDIST = 0
      needtp(3)=0
      needtp(5)=1
  700 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 800
      ic=id2n(icd)
      id=id2n(idd)
      if (ic.gt.id) then
        ia=id-ncore
        ir=ic
      else
        ia=ic-ncore
        ir=id
      endif
      ird=in2d(ir)
      do ib=1,nact
        ibd=in2d(ib+ncore)
        fmpb(ir)=fmpb(ir)-0.5*h2m(ird,ibd)*da(ib,ia)
          fock(ir,ir)=fock(ir,ir)-0.5*h2m(ird,ibd)*da(ib,ia)
         do is=ir+1,norb
          isd=in2d(is)
          fock(ir,is)=fock(ir,is)-0.5*h2m(isd,ibd)*da(ib,ia)
         enddo
      enddo
c----- (tj,ju) part of heff (t=r,u=a,j=b)
      iad=in2d(ia+ncore)
      do iu=ncore+1,ncore+nact
        iud=in2d(iu)
        iru=indice(ir,iu)
        heffd(iru)=heffd(iru)-h2m(iad,iud)
      enddo
      goto 700
  800 continue
c------ two-electron part with mulliken distr (ir,ii|**)
      KFREE=1
      IDIST = 0
      needtp(5)=0
      needtp(4)=1
  900 CALL NXTH2M(ICD,IDD,H2M,NEEDTP,WRK,KFREE,LFREE,IDIST)
      IF (IDIST .LT. 0) GO TO 1000
      ic=id2n(icd)
      id=id2n(idd)
      if (ic.gt.id) then
        ii=id
        ir=ic
      else
        ii=ic
        ir=id
      endif
      ird=in2d(ir)
      iid=in2d(ii)
      fmpb(ir)=fmpb(ir)-h2m(ird,iid)
      fock(ir,ir)=fock(ir,ir)-h2m(ird,iid)
      do is=ir+1,norb
        isd=in2d(is)
        fock(ir,is)=fock(ir,is)-h2m(isd,iid)
      enddo

c----- (tj,ju) part of heff (t=r,u=a)
      do iu=ncore+1,ncore+nact
        iud=in2d(iu)
        iru=indice(ir,iu)
        heffd(iru)=heffd(iru)-h2m(iid,iud)
      enddo
      goto 900
 1000 continue

C--- Cele Debug
C      write (lupri,*)'  i          fmpb '
C      do i=1,ncore
C         write (lupri,'(i3,f15.8)')i,fmpb(i)
C      enddo
C      do i=ncore+nact+1,norb
C         write (lupri,'(i3,f15.8)')i,fmpb(i)
C      enddo
C      write (lupri,*) 'HEFFD core-core'
C      do i=1,ncore
C       do j=i,ncore
C        ij=indice(i,j)
C        write (lupri,*)i,j,fock(i,j)
C        write (lupri,*)i,j,fock(j,i)
C       enddo
C      enddo
C      write (lupri,*) 'HEFFD core-virt'
C      do i=1,ncore
C       do j=ncore+nact+1,norb
C        ij=indice(i,j)
C        write (lupri,*)i,j,heffd(ij)
C       enddo
C      enddo
C      write (lupri,*) 'HEFFD core-act'
C      do i=1,ncore
C       do j=ncore+1,ncore+nact
C        ij=indice(i,j)
C        write (lupri,*)i,j,heffd(ij)
C       enddo
C      enddo
C      write (lupri,*) 'HEFFD act-virt'
C      do i=ncore+1,ncore+nact
C       do j=ncore+nact+1,norb
C        ij=indice(i,j)
C        write (lupri,*)i,j,heffd(ij)
C       enddo
C      enddo
C      write (lupri,*) 'HEFFD virt-virt'
C      do i=ncore+nact+1,norb
C       do j=i,norb
C        ij=indice(i,j)
C        write (lupri,*)i,j,fock(i,j)
C       enddo
C      enddo
C--- Cele Debug 
      im=0
      jm=0
      dum=0.d0
      do i=1,ncore
      do j=i+1,ncore
      if (abs(fock(i,j)).gt.dum)then
        im=i
        jm=j
        dum=abs(fock(i,j))
      endif
      enddo
      enddo
      write(lupri,'(a,d10.3)') 
     *  ' Maximum value of core F matrix is:',dum
      write(lupri,'(a,2i3)')
     *  ' for indices = ',in2d(im),in2d(jm)
      dum=0.d0
      im=0
      jm=0
      do i=ncore+nact+1,norb
      do j=i+1,norb
      if (abs(fock(i,j)).gt.dum)then
        im=i
        jm=j
        dum=abs(fock(i,j))
      endif
      enddo
      enddo
      write(lupri,'(a,d10.3)') 
     *  ' Maximum value of virt F matrix is:',dum
      write(lupri,'(a,2i3)')
     *  ' for indices = ',in2d(im),in2d(jm)
      deallocate (h2m)
      deallocate (fock)
      return
      end
! -- end of dypc.F --
