module cc_response_tools_module

  use fundamental
  use precision
  use typedef
  use typedeftype
  use dec_typedef_module
  use IntegralInterfaceMOD
  use tensor_interface_module

  ! DEC DEPENDENCIES (within deccc directory)   
  ! *****************************************
  use array4_simple_operations
  use cc_tools_module
  use ccintegrals
  use dec_fragment_utils
  use dec_tools_module

  public get_ccsd_multipliers_simple,cc_jacobian_rhtr,& 
       & ccsd_eigenvalue_solver
  private

contains

   !> \author Patrick Ettenhuber
   !> \Date September 2013
   !> \brief This debug routine calculates the residual of the CCSD left
   !transformations rho1 and rho2, on input, the converged CCSD amplitudes t1f
   !and t2f must be given, furthermore a first guess for the multipliers m1 and
   !m2, the full ao integrals as array4 the t1 transformation matrices xo yo xv
   !and yv being the occupied and virtual particle and hole matrices, the number
   !of occupied orbitals no, the number of virtual orbitals nv and the number of
   !basis functions nb, MyLsItem is required to calculate the Fock matrix in
   !here
   subroutine get_ccsd_multipliers_simple(rho1,rho2,t1f,t2f,m1,m2,fo,xo,yo,xv,yv,no,nv,nb,MyLsItem,JacobianLT,gao_ex)
     implicit none

     type(lsitem), intent(inout) :: MyLsItem
     real(realk),intent(inout) :: rho1(:,:),rho2(:,:,:,:)
     real(realk),intent(inout) :: t1f(:,:),t2f(:,:,:,:),m1(:,:),m2(:,:,:,:),fo(:,:)
     real(realk),intent(in)    :: xo(:,:),yo(:,:),xv(:,:),yv(:,:)
     integer, intent(in)       :: no,nv,nb
     !> Calculate Jacobian left transformation.
     !> This implies that the RHS is not added at the end
     !> If JacobianLT is not present it is effectively set to false
     logical,intent(in),optional :: JacobianLT
     type(array4),intent(inout),optional:: gao_ex
     real(realk), pointer      :: w1(:), w2(:), w3(:), w4(:)
     real(realk), pointer      :: Lovov(:),Lvoov(:), Lovoo(:)
     real(realk), pointer      :: gvovv(:), gvooo(:), govvv(:), gooov(:), govov(:)
     real(realk), pointer      :: goovv(:)
     real(realk), pointer      :: oof(:),ovf(:),vof(:),vvf(:)
     character(tensor_MSG_LEN)    :: msg
     integer                   :: v4,o4,o3v,ov3,o2v2,ov,b2,v2,o2
     type(matrix)              :: iFock, Dens
     integer                   :: i,a,j,b,ctr
     real(realk)               :: norm,nrmt2,nrmm2
     type(array4)              :: gao
     logical :: JacLT

     if(present(gao_ex))then
        gao = gao_ex
     else
        call get_full_eri(mylsitem,nb,gao)
     endif

     if(present(JacobianLT)) then
        JacLT = JacobianLT
     else
        JacLT=.false.
     end if

     b2   = nb*nb
     v2   = nv*nv
     o2   = no*no
     ov   = no*nv
     o2v2 = ov*ov
     ov3  = ov*v2
     o3v  = ov*o2
     o4   = o2*o2
     v4   = v2*v2

     if( DECinfo%PL>2 )then 

        write (msg,*)"t1 n**2 n"
        call print_norm(t1f,int(ov,kind=8),norm,.true.)
        print *,msg,norm,sqrt(norm)

        write (msg,*)"z1 n**2 n"
        call print_norm(m1,int(ov,kind=8),norm,.true.)
        print *,msg,norm,sqrt(norm)

        write (msg,*)"t2 n**2 n pn**2"
        call print_norm(t2f,int(o2v2,kind=8),norm,.true.)
        nrmt2 = 0.0E0_realk
        nrmm2 = 0.0E0_realk
        do j = 1, no
          do b = 1, nv
            do i = 1, no
              do a = 1,nv
                if(a+(i-1)*nv<=b+(j-1)*nv)then
                  nrmt2 = nrmt2 + t2f(a,i,b,j)**2
                  nrmm2 = nrmm2 + m2(a,i,b,j)**2
                endif
              enddo
            enddo
          enddo
        enddo
        print *,msg,norm,sqrt(norm),nrmt2
        
        write (msg,*)"z2 n**2 n pn**2 pn"
        call print_norm(m2,int(o2v2,kind=8),norm,.true.)
        print *,msg,norm,sqrt(norm),nrmm2,sqrt(nrmm2)
     endif

     call mem_alloc(w2,max(max(no,nv)*nb**3,max(max(max(max(o2v2,ov3),v4),o2*v2),o4)))

     call mem_alloc(oof,o2)
     call mem_alloc(ovf,ov)
     call mem_alloc(vof,ov)
     call mem_alloc(vvf,v2)

     call mem_alloc(govov,o2v2)
     call mem_alloc(goovv,o2v2)
     call mem_alloc(Lovov,o2v2)
     call mem_alloc(Lvoov,o2v2)


     call mem_alloc(gvovv,ov3)
     call mem_alloc(govvv,ov3)

     call mem_alloc(Lovoo,o3v)
     call mem_alloc(gvooo,o3v)
     call mem_alloc(gooov,o3v)


     !construct Ls from gs

     !govov
     call array4_read(gao)
     call successive_4ao_mo_trafo(nb,gao%val,xo,no,yv,nv,xo,no,yv,nv,w2)
     call dcopy(o2v2,gao%val,1,govov,1)

     call array_reorder_4d(2.0E0_realk,gao%val,no,nv,no,nv,[1,2,3,4],0.0E0_realk,Lovov)
     call array_reorder_4d(-1.0E0_realk,gao%val,no,nv,no,nv,[1,4,3,2],1.0E0_realk,Lovov)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvoov
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yo,no,xo,no,yv,nv,w2)
     !Lvoov (a i j b) += gvoov (a i j b)
     call array_reorder_4d(2.0E0_realk,gao%val,nv,no,no,nv,[1,2,3,4],0.0E0_realk,Lvoov)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvvoo
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yv,nv,xo,no,yo,no,w2)
     !write (msg,*)"gvvoo"
     !call print_norm(gao%val,int(no**2*nv**2,kind=8),msg)
     !Lvoov (a i j b) += gvvoo (a b i j)
     call array_reorder_4d(-1.0E0_realk,gao%val,nv,nv,no,no,[1,4,3,2],1.0E0_realk,Lvoov)

     call array_reorder_4d(1.0E0_realk,gao%val,nv,nv,no,no,[3,4,1,2],0.0E0_realk,goovv)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvvov
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yv,nv,xo,no,yv,nv,w2)
     call array_reorder_4d(1.0E0_realk,gao%val,nv,nv,no,nv,[3,4,1,2],0.0E0_realk,govvv)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvovv
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yo,no,xv,nv,yv,nv,w2)
     call dcopy(no*nv**3,gao%val,1,gvovv,1)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gooov
     call successive_4ao_mo_trafo(nb,gao%val,xo,no,yo,no,xo,no,yv,nv,w2)

     call dcopy(nv*no**3,gao%val,1,gooov,1)


     call array_reorder_4d(2.0E0_realk,gao%val,no,no,no,nv,[3,4,1,2],0.0E0_realk,Lovoo)
     call array_reorder_4d(-1.0E0_realk,gao%val,no,no,no,nv,[1,4,3,2],1.0E0_realk,Lovoo)

     call array4_dealloc(gao)
     call array4_read(gao)

     !gvooo
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yo,no,xo,no,yo,no,w2)
     call dcopy(nv*no**3,gao%val,1,gvooo,1)

     call array4_dealloc(gao)

     call mem_alloc(w1,max(max(max(max(o2v2,ov3),v4),o2*v2),o4))

     !Transform inactive Fock matrix into the different mo subspaces
     ! -> Foo
     call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,fo,nb,0.0E0_realk,w1,no)
     call dgemm('n','n',no,no,nb,1.0E0_realk,w1,no,yo,nb,0.0E0_realk,oof,no)
     ! -> Fov
     call dgemm('n','n',no,nv,nb,1.0E0_realk,w1,no,yv,nb,0.0E0_realk,ovf,no)
     ! -> Fvo
     call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,fo,nb,0.0E0_realk,w1,nv)
     call dgemm('n','n',nv,no,nb,1.0E0_realk,w1,nv,yo,nb,0.0E0_realk,vof,nv)
     ! -> Fvv
     call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1,nv,yv,nb,0.0E0_realk,vvf,nv)

     call mem_alloc(w3,max(max(max(max(o2v2,ov3),v4),o2*v2),o4))
     call mem_alloc(w4,max(max(ov3,o3v),o2v2))


     !The notation in this routine is according to Halkier et. al. J. chem. phys.,
     !Vol 107. No.3 15 july 1997, and the order of terms is according to the
     !DALTON program  for simple comparison of the terms

     !SINGLES EXPRESSIONS
     !*******************

     !rho f
     !-----
     rho1 = 0.0E0_realk

     ! part1
     !sort amps(dkfj) -> dkjf
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     ! w3 : \sum_{fj} t^{df}_{kj} (d k f j) Lovvv (j f e a)

     call array_reorder_4d(2.0E0_realk,govvv,no,nv,nv,nv,[1,2,3,4],0.0E0_realk,w1)
     call array_reorder_4d(-1.0E0_realk,govvv,no,nv,nv,nv,[1,4,3,2],1.0E0_realk,w1)
     call dgemm('n','n',ov,v2,ov,1.0E0_realk,w2,ov,w1,ov,0.0E0_realk,w3,ov)
     !write(msg,*)"rho f 1"
     !call print_norm(w3,int(ov3,kind=8),msg)
     ! part2

     ! sort t2f (e j f k) -[1,4,2,3]> t2f (e k j f)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w2)
     ! sort govvv(j a d f) -[1,4,2,3]> govvv (j f a d) 
     call array_reorder_4d(1.0E0_realk,govvv,no,nv,nv,nv,[1,4,2,3],0.0E0_realk,w1)
     ! w1 : \sum_{jf} t^{ef}_{jk}(e k j f) govvv_{jadf}(j f a d)
     call dgemm('n','n',ov,v2,ov,1.0E0_realk,w2,ov,w1,ov,0.0E0_realk,w4,ov)
     !write(msg,*)"rho f 2"
     !call print_norm(w4,int(ov3,kind=8),msg)
     ! sort result w4(e k a d) -[4,2,1,3]+> w3 (d k e a)
     call array_reorder_4d(-1.0E0_realk,w4,nv,no,nv,nv,[4,2,1,3],1.0E0_realk,w3)
     !write(msg,*)"rho f 2 - added"
     !call print_norm(w3,int(ov3,kind=8),msg)
     ! part3

     ! sort t2f (d j f k) -[1,4,2,3]> t2f (d k j f)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w2)
     ! \sum_{jf} t^{df}_{jk} (d k j f)(still in w2) govvv(j f e a)
     call dgemm('n','n',ov,v2,ov,-1.0E0_realk,w2,ov,govvv,ov,1.0E0_realk,w3,ov)
     !write(msg,*)"rho f 3"
     !call print_norm(w3,int(ov3,kind=8),msg)

     !\sum_{dke} \zeta^{d e}_{k i}(d k e , i)^T w3(d k e a)
     call dgemm('t','n',no,nv,v2*no,1.0E0_realk,m2,v2*no,w3,v2*no,0.0E0_realk,w2,no)
     !write(msg,*)"rho f1-3(LT21I)"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T

     call mat_transpose(no,nv,1.0E0_realk,w2,1.0E0_realk,rho1)

     !rho c - 1
     !-----
     ! part1
     ! sort \zeta^{df}_{kl} (d k f l) -[3 1 2 4]> w1(f d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! sort t^{de}_{kl} (d k e l) -[1,2,4,3]> w2(d k l e)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     ! w3 : \sum_{dkl} w1 (f d k l) w2 (d k l e)
     call dgemm('n','n',nv,nv,no*no*nv,1.0E0_realk,w1,nv,w2,nv*no*no,0.0E0_realk,w3,nv)
     ! w1 : \sum_{fe} w3 (fe) L_{feia} (f e i a)
     call array_reorder_4d(2.0E0_realk,govvv,no,nv,nv,nv,[3,4,1,2],0.0E0_realk,w1)
     call array_reorder_4d(-1.0E0_realk,govvv,no,nv,nv,nv,[3,2,1,4],1.0E0_realk,w1)
     call dgemv('t',v2,ov,1.0E0_realk,w1,v2,w3,1,0.0E0_realk,w2,1)
     !write(msg,*)"rho c - 1(LT21A)"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,1.0E0_realk,w2,1.0E0_realk,rho1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho1 after (LT21A)"
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho g 1-3
     !-----
     ! part1
     ! sort Lovoo (j f i l) -[2,1,3,4]> (f j i l)
     call array_reorder_4d(1.0E0_realk,Lovoo,no,nv,no,no,[2,1,3,4],0.0E0_realk,w2)
     ! -\sum_{jf} t^{df}_{kj} (d k f j) Lovoo (f j i l)
     call dgemm('n','n',ov,o2,ov,-1.0E0_realk,t2f,ov,w2,ov,0.0E0_realk,w3,ov)
     ! part2
     ! sort gooov (jkif) -[1,4,2,3]> (jfki)
     call array_reorder_4d(1.0E0_realk,gooov,no,no,no,nv,[1,4,2,3],0.0E0_realk,w2)
     ! sort t^{df}_{jl} (d j f l) -[1,4,2,3]> (dljf)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],0.0E0_realk,w1)
     ! \sum_{jf} t^{df}_{jl} gooov_{jkif}(jfki)
     call dgemm('n','n',ov,o2,ov,1.0E0_realk,w1,ov,w2,ov,0.0E0_realk,w4,ov)
     ! sort result w4(d l k i) -[1,3,4,2]+> w3 (d k i l)
     call array_reorder_4d(1.0E0_realk,w4,nv,no,no,no,[1,3,4,2],1.0E0_realk,w3)
     ! part3
     ! govoo from gooov through 3,4,1,2
     call array_reorder_4d(1.0E0_realk,gooov,no,no,no,nv,[3,4,1,2],0.0E0_realk,w2)
     !\sum_{jf} t^{df}_{jk}(d k j f)in w1 govoo_{jfil}(jfil)
     call dgemm('n','n',ov,o2,ov,1.0E0_realk,w1,ov,w2,ov,1.0E0_realk,w3,ov)

     ! sort w3 (d k i l) -[1,2,4,3]> w2 (d k l i)
     call array_reorder_4d(1.0E0_realk,w3,nv,no,no,no,[1,2,4,3],0.0E0_realk,w2)
     ! sort \zeta^{da}_{kl} (d k a l} -[3,1,2,4]> w1(a d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! \sum_{dkl} \zeta^{da}_{kl} w2(dkli)
     call dgemm('n','n',nv,no,o2*nv,1.0E0_realk,w1,nv,w2,o2*nv,0.0E0_realk,w4,nv)
     !write(msg,*)"rho g terms 1-3 (21H)"
     !call print_norm(w4,int(ov,kind=8),msg)
     !rho1 += w2(i a)
     call daxpy(ov,1.0E0_realk,w4,1,rho1,1)


     !rho e
     !-----
     ! part1
     ! \sum_{dke} \zeta^{de}_{ki}(d k e, i)^T g_{dkea} (d k e a)
     call dgemm('t','n',no,nv,v2*no,1.0E0_realk,m2,v2*no,gvovv,v2*no,0.0E0_realk,w2,no)
     !write(msg,*)"rho e - 1 (21DC)" 
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,1.0E0_realk,w2,1.0E0_realk,rho1)

     if( DECinfo%PL > 2) then
        write(msg,*)"singles after loop" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif


     ! f-term
     ! part4
     ! sort t^{de}_{jl} (d j e l) -[1,3,2,4]> (d e j l)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,3,2,4],0.0E0_realk,w1)
     ! sort m^{de}_{ki} (d k e i) -[2,4,1,3]> (k i d e)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[2,4,1,3],0.0E0_realk,w2)
     ! \sum_{de} m2(k i d e)  t^{de}_{jl} (d e j l)
     call dgemm('n','n',o2,o2,v2,1.0E0_realk,w2,o2,w1,v2,0.0E0_realk,w3,o2)
     ! w3 (k i j l) -> w2 (i j l k)
     call array_reorder_4d(1.0E0_realk,w3,no,no,no,no,[2,3,4,1],0.0E0_realk,w2)
    
     ! sort gooov(jkla) -[1,3,2,4]> (j l k a)
     call array_reorder_4d(1.0E0_realk,gooov,no,no,no,nv,[1,3,2,4],0.0E0_realk,w1)
     ! \sum_{jl} w3(i j l k) gooov_{jkla} (j l k a)
     call dgemm('n','n',no,nv,o2*no,1.0E0_realk,w2,no,w1,no*o2,0.0E0_realk,w3,no)
     !write(msg,*)"rho f 4(LT21G)"
     !call print_norm(w3,int(ov,kind=8),msg)
     call mat_transpose(no,nv,1.0E0_realk,w3,1.0E0_realk,rho1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho 4.f(21G):" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho g
     !-----
     ! part 4
     ! sort t^{fe}_{kl} (f k e l) -[2,4,1,3]> (k l f e)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[2,4,1,3],0.0E0_realk,w1)
     ! sort govvv_{iedf} -[3,4,1,2]> gvvov_{dfie} (d f i e) -[2,4,1,3]> (f e d i)
     ! : [4,2,3,1]
     call array_reorder_4d(1.0E0_realk,govvv,no,nv,nv,nv,[4,2,3,1],0.0E0_realk,w2)
     ! \sum_{ef} t^{fe}_{kl}(k l f e)  gvvov_{dfie} (f e d i)
     call dgemm('n','n',o2,ov,v2,1.0E0_realk,w1,o2,w2,v2,0.0E0_realk,w4,o2)
     ! sort result w4(k l d i) -[3,1,2,4]> w3 (d k l i)
     call array_reorder_4d(1.0E0_realk,w4,no,no,nv,no,[3,1,2,4],0.0E0_realk,w3)
     ! sort \zeta^{da}_{kl} (d k a l} -[3,1,2,4]> w1(a d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! \sum_{dkl} \zeta^{da}_{kl} (a d k l) w3(dkli)
     call dgemm('n','n',nv,no,o2*nv,1.0E0_realk,w1,nv,w3,o2*nv,0.0E0_realk,w4,nv)

     !rho e
     !-----
     ! part2
     ! sort: \zeta^{da}_{kl} (d k a l) -[3,1,2,4]> w1 (a, d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! sort: gvooo_{dkil} -[1,2,4,3]> (dkli)
     call array_reorder_4d(1.0E0_realk,gvooo,nv,no,no,no,[1,2,4,3],0.0E0_realk,w3)
     ! \sum_{dke} \zeta^{da}_{kl}(a, dkl) g_{dkil} (d k l i)
     call dgemm('n','n',nv,no,o2*nv,1.0E0_realk,w1,nv,w3,o2*nv,0.0E0_realk,w2,nv)
     !w4 += w2(i a)^T
     call daxpy(ov,1.0E0_realk,w2,1,w4,1)
     !write(msg,*)"4.g+2.e 21BF:"
     !call print_norm(w4,int(ov,kind=8),msg)
     call daxpy(ov,-1.0E0_realk,w4,1,rho1,1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho 4.g +2.e(21BF):" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho a
     !-----
     ! sort \hat{L}_{bkia} (b k i a) -> w1
     call dcopy(o2v2,Lvoov,1,w1,1)
     ! sort u2^{bd}_{kl} (b k d l) -[1,2,4,3]> (b k l d)
     call array_reorder_4d(2.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     call array_reorder_4d(-1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],1.0E0_realk,w2)
     ! w1 : \sum_{dl} u^{bd}_{kl}(b k l d) \hat{L}_{ldia} (l d i a) + \hat{L}_{bkia} (b k i a)
     call dgemm('n','n',ov,ov,ov,1.0E0_realk,w2,ov,Lovov,ov,1.0E0_realk,w1,ov)
     ! w2 : \sum_{bk}  w1(b k , i a)^T \zeta_k^b (b k)
     call dgemv('t',ov,ov,1.0E0_realk,w1,ov,m1,1,0.0E0_realk,w2,1)
     !write(msg,*)"rho a (11A)"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,1.0E0_realk,w2,1.0E0_realk,rho1)
     
     if( DECinfo%PL > 2) then
        write(msg,*)"rho a (11A):" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho b
     !-----
     !part2
     ! w2 : sort t(d j b k) -[3,1,4,2]> t (b d k j)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[3,1,4,2],0.0E0_realk,w2)
     ! w3 : sort \hat{L}_{kbid} (k b i d} -[3,2,4,1]> \hat{L}_{kbid} (i b d k)
     call array_reorder_4d(1.0E0_realk,Lovov,no,nv,no,nv,[3,2,4,1],0.0E0_realk,w3)
     ! w1 : sort  F_{ij} -> w1
     call dcopy(o2,oof,1,w1,1)
     ! w1 : \sum_{bdk} \hat{L}_{kbid}(i b d k) t^{db}_{jk} (b d k j) + \hat{F}_{ij}
     call dgemm('n','n',no,no,no*v2,1.0E0_realk,w3,no,w2,no*v2,1.0E0_realk,w1,no)
     ! w2 : \sum_{j} \zeta_{j}^{a} (a j) w1(ij)^T = w2
     call dgemm('n','t',nv,no,no,1.0E0_realk,m1,nv,w1,no,0.0E0_realk,w2,nv) 
     !write(msg,*)"rho b - 2"
     !call print_norm(w2,int(ov,kind=8),msg)
     call daxpy(ov,-1.0E0_realk,w2,1,rho1,1)
     !part1
     ! w2 : sort t(d l b k) -[3,2,1,4]> t (b l d k)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[3,2,1,4],0.0E0_realk,w2)
     ! w1 : sort  F_{ba} -> w1
     call dcopy(v2,vvf,1,w1,1)
     ! w1 : \sum_{bdk} t^{db}_{lk} (b l d k)\hat{L}_{ldka}(l d k a)  + \hat{F}_{ba}
     call dgemm('n','n',nv,nv,o2*nv,-1.0E0_realk,w2,nv,Lovov,o2*nv,1.0E0_realk,w1,nv)
     ! w2 : \sum_{j} w1(b a)^T \zeta_{j}^{b} (b j) = w2
     call dgemm('t','n',nv,no,nv,1.0E0_realk,w1,nv,m1,nv,0.0E0_realk,w2,nv) 
     !write(msg,*)"rho b - 1"
     !call print_norm(w2,int(ov,kind=8),msg)
     call daxpy(ov,1.0E0_realk,w2,1,rho1,1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho b (11B):"
        call print_norm(rho1,int(ov,kind=8),msg)
     endif
     
     !rho c - part2
     !-------------
     ! w3 : \sum_{dke} t^{de}_{kl} (d k e ,l)^T \zeta^{de}_{kj} (d k e j)
     call dgemm('t','n',no,no,no*v2,1.0E0_realk,t2f,v2*no,m2,v2*no,0.0E0_realk,w3,no)
     ! w1 : \sum_{lj} w3 (lj) L_{l j i a} (l j i a)
     call array_reorder_4d(2.0E0_realk,gooov,no,no,no,nv,[1,2,3,4],0.0E0_realk,w1)
     call array_reorder_4d(-1.0E0_realk,gooov,no,no,no,nv,[3,2,1,4],1.0E0_realk,w1)
     call dgemv('t',o2,ov,1.0E0_realk,w1,o2,w3,1,0.0E0_realk,w2,1)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,-1.0E0_realk,w2,1.0E0_realk,rho1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho c - 2(LT21B)"
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     !rho d
     !-----
     ! part1
     ! sort \zeta^{da}_{kl} (d k a l) -[3 1 2 4]> w1(a d k l)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! sort t^{de}_{kl} (d k e l) -[1,2,4,3]> w2(d k l e)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     ! w3 : \sum_{dkl} w1 (a d k l) w2 (d k l e)
     call dgemm('n','n',nv,nv,o2*nv,1.0E0_realk,w1,nv,w2,nv*o2,0.0E0_realk,w3,nv)
     ! w1 : \sum_{e} w3 (ae) F_{ie} (i e)^T
     call dgemm('n','t',nv,no,nv,1.0E0_realk,w3,nv,ovf,no,0.0E0_realk,w2,nv)
     !write(msg,*)"rho d - 1"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call daxpy(ov,-1.0E0_realk,w2,1,rho1,1)
     ! part2
     ! w3 : \sum_{dke} zeta^{de}_{ki} (d k e ,i)^T t^{de}_{kl} (d k e ,l)
     call dgemm('t','n',no,no,no*v2,1.0E0_realk,m2,v2*no,t2f,v2*no,0.0E0_realk,w3,no)
     ! w1 : \sum_{e} w3 (il) F_{la} (l a)
     call dgemm('n','n',no,nv,no,1.0E0_realk,w3,no,ovf,no,0.0E0_realk,w2,no)
     !write(msg,*)"rho d - 2"
     !call print_norm(w2,int(ov,kind=8),msg)
     !rho1 += w2(i a)^T
     call mat_transpose(no,nv,-1.0E0_realk,w2,1.0E0_realk,rho1)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho d(LT2EFM)"
        call print_norm(rho1,int(ov,kind=8),msg)

        write(msg,*)"singles end" 
        call print_norm(rho1,int(ov,kind=8),msg)
     endif

     call mem_dealloc(w4)
     
   
     !DOUBLES EXPRESSIONS
     !*******************

     !rho B
     !-----
     rho2 = 0.0E0_realk
     !gvvvv
     call array4_read(gao)
     call successive_4ao_mo_trafo(nb,gao%val,xv,nv,yv,nv,xv,nv,yv,nv,w2)
     ! sort gvvvv_{cadb} (cadb) -> (cd ab)
     call array_reorder_4d(1.0E0_realk,gao%val,nv,nv,nv,nv,[1,3,2,4],0.0E0_realk,w1)
     call array4_dealloc(gao)
     ! sort \zeta^{cd}_{ij} (c i d j) -> (i j c d)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[2,4,1,3],0.0E0_realk,w2)
     ! \sum_{cd} \zeta{cd}_{ij} (ij cd) gvvvv_{cadb} (cd ab)
     call dgemm('n','n',o2,v2,v2,1.0E0_realk,w2,o2,w1,v2,0.0E0_realk,w3,o2)
     !write(msg,*)"rho B"
     !call print_norm(w3,int(o2v2,kind=8),msg)
     ! order w3(i j a b) and add to residual rho2 (a i b j)
     call array_reorder_4d(0.5E0_realk,w3,no,no,nv,nv,[3,1,4,2],1.0E0_realk,rho2)

     !rho 2. H part 
     ! \sum_{c} \zeta_{j}^{c} (c,j)^T Lvvov_{cbia} (cbia)
     call array_reorder_4d(2.0E0_realk,govvv,no,nv,nv,nv,[3,4,1,2],0.0E0_realk,w2)
     call array_reorder_4d(-1.0E0_realk,govvv,no,nv,nv,nv,[3,2,1,4],1.0E0_realk,w2)
     call dgemm('t','n',no,v2*no,nv,1.0E0_realk,m1,nv,w2,nv,0.0E0_realk,w1,no)
     !write(msg,*)"rho H - 2"
     !call print_norm(w1,int(o2v2,kind=8),msg)
     ! order w1(j b i a) and add to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w1,no,nv,no,nv,[4,3,2,1],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"H(B-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif

     !rho E
     !----- 
     ! sort t^{cd}_{mn} (c m d n) -[1,3,2,4]> (c d m n)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,3,2,4],0.0E0_realk,w1)
     ! sort govov_{manb} (m a n b) -[1,3,2,4]> (m n a b)
     call array_reorder_4d(1.0E0_realk,govov,no,nv,no,nv,[1,3,2,4],0.0E0_realk,w2)
     ! \sum_{mn} t^{cd}_{mn} (c d m n) govov_{manb} (m n a b)
     call dgemm('n','n',v2,v2,o2,1.0E0_realk,w1,v2,w2,o2,0.0E0_realk,w3,v2)
     ! sort \zeta^{cd}_{ij} (c i d j) -[2,4,1,3]> (i j c d)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[2,4,1,3],0.0E0_realk,w2)
     ! \sum_{cd} \zeta^{cd}_{ij} (i j c d) w3 (c d a b)
     call dgemm('n','n',o2,v2,v2,0.5E0_realk,w2,o2,w3,v2,0.0E0_realk,w1,o2)
     !write(msg,*)"rho E"
     !call print_norm(w1,int(o2v2,kind=8),msg)
     ! order w1(i j  a b) and add to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w1,no,no,nv,nv,[3,1,4,2],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"E(AM-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !rho F
     !-----
     !part1
     ! \sum_{efm} \zeta^{ef}_{mj} (e m f , j)^T t^{ef}_{mn} (e m f , n)
     call dgemm('t','n',no,no,v2*no,1.0E0_realk,m2,v2*no,t2f,v2*no,0.0E0_realk,w1,no)
     ! sort Lovov{ianb} (i a n b) -[2,1,4,3]> (a i b n)
     call array_reorder_4d(1.0E0_realk,Lovov,no,nv,no,nv,[2,1,4,3],0.0E0_realk,w2)
     ! \sum_{n} Lovov_{ianb} (aibn) w1 (j n)^T
     call dgemm('n','t',v2*no,no,no,1.0E0_realk,w2,v2*no,w1,no,0.0E0_realk,w3,v2*no)
     !write(msg,*)"rho F - 1"
     !call print_norm(w3,int(o2v2,kind=8),msg)
     !rho2 += w1(a i b j)
     call daxpy(o2v2,-1.0E0_realk,w3,1,rho2,1)
     !part2
     ! sort \zeta^{ea}_{mn} (e m a n) -[3,1,2,4]> (a e m n)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[3,1,2,4],0.0E0_realk,w1)
     ! sort \t^{ef}_{mn} (e m f n) -[1,2,4,3]> (e m n f)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w2)
     ! \sum_{emn} \zeta^{ea}_{mn}(a e m n) t^{ef}_{mn} (e m n f)
     call dgemm('n','n',nv,nv,o2*nv,1.0E0_realk,w1,nv,w2,o2*nv,0.0E0_realk,w3,nv)
     ! sort Lovov{ifjb} (i f j b) -> (f i b j)
     call array_reorder_4d(1.0E0_realk,Lovov,no,nv,no,nv,[2,1,4,3],0.0E0_realk,w2)
     ! \sum_{f} w3 (a f) Lovov_{ifjb} (f i b j) 
     call dgemm('n','n',nv,o2*nv,nv,1.0E0_realk,w3,nv,w2,nv,0.0E0_realk,w1,nv)
     !write(msg,*)"rho F - 2"
     !call print_norm(w1,int(o2v2,kind=8),msg)
     !rho2 += w1(a i b j)
     call daxpy(o2v2,-1.0E0_realk,w1,1,rho2,1)

     if( DECinfo%PL > 2) then
        write (msg,*)"F(EM-TRM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !rho G
     !-----
     ! part1
     ! copy vv F to w3
     call dcopy(v2,vvf,1,w3,1)
     ! sort t^{fe}_{nm} (f n e m) -[3,2,1,4]> (e n f m)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[3,2,1,4],0.0E0_realk,w1)
     ! vv F - \sum_{fmn} t^{fe}_{nm}(e n f m) Lovov_{nfmb} (n f m b)
     call dgemm('n','n',nv,nv,o2*nv,-1.0E0_realk,w1,nv,Lovov,o2*nv,1.0E0_realk,w3,nv)
     ! sort \zeta^{ae}_{ij} (a i e j) -> (a i j e)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w1)
     ! \sum_{e} \zeta^{ae}_{ij} (a i j e) w3 (e a)
     call dgemm('n','n',o2*nv,nv,nv,1.0E0_realk,w1,o2*nv,w3,nv,0.0E0_realk,w2,o2*nv)
     !write(msg,*)"rho G - 1"
     !call print_norm(w2,int(o2v2,kind=8),msg)
     ! order w1(a i j b) and add to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w2,nv,no,no,nv,[1,2,4,3],1.0E0_realk,rho2)
     !part2
     ! copy oo F to w3
     call dcopy(o2,oof,1,w3,1)
     ! sort t^{fe}_{nm} (f n e m) -[4,3,1,2]> (m e f n)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[4,3,1,2],0.0E0_realk,w2)
     ! sort Lovov_{mejf} (mejf} -[3,1,2,4]> (j m e f)
     call array_reorder_4d(1.0E0_realk,Lovov,no,nv,no,nv,[3,1,2,4],0.0E0_realk,w1)
     ! oo F + \sum_{efm} Lovov_{mejf} (j mef) t^{fe}_{nm}(mef n) 
     call dgemm('n','n',no,no,v2*no,1.0E0_realk,w1,no,w2,v2*no,1.0E0_realk,w3,no)
     ! \sum_{n} \zeta^{ab}_{in} (a i b n) w3 (j n)^T
     call dgemm('n','t',v2*no,no,no,1.0E0_realk,m2,v2*no,w3,no,0.0E0_realk,w2,v2*no)
     !write(msg,*)"rho G - 2"
     !call print_norm(w2,int(o2v2,kind=8),msg)
     ! order w1(a i b j) and add to residual rho2 (a i b j)
     call array_reorder_4d(-1.0E0_realk,w2,nv,no,nv,no,[1,2,3,4],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"G(22E-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !rho A
     !-----
     ! copy goooo[j n i m ]  -[2,4,3,1]> w1[n m i j]
     !goooo
     call array4_read(gao)
     call successive_4ao_mo_trafo(nb,gao%val,xo,no,yo,no,xo,no,yo,no,w2)
     call array_reorder_4d(1.0E0_realk,gao%val,no,no,no,no,[2,4,3,1],0.0E0_realk,w1)
     call array4_dealloc(gao)
     ! sort govov{jfie}(jfie) -[4,2,3,1]> (e f i j)
     call array_reorder_4d(1.0E0_realk,govov,no,nv,no,nv,[4,2,3,1],0.0E0_realk,w2)
     ! sort t^{fe}_{nm} (f n e m) -[2,4,3,1]> (n m e f)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[2,4,3,1],0.0E0_realk,w3)
     !\sum_{ef} t^{fe}_{nm} (n m e f) govov_{jfie} (ef ij) + goooo
     call dgemm('n','n',o2,o2,v2,1.0E0_realk,w3,o2,w2,v2,1.0E0_realk,w1,o2)
     ! sort \zeta^{ab}_{mn} (a m b n) -[1,3,4,2]> (a b n m)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[1,3,4,2],0.0E0_realk,w2)
     !\frac{1}{2} \sum_{mn} w2:\zeta^{ab}_{mn} (a b n m} w1(n m i j)
     call dgemm('n','n',v2,o2,o2,0.5E0_realk,w2,v2,w1,o2,0.0E0_realk,w3,v2)
     !write(msg,*)"rho A"
     !call print_norm(w3,int(o2v2,kind=8),msg)
     ! order and add (a b i j) to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w3,nv,nv,no,no,[1,3,2,4],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"A(22A-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif

     !rho D
     !-----
     ! copy Lvoov in w3 
     call dcopy(o2v2,Lvoov,1,w3,1)
     ! sort u2^{ef}_{mn} (emfn) -> em nf
     call array_reorder_4d(2.0E0_realk,t2f,nv,no,nv,no,[1,2,4,3],0.0E0_realk,w1)
     call array_reorder_4d(-1.0E0_realk,t2f,nv,no,nv,no,[1,4,2,3],1.0E0_realk,w1)
     ! \sum_{fn} u2^{ef}_{mn} (emnf) Lovov_{nfjb} (nfjb) + Lvoov_{emjb} (emjb)
     call dgemm('n','n',ov,ov,ov,1.0E0_realk,w1,ov,Lovov,ov,1.0E0_realk,w3,ov)
     ! \sum_{em} \zeta^{ae}_{im} (a i e m) w3 (e m j b)
     call dgemm('n','n',ov,ov,ov,1.0E0_realk,m2,ov,w3,ov,0.0E0_realk,w2,ov)
     !write(msg,*)"rho D"
     !call print_norm(w2,int(o2v2,kind=8),msg)
     ! order w1(a i j b) and add to residual rho2 (a i b j)
     call array_reorder_4d(1.0E0_realk,w2,nv,no,no,nv,[1,2,4,3],1.0E0_realk,rho2)
     ! order w1(a j i b) and add to residual rho2 (a i b j)
     call array_reorder_4d(-0.5E0_realk,w2,nv,no,no,nv,[1,3,4,2],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"D(22D-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif

     !rho C
     !-----
     ! sort govov(nbif) -[4,1,2,3]> (fnbi)
     call array_reorder_4d(1.0E0_realk,govov,no,nv,no,nv,[4,1,2,3],0.0E0_realk,w1)
     ! sort t^{ef}_{nm} (enfm) -[1,4,3,2]> (emfn)
     call array_reorder_4d(1.0E0_realk,t2f,nv,no,nv,no,[1,4,3,2],0.0E0_realk,w2)
     ! sort goovv(imeb) -[3,2,4,1]> (e m b i)
     call array_reorder_4d(1.0E0_realk,goovv,no,no,nv,nv,[3,2,4,1],0.0E0_realk,w3)
     ! \sum_{fn} t^{ef)_{nm}(em fn) govov_{nbif} (fn bi) + goovv
     call dgemm('n','n',ov,ov,ov,-1.0E0_realk,w2,ov,w1,ov,1.0E0_realk,w3,ov)
     ! sort make anti-u2 analog of \zeta^{ae}_{mj} (a m e j) as 2 (a j e m) + (a m e j)
     call array_reorder_4d(2.0E0_realk,m2,nv,no,nv,no,[1,4,3,2],0.0E0_realk,w2)
     call array_reorder_4d(1.0E0_realk,m2,nv,no,nv,no,[1,2,3,4],1.0E0_realk,w2)
     ! \sum_{em} w2(a j e m) w3(em bi)
     call dgemm('n','n',ov,ov,ov,1.0E0_realk,w2,ov,w3,ov,0.0E0_realk,w1,ov)
     !write(msg,*)"rho C"
     !call print_norm(w1,int(o2v2,kind=8),msg)
     ! order w1(a j b i) and add to residual rho2 (a i b j)
     call array_reorder_4d(-0.5E0_realk,w1,nv,no,nv,no,[1,4,3,2],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write (msg,*)"C(22C-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif

     !rho H
     !-----
     ! 2*\zeta_i^a * F_{jb} - \zeta_j^b F_{ib}
     do j = 1, no
       do b = 1, nv
         do i = 1, no
           do a = 1 ,nv
             rho2(a,i,b,j) = rho2(a,i,b,j) + 2*m1(a,i)*ovf(j+(b-1)*no) - m1(a,j)*ovf(i+(b-1)*no)
           enddo
         enddo
       enddo
     enddo

     if( DECinfo%PL > 2) then
        write (msg,*)"H(A12-TERM)"
        print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !\sum_{k} Lovoo_{jbik} \zeta_{k}^{a} (a k)^T
     call dgemm('n','t',o2*nv,nv,no,1.0E0_realk,Lovoo,o2*nv,m1,nv,0.0E0_realk,w1,o2*nv)
     ! order w1(j b i a) and add to residual rho2 (a i b j)
     call array_reorder_4d(-1.0E0_realk,w1,no,nv,no,nv,[4,3,2,1],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
       write (msg,*)"H(C12-TERM)"
       print *,msg,symmetrized_packed_1norm(rho2,nv,no)
     endif


     !END OF CONTRIBUTIONS
     !********************


     !PERMUTE RHO 2
     call array_reorder_4d(1.0E0_realk,rho2,nv,no,nv,no,[1,2,3,4],0.0E0_realk,w1)
     call array_reorder_4d(1.0E0_realk,w1,nv,no,nv,no,[3,4,1,2],1.0E0_realk,rho2)

     if( DECinfo%PL > 2) then
        write(msg,*)"rho end"
        call print_norm(rho2,int(o2v2,kind=8),msg)
     endif

     !ADD RIGHT HAND SIDES  (not for Jacobian left transformation)
     if(.not. JacLT) then
        call array_reorder_4d(2.0E0_realk,Lovov,no,nv,no,nv,[2,1,4,3],1.0E0_realk,rho2)
        call mat_transpose(no,nv,2.0E0_realk,ovf,1.0E0_realk,rho1)
     end if

     if( DECinfo%PL > 2) then
        write(*,*)"rho + RHS"
        call print_norm(rho1,int(no*nv,kind=8),"R1+RHS")
        call print_norm(rho2,int(o2v2,kind=8),"R2+RHS")
        write(*,*)"---------"
     endif

     call mem_dealloc(oof)
     call mem_dealloc(ovf)
     call mem_dealloc(vof)
     call mem_dealloc(vvf)

     call mem_dealloc(govov)
     call mem_dealloc(goovv)
     call mem_dealloc(Lovov)
     call mem_dealloc(Lvoov)


     call mem_dealloc(gvovv)
     call mem_dealloc(govvv)

     call mem_dealloc(Lovoo)
     call mem_dealloc(gvooo)
     call mem_dealloc(gooov)

     call mem_dealloc(w1)
     call mem_dealloc(w2)
     call mem_dealloc(w3)

     if(.not.present(gao_ex))then
        call array4_free(gao)
     endif

    end subroutine get_ccsd_multipliers_simple

    !This function calculates the norm according to how it should be compared to
    !DALTON, because they construct the symmetrization for each term
    !individually and store the packed residual, furthermore they only print the
    !square of the 2 norm
    function symmetrized_packed_1norm(m,nv,no)result(nrm)
      implicit none
      integer,intent(in)     :: nv,no
      real(realk),intent(in) :: m(nv,no,nv,no)
      integer                :: a,i,b,j
      real(realk)            :: nrm

      nrm = 0.0E0_realk
      do b = 1, nv
        do j = 1, no
          do a = 1,nv
            do i = 1, no
              if(a+(i-1)*nv<=b+(j-1)*nv)then
                nrm = nrm + (m(a,i,b,j) + m(b,j,a,i) )**2
              endif
            enddo
          enddo
        enddo
      enddo

    end function symmetrized_packed_1norm
    subroutine get_ccsd_lhtr_integral_driven(ccmodel,deltafock,omega2,t2,fock,govov,no,nv,&
          ppfock,qqfock,pqfock,qpfock,xo,xv,yo,yv,nb,MyLsItem, omega1,iter,local,rest)
       implicit none

       !> CC model
       integer,intent(in) :: ccmodel
       !> Number of basis functions
       integer,intent(in) :: nb
       !> Number of occupied orbitals
       integer,intent(in) :: no
       !> Number of virtual orbitals
       integer,intent(in) :: nv

       ! derived types needed for the calculation
       !> Long-range correction to Fock matrix
       real(realk),intent(in)    :: deltafock(nb*nb)
       !> the zeroed doubles residual vector on input, on output this contains the
       !full doubles residual
       !real(realk),intent(inout) :: omega2(nv*nv*no*no)
       type(tensor),intent(inout) :: omega2
       !> the current guess amplitudes
       !real(realk),pointer,intent(in) ::t2(:)
       type(tensor),intent(inout) :: t2
       !> on output this contains the occupied-occupied block of the t1-fock matrix
       real(realk),intent(inout) :: ppfock(no*no)
       !> on output this contains the virtual-virtual block of the t1-fock matrix
       real(realk),intent(inout) :: qqfock(nv*nv)
       !> on output this contains the occupied-virtual block of the t1-fock matrix
       real(realk),intent(inout) :: pqfock(no*nv)
       !> on output this contains the virtual-occupied block of the t1-fock matrix
       real(realk),intent(inout) :: qpfock(nv*no)
       !> the ao-fock-matrix
       real(realk),intent(inout) :: fock(nb*nb)
       !> zeroed on input, on output this contains the singles residual
       real(realk),intent(inout) :: omega1(no*nv)
       !> on input this contains the transformation matrices for t1-transformed
       !integrals
       real(realk),pointer :: xo(:),xv(:),yo(:),yv(:)

       !> LS item with information needed for integrals
       type(lsitem), intent(inout) :: MyLsItem

       !govov is passed from the dec ccsd driver and is returned,
       !only it is used in another shape here in the subroutine
       !real(realk), dimension(nvirt*nocc*nvirt*nocc) :: govov
       ! for the master iter is the iteration, mdimg = 0
       ! for the slaves iter contains maximum allowed dim alpha
       ! and mdimg maximum allowed dim gamma
       integer,intent(in) :: iter
       !real(realk),intent(inout) :: govov(nv*no*nv*no)
       type(tensor),intent(inout) :: govov
       logical, intent(in) :: local
       !logical that specifies whether the amplitudes were read
       logical, optional, intent(inout) :: rest

       ! elementary types needed for the calculation
       type(mpi_realk)      :: gvvoo,gvoov,w0,w1,w2,w3,uigcj
       real(realk), pointer :: Had(:), t2_d(:,:,:,:), Gbi(:)
       type(c_ptr) :: Hadc,t2_dc, Gbic
       integer(kind=ls_mpik) :: Hadw,t2_dw,Gbiw,gvvoow,gvoovw
       type(tensor) :: tpl, tmi

       integer(kind=8) :: w0size,w1size,w2size,w3size,neloc

       ! Variables for mpi
       logical :: master,lg_master,parent
       integer :: fintel,nintel,fe,ne
       integer(kind=ls_mpik) :: nnod
       real(realk) :: startt, stopp
       integer(kind=ls_mpik) :: ierr

       integer :: sio4_mode, sio4_dims(4),sio4_tdim(4) 
       type(tensor) :: u2, sio4
       type(tensor) :: gvoova,gvvooa
       !special arrays for scheme=1
       type(tensor) :: t2jabi,u2kcjb
       integer,pointer       :: tasks(:)
       type(c_ptr)           :: tasksc
       integer(kind=ls_mpik) :: tasksw,taskslw
       integer(kind=ls_mpik) :: lg_me,lg_nnod
       integer(kind=8)       :: len81,len82
       integer(kind=4)       :: len41,len42
       integer               :: lenI1,lenI2
       integer,parameter :: inflen=5
       real(realk)       :: inf(inflen)
       logical :: lock_outside

       ! CHECKING and MEASURING variables
       integer(kind=long) :: maxsize64,dummy64
       integer :: myload,nelms,n4
       real(realk) :: tcpu, twall,tcpu1,twall1,tcpu2,twall2, deb1,deb2,ActuallyUsed
       real(realk) :: MemFree,MemFree2,MemFree3,MemFree4
       real(realk) :: tcpu_end,twall_end,time_a, time_c, time_d,time_singles
       real(realk) :: time_doubles,timewall_start,wait_time,max_wait_time,min_wait_time,ave_wait_time
       integer     :: scheme
       integer(kind=8) :: els2add
       logical :: memfound

       ! variables used for BATCH construction and INTEGRAL calculation
       integer :: alphaB,gammaB,dimAlpha,dimGamma
       integer :: dim1,dim2,dim3,K,MinAObatch
       integer :: GammaStart, GammaEnd, AlphaStart, AlphaEnd
       integer :: iorb,nthreads,idx,residual_nr
       type(DecAObatchinfo),pointer :: AOGammabatchinfo(:)
       type(DecAObatchinfo),pointer :: AOAlphabatchinfo(:)
       integer :: iAO,nAObatches,AOGammaStart,AOGammaEnd,AOAlphaStart,AOAlphaEnd,iprint
       logical :: MoTrans, NoSymmetry,SameMol
       Character(80)        :: FilenameCS,FilenamePS
       Character(80),pointer:: BatchfilenamesCS(:,:)
       Character(80),pointer:: BatchfilenamesPS(:,:)
       logical :: FoundInMem,doscreen
       integer, pointer :: orb2batchAlpha(:), batchsizeAlpha(:), batchindexAlpha(:)
       integer, pointer :: orb2batchGamma(:), batchsizeGamma(:), batchindexGamma(:)
       TYPE(DECscreenITEM)    :: DecScreen
       integer, pointer :: batchdimAlpha(:), batchdimGamma(:)     
       type(batchtoorb), pointer :: batch2orbAlpha(:)
       type(batchtoorb), pointer :: batch2orbGamma(:)
       Character            :: INTSPEC(5)
       logical :: fullRHS
       integer :: MaxAllowedDimAlpha,MaxActualDimAlpha,nbatchesAlpha
       integer :: MaxAllowedDimGamma,MaxActualDimGamma,nbatchesGamma

       integer :: a,b,i,j,l,m,n,c,d,fa,fg,la,lg,worksize
       integer :: nb2,nb3,nv2,no2,b2v,o2v,v2o,no3,vs,os
       integer(kind=8) :: nb4,o2v2,no4,buf_size
       integer :: tlen,tred,nor,nvr,goffs,aoffs
       integer :: prev_alphaB,mpi_buf,ccmodel_copy
       logical :: jobtodo,first_round,dynamic_load,restart,print_debug
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !TEST AND DEVELOPMENT VARIABLES!!!!!
       real(realk) :: op_start,op_stop, dt_last_phase
       real(realk) :: time_init,    time_init_work,    time_init_comm, time_cndonly
       real(realk) :: time_intloop, time_intloop_work, time_intloop_comm, time_intloop_idle, time_intloop_int
       real(realk) :: time_intloop_B1work, time_intloop_B1comm, time_intloop_stop
       real(realk) :: time_fock_mat,commtime,time_reduction2,time_Bonly
       real(realk) :: time_Bcnd,time_Bcnd_work,time_Bcnd_comm, time_Bcnd_idle
       real(realk) :: time_cnd,time_cnd_work,time_cnd_comm,time_get_ao_fock, time_get_mo_fock
       real(realk) :: time_Esing,time_Esing_work,time_Esing_comm, time_Esing_idle
       real(realk) :: unlock_time, waiting_time, flushing_time
       real(realk) :: phase_counters_int_dir(nphases)
       integer :: testmode(4)
       integer(kind=long) :: xyz,zyx1,zyx2,mem_allocated,HeapMemoryUsage
       logical :: debug
       character(tensor_MSG_LEN) :: msg
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef VAR_OMP
       integer, external :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
#endif
       character(4) :: def_atype
       integer, parameter :: bs = 1
       !init timing variables
       call time_start_phase(PHASE_WORK, twall = twall)
       call DetectHeapMemoryInit()
       time_init_work       = 0.0E0_realk
       time_init_comm       = 0.0E0_realk
       time_intloop_work    = 0.0E0_realk
       time_intloop_comm    = 0.0E0_realk
       time_intloop_int     = 0.0E0_realk
       time_intloop_B1work  = 0.0E0_realk
       time_intloop_B1comm  = 0.0E0_realk
       time_intloop_idle    = 0.0E0_realk
       time_cnd_work        = 0.0E0_realk
       time_cnd_comm        = 0.0E0_realk
       time_Bcnd_work       = 0.0E0_realk
       time_Bcnd_comm       = 0.0E0_realk
       time_Esing_work      = 0.0E0_realk
       time_Esing_comm      = 0.0E0_realk
       commtime             = 0.0E0_realk
       unlock_time          = 0.0E0_realk 
       waiting_time         = 0.0E0_realk
       flushing_time        = 0.0E0_realk


       ! Set default values for the path throug the routine
       ! **************************************************
       restart                  = .false.
       if(present(rest))restart = rest
       scheme                   = 0
       dynamic_load             = DECinfo%dyn_load
       startt                   = 0.0E0_realk
       stopp                    = 0.0E0_realk
       print_debug              = (DECinfo%PL>3)
       debug                    = .false.

       ! Set some shorthand notations
       ! ****************************
       nb2                      = nb*nb
       nb3                      = nb2*nb
       nb4                      = int((i8*nb3)*nb,kind=long)
       nv2                      = nv*nv
       no2                      = no*no
       no3                      = no2*no
       no4                      = int((i8*no3)*no,kind=long)
       b2v                      = nb2*nv
       o2v                      = no2*nv
       v2o                      = nv2*no
       o2v2                     = int((i8*nv2)*no2,kind=long)
       nor                      = no*(no+1)/2
       nvr                      = nv*(nv+1)/2
       vs                       = t2%tdim(1)
       os                       = t2%tdim(3)
       
       ! Memory info
       ! ***********
       call get_currently_available_memory(MemFree)

       ! Set integral info
       ! *****************
       INTSPEC(1)               = 'R' !R = Regular Basis set on the 1th center 
       INTSPEC(2)               = 'R' !R = Regular Basis set on the 2th center 
       INTSPEC(3)               = 'R' !R = Regular Basis set on the 3th center 
       INTSPEC(4)               = 'R' !R = Regular Basis set on the 4th center 
       INTSPEC(5)               = 'C' !C = Coulomb operator
       IF(DECinfo%useIchor)THEN
          iprint = 0           !print level for Ichor Integral code
          MoTrans = .FALSE.    !Do not transform to MO basis! 
          NoSymmetry = .FALSE. !Use Permutational Symmetry! 
          SameMol = .TRUE.     !Same molecule on all centers of the 4 center 2 electron integral
          !Determine the full number of AO batches - not to be confused with the batches of AOs
          !Required by the MAIN_ICHORERI_DRIVER unless all four dimensions are batched 
          iAO = 1
          call determine_Ichor_nAObatches(mylsitem%setting,iAO,'R',nAObatches,DECinfo%output)
       ELSE
          doscreen = MyLsItem%setting%scheme%cs_screen.OR.MyLsItem%setting%scheme%ps_screen
       ENDIF

       ! Set MPI related info
       ! ********************
       master                   = .true.
       lg_master                = .true.
       parent                   = .true.
       lg_me                    = int(0,kind=ls_mpik)
       lg_nnod                  = 1


       if(master.and.print_debug)then
          call print_norm(xo,int(nb*no,kind=8)," NORM(xo)    :")
          call print_norm(xv,int(nb*nv,kind=8)," NORM(xv)    :")
          call print_norm(yo,int(nb*no,kind=8)," NORM(yo)    :")
          call print_norm(yv,int(nb*nv,kind=8)," NORM(yv)    :")
       endif

       ! Initialize stuff
       IF(DECinfo%useIchor)THEN
          nullify(AOGammabatchinfo)
          nullify(AOalphabatchinfo)    
       ELSE
          nullify(orb2batchAlpha)
          nullify(batchdimAlpha)
          nullify(batchsizeAlpha)
          nullify(batch2orbAlpha)
          nullify(batchindexAlpha)
          nullify(orb2batchGamma)
          nullify(batchdimGamma)
          nullify(batchsizeGamma)
          nullify(batch2orbGamma)
          nullify(batchindexGamma)
       ENDIF
       nullify(Had)
       nullify(Gbi)
!#ifdef VAR_MPI
!       nullify(tasks)
!#endif

       if(master) then
          !==================================================
          !                  Batch construction             !
          !==================================================

          ! Get free memory and determine maximum batch sizes
          ! -------------------------------------------------
          IF(DECinfo%useIchor)THEN
             !Determine the minimum allowed AObatch size MinAObatch
             !In case of pure Helium atoms in cc-pVDZ ((4s,1p) -> [2s,1p]) MinAObatch = 3 (Px,Py,Pz)
             !In case of pure Carbon atoms in cc-pVDZ ((9s,4p,1d) -> [3s,2p,1d]) MinAObatch = 6 (the 2*(Px,Py,Pz))
             !In case of pure Carbon atoms in 6-31G   ((10s,4p) -> [3s,2p]) MinAObatch = 3 (Px,Py,Pz) 
             !'R'  !Specifies that it is the Regular AO basis that should be batched
             iAO = 4 !the center that the batching should occur on (they are all the same in this case)  
             call determine_MinimumAllowedAObatchSize(MyLsItem%setting,iAO,'R',MinAObatch)
          ELSE
             call determine_maxBatchOrbitalsize(DECinfo%output,MyLsItem%setting,MinAObatch,'R')
          ENDIF
          ! INSERT AUTOMATIC BATCH DETERMINATION
          MaxAllowedDimGamma = DECinfo%ccsdGbatch
          MaxAllowedDimAlpha = DECinfo%ccsdAbatch


       endif


       call tensor_zero(omega2)


       ! ************************************************
       ! * Determine batch information for Gamma batch  *
       ! ************************************************
       IF(DECinfo%useIchor)THEN
          iAO = 4 !Gamma is the 4. Center of the 4 center two electron coulomb integral
          !Determine how many batches of AOS based on the MaxAllowedDimGamma, the requested
          !size of the AO batches. iAO is the center that the batching should occur on. 
          !'R'  !Specifies that it is the Regular AO basis that should be batched 
          call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
               & nbatchesGamma,DECinfo%output)
          call mem_alloc(AOGammabatchinfo,nbatchesGamma)
          !Construct the batches of AOS based on the MaxAllowedDimGamma, the requested
          !size of the AO batches - MaxAllowedDimGamma must be unchanged since the call 
          !to determine_Ichor_nbatchesofAOS
          !MaxActualDimGamma is an output parameter indicating How big the biggest batch was, 
          !So MaxActualDimGamma must be less og equal to MaxAllowedDimGamma
          call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimGamma,&
               & nbatchesGamma,AOGammabatchinfo,MaxActualDimGamma,DECinfo%output)
       ELSE
          ! Orbital to batch information
          ! ----------------------------
          call mem_alloc(orb2batchGamma,nb)
          call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimGamma,&
               & nb,MaxActualDimGamma,batchsizeGamma,batchdimGamma,batchindexGamma,&
               &nbatchesGamma,orb2BatchGamma,'R')
       ENDIF
       if(master.and.DECinfo%PL>1)write(DECinfo%output,*) 'BATCH: Number of Gamma batches   = ', nbatchesGamma,&
          & 'with maximum size',MaxActualDimGamma

       IF(.NOT.DECinfo%useIchor)THEN
          ! Translate batchindex to orbital index
          ! -------------------------------------
          call mem_alloc(batch2orbGamma,nbatchesGamma)
          do idx=1,nbatchesGamma
             call mem_alloc(batch2orbGamma(idx)%orbindex,batchdimGamma(idx))
             batch2orbGamma(idx)%orbindex = 0
             batch2orbGamma(idx)%norbindex = 0
          end do
          do iorb=1,nb
             idx = orb2batchGamma(iorb)
             batch2orbGamma(idx)%norbindex = batch2orbGamma(idx)%norbindex+1
             K = batch2orbGamma(idx)%norbindex
             batch2orbGamma(idx)%orbindex(K) = iorb
          end do
       ENDIF
       ! ************************************************
       ! * Determine batch information for Alpha batch  *
       ! ************************************************

       IF(DECinfo%useIchor)THEN
          iAO = 3 !Alpha is the 3. Center of the 4 center two electron coulomb integral
          !Determine how many batches of AOS based on the MaxAllowedDimAlpha, the requested
          !size of the AO batches. iAO is the center that the batching should occur on. 
          !'R'  !Specifies that it is the Regular AO basis that should be batched 
          call determine_Ichor_nbatchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
               & nbatchesAlpha,DECinfo%output)
          call mem_alloc(AOAlphabatchinfo,nbatchesAlpha)
          !Construct the batches of AOS based on the MaxAllowedDimAlpha, the requested
          !size of the AO batches - MaxAllowedDimAlpha must be unchanged since the call 
          !to determine_Ichor_nbatchesofAOS
          !MaxActualDimAlpha is an output parameter indicating How big the biggest batch was, 
          !So MaxActualDimAlpha must be less og equal to MaxAllowedDimAlpha
          call determine_Ichor_batchesofAOS(mylsitem%setting,iAO,'R',MaxAllowedDimAlpha,&
               & nbatchesAlpha,AOAlphabatchinfo,MaxActualDimAlpha,DECinfo%output)
       ELSE
          ! Orbital to batch information
          ! ----------------------------
          call mem_alloc(orb2batchAlpha,nb)
          call build_batchesofAOS(DECinfo%output,mylsitem%setting,MaxAllowedDimAlpha,&
               & nb,MaxActualDimAlpha,batchsizeAlpha,batchdimAlpha,batchindexAlpha,nbatchesAlpha,orb2BatchAlpha,'R')
       ENDIF

       if(master.and.DECinfo%PL>1)write(DECinfo%output,*) 'BATCH: Number of Alpha batches   = ', nbatchesAlpha&
          &, 'with maximum size',MaxActualDimAlpha

       IF(.NOT.DECinfo%useIchor)THEN
          ! Translate batchindex to orbital index
          ! -------------------------------------
          call mem_alloc(batch2orbAlpha,nbatchesAlpha)
          do idx=1,nbatchesAlpha
             call mem_alloc(batch2orbAlpha(idx)%orbindex,batchdimAlpha(idx) )
             batch2orbAlpha(idx)%orbindex = 0
             batch2orbAlpha(idx)%norbindex = 0
          end do
          do iorb=1,nb
             idx = orb2batchAlpha(iorb)
             batch2orbAlpha(idx)%norbindex = batch2orbAlpha(idx)%norbindex+1
             K = batch2orbAlpha(idx)%norbindex
             batch2orbAlpha(idx)%orbindex(K) = iorb
          end do
       ENDIF

       ! ************************************************
       ! *  Allocate matrices used in the batched loop  *
       ! ************************************************

       ! Use the dense amplitudes
       ! ------------------------

       !get the t+ and t- for the Kobayshi-like B2 term
       call tensor_ainit(tpl,[nor,nvr],2,atype=def_atype)                                                
       call tensor_ainit(tmi,[nor,nvr],2,atype=def_atype)


       !if I am the working process, then
       call get_tpl_and_tmi(t2,tpl,tmi)

       if(master.and.print_debug)then
          call print_norm(tpl," NORM(tpl)   :")
          call print_norm(tmi," NORM(tmi)    :")
       endif


       !get u2 in pdm or local
       call tensor_ainit(u2, [nv,nv,no,no], 4, local=local, atype='LDAR' )
       !calculate u matrix: t[c d i j] -> t[d c i j], 2t[d c i j] - t[d c j i] = u [d c i j]
       call array_reorder_4d(  2.0E0_realk, t2%elm1,nv,nv,no,no,[2,1,3,4],0.0E0_realk,u2%elm1)
       call array_reorder_4d( -1.0E0_realk, t2%elm1,nv,nv,no,no,[2,1,4,3],1.0E0_realk,u2%elm1)

       if(print_debug) call print_norm(u2," NORM(u2)    :",print_on_rank=0)

        
       call tensor_ainit(sio4, [no,no,nor], 3, local=local, atype='LDAR' )
       call tensor_zero(sio4)

       !allocate working arrays --> adapt this to memory requirements
       call mem_alloc(w0,i8*nb**4)
       call mem_alloc(w1,i8*nb**4)
       call mem_alloc(w2,i8*nb**4)
       call mem_alloc(w3,i8*nb**4)

       IF(DECinfo%useIchor)THEN
          !Calculate Screening integrals 
          SameMOL = .TRUE. !Specifies same molecule on all centers 
          call SCREEN_ICHORERI_DRIVER(DECinfo%output,iprint,mylsitem%setting,INTSPEC,SameMOL)
       ELSE
          ! This subroutine builds the full screening matrix.
          call II_precalc_DECScreenMat(DECscreen,DECinfo%output,6,mylsitem%setting,&
               & nbatchesAlpha,nbatchesGamma,INTSPEC,DECinfo%IntegralThreshold)
          IF(mylsitem%setting%scheme%cs_screen .OR. mylsitem%setting%scheme%ps_screen)THEN
             call II_getBatchOrbitalScreen(DecScreen,mylsitem%setting,&
                  & nb,nbatchesAlpha,nbatchesGamma,&
                  & batchsizeAlpha,batchsizeGamma,batchindexAlpha,batchindexGamma,&
                  & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
             call II_getBatchOrbitalScreenK(DecScreen,mylsitem%setting,&
                  & nb,nbatchesAlpha,nbatchesGamma,batchsizeAlpha,batchsizeGamma,&
                  & batchindexAlpha,batchindexGamma,&
                  & batchdimAlpha,batchdimGamma,INTSPEC,DECinfo%output,DECinfo%output)
          ENDIF
          !setup LHS screening - the full AO basis is used so we can use the
          !                      full matrices:        FilenameCS and FilenamePS
          !Note that it is faster to calculate the integrals in the form
          !(dimAlpha,dimGamma,nbasis,nbasis) so the full AO basis is used on the RHS
          !but the integrals is stored and returned in (nbasis,nbasis,dimAlpha,dimGamma)
       ENDIF

#ifdef VAR_OMP
       nthreads=OMP_GET_MAX_THREADS()
       if(master.and.DECinfo%PL>2)write(DECinfo%output,*)&
          & 'Starting CCSD residuals - OMP. Number of threads: ', OMP_GET_MAX_THREADS()
#else
       nthreads=1
       if(master.and.DECinfo%PL>2)write(DECinfo%output,*) &
          &'Starting CCSD integral/amplitudes - NO OMP!'
#endif


!#ifdef VAR_MPI
!       if(.not.dynamic_load)then
!
!          ! Calculate the batches for a good load balance
!          lenI2 = nbatchesAlpha*nbatchesGamma
!          !call mem_alloc( tasks, tasksc, lenI2 )
!          call mem_alloc( tasks, lenI2 )
!
!          myload = 0
!          tasks  = 0
!
!          IF(DECinfo%useIchor)THEN
!             call mem_alloc(batchdimAlpha,nbatchesAlpha)
!             do idx=1,nbatchesAlpha
!                batchdimAlpha(idx) = AOAlphabatchinfo(idx)%dim 
!             enddo
!             call mem_alloc(batch2orbAlpha,nbatchesAlpha)
!             do idx=1,nbatchesAlpha
!                call mem_alloc(batch2orbAlpha(idx)%orbindex,1)
!                batch2orbAlpha(idx)%orbindex(1) = AOAlphabatchinfo(idx)%orbstart
!                batch2orbAlpha(idx)%norbindex = 1
!             end do
!             call mem_alloc(batchdimGamma,nbatchesGamma)
!             do idx=1,nbatchesGamma
!                batchdimGamma(idx) = AOGammabatchinfo(idx)%dim 
!             enddo
!             call mem_alloc(batch2orbGamma,nbatchesGamma)
!             do idx=1,nbatchesGamma
!                call mem_alloc(batch2orbGamma(idx)%orbindex,1)
!                batch2orbGamma(idx)%orbindex(1) = AOGammabatchinfo(idx)%orbstart
!                batch2orbGamma(idx)%norbindex = 1
!             end do
!             call distribute_mpi_jobs(tasks,nbatchesAlpha,nbatchesGamma,&
!                  & batchdimAlpha,batchdimGamma,myload,lg_nnod,lg_me,scheme,&
!                  & no,nv,nb,batch2orbAlpha,batch2orbGamma)
!             call mem_dealloc(batchdimAlpha)
!             call mem_dealloc(batchdimGamma)
!             do idx=1,nbatchesAlpha
!                call mem_dealloc(batch2orbAlpha(idx)%orbindex)
!             end do
!             call mem_dealloc(batch2orbAlpha)
!             do idx=1,nbatchesGamma
!                call mem_dealloc(batch2orbGamma(idx)%orbindex)
!             end do
!             call mem_dealloc(batch2orbGamma)
!          ELSE
!             call distribute_mpi_jobs(tasks,nbatchesAlpha,nbatchesGamma,batchdimAlpha,&
!                  &batchdimGamma,myload,lg_nnod,lg_me,scheme,no,nv,nb,batch2orbAlpha,&
!                  &batch2orbGamma)
!          ENDIF
!       else
!
!          lenI2 = nbatchesGamma
!
!          call mem_alloc( tasks, tasksc, lenI2 ) 
!
!          tasks = 0
!          if(lg_me == 0) tasks(1) = lg_nnod
!
!          call lsmpi_win_create(tasks,tasksw,nbatchesGamma,infpar%lg_comm)
!#ifdef VAR_HAVE_MPI3
!          call lsmpi_win_lock_all(tasksw,ass=mode)
!#endif
!
!       endif
!
!       
!       if(master.and.DECinfo%PL>2)then
!          write(*,'("CCSD time in lsmpi_win_unlock phase A",g10.3)') time_lsmpi_win_unlock - unlock_time
!          write(*,'("CCSD time in lsmpi_wait       phase A",g10.3)') time_lsmpi_wait       - waiting_time
!          write(*,'("CCSD time in lsmpi_win_flush  phase A",g10.3)') time_lsmpi_win_flush  - flushing_time
!       endif
!       unlock_time   = time_lsmpi_win_unlock 
!       waiting_time  = time_lsmpi_wait
!       flushing_time = time_lsmpi_win_flush
!#endif

       myload = 0

       fullRHS=(nbatchesGamma.EQ.1).AND.(nbatchesAlpha.EQ.1)

       !**********************************
       ! Begin the loop over gamma batches
       !**********************************

       first_round=.false.
       if(dynamic_load)first_round=.true.

       BatchGamma: do gammaB = 1,nbatchesGamma  ! AO batches
          IF(DECinfo%useIchor)THEN
             dimGamma     = AOGammabatchinfo(gammaB)%dim      ! Dimension of gamma batch
             GammaStart   = AOGammabatchinfo(gammaB)%orbstart ! First orbital index in gamma batch
             GammaEnd     = AOGammabatchinfo(gammaB)%orbEnd   ! Last orbital index in gamma batch
             AOGammaStart = AOGammabatchinfo(gammaB)%AOstart  ! First AO batch index in gamma batch
             AOGammaEnd   = AOGammabatchinfo(gammaB)%AOEnd    ! Last AO batch index in gamma batch
          ELSE
             dimGamma     = batchdimGamma(gammaB)                         ! Dimension of gamma batch
             GammaStart   = batch2orbGamma(gammaB)%orbindex(1)            ! First index in gamma batch
             GammaEnd     = batch2orbGamma(gammaB)%orbindex(dimGamma)     ! Last index in gamma batch
          ENDIF
          !short hand notation
          fg         = GammaStart
          lg         = dimGamma


          alphaB=0
             

          !**********************************
          ! Begin the loop over alpha batches
          !**********************************


          BatchAlpha: do while(alphaB<=nbatchesAlpha) ! AO batches

             !check if the current job is to be done by current node
             !call check_job(scheme,first_round,dynamic_load,alphaB,gammaB,nbatchesAlpha,&
             !   &nbatchesGamma,tasks,tasksw,print_debug)

             !break the loop if alpha become too large, necessary to account for all
             !of the mpi and non mpi schemes, this is accounted for, because static,
             !and dynamic load balancing are enabled
             if(alphaB>nbatchesAlpha) exit
             
             IF(DECinfo%useIchor)THEN
                dimAlpha = AOAlphabatchinfo(alphaB)%dim         ! Dimension of alpha batch
                AlphaStart = AOAlphabatchinfo(alphaB)%orbstart  ! First orbital index in alpha batch
                AlphaEnd = AOAlphabatchinfo(alphaB)%orbEnd      ! Last orbital index in alpha batch
                AOAlphaStart = AOAlphabatchinfo(alphaB)%AOstart ! First AO batch index in alpha batch
                AOAlphaEnd = AOAlphabatchinfo(alphaB)%AOEnd     ! Last AO batch index in alpha batch
             ELSE
                dimAlpha   = batchdimAlpha(alphaB)                              ! Dimension of alpha batch
                AlphaStart = batch2orbAlpha(alphaB)%orbindex(1)                 ! First index in alpha batch
                AlphaEnd   = batch2orbAlpha(alphaB)%orbindex(dimAlpha)          ! Last index in alpha batch
             ENDIF

             !short hand notation
             fa         = AlphaStart
             la         = dimAlpha
             myload     = myload + la * lg


             !setup RHS screening - here we only have a set of AO basisfunctions
             !                      so we use the batchscreening matrices.
             !                      like BatchfilenamesCS(alphaB,gammaB)
             !Note that it is faster to calculate the integrals in the form
             !(dimAlpha,dimGamma,nbasis,nbasis) so the subset of the AO basis is used on the LHS
             !but the integrals is stored and returned in (nbasis,nbasis,dimAlpha,dimGamma)
!             call time_start_phase(PHASE_WORK, at = time_intloop_work)
!             IF(DECinfo%useIchor)THEN
!                dim1 = nb*nb*dimAlpha*dimGamma   ! dimension for integral array
!                call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,nb,nb,dimAlpha,dimGamma,&
!                     & w1%d,INTSPEC,FULLRHS,1,nAObatches,1,nAObatches,AOAlphaStart,AOAlphaEnd,&
!                     & AOGammaStart,AOGammaEnd,MoTrans,nb,nb,dimAlpha,dimGamma,NoSymmetry,DECinfo%IntegralThreshold)
!             ELSE
!                IF(doscreen) Mylsitem%setting%LST_GAB_LHS => DECSCREEN%masterGabLHS
!                IF(doscreen) mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGab(alphaB,gammaB)%p
!                ! Get (beta delta | alphaB gammaB) integrals using (beta,delta,alphaB,gammaB) ordering
!                ! ************************************************************************************
!                dim1 = nb*nb*dimAlpha*dimGamma   ! dimension for integral array
!                ! Store integral in tmp1(1:dim1) array in (beta,delta,alphaB,gammaB) order
!                call LSTIMER('START',tcpu1,twall1,DECinfo%output)
!                !Mylsitem%setting%scheme%intprint=6
!                call II_GET_DECPACKED4CENTER_J_ERI(DECinfo%output,DECinfo%output, Mylsitem%setting, w1%d,batchindexAlpha(alphaB),&
!                     &batchindexGamma(gammaB),&
!                     &batchsizeAlpha(alphaB),batchsizeGamma(gammaB),nb,nb,dimAlpha,dimGamma,fullRHS,INTSPEC,DECinfo%IntegralThreshold)
!                !Mylsitem%setting%scheme%intprint=0
!             ENDIF
             call LSTIMER('START',tcpu2,twall2,DECinfo%output)

             !GET EXCHANGE DISTRIBUTION TYPE ERI
             call time_start_phase(PHASE_WORK, at = time_intloop_work)
             IF(DECinfo%useIchor)THEN
                !Build (batchA,full,batchC,full)
                call MAIN_ICHORERI_DRIVER(DECinfo%output,iprint,Mylsitem%setting,dimAlpha,nb,dimGamma,nb,&
                     & w1%d,INTSPEC,FULLRHS,AOAlphaStart,AOAlphaEnd,1,nAObatches,AOGammaStart,AOGammaEnd,&
                     & 1,nAObatches,MoTrans,dimAlpha,nb,dimGamma,nb,NoSymmetry,DECinfo%IntegralThreshold)
             ELSE
                IF(doscreen)Mylsitem%setting%LST_GAB_LHS => DECSCREEN%batchGabKLHS(alphaB)%p
                IF(doscreen)Mylsitem%setting%LST_GAB_RHS => DECSCREEN%batchGabKRHS(gammaB)%p
                
                call II_GET_DECPACKED4CENTER_K_ERI(DECinfo%output,DECinfo%output, &
                     & Mylsitem%setting,w1%d,batchindexAlpha(alphaB),batchindexGamma(gammaB),&
                     & batchsizeAlpha(alphaB),batchsizeGamma(gammaB),dimAlpha,nb,dimGamma,nb,&
                     & INTSPEC,fullRHS,DECinfo%IntegralThreshold)
             ENDIF
             call lsmpi_poke()


             if( Ccmodel > MODEL_CC2 )then
                if(fa<=fg+lg-1)then
                   !CHECK WHETHER THE TERM HAS TO BE DONE AT ALL, i.e. when the first
                   !element in the alpha batch has a smaller index as the last element in
                   !the gamma batch, chose the trafolength as minimum of alpha batch-length
                   !and the difference between first element of alpha batch and last element
                   !of gamma batch
                   call get_a22_and_prepb22_terms_ex(w0%d,w1%d,w2%d,w3%d,tpl,tmi,no,nv,nb,fa,fg,la,lg,&
                      &xo,yo,xv,yv,omega2,sio4,scheme,[w0%n,w1%n,w2%n,w3%n],lock_outside,&
                      &time_intloop_B1work, time_intloop_B1comm, scal=0.5E0_realk  )


                   !start a new timing phase after these terms
                   call time_start_phase(PHASE_WORK)

                endif
             endif



          end do BatchAlpha
       end do BatchGamma


       ! Free integral stuff
       ! *******************
       IF(DECinfo%useIchor)THEN
          call FREE_SCREEN_ICHORERI()
          call mem_dealloc(AOGammabatchinfo)
          call mem_dealloc(AOAlphabatchinfo)
       ELSE
          nullify(Mylsitem%setting%LST_GAB_LHS)
          nullify(Mylsitem%setting%LST_GAB_RHS)
          call free_decscreen(DECSCREEN)
          
          ! Free gamma stuff
          call mem_dealloc(orb2batchGamma)
          call mem_dealloc(batchdimGamma)
          call mem_dealloc(batchsizeGamma)
          call mem_dealloc(batchindexGamma)
          do idx=1,nbatchesGamma
             call mem_dealloc(batch2orbGamma(idx)%orbindex)
             batch2orbGamma(idx)%orbindex => null()
          end do
          call mem_dealloc(batch2orbGamma)
          
          ! Free alpha stuff
          call mem_dealloc(orb2batchAlpha)
          call mem_dealloc(batchdimAlpha)
          call mem_dealloc(batchsizeAlpha)
          call mem_dealloc(batchindexAlpha)
          do idx=1,nbatchesAlpha
             call mem_dealloc(batch2orbAlpha(idx)%orbindex)
             batch2orbAlpha(idx)%orbindex => null()
          end do
          call mem_dealloc(batch2orbAlpha)
       ENDIF

       !FREE STUFF DEPENDING ON WHAT YOU NEED

!       !CALC T_1 transformed density matrix
!       call dgemm('n','t',nb,nb,no,1.0E0_realk,yo,nb,xo,nb,0.0E0_realk,w0%d,nb)
!
!       !GET 2electron part of inactive T_1 fock matrix
!       call II_get_fock_mat_full(DECinfo%output,DECinfo%output,MyLsItem%setting,nb,w0%d,.false.,w2%d)
!
!       !Get 1 electron part of ...
!       call ii_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsItem%setting,&
!          & w0%d,nb,nb,AORdefault,AORdefault)
!
!       ! Add one- and two-electron contributions to Fock matrix -> maybe  we can get rid of this
!       call daxpy(nb2,1.0E0_realk,w0%d,1,w2%d,1)
!
!       !Free the density matrix
!       call mem_dealloc(w0)
!
!       ! KK: Add long-range Fock correction
!       call daxpy(nb2,1.0E0_realk,deltafock,1,w2%d,1)
!
!       call time_start_phase(PHASE_WORK, ttot = time_get_ao_fock, twall = time_get_mo_fock)
!
!       if(print_debug)then
!          call print_norm(deltafock,int((i8*nb)*nb,kind=8), " NORM(deltafock):")
!          call print_norm(w2%d,int((i8*nb)*nb,kind=8)," NORM(iFock)    :")
!       endif
!
!
!
!       !Transform inactive Fock matrix into the different mo subspaces
!       if (Ccmodel>MODEL_CC2) then
!          ! -> Foo
!          call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,w2%d,nb,0.0E0_realk,w1%d,no)
!          call dgemm('n','n',no,no,nb,1.0E0_realk,w1%d,no,yo,nb,0.0E0_realk,ppfock,no)
!          ! -> Fov
!          call dgemm('n','n',no,nv,nb,1.0E0_realk,w1%d,no,yv,nb,0.0E0_realk,pqfock,no)
!          ! -> Fvo
!          call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,w2%d,nb,0.0E0_realk,w1%d,nv)
!          call dgemm('n','n',nv,no,nb,1.0E0_realk,w1%d,nv,yo,nb,0.0E0_realk,qpfock,nv)
!          ! -> Fvv
!          call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1%d,nv,yv,nb,0.0E0_realk,qqfock,nv)
!       else
!          ! -> Foo
!          call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,fock,nb,0.0E0_realk,w1%d,no)
!          call dgemm('n','n',no,no,nb,1.0E0_realk,w1%d,no,yo,nb,0.0E0_realk,ppfock,no)
!          ! -> Fov
!          call dgemm('t','n',no,nb,nb,1.0E0_realk,xo,nb,w2%d,nb,0.0E0_realk,w1%d,no)
!          call dgemm('n','n',no,nv,nb,1.0E0_realk,w1%d,no,yv,nb,0.0E0_realk,pqfock,no)
!          ! -> Fvo
!          call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,w2%d,nb,0.0E0_realk,w1%d,nv)
!          call dgemm('n','n',nv,no,nb,1.0E0_realk,w1%d,nv,yo,nb,0.0E0_realk,qpfock,nv)
!          ! -> Fvv
!          call dgemm('t','n',nv,nb,nb,1.0E0_realk,xv,nb,fock,nb,0.0E0_realk,w1%d,nv)
!          call dgemm('n','n',nv,nv,nb,1.0E0_realk,w1%d,nv,yv,nb,0.0E0_realk,qqfock,nv)
!       endif



!       if(print_debug)then
!          call print_norm(ppfock,int((i8*no)*no,kind=8)," NORM(ppfock)   :")
!          call print_norm(pqfock,int((i8*no)*nv,kind=8)," NORM(pqfock)   :")
!          call print_norm(qpfock,int((i8*no)*nv,kind=8)," NORM(qpfock)   :")
!          call print_norm(qqfock,int((i8*nv)*nv,kind=8)," NORM(qqfock)   :")
!       endif

!       !Free the AO fock matrix
!       call mem_dealloc(w2)

!       if(master.and.DECinfo%PL>1) then
!          call DetectHeapMemory(mem_allocated,HeapMemoryUsage,&
!               & 'get_ccsd_residual_integral_driven',DECinfo%output)
!          WRITE(DECinfo%output,'(A,f9.3,A)')'Expected Memory Used ',ActuallyUsed,' GB'
!       endif
    end subroutine get_ccsd_lhtr_integral_driven


    !> \brief Wrapper to calculate Calculate CCSD Jacobian left-hand transformation.
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine cc_jacobian_lhtr(mylsitem,Co,Cv,xo,xv,yo,yv,t1fock,&
         & t1,t2,L1,L2,rho1,rho2)
      implicit none
      !> LS item structure
      type(lsitem), intent(inout) :: MyLsItem      
      !> Occ and virt MO coefficients. Note: Co=xo (particle) and Cv=yv (hole)
      !> (in practice intent(in))
      type(tensor),intent(inout) :: Co,Cv
      !> Particle (x) and hole (y) transformation matrices
      !> (in practice intent(in))
      type(tensor),intent(inout) :: xo,xv,yo,yv
      ! T1-transformed Fock matrix in AO basis
      !> (in practice intent(in))
      type(tensor) :: t1fock
      !> Singles and doubles amplitudes  (in practice these are intent(in))
      !> (in practice intent(in))
      type(tensor),intent(inout) :: t1,t2
      !> Singles (L1) and doubles (L2) components of trial vector
      !> (in practice intent(in))
      type(tensor),intent(inout) :: L1,L2
      !> Output: Singles (rho1) and doubles (rho2) components of Jacobian transformation
      !> on trial vector (L1,L2).
      type(tensor),intent(inout) :: rho1,rho2
      integer :: nbasis,nocc,nvirt

      nbasis = Co%dims(1)
      nocc = Co%dims(2)
      nvirt = Cv%dims(2)

      ! Calculate left-hand transformation for Jacobian
      call get_ccsd_multipliers_simple(rho1%elm2,rho2%elm4,t1%elm2,t2%elm4,&
           & L1%elm2,L2%elm4,&
           & t1fock%elm2,xo%elm2,yo%elm2,xv%elm2,yv%elm2,&
           & nocc,nvirt,nbasis,MyLsItem,JacobianLT=.true.)

      ! Add contribution to get LW1 model
      if(DECinfo%ccModel==MODEL_MP2 .and. DECinfo%LW1) then
         call add_lw1_contribution(nbasis,mylsitem,xo,&
              & yv,t2,L1,L2,rho1,rho2,.false.)
      end if

    end subroutine cc_jacobian_lhtr



    !> \brief Calculate CCSD Jacobian right-hand transformation.
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine cc_jacobian_rhtr(mylsitem,xo,xv,yo,&
         & yv,t1,t2,R1,R2,rho1,rho2)
      implicit none
      !> LS item structure
      type(lsitem), intent(inout) :: MyLsItem      
      !> Particle (x) and hole (y) transformation matrices for occ (dimension nbasis,nocc)
      !> and virt (dimension nbasis,nvirt) transformations
      type(tensor),intent(in) :: xo,xv,yo,yv
      !> Singles and doubles amplitudes
      type(tensor),intent(in) :: t1,t2 
      !> Singles (R1) and doubles (R2) components of trial vector
      type(tensor),intent(in) :: R1,R2
      !> Singles (rho1) and doubles (rho2) components of Jacobian transformation on trial vector.
      type(tensor),intent(inout) :: rho1,rho2
      type(tensor) :: rho11,rho12,rho21,rho22
      integer :: nbasis,nocc,nvirt,whattodo

      ! Dimensions
      nbasis = xo%dims(1)
      nocc = xo%dims(2)
      nvirt = xv%dims(2)

      call tensor_minit(rho11,[nvirt,nocc],2)
      call tensor_minit(rho12,[nvirt,nocc,nvirt,nocc],4)
      call tensor_minit(rho21,[nvirt,nocc],2)
      call tensor_minit(rho22,[nvirt,nocc,nvirt,nocc],4)


      ! Calculate 1^rho components for Jacobian RHS transformation,
      ! see Eqs. 55 and 57 in JCP 105, 6921 (1996)     
      ! Singles component: rho11  (Eq. 55)
      ! Doubles component: rho12  (Eq. 57)
      whattodo=1
      call jacobian_rhtr_workhorse(mylsitem,xo,xv,yo,&
           & yv,t2,R1,R2,rho11,rho12,whattodo)

      ! Calculate 2^rho components for Jacobian RHS transformation,
      ! see Eqs. 56 and 58 in JCP 105, 6921 (1996)     
      ! Singles component: rho21  (Eq. 56)
      ! Doubles component: rho22  (Eq. 58)
      whattodo=2
      call jacobian_rhtr_workhorse(mylsitem,xo,xv,yo,&
           & yv,t2,R1,R2,rho21,rho22,whattodo)

      ! Add rho contributions (Eq. 34 in JCP 105, 6921 (1996))
      call tensor_zero(rho1)
      call tensor_zero(rho2)
      call tensor_add(rho1,1.0_realk,rho11)
      call tensor_add(rho1,1.0_realk,rho21)
      call tensor_add(rho2,1.0_realk,rho12)
      call tensor_add(rho2,1.0_realk,rho22)

      call tensor_free(rho11)
      call tensor_free(rho12)
      call tensor_free(rho21)
      call tensor_free(rho22)


    end subroutine cc_jacobian_rhtr


    !> \brief Noddy implementation of CCSD residual made more general such that the components of 
    !> the Jacobian right-hand transformation can also be calculated.
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine jacobian_rhtr_workhorse(mylsitem,xo,xv,yo,&
         & yv,t2,R1,R2,rho1,rho2,whattodo)
      implicit none
      !> LS item structure
      type(lsitem), intent(inout) :: MyLsItem      
      !> Particle (x) and hole (y) transformation matrices for occ (dimension nbasis,nocc)
      !> and virt (dimension nbasis,nvirt) transformations
      type(tensor),intent(in) :: xo,xv,yo,yv
      !> Doubles amplitudes
      type(tensor),intent(in),target :: t2 
      !> Singles (R1) and doubles (R2) components of trial vector
      type(tensor),intent(in),target :: R1,R2
      !> Singles (rho1) and doubles (rho2) components of Jacobian transformation on trial vector.
      type(tensor),intent(inout) :: rho1,rho2
      !> What to do?
      !> 1. Calculate 1^rho components for Jacobian RH transformation (singles and doubles), 
      !>    see Eqs. 55 and 57 in JCP 105, 6921 (1996)
      !> 2. Calculate 2^rho components for Jacobian RH transformation (singles and doubles), 
      !>    see Eqs. 56 and 58 in JCP 105, 6921 (1996)
      !> 3. Calculate standard CCSD residual (then R1=CCSD singles, R2=CCSD doubles=t2), 
      !>    see Eqs. 13.7.80-83 and 13.7.101-105 in THE BOOK 
      integer,intent(in) :: whattodo
      type(array2) :: fvo,fov,fvv,foo
      type(array4) :: gvvov,gooov,gvovo,gvvvv,goooo,govov,goovv,gvoov
      integer :: nbasis,nocc,nvirt,i,j,k,l,m,a,b,c,d,e
      real(realk) :: uA,uB,Lldkc,fac
      real(realk),pointer :: tmp(:,:,:,:),tmp2(:,:,:,:),rho2CDE(:,:,:,:),tmpvv(:,:),tmpoo(:,:)
      logical :: somethingwrong,includeterm
      real(realk),pointer :: A2(:,:,:,:), B2(:,:,:,:)
      type(tensor) :: ttmp


      ! Dimensions
      nbasis = xo%dims(1)
      nocc = xo%dims(2)
      nvirt = xv%dims(2)

      ! Sanity check 1
      if(whattodo/=1 .and. whattodo/=2 .and. whattodo/=3) then
         print *, 'whattodo ', whattodo
         call lsquit('jacobian_rhtr_workhorse: Ill-defined whattodo!',-1)
      end if

      ! Sanity check 2: Dimensions
      somethingwrong=.false.
      if(t2%dims(1)/=nvirt .or. t2%dims(3)/=nvirt) somethingwrong=.true.
      if(t2%dims(2)/=nocc .or. t2%dims(4)/=nocc) somethingwrong=.true.
      if(R2%dims(1)/=nvirt .or. R2%dims(3)/=nvirt) somethingwrong=.true.
      if(R2%dims(2)/=nocc .or. R2%dims(4)/=nocc) somethingwrong=.true.
      if(rho2%dims(1)/=nvirt .or. rho2%dims(3)/=nvirt) somethingwrong=.true.
      if(rho2%dims(2)/=nocc .or. rho2%dims(4)/=nocc) somethingwrong=.true.
      if(R1%dims(1)/=nvirt .or. R1%dims(2)/=nocc) somethingwrong=.true.
      if(rho1%dims(1)/=nvirt .or. rho1%dims(2)/=nocc) somethingwrong=.true.
      if(somethingwrong) then
         print *, 't2   ', t2%dims
         print *, 'R2   ', R2%dims
         print *, 'rho2 ', rho2%dims
         print *, 'rho1 ', rho1%dims
         print *, 'R1   ', R1%dims
         call lsquit('jacobian_rhtr_workhorse: Dimension error!',-1)
      end if


      ! Calculate all integrals 
      call jacobian_rhtr_workhorse_integrals(mylsitem,xo,xv,yo,&
           & yv,R1,gvvov,gooov,gvovo,gvvvv,goooo,govov,goovv,gvoov, &
           & fvo,fov,fvv,foo,whattodo)


      ! Calculate terms in CCSD residual
      ! ********************************
      ! Notation according to THE BOOK
      ! Corresponding notation in Table I of JCP 105, 6921 (1996) is given in parenthesis
      call tensor_zero(rho1)
      call tensor_zero(rho2)


      ! SINGLES RESIDUAL
      ! ================

      ! We use generic doubles array A2, see Table I of JCP 105, 6921 (1996):
      ! whattodo = 1: Doubles array is t2
      ! whattodo = 2: Doubles array is R2
      ! whattodo = 3: Doubles array is t2
      if(whattodo==1 .or. whattodo==3) A2 => t2%elm4
      if(whattodo==2) A2 => R2%elm4

      do i=1,nocc
         do a=1,nvirt

            ! A1 term (G term in JCP 105, 6921 (1996))
            do c=1,nvirt
               do d=1,nvirt
                  do k=1,nocc
                     uA = 2.0_realk*A2(c,k,d,i) - A2(c,i,d,k)
                     rho1%elm2(a,i) = rho1%elm2(a,i) + uA*gvvov%val(a,d,k,c)
                  end do
               end do
            end do

            ! B1 term (H)
            do c=1,nvirt
               do k=1,nocc
                  do l=1,nocc
                     uA = 2.0_realk*A2(a,k,c,l) - A2(a,l,c,k)
                     rho1%elm2(a,i) = rho1%elm2(a,i) - uA*gooov%val(k,i,l,c)
                  end do
               end do
            end do

            ! C1 term (I)
            do c=1,nvirt
               do k=1,nocc
                  uA = 2.0_realk*A2(a,i,c,k) - A2(a,k,c,i)
                  rho1%elm2(a,i) = rho1%elm2(a,i) + uA*fov%val(k,c)
               end do
            end do

         end do
      end do

      ! D1 term (J)
      if(whattodo/=2) then
         ! Not include zero-order term for 2^rho singles component 
         ! (Eq. 37 in JCP 105, 6921 (1996))
         do i=1,nocc
            do a=1,nvirt
               rho1%elm2(a,i) = rho1%elm2(a,i) + fvo%val(a,i)
            end do
         end do
      end if


      ! DOUBLES RESIDUAL
      ! ================

      ! We use generic doubles arrays A2 and B2, see Table I of JCP 105, 6921 (1996):
      ! whattodo = 1: Doubles arrays A2 and B2 are both t2
      ! whattodo = 2: Doubles arrays A2 and B2 are R2 and t2, respectively
      ! whattodo = 3: Doubles arrays A2 and B2 are both t2
      if(whattodo==1 .or. whattodo==3) then
         A2 => t2%elm4
         B2 => t2%elm4
      else
         A2 => R2%elm4
         B2 => t2%elm4
      end if


      ! A2 (F and B)
      ! ************

      ! A2.1 term (F term in JCP 105, 6921 (1996))
      if(whattodo/=2) then
         ! Not include zero-order term for 2^rho doubles component 
         ! (Eq. 38 in JCP 105, 6921 (1996))
         do i=1,nocc
            do a=1,nvirt
               do j=1,nocc
                  do b=1,nvirt

                     rho2%elm4(a,i,b,j) = gvovo%val(a,i,b,j)

                  end do
               end do
            end do
         end do
      end if

      ! For whattodo=1 and Hald approximation all remaining terms are skipped
      Hald1: if( whattodo/=1 .or. (.not. DECinfo%haldapprox) ) then

         ! For whattodo=2 all terms involving doubles amplitudes (B2 pointer)
         ! are not included. This is controlled by "includeterm" logical.
         includeterm=.true.
         if(whattodo==2 .and. DECinfo%haldapprox) includeterm=.false.


         ! A2.2 term (B)
         do i=1,nocc
            do a=1,nvirt
               do j=1,nocc
                  do b=1,nvirt

                     do c=1,nvirt
                        do d=1,nvirt
                           rho2%elm4(a,i,b,j) = rho2%elm4(a,i,b,j) + &
                                & A2(c,i,d,j)*gvvvv%val(a,c,b,d)
                        end do
                     end do

                  end do
               end do
            end do
         end do


         ! B2 (A)
         ! ******
         ! Scaling of quadratic term (not present for 1^rho)
         if(whattodo==1) fac=0.0_realk
         if(whattodo==2 .or. whattodo==3) fac=1.0_realk

         call mem_alloc(tmp,nocc,nocc,nocc,nocc)
         do i=1,nocc
            do j=1,nocc

               do k=1,nocc
                  do l=1,nocc

                     ! Temporary quantity: 
                     ! tmp(k,i,l,j) = g(k,i,l,j) + fac * sum_{cd} B2(c,i,d,j) * g(k,c,l,d)
                     tmp(k,i,l,j) = goooo%val(k,i,l,j)

                     if(includeterm) then
                        do c=1,nvirt
                           do d=1,nvirt
                              tmp(k,i,l,j) = tmp(k,i,l,j) + fac*B2(c,i,d,j)*govov%val(k,c,l,d)
                           end do
                        end do
                     end if

                     ! rho2(a,b,i,j) += sum_{kl} A2(a,k,b,l) tmp(k,i,l,j)
                     ! Note that the linear term is then included with A2, not B2:
                     ! sum_{kl} A2(a,k,b,l) g(k,i,l,j)
                     do a=1,nvirt
                        do b=1,nvirt
                           rho2%elm4(a,i,b,j) = rho2%elm4(a,i,b,j) + A2(a,k,b,l)*tmp(k,i,l,j) 
                        end do
                     end do

                  end do
               end do

            end do
         end do



         ! For 2^rho, add quadratic term with B2 and A2 switched around,
         ! see Table I in JCP 105, 6921 (1996)
         Bquadratic: if(whattodo==2 .and. includeterm) then
            do i=1,nocc
               do j=1,nocc

                  do k=1,nocc
                     do l=1,nocc

                        ! Temporary quantity: tmp(k,i,l,j) = sum_{cd} A2(c,i,d,j) * g(k,c,l,d)
                        tmp(k,i,l,j) = 0.0_realk
                        do c=1,nvirt
                           do d=1,nvirt
                              tmp(k,i,l,j) = tmp(k,i,l,j) + fac*A2(c,i,d,j)*govov%val(k,c,l,d)
                           end do
                        end do

                        ! rho2(a,i,b,j) += sum_{kl} B2(a,k,b,l) tmp(k,i,l,j)
                        do a=1,nvirt
                           do b=1,nvirt
                              rho2%elm4(a,i,b,j) = rho2%elm4(a,i,b,j) + B2(a,k,b,l)*tmp(k,i,l,j) 
                           end do
                        end do

                     end do
                  end do

               end do
            end do
         end if Bquadratic
         call mem_dealloc(tmp)


         ! C2 (C)
         ! ******
         ! NOTE: Here the quadratic term is scaled by 0.5 for CCSD residual, 
         ! by 1.0 for 2^rho, and it is not present for 1^rho, see Table I in JCP 105, 6921 (1996)
         ! 
         if(whattodo==1) fac=0.0_realk
         if(whattodo==2) fac=1.0_realk
         if(whattodo==3) fac=0.5_realk
         call mem_alloc(tmp,nocc,nocc,nvirt,nvirt)
         call mem_alloc(tmp2,nvirt,nvirt,nocc,nocc)
         tmp2=0.0_realk
         do i=1,nocc
            do a=1,nvirt

               do k=1,nocc
                  do c=1,nvirt

                     ! tmp(k,i,a,c) = g(k,i,a,c) - fac * sum_{dl} B2(a,l,d,i) * govov(k,d,l,c)
                     tmp(k,i,a,c) = goovv%val(k,i,a,c)
                     if(includeterm) then
                        do l=1,nocc
                           do d=1,nvirt
                              tmp(k,i,a,c) = tmp(k,i,a,c) - fac*B2(a,l,d,i)*govov%val(k,d,l,c)
                           end do
                        end do
                     end if

                     ! tmp2(a,b,i,j) += - sum_{kc} A2(b,k,c,j) * tmp(k,i,a,c)
                     do j=1,nocc
                        do b=1,nvirt
                           tmp2(a,b,i,j) = tmp2(a,b,i,j) - A2(b,k,c,j) * tmp(k,i,a,c)
                        end do
                     end do

                  end do
               end do

            end do
         end do

         ! Add component with i<-->j scaled properly to get final C2 component
         ! (although still missing final symmetrization which is done commonly for C,D, and E
         !  terms at the very end)
         ! rho2CDE(a,i,b,j) = (0.5 + P_{ij}) * tmp2(a,b,i,j)
         call mem_alloc(rho2CDE,nvirt,nocc,nvirt,nocc)
         do j=1,nocc
            do i=1,nocc
               do a=1,nvirt
                  do b=1,nvirt
                     rho2CDE(a,i,b,j) = 0.5_realk*tmp2(a,b,i,j) + tmp2(a,b,j,i)
                  end do
               end do
            end do
         end do
         call mem_dealloc(tmp)
         call mem_dealloc(tmp2)


         ! D2 (D)
         ! ******  
         call mem_alloc(tmp,nvirt,nocc,nocc,nvirt)
         do i=1,nocc
            do a=1,nvirt


               do c=1,nvirt
                  do k=1,nocc

                     ! tmp(a,i,k,c) = L(a,i,k,c) + fac * sum_{dl} uB(a,i,d,l)*L(l,d,k,c)
                     ! ----------------------------------------------------------------

                     ! L(a,i,k,c) = 2*g(a,i,k,c) - g(a,c,k,i) = 2*g(a,i,k,c) - g(k,i,a,c) 
                     tmp(a,i,k,c) = 2.0_realk*gvoov%val(a,i,k,c) - goovv%val(k,i,a,c)
                     if(includeterm) then
                        do d=1,nvirt
                           do l=1,nocc
                              ! L(l,d,k,c) = 2*g(l,d,k,c) - g(l,c,k,d)
                              Lldkc = 2.0_realk*govov%val(l,d,k,c) - govov%val(l,c,k,d)
                              ! uB(a,d,i,l) = 2*B2(a,i,d,l) - B2(a,l,d,i)
                              uB = 2.0_realk*B2(a,i,d,l) - B2(a,l,d,i)
                              ! tmp(a,i,k,c) = L(a,i,k,c) + fac * sum_{dl} uB(a,i,d,l)*L(l,d,k,c) 
                              tmp(a,i,k,c) = tmp(a,i,k,c) + fac*uB*Lldkc
                           end do
                        end do
                     end if

                     ! Update:
                     ! rho2(a,b,i,j) += rho2(a,b,i,j) + 0.5 * sum_{ck} uA(b,j,c,k) * tmp(a,i,k,c)
                     do b=1,nvirt
                        do j=1,nocc
                           uA = 2.0_realk*A2(b,j,c,k) - A2(b,k,c,j)
                           rho2CDE(a,i,b,j) = rho2CDE(a,i,b,j) + 0.5_realk*uA*tmp(a,i,k,c)
                        end do
                     end do

                  end do
               end do


            end do
         end do
         call mem_dealloc(tmp)


         ! E2 (E1 and E2)
         ! **************
         ! Scaling of quadratic term (not present for 1^rho)
         if(whattodo==1) fac=0.0_realk
         if(whattodo==2 .or. whattodo==3) fac=1.0_realk

         ! Construct 2-dimensional intermediates
         call mem_alloc(tmpoo,nocc,nocc)
         call mem_alloc(tmpvv,nvirt,nvirt)

         ! tmpoo(k,j) = F(k,j) + sum_{cdl} uB(c,l,d,j) * g(k,d,l,c)
         do j=1,nocc
            do k=1,nocc
               tmpoo(k,j) = foo%val(k,j)

               if(includeterm) then
                  do l=1,nocc
                     do d=1,nvirt
                        do c=1,nvirt
                           uB = 2.0_realk*B2(c,l,d,j) - B2(c,j,d,l)
                           tmpoo(k,j) = tmpoo(k,j) + fac * uB * govov%val(k,d,l,c)
                        end do
                     end do
                  end do
               end if

            end do
         end do

         ! tmpvv(b,c) = F(b,c) - sum_{dkl} uB(b,k,d,l) * g(l,d,k,c)
         do b=1,nvirt
            do c=1,nvirt
               tmpvv(b,c) = fvv%val(b,c)

               if(includeterm) then
                  do l=1,nocc
                     do d=1,nvirt
                        do k=1,nocc
                           uB = 2.0_realk*B2(b,k,d,l) - B2(b,l,d,k)
                           tmpvv(b,c) = tmpvv(b,c) - fac * uB * govov%val(l,d,k,c)
                        end do
                     end do
                  end do
               end if

            end do
         end do

         ! E2.1 (E1)
         do j=1,nocc
            do c=1,nvirt
               do b=1,nvirt
                  do i=1,nocc
                     do a=1,nvirt
                        rho2CDE(a,i,b,j) = rho2CDE(a,i,b,j) + A2(a,i,c,j)*tmpvv(b,c)
                     end do
                  end do
               end do
            end do
         end do

         ! E2.2 (E2)
         do j=1,nocc
            do k=1,nocc
               do b=1,nvirt
                  do i=1,nocc
                     do a=1,nvirt
                        rho2CDE(a,i,b,j) = rho2CDE(a,i,b,j) - A2(a,i,b,k)*tmpoo(k,j)
                     end do
                  end do
               end do
            end do
         end do



         ! Add quadratic E terms with B2 and A2 switched around
         ! ----------------------------------------------------
         Equadratic: if(whattodo==2 .and. includeterm) then
            ! tmpoo(k,j) = sum_{cdl} uA(c,l,d,j) * g(k,d,l,c)
            do j=1,nocc
               do k=1,nocc
                  tmpoo(k,j) = 0.0_realk
                  do l=1,nocc
                     do d=1,nvirt
                        do c=1,nvirt
                           uA = 2.0_realk*A2(c,l,d,j) - A2(c,j,d,l)
                           tmpoo(k,j) = tmpoo(k,j) + fac * uA * govov%val(k,d,l,c)
                        end do
                     end do
                  end do
               end do
            end do

            ! tmpvv(b,c) = - sum_{dkl} uA(b,k,d,l) * g(l,d,k,c)
            do b=1,nvirt
               do c=1,nvirt
                  tmpvv(b,c) = 0.0_realk

                  do l=1,nocc
                     do d=1,nvirt
                        do k=1,nocc
                           uA = 2.0_realk*A2(b,k,d,l) - A2(b,l,d,k)
                           tmpvv(b,c) = tmpvv(b,c) - fac * uA * govov%val(l,d,k,c)
                        end do
                     end do
                  end do
               end do
            end do

            ! E2.1 (E1)
            do j=1,nocc
               do c=1,nvirt
                  do b=1,nvirt
                     do i=1,nocc
                        do a=1,nvirt
                           rho2CDE(a,i,b,j) = rho2CDE(a,i,b,j) + B2(a,i,c,j)*tmpvv(b,c)
                        end do
                     end do
                  end do
               end do
            end do

            ! E2.2 (E2)
            do j=1,nocc
               do k=1,nocc
                  do b=1,nvirt
                     do i=1,nocc
                        do a=1,nvirt
                           rho2CDE(a,i,b,j) = rho2CDE(a,i,b,j) - B2(a,i,b,k)*tmpoo(k,j)
                        end do
                     end do
                  end do
               end do
            end do

         end if Equadratic


         ! Symmetrize C, D, and E terms and add to A and B terms
         ! *****************************************************
         do j=1,nocc
            do b=1,nvirt
               do i=1,nocc
                  do a=1,nvirt
                     rho2%elm4(a,i,b,j) = rho2%elm4(a,i,b,j) + rho2CDE(a,i,b,j) + rho2CDE(b,j,a,i)
                  end do
               end do
            end do
         end do
         call mem_dealloc(rho2CDE)
         call mem_dealloc(tmpoo)
         call mem_dealloc(tmpvv)

      end if Hald1


      ! Add contribution to get LW1 model
      if(DECinfo%ccModel==MODEL_MP2 .and. whattodo==1 .and. DECinfo%LW1) then
         call add_lw1_contribution(nbasis,mylsitem,xo,&
              & yv,t2,R1,R2,rho1,rho2,.true.)
      end if


      ! Clean up
      call array4_free(gvvov)
      call array4_free(gooov)
      call array4_free(gvovo)
      call array4_free(gvvvv)
      call array4_free(goooo)
      call array4_free(govov)
      call array4_free(goovv)
      call array4_free(gvoov)
      call array2_free(foo)
      call array2_free(fov)
      call array2_free(fvo)
      call array2_free(fvv)
      nullify(A2)
      nullify(B2)

      
    end subroutine jacobian_rhtr_workhorse


    !> \brief Calculate integrals for jacobian_rhtr_workhorse.
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine jacobian_rhtr_workhorse_integrals(mylsitem,xo_tensor,xv_tensor,yo_tensor,&
         & yv_tensor,R1_tensor,gvvov,gooov,gvovo,gvvvv,goooo,govov,goovv,gvoov, &
         & fvo,fov,fvv,foo,whattodo)
      implicit none
      !> LS item structure
      type(lsitem), intent(inout) :: MyLsItem      
      !> Particle (x) and hole (y) transformation matrices for occ (dimension nbasis,nocc)
      !> and virt (dimension nbasis,nvirt) transformations
      type(tensor),intent(in) :: xo_tensor,xv_tensor,yo_tensor,yv_tensor
      !> Singles (R1) component of trial vector
      type(tensor),intent(in) :: R1_tensor
      !> Two-electron integrals in obvious notation (e.g. gvvov = (a b | i c) integrals)
      !> These are allocated inside this subroutine!
      type(array4),intent(inout) :: gvvov,gooov,gvovo,gvvvv,goooo,govov,goovv,gvoov
      !> Fock matrix elements in MO basis, e.g. fvo = < a | F | i >,
      !> using either T1-transformed integrals or trial-T1 transformed integrals (see below).
      !> These are allocated inside this subroutine!
      type(array2),intent(inout) :: fvo,fov,fvv,foo
      !> What to do? See jacobian_rhtr_workhorse.
      !> whattodo=1: Use trial-T1 transformed integrals 
      !>             (Hamiltonian is Eq. 45 in JCP 105, 6921 (1996))
      !> whattodo=2: Use T1-transformed integrals
      !> whattodo=3: Use T1-transformed integrals
      integer,intent(in) :: whattodo
      type(array2) :: xo,xv,yo,yv,fao,Dt1,hao,hoo,hvo,hov,hvv,xoR,xvR,yoR,yvR,R1
      type(array4) :: gao,gvooo
      integer :: nbasis,nocc,nvirt,i,j,a,b,k
      integer,dimension(2) :: bo,bv,bb,oo,ov,vo,vv
      real(realk) :: s

      ! Dimensions
      nbasis = xo_tensor%dims(1)
      nocc = xo_tensor%dims(2)
      nvirt = xv_tensor%dims(2)
      bo(1) = nbasis; bo(2) = nocc
      bv(1) = nbasis; bv(2) = nvirt
      bb(1) = nbasis; bb(2)=nbasis
      oo(1)=nocc; oo(2)=nocc
      ov(1)=nocc; ov(2)=nvirt
      vo(1)=nvirt; vo(2)=nocc
      vv(1)=nvirt; vv(2)=nvirt

      ! Transformation matrices in array2 format
      xo = array2_init(bo,xo_tensor%elm2)
      xv = array2_init(bv,xv_tensor%elm2)
      yo = array2_init(bo,yo_tensor%elm2)
      yv = array2_init(bv,yv_tensor%elm2)
      R1 = array2_init(vo,R1_tensor%elm2)

      ! Calculate full AO integrals
      call get_full_eri(mylsitem,nbasis,gao)

      ! Get MO integrals
      ! ****************
      ResidualT1_twoel: if(whattodo==1) then  ! Residual-T1 transformed integrals

         ! See Eq. 47 in JCP 105, 6921 (1996) for two-electron trial-T1 transformed integrals

         ! Residual-transformed coefficient matrices (Eqs. 49 and 50 in JCP 105, 6921 (1996))
         xoR = array2_init(bo)   ! zero, see Eqs. 49-51 in JCP 105, 6921 (1996) 
         xvR = array2_init(bv)
         yoR = array2_init(bo)   
         yvR = array2_init(bv)   ! zero, see Eqs. 49-51 in JCP 105, 6921 (1996) 
         call array2_matmul(xo,R1,xvR,'N','T',-1.0_realk,0.0_realk)
         call array2_matmul(yv,R1,yoR,'N','N',1.0_realk,0.0_realk)

         ! Calculate integrals according to Eq. 47 in JCP 105, 6921 (1996)
         call two_electron_trial_T1_integrals(gao,xv,yv,xo,yv,xvR,yvR,xoR,yvR,gvvov)
         call two_electron_trial_T1_integrals(gao,xo,yo,xo,yv,xoR,yoR,xoR,yvR,gooov)
         call two_electron_trial_T1_integrals(gao,xv,yo,xv,yo,xvR,yoR,xvR,yoR,gvovo)
         call two_electron_trial_T1_integrals(gao,xv,yv,xv,yv,xvR,yvR,xvR,yvR,gvvvv)
         call two_electron_trial_T1_integrals(gao,xo,yo,xo,yo,xoR,yoR,xoR,yoR,goooo)
         call two_electron_trial_T1_integrals(gao,xo,yv,xo,yv,xoR,yvR,xoR,yvR,govov)
         call two_electron_trial_T1_integrals(gao,xo,yo,xv,yv,xoR,yoR,xvR,yvR,goovv)
         call two_electron_trial_T1_integrals(gao,xv,yo,xo,yv,xvR,yoR,xoR,yvR,gvoov)
         call two_electron_trial_T1_integrals(gao,xv,yo,xo,yo,xvR,yoR,xoR,yoR,gvooo)

      else ! T1-transformed integrals

         gvvov = get_gmo_simple(gao,xv,yv,xo,yv)
         gooov = get_gmo_simple(gao,xo,yo,xo,yv)
         gvovo = get_gmo_simple(gao,xv,yo,xv,yo)
         gvvvv = get_gmo_simple(gao,xv,yv,xv,yv)
         goooo = get_gmo_simple(gao,xo,yo,xo,yo)
         govov = get_gmo_simple(gao,xo,yv,xo,yv)
         goovv = get_gmo_simple(gao,xo,yo,xv,yv)
         gvoov = get_gmo_simple(gao,xv,yo,xo,yv)
         gvooo = get_gmo_simple(gao,xv,yo,xo,yo)

      end if ResidualT1_twoel

      call array4_free(gao)


      ! Fock matrix
      ! ***********

      ! One-electron contribution
      hao = array2_init(bb)
      call ii_get_h1_mixed_full(DECinfo%output,DECinfo%output,MyLsItem%setting,&
           & hao%val,nbasis,nbasis,AORdefault,AORdefault)

      ResidualT1_oneel: if(whattodo==1) then  ! Residual-T1 transformed integrals
         call one_electron_trial_T1_integrals(hao,xo,yo,xoR,yoR,hoo)
         call one_electron_trial_T1_integrals(hao,xo,yv,xoR,yvR,hov)
         call one_electron_trial_T1_integrals(hao,xv,yo,xvR,yoR,hvo)
         call one_electron_trial_T1_integrals(hao,xv,yv,xvR,yvR,hvv)
      else
         hoo = array2_similarity_transformation(xo,hao,yo,oo)
         hov = array2_similarity_transformation(xo,hao,yv,ov)
         hvo = array2_similarity_transformation(xv,hao,yo,vo)
         hvv = array2_similarity_transformation(xv,hao,yv,vv)
      end if ResidualT1_oneel

      ! Occ-occ Fock matrix
      s = 2.0_realk
      foo = array2_init(oo,hoo%val)
      do j=1,nocc
         do i=1,nocc
            do k=1,nocc
               foo%val(i,j) = foo%val(i,j) + (s*goooo%val(i,j,k,k) - goooo%val(i,k,k,j))
            end do
         end do
      end do

      ! Occ-virt Fock matrix
      fov = array2_init(ov,hov%val)
      do a=1,nvirt
         do i=1,nocc
            do k=1,nocc
               fov%val(i,a) = fov%val(i,a) + (s*gooov%val(k,k,i,a) - gooov%val(i,k,k,a))
            end do
         end do
      end do

      ! Virt-occ Fock matrix
      fvo = array2_init(vo,hvo%val)
      do i=1,nocc
         do a=1,nvirt
            do k=1,nocc
               fvo%val(a,i) = fvo%val(a,i) + (s*gvooo%val(a,i,k,k) - gvooo%val(a,k,k,i))
            end do
         end do
      end do

      ! Virt-virt Fock matrix
      fvv = array2_init(vv,hvv%val)
      do b=1,nvirt
         do a=1,nvirt
            do k=1,nocc
               fvv%val(a,b) = fvv%val(a,b) + (s*goovv%val(k,k,a,b) - gvoov%val(a,k,k,b))
            end do
         end do
      end do

      call array2_free(hao)
      call array2_free(hoo)
      call array2_free(hov)
      call array2_free(hvo)
      call array2_free(hvv)
      if(whattodo==1) then  ! Residual-T1 transformed integrals
         call array2_free(xoR)
         call array2_free(xvR)
         call array2_free(yoR)
         call array2_free(yvR)
      end if


!!$         ! Fock matrix from T1-transformed density (not needed now but keep it commented out)
!!$
!!$         ! T1-transformed density: 
!!$         ! Dt1(rho,sigma) = sum_i Y_{rho i} X_{sigma i}
!!$         ! (some might say that this is the transposed of the T1-transposed density matrix,
!!$         !  however, this convention is chosen in accordance with the way Fock matrices
!!$         !  are built for nonsymmetric density matrices)
!!$         Dt1 = array2_init(bb)
!!$         call array2_matmul(yo,xo,Dt1,'N','T',1.0_realk,0.0_realk)
!!$         ! T1-transformed Fock matrix in AO basis
!!$         fao = array2_init(bb)
!!$         call dec_fock_transformation_fortran_array(nbasis,fao%val,Dt1%val,MyLsitem,.false.,incl_h=.true.)
!!$
!!$         ! T1-transformed Fock matrix in MO basis partioned into the four blocks
!!$         foo = array2_similarity_transformation(xo,fao,yo,oo)
!!$         fov = array2_similarity_transformation(xo,fao,yv,ov)
!!$         fvo = array2_similarity_transformation(xv,fao,yo,vo)
!!$         fvv = array2_similarity_transformation(xv,fao,yv,vv)
!!$
!!$         call array2_free(Dt1)
!!$         call array2_free(fao)


      ! Clean up
      call array2_free(xo)
      call array2_free(xv)
      call array2_free(yo)
      call array2_free(yv)
      call array2_free(R1)
      call array4_free(gvooo)

    end subroutine jacobian_rhtr_workhorse_integrals


    !> \brief Calculate two-electron trial-T1 transformed integrals, 
    !>        see Eq. 47 in JCP 105, 6921 (1996):
    !>
    !> g(p,q,r,s) = g1(p',q,r,s) + g2(p,q',r,s) + g3(p,q,r',s) + g4(p,q,r,s') 
    !>
    !> where primed index "n" is transformed with the znR transformation matrix,
    !> while the nonprimed indices are transformed with the z1,z2,z3,z4 input matrices.
    !>
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine two_electron_trial_T1_integrals(gao,z1,z2,z3,z4,z1R,z2R,z3R,z4R,g)

      implicit none
      !> Two-electron integrals in AO basis (not modified although intent(inout))
      type(array4),intent(inout) :: gao
      !> Transformation matrices as described above
      type(array2),intent(in) :: z1,z2,z3,z4,z1R,z2R,z3R,z4R
      !> Output two-electron integrals as described above
      !> Note: Also initialized here!
      type(array4),intent(inout) :: g
      type(array4) :: g1,g2,g3,g4
      integer :: i,j,k,l
      logical :: somethingwrong


      ! Sanity check
      somethingwrong=.false.
      if(z1%dims(1)/=z1R%dims(1) .or. z1%dims(2)/=z1R%dims(2)) somethingwrong=.true.
      if(z2%dims(1)/=z2R%dims(1) .or. z2%dims(2)/=z2R%dims(2)) somethingwrong=.true.
      if(z3%dims(1)/=z3R%dims(1) .or. z3%dims(2)/=z3R%dims(2)) somethingwrong=.true.
      if(z4%dims(1)/=z4R%dims(1) .or. z4%dims(2)/=z4R%dims(2)) somethingwrong=.true.
      if(somethingwrong) then
         call lsquit('two_electron_trial_T1_integrals: dimension mismatch!',-1)
      end if

      ! Calculate g1,g2,g3,g4 as described above
      g1 = get_gmo_simple(gao,z1R,z2,z3,z4)
      g2 = get_gmo_simple(gao,z1,z2R,z3,z4)
      g3 = get_gmo_simple(gao,z1,z2,z3R,z4)
      g4 = get_gmo_simple(gao,z1,z2,z3,z4R)

      ! Add up contributions to construct final two-electron integrals
      g = array4_init([z1%dims(2),z2%dims(2),z3%dims(2),z4%dims(2)])
      do l=1,z4%dims(2)
         do k=1,z3%dims(2)
            do j=1,z2%dims(2)
               do i=1,z1%dims(2)
                  g%val(i,j,k,l) = g1%val(i,j,k,l) + g2%val(i,j,k,l) + &
                       & g3%val(i,j,k,l) + g4%val(i,j,k,l)
               end do
            end do
         end do
      end do
      call array4_free(g1)
      call array4_free(g2)
      call array4_free(g3)
      call array4_free(g4)

    end subroutine two_electron_trial_T1_integrals



    !> \brief Calculate one-electron trial-T1 transformed integrals, 
    !>        see Eq. 46 in JCP 105, 6921 (1996):
    !>
    !> h(p,q) = h1(p',q) + h2(p,q') 
    !>
    !> where primed index "n" is transformed with the znR transformation matrix,
    !> while the nonprimed indices are transformed with the z1 and z2 input matrices.
    !>
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine one_electron_trial_T1_integrals(hao,z1,z2,z1R,z2R,h)

      implicit none
      !> One-electron integrals in AO basis (not modified although intent(inout))
      type(array2),intent(inout) :: hao
      !> Transformation matrices as described above (not modified although intent(inout))
      type(array2),intent(inout) :: z1,z2,z1R,z2R
      !> Output one-electron integrals as described above
      !> Note: Also initialized here!
      type(array2),intent(inout) :: h
      type(array2) :: h1,h2
      integer :: i,j
      logical :: somethingwrong


      ! Sanity check
      somethingwrong=.false.
      if(z1%dims(1)/=z1R%dims(1) .or. z1%dims(2)/=z1R%dims(2)) somethingwrong=.true.
      if(z2%dims(1)/=z2R%dims(1) .or. z2%dims(2)/=z2R%dims(2)) somethingwrong=.true.
      if(somethingwrong) then
         call lsquit('one_electron_trial_T1_integrals: dimension mismatch!',-1)
      end if

      ! Calculate h1 and h2 as described above
      h1 = array2_similarity_transformation(z1R,hao,z2,[z1%dims(2),z2%dims(2)])
      h2 = array2_similarity_transformation(z1,hao,z2R,[z1%dims(2),z2%dims(2)])

      ! Add up contributions to construct final one-electron integrals
      h = array2_add(h1,h2)

      call array2_free(h1)
      call array2_free(h2)

    end subroutine one_electron_trial_T1_integrals


    !> \brief Noddy code implementation of Davidson eigenvalue solver (J. Comp. Phys. 17, 87 (1975))
    !> for Jacobian right-hand side eigenvalue problem to get CCSD excitation energies.
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine ccsd_eigenvalue_solver(nbasis,nocc,nvirt,fAO,Co,Cv,mylsitem,t1,t2)
      implicit none
      !> Number of atomic/occupied/virtual orbitals
      integer,intent(in) :: nbasis,nocc,nvirt
      !> Fock matrix in AO basis
      real(realk),intent(in) :: fAO(nbasis,nbasis)
      !> LS item structure
      type(lsitem), intent(inout) :: MyLsItem
      !> Occupied and virtual MO coefficients
      !> (in practice these are intent(in))
      real(realk),intent(inout) :: Co(nbasis,nocc),Cv(nbasis,nvirt)
      !> Singles and doubles amplitudes
      !> (effectively intent(in) but need to be intent(inout) for practical purposes)
      type(tensor),intent(inout) :: t1,t2
      real(realk),pointer :: lambdaREAL(:),Arbig(:,:),alphaR(:,:),alphaL(:,:),Ar(:,:),alpha(:,:)
      integer :: Mold,k,M,p,i,j,iter
      integer(kind=long) :: maxnumeival,O,V
      type(tensor),pointer :: b1(:), b2(:),Ab1(:),Ab2(:)
      type(tensor) :: q1,zeta1,q2,zeta2,bopt1,bopt2
      type(tensor) :: xo,xv,yo,yv,t1fock,Co_tensor,Cv_tensor
      real(realk) :: tmp1,tmp2,fac,bTzeta,bnorm,sc
      real(realk),pointer :: res(:),lambdaIMAG(:),eival(:),singlescomp(:),foo(:,:),fvv(:,:)
      logical,pointer :: conv(:),singleex(:)
      logical :: allconv
      
      ! Sanity check
      if(DECinfo%haldapprox .and. DECinfo%JacobianLhtr) then
         call lsquit('Hald approximation not implemented for &
              & Jacobian left-hand transformation!',-1)
      end if

      ! Notation and implementation is exactly like p. 91 of J. Comp. Phys. 17, 87 (1975), 
      ! except that we may solve for more than one eigenvalue at the same time.
      ! Also the convergence check is different and simply based on the residual.

      ! How many eigenvalues k?
      k=DECinfo%JacobianNumEival

      ! Size of zero-order orthonormal subspace M (M must be equal to or larger than k)
      M = DECinfo%JacobianInitialSubspace
      if(M==0) then
         ! M was not set by input, set it equal to the number of requested excitation energies
         M=k
      else
         ! M set by input, only modify it if it is too small
         if(M < k) then
            write(DECinfo%output,*) 'WARNING! Initial Jacobian subspace set to', &
                 & M
            write(DECinfo%output,*) 'WARNING! This is smaller than the number of excitation energies',k
            write(DECinfo%output,*) 'WARNING! I overrule input and set initial Jacobian subspace to',k
            M=k
         end if
      end if


      ! Sanity check: Requested number of initial eigenvalues/start vectors is not
      ! larger than the maximum possible number of singles+doubles excitations.
      O = int(nocc,kind=8)
      V = int(nvirt,kind=8)
      ! Different unique exciations (i/=j and a/=b):
      ! Singles: t_i^a  --> O*V
      ! Doubles: t_ii^aa --> O*V
      !          t_ij^aa --> O*(O-1)*V / 2
      !          t_ii^ab --> O*V*(V-1) / 2
      !          t_ij^ab --> O*(O-1)*V*(V-1) / 2
      !          
      ! Division by 2 is due to symmetry: t_ij^ab = t_ji^ba
      ! 
      maxnumeival = O*V + O*V + O*(O-1)*V/2 + O*V*(V-1)/2 + O*(O-1)*V*(V-1)/2

      if( int(M,kind=8) > maxnumeival) then
         print *, 'Number of startvectors ', M
         print *, 'Max number of eigenvalues ',maxnumeival
         call lsquit('ccsd_eigenvalue_solver: Number of requested &
              & start vectors is too large!',-1)
      end if


      ! ********************************************************************************
      !                         INITIALIZATION BEFORE SOLVER LOOP                      !
      ! ********************************************************************************

      ! Occ-occ and virt-virt Fock matrix blocks
      call mem_alloc(foo,nocc,nocc)
      call mem_alloc(fvv,nvirt,nvirt)
      call dec_simple_basis_transform1(nbasis,nocc,Co,fAO,foo)
      call dec_simple_basis_transform1(nbasis,nvirt,Cv,fAO,fvv)

      ! Quick workaround. If the model does not use singles, we simply initialize the singles
      ! are and set the singles to zero
      if(.not. DECinfo%use_singles) then
         call tensor_minit(t1,[nvirt,nocc],2)
         call tensor_zero(t1)
      end if

      ! T1 transformed quantities
      call get_T1_transformed_for_eigenvalue_solver(nbasis,nocc,nvirt,mylsitem,Co,Cv,fAO,t1,&
           & xo,xv,yo,yv,t1fock)

      ! It's a mess but we need MO coefficients also in tensor format
      call tensor_minit(Co_tensor,[nbasis,nocc],2)
      call tensor_minit(Cv_tensor,[nbasis,nvirt],2)
      call tensor_convert(Co,Co_tensor)
      call tensor_convert(Cv,Cv_tensor)


      ! Allocate space for residual-related stuff
      ! -----------------------------------------
      call tensor_minit(q1,[nvirt,nocc],2)
      call tensor_minit(zeta1,[nvirt,nocc],2)
      call tensor_minit(q2,[nvirt,nocc,nvirt,nocc],4)
      call tensor_minit(zeta2,[nvirt,nocc,nvirt,nocc],4)
      call tensor_minit(bopt1,[nvirt,nocc],2)
      call tensor_minit(bopt2,[nvirt,nocc,nvirt,nocc],4)
      call mem_alloc(res,k)
      call mem_alloc(conv,k)
      call mem_alloc(singleex,M)
      call mem_alloc(eival,k)
      call mem_alloc(singlescomp,k)
      conv=.false.
      allconv=.false.


      ! Get start vectors and associated eigenvalues
      ! --------------------------------------------
      call mem_alloc(b1,DECinfo%JacobianMaxSubspace)
      call mem_alloc(b2,DECinfo%JacobianMaxSubspace)
      call mem_alloc(Ab1,DECinfo%JacobianMaxSubspace)
      call mem_alloc(Ab2,DECinfo%JacobianMaxSubspace)
      do i=1,M
         call tensor_minit(b1(i),[nvirt,nocc],2)
         call tensor_minit(b2(i),[nvirt,nocc,nvirt,nocc],4)
         call tensor_minit(Ab1(i),[nvirt,nocc],2)
         call tensor_minit(Ab2(i),[nvirt,nocc,nvirt,nocc],4)
      end do
      call mem_alloc(lambdaREAL,M)
      call ccsd_eigenvalue_solver_startguess(nbasis,mylsitem,xo,xv,&
           & yo,yv,M,nocc,nvirt,foo,fvv,b1(1:M),b2(1:M),lambdaREAL,singleex)
      eival(1:k) = lambdaREAL(1:k)


      ! Form Jacobian right-hand transformations A b on initial M trial vectors b
      ! -------------------------------------------------------------------------
      do i=1,M
         if(DECinfo%JacobianLhtr) then
            call cc_jacobian_lhtr(mylsitem,Co_tensor,Cv_tensor,xo,xv,yo,yv,t1fock,&
                 & t1,t2,b1(i),b2(i),Ab1(i),Ab2(i))
         else
            call cc_jacobian_rhtr(mylsitem,xo,xv,yo,yv,t1,t2,b1(i),b2(i),Ab1(i),Ab2(i))
         end if
      end do


      ! Form Jacobian in reduced space
      ! -------------------------------
      ! Note: Arbig always contains the components of the Jacobian constructed at the given time,
      !       and it is allocated such that there is room for more elements.
      !       Ar below will contain the same elements but with the current dimensions.

      call mem_alloc(Arbig,DECinfo%JacobianMaxSubspace,DECinfo%JacobianMaxSubspace)
      Arbig=0.0_realk
      do j=1,M
         do i=1,M
            ! Note that b(i) and Ab(i) both have structures (singles, doubles).
            ! The reduced matrix is the dotproduct <b(i),Ab(j)>,
            ! and it therefore becomes a sum of dot products of the singles
            ! and doubles components
            if(DECinfo%JacobianLhtr) then
               ! For left-hand transformation, Ab is really b^T A where
               ! b^T is a row vector and A is a matrix. We therefore effectively
               ! get entry (j,i) of the reduced Jacobian when we calculate b(i) "dot" Ab(j) which is
               ! really equal to {b^T A}(j) b(i).
               ! (***)
               Arbig(j,i) = tensor_ddot(b1(i),Ab1(j)) + tensor_ddot(b2(i),Ab2(j))
            else
               Arbig(i,j) = tensor_ddot(b1(i),Ab1(j)) + tensor_ddot(b2(i),Ab2(j))
            end if
         end do
      end do



      ! ***************************************************************************************
      ! *                                MAIN SOLVER ITERATIONS                               *
      ! ***************************************************************************************

      write(DECinfo%output,'(1X,a)') 'JAC ********************************************************'
      write(DECinfo%output,'(1X,a)') 'JAC        Information for Jacobian eigenvalue solver       '
      write(DECinfo%output,'(1X,a)') 'JAC        ------------------------------------------       '
      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,*) 'JAC Number of eigenvalues      ',k
      write(DECinfo%output,*) 'JAC Initial subspace dimension ',M
      write(DECinfo%output,*) 'JAC Maximum subspace dimension ',maxnumeival
      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,'(1X,a)') 'JAC Start guess for eigenvalues'
      write(DECinfo%output,'(1X,a)') 'JAC ---------------------------'
      do i=1,k
         if(singleex(i)) then
            write(DECinfo%output,'(1X,a,i7,g20.10,1X,a)') 'JAC ',i,lambdaREAL(i), '  (singles ex.)'
            singlescomp(i) = 100.0_realk
         else
            write(DECinfo%output,'(1X,a,i7,g20.10,1X,a)') 'JAC ',i,lambdaREAL(i), '  (doubles ex.)'
            singlescomp(i) = 0.0_realk
         end if
      end do
      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,'(1X,a)') 'JAC ********************************************************'

      call mem_dealloc(lambdaREAL)

      write(DECinfo%output,'(1X,a)') 'JAC'
      write(DECinfo%output,'(1X,a)') 'JAC Jacobian eigenvalue solver'
      write(DECinfo%output,'(1X,a)') 'JAC'
      write(DECinfo%output,'(1X,a)') 'JAC Which   Subspace     Eigenvalue           Residual       Conv?  '

      SolverLoop: do iter=1,DECinfo%JacobianMaxIter

         ! Copy elements of reduced Jacobian into array having the proper dimensions
         ! -------------------------------------------------------------------------
         call mem_alloc(Ar,M,M)
         do j=1,M
            do i=1,M
               Ar(i,j) = Arbig(i,j)
            end do
         end do


         ! Diagonalize reduced matrix Ar
         ! -----------------------------
         ! A alphaR = lambda alphaR
         ! alphaL A = alphaL lambda
         ! lambda is a diagonal matrix with eigenvalues on the diagonal
         ! 
         ! Note: We only consider real part of lambda! Imaginary part is simply removed
         ! but a warning is issued if the imaginary part is sizable.
         call mem_alloc(alphaR,M,M)
         call mem_alloc(alphaL,M,M)
         call mem_alloc(lambdaREAL,M)
         call mem_alloc(lambdaIMAG,M)
         call solve_nonsymmetric_eigenvalue_problem_unitoverlap(M,Ar,lambdaREAL,&
              & lambdaIMAG,alphaR,alphaL)
         eival(1:k) = lambdaREAL(1:k)
         call mem_dealloc(Ar)

         ! Any imaginary eigenvalue component? Print a warning if this is the case.
         call ccsd_eigenvalue_check(k,M,lambdaREAL,lambdaIMAG)


         ! NOTE: In general, for the Jacobian A we have:
         !
         ! A R = lambda A R 
         ! L A = lambda L A
         !
         ! where A is a matrix, R/L is the right/left column/row eigenvector for eigenvalue lambda. 
         ! For left/right eigenvectors, we need to work with left/right eigenvalues.
         call mem_alloc(alpha,M,M)
         do j=1,M
            do i=1,M
               if(DECinfo%JacobianLhtr) then
                  alpha(i,j) = alphaL(i,j)
               else
                  alpha(i,j) = alphaR(i,j)
               end if
            end do
         end do
         call mem_dealloc(alphaR)
         call mem_dealloc(alphaL)


         ! Save current subspace dimension
         Mold = M

         ! Residual, preconditioning and new trial vectors
         ! -----------------------------------------------
         ! Note. Here we check for k residuals (p=1,k), while in J. Comp. Phys. 17, 87 (1975)
         ! only the k'th eigenvalue is considered.
         ExcitationEnergyLoop: do p=1,k


            ! Current optimal eigenvector bopt for eigenvalue p:
            ! ''''''''''''''''''''''''''''''''''''''''''''''''''
            ! bopt = sum_i^M alpha(i,p) b(i) 
            ! 
            ! Residual q:
            ! '''''''''''
            ! q = A bopt - eival bopt
            ! 
            ! Residual contains two components: singles (q1) and doubles (q2)
            ! 
            call tensor_zero(q1)
            call tensor_zero(q2)
            call tensor_zero(bopt1)
            call tensor_zero(bopt2)

            do i=1,Mold
               ! Update A bopt part
               call tensor_add(q1,alpha(i,p),Ab1(i))
               call tensor_add(q2,alpha(i,p),Ab2(i))

               ! Calculate bopt
               call tensor_add(bopt1,alpha(i,p),b1(i))
               call tensor_add(bopt2,alpha(i,p),b2(i))
            end do


            ! Residual: A bopt - eival bopt  
            fac = - lambdaREAL(p)
            call tensor_add(q1,fac,bopt1)
            call tensor_add(q2,fac,bopt2)

            ! Singles component of optimal vector: |bopt1| / |bopt|
            call print_norm(bopt1,nrm=tmp1,returnsquared=.true.)
            tmp2 = SD_dotproduct(bopt1,bopt2)
            singlescomp(p) = (sqrt(tmp1) / sqrt(tmp2) )*100.0_realk

            ! Residual norm
            res(p) = SD_dotproduct(q1,q2)
            res(p) = sqrt(res(p))

            ! Check for convergence
            conv(p) = (res(p)<DECinfo%JacobianThr) 
            write(DECinfo%output,'(1X,a,i6,2X,i6,4X,g18.8,1X,g18.8,3X,L2)') &
                 & 'JAC',p,M,lambdaREAL(p),res(p),conv(p)

            if(conv(p)) then
               ! Do not include new vectors for eigenvalue p if it is already converged
               cycle ExcitationEnergyLoop
            else
               ! Update M
               M=M+1
            end if

            ! Sanity check for dimensions
            if(M > DECinfo%JacobianMaxSubspace) then
               call lsquit('ccsd_eigenvalue_solver: M is too large! Not possible to extend subspace!',-1)
            end if

            ! Preconditioning
            if(DECinfo%JacobianPrecond) then
               call precondition_jacobian_residual(nocc,nvirt,foo,fvv,lambdaREAL(p),q1,q2,zeta1,zeta2)
            else
               ! Simply copy elements
               zeta1%elm2 = q1%elm2
               zeta2%elm4 = q2%elm4
            end if

            ! Component of precondioned residual orthogonal to current trial vectors:
            !
            ! b(M) = [ 1 - sum_{i=1}^{M-1} b(i) b(i)^T ] zeta
            !
            call tensor_minit(b1(M),[nvirt,nocc],2)
            call tensor_minit(b2(M),[nvirt,nocc,nvirt,nocc],4)
            call tensor_zero(b1(M))
            call tensor_zero(b2(M))
            ! b(M+p) = zeta
            call tensor_add(b1(M),1.0_realk,zeta1)
            call tensor_add(b2(M),1.0_realk,zeta2)

            do i=1,M-1
               ! b(i)^T zeta  (just a number when b(i) and zeta are considered as vectors)
               bTzeta = - ( tensor_ddot(b1(i),zeta1) + tensor_ddot(b2(i),zeta2) )
               ! b(M+p) += - b(i) { b(i)^T zeta }
               call tensor_add(b1(M),bTzeta,b1(i))
               call tensor_add(b2(M),bTzeta,b2(i))
            end do

            ! Normalize b(M)
            bnorm = SD_dotproduct(b1(M),b2(M))
            bnorm = sqrt(bnorm)
            if(bnorm < 1.0e-15_realk) then
               ! Retreat! We do not want to include the current b vector anyways.
               ! The components spanned by the current b vector were 
               ! probably included by the b vector of one of the other excitation energies.
               ! Reset M and deallocate b1 and b2 again.
               call tensor_free(b1(M))
               call tensor_free(b2(M))
               M = M-1
               write(DECinfo%output,'(1X,a,2i7)') &
                    & 'WARNING: Zero norm after projection for exci/iter', p,iter
               cycle ExcitationEnergyLoop
            end if

            sc = 1.0_realk/bnorm
            call tensor_scale(b1(M),sc)
            call tensor_scale(b2(M),sc)

            ! Form Ab(M)
            call tensor_minit(Ab1(M),[nvirt,nocc],2)
            call tensor_minit(Ab2(M),[nvirt,nocc,nvirt,nocc],4)
            if(DECinfo%JacobianLhtr) then
               call cc_jacobian_lhtr(mylsitem,Co_tensor,Cv_tensor,xo,xv,yo,yv,t1fock,&
                    & t1,t2,b1(M),b2(M),Ab1(M),Ab2(M))
            else
               call cc_jacobian_rhtr(mylsitem,xo,xv,yo,&
                    & yv,t1,t2,b1(M),b2(M),Ab1(M),Ab2(M))
            end if

         end do ExcitationEnergyLoop

         ! Free stuff
         call mem_dealloc(alpha)
         call mem_dealloc(lambdaREAL)
         call mem_dealloc(lambdaIMAG)

         ! Are all eigenvalues converged?
         allconv = all(conv)
         if(allconv) then
            exit SolverLoop
         end if

         ! Calculate new blocks of reduced Jacobian matrix A with new M. A is:
         !
         !      |   A(1:Mold ,   1:Mold)       A(1:Mold ,   Mold+1:M)    |
         ! A =  |   A(Mold+1:M , 1:Mold)       A(Mold+1:M , Mold+1:M)    |
         ! 
         ! A(1:Mold , 1:Mold) has already been calculated, while the other blocks are
         ! calculated here.
         ! 
         ! Difference between right- and left transformations, see comment (***) above!
         do i=1,M
            do j=Mold+1,M
               if(DECinfo%JacobianLhtr) then
                  Arbig(j,i) = tensor_ddot(b1(i),Ab1(j)) + tensor_ddot(b2(i),Ab2(j))
                  if(i/=j) Arbig(i,j) = tensor_ddot(b1(j),Ab1(i)) + tensor_ddot(b2(j),Ab2(i))
               else
                  Arbig(i,j) = tensor_ddot(b1(i),Ab1(j)) + tensor_ddot(b2(i),Ab2(j))
                  if(i/=j) Arbig(j,i) = tensor_ddot(b1(j),Ab1(i)) + tensor_ddot(b2(j),Ab2(i))
               end if
            end do
         end do

      end do SolverLoop

      if(.not. allconv) then
         call lsquit('Jacobian eigenvalue solver did not converge!',-1)
      end if


      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,'(1X,a)') 'JAC ******************************************************************************'
      write(DECinfo%output,'(1X,a)') 'JAC CC excitation energies'
      write(DECinfo%output,'(1X,a)') 'JAC ******************************************************************************'
      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,'(1X,a)') 'JAC      Exci.    Hartree           eV            cm-1           ||T1|| (%)'
      do p=1,k
         write(DECinfo%output,'(1X,a,i7,4X,3g15.7,F10.2)') 'JAC ',p,eival(p),eival(p)*hartree_to_eV,&
              & eival(p)*hartree_to_cm1, singlescomp(p)
      end do
      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,'(1X,a)') 'JAC ******************************************************************************'
      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,'(1X,a)') 'JAC '
      write(DECinfo%output,'(1X,a)') 'JAC '


      do i=1,M
         call tensor_free(b1(i))
         call tensor_free(b2(i))
         call tensor_free(Ab1(i))
         call tensor_free(Ab2(i))
      end do
      call mem_dealloc(b1)
      call mem_dealloc(b2)
      call mem_dealloc(Ab1)
      call mem_dealloc(Ab2)

      call mem_dealloc(foo)
      call mem_dealloc(fvv)
      call mem_dealloc(eival)
      call mem_dealloc(singlescomp)
      call mem_dealloc(singleex)
      call mem_dealloc(res)
      call mem_dealloc(conv)
      call mem_dealloc(Arbig)
      call tensor_free(zeta1)
      call tensor_free(zeta2)
      call tensor_free(q1)
      call tensor_free(q2)
      call tensor_free(bopt1)
      call tensor_free(bopt2)
      call tensor_free(xo)
      call tensor_free(xv)
      call tensor_free(yo)
      call tensor_free(yv)
      call tensor_free(t1fock)
      call tensor_free(Co_tensor)
      call tensor_free(Cv_tensor)
      if(.not. DECinfo%use_singles) then
         call tensor_free(t1)
      end if

    end subroutine ccsd_eigenvalue_solver



    !> Precondition singles and doubles block for Jacobian trial vectors.
    !> This corresponds to step D of J. Comp. Phys. 17, 87 (1975) where we simply use
    !> Fock matrix differences to represent Jacobian diagonal.
    subroutine precondition_jacobian_residual(nocc,nvirt,foo,fvv,lambda,q1,q2,zeta1,zeta2)
      implicit none
      !> Number of occupied/virtual orbitals in molecule
      integer,intent(in) :: nocc,nvirt
      !> Occ-occ and virt-virt Fock matrix blocks in MO basis
      real(realk),intent(in) ::  foo(nocc,nocc),fvv(nvirt,nvirt)
      !> Eigenvalue lambda
      real(realk),intent(in) :: lambda
      !> Residual before preconditioning (singles "1" and doubles "2" components)
      type(tensor),intent(in) :: q1,q2
      !> Preconditioned residual (singles "1" and doubles "2" components)
      type(tensor),intent(inout) :: zeta1,zeta2
      integer :: i,j,a,b
      real(realk) :: Aelm

      ! Singles preconditioning
      do i=1,nocc
         do a=1,nvirt
            Aelm = fvv(a,a) - foo(i,i)
            if(abs(Aelm-lambda)<1.0e-10) then  
               print *, 'a,i,A,lambda',a,i,Aelm,lambda
               call lsquit('singles: unstable Jacobian residual',-1) 
            end if
            zeta1%elm2(a,i) = ( 1.0_realk/(lambda-Aelm) ) * q1%elm2(a,i)
         end do
      end do

      ! Doubles preconditioning
      do j=1,nocc
         do b=1,nvirt
            do i=1,nocc
               do a=1,nvirt
                  Aelm = fvv(a,a) + fvv(b,b) - foo(i,i) - foo(j,j)
                  if(abs(Aelm-lambda)<1.0e-10) then 
                     print *, 'a,i,b,j,A,lambda',a,i,b,j,Aelm,lambda
                     call lsquit('doubles: unstable Jacobian residual',-1) 
                  end if
                  zeta2%elm4(a,i,b,j) = ( 1.0_realk/(lambda-Aelm) ) * q2%elm4(a,i,b,j)
               end do
            end do
         end do
      end do

    end subroutine precondition_jacobian_residual


    !> \brief Calculate start vector and associated eigenvalues for CCSD eigenvalue solver.
    !> We simply take the lowest M orbital energies differences in singles and doubles spaces. Correction - also add dominating two-electron integrals.
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine ccsd_eigenvalue_solver_startguess(nbasis,mylsitem,xo_tensor,xv_tensor,yo_tensor,yv_tensor,&
         & M,nocc,nvirt,foo,fvv,b1,b2,eival,singleex)
      implicit none
      integer,intent(in) :: nbasis
      type(lsitem) :: mylsitem
      type(tensor) :: xo_tensor,xv_tensor,yo_tensor,yv_tensor
      !> Number of start vectors requested
      integer,intent(in) :: M
      !> Number of occupied/virtual orbitals in molecule
      integer,intent(in) :: nocc,nvirt
      !> Occ-occ and virt-virt Fock matrix blocks in MO basis
      real(realk),intent(in) ::  foo(nocc,nocc),fvv(nvirt,nvirt)
      !> Singles and Doubles components of start vectors
      type(tensor),intent(inout) :: b1(M), b2(M)
      !> Eigenvalues (orbital energy differences) for start vectors
      real(realk),intent(inout) :: eival(M)
      !> Is the initial eigenvalue a single excitation (true) or double excitation (false)
      logical,intent(inout) :: singleex(M)
      integer :: i,j,a,b,maxidx(1),MS,MD,idxS,idxD
      real(realk) :: ex,maxeival,eivalS(M),eivalD(M),mineivalD,mineivalS
      integer,dimension(M) :: s1,s2,d1,d2,d3,d4,sd
      logical,dimension(M) :: inclS, inclD
      type(array4) :: gao,gvoov,gvvoo
      type(array2) :: xo,xv,yo,yv

      singleex=.true.

      xo = array2_init([nbasis,nocc],xo_tensor%elm2)
      xv = array2_init([nbasis,nvirt],xv_tensor%elm2)
      yo = array2_init([nbasis,nocc],yo_tensor%elm2)
      yv = array2_init([nbasis,nvirt],yv_tensor%elm2)



      ! Calculate full AO integrals
      call get_full_eri(mylsitem,nbasis,gao)
      gvoov = get_gmo_simple(gao,xv,yo,xo,yv)
      gvvoo = get_gmo_simple(gao,xv,yv,xo,yo)
      call array4_free(gao)



      ! M lowest singles excitation energies
      ! ************************************
      eivalS = huge(1.0)
      maxeival = huge(1.0)
      maxidx=1
      s1=0; s2=0
      do i=1,nocc
         do a=1,nvirt

            ! Zero-order excitation energy (orbital difference)
            ex = fvv(a,a) - foo(i,i) + 2.0_realk*gvoov%val(a,i,i,a) - gvvoo%val(a,a,i,i)

            ! Is excitation energy lower than current maximum eigenvalue?
            if(ex < maxeival) then
               ! Replace largest excitation energy by excitation energy under consideration
               ! and also save indices
               eivalS( maxidx(1) ) = ex  
               s1( maxidx(1) ) = a
               s2( maxidx(1) ) = i

               ! New maximum eigenvalue and position of it
               maxeival = maxval(eivalS)
               maxidx = maxloc(eivalS)
            end if

         end do
      end do

      ! M lowest doubles excitation energies: Same strategy as for singles
      ! ******************************************************************
      eivalD = huge(1.0)
      maxeival = huge(1.0)
      maxidx=1
      d1=0; d2=0; d3=0; d4=0
      iloop: do i=1,nocc
         jloop: do j=1,nocc
            aloop: do a=1,nvirt
               bloop: do b=a,nvirt

                  ! avoid taking (a,i,a,j) combination twice
                  if(a==b .and. i>j) cycle bloop  


                  ex = fvv(a,a) + fvv(b,b) - foo(i,i) - foo(j,j) &
                       & + 2.0_realk*gvoov%val(a,i,i,a) - gvvoo%val(a,a,i,i) &
                       & + 2.0_realk*gvoov%val(b,j,j,a) - gvvoo%val(b,b,j,j) 

                  if(ex < maxeival) then
                     eivalD( maxidx(1) ) = ex
                     d1( maxidx(1) ) = a
                     d2( maxidx(1) ) = i
                     d3( maxidx(1) ) = b
                     d4( maxidx(1) ) = j

                     maxeival = maxval(eivalD)
                     maxidx = maxloc(eivalD)
                  end if

               end do bloop
            end do aloop
         end do jloop
      end do iloop


      ! M lowest excitation energies - singles OR doubles
      ! *************************************************
      inclS = .false.
      inclD = .false.

      do i=1,M  
         ! Find i'th lowest eigenvalue, regardless of whether it's singles or doubles
         mineivalS = huge(1.0)
         mineivalD = huge(1.0)


         do j=1,M

            ! Find minimum singles eigenvalue and corresponding index not already included
            if(.not. inclS(j)) then
               if(eivalS(j) < mineivalS) then
                  idxS = j   
                  mineivalS = eivalS(j)
               end if
            end if

            ! Find minimum doubles eigenvalue and corresponding index not already included
            if(.not. inclD(j)) then
               if(eivalD(j) < mineivalD) then
                  idxD = j
                  mineivalD = eivalD(j)
               end if
            end if

         end do

         ! Initialize singles and doubles start vectors
         call tensor_zero(b1(i))
         call tensor_zero(b2(i))

         ! Set i'th eigenvalue
         ! -------------------
         if(mineivalS < mineivalD) then ! the singles eigenvalue is the lowest

            ! this eigenvalue is now included and cannot be considered again
            inclS(idxS) = .true.  
            eival(i) = mineivalS
            singleex(i) = .true.

            ! Singles vector - 1 in position of eigenvalue, zero elsewhere
            ! Doubles vector - zero 
            b1(i)%elm2(s1(idxS),s2(idxS)) = 1.0_realk

         else ! the doubles eigenvalue is the lowest

            ! this eigenvalue is now included and cannot be considered again
            inclD(idxD) = .true.
            eival(i) = mineivalD
            singleex(i) = .false.

            ! Singles vector - zero
            ! Doubles vector - 1 in position of eigenvalue, zero elsewhere
            ! However, since doubles vector must have symmetry:
            ! X_bjai = X_aibj
            ! we need to set X_bjai = X_aibj = 1/sqrt(2)
            ! unless a=b and i=j
            if(d1(idxD)==d3(idxD) .and. d2(idxD)==d4(idxD)) then
               b2(i)%elm4(d1(idxD),d2(idxD),d3(idxD),d4(idxD)) = 1.0_realk
            else
               b2(i)%elm4(d1(idxD),d2(idxD),d3(idxD),d4(idxD)) = 1.0_realk/sqrt(2.0_realk)
               b2(i)%elm4(d3(idxD),d4(idxD),d1(idxD),d2(idxD)) = 1.0_realk/sqrt(2.0_realk)
            end if

         end if
         

      end do

      call array2_free(xo)
      call array2_free(xv)
      call array2_free(yo)
      call array2_free(yv)
      call array4_free(gvoov)
      call array4_free(gvvoo)

    end subroutine ccsd_eigenvalue_solver_startguess


    !> The "reality" of CCSD eigenvalues and print warnings if any eigenvalue contains a significant
    !> imaginary component.
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine ccsd_eigenvalue_check(k,M,lambdaREAL,lambdaIMAG)
      implicit none
      !> Number of target excitation energies
      integer,intent(in) :: k
      !> Dimension of current subspace
      integer,intent(in) :: M
      !> Real/imaginary components of the M eigenvalues
      real(realk),dimension(M),intent(in) :: lambdaREAL, lambdaIMAG
      !> If the imaginary eigenvalue component is above this threshold we print a warning
      real(realk),parameter :: eivalthr=1.0e-9
      integer :: i

      ! Target eigenvalues
      do i=1,k
         if( abs(lambdaIMAG(i))>eivalthr) then
            write(DECinfo%output,*) 'WARNING! Target eigenvalue contains imaginary component!'
            write(DECinfo%output,*) 'WARNING! We ignore the imaginary component and proceed.'
            write(DECinfo%output,*) 'WARNING! Be very careful when interpreting results!'
            write(DECinfo%output,'(1X,a,i7)') 'Eigenvalue number:     ',i
            write(DECinfo%output,'(1X,a,2g20.10)') 'Real/imag components:  ',&
                 & lambdaREAL(i),lambdaIMAG(i)
         end if
      end do

      ! Eigenvalues also present in subspace but which are not a target of calculation
      do i=k+1,M
         if( abs(lambdaIMAG(i))>eivalthr) then
            write(DECinfo%output,*) 'WARNING! Eigenvalue in subspace contains imaginary component!'
            write(DECinfo%output,*) 'WARNING! However, this is not a target eigenvalue'
            write(DECinfo%output,*) 'WARNING! so it is probably not a problem...'
            write(DECinfo%output,'(1X,a,i7)') 'Eigenvalue number:     ',i
            write(DECinfo%output,'(1X,a,2g20.10)') 'Real/imag components:  ',&
                 & lambdaREAL(i),lambdaIMAG(i)
         end if
      end do

    end subroutine ccsd_eigenvalue_check


    !> \brief Wrapper to calculate the necessary T1-transformed quantity used for
    !> CCSD Jacobian eigenvalue problem.
    !> Note: The tensors are also initialized here!
    !> \author Kasper Kristensen
    !> \date June 2015
    subroutine get_T1_transformed_for_eigenvalue_solver(nbasis,nocc,nvirt,mylsitem,Co,Cv,fAO,t1,&
         & xo,xv,yo,yv,t1fock)
      implicit none
      !> Number of atomic/occupied/virtual orbitals
      integer,intent(in) :: nbasis,nocc,nvirt
      !> Fock matrix in AO basis
      real(realk),intent(in) :: fAO(nbasis,nbasis)
      !> LS item structure
      type(lsitem), intent(inout) :: MyLsItem      
      !> Occ and virt MO coefficients. Note: Co=xo (particle) and Cv=yv (hole)
      !> (in practice these are intent(in))
      real(realk),intent(inout) :: Co(nbasis,nocc), Cv(nbasis,nvirt)
      !> Singles amplitudes  (in practice this is intent(in))
      type(tensor),intent(inout) :: t1
      !> Particle (x) and hole (y) matrices for occ and virt spaces (INITIALIZED HERE!)
      type(tensor),intent(inout) :: xo,xv,yo,yv
      !> T1 transformed Fock matrix in AO basis (INITIALIZED HERE!)
      type(tensor),intent(inout) :: t1fock
      type(tensor) :: fAO_tensor,Co_tensor,Cv_tensor

      ! Init tensors for hole and particle transformation coefficients + T1-transformed Fock
      call tensor_minit(xo,[nbasis,nocc],2)
      call tensor_minit(xv,[nbasis,nvirt],2)
      call tensor_minit(yo,[nbasis,nocc],2)
      call tensor_minit(yv,[nbasis,nvirt],2)
      call tensor_minit(t1fock,[nbasis,nbasis],2)

      ! Dirty workaround due to input in get_t1_matrices
      call tensor_minit(fAO_tensor,[nbasis,nbasis],2)
      call tensor_minit(Co_tensor,[nbasis,nocc],2)
      call tensor_minit(Cv_tensor,[nbasis,nvirt],2)
      call tensor_convert(fAO,fAO_tensor)
      call tensor_convert(Co,Co_tensor)
      call tensor_convert(Cv,Cv_tensor)

      ! Calculate T1-transformed Fock matrix
      call get_t1_matrices(MyLsitem,t1,Co_tensor,Cv_tensor,xo,yo,xv,yv,fAO_tensor,t1fock,.false.)

      call tensor_free(fAO_tensor)
      call tensor_free(Co_tensor)
      call tensor_free(Cv_tensor)

    end subroutine get_T1_transformed_for_eigenvalue_solver


    !> \brief Add contribution to Jacobian trial vector to get linear LW1 model 
    !> rather than exponential EW1 model.
    !> Note: For the Jacobian A only the doubles-singles block
    !> differs between LW1 and EW1.
    !> For right transformation, we get:
    !>
    !> (rho1) = ( A11 A12 )  (V1)     ( A11 V1 + A12 V2 )
    !> (rho2) = ( A21 A22 )  (V2)  =  ( A21 V1 + A22 V2 ) 
    !> 
    !> Here, only the A21 V1 contribution differs between EW1 and LW1. 
    !> The input rho2 is updated with this difference, i.e.:
    !>
    !> rho2 --> rho2 + (A21 V1)_{LW1} - (A21 V1)_{EW1}
    !> 
    !> rho1 is not modified for right transformation. 
    !> 
    !> For left transformation we only change V2 A21 contribution.
    !> Thus, rho1 is changed, while rho2 is not.
    !> 
    !> \author Kasper Kristensen
    !> \date July 2015
    subroutine add_lw1_contribution(nbasis,mylsitem,Co_tensor,&
         & Cv_tensor,t2,V1,V2,rho1,rho2,righttrans)
      implicit none
      !> Number of basis functions
      integer,intent(in) :: nbasis
      !> LSitem
      type(lsitem) :: mylsitem
      !> Occ and virt MO coefficients
      type(tensor) :: Co_tensor,Cv_tensor
      !> Doubles amplitudes
      type(tensor),intent(in) :: t2 
      !> Singles (V1) and doubles (V2) component of trial vector
      type(tensor),intent(in) :: V1,V2
      !> rho1 and rho2 vectors corresponding to Jacobian transf. on trial vectors
      type(tensor),intent(inout) :: rho1,rho2
      !> Right transformation (true) or left (false)
      logical,intent(in) :: righttrans
      integer :: nocc,nvirt,i,j,k,l,a,b,c,d
      real(realk),pointer :: Q(:,:),Y(:,:),deltaNS(:,:,:,:)
      real(realk) :: Lljkc, Lbdkc, Rai
      type(array4) :: gao,gooov,gvvov
      type(array2) :: Co,Cv

      ! Never do this for Hald approximation
      if(DECinfo%haldapprox) return

      ! Sanity check - LW1 only meaningful for MP2 model
      if(DECinfo%ccModel/=MODEL_MP2) then
         call lsquit('add_lw1_contribution: Only implemented for MP2 model!',-1)
      end if

      ! Dimensions
      nvirt = t2%dims(1)
      nocc = t2%dims(2)

      ! Calculate necessary MO integrals
      Co = array2_init([nbasis,nocc],Co_tensor%elm2)
      Cv = array2_init([nbasis,nvirt],Cv_tensor%elm2)
      call get_full_eri(mylsitem,nbasis,gao)
      gooov = get_gmo_simple(gao,Co,Co,Co,Cv)
      gvvov = get_gmo_simple(gao,Cv,Cv,Co,Cv)
      call array4_free(gao)


      ! For right-hand vector correction delta (LW1 - EW1) is:
      !
      ! delta_{aibj} = P_{ij}^{ab} [ - R_ai sum_{klc} t_kl^cb L_ljkc
      !                              + R_ai sum_{cdk} t_{kj}^{cd} L_{bdkc} ]
      !

      !  Q_bj = sum_{klc} t_kl^cb L_ljkc  
      call mem_alloc(Q,nvirt,nocc)
      Q = 0.0_realk
      do c=1,nvirt
         do k=1,nocc
            do j=1,nocc
               do l=1,nocc
                  Lljkc = 2.0_realk*gooov%val(l,j,k,c) - gooov%val(k,j,l,c)
                  do b=1,nvirt
                     Q(b,j) = Q(b,j) + t2%elm4(c,k,b,l) * Lljkc
                  end do
               end do
            end do
         end do
      end do

      !  Y_jb = sum_{cdk} t_{kj}^{cd} L_{bdkc}
      call mem_alloc(Y,nocc,nvirt)
      Y = 0.0_realk
      do k=1,nocc
         do d=1,nvirt
            do c=1,nvirt
               do b=1,nvirt
                  Lbdkc = 2.0_realk*gvvov%val(b,d,k,c) - gvvov%val(b,c,k,d)
                  do j=1,nocc
                     Y(j,b) = Y(j,b) + t2%elm4(c,k,d,j) * Lbdkc
                  end do
               end do
            end do
         end do
      end do

      RightOrLeft: if(righttrans) then  ! right-hand side

         ! Non-symmetrized deltaNS_{bjai} = R_ai [ - sum_{klc} t_kl^cb L_jlkc
         !                                         + sum_{cdk} t_{kj}^{cd} L_{bdkc} ]
         call mem_alloc(deltaNS,nvirt,nocc,nvirt,nocc)
         do a=1,nvirt
            do i=1,nocc
               Rai = V1%elm2(a,i)
               do b=1,nvirt
                  do j=1,nocc
                     deltaNS(b,j,a,i) = Rai * ( -Q(b,j) + Y(j,b) )
                  end do
               end do
            end do
         end do

         ! Update rho with symmetrized delta: 
         ! rho2_{aibj} --> rho2_{aibj} + P_{ij}^{ab} [ deltaNS(a,i,b,j) + deltaNS(b,j,a,i) ]
         do j=1,nocc
            do b=1,nvirt
               do i=1,nocc
                  do a=1,nvirt
                     rho2%elm4(a,i,b,j) = rho2%elm4(a,i,b,j) + deltaNS(a,i,b,j) + deltaNS(b,j,a,i) 
                  end do
               end do
            end do
         end do
         call mem_dealloc(deltaNS)

      else

         ! Left-hand side
         do i=1,nocc
            do a=1,nvirt
               do j=1,nocc
                  do b=1,nvirt
                     rho1%elm2(a,i) = rho1%elm2(a,i) &
                          & + ( Y(j,b) - Q(b,j) ) * V2%elm4(b,j,a,i)
                  end do
               end do
            end do
         end do

      end if RightOrLeft

      call mem_dealloc(Q)
      call mem_dealloc(Y)

      call array2_free(Co)
      call array2_free(Cv)
      call array4_free(gooov)
      call array4_free(gvvov)

    end subroutine add_lw1_contribution


  end module cc_response_tools_module
