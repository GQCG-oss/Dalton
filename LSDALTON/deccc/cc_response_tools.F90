module ccsd_lhtr_module

   use precision
   use typedef
   use typedeftype
   use dec_typedef_module

   ! DEC DEPENDENCIES (within deccc directory)   
   ! *****************************************
   use array4_simple_operations
   use cc_tools_module
   use ccintegrals

   public get_ccsd_multipliers_simple
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
   subroutine get_ccsd_multipliers_simple(rho1,rho2,t1f,t2f,m1,m2,fo,xo,yo,xv,yv,no,nv,nb,MyLsItem,gao_ex)
     implicit none

     type(lsitem), intent(inout) :: MyLsItem
     real(realk),intent(inout) :: rho1(:,:),rho2(:,:,:,:)
     real(realk),intent(inout) :: t1f(:,:),t2f(:,:,:,:),m1(:,:),m2(:,:,:,:),fo(:,:)
     real(realk),intent(in)    :: xo(:,:),yo(:,:),xv(:,:),yv(:,:)
     integer, intent(in)       :: no,nv,nb
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

     if(present(gao_ex))then
        gao = gao_ex
     else
        call get_full_eri(mylsitem,nb,gao)
     endif

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

     !ADD RIGHT HAND SIDES
     call array_reorder_4d(2.0E0_realk,Lovov,no,nv,no,nv,[2,1,4,3],1.0E0_realk,rho2)
     call mat_transpose(no,nv,2.0E0_realk,ovf,1.0E0_realk,rho1)

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
end module ccsd_lhtr_module
