!Simple tools common for cc routines
module cc_tools_module

   use precision
   use tensor_interface_module
   
   
   contains
   
   subroutine mo_work_dist(m,fai,tl,have_part,nod)
      implicit none
      integer,intent(in) :: m
      integer,intent(inout)::fai
      integer,intent(inout)::tl
      logical,intent(out)  :: have_part
      integer(kind=ls_mpik),optional,intent(inout)::nod
      integer(kind=ls_mpik) :: nnod, me
      integer :: l,ml

      me   = 0
      nnod = 1
#ifdef VAR_MPI
      nnod = infpar%lg_nodtot
      me   = infpar%lg_mynum
#endif

      if(present(nod))me=nod

      !Setting transformation variables for each rank
      !**********************************************
      l   = (m) / nnod
      ml  = mod(m,nnod)
      fai = me * l + 1
      tl  = l

      if(ml>0)then
         if(me<ml)then
            fai = fai + me
            tl  = l + 1
         else
            fai = fai + ml
            tl  = l
         endif
      endif

      !If too many nodes are used set the first element to an invalid counter
      !**********************************************************************
      if(fai > m)then
         fai = -1
         have_part = .false.
      else
         have_part = .true.
      endif

   end subroutine mo_work_dist

   !> \brief Reorder t to use symmetry in both occupied and virtual indices,
   !> thereby restricting the first virtual to be less equal to the second and
   !> make symmetric and antisymmetric combinations of these 
   !> \author Patrick Ettenhuber
   !> \date October 2012
   subroutine get_tpl_and_tmi(t2,nv,no,tpl,tmi)
      implicit none
      !> number of virtual and occupied orbitals
      integer, intent(in) :: nv, no
      !>input the full doubles ampitudes in the order nv nv no no
      real(realk),intent(in) :: t2(nv,nv,no,no)
      !> output amplitudes with indices restricted c<d,i<j
      real(realk), intent(inout) :: tpl(:)
      !> output amplitudes with indices restricted c>d,i<j
      real(realk), intent(inout) :: tmi(:)
      integer :: dims(4)
      integer ::i,j,d, ilej,cged,crvd
      integer :: counter,counter2
      dims(1)=nv
      dims(2)=nv
      dims(3)=no
      dims(4)=no

      !get the combined reduced virtual dimension
      crvd = dims(1) * (dims(1)+1) / 2
      !OMP PARALLEL PRIVATE(i,j,d,cged,ilej) SHARED(crvd,t2,tmi,tpl,dims) DEFAULT(NONE)
      !OMP DO
      do j=1,dims(3)
         do i=1,j
            cged=1
            ilej = i + (j-1) * j / 2
            do d=1,dims(2)
               !copy the elements c>d into tpl and tmi
               call dcopy(dims(2)-d+1,t2(d,d,i,j),1,tpl(cged+(ilej-1)*crvd),1)
               call dcopy(dims(2)-d+1,t2(d,d,i,j),1,tmi(cged+(ilej-1)*crvd),1)

               !add and subtract the counterparts for tpl and tmi respectively
               call daxpy(dims(1)-d+1,-1.0E0_realk,t2(d,d,i,j),dims(1),tmi(cged+(ilej-1)*crvd),1)
               call daxpy(dims(1)-d+1,1.0E0_realk,t2(d,d,i,j),dims(1),tpl(cged+(ilej-1)*crvd),1)

               !take care of plus combination diagonal elements to be only half
               tpl(cged+(ilej-1)*crvd)=0.5E0_realk*tpl(cged+(ilej-1)*crvd)

               !count in the triangular part of the matrix
               cged = cged+dims(2)-d+1
            enddo
         enddo
      enddo
      !OMP END DO
      !OMP END PARALLEL
   end subroutine get_tpl_and_tmi

   !> \brief calculate a and b terms in a kobayashi fashion
   !> \author Patrick Ettenhuber
   !> \date December 2012
   subroutine get_a22_and_prepb22_terms_ex(w0,w1,w2,w3,tpl,tmi,no,nv,nb,fa,fg,la,lg,&
         &xo,yo,xv,yv,om2,sio4,s,wszes,lo,twork,tcomm,order,rest_occ_om2,scal)
      implicit none
      !> workspace with exchange integrals
      real(realk),intent(inout) :: w1(:)
      !> empty workspace of correct sizes
      real(realk),intent(inout) :: w0(:),w2(:),w3(:)
      !> the t+ and t- combinations with a value of the amplitudes with the
      !diagonal elements divided by two
      real(realk),intent(inout) :: tpl(:),tmi(:)
      !> number of occupied, virutal and ao indices
      integer, intent(in) :: no,nv,nb
      !> first alpha and first gamma indices of the current loop
      integer, intent(in) :: fa,fg
      !> lengths of the alpha ang gamma batches in the currnet loop
      integer, intent(in) :: la,lg
      !> the doubles residual to update
      !real(realk), intent(inout) :: om2(:)
      type(array), intent(inout) :: om2
      !> the lambda transformation matrices
      real(realk), intent(in)    :: xo(:),yo(:),xv(:),yv(:)
      !> the sio4 matrix to calculate the b2.2 contribution
      real(realk),intent(inout) ::sio4(:)
      !> scheme
      integer,intent(in) :: s
      logical,intent(in) :: lo
      !timing information
      real(realk) :: twork,tcomm
      !> W0 SIZE
      integer(kind=8),intent(in)  :: wszes(4)
      integer,optional,intent(in) :: order(4)
      logical,optional,intent(in) :: rest_occ_om2
      real(realk),optional :: scal
      integer(kind=8)  :: wszes3(3)
      integer :: goffs,aoffs,tlen,tred,nor,nvr

      call time_start_phase(PHASE_WORK)

      nor=(no*(no+1))/2
      nvr=(nv*(nv+1))/2

      !Determine the offsets in the alpha and gamma indices, which arise from
      !non uniform batch distributions, e.g. consider:
      !Case 1:   gamma _____
      !               |\    |
      !               |_\___|
      !  alpha        |  \  |
      !               |___\_| here an offset in gamma for the second gamma batch is
      !               needed, since the elements to be considered do not
      !               begin with the first element in that batch
      goffs=0
      if(fa-fg>0)goffs=fa-fg
      !Case 2:   gamma ______
      !               |\ |   |
      !               | \|   |
      !  alpha        |  |\  |
      !               |__|_\_| here an offset in alpha for the second alpha batch is
      !               needed, since the elements to be considered do not
      !               begin with the first element in that batch
      aoffs=0
      if(fg-fa>0)aoffs=fg-fa

      !Determine the dimension of the triangular part in the current batch
      !and the total number of elements (tred) in the triangular plus
      !rectangular parts
      tred=0
      tlen=min(min(min(la,lg),fg+lg-fa),fa+la-fg)

      !calculate amount of triangular elements in the current batch
      if(fa+la-1>=fg)tred= tlen*(tlen+1)/2
      !add the rectangular contribution if lg is larger than tlen 
      if(fa>=fg.and.fg+lg-fa-tlen>0)tred=tred+(fg+lg-fa-tlen)*la
      if(fa<fg)tred=tred+lg*aoffs
      if(fa<fg.and.fa+la>fg) tred=tred+(lg-tlen)*(la-aoffs)
      !if only rectangular
      if(tlen<=0)then
         tlen=0
         aoffs=0
         goffs=0
         tred=la*lg
      endif
      !!SYMMETRIC COMBINATION
      !(w0):I+ [delta alpha<=gamma beta] <= (w1):I [alpha beta gamma delta] + (w1):I[alpha delta gamma beta]
      call get_I_plusminus_le(w0,w1,w2,'p',fa,fg,la,lg,nb,tlen,tred,goffs)
      !(w2):I+ [delta alpha<=gamma c] = (w0):I+ [delta alpha<=gamma beta] * Lambda^h[beta c]
      call dgemm('n','n',nb*tred,nv,nb,1.0E0_realk,w0,nb*tred,yv,nb,0.0E0_realk,w2,nb*tred)
      !(w0):I+ [alpha<=gamma c d] = (w2):I+ [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
      call dgemm('t','n',tred*nv,nv,nb,1.0E0_realk,w2,nb,yv,nb,0.0E0_realk,w0,nv*tred)
      !(w2):I+ [alpha<=gamma c>=d] <= (w0):I+ [alpha<=gamma c d] 
      call get_I_cged(w2,w0,tred,nv)
      !(w3.1):sigma+ [alpha<=gamma i>=j] = (w2):I+ [alpha<=gamma c>=d] * (w0):t+ [c>=d i>=j]
      call dgemm('n','n',tred,nor,nvr,0.5E0_realk,w2,tred,tpl,nvr,0.0E0_realk,w3,tred)


      !!ANTI-SYMMETRIC COMBINATION
      !(w0):I- [delta alpha<=gamma beta] <= (w1):I [alpha beta gamma delta] + (w1):I[alpha delta gamma beta]
      call get_I_plusminus_le(w0,w1,w2,'m',fa,fg,la,lg,nb,tlen,tred,goffs)
      !(w2):I- [delta alpha<=gamma c] = (w0):I- [delta alpha<=gamma beta] * Lambda^h[beta c]
      call dgemm('n','n',nb*tred,nv,nb,1.0E0_realk,w0,nb*tred,yv,nb,0.0E0_realk,w2,nb*tred)
      !(w0):I- [alpha<=gamma c d] = (w2):I- [delta, alpha<=gamma c] ^T * Lambda^h[delta d]
      call dgemm('t','n',tred*nv,nv,nb,1.0E0_realk,w2,nb,yv,nb,0.0E0_realk,w0,nv*tred)
      !(w2):I- [alpha<=gamma c<=d] <= (w0):I- [alpha<=gamma c d] 
      call get_I_cged(w2,w0,tred,nv)
      !(w3.2):sigma- [alpha<=gamma i<=j] = (w2):I- [alpha<=gamma c>=d] * (w0):t- [c>=d i>=j]
      call dgemm('n','n',tred,nor,nvr,0.5E0_realk,w2,tred,tmi,nvr,0.0E0_realk,w3(tred*nor+1),tred)

      !COMBINE THE TWO SIGMAS OF W3 IN W2
      !(w2):sigma[alpha<=gamma i<=j]=0.5*(w3.1):sigma+ [alpha<=gamma i<=j] + 0.5*(w3.2):sigma- [alpha <=gamm i<=j]
      !(w2):sigma[alpha>=gamma i<=j]=0.5*(w3.1):sigma+ [alpha<=gamma i<=j] - 0.5*(w3.2):sigma- [alpha <=gamm i<=j]
      wszes3 = [wszes(1),wszes(3),wszes(4)]
      call combine_and_transform_sigma(om2,w0,w2,w3,xv,xo,sio4,nor,tlen,tred,fa,fg,la,lg,&
         &no,nv,nb,goffs,aoffs,s,wszes3,lo,twork,tcomm,order=order, &
         &rest_occ_om2=rest_occ_om2,scal=scal)  

      call time_start_phase(PHASE_WORK, at=twork)
   end subroutine get_a22_and_prepb22_terms_ex

   !> \brief Combine sigma matrixes in symmetric and antisymmetric combinations 
   !> \author Patrick Ettenhuber
   !> \date October 2012
   subroutine combine_and_transform_sigma(omega,w0,w2,w3,xvirt,xocc,sio4,nor, tlen,tred,fa,fg,&
         & la,lg,no,nv,nb,goffs,aoffs,s,wszes,lock_outside,twork,tcomm, order,rest_occ_om2,scal,act_no, sio4_ilej, query )
      implicit none
      !\> omega should be the residual matrix which contains the second parts
      !of the A2 and B2 term
      !real(realk),intent(inout) :: omega(nv*nv*no*no)
      type(array),intent(inout) :: omega
      !> w0 is just some workspace on input
      real(realk),intent(inout) :: w0(:)
      !> w2 is just some workspace on input
      real(realk),intent(inout) :: w2(:)
      !> w3 contains the symmetric and antisymmetric combinations 
      real(realk),intent(inout) :: w3(:)
      !> sio4 are the reduced o4 integrals whic are used to calculate the B2.2
      !contribution after the loop, update them in the loops
      real(realk),intent(inout) :: sio4(:)
      !> Lambda p virutal part
      real(realk),intent(in) :: xvirt(:)
      !> Lambda p occupied part
      real(realk),intent(in) :: xocc(:)
      !> number of reduced occupied indices 
      integer,intent(in) :: nor
      !> total number of upper triangular elements in the batch
      integer,intent(in) :: tlen
      !> first element of alpha and gamma in the current batch
      integer,intent(in) :: fa,fg
      !> length of alpha and gamma in the current batch
      integer,intent(in) :: la, lg
      !> number of triangular elements in the batch
      integer,intent(in) :: tred
      !>number of occupied, virtual and ao indices
      integer,intent(in) :: no, nv, nb
      !>offsets in the currnt alpha and gamma batches to get to the first upper
      !triangular element
      integer,intent(in) :: goffs, aoffs
      !> scheme
      integer,intent(in) :: s
      logical,intent(in) :: lock_outside
      !> size of w0
      integer(kind=8),intent(inout)   :: wszes(3)
      integer,optional,intent(in)  :: order(4),act_no
      !restricted i<=j in the omega2 and or sio4
      logical,optional, intent(in) :: rest_occ_om2, sio4_ilej, query
      real(realk), optional, intent(in) :: scal
      !timing information
      real(realk) :: twork,tcomm
      !> the doubles amplitudes
      !real(realk),intent(in) :: amps(nv*nv*no*no)
      !type(array),intent(in) :: amps
      real(realk) :: scaleitby
      integer(kind=8)       :: pos,pos2,pos21,i,j,k,dim_big,dim_small,ttri,tsq,nel2cp,ncph,pos1
      integer ::occ,gamm,alpha,case_sel,full1,full2,offset1,offset2
      integer :: l1,l2,lsa,lsg,gamm_i_b,a,b,full1T,full2T,jump,ft1,ft2
      logical               :: second_trafo_step
      real(realk),pointer   :: dumm(:)
      integer               :: mv((nv*nv)/2),st,dims(2),no2
      real(realk),pointer   :: source(:,:),drain(:,:)
      integer(kind=ls_mpik) :: mode
      integer(kind=long)    :: o2v2
      logical               :: rest_o2_occ, rest_sio4,qu
      real(realk), pointer  :: h1(:,:,:,:), t1(:,:,:)
      !$ integer, external  :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads


      rest_o2_occ   = .false.
      if(present(rest_occ_om2 ))rest_o2_occ   = rest_occ_om2
      no2 = no
      if(present(act_no))no2=act_no
      rest_sio4 = .true.
      if(present(sio4_ilej))rest_sio4 = sio4_ilej
      qu = .false.
      if(present(query))qu = query

      o2v2 = int((i8*no)*no*nv*nv,kind=long)
#ifdef VAR_MPI
      mode = MPI_MODE_NOCHECK
#endif

      scaleitby=1.0E0_realk
      if(present(scal)) scaleitby = scal
      second_trafo_step=.false.
      lsa=fa+la-1
      lsg=fg+lg-1

      !building the sigma matrices is a complicated matter, due to non-uniform
      !batchsizes. A distiction between the different possible cases has to be
      !made and the elements copied respectively
      case_sel=0
      if(fa>=fg)then
         if(lsa>=lsg)then
            !case 0 is a triagular submatrix:
            !                                full2=tlen
            ! ---|                          |----|
            ! \  |                          | \  |  full1=tlen
            !  \ |   build full from that   |  \ |
            !   \|                          |---\|
            case_sel=1
            full1=tlen
            full2=tlen
            full1T=tlen
            full2T=tlen
            l1=0
            l2=0
         else if(lsa<lsg)then
            !case 1 is a  submatrix of the type:
            !                                   full2=lg-goffs
            ! ------|                          |------|
            ! \     |                          | \    |
            !  \    |   build full from that   |  \   |full1=la
            !   \---|                          |___\__|
            !                                  |   |0 |
            !            full1T=lg-goffs-tlen  |___|__|
            case_sel=2                         !full2T=la 
            full1=la
            full2=lg-goffs
            full1T=lg-goffs-tlen
            full2T=la
            l1=la
            l2=0
            second_trafo_step=.true.
         endif
      else
         if(lsa>=lsg)then
            !case 2 is a  submatrix of the type:
            !                                 full2=lg
            ! |---|                         |------|
            ! |   |                         |0 |   |
            ! |   |  build full from that   |__|   | full1=aoffs+tlen
            !  \  |                         |  |\  |
            !   \ |                         |  | \ |
            !    \|     full1T=lg           |__|__\|
            !                               full2T=aoffs
            case_sel=3
            full1=aoffs+tlen
            full2=tlen
            full1T=tlen
            full2T=aoffs
            l1=aoffs
            l2=aoffs
            second_trafo_step=.true.
         else if(lsa<lsg)then
            if(lsa>=fg)then
               !case 3 is a  submatrix of the type: 
               !                                       full2=lg
               ! |-----|                         |--------|
               ! |     |                         |0 |     |
               ! |     |  build full from that   |__|     |
               !  \    |                         |  |\    |  full1 = la
               !   \   |                         |  | \   |
               !    \--|                         |  |__\__|
               !                     full1T=lg   |     | 0|
               !                                 |_____|__|
               !                                    full2T = la
               case_sel=4
               full1=la
               full2=lg
               full1T=lg
               full2T=la
               l1=aoffs
               l2=aoffs
               second_trafo_step=.true.
            else
               !case 4 is a  submatrix of dimensions la lg, and can thus be transformed in
               !full anyways
               case_sel=5
               full1=la
               full2=lg
               full1T=lg
               full2T=la
               l1=la
               l2=0
               second_trafo_step=.true.
            endif
         endif
      endif
      if(case_sel==0)call lsquit("ERROR(combine_and_transform_sigma):case not known",DECinfo%output)

      !print *,"-------------------------------------------------------------------------------------"
      !print *,"       case        dim1        dim2            la          lg        tlen        tred"
      !print *,case_sel,full1,full2,la,lg,tlen,tred
      !print *,"-------------------------------------------------------------------------------------"


      if( qu )then

         !ATTENTION: KEEP UP TO DATE

         !w0:
         if(second_trafo_step)then
            wszes(1) = max(wszes(1),(i8*nor)*full1*full2+(i8*nor)*full1T*full2T)
         else
            wszes(1) = max(wszes(1),(i8*nor)*full1*full2)
         endif

         !w2:
         wszes(2) = max(wszes(2),(i8*nor)*full1*full2)
         wszes(2) = max(wszes(2),(i8*nv*nv)*nor)
         if(.not.rest_o2_occ)then
            wszes(2) = max(wszes(2),(i8*nv*nv)*no*no)
         endif
         if(second_trafo_step)then
            wszes(2) = max(wszes(2),(i8*full1T)*full2T*nor)
            wszes(2) = max(wszes(2),(i8*nv)*nv*nor)
         endif
         if( .not. rest_sio4 )then
            wszes(2) = max(wszes(2),(i8*no2*no2)*nor)
         endif

         !w3:
         wszes(3) = max(wszes(3),(i8*nor)*full1*full2)
         wszes(3) = max(wszes(3),(i8*nv*nor*full1)*full2)
         wszes(3) = max(wszes(3),(i8*no2)*nor*full1)
         if(second_trafo_step)then
            wszes(3) = max(wszes(3),(i8*nv*nor)*full1T)
        endif

      else

         !Zero the elements to update for testing, not needed in a performance
         !implementation
         pos=(i8*nor)*full1*full2
         if(second_trafo_step)pos=pos+(i8*nor)*full1T*full2T
         w0(1_long:pos)=0.0E0_realk

         !set required variables
         ttri      = tlen * (tlen+1)/2
         tsq       = tlen *  tlen
         pos       = 1
         occ       = 1
         dim_big   = (i8 * full1 ) * full2
         dim_small = (i8 * full1T) * full2T

         !$OMP PARALLEL DEFAULT(NONE)&
         !$OMP SHARED(w0,w3,case_sel,nor,goffs,lg,la,full1,full1T,ttri,tred,&
         !$OMP full2,full2T,tlen,l1,second_trafo_step,aoffs,dim_big,dim_small,l2)&
         !$OMP PRIVATE(occ,gamm,gamm_i_b,pos,nel2cp,pos2,jump,ft1,ft2,ncph,pos21,&
         !$OMP dims,i)
         !$OMP DO
         do occ=1,nor
            do gamm=1,lg-goffs
               gamm_i_b=gamm+goffs
               !SYMMETRIC COMBINATION OF THE SIGMAS

               !calculate the old position
               !**************************
               if(case_sel==3.or.case_sel==4)then
                  pos=1+(gamm      -1)*full1+(occ-1)*dim_big
               else
                  pos=1+(gamm+aoffs-1)*full1+(occ-1)*dim_big
               endif

               !calculate the new position
               !**************************
               if(gamm>tlen)then

                  !get the elements from the rectangular part of the batch

                  nel2cp = l1
                  pos2   = 1_long + ttri+(gamm-tlen-1)*(la-aoffs)+(occ-1)*((la-aoffs)*(lg-goffs-tlen)+ttri)

                  if(case_sel==4)then
                     nel2cp = nel2cp + tlen
                     pos2   = pos2   + tlen * aoffs + (gamm-tlen-1) * (la-tlen)+(occ-1)*(aoffs*lg)
                  endif

                  if(second_trafo_step)then
                     jump = full1T
                     ft1  = full1T
                     ft2  = full2T
                  else
                     jump = full1
                     ft1  = full1
                     ft2  = full2
                  endif

               else
                  !get the elements from the triangular part of the batch
                  pos2   = 1  + (gamm*(gamm-1)/2)+(gamm-1)*aoffs+(occ-1)*tred
                  nel2cp = l2 + gamm
                  jump   = full1
                  ft1    = full1
                  ft2    = full2
               endif

               !call dcopy(nel2cp,w3(pos2),1,w0(pos),1)
               w0(pos:pos+nel2cp-1) = w3(pos2:pos2+nel2cp-1)

               !get corresponding position in sigma- and add to output
               pos21=pos2+tred*nor
               !call daxpy(nel2cp,1.0E0_realk,w3(pos21),1,w0(pos),1)
               w0(pos:pos+nel2cp-1) =w0(pos:pos+nel2cp-1) + w3(pos21:pos21+nel2cp-1)    

               !ANTI-SYMMETRIC COMBINATION OF THE SIGMAS
               pos = gamm+aoffs+(occ-1)*ft1*ft2

               if(second_trafo_step.and.gamm>tlen) pos = pos+full1*full2*nor-tlen

               if(case_sel==3.or.case_sel==4)then
                  !fill diagonal part
                  if(gamm>tlen)then
                     ncph=0
                  else
                     ncph=gamm
                  endif
                  !call daxpy(ncph,1.0E0_realk,w3(pos2+aoffs),1,w0(aoffs+gamm+(occ-1)*full1*full2),full1)
                  !call daxpy(ncph,-1.0E0_realk,w3(pos21+aoffs),1,w0(aoffs+gamm+(occ-1)*dim_big),full1)

                  !because of the intrinsic omp-parallelizaton of daxpy the following
                  !lines replace the daxpy calls
                  pos = aoffs+gamm+(occ-1)*dim_big
                  do i=0,ncph-1
                     w0(pos+i*full1)=w0(pos+i*full1) + w3(pos2 +aoffs+i)
                     w0(pos+i*full1)=w0(pos+i*full1) - w3(pos21+aoffs+i)
                  enddo

                  !fill small matrix
                  if(gamm>tlen)then
                     ncph=nel2cp
                  else
                     ncph=nel2cp-gamm
                  endif
                  !call daxpy(ncph, 1.0E0_realk,w3(pos2 ),1,w0(gamm+(occ-1)*full1T*full2T+dim_big*nor),full1T)
                  !call daxpy(ncph,-1.0E0_realk,w3(pos21),1,w0(gamm+(occ-1)*full1T*full2T+dim_big*nor),full1T)
                  pos = gamm+(occ-1)*full1T*full2T+dim_big*nor
                  do i=0,ncph-1
                     w0(pos+i*full1T) = w0(pos+i*full1T) + w3(pos2 +i)
                     w0(pos+i*full1T) = w0(pos+i*full1T) - w3(pos21+i)
                  enddo
               else
                  !call daxpy(nel2cp,1.0E0_realk,w3(pos2),1,w0(pos),jump)
                  !call daxpy(nel2cp,-1.0E0_realk,w3(pos21),1,w0(pos),jump)
                  do i=0,nel2cp-1
                     w0(pos+i*jump) = w0(pos+i*jump) + w3(pos2 +i)
                     w0(pos+i*jump) = w0(pos+i*jump) - w3(pos21+i)
                  enddo
               endif
            enddo
         enddo
         !$OMP END DO
         !$OMP END PARALLEL


         !Print the individual contributions
         !if(.false.)then
         !  do occ=1,nor
         !    print *,"symmetric"
         !    print *,w3(1+(occ-1)*tred:tred+(occ-1)*tred)
         !    !print *,"antisymmetric"
         !    !print *,w3(nor*tred+1+(occ-1)*tred:nor*tred+tred+(occ-1)*tred)
         !    print *,"a g -> 1. and l.",1+(occ-1)*full1*full2,full1*full2+(occ-1)*full1*full2
         !    do alpha=1,full1
         !      print *,(w0(alpha+(gamm-1)*full1+(occ-1)*full1*full2),gamm=1,full2)
         !    enddo
         !    print *,""
         !    if(second_trafo_step)then
         !      print *,"g a -> 1. and l.",dim_big*nor+1+(occ-1)*full1T*full2T,dim_big*nor+full1T*full2T+(occ-1)*full1T*full2T
         !      do gamm=1,full1T
         !        print *,(w0(dim_big*nor+gamm+(alpha-1)*full1T+(occ-1)*full1T*full2T),alpha=1,full2T)
         !      enddo
         !      print *,""
         !    endif
         !    print *,""
         !  enddo
         !endif

         !add up the contributions for the sigma [ alpha gamma ] contributions
         ! get the order sigma[ gamma i j alpha ]
         call mat_transpose(full1,full2*nor,1.0E0_realk,w0,0.0E0_realk,w2)
         !transform gamma -> b
         call dgemm('t','n',nv,nor*full1,full2,1.0E0_realk,xvirt(fg+goffs),nb,w2,full2,0.0E0_realk,w3,nv)
         !transform alpha -> a , order is now sigma [ a b i j]
         call dgemm('t','t',nv,nv*nor,full1,1.0E0_realk,xvirt(fa),nb,w3,nor*nv,0.0E0_realk,w2,nv)


         if(.not.rest_o2_occ)then
            !square up the contributions if the residual itself has no restricions
            !in the indices i and j

            ! add up contributions in the residual with keeping track of i<j
            !$OMP PARALLEL DEFAULT(NONE) SHARED(no,w2,nv)&
            !$OMP PRIVATE(i,j,k,pos1,pos2) IF( no > 2 )
            do j=no,2,-1
               do i=j,1,-1
                  pos1=((i+j*(j-1)/2)-1)*nv*nv
                  pos2=(i-1)*nv*nv+(j-1)*no*nv*nv
                  !$OMP DO 
                  do k=1,nv*nv
                     w2(k+pos2) = w2(k+pos1)
                  enddo
                  !$OMP END DO
               enddo
            enddo
            !$OMP BARRIER
            !$OMP DO 
            do j=no,1,-1
               do i=j,1,-1
                  pos1=1+(i-1)*nv*nv+(j-1)*no*nv*nv
                  pos2=1+(j-1)*nv*nv+(i-1)*no*nv*nv
                  if(i/=j) w2(pos2:pos2+nv*nv-1) = w2(pos1:pos1+nv*nv-1)
               enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL

            do j=no,1,-1
               do i=j,1,-1
                  pos1=1+(i-1)*nv*nv+(j-1)*no*nv*nv
                  call alg513(w2(pos1:nv*nv+pos1-1),nv,nv,nv*nv,mv,(nv*nv)/2,st)
               enddo
            enddo

            if(s==4.or.s==3)then
               if( present(order) )then
                  call array_reorder_4d(scaleitby,w2,nv,nv,no,no,order,1.0E0_realk,omega%elm1)
               else
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
                  call assign_in_subblocks(omega%elm1,'+',w2,o2v2,scal2=scaleitby)
#else
                  !OMP WORKSHARE
                  omega%elm1(1_long:o2v2) = omega%elm1(1_long:o2v2) + scaleitby * w2(1_long:o2v2)
                  !OMP END WORKSHARE
#endif
               endif
            else if(s==2)then
#ifdef VAR_MPI
               if( lock_outside )call arr_lock_wins(omega,'s',mode)
               !$OMP WORKSHARE
               w2(1_long:o2v2) = scaleitby*w2(1_long:o2v2)
               !$OMP END WORKSHARE
               call time_start_phase(PHASE_COMM, at=twork)
               call array_add(omega,1.0E0_realk,w2,wrk=w3,iwrk=wszes(3))
               call time_start_phase(PHASE_WORK, at=tcomm)
#endif
            endif
         else
#ifdef VAR_LSDEBUG
            if(s/=4)call lsquit("ERROR(combine_and_transform_sigma)not implemented&
               & for other schemes than 4",-1)
            if(s/=4)call lsquit("ERROR(combine_and_transform_sigma)only implemented&
               & for pairs for now",-1)
#endif
            pos1 = 1
            call array_reorder_2d(scaleitby,w2(pos1:nv*nv+pos1-1),nv,nv,[2,1],1.0E0_realk,omega%elm1)

         endif

         !If the contributions are split in terms of the sigma matrix this adds the
         !sigma [gamma alpha] contributions
         if(second_trafo_step)then
            !get the order sigma [gamma i j alpha]
            pos2 = 1+full1*full2*nor
            if(case_sel==3.or.case_sel==4)then
               l1=fa
               l2=fg
            else
               l1=fa
               l2=fg+goffs+tlen
            endif
            call mat_transpose(full1T,full2T*nor,1.0E0_realk,&
               &w0(pos2:full1T*full2T*nor+pos2-1),0.0E0_realk,w2)

#ifdef  VAR_MPI
            if( lock_outside .and. s==2 )then
               call time_start_phase(PHASE_COMM, at=twork)
               call arr_unlock_wins(omega,.true.)
               call time_start_phase(PHASE_WORK, at=tcomm)
            endif
#endif 

            !transform gamma -> a
            call dgemm('t','n',nv,nor*full1T,full2T,1.0E0_realk,xvirt(l1),nb,w2,full2T,0.0E0_realk,w3,nv)
            !transform alpha -> , order is now sigma[b a i j]
            call dgemm('t','t',nv,nv*nor,full1T,1.0E0_realk,xvirt(l2),nb,w3,nor*nv,0.0E0_realk,w2,nv)

            if(.not.rest_o2_occ)then
               !$OMP PARALLEL DEFAULT(NONE) SHARED(no,w2,nv)&
               !$OMP PRIVATE(i,j,k,pos1,pos2) IF( no > 2 )
               do j=no,2,-1
                  do i=j,1,-1
                     pos1=1+((i+j*(j-1)/2)-1)*nv*nv
                     pos2=1+(i-1)*nv*nv+(j-1)*no*nv*nv
                     !$OMP DO 
                     do k=1,nv*nv
                        w2(k+pos2) = w2(k+pos1)
                     enddo
                     !$OMP END DO
                  enddo
               enddo
               !$OMP BARRIER
               !$OMP DO 
               do j=no,1,-1
                  do i=j,1,-1
                     pos1=1+(i-1)*nv*nv+(j-1)*no*nv*nv
                     pos2=1+(j-1)*nv*nv+(i-1)*no*nv*nv
                     if(i/=j) w2(pos2:pos2+nv*nv-1) = w2(pos1:pos1+nv*nv-1)
                  enddo
               enddo
               !$OMP END DO
               !$OMP END PARALLEL

               do j=no,1,-1
                  do i=j,1,-1
                     pos1=1+(i-1)*nv*nv+(j-1)*no*nv*nv
                     call alg513(w2(pos1:nv*nv+pos1-1),nv,nv,nv*nv,mv,(nv*nv)/2,st)
                  enddo
               enddo
               if(s==4.or.s==3)then
                  if( present(order) )then
                     call array_reorder_4d(scaleitby,w2,nv,nv,no,no,order,1.0E0_realk,omega%elm1)
                  else
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
                     call assign_in_subblocks(omega%elm1,'+',w2,o2v2,scal2=scaleitby)
#else
                     !OMP WORKSHARE
                     omega%elm1(1_long:o2v2) = omega%elm1(1_long:o2v2) + scaleitby * w2(1_long:o2v2)
                     !OMP END WORKSHARE
#endif
                  endif
               else if(s==2)then
                  !$OMP WORKSHARE
                  w2(1_long:o2v2) = scaleitby*w2(1_long:o2v2)
                  !$OMP END WORKSHARE
                  call time_start_phase(PHASE_COMM, at=twork)
                  call array_add(omega,1.0E0_realk,w2,wrk=w3,iwrk=wszes(3))
                  call time_start_phase(PHASE_WORK, at=tcomm)
               endif
            else
#ifdef VAR_LSDEBUG
               if(s/=4)call lsquit("ERROR(combine_and_transform_sigma)not implemented&
                  & for other schemes than 4",-1)
               if(s/=4)call lsquit("ERROR(combine_and_transform_sigma)only implemented&
                  & for pairs for now",-1)
#endif
               pos1 = 1
               call array_reorder_2d(scaleitby,w2(pos1:nv*nv+pos1-1),nv,nv,[2,1],1.0E0_realk,omega%elm1)
            endif
         endif


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !Construct the B2 term from the intermediates in w0
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !add up the contributions for the sigma [ alpha gamma ] contributions
         ! get the order sigma[ gamma i j alpha ]
         call mat_transpose(full1,full2*nor,1.0E0_realk,w0,0.0E0_realk,w2)
#ifdef VAR_MPI
         if(lock_outside.and.s==2)then
            call arr_unlock_wins(omega,.true.)
         endif
#endif
         !transform gamma -> l
         call dgemm('t','n',no2,nor*full1,full2,1.0E0_realk,xocc(fg+goffs),nb,w2,full2,0.0E0_realk,w3,no2)
         !transform alpha -> a , order is now sigma [ k l i j]
         if( rest_sio4 )then
            call dgemm('t','t',no2,no2*nor,full1,1.0E0_realk,xocc(fa),nb,w3,nor*no2,1.0E0_realk,sio4,no2)
         else

            call dgemm('t','t',no2,no2*nor,full1,1.0E0_realk,xocc(fa),nb,w3,nor*no2,0.0E0_realk,w2,no2)
            call ass_D1to3(w2,t1,[no2,no2,nor])
            call ass_D1to4(sio4,h1,[no,no,no2,no2])
            do j=no,1,-1
               do i=j,1,-1
                  call array_reorder_2d(1.0E0_realk,t1(:,:,i+j*(j-1)/2),no2,no2,[2,1],1.0E0_realk,h1(i,j,:,:))
                  if(i /= j)then
                     h1(j,i,:,:) =h1(j,i,:,:) +  t1(:,:,i+j*(j-1)/2)
                  endif
               enddo
            enddo
         endif




         !If the contributions are split in terms of the sigma matrix this adds the
         !sigma [gamma alpha i j] contributions
         if(second_trafo_step)then
            !get the order sigma [gamma i j alpha]

            pos2 = 1+full1*full2*nor
            if(case_sel==3.or.case_sel==4)then
               l1=fa
               l2=fg
            else
               l1=fa
               l2=fg+goffs+tlen
            endif
            call mat_transpose(full1T,full2T*nor,1.0E0_realk,&
               &w0(pos2:full1T*full2T*nor+pos2-1),0.0E0_realk,w2)
            !transform gamma -> l
            call dgemm('t','n',no2,nor*full1T,full2T,1.0E0_realk,xocc(l1),nb,w2,full2T,0.0E0_realk,w3,no2)
            !transform alpha -> k, order is now sigma[k l i j]
            if( rest_sio4 )then
               call dgemm('t','t',no2,no2*nor,full1T,1.0E0_realk,xocc(l2),nb,w3,nor*no2,1.0E0_realk,sio4,no2)
            else
               call dgemm('t','t',no2,no2*nor,full1T,1.0E0_realk,xocc(l2),nb,w3,nor*no2,0.0E0_realk,w2,no2)
               call ass_D1to3(w2,t1,[no2,no2,nor])
               call ass_D1to4(sio4,h1,[no,no,no2,no2])
               do j=no,1,-1
                  do i=j,1,-1
                     call array_reorder_2d(1.0E0_realk,t1(:,:,i+j*(j-1)/2),no2,no2,[2,1],1.0E0_realk,h1(i,j,:,:))
                     if(i /= j)then
                        h1(j,i,:,:) =h1(j,i,:,:) + t1(:,:,i+j*(j-1)/2)
                     endif
                  enddo
               enddo
            endif

         endif
      endif


   end subroutine combine_and_transform_sigma


   !> \brief Construct integral matrix which uses symmetries in the virtual
   !> indices
   !> \author Patrick Ettenhuber
   !> \date October 2012
   subroutine get_I_cged(Int_out,Int_in,m,nv)
      implicit none
      !> integral output with indices reduced to c>=d
      real(realk),intent(inout)::Int_out(:)
      !>full integral input m*c*d
      real(realk),intent(in) :: Int_in(:)
      !> leading dimension m and virtual dimension
      integer,intent(in)::m,nv
      integer ::d,pos,pos2,a,b,c,cged
      logical :: doit
#ifdef VAR_OMP
      integer :: tid,nthr
      integer, external :: omp_get_thread_num,omp_get_max_threads
      !    nthr = omp_get_max_threads()
      !    nthr = min(nthr,nv)
      !    call omp_set_num_threads(nthr)
#endif
      !OMP PARALLEL DEFAULT(NONE) SHARED(int_in,int_out,m,nv,nthr)&
      !OMP PRIVATE(pos,pos2,d,tid,doit)
#ifdef VAR_OMP
      !    tid = omp_get_thread_num()
#else 
      !    doit = .true.
#endif
      !    doit = .true.
      pos =1
      do d=1,nv
#ifdef VAR_OMP
         !      doit = (mod(d,nthr) == tid)
#endif
         !      if(doit) then
         pos2=1+(d-1)*m+(d-1)*nv*m
         !        call dcopy(m*(nv-d+1),Int_in(pos2),1,Int_out(pos),1)
         Int_out(pos:pos+m*(nv-d+1)-1) = Int_in(pos2:pos2+m*(nv-d+1)-1)
         !      endif
         pos=pos+m*(nv-d+1)
      enddo
      !OMP END PARALLEL
#ifdef VAR_OMP
      !call omp_set_num_threads(omp_get_max_threads())
#endif
   end subroutine get_I_cged


   !> \brief Construct symmetric and antisymmentric combinations of an itegral
   !matrix 
   !> \author Patrick Ettenhuber
   !> \date October 2012
   subroutine get_I_plusminus_le(w0,w1,w2,op,fa,fg,la,lg,nb,tlen,tred,goffs,qu,quarry)
      implicit none
      !> blank workspace
      real(realk),intent(inout) :: w0(:),w2(:)
      !> workspace containing the integrals
      real(realk),intent(in) :: w1(:)
      !> integer specifying the first element in alpha and gamma batch
      integer,intent(in) :: fa,fg
      !> integer specifying the length of alpha and gamma batches
      integer,intent(in) :: la,lg
      !> number of upper triangular elements in the current batch
      integer,intent(in) :: tred
      !> number of elements in triangular part
      integer,intent(in) :: tlen
      !> number of ao indices
      integer,intent(in) :: nb
      !> offset in the gamma batch
      integer,intent(in) :: goffs
      !> character specifying which operation to use
      character, intent(in) :: op
      !> query the array sizes
      logical, intent(in), optional :: qu
      integer(kind=8), intent(inout), optional :: quarry(3)
      integer :: i,alpha,beta,gamm,delta,cagc,cagi,bs,bctr
      integer :: alpha_b,beta_b,gamm_b,delta_b,elements,aleg
      integer :: globg, globa, loca,eldiag,elsqre,nbnb,nrnb,omp_level
      real(realk) ::chk,chk2,el
      real(realk),pointer :: trick(:,:,:)
      logical :: modb,query

      query = .false.
      if(present(qu).and.present(quarry)) query = qu

!#ifdef VAR_OMP
!      integer, external :: omp_get_level
!#endif
      omp_level = 0
!#ifdef VAR_OMP
!      omp_level = omp_get_level()
!#endif
      bs=int(sqrt(((8.0E6_realk)/1.6E1_realk)))
      !bs=5
      !print *,"block size",bs,(bs*bs*8)/1024.0E0_realk
      nbnb=(nb/bs)*bs
      modb=(mod(nb,bs)>0)
      bctr = bs-1
      cagi=tred

      if(query)then
         !IMPORTANT:
         !Make sure these sizes are up to date
         quarry(1) = max(quarry(1),(i8*la*nb)*lg*nb)
         quarry(2) = max(quarry(2),(i8*la*nb)*lg*nb)
         quarry(3) = max(quarry(3),(i8*nb*nb)*cagi)
      else
         call ass_D1to3(w2,trick,[nb,nb,cagi])
         if(op=='p')then
            call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
            aleg=0
            do gamm=0,lg-1
               do alpha=0,la-1
                  if(fa+alpha<=fg+gamm)then
                     !aleg = (alpha+(gamm*(gamm+1))/2) 
                     eldiag = aleg*nb*nb
                     elsqre = alpha*nb*nb+gamm*nb*nb*la
                     !print *,alpha,gamm,1+eldiag,nb*nb+eldiag,1+elsqre,nb*nb+elsqre,aleg,cagi,nb*nb
                     if(fa+alpha==fg+gamm)   call dscal(nb*nb,0.5E0_realk,w2(1+elsqre),1)
                     call dcopy(nb*nb,w2(1+elsqre),1,w2(1+eldiag),1)
                     !$OMP PARALLEL PRIVATE(el,delta_b,beta_b,beta,delta)&
                     !$OMP SHARED(bs,bctr,trick,nb,aleg,nbnb,modb)&
                     !$OMP DEFAULT(NONE)
                     if(nbnb>0)then
                        !$OMP DO
                        do delta_b=1,nbnb,bs
                           do beta_b=delta_b+bs,nbnb,bs
                              do delta=0,bctr
                                 do beta=0,bctr
                                    el=trick(beta+beta_b,delta_b+delta,aleg+1)
                                    trick(beta+beta_b,delta_b+delta,aleg+1)=&
                                       &trick(beta+beta_b,delta_b+delta,aleg+1) + &
                                       &trick(delta_b+delta,beta+beta_b,aleg+1)
                                    trick(delta_b+delta,beta+beta_b,aleg+1)=&
                                       &trick(delta_b+delta,beta+beta_b,aleg+1) + &
                                       &el
                                 enddo
                              enddo
                           enddo
                        enddo
                        !$OMP END DO NOWAIT
                     endif
                     if(nbnb>0.and.modb)then
                        !$OMP DO
                        do delta_b=1,nbnb,bs
                           do delta=0,bctr
                              do beta=nbnb+1,nb
                                 el=trick(beta,delta_b+delta,aleg+1)
                                 trick(beta,delta_b+delta,aleg+1)=&
                                    &trick(beta,delta_b+delta,aleg+1) + &
                                    &trick(delta_b+delta,beta,aleg+1)
                                 trick(delta_b+delta,beta,aleg+1)=&
                                    &trick(delta_b+delta,beta,aleg+1) + &
                                    &el
                              enddo
                           enddo
                        enddo
                        !$OMP END DO NOWAIT
                     endif
                     if(nbnb>0)then
                        !$OMP DO
                        do delta_b=1,nbnb,bs
                           do delta=0,bctr
                              do beta=delta+1,bctr
                                 el=trick(beta+delta_b,delta_b+delta,aleg+1)
                                 trick(beta+delta_b,delta_b+delta,aleg+1)=&
                                    &trick(beta+delta_b,delta_b+delta,aleg+1) + &
                                    &trick(delta_b+delta,beta+delta_b,aleg+1)
                                 trick(delta_b+delta,beta+delta_b,aleg+1)=&
                                    &trick(delta_b+delta,beta+delta_b,aleg+1) + &
                                    &el
                              enddo
                              trick(delta+delta_b,delta_b+delta,aleg+1)=&
                                 &trick(delta+delta_b,delta_b+delta,aleg+1) + &
                                 &trick(delta_b+delta,delta+delta_b,aleg+1)
                           enddo
                        enddo
                        !$OMP END DO NOWAIT
                     endif
                     !$OMP END PARALLEL 
                     if(modb)then
                        do delta=nbnb+1,nb
                           do beta=delta+1,nb
                              el=trick(beta,delta,aleg+1)
                              trick(beta,delta,aleg+1)=&
                                 &trick(beta,delta,aleg+1) + &
                                 &trick(delta,beta,aleg+1)
                              trick(delta,beta,aleg+1)=&
                                 &trick(delta,beta,aleg+1) + &
                                 &el
                           enddo
                           trick(delta,delta,aleg+1)=&
                              &trick(delta,delta,aleg+1) + &
                              &trick(delta,delta,aleg+1)
                        enddo
                     endif
                     aleg=aleg+1
                  endif
               enddo
            enddo
            call array_reorder_3d(1.0E0_realk,w2,nb,nb,cagi,[2,3,1],0.0E0_realk,w0)
         endif


         !Calculate the antisymmetric combination in the same way
         if(op=='m')then
            call array_reorder_4d(1.0E0_realk,w1,la,nb,lg,nb,[2,4,1,3],0.0E0_realk,w2)
            aleg=0
            do gamm=0,lg-1
               do alpha=0,la-1
                  if(fa+alpha<=fg+gamm)then
                     !aleg = (alpha+(gamm*(gamm+1))/2) 
                     eldiag = aleg*nb*nb
                     elsqre = alpha*nb*nb+gamm*nb*nb*la
                     call dcopy(nb*nb,w2(1+elsqre),1,w2(1+eldiag),1)
                     !$OMP PARALLEL PRIVATE(el,delta_b,beta_b,beta,delta)&
                     !$OMP SHARED(bctr,bs,trick,nb,aleg,nbnb,modb)&
                     !$OMP DEFAULT(NONE)
                     if(nbnb>0)then
                        !$OMP DO
                        do delta_b=1,nbnb,bs
                           do beta_b=delta_b+bs,nbnb,bs
                              do delta=0,bctr
                                 do beta=0,bctr
                                    el=trick(beta+beta_b,delta_b+delta,aleg+1)
                                    trick(beta+beta_b,delta_b+delta,aleg+1)=&
                                       &trick(delta_b+delta,beta+beta_b,aleg+1)- &
                                       &trick(beta+beta_b,delta_b+delta,aleg+1) 
                                    trick(delta_b+delta,beta+beta_b,aleg+1)=&
                                       &el - trick(delta_b+delta,beta+beta_b,aleg+1) 
                                 enddo
                              enddo
                           enddo
                        enddo
                        !$OMP END DO NOWAIT
                     endif
                     if(nbnb>0.and.modb)then
                        !$OMP DO
                        do delta_b=1,nbnb,bs
                           do delta=0,bctr
                              do beta=nbnb+1,nb
                                 el=trick(beta,delta_b+delta,aleg+1)
                                 trick(beta,delta_b+delta,aleg+1)=&
                                    &trick(delta_b+delta,beta,aleg+1)- &
                                    &trick(beta,delta_b+delta,aleg+1) 
                                 trick(delta_b+delta,beta,aleg+1)=&
                                    &el -trick(delta_b+delta,beta,aleg+1)
                              enddo
                           enddo
                        enddo
                        !$OMP END DO NOWAIT
                     endif
                     if(nbnb>0)then
                        !$OMP DO
                        do delta_b=1,nbnb,bs
                           do delta=0,bctr
                              do beta=delta+1,bctr
                                 el=trick(beta+delta_b,delta_b+delta,aleg+1)
                                 trick(beta+delta_b,delta_b+delta,aleg+1)=&
                                    &trick(delta_b+delta,beta+delta_b,aleg+1) - &
                                    &trick(beta+delta_b,delta_b+delta,aleg+1)
                                 trick(delta_b+delta,beta+delta_b,aleg+1)=&
                                    &el - trick(delta_b+delta,beta+delta_b,aleg+1)
                              enddo
                              trick(delta_b+delta,delta+delta_b,aleg+1)=&
                                 &trick(delta_b+delta,delta+delta_b,aleg+1) - &
                                 &trick(delta+delta_b,delta_b+delta,aleg+1)
                           enddo
                        enddo
                        !$OMP END DO NOWAIT
                     endif
                     !$OMP END PARALLEL 
                     if(modb)then
                        do delta=nbnb+1,nb
                           do beta=delta+1,nb
                              el=trick(beta,delta,aleg+1)
                              trick(beta,delta,aleg+1)=&
                                 &trick(delta,beta,aleg+1) - &
                                 &trick(beta,delta,aleg+1)
                              trick(delta,beta,aleg+1)=&
                                 &el- trick(delta,beta,aleg+1)
                           enddo
                           trick(delta,delta,aleg+1)=&
                              &trick(delta,delta,aleg+1) - &
                              &trick(delta,delta,aleg+1)
                        enddo
                     endif
                     aleg=aleg+1
                  endif
               enddo
            enddo
            call array_reorder_3d(1.0E0_realk,w2,nb,nb,cagi,[2,3,1],0.0E0_realk,w0)
         endif
         nullify(trick)
      endif
   end subroutine get_I_plusminus_le

   !> \brief subroutine to add contributions to the sio4 matrix which enters the
   !B2.2 term in the "non"-parallel region
   !> \author Patrick Ettenhuber
   !> \Date September 2012
   subroutine add_int_to_sio4(w0,w2,w3,nor,no,nv,nb,fa,fg,la,lg,xo,sio4,act_no)
      implicit none
      !> workspace containing the paritially transformed integrals ordered as I(i
      !gamma alpha j)
      real(realk),pointer :: w0(:)
      !> arbitrary workspace of correct size
      real(realk),pointer :: w2(:),w3(:)
      !> number of occupied, virutal and ao indices
      integer, intent(in) :: nor,no,nv,nb
      !> first alpha and first gamma indices of the current loop
      integer, intent(in) :: fa,fg
      !> lengths of the alpha ang gamma batches in the currnet loop
      integer, intent(in) :: la,lg
      !> transformation matrix for t1 transformed integrals "Lambda p"
      real(realk),intent(in) :: xo(nb*no)
      integer, intent(in), optional :: act_no
      !> sio4 storage space to update during the batched loops
      real(realk),pointer :: sio4(:)
      integer :: pos,i,j,no2
      integer(kind=8) :: pos1, pos2

      no2 = no
      if(present(act_no))no2 = act_no

      ! (w3):I[alpha gamma i j] <- (w0):I[i gamma alpha j]
      call array_reorder_4d(1.0E0_realk,w0,no,lg,la,no,[2,3,1,4],0.0E0_realk,w2)
      ! (w2):I[alpha gamma i <= j] <- (w3):I[alpha gamma i j]
      do j=1,no
         do i=1,j
            pos1=1_long+((i+j*(j-1)/2)-1)*la*(lg*i8)
            pos2=1_long+(i-1)*la*lg+(j-1)*la*lg*(no*i8)
            !call dcopy(la*lg,w2(1+(i-1)*la*lg+(j-1)*la*lg*no),1,w3(pos),1)
            w3(pos1:pos1+la*lg-1) = w2(pos2:pos2+la*lg-1)
         enddo
      enddo
      ! (w3):I[ gamma i <= j alpha] <- (w2):I[alpha gamma i <= j]
      call array_reorder_3d(1.0E0_realk,w3,lg,la,nor,[2,3,1],0.0E0_realk,w2)
      ! (w2):I[ l i <= j alpha] <- (w3):Lambda^p [gamma l ]^T I[gamma i <= j alpha]
      call dgemm('t','n',no2,nor*lg,la,1.0E0_realk,xo(fa),nb,w2,la,0.0E0_realk,w3,no2)
      ! (sio4):I[ k l i <= j] <-+ (w2):Lambda^p [alpha k ]^T I[l i <= j alpha]^T
      call dgemm('t','t',no2,nor*no2,lg,1.0E0_realk,xo(fg),nb,w3,nor*no2,1.0E0_realk,sio4,no2)

   end subroutine add_int_to_sio4

   !> \brief Get the b2.2 contribution constructed in the kobayashi scheme after
   !the loop to avoid steep scaling ste  !> \author Patrick Ettenhuber
   !> \date December 2012
   subroutine get_B22_contrib_mo(sio4,t2,w1,w2,no,nv,om2,s,lock_outside,tw,tc,no_par,order)
      implicit none
      !> the sio4 matrix from the kobayashi terms on input
      real(realk), intent(in) :: sio4(:)
      !> amplitudes
      !real(realk), intent(in) :: t2(*)
      type(array), intent(inout) :: t2
      !> some workspave
      real(realk), intent(inout) :: w1(:)
      real(realk), pointer :: w2(:)
      !> number of occupied, virutal and ao indices
      integer, intent(in) :: no,nv
      !> residual to be updated
      !real(realk), intent(inout) :: om2(*)
      type(array), intent(inout) :: om2
      !> integer specifying the calc-scheme
      integer, intent(in) :: s
      logical, intent(in) :: lock_outside
      !> work and communication time in B2 term
      real(realk), intent(inout) :: tw,tc
      logical, intent(in),optional :: no_par
      integer, intent(in),optional :: order(4)
      integer :: nor
      integer :: ml,l,tl,fai,lai
      integer :: tri,fri
      integer(kind=ls_mpik) :: nod,me,nnod,massa,mode
      real(realk) :: nrm1,nrm2,nrm3,nrm4
      integer ::  mv((nv*nv)/2),st
      integer(kind=8) :: o2v2,pos1,pos2,i,j,k,pos
      logical :: traf,np
      integer :: o(4)

      call time_start_phase(PHASE_WORK)

      np = .false.
      if(present(no_par))np = no_par
      o = [1,2,3,4]
      if(present(order)) o  = order

      me    = 0
      massa = 0
      nnod  = 1
#ifdef VAR_MPI
      massa = infpar%master
      nnod  = infpar%lg_nodtot
      me    = infpar%lg_mynum
      mode  = int(MPI_MODE_NOCHECK,kind=ls_mpik)
#endif
      o2v2=(i8*no)*no*nv*nv

      !Setting transformation variables for each rank
      !**********************************************
      if(np)then
         fai = 1
         tl  = nv*nv
      else
         call mo_work_dist(nv*nv,fai,tl,traf)
      endif

      if(DECinfo%PL>3.and.me==0)then
         write(DECinfo%output,'("Trafolength in striped B2:",I5)')tl
      endif

      nor=no*(no+1)/2

      ! do contraction
      if((s==4.or.s==3).and.traf)then


         w1=0.0E0_realk
         call dgemm('n','n',tl,nor,no*no,0.5E0_realk,t2%elm1(fai),nv*nv,sio4,no*no,0.0E0_realk,w1(fai),nv*nv)


      else if(s==2.and.traf)then

#ifdef VAR_MPI
         call mem_alloc(w2,tl*no*no)

         call time_start_phase(PHASE_COMM, at = tw )

         if(lock_outside)call arr_lock_wins(t2,'s',mode)
         call array_two_dim_1batch(t2,[1,2,3,4],'g',w2,2,fai,tl,lock_outside)
         if(lock_outside)call arr_unlock_wins(t2,.true.)

         call time_start_phase(PHASE_WORK, at = tc )

         w1=0.0E0_realk
         call dgemm('n','n',tl,nor,no*no,0.5E0_realk,w2,tl,sio4,no*no,0.0E0_realk,w1(fai),nv*nv)
         call mem_dealloc(w2)
#endif

      endif


      !$OMP PARALLEL DEFAULT(NONE) SHARED(no,w1,nv)&
      !$OMP PRIVATE(i,j,k,pos1,pos2)
      do j=no,2,-1
         do i=j,1,-1
            pos1=((i+j*(j-1)/2)-1)*nv*nv
            pos2=(i-1)*nv*nv+(j-1)*no*nv*nv
            !$OMP DO 
            do k=1,nv*nv
               w1(k+pos2) = w1(k+pos1)
            enddo
            !$OMP END DO
         enddo
      enddo
      !$OMP BARRIER
      !$OMP DO 
      do j=no,1,-1
         do i=j,1,-1
            pos1=1+(i-1)*nv*nv+(j-1)*no*nv*nv
            pos2=1+(j-1)*nv*nv+(i-1)*no*nv*nv
            if(i/=j) w1(pos2:pos2+nv*nv-1) = w1(pos1:pos1+nv*nv-1)
         enddo
      enddo
      !$OMP END DO
      !$OMP END PARALLEL

      do j=no,1,-1
         do i=j,1,-1
            pos1=1+(i-1)*nv*nv+(j-1)*no*nv*nv
            call alg513(w1(pos1:nv*nv+pos1-1),nv,nv,nv*nv,mv,(nv*nv)/2,st)
         enddo
      enddo

      if((s==4.or.s==3).and.traf)then

         if(present(order))then
            call array_reorder_4d(1.0E0_realk,w1,nv,nv,no,no,order,1.0E0_realk,om2%elm1)
         else
#ifdef VAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN
         call assign_in_subblocks(om2%elm1,'+',w1,o2v2)
#else
         !$OMP WORKSHARE
         om2%elm1(1:o2v2) = om2%elm1(1:o2v2) + w1(1:o2v2)
         !$OMP END WORKSHARE
#endif
         endif

      else if(s==2.and.traf)then

#ifdef VAR_MPI
         call time_start_phase(PHASE_COMM, at = tw )

         if(lock_outside)call arr_lock_wins(om2,'s',mode)
         call array_two_dim_1batch(om2,o,'a',w1,2,1,nv*nv,lock_outside) 
         call time_start_phase(PHASE_WORK, at = tc )
#endif

      endif

      call time_start_phase(PHASE_COMM, at = tw )
   end subroutine get_B22_contrib_mo



   subroutine print_tensor_unfolding_with_labels(mat,d1,l1,m1,d2,l2,m2,label)
      implicit none
      integer, intent(in) :: m1,m2
      integer, intent(in) :: d1(m1),d2(m2)
      character,intent(in) :: l1(m1),l2(m2)
      real(realk), intent(in) :: mat(*)
      character*(*), intent(in) :: label
      integer :: i,a,b,id1(m1),id2(m2),cd1,cd2

      print *,label(1:len(label))
      write (*,'("(")',advance='no') 
      cd1 = 1
      do i=1,m1
         if(i>1)write (*,'(X)',advance='no')
         write (*,'(A)',advance='no') l1(i)
         cd1 = cd1 * d1(i)
      enddo
      write (*,'(")  /  (")',advance='no') 
      cd2 = 1
      do i=1,m2
         if(i>1)write (*,'(X)',advance='no')
         write (*,'(A)',advance='no') l2(i)
         cd2 = cd2 * d2(i)
      enddo
      write (*,'(")  ")',advance='no') 

      do b = 1, cd2
         call get_midx(b,id2,d2,m2)
         write (*,'(12X,"(")',advance='no')
         do i = 1, m2
            if(i>1)write (*,'("/")',advance='no')
            write (*,'(I4)',advance='no') id2(i)
         enddo
         write (*,'(")")',advance='no')
      enddo
      write(*,*)

      do a = 1, cd1
         call get_midx(a,id1,d1,m1)
         write (*,'(" (")',advance='no') 
         do i = 1, m1
            if(i>1)write (*,'("/")',advance='no')
            write (*,'(I4)',advance='no') id1(i)
         enddo
         write (*,'(") ")',advance='no') 
         do b = 1, cd2
            write (*,'(g18.5)',advance='no') mat(a+(b-1)*cd1)
            if(b<cd2)write (*,'(",",X)',advance='no')
         enddo
         if(a<cd1)then
            write (*,'(";")')
         else
            write (*,*)
         endif
      enddo



   end subroutine print_tensor_unfolding_with_labels

end module cc_tools_module
