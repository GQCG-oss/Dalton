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

end module cc_tools_module
