      subroutine two2mean34b(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
      implicit double precision (a-h,o-z) 
      include 'para.h'
      logical*4 sameorb 
      dimension 
     *carteSO(norbsum,ncontmf,norbsum,ncontmf),
     *carteOO(norbsum,ncontmf,norbsum,ncontmf),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) then 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5D0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *carteSO(icartleft,jrun,icartright,irun)
      enddo
      enddo
      enddo
      enddo
      else 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb   
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5D0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *(carteSO(icartleft,jrun,icartright,irun)+
     *2d0*carteOO(icartleft,jrun,icartright,irun))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
      
