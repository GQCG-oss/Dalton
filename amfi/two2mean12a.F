      subroutine two2mean12a(carteSO,carteOO,occup,AOcoeffs,onecart,
     *ncontmf,norbsum,noccorb,sameorb)
      implicit double precision (a-h,o-z) 
      include 'para.h'
      logical*4 sameorb 
      dimension 
     *carteSO(ncontmf,norbsum,ncontmf,norbsum),
     *carteOO(ncontmf,norbsum,ncontmf,norbsum),
     *occup(*),AOcoeffs(MxcontL,*),onecart(MxcontL,MxcontL)
      if (sameorb) THEN 
      do icartleft=1,norbsum
      do icartright=1,norbsum
      coeff=0d0
      do Mrun=1,noccorb
      coeff=coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*
     *      AOcoeffs(icartright,Mrun)
      enddo 
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *carteSO(irun,icartleft,jrun,icartright)
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
      coeff=0.5d0*coeff  
      do irun=1,ncontmf
      do jrun=1,ncontmf
      onecart(irun,jrun)=onecart(irun,jrun)-coeff*
     *(carteSO(irun,icartleft,jrun,icartright)+
     *2d0*carteOO(irun,icartleft,jrun,icartright))
      enddo
      enddo
      enddo
      enddo
      endif 
      return 
      end 
      
