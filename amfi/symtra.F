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
#include <single.h>
C
      Subroutine symtrafo(LUPROP,LUPROP2,PROPRI,IPRINT,WRK,LFREE)
cbs   
cbs   Purpose: combine SO-integrals from amfi to symmetry-addapted   
cbs   integrals on one file AOPROPER_MF_SYM   
cbs   
cbs   LUPROP has the atomic integrals, LUPROP2 will have the molecular 
cbs   integrals    
cbs   
#include <implicit.h>
#include <dummy.h>
#include <priunit.h>
#include <para.h>
      parameter(maxorbs=700)
      parameter(maxcent=40)
      character*8 xa,ya,za,xa2
      character*3 END                     
CBS   character*20 filename               
      logical  EX, PROPRI
CBS   namelist /SYMTRA/ none
      dimension xa(4),ya(4),za(4)
      dimension xa2(4)
      dimension WRK(LFREE)
      dimension ncent(maxorbs), Lval(maxorbs),mval(maxorbs), 
     *nadpt(maxorbs),nphase(8,maxorbs),jdummy(8),    
     *Lhighcent(maxcent),Lcent(maxorbs,maxcent),Mcent(maxorbs,maxcent), 
     *ncontcent(0:Lmax,maxcent),ipC(3,Maxcent),
     *numballcart(maxcent),ifirstLM(0:Lmax,-Lmax:Lmax,maxcent)  ,
     *ncount(maxorbs)
      Character Label*8
c###############################################################################
      IPNT(I,J)=(max(i,j)*max(i,j)-max(i,j))/2 +min(i,j) 
      zero=0d0
      xa(1)='********'
      ya(1)='********'
      za(1)='********'
      xa(2)='        '
      ya(2)='        '
      Za(2)='        '
      xa(3)='ANTISYMM'
      ya(3)='ANTISYMM'
      Za(3)='ANTISYMM'
      xa(4)='X1MNF-SO'
      ya(4)='Y1MNF-SO'
      ZA(4)='Z1MNF-SO'  
c
c     read information from SYMINFO   
      ISYMUNIT = -1
      inquire(FILE='SYMINFO',EXIST=EX) 
      if (.not.EX) CALL QUIT('SYMINFO not present. Sorry')
      CALL GPOPEN(ISYMUNIT,'SYMINFO','UNKNOWN',' ','FORMATTED',IDUMMY,
     &            .FALSE.)
      rewind(isymunit)  
CBS   write(LUPRI,*) 'Symmetry adapation of the SO-integrals'
      read(isymunit,*)   
      read(isymunit,*)   
      read(isymunit,*)   
      numboffunct=0
      END = '   '
      do while(END.ne.'END') 
      numboffunct=numboffunct+1   
      read(isymunit,'(A3)') END   
      enddo    
CBS   write(LUPRI,*) 'there are totally ',numboffunct,' functions' 
      if (numboffunct.gt.maxorbs) THEN
         write(LUPRI,*) 'there are totally ',numboffunct,' functions'
         write(LUPRI,*) 'but maxorbs=',maxorbs
         CALL QUIT('increase maxorbs in symtrafo')
      endif
      rewind isymunit 
      read(isymunit,*)   
      read(isymunit,*)
      numbofcent=0
CBS   read(5,SYMTRA,END=4711)
4711  do irun=1,numboffunct
      read(isymunit,*) index,ncent(irun),lval(irun),
     *mval(irun),ncount(irun),nadpt(irun)
     *,(nphase(I,irun),I=1,nadpt(irun)) 
      numbofcent=max(numbofcent,ncent(irun)) 
Chj-start
      IF (LVAL(IRUN) .LT. 0) CALL QUIT('Negative L value on SYMINFO')
C     hjaaj nov 2001: I discovered that only s,p,d,f,g orbitals
C     was programmed in SYMOUT in herrdn.F, thus this consistency test.
C     Old code gave wrong lval for h orbitals etc.
C      -- modify SYMOUT if you need orbitals with higher L
C
Chj-end
      if (index.ne.irun) CALL QUIT('weird numbering  on SYMINFO')
      enddo 
      CALL GPCLOSE(ISYMUNIT,'KEEP')
CBS   write(LUPRI,*) 'number of unique centres' , numbofcent 
c
c     clean up arrays for new integrals  
      numboffunct3=(numboffunct*numboffunct+numboffunct)/2 
*
      nSCR=numboffunct3*3
      IPX = 1
      ipY=ipX+numboffunct3
      ipZ=ipY+numboffunct3
      IPSCR = IPZ + NUMBOFFUNCT3
      KLAST = IPZ + NSCR
      IF (KLAST .GT. LFREE) CALL STOPIT('SYMTRA',' ',KLAST,LFREE)
      Call DCopy(numboffunct3*3,Zero,0,WRK(ipx),1)
*
*
c     
c     loop over unique centres to read integrals and information   
c
      iunit=LUPROP
      nSCR=numboffunct3*3
      ipSCRX=ipSCR
      ipSCRY=ipSCRX+numboffunct3
      ipSCRZ=ipSCRY+numboffunct3    
      KLAST = IPZ + NSCR
      IF (KLAST .GT. LFREE) CALL STOPIT('SYMTRA',' ',KLAST,LFREE)
      length3_tot=0
*
      do icent=1,numbofcent    
*
CBS   write(LUPRI,*) 'read integrals and info for centre ',icent   
              read(iunit)  xa2
     *        ,numbofsym,(jdummy(I),
     *        i=1,numbofsym),
     *        numballcart(icent),(Lcent(i,icent),
     *        I=1,numballcart(icent)),
     *        (mcent(i,icent),I=1,numballcart(icent)),
     *        Lhighcent(icent),(ncontcent(I,icent),I=0,Lhighcent(icent))
CBS           write(LUPRI,*) numballcart(icent) , 
CBS  *        'functions on centre ',icent
              length3=ipnt(numballcart(icent),numballcart(icent)) 
              ipC(1,iCent)=ipSCRX
              ipC(2,iCent)=ipSCRY
              ipC(3,iCent)=ipSCRZ
              read(iunit) (Wrk(i),i=ipSCRX,ipSCRX+length3-1)
              read(iunit)  Ya
              read(iunit) (Wrk(i),i=ipSCRY,ipSCRY+length3-1)
              read(iunit)  Za
              read(iunit) (Wrk(i),i=ipSCRZ,ipSCRZ+length3-1)
              ipSCRX=ipSCRX+length3
              ipSCRY=ipSCRY+length3
              ipSCRZ=ipSCRZ+length3
              length3_tot=length3_tot+length3
culf
c      check if any L-value is missing
      LLhigh = Lhighcent(icent)
      do i=1,Lhighcent(icent)
         if(ncontcent(I,icent).eq.0) LLhigh=LLhigh-1
      enddo
      Lhighcent(icent)=LLhigh                                  
cbs   determize where the first function of a special type is..
      do Lrun=0,Lhighcent(icent) 
      do Mrun=-Lrun,Lrun 
      ifirstLM(Lrun,Mrun,icent)=ipnt(maxorbs,maxorbs)+1 
      enddo 
      enddo 
      do iorb=1,numballcart(icent)  
      Lrun=Lcent(iorb,icent) 
      Mrun=Mcent(iorb,icent) 
c     write(LUPRI,*) 'iorb,Lrun,mrun',iorb,Lrun,mrun 
      ifirstLM(Lrun,Mrun,icent)=min(iorb,ifirstLM(Lrun,Mrun,icent))   
      enddo 
cbs   determined..     
cbs   check if all of them were found 
      do Lrun=0,Lhighcent(icent) 
      do Mrun=-Lrun,Lrun 
      if(ifirstLM(Lrun,Mrun,icent).eq.(ipnt(maxorbs,maxorbs)+1)) then 
      write(LUPRI,*) 'problems for centre,L,M ',icent,Lrun,Mrun  
      CALL QUIT('problems with L- and M-values')
      endif 
      enddo 
      enddo 
      enddo    !end of loop over centres   
cbs
cbs   Finally the transformation!!!! 
cbs
cbs
      lauf=0
      do irun=1,numboffunct 
      icent=ncent(irun)
      ilval=lval(irun)
      imval=mval(irun)
      do jrun=1,irun        
      lauf=lauf+1
      jcent=ncent(jrun)
      jlval=lval(jrun)
      jmval=mval(jrun)
cbs   check for same centers 
      if (ncent(irun).eq.ncent(jrun)) then 
      if (lval(irun).eq.lval(jrun).and.lval(irun).gt.0) then 
      if (iabs(iabs(mval(irun))-iabs(mval(jrun))).le.1) then 
      if (iabs(iabs(mval(irun))+iabs(mval(jrun))).ne.0) then 
      if (irun.ne.jrun) then 
*
*
*
cbs   the only cases  where non-zero integrals occur    
      if (nadpt(irun).eq.1) then 
      coeff=1d0
      else 
      icoeff=0
      do icc=1,nadpt(irun)
      icoeff=icoeff+nphase(icc,irun)*nphase(icc,jrun)
      enddo 
      coeff=icoeff
CBS   coeff=coeff/(nadpt(irun)*nadpt(irun))  ! this is HERMIT  
      endif  
cbs   determine indices of atomic integrals 
      indexi=ifirstLM(ilval,imval,icent)+ncount(irun)-1
      indexj=ifirstLM(jlval,jmval,jcent)+ncount(jrun)-1
      laufalt=ipnt(indexi,indexj) 
*
      ipSCRX=ipC(1,ncent(irun))-1+laufalt
      ipSCRY=ipC(2,ncent(irun))-1+laufalt
      ipSCRZ=ipC(3,ncent(irun))-1+laufalt
      if (indexi.gt.indexj) then 
         Wrk(ipX+lauf-1)=coeff*Wrk(ipSCRX)
         Wrk(ipY+lauf-1)=coeff*Wrk(ipSCRY)
         Wrk(ipZ+lauf-1)=coeff*Wrk(ipSCRZ)
      else 
         Wrk(ipX+lauf-1)=-coeff*Wrk(ipSCRX)
         Wrk(ipY+lauf-1)=-coeff*Wrk(ipSCRY)
         Wrk(ipZ+lauf-1)=-coeff*Wrk(ipSCRZ)
      endif 
CBS
CBS  the integrals are now on on the WRK array 
CBS  WRK(ipx+i) (i=1,numboffunct3)
CBS  WRK(ipy+i) (i=1,numboffunct3)
CBS  WRK(ipz+i) (i=1,numboffunct3)
CBS
*
*   
      endif 
      endif 
      endif 
      endif 
      endif 
      enddo   
      enddo   
      if (lauf.ne.numboffunct3) CALL QUIT('error in numbering ')
*
COV   IPRINT=5
      IF (PROPRI .OR. IPRINT .GT. 4) THEN
         CALL AROUND('Integrals of operator: '//XA(4))
         CALL OUTPAK(WRK(IPX),NUMBOFFUNCT,1,LUPRI)
         CALL AROUND('Integrals of operator: '//YA(4))
         CALL OUTPAK(WRK(IPY),NUMBOFFUNCT,1,LUPRI)
         CALL AROUND('Integrals of operator: '//ZA(4))
         CALL OUTPAK(WRK(IPZ),NUMBOFFUNCT,1,LUPRI)
      END IF
      IPRINT=0
      CALL WRTPRO(WRK(IPX),NUMBOFFUNCT3,XA(4),XA(2))
      CALL WRTPRO(WRK(IPY),NUMBOFFUNCT3,YA(4),YA(2))
      CALL WRTPRO(WRK(IPZ),NUMBOFFUNCT3,ZA(4),ZA(2))
      Return
      End   
