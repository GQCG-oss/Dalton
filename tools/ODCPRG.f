      PROGRAM ODCPRG
C  
C   Optical Data Calculations
C   program for hyperpolarizabilities and
C   shielding polarizabilities via Finite Field 
C   on analytical magnetizabilities and nuclear
C   shieldings
C   Version March 1, 1997  by A. Rizzo, S.Coriani 
C   
      implicit double precision (a-h,o-z)
      logical*1 zprt
      logical*1 case1,case2,case3,case4,case5,case6,case7
      character*2 satom
      character*3 cpg
      character*80 filenam,Title
      character*9 casl
      dimension satom(10),chi(9),shield(9,10),field(3)
      dimension en(27),fie(27,3),cmag(27,9),cshi(27,9,10)
      dimension casl(27),ifound(27),cpg(6)
      data casl/' -x -y -z', ' -x -y  0', ' -x -y +z',
     1          ' -x  0 -z', ' -x  0  0', ' -x  0 +z',
     2          ' -x +y -z', ' -x +y  0', ' -x +y +z',
     3          '  0 -y -z', '  0 -y  0', '  0 -y +z',
     4          '  0  0 -z', '  0  0  0', '  0  0 +z',
     5          '  0 +y -z', '  0 +y  0', '  0 +y +z',
     6          ' +x -y -z', ' +x -y  0', ' +x -y +z',
     7          ' +x  0 -z', ' +x  0  0', ' +x  0 +z',
     8          ' +x +y -z', ' +x +y  0', ' +x +y +z'/
      data cpg/' Td','Civ','D2h','C2v','C3v','Dih'/
      data ifound/27*0/
      iunit=22
c
c Script file file Unit 5
c
      open (unit=5,file='readmg.dat')
c
c Input
c
 9999 read (5,'(a)',end=9998) Title
      read (5,*,end=9997) nfiles,zprt,ipg,issym,nist,nlast,nstep
      write (6,'(//A//)') Title
      write (6,'(3x,A,I3,A)') '...Reading ',
     #nfiles,' data files from DALTON...'
c
      do 123 ij=1,27
         ifound(ij)=0.0
  123 continue  
      do 101 j=1,nfiles
        read (5,'(a)',end=9996) filenam
        write (6,'(///3x,3a)') '...Now opening file ',filenam,'...'
        call rdfile (iunit,filenam,ish,satom,chi,shield,field,energy,
     #  zprt,imag)
c
c Find out field setup
c      
        i=1
        if (field(1)) 2,3,4
    3   i=i+1
        go to 2
    4   i=i+2
    2   ii=1
        if (field(2)) 20,30,40
   30   ii=ii+1
        go to 20
   40   ii=ii+2
   20   iii=1
        if (field(3)) 200,300,400
  300   iii=iii+1
        go to 200
  400   iii=iii+2
  200   index=(i-1)*9+(ii-1)*3+iii
        ifound(index)=1
c
        en(index)=energy
        do jj=1,3
          fie(index,jj)=field(jj)
        end do
        if (imag.eq.0) then
          do jj=1,9
            cmag(index,jj)=chi(jj)
          end do
        end if
        if (ish.gt.0) then
         do jj=1,9
          do jjj=1,ish
           cshi(index,jj,jjj)=shield(jj,jjj)
          end do
         end do
        end if
  101 continue
c
      write (6,'(///a)') 
     #    '... I found data for the following field setups:'
      do jota=1,27
       if (ifound(jota).eq.1) write (6,'(10x,a)') casl(jota)
      end do
c
      write (6,'(/A/)') Title
      write (6,'(///10x,a)')     ' Summary of our Results:'
c
      ipgr=ipg
   91 continue
      write(6,'(10x,a,i)')'Warning:  ipgr=', ipgr
      go to (109,209,309,409,509,609) ipgr
c   --------------------------------------------------
c    Case Td    (No symmetry - 3 calculations):  CH4
c 
c     Need     0  0  0
c     Need     0  0  z
c     Need     x  0  z
c
  109 case1=ifound(14).eq.1                         
      case2=ifound(15).eq.1 .or. ifound(13).eq.1     
      case3=ifound(24).eq.1 .or. ifound(22).eq.1.or.
     #      ifound(6).eq.1 .or. ifound(4).eq.1
      if (case1.and.case2.and.case3) then
        ifirst=14        
        isecnd=13 
        if (ifound(15).eq.1) isecnd=15
        ithird=4
        if (ifound(24).eq.1) ithird=24
        if (ifound(22).eq.1) ithird=22  
        if (ifound(6).eq.1) ithird=6   
c      
        if (imag.eq.0) then        
c
          cxyz=(cmag(isecnd,2)-cmag(ifirst,2))/
     #                fie(isecnd,3)
          ezzzz=2.d0*(cmag(isecnd,9)-cmag(ifirst,9))/
     #               (fie(isecnd,3)*fie(isecnd,3))
          exxzz=2.d0*(cmag(isecnd,1)-cmag(ifirst,1))/
     #               (fie(isecnd,3)*fie(isecnd,3))
          exzxz=(cmag(ithird,3)-cmag(ifirst,3))/
     #               (fie(ithird,1)*fie(ithird,3))
c      
          deta=2.d0*(ezzzz-exxzz+3.d0*exzxz)/5.d0
c
          write (6,'(10x,a,/)') ' (Hyper)magnetizabilities'
          write (6,'(10x,a,/)') 'eta(ab,gd): ab=BaBb gd=FgFd'
          write (6,'(10x,a,g20.5,/)') '...Point group...', ipg
          write (6,'(10x,a,g20.5,/)') '...Anysotropy...    ',deta
c
          write (6,'(10x,a,g20.5,/)') '...Chi(zz)...',cmag(ifirst,1)
          write (6,'(10x,a,g20.5,/)') '...CSI(xyz)      ',cxyz
          write (6,'(10x,a,g20.5)') '...ETA(zzzz)     ',ezzzz
          write (6,'(10x,a,g20.5)') '...ETA(xxzz)     ',exxzz
          write (6,'(10x,a,g20.5,/)') '...ETA(xzxz)     ',exzxz
c
        end if                                       
        if (ish.gt.0) then            
         if (issym.eq.ipgr) then
          do 119 iatom=nist,nlast,nstep
c
            Sxyz=(cshi(isecnd,2,iatom)-cshi(ifirst,2,iatom))/
     #                  fie(isecnd,3)
            Szzzz=2.d0*(cshi(isecnd,9,iatom)-cshi(ifirst,9,iatom))/
     #                 (fie(isecnd,3)*fie(isecnd,3))
            Sxxzz=2.d0*(cshi(isecnd,1,iatom)-cshi(ifirst,1,iatom))/
     #                 (fie(isecnd,3)*fie(isecnd,3))
            Sxzxz=(cshi(ithird,3,iatom)-cshi(ifirst,3,iatom))/
     #                 (fie(ithird,1)*fie(ithird,3))
            B=-(2.d0*Sxxzz+Szzzz)/6.d0
c
            write (6,'(//10x,a,a,/)')
     1 ' Shielding (Polarizabilities) for atom ',satom(iatom)
            write (6,'(10x,a,i,/)') 'Site Symmetry  ', issym
            write (6,'(10x,a,g20.5,/)')'SIGMA(zz)  ',cshi(ifirst,9,1)
            write (6,'(10x,a,g20.5,/)') '...SIGMA(xyz)      ',sxyz
            write (6,'(10x,a,g20.5)') '...SIGMA(zzzz)     ',szzzz
            write (6,'(10x,a,g20.5)') '...SIGMA(xxzz)     ',sxxzz
            write (6,'(10x,a,g20.5,/)') '...SIGMA(xzxz)     ',sxzxz
            write (6,'(10x,a,g20.5)') '...B     ',B
  119     continue
         else
          ipgr=issym
          imag=1
          goto 91
         end if
        end if         
      else
        write(6,*)'Insufficient # of fields for required symmetry !!!'
        write(6,*)'Check value set in readmg.dat'
      end if         
      go to 9999
c --------------------------------------------------
c  Case Civ (CO as ex. No symmetry - 5 calculations)
c
c     Need     0  0  0
c     Need     x  0  0
c     Need     0  0  z
c     Need     0  0 -z
c     Need     x  0  z
c
  209 case1=ifound(14).eq.1                         
      case2=ifound(13).eq.1.and.ifound(15).eq.1    
      case3=ifound(23).eq.1.or.ifound(5).eq.1     
      case4=ifound(24).eq.1.or.ifound(22).eq.1.or.
     #      ifound(6).eq.1.or.ifound(4).eq.1
      if (case1.and.case2.and.case3.and.case4) then 
        ifirst=14
        isecnd1=15
        isecnd2=13
        ithird=5
        if (ifound(23).eq.1) ithird=23
        ifourth=4
        if (ifound(24).eq.1) ifourth=24
        if (ifound(22).eq.1) ifourth=22
        if (ifound(6).eq.1) ifourth=6
c
        if (imag.eq.0) then   
c
          czzz=
     #((cmag(isecnd1,9)-cmag(ifirst,9))*fie(isecnd2,3)*fie(isecnd2,3)-
     # (cmag(isecnd2,9)-cmag(ifirst,9))*fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
          cxxz=
     #((cmag(isecnd1,1)-cmag(ifirst,1))*fie(isecnd2,3)*fie(isecnd2,3)-
     # (cmag(isecnd2,1)-cmag(ifirst,1))*fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
          cxzx=(cmag(ithird,3)-cmag(ifirst,3))/fie(ithird,1)
c
          exxxx=2.d0*(cmag(ithird,1)-cmag(ifirst,1))/
     #             (fie(ithird,1)*fie(ithird,1))
c
          ezzzz=2.d0*((cmag(isecnd2,9)-cmag(ifirst,9))*fie(isecnd1,3)-
     # (cmag(isecnd1,9)-cmag(ifirst,9))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
          ezzxx=2.d0*(cmag(ithird,9)-cmag(ifirst,9))/
     #           (fie(ithird,1)*fie(ithird,1))
c
          exxzz=2.d0*((cmag(isecnd2,1)-cmag(ifirst,1))*fie(isecnd1,3)-
     # (cmag(isecnd1,1)-cmag(ifirst,1))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
          eyyxx=2.d0*(cmag(ithird,5)-cmag(ifirst,5))/
     #           (fie(ithird,1)*fie(ithird,1))
c
          exzxz=(cmag(ifourth,3)-cmag(ifirst,3)-cxzx*
     #            fie(ifourth,1))/
     #           (fie(ifourth,1)*fie(ifourth,3))
c
          deta=(exxxx*7.d0-eyyxx*5.d0+ezzzz*2.d0-exxzz*2.d0
     #         -ezzxx*2.d0+exzxz*12.d0)/15.d0
c
          cmagav=(cmag(ifirst,1)+cmag(ifirst,5)+cmag(ifirst,9))/3.d0
c
          dchi=cmag(ifirst,9)-cmag(ifirst,1)
c
          write (6,'(10x,a,/)') ' (Hyper)magnetizabilities'
          write (6,'(10x,a,/)') 'eta(ab,gd): ab=BaBb gd=FgFd'
          write (6,'(10x,a,i,/)') '...Point group...', ipgr
          write (6,'(10x,a,g20.5,/)') '...Anysotropy...    ',deta
          write (6,'(10x,a,g20.5)') '...Chi(xx)     ',cmag(ifirst,1)
          write (6,'(10x,a,g20.5)') '...Chi(xx)     ',cmag(ifirst,9)
          write (6,'(10x,a,g20.5)') '...Chi(av)     ',cmagav
          write (6,'(10x,a,g20.5,/)') '...DeltaChi     ',dchi
          write (6,'(10x,a,g20.5)') '...CSI(xxz)      ',cxxz
          write (6,'(10x,a,g20.5)') '...CSI(zzz)      ',czzz
          write (6,'(10x,a,g20.5,/)') '...CSI(xzx)      ',cxzx
          write (6,'(10x,a,g20.5)') '...ETA(xxxx)     ',exxxx
          write (6,'(10x,a,g20.5)') '...ETA(zzzz)     ',ezzzz
          write (6,'(10x,a,g20.5)') '...ETA(xxzz)     ',exxzz
          write (6,'(10x,a,g20.5)') '...ETA(yyxx)     ',eyyxx
          write (6,'(10x,a,g20.5)') '...ETA(zzxx)     ',ezzxx
          write (6,'(10x,a,g20.5,/)') '...ETA(xzxz)     ',exzxz
        end if
c
        if (ish.gt.0) then
         if (issym.eq.ipgr) then
          do 219 iatom=nist,nlast,nstep
c
             Sxxxx=2.d0*(cshi(ithird,1,iatom)-cshi(ifirst,1,iatom))/
     #           (fie(ithird,1)*fie(ithird,1))
c
             Szzzz=2.d0*
     # ((cshi(isecnd2,9,iatom)-cshi(ifirst,9,iatom))*fie(isecnd1,3)-
     # (cshi(isecnd1,9,iatom)-cshi(ifirst,9,iatom))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
             Szzxx=2.d0*(cshi(ithird,9,iatom)-cshi(ifirst,9,iatom))/
     #           (fie(ithird,1)*fie(ithird,1))
c
             Sxxzz=2.d0*
     # ((cshi(isecnd2,1,iatom)-cshi(ifirst,1,iatom))*fie(isecnd1,3)-
     # (cshi(isecnd1,1,iatom)-cshi(ifirst,1,iatom))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
             Syyxx=2.d0*(cshi(ithird,5,iatom)-cshi(ifirst,5,iatom))/
     #           (fie(ithird,1)*fie(ithird,1))
c
             Szxx=(cshi(ithird,7,iatom)-cshi(ifirst,7,iatom))/
     #              fie(ithird,1)
c
             Sxzx=(cshi(ithird,3,iatom)-cshi(ifirst,3,iatom))/
     #              fie(ithird,1)
c
             Sxzxz=(cshi(ifourth,3,iatom)-cshi(ifirst,3,iatom)-
     #               sxzx*fie(ifourth,1))/
     #               (fie(ifourth,1)*fie(ifourth,3))
c
             Szxxz=(cshi(ifourth,7,iatom)-cshi(ifirst,7,iatom)-
     #               szxx*fie(ifourth,1))/
     #               (fie(ifourth,1)*fie(ifourth,3))
c
             Szzz=
     # ((cshi(isecnd1,9,iatom)-cshi(ifirst,9,iatom))*
     # fie(isecnd2,3)*fie(isecnd2,3)-
     # (cshi(isecnd2,9,iatom)-cshi(ifirst,9,iatom))*
     # fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
             Sxxz=
     # ((cshi(isecnd1,1,iatom)-cshi(ifirst,1,iatom))*
     # fie(isecnd2,3)*fie(isecnd2,3)-
     # (cshi(isecnd2,1,iatom)-cshi(ifirst,1,iatom))*
     # fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
             Az=-(2.d0*Sxxz+Szzz)/3.d0
             Bzz=-(2.d0*Sxxzz+Szzzz)/6.d0
             Bxx=-(Sxxxx+Syyxx+Szzxx)/6.d0
c
             sigave=(cshi(ifirst,1,iatom)+cshi(ifirst,5,iatom)+
     #                cshi(ifirst,9,iatom))/3.d0
c     
             write (6,'(//10x,a,a,/)')
     1 ' Shielding (Polarizabilities) for atom ',satom(iatom)
             write (6,'(10x,a,i)') 'Site Symmetry   ', issym
             write (6,'(10x,a,g20.5)')'SIGMA(xx)  ',cshi(ifirst,1,iatom)
             write (6,'(10x,a,g20.5)')'SIGMA(zz)  ',cshi(ifirst,9,iatom)
             write (6,'(10x,a,g20.5,/)') '...SIGMA(av.)     ',sigave
             write (6,'(10x,a,g20.5)') '...SIGMA(xxz)      ',sxxz
             write (6,'(10x,a,g20.5)') '...SIGMA(xzx)      ',sxzx
             write (6,'(10x,a,g20.5)') '...SIGMA(zxx)      ',szxx
             write (6,'(10x,a,g20.5,/)') '...SIGMA(zzz)      ',szzz
             write (6,'(10x,a,g20.5)') '...SIGMA(xxxx)     ',sxxxx
             write (6,'(10x,a,g20.5)') '...SIGMA(zzzz)     ',szzzz
             write (6,'(10x,a,g20.5)') '...SIGMA(xxzz)     ',sxxzz
             write (6,'(10x,a,g20.5)') '...SIGMA(zzxx)     ',szzxx
             write (6,'(10x,a,g20.5)') '...SIGMA(yyxx)     ',syyxx
             write (6,'(10x,a,g20.5)') '...SIGMA(xzxz)     ',sxzxz
             write (6,'(10x,a,g20.5,/)') '...SIGMA(zxxz)     ',szxxz
             write (6,'(10x,a,g20.5)') '...Az     ',Az
             write (6,'(10x,a,g20.5)') '...Bxx     ',Bxx
             write (6,'(10x,a,g20.5)') '...Bzz     ',Bzz
  219     continue
         else
          ipgr=issym
          imag=1
          goto 91
         endif
        end if
      else
        write(6,*)'Insufficient # of fields for required symmetry !!!'
        write(6,*)'Check value set in readmg.dat'
      end if
      go to 9999
c  -----------------------------------------------
c  Case D2h (C2H4  No symmetry - 8 calculations)
c
c     Need     0  0  0
c     Need     0  0  z
c     Need     x  0  0
c     Need     0  y  0
c     Need     x  y  0
c     Need     x  0  z
c     Need     0  y  z
c
  309 case1=ifound(14).eq.1                      
      case2=ifound(13).eq.1.or.ifound(15).eq.1 
      case3=ifound(23).eq.1.or.ifound(5).eq.1  
      case4=ifound(11).eq.1.or.ifound(17).eq.1
      case5=ifound(2).eq.1.or.ifound(8).eq.1.or.
     #       ifound(20).eq.1.or.ifound(26).eq.1
      case6=ifound(4).eq.1.or.ifound(6).eq.1.or.
     #       ifound(22).eq.1.or.ifound(24).eq.1
      case7=ifound(18).eq.1.or.ifound(10).eq.1
     #      .or.ifound(12).eq.1.or.ifound(16).eq.1
      if (case1.and.case2.and.case3.and.case4
     #         .and.case5.and.case6.and.case7) then     
        ifirst=14
        isecnd=15
        if (ifound(13).eq.1) isecnd=13  
        ithird=5 
        if (ifound(23).eq.1) ithird=23  
        ifourth=11                      
        if (ifound(17).eq.1) ifourth=17
        ififth=26                      
        if (ifound(2).eq.1) ififth=2   
        if (ifound(8).eq.1) ififth=8   
        if (ifound(20).eq.1) ififth=20 
        isixth=24                     
        if (ifound(4).eq.1) isixth=4  
        if (ifound(6).eq.1) isixth=6  
        if (ifound(22).eq.1) isixth=22
        isevth=18                    
        if (ifound(10).eq.1) isevth=10
        if (ifound(12).eq.1) isevth=12 
        if (ifound(16).eq.1) isevth=16
c
        if (imag.eq.0) then    
c
          Exxxx=2.d0*(cmag(ithird,1)-cmag(ifirst,1))/
     #           (fie(ithird,1)*fie(ithird,1))
c
          Eyyyy=2.d0*(cmag(ifourth,5)-cmag(ifirst,5))/
     #           (fie(ifourth,2)*fie(ifourth,2))
c      
          Ezzzz=2.d0*(cmag(isecnd,9)-cmag(ifirst,9))/
     #           (fie(isecnd,3)*fie(isecnd,3))
c
          Eyyxx=2.d0*(cmag(ithird,5)-cmag(ifirst,5))/
     #           (fie(ithird,1)*fie(ithird,1))
c
          Ezzxx=2.d0*(cmag(ithird,9)-cmag(ifirst,9))/
     #           (fie(ithird,1)*fie(ithird,1))
c
          Exxyy=2.d0*(cmag(ifourth,1)-cmag(ifirst,1))/
     #           (fie(ifourth,2)*fie(ifourth,2))
c
          Ezzyy=2.d0*(cmag(ifourth,9)-cmag(ifirst,9))/
     #           (fie(ifourth,2)*fie(ifourth,2))
c
          Exxzz=2.d0*(cmag(isecnd,1)-cmag(ifirst,1))/
     #           (fie(isecnd,3)*fie(isecnd,3))
c
          Eyyzz=2.d0*(cmag(isecnd,5)-cmag(ifirst,5))/
     #           (fie(isecnd,3)*fie(isecnd,3))
c      
          Exyxy=(cmag(ififth,2)-cmag(ifirst,2))/
     #           (fie(ififth,1)*fie(ififth,2))
c      
          Exzxz=(cmag(isixth,3)-cmag(ifirst,3))/
     #           (fie(isixth,1)*fie(isixth,3))
c      
          Eyzyz=(cmag(isevth,6)-cmag(ifirst,6))/
     #           (fie(isevth,2)*fie(isevth,3))
c      
          deta=(2.d0*exxxx+2.d0*eyyyy+2.d0*ezzzz
     #         +6.d0*exyxy+6.d0*exzxz+6.d0*eyzyz
     #         -exxyy-eyyxx-exxzz-ezzxx-eyyzz-ezzyy)/15.d0
c
          cmagav=(cmag(ifirst,1)+cmag(ifirst,5)+cmag(ifirst,9))/3.d0
c
          write (6,'(10x,a,/)') ' (Hyper)magnetizabilities'
          write (6,'(10x,a,/)') 'eta(ab,gd): ab=BaBb gd=FgFd'
          write (6,'(10x,a,i,/)') '...Point group...', ipgr
          write (6,'(10x,a,g20.5,/)') '...Anysotropy...    ',deta
          write (6,'(10x,a,g20.5)') '...Chi(xx)     ',cmag(ifirst,1)
          write (6,'(10x,a,g20.5)') '...Chi(yy)     ',cmag(ifirst,5)
          write (6,'(10x,a,g20.5)') '...Chi(zz)     ',cmag(ifirst,9)
          write (6,'(10x,a,g20.5,/)') '...Chi(av)     ',cmagav
          write (6,'(10x,a,g20.5)') '...ETA(xxxx)     ',exxxx
          write (6,'(10x,a,g20.5)') '...ETA(yyxx)     ',eyyxx
          write (6,'(10x,a,g20.5,/)') '...ETA(zzxx)     ',ezzxx
          write (6,'(10x,a,g20.5)') '...ETA(xxyy)     ',exxyy
          write (6,'(10x,a,g20.5)') '...ETA(yyyy)     ',eyyyy
          write (6,'(10x,a,g20.5,/)') '...ETA(zzyy)     ',ezzyy
          write (6,'(10x,a,g20.5)') '...ETA(xxzz)     ',exxzz
          write (6,'(10x,a,g20.5)') '...ETA(yyzz)     ',eyyzz
          write (6,'(10x,a,g20.5,/)') '...ETA(zzzz)     ',ezzzz
          write (6,'(10x,a,g20.5)') '...ETA(xyxy)     ',exyxy
          write (6,'(10x,a,g20.5)') '...ETA(xzxz)     ',exzxz
          write (6,'(10x,a,g20.5,/)') '...ETA(yzyz)     ',eyzyz
c
        end if
c
        if (ish.gt.0) then
         if(issym.eq.ipgr) then
          do 319 iatom=nist,nlast,nstep
c
            Sxxxx=2.d0*(cshi(ithird,1,iatom)-cshi(ifirst,1,iatom))/
     #                 (fie(ithird,1)*fie(ithird,1))
c
            Syyxx=2.d0*(cshi(ithird,5,iatom)-cshi(ifirst,5,iatom))/
     #                 (fie(ithird,1)*fie(ithird,1))
c
            Szzxx=2.d0*(cshi(ithird,9,iatom)-cshi(ifirst,9,iatom))/
     #                 (fie(ithird,1)*fie(ithird,1))
c
            Sxxyy=2.d0*(cshi(ifourth,1,iatom)-cshi(ifirst,1,iatom))/
     #                 (fie(ifourth,2)*fie(ifourth,2))
c
            Syyyy=2.d0*(cshi(ifourth,5,iatom)-cshi(ifirst,5,iatom))/
     #                 (fie(ifourth,2)*fie(ifourth,2))
c
            Szzyy=2.d0*(cshi(ifourth,9,iatom)-cshi(ifirst,9,iatom))/
     #                 (fie(ifourth,2)*fie(ifourth,2))
c
            Szzzz=2.d0*(cshi(isecnd,9,iatom)-cshi(ifirst,9,iatom))/
     #                 (fie(isecnd,3)*fie(isecnd,3))
c
            Sxxzz=2.d0*(cshi(isecnd,1,iatom)-cshi(ifirst,1,iatom))/
     #                 (fie(isecnd,3)*fie(isecnd,3))
c
            Syyzz=2.d0*(cshi(isecnd,5,iatom)-cshi(ifirst,5,iatom))/
     #                 (fie(isecnd,3)*fie(isecnd,3))
c
            Sxyxy=(cshi(ififth,2,iatom)-cshi(ifirst,2,iatom))/
     #                 (fie(ififth,1)*fie(ififth,2))
c
            Syxxy=(cshi(ififth,4,iatom)-cshi(ifirst,4,iatom))/
     #                 (fie(ififth,1)*fie(ififth,2))
c
            Sxzxz=(cshi(isixth,3,iatom)-cshi(ifirst,3,iatom))/
     #                 (fie(isixth,1)*fie(isixth,3))
c
            Szxxz=(cshi(isixth,7,iatom)-cshi(ifirst,7,iatom))/
     #                 (fie(isixth,1)*fie(isixth,3))
c
            Syzyz=(cshi(isevth,6,iatom)-cshi(ifirst,6,iatom))/
     #                 (fie(isevth,2)*fie(isevth,3))
c
            Szyyz=(cshi(isevth,8,iatom)-cshi(ifirst,8,iatom))/
     #                 (fie(isevth,2)*fie(isevth,3))
c      
            Bzz=-(Sxxzz+Syyzz+Szzzz)/6.d0
            Byy=-(Sxxyy+Syyyy+Szzyy)/6.d0
            Bxx=-(Sxxxx+Syyxx+Szzxx)/6.d0
c
            sigave=(cshi(ifirst,1,iatom)+cshi(ifirst,5,iatom)+
     #               cshi(ifirst,9,iatom))/3.d0
c
             write (6,'(//10x,a,a,/)')
     # ' Shielding (Polarizabilities) for atom ',satom(iatom)
             write (6,'(10x,a,i,/)') 'Site Symmetry  ', issym
             write (6,'(10x,a,g20.5)')'SIGMA(xx) ',cshi(ifirst,1,iatom)
             write (6,'(10x,a,g20.5)')'SIGMA(yy) ',cshi(ifirst,5,iatom)
             write (6,'(10x,a,g20.5)')'SIGMA(zz) ',cshi(ifirst,9,iatom)
             write (6,'(10x,a,g20.5,/)')'SIGMA(av.)     ',sigave
             write (6,'(10x,a,g20.5)') '...SIGMA(xxz)      ',sxxz
             write (6,'(10x,a,g20.5)') '...SIGMA(xzx)      ',sxzx
             write (6,'(10x,a,g20.5)') '...SIGMA(zxx)      ',szxx
             write (6,'(10x,a,g20.5)') '...SIGMA(yyz)      ',syyz
             write (6,'(10x,a,g20.5)') '...SIGMA(yzy)      ',syzy
             write (6,'(10x,a,g20.5)') '...SIGMA(zyy)      ',szyy
             write (6,'(10x,a,g20.5,/)') '...SIGMA(zzz)      ',szzz
             write (6,'(10x,a,g20.5)') '...SIGMA(xxxx)     ',sxxxx
             write (6,'(10x,a,g20.5)') '...SIGMA(yyxx)     ',syyxx
             write (6,'(10x,a,g20.5)') '...SIGMA(zzxx)     ',szzxx
             write (6,'(10x,a,g20.5)') '...SIGMA(xxyy)     ',sxxyy
             write (6,'(10x,a,g20.5)') '...SIGMA(yyyy)     ',syyyy
             write (6,'(10x,a,g20.5)') '...SIGMA(zzyy)     ',szzyy
             write (6,'(10x,a,g20.5)') '...SIGMA(xxzz)     ',sxxzz
             write (6,'(10x,a,g20.5)') '...SIGMA(yyzz)     ',syyzz
             write (6,'(10x,a,g20.5)') '...SIGMA(zzzz)     ',szzzz
             write (6,'(10x,a,g20.5)') '...SIGMA(xyxy)     ',sxyxy
             write (6,'(10x,a,g20.5)') '...SIGMA(yxxy)     ',syxxy
             write (6,'(10x,a,g20.5)') '...SIGMA(xzxz)     ',sxzxz
             write (6,'(10x,a,g20.5)') '...SIGMA(zxxz)     ',szxxz
             write (6,'(10x,a,g20.5)') '...SIGMA(yzyz)     ',syzyz
             write (6,'(10x,a,g20.5,/)') '...SIGMA(zyyz)     ',szyyz
             write (6,'(10x,a,g20.5)') '...Az     ',Az
             write (6,'(10x,a,g20.5)') '...Bxx     ',Bxx
             write (6,'(10x,a,g20.5)') '...Byy     ',Byy
             write (6,'(10x,a,g20.5)') '...Bzz     ',Bzz
  319     continue
         else
          ipgr=issym
          imag=1
          goto 91
         end if
        end if
      else
        write(6,*)'Insufficient # of fields for required symmetry !!!'
        write(6,*)'Check value set in readmg.dat'
      end if
      go to 9999
c -------------------------------------------
c  Case C2v, No symmetry - 8 calculations)
c
c     Need     0  0  0
c     Need     0  0  z
c     Need     0  0 -z  
c     Need     x  0  0
c     Need     0  y  0
c     Need     x  y  0
c     Need     x  0  z
c     Need     0  y  z
c
  409 case1=ifound(14).eq.1                      
      case2=ifound(13).eq.1.and.ifound(15).eq.1 
      case3=ifound(23).eq.1.or.ifound(5).eq.1   
      case4=ifound(11).eq.1.or.ifound(17).eq.1  
      case5=ifound(2).eq.1.or.ifound(8).eq.1.or.
     #       ifound(20).eq.1.or.ifound(26).eq.1
      case6=ifound(4).eq.1.or.ifound(6).eq.1.or. 
     #       ifound(22).eq.1.or.ifound(24).eq.1
      case7=ifound(18).eq.1.or.ifound(10).eq.1
     #      .or.ifound(12).eq.1.or.ifound(16).eq.1 
      if (case1.and.case2.and.case3.and.case4
     #   .and.case5.and.case6.and.case7) then     
        ifirst=14                                   
        isecnd1=15                                 
        isecnd2=13                                
        ithird=5                                 
        if (ifound(23).eq.1) ithird=23          
        ifourth=11                             
        if (ifound(17).eq.1) ifourth=17       
        ififth=26                            
        if (ifound(2).eq.1) ififth=2        
        if (ifound(8).eq.1) ififth=8       
        if (ifound(20).eq.1) ififth=20    
        isixth=24                        
        if (ifound(4).eq.1) isixth=4    
        if (ifound(6).eq.1) isixth=6   
        if (ifound(22).eq.1) isixth=22
        isevth=18                     
        if (ifound(10).eq.1) isevth=10
        if (ifound(12).eq.1) isevth=12
        if (ifound(16).eq.1) isevth=16
c
        if (imag.eq.0) then
c
          cxzx=(cmag(ithird,3)-cmag(ifirst,3))/fie(ithird,1)
c      
          cyzy=(cmag(ifourth,6)-cmag(ifirst,6))/fie(ifourth,2)
c      
          czzz=
     #((cmag(isecnd1,9)-cmag(ifirst,9))*fie(isecnd2,3)*fie(isecnd2,3)-
     # (cmag(isecnd2,9)-cmag(ifirst,9))*fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c      
          cxxz=
     #((cmag(isecnd1,1)-cmag(ifirst,1))*fie(isecnd2,3)*fie(isecnd2,3)-
     # (cmag(isecnd2,1)-cmag(ifirst,1))*fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
          cyyz=
     #((cmag(isecnd1,5)-cmag(ifirst,5))*fie(isecnd2,3)*fie(isecnd2,3)-
     # (cmag(isecnd2,5)-cmag(ifirst,5))*fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
          exxxx=2.d0*(cmag(ithird,1)-cmag(ifirst,1))/
     #           (fie(ithird,1)*fie(ithird,1))
c      
          eyyxx=2.d0*(cmag(ithird,5)-cmag(ifirst,5))/
     #           (fie(ithird,1)*fie(ithird,1))
c      
          ezzxx=2.d0*(cmag(ithird,9)-cmag(ifirst,9))/
     #           (fie(ithird,1)*fie(ithird,1))
c
          exxyy=2.d0*(cmag(ifourth,1)-cmag(ifirst,1))/
     #           (fie(ifourth,2)*fie(ifourth,2))
c      
          eyyyy=2.d0*(cmag(ifourth,5)-cmag(ifirst,5))/
     #           (fie(ifourth,2)*fie(ifourth,2))
c      
          ezzyy=2.d0*(cmag(ifourth,9)-cmag(ifirst,9))/
     #           (fie(ifourth,2)*fie(ifourth,2))
c      
          exxzz=2.d0*((cmag(isecnd2,1)-cmag(ifirst,1))*fie(isecnd1,3)-
     # (cmag(isecnd1,1)-cmag(ifirst,1))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c      
          eyyzz=2.d0*((cmag(isecnd2,5)-cmag(ifirst,5))*fie(isecnd1,3)-
     # (cmag(isecnd1,5)-cmag(ifirst,5))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c     
          ezzzz=2.d0*((cmag(isecnd2,9)-cmag(ifirst,9))*fie(isecnd1,3)-
     # (cmag(isecnd1,9)-cmag(ifirst,9))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c     
          exyxy=(cmag(ififth,2)-cmag(ifirst,2)-Cxyx*fie(ififth,1))/
     #           (fie(ififth,1)*fie(ififth,2))
c      
          exzxz=(cmag(isixth,3)-cmag(ifirst,3)-Cxzx*fie(isixth,1))/
     #           (fie(isixth,1)*fie(isixth,3))
c      
          eyzyz=(cmag(isevth,6)-cmag(ifirst,6)-Cyzy*fie(isevth,2))/
     #           (fie(isevth,2)*fie(isevth,3))
c      
          deta=(2.d0*exxxx+2.d0*eyyyy+2.d0*ezzzz
     #         +6.d0*exyxy+6.d0*exzxz+6.d0*eyzyz
     #         -exxyy-eyyxx-exxzz-ezzxx-eyyzz-ezzyy)/15.d0
c
          cmagav=(cmag(ifirst,1)+cmag(ifirst,5)+cmag(ifirst,9))/3.d0
          dcxxyy=cmag(ifirst,1)-cmag(ifirst,5)
          dcyyxx=cmag(ifirst,5)-cmag(ifirst,9)
          dcxxzz=cmag(ifirst,1)-cmag(ifirst,9)
c
          write (6,'(10x,a,/)') ' (Hyper)magnetizabilities'
          write (6,'(10x,a,/)') 'eta(ab,gd): ab=BaBb gd=FgFd'
          write (6,'(10x,a,i,/)') '...Point group...', ipgr
          write (6,'(10x,a,g20.5,/)') '...Anysotropy...    ',deta
          write (6,'(10x,a,g20.5)') '...Chi(xx)     ',cmag(ifirst,1)
          write (6,'(10x,a,g20.5)') '...Chi(yy)     ',cmag(ifirst,5)
          write (6,'(10x,a,g20.5)') '...Chi(zz)     ',cmag(ifirst,9)
          write (6,'(10x,a,g20.5,/)') '...Chi(av.)     ',cmagav
          write (6,'(10x,a,g20.5)') '...DeltaChi(xx,yy)',dcxxyy
          write (6,'(10x,a,g20.5)') '...DeltaChi(yy,xx)',dcyyxx
          write (6,'(10x,a,g20.5,/)') '...DeltaChi(xx,zz)',dcxxzz
          write (6,'(10x,a,g20.5)') '...Csi(xxz)     ',Cxxz
          write (6,'(10x,a,g20.5)') '...Csi(yyz)     ',Cyyz
          write (6,'(10x,a,g20.5)') '...Csi(zzz)     ',Czzz
          write (6,'(10x,a,g20.5)') '...Csi(xzx)     ',Cxzx
          write (6,'(10x,a,g20.5,/)') '...Csi(yzy)     ',Cyzy
          write (6,'(10x,a,g20.5)') '...ETA(xxxx)     ',exxxx
          write (6,'(10x,a,g20.5)') '...ETA(yyxx)     ',eyyxx
          write (6,'(10x,a,g20.5,/)') '...ETA(zzxx)     ',ezzxx
          write (6,'(10x,a,g20.5)') '...ETA(xxyy)     ',exxyy
          write (6,'(10x,a,g20.5)') '...ETA(yyyy)     ',eyyyy
          write (6,'(10x,a,g20.5,/)') '...ETA(zzyy)     ',ezzyy
          write (6,'(10x,a,g20.5)') '...ETA(xxzz)     ',exxzz
          write (6,'(10x,a,g20.5)') '...ETA(yyzz)     ',eyyzz
          write (6,'(10x,a,g20.5,/)') '...ETA(zzzz)     ',ezzzz
          write (6,'(10x,a,g20.5)') '...ETA(xyxy)     ',exyxy
          write (6,'(10x,a,g20.5)') '...ETA(xzxz)     ',exzxz
          write (6,'(10x,a,g20.5,/)') '...ETA(yzyz)     ',eyzyz
        end if
c
        if (ish.gt.0) then
         if(issym.eq.ipgr) then
          do 419 iatom=nist,nlast,nstep
c
             Szxx=(cshi(ithird,7,iatom)-cshi(ifirst,7,iatom))/
     #             fie(ithird,1)
c
             Sxzx=(cshi(ithird,3,iatom)-cshi(ifirst,3,iatom))/
     #             fie(ithird,1)
c
             Szyy=(cshi(ifourth,8,iatom)-cshi(ifirst,8,iatom))/
     #             fie(ifourth,2)
c      
             Syzy=(cshi(ifourth,6,iatom)-cshi(ifirst,6,iatom))/
     #             fie(ifourth,2)
c
             Szzz=
     # ((cshi(isecnd1,9,iatom)-cshi(ifirst,9,iatom))*
     # fie(isecnd2,3)*fie(isecnd2,3)-
     # (cshi(isecnd2,9,iatom)-cshi(ifirst,9,iatom))*
     # fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c      
             Sxxz=
     # ((cshi(isecnd1,1,iatom)-cshi(ifirst,1,iatom))*
     # fie(isecnd2,3)*fie(isecnd2,3)-
     # (cshi(isecnd2,1,iatom)-cshi(ifirst,1,iatom))*
     # fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c      
             Syyz=
     # ((cshi(isecnd1,5,iatom)-cshi(ifirst,5,iatom))*
     # fie(isecnd2,3)*fie(isecnd2,3)-
     # (cshi(isecnd2,5,iatom)-cshi(ifirst,5,iatom))*
     # fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
             Sxxxx=2.d0*(cshi(ithird,1,iatom)-cshi(ifirst,1,iatom))/
     #           (fie(ithird,1)*fie(ithird,1))
c
             Syyxx=2.d0*(cshi(ithird,5,iatom)-cshi(ifirst,5,iatom))/
     #           (fie(ithird,1)*fie(ithird,1))
c
             Szzxx=2.d0*(cshi(ithird,9,iatom)-cshi(ifirst,9,iatom))/
     #           (fie(ithird,1)*fie(ithird,1))
c
             Sxxyy=2.d0*(cshi(ifourth,1,iatom)-cshi(ifirst,1,iatom))/
     #          (fie(ifourth,2)*fie(ifourth,2))
c
             Syyyy=2.d0*(cshi(ifourth,5,iatom)-cshi(ifirst,5,iatom))/
     #           (fie(ifourth,2)*fie(ifourth,2))
c
             Szzyy=2.d0*(cshi(ifourth,9,iatom)-cshi(ifirst,9,iatom))/
     #           (fie(ifourth,2)*fie(ifourth,2))
c
             Szzzz=2.d0*
     # ((cshi(isecnd2,9,iatom)-cshi(ifirst,9,iatom))*fie(isecnd1,3)-
     # (cshi(isecnd1,9,iatom)-cshi(ifirst,9,iatom))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c      
             Sxxzz=2.d0*
     # ((cshi(isecnd2,1,iatom)-cshi(ifirst,1,iatom))*fie(isecnd1,3)-
     # (cshi(isecnd1,1,iatom)-cshi(ifirst,1,iatom))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c      
             Syyzz=2.d0*
     # ((cshi(isecnd2,5,iatom)-cshi(ifirst,5,iatom))*fie(isecnd1,3)-
     # (cshi(isecnd1,5,iatom)-cshi(ifirst,5,iatom))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
             Sxyxy=(cshi(ififth,2,iatom)-cshi(ifirst,2,iatom))/
     #           (fie(ififth,1)*fie(ififth,2))
c
             Syxxy=(cshi(ififth,4,iatom)-cshi(ifirst,4,iatom))/
     #           (fie(ififth,1)*fie(ififth,2))
c      
             Sxzxz=(cshi(isixth,3,iatom)-cshi(ifirst,3,iatom)
     #           -sxzx*fie(isixth,1))/(fie(isixth,1)*fie(isixth,3))
c      
             Szxxz=(cshi(isixth,7,iatom)-cshi(ifirst,7,iatom)
     #           -szxx*fie(isixth,1))/(fie(isixth,1)*fie(isixth,3))
c
             Syzyz=(cshi(isevth,6,iatom)-cshi(ifirst,6,iatom)
     #           -syzy*fie(isevth,2))/(fie(isevth,2)*fie(isevth,3))
c      
             Szyyz=(cshi(isevth,8,iatom)-cshi(ifirst,8,iatom)
     #           -szyy*fie(isevth,2))/(fie(isevth,2)*fie(isevth,3))
c      
             Az=-(Sxxz+Syyz+Szzz)/3.d0
             Bzz=-(Sxxzz+Syyzz+Szzzz)/6.d0
             Byy=-(Sxxyy+Syyyy+Szzyy)/6.d0
             Bxx=-(Sxxxx+Syyxx+Szzxx)/6.d0
c
             sigave=(cshi(ifirst,1,iatom)+cshi(ifirst,5,iatom)+
     #               cshi(ifirst,9,iatom))/3.d0
c
             write (6,'(//10x,a,a,/)')
     # ' Shielding (Polarizabilities) for atom ',satom(iatom)
             write (6,'(10x,a,i,/)') 'Site Symmetry   ',  issym
             write (6,'(10x,a,g20.5)')'SIGMA(xx)  ',cshi(ifirst,1,iatom)
             write (6,'(10x,a,g20.5)')'SIGMA(yy)  ',cshi(ifirst,5,iatom)
             write (6,'(10x,a,g20.5)')'SIGMA(zz)  ',cshi(ifirst,9,iatom)
             write (6,'(10x,a,g20.5,/)') '...SIGMA(av.)     ',sigave
             write (6,'(10x,a,g20.5)') '...SIGMA(xxz)      ',sxxz
             write (6,'(10x,a,g20.5)') '...SIGMA(xzx)      ',sxzx
             write (6,'(10x,a,g20.5)') '...SIGMA(zxx)      ',szxx
             write (6,'(10x,a,g20.5)') '...SIGMA(yyz)      ',syyz
             write (6,'(10x,a,g20.5)') '...SIGMA(yzy)      ',syzy
             write (6,'(10x,a,g20.5)') '...SIGMA(zyy)      ',szyy
             write (6,'(10x,a,g20.5,/)') '...SIGMA(zzz)      ',szzz
             write (6,'(10x,a,g20.5)') '...SIGMA(xxxx)     ',sxxxx
             write (6,'(10x,a,g20.5)') '...SIGMA(yyxx)     ',syyxx
             write (6,'(10x,a,g20.5)') '...SIGMA(zzxx)     ',szzxx
             write (6,'(10x,a,g20.5)') '...SIGMA(xxyy)     ',sxxyy
             write (6,'(10x,a,g20.5)') '...SIGMA(yyyy)     ',syyyy
             write (6,'(10x,a,g20.5)') '...SIGMA(zzyy)     ',szzyy
             write (6,'(10x,a,g20.5)') '...SIGMA(xxzz)     ',sxxzz
             write (6,'(10x,a,g20.5)') '...SIGMA(yyzz)     ',syyzz
             write (6,'(10x,a,g20.5)') '...SIGMA(zzzz)     ',szzzz
	     write (6,'(10x,a,g20.5)') '...SIGMA(xyxy)     ',sxyxy
	     write (6,'(10x,a,g20.5)') '...SIGMA(yxxy)     ',syxxy
             write (6,'(10x,a,g20.5)') '...SIGMA(xzxz)     ',sxzxz
             write (6,'(10x,a,g20.5)') '...SIGMA(zxxz)     ',szxxz
             write (6,'(10x,a,g20.5)') '...SIGMA(yzyz)     ',syzyz
             write (6,'(10x,a,g20.5,/)') '...SIGMA(zyyz)     ',szyyz
             write (6,'(10x,a,g20.5)') '...Az     ',Az
             write (6,'(10x,a,g20.5)') '...Bxx     ',Bxx
             write (6,'(10x,a,g20.5)') '...Byy     ',Byy
             write (6,'(10x,a,g20.5)') '...Bzz     ',Bzz
  419     continue 
         else
          ipgr=issym
          imag=1
          goto 91
         end if
        end if
      else
        write(6,*)'Insufficient # of fields for required symmetry !!!'
        write(6,*)'Check value set in readmg.dat'
      end if
      go to 9999
c-------------------------------------
c  Case C3v ( Shielding H1 in CH4)
c
c     Need     0  0  0
c     Need     0  0  z
c     Need     0  0 -z
c     Need     x  0  0
c     Need    -x  0  0
c     Need     x  0  z
c
  509 case1=ifound(14).eq.1                       
      case2=ifound(13).eq.1.and.ifound(15).eq.1  
      case3=ifound(23).eq.1.and.ifound(5).eq.1  
      case4=ifound(24).eq.1.or.ifound(22).eq.1.or.
     #      ifound(6).eq.1.or.ifound(4).eq.1
      if (case1.and.case2.and.case3.and.case4) then 
        ifirst=14
        isecnd1=15
        isecnd2=13
        ithird1=23
        ithird2=5 
        if (ifound(24).eq.1) ifourth=24
        if (ifound(22).eq.1) ifourth=22
        if (ifound(6).eq.1) ifourth=6
        if (ifound(4).eq.1) ifourth=4
        if (imagn.eq.0) then 
           write(6,*)'WARNING:  hyperpol for C3v not implemented'
        end if
        if (ish.gt.0) then
         if(issym.eq.ipgr) then
          do 519 iatom=nist,nlast,step
c
            Sxxx=
     # ((cshi(ithird1,1,iatom)-cshi(ifirst,1,iatom))*
     # fie(ithird2,1)*fie(ithird2,1)-
     # (cshi(ithird2,1,iatom)-cshi(ifirst,1,iatom))*
     # fie(ithird1,1)*fie(ithird1,1))/
     # (fie(ithird2,1)*fie(ithird1,1)*(fie(ithird2,1)-fie(ithird1,1)))
c      
            Sxxxx=2.d0*
     # ((cshi(ithird2,1,iatom)-cshi(ifirst,1,iatom))*fie(ithird1,1)-
     # (cshi(ithird1,1,iatom)-cshi(ifirst,1,iatom))*fie(ithird2,1))/
     # (fie(ithird2,1)*fie(ithird1,1)*(fie(ithird2,1)-fie(ithird1,1)))
c      
            Sxxyy=2.d0*(cshi(ithird1,5,iatom)-cshi(ifirst,5,iatom)+
     #         Sxxx*fie(ithird1,1))/(fie(ithird1,1)*fie(ithird1,1))
c      
            Szzxx=2.d0*(cshi(ithird1,9,iatom)-cshi(ifirst,9,iatom))/
     #           (fie(ithird1,1)*fie(ithird1,1))
c      
            Szxx=
     # ((cshi(ithird1,7,iatom)-cshi(ifirst,7,iatom))*
     # fie(ithird2,1)*fie(ithird2,1)-
     # (cshi(ithird2,7,iatom)-cshi(ifirst,7,iatom))*
     # fie(ithird1,1)*fie(ithird1,1))/
     # (fie(ithird2,1)*fie(ithird1,1)*(fie(ithird2,1)-fie(ithird1,1)))
c      
            Sxzx=
     # ((cshi(ithird1,3,iatom)-cshi(ifirst,3,iatom))*
     # fie(ithird2,1)*fie(ithird2,1)-
     # (cshi(ithird2,3,iatom)-cshi(ifirst,3,iatom))*
     # fie(ithird1,1)*fie(ithird1,1))/
     # (fie(ithird2,1)*fie(ithird1,1)*(fie(ithird2,1)-fie(ithird1,1)))
c
            Szzzz=
     #2.d0*((cshi(isecnd2,9,iatom)-cshi(ifirst,9,iatom))*fie(isecnd1,3)-
     # (cshi(isecnd1,9,iatom)-cshi(ifirst,9,iatom))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c      
            Sxxzz=2.d0*
     # ((cshi(isecnd2,1,iatom)-cshi(ifirst,1,iatom))*fie(isecnd1,3)-
     # (cshi(isecnd1,1,iatom)-cshi(ifirst,1,iatom))*fie(isecnd2,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c
            Szzz=
     # ((cshi(isecnd1,9,iatom)-cshi(ifirst,9,iatom))*
     # fie(isecnd2,3)*fie(isecnd2,3)-
     # (cshi(isecnd2,9,iatom)-cshi(ifirst,9,iatom))*
     # fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c      
            Sxxz=
     # ((cshi(isecnd1,1,iatom)-cshi(ifirst,1,iatom))*
     # fie(isecnd2,3)*fie(isecnd2,3)-
     # (cshi(isecnd2,1,iatom)-cshi(ifirst,1,iatom))*
     # fie(isecnd1,3)*fie(isecnd1,3))/
     # (fie(isecnd2,3)*fie(isecnd1,3)*(fie(isecnd2,3)-fie(isecnd1,3)))
c      
            Sxyxy=(Sxxxx-Sxxyy)/2.d0
c
            Sxzxx=2.d0*
     # ((cshi(ithird2,3,iatom)-cshi(ifirst,3,iatom))*fie(ithird1,1)-
     # (cshi(ithird1,3,iatom)-cshi(ifirst,3,iatom))*fie(ithird2,1))/
     # (fie(ithird2,1)*fie(ithird1,1)*(fie(ithird2,1)-fie(ithird1,1)))   
c
            Szxxx=
     #2.d0*((cshi(ithird2,9,iatom)-cshi(ifirst,9,iatom))*fie(ithird1,1)-    
     # (cshi(ithird1,9,iatom)-cshi(ifirst,9,iatom))*fie(ithird2,1))/
     # (fie(ithird2,1)*fie(ithird1,1)*(fie(ithird2,1)-fie(ithird1,1)))
c
            Sxzxz=
     #      2.d0*(cshi(ifourth,3,iatom)-cshi(ifirst,3,iatom)-sxzx*
     #      fie(ifourth,1)-5.d-1*Sxzxx*fie(ifourth,1)*fie(ifourth1,1))
     #      /(fie(ifourth,1)*fie(ifourth,3))
c      
            Szxxz=
     #      2.d0*(cshi(ifourth,9,iatom)-cshi(ifirst,9,iatom)-szxx*
     #      fie(ifourth,1)-5.d-1*Szxxx*fie(ifourth,1)*fie(ifourth1,1))
     #      /(fie(ifourth,1)*fie(ifourth,3))
c
            Az=-(2.d0*Sxxz+Szzz)/3.d0
            Bzz=-(2.d0*Sxxzz+Szzzz)/6.d0
c    
c     note that Bxx=Byy, Sxxyy=Syyxx
c
            Bxx=-(Sxxxx+Sxxyy+Szzxx)/6.d0
c
            sigave=(cshi(ifirst,1,iatom)+cshi(ifirst,5,iatom)+
     #               cshi(ifirst,9,iatom))/3.d0
c
            write (6,'(//10x,a,a,/)')
     # ' Shielding (Polarizabilities) for atom ',satom(iatom)
            write (6,'(10x,a,i,/)') 'Site Symmetry    ',issym
            write (6,'(10x,a,g20.5)') 'SIGMA(xx) ',cshi(ifirst,1,iatom)
            write (6,'(10x,a,g20.5)') 'SIGMA(yy) ',cshi(ifirst,5,iatom)
            write (6,'(10x,a,g20.5)') 'SIGMA(zz) ',cshi(ifirst,9,iatom)
            write (6,'(10x,a,g20.5,/)') '...SIGMA(av.)     ',sigave
            write (6,'(10x,a,g20.5)') '...SIGMA(xxx)      ',sxxx
            write (6,'(10x,a,g20.5)') '...SIGMA(xxz)      ',sxxz
            write (6,'(10x,a,g20.5)') '...SIGMA(xzx)      ',sxzx
            write (6,'(10x,a,g20.5)') '...SIGMA(zxx)      ',szxx
            write (6,'(10x,a,g20.5,/)') '...SIGMA(zzz)      ',szzz
            write (6,'(10x,a,g20.5)') '...SIGMA(xxxx)     ',sxxxx
            write (6,'(10x,a,g20.5)') '...SIGMA(zzxx)     ',szzxx
            write (6,'(10x,a,g20.5)') '...SIGMA(xxyy)     ',sxxyy
            write (6,'(10x,a,g20.5)') '...SIGMA(xxzz)     ',sxxzz
            write (6,'(10x,a,g20.5)') '...SIGMA(zzzz)     ',szzzz
            write (6,'(10x,a,g20.5)') '...SIGMA(xyxy)     ',sxyxy
            write (6,'(10x,a,g20.5)') '...SIGMA(xzxx)     ',sxzxx
            write (6,'(10x,a,g20.5)') '...SIGMA(zxxx)     ',szxxx
            write (6,'(10x,a,g20.5)') '...SIGMA(xzxz)     ',sxzxz
            write (6,'(10x,a,g20.5)') '...SIGMA(zxxz)     ',szxxz
            write (6,'(10x,a,g20.5)') '...Az     ',Az
            write (6,'(10x,a,g20.5)') '...Bxx     ',Bxx
            write (6,'(10x,a,g20.5)') '...Bzz     ',Bzz
  519     continue
         else
          ipgr=issym
          imag=1
          goto 91
         endif
        end if
      else
        write(6,*)'Insufficient # of fields for required symmetry !!!'
        write(6,*)'Check value set in readmg.dat'
      end if
      go to 9999
C --------------------------------
c Case Dih (N2  4 calculations)
c
c     Need     0  0  0
c     Need     0  0  z
c     Need     x  0  0
c     Need     x  0  z
c
  609 case1=ifound(14).eq.1                         
      case2=ifound(5).eq.1.or.ifound(23).eq.1      
      case3=ifound(15).eq.1.or.ifound(13).eq.1     
      case4=ifound(24).eq.1.or.ifound(22).eq.1.or. 
     #      ifound(6).eq.1.or.ifound(4).eq.1
      if (case1.and.case2.and.case3.and.case4) then
        ifirst=14                                 
        isecnd=5                                 
        if (ifound(23).eq.1) isecnd=23          
        ithird=13                              
        if (ifound(15).eq.1) ithird=15        
        ifourth=4                            
        if (ifound(24).eq.1) ifourth=24     
        if (ifound(22).eq.1) ifourth=22    
        if (ifound(6).eq.1) ifourth=6     
c      
        if (imag.eq.0) then  
c
          Ezzzz=2.d0*(cmag(ithird,9)-cmag(ifirst,9))/
     #           (fie(ithird,3)*fie(ithird,3))
          Exxxx=2.d0*(cmag(isecnd,1)-cmag(ifirst,1))/
     #           (fie(isecnd,1)*fie(isecnd,1))
          Exxzz=2.d0*(cmag(ithird,1)-cmag(ifirst,1))/
     #           (fie(ithird,3)*fie(ithird,3))
          Ezzxx=2.d0*(cmag(isecnd,9)-cmag(ifirst,9))/
     #           (fie(isecnd,1)*fie(isecnd,1))
          Eyyxx=2.d0*(cmag(isecnd,5)-cmag(ifirst,5))/
     #           (fie(isecnd,1)*fie(isecnd,1))
          Exzxz=(cmag(ifourth,3)-cmag(ifirst,3))/
     #           (fie(ifourth,1)*fie(ifourth,3))
c      
          Deta=(Exxxx*7.d0-Eyyxx*5.d0+Ezzzz*2.d0-Exxzz*2.d0
     #         -Ezzxx*2.d0+Exzxz*12.d0)/15.d0
c
          write (6,'(10x,a,/)') ' (Hyper)magnetizabilities'
          write (6,'(10x,a,/)') 'eta(ab,gd): ab=BaBb gd=FgFd'
          write (6,'(10x,a,i,/)') '...Point group...', ipgr
          write (6,'(10x,a,g20.5,/)') '...Anysotropy...    ',deta
          write (6,'(10x,a,g20.5)') '...Chi(xx)   ',cmag(ifirst,1)
          write (6,'(10x,a,g20.5,/)') '...Chi(zz)   ',cmag(ifirst,9)
          write (6,'(10x,a,g20.5)') '...ETA(zzzz)     ',ezzzz
          write (6,'(10x,a,g20.5)') '...ETA(xxxx)     ',exxxx
          write (6,'(10x,a,g20.5)') '...ETA(xxzz)     ',exxzz
          write (6,'(10x,a,g20.5)') '...ETA(yyxx)     ',eyyxx
          write (6,'(10x,a,g20.5)') '...ETA(zzxx)     ',ezzxx
          write (6,'(10x,a,g20.5,/)') '...ETA(xzxz)     ',exzxz
        end if                                       
        if (ish.gt.0) then
         if (issym.eq.ipgr) then
          write(6,*)' Case Dih not implemeted for shieldpol'
         else
          ipgr=issym
          imag=1
          goto 91
         endif
        end if
      else
        write(6,*)'Insufficient # of fields for required symmetry !!!'
        write(6,*)'Check value set in readmg.dat'
      end if
      go to 9999
c
 9998 stop
 9997 stop '...ONLY TITLE IN INPUT'
 9996 stop '...NOT ENOUGH FILENAMES IN INPUT'
      end
c--------------------------------------------------------
      SUBROUTINE GLENG (IBEGIN,STRINGA,FORM,ISFORM,IEND)
c
      implicit double precision (a-h,o-z)
      character*80 stringa
      character*31 form
      character*1 f1
c
      f1='F'
      i=ibegin
    1 i=i+1
      if (stringa(i:i).eq.' ') go to 1
      init=i
      iwhat=init-ibegin-1
      if (iwhat.eq.0) form(isform:isform+3)='    '
      if (iwhat.lt.10.and.iwhat.gt.0) 
     #                 write (form(isform+1:isform+1),'(I1)') iwhat
      if (iwhat.ge.10) write (form(isform:isform+1),'(I2)') iwhat
    2 i=i+1
      if (stringa(i:i).eq.'d'.or.stringa(i:i).eq.'D') f1='D'
      if (stringa(i:i).eq.'e'.or.stringa(i:i).eq.'E') f1='E'
      if (stringa(i:i).eq.'g'.or.stringa(i:i).eq.'G') f1='G'
      if (stringa(i:i).eq.' ') istop=i-1
      if (stringa(i:i).ne.' ') go to 2
      if (istop-init+1.lt.10) then
      form(isform+4:isform+4)=' '
      write (form(isform+5:isform+6),'(a1,i1)') f1,istop-init+1
      else
      write (form(isform+4:isform+6),'(a1,i2)') f1,istop-init+1
      end if
      iend=i-1
      return
      end
c--------------------------------------------------------------
      SUBROUTINE RDFILE (iunit,filenam,ish,satom,chi,shield,
     1                   field,energy,zprt,imag)
      implicit double precision (a-h,o-z)
      logical*1 zprt
      character*2 satom
      character*80 filenam
      dimension satom(*),chi(*),shield(9,10),field(*)
c
      energy=0.d0
      do jj=1,3
       field(jj)=0.d0
      end do
      do jj=1,9
       chi(jj)=0.d0
      end do
      do jj1=1,9
       do jj2=1,10
        shield(jj1,jj2)=0.d0
       end do
      end do
C
C Open file
C     
      open (unit=iunit,file=filenam,status='old')
C
C Read Fields
C      
      call readmg (iunit,'XDIPLEN',len('XDIPLEN'),value,iretx)
      if (iretx.eq.0) field(1)=value
      call readmg (iunit,'YDIPLEN',len('YDIPLEN'),value,irety)
      if (irety.eq.0) field(2)=value
      call readmg (iunit,'ZDIPLEN',len('ZDIPLEN'),value,iretz)
      if (iretz.eq.0) field(3)=value
      if (iretx+irety+iretz.eq.3) then
        write (6,'(////10x,A)') '...   ABORTING RUN  ...'
        write (6,'(1x,a)') filenam
        write (6,'(10x,a)') '.........  NO FIELDS IN INPUT'
        stop
      end if
C
C Read Energy
C      
      call readmg(iunit,'Total energy',len('Total energy'),value,iret)
      if (iret.eq.0) energy=value
C
C Read Magnetizability tensor
C
      call readte (iunit,chi,imag)
C
C Find how many shieldings
C
      call sumshi (iunit,ish,satom)
      if (ish.gt.0) call rdshi (iunit,ish,shield)
C
C Close file
C     
      close (unit=iunit)
c
      if (zprt) then
        write (6,'(3x,3a)') '...Data read from file ',filenam,':'
        write (6,*) '   ...X FIELD: ',field(1)
        write (6,*) '   ...Y FIELD: ',field(2)
        write (6,*) '   ...Z FIELD: ',field(3)
        write (6,*) '   ...Energy: ',energy
        if (imag.eq.0) then
          write (6,*) '   ...Magnetizability tensor:'
          write (6,*)  chi(1),chi(2),chi(3)
          write (6,*)  chi(4),chi(5),chi(6)
          write (6,*)  chi(7),chi(8),chi(9)
        end if
        if (ish.gt.0) then
          write (6,'(3x,a,i2,a)') '..We read ',ish,' Shield. tensor..'
          do j=1,ish
            write (6,'(3x,2a)') '...Shielding tensor for ',satom(j)
            write (6,*)  shield(1,j),shield(2,j),shield(3,j)
            write (6,*)  shield(4,j),shield(5,j),shield(6,j)
            write (6,*)  shield(7,j),shield(8,j),shield(9,j)
          end do
        end if
      end if
c
      return
      end
c-----------------------------------------------------------------
      SUBROUTINE RDSHI (IUNIT,ISHI,SHIELD)
c
      implicit double precision (a-h,o-z)   
      character*80 string
      dimension shield(9,10)
      rewind (unit=iunit)
      do j=1,ishi
    1   read (iunit,'(a80)',end=9) string
        ll=index(string,'Total shielding tensor (ppm):')
        if (ll.eq.0) go to 1
        read (iunit,'(a80)',end=9) string
        string(1:8)='        '
        call threere (string,shield(1,j),shield(2,j),shield(3,j))
        read (iunit,'(a80)',end=9) string
        string(1:8)='        '
        call threere (string,shield(4,j),shield(5,j),shield(6,j))
        read (iunit,'(a80)',end=9) string
        string(1:8)='        '
        call threere (string,shield(7,j),shield(8,j),shield(9,j))
      end do
    9 continue
      return 
      end
c --------------------------------------------------------------
      SUBROUTINE READMG (IUNIT,WHAT,LLENG,VALUE,IRET)
c
      implicit double precision (a-h,o-z)
      character*80 string
      character*7 what,sfor
      rewind (unit=iunit)
    1 read (iunit,'(a80)',end=9) string
      ll=index(string,what)
      if (ll.eq.0) go to 1
      do i=ll,ll+lleng
        string(i:i)=' '
      end do
      call snumber (string,sfor)
      read (string,sfor) value
      iret=0
      return
    9 iret=1
      return
      end
c  ---------------------------------------
      SUBROUTINE READTE (IUNIT,CHI,IRET)
c
      implicit double precision (a-h,o-z)
      dimension chi(9)
      character*80 string
      iret=0
      rewind (unit=iunit)
    1 read (iunit,'(a80)',end=9) string
      ll=index(string,'magneti')
      if (ll.eq.0) go to 1
      read (iunit,'(a80)',end=9) string
      string(1:4)='    '
      call threere (string,chi(1),chi(2),chi(3))
      read (iunit,'(a80)',end=9) string
      string(1:4)='    '
      call threere (string,chi(4),chi(5),chi(6))
      read (iunit,'(a80)',end=9) string
      string(1:4)='    '
      call threere (string,chi(7),chi(8),chi(9))
      return
    9 iret=1
      return
      end
c ------------------------------------------------
      SUBROUTINE SNUMBER (STRING,SFOR)
c
      implicit double precision (a-h,o-z)
      character*80 string
      character*7 sfor
      sfor='(g  .0)'
      i=0
    1 i=i+1
      if (string(i:i).eq.' ') go to 1
      istart=i
    2 i=i+1
      if (string(i:i).ne.' ') go to 2
      istop=i
      ilen=istop-istart+1
      write (sfor(3:4),'(i2)') ilen     
      i=0
    3 i=i+1
      string(i:i)=string(istart+i-1:istart+i-1)
      if (istart+i-1.lt.80) go to 3
      do j=i+1,80
        string(j:j)=' '
      end do
      return
      end
c ----------------------------------------
      SUBROUTINE SUMSHI (IUNIT,ISHI,SATOM)
c
      implicit double precision (a-h,o-z)   
      character*80 string
      character*2 satom
      dimension satom(10)
      ishi=0
      rewind (unit=iunit)
    1 read (iunit,'(a80)',end=9) string
      ll=index(string,'Chemical shielding for')
      if (ll.eq.0) go to 1
      ishi=ishi+1
      satom(ishi)=string(24:25)
      go to 1
    9 continue
      return 
      end
c ---------------------------------------------
      SUBROUTINE THREERE (STRINGA,A1,A2,A3)
c
      implicit double precision (a-h,o-z)
      character*80 stringa
      character*31 form
c
      form   ='(  X, cc.0,  X, ff.0,  X, ii.0)'
c      
      call gleng (0,stringa,form,2,iend)
      call gleng (iend,stringa,form,12,iend)
      call gleng (iend,stringa,form,22,iend)      
c
      read (stringa,form) a1,a2,a3
      return
      end
