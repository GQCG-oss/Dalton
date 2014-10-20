MODULE IchorEriCoulombintegralCPUOBSGeneralModSegSize
!Automatic Generated Code (AGC) by runOBSdriver.f90 in tools directory
!Contains routines for Segmented contracted Basisset 
use IchorCommonMod
use IchorEriCoulombintegralCPUMcMGeneralMod
  
private   
public :: ICI_CPU_OBS_general_sizeSeg  
  
CONTAINS
  subroutine ICI_CPU_OBS_general_sizeSeg(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP)
    implicit none
    integer,intent(inout) :: TMParray1maxsize,TMParray2maxsize
    integer,intent(in) :: AngmomA,AngmomB,AngmomC,AngmomD
    integer,intent(in) :: nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP
    integer,intent(in) :: nContA,nContB,nContC,nContD
    integer,intent(in) :: nPrimA,nPrimB,nPrimC,nPrimD
    ! local variables
    integer :: AngmomID
    
    AngmomID = 1000*AngmomA+100*AngmomB+10*AngmomC+AngmomD
    IF(UseGeneralCode) AngmomID = AngmomID + 10000 !force to use general code
    TMParray2maxSize = 1
    TMParray1maxSize = 1
    SELECT CASE(AngmomID)
    CASE(   0)  !Angmom(A= 0,B= 0,C= 0,D= 0) combi
    TMParray2maxSize = MAX(TMParray2maxSize,1)
    CASE(   1)  !Angmom(A= 0,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
    CASE(   2)  !Angmom(A= 0,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE(  10)  !Angmom(A= 0,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
    CASE(  11)  !Angmom(A= 0,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
    CASE(  12)  !Angmom(A= 0,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(  20)  !Angmom(A= 0,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE(  21)  !Angmom(A= 0,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(  22)  !Angmom(A= 0,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35)
       TMParray2maxSize = MAX(TMParray2maxSize,36)
    CASE( 100)  !Angmom(A= 0,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
    CASE( 101)  !Angmom(A= 0,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
    CASE( 102)  !Angmom(A= 0,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE( 110)  !Angmom(A= 0,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
    CASE( 111)  !Angmom(A= 0,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE( 112)  !Angmom(A= 0,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE( 120)  !Angmom(A= 0,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE( 121)  !Angmom(A= 0,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE( 122)  !Angmom(A= 0,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140)
       TMParray1maxSize = MAX(TMParray1maxSize,105)
       TMParray2maxSize = MAX(TMParray2maxSize,108)
    CASE( 200)  !Angmom(A= 0,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE( 201)  !Angmom(A= 0,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,24)
       TMParray2maxSize = MAX(TMParray2maxSize,20)
    CASE( 202)  !Angmom(A= 0,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE( 210)  !Angmom(A= 0,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,24)
       TMParray2maxSize = MAX(TMParray2maxSize,20)
    CASE( 211)  !Angmom(A= 0,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
    CASE( 212)  !Angmom(A= 0,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,120)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE( 220)  !Angmom(A= 0,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE( 221)  !Angmom(A= 0,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,120)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE( 222)  !Angmom(A= 0,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,210)
       TMParray2maxSize = MAX(TMParray2maxSize,175)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
    CASE(1000)  !Angmom(A= 1,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4)
    CASE(1001)  !Angmom(A= 1,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
    CASE(1002)  !Angmom(A= 1,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(1010)  !Angmom(A= 1,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,16)
       TMParray1maxSize = MAX(TMParray1maxSize,12)
    CASE(1011)  !Angmom(A= 1,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE(1012)  !Angmom(A= 1,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1020)  !Angmom(A= 1,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(1021)  !Angmom(A= 1,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1022)  !Angmom(A= 1,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140)
       TMParray1maxSize = MAX(TMParray1maxSize,105)
       TMParray2maxSize = MAX(TMParray2maxSize,108)
    CASE(1100)  !Angmom(A= 1,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
    CASE(1101)  !Angmom(A= 1,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,36)
    CASE(1102)  !Angmom(A= 1,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1110)  !Angmom(A= 1,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,36)
    CASE(1111)  !Angmom(A= 1,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(1112)  !Angmom(A= 1,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,162)
    CASE(1120)  !Angmom(A= 1,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
       TMParray2maxSize = MAX(TMParray2maxSize,54)
    CASE(1121)  !Angmom(A= 1,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,162)
    CASE(1122)  !Angmom(A= 1,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,315)
       TMParray2maxSize = MAX(TMParray2maxSize,324)
    CASE(1200)  !Angmom(A= 1,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(1201)  !Angmom(A= 1,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,72)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(1202)  !Angmom(A= 1,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(1210)  !Angmom(A= 1,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,72)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(1211)  !Angmom(A= 1,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
    CASE(1212)  !Angmom(A= 1,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
    CASE(1220)  !Angmom(A= 1,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(1221)  !Angmom(A= 1,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
    CASE(1222)  !Angmom(A= 1,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700)
       TMParray1maxSize = MAX(TMParray1maxSize,630)
       TMParray2maxSize = MAX(TMParray2maxSize,525)
       TMParray1maxSize = MAX(TMParray1maxSize,540)
    CASE(2000)  !Angmom(A= 2,B= 0,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,3*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,10)
       TMParray2maxSize = MAX(TMParray2maxSize,6)
    CASE(2001)  !Angmom(A= 2,B= 0,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,24)
       TMParray2maxSize = MAX(TMParray2maxSize,20)
    CASE(2002)  !Angmom(A= 2,B= 0,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE(2010)  !Angmom(A= 2,B= 0,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,40)
       TMParray1maxSize = MAX(TMParray1maxSize,24)
       TMParray2maxSize = MAX(TMParray2maxSize,20)
    CASE(2011)  !Angmom(A= 2,B= 0,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
    CASE(2012)  !Angmom(A= 2,B= 0,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,120)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2020)  !Angmom(A= 2,B= 0,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,60)
       TMParray2maxSize = MAX(TMParray2maxSize,50)
       TMParray1maxSize = MAX(TMParray1maxSize,30)
    CASE(2021)  !Angmom(A= 2,B= 0,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,120)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2022)  !Angmom(A= 2,B= 0,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,210)
       TMParray2maxSize = MAX(TMParray2maxSize,175)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
    CASE(2100)  !Angmom(A= 2,B= 1,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,4*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,20)
       TMParray2maxSize = MAX(TMParray2maxSize,18)
    CASE(2101)  !Angmom(A= 2,B= 1,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,72)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(2102)  !Angmom(A= 2,B= 1,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2110)  !Angmom(A= 2,B= 1,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,80)
       TMParray1maxSize = MAX(TMParray1maxSize,72)
       TMParray2maxSize = MAX(TMParray2maxSize,60)
    CASE(2111)  !Angmom(A= 2,B= 1,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
    CASE(2112)  !Angmom(A= 2,B= 1,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
    CASE(2120)  !Angmom(A= 2,B= 1,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,200)
       TMParray1maxSize = MAX(TMParray1maxSize,180)
       TMParray2maxSize = MAX(TMParray2maxSize,150)
       TMParray1maxSize = MAX(TMParray1maxSize,90)
    CASE(2121)  !Angmom(A= 2,B= 1,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,400)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,300)
       TMParray1maxSize = MAX(TMParray1maxSize,270)
    CASE(2122)  !Angmom(A= 2,B= 1,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700)
       TMParray1maxSize = MAX(TMParray1maxSize,630)
       TMParray2maxSize = MAX(TMParray2maxSize,525)
       TMParray1maxSize = MAX(TMParray1maxSize,540)
    CASE(2200)  !Angmom(A= 2,B= 2,C= 0,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,5*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,35)
       TMParray2maxSize = MAX(TMParray2maxSize,36)
    CASE(2201)  !Angmom(A= 2,B= 2,C= 0,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140)
       TMParray1maxSize = MAX(TMParray1maxSize,144)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
    CASE(2202)  !Angmom(A= 2,B= 2,C= 0,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,250)
       TMParray1maxSize = MAX(TMParray1maxSize,150)
    CASE(2210)  !Angmom(A= 2,B= 2,C= 1,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,6*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,56*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,140)
       TMParray1maxSize = MAX(TMParray1maxSize,144)
       TMParray2maxSize = MAX(TMParray2maxSize,100)
    CASE(2211)  !Angmom(A= 2,B= 2,C= 1,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,250)
    CASE(2212)  !Angmom(A= 2,B= 2,C= 1,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700)
       TMParray1maxSize = MAX(TMParray1maxSize,720)
       TMParray2maxSize = MAX(TMParray2maxSize,500)
       TMParray1maxSize = MAX(TMParray1maxSize,450)
    CASE(2220)  !Angmom(A= 2,B= 2,C= 2,D= 0) combi
       TMParray2maxSize = MAX(TMParray2maxSize,7*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,84*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,350)
       TMParray1maxSize = MAX(TMParray1maxSize,360)
       TMParray2maxSize = MAX(TMParray2maxSize,250)
       TMParray1maxSize = MAX(TMParray1maxSize,150)
    CASE(2221)  !Angmom(A= 2,B= 2,C= 2,D= 1) combi
       TMParray2maxSize = MAX(TMParray2maxSize,8*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,120*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,700)
       TMParray1maxSize = MAX(TMParray1maxSize,720)
       TMParray2maxSize = MAX(TMParray2maxSize,500)
       TMParray1maxSize = MAX(TMParray1maxSize,450)
    CASE(2222)  !Angmom(A= 2,B= 2,C= 2,D= 2) combi
       TMParray2maxSize = MAX(TMParray2maxSize,9*nPrimQP)
       TMParray1maxSize = MAX(TMParray1maxSize,165*nPrimQP)
       TMParray2maxSize = MAX(TMParray2maxSize,1225)
       TMParray1maxSize = MAX(TMParray1maxSize,1260)
       TMParray2maxSize = MAX(TMParray2maxSize,875)
       TMParray1maxSize = MAX(TMParray1maxSize,900)
    CASE DEFAULT
     call ICI_CPU_McM_general_size(TMParray1maxsize,&
         & TMParray2maxsize,AngmomA,AngmomB,AngmomC,AngmomD,&
         & nContA,nContB,nContC,nContD,&
         & nPrimA,nPrimB,nPrimC,nPrimD,&
         & nPrimP,nPrimQ,nContP,nContQ,nPrimQP,nContQP,&
         & .TRUE.,.TRUE.)
    END SELECT
  end subroutine ICI_CPU_OBS_general_sizeSeg
  
END MODULE IchorEriCoulombintegralCPUOBSGeneralModSegSize
