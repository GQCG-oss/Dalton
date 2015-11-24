!
!...   Copyright (c) 2015 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2015 (2015), see http://daltonprogram.org"
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
C   /* Deck angkin */
      SUBROUTINE ANGKIN(xmos,WORK,llwork)
C
C     J.I.Melo, July. 2005
C Linux version, Melo, (esta la hizo MELO solo !!)
C
C----------------------------------------
C
Cc      dalton.cm = lucme
Cc      dalton.out = lupri
Cc      aca esta nasht = 0.5 naos
Cc      aca esta norbt = nmos
Cc      aca estan los print levels =  iprsir
#include "priunit.h"
#include "infpri.h"
#include "inftap.h"
#include "inforb.h"
#include "chrnos.h"
#include "lrescinf.h"

C
C ------------------------------------------------------
C  PARAMETOS :		LEO SIRIUS.RST
C ------------------------------------------------------
c      INTEGER  NSYM, norb, nrhf
      INTEGER  number, k
      CHARACTER*8  B(4), Labl(10), LName(6)
      CHARACTER*8 RTNLBL(2),RTNLBLa(2),RTNLBLb(2)
      CHARACTER*12 char
C ------------------------------------------------------
C  PARAMETROS : 	MATRICES DE TRABAJO 
C ------------------------------------------------------
       integer  naos, nmos, nelec, lineas
       double precision Xmos, Cmos, Vec1, Vec2, Vec3,
     &    Mat1, Mat2, Mat3, MatKin, VecKin
!      the following become "AUTOMATIC ARRAYS"
!      better to allocate from WORK /hjaaj Sep 2013
       dimension xmos(NNBASX), CCMOS(NBAST,NBAST), 
     &    Vec1(NNBASX), Vec2(NNBASX), Vec3(NNBASX),
     &    Mat1(NBAST,NBAST), Mat2(NBAST,NBAST), Mat3(NBAST,NBAST),
     &    MatKin(NBAST,NBAST), VecKin(NNBASX)
C ------------------------------------------------------
C  TITLUOS :
C ------------------------------------------------------
      write(LUPRI,'(/721A1/)')('*',I=1,72)
      WRITE (LUPRI,*)' This is output from A.1.B :'
      WRITE(LUPRI,*) ' PSOK will be calculated as PSO.1.KIN '
      WRITE(LUPRI,*)'  and written into AOPROPER'
      write(lupri,*) 
      write(LUPRI,'(/72A1/)')('*',I=1,72)
C
C ======================================================
C ======================================================
C      -----------------------------------------------------
C      busco el numero de electrones en capa cerrada (nelec)
C      -----------------------------------------------------
        nelec = 0
        nelec = 2*NISHT
C      --------------------------------------
C      busco los orbitales molecualres (nmos)
C      --------------------------------------
        nmos = 0
        nmos = norbt
        
C      ----------------------------------------------------------
C      busco los orbitales atomicos < numero de orb. molec (naos)
C      ----------------------------------------------------------      

cz       if (JIMPRT .gt. 2) then 

        write(lupri,*) 'NBAST = NAOS ?? ', NBAST
        write(lupri,*) 'Miro variables comons.'
        write(lupri,*) 'LRATOM :', LRATOM
        write(lupri,*) 'JIMPRT :', JIMPRT
        write(lupri,*) 'iprsir = ', iprsir
        WRITE(lupri,*) 'nelec = 2*nisht?, nisht :', nisht
        WRITE(lupri,*) 'nelec = 2*NOCCT?, nocct :', nocct
        WRITE(lupri,*) 'nmos = norbt?, norbt :', norbt
        WRITE(lupri,*) 'naos = nbasis(iSYM) ?, nbas(1)', nbas(1)
cz       endif
C
C ------------------------------------------------------
C	LEO SIRIUS.RST PARA LOS CMOS
C ------------------------------------------------------
        NAOS = NBAST
        NMOS = NAOS 
        nelec = NOCCT
        number = NBAST*(NBAST+1)/2
        write(LUPRI,*)'nelec : nocct, nelec :', nelec 
        nelemsmol = NAOS*NMOS
         k = 1
         do i=1,naos
            do j=1,naos
            CCMOS(i,j) = xmos(k)
c            write(lupri,*) 'i,j', i,j, xmos(k)
            k = k + 1
            end do
         end do
         if (JIMPRT .gt. 3) then
            write(lupri,*)
            write(lupri,*) ' *************** CCMOS *****************'
            write(lupri,*) '                  OJO !! '
            write(lupri,*) '      en DALTON.OUT esta traspuesta '
            write(lupri,*) '              CCMOS(nelec,naos)   '
            write(lupri,*) ' ***************************************'
            do i=1,nelec
               write(lupri,*)
               write(lupri,*) 'Orbital Molecular Occ:',i
               do j=1,nelec
                  write(lupri,*) j,CCMOS(i,j)
               end do
            enddo
         endif
C ---------------------------------------------------------
C	SECCION CUENTA 
C ---------------------------------------------------------
C    Tengo :
C    naos        :  numero de orbitales de base atomica
C    nelec       :  numero de closed shell electrons
C    CCMOS(i,j)   :  Matrix coef. moleculares de SIRIUS.RST
C ---------------------------------------------------------
      write(lupri,'(/72A1/)')('*',I=1,72)
      write(lupri,*) 'empieza la cuenta'
      write(lupri,*) 'LRATOM is ', LRATOM
C
C Integral names to be calculated
C --------------------------------

C     ANGKIN for {Lkin,PSO}, if desired but not done here today
C
      LName(1)='XANGKINE'
      LName(2)='YANGKINE'
      LName(3)='ZANGKINE'
C
C     PSOK {L,PSOK}
C
      IJ = 4
      npos = 3*LRATOM-2
      DO I=npos,npos+2
         LName(IJ)= 'PSOKI'//CHRNOS(I/100)//CHRNOS(I/10)//
     &               CHRNOS(MOD(I,10))
         IJ = IJ + 1
      END DO
C
C     Print label names 
      write(LUPRI,*)' NAME of labels to be written on LUPROP :'
      DO I=1,6
         write(LUPRI,*)' #', I, LName(I)
      enddo
C 
C  Integrals needed for the calc.
C
      Labl(10)='KINENERG'
      Labl(1)='XANGMOM '
      Labl(2)='YANGMOM '
      Labl(3)='ZANGMOM '
      IJ=4
      DO I=npos,npos+2
        Labl(IJ)= 'PSO '//CHRNOS(I/100)//CHRNOS(I/10)//
     &         CHRNOS(MOD(I,10))
        IJ = IJ + 1
      END DO
C
C     Print label names 
      write(LUPRI,*)' Labels to be USED in calc :'
      DO I=1,10
         write(LUPRI,*)' #', I, Labl(I)
      enddo
C
C  Init of some matrices
C  -----------------------
      do mu=1,number
        vec1(mu)  = 0.0
        vec2(mu)  = 0.0
        vec3(mu)  = 0.0
      enddo
      Do j=1, naos
      Do i=1, naos
        Mat1(i,j) = 0.0
        Mat2(i,j) = 0.0
        Mat3(i,j) = 0.0
        MatKin(i,j)  = 0.0
      enddo
      enddo
C
      write(lupri,*)
c      CALL getmatrix(Labl(10),RTNLBLa,VecKin,number)
      CALL getmatrix('KINENERG',RTNLBLa,VecKin,number)
      CALL GPCLOSE(LUPROP,'KEEP')
      CALL Matform(VecKin,MatKin,naos,RTNLBLa(2),number,lupri)
C
C   Doing AngKin or L-Kin 
C
cxz      DO nlab=1, 3
cxz        write(lupri,*)
cxz        CALL Getmatrix(Labl(nlab),RTNLBLb,Vec1,number)
cxz        CALL Matform(Vec1,Mat1,naos,RTNLBLb(2),number,lupri)
cxzC
cxz        CALL CuentaAB(MatKin,Mat1,CCMOS,Mat2,nelec,naos)
cxzcz        if (JIMPRT.gt.2) then
cxz             write(lupri,*) '  Despues de cuenta AB = A.1.B : '
cxzc j2009            Call MatOut(Mat2,'        ', naos, lupri)
cxz             write(lupri,*)
cxzcz        endif
cxz        CALL Vecform(Vec2, Mat2,naos,RTNLBLb(2),number,lupri)
cxz        if (JIMPRT.gt.3) then
cxz           CALL HEADER(Labl(nlab),-1)
cxz           CALL OUTPAK(Vec2,Nbast,1,LUPRI)
cxz        endif
cxzC
cxz        Call Add(Vec2,Vec2,Vec3,number,1.0)
cxz        CALL GETDAT(RTNLBL(1),RTNLBL(2))
cxz        if (RTNLBLa(2).eq.RTNLBLb(2)) then 
cxz            RTNLBL(2) = RTNLBLa(2)
cxz          else 
cxz            RTNLBL(2)='ANTISYMM'
cxz        endif
cxz        write(lupri,*)'  Writing A.B Op.: ', LName(nlab),' ',RTNLBL
cxz        if (JIMPRT.gt.3) then
cxz           CALL HEADER(LName(nlab),-1)
cxz           CALL OUTPAK(Vec2,Nbast,1,LUPRI)
cxz           write(lupri,'(/72A1/)')('=',I=1,72)   
cxz        endif
cxz        write(*,*) ' about to write :' , LName(nlab)
cxz        CALL wrtpro(Vec2,number,LName(nlab),RTNLBL)
cxz        CALL GPCLOSE(LUPROP,'KEEP')
cxz      ENDDO
      write(lupri,'(/72A1/)')('*',I=1,72)
C
C   Doing Pso-Kin or Pso-K
C
      Do  nlab = 4, 6
        write(lupri,*)
C
        CALL Getmatrix(Labl(nlab),RTNLBLb,Vec1,number)
        CALL Matform(Vec1,Mat1,naos,RTNLBLb(2),number,lupri)
C  Do : A.B
        CALL CuentaAB(MatKin,Mat1,CCMOS,Mat2,nelec,naos,NBAST)
        CALL Vecform(Vec1,Mat2,naos,RTNLBLb(2),number,lupri)
C
        if (JIMPRT.gt.2) then
           write(lupri,*) '  Despues de cuenta AB = A.1.B : '
           char=Labl(nlab)//'.p2'
           CALL HEADER(char,-1)
           write(LUPRI,*) 
           if (JIMPRT.gt.3) then
              CALL OUTPAK(Vec1,Nbast,1,LUPRI)
           endif
        endif
C
C
C  Do : B.A
        CALL CuentaAB(Mat1,MatKin,CCMOS,Mat2,nelec,naos,NBAST)
        CALL Vecform(Vec2,Mat2,naos,RTNLBLb(2),number,lupri)
C
        if (JIMPRT.gt.2) then
           write(lupri,*) '  Despues de cuenta BA = B.1.A : '
           char='p2.'//Labl(nlab)
           CALL HEADER(char,-1)
           write(lupri,*)
           if (JIMPRT.gt.3) then
              CALL OUTPAK(Vec2,Nbast,1,LUPRI)
           endif
        endif
C
C  Add : A.B + B.A 
        Do i=1,Number
           Vec3(i) = Vec1(i) + Vec2(i) 
        enddo
c       Call Add(Vec1,Vec2,Vec3,number,1.0)
        if (JIMPRT.gt.2) then
           char='Total PSOKI : '//LName(nlab)
           CALL HEADER(char,-1)
           if (JIMPRT.gt.3) then
              CALL OUTPAK(Vec3,Nbast,1,LUPRI)
           write(lupri,'(/72A1/)')('=',I=1,72)
           endif
        endif
C
C  Write Section to AOPROPER
C
        CALL GETDAT(RTNLBL(1),RTNLBL(2))
        if (RTNLBLa(2).eq.RTNLBLb(2)) then
            RTNLBL(2) = RTNLBLa(2)
          else
             RTNLBL(2)='ANTISYMM'
c            RTNLBL(2)='SYMMETRI'
        endif
        write(lupri,*)'  Writing A.B Op.: ', LName(nlab),' ',RTNLBL
        if (JIMPRT.gt.3) then
           CALL HEADER(LName(nlab),-1)
           CALL OUTPAK(Vec3,Nbast,1,LUPRI)
           write(lupri,'(/72A1/)')('=',I=1,72)
        endif
        CALL wrtpro(Vec3,number,LName(nlab),RTNLBL)
        CALL GPCLOSE(LUPROP,'KEEP')
      ENDDO
      write(lupri,'(/72A1/)')('*',I=1,72)

  
      END
C =================================================================
C =================================================================
C
        SUBROUTINE Matform(VEC,MAT,naos,label,number,lupri)
C
C ***************************************************
C       PARAMETROS : PASA A MATRIZ UN VECTOR DADO
C   VEC   : Matriz in vector form
C   Mat   : Matriz in Matrix form
C   naos  : Numero de orbitales en base atomicas
C
C   nota : las matrices vienen en formato LT
C
C -------------------------------------------------
       double precision VEC, Mat
       character*8 label
       INTEGER   naos, number, lupri
       DIMENSION VEC(*), MAT(NAOS,NAOS)
C
C       write(*,*) '    Empieza VEC --> MAT'
C
      k = 1
      Do i=1,naos
         Do j=1,i
         mat(i,j) = vec(k)
         mat(j,i) = vec(k)
         if (label.eq.'ANTISYMM') then 
            mat(j,i) = -1.0*vec(k)
         endif
         k=k+1
         enddo
      enddo
cd      write(lupri,*) '     At MatForm : ', label
cd      write(lupri,*) 
cd      Call MatOut(Mat,label, naos, lupri)
cd      write(lupri,'(/72A1/)')('-',I=1,72)
C
      RETURN
      END
C =================================================================
C =================================================================
C
        SUBROUTINE Vecform(Vec,Mat,naos,label,number,lupri)
C
C ***************************************************
C       PARAMETROS
C   VEC   : Matriz in vector form
C   Mat   : Matriz in Matrix form
C   naos  : Numero de orbitales en base atomicas
C   label : 0:SYMM, 1: ANTYSIMM
C
C   nota : las matrices vienen en formato....
C
C -------------------------------------------------
       double precision VEC, Mat
       character*8 label
       INTEGER  naos, number, lupri, i, j
       DIMENSION VEC(*), MAT(naos,naos)
C
C       write(*,*) '    Empieza VEC --> MAT'
C
      k = 1
      Do i = 1, naos
         Do j = 1, i
         vec(k) = mat(i,j)
         k = k + 1
         enddo
      enddo
cd      write(lupri,*) '  At VecForm Symetri debagear aca: ', label
cd         do k=1,number
cd         write(lupri,*) '    ', k, vec(k)
cd         enddo
cd      write(lupri,*)
cd      write(lupri,'(/72A1/)')('-',I=1,72)
C
      RETURN
      END
C
C =================================================================
C =================================================================
C
        SUBROUTINE CuentaAB(A,B,Cmos,Cout,nelec,naos,NBAST)
C
C ***************************************************
C       PARAMETROS : A.B  --> A.1.B = Cout
C   A     : Matriz A vector form
C   B     : Matriz B vector form
C   Ccmos  : Matriz de coeficientes moleculares (Cmos), vec
C   nelec : Numero de electrones   !!! NOT USED !!!
C   naos : Numero de orbitales en base atomicas
C   -------------------------------------------------
       double precision A, B, Cmos , aux1, aux2, Cout
       double precision g11, g22, Sum, ge1, ge2
       INTEGER   nelec, naos, mu, n , i, j
       DIMENSION A(NBAST,NBAST),B(NBAST,NBAST), Cmos(NBAST,NBAST)
       DIMENSION aux1(NBAST,NBAST), aux2(NBAST,NBAST), Cout(NBAST,NBAST)
C ***************************************************
C       write(*,*) '          Empieza la cuenta AB = A.1.B'
C                                               
       Do i = 1,naos
          Do j = 1, naos
             aux1(i,j)=0.0
             aux2(i,j) =0.0
             Cout(i,j) =0.0
          enddo
       enddo
       Do mu = 1, naos
          Do n = 1, naos
          aux2(mu,n) = g22(mu,n,naos,A,Cmos,NBAST)
          aux1(n,mu) = g11(n,mu,naos,Cmos,B,NBAST)
          enddo
       enddo
C    
       Do mu = 1, naos
          Do nu = 1, naos
             Sum = 0.0
             Do k = 1, naos
                Sum = Sum + aux2(mu,k)*aux1(k,nu)
             enddo
             Cout(mu,nu) = Sum
          enddo
       enddo
C
cd       write(lupri,*) '  At Cuenta AB = A.1.B '
cd       Call MatOut(A,'A :     ', naos, lupri)
cd       Call MatOut(B,'B :     ', naos, lupri)
cd       Call MatOut(Cmos,'Cmos :  ', naos, lupri)
cd       Call MatOut(aux1,'g1(n,mu)', naos, lupri)
cd       Call MatOut(aux2,'g2(mu,n)', naos, lupri)
cd       Call MatOut(Cout,'Cout  : ', naos, lupri)
cd       write(lupri,'(/72A1/)')('-',I=1,72)
cd       write(lupri,*)
       return
       end
                                               
C ***********************************************************
       FUNCTION g22(i, mu, lim, mat1, mat2,nbas)
       integer*4 i, mu, lim, k
       double precision mat1, mat2, S, g22
       dimension mat1(nbas,nbas), mat2(nbas,nbas)
       S = 0.0
       do k=1,lim
          S = S + mat1(i,k)*mat2(mu,k)
       enddo
       g22 = S
       return
       end
C *********************************************
       FUNCTION g11(i,mu,lim,mat1,mat2,nbas)
       double precision  mat1, mat2, S, g11
       integer*4 i, mu, lim, k
       dimension mat1(nbas,nbas), mat2(nbas,nbas)
C OJO q estoy usando la mat 2 TRASPUESTA DE como viene !!!!
       S = 0.0
       do k=1,lim
          S = S + mat1(i,k)*mat2(k,mu)
       enddo
       g11 = S
       return
       end
C
C
C =================================================================
C =================================================================
       SUBROUTINE getmatrix(label,RTNLBL,TOTOVL,nelmnt)
#include "inforb.h"
#include "priunit.h"
#include "inftap.h"
#include "lrescinf.h"

       LOGICAL ANTSYM, FNDLB2
       CHARACTER*8 RTNLBL(2),label
       INTEGER  NELEMNT, NROW
       DIMENSION TOTOVL(NELMNT)
C
C
C     Read AO LABEL*8  matrix from AOPROPER file, 
C     devuelve VECTOR y RTNLABEL *2 date-symm
C
       NROW = Nbast
       NELMNT = nbast*(nbast + 1)/2
       CALL GPOPEN(LUPROP,'AOPROPER','OLD',' ','UNFORMATTED',IDUMMY,
     &  .FALSE.)
       REWIND(LUPROP)
       IF (FNDLB2(label,RTNLBL,LUPROP)) THEN
          CALL READT(LUPROP,NELMNT,TOTOVL)
          write(lupri,*)
          write(lupri,*) '  At getmatrix label: ', label,' ', RTNLBL
          if (JIMPRT .gt. 5) then 
             CALL HEADER('Read from AOPROPER',-1)
             CALL OUTPAK(TOTOVL,NROW,1,LUPRI)
             write(lupri,'(/72A1/)')('-',I=1,72)
          endif
       ELSE
       WRITE (*,'(//3A)') ' DALTON error: "',label,'" matrix'//
     *      ' not found on LUPROP.'
       CALL QUIT('@GETMATRIX error: property not found on LUPROP')
       END IF
cj       CALL GPCLOSE(LUPROP,'KEEP')
       RETURN
       END
C
C =============================================================
C
      Subroutine Add(A,B,C,number,sclfac)
                                             
            INTEGER Number
            double precision A, B, C, aux, sclfac
            dimension A(*),B(*),C(*)
C
C   deberia sumar dos vectores 
C
       Do i=1,Number
          C(i) = sclfac* ( A(i) + B(i) )
       enddo
       return
       end
C
C =============================================================
C
cc      Subroutine MatOut(mat,label,naos,lupri)
cc
cc            INTEGER naos, lupri
cc            character*8 label
cc            double precision mat 
cc            dimension mat(700,700)
ccC
cc        write(lupri,*) label
cc             do i=1,naos
cc             write(lupri,*) '   ',i, (Mat(i,j),j=1,naos)
cc             enddo
cc             write(lupri,*)
cc       return
cc       end
