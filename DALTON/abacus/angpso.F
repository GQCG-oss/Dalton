!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!

      SUBROUTINE ANGPSO(xmos,llwork,PRT,char1,char2,Av)
C
C Linux version, Melo, azua, giribet, alejandrito, 
C  aucar, romero, homero,bart y lisa
C
#include "inftap.h"
#include "inforb.h"
#include "priunit.h"
c#include "symmet.h"

C ------------------------------------------------------
C        DOUBLE PRECISION matrix
C        dimension matrix(*,*)
C ------------------------------------------------------
C  PARAMETROS :       LEO DALTON.INP
C ------------------------------------------------------
      DOUBLE PRECISION TIME, SECOND
      CHARACTER*8 B0(2), LABEL(2),char1, char2
      integer PRT,nlab
C ------------------------------------------------------
C  PARAMETROS :       LEO DALTON.CM
C ------------------------------------------------------

C included inforb to have bas info

C ------------------------------------------------------
C  PARAMETOS :        LEO SIRIUS.RST
C ------------------------------------------------------
      CHARACTER*8  STARS, B(4)
      PARAMETER (STARS = '********')
      double precision cmos(NBAST,NBAST), xmos(5000000)
      integer*8 k

C ------------------------------------------------------
C  PARAMETROS :       WORKING ARRAYS  
C ------------------------------------------------------
C
      double precision A1(NBAST,NBAST), A2(NBAST,NBAST), ORBC(NBAST)
      double precision OVERLP(NBAST,NBAST), Av, uno , p(NBAST,NBAST), S1
C
C
C ======================================================
C ======================================================
C ------------------------------------------------------
C       READ DALTON.INP  and take out INTEGRALS
C ------------------------------------------------------
        B0(1)='12345678'
        B0(2)='12345678'
c       nlab=2 because we do in pairs from outside de loops obver
c       all teh integrals, to avoid extra calculations
c       LABEL(1)='XANGMOM'
c       LABEL(2)='YANGMOM'
c       LABEL(3)='ZANGMOM'
c       LABEL(4)='PSO 001'
c       LABEL(5)='PSO 002'
c       LABEL(6)='PSO 003'
C ------------------------------------------------------
C       READ DALTON.CM  for nbas(naos). Deprecated. not read anymore
C ------------------------------------------------------
       naos = NBAST 
       nmos= naos
c       nelec=10
c       write(*,*) 'nelec = NOCCT ?. elec. cpa closed: ', NOCCT, nelec
c       write(*,*) ' NBAST', NBAST
 
c      write(lupri,*)
c      write(lupri,'(/72A1/)')('*',I=1,72)
C
C ------------------------------------------------------
C     READ SIRIUS.RST for CMOS
C ------------------------------------------------------

cz    write(*,*) 'Opening SIRIUS.RST file'
cz    write(*,*)

            nelemsmol = naos*nmos
c            write(lupri,*) 'nelemsmol =',naos,'*',nmos
c            write(*,*)
            k = 1
c            write(lupri,*) '* Write  xmos vector'
c            write(lupri,*) '** write i,j,k,cmos(i,j)=xmos(k)'
c            write(lupri,*) 'Chq  DALTON.OUT warning it has transopsed '
c            write(lupri,*) xmos
c            write(lupri,*)
c            write(lupri,*)
C *jmelo*            do i=1,nelec/2(=NOCCT)  in case of occ orbs
            do i=1,naos
               do j=1,naos
               cmos(i,j) = xmos(k)
C               write(lupri,*) i,j,k, cmos(i,j)
               k = k + 1
              end do
            end do
C
C
C ----------------------------------------------------------
C       We do the matrix prods here  
C ----------------------------------------------------------
C    Tengo : most of this are deprecated, we now use an include
C    naos        :  numero de orbitales de base atomica
C    nmols       :  numero de orbitales de base molec 
C    nelec(NOCCT):  Orbitales Ocupados, numero de closed shell electrons
C    cmos(i,j)   :  Matrix coef. moleculares de SIRIUS.RST
C    LABEL(nlab) :  Leo las integrales a usar en DALTON.INP
C ----------------------------------------------------------
      nlab=2
      LABEL(1)=char1
      LABEL(2)=char2
c      write(lupri,*) 
c      write(lupri,*) '*******************************************' 
c      write(lupri,*) '***** Debug : Starting Calculations *******' 
c      write(lupri,*) '*******************************************' 
c      write(lupri,*) 
      CALL READAOP( OVERLP, naos, 'OVERLAP ')
c      write(lupri,*) ' ****** OVERLAP *******'
C      write(lupri,*) ((OVELP(i,j),j=1,naos), i=1,naos)
      IF (B0(1).Eq.'ORBITAL') then
      write(lupri,*)' ****** ORBITAL CONTRIBUTUIONS DETECTED !  *******'
      ENDIF 

cjim-dbg       DO j = 1, nlab
cjim-dbg       DO k = j, nlab
C  we keep the general form in case we use it somewhere else.
        j=1
        k=2
        do i1=1,NBAST
          ORBC(i1) = 0.0
          do i2=1,NBAST  
            A2(i1,i2) = 0.0
          enddo
        enddo
        CALL READAOP( A1, naos, LABEL(j))
        CALL READAOP( A2, naos, LABEL(k))
c       write(lupri,*) 'Llamando Cuenta para : ', LABEL(j),'  ', LABEL(k)
        CALL cuenta(ORBC,A1, A2,Cmos,NOCCT,naos,nmos,NBAST)
        IF (B0(1).Eq.'ORBITAL') then
            write(lupri,*) '< ',LABEL(j),' >< ',LABEL(k),'>' 
        ENDIF 
        S1 = 0.0
c jix : on NaN issue . should be  nelec/2
        do i1 = 1,NOCCT
           IF (B0(1).Eq.'ORBITAL') then
C              write(lupri,*) 'Orbital : ',i1, '  Valor', ORBC(i1)
               write(lupri,*) ORBC(i1)
           ENDIF 
           S1 = S1 + ORBC(i1)
        enddo
        Av = 2.0*S1
        IF (PRT.GT.2) THEN 
        write(lupri,*) '@ANGPSO< ',LABEL(j),' >< ',LABEL(k),' >  = ',Av 
c           write(*,*) '< ',LABEL(j),' >< ',LABEL(k),' >  = ', Av 
           if (B0(1).Eq.'ORBITAL') then
              write(lupri,*) '======================================= '
           endif
        ENDIF
cjim-dbg       ENDDO
C   *****************************
C      Writing the Matrices 
C   *****************************
      IF (B0(2).Eq.'INTEGRAL') then
         write(lupri,*) ' ****** INTEGRALS OUTPUT   DETECTED !  *******'
          write(lupri,*) '        in atomic basis set  '
          write(lupri,*) LABEL(j)
          do mu=1,naos
          do nu=mu,naos
             write(lupri,*) mu, nu, A1(mu,nu)
          enddo
          enddo
      ENDIF
C
cjim-dbg       ENDDO
c      write(lupri,*) 
c      write(lupri,*) 
C
C ----------------------------------------------------------
C         SECCION CHEQUEO DE LA BASE
C ----------------------------------------------------------
      IF (PRT.GT.3) THEN 
         write(lupri,*) '*******************************************' 
         write(lupri,*) '********** Basis set check    *************'
         write(lupri,*) '*******************************************' 
         uno = 0.0
         do i = 1,NBAST
            do j = 1,NBAST
            p(i,j) = 0.0
            A2(i,j)= 0.0
            enddo
         enddo
C        write(lupri,*) '********* OVERLAP **********' 
C        write(lupri,*) ((OVELP(i,j),j=1,naos), i=1,naos)
C        write(lupri,*) 
         IF (B0(2).Eq.'INTEGRAL') then
            write(lupri,*) ' **** INTEGRALS OUTPUT   DETECTED ! *******'
            write(lupri,*) ' ************** Cmos OJO !! **************' 
            write(lupri,*) ' ***** en DALTON.OUT has the transposed **' 
            write(lupri,*) ' ******      Cmos(NOCCT,naos)       ******' 
            write(lupri,*) ' *****************************************' 
C           write(lupri,*) ((Cmos(i,j),j=1,naos), i=1,NOCCT)
            do i=1,NOCCT
               write(lupri,*) (Cmos(i,j),j=1,naos)
            enddo
         ENDIF
C        ---------------------------------------------
C          CHECK  Cmos has norm = 1 in atomic basis !!
C        ---------------------------------------------
         write(lupri,*) 
         write(lupri,*) 'Norm of  Cmos by rows in atomic basis :' 
         Avv = 0.0
         do i = 1,NOCCT
            do j = 1,naos
               uno = uno + Cmos(i,j)*Cmos(i,j)
               uno = sqrt(uno)
            enddo
            write(lupri,*) '  Orbital', i, ' : ',uno
            Avv = Avv + uno 
         enddo
         write(lupri,*) '** Total summ of row norms :', Avv
         write(lupri,*) 
C     ------------------------
C     CHEQUEO Tr(PS) Mulliken 
C     ------------------------
         Avv = 0.0
         do i = 1,naos
            do j = 1,naos
               do k=1,NOCCT
               p(i,j) = p(i,j) + Cmos(k,i)*Cmos(k,j)
               enddo
            enddo
         enddo
         do i = 1,naos
            do j = 1,naos
               do k=1,naos
               A2(i,j) = A2(i,j) + P(i,k)*OVERLP(k,j)
               enddo
            enddo
         enddo
         do k=1,naos
            Avv = Avv + A2(k,k)
         enddo
         write(lupri,*) '  Tr(P.S) = ',Avv
      ENDIF
      RETURN
      END
C ****************************************************************

        SUBROUTINE READAOP(A, n, LABEL)
C
C Linux version, Melo, azua, giribet, alejandrito, aucar, romero, homero,
C  bart y lisa
C
C       A : Devuelve la matriz de los coeficientes buscados
C       n : Dim(A)
C   LABEL : Etiqueta de la matriz buscada
C
C -------------------------------------------------
C     PARAMETOS DE LEO AOPROPER
C -------------------------------------------------
      CHARACTER*8  STARS, LABEL, B(4)
      double precision A(n,n)
      integer*4 nelem,  mu, nu
      PARAMETER (STARS = '********')
      PARAMETER (LUINP222 = 222)
      real(8), allocatable :: xx(:)
C ===================================================
C ------------------------------------------------
C       LEO AOPROPER Y ARMO EL VECTOR xx 
C                  con las integrales 
C         
C ------------------------------------------------
C
       open(222,file='AOPROPER',status='OLD',FORM='UNFORMATTED')
C       write(*,*) 'Opening AOPROPER file to look for :', LABEL
C       write(*,*) '*******************************************'
C       write(*,*)
       REWIND LUINP222
       IREC = 0
       IERR = 0
  1    READ (LUINP222,END=3,ERR=2) B
       IREC = IREC + 1
       IF (B(1) .NE. STARS) GO TO 1
C          WRITE (*, '(5X,I5,3X,4(2X,A8))')  IREC, B
           IF (B(4).EQ.LABEL) then
              IF (B(3).EQ.'SYMMETRI') then
                 nelem = n*(n+1)/2
                 allocate(xx(nelem))
                 read (LUINP222,end=3,err=2) (xx(i),i=1,nelem)
C                 write(*,*) 'Found : ', B ,' on FILE AOPPROPER'
                 k = 1
                 do nu=1,n
                    do mu=1,nu
                     A(mu,nu) = xx(k)
                     A(nu,mu) = A(mu,nu)
C                      write(*,*) mu, nu, a(mu,nu)
                     k = k + 1 
                     enddo
                 enddo
                 deallocate(xx)
              ELSE
                 nelem = n*(n+1)/2
                 allocate(xx(nelem))
                 read (LUINP222,end=3,err=2) (xx(i),i=1,nelem)
C                 write(*,*) 'Found : ', B ,' on FILE AOPPROPER'
                 k = 1
                 do nu = 1,n
                    do mu = 1,nu
C                      IF (mu.EQ.nu) then
C                        A(mu,nu) = 0.0 
C                      ELSE
                        A(mu,nu) = xx(k)
                        A(nu,mu) = -1.0*A(mu,nu)
                        k = k + 1 
C                        write(*,*) mu, nu, a(mu,nu)
C                      ENDIF
                    enddo
                 enddo   
                 deallocate(xx)
              ENDIF  
           ENDIF   
      GO TO 1
C
   2  CONTINUE
      IREC = IREC + 1
      IERR = IERR + 1
      WRITE (*, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
      REWIND LUINP222
      DO 102 I = 1,IREC
         READ (LUINP222) J
  102 CONTINUE
      IF (IERR .LE. 2) GO TO 1
  202 CONTINUE
         READ (LUINP222,END=3) J
         IREC = IREC + 1
      GO TO 202
C
   3  CONTINUE
C     WRITE (*,'(/I10,A)') IREC,
C    *   ' records read before EOF on file. AOPROPER'
      close(LUINP222)
C      write(*,*) 'Closing File AOPROPER '
C      write(*,*) '*****************************************'
      RETURN
      END
C
C ***********************************************************************
C
      SUBROUTINE cuenta(ORBC,A1,A2,C,noc,naos,nmos,NBAST)
C
C     Melo, azua, giribet, alejandrito, aucar, romero, homero,
C  bart and lisa
C
C ***********************************************************************
C   PARAMETROS  
C   ORBC(a) : Orbital Contribution of <A.B>
C    A1, A2 : Matrices con integrales en base atomica
C         C : Matriz de coeficientes moleculares (cmos)
C     nelec : Numero de electrones now : noc
C      naos : Numero de orbitales en base atomicas 
C      nmos : Numero de obitales moleculares
C       a   : Ocupados = nelec/2, AHORA 'noc'
C       i   : Todos  = nmos
C -------------------------------------------------
       DIMENSION A1(NBAST,NBAST), A2(NBAST,NBAST), C(NBAST,NBAST)
       DIMENSION CT(NBAST,NBAST),gama(NBAST,NBAST)
       DIMENSION alfa(NBAST,NBAST), beta(NBAST,NBAST)
       DIMENSION delta(NBAST,NBAST)
       double precision A1, A2, C, CT, ORBC(5000000)
       double precision alfa, beta, gama, delta, g1, g2
C       double precision mat1, mat2
       INTEGER*4  naos, nmos, noc, a, i, mu, sigma
C ***************************************************
C ***************************************************
C ------------------------------------------------   
C      GENERO C TRANSPUESTA
          Do i = 1, nmos
             Do j = 1, naos
             CT(j,i)= C(i,j)
             enddo
          enddo
C ------------------------------------------------   

C       write(*,*) 'Empieza la cuenta'
       S1  = 0.0
       Av  = 0.0
       aux = 0.0
C ***********************************************************
C
C       write(*,*) 'alfa'     
       Do mu = 1, naos
          Do i = 1, nmos
          alfa(mu,i)= g2(mu,i,naos,A1,C,NBAST)
C         write(*,*) alfa(mu,i)
          enddo
       enddo
C       write(*,*) '****'
C       write(*,*) 'beta'     
       Do sigma=1, naos
          Do a =1, noc
          beta(sigma,a)= g2(sigma,a,naos,A2,C,NBAST)
C         write(*,*) beta(sigma,a)
          enddo
       enddo
C       write(*,*) 'gama'     
       Do a = 1, noc
          Do i = 1, nmos
          gama(a,i)= g1(a,i,naos,CT,alfa,NBAST)
C          write(*,*) gama(a,i)
          enddo
       enddo
C       write(*,*) 'delta'     
       Do i = 1, nmos
          Do a = 1, noc
          delta(i,a)=  g1(i,a,naos,CT,beta,NBAST)
C          write(*,*) delta(i,a)
          enddo
       enddo
C jix : on NaN issue : should be  nelec/2 
       do a = 1,noc
          do i = 1,nmos
          ORBC(a) = ORBC(a) + gama(a,i)*delta(i,a)
          enddo
       enddo
       return
       end
C ***********************************************************
       FUNCTION g2(i, mu, lim, mat1, mat2,NBAST)
       dimension mat1(NBAST,NBAST), mat2(NBAST,NBAST)
       double precision mat1, mat2, S, g2
       integer*4 i, mu, lim, k
       S = 0.0
C        write(*,*)   mat1
C       write(*,*) 'Function g2'
C       write(*,*) '****'
       do k=1,lim
          S = S + mat1(i,k)*mat2(mu,k)
       enddo
       g2 = S
C       write(*,*) 'FIN g2 : ', S
       return
       end
C *********************************************
       FUNCTION g1(i,mu,lim,mat1,mat2,NBAST)
       dimension mat1(NBAST,NBAST), mat2(NBAST,NBAST)
       double precision  mat1, mat2, S, g1
       integer*4 i, mu, lim, k
       S = 0.0
       do k=1,lim
          S = S+ mat1(k,i)*mat2(k,mu)
       enddo
       g1 = S
       return
       end

