      Subroutine Jacobi(A,C,N,LdC)
************************************************************************
*                                                                      *
*     (c) Copyright. All rights reserved                               *
*                                                                      *
*     No part of this code may be copied, redestributed or included    *
*     into any commercial product without the written permission of    *
*     the author. The use is restriced to research purposes only.      *
*                                                                      *
************************************************************************
*                                                                      *
*     Jacobi diagonalization                                           *
*                                                                      *
*     calling parameters:                                              *
*     A       : array of double precision real, input/output           *
*               On input this is the matrix to be diagonalized.        *
*               On ouput the diagonal matrix is returned.              *
*     C       : array of double precision real                         *
*               On input this are the basis vectors.                   *
*               The rotated basis vectors are returned.                *
*     N       : Integer, input                                         *
*               Dimension of the matrix to be diagonalized             *
*     LdC     : Integer, input                                         *
*               Leading dimension of matrix C.                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher, University of Lund, Sweden, 1993                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
*
      Dimension A(*)
      Dimension C(LdC,*)
*
      Parameter ( Eps=1.D-12 )
      Parameter ( EpsSqr=Eps*Eps )
*----------------------------------------------------------------------*
*     initialisation step:                                             *
*     - compute the largest subdiagonal element                        *
*----------------------------------------------------------------------*
      ASqrMax=0.0D0
      Do i=2,N
         ii=i*(i-1)/2
         Do j=1,i-1
            Temp=A(ii+j)*A(ii+j)
            ASqrMax=Max(ASqrMax,Temp)
         End Do
      End Do
      If ( ASqrMax.le.EpsSqr ) Return  
      Thrs=0.1*ASqrMax
*----------------------------------------------------------------------*
*     start loop over sweeps:                                          *
*     - no pivot elementsearch is implemented. However, the rotation   *
*       is applied only to those submatrix elements which are bigger   *
*       than a threshold. The latter is based on the value of the      *
*       Largest subdiagonal element.                                   *
*----------------------------------------------------------------------*
      iter=0
100   iter=iter+1
      Do i=2,N
         ii=i*(i-1)/2
         Do j=1,i-1
            jj=j*(j-1)/2
            Aij=A(ii+j)
            AijSqr=Aij*Aij
            If ( AijSqr.ge.Thrs ) then
*---  compute the rotation angle --------------------------------------*
               Aii=A(ii+i)
               Ajj=A(jj+j)
               Diff=Aii-Ajj
               SigRot=1.0
               If ( Diff.lt.0.0D0 ) then
                  SigRot=-SigRot
                  Diff=-Diff
               End If
               Temp=Diff+SQRT(Diff*Diff+4.0*AijSqr)
               TanA=2.0*SigRot*Aij/Temp
               CosA=1.0/SQRT(1.0+TanA*TanA)
               SinA=CosA*TanA
*---  update rows/columnes of the matrix A ----------------------------*
               kmin=1  
               kmax=j-1
               kleft=Mod(kmax-kmin+1,4)
               If ( kleft.eq.1 ) then
                  Aii1=A(ii+1)
                  Ajj1=A(jj+1)
                  A(ii+1)=SinA*Ajj1+CosA*Aii1
                  A(jj+1)=CosA*Ajj1-SinA*Aii1
               Else If ( kleft.eq.2 ) then
                  Ajj1=A(jj+1)
                  Ajj2=A(jj+2)
                  Aii1=A(ii+1)
                  Aii2=A(ii+2)
                  A(ii+1)=SinA*Ajj1+CosA*Aii1
                  A(ii+2)=SinA*Ajj2+CosA*Aii2
                  A(jj+1)=CosA*Ajj1-SinA*Aii1
                  A(jj+2)=CosA*Ajj2-SinA*Aii2
               Else If ( kleft.eq.3 ) then
                  Ajj1=A(jj+1)
                  Ajj2=A(jj+2)
                  Ajj3=A(jj+3)
                  Aii1=A(ii+1)
                  Aii2=A(ii+2)
                  Aii3=A(ii+3)
                  A(ii+1)=SinA*Ajj1+CosA*Aii1
                  A(ii+2)=SinA*Ajj2+CosA*Aii2
                  A(ii+3)=SinA*Ajj3+CosA*Aii3
                  A(jj+1)=CosA*Ajj1-SinA*Aii1
                  A(jj+2)=CosA*Ajj2-SinA*Aii2
                  A(jj+3)=CosA*Ajj3-SinA*Aii3
               End If
               kmin=kmin+kleft
               Do k=kmin,kmax,4
                  Ajj0=A(jj+k+0)
                  Ajj1=A(jj+k+1)
                  Ajj2=A(jj+k+2)
                  Ajj3=A(jj+k+3)
                  Aii0=A(ii+k+0)
                  Aii1=A(ii+k+1)
                  Aii2=A(ii+k+2)
                  Aii3=A(ii+k+3)
                  A(ii+k+0)=SinA*Ajj0+CosA*Aii0
                  A(ii+k+1)=SinA*Ajj1+CosA*Aii1
                  A(ii+k+2)=SinA*Ajj2+CosA*Aii2
                  A(ii+k+3)=SinA*Ajj3+CosA*Aii3
                  A(jj+k+0)=CosA*Ajj0-SinA*Aii0
                  A(jj+k+1)=CosA*Ajj1-SinA*Aii1
                  A(jj+k+2)=CosA*Ajj2-SinA*Aii2
                  A(jj+k+3)=CosA*Ajj3-SinA*Aii3
               End Do
               kmin=j+1
               kmax=i-1
               kleft=Mod(kmax-kmin+1,4)
               kk=jj+j
               If ( kleft.eq.1 ) then
                  k0=kk
                  Ak0j=A(k0+j)
                  Aii0=A(ii+kmin+0)
                  A(ii+kmin+0)=SinA*Ak0j+CosA*Aii0
                  A(kk+j)=CosA*Ak0j-SinA*Aii0
                  kk=k0+kmin
               Else If ( kleft.eq.2 ) then
                  k0=kk
                  k1=k0+kmin
                  Ak0j=A(k0+j)
                  Ak1j=A(k1+j)
                  Aii0=A(ii+kmin+0)
                  Aii1=A(ii+kmin+1)
                  A(k0+j)=CosA*Ak0j-SinA*Aii0
                  A(k1+j)=CosA*Ak1j-SinA*Aii1
                  A(ii+kmin+0)=SinA*Ak0j+CosA*Aii0
                  A(ii+kmin+1)=SinA*Ak1j+CosA*Aii1
                  kk=k1+kmin+1
               Else If ( kleft.eq.3 ) then
                  k0=kk
                  k1=k0+kmin
                  k2=k1+kmin+1
                  Ak0j=A(k0+j)
                  Ak1j=A(k1+j)
                  Ak2j=A(k2+j)
                  Aii0=A(ii+kmin+0)
                  Aii1=A(ii+kmin+1)
                  Aii2=A(ii+kmin+2)
                  A(k0+j)=CosA*Ak0j-SinA*Aii0
                  A(k1+j)=CosA*Ak1j-SinA*Aii1
                  A(k2+j)=CosA*Ak2j-SinA*Aii2
                  A(ii+kmin+0)=SinA*Ak0j+CosA*Aii0
                  A(ii+kmin+1)=SinA*Ak1j+CosA*Aii1
                  A(ii+kmin+2)=SinA*Ak2j+CosA*Aii2
                  kk=k2+kmin+2
               End If
               kmin=kmin+kleft
               Do k=kmin,kmax,4
                  k0=kk
                  k1=k0+k
                  k2=k1+k+1
                  k3=k2+k+2
                  Ak0j=A(k0+j)
                  Ak1j=A(k1+j)
                  Ak2j=A(k2+j)
                  Ak3j=A(k3+j)
                  Aii0=A(ii+k+0)
                  Aii1=A(ii+k+1)
                  Aii2=A(ii+k+2)
                  Aii3=A(ii+k+3)
                  A(k0+j)=CosA*Ak0j-SinA*Aii0
                  A(k1+j)=CosA*Ak1j-SinA*Aii1
                  A(k2+j)=CosA*Ak2j-SinA*Aii2
                  A(k3+j)=CosA*Ak3j-SinA*Aii3
                  A(ii+k+0)=SinA*Ak0j+CosA*Aii0
                  A(ii+k+1)=SinA*Ak1j+CosA*Aii1
                  A(ii+k+2)=SinA*Ak2j+CosA*Aii2
                  A(ii+k+3)=SinA*Ak3j+CosA*Aii3
                  kk=k3+k+3
               End Do
               kmin=i+1
               kmax=N
               kleft=Mod(kmax-kmin+1,4)
               kk=ii+i
               If ( kleft.eq.1 ) then
                  k0=kk
                  Ak0j=A(k0+j)
                  Ak0i=A(k0+i)
                  A(k0+j)=CosA*Ak0j-SinA*Ak0i
                  A(k0+i)=SinA*Ak0j+CosA*Ak0i
                  kk=k0+kmin
               Else If ( kleft.eq.2 ) then
                  k0=kk
                  k1=k0+kmin
                  Ak0j=A(k0+j)
                  Ak0i=A(k0+i)
                  Ak1j=A(k1+j)
                  Ak1i=A(k1+i)
                  A(k0+j)=CosA*Ak0j-SinA*Ak0i
                  A(k0+i)=SinA*Ak0j+CosA*Ak0i
                  A(k1+j)=CosA*Ak1j-SinA*Ak1i
                  A(k1+i)=SinA*Ak1j+CosA*Ak1i
                  kk=k1+kmin+1
               Else If ( kleft.eq.3 ) then
                  k0=kk
                  k1=k0+kmin
                  k2=k1+kmin+1
                  Ak0j=A(k0+j)
                  Ak0i=A(k0+i)
                  Ak1j=A(k1+j)
                  Ak1i=A(k1+i)
                  Ak2j=A(k2+j)
                  Ak2i=A(k2+i)
                  A(k0+j)=CosA*Ak0j-SinA*Ak0i
                  A(k0+i)=SinA*Ak0j+CosA*Ak0i
                  A(k1+j)=CosA*Ak1j-SinA*Ak1i
                  A(k1+i)=SinA*Ak1j+CosA*Ak1i
                  A(k2+j)=CosA*Ak2j-SinA*Ak2i
                  A(k2+i)=SinA*Ak2j+CosA*Ak2i
                  kk=k2+kmin+2
               End If
               kmin=kmin+kleft
               Do k=kmin,kmax,4
                  k0=kk
                  k1=k0+k
                  k2=k1+k+1
                  k3=k2+k+2
                  Ak0j=A(k0+j)
                  Ak0i=A(k0+i)
                  Ak1j=A(k1+j)
                  Ak1i=A(k1+i)
                  Ak2j=A(k2+j)
                  Ak2i=A(k2+i)
                  Ak3j=A(k3+j)
                  Ak3i=A(k3+i)
                  A(k0+j)=CosA*Ak0j-SinA*Ak0i
                  A(k0+i)=SinA*Ak0j+CosA*Ak0i
                  A(k1+j)=CosA*Ak1j-SinA*Ak1i
                  A(k1+i)=SinA*Ak1j+CosA*Ak1i
                  A(k2+j)=CosA*Ak2j-SinA*Ak2i
                  A(k2+i)=SinA*Ak2j+CosA*Ak2i
                  A(k3+j)=CosA*Ak3j-SinA*Ak3i
                  A(k3+i)=SinA*Ak3j+CosA*Ak3i
                  kk=k3+k+3
               End Do
*---  update the diagonal elements of A -------------------------------*
               Temp=2.0*CosA*SinA*Aij
               CosA2=CosA*CosA
               SinA2=SinA*SinA
               A(jj+j)=SinA2*Aii+CosA2*Ajj-Temp
               A(ii+j)=0.0D0
               A(ii+i)=CosA2*Aii+SinA2*Ajj+Temp
*---  update rows/columnes of the eigenvectors C ----------------------*
               kmin=1  
               kmax=LdC
               kleft=Mod(kmax-kmin+1,4)
               If ( kleft.eq.1 ) then
                  C1j=C(1,j)
                  C1i=C(1,i)
                  C(1,i)=SinA*C1j+CosA*C1i
                  C(1,j)=CosA*C1j-SinA*C1i
               Else If ( kleft.eq.2 ) then
                  C1j=C(1,j)
                  C2j=C(2,j)
                  C1i=C(1,i)
                  C2i=C(2,i)
                  C(1,i)=SinA*C1j+CosA*C1i
                  C(2,i)=SinA*C2j+CosA*C2i
                  C(1,j)=CosA*C1j-SinA*C1i
                  C(2,j)=CosA*C2j-SinA*C2i
               Else If ( kleft.eq.3 ) then
                  C1j=C(1,j)
                  C2j=C(2,j)
                  C3j=C(3,j)
                  C1i=C(1,i)
                  C2i=C(2,i)
                  C3i=C(3,i)
                  C(1,i)=SinA*C1j+CosA*C1i
                  C(2,i)=SinA*C2j+CosA*C2i
                  C(3,i)=SinA*C3j+CosA*C3i
                  C(1,j)=CosA*C1j-SinA*C1i
                  C(2,j)=CosA*C2j-SinA*C2i
                  C(3,j)=CosA*C3j-SinA*C3i
               End If
               kmin=kmin+kleft
               Do k=kmin,kmax,4
                  C0j=C(k+0,j)
                  C1j=C(k+1,j)
                  C2j=C(k+2,j)
                  C3j=C(k+3,j)
                  C0i=C(k+0,i)
                  C1i=C(k+1,i)
                  C2i=C(k+2,i)
                  C3i=C(k+3,i)
                  C(k+0,i)=SinA*C0j+CosA*C0i
                  C(k+1,i)=SinA*C1j+CosA*C1i
                  C(k+2,i)=SinA*C2j+CosA*C2i
                  C(k+3,i)=SinA*C3j+CosA*C3i
                  C(k+0,j)=CosA*C0j-SinA*C0i
                  C(k+1,j)=CosA*C1j-SinA*C1i
                  C(k+2,j)=CosA*C2j-SinA*C2i
                  C(k+3,j)=CosA*C3j-SinA*C3i
               End Do
            End If
         End Do
      End Do
*---  find the value of the largest subdiagonal element ---------------*
      ASqrMax=0.0D0
      Do i=2,N
         ii=i*(i-1)/2
         Do j=1,i-1
            Temp=A(ii+j)*A(ii+j)
            ASqrMax=Max(ASqrMax,Temp)
         End Do
      End Do
*---  update the threshold --------------------------------------------*
      Thrs=0.1*ASqrMax
*---  check the accuracy ----------------------------------------------*
      If ( ASqrMax.gt.EpsSqr ) Goto 100
*----------------------------------------------------------------------*
*     Thats the end                                                    *
*----------------------------------------------------------------------*
200   Return
      End 
