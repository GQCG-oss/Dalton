C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
      subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      dimension a(1), b(1), c(1)
      call sgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      return
      end
      subroutine dgemv (trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      dimension a(1),x(1),y(1)
      call sgemv (trans,m,n,alpha,a,lda,x,incx,beta,y,incy)
      return
      end
      function ddot(n,x,nx,y,ny)
      dimension x(1),y(1)
      ddot=sdot(n,x,nx,y,ny)
      return
      end
      subroutine dswap(n,x,nx,y,ny)
      dimension x(1),y(1)
      call sswap(n,x,nx,y,ny)
      return
      end
      subroutine dcopy(n,x,nx,y,ny)
      dimension x(1),y(1)
      call scopy(n,x,nx,y,ny)
      return
      end
      subroutine dscal(n,a,x,nx)
      dimension x(1)
      call sscal(n,a,x,nx)
      return
      end
      subroutine daxpy(n,a,x,nx,y,ny)
      dimension x(1),y(1)
      call saxpy(n,a,x,nx,y,ny)
      return
      end
      integer function idamax(n,x,nx)
      dimension x(1)
      idamax=isamax(n,x,nx)
      return
      end
      function dnrm2(n,x,incx)
      dimension x(1)
      dnrm2=snrm2(n,x,incx)
      return
      end
      function dasum(n,x,incx)
      dimension x(1)
      dasum=sasum(n,x,incx)
      return
      end

