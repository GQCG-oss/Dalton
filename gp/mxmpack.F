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
C
C$Id: mxmpack.F,v 1.3 2001-05-31 06:38:26 vebjornb Exp $
C
#include <single.h>
C  /* Deck mxm */
      SUBROUTINE MXMDUM
      END
#if defined (SYS_T90)
C...    Spastic T90 (at least the IEEE ones) has no MXM of its own...
C...    Call to SGEMM.     PRT
C
      SUBROUTINE MXM(A,NAR,B,NAC,C,NBC)
      DIMENSION A(NAR,NAC), B(NAC,NBC), C(NAR,NBC)
      CALL SGEMM('N', 'N', NAR, NBC, NAC, 1.0,
     &            A, NAR, B, NAC, 0.0, C, NAR)
      RETURN
      END
#endif
#if !defined (SYS_CRAY) && !defined (SYS_T3D)  && !defined (SYS_T90)
      SUBROUTINE MXM(A,NAR,B,NAC,C,NBC)
#include <implicit.h>
      DIMENSION A(NAR,NAC), B(NAC,NBC), C(NAR,NBC)
#if defined (VAR_F90)
C
C     tuh Mar 1988 - Alliant spends a lot of time setting up this call,
C     it might therefore be best to avoid calling MATMUL
C
      C=MATMUL(A,B)
#else
#if defined (VAR_ESSL)
C
C ESSL DGEMUL uses same cpu time as DGEMM /930330-hjaaj
C     but has no problem with NaN /931022-hjaaj
      CALL DGEMUL(A,NAR,'N',B,NAC,'N',C,NAR,NAR,NAC,NBC)
#endif
#endif
#if !defined (VAR_ESSL) && !defined (VAR_F90)
#if defined (SYS_CONVEX) || defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_AIX) || defined (VAR_BLAS3) || defined (SYS_SUN) || defined (SYS_LINUX) || defined (SYS_HPUX) || defined (SYS_HAL)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
C
C     Necessary to zero out because of NaN, tuh
C
      CALL DZERO(C,NAR*NBC)
      CALL DGEMM('N','N',NAR,NBC,NAC,D1,A,NAR,B,NAC,D0,C,NAR)
#else
C
C     C(I,J) = A(I,K) * B(K,J)
C
      DO 100 J = 1,NBC
         DO 100 I = 1,NAR
  100       C(I,J) = A(I,1)*B(1,J)
C
      DO 300 K = 2,NAC
         DO 300 J = 1,NBC
            BKJ = B(K,J)
            DO 300 I = 1,NAR
  300          C(I,J) = C(I,J) + A(I,K) * BKJ
C
#endif
#endif
      RETURN
      END
#endif
C  /* Deck mxmd */
#if defined (VAR_NOTUSED  )
C930522-hjaaj: MXMD not used in Hermit/Abacus any more
#if defined (VAR_F90)
      SUBROUTINE MXMD(A,NAR,NARA,B,NAC,NACA,C,NBC)
C
C     tuh Mar 1988 - Alliant spends a lot of time setting up this call,
C     it might therefore be best to avoid calling MATMUL
C
#include <implicit.h>
      DIMENSION A(NARA,NACA), B(NAC,NBC), C(NAR,NBC)
      C(1:NARA,:)=MATMUL(A,B(1:NACA,:))
      RETURN
      END
#endif
#if defined (SYS_CONVEX) || defined (VAR_ESSL) || defined (VAR_BLAS3) || defined (SYS_HAL)
      SUBROUTINE MXMD(A,NAR,NARA,B,NAC,NACA,C,NBC)
C     890905 tuh
#include <implicit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      DIMENSION A(NARA,NACA), B(NAC,NBC), C(NAR,NBC)
      CALL DZERO(C,NAR*NBC)
      CALL DGEMM('N','N',NARA,NBC,NACA,D1,A,NARA,B,NAC,D0,C,NAR)
      RETURN
      END
#endif
#if !defined (SYS_CONVEX) && !defined (VAR_ESSL) && !defined (VAR_F90) && !defined (VAR_BLAS3) && !defined (SYS_HAL)
      SUBROUTINE MXMD(A,NAR,NARA,B,NAC,NACA,C,NBC)
C
C     Modified MXM by Trygve U. Helgaker
C
#include <implicit.h>
      DIMENSION A(NARA,NACA), B(NAC,NBC), C(NAR,NBC)
C
C     C(1:NARA,:)=MATMUL(A,B(1:NACA,:))
C     C(I,J) = A(I,K) * B(K,J)
C
      DO 100 J = 1,NBC
         DO 100 I = 1,NARA
  100       C(I,J) = A(I,1)*B(1,J)
C
      DO 300 K = 2,NACA
         DO 300 J = 1,NBC
            BKJ = B(K,J)
            DO 300 I = 1,NARA
  300          C(I,J) = C(I,J) + A(I,K) * BKJ
C
      RETURN
      END
#endif
#endif
