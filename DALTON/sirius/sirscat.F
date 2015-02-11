C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Agren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
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
C
C
C  /* Deck sca_vw */
      SUBROUTINE SCA_VW(NIJ,IJLIST,VWMAT,CMO,WORK,KFREE,LFREE)
C
C Copyright 5-Nov-1997 Hans Joergen Aa. Jensen.
C
C This subroutines returns (ij = 1,NIJ):
#ifndef SCA_VW
C    VWMAT(a,b,1,ij) = (IJ | ab)
C    VWMAT(a,b,2,ij) = (Ia | Jb)
#else
C    VWMAT(a,b,1,ij) = <Ia || Jb> = (IJ|ab) - (Ib|Ja)
C    VWMAT(a,b,2,ij) = <IJ || ab> = (Ia | Jb) - (Ib | Ja)
#endif
C where
C    I = IJLIST(1,ij)
C    J = IJLIST(2,ij)
C
C
#include "implicit.h"
      DIMENSION IJLIST(2,NIJ)
      DIMENSION VWMAT(N2BASX,2,NIJ), CMO(NCMOT), WORK(*)
C
C Used from common blocks:
C  INFINP : DIRFCK
C  INFORB : N2BASX, NCMOT
C
#include "maxorb.h"
#include "infinp.h"
#include "inforb.h"
C
C Local arrays:
C
      PARAMETER (MXVW = 100)
      DIMENSION ISYMDM(2*MXVW), IFCTYP(2*MXVW)
C
C     Calculate how many V,W matrices we can calculate
C     in one call to SIRFCK.
C     > LVW is work memory needed for calculating on V,W pair.
C     > 1mw is estimate of what is needed in TWOINT(hermit) and the rest
      LVW = 5*N2BASX
C     (5 from 2 DMAT in SCA_FSET, 2 DMAT_AO in HERFCK(hermit)
C      and 2 FMAT in SKLFC1(hermit)).
      NVWMX = (LFREE - 1000000) / LVW
      NVWMX = MIN(NVWMX,NIJ,MXVW)
      IF (NVWMX .LT. 1)
     &    CALL STOPIT('SCA_VW','SIRFCK',(1000000+LVW),LFREE)
C
      KFRSAV = KFREE
      CALL MEMGET('REAL',KDMAT,2*NVWMX*N2BASX,WORK,KFREE,LFREE)
C
C     Zero VWMAT before calculation
C
      CALL DZERO(VWMAT,NIJ*2*N2BASX)
C
      DO IJBLOK = 1,NIJ,NVWMX
         IJEND = MIN(NIJ,IJBLOK+NVWMX-1)
         NVW = IJEND + 1 - IJBLOK
C
C        Construct DMAT_IJ(a,b) density matrices
C
         JDMAT = KDMAT
         JIND  = 1
         DO IJ = IJBLOK,IJEND
            I = IJLIST(1,IJ)
            J = IJLIST(2,IJ)
            CALL SCA_FSET(I,J,WORK(JDMAT),ISYMDM(JIND),IFCTYP(JIND),CMO)
            JDMAT = JDMAT + 2*N2BASX
            JIND  = JIND  + 2
         END DO
C
C        Calculate coulomb + exchange Fock matrices for
C        the current block of (I,J) densities.
C
C        if (dirfck) construct Fock matrices directly
C        else construct Fock matrix from LUINTA
C
         CALL SIRFCK(VWMAT(1,1,IJBLOK),WORK(KDMAT),2*NVW,ISYMDM,IFCTYP,
     &               DIRFCK,WORK(KFREE),LFREE)
C        CALL SIRFCK(FMAT,DMAT,NDMAT,ISYMDM,IFCTYP,DIRECT,WRK,LWORK)
C
#ifdef SCA_VW
C        combine coulomb + exchange Fock matrices to V and W
C
         CALL SCA_F2VW(NVW,VWMAT(1,1,IJBLOK))
#endif
      END DO
C
      CALL MEMREL('SCA_DRV',WORK,1,KFRSAV,KFREE,LFREE)
      RETURN
      END
C  /* Deck sca_fset */
      SUBROUTINE SCA_FSET(I,J,DMAT,ISYMDM,IFCTYP,CMO)
C
C Copyright 3-Nov-1997 Hans Joergen Aa. Jensen.
C
C This subroutine generates information needed
C for calculation of the Coulomb and exchange matrices
C with SIRFCK routine:
#ifndef SCA_VW
C    FMAT(a,b,1) = (IJ | ab)
C    FMAT(a,b,2) = (Ia | Jb)
#else
C    FMAT(a,b,1) = (IJ | ab) = (JI | ab)
C    FMAT(a,b,2) = (Ib | Ja) = (Ja | Ib)
C
C These will in subroutine SCA_F2VW be recombined to give:
C VIJ(a,b) = <IJ || ab> and WIJ(a,b) = <Ia || Jb>.
#endif
C
#include "implicit.h"
      DIMENSION DMAT(N2BASX,2), ISYMDM(2), IFCTYP(2), CMO(NCMOT)
      PARAMETER (DM2 = -2.0D0)
C
C Used from common blocks:
C  INFORB : MULD2H(), N2BASX, NCMOT, ICMO(),NBAS()
C  INFIND : ISMO()
C
#include "maxorb.h"
#include "maxash.h"
#include "inforb.h"
#include "infind.h"
C
      IF (I .EQ. J) THEN
C        DMAT(a,b,*) = CMO(a,I)*CMO(b,I) = DMAT(b,a,*)
         ISYM = ISMO(I)
         IFCTYP(1) = 11
C        symmetric density matrix, coulomb only
         IFCTYP(2) = 12
C        symmetric density matrix, exchange only
         ISYMDM(1) = 1
         ISYMDM(2) = 1
C
         NBASI = NBAS(ISYM)
         JCMOI = ICMO(ISYM) + (I-IORB(ISYM)-1)*NBASI + 1
         JDMAT = NBAST*IBAS(ISYM) + IBAS(ISYM) + 1
         CALL DZERO(DMAT(1,2),N2BASX)
         CALL DGEMM('N','T',NBASI,NBASI,1,1.D0,
     &              CMO(JCMOI),NBASI,
     &              CMO(JCMOI),NBASI,0.D0,
     &              DMAT(JDMAT,2),NBAST)
         CALL DCOPY(N2BASX,DMAT(1,2),1,DMAT(1,1),1)
      ELSE
         ISYM = ISMO(I)
         JSYM = ISMO(J)
         IFCTYP(1) = 11
C        symmetric density matrix, coulomb only
         IFCTYP(2) =  2
C        general density matrix, exchange only
         ISYMDM(1) = MULD2H(ISYM,JSYM)
         ISYMDM(2) = MULD2H(ISYM,JSYM)
C
         NBASI = NBAS(ISYM)
         JCMOI = ICMO(ISYM) + (I-IORB(ISYM)-1)*NBASI + 1
         NBASJ = NBAS(JSYM)
         JCMOJ = ICMO(JSYM) + (J-IORB(JSYM)-1)*NBASJ + 1
#ifndef SCA_VW
         JDMAT = NBAST*IBAS(JSYM) + IBAS(ISYM) + 1
C        DMAT(a,b,2) = CMO(a,I)*CMO(b,J) .ne. DMAT(b,a,2)
         CALL DZERO(DMAT(1,2),N2BASX)
         CALL DGEMM('N','T',NBASI,NBASJ,1,1.D0,
     &              CMO(JCMOI),NBASI,
     &              CMO(JCMOJ),NBASJ,0.D0,
     &              DMAT(JDMAT,2),NBAST)
#else
         JDMAT = NBAST*IBAS(ISYM) + IBAS(JSYM) + 1
C        DMAT(a,b,2) = CMO(a,J)*CMO(b,I) .ne. DMAT(b,a,2)
         CALL DZERO(DMAT(1,2),N2BASX)
         CALL DGEMM('N','T',NBASJ,NBASI,1,1.D0,
     &              CMO(JCMOJ),NBASJ,
     &              CMO(JCMOI),NBASI,0.D0,
     &              DMAT(JDMAT,2),NBAST)
#endif
C        DMAT(a,b,1) is DMAT(a,b,2) symmetrized for Coulomb matrix
         CALL DGETSI(NBAST,DMAT(1,2),DMAT(1,1))
      END IF
      CALL DSCAL(N2BASX,DM2,DMAT(1,2),1)
C     ... IFCTYP(2)=12 calculates -0.5*(I*|*J)
      RETURN
      END
C  /* Deck sca_f2vw */
      SUBROUTINE SCA_F2VW(NFMAT,FMAT)
C
C Copyright 3-Nov-1997 Hans Joergen Aa. Jensen.
C
C On input (calculated in SIRFCK after SCA_FSET):
C    FMAT(a,b,1) = (IJ | ab) = (JI | ab)
C    FMAT(a,b,2) = (Ib | Ja) = (Ja | Ib)
C
C On output (matrices for scattering):
C    FMAT(a,b,1) = <Ia || Jb> = (IJ|ab) - (Ib|Ja)
C                = FMAT(a,b,1) - FMAT(a,b,2)
C    FMAT(a,b,2) = <IJ || ab> = (Ia | Jb) - (Ib | Ja)
C                = FMAT(b,a,2) - FMAT(a,b,2)
C
#include "implicit.h"
C
      DIMENSION FMAT(NBAST,NBAST,2,NFMAT)
C
      PARAMETER (D0 = 0.0D0)
C
C Used from common blocks:
C  INFORB : NBAST
C
#include "inforb.h"
C
      DO IFMAT = 1,NFMAT
C
         DO IB = 1,NBAST
            DO IA = 1,NBAST
               FMAT(IA,IB,1,IFMAT) =
     &         FMAT(IA,IB,1,IFMAT) - FMAT(IA,IB,2,IFMAT)
            END DO
         END DO
C
         DO IB = 1,NBAST
            FMAT(IB,IB,2,IFMAT) = D0
            DO IA = 1,IB-1
               FMAT(IA,IB,2,IFMAT) =
     &         FMAT(IB,IA,2,IFMAT) - FMAT(IA,IB,2,IFMAT)
               FMAT(IB,IA,2,IFMAT) = -FMAT(IA,IB,2,IFMAT)
            END DO
         END DO
C
      END DO
      RETURN
      END
