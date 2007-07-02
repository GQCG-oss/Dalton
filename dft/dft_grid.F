C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
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
      SUBROUTINE RADLMG(RADNDE,RADWT,NR,RADERR,NRADPT,NUCORB,AA,IPRINT)
#include <implicit.h>
      PARAMETER(D0 = 0D0, D1 = 1D0,D2 = 2D0, D3 = 3D0)
#include <dummy.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <ccom.h>
#include <nuclei.h>
#include <priunit.h>
      DIMENSION NUCORB(NHTYP,2),AA(2,NHTYP,2)
      DIMENSION RADWT(NRADPT), RADNDE(NRADPT)
      CHARACTER SPDCAR*1
c     
C     Grid spacing to H and inner grid point to AH
      IF(IPRINT.GE.3)WRITE(LUPRI,'(A)') '* Grid spacing'
      NR = 0
      H  = DUMMY
      AH = 0D0
      DO LL = 1,NHTYP
         L = LL-1
         NBAS=NUCORB(LL,1)+NUCORB(LL,2)
         IF(NBAS.GT.0) THEN
            HTMP = DISERR(RADERR,L)
            H = MIN(H,HTMP)
            IF(IPRINT.GE.3) THEN
               WRITE(LUPRI,'(3X,A1,A,F6.3)')
     &              SPDCAR(L),'-orbitals --> ',HTMP
            ENDIF
         ENDIF
c        AH = MAX(AA(1,LL,1),AA(1,LL,2)) 
         AH = MAX(AH,AA(1,LL,1))
      ENDDO
      IF(AH .EQ. 0d0) RETURN
      EPH = EXP(H)
      IF(IPRINT.GE.3)WRITE(LUPRI,'(1X,A,F6.3)') 'Value chosen: ',H
C...  Inner grid point AA->R transformation.
      AH = D2*AH
      IF(IPRINT.GE.3)WRITE(LUPRI,*) 'AH = ',AH
      RL = ((1.9D0+LOG(RADERR))/D3)-(LOG(AH)/D2)
      RL = EXP(RL)
      IF(IPRINT.GE.3)
     &     WRITE(LUPRI,'(A,1P,E12.5)') '* Inner grid point:',RL
C...  Outer point
      IF(IPRINT.GE.3) WRITE(LUPRI,'(A)') '* Outer point:'
      RH = D0
      DO LL = 1,NHTYP
         L = LL-1
         AL=DUMMY
         IF(NUCORB(LL,1).GT.0) AL=AA(2,LL,1)
         IF(NUCORB(LL,2).GT.0) AL=MIN(AL,AA(2,LL,2))
         IF(AL.LT.DUMMY) THEN
            AL = AL+AL
            RHTMP = OUTERR(AL,L,RADERR)           
            RH=MAX(RH,RHTMP)
            IF(IPRINT.GE.3) THEN
               WRITE(LUPRI,'(3X,A1,A,F6.3)')
     &              SPDCAR(L),'-orbitals --> ',RHTMP
            ENDIF
         ENDIF
      ENDDO
      IF(IPRINT.GE.3)WRITE(LUPRI,'(1X,A,F6.3)') 'Value chosen: ',RH
      GRDC = RL/(EPH-D1)
      IF(IPRINT.GE.3)WRITE(LUPRI,'(A,E12.5)') 'Constant c:',GRDC
      NR = NINT(LOG(D1+(RH/GRDC))/H)
      IF(IPRINT.GE.3)WRITE(LUPRI,'(A,I5)') 'Number of points:',NR
      IF(NR.GT.NRADPT) CALL QUIT('Too many radial points.')
      RADNDE(NR) = RL
      RADWT(NR)  = (RL+GRDC)*RL*RL*H
      DO IR = NR-1,1,-1
         RADNDE(IR) = (RADNDE(IR+1)+GRDC)*EPH-GRDC
         RADWT(IR) = (RADNDE(IR)+GRDC)*RADNDE(IR)*RADNDE(IR)*H
      ENDDO
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
      FUNCTION DISERR(RD,L)
C                                                                      C
C     Provide grid spacing h for given angular momentum L              C
C     and discretization error RD                                      C
C                                                                      C
C     Based on eqs. (17) and (18) of                                   C
C       R. Lindh, P.-Aa. Malmqvist and L. Gagliardi                    C 
C       "Molecular integrals by numerical quadrature",                 C
C       Theor. Chem. Acc. 106 (2001) 178-187                           C
C                                                                      C
C     The array CF(4,L) contains coefficients of a 3rd order           C
C     polynomial fit to provide start values for the                   C
C     determination of H by a Newton-Raphson search.                   C
C                                                                      C
C     Written by T. Saue July 2002                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(ACC=1.0D-5)
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,TI=1.0D1)
      PARAMETER (PI = 3.14159 26535 89793 D00)
      PARAMETER(MXIT=20)
      DIMENSION CF(4,0:4)
      DATA CF/0.91570D0,0.78806D-1,0.28056D-2,3.4197D-05,
     &        0.74912D0,0.61502D-1,0.21558D-2,2.6100D-05,
     &        0.65449D0,0.52322D-1,0.18217D-2,2.2004D-05,
     &        0.59321D0,0.46769D-1,0.16261D-2,1.9649D-05,
     &        0.55125D0,0.43269D-1,0.15084D-2,1.8270D-05/
C
C     Initialization
C
      FAC  = SQRT(D2)*D2*D2
      IFAC = 1
      DO I = 1,L
        FAC   = FAC*D2
        IFAC  = IFAC*(2*I+1)
      ENDDO
      FAC = FAC/IFAC
      LM = MIN(L,4)
      RDLOG = LOG(RD)
      DISERR = POLVAL(3,CF(1,LM),RDLOG)
      HTLOG = LOG(DISERR)
C     Newton-Raphson search
      DO IT = 1,MXIT
        PIH  = PI/DISERR
        PIHL = PIH
        PIEX = PI*PIH/D2
        DO I = 1,L
          PIHL = PIHL*PIH
        ENDDO
        U0   = FAC*PIHL*EXP(-PIEX)
        U1   = U0*((PIEX/DISERR)-(L+1)/PIH)
        F0   = LOG(U0)-RDLOG
        F1   = DISERR*U1/U0
        DX = F0/F1
        HTLOG = HTLOG - DX
        DISERR = EXP(HTLOG)
        IF(ABS(DX).LT.ACC) RETURN
      ENDDO
      STOP 'Error in DISERR'

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
      FUNCTION OUTERR(AL,L,RD)
C                                                                      C
C     Provide outer grid point for given angular momentum L            C
C     outer exponent AL and discretization error RD                    C
C                                                                      C
C     Based on eq. (19) of                                             C
C       R. Lindh, P.-Aa. Malmqvist and L. Gagliardi                    C
C       "Molecular integrals by numerical quadrature",                 C
C       Theor. Chem. Acc. 106 (2001) 178-187                           C
C                                                                      C
C     The variable U = AL*R*R is found by a Newton-Raphson search.     C 
C                                                                      C
C     Written by T. Saue July 2002                                     C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(ACC=1.0D-6)
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,TI=1.0D1)
      PARAMETER (PI = 3.14159 26535 89793 D00)
      PARAMETER(MXIT=8)
C
C     Initialization
C
      TOLEN = D2
      FAC = D1
      DO I = 1,L
        TOLEN = TOLEN*D2
        FAC   = FAC*(2*I+1)
      ENDDO
      EXPL = (2*L+1)/D2
      A = SQRT(PI)*FAC/TOLEN
      ALN = LOG(A)
      RLN = LOG(RD)
      U = 35.0D0
C     Newton-Raphson search
      DO IT = 1,MXIT
        F0HLN = ALN+EXPL*LOG(U)-U-RLN
        F1HLN = EXPL/U-D1
        DX = F0HLN/F1HLN
        U = U - DX
        IF(ABS(DX).LT.ACC) THEN
          OUTERR = SQRT(U/AL)
          RETURN
        ENDIF
      ENDDO
      STOP 'Error in OUTERR'

      RETURN
      END
      FUNCTION GETMAXL()
#include <implicit.h>
#include <maxaqn.h>
#include <ccom.h>
      INTEGER GETMAXL
      GETMAXL = NHTYP
      END
      SUBROUTINE NUCBAS(NUCORB,AA,IPRINT)
C***********************************************************************
C
C     Extract basis information for all centers
C
C     Written by T.Saue March 12 2001
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0=0.0D0)
C
#include <maxaqn.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <dummy.h>
C
#include <nuclei.h>
#include <ccom.h>
#include <shells.h>
#include <primit.h>
      CHARACTER SPDCAR*1
      DIMENSION NUCORB(NHTYP,2,NUCIND),AA(2,NHTYP,2,NUCIND)
C
C     Initialize
C
      NDIM = 2*(NHTYP)*NUCIND
      CALL IZERO(NUCORB,NDIM)
C
      JCENT = 0
      JPRIM = -1
      JC = -1
      JLVAL = -1
      DO ISHELL = 1,KMAX
        ICENT = NCENT(ISHELL)
        IF(ICENT.NE.JCENT) THEN
          JCENT = ICENT
          JLVAL = 0
        ENDIF
        IC = LCLASS(ISHELL)
        IF(IC.NE.JC) THEN
          JC    = IC
          JLVAL = 0
        ENDIF
        ILVAL = NHKT(ISHELL)
        IF(ILVAL.NE.JLVAL) THEN
          JLVAL = ILVAL
          NUCORB(ILVAL,IC,ICENT) = 0
          AA(1,ILVAL,IC,ICENT)=D0
          AA(2,ILVAL,IC,ICENT)=DUMMY
        ENDIF
        NUCORB(ILVAL,IC,ICENT)=NUCORB(ILVAL,IC,ICENT)+1
        IPRIM = JSTRT(ISHELL)
        IF(IPRIM.NE.JPRIM) THEN
          JPRIM = IPRIM
          NPRIM = NUCO(ISHELL)
          DO IEXP = 1,NPRIM
            A=PRIEXP(IPRIM+IEXP)
            AA(1,ILVAL,IC,ICENT)=MAX(AA(1,ILVAL,IC,ICENT),A)
            AA(2,ILVAL,IC,ICENT)=MIN(AA(2,ILVAL,IC,ICENT),A)
          ENDDO
        ENDIF
      ENDDO
C
      IF(IPRINT.GE.2) THEN
         CALL HEADER('NUCORB:Basis set information:',-1)
        DO I = 1,NUCIND
          WRITE(LUPRI,'(/A,A4,A/)') '*** Center: ',NAMN(I),' ***'
          WRITE(LUPRI,'(2X,A)') '* Large components:'
          IC = 1
          DO LL = 1,NHTYP
          IF(NUCORB(LL,IC,I).GT.0) THEN
            L=LL-1
            WRITE(LUPRI,'(3X,A1,A,I6,2(3X,A,E12.5))')
     &      SPDCAR(L),'-orbitals:',NUCORB(LL,IC,I),
     &      'Alpha_H :',AA(1,LL,IC,I),
     &      'Alpha_L :',AA(2,LL,IC,I)
          ENDIF
          ENDDO
          WRITE(LUPRI,'(2X,A)') '* Small components:'
          IC = 2
          DO LL = 1,NHTYP
          IF(NUCORB(LL,IC,I).GT.0) THEN
            L=LL-1
            WRITE(LUPRI,'(3X,A1,A,I6,2(3X,A,E12.5))')
     &      SPDCAR(L),'-orbitals:',NUCORB(LL,IC,I),
     &      'Alpha_H: ',AA(1,LL,IC,I),
     &      'Alpha_L: ',AA(2,LL,IC,I)
          ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
c     Trond: polynomial evaluation.
      FUNCTION POLVAL(NORDER,B,XVAL)
#include <implicit.h>
      DIMENSION B(NORDER+1)
      POLVAL = B(1)
      XBUF = 1
      DO 10 I = 2,(NORDER+1)
        XBUF = XBUF*XVAL
        POLVAL = POLVAL + XBUF*B(I)
   10 CONTINUE     
      END
