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
      SUBROUTINE DENFIT_NODDY_INT
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     This routine is a driver for calculating two- and three-center
*     integrals, which are then put on file.
*     
*     (alpha|beta) integrals are needed to solve linear equations for 
*                  the density-fitting coefficeints
*
*     (ab|alpha)   integrals are needed both to solve linear equations 
*                  for the density-fitting coefficeints and to construct
*                  the Coulomb contribution to Fock matrix, J_{ab}
*
*     a,b         - usual basis set
*     alpha,beta  - auxilliary basis for density fitting
*
* 
*     Filip Pawlowski, Oslo, 13.Oct 2004
*
****************************************************************************
*
C
C FIXME: Clean up calling list and declarations after things are working!!!
C        At the moment it's a copy of header files from the REAPRI routine.
C
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
C
#include <abainf.h>
#include <cbirea.h>
C
#include <molinp.h>
#include <ccom.h>
#include <nuclei.h>
#include <frame.h>
#include <primit.h>
#include <shells.h>
#include <symmet.h>
#include <aosotr.h>
#include <chrsgn.h>
#include <inforb.h>

C
      CHARACTER TITLER*14, SPACER*14
      PARAMETER (TITLER = 'NODDY_DF_INT> ')
      PARAMETER (SPACER = '              ')
C
      LOGICAL LDEBUG
      PARAMETER (LDEBUG = .FALSE.)
C
      !FIXME:Current implementation only for s functions
      DO I = 1,KMAX
         IF (NHKT(I).GT.1) CALL QUIT('DenFit-noddy '
     &                          //'implemented for s only')
      END DO
C
      IF (NSYM.GT.1) CALL QUIT('No symmetry in DenFit-noddy')

C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) TITLER//'Entered DENFIT_NODDY_INT routine'
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) TITLER//'Parameters of primitives'
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)' # |        exponent  |      x '
     &                       //'     |      y      |      z '
     &                       //'     | CCF (norm)  '
         DO I = 1,NPBAS
            WRITE(LUPRI,*)'---------------------------------'
     &                          //'-------------------------------'
     &                          //'--------------'
            WRITE(LUPRI,'(I3,A2,F17.10,A3,F11.8,A3,F11.8,A3,F11.8,
     &                    A3,F12.8)')
     &       I,' |',PRIEXP(I),' | ',PRICRX(I),' | ',PRICRY(I),
     &      ' | ',PRICRZ(I),' | ',PRICCF(I,numcft(i))
!FIXME:PRICCF shall refer to actual number of contracted, not "1"
         END DO
*        WRITE(LUPRI,*)SPACER//' # |        exponent  |      x '
*    &                       //'     |      y      |      z '
*        DO I = 1,NPBAS
*           WRITE(LUPRI,*)SPACER//'---------------------------------'
*    &                          //'-------------------------------'
*           WRITE(LUPRI,'(A14,I3,A2,F17.10,A3,F11.8,A3,F11.8,A3,F11.8)')
*    &      SPACER,I,' |',PRIEXP(I),' | ',PRICRX(I),' | ',PRICRY(I),
*    &      ' | ',PRICRZ(I)
*        END DO
      END IF
C
      !Calculate (alpha|beta) integrals and put on file 
      !(needed to solve linear equations for the density-fitting coefficeints)
      CALL DENFIT_NODDY_INT2
C
      !Calculate (ab|alpha) integrals and put on file 
      !(needed both to solve linear equations for the density-fitting 
      !coefficeints and to construct the Coulomb contribution to 
      !Fock matrix, J_{ab})
      CALL DENFIT_NODDY_INT3
C
*     !Calculate (ab|cd) integrals - just to check that my integral
*     !routine works correctly...
      call denfit_noddy_int4
C
      RETURN
      END
C
C  /* Deck denfit_noddy_int2 */
      SUBROUTINE DENFIT_NODDY_INT2
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     This routine calculates two-center (alpha|beta) integrals, which are 
*     then put on file.
*     
*     The integrals are needed to solve linear equations for 
*     the density-fitting coefficeints
*
*     alpha,beta  - auxilliary basis functions for density fitting
*
* 
*     Filip Pawlowski, Oslo, 13.Oct 2004
*
****************************************************************************
*
C
C FIXME: Clean up calling list and declarations after things are working!!!
C        At the moment it's a copy of header files from the REAPRI routine.
C
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
C
#include <abainf.h>
#include <cbirea.h>
C
#include <molinp.h>
#include <ccom.h>
#include <nuclei.h>
#include <frame.h>
#include <primit.h>
#include <shells.h>
#include <symmet.h>
#include <aosotr.h>
#include <chrsgn.h>
#include <hertop.h>

C
      CHARACTER*17 FNAB 
      PARAMETER(FNAB = 'DENFIT_NODDY_INT2')
C
      CHARACTER TITLER*15, SPACER*15
      PARAMETER (TITLER = 'NODDY_DF_INT2> ')
      PARAMETER (SPACER = '               ')
C
      LOGICAL LDEBUG
      PARAMETER (LDEBUG = .FALSE.)
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0)
      PARAMETER (DSMMY = 1.0D-10)
C
      !FIXME You don't need to keep them unless for debugging print!!!
      !      So FAC could be dimensionless
      DIMENSION FAC(NPBAS,NPBAS)
C
      !FIXME: It would probably be enough to give
      !       (JTOP+1)(JTOP+2)/2 which is still fitted only for the highest L
      !       so maybe dimension statment caould be put inside the loops????
      DIMENSION IKM1((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3) 
      DIMENSION IKM2((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3)
C
C FIXME:NPBAS should be replaced by NPBASAUX or sth
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) TITLER//'Entered NODDY_DF_INT2 routine'
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) TITLER//'Number of primitives :',NPBAS
         WRITE(LUPRI,*) TITLER//'Highest angular momentum :',JTOP
         WRITE(LUPRI,*) TITLER//'Number of (alpha|beta) integrals '
     &                        //'to calculate :',NPBAS*NPBAS
         WRITE(LUPRI,*) 
      END IF
C
      !Open the file for (alpha|beta) integrals
      LUAB = -1
      CALL WOPEN2(LUAB,FNAB,64,0)
C
*     write(lupri,*)'(alpha|beta) ints in noddy'
      DO IALPHA = 1,NPBAS
        DO IBETA = 1,NPBAS
         !initialize the integral
         XINT = D0
         !get 2*pi^{5/2}/(pq*sqrt(p+q)) factor, 
         !where p - exp for alpha function, q -exponent for beta function
         CALL GET_FACTOR(FAC(IALPHA,IBETA),PRIEXP(IALPHA),PRIEXP(IBETA))
C
         !Get reduced exponent alpha = p*q/(p+q)
         !where p = a + b ,  q = c + d
         EXPP = PRIEXP(IALPHA) + 0.0D0
         EXPQ = PRIEXP(IBETA)  + 0.0D0
         EXPALPHA = (EXPP*EXPQ)/(EXPP+EXPQ)
C
         !Get R_{PQ}^2 distance between the two overlap distributions
         PX = (PRIEXP(IALPHA)*PRICRX(IALPHA) + 0.0D0*0.0D0)/EXPP
         QX = (PRIEXP(IBETA) *PRICRX(IBETA)  + 0.0D0*0.0D0)/EXPQ
         PY = (PRIEXP(IALPHA)*PRICRY(IALPHA) + 0.0D0*0.0D0)/EXPP
         QY = (PRIEXP(IBETA) *PRICRY(IBETA)  + 0.0D0*0.0D0)/EXPQ
         PZ = (PRIEXP(IALPHA)*PRICRZ(IALPHA) + 0.0D0*0.0D0)/EXPP
         QZ = (PRIEXP(IBETA) *PRICRZ(IBETA)  + 0.0D0*0.0D0)/EXPQ
         XPQ = PX - QX
         IF (ABS(XPQ).LT.DSMMY) XPQ = D0
         YPQ = PY - QY
         IF (ABS(YPQ).LT.DSMMY) YPQ = D0
         ZPQ = PZ - QZ
         IF (ABS(ZPQ).LT.DSMMY) ZPQ = D0
         RPQ2 = XPQ**2 + YPQ**2 + ZPQ**2
*        write(*,*)ialpha,ibeta
*        write(*,*)'pz,qz ',pz,qz
*        write(*,*)'zpq  = ',zpq 
*        write(*,*)'rpq2 = ',rpq2
*        write(*,*)'erf =',erf(sqrt(EXPALPHA*RPQ2))
*        write(*,*)'R0 =',sqrt(3.14/(4*EXPALPHA*RPQ2))*
*    & erf(sqrt(EXPALPHA*RPQ2))
*        write(*,*)'-------------'
C
C        FIXME This part has to be fixed to loop through the different
C        cartesian components. 
C        Instead of "0" pass a general angular momentum.
C        Get out IINDEX from GET_IKM and loop over that!
C        For each IINDEX we should get a different integral
C        (because it's a different primitive)
C
         !Get angular momentum for alpha and beta functions 
         CALL GET_IKM(0,IKM1)
         CALL GET_IKM(0,IKM2)
C
         !Set TUV loop
         NT = IKM1(1,1) + 0
         NU = IKM1(1,2) + 0
         NV = IKM1(1,3) + 0
         !Set TauNuPhi loop
         NTAU = IKM2(1,1) + 0
         NNU  = IKM2(1,2) + 0
         NPHI = IKM2(1,3) + 0
C
         DO IT = 0,NT
          !Get E^{i0}_t
          CALL GET_EIJT(IKM1(1,1),PRIEXP(IALPHA),PRICRX(IALPHA),
     &                  0,0.0D0,0.0D0,IT,EIJT)
C
          DO IU = 0,NU
           !Get E^{k0}_u
           CALL GET_EIJT(IKM1(1,2),PRIEXP(IALPHA),PRICRY(IALPHA),
     &                   0,0.0D0,0.0D0,IU,EKLU)
C
           DO IV = 0,NV
            !Get E^{m0}_v
            CALL GET_EIJT(IKM1(1,3),PRIEXP(IALPHA),PRICRZ(IALPHA),
     &                   0,0.0D0,0.0D0,IV,EMNV)
C
            ETUV = EIJT*EKLU*EMNV
C
            DO ITAU = 0,NTAU
             !Get E^{i0}_t
             CALL GET_EIJT(IKM2(1,1),PRIEXP(IBETA),PRICRX(IBETA),
     &                     0,0.0D0,0.0D0,ITAU,EIJTAU)
C
             DO INU = 0,NNU
              !Get E^{k0}_u
              CALL GET_EIJT(IKM2(1,2),PRIEXP(IBETA),PRICRY(IBETA),
     &                      0,0.0D0,0.0D0,INU,EKLNU)
C
              DO IPHI = 0,NPHI
               !Get E^{m0}_v
               CALL GET_EIJT(IKM2(1,3),PRIEXP(IBETA),PRICRZ(IBETA),
     &                      0,0.0D0,0.0D0,IPHI,EMNPHI)
C
               ETAUNUPHI = EIJTAU*EKLNU*EMNPHI
C
               !Calculate Hermite Coulomb Integrals
               CALL GET_RTUV(IT+ITAU,IU+INU,IV+IPHI,EXPALPHA,RPQ2,RTUV)
*     write(*,*)ialpha,ibeta,rtuv
C
               !Finally... calculate the integral g_abcd
               XINT = XINT
     &              + (-D1)**(ITAU+INU+IPHI)
     &                *ETUV*ETAUNUPHI
     &                *RTUV
C
              END DO !IPHI
             END DO !INU
            END DO !ITAU
C
           END DO !IV
          END DO !IU
         END DO !IT
C
         !Don't forget about the factor in front...
         XINT = FAC(IALPHA,IBETA)*XINT
C
         !You also need to multiply by the contraction coeffs
         !(though this noddy code works at the moment only
         !with uncontracted basis, but still there is some
         !renormalization issue...)
         !FIXME:PRICCF shall refer to actual number of contracted, 
         !not "1"
         XINT = XINT*PRICCF(IALPHA,numcft(ialpha))
     &              *PRICCF(IBETA,numcft(ibeta))
C
        !Write (alpha|beta) integrals to file
        IADR = (IBETA + NPBAS*(IALPHA-1))
        LENGTH = 1
        CALL PUTWA2(LUAB,FNAB,XINT,IADR,LENGTH)
*       write(lupri,*)'XINT(',IALPHA,',',IBETA,') = ',XINT
        END DO !IBETA
      END DO !IALPHA
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) TITLER//'PI-exponents factors:'
         WRITE(LUPRI,*)
         WRITE(LUPRI,*)SPACER//' alpha | beta  | factor '
         DO I = 1,NPBAS
            DO J = 1,NPBAS
               WRITE(LUPRI,*)SPACER//'--------------------------------'
               WRITE(LUPRI,'(A14,I3,A2,I3,A2,F20.10)')
     &         SPACER,I,' | ',J,' | ',FAC(I,J)
            END DO
         END DO
      END IF  
cfp
*     do i =1,npbas
*     do j =1,npbas
*     iadr = (j+npbas*(i-1))
*     call getwa2(luab,fnab,tempoo,iadr,1)
*     write(*,*)'XINTread(',i,',',j,') = ',tempoo
*     end do
*     end do
cfp
C
      !Close and keep the file with (alpha|beta) inmtegrals
      CALL WCLOSE2(LUAB,FNAB,'KEEP')
C
      RETURN
      END
C  /* Deck denfit_noddy_int3 */
      SUBROUTINE DENFIT_NODDY_INT3
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     This routine calculates two-center (ab|alpha) integrals, which are 
*     then put on file.
*     
*     The integrals are needed to solve linear equations for 
*     the density-fitting coefficeints and to construct the J_ab
*     part of the Fock operator
*
*     a,b    - standard basis functions 
*     alpha  - auxilliary basis functions for density fitting
*
* 
*     Filip Pawlowski, Oslo, 14.Oct 2004
*
****************************************************************************
*
C
C FIXME: Clean up calling list and declarations after things are working!!!
C        At the moment it's a copy of header files from the REAPRI routine.
C
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
C
#include <abainf.h>
#include <cbirea.h>
C
#include <molinp.h>
#include <ccom.h>
#include <nuclei.h>
#include <frame.h>
#include <primit.h>
#include <shells.h>
#include <symmet.h>
#include <aosotr.h>
#include <chrsgn.h>
#include <hertop.h>

C
      CHARACTER*17 FNABAL 
      PARAMETER(FNABAL = 'DENFIT_NODDY_INT3')
C
      CHARACTER TITLER*15, SPACER*15
      PARAMETER (TITLER = 'NODDY_DF_INT3> ')
      PARAMETER (SPACER = '               ')
C
      LOGICAL LDEBUG
      PARAMETER (LDEBUG = .FALSE.)
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0)
      PARAMETER (DSMMY = 1.0D-10)
C
      !FIXME: It would probably be enough to give
      !       (JTOP+1)(JTOP+2)/2 which is still fitted only for the highest L
      !       so maybe dimension statment caould be put inside the loops????
      DIMENSION IKM1((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3) 
      DIMENSION JLN1((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3) 
      DIMENSION IKM2((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3)
C
      NPBASAUX = NPBAS
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) TITLER//'Entered NODDY_DF_INT3 routine'
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) TITLER//'Number of primitives in standard '
     &                        //'basis:',NPBAS
         WRITE(LUPRI,*) TITLER//'Number of primitives in auxillary '
     &                        //'basis:',NPBASAUX
         WRITE(LUPRI,*) TITLER//'Highest angular momentum :',JTOP
         WRITE(LUPRI,*) TITLER//'Number of (ab|alpha) integrals '
     &                        //'to calculate :',NPBAS*NPBAS*NPBASAUX
         WRITE(LUPRI,*) 
      END IF
C
      !Open the file for (ab|alpha) integrals
      LUABAL = -1
      CALL WOPEN2(LUABAL,FNABAL,64,0)
C
      DO IA = 1,NPBAS
        DO IB = 1,NPBAS
         DO IALPHA = 1,NPBASAUX
C
          !initialize the integral
          XINT = D0
C
          !get sum of exponents: p = a + b ,  q = alpha + 0
          EXPP = PRIEXP(IA) + PRIEXP(IB)
          EXPQ = PRIEXP(IALPHA)  + 0.0D0
C
          !get 2*pi^{5/2}/(pq*sqrt(p+q)) factor, 
          CALL GET_FACTOR(FAC,EXPP,EXPQ)
C
          !Get reduced exponent expalpha = p*q/(p+q)
          EXPALPHA = (EXPP*EXPQ)/(EXPP+EXPQ)
C
          !Get coordinates of P (=A+B)  overlap distribution
          PX = (PRIEXP(IA)*PRICRX(IA) + PRIEXP(IB)*PRICRX(IB))/EXPP
          PY = (PRIEXP(IA)*PRICRY(IA) + PRIEXP(IB)*PRICRY(IB))/EXPP
          PZ = (PRIEXP(IA)*PRICRZ(IA) + PRIEXP(IB)*PRICRZ(IB))/EXPP
C
          !Get coordinates of Q=(ALPHA+0) overlap distribution
          QX = (PRIEXP(IALPHA)*PRICRX(IALPHA)  + 0.0D0*0.0D0)/EXPQ
          QY = (PRIEXP(IALPHA)*PRICRY(IALPHA)  + 0.0D0*0.0D0)/EXPQ
          QZ = (PRIEXP(IALPHA)*PRICRZ(IALPHA)  + 0.0D0*0.0D0)/EXPQ
C
          !Get components of distance vector between P and Q 
          XPQ = PX - QX
          IF (ABS(XPQ).LT.DSMMY) XPQ = D0
          YPQ = PY - QY
          IF (ABS(YPQ).LT.DSMMY) YPQ = D0
          ZPQ = PZ - QZ
          IF (ABS(ZPQ).LT.DSMMY) ZPQ = D0
C
          !Get R_{PQ}^2 distance between the two overlap distributions (P & Q)
          RPQ2 = XPQ**2 + YPQ**2 + ZPQ**2
*
*         write(*,*)ialpha,ibeta
*         write(*,*)'pz,qz ',pz,qz
*         write(*,*)'zpq  = ',zpq 
*         write(*,*)'rpq2 = ',rpq2
*         write(*,*)'erf =',erf(sqrt(EXPALPHA*RPQ2))
*         write(*,*)'R0 =',sqrt(3.14/(4*EXPALPHA*RPQ2))*
*    & er f(sqrt(EXPALPHA*RPQ2))
*         write(*,*)'-------------'
C
C         FIXME This part has to be fixed to loop through the different
C         cartesian components. 
C         Instead of "0" pass a general angular momentum.
C         Get out IINDEX from GET_IKM and loop over that!
C         For each IINDEX we should get a different integral
C         (because it's a different primitive)
C
          !Get i1, k1, m1 for function a (i1+k1+m1=l_a)
          CALL GET_IKM(0,IKM1)
          !Get j1, l1, n1 for function b (j1+l1+n1=l_b)
          CALL GET_IKM(0,JLN1)
          !Get i2, k2, m2 for function alpha (i2+k2+m2=l_alpha)
          CALL GET_IKM(0,IKM2)
C
          !Set TUV loop
          NT = IKM1(1,1) + JLN1(1,1)
          NU = IKM1(1,2) + JLN1(1,2)
          NV = IKM1(1,3) + JLN1(1,3)
          !Set TauNuPhi loop
          NTAU = IKM2(1,1) + 0
          NNU  = IKM2(1,2) + 0
          NPHI = IKM2(1,3) + 0
C
          DO IT = 0,NT
           !Get E^{i0}_t
           CALL GET_EIJT(IKM1(1,1),PRIEXP(IA),PRICRX(IA),
     &                   JLN1(1,1),PRIEXP(IB),PRICRX(IB),
     &                   IT,EIJT)
C
           DO IU = 0,NU
            !Get E^{k0}_u
            CALL GET_EIJT(IKM1(1,2),PRIEXP(IA),PRICRY(IA),
     &                    JLN1(1,2),PRIEXP(IB),PRICRY(IB),
     &                    IU,EKLU)
C
            DO IV = 0,NV
             !Get E^{m0}_v
             CALL GET_EIJT(IKM1(1,3),PRIEXP(IA),PRICRZ(IA),
     &                     JLN1(1,3),PRIEXP(IB),PRICRZ(IB),
     &                     IV,EMNV)
C
             ETUV = EIJT*EKLU*EMNV
C
             DO ITAU = 0,NTAU
              !Get E^{i0}_t
              CALL GET_EIJT(IKM2(1,1),PRIEXP(IALPHA),PRICRX(IALPHA),
     &                      0,0.0D0,0.0D0,ITAU,EIJTAU)
C
              DO INU = 0,NNU
               !Get E^{k0}_u
               CALL GET_EIJT(IKM2(1,2),PRIEXP(IALPHA),PRICRY(IALPHA),
     &                       0,0.0D0,0.0D0,INU,EKLNU)
C
               DO IPHI = 0,NPHI
                !Get E^{m0}_v
                CALL GET_EIJT(IKM2(1,3),PRIEXP(IALPHA),PRICRZ(IALPHA),
     &                       0,0.0D0,0.0D0,IPHI,EMNPHI)
C
                ETAUNUPHI = EIJTAU*EKLNU*EMNPHI
C
                !Calculate Hermite Coulomb Integrals
                CALL GET_RTUV(IT+ITAU,IU+INU,IV+IPHI,EXPALPHA,RPQ2,RTUV)

*               write(*,*)ialpha,ibeta,rtuv
C
                !Finally... calculate the integral g_abcd
                XINT = XINT
     &               + (-D1)**(ITAU+INU+IPHI)
     &                 *ETUV*ETAUNUPHI
     &                 *RTUV
C
               END DO !IPHI
              END DO !INU
             END DO !ITAU
C
            END DO !IV
           END DO !IU
          END DO !IT
C
          !Don't forget about the factor in front...
          XINT = FAC*XINT

          !You also need to multiply by the contraction coeffs
          !(though this noddy code works at the moment only
          !with uncontracted basis, but still there is some
          !renormalization issue...)
          !FIXME:PRICCF shall refer to actual number of contracted, 
          !not "1"
          XINT = XINT*PRICCF(IA,numcft(ia))*PRICCF(IB,numcft(ib))
     &                *PRICCF(IALPHA,numcft(ialpha))

C
          !Write (ab|alpha) integrals to file
          IADR = (IALPHA + NPBASAUX*(IB-1)+NPBASAUX*NPBAS*(IA-1))
          LENGTH = 1
          CALL PUTWA2(LUABAL,FNABAL,XINT,IADR,LENGTH)
*         write(*,*)'XINT(',IA,',',IB,'|',IALPHA,') = ',XINT
         END DO !IALPHA
        END DO !IB
      END DO !IA
C
cfp
*     do i =1,npbas
*     do j =1,npbas
*     do k =1,npbasaux
*     iadr = (k+npbasaux*(j-1)+npbasaux*npbas*(i-1))
*     call getwa2(LUABAL,FNABAL,tempoo,iadr,1)
*     write(*,*)'XINTread(',i,',',j,'|',k,') = ',tempoo
*     end do
*     end do
*     end do
cfp
C
      !Close and keep the file with (ab|alpha) inmtegrals
      CALL WCLOSE2(LUABAL,FNABAL,'KEEP')
C
      RETURN
      END
C  /* Deck denfit_noddy_int4 */
      SUBROUTINE DENFIT_NODDY_INT4
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     This routine calculates two-center (ab|cd) integrals, which can 
*     then be put on file or printed.
*     
*     The integrals are useful for debugging purposes and SHOULD NOT
*     be calculated in a normal run of the density fitting noddy code.
*
*     Filip Pawlowski, Oslo, 14.Oct 2004
*
****************************************************************************
*
C
C FIXME: Clean up calling list and declarations after things are working!!!
C        At the moment it's a copy of header files from the REAPRI routine.
C
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
C
#include <abainf.h>
#include <cbirea.h>
C
#include <molinp.h>
#include <ccom.h>
#include <nuclei.h>
#include <frame.h>
#include <primit.h>
#include <shells.h>
#include <symmet.h>
#include <aosotr.h>
#include <chrsgn.h>
#include <hertop.h>

C
      CHARACTER*17 FNABCD 
      PARAMETER(FNABCD = 'DENFIT_NODDY_INT4')
C
      CHARACTER TITLER*15, SPACER*15
      PARAMETER (TITLER = 'NODDY_DF_INT4> ')
      PARAMETER (SPACER = '               ')
C
      LOGICAL LDEBUG
      PARAMETER (LDEBUG = .FALSE.)
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0)
      PARAMETER (DSMMY = 1.0D-10)
C
      !FIXME: It would probably be enough to give
      !       (JTOP+1)(JTOP+2)/2 which is still fitted only for the highest L
      !       so maybe dimension statment could be put inside the loops????
      DIMENSION IKM1((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3) 
      DIMENSION JLN1((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3) 
      DIMENSION IKM2((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3)
      DIMENSION JLN2((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3) 
C
*     NPBASAUX = NPBAS
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) TITLER//'Entered NODDY_DF_INT4 routine'
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) TITLER//'Number of primitives in standard '
     &                        //'basis:',NPBAS
*        WRITE(LUPRI,*) TITLER//'Number of primitives in auxillary '
*    &                        //'basis:',NPBASAUX
         WRITE(LUPRI,*) TITLER//'Highest angular momentum :',JTOP
         WRITE(LUPRI,*) TITLER//'Number of (ab|cd) integrals '
     &                        //'to calculate :',NPBAS**4
         WRITE(LUPRI,*) 
      END IF
C
      !Open the file for (ab|cd) integrals
      LUABCD = -1
      CALL WOPEN2(LUABCD,FNABCD,64,0)
C
cfp
*     icounter = 0
cfp
      DO IA = 1,NPBAS
        DO IB = 1,NPBAS
         DO IC = 1,NPBAS
          DO ID = 1,NPBAS
cfp
*     icounter = icounter + 1
cfp
C
           !initialize the integral
           XINT = D0
C
           !get sum of exponents: p = a + b ,  q = c + d
           EXPP = PRIEXP(IA) + PRIEXP(IB)
           EXPQ = PRIEXP(IC) + PRIEXP(ID)
C
           !get 2*pi^{5/2}/(pq*sqrt(p+q)) factor, 
           CALL GET_FACTOR(FAC,EXPP,EXPQ)
C
           !Get reduced exponent expalpha = p*q/(p+q)
           EXPALPHA = (EXPP*EXPQ)/(EXPP+EXPQ)
C
           !Get coordinates of P (=A+B)  overlap distribution
           PX = (PRIEXP(IA)*PRICRX(IA) + PRIEXP(IB)*PRICRX(IB))/EXPP
           PY = (PRIEXP(IA)*PRICRY(IA) + PRIEXP(IB)*PRICRY(IB))/EXPP
           PZ = (PRIEXP(IA)*PRICRZ(IA) + PRIEXP(IB)*PRICRZ(IB))/EXPP
C
           !Get coordinates of Q=(C+D) overlap distribution
           QX = (PRIEXP(IC)*PRICRX(IC) + PRIEXP(ID)*PRICRX(ID))/EXPQ
           QY = (PRIEXP(IC)*PRICRY(IC) + PRIEXP(ID)*PRICRY(ID))/EXPQ
           QZ = (PRIEXP(IC)*PRICRZ(IC) + PRIEXP(ID)*PRICRZ(ID))/EXPQ
C
           !Get components of distance vector between P and Q 
           XPQ = PX - QX
           IF (ABS(XPQ).LT.DSMMY) XPQ = D0
           YPQ = PY - QY
           IF (ABS(YPQ).LT.DSMMY) YPQ = D0
           ZPQ = PZ - QZ
           IF (ABS(ZPQ).LT.DSMMY) ZPQ = D0
C
           !Get R_{PQ}^2 distance between the two overlap distributions (P & Q)
           RPQ2 = XPQ**2 + YPQ**2 + ZPQ**2
*
*          write(*,*)ialpha,ibeta
*          write(*,*)'pz,qz ',pz,qz
*          write(*,*)'zpq  = ',zpq 
*          write(*,*)'rpq2 = ',rpq2
*          write(*,*)'erf =',erf(sqrt(EXPALPHA*RPQ2))
*          write(*,*)'R0 =',sqrt(3.14/(4*EXPALPHA*RPQ2))*
*    & er  f(sqrt(EXPALPHA*RPQ2))
*          write(*,*)'-------------'
C
C          FIXME This part has to be fixed to loop through the different
C          cartesian components. 
C          Instead of "0" pass a general angular momentum.
C          Get out IINDEX from GET_IKM and loop over that!
C          For each IINDEX we should get a different integral
C          (because it's a different primitive)
C
           !Get i1, k1, m1 for function a (i1+k1+m1=l_a)
           CALL GET_IKM(0,IKM1)
           !Get j1, l1, n1 for function b (j1+l1+n1=l_b)
           CALL GET_IKM(0,JLN1)
           !Get i2, k2, m2 for function c (i2+k2+m2=l_c)
           CALL GET_IKM(0,IKM2)
           !Get j2, l2, n2 for function d (j2+l2+n2=l_d)
           CALL GET_IKM(0,JLN2)
C
           !Set TUV loop
           NT = IKM1(1,1) + JLN1(1,1)
           NU = IKM1(1,2) + JLN1(1,2)
           NV = IKM1(1,3) + JLN1(1,3)
           !Set TauNuPhi loop
           NTAU = IKM2(1,1) + JLN2(1,1)
           NNU  = IKM2(1,2) + JLN2(1,2)
           NPHI = IKM2(1,3) + JLN2(1,3)
C
           DO IT = 0,NT
            !Get E^{i0}_t
            CALL GET_EIJT(IKM1(1,1),PRIEXP(IA),PRICRX(IA),
     &                    JLN1(1,1),PRIEXP(IB),PRICRX(IB),
     &                    IT,EIJT)
C
            DO IU = 0,NU
             !Get E^{k0}_u
             CALL GET_EIJT(IKM1(1,2),PRIEXP(IA),PRICRY(IA),
     &                     JLN1(1,2),PRIEXP(IB),PRICRY(IB),
     &                     IU,EKLU)
C
             DO IV = 0,NV
              !Get E^{m0}_v
              CALL GET_EIJT(IKM1(1,3),PRIEXP(IA),PRICRZ(IA),
     &                      JLN1(1,3),PRIEXP(IB),PRICRZ(IB),
     &                      IV,EMNV)
C
              ETUV = EIJT*EKLU*EMNV
C
              DO ITAU = 0,NTAU
               !Get E^{i0}_t
               CALL GET_EIJT(IKM2(1,1),PRIEXP(IC),PRICRX(IC),
     &                       JLN2(1,1),PRIEXP(ID),PRICRX(ID),
     &                       ITAU,EIJTAU)
C
               DO INU = 0,NNU
                !Get E^{k0}_u
                CALL GET_EIJT(IKM2(1,2),PRIEXP(IC),PRICRY(IC),
     &                        JLN2(1,2),PRIEXP(ID),PRICRY(ID),
     &                        INU,EKLNU)
C
                DO IPHI = 0,NPHI
                 !Get E^{m0}_v
                 CALL GET_EIJT(IKM2(1,3),PRIEXP(IC),PRICRZ(IC),
     &                         JLN2(1,3),PRIEXP(ID),PRICRZ(ID),
     &                         IPHI,EMNPHI)
C
                 ETAUNUPHI = EIJTAU*EKLNU*EMNPHI
C
                 !Calculate Hermite Coulomb Integrals
                CALL GET_RTUV(IT+ITAU,IU+INU,IV+IPHI,EXPALPHA,RPQ2,RTUV)
*               if ((ia.eq.1).and.(ib.eq.1).and.(ic.eq.1).and.(id.eq.4))
*    &          then
*               icounter = icounter + 1
*               write(*,*) EXPALPHA*RPQ2
*               write(*,*)'scaledRTUV(',IA,',',IB,'|',IC,',',ID,') = ',
*    &          FAC*RTUV
*               write(*,*)'ETUV = ',ETUV
*               write(*,*)'EMNPHI = ', EMNPHI
*               write(*,*)'ETAUNUPHI = ', ETAUNUPHI
*               write(*,*)'rpq2 = ', RPQ2
*               xmuab = PRIEXP(ia)*PRIEXP(ib)/(PRIEXP(ia)
*    &                               +PRIEXP(ib))
*               write(*,*)'muab = ', xmuab
*               write(*,*)'xab = ', PRICRX(IA) - PRICRX(Ib)
*               write(*,*)'yab = ', PRICRy(IA) - PRICRy(Ib)
*               zab = PRICRz(IA) - PRICRz(Ib)
*               write(*,*)'zab = ', zab
*               xmucd = PRIEXP(ic)*PRIEXP(id)/(PRIEXP(ic)
*    &                               +PRIEXP(id))
*               write(*,*)'mucd = ', xmucd
*               write(*,*)'xcd = ', PRICRX(Ic) - PRICRX(Id)
*               write(*,*)'ycd = ', PRICRy(Ic) - PRICRy(Id)
*               zcd = PRICRz(Ic) - PRICRz(Id)
*               write(*,*)'zcd = ', zcd
*               PI = 2.0D0*ACOS(0.0D0)
*               write(*,*)'expalpha = ', EXPALPHA
*               X4alpi = sqrt((4.0d0*EXPALPHA)/pi)
*               write(*,*)'sqrt(4al/pi) = ', X4alpi
*               write(*,*)
*            write(*,*)'pi/p = ',(pi/expp)
*            write(*,*)'pi/q = ',(pi/expq)
*            write(*,*)'(pi/p) = ',(pi/expp)**(3.0d0/2.0d0)
*            write(*,*)'(pi/q) = ',(pi/expq)**(3.0d0/2.0d0)
*     Sab =(pi/expp)**(3.0d0/2.0d0)*exp(-1.0d0*xmuab*zab**2.0d0)
*     write(*,*)'Sab =',Sab
*     Scd =(pi/expq)**(3.0d0/2.0d0)*exp(-1.0d0*xmucd*zcd**2.0d0)
*     write(*,*)'Scd =',Scd
*     SSab =exp(-1.0d0*xmuab*zab**2.0d0)
*     write(*,*)'Sab/N =',SSab
*     SScd =exp(-1.0d0*xmucd*zcd**2.0d0)
*     write(*,*)'Scd/N =',SScd
*     write(*,*)'X4alpi*Sab*Scd = ',X4alpi*Sab*Scd
*     XR000=SQRT(PI/(4.0d0*EXPALPHA*RPQ2))*ERF(SQRT(EXPALPHA*RPQ2))
*               write(*,*)'XR000 = ',XR000
*               end if

*                write(*,*)ialpha,ibeta,rtuv
C
                 !Finally... calculate the integral g_abcd
                 XINT = XINT
     &                + (-D1)**(ITAU+INU+IPHI)
     &                  *ETUV*ETAUNUPHI
     &                  *RTUV
C
                END DO !IPHI
               END DO !INU
              END DO !ITAU
C
             END DO !IV
            END DO !IU
           END DO !IT
C
           !Don't forget about the factor in front...
           XINT = FAC*XINT
C
           !You also need to multiply by the contraction coeffs
           !(though this noddy code works at the moment only
           !with uncontracted basis, but still there is some
           !renormalization issue...)
           !FIXME:PRICCF shall refer to actual number of contracted, 
           !not "1"
*          write(*,*)'XINT(',IA,',',IB,'|',IC,',',ID,') = ',XINT
           XINT = XINT*PRICCF(IA,numcft(ia))*PRICCF(IB,numcft(ib))
     &                *PRICCF(IC,numcft(ic))*PRICCF(ID,numcft(id))
C
*          !Write (ab|cd) integrals to file
           IADR = (ID + NPBAS*(IC-1)+NPBAS*NPBAS*(IB-1)
     &                + NPBAS*NPBAS*NPBAS*(IA-1))
           LENGTH = 1
           CALL PUTWA2(LUABCD,FNABCD,XINT,IADR,LENGTH)
*          if ((ia.eq.1).and.(ib.eq.1).and.(ic.eq.1).and.(id.eq.4))
*    &          then
*          write(*,*)'XINT/FAC(',IA,',',IB,'|',IC,',',ID,') = '
*    &                ,XINT/fac
*          write(*,*)'XINT(',IA,',',IB,'|',IC,',',ID,') = ',XINT
*          end if
          END DO !ID
         END DO !IC
        END DO !IB
      END DO !IA
cfp
*     write(*,*)'icounter = ',icounter
cfp
C
cfp
*     do i =1,npbas
*     do j =1,npbas
*     do k =1,npbas
*     do l =1,npbas
*     iadr = (l+npbas*(k-1)+npbas*npbas*(j-1)+npbas*npbas*npbas*(i-1))
*     call getwa2(LUABCD,FNABCD,tempoo,iadr,1)
*     write(*,*)'XINTread(',i,',',j,'|',k,',',l,') = ',tempoo
*     end do
*     end do
*     end do
*     end do
cfp
C
      !Close and keep the file with (ab|cd) integrals
      CALL WCLOSE2(LUABCD,FNABCD,'KEEP')
C
      RETURN
      END
C  /* Deck get_factor */
      SUBROUTINE GET_FACTOR(FAC,P,Q)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     This routine calculates the factor in front of the g_{abcd} integral:
*
*     FAC = 2*PI^{5/2} / (PQ*SQRT(P+Q)),
*
*     where P = A + B
*           Q = C + D
*
*           A,B,C,D are exponents of a,b,c,d functions
*
* 
*     Filip Pawlowski, Oslo, 13.Oct 2004
*
****************************************************************************
*
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
C
!#if defined (SYS_CRAY)
!      REAL P,Q,FAC
!      REAL PI,R2PI52,D2,D0
!#else
!      DOUBLE PRECISION P,Q,FAC
!      DOUBLE PRECISION PI,R2PI52,D2,D0
!#endif
      PARAMETER (D0 = 0.0D0, D2 = 2.0D0,D5 = 5.0D0)

      !get PI and 2*PI^{5/2} 
      PI = D2*ACOS(D0)
      R2PI52 = D2*(SQRT(PI)**D5)
C
      !FAC = 2*PI^{5/2} / (PQ*SQRT(P+Q))
      FAC = R2PI52/(P*Q*SQRT(P+Q))
C
      RETURN
      END 
C  /* Deck get_ikm */
      SUBROUTINE GET_IKM(LMAX,IKM)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     For a given angular momentum LMAX get all the possible combinations
*     of cartisian exponents I,K,M:
*
*     LMAX = I + K + M
* 
*     Filip Pawlowski, Oslo, 13.Oct 2004
*
****************************************************************************
*
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
#include <hertop.h>
C
         DIMENSION IKM((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3)
C
         IINDEX = 0
C        IKM(J+1,1) = J
         DO I = LMAX, 0, -1
            DO K = LMAX - I, 0, -1
                  M = LMAX - I - K
                  IINDEX = IINDEX + 1
                  IKM(IINDEX,1) = I
                  IKM(IINDEX,2) = K
                  IKM(IINDEX,3) = M
            END DO
         END DO

C
      RETURN
      END 
C
C  /* Deck get_eijt */
      SUBROUTINE GET_EIJT(I,A,AX,J,B,BX,IT,EIJT)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     Calculate E^{ij}_t using the recurrence relations
* 
*     Filip Pawlowski, Oslo, 13.Oct 2004
*
****************************************************************************
*
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D1 = 1.0D0,D2 = 2.0D0)
C
      !FIXME: At the moment density fitting implemented only for s functions
      IF ( (I.GT.0).OR.(J.GT.0).OR.(IT.GT.0) ) THEN
         WRITE(LUPRI,*) 'I = ', I, '  ,J = ',J,' , T = ',IT
         WRITE(LUPRI,*) 'GET_EIJT only implemented for E000'
         CALL QUIT('Case not implemented in GET_EIJT')
      END IF
C
      !calculate the reduced exponent 
      XMU  = (A*B)/(A+B)
      !calculate the Gaussian separation
      XAB = AX - BX
      !K^{x}_{ab}
      EIJT = EXP(-D1*XMU*XAB**D2)
      RETURN
      END
C  /* Deck get_rtuv */
      SUBROUTINE GET_RTUV(IT,IU,IV,EXPALPHA,RPQ2,RTUV)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     Calculate Hermite Coulomb Integrals R_{t,u,v}(alpha,R_pq)
*     based on error function:
*
*    R_{t,u,v}(alpha,R_pq) = (d/dPx)^t (d/dPy)^u (d/dPz)^v R_{0,0,0}(alpha,R_pq)
*
*     where
*     R_{0,0,0}(alpha,R_pq) = F_0(alpha*R_pq^2)
*
*     F_0 is the Boys function and may be caluclated using the error function
*
*     The derivatives d/dPx, etc are incorporated exploiting recurence
*     relations.
*
!!!!!!!!!!!!!!!!!!!!!!! OBS ! OBS ! OBS ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  RPQ2 is assumed to come in as ALREADY squered R_pq, i.e.:
!
!              RPQ2 = R_pq^2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
* 
*     Filip Pawlowski, Oslo, 13.Oct 2004
*
****************************************************************************
*
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
C
      !FIXME: At the moment density fitting implemented only for s functions
      IF ( (IT.GT.0).OR.(IU.GT.0).OR.(IV.GT.0) ) THEN
         WRITE(LUPRI,*) 'T = ', IT, '  ,U = ',IU,' , V = ',IV
         WRITE(LUPRI,*) 'GET_RTUV only implemented for R000'
         CALL QUIT('Case not implemented in GET_RTUV')
      END IF
C
      !get PI 
      PI = D2*ACOS(D0)
C
      !Calculate the starting point - R000 - based on the error function
      IF ((EXPALPHA*RPQ2).EQ.D0) THEN
         R000 = 1 !F_n(0) = 1/(2n +1); here n=0
      ELSE
         CALL QUIT('Noddy code temporarily disabled')
cluuk    Error function call disabled since not all linkers will find this c-routine
cluuk    R000 =SQRT(PI/(D4*EXPALPHA*RPQ2))*ERF(SQRT(EXPALPHA*RPQ2))
      END IF
      RTUV = R000
C
      RETURN
      END
C  /* Deck denfit_coeffs_noddy */
      SUBROUTINE DENFIT_COEFFS_NODDY(N,NAUX,DMAT,XC)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     
*--------------------------------------------------------------------------
*     
*     PURPOSE:
*     ========
*
*     Calculate coeffcients c_alpha of the densitty-fitting expansion:
*
*       rho^{tilde} = sum_{alpha} c_{alpha} * alpha,
*
*     where 
*       rho^{tilde} - the fitted density 
*       alpha             - auxilliary basis functions
*
*     The coefficeints are calculated by solving the set of linear eqs:
*
*       sum_{alpha,beta} (alpha|beta) * c_{beta} = gamma_{alpha}    (1)
*
*     where
*       gamma_{alpha} = sum_{ab} D_ab (ab|alpha)
*
*       D_ab          - AO density matrix (DMAT(N,N))
*       (alpha|beta)  - integrals calculated in DENFIT_NODDY_INT2 routine
*                       (sitting on 'DENFIT_NODDY_INT2' file)
*       (ab|alpha)    - integrals calculated in DENFIT_NODDY_INT3 routine
*                       (sitting on 'DENFIT_NODDY_INT3' file)
*
*--------------------------------------------------------------------------
*
*     IMPLEMENTATION:
*     ===============
*
*     Eq. (1) above may be written in matrix form as:
*
*       M * C = GAMMA 
*
*     where
*       C     - the solution vector of size NAUX with expansion coeffs c_{alpha}
*       M     - matrix NAUX x NAUX with (alpha|beta) elements
*       GAMMA - vector of size NAUX containing right-hand sides gamma_{alpha}
*
*     The C solution vector is found in this noddy routine in a very
*     primitive and unefficient manner:
*
*       C = M^{-1} * GAMMA
*
*     The real implementation should call a lapack routine (e.g., dposv routine)
*
*            
* 
*     Filip Pawlowski, Oslo, 16.Oct 2004
*
****************************************************************************
*
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
C
      CHARACTER TITLER*21, SPACER*21
      PARAMETER (TITLER = 'DENFIT_COEFFS_NODDY> ')
      PARAMETER (SPACER = '                     ')
C
      LOGICAL LDEBUG
      PARAMETER (LDEBUG = .FALSE.)
C
      DIMENSION DMAT(N,N),XGAMMA(NAUX),XM(NAUX,NAUX),XMINV(NAUX,NAUX)
      DIMENSION XC(NAUX)
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SPACER//'***********************************'
         WRITE(LUPRI,*) TITLER//'Entered DENFIT_COEFFS_NODDY routine'
         WRITE(LUPRI,*) SPACER//'***********************************'
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) TITLER//'Number of primitives in standard '
     &                        //'basis:',N
         WRITE(LUPRI,*) TITLER//'Number of primitives in auxillary '
     &                        //'basis:',NAUX
         WRITE(LUPRI,*)
      END IF
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) SPACER//'------------------'
         WRITE(LUPRI,*) TITLER//'AO Density matirx:'
         WRITE(LUPRI,*) SPACER//'------------------'
         CALL OUTPUT(DMAT,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      !Calculate GAMMA vector
      ! gamma_{alpha} = sum_{ab} D_ab (ab|alpha)
      CALL GET_GAMMA_VEC(N,NAUX,DMAT,XGAMMA)
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) SPACER//'-------------'
         WRITE(LUPRI,*) TITLER//'GAMMA vector:'
         WRITE(LUPRI,*) SPACER//'-------------'
         CALL OUTPUT(XGAMMA,1,NAUX,1,1,NAUX,1,1,LUPRI)
      END IF
C
      !Calculate M matrix
      !M_{alpha,beta} = (alpha|beta)
      CALL GET_M_MAT(XM,NAUX)
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) SPACER//'---------'
         WRITE(LUPRI,*) TITLER//'M matrix:'
         WRITE(LUPRI,*) SPACER//'---------'
         CALL OUTPUT(XM,1,NAUX,1,NAUX,NAUX,NAUX,1,LUPRI)
      END IF
C
      !Invert M matrix: MINV = M^{-1}
      CALL DF_MATINV(XM,XMINV,NAUX)
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) SPACER//'--------------'
         WRITE(LUPRI,*) TITLER//'M^{-1} matrix:'
         WRITE(LUPRI,*) SPACER//'--------------'
         CALL OUTPUT(XMINV,1,NAUX,1,NAUX,NAUX,NAUX,1,LUPRI)
      END IF
C
      !So... we are ready now to find the density fitting coeffs vector:
      !C = M^{-1} * GAMMA
      CALL GET_COEFFS(XMINV,XGAMMA,XC,NAUX)
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) SPACER//'--------------------'
         WRITE(LUPRI,*) TITLER//'Coefficients vector:'
         WRITE(LUPRI,*) SPACER//'--------------------'
         CALL OUTPUT(XC,1,NAUX,1,1,NAUX,1,1,LUPRI)
      END IF
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,*) SPACER//'***********************************'
         WRITE(LUPRI,*) TITLER//'Leaving DENFIT_COEFFS_NODDY routine'
         WRITE(LUPRI,*) SPACER//'***********************************'
      END IF
C
      RETURN
      END
C  /* Deck get_gamma_vec */
      SUBROUTINE GET_GAMMA_VEC(N,NAUX,DMAT,XGAMMA)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     
*--------------------------------------------------------------------------
*     
*     Calculate GAMMA vector needed for solving the equations for
*     density fitting coeffs (see DENFIT_COEFFS_NODDY routine).
*
*     The elements of the vector are:
*
*       gamma_{alpha} = sum_{ab} D_ab (ab|alpha)
*
*     where
*       D_ab          - AO density matrix (DMAT(N,N))
*       (ab|alpha)    - integrals calculated in DENFIT_NODDY_INT3 routine
*                       (sitting on 'DENFIT_NODDY_INT3' file)
*
* 
*     Filip Pawlowski, Oslo, 16.Oct 2004
*
****************************************************************************
*
#include <implicit.h>
#include <priunit.h>
C
      CHARACTER*17 FNABAL
      PARAMETER(FNABAL = 'DENFIT_NODDY_INT3')
C
      LOGICAL LDEBUG
      PARAMETER(LDEBUG = .FALSE.)
C
      DIMENSION DMAT(N,N),XGAMMA(NAUX)
C
      !Open the file for (ab|alpha) integrals
      LUABAL = -1
      CALL WOPEN2(LUABAL,FNABAL,64,0)
*
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) 'AO density matrix in GET_GAMMA_VEC:'
        CALL OUTPUT(DMAT,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      !Initialize XGAMMA vector
      CALL DZERO(XGAMMA,NAUX)
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) 'XGAMMA at the beginning of GET_GAMMA_VEC:'
        CALL OUTPUT(XGAMMA,1,NAUX,1,1,NAUX,1,1,LUPRI)
      END IF
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) '(ab|alpha) integrals in GET_GAMMA_VEC:'
      END IF
C
      !Calculate gamma_{alpha} = sum_{ab} D_ab (ab|alpha)
      DO IA =1,N
         DO IB =1,N
            DO IALPHA =1,NAUX
C
               !Get (ab|alpha) integral
               LENGTH = 1
               IADR = (IALPHA+NAUX*(IB-1)+NAUX*N*(IA-1))
               CALL GETWA2(LUABAL,FNABAL,XINT3,IADR,LENGTH)
C
               IF (LDEBUG) THEN
                 WRITE(LUPRI,*)'(',IA,',',IB,'|',IALPHA,') = ',XINT3
               END IF
C
               !Contract: gamma_{alpha} = gamma_{alpha} 
               !                        + D_ab * (ab|alpha)
               XGAMMA(IALPHA) = XGAMMA(IALPHA) + XINT3*DMAT(IA,IB)
C
            END DO !IALPHA
         END DO !IB
      END DO !IA
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 'XGAMMA at the end of GET_GAMMA_VEC:'
         CALL OUTPUT(XGAMMA,1,NAUX,1,1,NAUX,1,1,LUPRI)
      END IF
C
      !Close and keep the file with (ab|alpha) integrals
      CALL WCLOSE2(LUABAL,FNABAL,'KEEP')

      RETURN
      END
C  /* Deck get_m_mat */
      SUBROUTINE GET_M_MAT(XM,NAUX)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     
*--------------------------------------------------------------------------
*     
*     Calculate M matrix needed for solving the equations for
*     density fitting coeffs (see DENFIT_COEFFS_NODDY routine).
*
*     The elements of the matrix are:
*
*       M_{alpha,beta} = (alpha|beta)
*
*     so we simply read in the (alpha|beta) integrals calculated 
*     in DENFIT_NODDY_INT2 routine (sitting on 'DENFIT_NODDY_INT2' file)
*
* 
*     Filip Pawlowski, Oslo, 16.Oct 2004
*
****************************************************************************
*
#include <implicit.h>
#include <priunit.h>
C
      CHARACTER*17 FNAB
      PARAMETER(FNAB = 'DENFIT_NODDY_INT2')
C
      LOGICAL LDEBUG
      PARAMETER (LDEBUG = .FALSE.)
C
      DIMENSION XM(NAUX,NAUX)
C
      !Open the file for (alpha|beta) integrals
      LUAB = -1
      CALL WOPEN2(LUAB,FNAB,64,0)
C
      !Construct M matrix 
      ! ( M_{alpha,beta} = (alpha|beta) )
      DO IALPHA = 1,NAUX
         DO IBETA = 1,NAUX
            LENGTH = 1
            IADR = (IBETA+NAUX*(IALPHA-1))
            CALL GETWA2(LUAB,FNAB,XM(IALPHA,IBETA),IADR,LENGTH)
         END DO
      END DO
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)'M matrix in GET_M_MAT:'
         CALL OUTPUT(XM,1,NAUX,1,NAUX,NAUX,NAUX,1,LUPRI)
      END IF
C
      RETURN
      END 
C  /* Deck matinv */
      SUBROUTINE DF_MATINV(A,AINV,N)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     
*
*     Invert matrix A(N,N):
*
*        AINV = A^{-1}
*
*     (Refers to DGEINV routine from pdpack/linextra.F)
*
*     OBS!!! Note that it could easily be a real (not noddy) routine
*            IF the arrays were allocated dynamically (using the 
*            workspace WORK).
* 
*     Filip Pawlowski, Oslo, 16.Oct 2004
*
****************************************************************************
*
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D1 = 1.0D0,D0 = 0.0D0)
C
      LOGICAL LDEBUG
      PARAMETER (LDEBUG = .FALSE.)
C
      DIMENSION A(N,N),AINV(N,N),IPVT(N),WRKV(N)
      DIMENSION C(N,N)
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)'Input matrix in DF_MATINV:'
         CALL OUTPUT(A,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      !INVERT
      CALL DGEINV(N,A,AINV,IPVT,WRKV,INFO)
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)'Output matrix in DF_MATINV:'
         CALL OUTPUT(AINV,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      IF (LDEBUG) THEN
         !Check that A * AINV = I
         DO I = 1,N
            DO J = 1,N
               C(I,J) = D0
            END DO
         END DO
*        write(lupri,*)'c after initialization...'
*        call output(c,1,n,1,n,n,n,1,lupri)
         CALL DGEMM('N','N',N,N,N,
     &              D1,A,MAX(1,N),
     &              AINV,MAX(1,N),
     &              D1,C,MAX(1,N))
         WRITE(LUPRI,*)'DF_MATINV: Checking that A * AINV = I'
         CALL OUTPUT(C,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck get_coeffs */
      SUBROUTINE GET_COEFFS(XMINV,XGAMMA,XC,NAUX)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     
*     Calculate densitty fitting coefficients vector as:
*
*       C(NAUX) = MINV(NAUX,NAUX) * XGAMMA(NAUX)
*
*
*     Filip Pawlowski, Oslo, 16.Oct 2004
*
****************************************************************************
*
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (D1 = 1.0D0)
C
      LOGICAL LDEBUG
      PARAMETER (LDEBUG = .FALSE.)
C
      DIMENSION XMINV(NAUX,NAUX),XGAMMA(NAUX),XC(NAUX)
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)'XMINV matrix in GET_COEFFS:'
         CALL OUTPUT(XMINV,1,NAUX,1,NAUX,NAUX,NAUX,1,LUPRI)
      END IF
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)'XGAMMA vector in GET_COEFFS:'
         CALL OUTPUT(XGAMMA,1,NAUX,1,1,NAUX,1,1,LUPRI)
      END IF
C
      !Initialize the C vector (using beta = 0.0d0 in DGEMV
      !                         leads sometimes to problems!!!)
      CALL DZERO(XC,NAUX)
*     write(lupri,*)'coeffs vector after initialization:'
*     call output(xc,1,naux,1,1,naux,1,1,lupri)
C     
      !Multiply C(NAUX) = MINV(NAUX,NAUX) * XGAMMA(NAUX)
      CALL DGEMV('N',NAUX,NAUX,
     &           D1,XMINV,MAX(1,NAUX),
     &           XGAMMA,1,
     &           D1,XC,1)
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*)'Result (XC) vector in GET_COEFFS:'
         CALL OUTPUT(XC,1,NAUX,1,1,NAUX,1,1,LUPRI)
      END IF
C
      RETURN
      END
C  /* Deck jab_noddy */
      SUBROUTINE JAB_NODDY(N,DMAT,XJAB)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     
*--------------------------------------------------------------------------
*     
*     Calculate two-electron part (Coulomb only) of the Fock operator:
*
*       J_ab = sum_cd (ab|cd) * D_cd
*
*     where
*       D_cd     - AO density matrix (DMAT(N,N))
*       (ab|cd)  - integrals calculated in DENFIT_NODDY_INT4 routine
*                  (sitting on 'DENFIT_NODDY_INT4' file)
*
*     The goal of this routine is to debug the noddy code.
*     It serves two main purposes:
*
*      - checking (by comparison to the J_ab from the real code)
*        that the (ab|cd) integrals are calculated correctly
*        in the DENFIT_NODDY_INT4 routine (then these integrals
*        may be used to further check the integrals coming out
*        from DENFIT_NODDY_INT2 and DENFIT_NODDY_INT3 routines)
*
*      - checking how good J_ab is approximated by J_{ab}^{tilde}
*        from JAB_DENFIT_NODDY routine.
*
* 
*     Filip Pawlowski, Oslo, 17.Oct 2004
*
****************************************************************************
*
#include <implicit.h>
#include <priunit.h>
C
      CHARACTER*17 FNABCD
      PARAMETER(FNABCD = 'DENFIT_NODDY_INT4')
C
      PARAMETER (D0 = 0.0D0)
C
      LOGICAL LDEBUG
      PARAMETER(LDEBUG = .FALSE.)
C
      DIMENSION DMAT(N,N),XJAB(N,N)
C
      !Open the file for (ab|cd) integrals
      LUABCD = -1
      CALL WOPEN2(LUABCD,FNABCD,64,0)
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) 'AO density matrix in JAB_NODDY:'
        CALL OUTPUT(DMAT,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      !Initialize XJAB vector
      DO IA = 1,N
         DO IB = 1,N
            XJAB(IA,IB) = D0
         END DO
      END DO
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) 'XJAB at the beginning of JAB_NODDY:'
        CALL OUTPUT(XJAB,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) '(ab|cd) integrals in JAB_NODDY:'
      END IF
C
      !Calculate J_ab = sum_{cd} (ab|cd) * D_cd 
      DO IA =1,N
         DO IB =1,N
            DO IC =1,N
             DO ID =1,N
C
                 !Get (ab|cd) integral
                 LENGTH = 1
                 IADR = (ID+N*(IC-1)+N*N*(IB-1)+N*N*N*(IA-1))
                 CALL GETWA2(LUABCD,FNABCD,XINT4,IADR,LENGTH)
C
                 IF (LDEBUG) THEN
                  WRITE(LUPRI,*)'(',IA,',',IB,'|',IC,',',ID,') = ',XINT4
                 END IF
C
                 !Contract: J_ab = J_ab + (ab|cd) * D_cd
                 XJAB(IA,IB) = XJAB(IA,IB) + XINT4*DMAT(IC,ID)
C
C
             END DO !ID
            END DO !IC
         END DO !IB
      END DO !IA
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) 'XJAB at the end of JAB_NODDY:'
        CALL OUTPUT(XJAB,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      !Close and keep the file with (ab|cd) integrals
      CALL WCLOSE2(LUABCD,FNABCD,'KEEP')

      RETURN
      END
C  /* Deck jab_denfit_noddy */
      SUBROUTINE JAB_DENFIT_NODDY(N,NAUX,C,XJAB)
*
****************************************************************************
*
*     This is a part of the density-fitting noddy code. 
*     
*--------------------------------------------------------------------------
*     
*     Calculate two-electron part (Coulomb only) of the Fock operator
*     as defined by density fitting:
*
*       J_ab = sum_{alpha} (ab|alpha) * c_{alpha}
*
*     where
*       c_{alpha}    - density fitting coefficients (calculated in
*                      GET_COEFFS routine): C(NAUX)
*       (ab|alpha)   - integrals calculated in DENFIT_NODDY_INT3 routine
*                      (sitting on 'DENFIT_NODDY_INT3' file)
*
*       a,b   - regular basis-set functions
*       alpha - auxiliary basis-set functions (fitting functions)
*
* 
*     Filip Pawlowski, Oslo, 18.Oct 2004
*
****************************************************************************
*
#include <implicit.h>
#include <priunit.h>
C
      CHARACTER*17 FNABAL
      PARAMETER(FNABAL = 'DENFIT_NODDY_INT3')
C
      PARAMETER (D0 = 0.0D0)
C
      LOGICAL LDEBUG
      PARAMETER(LDEBUG = .FALSE.)
C
      DIMENSION XJAB(N,N),C(NAUX)
C
      !Open the file for (ab|alpha) integrals
      LUABAL = -1
      CALL WOPEN2(LUABAL,FNABAL,64,0)
C
      !Initialize XJAB vector
      DO IA = 1,N
         DO IB = 1,N
            XJAB(IA,IB) = D0
         END DO
      END DO
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) 'XJAB at the beginning of JAB_DENFIT_NODDY:'
        CALL OUTPUT(XJAB,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) 'C coeffs at the beginning of JAB_DENFIT_NODDY:'
        CALL OUTPUT(C,1,NAUX,1,1,NAUX,1,1,LUPRI)
      END IF
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) '(ab|alpha) integrals in JAB_DENFIT_NODDY:'
      END IF
C
      !Calculate J_ab = sum_{alpha} (ab|alpha) * c_{alpha} 
      DO IA =1,N
         DO IB =1,N
            DO IALPHA =1,NAUX
C
                 !Get (ab|cd) integral
                 LENGTH = 1
                 IADR = (IALPHA+NAUX*(IB-1)+NAUX*N*(IA-1))
                 CALL GETWA2(LUABAL,FNABAL,XINT3,IADR,LENGTH)
C
                 IF (LDEBUG) THEN
                  WRITE(LUPRI,*)'(',IA,',',IB,'|',IALPHA,') = ',XINT3
                 END IF
C
                 !Contract: J_ab = J_ab + (ab|cd) * D_cd
                 XJAB(IA,IB) = XJAB(IA,IB) + XINT3*C(IALPHA)
C
C
            END DO !IALPHA
         END DO !IB
      END DO !IA
C
      IF (LDEBUG) THEN
        WRITE(LUPRI,*) 'XJAB at the end of JAB_DENFIT_NODDY:'
        CALL OUTPUT(XJAB,1,N,1,N,N,N,1,LUPRI)
      END IF
C
      !Close and keep the file with (ab|alpha) integrals
      CALL WCLOSE2(LUABAL,FNABAL,'KEEP')

      RETURN
      END
C  /* Deck denfit_noddy_int_debugger */
      SUBROUTINE DENFIT_NODDY_INT_DEBUGGER
*
****************************************************************************
*ACTUALLY... This is a routine to debug 2- and 3-index integrals.
*
* It is based on DENFIT_NODDY_INT4 routine.
*
*
*--------------------------------------------------------------------------
*     This is a part of the density-fitting noddy code. 
*     This routine calculates two-center (ab|cd) integrals, which can 
*     then be put on file or printed.
*     
*     The integrals are useful for debugging purposes and SHOULD NOT
*     be calculated in a normal run of the density fitting noddy code.
*
*     Filip Pawlowski, Oslo, 14.Oct 2004
*
****************************************************************************
*
C
C FIXME: Clean up calling list and declarations after things are working!!!
C        At the moment it's a copy of header files from the REAPRI routine.
C
C     IMPLICIT NONE
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
#include <maxorb.h>
#include <aovec.h>
#include <maxaqn.h>
C
#include <abainf.h>
#include <cbirea.h>
C
#include <molinp.h>
#include <ccom.h>
#include <nuclei.h>
#include <frame.h>
#include <primit.h>
#include <shells.h>
#include <symmet.h>
#include <aosotr.h>
#include <chrsgn.h>
#include <hertop.h>

C
      CHARACTER*17 FNABCD 
      PARAMETER(FNABCD = 'DENFIT_NODDY_INT_DEBUGGER')
C
      CHARACTER TITLER*17, SPACER*17
      PARAMETER (TITLER = 'NODDY_DF_INT_DB> ')
      PARAMETER (SPACER = '               ')
C
      LOGICAL LDEBUG
      PARAMETER (LDEBUG = .FALSE.)
C
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0)
      PARAMETER (DSMMY = 1.0D-10)
C
      !FIXME: It would probably be enough to give
      !       (JTOP+1)(JTOP+2)/2 which is still fitted only for the highest L
      !       so maybe dimension statment could be put inside the loops????
      DIMENSION IKM1((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3) 
      DIMENSION JLN1((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3) 
      DIMENSION IKM2((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3)
      DIMENSION JLN2((JTOP+1)*(JTOP+2)*(JTOP+3)/6,3) 
C
*     NPBASAUX = NPBAS
C
      IF (LDEBUG) THEN
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) TITLER//'Entered NODDY_DF_INT_DEBUGER routine'
         WRITE(LUPRI,*) 
         WRITE(LUPRI,*) TITLER//'Number of primitives in standard '
     &                        //'basis:',NPBAS
*        WRITE(LUPRI,*) TITLER//'Number of primitives in auxillary '
*    &                        //'basis:',NPBASAUX
         WRITE(LUPRI,*) TITLER//'Highest angular momentum :',JTOP
         WRITE(LUPRI,*) TITLER//'Number of (ab|cd) integrals '
     &                        //'to calculate :',NPBAS**4
         WRITE(LUPRI,*) 
      END IF
C
      !Open the file for (ab|cd) integrals
      LUABCD = -1
      CALL WOPEN2(LUABCD,FNABCD,64,0)
C
cfp
*     icounter = 0
cfp
      DO IA = 1,NPBAS
cint2
*       DO IB = 1,NPBAS
cint2
         DO IC = 1,NPBAS
cint3
*         DO ID = 1,NPBAS
cint3
cfp
*     icounter = icounter + 1
cfp
C
           !initialize the integral
           XINT = D0
C
           !get sum of exponents: p = a + b ,  q = c + d
cint2
*          EXPP = PRIEXP(IA) + PRIEXP(IB)
           EXPP = PRIEXP(IA) 
cint2
cint3
*          EXPQ = PRIEXP(IC) + PRIEXP(ID)
           EXPQ = PRIEXP(IC) 
cint3
C
           !get 2*pi^{5/2}/(pq*sqrt(p+q)) factor, 
           CALL GET_FACTOR(FAC,EXPP,EXPQ)
C
           !Get reduced exponent expalpha = p*q/(p+q)
           EXPALPHA = (EXPP*EXPQ)/(EXPP+EXPQ)
C
           !Get coordinates of P (=A+B)  overlap distribution
cint2
*          PX = (PRIEXP(IA)*PRICRX(IA) + PRIEXP(IB)*PRICRX(IB))/EXPP
*          PY = (PRIEXP(IA)*PRICRY(IA) + PRIEXP(IB)*PRICRY(IB))/EXPP
*          PZ = (PRIEXP(IA)*PRICRZ(IA) + PRIEXP(IB)*PRICRZ(IB))/EXPP
           PX = (PRIEXP(IA)*PRICRX(IA) )/EXPP
           PY = (PRIEXP(IA)*PRICRY(IA) )/EXPP
           PZ = (PRIEXP(IA)*PRICRZ(IA) )/EXPP
cint2
C
           !Get coordinates of Q=(C+D) overlap distribution
cint3
*          QX = (PRIEXP(IC)*PRICRX(IC) + PRIEXP(ID)*PRICRX(ID))/EXPQ
*          QY = (PRIEXP(IC)*PRICRY(IC) + PRIEXP(ID)*PRICRY(ID))/EXPQ
*          QZ = (PRIEXP(IC)*PRICRZ(IC) + PRIEXP(ID)*PRICRZ(ID))/EXPQ
           QX = (PRIEXP(IC)*PRICRX(IC))/EXPQ
           QY = (PRIEXP(IC)*PRICRY(IC))/EXPQ
           QZ = (PRIEXP(IC)*PRICRZ(IC))/EXPQ
cint3
C
           !Get components of distance vector between P and Q 
           XPQ = PX - QX
           IF (ABS(XPQ).LT.DSMMY) XPQ = D0
           YPQ = PY - QY
           IF (ABS(YPQ).LT.DSMMY) YPQ = D0
           ZPQ = PZ - QZ
           IF (ABS(ZPQ).LT.DSMMY) ZPQ = D0
C
           !Get R_{PQ}^2 distance between the two overlap distributions (P & Q)
           RPQ2 = XPQ**2 + YPQ**2 + ZPQ**2
*
*          write(*,*)ialpha,ibeta
*          write(*,*)'pz,qz ',pz,qz
*          write(*,*)'zpq  = ',zpq 
*          write(*,*)'rpq2 = ',rpq2
*          write(*,*)'erf =',erf(sqrt(EXPALPHA*RPQ2))
*          write(*,*)'R0 =',sqrt(3.14/(4*EXPALPHA*RPQ2))*
*    & er  f(sqrt(EXPALPHA*RPQ2))
*          write(*,*)'-------------'
C
C          FIXME This part has to be fixed to loop through the different
C          cartesian components. 
C          Instead of "0" pass a general angular momentum.
C          Get out IINDEX from GET_IKM and loop over that!
C          For each IINDEX we should get a different integral
C          (because it's a different primitive)
C
           !Get i1, k1, m1 for function a (i1+k1+m1=l_a)
           CALL GET_IKM(0,IKM1)
           !Get j1, l1, n1 for function b (j1+l1+n1=l_b)
           CALL GET_IKM(0,JLN1)
           !Get i2, k2, m2 for function c (i2+k2+m2=l_c)
           CALL GET_IKM(0,IKM2)
           !Get j2, l2, n2 for function d (j2+l2+n2=l_d)
           CALL GET_IKM(0,JLN2)
C
           !Set TUV loop
           NT = IKM1(1,1) + JLN1(1,1)
           NU = IKM1(1,2) + JLN1(1,2)
           NV = IKM1(1,3) + JLN1(1,3)
           !Set TauNuPhi loop
           NTAU = IKM2(1,1) + JLN2(1,1)
           NNU  = IKM2(1,2) + JLN2(1,2)
           NPHI = IKM2(1,3) + JLN2(1,3)
C
           DO IT = 0,NT
            !Get E^{i0}_t
cint2
*           CALL GET_EIJT(IKM1(1,1),PRIEXP(IA),PRICRX(IA),
*    &                    JLN1(1,1),PRIEXP(IB),PRICRX(IB),
*    &                    IT,EIJT)
            CALL GET_EIJT(IKM1(1,1),PRIEXP(IA),PRICRX(IA),
     &                    JLN1(1,1),0.0d0,0.0d0,
     &                    IT,EIJT)
cint2
C
            DO IU = 0,NU
             !Get E^{k0}_u
cint2
*            CALL GET_EIJT(IKM1(1,2),PRIEXP(IA),PRICRY(IA),
*    &                     JLN1(1,2),PRIEXP(IB),PRICRY(IB),
*    &                     IU,EKLU)
             CALL GET_EIJT(IKM1(1,2),PRIEXP(IA),PRICRY(IA),
     &                     JLN1(1,2),0.0d0,0.0d0,
     &                     IU,EKLU)
cint2
C
             DO IV = 0,NV
              !Get E^{m0}_v
cint2
*             CALL GET_EIJT(IKM1(1,3),PRIEXP(IA),PRICRZ(IA),
*    &                      JLN1(1,3),PRIEXP(IB),PRICRZ(IB),
*    &                      IV,EMNV)
              CALL GET_EIJT(IKM1(1,3),PRIEXP(IA),PRICRZ(IA),
     &                      JLN1(1,3),0.0d0,0.0d0,
     &                      IV,EMNV)
cint2
C
              ETUV = EIJT*EKLU*EMNV
C
              DO ITAU = 0,NTAU
               !Get E^{i0}_t
cint3
*              CALL GET_EIJT(IKM2(1,1),PRIEXP(IC),PRICRX(IC),
*    &                       JLN2(1,1),PRIEXP(ID),PRICRX(ID),
*    &                       ITAU,EIJTAU)
               CALL GET_EIJT(IKM2(1,1),PRIEXP(IC),PRICRX(IC),
     &                       JLN2(1,1),0.0d0,0.0d0,
     &                       ITAU,EIJTAU)
cint3
C
               DO INU = 0,NNU
                !Get E^{k0}_u
cint3
*               CALL GET_EIJT(IKM2(1,2),PRIEXP(IC),PRICRY(IC),
*    &                        JLN2(1,2),PRIEXP(ID),PRICRY(ID),
*    &                        INU,EKLNU)
                CALL GET_EIJT(IKM2(1,2),PRIEXP(IC),PRICRY(IC),
     &                        JLN2(1,2),0.0d0,0.0d0,
     &                        INU,EKLNU)
cint3
C
                DO IPHI = 0,NPHI
                 !Get E^{m0}_v
cint3
*                CALL GET_EIJT(IKM2(1,3),PRIEXP(IC),PRICRZ(IC),
*    &                         JLN2(1,3),PRIEXP(ID),PRICRZ(ID),
*    &                         IPHI,EMNPHI)
                 CALL GET_EIJT(IKM2(1,3),PRIEXP(IC),PRICRZ(IC),
     &                         JLN2(1,3),0.0d0,0.0d0,
     &                         IPHI,EMNPHI)
cint3
C
                 ETAUNUPHI = EIJTAU*EKLNU*EMNPHI
C
                 !Calculate Hermite Coulomb Integrals
                CALL GET_RTUV(IT+ITAU,IU+INU,IV+IPHI,EXPALPHA,RPQ2,RTUV)
*               if ((ia.eq.1).and.(ib.eq.1).and.(ic.eq.1).and.(id.eq.4))
*    &          then
*               icounter = icounter + 1
*               write(*,*) EXPALPHA*RPQ2
*               write(*,*)'scaledRTUV(',IA,',',IB,'|',IC,',',ID,') = ',
*    &          FAC*RTUV
*               write(*,*)'ETUV = ',ETUV
*               write(*,*)'EMNPHI = ', EMNPHI
*               write(*,*)'ETAUNUPHI = ', ETAUNUPHI
*               write(*,*)'rpq2 = ', RPQ2
*               xmuab = PRIEXP(ia)*PRIEXP(ib)/(PRIEXP(ia)
*    &                               +PRIEXP(ib))
*               write(*,*)'muab = ', xmuab
*               write(*,*)'xab = ', PRICRX(IA) - PRICRX(Ib)
*               write(*,*)'yab = ', PRICRy(IA) - PRICRy(Ib)
*               zab = PRICRz(IA) - PRICRz(Ib)
*               write(*,*)'zab = ', zab
*               xmucd = PRIEXP(ic)*PRIEXP(id)/(PRIEXP(ic)
*    &                               +PRIEXP(id))
*               write(*,*)'mucd = ', xmucd
*               write(*,*)'xcd = ', PRICRX(Ic) - PRICRX(Id)
*               write(*,*)'ycd = ', PRICRy(Ic) - PRICRy(Id)
*               zcd = PRICRz(Ic) - PRICRz(Id)
*               write(*,*)'zcd = ', zcd
*               PI = 2.0D0*ACOS(0.0D0)
*               write(*,*)'expalpha = ', EXPALPHA
*               X4alpi = sqrt((4.0d0*EXPALPHA)/pi)
*               write(*,*)'sqrt(4al/pi) = ', X4alpi
*               write(*,*)
*            write(*,*)'pi/p = ',(pi/expp)
*            write(*,*)'pi/q = ',(pi/expq)
*            write(*,*)'(pi/p) = ',(pi/expp)**(3.0d0/2.0d0)
*            write(*,*)'(pi/q) = ',(pi/expq)**(3.0d0/2.0d0)
*     Sab =(pi/expp)**(3.0d0/2.0d0)*exp(-1.0d0*xmuab*zab**2.0d0)
*     write(*,*)'Sab =',Sab
*     Scd =(pi/expq)**(3.0d0/2.0d0)*exp(-1.0d0*xmucd*zcd**2.0d0)
*     write(*,*)'Scd =',Scd
*     SSab =exp(-1.0d0*xmuab*zab**2.0d0)
*     write(*,*)'Sab/N =',SSab
*     SScd =exp(-1.0d0*xmucd*zcd**2.0d0)
*     write(*,*)'Scd/N =',SScd
*     write(*,*)'X4alpi*Sab*Scd = ',X4alpi*Sab*Scd
*     XR000=SQRT(PI/(4.0d0*EXPALPHA*RPQ2))*ERF(SQRT(EXPALPHA*RPQ2))
*               write(*,*)'XR000 = ',XR000
*               end if

*                write(*,*)ialpha,ibeta,rtuv
C
                 !Finally... calculate the integral g_abcd
                 XINT = XINT
     &                + (-D1)**(ITAU+INU+IPHI)
     &                  *ETUV*ETAUNUPHI
     &                  *RTUV
C
                END DO !IPHI
               END DO !INU
              END DO !ITAU
C
             END DO !IV
            END DO !IU
           END DO !IT
C
           !Don't forget about the factor in front...
           XINT = FAC*XINT
C
           !You also need to multiply by the contraction coeffs
           !(though this noddy code works at the moment only
           !with uncontracted basis, but still there is some
           !renormalization issue...)
           !FIXME:PRICCF shall refer to actual number of contracted, 
           !not "1"
*          write(*,*)'XINT_db(',IA,',',IB,'|',IC,',',0,') = ',XINT
*          write(*,*)'XINT_db(',IA,',',0,'|',IC,',',0,') = ',XINT
cint3
*          XINT = XINT*PRICCF(IA,numcft(ia))*PRICCF(IB,numcft(ib))
*    &                *PRICCF(IC,numcft(ic))*PRICCF(ID,numcft(id))
cint2
*          XINT = XINT*PRICCF(IA,numcft(ia))*PRICCF(IB,numcft(ib))
*    &                *PRICCF(IC,numcft(ic))
           XINT = XINT*PRICCF(IA,numcft(ia))*PRICCF(IC,numcft(ic))
cint2
cint3
C
*          !Write (ab|cd) integrals to file
*          IADR = (ID + NPBAS*(IC-1)+NPBAS*NPBAS*(IB-1)
*    &                + NPBAS*NPBAS*NPBAS*(IA-1))
*          LENGTH = 1
*          CALL PUTWA2(LUABCD,FNABCD,XINT,IADR,LENGTH)
*          if ((ia.eq.1).and.(ib.eq.1).and.(ic.eq.1).and.(id.eq.4))
*    &          then
*          write(*,*)'XINT/FAC(',IA,',',IB,'|',IC,',',ID,') = '
*    &                ,XINT/fac
*          write(*,*)'XINTcont_db(',IA,',',IB,'|',IC,',',0,') = ',XINT
*          write(*,*)'XINTcont_db(',IA,',',0,'|',IC,',',0,') = ',XINT
*          end if
cint3
*         END DO !ID
cint3
         END DO !IC
cint2
*       END DO !IB
cint2
      END DO !IA
cfp
*     write(*,*)'icounter = ',icounter
cfp
C
cfp
*     do i =1,npbas
*     do j =1,npbas
*     do k =1,npbas
*     do l =1,npbas
*     iadr = (l+npbas*(k-1)+npbas*npbas*(j-1)+npbas*npbas*npbas*(i-1))
*     call getwa2(LUABCD,FNABCD,tempoo,iadr,1)
*     write(*,*)'XINTread(',i,',',j,'|',k,',',l,') = ',tempoo
*     end do
*     end do
*     end do
*     end do
cfp
C
      !Close and keep the file with (ab|cd) integrals
      CALL WCLOSE2(LUABCD,FNABCD,'KEEP')
C
      RETURN
      END
